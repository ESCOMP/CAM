!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_cc_emissions
!
! !DESCRIPTION: Module hco\_cc\_emissions provides emissions to CAM-chem in CESM.
!               This module replaces mo\_extfrc and mo\_srf\_emissions in CAM-chem
!               when HEMCO is enabled, and are otherwise stubs if HEMCO-CESM is not
!               enabled at runtime.
!
!               These subroutines emulate the behavior of extfrc\_set (3-D emissions)
!               and set\_srf\_emissions in CAM-chem for API compatibility.
!\\
!\\
! !INTERFACE:
!
module hco_cc_emissions
!
! !USES:
!
	! CESM types
	use shr_kind_mod,     only: r8 => shr_kind_r8

	! Run control
	use spmd_utils,       only: masterproc
	use cam_abortutils,   only: endrun
	use cam_logfile,      only: iulog

	! Grid information
	use ppgrid,           only: pver

	! Chemistry mechanism properties
    use chem_mods,        only: gas_pcnst
    use mo_tracname,      only: solsym

    ! Physics buffer operations
	use physics_buffer, only: physics_buffer_desc
    use physics_buffer, only : pbuf_get_field, pbuf_get_index

	implicit none
	private
!
! !PUBLIC MEMBER FUNCTIONS:
!
	public :: hco_set_srf_emissions, hco_set_extfrc
!
! !REMARKS:
!  None
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_set_srf_emissions
!
! !DESCRIPTION: Sets surface level emissions (note, TOA is level 1, so this is level pver)
!               for this latitude slice (chunk).
!\\
!\\
! !INTERFACE:
!
    subroutine hco_set_srf_emissions( lchnk, ncol, sflx, pbuf )
!
! !USES:
!
    implicit none
!
! !INPUT PARAMETERS:
!
    integer,  intent(in)               :: ncol        ! columns in chunk
    integer,  intent(in)               :: lchnk       ! chunk index
    real(r8), intent(out)              :: sflx(:,:)   ! surface emissions ( kg/m^2/s )
    type(physics_buffer_desc), pointer :: pbuf        ! pbuf in chunk
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(r8), pointer     :: pbuf_ik(:,:)              ! ptr to pbuf data (/pcols,pver/)
    integer               :: tmpIdx                    ! pbuf field id
    character(len=255)    :: fldname_ns                ! field name HCO_NH3
    integer               :: RC                        ! return code (dummy)

 	!--------------------------------------------------------
    ! ... set HEMCO_CESM emissions
    ! hplin 7/19/20
    !--------------------------------------------------------

    ! for every species index retrieve the species name, compute the pbuf name,
    ! and write it into sflx(col, n)
    ! where n is spc_ndx

    ! if the pbuf exists, set has_emis(1:gas_pcnst) to .true.
    ! this process is supposed to be set by srf_emissions_inti but it is just
    ! used below, so we shunt it here and decide later

    ! ncol: # of columns in chunk
    ! lchnk: chunk number

    ! sflx is given in (pcols, gas_pcnst) so it is a in-chunk slice of the
    ! srf flux specifier. maybe this loop needs to be done higher up so
    ! we loop over the pbuf to prevent inquiries. tbd hplin 7/19/20

    do n = 1, gas_pcnst
        ! species name: solsym(n)
        fldname_ns = 'HCO_' // trim(solsym(n))
        tmpIdx = pbuf_get_index(fldname_ns, RC)
        if ( tmpIdx < 0 ) then
            if ( masterproc ) write(iulog,*) "mo_srf_emissions hemco: Field not found ", TRIM(fldname_ns)
        else
            ! this is already in chunk, retrieve it
            call pbuf_get_field(pbuf, tmpIdx, pbuf_ik)

            if(.not. associated(pbuf_ik)) then ! sanity check
              call endrun("mo_srf_emissions hemco: FATAL - tmpIdx > 0 but pbuf_ik unassoc")
            endif

            ! for each col retrieve data from pbuf_ik(I, K)
            sflx(1:ncol,n) = pbuf_ik(1:ncol,pver) ! only surface emissions for now, remember vertical is inverted

            has_emis(n) = .true.

            ! if(masterproc) write(iulog,*) "mo_srf_emissions hemco: debug added emiss for", solsym(n), maxval(sflx(1:ncol, n))
        endif
    enddo

	end subroutine hco_set_srf_emissions
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_set_extfrc
!
! !DESCRIPTION: Set 3-D emissions apart from surface level.
!\\
!\\
! !INTERFACE:
!
    subroutine hco_set_extfrc( lchnk, zint, frcing, ncol, pbuf )
!
! !USES:
!
		implicit none
		use mo_chem_utls,   only: get_spc_ndx
	    use mo_chem_utls,   only: get_extfrc_ndx

	    ! Check list whether this species has external forcing from dataset - this is a CAM-chem flag
	    ! and this is CAM-chem specific.
	    use chem_mods,      only: frc_from_dataset, extfrc_lst
	    use chem_mods,      only: extcnt, adv_mass

	    use mo_constants,   only: avogadro
!
! !INPUT PARAMETERS:
!
	    integer,  intent(in)               :: ncol                       ! columns in chunk
	    integer,  intent(in)               :: lchnk                      ! chunk index
	    real(r8), intent(in)               :: zint(ncol, pverp)          ! interface geopot above surface (km)
	    real(r8), intent(inout)            :: frcing(ncol,pver,extcnt)   ! insitu forcings (molec/cm^3/s)
        type(physics_buffer_desc), pointer :: pbuf                       ! pbuf in chunk
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version based on original from 14 Nov 2020
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(r8), parameter   :: cm2_to_m2 = 1.e4_r8
    real(r8), parameter   :: kg_to_g   = 1.e-3_r8
    real(r8), parameter   :: km_to_cm  = 1.e5_r8

    ! Loop idxs
    integer               :: m, n, k

    ! For compatibility with XFRC_ diagnostic in CAM-chem
    real(r8)              :: frcing_col(1:ncol), frcing_col_kg(1:ncol)
    character(len=max_fieldname_len) :: xfcname
    real(r8)              :: molec_to_kg


    real(r8), pointer     :: pbuf_ik(:,:)              ! ptr to pbuf data (/pcols,pver/)
    integer               :: tmpIdx                    ! pbuf field id
    character(len=255)    :: fldname_ns                ! field name HCO_NH3
    integer               :: RC                        ! return code (dummy)
    real(r8)              :: kg_to_molec

    ! for every species index retrieve the species name, compute the pbuf name,
    ! and write it into frcing(col, lev, n)
    ! where n is FRC_IDX --
    !

    !******************************************************************************************************
    ! WARNING: ONLY SPECIES THAT ARE EXTERNALLY FORCED AND SPECIFIED IN mo_sim_dat.F90
    ! CAN HAVE 3-D EMISSIONS, OTHERWISE THEY WILL BE IGNORED!!!
    !
    !  the n = frc_idx comes from get_extfrc_ndx( spc_name ).
    !  it is too computationally expensive to check all fields to see if they are 3-d emitted
    !
    !  so PLEASE verify that all your species are in mo_sim_dat.F90:: extfrc_lst before attempting
    !  to have 3-D emissions for them.
    !
    ! I will still loop through all species, get its symbol and attempt to inject extfrc emissions
    ! for it, but it may not be guaranteed to be done.
    !******************************************************************************************************

    ! if the pbuf exists, set has_emis(1:gas_pcnst) to .true.
    ! this process is supposed to be set by srf_emissions_inti but it is just
    ! used below, so we shunt it here and decide later

    ! ncol: # of columns in chunk
    ! lchnk: chunk number

    do n = 1, gas_pcnst
      ! check if extfrc available?
      m = get_extfrc_ndx(trim(solsym(n)))

      if(m > 0) then
        ! confirm extfrc present
        has_emis(m) = .true.

        ! add extfrc
        ! "external insitu forcing" (1/cm^3/s) -- NOTE UNITS COMING OUT OF HEMCO are
        ! kg/m2/s, so unit conversion must be done
        !
        ! using species factor...
        ! (kg_to_g is actually kg/g...)
        !
        !     1 / (kg/molec cm2/m2) = molec/kg m2/cm2
        !
        ! kg/m2/s * molec/kg m2/cm2 = molec/cm2/s
        ! now divide by z-interface height (in CM!) for each height to get the right answer!
        ! (hplin, 11/14/20)
        kg_to_molec = 1/(adv_mass(n) / avogadro * cm2_to_m2 * kg_to_g)

        ! species name: solsym(n)
        fldname_ns = 'HCO_' // trim(solsym(n))
        ! if(masterproc) write(iulog,*) "mo_extfrc hemco: Adding extfrc for", fldname_ns

        tmpIdx = pbuf_get_index(fldname_ns, RC)
        if(tmpIdx < 0) then
          if(masterproc) then
            write(iulog,*) "mo_extfrc hemco: Field not found ", TRIM(fldname_ns)
          endif
        else
          ! this is already in chunk, retrieve it
          call pbuf_get_field(pbuf, tmpIdx, pbuf_ik)

          if(.not. associated(pbuf_ik)) then ! sanity check
            call endrun("mo_extfrc hemco: FATAL - tmpIdx > 0 but pbuf_ik unassoc")
          endif

          ! for each col retrieve data from pbuf_ik(I, K)
          do k = 1, pver-1
            frcing(:ncol,k,m) = pbuf_ik(1:ncol,k) * kg_to_molec / ((zint(:ncol,k)-zint(:ncol,k+1)) * km_to_cm)
          enddo
          ! remember vertical is inverted - REMOVE the top level as it is injected in mo_srf_emissions instead

          if ( frc_from_dataset(m) ) then
             xfcname = trim(extfrc_lst(m))//'_XFRC'
             call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )
             spc_ndx = get_spc_ndx( extfrc_lst(m) )
             molec_to_kg = adv_mass( spc_ndx ) / avogadro * cm2_to_m2 * kg_to_g

             frcing_col(:ncol) = 0._r8
             frcing_col_kg(:ncol) = 0._r8
             do k = 1, pver
                frcing_col(:ncol) = frcing_col(:ncol) + frcing(:ncol,k,m)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm
                frcing_col_kg(:ncol) = frcing_col_kg(:ncol) + frcing(:ncol,k,m)*(zint(:ncol,k)-zint(:ncol,k+1))*km_to_cm*molec_to_kg
             enddo

             xfcname = trim(extfrc_lst(m))//'_CLXF'
             call outfld( xfcname, frcing_col(:ncol), ncol, lchnk )
             xfcname = trim(extfrc_lst(m))//'_CMXF'
             call outfld( xfcname, frcing_col_kg(:ncol), ncol, lchnk )
             if ( masterproc ) then
                 write(iulog,*) "mo_extfrc hemco: debug added 3D emiss for ", TRIM(solsym(n)), maxval(frcing(:ncol,:,m))
             endif
          endif
        endif
      endif
    enddo

	end subroutine hco_set_extfrc
!EOC
end module hco_cc_emissions