!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_cc_emissions
!
! !DESCRIPTION: Module hco\_cc\_emissions provides emissions to CAM-chem in CESM.
!               This module replaces mo\_extfrc and mo\_srf\_emissions in CAM-chem
!               when HEMCO is enabled, otherwise they are not called if HEMCO-CESM is not
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
    use ppgrid,           only: pver, pverp

    ! Chemistry mechanism properties
    use chem_mods,        only: gas_pcnst
    use mo_tracname,      only: solsym

    ! Physics buffer operations
    use physics_buffer,   only: physics_buffer_desc
    use physics_buffer,   only: pbuf_get_field, pbuf_get_index

    ! Compat with XFRC diagn
    use cam_history,         only: outfld
    use cam_history_support, only: max_fieldname_len

    implicit none
    private
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public :: hco_extfrc_inti
    public :: hco_set_srf_emissions, hco_set_extfrc
!
! !REMARKS:
!  None
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version
!  04 Nov 2022 - H.P. Lin    - Now initialize extfrc fields in here
!EOP
!------------------------------------------------------------------------------
!BOC
    logical :: pcnst_is_extfrc(gas_pcnst)              ! Is external forcing? (3-D data)
    integer :: pcnst_extfrc_ndx(gas_pcnst)             ! External forcing index from get_extfrc_ndx
    integer :: hco_pbuf_idx(gas_pcnst)                 ! Physics buffer indices for HCO_* fields from HEMCO
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
    type(physics_buffer_desc), pointer :: pbuf(:)     ! pbuf in chunk
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version
!  12 Jan 2023 - H.P. Lin    - Check if pbuf is 2-D or 3-D first
!  09 Feb 2023 - H.P. Lin    - For 3-D pbuf, no longer set cflx and use 3-D forcing only.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer               :: n

    real(r8), pointer     :: pbuf_ptr_3d(:,:)          ! ptr to pbuf data (/pcols,pver/)
    real(r8), pointer     :: pbuf_ptr_2d(:)            ! ptr to pbuf data (/pcols/)
    integer               :: tmpIdx                    ! pbuf field id

    ! reset sflx here. (same as mo_srf_emissions.F90)
    ! sflx is defined in chem_emissions (chemistry.F90) but without default values, and is
    ! later added to cam_in%cflx. it must be initialized in this subroutine.
    sflx(:,:) = 0._r8

    !--------------------------------------------------------
    ! ... set HEMCO emissions
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
        tmpIdx = hco_pbuf_idx(n)
        if(tmpIdx > 0) then
            if(pcnst_is_extfrc(n)) then ! 3-D data
                ! if species is 3-D data, then all forcings set through 3-D. no longer process
                ! their emissions here.
            else ! 2-D data
                call pbuf_get_field(pbuf, tmpIdx, pbuf_ptr_2d)

                ! for each col retrieve data from pbuf_ptr(I, K)
                sflx(1:ncol,n) = pbuf_ptr_2d(1:ncol)

                pbuf_ptr_2d => null()
            endif
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
        use mo_chem_utls,   only: get_spc_ndx

        ! Check list whether this species has external forcing from dataset - this is a CAM-chem flag
        ! and this is CAM-chem specific.
        use chem_mods,      only: frc_from_dataset, extfrc_lst
        use chem_mods,      only: extcnt, adv_mass

        use mo_constants,   only: avogadro
        implicit none
!
! !INPUT PARAMETERS:
!
        integer,  intent(in)               :: ncol                       ! columns in chunk
        integer,  intent(in)               :: lchnk                      ! chunk index
        real(r8), intent(in)               :: zint(ncol, pverp)          ! interface geopot above surface (km)
        real(r8), intent(inout)            :: frcing(ncol,pver,extcnt)   ! insitu forcings (molec/cm^3/s)
        type(physics_buffer_desc), pointer :: pbuf(:)                    ! pbuf in chunk
!
! !REVISION HISTORY:
!  08 Aug 2022 - H.P. Lin    - Initial version based on original from 14 Nov 2020
!  09 Feb 2023 - H.P. Lin    - Use full 3-D emissions, including surface, if available
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
    integer               :: spc_ndx


    real(r8), pointer     :: pbuf_ik(:,:)              ! ptr to pbuf data (/pcols,pver/)
    integer               :: tmpIdx                    ! pbuf field id
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

    ! Zero out frcing to be consistent with mo_extfrc
    frcing(:,:,:) = 0._r8

    do n = 1, gas_pcnst
      ! check if extfrc available?
      if(pcnst_is_extfrc(n)) then
        ! add extfrc
        ! "external insitu forcing" (1/cm^3/s)
        m = pcnst_extfrc_ndx(n)
        tmpIdx = hco_pbuf_idx(n)
        if(tmpIdx > 0) then
          ! Note: units coming out of HEMCO are in kg/m2/s, so unit conversion must be done
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

          ! this is already in chunk, retrieve it.
          ! if the field does not exist, pbuf_get_field will return an error, so sanity check for pbuf_ik is not needed.
          call pbuf_get_field(pbuf, tmpIdx, pbuf_ik)

          ! for each col retrieve data from pbuf_ik(I, K)
          ! this includes surface layer.
          do k = 1, pver
            frcing(:ncol,k,m) = frcing(:ncol,k,m) + pbuf_ik(1:ncol,k) * kg_to_molec / ((zint(:ncol,k)-zint(:ncol,k+1)) * km_to_cm)
          enddo

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
          endif
        endif
      endif
    enddo

    end subroutine hco_set_extfrc
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_extfrc_inti
!
! !DESCRIPTION: Initialize external forcing related diagnostic fields
!\\
!\\
! !INTERFACE:
!
    subroutine hco_extfrc_inti( )
!
! !USES:
!
        use chem_mods,    only: frc_from_dataset, extcnt, extfrc_lst
        use cam_history,  only: addfld, add_default, horiz_only
        use phys_control, only: phys_getopts
        use mo_chem_utls, only: get_extfrc_ndx
        implicit none
!
! !REVISION HISTORY:
!  04 Nov 2022 - H.P. Lin    - Initial version based on extfrc_inti
!  10 Apr 2023 - H.P. Lin    - Now move pcnst_is_extfrc initialization here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    logical  :: history_aerosol
    logical  :: history_chemistry
    logical  :: history_cesm_forcing
    character(len=16)  :: spc_name
    integer  :: n

    character(len=255)    :: fldname_ns                ! field name HCO_NH3
    integer               :: RC                        ! return code (dummy)

    ! for first run, cache results of 3-D or 2-D scan within pcnst_is_extfrc
    ! to avoid lengthy lookups in future timesteps. hplin 1/12/23
    do n = 1, gas_pcnst
        pcnst_extfrc_ndx(n) = get_extfrc_ndx(trim(solsym(n)))
        pcnst_is_extfrc(n) = (pcnst_extfrc_ndx(n) > 0)

        ! construct information about HCO_* corresponding pbuf location
        fldname_ns = 'HCO_' // trim(solsym(n))
        hco_pbuf_idx(n) = pbuf_get_index(fldname_ns, RC)
    enddo

    if(masterproc) then
        write(iulog,*) "hco_set_srf_emissions: first run pcnst_is_extfrc cache, extfrc_ndx:"
        do n = 1, gas_pcnst
           write(iulog,*) trim(solsym(n)), ' : ', pcnst_is_extfrc(n), pcnst_extfrc_ndx(n)
        end do
     endif

    ! Replicate functionality in extfrc_inti to create _XFRC... diagnostics
    call phys_getopts( &
         history_aerosol_out = history_aerosol, &
         history_chemistry_out = history_chemistry, &
         history_cesm_forcing_out = history_cesm_forcing )

    do n= 1,extcnt
       if (frc_from_dataset(n)) then
          spc_name = extfrc_lst(n)
          call addfld( trim(spc_name)//'_XFRC', (/ 'lev' /), 'A',  'molec/cm3/s', &
               'external forcing for '//trim(spc_name) )
          call addfld( trim(spc_name)//'_CLXF', horiz_only,  'A',  'molec/cm2/s', &
               'vertically intergrated external forcing for '//trim(spc_name) )
          call addfld( trim(spc_name)//'_CMXF', horiz_only,  'A',  'kg/m2/s', &
               'vertically intergrated external forcing for '//trim(spc_name) )
          if ( history_aerosol .or. history_chemistry ) then
             call add_default( trim(spc_name)//'_CLXF', 1, ' ' )
             call add_default( trim(spc_name)//'_CMXF', 1, ' ' )
          endif
          if ( history_cesm_forcing .and. spc_name == 'NO2' ) then
             call add_default( trim(spc_name)//'_CLXF', 1, ' ' )
             call add_default( trim(spc_name)//'_CMXF', 1, ' ' )
          endif
       endif
    enddo

    end subroutine hco_extfrc_inti
!EOC
end module hco_cc_emissions
