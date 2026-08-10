!-----------------------------------------------------------------------
! CAM wrapper for mo_setsox (portable aqueous sulfur chemistry).
!-----------------------------------------------------------------------
module mo_setsox_cam

   use shr_kind_mod, only : r8 => shr_kind_r8
   use cam_logfile,  only : iulog
   use physics_types,only : physics_state
   use aerosol_state_mod, only: aerosol_state

   implicit none

   private
   public :: sox_inti, setsox
   public :: has_sox

   logical :: has_sox = .true.

contains
   subroutine sox_inti(aero_props)
      !-----------------------------------------------------------------------
      !        ... initialize the hetero sox routine
      !-----------------------------------------------------------------------

      use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
      use spmd_utils,   only : masterproc
      use phys_control, only : phys_getopts, cam_chempkg_is
      use carma_flags_mod, only : carma_do_cloudborne
      use mo_setsox,       only : setsox_init
      use sox_cldaero_mod, only : sox_cldaero_init
      use aerosol_properties_mod, only : aerosol_properties
      use mo_constants,   only : pi
      use cam_abortutils, only : endrun

      class(aerosol_properties), target, intent(in) :: aero_props

      logical :: modal_aerosols

      logical :: cloud_borne
      integer :: id_msa
      integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
      integer :: id_so4, id_h2so4
      logical :: inv_so2, inv_nh3, inv_hno3, inv_h2o2, inv_ho2, inv_o3

      ! Indices for species in the shared array of Henry's Law constant parameters
      integer :: heff_id_hno3, heff_id_so2, heff_id_nh3, heff_id_co2, heff_id_h2o2, heff_id_o3

      character(len=512) :: errmsg
      integer            :: errflg

      id_so4 = -1
      id_h2so4 = -1

      call phys_getopts( prog_modal_aero_out=modal_aerosols )
      cloud_borne = modal_aerosols .or. carma_do_cloudborne

      !-----------------------------------------------------------------
      !       ... get species indicies
      !-----------------------------------------------------------------

      if (cloud_borne) then
         id_h2so4 = get_spc_ndx( 'H2SO4' )
      else
         id_so4 = get_spc_ndx( 'SO4' )
      endif
      id_msa = get_spc_ndx( 'MSA' )

      inv_so2 = .false.
      id_so2 = get_inv_ndx( 'SO2' )
      inv_so2 = id_so2 > 0
      if ( .not. inv_so2 ) then
         id_so2 = get_spc_ndx( 'SO2' )
      endif

      inv_NH3 = .false.
      id_NH3 = get_inv_ndx( 'NH3' )
      inv_NH3 = id_NH3 > 0
      if ( .not. inv_NH3 ) then
         id_NH3 = get_spc_ndx( 'NH3' )
      endif

      inv_HNO3 = .false.
      id_HNO3 = get_inv_ndx( 'HNO3' )
      inv_HNO3 = id_hno3 > 0
      if ( .not. inv_HNO3 ) then
         id_HNO3 = get_spc_ndx( 'HNO3' )
      endif

      inv_H2O2 = .false.
      id_H2O2 = get_inv_ndx( 'H2O2' )
      inv_H2O2 = id_H2O2 > 0
      if ( .not. inv_H2O2 ) then
         id_H2O2 = get_spc_ndx( 'H2O2' )
      endif

      inv_HO2 = .false.
      id_HO2 = get_inv_ndx( 'HO2' )
      inv_HO2 = id_HO2 > 0
      if ( .not. inv_HO2 ) then
         id_HO2 = get_spc_ndx( 'HO2' )
      endif

      inv_o3 = get_inv_ndx( 'O3' ) > 0
      if (inv_o3) then
         id_o3 = get_inv_ndx( 'O3' )
      else
         id_o3 = get_spc_ndx( 'O3' )
      endif
      inv_ho2 = get_inv_ndx( 'HO2' ) > 0
      if (inv_ho2) then
         id_ho2 = get_inv_ndx( 'HO2' )
      else
         id_ho2 = get_spc_ndx( 'HO2' )
      endif

      has_sox = (id_so2>0) .and. (id_h2o2>0) .and. (id_o3>0) .and. (id_ho2>0)
      if (cloud_borne) then
         has_sox = has_sox .and. (id_h2so4>0)
      else
         has_sox = has_sox .and. (id_so4>0) .and. (id_nh3>0)
      endif

      ! Lookup Effective Henry's Law Constant parameters from the common
      ! data file read in the shared code.
      heff_id_hno3 = get_heff_index( 'HNO3' )
      heff_id_so2  = get_heff_index( 'SO2'  )
      heff_id_nh3  = get_heff_index( 'NH3'  )
      heff_id_co2  = get_heff_index( 'CO2'  )
      heff_id_h2o2 = get_heff_index( 'H2O2' )
      heff_id_o3   = get_heff_index( 'OX'   )

      has_sox = has_sox .and. (heff_id_hno3 > 0) .and. (heff_id_so2 > 0) &
                .and. (heff_id_nh3 > 0) .and. (heff_id_co2 > 0) &
                .and. (heff_id_h2o2 > 0) .and. (heff_id_o3 > 0)

      if (masterproc) then
         write(iulog,*) 'sox_inti: has_sox = ',has_sox
      endif

      if( has_sox ) then
         if (masterproc) then
            write(iulog,*) '-----------------------------------------'
            write(iulog,*) ' mo_setsox will do sox aerosols'
            write(iulog,*) '-----------------------------------------'
         endif
      else
         if (masterproc) then
            write(iulog,*) '-----------------------------------------'
            write(iulog,*) ' mo_setsox will not do sox aerosols'
            write(iulog,*) '-----------------------------------------'
         endif
         return
      end if

      ! call the portable init subroutines:
      call setsox_init( cloud_borne_in=cloud_borne, &
                       id_so2_in=id_so2,     inv_so2_in=inv_so2,   &
                       id_nh3_in=id_nh3,     inv_nh3_in=inv_nh3,   &
                       id_hno3_in=id_hno3,   inv_hno3_in=inv_hno3, &
                       id_h2o2_in=id_h2o2,   inv_h2o2_in=inv_h2o2, &
                       id_ho2_in=id_ho2,     inv_ho2_in=inv_ho2,   &
                       id_o3_in=id_o3,       inv_o3_in=inv_o3,     &
                       id_h2so4_in=id_h2so4, id_so4_in=id_so4,     id_msa_in=id_msa, &
                       heff_id_hno3_in=heff_id_hno3, heff_id_so2_in=heff_id_so2,   &
                       heff_id_nh3_in=heff_id_nh3,   heff_id_co2_in=heff_id_co2,   &
                       heff_id_h2o2_in=heff_id_h2o2, heff_id_o3_in=heff_id_o3 )

      call sox_cldaero_init(aero_props, &
                            id_msa_in=id_msa, id_h2so4_in=id_h2so4, id_so2_in=id_so2, &
                            id_h2o2_in=id_h2o2, id_nh3_in=id_nh3, pi_in=pi, &
                            ! sulfur oxidation is performed internally to GEOS-Chem, so the aerosol
                            ! and gas updates are not applied here (avoids double counting)
                            do_aqueous_sulfur_chemistry_aerosol_update_in=.not. cam_chempkg_is('geoschem_mam4'), &
                            errmsg=errmsg, errflg=errflg)
      if (errflg /= 0) then
         call endrun(trim(errmsg))
      end if

   end subroutine sox_inti

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine setsox( aero_state, state, &
                     pbuf,   &
                     ncol,   &
                     dtime,  &
                     press,  &
                     pdel,   &
                     tfld,   &
                     mbar,   &
                     lwc,    &
                     cldfrc, &
                     cldnum, &
                     invariants, &
                     qcw,    &
                     qin,    &
                     xphlwc, &
                     aqso4,  &
                     aqh2so4,&
                     aqso4_h2o2, &
                     aqso4_o3,   &
                     yph_in,  &
                     aqso4_h2o2_3d, &
                     aqso4_o3_3d &
                     )

      use physconst,    only : avogad, boltz, r_universal, mwco2, mwdry, gravit
      use ppgrid,       only : pver
      use shr_drydep_mod,  only : dheff
      use physics_buffer,  only : physics_buffer_desc
      use rad_constituents, only : rad_cnst_get_gas
      use mo_setsox,       only : setsox_sub
      use cam_abortutils,  only : endrun

      !-----------------------------------------------------------------------
      !      ... Dummy arguments
      !-----------------------------------------------------------------------
      class(aerosol_state), intent(in) :: aero_state
      type(physics_state),                intent(in)    :: state   ! Physics state variables
      type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:) ! Physics buffer
      integer,          intent(in)    :: ncol              ! num of columns in chunk
      real(r8),         intent(in)    :: dtime             ! time step (sec)
      real(r8),         intent(in)    :: press(:,:)        ! midpoint pressure ( Pa )
      real(r8),         intent(in)    :: pdel(:,:)         ! pressure thickness of levels (Pa)
      real(r8),         intent(in)    :: tfld(:,:)         ! temperature
      real(r8),         intent(in)    :: mbar(:,:)         ! mean wet atmospheric mass ( amu )
      real(r8), target, intent(in)    :: lwc(:,:)          ! cloud liquid water content (kg/kg)
      real(r8), target, intent(in)    :: cldfrc(:,:)       ! cloud fraction
      real(r8),         intent(in)    :: cldnum(:,:)       ! droplet number concentration (#/kg)
      real(r8),         intent(in)    :: invariants(:,:,:)
      real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
      real(r8),         intent(inout) :: qin(:,:,:)        ! transported species ( vmr )
      real(r8),         intent(out)   :: xphlwc(:,:)       ! pH value multiplied by lwc

      real(r8),         intent(out)   :: aqso4(:,:)                   ! aqueous phase chemistry
      real(r8),         intent(out)   :: aqh2so4(:,:)                 ! aqueous phase chemistry
      real(r8),         intent(out)   :: aqso4_h2o2(:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
      real(r8),         intent(out)   :: aqso4_o3(:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)
      real(r8),         intent(in), optional :: yph_in                ! ph value
      real(r8),         intent(out), optional :: aqso4_h2o2_3d(:, :)  ! 3D SO4 aqueous phase chemistry due to H2O2 (kg/m2)
      real(r8),         intent(out), optional :: aqso4_o3_3d(:, :)    ! 3D SO4 aqueous phase chemistry due to O3 (kg/m2)

      !-----------------------------------------------------------------------
      !      ... Local variables
      !-----------------------------------------------------------------------
      real(r8), pointer :: co2_mass_mixing_ratio(:,:) ! kg kg-1

      character(len=512) :: errmsg
      integer            :: errflg

      call rad_cnst_get_gas(0, 'CO2', state, pbuf, co2_mass_mixing_ratio)

      ! call the portable subroutine:
      call setsox_sub(
      aero_state = aero_state, &
         ncol       = ncol,       &
         pver       = pver,       &
         dtime      = dtime,      &
         press      = press,      &
         pdel       = pdel,       &
         tfld       = tfld,       &
         mbar       = mbar,       &
         lwc        = lwc,        &
         cldfrc     = cldfrc,     &
         cldnum     = cldnum,     &
         invariants = invariants, &
         co2_mass_mixing_ratio          = co2_mass_mixing_ratio, &
         dheff                          = dheff,       &
         AVOGADRO_KMOL                  = avogad,      &
         BOLTZMANN                      = boltz,       &
         GAS_CONSTANT_KMOL              = r_universal, &
         MOLECULAR_WEIGHT_CO2_G_MOL     = mwco2,       &
         MOLECULAR_WEIGHT_DRY_AIR_G_MOL = mwdry,       &
         gravit     = gravit,     &
         qcw        = qcw,        &
         qin        = qin,        &
         xphlwc     = xphlwc,     &
         aqso4      = aqso4,      &
         aqh2so4    = aqh2so4,    &
         aqso4_h2o2 = aqso4_h2o2, &
         aqso4_o3   = aqso4_o3,   &
         errmsg     = errmsg,     &
         errflg     = errflg,     &
         yph_in        = yph_in,        &
         aqso4_h2o2_3d = aqso4_h2o2_3d, &
         aqso4_o3_3d   = aqso4_o3_3d )

      if (errflg /= 0) then
         call endrun(trim(errmsg))
      end if

   end subroutine setsox

   !-----------------------------------------------------------------
   !       ... looks up Effective Henry's Law Constant parameters
   !-----------------------------------------------------------------
   pure integer function get_heff_index(species_name) result(index)
      use shr_drydep_mod, only: species_name_table

      character(len=*), intent(in) :: species_name

      do index = 1, size(species_name_table)
         if (trim(adjustl(species_name)) == &
             trim(adjustl(species_name_table(index)))) return
      end do
      index = -1
   end function get_heff_index

end module mo_setsox_cam
