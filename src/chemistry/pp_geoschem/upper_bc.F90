
module upper_bc

!---------------------------------------------------------------------------------
! Module to compute the upper boundary condition for temperature (dry static energy)
! and trace gases. Uses the MSIS model, and SNOE and TIME GCM data.
!
! original code by Stacy Walters
! adapted by B. A. Boville
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod,only: grav   => shr_const_g,     &   ! gravitational constant (m/s^2)
                          kboltz => shr_const_boltz, &   ! Boltzmann constant
                          pi => shr_const_pi,        &   ! pi
                          rEarth => shr_const_rearth     ! Earth radius 
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use ref_pres,     only: ptop_ref

  implicit none
  private
  save
!
! Public interfaces
!
  public :: ubc_defaultopts    ! set default values of namelist variables
  public :: ubc_setopts        ! get namelist input
  public :: ubc_init           ! global initialization
  public :: ubc_timestep_init  ! time step initialization
  public :: ubc_get_vals       ! get ubc values for this step

! Namelist variables
  character(len=256) :: snoe_ubc_file = ' '
  real(r8)           :: t_pert_ubc  = 0._r8
  real(r8)           :: no_xfac_ubc = 1._r8

  character(len=256) :: tgcm_ubc_file = ' '
  integer            :: tgcm_ubc_cycle_yr = 0
  integer            :: tgcm_ubc_fixed_ymd = 0
  integer            :: tgcm_ubc_fixed_tod = 0
  integer            :: f_ndx, hf_ndx
  character(len=32)  :: tgcm_ubc_data_type = 'CYCLICAL'

  logical :: apply_upper_bc = .true.

!================================================================================================
contains
!================================================================================================

subroutine ubc_defaultopts(tgcm_ubc_file_out, tgcm_ubc_data_type_out, tgcm_ubc_cycle_yr_out, tgcm_ubc_fixed_ymd_out, &
     tgcm_ubc_fixed_tod_out, snoe_ubc_file_out, t_pert_ubc_out, no_xfac_ubc_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   real(r8), intent(out), optional         :: t_pert_ubc_out
   real(r8), intent(out), optional         :: no_xfac_ubc_out
   character(len=*), intent(out), optional :: tgcm_ubc_file_out
   character(len=*), intent(out), optional :: snoe_ubc_file_out
   integer         , intent(out), optional :: tgcm_ubc_cycle_yr_out
   integer         , intent(out), optional :: tgcm_ubc_fixed_ymd_out
   integer         , intent(out), optional :: tgcm_ubc_fixed_tod_out
   character(len=*), intent(out), optional :: tgcm_ubc_data_type_out

!-----------------------------------------------------------------------

   if ( present(tgcm_ubc_file_out) ) then
      tgcm_ubc_file_out = tgcm_ubc_file
   endif
   if ( present(tgcm_ubc_data_type_out) ) then
      tgcm_ubc_data_type_out = tgcm_ubc_data_type
   endif
   if ( present(tgcm_ubc_cycle_yr_out) ) then
      tgcm_ubc_cycle_yr_out = tgcm_ubc_cycle_yr
   endif
   if ( present(tgcm_ubc_fixed_ymd_out) ) then
      tgcm_ubc_fixed_ymd_out = tgcm_ubc_fixed_ymd
   endif
   if ( present(tgcm_ubc_fixed_tod_out) ) then
      tgcm_ubc_fixed_tod_out = tgcm_ubc_fixed_tod
   endif
   if ( present(snoe_ubc_file_out) ) then
      snoe_ubc_file_out = snoe_ubc_file
   endif
   if ( present(t_pert_ubc_out) ) then
      t_pert_ubc_out = t_pert_ubc
   endif
   if ( present(no_xfac_ubc_out) ) then
      no_xfac_ubc_out = no_xfac_ubc
   endif

end subroutine ubc_defaultopts

!================================================================================================

subroutine ubc_setopts(tgcm_ubc_file_in, tgcm_ubc_data_type_in, tgcm_ubc_cycle_yr_in, tgcm_ubc_fixed_ymd_in, &
     tgcm_ubc_fixed_tod_in, snoe_ubc_file_in, t_pert_ubc_in, no_xfac_ubc_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
!-----------------------------------------------------------------------

   use cam_abortutils, only : endrun

   real(r8), intent(in), optional         :: t_pert_ubc_in
   real(r8), intent(in), optional         :: no_xfac_ubc_in
   character(len=*), intent(in), optional :: tgcm_ubc_file_in
   character(len=*), intent(in), optional :: snoe_ubc_file_in
   integer         , intent(in), optional :: tgcm_ubc_cycle_yr_in
   integer         , intent(in), optional :: tgcm_ubc_fixed_ymd_in
   integer         , intent(in), optional :: tgcm_ubc_fixed_tod_in
   character(len=*), intent(in), optional :: tgcm_ubc_data_type_in

!-----------------------------------------------------------------------

   if ( present(tgcm_ubc_file_in) ) then
      tgcm_ubc_file = tgcm_ubc_file_in
   endif
   if ( present(tgcm_ubc_data_type_in) ) then
      tgcm_ubc_data_type = tgcm_ubc_data_type_in
   endif
   if ( present(tgcm_ubc_cycle_yr_in) ) then
      tgcm_ubc_cycle_yr = tgcm_ubc_cycle_yr_in
   endif
   if ( present(tgcm_ubc_fixed_ymd_in) ) then
      tgcm_ubc_fixed_ymd = tgcm_ubc_fixed_ymd_in
   endif
   if ( present(tgcm_ubc_fixed_tod_in) ) then
      tgcm_ubc_fixed_tod = tgcm_ubc_fixed_tod_in
   endif
   if ( present(snoe_ubc_file_in) ) then
      snoe_ubc_file = snoe_ubc_file_in
   endif
   if ( present(t_pert_ubc_in) ) then
      t_pert_ubc = t_pert_ubc_in
   endif
   if ( present(no_xfac_ubc_in) ) then
      no_xfac_ubc = no_xfac_ubc_in
      if( no_xfac_ubc < 0._r8 ) then
         write(iulog,*) 'ubc_setopts: no_xfac_ubc = ',no_xfac_ubc,' must be >= 0'
         call endrun
      end if
   endif

end subroutine ubc_setopts

!===============================================================================

  subroutine ubc_init()
!-----------------------------------------------------------------------
! Initialization of time independent fields for the upper boundary condition
! Calls initialization routine for MSIS, TGCM and SNOE
!-----------------------------------------------------------------------

    ! Assume we are running in a simulation with ptop >= 1 Pa
    apply_upper_bc = .false.

    if (.not.apply_upper_bc) return

  end subroutine ubc_init

!===============================================================================

  subroutine ubc_timestep_init(pbuf2d, state)
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

    use solar_parms_data, only: kp=>solar_parms_kp, ap=>solar_parms_ap, f107=>solar_parms_f107
    use solar_parms_data, only: f107a=>solar_parms_f107a, f107p=>solar_parms_f107p
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.apply_upper_bc) return

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, pint, zi, t, q, omega, phis, &
                           msis_temp, ubc_mmr, ubc_flux)

!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------
    use cam_abortutils,   only: endrun
    use physconst,        only: avogad, rairv, mbarv, rga ! Avogadro, gas constant, mean mass, universal gas constant
    use phys_control,     only: waccmx_is
    use constituents,     only: cnst_get_ind, cnst_mw, cnst_fixed_ubc  ! Needed for ubc_flux

!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc
    real(r8), intent(in)  :: t(pcols,pver)         ! midpoint temperature
    real(r8), intent(in),target :: q(pcols,pver,pcnst)   ! contituent mixing ratios (kg/kg)
    real(r8), intent(in)  :: omega(pcols,pver)     ! Vertical pressure velocity (Pa/s)
    real(r8), intent(in)  :: phis(pcols)           ! Surface geopotential (m2/s2)

    real(r8), intent(out) :: msis_temp(pcols)      ! upper bndy temperature (K)
    real(r8), intent(out) :: ubc_mmr(pcols,pcnst)  ! upper bndy mixing ratios (kg/kg)
    real(r8), intent(out) :: ubc_flux(pcols,pcnst) ! upper bndy flux (kg/s/m^2)

!---------------------------Local storage-------------------------------
    integer :: m                                   ! constituent index
    integer :: ierr                                ! error flag for allocates
    integer :: indx_H                              ! cnst index for H
    integer :: indx_HE                             ! cnst index for He
    integer :: iCol                                ! column loop counter

    real(r8), parameter :: m2km = 1.e-3_r8         ! meter to km
    real(r8) :: rho_top(pcols)                     ! density at top interface
    real(r8) :: z_top(pcols)                       ! height of top interface (km)

    real(r8), parameter :: hfluxlimitfac = 0.72_r8 ! Hydrogen upper boundary flux limiting factor

    real(r8) :: nmbartop                           ! Top level density (rho)
    real(r8) :: zkt                                ! Factor for H Jean's escape flux calculation
    real(r8) :: nDensHETop                         ! Helium number density (kg/m3)
    real(r8) :: pScaleHeight                       ! Scale height (m)
    real(r8) :: wN2                                ! Neutral vertical velocity second level (m/s)
    real(r8) :: wN3                                ! Neutral vertical velocity at third level (m/s)
    real(r8) :: wNTop                              ! Neutral vertical velocity at top level (m/s)

    real(r8), pointer :: qh_top(:)         ! Top level hydrogen mixing ratio (kg/kg)
!-----------------------------------------------------------------------

    ubc_mmr(:,:) = 0._r8
    ubc_flux(:,:) = 0._r8
    msis_temp(:) = 0._r8

    if (.not. apply_upper_bc) return

  end subroutine ubc_get_vals

end module upper_bc
