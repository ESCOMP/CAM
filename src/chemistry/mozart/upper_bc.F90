
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

  logical :: apply_upper_bc = .false.

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
    use mo_tgcm_ubc, only: tgcm_ubc_inti
    use mo_snoe,     only: snoe_inti
    use mo_msis_ubc, only: msis_ubc_inti
    use constituents,only: cnst_get_ind

!---------------------------Local workspace-----------------------------
    logical :: zonal_avg
!-----------------------------------------------------------------------
    apply_upper_bc = ptop_ref<1._r8 ! Pa

    if (.not.apply_upper_bc) return

    call cnst_get_ind('F', f_ndx, abort=.false.)
    call cnst_get_ind('HF', hf_ndx, abort=.false.)

    zonal_avg     = .false.

!-----------------------------------------------------------------------
!       ... initialize the tgcm upper boundary module
!-----------------------------------------------------------------------
    call tgcm_ubc_inti( tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod)
    if (masterproc) write(iulog,*) 'ubc_init: after tgcm_ubc_inti'

!-----------------------------------------------------------------------
!       ... initialize the snoe module
!-----------------------------------------------------------------------
    call snoe_inti(snoe_ubc_file)
    if (masterproc) write(iulog,*) 'ubc_init: after snoe_inti'

!-----------------------------------------------------------------------
!       ... initialize the msis module
!-----------------------------------------------------------------------
    call msis_ubc_inti( zonal_avg )
    if (masterproc) write(iulog,*) 'ubc_init: after msis_ubc_inti'

  end subroutine ubc_init

!===============================================================================

  subroutine ubc_timestep_init(pbuf2d, state)
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

    use solar_parms_data, only: kp=>solar_parms_kp, ap=>solar_parms_ap, f107=>solar_parms_f107
    use solar_parms_data, only: f107a=>solar_parms_f107a, f107p=>solar_parms_f107p
    use mo_msis_ubc,      only: msis_timestep_init
    use mo_tgcm_ubc,      only: tgcm_timestep_init
    use mo_snoe,          only: snoe_timestep_init
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.apply_upper_bc) return

    call msis_timestep_init( ap, f107p, f107a )
    call tgcm_timestep_init( pbuf2d, state )
    call snoe_timestep_init( kp, f107 )

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, pint, zi, t, q, omega, phis, &
                           msis_temp, ubc_mmr, ubc_flux)

!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------
    use mo_msis_ubc,      only: get_msis_ubc
    use mo_snoe,          only: set_no_ubc, ndx_no
    use mo_tgcm_ubc,      only: set_tgcm_ubc
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

    call get_msis_ubc( lchnk, ncol, msis_temp, ubc_mmr )
    if( t_pert_ubc /= 0._r8 ) then
       msis_temp(:ncol) = msis_temp(:ncol) + t_pert_ubc
       if( any( msis_temp(:ncol) < 0._r8 ) ) then
          write(iulog,*) 'ubc_get_vals: msis temp < 0 after applying offset = ',t_pert_ubc
          call endrun
       end if
    end if
   
    !--------------------------------------------------------------------------------------------
    ! For WACCM-X, calculate upper boundary H flux
    !--------------------------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 

      call cnst_get_ind('H',  indx_H)
      qh_top => q(:,1,indx_H)

      do iCol = 1, ncol
        !--------------------------------------------------
        ! Get total density (rho) at top level
        !--------------------------------------------------
        nmbartop = 0.5_r8 * (pint(iCol,1) + pint(iCol,2)) / ( rairv(iCol,1,lchnk) * t(iCol,1) ) 

        !---------------------------------------------------------------------
        ! Calculate factor for Jean's escape flux once here, used twice below
        !---------------------------------------------------------------------
        zkt = (rEarth + ( 0.5_r8 * ( zi(iCol,1) + zi(iCol,2) ) + rga * phis(iCol) ) ) * &
                                                   cnst_mw(indx_H) / avogad * grav / ( kboltz * t(iCol,1) )
      
        ubc_flux(iCol,indx_H) = hfluxlimitfac * SQRT(kboltz/(2.0_r8 * pi * cnst_mw(indx_H) / avogad)) * &
                                qh_top(iCol) * nmbartop * &
                                SQRT(t(iCol,1)) * (1._r8 + zkt) * EXP(-zkt)
                                
        ubc_flux(iCol,indx_H) = ubc_flux(iCol,indx_H) * &
                                (2.03E-13_r8 * qh_top(iCol) * nmbartop / (cnst_mw(indx_H) / avogad) * t(iCol,1))
                                       
        !--------------------------------------------------------------------------------------------------------------
        !  Need to get helium number density (SI units) from mass mixing ratio.  mbarv is kg/mole, same as rMass units
        !  kg/kg * (kg/mole)/(kg/mole) * (Pa or N/m*m)/((Joules/K or N*m/K) * (K)) = m-3
        !--------------------------------------------------------------------------------------------------------------- 
!        nDensHETop  = qhe_top(iCol) * mbarv(iCol,1,lchnk) / cnst_mw(indx_HE) * &
!                                   0.5_r8 * (pint(iCol,1) + pint(iCol,2)) / (kboltz * t(iCol,1))
!       
!        !------------------------------------------------------------------------------------------------------
!        !  Get midpoint vertical velocity for top level by extrapolating from two levels below top (Pa/s)*P
!        !------------------------------------------------------------------------------------------------------                 
!
!        pScaleHeight = .5_r8*(rairv(iCol,2,lchnk)*t(iCol,1) + rairv(iCol,1,lchnk)*t(iCol,1)) / grav      
!        wN2 = -omega(iCol,2) / 0.5_r8 * (pint(iCol,1) + pint(iCol,2)) * pScaleHeight 
! 
!        pScaleHeight = .5_r8 * (rairv(iCol,3,lchnk)*t(iCol,2) + rairv(iCol,2,lchnk)*t(iCol,2)) / grav      
!        wN3 = -omega(iCol,3) / 0.5_r8 * (pint(iCol,1) + pint(iCol,2)) * pScaleHeight 
!  
!        !----------------------------------------------------
!        !  Get top midpoint level vertical velocity
!        !----------------------------------------------------
!        wNTop = 1.5_r8 * wN2 - 0.5_r8 * wN3 
!
!        !-----------------------------------------------------------------------------------------------------------------
!        ! Helium upper boundary flux is just helium density multiplied by vertical velocity (kg*/m3)*(m/s) = kg/s/m^2) 
!        !-----------------------------------------------------------------------------------------------------------------
!        ubc_flux(iCol,indx_HE) = -ndensHETop * wNTop 
!                         
      enddo

      ubc_mmr(:ncol,ndx_no) = 0.0_r8

    else ! for waccm

       rho_top(:ncol) = pint(:ncol,1) / (rairv(:ncol,1,lchnk)*msis_temp(:ncol))
       z_top(:ncol)   = m2km * zi(:ncol,1)

       call set_no_ubc  ( lchnk, ncol, z_top, ubc_mmr, rho_top )
       if( ndx_no > 0 .and. no_xfac_ubc /= 1._r8 ) then
          ubc_mmr(:ncol,ndx_no) = no_xfac_ubc * ubc_mmr(:ncol,ndx_no)
       end if

    endif

    call set_tgcm_ubc( lchnk, ncol, ubc_mmr, mbarv(:,1,lchnk))

    if (f_ndx .GT. 0) then
      ubc_mmr(:ncol, f_ndx) = 1.0e-15_r8
    endif
    if (hf_ndx .GT. 0) then
      ubc_mmr(:ncol, hf_ndx) = 1.0e-15_r8
    endif

    ! Zero out constituent ubc's that are not used.
    do m = 1, pcnst
       if (.not. cnst_fixed_ubc(m)) then
          ubc_mmr(:,m) = 0._r8
       end if
    end do

  end subroutine ubc_get_vals

end module upper_bc
