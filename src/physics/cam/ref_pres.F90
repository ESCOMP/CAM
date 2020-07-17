module ref_pres
!--------------------------------------------------------------------------
! 
! Provides access to reference pressures for use by the physics
! parameterizations.  The pressures are provided by the dynamical core
! since it determines the grid used by the physics.
! 
! Note that the init method for this module is called before the init
! method in physpkg; therefore, most physics modules can use these
! reference pressures during their init phases.
! 
!--------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use ppgrid,         only: pver, pverp

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun

implicit none
public
save

! Reference pressures (Pa)
real(r8), protected :: pref_edge(pverp)     ! Layer edges
real(r8), protected :: pref_mid(pver)       ! Layer midpoints
real(r8), protected :: pref_mid_norm(pver)  ! Layer midpoints normalized by
                                            ! surface pressure ('eta' coordinate)

real(r8), protected :: ptop_ref             ! Top of model
real(r8), protected :: psurf_ref            ! Surface pressure

! Number of top levels using pure pressure representation
integer, protected :: num_pr_lev

! Pressure used to set troposphere cloud physics top (Pa)
real(r8), protected :: trop_cloud_top_press = 0._r8
! Top level for troposphere cloud physics
integer, protected :: trop_cloud_top_lev

! Pressure used to set MAM process top (Pa)
real(r8), protected :: clim_modal_aero_top_press = 0._r8
! Top level for MAM processes that impact climate
integer, protected :: clim_modal_aero_top_lev

! Molecular diffusion is calculated only if the model top is below this
! pressure (Pa).
real(r8), protected :: do_molec_press = 0.1_r8
! Pressure used to set bottom of molecular diffusion region (Pa).
real(r8), protected :: molec_diff_bot_press = 50._r8
! Flag for molecular diffusion, and molecular diffusion level index.
logical, protected :: do_molec_diff = .false.
integer, protected :: nbot_molec = 0

!====================================================================================
contains
!====================================================================================

subroutine ref_pres_readnl(nlfile)

   use spmd_utils,      only: masterproc
   use cam_abortutils,  only: endrun
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'ref_pres_readnl'

   namelist /ref_pres_nl/ trop_cloud_top_press, clim_modal_aero_top_press,&
        do_molec_press, molec_diff_bot_press
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'ref_pres_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, ref_pres_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! Check that top for modal aerosols is not lower than
      ! top for clouds.
      if (clim_modal_aero_top_press > trop_cloud_top_press) &
           call endrun("ERROR: clim_modal_aero_top press must be less &
           &than or equal to trop_cloud_top_press.")
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(trop_cloud_top_press,            1 , mpir8,   0, mpicom)
   call mpibcast(clim_modal_aero_top_press,       1 , mpir8,   0, mpicom)
   call mpibcast(do_molec_press,                  1 , mpir8,   0, mpicom)
   call mpibcast(molec_diff_bot_press,            1 , mpir8,   0, mpicom)
#endif

end subroutine ref_pres_readnl

!====================================================================================

subroutine ref_pres_init(pref_edge_in, pref_mid_in, num_pr_lev_in)

   ! Initialize reference pressures

   ! arguments
   real(r8), intent(in) :: pref_edge_in(:) ! reference pressure at layer edges (Pa)
   real(r8), intent(in) :: pref_mid_in(:)  ! reference pressure at layer midpoints (Pa)
   integer,  intent(in) :: num_pr_lev_in   ! number of top levels using pure pressure representation
   !---------------------------------------------------------------------------

   pref_edge = pref_edge_in
   pref_mid  = pref_mid_in
   num_pr_lev = num_pr_lev_in

   ptop_ref = pref_edge(1)
   psurf_ref = pref_edge(pverp)

   pref_mid_norm = pref_mid/psurf_ref

   ! Find level corresponding to the top of troposphere clouds.
   trop_cloud_top_lev = press_lim_idx(trop_cloud_top_press, &
      top=.true.)

   ! Find level corresponding to the top for MAM processes.
   clim_modal_aero_top_lev = press_lim_idx(clim_modal_aero_top_press, &
      top=.true.)

   ! Find level corresponding to the molecular diffusion bottom.
   do_molec_diff = (ptop_ref < do_molec_press)
   if (do_molec_diff) then
      nbot_molec = press_lim_idx(molec_diff_bot_press, &
         top=.false.)
   end if

end subroutine ref_pres_init

!====================================================================================

! Convert pressure limiters to the appropriate level.
pure function press_lim_idx(p, top) result(k_lim)
  ! Pressure
  real(r8), intent(in) :: p
  ! Is this a top or bottom limit?
  logical,  intent(in) :: top
  integer :: k_lim, k

  if (top) then
     k_lim = pver+1
     do k = 1, pver
        if (pref_mid(k) > p) then
           k_lim = k
           exit
        end if
     end do
  else
     k_lim = 0
     do k = pver, 1, -1
        if (pref_mid(k) < p) then
           k_lim = k
           exit
        end if
     end do
  end if

end function press_lim_idx

!====================================================================================

subroutine std_atm_pres(height, pstd)

   ! Use barometric formula for U.S. Standard Atmosphere to convert heights to pressures.
   ! This formula is valid up to 86 km.
   ! https://en.wikipedia.org/wiki/Barometric_formula

   ! arguments
   real(r8), intent(in)  :: height(:) ! height above sea level in meters
   real(r8), intent(out) :: pstd(:)   ! std pressure in Pa

   ! local vars
   integer, parameter  :: nreg = 7  ! number of regions
   real(r8), parameter :: hb(nreg) = & ! height a bottom of layer (m)
      (/0.0_r8, 1.1e4_r8, 2.0e4_r8, 3.2e4_r8, 4.7e4_r8, 5.1e4_r8, 7.1e4_r8/)
   real(r8), parameter :: pb(nreg) = & ! standard pressure (Pa)
      (/101325._r8, 22632.1_r8, 5474.89_r8, 868.02_r8, 110.91_r8, 66.94_r8, 3.96_r8/)
   real(r8), parameter :: tb(nreg) = & ! standard temperature (K)
      (/288.15_r8, 216.65_r8, 216.65_r8, 228.65_r8, 270.65_r8, 270.65_r8, 214.65_r8/)
   real(r8), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
      (/-0.0065_r8, 0.0_r8, 0.001_r8, 0.0028_r8, 0.0_r8, -0.0028_r8, -0.002_r8/)
   real(r8), parameter :: rg = 8.3144598_r8 ! universal gas constant (J/mol/K)
   real(r8), parameter :: g0 = 9.80665_r8   ! gravitational acceleration (m/s^2)
   real(r8), parameter :: mw = 0.0289644_r8 ! molar mass of dry air (kg/mol)
   real(r8), parameter :: c1 = g0*mw/rg

   integer :: i, ii, k, nlev
   logical :: found_region
   character(len=*), parameter :: routine = 'ref_pres::std_atm_pres'
   !---------------------------------------------------------------------------

   nlev = size(height)
   do k = 1, nlev

      ! find region containing height
      found_region = .false.
      find_region: do i = nreg, 1, -1
         if (height(k) >= hb(i)) then
            ii = i
            found_region = .true.
            exit find_region
         end if
      end do find_region

      if (.not. found_region) then
         write(iulog,*) routine, ': illegal height: ', height(k)
         call endrun(routine// ': illegal height < 0. ')
      end if

      if (lb(ii) /= 0._r8) then
         pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
      else
         pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
      end if

   end do

end subroutine std_atm_pres

!====================================================================================

end module ref_pres
