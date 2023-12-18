module gw_movmtn

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use gw_utils, only: r8

implicit none
private
save

public :: MovMtnSourceDesc
public :: gw_movmtn_src

type :: MovMtnSourceDesc
   ! Whether wind speeds are shifted to be relative to storm cells.
   logical :: storm_shift
   ! Index for level where wind speed is used as the source speed. ->700hPa
   integer :: k
   ! Heating depths below this value [m] will be ignored.
   real(r8) :: min_hdepth
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh !-bounds of the lookup table heating depths 
   integer :: maxuh ! bounds of the lookup table wind
   ! Heating depths [m].
   real(r8), allocatable :: hd(:), uh(:)
   ! Table of source spectra.
   real(r8), allocatable :: mfcc(:,:,:)  !is the lookup table f(depth, wind, phase speed)
end type MovMtnSourceDesc

contains

!==========================================================================

subroutine gw_movmtn_src(ncol, band, desc, u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, pver, qbo_hdepth_scaling

!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Wavelengths triggered by convection.
  type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  type(MovMtnSourceDesc), intent(in) :: desc

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Heating rate due to convection.
  real(r8), intent(in) :: netdt(:,:)  !from Zhang McFarlane
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1) !tau = momentum flux (m2/s2) at interface level ngwv = band of phase speeds
  ! Projectin of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol) !determined by vector direction of wind at 700hPa
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real(r8), intent(out) :: hdepth(ncol)    !calculated here in this code

!---------------------------Local Storage-------------------------------
  ! Column and (vertical) level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real(r8) :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real(r8) :: q0(ncol), qj(ncol)
  ! determined by vector direction of wind at 700hPa
  real(r8) :: xv700(ncol), yv700(ncol), u700(ncol) 
  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real(r8) :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real(r8) :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real(r8) :: CS(ncol),CS1(ncol)
  ! Wind speeds in wave direction
  real(r8) :: udiff(ncol),vdiff(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift
  ! Other wind quantities
  real(r8) :: ut(ncol),uc(ncol),umm(ncol)
  ! Tau from moving mountain lookup table
  real(r8) :: taumm(ncol)
  ! Heating rate conversion factor.  -> tuning factors
  real(r8), parameter :: CF = 20._r8  !(1/ (5%))  -> 5% of grid cell is covered with convection
  ! Averaging length.
  real(r8), parameter :: AL = 1.0e5_r8
  ! Index for moving mountain lookuptable
  real(r8) :: hdmm_idx(ncol), uhmm_idx(ncol)
  ! Index for ground based phase speed bin
  real(r8) :: c0(ncol,-band%ngwv:band%ngwv), c_idx(ncol,-band%ngwv:band%ngwv)
  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

  tau = 0.0_r8
  hdepth = 0.0_r8
  q0 = 0.0_r8
  tau0 = 0.0_r8

    !write(*,*) " DUMMY call ... of gw_movmtn_src "
    !RETURN

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source (steering level), then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction. (at 700hPa)
  uconv = u(:,desc%k)  !k defined in line21 (at specified altitude)
  vconv = v(:,desc%k)
  
  ! all GW calculations on a plane, which in our case is the wind at 700hPa  source level -> ubi is wind in this plane
  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv700, yv700, u700)
  !u700 = dot_2d(uconv, vconv, xv700,yv700)
  
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.  !nedt is heating profile from Zhang McFarlane (it's pressure coordinates, therefore k=0 is the top)
  boti = 0 !bottom
  topi = 0  !top
  do k = pver, 1, -1 !start at surface
     do i = 1, ncol
        if (boti(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0_r8) boti(i) = k
           end if
        else if (topi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000._r8) then
              topi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (.not. (netdt(i,k) > 0.0_r8)) topi(i) = k
           end if
        end if
     end do
     ! When all done, exit.
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.  (top-bottom altitudes)
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
!  hdepth = hdepth*qbo_hdepth_scaling
! where in the lookup table do I find this heating depth
  hd_idx = index_of_nearest(hdepth, desc%hd)
  !write(111) hd_idx
  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though. 
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0
 
  ! Maximum heating rate.
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:,k))
     end where
  end do

  ! Multipy by conversion factor (now 20* larger than what Zahng McFarlane said as they try to describe heating over 100km grid cell)
  q0 = q0 * CF
  qj = 9.81/285*q0 ! unit conversion to m/s3
  !write(117) qj
 
  ! Moving Mountain wind speeds
  ! Storm speed is taken to be the source wind speed.
  CS1 = max(sqrt(uconv**2._r8 + vconv**2._r8)-10._r8, 0._r8)
  CS = CS1*xv700 + CS1*yv700 
  ! calculate the xv and yv in the frame of reference of the wave 
  
  do i=1,ncol
     udiff(i) = u(i,topi(i)) - uconv(i)
     vdiff(i) = v(i,topi(i)) - vconv(i)
    ! xv(i) = -1._r8*udiff(i)/abs(udiff(i))
    ! yv(i) = -1._r8*vdiff(i)/abs(vdiff(i))
  end do
  call get_unit_vector(udiff, vdiff, xv, yv, ubi(:,desc%k))
  
  ! Project the local wind at midpoints onto the source wind. looping through altitudes ubm is a profile projected on to the steering level
  do k = 1, pver
     ubm(:,k) = -1._r8*dot_2d(u(:,k), v(:,k), xv, yv)
  end do
  
  ! Compute the interface wind projection by averaging the midpoint winds. (both same wind profile, just at different points of the grid)
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! determine wind for lookup table
  ! need wind speed at the top of the convecitve cell and at the steering level (700hPa)
  uh = 0._r8
  do i=1,ncol
     ut(i) = ubm(i,topi(i))
     uh(i) = ut(i) - CS(i) ! wind at top in the frame moving with the cell (at700hPa)
  end do

  ! Set phase speeds; just use reference speeds.
  c = spread(band%cref, 1, ncol)
  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else leave
     ! tau = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then
        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------
        !hdmm_idx = index_of_nearest(hdepth, desc%hd)
        !write(110) hdmm_idx
        uhmm_idx = index_of_nearest(uh, desc%uh)
        taumm(i) = abs(desc%mfcc(uhmm_idx(i),hd_idx(i),0))
        !write(1111) taumm(i)
        taumm(i) = taumm(i)*qj(i)*qj(i)/AL/1000._r8
        ! assign sign to MF based on the ground based phase speed, ground based phase speed = CS
        taumm(i) = -1._r8*sign(taumm(i),CS(i))
        !find the right phase speed bin
        c0(i,:) = CS(i) 
        c_idx(i,:) = index_of_nearest(c0(i,:),c(i,:))
        !write(1112) taumm(i)
        tau(i,c_idx(i,:),topi(i):topi(i)+1) = taumm(i) !input tau to top +1 level, interface level just below top of heating, remember it's in pressure - everything is upside down (source level of GWs, level where GWs are launched)
        
     end if ! heating depth above min and not at the pole
     
  enddo
  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi


end subroutine gw_movmtn_src

! Short routine to get the indices of a set of values rounded to their
! nearest points on a grid.
function index_of_nearest(x, grid) result(idx)
  real(r8), intent(in) :: x(:)
  real(r8), intent(in) :: grid(:)

  integer :: idx(size(x))

  real(r8) :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2._r8

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_movmtn
