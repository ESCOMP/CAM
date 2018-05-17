module gw_oro

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use gw_utils, only: r8
use coords_1d, only: Coords1D

implicit none
private
save

! Public interface
public :: gw_oro_src

contains

!==========================================================================

subroutine gw_oro_src(ncol, band, p, &
     u, v, t, sgh, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_common, only: GWBand, pver, rair
  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(ncol)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)
  ! Midpoint Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real(r8) :: hdsp(ncol)
  ! Max orographic standard deviation to use.
  real(r8) :: sghmax
  ! c=0 stress from orography.
  real(r8) :: tauoro(ncol)
  ! Averages over source region.
  real(r8) :: nsrc(ncol) ! B-V frequency.
  real(r8) :: rsrc(ncol) ! Density.
  real(r8) :: usrc(ncol) ! Zonal wind.
  real(r8) :: vsrc(ncol) ! Meridional wind.

  ! Difference in interface pressure across source region.
  real(r8) :: dpsrc(ncol)

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves
  real(r8), parameter :: orohmin = 10._r8
  ! min wind speed for orographic waves
  real(r8), parameter :: orovmin = 2._r8

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp = 2.0_r8 * sgh

  k = pver
  src_level = k-1
  rsrc = p%mid(:,k)/(rair*t(:,k)) * p%del(:,k)
  usrc = u(:,k) * p%del(:,k)
  vsrc = v(:,k) * p%del(:,k)
  nsrc = nm(:,k)* p%del(:,k)

  do k = pver-1, 1, -1
     do i = 1, ncol
        if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
           src_level(i) = k-1
           rsrc(i) = rsrc(i) + &
                p%mid(i,k) / (rair*t(i,k)) * p%del(i,k)
           usrc(i) = usrc(i) + u(i,k) * p%del(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * p%del(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* p%del(i,k)
        end if
     end do
     ! Break the loop when all source levels found.
     if (all(src_level >= k)) exit
  end do

  do i = 1, ncol
     dpsrc(i) = p%ifc(i,pver+1) - p%ifc(i,src_level(i)+1)
  end do

  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc

  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,pver+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  do i = 1, ncol
     if ((ubi(i,pver+1) > orovmin) .and. (hdsp(i) > orohmin)) then
        sghmax = band%fcrit2 * (ubi(i,pver+1) / nsrc(i))**2
        tauoro(i) = 0.5_r8 * band%kwv * min(hdsp(i)**2, sghmax) * &
             rsrc(i) * nsrc(i) * ubi(i,pver+1)
     else
        tauoro(i) = 0._r8
        src_level(i) = pver
     end if
  end do

  ! Set the phase speeds and wave numbers in the direction of the source
  ! wind. Set the source stress magnitude (positive only, note that the
  ! sign of the stress is the same as (c-u).
  tau = 0._r8
  do k = pver, minval(src_level), -1
     where (src_level <= k) tau(:,0,k+1) = tauoro
  end do

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0._r8

end subroutine gw_oro_src

end module gw_oro
