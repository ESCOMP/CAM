module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

use gw_utils, only: r8

implicit none
private
save

public :: BeresSourceDesc
public :: gw_beres_src

type :: BeresSourceDesc
   ! Whether wind speeds are shifted to be relative to storm cells.
   logical :: storm_shift
   ! Index for level where wind speed is used as the source speed.
   integer :: k
   ! Heating depths below this value [m] will be ignored.
   real(r8) :: min_hdepth
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh
   integer :: maxuh
   ! Heating depths [m].
   real(r8), allocatable :: hd(:)
   ! Table of source spectra.
   real(r8), allocatable :: mfcc(:,:,:)
end type BeresSourceDesc

contains

!==========================================================================

subroutine gw_beres_src(ncol, band, desc, u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth, maxq0)
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
  type(BeresSourceDesc), intent(in) :: desc

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Heating rate due to convection.
  real(r8), intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)

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

  ! Heating depth [m] and maximum heating in each column.
  real(r8), intent(out) :: hdepth(ncol), maxq0(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real(r8) :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real(r8) :: q0(ncol)

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
  real(r8) :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Heating rate conversion factor.
  real(r8), parameter :: CF = 20._r8
  ! Averaging length.
  real(r8), parameter :: AL = 1.0e5_r8

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------

  tau = 0.0_r8
  hdepth = 0.0_r8
  q0 = 0.0_r8
  tau0 = 0.0_r8

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source level, then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction.
  uconv = u(:,desc%k)
  vconv = v(:,desc%k)

  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv, yv, ubi(:,desc%k+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  boti = 0
  topi = 0
  do k = pver, 1, -1
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

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc%hd)

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

  !output max heating rate in K/day
  maxq0 = q0*24._r8*3600._r8

  ! Multipy by conversion factor
  q0 = q0 * CF

  if (desc%storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     CS = sign(max(abs(ubm(:,desc%k))-10._r8, 0._r8), ubm(:,desc%k))

     ! Average wind in heating region, relative to storm cells.
     uh = 0._r8
     do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
           uh = uh + ubm(:,k)/(boti-topi+1)
        end where
     end do

     uh = uh - CS

  else

     ! For shallow convection, wind is relative to ground, and "heating
     ! region" wind is just the source level wind.
     uh = ubm(:,desc%k)

  end if

  ! Limit uh to table range.
  uh = min(uh, real(desc%maxuh, r8))
  uh = max(uh, -real(desc%maxuh, r8))

  ! Speeds for critical level filtering.
  Umini =  band%ngwv
  Umaxi = -band%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band%dc))
     end where
  end do

  Umini = max(Umini, -band%ngwv)
  Umaxi = min(Umaxi, band%ngwv)

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else set
     ! tau0 = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0_r8
 
        tau(i,:,topi(i)+1) = tau0

     end if ! heating depth above min and not at the pole

  enddo

  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi

  ! Set phase speeds; just use reference speeds.
  c = spread(band%cref, 1, ncol)

end subroutine gw_beres_src

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

end module gw_convect
