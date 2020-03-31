module prim_si_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8

  implicit none
  private

  public :: geopotential_t
contains

!
!  The hydrostatic routine from CAM physics.
!  (FV stuff removed)
!  t,q input changed to take t_v
  !  removed gravit, so this routine returns PHI, not zm

  !xxx this subroutine is not used
subroutine geopotential_t(                                 &
       pmid   , pdel   ,  tv      , rair   ,  zm)

!-----------------------------------------------------------------------
!
! Purpose:
! Compute the geopotential height (above the surface) at the midpoints and
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------
    use dimensions_mod,     only : nlev, nlevp, np
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments



    real(r8), intent(in) :: pmid (np*np,nlev)    ! Midpoint pressures
    real(r8), intent(in) :: pdel (np*np,nlev)    ! layer thickness
    real(r8), intent(in) :: tv    (np*np,nlev)    ! temperature
    real(r8), intent(in) :: rair                 ! Gas constant for dry air
    ! real(r8), intent(in) :: gravit               ! Acceleration of gravity
    ! real(r8), intent(in) :: zvir                 ! rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zm(np*np,nlev)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
    integer :: ncol=np*np             ! Number of longitudes

    integer  :: i,k                ! Lon, level indices
    real(r8) :: hkk(np*np)         ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(np*np)         ! off-diagonal element
    real(r8) :: rog                ! Rair / gravit
    real(r8) :: zi(np*np,nlevp)     ! Height above surface at interfaces
!
!-----------------------------------------------------------------------
!
!    rog = rair/gravit
    rog = rair

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,nlevp) = 0.0_r8
    end do

! Compute zi, zm from bottom up.
! Note, zi(i,k) is the interface above zm(i,k)
    do k = nlev, 1, -1
! First set hydrostatic elements consistent with dynamics
       do i = 1,ncol
          hkl(i) = pdel(i,k) / pmid(i,k)
          hkk(i) = 0.5_r8 * hkl(i)
       end do

! Now compute tv, zm, zi
       do i = 1,ncol
          ! tvfac   = 1._r8 + zvir * q(i,k)
          ! tv      = t(i,k) * tvfac
          zm(i,k) = zi(i,k+1) + rog * tv(i,k) * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv(i,k) * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t
end module prim_si_mod
