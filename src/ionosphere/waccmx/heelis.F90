module heelis
  use shr_kind_mod  ,only: r8 => shr_kind_r8 ! 8-byte reals
  use edyn_maggrid  ,only: nmlon,nmlonp1,nmlat,ylonm,ylatm
  use edyn_geogrid  ,only: nlat
  use heelis_mod    ,only: heelis_update, heelis_flwv32
!
! phihm and pfrac are output of this module:
!
  use edyn_solve, only: phihm ! output high-latitude potential (nmlonp1,nmlat)

  implicit none
  save
  private

  public :: heelis_model
!
  real(r8), parameter :: h2deg = 15._r8   ! hour to degree
  integer, parameter  :: isouth = 1
  integer, parameter  :: inorth = 2

  contains
!-----------------------------------------------------------------------
  subroutine heelis_model(sunlons)
    use aurora_params, only: aurora_params_set

! Driver for Heelis empirical model to calculate high-latitude potential.
!
! Args:
    real(r8),intent(in) :: sunlons(nlat)  ! sun's location

!
! Set auroral parameters:
!
    call heelis_update()
    aurora_params_set = .true. ! this prevents unnecessary update in column physics

!
! Calculate  the heelis potential phihm in geomagnetic coordinates:
! (potm calls sub flwv32)
!
    call potm(sunlons)

  end subroutine heelis_model

!-----------------------------------------------------------------------
  subroutine potm(sunlons)
    use edyn_params, only: pi_dyn ! pi used in dynamo calculations
!
! Calculate heelis potential in geomagnetic coordinates.
!
! Args:
    real(r8),intent(in) :: sunlons(nlat)
!
! Local:
    integer :: j
    real(r8),dimension(nmlon) :: dlat,dlon,ratio
    integer,dimension(nmlon) :: iflag
!
    ratio(:) = 1._r8
    do j=1,nmlat
      iflag(:) = 1 ! must be updated at each j
      dlat(:) = ylatm(j)
      dlon(:) = ylonm(1:nmlon)-sunlons(1)
!
! flwv32 returns single-level Heelis potential in geomag coords:
!
      if (abs(ylatm(j)) > pi_dyn/6._r8) then
        call heelis_flwv32(dlat,dlon,ratio,pi_dyn,iflag,nmlon,phihm(:,j))
      else
        phihm(1:nmlon,j) = 0._r8
      endif
    enddo ! j=1,nmlat
!
! Periodic point:
    do j=1,nmlat
      phihm(nmlonp1,j) = phihm(1,j)
    enddo ! j=1,nmlat
  end subroutine potm
!-----------------------------------------------------------------------
end module heelis
