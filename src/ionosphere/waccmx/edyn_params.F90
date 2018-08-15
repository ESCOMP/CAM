module edyn_params
!
! Constants for edynamo.
!
  use shr_kind_mod, only: r8 => shr_kind_r8            ! 8-byte reals
  use physconst, only: pi

  implicit none
  save

  private

  public :: pi, pi_dyn, re_dyn, r0, re, rtd, dtr, finit, h0, hs
  public :: kbotdyn, pbotdyn, cm2km

  real(r8),parameter :: &
    finit = 0._r8,      & ! initialization value
    re = 6.37122e8_r8,  & ! earth radius (cm)
    h0 = 9.7e6_r8,      & ! minimum height (cm) 
    r0 = re+h0,         & ! min height from earth center
    hs = 1.3e7_r8,      &
    cm2km = 1.e-5_r8      ! cm to km conversion
!
! Special pi for mag field calculations. If pi=4.*atan(1.) and code is
! linked with -lmass lib, then the last 2 digits (16th and 17th) of pi
! are different (56 instead of 12), resulting in theta0(j=49)==0., which
! is wrong (should be .1110e-15).
!   
  real(r8),parameter :: pi_dyn = 3.14159265358979312_r8 ! pi for dynamo
  real(r8),parameter :: re_dyn = 6.378165e8_r8          ! earth radius (cm) for dynamo
!
  real(r8),parameter :: dtr = pi/180._r8 ! degrees to radians
  real(r8),parameter :: rtd = 180._r8/pi ! radians to degrees
! 
! kbotdyn is the column index at which upward dynamo integrals begin. 
! This should correspond to about 85 km (zbotdyn). The index is determined
! by function find_kbotdyn (edynamo.F90) at every step (called by sub
! dynamo_input). The function insures that all processors use the same
! (minimum) kbotdyn.
!
  real(r8),parameter :: pbotdyn = 1.0_r8 ! Pa pressure (~80 km) at which to set kbotdyn
  integer :: kbotdyn = -1                     

end module edyn_params
