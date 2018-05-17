module eul_control_mod

! Eulerian dynamics shared data

use shr_kind_mod, only: r8=>shr_kind_r8
use pmgrid,       only: plat, plon, plev
use spmd_utils,   only: masterproc
use pspect,       only: pnmax

implicit none
private
save

real(r8) ,public ::  tmass(plat)  ! Mass integral for each latitude pair
real(r8) ,public ::  tmass0       ! Specified dry mass of atmosphere
real(r8) ,public ::  tmassf       ! Global mass integral
real(r8) ,public ::  qmassf       ! Global moisture integral
real(r8) ,public ::  fixmas       ! Proportionality factor for ps in dry mass fixer
real(r8) ,public ::  qmass1       ! Contribution to global moisture integral (mass
                                  !  weighting is based upon the "A" part of the hybrid grid)
real(r8) ,public ::  qmass2       ! Contribution to global moisture integral (mass
                                     !  weighting is based upon the "B" part of the hybrid grid)
real(r8) ,public ::  pdela(plon,plev)! pressure difference between interfaces (pressure
                                     !  defined using the "A" part of hybrid grid only)
real(r8) ,public ::  zgsint       ! global integral of geopotential height

integer  ,public :: pcray                   ! length of vector register (words) for FFT workspace
parameter (pcray=64)

real(r8) ,public :: trig (3*plon/2+1,plat)  ! trigonometric funct values used by fft
integer  ,public :: ifax(19,plat)           ! fft factorization of plon/2
real(r8), public :: cnfac                   ! Courant num factor(multiply by max |V|)
real(r8), public :: cnlim                   ! Maximum allowable courant number
real(r8), public :: hdfsd2(pnmax)   	       ! Del^2 mult. for each wave (vort-div)
real(r8), public :: hdfst2(pnmax)   	       ! Del^2 multiplier for each wave (t-q)
real(r8), public :: hdfsdn(pnmax)   	       ! Del^N mult. for each wave (vort-div)
real(r8), public :: hdfstn(pnmax)   	       ! Del^N multiplier for each wave (t-q)
real(r8), public :: hdiftq(pnmax,plev)      ! Temperature-tracer diffusion factors
real(r8), public :: hdifzd(pnmax,plev)      ! Vorticity-divergence diffusion factors
integer, parameter, public :: kmxhd2 = 2    ! Bottom level for increased del^2 diffusion
integer,  public :: nindex(plev)            ! Starting index for spectral truncation
integer,  public :: nmaxhd                  ! Maximum two dimensional wave number

! Variables set by namelist
real(r8), public :: dif2             ! del2 horizontal diffusion coeff.
integer,  public :: hdif_order       ! Order of horizontal diffusion operator
integer,  public :: kmnhdn           ! Nth order diffusion applied at and below layer kmnhdn.
                                     ! 2nd order diffusion is applied above layer kmnhdn.
real(r8), public :: hdif_coef        ! Nth order horizontal diffusion coefficient.
real(r8), public :: divdampn         ! Number of days (from nstep 0) to run divergence
real(r8), public :: eps              ! time filter coefficient. Defaults to 0.06.
integer,  public :: kmxhdc           ! number of levels (starting from model top) to apply Courant limiter.
integer,  public :: eul_nsplit       ! Intended number of dynamics timesteps per physics timestep
   
end module eul_control_mod
