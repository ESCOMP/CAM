module pbl_utils
!-----------------------------------------------------------------------!
! Module to hold PBL-related subprograms that may be used with multiple !
! different vertical diffusion schemes.                                 !
!                                                                       !
! Public subroutines:                                                   !
!
!     calc_obklen                                                       !
!                                                                       !
!------------------ History --------------------------------------------!
! Created: Apr. 2012, by S. Santos                                      !
!-----------------------------------------------------------------------!

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none
private

! Public Procedures
!----------------------------------------------------------------------!
! Excepting the initialization procedure, these are elemental
! procedures, so they can accept scalars or any dimension of array as
! arguments, as long as all arguments have the same number of
! elements.
public pbl_utils_init
public calc_ustar
public calc_obklen
public virtem
public compute_radf
public austausch_atm

real(r8), parameter :: ustar_min = 0.01_r8

real(r8) :: g         ! acceleration of gravity
real(r8) :: vk        ! Von Karman's constant
real(r8) :: cpair     ! specific heat of dry air
real(r8) :: rair      ! gas constant for dry air
real(r8) :: zvir      ! rh2o/rair - 1


!------------------------------------------------------------------------!
! Purpose: Compilers aren't creating optimized vector versions of        !
!          elemental routines, so we'll explicitly create them and bind  !
!          them via an interface for transparent use                     !
!------------------------------------------------------------------------!
interface calc_ustar
  module procedure calc_ustar_scalar
  module procedure calc_ustar_vector
end interface 

interface calc_obklen
  module procedure calc_obklen_scalar
  module procedure calc_obklen_vector
end interface

interface virtem
  module procedure virtem_vector1D
  module procedure virtem_vector2D  ! Used in hb_diff.F90
end interface



contains

subroutine pbl_utils_init(g_in,vk_in,cpair_in,rair_in,zvir_in)

  !-----------------------------------------------------------------------!
  ! Purpose: Set constants to be used in calls to later functions         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: g_in       ! acceleration of gravity
  real(r8), intent(in) :: vk_in      ! Von Karman's constant
  real(r8), intent(in) :: cpair_in   ! specific heat of dry air
  real(r8), intent(in) :: rair_in    ! gas constant for dry air
  real(r8), intent(in) :: zvir_in    ! rh2o/rair - 1

  g = g_in
  vk = vk_in
  cpair = cpair_in
  rair = rair_in
  zvir = zvir_in

end subroutine pbl_utils_init

subroutine calc_ustar_scalar( t,    pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!

  real(r8), intent(in) :: t         ! surface temperature
  real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

  real(r8), intent(out) :: rrho     ! 1./bottom level density
  real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

  rrho = rair * t / pmid
  ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )

end subroutine calc_ustar_scalar

subroutine calc_ustar_vector(n, t, pmid, taux, tauy, &
                                 rrho, ustar)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate ustar and bottom level density (necessary for      !
  !  Obukhov length calculation).                                         !
  !-----------------------------------------------------------------------!
  integer, intent(in) :: n             ! Length of vectors

  real(r8), intent(in) :: t(n)         ! surface temperature
  real(r8), intent(in) :: pmid(n)      ! midpoint pressure (bottom level)
  real(r8), intent(in) :: taux(n)      ! surface u stress [N/m2]
  real(r8), intent(in) :: tauy(n)      ! surface v stress [N/m2]


  real(r8), intent(out) :: rrho(n)     ! 1./bottom level density
  real(r8), intent(out) :: ustar(n)    ! surface friction velocity [m/s]


  rrho = rair * t / pmid
  ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )

end subroutine calc_ustar_vector

subroutine calc_obklen_scalar( ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!

  real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]

  real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen        ! Obukhov length

  ! Need kinematic fluxes for Obukhov:
  khfs = shflx*rrho/cpair
  kqfs = qflx*rrho
  kbfs = khfs + zvir*ths*kqfs

  ! Compute Obukhov length:
  obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

end subroutine calc_obklen_scalar

subroutine calc_obklen_vector(n, ths,  thvs, qflx, shflx, rrho, ustar, &
                                  khfs, kqfs, kbfs, obklen)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate Obukhov length and kinematic fluxes.               !
  !-----------------------------------------------------------------------!
  integer, intent(in) :: n                  ! Length of vectors

  real(r8), intent(in)  :: ths(n)           ! potential temperature at surface [K]
  real(r8), intent(in)  :: thvs(n)          ! virtual potential temperature at surface
  real(r8), intent(in)  :: qflx(n)          ! water vapor flux (kg/m2/s)
  real(r8), intent(in)  :: shflx(n)         ! surface heat flux (W/m2)

  real(r8), intent(in)  :: rrho(n)          ! 1./bottom level density [ m3/kg ]
  real(r8), intent(in)  :: ustar(n)         ! Surface friction velocity [ m/s ]

  real(r8), intent(out) :: khfs(n)          ! sfc kinematic heat flux [mK/s]
  real(r8), intent(out) :: kqfs(n)          ! sfc kinematic water vapor flux [m/s]
  real(r8), intent(out) :: kbfs(n)          ! sfc kinematic buoyancy flux [m^2/s^3]
  real(r8), intent(out) :: obklen(n)        ! Obukhov length


  ! Need kinematic fluxes for Obukhov:
  khfs = shflx*rrho/cpair
  kqfs = qflx*rrho
  kbfs = khfs + zvir*ths*kqfs

  ! Compute Obukhov length:
  obklen = -thvs * ustar**3 / (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

end subroutine calc_obklen_vector

subroutine virtem_vector1D(n, t,q, virtem)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  integer,  intent(in) :: n              ! vector length

  real(r8), intent(in) :: t(n), q(n)
  real(r8), intent(out):: virtem(n)

  virtem = t * (1.0_r8 + zvir*q)

end subroutine virtem_vector1D

subroutine virtem_vector2D(n, m, t, q, virtem)

  !-----------------------------------------------------------------------!
  ! Purpose: Calculate virtual temperature from temperature and specific  !
  !  humidity.                                                            !
  !-----------------------------------------------------------------------!

  integer,  intent(in) :: n, m            ! vector lengths

  real(r8), intent(in) :: t(n,m), q(n,m)
  real(r8), intent(out):: virtem(n,m)

  virtem = t * (1.0_r8 + zvir*q)

end subroutine virtem_vector2D


subroutine compute_radf( choice_radf, i, pcols, pver, ncvmax, ncvfin, ktop, qmin, &
                         ql, pi, qrlw, g, cldeff, zi, chs, lwp_CL, opt_depth_CL,  &
                         radinvfrac_CL, radf_CL )
  ! -------------------------------------------------------------------------- !
  ! Purpose:                                                                   !
  ! Calculate cloud-top radiative cooling contribution to buoyancy production. !
  ! Here,  'radf' [m2/s3] is additional buoyancy flux at the CL top interface  !
  ! associated with cloud-top LW cooling being mainly concentrated near the CL !
  ! top interface ( just below CL top interface ).  Contribution of SW heating !
  ! within the cloud is not included in this radiative buoyancy production     !
  ! since SW heating is more broadly distributed throughout the CL top layer.  !
  ! -------------------------------------------------------------------------- !

  !-----------------!
  ! Input variables !
  !-----------------!
  character(len=6), intent(in) :: choice_radf  ! Method for calculating radf
  integer,  intent(in)  :: i                   ! Index of current column
  integer,  intent(in)  :: pcols               ! Number of atmospheric columns
  integer,  intent(in)  :: pver                ! Number of atmospheric layers
  integer,  intent(in)  :: ncvmax              ! Max numbers of CLs (perhaps equal to pver)
  integer,  intent(in)  :: ncvfin(pcols)       ! Total number of CL in column
  integer,  intent(in)  :: ktop(pcols, ncvmax) ! ktop for current column
  real(r8), intent(in)  :: qmin                ! Minimum grid-mean LWC counted as clouds [kg/kg]
  real(r8), intent(in)  :: ql(pcols, pver)     ! Liquid water specific humidity [ kg/kg ]
  real(r8), intent(in)  :: pi(pcols, pver+1)   ! Interface pressures [ Pa ]
  real(r8), intent(in)  :: qrlw(pcols, pver)   ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
  real(r8), intent(in)  :: g                   ! Gravitational acceleration
  real(r8), intent(in)  :: cldeff(pcols,pver)  ! Effective Cloud Fraction [fraction]
  real(r8), intent(in)  :: zi(pcols, pver+1)   ! Interface heights [ m ]
  real(r8), intent(in)  :: chs(pcols, pver+1)  ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces.

  !------------------!
  ! Output variables !
  !------------------!
  real(r8), intent(out) :: lwp_CL(ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
  real(r8), intent(out) :: opt_depth_CL(ncvmax)   ! Optical depth of the CL top layer
  real(r8), intent(out) :: radinvfrac_CL(ncvmax)  ! Fraction of LW radiative cooling confined in the top portion of CL
  real(r8), intent(out) :: radf_CL(ncvmax)        ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]

  !-----------------!
  ! Local variables !
  !-----------------!
  integer :: kt, ncv
  real(r8) :: lwp, opt_depth, radinvfrac, radf


  !-----------------!
  ! Begin main code !
  !-----------------!
  lwp_CL        = 0._r8
  opt_depth_CL  = 0._r8
  radinvfrac_CL = 0._r8
  radf_CL       = 0._r8

  ! ---------------------------------------- !
  ! Perform do loop for individual CL regime !
  ! ---------------------------------------- !
  do ncv = 1, ncvfin(i)
    kt = ktop(i,ncv)
    !-----------------------------------------------------!
    ! Compute radf for each CL regime and for each column !
    !-----------------------------------------------------!
    if( choice_radf .eq. 'orig' ) then
      if( ql(i,kt) .gt. qmin .and. ql(i,kt-1) .lt. qmin ) then
        lwp       = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        ! Approximate LW cooling fraction concentrated at the inversion by using
        ! polynomial approx to exact formula 1-2/opt_depth+2/(exp(opt_depth)-1))

        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling = [ W/kg ]
        radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)
        ! We can disable cloud LW cooling contribution to turbulence by uncommenting:
        ! radf = 0._r8
      end if

    elseif( choice_radf .eq. 'ramp' ) then

      lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
      opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
      radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
      radinvfrac  = max(cldeff(i,kt)-cldeff(i,kt-1),0._r8) * radinvfrac
      radf        = qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) ! Cp*radiative cooling [W/kg]
      radf        = max( radinvfrac * radf * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 ) * chs(i,kt)

    elseif( choice_radf .eq. 'maxi' ) then

      ! Radiative flux divergence both in 'kt' and 'kt-1' layers are included
      ! 1. From 'kt' layer
        lwp         = ql(i,kt) * ( pi(i,kt+1) - pi(i,kt) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
        radf        = max( radinvfrac * qrlw(i,kt) / ( pi(i,kt) - pi(i,kt+1) ) * ( zi(i,kt) - zi(i,kt+1) ), 0._r8 )
      ! 2. From 'kt-1' layer and add the contribution from 'kt' layer
        lwp         = ql(i,kt-1) * ( pi(i,kt) - pi(i,kt-1) ) / g
        opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
        radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth) + opt_depth**2 )
        radf        = radf + max( radinvfrac * qrlw(i,kt-1) / ( pi(i,kt-1) - pi(i,kt) ) * ( zi(i,kt-1) - zi(i,kt) ), 0._r8 )
        radf        = max( radf, 0._r8 ) * chs(i,kt)

    endif

    lwp_CL(ncv)        = lwp
    opt_depth_CL(ncv)  = opt_depth
    radinvfrac_CL(ncv) = radinvfrac
    radf_CL(ncv)       = radf
  end do ! ncv = 1, ncvfin(i)
end subroutine compute_radf

subroutine austausch_atm(pcols, ncol, pver, ntop, nbot, ml2, ri, s2, kvf)

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Purpose: Computes exchange coefficients for free turbulent flows.     !
  !                                                                       !
  ! Method:                                                               !
  !                                                                       !
  ! The free atmosphere diffusivities are based on standard mixing length !
  ! forms for the neutral diffusivity multiplied by functns of Richardson !
  ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for    !
  ! momentum, potential temperature, and constitutents.                   !
  !                                                                       !
  ! The stable Richardson num function (Ri>0) is taken from Holtslag and  !
  ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))    !
  ! The unstable Richardson number function (Ri<0) is taken from  CCM1.   !
  ! f = sqrt(1 - 18*Ri)                                                   !
  !                                                                       !
  ! Author: B. Stevens (rewrite, August 2000)                             !
  !                                                                       !
  !---------------------------------------------------------------------- !

  ! --------------- !
  ! Input arguments !
  ! --------------- !

  integer,  intent(in)  :: pcols                ! Atmospheric columns dimension size
  integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
  integer,  intent(in)  :: pver                 ! Number of atmospheric layers
  integer,  intent(in)  :: ntop                 ! Top layer for calculation
  integer,  intent(in)  :: nbot                 ! Bottom layer for calculation

  real(r8), intent(in)  :: ml2(pver+1)          ! Mixing lengths squared
  real(r8), intent(in)  :: s2(pcols,pver)       ! Shear squared
  real(r8), intent(in)  :: ri(pcols,pver)       ! Richardson no

  ! ---------------- !
  ! Output arguments !
  ! ---------------- !

  real(r8), intent(out) :: kvf(pcols,pver+1)    ! Eddy diffusivity for heat and tracers

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  real(r8)              :: fofri                ! f(ri)
  real(r8)              :: kvn                  ! Neutral Kv

  integer               :: i                    ! Longitude index
  integer               :: k                    ! Vertical index

  real(r8), parameter :: zkmin =  0.01_r8       ! Minimum kneutral*f(ri).

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  kvf(:ncol,:)           = 0.0_r8

  ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm.

  do k = ntop, nbot - 1
     do i = 1, ncol
        if( ri(i,k) < 0.0_r8 ) then
           fofri = sqrt( max( 1._r8 - 18._r8 * ri(i,k), 0._r8 ) )
        else
           fofri = 1.0_r8 / ( 1.0_r8 + 10.0_r8 * ri(i,k) * ( 1.0_r8 + 8.0_r8 * ri(i,k) ) )
        end if
        kvn = ml2(k) * sqrt(s2(i,k))
        kvf(i,k+1) = max( zkmin, kvn * fofri )
     end do
  end do

end subroutine austausch_atm

end module pbl_utils
