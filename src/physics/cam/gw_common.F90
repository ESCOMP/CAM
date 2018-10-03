module gw_common

!
! This module contains code common to different gravity wave
! parameterizations.
!
use gw_utils, only: r8
use coords_1d, only: Coords1D


implicit none
private
save

! Public interface.

public :: GWBand

public :: gw_common_init
public :: gw_prof
public :: gw_drag_prof
public :: qbo_hdepth_scaling
public :: calc_taucd, momentum_flux, momentum_fixer
public :: energy_change, energy_fixer
public :: coriolis_speed, adjust_inertial

public :: pver
public :: west, east, north, south
public :: pi
public :: gravit
public :: rair

! Number of levels in the atmosphere.
integer, protected :: pver = 0

! Whether or not to enforce an upper boundary condition of tau = 0.
logical :: tau_0_ubc = .false.

! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4

! Scaling factor for generating QBO
real(r8), protected :: qbo_hdepth_scaling

! 3.14159...
real(r8), parameter :: pi = acos(-1._r8)

! Acceleration due to gravity.
real(r8), protected :: gravit = huge(1._r8)

! Gas constant for dry air.
real(r8), protected :: rair = huge(1._r8)

!
! Private variables
!

! Interface levels for gravity wave sources.
integer :: ktop = huge(1)

! Background diffusivity.
real(r8), parameter :: dback = 0.05_r8

! rair/gravit
real(r8) :: rog = huge(1._r8)

! Newtonian cooling coefficients.
real(r8), allocatable :: alpha(:)

! Inverse Prandtl number.
real(r8) :: prndl

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(r8), parameter :: taumin = 1.e-10_r8
! Maximum wind tendency from stress divergence (before efficiency applied).
! 400 m/s/day
real(r8), parameter :: tndmax = 400._r8 / 86400._r8
! Maximum allowed change in u-c (before efficiency applied).
real(r8), parameter :: umcfac = 0.5_r8
! Minimum value of (u-c)**2.
real(r8), parameter :: ubmc2mn = 0.01_r8

! Type describing a band of wavelengths into which gravity waves can be
! emitted.
! Currently this has to have uniform spacing (i.e. adjacent elements of
! cref are exactly dc apart).
type :: GWBand
   ! Dimension of the spectrum.
   integer :: ngwv
   ! Delta between nearest phase speeds [m/s].
   real(r8) :: dc
   ! Reference speeds [m/s].
   real(r8), allocatable :: cref(:)
   ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
   real(r8) :: fcrit2
   ! Horizontal wave number [1/m].
   real(r8) :: kwv
   ! Effective horizontal wave number [1/m] (fcrit2*kwv).
   real(r8) :: effkwv
end type GWBand

interface GWBand
   module procedure new_GWBand
end interface

contains

!==========================================================================

! Constructor for a GWBand that calculates derived components.
function new_GWBand(ngwv, dc, fcrit2, wavelength) result(band)
  ! Used directly to set the type's components.
  integer, intent(in) :: ngwv
  real(r8), intent(in) :: dc
  real(r8), intent(in) :: fcrit2
  ! Wavelength in meters.
  real(r8), intent(in) :: wavelength

  ! Output.
  type(GWBand) :: band

  ! Wavenumber index.
  integer :: l

  ! Simple assignments.
  band%ngwv = ngwv
  band%dc = dc
  band%fcrit2 = fcrit2

  ! Uniform phase speed reference grid.
  allocate(band%cref(-ngwv:ngwv))
  band%cref = [( dc * l, l = -ngwv, ngwv )]

  ! Wavenumber and effective wavenumber come from the wavelength.
  band%kwv = 2._r8*pi / wavelength
  band%effkwv = band%fcrit2 * band%kwv

end function new_GWBand

!==========================================================================

subroutine gw_common_init(pver_in, &
     tau_0_ubc_in, ktop_in, gravit_in, rair_in, alpha_in, & 
     prndl_in, qbo_hdepth_scaling_in, errstring)

  integer,  intent(in) :: pver_in
  logical,  intent(in) :: tau_0_ubc_in
  integer,  intent(in) :: ktop_in
  real(r8), intent(in) :: gravit_in
  real(r8), intent(in) :: rair_in
  real(r8), intent(in) :: alpha_in(:)
  real(r8), intent(in) :: prndl_in
  real(r8), intent(in) :: qbo_hdepth_scaling_in
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  integer :: ierr

  errstring = ""

  pver = pver_in
  tau_0_ubc = tau_0_ubc_in
  ktop = ktop_in
  gravit = gravit_in
  rair = rair_in
  allocate(alpha(pver+1), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return
  alpha = alpha_in
  prndl = prndl_in
  qbo_hdepth_scaling = qbo_hdepth_scaling_in

  rog = rair/gravit

end subroutine gw_common_init

!==========================================================================

subroutine gw_prof (ncol, p, cpair, t, rhoi, nm, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  use gw_utils, only: midpoint_interp
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Specific heat of dry air, constant pressure.
  real(r8), intent(in) :: cpair
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)

  ! Interface density.
  real(r8), intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real(r8) :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real(r8) :: n2

  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real(r8), parameter :: n2min = 5.e-5_r8

  !------------------------------------------------------------------------
  ! Determine the interface densities and Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------

  ! The top interface values are calculated assuming an isothermal
  ! atmosphere above the top level.
  k = 1
  do i = 1, ncol
     ti(i,k) = t(i,k)
     rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
     ni(i,k) = sqrt(gravit*gravit / (cpair*ti(i,k)))
  end do

  ! Interior points use centered differences.
  ti(:,2:pver) = midpoint_interp(t)
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) * p%rdst(i,k-1)
        n2 = gravit*gravit/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do

  ! Bottom interface uses bottom level temperature, density; next interface
  ! B-V frequency.
  k = pver+1
  do i = 1, ncol
     ti(i,k) = t(i,k-1)
     rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
     ni(i,k) = ni(i,k-1)
  end do

  !------------------------------------------------------------------------
  ! Determine the midpoint Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------
  nm = midpoint_interp(ni)

end subroutine gw_prof

!==========================================================================

subroutine gw_drag_prof(ncol, band, p, src_level, tend_level, dt, &
     t, vramp,   &
     piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   &
     effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, &
     ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, ro_adjust, &
     kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff )

  !-----------------------------------------------------------------------
  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !-----------------------------------------------------------------------

  use gw_diffusion, only: gw_ediff, gw_diff_tend
  use linear_1d_operators, only: TriDiagDecomp

  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real(r8), intent(in) :: dt

  ! Midpoint and interface temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Log of interface pressures.
  real(r8), intent(in) :: piln(ncol,pver+1)
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(r8), intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Molecular thermal diffusivity.
  real(r8), intent(in) :: kvtt(ncol,pver+1)
  ! Constituent array.
  real(r8), intent(in) :: q(:,:,:)
  ! Dry static energy.
  real(r8), intent(in) :: dse(ncol,pver)
  ! Coefficient to ramp down diffusion coeff.
  real(r8), pointer, intent(in) :: vramp(:)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Zonal/meridional wind tendencies.
  real(r8), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(r8), intent(out) :: ttgw(ncol,pver)
  ! Gravity wave constituent tendency.
  real(r8), intent(out) :: qtgw(:,:,:)

  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pver+1)

  ! Gravity wave wind tendency for each wave.
  real(r8), intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8), intent(out) :: dttdf(ncol,pver)
  real(r8), intent(out) :: dttke(ncol,pver)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  ! Diagnosed horizontal wavenumber for ridges.
  real(r8), intent(in), optional :: &
       kwvrdg(ncol)

  ! Factor for saturation calculation. Here backwards 
  ! compatibility. I believe it should be 1.0 (jtb). 
  ! Looks like it has been 2.0 for a while in CAM.
  real(r8), intent(in), optional :: &
       satfac_in

  logical, intent(in), optional :: lapply_effgw_in, lapply_vdiff

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, m, i

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! "Total" and saturation diffusivity.
  real(r8) :: d(ncol)
  ! Imaginary part of vertical wavenumber.
  real(r8) :: mi(ncol)
  ! Stress after damping.
  real(r8) :: taudmp(ncol)
  ! Saturation stress.
  real(r8) :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real(r8) :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real(r8) :: ubt(ncol,pver), ubtl(ncol)
  real(r8) :: wrk(ncol)
  ! Ratio used for ubt tndmax limiting.
  real(r8) :: ubt_lim_ratio(ncol)

  ! saturation factor. Defaults to 2.0
  ! unless overidden by satfac_in
  real(r8) :: satfac

  logical :: lapply_effgw,do_vertical_diffusion

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp

  !------------------------------------------------------------------------

  if (present(satfac_in)) then
     satfac = satfac_in
  else
     satfac = 2._r8
  endif

  ! Default behavior is to apply vertical diffusion.
  ! The user has the option to turn off vert diffusion
  do_vertical_diffusion = .true.
  if (present(lapply_vdiff)) then
     do_vertical_diffusion = lapply_vdiff
  endif

  ! Default behavior is to apply effgw and
  ! tendency limiters as designed by Sean
  ! Santos (lapply_effgw=.TRUE.). However,
  ! WACCM non-oro GW need to be retuned before
  ! this can done to them. --jtb 03/02/16
  if (present(lapply_effgw_in)) then
      lapply_effgw = lapply_effgw_in
  else
      lapply_effgw = .TRUE.
  endif

  
  ! Lowest levels that loops need to iterate over.
  kbot_tend = maxval(tend_level)
  kbot_src = maxval(src_level)

  ! Initialize gravity wave drag tendencies to zero.

  utgw = 0._r8
  vtgw = 0._r8

  gwut = 0._r8

  dttke = 0._r8
  ttgw = 0._r8

  dttdf = 0._r8
  qtgw = 0._r8

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi = 0._r8
  taudmp = 0._r8
  tausat = 0._r8
  ubmc = 0._r8
  ubmc2 = 0._r8
  wrk = 0._r8

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
  ! do k = kbot_src-1, ktop, -1 !++jtb I think this is right 
  do k = kbot_src, ktop, -1  !++ but this is in model now 
     
     ! Determine the diffusivity for each column.

     d = dback + kvtt(:,k)

     do l = -band%ngwv, band%ngwv

        ! Determine the absolute value of the saturation stress.
        ! Define critical levels where the sign of (u-c) changes between
        ! interfaces.
        ubmc = ubi(:,k) - c(:,l)

        tausat = 0.0_r8

        if (present(kwvrdg)) then
           where (src_level >= k)
              ! Test to see if u-c has the same sign here as the level below.
              where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
                 tausat = abs(  kwvrdg  * rhoi(:,k) * ubmc**3 / &
                    (satfac*ni(:,k)))
              end where
           end where
        else
           where (src_level >= k)
              ! Test to see if u-c has the same sign here as the level below.
              where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
                 tausat = abs(band%effkwv * rhoi(:,k) * ubmc**3 / &
                    (satfac*ni(:,k)))
              end where
           end where
        end if

        if (present(ro_adjust)) then
           where (src_level >= k)
              tausat = tausat * sqrt(ro_adjust(:,l,k))
           end where
        end if

        if (present(kwvrdg)) then
           where (src_level >= k)
              ! Compute stress for each wave. The stress at this level is the
              ! min of the saturation stress and the stress at the level below
              ! reduced by damping. The sign of the stress must be the same as
              ! at the level below.

              ubmc2 = max(ubmc**2, ubmc2mn)
              mi = ni(:,k) / (2._r8 *   kwvrdg * ubmc2) * &  ! Is this 2._r8 related to satfac?
                 (alpha(k) + ni(:,k)**2/ubmc2 * d)
              wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))

              taudmp = tau(:,l,k+1)

              ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
              ! we limit tau, so instead limit the arrays used to set it.
              where (tausat <= taumin) tausat = 0._r8
              where (taudmp <= taumin) taudmp = 0._r8

              tau(:,l,k) = min(taudmp, tausat)
           end where

        else

           where (src_level >= k)

              ! Compute stress for each wave. The stress at this level is the
              ! min of the saturation stress and the stress at the level below
              ! reduced by damping. The sign of the stress must be the same as
              ! at the level below.

              ubmc2 = max(ubmc**2, ubmc2mn)
              mi = ni(:,k) / (2._r8 * band%kwv * ubmc2) * &
                 (alpha(k) + ni(:,k)**2/ubmc2 * d)
              wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))

              taudmp = tau(:,l,k+1) * exp(wrk)

              ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
              ! we limit tau, so instead limit the arrays used to set it.
              where (tausat <= taumin) tausat = 0._r8
              where (taudmp <= taumin) taudmp = 0._r8

              tau(:,l,k) = min(taudmp, tausat)
           end where
        endif

     end do
  end do

  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0._r8

  ! Apply efficiency to completed stress profile.
  if (lapply_effgw) then
     do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           where (k-1 <= tend_level)
              tau(:,l,k) = tau(:,l,k) * effgw
           end where
        end do
     end do
  end if

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------

  ! Loop over levels from top to bottom
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_r8

     do l = -band%ngwv, band%ngwv    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit * (tau(:,l,k+1)-tau(:,l,k)) * p%rdel(:,k)

        ! Apply first tendency limit to maintain numerical stability.
        ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! The limiter is somewhat stricter, so that we don't come anywhere
        ! near reversing c-u.
        ubtl = min(ubtl, umcfac * abs(c(:,l)-ubm(:,k)) / dt)

        if (.not. lapply_effgw) ubtl = min(ubtl, tndmax)
        
        where (k <= tend_level)

           ! Save tendency for each wave (for later computation of kzz).
           ! sign function returns magnitude of ubtl with sign of c-ubm 
           ! Renders ubt/ubm check for mountain waves unecessary
           gwut(:,k,l) = sign(ubtl, c(:,l)-ubm(:,k))
           ubt(:,k) = ubt(:,k) + gwut(:,k,l)

        end where

     end do

     if (lapply_effgw) then
        ! Apply second tendency limit to maintain numerical stability.
        ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
        ! permitted.
        ! This can only happen above tend_level, so don't bother checking the
        ! level explicitly.
        where (abs(ubt(:,k)) > tndmax)
           ubt_lim_ratio = tndmax/abs(ubt(:,k))
           ubt(:,k) = ubt_lim_ratio * ubt(:,k)
        elsewhere
           ubt_lim_ratio = 1._r8
        end where
     else
        ubt_lim_ratio = 1._r8
     end if
     
     do l = -band%ngwv, band%ngwv
        gwut(:,k,l) = ubt_lim_ratio*gwut(:,k,l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        where (k <= tend_level)
           tau(:,l,k+1) = tau(:,l,k) + & 
                abs(gwut(:,k,l)) * p%del(:,k) / gravit 
        end where
     end do

     ! Project the mean wind tendency onto the components.
     where (k <= tend_level)
        utgw(:,k) = ubt(:,k) * xv
        vtgw(:,k) = ubt(:,k) * yv
     end where

     if (associated(vramp)) then
        utgw(:,k) = utgw(:,k) * vramp(k)
        vtgw(:,k) = vtgw(:,k) * vramp(k)
     endif

     ! End of level loop.
  end do


  ! Block to undo Sean Santos mods to effgw and limiters.
  ! Here because non-oro GW in WACCM need extensive re-tuning
  ! before Sean's mods can be adopted. --jtb 03/02/16
  !==========================================
  if (.not.(lapply_effgw)) then
     do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           where (k-1 <= tend_level)
              tau(:,l,k) = tau(:,l,k) * effgw
           end where
        end do
     end do
     do k = ktop, kbot_tend
        do l = -band%ngwv, band%ngwv
           gwut(:,k,l) = gwut(:,k,l) * effgw
        end do
        utgw(:,k) = utgw(:,k) * effgw
        vtgw(:,k) = vtgw(:,k) * effgw
     end do
  end if
  !===========================================

  if (do_vertical_diffusion) then

     ! Calculate effective diffusivity and LU decomposition for the
     ! vertical diffusion solver.
     call gw_ediff (ncol, pver, band%ngwv, kbot_tend, ktop, tend_level, &
          gwut, ubm, nm, rhoi, dt, prndl, gravit, p, c, vramp, &
          egwdffi, decomp, ro_adjust=ro_adjust)

     ! Calculate tendency on each constituent.
     do m = 1, size(q,3)

        call gw_diff_tend(ncol, pver, kbot_tend, ktop, q(:,:,m), &
             dt, decomp, qtgw(:,:,m))

     enddo

     ! Calculate tendency from diffusing dry static energy (dttdf).
     call gw_diff_tend(ncol, pver, kbot_tend, ktop, dse, dt, decomp, dttdf)

  endif

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
  do l = -band%ngwv, band%ngwv
     do k = ktop, kbot_tend
        dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke + dttdf

  if (associated(vramp)) then
     do k = ktop, kbot_tend
        ttgw(:,k) = ttgw(:,k) * vramp(k)
     enddo
  endif

  ! Deallocate decomp.
  call decomp%finalize()

end subroutine gw_drag_prof

!==========================================================================

! Calculate Reynolds stress for waves propagating in each cardinal
! direction.

function calc_taucd(ncol, ngwv, tend_level, tau, c, xv, yv, ubi) &
     result(taucd)

  ! Column and gravity wave wavenumber dimensions.
  integer, intent(in) :: ncol, ngwv
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(:)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(:,-ngwv:,:)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(:,-ngwv:)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(:), yv(:)
  ! Projection of wind at interfaces.
  real(r8), intent(in) :: ubi(:,:)

  real(r8) :: taucd(ncol,pver+1,4)

  ! Indices.
  integer :: i, k, l

  ! ubi at tend_level.
  real(r8) :: ubi_tend(ncol)

  ! Signed wave Reynolds stress.
  real(r8) :: tausg(ncol)

  ! Reynolds stress for waves propagating behind and forward of the wind.
  real(r8) :: taub(ncol)
  real(r8) :: tauf(ncol)

  taucd = 0._r8
  tausg = 0._r8

  ubi_tend = (/ (ubi(i,tend_level(i)+1), i = 1, ncol) /)

  do k = ktop, maxval(tend_level)+1

     taub = 0._r8
     tauf = 0._r8

     do l = -ngwv, ngwv
        where (k-1 <= tend_level)

           tausg = sign(tau(:,l,k), c(:,l)-ubi(:,k))

           where ( c(:,l) < ubi_tend )
              taub = taub + tausg
           elsewhere
              tauf = tauf + tausg
           end where

        end where
     end do

     where (k-1 <= tend_level)
        where (xv > 0._r8)
           taucd(:,k,east) = tauf * xv
           taucd(:,k,west) = taub * xv
        elsewhere
           taucd(:,k,east) = taub * xv
           taucd(:,k,west) = tauf * xv
        end where

        where ( yv > 0._r8)
           taucd(:,k,north) = tauf * yv
           taucd(:,k,south) = taub * yv
        elsewhere
           taucd(:,k,north) = taub * yv
           taucd(:,k,south) = tauf * yv
        end where
     end where

  end do

end function calc_taucd

!==========================================================================

! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real(r8), intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real(r8), intent(out) :: um_flux(:), vm_flux(:)

  integer :: i

  ! Tendency for U & V below source level.
  do i = 1, size(tend_level)
     um_flux(i) = taucd(i,tend_level(i)+1, east) + &
                  taucd(i,tend_level(i)+1, west)
     vm_flux(i) = taucd(i,tend_level(i)+1,north) + &
                  taucd(i,tend_level(i)+1,south)
  end do

end subroutine momentum_flux

!==========================================================================

! Subtracts a change in momentum in the gravity wave levels from wind
! tendencies in lower levels, ensuring momentum conservation.
subroutine momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Components of momentum change sourced from the bottom.
  real(r8), intent(in) :: um_flux(:), vm_flux(:)
  ! Wind tendencies.
  real(r8), intent(inout) :: utgw(:,:), vtgw(:,:)

  ! Indices.
  integer :: i, k
  ! Reciprocal of total mass.
  real(r8) :: rdm(size(tend_level))
  ! Average changes in velocity from momentum change being spread over
  ! total mass.
  real(r8) :: du(size(tend_level)), dv(size(tend_level))

  ! Total mass from ground to source level: rho*dz = dp/gravit
  do i = 1, size(tend_level)
     rdm(i) = gravit/(p%ifc(i,pver+1)-p%ifc(i,tend_level(i)+1))
  end do

  ! Average velocity changes.
  du = -um_flux*rdm
  dv = -vm_flux*rdm

  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        utgw(:,k) = utgw(:,k) + du
        vtgw(:,k) = vtgw(:,k) + dv
     end where
  end do
  
end subroutine momentum_fixer

!==========================================================================

! Calculate the change in total energy from tendencies up to this point.
subroutine energy_change(dt, p, u, v, dudt, dvdt, dsdt, de)

  ! Time step.
  real(r8), intent(in) :: dt
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Winds at start of time step.
  real(r8), intent(in) :: u(:,:), v(:,:)
  ! Wind tendencies.
  real(r8), intent(in) :: dudt(:,:), dvdt(:,:)
  ! Heating tendency.
  real(r8), intent(in) :: dsdt(:,:)
  ! Change in energy.
  real(r8), intent(out) :: de(:)

  ! Level index.
  integer :: k

  ! Net gain/loss of total energy in the column.
  de = 0.0_r8
  do k = 1, pver
     de = de + p%del(:,k)/gravit * (dsdt(:,k) + &
          dudt(:,k)*(u(:,k)+dudt(:,k)*0.5_r8*dt) + &
          dvdt(:,k)*(v(:,k)+dvdt(:,k)*0.5_r8*dt) )
  end do

end subroutine energy_change

!==========================================================================

! Subtract change in energy from the heating tendency in the levels below
! the gravity wave region.
subroutine energy_fixer(tend_level, p, de, ttgw)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Change in energy.
  real(r8), intent(in) :: de(:)
  ! Heating tendency.
  real(r8), intent(inout) :: ttgw(:,:)

  ! Column/level indices.
  integer :: i, k
  ! Energy change to apply divided by all the mass it is spread across.
  real(r8) :: de_dm(size(tend_level))

  do i = 1, size(tend_level)
     de_dm(i) = -de(i)*gravit/(p%ifc(i,pver+1)-p%ifc(i,tend_level(i)+1))
  end do

  ! Subtract net gain/loss of total energy below tend_level.
  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        ttgw(:,k) = ttgw(:,k) + de_dm
     end where
  end do

end subroutine energy_fixer

!==========================================================================

! Calculates absolute value of the local Coriolis frequency divided by the
! spatial frequency kwv, which gives a characteristic speed in m/s.
function coriolis_speed(band, lat)
  ! Inertial gravity wave lengths.
  type(GWBand), intent(in) :: band
  ! Latitude in radians.
  real(r8), intent(in) :: lat(:)

  real(r8) :: coriolis_speed(size(lat))

  ! 24*3600 = 86400 seconds in a day.
  real(r8), parameter :: omega_earth = 2._r8*pi/86400._r8

  coriolis_speed = abs(sin(lat) * 2._r8 * omega_earth / band%kwv)

end function coriolis_speed

!==========================================================================

subroutine adjust_inertial(band, tend_level, &
     u_coriolis, c, ubi, tau, ro_adjust)
  ! Inertial gravity wave lengths.
  type(GWBand), intent(in) :: band
  ! Levels above which tau is calculated.
  integer, intent(in) :: tend_level(:)
  ! Absolute value of the Coriolis frequency for each column,
  ! divided by kwv [m/s].
  real(r8), intent(in) :: u_coriolis(:)
  ! Wave propagation speed.
  real(r8), intent(in) :: c(:,-band%ngwv:)
  ! Wind speed in the direction of wave propagation.
  real(r8), intent(in) :: ubi(:,:)

  ! Tau will be adjusted by blocking wave propagation through cells where
  ! the Coriolis effect prevents it.
  real(r8), intent(inout) :: tau(:,-band%ngwv:,:)
  ! Dimensionless Coriolis term used to reduce gravity wave strength.
  ! Equal to max(0, 1 - (1/ro)^2), where ro is the Rossby number of the
  ! wind with respect to inertial waves.
  real(r8), intent(out) :: ro_adjust(:,-band%ngwv:,:)

  ! Column/level/wavenumber indices.
  integer :: i, k, l

  ! For each column and wavenumber, are we clear of levels that block
  ! upward propagation?
  logical :: unblocked_mask(size(tend_level),-band%ngwv:band%ngwv)

  unblocked_mask = .true.
  ro_adjust = 0._r8

  ! Iterate from the bottom up, through every interface level where tau is
  ! set.
  do k = maxval(tend_level)+1, ktop, -1
     do l = -band%ngwv, band%ngwv
        do i = 1, size(tend_level)
           ! Only operate on valid levels for this column.
           if (k <= tend_level(i) + 1) then
              ! Block waves if Coriolis is too strong.
              ! By setting the mask in this way, we avoid division by zero.
              unblocked_mask(i,l) = unblocked_mask(i,l) .and. &
                   (abs(ubi(i,k) - c(i,l)) > u_coriolis(i))
              if (unblocked_mask(i,l)) then
                 ro_adjust(i,l,k) = &
                      1._r8 - (u_coriolis(i)/(ubi(i,k)-c(i,l)))**2
              else
                 tau(i,l,k) = 0._r8
              end if
           end if
        end do
     end do
  end do

end subroutine adjust_inertial

end module gw_common
