module TJ2016
  !------------------------------------------------------------------------------------
  !
  ! Purpose: Implement idealized moist Held-Suarez forcings described in the TJ16 paper
  !          Thatcher, D. R. and C. Jablonowski (2016),
  !          "A moist aquaplanet variant of the Held-Suarez test
  !          for atmospheric model dynamical cores",
  !          Geosci. Model Dev., Vol. 9, 1263-1292,
  !          doi:10.5194/gmd-9-1263-2016
  !
  !          The moist simplified physics processes are based on the paper by
  !          Reed, K. A. and C. Jablonowski (2012), "Idealized tropical 
  !          cyclone simulations of intermediate complexity: A test case 
  !          for AGCMs", J. Adv. Model. Earth Syst., Vol. 4, M04001,
  !          doi:10.1029/2011MS000099 
  !
  !          The default configuration of this routine selects the 
  !          moist Held-Suarez forcing (TJ16_moist_HS). The routine can also be changed
  !          to select the Reed-Jablonowski (RJ) "simple-physics" forcing for e.g. an 
  !          idealized tropical cyclone simulation.
  !          The switch is implemented via the variable:
  !          simple_physics_option = "TJ16"    (default, moist Held-Suarez)
  !          or
  !          simple_physics_option = "RJ12"    (optional, alternative setting)
  !-----------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: pi => shr_const_pi

  implicit none
  private
  save

  public :: Thatcher_Jablonowski_set_const  ! Store constants
  public :: Thatcher_Jablonowski_precip     ! Moist physics
  public :: Thatcher_Jablonowski_sfc_pbl_hs ! Surface, PBL and Held-Suarez

  ! Private data
  real(r8)              :: gravit    ! g: gravitational acceleration (m/s2)
  real(r8)              :: cappa     ! Rd/cp
  real(r8)              :: rair      ! Rd: dry air gas constant (J/K/kg)
  real(r8)              :: cpair     ! cp: specific heat of dry air (J/K/kg)
  real(r8)              :: latvap    ! L: latent heat of vaporization (J/kg)
  real(r8)              :: rh2o      ! Rv: water vapor gas constant (J/K/kg)
  real(r8)              :: epsilo    ! Rd/Rv: ratio of h2o to dry air molecular weights
  real(r8)              :: rhoh2o    ! density of liquid water (kg/m3)
  real(r8)              :: zvir      ! (rh2o/rair) - 1, needed for virtual temperaturr
  real(r8)              :: ps0       ! Base state surface pressure (Pa)
  real(r8), allocatable :: etamid(:) ! hybrid coordinate - midpoints

CONTAINS

  subroutine Thatcher_Jablonowski_set_const(gravit_in, cappa_in, rair_in,    &
       cpair_in, latvap_in, rh2o_in, epsilo_in, rhoh2o_in, zvir_in, ps0_in, etamid_in)
    real(r8), intent(in) :: gravit_in
    real(r8), intent(in) :: cappa_in
    real(r8), intent(in) :: rair_in
    real(r8), intent(in) :: cpair_in
    real(r8), intent(in) :: latvap_in
    real(r8), intent(in) :: rh2o_in
    real(r8), intent(in) :: epsilo_in
    real(r8), intent(in) :: rhoh2o_in
    real(r8), intent(in) :: zvir_in
    real(r8), intent(in) :: ps0_in
    real(r8), intent(in) :: etamid_in(:)

    gravit = gravit_in
    cappa  = cappa_in
    rair   = rair_in
    cpair  = cpair_in
    latvap = latvap_in
    rh2o   = rh2o_in
    epsilo = epsilo_in
    rhoh2o = rhoh2o_in
    zvir   = zvir_in
    ps0    = ps0_in

    allocate(etamid(size(etamid_in)))
    etamid = etamid_in

  end subroutine Thatcher_Jablonowski_set_const


!=======================================================================
! Moist processes   
!=======================================================================
  subroutine Thatcher_Jablonowski_precip(ncol, pver, dtime,                &
       pmid, pdel, T, qv, relhum, precl, precc)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------

    integer,  intent(in)    :: ncol                    ! number of columns
    integer,  intent(in)    :: pver                    ! number of vertical levels
    real(r8), intent(in)    :: dtime                   ! time step (s)
    real(r8), intent(in)    :: pmid(ncol,pver)         ! mid-point pressure (Pa)
    real(r8), intent(in)    :: pdel(ncol,pver)         ! layer thickness (Pa)

    real(r8), intent(inout) :: T(ncol,pver)            ! temperature (K)
    real(r8), intent(inout) :: qv(ncol,pver)           ! specific humidity Q (kg/kg)

    real(r8), intent(out)   :: relhum(ncol,pver)       ! relative humidity
    real(r8), intent(out)   :: precl(ncol)             ! large-scale precipitation rate (m/s)
    real(r8), intent(out)   :: precc(ncol)             ! convective precipitation (m/s) (optional)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------

    ! Simple physics specific constants and variables

    real(r8), parameter :: T0=273.16_r8    ! control temperature (K) for calculation of qsat
    real(r8), parameter :: e0=610.78_r8    ! saturation vapor pressure (Pa) at T0 for calculation of qsat

    ! Variables for condensation and precipitation
    real(r8) :: qsat                       ! saturation value for Q (kg/kg)
    real(r8) :: tmp, tmp_t, tmp_q
    ! Loop variables
    integer  :: i, k

    !==========================================================================
    ! Set intial total, convective, and large scale precipitation rates to zero
    !==========================================================================
    precc   = 0.0_r8
    precl   = 0.0_r8

    !=========================================================================
    ! Placeholder location for an optional deep convection parameterization (not included here)
    !=========================================================================
    ! An example could be the simplified Betts-Miller (SBM) convection
    ! parameterization described in Frierson (JAS, 2007).
    ! The parameterization is expected to update 
    ! the convective precipitation rate precc and the temporary state variables
    ! T and qv. T and qv will then be updated again with the 
    ! large-scale condensation process below.

    !=========================================================================
    ! Large-Scale Condensation and Precipitation without cloud stage 
    !=========================================================================
    do k = 1, pver
      do i = 1, ncol
        qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/T(i,k))-1._r8/T0)) ! saturation value for Q
        if (qv(i,k) > qsat) then
          ! if > 100% relative humidity rain falls out
          tmp         = 1._r8/dtime*(qv(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*T(i,k)**2))) ! condensation rate
          tmp_t       = latvap/cpair*tmp       ! dT/dt tendency from large-scale condensation
          tmp_q       = -tmp                   ! dqv/dt tendency from large-scale condensation
          precl(i)    = precl(i) + tmp*pdel(i,k)/(gravit*rhoh2o) ! large-scale precipitation rate (m/s)
          T(i,k)      = T(i,k)   + tmp_t*dtime ! update T (temperature)
          qv(i,k)     = qv(i,k)  + tmp_q*dtime ! update qv (specific humidity)
          ! recompute qsat with updated T
          qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/T(i,k))-1._r8/T0)) ! saturation value for Q
        end if

        relhum(i,k) = qv(i,k) / qsat * 100._r8 ! in percent

      end do
    end do

  end subroutine Thatcher_Jablonowski_precip


!=======================================================================
! Surface fluxes and planetary boundary layer parameterization  
!=======================================================================
  subroutine Thatcher_Jablonowski_sfc_pbl_hs(ncol, pver, dtime, clat,   &
       PS, pmid, pint, lnpint, rpdel, T, U, V, qv, shflx, lhflx, taux, tauy, &
       evap, dqdt_vdiff, dtdt_vdiff, dtdt_heating, Km, Ke, Tsurf)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------

    integer,  intent(in)    :: ncol                      ! number of columns
    integer,  intent(in)    :: pver                      ! number of vertical levels
    real(r8), intent(in)    :: dtime                     ! time step (s)
    real(r8), intent(in)    :: clat(ncol)                ! latitude
    real(r8), intent(in)    :: PS(ncol)                  ! surface pressure (Pa)
    real(r8), intent(in)    :: pmid(ncol,pver)           ! mid-point pressure (Pa)
    real(r8), intent(in)    :: pint(ncol,pver+1)         ! interface pressure (Pa)
    real(r8), intent(in)    :: lnpint(ncol,2)            ! ln(interface pressure (Pa)) at and above the surface 
    real(r8), intent(in)    :: rpdel(ncol,pver)          ! reciprocal of layer thickness (Pa)

    real(r8), intent(inout) :: T(ncol,pver)              ! temperature (K)
    real(r8), intent(inout) :: U(ncol,pver)              ! zonal wind (m/s)
    real(r8), intent(inout) :: V(ncol,pver)              ! meridional wind (m/s)
    real(r8), intent(inout) :: qv(ncol,pver)             ! moisture variable (vapor form) Q (kg/kg)

    real(r8), intent(out)   :: shflx(ncol)               ! surface sensible heat flux (W/m2)
    real(r8), intent(out)   :: lhflx(ncol)               ! surface latent heat flux   (W/m2)
    real(r8), intent(out)   :: taux(ncol)                ! surface momentum flux in the zonal direction (N/m2)
    real(r8), intent(out)   :: tauy(ncol)                ! surface momentum flux in the meridional direction (N/m2)
    real(r8), intent(out)   :: evap(ncol)                ! surface water flux (kg/m2/s)
    real(r8), intent(out)   :: dqdt_vdiff(ncol,pver)     ! Q tendency due to vertical diffusion (PBL) (kg/kg/s)
    real(r8), intent(out)   :: dtdt_vdiff(ncol,pver)     ! T tendency due to vertical diffusion (PBL) in K/s
    real(r8), intent(out)   :: dtdt_heating(ncol,pver)   ! temperature tendency in K/s from relaxation
    real(r8), intent(out)   :: Km(ncol,pver+1)           ! Eddy diffusivity for boundary layer calculations
    real(r8), intent(out)   :: Ke(ncol,pver+1)           ! Eddy diffusivity for boundary layer calculations
    real(r8), intent(out)   :: Tsurf(ncol)               ! sea surface temperature K (varied by latitude)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------

    ! Constants and variables for the modified Held-Suarez forcing
    real(r8), parameter :: sec_per_day = 86400._r8       ! number of seconds per day
    real(r8), parameter :: kf=1._r8/( 1._r8*sec_per_day) ! 1./efolding_time for wind dissipation  (1/s)
    real(r8), parameter :: ka=1._r8/(40._r8*sec_per_day) ! 1./efolding_time for temperature diss. (1/s)
    real(r8), parameter :: ks=1._r8/( 4._r8*sec_per_day) ! 1./efolding_time for temperature diss. (1/s)
    real(r8), parameter :: sigmab=0.7_r8                 ! threshold sigma level (PBL level)
    real(r8), parameter :: onemsig=1._r8-sigmab          ! 1. - sigma_reference
    real(r8), parameter :: t00 = 200._r8                 ! minimum reference temperature (K)
    real(r8), parameter :: t_max=294._r8                 ! modified maximum HS equilibrium temperature (HS original is 315 K)
    real(r8), parameter :: delta_T=65._r8                ! difference in eq-polar HS equilibrium temperature (HS original is 60 K)
    real(r8), parameter :: delta_theta=10._r8            ! parameter for vertical temperature gradient (K)
    real(r8)            :: kv                            ! 1./efolding_time (normalized) for wind (1/s)
    real(r8)            :: kt                            ! 1./efolding_time for temperature diss. (1/s)
    real(r8)            :: trefa                         ! "radiative equilibrium" T (K)
    real(r8)            :: trefc                         ! used in calc of "radiative equilibrium" T

    ! Trig functions
    real(r8) :: cossq(ncol)                              ! coslat**2
    real(r8) :: cossqsq(ncol)                            ! coslat**4
    real(r8) :: sinsq(ncol)                              ! sinlat**2
    real(r8) :: coslat(ncol)                             ! cosine(latitude)

    ! Simplified physics: constants
    real(r8), parameter :: T_min      = 271._r8             ! Minimum sea surface temperature (K)
    real(r8), parameter :: del_T      = 29._r8              ! difference in eq-polar sea surface temperature (K)
    real(r8), parameter :: T_width    = 26.0_r8*pi/180.0_r8 ! width parameter for sea surface temperature (C)
    real(r8), parameter :: Tsurf_RJ12 = 302.15_r8           ! constant sea surface temperature (K) for RJ12

    real(r8), parameter :: T0=273.16_r8       ! Control temperature (K) for calculation of qsat
    real(r8), parameter :: e0=610.78_r8       ! Saturation vapor pressure (Pa) at T0 for calculation of qsat
    real(r8), parameter :: Cd0=0.0007_r8      ! Constant for calculating Cd from Smith and Vogl (2008)
    real(r8), parameter :: Cd1=0.000065_r8    ! Constant for calculating Cd from Smith and Vogl (2008)
    real(r8), parameter :: Cm=0.002_r8        ! Constant for calculating Cd from Smith and Vogl (2008)
    real(r8), parameter :: v20=20.0_r8        ! Threshold wind speed (m/s) for calculating Cd from Smith and Vogl (2008)
    real(r8)            :: C                  ! Surface exchange coefficient for sensible and latent heat, depends on simple_physics_option
    real(r8), parameter :: pbltop=85000._r8   ! Pressure (Pa) at the top of boundary layer
    real(r8), parameter :: pblconst=10000._r8 ! Constant (Pa) for the calculation of the decay of diffusivity

    ! Variables for the simple-physics and moist HS boundary layer turbulence calculation
    real(r8)            :: wind(ncol)         ! wind speed at the lowest model level (m/s)
    real(r8)            :: rho(ncol)          ! Air density near the ground (kg/m3)
    real(r8)            :: Cd(ncol)           ! Drag coefficient for momentum
    real(r8)            :: za(ncol)           ! Height at midpoint of the lowest model level (m)
    real(r8)            :: dlnpint            ! Used for calculation of heights

    ! Variables for the simple-physics and moist HS boundary layer turbulence calculation (for T and qv)
    real(r8)            :: CA(ncol,pver)      ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CC(ncol,pver)      ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CE(ncol,pver+1)    ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CFt(ncol,pver+1)   ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CFq(ncol,pver+1)   ! Matrix Coefficents for PBL Scheme

    ! Variables for the simple-physics boundary layer turbulence calculation for u and v, not used by JT16, only by RJ12 
    real(r8)            :: CAm(ncol,pver)     ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CCm(ncol,pver)     ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CEm(ncol,pver+1)   ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CFu(ncol,pver+1)   ! Matrix Coefficents for PBL Scheme
    real(r8)            :: CFv(ncol,pver+1)   ! Matrix Coefficents for PBL Scheme

    ! Variable for surface flux calculation
    real(r8)            :: qsat               ! saturation value for Q (kg/kg)

    ! Temporary storage variable
    real(r8)            :: tmp

    ! Loop variables
    integer             :: i, k

    ! Define simple_physics_option to either "TJ16" (moist HS) or "RJ12" (simple-physics)
    character(LEN=4)    :: simple_physics_option

    ! Set the simple_physics_option "TJ16" (default, moist HS)
    simple_physics_option = "TJ16"
    ! simple_physics_option = "RJ12"   ! alternative simple-physics forcing, Reed and Jablonowski (2012)

    !==========================================================================
    ! Calculate Sea Surface Temperature and set exchange coefficient
    !==========================================================================
    if (simple_physics_option == "TJ16") then
      C=0.0044_r8        ! Surface exchange coefficient for sensible and latent heat for moist HS
      do i = 1, ncol     ! set SST profile
        Tsurf(i) = del_T*exp(-(((clat(i))**2.0_r8)/(2.0_r8*(T_width**2.0_r8)))) + T_min
      end do
    else                 ! settings for RJ12
      C     = 0.0011_r8  ! Surface exchange coefficient for sensible and latent heat for simple-physics
      Tsurf = Tsurf_RJ12 ! constant SST
    endif

    !==========================================================================
    ! Pre-calculate trig functions
    !==========================================================================
    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    !==========================================================================
    ! Initialize accumulated tendencies due to Eddy diffusion
    !==========================================================================
    dqdt_vdiff = 0.0_r8
    dtdt_vdiff = 0.0_r8

    !==========================================================================
    ! Calculate hydrostatic height za of the lowermost model level
    !==========================================================================
    do i = 1, ncol
      dlnpint = (lnpint(i,2) - lnpint(i,1))
      za(i) = rair/gravit*T(i,pver)*(1._r8+zvir*qv(i,pver))*0.5_r8*dlnpint
    end do

    !==========================================================================
    ! Simple-physics surface fluxes and turbulence scheme for heat and moisture
    !
    ! The PBL parameterization is based on a simplified Ekman
    ! theory (constant Ke below 850 hPa). Ke is updated at each time step
    ! and is linked to surface conditions. First, T and Q are updated with the
    ! surface flux at the lowermost model level and then the semi-implicit
    ! PBL scheme is applied.
    !
    ! Details of the surface flux and PBL implementation can be found in:
    ! Thatcher and Jablonowski (GMD, 2016) and Reed and Jablonowski (JAMES, 2012).
    !
    ! Note that the exchange coefficient C is set to a different constant 
    ! in TJ16 and RJ12.
    !==========================================================================

    !--------------------------------------------------------------------------
    ! Compute magnitude of the low-level wind, and diffusion coeffients (Ke and Km)
    ! for PBL turbulence scheme (Eddy diffusivity), 
    ! Ke is used for heat and moisture (used by TJ16 and RJ12)
    ! Km is used for momentum (not used by TJ16, only RJ12)
    !--------------------------------------------------------------------------
    do i = 1, ncol
      wind(i) = sqrt(U(i,pver)**2 + V(i,pver)**2)   ! wind speed closest to the surface
    end do
    do i = 1, ncol
      Ke(i,pver+1) = C*wind(i)*za(i)
      if (wind(i) < v20) then                       ! if wind speed is less than 20 m/s
        Cd(i)        = Cd0+Cd1*wind(i)
        Km(i,pver+1) = Cd(i)*wind(i)*za(i)
      else
        Cd(i)        = Cm
        Km(i,pver+1) = Cm*wind(i)*za(i)
      end if
    end do

    do k = 1, pver
      do i = 1, ncol
        if( pint(i,k) >= pbltop) then
          ! keep diffusion coefficients constant below pbltop 
          Km(i,k) = Km(i,pver+1)
          Ke(i,k) = Ke(i,pver+1)
        else
          ! PBL diffusion coefficients are dragged to zero above pbltop
          Km(i,k) = Km(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)
          Ke(i,k) = Ke(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)
        end if
      end do
    end do

    !--------------------------------------------------------------------------
    ! Compute sensible and latent heat surface fluxes using an implicit approach
    ! and update the variables T and qv
    ! note: this only occurs in the lowermost model level
    !--------------------------------------------------------------------------
    do i = 1, ncol
      qsat   = epsilo*e0/PS(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))     ! saturation value for Q at the surface
      rho(i) = pmid(i,pver)/(rair * T(i,pver) *(1._r8+zvir*qv(i,pver)))          ! air density at the lowest level rho = p/(Rd Tv)

      tmp                = (T(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i)) ! new T
      dtdt_vdiff(i,pver) = (tmp-T(i,pver))/dtime                                 ! T tendency due to surface flux 
      shflx(i)           = rho(i) * cpair * C*wind(i)*(Tsurf(i)-T(i,pver))       ! sensible heat flux (W/m2)
      T(i,pver)          = tmp                                                   ! update T
 
      tmp                = (qv(i,pver)+C*wind(i)*qsat*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i)) ! new Q 
      dqdt_vdiff(i,pver) = (tmp-qv(i,pver))/dtime                                ! Q tendency due to surface flux
      lhflx(i)           = rho(i) * latvap * C*wind(i)*(qsat-qv(i,pver))         ! latent heat flux (W/m2) 
      evap(i)            = rho(i) * C*wind(i)*(qsat-qv(i,pver))                  ! surface water flux (kg/m2/s)
      qv(i,pver)         = tmp                                                   ! update Q
    end do

    if (simple_physics_option == "RJ12") then
      !--------------------------------------------------------------------------
      ! If the configuration is set to the simple-physics package by RJ12 compute
      ! surface momentum fluxes using an implicit approach and update the variables u and v 
      ! note: this only occurs in the lowermost model level and the density field rho from
      ! above is used 
      !--------------------------------------------------------------------------
      do i = 1, ncol
        tmp          = Cd(i) * wind(i)
        taux(i)      = -rho(i) * tmp * U(i,pver)                      ! zonal surface momentum flux (N/m2) 
        U(i,pver)    = U(i,pver)/(1._r8+tmp*dtime/za(i))              ! new U
        tauy(i)      = -rho(i) * tmp * V(i,pver)                      ! meridional surface momentum flux (N/m2) 
        V(i,pver)    = V(i,pver)/(1._r8+tmp*dtime/za(i))              ! new V
      enddo
    endif

    !--------------------------------------------------------------------------
    ! Calculate Diagonal Variables for PBL Scheme (semi-implicit technique follows the CESM PBL implementation)
    !--------------------------------------------------------------------------
    do k = 1, pver-1
      do i = 1, ncol
        rho(i)     = (pint(i,k+1)/(rair*(T(i,k+1)*(1._r8+zvir*qv(i,k+1))+T(i,k)*(1._r8+zvir*qv(i,k)))/2.0_r8))
        CA(i,k)    = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho(i)*rho(i)/(pmid(i,k+1)-pmid(i,k))
        CC(i,k+1)  = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho(i)*rho(i)/(pmid(i,k+1)-pmid(i,k))
        ! the next two PBL variables are initialized here for the potential use of RJ12 instead of TJ16
        ! since they need to use the same density field rho
        CAm(i,k)   = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho(i)*rho(i)/(pmid(i,k+1)-pmid(i,k))
        CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho(i)*rho(i)/(pmid(i,k+1)-pmid(i,k))
      end do
    end do
    do i = 1, ncol
      CA(i,pver)    = 0._r8
      CC(i,1)       = 0._r8
      CE(i,pver+1)  = 0._r8
      CFt(i,pver+1) = 0._r8
      CFq(i,pver+1) = 0._r8
    end do
    do i = 1, ncol
      do k = pver, 1, -1
        CE(i,k)  = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
        CFt(i,k) = ((ps0/pmid(i,k))**cappa*T(i,k)+CA(i,k)*CFt(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
        CFq(i,k) = (qv(i,k)+CA(i,k)*CFq(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
      end do
    end do

    !--------------------------------------------------------------------------
    ! Calculate the updated temperature T and moisture Q fields
    !--------------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! First: calculate the PBL mixing tendencies at the top model level
    !---------------------------------------------------------------------
    do i = 1, ncol
      tmp             = CFt(i,1)*(pmid(i,1)/ps0)**cappa         ! new T at the model top
      dtdt_vdiff(i,1) = (tmp-T(i,1))/dtime                      ! T tendency due to PBL diffusion (model top)
      T(i,1)          = tmp                                     ! update T at the model top

      dqdt_vdiff(i,1) = (CFq(i,1)-qv(i,1))/dtime                ! Q tendency due to PBL diffusion (model top)
      qv(i,1)         = CFq(i,1)                                ! update Q at the model top
    end do

    !-----------------------------------------
    ! PBL mixing at all other model levels
    !-----------------------------------------
    do i = 1, ncol
      do k = 2, pver
        tmp             = (CE(i,k)*T(i,k-1)*(ps0/pmid(i,k-1))**cappa+CFt(i,k))*(pmid(i,k)/ps0)**cappa  ! new T
        dtdt_vdiff(i,k) = dtdt_vdiff(i,k) + (tmp-T(i,k))/dtime  ! update the T tendency due to surface fluxes and the PBL diffusion
        T(i,k)          = tmp                                   ! update T

        tmp             = CE(i,k)*qv(i,k-1)+CFq(i,k)            ! new Q
        dqdt_vdiff(i,k) = dqdt_vdiff(i,k) + (tmp-qv(i,k))/dtime ! update the Q tendency due to surface fluxes and the PBL diffusion
        qv(i,k)         = tmp                                   ! update Q
      end do
    end do

    if (simple_physics_option == "TJ16") then
      !==========================================================================
      ! modified HS forcing (see Thatcher and Jablonowski (GMD, 2016))
      !--------------------------------------------------------------------------
      ! The original Held-Suarez (HS) physics algorithm is described in
      !
      !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
      !   intercomparison of the dynamical cores of atmospheric general
      !   circulation models.
      !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830
      !
      ! The modified version uses the redefined parameters: trefc, delta_T  
      !==========================================================================

      !--------------------------------------------------------------------------
      !  Compute frictional tendency from HS Rayleigh Friction (RF) at the lowest
      !  level as a diagnostic (surface momentum fluxes)
      !--------------------------------------------------------------------------
       kv  = kf*(etamid(pver) - sigmab)/onemsig                                 ! RF coefficient at the lowest level
       do i = 1, ncol
         dlnpint = (lnpint(i,2) - lnpint(i,1))
         za(i)   = rair/gravit*T(i,pver)*(1._r8+zvir*qv(i,pver))*0.5_r8*dlnpint ! height of lowest full model level
         rho(i)  = pmid(i,pver)/(rair * T(i,pver) *(1._r8+zvir*qv(i,pver)))     ! air density at the lowest level rho = p/(Rd Tv)
         taux(i) = -kv * rho(i) * U(i,pver) * za(i)                             ! U surface momentum flux in N/m2
         tauy(i) = -kv * rho(i) * V(i,pver) * za(i)                             ! V surface momentum flux in N/m2
       end do

       !--------------------------------------------------------------------------
       ! Apply HS Rayleigh Friction (RF) near the surface (below eta=0.7):
       ! represents surface stresses and PBL diffusion for U and V
       !--------------------------------------------------------------------------
       do k = 1, pver
         if (etamid(k) > sigmab) then
           kv  = kf*(etamid(k) - sigmab)/onemsig                         ! RF coefficient 
           do i=1,ncol
             U(i,k) = U(i,k) -kv*U(i,k)*dtime                            ! apply RF to U
             V(i,k) = V(i,k) -kv*V(i,k)*dtime                            ! apply RF to V
           end do
         end if
       end do

       !-----------------------------------------------------------------------
       ! Compute idealized radiative heating rates (with modified HS equilibrium temperature)
       ! mimics radiation
       !-----------------------------------------------------------------------
       do k = 1, pver
         if (etamid(k) > sigmab) then                                    ! lower atmosphere
           do i = 1, ncol
             kt = ka + (ks - ka)*cossqsq(i)*(etamid(k) - sigmab)/onemsig ! relaxation coefficent varies in the vertical
             trefc             = T_max - delta_T*sinsq(i)
             trefa             = (trefc - delta_theta*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
             trefa             = max(t00,trefa)                          ! relaxation temperature
             dtdt_heating(i,k) = (trefa - T(i,k))*kt                     ! temperature forcing due to relaxation
             T(i,k)            = T(i,k) + dtdt_heating(i,k)*dtime        ! update T
           end do
         else
           do i=1,ncol
             trefc             = T_max - delta_T*sinsq(i)
             trefa             = (trefc - delta_theta*cossq(i)*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**cappa
             trefa             = max(t00,trefa)                          ! relaxation temperature
             dtdt_heating(i,k) = (trefa - T(i,k))*ka                     ! temperature forcing due to relaxation
             T(i,k)            = T(i,k) + dtdt_heating(i,k)*dtime        ! update T
           end do
         end if
       end do

    else 
      !==========================================================================
      ! RJ12: Surface flux and PBL forcing of u and v follows the Reed-Jablonowski simple-physics configuration
      !       no HS temperature relaxation is used which limits this configuration to
      !       short simulation periods (under 30 days)
      !--------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      ! Calculate Diagonal Variables for PBL Scheme (semi-implicit technique follows the CESM PBL implementation)
      ! The fields CAm and CCm are also initialized above to guarantee the use of the same density.
      !--------------------------------------------------------------------------
      do i = 1, ncol
        CAm(i,pver)   = 0._r8
        CCm(i,1)      = 0._r8
        CEm(i,pver+1) = 0._r8
        CFu(i,pver+1) = 0._r8
        CFv(i,pver+1) = 0._r8
      end do
      do i = 1, ncol
        do k = pver, 1, -1
          CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
          CFu(i,k) = (U(i,k)+CAm(i,k)*CFu(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
          CFv(i,k) = (V(i,k)+CAm(i,k)*CFv(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
        end do
      end do

      !--------------------------------------------------------------------------
      ! Calculate the updated velocity fields U and V
      !--------------------------------------------------------------------------

      !---------------------------------------------------------------------
      ! First: calculate the PBL diffusive tendencies at the top model level
      !---------------------------------------------------------------------
      do i = 1, ncol
        U(i,1)    = CFu(i,1)                                 ! new U at the model top
        V(i,1)    = CFv(i,1)                                 ! new V at the model top
      end do

      !-----------------------------------------
      ! PBL diffusion of U and V at all other model levels
      !-----------------------------------------
      do i = 1, ncol
        do k = 2, pver
          U(i,k) = CEm(i,k)*U(i,k-1) + CFu(i,k)     ! new U
          V(i,k) = CEm(i,k)*V(i,k-1) + CFv(i,k)     ! new V
        end do
      end do
    endif

  end subroutine Thatcher_Jablonowski_sfc_pbl_hs

  !=======================================================================

end module TJ2016
