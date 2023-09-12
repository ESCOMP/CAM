module frierson
!------------------------------------------------------------------------------------
!
! Purpose: Implement idealized forcings described in
!          Frierson, et al. (2006), "A Gray-Radiation
!          Aquaplanet Moist GCM. Part I. Static Stability and Eddy Scale"
!          J. Atmos. Sci., Vol. 63, 2548-2566.
!
!          DOI: https://doi.org/10.1175/JAS3753.1
!
!====================================================================================
  !
  ! The only modules that are permitted
  !--------------------------------------
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: pi => shr_const_pi

  ! Set all Global values and routine to private by default
  ! and then explicitly set their exposure
  !---------------------------------------------------------
  implicit none
  private
  save

  public:: frierson_set_const
  public:: frierson_condensate_NONE
  public:: frierson_condensate
  public:: frierson_condensate_TJ16
  public:: frierson_condensate_USER
  public:: frierson_pbl
  public:: frierson_pbl_USER
  public:: frierson_radiation
  public:: frierson_radiation_USER

  ! Global Tuning Parameters:
  !   T0 and E0  are the temperature and saturation vapor pressure used
  !   to calculate qsat values, the saturation value for Q (kg/kg)
  !--------------------------------------------------------------------
  real(r8):: T0
  real(r8):: E0
  real(r8):: Erad
  real(r8):: Wind_min
  real(r8):: Z0
  real(r8):: Ri_c
  real(r8):: Karman
  real(r8):: Fb
  real(r8):: Rs0
  real(r8):: DeltaS
  real(r8):: Tau_eqtr
  real(r8):: Tau_pole
  real(r8):: LinFrac
  real(r8):: Boltz
  real(r8):: C_ocn

  ! Private data
  !----------------------
  real(r8),private :: gravit    ! g: gravitational acceleration (m/s2)
  real(r8),private :: cappa     ! Rd/cp
  real(r8),private :: rair      ! Rd: dry air gas constant (J/K/kg)
  real(r8),private :: cpair     ! cp: specific heat of dry air (J/K/kg)
  real(r8),private :: latvap    ! L: latent heat of vaporization (J/kg)
  real(r8),private :: rh2o      ! Rv: water vapor gas constant (J/K/kg)
  real(r8),private :: epsilo    ! Rd/Rv: ratio of h2o to dry air molecular weights
  real(r8),private :: rhoh2o    ! density of liquid water (kg/m3)
  real(r8),private :: zvir      ! (rh2o/rair) - 1, needed for virtual temperature
  real(r8),private :: ps0       ! Base state surface pressure (Pa)

  real(r8),private :: latvap_div_cpair ! latvap/cpair
  real(r8),private :: latvap_div_rh2o  ! latvap/rh2o

  real(r8),private,allocatable:: etamid(:) ! hybrid coordinate - midpoints


contains
  !=======================================================================
  subroutine frierson_set_const(I_gravit,I_cappa   ,I_rair    ,I_cpair  ,I_latvap  , &
                                I_rh2o  ,I_epsilo  ,I_rhoh2o  ,I_zvir   ,I_ps0     , &
                                I_etamid,I_T0      ,I_E0      ,I_Erad   ,I_Wind_min, &
                                I_Z0    ,I_Ri_c    ,I_Karman  ,I_Fb     ,I_Rs0     , &
                                I_DeltaS,I_Tau_eqtr,I_Tau_pole,I_LinFrac,I_Boltz   , &
                                I_Cocn                                               )
    !
    ! frierson_set_const: Set parameters and constants for the Frierson
    !                     Model fomulation. Optional inputs can be provided
    !                     to over-ride the model defaults.
    !=====================================================================

  use cam_abortutils,      only: handle_allocate_error

    !
    ! Passed Variables
    !-------------------
    real(r8),intent(in):: I_gravit
    real(r8),intent(in):: I_cappa
    real(r8),intent(in):: I_rair
    real(r8),intent(in):: I_cpair
    real(r8),intent(in):: I_latvap
    real(r8),intent(in):: I_rh2o
    real(r8),intent(in):: I_epsilo
    real(r8),intent(in):: I_rhoh2o
    real(r8),intent(in):: I_zvir
    real(r8),intent(in):: I_ps0
    real(r8),intent(in):: I_etamid(:)

    real(r8),intent(in) :: I_T0
    real(r8),intent(in) :: I_E0
    real(r8),intent(in) :: I_Erad
    real(r8),intent(in) :: I_Wind_min
    real(r8),intent(in) :: I_Z0
    real(r8),intent(in) :: I_Ri_c
    real(r8),intent(in) :: I_Karman
    real(r8),intent(in) :: I_Fb
    real(r8),intent(in) :: I_Rs0
    real(r8),intent(in) :: I_DeltaS
    real(r8),intent(in) :: I_Tau_eqtr
    real(r8),intent(in) :: I_Tau_pole
    real(r8),intent(in) :: I_LinFrac
    real(r8),intent(in) :: I_Boltz
    real(r8),intent(in) :: I_Cocn

    integer :: ierr

    ! Set global constants for later use
    !------------------------------------
    gravit   = I_gravit
    cappa    = I_cappa
    rair     = I_rair
    cpair    = I_cpair
    latvap   = I_latvap
    rh2o     = I_rh2o
    epsilo   = I_epsilo
    rhoh2o   = I_rhoh2o
    zvir     = I_zvir
    ps0      = I_ps0
    T0       = I_T0
    E0       = I_E0
    Erad     = I_Erad
    Wind_min = I_Wind_min
    Z0       = I_Z0
    Ri_c     = I_Ri_c
    Karman   = I_Karman
    Fb       = I_Fb
    Rs0      = I_Rs0
    DeltaS   = I_DeltaS
    Tau_eqtr = I_Tau_eqtr
    Tau_pole = I_Tau_pole
    LinFrac  = I_LinFrac
    Boltz    = I_Boltz
    C_ocn    = I_Cocn

    latvap_div_cpair = latvap/cpair
    latvap_div_rh2o  = latvap/rh2o

    ! allocate space and set the level information
    !----------------------------------------------
    allocate(etamid(size(I_etamid)),stat=ierr)
    if (ierr /= 0) then
      call handle_allocate_error(ierr, 'frierson_set_const', 'etamid')
    end if

    etamid = I_etamid

  end subroutine frierson_set_const
  !=======================================================================


  !=======================================================================
  subroutine frierson_condensate_NONE(ncol,pver,pmid,T,qv,relhum,precl)
    !
    ! Precip_process: Implement NO large-scale condensation/precipitation
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)
    real(r8),intent(out)  :: relhum(ncol,pver) ! relative humidity
    real(r8),intent(out)  :: precl (ncol)      ! large-scale precipitation rate (m/s)
    !
    ! Local Values
    !-------------
    real(r8):: qsat
    integer :: i, k

    ! Set large scale precipitation rates to zero
    !--------------------------------------------------------------------------
    precl(:) = 0.0_r8

    ! Large-Scale Condensation and Precipitation without cloud stage
    !---------------------------------------------------------------
    do k = 1, pver
       do i = 1, ncol
         ! calculate saturation value for Q
         !----------------------------------
         qsat = epsilo*E0/pmid(i,k)*exp(-latvap_div_rh2o*((1._r8/T(i,k))-1._r8/T0))

         ! Set percent relative humidity
         !-------------------------------
         relhum(i,k) = (qv(i,k)/qsat)*100._r8
       end do
    end do

  end subroutine frierson_condensate_NONE
  !=======================================================================


  !=======================================================================
  subroutine frierson_condensate(ncol,pver,dtime,pmid,pdel,T,qv,relhum,precl,evapdt,evapdq)
    !
    ! Precip_process: Implement large-scale condensation/precipitation
    !                 from Frierson 2006.
    !
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: dtime             ! time step (s)
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pdel  (ncol,pver) ! layer thickness (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)
    real(r8),intent(out)  :: relhum(ncol,pver) ! relative humidity
    real(r8),intent(out)  :: precl (ncol)      ! large-scale precipitation rate (m/s)
    real(r8),intent(out)  :: evapdt(ncol,pver) ! T tendency due to re-evaporation
    real(r8),intent(out)  :: evapdq(ncol,pver) ! Q tendency due to re-evaporation
    !
    ! Local Values
    !-------------
    logical,parameter :: do_evap = .true.

    real(r8):: esat (ncol,pver)
    real(r8):: qsat (ncol,pver)
    real(r8):: dqsat(ncol,pver)
    real(r8):: qdel (ncol,pver)
    real(r8):: tdel (ncol,pver)
    real(r8):: qnew (ncol,pver)
    real(r8):: tnew (ncol,pver)
    real(r8):: qext (ncol)
    real(r8):: qdef (ncol)

    integer :: i, k

    ! Large-Scale Condensation and Precipitation
    !--------------------------------------------
    do k = 1,pver

      ! calculate saturation vapor pressure
      !-------------------------------------
      esat(:,k) = E0*exp(-(latvap_div_rh2o)*((1._r8/T(:,k))-1._r8/T0))

      ! calculate saturation value for Q
      !----------------------------------
      do i = 1,ncol
        if(pmid(i,k) > (1._r8-epsilo)*esat(i,k)) then
          qsat (i,k) = epsilo*esat(i,k)/pmid(i,k)
          dqsat(i,k) = (latvap_div_rh2o)*qsat(i,k)/(T(i,k)**2)
        else
          qsat (i,k) = 0._r8
          dqsat(i,k) = 0._r8
        endif
      end do

      ! if > 100% relative humidity, rain falls out
      !---------------------------------------------
      where(((qv(:,k)-qsat(:,k))*qsat(:,k)) > 0._r8)
        qdel (:,k) = (qsat(:,k)-qv(:,k))/(1._r8+(latvap_div_cpair)*dqsat(:,k))
        tdel (:,k) = -(latvap_div_cpair)*qdel(:,k)
      else where
        qdel (:,k) = 0._r8
        tdel (:,k) = 0._r8
      end where

      ! Update temperature and water vapor
      !-----------------------------------
      qnew(:,k) = qv(:,k) + qdel(:,k)
      tnew(:,k) =  T(:,k) + tdel(:,k)
    end do

    ! optionally allow for re-evaporation of falling precip
    !-------------------------------------------------------
    if(do_evap) then
      ! Initialize work array for excess Q
      !--------------------------------------
      qext(:) = 0._r8

      ! Loop down thru the model levels
      !---------------------------------
      do k = 1, pver

        ! Add qdel for the current level to the excess
        !----------------------------------------------
        where(qdel(:,k) < 0._r8) qext(:) = qext(:) - qdel(:,k)*pdel(:,k)/gravit

        ! Evaporate excess Q where needed
        !----------------------------------
        qdef(:) = 0._r8
        where((qdel(:,k) >= 0._r8).and.(qext(:) > 0._r8))
          qext(:)   = qext(:)*gravit/pdel(:,k)
          qdef(:)   = (qsat(:,k)-qv(:,k))/(1._r8+(latvap_div_cpair)*dqsat(:,k))
          qdef(:)   = min(qext(:),max(qdef(:),0._r8))
          qdel(:,k) = qdel(:,k) + qdef(:)
          tdel(:,k) = tdel(:,k) -(latvap_div_cpair)*qdef(:)
          qext(:)   = (qext(:)-qdef(:))*pdel(:,k)/gravit

          ! Update temperature and water vapor
          !-----------------------------------
          qnew(:,k) = qv(:,k) + qdel(:,k)
          tnew(:,k) =  T(:,k) + tdel(:,k)
        end where

        ! Save T/Q tendencies due to re-evaporation
        !--------------------------------------------
        evapdq(:,k) =  qdef(:)/dtime
        evapdt(:,k) = -qdef(:)*latvap_div_cpair/dtime
      end do ! k = 1, pver
    else
      ! Set T/Q re-evaporation tendencies to 0
      !--------------------------------------------
      evapdt(:,:) = 0._r8
      evapdq(:,:) = 0._r8
    endif

    ! Set large scale precipitation rates to zero
    !--------------------------------------------------------------------------
    precl(:) = 0.0_r8

    ! Calculate resulting precip value and relative humidity
    !--------------------------------------------------------
    do k = 1, pver
      precl (:)   = precl(:) - (qdel(:,k)*pdel(:,k))/(gravit*rhoh2o)
      qsat  (:,k) = (epsilo/pmid(:,k))*E0*exp(-latvap_div_rh2o*((1._r8/tnew(:,k))-1._r8/T0))
      relhum(:,k) = (qnew(:,k)/qsat (:,k))*100._r8
    end do
    precl(:) = max(precl(:),0._r8)/dtime

    ! Update T and qv values due to precipitation
    !--------------------------------------------
    qv(:,:) = qnew(:,:)
    T (:,:) = tnew(:,:)

  end subroutine frierson_condensate
  !=======================================================================


  !=======================================================================
  subroutine frierson_condensate_TJ16(ncol,pver,dtime,pmid,pdel,T,qv,relhum,precl)
    !
    ! Precip_process: Implement large-scale condensation/precipitation
    !                 from TJ16.
    !
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: dtime             ! time step (s)
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pdel  (ncol,pver) ! layer thickness (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)
    real(r8),intent(out)  :: relhum(ncol,pver) ! relative humidity
    real(r8),intent(out)  :: precl (ncol)      ! large-scale precipitation rate (m/s)
    !
    ! Local Values
    !-------------
    real(r8):: qsat
    real(r8):: Crate
    integer :: i, k

    ! Set large scale precipitation rates to zero
    !--------------------------------------------------------------------------
    precl(:) = 0.0_r8

    ! Large-Scale Condensation and Precipitation without cloud stage
    !---------------------------------------------------------------
    do k = 1, pver
       do i = 1, ncol
         ! calculate saturation value for Q
         !----------------------------------
         qsat = epsilo*E0/pmid(i,k)*exp(-latvap_div_rh2o*((1._r8/T(i,k))-1._r8/T0))

         ! if > 100% relative humidity rain falls out
         !-------------------------------------------
         if(qv(i,k) > qsat) then
           ! calc the condensation and large-scale precipitation(m/s) rates
           !-------------------------------------------------------------------
           Crate    =  ((qv(i,k)-qsat)/dtime)                                     &
                      /(1._r8+(latvap_div_cpair)*(epsilo*latvap*qsat/(rair*T(i,k)**2)))
           precl(i) = precl(i) + (Crate*pdel(i,k))/(gravit*rhoh2o)

           ! Update T and qv values due to precipitation
           !--------------------------------------------
           T (i,k) = T (i,k) + Crate*(latvap_div_cpair)*dtime
           qv(i,k) = qv(i,k) - Crate*dtime

           ! recompute qsat with updated T
           !-------------------------------
           qsat = epsilo*E0/pmid(i,k)*exp(-latvap_div_rh2o*((1._r8/T(i,k))-1._r8/T0))
         endif

         ! Set percent relative humidity
         !-------------------------------
         relhum(i,k) = (qv(i,k)/qsat)*100._r8
       end do
    end do

  end subroutine frierson_condensate_TJ16
  !=======================================================================


  !=======================================================================
  subroutine frierson_condensate_USER(ncol,pver,dtime,pmid,pdel,T,qv,relhum,precl)
    !
    ! frierson_condensate_USER: This routine is a stub which users can use
    !                           to develop and test their own large scale
    !                           condensation scheme
    !=======================================================================
    !
    ! Passed Variables
    !---------------------
    integer ,intent(in)   :: ncol              ! number of columns
    integer ,intent(in)   :: pver              ! number of vertical levels
    real(r8),intent(in)   :: dtime             ! time step (s)
    real(r8),intent(in)   :: pmid  (ncol,pver) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pdel  (ncol,pver) ! layer thickness (Pa)
    real(r8),intent(inout):: T     (ncol,pver) ! temperature (K)
    real(r8),intent(inout):: qv    (ncol,pver) ! specific humidity Q (kg/kg)
    real(r8),intent(out)  :: relhum(ncol,pver) ! relative humidity
    real(r8),intent(out)  :: precl (ncol)      ! large-scale precipitation rate (m/s)
    !
    ! Local Values
    !-------------
    real(r8):: qsat
    integer :: i, k

    ! Set large scale precipitation rates to zero
    !--------------------------------------------------------------------------
    precl(:) = 0.0_r8

    ! Large-Scale Condensation and Precipitation without cloud stage
    !---------------------------------------------------------------
    do k = 1, pver
      do i = 1, ncol
        ! calculate saturation value for Q
        !----------------------------------
        qsat = epsilo*E0/pmid(i,k)*exp(-latvap_div_rh2o*((1._r8/T(i,k))-1._r8/T0))

        ! Set percent relative humidity
        !-------------------------------
        relhum(i,k) = (qv(i,k)/qsat)*100._r8
      end do
    end do

  end subroutine frierson_condensate_USER
  !=======================================================================


  !=======================================================================
  subroutine frierson_pbl(ncol, pver, dtime, pmid, pint, Zm, Zi,             &
                          Psfc, Tsfc, Qsfc, T, U, V, Q,  Fsw, Fdn,           &
                          Cdrag, Km, Ke, VSE, Z_pbl, Rf, dQa, dTa, dUa, dVa, &
                          LHflux, SHflux, TUflux, TVflux                     )
    !
    ! The implicit PBL parameterization based on Frierson, et al. 2006.
    !
    ! frierson_pbl: This is an implementation of the implicit computation
    !               derived from the code of the Frierson model. The
    !               calculations are roughly divided up into sections of
    !               the model where they should be carried out.
    !
    !==========================================================================
    !
    ! Passed Variables
    !------------------
    integer ,intent(in)   :: ncol                 ! Number of columns
    integer ,intent(in)   :: pver                 ! Number of levels
    real(r8),intent(in)   :: dtime                ! Time Step
    real(r8),intent(in)   :: pmid  (ncol,pver)    ! Pressure at model levels
    real(r8),intent(in)   :: pint  (ncol,pver+1)  ! Pressure at interface levels.
    real(r8),intent(in)   :: Zm    (ncol,pver)    ! Height at model levels.
    real(r8),intent(in)   :: Zi    (ncol,pver)    ! Height at interface levels.
    real(r8),intent(in)   :: Psfc  (ncol)         ! Surface Pressure.
    real(r8),intent(inout):: Tsfc  (ncol)         ! SST temperature K
    real(r8),intent(inout):: Qsfc  (ncol)         ! sea surface water vapor (kg/kg)
    real(r8),intent(inout):: T     (ncol,pver)    ! ATM Temperature values.
    real(r8),intent(inout):: U     (ncol,pver)    ! ATM Zonal Wind values.
    real(r8),intent(inout):: V     (ncol,pver)    ! ATM Meridional Wind values.
    real(r8),intent(inout):: Q     (ncol,pver)    ! ATM Water vapor values.
    real(r8),intent(in)   :: Fdn   (ncol)         ! Downward LW flux at surface
    real(r8),intent(in)   :: Fsw   (ncol)         ! Net SW flux at surface from gray radiation
    real(r8),intent(out)  :: Cdrag (ncol)         ! Surface drage coef.
    real(r8),intent(out)  :: Km    (ncol,pver+1)  ! Eddy diffusivity for PBL (momentum)
    real(r8),intent(out)  :: Ke    (ncol,pver+1)  ! Eddy diffusivity for PBL
    real(r8),intent(out)  :: VSE   (ncol,pver)    ! Virtual-Dry Static energy
    real(r8),intent(out)  :: Z_pbl (ncol)         ! Height of PBL layer.
    real(r8),intent(out)  :: Rf    (ncol,pver)
    real(r8),intent(out)  :: dTa   (ncol,pver)
    real(r8),intent(out)  :: dQa   (ncol,pver)
    real(r8),intent(out)  :: dUa   (ncol,pver)
    real(r8),intent(out)  :: dVa   (ncol,pver)
    real(r8),intent(out)  :: LHflux(ncol)         ! Latent Heat Flux
    real(r8),intent(out)  :: SHflux(ncol)         ! Sensible Heat Flux
    real(r8),intent(out)  :: TUflux(ncol)         ! U Surface stress
    real(r8),intent(out)  :: TVflux(ncol)         ! V Surface stress
    !
    ! Local Values
    !---------------
    real(r8):: Tv    (ncol,pver)
    real(r8):: Thv   (ncol,pver)
    real(r8):: Ws    (ncol,pver)
    real(r8):: rho   (ncol)  ! Air density near the ground (kg/m3)
    real(r8):: Z_sfc (ncol)
    real(r8):: Rf_sfc(ncol)
    real(r8):: Ri_a  (ncol)
    real(r8):: Ri    (ncol,pver)
    integer :: K_sfc (ncol)
    integer :: K_pbl (ncol)
    real(r8):: Ke_pbl(ncol)
    real(r8):: Km_pbl(ncol)
    real(r8):: Z_a   (ncol)  ! Height at midpoint of the lowest model level (m)
    real(r8):: Ws_a  (ncol)  ! wind speed at the lowest model level (m/s)
    real(r8):: Thv_a (ncol)
    real(r8):: Thv_s (ncol)
    real(r8):: Ustar (ncol)
    real(r8):: Bstar (ncol)
    real(r8):: ZETA,PHI

    real(r8):: MU (ncol,pver)
    real(r8):: NUe(ncol,pver)
    real(r8):: NUm(ncol,pver)
    real(r8):: Am (ncol,pver)
    real(r8):: Bm (ncol,pver)
    real(r8):: Cm (ncol,pver)
    real(r8):: Ae (ncol,pver)
    real(r8):: Be (ncol,pver)
    real(r8):: Ce (ncol,pver)
    real(r8):: FLu(ncol,pver)
    real(r8):: FLv(ncol,pver)
    real(r8):: FLq(ncol,pver)
    real(r8):: FLt(ncol,pver)
    real(r8):: Et (ncol,pver)
    real(r8):: Eq (ncol,pver)
    real(r8):: Eu (ncol,pver)
    real(r8):: Ev (ncol,pver)

    real(r8):: Fval_t(ncol,pver)
    real(r8):: Fval_q(ncol,pver)
    real(r8):: Fval_u(ncol,pver)
    real(r8):: Fval_v(ncol,pver)
    real(r8):: Eval_m(ncol,pver)
    real(r8):: Eval_e(ncol,pver)
    integer :: i, k

    real(r8):: Su(ncol,pver)
    real(r8):: Sv(ncol,pver)
    real(r8):: St(ncol,pver)
    real(r8):: Sq(ncol,pver)

    real(r8)::  Estar_u(ncol)
    real(r8)::  Estar_v(ncol)
    real(r8)::  Estar_q(ncol)
    real(r8)::  Estar_t(ncol)
    real(r8)::  dFa_dTa(ncol)
    real(r8)::  dFa_dQa(ncol)
    real(r8)::  dFa_dUa(ncol)
    real(r8)::  dFa_dVa(ncol)

    real(r8)::  Th_a    (ncol)
    real(r8)::  Th_s    (ncol)
    real(r8)::  Ft      (ncol)
    real(r8)::  dFt_dTa (ncol)
    real(r8)::  dFt_dTs (ncol)
    real(r8)::  Fq      (ncol)
    real(r8)::  dFq_dQa (ncol)
    real(r8)::  dFq_dTs (ncol)
    real(r8)::  Fu      (ncol)
    real(r8)::  dFu_dUa (ncol)
    real(r8)::  Fv      (ncol)
    real(r8)::  dFv_dVa (ncol)
    real(r8)::  Fup     (ncol)
    real(r8)::  dFup_dTs(ncol)

    real(r8)::  FN_u (ncol)
    real(r8)::  FN_v (ncol)
    real(r8)::  EN_t (ncol)
    real(r8)::  FN_t (ncol)
    real(r8)::  EN_q (ncol)
    real(r8)::  FN_q (ncol)
    real(r8)::  Flux (ncol)
    real(r8)::  dFlux(ncol)
    real(r8)::  dTs  (ncol)

    real(r8):: Tsfc_bc(ncol)

  !============================================================================
  ! <PHYSICS> tphysbc():
  !
  !     Required Values:
  !       T(:,:),Q(:,:),U(:,:),V(:,:)
  !       Pmid(:,:),Pint(:,:),Zm(:,:),Zi(:,:)
  !       Tsfc(:),Qsfc(:),Psfc(:)
  !============================================================================

    ! Sx() values allow for explicit source tendencies to be passed to
    ! implicit PBL calculation.  Set all values to 0. for now.
    !-------------------------------------------------------------------------
    Su(:,:) = 0._r8
    Sv(:,:) = 0._r8
    St(:,:) = 0._r8
    Sq(:,:) = 0._r8

    ! Calc some values we will need later on
    !------------------------------------------
    do k = 1, pver
      Ws (:,k) = sqrt(U(:,k)**2 + V(:,k)**2 + Wind_min)
      Tv (:,k) = T (:,k)*(1._r8+zvir*Q(:,k))
      Thv(:,k) = Tv(:,k)*((ps0/pmid(:,k))**cappa)
      VSE(:,k) = Tv(:,k)+gravit*Zm(:,k)/cpair
    end do

    ! Calculate Drag Coef and related values
    !-----------------------------------------
    do i = 1,ncol
      Z_a  (i) = Zm (i,pver)
      Ws_a (i) = Ws (i,pver)
      Thv_a(i) = Thv(i,pver)
      Thv_s(i) = Tsfc(i)*(1._r8+zvir*Qsfc(i)  )*((ps0/Psfc(i))**cappa)
      Ri_a (i) = (gravit*Z_a(i)/(Ws_a(i)**2))*(Thv_a(i)-Thv_s(i))/Thv_s(i)
      if(Ri_a(i) <= 0._r8) then
        Cdrag(i) = (Karman/log((Z_a(i)/Z0)))**2
      elseif(Ri_a(i) >= Ri_c) then
        Cdrag(i) = 0._r8
      else
        Cdrag(i) = ((1._r8-(Ri_a(i)/Ri_c))*Karman/log((Z_a(i)/Z0)))**2
      endif
      Ustar(i) = sqrt(Cdrag(i))*Ws_a(i)
      Bstar(i) = sqrt(Cdrag(i))*(gravit*(Thv_a(i)-Thv_s(i))/Thv_s(i))
    end do

    ! Calculate a bulk Richardson number and determine
    ! depths of boundary/surface layers.
    !----------------------------------------------------
    do k = 1,pver
     Ri(:,k) = (gravit*Zm(:,k)/(Ws(:,k)**2))*(VSE(:,k)-VSE(:,pver))/VSE(:,pver)
     Rf(:,k) = Ri(:,k)/Ri_c
    end do

    do i =1,ncol
      Z_pbl(i) = Zm(i,pver)
      K_pbl(i) = pver
      do k = (pver-1),1,-1
        if(Rf(i,k) > 1._r8) then
          K_pbl(i) = k + 1
          Z_pbl(i) = (Zm(i,k+1)*(Rf(i,k)- 1._r8    )                    &
                     +Zm(i,k  )*( 1._r8 - Rf(i,k+1)))/(Rf(i,k)-Rf(i,k+1))
          exit
        endif
      end do

      ! surface layer height is a fixed fraction of the PBL
      ! determine the corresponding level index and Rf value
      !-----------------------------------------------------
      Z_sfc(i) = Fb*Z_pbl(i)
      K_sfc(i) = pver
      do k = (pver-1),1,-1
        if(Zm(i,k) > Z_sfc(i)) then
          K_sfc (i) = k + 1
          Rf_sfc(i) = (Rf(i,k+1)*(Zm(i,k)  - Z_sfc(i) )                    &
                     + Rf(i,k  )*(Z_sfc(i) - Zm(i,k+1)))/(Zm(i,k)-Zm(i,k+1))
          exit
        endif
      end do
    end do ! i =1,ncol

    ! Compute diffusion coefs
    !-------------------------
    Ke(:,:)   = 0._r8
    Ke_pbl(:) = 0._r8
    do i = 1,ncol
      if (Cdrag(i) /= 0._r8) then
        do k = pver,K_pbl(i),-1
          ZETA = Zi(i,k)*Karman*Bstar(i)/(Ustar(i)*Ustar(i))
          if(ZETA < 0._r8) then
            if(k >= K_sfc(i)) then
              Ke(i,k) = Karman*Ustar(i)*Zi(i,k)
            else
              Ke(i,k) = Karman*Ustar(i)*Zi(i,k)                  &
                              *(((Z_pbl(i)-Zi(i,k))/(Z_pbl(i)-Z_sfc(i)))**2)
            endif
          elseif (ZETA < Ri_c) then
            PHI = 1._r8 + ZETA
            if(k >= K_sfc(i)) then
              Ke(i,k) = Karman*Ustar(i)*Zi(i,k)/PHI
            else
              Ke(i,k) = Karman*Ustar(i)*Zi(i,k)                  &
                              *(((Z_pbl(i)-Zi(i,k))/(Z_pbl(i)-Z_sfc(i)))**2)/PHI
            endif
          endif
        end do
        Ke_pbl(i) = Ke(i,K_sfc(i))*Z_sfc(i)/Zi(i,K_sfc(i))
      endif
    end do

    ! The Same coefs used for momentum
    !-----------------------------------
    Km(:,:)   = Ke(:,:)
    Km_pbl(:) = Ke_pbl(:)

    ! Compute downward values for the implicit PBL scheme
    !-----------------------------------------------------
    do k = 1,pver
      MU (:,k) = gravit*dtime/(Pint(:,k+1) - Pint(:,k))
    end do

    NUe(:,:) = 0._r8
    NUm(:,:) = 0._r8
    do k = 2,pver
      rho(:)   = 2._r8*Pint(:,k)/(rair*(Tv(:,k)+Tv(:,k-1)))
      NUe(:,k) = rho(:)*Ke(:,k)/(Zm(:,k)-Zm(:,k-1))
      NUm(:,k) = rho(:)*Km(:,k)/(Zm(:,k)-Zm(:,k-1))
    end do

    Am(:,1   ) = MU(:,1)*NUm(:,2)
    Cm(:,1   ) = 0._r8
    Am(:,pver) = 0._r8
    Cm(:,pver) = MU(:,pver)*NUm(:,pver)
    Ae(:,1   ) = MU(:,1   )*NUe(:,2)
    Ce(:,1   ) = 0._r8
    Ae(:,pver) = 0._r8
    Ce(:,pver) = MU(:,pver)*NUe(:,pver)
    do k = 2,(pver-1)
      Am(:,k) = MU(:,k)*NUm(:,k+1)
      Cm(:,k) = MU(:,k)*NUm(:,k  )
      Ae(:,k) = MU(:,k)*NUe(:,k+1)
      Ce(:,k) = MU(:,k)*NUe(:,k  )
    end do
    Bm(:,:) = 1._r8 - Am(:,:) - Cm(:,:)
    Be(:,:) = 1._r8 - Ae(:,:) - Ce(:,:)

    FLu(:,1) = 0._r8
    FLv(:,1) = 0._r8
    FLq(:,1) = 0._r8
    FLt(:,1) = 0._r8
    do k = 2,pver
      FLu(:,k) = NUm(:,k)*(U  (:,k)-U  (:,k-1))
      FLv(:,k) = NUm(:,k)*(V  (:,k)-V  (:,k-1))
      FLq(:,k) = NUe(:,k)*(Q  (:,k)-Q  (:,k-1))
      FLt(:,k) = NUe(:,k)*(VSE(:,k)-VSE(:,k-1))
    end do
    do k = 1,(pver-1)
      Eu(:,k) = Su(:,k) + MU(:,k)*(FLu(:,k)-FLu(:,k+1))
      Ev(:,k) = Sv(:,k) + MU(:,k)*(FLv(:,k)-FLv(:,k+1))
      Eq(:,k) = Sq(:,k) + MU(:,k)*(FLq(:,k)-FLq(:,k+1))
      Et(:,k) = St(:,k) + MU(:,k)*(FLt(:,k)-FLt(:,k+1))
    end do
    Eu(:,pver) = Su(:,pver) + MU(:,pver)*FLu(:,pver)
    Ev(:,pver) = Sv(:,pver) + MU(:,pver)*FLv(:,pver)
    Eq(:,pver) = Sq(:,pver) + MU(:,pver)*FLq(:,pver)
    Et(:,pver) = St(:,pver) + MU(:,pver)*FLt(:,pver)

    Eval_m(:,1) = -Am(:,1)/Bm(:,1)
    Eval_e(:,1) = -Ae(:,1)/Be(:,1)
    Fval_u(:,1) =  Eu(:,1)/Bm(:,1)
    Fval_v(:,1) =  Ev(:,1)/Bm(:,1)
    Fval_q(:,1) =  Eq(:,1)/Be(:,1)
    Fval_t(:,1) =  Et(:,1)/Be(:,1)
    do k = 2,(pver-1)
      Eval_m(:,k) = -Am(:,k)/(Bm(:,k)+Cm(:,k)*Eval_m(:,k-1))
      Eval_e(:,k) = -Ae(:,k)/(Be(:,k)+Ce(:,k)*Eval_e(:,k-1))
      Fval_u(:,k) =  (Eu(:,k)-Cm(:,k)*Fval_u(:,k-1)) &
                    /(Bm(:,k)+Cm(:,k)*Eval_m(:,k-1))
      Fval_v(:,k) =  (Ev(:,k)-Cm(:,k)*Fval_v(:,k-1)) &
                    /(Bm(:,k)+Cm(:,k)*Eval_m(:,k-1))
      Fval_q(:,k) =  (Eq(:,k)-Ce(:,k)*Fval_q(:,k-1)) &
                    /(Be(:,k)+Ce(:,k)*Eval_e(:,k-1))
      Fval_t(:,k) =  (Et(:,k)-Ce(:,k)*Fval_t(:,k-1)) &
                    /(Be(:,k)+Ce(:,k)*Eval_e(:,k-1))
    end do
    Eval_m(:,pver) = 0._r8
    Eval_e(:,pver) = 0._r8
    Fval_u(:,pver) = 0._r8
    Fval_v(:,pver) = 0._r8
    Fval_q(:,pver) = 0._r8
    Fval_t(:,pver) = 0._r8

    Estar_u(:) = (Eu(:,pver)-Cm(:,pver)*Fval_u(:,pver-1))
    Estar_v(:) = (Ev(:,pver)-Cm(:,pver)*Fval_v(:,pver-1))
    Estar_q(:) = (Eq(:,pver)-Ce(:,pver)*Fval_q(:,pver-1))
    Estar_t(:) = (Et(:,pver)-Ce(:,pver)*Fval_t(:,pver-1))

    dFa_dTa(:) = NUe(:,pver)*(1._r8-Eval_e(:,pver-1))
    dFa_dQa(:) = NUe(:,pver)*(1._r8-Eval_e(:,pver-1))
    dFa_dUa(:) = NUm(:,pver)*(1._r8-Eval_m(:,pver-1))
    dFa_dVa(:) = NUm(:,pver)*(1._r8-Eval_m(:,pver-1))

  !============================================================================
  ! <SHARE> flux calculation():
  !
  !     Required Values:
  !       Passed from <PHYSICS>:
  !           T(:,pver),Q(:,pver),U(:,pver),V(:,pver),Cdrag(:),Pmid(:,pver)
  !           MU(:,pver),dFa_dTa(:),dFa_dQa(:),dFa_dUa(:),dFa_dVa(:)
  !           Estar_t(:),Estar_q(:),Estar_u(:),Estar_v(:)
  !       Passed from <DOCN>:
  !           Tsfc(:),Qsfc(:),Psfc(:)
  !
  !============================================================================

    ! Calculate Surface flux values and their derivatives
    !--------------------------------------------------------
    do i = 1, ncol
      Th_a(i) = T    (i,pver)*((ps0/pmid(i,pver))**cappa)
      Th_s(i) = Tsfc(i)     *((ps0/Psfc  (i)  )**cappa)
      rho (i) = pmid (i,pver)/(rair*Tv(i,pver))

      Ft     (i) =  rho(i)*Cdrag(i)*Ws_a(i)*(Th_s (i) - Th_a(i))
      Fq     (i) =  rho(i)*Cdrag(i)*Ws_a(i)*(Qsfc(i) - Q(i,pver))
      Fu     (i) = -rho(i)*Cdrag(i)*Ws_a(i)*U(i,pver)
      Fv     (i) = -rho(i)*Cdrag(i)*Ws_a(i)*V(i,pver)
      Fup    (i) =  Boltz*Tsfc(i)**4

      dFt_dTa(i) = -rho(i)*Cdrag(i)*Ws_a(i)*((ps0/pmid(i,pver))**cappa)
      dFq_dQa(i) = -rho(i)*Cdrag(i)*Ws_a(i)
      dFu_dUa(i) = -rho(i)*Cdrag(i)*Ws_a(i)
      dFv_dVa(i) = -rho(i)*Cdrag(i)*Ws_a(i)

      dFt_dTs(i) =  rho(i)*Cdrag(i)*Ws_a(i)*((ps0/Psfc  (i)  )**cappa)
      dFq_dTs(i) =  rho(i)*Cdrag(i)*Ws_a(i)*Qsfc(i)*latvap/(rh2o*(Tsfc(i)**2))
      dFup_dTs(i) = 4._r8*Boltz*Tsfc(i)**3
    end do

    ! Incorporate surface fluxes into implicit scheme, then
    ! update flux values and derivatives
    !------------------------------------------
    FN_u   (:) =  (Estar_u(:) + MU(:,pver)*Fu(:))          &
                 /(1._r8-MU(:,pver)*(dFa_dUa(:)+dFu_dUa(:)))
    FN_v   (:) =  (Estar_v(:) + MU(:,pver)*Fv(:))          &
                 /(1._r8-MU(:,pver)*(dFa_dVa(:)+dFv_dVa(:)))
    FN_t   (:) =  (Estar_t(:) + MU(:,pver)*Ft(:))          &
                 /(1._r8-MU(:,pver)*(dFa_dTa(:)+dFt_dTa(:)))
    FN_q   (:) =  (Estar_q(:) + MU(:,pver)*Fq(:))          &
                 /(1._r8-MU(:,pver)*(dFa_dQa(:)+dFq_dQa(:)))

    EN_t   (:) =  (             MU(:,pver)*dFt_dTs(:) )    &
                 /(1._r8-MU(:,pver)*(dFa_dTa(:)+dFt_dTa(:)))
    EN_q   (:) =  (             MU(:,pver)*dFq_dTs(:) )    &
                 /(1._r8-MU(:,pver)*(dFa_dQa(:)+dFq_dQa(:)))

    Ft     (:) =      Ft(:) + dFt_dTa(:)*FN_t(:)
    Fq     (:) =      Fq(:) + dFq_dQa(:)*FN_q(:)

    dFt_dTs(:) = dFt_dTs(:) + dFt_dTa(:)*EN_t(:)
    dFq_dTs(:) = dFq_dTs(:) + dFq_dQa(:)*EN_q(:)

  !============================================================================
  ! <DOCN> surface calculation():
  !
  !     Required Values:
  !       Passed from <SHARE>:
  !           Fup(:),Ft(:),Fq(:)
  !           dFup_dTs(:),dFt_dTs(:),dFq_dTs(:)
  !           Fsw(:)
  !       Passed from <PHYSICS>:
  !           Fdn(:)
  !
  !============================================================================

    ! Update surface values
    !-----------------------
    Tsfc_bc(:) = Tsfc(:)

    Flux (:) = (dtime/C_ocn)*(Fdn(:) -     Fup(:) +      Fsw(:)  &
                                     -cpair*Ft(:) -latvap*Fq(:)  )
    dFlux(:) = (dtime/C_ocn)*(-dFup_dTs(:) -cpair*dFt_dTs(:) -latvap*dFq_dTs(:))
    Tsfc (:) = Tsfc(:) + (Flux(:)/(1._r8-dFlux(:)))
    Qsfc (:) = epsilo*E0/Psfc(:)*exp(-latvap_div_rh2o*((1._r8/Tsfc(:))-1._r8/T0))
    dTs  (:) = Tsfc(:) - Tsfc_bc(:)

    LHflux(:) = latvap*Fq(:)
    SHflux(:) = cpair *Ft(:)
    TUflux(:) = Fu(:)
    TVflux(:) = Fv(:)

  !============================================================================
  ! <PHYSICS> tphysac():
  !
  !     Required Values:
  !       Passed from <SHARE>:
  !           FN_t(:),FN_q(:),FN_u(:),FN_v(:)
  !           EN_t(:),EN_q(:),dTs(:)
  !       Passed from <PHYSICS>:
  !           Fval_t(:),Fval_q(:),Fval_u(:),Fval_v(:)
  !           Eval_e(:),Eval_m(:)
  !
  !============================================================================

    ! Compute upward values for the implicit PBL scheme
    !-----------------------------------------------------
    dTa(:,pver) = FN_t(:) + EN_t(:)*dTs(:)
    dQa(:,pver) = FN_q(:) + EN_q(:)*dTs(:)
    dUa(:,pver) = FN_u(:)
    dVa(:,pver) = FN_v(:)
    do k=(pver-1),1,-1
      dTa(:,k) = Fval_t(:,k) + Eval_e(:,k)*dTa(:,k+1)
      dQa(:,k) = Fval_q(:,k) + Eval_e(:,k)*dQa(:,k+1)
      dUa(:,k) = Fval_u(:,k) + Eval_m(:,k)*dUa(:,k+1)
      dVa(:,k) = Fval_v(:,k) + Eval_m(:,k)*dVa(:,k+1)
    end do

    ! Update atmosphere values
    !--------------------------
    U(:,:) = U(:,:) + dUa(:,:)
    V(:,:) = V(:,:) + dVa(:,:)
    Q(:,:) = Q(:,:) + dQa(:,:)
    T(:,:) = T(:,:) + dTa(:,:)

    ! Return resulting Tendency values
    !----------------------------------
    dUa(:,:) = dUa(:,:)/dtime
    dVa(:,:) = dVa(:,:)/dtime
    dQa(:,:) = dQa(:,:)/dtime
    dTa(:,:) = dTa(:,:)/dtime

  end subroutine frierson_pbl
  !=======================================================================


  !=======================================================================
  subroutine frierson_pbl_USER(ncol, pver, dtime, pmid, pint, Zm, Zi,             &
                               Psfc, Tsfc, Qsfc, T, U, V, Q,  Fsw, Fdn,           &
                               Cdrag, Km, Ke, VSE, Z_pbl, Rf, dQa, dTa, dUa, dVa, &
                               LHflux, SHflux, TUflux, TVflux                     )
    !
    ! frierson_pbl_USER: This routine is a stub which users can use
    !                    to develop and test their own PBL scheme
    !==========================================================================
    !
    ! Passed Variables
    !------------------
    integer ,intent(in)   :: ncol                 ! Number of columns
    integer ,intent(in)   :: pver                 ! Number of levels
    real(r8),intent(in)   :: dtime                ! Time Step
    real(r8),intent(in)   :: pmid  (ncol,pver)    ! Pressure at model levels
    real(r8),intent(in)   :: pint  (ncol,pver+1)  ! Pressure at interface levels.
    real(r8),intent(in)   :: Zm    (ncol,pver)    ! Height at model levels.
    real(r8),intent(in)   :: Zi    (ncol,pver)    ! Height at interface levels.
    real(r8),intent(in)   :: Psfc  (ncol)         ! Surface Pressure.
    real(r8),intent(inout):: Tsfc  (ncol)         ! SST temperature K
    real(r8),intent(inout):: Qsfc  (ncol)         ! sea surface water vapor (kg/kg)
    real(r8),intent(inout):: T     (ncol,pver)    ! ATM Temperature values.
    real(r8),intent(inout):: U     (ncol,pver)    ! ATM Zonal Wind values.
    real(r8),intent(inout):: V     (ncol,pver)    ! ATM Meridional Wind values.
    real(r8),intent(inout):: Q     (ncol,pver)    ! ATM Water vapor values.
    real(r8),intent(in)   :: Fdn   (ncol)         ! Downward LW flux at surface
    real(r8),intent(in)   :: Fsw   (ncol)         ! Net SW flux at surface from gray radiation
    real(r8),intent(out)  :: Cdrag (ncol)         ! Surface drage coef.
    real(r8),intent(out)  :: Km    (ncol,pver+1)  ! Eddy diffusivity for PBL
    real(r8),intent(out)  :: Ke    (ncol,pver+1)  ! Eddy diffusivity for PBL
    real(r8),intent(out)  :: VSE   (ncol,pver)    ! Virtual-Dry Static energy.(huh?)
    real(r8),intent(out)  :: Z_pbl (ncol)         ! Height of PBL layer.
    real(r8),intent(out)  :: Rf    (ncol,pver)
    real(r8),intent(out)  :: dTa   (ncol,pver)
    real(r8),intent(out)  :: dQa   (ncol,pver)
    real(r8),intent(out)  :: dUa   (ncol,pver)
    real(r8),intent(out)  :: dVa   (ncol,pver)
    real(r8),intent(out)  :: LHflux(ncol)         ! Latent Heat Flux
    real(r8),intent(out)  :: SHflux(ncol)         ! Sensible Heat Flux
    real(r8),intent(out)  :: TUflux(ncol)         ! U Surface stress
    real(r8),intent(out)  :: TVflux(ncol)         ! V Surface stress
    !
    ! Local Values
    !---------------


    Cdrag  = 0._r8
    Km     = 0._r8
    Ke     = 0._r8
    VSE    = 0._r8
    Z_pbl  = 0._r8
    Rf     = 0._r8
    dTa    = 0._r8
    dQa    = 0._r8
    dUa    = 0._r8
    dVa    = 0._r8
    LHflux = 0._r8
    SHflux = 0._r8
    TUflux = 0._r8
    TVflux = 0._r8

  end subroutine frierson_pbl_USER
  !=======================================================================


  !=======================================================================
  subroutine frierson_radiation(ncol,pver,dtime,clat,pint,pmid,  &
                                Psfc,Tsfc,Qsfc,T,qv,dtdt_rad, &
                                Fsolar,Fup_s,Fdown_s,Fup_toa,Fdown_toa)
    !
    ! The gray radiation parameterization based on Frierson, et al. 2006.
    !
    ! frierson_radiation: This is an implementation of the gray radiation
    !                     scheme used in the Frierson model.
    !==========================================================================
    !
    ! Passed Variables
    !-------------------
    integer ,intent(in)   :: ncol                  ! number of columns
    integer ,intent(in)   :: pver                  ! number of vertical levels
    real(r8),intent(in)   :: dtime                 ! time step (s)
    real(r8),intent(in)   :: clat    (ncol)        ! latitude
    real(r8),intent(in)   :: pint    (ncol,pver+1) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pmid    (ncol,pver)   ! mid-point pressure (Pa)
    real(r8),intent(in)   :: Psfc    (ncol)        ! surface pressure
    real(r8),intent(in)   :: Tsfc    (ncol)        ! surface temperature (K)
    real(r8),intent(in)   :: Qsfc    (ncol)
    real(r8),intent(inout):: T       (ncol,pver)   ! temperature (K)
    real(r8),intent(in)   :: qv      (ncol,pver)   ! Q (kg/kg)
    real(r8),intent(out)  :: dtdt_rad(ncol,pver)   ! temperature tendency in K/s from relaxation
    real(r8),intent(out)  :: Fsolar  (ncol)        !
    real(r8),intent(out)  :: Fup_s   (ncol)        !
    real(r8),intent(out)  :: Fdown_toa(ncol)       !
    real(r8),intent(out)  :: Fup_toa  (ncol)       !
    real(r8),intent(out)  :: Fdown_s (ncol)        !
    !
    ! Local Values
    !-------------
    real(r8):: sinsq  (ncol)    ! sinlat**2
    real(r8):: Tv_srf (ncol)
    real(r8):: Tv     (ncol,pver  )
    real(r8):: Zm     (ncol,pver  )
    real(r8):: Zscl   (ncol)
    real(r8):: Tbar   (ncol)
    real(r8):: Tdif   (ncol)
    real(r8):: Zfrac  (ncol)
    real(r8):: Pfrac  (ncol)
    real(r8):: Tau_lat(ncol)
    real(r8):: Tau    (ncol,pver+1)
    real(r8):: Zi     (ncol,pver+1)
    real(r8):: Fup    (ncol,pver+1)
    real(r8):: Fdown  (ncol,pver+1)
    real(r8):: Bval   (ncol,pver)
    real(r8):: Etau   (ncol,pver)
    integer :: k

    ! Calc current Tv values
    !---------------------------------
    Tv_srf(:)  = Tsfc(:)*(1._r8+zvir*Qsfc(:))
    do k = 1, pver
      Tv(:,k) = T(:,k)*(1._r8+zvir*qv(:,k))
    end do

    ! Calc Geopotential Heights at model interface
    ! levels and at model layer levels
    !----------------------------------------------
    Zi(:,pver+1) = 0._r8
    Zm(:,pver  ) = rair*((Tv(:,pver)+Tv_srf(:))/2._r8)*log(Psfc(:)/pmid(:,pver))/gravit
    do k = pver-1,1,-1
      Zscl  (:) = rair*log(pmid(:,k+1)/pmid(:,k))/gravit
      Tbar  (:) = (Tv(:,k)+Tv(:,k+1))/2._r8
      Tdif  (:) = (Tv(:,k)-Tv(:,k+1))/2._r8
      Zfrac (:) = log(pint(:,k+1)/pmid(:,k+1))/log(pmid(:,k)/pmid(:,k+1))
      Zm(:,k  ) = Zm(:,k+1) + Zscl(:)*Tbar(:)
      Zi(:,k+1) = Zm(:,k+1) + Zscl(:)*((Tv(:,k+1)-2._r8*Tdif(:))*Zfrac(:)   &
                                       +                Tdif(:) *Zfrac(:)**2)
    end do
    Zfrac(:) = log(pint(:,1)/pmid(:,2))/log(pmid(:,1)/pmid(:,2))
    Zi(:,1)  = Zm(:,2) + Zscl(:)*((Tv(:,2)-2._r8*Tdif(:))*Zfrac(:)   &
                                       +         Tdif(:) *Zfrac(:)**2)

    ! Set Solar flux
    !------------------------
    sinsq (:) = sin(clat(:))*sin(clat(:))
    Fsolar(:) = (Rs0/4._r8)*(1._r8 + DeltaS*(1._r8 - 3._r8*sinsq(:))/4._r8)

    ! Calc optical depths
    !------------------------
    Tau_lat(:) = Tau_eqtr + (Tau_pole-Tau_eqtr)*sinsq(:)
    do k = 1,pver+1
      Pfrac(:) = pint(:,k)/Psfc(:)
      Tau(:,k) = Tau_lat(:)*(LinFrac*Pfrac(:) + (1._r8-LinFrac)*Pfrac(:)**4)
    end do

    ! Lowest order solution for up/down flux assumes B is constant for the layer
    !----------------------------------------------------------------------------
    do k=1,pver
      Bval(:,k) = Boltz*T(:,k)**4
      Etau(:,k) = exp(Tau(:,k)-Tau(:,k+1))
    end do

    Fup(:,pver+1) = Boltz*Tsfc(:)**4
    do k=pver,1,-1
      Fup(:,k) = Fup(:,k+1)*Etau(:,k) + Bval(:,k)*(1._r8-Etau(:,k))
    end do

    Fdown(:,1) = 0._r8
    do k=1,pver
      Fdown(:,k+1) = Fdown(:,k)*Etau(:,k) + Bval(:,k)*(1._r8-Etau(:,k))
    end do

    ! Calc Radiative heating
    !-------------------------
    do k=1,pver
      dtdt_rad(:,k) = -(cappa*Tv(:,k)/pmid(:,k))                           &
                      *((Fup(:,k+1)-Fdown(:,k+1)) - (Fup(:,k)-Fdown(:,k))) &
                       /(              Zi(:,k+1)  -   Zi(:,k)            )
    end do

    ! Return Upward/Downward long wave radiation at Surface and TOA
    !----------------------------------------------------------
    Fup_s  (:) = Fup  (:,pver+1)
    Fdown_s(:) = Fdown(:,pver+1)
    Fup_toa  (:) = Fup  (:,1)
    Fdown_toa(:) = Fdown(:,1)

    ! Update T values
    !-------------------
    do k=1,pver
      T(:,k) =  T(:,k) + dtdt_rad(:,k)*dtime
    end do

  end subroutine frierson_radiation
  !=======================================================================


  !=======================================================================
  subroutine frierson_radiation_USER(ncol,pver,dtime,clat,pint,pmid,  &
                                     Psfc,Tsfc,Qsfc,T,qv,dtdt_rad, &
                                     Fsolar,Fup_s,Fdown_s,Fup_toa,Fdown_toa)
    !
    ! frierson_radiation_USER: This routine is a stub which users can use
    !                          to develop and test their own radiation scheme
    !==========================================================================
    !
    ! Passed Variables
    !-------------------
    integer ,intent(in)   :: ncol                  ! number of columns
    integer ,intent(in)   :: pver                  ! number of vertical levels
    real(r8),intent(in)   :: dtime                 ! time step (s)
    real(r8),intent(in)   :: clat    (ncol)        ! latitude
    real(r8),intent(in)   :: pint    (ncol,pver+1) ! mid-point pressure (Pa)
    real(r8),intent(in)   :: pmid    (ncol,pver)   ! mid-point pressure (Pa)
    real(r8),intent(in)   :: Psfc    (ncol)        ! surface pressure
    real(r8),intent(in)   :: Tsfc    (ncol)        ! surface temperature (K)
    real(r8),intent(in)   :: Qsfc    (ncol)
    real(r8),intent(inout):: T       (ncol,pver)   ! temperature (K)
    real(r8),intent(in)   :: qv      (ncol,pver)   ! Q (kg/kg)
    real(r8),intent(out)  :: dtdt_rad(ncol,pver)   ! temperature tendency in K/s from relaxation
    real(r8),intent(out)  :: Fsolar  (ncol)        !
    real(r8),intent(out)  :: Fup_s   (ncol)        !
    real(r8),intent(out)  :: Fdown_s (ncol)        !
    real(r8),intent(out)  :: Fup_toa  (ncol)       !
    real(r8),intent(out)  :: Fdown_toa(ncol)      !
    !
    ! Local Values
    !-------------

    dtdt_rad = 0._r8
    Fsolar   = 0._r8
    Fup_s    = 0._r8
    Fdown_s  = 0._r8
    Fup_toa    = 0._r8
    Fdown_toa  = 0._r8


  end subroutine frierson_radiation_USER
  !=======================================================================

end module frierson
