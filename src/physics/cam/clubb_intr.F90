module clubb_intr

  !----------------------------------------------------------------------------------------------------- !
  ! Module to interface CAM with Cloud Layers Unified by Bi-normals (CLUBB), developed                   !
  !    by the University of Wisconsin Milwaukee Group (UWM).                                             !
  !                                                                                                      !
  ! CLUBB replaces the exisiting turbulence, shallow convection, and macrophysics in CAM5                !  
  !                                                                                                      !  
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      ! 
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! Authors:  P. Bogenschutz, C. Craig, A. Gettelman                                                     ! 
  ! Modified by: K Thayer-Calder                                                     ! 
  !                                                                                                      ! 
  !----------------------------------------------------------------------------------------------------- !

  use shr_kind_mod,        only: r8=>shr_kind_r8                                                                  
  use shr_const_mod,       only: shr_const_pi
  use ppgrid,              only: pver, pverp, pcols, begchunk, endchunk
  use phys_control,        only: phys_getopts
  use physconst,           only: cpair, gravit, rga, latvap, latice, zvir, rh2o, karman, rhoh2o
  use air_composition,     only: rairv, cpairv
  use cam_history_support, only: fillvalue, max_fieldname_len

  use spmd_utils,       only: masterproc 
  use constituents,     only: pcnst, cnst_add
  use pbl_utils,        only: calc_ustar, calc_obklen
  use ref_pres,         only: top_lev => trop_cloud_top_lev
  use zm_conv_intr,     only: zmconv_microp
#ifdef CLUBB_SGS
  use clubb_api_module,    only: pdf_parameter, implicit_coefs_terms
  use clubb_api_module,    only: clubb_config_flags_type, grid, stats, nu_vertical_res_dep
  use clubb_api_module,    only: nparams
  use clubb_mf,            only: do_clubb_mf, do_clubb_mf_diag
  use cloud_fraction,      only: dp1, dp2
  use clubb_mf,            only: do_clubb_mf, do_clubb_mf_diag, clubb_mf_nup, do_clubb_mf_rad, clubb_mf_Lopt, &
                                 clubb_mf_ddalph, clubb_mf_up_ndt, clubb_mf_cp_ndt, do_clubb_mf_cmt
  use cam_history_support, only: add_hist_coord
#endif

  implicit none
#ifdef CLUBB_SGS
  ! Variables that contains all the statistics

  type (stats), target, save :: stats_zt(pcols),      & ! stats_zt grid
                                stats_zm(pcols),      & ! stats_zm grid
                                stats_rad_zt(pcols),  & ! stats_rad_zt grid
                                stats_rad_zm(pcols),  & ! stats_rad_zm grid
                                stats_sfc(pcols)        ! stats_sfc
                                
!$omp threadprivate(stats_zt, stats_zm, stats_rad_zt, stats_rad_zm, stats_sfc)

#endif

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, clubb_emissions_cam, &
#ifdef CLUBB_SGS
            ! This utilizes CLUBB specific variables in its interface
            stats_init_clubb, &
            stats_zt, stats_zm, stats_sfc, &
            stats_rad_zt, stats_rad_zm, &
            stats_end_timestep_clubb, & 
#endif
            clubb_readnl, &
            clubb_init_cnst, &
            clubb_implements_cnst

#ifdef CLUBB_SGS
  ! Both of these utilize CLUBB specific variables in their interface
  private :: stats_zero, stats_avg
#endif

  logical, public :: do_cldcool
  logical         :: clubb_do_icesuper
#ifdef CLUBB_SGS
  type(clubb_config_flags_type), public :: clubb_config_flags
  real(r8), dimension(nparams), public :: clubb_params    ! Adjustable CLUBB parameters (C1, C2 ...)
#endif

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer, parameter :: &
      grid_type    = 3, &               ! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0                  ! The hydromet array in SAM-CLUBB is currently 0 elements
   
  ! Even though sclr_dim is set to 0, the dimension here is set to 1 to prevent compiler errors
  ! See github ticket larson-group/cam#133 for details
  real(r8), parameter, dimension(1) :: &
      sclr_tol = 1.e-8_r8               ! Total water in kg/kg

  character(len=6) :: saturation_equation

  real(r8), parameter :: &
      theta0   = 300._r8, &             ! Reference temperature                     [K]
      ts_nudge = 86400._r8, &           ! Time scale for u/v nudging (not used)     [s]
      p0_clubb = 100000._r8
      
  integer, parameter :: & 
    sclr_dim = 0                        ! Higher-order scalars, set to zero

  real(r8), parameter :: &
    wp3_const = 1._r8                   ! Constant to add to wp3 when moments are advected
    
  real(r8), parameter :: & 
    wpthlp_const = 10.0_r8              ! Constant to add to wpthlp when moments are advected
    
  real(r8), parameter :: & 
    wprtp_const = 0.01_r8               ! Constant to add to wprtp when moments are advected
    
  real(r8), parameter :: & 
    rtpthlp_const = 0.01_r8             ! Constant to add to rtpthlp when moments are advected
    
  real(r8), parameter :: unset_r8 = huge(1.0_r8)
  integer, parameter  :: unset_i = huge(1)

  ! Commonly used temperature for the melting temp of ice crystals [K] 
  real(r8), parameter :: meltpt_temp = 268.15_r8  
    
  real(r8) :: clubb_timestep = unset_r8  ! Default CLUBB timestep, unless overwriten by namelist
  real(r8) :: clubb_rnevap_effic = unset_r8

  real(r8) :: clubb_c1 = unset_r8
  real(r8) :: clubb_c1b = unset_r8
  real(r8) :: clubb_C2rt = unset_r8
  real(r8) :: clubb_C2thl = unset_r8
  real(r8) :: clubb_C2rtthl = unset_r8
  real(r8) :: clubb_C4 = unset_r8
  real(r8) :: clubb_C6rt = unset_r8
  real(r8) :: clubb_c6rtb = unset_r8
  real(r8) :: clubb_c6rtc = unset_r8
  real(r8) :: clubb_c6thl = unset_r8
  real(r8) :: clubb_c6thlb = unset_r8
  real(r8) :: clubb_c6thlc = unset_r8
  real(r8) :: clubb_C8 = unset_r8
  real(r8) :: clubb_C8b = unset_r8
  real(r8) :: clubb_C7 = unset_r8
  real(r8) :: clubb_C7b = unset_r8
  real(r8) :: clubb_c11 = unset_r8
  real(r8) :: clubb_c11b = unset_r8
  real(r8) :: clubb_c14 = unset_r8
  real(r8) :: clubb_C_wp3_pr_turb = unset_r8
  real(r8) :: clubb_c_K1 = unset_r8
  real(r8) :: clubb_c_K2 = unset_r8
  real(r8) :: clubb_nu2 = unset_r8
  real(r8) :: clubb_c_K8 = unset_r8
  real(r8) :: clubb_c_K9 = unset_r8
  real(r8) :: clubb_nu9 = unset_r8
  real(r8) :: clubb_c_K10 = unset_r8
  real(r8) :: clubb_c_K10h = unset_r8
  real(r8) :: clubb_C_invrs_tau_bkgnd = unset_r8
  real(r8) :: clubb_C_invrs_tau_sfc = unset_r8
  real(r8) :: clubb_C_invrs_tau_shear = unset_r8
  real(r8) :: clubb_C_invrs_tau_N2 = unset_r8
  real(r8) :: clubb_C_invrs_tau_N2_wp2 = unset_r8
  real(r8) :: clubb_C_invrs_tau_N2_xp2 = unset_r8
  real(r8) :: clubb_C_invrs_tau_N2_wpxp = unset_r8
  real(r8) :: clubb_C_invrs_tau_N2_clear_wp3 = unset_r8
  real(r8) :: clubb_C_uu_shr = unset_r8
  real(r8) :: clubb_C_uu_buoy = unset_r8
  real(r8) :: clubb_gamma_coef = unset_r8
  real(r8) :: clubb_gamma_coefb = unset_r8
  real(r8) :: clubb_beta = unset_r8
  real(r8) :: clubb_lambda0_stability_coef = unset_r8
  real(r8) :: clubb_lmin_coef = unset_r8
  real(r8) :: clubb_mult_coef = unset_r8
  real(r8) :: clubb_Skw_denom_coef = unset_r8
  real(r8) :: clubb_skw_max_mag = unset_r8
  real(r8) :: clubb_up2_sfc_coef = unset_r8
  real(r8) :: clubb_C_wp2_splat = unset_r8
  real(r8) :: clubb_wpxp_L_thresh = unset_r8
  real(r8) :: clubb_detliq_rad = unset_r8
  real(r8) :: clubb_detice_rad = unset_r8
  real(r8) :: clubb_detphase_lowtemp = unset_r8
  
  integer :: &
    clubb_iiPDF_type,          & ! Selected option for the two-component normal
                                 ! (double Gaussian) PDF type to use for the w, rt,
                                 ! and theta-l (or w, chi, and eta) portion of
                                 ! CLUBB's multivariate, two-component PDF.
    clubb_ipdf_call_placement = unset_i, & ! Selected option for the placement of the call to
                                           ! CLUBB's PDF.
    clubb_penta_solve_method = unset_i,  & ! Specifier for method to solve the penta-diagonal system
    clubb_tridiag_solve_method = unset_i   ! Specifier for method to solve tri-diagonal systems


   
  logical :: &
    clubb_l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                          ! precipitation fraction is automatically set to 1 when this
                                          ! flag is turned off.
    clubb_l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                          ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                          ! <w'thl'>, <sclr>, and <w'sclr'> in subroutine
                                          ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                          ! approximated by eddy diffusivity when <u> and <v> are
                                          ! advanced in subroutine advance_windm_edsclrm.
    clubb_l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                          ! the overall correlation of w and x (w and rt, as well as w
                                          ! and theta-l) within the limits of -max_mag_correlation_flux
                                          ! to max_mag_correlation_flux.
    clubb_l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                          ! thlp2) on keeping the overall correlation of w and x within
                                          ! the limits of -max_mag_correlation_flux to
                                          ! max_mag_correlation_flux.
    clubb_l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                          ! turbulent dissipation coefficient, C2.
    clubb_l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
    clubb_l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
    clubb_l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
    clubb_l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                          ! differencing approximation rather than a centered
                                          ! differencing for turbulent or mean advection terms. It
                                          ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                          ! sclrpthlp.
    clubb_l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                          ! differencing approximation rather than a centered
                                          ! differencing for turbulent or mean advection terms. It
                                          ! affects rtm, thlm, sclrm, um and vm.
    clubb_l_uv_nudge,                   & ! For wind speed nudging.
    clubb_l_rtm_nudge,                  & ! For rtm nudging
    clubb_l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                          ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
    clubb_l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                          ! compute the varibles that are output from high order
                                          ! closure
    clubb_l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                          ! thermodynamic-level variables output from pdf_closure.
    clubb_l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                          ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                          ! output from pdf_closure.
    clubb_l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                          ! subroutine pdf_closure twice.  If true, pdf_closure is
                                          ! called first on thermodynamic levels and then on momentum
                                          ! levels so that each variable is computed on its native
                                          ! level.  If false, pdf_closure is only called on
                                          ! thermodynamic levels, and variables which belong on
                                          ! momentum levels are interpolated.
    clubb_l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                          ! terms.  Setting to .false. means that a_1 and a_3 are
                                          ! pulled outside of the derivative in
                                          ! advance_wp2_wp3_module.F90 and in
                                          ! advance_xp2_xpyp_module.F90.
    clubb_l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                          ! than a centered discretization for the portion
                                          ! of the wp3 turbulent advection term for ADG1
                                          ! that is linearized in terms of wp3<t+1>.
                                          ! (Requires ADG1 PDF and clubb_l_standard_term_ta).
    clubb_l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                          ! differencing approximation rather than a centered 
                                          ! differencing for turbulent advection terms. 
                                          ! It affects  wpxp only.
    clubb_l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                          ! differencing approximation rather than a centered
                                          ! differencing for turbulent advection terms. It affects
                                          ! xpyp only.
    clubb_l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                          ! and rcm to help increase cloudiness at coarser grid
                                          ! resolutions.
    clubb_l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
    clubb_l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
    clubb_l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
    clubb_l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
    clubb_l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                          ! correction
    clubb_l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                          ! -(2/3)*em/tau_zm, as in Bougeault (1981)
    clubb_l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
    clubb_l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
    clubb_l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                          ! mixing length scale as Lscale = tau * tke
    clubb_l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
    clubb_l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
    clubb_l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
    clubb_l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                          ! saturated atmospheres (from Durran and Klemp, 1982)
    clubb_l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
    clubb_l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
    clubb_l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
    clubb_l_e3sm_config,                & ! Run model with E3SM settings
    clubb_l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                          ! a variable estimate of layer depth based on the depth
                                          ! over which wpthlp is positive near the ground when true
                                          ! More information can be found by
                                          ! Looking at issue #905 on the clubb repo
    clubb_l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
    clubb_l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
    clubb_l_smooth_Heaviside_tau_wpxp,  & ! Use smooth Heaviside 'Peskin' in computation of invrs_tau
    clubb_l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
    clubb_l_linearize_pbl_winds,        & ! Flag to turn on code to linearize PBL winds
    clubb_l_single_C2_Skw,              & ! Use a single Skewness dependent C2 for rtp2, thlp2, and
                                          ! rtpthlp
    clubb_l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
    clubb_l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
    clubb_l_update_pressure,            & ! Flag for having CLUBB update pressure and exner
    clubb_l_mono_flux_lim_thlm,         & ! Flag to turn on monotonic flux limiter for thlm
    clubb_l_mono_flux_lim_rtm,          & ! Flag to turn on monotonic flux limiter for rtm
    clubb_l_mono_flux_lim_um,           & ! Flag to turn on monotonic flux limiter for um
    clubb_l_mono_flux_lim_vm,           & ! Flag to turn on monotonic flux limiter for vm
    clubb_l_mono_flux_lim_spikefix,     & ! Flag to implement monotonic flux limiter code that
                                          ! eliminates spurious drying tendencies at model top  
    clubb_l_intr_sfc_flux_smooth = .false.! Add a locally calculated roughness to upwp and vpwp sfc fluxes

!  Constant parameters
  logical, parameter, private :: &
    l_implemented    = .true.,        &  ! Implemented in a host model (always true)
    l_host_applies_sfc_fluxes = .false.  ! Whether the host model applies the surface fluxes
    
  logical, parameter, private :: &
    apply_to_heat    = .false.           ! Apply WACCM energy fixer to heat or not (.true. = yes (duh))  

  logical            :: lq(pcnst)
  logical            :: prog_modal_aero
  logical            :: do_rainturb
  logical            :: clubb_do_adv
  logical            :: clubb_do_liqsupersat = .false.
  logical            :: clubb_do_energyfix   = .true.
  logical            :: history_budget
  logical            :: do_hb_above_clubb    = .false.
  integer            :: history_budget_histfile_num
  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB
  integer            :: offset
 
!  define physics buffer indicies here       
  integer :: &
    wp2_idx, &         	! vertical velocity variances
    wp3_idx, &         	! third moment of vertical velocity
    wpthlp_idx, &      	! turbulent flux of thetal
    wprtp_idx, &       	! turbulent flux of total water
    rtpthlp_idx, &     	! covariance of thetal and rt
    rtp2_idx, &        	! variance of total water
    thlp2_idx, &       	! variance of thetal
    rtp3_idx, &        	! total water 3rd order
    thlp3_idx, &       	! thetal 3rd order
    up2_idx, &         	! variance of east-west wind
    vp2_idx, &         	! variance of north-south wind
    up3_idx, &         	! east-west wind 3rd order
    vp3_idx, &         	! north-south wind 3rd order
    upwp_idx, &        	! east-west momentum flux
    vpwp_idx, &        	! north-south momentum flux
    thlm_idx, &        	! mean thetal
    rtm_idx, &         	! mean total water mixing ratio
    um_idx, &         	! mean of east-west wind
    vm_idx, &           ! mean of north-south wind
    wpthvp_idx, &       ! buoyancy flux
    wp2thvp_idx, &      ! second order buoyancy term
    rtpthvp_idx, &      ! moisture buoyancy correlation
    thlpthvp_idx, &     ! temperature buoyancy correlation
    sclrpthvp_idx, &    ! passive scalar buoyancy correlation
    wp2rtp_idx, &       ! w'^2 rt'
    wp2thlp_idx, &      ! w'^2 thl'
    uprcp_idx, &        ! < u' r_c' >
    vprcp_idx, &        ! < v' r_c' >
    rc_coef_idx, &      ! Coefficient of X'r_c' in Eq. (34)
    wp4_idx, &          ! w'^4
    wpup2_idx, &        ! w'u'^2
    wpvp2_idx, &        ! w'v'^2
    wp2up2_idx, &       ! w'^2 u'^2
    wp2vp2_idx, &       ! w'^2 v'^2
    cloud_frac_idx, &   ! CLUBB's cloud fraction
    cld_idx, &         	! Cloud fraction
    concld_idx, &       ! Convective cloud fraction
    ast_idx, &          ! Stratiform cloud fraction
    alst_idx, &         ! Liquid stratiform cloud fraction
    aist_idx, &         ! Ice stratiform cloud fraction
    qlst_idx, &         ! Physical in-cloud LWC
    qist_idx, &         ! Physical in-cloud IWC
    dp_frac_idx, &      ! deep convection cloud fraction
    sh_frac_idx, &      ! shallow convection cloud fraction
    kvh_idx, &		! CLUBB eddy diffusivity on thermo levels
    pblh_idx, &         ! PBL pbuf
    icwmrdp_idx, &	! In cloud mixing ratio for deep convection
    icwmrsh_idx, &      ! In cloud mixing ratio for shallow convection (MF)
    tke_idx, &          ! turbulent kinetic energy
    tpert_idx, &        ! temperature perturbation from PBL
    fice_idx, &         ! fice_idx index in physics buffer
    cmeliq_idx, &       ! cmeliq_idx index in physics buffer
    relvar_idx, &       ! relative cloud water variance
    accre_enhan_idx, &  ! optional accretion enhancement factor for MG
    npccn_idx, &        ! liquid ccn number concentration
    naai_idx, &         ! ice number concentration
    prer_evap_idx, &    ! rain evaporation rate
    qrl_idx, &          ! longwave cooling rate
    radf_idx, & 
    qsatfac_idx, &      ! subgrid cloud water saturation scaling factor 
    ice_supersat_idx, & ! ice cloud fraction for SILHS
    rcm_idx, &          ! Cloud water mixing ratio for SILHS
    ztodt_idx,&         ! physics timestep for SILHS
    clubbtop_idx        ! level index for CLUBB top

  ! Indices for microphysical covariance tendencies
  integer :: &
    rtp2_mc_zt_idx,   &
    thlp2_mc_zt_idx,  &
    wprtp_mc_zt_idx,  &
    wpthlp_mc_zt_idx, &
    rtpthlp_mc_zt_idx

  integer :: &          ! added pbuf fields for clubb to have restart bfb when ipdf_call_placement=2
    pdf_zm_w_1_idx, &
    pdf_zm_w_2_idx, &
    pdf_zm_varnce_w_1_idx, &
    pdf_zm_varnce_w_2_idx, &
    pdf_zm_mixt_frac_idx

  integer, public :: & 
    ixthlp2 = 0, &
    ixwpthlp = 0, &
    ixwprtp = 0, &
    ixwp2 = 0, &
    ixwp3 = 0, &
    ixrtpthlp = 0, &
    ixrtp2 = 0, &
    ixup2 = 0, &
    ixvp2 = 0

  integer :: cmfmc_sh_idx = 0

  integer :: &
    dlfzm_idx  = -1,    & ! ZM detrained convective cloud water mixing ratio.
    difzm_idx  = -1,    & ! ZM detrained convective cloud ice mixing ratio.
    dnlfzm_idx = -1,    & ! ZM detrained convective cloud water num concen.
    dnifzm_idx = -1       ! ZM detrained convective cloud ice num concen.

  integer :: &
    qtm_macmic1_idx, &
    qtm_macmic2_idx, &
    thlm_macmic1_idx, &
    thlm_macmic2_idx, &
    rcm_macmic_idx, &
    cldfrac_macmic_idx, &
    wpthlp_macmic_idx, &
    wprtp_macmic_idx, &
    wpthvp_macmic_idx, &
    mf_wpthlp_macmic_idx, &
    mf_wprtp_macmic_idx, &
    mf_wpthvp_macmic_idx, &
    up_macmic1_idx, &
    up_macmic2_idx, &
    dn_macmic1_idx, &
    dn_macmic2_idx, &
    upa_macmic1_idx, &
    upa_macmic2_idx, &
    dna_macmic1_idx, &
    dna_macmic2_idx, &
    thlu_macmic1_idx, &
    thlu_macmic2_idx, &
    qtu_macmic1_idx, &
    qtu_macmic2_idx, &
    thld_macmic1_idx, &
    thld_macmic2_idx, &
    qtd_macmic1_idx, &
    qtd_macmic2_idx, &
    dthl_macmic1_idx, &
    dthl_macmic2_idx, &
    dqt_macmic1_idx, &
    dqt_macmic2_idx, &
    dthlu_macmic1_idx, &
    dthlu_macmic2_idx, &
    dqtu_macmic1_idx, &
    dqtu_macmic2_idx, &
    dthld_macmic1_idx, &
    dthld_macmic2_idx, &
    dqtd_macmic1_idx, &
    dqtd_macmic2_idx, &
    ztop_macmic1_idx, &
    ztop_macmic2_idx, &
    ddcp_macmic1_idx, &
    ddcp_macmic2_idx

  integer :: &
    prec_sh_idx, &
    snow_sh_idx

  integer :: ztopmn_idx
  integer :: ztopma_idx
  integer :: ztopm1_macmic_idx
  integer :: ddcp_idx
  integer :: ddcp_macmic_idx
  integer :: ddcpmn_idx
  integer :: cbm1_idx
  integer :: cbm1_macmic_idx

  !  Output arrays for CLUBB statistics    
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90
  character(len=16)  :: deep_scheme      ! Default set in phys_control.F90
  character(len=16)  :: subcol_scheme

  integer, parameter :: ncnst=9
  character(len=8)   :: cnst_names(ncnst)
  logical            :: do_cnst=.false.

#ifdef CLUBB_SGS
  type(pdf_parameter), target, allocatable, public, protected :: &
                              pdf_params_chnk(:)    ! PDF parameters (thermo. levs.) [units vary]
                                              
  type(pdf_parameter), target, allocatable :: pdf_params_zm_chnk(:) ! PDF parameters on momentum levs. [units vary]
  
  type(implicit_coefs_terms), target, allocatable :: pdf_implicit_coefs_terms_chnk(:) ! PDF impl. coefs. & expl. terms      [units vary]
#endif

  contains
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_register_cam( )
!-------------------------------------------------------------------------------
! Description:
!   Register the constituents and fields in the physics buffer
! Author: P. Bogenschutz, C. Craig, A. Gettelman
! Modified: 7/2013 by K Thayer-Calder to include support for SILHS/subcolumns
!
!-------------------------------------------------------------------------------
#ifdef CLUBB_SGS

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    !  Add CLUBB fields to pbuf 
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dtype_i4, dyn_time_lvls
    use subcol_utils,    only: subcol_get_scheme

    integer :: cld_macmic_num_steps

    !----- Begin Code -----    
    call phys_getopts( eddy_scheme_out                 = eddy_scheme, &
                       deep_scheme_out                 = deep_scheme, & 
                       history_budget_out              = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num , &
                       do_hb_above_clubb_out           = do_hb_above_clubb, &
                       cld_macmic_num_steps_out        = cld_macmic_num_steps)

    subcol_scheme = subcol_get_scheme()

    if (trim(subcol_scheme) == 'SILHS') then
      saturation_equation = "flatau"
    else
      saturation_equation = "gfdl"       ! Goff & Gratch (1946) approximation for SVP
    end if

    if (clubb_do_adv) then
       cnst_names =(/'THLP2  ','RTP2   ','RTPTHLP','WPTHLP ','WPRTP  ','WP2    ','WP3    ','UP2    ','VP2    '/)
       do_cnst=.true.
       !  If CLUBB moments are advected, do not output them automatically which is typically done.  Some moments
       !    need a constant added to them before they are advected, thus this would corrupt the output.  
       !    Users should refer to the "XXXX_CLUBB" (THLP2_CLUBB for instance) output variables for these moments
       call cnst_add(trim(cnst_names(1)),0._r8,0._r8,0._r8,ixthlp2,longname='second moment vertical velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(2)),0._r8,0._r8,0._r8,ixrtp2,longname='second moment rtp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(3)),0._r8,0._r8,-999999._r8,ixrtpthlp,longname='covariance rtp thlp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(4)),0._r8,0._r8,-999999._r8,ixwpthlp,longname='CLUBB heat flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(5)),0._r8,0._r8,-999999._r8,ixwprtp,longname='CLUBB moisture flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(6)),0._r8,0._r8,0._r8,ixwp2,longname='CLUBB wp2',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(7)),0._r8,0._r8,-999999._r8,ixwp3,longname='CLUBB 3rd moment vert velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(8)),0._r8,0._r8,0._r8,ixup2,longname='CLUBB 2nd moment u wind',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(9)),0._r8,0._r8,0._r8,ixvp2,longname='CLUBB 2nd moment v wind',cam_outfld=.false.)
    end if
    if (do_hb_above_clubb) then
      call pbuf_add_field('clubbtop', 'physpkg', dtype_i4, (/pcols/), clubbtop_idx)
    endif

    !  put pbuf_add calls here (see macrop_driver.F90 for sample) use indicies defined at top
    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/),                    pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/),             tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/),             kvh_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/),                    tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    cld_idx)
    call pbuf_add_field('FICE',       'global',dtype_r8, (/pcols,pver/),               fice_idx)
    call pbuf_add_field('RAD_CLUBB',  'global', dtype_r8, (/pcols,pver/),               radf_idx)
    call pbuf_add_field('CMELIQ',     'global',dtype_r8, (/pcols,pver/),                  cmeliq_idx)
    call pbuf_add_field('QSATFAC',    'global',dtype_r8, (/pcols,pver/),                qsatfac_idx)

    call pbuf_add_field('WP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp2_idx)
    call pbuf_add_field('WP3_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp3_idx)
    call pbuf_add_field('WPTHLP_nadv',     'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wpthlp_idx)
    call pbuf_add_field('WPRTP_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wprtp_idx)
    call pbuf_add_field('RTPTHLP_nadv',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtpthlp_idx)
    call pbuf_add_field('RTP2_nadv',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtp2_idx)
    call pbuf_add_field('THLP2_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlp2_idx)
    call pbuf_add_field('UP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), up2_idx)
    call pbuf_add_field('VP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vp2_idx)    

    call pbuf_add_field('RTP3',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtp3_idx)
    call pbuf_add_field('THLP3',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlp3_idx)
    call pbuf_add_field('UP3',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), up3_idx)
    call pbuf_add_field('VP3',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vp3_idx)

    call pbuf_add_field('UPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), upwp_idx)
    call pbuf_add_field('VPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vpwp_idx)
    call pbuf_add_field('THLM',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlm_idx)
    call pbuf_add_field('RTM',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtm_idx)
    call pbuf_add_field('UM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), um_idx)
    call pbuf_add_field('VM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vm_idx)

    call pbuf_add_field('WPTHVP',     'global', dtype_r8, (/pcols,pverp/), wpthvp_idx)
    call pbuf_add_field('WP2THVP',    'global', dtype_r8, (/pcols,pverp/), wp2thvp_idx)
    call pbuf_add_field('RTPTHVP',    'global', dtype_r8, (/pcols,pverp/), rtpthvp_idx)
    call pbuf_add_field('THLPTHVP',   'global', dtype_r8, (/pcols,pverp/), thlpthvp_idx)
    call pbuf_add_field('CLOUD_FRAC', 'global', dtype_r8, (/pcols,pverp/), cloud_frac_idx)
    call pbuf_add_field('ISS_FRAC',   'global',  dtype_r8, (/pcols,pverp/), ice_supersat_idx)
    call pbuf_add_field('RCM',        'physpkg', dtype_r8, (/pcols,pverp/), rcm_idx)
    call pbuf_add_field('ZTODT',      'physpkg', dtype_r8, (/pcols/),       ztodt_idx)
    call pbuf_add_field('WP2RTP',     'global', dtype_r8, (/pcols,pverp/), wp2rtp_idx)
    call pbuf_add_field('WP2THLP',    'global', dtype_r8, (/pcols,pverp/), wp2thlp_idx)
    call pbuf_add_field('UPRCP',      'global', dtype_r8, (/pcols,pverp/), uprcp_idx)
    call pbuf_add_field('VPRCP',      'global', dtype_r8, (/pcols,pverp/), vprcp_idx)
    call pbuf_add_field('RC_COEF',    'global', dtype_r8, (/pcols,pverp/), rc_coef_idx)
    call pbuf_add_field('WP4',        'global', dtype_r8, (/pcols,pverp/), wp4_idx)
    call pbuf_add_field('WPUP2',      'global', dtype_r8, (/pcols,pverp/), wpup2_idx)
    call pbuf_add_field('WPVP2',      'global', dtype_r8, (/pcols,pverp/), wpvp2_idx)
    call pbuf_add_field('WP2UP2',     'global', dtype_r8, (/pcols,pverp/), wp2up2_idx)
    call pbuf_add_field('WP2VP2',     'global', dtype_r8, (/pcols,pverp/), wp2vp2_idx)

    ! For SILHS microphysical covariance contributions
    call pbuf_add_field('rtp2_mc_zt', 'global', dtype_r8, (/pcols,pverp/), rtp2_mc_zt_idx)
    call pbuf_add_field('thlp2_mc_zt','global', dtype_r8, (/pcols,pverp/), thlp2_mc_zt_idx)
    call pbuf_add_field('wprtp_mc_zt','global', dtype_r8, (/pcols,pverp/), wprtp_mc_zt_idx)
    call pbuf_add_field('wpthlp_mc_zt','global',dtype_r8, (/pcols,pverp/), wpthlp_mc_zt_idx)
    call pbuf_add_field('rtpthlp_mc_zt','global',dtype_r8,(/pcols,pverp/), rtpthlp_mc_zt_idx)

    call pbuf_add_field('pdf_zm_w_1',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_1_idx)
    call pbuf_add_field('pdf_zm_w_2',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_2_idx)
    call pbuf_add_field('pdf_zm_var_w_1', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_1_idx)
    call pbuf_add_field('pdf_zm_var_w_2', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_2_idx)
    call pbuf_add_field('pdf_zm_mixt_frac',  'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_mixt_frac_idx)

    ! these extrra coord vars don't seem to work for interpolate_output=.true.
    call add_hist_coord('ncyc', cld_macmic_num_steps, 'macro/micro cycle index')
    call add_hist_coord('nens', clubb_mf_nup, 'clubb+mf ensemble size')

    call pbuf_add_field('qtm_macmic1'           ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), qtm_macmic1_idx)
    call pbuf_add_field('qtm_macmic2'           ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), qtm_macmic2_idx)
    call pbuf_add_field('thlm_macmic1'       ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), thlm_macmic1_idx)
    call pbuf_add_field('thlm_macmic2'       ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), thlm_macmic2_idx)
    call pbuf_add_field('RCM_CLUBB_macmic'    ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), rcm_macmic_idx)
    call pbuf_add_field('CLDFRAC_CLUBB_macmic','physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), cldfrac_macmic_idx)
    call pbuf_add_field('WPTHLP_CLUBB_macmic' ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), wpthlp_macmic_idx)
    call pbuf_add_field('WPRTP_CLUBB_macmic'  ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), wpthvp_macmic_idx)
    call pbuf_add_field('WPTHVP_CLUBB_macmic' ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), wprtp_macmic_idx)
    if (do_clubb_mf) then
      call pbuf_add_field('edmf_thlflx_macmic' ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), mf_wpthlp_macmic_idx)
      call pbuf_add_field('edmf_qtflx_macmic'  ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), mf_wprtp_macmic_idx)
      call pbuf_add_field('edmf_thvflx_macmic' ,'physpkg',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), mf_wpthvp_macmic_idx)
      call pbuf_add_field('ZTOPMN'             ,'global' ,  dtype_r8, (/clubb_mf_up_ndt,pcols/), ztopmn_idx)
      call pbuf_add_field('ZTOPMA'             ,'global' ,  dtype_r8, (/pcols/), ztopma_idx)
      call pbuf_add_field('ZTOP_MACMIC'        ,'physpkg',  dtype_r8, (/pcols/), ztopm1_macmic_idx)
      call pbuf_add_field('DDCP'               ,'global' ,  dtype_r8, (/pcols/), ddcp_idx)
      call pbuf_add_field('DDCP_MACMIC'        ,'physpkg',  dtype_r8, (/pcols/), ddcp_macmic_idx)
      call pbuf_add_field('DDCPMN'             ,'global' ,  dtype_r8, (/clubb_mf_cp_ndt,pcols/), ddcpmn_idx)
      call pbuf_add_field('CBM1'               ,'global' ,  dtype_r8, (/pcols/), cbm1_idx)
      call pbuf_add_field('CBM1_MACMIC'        ,'physpkg',  dtype_r8, (/pcols/), cbm1_macmic_idx)
      call pbuf_add_field('up_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), up_macmic1_idx)
      call pbuf_add_field('up_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), up_macmic2_idx)
      call pbuf_add_field('dn_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dn_macmic1_idx)
      call pbuf_add_field('dn_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dn_macmic2_idx)
      call pbuf_add_field('upa_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), upa_macmic1_idx)
      call pbuf_add_field('upa_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), upa_macmic2_idx)
      call pbuf_add_field('dna_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dna_macmic1_idx)
      call pbuf_add_field('dna_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dna_macmic2_idx)
      call pbuf_add_field('thlu_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), thlu_macmic1_idx)
      call pbuf_add_field('thlu_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), thlu_macmic2_idx)
      call pbuf_add_field('qtu_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), qtu_macmic1_idx)
      call pbuf_add_field('qtu_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), qtu_macmic2_idx)
      call pbuf_add_field('thld_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), thld_macmic1_idx)
      call pbuf_add_field('thld_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), thld_macmic2_idx)
      call pbuf_add_field('qtd_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), qtd_macmic1_idx)
      call pbuf_add_field('qtd_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), qtd_macmic2_idx)
      call pbuf_add_field('dthl_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dthl_macmic1_idx)
      call pbuf_add_field('dthl_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dthl_macmic2_idx)
      call pbuf_add_field('dqt_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dqt_macmic1_idx)
      call pbuf_add_field('dqt_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps*clubb_mf_nup/), dqt_macmic2_idx)
      call pbuf_add_field('dthlu_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dthlu_macmic1_idx)
      call pbuf_add_field('dthlu_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dthlu_macmic2_idx)
      call pbuf_add_field('dqtu_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dqtu_macmic1_idx)
      call pbuf_add_field('dqtu_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dqtu_macmic2_idx)
      call pbuf_add_field('dthld_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dthld_macmic1_idx)
      call pbuf_add_field('dthld_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dthld_macmic2_idx)
      call pbuf_add_field('dqtd_macmic1' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dqtd_macmic1_idx)
      call pbuf_add_field('dqtd_macmic2' ,'global',  dtype_r8, (/pcols,pverp*cld_macmic_num_steps/), dqtd_macmic2_idx)
      call pbuf_add_field('ztop_macmic1' ,'global',  dtype_r8, (/pcols,cld_macmic_num_steps/), ztop_macmic1_idx)
      call pbuf_add_field('ztop_macmic2' ,'global',  dtype_r8, (/pcols,cld_macmic_num_steps/), ztop_macmic2_idx)
      call pbuf_add_field('ddcp_macmic1' ,'global',  dtype_r8, (/pcols,cld_macmic_num_steps/), ddcp_macmic1_idx)
      call pbuf_add_field('ddcp_macmic2' ,'global',  dtype_r8, (/pcols,cld_macmic_num_steps/), ddcp_macmic2_idx)
    end if

#endif 

  end subroutine clubb_register_cam
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

function clubb_implements_cnst(name)

  !----------------------------------------------------------------------------- !
  !                                                                              !
  ! Return true if specified constituent is implemented by this package          !
  !                                                                              !
  !----------------------------------------------------------------------------- !

   character(len=*), intent(in) :: name      ! constituent name
   logical :: clubb_implements_cnst     ! return value

   !-----------------------------------------------------------------------

   clubb_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function clubb_implements_cnst

 
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

subroutine clubb_init_cnst(name, latvals, lonvals, mask, q)
#ifdef CLUBB_SGS
    use clubb_api_module,        only: w_tol_sqd, rt_tol, thl_tol  
#endif

   !----------------------------------------------------------------------- !
   !                                                                        !
   ! Initialize the state if clubb_do_adv                                   !
   !                                                                        !
   !----------------------------------------------------------------------- !

   character(len=*), intent(in)  :: name       ! constituent name
   real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev

   !-----------------------------------------------------------------------
   integer :: k, nlev

#ifdef CLUBB_SGS
   if (clubb_do_adv) then
      nlev = size(q, 2)
      do k = 1, nlev
         if (trim(name) == trim(cnst_names(1))) then
            where(mask)
               q(:,k) = thl_tol**2
            end where
         end if
         if (trim(name) == trim(cnst_names(2))) then
            where(mask)
               q(:,k) = rt_tol**2
            end where
         end if
         if (trim(name) == trim(cnst_names(3))) then
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
         if (trim(name) == trim(cnst_names(4))) then
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
         if (trim(name) == trim(cnst_names(5))) then
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
         if (trim(name) == trim(cnst_names(6))) then
            where(mask)
               q(:,k) = w_tol_sqd
            end where
         end if
         if (trim(name) == trim(cnst_names(7))) then
            where(mask)
               q(:,k) = 0.0_r8
            end where
         end if
         if (trim(name) == trim(cnst_names(8))) then
            where(mask)
               q(:,k) = w_tol_sqd
            end where
         end if
         if (trim(name) == trim(cnst_names(9))) then
            where(mask)
               q(:,k) = w_tol_sqd
            end where
         end if
      end do
   end if
#endif

end subroutine clubb_init_cnst

  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_readnl(nlfile)

#ifdef CLUBB_SGS
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical, mpi_real8, &
                               mpi_integer
    use clubb_mf,        only: clubb_mf_readnl
    
    use clubb_api_module, only: &
      set_default_clubb_config_flags_api, & ! Procedure(s)
      initialize_clubb_config_flags_type_api, &
      l_stats, l_output_rad_files
#endif

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

#ifdef CLUBB_SGS

    character(len=*), parameter :: sub = 'clubb_readnl'

    logical :: clubb_history = .false., clubb_rad_history = .false.  ! Stats enabled (T/F) 
    logical :: clubb_cloudtop_cooling = .false., clubb_rainevap_turb = .false.

    integer :: iunit, read_status, ierr

    namelist /clubb_his_nl/ clubb_history, clubb_rad_history
    namelist /clubbpbl_diff_nl/ clubb_cloudtop_cooling, clubb_rainevap_turb, &
                                clubb_do_adv, clubb_timestep,  &
                                clubb_rnevap_effic,clubb_do_icesuper
    namelist /clubb_params_nl/ clubb_beta, &
         clubb_c1, &
         clubb_c1b, &
         clubb_c11, &
         clubb_c11b, &
         clubb_c14, &
         clubb_C2rt, &
         clubb_C2rtthl, &
         clubb_C2thl, &
         clubb_C4, &
         clubb_c6rt, &
         clubb_c6rtb, &
         clubb_c6rtc, &
         clubb_c6thl, &
         clubb_c6thlb, &
         clubb_c6thlc, &
         clubb_C7, &
         clubb_C7b, &
         clubb_C8, &
         clubb_C8b, &
         clubb_C_invrs_tau_bkgnd, &
         clubb_C_invrs_tau_sfc, &
         clubb_C_invrs_tau_shear, &
         clubb_C_invrs_tau_N2, &
         clubb_C_invrs_tau_N2_clear_wp3, &
         clubb_C_invrs_tau_N2_wp2, &
         clubb_C_invrs_tau_N2_wpxp, &
         clubb_C_invrs_tau_N2_xp2, &
         clubb_c_K1, &
         clubb_c_K10, & 
         clubb_c_K10h, &
         clubb_c_K2, &
         clubb_c_K8, &
         clubb_c_K9, &
         clubb_C_uu_shr, &
         clubb_C_uu_buoy, &
         clubb_C_wp2_splat, &
         clubb_C_wp3_pr_turb, &
         clubb_detice_rad, &
         clubb_detliq_rad, &
         clubb_detphase_lowtemp, &
         clubb_do_energyfix, &
         clubb_do_liqsupersat, &
         clubb_gamma_coef, &
         clubb_gamma_coefb, &
         clubb_ipdf_call_placement, &
         clubb_lambda0_stability_coef, &
         clubb_lmin_coef, &
         clubb_l_brunt_vaisala_freq_moist, &
         clubb_l_call_pdf_closure_twice, &
         clubb_l_damp_wp2_using_em, &
         clubb_l_damp_wp3_Skw_squared, &
         clubb_l_diag_Lscale_from_tau, &
         clubb_l_do_expldiff_rtm_thlm, &
         clubb_l_e3sm_config, &
         clubb_l_enable_relaxed_clipping, &
         clubb_l_godunov_upwind_wpxp_ta, &
         clubb_l_godunov_upwind_xpyp_ta, &
         clubb_l_intr_sfc_flux_smooth, &
         clubb_l_lmm_stepping, &
         clubb_l_lscale_plume_centered, &
         clubb_l_min_wp2_from_corr_wx, &
         clubb_l_min_xp2_from_corr_wx, &
         clubb_l_mono_flux_lim_rtm, &
         clubb_l_mono_flux_lim_spikefix, &
         clubb_l_mono_flux_lim_thlm, &
         clubb_l_mono_flux_lim_um, &
         clubb_l_mono_flux_lim_vm, &
         clubb_l_partial_upwind_wp3, &
         clubb_l_predict_upwp_vpwp, &
         clubb_l_rcm_supersat_adj, &
         clubb_l_smooth_Heaviside_tau_wpxp, &
         clubb_l_stability_correct_tau_zm, &
         clubb_l_standard_term_ta, &
         clubb_l_trapezoidal_rule_zm, &
         clubb_l_trapezoidal_rule_zt, &
         clubb_l_upwind_xpyp_ta, &
         clubb_l_use_C11_Richardson, &
         clubb_l_use_C7_Richardson, &
         clubb_l_use_cloud_cover, &
         clubb_l_use_shear_Richardson, &
         clubb_l_use_thvm_in_bv_freq, &
         clubb_l_use_tke_in_wp2_wp3_K_dfsn, &
         clubb_l_use_tke_in_wp3_pr_turb_term, &
         clubb_l_vary_convect_depth, &
         clubb_l_vert_avg_closure, &
         clubb_mult_coef, &
         clubb_nu2, &
         clubb_nu9, &
         clubb_penta_solve_method, &
         clubb_Skw_denom_coef, &
         clubb_skw_max_mag, &
         clubb_tridiag_solve_method, &
         clubb_up2_sfc_coef, &
         clubb_wpxp_L_thresh

    !----- Begin Code -----

    !  Determine if we want clubb_history to be output  
    clubb_history      = .false.   ! Initialize to false
    l_stats            = .false.   ! Initialize to false
    l_output_rad_files = .false.   ! Initialize to false
    do_cldcool         = .false.   ! Initialize to false
    do_rainturb        = .false.   ! Initialize to false
    
    ! Initialize namelist variables to clubb defaults
    call set_default_clubb_config_flags_api( clubb_iiPDF_type, & ! Out
                                             clubb_ipdf_call_placement, & ! Out
                                             clubb_penta_solve_method, & ! Out
                                             clubb_tridiag_solve_method, & ! Out
                                             clubb_l_use_precip_frac, & ! Out
                                             clubb_l_predict_upwp_vpwp, & ! Out
                                             clubb_l_min_wp2_from_corr_wx, & ! Out
                                             clubb_l_min_xp2_from_corr_wx, & ! Out
                                             clubb_l_C2_cloud_frac, & ! Out
                                             clubb_l_diffuse_rtm_and_thlm, & ! Out
                                             clubb_l_stability_correct_Kh_N2_zm, & ! Out
                                             clubb_l_calc_thlp2_rad, & ! Out
                                             clubb_l_upwind_xpyp_ta, & ! Out
                                             clubb_l_upwind_xm_ma, & ! Out
                                             clubb_l_uv_nudge, & ! Out
                                             clubb_l_rtm_nudge, & ! Out
                                             clubb_l_tke_aniso, & ! Out
                                             clubb_l_vert_avg_closure, & ! Out
                                             clubb_l_trapezoidal_rule_zt, & ! Out
                                             clubb_l_trapezoidal_rule_zm, & ! Out
                                             clubb_l_call_pdf_closure_twice, & ! Out
                                             clubb_l_standard_term_ta, & ! Out
                                             clubb_l_partial_upwind_wp3, & ! Out
                                             clubb_l_godunov_upwind_wpxp_ta, & ! Out
                                             clubb_l_godunov_upwind_xpyp_ta, & ! Out
                                             clubb_l_use_cloud_cover, & ! Out
                                             clubb_l_diagnose_correlations, & ! Out
                                             clubb_l_calc_w_corr, & ! Out
                                             clubb_l_const_Nc_in_cloud, & ! Out
                                             clubb_l_fix_w_chi_eta_correlations, & ! Out
                                             clubb_l_stability_correct_tau_zm, & ! Out
                                             clubb_l_damp_wp2_using_em, & ! Out
                                             clubb_l_do_expldiff_rtm_thlm, & ! Out
                                             clubb_l_Lscale_plume_centered, & ! Out
                                             clubb_l_diag_Lscale_from_tau, & ! Out
                                             clubb_l_use_C7_Richardson, & ! Out
                                             clubb_l_use_C11_Richardson, & ! Out
                                             clubb_l_use_shear_Richardson, & ! Out
                                             clubb_l_brunt_vaisala_freq_moist, & ! Out
                                             clubb_l_use_thvm_in_bv_freq, & ! Out
                                             clubb_l_rcm_supersat_adj, & ! Out
                                             clubb_l_damp_wp3_Skw_squared, & ! Out
                                             clubb_l_prescribed_avg_deltaz, & ! Out
                                             clubb_l_lmm_stepping, & ! Out
                                             clubb_l_e3sm_config, & ! Out
                                             clubb_l_vary_convect_depth, & ! Out
                                             clubb_l_use_tke_in_wp3_pr_turb_term, & ! Out
                                             clubb_l_use_tke_in_wp2_wp3_K_dfsn, & ! Out
                                             clubb_l_smooth_Heaviside_tau_wpxp, & ! Out
                                             clubb_l_enable_relaxed_clipping, & ! Out
                                             clubb_l_linearize_pbl_winds, & ! Out
                                             clubb_l_mono_flux_lim_thlm, & ! Out
                                             clubb_l_mono_flux_lim_rtm, & ! Out
                                             clubb_l_mono_flux_lim_um, & ! Out
                                             clubb_l_mono_flux_lim_vm, & ! Out
                                             clubb_l_mono_flux_lim_spikefix ) ! Out

    !  Call CLUBB+MF namelist
    call clubb_mf_readnl(nlfile)

    !  Read namelist to determine if CLUBB history should be called
    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' )

      call find_group_name(iunit, 'clubb_his_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_his_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist')
         end if
      end if

      call find_group_name(iunit, 'clubb_params_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_params_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist')
         end if
      else
         call endrun('clubb_readnl:  error reading namelist')
      end if

      call find_group_name(iunit, 'clubbpbl_diff_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubbpbl_diff_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist')
         end if
      end if

      close(unit=iunit)
      call freeunit(iunit)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(clubb_history,                1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_history")
    call mpi_bcast(clubb_rad_history,            1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_rad_history")
    call mpi_bcast(clubb_do_icesuper,            1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_do_icesuper")
    call mpi_bcast(clubb_cloudtop_cooling,       1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_cloudtop_cooling")
    call mpi_bcast(clubb_rainevap_turb,          1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_rainevap_turb")
    call mpi_bcast(clubb_do_adv,                 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_do_adv")
    call mpi_bcast(clubb_timestep,               1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_timestep")
    call mpi_bcast(clubb_rnevap_effic,           1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_rnevap_effic")

    call mpi_bcast(clubb_c1,                    1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c1")
    call mpi_bcast(clubb_c1b,                    1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c1b")
    call mpi_bcast(clubb_c11,                    1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c11")
    call mpi_bcast(clubb_c11b,                   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c11b")
    call mpi_bcast(clubb_c14,                    1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c14")
    call mpi_bcast(clubb_C_wp3_pr_turb,          1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_wp3_pr_turb")
    call mpi_bcast(clubb_c6rt,                   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6rt")
    call mpi_bcast(clubb_c6rtb,                  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6rtb")
    call mpi_bcast(clubb_c6rtc,                  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6rtc")
    call mpi_bcast(clubb_c6thl,                 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6thl")
    call mpi_bcast(clubb_c6thlb,                 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6thlb")
    call mpi_bcast(clubb_c6thlc,                 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c6thlc")
    call mpi_bcast(clubb_wpxp_L_thresh,          1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_wpxp_L_thresh")
    call mpi_bcast(clubb_mult_coef,              1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mult_coef")
    call mpi_bcast(clubb_gamma_coef,             1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_gamma_coef")
    call mpi_bcast(clubb_c_K10,                  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K10")
    call mpi_bcast(clubb_c_K10h,                  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K10h")    
    call mpi_bcast(clubb_beta,                   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_beta")
    call mpi_bcast(clubb_C2rt,                   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C2rt")
    call mpi_bcast(clubb_C2thl,                  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C2thl")
    call mpi_bcast(clubb_C2rtthl,                1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C2rtthl")
    call mpi_bcast(clubb_C8,                     1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C8")
    call mpi_bcast(clubb_C8b,                     1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C8b")
    call mpi_bcast(clubb_C7,                     1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C7")
    call mpi_bcast(clubb_C7b,                    1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C7b")
    call mpi_bcast(clubb_Skw_denom_coef,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_Skw_denom_coef")
    call mpi_bcast(clubb_C4,                     1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C4")
    call mpi_bcast(clubb_C_uu_shr,               1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_uu_shr")
    call mpi_bcast(clubb_C_uu_buoy,              1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_uu_buoy")
    call mpi_bcast(clubb_c_K1,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K1")
    call mpi_bcast(clubb_c_K2,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K2")
    call mpi_bcast(clubb_nu2,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_nu2")
    call mpi_bcast(clubb_c_K8,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K8")
    call mpi_bcast(clubb_c_K9,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_c_K9")
    call mpi_bcast(clubb_nu9,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_nu9")
    call mpi_bcast(clubb_C_wp2_splat,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_wp2_splat")
    call mpi_bcast(clubb_lambda0_stability_coef, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_lambda0_stability_coef")
    call mpi_bcast(clubb_l_lscale_plume_centered,1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_lscale_plume_centered")
    call mpi_bcast(clubb_do_liqsupersat,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_do_liqsupersat")
    call mpi_bcast(clubb_do_energyfix,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_do_energyfix")
    call mpi_bcast(clubb_C_invrs_tau_bkgnd,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_bkgnd")    
    call mpi_bcast(clubb_C_invrs_tau_sfc,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_sfc")    
    call mpi_bcast(clubb_C_invrs_tau_shear,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_shear")    
    call mpi_bcast(clubb_C_invrs_tau_N2,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_N2")    
    call mpi_bcast(clubb_C_invrs_tau_N2_wp2,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_N2_wp2")    
    call mpi_bcast(clubb_C_invrs_tau_N2_xp2,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_N2_xp2")    
    call mpi_bcast(clubb_C_invrs_tau_N2_wpxp,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_N2_wpxp")    
    call mpi_bcast(clubb_C_invrs_tau_N2_clear_wp3,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_C_invrs_tau_N2_clear_wp3")    
    call mpi_bcast(clubb_lmin_coef, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_lmin_coef")
    call mpi_bcast(clubb_skw_max_mag, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_skw_max_mag")
    call mpi_bcast(clubb_l_stability_correct_tau_zm, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_stability_correct_tau_zm")
    call mpi_bcast(clubb_gamma_coefb, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_gamma_coefb")
    call mpi_bcast(clubb_up2_sfc_coef, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_up2_sfc_coef")
    call mpi_bcast(clubb_detliq_rad, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_detliq_rad")
    call mpi_bcast(clubb_detice_rad, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_detice_rad")
    call mpi_bcast(clubb_detphase_lowtemp, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_detphase_lowtemp")

    call mpi_bcast(clubb_l_use_C7_Richardson,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_C7_Richardson")
    call mpi_bcast(clubb_l_use_C11_Richardson,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_C11_Richardson")
    call mpi_bcast(clubb_l_use_shear_Richardson,       1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_shear_Richardson")
    call mpi_bcast(clubb_l_brunt_vaisala_freq_moist,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_brunt_vaisala_freq_moist")
    call mpi_bcast(clubb_l_use_thvm_in_bv_freq,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_thvm_in_bv_freq")
    call mpi_bcast(clubb_l_rcm_supersat_adj,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_rcm_supersat_adj")
    call mpi_bcast(clubb_l_damp_wp3_Skw_squared,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_damp_wp3_Skw_squared")
    call mpi_bcast(clubb_l_predict_upwp_vpwp,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_predict_upwp_vpwp")
    call mpi_bcast(clubb_l_min_wp2_from_corr_wx,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_min_wp2_from_corr_wx")
    call mpi_bcast(clubb_l_min_xp2_from_corr_wx,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_min_xp2_from_corr_wx")
    call mpi_bcast(clubb_l_upwind_xpyp_ta,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_upwind_xpyp_ta")
    call mpi_bcast(clubb_l_godunov_upwind_wpxp_ta,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_godunov_upwind_wpxp_ta")
    call mpi_bcast(clubb_l_godunov_upwind_xpyp_ta,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_godunov_upwind_xpyp_ta")
    call mpi_bcast(clubb_l_vert_avg_closure,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_vert_avg_closure")
    call mpi_bcast(clubb_l_trapezoidal_rule_zt,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_trapezoidal_rule_zt")
    call mpi_bcast(clubb_l_trapezoidal_rule_zm,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_trapezoidal_rule_zm")
    call mpi_bcast(clubb_l_call_pdf_closure_twice,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_call_pdf_closure_twice")
    call mpi_bcast(clubb_l_use_cloud_cover,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_cloud_cover")
    call mpi_bcast(clubb_l_diag_Lscale_from_tau,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_diag_Lscale_from_tau")
    call mpi_bcast(clubb_l_damp_wp2_using_em,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_damp_wp2_using_em")
    call mpi_bcast(clubb_l_do_expldiff_rtm_thlm,      1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_do_expldiff_rtm_thlm")
    call mpi_bcast(clubb_l_lmm_stepping,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_lmm_stepping")
    call mpi_bcast(clubb_l_e3sm_config,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_e3sm_config")
    call mpi_bcast(clubb_l_enable_relaxed_clipping,       1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_enable_relaxed_clipping")
    call mpi_bcast(clubb_l_use_tke_in_wp3_pr_turb_term,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_tke_in_wp3_pr_turb_term")
    call mpi_bcast(clubb_l_use_tke_in_wp2_wp3_K_dfsn,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_tke_in_wp2_wp3_K_dfsn")
    call mpi_bcast(clubb_l_smooth_Heaviside_tau_wpxp,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_smooth_Heaviside_tau_wpxp")
    call mpi_bcast(clubb_ipdf_call_placement,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_ipdf_call_placement")
    call mpi_bcast(clubb_l_mono_flux_lim_thlm,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_mono_flux_lim_thlm")
    call mpi_bcast(clubb_l_mono_flux_lim_rtm,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_mono_flux_lim_rtm")
    call mpi_bcast(clubb_l_mono_flux_lim_um,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_mono_flux_lim_um")
    call mpi_bcast(clubb_l_mono_flux_lim_vm,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_mono_flux_lim_vm")
    call mpi_bcast(clubb_l_mono_flux_lim_spikefix,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_mono_flux_lim_spikefix")
    call mpi_bcast(clubb_penta_solve_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_penta_solve_method")
    call mpi_bcast(clubb_tridiag_solve_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_tridiag_solve_method")
    call mpi_bcast(clubb_l_intr_sfc_flux_smooth,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_intr_sfc_flux_smooth")
    call mpi_bcast(clubb_l_vary_convect_depth,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_vary_convect_depth")
    call mpi_bcast(clubb_l_standard_term_ta,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_standard_term_ta")
    call mpi_bcast(clubb_l_partial_upwind_wp3,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_partial_upwind_wp3")

    !  Overwrite defaults if they are true
    if (clubb_history) l_stats = .true.
    if (clubb_rad_history) l_output_rad_files = .true. 
    if (clubb_cloudtop_cooling) do_cldcool = .true.
    if (clubb_rainevap_turb) do_rainturb = .true.
    
    ! Check that all namelists have been set
    if(clubb_timestep == unset_r8) call endrun(sub//": FATAL: clubb_timestep is not set")
    if(clubb_rnevap_effic == unset_r8) call endrun(sub//": FATAL:clubb_rnevap_effic  is not set")

    if(clubb_c1 == unset_r8) call endrun(sub//": FATAL: clubb_c1 is not set")
    if(clubb_c1b == unset_r8) call endrun(sub//": FATAL: clubb_c1b is not set")
    if(clubb_C2rt == unset_r8) call endrun(sub//": FATAL: clubb_C2rt is not set")
    if(clubb_C2thl == unset_r8) call endrun(sub//": FATAL: clubb_C2thl is not set")
    if(clubb_C2rtthl == unset_r8) call endrun(sub//": FATAL: clubb_C2rtthl is not set")
    if(clubb_C4 == unset_r8) call endrun(sub//": FATAL: clubb_C4 is not set")
    if(clubb_C_uu_shr == unset_r8) call endrun(sub//": FATAL: clubb_C_uu_shr is not set")
    if(clubb_C_uu_buoy == unset_r8) call endrun(sub//": FATAL: clubb_C_uu_buoy is not set")
    if(clubb_c6rt == unset_r8) call endrun(sub//": FATAL: clubb_c6rt is not set")
    if(clubb_c6rtb == unset_r8) call endrun(sub//": FATAL: clubb_c6rtb is not set")
    if(clubb_c6rtc == unset_r8) call endrun(sub//": FATAL: clubb_c6rtc is not set")
    if(clubb_c6thl == unset_r8) call endrun(sub//": FATAL: clubb_c6thl is not set")
    if(clubb_c6thlb == unset_r8) call endrun(sub//": FATAL: clubb_c6thlb is not set")
    if(clubb_c6thlc == unset_r8) call endrun(sub//": FATAL: clubb_c6thlc is not set")
    if(clubb_wpxp_L_thresh == unset_r8) call endrun(sub//": FATAL: clubb_wpxp_L_thresh is not set")
    if(clubb_C8 == unset_r8) call endrun(sub//": FATAL: clubb_C8 is not set")
    if(clubb_C8b == unset_r8) call endrun(sub//": FATAL: clubb_C8b is not set")
    if(clubb_C7 == unset_r8) call endrun(sub//": FATAL: clubb_C7 is not set")
    if(clubb_C7b == unset_r8) call endrun(sub//": FATAL: clubb_C7b is not set")
    if(clubb_c11 == unset_r8) call endrun(sub//": FATAL: clubb_c11 is not set")
    if(clubb_c11b == unset_r8) call endrun(sub//": FATAL: clubb_c11b is not set")
    if(clubb_c14 == unset_r8) call endrun(sub//": FATAL: clubb_c14 is not set")
    if(clubb_C_wp3_pr_turb == unset_r8) call endrun(sub//": FATAL: clubb_C_wp3_pr_turb is not set")
    if(clubb_c_K1 == unset_r8) call endrun(sub//": FATAL: clubb_c_K1 is not set")
    if(clubb_c_K2 == unset_r8) call endrun(sub//": FATAL: clubb_c_K2 is not set")
    if(clubb_nu2 == unset_r8) call endrun(sub//": FATAL: clubb_nu2 is not set")
    if(clubb_c_K8 == unset_r8) call endrun(sub//": FATAL: clubb_c_K8 is not set")
    if(clubb_c_K9 == unset_r8) call endrun(sub//": FATAL: clubb_c_K9 is not set")
    if(clubb_nu9 == unset_r8) call endrun(sub//": FATAL: clubb_nu9 is not set")
    if(clubb_c_K10 == unset_r8) call endrun(sub//": FATAL: clubb_c_K10 is not set")
    if(clubb_c_K10h == unset_r8) call endrun(sub//": FATAL: clubb_c_K10h is not set")
    if(clubb_C_invrs_tau_bkgnd == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_bkgnd is not set")
    if(clubb_C_invrs_tau_sfc == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_sfc is not set")
    if(clubb_C_invrs_tau_shear == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_shear is not set")
    if(clubb_C_invrs_tau_N2 == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_N2 is not set")
    if(clubb_C_invrs_tau_N2_wp2 == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_N2_wp2 is not set")
    if(clubb_C_invrs_tau_N2_xp2 == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_N2_xp2 is not set")
    if(clubb_C_invrs_tau_N2_wpxp == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_N2_wpxp is not set")
    if(clubb_C_invrs_tau_N2_clear_wp3 == unset_r8) call endrun(sub//": FATAL: clubb_C_invrs_tau_N2_clear_wp3 is not set")
    if(clubb_gamma_coef == unset_r8) call endrun(sub//": FATAL: clubb_gamma_coef is not set")
    if(clubb_gamma_coefb == unset_r8) call endrun(sub//": FATAL: clubb_gamma_coefb is not set")
    if(clubb_beta == unset_r8) call endrun(sub//": FATAL: clubb_beta is not set")
    if(clubb_lambda0_stability_coef == unset_r8) call endrun(sub//": FATAL: clubb_lambda0_stability_coef is not set")
    if(clubb_lmin_coef == unset_r8) call endrun(sub//": FATAL: clubb_lmin_coef is not set")
    if(clubb_mult_coef == unset_r8) call endrun(sub//": FATAL: clubb_mult_coef is not set")
    if(clubb_Skw_denom_coef == unset_r8) call endrun(sub//": FATAL: clubb_Skw_denom_coef is not set")
    if(clubb_skw_max_mag == unset_r8) call endrun(sub//": FATAL: clubb_skw_max_mag is not set")
    if(clubb_up2_sfc_coef == unset_r8) call endrun(sub//": FATAL: clubb_up2_sfc_coef is not set")
    if(clubb_C_wp2_splat == unset_r8) call endrun(sub//": FATAL: clubb_C_wp2_splatis not set")
    if(clubb_detliq_rad == unset_r8) call endrun(sub//": FATAL: clubb_detliq_rad not set")
    if(clubb_detice_rad == unset_r8) call endrun(sub//": FATAL: clubb_detice_rad not set")
    if(clubb_ipdf_call_placement == unset_i) call endrun(sub//": FATAL: clubb_ipdf_call_placement not set")
    if(clubb_detphase_lowtemp == unset_r8) call endrun(sub//": FATAL: clubb_detphase_lowtemp not set")
    if(clubb_penta_solve_method == unset_i) call endrun(sub//": FATAL: clubb_penta_solve_method not set")
    if(clubb_tridiag_solve_method == unset_i) call endrun(sub//": FATAL: clubb_tridiag_solve_method not set")
    if(clubb_detphase_lowtemp >= meltpt_temp) & 
    call endrun(sub//": ERROR: clubb_detphase_lowtemp must be less than 268.15 K")

    call initialize_clubb_config_flags_type_api( clubb_iiPDF_type, & ! In
                                                 clubb_ipdf_call_placement, & ! In
                                                 clubb_penta_solve_method, & ! In
                                                 clubb_tridiag_solve_method, & ! In
                                                 clubb_l_use_precip_frac, & ! In
                                                 clubb_l_predict_upwp_vpwp, & ! In
                                                 clubb_l_min_wp2_from_corr_wx, & ! In
                                                 clubb_l_min_xp2_from_corr_wx, & ! In
                                                 clubb_l_C2_cloud_frac, & ! In
                                                 clubb_l_diffuse_rtm_and_thlm, & ! In
                                                 clubb_l_stability_correct_Kh_N2_zm, & ! In
                                                 clubb_l_calc_thlp2_rad, & ! In
                                                 clubb_l_upwind_xpyp_ta, & ! In
                                                 clubb_l_upwind_xm_ma, & ! In
                                                 clubb_l_uv_nudge, & ! In
                                                 clubb_l_rtm_nudge, & ! In
                                                 clubb_l_tke_aniso, & ! In
                                                 clubb_l_vert_avg_closure, & ! In
                                                 clubb_l_trapezoidal_rule_zt, & ! In
                                                 clubb_l_trapezoidal_rule_zm, & ! In
                                                 clubb_l_call_pdf_closure_twice, & ! In
                                                 clubb_l_standard_term_ta, & ! In
                                                 clubb_l_partial_upwind_wp3, & ! In
                                                 clubb_l_godunov_upwind_wpxp_ta, & ! In
                                                 clubb_l_godunov_upwind_xpyp_ta, & ! In
                                                 clubb_l_use_cloud_cover, & ! In
                                                 clubb_l_diagnose_correlations, & ! In
                                                 clubb_l_calc_w_corr, & ! In
                                                 clubb_l_const_Nc_in_cloud, & ! In
                                                 clubb_l_fix_w_chi_eta_correlations, & ! In
                                                 clubb_l_stability_correct_tau_zm, & ! In
                                                 clubb_l_damp_wp2_using_em, & ! In
                                                 clubb_l_do_expldiff_rtm_thlm, & ! In
                                                 clubb_l_Lscale_plume_centered, & ! In
                                                 clubb_l_diag_Lscale_from_tau, & ! In
                                                 clubb_l_use_C7_Richardson, & ! In
                                                 clubb_l_use_C11_Richardson, & ! In
                                                 clubb_l_use_shear_Richardson, & ! In
                                                 clubb_l_brunt_vaisala_freq_moist, & ! In
                                                 clubb_l_use_thvm_in_bv_freq, & ! In
                                                 clubb_l_rcm_supersat_adj, & ! In
                                                 clubb_l_damp_wp3_Skw_squared, & ! In
                                                 clubb_l_prescribed_avg_deltaz, & ! In
                                                 clubb_l_lmm_stepping, & ! In
                                                 clubb_l_e3sm_config, & ! In
                                                 clubb_l_vary_convect_depth, & ! In
                                                 clubb_l_use_tke_in_wp3_pr_turb_term, & ! In
                                                 clubb_l_use_tke_in_wp2_wp3_K_dfsn, & ! In 
                                                 clubb_l_smooth_Heaviside_tau_wpxp, & ! In
                                                 clubb_l_enable_relaxed_clipping, & ! In
                                                 clubb_l_linearize_pbl_winds, & ! In
                                                 clubb_l_mono_flux_lim_thlm, & ! In
                                                 clubb_l_mono_flux_lim_rtm, & ! In
                                                 clubb_l_mono_flux_lim_um, & ! In
                                                 clubb_l_mono_flux_lim_vm, & ! In
                                                 clubb_l_mono_flux_lim_spikefix, & ! In
                                                 clubb_config_flags ) ! Out

#endif
  end subroutine clubb_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf2d)
!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Modifications: K Thayer-Calder 2013 July and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------



#ifdef CLUBB_SGS

    !  From CAM libraries
    use cam_history,            only: addfld, add_default, horiz_only
    use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx
    use cam_abortutils,         only: endrun

    ! These are needed to set parameters
    use clubb_api_module, only: &
         core_rknd, em_min, &
         ilambda0_stability_coef, ic_K10, ic_K10h, iC7, iC7b, iC8, iC8b, iC11, iC11b, iC4, iC_uu_shr, iC_uu_buoy, &
         iC1, iC1b, iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, iup2_sfc_coef, iwpxp_L_thresh, &
         iC14, iC_wp3_pr_turb, igamma_coef, igamma_coefb, imult_coef, ilmin_coef, &
         iSkw_denom_coef, ibeta, iskw_max_mag, &
         iC_invrs_tau_bkgnd,iC_invrs_tau_sfc,iC_invrs_tau_shear,iC_invrs_tau_N2,iC_invrs_tau_N2_wp2, &
         iC_invrs_tau_N2_xp2,iC_invrs_tau_N2_wpxp,iC_invrs_tau_N2_clear_wp3,iC_uu_shr,iC_uu_buoy, &
         iC2rt, iC2thl, iC2rtthl, ic_K1, ic_K2, inu2, ic_K8, ic_K9, inu9, iC_wp2_splat, params_list

    use clubb_api_module, only: &
         print_clubb_config_flags_api, &
         setup_clubb_core_api, &
         init_pdf_params_api, &
         time_precision, &
         core_rknd, &
         set_clubb_debug_level_api, &
         clubb_fatal_error, &     ! Error code value to indicate a fatal error
         nparams, &
         set_default_parameters_api, &
         read_parameters_api, &
         l_stats, &
         l_stats_samp, &
         l_grads, &
         w_tol_sqd, &
         rt_tol, &
         thl_tol

    !  These are only needed if we're using a passive scalar
    use clubb_api_module, only: &
         iisclr_rt, &
         iisclr_thl, &
         iisclr_CO2, &
         iiedsclr_rt, &
         iiedsclr_thl, &
         iiedsclr_CO2

    use time_manager,           only: is_first_step
    use clubb_api_module,       only: hydromet_dim
    use constituents,           only: cnst_get_ind
    use phys_control,           only: phys_getopts
    use spmd_utils,             only: iam
    use cam_logfile,            only: iulog
#endif

    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CLUBB_SGS

    real(kind=time_precision) :: dum1, dum2, dum3
    
    ! The similar name to clubb_history is unfortunate...
    logical :: history_amwg, history_clubb

    integer :: cld_macmic_num_steps

    integer :: err_code                   ! Code for when CLUBB fails
    integer :: i, j, k, l                    ! Indices
    integer :: nmodes, nspec, m
    integer :: ixq, ixcldice, ixcldliq, ixnumliq, ixnumice
    integer :: lptr

    logical, parameter :: l_input_fields = .false. ! Always false for CAM-CLUBB.
    logical, parameter :: l_update_pressure = .false. ! Always false for CAM-CLUBB.

    integer :: nlev

    real(r8) :: &
      C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
      C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
      C7, C7b, C7c, C8, C8b, C10, &
      C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  &
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
      gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
      Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
      rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
      C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, a3_coef_min

    !----- Begin Code -----

    nlev = pver + 1 - top_lev

    if (core_rknd /= r8) then
      call endrun('clubb_ini_cam:  CLUBB library core_rknd must match CAM r8 and it does not')
    end if

    ! Allocate PDF parameters across columns and chunks
    allocate( &
       pdf_params_chnk(begchunk:endchunk),   &
       pdf_params_zm_chnk(begchunk:endchunk), &
       pdf_implicit_coefs_terms_chnk(begchunk:endchunk) )

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that  
    ! CLUBB does not transport aerosol consituents.  Therefore, need to 
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents) 
    ! ----------------------------------------------------------------- !

    call phys_getopts(prog_modal_aero_out=prog_modal_aero, &
                      history_amwg_out=history_amwg, &
                      history_clubb_out=history_clubb, &
                      do_hb_above_clubb_out=do_hb_above_clubb, &
                      cld_macmic_num_steps_out=cld_macmic_num_steps)

    !  Select variables to apply tendencies back to CAM
 
    ! Initialize all consituents to true to start
    lq(1:pcnst) = .true.
    edsclr_dim  = pcnst

    call cnst_get_ind('Q',ixq)
    call cnst_get_ind('NUMICE',ixnumice)
    call cnst_get_ind('NUMLIQ',ixnumliq)
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)

    if (prog_modal_aero) then
       ! Turn off modal aerosols and decrement edsclr_dim accordingly
       call rad_cnst_get_info(0, nmodes=nmodes)
 
       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, lptr)
          lq(lptr)=.false.
          edsclr_dim = edsclr_dim-1
 
          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             lq(lptr)=.false.
             edsclr_dim = edsclr_dim-1
          end do
       end do
 
       !  In addition, if running with MAM, droplet number is transported
       !  in dropmixnuc, therefore we do NOT want CLUBB to apply transport
       !  tendencies to avoid double counted.  Else, we apply tendencies.
       lq(ixnumliq) = .false.
       edsclr_dim = edsclr_dim-1
    endif

    ! ----------------------------------------------------------------- !
    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    ! ----------------------------------------------------------------- !
    call set_clubb_debug_level_api( 0 )

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !


    !  Defaults
    l_stats_samp = .false.
    l_grads = .false.

    !  Overwrite defaults if needbe     
    if (l_stats) l_stats_samp = .true.

    !  Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC 
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    icwmrsh_idx = pbuf_get_index('ICWMRSH')     ! In-cloud shallow convective mixing ratio (EDMF)
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')
    naai_idx        = pbuf_get_index('NAAI')
    npccn_idx       = pbuf_get_index('NPCCN')

    ! CLUBB+MF
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

    iisclr_rt  = -1
    iisclr_thl = -1
    iisclr_CO2 = -1

    iiedsclr_rt  = -1
    iiedsclr_thl = -1
    iiedsclr_CO2 = -1

    if (zmconv_microp) then
       dlfzm_idx = pbuf_get_index('DLFZM')
       difzm_idx = pbuf_get_index('DIFZM')
       dnlfzm_idx = pbuf_get_index('DNLFZM')
       dnifzm_idx = pbuf_get_index('DNIFZM')
    end if

    ! ----------------------------------------------------------------- !
    ! Define number of tracers for CLUBB to diffuse
    ! ----------------------------------------------------------------- !    
    
    if (clubb_l_do_expldiff_rtm_thlm) then
       offset = 2 ! diffuse temperature and moisture explicitly
       edsclr_dim = edsclr_dim + offset 
    endif
    
    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !
    
    !  Read in parameters for CLUBB.  Just read in default values
    call set_default_parameters_api( &
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
               C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
               lambda0_stability_coef, mult_coef, taumin, taumax, &
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, &
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
               Cx_min, Cx_max, Richardson_num_min, &
               Richardson_num_max, a3_coef_min )

    call read_parameters_api( -99, "", &
                              C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
                              C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
                              C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
                              C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
                              C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
                              C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                              c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
                              c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
                              slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
                              coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
                              gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
                              omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
                              lambda0_stability_coef, mult_coef, taumin, taumax, &
                              Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
                              Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                              thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
                              Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
                              altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
                              C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, &
                              C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
                              C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
                              C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
                              Cx_min, Cx_max, Richardson_num_min, &
                              Richardson_num_max, a3_coef_min, &
                              clubb_params )

    clubb_params(iC2rtthl) = clubb_C2rtthl
    clubb_params(iC8) = clubb_C8
    clubb_params(iC11) = clubb_c11
    clubb_params(iC11b) = clubb_c11b
    clubb_params(iC14) = clubb_c14
    clubb_params(iC_wp3_pr_turb) = clubb_C_wp3_pr_turb
    clubb_params(ic_K10) = clubb_c_K10
    clubb_params(imult_coef) = clubb_mult_coef
    clubb_params(iSkw_denom_coef) = clubb_Skw_denom_coef
    clubb_params(iC2rt) = clubb_C2rt
    clubb_params(iC2thl) = clubb_C2thl
    clubb_params(ibeta) = clubb_beta
    clubb_params(iC6rt) = clubb_c6rt
    clubb_params(iC6rtb) = clubb_c6rtb
    clubb_params(iC6rtc) = clubb_c6rtc
    clubb_params(iC6thl) = clubb_c6thl
    clubb_params(iC6thlb) = clubb_c6thlb
    clubb_params(iC6thlc) = clubb_c6thlc
    clubb_params(iwpxp_L_thresh) = clubb_wpxp_L_thresh
    clubb_params(iC7) = clubb_C7
    clubb_params(iC7b) = clubb_C7b
    clubb_params(igamma_coef) = clubb_gamma_coef
    clubb_params(ic_K10h) = clubb_c_K10h
    clubb_params(ilambda0_stability_coef) = clubb_lambda0_stability_coef
    clubb_params(ilmin_coef) = clubb_lmin_coef
    clubb_params(iC8b) = clubb_C8b
    clubb_params(iskw_max_mag) = clubb_skw_max_mag
    clubb_params(iC1)  = clubb_C1
    clubb_params(iC1b) = clubb_C1b
    clubb_params(igamma_coefb) = clubb_gamma_coefb
    clubb_params(iup2_sfc_coef) = clubb_up2_sfc_coef
    clubb_params(iC4) = clubb_C4
    clubb_params(iC_uu_shr) = clubb_C_uu_shr
    clubb_params(iC_uu_buoy) = clubb_C_uu_buoy
    clubb_params(ic_K1) = clubb_c_K1
    clubb_params(ic_K2) = clubb_c_K2
    clubb_params(inu2)  = clubb_nu2
    clubb_params(ic_K8) = clubb_c_K8
    clubb_params(ic_K9) = clubb_c_K9
    clubb_params(inu9)  = clubb_nu9
    clubb_params(iC_wp2_splat) = clubb_C_wp2_splat
    clubb_params(iC_invrs_tau_bkgnd) = clubb_C_invrs_tau_bkgnd
    clubb_params(iC_invrs_tau_sfc) = clubb_C_invrs_tau_sfc
    clubb_params(iC_invrs_tau_shear) = clubb_C_invrs_tau_shear
    clubb_params(iC_invrs_tau_N2) = clubb_C_invrs_tau_N2
    clubb_params(iC_invrs_tau_N2_wp2) = clubb_C_invrs_tau_N2_wp2
    clubb_params(iC_invrs_tau_N2_xp2) = clubb_C_invrs_tau_N2_xp2
    clubb_params(iC_invrs_tau_N2_wpxp) = clubb_C_invrs_tau_N2_wpxp
    clubb_params(iC_invrs_tau_N2_clear_wp3) = clubb_C_invrs_tau_N2_clear_wp3
   
    !  Set up CLUBB core.  Note that some of these inputs are overwritten
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.
!$OMP PARALLEL
    call setup_clubb_core_api &
         ( nlev+1, theta0, ts_nudge, &                                ! In
           hydromet_dim,  sclr_dim, &                                 ! In
           sclr_tol, edsclr_dim, clubb_params, &                      ! In
           l_host_applies_sfc_fluxes, &                               ! In
           saturation_equation, &                                     ! In
           l_input_fields, &                                          ! In
           clubb_config_flags%iiPDF_type, &                           ! In
           clubb_config_flags%ipdf_call_placement, &                  ! In
           clubb_config_flags%l_predict_upwp_vpwp, &                  ! In
           clubb_config_flags%l_min_xp2_from_corr_wx, &               ! In
           clubb_config_flags%l_prescribed_avg_deltaz, &              ! In
           clubb_config_flags%l_damp_wp2_using_em, &                  ! In
           clubb_config_flags%l_stability_correct_tau_zm, &           ! In
           clubb_config_flags%l_enable_relaxed_clipping, &            ! In
           clubb_config_flags%l_diag_Lscale_from_tau, &               ! In
           err_code )    ! Out

    if ( err_code == clubb_fatal_error ) then
       call endrun('clubb_ini_cam:  FATAL ERROR CALLING SETUP_CLUBB_CORE')
    end if
!$OMP END PARALLEL

    ! Print the list of CLUBB parameters
    if ( masterproc ) then
       do j = 1, nparams, 1
          write(iulog,*) params_list(j), " = ", clubb_params(j)
       enddo
    endif

    ! Print configurable CLUBB flags
    if ( masterproc ) then
       write(iulog,'(a,i0,a)') " CLUBB configurable flags "
       call print_clubb_config_flags_api( iulog, clubb_config_flags ) ! Intent(in)
    end if
    
    ! ----------------------------------------------------------------- !
    ! Add output fields for the history files
    ! ----------------------------------------------------------------- !

    !  These are default CLUBB output.  Not the higher order history budgets
    call addfld ('RHO_CLUBB',        (/ 'lev' /),  'A', 'kg/m3',    'Air Density')
    call addfld ('UP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Zonal Velocity Variance')
    call addfld ('VP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Meridional Velocity Variance')
    call addfld ('WP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Vertical Velocity Variance')
    call addfld ('WP2_ZT_CLUBB',     (/ 'lev' /),  'A', 'm2/s2',    'Vert Vel Variance on zt grid')
    call addfld ('UPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',    'Zonal Momentum Flux')
    call addfld ('VPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',    'Meridional Momentum Flux')
    call addfld ('WP3_CLUBB',        (/ 'lev' /),  'A', 'm3/s3',    'Third Moment Vertical Velocity')
    call addfld ('WPTHLP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',     'Heat Flux')
    call addfld ('WPRTP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',     'Moisture Flux')
    call addfld ('RTP2_CLUBB',       (/ 'ilev' /), 'A', 'kg^2/kg^2', 'Moisture Variance')
    call addfld ('RTP2_ZT_CLUBB',    (/ 'lev' /),  'A', 'kg^2/kg^2','Moisture Variance on zt grid')
    call addfld ('PDFP_RTP2_CLUBB',  (/ 'ilev' /), 'A', 'kg^2/kg^2','PDF Rtot Variance')
    call addfld ('THLP2_CLUBB',      (/ 'ilev' /), 'A', 'K^2',      'Temperature Variance')
    call addfld ('THLP2_ZT_CLUBB',   (/ 'lev' /),  'A', 'K^2',      'Temperature Variance on zt grid')
    call addfld ('RTPTHLP_CLUBB',    (/ 'ilev' /), 'A', 'K kg/kg',   'Temp. Moist. Covariance')
    call addfld ('RCM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water Mixing Ratio')
    call addfld ('RTM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Total Water Mixing Ratio')
    call addfld ('THLM_CLUBB',       (/ 'lev' /),  'A', 'K',         'Liquid Water Potential Temperature')
    call addfld ('WPRCP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',     'Liquid Water Flux')
    call addfld ('CLOUDFRAC_CLUBB',  (/ 'lev' /),  'A', 'fraction', 'Cloud Fraction')
    call addfld ('RCMINLAYER_CLUBB', (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water in Layer')
    call addfld ('CLOUDCOVER_CLUBB', (/ 'lev' /),  'A', 'fraction', 'Cloud Cover') 
    call addfld ('WPTHVP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',     'Buoyancy Flux')
    call addfld ('RVMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Water vapor tendency')
    call addfld ('STEND_CLUBB',      (/ 'lev' /),  'A', 'J/(kg s)', 'Static energy tendency')
    call addfld ('RCMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Liquid Water Tendency')
    call addfld ('RIMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Ice Tendency')
    call addfld ('UTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',   'U-wind Tendency')
    call addfld ('VTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',   'V-wind Tendency')
    call addfld ('ZT_CLUBB',         (/ 'lev' /),  'A', 'm',        'Thermodynamic Heights')
    call addfld ('ZM_CLUBB',         (/ 'ilev' /), 'A', 'm',        'Momentum Heights')     
    call addfld ('UM_CLUBB',         (/ 'lev' /),  'A', 'm/s',      'Zonal Wind')
    call addfld ('VM_CLUBB',         (/ 'lev' /),  'A', 'm/s',      'Meridional Wind')
    call addfld ('WM_ZT_CLUBB',      (/ 'lev' /),  'A', 'm/s',      'Vertical Velocity')
    call addfld ('PBLH',             horiz_only,   'A', 'm',        'PBL height')
    call addfld ('CLDST',            (/ 'lev' /),  'A', 'fraction', 'Stratus cloud fraction')
    call addfld ('ZMDLF',            (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained liquid water from ZM convection')
    call addfld ('TTENDICE',         (/ 'lev' /),  'A', 'K/s',      'T tendency from Ice Saturation Adjustment')
    call addfld ('QVTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'Q tendency from Ice Saturation Adjustment')
    call addfld ('QITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'CLDICE tendency from Ice Saturation Adjustment')
    call addfld ('NITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'NUMICE tendency from Ice Saturation Adjustment')


    call addfld ('QCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'CLDICE tendency from Ice Saturation Adjustment')
    call addfld ('NCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'NUMICE tendency from Ice Saturation Adjustment')
    call addfld ('FQTENDICE',        (/ 'lev' /),  'A', 'fraction', 'Frequency of Ice Saturation Adjustment')

    call addfld ('DPDLFLIQ',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained liquid water from deep convection')
    call addfld ('DPDLFICE',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice from deep convection')  
    call addfld ('DPDLFT',           (/ 'lev' /),  'A', 'K/s',      'T-tendency due to deep convective detrainment') 
    call addfld ('RELVAR',           (/ 'lev' /),  'A', '-',        'Relative cloud water variance')
    call addfld ('CLUBB_GRID_SIZE',  horiz_only,   'A', 'm',        'Horizontal grid box size seen by CLUBB')


    call addfld ('ZMDLFI',           (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice water from ZM convection')
    call addfld ('CONCLD',           (/ 'lev' /),  'A', 'fraction', 'Convective cloud cover')
    call addfld ('CMELIQ',           (/ 'lev' /),  'A', 'kg/kg/s',  'Rate of cond-evap of liq within the cloud')
    call addfld ('DETNLIQTND',       (/ 'lev' /),  'A', '1/kg/s',   'CLDNUM tendency in detrained water')

    call addfld ('QSATFAC',          (/ 'lev' /),  'A', '-', 'Subgrid cloud water saturation scaling factor')
    call addfld ('KVH_CLUBB',        (/ 'ilev' /), 'A', 'm2/s', 'CLUBB vertical diffusivity of heat/moisture on interface levels')
    call addfld ('ELEAK_CLUBB',      horiz_only,   'A', 'W/m2', 'CLUBB energy leak')
    call addfld ('TFIX_CLUBB',       horiz_only,   'A', 'K', 'Temperature increment to conserve energy')

    call addfld ('TKE_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2', 'CLUBB tke on interface levels')
    ! ---------------------------------------------------------------------------- !
    ! Below are for detailed analysis of EDMF Scheme                               !
    ! ---------------------------------------------------------------------------- !
    if (do_clubb_mf) then
      call addfld ( 'edmf_DRY_A'    , (/ 'ilev' /), 'A', 'fraction', 'Dry updraft area fraction (EDMF)' )
      call addfld ( 'edmf_MOIST_A'  , (/ 'ilev' /), 'A', 'fraction', 'Moist updraft area fraction (EDMF)' )
      call addfld ( 'edmf_DRY_W'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft vertical velocity (EDMF)' )
      call addfld ( 'edmf_MOIST_W'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft vertical velocity (EDMF)' )
      call addfld ( 'edmf_DRY_QT'   , (/ 'ilev' /), 'A', 'kg/kg'   , 'Dry updraft total water mixing ratio (EDMF)' )
      call addfld ( 'edmf_MOIST_QT' , (/ 'ilev' /), 'A', 'kg/kg'   , 'Moist updraft total water mixing ratio (EDMF)' )
      call addfld ( 'edmf_DRY_THL'  , (/ 'ilev' /), 'A', 'K'       , 'Dry updraft liquid-ice potential temperature (EDMF)' )
      call addfld ( 'edmf_MOIST_THL', (/ 'ilev' /), 'A', 'K'       , 'Moist updraft liquid-ice potential temperature (EDMF)' )
      call addfld ( 'edmf_DRY_U'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft zonal velocity (EDMF)' )
      call addfld ( 'edmf_MOIST_U'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft zonal velocity (EDMF)' )
      call addfld ( 'edmf_DRY_V'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft meridional velocity (EDMF)' )
      call addfld ( 'edmf_MOIST_V'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft meridional velocity (EDMF)' )
      call addfld ( 'edmf_MOIST_QC' , (/ 'ilev' /), 'A', 'kg/kg'   , 'Moist updraft condensate mixing ratio (EDMF)' )
      call addfld ( 'edmf_precc'    , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft precipitation rate (EDMF)' )
      call addfld ( 'edmf_S_AE'     , (/ 'ilev' /), 'A', 'fraction', '1 minus sum of a_i*w_i (EDMF)' )
      call addfld ( 'edmf_S_AW'     , (/ 'ilev' /), 'A', 'm/s'     , 'Sum of a_i*w_i (EDMF)' )
      call addfld ( 'edmf_S_AWW'    , (/ 'ilev' /), 'A', 'm2/s2'   , 'Sum of a_i*w_i*w_i (EDMF)' )
      call addfld ( 'edmf_S_AWTHL'  , (/ 'ilev' /), 'A', 'K m/s'   , 'Sum of a_i*w_i*thl_i (EDMF)' )
      call addfld ( 'edmf_S_AWQT'   , (/ 'ilev' /), 'A', 'kgm/kgs' , 'Sum of a_i*w_i*q_ti (EDMF)' )
      call addfld ( 'edmf_S_AWTH'   , (/ 'ilev' /), 'A', 'K m/s'   , 'Sum of a_i*w_i*th_i (EDMF)' )
      call addfld ( 'edmf_S_AWQV'   , (/ 'ilev' /), 'A', 'kgm/kgs' , 'Sum of a_i*w_i*q_vi (EDMF)' )
      call addfld ( 'edmf_S_AWU'    , (/ 'ilev' /), 'A', 'm2/s2'   , 'Sum of a_i*w_i*u_i (EDMF)' )
      call addfld ( 'edmf_S_AWV'    , (/ 'ilev' /), 'A', 'm2/s2'   , 'Sum of a_i*w_i*v_i (EDMF)' )
      call addfld ( 'edmf_thlforcup', (/ 'lev' /),  'A', 'K/s'     , 'thl updraft forcing (EDMF)' )
      call addfld ( 'edmf_qtforcup' , (/ 'lev' /),  'A', 'kg/kg/s' , 'qt updraft forcing (EDMF)' )
      call addfld ( 'edmf_thlforcdn', (/ 'lev' /),  'A', 'K/s'     , 'thl downdraft forcing (EDMF)' )
      call addfld ( 'edmf_qtforcdn' , (/ 'lev' /),  'A', 'kg/kg/s' , 'qt downdraft forcing (EDMF)' )
      call addfld ( 'edmf_thlforc'  , (/ 'lev' /),  'A', 'K/s'     , 'thl forcing (EDMF)' )
      call addfld ( 'edmf_qtforc'   , (/ 'lev' /),  'A', 'kg/kg/s' , 'qt forcing (EDMF)' )
      call addfld ( 'edmf_thlflxup' , (/ 'ilev' /), 'A', 'K m/s'    , 'thl updraft flux (EDMF)' )
      call addfld ( 'edmf_qtflxup'  , (/ 'ilev' /), 'A', 'kg/kg m/s', 'qt updraft flux (EDMF)' )
      call addfld ( 'edmf_thlflxdn' , (/ 'ilev' /), 'A', 'K m/s'    , 'thl downdraft flux (EDMF)' )
      call addfld ( 'edmf_qtflxdn'  , (/ 'ilev' /), 'A', 'kg/kg m/s', 'qt downdraft flux (EDMF)' )
      call addfld ( 'edmf_thlflx'   , (/ 'ilev' /), 'A', 'K m/s'    , 'thl flux (EDMF)' )
      call addfld ( 'edmf_qtflx'    , (/ 'ilev' /), 'A', 'kg/kg m/s', 'qt flux (EDMF)' )
      call addfld ( 'edmf_thvflx'   , (/ 'ilev' /), 'A', 'K m/s'    , 'thv flux (EDMF)' )
      call addfld ( 'edmf_uflxup'   , (/ 'ilev' /), 'A', 'm2/s2'    , 'u updraft flux (EDMF)' )
      call addfld ( 'edmf_vflxup'   , (/ 'ilev' /), 'A', 'm2/s2'    , 'v updraft flux (EDMF)' )
      call addfld ( 'edmf_uflxdn'   , (/ 'ilev' /), 'A', 'm2/s2'    , 'u downdraft flux (EDMF)' )
      call addfld ( 'edmf_vflxdn'   , (/ 'ilev' /), 'A', 'm2/s2'    , 'v downdraft flux (EDMF)' )
      call addfld ( 'edmf_uflx'     , (/ 'ilev' /), 'A', 'm2/s2'    , 'u flux (EDMF)' )
      call addfld ( 'edmf_vflx'     , (/ 'ilev' /), 'A', 'm2/s2'    , 'v flux (EDMF)' )
      call addfld ( 'edmf_sqtup'    , (/ 'lev' /), 'A', 'kg/kg/s' , 'Plume updraft microphysics tendency (EDMF)' )
      call addfld ( 'edmf_sqtdn'    , (/ 'lev' /), 'A', 'kg/kg/s' , 'Plume downdraft microphysics tendency (EDMF)' )
      call addfld ( 'edmf_rcm'      , (/ 'ilev' /), 'A', 'kg/kg'   , 'grid mean cloud (EDMF)' )
      call addfld ( 'edmf_cloudfrac', (/ 'lev' /),  'A', 'fraction', 'grid mean cloud fraction (EDMF)' )
      call addfld ( 'edmf_ent'      , (/ 'lev' /),  'A', '1/m'     , 'ensemble mean entrainment (EDMF)' )
      call addfld ( 'edmf_ztop'     ,  horiz_only,  'A', 'm'       , 'edmf ztop',       flag_xyfill=.True.)
      call addfld ( 'edmf_ddcp'     ,  horiz_only,  'A', 'm/s'     , 'edmf ddcp',       flag_xyfill=.True.)
      call addfld ( 'edmf_L0'       ,  horiz_only,  'A', 'm'       , 'edmf dynamic L0', flag_xyfill=.True.)
      call addfld ( 'edmf_freq'       ,  horiz_only,  'A', 'unitless', 'edmf frequency mf is active', flag_xyfill=.True.)
      call addfld ( 'edmf_cfl'      ,  horiz_only,  'A', 'unitless', 'max flux cfl number (EDMF)' )
      call addfld ( 'edmf_cape'     ,  horiz_only,  'A', 'J/kg'    , 'ensemble mean CAPE (EDMF)' )
      call addfld ( 'edmf_upa'      , (/ 'ilev', 'nens' /), 'A', 'fraction', 'Plume updraft area fraction (EDMF)' )
      call addfld ( 'edmf_upw'      , (/ 'ilev', 'nens' /), 'A', 'm/s'     , 'Plume updraft vertical velocity (EDMF)' )
      call addfld ( 'edmf_upmf'     , (/ 'ilev', 'nens' /), 'A', 'kg/m2/s' , 'Plume updraft mass flux (EDMF)' )
      call addfld ( 'edmf_upqt'     , (/ 'ilev', 'nens' /), 'A', 'kg/kg'   , 'Plume updraft total water mixing ratio (EDMF)' )
      call addfld ( 'edmf_upthl'    , (/ 'ilev', 'nens' /), 'A', 'K'       , 'Plume updraft liquid potential temperature (EDMF)' )
      call addfld ( 'edmf_upthv'    , (/ 'ilev', 'nens' /), 'A', 'm/s'     , 'Plume updraft virtual potential temperature (EDMF)' )
      call addfld ( 'edmf_upth'     , (/ 'ilev', 'nens' /), 'A', 'm/s'     , 'Plume updraft potential temperature (EDMF)' )
      call addfld ( 'edmf_upqc'     , (/ 'ilev', 'nens' /), 'A', 'kg/kg'   , 'Plume updraft condensate mixing ratio (EDMF)' )
      call addfld ( 'edmf_upent'    , (/ 'ilev', 'nens' /), 'A', '1/m'     , 'Plume updraft entrainment rate (EDMF)' )
      call addfld ( 'edmf_updet'    , (/ 'ilev', 'nens' /), 'A', '1/m'     , 'Plume updraft dettrainment rate (EDMF)' )
      call addfld ( 'edmf_upbuoy'   , (/  'ilev', 'nens' /), 'A', 'm/s2'   , 'Plume updraft buoyancy (EDMF)' )
      call addfld ( 'edmf_dnw'      , (/ 'ilev', 'nens' /), 'A', 'm/s'     , 'Plume downdraft vertical velocity (EDMF)' )
      call addfld ( 'edmf_dnthl'    , (/ 'ilev', 'nens' /), 'A', 'K'       , 'Plume downdraft liquid potential temperature (EDMF)' )
      call addfld ( 'edmf_dnqt'     , (/ 'ilev', 'nens' /), 'A', 'kg/kg'   , 'Plume downdraft total water mixing ratio (EDMF)' )
    end if 
    call addfld ('qtm_macmic1'           , (/ 'ilev', 'ncyc' /), 'A', 'kg/kg'   , 'QT at macro/micro substep')
    call addfld ('qtm_macmic2'           , (/ 'ilev', 'ncyc' /), 'A', 'kg/kg'   , 'QT at macro/micro substep')
    call addfld ('thlm_macmic1'       , (/ 'ilev', 'ncyc' /), 'A', 'K'       , 'THETAL at macro/micro substep')
    call addfld ('thlm_macmic2'       , (/ 'ilev', 'ncyc' /), 'A', 'K'       , 'THETAL at macro/micro substep')
    call addfld ('RCM_CLUBB_macmic'    , (/ 'ilev', 'ncyc' /), 'A', 'kg/kg'   , 'RCM CLUBB at macro/micro substep')
    call addfld ('CLDFRAC_CLUBB_macmic', (/ 'ilev', 'ncyc' /), 'A', 'fraction', 'CLDFRAC CLUBB at macro/micro substep')
    call addfld ('WPTHLP_CLUBB_macmic' , (/ 'ilev', 'ncyc' /), 'A', 'W/m2'    , 'Heat Flux at macro/micro substep')
    call addfld ('WPRTP_CLUBB_macmic'  , (/ 'ilev', 'ncyc' /), 'A', 'W/m2'    , 'Moisture Flux at macro/micro substep')
    call addfld ('WPTHVP_CLUBB_macmic' , (/ 'ilev', 'ncyc' /), 'A', 'W/m2'    , 'Buoyancy Flux at macro/micro substep')
    if (do_clubb_mf) then
      call addfld ( 'edmf_thlflx_macmic', (/ 'ilev', 'ncyc' /), 'A', 'K m/s'  , 'thl flux (EDMF) at macro/micro substep' )
      call addfld ( 'edmf_thvflx_macmic', (/ 'ilev', 'ncyc' /), 'A', 'K m/s'  , 'thv flux (EDMF) at macro/micro substep' )
      call addfld ( 'edmf_qtflx_macmic' , (/ 'ilev', 'ncyc' /), 'A', 'kg/kg m/s' , 'qt flux (EDMF) at macro/micro substep' )
      call addfld ( 'up_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'up' )
      call addfld ( 'up_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'up' )
      call addfld ( 'dn_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'dn' )
      call addfld ( 'dn_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'dn' )
      call addfld ( 'upa_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'fraction'  , 'upa' )
      call addfld ( 'upa_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'fraction'  , 'upa' )
      call addfld ( 'dna_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'fraction'  , 'dna' )
      call addfld ( 'dna_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'fraction'  , 'dna' )
      call addfld ( 'thlu_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'thl up' )
      call addfld ( 'thlu_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'thl up' )
      call addfld ( 'qtu_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'qt up' )
      call addfld ( 'qtu_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'qt up' )
      call addfld ( 'thld_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'thl dn' )
      call addfld ( 'thld_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'thl dn' )
      call addfld ( 'qtd_macmic1', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'qt dn' )
      call addfld ( 'qtd_macmic2', (/ 'ilev', 'nens', 'ncyc' /), 'I', 'm/s'  , 'qt dn' )
      call addfld ( 'dthl_macmic1', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl tend' )
      call addfld ( 'dthl_macmic2', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl tend' )
      call addfld ( 'dqt_macmic1', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt tend' )
      call addfld ( 'dqt_macmic2', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt tend' )
      call addfld ( 'dthlu_macmic1', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl up tend' )
      call addfld ( 'dthlu_macmic2', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl up tend' )
      call addfld ( 'dqtu_macmic1', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt up tend' )
      call addfld ( 'dqtu_macmic2', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt up tend' )
      call addfld ( 'dthld_macmic1', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl dn tend' )
      call addfld ( 'dthld_macmic2', (/ 'ilev', 'ncyc' /), 'I', 'm/s'  , 'thl dn tend' )
      call addfld ( 'dqtd_macmic1', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt dn tend' )
      call addfld ( 'dqtd_macmic2', (/ 'ilev',  'ncyc' /), 'I', 'm/s'  , 'qt dn tend' )
      call addfld ( 'ztop_macmic1', (/ 'ncyc' /), 'I', 'm/s'  , 'ztop' )
      call addfld ( 'ztop_macmic2', (/ 'ncyc' /), 'I', 'm/s'  , 'ztop' )
      call addfld ( 'ddcp_macmic1', (/ 'ncyc' /), 'I', 'm/s'  , 'ddcp' )
      call addfld ( 'ddcp_macmic2', (/ 'ncyc' /), 'I', 'm/s'  , 'ddcp' )
    end if

    !  Initialize statistics, below are dummy variables
    dum1 = 300._r8
    dum2 = 1200._r8
    dum3 = 300._r8

    if (l_stats) then
      
      do i=1, pcols
        call stats_init_clubb( .true., dum1, dum2, &
                               nlev+1, nlev+1, nlev+1, dum3, &
                               stats_zt(i), stats_zm(i), stats_sfc(i), &
                               stats_rad_zt(i), stats_rad_zm(i))
      end do             

       allocate(out_zt(pcols,pverp,stats_zt(1)%num_output_fields))
       allocate(out_zm(pcols,pverp,stats_zm(1)%num_output_fields))
       allocate(out_sfc(pcols,1,stats_sfc(1)%num_output_fields))

       allocate(out_radzt(pcols,pverp,stats_rad_zt(1)%num_output_fields))
       allocate(out_radzm(pcols,pverp,stats_rad_zm(1)%num_output_fields))

    endif
  
    ! ----------------------------------------------------------------- !
    ! Make all of this output default, this is not CLUBB history
    ! ----------------------------------------------------------------- !
  
    if (clubb_do_adv .or. history_clubb) then 
       call add_default('RELVAR',           1, ' ')
       call add_default('RHO_CLUBB',        1, ' ')
       call add_default('UP2_CLUBB',        1, ' ')
       call add_default('VP2_CLUBB',        1, ' ')
       call add_default('WP2_CLUBB',        1, ' ')
       call add_default('WP2_ZT_CLUBB',     1, ' ')
       call add_default('WP3_CLUBB',        1, ' ')
       call add_default('UPWP_CLUBB',       1, ' ')
       call add_default('VPWP_CLUBB',       1, ' ')
       call add_default('WPTHLP_CLUBB',     1, ' ')
       call add_default('WPRTP_CLUBB',      1, ' ')
       call add_default('RTP2_CLUBB',       1, ' ')
       call add_default('RTP2_ZT_CLUBB',    1, ' ')
       call add_default('PDFP_RTP2_CLUBB',  1, ' ')
       call add_default('THLP2_CLUBB',      1, ' ')
       call add_default('THLP2_ZT_CLUBB',   1, ' ')
       call add_default('RTPTHLP_CLUBB',    1, ' ')
       call add_default('RCM_CLUBB',        1, ' ')
       call add_default('RTM_CLUBB',        1, ' ')
       call add_default('THLM_CLUBB',       1, ' ')
       call add_default('WPRCP_CLUBB',      1, ' ')
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('RCMINLAYER_CLUBB', 1, ' ')
       call add_default('CLOUDCOVER_CLUBB', 1, ' ')
       call add_default('WPTHVP_CLUBB',     1, ' ')
       call add_default('RVMTEND_CLUBB',    1, ' ')
       call add_default('STEND_CLUBB',      1, ' ')
       call add_default('RCMTEND_CLUBB',    1, ' ')
       call add_default('RIMTEND_CLUBB',    1, ' ')
       call add_default('UTEND_CLUBB',      1, ' ')
       call add_default('VTEND_CLUBB',      1, ' ')
       call add_default('ZT_CLUBB',         1, ' ')
       call add_default('ZM_CLUBB',         1, ' ')   
       call add_default('UM_CLUBB',         1, ' ')
       call add_default('VM_CLUBB',         1, ' ')
       call add_default('WM_ZT_CLUBB',      1, ' ')
       call add_default('PBLH',             1, ' ')
       call add_default('CONCLD',           1, ' ')
       call add_default('TKE_CLUBB',        1, ' ')
       call add_default('ELEAK_CLUBB',      1, ' ')
       call add_default('TFIX_CLUBB',       1, ' ')
    end if
      
     if (history_amwg) then
        call add_default('PBLH',             1, ' ')
     end if

    if (history_clubb) then

       call add_default('RELVAR',           1, ' ')
       call add_default('RHO_CLUBB',        1, ' ')
       call add_default('UPWP_CLUBB',       1, ' ')
       call add_default('VPWP_CLUBB',       1, ' ')
       call add_default('RCM_CLUBB',        1, ' ')
       call add_default('WPRCP_CLUBB',      1, ' ')
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('RCMINLAYER_CLUBB', 1, ' ')
       call add_default('CLOUDCOVER_CLUBB', 1, ' ')
       call add_default('WPTHVP_CLUBB',     1, ' ')
       call add_default('RVMTEND_CLUBB',    1, ' ')
       call add_default('STEND_CLUBB',      1, ' ')
       call add_default('RCMTEND_CLUBB',    1, ' ')
       call add_default('RIMTEND_CLUBB',    1, ' ')
       call add_default('UTEND_CLUBB',      1, ' ')
       call add_default('VTEND_CLUBB',      1, ' ')
       call add_default('ZT_CLUBB',         1, ' ')
       call add_default('ZM_CLUBB',         1, ' ')   
       call add_default('UM_CLUBB',         1, ' ')
       call add_default('VM_CLUBB',         1, ' ')
       !call add_default('SL',               1, ' ')
       !call add_default('QT',               1, ' ')
       call add_default('CONCLD',           1, ' ')

       if (do_clubb_mf_diag) then
         call add_default( 'edmf_DRY_A'    , 1, ' ')
         call add_default( 'edmf_MOIST_A'  , 1, ' ')
         call add_default( 'edmf_DRY_W'    , 1, ' ')
         call add_default( 'edmf_MOIST_W'  , 1, ' ')
         call add_default( 'edmf_DRY_QT'   , 1, ' ')
         call add_default( 'edmf_MOIST_QT' , 1, ' ')
         call add_default( 'edmf_DRY_THL'  , 1, ' ')
         call add_default( 'edmf_MOIST_THL', 1, ' ')
         call add_default( 'edmf_DRY_U'    , 1, ' ')
         call add_default( 'edmf_MOIST_U'  , 1, ' ')
         call add_default( 'edmf_DRY_V'    , 1, ' ')
         call add_default( 'edmf_MOIST_V'  , 1, ' ')
         call add_default( 'edmf_MOIST_QC' , 1, ' ')
         call add_default( 'edmf_precc'    , 1, ' ')
         call add_default( 'edmf_S_AE'     , 1, ' ')
         call add_default( 'edmf_S_AW'     , 1, ' ')
         call add_default( 'edmf_S_AWW'    , 1, ' ')
         call add_default( 'edmf_S_AWTH'   , 1, ' ')
         call add_default( 'edmf_S_AWTHL'  , 1, ' ')
         call add_default( 'edmf_S_AWQT'   , 1, ' ')
         call add_default( 'edmf_S_AWU'    , 1, ' ')
         call add_default( 'edmf_S_AWV'    , 1, ' ')
         call add_default( 'edmf_thlflxup' , 1, ' ')
         call add_default( 'edmf_qtflxup'  , 1, ' ')
         call add_default( 'edmf_thlflxdn' , 1, ' ')
         call add_default( 'edmf_qtflxdn'  , 1, ' ')
         call add_default( 'edmf_thlflx'   , 1, ' ')
         call add_default( 'edmf_thvflx'   , 1, ' ')
         call add_default( 'edmf_uflxup'   , 1, ' ')
         call add_default( 'edmf_vflxup'   , 1, ' ')
         call add_default( 'edmf_uflxdn'   , 1, ' ')
         call add_default( 'edmf_vflxdn'   , 1, ' ')
         call add_default( 'edmf_uflx'     , 1, ' ')
         call add_default( 'edmf_vflx'     , 1, ' ')
         call add_default( 'edmf_qtflx'    , 1, ' ')

         call add_default( 'edmf_thlforcup', 1, ' ')
         call add_default( 'edmf_qtforcup' , 1, ' ')
         call add_default( 'edmf_thlforcdn', 1, ' ')
         call add_default( 'edmf_qtforcdn' , 1, ' ')

         call add_default( 'edmf_thlforc'  , 1, ' ')
         call add_default( 'edmf_qtforc'   , 1, ' ')
         call add_default( 'edmf_sqtup'    , 1, ' ')
         call add_default( 'edmf_sqtdn'    , 1, ' ')
         call add_default( 'edmf_rcm'      , 1, ' ')
         call add_default( 'edmf_cloudfrac', 1, ' ')
         call add_default( 'edmf_ent'      , 1, ' ')
         call add_default( 'edmf_ztop'     , 1, ' ')
         call add_default( 'edmf_ddcp'     , 1, ' ')
         call add_default( 'edmf_L0'       , 1, ' ')
         call add_default( 'edmf_freq'       , 1, ' ')
         call add_default( 'edmf_cape'     , 1, ' ')
         call add_default( 'edmf_cfl'     , 1, ' ')
       end if
       !call add_default( 'qtm_macmic1'           , 1, ' ')
       !call add_default( 'qtm_macmic2'           , 1, ' ')
       !call add_default( 'thlm_macmic1'       , 1, ' ')
       !call add_default( 'thlm_macmic2'       , 1, ' ')
       call add_default( 'RCM_CLUBB_macmic'    , 1, ' ')
       call add_default( 'CLDFRAC_CLUBB_macmic', 1, ' ')
       call add_default( 'WPTHLP_CLUBB_macmic' , 1, ' ')
       call add_default( 'WPRTP_CLUBB_macmic'  , 1, ' ')
       call add_default( 'WPTHVP_CLUBB_macmic' , 1, ' ')
       if (do_clubb_mf_diag) then
         call add_default( 'edmf_thlflx_macmic' , 1, ' ')
         call add_default( 'edmf_qtflx_macmic'  , 1, ' ')
         call add_default( 'edmf_thvflx_macmic' , 1, ' ')
         !call add_default( 'up_macmic1' , 1, ' ')
         !call add_default( 'up_macmic2' , 1, ' ')
         !call add_default( 'dn_macmic1' , 1, ' ')
         !call add_default( 'dn_macmic2' , 1, ' ')
         !call add_default( 'upa_macmic1' , 1, ' ')
         !call add_default( 'upa_macmic2' , 1, ' ')
         !call add_default( 'dna_macmic1' , 1, ' ')
         !call add_default( 'dna_macmic2' , 1, ' ')
         !call add_default( 'thlu_macmic1' , 1, ' ')
         !call add_default( 'thlu_macmic2' , 1, ' ')
         !call add_default( 'qtu_macmic1' , 1, ' ')
         !call add_default( 'qtu_macmic2' , 1, ' ')
         !call add_default( 'thld_macmic1' , 1, ' ')
         !call add_default( 'thld_macmic2' , 1, ' ')
         !call add_default( 'qtd_macmic1' , 1, ' ')
         !call add_default( 'qtd_macmic2' , 1, ' ')
         !call add_default( 'dthl_macmic1' , 1, ' ')
         !call add_default( 'dthl_macmic2' , 1, ' ')
         !call add_default( 'dqt_macmic1' , 1, ' ')
         !call add_default( 'dqt_macmic2' , 1, ' ')
         !call add_default( 'dthlu_macmic1' , 1, ' ')
         !call add_default( 'dthlu_macmic2' , 1, ' ')
         !call add_default( 'dqtu_macmic1' , 1, ' ')
         !call add_default( 'dqtu_macmic2' , 1, ' ')
         !call add_default( 'dthld_macmic1' , 1, ' ')
         !call add_default( 'dthld_macmic2' , 1, ' ')
         !call add_default( 'dqtd_macmic1' , 1, ' ')
         !call add_default( 'dqtd_macmic2' , 1, ' ')
         !call add_default( 'ztop_macmic1' , 1, ' ')
         !call add_default( 'ztop_macmic2' , 1, ' ')
         !call add_default( 'ddcp_macmic1' , 1, ' ')
         !call add_default( 'ddcp_macmic2' , 1, ' ')
       end if
    end if

    if (history_amwg) then
       call add_default('PBLH',             1, ' ')
    end if

    if (do_clubb_mf_diag) then
       call add_default( 'edmf_DRY_A'    , 1, ' ')
       call add_default( 'edmf_MOIST_A'  , 1, ' ')
       call add_default( 'edmf_DRY_W'    , 1, ' ')
       call add_default( 'edmf_MOIST_W'  , 1, ' ')
       call add_default( 'edmf_DRY_QT'   , 1, ' ')
       call add_default( 'edmf_MOIST_QT' , 1, ' ')
       call add_default( 'edmf_DRY_THL'  , 1, ' ')
       call add_default( 'edmf_MOIST_THL', 1, ' ')
       call add_default( 'edmf_DRY_U'    , 1, ' ')
       call add_default( 'edmf_MOIST_U'  , 1, ' ')
       call add_default( 'edmf_DRY_V'    , 1, ' ')
       call add_default( 'edmf_MOIST_V'  , 1, ' ')
       call add_default( 'edmf_MOIST_QC' , 1, ' ')
       call add_default( 'edmf_S_AE'     , 1, ' ')
       call add_default( 'edmf_S_AW'     , 1, ' ')
       call add_default( 'edmf_S_AWTHL'  , 1, ' ')
       call add_default( 'edmf_S_AWQT'   , 1, ' ')
       call add_default( 'edmf_S_AWU'    , 1, ' ')
       call add_default( 'edmf_S_AWV'    , 1, ' ')
       call add_default( 'edmf_thlflx'   , 1, ' ')
       call add_default( 'edmf_qtflx'    , 1, ' ')
    end if

    if (history_budget) then       
       call add_default('DPDLFLIQ',         history_budget_histfile_num, ' ')
       call add_default('DPDLFICE',         history_budget_histfile_num, ' ')
       call add_default('DPDLFT',           history_budget_histfile_num, ' ') 
       call add_default('STEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('RCMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RIMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RVMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('UTEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('VTEND_CLUBB',      history_budget_histfile_num, ' ')
    endif
     

    ! --------------- !
    ! First step?     !
    ! Initialization  !
    ! --------------- !

    !  Is this the first time step?  If so then initialize CLUBB variables as follows
    if (is_first_step()) then

       call pbuf_set_field(pbuf2d, wp2_idx,     w_tol_sqd)
       call pbuf_set_field(pbuf2d, wp3_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wpthlp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wprtp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, rtp2_idx,    rt_tol**2)
       call pbuf_set_field(pbuf2d, thlp2_idx,   thl_tol**2)
       call pbuf_set_field(pbuf2d, up2_idx,     w_tol_sqd)
       call pbuf_set_field(pbuf2d, vp2_idx,     w_tol_sqd)
      
       call pbuf_set_field(pbuf2d, rtp3_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, thlp3_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, up3_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, vp3_idx,     0.0_r8)
      
       call pbuf_set_field(pbuf2d, upwp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, vpwp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, wpthvp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wp2thvp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthvp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, thlpthvp_idx,0.0_r8)
       call pbuf_set_field(pbuf2d, rcm_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, cloud_frac_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, tke_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, kvh_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, radf_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, wp2rtp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wp2thlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, uprcp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, vprcp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rc_coef_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, wp4_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wpup2_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wpvp2_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wp2up2_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wp2vp2_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, ice_supersat_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, thlm_macmic1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, thlm_macmic2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, qtm_macmic1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, qtm_macmic2_idx, 0.0_r8)

       ! Initialize SILHS covariance contributions
       call pbuf_set_field(pbuf2d, rtp2_mc_zt_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, thlp2_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wprtp_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wpthlp_mc_zt_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthlp_mc_zt_idx, 0.0_r8)

       if (do_clubb_mf) then
         call pbuf_set_field(pbuf2d, ztopmn_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ztopma_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ztopm1_macmic_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ddcp_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ddcp_macmic_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ddcpmn_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, cbm1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, cbm1_macmic_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, up_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, up_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dn_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dn_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, upa_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, upa_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dna_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dna_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, thlu_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, thlu_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, qtu_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, qtu_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, thld_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, thld_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, qtd_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, qtd_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthl_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthl_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqt_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqt_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthlu_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthlu_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqtu_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqtu_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthld_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dthld_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqtd_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, dqtd_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ztop_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ztop_macmic2_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ddcp_macmic1_idx, 0.0_r8)
         call pbuf_set_field(pbuf2d, ddcp_macmic2_idx, 0.0_r8)
       end if
       call pbuf_set_field(pbuf2d, pdf_zm_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_mixt_frac_idx, 0.0_r8)

    endif
  
    ! The following is physpkg, so it needs to be initialized every time
    call pbuf_set_field(pbuf2d, fice_idx,    0.0_r8)
   
    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

#endif
    end subroutine clubb_ini_cam
    
    
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  subroutine clubb_tend_cam( state,   ptend_all,   pbuf,     hdtime, &
                             cmfmc,   cam_in,                        & 
                             macmic_it, cld_macmic_num_steps,dlf, det_s, det_ice)

  !-------------------------------------------------------------------------------
  ! Description: Provide tendencies of shallow convection, turbulence, and 
  !              macrophysics from CLUBB to CAM
  !   
  ! Author: Cheryl Craig, March 2011
  ! Modifications: Pete Bogenschutz, March 2011 and onward
  ! Origin: Based heavily on UWM clubb_init.F90
  ! References:
  !   None
  !-------------------------------------------------------------------------------

    use physics_types,  only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init, &
                              physics_ptend_sum, physics_update, set_wet_to_dry

    use physics_buffer, only: pbuf_old_tim_idx, pbuf_get_field, physics_buffer_desc
    use physics_buffer, only: pbuf_set_field

    use constituents,   only: cnst_get_ind, cnst_type
    use camsrfexch,     only: cam_in_t
    use time_manager,   only: is_first_step   
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
    use tropopause,     only: tropopause_findChemTrop
    use time_manager,   only: get_nstep, is_first_restart_step

    use wv_saturation,   only: qsat
    use interpolate_data,only: vertinterp

#ifdef CLUBB_SGS
    use hb_diff,                   only: pblintd
    use scamMOD,                   only: single_column,scm_clubb_iop_name
    use clubb_api_module, only: &
      nparams, &
      setup_parameters_api, &
      time_precision, &
      advance_clubb_core_api, &
      zt2zm_api, zm2zt_api, &
      setup_grid_heights_api, &
      em_min, &
      w_tol_sqd, &
      rt_tol, &
      thl_tol, &
      l_stats, &
      stats_tsamp, &
      stats_tout, &
      l_output_rad_files, &
      stats_begin_timestep_api, &
      hydromet_dim, calculate_thlp2_rad_api, update_xp2_mc_api, &
      sat_mixrat_liq_api, &
      fstderr, &
      ipdf_post_advance_fields, &
      copy_single_pdf_params_to_multi, &
      copy_multi_pdf_params_to_single, &
      pdf_parameter, &
      init_pdf_params_api, &
      init_pdf_implicit_coefs_terms_api, &
      setup_grid_api

    use clubb_api_module, only: &
      clubb_fatal_error    ! Error code value to indicate a fatal error

    use cldfrc2m,                  only: aist_vector, rhmini_const, rhmaxi_const, rhminis_const, rhmaxis_const
    use cam_history,               only: outfld

    use macrop_driver,             only: liquid_macro_tend
    use clubb_mf,                  only: integrate_mf
   
    use perf_mod

#endif

    implicit none
   
    ! ---------------------------------------------------- !
    !                    Input Auguments                   !
    ! ---------------------------------------------------- !

    type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
    type(cam_in_t),      intent(in)    :: cam_in
    real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]
    real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
    real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
    integer,             intent(in)    :: cld_macmic_num_steps     ! number of mac-mic iterations
    integer,             intent(in)    :: macmic_it                ! number of mac-mic iterations
    
    ! ---------------------------------------------------- !
    !                Input-Output Auguments                !
    ! ---------------------------------------------------- !
    
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! ---------------------------------------------------- !
    !                   Output Auguments                   !
    ! ---------------------------------------------------- !

    type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

    ! These two variables are needed for energy check    
    real(r8),            intent(out)   :: det_s(pcols)              ! Integral of detrained static energy from ice
    real(r8),            intent(out)   :: det_ice(pcols)            ! Integral of detrained ice for energy check

        
    ! ---------------------------------------------------- !
    !                   Local Variables                    !
    ! ---------------------------------------------------- !

#ifdef CLUBB_SGS

    type(physics_state) :: state1                ! Local copy of state variable
    type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all
   
    integer :: i, j,  k, t, ixind, nadv
    integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
    integer :: itim_old
    integer :: ncol, lchnk                       ! # of columns, and chunk identifier
    integer :: err_code                          ! Diagnostic, for if some calculation goes amiss.
    integer :: icnt  
    logical :: lq2(pcnst)

    integer :: iter

    integer :: clubbtop(pcols)

    real(r8) :: frac_limit, ic_limit

    real(r8) :: dtime                                   ! CLUBB time step                               [s]
    real(r8) :: zt_out(pcols,pverp)                        ! output for the thermo CLUBB grid                   [m]
    real(r8) :: zi_out(pcols,pverp)                        ! output for momentum CLUBB grid                     [m]
    real(r8) :: ubar                                      ! surface wind                                [m/s]
    real(r8) :: ustar                                     ! surface stress                              [m/s]
    real(r8) :: z0                                        ! roughness height                            [m]
    real(r8) :: bflx22(pcols)                          ! Variable for buoyancy flux for pbl            [K m/s]
    real(r8) :: qclvar(pcols,pverp)              ! cloud water variance                          [kg^2/kg^2]
    real(r8) :: zo(pcols)                               ! roughness height                              [m]
    real(r8) :: dz_g(pcols,pver)                       ! thickness of layer                            [m]
    real(r8) :: relvarmax

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol) :: &
      fcor, &                             ! Coriolis forcing                            [s^-1]
      sfc_elevation, &                    ! Elevation of ground                         [m AMSL][m]
      wpthlp_sfc, &                       ! w' theta_l' at surface                      [(m K)/s]
      wprtp_sfc, &                        ! w' r_t' at surface                          [(kg m)/( kg s)]
      upwp_sfc, &                         ! u'w' at surface                             [m^2/s^2]
      vpwp_sfc, &                         ! v'w' at surface                             [m^2/s^2]
      upwp_sfc_pert, &                    ! perturbed u'w' at surface                   [m^2/s^2]
      vpwp_sfc_pert, &                    ! perturbed v'w' at surface                   [m^2/s^2]
      grid_dx, grid_dy                    ! CAM grid [m]

    real(r8), dimension(state%ncol,sclr_dim) :: &
      wpsclrp_sfc            ! Scalar flux at surface                        [{units vary} m/s]

    real(r8), dimension(state%ncol,edsclr_dim) :: &
      wpedsclrp_sfc        ! Eddy-scalar flux at surface                   [{units vary} m/s]

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,pverp+1-top_lev) :: &
      thlm_forcing,             & ! theta_l forcing (thermodynamic levels)      [K/s]
      rtm_forcing,              & ! r_t forcing (thermodynamic levels)          [(kg/kg)/s]
      um_forcing,               & ! u wind forcing (thermodynamic levels)       [m/s/s]
      vm_forcing,               & ! v wind forcing (thermodynamic levels)       [m/s/s]
      wprtp_forcing,            &
      wpthlp_forcing,           &
      rtp2_forcing,             &
      thlp2_forcing,            &
      rtpthlp_forcing,          &
      wm_zm,                    & ! w mean wind component on momentum levels    [m/s]
      wm_zt,                    & ! w mean wind component on thermo. levels     [m/s]
      rtm_ref,                  & ! Initial profile of rtm                      [kg/kg]
      thlm_ref,                 & ! Initial profile of thlm                     [K]
      um_ref,                   & ! Initial profile of um                       [m/s]
      vm_ref,                   & ! Initial profile of vm                       [m/s]
      ug,                       & ! U geostrophic wind                          [m/s]
      vg,                       & ! V geostrophic wind                          [m/s]
      p_in_Pa,                  & ! Air pressure (thermodynamic levels)         [Pa]
      rho_zm,                   & ! Air density on momentum levels              [kg/m^3]
      rho_zt,                   & ! Air density on thermo levels                [kg/m^3]
      exner,                    & ! Exner function (thermodynamic levels)       [-]
      rho_ds_zm,                & ! Dry, static density on momentum levels              [kg/m^3]
      rho_ds_zt,                & ! Dry, static density on thermodynamic levels         [kg/m^3]
      invrs_rho_ds_zm,          & ! Inv. dry, static density on momentum levels         [m^3/kg]
      invrs_rho_ds_zt,          & ! Inv. dry, static density on thermo. levels          [m^3/kg]
      thv_ds_zm,                & ! Dry, base-state theta_v on momentum levels          [K]
      thv_ds_zt,                & ! Dry, base-state theta_v on thermo. levels           [K]
      rfrzm,                    &
      radf,                     &
      um_in,                    & ! meridional wind                             [m/s]
      vm_in,                    & ! zonal wind                                  [m/s]
      upwp_in,                  & ! meridional wind flux                                [m^2/s^2]
      vpwp_in,                  & ! zonal wind flux                             [m^2/s^2]
      up2_in,                   & ! meridional wind variance                    [m^2/s^2]
      vp2_in,                   & ! zonal wind variance                         [m^2/s^2]
      up3_in,                   & ! meridional wind third-order                   [m^3/s^3]
      vp3_in,                   & ! zonal wind third-order                        [m^3/s^3]
      thlm_in,                  & ! liquid water potential temperature (thetal) [K]
      rvm_in,                   & ! water vapor mixing ratio                      [kg/kg]
      rtm_in,                   & ! total water mixing ratio                    [kg/kg]
      wprtp_in,                 & ! turbulent flux of total water                       [kg/kg m/s]
      wpthlp_in,                & ! turbulent flux of thetal                    [K m/s]
      tke_in,                   & ! TKE                                         [m^2/s^2]
      wp2_in,                   & ! vertical velocity variance (CLUBB)          [m^2/s^2]
      wp3_in,                   & ! third moment vertical velocity              [m^3/s^3]
      rtp2_in,                  & ! total water variance                                [kg^2/kg^2]
      rtp2_zt,                  & ! CLUBB R-tot variance on thermo levs
      thl2_zt,                  & ! CLUBB Theta-l variance on thermo levs         [K^2]
      wp2_zt,                   & ! CLUBB W variance on theromo levs              [m^2/s^2]
      rtp3_in,                  & ! total water 3rd order                               [kg^3/kg^3]
      thlp2_in,                 & ! thetal variance                             [K^2]
      thlp3_in,                 & ! thetal 3rd order                            [K^3]
      rtpthlp_in,               & ! covariance of thetal and qt                 [kg/kg K]
      rcm_inout,                & ! CLUBB output of liquid water mixing ratio   [kg/kg]
      rcm_out_zm,               &
      cloud_frac_inout,         & ! CLUBB output of cloud fraction              [fraction]
      wpthvp_in,                & ! w'th_v' (momentum levels)                   [m/s K]
      wp2thvp_in,               & ! w'^2 th_v' (thermodynamic levels)           [m^2/s^2 K]
      rtpthvp_in,               & ! r_t'th_v' (momentum levels)                 [kg/kg K]
      thlpthvp_in,              & ! th_l'th_v' (momentum levels)                        [K^2]
      ice_supersat_frac_inout,  &
      um_pert_inout,            & ! Perturbed U wind                          [m/s]
      vm_pert_inout,            & ! Perturbed V wind                          [m/s]
      upwp_pert_inout,          & ! Perturbed u'w'                            [m^2/s^2]
      vpwp_pert_inout,          & ! Perturbed v'w'                            [m^2/s^2]
      khzm_out,                 & ! Eddy diffusivity of heat/moisture on momentum (i.e. interface) levels  [m^2/s]
      khzt_out,                 & ! eddy diffusivity on thermo grids              [m^2/s]
      qclvar_out,               & ! cloud water variance                          [kg^2/kg^2]
      thlprcp_out,              &
      wprcp_out,                & ! CLUBB output of flux of liquid water                [kg/kg m/s]
      w_up_in_cloud_out,        &
      w_down_in_cloud_out,      &
      cloudy_updraft_frac_out,  &
      cloudy_downdraft_frac_out,&
      rcm_in_layer_out,         & ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
      cloud_cover_out,          & ! CLUBB output of in-cloud cloud fraction     [fraction]
      invrs_tau_zm_out,         & ! CLUBB output of 1 divided by time-scale     [1/s]
      rtp2_mc_out,              & ! total water tendency from rain evap
      thlp2_mc_out,             & ! thetal tendency from rain evap
      wprtp_mc_out,             &
      wpthlp_mc_out,            &
      rtpthlp_mc_out,           &
      pre_in,                   & ! input for precip evaporation
      qrl_clubb,                &
      qrl_zm,                   &
      wp2rtp_inout,             & ! w'^2 rt' (thermodynamic levels)
      wp2thlp_inout,            & ! w'^2 thl' (thermodynamic levels)
      uprcp_inout,              & ! < u' r_c' > (momentum levels)
      vprcp_inout,              & ! < v' r_c' > (momentum levels)
      rc_coef_inout,            & ! Coef. of X'r_c' in Eq. (34) (t-levs.)
      wp4_inout,                & ! w'^4 (momentum levels
      wpup2_inout,              & ! w'u'^2 (thermodynamic levels)
      wpvp2_inout,              & ! w'v'^2 (thermodynamic levels)
      wp2up2_inout,             & ! w'^2 u'^2 (momentum levels)
      wp2vp2_inout,             & ! w'^2 v'^2 (momentum levels)
      zt_g,                     & ! Thermodynamic grid of CLUBB                 [m]
      zi_g                                       ! Momentum grid of CLUBB                       [m]

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,pverp+1-top_lev,sclr_dim) :: &
      sclrm_forcing,  & ! Passive scalar forcing              [{units vary}/s]
      sclrm,          & ! Passive scalar mean (thermo. levels)          [units vary]
      sclrp2,         & ! sclr'^2 (momentum levels)                     [{units vary}^2]
      sclrp3,         & ! sclr'^3 (thermo. levels)                      [{units vary}^3]
      sclrprtp,       & ! sclr'rt' (momentum levels)          [{units vary} (kg/kg)]
      sclrpthlp,      & ! sclr'thlp' (momentum levels)        [{units vary} (K)]
      wpsclrp           ! w'sclr' (momentum levels)                     [{units vary} m/s]

    real(r8), dimension(state%ncol,pverp,sclr_dim) :: &
      sclrpthvp_inout  ! sclr'th_v' (momentum levels)                  [{units vary} (K)]

    real(r8), dimension(state%ncol,pverp+1-top_lev,edsclr_dim) :: &
      edsclrm_forcing,  & ! Eddy passive scalar forcing         [{units vary}/s]
      edsclr_in           ! Scalars to be diffused through CLUBB                [units vary]

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,pverp+1-top_lev,hydromet_dim) :: &
      hydromet,     &
      wphydrometp,  &
      wp2hmp,       &
      rtphmp_zt,    &
      thlphmp_zt
 
    ! Variables below are needed to compute energy integrals for conservation
    ! NOTE: Arrays of size PCOLS (all possible columns) can be used to access State, PBuf and History Subroutines
    real(r8) :: ke_a(pcols), ke_b(pcols), te_a(pcols), te_b(pcols)
    real(r8) :: wv_a(pcols), wv_b(pcols), wl_b(pcols), wl_a(pcols)
    real(r8) :: se_dis(pcols), se_a(pcols), se_b(pcols), clubb_s(pcols,pver)
    real(r8) :: eleak(pcols)

    real(r8) :: inv_exner_clubb(pcols,pverp)     ! Inverse exner function consistent with CLUBB  [-]
    real(r8) :: inv_exner_clubb_surf(pcols)      ! Inverse exner function at the surface
    real(r8) :: wpthlp_output(pcols,pverp)       ! Heat flux output variable                     [W/m2]
    real(r8) :: wprtp_output(pcols,pverp)        ! Total water flux output variable              [W/m2]
    real(r8) :: wp3_output(pcols,pverp)          ! wp3 output                                    [m^3/s^3]
    real(r8) :: rtpthlp_output(pcols,pverp)      ! rtpthlp ouptut                                [K kg/kg]
    real(r8) :: qt_output(pcols,pver)            ! Total water mixing ratio for output           [kg/kg]
    real(r8) :: thetal_output(pcols,pver)        ! Liquid water potential temperature output     [K]
    real(r8) :: sl_output(pcols,pver)            ! Liquid water static energy                    [J/kg]
    real(r8) :: ustar2(pcols)                    ! Surface stress for PBL height                 [m2/s2]
    real(r8) :: rho(pcols,pverp)     		! Midpoint density in CAM      			[kg/m^3]
    real(r8) :: thv(pcols,pverp)   		! virtual potential temperature			[K]
    real(r8) :: edsclr_out(pcols,pverp,edsclr_dim)     ! Scalars to be diffused through CLUBB		[units vary]
    real(r8) :: rcm_in_layer(pcols,pverp)	! CLUBB in-cloud liquid water mixing ratio	[kg/kg]
    real(r8) :: cloud_cover(pcols,pverp)		! CLUBB in-cloud cloud fraction			[fraction]
    real(r8) :: wprcp(pcols,pverp)		! CLUBB liquid water flux			[m/s kg/kg]
    real(r8) :: wpthvp_diag(pcols,pverp)		! CLUBB buoyancy flux				[W/m^2]
    real(r8) :: rvm(pcols,pverp)
    real(r8) :: pdfp_rtp2(pcols, pverp)          ! Calculated R-tot variance from pdf_params     [kg^2/kg^2]
    real(r8) :: rtp2_zt_out(pcols, pverp)        ! CLUBB R-tot variance on thermo levs           [kg^2/kg^2]
    real(r8) :: thl2_zt_out(pcols, pverp)        ! CLUBB Theta-l variance on thermo levs
    real(r8) :: wp2_zt_out(pcols, pverp) 
    real(r8) :: dlf_liq_out(pcols, pverp)        ! Detrained liquid water from ZM                [kg/kg/s]
    real(r8) :: dlf_ice_out(pcols, pverp)        ! Detrained ice water from ZM                   [kg/kg/s]
    real(r8) :: wm_zt_out(pcols, pverp)          ! CLUBB mean W on thermo levs output            [m/s]
    real(r8) :: mean_rt                          ! Calculated R-tot mean from pdf_params (temp)  [kg/kg]
    real(r8) :: dlf2(pcols,pver)                 ! Detraining cld H20 from shallow convection    [kg/kg/day]
    real(r8) :: eps                              ! Rv/Rd                                         [-]
    real(r8) :: dum1                             ! dummy variable                                [units vary]
    real(r8) :: obklen(pcols)                    ! Obukov length                                 [m]
    real(r8) :: kbfs(pcols)                      ! Kinematic Surface heat flux                   [K m/s]
    real(r8) :: th(pcols,pver)                   ! potential temperature                         [K]
    real(r8) :: dummy2(pcols)                    ! dummy variable                                [units vary]
    real(r8) :: dummy3(pcols)                    ! dummy variable                                [units vary]
    real(r8) :: kinheat(pcols)                   ! Kinematic Surface heat flux                   [K m/s]
    real(r8) :: rrho(pcols)                      ! Inverse of air density                        [1/kg/m^3]
    real(r8) :: kinwat(pcols)                    ! Kinematic water vapor flux                    [m/s]
    real(r8) :: latsub
    real(r8) :: thlp2_rad_out(pcols,pverp+1-top_lev)
    real(r8) :: apply_const, rtm_test
    real(r8) :: dl_rad, di_rad, dt_low

    character(len=200) :: temp1, sub             ! Strings needed for CLUBB output
    real(kind=time_precision)                 :: time_elapsed                ! time keep track of stats          [s]
    integer :: stats_nsamp, stats_nout           ! Stats sampling and output intervals for CLUBB [timestep]

    real(r8) :: rtm_integral_1, rtm_integral_update, rtm_integral_forcing

    ! ---------------------------------------------------- !
    !                    Pointers                          !
    ! ---------------------------------------------------- !
    
    real(r8), pointer, dimension(:,:) :: wp2      ! vertical velocity variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wp3      ! third moment of vertical velocity		[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: wpthlp   ! turbulent flux of thetal			[m/s K]
    real(r8), pointer, dimension(:,:) :: wprtp    ! turbulent flux of moisture			[m/s kg/kg]
    real(r8), pointer, dimension(:,:) :: rtpthlp  ! covariance of thetal and qt			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: rtp2     ! moisture variance				[kg^2/kg^2]
    real(r8), pointer, dimension(:,:) :: thlp2    ! temperature variance				[K^2]
    real(r8), pointer, dimension(:,:) :: rtp3     ! moisture 3rd order				[kg^3/kg^3]
    real(r8), pointer, dimension(:,:) :: thlp3    ! temperature 3rd order			[K^3]
    real(r8), pointer, dimension(:,:) :: up2      ! east-west wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vp2      ! north-south wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: up3      ! east-west wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: vp3      ! north-south wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: upwp     ! east-west momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vpwp     ! north-south momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wpthvp   ! w'th_v' (momentum levels)			[m/s K]
    real(r8), pointer, dimension(:,:) :: wp2thvp  ! w'^2 th_v' (thermodynamic levels)		[m^2/s^2 K]
    real(r8), pointer, dimension(:,:) :: rtpthvp  ! r_t'th_v' (momentum levels)			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: thlpthvp ! th_l'th_v' (momentum levels)			[K^2]
    real(r8), pointer, dimension(:,:) :: cloud_frac ! Cloud fraction (thermodynamic levels)	[K^2]
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_1        !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_2        !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_1 !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_2 !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_mixt_frac  !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: wp2rtp    ! w'^2 rt' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2thlp   ! w'^2 thl' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: uprcp     ! < u' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: vprcp     ! < v' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: rc_coef   ! Coef. of X'r_c' in Eq. (34) (t-levs.)
    real(r8), pointer, dimension(:,:) :: wp4       ! w'^4 (momentum levels
    real(r8), pointer, dimension(:,:) :: wpup2     ! w'u'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wpvp2     ! w'v'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2up2    ! w'^2 u'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: wp2vp2    ! w'^2 v'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: thlm     ! mean temperature				[K]
    real(r8), pointer, dimension(:,:) :: rtm      ! mean moisture mixing ratio			[kg/kg]
    real(r8), pointer, dimension(:,:) :: rcm      ! CLUBB cloud water mixing ratio               [kg/kg]
    real(r8), pointer, dimension(:)   :: ztodtptr ! timestep to send to SILHS
    real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind				[m/s]
    real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind			[m/s]
    real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction 				[fraction]
    real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction 		[fraction]    
    real(r8), pointer, dimension(:,:) :: khzm     ! CLUBB's eddy diffusivity of heat/moisture on momentum (i.e. interface) levels          [m^2/s]
    real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height                [m]
    real(r8), pointer, dimension(:,:) :: tke      ! turbulent kinetic energy                     [m^2/s^2]
    real(r8), pointer, dimension(:,:) :: dp_icwmr ! deep convection in cloud mixing ratio        [kg/kg]
    real(r8), pointer, dimension(:,:) :: sh_icwmr ! shallow convection (EDMF) in cloud mixing ratio [kg/kg]
    real(r8), pointer, dimension(:,:) :: ice_supersat_frac ! Cloud fraction of ice clouds (pverp)[fraction] 
    real(r8), pointer, dimension(:,:) :: relvar   ! relative cloud water variance                [-]
    real(r8), pointer, dimension(:,:) :: accre_enhan ! accretion enhancement factor              [-]
    real(r8), pointer, dimension(:,:) :: naai
    real(r8), pointer, dimension(:,:) :: cmeliq 
    real(r8), pointer, dimension(:,:) :: cmfmc_sh ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]

    real(r8), pointer, dimension(:,:) :: qsatfac
    real(r8), pointer, dimension(:,:) :: npccn
    real(r8), pointer, dimension(:,:) :: prer_evap
    real(r8), pointer, dimension(:,:) :: qrl
    real(r8), pointer, dimension(:,:) :: radf_clubb

    ! SILHS covariance contributions
    real(r8), pointer, dimension(:,:) :: rtp2_mc_zt
    real(r8), pointer, dimension(:,:) :: thlp2_mc_zt
    real(r8), pointer, dimension(:,:) :: wprtp_mc_zt
    real(r8), pointer, dimension(:,:) :: wpthlp_mc_zt
    real(r8), pointer, dimension(:,:) :: rtpthlp_mc_zt

    real(r8)  qitend(pcols,pver)
    real(r8)  initend(pcols,pver)  ! Needed for ice supersaturation adjustment calculation

    ! ZM microphysics
    real(r8), pointer :: dlfzm(:,:)  ! ZM detrained convective cloud water mixing ratio.
    real(r8), pointer :: difzm(:,:)  ! ZM detrained convective cloud ice mixing ratio.
    real(r8), pointer :: dnlfzm(:,:) ! ZM detrained convective cloud water num concen.
    real(r8), pointer :: dnifzm(:,:) ! ZM detrained convective cloud ice num concen.

    real(r8)                          :: stend(pcols,pver)
    real(r8)                          :: qvtend(pcols,pver)
    real(r8)                          :: qctend(pcols,pver)
    real(r8)                          :: inctend(pcols,pver)
    real(r8)                          :: fqtend(pcols,pver)
    real(r8)                          :: rhmini(pcols)
    real(r8)                          :: rhmaxi(pcols)
    integer                           :: troplev(pcols)
    logical                           :: lqice(pcnst)
    logical                           :: apply_to_surface(pcols)

    ! CLUBB-MF pointers
    real(r8),pointer :: prec_sh(:)   ! total precipitation from MF
    real(r8),pointer :: snow_sh(:)   ! snow from MF
 
    real(r8), pointer :: ztopmn(:,:)
    real(r8), pointer :: ztopma(:)
    real(r8), pointer :: ztopm1_macmic(:)
 
    real(r8), pointer :: ddcp(:)
    real(r8), pointer :: ddcp_macmic(:)
    real(r8), pointer :: ddcpmn(:,:)
 
    real(r8), pointer :: cbm1(:)
    real(r8), pointer :: cbm1_macmic(:)
 
    real(r8), pointer :: qtm_macmic1(:,:)
    real(r8), pointer :: qtm_macmic2(:,:)
    real(r8), pointer :: thlm_macmic1(:,:)
    real(r8), pointer :: thlm_macmic2(:,:)
    real(r8), pointer :: rcm_macmic(:,:)
    real(r8), pointer :: cldfrac_macmic(:,:)
    real(r8), pointer :: wpthlp_macmic(:,:)
    real(r8), pointer :: wprtp_macmic(:,:)
    real(r8), pointer :: wpthvp_macmic(:,:)
    real(r8), pointer :: mf_thlflx_macmic(:,:)
    real(r8), pointer :: mf_qtflx_macmic(:,:)
    real(r8), pointer :: mf_thvflx_macmic(:,:)
    real(r8), pointer :: up_macmic1(:,:)
    real(r8), pointer :: up_macmic2(:,:)
    real(r8), pointer :: dn_macmic1(:,:)
    real(r8), pointer :: dn_macmic2(:,:)
    real(r8), pointer :: upa_macmic1(:,:)
    real(r8), pointer :: upa_macmic2(:,:)
    real(r8), pointer :: dna_macmic1(:,:)
    real(r8), pointer :: dna_macmic2(:,:)
    real(r8), pointer :: thlu_macmic1(:,:)
    real(r8), pointer :: thlu_macmic2(:,:)
    real(r8), pointer :: qtu_macmic1(:,:)
    real(r8), pointer :: qtu_macmic2(:,:)
    real(r8), pointer :: thld_macmic1(:,:)
    real(r8), pointer :: thld_macmic2(:,:)
    real(r8), pointer :: qtd_macmic1(:,:)
    real(r8), pointer :: qtd_macmic2(:,:)
    real(r8), pointer :: dthl_macmic1(:,:)
    real(r8), pointer :: dthl_macmic2(:,:)
    real(r8), pointer :: dqt_macmic1(:,:)
    real(r8), pointer :: dqt_macmic2(:,:)
    real(r8), pointer :: dthlu_macmic1(:,:)
    real(r8), pointer :: dthlu_macmic2(:,:)
    real(r8), pointer :: dqtu_macmic1(:,:)
    real(r8), pointer :: dqtu_macmic2(:,:)
    real(r8), pointer :: dthld_macmic1(:,:)
    real(r8), pointer :: dthld_macmic2(:,:)
    real(r8), pointer :: dqtd_macmic1(:,:)
    real(r8), pointer :: dqtd_macmic2(:,:)
    real(r8), pointer :: ztop_macmic1(:,:)
    real(r8), pointer :: ztop_macmic2(:,:)
    real(r8), pointer :: ddcp_macmic1(:,:)
    real(r8), pointer :: ddcp_macmic2(:,:)

    real(r8), dimension(pcols)           :: mf_ztop_output,    mf_L0_output,        &
                                            mf_cape_output,    mf_cfl_output,       &
                                            mf_ddcp_output,    mf_freq_output
    !
    ! MF outputs to outfld
    real(r8), dimension(pcols,pver)      :: mf_thlforcup_output, mf_qtforcup_output,  & ! thermodynamic grid
                                            mf_thlforcdn_output, mf_qtforcdn_output,  & ! thermodynamic grid
                                            mf_thlforc_output,   mf_qtforc_output,    & ! thermodynamic grid
                                            mf_ent_output,                            & ! thermodynamic grid
                                            mf_sqtup_output,     mf_sqtdn_output,     & ! thermodynamic grid
                                            mf_qc_output,        mf_cloudfrac_output    ! thermodynamic grid
 
    ! MF plume level outputs
    real(r8), dimension(pcols,pverp,clubb_mf_nup) ::           mf_upa_flip,         &
                                                               mf_upw_flip,         &
                                                               mf_upmf_flip,        &
                                                               mf_upqt_flip,        &
                                                               mf_upthl_flip,       &
                                                               mf_upthv_flip,       &
                                                               mf_upth_flip,        &
                                                               mf_upqc_flip,        &
                                                               mf_upbuoy_flip,      &
                                                               mf_upent_flip,       &
                                                               mf_updet_flip
    ! MF plume level outputs to outfld
    real(r8), dimension(pcols,pverp*clubb_mf_nup) ::           mf_upa_output,       &
                                                               mf_upw_output,       &
                                                               mf_upmf_output,      &
                                                               mf_upqt_output,      &
                                                               mf_upthl_output,     &
                                                               mf_upthv_output,     &
                                                               mf_upth_output,      &
                                                               mf_upqc_output,      &
                                                               mf_upent_output,     &
                                                               mf_updet_output,     &
                                                               mf_upbuoy_output
    ! MF plume level outputs
    real(r8), dimension(pcols,pverp,clubb_mf_nup) ::           mf_dnw_flip,         &
                                                               mf_dnthl_flip,       &
                                                               mf_dnqt_flip
 
    ! MF plume level outputs to outfld
    real(r8), dimension(pcols,pverp*clubb_mf_nup) ::           mf_dnw_output,       &
                                                               mf_dnthl_output,     &
                                                               mf_dnqt_output
 
    ! MF Plume
    real(r8), pointer                    :: tpert(:)
 
    ! MF outputs to outfld
    ! NOTE: Arrays of size PCOLS (all possible columns) can be used to access State, PBuf and History Subroutines
    real(r8), dimension(pcols,pverp)     :: mf_dry_a_output,   mf_moist_a_output,   &
                                            mf_dry_w_output,   mf_moist_w_output,   &
                                            mf_dry_qt_output,  mf_moist_qt_output,  &
                                            mf_dry_thl_output, mf_moist_thl_output, &
                                            mf_dry_u_output,   mf_moist_u_output,   &
                                            mf_dry_v_output,   mf_moist_v_output,   &
                                                               mf_moist_qc_output,  &
                                            s_ae_output,       s_aw_output,         &
                                            s_awthlup_output,  s_awqtup_output, s_awuup_output, s_awvup_output,        &
                                            s_awthldn_output,  s_awqtdn_output, s_awudn_output, s_awvdn_output,        &
                                            s_awthl_output,    s_awqt_output,       &
                                            s_awu_output,      s_awv_output,        &
                                            s_aww_output,                                                              &
                                            mf_thlflxup_output,mf_qtflxup_output, mf_uflxup_output, mf_vflxup_output,  &
                                            mf_thlflxdn_output,mf_qtflxdn_output, mf_uflxdn_output, mf_vflxdn_output,  &
                                            mf_thlflx_output,  mf_qtflx_output,   mf_uflx_output,   mf_vflx_output,    &
                                            mf_thvflx_output,                                                          &
                                            mf_rcm_output,     mf_precc_output
    ! MF work arrays (of size NCOL)
    real(r8), dimension(state%ncol,pverp):: mf_dry_a,   mf_moist_a,    &
                                            mf_dry_w,   mf_moist_w,    &
                                            mf_dry_qt,  mf_moist_qt,   &
                                            mf_dry_thl, mf_moist_thl,  &
                                            mf_dry_u,   mf_moist_u,    &
                                            mf_dry_v,   mf_moist_v,    &
                                                        mf_moist_qc,   &
                                            s_ae,       s_ac,          &
                                            s_aup,      s_adn,         &
                                            s_aw,                      &
                                            s_awup,     s_awdn,        &
                                            s_aww,                     &
                                            s_awwup,    s_awwdn,       &
                                            s_awthlup,  s_awqtup, s_awuup, s_awvup,      &
                                            s_awthldn,  s_awqtdn, s_awudn, s_awvdn,      &
                                            s_awthl,    s_awqt,        &
                                            s_awu,      s_awv,         &          
                                            mf_sqtup,   mf_sthlup,     &
                                            mf_sqtdn,   mf_sthldn,     &
                                            mf_sqt,     mf_sthl,       &
                                                        mf_precc       

   real(r8), dimension(state%ncol,pverp) :: mf_thlflxup,      mf_qtflxup,       mf_uflxup,  mf_vflxup,  &
                                            mf_thlflxdn,      mf_qtflxdn,       mf_uflxdn,  mf_vflxdn,  &
                                            mf_thlflx,        mf_qtflx,         mf_uflx,    mf_vflx,    &
                                            mf_thvflx,                                                  &
                                            mf_thlforcup,     mf_qtforcup,                              &
                                            mf_thlforcdn,     mf_qtforcdn,                              &
                                            mf_thlforcup_nadv,mf_qtforcup_nadv,                         &
                                            mf_thlforcdn_nadv,mf_qtforcdn_nadv,                         &
                                            mf_thlforc_nadv,  mf_qtforc_nadv,                           &
                                            mf_qc,            mf_cloudfrac,                             &
                                            mf_qc_nadv,       mf_cloudfrac_nadv,                        &
                                            mf_qc_zt,         mf_cloudfrac_zt,                          &
                                            mf_rcm,           mf_rcm_nadv,                              &
                                            mf_ent_nadv

   real(r8), dimension(state%ncol,pverp,clubb_mf_nup) :: mf_upa,    mf_dna,       &
                                                         mf_upw,    mf_dnw,       &
                                                         mf_upmf,                 &
                                                         mf_upqt,   mf_dnqt,      &
                                                         mf_upthl,  mf_dnthl,     &
                                                         mf_upthv,  mf_dnthv,     &
                                                         mf_upth,   mf_dnth,      &
                                                         mf_upqc,   mf_dnqc,      &
                                                         mf_upbuoy,               &
                                                         mf_updet,                &
                                                         mf_upent

    real(r8) :: inv_rh2o ! To reduce the number of divisions in clubb_tend

    real(r8), dimension(state%ncol,pverp,clubb_mf_nup) :: flip
    real(r8), dimension(state%ncol,pverp) :: lilflip

    ! CFL limiter vars
    real(r8), parameter                  :: cflval = 1._r8
    real(r8)                             :: lambda
    real(r8), dimension(state%ncol)      :: cflfac,     max_cfl,        &
                                            th_sfc,     max_cfl_nadv
    logical                              :: cfllim

    real(r8), dimension(state%ncol)      :: mf_ztop,    mf_ztop_nadv,   &
                                            mf_ztopm1,  mf_ztopm1_nadv, &
                                            mf_precc_nadv, mf_snow_nadv,&
                                            mf_L0,      mf_L0_nadv,     &
                                            mf_ddcp,    mf_ddcp_nadv,   &
                                            mf_cbm1,    mf_cbm1_nadv,   &
                                                        mf_freq_nadv

    real(r8), dimension(state%ncol,pver) :: esat,      rh
    real(r8), dimension(state%ncol,pver) :: mq,        mqsat
    real(r8), dimension(state%ncol)      :: rhlev,     rhinv

    real(r8), dimension(state%ncol,pverp):: rtm_zm_in,  thlm_zm_in,    & ! momentum grid
                                            dzt,        invrs_dzt,     & ! thermodynamic grid
                                                        invrs_exner_zt,& ! thermodynamic grid
                                            kappa_zt,   qc_zt,         & ! thermodynamic grid
                                            th_zt,      qv_zt,         & ! momentum grid
                                            th_zm,      qv_zm,         & ! momentum grid
                                                        qc_zm,         & ! momentum grid
                                            kappa_zm,   p_in_Pa_zm,    & ! momentum grid
                                            dzm,        invrs_exner_zm   ! momentum grid

    real(r8) :: temp2d(pcols,pver), temp2dp(pcols,pverp)  ! temporary array for holding scaled outputs

    integer :: nlev

    intrinsic :: max

    character(len=*), parameter :: subr='clubb_tend_cam'
                          
    type(grid) :: gr
    integer :: begin_height, end_height
    
    type(nu_vertical_res_dep) :: nu_vert_res_dep   ! Vertical resolution dependent nu values
    real(r8) :: lmin
 
#endif
    det_s(:)   = 0.0_r8
    det_ice(:) = 0.0_r8

#ifdef CLUBB_SGS

    !-----------------------------------------------------------------------------------!
    !                           MAIN COMPUTATION BEGINS HERE                            !                                              !
    !-----------------------------------------------------------------------------------!

    call t_startf("clubb_tend_cam")

    nlev = pver + 1 - top_lev

    rtp2_zt_out = 0._r8
    thl2_zt_out = 0._r8
    wp2_zt_out  = 0._r8
    pdfp_rtp2   = 0._r8
    wm_zt_out   = 0._r8

    temp2d      = 0._r8 
    temp2dp     = 0._r8

    dl_rad = clubb_detliq_rad
    di_rad = clubb_detice_rad
    dt_low = clubb_detphase_lowtemp

    frac_limit = 0.01_r8
    ic_limit   = 1.e-12_r8
    inv_rh2o = 1._r8/rh2o

    if (clubb_do_adv) then 
      apply_const = 1._r8  ! Initialize to one, only if CLUBB's moments are advected
    else
      apply_const = 0._r8  ! Never want this if CLUBB's moments are not advected
    endif

    !  Get indicees for cloud and ice mass and cloud and ice number
    call cnst_get_ind('Q',ixq)
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)
    call cnst_get_ind('NUMLIQ',ixnumliq)
    call cnst_get_ind('NUMICE',ixnumice)

    if (clubb_do_icesuper) then
      call pbuf_get_field(pbuf, naai_idx, naai)
    end if

    !  Initialize physics tendency arrays, copy the state to state1 array to use in this routine
    call physics_ptend_init(ptend_all, state%psetcols, 'clubb')

    ! Copy the state to state1 array to use in this routine
    call physics_state_copy(state, state1)

    ! constituents are all treated as dry mmr by clubb
    call set_wet_to_dry(state1)

    if (clubb_do_liqsupersat) then
      call pbuf_get_field(pbuf, npccn_idx, npccn)
    endif

    !  Determine number of columns and which chunk computation is to be performed on
    ncol = state%ncol
    lchnk = state%lchnk

    !  Determine time step of physics buffer
    itim_old = pbuf_old_tim_idx()

    !  Establish associations between pointers and physics buffer fields   
    call pbuf_get_field(pbuf, wp2_idx,     wp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, wp3_idx,     wp3,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, wpthlp_idx,  wpthlp,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, wprtp_idx,   wprtp,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, rtpthlp_idx, rtpthlp, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, rtp2_idx,    rtp2,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, thlp2_idx,   thlp2,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, up2_idx,     up2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, vp2_idx,     vp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

    call pbuf_get_field(pbuf, rtp3_idx,    rtp3,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, thlp3_idx,   thlp3,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, up3_idx,     up3,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, vp3_idx,     vp3,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

    call pbuf_get_field(pbuf, upwp_idx,    upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, vpwp_idx,    vpwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, wpthvp_idx,  wpthvp)
    call pbuf_get_field(pbuf, wp2thvp_idx, wp2thvp)
    call pbuf_get_field(pbuf, rtpthvp_idx, rtpthvp)
    call pbuf_get_field(pbuf, thlpthvp_idx,thlpthvp)
    call pbuf_get_field(pbuf, rcm_idx,     rcm)
    call pbuf_get_field(pbuf, cloud_frac_idx, cloud_frac)

    call pbuf_get_field(pbuf, pdf_zm_w_1_idx, pdf_zm_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, pdf_zm_w_2_idx, pdf_zm_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_1_idx, pdf_zm_varnce_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_2_idx, pdf_zm_varnce_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, pdf_zm_mixt_frac_idx, pdf_zm_mixt_frac, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

    call pbuf_get_field(pbuf, wp2rtp_idx, wp2rtp)
    call pbuf_get_field(pbuf, wp2thlp_idx, wp2thlp)
    call pbuf_get_field(pbuf, uprcp_idx, uprcp)
    call pbuf_get_field(pbuf, vprcp_idx, vprcp)
    call pbuf_get_field(pbuf, rc_coef_idx, rc_coef)
    call pbuf_get_field(pbuf, wp4_idx, wp4)
    call pbuf_get_field(pbuf, wpup2_idx, wpup2)
    call pbuf_get_field(pbuf, wpvp2_idx, wpvp2)
    call pbuf_get_field(pbuf, wp2up2_idx, wp2up2)
    call pbuf_get_field(pbuf, wp2vp2_idx, wp2vp2)
    call pbuf_get_field(pbuf, thlm_idx,    thlm,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, rtm_idx,     rtm,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, um_idx,      um,      start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
    call pbuf_get_field(pbuf, vm_idx,      vm,      start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   
    call pbuf_get_field(pbuf, tke_idx,     tke)
    call pbuf_get_field(pbuf, qrl_idx,     qrl)
    call pbuf_get_field(pbuf, radf_idx,    radf_clubb)

    call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

    call pbuf_get_field(pbuf, qsatfac_idx, qsatfac)

    call pbuf_get_field(pbuf, prer_evap_idx, prer_evap)
    call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
    call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
    call pbuf_get_field(pbuf, ice_supersat_idx, ice_supersat_frac)
    call pbuf_get_field(pbuf, ztodt_idx,   ztodtptr)
    call pbuf_get_field(pbuf, relvar_idx,  relvar)
    call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
    call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
    call pbuf_get_field(pbuf, kvh_idx,     khzm)
    call pbuf_get_field(pbuf, pblh_idx,    pblh)
    call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
    call pbuf_get_field(pbuf, icwmrsh_idx, sh_icwmr)
    call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

    ! SILHS covariance contributions
    call pbuf_get_field(pbuf, rtp2_mc_zt_idx,    rtp2_mc_zt)
    call pbuf_get_field(pbuf, thlp2_mc_zt_idx,   thlp2_mc_zt)
    call pbuf_get_field(pbuf, wprtp_mc_zt_idx,   wprtp_mc_zt)
    call pbuf_get_field(pbuf, wpthlp_mc_zt_idx,  wpthlp_mc_zt)
    call pbuf_get_field(pbuf, rtpthlp_mc_zt_idx, rtpthlp_mc_zt)
   
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )

   call pbuf_get_field(pbuf, qtm_macmic1_idx, qtm_macmic1)
   call pbuf_get_field(pbuf, qtm_macmic2_idx, qtm_macmic2)
   call pbuf_get_field(pbuf, thlm_macmic1_idx, thlm_macmic1)
   call pbuf_get_field(pbuf, thlm_macmic2_idx, thlm_macmic2)

   call pbuf_get_field(pbuf, rcm_macmic_idx, rcm_macmic)
   call pbuf_get_field(pbuf, cldfrac_macmic_idx, cldfrac_macmic)
   call pbuf_get_field(pbuf, wpthlp_macmic_idx, wpthlp_macmic)
   call pbuf_get_field(pbuf, wprtp_macmic_idx, wprtp_macmic)
   call pbuf_get_field(pbuf, wpthvp_macmic_idx, wpthvp_macmic)
   if (do_clubb_mf) then
     call pbuf_get_field(pbuf, mf_wpthlp_macmic_idx, mf_thlflx_macmic)
     call pbuf_get_field(pbuf, mf_wprtp_macmic_idx, mf_qtflx_macmic)
     call pbuf_get_field(pbuf, mf_wpthvp_macmic_idx, mf_thvflx_macmic)
     call pbuf_get_field(pbuf, tpert_idx, tpert)

     call pbuf_get_field(pbuf, ztopmn_idx, ztopmn)
     call pbuf_get_field(pbuf, ztopma_idx, ztopma)
     call pbuf_get_field(pbuf, ztopm1_macmic_idx, ztopm1_macmic)

     call pbuf_get_field(pbuf, ddcp_idx, ddcp)
     call pbuf_get_field(pbuf, ddcp_macmic_idx, ddcp_macmic)
     call pbuf_get_field(pbuf, ddcpmn_idx, ddcpmn)

     call pbuf_get_field(pbuf, cbm1_idx, cbm1)
     call pbuf_get_field(pbuf, cbm1_macmic_idx, cbm1_macmic)

     call pbuf_get_field(pbuf, up_macmic1_idx, up_macmic1)
     call pbuf_get_field(pbuf, up_macmic2_idx, up_macmic2)
     call pbuf_get_field(pbuf, dn_macmic1_idx, dn_macmic1)
     call pbuf_get_field(pbuf, dn_macmic2_idx, dn_macmic2)
     call pbuf_get_field(pbuf, upa_macmic1_idx, upa_macmic1)
     call pbuf_get_field(pbuf, upa_macmic2_idx, upa_macmic2)
     call pbuf_get_field(pbuf, dna_macmic1_idx, dna_macmic1)
     call pbuf_get_field(pbuf, dna_macmic2_idx, dna_macmic2)
     call pbuf_get_field(pbuf, thlu_macmic1_idx, thlu_macmic1)
     call pbuf_get_field(pbuf, thlu_macmic2_idx, thlu_macmic2)
     call pbuf_get_field(pbuf, qtu_macmic1_idx, qtu_macmic1)
     call pbuf_get_field(pbuf, qtu_macmic2_idx, qtu_macmic2)
     call pbuf_get_field(pbuf, thld_macmic1_idx, thld_macmic1)
     call pbuf_get_field(pbuf, thld_macmic2_idx, thld_macmic2)
     call pbuf_get_field(pbuf, qtd_macmic1_idx, qtd_macmic1)
     call pbuf_get_field(pbuf, qtd_macmic2_idx, qtd_macmic2)
     call pbuf_get_field(pbuf, dthl_macmic1_idx, dthl_macmic1)
     call pbuf_get_field(pbuf, dthl_macmic2_idx, dthl_macmic2)
     call pbuf_get_field(pbuf, dqt_macmic1_idx, dqt_macmic1)
     call pbuf_get_field(pbuf, dqt_macmic2_idx, dqt_macmic2)
     call pbuf_get_field(pbuf, dthlu_macmic1_idx, dthlu_macmic1)
     call pbuf_get_field(pbuf, dthlu_macmic2_idx, dthlu_macmic2)
     call pbuf_get_field(pbuf, dqtu_macmic1_idx, dqtu_macmic1)
     call pbuf_get_field(pbuf, dqtu_macmic2_idx, dqtu_macmic2)
     call pbuf_get_field(pbuf, dthld_macmic1_idx, dthld_macmic1)
     call pbuf_get_field(pbuf, dthld_macmic2_idx, dthld_macmic2)
     call pbuf_get_field(pbuf, dqtd_macmic1_idx, dqtd_macmic1)
     call pbuf_get_field(pbuf, dqtd_macmic2_idx, dqtd_macmic2)
     call pbuf_get_field(pbuf, ztop_macmic1_idx, ztop_macmic1)
     call pbuf_get_field(pbuf, ztop_macmic2_idx, ztop_macmic2)
     call pbuf_get_field(pbuf, ddcp_macmic1_idx, ddcp_macmic1)
     call pbuf_get_field(pbuf, ddcp_macmic2_idx, ddcp_macmic2)

     ! SVP
     do k = 1, pver
        call qsat(state%t(1:ncol,k), state%pmid(1:ncol,k), esat(1:ncol,k), rh(1:ncol,k), ncol)
     end do

     rhlev(:ncol) = 0._r8
     if (clubb_mf_Lopt==7 .or. clubb_mf_Lopt==6) then
       ! Interpolate RH to 500 hPa
       rh(:ncol,:) = state%q(:ncol,:,1)/rh(:ncol,:)
       call vertinterp(ncol, ncol, pver, state%pmid(:ncol,:), 50000._r8, rh, rhlev, &
           extrapolate='Z', ln_interp=.true., ps=state%ps(:ncol), phis=state%phis(:ncol), tbot=state%t(:ncol,pver))
     else if (clubb_mf_Lopt==8) then
       ! Mass of q, by layer and vertically integrated
       mq(:ncol,:) = state%q(:ncol,:,1) * state%pdel(:ncol,:) * rga
       mqsat(:ncol,:) = rh(:ncol,:) * state%pdel(:ncol,:) * rga
       do k=2,pver
         mq(:ncol,1) = mq(:ncol,1) + mq(:ncol,k)
         mqsat(:ncol,1) = mqsat(:ncol,1) + mqsat(:ncol,k)
       end do
       rhlev(:ncol) = mq(:ncol,1)/mqsat(:ncol,1) 
     end if     
   end if

    ! Allocate pdf_params only if they aren't allocated already.
    if ( .not. allocated(pdf_params_chnk(lchnk)%mixt_frac) ) then
      call init_pdf_params_api( pverp+1-top_lev, ncol, pdf_params_chnk(lchnk) )
      call init_pdf_params_api( pverp+1-top_lev, ncol, pdf_params_zm_chnk(lchnk) )
    end if
    
    if ( .not. allocated(pdf_implicit_coefs_terms_chnk(lchnk)%coef_wp4_implicit) ) then
      call init_pdf_implicit_coefs_terms_api( pverp+1-top_lev, ncol, sclr_dim, &
                                               pdf_implicit_coefs_terms_chnk(lchnk) )
    end if

    ! Initialize the apply_const variable (note special logic is due to eularian backstepping)
    if (clubb_do_adv .and. (is_first_step() .or. all(wpthlp(1:ncol,1:pver)  ==  0._r8))) then
      apply_const = 0._r8  ! On first time through do not remove constant 
                           !  from moments since it has not been added yet 
    endif

    !  Set the ztodt timestep in pbuf for SILHS
    ztodtptr(:) = 1.0_r8*hdtime

    ! Define the grid box size.  CLUBB needs this information to determine what
    !  the maximum length scale should be.  This depends on the column for 
    !  variable mesh grids and lat-lon grids
    if (single_column) then
      ! If single column specify grid box size to be something
      !  similar to a GCM run
      grid_dx(:) = 100000._r8
      grid_dy(:) = 100000._r8
    else
     
      call grid_size(state1, grid_dx, grid_dy)

    end if

    if (clubb_do_icesuper) then

      ! -------------------------------------- !
      ! Ice Saturation Adjustment Computation  !
      ! -------------------------------------- !

      lq2(:)  = .FALSE.
      lq2(1)  = .TRUE.
      lq2(ixcldice) = .TRUE.
      lq2(ixnumice) = .TRUE.
   
      latsub = latvap + latice
    
      call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq2 )
   
      stend(:ncol,:)=0._r8
      qvtend(:ncol,:)=0._r8
      qitend(:ncol,:)=0._r8
      initend(:ncol,:)=0._r8

      call ice_macro_tend(naai(1:ncol,top_lev:pver), state1%t(1:ncol,top_lev:pver),                       &
                          state1%pmid(1:ncol,top_lev:pver), state1%q(1:ncol,top_lev:pver,1),              &
                          state1%q(1:ncol,top_lev:pver,ixcldice), state1%q(1:ncol,top_lev:pver,ixnumice), &
                          latsub, hdtime, stend(1:ncol,top_lev:pver), qvtend(1:ncol,top_lev:pver),        &
                          qitend(1:ncol,top_lev:pver), initend(1:ncol,top_lev:pver), ncol*(pver-top_lev+1))

      ! update local copy of state with the tendencies
      ptend_loc%q(:ncol,top_lev:pver,1)=qvtend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixcldice)=qitend(:ncol,top_lev:pver)  
      ptend_loc%q(:ncol,top_lev:pver,ixnumice)=initend(:ncol,top_lev:pver)
      ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver) 

      ! Add the ice tendency to the output tendency
      call physics_ptend_sum(ptend_loc, ptend_all, ncol)

      ! ptend_loc is reset to zero by this call
      call physics_update(state1, ptend_loc, hdtime)

      !Write output for tendencies:
      temp2d(:ncol,:pver) =  stend(:ncol,:pver)/cpairv(:ncol,:pver,lchnk)
      call outfld( 'TTENDICE',  temp2d, pcols, lchnk )
      call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
      call outfld( 'QITENDICE', qitend, pcols, lchnk )
      call outfld( 'NITENDICE', initend, pcols, lchnk )
   
    endif


    !  Determine CLUBB time step and make it sub-step friendly
    !  For now we want CLUBB time step to be 5 min since that is 
    !  what has been scientifically validated.  However, there are certain
    !  instances when a 5 min time step will not be possible (based on 
    !  host model time step or on macro-micro sub-stepping
    dtime = clubb_timestep 
   
    !  Now check to see if dtime is greater than the host model 
    !    (or sub stepped) time step.  If it is, then simply 
    !    set it equal to the host (or sub step) time step.  
    !    This section is mostly to deal with small host model
    !    time steps (or small sub-steps) 
    if (dtime > hdtime) then
      dtime = hdtime
    endif
   
    !  Now check to see if CLUBB time step divides evenly into
    !    the host model time step.  If not, force it to divide evenly.
    !    We also want it to be 5 minutes or less.  This section is
    !    mainly for host model time steps that are not evenly divisible
    !    by 5 minutes   
    if (mod(hdtime,dtime) .ne. 0) then
      dtime = hdtime/2._r8
      do while (dtime > clubb_timestep) 
        dtime = dtime/2._r8
      end do
    endif   

    !  If resulting host model time step and CLUBB time step do not divide evenly
    !    into each other, have model throw a fit.  
    if (mod(hdtime,dtime) .ne. 0) then
      call endrun(subr//':  CLUBB time step and HOST time step NOT compatible')
    endif
   
    !  determine number of timesteps CLUBB core should be advanced, 
    !  host time step divided by CLUBB time step  
    nadv = max(hdtime/dtime,1._r8)
  
    !  Initialize forcings for transported scalars to zero 
    sclrm_forcing(:,:,:)   = 0._r8
    edsclrm_forcing(:,:,:) = 0._r8
    sclrm(:,:,:)           = 0._r8
   
    !  Compute inverse exner function consistent with CLUBB's definition, which uses a constant
    !  surface pressure.  CAM's exner (in state) does not.  Therefore, for consistent 
    !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables 
    !  (such as thlm), use "inv_exner_clubb" otherwise use the exner in state
    do k=1,pver
      do i=1,ncol
        inv_exner_clubb(i,k) = 1._r8/((state1%pmid(i,k)/p0_clubb)**(rairv(i,k,lchnk)/cpairv(i,k,lchnk)))
      enddo
    enddo

    !  Compute exner at the surface for converting the sensible heat fluxes
    !  to a flux of potential temperature for use as clubb's boundary conditions
    do i=1,ncol
       inv_exner_clubb_surf(i) = 1._r8/((state1%pmid(i,pver)/p0_clubb)**(rairv(i,pver,lchnk)/cpairv(i,pver,lchnk)))
    enddo

    !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state 
    !  from the CAM state 
    do k=1,pver   ! loop over levels
      do i=1,ncol ! loop over columns

        rtm(i,k)  = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)
        rvm(i,k)  = state1%q(i,k,ixq)
        um(i,k)   = state1%u(i,k)
        vm(i,k)   = state1%v(i,k)
        thlm(i,k) = ( state1%t(i,k) &
                      - (latvap/cpairv(i,k,lchnk))*state1%q(i,k,ixcldliq) ) &
                    * inv_exner_clubb(i,k)

        if (clubb_do_adv) then
          if (macmic_it  ==  1) then 

            !  Note that some of the moments below can be positive or negative.  
            !    Remove a constant that was added to prevent dynamics from clipping                                                                   
            !    them to prevent dynamics from making them positive.  
            thlp2(i,k)   = state1%q(i,k,ixthlp2)
            rtp2(i,k)    = state1%q(i,k,ixrtp2)
            rtpthlp(i,k) = state1%q(i,k,ixrtpthlp) - (rtpthlp_const*apply_const)
            wpthlp(i,k)  = state1%q(i,k,ixwpthlp) - (wpthlp_const*apply_const)
            wprtp(i,k)   = state1%q(i,k,ixwprtp) - (wprtp_const*apply_const)
            wp2(i,k)     = state1%q(i,k,ixwp2)
            wp3(i,k)     = state1%q(i,k,ixwp3) - (wp3_const*apply_const)
            up2(i,k)     = state1%q(i,k,ixup2)
            vp2(i,k)     = state1%q(i,k,ixvp2)
          endif
        endif

      enddo
    enddo
   
    if (clubb_do_adv) then
      ! If not last step of macmic loop then set apply_const back to 
      !   zero to prevent output from being corrupted.  
      if (macmic_it  ==  cld_macmic_num_steps) then    
        apply_const = 1._r8 
      else
        apply_const = 0._r8
      endif
    endif   

    rtm(1:ncol,pverp)  = rtm(1:ncol,pver)
    um(1:ncol,pverp)   = state1%u(1:ncol,pver)
    vm(1:ncol,pverp)   = state1%v(1:ncol,pver)
    thlm(1:ncol,pverp) = thlm(1:ncol,pver)

    if (clubb_do_adv) then
      thlp2(1:ncol,pverp)   = thlp2(1:ncol,pver)
      rtp2(1:ncol,pverp)    = rtp2(1:ncol,pver)
      rtpthlp(1:ncol,pverp) = rtpthlp(1:ncol,pver)
      wpthlp(1:ncol,pverp)  = wpthlp(1:ncol,pver)
      wprtp(1:ncol,pverp)   = wprtp(1:ncol,pver)
      wp2(1:ncol,pverp)     = wp2(1:ncol,pver)
      wp3(1:ncol,pverp)     = wp3(1:ncol,pver)
      up2(1:ncol,pverp)     = up2(1:ncol,pver)
      vp2(1:ncol,pverp)     = vp2(1:ncol,pver)
    endif

    !  Compute virtual potential temperature, which is needed for CLUBB
    do k=1,pver
      do i=1,ncol
        thv(i,k) = state1%t(i,k)*inv_exner_clubb(i,k)*(1._r8+zvir*state1%q(i,k,ixq)&
                   -state1%q(i,k,ixcldliq))
      enddo
    enddo

   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq)

   call tropopause_findChemTrop(state, troplev)

   ! Initialize EDMF outputs
   if (do_clubb_mf) then 
     mf_dry_a_output(:,:)     = 0._r8
     mf_moist_a_output(:,:)   = 0._r8
     mf_dry_w_output(:,:)     = 0._r8
     mf_moist_w_output(:,:)   = 0._r8
     mf_dry_qt_output(:,:)    = 0._r8
     mf_moist_qt_output(:,:)  = 0._r8
     mf_dry_thl_output(:,:)   = 0._r8
     mf_moist_thl_output(:,:) = 0._r8
     mf_dry_u_output(:,:)     = 0._r8
     mf_moist_u_output(:,:)   = 0._r8
     mf_dry_v_output(:,:)     = 0._r8
     mf_moist_v_output(:,:)   = 0._r8
     mf_moist_qc_output(:,:)  = 0._r8
     mf_precc_output(:,:)     = 0._r8
     s_ae_output(:,:)         = 0._r8
     s_aw_output(:,:)         = 0._r8
     s_awthlup_output(:,:)    = 0._r8
     s_awqtup_output(:,:)     = 0._r8
     s_awthldn_output(:,:)    = 0._r8
     s_awqtdn_output(:,:)     = 0._r8
     s_awthl_output(:,:)      = 0._r8
     s_awqt_output(:,:)       = 0._r8
     s_awuup_output(:,:)      = 0._r8
     s_awvup_output(:,:)      = 0._r8
     s_awudn_output(:,:)      = 0._r8
     s_awvdn_output(:,:)      = 0._r8
     s_awu_output(:,:)        = 0._r8
     s_awv_output(:,:)        = 0._r8
     s_aww_output(:,:)        = 0._r8
     mf_upa_output(:,:)       = 0._r8
     mf_upw_output(:,:)       = 0._r8
     mf_upmf_output(:,:)      = 0._r8
     mf_upqt_output(:,:)      = 0._r8
     mf_upthl_output(:,:)     = 0._r8
     mf_upthv_output(:,:)     = 0._r8
     mf_upth_output(:,:)      = 0._r8
     mf_upqc_output(:,:)      = 0._r8
     mf_upbuoy_output(:,:)    = 0._r8
     mf_upent_output(:,:)     = 0._r8
     mf_updet_output(:,:)     = 0._r8
     mf_upa_flip(:,:,:)       = 0._r8
     mf_upw_flip(:,:,:)       = 0._r8
     mf_upmf_flip(:,:,:)      = 0._r8
     mf_upqt_flip(:,:,:)      = 0._r8
     mf_upthl_flip(:,:,:)     = 0._r8
     mf_upthv_flip(:,:,:)     = 0._r8
     mf_upth_flip(:,:,:)      = 0._r8
     mf_upqc_flip(:,:,:)      = 0._r8
     mf_upbuoy_flip(:,:,:)    = 0._r8
     mf_upent_flip(:,:,:)     = 0._r8
     mf_updet_flip(:,:,:)     = 0._r8
     mf_thlflxup_output(:,:)  = 0._r8
     mf_qtflxup_output(:,:)   = 0._r8
     mf_thlflxdn_output(:,:)  = 0._r8
     mf_qtflxdn_output(:,:)   = 0._r8
     mf_thlflx_output(:,:)    = 0._r8
     mf_qtflx_output(:,:)     = 0._r8
     mf_thvflx_output(:,:)    = 0._r8
     mf_uflxup_output(:,:)    = 0._r8
     mf_vflxup_output(:,:)    = 0._r8
     mf_uflxdn_output(:,:)    = 0._r8
     mf_vflxdn_output(:,:)    = 0._r8
     mf_uflx_output(:,:)      = 0._r8
     mf_vflx_output(:,:)      = 0._r8
     mf_thlforcup_output(:,:) = 0._r8
     mf_qtforcup_output(:,:)  = 0._r8
     mf_thlforcdn_output(:,:) = 0._r8
     mf_qtforcdn_output(:,:)  = 0._r8
     mf_thlforc_output(:,:)   = 0._r8
     mf_qtforc_output(:,:)    = 0._r8
     mf_sqtup_output(:,:)     = 0._r8
     mf_sqtdn_output(:,:)     = 0._r8
     mf_rcm_output(:,:)       = 0._r8
     mf_cloudfrac_output(:,:) = 0._r8
     mf_ent_output(:,:)       = 0._r8
     mf_qc_output(:,:)        = 0._r8
     mf_ztop_output(:)        = 0._r8
     mf_ddcp_output(:)        = 0._r8
     mf_L0_output(:)          = 0._r8
     mf_freq_output(:)        = 0._r8
     mf_cape_output(:)        = 0._r8
     mf_cfl_output(:)         = 0._r8
     mf_dnw_output(:,:)       = 0._r8
     mf_dnthl_output(:,:)     = 0._r8
     mf_dnqt_output(:,:)      = 0._r8
     mf_dnw_flip(:,:,:)       = 0._r8
     mf_dnthl_flip(:,:,:)     = 0._r8
     mf_dnqt_flip(:,:,:)      = 0._r8
    end if

    call t_startf("clubb_tend_cam_i_loop")

    ! Determine Coriolis force at given latitude. This is never used
    ! when CLUBB is implemented in a host model, therefore just set
    ! to zero.
    fcor(:) = 0._r8

    ! Define the CLUBB momentum grid (in height, units of m)
    do k=1, nlev+1
      do i=1, ncol
        zi_g(i,k) = state1%zi(i,pverp-k+1)-state1%zi(i,pver+1)
      end do
    end do

    ! Define the CLUBB thermodynamic grid (in units of m)
    do k=1, nlev
      do i=1, ncol
        zt_g(i,k+1) = state1%zm(i,pver-k+1)-state1%zi(i,pver+1)
      end do
    end do

    do k=1, pver
      do i=1, ncol
        dz_g(i,k) = state1%zi(i,k)-state1%zi(i,k+1)  ! compute thickness
      end do
    end do

    !  Thermodynamic ghost point is below surface
    do i=1, ncol
      zt_g(i,1) = -1._r8*zt_g(i,2)
    end do

    do i=1, ncol
      !  Set the elevation of the surface
      sfc_elevation(i) = state1%zi(i,pver+1)
    end do

    !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
    !  Inputs for the momentum levels are set below setup_clubb core
    do k=1,nlev
      do i=1, ncol
        ! base state (dry) variables
        rho_ds_zt(i,k+1)       = rga*(state1%pdeldry(i,pver-k+1)/dz_g(i,pver-k+1))
        invrs_rho_ds_zt(i,k+1) = 1._r8/(rho_ds_zt(i,k+1))

        ! full state (moist) variables
        p_in_Pa(i,k+1)         = state1%pmid(i,pver-k+1)
        exner(i,k+1)           = 1._r8/inv_exner_clubb(i,pver-k+1)
        thv(i,k+1)             = state1%t(i,pver-k+1)*inv_exner_clubb(i,pver-k+1)*(1._r8+zvir*state1%q(i,pver-k+1,ixq) &
                                   -state1%q(i,pver-k+1,ixcldliq))
        rho_zt(i,k+1)          = rga*state1%pdel(i,pver-k+1)/dz_g(i,pver-k+1)

        ! exception - setting this to moist thv
        thv_ds_zt(i,k+1)       = thv(i,k+1)

        rfrzm(i,k+1)           = state1%q(i,pver-k+1,ixcldice)
        radf(i,k+1)            = radf_clubb(i,pver-k+1)
        qrl_clubb(i,k+1)       = qrl(i,pver-k+1)/(cpairv(i,k,lchnk)*state1%pdeldry(i,pver-k+1))
      end do
    end do

    !  Compute mean w wind on thermo grid, convert from omega to w
    do k=1,nlev
      do i=1,ncol
        wm_zt(i,k+1) = -1._r8*(state1%omega(i,pver-k+1)-state1%omega(i,pver))/(rho_zt(i,k+1)*gravit)
      end do
    end do

    !  Below computes the same stuff for the ghost point.  May or may
    !  not be needed, just to be safe to avoid NaN's
    thv_ds_zt(:ncol,1)       = thv_ds_zt(:ncol,2)
    rho_ds_zt(:ncol,1)       = rho_ds_zt(:ncol,2)
    invrs_rho_ds_zt(:ncol,1) = invrs_rho_ds_zt(:ncol,2)
    p_in_Pa(:ncol,1)         = p_in_Pa(:ncol,2)
    exner(:ncol,1)           = exner(:ncol,2)
    thv(:ncol,1)             = thv(:ncol,2)
    rho_zt(:ncol,1)          = rho_zt(:ncol,2)
    rfrzm(:ncol,1)           = rfrzm(:ncol,2)
    radf(:ncol,1)            = radf(:ncol,2)
    qrl_clubb(:ncol,1)       = qrl_clubb(:ncol,2)
    wm_zt(:ncol,1)           = wm_zt(:ncol,2)

    
    ! ------------------------------------------------- !
    ! Begin case specific code for SCAM cases.          !
    ! This section of code block is NOT called in       !
    ! global simulations                                !
    ! ------------------------------------------------- !
    if (single_column) then

      !  Initialize zo if variable ustar is used
      if (cam_in%landfrac(1) >= 0.5_r8) then
        zo(1) = 0.035_r8
      else
        zo(1) = 0.0001_r8
      endif

      !  Compute surface wind (ubar)
      ubar = sqrt(um(1,pver)**2+vm(1,pver)**2)
      if (ubar <  0.25_r8) ubar = 0.25_r8
  
      !  Below denotes case specifics for surface momentum
      !  and thermodynamic fluxes, depending on the case

      !  Define ustar (based on case, if not variable)
      ustar = 0.25_r8   ! Initialize ustar in case no case

      if(trim(scm_clubb_iop_name)  ==  'BOMEX_5day') then
        ustar = 0.28_r8
      endif

      if(trim(scm_clubb_iop_name)  ==  'ATEX_48hr') then
        ustar = 0.30_r8
      endif

      if(trim(scm_clubb_iop_name)  ==  'RICO_3day') then
        ustar      = 0.28_r8
      endif

      if(trim(scm_clubb_iop_name)  ==  'arm97' .or. trim(scm_clubb_iop_name)  ==  'gate' .or. &
         trim(scm_clubb_iop_name)  ==  'toga' .or. trim(scm_clubb_iop_name)  ==  'mpace' .or. &
         trim(scm_clubb_iop_name)  ==  'ARM_CC') then

          bflx22(1) = (gravit/theta0)*wpthlp_sfc(1)
          ustar  = diag_ustar(zt_g(1,2),bflx22(1),ubar,zo(1))
      endif

      !  Compute the surface momentum fluxes, if this is a SCAM simulation
      upwp_sfc(1) = -um(1,pver)*ustar**2/ubar
      vpwp_sfc(1) = -vm(1,pver)*ustar**2/ubar

    end if

    !  Define surface sources for transported variables for diffusion, will
    !  be zero as these tendencies are done in vertical_diffusion
    do ixind=1,edsclr_dim
      do i=1,ncol
        wpedsclrp_sfc(i,ixind) = 0._r8
      end do
    end do

    !  Set stats output and increment equal to CLUBB and host dt
    stats_tsamp = dtime
    stats_tout  = hdtime

    stats_nsamp = nint(stats_tsamp/dtime)
    stats_nout = nint(stats_tout/dtime)

    !  Heights need to be set at each timestep.  Therefore, recall
    !  setup_grid and setup_parameters for this.

    !  Set-up CLUBB core at each CLUBB call because heights can change
    !  Important note:  do not make any calls that use CLUBB grid-height
    !                   operators (such as zt2zm_api, etc.) until AFTER the
    !                   call to setup_grid_heights_api.
    call setup_grid_api( nlev+1, ncol, sfc_elevation, l_implemented,      & ! intent(in)
                         grid_type, zi_g(:,2), zi_g(:,1), zi_g(:,nlev+1), & ! intent(in)
                         zi_g, zt_g,                                      & ! intent(in)
                         gr, begin_height, end_height )                     ! intent(out)

    call setup_parameters_api( zi_g(:,2), clubb_params, nlev+1, ncol, grid_type,  & ! intent(in)
                               zi_g, zt_g,                                        & ! intent(in)
                               clubb_config_flags%l_prescribed_avg_deltaz,        & ! intent(in)
                               lmin, nu_vert_res_dep, err_code )                    ! intent(out)
    if ( err_code == clubb_fatal_error ) then
       call endrun(subr//':  Fatal error in CLUBB setup_parameters')
    end if


    !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
    !  forcings already applied through CAM
    thlm_forcing(:,:) = 0._r8
    rtm_forcing(:,:) = 0._r8
    um_forcing(:,:) = 0._r8
    vm_forcing(:,:) = 0._r8


    rtm_ref(:,:) = 0.0_r8
    thlm_ref(:,:) = 0.0_r8
    um_ref(:,:) = 0.0_r8
    vm_ref(:,:) = 0.0_r8
    ug(:,:) = 0.0_r8
    vg(:,:) = 0.0_r8
    
    ! Add forcings for SILHS covariance contributions
    rtp2_forcing    = zt2zm_api( pverp+1-top_lev, ncol, gr, rtp2_mc_zt(1:ncol,:) )
    thlp2_forcing   = zt2zm_api( pverp+1-top_lev, ncol, gr, thlp2_mc_zt(1:ncol,:) )
    wprtp_forcing   = zt2zm_api( pverp+1-top_lev, ncol, gr, wprtp_mc_zt(1:ncol,:) )
    wpthlp_forcing  = zt2zm_api( pverp+1-top_lev, ncol, gr, wpthlp_mc_zt(1:ncol,:) )
    rtpthlp_forcing = zt2zm_api( pverp+1-top_lev, ncol, gr, rtpthlp_mc_zt(1:ncol,:) )

    ! Zero out SILHS covariance contribution terms
    rtp2_mc_zt(:,:) = 0.0_r8
    thlp2_mc_zt(:,:) = 0.0_r8
    wprtp_mc_zt(:,:) = 0.0_r8
    wpthlp_mc_zt(:,:) = 0.0_r8
    rtpthlp_mc_zt(:,:) = 0.0_r8


    !  Compute some inputs from the thermodynamic grid
    !  to the momentum grid
    rho_ds_zm       = zt2zm_api( pverp+1-top_lev, ncol, gr, rho_ds_zt(1:ncol,:) )
    rho_zm          = zt2zm_api( pverp+1-top_lev, ncol, gr, rho_zt(1:ncol,:) )
    invrs_rho_ds_zm = zt2zm_api( pverp+1-top_lev, ncol, gr, invrs_rho_ds_zt(1:ncol,:) )
    thv_ds_zm       = zt2zm_api( pverp+1-top_lev, ncol, gr, thv_ds_zt(1:ncol,:) )
    wm_zm           = zt2zm_api( pverp+1-top_lev, ncol, gr, wm_zt(1:ncol,:) )

    !  Surface fluxes provided by host model
    do i=1,ncol
      wpthlp_sfc(i) = cam_in%shf(i)/(cpairv(i,pver,lchnk)*rho_ds_zm(i,1)) ! Sensible heat flux
      wpthlp_sfc(i) = wpthlp_sfc(i)*inv_exner_clubb_surf(i)   ! Potential temperature flux
      wprtp_sfc(i)  = cam_in%cflx(i,1)/rho_ds_zm(i,1)         ! Moisture flux
    end do

    ! Implementation after Thomas Toniazzo (NorESM) and Colin Zarzycki (PSU)
    !  Other Surface fluxes provided by host model
    if( (cld_macmic_num_steps > 1) .and. clubb_l_intr_sfc_flux_smooth ) then
       ! Adjust surface stresses using winds from the prior macmic iteration
       do i=1,ncol
          ubar = sqrt(state1%u(i,pver)**2+state1%v(i,pver)**2)
          if (ubar <  0.25_r8) ubar = 0.25_r8

          call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
               rrho(i), ustar )

          upwp_sfc(i) = -state1%u(i,pver)*ustar**2/ubar
          vpwp_sfc(i) = -state1%v(i,pver)*ustar**2/ubar
       end do
    else
       do i=1,ncol
          upwp_sfc(i)   = cam_in%wsx(i)/rho_ds_zm(i,1)               ! Surface meridional momentum flux
          vpwp_sfc(i)   = cam_in%wsy(i)/rho_ds_zm(i,1)               ! Surface zonal momentum flux
       end do
    endif

    ! Perturbed winds are not used in CAM
    upwp_sfc_pert = 0.0_r8
    vpwp_sfc_pert = 0.0_r8

    !  Need to flip arrays around for CLUBB core
    do k=1,nlev+1
      do i=1,ncol
        um_in(i,k)      = um(i,pverp-k+1)
        vm_in(i,k)      = vm(i,pverp-k+1)
        upwp_in(i,k)    = upwp(i,pverp-k+1)
        vpwp_in(i,k)    = vpwp(i,pverp-k+1)
        wpthvp_in(i,k)  = wpthvp(i,pverp-k+1)
        wp2thvp_in(i,k) = wp2thvp(i,pverp-k+1)
        rtpthvp_in(i,k) = rtpthvp(i,pverp-k+1)
        thlpthvp_in(i,k)= thlpthvp(i,pverp-k+1)
        up2_in(i,k)     = up2(i,pverp-k+1)
        vp2_in(i,k)     = vp2(i,pverp-k+1)
        up3_in(i,k)     = up3(i,pverp-k+1)
        vp3_in(i,k)     = vp3(i,pverp-k+1)
        wp2_in(i,k)     = wp2(i,pverp-k+1)
        wp3_in(i,k)     = wp3(i,pverp-k+1)
        rtp2_in(i,k)    = rtp2(i,pverp-k+1)
        thlp2_in(i,k)   = thlp2(i,pverp-k+1)
        rtp3_in(i,k)    = rtp3(i,pverp-k+1)
        thlp3_in(i,k)   = thlp3(i,pverp-k+1)
        thlm_in(i,k)    = thlm(i,pverp-k+1)
        rtm_in(i,k)     = rtm(i,pverp-k+1)
        rvm_in(i,k)     = rvm(i,pverp-k+1)
        wprtp_in(i,k)   = wprtp(i,pverp-k+1)
        wpthlp_in(i,k)  = wpthlp(i,pverp-k+1)
        tke_in(i,k)     = tke(i,pverp-k+1)
        rtpthlp_in(i,k) = rtpthlp(i,pverp-k+1)
        cloud_frac_inout(i,k) = cloud_frac(i,pverp-k+1)
        if (k>1) then
          rcm_inout(i,k) = state1%q(i,pverp-k+1,ixcldliq)
        end if

        ! We only need to copy pdf_params from pbuf if this is a restart and
        ! we're calling pdf_closure at the end of advance_clubb_core
        if ( is_first_restart_step() &
             .and. clubb_config_flags%ipdf_call_placement .eq. ipdf_post_advance_fields ) then
          pdf_params_zm_chnk(lchnk)%w_1(i,k)        = pdf_zm_w_1(i,pverp-k+1)
          pdf_params_zm_chnk(lchnk)%w_2(i,k)        = pdf_zm_w_2(i,pverp-k+1)
          pdf_params_zm_chnk(lchnk)%varnce_w_1(i,k) = pdf_zm_varnce_w_1(i,pverp-k+1)
          pdf_params_zm_chnk(lchnk)%varnce_w_2(i,k) = pdf_zm_varnce_w_2(i,pverp-k+1)
          pdf_params_zm_chnk(lchnk)%mixt_frac(i,k)  = pdf_zm_mixt_frac(i,pverp-k+1)
        end if

        sclrpthvp_inout(i,k,:) = 0._r8
        wp2rtp_inout(i,k)  = wp2rtp(i,pverp-k+1)
        wp2thlp_inout(i,k) = wp2thlp(i,pverp-k+1)
        uprcp_inout(i,k)   = uprcp(i,pverp-k+1)
        vprcp_inout(i,k)   = vprcp(i,pverp-k+1)
        rc_coef_inout(i,k) = rc_coef(i,pverp-k+1)
        wp4_inout(i,k)     = wp4(i,pverp-k+1)
        wpup2_inout(i,k)   = wpup2(i,pverp-k+1)
        wpvp2_inout(i,k)   = wpvp2(i,pverp-k+1)
        wp2up2_inout(i,k)  = wp2up2(i,pverp-k+1)
        wp2vp2_inout(i,k)  = wp2vp2(i,pverp-k+1)
        ice_supersat_frac_inout(i,k) = ice_supersat_frac(i,pverp-k+1)
      end do
    end do

    ! Perturbed winds are not used in CAM
    um_pert_inout   = 0.0_r8
    vm_pert_inout   = 0.0_r8
    upwp_pert_inout = 0.0_r8
    vpwp_pert_inout = 0.0_r8

    do k=2,nlev+1
      do i=1,ncol
        pre_in(i,k) = prer_evap(i,pverp-k+1)
      end do
    end do

    do i=1,ncol
      pre_in(i,1) = pre_in(i,2)
    end do

    do i=1,ncol
      rcm_inout(i,1) = rcm_inout(i,2)
    end do

    !  Initialize these to prevent crashing behavior
    do k=1,nlev+1
      do i=1,ncol
        wprcp_out(i,k)        = 0._r8
        rcm_in_layer_out(i,k) = 0._r8
        cloud_cover_out(i,k)  = 0._r8
        edsclr_in(i,k,:)      = 0._r8
        khzm_out(i,k)         = 0._r8
        khzt_out(i,k)         = 0._r8
      end do
    end do

    !  higher order scalar stuff, put to zero
    do ixind=1, sclr_dim
      do k=1, nlev+1
        do i=1, ncol
          sclrm(i,k,ixind)     = 0._r8
          wpsclrp(i,k,ixind)   = 0._r8
          sclrp2(i,k,ixind)    = 0._r8
          sclrp3(i,k,ixind)    = 0._r8
          sclrprtp(i,k,ixind)  = 0._r8
          sclrpthlp(i,k,ixind) = 0._r8
          wpsclrp_sfc(i,ixind) = 0._r8
        end do
      end do
    end do

    do ixind=1, hydromet_dim
      do k=1, nlev+1
        do i=1, ncol
          hydromet(i,k,ixind)    = 0._r8
          wphydrometp(i,k,ixind) = 0._r8
          wp2hmp(i,k,ixind)      = 0._r8
          rtphmp_zt(i,k,ixind)   = 0._r8
          thlphmp_zt(i,k,ixind)  = 0._r8
        end do
      end do
    end do

    ! pressure,exner on momentum grid needed for mass flux calc.
    if (do_clubb_mf) then

      kappa_zt(:,:) = 0._r8
      qc_zt(:,:) = 0._r8
      invrs_exner_zt(:,:) = 0._r8
      do k=1,pver
        do i=1,ncol
          kappa_zt(i,k+1) = (rairv(i,pver-k+1,lchnk)/cpairv(i,pver-k+1,lchnk))
          qc_zt(i,k+1) = state1%q(i,pver-k+1,ixcldliq)
          invrs_exner_zt(i,k+1) = inv_exner_clubb(i,pver-k+1)
        end do
      end do

      kappa_zt(:ncol,1) = kappa_zt(:ncol,2)
      qc_zt(:ncol,1) = qc_zt(:ncol,2)
      invrs_exner_zt(:ncol,1) = invrs_exner_zt(:ncol,2)

      kappa_zm(:,:) = 0._r8
      kappa_zm(1:ncol,:) = zt2zm_api(pverp+1-top_lev, ncol, gr, kappa_zt(1:ncol,:))

      p_in_Pa_zm(:,:) = 0._r8
      invrs_exner_zm(:,:) = 0._r8
      do k=1,pverp
        do i=1,ncol
          p_in_Pa_zm(i,k) = state1%pint(i,pverp-k+1)
          invrs_exner_zm(i,k) = 1._r8/((p_in_Pa_zm(i,k)/p0_clubb)**(kappa_zm(i,k)))
        end do
      end do

      th_sfc(:) = 0._r8
      th_sfc(1:ncol) = cam_in%ts(1:ncol)*invrs_exner_zm(1:ncol,1)

      call calc_ustar( ncol, state1%t(:ncol,pver), state1%pmid(:ncol,pver), cam_in%wsx(:ncol), cam_in%wsy(:ncol), &
                       rrho(:ncol), ustar2(:ncol) )
    end if

    if (clubb_do_adv) then
      if (macmic_it  ==  1) then

        wp2_in     = zt2zm_api(pverp+1-top_lev, ncol, gr, wp2_in )
        wpthlp_in  = zt2zm_api(pverp+1-top_lev, ncol, gr, wpthlp_in )
        wprtp_in   = zt2zm_api(pverp+1-top_lev, ncol, gr, wprtp_in )
        up2_in     = zt2zm_api(pverp+1-top_lev, ncol, gr, up2_in )
        vp2_in     = zt2zm_api(pverp+1-top_lev, ncol, gr, vp2_in )
        thlp2_in   = zt2zm_api(pverp+1-top_lev, ncol, gr, thlp2_in )
        rtp2_in    = zt2zm_api(pverp+1-top_lev, ncol, gr, rtp2_in )
        rtpthlp_in = zt2zm_api(pverp+1-top_lev, ncol, gr, rtpthlp_in )

        do k=1,nlev+1
          do i=1,ncol
            thlp2_in(i,k) = max(thl_tol**2,thlp2_in(i,k))
            rtp2_in(i,k)  = max(rt_tol**2,rtp2_in(i,k))
            wp2_in(i,k)   = max(w_tol_sqd,wp2_in(i,k))
            up2_in(i,k)   = max(w_tol_sqd,up2_in(i,k))
            vp2_in(i,k)   = max(w_tol_sqd,vp2_in(i,k))
          end do
        end do

      end if
    end if

    !  Do the same for tracers
    icnt=0
    do ixind=1,pcnst
      if (lq(ixind))  then

        icnt = icnt+1

        do k=1,nlev
          do i=1,ncol
            edsclr_in(i,k+1,icnt) = state1%q(i,pver-k+1,ixind)
          end do
        end do

        do i=1,ncol
          edsclr_in(i,1,icnt) = edsclr_in(i,2,icnt)
        end do

      end if
    end do


    if (clubb_l_do_expldiff_rtm_thlm) then
      do k=1,nlev
        do i=1, ncol
          edsclr_in(i,k+1,icnt+1) = thlm(i,pver-k+1)
          edsclr_in(i,k+1,icnt+2) = rtm(i,pver-k+1)
        end do
      end do

      edsclr_in(:ncol,1,icnt+1) = edsclr_in(:ncol,2,icnt+1)
      edsclr_in(:ncol,1,icnt+2) = edsclr_in(:ncol,2,icnt+2)

    endif

    if (macmic_it==1) thlm_macmic1(:ncol,:) = 0._r8
    if (macmic_it==1) thlm_macmic2(:ncol,:) = 0._r8
    if (macmic_it==1) qtm_macmic1(:ncol,:) = 0._r8
    if (macmic_it==1) qtm_macmic2(:ncol,:) = 0._r8

    if (do_clubb_mf) then
      mf_L0          = 0._r8
      mf_L0_nadv     = 0._r8
      mf_ztop        = 0._r8
      mf_ztop_nadv   = 0._r8
      mf_ztopm1      = 0._r8
      mf_ztopm1_nadv = 0._r8
      mf_ddcp_nadv   = 0._r8
      mf_cbm1        = 0._r8
      mf_cbm1_nadv   = 0._r8
      mf_freq_nadv   = 0._r8

      if (macmic_it==1) ztopm1_macmic(:ncol) = 0._r8
      if (macmic_it==1) ddcp_macmic(:ncol) = 0._r8
      if (macmic_it==1) cbm1_macmic(:ncol) = 0._r8

      if (macmic_it==1) up_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) up_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dn_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dn_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) upa_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) upa_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dna_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dna_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) thlu_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) thlu_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) qtu_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) qtu_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) thld_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) thld_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) qtd_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) qtd_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dthl_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dthl_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dqt_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dqt_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dthlu_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dthlu_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dqtu_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dqtu_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dthld_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dthld_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) dqtd_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) dqtd_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) ztop_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) ztop_macmic2(:ncol,:) = 0._r8
      if (macmic_it==1) ddcp_macmic1(:ncol,:) = 0._r8
      if (macmic_it==1) ddcp_macmic2(:ncol,:) = 0._r8

!+++ARH - Temporary hack - pbuf_set_field is apparently not taking?
      if (is_first_step() .and. macmic_it==1) then
        ddcp(:ncol) = 0._r8
      end if

      mf_precc_nadv(:ncol)      = 0._r8
      mf_snow_nadv(:ncol)       = 0._r8

      mf_qc(:ncol,:pverp)             = 0._r8
      mf_rcm(:ncol,:pverp)            = 0._r8
      mf_cloudfrac(:ncol,:pverp)      = 0._r8
      mf_qc_nadv(:ncol,:pverp)        = 0._r8
      mf_rcm_nadv(:ncol,:pverp)       = 0._r8
      mf_cloudfrac_nadv(:ncol,:pverp) = 0._r8

      mf_thlforcup_nadv(:ncol,:pverp) = 0._r8
      mf_qtforcup_nadv(:ncol,:pverp)  = 0._r8
      mf_thlforcdn_nadv(:ncol,:pverp) = 0._r8
      mf_qtforcdn_nadv(:ncol,:pverp)  = 0._r8
      mf_thlforc_nadv(:ncol,:pverp)   = 0._r8
      mf_qtforc_nadv(:ncol,:pverp)    = 0._r8
      mf_ent_nadv(:ncol,:pverp)       = 0._r8

      max_cfl_nadv(:ncol) = 0._r8
    end if

    do t=1,nadv    ! do needed number of "sub" timesteps for each CAM step
    
      !  Increment the statistics then being stats timestep
      if (l_stats) then
         call stats_begin_timestep_api(t, stats_nsamp, stats_nout)
      endif

      !#######################################################################
      !###################### CALL MF DIAGNOSTIC PLUMES ######################
      !#######################################################################
      if (do_clubb_mf) then


        dzt(:,:) = 0._r8
        do k=2,pverp
          do i=1, ncol
            dzt(i,k) = zi_g(i,k) - zi_g(i,k-1)
          end do
        end do

        dzt(:ncol,1) = dzt(:ncol,2)

        invrs_dzt(:,:) = 0._r8
        invrs_dzt(:ncol,:) = 1._r8/dzt(:ncol,:)

        rtm_zm_in(:,:) = 0._r8
        thlm_zm_in(:,:) = 0._r8
        th_zm(:,:) = 0._r8
        qv_zm(:,:) = 0._r8
        qc_zm(:,:) = 0._r8
        rtm_zm_in(1:ncol,:)  = zt2zm_api( pverp+1-top_lev, ncol, gr, rtm_in(1:ncol,:) )
        thlm_zm_in(1:ncol,:) = zt2zm_api( pverp+1-top_lev, ncol, gr, thlm_in(1:ncol,:) )
        th_zm(1:ncol,:)      = zt2zm_api( pverp+1-top_lev, ncol, gr, th_zt(1:ncol,:) )
        qv_zm(1:ncol,:)      = zt2zm_api( pverp+1-top_lev, ncol, gr, qv_zt(1:ncol,:) )
        qc_zm(1:ncol,:)      = zt2zm_api( pverp+1-top_lev, ncol, gr, qc_zt(1:ncol,:) )

        if (t>1) then
          ! update thv if clubb is subcycled
          do i=1, ncol
            thv_ds_zt(i,:pverp) = thlm_in(i,:pverp) &
                                + latvap*rcm_inout(i,:pverp)*invrs_exner_zt(i,:pverp)/cpair
            thv_ds_zt(i,:pverp) = thv_ds_zt(i,:pverp) &
                                * (1._r8+zvir*(rtm_in(i,:pverp)-rcm_inout(i,:pverp)) - rcm_inout(i,:pverp))
          end do
          thv_ds_zm(1:ncol,:) = zt2zm_api( pverp+1-top_lev, ncol, gr, thv_ds_zt(1:ncol,:) )
        end if

        mf_ztopm1(1:ncol) = ztopma(1:ncol) 
        mf_ddcp(1:ncol) = ddcp(1:ncol)
        mf_cbm1(1:ncol) = cbm1(1:ncol)

        rhinv(1:ncol) = 0._r8
        do i=1, ncol
          if (rhlev(i) >= 1._r8) rhlev(i) = 0.990_r8
          if (rhlev(i) > 0._r8) rhinv(i) = 1._r8 / ( (1._r8/rhlev(i)) - 1._r8 )
        end do

        do i=1, ncol
          call integrate_mf( pverp,                                                                                   & ! input
                             rho_zm(i,:),    dzm(i,:),         zi_g(i,:),       p_in_Pa_zm(i,:), invrs_exner_zm(i,:), & ! input
                             rho_zt(i,:),    dzt(i,:),         zt_g(i,:),       p_in_Pa(i,:),    invrs_exner_zt(i,:), & ! input
                             um_in(i,:),     vm_in(i,:),       thlm_in(i,:),    rtm_in(i,:),     thv_ds_zt(i,:),      & ! input
                                                               th_zt(i,:),      qv_zt(i,:),      qc_zt(i,:),          & ! input
                                                               thlm_zm_in(i,:), rtm_zm_in(i,:),  thv_ds_zm(i,:),      & ! input
                                                               th_zm(i,:),      qv_zm(i,:),      qc_zm(i,:),          & ! input
                             ustar2(i),      th_sfc(i),        wpthlp_sfc(i),   wprtp_sfc(i),    pblh(i),             & ! input
                             wpthlp_in(i,:), tke_in(i,:),      tpert(i),        mf_ztopm1(i),    rhinv(i),            & ! input                     
                                             mf_cape_output(i),mf_ddcp(i),      mf_cbm1(i),                           & ! output - plume diagnostics
                             mf_upa(i,:,:),    mf_dna(i,:,:),                                                         & ! output - plume diagnostics
                             mf_upw(i,:,:),    mf_dnw(i,:,:),                                                         & ! output - plume diagnostics
                             mf_upmf(i,:,:),                                                                          & ! output - plume diagnostics
                             mf_upqt(i,:,:),   mf_dnqt(i,:,:),                                                        & ! output - plume diagnostics
                             mf_upthl(i,:,:),  mf_dnthl(i,:,:),                                                       & ! output - plume diagnostics
                             mf_upthv(i,:,:),  mf_dnthv(i,:,:),                                                       & ! output - plume diagnostics
                             mf_upth(i,:,:),   mf_dnth(i,:,:),                                                        & ! output - plume diagnostics
                             mf_upqc(i,:,:),   mf_dnqc(i,:,:),                                                        & ! output - plume diagnostics
                             mf_upbuoy(i,:,:),                                                                        & ! output - plume diagnostics
                             mf_upent(i,:,:),                                                                         & ! output - plume diagnostics
                             mf_updet(i,:,:),                                                                         & ! output - plume diagnostics
                             mf_dry_a(i,:),    mf_moist_a(i,:),                                                       & ! output - plume diagnostics
                             mf_dry_w(i,:),    mf_moist_w(i,:),                                                       & ! output - plume diagnostics
                             mf_dry_qt(i,:),   mf_moist_qt(i,:),                                                      & ! output - plume diagnostics
                             mf_dry_thl(i,:),  mf_moist_thl(i,:),                                                     & ! output - plume diagnostics
                             mf_dry_u(i,:),    mf_moist_u(i,:),                                                       & ! output - plume diagnostics
                             mf_dry_v(i,:),    mf_moist_v(i,:),                                                       & ! output - plume diagnostics
                                               mf_moist_qc(i,:),                                                      & ! output - plume diagnostics
                             s_ae(i,:),                                                                               & ! output - plume diagnostics
                             s_ac(i,:),        s_aup(i,:),      s_adn(i,:),                                           & ! output - plume diagnostics
                             s_aw(i,:),        s_awup(i,:),     s_awdn(i,:),                                          & ! output - plume diagnostics
                             s_aww(i,:),       s_awwup(i,:),    s_awwdn(i,:),                                         & ! output - plume diagnostics
                             s_awthlup(i,:),   s_awqtup(i,:),   s_awuup(i,:),   s_awvup(i,:),                         & ! output - plume diagnostics
                             s_awthldn(i,:),   s_awqtdn(i,:),   s_awudn(i,:),   s_awvdn(i,:),                         & ! output - plume diagnostics
                             s_awthl(i,:),     s_awqt(i,:),                                                           & ! output - plume diagnostics
                             s_awu(i,:),       s_awv(i,:),                                                            & ! output - plume diagnostics
                             mf_thlflxup(i,:), mf_qtflxup(i,:), mf_uflxup(i,:), mf_vflxup(i,:),                       & ! output - plume diagnostics
                             mf_thlflxdn(i,:), mf_qtflxdn(i,:), mf_uflxdn(i,:), mf_vflxdn(i,:),                       & ! output - plume diagnostics
                             mf_thlflx(i,:),   mf_qtflx(i,:),   mf_uflx(i,:),   mf_vflx(i,:),                         & ! output - variables needed for solver
                             mf_thvflx(i,:),                                                                          & ! output - plume diagnostics
                             mf_sqtup(i,:),    mf_sthlup(i,:),                                                        & ! output - plume diagnostics
                             mf_sqtdn(i,:),    mf_sthldn(i,:),                                                        & ! output - plume diagnostics
                             mf_sqt(i,:),      mf_sthl(i,:),                                                          & ! output - variables needed for solver
                             mf_precc(i,:),                                                                           & ! output - plume diagnostics
                             mf_ztop(i),       mf_L0(i) )
        end do

        ! CFL limiter
        cfllim = .true.
        cflfac(:ncol) = 1._r8
        s_aw(:ncol,1) = 0._r8
        max_cfl(:ncol)= 0._r8
        do i=1,ncol
          do k=2,pverp
             max_cfl(i) = max(max_cfl(i),dtime*invrs_dzt(i,k)*max(abs(s_aw(i,k-1)),abs(s_aw(i,k))))
          end do
          if (max_cfl(i).gt.cflval.and.cfllim) cflfac(i) = cflval/max_cfl(i)
        end do

           ! Scale microphys so it can't drive qt negative
        do k=2,pverp
          do i=1,ncol
            if ((-1._r8*mf_sqt(i,k)*dtime) > rtm_in(i,k)) then
              lambda = -1._r8*rtm_in(i,k)/(mf_sqt(i,k)*dtime)
              mf_sqt(i,k) = lambda*mf_sqt(i,k)
              mf_sthl(i,k) = lambda*mf_sthl(i,k)
              mf_sqtup(i,k) = lambda*mf_sqtup(i,k)
              mf_sthlup(i,k) = lambda*mf_sthlup(i,k)
              mf_sqtdn(i,k) = lambda*mf_sqtdn(i,k)
              mf_sthldn(i,k) = lambda*mf_sthldn(i,k)
            end if
          end do
        end do

        ! Recalculate precip using new microphys forcing
        mf_precc(:ncol,:pverp) = 0._r8
        do k=pverp,2,-1
          do i=1,ncol
            mf_precc(i,k-1) = mf_precc(i,k) - rho_zt(i,k)*dzt(i,k)*mf_sqt(i,k)
          end do
        end do



        ! pass MF turbulent advection term as CLUBB explicit forcing term
        rtm_forcing(:ncol,:)  = 0._r8
        thlm_forcing(:ncol,:) = 0._r8
        mf_qtforcup(:ncol,:)  = 0._r8
        mf_thlforcup(:ncol,:) = 0._r8
        mf_qtforcdn(:ncol,:)  = 0._r8
        mf_thlforcdn(:ncol,:) = 0._r8

        do k=2,pverp
          do i=1, ncol
            rtm_forcing(i,k)  = rtm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                              ((rho_ds_zm(i,k) * mf_qtflx(i,k)) - (rho_ds_zm(i,k-1) * mf_qtflx(i,k-1))) &
                               + mf_sqt(i,k)

            thlm_forcing(i,k) = thlm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                               ((rho_ds_zm(i,k) * mf_thlflx(i,k)) - (rho_ds_zm(i,k-1) * mf_thlflx(i,k-1))) &
                               + mf_sthl(i,k)

            mf_qtforcup(i,k)  = mf_qtforcup(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                               ((rho_ds_zm(i,k) * mf_qtflxup(i,k)) - (rho_ds_zm(i,k-1) * mf_qtflxup(i,k-1))) &
                               + mf_sqtup(i,k)

            mf_thlforcup(i,k) = mf_thlforcup(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                               ((rho_ds_zm(i,k) * mf_thlflxup(i,k)) - (rho_ds_zm(i,k-1) * mf_thlflxup(i,k-1))) &
                               + mf_sthlup(i,k)

            mf_qtforcdn(i,k)  = mf_qtforcdn(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                               ((rho_ds_zm(i,k) * mf_qtflxdn(i,k)) - (rho_ds_zm(i,k-1) * mf_qtflxdn(i,k-1))) &
                               + mf_sqtdn(i,k)

            mf_thlforcdn(i,k) = mf_thlforcdn(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                               ((rho_ds_zm(i,k) * mf_thlflxdn(i,k)) - (rho_ds_zm(i,k-1) * mf_thlflxdn(i,k-1))) &
                               + mf_sthldn(i,k)

          end do
        end do

        if (do_clubb_mf_cmt) then
          ! convective momentum transport
          um_forcing(:ncol,:) = 0._r8
          vm_forcing(:ncol,:) = 0._r8
          do k=2,pverp
            do i=1, ncol
              um_forcing(i,k)   = um_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                                 ((rho_ds_zm(i,k) * mf_uflx(i,k)) - (rho_ds_zm(i,k-1) * mf_uflx(i,k-1)))

              vm_forcing(i,k)   = vm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * cflfac(i) * &
                                 ((rho_ds_zm(i,k) * mf_vflx(i,k)) - (rho_ds_zm(i,k-1) * mf_vflx(i,k-1)))
            end do
          end do
        end if

        do i=1,ncol
           ! compute ensemble cloud properties
           mf_qc_nadv(i,:pverp)        = mf_qc_nadv(i,:pverp) + mf_moist_qc(i,:pverp)
           mf_rcm_nadv(i,:pverp)       = mf_rcm_nadv(i,:pverp) + mf_moist_a(i,:pverp)*mf_moist_qc(i,:pverp)
           mf_cloudfrac_nadv(i,:pverp) = mf_cloudfrac_nadv(i,:pverp) + mf_moist_a(i,:pverp)

           ! [kg/m2/s]->[m/s]
           mf_precc_nadv(i) = mf_precc_nadv(i) + mf_precc(i,1)/1000._r8
           mf_snow_nadv(i)  = 0._r8

           ! accumulate over nadv subcycles
           mf_L0_nadv(i)     = mf_L0_nadv(i) + mf_L0(i)
           mf_ztop_nadv(i)   = mf_ztop_nadv(i) + mf_ztop(i)
           mf_ztopm1_nadv(i) = mf_ztopm1_nadv(i) + mf_ztopm1(i) 
           mf_ddcp_nadv(i)   = mf_ddcp_nadv(i) + mf_ddcp(i)
           mf_cbm1_nadv(i)   = mf_cbm1_nadv(i) + mf_cbm1(i)

           if (mf_ztop(i) > 0._r8) mf_freq_nadv(i) = mf_freq_nadv(i) + 1._r8

           mf_thlforcup_nadv(i,:pverp) = mf_thlforcup_nadv(i,:pverp) + mf_thlforcup(i,:pverp)
           mf_qtforcup_nadv(i,:pverp)  = mf_qtforcup_nadv(i,:pverp) + mf_qtforcup(i,:pverp)
           mf_thlforcdn_nadv(i,:pverp) = mf_thlforcdn_nadv(i,:pverp) + mf_thlforcdn(i,:pverp)
           mf_qtforcdn_nadv(i,:pverp)  = mf_qtforcdn_nadv(i,:pverp) + mf_qtforcdn(i,:pverp)
           mf_thlforc_nadv(i,:pverp) = mf_thlforc_nadv(i,:pverp) + thlm_forcing(i,:pverp)
           mf_qtforc_nadv(i,:pverp)  = mf_qtforc_nadv(i,:pverp) + rtm_forcing(i,:pverp)

           mf_ent_nadv(i,:pverp)     = mf_ent_nadv(i,:pverp) + s_awu(i,:pverp)
           
           max_cfl_nadv(i) = MAX(max_cfl(i),max_cfl_nadv(i))
         end do

         if (t==1) then

           ztop_macmic1(:ncol,macmic_it) = mf_ztopm1(:ncol)
           ddcp_macmic1(:ncol,macmic_it) = mf_ddcp(:ncol)

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upw(:ncol,k,:clubb_mf_nup)
           end do
           
           do k=1,clubb_mf_nup
             up_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnw(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             dn_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upa(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             upa_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dna(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             dna_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upthl(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             thlu_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upqt(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             qtu_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnthl(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             thld_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnqt(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             qtd_macmic1(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = thlm_forcing(:ncol,k)
           end do
           dthl_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = rtm_forcing(:ncol,k)
           end do
           dqt_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_thlforcup(:ncol,k)
           end do
           dthlu_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_qtforcup(:ncol,k)
           end do
           dqtu_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_thlforcdn(:ncol,k)
           end do
           dthld_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_qtforcdn(:ncol,k)
           end do
           dqtd_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

         else if (t==2) then

           ztop_macmic2(:ncol,macmic_it) = mf_ztopm1(:ncol)
           ddcp_macmic2(:ncol,macmic_it) = mf_ddcp(:ncol)

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upw(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             up_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnw(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             dn_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upa(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             upa_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dna(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             dna_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upthl(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             thlu_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_upqt(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             qtu_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnthl(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             thld_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             flip(:ncol,pverp-k+1,:clubb_mf_nup) = mf_dnqt(:ncol,k,:clubb_mf_nup)
           end do

           do k=1,clubb_mf_nup
             qtd_macmic2(:ncol, 1+pverp*(clubb_mf_nup*(macmic_it-1)+k-1):pverp*(clubb_mf_nup*(macmic_it-1)+k) ) = flip(:ncol,:pverp,k)
           end do

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = thlm_forcing(:ncol,k)
           end do
           dthl_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = rtm_forcing(:ncol,k)
           end do
           dqt_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_thlforcup(:ncol,k)
           end do
           dthlu_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_qtforcup(:ncol,k)
           end do
           dqtu_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_thlforcdn(:ncol,k)
           end do
           dthld_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

           do k=1,nlev+1
             lilflip(:ncol,pverp-k+1) = mf_qtforcdn(:ncol,k)
           end do
           dqtd_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

         end if 

      end if !clubbmf

      if (t==1) then

        do k=1,nlev+1
          lilflip(:ncol,pverp-k+1) = thlm_in(:ncol,k)
        end do
        thlm_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

        do k=1,nlev+1
          lilflip(:ncol,pverp-k+1) = rtm_in(:ncol,k)
        end do
        qtm_macmic1(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

      else if (t==2) then

        do k=1,nlev+1
          lilflip(:ncol,pverp-k+1) = thlm_in(:ncol,k)
        end do
        thlm_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

        do k=1,nlev+1
          lilflip(:ncol,pverp-k+1) = rtm_in(:ncol,k)
        end do
        qtm_macmic2(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = lilflip(:ncol,:pverp)

      end if

      !  Advance CLUBB CORE one timestep in the future
      call advance_clubb_core_api( gr, pverp+1-top_lev, ncol, &
          l_implemented, dtime, fcor, sfc_elevation, hydromet_dim, &
          thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
          sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
          wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
          rtpthlp_forcing, wm_zm, wm_zt, &
          wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
          wpsclrp_sfc, wpedsclrp_sfc, &
          upwp_sfc_pert, vpwp_sfc_pert, &
          rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
          p_in_Pa, rho_zm, rho_zt, exner, &
          rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
          invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &
          rfrzm, radf, &
          wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
          grid_dx, grid_dy, &
          clubb_params, nu_vert_res_dep, lmin, &
          clubb_config_flags, &
          stats_zt(:ncol), stats_zm(:ncol), stats_sfc(:ncol), &
          um_in, vm_in, upwp_in, vpwp_in, up2_in, vp2_in, up3_in, vp3_in, &
          thlm_in, rtm_in, wprtp_in, wpthlp_in, &
          wp2_in, wp3_in, rtp2_in, rtp3_in, thlp2_in, thlp3_in, rtpthlp_in, &
          sclrm, &
          sclrp2, sclrp3, sclrprtp, sclrpthlp, &
          wpsclrp, edsclr_in, err_code, &
          rcm_inout, cloud_frac_inout, &
          wpthvp_in, wp2thvp_in, rtpthvp_in, thlpthvp_in, &
          sclrpthvp_inout, &
          wp2rtp_inout, wp2thlp_inout, uprcp_inout, &
          vprcp_inout, rc_coef_inout, &
          wp4_inout, wpup2_inout, wpvp2_inout, &
          wp2up2_inout, wp2vp2_inout, ice_supersat_frac_inout, &
          um_pert_inout, vm_pert_inout, upwp_pert_inout, vpwp_pert_inout, &
          pdf_params_chnk(lchnk), pdf_params_zm_chnk(lchnk), &
          pdf_implicit_coefs_terms_chnk(lchnk), &
          khzm_out, khzt_out, &
          qclvar_out, thlprcp_out, &
          wprcp_out, w_up_in_cloud_out, w_down_in_cloud_out,  &
          cloudy_updraft_frac_out, cloudy_downdraft_frac_out, &
          rcm_in_layer_out, cloud_cover_out, invrs_tau_zm_out )
      
      ! Note that CLUBB does not produce an error code specific to any column, and
      ! one value only for the entire chunk
      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Fatal error in CLUBB: at timestep ", get_nstep()
        write(fstderr,*) "LAT Range: ", state1%lat(1), " -- ", state1%lat(ncol)
        write(fstderr,*) "LON: Range:", state1%lon(1), " -- ", state1%lon(ncol)
        call endrun(subr//':  Fatal error in CLUBB library')
      end if

!+++ARH  ! clip negative water
      do k=1,pverp
        do i=1,ncol
          if (rtm_in(i,k) < rcm_inout(i,k)) rtm_in(i,k) = rcm_inout(i,k)
        end do
      end do
      
      if (do_rainturb) then
        
        do k=1,nlev+1
          do i=1,ncol
            rvm_in(i,k) = rtm_in(i,k) - rcm_inout(i,k) 
          end do
        end do
        
        call update_xp2_mc_api( gr, nlev+1, ncol, dtime, cloud_frac_inout, &
          rcm_inout, rvm_in, thlm_in, wm_zt, &
          exner, pre_in, pdf_params_chnk(lchnk), &
          rtp2_mc_out, thlp2_mc_out, &
          wprtp_mc_out, wpthlp_mc_out, &
          rtpthlp_mc_out)

        do k=1,nlev+1
          do i=1,ncol
            dum1 = (1._r8 - cam_in%landfrac(i))

            ! update turbulent moments based on rain evaporation  
            rtp2_in(i,k)   = rtp2_in(i,k)   + clubb_rnevap_effic * dum1 * rtp2_mc_out(i,k)   * dtime
            thlp2_in(i,k)  = thlp2_in(i,k)  + clubb_rnevap_effic * dum1 * thlp2_mc_out(i,k)  * dtime  
            wprtp_in(i,k)  = wprtp_in(i,k)  + clubb_rnevap_effic * dum1 * wprtp_mc_out(i,k)  * dtime
            wpthlp_in(i,k) = wpthlp_in(i,k) + clubb_rnevap_effic * dum1 * wpthlp_mc_out(i,k) * dtime
          end do
        end do
        
      end if  
    
      if (do_cldcool) then
        
        rcm_out_zm = zt2zm_api(pverp+1-top_lev, ncol, gr, rcm_inout )
        qrl_zm     = zt2zm_api(pverp+1-top_lev, ncol, gr, qrl_clubb )
        thlp2_rad_out(:,:) = 0._r8
        
        do i=1, ncol
          call calculate_thlp2_rad_api(nlev+1, rcm_out_zm(i,:), thlprcp_out(i,:), qrl_zm(i,:), clubb_params, &
                                       thlp2_rad_out(i,:))
        end do
        
        do i=1, ncol
          thlp2_in(i,:) = thlp2_in(i,:) + thlp2_rad_out(i,:) * dtime
          thlp2_in(i,:) = max(thl_tol**2,thlp2_in(i,:))
        end do
          
      end if
      
      !  Check to see if stats should be output, here stats are read into
      !  output arrays to make them conformable to CAM output
      if (l_stats) then
        do i=1, ncol
          call stats_end_timestep_clubb(i, stats_zt(i), stats_zm(i), stats_rad_zt(i), stats_rad_zm(i), stats_sfc(i), &
                                        out_zt, out_zm, out_radzt, out_radzm, out_sfc)
        end do
      end if

    enddo  ! end time loop

    if (do_clubb_mf) then
      ! average over nadv
      mf_L0_nadv   = mf_L0_nadv/REAL(nadv)
      mf_ztop_nadv = mf_ztop_nadv/REAL(nadv)
      mf_ztopm1_nadv = mf_ztopm1_nadv/REAL(nadv)
      mf_ddcp_nadv = mf_ddcp_nadv/REAL(nadv)
      mf_cbm1_nadv = mf_cbm1_nadv/REAL(nadv)
      mf_freq_nadv = mf_freq_nadv/REAL(nadv)

      ! accumulate in buffer
      ztopm1_macmic(:ncol) = ztopm1_macmic(:ncol) + mf_ztopm1_nadv(:ncol)
      ddcp_macmic(:ncol) = ddcp_macmic(:ncol) + mf_ddcp_nadv(:ncol)
      cbm1_macmic(:ncol) = cbm1_macmic(:ncol) + mf_cbm1_nadv(:ncol)

      if (macmic_it == cld_macmic_num_steps) then

        cbm1(:ncol) = cbm1_macmic(:ncol)/REAL(cld_macmic_num_steps)

        if (clubb_mf_up_ndt == 1) then
          ztopma(:ncol) = ztopm1_macmic(:ncol)/REAL(cld_macmic_num_steps)
        else
          ztopmn(2:clubb_mf_up_ndt,:ncol) = ztopmn(1:clubb_mf_up_ndt-1,:ncol)
          ztopmn(1,:ncol) = ztopm1_macmic(:ncol)/REAL(cld_macmic_num_steps)
          ztopma(:ncol) = 0._r8
          do t=1,clubb_mf_up_ndt
            ztopma(:ncol) = ztopma(:ncol) + ztopmn(t,:ncol)
          end do
          ztopma(:ncol) = ztopma(:ncol)/REAL(clubb_mf_up_ndt)
        end if

        if (clubb_mf_cp_ndt == 1) then
          ddcp(:ncol) = ddcp_macmic(:ncol)/REAL(cld_macmic_num_steps)
        else
          ddcpmn(2:clubb_mf_cp_ndt,:ncol) = ddcpmn(1:clubb_mf_cp_ndt-1,:ncol)
          ddcpmn(1,:ncol) = ddcp_macmic(:ncol)/REAL(cld_macmic_num_steps)
          ddcp(:ncol) = 0._r8
          do t=1,clubb_mf_cp_ndt
            ddcp(:ncol) = ddcp(:ncol) + ddcpmn(t,:ncol)
          end do
          ddcp(:ncol) = ddcp(:ncol)/REAL(clubb_mf_cp_ndt)
        end if

        ddcp(:ncol) = clubb_mf_ddalph*ddcp(:ncol)

      end if

      mf_qc(:ncol,:pverp)        = mf_qc_nadv(:ncol,:pverp)/REAL(nadv)
      mf_rcm(:ncol,:pverp)       = mf_rcm_nadv(:ncol,:pverp)/REAL(nadv)
      mf_cloudfrac(:ncol,:pverp) = mf_cloudfrac_nadv(:ncol,:pverp)/REAL(nadv)
      prec_sh(:ncol)             = mf_precc_nadv(:ncol)/REAL(nadv)
      snow_sh(:ncol)             = mf_snow_nadv(:ncol)/REAL(nadv)

      mf_thlforcup_nadv(:ncol,:pverp) = mf_thlforcup_nadv(:ncol,:pverp)/REAL(nadv)
      mf_qtforcup_nadv(:ncol,:pverp)  = mf_qtforcup_nadv(:ncol,:pverp)/REAL(nadv)
      mf_thlforcdn_nadv(:ncol,:pverp) = mf_thlforcdn_nadv(:ncol,:pverp)/REAL(nadv)
      mf_qtforcdn_nadv(:ncol,:pverp)  = mf_qtforcdn_nadv(:ncol,:pverp)/REAL(nadv)
      mf_thlforc_nadv(:ncol,:pverp)   = mf_thlforc_nadv(:ncol,:pverp)/REAL(nadv)
      mf_qtforc_nadv(:ncol,:pverp)    = mf_qtforc_nadv(:ncol,:pverp)/REAL(nadv)

      mf_ent_nadv(:ncol,:pverp) = mf_ent_nadv(:ncol,:pverp)/REAL(nadv)

    end if !clubbmf

    if (clubb_do_adv) then
      if (macmic_it  ==  cld_macmic_num_steps) then 
        
        wp2_in     = zm2zt_api( pverp+1-top_lev, ncol, gr, wp2_in )   
        wpthlp_in  = zm2zt_api( pverp+1-top_lev, ncol, gr, wpthlp_in )
        wprtp_in   = zm2zt_api( pverp+1-top_lev, ncol, gr, wprtp_in )
        up2_in     = zm2zt_api( pverp+1-top_lev, ncol, gr, up2_in )
        vp2_in     = zm2zt_api( pverp+1-top_lev, ncol, gr, vp2_in )
        thlp2_in   = zm2zt_api( pverp+1-top_lev, ncol, gr, thlp2_in )
        rtp2_in    = zm2zt_api( pverp+1-top_lev, ncol, gr, rtp2_in )
        rtpthlp_in = zm2zt_api( pverp+1-top_lev, ncol, gr, rtpthlp_in ) 

        do k=1,nlev+1
          do i=1, ncol
            thlp2_in(i,k) = max(thl_tol**2, thlp2_in(i,k))
            rtp2_in(i,k)  = max(rt_tol**2, rtp2_in(i,k))
            wp2_in(i,k)   = max(w_tol_sqd, wp2_in(i,k))
            up2_in(i,k)   = max(w_tol_sqd, up2_in(i,k))
            vp2_in(i,k)   = max(w_tol_sqd, vp2_in(i,k))
          end do
        end do
        
      end if
    end if
    
    ! Convert RTP2 and THLP2 to thermo grid for output
    rtp2_zt = zm2zt_api( pverp+1-top_lev, ncol, gr, rtp2_in )
    thl2_zt = zm2zt_api( pverp+1-top_lev, ncol, gr, thlp2_in )
    wp2_zt  = zm2zt_api( pverp+1-top_lev, ncol, gr, wp2_in )

    ! Need moist_qc and cloudfrac on thermo grid for output
    mf_qc_zt(:,:) = 0._r8
    mf_cloudfrac_zt(:,:) = 0._r8
    mf_qc_zt = zm2zt_api( pverp+1-top_lev, ncol, gr, mf_qc)
    mf_cloudfrac_zt = zm2zt_api( pverp+1-top_lev, ncol, gr, mf_cloudfrac)
    !mf_qc_zt(1:ncol,:) = zm2zt_api( pverp+1-top_lev, ncol, gr, mf_qc(1:ncol,:))
    !mf_cloudfrac_zt(1:ncol,:) = zm2zt_api( pverp+1-top_lev, ncol, gr, mf_cloudfrac(1:ncol,:))

    !  Arrays need to be "flipped" to CAM grid 
    do k=1, nlev+1
      do i=1, ncol
        um(i,pverp-k+1)           = um_in(i,k)
        vm(i,pverp-k+1)           = vm_in(i,k)
        upwp(i,pverp-k+1)         = upwp_in(i,k)
        vpwp(i,pverp-k+1)         = vpwp_in(i,k)
        wpthvp(i,pverp-k+1)       = wpthvp_in(i,k)
        wp2thvp(i,pverp-k+1)      = wp2thvp_in(i,k)
        rtpthvp(i,pverp-k+1)      = rtpthvp_in(i,k)
        thlpthvp(i,pverp-k+1)     = thlpthvp_in(i,k)
        up2(i,pverp-k+1)          = up2_in(i,k)
        vp2(i,pverp-k+1)          = vp2_in(i,k)
        up3(i,pverp-k+1)          = up3_in(i,k)
        vp3(i,pverp-k+1)          = vp3_in(i,k)
        thlm(i,pverp-k+1)         = thlm_in(i,k)
        rtm(i,pverp-k+1)          = rtm_in(i,k)
        wprtp(i,pverp-k+1)        = wprtp_in(i,k)
        wpthlp(i,pverp-k+1)       = wpthlp_in(i,k)
        wp2(i,pverp-k+1)          = wp2_in(i,k)
        wp3(i,pverp-k+1)          = wp3_in(i,k)
        rtp2(i,pverp-k+1)         = rtp2_in(i,k)
        thlp2(i,pverp-k+1)        = thlp2_in(i,k)
        rtp3(i,pverp-k+1)         = rtp3_in(i,k)
        thlp3(i,pverp-k+1)        = thlp3_in(i,k)
        rtpthlp(i,pverp-k+1)      = rtpthlp_in(i,k)
        rcm(i,pverp-k+1)          = rcm_inout(i,k)
        wprcp(i,pverp-k+1)        = wprcp_out(i,k)
        cloud_frac(i,pverp-k+1)   = min(cloud_frac_inout(i,k),1._r8)
        pdf_zm_w_1(i,pverp-k+1)   = pdf_params_zm_chnk(lchnk)%w_1(i,k)
        pdf_zm_w_2(i,pverp-k+1)   = pdf_params_zm_chnk(lchnk)%w_2(i,k)
        pdf_zm_varnce_w_1(i,pverp-k+1) = pdf_params_zm_chnk(lchnk)%varnce_w_1(i,k)
        pdf_zm_varnce_w_2(i,pverp-k+1) = pdf_params_zm_chnk(lchnk)%varnce_w_2(i,k)
        pdf_zm_mixt_frac(i,pverp-k+1)  = pdf_params_zm_chnk(lchnk)%mixt_frac(i,k)
        rcm_in_layer(i,pverp-k+1) = rcm_in_layer_out(i,k)
        cloud_cover(i,pverp-k+1)  = min(cloud_cover_out(i,k),1._r8)
        zt_out(i,pverp-k+1)       = zt_g(i,k)
        zi_out(i,pverp-k+1)       = zi_g(i,k)
        khzm(i,pverp-k+1)         = khzm_out(i,k)
        qclvar(i,pverp-k+1)       = min(1._r8,qclvar_out(i,k))
        wm_zt_out(i,pverp-k+1)    = wm_zt(i,k)
        wp2rtp(i,pverp-k+1)       = wp2rtp_inout(i,k)
        wp2thlp(i,pverp-k+1)      = wp2thlp_inout(i,k)
        uprcp(i,pverp-k+1)        = uprcp_inout(i,k)
        vprcp(i,pverp-k+1)        = vprcp_inout(i,k)
        rc_coef(i,pverp-k+1)      = rc_coef_inout(i,k)
        wp4(i,pverp-k+1)          = wp4_inout(i,k)
        wpup2(i,pverp-k+1)        = wpup2_inout(i,k)
        wpvp2(i,pverp-k+1)        = wpvp2_inout(i,k)
        wp2up2(i,pverp-k+1)       = wp2up2_inout(i,k)
        wp2vp2(i,pverp-k+1)       = wp2vp2_inout(i,k)
        ice_supersat_frac(i,pverp-k+1) = ice_supersat_frac_inout(i,k)

        rtp2_zt_out(i,pverp-k+1)  = rtp2_zt(i,k)
        thl2_zt_out(i,pverp-k+1)  = thl2_zt(i,k)
        wp2_zt_out(i,pverp-k+1)   = wp2_zt(i,k)
        
      end do
    end do

    do k=1, nlev+1
      do i=1, ncol
        
        mean_rt = pdf_params_chnk(lchnk)%mixt_frac(i,k) &
                  * pdf_params_chnk(lchnk)%rt_1(i,k) &
                  + ( 1.0_r8 - pdf_params_chnk(lchnk)%mixt_frac(i,k) ) &
                    * pdf_params_chnk(lchnk)%rt_2(i,k)
        
        pdfp_rtp2(i,pverp-k+1) = pdf_params_chnk(lchnk)%mixt_frac(i,k) &
                                 * ( ( pdf_params_chnk(lchnk)%rt_1(i,k) - mean_rt )**2 &
                                     + pdf_params_chnk(lchnk)%varnce_rt_1(i,k) ) &
                                 + ( 1.0_r8 - pdf_params_chnk(lchnk)%mixt_frac(i,k) ) &
                                   * ( ( pdf_params_chnk(lchnk)%rt_2(i,k) - mean_rt )**2 &
                                       + pdf_params_chnk(lchnk)%varnce_rt_2(i,k) )
      end do
    end do

    do ixind=1,edsclr_dim
      do k=1, nlev+1
        do i=1, ncol
          edsclr_out(i,pverp-k+1,ixind) = edsclr_in(i,k,ixind)
        end do
      end do
    end do

    if (do_clubb_mf) then
      do k=1, nlev+1
        do i=1, ncol
          mf_dry_a_output(i,pverp-k+1)     = mf_dry_a(i,k)
          mf_moist_a_output(i,pverp-k+1)   = mf_moist_a(i,k)
          mf_dry_w_output(i,pverp-k+1)     = mf_dry_w(i,k)
          mf_moist_w_output(i,pverp-k+1)   = mf_moist_w(i,k)
          mf_dry_qt_output(i,pverp-k+1)    = mf_dry_qt(i,k)
          mf_moist_qt_output(i,pverp-k+1)  = mf_moist_qt(i,k)
          mf_dry_thl_output(i,pverp-k+1)   = mf_dry_thl(i,k)
          mf_moist_thl_output(i,pverp-k+1) = mf_moist_thl(i,k)
          mf_dry_u_output(i,pverp-k+1)     = mf_dry_u(i,k)
          mf_moist_u_output(i,pverp-k+1)   = mf_moist_u(i,k)
          mf_dry_v_output(i,pverp-k+1)     = mf_dry_v(i,k)
          mf_moist_v_output(i,pverp-k+1)   = mf_moist_v(i,k)
          mf_moist_qc_output(i,pverp-k+1)  = mf_moist_qc(i,k)

          s_ae_output(i,pverp-k+1)         = s_ae(i,k)
          s_aw_output(i,pverp-k+1)         = s_aw(i,k)

          s_awthlup_output(i,pverp-k+1)    = s_awthlup(i,k)
          s_awqtup_output(i,pverp-k+1)     = s_awqtup(i,k)
          s_awthldn_output(i,pverp-k+1)    = s_awthldn(i,k)
          s_awqtdn_output(i,pverp-k+1)     = s_awqtdn(i,k)
          s_awthl_output(i,pverp-k+1)      = s_awthl(i,k)
          s_awqt_output(i,pverp-k+1)       = s_awqt(i,k)

          s_awuup_output(i,pverp-k+1)      = s_awuup(i,k)
          s_awvup_output(i,pverp-k+1)      = s_awvup(i,k)
          s_awudn_output(i,pverp-k+1)      = s_awudn(i,k)
          s_awvdn_output(i,pverp-k+1)      = s_awvdn(i,k)

          s_awu_output(i,pverp-k+1)        = s_awu(i,k)
          s_awv_output(i,pverp-k+1)        = s_awv(i,k)
          s_aww_output(i,pverp-k+1)        = s_aww(i,k)

          mf_thlflxup_output(i,pverp-k+1)  = mf_thlflxup(i,k)
          mf_qtflxup_output(i,pverp-k+1)   = mf_qtflxup(i,k)
          mf_thlflxdn_output(i,pverp-k+1)  = mf_thlflxdn(i,k)
          mf_qtflxdn_output(i,pverp-k+1)   = mf_qtflxdn(i,k)
          mf_thlflx_output(i,pverp-k+1)    = mf_thlflx(i,k)
          mf_qtflx_output(i,pverp-k+1)     = mf_qtflx(i,k)
          mf_thvflx_output(i,pverp-k+1)    = mf_thvflx(i,k)

          mf_uflxup_output(i,pverp-k+1)    = mf_uflxup(i,k)
          mf_vflxup_output(i,pverp-k+1)    = mf_vflxup(i,k)
          mf_uflxdn_output(i,pverp-k+1)    = mf_uflxdn(i,k)
          mf_vflxdn_output(i,pverp-k+1)    = mf_vflxdn(i,k)
          mf_uflx_output(i,pverp-k+1)      = mf_uflx(i,k)
          mf_vflx_output(i,pverp-k+1)      = mf_vflx(i,k)

          mf_rcm_output(i,pverp-k+1)       = mf_rcm(i,k)
          mf_precc_output(i,pverp-k+1)     = mf_precc(i,k)

          if (k.ne.1) then
            mf_thlforcup_output(i,pverp-k+1)          = mf_thlforcup_nadv(i,k)
            mf_qtforcup_output(i,pverp-k+1)           = mf_qtforcup_nadv(i,k)
            mf_thlforcdn_output(i,pverp-k+1)          = mf_thlforcdn_nadv(i,k)
            mf_qtforcdn_output(i,pverp-k+1)           = mf_qtforcdn_nadv(i,k)
            mf_thlforc_output(i,pverp-k+1)            = mf_thlforc_nadv(i,k)
            mf_qtforc_output(i,pverp-k+1)             = mf_qtforc_nadv(i,k)
            mf_sqtup_output(i,pverp-k+1)              = mf_sqtup(i,k)
            mf_sqtdn_output(i,pverp-k+1)              = mf_sqtdn(i,k)

            mf_cloudfrac_output(i,pverp-k+1)          = mf_cloudfrac_zt(i,k)
            mf_ent_output(i,pverp-k+1)                = mf_ent_nadv(i,k)
            mf_qc_output(i,pverp-k+1)                 = mf_qc_zt(i,k)
          end if

          mf_upa_flip(i,pverp-k+1,:clubb_mf_nup)       = mf_upa(i,k,:clubb_mf_nup)
          mf_upw_flip(i,pverp-k+1,:clubb_mf_nup)       = mf_upw(i,k,:clubb_mf_nup)
          mf_upmf_flip(i,pverp-k+1,:clubb_mf_nup)      = mf_upmf(i,k,:clubb_mf_nup)
          mf_upqt_flip(i,pverp-k+1,:clubb_mf_nup)      = mf_upqt(i,k,:clubb_mf_nup)
          mf_upthl_flip(i,pverp-k+1,:clubb_mf_nup)     = mf_upthl(i,k,:clubb_mf_nup)
          mf_upthv_flip(i,pverp-k+1,:clubb_mf_nup)     = mf_upthv(i,k,:clubb_mf_nup)
          mf_upth_flip(i,pverp-k+1,:clubb_mf_nup)      = mf_upth(i,k,:clubb_mf_nup)
          mf_upqc_flip(i,pverp-k+1,:clubb_mf_nup)      = mf_upqc(i,k,:clubb_mf_nup)
          mf_upent_flip(i,pverp-k+1,:clubb_mf_nup)     = mf_upent(i,k,:clubb_mf_nup)
          mf_updet_flip(i,pverp-k+1,:clubb_mf_nup)     = mf_updet(i,k,:clubb_mf_nup)
          mf_upbuoy_flip(i,pverp-k+1,:clubb_mf_nup)    = mf_upbuoy(i,k,:clubb_mf_nup)
          mf_dnw_flip(i,pverp-k+1,:clubb_mf_nup)       = mf_dnw(i,k,:clubb_mf_nup)
          mf_dnthl_flip(i,pverp-k+1,:clubb_mf_nup)     = mf_dnthl(i,k,:clubb_mf_nup)
          mf_dnqt_flip(i,pverp-k+1,:clubb_mf_nup)      = mf_dnqt(i,k,:clubb_mf_nup)
        end do
      end do
    end if !clubbmf

    if (do_clubb_mf) then
      ! these fillvalues won't average correctly
      !if (mf_ztop_nadv == 0._r8) mf_ztop_nadv = fillvalue
      !if (mf_L0_nadv == 0._r8) mf_L0_nadv = fillvalue

      !mf_ztop_output(i) = ztopma(i) !mf_ztop_nadv
      mf_ztop_output(:ncol) = mf_ztop_nadv(:ncol)

      mf_L0_output(:ncol)   = mf_L0_nadv(:ncol)
      mf_cfl_output(:ncol)  = max_cfl_nadv(:ncol)

      !mf_ddcp_output(i) = ddcp(i) !mf_ddcp_nadv !ddcp(i)
      mf_ddcp_output(:ncol) = mf_ddcp_nadv(:ncol)
      mf_freq_output(:ncol) = mf_freq_nadv(:ncol)

      do k=1,clubb_mf_nup
        mf_upa_output(:ncol,pverp*(k-1)+1:pverp*k)   = mf_upa_flip(:ncol,:pverp,k)
        mf_upw_output(:ncol,pverp*(k-1)+1:pverp*k)   = mf_upw_flip(:ncol,:pverp,k)
        mf_upmf_output(:ncol,pverp*(k-1)+1:pverp*k)  = mf_upmf_flip(:ncol,:pverp,k)
        mf_upqt_output(:ncol,pverp*(k-1)+1:pverp*k)  = mf_upqt_flip(:ncol,:pverp,k)
        mf_upthl_output(:ncol,pverp*(k-1)+1:pverp*k) = mf_upthl_flip(:ncol,:pverp,k)
        mf_upthv_output(:ncol,pverp*(k-1)+1:pverp*k) = mf_upthv_flip(:ncol,:pverp,k)
        mf_upth_output(:ncol,pverp*(k-1)+1:pverp*k)  = mf_upth_flip(:ncol,:pverp,k)
        mf_upqc_output(:ncol,pverp*(k-1)+1:pverp*k)  = mf_upqc_flip(:ncol,:pverp,k)
        mf_upent_output(:ncol,pverp*(k-1)+1:pverp*k) = mf_upent_flip(:ncol,:pverp,k)
        mf_updet_output(:ncol,pverp*(k-1)+1:pverp*k) = mf_updet_flip(:ncol,:pverp,k)
        mf_upbuoy_output(:ncol,pverp*(k-1)+1:pverp*k)= mf_upbuoy_flip(:ncol,:pverp,k)
        mf_dnw_output(:ncol,pverp*(k-1)+1:pverp*k)   = mf_dnw_flip(:ncol,:pverp,k)
        mf_dnthl_output(:ncol,pverp*(k-1)+1:pverp*k) = mf_dnthl_flip(:ncol,:pverp,k)
        mf_dnqt_output(:ncol,pverp*(k-1)+1:pverp*k)  = mf_dnqt_flip(:ncol,:pverp,k)
      end do
    end if !clubbmf

      ! Values to use above top_lev, for variables that have not already been
      ! set up there. These are mostly fill values that should not actually be
      ! used in the run, but may end up in diagnostic output.
    do k=1, top_lev-1
      do i=1, ncol
        upwp(i,k)         = 0._r8
        vpwp(i,k)         = 0._r8
        rcm(i,k)          = 0._r8
        wprcp(i,k)        = 0._r8
        cloud_frac(i,k)   = 0._r8
        rcm_in_layer(i,k) = 0._r8
        zt_out(i,k)       = 0._r8
        zi_out(i,k)       = 0._r8
        khzm(i,k)         = 0._r8
        qclvar(i,k)       = 2._r8
      end do
    end do 

    ! enforce zero tracer tendencies above the top_lev level -- no change
    icnt=0
    do ixind=1,pcnst
      if (lq(ixind)) then 
        icnt=icnt+1
        
        do i=1, ncol
          edsclr_out(i,:top_lev-1,icnt) = state1%q(i,:top_lev-1,ixind)
        end do
          
      end if
    end do

    !  Fill up arrays needed for McICA.  Note we do not want the ghost point,
    !   thus why the second loop is needed.
    zi_out(:,1) = 0._r8

    ! Compute static energy using CLUBB's variables
    do k=1,pver
      do i=1, ncol
        clubb_s(i,k) = cpairv(i,k,lchnk) * thlm(i,k) / inv_exner_clubb(i,k) &
                       + latvap * rcm(i,k) &
                       + gravit * state1%zm(i,k) + state1%phis(i)
      end do
    end do

    ! Section below is concentrated on energy fixing for conservation.
    !   because CLUBB and CAM's thermodynamic variables are different.  

    ! Initialize clubbtop to top_lev, for finding the highlest level CLUBB is
    !  active for informing where to apply the energy fixer.
    do i=1, ncol
      clubbtop(i) = top_lev
      do while ((rtp2(i,clubbtop(i)) <= 1.e-15_r8 .and. rcm(i,clubbtop(i))  ==  0._r8) .and. clubbtop(i) <  pver)
        clubbtop(i) = clubbtop(i) + 1
      end do    
    end do
    !
    ! set pbuf field so that HB scheme is only applied above CLUBB top
    !
    if (do_hb_above_clubb) then
      call pbuf_set_field(pbuf, clubbtop_idx, clubbtop)
    endif


    ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
    ! after CLUBB is called.  This is for energy conservation purposes.
    se_a(:) = 0._r8
    ke_a(:) = 0._r8
    wv_a(:) = 0._r8
    wl_a(:) = 0._r8
       
    do k=1,pver
      do i=1, ncol
        se_a(i) = se_a(i) + clubb_s(i,k)*state1%pdel(i,k)*rga
        ke_a(i) = ke_a(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*rga
        wv_a(i) = wv_a(i) + (rtm(i,k)-rcm(i,k))*state1%pdeldry(i,k)*rga
        wl_a(i) = wl_a(i) + (rcm(i,k))*state1%pdeldry(i,k)*rga
      end do
    end do   
    
    ! Do the same as above, but for before CLUBB was called.
    se_b(:) = 0._r8
    ke_b(:) = 0._r8
    wv_b(:) = 0._r8
    wl_b(:) = 0._r8   
    
    do k=1, pver
      do i=1, ncol
        se_b(i) = se_b(i) + state1%s(i,k)*state1%pdel(i,k)*rga
        ke_b(i) = ke_b(i) + 0.5_r8*(state1%u(i,k)**2+state1%v(i,k)**2)*state1%pdel(i,k)*rga
        wv_b(i) = wv_b(i) + state1%q(i,k,ixq)*state1%pdeldry(i,k)*rga
        wl_b(i) = wl_b(i) + state1%q(i,k,ixcldliq)*state1%pdeldry(i,k)*rga
      end do
    end do
      
     
    do i=1, ncol
      ! Based on these integrals, compute the total energy before and after CLUBB call
      te_a(i) = se_a(i) + ke_a(i) + (latvap+latice) * wv_a(i) + latice * wl_a(i)
      te_b(i) = se_b(i) + ke_b(i) + (latvap+latice) * wv_b(i) + latice * wl_b(i)
      
      ! Take into account the surface fluxes of heat and moisture
      !  Use correct qflux from cam_in, not lhf/latvap as was done previously
      te_b(i) = te_b(i) + (cam_in%shf(i)+cam_in%cflx(i,1)*(latvap+latice)) * hdtime      

      ! subtract enthalpy of falling precip from tb
      te_b(i) = te_b(i) - prec_sh(i)*1000._r8*latice*hdtime

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis(i) = (te_a(i) - te_b(i))/(state1%pint(i,pverp)-state1%pint(i,clubbtop(i)))
    end do


    ! Fix the total energy coming out of CLUBB so it achieves energy conservation.
    ! Apply this fixer throughout the column evenly, but only at layers where 
    ! CLUBB is active.
    !
    ! NOTE: The energy fixer seems to cause the climate to change significantly
    ! when using specified dynamics, so allow this to be turned off via a namelist
    ! variable.
    if (clubb_do_energyfix) then
      do i=1, ncol
        do k=clubbtop(i),pver
          clubb_s(i,k) = clubb_s(i,k) - se_dis(i)*gravit
        end do
        ! convert to units of +ve [K]
        se_dis(i) = -1._r8*se_dis(i)*gravit/cpairv(i,pver,lchnk)
      end do
    endif

      
    !  Now compute the tendencies of CLUBB to CAM, note that pverp is the ghost point
    !  for all variables and therefore is never called in this loop
    do k=1, pver
      do i=1, ncol

        ptend_loc%u(i,k)          = (um(i,k) - state1%u(i,k))               / hdtime ! east-west wind
        ptend_loc%v(i,k)          = (vm(i,k) - state1%v(i,k))               / hdtime ! north-south wind
        ptend_loc%q(i,k,ixq)      = (rtm(i,k) - rcm(i,k)-state1%q(i,k,ixq)) / hdtime ! water vapor
        ptend_loc%q(i,k,ixcldliq) = (rcm(i,k) - state1%q(i,k,ixcldliq))     / hdtime ! Tendency of liquid water
        ptend_loc%s(i,k)          = (clubb_s(i,k) - state1%s(i,k))          / hdtime ! Tendency of static energy

      end do
    end do
    
    
    if (clubb_do_adv) then
      if (macmic_it == cld_macmic_num_steps) then
        
        do k=1, pver
          do i=1, ncol

            ! Here add a constant to moments which can be either positive or 
            !  negative.  This is to prevent clipping when dynamics tries to
            !  make all constituents positive 
            wp3(i,k)     = wp3(i,k)     + wp3_const
            rtpthlp(i,k) = rtpthlp(i,k) + rtpthlp_const
            wpthlp(i,k)  = wpthlp(i,k)  + wpthlp_const
            wprtp(i,k)   = wprtp(i,k)   + wprtp_const

            ptend_loc%q(i,k,ixthlp2)   = (thlp2(i,k)   - state1%q(i,k,ixthlp2))   / hdtime ! THLP Variance
            ptend_loc%q(i,k,ixrtp2)    = (rtp2(i,k)    - state1%q(i,k,ixrtp2))    / hdtime ! RTP Variance
            ptend_loc%q(i,k,ixrtpthlp) = (rtpthlp(i,k) - state1%q(i,k,ixrtpthlp)) / hdtime ! RTP THLP covariance
            ptend_loc%q(i,k,ixwpthlp)  = (wpthlp(i,k)  - state1%q(i,k,ixwpthlp))  / hdtime ! WPTHLP 
            ptend_loc%q(i,k,ixwprtp)   = (wprtp(i,k)   - state1%q(i,k,ixwprtp))   / hdtime ! WPRTP
            ptend_loc%q(i,k,ixwp2)     = (wp2(i,k)     - state1%q(i,k,ixwp2))     / hdtime ! WP2
            ptend_loc%q(i,k,ixwp3)     = (wp3(i,k)     - state1%q(i,k,ixwp3))     / hdtime ! WP3
            ptend_loc%q(i,k,ixup2)     = (up2(i,k)     - state1%q(i,k,ixup2))     / hdtime ! UP2
            ptend_loc%q(i,k,ixvp2)     = (vp2(i,k)     - state1%q(i,k,ixvp2))     / hdtime ! VP2
            
          end do
        end do
        
      else
        
        do k=1, pver
          do i=1, ncol
            ptend_loc%q(i,k,ixthlp2)   = 0._r8
            ptend_loc%q(i,k,ixrtp2)    = 0._r8
            ptend_loc%q(i,k,ixrtpthlp) = 0._r8
            ptend_loc%q(i,k,ixwpthlp)  = 0._r8
            ptend_loc%q(i,k,ixwprtp)   = 0._r8
            ptend_loc%q(i,k,ixwp2)     = 0._r8
            ptend_loc%q(i,k,ixwp3)     = 0._r8
            ptend_loc%q(i,k,ixup2)     = 0._r8
            ptend_loc%q(i,k,ixvp2)     = 0._r8 
          end do
        end do
        
      end if
    end if

    !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
    !  Loading up this array doesn't mean the tendencies are applied.  
    ! edsclr_out is compressed with just the constituents being used, ptend and state are not compressed
    icnt=0
    do ixind=1,pcnst
      if (lq(ixind)) then
        icnt=icnt+1
        if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
            (ixind /= ixthlp2)   .and. (ixind /= ixrtp2)   .and.&
            (ixind /= ixrtpthlp) .and. (ixind /= ixwpthlp) .and.&
            (ixind /= ixwprtp)   .and. (ixind /= ixwp2)    .and.&
            (ixind /= ixwp3)     .and. (ixind /= ixup2)    .and. (ixind /= ixvp2) ) then
            
          do k=1, pver
            do i=1, ncol
              ptend_loc%q(i,k,ixind) = (edsclr_out(i,k,icnt)-state1%q(i,k,ixind))/hdtime ! transported constituents 
            end do
          end do
              
        end if
      end if
    end do
  
    call t_stopf("clubb_tend_cam_i_loop")

    call outfld('KVH_CLUBB', khzm, pcols, lchnk)

    ! dte / hdtime = [kg/s2]/[s] = W/m2
    eleak(:ncol) = (te_a(:ncol) - te_b(:ncol))/hdtime
    call outfld('ELEAK_CLUBB', eleak, pcols, lchnk)
    call outfld('TFIX_CLUBB', se_dis, pcols, lchnk)

    ! Add constant to ghost point so that output is not corrupted 
    if (clubb_do_adv) then
      if (macmic_it == cld_macmic_num_steps) then
        wp3(:,pverp)     = wp3(:,pverp)     + wp3_const
        rtpthlp(:,pverp) = rtpthlp(:,pverp) + rtpthlp_const
        wpthlp(:,pverp)  = wpthlp(:,pverp)  + wpthlp_const
        wprtp(:,pverp)   = wprtp(:,pverp)   + wprtp_const
      end if
    end if   

    ! ------------------------------------------------- !
    ! End column computation of CLUBB, begin to apply   !
    ! and compute output, etc                           !
    ! ------------------------------------------------- !

    !  Output CLUBB tendencies (convert dry basis to wet for consistency with  history variable definition)
    temp2d(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    call outfld( 'RVMTEND_CLUBB', temp2d, pcols, lchnk)

    temp2d(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    call outfld( 'RCMTEND_CLUBB', temp2d, pcols, lchnk)

    temp2d(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldice)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    call outfld( 'RIMTEND_CLUBB', temp2d, pcols, lchnk)

    call outfld( 'STEND_CLUBB', ptend_loc%s,pcols, lchnk)
    call outfld( 'UTEND_CLUBB', ptend_loc%u,pcols, lchnk)
    call outfld( 'VTEND_CLUBB', ptend_loc%v,pcols, lchnk)

    cmeliq(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    call outfld( 'CMELIQ', cmeliq, pcols, lchnk)

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state1,ptend_loc,hdtime)
 
    ! Due to the order of operation of CLUBB, which closes on liquid first, 
    ! then advances it's predictive equations second, this can lead to 
    ! RHliq > 1 directly before microphysics is called.  Therefore, we use 
    ! ice_macro_tend to enforce RHliq <= 1 everywhere before microphysics is called. 
 
    if (clubb_do_liqsupersat) then
 
      ! -------------------------------------- !
      ! Ice Saturation Adjustment Computation  !
      ! -------------------------------------- !
 
      latsub = latvap + latice

      lq2(:)        = .FALSE.
      lq2(ixq)      = .TRUE.
      lq2(ixcldliq) = .TRUE.
      lq2(ixnumliq) = .TRUE.
 
      call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq2 )
 
      stend(:ncol,:)=0._r8
      qvtend(:ncol,:)=0._r8
      qctend(:ncol,:)=0._r8
      inctend(:ncol,:)=0._r8
 
      call liquid_macro_tend(npccn(1:ncol,top_lev:pver), state1%t(1:ncol,top_lev:pver),                      &
                             state1%pmid(1:ncol,top_lev:pver), state1%q(1:ncol,top_lev:pver,ixq),            &
                             state1%q(1:ncol,top_lev:pver,ixcldliq), state1%q(1:ncol,top_lev:pver,ixnumliq), &
                             latvap, hdtime, stend(1:ncol,top_lev:pver),qvtend(1:ncol,top_lev:pver),         &
                             qctend(1:ncol,top_lev:pver), inctend(1:ncol,top_lev:pver), ncol*(pver-top_lev+1))

      ! update local copy of state with the tendencies
      ptend_loc%q(:ncol,top_lev:pver,ixq)=qvtend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixcldliq)=qctend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixnumliq)=inctend(:ncol,top_lev:pver)
      ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver)
 
      ! Add the ice tendency to the output tendency
      call physics_ptend_sum(ptend_loc, ptend_all, ncol)
 
      ! ptend_loc is reset to zero by this call
      call physics_update(state1, ptend_loc, hdtime)
 
      ! Write output for tendencies:
      !        oufld: QVTENDICE,QCTENDICE,NCTENDICE,FQTENDICE
      temp2d(:ncol,:pver) =  stend(:ncol,:pver)/cpairv(:ncol,:pver,lchnk)
      call outfld( 'TTENDICE', temp2d, pcols, lchnk )
      call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
      call outfld( 'QCTENDICE', qctend, pcols, lchnk )
      call outfld( 'NCTENDICE', inctend, pcols, lchnk )
 
      where(qctend .ne. 0._r8)
        fqtend = 1._r8
      elsewhere
        fqtend = 0._r8
      end where
 
      call outfld( 'FQTENDICE', fqtend, pcols, lchnk )
    end if
 
    ! ------------------------------------------------------------ !
    ! The rest of the code deals with diagnosing variables         !
    ! for microphysics/radiation computation and macrophysics      !
    ! ------------------------------------------------------------ !
     
    ! --------------------------------------------------------------------------------- !  
    !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
    !  Detrainment of convective condensate into the environment or stratiform cloud    !
    ! --------------------------------------------------------------------------------- ! 

    !  Initialize the shallow convective detrainment rate, will always be zero
    dlf2(:,:) = 0.0_r8
    dlf_liq_out(:,:) = 0.0_r8
    dlf_ice_out(:,:) = 0.0_r8

    lqice(:)        = .false.
    lqice(ixcldliq) = .true.
    lqice(ixcldice) = .true.
    lqice(ixnumliq) = .true.
    lqice(ixnumice) = .true.

    call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lq=lqice)

    if (zmconv_microp) then
      call pbuf_get_field(pbuf, dlfzm_idx, dlfzm)
      call pbuf_get_field(pbuf, difzm_idx, difzm)
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlfzm)
      call pbuf_get_field(pbuf, dnifzm_idx, dnifzm)
    end if
   
    do k=1,pver
      do i=1,ncol
        if( state1%t(i,k) > meltpt_temp ) then
          dum1 = 0.0_r8
        elseif ( state1%t(i,k) < dt_low ) then
          dum1 = 1.0_r8
        else
          dum1 = ( meltpt_temp - state1%t(i,k) ) / ( meltpt_temp - dt_low ) 
        endif

        if (zmconv_microp) then
          ptend_loc%q(i,k,ixcldliq) = dlfzm(i,k) + dlf2(i,k) * ( 1._r8 - dum1 )
          ptend_loc%q(i,k,ixcldice) = difzm(i,k) + dlf2(i,k) * dum1

          ptend_loc%q(i,k,ixnumliq) = dnlfzm(i,k) + 3._r8 * ( dlf2(i,k) * ( 1._r8 - dum1 ) )   &
                                                 / (4._r8*3.14_r8*dl_rad**3*997._r8)      ! Shallow Convection
          ptend_loc%q(i,k,ixnumice) = dnifzm(i,k) + 3._r8 * ( dlf2(i,k) * dum1 ) &
                                                 / (4._r8*3.14_r8*di_rad**3*500._r8)      ! Shallow Convection
          ptend_loc%s(i,k)          = dlf2(i,k) * dum1 * latice
        else       

          ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
          ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
          ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
                                   / (4._r8*3.14_r8*dl_rad**3*997._r8) + & ! Deep    Convection
                                   3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
                                   / (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection 
          ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
                                   / (4._r8*3.14_r8*di_rad**3*500._r8) + & ! Deep    Convection
                                   3._r8 * (                         dlf2(i,k)    *  dum1 ) &
                                   / (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
          ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

          dlf_liq_out(i,k) = dlf(i,k) * ( 1._r8 - dum1 ) 
          dlf_ice_out(i,k) = dlf(i,k) * dum1
        end if
        ! convert moist dlf tendencies to dry
        ptend_loc%q(i,k,ixcldliq) = ptend_loc%q(i,k,ixcldliq)*state1%pdel(i,k)/state1%pdeldry(i,k)
        ptend_loc%q(i,k,ixcldice) = ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)/state1%pdeldry(i,k)

        ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
        !   track of the integrals of ice and static energy that is effected from conversion to ice
        !   so that the energy checker doesn't complain.
        det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)*rga
        det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdeldry(i,k)*rga
      enddo
    enddo
   
    det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water

    ! output moist basis to be consistent with history variable definition
    temp2d(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)    
    call outfld( 'DPDLFLIQ', temp2d, pcols, lchnk)

    ! output moist basis to be consistent with history variable definition
    temp2d(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldice)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    call outfld( 'DPDLFICE', temp2d, pcols, lchnk)
   
    temp2d(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpairv(:ncol,:pver, lchnk)
    call outfld( 'DPDLFT',   temp2d, pcols, lchnk)
  
    call outfld( 'DETNLIQTND', ptend_loc%q(:,:,ixnumliq),pcols, lchnk )

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state1,ptend_loc,hdtime)

    ! ptend_all now has all accumulated tendencies.  Convert the tendencies for the
    ! wet constituents to wet air basis.
    do ixind = 1, pcnst
      if (lq(ixind) .and. cnst_type(ixind) == 'wet') then
        do k = 1, pver
          do i = 1, ncol
            ptend_all%q(i,k,ixind) = ptend_all%q(i,k,ixind)*state1%pdeldry(i,k)/state1%pdel(i,k)
          end do
        end do
      end if
    end do

    ! ------------------------------------------------- !
    ! Diagnose relative cloud water variance            !
    ! ------------------------------------------------- !

    if (deep_scheme  ==  'CLUBB_SGS') then
      relvarmax = 2.0_r8
    else
      relvarmax = 10.0_r8
    endif
   
    relvar(:,:) = relvarmax  ! default

    if (deep_scheme .ne. 'CLUBB_SGS') then    
      where (rcm(:ncol,:pver) /= 0 .and. qclvar(:ncol,:pver) /= 0) &
          relvar(:ncol,:pver) = min(relvarmax,max(0.001_r8,rcm(:ncol,:pver)**2/qclvar(:ncol,:pver)))
    endif
   
    ! ------------------------------------------------- !
    ! Optional Accretion enhancement factor             !
    ! ------------------------------------------------- !   
     accre_enhan(:ncol,:pver) = 1._r8
   
    ! ------------------------------------------------- !
    ! Diagnose some output variables                    !
    ! ------------------------------------------------- !

    !  density
    rho(1:ncol,1:pver) = rga*state1%pdel(1:ncol,1:pver)/(state1%zi(1:ncol,1:pver)-state1%zi(1:ncol,2:pverp))
    rho(1:ncol,pverp) = rho(1:ncol,pver)

    wpthvp_diag(:,:) = 0.0_r8
    do k=1,pver
      do i=1,ncol
        eps = rairv(i,k,lchnk)*inv_rh2o
        !  buoyancy flux
        wpthvp_diag(i,k) = (wpthlp(i,k)-(apply_const*wpthlp_const))+((1._r8-eps)/eps)*theta0* &
                       (wprtp(i,k)-(apply_const*wprtp_const))+((latvap/cpairv(i,k,lchnk))* &
                       state1%exner(i,k)-(1._r8/eps)*theta0)*wprcp(i,k)

        !  total water mixing ratio
        qt_output(i,k) = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)+state1%q(i,k,ixcldice)
        !  liquid water potential temperature
        thetal_output(i,k) = (state1%t(i,k)*state1%exner(i,k))-(latvap/cpairv(i,k,lchnk))*state1%q(i,k,ixcldliq)
        !  liquid water static energy
        sl_output(i,k) = cpairv(i,k,lchnk)*state1%t(i,k)+gravit*state1%zm(i,k)-latvap*state1%q(i,k,ixcldliq)
      enddo
    enddo
   
    do k=1,pverp
      do i=1,ncol
         !wpthlp_output(i,k)  = (wpthlp(i,k)-(apply_const*wpthlp_const))*rho(i,k)*cpair !  liquid water potential temperature flux
         !wprtp_output(i,k)   = (wprtp(i,k)-(apply_const*wprtp_const))*rho(i,k)*latvap  !  total water mixig ratio flux
         wpthlp_output(i,k)  = wpthlp(i,k)
         wprtp_output(i,k)   = wprtp(i,k)
         rtpthlp_output(i,k) = rtpthlp(i,k)-(apply_const*rtpthlp_const)                !  rtpthlp output
         wp3_output(i,k)     = wp3(i,k) - (apply_const*wp3_const)                      !  wp3 output
         tke(i,k)            = 0.5_r8*(up2(i,k)+vp2(i,k)+wp2(i,k))                     !  turbulent kinetic energy
         if (do_clubb_mf) then
           ! comment out for kinemtatic fluxes 
           !mf_thlflx_output(i,k) = mf_thlflx_output(i,k)*rho(i,k)*cpair
           !mf_thvflx_output(i,k) = mf_thvflx_output(i,k)*rho(i,k)*cpair
           !mf_qtflx_output(i,k)  = mf_qtflx_output(i,k)*rho(i,k)*latvap

           mf_precc_output(i,k)  = mf_precc_output(i,k)/rhoh2o
         end if
      enddo
    enddo
   
   rcm_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = rcm(:ncol,:pverp)
   cldfrac_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = cloud_frac(:ncol,:pverp)
   wpthlp_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = wpthlp_output(:ncol,:pverp)
   wprtp_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = wprtp_output(:ncol,:pverp)
   wpthvp_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = wpthvp(:ncol,:pverp)
   if (do_clubb_mf) then
     mf_thlflx_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = mf_thlflx_output(:ncol,:pverp)
     mf_qtflx_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = mf_qtflx_output(:ncol,:pverp)
     mf_thvflx_macmic(:ncol,pverp*(macmic_it-1)+1:pverp*macmic_it) = mf_thvflx_output(:ncol,:pverp)
   end if

    ! --------------------------------------------------------------------------------- ! 
    !  Diagnose some quantities that are computed in macrop_tend here.                  !
    !  These are inputs required for the microphysics calculation.                      !
    !                                                                                   !
    !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
    ! --------------------------------------------------------------------------------- ! 
 
    !  initialize variables 
    alst(:,:) = 0.0_r8
    qlst(:,:) = 0.0_r8 
 
    do k=1,pver
      do i=1,ncol
        alst(i,k) = cloud_frac(i,k)   
        qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
      enddo
    enddo
 
    ! --------------------------------------------------------------------------------- !  
    !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
    ! --------------------------------------------------------------------------------- ! 
 
    deepcu(:,:) = 0.0_r8
    shalcu(:,:) = 0.0_r8
    sh_icwmr(:,:) = 0.0_r8

    do k=1,pver-1
      do i=1,ncol
         !  diagnose the deep convective cloud fraction, as done in macrophysics based on the 
         !  deep convective mass flux, read in from pbuf.  Since shallow convection is never 
         !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
         !  fraction is purely from deep convection scheme.  
         deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(i,k+1)-cmfmc_sh(i,k+1))),0.6_r8))

         if (do_clubb_mf_rad) then
           shalcu(i,k) = mf_cloudfrac_output(i,k)
           sh_icwmr(i,k) = mf_qc_output(i,k)
         end if
       
         if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
            deepcu(i,k) = 0._r8
         endif

         if (shalcu(i,k) <= frac_limit .or. sh_icwmr(i,k) < ic_limit) then
            shalcu(i,k) = 0._r8
         endif
             
         !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable 
         !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
         !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud 
         !  from CLUBB plus the deep convective cloud fraction
         concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+deepcu(i,k)+shalcu(i,k),0.80_r8)
      enddo
    enddo
   
    if (single_column) then
      if (trim(scm_clubb_iop_name)  ==  'ATEX_48hr'       .or. &
          trim(scm_clubb_iop_name)  ==  'BOMEX_5day'      .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf01_4day' .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf02_06hr' .or. &
          trim(scm_clubb_iop_name)  ==  'RICO_3day'       .or. &
          trim(scm_clubb_iop_name)  ==  'ARM_CC') then
       
         deepcu(:,:) = 0.0_r8
         concld(:,:) = 0.0_r8
   
      endif       
    endif
   
    ! --------------------------------------------------------------------------------- !  
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !  

    aist(:,:top_lev-1) = 0._r8
    qsatfac(:, :) = 0._r8 ! Zero out entire profile in case qsatfac is left undefined in aist_vector below 

    do k = top_lev, pver

      ! For Type II PSC and for thin cirrus, the clouds can be thin, but
      ! extensive and they should start forming when the gridbox mean saturation
      ! reaches 1.0.
      !
      ! For now, use the tropopause diagnostic to determine where the Type II
      ! PSC should be, but in the future wold like a better metric that can also
      ! identify the level for thin cirrus. Include the tropopause level so that
      ! the cold point tropopause will use the stratospheric values.
      where (k <= troplev)
        rhmini = rhminis_const
        rhmaxi = rhmaxis_const
      elsewhere
        rhmini = rhmini_const
        rhmaxi = rhmaxi_const
      end where

      if ( trim(subcol_scheme) == 'SILHS' ) then
        call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
             state1%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol )
      else
        call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
              state1%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol,&
              qsatfac_out=qsatfac(:,k), rhmini_in=rhmini, rhmaxi_in=rhmaxi)
      endif
    enddo
  
    ! --------------------------------------------------------------------------------- !  
    !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
    !                                                                                   !
    !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
    !  because CLUBB does nothing with ice.  Here I simply overwrite the liquid stratus ! 
    !  fraction that was coded in macrop_driver                                         !
    ! --------------------------------------------------------------------------------- !  
 
    !  Recompute net stratus fraction using maximum over-lapping assumption, as done
    !  in macrophysics code, using alst computed above and aist read in from physics buffer
   
    do k=1,pver
      do i=1,ncol
        ast(i,k) = max(alst(i,k),aist(i,k))
        qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k)) 
      enddo
    enddo
   
   !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just 
   !  be outputting the shallow convective cloud fraction 
    do k=1,pver
      do i=1,ncol
         cloud_frac(i,k) = min(ast(i,k)+deepcu(i,k)+shalcu(i,k),1.0_r8)
      enddo
    enddo
   
    ! --------------------------------------------------------------------------------- !  
    !  DIAGNOSE THE PBL DEPTH                                                           !
    !  this is needed for aerosol code                                                  !
    ! --------------------------------------------------------------------------------- !       
    do i=1,ncol
      do k=1,pver
         !use local exner since state%exner is not a proper exner
         th(i,k) = state1%t(i,k)*inv_exner_clubb(i,k)
         !thv should have condensate loading to be consistent with earlier def's in this module
         thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq) - state1%q(i,k,ixcldliq))
      enddo
    enddo
 
    ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
    rrho(1:ncol) = (rga)*(state1%pdel(1:ncol,pver)/dz_g(1:ncol,pver))
    call calc_ustar( ncol, state1%t(1:ncol,pver), state1%pmid(1:ncol,pver), cam_in%wsx(1:ncol), cam_in%wsy(1:ncol), &
                    rrho(1:ncol), ustar2(1:ncol))
    ! use correct qflux from coupler
    call calc_obklen( ncol, th(1:ncol,pver), thv(1:ncol,pver), cam_in%cflx(1:ncol,1), cam_in%shf(1:ncol), &
                      rrho(1:ncol), ustar2(1:ncol), kinheat(1:ncol), kinwat(1:ncol), kbfs(1:ncol), &
                      obklen(1:ncol))
   
    dummy2(:) = 0._r8
    dummy3(:) = 0._r8
   
    where (kbfs(:ncol)  ==  -0.0_r8) kbfs(:ncol) = 0.0_r8

    !  Compute PBL depth according to Holtslag-Boville Scheme
    call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kbfs, pblh, dummy2, &
                state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)

    !  Output the PBL depth
    call outfld('PBLH', pblh, pcols, lchnk)
 
    ! Assign the first pver levels of cloud_frac back to cld
    cld(:,1:pver) = cloud_frac(:,1:pver)

    ! --------------------------------------------------------------------------------- !   
    !  END CLOUD FRACTION DIAGNOSIS, begin to store variables back into buffer          !
    ! --------------------------------------------------------------------------------- !  
 
    !  Output calls of variables goes here
    call outfld( 'RELVAR',           relvar,                pcols, lchnk )
    call outfld( 'RHO_CLUBB',        rho(:,1:pver),         pcols, lchnk )
    call outfld( 'WP2_CLUBB',        wp2,                   pcols, lchnk )
    call outfld( 'UP2_CLUBB',        up2,                   pcols, lchnk )
    call outfld( 'VP2_CLUBB',        vp2,                   pcols, lchnk )
    call outfld( 'WP3_CLUBB',        wp3_output(:,1:pver),  pcols, lchnk )
    call outfld( 'UPWP_CLUBB',       upwp,                  pcols, lchnk )
    call outfld( 'VPWP_CLUBB',       vpwp,                  pcols, lchnk )
    call outfld( 'WPTHLP_CLUBB',     wpthlp_output,         pcols, lchnk )
    call outfld( 'WPRTP_CLUBB',      wprtp_output,          pcols, lchnk )
    call outfld( 'RTP2_CLUBB',       rtp2,                  pcols, lchnk )
    call outfld( 'RTPTHLP_CLUBB',    rtpthlp_output,        pcols, lchnk )
    call outfld( 'RCM_CLUBB',        rcm(:,1:pver),         pcols, lchnk )
    call outfld( 'RTM_CLUBB',        rtm(:,1:pver),         pcols, lchnk )
    call outfld( 'THLM_CLUBB',       thlm(:,1:pver),        pcols, lchnk )

    temp2dp(:ncol,:) = wprcp(:ncol,:) * latvap
    call outfld( 'WPRCP_CLUBB',      temp2dp,                 pcols, lchnk )

    temp2dp(:ncol,:) = wpthvp(:ncol,:) * cpair
    call outfld( 'WPTHVP_CLUBB',     temp2dp,                 pcols, lchnk )

    call outfld( 'RTP2_ZT_CLUBB',    rtp2_zt_out(:,1:pver),   pcols, lchnk )
    call outfld( 'THLP2_ZT_CLUBB',   thl2_zt_out(:,1:pver),   pcols, lchnk )
    call outfld( 'WP2_ZT_CLUBB',     wp2_zt_out(:,1:pver),    pcols, lchnk )
    call outfld( 'PDFP_RTP2_CLUBB',  pdfp_rtp2,               pcols, lchnk )
    call outfld( 'THLP2_CLUBB',      thlp2,                   pcols, lchnk )
    call outfld( 'RCMINLAYER_CLUBB', rcm_in_layer(:,1:pver),  pcols, lchnk )
    call outfld( 'CLOUDFRAC_CLUBB',  alst,                    pcols, lchnk )
    call outfld( 'CLOUDCOVER_CLUBB', cloud_frac(:,1:pver),    pcols, lchnk )
    call outfld( 'ZT_CLUBB',         zt_out(:,1:pver),        pcols, lchnk )
    call outfld( 'ZM_CLUBB',         zi_out,                  pcols, lchnk )
    call outfld( 'UM_CLUBB',         um(:,1:pver),            pcols, lchnk )
    call outfld( 'VM_CLUBB',         vm(:,1:pver),            pcols, lchnk )
    call outfld( 'WM_ZT_CLUBB',      wm_zt_out(:,1:pver),     pcols, lchnk )
    call outfld( 'CONCLD',           concld,                  pcols, lchnk )
    call outfld( 'DP_CLD',           deepcu,                  pcols, lchnk )
    call outfld( 'ZMDLF',            dlf_liq_out,             pcols, lchnk )
    call outfld( 'ZMDLFI',           dlf_ice_out,             pcols, lchnk )
    call outfld( 'CLUBB_GRID_SIZE',  grid_dx,                 pcols, lchnk )
    call outfld( 'QSATFAC',          qsatfac,                 pcols, lchnk )

    call outfld( 'PRECSH' ,          prec_sh(:ncol),          pcols, lchnk )
    call outfld( 'TKE_CLUBB',        tke,                     pcols, lchnk )
 
   ! --------------------------------------------------------------- !
   ! Writing state variables after EDMF scheme for detailed analysis !
   ! --------------------------------------------------------------- !
    if (do_clubb_mf) then
      call outfld( 'edmf_DRY_A'    , mf_dry_a_output,           pcols, lchnk )
      call outfld( 'edmf_MOIST_A'  , mf_moist_a_output,         pcols, lchnk )
      call outfld( 'edmf_DRY_W'    , mf_dry_w_output,           pcols, lchnk )
      call outfld( 'edmf_MOIST_W'  , mf_moist_w_output,         pcols, lchnk )
      call outfld( 'edmf_DRY_QT'   , mf_dry_qt_output,          pcols, lchnk )
      call outfld( 'edmf_MOIST_QT' , mf_moist_qt_output,        pcols, lchnk )
      call outfld( 'edmf_DRY_THL'  , mf_dry_thl_output,         pcols, lchnk )
      call outfld( 'edmf_MOIST_THL', mf_moist_thl_output,       pcols, lchnk )
      call outfld( 'edmf_DRY_U'    , mf_dry_u_output,           pcols, lchnk )
      call outfld( 'edmf_MOIST_U'  , mf_moist_u_output,         pcols, lchnk )
      call outfld( 'edmf_DRY_V'    , mf_dry_v_output,           pcols, lchnk )
      call outfld( 'edmf_MOIST_V'  , mf_moist_v_output,         pcols, lchnk )
      call outfld( 'edmf_MOIST_QC' , mf_moist_qc_output,        pcols, lchnk )
      call outfld( 'edmf_precc'    , mf_precc_output,           pcols, lchnk )
      call outfld( 'edmf_S_AE'     , s_ae_output,               pcols, lchnk )
      call outfld( 'edmf_S_AW'     , s_aw_output,               pcols, lchnk )
      call outfld( 'edmf_S_AWW'    , s_aww_output,              pcols, lchnk )
      call outfld( 'edmf_S_AWTHL'  , s_awthl_output,            pcols, lchnk )
      call outfld( 'edmf_S_AWQT'   , s_awqt_output,             pcols, lchnk )
      call outfld( 'edmf_S_AWU'    , s_awu_output,              pcols, lchnk )
      call outfld( 'edmf_S_AWV'    , s_awv_output,              pcols, lchnk )
      call outfld( 'edmf_thlforcup', mf_thlforcup_output,       pcols, lchnk )
      call outfld( 'edmf_qtforcup' , mf_qtforcup_output,        pcols, lchnk )
      call outfld( 'edmf_thlforcdn', mf_thlforcdn_output,       pcols, lchnk )
      call outfld( 'edmf_qtforcdn' , mf_qtforcdn_output,        pcols, lchnk )
      call outfld( 'edmf_thlforc'  , mf_thlforc_output,         pcols, lchnk )
      call outfld( 'edmf_qtforc'   , mf_qtforc_output,          pcols, lchnk )
      call outfld( 'edmf_thlflxup' , mf_thlflxup_output,        pcols, lchnk )
      call outfld( 'edmf_qtflxup'  , mf_qtflxup_output,         pcols, lchnk )
      call outfld( 'edmf_thlflxdn' , mf_thlflxdn_output,        pcols, lchnk )
      call outfld( 'edmf_qtflxdn'  , mf_qtflxdn_output,         pcols, lchnk )
      call outfld( 'edmf_thlflx'   , mf_thlflx_output,          pcols, lchnk )
      call outfld( 'edmf_qtflx'    , mf_qtflx_output,           pcols, lchnk )
      call outfld( 'edmf_thvflx'   , mf_thvflx_output,          pcols, lchnk )
      call outfld( 'edmf_rcm'      , mf_rcm_output,             pcols, lchnk )
      call outfld( 'edmf_uflxup'   , mf_uflxup_output,          pcols, lchnk )
      call outfld( 'edmf_vflxup'   , mf_vflxup_output,          pcols, lchnk )
      call outfld( 'edmf_uflxdn'   , mf_uflxdn_output,          pcols, lchnk )
      call outfld( 'edmf_vflxdn'   , mf_vflxdn_output,          pcols, lchnk )
      call outfld( 'edmf_uflx'     , mf_uflx_output,            pcols, lchnk )
      call outfld( 'edmf_vflx'     , mf_vflx_output,            pcols, lchnk )
      call outfld( 'edmf_cloudfrac', mf_cloudfrac_output,       pcols, lchnk )
      call outfld( 'edmf_ent'      , mf_ent_output,             pcols, lchnk )
      call outfld( 'edmf_upa'      , mf_upa_output,             pcols, lchnk )
      call outfld( 'edmf_upw'      , mf_upw_output,             pcols, lchnk )
      call outfld( 'edmf_upmf'     , mf_upmf_output,            pcols, lchnk )
      call outfld( 'edmf_upqt'     , mf_upqt_output,            pcols, lchnk )
      call outfld( 'edmf_upthl'    , mf_upthl_output,           pcols, lchnk )
      call outfld( 'edmf_upthv'    , mf_upthv_output,           pcols, lchnk )
      call outfld( 'edmf_upth'     , mf_upth_output,            pcols, lchnk )
      call outfld( 'edmf_upqc'     , mf_upqc_output,            pcols, lchnk )
      call outfld( 'edmf_upbuoy'   , mf_upbuoy_output,          pcols, lchnk )
      call outfld( 'edmf_upent'    , mf_upent_output,           pcols, lchnk )
      call outfld( 'edmf_updet'    , mf_updet_output,           pcols, lchnk )
      call outfld( 'edmf_dnw'      , mf_dnw_output,             pcols, lchnk )
      call outfld( 'edmf_dnthl'    , mf_dnthl_output,           pcols, lchnk )
      call outfld( 'edmf_dnqt'     , mf_dnqt_output,            pcols, lchnk )
      call outfld( 'edmf_sqtup'    , mf_sqtup_output,           pcols, lchnk )
      call outfld( 'edmf_sqtdn'    , mf_sqtdn_output,           pcols, lchnk )

      ! macmic_it==1 ensures that this is ddcp aeraged over the prior time-steps
      if (macmic_it==1) call outfld( 'edmf_ztop'     , mf_ztop_output,            pcols, lchnk )
      if (macmic_it==1) call outfld( 'edmf_ddcp'     , mf_ddcp_output,            pcols, lchnk )

      call outfld( 'edmf_L0'       , mf_L0_output,              pcols, lchnk )
      call outfld( 'edmf_freq'     , mf_freq_output,            pcols, lchnk )
      call outfld( 'edmf_cape'     , mf_cape_output,            pcols, lchnk )
      call outfld( 'edmf_cfl'      , mf_cfl_output,             pcols, lchnk )
      call outfld( 'ICWMRSH'       , sh_icwmr,                  pcols, lchnk )
    end if

    if (macmic_it==cld_macmic_num_steps) then
      call outfld( 'qtm_macmic1'         , qtm_macmic1,            pcols, lchnk )
      call outfld( 'qtm_macmic2'         , qtm_macmic2,            pcols, lchnk )
      call outfld( 'thlm_macmic1'        , thlm_macmic1,           pcols, lchnk )
      call outfld( 'thlm_macmic2'        , thlm_macmic2,           pcols, lchnk )
 
      call outfld( 'RCM_CLUBB_macmic'     , rcm_macmic,           pcols, lchnk )
      call outfld( 'CLDFRAC_CLUBB_macmic' , cldfrac_macmic,       pcols, lchnk )
 
      call outfld( 'WPTHLP_CLUBB_macmic'  , wpthlp_macmic,        pcols, lchnk )
      call outfld( 'WPRTP_CLUBB_macmic'   , wprtp_macmic,         pcols, lchnk )
      call outfld( 'WPTHVP_CLUBB_macmic'  , wpthvp_macmic,        pcols, lchnk )
      if (do_clubb_mf) then
        call outfld( 'edmf_thlflx_macmic'   , mf_thlflx_macmic,     pcols, lchnk )
        call outfld( 'edmf_qtflx_macmic'    , mf_qtflx_macmic,      pcols, lchnk )
        call outfld( 'edmf_thvflx_macmic'   , mf_thvflx_macmic,     pcols, lchnk )
        call outfld( 'up_macmic1'   , up_macmic1,     pcols, lchnk )
        call outfld( 'up_macmic2'   , up_macmic2,     pcols, lchnk )
        call outfld( 'dn_macmic1'   , dn_macmic1,     pcols, lchnk )
        call outfld( 'dn_macmic2'   , dn_macmic2,     pcols, lchnk )
        call outfld( 'upa_macmic1'  , upa_macmic1,    pcols, lchnk )
        call outfld( 'upa_macmic2'  , upa_macmic2,    pcols, lchnk )
        call outfld( 'dna_macmic1'  , dna_macmic1,    pcols, lchnk )
        call outfld( 'dna_macmic2'  , dna_macmic2,    pcols, lchnk )
        call outfld( 'thlu_macmic1' , thlu_macmic1,   pcols, lchnk )
        call outfld( 'thlu_macmic2' , thlu_macmic2,   pcols, lchnk )
        call outfld( 'qtu_macmic1'  , qtu_macmic1,    pcols, lchnk )
        call outfld( 'qtu_macmic2'  , qtu_macmic2,    pcols, lchnk )
        call outfld( 'thld_macmic1' , thld_macmic1,   pcols, lchnk )
        call outfld( 'thld_macmic2' , thld_macmic2,   pcols, lchnk )
        call outfld( 'qtd_macmic1'  , qtd_macmic1,    pcols, lchnk )
        call outfld( 'qtd_macmic2'  , qtd_macmic2,    pcols, lchnk )
        call outfld( 'dthl_macmic1' , dthl_macmic1,   pcols, lchnk )
        call outfld( 'dthl_macmic2' , dthl_macmic2,   pcols, lchnk )
        call outfld( 'dqt_macmic1'  , dqt_macmic1,    pcols, lchnk )
        call outfld( 'dqt_macmic2'  , dqt_macmic2,    pcols, lchnk )
        call outfld( 'dthlu_macmic1' , dthlu_macmic1,   pcols, lchnk )
        call outfld( 'dthlu_macmic2' , dthlu_macmic2,   pcols, lchnk )
        call outfld( 'dqtu_macmic1'  , dqtu_macmic1,    pcols, lchnk )
        call outfld( 'dqtu_macmic2'  , dqtu_macmic2,    pcols, lchnk )
        call outfld( 'dthld_macmic1' , dthld_macmic1,   pcols, lchnk )
        call outfld( 'dthld_macmic2' , dthld_macmic2,   pcols, lchnk )
        call outfld( 'dqtd_macmic1'  , dqtd_macmic1,    pcols, lchnk )
        call outfld( 'dqtd_macmic2'  , dqtd_macmic2,    pcols, lchnk )
        call outfld( 'ztop_macmic1'  , ztop_macmic1,    pcols, lchnk )
        call outfld( 'ztop_macmic2'  , ztop_macmic2,    pcols, lchnk )
        call outfld( 'ddcp_macmic1'  , ddcp_macmic1,    pcols, lchnk )
        call outfld( 'ddcp_macmic2'  , ddcp_macmic2,    pcols, lchnk )
      end if
    end if

    !  Output CLUBB history here
    if (l_stats) then 
      
      do j=1,stats_zt(1)%num_output_fields
   
        temp1 = trim(stats_zt(1)%file%grid_avg_var(j)%name)
        sub   = temp1
        if (len(temp1) >  max_fieldname_len) sub = temp1(1:max_fieldname_len)
 
        call outfld(trim(sub), out_zt(:,:,j), pcols, lchnk )
      enddo
   
      do j=1,stats_zm(1)%num_output_fields
   
        temp1 = trim(stats_zm(1)%file%grid_avg_var(j)%name)
        sub   = temp1
        if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
   
        call outfld(trim(sub),out_zm(:,:,j), pcols, lchnk)
      enddo

      if (l_output_rad_files) then  
        do j=1,stats_rad_zt(1)%num_output_fields
          call outfld(trim(stats_rad_zt(1)%file%grid_avg_var(j)%name), out_radzt(:,:,j), pcols, lchnk)
        enddo
   
        do j=1,stats_rad_zm(1)%num_output_fields
          call outfld(trim(stats_rad_zm(1)%file%grid_avg_var(j)%name), out_radzm(:,:,j), pcols, lchnk)
        enddo
      endif
   
      do j=1,stats_sfc(1)%num_output_fields
        call outfld(trim(stats_sfc(1)%file%grid_avg_var(j)%name), out_sfc(:,:,j), pcols, lchnk)
      enddo
   
    endif
   
    call t_stopf("clubb_tend_cam")
   
    return
#endif
  end subroutine clubb_tend_cam

  subroutine clubb_emissions_cam (state, cam_in, ptend)

  !-------------------------------------------------------------------------------
  ! Description: Apply surface fluxes of constituents to lowest model level
  !              except water vapor (applied in clubb_tend_cam)
  !
  ! Author: Adam Herrington, November 2022
  ! Origin: Based on E3SM's clubb_surface subroutine
  ! References:
  !   None
  !-------------------------------------------------------------------------------
  use physics_types,      only: physics_ptend, physics_ptend_init, physics_state
  use constituents,       only: cnst_type
  use camsrfexch,         only: cam_in_t

  ! --------------- !
  ! Input Arguments !
  ! --------------- !
  type(physics_state), intent(in)  :: state                     ! Physics state variables
  type(cam_in_t),      intent(in)  :: cam_in                    ! Surface inputs

  ! ---------------------- !
  ! Output Arguments       !
  ! ---------------------- !
  type(physics_ptend), intent(out) :: ptend                      ! Individual parameterization tendencies

  ! --------------- !
  ! Local Variables !
  ! --------------- !
  integer  :: m, ncol                                 
  logical  :: lq(pcnst)

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !
  ncol = state%ncol

  lq(1) = .false.
  lq(2:) = .true.
  call physics_ptend_init(ptend,state%psetcols, "clubb emissions", lq=lq)

  ! Apply tracer fluxes to lowest model level (except water vapor)
  do m = 2,pcnst
    ptend%q(:ncol,pver,m) = cam_in%cflx(:ncol,m)*state%rpdel(:ncol,pver)*gravit
  end do

  ! Convert tendencies of dry constituents to dry basis.
  do m = 2,pcnst
     if (cnst_type(m).eq.'dry') then
        ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m)*state%pdel(:ncol,pver)*state%rpdeldry(:ncol,pver)
     endif
  end do

  end subroutine clubb_emissions_cam  

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

! Saturation adjustment for ice
! Add ice mass if supersaturated
subroutine ice_macro_tend(naai,t,p,qv,qi,ni,xxls,deltat,stend,qvtend,qitend,nitend,vlen) 

  use wv_sat_methods, only: wv_sat_qsat_ice

  integer,                   intent(in)  :: vlen
  real(r8), dimension(vlen), intent(in)  :: naai   !Activated number of ice nuclei 
  real(r8), dimension(vlen), intent(in)  :: t      !temperature (k)
  real(r8), dimension(vlen), intent(in)  :: p      !pressure (pa)
  real(r8), dimension(vlen), intent(in)  :: qv     !water vapor mixing ratio
  real(r8), dimension(vlen), intent(in)  :: qi     !ice mixing ratio
  real(r8), dimension(vlen), intent(in)  :: ni     !ice number concentration
  real(r8),                  intent(in)  :: xxls   !latent heat of freezing
  real(r8),                  intent(in)  :: deltat !timestep
  real(r8), dimension(vlen), intent(out) :: stend  ! 'temperature' tendency 
  real(r8), dimension(vlen), intent(out) :: qvtend !vapor tendency
  real(r8), dimension(vlen), intent(out) :: qitend !ice mass tendency
  real(r8), dimension(vlen), intent(out) :: nitend !ice number tendency  
 
  real(r8) :: ESI(vlen)
  real(r8) :: QSI(vlen)
  integer  :: i

  do i = 1, vlen
     stend(i)  = 0._r8
     qvtend(i) = 0._r8
     qitend(i) = 0._r8
     nitend(i) = 0._r8
  end do

! calculate qsati from t,p,q
  do i = 1, vlen
     call wv_sat_qsat_ice(t(i), p(i), ESI(i), QSI(i))
  end do

  do i = 1, vlen
     if (naai(i) > 1.e-18_r8 .and. qv(i) > QSI(i)) then

        qitend(i) = (qv(i)-QSI(i))/deltat
        qvtend(i) = 0._r8 - qitend(i)
        stend(i)  = qitend(i) * xxls      ! moist static energy tend...[J/kg/s] !
   
        ! if ice exists (more than 1 L-1) and there is condensation, do not add to number (= growth), else, add 10um ice
        if (ni(i) < 1.e3_r8 .and. (qi(i)+qitend(i)*deltat) > 1.e-18_r8) then
           nitend(i) = nitend(i) + 3._r8 * qitend(i)/(4._r8*3.14_r8* 10.e-6_r8**3*997._r8)
        end if

     end if
  end do

end subroutine ice_macro_tend

#ifdef CLUBB_SGS
! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!

real(r8) function diag_ustar( z, bflx, wnd, z0 ) 

use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

implicit none

real(r8), parameter      :: am   =  4.8_r8   !   "          "         "
real(r8), parameter      :: bm   = 19.3_r8  !   "          "         "

real(r8), parameter      :: grav = shr_const_g
real(r8), parameter      :: vonk = shr_const_karman
real(r8), parameter      :: pi   = shr_const_pi

real(r8), intent (in)    :: z             ! height where u locates
real(r8), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real(r8), intent (in)    :: wnd           ! wind speed at z
real(r8), intent (in)    :: z0            ! momentum roughness height


integer :: iterate
real(r8)    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0_r8 - 3.0_r8*log( 2.0_r8 )

ustar =  wnd*klnz
if (abs(bflx) > 1.e-6_r8) then
   do iterate=1,4

      if (ustar > 1.e-6_r8) then
         lmo   = -ustar**3 / ( vonk * bflx )
         zeta  = z/lmo
         if (zeta > 0._r8) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
         else
            x     = sqrt( sqrt( 1.0_r8 - bm*zeta ) )
            psi1  = 2._r8*log( 1.0_r8+x ) + log( 1.0_r8+x*x ) - 2._r8*atan( x ) + c1
            ustar = wnd*vonk/(lnz - psi1)
         end if

      endif

   end do
end if


diag_ustar = ustar

return


end function diag_ustar
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS

  subroutine stats_init_clubb( l_stats_in, stats_tsamp_in, stats_tout_in, &
                               nnzp, nnrad_zt,nnrad_zm, delt, &
                               stats_zt, stats_zm, stats_sfc, &
                               stats_rad_zt, stats_rad_zm)
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.   
    
    !-----------------------------------------------------------------------


    use clubb_api_module, only: &
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      ztscr07, & 
      ztscr08, & 
      ztscr09, & 
      ztscr10, & 
      ztscr11, & 
      ztscr12, & 
      ztscr13, & 
      ztscr14, & 
      ztscr15, & 
      ztscr16, & 
      ztscr17, & 
      ztscr18, & 
      ztscr19, & 
      ztscr20, & 
      ztscr21

    use clubb_api_module, only: &
      zmscr01, & 
      zmscr02, & 
      zmscr03, & 
      zmscr04, & 
      zmscr05, & 
      zmscr06, & 
      zmscr07, & 
      zmscr08, & 
      zmscr09, & 
      zmscr10, & 
      zmscr11, & 
      zmscr12, & 
      zmscr13, & 
      zmscr14, & 
      zmscr15, &
      zmscr16, &
      zmscr17, &
      l_stats, &
      l_output_rad_files, & 
      stats_tsamp,   & 
      stats_tout,    & 
      l_stats_samp,  & 
      l_stats_last, & 
      l_netcdf, & 
      l_grads

    use clubb_api_module, only:        time_precision, &   !
                                       nvarmax_zm, stats_init_zm_api, & !
                                       nvarmax_zt, stats_init_zt_api, & !
                                       nvarmax_rad_zt, stats_init_rad_zt_api, & !
                                       nvarmax_rad_zm, stats_init_rad_zm_api, & !
                                       nvarmax_sfc, stats_init_sfc_api, & !
                                       fstderr, var_length !
    use cam_abortutils,         only: endrun
    use cam_history,            only: addfld, horiz_only
    use namelist_utils,         only: find_group_name
    use units,                  only: getunit, freeunit
    use spmd_utils,             only: mpicom, mstrid=>masterprocid, mpi_character

    implicit none

    ! Input Variables

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count] 
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::   delt         ! Timestep (dtmain in CLUBB)         [s]
    
    ! Output Variables
    type (stats), intent(out) :: stats_zt,      & ! stats_zt grid
                                 stats_zm,      & ! stats_zm grid
                                 stats_rad_zt,  & ! stats_rad_zt grid
                                 stats_rad_zm,  & ! stats_rad_zm grid
                                 stats_sfc        ! stats_sfc


    !  Local Variables

    !  Namelist Variables

    character(len=*), parameter :: subr = 'stats_init_clubb'

    character(len=var_length), dimension(nvarmax_zt)     ::   clubb_vars_zt      ! Variables on the thermodynamic levels
    character(len=var_length), dimension(nvarmax_zm)     ::   clubb_vars_zm      ! Variables on the momentum levels
    character(len=var_length), dimension(nvarmax_rad_zt) ::   clubb_vars_rad_zt  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_rad_zm) ::   clubb_vars_rad_zm  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_sfc)    ::   clubb_vars_sfc     ! Variables at the model surface

    namelist /clubb_stats_nl/ & 
      clubb_vars_zt, & 
      clubb_vars_zm, &
      clubb_vars_rad_zt, &
      clubb_vars_rad_zm, & 
      clubb_vars_sfc

    !  Local Variables

    logical :: l_error, &
               first_call = .false.

    character(len=200) :: temp1, sub

    integer :: i, ntot, read_status
    integer :: iunit, ierr

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in
    
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in

    if ( .not. l_stats ) then
       l_stats_samp  = .false.
       l_stats_last  = .false.
       return
    end if

    !  Initialize namelist variables

    clubb_vars_zt     = ''
    clubb_vars_zm     = ''
    clubb_vars_rad_zt = ''
    clubb_vars_rad_zm = ''
    clubb_vars_sfc    = ''

    !  Read variables to compute from the namelist    
    if (masterproc) then
       iunit= getunit()
       open(unit=iunit,file="atm_in",status='old')
       call find_group_name(iunit, 'clubb_stats_nl', status=read_status)
       if (read_status == 0) then
          read(unit=iunit, nml=clubb_stats_nl, iostat=read_status)
          if (read_status /= 0) then
             call endrun('stats_init_clubb:  error reading namelist')
          end if
       end if
       close(unit=iunit)
       call freeunit(iunit)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(clubb_vars_zt,      var_length*nvarmax_zt,       mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_zt")
    call mpi_bcast(clubb_vars_zm,      var_length*nvarmax_zm,       mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_zm")
    call mpi_bcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_rad_zt")
    call mpi_bcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_rad_zm")
    call mpi_bcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_sfc")

    !  Hardcode these for use in CAM-CLUBB, don't want either
    l_netcdf = .false.
    l_grads  = .false.

    !  Check sampling and output frequencies

    !  The model time step length, delt (which is dtmain), should multiply
    !  evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) ) > 1.e-8_r8 ) then
       l_error = .true.  ! This will cause the run to stop.
       write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                        'delt (which is dtmain).  Check the appropriate ',  &
                        'model.in file.'
       write(fstderr,*) 'stats_tsamp = ', stats_tsamp
       write(fstderr,*) 'delt = ', delt
    endif

    !  Initialize zt (mass points)

    i = 1
    do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0 .and. & 
               len_trim(clubb_vars_zt(i))   /= 0 .and. & 
               i <= nvarmax_zt )
       i = i + 1
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zt than allowed for by nvarmax_zt."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                        "in the stats namelist, or change nvarmax_zt."
       write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
       call endrun ("stats_init_clubb:  number of zt statistical variables exceeds limit")
    endif

    stats_zt%num_output_fields = ntot
    stats_zt%kk = nnzp

    allocate( stats_zt%z( stats_zt%kk ) )

    allocate( stats_zt%accum_field_values( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%accum_num_samples( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%l_in_update( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )

    allocate( stats_zt%file%grid_avg_var( stats_zt%num_output_fields ) )
    allocate( stats_zt%file%z( stats_zt%kk ) )

    first_call = (.not. allocated(ztscr01))

    !  Allocate scratch space
    if (first_call) allocate( ztscr01(stats_zt%kk) )
    if (first_call) allocate( ztscr02(stats_zt%kk) )
    if (first_call) allocate( ztscr03(stats_zt%kk) )
    if (first_call) allocate( ztscr04(stats_zt%kk) )
    if (first_call) allocate( ztscr05(stats_zt%kk) )
    if (first_call) allocate( ztscr06(stats_zt%kk) )
    if (first_call) allocate( ztscr07(stats_zt%kk) )
    if (first_call) allocate( ztscr08(stats_zt%kk) )
    if (first_call) allocate( ztscr09(stats_zt%kk) )
    if (first_call) allocate( ztscr10(stats_zt%kk) )
    if (first_call) allocate( ztscr11(stats_zt%kk) )
    if (first_call) allocate( ztscr12(stats_zt%kk) )
    if (first_call) allocate( ztscr13(stats_zt%kk) )
    if (first_call) allocate( ztscr14(stats_zt%kk) )
    if (first_call) allocate( ztscr15(stats_zt%kk) )
    if (first_call) allocate( ztscr16(stats_zt%kk) )
    if (first_call) allocate( ztscr17(stats_zt%kk) )
    if (first_call) allocate( ztscr18(stats_zt%kk) )
    if (first_call) allocate( ztscr19(stats_zt%kk) )
    if (first_call) allocate( ztscr20(stats_zt%kk) )
    if (first_call) allocate( ztscr21(stats_zt%kk) )

    ztscr01 = 0.0_r8
    ztscr02 = 0.0_r8
    ztscr03 = 0.0_r8
    ztscr04 = 0.0_r8
    ztscr05 = 0.0_r8
    ztscr06 = 0.0_r8
    ztscr07 = 0.0_r8
    ztscr08 = 0.0_r8
    ztscr09 = 0.0_r8
    ztscr10 = 0.0_r8
    ztscr11 = 0.0_r8
    ztscr12 = 0.0_r8
    ztscr13 = 0.0_r8
    ztscr14 = 0.0_r8
    ztscr15 = 0.0_r8
    ztscr16 = 0.0_r8
    ztscr17 = 0.0_r8
    ztscr18 = 0.0_r8
    ztscr19 = 0.0_r8
    ztscr20 = 0.0_r8
    ztscr21 = 0.0_r8

    !  Default initialization for array indices for zt
    if (first_call) then
      call stats_init_zt_api( clubb_vars_zt, l_error, &
                              stats_zt )
    end if

    !  Initialize zm (momentum points)

    i = 1
    do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  .and. & 
               len_trim(clubb_vars_zm(i)) /= 0    .and. & 
               i <= nvarmax_zm )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zm than allowed for by nvarmax_zm."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                        "in the stats namelist, or change nvarmax_zm."
       write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
       call endrun ("stats_init_clubb:  number of zm statistical variables exceeds limit")
    endif

    stats_zm%num_output_fields = ntot
    stats_zm%kk = nnzp

    allocate( stats_zm%z( stats_zm%kk ) )

    allocate( stats_zm%accum_field_values( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%accum_num_samples( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%l_in_update( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )

    allocate( stats_zm%file%grid_avg_var( stats_zm%num_output_fields ) )
    allocate( stats_zm%file%z( stats_zm%kk ) )

    !  Allocate scratch space

    if (first_call) allocate( zmscr01(stats_zm%kk) )
    if (first_call) allocate( zmscr02(stats_zm%kk) )
    if (first_call) allocate( zmscr03(stats_zm%kk) )
    if (first_call) allocate( zmscr04(stats_zm%kk) )
    if (first_call) allocate( zmscr05(stats_zm%kk) )
    if (first_call) allocate( zmscr06(stats_zm%kk) )
    if (first_call) allocate( zmscr07(stats_zm%kk) )
    if (first_call) allocate( zmscr08(stats_zm%kk) )
    if (first_call) allocate( zmscr09(stats_zm%kk) )
    if (first_call) allocate( zmscr10(stats_zm%kk) )
    if (first_call) allocate( zmscr11(stats_zm%kk) )
    if (first_call) allocate( zmscr12(stats_zm%kk) )
    if (first_call) allocate( zmscr13(stats_zm%kk) )
    if (first_call) allocate( zmscr14(stats_zm%kk) )
    if (first_call) allocate( zmscr15(stats_zm%kk) )
    if (first_call) allocate( zmscr16(stats_zm%kk) )
    if (first_call) allocate( zmscr17(stats_zm%kk) )

    zmscr01 = 0.0_r8
    zmscr02 = 0.0_r8
    zmscr03 = 0.0_r8
    zmscr04 = 0.0_r8
    zmscr05 = 0.0_r8
    zmscr06 = 0.0_r8
    zmscr07 = 0.0_r8
    zmscr08 = 0.0_r8
    zmscr09 = 0.0_r8
    zmscr10 = 0.0_r8
    zmscr11 = 0.0_r8
    zmscr12 = 0.0_r8
    zmscr13 = 0.0_r8
    zmscr14 = 0.0_r8
    zmscr15 = 0.0_r8
    zmscr16 = 0.0_r8
    zmscr17 = 0.0_r8

    if (first_call) then
      call stats_init_zm_api( clubb_vars_zm, l_error, &
                              stats_zm )
    end if

    !  Initialize rad_zt (radiation points)

    if (l_output_rad_files) then
    
       i = 1
       do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  .and. & 
                  len_trim(clubb_vars_rad_zt(i))   /= 0  .and. & 
                  i <= nvarmax_rad_zt )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zt ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                           "in the stats namelist, or change nvarmax_rad_zt."
          write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
          call endrun ("stats_init_clubb:  number of rad_zt statistical variables exceeds limit")
       endif

      stats_rad_zt%num_output_fields = ntot
      stats_rad_zt%kk = nnrad_zt

      allocate( stats_rad_zt%z( stats_rad_zt%kk ) )

      allocate( stats_rad_zt%accum_field_values( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%accum_num_samples( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%l_in_update( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )

      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                     stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )

      allocate( stats_rad_zt%file%grid_avg_var( stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%file%z( stats_rad_zt%kk ) )

       call stats_init_rad_zt_api( clubb_vars_rad_zt, l_error, &
                                   stats_rad_zt )

       !  Initialize rad_zm (radiation points)
 
       i = 1
       do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0 .and. & 
                  len_trim(clubb_vars_rad_zm(i))   /= 0 .and. & 
                  i <= nvarmax_rad_zm )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zm ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                           "in the stats namelist, or change nvarmax_rad_zm."
          write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
          call endrun ("stats_init_clubb:  number of rad_zm statistical variables exceeds limit")
       endif

       stats_rad_zm%num_output_fields = ntot
       stats_rad_zm%kk = nnrad_zm

       allocate( stats_rad_zm%z( stats_rad_zm%kk ) )

       allocate( stats_rad_zm%accum_field_values( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%accum_num_samples( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%l_in_update( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )

       call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                     stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )

       allocate( stats_rad_zm%file%grid_avg_var( stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%file%z( stats_rad_zm%kk ) )
   
       call stats_init_rad_zm_api( clubb_vars_rad_zm, l_error, &
                                   stats_rad_zm )
    end if ! l_output_rad_files


    !  Initialize sfc (surface point)

    i = 1
    do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0 .and. & 
               len_trim(clubb_vars_sfc(i))   /= 0 .and. & 
               i <= nvarmax_sfc )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_sfc than allowed for by nvarmax_sfc."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                        "in the stats namelist, or change nvarmax_sfc."
       write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
       call endrun ("stats_init_clubb:  number of sfc statistical variables exceeds limit")
    endif

    stats_sfc%num_output_fields = ntot
    stats_sfc%kk = 1

    allocate( stats_sfc%z( stats_sfc%kk ) )

    allocate( stats_sfc%accum_field_values( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%accum_num_samples( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%l_in_update( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )

    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    allocate( stats_sfc%file%grid_avg_var( stats_sfc%num_output_fields ) )
    allocate( stats_sfc%file%z( stats_sfc%kk ) )

    if (first_call) then
      call stats_init_sfc_api( clubb_vars_sfc, l_error, &
                               stats_sfc )
    end if

    ! Check for errors

    if ( l_error ) then
       call endrun ('stats_init:  errors found')
    endif

!   Now call add fields
    if (first_call) then
      
      do i = 1, stats_zt%num_output_fields
      
        temp1 = trim(stats_zt%file%grid_avg_var(i)%name)
        sub   = temp1
        if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
       
          call addfld(trim(sub),(/ 'ilev' /),&
               'A',trim(stats_zt%file%grid_avg_var(i)%units),trim(stats_zt%file%grid_avg_var(i)%description))
      enddo
      
      do i = 1, stats_zm%num_output_fields
      
        temp1 = trim(stats_zm%file%grid_avg_var(i)%name)
        sub   = temp1
        if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
      
         call addfld(trim(sub),(/ 'ilev' /),&
              'A',trim(stats_zm%file%grid_avg_var(i)%units),trim(stats_zm%file%grid_avg_var(i)%description))
      enddo

      if (l_output_rad_files) then     

         do i = 1, stats_rad_zt%num_output_fields
            call addfld(trim(stats_rad_zt%file%grid_avg_var(i)%name),(/ 'ilev' /),&
               'A',trim(stats_rad_zt%file%grid_avg_var(i)%units),trim(stats_rad_zt%file%grid_avg_var(i)%description))
         enddo
      
         do i = 1, stats_rad_zm%num_output_fields
            call addfld(trim(stats_rad_zm%file%grid_avg_var(i)%name),(/ 'ilev' /),&
               'A',trim(stats_rad_zm%file%grid_avg_var(i)%units),trim(stats_rad_zm%file%grid_avg_var(i)%description))
         enddo
      endif
      
      do i = 1, stats_sfc%num_output_fields
         call addfld(trim(stats_sfc%file%grid_avg_var(i)%name),horiz_only,&
              'A',trim(stats_sfc%file%grid_avg_var(i)%units),trim(stats_sfc%file%grid_avg_var(i)%description))
      enddo
      
    end if

    return

  end subroutine stats_init_clubb 
  
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS  
  subroutine stats_end_timestep_clubb(thecol, stats_zt, stats_zm, stats_rad_zt, stats_rad_zm, stats_sfc, &
                                      out_zt, out_zm, out_radzt, out_radzm, out_sfc)
    !-----------------------------------------------------------------------
    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------



    use shr_infnan_mod, only: is_nan => shr_infnan_isnan

    use clubb_api_module, only: &
        fstderr, & ! Constant(s)
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout, &
        l_output_rad_files, &
        clubb_at_least_debug_level_api ! Procedure(s)

    use cam_abortutils,  only: endrun

    implicit none

    integer :: thecol
    
    ! Input Variables
    type (stats), intent(inout) :: stats_zt,      & ! stats_zt grid
                                   stats_zm,      & ! stats_zm grid
                                   stats_rad_zt,  & ! stats_rad_zt grid
                                   stats_rad_zm,  & ! stats_rad_zm grid
                                   stats_sfc        ! stats_sfc
    
    ! Inout variables
    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pverp,stats_zt%num_output_fields)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,stats_zt%num_output_fields)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pverp,stats_rad_zt%num_output_fields)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%num_output_fields)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%num_output_fields)

    ! Local Variables

    integer :: i, k
    logical :: l_error

    !  Check if it is time to write to file

    if ( .not. l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Compute averages
    call stats_avg( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, stats_zt%accum_num_samples )
    call stats_avg( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, stats_zm%accum_num_samples )
    if (l_output_rad_files) then
      call stats_avg( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                      stats_rad_zt%accum_num_samples )
      call stats_avg( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                      stats_rad_zm%accum_num_samples )
    end if
    call stats_avg( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, stats_sfc%accum_num_samples )

   !  Here we are not outputting the data, rather reading the stats into 
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.    
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk 
         out_zt(thecol,pverp-k+1,i) = stats_zt%accum_field_values(1,1,k,i)
         if(is_nan(out_zt(thecol,k,i))) out_zt(thecol,k,i) = 0.0_r8
      enddo   
    enddo

    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zt%kk 
         out_zm(thecol,pverp-k+1,i) = stats_zm%accum_field_values(1,1,k,i)
         if(is_nan(out_zm(thecol,k,i))) out_zm(thecol,k,i) = 0.0_r8
      enddo   
    enddo

    if (l_output_rad_files) then 
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk 
          out_radzt(thecol,pverp-k+1,i) = stats_rad_zt%accum_field_values(1,1,k,i)
          if(is_nan(out_radzt(thecol,k,i))) out_radzt(thecol,k,i) = 0.0_r8
        enddo   
      enddo
    
      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk 
          out_radzm(thecol,pverp-k+1,i) = stats_rad_zm%accum_field_values(1,1,k,i)
          if(is_nan(out_radzm(thecol,k,i))) out_radzm(thecol,k,i) = 0.0_r8
        enddo   
      enddo

      ! Fill in values above the CLUBB top.
      out_zt(thecol,:top_lev-1,:) = 0.0_r8
      out_zm(thecol,:top_lev-1,:) = 0.0_r8
      out_radzt(thecol,:top_lev-1,:) = 0.0_r8
      out_radzm(thecol,:top_lev-1,:) = 0.0_r8

    endif ! l_output_rad_files
    
    do i = 1, stats_sfc%num_output_fields
      out_sfc(thecol,1,i) = stats_sfc%accum_field_values(1,1,1,i)   
      if(is_nan(out_sfc(thecol,1,i))) out_sfc(thecol,1,i) = 0.0_r8
    enddo

    !  Reset sample fields
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                       stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )
      call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                       stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )
    end if
    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    return

  end subroutine stats_end_timestep_clubb
#endif  
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS
  
    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, num_output_fields, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_api_module, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, num_output_fields

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,num_output_fields), intent(out)    :: x
    integer(kind=stat_nknd), dimension(1,1,kk,num_output_fields), intent(out) :: n
    logical, dimension(1,1,kk,num_output_fields), intent(out)                 :: l_in_update

    !  Zero out arrays

    if ( num_output_fields > 0 ) then
       x(:,:,:,:) = 0.0_r8
       n(:,:,:,:) = 0
       l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero
  
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  
#ifdef CLUBB_SGS
    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, num_output_fields, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_api_module, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: num_output_fields, kk
    integer(kind=stat_nknd), dimension(1,1,kk,num_output_fields), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,num_output_fields), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m=1,num_output_fields
       do k=1,kk

          if ( n(1,1,k,m) > 0 ) then
             x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
          end if

       end do
    end do

    return

  end subroutine stats_avg

  subroutine grid_size(state, grid_dx, grid_dy)
  ! Determine the size of the grid for each of the columns in state

  use phys_grid,       only: get_area_p
  use shr_const_mod,   only: shr_const_pi
  use physics_types,   only: physics_state

   
  type(physics_state), intent(in) :: state
  real(r8), intent(out)           :: grid_dx(state%ncol), grid_dy(state%ncol)   ! CAM grid [m]

  real(r8), parameter :: earth_ellipsoid1 = 111132.92_r8 ! first coefficient, meters per degree longitude at equator
  real(r8), parameter :: earth_ellipsoid2 = 559.82_r8 ! second expansion coefficient for WGS84 ellipsoid
  real(r8), parameter :: earth_ellipsoid3 = 1.175_r8 ! third expansion coefficient for WGS84 ellipsoid

  real(r8) :: mpdeglat, column_area, degree
  integer  :: i

  ! determine the column area in radians
  do i=1,state%ncol
      column_area = get_area_p(state%lchnk,i)
      degree = sqrt(column_area)*(180._r8/shr_const_pi)
       
      ! Now find meters per degree latitude
      ! Below equation finds distance between two points on an ellipsoid, derived from expansion
      !  taking into account ellipsoid using World Geodetic System (WGS84) reference 
      mpdeglat = earth_ellipsoid1 - earth_ellipsoid2 * cos(2._r8*state%lat(i)) + earth_ellipsoid3 * cos(4._r8*state%lat(i))
      grid_dx(i) = mpdeglat * degree
      grid_dy(i) = grid_dx(i) ! Assume these are the same
  enddo   

  end subroutine grid_size     

#endif
  
end module clubb_intr
