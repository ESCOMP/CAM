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
  ! Modified by: K Thayer-Calder                                                                         !
  !                                                                                                      !
  !----------------------------------------------------------------------------------------------------- !

  use ref_pres,            only: trop_cloud_top_press

  use shr_kind_mod,        only: r8=>shr_kind_r8
  use ppgrid,              only: pver, pverp, pcols, begchunk, endchunk
  use phys_control,        only: phys_getopts
  use physconst,           only: cpair, gravit, rga, latvap, latice, zvir, rh2o, karman, pi, rair
  use air_composition,     only: rairv, cpairv
  use cam_history_support, only: max_fieldname_len

  use spmd_utils,          only: masterproc
  use constituents,        only: pcnst, cnst_add, cnst_ndropmixed
  use atmos_phys_pbl_utils,only: calc_friction_velocity, calc_kinematic_heat_flux, calc_ideal_gas_rrho, &
                                 calc_kinematic_water_vapor_flux, calc_kinematic_buoyancy_flux, calc_obukhov_length
  use ref_pres,            only: top_lev => trop_cloud_top_lev
  use scamMOD,             only: single_column, scm_clubb_iop_name, scm_cambfb_mode

#ifdef CLUBB_SGS
  use clubb_api_module,    only: pdf_parameter, implicit_coefs_terms, &
                                 clubb_config_flags_type, grid, stats, &
                                 nu_vertical_res_dep, stats_metadata_type, &
                                 hm_metadata_type, sclr_idx_type, &
                                 nparams

  use clubb_mf,            only: do_clubb_mf, do_clubb_mf_diag
  use cloud_fraction,      only: dp1, dp2
#endif

  implicit none

#ifdef CLUBB_SGS

#endif

  private

  save

  ! Subroutines to make public
  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, clubb_emissions_cam, &
            clubb_readnl, clubb_init_cnst, clubb_implements_cnst

#ifdef CLUBB_SGS
  
  ! NOTE: the only reason for anything in this section being set to public is for use with SILHS

  public :: stats_init_clubb, stats_end_timestep_clubb

  type(clubb_config_flags_type), public  :: &
    clubb_config_flags

  real(r8), dimension(1,nparams), public :: &
    clubb_params_single_col    ! Adjustable CLUBB parameters (C1, C2 ...)

  ! Variables that contains all the statistics
  type (stats), save, public :: &
    stats_zt(pcols),      & ! stats_zt grid
    stats_zm(pcols),      & ! stats_zm grid
    stats_rad_zt(pcols),  & ! stats_rad_zt grid
    stats_rad_zm(pcols),  & ! stats_rad_zm grid
    stats_sfc(pcols)        ! stats_sfc

  type (hm_metadata_type), public :: &
    hm_metadata

  type (stats_metadata_type), public :: &
    stats_metadata

  type (sclr_idx_type), public :: &
    sclr_idx

  logical, public, parameter :: &
    l_ascending_grid = .true. ! Set clubb to ascending mode, which is opposite of the 
                               ! cam grid the rest of this code uses, thus it requires
                               ! an expensive array flipping step before calling clubb

  integer, public :: &
    nzm_clubb,          & ! Number of vertical levels used by CLUBB momentum variables
    nzt_clubb             ! Number of vertical levels used by CLUBB thermodynamic variables

  ! These are zero by default, but will be set by SILHS before they are used by subcolumns
  integer, public :: &
    hydromet_dim = 0, &
    pdf_dim      = 0

  type(pdf_parameter), allocatable, public :: &
    pdf_params_chnk(:)                ! PDF parameters (thermo. levs.) [units vary]

  type(pdf_parameter), allocatable :: &
    pdf_params_zm_chnk(:)             ! PDF parameters on momentum levs. [units vary]

  type(implicit_coefs_terms), allocatable :: &
    pdf_implicit_coefs_terms_chnk(:)  ! PDF impl. coefs. & expl. terms      [units vary]

  real(r8), public :: &
    ztodt  ! model timestep
#endif

  ! ------------------------------------------------------------ !
  !                           CONSTANTS                          !
  ! ------------------------------------------------------------ !

  integer, parameter :: &
      grid_type    = 3, &               ! The 2 option specifies stretched thermodynamic levels
      sclr_dim     = 0                  ! Higher-order scalars, set to zero

  ! Even though sclr_dim is set to 0, the dimension here is set to 1 to prevent compiler errors
  ! See github ticket larson-group/cam#133 for details
  real(r8), parameter, dimension(1) :: &
      sclr_tol = 1.e-8_r8               ! Total water in kg/kg

  real(r8), parameter :: &
      rtm_min                 = epsilon( rtm_min ),   & ! Value below which rtm will be nudged [kg/kg]
      rtm_nudge_max_altitude  = 10000._r8,            & ! Highest altitude at which to nudge rtm [m]
      theta0                  = 300._r8,              & ! Reference temperature                     [K]
      ts_nudge                = 86400._r8,            & ! Time scale for u/v nudging (not used)     [s]
      p0_clubb                = 100000._r8

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

  logical, parameter, private :: &
    l_implemented    = .true.            ! Implemented in a host model (always true)

  ! ----------------------------------------------------------------- !
  !   Things shared between subroutines: generally because they are   !
  !   set by an initialization routine, then used by clubb_tend_cam   !
  ! ----------------------------------------------------------------- !

  logical :: do_cldcool
  logical :: clubb_do_icesuper

  logical :: &
    clubb_l_intr_sfc_flux_smooth = .false. ! Add a locally calculated roughness to upwp and vpwp sfc fluxes

  logical            :: lq(pcnst)
  logical            :: do_rainturb
  logical            :: clubb_do_adv
  logical            :: clubb_do_liqsupersat = .false.
  logical            :: clubb_do_energyfix   = .true.
  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB

  integer :: &
    ixthlp2 = 0, &
    ixwpthlp = 0, &
    ixwprtp = 0, &
    ixwp2 = 0, &
    ixwp3 = 0, &
    ixrtpthlp = 0, &
    ixrtp2 = 0, &
    ixup2 = 0, &
    ixvp2 = 0

  !  Output arrays for CLUBB statistics
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  ! Outputs from phys_getopts
  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90
  character(len=16)  :: deep_scheme      ! Default set in phys_control.F90
  logical            :: history_budget
  integer            :: history_budget_histfile_num
  logical            :: do_hb_above_clubb    = .false.

  character(len=16)  :: subcol_scheme

  ! For clubb_do_adv
  integer, parameter :: ncnst=9
  character(len=8)   :: cnst_names(ncnst)
  logical            :: do_cnst=.false.


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
  real(r8) :: clubb_bv_efold = unset_r8
  real(r8) :: clubb_wpxp_Ri_exp = unset_r8
  real(r8) :: clubb_z_displace = unset_r8

  integer :: &
    clubb_iiPDF_type,          & ! Selected option for the two-component normal
                                 ! (double Gaussian) PDF type to use for the w, rt,
                                 ! and theta-l (or w, chi, and eta) portion of
                                 ! CLUBB's multivariate, two-component PDF.
    clubb_ipdf_call_placement = unset_i, & ! Selected option for the placement of the call to
                                           ! CLUBB's PDF.
    clubb_penta_solve_method = unset_i,  & ! Specifier for method to solve the penta-diagonal system
    clubb_tridiag_solve_method = unset_i,& ! Specifier for method to solve tri-diagonal systems
    clubb_saturation_equation = unset_i, & ! Specifier for which saturation formula to use
    clubb_grid_remap_method = unset_i, & ! Specifier for which method should be used to
                                         ! map values from one grid to another
                                         ! (starts at 1, so 0 is an invalid option for this flag)
    clubb_grid_adapt_in_time_method = unset_i,     & ! Specifier for how the grid density method should
                                                     ! be constructed if the grid should be adapted over time
                                                     ! (set to 0 for no adaptation)
    clubb_fill_holes_type = unset_i                  ! Option for which type of hole filler to use in the 
                                                     ! fill_holes_vertical procedure


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
    clubb_l_use_wp3_lim_with_smth_Heaviside, & ! Flag to activate mods on wp3 limiters for conv test
    clubb_l_smooth_Heaviside_tau_wpxp,  & ! Use smooth Heaviside 'Peskin' in computation of invrs_tau
    clubb_l_modify_limiters_for_cnvg_test, & ! Flag to activate mods on limiters for conv test
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
    clubb_l_host_applies_sfc_fluxes,    & ! Whether the host model applies the surface fluxes
    clubb_l_wp2_fill_holes_tke,         & ! Whether TKE is taken from up2 and vp2 to fill holes in wp2
    clubb_l_add_dycore_grid               ! Flag to remap values from dycore grid

  ! ------------------------------------------------------------ !
  !             Indices for physics buffer (pbuf)                !
  ! ------------------------------------------------------------ !
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
    wp2rtp_idx, &       ! w'^2 rt'
    wp2thlp_idx, &      ! w'^2 thl'
    uprcp_idx, &        ! < u' r_c' >
    vprcp_idx, &        ! < v' r_c' >
    rc_coef_zm_idx, &   ! Coefficient of X'r_c' in Eq. (34)
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
    qsatfac_idx, &      ! subgrid cloud water saturation scaling factor
    ice_supersat_idx, & ! ice cloud fraction for SILHS
    rcm_idx, &          ! Cloud water mixing ratio for SILHS
    clubbtop_idx        ! level index for CLUBB top

  ! For Gravity Wave code
  integer :: &
    ttend_clubb_idx, &
    ttend_clubb_mc_idx, &
    upwp_clubb_gw_idx, &
    upwp_clubb_gw_mc_idx, &
    vpwp_clubb_gw_idx, &
    vpwp_clubb_gw_mc_idx, &
    thlp2_clubb_gw_idx, &
    thlp2_clubb_gw_mc_idx, &
    wpthlp_clubb_gw_idx, &
    wpthlp_clubb_gw_mc_idx

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

  integer :: &
    cmfmc_sh_idx = 0

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

    !----- Begin Code -----
    call phys_getopts( eddy_scheme_out                 = eddy_scheme, &
                       deep_scheme_out                 = deep_scheme, &
                       history_budget_out              = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       do_hb_above_clubb_out           = do_hb_above_clubb)

    subcol_scheme = subcol_get_scheme()

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
    
    ! Determine number of vertical levels used in clubb, thermo variables are nzt_clubb
    ! and momentum variables are nzm_clubb
    nzt_clubb = pver  + 1 - top_lev
    nzm_clubb = pverp + 1 - top_lev

    if (do_hb_above_clubb) then
      call pbuf_add_field('clubbtop', 'physpkg', dtype_i4, (/pcols/), clubbtop_idx)
    endif

    !  put pbuf_add calls here (see macrop_driver.F90 for sample) use indicies defined at top
    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/),                      pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/),               tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/),               kvh_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/),                      tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),   cld_idx)
    call pbuf_add_field('FICE',       'physpkg',dtype_r8, (/pcols,pver/),                 fice_idx)
    call pbuf_add_field('CMELIQ',     'physpkg',dtype_r8, (/pcols,pver/),                 cmeliq_idx)
    call pbuf_add_field('QSATFAC',    'physpkg',dtype_r8, (/pcols,pver/),                 qsatfac_idx)

    ! pbuf fields for Gravity Wave scheme
    call pbuf_add_field('TTEND_CLUBB',     'physpkg', dtype_r8, (/pcols,pver/),      ttend_clubb_idx )
    call pbuf_add_field('UPWP_CLUBB_GW',   'physpkg', dtype_r8, (/pcols,pverp/),   upwp_clubb_gw_idx )
    call pbuf_add_field('VPWP_CLUBB_GW',   'physpkg', dtype_r8, (/pcols,pverp/),   vpwp_clubb_gw_idx )
    call pbuf_add_field('THLP2_CLUBB_GW',  'physpkg', dtype_r8, (/pcols,pverp/),  thlp2_clubb_gw_idx )
    call pbuf_add_field('WPTHLP_CLUBB_GW', 'physpkg', dtype_r8, (/pcols,pverp/), wpthlp_clubb_gw_idx )


    ! For SILHS microphysical covariance contributions
    call pbuf_add_field('rtp2_mc_zt',     'global', dtype_r8, (/pcols,nzt_clubb/), rtp2_mc_zt_idx)
    call pbuf_add_field('thlp2_mc_zt',    'global', dtype_r8, (/pcols,nzt_clubb/), thlp2_mc_zt_idx)
    call pbuf_add_field('wprtp_mc_zt',    'global', dtype_r8, (/pcols,nzt_clubb/), wprtp_mc_zt_idx)
    call pbuf_add_field('wpthlp_mc_zt',   'global', dtype_r8, (/pcols,nzt_clubb/), wpthlp_mc_zt_idx)
    call pbuf_add_field('rtpthlp_mc_zt',  'global', dtype_r8, (/pcols,nzt_clubb/), rtpthlp_mc_zt_idx)


    ! Only in clubb_intr.F90, these are safe to dimensions (ngrdcol,nzm_clubb) or (ngrdcol,nzt_clubb)
    call pbuf_add_field('pdf_zm_w_1',         'global', dtype_r8, (/pcols,nzm_clubb/), pdf_zm_w_1_idx)
    call pbuf_add_field('pdf_zm_w_2',         'global', dtype_r8, (/pcols,nzm_clubb/), pdf_zm_w_2_idx)
    call pbuf_add_field('pdf_zm_var_w_1',     'global', dtype_r8, (/pcols,nzm_clubb/), pdf_zm_varnce_w_1_idx)
    call pbuf_add_field('pdf_zm_var_w_2',     'global', dtype_r8, (/pcols,nzm_clubb/), pdf_zm_varnce_w_2_idx)
    call pbuf_add_field('pdf_zm_mixt_frac',   'global', dtype_r8, (/pcols,nzm_clubb/), pdf_zm_mixt_frac_idx)

    call pbuf_add_field('WPTHVP',     'global', dtype_r8, (/pcols,nzm_clubb/), wpthvp_idx)
    call pbuf_add_field('RTPTHVP',    'global', dtype_r8, (/pcols,nzm_clubb/), rtpthvp_idx)
    call pbuf_add_field('THLPTHVP',   'global', dtype_r8, (/pcols,nzm_clubb/), thlpthvp_idx)
    call pbuf_add_field('UPRCP',      'global', dtype_r8, (/pcols,nzm_clubb/), uprcp_idx)
    call pbuf_add_field('VPRCP',      'global', dtype_r8, (/pcols,nzm_clubb/), vprcp_idx)
    call pbuf_add_field('RC_COEF_ZM', 'global', dtype_r8, (/pcols,nzm_clubb/), rc_coef_zm_idx)
    call pbuf_add_field('WP4',        'global', dtype_r8, (/pcols,nzm_clubb/), wp4_idx)
    call pbuf_add_field('WP2UP2',     'global', dtype_r8, (/pcols,nzm_clubb/), wp2up2_idx)
    call pbuf_add_field('WP2VP2',     'global', dtype_r8, (/pcols,nzm_clubb/), wp2vp2_idx)

    call pbuf_add_field('UPWP',            'global', dtype_r8, (/pcols,nzm_clubb/), upwp_idx)
    call pbuf_add_field('VPWP',            'global', dtype_r8, (/pcols,nzm_clubb/), vpwp_idx)
    call pbuf_add_field('WPTHLP_nadv',     'global', dtype_r8, (/pcols,nzm_clubb/), wpthlp_idx)
    call pbuf_add_field('WPRTP_nadv',      'global', dtype_r8, (/pcols,nzm_clubb/), wprtp_idx)
    call pbuf_add_field('RTPTHLP_nadv',    'global', dtype_r8, (/pcols,nzm_clubb/), rtpthlp_idx)
    call pbuf_add_field('RTP2_nadv',       'global', dtype_r8, (/pcols,nzm_clubb/), rtp2_idx)
    call pbuf_add_field('THLP2_nadv',      'global', dtype_r8, (/pcols,nzm_clubb/), thlp2_idx)

    call pbuf_add_field('TTEND_CLUBB_MC',     'physpkg', dtype_r8, (/pcols,nzt_clubb/), ttend_clubb_mc_idx)
    call pbuf_add_field('UPWP_CLUBB_GW_MC',   'physpkg', dtype_r8, (/pcols,nzm_clubb/), upwp_clubb_gw_mc_idx)
    call pbuf_add_field('VPWP_CLUBB_GW_MC',   'physpkg', dtype_r8, (/pcols,nzm_clubb/), vpwp_clubb_gw_mc_idx)
    call pbuf_add_field('THLP2_CLUBB_GW_MC',  'physpkg', dtype_r8, (/pcols,nzm_clubb/), thlp2_clubb_gw_mc_idx)
    call pbuf_add_field('WPTHLP_CLUBB_GW_MC', 'physpkg', dtype_r8, (/pcols,nzm_clubb/), wpthlp_clubb_gw_mc_idx)

    call pbuf_add_field('WP2THVP',    'global',  dtype_r8, (/pcols,nzt_clubb/), wp2thvp_idx)
    call pbuf_add_field('RCM',        'physpkg', dtype_r8, (/pcols,nzt_clubb/), rcm_idx)
    call pbuf_add_field('THLM',       'global',  dtype_r8, (/pcols,nzt_clubb/), thlm_idx)
    call pbuf_add_field('WP2RTP',     'global',  dtype_r8, (/pcols,nzt_clubb/), wp2rtp_idx)
    call pbuf_add_field('WP2THLP',    'global',  dtype_r8, (/pcols,nzt_clubb/), wp2thlp_idx)
    call pbuf_add_field('WPUP2',      'global',  dtype_r8, (/pcols,nzt_clubb/), wpup2_idx)
    call pbuf_add_field('WPVP2',      'global',  dtype_r8, (/pcols,nzt_clubb/), wpvp2_idx)

    call pbuf_add_field('UM',         'global', dtype_r8, (/pcols,nzt_clubb/), um_idx)
    call pbuf_add_field('VM',         'global', dtype_r8, (/pcols,nzt_clubb/), vm_idx)
    call pbuf_add_field('RTP3',       'global', dtype_r8, (/pcols,nzt_clubb/), rtp3_idx)
    call pbuf_add_field('THLP3',      'global', dtype_r8, (/pcols,nzt_clubb/), thlp3_idx)
    call pbuf_add_field('UP3',        'global', dtype_r8, (/pcols,nzt_clubb/), up3_idx)
    call pbuf_add_field('VP3',        'global', dtype_r8, (/pcols,nzt_clubb/), vp3_idx)
    call pbuf_add_field('WP3_nadv',   'global', dtype_r8, (/pcols,nzt_clubb/), wp3_idx)

    call pbuf_add_field('UP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), up2_idx)
    call pbuf_add_field('VP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), vp2_idx)
    call pbuf_add_field('WP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), wp2_idx)

    call pbuf_add_field('CLOUD_FRAC', 'global', dtype_r8, (/pcols,pver/), cloud_frac_idx)

    ! Only in clubb_intr.F90 or SILHS
    call pbuf_add_field('ISS_FRAC',   'global', dtype_r8, (/pcols,nzt_clubb/), ice_supersat_idx)
    call pbuf_add_field('RTM',        'global', dtype_r8, (/pcols,nzt_clubb/), rtm_idx)

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
      initialize_clubb_config_flags_type_api
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
         clubb_bv_efold, &
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
         clubb_grid_adapt_in_time_method, &
         clubb_fill_holes_type, &
         clubb_grid_remap_method, &
         clubb_iiPDF_type, &
         clubb_ipdf_call_placement, &
         clubb_lambda0_stability_coef, &
         clubb_lmin_coef, &
         clubb_l_add_dycore_grid, &
         clubb_l_brunt_vaisala_freq_moist, &
         clubb_l_C2_cloud_frac, &
         clubb_l_calc_thlp2_rad, &
         clubb_l_calc_w_corr, &
         clubb_l_call_pdf_closure_twice, &
         clubb_l_const_Nc_in_cloud, &
         clubb_l_damp_wp2_using_em, &
         clubb_l_damp_wp3_Skw_squared, &
         clubb_l_diag_Lscale_from_tau, &
         clubb_l_diagnose_correlations, &
         clubb_l_diffuse_rtm_and_thlm, &
         clubb_l_do_expldiff_rtm_thlm, &
         clubb_l_e3sm_config, &
         clubb_l_enable_relaxed_clipping, &
         clubb_l_fix_w_chi_eta_correlations, &
         clubb_l_godunov_upwind_wpxp_ta, &
         clubb_l_godunov_upwind_xpyp_ta, &
         clubb_l_intr_sfc_flux_smooth, &
         clubb_l_lmm_stepping, &
         clubb_l_lscale_plume_centered, &
         clubb_l_min_wp2_from_corr_wx, &
         clubb_l_min_xp2_from_corr_wx, &
         clubb_l_modify_limiters_for_cnvg_test, &
         clubb_l_mono_flux_lim_rtm, &
         clubb_l_mono_flux_lim_spikefix, &
         clubb_l_mono_flux_lim_thlm, &
         clubb_l_mono_flux_lim_um, &
         clubb_l_mono_flux_lim_vm, &
         clubb_l_partial_upwind_wp3, &
         clubb_l_predict_upwp_vpwp, &
         clubb_l_prescribed_avg_deltaz, &
         clubb_l_rcm_supersat_adj, &
         clubb_l_rtm_nudge, &
         clubb_l_smooth_Heaviside_tau_wpxp, &
         clubb_l_stability_correct_Kh_N2_zm, &
         clubb_l_stability_correct_tau_zm, &
         clubb_l_standard_term_ta, &
         clubb_l_tke_aniso, &
         clubb_l_trapezoidal_rule_zm, &
         clubb_l_trapezoidal_rule_zt, &
         clubb_l_upwind_xm_ma, &
         clubb_l_upwind_xpyp_ta, &
         clubb_l_use_C11_Richardson, &
         clubb_l_use_C7_Richardson, &
         clubb_l_use_cloud_cover, &
         clubb_l_use_precip_frac, &
         clubb_l_use_shear_Richardson, &
         clubb_l_use_thvm_in_bv_freq, &
         clubb_l_use_tke_in_wp2_wp3_K_dfsn, &
         clubb_l_use_tke_in_wp3_pr_turb_term, &
         clubb_l_use_wp3_lim_with_smth_Heaviside, &
         clubb_l_uv_nudge, &
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
         clubb_wpxp_L_thresh, &
         clubb_wpxp_Ri_exp, &
         clubb_z_displace

    !----- Begin Code -----

    !  Determine if we want clubb_history to be output
    clubb_history                     = .false.   ! Initialize to false
    stats_metadata%l_stats            = .false.   ! Initialize to false
    stats_metadata%l_output_rad_files = .false.   ! Initialize to false
    do_cldcool                        = .false.   ! Initialize to false
    do_rainturb                       = .false.   ! Initialize to false

    ! Initialize namelist variables to clubb defaults
    call set_default_clubb_config_flags_api( clubb_iiPDF_type, & ! Out
                                             clubb_ipdf_call_placement, & ! Out
                                             clubb_penta_solve_method, & ! Out
                                             clubb_tridiag_solve_method, & ! Out
                                             clubb_saturation_equation, & ! Out
                                             clubb_grid_remap_method, & ! Out
                                             clubb_grid_adapt_in_time_method, & ! Out
                                             clubb_fill_holes_type, & ! Out
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
                                             clubb_l_use_wp3_lim_with_smth_Heaviside, & ! Out
                                             clubb_l_smooth_Heaviside_tau_wpxp, & ! Out
                                             clubb_l_modify_limiters_for_cnvg_test, & ! Out
                                             clubb_l_enable_relaxed_clipping, & ! Out
                                             clubb_l_linearize_pbl_winds, & ! Out
                                             clubb_l_mono_flux_lim_thlm, & ! Out
                                             clubb_l_mono_flux_lim_rtm, & ! Out
                                             clubb_l_mono_flux_lim_um, & ! Out
                                             clubb_l_mono_flux_lim_vm, & ! Out
                                             clubb_l_mono_flux_lim_spikefix, &  ! Out
                                             clubb_l_host_applies_sfc_fluxes, & ! Out
                                             clubb_l_wp2_fill_holes_tke, & ! Out
                                             clubb_l_add_dycore_grid ) ! Out

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
    call mpi_bcast(clubb_bv_efold,         1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_bv_efold")
    call mpi_bcast(clubb_wpxp_Ri_exp,      1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_wpxp_Ri_exp")
    call mpi_bcast(clubb_z_displace,       1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_z_displace")
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
    call mpi_bcast(clubb_iiPDF_type, 1, mpi_integer,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_iiPDF_type")

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
    call mpi_bcast(clubb_l_use_wp3_lim_with_smth_Heaviside, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_wp3_lim_with_smth_Heaviside")
    call mpi_bcast(clubb_l_smooth_Heaviside_tau_wpxp,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_smooth_Heaviside_tau_wpxp")
    call mpi_bcast(clubb_l_modify_limiters_for_cnvg_test, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_modify_limiters_for_cnvg_test")
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
    call mpi_bcast(clubb_l_host_applies_sfc_fluxes,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_host_applies_sfc_fluxes")
    call mpi_bcast(clubb_l_wp2_fill_holes_tke,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_wp2_fill_holes_tke")
    call mpi_bcast(clubb_l_add_dycore_grid,   1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_add_dycore_grid")
    call mpi_bcast(clubb_penta_solve_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_penta_solve_method")
    call mpi_bcast(clubb_tridiag_solve_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_tridiag_solve_method")
    call mpi_bcast(clubb_saturation_equation,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_saturation_equation")
    call mpi_bcast(clubb_grid_remap_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_grid_remap_method")
    call mpi_bcast(clubb_grid_adapt_in_time_method,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_grid_adapt_in_time_method")
    call mpi_bcast(clubb_fill_holes_type,    1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_fill_holes_type")
    call mpi_bcast(clubb_l_intr_sfc_flux_smooth,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_intr_sfc_flux_smooth")
    call mpi_bcast(clubb_l_vary_convect_depth,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_vary_convect_depth")
    call mpi_bcast(clubb_l_standard_term_ta,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_standard_term_ta")
    call mpi_bcast(clubb_l_partial_upwind_wp3,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_partial_upwind_wp3")
    call mpi_bcast(clubb_l_C2_cloud_frac,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_C2_cloud_frac")
    call mpi_bcast(clubb_l_calc_thlp2_rad,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_calc_thlp2_rad")
    call mpi_bcast(clubb_l_calc_w_corr,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_calc_w_corr")
    call mpi_bcast(clubb_l_const_Nc_in_cloud,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_const_Nc_in_cloud")
    call mpi_bcast(clubb_l_diagnose_correlations,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_diagnose_correlations")
    call mpi_bcast(clubb_l_diffuse_rtm_and_thlm,    1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_diffuse_rtm_and_thlm")
    call mpi_bcast(clubb_l_fix_w_chi_eta_correlations, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_fix_w_chi_eta_correlations")
    call mpi_bcast(clubb_l_prescribed_avg_deltaz, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_prescribed_avg_deltaz")
    call mpi_bcast(clubb_l_rtm_nudge, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_rtm_nudge")
    call mpi_bcast(clubb_l_stability_correct_Kh_N2_zm, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_stability_correct_Kh_N2_zm")
    call mpi_bcast(clubb_l_tke_aniso, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_tke_aniso")
    call mpi_bcast(clubb_l_upwind_xm_ma, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_upwind_xm_ma")
    call mpi_bcast(clubb_l_use_precip_frac, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_use_precip_frac")
    call mpi_bcast(clubb_l_uv_nudge, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_uv_nudge")

    !  Overwrite defaults if they are true
    if (clubb_history) stats_metadata%l_stats = .true.
    if (clubb_rad_history) stats_metadata%l_output_rad_files = .true.
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
    if(clubb_C_wp2_splat == unset_r8) call endrun(sub//": FATAL: clubb_C_wp2_splat is not set")
    if(clubb_bv_efold == unset_r8) call endrun(sub//": FATAL: clubb_bv_efold is not set")
    if(clubb_wpxp_Ri_exp == unset_r8) call endrun(sub//": FATAL: clubb_wpxp_Ri_exp is not set")
    if(clubb_z_displace == unset_r8) call endrun(sub//": FATAL: clubb_z_displace is not set")
    if(clubb_detliq_rad == unset_r8) call endrun(sub//": FATAL: clubb_detliq_rad not set")
    if(clubb_detice_rad == unset_r8) call endrun(sub//": FATAL: clubb_detice_rad not set")
    if(clubb_ipdf_call_placement == unset_i) call endrun(sub//": FATAL: clubb_ipdf_call_placement not set")
    if(clubb_detphase_lowtemp == unset_r8) call endrun(sub//": FATAL: clubb_detphase_lowtemp not set")
    if(clubb_penta_solve_method == unset_i) call endrun(sub//": FATAL: clubb_penta_solve_method not set")
    if(clubb_tridiag_solve_method == unset_i) call endrun(sub//": FATAL: clubb_tridiag_solve_method not set")
    if(clubb_saturation_equation == unset_i) call endrun(sub//": FATAL: clubb_saturation_equation not set")
    if(clubb_grid_remap_method == unset_i) call endrun(sub//": FATAL: clubb_grid_remap_method not set")
    if(clubb_grid_adapt_in_time_method == unset_i) call endrun(sub//": FATAL: clubb_grid_adapt_in_time_method not set")
    if(clubb_fill_holes_type == unset_i) call endrun(sub//": FATAL: clubb_fill_holes_type not set")
    if(clubb_detphase_lowtemp >= meltpt_temp) &
    call endrun(sub//": ERROR: clubb_detphase_lowtemp must be less than 268.15 K")

    call initialize_clubb_config_flags_type_api( clubb_iiPDF_type, & ! In
                                                 clubb_ipdf_call_placement, & ! In
                                                 clubb_penta_solve_method, & ! In
                                                 clubb_tridiag_solve_method, & ! In
                                                 clubb_saturation_equation, & ! In
                                                 clubb_grid_remap_method, & ! In
                                                 clubb_grid_adapt_in_time_method, & ! In
                                                 clubb_fill_holes_type, & ! In
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
                                                 clubb_l_use_wp3_lim_with_smth_Heaviside, & ! In
                                                 clubb_l_smooth_Heaviside_tau_wpxp, & ! In
                                                 clubb_l_modify_limiters_for_cnvg_test, & ! In
                                                 clubb_l_enable_relaxed_clipping, & ! In
                                                 clubb_l_linearize_pbl_winds, & ! In
                                                 clubb_l_mono_flux_lim_thlm, & ! In
                                                 clubb_l_mono_flux_lim_rtm, & ! In
                                                 clubb_l_mono_flux_lim_um, & ! In
                                                 clubb_l_mono_flux_lim_vm, & ! In
                                                 clubb_l_mono_flux_lim_spikefix, & ! In
                                                 clubb_l_host_applies_sfc_fluxes, & ! In
                                                 clubb_l_wp2_fill_holes_tke, & ! In
                                                 clubb_l_add_dycore_grid, & ! In
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
         iC_invrs_tau_N2_xp2,iC_invrs_tau_N2_wpxp,iC_invrs_tau_N2_clear_wp3, &
         iC2rt, iC2thl, iC2rtthl, ic_K1, ic_K2, inu2, ic_K8, ic_K9, inu9, iC_wp2_splat, ibv_efold, &
         iwpxp_Ri_exp, iz_displace, &
         params_list

    use clubb_api_module, only: &
         print_clubb_config_flags_api, &
         check_clubb_settings_api, &
         init_pdf_params_api, &
         time_precision, &
         core_rknd, &
         set_clubb_debug_level_api, &
         clubb_fatal_error, &     ! Error code value to indicate a fatal error
         err_info_type, &
         init_default_err_info_api, &
         cleanup_err_info_api, &
         nparams, &
         init_clubb_params_api, &
         w_tol_sqd, &
         rt_tol, &
         thl_tol, &
         saturation_bolton, & ! Constant for Bolton approximations of saturation
         saturation_gfdl,   & ! Constant for the GFDL approximation of saturation
         saturation_flatau, & ! Constant for Flatau approximations of saturation
         saturation_lookup    ! Use a lookup table for mixing length

    use time_manager,              only: is_first_step
    use constituents,           only: cnst_get_ind
    use phys_control,           only: phys_getopts
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

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header
      
    integer :: i, j, k, l                    ! Indices
    integer :: nmodes, nspec, m
    integer :: ixq, ixcldice, ixcldliq, ixnumliq, ixnumice
    integer :: lptr

    logical, parameter :: l_input_fields = .false. ! Always false for CAM-CLUBB.
    logical, parameter :: l_update_pressure = .false. ! Always false for CAM-CLUBB.

    integer :: ierr=0

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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, wpxp_Ri_exp, &
      a3_coef_min, a_const, bv_efold, z_displace

    !----- Begin Code -----

    if (core_rknd /= r8) then
      call endrun('clubb_ini_cam:  CLUBB library core_rknd must match CAM r8 and it does not')
    end if

    ! Allocate PDF parameters across columns and chunks
    allocate( &
       pdf_params_chnk(begchunk:endchunk),   &
       pdf_params_zm_chnk(begchunk:endchunk), &
       pdf_implicit_coefs_terms_chnk(begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) call endrun(' clubb_ini_cam: failed to allocate pdf_params')

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that
    ! CLUBB does not transport aerosol consituents.  Therefore, need to
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents)
    ! ----------------------------------------------------------------- !

    call phys_getopts(history_amwg_out=history_amwg, &
                      history_clubb_out=history_clubb, &
                      do_hb_above_clubb_out=do_hb_above_clubb)

    !  Select variables to apply tendencies back to CAM

    ! Initialize all consituents to true to start
    lq(1:pcnst) = .true.
    edsclr_dim  = pcnst

    call cnst_get_ind('Q',ixq)
    call cnst_get_ind('NUMICE',ixnumice)
    call cnst_get_ind('NUMLIQ',ixnumliq)
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)

    do m = 1, pcnst
       if (cnst_ndropmixed(m)) then
          lq(m)=.false.
          !  Droplet number is transported in dropmixnuc, therefore we
          !  do NOT want CLUBB to apply transport tendencies to avoid double
          !  counting.  Else, we apply tendencies.
          edsclr_dim = edsclr_dim-1
       endif
    enddo

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
    stats_metadata%l_stats_samp = .false.
    stats_metadata%l_grads = .false.

    !  Overwrite defaults if needed
    if (stats_metadata%l_stats) stats_metadata%l_stats_samp = .true.

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
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')
    naai_idx        = pbuf_get_index('NAAI')
    npccn_idx       = pbuf_get_index('NPCCN')


    sclr_idx%iisclr_rt  = -1
    sclr_idx%iisclr_thl = -1
    sclr_idx%iisclr_CO2 = -1

    sclr_idx%iiedsclr_rt  = -1
    sclr_idx%iiedsclr_thl = -1
    sclr_idx%iiedsclr_CO2 = -1

    ! ----------------------------------------------------------------- !
    ! Define number of tracers for CLUBB to diffuse
    ! ----------------------------------------------------------------- !

    if (clubb_l_do_expldiff_rtm_thlm) then
      ! add 2 since we want to diffuse temperature and moisture explicitly as well
      edsclr_dim = edsclr_dim + 2
    endif

    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !

     call init_clubb_params_api( 1, -99, "", &
                                     clubb_params_single_col )

    clubb_params_single_col(:,iC2rtthl) = clubb_C2rtthl
    clubb_params_single_col(:,iC8) = clubb_C8
    clubb_params_single_col(:,iC11) = clubb_c11
    clubb_params_single_col(:,iC11b) = clubb_c11b
    clubb_params_single_col(:,iC14) = clubb_c14
    clubb_params_single_col(:,iC_wp3_pr_turb) = clubb_C_wp3_pr_turb
    clubb_params_single_col(:,ic_K10) = clubb_c_K10
    clubb_params_single_col(:,imult_coef) = clubb_mult_coef
    clubb_params_single_col(:,iSkw_denom_coef) = clubb_Skw_denom_coef
    clubb_params_single_col(:,iC2rt) = clubb_C2rt
    clubb_params_single_col(:,iC2thl) = clubb_C2thl
    clubb_params_single_col(:,ibeta) = clubb_beta
    clubb_params_single_col(:,iC6rt) = clubb_c6rt
    clubb_params_single_col(:,iC6rtb) = clubb_c6rtb
    clubb_params_single_col(:,iC6rtc) = clubb_c6rtc
    clubb_params_single_col(:,iC6thl) = clubb_c6thl
    clubb_params_single_col(:,iC6thlb) = clubb_c6thlb
    clubb_params_single_col(:,iC6thlc) = clubb_c6thlc
    clubb_params_single_col(:,iwpxp_L_thresh) = clubb_wpxp_L_thresh
    clubb_params_single_col(:,iC7) = clubb_C7
    clubb_params_single_col(:,iC7b) = clubb_C7b
    clubb_params_single_col(:,igamma_coef) = clubb_gamma_coef
    clubb_params_single_col(:,ic_K10h) = clubb_c_K10h
    clubb_params_single_col(:,ilambda0_stability_coef) = clubb_lambda0_stability_coef
    clubb_params_single_col(:,ilmin_coef) = clubb_lmin_coef
    clubb_params_single_col(:,iC8b) = clubb_C8b
    clubb_params_single_col(:,iskw_max_mag) = clubb_skw_max_mag
    clubb_params_single_col(:,iC1)  = clubb_C1
    clubb_params_single_col(:,iC1b) = clubb_C1b
    clubb_params_single_col(:,igamma_coefb) = clubb_gamma_coefb
    clubb_params_single_col(:,iup2_sfc_coef) = clubb_up2_sfc_coef
    clubb_params_single_col(:,iC4) = clubb_C4
    clubb_params_single_col(:,iC_uu_shr) = clubb_C_uu_shr
    clubb_params_single_col(:,iC_uu_buoy) = clubb_C_uu_buoy
    clubb_params_single_col(:,ic_K1) = clubb_c_K1
    clubb_params_single_col(:,ic_K2) = clubb_c_K2
    clubb_params_single_col(:,inu2)  = clubb_nu2
    clubb_params_single_col(:,ic_K8) = clubb_c_K8
    clubb_params_single_col(:,ic_K9) = clubb_c_K9
    clubb_params_single_col(:,inu9)  = clubb_nu9
    clubb_params_single_col(:,iC_wp2_splat) = clubb_C_wp2_splat
    clubb_params_single_col(:,iC_invrs_tau_bkgnd) = clubb_C_invrs_tau_bkgnd
    clubb_params_single_col(:,iC_invrs_tau_sfc) = clubb_C_invrs_tau_sfc
    clubb_params_single_col(:,iC_invrs_tau_shear) = clubb_C_invrs_tau_shear
    clubb_params_single_col(:,iC_invrs_tau_N2) = clubb_C_invrs_tau_N2
    clubb_params_single_col(:,iC_invrs_tau_N2_wp2) = clubb_C_invrs_tau_N2_wp2
    clubb_params_single_col(:,iC_invrs_tau_N2_xp2) = clubb_C_invrs_tau_N2_xp2
    clubb_params_single_col(:,iC_invrs_tau_N2_wpxp) = clubb_C_invrs_tau_N2_wpxp
    clubb_params_single_col(:,iC_invrs_tau_N2_clear_wp3) = clubb_C_invrs_tau_N2_clear_wp3
    clubb_params_single_col(:,ibv_efold) = clubb_bv_efold
    clubb_params_single_col(:,iwpxp_Ri_exp) = clubb_wpxp_Ri_exp
    clubb_params_single_col(:,iz_displace) = clubb_z_displace

    ! Override clubb default
    if ( trim(subcol_scheme) == 'SILHS' ) then
      clubb_config_flags%saturation_formula = saturation_flatau
    else
      clubb_config_flags%saturation_formula = saturation_gfdl     ! Goff & Gratch (1946) approximation for SVP
    end if

    !  Set up CLUBB core.  Note that some of these inputs are overwritten
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.
    !! Initialize err_info with default values since info is not available here
    call init_default_err_info_api(1, err_info)
!$OMP PARALLEL
    call check_clubb_settings_api( 1, clubb_params_single_col,  & ! Intent(in)
                                   l_implemented,               & ! Intent(in)
                                   l_input_fields,              & ! Intent(in)
                                   clubb_config_flags,          & ! intent(in)
                                   err_info )                     ! Intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       call endrun('clubb_ini_cam: FATAL ERROR CALLING CHECK_CLUBB_SETTINGS_API')
    end if
!$OMP END PARALLEL

    ! Cleanup err_info since it is not needed anymore
    call cleanup_err_info_api(err_info)

    ! Print the list of CLUBB parameters
    if ( masterproc ) then
       do j = 1, nparams, 1
          write(iulog,*) params_list(j), " = ", clubb_params_single_col(1,j)
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
    call addfld ('RHO_CLUBB',        (/ 'lev' /),  'A', 'kg/m3',    'Air Density', sampled_on_subcycle=.true.)
    call addfld ('UP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Zonal Velocity Variance', sampled_on_subcycle=.true.)
    call addfld ('VP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Meridional Velocity Variance', sampled_on_subcycle=.true.)
    call addfld ('WP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',    'Vertical Velocity Variance', sampled_on_subcycle=.true.)
    call addfld ('WP2_ZT_CLUBB',     (/ 'lev' /),  'A', 'm2/s2',    'Vert Vel Variance on zt grid', sampled_on_subcycle=.true.)
    call addfld ('UPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',    'Zonal Momentum Flux', sampled_on_subcycle=.true.)
    call addfld ('VPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',    'Meridional Momentum Flux', sampled_on_subcycle=.true.)
    call addfld ('WP3_CLUBB',        (/ 'lev' /),  'A', 'm3/s3',    'Third Moment Vertical Velocity', sampled_on_subcycle=.true.)
    call addfld ('WPTHLP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',     'Heat Flux', sampled_on_subcycle=.true.)
    call addfld ('WPRTP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',     'Moisture Flux', sampled_on_subcycle=.true.)
    call addfld ('RTP2_CLUBB',       (/ 'ilev' /), 'A', 'kg^2/kg^2', 'Moisture Variance', sampled_on_subcycle=.true.)
    call addfld ('RTP2_ZT_CLUBB',    (/ 'lev' /),  'A', 'kg^2/kg^2','Moisture Variance on zt grid', sampled_on_subcycle=.true.)
    call addfld ('pdfp_rtp2_output_CLUBB',  (/ 'lev' /),  'A', 'kg^2/kg^2','PDF Rtot Variance', sampled_on_subcycle=.true.)
    call addfld ('THLP2_CLUBB',      (/ 'ilev' /), 'A', 'K^2',      'Temperature Variance', sampled_on_subcycle=.true.)
    call addfld ('THLP2_ZT_CLUBB',   (/ 'lev' /),  'A', 'K^2',      'Temperature Variance on zt grid', sampled_on_subcycle=.true.)
    call addfld ('RTPTHLP_CLUBB',    (/ 'ilev' /), 'A', 'K kg/kg',   'Temp. Moist. Covariance', sampled_on_subcycle=.true.)
    call addfld ('RCM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water Mixing Ratio', sampled_on_subcycle=.true.)
    call addfld ('RTM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Total Water Mixing Ratio', sampled_on_subcycle=.true.)
    call addfld ('THLM_CLUBB',       (/ 'lev' /),  'A', 'K',         'Liquid Water Potential Temperature', sampled_on_subcycle=.true.)
    call addfld ('WPRCP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',     'Liquid Water Flux', sampled_on_subcycle=.true.)
    call addfld ('CLOUDFRAC_CLUBB',  (/ 'lev' /),  'A', 'fraction', 'Cloud Fraction', sampled_on_subcycle=.true.)
    call addfld ('RCMINLAYER_CLUBB', (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water in Layer', sampled_on_subcycle=.true.)
    call addfld ('CLOUDCOVER_CLUBB', (/ 'lev' /),  'A', 'fraction', 'Cloud Cover', sampled_on_subcycle=.true.)
    call addfld ('WPTHVP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',     'Buoyancy Flux', sampled_on_subcycle=.true.)
    call addfld ('RVMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Water vapor tendency', sampled_on_subcycle=.true.)
    call addfld ('STEND_CLUBB',      (/ 'lev' /),  'A', 'J/(kg s)', 'Static energy tendency', sampled_on_subcycle=.true.)
    call addfld ('RCMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Liquid Water Tendency', sampled_on_subcycle=.true.)
    call addfld ('RIMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Ice Tendency', sampled_on_subcycle=.true.)
    call addfld ('UTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',   'U-wind Tendency', sampled_on_subcycle=.true.)
    call addfld ('VTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',   'V-wind Tendency', sampled_on_subcycle=.true.)
    call addfld ('ZT_CLUBB',         (/ 'lev' /),  'A', 'm',        'Thermodynamic Heights', sampled_on_subcycle=.true.)
    call addfld ('ZM_CLUBB',         (/ 'ilev' /), 'A', 'm',        'Momentum Heights', sampled_on_subcycle=.true.)
    call addfld ('UM_CLUBB',         (/ 'lev' /),  'A', 'm/s',      'Zonal Wind', sampled_on_subcycle=.true.)
    call addfld ('VM_CLUBB',         (/ 'lev' /),  'A', 'm/s',      'Meridional Wind', sampled_on_subcycle=.true.)
    call addfld ('WM_ZT_CLUBB',      (/ 'lev' /),  'A', 'm/s',      'Vertical Velocity', sampled_on_subcycle=.true.)
    call addfld ('PBLH',             horiz_only,   'A', 'm',        'PBL height', sampled_on_subcycle=.true.)
    call addfld ('CLDST',            (/ 'lev' /),  'A', 'fraction', 'Stratus cloud fraction', sampled_on_subcycle=.true.)
    call addfld ('ZMDLF',            (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained liquid water from ZM convection', sampled_on_subcycle=.true.)
    call addfld ('TTENDICE',         (/ 'lev' /),  'A', 'K/s',      'T tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('QVTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'Q tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('QITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'CLDICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('NITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'NUMICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)


    call addfld ('QCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'CLDICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('NCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'NUMICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('FQTENDICE',        (/ 'lev' /),  'A', 'fraction', 'Frequency of Ice Saturation Adjustment', sampled_on_subcycle=.true.)

    call addfld ('DPDLFLIQ',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained liquid water from deep convection', sampled_on_subcycle=.true.)
    call addfld ('DPDLFICE',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice from deep convection', sampled_on_subcycle=.true.)
    call addfld ('DPDLFT',           (/ 'lev' /),  'A', 'K/s',      'T-tendency due to deep convective detrainment', sampled_on_subcycle=.true.)
    call addfld ('RELVAR',           (/ 'lev' /),  'A', '-',        'Relative cloud water variance', sampled_on_subcycle=.true.)
    call addfld ('CLUBB_GRID_SIZE',  horiz_only,   'A', 'm',        'Horizontal grid box size seen by CLUBB', sampled_on_subcycle=.true.)


    call addfld ('ZMDLFI',           (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice water from ZM convection', sampled_on_subcycle=.true.)
    call addfld ('CONCLD',           (/ 'lev' /),  'A', 'fraction', 'Convective cloud cover', sampled_on_subcycle=.true.)
    call addfld ('CMELIQ',           (/ 'lev' /),  'A', 'kg/kg/s',  'Rate of cond-evap of liq within the cloud', sampled_on_subcycle=.true.)
    call addfld ('DETNLIQTND',       (/ 'lev' /),  'A', '1/kg/s',   'CLDNUM tendency in detrained water', sampled_on_subcycle=.true.)

    call addfld ('QSATFAC',          (/ 'lev' /),  'A', '-', 'Subgrid cloud water saturation scaling factor', sampled_on_subcycle=.true.)
    call addfld ('KVH_CLUBB',        (/ 'ilev' /), 'A', 'm2/s', 'CLUBB vertical diffusivity of heat/moisture on interface levels', sampled_on_subcycle=.true.)
    call addfld ('ELEAK_CLUBB',      horiz_only,   'A', 'W/m2', 'CLUBB energy leak', sampled_on_subcycle=.true.)
    call addfld ('TFIX_CLUBB',       horiz_only,   'A', 'K', 'Temperature increment to conserve energy', sampled_on_subcycle=.true.)

    ! ---------------------------------------------------------------------------- !
    ! Below are for detailed analysis of EDMF Scheme                               !
    ! ---------------------------------------------------------------------------- !
    if (do_clubb_mf) then
      call addfld ( 'edmf_DRY_A'    , (/ 'ilev' /), 'A', 'fraction', 'Dry updraft area fraction (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_A'  , (/ 'ilev' /), 'A', 'fraction', 'Moist updraft area fraction (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_DRY_W'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft vertical velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_W'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft vertical velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_DRY_QT'   , (/ 'ilev' /), 'A', 'kg/kg'   , 'Dry updraft total water mixing ratio (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_QT' , (/ 'ilev' /), 'A', 'kg/kg'   , 'Moist updraft total water mixing ratio (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_DRY_THL'  , (/ 'ilev' /), 'A', 'K'       , 'Dry updraft liquid-ice potential temperature (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_THL', (/ 'ilev' /), 'A', 'K'       , 'Moist updraft liquid-ice potential temperature (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_DRY_U'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft zonal velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_U'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft zonal velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_DRY_V'    , (/ 'ilev' /), 'A', 'm/s'     , 'Dry updraft meridional velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_V'  , (/ 'ilev' /), 'A', 'm/s'     , 'Moist updraft meridional velocity (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_MOIST_QC' , (/ 'ilev' /), 'A', 'kg/kg'   , 'Moist updraft condensate mixing ratio (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AE'     , (/ 'ilev' /), 'A', 'fraction', '1 minus sum of a_i*w_i (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AW'     , (/ 'ilev' /), 'A', 'm/s'     , 'Sum of a_i*w_i (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AWTHL'  , (/ 'ilev' /), 'A', 'K m/s'   , 'Sum of a_i*w_i*thl_i (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AWQT'   , (/ 'ilev' /), 'A', 'kgm/kgs' , 'Sum of a_i*w_i*q_ti (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AWU'    , (/ 'ilev' /), 'A', 'm2/s2'   , 'Sum of a_i*w_i*u_i (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_S_AWV'    , (/ 'ilev' /), 'A', 'm2/s2'   , 'Sum of a_i*w_i*v_i (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_thlflx'   , (/ 'ilev' /), 'A', 'W/m2'    , 'thl flux (EDMF)', sampled_on_subcycle=.true.)
      call addfld ( 'edmf_qtflx'    , (/ 'ilev' /), 'A', 'W/m2'    , 'qt flux (EDMF)', sampled_on_subcycle=.true.)
    end if

    if ( trim(subcol_scheme) /= 'SILHS' ) then
       ! hm_metadata is set up by calling init_pdf_hydromet_arrays_api in subcol_init_SILHS.
       ! So if we are not using silhs, we allocate the parts of hm_metadata that need allocating
       ! in order to making intel debug tests happy.
       allocate( hm_metadata%hydromet_list(1), stat=ierr)
       if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate hm_metadata%hydromet_list' )
       allocate( hm_metadata%l_mix_rat_hm(1), stat=ierr)
       if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate hm_metadata%l_mix_rat_hm' )
    end if

    !  Initialize statistics, below are dummy variables
    dum1 = 300._r8
    dum2 = 1200._r8
    dum3 = 300._r8

    if (stats_metadata%l_stats) then

      call stats_init_clubb( .true., dum1, dum2, &
                             nzm_clubb, nzt_clubb, nzm_clubb, dum3, &
                             stats_zt(:), stats_zm(:), stats_sfc(:), &
                             stats_rad_zt(:), stats_rad_zm(:))

       allocate(out_zt(pcols,pver,stats_zt(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate out_zt' )
       allocate(out_zm(pcols,pverp,stats_zm(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate out_zm' )
       allocate(out_sfc(pcols,1,stats_sfc(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate out_sfc' )

       if ( stats_metadata%l_output_rad_files ) then
          allocate(out_radzt(pcols,pver,stats_rad_zt(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate out_radzt' )
          allocate(out_radzm(pcols,pverp,stats_rad_zm(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) call endrun( 'clubb_ini_cam: Unable to allocate out_radzm' )
       end if

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
       call add_default('pdfp_rtp2_output_CLUBB',  1, ' ')
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
    endif

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
       call pbuf_set_field(pbuf2d, wp2rtp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wp2thlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, uprcp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, vprcp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rc_coef_zm_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, wp4_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wpup2_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wpvp2_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wp2up2_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wp2vp2_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, ice_supersat_idx, 0.0_r8)

       ! Initialize SILHS covariance contributions
       call pbuf_set_field(pbuf2d, rtp2_mc_zt_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, thlp2_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wprtp_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, wpthlp_mc_zt_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthlp_mc_zt_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, pdf_zm_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_mixt_frac_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d,  ttend_clubb_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  upwp_clubb_gw_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  vpwp_clubb_gw_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  thlp2_clubb_gw_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  wpthlp_clubb_gw_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d,  ttend_clubb_mc_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  upwp_clubb_gw_mc_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  vpwp_clubb_gw_mc_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  thlp2_clubb_gw_mc_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d,  wpthlp_clubb_gw_mc_idx, 0.0_r8)


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
    use perf_mod,       only: t_startf, t_stopf

#ifdef CLUBB_SGS
    use holtslag_boville_diff, only: hb_pbl_dependent_coefficients_run
    use spmd_utils, only: iam
    use clubb_api_module, only: &
      nparams, &
      calc_derrived_params_api, &
      check_parameters_api, &
      time_precision, &
      advance_clubb_core_api, &
      zt2zm_api, zm2zt_api, &
      setup_grid_heights_api, &
      em_min, &
      w_tol_sqd, &
      rt_tol, &
      thl_tol, &
      stats_begin_timestep_api, &
      calculate_thlp2_rad_api, update_xp2_mc_api, &
      sat_mixrat_liq_api, &
      fstderr, &
      ipdf_post_advance_fields, &
      copy_single_pdf_params_to_multi, &
      copy_multi_pdf_params_to_single, &
      pdf_parameter, &
      init_pdf_params_api, &
      init_pdf_implicit_coefs_terms_api, &
      setup_grid_api, &
      iiPDF_new, &
      iiPDF_new_hybrid

    ! Import setup for CLUBB error messaging
    use clubb_api_module, only: &
      clubb_fatal_error,    & ! Error code value to indicate a fatal error
      err_info_type,        &
      init_err_info_api,    &
      cleanup_err_info_api

    use cldfrc2m,                  only: aist_vector, rhmini_const, rhmaxi_const, rhminis_const, rhmaxis_const
    use cam_history,               only: outfld

    use macrop_driver,             only: liquid_macro_tend
    use clubb_mf,                  only: integrate_mf

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
    !                Pointers for pbuf                     !
    ! ---------------------------------------------------- !

    real(r8), pointer, dimension(:,:) :: wp2_pbuf      ! vertical velocity variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wp3_pbuf      ! third moment of vertical velocity		[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: wpthlp_pbuf   ! turbulent flux of thetal			[m/s K]
    real(r8), pointer, dimension(:,:) :: wprtp_pbuf    ! turbulent flux of moisture			[m/s kg/kg]
    real(r8), pointer, dimension(:,:) :: rtpthlp_pbuf  ! covariance of thetal and qt			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: rtp2_pbuf     ! moisture variance				[kg^2/kg^2]
    real(r8), pointer, dimension(:,:) :: thlp2_pbuf    ! temperature variance				[K^2]
    real(r8), pointer, dimension(:,:) :: rtp3_pbuf     ! moisture 3rd order				[kg^3/kg^3]
    real(r8), pointer, dimension(:,:) :: thlp3_pbuf    ! temperature 3rd order			[K^3]
    real(r8), pointer, dimension(:,:) :: up2_pbuf      ! east-west wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vp2_pbuf      ! north-south wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: up3_pbuf      ! east-west wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: vp3_pbuf      ! north-south wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: upwp_pbuf     ! east-west momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vpwp_pbuf     ! north-south momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wpthvp_pbuf   ! w'th_v' (momentum levels)			[m/s K]
    real(r8), pointer, dimension(:,:) :: wp2thvp_pbuf  ! w'^2 th_v' (thermodynamic levels)		[m^2/s^2 K]
    real(r8), pointer, dimension(:,:) :: rtpthvp_pbuf  ! r_t'th_v' (momentum levels)			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: thlpthvp_pbuf ! th_l'th_v' (momentum levels)			[K^2]
    real(r8), pointer, dimension(:,:) :: cloud_frac_pbuf ! Cloud fraction (thermodynamic levels)	[K^2]
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_1_pbuf        !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_2_pbuf        !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_1_pbuf !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_2_pbuf !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_mixt_frac_pbuf  !work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: wp2rtp_pbuf    ! w'^2 rt' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2thlp_pbuf   ! w'^2 thl' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: uprcp_pbuf     ! < u' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: vprcp_pbuf     ! < v' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: rc_coef_zm_pbuf   ! Coef. of X'r_c' in Eq. (34) (t-levs.)
    real(r8), pointer, dimension(:,:) :: wp4_pbuf       ! w'^4 (momentum levels
    real(r8), pointer, dimension(:,:) :: wpup2_pbuf     ! w'u'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wpvp2_pbuf     ! w'v'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2up2_pbuf    ! w'^2 u'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: wp2vp2_pbuf    ! w'^2 v'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: thlm_pbuf     ! mean temperature				[K]
    real(r8), pointer, dimension(:,:) :: rtm_pbuf      ! mean moisture mixing ratio			[kg/kg]
    real(r8), pointer, dimension(:,:) :: rcm_pbuf      ! CLUBB cloud water mixing ratio               [kg/kg]
    real(r8), pointer, dimension(:,:) :: um_pbuf       ! mean east-west wind				[m/s]
    real(r8), pointer, dimension(:,:) :: vm_pbuf       ! mean north-south wind			[m/s]
    real(r8), pointer, dimension(:,:) :: cld_pbuf      ! cloud fraction 				[fraction]
    real(r8), pointer, dimension(:,:) :: concld_pbuf   ! convective cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: ast_pbuf      ! stratiform cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: alst_pbuf     ! liquid stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: aist_pbuf     ! ice stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: qlst_pbuf     ! Physical in-stratus LWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: qist_pbuf     ! Physical in-stratus IWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: deepcu_pbuf   ! deep convection cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: shalcu_pbuf   ! shallow convection cloud fraction 		[fraction]
    real(r8), pointer, dimension(:,:) :: khzm_pbuf     ! CLUBB's eddy diffusivity of heat/moisture on momentum (i.e. interface) levels          [m^2/s]
    real(r8), pointer, dimension(:)   :: pblh_pbuf     ! planetary boundary layer height                [m]
    real(r8), pointer, dimension(:,:) :: tke_pbuf      ! turbulent kinetic energy                     [m^2/s^2]
    real(r8), pointer, dimension(:,:) :: dp_icwmr_pbuf ! deep convection in cloud mixing ratio        [kg/kg]
    real(r8), pointer, dimension(:,:) :: ice_supersat_frac_pbuf ! Cloud fraction of ice clouds (pver)[fraction]
    real(r8), pointer, dimension(:,:) :: relvar_pbuf   ! relative cloud water variance                [-]
    real(r8), pointer, dimension(:,:) :: accre_enhan_pbuf ! accretion enhancement factor              [-]
    real(r8), pointer, dimension(:,:) :: naai_pbuf
    real(r8), pointer, dimension(:,:) :: cmeliq_pbuf
    real(r8), pointer, dimension(:,:) :: cmfmc_sh_pbuf ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]

    real(r8), pointer, dimension(:,:) :: qsatfac_pbuf
    real(r8), pointer, dimension(:,:) :: npccn_pbuf
    real(r8), pointer, dimension(:,:) :: prer_evap_pbuf
    real(r8), pointer, dimension(:,:) :: qrl_pbuf

    ! SILHS covariance contributions
    real(r8), pointer, dimension(:,:) :: rtp2_mc_zt_pbuf
    real(r8), pointer, dimension(:,:) :: thlp2_mc_zt_pbuf
    real(r8), pointer, dimension(:,:) :: wprtp_mc_zt_pbuf
    real(r8), pointer, dimension(:,:) :: wpthlp_mc_zt_pbuf
    real(r8), pointer, dimension(:,:) :: rtpthlp_mc_zt_pbuf

    ! Connections to Gravity Wave parameterization
    real(r8), pointer, dimension(:,:) :: ttend_clubb_pbuf
    real(r8), pointer, dimension(:,:) :: upwp_clubb_gw_pbuf
    real(r8), pointer, dimension(:,:) :: vpwp_clubb_gw_pbuf
    real(r8), pointer, dimension(:,:) :: thlp2_clubb_gw_pbuf
    real(r8), pointer, dimension(:,:) :: wpthlp_clubb_gw_pbuf

    real(r8), pointer, dimension(:,:) :: ttend_clubb_mc_pbuf
    real(r8), pointer, dimension(:,:) :: upwp_clubb_gw_mc_pbuf
    real(r8), pointer, dimension(:,:) :: vpwp_clubb_gw_mc_pbuf
    real(r8), pointer, dimension(:,:) :: thlp2_clubb_gw_mc_pbuf
    real(r8), pointer, dimension(:,:) :: wpthlp_clubb_gw_mc_pbuf

    ! ---------------------------------------------------- !
    !                   Local Variables                    !
    ! ---------------------------------------------------- !

    integer :: i !Must be delcared outside "CLUBB_SGS" ifdef for det_s and det_ice zero-ing loops

#ifdef CLUBB_SGS

    real(r8), parameter :: &
      rad2deg=180.0_r8/pi

    character(len=*), parameter :: subr='clubb_tend_cam'

    type(physics_state) :: state1                ! Local copy of state variable
    type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header

    type(grid) :: &
      gr, gr_a          ! CLUBB grid data structure

    type(nu_vertical_res_dep) :: &
      nu_vert_res_dep   ! Vertical resolution dependent nu values

    real(r8), dimension(state%ncol,nparams) :: &
      clubb_params    ! Adjustable CLUBB parameters (C1, C2 ...)

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLES SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol) :: &
      deltaz, &
      fcor, &                             ! Coriolis forcing 			      	[s^-1]
      sfc_elevation, &    		            ! Elevation of ground			      	[m AMSL][m]
      wpthlp_sfc, &                       ! w' theta_l' at surface                      [(m K)/s]
      wprtp_sfc, &                        ! w' r_t' at surface                          [(kg m)/( kg s)]
      upwp_sfc, &                         ! u'w' at surface                             [m^2/s^2]
      vpwp_sfc, &                         ! v'w' at surface                             [m^2/s^2]
      p_sfc, &                            ! pressure at surface                         [Pa]
      upwp_sfc_pert, &                    ! perturbed u'w' at surface                   [m^2/s^2]
      vpwp_sfc_pert, &                    ! perturbed v'w' at surface                   [m^2/s^2]
      grid_dx, grid_dy                    ! CAM grid [m]

    real(r8), dimension(state%ncol,sclr_dim) :: &
      wpsclrp_sfc            ! Scalar flux at surface                        [{units vary} m/s]

    real(r8), dimension(state%ncol,edsclr_dim) :: &
      wpedsclrp_sfc        ! Eddy-scalar flux at surface                   [{units vary} m/s]


    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzt_clubb) :: &
      thlm_forcing,             & ! theta_l forcing (thermodynamic levels)      [K/s]
      rtm_forcing,              & ! r_t forcing (thermodynamic levels)          [(kg/kg)/s]
      um_forcing,               & ! u wind forcing (thermodynamic levels)     	[m/s/s]
      vm_forcing,               & ! v wind forcing (thermodynamic levels)     	[m/s/s]
      wm_zt,                    & ! w mean wind component on thermo. levels   	[m/s]
      rtm_ref,                  & ! Initial profile of rtm                      [kg/kg]
      thlm_ref,                 & ! Initial profile of thlm                     [K]
      um_ref,                   & ! Initial profile of um                       [m/s]
      vm_ref,                   & ! Initial profile of vm                       [m/s]
      ug,                       & ! U geostrophic wind                          [m/s]
      vg,                       & ! V geostrophic wind                          [m/s]
      p_in_Pa,                  & ! Air pressure (thermodynamic levels)       	[Pa]
      rho_zt,                   & ! Air density on thermo levels                [kg/m^3]
      exner,                    & ! Exner function (thermodynamic levels)       [-]
      rho_ds_zt,                & ! Dry, static density on thermodynamic levels 	[kg/m^3]
      invrs_rho_ds_zt,          & ! Inv. dry, static density on thermo. levels  	[m^3/kg]
      thv_ds_zt,                & ! Dry, base-state theta_v on thermo. levels   	[K]
      rfrzm,                    &
      rvm_in,                   & ! water vapor mixing ratio                      [kg/kg]
      rtp2_zt,                  & ! CLUBB R-tot variance on thermo levs
      thl2_zt,                  & ! CLUBB Theta-l variance on thermo levs         [K^2]
      wp2_zt,                   & ! CLUBB W variance on theromo levs              [m^2/s^2]
      cloud_frac_inout,         & ! CLUBB output of cloud fraction		[fraction]
      um_pert_inout,            & ! Perturbed U wind                          [m/s]
      vm_pert_inout,            & ! Perturbed V wind                          [m/s]
      khzt_out,                 & ! eddy diffusivity on thermo grids              [m^2/s]
      w_up_in_cloud_out,        &
      w_down_in_cloud_out,      &
      cloudy_updraft_frac_out,  &
      cloudy_downdraft_frac_out,&
      rcm_in_layer,         & ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
      cloud_cover_out,          & ! CLUBB output of in-cloud cloud fraction	[fraction]
      pre_in,                   & ! input for precip evaporation
      qrl_clubb,                &
      qclvar_out,               & ! cloud water variance                          [kg^2/kg^2]
      zt_g,                     & ! Thermodynamic grid of CLUBB		      	[m]
      Lscale,                   &

      ! MF local thermodynamic vars
                  invrs_dzt,     & ! thermodynamic grid
                  invrs_exner_zt,& ! thermodynamic grid
      kappa_zt,   qc_zt            ! thermodynamic grid

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLES SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzm_clubb) :: &
      thlp2_rad,                &
      wprtp_forcing,            &
      wpthlp_forcing,           &
      rtp2_forcing,             &
      thlp2_forcing,            &
      rtpthlp_forcing,          &
      wm_zm,                    & ! w mean wind component on momentum levels  	[m/s]
      rho_zm,                   & ! Air density on momentum levels              [kg/m^3]
      rho_ds_zm,                & ! Dry, static density on momentum levels      	[kg/m^3]
      invrs_rho_ds_zm,          & ! Inv. dry, static density on momentum levels 	[m^3/kg]
      thv_ds_zm,                & ! Dry, base-state theta_v on momentum levels  	[K]
      upwp_pert_inout,          & ! Perturbed u'w'                            [m^2/s^2]
      vpwp_pert_inout,          & ! Perturbed v'w'                            [m^2/s^2]
      khzm_out,                 & ! Eddy diffusivity of heat/moisture on momentum (i.e. interface) levels  [m^2/s]
      thlprcp_out,              &
      wprcp_out,                & ! CLUBB output of flux of liquid water		[kg/kg m/s]
      invrs_tau_zm_out,         & ! CLUBB output of 1 divided by time-scale	[1/s]
      rtp2_mc_out,              & ! total water tendency from rain evap
      thlp2_mc_out,             & ! thetal tendency from rain evap
      wprtp_mc_out,             &
      wpthlp_mc_out,            &
      rtpthlp_mc_out,           &
      uprcp_inout,              & ! < u' r_c' > (momentum levels)
      vprcp_inout,              & ! < v' r_c' > (momentum levels)
      zi_g,                     & ! Momentum grid of CLUBB		      	[m]

      ! MF Plume
      mf_dry_a,   mf_moist_a,    &
      mf_dry_w,   mf_moist_w,    &
      mf_dry_qt,  mf_moist_qt,   &
      mf_dry_thl, mf_moist_thl,  &
      mf_dry_u,   mf_moist_u,    &
      mf_dry_v,   mf_moist_v,    &
                  mf_moist_qc,   &
      s_ae,       s_aw,          &
      s_awthl,    s_awqt,        &
      s_awql,     s_awqi,        &
      s_awu,      s_awv,         &
      mf_thlflx,  mf_qtflx,      &

      ! MF local momentum vars
      rtm_zm_in,  thlm_zm_in,    & ! momentum grid
      kappa_zm,   p_in_Pa_zm,    & ! momentum grid
                  invrs_exner_zm   ! momentum grid
      
    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzt_clubb,sclr_dim) :: &
      sclrm_forcing,  & ! Passive scalar forcing              [{units vary}/s]
      sclrm,          & ! Passive scalar mean (thermo. levels)          [units vary]
      sclrp3            ! sclr'^3 (thermo. levels)                      [{units vary}^3]
      
    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzm_clubb,sclr_dim) :: &
      sclrp2,         & ! sclr'^2 (momentum levels)                     [{units vary}^2]
      sclrprtp,       & ! sclr'rt' (momentum levels)          [{units vary} (kg/kg)]
      sclrpthlp,      & ! sclr'thlp' (momentum levels)        [{units vary} (K)]
      wpsclrp,        & ! w'sclr' (momentum levels)                     [{units vary} m/s]
      sclrpthvp_inout   ! sclr'th_v' (momentum levels)                  [{units vary} (K)]

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzt_clubb,edsclr_dim) :: &
      edsclrm_forcing,  & ! Eddy passive scalar forcing         [{units vary}/s]
      edsclr_in           ! Scalars to be diffused through CLUBB 		[units vary]

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzt_clubb,hydromet_dim) :: &
      wp2hmp,       &
      rtphmp_zt,    &
      thlphmp_zt

    ! Local CLUBB variables dimensioned as NCOL (only useful columns) to be sent into the clubb run api
    ! NOTE: THESE VARIABLS SHOULD NOT BE USED IN PBUF OR OUTFLD (HISTORY) SUBROUTINES
    real(r8), dimension(state%ncol,nzm_clubb,hydromet_dim) :: &
      wphydrometp


    ! Variables used for output (zm)
    real(r8), dimension(pcols,pverp) :: &
      zi_output, & ! output for momentum CLUBB grid                    [m]
      wpthlp_output, &  ! Heat flux output variable                     [W/m2]
      rtpthlp_output, & ! rtpthlp ouptut                                [K kg/kg]
      wprtp_output, &       ! Total water flux output variable              [W/m2]
      wp2_output, &
      up2_output, &
      vp2_output, &
      upwp_output, &
      vpwp_output, &
      rtp2_output, &
      wprcp_clubb_output, &
      wpthvp_clubb_output, &
      thlp2_output,  &
      dlf_liq_out, & ! Detrained liquid water from ZM                [kg/kg/s]
      dlf_ice_out, & ! Detrained ice water from ZM                   [kg/kg/s]

      ! MF outputs to outfld
      ! NOTE: Arrays of size PCOLS (all possible columns) can be used to access State, PBuf and History Subroutines
      mf_dry_a_output,   mf_moist_a_output,   &
      mf_dry_w_output,   mf_moist_w_output,   &
      mf_dry_qt_output,  mf_moist_qt_output,  &
      mf_dry_thl_output, mf_moist_thl_output, &
      mf_dry_u_output,   mf_moist_u_output,   &
      mf_dry_v_output,   mf_moist_v_output,   &
                          mf_moist_qc_output,  &
      s_ae_output,       s_aw_output,         &
      s_awthl_output,    s_awqt_output,       &
      s_awql_output,     s_awqi_output,       &
      s_awu_output,      s_awv_output,        &
      mf_thlflx_output,  mf_qtflx_output

    ! Variables used for output (zt)
    real(r8), dimension(pcols,pver) :: &
      rvmtend_clubb_output, &
      rcmtend_clubb_output, &
      rimtend_clubb_output, &
      stend_clubb_output, &
      utend_clubb_output, &
      vtend_clubb_output, &
      dpdlfliq_output, &
      dpdlfice_output, &
      dpdlft_output, &
      detnliquid_output, &
      zt_output, &  ! output for the thermo CLUBB grid                   [m]
      rtp2_zt_output, &         ! CLUBB R-tot variance on thermo levs           [kg^2/kg^2]
      wp3_output, &     ! wp3 output                                    [m^3/s^3]
      thl2_zt_output, &     ! CLUBB Theta-l variance on thermo levs
      wp2_zt_output, & 
      rcm_in_layer_output, &  	! CLUBB in-cloud liquid water mixing ratio	[kg/kg]
      pdfp_rtp2_output, &      ! Calculated R-tot variance from pdf_params     [kg^2/kg^2]
      wm_zt_output, &        ! CLUBB mean W on thermo levs output            [m/s]
      rcm_output, &
      rtm_output, &
      thlm_output, &
      um_output, &
      vm_output

    real(r8), dimension(pcols) :: &
      rhmini, &
      rhmaxi, &
      rtm_integral_vtend, &
      rtm_integral_ltend, &
      se_dis, &
      eleak, &
      ustar2, &                    ! Surface stress for PBL height                 [m2/s2]
      obklen, &                    ! Obukov length                                 [m]
      kbfs, &                      ! Kinematic Surface heat flux                   [K m/s]
      kinheat, &                   ! Kinematic Surface heat flux                   [K m/s]
      rrho, &                      ! Inverse of air density                        [1/kg/m^3]
      kinwat, &                       ! Kinematic water vapor flux                    [m/s]
      dummy2, &                    ! dummy variable                                [units vary]
      dummy3                    ! dummy variable                                [units vary]

    real(r8), dimension(pcols,pver) :: &
      temp2d, &  ! temporary array for holding scaled outputs
      qitend, &
      initend, &  ! Needed for ice supersaturation adjustment calculation
      stend, &
      qvtend, &
      qctend, &
      inctend, &
      qclvar, & ! cloud water variance                          [kg^2/kg^2]
      dz_g, & ! thickness of layer                            [m]
      clubb_s, &
      inv_exner_clubb, &  ! Inverse exner function consistent with CLUBB  [-]
      thv, &     ! virtual potential temperature			[K]
      th     ! potential temperature                         [K]

    real(r8), dimension(pcols,pverp) :: &
      rho     		! Midpoint density in CAM      			[kg/m^3]

    real(r8) :: &
      dlf2, & ! Detraining cld H20 from shallow convection    [kg/kg/day]
      dum1, &   ! dummy variable                                [units vary]
      invrs_hdtime, &
      lmin, &
      mixt_frac_max_mag, &
      dtime, &				        ! CLUBB time step                               [s]
      ubar, &				          ! surface wind                                [m/s]
      ustar, &				          ! surface stress				[m/s]
      bflx22, &                          ! Variable for buoyancy flux for pbl            [K m/s]
      zo, &                               ! roughness height                              [m]
      relvarmax, &
      frac_limit, ic_limit, &
      mean_rt, &                          ! Calculated R-tot mean from pdf_params (temp)  [kg/kg]
      latsub, &
      apply_const, &
      dl_rad, di_rad, dt_low, &
      ! Variables below are needed to compute energy integrals for conservation
      te_a, se_a, ke_a, wv_a, wl_a, &
      te_b, se_b, ke_b, wv_b, wl_b

    intrinsic :: max

    logical, dimension(pcnst) :: &
      lq2, &
      lqice

    character(len=200) :: temp1, sub             ! Strings needed for CLUBB output
    character(len=512) :: errmsg

    integer, dimension(pcols) :: &
      clubbtop, &
      troplev

    integer :: &
      errflg, &
      j, k, t, ixind, nadv, n,      & ! Loop variables
      k_cam, k_clubb, sclr, edsclr, & ! Loop variables
      ixcldice, ixcldliq, ixnumliq, ixnumice, ixq, &
      itim_old, &
      ncol, lchnk, &  ! # of columns, and chunk identifier
      icnt, &
      stats_nsamp, stats_nout           ! Stats sampling and output intervals for CLUBB [timestep]

#endif

  call t_startf('clubb_tend_cam')

  do i = 1, pcols
    det_s(i)   = 0.0_r8
    det_ice(i) = 0.0_r8
  end do

#ifdef CLUBB_SGS

#ifdef _OPENACC
    ! These options have not been GPUized
    if ( do_clubb_mf ) call endrun(subr//': do_clubb_mf=.true. not available when compiling with OpenACC')
    if ( do_rainturb ) call endrun(subr//': do_rainturb=.true. not available when compiling with OpenACC')
    if ( do_cldcool )  call endrun(subr//': do_cldcool=.true. not available when compiling with OpenACC')
    if ( clubb_do_icesuper )  call endrun(subr//': clubb_do_icesuper=.true. not available when compiling with OpenACC')
    if ( single_column .and. .not. scm_cambfb_mode )  then
      call endrun(subr//': (single_column && !scm_cambfb_mode)=.true. not available when compiling with OpenACC')
    end if
#endif

    !-----------------------------------------------------------------------------------!
    !                           MAIN COMPUTATION BEGINS HERE                            !
    !-----------------------------------------------------------------------------------!
    print *, "top_lev = ", top_lev 
    print *, "trop_cloud_top_press = ", trop_cloud_top_press 

    call t_startf('clubb_tend_cam:NAR')

    !  Get indicees for cloud and ice mass and cloud and ice number
    call cnst_get_ind('Q',ixq)
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)
    call cnst_get_ind('NUMLIQ',ixnumliq)
    call cnst_get_ind('NUMICE',ixnumice)

    !  Determine time step of physics buffer
    itim_old = pbuf_old_tim_idx()

    !  Establish associations between pointers and physics buffer fields
    call pbuf_get_field(pbuf, wp2_idx,     wp2_pbuf )
    call pbuf_get_field(pbuf, wp3_idx,     wp3_pbuf )
    call pbuf_get_field(pbuf, wpthlp_idx,  wpthlp_pbuf )
    call pbuf_get_field(pbuf, wprtp_idx,   wprtp_pbuf )
    call pbuf_get_field(pbuf, rtpthlp_idx, rtpthlp_pbuf )
    call pbuf_get_field(pbuf, rtp2_idx,    rtp2_pbuf )
    call pbuf_get_field(pbuf, thlp2_idx,   thlp2_pbuf )
    call pbuf_get_field(pbuf, up2_idx,     up2_pbuf )
    call pbuf_get_field(pbuf, vp2_idx,     vp2_pbuf )

    call pbuf_get_field(pbuf, rtp3_idx,    rtp3_pbuf )
    call pbuf_get_field(pbuf, thlp3_idx,   thlp3_pbuf )
    call pbuf_get_field(pbuf, up3_idx,     up3_pbuf )
    call pbuf_get_field(pbuf, vp3_idx,     vp3_pbuf )

    call pbuf_get_field(pbuf, upwp_idx,    upwp_pbuf )
    call pbuf_get_field(pbuf, vpwp_idx,    vpwp_pbuf )
    call pbuf_get_field(pbuf, wpthvp_idx,  wpthvp_pbuf)
    call pbuf_get_field(pbuf, wp2thvp_idx, wp2thvp_pbuf)
    call pbuf_get_field(pbuf, rtpthvp_idx, rtpthvp_pbuf)
    call pbuf_get_field(pbuf, thlpthvp_idx,thlpthvp_pbuf)
    call pbuf_get_field(pbuf, rcm_idx,     rcm_pbuf)
    call pbuf_get_field(pbuf, cloud_frac_idx, cloud_frac_pbuf)

    call pbuf_get_field(pbuf, pdf_zm_w_1_idx,         pdf_zm_w_1_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_w_2_idx,         pdf_zm_w_2_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_1_idx,  pdf_zm_varnce_w_1_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_2_idx,  pdf_zm_varnce_w_2_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_mixt_frac_idx,   pdf_zm_mixt_frac_pbuf )

    call pbuf_get_field(pbuf, wp2rtp_idx, wp2rtp_pbuf)
    call pbuf_get_field(pbuf, wp2thlp_idx, wp2thlp_pbuf)
    call pbuf_get_field(pbuf, uprcp_idx, uprcp_pbuf)
    call pbuf_get_field(pbuf, vprcp_idx, vprcp_pbuf)
    call pbuf_get_field(pbuf, rc_coef_zm_idx, rc_coef_zm_pbuf)
    call pbuf_get_field(pbuf, wp4_idx, wp4_pbuf)
    call pbuf_get_field(pbuf, wpup2_idx, wpup2_pbuf)
    call pbuf_get_field(pbuf, wpvp2_idx, wpvp2_pbuf)
    call pbuf_get_field(pbuf, wp2up2_idx, wp2up2_pbuf)
    call pbuf_get_field(pbuf, wp2vp2_idx, wp2vp2_pbuf)
    call pbuf_get_field(pbuf, thlm_idx,    thlm_pbuf )
    call pbuf_get_field(pbuf, rtm_idx,     rtm_pbuf )
    call pbuf_get_field(pbuf, um_idx,      um_pbuf )
    call pbuf_get_field(pbuf, vm_idx,      vm_pbuf )

    call pbuf_get_field(pbuf, tke_idx,     tke_pbuf)
    call pbuf_get_field(pbuf, qrl_idx,     qrl_pbuf)

    call pbuf_get_field(pbuf, cld_idx,     cld_pbuf,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, concld_idx,  concld_pbuf,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, ast_idx,     ast_pbuf,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, alst_idx,    alst_pbuf,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, aist_idx,    aist_pbuf,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, qlst_idx,    qlst_pbuf,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, qist_idx,    qist_pbuf,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

    call pbuf_get_field(pbuf, qsatfac_idx, qsatfac_pbuf)

    call pbuf_get_field(pbuf, prer_evap_idx, prer_evap_pbuf)
    call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan_pbuf)
    call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq_pbuf)
    call pbuf_get_field(pbuf, ice_supersat_idx, ice_supersat_frac_pbuf)
    call pbuf_get_field(pbuf, relvar_idx,  relvar_pbuf)
    call pbuf_get_field(pbuf, dp_frac_idx, deepcu_pbuf)
    call pbuf_get_field(pbuf, sh_frac_idx, shalcu_pbuf)
    call pbuf_get_field(pbuf, kvh_idx,     khzm_pbuf)
    call pbuf_get_field(pbuf, pblh_idx,    pblh_pbuf)
    call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr_pbuf)
    call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh_pbuf)

    ! SILHS covariance contributions
    call pbuf_get_field(pbuf, rtp2_mc_zt_idx,    rtp2_mc_zt_pbuf)
    call pbuf_get_field(pbuf, thlp2_mc_zt_idx,   thlp2_mc_zt_pbuf)
    call pbuf_get_field(pbuf, wprtp_mc_zt_idx,   wprtp_mc_zt_pbuf)
    call pbuf_get_field(pbuf, wpthlp_mc_zt_idx,  wpthlp_mc_zt_pbuf)
    call pbuf_get_field(pbuf, rtpthlp_mc_zt_idx, rtpthlp_mc_zt_pbuf)

    ! For Gravity Wave
    call pbuf_get_field(pbuf, ttend_clubb_idx,       ttend_clubb_pbuf )
    call pbuf_get_field(pbuf, thlp2_clubb_gw_idx,    thlp2_clubb_gw_pbuf )
    call pbuf_get_field(pbuf, upwp_clubb_gw_idx,     upwp_clubb_gw_pbuf )
    call pbuf_get_field(pbuf, vpwp_clubb_gw_idx,     vpwp_clubb_gw_pbuf )
    call pbuf_get_field(pbuf, wpthlp_clubb_gw_idx,   wpthlp_clubb_gw_pbuf )

    call pbuf_get_field(pbuf, ttend_clubb_mc_idx,     ttend_clubb_mc_pbuf )
    call pbuf_get_field(pbuf, thlp2_clubb_gw_mc_idx,  thlp2_clubb_gw_mc_pbuf )
    call pbuf_get_field(pbuf, upwp_clubb_gw_mc_idx,   upwp_clubb_gw_mc_pbuf )
    call pbuf_get_field(pbuf, vpwp_clubb_gw_mc_idx,   vpwp_clubb_gw_mc_pbuf )
    call pbuf_get_field(pbuf, wpthlp_clubb_gw_mc_idx, wpthlp_clubb_gw_mc_pbuf )

    if (clubb_do_icesuper) then
      call pbuf_get_field(pbuf, naai_idx, naai_pbuf)
    end if

    !  Initialize physics tendency arrays
    call physics_ptend_init(ptend_all, state%psetcols, 'clubb')

    ! Copy the state to state1 array to use in this routine
    call physics_state_copy(state, state1)

    ! Constituents are all treated as dry mmr by clubb.  Convert the water species to
    ! a dry basis.
    call set_wet_to_dry(state1, convert_cnst_type='wet')

    if (clubb_do_liqsupersat) then
      call pbuf_get_field(pbuf, npccn_idx, npccn_pbuf)
    endif

    ! Define the grid box size.  CLUBB needs this information to determine what
    !  the maximum length scale should be.  This depends on the column for
    !  variable mesh grids and lat-lon grids
    call grid_size(state1, grid_dx, grid_dy)

    ! Determine number of columns and which chunk computation is to be performed on
    ncol = state%ncol
    lchnk = state%lchnk

    ! Allocate pdf_params only if they aren't allocated already.
    if ( .not. allocated(pdf_params_chnk(lchnk)%mixt_frac) ) then
      call init_pdf_params_api( nzt_clubb, ncol, pdf_params_chnk(lchnk) )
    end if

    ! pdf_params_zm are only used if l_call_pdf_closure_twice=.true.
    if ( clubb_config_flags%l_call_pdf_closure_twice ) then
      if ( .not. allocated(pdf_params_zm_chnk(lchnk)%mixt_frac) ) then
        call init_pdf_params_api( nzm_clubb, ncol, pdf_params_zm_chnk(lchnk) )
      end if
    end if

    ! pdf_implicit_coefs_terms are only used if iiPDF_type = iiPDF_new or iiPDF_new_hybrid
    if ( clubb_config_flags%iiPDF_type == iiPDF_new         .or. &
         clubb_config_flags%iiPDF_type == iiPDF_new_hybrid       ) then

      if ( .not. allocated(pdf_implicit_coefs_terms_chnk(lchnk)%coef_wp4_implicit) ) then
        call init_pdf_implicit_coefs_terms_api( nzt_clubb, ncol, sclr_dim, &
                                                pdf_implicit_coefs_terms_chnk(lchnk) )
      end if

    end if

    ! Initialize err_info with parallelization and geographical info
    call init_err_info_api(ncol, lchnk, iam, state1%lat*rad2deg, state1%lon*rad2deg, err_info)

    !--------------------- Scalar Setting --------------------

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

    ! Precalculte the hdtime inverse
    invrs_hdtime = 1._r8 / hdtime


    !  Set stats output and increment equal to CLUBB and host dt
    stats_metadata%stats_tsamp = dtime
    stats_metadata%stats_tout  = hdtime

    stats_nsamp = nint(stats_metadata%stats_tsamp/dtime)
    stats_nout = nint(stats_metadata%stats_tout/dtime)


    if (clubb_do_adv) then
      apply_const = 1._r8  ! Initialize to one, only if CLUBB's moments are advected
    else
      apply_const = 0._r8  ! Never want this if CLUBB's moments are not advected
    endif

    ! Initialize the apply_const variable (note special logic is due to eulerian backstepping)
    if (clubb_do_adv .and. (is_first_step() .or. all(wpthlp_pbuf(1:ncol,:)  ==  0._r8))) then
      apply_const = 0._r8  ! On first time through do not remove constant
                           !  from moments since it has not been added yet
    endif

    !--------------------- Initializations --------------------

    !  Set the ztodt timestep in pbuf for SILHS, this is needed because hdtime is not input to silhs
    ztodt = 1.0_r8 * hdtime

    call t_stopf('clubb_tend_cam:NAR')
    call t_startf('clubb_tend_cam:acc_copyin')
    !$acc data copyin( sclr_idx, clubb_params_single_col, grid_dx, grid_dy, rairv, cpairv, qrl_pbuf, &
    !$acc              pdf_params_chnk(lchnk), pdf_params_zm_chnk(lchnk), &
    !$acc              state1, state1%q, state1%u, state1%v, state1%t, state1%pmid, state1%s, state1%pint, &
    !$acc              state1%zm, state1%zi, state1%pdeldry, state1%pdel, state1%omega, state1%phis, &
    !$acc              cam_in, cam_in%shf, cam_in%wsx, cam_in%wsy, cam_in%cflx, &
    !$acc              rrho, prer_evap_pbuf, rtp2_mc_zt_pbuf, thlp2_mc_zt_pbuf, wprtp_mc_zt_pbuf, wpthlp_mc_zt_pbuf, rtpthlp_mc_zt_pbuf, &
    !$acc              err_info, err_info%err_header ) &
    !$acc        copy( um_pbuf, vm_pbuf, upwp_pbuf, vpwp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, rtpthvp_pbuf, thlpthvp_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, &
    !$acc              wp2_pbuf, wp3_pbuf, rtp2_pbuf, thlp2_pbuf, rtp3_pbuf, thlp3_pbuf, thlm_pbuf, rtm_pbuf, wprtp_pbuf, wpthlp_pbuf, rtpthlp_pbuf, &
    !$acc              pdf_zm_w_1_pbuf, pdf_zm_w_2_pbuf, pdf_zm_varnce_w_1_pbuf, pdf_zm_varnce_w_2_pbuf, pdf_zm_mixt_frac_pbuf, &
    !$acc              cloud_frac_pbuf, wp2rtp_pbuf, wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, &
    !$acc              wp2up2_pbuf, wp2vp2_pbuf, ice_supersat_frac_pbuf, &
    !$acc              pdf_params_zm_chnk(lchnk)%w_1, pdf_params_zm_chnk(lchnk)%w_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%varnce_w_1, pdf_params_zm_chnk(lchnk)%varnce_w_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%mixt_frac ) &
    !$acc     copyout( temp2d, inv_exner_clubb, &
    !$acc              rcm_pbuf, khzm_pbuf, qclvar, thv, dz_g, &
    !$acc              clubbtop, se_dis, eleak, clubb_s ) &
    !$acc      create( upwp_sfc_pert, vpwp_sfc_pert, khzt_out, khzm_out, &
    !$acc              fcor, &
    !$acc              thlp3_in, rvm_in, cloud_frac_inout, &
    !$acc              uprcp_inout, vprcp_inout, &
    !$acc              pre_in, kappa_zt, qc_zt, invrs_exner_zt, kappa_zm, p_in_Pa_zm, &
    !$acc              invrs_exner_zm, cloud_cover_out, rcm_in_layer, wprcp_out, &
    !$acc              qclvar_out, w_up_in_cloud_out, cloudy_downdraft_frac_out, &
    !$acc              w_down_in_cloud_out, invrs_tau_zm_out, vm_pert_inout, upwp_pert_inout, vpwp_pert_inout, &
    !$acc              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
    !$acc              wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    !$acc              rtpthlp_forcing, wm_zm, wm_zt, rho_zm, rho_zt, rho_ds_zm, rho_ds_zt, &
    !$acc              invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, rfrzm, &
    !$acc              wpthlp_sfc, clubb_params, sfc_elevation, wprtp_sfc, upwp_sfc, vpwp_sfc, &
    !$acc              rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, p_in_Pa, exner, um_pert_inout, &
    !$acc              thlprcp_out, deltaz, zi_g, zt_g, qrl_clubb, p_sfc, &
    !$acc              err_info%err_code, &
    !$acc              pdf_params_chnk(lchnk)%w_1, pdf_params_chnk(lchnk)%w_2, &
    !$acc              pdf_params_chnk(lchnk)%varnce_w_1, pdf_params_chnk(lchnk)%varnce_w_2, &
    !$acc              pdf_params_chnk(lchnk)%rt_1, pdf_params_chnk(lchnk)%rt_2, &
    !$acc              pdf_params_chnk(lchnk)%varnce_rt_1, pdf_params_chnk(lchnk)%varnce_rt_2,  &
    !$acc              pdf_params_chnk(lchnk)%thl_1, pdf_params_chnk(lchnk)%thl_2, &
    !$acc              pdf_params_chnk(lchnk)%varnce_thl_1, pdf_params_chnk(lchnk)%varnce_thl_2, &
    !$acc              pdf_params_chnk(lchnk)%corr_w_rt_1, pdf_params_chnk(lchnk)%corr_w_rt_2,  &
    !$acc              pdf_params_chnk(lchnk)%corr_w_thl_1, pdf_params_chnk(lchnk)%corr_w_thl_2, &
    !$acc              pdf_params_chnk(lchnk)%corr_rt_thl_1, pdf_params_chnk(lchnk)%corr_rt_thl_2,&
    !$acc              pdf_params_chnk(lchnk)%alpha_thl, pdf_params_chnk(lchnk)%alpha_rt, &
    !$acc              pdf_params_chnk(lchnk)%crt_1, pdf_params_chnk(lchnk)%crt_2, pdf_params_chnk(lchnk)%cthl_1, &
    !$acc              pdf_params_chnk(lchnk)%cthl_2, pdf_params_chnk(lchnk)%chi_1, &
    !$acc              pdf_params_chnk(lchnk)%chi_2, pdf_params_chnk(lchnk)%stdev_chi_1, &
    !$acc              pdf_params_chnk(lchnk)%stdev_chi_2, pdf_params_chnk(lchnk)%stdev_eta_1, &
    !$acc              pdf_params_chnk(lchnk)%stdev_eta_2, pdf_params_chnk(lchnk)%covar_chi_eta_1, &
    !$acc              pdf_params_chnk(lchnk)%covar_chi_eta_2, pdf_params_chnk(lchnk)%corr_w_chi_1, &
    !$acc              pdf_params_chnk(lchnk)%corr_w_chi_2, pdf_params_chnk(lchnk)%corr_w_eta_1, &
    !$acc              pdf_params_chnk(lchnk)%corr_w_eta_2, pdf_params_chnk(lchnk)%corr_chi_eta_1, &
    !$acc              pdf_params_chnk(lchnk)%corr_chi_eta_2, pdf_params_chnk(lchnk)%rsatl_1, &
    !$acc              pdf_params_chnk(lchnk)%rsatl_2, pdf_params_chnk(lchnk)%rc_1, pdf_params_chnk(lchnk)%rc_2, &
    !$acc              pdf_params_chnk(lchnk)%cloud_frac_1, pdf_params_chnk(lchnk)%cloud_frac_2,  &
    !$acc              pdf_params_chnk(lchnk)%mixt_frac, pdf_params_chnk(lchnk)%ice_supersat_frac_1, &
    !$acc              pdf_params_chnk(lchnk)%ice_supersat_frac_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%rt_1, pdf_params_zm_chnk(lchnk)%rt_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%varnce_rt_1, pdf_params_zm_chnk(lchnk)%varnce_rt_2,  &
    !$acc              pdf_params_zm_chnk(lchnk)%thl_1, pdf_params_zm_chnk(lchnk)%thl_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%varnce_thl_1, pdf_params_zm_chnk(lchnk)%varnce_thl_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_w_rt_1, pdf_params_zm_chnk(lchnk)%corr_w_rt_2,  &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_w_thl_1, pdf_params_zm_chnk(lchnk)%corr_w_thl_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_rt_thl_1, pdf_params_zm_chnk(lchnk)%corr_rt_thl_2,&
    !$acc              pdf_params_zm_chnk(lchnk)%alpha_thl, pdf_params_zm_chnk(lchnk)%alpha_rt, &
    !$acc              pdf_params_zm_chnk(lchnk)%crt_1, pdf_params_zm_chnk(lchnk)%crt_2, pdf_params_zm_chnk(lchnk)%cthl_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%cthl_2, pdf_params_zm_chnk(lchnk)%chi_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%chi_2, pdf_params_zm_chnk(lchnk)%stdev_chi_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%stdev_chi_2, pdf_params_zm_chnk(lchnk)%stdev_eta_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%stdev_eta_2, pdf_params_zm_chnk(lchnk)%covar_chi_eta_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%covar_chi_eta_2, pdf_params_zm_chnk(lchnk)%corr_w_chi_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_w_chi_2, pdf_params_zm_chnk(lchnk)%corr_w_eta_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_w_eta_2, pdf_params_zm_chnk(lchnk)%corr_chi_eta_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%corr_chi_eta_2, pdf_params_zm_chnk(lchnk)%rsatl_1, &
    !$acc              pdf_params_zm_chnk(lchnk)%rsatl_2, pdf_params_zm_chnk(lchnk)%rc_1, pdf_params_zm_chnk(lchnk)%rc_2, &
    !$acc              pdf_params_zm_chnk(lchnk)%cloud_frac_1, pdf_params_zm_chnk(lchnk)%cloud_frac_2,  &
    !$acc              pdf_params_zm_chnk(lchnk)%ice_supersat_frac_1, pdf_params_zm_chnk(lchnk)%ice_supersat_frac_2 )

    !$acc data if( sclr_dim > 0 ) &
    !$acc      create( wpsclrp_sfc, sclrm_forcing, sclrm, wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp_inout) &
    !$acc      copyin( sclr_tol )

    !$acc data if( edsclr_dim > 0 ) &
    !$acc      create( wpedsclrp_sfc, edsclrm_forcing ) &
    !$acc     copy( edsclr_in )

    !$acc data if( hydromet_dim > 0 ) &
    !$acc      create( wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt ) &
    !$acc      copyin( hm_metadata, hm_metadata%l_mix_rat_hm )
    call t_stopf('clubb_tend_cam:acc_copyin')
    call t_startf('clubb_tend_cam:ACCR')

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, pver
      do i = 1, pcols
        temp2d(i,k)      = 0._r8
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
        !  forcings already applied through CAM
        thlm_forcing(i,k) = 0._r8
        rtm_forcing(i,k) = 0._r8
        um_forcing(i,k) = 0._r8
        vm_forcing(i,k) = 0._r8

        rtm_ref(i,k) = 0.0_r8
        thlm_ref(i,k) = 0.0_r8
        um_ref(i,k) = 0.0_r8
        vm_ref(i,k) = 0.0_r8
        ug(i,k) = 0.0_r8
        vg(i,k) = 0.0_r8
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol
        ! Perturbed winds are not used in CAM
        um_pert_inout(i,k) = 0.0_r8
        vm_pert_inout(i,k) = 0.0_r8
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        ! Perturbed winds are not used in CAM
        upwp_pert_inout(i,k) = 0.0_r8
        vpwp_pert_inout(i,k) = 0.0_r8
      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol
      ! Perturbed winds are not used in CAM
      upwp_sfc_pert(i) = 0.0_r8
      vpwp_sfc_pert(i) = 0.0_r8

      ! Determine Coriolis force at given latitude. This is never used
      ! when CLUBB is implemented in a host model, therefore just set
      ! to zero.
      fcor(i) = 0._r8
    end do

    if ( sclr_dim > 0 ) then
      !  higher order scalar stuff, put to zero
      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, nzt_clubb
          do i=1, ncol
            sclrm(i,k,sclr)           = 0._r8
            sclrp3(i,k,sclr)          = 0._r8
            sclrm_forcing(i,k,sclr)   = 0._r8
          end do
        end do
      end do

      !  higher order scalar stuff, put to zero
      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, nzm_clubb
          do i=1, ncol
            wpsclrp(i,k,sclr)         = 0._r8
            sclrp2(i,k,sclr)          = 0._r8
            sclrprtp(i,k,sclr)        = 0._r8
            sclrpthlp(i,k,sclr)       = 0._r8
            sclrpthvp_inout(i,k,sclr) = 0._r8
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do sclr = 1, sclr_dim
        do i=1, ncol
          wpsclrp_sfc(i,sclr) = 0._r8
        end do
      end do
    end if

    if ( edsclr_dim > 0 ) then
      !$acc parallel loop gang vector collapse(3) default(present)
      do edsclr = 1, edsclr_dim
        do i = 1, ncol
          do k = 1, nzt_clubb
            edsclrm_forcing(i,k,edsclr) = 0._r8
          end do
        end do
      end do

      !  Define surface sources for transported variables for diffusion, will
      !  be zero as these tendencies are done in vertical_diffusion
      !$acc parallel loop gang vector collapse(2) default(present)
      do edsclr = 1, edsclr_dim
        do i = 1, ncol
          wpedsclrp_sfc(i,edsclr) = 0._r8
        end do
      end do
    end if

    if ( hydromet_dim > 0 ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do ixind=1, hydromet_dim
        do k=1, nzt_clubb
          do i=1, ncol
            wp2hmp(i,k,ixind)      = 0._r8
            rtphmp_zt(i,k,ixind)   = 0._r8
            thlphmp_zt(i,k,ixind)  = 0._r8
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(3) default(present)
      do ixind=1, hydromet_dim
        do k=1, nzm_clubb
          do i=1, ncol
            wphydrometp(i,k,ixind) = 0._r8
          end do
        end do
      end do

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

      do i = 1, ncol
        do k = 1, pver
          stend(i,k)    = 0._r8
          qvtend(i,k)   = 0._r8
          qitend(i,k)   = 0._r8
          initend(i,k)  = 0._r8
        end do
      end do

      call t_startf('clubb_tend_cam:ice_macro_tend')
      call ice_macro_tend(naai_pbuf(1:ncol,top_lev:pver), state1%t(1:ncol,top_lev:pver),                       &
                          state1%pmid(1:ncol,top_lev:pver), state1%q(1:ncol,top_lev:pver,1),              &
                          state1%q(1:ncol,top_lev:pver,ixcldice), state1%q(1:ncol,top_lev:pver,ixnumice), &
                          latsub, hdtime, stend(1:ncol,top_lev:pver), qvtend(1:ncol,top_lev:pver),        &
                          qitend(1:ncol,top_lev:pver), initend(1:ncol,top_lev:pver), ncol * nzt_clubb )
      call t_stopf('clubb_tend_cam:ice_macro_tend')

      ! update local copy of state with the tendencies
      do i = 1, ncol
        do k = top_lev, pver
          ptend_loc%q(i,k,1)         = qvtend(i,k)
          ptend_loc%q(i,k,ixcldice)  = qitend(i,k)
          ptend_loc%q(i,k,ixnumice)  = initend(i,k)
          ptend_loc%s(i,k)           = stend(i,k)
        end do
      end do

      ! Add the ice tendency to the output tendency
      call physics_ptend_sum(ptend_loc, ptend_all, ncol)

      ! ptend_loc is reset to zero by this call
      call physics_update(state1, ptend_loc, hdtime)

      ! Write output for tendencies:
      do i = 1, ncol
        do k = 1, pver
          temp2d(i,k) =  stend(i,k) / cpairv(i,k,lchnk)
        end do
      end do

      call outfld( 'TTENDICE',  temp2d, pcols, lchnk )
      call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
      call outfld( 'QITENDICE', qitend, pcols, lchnk )
      call outfld( 'NITENDICE', initend, pcols, lchnk )

    endif

    if ( clubb_do_adv ) then

      if (macmic_it  ==  1) then

        !  Note that some of the moments below can be positive or negative.
        !    Remove a constant that was added to prevent dynamics from clipping
        !    them to prevent dynamics from making them positive.
        do k = 1, nzm_clubb-1
          do i = 1, ncol
            k_cam = top_lev - 1 + k  
            thlp2_pbuf(i,k)   = state1%q(i,k_cam,  ixthlp2)
            rtp2_pbuf(i,k)    = state1%q(i,k_cam,   ixrtp2)
            rtpthlp_pbuf(i,k) = state1%q(i,k_cam,ixrtpthlp) - ( rtpthlp_const * apply_const )
            wpthlp_pbuf(i,k)  = state1%q(i,k_cam, ixwpthlp) - ( wpthlp_const  * apply_const )
            wprtp_pbuf(i,k)   = state1%q(i,k_cam,  ixwprtp) - ( wprtp_const   * apply_const )
            wp2_pbuf(i,k)     = state1%q(i,k_cam,    ixwp2)
            wp3_pbuf(i,k)     = state1%q(i,k_cam,    ixwp3) - ( wp3_const     * apply_const )
            up2_pbuf(i,k)     = state1%q(i,k_cam,    ixup2)
            vp2_pbuf(i,k)     = state1%q(i,k_cam,    ixvp2)
          enddo
        enddo

      endif

      ! If not last step of macmic loop then set apply_const back to
      !   zero to prevent output from being corrupted.
      if (macmic_it  ==  cld_macmic_num_steps) then
        apply_const = 1._r8
      else
        apply_const = 0._r8
      endif

      do i = 1, ncol
        thlp2_pbuf(i,nzm_clubb)   = thlp2_pbuf(i,nzm_clubb-1)
        rtp2_pbuf(i,nzm_clubb)    = rtp2_pbuf(i,nzm_clubb-1)
        rtpthlp_pbuf(i,nzm_clubb) = rtpthlp_pbuf(i,nzm_clubb-1)
        wpthlp_pbuf(i,nzm_clubb)  = wpthlp_pbuf(i,nzm_clubb-1)
        wprtp_pbuf(i,nzm_clubb)   = wprtp_pbuf(i,nzm_clubb-1)
        wp2_pbuf(i,nzm_clubb)     = wp2_pbuf(i,nzm_clubb-1)
        up2_pbuf(i,nzm_clubb)     = up2_pbuf(i,nzm_clubb-1)
        vp2_pbuf(i,nzm_clubb)     = vp2_pbuf(i,nzm_clubb-1)
      end do

    endif

    !$acc parallel loop gang vector collapse(2) default(present)
    do n = 1, nparams
      do i = 1, ncol
        clubb_params(i,n) = clubb_params_single_col(1,n)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1, pver
      do i=1, ncol

        !  Compute inverse exner function consistent with CLUBB's definition, which uses a constant
        !  surface pressure.  CAM's exner (in state) does not.  Therefore, for consistent
        !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
        !  (such as thlm), use "inv_exner_clubb" otherwise use the exner in state
        inv_exner_clubb(i,k) = 1._r8 / ( ( state1%pmid(i,k) / p0_clubb )**( rairv(i,k,lchnk) / cpairv(i,k,lchnk) ) )

        !  Compute virtual potential temperature, which is needed for CLUBB
        thv(i,k) = state1%t(i,k) * inv_exner_clubb(i,k) &
                   * ( 1._r8 + zvir * state1%q(i,k,ixq) - state1%q(i,k,ixcldliq) )

        dz_g(i,k) = state1%zi(i,k) - state1%zi(i,k+1)  ! compute thickness

      enddo
    enddo

    !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
    !  Inputs for the momentum levels are set below setup_clubb core
    ! Flipped grid calcs
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k  

        !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state
        !  from the CAM state
        rtm_pbuf(i,k)  = state1%q(i,k_cam,ixq) + state1%q(i,k_cam,ixcldliq)

        thlm_pbuf(i,k) = ( state1%t(i,k_cam) - ( latvap / cpairv(i,k_cam,lchnk) ) * state1%q(i,k_cam,ixcldliq) ) &
                         * inv_exner_clubb(i,k_cam)

        um_pbuf(i,k)   = state1%u(i,k_cam)
        vm_pbuf(i,k)   = state1%v(i,k_cam)

        ! Define the CLUBB thermodynamic grid (in units of m)
        zt_g(i,k) = state1%zm(i,k_cam) - state1%zi(i,pverp)

        ! base state (dry) variables
        rho_ds_zt(i,k)       = rga * ( state1%pdeldry(i,k_cam) / dz_g(i,k_cam) )
        invrs_rho_ds_zt(i,k) = 1._r8 / rho_ds_zt(i,k)

        ! full state (moist) variables
        p_in_Pa(i,k)         = state1%pmid(i,k_cam)
        exner(i,k)           = 1._r8/inv_exner_clubb(i,k_cam)
        thv(i,k)             = state1%t(i,k_cam) * inv_exner_clubb(i,k_cam)  &
                               * (1._r8 + zvir * state1%q(i,k_cam,ixq) - state1%q(i,k_cam,ixcldliq))
        rho_zt(i,k)          = rga*state1%pdel(i,k_cam)/dz_g(i,k_cam)

        ! exception - setting this to moist thv
        thv_ds_zt(i,k)       = thv(i,k)

        rfrzm(i,k)           = state1%q(i,k_cam,ixcldice)

        !  Compute mean w wind on thermo grid, convert from omega to w
        wm_zt(i,k) = -1._r8*(state1%omega(i,k_cam)-state1%omega(i,pver))/(rho_zt(i,k)*gravit)
      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol

      deltaz(i) = state1%zi(i,pverp-1) - state1%zi(i,pverp)

      ! For consistency, set surface pressure to p_in_Pa at the first thermo.
      ! level above the surface (according to the CLUBB ascending grid).
      p_sfc(i) = state1%pmid(i,pver)

      !  Set the elevation of the surface
      sfc_elevation(i) = state1%zi(i,pverp)

    end do



    ! Define the CLUBB momentum grid (in height, units of m)
    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1, nzm_clubb
      do i=1, ncol
        k_cam = top_lev - 1 + k  
        zi_g(i,k) = state1%zi(i,k_cam) - state1%zi(i,pverp)
      end do
    end do


    !  Heights need to be set at each timestep.  Therefore, recall
    !  setup_grid and calc_derrived_params for this.
    !  IMPORTANT NOTE:  do not make any calls that use CLUBB grid-height
    !                   operators (such as zt2zm_api, etc.) until AFTER the
    !                   call to setup_grid_heights_api.
    call t_stopf('clubb_tend_cam:ACCR')
    call t_startf('clubb_tend_cam:NAR')
    !$acc update host( deltaz, zi_g, zt_g, clubb_params, sfc_elevation )

    ! Calculate grid assuming a descending grid (cam grid), since we want to
    ! confine ascending behavior to advance_clubb_core
    call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,   & ! intent(in)
                         .false., grid_type,                     & ! intent(in)
                         deltaz, zi_g(:,nzm_clubb), zi_g(:,1),            & ! intent(in)
                         zi_g, zt_g,                                      & ! intent(in)
                         gr, err_info )                                     ! intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       call endrun(subr//':  '//err_info%err_header_global//NEW_LINE('a')// &
                   'in CLUBB setup_grid')
    end if

    call calc_derrived_params_api( gr, ncol, grid_type, deltaz,                 & ! Intent(in)
                                   clubb_params,                                & ! Intent(in)
                                   clubb_config_flags%l_prescribed_avg_deltaz,  & ! Intent(in)
                                   nu_vert_res_dep, lmin,                       & ! intent(inout)
                                   mixt_frac_max_mag )                            ! intent(inout)

    call check_parameters_api( ncol, clubb_params, lmin, & ! Intent(in)
                               err_info )                  ! Intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       call endrun(subr//': '//err_info%err_header_global//NEW_LINE('a')// &
                   'in CLUBB check_parameters_api')
    end if

    ! CLUBB's grid data structure (gr) and nu_vert_res_dep contain arrays that need to
    ! be copied to the GPU
    call t_stopf('clubb_tend_cam:NAR')
    call t_startf('clubb_tend_cam:acc_copyin')
    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6)
    call t_stopf('clubb_tend_cam:acc_copyin')
    call t_startf('clubb_tend_cam:ACCR')
    
#ifdef SILHS
    ! Add forcings for SILHS covariance contributions
    rtp2_forcing    = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,    rtp2_mc_zt_pbuf(1:ncol,:) )
    thlp2_forcing   = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,   thlp2_mc_zt_pbuf(1:ncol,:) )
    wprtp_forcing   = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,   wprtp_mc_zt_pbuf(1:ncol,:) )
    wpthlp_forcing  = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,  wpthlp_mc_zt_pbuf(1:ncol,:) )
    rtpthlp_forcing = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rtpthlp_mc_zt_pbuf(1:ncol,:) )

    ! Zero out SILHS covariance contribution terms
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, pcols
        rtp2_mc_zt_pbuf(i,k)     = 0.0_r8
        thlp2_mc_zt_pbuf(i,k)    = 0.0_r8
        wprtp_mc_zt_pbuf(i,k)    = 0.0_r8
        wpthlp_mc_zt_pbuf(i,k)   = 0.0_r8
        rtpthlp_mc_zt_pbuf(i,k)  = 0.0_r8
      end do
    end do
#else
    ! Set forcings to zero if not using SILHS
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        rtp2_forcing(i,k)    = 0._r8
        thlp2_forcing(i,k)   = 0._r8
        wprtp_forcing(i,k)   = 0._r8
        wpthlp_forcing(i,k)  = 0._r8
        rtpthlp_forcing(i,k) = 0._r8
      end do
    end do
#endif

    ! Compute some inputs from the thermodynamic grid to the momentum grid
    rho_ds_zm       = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rho_ds_zt )
    invrs_rho_ds_zm = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, invrs_rho_ds_zt )
    rho_zm          = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rho_zt )
    thv_ds_zm       = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, thv_ds_zt )
    wm_zm           = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, wm_zt )

    ! Surface fluxes provided by host model
    !$acc parallel loop gang vector default(present)
    do i=1,ncol
      wpthlp_sfc(i) = cam_in%shf(i) / ( cpairv(i,pver,lchnk) * rho_ds_zm(i,nzm_clubb) ) ! Sensible heat flux
      wpthlp_sfc(i) = wpthlp_sfc(i) * inv_exner_clubb(i,pver)                          ! Potential temperature flux
      wprtp_sfc(i)  = cam_in%cflx(i,1) / rho_ds_zm(i,nzm_clubb)                         ! Moisture flux
    end do


    ! ------------------------------------------------- !
    ! Begin case specific code for SCAM cases.          !
    ! This section of code block is NOT called in       !
    ! global simulations                                !
    ! ------------------------------------------------- !
    if (single_column .and. .not. scm_cambfb_mode) then

      !  Initialize zo if variable ustar is used
      if (cam_in%landfrac(1) >= 0.5_r8) then
        zo = 0.035_r8
      else
        zo = 0.0001_r8
      endif

      !  Compute surface wind (ubar)
      ubar = sqrt(um_pbuf(1,nzt_clubb)**2+vm_pbuf(1,nzt_clubb)**2)
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

          bflx22 = (gravit/theta0)*wpthlp_sfc(1)
          ustar  = diag_ustar(zt_g(1,nzt_clubb-1),bflx22,ubar,zo)
      endif

      !  Compute the surface momentum fluxes, if this is a SCAM simulation
      upwp_sfc(1) = -um_pbuf(1,nzt_clubb)*ustar**2/ubar
      vpwp_sfc(1) = -vm_pbuf(1,nzt_clubb)*ustar**2/ubar

    end if
    
    ! Implementation after Thomas Toniazzo (NorESM) and Colin Zarzycki (PSU)
    !  Other Surface fluxes provided by host model
    if( (cld_macmic_num_steps > 1) .and. clubb_l_intr_sfc_flux_smooth ) then

      call t_stopf('clubb_tend_cam:ACCR')
      call t_startf('clubb_tend_cam:NAR')
      !$acc update host( state1%u, state1%v, state1%t, state1%pmid, cam_in%wsx, cam_in%wsy, rrho )

      ! Adjust surface stresses using winds from the prior macmic iteration
      do i=1,ncol
        ubar = sqrt(state1%u(i,pver)**2+state1%v(i,pver)**2)
        if (ubar <  0.25_r8) ubar = 0.25_r8

        rrho(i) = calc_ideal_gas_rrho(rair, state1%t(i,pver), state1%pmid(i,pver))
        ustar   = calc_friction_velocity(cam_in%wsx(i), cam_in%wsy(i), rrho(i))

        upwp_sfc(i) = -state1%u(i,pver)*ustar**2/ubar
        vpwp_sfc(i) = -state1%v(i,pver)*ustar**2/ubar
      end do

      !$acc update device( upwp_sfc, vpwp_sfc )
      call t_stopf('clubb_tend_cam:NAR')
      call t_startf('clubb_tend_cam:ACCR')

    else

      !$acc parallel loop gang vector default(present)
      do i=1,ncol
        upwp_sfc(i)   = cam_in%wsx(i) / rho_ds_zm(i,nzm_clubb)               ! Surface meridional momentum flux
        vpwp_sfc(i)   = cam_in%wsy(i) / rho_ds_zm(i,nzm_clubb)               ! Surface zonal momentum flux
      end do

    endif

    call t_startf('clubb_tend_cam:flip-index')

    !  Need to flip zt arrays around for CLUBB core
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol
        k_cam = top_lev - 1 + k
        cloud_frac_inout(i,k) = cloud_frac_pbuf(i,k_cam)
        pre_in(i,k)           = prer_evap_pbuf(i,k_cam)
        rcm_pbuf(i,k)         = state1%q(i,k_cam,ixcldliq)
        rvm_in(i,k)           = state1%q(i,k_cam,ixq)
      end do
    end do

    ! We only need to copy pdf_params from pbuf if this is a restart, we're calling pdf_closure 
    ! at the end of advance_clubb_core, and calling it twice for pdf_params_zm as well
    if ( is_first_restart_step() &
         .and. clubb_config_flags%l_call_pdf_closure_twice &
         .and. clubb_config_flags%ipdf_call_placement .eq. ipdf_post_advance_fields ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm_clubb
        do i = 1, ncol
          pdf_params_zm_chnk(lchnk)%w_1(i,k)        = pdf_zm_w_1_pbuf(i,k)
          pdf_params_zm_chnk(lchnk)%w_2(i,k)        = pdf_zm_w_2_pbuf(i,k)
          pdf_params_zm_chnk(lchnk)%varnce_w_1(i,k) = pdf_zm_varnce_w_1_pbuf(i,k)
          pdf_params_zm_chnk(lchnk)%varnce_w_2(i,k) = pdf_zm_varnce_w_2_pbuf(i,k)
          pdf_params_zm_chnk(lchnk)%mixt_frac(i,k)  = pdf_zm_mixt_frac_pbuf(i,k)
        end do
      end do

    end if

    ! pressure,exner on momentum grid needed for mass flux calc.
    if (do_clubb_mf) then

      do k=1,nzt_clubb
        do i=1,ncol
          k_cam = top_lev - 1 + k
          kappa_zt(i,k)       = rairv(i,k_cam,lchnk) / cpairv(i,k_cam,lchnk)
          qc_zt(i,k)          = state1%q(i,k_cam,ixcldliq)
          invrs_exner_zt(i,k) = inv_exner_clubb(i,k_cam)
        end do
      end do

      kappa_zm = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, kappa_zt )

      do k=1,nzm_clubb
        do i=1,ncol
          k_cam = top_lev - 1 + k
          p_in_Pa_zm(i,k)     = state1%pint(i,k_cam)
          invrs_exner_zm(i,k) = 1._r8 / ( (p_in_Pa_zm(i,k)/p0_clubb)**kappa_zm(i,k) )
        end do
      end do

    end if

    if ( clubb_do_adv .and. macmic_it == 1 ) then
      do k = 1, nzm_clubb
        do i = 1, ncol
          thlp2_pbuf(i,k) = max( thl_tol**2, thlp2_pbuf(i,k) )
          rtp2_pbuf(i,k)  = max(  rt_tol**2,  rtp2_pbuf(i,k) )
          wp2_pbuf(i,k)   = max(  w_tol_sqd,   wp2_pbuf(i,k) )
          up2_pbuf(i,k)   = max(  w_tol_sqd,   up2_pbuf(i,k) )
          vp2_pbuf(i,k)   = max(  w_tol_sqd,   vp2_pbuf(i,k) )
        end do
      end do
    end if

    if ( edsclr_dim > 0 ) then

      !  Do the same for tracers
      icnt=0
      do ixind=1,pcnst
        if (lq(ixind))  then

          icnt = icnt+1

          !$acc parallel loop gang vector collapse(2) default(present)
          do k=1,nzt_clubb
            do i=1,ncol
              k_cam = top_lev - 1 + k
              edsclr_in(i,k,icnt) = state1%q(i,k_cam,ixind)
            end do
          end do

        end if
      end do

      if (clubb_l_do_expldiff_rtm_thlm) then

        !$acc parallel loop gang vector collapse(2) default(present)
        do k=1,nzt_clubb
          do i=1, ncol
            k_cam = top_lev - 1 + k
            edsclr_in(i,k,icnt+1) = thlm_pbuf(i,k)
            edsclr_in(i,k,icnt+2) = rtm_pbuf(i,k)
          end do
        end do

      endif

    end if

    call t_stopf('clubb_tend_cam:flip-index')

    do t=1,nadv    ! do needed number of "sub" timesteps for each CAM step

      !  Increment the statistics then begin stats timestep
      if (stats_metadata%l_stats) then
        call stats_begin_timestep_api( t, stats_nsamp, stats_nout, &
                                       stats_metadata )
      endif

      !#######################################################################
      !###################### CALL MF DIAGNOSTIC PLUMES ######################
      !#######################################################################
      if (do_clubb_mf) then
        call t_startf('clubb_tend_cam:do_clubb_mf')

        do k = 1, nzt_clubb
          do i=1, ncol
            invrs_dzt(i,k)  = 1._r8 / dz_g(i,k)
          end do
        end do

        rtm_zm_in  = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,  rtm_pbuf(:ncol,:) )
        thlm_zm_in = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, thlm_pbuf(:ncol,:) )

        !--------------------------------------- integrate_mf call and flip ---------------------------------------
        ! integrate_mf assumes an ascending grid, which is the opposide of the cam grid that
        ! clubb_intr now mainly uses, so we need to flip the fields before calling integrate_mf
        !
        ! Ideally, integrate_mf would operate in descending mode, then we could remove the flipping.
        ! If the column loop gets pushed into it, we can also avoid the array slicing.

        p_in_Pa         = p_in_Pa(:,nzt_clubb:1:-1)
        invrs_exner_zt  = invrs_exner_zt(:,nzt_clubb:1:-1)
        um_pbuf           = um_pbuf(:,nzt_clubb:1:-1)
        vm_pbuf           = vm_pbuf(:,nzt_clubb:1:-1)
        thlm_pbuf         = thlm_pbuf(:,nzt_clubb:1:-1)
        rtm_pbuf          = rtm_pbuf(:,nzt_clubb:1:-1)

        dz_g(:,top_lev:pver)  = dz_g(:,pver:top_lev:-1)
        thv(:,top_lev:pver)   = thv(:,pver:top_lev:-1)

        ! Flip zm inputs
        zi_g            = zi_g(:,nzm_clubb:1:-1)
        p_in_Pa_zm      = p_in_Pa_zm(:,nzm_clubb:1:-1)
        invrs_exner_zm  = invrs_exner_zm(:,nzm_clubb:1:-1)
        thlm_zm_in      = thlm_zm_in(:,nzm_clubb:1:-1)
        rtm_zm_in       = rtm_zm_in(:,nzm_clubb:1:-1)

        do i=1, ncol
          call integrate_mf( nzm_clubb, nzt_clubb, dz_g(i,:), zi_g(i,:), p_in_Pa_zm(i,:), invrs_exner_zm(i,:), & ! input
                                                                         p_in_Pa(i,:),    invrs_exner_zt(i,:), & ! input
                            um_pbuf(i,:), vm_pbuf(i,:), thlm_pbuf(i,:),    rtm_pbuf(i,:), thv(i,1:nzt_clubb),    & ! input
                                                    thlm_zm_in(i,:), rtm_zm_in(i,:),                  & ! input
                                                    wpthlp_sfc(i), wprtp_sfc(i),  pblh_pbuf(i),            & ! input
                            mf_dry_a(i,:),    mf_moist_a(i,:),                                        & ! output - plume diagnostics
                            mf_dry_w(i,:),    mf_moist_w(i,:),                                        & ! output - plume diagnostics
                            mf_dry_qt(i,:),   mf_moist_qt(i,:),                                       & ! output - plume diagnostics
                            mf_dry_thl(i,:),  mf_moist_thl(i,:),                                      & ! output - plume diagnostics
                            mf_dry_u(i,:),    mf_moist_u(i,:),                                        & ! output - plume diagnostics
                            mf_dry_v(i,:),    mf_moist_v(i,:),                                        & ! output - plume diagnostics
                                              mf_moist_qc(i,:),                                       & ! output - plume diagnostics
                              s_ae(i,:),      s_aw(i,:),                                              & ! output - plume diagnostics
                              s_awthl(i,:),   s_awqt(i,:),                                            & ! output - plume diagnostics
                              s_awql(i,:),    s_awqi(i,:),                                            & ! output - plume diagnostics
                              s_awu(i,:),     s_awv(i,:),                                             & ! output - plume diagnostics
                              mf_thlflx(i,:), mf_qtflx(i,:) )                                 ! output - variables needed for solver
        end do
        
        ! Flip zt inputs back
        p_in_Pa         = p_in_Pa(:,nzt_clubb:1:-1)
        invrs_exner_zt  = invrs_exner_zt(:,nzt_clubb:1:-1)
        um_pbuf           = um_pbuf(:,nzt_clubb:1:-1)
        vm_pbuf           = vm_pbuf(:,nzt_clubb:1:-1)
        thlm_pbuf         = thlm_pbuf(:,nzt_clubb:1:-1)
        rtm_pbuf          = rtm_pbuf(:,nzt_clubb:1:-1)

        dz_g(:,top_lev:pver)  = dz_g(:,pver:top_lev:-1)
        thv(:,top_lev:pver)   = thv(:,pver:top_lev:-1)

        ! Flip zm inputs back
        zi_g            = zi_g(:,nzm_clubb:1:-1)
        p_in_Pa_zm      = p_in_Pa_zm(:,nzm_clubb:1:-1)
        invrs_exner_zm  = invrs_exner_zm(:,nzm_clubb:1:-1)
        thlm_zm_in      = thlm_zm_in(:,nzm_clubb:1:-1)
        rtm_zm_in       = rtm_zm_in(:,nzm_clubb:1:-1)
        
        ! Flip clubb_mf output, since it assumes an ascending grid currently
        mf_dry_a     = mf_dry_a(:,nzm_clubb:1:-1)
        mf_moist_a   = mf_moist_a(:,nzm_clubb:1:-1)
        mf_dry_w     = mf_dry_w(:,nzm_clubb:1:-1)
        mf_moist_w   = mf_moist_w(:,nzm_clubb:1:-1)
        mf_dry_qt    = mf_dry_qt(:,nzm_clubb:1:-1)
        mf_moist_qt  = mf_moist_qt(:,nzm_clubb:1:-1)
        mf_dry_thl   = mf_dry_thl(:,nzm_clubb:1:-1)
        mf_moist_thl = mf_moist_thl(:,nzm_clubb:1:-1)
        mf_dry_u     = mf_dry_u(:,nzm_clubb:1:-1)
        mf_moist_u   = mf_moist_u(:,nzm_clubb:1:-1)
        mf_dry_v     = mf_dry_v(:,nzm_clubb:1:-1)
        mf_moist_v   = mf_moist_v(:,nzm_clubb:1:-1)
        mf_moist_qc  = mf_moist_qc(:,nzm_clubb:1:-1)
        mf_thlflx    = mf_thlflx(:,nzm_clubb:1:-1)
        mf_qtflx     = mf_qtflx(:,nzm_clubb:1:-1)
        s_ae         = s_ae(:,nzm_clubb:1:-1)
        s_aw         = s_aw(:,nzm_clubb:1:-1)
        s_awthl      = s_awthl(:,nzm_clubb:1:-1)
        s_awqt       = s_awqt(:,nzm_clubb:1:-1)
        s_awql       = s_awql(:,nzm_clubb:1:-1)
        s_awqi       = s_awqi(:,nzm_clubb:1:-1)
        s_awu        = s_awu(:,nzm_clubb:1:-1)
        s_awv        = s_awv(:,nzm_clubb:1:-1)
        mf_thlflx    = mf_thlflx(:,nzm_clubb:1:-1)
        mf_qtflx     = mf_qtflx(:,nzm_clubb:1:-1)
        !--------------------------------------- END integrate_mf call and flip ---------------------------------------

        ! pass MF turbulent advection term as CLUBB explicit forcing term
        do k = 1, nzt_clubb
          do i = 1, ncol
            rtm_forcing(i,k)  = rtm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * &
                              ((rho_ds_zm(i,k) * mf_qtflx(i,k)) - (rho_ds_zm(i,k+1) * mf_qtflx(i,k+1)))

            thlm_forcing(i,k) = thlm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dzt(i,k) * &
                               ((rho_ds_zm(i,k) * mf_thlflx(i,k)) - (rho_ds_zm(i,k+1) * mf_thlflx(i,k+1)))
          end do
        end do
        call t_stopf('clubb_tend_cam:do_clubb_mf')

      end if

      !  Advance CLUBB CORE one timestep in the future
      call t_startf('clubb_tend_cam:advance_clubb_core_api')
      
      if ( l_ascending_grid ) then

        ! CLUBB is to be run in ascending mode, which has the surface at k=1, which is 
        ! the opposite of the cam grid that the rest of clubb_intr uses, so
        ! we need to flip the fields (in the vertical dimensions) before calling advance_clubb_core
        !
        ! NOTE: We do not neccesarily flip all arrays, only ones that are used within this
        !       subroutine (clubb_tend_cam). For example, only the pdf_params fields that 
        !       are used within this subroutine (or used in a subroutine we call) need to
        !       be flipped. 

        thlm_forcing              =              thlm_forcing(:,nzt_clubb:1:-1)
        rtm_forcing               =               rtm_forcing(:,nzt_clubb:1:-1)
        um_forcing                =                um_forcing(:,nzt_clubb:1:-1)
        vm_forcing                =                vm_forcing(:,nzt_clubb:1:-1)
        wm_zt                     =                     wm_zt(:,nzt_clubb:1:-1)
        rho_zt                    =                    rho_zt(:,nzt_clubb:1:-1)
        rho_ds_zt                 =                 rho_ds_zt(:,nzt_clubb:1:-1)
        invrs_rho_ds_zt           =           invrs_rho_ds_zt(:,nzt_clubb:1:-1)
        thv_ds_zt                 =                 thv_ds_zt(:,nzt_clubb:1:-1)
        rtm_ref                   =                   rtm_ref(:,nzt_clubb:1:-1)
        thlm_ref                  =                  thlm_ref(:,nzt_clubb:1:-1)
        um_ref                    =                    um_ref(:,nzt_clubb:1:-1)
        vm_ref                    =                    vm_ref(:,nzt_clubb:1:-1)
        ug                        =                        ug(:,nzt_clubb:1:-1)
        vg                        =                        vg(:,nzt_clubb:1:-1)
        p_in_Pa                   =                   p_in_Pa(:,nzt_clubb:1:-1)
        exner                     =                     exner(:,nzt_clubb:1:-1)
        rfrzm                     =                     rfrzm(:,nzt_clubb:1:-1)
        um_pbuf                     =                     um_pbuf(:,nzt_clubb:1:-1)
        vm_pbuf                     =                     vm_pbuf(:,nzt_clubb:1:-1)
        up3_pbuf                    =                    up3_pbuf(:,nzt_clubb:1:-1)
        vp3_pbuf                    =                    vp3_pbuf(:,nzt_clubb:1:-1)
        wp3_pbuf                    =                    wp3_pbuf(:,nzt_clubb:1:-1)
        rtp3_pbuf                   =                   rtp3_pbuf(:,nzt_clubb:1:-1)
        thlp3_pbuf                  =                  thlp3_pbuf(:,nzt_clubb:1:-1)
        rcm_pbuf                 =                 rcm_pbuf(:,nzt_clubb:1:-1)
        cloud_frac_inout          =          cloud_frac_inout(:,nzt_clubb:1:-1)
        wpup2_pbuf               =               wpup2_pbuf(:,nzt_clubb:1:-1)
        wpvp2_pbuf               =               wpvp2_pbuf(:,nzt_clubb:1:-1)
        wp2rtp_pbuf              =              wp2rtp_pbuf(:,nzt_clubb:1:-1)
        wp2thlp_pbuf             =             wp2thlp_pbuf(:,nzt_clubb:1:-1)
        ice_supersat_frac_pbuf   =   ice_supersat_frac_pbuf(:,nzt_clubb:1:-1)
        um_pert_inout             =             um_pert_inout(:,nzt_clubb:1:-1)
        vm_pert_inout             =             vm_pert_inout(:,nzt_clubb:1:-1)
        wp2thvp_pbuf                =                wp2thvp_pbuf(:,nzt_clubb:1:-1)
        rtm_pbuf                    =                    rtm_pbuf(:,nzt_clubb:1:-1)
        thlm_pbuf                   =                   thlm_pbuf(:,nzt_clubb:1:-1)

        wprtp_forcing     =    wprtp_forcing(:,nzm_clubb:1:-1)
        wpthlp_forcing    =   wpthlp_forcing(:,nzm_clubb:1:-1)
        rtp2_forcing      =     rtp2_forcing(:,nzm_clubb:1:-1)
        thlp2_forcing     =    thlp2_forcing(:,nzm_clubb:1:-1)
        rtpthlp_forcing   =  rtpthlp_forcing(:,nzm_clubb:1:-1)
        wm_zm             =            wm_zm(:,nzm_clubb:1:-1)
        rho_zm            =           rho_zm(:,nzm_clubb:1:-1)
        rho_ds_zm         =        rho_ds_zm(:,nzm_clubb:1:-1)
        invrs_rho_ds_zm   =  invrs_rho_ds_zm(:,nzm_clubb:1:-1)
        thv_ds_zm         =        thv_ds_zm(:,nzm_clubb:1:-1)
        upwp_pbuf           =          upwp_pbuf(:,nzm_clubb:1:-1)
        vpwp_pbuf           =          vpwp_pbuf(:,nzm_clubb:1:-1)
        up2_pbuf            =           up2_pbuf(:,nzm_clubb:1:-1)
        vp2_pbuf            =           vp2_pbuf(:,nzm_clubb:1:-1)
        wprtp_pbuf          =         wprtp_pbuf(:,nzm_clubb:1:-1)
        wpthlp_pbuf         =        wpthlp_pbuf(:,nzm_clubb:1:-1)
        wp2_pbuf            =           wp2_pbuf(:,nzm_clubb:1:-1)
        rtp2_pbuf           =          rtp2_pbuf(:,nzm_clubb:1:-1)
        thlp2_pbuf          =         thlp2_pbuf(:,nzm_clubb:1:-1)
        rtpthlp_pbuf        =       rtpthlp_pbuf(:,nzm_clubb:1:-1)
        wpthvp_pbuf         =        wpthvp_pbuf(:,nzm_clubb:1:-1)
        rtpthvp_pbuf        =       rtpthvp_pbuf(:,nzm_clubb:1:-1)
        thlpthvp_pbuf       =      thlpthvp_pbuf(:,nzm_clubb:1:-1)
        uprcp_pbuf       =      uprcp_pbuf(:,nzm_clubb:1:-1)
        vprcp_pbuf       =      vprcp_pbuf(:,nzm_clubb:1:-1)
        rc_coef_zm_pbuf  = rc_coef_zm_pbuf(:,nzm_clubb:1:-1)
        wp4_pbuf         =        wp4_pbuf(:,nzm_clubb:1:-1)
        wp2up2_pbuf      =     wp2up2_pbuf(:,nzm_clubb:1:-1)
        wp2vp2_pbuf      =     wp2vp2_pbuf(:,nzm_clubb:1:-1)
        upwp_pert_inout   =  upwp_pert_inout(:,nzm_clubb:1:-1)
        vpwp_pert_inout   =  vpwp_pert_inout(:,nzm_clubb:1:-1)

        if ( edsclr_dim > 0 ) then
          edsclr_in = edsclr_in(:,nzt_clubb:1:-1,:)
          edsclrm_forcing = edsclrm_forcing(:,nzt_clubb:1:-1,:)
        end if

        if ( sclr_dim > 0 ) then
          
          sclrm_forcing = sclrm_forcing(:,nzt_clubb:1:-1,:)
          sclrm         = sclrm(:,nzt_clubb:1:-1,:)
          sclrp3        = sclrp3(:,nzt_clubb:1:-1,:)

          sclrp2          = sclrp2(:,nzm_clubb:1:-1,:)
          sclrprtp        = sclrprtp(:,nzm_clubb:1:-1,:)
          sclrpthlp       = sclrpthlp(:,nzm_clubb:1:-1,:)
          wpsclrp         = wpsclrp(:,nzm_clubb:1:-1,:)
          sclrpthvp_inout = sclrpthvp_inout(:,nzm_clubb:1:-1,:)
        end if
    
        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid
        ! only because these are need to be stored for restarts
        if ( clubb_config_flags%l_call_pdf_closure_twice ) then
          pdf_params_zm_chnk(lchnk)%w_1        = pdf_params_zm_chnk(lchnk)%w_1       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%w_2        = pdf_params_zm_chnk(lchnk)%w_2       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%varnce_w_1 = pdf_params_zm_chnk(lchnk)%varnce_w_1(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%varnce_w_2 = pdf_params_zm_chnk(lchnk)%varnce_w_2(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%mixt_frac  = pdf_params_zm_chnk(lchnk)%mixt_frac (:,nzm_clubb:1:-1)
        end if
        
        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid 
        ! only for pdfp_rtp2_output calc 
        pdf_params_chnk(lchnk)%mixt_frac    = pdf_params_chnk(lchnk)%mixt_frac(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rt_1         = pdf_params_chnk(lchnk)%rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rt_2         = pdf_params_chnk(lchnk)%rt_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_rt_1  = pdf_params_chnk(lchnk)%varnce_rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_rt_2  = pdf_params_chnk(lchnk)%varnce_rt_2(:,nzt_clubb:1:-1)

        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid 
        ! only for update_xp2_mc_api call
        pdf_params_chnk(lchnk)%w_1          = pdf_params_chnk(lchnk)%w_1         (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%w_2          = pdf_params_chnk(lchnk)%w_2         (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_w_1   = pdf_params_chnk(lchnk)%varnce_w_1  (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_w_2   = pdf_params_chnk(lchnk)%varnce_w_2  (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%thl_1        = pdf_params_chnk(lchnk)%thl_1            (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%thl_2        = pdf_params_chnk(lchnk)%thl_2            (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_thl_1 = pdf_params_chnk(lchnk)%varnce_thl_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_thl_2 = pdf_params_chnk(lchnk)%varnce_thl_2(:,nzt_clubb:1:-1)
  
        ! These are flipped for silhs, which uses a cam grid
        pdf_params_chnk(lchnk)%rc_1                 = pdf_params_chnk(lchnk)%rc_1               (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rc_2                 = pdf_params_chnk(lchnk)%rc_2               (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cloud_frac_1         = pdf_params_chnk(lchnk)%cloud_frac_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cloud_frac_2         = pdf_params_chnk(lchnk)%cloud_frac_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%chi_1                = pdf_params_chnk(lchnk)%chi_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%chi_2                = pdf_params_chnk(lchnk)%chi_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%stdev_chi_1          = pdf_params_chnk(lchnk)%stdev_chi_1        (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%stdev_chi_2          = pdf_params_chnk(lchnk)%stdev_chi_2        (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%crt_1                = pdf_params_chnk(lchnk)%crt_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%crt_2                = pdf_params_chnk(lchnk)%crt_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cthl_1               = pdf_params_chnk(lchnk)%cthl_1             (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cthl_2               = pdf_params_chnk(lchnk)%cthl_2             (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%ice_supersat_frac_1  = pdf_params_chnk(lchnk)%ice_supersat_frac_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%ice_supersat_frac_2  = pdf_params_chnk(lchnk)%ice_supersat_frac_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_chi_eta_1       = pdf_params_chnk(lchnk)%corr_chi_eta_1     (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_chi_eta_2       = pdf_params_chnk(lchnk)%corr_chi_eta_2     (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_w_chi_1         = pdf_params_chnk(lchnk)%corr_w_chi_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_w_chi_2         = pdf_params_chnk(lchnk)%corr_w_chi_2       (:,nzt_clubb:1:-1)
          
        if ( t == 1 ) then

          ! we are in ascending mode, need to calculate ascending grid
          call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,   & ! intent(in)
                              l_ascending_grid, grid_type,                     & ! intent(in)
                              deltaz, zi_g(:,1), zi_g(:,nzm_clubb),            & ! intent(in)
                              zi_g(:,nzm_clubb:1:-1), zt_g(:,nzt_clubb:1:-1),   & ! intent(in)
                              gr_a, err_info )                                     ! intent(inout)
        end if
      else
  
        if ( t == 1 ) then
          ! not in ascending mode, so we calculate gr_a the same as gr
          call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,   & ! intent(in)
                              l_ascending_grid, grid_type,                     & ! intent(in)
                              deltaz, zi_g(:,nzm_clubb), zi_g(:,1),            & ! intent(in)
                              zi_g, zt_g,                                      & ! intent(in)
                              gr_a, err_info )                                     ! intent(inout)
        end if
  
      end if


      call advance_clubb_core_api( gr_a, nzm_clubb, nzt_clubb, ncol, &        ! ins
          l_implemented, dtime, fcor(:ncol), sfc_elevation(:ncol), &
          hydromet_dim, &
          sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &
          thlm_forcing(:ncol,:), rtm_forcing(:ncol,:), um_forcing(:ncol,:), vm_forcing(:ncol,:), &
          sclrm_forcing(:ncol,:,:), edsclrm_forcing(:ncol,:,:), wprtp_forcing(:ncol,:), &
          wpthlp_forcing(:ncol,:), rtp2_forcing(:ncol,:), thlp2_forcing(:ncol,:), &
          rtpthlp_forcing(:ncol,:), wm_zm(:ncol,:), wm_zt(:ncol,:), &
          wpthlp_sfc(:ncol), wprtp_sfc(:ncol), upwp_sfc(:ncol), vpwp_sfc(:ncol), p_sfc(:ncol), &
          wpsclrp_sfc(:ncol,:), wpedsclrp_sfc(:ncol,:), &
          upwp_sfc_pert(:ncol), vpwp_sfc_pert(:ncol), &
          rtm_ref(:ncol,:), thlm_ref(:ncol,:), um_ref(:ncol,:), vm_ref(:ncol,:), ug(:ncol,:), vg(:ncol,:), &
          p_in_Pa(:ncol,:), rho_zm(:ncol,:), rho_zt(:ncol,:), exner(:ncol,:), &
          rho_ds_zm(:ncol,:), rho_ds_zt(:ncol,:), invrs_rho_ds_zm(:ncol,:), &
          invrs_rho_ds_zt(:ncol,:), thv_ds_zm(:ncol,:), thv_ds_zt(:ncol,:), &
          hm_metadata%l_mix_rat_hm, &
          rfrzm(:ncol,:), &
          wphydrometp(:ncol,:,:), wp2hmp(:ncol,:,:), rtphmp_zt(:ncol,:,:), thlphmp_zt(:ncol,:,:), &
          grid_dx(:ncol), grid_dy(:ncol), &
          clubb_params(:ncol,:), nu_vert_res_dep, lmin, &
          mixt_frac_max_mag, theta0, ts_nudge, &
          rtm_min, rtm_nudge_max_altitude, &
          clubb_config_flags, &
          stats_metadata, &
          stats_zt(:ncol), stats_zm(:ncol), stats_sfc(:ncol), &                 ! inouts
          um_pbuf(:ncol,:), vm_pbuf(:ncol,:), upwp_pbuf(:ncol,:), vpwp_pbuf(:ncol,:), up2_pbuf(:ncol,:), vp2_pbuf(:ncol,:), up3_pbuf(:ncol,:), vp3_pbuf(:ncol,:), &
          thlm_pbuf(:ncol,:), rtm_pbuf(:ncol,:), wprtp_pbuf(:ncol,:), wpthlp_pbuf(:ncol,:), &
          wp2_pbuf(:ncol,:), wp3_pbuf(:ncol,:), rtp2_pbuf(:ncol,:), rtp3_pbuf(:ncol,:), thlp2_pbuf(:ncol,:), thlp3_pbuf(:ncol,:), rtpthlp_pbuf(:ncol,:), &
          sclrm(:ncol,:,:), &
          sclrp2(:ncol,:,:), sclrp3(:ncol,:,:), sclrprtp(:ncol,:,:), sclrpthlp(:ncol,:,:), &
          wpsclrp(:ncol,:,:), edsclr_in(:ncol,:,:), err_info, &
          rcm_pbuf(:ncol,:), cloud_frac_inout(:ncol,:), &
          wpthvp_pbuf(:ncol,:), wp2thvp_pbuf(:ncol,:), rtpthvp_pbuf(:ncol,:), thlpthvp_pbuf(:ncol,:), &
          sclrpthvp_inout(:ncol,:,:), &
          wp2rtp_pbuf(:ncol,:), wp2thlp_pbuf(:ncol,:), uprcp_pbuf(:ncol,:), &
          vprcp_pbuf(:ncol,:), rc_coef_zm_pbuf(:ncol,:), &
          wp4_pbuf(:ncol,:), wpup2_pbuf(:ncol,:), wpvp2_pbuf(:ncol,:), &
          wp2up2_pbuf(:ncol,:), wp2vp2_pbuf(:ncol,:), ice_supersat_frac_pbuf(:ncol,:), &
          um_pert_inout(:ncol,:), vm_pert_inout(:ncol,:), upwp_pert_inout(:ncol,:), vpwp_pert_inout(:ncol,:), &
          pdf_params_chnk(lchnk), pdf_params_zm_chnk(lchnk), &
          pdf_implicit_coefs_terms_chnk(lchnk), &
          khzm_out(:ncol,:), khzt_out(:ncol,:), &                                                  ! outs
          qclvar_out(:ncol,:), thlprcp_out(:ncol,:), &
          wprcp_out(:ncol,:), w_up_in_cloud_out(:ncol,:), w_down_in_cloud_out(:ncol,:),  &
          cloudy_updraft_frac_out(:ncol,:), cloudy_downdraft_frac_out(:ncol,:), &
          rcm_in_layer(:ncol,:), cloud_cover_out(:ncol,:), invrs_tau_zm_out(:ncol,:), &
          Lscale(:ncol,:) )
          

      if ( l_ascending_grid ) then

        ! If running in ascending mode, we flip the arrays before calling advance_clubb_core
        ! so we need to flip them back. This section should flip every array that was flipped 
        ! before the advance_clubb_core call.

        thlm_forcing              =              thlm_forcing(:,nzt_clubb:1:-1)
        rtm_forcing               =               rtm_forcing(:,nzt_clubb:1:-1)
        um_forcing                =                um_forcing(:,nzt_clubb:1:-1)
        vm_forcing                =                vm_forcing(:,nzt_clubb:1:-1)
        wm_zt                     =                     wm_zt(:,nzt_clubb:1:-1)
        rho_zt                    =                    rho_zt(:,nzt_clubb:1:-1)
        rho_ds_zt                 =                 rho_ds_zt(:,nzt_clubb:1:-1)
        invrs_rho_ds_zt           =           invrs_rho_ds_zt(:,nzt_clubb:1:-1)
        thv_ds_zt                 =                 thv_ds_zt(:,nzt_clubb:1:-1)
        khzt_out                  =                  khzt_out(:,nzt_clubb:1:-1)
        rtm_ref                   =                   rtm_ref(:,nzt_clubb:1:-1)
        thlm_ref                  =                  thlm_ref(:,nzt_clubb:1:-1)
        um_ref                    =                    um_ref(:,nzt_clubb:1:-1)
        vm_ref                    =                    vm_ref(:,nzt_clubb:1:-1)
        ug                        =                        ug(:,nzt_clubb:1:-1)
        vg                        =                        vg(:,nzt_clubb:1:-1)
        p_in_Pa                   =                   p_in_Pa(:,nzt_clubb:1:-1)
        exner                     =                     exner(:,nzt_clubb:1:-1)
        rfrzm                     =                     rfrzm(:,nzt_clubb:1:-1)
        um_pbuf                     =                     um_pbuf(:,nzt_clubb:1:-1)
        vm_pbuf                     =                     vm_pbuf(:,nzt_clubb:1:-1)
        up3_pbuf                    =                    up3_pbuf(:,nzt_clubb:1:-1)
        vp3_pbuf                    =                    vp3_pbuf(:,nzt_clubb:1:-1)
        wp3_pbuf                    =                    wp3_pbuf(:,nzt_clubb:1:-1)
        rtp3_pbuf                   =                   rtp3_pbuf(:,nzt_clubb:1:-1)
        thlp3_pbuf                  =                  thlp3_pbuf(:,nzt_clubb:1:-1)
        rcm_pbuf                 =                 rcm_pbuf(:,nzt_clubb:1:-1)
        cloud_frac_inout          =          cloud_frac_inout(:,nzt_clubb:1:-1)
        wpup2_pbuf               =               wpup2_pbuf(:,nzt_clubb:1:-1)
        wpvp2_pbuf               =               wpvp2_pbuf(:,nzt_clubb:1:-1)
        wp2rtp_pbuf              =              wp2rtp_pbuf(:,nzt_clubb:1:-1)
        wp2thlp_pbuf             =             wp2thlp_pbuf(:,nzt_clubb:1:-1)
        qclvar_out                =                qclvar_out(:,nzt_clubb:1:-1)
        cloud_cover_out           =           cloud_cover_out(:,nzt_clubb:1:-1)
        w_up_in_cloud_out         =         w_up_in_cloud_out(:,nzt_clubb:1:-1)
        w_down_in_cloud_out       =       w_down_in_cloud_out(:,nzt_clubb:1:-1)
        cloudy_updraft_frac_out   =   cloudy_updraft_frac_out(:,nzt_clubb:1:-1)
        cloudy_downdraft_frac_out = cloudy_downdraft_frac_out(:,nzt_clubb:1:-1)
        rcm_in_layer          =          rcm_in_layer(:,nzt_clubb:1:-1)
        ice_supersat_frac_pbuf   =   ice_supersat_frac_pbuf(:,nzt_clubb:1:-1)
        um_pert_inout             =             um_pert_inout(:,nzt_clubb:1:-1)
        vm_pert_inout             =             vm_pert_inout(:,nzt_clubb:1:-1)
        wp2thvp_pbuf                =                wp2thvp_pbuf(:,nzt_clubb:1:-1)
        rtm_pbuf                    =                    rtm_pbuf(:,nzt_clubb:1:-1)
        thlm_pbuf                   =                   thlm_pbuf(:,nzt_clubb:1:-1)
        Lscale                    =                    Lscale(:,nzt_clubb:1:-1)

        wprtp_forcing     =    wprtp_forcing(:,nzm_clubb:1:-1)
        wpthlp_forcing    =   wpthlp_forcing(:,nzm_clubb:1:-1)
        rtp2_forcing      =     rtp2_forcing(:,nzm_clubb:1:-1)
        thlp2_forcing     =    thlp2_forcing(:,nzm_clubb:1:-1)
        rtpthlp_forcing   =  rtpthlp_forcing(:,nzm_clubb:1:-1)
        wm_zm             =            wm_zm(:,nzm_clubb:1:-1)
        rho_zm            =           rho_zm(:,nzm_clubb:1:-1)
        rho_ds_zm         =        rho_ds_zm(:,nzm_clubb:1:-1)
        invrs_rho_ds_zm   =  invrs_rho_ds_zm(:,nzm_clubb:1:-1)
        thv_ds_zm         =        thv_ds_zm(:,nzm_clubb:1:-1)
        upwp_pbuf           =          upwp_pbuf(:,nzm_clubb:1:-1)
        vpwp_pbuf           =          vpwp_pbuf(:,nzm_clubb:1:-1)
        up2_pbuf            =           up2_pbuf(:,nzm_clubb:1:-1)
        vp2_pbuf            =           vp2_pbuf(:,nzm_clubb:1:-1)
        wprtp_pbuf          =         wprtp_pbuf(:,nzm_clubb:1:-1)
        wpthlp_pbuf         =        wpthlp_pbuf(:,nzm_clubb:1:-1)
        wp2_pbuf            =           wp2_pbuf(:,nzm_clubb:1:-1)
        rtp2_pbuf           =          rtp2_pbuf(:,nzm_clubb:1:-1)
        thlp2_pbuf          =         thlp2_pbuf(:,nzm_clubb:1:-1)
        rtpthlp_pbuf        =       rtpthlp_pbuf(:,nzm_clubb:1:-1)
        wpthvp_pbuf         =        wpthvp_pbuf(:,nzm_clubb:1:-1)
        rtpthvp_pbuf        =       rtpthvp_pbuf(:,nzm_clubb:1:-1)
        thlpthvp_pbuf       =      thlpthvp_pbuf(:,nzm_clubb:1:-1)
        uprcp_pbuf       =      uprcp_pbuf(:,nzm_clubb:1:-1)
        vprcp_pbuf       =      vprcp_pbuf(:,nzm_clubb:1:-1)
        rc_coef_zm_pbuf  = rc_coef_zm_pbuf(:,nzm_clubb:1:-1)
        wp4_pbuf         =        wp4_pbuf(:,nzm_clubb:1:-1)
        wp2up2_pbuf      =     wp2up2_pbuf(:,nzm_clubb:1:-1)
        wp2vp2_pbuf      =     wp2vp2_pbuf(:,nzm_clubb:1:-1)
        upwp_pert_inout   =  upwp_pert_inout(:,nzm_clubb:1:-1)
        vpwp_pert_inout   =  vpwp_pert_inout(:,nzm_clubb:1:-1)
        khzm_out          =         khzm_out(:,nzm_clubb:1:-1)
        thlprcp_out       =      thlprcp_out(:,nzm_clubb:1:-1)
        wprcp_out         =        wprcp_out(:,nzm_clubb:1:-1)
        invrs_tau_zm_out  = invrs_tau_zm_out(:,nzm_clubb:1:-1)

        if ( edsclr_dim > 0 ) then
          edsclr_in = edsclr_in(:,nzt_clubb:1:-1,:)
          edsclrm_forcing = edsclrm_forcing(:,nzt_clubb:1:-1,:)
        end if

        if ( sclr_dim > 0 ) then
          
          sclrm_forcing = sclrm_forcing(:,nzt_clubb:1:-1,:)
          sclrm         = sclrm(:,nzt_clubb:1:-1,:)
          sclrp3        = sclrp3(:,nzt_clubb:1:-1,:)

          sclrp2          = sclrp2(:,nzm_clubb:1:-1,:)
          sclrprtp        = sclrprtp(:,nzm_clubb:1:-1,:)
          sclrpthlp       = sclrpthlp(:,nzm_clubb:1:-1,:)
          wpsclrp         = wpsclrp(:,nzm_clubb:1:-1,:)
          sclrpthvp_inout = sclrpthvp_inout(:,nzm_clubb:1:-1,:)
        end if
    
        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid
        ! only because these are need to be stored for restarts
        if ( clubb_config_flags%l_call_pdf_closure_twice ) then
          pdf_params_zm_chnk(lchnk)%w_1        = pdf_params_zm_chnk(lchnk)%w_1       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%w_2        = pdf_params_zm_chnk(lchnk)%w_2       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%varnce_w_1 = pdf_params_zm_chnk(lchnk)%varnce_w_1(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%varnce_w_2 = pdf_params_zm_chnk(lchnk)%varnce_w_2(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk(lchnk)%mixt_frac  = pdf_params_zm_chnk(lchnk)%mixt_frac (:,nzm_clubb:1:-1)
        end if
        
        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid 
        ! only for pdfp_rtp2_output calc 
        pdf_params_chnk(lchnk)%mixt_frac    = pdf_params_chnk(lchnk)%mixt_frac(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rt_1         = pdf_params_chnk(lchnk)%rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rt_2         = pdf_params_chnk(lchnk)%rt_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_rt_1  = pdf_params_chnk(lchnk)%varnce_rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_rt_2  = pdf_params_chnk(lchnk)%varnce_rt_2(:,nzt_clubb:1:-1)

        ! These are flipped, ensuring these are stored in descending mode, regardless of l_ascending_grid 
        ! only for update_xp2_mc_api call
        pdf_params_chnk(lchnk)%w_1        = pdf_params_chnk(lchnk)%w_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%w_2        = pdf_params_chnk(lchnk)%w_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_w_1 = pdf_params_chnk(lchnk)%varnce_w_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_w_2 = pdf_params_chnk(lchnk)%varnce_w_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%thl_1        = pdf_params_chnk(lchnk)%thl_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%thl_2        = pdf_params_chnk(lchnk)%thl_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_thl_1 = pdf_params_chnk(lchnk)%varnce_thl_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%varnce_thl_2 = pdf_params_chnk(lchnk)%varnce_thl_2(:,nzt_clubb:1:-1)
  
        ! These are flipped for silhs, which uses a cam grid
        pdf_params_chnk(lchnk)%rc_1                 = pdf_params_chnk(lchnk)%rc_1               (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%rc_2                 = pdf_params_chnk(lchnk)%rc_2               (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cloud_frac_1         = pdf_params_chnk(lchnk)%cloud_frac_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cloud_frac_2         = pdf_params_chnk(lchnk)%cloud_frac_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%chi_1                = pdf_params_chnk(lchnk)%chi_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%chi_2                = pdf_params_chnk(lchnk)%chi_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%stdev_chi_1          = pdf_params_chnk(lchnk)%stdev_chi_1        (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%stdev_chi_2          = pdf_params_chnk(lchnk)%stdev_chi_2        (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%crt_1                = pdf_params_chnk(lchnk)%crt_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%crt_2                = pdf_params_chnk(lchnk)%crt_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cthl_1               = pdf_params_chnk(lchnk)%cthl_1             (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%cthl_2               = pdf_params_chnk(lchnk)%cthl_2             (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%ice_supersat_frac_1  = pdf_params_chnk(lchnk)%ice_supersat_frac_1(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%ice_supersat_frac_2  = pdf_params_chnk(lchnk)%ice_supersat_frac_2(:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_chi_eta_1       = pdf_params_chnk(lchnk)%corr_chi_eta_1     (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_chi_eta_2       = pdf_params_chnk(lchnk)%corr_chi_eta_2     (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_w_chi_1         = pdf_params_chnk(lchnk)%corr_w_chi_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk(lchnk)%corr_w_chi_2         = pdf_params_chnk(lchnk)%corr_w_chi_2       (:,nzt_clubb:1:-1)
      end if

      call t_stopf('clubb_tend_cam:advance_clubb_core_api')

      ! Note that CLUBB does not produce an error code specific to any column, and
      ! one value only for the entire chunk
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Fatal error in CLUBB advance_clubb_core: at timestep ", get_nstep()
        call endrun(subr//': '//err_info%err_header_global//NEW_LINE('a')//'Fatal error in CLUBB advance_clubb_core')
      end if

      if ( do_rainturb ) then

        call t_startf('clubb_tend_cam:do_rainturb')

        do k=1,nzt_clubb
          do i=1,ncol
            rvm_in(i,k) = rtm_pbuf(i,k) - rcm_pbuf(i,k)
          end do
        end do

        call update_xp2_mc_api( gr, nzm_clubb, nzt_clubb, ncol, dtime, cloud_frac_inout, &
                                rcm_pbuf(:ncol,:), rvm_in, thlm_pbuf(:ncol,:), wm_zt, &
                                exner, pre_in, pdf_params_chnk(lchnk), &
                                rtp2_mc_out, thlp2_mc_out, &
                                wprtp_mc_out, wpthlp_mc_out, &
                                rtpthlp_mc_out)

        do k=1,nzm_clubb
          do i=1,ncol
            dum1 = (1._r8 - cam_in%landfrac(i))

            ! update turbulent moments based on rain evaporation
            rtp2_pbuf(i,k)   = rtp2_pbuf(i,k)   + clubb_rnevap_effic * dum1 * rtp2_mc_out(i,k)   * dtime
            thlp2_pbuf(i,k)  = thlp2_pbuf(i,k)  + clubb_rnevap_effic * dum1 * thlp2_mc_out(i,k)  * dtime
            wprtp_pbuf(i,k)  = wprtp_pbuf(i,k)  + clubb_rnevap_effic * dum1 * wprtp_mc_out(i,k)  * dtime
            wpthlp_pbuf(i,k) = wpthlp_pbuf(i,k) + clubb_rnevap_effic * dum1 * wpthlp_mc_out(i,k) * dtime
          end do
        end do

        call t_stopf('clubb_tend_cam:do_rainturb')

      end if

      if (do_cldcool) then

        call t_startf('clubb_tend_cam:do_cldcool')

        thlp2_rad(:,:) = 0._r8
        
        do k = 1, nzt_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            qrl_clubb(i,k) = qrl_pbuf(i,k_cam) / ( cpairv(i,k_cam,lchnk) * state1%pdeldry(i,k_cam) )
          end do
        end do

        call calculate_thlp2_rad_api( ncol, nzm_clubb, nzt_clubb, gr, &
                                      rcm_pbuf(:ncol,:), thlprcp_out, qrl_clubb, clubb_params, &
                                      thlp2_rad )

        do k=1,nzm_clubb
          do i=1, ncol
            thlp2_pbuf(i,k) = max( thl_tol**2, thlp2_pbuf(i,k) + thlp2_rad(i,k) * dtime )
          end do
        end do

        call t_stopf('clubb_tend_cam:do_cldcool')

      end if

      !  Check to see if stats should be output, here stats are read into
      !  output arrays to make them conformable to CAM output
      if (stats_metadata%l_stats) then
        call t_startf('clubb_tend_cam:stats_end_timestep_clubb')
        do i=1, ncol
          call stats_end_timestep_clubb(i, stats_zt(i), stats_zm(i), stats_rad_zt(i), stats_rad_zm(i), stats_sfc(i), &
                                        out_zt, out_zm, out_radzt, out_radzm, out_sfc)
        end do
        call t_stopf('clubb_tend_cam:stats_end_timestep_clubb')
      end if

    enddo  ! end time loop

    call t_startf('clubb_tend_cam:flip-index')

    if ( clubb_do_adv .and. macmic_it  ==  cld_macmic_num_steps ) then
      do k=1,nzm_clubb
        do i=1, ncol
          thlp2_pbuf(i,k) = max( thl_tol**2, thlp2_pbuf(i,k) )
          rtp2_pbuf(i,k)  = max(  rt_tol**2,  rtp2_pbuf(i,k) )
          wp2_pbuf(i,k)   = max(  w_tol_sqd,   wp2_pbuf(i,k) )
          up2_pbuf(i,k)   = max(  w_tol_sqd,   up2_pbuf(i,k) )
          vp2_pbuf(i,k)   = max(  w_tol_sqd,   vp2_pbuf(i,k) )
        end do
      end do
    end if

    !  Arrays need to be "flipped" to CAM grid
    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1, nzt_clubb
      do i=1, ncol
        k_cam = top_lev - 1 + k
        cloud_frac_pbuf(i,k_cam)  = cloud_frac_inout(i,k)
        qclvar(i,k_cam)           = min( 1._r8, qclvar_out(i,k) )

      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1, nzm_clubb
      do i=1, ncol
        k_cam = top_lev - 1 + k
        khzm_pbuf(i,k_cam)         = khzm_out(i,k)

        ! pdf_params_zm_chnk is already persistent across calls, but we 
        ! save a pbuf version for restarts
        if ( clubb_config_flags%l_call_pdf_closure_twice ) then
          pdf_zm_w_1_pbuf(i,k)        = pdf_params_zm_chnk(lchnk)%w_1(i,k)
          pdf_zm_w_2_pbuf(i,k)        = pdf_params_zm_chnk(lchnk)%w_2(i,k)
          pdf_zm_varnce_w_1_pbuf(i,k) = pdf_params_zm_chnk(lchnk)%varnce_w_1(i,k)
          pdf_zm_varnce_w_2_pbuf(i,k) = pdf_params_zm_chnk(lchnk)%varnce_w_2(i,k)
          pdf_zm_mixt_frac_pbuf(i,k)  = pdf_params_zm_chnk(lchnk)%mixt_frac(i,k)
        end if
      end do
    end do

    call t_stopf('clubb_tend_cam:flip-index')

    ! Values to use above top_lev, for variables that have not already been
    ! set up there. These are mostly fill values that should not actually be
    ! used in the run, but may end up in diagnostic output.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1, top_lev-1
      do i=1, ncol
        cloud_frac_pbuf(i,k)   = 0._r8
        khzm_pbuf(i,k)         = 0._r8
        qclvar(i,k)            = 2._r8
      end do
    end do

!---- TODO: there seems to a an above top_lev interaction here that changes answers.
!           The error occured because we were zeroing out the [1:top_lev-1] values in
!           in rcm_pbuf (when it still contained those levels), when it should've been
!           set to state1%q(:,:,ixcldliq) and unchanged by clubb.
    ! Compute static energy using CLUBB's variables
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, top_lev-1
      do i = 1, ncol
        clubb_s(i,k) = cpairv(i,k,lchnk) * ( ( state1%t(i,k) - ( latvap / cpairv(i,k,lchnk) ) * state1%q(i,k,ixcldliq) ) &
                                             * inv_exner_clubb(i,k) ) / inv_exner_clubb(i,k) &
                       + latvap * 0._r8 &                           ! error kept for BFBness
                       !+ latvap * state1%q(i,k,ixcldliq) &         ! correct line
                       + gravit * state1%zm(i,k) + state1%phis(i)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = top_lev, pver
      do i = 1, ncol
        k_clubb = k + 1 - top_lev
        clubb_s(i,k) = cpairv(i,k,lchnk) * thlm_pbuf(i,k_clubb) / inv_exner_clubb(i,k) &
                       + latvap * rcm_pbuf(i,k_clubb) &
                       + gravit * state1%zm(i,k) + state1%phis(i)
      end do
    end do
!--------------------------------- END TODO ---------------------------------

    ! Section below is concentrated on energy fixing for conservation.
    !   because CLUBB and CAM's thermodynamic variables are different.

    ! Initialize clubbtop to top_lev, for finding the highlest level CLUBB is
    !  active for informing where to apply the energy fixer.
    !$acc parallel loop gang vector default(present)
    do i=1, ncol
      clubbtop(i) = top_lev
      k_clubb     = clubbtop(i) + 1 - top_lev
      do while ((rtp2_pbuf(i,k_clubb) <= 1.e-15_r8 .and. rcm_pbuf(i,k_clubb)  ==  0._r8) .and. clubbtop(i) <  pver)
        clubbtop(i) = clubbtop(i) + 1
        k_clubb     = clubbtop(i) + 1 - top_lev
      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i=1, ncol

      se_a = 0._r8
      ke_a = 0._r8
      wv_a = 0._r8
      wl_a = 0._r8

      se_b = 0._r8
      ke_b = 0._r8
      wv_b = 0._r8
      wl_b = 0._r8

!---- TODO: there seems to a an above top_lev interaction here that changes answers.
!           The error occured because we were zeroing out the [1:top_lev-1] values in
!           in rcm_pbuf (when it still contained those levels), when it should've been
!           set to state1%q(:,:,ixcldliq) and unchanged by clubb.
      do k = 1, top_lev - 1
        ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
        ! after CLUBB is called.  This is for energy conservation purposes.
        se_a = se_a + clubb_s(i,k)*state1%pdel(i,k)*rga
        ke_a = ke_a + 0.5_r8*(state1%u(i,k)**2+state1%v(i,k)**2)*state1%pdel(i,k)*rga
        wv_a = wv_a + ( state1%q(i,k,ixq) + state1%q(i,k,ixcldliq) ) * state1%pdeldry(i,k) * rga   ! error kept for BFBness
        wl_a = wl_a + 0.0_r8                                                                       ! error kept for BFBness
        !wv_a = wv_a + state1%q(i,k,ixq)*state1%pdeldry(i,k)*rga                                   ! correct way
        !wl_a = wl_a + state1%q(i,k,ixcldliq)*state1%pdeldry(i,k)*rga                              ! correct way
      end do

      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.
      do k = top_lev, pver
        k_clubb     = k + 1 - top_lev
        se_a = se_a + clubb_s(i,k)*state1%pdel(i,k)*rga
        ke_a = ke_a + 0.5_r8*(um_pbuf(i,k_clubb)**2+vm_pbuf(i,k_clubb)**2)*state1%pdel(i,k)*rga
        wv_a = wv_a + (rtm_pbuf(i,k_clubb)-rcm_pbuf(i,k_clubb))*state1%pdeldry(i,k)*rga
        wl_a = wl_a + (rcm_pbuf(i,k_clubb))*state1%pdeldry(i,k)*rga
      end do
!--------------------------------- END TODO ---------------------------------

      ! Based on these integrals, compute the total energy after CLUBB call
      te_a = se_a + ke_a + (latvap+latice) * wv_a + latice * wl_a

      do k=1, pver
        ! Do the same as above, but for before CLUBB was called.
        se_b = se_b + state1%s(i,k)*state1%pdel(i,k)*rga
        ke_b = ke_b + 0.5_r8*(state1%u(i,k)**2+state1%v(i,k)**2)*state1%pdel(i,k)*rga
        wv_b = wv_b + state1%q(i,k,ixq)*state1%pdeldry(i,k)*rga
        wl_b = wl_b + state1%q(i,k,ixcldliq)*state1%pdeldry(i,k)*rga
      end do

      ! Based on these integrals, compute the total energy before CLUBB call
      te_b = se_b + ke_b + (latvap+latice) * wv_b + latice * wl_b

      ! Take into account the surface fluxes of heat and moisture
      !  Use correct qflux from cam_in, not lhf/latvap as was done previously
      te_b = te_b + (cam_in%shf(i)+cam_in%cflx(i,1)*(latvap+latice)) * hdtime

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis(i) = ( te_a - te_b ) / ( state1%pint(i,pverp) - state1%pint(i,clubbtop(i)) )

      eleak(i) = ( te_a - te_b ) * invrs_hdtime

    end do

    ! Fix the total energy coming out of CLUBB so it achieves energy conservation.
    ! Apply this fixer throughout the column evenly, but only at layers where
    ! CLUBB is active.
    !
    ! NOTE: The energy fixer seems to cause the climate to change significantly
    ! when using specified dynamics, so allow this to be turned off via a namelist
    ! variable.
    if (clubb_do_energyfix) then

      !$acc parallel loop gang vector default(present)
      do i=1, ncol

        do k=clubbtop(i),pver
          clubb_s(i,k) = clubb_s(i,k) - se_dis(i)*gravit
        end do
        ! convert to units of +ve [K]
        se_dis(i) = -1._r8*se_dis(i)*gravit/cpairv(i,pver,lchnk)

      end do

    endif


    call t_stopf('clubb_tend_cam:ACCR')

    call t_startf('clubb_tend_cam:acc_copyout')
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    call t_stopf('clubb_tend_cam:acc_copyout')

    call t_startf('clubb_tend_cam:NAR')


    call physics_ptend_init( ptend_loc, state%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq )

!---- TODO: there seems to a an above top_lev interaction here that changes answers.
!           The error occured because we were zeroing out the [1:top_lev-1] values in
!           in rcm_pbuf (when it still contained those levels), when it should've been
!           set to state1%q(:,:,ixcldliq) and unchanged by clubb. Had it been set correctly
!           the ptend_loc%q terms would simplify to zero. I've left the (erroneous?) interaction
!           for now to maintain BFBness
    do k = 1, top_lev-1
      do i=1, ncol
        ptend_loc%u(i,k)          = 0.0_r8
        ptend_loc%v(i,k)          = 0.0_r8
        ptend_loc%q(i,k,ixq)      = (   state1%q(i,k,ixcldliq))     * invrs_hdtime ! error kept for BFBness
        ptend_loc%q(i,k,ixcldliq) = ( - state1%q(i,k,ixcldliq))     * invrs_hdtime ! error kept for BFBness
        ! ptend_loc%q(i,k,ixq)      = 0.0_r8                                       ! correct line
        ! ptend_loc%q(i,k,ixcldliq) = 0.0_r8                                       ! correct line
        ptend_loc%s(i,k)          = (clubb_s(i,k) - state1%s(i,k))  * invrs_hdtime ! Tendency of static energy
      end do
    end do

    do k = top_lev, pver
      do i=1, ncol
        k_clubb     = k + 1 - top_lev
        ptend_loc%u(i,k)          = (um_pbuf(i,k_clubb) - state1%u(i,k))               * invrs_hdtime ! east-west wind
        ptend_loc%v(i,k)          = (vm_pbuf(i,k_clubb) - state1%v(i,k))               * invrs_hdtime ! north-south wind
        ptend_loc%q(i,k,ixq)      = (rtm_pbuf(i,k_clubb) - rcm_pbuf(i,k_clubb)-state1%q(i,k,ixq)) * invrs_hdtime ! water vapor
        ptend_loc%q(i,k,ixcldliq) = (rcm_pbuf(i,k_clubb) - state1%q(i,k,ixcldliq))     * invrs_hdtime ! Tendency of liquid water
        ptend_loc%s(i,k)          = (clubb_s(i,k) - state1%s(i,k))          * invrs_hdtime ! Tendency of static energy
      end do
    end do
!--------------------------------- END TODO ---------------------------------

    do i=1, ncol

      !  Now compute the tendencies of CLUBB to CAM
      rtm_integral_vtend(i) = 0._r8
      rtm_integral_ltend(i) = 0._r8

      do k=1, pver
        rtm_integral_ltend(i) = rtm_integral_ltend(i) + ptend_loc%q(i,k,ixcldliq) * state1%pdel(i,k)
        rtm_integral_vtend(i) = rtm_integral_vtend(i) + ptend_loc%q(i,k,     ixq) * state1%pdel(i,k)
      end do

      rtm_integral_ltend(i) = rtm_integral_ltend(i)/gravit
      rtm_integral_vtend(i) = rtm_integral_vtend(i)/gravit

    end do


    ! need to initialize macmic coupling to zero
    if ( macmic_it == 1 ) then
      ttend_clubb_mc_pbuf(:,:)     = 0._r8
      upwp_clubb_gw_mc_pbuf(:,:)   = 0._r8
      vpwp_clubb_gw_mc_pbuf(:,:)   = 0._r8
      thlp2_clubb_gw_mc_pbuf(:,:)  = 0._r8
      wpthlp_clubb_gw_mc_pbuf(:,:) = 0._r8
    end if

    !  Accumulate vars through macmic subcycle for Gravity Wave parameterization
    ttend_clubb_mc_pbuf    (1:ncol,1:nzt_clubb) =     ttend_clubb_mc_pbuf(1:ncol,1:nzt_clubb) + ptend_loc%s(1:ncol,top_lev:pver) / cpair
    upwp_clubb_gw_mc_pbuf  (1:ncol,1:nzm_clubb) =   upwp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) + upwp_pbuf  (1:ncol,1:nzm_clubb)
    vpwp_clubb_gw_mc_pbuf  (1:ncol,1:nzm_clubb) =   vpwp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) + vpwp_pbuf  (1:ncol,1:nzm_clubb)
    thlp2_clubb_gw_mc_pbuf (1:ncol,1:nzm_clubb) =  thlp2_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) + thlp2_pbuf (1:ncol,1:nzm_clubb)
    wpthlp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) = wpthlp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) + wpthlp_pbuf(1:ncol,1:nzm_clubb)

    ! And average at last macmic step
    if (macmic_it == cld_macmic_num_steps) then
      ttend_clubb_pbuf    (1:ncol,top_lev:pver )  =    ttend_clubb_mc_pbuf(1:ncol,1:nzt_clubb) / REAL(cld_macmic_num_steps,r8)
      upwp_clubb_gw_pbuf  (1:ncol,top_lev:pverp) =   upwp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) / REAL(cld_macmic_num_steps,r8)
      vpwp_clubb_gw_pbuf  (1:ncol,top_lev:pverp) =   vpwp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) / REAL(cld_macmic_num_steps,r8)
      thlp2_clubb_gw_pbuf (1:ncol,top_lev:pverp) =  thlp2_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) / REAL(cld_macmic_num_steps,r8)
      wpthlp_clubb_gw_pbuf(1:ncol,top_lev:pverp) = wpthlp_clubb_gw_mc_pbuf(1:ncol,1:nzm_clubb) / REAL(cld_macmic_num_steps,r8)
    end if

    if (clubb_do_adv) then
      if (macmic_it == cld_macmic_num_steps) then

        ! Zero above top_lev 
        do k = 1, top_lev - 1
          do i=1, ncol
            ptend_loc%q(i,k,ixthlp2)   = 0.0_r8
            ptend_loc%q(i,k,ixrtp2)    = 0.0_r8
            ptend_loc%q(i,k,ixrtpthlp) = 0.0_r8
            ptend_loc%q(i,k,ixwpthlp)  = 0.0_r8
            ptend_loc%q(i,k,ixwprtp)   = 0.0_r8
            ptend_loc%q(i,k,ixwp2)     = 0.0_r8
            ptend_loc%q(i,k,ixwp3)     = 0.0_r8
            ptend_loc%q(i,k,ixup2)     = 0.0_r8
            ptend_loc%q(i,k,ixvp2)     = 0.0_r8

          end do
        end do

        do k = top_lev, pver
          do i=1, ncol

            k_clubb = k + 1 - top_lev

            ! Here add a constant to moments which can be either positive or
            !  negative.  This is to prevent clipping when dynamics tries to
            !  make all constituents positive
            wp3_pbuf(i,k_clubb)     = wp3_pbuf(i,k_clubb)     + wp3_const
            rtpthlp_pbuf(i,k_clubb) = rtpthlp_pbuf(i,k_clubb) + rtpthlp_const
            wpthlp_pbuf(i,k_clubb)  = wpthlp_pbuf(i,k_clubb)  + wpthlp_const
            wprtp_pbuf(i,k_clubb)   = wprtp_pbuf(i,k_clubb)   + wprtp_const

            ptend_loc%q(i,k,ixthlp2)   = (thlp2_pbuf(i,k_clubb)   - state1%q(i,k,ixthlp2))   * invrs_hdtime ! THLP Variance
            ptend_loc%q(i,k,ixrtp2)    = (rtp2_pbuf(i,k_clubb)    - state1%q(i,k,ixrtp2))    * invrs_hdtime ! RTP Variance
            ptend_loc%q(i,k,ixrtpthlp) = (rtpthlp_pbuf(i,k_clubb) - state1%q(i,k,ixrtpthlp)) * invrs_hdtime ! RTP THLP covariance
            ptend_loc%q(i,k,ixwpthlp)  = (wpthlp_pbuf(i,k_clubb)  - state1%q(i,k,ixwpthlp))  * invrs_hdtime ! WPTHLP
            ptend_loc%q(i,k,ixwprtp)   = (wprtp_pbuf(i,k_clubb)   - state1%q(i,k,ixwprtp))   * invrs_hdtime ! WPRTP
            ptend_loc%q(i,k,ixwp2)     = (wp2_pbuf(i,k_clubb)     - state1%q(i,k,ixwp2))     * invrs_hdtime ! WP2
            ptend_loc%q(i,k,ixwp3)     = (wp3_pbuf(i,k_clubb)     - state1%q(i,k,ixwp3))     * invrs_hdtime ! WP3
            ptend_loc%q(i,k,ixup2)     = (up2_pbuf(i,k_clubb)     - state1%q(i,k,ixup2))     * invrs_hdtime ! UP2
            ptend_loc%q(i,k,ixvp2)     = (vp2_pbuf(i,k_clubb)     - state1%q(i,k,ixvp2))     * invrs_hdtime ! VP2

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
    ! edsclr_in is compressed with just the constituents being used, ptend and state are not compressed
    icnt=0
    do ixind=1,pcnst
      if (lq(ixind)) then
        icnt=icnt+1
        if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
            (ixind /= ixthlp2)   .and. (ixind /= ixrtp2)   .and.&
            (ixind /= ixrtpthlp) .and. (ixind /= ixwpthlp) .and.&
            (ixind /= ixwprtp)   .and. (ixind /= ixwp2)    .and.&
            (ixind /= ixwp3)     .and. (ixind /= ixup2)    .and. (ixind /= ixvp2) ) then

          ! Copy CLUBB's edsclr values
          do k=top_lev, pver
            do i=1, ncol
              k_clubb = k + 1 - top_lev
              ptend_loc%q(i,k,ixind) = (edsclr_in(i,k_clubb,icnt)-state1%q(i,k,ixind)) / hdtime ! transported constituents
            end do
          end do

          ! Zero out levels above top_lev
          do k=1, top_lev-1
            do i=1, ncol
              ptend_loc%q(i,k,ixind) = 0._r8
            end do
          end do

        end if
      end if
    end do

    rvmtend_clubb_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    rcmtend_clubb_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    rimtend_clubb_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldice)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    stend_clubb_output(:ncol,:pver) = ptend_loc%s(:ncol,:pver)
    utend_clubb_output(:ncol,:pver) = ptend_loc%u(:ncol,:pver)
    vtend_clubb_output(:ncol,:pver) = ptend_loc%v(:ncol,:pver)
    cmeliq_pbuf(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)

    !
    ! set pbuf field so that HB scheme is only applied above CLUBB top
    !
    if (do_hb_above_clubb) then
      call pbuf_set_field(pbuf, clubbtop_idx, clubbtop)
    endif

    ! ------------------------------------------------- !
    ! End column computation of CLUBB, begin to apply   !
    ! and compute output, etc                           !
    ! ------------------------------------------------- !

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state1,ptend_loc,hdtime)

    ! Due to the order of operation of CLUBB, which closes on liquid first,
    ! then advances it's predictive equations second, this can lead to
    ! RHliq > 1 directly before microphysics is called.  Therefore, we use
    ! ice_macro_tend to enforce RHliq <= 1 everywhere before microphysics is called.

    if (clubb_do_liqsupersat) then

      call t_startf('clubb_cam_tend:do_liqsupersat')
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

      call liquid_macro_tend(npccn_pbuf(1:ncol,top_lev:pver), state1%t(1:ncol,top_lev:pver),                      &
                             state1%pmid(1:ncol,top_lev:pver), state1%q(1:ncol,top_lev:pver,ixq),            &
                             state1%q(1:ncol,top_lev:pver,ixcldliq), state1%q(1:ncol,top_lev:pver,ixnumliq), &
                             latvap, hdtime, stend(1:ncol,top_lev:pver),qvtend(1:ncol,top_lev:pver),         &
                             qctend(1:ncol,top_lev:pver), inctend(1:ncol,top_lev:pver), ncol * nzt_clubb )

      ! update local copy of state with the tendencies
      ptend_loc%q(:ncol,top_lev:pver,ixq)       =  qvtend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixcldliq)  =  qctend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixnumliq)  = inctend(:ncol,top_lev:pver)
      ptend_loc%s(:ncol,top_lev:pver)           =   stend(:ncol,top_lev:pver)

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
        temp2d = 1._r8
      elsewhere
        temp2d = 0._r8
      end where

      call outfld( 'FQTENDICE', temp2d, pcols, lchnk )
      call t_stopf('clubb_cam_tend:do_liqsupersat')
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
    dlf2 = 0.0_r8
    dlf_liq_out(:,:) = 0.0_r8
    dlf_ice_out(:,:) = 0.0_r8

    lqice(:)        = .false.
    lqice(ixcldliq) = .true.
    lqice(ixcldice) = .true.
    lqice(ixnumliq) = .true.
    lqice(ixnumice) = .true.

    dl_rad = clubb_detliq_rad
    di_rad = clubb_detice_rad
    dt_low = clubb_detphase_lowtemp

    call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lq=lqice)

    do k=1,pver
      do i=1,ncol

        if( state1%t(i,k) > meltpt_temp ) then
          dum1 = 0.0_r8
        elseif ( state1%t(i,k) < dt_low ) then
          dum1 = 1.0_r8
        else
          dum1 = ( meltpt_temp - state1%t(i,k) ) / ( meltpt_temp - dt_low )
        endif

        ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
        ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
        ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2 )) * ( 1._r8 - dum1 ) ) &
                                   / (4._r8*3.14_r8*dl_rad**3*997._r8) + & ! Deep    Convection
                                   3._r8 * (                         dlf2    * ( 1._r8 - dum1 ) ) &
                                   / (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection
        ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2 )) *  dum1 ) &
                                   / (4._r8*3.14_r8*di_rad**3*500._r8) + & ! Deep    Convection
                                   3._r8 * (                         dlf2    *  dum1 ) &
                                   / (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
        ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

        dlf_liq_out(i,k) = dlf(i,k) * ( 1._r8 - dum1 )
        dlf_ice_out(i,k) = dlf(i,k) * dum1

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
    dpdlfliq_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldliq)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    dpdlfice_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixcldice)*state1%pdeldry(:ncol,:pver)/state1%pdel(:ncol,:pver)
    dpdlft_output(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpairv(:ncol,:pver, lchnk)
    detnliquid_output(:ncol,:pver) = ptend_loc%q(:ncol,:pver,ixnumliq)

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

    do i = 1, ncol
      do k = 1, pver
        relvar_pbuf(i,k) = relvarmax  ! default
      end do
    end do

    if (deep_scheme .ne. 'CLUBB_SGS') then
      do i = 1, ncol
        do k = top_lev, pver
          k_clubb = k + 1 - top_lev
          if ( rcm_pbuf(i,k_clubb) /= 0 .and. qclvar(i,k) /= 0 ) then
            relvar_pbuf(i,k) = min( relvarmax, max(0.001_r8, rcm_pbuf(i,k_clubb)**2 / qclvar(i,k) ) )
          end if
        end do
      end do
    endif

    ! ------------------------------------------------- !
    ! Optional Accretion enhancement factor             !
    ! ------------------------------------------------- !
    accre_enhan_pbuf(:ncol,:pver) = 1._r8
    
    !  turbulent kinetic energy
    do k =  top_lev, pverp
      do i = 1, ncol
        k_clubb     = k + 1 - top_lev
        tke_pbuf(i,k) = 0.5_r8 * ( up2_pbuf(i,k_clubb) + vp2_pbuf(i,k_clubb) + wp2_pbuf(i,k_clubb) ) 
      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  Diagnose some quantities that are computed in macrop_tend here.                  !
    !  These are inputs required for the microphysics calculation.                      !
    !                                                                                   !
    !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
    ! --------------------------------------------------------------------------------- !

    !  initialize variables
    alst_pbuf(:,:) = 0.0_r8
    qlst_pbuf(:,:) = 0.0_r8

    do k = top_lev, pver
      do i = 1, ncol
        k_clubb     = k + 1 - top_lev
        alst_pbuf(i,k) = cloud_frac_pbuf(i,k)
        qlst_pbuf(i,k) = rcm_pbuf(i,k_clubb) / max( 0.01_r8, alst_pbuf(i,k) )  ! Incloud stratus condensate mixing ratio
      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
    ! --------------------------------------------------------------------------------- !

    frac_limit = 0.01_r8
    ic_limit   = 1.e-12_r8
    deepcu_pbuf(:,:) = 0.0_r8
    shalcu_pbuf(:,:) = 0.0_r8

    do k=1,pver-1
      do i=1,ncol
        !  diagnose the deep convective cloud fraction, as done in macrophysics based on the
        !  deep convective mass flux, read in from pbuf.  Since shallow convection is never
        !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
        !  fraction is purely from deep convection scheme.
        deepcu_pbuf(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(i,k+1)-cmfmc_sh_pbuf(i,k+1))),0.6_r8))
        shalcu_pbuf(i,k) = 0._r8

        if (deepcu_pbuf(i,k) <= frac_limit .or. dp_icwmr_pbuf(i,k) < ic_limit) then
          deepcu_pbuf(i,k) = 0._r8
        endif

        !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable
        !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
        !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud
        !  from CLUBB plus the deep convective cloud fraction
        concld_pbuf(i,k) = min(cloud_frac_pbuf(i,k)-alst_pbuf(i,k)+deepcu_pbuf(i,k),0.80_r8)
      enddo
    enddo

    if (single_column .and. .not. scm_cambfb_mode) then
      if (trim(scm_clubb_iop_name)  ==  'ATEX_48hr'       .or. &
          trim(scm_clubb_iop_name)  ==  'BOMEX_5day'      .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf01_4day' .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf02_06hr' .or. &
          trim(scm_clubb_iop_name)  ==  'RICO_3day'       .or. &
          trim(scm_clubb_iop_name)  ==  'ARM_CC') then

         deepcu_pbuf(:,:) = 0.0_r8
         concld_pbuf(:,:) = 0.0_r8

      endif
    endif

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !

    !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
    troplev(:) = 0
    !REMOVECAM_END
    call tropopause_findChemTrop( state, troplev )

    aist_pbuf(:,:top_lev-1) = 0._r8
    qsatfac_pbuf(:, :) = 0._r8 ! Zero out entire profile in case qsatfac is left undefined in aist_vector below

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
             state1%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist_pbuf(:,k),ncol )
      else
        call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
              state1%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist_pbuf(:,k),ncol,&
              qsatfac_out=qsatfac_pbuf(:,k), rhmini_in=rhmini, rhmaxi_in=rhmaxi)
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
        ast_pbuf(i,k) = max(alst_pbuf(i,k),aist_pbuf(i,k))
        qist_pbuf(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist_pbuf(i,k))
      enddo
    enddo

   !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just
   !  be outputting the shallow convective cloud fraction
    do k=1,pver
      do i=1,ncol
        cloud_frac_pbuf(i,k) = min(ast_pbuf(i,k)+deepcu_pbuf(i,k),1.0_r8)
      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  DIAGNOSE THE PBL DEPTH                                                           !
    !  this is needed for aerosol code                                                  !
    ! --------------------------------------------------------------------------------- !
    do i=1,ncol
      do k=1,pver
         !subroutine pblind expects "Stull" definition of Exner
         th(i,k) = state1%t(i,k)*state1%exner(i,k)
         !thv should have condensate loading to be consistent with earlier def's in this module
         thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq) - state1%q(i,k,ixcldliq))
      enddo
    enddo

    ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
    rrho   (1:ncol) = calc_ideal_gas_rrho(rair, state1%t(1:ncol,pver), state1%pmid(1:ncol,pver))
    ustar2 (1:ncol) = calc_friction_velocity(cam_in%wsx(1:ncol), cam_in%wsy(1:ncol), rrho(1:ncol))
    ! use correct qflux from coupler
    kinheat(1:ncol) = calc_kinematic_heat_flux(cam_in%shf(1:ncol), rrho(1:ncol), cpair)
    kinwat (1:ncol) = calc_kinematic_water_vapor_flux(cam_in%cflx(1:ncol,1), rrho(1:ncol))
    kbfs   (1:ncol) = calc_kinematic_buoyancy_flux(kinheat(1:ncol), zvir, th(1:ncol,pver), kinwat(1:ncol))
    obklen (1:ncol) = calc_obukhov_length(thv(1:ncol,pver), ustar2(1:ncol), gravit, karman, kbfs(1:ncol))


    where (kbfs(:ncol)  ==  -0.0_r8) kbfs(:ncol) = 0.0_r8

    ! Compute PBL depth according to Holtslag-Boville Scheme -- only pblh is needed here
    ! and other outputs are discarded
    !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
    pblh_pbuf(:) = 0._r8
    dummy2(:) = 0._r8
    dummy3(:) = 0._r8
    !REMOVECAM_END
    call hb_pbl_dependent_coefficients_run( &
      ncol      = ncol,                                      &
      pver      = pver,                                      &
      pverp     = pverp,                                     &
      gravit    = gravit,                                    &
      z         = state1%zm(:ncol,:pver),                    &
      zi        = state1%zi(:ncol,:pverp),                   &
      u         = state1%u(:ncol,:pver),                     &
      v         = state1%v(:ncol,:pver),                     &
      cldn      = cloud_frac_pbuf(:ncol,:pver),                   &
      ! Inputs from CLUBB (not HB coefficients)
      thv       = thv(:ncol,:pver),                          &
      ustar     = ustar2(:ncol),                             &
      kbfs      = kbfs(:ncol),                               &
      obklen    = obklen(:ncol),                             &
      ! Output variables
      pblh      = pblh_pbuf(:ncol),                               &
      wstar     = dummy2(:ncol),                             &
      bge       = dummy3(:ncol),                             &
      errmsg    = errmsg,                                    &
      errflg    = errflg)

    ! Assign the first pver levels of cloud_frac back to cld
    cld_pbuf(:,1:pver) = cloud_frac_pbuf(:,1:pver)

    ! --------------------------------------------------------------------------------- !
    !  END CLOUD FRACTION DIAGNOSIS, begin to store variables back into buffer          !
    ! --------------------------------------------------------------------------------- !

    call outfld( 'DETNLIQTND', detnliquid_output,pcols, lchnk )

    !  Output CLUBB tendencies (convert dry basis to wet for consistency with  history variable definition)
    call outfld( 'RVMTEND_CLUBB', rvmtend_clubb_output, pcols, lchnk)
    call outfld( 'RCMTEND_CLUBB', rcmtend_clubb_output, pcols, lchnk)
    call outfld( 'RIMTEND_CLUBB', rimtend_clubb_output, pcols, lchnk)
    call outfld( 'STEND_CLUBB', stend_clubb_output, pcols, lchnk)
    call outfld( 'UTEND_CLUBB', utend_clubb_output, pcols, lchnk)
    call outfld( 'VTEND_CLUBB', vtend_clubb_output, pcols, lchnk)

    call outfld( 'CMELIQ', cmeliq_pbuf, pcols, lchnk)

    ! output moist basis to be consistent with history variable definition
    call outfld( 'DPDLFLIQ', dpdlfliq_output, pcols, lchnk)
    call outfld( 'DPDLFICE', dpdlfice_output, pcols, lchnk)
    call outfld( 'DPDLFT',   dpdlft_output, pcols, lchnk)

    !  Output the PBL depth
    call outfld('PBLH', pblh_pbuf, pcols, lchnk)

    call outfld('KVH_CLUBB', khzm_pbuf, pcols, lchnk)
    call outfld('ELEAK_CLUBB', eleak, pcols, lchnk)
    call outfld('TFIX_CLUBB', se_dis, pcols, lchnk)

    ! ------------------------------------------------- !
    ! Diagnose some output variables                    !
    ! ------------------------------------------------- !

    !  density
    rho(1:ncol,1:pver) = rga*state1%pdel(1:ncol,1:pver)/(state1%zi(1:ncol,1:pver)-state1%zi(1:ncol,2:pverp))
    rho(1:ncol,pverp) = rho(1:ncol,pver)

    do k = top_lev, pverp
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        zi_output(i,k)       = zi_g(i,k_clubb)
        
        wp2_output(i,k) = wp2_pbuf(i,k_clubb)
        up2_output(i,k) = up2_pbuf(i,k_clubb)
        vp2_output(i,k) = vp2_pbuf(i,k_clubb)
        upwp_output(i,k) = upwp_pbuf(i,k_clubb)
        vpwp_output(i,k) = vpwp_pbuf(i,k_clubb)
        rtp2_output(i,k) = rtp2_pbuf(i,k_clubb)
        wprcp_clubb_output(i,k) = wprcp_out(i,k_clubb) * latvap
        wpthvp_clubb_output(i,k) = wpthvp_pbuf(i,k_clubb) * cpair
        thlp2_output(i,k) = thlp2_pbuf(i,k_clubb)

        wpthlp_output(i,k)  = (wpthlp_pbuf(i,k_clubb)-(apply_const*wpthlp_const))*rho(i,k)*cpair !  liquid water potential temperature flux
        wprtp_output(i,k)   = (wprtp_pbuf(i,k_clubb)-(apply_const*wprtp_const))*rho(i,k)*latvap  !  total water mixig ratio flux
        rtpthlp_output(i,k) = rtpthlp_pbuf(i,k_clubb)-(apply_const*rtpthlp_const)                !  rtpthlp output

      end do
    end do

    ! Convert RTP2 and THLP2 to thermo grid for output
    rtp2_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr,  rtp2_pbuf(:ncol,:) )
    thl2_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr, thlp2_pbuf(:ncol,:) )
    wp2_zt  = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr,   wp2_pbuf(:ncol,:) )

    do k = top_lev, pver
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        rcm_output(i,k) = rcm_pbuf(i,k_clubb)
        rtm_output(i,k) = rtm_pbuf(i,k_clubb)
        thlm_output(i,k) = thlm_pbuf(i,k_clubb)
        um_output(i,k) = um_pbuf(i,k_clubb)
        vm_output(i,k) = vm_pbuf(i,k_clubb)

        rcm_in_layer_output(i,k) = rcm_in_layer(i,k_clubb)
        zt_output(i,k)       = zt_g(i,k_clubb)
        wm_zt_output(i,k)    = wm_zt(i,k_clubb)

        rtp2_zt_output(i,k)  = rtp2_zt(i,k_clubb)
        thl2_zt_output(i,k)  = thl2_zt(i,k_clubb)
        wp2_zt_output(i,k)   = wp2_zt(i,k_clubb)

        wp3_output(i,k) = wp3_pbuf(i,k_clubb) - (apply_const*wp3_const)                      !  wp3 output

      end do
    end do

    do k=1, nzt_clubb
      do i=1, ncol

        mean_rt = pdf_params_chnk(lchnk)%mixt_frac(i,k) &
                  * pdf_params_chnk(lchnk)%rt_1(i,k) &
                  + ( 1.0_r8 - pdf_params_chnk(lchnk)%mixt_frac(i,k) ) &
                    * pdf_params_chnk(lchnk)%rt_2(i,k)

        k_cam = top_lev - 1 + k

        pdfp_rtp2_output(i,k_cam) = pdf_params_chnk(lchnk)%mixt_frac(i,k) &
                                * ( ( pdf_params_chnk(lchnk)%rt_1(i,k) - mean_rt )**2 &
                                    + pdf_params_chnk(lchnk)%varnce_rt_1(i,k) ) &
                                + ( 1.0_r8 - pdf_params_chnk(lchnk)%mixt_frac(i,k) ) &
                                  * ( ( pdf_params_chnk(lchnk)%rt_2(i,k) - mean_rt )**2 &
                                      + pdf_params_chnk(lchnk)%varnce_rt_2(i,k) )
      end do
    end do

    do k = 1, top_lev-1
      do i = 1, ncol
        wp2_output(i,k) = 0._r8
        up2_output(i,k) = 0._r8
        vp2_output(i,k) = 0._r8
        rtp2_output(i,k) = 0._r8
        thlp2_output(i,k) = 0._r8
        zt_output(i,k) = 0._r8
        rtp2_zt_output(i,k) = 0._r8
        wp3_output(i,k) = 0._r8
        thl2_zt_output(i,k) = 0._r8
        wp2_zt_output(i,k) = 0._r8
        rcm_in_layer_output(i,k) = 0._r8
        pdfp_rtp2_output(i,k) = 0._r8
        wm_zt_output(i,k) = 0._r8
        rcm_output(i,k) = 0._r8
        rtm_output(i,k) = 0._r8
        thlm_output(i,k) = 0._r8
        um_output(i,k) = 0._r8
        vm_output(i,k) = 0._r8
        zi_output(i,k) = 0._r8
        wpthlp_output(i,k) = 0._r8
        rtpthlp_output(i,k) = 0._r8
        wprtp_output(i,k) = 0._r8
        upwp_output(i,k) = 0._r8
        vpwp_output(i,k) = 0._r8
        wprcp_clubb_output(i,k) = 0._r8
        wpthvp_clubb_output(i,k) = 0._r8
      end do
    end do

    !  Output calls of variables goes here
    call outfld( 'WP2_CLUBB',        wp2_output,                     pcols, lchnk )
    call outfld( 'UP2_CLUBB',        up2_output,                     pcols, lchnk )
    call outfld( 'VP2_CLUBB',        vp2_output,                     pcols, lchnk )
    call outfld( 'WP3_CLUBB',        wp3_output,              pcols, lchnk )
    call outfld( 'UPWP_CLUBB',       upwp_output,                    pcols, lchnk )
    call outfld( 'VPWP_CLUBB',       vpwp_output,                    pcols, lchnk )
    call outfld( 'WPTHLP_CLUBB',     wpthlp_output,           pcols, lchnk )
    call outfld( 'WPRTP_CLUBB',      wprtp_output,            pcols, lchnk )
    call outfld( 'RTP2_CLUBB',       rtp2_output,                    pcols, lchnk )
    call outfld( 'RTPTHLP_CLUBB',    rtpthlp_output,          pcols, lchnk )
    call outfld( 'RCM_CLUBB',        rcm_output,                     pcols, lchnk )
    call outfld( 'RTM_CLUBB',        rtm_output,                     pcols, lchnk )
    call outfld( 'THLM_CLUBB',       thlm_output,                    pcols, lchnk )
    call outfld( 'WPRCP_CLUBB',      wprcp_clubb_output,             pcols, lchnk )
    call outfld( 'WPTHVP_CLUBB',     wpthvp_clubb_output,            pcols, lchnk )
    call outfld( 'RTP2_ZT_CLUBB',    rtp2_zt_output,             pcols, lchnk )
    call outfld( 'THLP2_ZT_CLUBB',   thl2_zt_output,             pcols, lchnk )
    call outfld( 'WP2_ZT_CLUBB',     wp2_zt_output,              pcols, lchnk )
    call outfld( 'pdfp_rtp2_output_CLUBB',  pdfp_rtp2_output,               pcols, lchnk )
    call outfld( 'THLP2_CLUBB',      thlp2_output,                   pcols, lchnk )
    call outfld( 'RCMINLAYER_CLUBB', rcm_in_layer_output,            pcols, lchnk )
    call outfld( 'CLOUDCOVER_CLUBB', cloud_frac_pbuf,              pcols, lchnk )
    call outfld( 'ZT_CLUBB',         zt_output,                  pcols, lchnk )
    call outfld( 'ZM_CLUBB',         zi_output,                  pcols, lchnk )
    call outfld( 'UM_CLUBB',         um_output,                      pcols, lchnk )
    call outfld( 'VM_CLUBB',         vm_output,                      pcols, lchnk )
    call outfld( 'WM_ZT_CLUBB',      wm_zt_output,               pcols, lchnk )

    call outfld( 'RELVAR',           relvar_pbuf,                  pcols, lchnk )
    call outfld( 'RHO_CLUBB',        rho(:,1:pver),           pcols, lchnk )
    call outfld( 'CLOUDFRAC_CLUBB',  alst_pbuf,                    pcols, lchnk )
    call outfld( 'CONCLD',           concld_pbuf,                  pcols, lchnk )
    call outfld( 'DP_CLD',           deepcu_pbuf,                  pcols, lchnk )
    call outfld( 'ZMDLF',            dlf_liq_out,             pcols, lchnk )
    call outfld( 'ZMDLFI',           dlf_ice_out,             pcols, lchnk )
    call outfld( 'CLUBB_GRID_SIZE',  grid_dx,                 pcols, lchnk )
    call outfld( 'QSATFAC',          qsatfac_pbuf,                 pcols, lchnk)


    ! --------------------------------------------------------------- !
    ! Writing state variables after EDMF scheme for detailed analysis !
    ! --------------------------------------------------------------- !
    if (do_clubb_mf) then

      do k = top_lev, pverp
        do i = 1, ncol
          k_clubb = k + 1 - top_lev
          mf_dry_a_output(i,k)     = mf_dry_a(i,k_clubb)
          mf_moist_a_output(i,k)   = mf_moist_a(i,k_clubb)
          mf_dry_w_output(i,k)     = mf_dry_w(i,k_clubb)
          mf_moist_w_output(i,k)   = mf_moist_w(i,k_clubb)
          mf_dry_qt_output(i,k)    = mf_dry_qt(i,k_clubb)
          mf_moist_qt_output(i,k)  = mf_moist_qt(i,k_clubb)
          mf_dry_thl_output(i,k)   = mf_dry_thl(i,k_clubb)
          mf_moist_thl_output(i,k) = mf_moist_thl(i,k_clubb)
          mf_dry_u_output(i,k)     = mf_dry_u(i,k_clubb)
          mf_moist_u_output(i,k)   = mf_moist_u(i,k_clubb)
          mf_dry_v_output(i,k)     = mf_dry_v(i,k_clubb)
          mf_moist_v_output(i,k)   = mf_moist_v(i,k_clubb)
          mf_moist_qc_output(i,k)  = mf_moist_qc(i,k_clubb)
          s_ae_output(i,k)         = s_ae(i,k_clubb)
          s_aw_output(i,k)         = s_aw(i,k_clubb)
          s_awthl_output(i,k)      = s_awthl(i,k_clubb)
          s_awqt_output(i,k)       = s_awqt(i,k_clubb)
          s_awql_output(i,k)       = s_awql(i,k_clubb)
          s_awqi_output(i,k)       = s_awqi(i,k_clubb)
          s_awu_output(i,k)        = s_awu(i,k_clubb)
          s_awv_output(i,k)        = s_awv(i,k_clubb)
          mf_thlflx_output(i,k)    = mf_thlflx(i,k_clubb)*rho(i,k)*cpair
          mf_qtflx_output(i,k)     = mf_qtflx(i,k_clubb)*rho(i,k)*latvap
        end do
      end do

      do k = 1, top_lev-1
        do i = 1, ncol
          mf_dry_a_output(i,k)     = 0._r8
          mf_moist_a_output(i,k)   = 0._r8
          mf_dry_w_output(i,k)     = 0._r8
          mf_moist_w_output(i,k)   = 0._r8
          mf_dry_qt_output(i,k)    = 0._r8
          mf_moist_qt_output(i,k)  = 0._r8
          mf_dry_thl_output(i,k)   = 0._r8
          mf_moist_thl_output(i,k) = 0._r8
          mf_dry_u_output(i,k)     = 0._r8
          mf_moist_u_output(i,k)   = 0._r8
          mf_dry_v_output(i,k)     = 0._r8
          mf_moist_v_output(i,k)   = 0._r8
          mf_moist_qc_output(i,k)  = 0._r8
          s_ae_output(i,k)         = 0._r8
          s_aw_output(i,k)         = 0._r8
          s_awthl_output(i,k)      = 0._r8
          s_awqt_output(i,k)       = 0._r8
          s_awql_output(i,k)       = 0._r8
          s_awqi_output(i,k)       = 0._r8
          s_awu_output(i,k)        = 0._r8
          s_awv_output(i,k)        = 0._r8
          mf_thlflx_output(i,k)    = 0._r8
          mf_qtflx_output(i,k)     = 0._r8
        end do
      end do

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
      call outfld( 'edmf_S_AE'     , s_ae_output,               pcols, lchnk )
      call outfld( 'edmf_S_AW'     , s_aw_output,               pcols, lchnk )
      call outfld( 'edmf_S_AWTHL'  , s_awthl_output,            pcols, lchnk )
      call outfld( 'edmf_S_AWQT'   , s_awqt_output,             pcols, lchnk )
      call outfld( 'edmf_S_AWU'    , s_awu_output,              pcols, lchnk )
      call outfld( 'edmf_S_AWV'    , s_awv_output,              pcols, lchnk )
      call outfld( 'edmf_thlflx'   , mf_thlflx_output,          pcols, lchnk )
      call outfld( 'edmf_qtflx'    , mf_qtflx_output,           pcols, lchnk )

    end if

    !  Output CLUBB history here
    if (stats_metadata%l_stats) then

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

      if (stats_metadata%l_output_rad_files) then
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
    call t_stopf('clubb_tend_cam:NAR')
#endif

    ! Cleanup err_info
    call cleanup_err_info_api(err_info)

    call t_stopf('clubb_tend_cam')

    return

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
subroutine ice_macro_tend(naai_pbuf,t,p,qv,qi,ni,xxls,deltat,stend,qvtend,qitend,nitend,vlen)

  use wv_sat_methods, only: wv_sat_qsat_ice

  integer,                   intent(in)  :: vlen
  real(r8), dimension(vlen), intent(in)  :: naai_pbuf   !Activated number of ice nuclei
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
     if (naai_pbuf(i) > 1.e-18_r8 .and. qv(i) > QSI(i)) then

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

    !----------------------- Input Variables -----------------------

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    real(kind=time_precision), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::   delt         ! Timestep (dtmain in CLUBB)         [s]

    !----------------------- Output Variables -----------------------
    type (stats), intent(out), dimension(pcols) :: &
      stats_zt,      & ! stats_zt grid
      stats_zm,      & ! stats_zm grid
      stats_rad_zt,  & ! stats_rad_zt grid
      stats_rad_zm,  & ! stats_rad_zm grid
      stats_sfc        ! stats_sfc


    !----------------------- Local Variables -----------------------

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

    logical :: l_error

    character(len=200) :: temp1, sub

    integer :: i, ntot, read_status, j
    integer :: iunit, ierr

    !----------------------- Begin Code -----------------------

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    stats_metadata%l_stats = l_stats_in

    stats_metadata%stats_tsamp = stats_tsamp_in
    stats_metadata%stats_tout  = stats_tout_in

    if ( .not. stats_metadata%l_stats ) then
       stats_metadata%l_stats_samp  = .false.
       stats_metadata%l_stats_last  = .false.
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
    stats_metadata%l_netcdf = .false.
    stats_metadata%l_grads  = .false.

    !  Check sampling and output frequencies
    do j = 1, pcols

      !  The model time step length, delt (which is dtmain), should multiply
      !  evenly into the statistical sampling time step length, stats_tsamp.
      if ( abs( stats_metadata%stats_tsamp/delt - floor(stats_metadata%stats_tsamp/delt) ) > 1.e-8_r8 ) then
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                          'the clubb time step (delt below)'
         write(fstderr,*) 'stats_tsamp = ', stats_metadata%stats_tsamp
         write(fstderr,*) 'delt = ', delt
         call endrun ("stats_init_clubb:  CLUBB stats_tsamp must be an even multiple of the timestep")
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
         l_error = .true.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_zt than allowed for by nvarmax_zt."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                          "in the stats namelist, or change nvarmax_zt."
         write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
         call endrun ("stats_init_clubb:  number of zt statistical variables exceeds limit")
      endif

      stats_zt(j)%num_output_fields = ntot
      stats_zt(j)%kk = nnzp - 1

      allocate( stats_zt(j)%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%z")

      allocate( stats_zt(j)%accum_field_values( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%accum_field_values")
      allocate( stats_zt(j)%accum_num_samples( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%accum_num_samples")
      allocate( stats_zt(j)%l_in_update( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%l_in_update")
      call stats_zero( stats_zt(j)%kk, stats_zt(j)%num_output_fields, stats_zt(j)%accum_field_values, &
                       stats_zt(j)%accum_num_samples, stats_zt(j)%l_in_update )

      allocate( stats_zt(j)%file%grid_avg_var( stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%file%grid_avg_var")
      allocate( stats_zt(j)%file%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) call endrun("stats_init_clubb: Failed to allocate stats_zt%file%z")

      !  Default initialization for array indices for zt
      call stats_init_zt_api( hydromet_dim, sclr_dim, edsclr_dim, &
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, &
                              clubb_vars_zt, &
                              l_error, &
                              stats_metadata, stats_zt(j) )

      !  Initialize zm (momentum points)

      i = 1
      do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  .and. &
                 len_trim(clubb_vars_zm(i)) /= 0    .and. &
                 i <= nvarmax_zm )
         i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_zm ) then
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_zm than allowed for by nvarmax_zm."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                          "in the stats namelist, or change nvarmax_zm."
         write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
         call endrun ("stats_init_clubb:  number of zm statistical variables exceeds limit")
      endif

      stats_zm(j)%num_output_fields = ntot
      stats_zm(j)%kk = nnzp

      allocate( stats_zm(j)%z( stats_zm(j)%kk ) )

      allocate( stats_zm(j)%accum_field_values( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%accum_num_samples( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%l_in_update( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      call stats_zero( stats_zm(j)%kk, stats_zm(j)%num_output_fields, stats_zm(j)%accum_field_values, &
                       stats_zm(j)%accum_num_samples, stats_zm(j)%l_in_update )

      allocate( stats_zm(j)%file%grid_avg_var( stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%file%z( stats_zm(j)%kk ) )

      call stats_init_zm_api( hydromet_dim, sclr_dim, edsclr_dim, &
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, &
                              clubb_vars_zm, &
                              l_error, &
                              stats_metadata, stats_zm(j) )

      !  Initialize rad_zt (radiation points)

      if (stats_metadata%l_output_rad_files) then

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

        stats_rad_zt(j)%num_output_fields = ntot
        stats_rad_zt(j)%kk = nnrad_zt

        allocate( stats_rad_zt(j)%z( stats_rad_zt(j)%kk ) )

        allocate( stats_rad_zt(j)%accum_field_values( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%accum_num_samples( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%l_in_update( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )

        call stats_zero( stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields, stats_rad_zt(j)%accum_field_values, &
                       stats_rad_zt(j)%accum_num_samples, stats_rad_zt(j)%l_in_update )

        allocate( stats_rad_zt(j)%file%grid_avg_var( stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%file%z( stats_rad_zt(j)%kk ) )

         call stats_init_rad_zt_api( clubb_vars_rad_zt, &
                                     l_error, &
                                     stats_metadata, stats_rad_zt(j) )

         !  Initialize rad_zm (radiation points)

         i = 1
         do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0 .and. &
                    len_trim(clubb_vars_rad_zm(i))   /= 0 .and. &
                    i <= nvarmax_rad_zm )
            i = i + 1
         end do
         ntot = i - 1
         if ( ntot == nvarmax_rad_zm ) then
            l_error = .true.  ! This will cause the run to stop.
            write(fstderr,*) "There are more statistical variables listed in ",  &
                             "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
            write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                             "in the stats namelist, or change nvarmax_rad_zm."
            write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
            call endrun ("stats_init_clubb:  number of rad_zm statistical variables exceeds limit")
         endif

         stats_rad_zm(j)%num_output_fields = ntot
         stats_rad_zm(j)%kk = nnrad_zm

         allocate( stats_rad_zm(j)%z( stats_rad_zm(j)%kk ) )

         allocate( stats_rad_zm(j)%accum_field_values( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%accum_num_samples( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%l_in_update( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )

         call stats_zero( stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields, stats_rad_zm(j)%accum_field_values, &
                       stats_rad_zm(j)%accum_num_samples, stats_rad_zm(j)%l_in_update )

         allocate( stats_rad_zm(j)%file%grid_avg_var( stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%file%z( stats_rad_zm(j)%kk ) )

         call stats_init_rad_zm_api( clubb_vars_rad_zm, &
                                     l_error, &
                                     stats_metadata, stats_rad_zm(j) )
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
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_sfc than allowed for by nvarmax_sfc."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                          "in the stats namelist, or change nvarmax_sfc."
         write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
         call endrun ("stats_init_clubb:  number of sfc statistical variables exceeds limit")
      endif

      stats_sfc(j)%num_output_fields = ntot
      stats_sfc(j)%kk = 1

      allocate( stats_sfc(j)%z( stats_sfc(j)%kk ) )

      allocate( stats_sfc(j)%accum_field_values( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%accum_num_samples( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%l_in_update( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )

      call stats_zero( stats_sfc(j)%kk, stats_sfc(j)%num_output_fields, stats_sfc(j)%accum_field_values, &
                       stats_sfc(j)%accum_num_samples, stats_sfc(j)%l_in_update )

      allocate( stats_sfc(j)%file%grid_avg_var( stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%file%z( stats_sfc(j)%kk ) )

      call stats_init_sfc_api( clubb_vars_sfc, &
                               l_error, &
                               stats_metadata, stats_sfc(j) )
    end do

    ! Check for errors

    if ( l_error ) then
       call endrun ('stats_init:  errors found')
    endif

    ! Now call add fields

    do i = 1, stats_zt(1)%num_output_fields

      temp1 = trim(stats_zt(1)%file%grid_avg_var(i)%name)
      sub   = temp1
      if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)

        call addfld( trim(sub), (/ 'ilev' /), 'A', &
                     trim(stats_zt(1)%file%grid_avg_var(i)%units), &
                     trim(stats_zt(1)%file%grid_avg_var(i)%description), &
                     sampled_on_subcycle=.true. )
    enddo

    do i = 1, stats_zm(1)%num_output_fields

      temp1 = trim(stats_zm(1)%file%grid_avg_var(i)%name)
      sub   = temp1
      if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)

       call addfld( trim(sub), (/ 'ilev' /), 'A', &
                    trim(stats_zm(1)%file%grid_avg_var(i)%units), &
                    trim(stats_zm(1)%file%grid_avg_var(i)%description), &
                    sampled_on_subcycle=.true. )
    enddo

    if (stats_metadata%l_output_rad_files) then

       do i = 1, stats_rad_zt(1)%num_output_fields
          temp1 = trim(stats_rad_zt(1)%file%grid_avg_var(i)%name)
          sub   = temp1
          if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
          call addfld( trim(sub), (/ 'ilev' /), 'A', &
                       trim(stats_rad_zt(1)%file%grid_avg_var(i)%units), &
                       trim(stats_rad_zt(1)%file%grid_avg_var(i)%description), &
                       sampled_on_subcycle=.true. )
       enddo

       do i = 1, stats_rad_zm(1)%num_output_fields
          temp1 = trim(stats_rad_zm(1)%file%grid_avg_var(i)%name)
          sub   = temp1
          if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
          call addfld( trim(sub), (/ 'ilev' /), 'A', &
                       trim(stats_rad_zm(1)%file%grid_avg_var(i)%units), &
                       trim(stats_rad_zm(1)%file%grid_avg_var(i)%description), &
                       sampled_on_subcycle=.true. )
       enddo
    endif

    do i = 1, stats_sfc(1)%num_output_fields
       temp1 = trim(stats_sfc(1)%file%grid_avg_var(i)%name)
       sub   = temp1
       if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
       call addfld( trim(sub), horiz_only, 'A', &
                    trim(stats_sfc(1)%file%grid_avg_var(i)%units), &
                    trim(stats_sfc(1)%file%grid_avg_var(i)%description), &
                    sampled_on_subcycle=.true. )
    enddo


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
    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pver,stats_zt%num_output_fields)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,stats_zt%num_output_fields)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pver,stats_rad_zt%num_output_fields)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%num_output_fields)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%num_output_fields)

    ! Local Variables

    integer :: i, k
    logical :: l_error

    !  Check if it is time to write to file

    if ( .not. stats_metadata%l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Compute averages
    call stats_avg( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, stats_zt%accum_num_samples )
    call stats_avg( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, stats_zm%accum_num_samples )
    if (stats_metadata%l_output_rad_files) then
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

        ! The data stored in stats types are ascending if l_ascending_grid = .true.
        if ( l_ascending_grid ) then
          out_zt(thecol,pver+1-k,i) = stats_zt%accum_field_values(1,1,k,i)
        else
          out_zt(thecol,top_lev-1+k,i) = stats_zt%accum_field_values(1,1,k,i)
        end if

        if(is_nan(out_zt(thecol,k,i))) out_zt(thecol,k,i) = 0.0_r8

      enddo
    enddo

    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zm%kk

        ! The data stored in stats types are ascending if l_ascending_grid = .true.
        if ( l_ascending_grid ) then
          out_zm(thecol,pverp+1-k,i) = stats_zm%accum_field_values(1,1,k,i)
        else
          out_zm(thecol,top_lev-1+k,i) = stats_zm%accum_field_values(1,1,k,i)
        end if

        if(is_nan(out_zm(thecol,k,i))) out_zm(thecol,k,i) = 0.0_r8

      enddo
    enddo

    if (stats_metadata%l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk

          ! The data stored in stats types are ascending if l_ascending_grid = .true.
          if ( l_ascending_grid ) then
            out_radzt(thecol,pver+1-k,i) = stats_rad_zt%accum_field_values(1,1,k,i)
          else
            out_radzt(thecol,top_lev-1+k,i) = stats_rad_zt%accum_field_values(1,1,k,i)
          end if

          if(is_nan(out_radzt(thecol,k,i))) out_radzt(thecol,k,i) = 0.0_r8

        enddo
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk

          ! The data stored in stats types are ascending if l_ascending_grid = .true.
          if ( l_ascending_grid ) then
            out_radzm(thecol,pverp+1-k,i) = stats_rad_zm%accum_field_values(1,1,k,i)
          else
            out_radzm(thecol,top_lev-1+k,i) = stats_rad_zm%accum_field_values(1,1,k,i)
          end if

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
    if (stats_metadata%l_output_rad_files) then
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
