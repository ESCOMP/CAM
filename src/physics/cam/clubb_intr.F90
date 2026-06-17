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

  use shr_kind_mod,        only: r8=>shr_kind_r8
  use ppgrid,              only: pver, pverp, pcols, begchunk, endchunk
  use phys_control,        only: phys_getopts
  use physconst,           only: cpair, gravit, rga, latvap, latice, zvir, rh2o, karman, pi, rair, omega
  use air_composition,     only: rairv, cpairv
  use cam_history_support, only: max_fieldname_len

  use spmd_utils,          only: masterproc
  use constituents,        only: pcnst, cnst_add, cnst_ndropmixed
  use ref_pres,            only: top_lev => trop_cloud_top_lev
  use scamMOD,             only: single_column, scm_clubb_iop_name, scm_cambfb_mode

#ifdef CLUBB_SGS
  use clubb,               only: clubb_init, clubb1_run, clubb2_run, clubb3_run, stats_zero
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

  type(clubb_config_flags_type), public  :: &
    clubb_config_flags

  real(r8), dimension(1,nparams), public :: &
    clubb_params_single_col    ! Adjustable CLUBB parameters (C1, C2 ...)

  ! Variables that contains all the statistics
  type (stats), public :: &
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
    ztodtptr ! model timestep
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
      p0_clubb                = 100000._r8,           &
      inv_p0_clubb            = 1._r8 / 100000._r8

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

  logical :: &
    clubb_l_intr_sfc_flux_smooth = .false. ! Add a locally calculated roughness to upwp and vpwp sfc fluxes

  logical :: &
    clubb_l_ascending_grid = .false.  ! Run clubb in ascending mode, which is opposite of the 
                                      ! cam grid the rest of this code uses, thus it requires
                                      ! an expensive array flipping step before calling advance_clubb_core.
                                      ! This is mainly for testing, it should not significantly change answers
  
  logical            :: lq(pcnst)
  logical            :: do_rainturb
  logical            :: clubb_do_adv
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
    clubb_iiPDF_type,                               & ! Selected option for the two-component normal
                                                      ! (double Gaussian) PDF type to use for the w, rt,
                                                      ! and theta-l (or w, chi, and eta) portion of
                                                      ! CLUBB's multivariate, two-component PDF.
    clubb_ipdf_call_placement = unset_i,            & ! Selected option for the placement of the call to
                                                      ! CLUBB's PDF.
    clubb_penta_solve_method = unset_i,             & ! Specifier for method to solve the penta-diagonal system
    clubb_tridiag_solve_method = unset_i,           & ! Specifier for method to solve tri-diagonal systems
    clubb_saturation_equation = unset_i,            & ! Specifier for which saturation formula to use
    clubb_grid_remap_method = unset_i,              & ! Specifier for which method should be used to
                                                      ! map values from one grid to another
                                                      ! (starts at 1, so 0 is an invalid option for this flag)
    clubb_grid_adapt_in_time_method = unset_i,      & ! Specifier for how the grid density method should
                                                      ! be constructed if the grid should be adapted over time
                                                      ! (set to 0 for no adaptation)
    clubb_fill_holes_type = unset_i                   ! Option for which type of hole filler to use in the 
                                                      ! fill_holes_vertical procedure


  logical :: &
    clubb_l_use_precip_frac,                        & ! Flag to use precipitation fraction in KK microphysics. The
                                                      ! precipitation fraction is automatically set to 1 when this
                                                      ! flag is turned off.
    clubb_l_predict_upwp_vpwp,                      & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                                      ! <w'thl'>, <sclr>, and <w'sclr'> in subroutine
                                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                                      ! approximated by eddy diffusivity when <u> and <v> are
                                                      ! advanced in subroutine advance_windm_edsclrm.
    clubb_l_ho_nontrad_coriolis,                    & ! Flag to implement the nontraditional Coriolis terms in the
                                                      ! prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
    clubb_l_ho_trad_coriolis,                       & ! Flag to implement the traditional Coriolis terms in the
                                                      ! prognostic equations of <v'w'> and <u'w'>.
    clubb_l_min_wp2_from_corr_wx,                   & ! Flag to base the threshold minimum value of wp2 on keeping
                                                      ! the overall correlation of w and x (w and rt, as well as w
                                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                                      ! to max_mag_correlation_flux.
    clubb_l_min_xp2_from_corr_wx,                   & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                                      ! thlp2) on keeping the overall correlation of w and x within
                                                      ! the limits of -max_mag_correlation_flux to
                                                      ! max_mag_correlation_flux.
    clubb_l_C2_cloud_frac,                          & ! Flag to use cloud fraction to adjust the value of the
                                                      ! turbulent dissipation coefficient, C2.
    clubb_l_diffuse_rtm_and_thlm,                   & ! Diffuses rtm and thlm
    clubb_l_stability_correct_Kh_N2_zm,             & ! Divides Kh_N2_zm by a stability factor
    clubb_l_calc_thlp2_rad,                         & ! Include the contribution of radiation to thlp2
    clubb_l_upwind_xpyp_ta,                         & ! This flag determines whether we want to use an upwind
                                                      ! differencing approximation rather than a centered
                                                      ! differencing for turbulent or mean advection terms. It
                                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                                      ! sclrpthlp.
    clubb_l_upwind_xm_ma,                           & ! This flag determines whether we want to use an upwind
                                                      ! differencing approximation rather than a centered
                                                      ! differencing for turbulent or mean advection terms. It
                                                      ! affects rtm, thlm, sclrm, um and vm.
    clubb_l_uv_nudge,                               & ! For wind speed nudging.
    clubb_l_rtm_nudge,                              & ! For rtm nudging
    clubb_l_tke_aniso,                              & ! For anisotropic turbulent kinetic energy, i.e.
                                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
    clubb_l_vert_avg_closure,                       & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                                      ! compute the varibles that are output from high order
                                                      ! closure
    clubb_l_trapezoidal_rule_zt,                    & ! If true, the trapezoidal rule is called for the
                                                      ! thermodynamic-level variables output from pdf_closure.
    clubb_l_trapezoidal_rule_zm,                    & ! If true, the trapezoidal rule is called for three
                                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                                      ! output from pdf_closure.
    clubb_l_call_pdf_closure_twice,                 & ! This logical flag determines whether or not to call
                                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                                      ! called first on thermodynamic levels and then on momentum
                                                      ! levels so that each variable is computed on its native
                                                      ! level.  If false, pdf_closure is only called on
                                                      ! thermodynamic levels, and variables which belong on
                                                      ! momentum levels are interpolated.
    clubb_l_standard_term_ta,                       & ! Use the standard discretization for the turbulent advection
                                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                                      ! pulled outside of the derivative in
                                                      ! advance_wp2_wp3_module.F90 and in
                                                      ! advance_xp2_xpyp_module.F90.
    clubb_l_partial_upwind_wp3,                     & ! Flag to use an "upwind" discretization rather
                                                      ! than a centered discretization for the portion
                                                      ! of the wp3 turbulent advection term for ADG1
                                                      ! that is linearized in terms of wp3<t+1>.
                                                      ! (Requires ADG1 PDF and clubb_l_standard_term_ta).
    clubb_l_godunov_upwind_wpxp_ta,                 & ! This flag determines whether we want to use an upwind
                                                      ! differencing approximation rather than a centered
                                                      ! differencing for turbulent advection terms.
                                                      ! It affects  wpxp only.
    clubb_l_godunov_upwind_xpyp_ta,                 & ! This flag determines whether we want to use an upwind
                                                      ! differencing approximation rather than a centered
                                                      ! differencing for turbulent advection terms. It affects
                                                      ! xpyp only.
    clubb_l_use_cloud_cover,                        & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                                      ! and rcm to help increase cloudiness at coarser grid
                                                      ! resolutions.
    clubb_l_diagnose_correlations,                  & ! Diagnose correlations instead of using fixed ones
    clubb_l_calc_w_corr,                            & ! Calculate the correlations between w and the hydrometeors
    clubb_l_const_Nc_in_cloud,                      & ! Use a constant cloud droplet conc. within cloud (K&K)
    clubb_l_fix_w_chi_eta_correlations,             & ! Use a fixed correlation for s and t Mellor(chi/eta)
    clubb_l_stability_correct_tau_zm,               & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                                      ! correction
    clubb_l_damp_wp2_using_em,                      & ! In wp2 equation, use a dissipation formula of
                                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
    clubb_l_do_expldiff_rtm_thlm,                   & ! Diffuse rtm and thlm explicitly
    clubb_l_Lscale_plume_centered,                  & ! Alternate that uses the PDF to compute the perturbed values
    clubb_l_diag_Lscale_from_tau,                   & ! First diagnose dissipation time tau, and then diagnose the
                                                      ! mixing length scale as Lscale = tau * tke
    clubb_l_use_C7_Richardson,                      & ! Parameterize C7 based on Richardson number
    clubb_l_use_C11_Richardson,                     & ! Parameterize C11 and C16 based on Richardson number
    clubb_l_use_shear_Richardson,                   & ! Use shear in the calculation of Richardson number
    clubb_l_brunt_vaisala_freq_moist,               & ! Use a different formula for the Brunt-Vaisala frequency in
                                                      ! saturated atmospheres (from Durran and Klemp, 1982)
    clubb_l_use_thvm_in_bv_freq,                    & ! Use thvm in the calculation of Brunt-Vaisala frequency
    clubb_l_rcm_supersat_adj,                       & ! Add excess supersaturated vapor to cloud water
    clubb_l_lmm_stepping,                           & ! Apply Linear Multistep Method (LMM) Stepping
    clubb_l_e3sm_config,                            & ! Run model with E3SM settings
    clubb_l_vary_convect_depth,                     & ! Flag used to calculate convective velocity using
                                                      ! a variable estimate of layer depth based on the depth
                                                      ! over which wpthlp is positive near the ground when true
                                                      ! More information can be found by
                                                      ! Looking at issue #905 on the clubb repo
    clubb_l_use_tke_in_wp3_pr_turb_term,            & ! Use TKE formulation for wp3 pr_turb term
    clubb_l_use_tke_in_wp2_wp3_K_dfsn,              & ! Use TKE in eddy diffusion for wp2 and wp3
    clubb_l_use_wp3_lim_with_smth_Heaviside,        & ! Flag to activate mods on wp3 limiters for conv test
    clubb_l_smooth_Heaviside_tau_wpxp,              & ! Use smooth Heaviside 'Peskin' in computation of invrs_tau
    clubb_l_modify_limiters_for_cnvg_test,          & ! Flag to activate mods on limiters for conv test
    clubb_l_enable_relaxed_clipping,                & ! Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats
    clubb_l_linearize_pbl_winds,                    & ! Flag to turn on code to linearize PBL winds
    clubb_l_single_C2_Skw,                          & ! Use a single Skewness dependent C2 for rtp2, thlp2, and
                                                      ! rtpthlp
    clubb_l_damp_wp3_Skw_squared,                   & ! Set damping on wp3 to use Skw^2 rather than Skw^4
    clubb_l_prescribed_avg_deltaz,                  & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
    clubb_l_update_pressure,                        & ! Flag for having CLUBB update pressure and exner
    clubb_l_mono_flux_lim_thlm,                     & ! Flag to turn on monotonic flux limiter for thlm
    clubb_l_mono_flux_lim_rtm,                      & ! Flag to turn on monotonic flux limiter for rtm
    clubb_l_mono_flux_lim_um,                       & ! Flag to turn on monotonic flux limiter for um
    clubb_l_mono_flux_lim_vm,                       & ! Flag to turn on monotonic flux limiter for vm
    clubb_l_mono_flux_lim_spikefix,                 & ! Flag to implement monotonic flux limiter code that
                                                      ! eliminates spurious drying tendencies at model top
    clubb_l_host_applies_sfc_fluxes,                & ! Whether the host model applies the surface fluxes
    clubb_l_wp2_fill_holes_tke,                     & ! Whether TKE is taken from up2 and vp2 to fill holes in wp2
    clubb_l_add_dycore_grid                           ! Flag to remap values from dycore grid

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
    wpthvp_idx, &       ! buoyancy flux
    wp2thvp_idx, &      ! second order buoyancy term
    wp2up_idx, &        ! w'^2 u'
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
    cld_idx, &         	! Cloud fraction
    concld_idx, &       ! Convective cloud fraction
    ast_idx, &          ! Stratiform cloud fraction
    alst_idx, &         ! Liquid stratiform cloud fraction
    aist_idx, &         ! Ice stratiform cloud fraction
    qlst_idx, &         ! Physical in-cloud LWC
    qist_idx, &         ! Physical in-cloud IWC
    dp_frac_idx, &      ! deep convection cloud fraction
    sh_frac_idx, &      ! shallow convection cloud fraction
    kvh_idx, &		      ! CLUBB eddy diffusivity on thermo levels
    pblh_idx, &         ! PBL pbuf
    icwmrdp_idx, &	    ! In cloud mixing ratio for deep convection
    tke_idx, &          ! turbulent kinetic energy
    tpert_idx, &        ! temperature perturbation from PBL
    fice_idx, &         ! fice_idx index in physics buffer
    cmeliq_idx, &       ! cmeliq_idx index in physics buffer
    relvar_idx, &       ! relative cloud water variance
    npccn_idx, &        ! liquid ccn number concentration
    naai_idx, &         ! ice number concentration
    prer_evap_idx, &    ! rain evaporation rate
    qrl_idx, &          ! longwave cooling rate
    qsatfac_idx, &      ! subgrid cloud water saturation scaling factor
    ice_supersat_idx, & ! ice cloud fraction for SILHS
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

  ! added pbuf fields for clubb to have restart bfb when ipdf_call_placement=2
  integer :: &          
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
    call pbuf_add_field('TTEND_CLUBB',     'physpkg', dtype_r8, (/pcols,pver /), ttend_clubb_idx )
    call pbuf_add_field('UPWP_CLUBB_GW',   'physpkg', dtype_r8, (/pcols,pverp/), upwp_clubb_gw_idx )
    call pbuf_add_field('VPWP_CLUBB_GW',   'physpkg', dtype_r8, (/pcols,pverp/), vpwp_clubb_gw_idx )
    call pbuf_add_field('THLP2_CLUBB_GW',  'physpkg', dtype_r8, (/pcols,pverp/), thlp2_clubb_gw_idx )
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
    call pbuf_add_field('WP2UP',      'global',  dtype_r8, (/pcols,nzt_clubb/), wp2up_idx)
    call pbuf_add_field('WP2RTP',     'global',  dtype_r8, (/pcols,nzt_clubb/), wp2rtp_idx)
    call pbuf_add_field('WP2THLP',    'global',  dtype_r8, (/pcols,nzt_clubb/), wp2thlp_idx)
    call pbuf_add_field('WPUP2',      'global',  dtype_r8, (/pcols,nzt_clubb/), wpup2_idx)
    call pbuf_add_field('WPVP2',      'global',  dtype_r8, (/pcols,nzt_clubb/), wpvp2_idx)

    call pbuf_add_field('RTP3',       'global', dtype_r8, (/pcols,nzt_clubb/), rtp3_idx)
    call pbuf_add_field('THLP3',      'global', dtype_r8, (/pcols,nzt_clubb/), thlp3_idx)
    call pbuf_add_field('UP3',        'global', dtype_r8, (/pcols,nzt_clubb/), up3_idx)
    call pbuf_add_field('VP3',        'global', dtype_r8, (/pcols,nzt_clubb/), vp3_idx)
    call pbuf_add_field('WP3_nadv',   'global', dtype_r8, (/pcols,nzt_clubb/), wp3_idx)

    call pbuf_add_field('UP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), up2_idx)
    call pbuf_add_field('VP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), vp2_idx)
    call pbuf_add_field('WP2_nadv',   'global', dtype_r8, (/pcols,nzm_clubb/), wp2_idx)

    ! Only in clubb_intr.F90 or SILHS
    call pbuf_add_field('ISS_FRAC',   'global', dtype_r8, (/pcols,nzt_clubb/), ice_supersat_idx)

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
                                clubb_rnevap_effic, &
                                clubb_l_ascending_grid
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
         clubb_l_ho_nontrad_coriolis, &
         clubb_l_ho_trad_coriolis, &
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
    clubb_l_ascending_grid            = .false.   ! Initialize to false

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
                                             clubb_l_ho_nontrad_coriolis, & ! Out
                                             clubb_l_ho_trad_coriolis, & ! Out
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
    call mpi_bcast(clubb_cloudtop_cooling,       1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_cloudtop_cooling")
    call mpi_bcast(clubb_rainevap_turb,          1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_rainevap_turb")
    call mpi_bcast(clubb_do_adv,                 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_do_adv")
    call mpi_bcast(clubb_l_ascending_grid,       1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_ascending_grid")
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
    call mpi_bcast(clubb_l_ho_nontrad_coriolis,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_ho_nontrad_coriolis")
    call mpi_bcast(clubb_l_ho_trad_coriolis,         1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_l_ho_trad_coriolis")
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
    if ( clubb_history          ) stats_metadata%l_stats            = .true.
    if ( clubb_rad_history      ) stats_metadata%l_output_rad_files = .true.
    if ( clubb_cloudtop_cooling ) do_cldcool                        = .true.
    if ( clubb_rainevap_turb    ) do_rainturb                       = .true.

    ! Check that all namelists have been set
    if ( clubb_timestep                   == unset_r8 ) call endrun( sub//": FATAL: clubb_timestep is not set")
    if ( clubb_rnevap_effic               == unset_r8 ) call endrun( sub//": FATAL:clubb_rnevap_effic  is not set")

    if ( clubb_c1                         == unset_r8 ) call endrun( sub//": FATAL: clubb_c1 is not set")
    if ( clubb_c1b                        == unset_r8 ) call endrun( sub//": FATAL: clubb_c1b is not set")
    if ( clubb_C2rt                       == unset_r8 ) call endrun( sub//": FATAL: clubb_C2rt is not set")
    if ( clubb_C2thl                      == unset_r8 ) call endrun( sub//": FATAL: clubb_C2thl is not set")
    if ( clubb_C2rtthl                    == unset_r8 ) call endrun( sub//": FATAL: clubb_C2rtthl is not set")
    if ( clubb_C4                         == unset_r8 ) call endrun( sub//": FATAL: clubb_C4 is not set")
    if ( clubb_C_uu_shr                   == unset_r8 ) call endrun( sub//": FATAL: clubb_C_uu_shr is not set")
    if ( clubb_C_uu_buoy                  == unset_r8 ) call endrun( sub//": FATAL: clubb_C_uu_buoy is not set")
    if ( clubb_c6rt                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c6rt is not set")
    if ( clubb_c6rtb                      == unset_r8 ) call endrun( sub//": FATAL: clubb_c6rtb is not set")
    if ( clubb_c6rtc                      == unset_r8 ) call endrun( sub//": FATAL: clubb_c6rtc is not set")
    if ( clubb_c6thl                      == unset_r8 ) call endrun( sub//": FATAL: clubb_c6thl is not set")
    if ( clubb_c6thlb                     == unset_r8 ) call endrun( sub//": FATAL: clubb_c6thlb is not set")
    if ( clubb_c6thlc                     == unset_r8 ) call endrun( sub//": FATAL: clubb_c6thlc is not set")
    if ( clubb_wpxp_L_thresh              == unset_r8 ) call endrun( sub//": FATAL: clubb_wpxp_L_thresh is not set")
    if ( clubb_C8                         == unset_r8 ) call endrun( sub//": FATAL: clubb_C8 is not set")
    if ( clubb_C8b                        == unset_r8 ) call endrun( sub//": FATAL: clubb_C8b is not set")
    if ( clubb_C7                         == unset_r8 ) call endrun( sub//": FATAL: clubb_C7 is not set")
    if ( clubb_C7b                        == unset_r8 ) call endrun( sub//": FATAL: clubb_C7b is not set")
    if ( clubb_c11                        == unset_r8 ) call endrun( sub//": FATAL: clubb_c11 is not set")
    if ( clubb_c11b                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c11b is not set")
    if ( clubb_c14                        == unset_r8 ) call endrun( sub//": FATAL: clubb_c14 is not set")
    if ( clubb_C_wp3_pr_turb              == unset_r8 ) call endrun( sub//": FATAL: clubb_C_wp3_pr_turb is not set")
    if ( clubb_c_K1                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K1 is not set")
    if ( clubb_c_K2                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K2 is not set")
    if ( clubb_nu2                        == unset_r8 ) call endrun( sub//": FATAL: clubb_nu2 is not set")
    if ( clubb_c_K8                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K8 is not set")
    if ( clubb_c_K9                       == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K9 is not set")
    if ( clubb_nu9                        == unset_r8 ) call endrun( sub//": FATAL: clubb_nu9 is not set")
    if ( clubb_c_K10                      == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K10 is not set")
    if ( clubb_c_K10h                     == unset_r8 ) call endrun( sub//": FATAL: clubb_c_K10h is not set")
    if ( clubb_C_invrs_tau_bkgnd          == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_bkgnd is not set")
    if ( clubb_C_invrs_tau_sfc            == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_sfc is not set")
    if ( clubb_C_invrs_tau_shear          == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_shear is not set")
    if ( clubb_C_invrs_tau_N2             == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_N2 is not set")
    if ( clubb_C_invrs_tau_N2_wp2         == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_N2_wp2 is not set")
    if ( clubb_C_invrs_tau_N2_xp2         == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_N2_xp2 is not set")
    if ( clubb_C_invrs_tau_N2_wpxp        == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_N2_wpxp is not set")
    if ( clubb_C_invrs_tau_N2_clear_wp3   == unset_r8 ) call endrun( sub//": FATAL: clubb_C_invrs_tau_N2_clear_wp3 is not set")
    if ( clubb_gamma_coef                 == unset_r8 ) call endrun( sub//": FATAL: clubb_gamma_coef is not set")
    if ( clubb_gamma_coefb                == unset_r8 ) call endrun( sub//": FATAL: clubb_gamma_coefb is not set")
    if ( clubb_beta                       == unset_r8 ) call endrun( sub//": FATAL: clubb_beta is not set")
    if ( clubb_lambda0_stability_coef     == unset_r8 ) call endrun( sub//": FATAL: clubb_lambda0_stability_coef is not set")
    if ( clubb_lmin_coef                  == unset_r8 ) call endrun( sub//": FATAL: clubb_lmin_coef is not set")
    if ( clubb_mult_coef                  == unset_r8 ) call endrun( sub//": FATAL: clubb_mult_coef is not set")
    if ( clubb_Skw_denom_coef             == unset_r8 ) call endrun( sub//": FATAL: clubb_Skw_denom_coef is not set")
    if ( clubb_skw_max_mag                == unset_r8 ) call endrun( sub//": FATAL: clubb_skw_max_mag is not set")
    if ( clubb_up2_sfc_coef               == unset_r8 ) call endrun( sub//": FATAL: clubb_up2_sfc_coef is not set")
    if ( clubb_C_wp2_splat                == unset_r8 ) call endrun( sub//": FATAL: clubb_C_wp2_splat is not set")
    if ( clubb_bv_efold                   == unset_r8 ) call endrun( sub//": FATAL: clubb_bv_efold is not set")
    if ( clubb_wpxp_Ri_exp                == unset_r8 ) call endrun( sub//": FATAL: clubb_wpxp_Ri_exp is not set")
    if ( clubb_z_displace                 == unset_r8 ) call endrun( sub//": FATAL: clubb_z_displace is not set")
    if ( clubb_detliq_rad                 == unset_r8 ) call endrun( sub//": FATAL: clubb_detliq_rad not set")
    if ( clubb_detice_rad                 == unset_r8 ) call endrun( sub//": FATAL: clubb_detice_rad not set")
    if ( clubb_ipdf_call_placement        == unset_i  ) call endrun( sub//": FATAL: clubb_ipdf_call_placement not set")
    if ( clubb_penta_solve_method         == unset_i  ) call endrun( sub//": FATAL: clubb_penta_solve_method not set")
    if ( clubb_tridiag_solve_method       == unset_i  ) call endrun( sub//": FATAL: clubb_tridiag_solve_method not set")
    if ( clubb_saturation_equation        == unset_i  ) call endrun( sub//": FATAL: clubb_saturation_equation not set")
    if ( clubb_grid_remap_method          == unset_i  ) call endrun( sub//": FATAL: clubb_grid_remap_method not set")
    if ( clubb_grid_adapt_in_time_method  == unset_i  ) call endrun( sub//": FATAL: clubb_grid_adapt_in_time_method not set")
    if ( clubb_fill_holes_type            == unset_i  ) call endrun( sub//": FATAL: clubb_fill_holes_type not set")

    if ( clubb_detphase_lowtemp           == unset_r8 ) call endrun( sub//": FATAL: clubb_detphase_lowtemp not set")
    if ( clubb_detphase_lowtemp        >= meltpt_temp ) call endrun( sub//": ERROR: clubb_detphase_lowtemp must be less than 268.15 K")

    call initialize_clubb_config_flags_type_api( clubb_iiPDF_type, &                        ! In        
                                                 clubb_ipdf_call_placement, &               ! In                
                                                 clubb_penta_solve_method, &                ! In                
                                                 clubb_tridiag_solve_method, &              ! In                  
                                                 clubb_saturation_equation, &               ! In                
                                                 clubb_grid_remap_method, &                 ! In              
                                                 clubb_grid_adapt_in_time_method, &         ! In                      
                                                 clubb_fill_holes_type, &                   ! In            
                                                 clubb_l_use_precip_frac, &                 ! In              
                                                 clubb_l_predict_upwp_vpwp, &               ! In   
                                                 clubb_l_ho_nontrad_coriolis, &             ! In
                                                 clubb_l_ho_trad_coriolis, &                ! In             
                                                 clubb_l_min_wp2_from_corr_wx, &            ! In                    
                                                 clubb_l_min_xp2_from_corr_wx, &            ! In                    
                                                 clubb_l_C2_cloud_frac, &                   ! In            
                                                 clubb_l_diffuse_rtm_and_thlm, &            ! In                    
                                                 clubb_l_stability_correct_Kh_N2_zm, &      ! In                          
                                                 clubb_l_calc_thlp2_rad, &                  ! In              
                                                 clubb_l_upwind_xpyp_ta, &                  ! In              
                                                 clubb_l_upwind_xm_ma, &                    ! In            
                                                 clubb_l_uv_nudge, &                        ! In        
                                                 clubb_l_rtm_nudge, &                       ! In        
                                                 clubb_l_tke_aniso, &                       ! In        
                                                 clubb_l_vert_avg_closure, &                ! In                
                                                 clubb_l_trapezoidal_rule_zt, &             ! In                  
                                                 clubb_l_trapezoidal_rule_zm, &             ! In                  
                                                 clubb_l_call_pdf_closure_twice, &          ! In                      
                                                 clubb_l_standard_term_ta, &                ! In                
                                                 clubb_l_partial_upwind_wp3, &              ! In                  
                                                 clubb_l_godunov_upwind_wpxp_ta, &          ! In                      
                                                 clubb_l_godunov_upwind_xpyp_ta, &          ! In                      
                                                 clubb_l_use_cloud_cover, &                 ! In              
                                                 clubb_l_diagnose_correlations, &           ! In                    
                                                 clubb_l_calc_w_corr, &                     ! In          
                                                 clubb_l_const_Nc_in_cloud, &               ! In                
                                                 clubb_l_fix_w_chi_eta_correlations, &      ! In                          
                                                 clubb_l_stability_correct_tau_zm, &        ! In                        
                                                 clubb_l_damp_wp2_using_em, &               ! In                
                                                 clubb_l_do_expldiff_rtm_thlm, &            ! In                    
                                                 clubb_l_Lscale_plume_centered, &           ! In                    
                                                 clubb_l_diag_Lscale_from_tau, &            ! In                    
                                                 clubb_l_use_C7_Richardson, &               ! In                
                                                 clubb_l_use_C11_Richardson, &              ! In                  
                                                 clubb_l_use_shear_Richardson, &            ! In                    
                                                 clubb_l_brunt_vaisala_freq_moist, &        ! In                        
                                                 clubb_l_use_thvm_in_bv_freq, &             ! In                  
                                                 clubb_l_rcm_supersat_adj, &                ! In                
                                                 clubb_l_damp_wp3_Skw_squared, &            ! In                    
                                                 clubb_l_prescribed_avg_deltaz, &           ! In                    
                                                 clubb_l_lmm_stepping, &                    ! In            
                                                 clubb_l_e3sm_config, &                     ! In          
                                                 clubb_l_vary_convect_depth, &              ! In                  
                                                 clubb_l_use_tke_in_wp3_pr_turb_term, &     ! In                          
                                                 clubb_l_use_tke_in_wp2_wp3_K_dfsn, &       ! In                        
                                                 clubb_l_use_wp3_lim_with_smth_Heaviside, & ! In                              
                                                 clubb_l_smooth_Heaviside_tau_wpxp, &       ! In                        
                                                 clubb_l_modify_limiters_for_cnvg_test, &   ! In                            
                                                 clubb_l_enable_relaxed_clipping, &         ! In                      
                                                 clubb_l_linearize_pbl_winds, &             ! In                  
                                                 clubb_l_mono_flux_lim_thlm, &              ! In                  
                                                 clubb_l_mono_flux_lim_rtm, &               ! In                
                                                 clubb_l_mono_flux_lim_um, &                ! In                
                                                 clubb_l_mono_flux_lim_vm, &                ! In                
                                                 clubb_l_mono_flux_lim_spikefix, &          ! In                      
                                                 clubb_l_host_applies_sfc_fluxes, &         ! In                      
                                                 clubb_l_wp2_fill_holes_tke, &              ! In                  
                                                 clubb_l_add_dycore_grid, &                 ! In              
                                                 clubb_config_flags )                       ! Out

#endif
  end subroutine clubb_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf_ini)
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
    use cam_abortutils,         only: endrun

    use time_manager,           only: is_first_step
    use phys_control,           only: phys_getopts
    use cam_logfile,            only: iulog
    use spmd_utils,             only: mpicom, mstrid=>masterprocid, mpi_character
    use clubb_api_module,       only: nvarmax_zm, &
         nvarmax_zt, &
         nvarmax_rad_zt, &
         nvarmax_rad_zm, &
         nvarmax_sfc, &
         var_length

#endif

    use namelist_utils,         only: find_group_name
    use units,                  only: getunit, freeunit
    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc

    implicit none

    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf_ini(:,:)

#ifdef CLUBB_SGS

    ! The similar name to clubb_history is unfortunate...
    logical :: history_amwg, history_clubb

    logical, parameter :: l_input_fields = .false. ! Always false for CAM-CLUBB.

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

    integer :: iunit, read_status, ierr

    integer :: errflg
    character(len=200) :: errmsg

    call phys_getopts(history_amwg_out=history_amwg, &
                      history_clubb_out=history_clubb, &
                      do_hb_above_clubb_out=do_hb_above_clubb)

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !

    !  Define physics buffers indexes
    cld_idx             = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx          = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx             = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx            = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx            = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx            = pbuf_get_index('QLST')        ! Physical in-stratus LWC
    qist_idx            = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx         = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx         = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx         = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx          = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    prer_evap_idx       = pbuf_get_index('PRER_EVAP')
    qrl_idx             = pbuf_get_index('QRL')
    cmfmc_sh_idx        = pbuf_get_index('CMFMC_SH')
    naai_idx            = pbuf_get_index('NAAI')
    npccn_idx           = pbuf_get_index('NPCCN')

    ! ----------------------------------------------------------------- !
    ! Add output fields for the history files
    ! ----------------------------------------------------------------- !

    !  These are default CLUBB output.  Not the higher order history budgets
    call addfld ('RHO_CLUBB',        (/ 'lev' /),  'A', 'kg/m3',     'Air Density',                                    sampled_on_subcycle = .true. )
    call addfld ('UP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',     'Zonal Velocity Variance',                        sampled_on_subcycle = .true. )
    call addfld ('VP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',     'Meridional Velocity Variance',                   sampled_on_subcycle = .true. )
    call addfld ('WP2_CLUBB',        (/ 'ilev' /), 'A', 'm2/s2',     'Vertical Velocity Variance',                     sampled_on_subcycle = .true. )
    call addfld ('WP2_ZT_CLUBB',     (/ 'lev' /),  'A', 'm2/s2',     'Vert Vel Variance on zt grid',                   sampled_on_subcycle = .true. )
    call addfld ('UPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',     'Zonal Momentum Flux',                            sampled_on_subcycle = .true. )
    call addfld ('VPWP_CLUBB',       (/ 'ilev' /), 'A', 'm2/s2',     'Meridional Momentum Flux',                       sampled_on_subcycle = .true. )
    call addfld ('WP3_CLUBB',        (/ 'lev' /),  'A', 'm3/s3',     'Third Moment Vertical Velocity',                 sampled_on_subcycle = .true. )
    call addfld ('WPTHLP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',      'Heat Flux',                                      sampled_on_subcycle = .true. )
    call addfld ('WPRTP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',      'Moisture Flux',                                  sampled_on_subcycle = .true. )
    call addfld ('RTP2_CLUBB',       (/ 'ilev' /), 'A', 'kg^2/kg^2', 'Moisture Variance',                              sampled_on_subcycle = .true. )
    call addfld ('RTP2_ZT_CLUBB',    (/ 'lev' /),  'A', 'kg^2/kg^2', 'Moisture Variance on zt grid',                   sampled_on_subcycle = .true. )
    call addfld ('THLP2_CLUBB',      (/ 'ilev' /), 'A', 'K^2',       'Temperature Variance',                           sampled_on_subcycle = .true. )
    call addfld ('THLP2_ZT_CLUBB',   (/ 'lev' /),  'A', 'K^2',       'Temperature Variance on zt grid',                sampled_on_subcycle = .true. )
    call addfld ('RTPTHLP_CLUBB',    (/ 'ilev' /), 'A', 'K kg/kg',   'Temp. Moist. Covariance',                        sampled_on_subcycle = .true. )
    call addfld ('RCM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water Mixing Ratio',                       sampled_on_subcycle = .true. )
    call addfld ('RTM_CLUBB',        (/ 'lev' /),  'A', 'kg/kg',     'Total Water Mixing Ratio',                       sampled_on_subcycle = .true. )
    call addfld ('THLM_CLUBB',       (/ 'lev' /),  'A', 'K',         'Liquid Water Potential Temperature',             sampled_on_subcycle = .true. )
    call addfld ('WPRCP_CLUBB',      (/ 'ilev' /), 'A', 'W/m2',      'Liquid Water Flux',                              sampled_on_subcycle = .true. )
    call addfld ('CLOUDFRAC_CLUBB',  (/ 'lev' /),  'A', 'fraction',  'Cloud Fraction',                                 sampled_on_subcycle = .true. )
    call addfld ('RCMINLAYER_CLUBB', (/ 'lev' /),  'A', 'kg/kg',     'Cloud Water in Layer',                           sampled_on_subcycle = .true. )
    call addfld ('CLOUDCOVER_CLUBB', (/ 'lev' /),  'A', 'fraction',  'Cloud Cover',                                    sampled_on_subcycle = .true. )
    call addfld ('WPTHVP_CLUBB',     (/ 'ilev' /), 'A', 'W/m2',      'Buoyancy Flux',                                  sampled_on_subcycle = .true. )
    call addfld ('RVMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Water vapor tendency',                           sampled_on_subcycle = .true. )
    call addfld ('STEND_CLUBB',      (/ 'lev' /),  'A', 'J/(kg s)',  'Static energy tendency',                         sampled_on_subcycle = .true. )
    call addfld ('RCMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Liquid Water Tendency',                    sampled_on_subcycle = .true. )
    call addfld ('RIMTEND_CLUBB',    (/ 'lev' /),  'A', 'kg/kg /s',  'Cloud Ice Tendency',                             sampled_on_subcycle = .true. )
    call addfld ('UTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',    'U-wind Tendency',                                sampled_on_subcycle = .true. )
    call addfld ('VTEND_CLUBB',      (/ 'lev' /),  'A', 'm/s /s',    'V-wind Tendency',                                sampled_on_subcycle = .true. )
    call addfld ('ZT_CLUBB',         (/ 'lev' /),  'A', 'm',         'Thermodynamic Heights',                          sampled_on_subcycle = .true. )
    call addfld ('ZM_CLUBB',         (/ 'ilev' /), 'A', 'm',         'Momentum Heights',                               sampled_on_subcycle = .true. )
    call addfld ('UM_CLUBB',         (/ 'lev' /),  'A', 'm/s',       'Zonal Wind',                                     sampled_on_subcycle = .true. )
    call addfld ('VM_CLUBB',         (/ 'lev' /),  'A', 'm/s',       'Meridional Wind',                                sampled_on_subcycle = .true. )
    call addfld ('WM_ZT_CLUBB',      (/ 'lev' /),  'A', 'm/s',       'Vertical Velocity',                              sampled_on_subcycle = .true. )
    call addfld ('CLDST',            (/ 'lev' /),  'A', 'fraction',  'Stratus cloud fraction',                         sampled_on_subcycle = .true. )
    call addfld ('ZMDLF',            (/ 'lev' /),  'A', 'kg/kg/s',   'Detrained liquid water from ZM convection',      sampled_on_subcycle = .true. )
    call addfld ('TTENDICE',         (/ 'lev' /),  'A', 'K/s',       'T tendency from Ice Saturation Adjustment',      sampled_on_subcycle = .true. )
    call addfld ('QVTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',   'Q tendency from Ice Saturation Adjustment',      sampled_on_subcycle = .true. )
    call addfld ('QITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',   'CLDICE tendency from Ice Saturation Adjustment', sampled_on_subcycle = .true. )
    call addfld ('NITENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',   'NUMICE tendency from Ice Saturation Adjustment', sampled_on_subcycle = .true. )

    call addfld ('PBLH',                    horiz_only,   'A', 'm',         'PBL height',         sampled_on_subcycle=.true.)
    call addfld ('PDFP_RTP2_CLUBB',  (/ 'lev' /),  'A', 'kg^2/kg^2', 'PDF Rtot Variance',  sampled_on_subcycle=.true.)

    call addfld ('QCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'CLDICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('NCTENDICE',        (/ 'lev' /),  'A', 'kg/kg/s',  'NUMICE tendency from Ice Saturation Adjustment', sampled_on_subcycle=.true.)
    call addfld ('FQTENDICE',        (/ 'lev' /),  'A', 'fraction', 'Frequency of Ice Saturation Adjustment',         sampled_on_subcycle=.true.)

    call addfld ('DPDLFLIQ',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained liquid water from deep convection',    sampled_on_subcycle=.true.)
    call addfld ('DPDLFICE',         (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice from deep convection',             sampled_on_subcycle=.true.)
    call addfld ('DPDLFT',           (/ 'lev' /),  'A', 'K/s',      'T-tendency due to deep convective detrainment',  sampled_on_subcycle=.true.)
    call addfld ('RELVAR',           (/ 'lev' /),  'A', '-',        'Relative cloud water variance',                  sampled_on_subcycle=.true.)
    call addfld ('CLUBB_GRID_SIZE',  horiz_only,   'A', 'm',        'Horizontal grid box size seen by CLUBB',         sampled_on_subcycle=.true.)


    call addfld ('ZMDLFI',           (/ 'lev' /),  'A', 'kg/kg/s',  'Detrained ice water from ZM convection',     sampled_on_subcycle=.true.)
    call addfld ('CONCLD',           (/ 'lev' /),  'A', 'fraction', 'Convective cloud cover',                     sampled_on_subcycle=.true.)
    call addfld ('CMELIQ',           (/ 'lev' /),  'A', 'kg/kg/s',  'Rate of cond-evap of liq within the cloud',  sampled_on_subcycle=.true.)
    call addfld ('DETNLIQTND',       (/ 'lev' /),  'A', '1/kg/s',   'CLDNUM tendency in detrained water',         sampled_on_subcycle=.true.)

    call addfld ('KVH_CLUBB',        (/ 'ilev' /), 'A', 'm2/s', 'CLUBB vertical diffusivity of heat/moisture on interface levels', sampled_on_subcycle=.true.)
    call addfld ('QSATFAC',          (/ 'lev' /),  'A', '-',    'Subgrid cloud water saturation scaling factor',      sampled_on_subcycle=.true.)
    call addfld ('ELEAK_CLUBB',      horiz_only,   'A', 'W/m2', 'CLUBB energy leak',                                  sampled_on_subcycle=.true.)
    call addfld ('TFIX_CLUBB',       horiz_only,   'A', 'K',    'Temperature increment to conserve energy',           sampled_on_subcycle=.true.)

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

    ! ----------------------------------------------------------------- !
    ! Make all of this output default, this is not CLUBB history
    ! ----------------------------------------------------------------- !

    if (clubb_do_adv .or. history_clubb) then
       call add_default('RELVAR',                   1, ' ')
       call add_default('RHO_CLUBB',                1, ' ')
       call add_default('UP2_CLUBB',                1, ' ')
       call add_default('VP2_CLUBB',                1, ' ')
       call add_default('WP2_CLUBB',                1, ' ')
       call add_default('WP2_ZT_CLUBB',             1, ' ')
       call add_default('WP3_CLUBB',                1, ' ')
       call add_default('UPWP_CLUBB',               1, ' ')
       call add_default('VPWP_CLUBB',               1, ' ')
       call add_default('WPTHLP_CLUBB',             1, ' ')
       call add_default('WPRTP_CLUBB',              1, ' ')
       call add_default('RTP2_CLUBB',               1, ' ')
       call add_default('RTP2_ZT_CLUBB',            1, ' ')
       call add_default('PDFP_RTP2_CLUBB',          1, ' ')
       call add_default('THLP2_CLUBB',              1, ' ')
       call add_default('THLP2_ZT_CLUBB',           1, ' ')
       call add_default('RTPTHLP_CLUBB',            1, ' ')
       call add_default('RCM_CLUBB',                1, ' ')
       call add_default('RTM_CLUBB',                1, ' ')
       call add_default('THLM_CLUBB',               1, ' ')
       call add_default('WPRCP_CLUBB',              1, ' ')
       call add_default('CLOUDFRAC_CLUBB',          1, ' ')
       call add_default('RCMINLAYER_CLUBB',         1, ' ')
       call add_default('CLOUDCOVER_CLUBB',         1, ' ')
       call add_default('WPTHVP_CLUBB',             1, ' ')
       call add_default('RVMTEND_CLUBB',            1, ' ')
       call add_default('STEND_CLUBB',              1, ' ')
       call add_default('RCMTEND_CLUBB',            1, ' ')
       call add_default('RIMTEND_CLUBB',            1, ' ')
       call add_default('UTEND_CLUBB',              1, ' ')
       call add_default('VTEND_CLUBB',              1, ' ')
       call add_default('ZT_CLUBB',                 1, ' ')
       call add_default('ZM_CLUBB',                 1, ' ')
       call add_default('UM_CLUBB',                 1, ' ')
       call add_default('VM_CLUBB',                 1, ' ')
       call add_default('WM_ZT_CLUBB',              1, ' ')
       call add_default('PBLH',                     1, ' ')
       call add_default('CONCLD',                   1, ' ')
    endif

    if (history_amwg) then
       call add_default('PBLH',           1, ' ')
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

       call pbuf_set_field(pbuf_ini, wp2_idx,     w_tol_sqd)
       call pbuf_set_field(pbuf_ini, wp3_idx,     0.0_r8)
       call pbuf_set_field(pbuf_ini, wpthlp_idx,  0.0_r8)
       call pbuf_set_field(pbuf_ini, wprtp_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini, rtpthlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf_ini, rtp2_idx,    rt_tol**2)
       call pbuf_set_field(pbuf_ini, thlp2_idx,   thl_tol**2)
       call pbuf_set_field(pbuf_ini, up2_idx,     w_tol_sqd)
       call pbuf_set_field(pbuf_ini, vp2_idx,     w_tol_sqd)

       call pbuf_set_field(pbuf_ini, rtp3_idx,    0.0_r8)
       call pbuf_set_field(pbuf_ini, thlp3_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini, up3_idx,     0.0_r8)
       call pbuf_set_field(pbuf_ini, vp3_idx,     0.0_r8)

       call pbuf_set_field(pbuf_ini, upwp_idx,          0.0_r8)
       call pbuf_set_field(pbuf_ini, vpwp_idx,          0.0_r8)
       call pbuf_set_field(pbuf_ini, wpthvp_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, wp2thvp_idx,       0.0_r8)
       call pbuf_set_field(pbuf_ini, rtpthvp_idx,       0.0_r8)
       call pbuf_set_field(pbuf_ini, thlpthvp_idx,      0.0_r8)
       call pbuf_set_field(pbuf_ini, tke_idx,           0.0_r8)
       call pbuf_set_field(pbuf_ini, kvh_idx,           0.0_r8)
       call pbuf_set_field(pbuf_ini, wp2rtp_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, wp2thlp_idx,       0.0_r8)
       call pbuf_set_field(pbuf_ini, uprcp_idx,         0.0_r8)
       call pbuf_set_field(pbuf_ini, vprcp_idx,         0.0_r8)
       call pbuf_set_field(pbuf_ini, rc_coef_zm_idx,    0.0_r8)
       call pbuf_set_field(pbuf_ini, wp4_idx,           0.0_r8)
       call pbuf_set_field(pbuf_ini, wpup2_idx,         0.0_r8)
       call pbuf_set_field(pbuf_ini, wpvp2_idx,         0.0_r8)
       call pbuf_set_field(pbuf_ini, wp2up2_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, wp2vp2_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, ice_supersat_idx,  0.0_r8)

       ! Initialize SILHS covariance contributions
       call pbuf_set_field(pbuf_ini, rtp2_mc_zt_idx,    0.0_r8)
       call pbuf_set_field(pbuf_ini, thlp2_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini, wprtp_mc_zt_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini, wpthlp_mc_zt_idx,  0.0_r8)
       call pbuf_set_field(pbuf_ini, rtpthlp_mc_zt_idx, 0.0_r8)

       call pbuf_set_field(pbuf_ini, pdf_zm_w_1_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, pdf_zm_w_2_idx,        0.0_r8)
       call pbuf_set_field(pbuf_ini, pdf_zm_varnce_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf_ini, pdf_zm_varnce_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf_ini, pdf_zm_mixt_frac_idx,  0.0_r8)

       call pbuf_set_field(pbuf_ini,  ttend_clubb_idx,      0.0_r8)
       call pbuf_set_field(pbuf_ini,  upwp_clubb_gw_idx,    0.0_r8)
       call pbuf_set_field(pbuf_ini,  vpwp_clubb_gw_idx,    0.0_r8)
       call pbuf_set_field(pbuf_ini,  thlp2_clubb_gw_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini,  wpthlp_clubb_gw_idx,  0.0_r8)

       call pbuf_set_field(pbuf_ini,  ttend_clubb_mc_idx,     0.0_r8)
       call pbuf_set_field(pbuf_ini,  upwp_clubb_gw_mc_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini,  vpwp_clubb_gw_mc_idx,   0.0_r8)
       call pbuf_set_field(pbuf_ini,  thlp2_clubb_gw_mc_idx,  0.0_r8)
       call pbuf_set_field(pbuf_ini,  wpthlp_clubb_gw_mc_idx, 0.0_r8)

    endif

    ! The following is physpkg, so it needs to be initialized every time
    call pbuf_set_field(pbuf_ini, fice_idx,    0.0_r8)

    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

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
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_zz")
    call mpi_bcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_rad_zt")
    call mpi_bcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_rad_zm")
    call mpi_bcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpi_character, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(subr//": FATAL: mpi_bcast: clubb_vars_sfc")

    ! Call the CAM-SIMA layer
    call clubb_init(pcols, pver, pverp, pcnst, begchunk, endchunk, & ! in
                    masterproc, mpicom, mstrid, mpi_character, & ! in
                    iulog, max_fieldname_len, & ! in
                    sclr_dim, hydromet_dim, nzt_clubb, nzm_clubb, & ! in
                    l_implemented, l_input_fields, clubb_l_do_expldiff_rtm_thlm, & ! in
                    cnst_ndropmixed, subcol_scheme, & ! in
                    clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, & ! in
                    clubb_vars_rad_zt, clubb_vars_rad_zm, & ! in
                    edsclr_dim, clubb_params_single_col, & ! inout
                    out_zt, out_zm, out_sfc, out_radzt, out_radzm, & ! inout
                    clubb_C2rtthl, clubb_C8, clubb_c11, clubb_c11b, clubb_c14, & ! inout
                    clubb_C_wp3_pr_turb, clubb_c_K10, clubb_mult_coef, & ! inout
                    clubb_Skw_denom_coef, clubb_C2rt, clubb_C2thl, clubb_beta, & ! inout
                    clubb_c6rt, clubb_c6rtb, clubb_c6rtc, clubb_c6thl, clubb_c6thlb, & ! inout
                    clubb_c6thlc, clubb_wpxp_L_thresh, clubb_C7, clubb_C7b, & ! inout
                    clubb_gamma_coef, clubb_c_K10h, clubb_lambda0_stability_coef, & ! inout
                    clubb_lmin_coef, clubb_C8b, clubb_skw_max_mag, clubb_C1, clubb_C1b, & ! inout
                    clubb_gamma_coefb, clubb_up2_sfc_coef, clubb_C4, clubb_C_uu_shr, & ! inout
                    clubb_C_uu_buoy, clubb_c_K1, clubb_c_K2, clubb_nu2, clubb_c_K8, & ! inout
                    clubb_c_K9, clubb_nu9, clubb_C_wp2_splat, clubb_C_invrs_tau_bkgnd, & ! inout
                    clubb_C_invrs_tau_sfc, clubb_C_invrs_tau_shear, clubb_C_invrs_tau_N2, & ! inout
                    clubb_C_invrs_tau_N2_wp2, clubb_C_invrs_tau_N2_xp2, clubb_C_invrs_tau_N2_wpxp, & ! inout
                    clubb_C_invrs_tau_N2_clear_wp3, clubb_bv_efold, clubb_wpxp_Ri_exp, clubb_z_displace, & ! inout
                    lq, stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, & ! inout
                    pdf_params_chnk, pdf_params_zm_chnk, pdf_implicit_coefs_terms_chnk, & ! inout
                    stats_metadata, hm_metadata, clubb_config_flags, sclr_idx, & ! inout
                    errmsg, errflg ) ! out

    if (errflg /= 0) then
      call endrun(errmsg)
    end if

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
    use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

    use constituents,   only: cnst_get_ind, cnst_type
    use camsrfexch,     only: cam_in_t
    use time_manager,   only: is_first_step
    use cam_abortutils, only: endrun
    use cam_logfile,    only: iulog
!BAS still using tropopause but prob need to switch over to SIMA side
    use tropopause,     only: tropopause_findChemTrop
    use time_manager,   only: get_nstep, is_first_restart_step, get_curr_calday, get_calday
    use perf_mod,       only: t_startf, t_stopf

#ifdef CLUBB_SGS
    use spmd_utils, only: iam
    use clubb_api_module, only: &
      zm2zt_api, &
      init_pdf_params_api, &
      init_pdf_implicit_coefs_terms_api, &
      iiPDF_new, &
      iiPDF_new_hybrid

    ! Import setup for CLUBB error messaging
    use clubb_api_module, only: &
      err_info_type,        &
      cleanup_err_info_api

    use cldfrc2m,                  only: aist_vector, rhmini_const, rhmaxi_const, rhminis_const, rhmaxis_const
    use cam_history,               only: outfld

    use macrop_driver,             only: liquid_macro_tend

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

    real(r8), pointer, dimension(:,:) :: wp2_pbuf                   ! vertical velocity variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wp3_pbuf                   ! third moment of vertical velocity		[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: wpthlp_pbuf                ! turbulent flux of thetal			[m/s K]
    real(r8), pointer, dimension(:,:) :: wprtp_pbuf                 ! turbulent flux of moisture			[m/s kg/kg]
    real(r8), pointer, dimension(:,:) :: rtpthlp_pbuf               ! covariance of thetal and qt			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: rtp2_pbuf                  ! moisture variance				[kg^2/kg^2]
    real(r8), pointer, dimension(:,:) :: thlp2_pbuf                 ! temperature variance				[K^2]
    real(r8), pointer, dimension(:,:) :: rtp3_pbuf                  ! moisture 3rd order				[kg^3/kg^3]
    real(r8), pointer, dimension(:,:) :: thlp3_pbuf                 ! temperature 3rd order			[K^3]
    real(r8), pointer, dimension(:,:) :: up2_pbuf                   ! east-west wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vp2_pbuf                   ! north-south wind variance			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: up3_pbuf                   ! east-west wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: vp3_pbuf                   ! north-south wind 3rd order			[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: upwp_pbuf                  ! east-west momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: vpwp_pbuf                  ! north-south momentum flux			[m^2/s^2]
    real(r8), pointer, dimension(:,:) :: wpthvp_pbuf                ! w'th_v' (momentum levels)			[m/s K]
    real(r8), pointer, dimension(:,:) :: wp2thvp_pbuf               ! w'^2 th_v' (thermodynamic levels)		[m^2/s^2 K]
    real(r8), pointer, dimension(:,:) :: wp2up_pbuf                 ! w'^2 u' (thermodynamic levels)		[m^3/s^3]
    real(r8), pointer, dimension(:,:) :: rtpthvp_pbuf               ! r_t'th_v' (momentum levels)			[kg/kg K]
    real(r8), pointer, dimension(:,:) :: thlpthvp_pbuf              ! th_l'th_v' (momentum levels)			[K^2]
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_1_pbuf            ! work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_w_2_pbuf            ! work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_1_pbuf     ! work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_varnce_w_2_pbuf     ! work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: pdf_zm_mixt_frac_pbuf      ! work pointer for pdf_params_zm
    real(r8), pointer, dimension(:,:) :: wp2rtp_pbuf                ! w'^2 rt' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2thlp_pbuf               ! w'^2 thl' (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: uprcp_pbuf                 ! < u' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: vprcp_pbuf                 ! < v' r_c' > (momentum levels)
    real(r8), pointer, dimension(:,:) :: rc_coef_zm_pbuf            ! Coef. of X'r_c' in Eq. (34) (t-levs.)
    real(r8), pointer, dimension(:,:) :: wp4_pbuf                   ! w'^4 (momentum levels
    real(r8), pointer, dimension(:,:) :: wpup2_pbuf                 ! w'u'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wpvp2_pbuf                 ! w'v'^2 (thermodynamic levels)
    real(r8), pointer, dimension(:,:) :: wp2up2_pbuf                ! w'^2 u'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: wp2vp2_pbuf                ! w'^2 v'^2 (momentum levels)
    real(r8), pointer, dimension(:,:) :: cld_pbuf                   ! cloud fraction 				[fraction]
    real(r8), pointer, dimension(:,:) :: concld_pbuf                ! convective cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: ast_pbuf                   ! stratiform cloud fraction			[fraction]
    real(r8), pointer, dimension(:,:) :: alst_pbuf                  ! liquid stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: aist_pbuf                  ! ice stratiform cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: qlst_pbuf                  ! Physical in-stratus LWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: qist_pbuf                  ! Physical in-stratus IWC			[kg/kg]
    real(r8), pointer, dimension(:,:) :: deepcu_pbuf                ! deep convection cloud fraction		[fraction]
    real(r8), pointer, dimension(:,:) :: shalcu_pbuf                ! shallow convection cloud fraction 		[fraction]
    real(r8), pointer, dimension(:,:) :: khzm_pbuf                  ! CLUBB's eddy diffusivity of heat/moisture on momentum  levels [m^2/s]
    real(r8), pointer, dimension(:)   :: pblh_pbuf                  ! planetary boundary layer height                [m]
    real(r8), pointer, dimension(:,:) :: tke_pbuf                   ! turbulent kinetic energy                     [m^2/s^2]
    real(r8), pointer, dimension(:,:) :: dp_icwmr_pbuf              ! deep convection in cloud mixing ratio        [kg/kg]
    real(r8), pointer, dimension(:,:) :: ice_supersat_frac_pbuf     ! Cloud fraction of ice clouds (pver)[fraction]
    real(r8), pointer, dimension(:,:) :: relvar_pbuf                ! relative cloud water variance                [-]
    real(r8), pointer, dimension(:,:) :: naai_pbuf
    real(r8), pointer, dimension(:,:) :: cmeliq_pbuf
    real(r8), pointer, dimension(:,:) :: cmfmc_sh_pbuf              ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]

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

!BAS to sima    real(r8), parameter :: &
!BAS to sima      rad2deg=180.0_r8/pi

    character(len=*), parameter :: subr='clubb_tend_cam'

    type(physics_state) :: state_loc                ! Local copy of state variable
    type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header

    type(grid) :: &
      gr          ! CLUBB grid data structure

    real(r8), dimension(state%ncol) :: &
      grid_dx, grid_dy                    ! CAM grid [m]

    real(r8), dimension(state%ncol,nzt_clubb) :: &
      rtm,                            & ! mean moisture mixing ratio			              [kg/kg]
      thlm,                           & ! mean temperature				                      [K]
      rcm,                            & ! CLUBB cloud water mixing ratio                [kg/kg]
      um,                             & ! mean east-west wind				                    [m/s]
      vm,                             & ! mean north-south wind			                    [m/s]
      wm_zt,                          & ! w mean wind component on thermo. levels   	  [m/s]
      rho_zt,                         & ! Air density on thermo levels                  [kg/m^3]
      exner,                          & ! Exner function (thermodynamic levels)         [-]
      rtp2_zt,                        & ! CLUBB R-tot variance on thermo levs
      thl2_zt,                        & ! CLUBB Theta-l variance on thermo levs         [K^2]
      wp2_zt,                         & ! CLUBB W variance on theromo levs              [m^2/s^2]
      cloud_frac,                     & ! CLUBB output of cloud fraction                [fraction]
      rcm_in_layer,                   & ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
      zt_g                              ! Thermodynamic grid of CLUBB		      	        [m]

    real(r8), dimension(state%ncol,nzm_clubb) :: &
      rho_zm,                   & ! Air density on momentum levels                        [kg/m^3]
      wprcp,                    & ! CLUBB output of flux of liquid water                  [kg/kg m/s]
      zi_g,                     & ! Momentum grid of CLUBB		      	                    [m]

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
      mf_thlflx,  mf_qtflx
      
    ! Variables used for output (zm)
    real(r8), dimension(pcols,pverp) :: &
      zi_output,                & ! output for momentum CLUBB grid                [m]
      wpthlp_output,            & ! Heat flux output variable                     [W/m2]
      rtpthlp_output,           & ! rtpthlp ouptut                                [K kg/kg]
      wprtp_output,             & ! Total water flux output variable              [W/m2]
      wp2_output,               &
      up2_output,               &
      vp2_output,               &
      upwp_output,              &
      vpwp_output,              &
      rtp2_output,              &
      wprcp_clubb_output,       &
      wpthvp_clubb_output,      &
      thlp2_output,             &
      dlf_liq_out,              & ! Detrained liquid water from ZM                [kg/kg/s]
      dlf_ice_out,              & ! Detrained ice water from ZM                   [kg/kg/s]

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
      rvmtend_clubb_output,           &
      rcmtend_clubb_output,           &
      rimtend_clubb_output,           &
      stend_clubb_output,             &
      utend_clubb_output,             &
      vtend_clubb_output,             &
      dpdlfliq_output,                &
      dpdlfice_output,                &
      dpdlft_output,                  &
      detnliquid_output,              &
      zt_output,                      & ! output for the thermo CLUBB grid              [m]
      rtp2_zt_output,                 & ! CLUBB R-tot variance on thermo levs           [kg^2/kg^2]
      wp3_output,                     & ! wp3 output                                    [m^3/s^3]
      thl2_zt_output,                 & ! CLUBB Theta-l variance on thermo levs
      wp2_zt_output,                  & 
      rcm_in_layer_output,            & ! CLUBB in-cloud liquid water mixing ratio	    [kg/kg]
      pdfp_rtp2_output,               & ! Calculated R-tot variance from pdf_params     [kg^2/kg^2]
      wm_zt_output,                   & ! CLUBB mean W on thermo levs output            [m/s]
      rcm_output,                     &
      rtm_output,                     &
      thlm_output,                    &
      um_output,                      &
      vm_output,                      &
      rho_output

    real(r8), dimension(pcols) :: &
      rhmini,           &
      rhmaxi,           &
      se_dis,           &
      eleak

    real(r8), dimension(pcols,pver) :: &
      invrs_cpairv

    real(r8) :: &
      inv_exner_tmp,            & ! Inverse exner function consistent with CLUBB  [-]
      mean_rt,                  & ! Calculated R-tot mean from pdf_params (temp)  [kg/kg]
      apply_const

    intrinsic :: max

    logical, dimension(pcnst) :: &
      lqice

    character(len=200) :: temp1, sub             ! Strings needed for CLUBB output
    character(len=512) :: errmsg

    integer, dimension(pcols) :: &
      clubbtop_pbuf, &
      troplev

    integer :: &
      errflg, &
      j, k, t, ixind, n,      & ! Loop variables
      k_cam, k_clubb, & ! Loop variables
      ixcldice, ixcldliq, ixnumliq, &
      ixnumice, ixq, &
      itim_old, &
      ncol, lchnk                  ! # of columns, and chunk identifier

#endif
!BAS for new tropopause calculation but not used yet
    real(r8) :: calday

    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)

    real(r8) :: tropp_days(12)
!end BAS

!BAS
    logical :: first_step, first_restart_step
    integer :: nstep
!end BAS

  call t_startf('clubb_tend_cam')

!BAS for new tropopause calc but not used yet
  calday = get_curr_calday()

  do n = 1,12
    tropp_days(n) = get_calday( dates(n), 0 )
  end do
!end BAS

  do i = 1, pcols
    det_s(i)   = 0.0_r8
    det_ice(i) = 0.0_r8
  end do

#ifdef CLUBB_SGS

#ifdef _OPENACC
    ! These options have not been GPUized
    if ( clubb_l_ascending_grid ) call endrun(subr//': clubb_l_ascending_grid=.true. not available when compiling with OpenACC')
    if ( do_clubb_mf )            call endrun(subr//': do_clubb_mf=.true. not available when compiling with OpenACC')
    if ( do_rainturb )            call endrun(subr//': do_rainturb=.true. not available when compiling with OpenACC')
    if ( do_cldcool )             call endrun(subr//': do_cldcool=.true. not available when compiling with OpenACC')
    if ( single_column .and. .not. scm_cambfb_mode )  then
      call endrun(subr//': (single_column && !scm_cambfb_mode)=.true. not available when compiling with OpenACC')
    end if
#endif

    !-----------------------------------------------------------------------------------!
    !                           MAIN COMPUTATION BEGINS HERE                            !
    !-----------------------------------------------------------------------------------!

    call t_startf('clubb_tend_cam:non_acc_region')

    !  Get indicees for cloud and ice mass and cloud and ice number
    call cnst_get_ind('Q',ixq)
    call cnst_get_ind('CLDLIQ',ixcldliq)
    call cnst_get_ind('CLDICE',ixcldice)
    call cnst_get_ind('NUMLIQ',ixnumliq)
    call cnst_get_ind('NUMICE',ixnumice)

    !  Determine time step of physics buffer
    itim_old = pbuf_old_tim_idx()

    !  Establish associations between pointers and physics buffer fields
    call pbuf_get_field(pbuf, wp2_idx,        wp2_pbuf )
    call pbuf_get_field(pbuf, wp3_idx,        wp3_pbuf )
    call pbuf_get_field(pbuf, wpthlp_idx,     wpthlp_pbuf )
    call pbuf_get_field(pbuf, wprtp_idx,      wprtp_pbuf )
    call pbuf_get_field(pbuf, rtpthlp_idx,    rtpthlp_pbuf )
    call pbuf_get_field(pbuf, rtp2_idx,       rtp2_pbuf )
    call pbuf_get_field(pbuf, thlp2_idx,      thlp2_pbuf )
    call pbuf_get_field(pbuf, up2_idx,        up2_pbuf )
    call pbuf_get_field(pbuf, vp2_idx,        vp2_pbuf )

    call pbuf_get_field(pbuf, rtp3_idx,       rtp3_pbuf )
    call pbuf_get_field(pbuf, thlp3_idx,      thlp3_pbuf )
    call pbuf_get_field(pbuf, up3_idx,        up3_pbuf )
    call pbuf_get_field(pbuf, vp3_idx,        vp3_pbuf )

    call pbuf_get_field(pbuf, upwp_idx,       upwp_pbuf )
    call pbuf_get_field(pbuf, vpwp_idx,       vpwp_pbuf )
    call pbuf_get_field(pbuf, wpthvp_idx,     wpthvp_pbuf)
    call pbuf_get_field(pbuf, wp2thvp_idx,    wp2thvp_pbuf)
    call pbuf_get_field(pbuf, wp2up_idx,      wp2up_pbuf)
    call pbuf_get_field(pbuf, rtpthvp_idx,    rtpthvp_pbuf)
    call pbuf_get_field(pbuf, thlpthvp_idx,   thlpthvp_pbuf)

    call pbuf_get_field(pbuf, pdf_zm_w_1_idx,         pdf_zm_w_1_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_w_2_idx,         pdf_zm_w_2_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_1_idx,  pdf_zm_varnce_w_1_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_varnce_w_2_idx,  pdf_zm_varnce_w_2_pbuf )
    call pbuf_get_field(pbuf, pdf_zm_mixt_frac_idx,   pdf_zm_mixt_frac_pbuf )

    call pbuf_get_field(pbuf, wp2rtp_idx,       wp2rtp_pbuf     )
    call pbuf_get_field(pbuf, wp2thlp_idx,      wp2thlp_pbuf    )
    call pbuf_get_field(pbuf, uprcp_idx,        uprcp_pbuf      )
    call pbuf_get_field(pbuf, vprcp_idx,        vprcp_pbuf      )
    call pbuf_get_field(pbuf, rc_coef_zm_idx,   rc_coef_zm_pbuf )
    call pbuf_get_field(pbuf, wp4_idx,          wp4_pbuf        )
    call pbuf_get_field(pbuf, wpup2_idx,        wpup2_pbuf      )
    call pbuf_get_field(pbuf, wpvp2_idx,        wpvp2_pbuf      )
    call pbuf_get_field(pbuf, wp2up2_idx,       wp2up2_pbuf     )
    call pbuf_get_field(pbuf, wp2vp2_idx,       wp2vp2_pbuf     )

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

    call pbuf_get_field(pbuf, prer_evap_idx,      prer_evap_pbuf)
    call pbuf_get_field(pbuf, cmeliq_idx,         cmeliq_pbuf)
    call pbuf_get_field(pbuf, ice_supersat_idx,   ice_supersat_frac_pbuf)
    call pbuf_get_field(pbuf, relvar_idx,         relvar_pbuf)
    call pbuf_get_field(pbuf, dp_frac_idx,        deepcu_pbuf)
    call pbuf_get_field(pbuf, sh_frac_idx,        shalcu_pbuf)
    call pbuf_get_field(pbuf, kvh_idx,            khzm_pbuf)
    call pbuf_get_field(pbuf, pblh_idx,           pblh_pbuf)
    call pbuf_get_field(pbuf, icwmrdp_idx,        dp_icwmr_pbuf)
    call pbuf_get_field(pbuf, cmfmc_sh_idx,       cmfmc_sh_pbuf)

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

    !  Initialize physics tendency arrays
    call physics_ptend_init(ptend_all, state%psetcols, 'clubb')

    ! Copy the state to state_loc array to use in this routine
    call physics_state_copy(state, state_loc)

    ! Constituents are all treated as dry mmr by clubb.  Convert the water species to
    ! a dry basis.
    call set_wet_to_dry(state_loc, convert_cnst_type='wet')

    ! Define the grid box size.  CLUBB needs this information to determine what
    !  the maximum length scale should be.  This depends on the column for
    !  variable mesh grids and lat-lon grids
    call grid_size(state_loc, grid_dx, grid_dy)

    ! Determine number of columns and which chunk computation is to be performed on
    ncol = state%ncol
    lchnk = state%lchnk

    first_step = is_first_step()
    first_restart_step = is_first_restart_step()

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

!BAS does this make sense to do? 
    !$acc data copyin( state_loc, cam_in )

!BAS moved up from below to keep on CAM side
    call physics_ptend_init( ptend_loc, state%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq )

    call clubb1_run(ncol, pcols, lchnk, iam, nstep, state_loc%lat, state_loc%lon, hdtime, ztodtptr, &
                        pver, pverp, pcnst, clubb_timestep, gr, apply_const, &
                        nzt_clubb, nzm_clubb, sclr_dim, edsclr_dim, hydromet_dim, &
                        stats_metadata, hm_metadata, clubb_do_adv, first_step, first_restart_step, &
                        single_column, scm_cambfb_mode, scm_clubb_iop_name, &
                        shr_const_karman, shr_const_pi, shr_const_g, omega, theta0, &
                        macmic_it, top_lev, rtpthlp_const, wpthlp_const, wprtp_const, sclr_tol, &
                        ts_nudge, rtm_min, rtm_nudge_max_altitude, &
                        wp3_const, cld_macmic_num_steps, clubb_params_single_col, &
                        cpair, cpairv, invrs_cpairv, rair, rga, inv_p0_clubb, rairv, zvir, latvap, latice, &
                        gravit, clubb_rnevap_effic, do_cldcool, do_rainturb, &
                        do_clubb_mf, l_implemented, grid_type, lq, deep_scheme, &
                        state_loc%q, state_loc%u, state_loc%v, state_loc%t, state_loc%pmid, &
                        state_loc%zm, state_loc%phis, state_loc%pdel, state_loc%pdeldry, state_loc%s, &
                        state_loc%pint, state_loc%zi, state_loc%omega, wprcp, &
                        cam_in%wsx, cam_in%wsy, cam_in%cflx, cam_in%shf, cam_in%landfrac, &
                        ptend_loc%q, ptend_loc%u, ptend_loc%v, ptend_loc%s, &
                        sclr_idx, clubb_l_ascending_grid, clubb_do_energyfix, &
                        ixq, ixcldliq, ixcldice, ixrtpthlp, ixwpthlp, &
                        ixwprtp, ixwp3, ixwp2, ixthlp2, ixrtp2, ixup2, ixvp2, &
                        clubbtop_pbuf, clubb_l_intr_sfc_flux_smooth, &
                        pdf_params_chnk, pdf_params_zm_chnk, pdf_implicit_coefs_terms_chnk, &
                        clubb_config_flags, &
                        eleak, se_dis, rho_zm, rho_zt, exner, cloud_frac, &
                        zi_g, zt_g, &
                        grid_dx, grid_dy, &
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
                        mf_thlflx,  mf_qtflx, &
                        thlm, rtm, um, vm, wm_zt, rcm, rcm_in_layer, &
                        wp2_pbuf, wp3_pbuf, wpthlp_pbuf, wprtp_pbuf, &
                        rtpthlp_pbuf, rtp2_pbuf, thlp2_pbuf, rtp3_pbuf, &
                        thlp3_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, &
                        upwp_pbuf, vpwp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, &
                        rtpthvp_pbuf, thlpthvp_pbuf, pdf_zm_w_1_pbuf, pdf_zm_w_2_pbuf, &
                        pdf_zm_varnce_w_1_pbuf, pdf_zm_varnce_w_2_pbuf, pdf_zm_mixt_frac_pbuf, &
                        wp2rtp_pbuf, wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, &
                        wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, wp2up2_pbuf, wp2vp2_pbuf, cld_pbuf, &
                        concld_pbuf, ast_pbuf, alst_pbuf, aist_pbuf, qlst_pbuf, qist_pbuf, &
                        deepcu_pbuf, shalcu_pbuf, khzm_pbuf, pblh_pbuf, tke_pbuf, dp_icwmr_pbuf, &
                        ice_supersat_frac_pbuf, relvar_pbuf, naai_pbuf, cmeliq_pbuf, &
                        cmfmc_sh_pbuf, qsatfac_pbuf, npccn_pbuf, prer_evap_pbuf, qrl_pbuf, &
                        rtp2_mc_zt_pbuf, thlp2_mc_zt_pbuf, wprtp_mc_zt_pbuf, &
                        wpthlp_mc_zt_pbuf, rtpthlp_mc_zt_pbuf, ttend_clubb_pbuf, &
                        upwp_clubb_gw_pbuf, vpwp_clubb_gw_pbuf, thlp2_clubb_gw_pbuf, &
                        wpthlp_clubb_gw_pbuf, ttend_clubb_mc_pbuf, upwp_clubb_gw_mc_pbuf, &
                        vpwp_clubb_gw_mc_pbuf, thlp2_clubb_gw_mc_pbuf, wpthlp_clubb_gw_mc_pbuf, &
                        stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, &
                        out_zt, out_zm, out_sfc, out_radzt, out_radzm, &
                        errmsg, errflg )

    if (errflg /= 0) then
      call endrun(errmsg)
    end if

    do k = 1, pver
      do i = 1, ncol
        rvmtend_clubb_output(i,k) = ptend_loc%q(i,k,ixq)      * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        rcmtend_clubb_output(i,k) = ptend_loc%q(i,k,ixcldliq) * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        rimtend_clubb_output(i,k) = ptend_loc%q(i,k,ixcldice) * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        cmeliq_pbuf         (i,k) = ptend_loc%q(i,k,ixcldliq) * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        stend_clubb_output  (i,k) = ptend_loc%s(i,k)
        utend_clubb_output  (i,k) = ptend_loc%u(i,k)
        vtend_clubb_output  (i,k) = ptend_loc%v(i,k)
      end do
    end do

    !
    ! set pbuf field so that HB scheme is only applied above CLUBB top
    !
    if (do_hb_above_clubb) then
      call pbuf_set_field(pbuf, clubbtop_idx, clubbtop_pbuf)
    endif

    ! ------------------------------------------------- !
    ! End column computation of CLUBB, begin to apply   !
    ! and compute output, etc                           !
    ! ------------------------------------------------- !

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state_loc,ptend_loc,hdtime)

    ! ------------------------------------------------------------ !
    ! The rest of the code deals with diagnosing variables         !
    ! for microphysics/radiation computation and macrophysics      !
    ! ------------------------------------------------------------ !

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
    !  Detrainment of convective condensate into the environment or stratiform cloud    !
    ! --------------------------------------------------------------------------------- !

    lqice(:)        = .false.
    lqice(ixcldliq) = .true.
    lqice(ixcldice) = .true.
    lqice(ixnumliq) = .true.
    lqice(ixnumice) = .true.

    call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lq=lqice)

    call clubb2_run(ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice, & ! in
                    clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp, &! in
                    meltpt_temp, latice, rga, & ! in
                    dlf, state_loc%t, state_loc%pdel, state_loc%pdeldry, & ! in
                    ptend_loc%q, ptend_loc%s, det_s, det_ice, & ! inout
                    dlf_liq_out, dlf_ice_out ) ! out

    do k = 1, pver
      do i = 1, ncol
        dpdlfliq_output(i,k)    = ptend_loc%q(i,k,ixcldliq) * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        dpdlfice_output(i,k)    = ptend_loc%q(i,k,ixcldice) * state_loc%pdeldry(i,k) / state_loc%pdel(i,k)
        dpdlft_output(i,k)      = ptend_loc%s(i,k) * invrs_cpairv(i,k)
        detnliquid_output(i,k)  = ptend_loc%q(i,k,ixnumliq)
      end do
    end do

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state_loc,ptend_loc,hdtime)

    !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
    troplev(:) = 0
    !REMOVECAM_END
    call tropopause_findChemTrop( state, troplev )

    aist_pbuf(:,:top_lev-1) = 0._r8
    qsatfac_pbuf(:, :) = 0._r8

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
        call aist_vector(state_loc%q(:,k,ixq),state_loc%t(:,k),state_loc%pmid(:,k),state_loc%q(:,k,ixcldice), &
             state_loc%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist_pbuf(:,k),ncol )
      else
        call aist_vector(state_loc%q(:,k,ixq),state_loc%t(:,k),state_loc%pmid(:,k),state_loc%q(:,k,ixcldice), &
              state_loc%q(:,k,ixnumice), cam_in%landfrac(:),cam_in%snowhland(:),aist_pbuf(:,k),ncol,&
              qsatfac_out=qsatfac_pbuf(:,k), rhmini_in=rhmini, rhmaxi_in=rhmaxi)
      endif
    enddo

    call clubb3_run(pcols, ncol, pver, pverp, pcnst, top_lev, & ! in
                    ixq, ixcldice, ixcldliq, ixnumice, & ! in
                    rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const, & ! in
                    dp1, dp2, zvir, rair, cpair, gravit, karman, & ! in
                    calday, tropp_days, & ! in
                    state_loc%lat, state_loc%phis, cam_in%landfrac, cam_in%snowhland, & ! in
                    cam_in%wsx, cam_in%wsy, cam_in%shf, & ! in
                    state_loc%pint, state_loc%pmid, state_loc%pdel, state_loc%pdeldry, & ! in
                    rcm, cloud_frac, state_loc%t, exner, & ! in
                    state_loc%exner, state_loc%zm, state_loc%zi, state_loc%u, & ! in
                    state_loc%v, cmfmc, cam_in%cflx, state_loc%q, & ! in
                    single_column, scm_cambfb_mode, lq, & ! in
                    cnst_type, scm_clubb_iop_name, subcol_scheme, & ! in
                    pblh_pbuf, alst_pbuf, qlst_pbuf, deepcu_pbuf, shalcu_pbuf, & ! inout
                    cmfmc_sh_pbuf, dp_icwmr_pbuf, concld_pbuf, aist_pbuf,  & ! inout
                    qsatfac_pbuf, ast_pbuf, qist_pbuf, cld_pbuf, ptend_all%q, troplev, & ! inout
                    errmsg, errflg ) ! out

    if (errflg /= 0) then
      call endrun(errmsg)
    end if

    !----------------------------------------- Output section -----------------------------------------

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

    do k = top_lev, pverp
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        zi_output(i,k)              =        zi_g(i,k_clubb)
        wp2_output(i,k)             =    wp2_pbuf(i,k_clubb)
        up2_output(i,k)             =    up2_pbuf(i,k_clubb)
        vp2_output(i,k)             =    vp2_pbuf(i,k_clubb)
        upwp_output(i,k)            =   upwp_pbuf(i,k_clubb)
        vpwp_output(i,k)            =   vpwp_pbuf(i,k_clubb)
        rtp2_output(i,k)            =   rtp2_pbuf(i,k_clubb)
        wprcp_clubb_output(i,k)     =       wprcp(i,k_clubb) * latvap
        wpthvp_clubb_output(i,k)    = wpthvp_pbuf(i,k_clubb) * cpair
        thlp2_output(i,k)           =  thlp2_pbuf(i,k_clubb)

        wpthlp_output(i,k)  = ( wpthlp_pbuf(i,k_clubb) - (apply_const *  wpthlp_const) ) &
                              * rho_zm(i,k_clubb) * cpair !  liquid water potential temperature flux

        wprtp_output(i,k)   = (  wprtp_pbuf(i,k_clubb) - (apply_const *   wprtp_const) ) &
                              * rho_zm(i,k_clubb) * latvap  !  total water mixig ratio flux

        rtpthlp_output(i,k) =  rtpthlp_pbuf(i,k_clubb) - (apply_const * rtpthlp_const)

      end do
    end do

    ! Convert RTP2 and THLP2 to thermo grid for output
    rtp2_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr,  rtp2_pbuf(:ncol,:) )
    thl2_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr, thlp2_pbuf(:ncol,:) )
    wp2_zt  = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr,   wp2_pbuf(:ncol,:) )

    do k = top_lev, pver
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        rho_output(i,k)               = rho_zt(i,k_clubb)
        rcm_output(i,k)               = rcm(i,k_clubb)
        rtm_output(i,k)               = rtm(i,k_clubb)
        thlm_output(i,k)              = thlm(i,k_clubb)
        um_output(i,k)                = um(i,k_clubb)
        vm_output(i,k)                = vm(i,k_clubb)
        rcm_in_layer_output(i,k)      = rcm_in_layer(i,k_clubb)
        zt_output(i,k)                = zt_g(i,k_clubb)
        wm_zt_output(i,k)             = wm_zt(i,k_clubb)
        rtp2_zt_output(i,k)           = rtp2_zt(i,k_clubb)
        thl2_zt_output(i,k)           = thl2_zt(i,k_clubb)
        wp2_zt_output(i,k)            = wp2_zt(i,k_clubb)
        wp3_output(i,k)               = wp3_pbuf(i,k_clubb) - (apply_const*wp3_const)

      end do
    end do

    do k = 1, nzt_clubb
      do i = 1, ncol

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
        rho_output(i,k)             = 0._r8
        wp2_output(i,k)             = 0._r8
        up2_output(i,k)             = 0._r8
        vp2_output(i,k)             = 0._r8
        rtp2_output(i,k)            = 0._r8
        thlp2_output(i,k)           = 0._r8
        zt_output(i,k)              = 0._r8
        rtp2_zt_output(i,k)         = 0._r8
        wp3_output(i,k)             = 0._r8
        thl2_zt_output(i,k)         = 0._r8
        wp2_zt_output(i,k)          = 0._r8
        rcm_in_layer_output(i,k)    = 0._r8
        pdfp_rtp2_output(i,k)       = 0._r8
        wm_zt_output(i,k)           = 0._r8
        rcm_output(i,k)             = 0._r8
        rtm_output(i,k)             = 0._r8
        thlm_output(i,k)            = 0._r8
        um_output(i,k)              = 0._r8
        vm_output(i,k)              = 0._r8
        zi_output(i,k)              = 0._r8
        wpthlp_output(i,k)          = 0._r8
        rtpthlp_output(i,k)         = 0._r8
        wprtp_output(i,k)           = 0._r8
        upwp_output(i,k)            = 0._r8
        vpwp_output(i,k)            = 0._r8
        wprcp_clubb_output(i,k)     = 0._r8
        wpthvp_clubb_output(i,k)    = 0._r8
      end do
    end do

    !  Output calls of variables goes here
    call outfld( 'WP2_CLUBB',        wp2_output,                     pcols, lchnk )
    call outfld( 'UP2_CLUBB',        up2_output,                     pcols, lchnk )
    call outfld( 'VP2_CLUBB',        vp2_output,                     pcols, lchnk )
    call outfld( 'WP3_CLUBB',        wp3_output,                     pcols, lchnk )
    call outfld( 'UPWP_CLUBB',       upwp_output,                    pcols, lchnk )
    call outfld( 'VPWP_CLUBB',       vpwp_output,                    pcols, lchnk )
    call outfld( 'WPTHLP_CLUBB',     wpthlp_output,                  pcols, lchnk )
    call outfld( 'WPRTP_CLUBB',      wprtp_output,                   pcols, lchnk )
    call outfld( 'RTP2_CLUBB',       rtp2_output,                    pcols, lchnk )
    call outfld( 'RTPTHLP_CLUBB',    rtpthlp_output,                 pcols, lchnk )
    call outfld( 'RCM_CLUBB',        rcm_output,                     pcols, lchnk )
    call outfld( 'RTM_CLUBB',        rtm_output,                     pcols, lchnk )
    call outfld( 'THLM_CLUBB',       thlm_output,                    pcols, lchnk )
    call outfld( 'WPRCP_CLUBB',      wprcp_clubb_output,             pcols, lchnk )
    call outfld( 'WPTHVP_CLUBB',     wpthvp_clubb_output,            pcols, lchnk )
    call outfld( 'RTP2_ZT_CLUBB',    rtp2_zt_output,                 pcols, lchnk )
    call outfld( 'THLP2_ZT_CLUBB',   thl2_zt_output,                 pcols, lchnk )
    call outfld( 'WP2_ZT_CLUBB',     wp2_zt_output,                  pcols, lchnk )
    call outfld( 'PDFP_RTP2_CLUBB',  pdfp_rtp2_output,               pcols, lchnk )
    call outfld( 'THLP2_CLUBB',      thlp2_output,                   pcols, lchnk )
    call outfld( 'RCMINLAYER_CLUBB', rcm_in_layer_output,            pcols, lchnk )
    call outfld( 'ZT_CLUBB',         zt_output,                      pcols, lchnk )
    call outfld( 'ZM_CLUBB',         zi_output,                      pcols, lchnk )
    call outfld( 'UM_CLUBB',         um_output,                      pcols, lchnk )
    call outfld( 'VM_CLUBB',         vm_output,                      pcols, lchnk )
    call outfld( 'WM_ZT_CLUBB',      wm_zt_output,                   pcols, lchnk )
    call outfld( 'RHO_CLUBB',        rho_output,                     pcols, lchnk )

    call outfld( 'RELVAR',           relvar_pbuf,                    pcols, lchnk )
    call outfld( 'CLOUDCOVER_CLUBB', cld_pbuf,                       pcols, lchnk )
    call outfld( 'CLOUDFRAC_CLUBB',  alst_pbuf,                      pcols, lchnk )
    call outfld( 'CONCLD',           concld_pbuf,                    pcols, lchnk )
    call outfld( 'DP_CLD',           deepcu_pbuf,                    pcols, lchnk )
    call outfld( 'ZMDLF',            dlf_liq_out,                    pcols, lchnk )
    call outfld( 'ZMDLFI',           dlf_ice_out,                    pcols, lchnk )
    call outfld( 'CLUBB_GRID_SIZE',  grid_dx,                        pcols, lchnk )
    call outfld( 'QSATFAC',          qsatfac_pbuf,                   pcols, lchnk )


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
          mf_thlflx_output(i,k)    = mf_thlflx(i,k_clubb) * rho_zm(i,k_clubb) * cpair
          mf_qtflx_output(i,k)     = mf_qtflx(i,k_clubb) * rho_zm(i,k_clubb) * latvap
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

      do j = 1, stats_zt(1)%num_output_fields

        temp1 = trim(stats_zt(1)%file%grid_avg_var(j)%name)
        sub   = temp1
        if (len(temp1) >  max_fieldname_len) sub = temp1(1:max_fieldname_len)

        call outfld(trim(sub), out_zt(:,:,j), pcols, lchnk )
      enddo

      do j = 1, stats_zm(1)%num_output_fields

        temp1 = trim(stats_zm(1)%file%grid_avg_var(j)%name)
        sub   = temp1
        if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)

        call outfld(trim(sub),out_zm(:,:,j), pcols, lchnk)
      enddo

      if (stats_metadata%l_output_rad_files) then
        do j = 1, stats_rad_zt(1)%num_output_fields
          call outfld(trim(stats_rad_zt(1)%file%grid_avg_var(j)%name), out_radzt(:,:,j), pcols, lchnk)
        enddo

        do j = 1, stats_rad_zm(1)%num_output_fields
          call outfld(trim(stats_rad_zm(1)%file%grid_avg_var(j)%name), out_radzm(:,:,j), pcols, lchnk)
        enddo
      endif

      do j = 1, stats_sfc(1)%num_output_fields
        call outfld(trim(stats_sfc(1)%file%grid_avg_var(j)%name), out_sfc(:,:,j), pcols, lchnk)
      enddo

    endif
    call t_stopf('clubb_tend_cam:non_acc_region')

    ! Cleanup err_info
    call cleanup_err_info_api(err_info)
#endif

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
  do m = 2, pcnst
    ptend%q(:ncol,pver,m) = cam_in%cflx(:ncol,m)*state%rpdel(:ncol,pver)*gravit
  end do

  ! Convert tendencies of dry constituents to dry basis.
  do m = 2, pcnst
     if (cnst_type(m).eq.'dry') then
        ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m)*state%pdel(:ncol,pver)*state%rpdeldry(:ncol,pver)
     endif
  end do

  end subroutine clubb_emissions_cam

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
#ifdef CLUBB_SGS
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
  do i = 1, state%ncol
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
