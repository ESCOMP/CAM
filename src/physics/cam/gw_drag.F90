module gw_drag

!--------------------------------------------------------------------------
! CAM and WACCM gravity wave parameterizations were merged by Sean Patrick
! Santos in Summer 2013, and at the same time, gw_drag was split into
! various modules. This is the CAM interface and driver module. The below
! notes are for the old CAM and WACCM versions of gw_drag.
!--------------------------------------------------------------------------
! This file came from wa17 and was modified by Fabrizio: 07-02-2004
! Standard gw_drag with modification (6) of latitude profile of gw spectrum
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the forcing due to parameterized gravity waves. Both an
! orographic and an internal source spectrum are considered.
!
! Author: Byron Boville
!
!--------------------------------------------------------------------------
  use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
  use shr_log_mod,    only: errMsg => shr_log_errMsg
  use shr_assert_mod, only: shr_assert

  use ppgrid,         only: pcols, pver, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,     only: masterproc
  use cam_history,    only: outfld
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  use ref_pres,       only: do_molec_diff, nbot_molec, press_lim_idx
  use physconst,      only: cpair

  ! These are the actual switches for different gravity wave sources.
  use phys_control,   only: use_gw_oro, use_gw_front, use_gw_front_igw, &
                            use_gw_convect_dp, use_gw_convect_sh,       &
                            use_simple_phys

  use gw_common,      only: GWBand
  use gw_convect,     only: BeresSourceDesc
  use gw_front,       only: CMSourceDesc

! Typical module header
  implicit none
  private
  save

!
! PUBLIC: interfaces
!
  public :: gw_drag_readnl           ! Read namelist
  public :: gw_init                  ! Initialization
  public :: gw_tend                  ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  real(r8), parameter :: unset_r8 = huge(1._r8)

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro
  ! Medium scale waves.
  type(GWBand) :: band_mid
  ! Long scale waves for IGWs.
  type(GWBand) :: band_long

  ! Top level for gravity waves.
  integer, parameter :: ktop = 1
  ! Bottom level for frontal waves.
  integer :: kbot_front

  ! Factor for SH orographic waves.
  real(r8) :: gw_oro_south_fac = 1._r8
  
  ! Frontogenesis function critical threshold.
  real(r8) :: frontgfc = unset_r8

  ! Tendency efficiencies.

  ! Ridge scheme.
  logical  :: use_gw_rdg_beta      = .false.
  integer  :: n_rdg_beta           = 1
  real(r8) :: effgw_rdg_beta       = unset_r8
  real(r8) :: effgw_rdg_beta_max   = unset_r8
  real(r8) :: rdg_beta_cd_llb      = unset_r8  ! Low-level obstacle drag coefficient Ridge scheme.
  logical  :: trpd_leewv_rdg_beta  = .false.

  logical  :: use_gw_rdg_gamma     = .false.
  integer  :: n_rdg_gamma          = -1
  real(r8) :: effgw_rdg_gamma      = unset_r8
  real(r8) :: effgw_rdg_gamma_max  = unset_r8
  real(r8) :: rdg_gamma_cd_llb     = unset_r8
  logical  :: trpd_leewv_rdg_gamma = .false.
  character(len=cl) :: bnd_rdggm   = 'bnd_rdggm' ! full pathname for meso-Gamma ridge dataset

  ! Orography.
  real(r8) :: effgw_oro = unset_r8
  ! C&M scheme.
  real(r8) :: effgw_cm = unset_r8
  ! C&M scheme (inertial waves).
  real(r8) :: effgw_cm_igw = unset_r8
  ! Beres (deep convection).
  real(r8) :: effgw_beres_dp = unset_r8
  ! Beres (shallow convection).
  real(r8) :: effgw_beres_sh = unset_r8

  ! Horzontal wavelengths [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  real(r8), parameter :: wavelength_long = 3.e5_r8

  ! Background stress source strengths.
  real(r8) :: taubgnd = unset_r8
  real(r8) :: taubgnd_igw = unset_r8

  ! Whether or not to use a polar taper for frontally generated waves.
  logical :: gw_polar_taper = .false.

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical :: tau_0_ubc = .false.

  ! Whether or not to limit tau *before* applying any efficiency factors.
  logical :: gw_limit_tau_without_eff = .false.

  ! Whether or not to apply tendency max
  logical :: gw_apply_tndmax = .true.

  ! Files to read Beres source spectra from.
  character(len=256) :: gw_drag_file = ""
  character(len=256) :: gw_drag_file_sh = ""

  ! Beres settings and table.
  type(BeresSourceDesc) :: beres_dp_desc
  type(BeresSourceDesc) :: beres_sh_desc

  ! Width of gaussian used to create frontogenesis tau profile [m/s].
  real(r8), parameter :: front_gaussian_width = 30._r8

  ! Frontogenesis wave settings.
  type(CMSourceDesc) :: cm_desc
  type(CMSourceDesc) :: cm_igw_desc

  ! Indices into pbuf
  integer :: kvt_idx      = -1
  integer :: ttend_dp_idx = -1
  integer :: ttend_sh_idx = -1
  integer :: frontgf_idx  = -1
  integer :: frontga_idx  = -1
  integer :: sgh_idx      = -1

  ! anisotropic ridge fields
  integer, parameter :: prdg = 16

  real(r8), allocatable, dimension(:,:), target :: &
     rdg_gbxar

     ! Meso Beta
  real(r8), allocatable, dimension(:,:,:), target :: &
     rdg_hwdth,  &
     rdg_clngt,  &
     rdg_mxdis,  &
     rdg_anixy,  &
     rdg_angll

     ! Meso Gamma
  real(r8), allocatable, dimension(:,:,:), target :: &
     rdg_hwdthg, &
     rdg_clngtg, &
     rdg_mxdisg, &
     rdg_anixyg, &
     rdg_angllg

  ! Water constituent indices for budget
  integer :: ixcldliq = -1
  integer :: ixcldice = -1

  ! Prefixes for history field names
  character(len=1), parameter :: cm_pf = " "
  character(len=1), parameter :: cm_igw_pf = "I"
  character(len=1), parameter :: beres_dp_pf = "B"
  character(len=1), parameter :: beres_sh_pf = "S"

  ! namelist 
  logical          :: history_amwg                   ! output the variables used by the AMWG diag package
  logical  :: gw_lndscl_sgh = .true. ! scale SGH by land frac
  real(r8) :: gw_prndl = 0.25_r8
  real(r8) :: gw_qbo_hdepth_scaling = 1._r8 ! heating depth scaling factor

  logical :: gw_top_taper=.false.
  real(r8), pointer :: vramp(:)=>null()

!==========================================================================
contains
!==========================================================================

subroutine gw_drag_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, &
                             mpi_character, mpi_logical, mpi_integer
  use gw_rdg,          only: gw_rdg_readnl

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_drag_readnl'

  ! Maximum wave number and width of spectrum bins.
  integer :: pgwv = -1
  real(r8) :: gw_dc = unset_r8
  integer :: pgwv_long = -1
  real(r8) :: gw_dc_long = unset_r8

  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = unset_r8   ! critical froude number squared

  namelist /gw_drag_nl/ pgwv, gw_dc, pgwv_long, gw_dc_long, tau_0_ubc, &
       effgw_beres_dp, effgw_beres_sh, effgw_cm, effgw_cm_igw, effgw_oro, &
       fcrit2, frontgfc, gw_drag_file, gw_drag_file_sh, taubgnd, &
       taubgnd_igw, gw_polar_taper, &
       use_gw_rdg_beta, n_rdg_beta, effgw_rdg_beta, effgw_rdg_beta_max, &
       rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
       use_gw_rdg_gamma, n_rdg_gamma, effgw_rdg_gamma, effgw_rdg_gamma_max, &
       rdg_gamma_cd_llb, trpd_leewv_rdg_gamma, bnd_rdggm, &
       gw_oro_south_fac, gw_limit_tau_without_eff, &
       gw_lndscl_sgh, gw_prndl, gw_apply_tndmax, gw_qbo_hdepth_scaling, &
       gw_top_taper
  !----------------------------------------------------------------------

  if (use_simple_phys) return

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'gw_drag_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, gw_drag_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(pgwv, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: pgwv")
  call mpi_bcast(gw_dc, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_dc")
  call mpi_bcast(pgwv_long, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: pgwv_long")
  call mpi_bcast(gw_dc_long, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_dc_long")
  call mpi_bcast(tau_0_ubc, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: tau_0_ubc")
  call mpi_bcast(effgw_beres_dp, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_beres_dp")
  call mpi_bcast(effgw_beres_sh, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_beres_sh")
  call mpi_bcast(effgw_cm, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_cm")
  call mpi_bcast(effgw_cm_igw, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_cm_igw")
  call mpi_bcast(effgw_oro, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_oro")

  call mpi_bcast(use_gw_rdg_beta, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_gw_rdg_beta")
  call mpi_bcast(n_rdg_beta, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: n_rdg_beta")
  call mpi_bcast(effgw_rdg_beta, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_beta")
  call mpi_bcast(effgw_rdg_beta_max, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_beta_max")
  call mpi_bcast(rdg_beta_cd_llb, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rdg_beta_cd_llb")
  call mpi_bcast(trpd_leewv_rdg_beta, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: trpd_leewv_rdg_beta")

  call mpi_bcast(use_gw_rdg_gamma, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_gw_rdg_gamma")
  call mpi_bcast(n_rdg_gamma, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: n_rdg_gamma")
  call mpi_bcast(effgw_rdg_gamma, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_gamma")
  call mpi_bcast(effgw_rdg_gamma_max, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_gamma_max")
  call mpi_bcast(rdg_gamma_cd_llb, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rdg_gamma_cd_llb")
  call mpi_bcast(trpd_leewv_rdg_gamma, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: trpd_leewv_rdg_gamma")
  call mpi_bcast(bnd_rdggm, len(bnd_rdggm), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: bnd_rdggm")

  call mpi_bcast(gw_oro_south_fac, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_oro_south_fac")
  call mpi_bcast(fcrit2, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fcrit2")
  call mpi_bcast(frontgfc, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frontgfc")
  call mpi_bcast(taubgnd, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: taubgnd")
  call mpi_bcast(taubgnd_igw, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: taubgnd_igw")

  call mpi_bcast(gw_polar_taper, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_polar_taper")
  call mpi_bcast(gw_limit_tau_without_eff, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_limit_tau_without_eff")
  call mpi_bcast(gw_apply_tndmax, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_apply_tndmax")

  call mpi_bcast(gw_top_taper, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_top_taper")

  call mpi_bcast(gw_lndscl_sgh, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_lndscl_sgh")
  call mpi_bcast(gw_prndl, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_prndl")
  call mpi_bcast(gw_qbo_hdepth_scaling, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_qbo_hdepth_scaling")

  call mpi_bcast(gw_drag_file, len(gw_drag_file), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file")
  call mpi_bcast(gw_drag_file_sh, len(gw_drag_file_sh), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file_sh")


  ! Check if fcrit2 was set.
  call shr_assert(fcrit2 /= unset_r8, &
       "gw_drag_readnl: fcrit2 must be set via the namelist."// &
       errMsg(__FILE__, __LINE__))

  ! Check if pgwv was set.
  call shr_assert(pgwv >= 0, &
       "gw_drag_readnl: pgwv must be set via the namelist and &
       &non-negative."// &
       errMsg(__FILE__, __LINE__))

  ! Check if gw_dc was set.
  call shr_assert(gw_dc /= unset_r8, &
       "gw_drag_readnl: gw_dc must be set via the namelist."// &
       errMsg(__FILE__, __LINE__))

  band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)
  band_mid = GWBand(pgwv, gw_dc, 1.0_r8, wavelength_mid)
  band_long = GWBand(pgwv_long, gw_dc_long, 1.0_r8, wavelength_long)

  if (use_gw_rdg_gamma .or. use_gw_rdg_beta) then
     call gw_rdg_readnl(nlfile)
  end if

end subroutine gw_drag_readnl

!==========================================================================

subroutine gw_init()
  !-----------------------------------------------------------------------
  ! Time independent initialization for multiple gravity wave
  ! parameterization.
  !-----------------------------------------------------------------------

  use cam_history,      only: addfld, add_default, horiz_only
  use cam_history,      only: register_vector_field
  use interpolate_data, only: lininterp
  use phys_control,     only: phys_getopts
  use physics_buffer,   only: pbuf_get_index
  use constituents,     only: cnst_get_ind

  use cam_initfiles,    only: topo_file_get_id

  ! temporary for restart with ridge scheme
  use cam_initfiles,    only: bnd_topo
  
  use cam_pio_utils,    only: cam_pio_openfile
  use cam_grid_support, only: cam_grid_check, cam_grid_id
  use cam_grid_support, only: cam_grid_get_dim_names
  use pio,              only: file_desc_t, pio_nowrite, pio_closefile
  use ncdio_atm,        only: infld
  use ioFileMod,        only: getfil

  use ref_pres,   only: pref_edge, pref_mid
  use physconst,  only: gravit, rair, rearth

  use gw_common,  only: gw_common_init
  use gw_front,   only: gaussian_cm_desc

  !---------------------------Local storage-------------------------------

  integer          :: i, l, k
  character(len=1) :: cn

  ! Index for levels at specific pressures.
  integer :: kfront

  ! output tendencies and state variables for CAM4 temperature,
  ! water vapor, cloud ice and cloud liquid budgets.
  logical :: history_budget
  ! output history file number for budget fields
  integer :: history_budget_histfile_num
  ! output variables of interest in WACCM runs
  logical :: history_waccm

  ! Interpolated Newtonian cooling coefficients.
  real(r8) :: alpha(pver+1)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer, parameter :: nalph = 71
  real(r8) :: alpha0(nalph) = [ &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.10133333_r8,  0.104_r8,       &
       0.108_r8,       0.112_r8,       0.116_r8,       0.12066667_r8,  &
       0.126_r8,       0.132_r8,       0.138_r8,       0.144_r8,       &
       0.15133333_r8,  0.16_r8,        0.17_r8,        0.18_r8,        &
       0.19_r8,        0.19933333_r8,  0.208_r8,       0.216_r8,       &
       0.224_r8,       0.232_r8,       0.23466667_r8,  0.232_r8,       &
       0.224_r8,       0.216_r8,       0.208_r8,       0.20133333_r8,  &
       0.196_r8,       0.192_r8,       0.188_r8,       0.184_r8,       &
       0.18266667_r8,  0.184_r8,       0.188_r8,       0.192_r8,       &
       0.196_r8,       0.19333333_r8,  0.184_r8,       0.168_r8,       &
       0.152_r8,       0.136_r8,       0.12133333_r8,  0.108_r8,       &
       0.096_r8,       0.084_r8,       0.072_r8,       0.061_r8,       &
       0.051_r8,       0.042_r8,       0.033_r8,       0.024_r8,       &
       0.017666667_r8, 0.014_r8,       0.013_r8,       0.012_r8,       &
       0.011_r8,       0.010333333_r8, 0.01_r8,        0.01_r8,        &
       0.01_r8,        0.01_r8,        0.01_r8                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(r8) :: palph(nalph) = [ &
       2.06115E-06_r8, 2.74280E-06_r8, 3.64988E-06_r8, 4.85694E-06_r8, &
       6.46319E-06_r8, 8.60065E-06_r8, 1.14450E-05_r8, 1.52300E-05_r8, &
       2.02667E-05_r8, 2.69692E-05_r8, 3.58882E-05_r8, 4.77568E-05_r8, &
       6.35507E-05_r8, 8.45676E-05_r8, 0.000112535_r8, 0.000149752_r8, &
       0.000199277_r8, 0.000265180_r8, 0.000352878_r8, 0.000469579_r8, &
       0.000624875_r8, 0.000831529_r8, 0.00110653_r8,  0.00147247_r8,  &
       0.00195943_r8,  0.00260744_r8,  0.00346975_r8,  0.00461724_r8,  &
       0.00614421_r8,  0.00817618_r8,  0.0108801_r8,   0.0144783_r8,   &
       0.0192665_r8,   0.0256382_r8,   0.0341170_r8,   0.0453999_r8,   &
       0.0604142_r8,   0.0803939_r8,   0.106981_r8,    0.142361_r8,    &
       0.189442_r8,    0.252093_r8,    0.335463_r8,    0.446404_r8,    &
       0.594036_r8,    0.790490_r8,    1.05192_r8,     1.39980_r8,     &
       1.86273_r8,     2.47875_r8,     3.29851_r8,     4.38936_r8,     &
       5.84098_r8,     7.77266_r8,     10.3432_r8,     13.7638_r8,     &
       18.3156_r8,     24.3728_r8,     32.4332_r8,     43.1593_r8,     &
       57.4326_r8,     76.4263_r8,     101.701_r8,     135.335_r8,     &
       180.092_r8,     239.651_r8,     318.907_r8,     424.373_r8,     &
       564.718_r8,     751.477_r8,     1000._r8                        &
       ]

  ! Read data from file
  type(file_desc_t), pointer :: fh_topo
  type(file_desc_t)  :: fh_rdggm
  integer            :: grid_id
  character(len=8)   :: dim1name, dim2name
  logical            :: found
  character(len=256) :: bnd_rdggm_loc   ! filepath of topo file on local disk

  ! Allow reporting of error messages.
  character(len=128) :: errstring
  character(len=*), parameter :: sub = 'gw_init'

  ! temporary workaround for restart w/ ridge scheme
  character(len=256) :: bnd_topo_loc   ! filepath of topo file on local disk

  integer :: botndx,topndx

  !-----------------------------------------------------------------------

  if (do_molec_diff) then
     kvt_idx     = pbuf_get_index('kvt')
  end if

  if (masterproc) then
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_mid%ngwv = ", band_mid%ngwv
     do l = -band_mid%ngwv, band_mid%ngwv
        write (iulog,'(A,I0,A,F7.2)') &
             "GW_DRAG: band_mid%cref(",l,") = ",band_mid%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_mid%kwv = ', band_mid%kwv
     write(iulog,*) 'GW_DRAG: band_mid%fcrit2 = ', band_mid%fcrit2
     write(iulog,*) ' '
     write(iulog,*) "GW_DRAG: band_long%ngwv = ", band_long%ngwv
     do l = -band_long%ngwv, band_long%ngwv
        write (iulog,'(A,I2,A,F7.2)') &
             "GW_DRAG: band_long%cref(",l,") = ",band_long%cref(l)
     enddo
     write(iulog,*) 'GW_DRAG: band_long%kwv = ', band_long%kwv
     write(iulog,*) 'GW_DRAG: band_long%fcrit2 = ', band_long%fcrit2
     write(iulog,*) ' '
  end if

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa

  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._r8
     alpha0(k) = max(alpha0(k), 1.e-6_r8)
     palph(k) = palph(k)*1.e2_r8
  end do

  ! interpolate to current vertical grid and obtain alpha

  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pver+1)
  if (masterproc) then
     write (iulog,*) 'gw_init: newtonian damping (1/day):'
     write (iulog,fmt='(a4,a12,a10)') ' k  ','  pref_edge      ', &
          '  alpha   '
     do k = 1, pver+1
        write (iulog,fmt='(i4,1e12.5,1f10.2)') k,pref_edge(k), &
             alpha(k)*86400._r8
     end do
  end if

  if (masterproc) then
     write(iulog,*) 'KTOP        =',ktop
  end if

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm, &
       history_amwg_out   = history_amwg  )

  ! Initialize subordinate modules.
  call gw_common_init(pver,&
       tau_0_ubc, ktop, gravit, rair, alpha, gw_prndl, & 
       gw_qbo_hdepth_scaling, errstring )
  call shr_assert(trim(errstring) == "", "gw_common_init: "//errstring// &
       errMsg(__FILE__, __LINE__))

  if ( use_gw_oro ) then

     if (effgw_oro == unset_r8) then
        call endrun("gw_drag_init: Orographic gravity waves enabled, &
             &but effgw_oro was not set.")
     end if
  end if
     
  if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then

     sgh_idx = pbuf_get_index('SGH')

     ! Declare history variables for orographic term
     call addfld ('TAUAORO',    (/ 'ilev' /), 'I','N/m2',  &
          'Total stress from original OGW scheme')
     call addfld ('TTGWORO',    (/ 'lev' /), 'A','K/s',  &
          'T tendency - orographic gravity wave drag')
     call addfld ('TTGWSDFORO', (/ 'lev' /), 'A','K/s',  &
          'T tendency - orographic gravity wave, diffusion.')
     call addfld ('TTGWSKEORO', (/ 'lev' /), 'A','K/s',  &
          'T tendency - orographic gravity wave, breaking KE.')
     call addfld ('UTGWORO',    (/ 'lev' /), 'A','m/s2', &
          'U tendency - orographic gravity wave drag')
     call addfld ('VTGWORO',    (/ 'lev' /), 'A','m/s2', &
          'V tendency - orographic gravity wave drag')
     call register_vector_field('UTGWORO', 'VTGWORO')
     call addfld ('TAUGWX',     horiz_only,  'A','N/m2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAUGWY',     horiz_only,  'A','N/m2', &
          'Meridional gravity wave surface stress')
     call register_vector_field('TAUGWX', 'TAUGWY')

     if (history_amwg) then
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

     if (history_waccm) then
        call add_default('UTGWORO ', 1, ' ')
        call add_default('VTGWORO ', 1, ' ')
        call add_default('TAUGWX  ', 1, ' ')
        call add_default('TAUGWY  ', 1, ' ')
     end if

  end if

  if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     grid_id = cam_grid_id('physgrid')
     if (.not. cam_grid_check(grid_id)) then
        call endrun(sub//': ERROR: no "physgrid" grid')
     end if
     call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
  end if

  if (use_gw_rdg_beta) then

     if (effgw_rdg_beta == unset_r8) then
        call endrun(sub//": ERROR: Anisotropic OGW enabled, &
                     &but effgw_rdg_beta was not set.")
     end if

     fh_topo => topo_file_get_id()
     bnd_topo_loc = ' '
     if (.not. associated(fh_topo)) then

        ! Try to open topo file here.  This workaround will not be needed
        ! once the refactored initialization sequence is on trunk.
        
        allocate(fh_topo)
        ! Error exit is from getfil if file not found.
        call getfil(bnd_topo, bnd_topo_loc)
        call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)

     end if

     ! Get beta ridge data
     allocate( &
        rdg_gbxar(pcols,begchunk:endchunk),      &
        rdg_hwdth(pcols,prdg,begchunk:endchunk), &
        rdg_clngt(pcols,prdg,begchunk:endchunk), &
        rdg_mxdis(pcols,prdg,begchunk:endchunk), &
        rdg_anixy(pcols,prdg,begchunk:endchunk), &
        rdg_angll(pcols,prdg,begchunk:endchunk)  )
   
     call infld('GBXAR', fh_topo, dim1name, dim2name, 1, pcols, &
                         begchunk, endchunk, rdg_gbxar, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: GBXAR not found on topo file')
     rdg_gbxar = rdg_gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2
 
     call infld('HWDTH', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_hwdth, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: HWDTH not found on topo file')

     call infld('CLNGT', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_clngt, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: CLNGT not found on topo file')

     call infld('MXDIS', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_mxdis, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: MXDIS not found on topo file')

     call infld('ANIXY', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_anixy, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANIXY not found on topo file')

     call infld('ANGLL', fh_topo, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_angll, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANGLL not found on topo file')

     ! close topo file only if it was opened here
     if (len_trim(bnd_topo_loc) > 0) then
        call pio_closefile(fh_topo)
     end if

     call addfld('UEGW',  (/ 'lev' /) , 'A'  ,'1/s' ,  &
          'Zonal wind profile-entry to GW ' )
     call addfld('VEGW',  (/ 'lev' /) , 'A'  ,'1/s' ,  &
          'Merdional wind profile-entry to GW ' )
     call register_vector_field('UEGW','VEGW')
     call addfld('TEGW',  (/ 'lev' /) , 'A'  ,'K' ,  &
          'Temperature profile-entry to GW ' )

     call addfld('TAU1RDGBETAM' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call addfld('UBM1BETA',  (/ 'lev' /) , 'A'  ,'1/s' ,  &
          'On-ridge wind profile           ' )
     call addfld('UBT1RDGBETA' , (/ 'lev' /) , 'I'  ,'m/s' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGBETAY' , (/ 'ilev' /), 'I', 'N/m2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGBETAX' , (/ 'ilev' /), 'I', 'N/m2', &
          'Ridge based momentum flux profile')
        call register_vector_field('TAU'//cn//'RDGBETAX','TAU'//cn//'RDGBETAY')
        call addfld('UT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm/s', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm/s', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGBETA','VT'//cn//'RDGBETA')
     end do

     call addfld('TAUARDGBETAY' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call addfld('TAUARDGBETAX' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGBETAX','TAUARDGBETAY')

     if (history_waccm) then
        call add_default('TAUARDGBETAX', 1, ' ')
        call add_default('TAUARDGBETAY  ', 1, ' ')
     end if

  end if

  if (use_gw_rdg_gamma) then

     if (effgw_rdg_gamma == unset_r8) then
        call endrun(sub//": ERROR: Anisotropic OGW enabled, but effgw_rdg_gamma was not set.")
     end if
     
     call getfil(bnd_rdggm, bnd_rdggm_loc, iflag=1, lexist=found)
     if (found) then
        call cam_pio_openfile(fh_rdggm, bnd_rdggm_loc, PIO_NOWRITE)
     else
        call endrun(sub//': ERROR: file for gamma ridges not found: bnd_rdggm='// &
                   trim(bnd_rdggm))
     end if

     if (.not. allocated(rdg_gbxar)) then
        allocate(rdg_gbxar(pcols,begchunk:endchunk))
        call infld('GBXAR', fh_rdggm, dim1name, dim2name, 1, pcols, &
                            begchunk, endchunk, rdg_gbxar, found, gridname='physgrid')
        if (.not. found) call endrun(sub//': ERROR: GBXAR not found on bnd_rdggm')
        rdg_gbxar = rdg_gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2
     end if

     ! Get meso-gamma ridge data
     allocate( &
        rdg_hwdthg(pcols,prdg,begchunk:endchunk), &
        rdg_clngtg(pcols,prdg,begchunk:endchunk), &
        rdg_mxdisg(pcols,prdg,begchunk:endchunk), &
        rdg_anixyg(pcols,prdg,begchunk:endchunk), &
        rdg_angllg(pcols,prdg,begchunk:endchunk)  )

     call infld('HWDTH', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_hwdthg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: HWDTH not found on bnd_rdggm')

     call infld('CLNGT', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_clngtg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: CLNGT not found on bnd_rdggm')

     call infld('MXDIS', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_mxdisg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: MXDIS not found on bnd_rdggm')

     call infld('ANIXY', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_anixyg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANIXY not found on bnd_rdggm')

     call infld('ANGLL', fh_rdggm, dim1name, 'nrdg', dim2name, 1, pcols, &
                1, prdg, begchunk, endchunk, rdg_angllg, found, gridname='physgrid')
     if (.not. found) call endrun(sub//': ERROR: ANGLL not found on bnd_rdggm')
 
     call pio_closefile(fh_rdggm)

     call addfld ('TAU1RDGGAMMAM' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call addfld ('UBM1GAMMA',  (/ 'lev' /) , 'A'  ,'1/s' ,  &
          'On-ridge wind profile           ' )
     call addfld ('UBT1RDGGAMMA' , (/ 'lev' /) , 'I'  ,'m/s' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGGAMMAY', (/ 'ilev' /), 'I', 'N/m2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGGAMMAX', (/ 'ilev' /), 'I', 'N/m2', &
          'Ridge based momentum flux profile')
        call addfld('UT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm/s', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm/s', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGGAMMA','VT'//cn//'RDGGAMMA')
     end do

     call addfld ('TAUARDGGAMMAY' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call addfld ('TAUARDGGAMMAX' , (/ 'ilev' /) , 'I'  ,'N/m2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGGAMMAX','TAUARDGGAMMAY')
     call addfld ('TAURDGGMX',     horiz_only,  'A','N/m2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAURDGGMY',     horiz_only,  'A','N/m2', &
          'Meridional gravity wave surface stress')
     call register_vector_field('TAURDGGMX','TAURDGGMY')
     call addfld ('UTRDGGM' , (/ 'lev' /) , 'I'  ,'m/s' , &
          'U wind tendency from ridge 6     ')
     call addfld ('VTRDGGM' , (/ 'lev' /) , 'I'  ,'m/s' , &
          'V wind tendency from ridge 6     ')
     call register_vector_field('UTRDGGM','VTRDGGM')
  end if

  if (use_gw_front .or. use_gw_front_igw) then

     frontgf_idx = pbuf_get_index('FRONTGF')
     frontga_idx = pbuf_get_index('FRONTGA')

     call shr_assert(unset_r8 /= frontgfc, &
          "gw_drag_init: Frontogenesis enabled, but frontgfc was &
          & not set!"// &
          errMsg(__FILE__, __LINE__))

     do k = 0, pver
        ! Check frontogenesis at 600 hPa.
        if (pref_edge(k+1) < 60000._r8) kfront = k+1
     end do

     ! Source waves from 500 hPa.
     kbot_front = maxloc(pref_edge, 1, (pref_edge < 50000._r8)) - 1

     if (masterproc) then
        write (iulog,*) 'KFRONT      =',kfront
        write (iulog,*) 'KBOT_FRONT  =',kbot_front
        write(iulog,*) ' '
     end if

     call addfld ('FRONTGF', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')
     call addfld ('FRONTGFA', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')

     if (history_waccm) then
        call add_default('FRONTGF', 1, ' ')
        call add_default('FRONTGFA', 1, ' ')
     end if

  end if

  if (use_gw_front) then

     call shr_assert(all(unset_r8 /= [ effgw_cm, taubgnd ]), &
          "gw_drag_init: Frontogenesis mid-scale waves enabled, but not &
          &all required namelist variables were set!"// &
          errMsg(__FILE__, __LINE__))

     if (masterproc) then
        write(iulog,*) 'gw_init: gw spectrum taubgnd, ', &
             'effgw_cm = ',taubgnd, effgw_cm
        write(iulog,*) ' '
     end if

     cm_desc = gaussian_cm_desc(band_mid, kbot_front, kfront, frontgfc, &
          taubgnd, front_gaussian_width)

     ! Output for gravity waves from frontogenesis.
     call gw_spec_addflds(prefix=cm_pf, scheme="C&M", band=band_mid, &
          history_defaults=history_waccm)

  end if

  if (use_gw_front_igw) then

     call shr_assert(all(unset_r8 /= [ effgw_cm_igw, taubgnd_igw ]), &
          "gw_drag_init: Frontogenesis inertial waves enabled, but not &
          &all required namelist variables were set!"// &
          errMsg(__FILE__, __LINE__))

     if (masterproc) then
        write(iulog,*) 'gw_init: gw spectrum taubgnd_igw, ', &
             'effgw_cm_igw = ',taubgnd_igw, effgw_cm_igw
        write(iulog,*) ' '
     end if

     cm_igw_desc = gaussian_cm_desc(band_long, kbot_front, kfront, frontgfc, &
          taubgnd_igw, front_gaussian_width)

     ! Output for gravity waves from frontogenesis.
     call gw_spec_addflds(prefix=cm_igw_pf, scheme="C&M IGW", &
          band=band_long, history_defaults=history_waccm)

  end if

  if (use_gw_convect_dp) then

     ttend_dp_idx    = pbuf_get_index('TTEND_DP')

     ! Set the deep scheme specification components.
     beres_dp_desc%storm_shift = .true.

     do k = 0, pver
        ! 700 hPa index
        if (pref_edge(k+1) < 70000._r8) beres_dp_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'Beres deep level =',beres_dp_desc%k
     end if

     ! Don't use deep convection heating depths below this limit.
     ! This is important for QBO. Bad result if left at 2.5 km.
     beres_dp_desc%min_hdepth = 1000._r8

     ! Read Beres file.

     call shr_assert(trim(gw_drag_file) /= "", &
          "gw_drag_init: No gw_drag_file provided for Beres deep &
          &scheme. Set this via namelist."// &
          errMsg(__FILE__, __LINE__))

     call gw_init_beres(gw_drag_file, band_mid, beres_dp_desc)

     ! Output for gravity waves from the Beres scheme (deep).
     call gw_spec_addflds(prefix=beres_dp_pf, scheme="Beres (deep)", &
          band=band_mid, history_defaults=history_waccm)

     call addfld ('NETDT',(/ 'lev' /), 'A','K/s', &
          'Net heating rate')
     call addfld ('MAXQ0',horiz_only  ,  'A','K/day', &
          'Max column heating rate')
     call addfld ('HDEPTH',horiz_only,    'A','km', &
          'Heating Depth')

     if (history_waccm) then
        call add_default('NETDT    ', 1, ' ')
        call add_default('HDEPTH   ', 1, ' ')
        call add_default('MAXQ0    ', 1, ' ')
     end if

  end if

  if (use_gw_convect_sh) then

     ttend_sh_idx    = pbuf_get_index('TTEND_SH')

     ! Set the shallow scheme specification components.
     beres_sh_desc%storm_shift = .false.

     do k = 0, pver
        ! 900 hPa index
        if (pref_edge(k+1) < 90000._r8) beres_sh_desc%k = k+1
     end do

     if (masterproc) then
        write (iulog,*) 'Beres shallow level =',beres_sh_desc%k
     end if

     ! Use all heating depths for shallow convection.
     beres_sh_desc%min_hdepth = 0._r8

     ! Read Beres file.

     call shr_assert(trim(gw_drag_file_sh) /= "", &
          "gw_drag_init: No gw_drag_file_sh provided for Beres shallow &
          &scheme. Set this via namelist."// &
          errMsg(__FILE__, __LINE__))

     call gw_init_beres(gw_drag_file_sh, band_mid, beres_sh_desc)

     ! Output for gravity waves from the Beres scheme (shallow).
     call gw_spec_addflds(prefix=beres_sh_pf, scheme="Beres (shallow)", &
          band=band_mid, history_defaults=history_waccm)

     call addfld ('SNETDT',(/ 'lev' /), 'A','K/s', &
          'Net heating rate')
     call addfld ('SMAXQ0',horiz_only  ,  'A','K/day', &
          'Max column heating rate')
     call addfld ('SHDEPTH',horiz_only,    'A','km', &
          'Heating Depth')

     if (history_waccm) then
        call add_default('SNETDT   ', 1, ' ')
        call add_default('SHDEPTH  ', 1, ' ')
        call add_default('SMAXQ0   ', 1, ' ')
     end if

  end if

  call addfld ('EKGW' ,(/ 'ilev' /), 'A','M2/S', &
       'Effective Kzz due to diffusion by gravity waves')

  if (history_waccm) then
     call add_default('EKGW', 1, ' ')
  end if

  call addfld ('UTGW_TOTAL',    (/ 'lev' /), 'A','m/s2', &
       'Total U tendency due to gravity wave drag')
  call addfld ('VTGW_TOTAL',    (/ 'lev' /), 'A','m/s2', &
       'Total V tendency due to gravity wave drag')
  call register_vector_field('UTGW_TOTAL', 'VTGW_TOTAL')

  ! Total temperature tendency output.
  call addfld ('TTGW', (/ 'lev' /), 'A', 'K/s',  &
       'T tendency - gravity wave drag')

  ! Water budget terms.
  call addfld('QTGW',(/ 'lev' /), 'A','kg/kg/s', &
       'Q tendency - gravity wave drag')
  call addfld('CLDLIQTGW',(/ 'lev' /), 'A','kg/kg/s', &
       'CLDLIQ tendency - gravity wave drag')
  call addfld('CLDICETGW',(/ 'lev' /), 'A','kg/kg/s', &
       'CLDICE tendency - gravity wave drag')

  if ( history_budget ) then
     call add_default('TTGW', history_budget_histfile_num, ' ')
     call add_default('QTGW', history_budget_histfile_num, ' ')
     call add_default('CLDLIQTGW', history_budget_histfile_num, ' ')
     call add_default('CLDICETGW', history_budget_histfile_num, ' ')
  end if

  ! Get indices to actually output the above.
  call cnst_get_ind("CLDLIQ", ixcldliq)
  call cnst_get_ind("CLDICE", ixcldice)

  if (gw_top_taper) then
     allocate(vramp(pver))
     vramp(:) = 1._r8
     topndx = 1
     botndx = press_lim_idx( 0.6E-02_r8, top=.true. )
     if (botndx>1) then
        do k=botndx,topndx,-1
           vramp(k) = vramp(k+1)/(pref_edge(k+1)/pref_edge(k))
        end do
        if (masterproc) then
           write(iulog,'(A)') 'GW taper coef (vramp):'
           do k=1,pver
              write(iulog,"('k: ',I4,' taper coef,press(Pa): ',F12.8,E12.4)") k, vramp(k), pref_mid(k)
           enddo
        endif
     endif
  end if

end subroutine gw_init

!==========================================================================

subroutine gw_init_beres(file_name, band, desc)

  use ioFileMod, only: getfil
  use pio, only: file_desc_t, pio_nowrite, pio_inq_varid, pio_get_var, &
       pio_closefile
  use cam_pio_utils, only: cam_pio_openfile

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band

  type(BeresSourceDesc), intent(inout) :: desc

  type(file_desc_t) :: gw_file_desc

  ! PIO variable ids and error code.
  integer :: mfccid, hdid, stat

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  ! Full path to gw_drag_file.
  character(len=256) :: file_path

  character(len=256) :: msg

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  call getfil(file_name, file_path)
  call cam_pio_openfile(gw_file_desc, file_path, pio_nowrite)

  ! Get HD (heating depth) dimension.

  desc%maxh = get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! Get MW (mean wind) dimension.

  desc%maxuh = get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.

  ngwv_file = get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2

  call shr_assert(ngwv_file >= band%ngwv, &
       "gw_beres_init: PS in lookup table file does not cover the whole &
       &spectrum implied by the model's ngwv.")

  ! Allocate hd and get data.

  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_init_beres: Allocation error (hd): "//msg// &
       errMsg(__FILE__, __LINE__))

  stat = pio_inq_varid(gw_file_desc,'HD',hdid)

  call handle_pio_error(stat, &
       'Error finding HD in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, hdid, start=[1], count=[desc%maxh], &
       ival=desc%hd)

  call handle_pio_error(stat, &
       'Error reading HD from: '//trim(file_path))

  ! While not currently documented in the file, it uses kilometers. Convert
  ! to meters.
  desc%hd = desc%hd*1000._r8

  ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
  ! model determines wavenumber dimension.

  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)

  call shr_assert(stat == 0, &
       "gw_init_beres: Allocation error (mfcc): "//msg// &
       errMsg(__FILE__, __LINE__))

  ! Get mfcc data.

  stat = pio_inq_varid(gw_file_desc,'mfcc',mfccid)

  call handle_pio_error(stat, &
       'Error finding mfcc in: '//trim(file_path))

  stat = pio_get_var(gw_file_desc, mfccid, &
       start=[1,1,ngwv_file-band%ngwv+1], count=shape(desc%mfcc), &
       ival=desc%mfcc)

  call handle_pio_error(stat, &
       'Error reading mfcc from: '//trim(file_path))

  call pio_closefile(gw_file_desc)

  if (masterproc) then

     write(iulog,*) "Read in source spectra from file."
     write(iulog,*) "mfcc max, min = ", &
          maxval(desc%mfcc),", ",minval(desc%mfcc)

  endif

end subroutine gw_init_beres

!==========================================================================

! Utility to reduce the repetitiveness of reads during initialization.
function get_pio_dimlen(file_desc, dim_name, file_path) result(dimlen)

  use pio, only: file_desc_t, pio_inq_dimid, pio_inq_dimlen

  type(file_desc_t), intent(in) :: file_desc
  character(len=*), intent(in) :: dim_name

  ! File path, for use in error messages only.
  character(len=*), intent(in) :: file_path

  integer :: dimlen

  integer :: dimid, stat

  stat = pio_inq_dimid(file_desc, dim_name, dimid)

  call handle_pio_error(stat, &
       "Error finding dimension "//dim_name//" in: "//file_path)

  stat = pio_inq_dimlen(file_desc, dimid, dimlen)

  call handle_pio_error(stat, &
       "Error reading dimension "//dim_name//" from: "//file_path)

end function get_pio_dimlen

!==========================================================================

! In fact, we'd usually expect PIO errors to abort the run before you can
! even check the error code. But just in case, use this little assert.
subroutine handle_pio_error(stat, message)
  use pio, only: pio_noerr
  integer, intent(in) :: stat
  character(len=*) :: message

  call shr_assert(stat == pio_noerr, &
       "PIO error in gw_init_beres: "//trim(message)// &
       errMsg(__FILE__, __LINE__))

end subroutine handle_pio_error

!==========================================================================

subroutine gw_tend(state, pbuf, dt, ptend, cam_in, flx_heat)
  !-----------------------------------------------------------------------
  ! Interface for multiple gravity wave drag parameterization.
  !-----------------------------------------------------------------------

  use physics_types,  only: physics_state_copy, set_dry_to_wet
  use constituents,   only: cnst_type
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field
  use camsrfexch, only: cam_in_t
  ! Location-dependent cpair
  use physconst,  only: cpairv, pi
  use coords_1d,  only: Coords1D
  use gw_common,  only: gw_prof, gw_drag_prof, calc_taucd, &
       momentum_flux, momentum_fixer, energy_change, energy_fixer, &
       coriolis_speed, adjust_inertial
  use gw_oro,     only: gw_oro_src
  use gw_front,   only: gw_cm_src
  use gw_convect, only: gw_beres_src

  !------------------------------Arguments--------------------------------
  type(physics_state), intent(in) :: state   ! physics state structure
  type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer
  real(r8), intent(in) :: dt                    ! time step
  ! Parameterization net tendencies.
  type(physics_ptend), intent(out):: ptend
  type(cam_in_t), intent(in) :: cam_in
  real(r8), intent(out) :: flx_heat(pcols)

  !---------------------------Local storage-------------------------------

  type(physics_state) :: state1     ! Local copy of state variable

  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns

  integer :: i, k                   ! loop indices

  type(Coords1D) :: p               ! Pressure coordinates

  real(r8) :: ttgw(state%ncol,pver) ! temperature tendency
  real(r8) :: utgw(state%ncol,pver) ! zonal wind tendency
  real(r8) :: vtgw(state%ncol,pver) ! meridional wind tendency

  real(r8) :: ni(state%ncol,pver+1) ! interface Brunt-Vaisalla frequency
  real(r8) :: nm(state%ncol,pver)   ! midpoint Brunt-Vaisalla frequency
  real(r8) :: rhoi(state%ncol,pver+1)     ! interface density
  real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
  real(r8) :: tau0x(state%ncol)     ! c=0 sfc. stress (zonal)
  real(r8) :: tau0y(state%ncol)     ! c=0 sfc. stress (meridional)
  real(r8) :: ubi(state%ncol,pver+1)! projection of wind at interfaces
  real(r8) :: ubm(state%ncol,pver)  ! projection of wind at midpoints
  real(r8) :: xv(state%ncol)        ! unit vector of source wind (x)
  real(r8) :: yv(state%ncol)        ! unit vector of source wind (y)

  integer :: m                      ! dummy integers
  real(r8) :: qtgw(state%ncol,pver,pcnst) ! constituents tendencies

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8) :: taucd(state%ncol,pver+1,4)

  ! gravity wave wind tendency for each wave
  real(r8), allocatable :: gwut(:,:,:)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(state%ncol,pver)
  real(r8) :: dttke(state%ncol,pver)

  ! Wave phase speeds for each column
  real(r8), allocatable :: c(:,:)

  ! Efficiency for a gravity wave source.
  real(r8) :: effgw(state%ncol)

  ! Coriolis characteristic speed.
  real(r8) :: u_coriolis(state%ncol)

  ! Adjustment for inertial gravity waves.
  real(r8), allocatable :: ro_adjust(:,:,:)

  ! pbuf fields
  ! Molecular diffusivity
  real(r8), pointer :: kvt_in(:,:)
  real(r8) :: kvtt(state%ncol,pver+1)

  ! Frontogenesis
  real(r8), pointer :: frontgf(:,:)
  real(r8), pointer :: frontga(:,:)

  ! Temperature change due to deep convection.
  real(r8), pointer :: ttend_dp(:,:)
  ! Temperature change due to shallow convection.
  real(r8), pointer :: ttend_sh(:,:)

  ! Standard deviation of orography.
  real(r8), pointer :: sgh(:)

  ! gridbox area
  real(r8), pointer :: gbxar(:)

     ! Beta ridges
  ! width of ridges.
  real(r8), pointer :: hwdth(:,:)
  ! length of ridges.
  real(r8), pointer :: clngt(:,:)
  ! Maximum deviations of ridges.
  real(r8), pointer :: mxdis(:,:)
  ! orientation of ridges.
  real(r8), pointer :: angll(:,:)
  ! anisotropy of ridges.
  real(r8), pointer :: anixy(:,:)

     ! Gamma ridges
  ! width of ridges.
  real(r8), pointer :: hwdthg(:,:)
  ! length of ridges.
  real(r8), pointer :: clngtg(:,:)
  ! Maximum deviations of ridges.
  real(r8), pointer :: mxdisg(:,:)
  ! orientation of ridges.
  real(r8), pointer :: angllg(:,:)
  ! anisotropy of ridges.
  real(r8), pointer :: anixyg(:,:)

  ! Indices of gravity wave source and lowest level where wind tendencies
  ! are allowed.
  integer :: src_level(state%ncol)
  integer :: tend_level(state%ncol)

  ! Convective source heating depth.
  ! heating depth
  real(r8) :: hdepth(state%ncol)
  ! maximum heating rate
  real(r8) :: maxq0(state%ncol)

  ! Scale sgh to account for landfrac.
  real(r8) :: sgh_scaled(state%ncol)

  ! Parameters for the IGW polar taper.
  real(r8), parameter :: degree2radian = pi/180._r8
  real(r8), parameter :: al0 = 82.5_r8 * degree2radian
  real(r8), parameter :: dlat0 = 5.0_r8 * degree2radian
  
  ! effective gw diffusivity at interfaces needed for output
  real(r8) :: egwdffi(state%ncol,pver+1)
  ! sum from the two types of spectral GW
  real(r8) :: egwdffi_tot(state%ncol,pver+1)

  ! Momentum fluxes used by fixer.
  real(r8) :: um_flux(state%ncol), vm_flux(state%ncol)
  ! Energy change used by fixer.
  real(r8) :: de(state%ncol)

  ! Which constituents are being affected by diffusion.
  logical  :: lq(pcnst)

  ! Contiguous copies of state arrays.
  real(r8) :: dse(state%ncol,pver)
  real(r8) :: t(state%ncol,pver)
  real(r8) :: u(state%ncol,pver)
  real(r8) :: v(state%ncol,pver)
  real(r8) :: q(state%ncol,pver,pcnst)
  real(r8) :: piln(state%ncol,pver+1)
  real(r8) :: zm(state%ncol,pver)
  real(r8) :: zi(state%ncol,pver+1)
  !------------------------------------------------------------------------

  ! Make local copy of input state.
  call physics_state_copy(state, state1)

  ! constituents are all treated as wet mmr
  call set_dry_to_wet(state1)

  lchnk = state1%lchnk
  ncol  = state1%ncol

  p = Coords1D(state1%pint(:ncol,:))

  dse = state1%s(:ncol,:)
  t = state1%t(:ncol,:)
  u = state1%u(:ncol,:)
  v = state1%v(:ncol,:)
  q = state1%q(:ncol,:,:)
  piln = state1%lnpint(:ncol,:)
  zm = state1%zm(:ncol,:)
  zi = state1%zi(:ncol,:)

  lq = .true.
  call physics_ptend_init(ptend, state1%psetcols, "Gravity wave drag", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  ! Profiles of background state variables
  call gw_prof(ncol, p, cpair, t, rhoi, nm, ni)

  if (do_molec_diff) then
     !--------------------------------------------------------
     ! Initialize and calculate local molecular diffusivity
     !--------------------------------------------------------

     call pbuf_get_field(pbuf, kvt_idx, kvt_in)  ! kvt_in(1:pcols,1:pver+1)

     ! Set kvtt from pbuf field; kvtt still needs a factor of 1/cpairv.
     kvtt = kvt_in(:ncol,:)

     ! Use linear extrapolation of cpairv to top interface.
     kvtt(:,1) = kvtt(:,1) / &
          (1.5_r8*cpairv(:ncol,1,lchnk) - &
          0.5_r8*cpairv(:ncol,2,lchnk))

     ! Interpolate cpairv to other interfaces.
     do k = 2, nbot_molec
        kvtt(:,k) = kvtt(:,k) / &
             (cpairv(:ncol,k+1,lchnk)+cpairv(:ncol,k,lchnk)) * 2._r8
     enddo

  else

     kvtt = 0._r8

  end if

  if (use_gw_front_igw) then
     u_coriolis = coriolis_speed(band_long, state1%lat(:ncol))
  end if

  ! Totals that accumulate over different sources.
  egwdffi_tot = 0._r8
  flx_heat = 0._r8
  
  if (use_gw_convect_dp) then
     !------------------------------------------------------------------
     ! Convective gravity waves (Beres scheme, deep).
     !------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1))
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv))
     allocate(c(ncol,-band_mid%ngwv:band_mid%ngwv))

     ! Set up heating
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._r8 - abs(state1%lat(:ncol)) >= 4*epsilon(1._r8))
        effgw = effgw_beres_dp
     elsewhere
        effgw = 0._r8
     end where

     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, band_mid, beres_dp_desc, &
          u, v, ttend_dp(:ncol,:), zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0)

     ! Solve for the drag profile with Beres source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
          t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   c,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, c, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Store constituents tendencies
     do m=1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! Add the momentum tendencies to the output tendency arrays.
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

     ! Change ttgw to a temperature tendency before outputing it.
     ttgw = ttgw / cpair
     call gw_spec_outflds(beres_dp_pf, lchnk, ncol, band_mid, c, u, v, &
          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
          taucd)

     ! Diagnostic outputs (convert hdepth to km).
     call outfld('NETDT', ttend_dp, pcols, lchnk)
     call outfld('HDEPTH', hdepth/1000._r8, ncol, lchnk)
     call outfld('MAXQ0', maxq0, ncol, lchnk)

     deallocate(tau, gwut, c)

  end if

  if (use_gw_convect_sh) then
     !------------------------------------------------------------------
     ! Convective gravity waves (Beres scheme, shallow).
     !------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1))
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv))
     allocate(c(ncol,-band_mid%ngwv:band_mid%ngwv))

     ! Set up heating
     call pbuf_get_field(pbuf, ttend_sh_idx, ttend_sh)

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2._r8 - abs(state1%lat(:ncol)) >= 4*epsilon(1._r8))
        effgw = effgw_beres_sh
     elsewhere
        effgw = 0._r8
     end where

     ! Determine wave sources for Beres shallow scheme
     call gw_beres_src(ncol, band_mid, beres_sh_desc, &
          u, v, ttend_sh(:ncol,:), zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0)

     ! Solve for the drag profile with Beres source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level,  dt, &
          t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   c,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, c, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Store constituents tendencies
     do m=1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Add the momentum tendencies to the output tendency arrays.
     ! Don't calculate fixers, since we are too close to the ground to
     ! spread momentum/energy differences across low layers.
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

     ! Calculate energy change for output to CAM's energy checker.
     ! This is sort of cheating; we don't have a good a priori idea of the
     ! energy coming from surface stress, so we just integrate what we and
     ! actually have so far and overwrite flx_heat with that.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
     flx_heat(:ncol) = de

     ! Change ttgw to a temperature tendency before outputing it.
     ttgw = ttgw / cpair
     call gw_spec_outflds(beres_sh_pf, lchnk, ncol, band_mid, c, u, v, &
          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
          taucd)

     ! Diagnostic outputs (convert SHDEPTH to km).
     call outfld ('SNETDT', ttend_sh, pcols, lchnk)
     call outfld ('SHDEPTH', hdepth/1000._r8, ncol, lchnk)
     call outfld ('SMAXQ0', maxq0, ncol, lchnk)

     deallocate(tau, gwut, c)

  end if

  if (use_gw_front .or. use_gw_front_igw) then
     ! Get frontogenesis physics buffer fields set by dynamics.
     call pbuf_get_field(pbuf, frontgf_idx, frontgf)
     call pbuf_get_field(pbuf, frontga_idx, frontga)

     ! Output for diagnostics.
     call outfld ('FRONTGF', frontgf, pcols, lchnk)
     call outfld ('FRONTGFA', frontga, pcols, lchnk)
  end if

  if (use_gw_front) then
     !------------------------------------------------------------------
     ! Frontally generated gravity waves
     !------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1))
     allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv))
     allocate(c(ncol,-band_mid%ngwv:band_mid%ngwv))

     ! Efficiency of gravity wave momentum transfer.
     effgw = effgw_cm
     ! Frontogenesis is too high at the poles (at least for the FV
     ! dycore), so introduce a polar taper.
     if (gw_polar_taper) effgw = effgw * cos(state1%lat(:ncol))

     ! Determine the wave source for C&M background spectrum
     call gw_cm_src(ncol, band_mid, cm_desc, u, v, frontgf(:ncol,:), &
          src_level, tend_level, tau, ubm, ubi, xv, yv, c)

     ! Solve for the drag profile with C&M source spectrum.
     call gw_drag_prof(ncol, band_mid, p, src_level, tend_level,  dt, &
          t, vramp,   &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   c,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, c, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     !Add the constituent tendencies
     do m=1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! add the momentum tendencies to the output tendency arrays
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

     ! Change ttgw to a temperature tendency before outputing it.
     ttgw = ttgw / cpair
     call gw_spec_outflds(cm_pf, lchnk, ncol, band_mid, c, u, v, &
          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
          taucd)

     deallocate(tau, gwut, c)

  end if

  if (use_gw_front_igw) then
     !------------------------------------------------------------------
     ! Frontally generated inertial gravity waves
     !------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,-band_long%ngwv:band_long%ngwv,pver+1))
     allocate(gwut(ncol,pver,-band_long%ngwv:band_long%ngwv))
     allocate(c(ncol,-band_long%ngwv:band_long%ngwv))
     allocate(ro_adjust(ncol,-band_long%ngwv:band_long%ngwv,pver+1))

     ! Efficiency of gravity wave momentum transfer.
     effgw = effgw_cm_igw

     ! Frontogenesis is too high at the poles (at least for the FV
     ! dycore), so introduce a polar taper.
     if (gw_polar_taper) then
        where (abs(state1%lat(:ncol)) <= 89._r8*degree2radian)
           effgw = effgw * 0.25_r8 * &
                 (1._r8+tanh((state1%lat(:ncol)+al0)/dlat0)) * &
                 (1._r8-tanh((state1%lat(:ncol)-al0)/dlat0))
        elsewhere
           effgw = 0._r8
        end where
     end if

     ! Determine the wave source for C&M background spectrum
     call gw_cm_src(ncol, band_long, cm_igw_desc, u, v, frontgf(:ncol,:), &
          src_level, tend_level, tau, ubm, ubi, xv, yv, c)

     call adjust_inertial(band_long, tend_level, u_coriolis, c, ubi, &
          tau, ro_adjust)

     ! Solve for the drag profile with C&M source spectrum.
     call gw_drag_prof(ncol, band_long, p, src_level, tend_level,  dt, &
          t, vramp,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   c,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke, ro_adjust=ro_adjust, &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_long%ngwv, tend_level, tau, c, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     !Add the constituent tendencies
     do m=1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! add the momentum tendencies to the output tendency arrays
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
     call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

     ! Change ttgw to a temperature tendency before outputing it.
     ttgw = ttgw / cpair
     call gw_spec_outflds(cm_igw_pf, lchnk, ncol, band_long, c, u, v, &
          xv, yv, gwut, dttdf, dttke, tau(:,:,2:), utgw, vtgw, ttgw, &
          taucd)

     deallocate(tau, gwut, c, ro_adjust)

  end if

  if (use_gw_oro) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Allocate wavenumber fields.
     allocate(tau(ncol,band_oro%ngwv:band_oro%ngwv,pver+1))
     allocate(gwut(ncol,pver,band_oro%ngwv:band_oro%ngwv))
     allocate(c(ncol,band_oro%ngwv:band_oro%ngwv))

     ! Efficiency of gravity wave momentum transfer.
     ! Take into account that wave sources are only over land.
     call pbuf_get_field(pbuf, sgh_idx, sgh)
 
     if (gw_lndscl_sgh) then
        where (cam_in%landfrac(:ncol) >= epsilon(1._r8))
           effgw = effgw_oro * cam_in%landfrac(:ncol)
           sgh_scaled = sgh(:ncol) / sqrt(cam_in%landfrac(:ncol))
        elsewhere
           effgw = 0._r8
           sgh_scaled = 0._r8
        end where

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band_oro, p, &
           u, v, t, sgh_scaled, zm, nm, &
           src_level, tend_level, tau, ubm, ubi, xv, yv, c)
     else
        effgw = effgw_oro

        ! Determine the orographic wave source
        call gw_oro_src(ncol, band_oro, p, &
             u, v, t, sgh(:ncol), zm, nm, &
             src_level, tend_level, tau, ubm, ubi, xv, yv, c)
     endif
     do i = 1, ncol
        if (state1%lat(i) < 0._r8) then
           tau(i,:,:) = tau(i,:,:) * gw_oro_south_fac
        end if
     end do

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, band_oro, p, src_level, tend_level,   dt,   &
          t, vramp,   &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,c,          kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! For orographic waves, don't bother with taucd, since there are no
     ! momentum conservation routines or directional diagnostics.

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     ! Add the orographic tendencies to the spectrum tendencies.
     ! Don't calculate fixers, since we are too close to the ground to
     ! spread momentum/energy differences across low layers.
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
        ! Convert to temperature tendency for output.
        ttgw(:,k) = ttgw(:,k) / cpairv(:ncol, k, lchnk)
     end do

     ! Calculate energy change for output to CAM's energy checker.
     ! This is sort of cheating; we don't have a good a priori idea of the
     ! energy coming from surface stress, so we just integrate what we and
     ! actually have so far and overwrite flx_heat with that.
     call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
     flx_heat(:ncol) = de

     do m = 1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Write output fields to history file
     call outfld('TAUAORO', tau(:,0,:),  ncol, lchnk)
     call outfld('UTGWORO', utgw,  ncol, lchnk)
     call outfld('VTGWORO', vtgw,  ncol, lchnk)
     call outfld('TTGWORO', ttgw,  ncol, lchnk)
     call outfld('TTGWSDFORO', dttdf / cpair,  ncol, lchnk)
     call outfld('TTGWSKEORO', dttke / cpair,  ncol, lchnk)
     tau0x = tau(:,0,pver+1) * xv
     tau0y = tau(:,0,pver+1) * yv
     call outfld('TAUGWX', tau0x, ncol, lchnk)
     call outfld('TAUGWY', tau0y, ncol, lchnk)

     deallocate(tau, gwut, c)

  end if

  if (use_gw_rdg_beta) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Efficiency of gravity wave momentum transfer.
     ! Take into account that wave sources are only over land.

     hwdth => rdg_hwdth(:ncol,:,lchnk)
     clngt => rdg_clngt(:ncol,:,lchnk)
     gbxar => rdg_gbxar(:ncol,lchnk)
     mxdis => rdg_mxdis(:ncol,:,lchnk)
     angll => rdg_angll(:ncol,:,lchnk)
     anixy => rdg_anixy(:ncol,:,lchnk)

     where(mxdis < 0._r8)
        mxdis = 0._r8
     end where

     ! Save state at top of routine
     ! Useful for unit testing checks
     call outfld('UEGW', u ,  ncol, lchnk)
     call outfld('VEGW', v ,  ncol, lchnk)
     call outfld('TEGW', t ,  ncol, lchnk)

     call gw_rdg_calc(&
        'BETA ', ncol, lchnk, n_rdg_beta, dt,     &
        u, v, t, p, piln, zm, zi,                 &
        nm, ni, rhoi, kvtt, q, dse,               &
        effgw_rdg_beta, effgw_rdg_beta_max,       &
        hwdth, clngt, gbxar, mxdis, angll, anixy, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta,     &
        ptend, flx_heat)

  end if

  if (use_gw_rdg_gamma) then
     !---------------------------------------------------------------------
     ! Orographic stationary gravity waves
     !---------------------------------------------------------------------

     ! Efficiency of gravity wave momentum transfer.
     ! Take into account that wave sources are only over land.

     hwdthg => rdg_hwdthg(:ncol,:,lchnk)
     clngtg => rdg_clngtg(:ncol,:,lchnk)
     gbxar  => rdg_gbxar(:ncol,lchnk)
     mxdisg => rdg_mxdisg(:ncol,:,lchnk)
     angllg => rdg_angllg(:ncol,:,lchnk)
     anixyg => rdg_anixyg(:ncol,:,lchnk)

     where(mxdisg < 0._r8)
        mxdisg = 0._r8
     end where

     call gw_rdg_calc(&
        'GAMMA', ncol, lchnk, n_rdg_gamma, dt,         &
        u, v, t, p, piln, zm, zi,                      &
        nm, ni, rhoi, kvtt, q, dse,                    &
        effgw_rdg_gamma, effgw_rdg_gamma_max,          &
        hwdthg, clngtg, gbxar, mxdisg, angllg, anixyg, &
        rdg_gamma_cd_llb, trpd_leewv_rdg_gamma,        &
        ptend, flx_heat)

  endif

  ! Convert the tendencies for the dry constituents to dry air basis.
  do m = 1, pcnst
     if (cnst_type(m).eq.'dry') then
        do k = 1, pver
           do i = 1, ncol
              ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
           end do
        end do
     end if
  end do

  ! Write totals to history file.
  call outfld('EKGW', egwdffi_tot , ncol, lchnk)
  call outfld('TTGW', ptend%s/cpairv(:,:,lchnk),  pcols, lchnk)
 
  call outfld('UTGW_TOTAL', ptend%u, pcols, lchnk)
  call outfld('VTGW_TOTAL', ptend%v, pcols, lchnk)

  call outfld('QTGW', ptend%q(:,:,1), pcols, lchnk)
  call outfld('CLDLIQTGW', ptend%q(:,:,ixcldliq), pcols, lchnk)
  call outfld('CLDICETGW', ptend%q(:,:,ixcldice), pcols, lchnk)

  ! Destroy objects.
  call p%finalize()

end subroutine gw_tend

!==========================================================================

subroutine gw_rdg_calc( &
   type, ncol, lchnk, n_rdg, dt, &
   u, v, t, p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, &
   effgw_rdg, effgw_rdg_max, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   rdg_cd_llb, trpd_leewv, &
   ptend, flx_heat)

   use coords_1d,  only: Coords1D
   use gw_rdg,     only: gw_rdg_src, gw_rdg_belowpeak, gw_rdg_break_trap, gw_rdg_do_vdiff
   use gw_common,  only: gw_drag_prof, energy_change

   character(len=5), intent(in) :: type         ! BETA or GAMMA
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   integer,          intent(in) :: n_rdg
   real(r8),         intent(in) :: dt           ! Time step.

   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   type(Coords1D),   intent(in) :: p               ! Pressure coordinates.
   real(r8),         intent(in) :: piln(ncol,pver+1)  ! Log of interface pressures.
   real(r8),         intent(in) :: zm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pver+1) ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: nm(ncol,pver)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: ni(ncol,pver+1) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8),         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
   real(r8),         intent(in) :: q(:,:,:)        ! Constituent array.
   real(r8),         intent(in) :: dse(ncol,pver)  ! Dry static energy.


   real(r8),         intent(in) :: effgw_rdg       ! Tendency efficiency.
   real(r8),         intent(in) :: effgw_rdg_max
   real(r8),         intent(in) :: hwdth(ncol,prdg) ! width of ridges.
   real(r8),         intent(in) :: clngt(ncol,prdg) ! length of ridges.
   real(r8),         intent(in) :: gbxar(ncol)      ! gridbox area

   real(r8),         intent(in) :: mxdis(ncol,prdg) ! Height estimate for ridge (m).
   real(r8),         intent(in) :: angll(ncol,prdg) ! orientation of ridges.
   real(r8),         intent(in) :: anixy(ncol,prdg) ! Anisotropy parameter.

   real(r8),         intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
   logical,          intent(in) :: trpd_leewv

   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies.

   real(r8),        intent(out) :: flx_heat(pcols)

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: c(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! horiz wavenumber [anisotropic orography].
   real(r8) :: kwvrdg(ncol)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(r8) :: ubm(ncol,pver)
   real(r8) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(r8) :: xv(ncol)
   real(r8) :: yv(ncol)

   ! Averages over source region.
   real(r8) :: ubmsrc(ncol) ! On-ridge wind.
   real(r8) :: usrc(ncol)   ! Zonal wind.
   real(r8) :: vsrc(ncol)   ! Meridional wind.
   real(r8) :: nsrc(ncol)   ! B-V frequency.
   real(r8) :: rsrc(ncol)   ! Density.

   ! normalized wavenumber
   real(r8) :: m2src(ncol)

   ! Top of low-level flow layer.
   real(r8) :: tlb(ncol)

   ! Bottom of linear wave region.
   real(r8) :: bwv(ncol)

   ! Froude numbers for flow/drag regimes
   real(r8) :: Fr1(ncol)
   real(r8) :: Fr2(ncol)
   real(r8) :: Frx(ncol)

   ! Wave Reynolds stresses at source level
   real(r8) :: tauoro(ncol)
   real(r8) :: taudsw(ncol)

   ! Surface streamline displacement height for linear waves.
   real(r8) :: hdspwv(ncol)

   ! Surface streamline displacement height for downslope wind regime.
   real(r8) :: hdspdw(ncol)

   ! Wave breaking level
   real(r8) :: wbr(ncol)

   real(r8) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(r8) :: taurx(ncol,pver+1)
   real(r8) :: taurx0(ncol,pver+1)
   real(r8) :: taury(ncol,pver+1)
   real(r8) :: taury0(ncol,pver+1)

   ! U,V tendency accumulators
   real(r8) :: utrdg(ncol,pver)
   real(r8) :: vtrdg(ncol,pver)
   real(r8) :: ttrdg(ncol,pver)

   ! Energy change used by fixer.
   real(r8) :: de(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,band_oro%ngwv:band_oro%ngwv,pver+1))
   allocate(gwut(ncol,pver,band_oro%ngwv:band_oro%ngwv))
   allocate(c(ncol,band_oro%ngwv:band_oro%ngwv))

   ! initialize accumulated momentum fluxes and tendencies
   taurx = 0._r8
   taury = 0._r8 
   ttrdg = 0._r8
   utrdg = 0._r8
   vtrdg = 0._r8

   do nn = 1, n_rdg
  
      kwvrdg  = 0.001_r8 / ( hwdth(:,nn) + 0.001_r8 ) ! this cant be done every time step !!!
      isoflag = 0   
      effgw   = effgw_rdg * ( hwdth(1:ncol,nn)* clngt(1:ncol,nn) ) / gbxar(1:ncol)
      effgw   = min( effgw_rdg_max , effgw )

      call gw_rdg_src(ncol, band_oro, p, &
         u, v, t, mxdis(:,nn), angll(:,nn), anixy(:,nn), kwvrdg, isoflag, zi, nm, &
         src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv,  & 
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c)

      call gw_rdg_belowpeak(ncol, band_oro, rdg_cd_llb, &
         t, mxdis(:,nn), anixy(:,nn), kwvrdg, & 
         zi, nm, ni, rhoi, &
         src_level, tau, & 
         ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, & 
         tauoro, taudsw, hdspwv, hdspdw)

      call gw_rdg_break_trap(ncol, band_oro, &
         zi, nm, ni, ubm, ubi, rhoi, kwvrdg , bwv, tlb, wbr, & 
         src_level, tlb_level, hdspwv, hdspdw,  mxdis(:,nn), & 
         tauoro, taudsw, tau, & 
         ldo_trapped_waves=trpd_leewv)
     
      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
         t, vramp,    &
         piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
         effgw, c, kvtt, q, dse, tau, utgw, vtgw, &
         ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, &
         kwvrdg=kwvrdg, & 
         satfac_in = 1._r8, lapply_vdiff=gw_rdg_do_vdiff )

      ! Add the tendencies from each ridge to the totals.
      do k = 1, pver
         ! diagnostics
         utrdg(:,k) = utrdg(:,k) + utgw(:,k)
         vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
         ttrdg(:,k) = ttrdg(:,k) + ttgw(:,k)
         ! physics tendencies
         ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
         ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
         ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
      end do

      do m = 1, pcnst
         do k = 1, pver
            ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
         end do
      end do

      do k = 1, pver+1
         taurx0(:,k) =  tau(:,0,k)*xv
         taury0(:,k) =  tau(:,0,k)*yv
         taurx(:,k)  =  taurx(:,k) + taurx0(:,k)
         taury(:,k)  =  taury(:,k) + taury0(:,k)
      end do

      if (nn == 1) then
         call outfld('TAU1RDG'//trim(type)//'M', tau(:,0,:),  ncol, lchnk)
         call outfld('UBM1'//trim(type),         ubm,         ncol, lchnk)
         call outfld('UBT1RDG'//trim(type),      gwut,        ncol, lchnk)
      end if

      if (nn <= 6) then
         write(cn, '(i1)') nn
         call outfld('TAU'//cn//'RDG'//trim(type)//'X', taurx0,  ncol, lchnk)
         call outfld('TAU'//cn//'RDG'//trim(type)//'Y', taury0,  ncol, lchnk)
         call outfld('UT'//cn//'RDG'//trim(type),       utgw,    ncol, lchnk)
         call outfld('VT'//cn//'RDG'//trim(type),       vtgw,    ncol, lchnk)
      end if

   end do ! end of loop over multiple ridges

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de

   call outfld('TAUARDG'//trim(type)//'X', taurx,  ncol, lchnk)
   call outfld('TAUARDG'//trim(type)//'Y', taury,  ncol, lchnk)

   if (trim(type) == 'BETA') then
      fname(1) = 'TAUGWX'
      fname(2) = 'TAUGWY'
      fname(3) = 'UTGWORO'
      fname(4) = 'VTGWORO'
   else if (trim(type) == 'GAMMA') then
      fname(1) = 'TAURDGGMX'
      fname(2) = 'TAURDGGMY'
      fname(3) = 'UTRDGGM'
      fname(4) = 'VTRDGGM'
   else
      call endrun('gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
                  //' type= '//type)
   end if

   call outfld(fname(1), taurx(:,pver+1), ncol, lchnk)
   call outfld(fname(2), taury(:,pver+1), ncol, lchnk)
   call outfld(fname(3), utrdg,  ncol, lchnk)
   call outfld(fname(4), vtrdg,  ncol, lchnk)
   call outfld('TTGWORO', ttrdg / cpair,  ncol, lchnk)

   deallocate(tau, gwut, c)

end subroutine gw_rdg_calc

!==========================================================================

! Add all history fields for a gravity wave spectrum source.
subroutine gw_spec_addflds(prefix, scheme, band, history_defaults)
  use cam_history, only: addfld, add_default, register_vector_field

  !------------------------------Arguments--------------------------------

  ! One character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Gravity wave scheme name prepended to output field descriptions.
  character(len=*), intent(in) :: scheme
  ! Wave speeds.
  type(GWBand), intent(in) :: band
  ! Whether or not to call add_default for fields output by WACCM.
  logical, intent(in) :: history_defaults

  !---------------------------Local storage-------------------------------

  integer :: l
  ! 7 chars is enough for "-100.00"
  character(len=7)  :: fnum
  ! 10 chars is enough for "BTAUXSn32"
  character(len=10) :: dumc1x, dumc1y
  ! Allow 80 chars for description
  character(len=80) dumc2

  !-----------------------------------------------------------------------

  ! Overall wind tendencies.
  call addfld (trim(prefix)//'UTGWSPEC',(/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency - gravity wave spectrum')
  call addfld (trim(prefix)//'VTGWSPEC',(/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' V tendency - gravity wave spectrum')
  call register_vector_field(trim(prefix)//'UTGWSPEC',trim(prefix)//'VTGWSPEC')

  call addfld (trim(prefix)//'TTGWSPEC',(/ 'lev' /), 'A','K/s', &
       trim(scheme)//' T tendency - gravity wave spectrum')

  ! Wind tendencies broken across five spectral bins.
  call addfld (trim(prefix)//'UTEND1',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   c < -40')
  call addfld (trim(prefix)//'UTEND2',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency  -40 < c < -15')
  call addfld (trim(prefix)//'UTEND3',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency  -15 < c <  15')
  call addfld (trim(prefix)//'UTEND4',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   15 < c <  40')
  call addfld (trim(prefix)//'UTEND5',  (/ 'lev' /), 'A','m/s2', &
       trim(scheme)//' U tendency   40 < c ')

  ! Reynold's stress toward each cardinal direction, and net zonal stress.
  call addfld (trim(prefix)//'TAUE' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Eastward Reynolds stress')
  call addfld (trim(prefix)//'TAUW' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Westward Reynolds stress')
  call addfld (trim(prefix)//'TAUNET' , (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' E+W Reynolds stress')
  call addfld (trim(prefix)//'TAUN' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Northward Reynolds stress')
  call addfld (trim(prefix)//'TAUS' ,   (/ 'ilev' /), 'A','Pa', &
       trim(scheme)//' Southward Reynolds stress')

  ! Momentum flux in each direction.
  call addfld (trim(prefix)//'EMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Eastward MF')
  call addfld (trim(prefix)//'WMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Westward MF')
  call addfld (trim(prefix)//'NMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Northward MF')
  call addfld (trim(prefix)//'SMF',       (/ 'lev' /), 'A','Pa', &
       trim(scheme)//' Southward MF')

  ! Temperature tendency terms.
  call addfld (trim(prefix)//'TTGWSDF' , (/ 'lev' /), 'A','K/s', &
       trim(scheme)//' t tendency - diffusion term')
  call addfld (trim(prefix)//'TTGWSKE' , (/ 'lev' /), 'A','K/s', &
       trim(scheme)//' t tendency - kinetic energy conversion term')

  ! Gravity wave source spectra by wave number.
  do l=-band%ngwv,band%ngwv
     ! String containing reference speed.
     write (fnum,fmt='(f7.2)') band%cref(l)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
     dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m/s"
     call addfld (trim(dumc1x),(/ 'lev' /), 'A','Pa',dumc2)
     call addfld (trim(dumc1y),(/ 'lev' /), 'A','Pa',dumc2)

  end do

  if (history_defaults) then
     call add_default(trim(prefix)//'UTGWSPEC', 1, ' ')
     call add_default(trim(prefix)//'VTGWSPEC', 1, ' ')
     call add_default(trim(prefix)//'TTGWSPEC', 1, ' ')
     call add_default(trim(prefix)//'TAUE', 1, ' ')
     call add_default(trim(prefix)//'TAUW', 1, ' ')
     call add_default(trim(prefix)//'TAUNET', 1, ' ')
     call add_default(trim(prefix)//'TAUN', 1, ' ')
     call add_default(trim(prefix)//'TAUS', 1, ' ')
  end if

end subroutine gw_spec_addflds

!==========================================================================

! Outputs for spectral waves.
subroutine gw_spec_outflds(prefix, lchnk, ncol, band, c, u, v, xv, yv, &
     gwut, dttdf, dttke, tau, utgw, vtgw, ttgw, taucd)

  use gw_common, only: west, east, south, north

  ! One-character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Chunk and number of columns in the chunk.
  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  ! Wave speeds.
  type(GWBand), intent(in) :: band
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Winds at cell midpoints.
  real(r8), intent(in) :: u(ncol,pver)
  real(r8), intent(in) :: v(ncol,pver)
  ! Unit vector in the direction of wind at source level.
  real(r8), intent(in) :: xv(ncol)
  real(r8), intent(in) :: yv(ncol)
  ! Wind tendency for each wave.
  real(r8), intent(in) :: gwut(ncol,pver,-band%ngwv:band%ngwv)
  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(ncol,pver)
  real(r8) :: dttke(ncol,pver)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver)
  ! Zonal and meridional total wind tendencies.
  real(r8), intent(in) :: utgw(ncol,pver)
  real(r8), intent(in) :: vtgw(ncol,pver)
  ! Temperature tendencies.
  real(r8), intent(in) :: ttgw(ncol,pver)
  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8), intent(in) :: taucd(ncol,pver+1,4)

  ! Indices
  integer :: i, k, l
  integer :: ix(ncol, -band%ngwv:band%ngwv), iy(ncol, -band%ngwv:band%ngwv)
  integer :: iu(ncol), iv(ncol)

  ! Zonal wind tendency, broken up into five bins.
  real(r8) :: utb(ncol, pver, 5)
  ! Definition of the bin boundaries.
  real(r8), parameter :: bounds(4) = (/ -40._r8, -15._r8, &
       15._r8, 40._r8 /)

  ! Momentum flux in the four cardinal directions.
  real(r8) :: mf(ncol, pver, 4)

  ! Wave stress in zonal/meridional direction
  real(r8) :: taux(ncol,-band%ngwv:band%ngwv,pver)
  real(r8) :: tauy(ncol,-band%ngwv:band%ngwv,pver)

  ! Temporaries for output
  real(r8) :: dummyx(ncol,pver)
  real(r8) :: dummyy(ncol,pver)
  ! Variable names
  character(len=10) :: dumc1x, dumc1y


  ! Accumulate wind tendencies binned according to phase speed.

  utb = 0._r8

  ! Find which output bin the phase speed corresponds to.
  ix = find_bin(c)

  ! Put the wind tendency in that bin.
  do l = -band%ngwv, band%ngwv
     do k = 1, pver
        do i = 1, ncol
           utb(i,k,ix(i,l)) = utb(i,k,ix(i,l)) + gwut(i,k,l)
        end do
     end do
  end do

  ! Find just the zonal part.
  do l = 1, 5
     do k = 1, pver
        utb(:, k, l) = utb(:, k, l) * xv
     end do
  end do

  call outfld(trim(prefix)//'UTEND1', utb(:,:,1), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND2', utb(:,:,2), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND3', utb(:,:,3), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND4', utb(:,:,4), ncol, lchnk)
  call outfld(trim(prefix)//'UTEND5', utb(:,:,5), ncol, lchnk)

  ! Output temperature tendencies due to diffusion and from kinetic energy.
  call outfld(trim(prefix)//'TTGWSDF', dttdf / cpair, ncol, lchnk)
  call outfld(trim(prefix)//'TTGWSKE', dttke / cpair, ncol, lchnk)


  ! Output tau broken down into zonal and meridional components.

  taux = 0._r8
  tauy = 0._r8

  ! Project c, and convert each component to a wavenumber index.
  ! These are mappings from the wavenumber index of tau to those of taux
  ! and tauy, respectively.
  do l=-band%ngwv,band%ngwv
     ix(:,l) = c_to_l(c(:,l)*xv)
     iy(:,l) = c_to_l(c(:,l)*yv)
  end do

  ! Find projection of tau.
  do k = 1, pver
     do l = -band%ngwv,band%ngwv
        do i = 1, ncol
           taux(i,ix(i,l),k) = taux(i,ix(i,l),k) &
                + abs(tau(i,l,k)*xv(i))
           tauy(i,iy(i,l),k) = tauy(i,iy(i,l),k) &
                + abs(tau(i,l,k)*yv(i))
        end do
     end do
  end do

  do l=-band%ngwv,band%ngwv

     dummyx = taux(:,l,:)
     dummyy = tauy(:,l,:)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)

     call outfld(dumc1x,dummyx,ncol,lchnk)
     call outfld(dumc1y,dummyy,ncol,lchnk)

  enddo


  ! Output momentum flux in each cardinal direction.
  mf = 0._r8

  do k = 1, pver

     ! Convert wind speed components to wavenumber indices.
     iu = c_to_l(u(:,k))
     iv = c_to_l(v(:,k))

     ! Sum tau components in each cardinal direction.
     ! Split west/east and north/south based on whether wave speed exceeds
     ! wind speed.
     do l = -band%ngwv, band%ngwv

        where (iu > l)
           mf(:,k,west) = mf(:,k,west) + taux(:,l,k)
        elsewhere
           mf(:,k,east) = mf(:,k,east) + taux(:,l,k)
        end where

        where (iv > l)
           mf(:,k,south) = mf(:,k,south) + tauy(:,l,k)
        elsewhere
           mf(:,k,north) = mf(:,k,north) + tauy(:,l,k)
        end where

     end do

  end do

  call outfld(trim(prefix)//'WMF',mf(:,:,west),ncol,lchnk)
  call outfld(trim(prefix)//'EMF',mf(:,:,east),ncol,lchnk)
  call outfld(trim(prefix)//'SMF',mf(:,:,south),ncol,lchnk)
  call outfld(trim(prefix)//'NMF',mf(:,:,north),ncol,lchnk)

  ! Simple output fields written to history file.
  ! Total wind tendencies.
  call outfld (trim(prefix)//'UTGWSPEC', utgw , ncol, lchnk)
  call outfld (trim(prefix)//'VTGWSPEC', vtgw , ncol, lchnk)
  call outfld (trim(prefix)//'TTGWSPEC', ttgw , ncol, lchnk)

  ! Tau in each direction.
  call outfld (trim(prefix)//'TAUE', taucd(:,:,east), ncol, lchnk)
  call outfld (trim(prefix)//'TAUW', taucd(:,:,west), ncol, lchnk)
  call outfld (trim(prefix)//'TAUN', taucd(:,:,north), ncol, lchnk)
  call outfld (trim(prefix)//'TAUS', taucd(:,:,south), ncol, lchnk)

  call outfld (trim(prefix)//'TAUNET', taucd(:,:,east)+taucd(:,:,west), &
       ncol, lchnk)

contains

  ! Given a value, finds which bin marked by "bounds" the value falls
  ! into.
  elemental function find_bin(val) result(idx)
    real(r8), intent(in) :: val

    integer :: idx

    ! We just have to count how many bounds are exceeded.
    if (val >= 0._r8) then
       idx = count(val > bounds) + 1
    else
       idx = count(val >= bounds) + 1
    end if

  end function find_bin

  ! Convert a speed to a wavenumber between -ngwv and ngwv.
  elemental function c_to_l(c) result(l)
    real(r8), intent(in) :: c

    integer :: l

    l = min( max(int(c/band%dc),-band%ngwv), band%ngwv )

  end function c_to_l

end subroutine gw_spec_outflds

!==========================================================================

! Generates names for tau output across the wave spectrum (e.g.
! BTAUXSn01 or TAUYSp05).
! Probably this should use a wavenumber dimension on one field rather
! than creating a ton of numbered fields.
character(len=9) pure function tau_fld_name(l, prefix, x_not_y)
  ! Wavenumber
  integer, intent(in) :: l
  ! Single-character prefix for output
  character(len=1), intent(in) :: prefix
  ! X or Y?
  logical, intent(in) :: x_not_y

  character(len=2) :: num_str

  tau_fld_name = trim(prefix)

  tau_fld_name = trim(tau_fld_name)//"TAU"

  if (x_not_y) then
     tau_fld_name = trim(tau_fld_name)//"XS"
  else
     tau_fld_name = trim(tau_fld_name)//"YS"
  end if

  if (l < 0) then
     tau_fld_name = trim(tau_fld_name)//"n"
  else
     tau_fld_name = trim(tau_fld_name)//"p"
  end if

  write(num_str,'(I2.2)') abs(l)

  tau_fld_name = trim(tau_fld_name)//num_str

end function tau_fld_name

!==========================================================================

end module gw_drag
