! CAM interface for the gravity waves drag parameterization
! now moved to CCPP
module gw_drag_cam
  use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
  use shr_log_mod,    only: shr_errMsg => shr_log_errMsg
  use shr_assert_mod, only: shr_assert

  use ppgrid,         only: pcols, pver, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,     only: masterproc
  use cam_history,    only: outfld
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use error_messages, only: alloc_err

  use ref_pres,       only: do_molec_diff, nbot_molec
  use physconst,      only: cpair

  ! These are the actual switches for different gravity wave sources.
  use phys_control,   only: use_gw_oro, use_gw_front, use_gw_front_igw, &
                            use_gw_convect_dp, use_gw_convect_sh,       &
                            use_simple_phys, use_gw_movmtn_pbl, phys_getopts

  use physics_buffer, only: pbuf_get_field

  implicit none
  private
  save

!
! PUBLIC: interfaces
!
  public :: gw_drag_cam_readnl           ! Read namelist
  public :: gw_drag_cam_init             ! Initialization
  public :: gw_drag_cam_tend             ! interface to actual parameterization

!
! PRIVATE: Rest of the data and interfaces are private to this module
!
  real(r8), parameter :: unset_r8 = huge(1._r8)

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
  ! PBL moving mtn
  real(r8) :: effgw_movmtn_pbl = unset_r8
  integer  :: movmtn_ksteer  = -1
  integer  :: movmtn_klaunch = -1
  integer  :: movmtn_source  = -1
  real(r8) :: movmtn_psteer  = unset_r8
  real(r8) :: movmtn_plaunch = unset_r8

  ! Parameters controlling isotropic residual
  ! orographic GW.
  logical :: use_gw_rdg_resid = .false.
  real(r8) :: effgw_rdg_resid = unset_r8

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
  character(len=cl) :: gw_drag_file = ""
  character(len=cl) :: gw_drag_file_sh = ""
  character(len=cl) :: gw_drag_file_mm = ""

  ! Indices into pbuf
  integer :: kvt_idx      = -1
  integer :: ttend_dp_idx = -1
  integer :: ttend_sh_idx = -1
  integer :: frontgf_idx  = -1
  integer :: frontga_idx  = -1

  integer :: vort4gw_idx  = -1

  integer :: sgh_idx      = -1

  ! From CLUBB
  integer :: ttend_clubb_idx  = -1
  integer :: upwp_clubb_gw_idx   = -1
  integer :: vpwp_clubb_gw_idx   = -1
  integer :: thlp2_clubb_gw_idx  = -1
  integer :: wpthlp_clubb_gw_idx  = -1

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

  ! Width of gaussian used to create frontogenesis tau profile [m s-1].
  real(r8) :: front_gaussian_width = -huge(1._r8)

  real(r8) :: alpha_gw_movmtn

  logical :: gw_top_taper=.false.

  ! Maximum wave number and width of spectrum bins.
  integer :: pgwv = -1
  real(r8) :: gw_dc = unset_r8
  integer :: pgwv_long = -1
  real(r8) :: gw_dc_long = unset_r8

  ! Temperature change due to deep convection.
  real(r8), pointer :: ttend_dp(:,:)
  ! Temperature change due to shallow convection.
  real(r8), pointer :: ttend_sh(:,:)

  !  New couplings from CLUBB
  real(r8), pointer :: ttend_clubb(:,:)
  real(r8), pointer :: thlp2_clubb_gw(:,:)
  real(r8), pointer :: wpthlp_clubb_gw(:,:)
  real(r8), pointer :: upwp_clubb_gw(:,:)
  real(r8), pointer :: vpwp_clubb_gw(:,:)
  real(r8), pointer :: vort4gw(:,:)

  ! Gravity wave Ridge scheme namelist
  logical  :: gw_rdg_do_divstream, gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
              gw_rdg_do_backward_compat

  logical  :: gw_rdg_do_vdiff = .true.

  real(r8) :: gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
              gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
              gw_rdg_orohmin, gw_rdg_orovmin, gw_rdg_orostratmin, gw_rdg_orom2min

  ! Gravity wave Ridge scheme data
  ! this is grid dependent and stored by chunk at initialization (last dimension is lchnk)
  ! the original variants are used for Meso Beta and stored in bnd_topo_file
  ! the "g" variants are used for Meso Gamma and stored in bnd_rdggm_file
  integer, parameter :: prdg = 16
  real(r8), pointer, dimension(:,:)   :: rdg_gbxar, rdg_gbxarg
  real(r8), pointer, dimension(:,:)   :: rdg_isovar!, rdg_isovarg
  real(r8), pointer, dimension(:,:)   :: rdg_isowgt!, rdg_isowgtg
  real(r8), pointer, dimension(:,:,:) :: rdg_hwdth, rdg_hwdthg
  real(r8), pointer, dimension(:,:,:) :: rdg_clngt, rdg_clngtg
  real(r8), pointer, dimension(:,:,:) :: rdg_mxdis, rdg_mxdisg
  real(r8), pointer, dimension(:,:,:) :: rdg_anixy, rdg_anixyg
  real(r8), pointer, dimension(:,:,:) :: rdg_angll, rdg_angllg

  ! State vramp
  real(r8), pointer :: vramp(:) => null()


!==========================================================================
contains
!==========================================================================

subroutine gw_drag_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, &
                             mpi_character, mpi_logical, mpi_integer

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gw_drag_cam_readnl'

  namelist /gw_drag_nl/ pgwv, gw_dc, pgwv_long, gw_dc_long, tau_0_ubc, &
       effgw_beres_dp, effgw_beres_sh, effgw_cm, effgw_cm_igw, effgw_oro, &
       frontgfc, gw_drag_file, gw_drag_file_sh, gw_drag_file_mm, taubgnd, &
       taubgnd_igw, gw_polar_taper, &
       use_gw_rdg_beta, n_rdg_beta, effgw_rdg_beta, effgw_rdg_beta_max, &
       rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
       use_gw_rdg_gamma, n_rdg_gamma, effgw_rdg_gamma, effgw_rdg_gamma_max, &
       rdg_gamma_cd_llb, trpd_leewv_rdg_gamma, bnd_rdggm, &
       gw_oro_south_fac, gw_limit_tau_without_eff, &
       gw_lndscl_sgh, gw_prndl, gw_apply_tndmax, gw_qbo_hdepth_scaling, &
       gw_top_taper, front_gaussian_width, alpha_gw_movmtn, use_gw_rdg_resid, &
       effgw_rdg_resid, effgw_movmtn_pbl, movmtn_source, movmtn_psteer, &
       movmtn_plaunch

   namelist /gw_rdg_nl/ gw_rdg_do_divstream, gw_rdg_C_BetaMax_DS, gw_rdg_C_GammaMax, &
                   gw_rdg_Frx0, gw_rdg_Frx1, gw_rdg_C_BetaMax_SM, gw_rdg_Fr_c, &
                   gw_rdg_do_smooth_regimes, gw_rdg_do_adjust_tauoro, &
                   gw_rdg_do_backward_compat, gw_rdg_orohmin, gw_rdg_orovmin, &
                   gw_rdg_orostratmin, gw_rdg_orom2min, gw_rdg_do_vdiff

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
  call mpi_bcast(gw_drag_file_mm, len(gw_drag_file_mm), mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_drag_file_mm")

  call mpi_bcast(front_gaussian_width, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: front_gaussian_width")

  call mpi_bcast(alpha_gw_movmtn, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: alpha_gw_movmtn")
  call mpi_bcast(effgw_movmtn_pbl, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_movmtn_pbl")
  call mpi_bcast(movmtn_source, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: movmtn_source")
  call mpi_bcast(movmtn_psteer, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: movmtn_psteer")
  call mpi_bcast(movmtn_plaunch, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: movmtn_plaunch")

  call mpi_bcast(use_gw_rdg_resid, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_gw_rdg_resid")
  call mpi_bcast(effgw_rdg_resid, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: effgw_rdg_resid")

  ! Check if pgwv was set.
  call shr_assert(pgwv >= 0, &
       "gw_drag_cam_readnl: pgwv must be set via the namelist and &
       &non-negative."// &
       shr_errMsg(__FILE__, __LINE__))

  ! Check if gw_dc was set.
  call shr_assert(gw_dc /= unset_r8, &
       "gw_drag_cam_readnl: gw_dc must be set via the namelist."// &
       shr_errMsg(__FILE__, __LINE__))

  if (use_gw_rdg_gamma .or. use_gw_rdg_beta) then
     if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'gw_rdg_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, gw_rdg_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(sub // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

      ! Broadcast the local variables

      call mpi_bcast(gw_rdg_do_divstream, 1, mpi_logical, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_divstream")
      call mpi_bcast(gw_rdg_do_smooth_regimes, 1, mpi_logical, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_smooth_regimes")
      call mpi_bcast(gw_rdg_do_adjust_tauoro, 1, mpi_logical, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_adjust_tauoro")
      call mpi_bcast(gw_rdg_do_backward_compat, 1, mpi_logical, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_backward_compat")

      call mpi_bcast(gw_rdg_C_BetaMax_DS, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_BetaMax_DS")
      call mpi_bcast(gw_rdg_C_GammaMax, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_GammaMax")
      call mpi_bcast(gw_rdg_Frx0, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Frx0")
      call mpi_bcast(gw_rdg_Frx1, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Frx1")
      call mpi_bcast(gw_rdg_C_BetaMax_SM, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_C_BetaMax_SM")
      call mpi_bcast(gw_rdg_Fr_c, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_Fr_c")
      call mpi_bcast(gw_rdg_orohmin, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orohmin")
      call mpi_bcast(gw_rdg_orovmin, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orovmin")
      call mpi_bcast(gw_rdg_orostratmin, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orostratmin")
      call mpi_bcast(gw_rdg_orom2min, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_orom2min")

      call mpi_bcast(gw_rdg_do_vdiff, 1, mpi_logical, mstrid, mpicom, ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: gw_rdg_do_vdiff")

      if (gw_rdg_Fr_c > 1.0_r8) call endrun(sub//": FATAL: gw_rdg_Fr_c must be <= 1")
  end if

end subroutine gw_drag_cam_readnl

!==========================================================================

subroutine gw_drag_cam_init()

  !-----------------------------------------------------------------------
  ! Time independent initialization for multiple gravity wave
  ! parameterization.
  !-----------------------------------------------------------------------

  use cam_history,      only: addfld, add_default, horiz_only
  use cam_history,      only: register_vector_field
  use air_composition,  only: cpairv
  use interpolate_data, only: lininterp
  use phys_control,     only: phys_getopts
  use physics_buffer,   only: pbuf_get_index
  use constituents,     only: cnst_get_ind

  use cam_initfiles,    only: topo_file_get_id

  ! temporary for restart with ridge scheme
  use cam_initfiles,    only: bnd_topo

  use cam_grid_support, only: cam_grid_check, cam_grid_id
  use cam_grid_support, only: cam_grid_get_dim_names
  use pio,              only: file_desc_t, PIO_NOWRITE, pio_closefile
  use cam_pio_utils,    only: cam_pio_openfile
  use ncdio_atm,        only: infld
  use ioFileMod,        only: getfil

  use ref_pres,   only: pref_edge, pref_mid
  use physconst,  only: gravit, rair, rearth, pi

  use ppgrid, only: pcols

  use gw_common,                         only: gravity_wave_drag_common_init

  use gravity_wave_drag_orographic,      only: gravity_wave_drag_orographic_init
  use gravity_wave_drag_top_taper,       only: gravity_wave_drag_top_taper_init
  use gravity_wave_drag_frontogenesis,   only: gravity_wave_drag_frontogenesis_init
  use gravity_wave_drag_ridge,           only: gw_rdg_init
  use gravity_wave_drag_moving_mountain, only: gravity_wave_drag_moving_mountain_init
  use gravity_wave_drag_convection,      only: gravity_wave_drag_convection_init

  !---------------------------Local storage-------------------------------

  integer          :: i, l, k, lchnk
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

  ! Read data from file
  type(file_desc_t), pointer :: fh_topo
  type(file_desc_t)  :: fh_rdggm
  integer            :: grid_id
  character(len=8)   :: dim1name, dim2name
  logical            :: found
  logical            :: found_rdggm
  character(len=cl) :: bnd_rdggm_loc   ! filepath of topo file on local disk

  ! Allow reporting of error messages.
  character(len=128) :: errstring
  character(len=*), parameter :: sub = 'gw_init'

  ! temporary workaround for restart w/ ridge scheme
  character(len=cl) :: bnd_topo_loc   ! filepath of topo file on local disk

  character(len=512)              :: errmsg
  integer                         :: errflg
  character(len=cl)              :: gw_drag_file_loc
  character(len=cl)              :: gw_drag_file_mm_loc
  character(len=cl)              :: gw_drag_file_sh_loc

  !-----------------------------------------------------------------------

  if (use_gw_convect_dp) then
     if (trim(gw_drag_file) /= "") then
        call getfil(gw_drag_file, gw_drag_file_loc)
     else
        call endrun('ERROR gw_drag_cam_init: use_gw_convect_dp is true but gw_drag_file is not specified in namelist')
     end if
  end if
  if (use_gw_convect_sh) then
     if (trim(gw_drag_file_sh) /= "") then
        call getfil(gw_drag_file_sh, gw_drag_file_sh_loc)
     else
        call endrun('ERROR gw_drag_cam_init: use_gw_convect_sh is true but gw_drag_file_sh is not specified in namelist')
     end if
  end if
  !movmtn files
  if (use_gw_movmtn_pbl) then
     if( trim(gw_drag_file_mm) /= "") then
        call getfil(gw_drag_file_mm, gw_drag_file_mm_loc)
     else
        call endrun('ERROR gw_drag_cam_init: use_gw_movmtn_pbl is true but gw_drag_file_mm is not specified in namelist')
     end if
  end if

  !rdg files
  if (use_gw_rdg_beta) then
     if (trim(bnd_topo) /= "") then
        call getfil(bnd_topo, bnd_topo_loc)
     else
        call endrun('ERROR gw_drag_cam_init: use_gw_rdg_beta is true but bnd_topo file is not specified in namelist')
     end if
  end if
  if (use_gw_rdg_gamma) then
     if ( trim(bnd_rdggm) /= "") then
        call getfil(bnd_rdggm, bnd_rdggm_loc)
     else
        call endrun('ERROR gw_drag_cam_init: use_gw_rdg_gamma is true but bnd_rdggm file is not specified in namelist')
     end if
  end if

  ! Sanity checks
  if(use_gw_oro) then
    if (effgw_oro == unset_r8) then
      call endrun("gw_init: Orographic gravity waves enabled, but effgw_oro was not set.")
    end if
  end if

  if(use_gw_front .or. use_gw_front_igw) then
    if(frontgfc == unset_r8) then
      call endrun("gw_init: Frontogenesis enabled, but frontgfc was not set!")
    endif
  endif

  ! pbuf index initialization
  kvt_idx = pbuf_get_index('kvt')

  if(use_gw_front .or. use_gw_front_igw) then
     frontgf_idx = pbuf_get_index('FRONTGF')
     frontga_idx = pbuf_get_index('FRONTGA')
  endif

  if(use_gw_movmtn_pbl) then
     ! get pbuf indices for CLUBB couplings
     ttend_clubb_idx     = pbuf_get_index('TTEND_CLUBB')
     thlp2_clubb_gw_idx  = pbuf_get_index('THLP2_CLUBB_GW')
     upwp_clubb_gw_idx   = pbuf_get_index('UPWP_CLUBB_GW')
     vpwp_clubb_gw_idx   = pbuf_get_index('VPWP_CLUBB_GW')
     wpthlp_clubb_gw_idx = pbuf_get_index('WPTHLP_CLUBB_GW')
     vort4gw_idx         = pbuf_get_index('VORT4GW')

     ttend_dp_idx        = pbuf_get_index('TTEND_DP')
  endif

  if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     sgh_idx = pbuf_get_index('SGH')
  endif

  if (use_gw_convect_dp) ttend_dp_idx    = pbuf_get_index('TTEND_DP')
  if (use_gw_convect_sh) ttend_sh_idx    = pbuf_get_index('TTEND_SH')

  ! Output initialization status
  if(masterproc) then
    write(iulog,*) "gw_drag_cam_init: use_gw_movmtn_pbl = ", use_gw_movmtn_pbl

    write(iulog,*) "gw_drag_cam_init: use_gw_convect_dp = ", use_gw_convect_dp
    write(iulog,*) "gw_drag_cam_init: use_gw_convect_sh = ", use_gw_convect_sh

    write(iulog,*) "gw_drag_cam_init: use_gw_front = ", use_gw_front
    write(iulog,*) "gw_drag_cam_init: use_gw_front_igw = ", use_gw_front_igw

    write(iulog,*) "gw_drag_cam_init: use_gw_oro = ", use_gw_oro
    write(iulog,*) "gw_drag_cam_init: use_gw_rdg_beta = ", use_gw_rdg_beta
    write(iulog,*) "gw_drag_cam_init: use_gw_rdg_gamma = ", use_gw_rdg_gamma
  endif

  ! Call the CCPPized initialization subroutines for the common module.
  call gravity_wave_drag_common_init( &
       pver_in                      = pver, &
       amIRoot                      = masterproc, &
       iulog                        = iulog, &
       pref_edge                    = pref_edge, &
       tau_0_ubc_in                 = tau_0_ubc, &
       ktop_in                      = ktop, &
       gravit_in                    = gravit, &
       rair_in                      = rair, &
       prndl_in                     = gw_prndl, &
       qbo_hdepth_scaling_in        = gw_qbo_hdepth_scaling, &
       errmsg                       = errmsg, &
       errflg                       = errflg)

  ! Call the CCPPized initialization subroutines for individual parameterizations.
  if(use_gw_movmtn_pbl) then
    call gravity_wave_drag_moving_mountain_init(pver = pver, &
                        masterproc = masterproc, &
                        iulog = iulog, &
                        file_path = gw_drag_file_mm_loc, &
                        gw_delta_c = gw_dc, & ! mid
                        pref_edge = pref_edge, &
                        movmtn_psteer = movmtn_psteer, &
                        movmtn_plaunch = movmtn_plaunch, &
                        movmtn_source_nl = movmtn_source, &
                        errmsg = errmsg, errflg = errflg)
    if (errflg /= 0) call endrun(errmsg)
  endif

  if(use_gw_convect_dp .or. use_gw_convect_sh) then
    call gravity_wave_drag_convection_init( &
      pver                = pver, &
      pi                  = pi, &
      masterproc          = masterproc, &
      iulog               = iulog, &
      gw_drag_file_sh     = gw_drag_file_sh_loc, &
      gw_drag_file_dp     = gw_drag_file_loc, &
      pref_edge           = pref_edge, &
      gw_delta_c          = gw_dc, &
      pgwv                = pgwv, &
      use_gw_convect_dp   = use_gw_convect_dp, &
      use_gw_convect_sh   = use_gw_convect_sh, &
      errmsg              = errmsg, &
      errflg              = errflg)
    if (errflg /= 0) call endrun(errmsg)
  end if

  if(use_gw_oro) then
    call gravity_wave_drag_orographic_init( &
      gw_delta_c          = gw_dc, &
      effgw_oro           = effgw_oro, &
      errmsg              = errmsg, &
      errflg              = errflg)
    if(errflg /= 0) call endrun(errmsg)
  endif

  ! note: gw_limit_tau_without_eff appears unused.
  if(errflg /= 0) then
    call endrun("gw_drag_init: " // errmsg)
  endif

  if(use_gw_front .or. use_gw_front_igw) then
    call gravity_wave_drag_frontogenesis_init( &
         pver                 = pver, &
         pi                   = pi, &
         masterproc           = masterproc, &
         iulog                = iulog, &
         pref_edge            = pref_edge, &
         frontgfc             = frontgfc, &
         gw_delta_c           = gw_dc, &
         pgwv                 = pgwv, &
         pgwv_long            = pgwv_long, &
         taubgnd              = taubgnd, &
         taubgnd_igw          = taubgnd_igw, &
         effgw_cm             = effgw_cm, &
         effgw_cm_igw         = effgw_cm_igw, &
         use_gw_front         = use_gw_front, &
         use_gw_front_igw     = use_gw_front_igw, &
         front_gaussian_width = front_gaussian_width, &
         errmsg               = errmsg, &
         errflg               = errflg)
    if(errflg /= 0) then
      call endrun("gravity_wave_drag_frontogenesis_init: " // errmsg)
    endif
  endif

  ! Initialize tapering
  allocate(vramp(pver), stat=errflg)
  call gravity_wave_drag_top_taper_init( &
       pver = pver, &
       amIRoot = masterproc, &
       iulog = iulog, &
       gw_top_taper = gw_top_taper, &
       pref_edge = pref_edge, &
       pref_mid = pref_mid, &
       vramp = vramp(:pver), &
       errmsg = errmsg, &
       errflg = errflg)
  if(errflg /= 0) then
    call endrun(errmsg)
  endif

  ! For gravity wave ridge scheme, initialize and read data into model state.
  ! This has to be at the host/CAM interface level because of chunking.

  ! First, initialize the data to zeros to populate default values for unused columns
  allocate(rdg_gbxar(pcols, begchunk:endchunk), stat=errflg)
  allocate(rdg_isovar(pcols, begchunk:endchunk), stat=errflg)
  allocate(rdg_isowgt(pcols, begchunk:endchunk), stat=errflg)
  allocate(rdg_hwdth(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_clngt(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_mxdis(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_anixy(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_angll(pcols, prdg, begchunk:endchunk), stat=errflg)

  rdg_gbxar(:,:) = 0._r8
  rdg_isovar(:,:) = 0._r8
  rdg_isowgt(:,:) = 0._r8
  rdg_hwdth(:,:,:) = 0._r8
  rdg_clngt(:,:,:) = 0._r8
  rdg_mxdis(:,:,:) = 0._r8
  rdg_anixy(:,:,:) = 0._r8
  rdg_angll(:,:,:) = 0._r8

  allocate(rdg_gbxarg(pcols, begchunk:endchunk), stat=errflg)
  allocate(rdg_hwdthg(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_clngtg(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_mxdisg(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_anixyg(pcols, prdg, begchunk:endchunk), stat=errflg)
  allocate(rdg_angllg(pcols, prdg, begchunk:endchunk), stat=errflg)

  rdg_gbxarg(:,:) = 0._r8
  rdg_hwdthg(:,:,:) = 0._r8
  rdg_clngtg(:,:,:) = 0._r8
  rdg_mxdisg(:,:,:) = 0._r8
  rdg_anixyg(:,:,:) = 0._r8
  rdg_angllg(:,:,:) = 0._r8

  if(use_gw_rdg_beta .or. use_gw_rdg_gamma) then
    grid_id = cam_grid_id('physgrid')
    if(.not. cam_grid_check(grid_id)) call endrun(sub // ': ERROR: no "physgrid" grid')
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
  endif

  if(use_gw_rdg_beta) then

      ! The CCPPized initialization routine, gravity_wave_drag_beta_ridge_init,
      ! uses the I/O reader which is not chunk aware. Thus, we have to lift this
      ! up in to the host model layer and use infld to read it in;
      ! SIMA can still use the "native" routine with ncol only.
      fh_topo => topo_file_get_id()
      if (.not. associated(fh_topo)) then
         ! Try to open topo file here.  This workaround will not be needed
         ! once the refactored initialization sequence is on trunk.
         allocate(fh_topo)
         call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)
      end if

      call infld('GBXAR', fh_topo, dim1name, dim2name, 1, pcols, &
                          begchunk, endchunk, rdg_gbxar, found, gridname='physgrid')
      if (.not. found) call endrun(sub//': ERROR: GBXAR not found on topo file')
      rdg_gbxar = rdg_gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2

      call infld('ISOVAR', fh_topo, dim1name, dim2name, 1, pcols, &
                          begchunk, endchunk, rdg_isovar, found, gridname='physgrid')
      ! ++jtb - Temporary fix until topo files contain this variable
      if (.not. found) rdg_isovar(:,:) = 0._r8

      call infld('ISOWGT', fh_topo, dim1name, dim2name, 1, pcols, &
                          begchunk, endchunk, rdg_isowgt, found, gridname='physgrid')
      ! ++jtb - Temporary fix until topo files contain this variable
      if (.not. found) rdg_isowgt(:,:) = 0._r8

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
  endif

  ! Gamma ridge data initialization
  if(use_gw_rdg_gamma) then

      call cam_pio_openfile(fh_rdggm, bnd_rdggm_loc, PIO_NOWRITE)
      if (.not. use_gw_rdg_beta) then
        call infld('GBXAR', fh_rdggm, dim1name, dim2name, 1, pcols, &
                            begchunk, endchunk, rdg_gbxarg, found, gridname='physgrid')
        if (.not. found) call endrun(sub//': ERROR: GBXAR not found on bnd_rdggm')
        rdg_gbxarg = rdg_gbxarg * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2
      else
        rdg_gbxarg = rdg_gbxar
      endif

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
  endif

  ! Call the underlying ridge scheme initialization to populate namelist variables
  ! into module
  if(use_gw_rdg_beta .or. use_gw_rdg_gamma) then
    call gw_rdg_init( &
      gw_delta_c             = gw_dc, &
      effgw_rdg_beta         = effgw_rdg_beta, &
      effgw_rdg_gamma        = effgw_rdg_gamma, &
      use_gw_rdg_beta_in     = use_gw_rdg_beta, &
      use_gw_rdg_gamma_in    = use_gw_rdg_gamma, &
      gw_rdg_do_divstream_nl = gw_rdg_do_divstream, &
      gw_rdg_C_BetaMax_DS_nl = gw_rdg_C_BetaMax_DS, &
      gw_rdg_C_GammaMax_nl   = gw_rdg_C_GammaMax, &
      gw_rdg_Frx0_nl         = gw_rdg_Frx0, &
      gw_rdg_Frx1_nl         = gw_rdg_Frx1, &
      gw_rdg_C_BetaMax_SM_nl = gw_rdg_C_BetaMax_SM, &
      gw_rdg_Fr_c_nl         = gw_rdg_Fr_c, &
      gw_rdg_do_smooth_regimes_nl  = gw_rdg_do_smooth_regimes, &
      gw_rdg_do_adjust_tauoro_nl   = gw_rdg_do_adjust_tauoro, &
      gw_rdg_do_backward_compat_nl = gw_rdg_do_backward_compat, &
      gw_rdg_orohmin_nl      = gw_rdg_orohmin, &
      gw_rdg_orovmin_nl      = gw_rdg_orovmin, &
      gw_rdg_orostratmin_nl  = gw_rdg_orostratmin, &
      gw_rdg_orom2min_nl     = gw_rdg_orom2min, &
      gw_rdg_do_vdiff_nl     = gw_rdg_do_vdiff, &
      errmsg = errmsg, &
      errflg = errflg)
    if (errflg /= 0) return
  endif

  !--------------------------------------------------
  ! CAM-specific diagnostic initialization code.
  ! Code below should not contain any physics logic, and should be for diagnostics only.
  ! Similar code is replicated in sima_diagnostics for CAM-SIMA.
  !--------------------------------------------------

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm, &
       history_amwg_out   = history_amwg  )

  ! Orographic terms.
  if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     call gw_drag_cam_oro_diag_init()
  end if

  ! Frontogenesis terms.
  if (use_gw_front .or. use_gw_front_igw) then
     ! Output for gravity waves from frontogenesis.
     call gw_drag_cam_front_diag_init(use_gw_front, use_gw_front_igw, &
                                      pgwv, gw_dc, pgwv_long, gw_dc_long)
  end if

  ! Moving mountain terms.
  if (use_gw_movmtn_pbl) then
     call gw_drag_cam_movmtn_diag_init(use_gw_movmtn_pbl, gw_drag_file_mm, movmtn_psteer, movmtn_plaunch, movmtn_source)
  end if

  ! Gravity wave convective terms.
  if (use_gw_convect_dp .or. use_gw_convect_sh) then
     call gw_drag_cam_beres_diag_init(use_gw_convect_dp, use_gw_convect_sh, gw_drag_file, gw_drag_file_sh, pgwv,gw_dc)
  end if

  ! Ridge meso-beta/gamma terms.
  if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     call gw_drag_cam_rdg_diag_init(use_gw_rdg_beta,use_gw_rdg_gamma)
  end if
  call addfld ('EKGW' ,(/ 'ilev' /), 'A','M2/S', 'Effective Kzz due to diffusion by gravity waves')

  if (history_waccm) then
     call add_default('EKGW', 1, ' ')
  end if

  ! Total tendencies for all gravity wave parameterizations.
  call addfld ('UTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', 'Total U tendency due to gravity wave drag')
  call addfld ('VTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', 'Total V tendency due to gravity wave drag')
  call register_vector_field('UTGW_TOTAL', 'VTGW_TOTAL')
  call addfld ('TTGW', (/ 'lev' /), 'A', 'K s-1', 'T tendency - gravity wave drag')
  call addfld('QTGW',(/ 'lev' /), 'A','kg/kg/s', 'Q tendency - gravity wave drag')
  call addfld('CLDLIQTGW',(/ 'lev' /), 'A','kg/kg/s', 'CLDLIQ tendency - gravity wave drag')
  call addfld('CLDICETGW',(/ 'lev' /), 'A','kg/kg/s', 'CLDICE tendency - gravity wave drag')

  if ( history_budget ) then
     call add_default('TTGW', history_budget_histfile_num, ' ')
     call add_default('QTGW', history_budget_histfile_num, ' ')
     call add_default('CLDLIQTGW', history_budget_histfile_num, ' ')
     call add_default('CLDICETGW', history_budget_histfile_num, ' ')
  end if

  ! Get indices to actually output the above.
  call cnst_get_ind("CLDLIQ", ixcldliq)
  call cnst_get_ind("CLDICE", ixcldice)

end subroutine gw_drag_cam_init

!==========================================================================

subroutine gw_drag_cam_front_diag_init(use_gw_front, use_gw_front_igw, pgwv, gw_dc, pgwv_long, gw_dc_long)
  use cam_history, only: addfld, add_default, register_vector_field

  logical, intent(in)  :: use_gw_front
  logical, intent(in)  :: use_gw_front_igw
  integer, intent(in)  :: pgwv
  real(r8), intent(in) :: gw_dc
  integer, intent(in)  :: pgwv_long
  real(r8), intent(in) :: gw_dc_long

  logical :: history_waccm

  ! Used to decide whether temperature tendencies should be output
  call phys_getopts(history_waccm_out = history_waccm)

  call addfld ('FRONTGF', (/ 'lev' /), 'A', 'K^2/M^2/S', 'Frontogenesis function at gws src level')
  call addfld ('FRONTGFA', (/ 'lev' /), 'A', 'K^2/M^2/S', 'Frontogenesis function at gws src level')

  if (history_waccm) then
     call add_default('FRONTGF', 1, ' ')
     call add_default('FRONTGFA', 1, ' ')
  end if

  if (use_gw_front) then
     ! Output for gravity waves from frontogenesis (C&M scheme)
     call gw_spec_addflds(pgwv, gw_dc, prefix=cm_pf, scheme="C&M", &
          history_defaults=history_waccm)
  end if

  if (use_gw_front_igw) then
     ! Output for inertial gravity waves from frontogenesis (C&M IGW scheme)
     call gw_spec_addflds(pgwv_long, gw_dc_long, prefix=cm_igw_pf, scheme="C&M IGW", &
          history_defaults=history_waccm)
  end if

end subroutine gw_drag_cam_front_diag_init

subroutine gw_drag_cam_oro_diag_init()
  use cam_history, only: addfld, add_default, register_vector_field
  use cam_history_support, only: horiz_only

  logical :: history_amwg, history_waccm

  ! Used to decide whether temperature tendencies should be output
  call phys_getopts(history_amwg_out = history_amwg, &
                    history_waccm_out = history_waccm)

  ! Declare history variables for orographic term
  call addfld ('TAUAORO',    (/ 'ilev' /), 'I','N m-2',  &
      'Total stress from original OGW scheme')
  call addfld ('TTGWORO',    (/ 'lev' /), 'A','K s-1',  &
      'T tendency - orographic gravity wave drag')
  call addfld ('TTGWSDFORO', (/ 'lev' /), 'A','K s-1',  &
      'T tendency - orographic gravity wave, diffusion.')
  call addfld ('TTGWSKEORO', (/ 'lev' /), 'A','K s-1',  &
      'T tendency - orographic gravity wave, breaking KE.')
  call addfld ('UTGWORO',    (/ 'lev' /), 'A','m s-2', &
      'U tendency - orographic gravity wave drag')
  call addfld ('VTGWORO',    (/ 'lev' /), 'A','m s-2', &
      'V tendency - orographic gravity wave drag')
  call register_vector_field('UTGWORO', 'VTGWORO')
  call addfld ('TAUGWX',     horiz_only,  'A','N m-2', &
      'Zonal gravity wave surface stress')
  call addfld ('TAUGWY',     horiz_only,  'A','N m-2', &
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

end subroutine gw_drag_cam_oro_diag_init

subroutine gw_drag_cam_rdg_diag_init(use_gw_rdg_beta,use_gw_rdg_gamma)
  use cam_history, only: addfld, add_default, register_vector_field
  use cam_history_support, only: horiz_only

  logical, intent(in) ::  use_gw_rdg_beta
  logical, intent(in) ::  use_gw_rdg_gamma

  character(len=1) :: cn
  integer          :: i
  logical :: history_waccm

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_waccm_out = history_waccm)

  if (use_gw_rdg_beta) then

     call addfld('UEGW',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'Zonal wind profile-entry to GW ' )
     call addfld('VEGW',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'Merdional wind profile-entry to GW ' )
     call register_vector_field('UEGW','VEGW')
     call addfld('TEGW',  (/ 'lev' /) , 'A'  ,'K' ,  &
          'Temperature profile-entry to GW ' )
     call addfld('ZEGW',  (/ 'ilev' /) , 'A'  ,'m' ,  &
          'interface geopotential heights in GW code ' )
     call addfld('ZMGW',  (/ 'lev' /) , 'A'  ,'m' ,  &
          'midlayer geopotential heights in GW code ' )

     call addfld('TAUARDGBETAY' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('TAUARDGBETAX' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGBETAX','TAUARDGBETAY')

     if (history_waccm) then
        call add_default('TAUARDGBETAX', 1, ' ')
        call add_default('TAUARDGBETAY  ', 1, ' ')
     end if
  end if

  if (use_gw_rdg_gamma) then

     call addfld ('TAUARDGGAMMAY' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld ('TAUARDGGAMMAX' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call register_vector_field('TAUARDGGAMMAX','TAUARDGGAMMAY')
     call addfld ('TAURDGGMX',     horiz_only,  'A','N m-2', &
          'Zonal gravity wave surface stress')
     call addfld ('TAURDGGMY',     horiz_only,  'A','N m-2', &
          'Meridional gravity wave surface stress')
     call register_vector_field('TAURDGGMX','TAURDGGMY')
     call addfld ('UTRDGGM' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'U wind tendency from ridge 6     ')
     call addfld ('VTRDGGM' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'V wind tendency from ridge 6     ')
     call register_vector_field('UTRDGGM','VTRDGGM')

  end if
end subroutine gw_drag_cam_rdg_diag_init
!==========================================================================

subroutine gw_drag_cam_beres_diag_init(use_gw_convect_dp,use_gw_convect_sh, gw_drag_file, gw_drag_file_sh, ngwv, dc)

  use physconst,       only: pi
  use ref_pres,        only: pref_edge
  use cam_history, only: addfld, add_default, register_vector_field
  use cam_history_support, only: horiz_only

  logical, intent(in)            :: use_gw_convect_dp,use_gw_convect_sh
  integer, intent(in)            :: ngwv
  real(r8), intent(in)           :: dc
  character(len=cl), intent(in)  :: gw_drag_file
  character(len=cl), intent(in)  :: gw_drag_file_sh
  ! output variables of interest in WACCM runs
  logical :: history_waccm
  character(len=512)              :: errmsg
  integer                         :: errflg

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_waccm_out = history_waccm)

  if (use_gw_convect_dp) then
     ! Output for gravity waves from the Beres scheme (deep).
     call gw_spec_addflds(ngwv, dc, prefix=beres_dp_pf, scheme="Beres (deep)", &
          history_defaults=history_waccm)

     call addfld ('NETDT',(/ 'lev' /), 'A','K s-1', &
          'Net heating rate')
     call addfld ('MAXQ0',horiz_only  ,  'A','K day-1', &
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
     ! Output for gravity waves from the Beres scheme (shallow).
     call gw_spec_addflds(ngwv, dc, prefix=beres_sh_pf, scheme="Beres (shallow)", &
          history_defaults=history_waccm)

     call addfld ('SNETDT',(/ 'lev' /), 'A','K s-1', &
          'Net heating rate')
     call addfld ('SMAXQ0',horiz_only  ,  'A','K day-1', &
          'Max column heating rate')
     call addfld ('SHDEPTH',horiz_only,    'A','km', &
          'Heating Depth')

     if (history_waccm) then
        call add_default('SNETDT   ', 1, ' ')
        call add_default('SHDEPTH  ', 1, ' ')
        call add_default('SMAXQ0   ', 1, ' ')
     end if
  end if

end subroutine gw_drag_cam_beres_diag_init

!==============================================================
subroutine gw_drag_cam_movmtn_diag_init(use_gw_movmtn_pbl, file_name, psteer, plaunch, source)

  use ioFileMod, only: getfil
  use ref_pres,   only: pref_edge
  use physics_buffer,   only: pbuf_get_index
  use cam_history, only: addfld, add_default, register_vector_field
  use cam_history_support, only: horiz_only

  character(len=*), intent(in) :: file_name
  logical, intent(in) :: use_gw_movmtn_pbl
  real(r8), intent(in) :: psteer
  real(r8), intent(in) :: plaunch
  integer, intent(in)  :: source

  ! PIO variable ids and error code.
  integer :: mfccid, uhid, hdid, stat

  ! Full path to gw_drag_file.
  character(len=cl) :: file_path
  character(len=cl) :: msg

  character(len=512):: errormsg
  integer           :: errorflg

  if (use_gw_movmtn_pbl) then
     call addfld ('VORT4GW', (/ 'lev' /), 'A', 's-1', &
          'Vorticity')
     call addfld ('GWUT_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - ubm component')
     call addfld ('UTGW_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - u component')
     call addfld ('VTGW_MOVMTN',(/ 'lev' /), 'I','m s-2', &
          'Mov Mtn dragforce - v component')
     call addfld('TAU_MOVMTN', (/ 'ilev' /), 'I', 'N m-2', &
          'Moving Mountain momentum flux profile')
     call addfld('U_MOVMTN_IN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint zonal input wind')
     call addfld('V_MOVMTN_IN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint meridional input wind')
     call addfld('UBI_MOVMTN', (/ 'ilev' /), 'I', 'm s-1', &
          'Moving Mountain - interface wind in direction of wave')
     call addfld('UBM_MOVMTN', (/ 'lev' /), 'I', 'm s-1', &
          'Moving Mountain - midpoint wind in direction of wave')
     call addfld ('HDEPTH_MOVMTN',horiz_only,'I','km', &
          'Heating Depth')
     call addfld ('UCELL_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - Source-level X-wind')
     call addfld ('VCELL_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - Source-level Y-wind')
     call addfld ('CS_MOVMTN',horiz_only,'I','m s-1', &
          'Gravity Wave Moving Mountain - phase speed in direction of wave')
     call addfld ('STEER_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - steering level for movmtn GW')
     call addfld ('SRC_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - launch level for movmtn GW')
     call addfld ('TND_LEVEL_MOVMTN',horiz_only,'I','1', &
          'Gravity Wave Moving Mountain - tendency lowest level for movmtn GW')
     call addfld ('NETDT_MOVMTN',(/ 'lev' /),'I','K s-1', &
          'Gravity Wave Moving Mountain - Net heating rate')
     call addfld ('TTEND_CLUBB',(/ 'lev' /),'A','K s-1', &
          'Gravity Wave Moving Mountain - CLUBB Net heating rate')
     call addfld ('THLP2_CLUBB_GW',(/ 'ilev' /),'A','K+2', &
          'Gravity Wave Moving Mountain - THLP variance from CLUBB to GW')
     call addfld ('WPTHLP_CLUBB_GW',(/ 'ilev' /),'A','Km s-2', &
          'Gravity Wave Moving Mountain - WPTHLP from CLUBB to GW')
     call addfld ('UPWP_CLUBB_GW',(/ 'ilev' /),'A','m+2 s-2', &
          'Gravity Wave Moving Mountain - X-momflux from CLUBB to GW')
     call addfld ('VPWP_CLUBB_GW',(/ 'ilev' /),'A','m+2 s-2', &
          'Gravity Wave Moving Mountain - Y-momflux from CLUBB to GW')
     call addfld ('XPWP_SRC_MOVMTN',horiz_only,'I','m+2 s-2', &
          'Gravity Wave Moving Mountain - flux source for moving mtn')
  end if

end subroutine gw_drag_cam_movmtn_diag_init
!==========================================================================


!==========================================================================

! In fact, we'd usually expect PIO errors to abort the run before you can
! even check the error code. But just in case, use this little assert.
subroutine handle_pio_error(stat, message)
  use pio, only: pio_noerr
  integer, intent(in) :: stat
  character(len=*) :: message

  call shr_assert(stat == pio_noerr, &
       "PIO error:"//trim(message)// &
       shr_errMsg(__FILE__, __LINE__))

end subroutine handle_pio_error

!==========================================================================

subroutine gw_drag_cam_tend(state, pbuf, dt, ptend, cam_in, flx_heat)
  !-----------------------------------------------------------------------
  ! Interface for multiple gravity wave drag parameterization.
  !-----------------------------------------------------------------------

  use physics_types,   only: physics_state_copy, set_dry_to_wet
  use constituents,    only: cnst_type
  use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
  use camsrfexch,      only: cam_in_t
  ! Location-dependent cpair
  use air_composition, only: cpairv
  use physconst,       only: pi, rair, gravit, cpair
  use time_manager,    only: get_step_size
  use coords_1d,       only: Coords1D

  ! CCPPized subroutines
  use gravity_wave_drag_interstitials, only: gravity_wave_drag_prepare_profiles_run
  use gravity_wave_drag_moving_mountain,  only: gravity_wave_drag_moving_mountain_run
  use gravity_wave_drag_convection,      only: gravity_wave_drag_convection_deep_run
  use gravity_wave_drag_convection,      only: gravity_wave_drag_convection_shallow_run
  use gravity_wave_drag_frontogenesis,   only: gravity_wave_drag_frontogenesis_run
  use gravity_wave_drag_frontogenesis,   only: gravity_wave_drag_frontogenesis_inertial_run
  use gravity_wave_drag_orographic,      only: gravity_wave_drag_orographic_run
  use gravity_wave_drag_ridge,           only: gravity_wave_drag_ridge_beta_run
  use gravity_wave_drag_ridge,           only: gravity_wave_drag_ridge_gamma_run
  use gravity_wave_drag_interstitials,   only: gravity_wave_drag_prepare_profiles_timestep_final

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
  integer :: istat

  integer :: i, k                   ! loop indices

  integer :: m                      ! dummy integers
  real(r8) :: qtgw(state%ncol,pver,pcnst) ! constituents tendencies

  ! gravity wave wind tendency for each wave
  real(r8), allocatable :: gwut(:,:,:)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(state%ncol,pver)
  real(r8) :: dttke(state%ncol,pver)

  ! Wave phase speeds for each column
  real(r8), allocatable :: phase_speeds(:,:)

  ! Efficiency for a gravity wave source.
  real(r8) :: effgw(state%ncol)

  ! pbuf fields
  ! Molecular diffusivity
  real(r8), pointer :: kvt_in(:,:)
  real(r8) :: kvtt(state%ncol,pver+1)
  real(r8) :: kvt_gw(state%ncol,pver+1)
  real(r8) :: sgharr(state%ncol)

  ! Frontogenesis
  real(r8), pointer :: frontgf(:,:)
  real(r8), pointer :: frontga(:,:)

  ! Standard deviation of orography.
  real(r8), pointer :: sgh(:)

  ! Indices of gravity wave source and lowest level where wind tendencies
  ! are allowed.
  integer :: src_level(state%ncol)
  integer :: tend_level(state%ncol)

  ! Convective source heating depth.
  ! heating depth
  real(r8) :: hdepth(state%ncol)
  ! maximum heating rate
  real(r8) :: maxq0(state%ncol)

  ! Temporaries for output from individual gw schemes.
  real(r8) :: ubi(state%ncol, pver+1)! projection of wind at interfaces
  real(r8) :: ubm(state%ncol, pver)  ! projection of wind at midpoints
  real(r8) :: xv(state%ncol)        ! unit vector of source wind (x)
  real(r8) :: yv(state%ncol)        ! unit vector of source wind (y)
  real(r8) :: ttgw(state%ncol, pver) ! temperature tendency
  real(r8) :: utgw(state%ncol, pver) ! zonal wind tendency
  real(r8) :: vtgw(state%ncol, pver) ! meridional wind tendency

  ! Scale sgh to account for landfrac.
  real(r8) :: sgh_scaled(state%ncol)

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

  ! Temporaries for profiles from interstitial
  type(Coords1D) :: p
  real(r8) :: rhoi(pcols, pver+1)
  real(r8) :: nm(pcols, pver)
  real(r8) :: ni(pcols, pver+1)

  ! Temporaries for diagnostic
  real(r8) :: taurx(state%ncol, pver + 1) ! wave stress in zonal direction
  real(r8) :: taury(state%ncol, pver + 1) ! wave stress in meridional direction
  real(r8) :: tauardgx(state%ncol, pver + 1) ! ridge based momentum profile
  real(r8) :: tauardgy(state%ncol, pver + 1) ! ridge based momentum profile
  real(r8) :: utrdg(state%ncol, pver) ! tendency accummulators
  real(r8) :: vtrdg(state%ncol, pver)
  real(r8) :: ttrdg(state%ncol, pver)
  real(r8) :: tau0x(pcols), tau0y(pcols)
  real(r8) :: taua(pcols, pver+1), tau0(pcols, pver+1)
  real(r8) :: gwut0(pcols, pver)
  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8) :: taucd_west (pcols, pver+1)
  real(r8) :: taucd_east (pcols, pver+1)
  real(r8) :: taucd_south(pcols, pver+1)
  real(r8) :: taucd_north(pcols, pver+1)
  real(r8) :: utend1(pcols, pver)
  real(r8) :: utend2(pcols, pver)
  real(r8) :: utend3(pcols, pver)
  real(r8) :: utend4(pcols, pver)
  real(r8) :: utend5(pcols, pver)

  character(len=64)               :: scheme_name
  character(len=512)              :: errmsg
  integer                         :: errflg

  real(r8) :: ttend_dp_arr(pcols, pver)
  real(r8) :: ttend_sh_arr(pcols, pver)

  !------------------------------------------------------------------------
  ! Make local copy of input state.
  call physics_state_copy(state, state1)

  lchnk = state1%lchnk
  ncol  = state1%ncol
  lq(:) = .true.
  call physics_ptend_init(ptend, state1%psetcols, "Gravity wave drag", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  if (do_molec_diff) then
     !--------------------------------------------------------
     ! Initialize and calculate local molecular diffusivity
     !--------------------------------------------------------
     call pbuf_get_field(pbuf, kvt_idx, kvt_in)  ! kvtt(1:pcols,1:pver+1)
     kvtt = kvt_in(:ncol,:)
  else
     kvtt = 0._r8
  end if

  ! Totals that accumulate over different sources, initialize over pcols.
  egwdffi_tot(:,:) = 0._r8
  flx_heat(:) = 0._r8

  if (use_gw_front .or. use_gw_front_igw) then
     ! Get frontogenesis physics buffer fields set by dynamics.
     call pbuf_get_field(pbuf, frontgf_idx, frontgf)
     call pbuf_get_field(pbuf, frontga_idx, frontga)

     ! Output for diagnostics.
     call outfld ('FRONTGF', frontgf, pcols, lchnk)
     call outfld ('FRONTGFA', frontga, pcols, lchnk)
  end if

  if(use_gw_movmtn_pbl .or. use_gw_convect_dp) then
    ! Set up heating
    call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
    ttend_dp_arr(:ncol, :pver) = ttend_dp(:ncol, :pver)
  else
    ttend_dp_arr(:,:) = 0._r8
  endif

  if (use_gw_movmtn_pbl) then
     ! New couplings from CLUBB
     call pbuf_get_field(pbuf, ttend_clubb_idx, ttend_clubb)
     call pbuf_get_field(pbuf, thlp2_clubb_gw_idx, thlp2_clubb_gw)
     call pbuf_get_field(pbuf, wpthlp_clubb_gw_idx, wpthlp_clubb_gw)
     call pbuf_get_field(pbuf, upwp_clubb_gw_idx, upwp_clubb_gw)
     call pbuf_get_field(pbuf, vpwp_clubb_gw_idx, vpwp_clubb_gw)
     call pbuf_get_field(pbuf, vort4gw_idx, vort4gw)
  end if

  if (use_gw_convect_sh) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_sh_idx, ttend_sh)
     ttend_sh_arr(:ncol,:pver) = ttend_sh(:ncol,:pver)
  else
     ttend_sh_arr(:,:) = 0._r8
  endif

  if (use_gw_oro) then
     call pbuf_get_field(pbuf, sgh_idx, sgh)
     sgharr(:ncol) = sgh(:ncol)
  else
     sgharr(:) = 0._r8
  end if

  rhoi(:,:) = 0._r8
  nm(:,:) = 0._r8
  ni(:,:) = 0._r8

  ! Call the CCPPized subroutine to compute necessary profiles for gravity wave drag parameterizations.
  call gravity_wave_drag_prepare_profiles_run( &
       ncol    = ncol, &
       pver    = pver, &
       cpair   = cpair, &
       rair    = rair, &
       gravit  = gravit, &
       pint    = state1%pint(:ncol,:pver+1), &
       t       = state1%t(:ncol,:pver), &
       do_molec_diff = do_molec_diff, &
       nbot_molec = nbot_molec, &
       cpairv  = cpairv(:ncol,:pver,lchnk), &
       kvt     = kvtt(:ncol,:pver+1), &
       ! below output
       p       = p, &
       rhoi    = rhoi(:ncol,:pver+1), &
       nm      = nm(:ncol,:pver), &
       ni      = ni(:ncol,:pver+1), &
       egwdffi = egwdffi_tot(:ncol,:pver+1), &
       kvt_gw  = kvt_gw(:ncol,:pver+1), &
       flx_heat= flx_heat(:ncol), &
       scheme_name = scheme_name, &
       errmsg  = errmsg, &
       errflg  = errflg)

  ! Call the CCPPized subroutine for the moving mountain gravity waves.
  if (use_gw_movmtn_pbl) then
    call outfld('U_MOVMTN_IN', state1%u, ncol, lchnk)
    call outfld('V_MOVMTN_IN', state1%v, ncol, lchnk)

    tau0(:,:) = 0._r8
    gwut0(:,:) = 0._r8

    call gravity_wave_drag_moving_mountain_run( &
      ncol                = ncol, &
      pver                = pver, &
      pcnst               = pcnst, &
      gravit              = gravit, &
      rair                = rair, &
      dt                  = dt, &
      p                   = p, &
      vramp               = vramp, &
      state_u             = state1%u(:ncol,:), &
      state_v             = state1%v(:ncol,:), &
      state_t             = state1%t(:ncol,:), &
      state_q             = state1%q(:ncol,:,:), &
      dse                 = state1%s(:ncol,:), &
      pint                = state1%pint(:ncol,:), &
      piln                = state1%lnpint(:ncol,:), &
      rhoi                = rhoi(:ncol,:), &
      nm                  = nm(:ncol,:), &
      ni                  = ni(:ncol,:), &
      kvt_gw              = kvt_gw(:ncol,:pver+1), &
      ttend_dp            = ttend_dp_arr(:ncol,:), &
      ttend_clubb         = ttend_clubb(:ncol,:), &
      upwp_clubb          = upwp_clubb_gw(:ncol,:), &
      vpwp_clubb          = vpwp_clubb_gw(:ncol,:), &
      vorticity           = vort4gw(:ncol,:), &           ! only in SE dycore.
      zm                  = state1%zm(:ncol,:), &
      alpha_gw_movmtn     = alpha_gw_movmtn, &
      effgw_movmtn_pbl    = effgw_movmtn_pbl, &
      gw_apply_tndmax     = gw_apply_tndmax, &
      use_gw_movmtn_pbl   = use_gw_movmtn_pbl, &
      ! Input/output arguments
      q_tend              = ptend%q(:ncol,:pver,:pcnst), &
      u_tend              = ptend%u(:ncol,:pver), &
      v_tend              = ptend%v(:ncol,:pver), &
      s_tend              = ptend%s(:ncol,:pver), &
      ! Output arguments
      src_level           = src_level(:ncol), &
      tend_level          = tend_level(:ncol), &
      ubm                 = ubm(:ncol,:pver), &
      ubi                 = ubi(:ncol,:pver+1), &
      xv                  = xv(:ncol), &
      yv                  = yv(:ncol), &
      hdepth              = hdepth(:ncol), &
      utgw                = utgw(:ncol,:pver), &
      vtgw                = vtgw(:ncol,:pver), &
      ttgw                = ttgw(:ncol,:pver), &
      qtgw                = qtgw(:ncol,:pver,:pcnst), &
      egwdffi_tot         = egwdffi_tot(:ncol,:pver+1), &
      dttdf               = dttdf(:ncol,:pver), &
      dttke               = dttke(:ncol,:pver), &
      tau0                = tau0(:ncol,:pver+1), &
      gwut0               = gwut0(:ncol,:pver), &
      errmsg              = errmsg, &
      errflg              = errflg)

    !-------------------------------------------------------------
    ! gw_movmtn_src returns wave-relative wind profiles ubm,ubi
    ! and unit vector components describing direction of wavevector
    ! and application of wave-drag force. I believe correct setting
    ! for c is c=0, since it is incorporated in ubm and (xv,yv)
    !--------------------------------------------------------------

    call outfld('SRC_LEVEL_MOVMTN', real(src_level,r8), ncol, lchnk)
    call outfld('TND_LEVEL_MOVMTN', real(tend_level,r8), ncol, lchnk)
    call outfld('UBI_MOVMTN', ubi, ncol, lchnk)
    call outfld('UBM_MOVMTN', ubm, ncol, lchnk)

    call outfld('TAU_MOVMTN', tau0, ncol, lchnk)
    call outfld('GWUT_MOVMTN', gwut0, ncol, lchnk)
    call outfld('VTGW_MOVMTN', vtgw, ncol, lchnk)
    call outfld('UTGW_MOVMTN', utgw, ncol, lchnk)
    call outfld('HDEPTH_MOVMTN', hdepth/1000._r8, ncol, lchnk)
    call outfld('NETDT_MOVMTN', ttend_dp, pcols, lchnk)
    call outfld('TTEND_CLUBB', ttend_clubb, pcols, lchnk)
    call outfld('THLP2_CLUBB_GW', thlp2_clubb_gw, pcols, lchnk)
    call outfld('WPTHLP_CLUBB_GW', wpthlp_clubb_gw, pcols, lchnk)
    call outfld('UPWP_CLUBB_GW', upwp_clubb_gw, pcols, lchnk)
    call outfld('VPWP_CLUBB_GW', vpwp_clubb_gw, pcols, lchnk)
    call outfld('VORT4GW', vort4gw, pcols, lchnk)
  end if

  ! Convective gravity waves (Beres scheme, deep).
  if (use_gw_convect_dp) then
    taucd_west(:,:) = 0._r8
    taucd_east(:,:) = 0._r8
    taucd_south(:,:) = 0._r8
    taucd_north(:,:) = 0._r8

    call gravity_wave_drag_convection_deep_run( &
          ncol            = ncol, &
          pver            = pver, &
          pcnst           = pcnst, &
          dt              = dt, &
          p               = p, &
          vramp           = vramp, &
          pi              = pi, &
          cpair           = cpair, &
          effgw_beres_dp  = effgw_beres_dp, &
          gw_apply_tndmax = gw_apply_tndmax, &
          u               = state1%u(:ncol,:pver), &
          v               = state1%v(:ncol,:pver), &
          t               = state1%t(:ncol,:pver), &
          q               = state1%q(:ncol,:pver,:), &
          dse             = state1%s(:ncol,:pver), &
          piln            = state1%lnpint(:ncol,:pver+1), &
          rhoi            = rhoi(:ncol,:pver+1), &
          nm              = nm(:ncol,:pver), &
          ni              = ni(:ncol,:pver+1), &
          kvt_gw          = kvt_gw(:ncol,:pver+1), &
          ttend_dp        = ttend_dp_arr(:ncol,:pver), &
          zm              = state1%zm(:ncol,:pver), &
          lat             = state1%lat(:ncol), &
          tend_q          = ptend%q(:ncol,:pver,:), &
          tend_u          = ptend%u(:ncol,:pver), &
          tend_v          = ptend%v(:ncol,:pver), &
          tend_s          = ptend%s(:ncol,:pver), &
          flx_heat        = flx_heat(:ncol), &
          src_level       = src_level(:ncol), &
          tend_level      = tend_level(:ncol), &
          ubm             = ubm(:ncol,:pver), &
          ubi             = ubi(:ncol,:pver+1), &
          xv              = xv(:ncol), &
          yv              = yv(:ncol), &
          hdepth          = hdepth(:ncol), &
          maxq0           = maxq0(:ncol), &
          utgw            = utgw(:ncol,:pver), &
          vtgw            = vtgw(:ncol,:pver), &
          ttgw            = ttgw(:ncol,:pver), &
          qtgw            = qtgw(:ncol,:pver,:), &
          egwdffi_tot     = egwdffi_tot(:ncol,:pver+1), &
          dttdf           = dttdf(:ncol,:pver), &
          dttke           = dttke(:ncol,:pver), &
          taucd_west      = taucd_west(:ncol,:pver+1), &
          taucd_east      = taucd_east(:ncol,:pver+1), &
          taucd_south     = taucd_south(:ncol,:pver+1), &
          taucd_north     = taucd_north(:ncol,:pver+1), &
          errmsg          = errmsg, &
          errflg          = errflg)

    call outfld(trim(beres_dp_pf) // 'TTGWSDF', dttdf / cpair, ncol, lchnk)
    call outfld(trim(beres_dp_pf) // 'TTGWSKE', dttke / cpair, ncol, lchnk)

    ! Simple output fields written to history file.
    ! Total wind tendencies.
    call outfld (trim(beres_dp_pf)//'UTGWSPEC', utgw , ncol, lchnk)
    call outfld (trim(beres_dp_pf)//'VTGWSPEC', vtgw , ncol, lchnk)
    call outfld (trim(beres_dp_pf)//'TTGWSPEC', ttgw , ncol, lchnk)

    ! Tau in each direction.
    call outfld (trim(beres_dp_pf)//'TAUE', taucd_east,  pcols, lchnk)
    call outfld (trim(beres_dp_pf)//'TAUW', taucd_west,  pcols, lchnk)
    call outfld (trim(beres_dp_pf)//'TAUN', taucd_north, pcols, lchnk)
    call outfld (trim(beres_dp_pf)//'TAUS', taucd_south, pcols, lchnk)
    call outfld (trim(beres_dp_pf)//'TAUNET', taucd_east + taucd_west, pcols, lchnk)

    ! Diagnostic outputs (convert hdepth to km).
    call outfld('NETDT', ttend_dp, pcols, lchnk)
    call outfld('HDEPTH', hdepth/1000._r8, ncol, lchnk)
    call outfld('MAXQ0', maxq0, ncol, lchnk)
  end if

  ! Convective gravity waves (Beres scheme, shallow).
  if (use_gw_convect_sh) then
    call gravity_wave_drag_convection_shallow_run( &
          ncol            = ncol, &
          pver            = pver, &
          pcnst           = pcnst, &
          dt              = dt, &
          p               = p, &
          vramp           = vramp, &
          pi              = pi, &
          cpair           = cpair, &
          effgw_beres_sh  = effgw_beres_sh, &
          gw_apply_tndmax = gw_apply_tndmax, &
          u               = state1%u(:ncol,:pver), &
          v               = state1%v(:ncol,:pver), &
          t               = state1%t(:ncol,:pver), &
          q               = state1%q(:ncol,:pver,:), &
          dse             = state1%s(:ncol,:pver), &
          piln            = state1%lnpint(:ncol,:pver+1), &
          rhoi            = rhoi(:ncol,:pver+1), &
          nm              = nm(:ncol,:pver), &
          ni              = ni(:ncol,:pver+1), &
          kvt_gw          = kvt_gw(:ncol,:pver+1), &
          ttend_sh        = ttend_sh_arr(:ncol,:pver), &
          zm              = state1%zm(:ncol,:pver), &
          lat             = state1%lat(:ncol), &
          tend_q          = ptend%q(:ncol,:pver,:), &
          tend_u          = ptend%u(:ncol,:pver), &
          tend_v          = ptend%v(:ncol,:pver), &
          tend_s          = ptend%s(:ncol,:pver), &
          flx_heat        = flx_heat(:ncol), &
          src_level       = src_level(:ncol), &
          tend_level      = tend_level(:ncol), &
          ubm             = ubm(:ncol,:pver), &
          ubi             = ubi(:ncol,:pver+1), &
          xv              = xv(:ncol), &
          yv              = yv(:ncol), &
          hdepth          = hdepth(:ncol), &
          maxq0           = maxq0(:ncol), &
          utgw            = utgw(:ncol,:pver), &
          vtgw            = vtgw(:ncol,:pver), &
          ttgw            = ttgw(:ncol,:pver), &
          qtgw            = qtgw(:ncol,:pver,:), &
          egwdffi_tot     = egwdffi_tot(:ncol,:pver+1), &
          dttdf           = dttdf(:ncol,:pver), &
          dttke           = dttke(:ncol,:pver), &
          errmsg          = errmsg, &
          errflg          = errflg)

    call outfld(trim(beres_sh_pf) // 'TTGWSDF', dttdf / cpair, ncol, lchnk)
    call outfld(trim(beres_sh_pf) // 'TTGWSKE', dttke / cpair, ncol, lchnk)

    ! Simple output fields written to history file.
    ! Total wind tendencies.
    call outfld (trim(beres_sh_pf)//'UTGWSPEC', utgw , ncol, lchnk)
    call outfld (trim(beres_sh_pf)//'VTGWSPEC', vtgw , ncol, lchnk)
    call outfld (trim(beres_sh_pf)//'TTGWSPEC', ttgw , ncol, lchnk)

    ! Diagnostic outputs (convert hdepth to km).
    call outfld('SNETDT',  ttend_sh, pcols, lchnk)
    call outfld('SHDEPTH', hdepth/1000._r8, ncol, lchnk)
    call outfld('SMAXQ0',  maxq0, ncol, lchnk)
  end if

  ! Call the CCPPized subroutine
  if(use_gw_front) then
    taucd_west(:,:) = 0._r8
    taucd_east(:,:) = 0._r8
    taucd_south(:,:) = 0._r8
    taucd_north(:,:) = 0._r8
    utend1(:,:) = 0._r8
    utend2(:,:) = 0._r8
    utend3(:,:) = 0._r8
    utend4(:,:) = 0._r8
    utend5(:,:) = 0._r8

    call gravity_wave_drag_frontogenesis_run( &
         ncol             = ncol, &
         pver             = pver, &
         pcnst            = pcnst, &
         dt               = dt, &
         cpair            = cpair, &
         p                = p, &
         vramp            = vramp, &
         piln             = state1%lnpint(:ncol,:), &
         rhoi             = rhoi(:ncol,:pver+1), &
         nm               = nm(:ncol,:pver), &
         ni               = ni(:ncol,:pver+1), &
         effgw_cm         = effgw_cm, &
         gw_polar_taper   = gw_polar_taper, &
         gw_apply_tndmax  = gw_apply_tndmax, &
         lat              = state1%lat(:ncol), &
         u                = state1%u(:ncol,:), &
         v                = state1%v(:ncol,:), &
         t                = state1%t(:ncol,:), &
         q                = state1%q(:ncol,:,:), &
         dse              = state1%s(:ncol,:), &
         frontgf          = frontgf(:ncol,:pver), &
         kvt_gw           = kvt_gw(:ncol,:pver+1), &
         ! below input/output (accummulated tendencies)
         tend_q           = ptend%q(:ncol,:,:), &
         tend_u           = ptend%u(:ncol,:), &
         tend_v           = ptend%v(:ncol,:), &
         tend_s           = ptend%s(:ncol,:), &
         egwdffi_tot      = egwdffi_tot(:ncol,:pver+1), &
         ! below output
         src_level        = src_level(:ncol), &
         tend_level       = tend_level(:ncol), &
         ubm              = ubm(:ncol,:pver), &
         ubi              = ubi(:ncol,:pver+1), &
         xv               = xv(:ncol), &
         yv               = yv(:ncol), &
         utgw             = utgw(:ncol,:pver), &
         vtgw             = vtgw(:ncol,:pver), &
         ttgw             = ttgw(:ncol,:pver), &
         qtgw             = qtgw(:ncol,:pver,:), &
         dttdf            = dttdf(:ncol,:pver), &
         dttke            = dttke(:ncol,:pver), &
         taucd_west       = taucd_west(:ncol,:pver+1), &
         taucd_east       = taucd_east(:ncol,:pver+1), &
         taucd_south      = taucd_south(:ncol,:pver+1), &
         taucd_north      = taucd_north(:ncol,:pver+1), &
         utend1           = utend1(:ncol,:pver), &
         utend2           = utend2(:ncol,:pver), &
         utend3           = utend3(:ncol,:pver), &
         utend4           = utend4(:ncol,:pver), &
         utend5           = utend5(:ncol,:pver), &
         flx_heat         = flx_heat(:ncol), &
         errmsg           = errmsg, &
         errflg           = errflg)

    if(errflg /= 0) then
      call endrun("gravity_wave_drag_frontogenesis_run: " // errmsg)
    endif

    ! Output wind tendencies binned by phase speed.
    call outfld(trim(cm_pf)//'UTEND1', utend1, pcols, lchnk)
    call outfld(trim(cm_pf)//'UTEND2', utend2, pcols, lchnk)
    call outfld(trim(cm_pf)//'UTEND3', utend3, pcols, lchnk)
    call outfld(trim(cm_pf)//'UTEND4', utend4, pcols, lchnk)
    call outfld(trim(cm_pf)//'UTEND5', utend5, pcols, lchnk)

    ! Output temperature tendencies due to diffusion and from kinetic energy.
    call outfld(trim(cm_pf)//'TTGWSDF', dttdf / cpair, ncol, lchnk)
    call outfld(trim(cm_pf)//'TTGWSKE', dttke / cpair, ncol, lchnk)

    ! Simple output fields written to history file.
    ! Total wind tendencies.
    call outfld (trim(cm_pf)//'UTGWSPEC', utgw , ncol, lchnk)
    call outfld (trim(cm_pf)//'VTGWSPEC', vtgw , ncol, lchnk)
    call outfld (trim(cm_pf)//'TTGWSPEC', ttgw , ncol, lchnk)

    ! Tau in each direction.
    call outfld (trim(cm_pf)//'TAUE', taucd_east,  pcols, lchnk)
    call outfld (trim(cm_pf)//'TAUW', taucd_west,  pcols, lchnk)
    call outfld (trim(cm_pf)//'TAUN', taucd_north, pcols, lchnk)
    call outfld (trim(cm_pf)//'TAUS', taucd_south, pcols, lchnk)
    call outfld (trim(cm_pf)//'TAUNET', taucd_east + taucd_west, pcols, lchnk)
  endif

  if(use_gw_front_igw) then
    call gravity_wave_drag_frontogenesis_inertial_run( &
         ncol             = ncol, &
         pver             = pver, &
         pcnst            = pcnst, &
         dt               = dt, &
         cpair            = cpair, &
         p                = p, &
         vramp            = vramp, &
         piln             = state1%lnpint(:ncol,:), &
         rhoi             = rhoi(:ncol,:pver+1), &
         nm               = nm(:ncol,:pver), &
         ni               = ni(:ncol,:pver+1), &
         effgw_cm_igw     = effgw_cm_igw, &
         gw_polar_taper   = gw_polar_taper, &
         gw_apply_tndmax  = gw_apply_tndmax, &
         lat              = state1%lat(:ncol), &
         u                = state1%u(:ncol,:), &
         v                = state1%v(:ncol,:), &
         t                = state1%t(:ncol,:), &
         q                = state1%q(:ncol,:,:), &
         dse              = state1%s(:ncol,:), &
         frontgf          = frontgf(:ncol,:pver), &
         kvt_gw           = kvt_gw(:ncol,:pver+1), &
         ! below input/output (accummulated tendencies)
         tend_q           = ptend%q(:ncol,:,:), &
         tend_u           = ptend%u(:ncol,:), &
         tend_v           = ptend%v(:ncol,:), &
         tend_s           = ptend%s(:ncol,:), &
         egwdffi_tot      = egwdffi_tot(:ncol,:pver+1), &
         ! below output
         src_level        = src_level(:ncol), &
         tend_level       = tend_level(:ncol), &
         ubm              = ubm(:ncol,:pver), &
         ubi              = ubi(:ncol,:pver+1), &
         xv               = xv(:ncol), &
         yv               = yv(:ncol), &
         utgw             = utgw(:ncol,:pver), &
         vtgw             = vtgw(:ncol,:pver), &
         ttgw             = ttgw(:ncol,:pver), &
         qtgw             = qtgw(:ncol,:pver,:), &
         dttdf            = dttdf(:ncol,:pver), &
         dttke            = dttke(:ncol,:pver), &
         flx_heat         = flx_heat(:ncol), &
         errmsg           = errmsg, &
         errflg           = errflg)

    if(errflg /= 0) then
      call endrun("gravity_wave_drag_frontogenesis_inertial_run: " // errmsg)
    endif
  endif

  if(use_gw_oro) then
    tau0x(:) = 0._r8
    tau0y(:) = 0._r8
    taua(:,:) = 0._r8

    call gravity_wave_drag_orographic_run( &
      ncol              = ncol, &
      pver              = pver, &
      pcnst             = pcnst, &
      dt                = dt, &
      rair              = rair, &
      p                 = p, &
      vramp             = vramp(:), &
      piln              = state1%lnpint(:ncol,:pver+1), &
      rhoi              = rhoi(:ncol,:pver+1), &
      nm                = nm(:ncol,:pver), &
      ni                = ni(:ncol,:pver+1), &
      effgw_oro         = effgw_oro, &
      gw_lndscl_sgh     = gw_lndscl_sgh, &
      gw_oro_south_fac  = gw_oro_south_fac, &
      gw_apply_tndmax   = gw_apply_tndmax, &
      landfrac          = cam_in%landfrac(:ncol), &
      lat               = state1%lat(:ncol), &
      u                 = state1%u(:ncol,:pver), &
      v                 = state1%v(:ncol,:pver), &
      t                 = state1%t(:ncol,:pver), &
      q                 = state1%q(:ncol,:pver,:pcnst), &
      dse               = state1%s(:ncol,:pver), &
      sgh               = sgharr(:ncol), &
      zm                = state1%zm(:ncol,:pver), &
      kvt_gw            = kvt_gw(:ncol,:pver+1), &
      cpairv            = cpairv(:ncol,:pver,lchnk), &
      tend_q            = ptend%q(:ncol,:pver,:pcnst), &
      tend_u            = ptend%u(:ncol,:pver), &
      tend_v            = ptend%v(:ncol,:pver), &
      tend_s            = ptend%s(:ncol,:pver), &
      src_level         = src_level(:ncol), &
      tend_level        = tend_level(:ncol), &
      ubm               = ubm(:ncol,:pver), &
      ubi               = ubi(:ncol,:pver+1), &
      xv                = xv(:ncol), &
      yv                = yv(:ncol), &
      utgw              = utgw(:ncol,:pver), &
      vtgw              = vtgw(:ncol,:pver), &
      ttgw              = ttgw(:ncol,:pver), &
      qtgw              = qtgw(:ncol,:pver,:), &
      dttdf             = dttdf(:ncol,:pver), &
      dttke             = dttke(:ncol,:pver), &
      egwdffi_tot       = egwdffi_tot(:ncol,:pver+1), &
      tau0x             = tau0x(:ncol), &
      tau0y             = tau0y(:ncol), &
      taua              = taua(:ncol,:pver+1), &
      flx_heat          = flx_heat(:ncol), &
      errmsg            = errmsg, &
      errflg            = errflg)
    if(errflg /= 0) then
      call endrun("gravity_wave_drag_orographic_run: " // errmsg)
    endif

    ! Write output fields to history file
    call outfld('TAUAORO', taua,  ncol, lchnk)
    call outfld('UTGWORO', utgw,  ncol, lchnk)
    call outfld('VTGWORO', vtgw,  ncol, lchnk)
    call outfld('TTGWORO', ttgw,  ncol, lchnk)
    call outfld('TTGWSDFORO', dttdf / cpair,  ncol, lchnk)
    call outfld('TTGWSKEORO', dttke / cpair,  ncol, lchnk)
    call outfld('TAUGWX', tau0x, ncol, lchnk)
    call outfld('TAUGWY', tau0y, ncol, lchnk)
  endif

  if (use_gw_rdg_beta) then
    ! Save state at top of routine
    ! Useful for unit testing checks
    call outfld('UEGW', state1%u ,  ncol, lchnk)
    call outfld('VEGW', state1%v ,  ncol, lchnk)
    call outfld('TEGW', state1%t ,  ncol, lchnk)
    call outfld('ZEGW', state1%zi , ncol, lchnk)
    call outfld('ZMGW', state1%zm , ncol, lchnk)

    call gravity_wave_drag_ridge_beta_run( &
      ncol                    = ncol, &
      pver                    = pver, &
      pcnst                   = pcnst, &
      dt                      = dt, &
      pi                      = pi, &
      cpair                   = cpair, &
      rair                    = rair, &
      vramp                   = vramp, &
      p                       = p, &
      n_rdg_beta              = n_rdg_beta, &
      u                       = state1%u(:ncol,:), &
      v                       = state1%v(:ncol,:), &
      t                       = state1%t(:ncol,:), &
      q                       = state1%q(:ncol,:,:), &
      dse                     = state1%s(:ncol,:), &
      piln                    = state1%lnpint(:ncol,:), &
      zm                      = state1%zm(:ncol,:), &
      zi                      = state1%zi(:ncol,:pver+1), &
      nm                      = nm(:ncol,:), &
      ni                      = ni(:ncol,:pver+1), &
      rhoi                    = rhoi(:ncol,:pver+1), &
      kvt_gw                  = kvt_gw(:ncol,:pver+1), &
      use_gw_rdg_resid        = use_gw_rdg_resid, &
      effgw_rdg_resid         = effgw_rdg_resid, &
      effgw_rdg_beta          = effgw_rdg_beta, &
      effgw_rdg_beta_max      = effgw_rdg_beta_max, &
      rdg_beta_cd_llb         = rdg_beta_cd_llb, &
      trpd_leewv_rdg_beta     = trpd_leewv_rdg_beta, &
      ! Input data for beta waves - stored for all chunks
      gbxar                   = rdg_gbxar(:ncol,lchnk), &
      isovar                  = rdg_isovar(:ncol,lchnk), &
      isowgt                  = rdg_isowgt(:ncol,lchnk), &
      hwdth                   = rdg_hwdth(:ncol,:,lchnk), &
      clngt                   = rdg_clngt(:ncol,:,lchnk), &
      mxdis                   = rdg_mxdis(:ncol,:,lchnk), &
      anixy                   = rdg_anixy(:ncol,:,lchnk), &
      angll                   = rdg_angll(:ncol,:,lchnk), &
      ! Input/output arguments
      s_tend                  = ptend%s(:ncol,:pver), &
      q_tend                  = ptend%q(:ncol,:pver,:pcnst), &
      u_tend                  = ptend%u(:ncol,:pver), &
      v_tend                  = ptend%v(:ncol,:pver), &
      ! Output arguments
      flx_heat                = flx_heat(:ncol), &
      taurx                   = taurx(:ncol,:pver+1), &
      taury                   = taury(:ncol,:pver+1), &
      tauardgx                = tauardgx(:ncol,:pver+1), &
      tauardgy                = tauardgy(:ncol,:pver+1), &
      utrdg                   = utrdg(:ncol,:pver), &
      vtrdg                   = vtrdg(:ncol,:pver), &
      ttrdg                   = ttrdg(:ncol,:pver), &
      errmsg                  = errmsg, &
      errflg                  = errflg)

     call outfld('TAUGWX', taurx(:,pver+1), ncol, lchnk)
     call outfld('TAUGWY', taury(:,pver+1), ncol, lchnk)
     call outfld('UTGWORO', utrdg,  ncol, lchnk)
     call outfld('VTGWORO', vtrdg,  ncol, lchnk)
     call outfld('TTGWORO', ttrdg,  ncol, lchnk)

     call outfld('TAUARDGBETAX', tauardgx, ncol, lchnk)
     call outfld('TAUARDGBETAY', tauardgy, ncol, lchnk)

  end if

  if (use_gw_rdg_gamma) then
    ! Save state at top of routine
    ! Useful for unit testing checks
    call outfld('UEGW', state1%u ,  ncol, lchnk)
    call outfld('VEGW', state1%v ,  ncol, lchnk)
    call outfld('TEGW', state1%t ,  ncol, lchnk)
    call outfld('ZEGW', state1%zi , ncol, lchnk)
    call outfld('ZMGW', state1%zm , ncol, lchnk)

    call gravity_wave_drag_ridge_gamma_run( &
      ncol                    = ncol, &
      pver                    = pver, &
      pcnst                   = pcnst, &
      dt                      = dt, &
      pi                      = pi, &
      cpair                   = cpair, &
      rair                    = rair, &
      vramp                   = vramp, &
      p                       = p, &
      n_rdg_gamma             = n_rdg_gamma, &
      u                       = state1%u(:ncol,:), &
      v                       = state1%v(:ncol,:), &
      t                       = state1%t(:ncol,:), &
      q                       = state1%q(:ncol,:,:), &
      dse                     = state1%s(:ncol,:), &
      piln                    = state1%lnpint(:ncol,:), &
      zm                      = state1%zm(:ncol,:), &
      zi                      = state1%zi(:ncol,:pver+1), &
      nm                      = nm(:ncol,:), &
      ni                      = ni(:ncol,:pver+1), &
      rhoi                    = rhoi(:ncol,:pver+1), &
      kvt_gw                  = kvt_gw(:ncol,:pver+1), &
      effgw_rdg_resid         = effgw_rdg_resid, &
      use_gw_rdg_resid        = use_gw_rdg_resid, &
      effgw_rdg_gamma         = effgw_rdg_gamma, &
      effgw_rdg_gamma_max     = effgw_rdg_gamma_max, &
      rdg_gamma_cd_llb        = rdg_gamma_cd_llb, &
      trpd_leewv_rdg_gamma    = trpd_leewv_rdg_gamma, &
      ! Input data for beta/gamma waves - stored for all chunks
      gbxarg                  = rdg_gbxarg(:ncol,lchnk), &
      hwdthg                  = rdg_hwdthg(:ncol,:,lchnk), &
      clngtg                  = rdg_clngtg(:ncol,:,lchnk), &
      mxdisg                  = rdg_mxdisg(:ncol,:,lchnk), &
      anixyg                  = rdg_anixyg(:ncol,:,lchnk), &
      angllg                  = rdg_angllg(:ncol,:,lchnk), &
      ! Input/output arguments
      s_tend                  = ptend%s(:ncol,:pver), &
      q_tend                  = ptend%q(:ncol,:pver,:pcnst), &
      u_tend                  = ptend%u(:ncol,:pver), &
      v_tend                  = ptend%v(:ncol,:pver), &
      ! Output arguments
      flx_heat                = flx_heat(:ncol), &
      taurx                   = taurx(:ncol,:pver+1), &
      taury                   = taury(:ncol,:pver+1), &
      tauardgx                = tauardgx(:ncol,:pver+1), &
      tauardgy                = tauardgy(:ncol,:pver+1), &
      utrdg                   = utrdg(:ncol,:pver), &
      vtrdg                   = vtrdg(:ncol,:pver), &
      ttrdg                   = ttrdg(:ncol,:pver), &
      errmsg                  = errmsg, &
      errflg                  = errflg)

     call outfld('TAURDGGMX', taurx(:,pver+1), ncol, lchnk)
     call outfld('TAURDGGMY', taury(:,pver+1), ncol, lchnk)
     call outfld('UTRDGGM', utrdg,  ncol, lchnk)
     call outfld('VTRDGGM', vtrdg,  ncol, lchnk)
     call outfld('TTGWORO', ttrdg,  ncol, lchnk)

     call outfld('TAUARDGGAMMAX', tauardgx, ncol, lchnk)
     call outfld('TAUARDGGAMMAY', tauardgy, ncol, lchnk)
  end if

  ! Call the CCPPized subroutine to clean up
  call gravity_wave_drag_prepare_profiles_timestep_final( &
    p      = p, &
    errmsg = errmsg, &
    errflg = errflg)

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


end subroutine gw_drag_cam_tend


!==========================================================================
!==========================================================================

! Add all history fields for a gravity wave spectrum source.
subroutine gw_spec_addflds(ngwv, dc, prefix, scheme, history_defaults)
  use cam_history, only: addfld, add_default, register_vector_field

  !------------------------------Arguments--------------------------------

  ! One character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Gravity wave scheme name prepended to output field descriptions.
  character(len=*), intent(in) :: scheme
  integer, intent(in)               :: ngwv
  real(r8), intent(in)              :: dc
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

  real(r8), allocatable :: cref(:)

  !-----------------------------------------------------------------------

  ! Uniform phase speed reference grid.
  allocate(cref(-ngwv:ngwv))
  cref = [( dc * l, l = -ngwv, ngwv )]

  ! Overall wind tendencies.
  call addfld (trim(prefix)//'UTGWSPEC',(/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency - gravity wave spectrum')
  call addfld (trim(prefix)//'VTGWSPEC',(/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' V tendency - gravity wave spectrum')
  call register_vector_field(trim(prefix)//'UTGWSPEC',trim(prefix)//'VTGWSPEC')

  call addfld (trim(prefix)//'TTGWSPEC',(/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' T tendency - gravity wave spectrum')

  ! Wind tendencies broken across five spectral bins.
  call addfld (trim(prefix)//'UTEND1',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency   c < -40')
  call addfld (trim(prefix)//'UTEND2',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency  -40 < c < -15')
  call addfld (trim(prefix)//'UTEND3',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency  -15 < c <  15')
  call addfld (trim(prefix)//'UTEND4',  (/ 'lev' /), 'A','m s-2', &
       trim(scheme)//' U tendency   15 < c <  40')
  call addfld (trim(prefix)//'UTEND5',  (/ 'lev' /), 'A','m s-2', &
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
  call addfld (trim(prefix)//'TTGWSDF' , (/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' t tendency - diffusion term')
  call addfld (trim(prefix)//'TTGWSKE' , (/ 'lev' /), 'A','K s-1', &
       trim(scheme)//' t tendency - kinetic energy conversion term')

  ! Gravity wave source spectra by wave number.
  do l=-ngwv,ngwv
     ! String containing reference speed.
     write (fnum,fmt='(f7.2)') cref(l)

     dumc1x = tau_fld_name(l, prefix, x_not_y=.true.)
     dumc1y = tau_fld_name(l, prefix, x_not_y=.false.)
     dumc2 = trim(scheme)//" tau at c= "//trim(fnum)//" m s-1"
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
  deallocate(cref)

end subroutine gw_spec_addflds

!==========================================================================

! Outputs for spectral waves.
subroutine gw_spec_outflds(prefix, lchnk, ncol, ngwv, phase_speeds, dc, u, v, xv, yv, &
     gwut, dttdf, dttke, tau, utgw, vtgw, ttgw, taucd)

  use gw_common, only: west, east, south, north

  ! One-character prefix prepended to output fields.
  character(len=1), intent(in) :: prefix
  ! Chunk and number of columns in the chunk.
  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  integer, intent(in) :: ngwv
  ! Wave speeds.
!jt  type(GWBand), intent(in) :: band
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: phase_speeds(ncol,-ngwv:ngwv)
  real(r8), intent(in) :: dc
  ! Winds at cell midpoints.
  real(r8), intent(in) :: u(ncol,pver)
  real(r8), intent(in) :: v(ncol,pver)
  ! Unit vector in the direction of wind at source level.
  real(r8), intent(in) :: xv(ncol)
  real(r8), intent(in) :: yv(ncol)
  ! Wind tendency for each wave.
  real(r8), intent(in) :: gwut(ncol,pver,-ngwv:ngwv)
  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8) :: dttdf(ncol,pver)
  real(r8) :: dttke(ncol,pver)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(ncol,-ngwv:ngwv,pver)
  ! Zonal and meridional total wind tendencies.
  real(r8), intent(in) :: utgw(ncol,pver)
  real(r8), intent(in) :: vtgw(ncol,pver)
  ! Temperature tendencies.
  real(r8), intent(in) :: ttgw(ncol,pver)
  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8), intent(in) :: taucd(ncol,pver+1,4)

  ! Indices
  integer :: i, k, l
  integer :: ix(ncol, -ngwv:ngwv), iy(ncol, -ngwv:ngwv)
  integer :: iu(ncol), iv(ncol)

  ! Zonal wind tendency, broken up into five bins.
  real(r8) :: utb(ncol, pver, 5)
  ! Definition of the bin boundaries.
  real(r8), parameter :: bounds(4) = (/ -40._r8, -15._r8, &
       15._r8, 40._r8 /)

  ! Momentum flux in the four cardinal directions.
  real(r8) :: mf(ncol, pver, 4)

  ! Wave stress in zonal/meridional direction
  real(r8) :: taux(ncol,-ngwv:ngwv,pver)
  real(r8) :: tauy(ncol,-ngwv:ngwv,pver)

  ! Temporaries for output
  real(r8) :: dummyx(ncol,pver)
  real(r8) :: dummyy(ncol,pver)
  ! Variable names
  character(len=10) :: dumc1x, dumc1y


  ! Accumulate wind tendencies binned according to phase speed.
  utb = 0._r8

  ! Find which output bin the phase speed corresponds to.
  ix = find_bin(phase_speeds)

  ! Put the wind tendency in that bin.
  do l = -ngwv, ngwv
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

  ! Project phase_speeds, and convert each component to a wavenumber index.
  ! These are mappings from the wavenumber index of tau to those of taux
  ! and tauy, respectively.
  do l=-ngwv,ngwv
     ix(:,l) = c_to_l(phase_speeds(:,l)*xv)
     iy(:,l) = c_to_l(phase_speeds(:,l)*yv)
  end do

  ! Find projection of tau.
  do k = 1, pver
     do l = -ngwv,ngwv
        do i = 1, ncol
           taux(i,ix(i,l),k) = taux(i,ix(i,l),k) &
                + abs(tau(i,l,k)*xv(i))
           tauy(i,iy(i,l),k) = tauy(i,iy(i,l),k) &
                + abs(tau(i,l,k)*yv(i))
        end do
     end do
  end do

  do l=-ngwv,ngwv

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
     do l = -ngwv, ngwv

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

    l = min( max(int(c/dc),-ngwv), ngwv )

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

end module gw_drag_cam
