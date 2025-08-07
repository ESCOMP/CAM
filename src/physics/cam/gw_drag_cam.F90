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

  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = unset_r8   ! critical froude number squared

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

  use gw_drag,        only: gw_drag_init
  use gravity_wave_drag_top_taper, only: gravity_wave_drag_top_taper_init

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

  call gw_drag_init( &
       iulog_in                     = iulog, &
       ktop_in                      = ktop, &
       masterproc_in                = masterproc, &
       ncol                         = pcols, &
       pver                         = pver, &
       gravit_in                    = gravit, &
       rair_in                      = rair, &
       pi_in                        = pi, &
       fcrit2_in                    = fcrit2, &
       rearth_in                    = rearth, &
       pref_edge                    = pref_edge, &
       pgwv_nl                      = pgwv, &
       gw_dc_nl                     = gw_dc, &
       pgwv_long_nl                 = pgwv_long, &
       gw_dc_long_nl                = gw_dc_long, &
       tau_0_ubc_nl                 = tau_0_ubc, &
       effgw_beres_dp_nl            = effgw_beres_dp, &
       effgw_beres_sh_nl            = effgw_beres_sh, &
       effgw_cm_nl                  = effgw_cm, &
       effgw_cm_igw_nl              = effgw_cm_igw, &
       effgw_oro_nl                 = effgw_oro, &
       frontgfc_nl                  = frontgfc, &
       gw_drag_file_nl              = gw_drag_file_loc, &
       gw_drag_file_sh_nl           = gw_drag_file_sh_loc, &
       gw_drag_file_mm_nl           = gw_drag_file_mm_loc, &
       taubgnd_nl                   = taubgnd, &
       taubgnd_igw_nl               = taubgnd_igw, &
       gw_polar_taper_nl            = gw_polar_taper, &
       use_gw_rdg_beta_nl           = use_gw_rdg_beta, &
       n_rdg_beta_nl                = n_rdg_beta, &
       effgw_rdg_beta_nl            = effgw_rdg_beta, &
       effgw_rdg_beta_max_nl        = effgw_rdg_beta_max, &
       rdg_beta_cd_llb_nl           = rdg_beta_cd_llb, &
       trpd_leewv_rdg_beta_nl       = trpd_leewv_rdg_beta, &
       use_gw_rdg_gamma_nl          = use_gw_rdg_gamma, &
       n_rdg_gamma_nl               = n_rdg_gamma, &
       effgw_rdg_gamma_nl           = effgw_rdg_gamma, &
       effgw_rdg_gamma_max_nl       = effgw_rdg_gamma_max, &
       rdg_gamma_cd_llb_nl          = rdg_gamma_cd_llb, &
       trpd_leewv_rdg_gamma_nl      = trpd_leewv_rdg_gamma, &
       gw_oro_south_fac_nl          = gw_oro_south_fac, &
       gw_limit_tau_without_eff_nl  = gw_limit_tau_without_eff, &
       gw_lndscl_sgh_nl             = gw_lndscl_sgh, &
       gw_prndl_nl                  = gw_prndl, &
       gw_apply_tndmax_nl           = gw_apply_tndmax, &
       gw_qbo_hdepth_scaling_nl     = gw_qbo_hdepth_scaling, &
       gw_top_taper_nl              = gw_top_taper, &
       front_gaussian_width_nl      = front_gaussian_width, &
       alpha_gw_movmtn_nl           = alpha_gw_movmtn, &
       use_gw_rdg_resid_in          = use_gw_rdg_resid, &
       effgw_rdg_resid_in           = effgw_rdg_resid, &
       effgw_movmtn_pbl_in          = effgw_movmtn_pbl, &
       movmtn_source_in             = movmtn_source, &
       movmtn_psteer_in             = movmtn_psteer, &
       movmtn_plaunch_in            = movmtn_plaunch, &
       gw_rdg_do_divstream_nl       = gw_rdg_do_divstream, &
       gw_rdg_do_smooth_regimes_nl  = gw_rdg_do_smooth_regimes, &
       gw_rdg_do_adjust_tauoro_nl   = gw_rdg_do_adjust_tauoro, &
       gw_rdg_do_backward_compat_nl = gw_rdg_do_backward_compat, &
       gw_rdg_do_vdiff_nl           = gw_rdg_do_vdiff, &
       gw_rdg_C_BetaMax_DS_nl       = gw_rdg_C_BetaMax_DS, &
       gw_rdg_C_GammaMax_nl         = gw_rdg_C_GammaMax, &
       gw_rdg_Frx0_nl               = gw_rdg_Frx0, &
       gw_rdg_Frx1_nl               = gw_rdg_Frx1, &
       gw_rdg_C_BetaMax_SM_nl       = gw_rdg_C_BetaMax_SM, &
       gw_rdg_Fr_c_nl               = gw_rdg_Fr_c, &
       gw_rdg_orohmin_nl            = gw_rdg_orohmin, &
       gw_rdg_orovmin_nl            = gw_rdg_orovmin, &
       gw_rdg_orostratmin_nl        = gw_rdg_orostratmin, &
       gw_rdg_orom2min_nl           = gw_rdg_orom2min, &
       use_gw_oro_in                = use_gw_oro, &
       use_gw_front_in              = use_gw_front, &
       use_gw_front_igw_in          = use_gw_front_igw, &
       use_gw_convect_dp_in         = use_gw_convect_dp, &
       use_gw_convect_sh_in         = use_gw_convect_sh, &
       use_simple_phys_in           = use_simple_phys, &
       use_gw_movmtn_pbl_in         = use_gw_movmtn_pbl, &
       do_molec_diff_in             = do_molec_diff, &
       nbot_molec_in                = nbot_molec, &
       wavelength_mid_in            = wavelength_mid, &
       wavelength_long_in           = wavelength_long, &
       errmsg                       = errmsg, &
       errflg                       = errflg)

  if(errflg /= 0) then
    call endrun("gw_init: " // errmsg)
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
    call endrun("gravity_wave_drag_top_taper_init: " // errmsg)
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


  !------- diagnostics code

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm, &
       history_amwg_out   = history_amwg  )

  if ( use_gw_oro ) then

     if (effgw_oro == unset_r8) then
        call endrun("gw_init: Orographic gravity waves enabled, &
             &but effgw_oro was not set.")
     end if
  end if

  if (use_gw_oro .or. use_gw_rdg_beta .or. use_gw_rdg_gamma) then

     sgh_idx = pbuf_get_index('SGH')

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

  end if

  ! ========= gw front initialization! ==========================
  if (use_gw_front .or. use_gw_front_igw) then

     frontgf_idx = pbuf_get_index('FRONTGF')
     frontga_idx = pbuf_get_index('FRONTGA')

     call shr_assert(unset_r8 /= frontgfc, &
          "gw_init: Frontogenesis enabled, but frontgfc was &
          & not set!"// &
          shr_errMsg(__FILE__, __LINE__))

     call addfld ('FRONTGF', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')
     call addfld ('FRONTGFA', (/ 'lev' /), 'A', 'K^2/M^2/S', &
          'Frontogenesis function at gws src level')

     if (history_waccm) then
        call add_default('FRONTGF', 1, ' ')
        call add_default('FRONTGFA', 1, ' ')
     end if

     if (use_gw_front) then
        if (masterproc) then
           write(iulog,*) 'gw_init: gw spectrum taubgnd, ', &
                'effgw_cm = ',taubgnd, effgw_cm
           write(iulog,*) ' '
        end if
!jt ******** Add this to a gw_drag_cam_front_diag_init subroutine
        ! Output for gravity waves from frontogenesis.
!jt ******** create bands on CAM side to pass to spec addflds
        !jt call gw_spec_addflds(prefix=cm_pf, scheme="C&M", band=band_mid, &
        !jt     history_defaults=history_waccm)
     end if
  end if
  ! ========= Moving Mountain initialization! ==========================
  if (use_gw_movmtn_pbl) then
     ! get pbuf indices for CLUBB couplings
     ttend_clubb_idx     = pbuf_get_index('TTEND_CLUBB')
     thlp2_clubb_gw_idx  = pbuf_get_index('THLP2_CLUBB_GW')
     upwp_clubb_gw_idx   = pbuf_get_index('UPWP_CLUBB_GW')
     vpwp_clubb_gw_idx   = pbuf_get_index('VPWP_CLUBB_GW')
     wpthlp_clubb_gw_idx = pbuf_get_index('WPTHLP_CLUBB_GW')
     vort4gw_idx         = pbuf_get_index('VORT4GW')

     if (masterproc) then
        write (iulog,*) 'Moving Mountain development code call init_movmtn'
     end if

     call gw_drag_cam_movmtn_diag_init(use_gw_movmtn_pbl, gw_drag_file_mm,movmtn_psteer, movmtn_plaunch, movmtn_source)
  end if
  ! ========= Convect diagnostics init! ==========================

  if (use_gw_convect_dp .or. use_gw_convect_sh) then

     if (use_gw_convect_dp) ttend_dp_idx    = pbuf_get_index('TTEND_DP')

     if (use_gw_convect_sh) ttend_sh_idx    = pbuf_get_index('TTEND_SH')

     call gw_drag_cam_beres_diag_init(use_gw_convect_dp,use_gw_convect_sh, gw_drag_file, gw_drag_file_sh, pgwv,gw_dc)
  end if
  ! ========= Convect Beta/Gamma initialization! ==========================
  if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     call gw_drag_cam_rdg_diag_init(use_gw_rdg_beta,use_gw_rdg_gamma)
  end if
  call addfld ('EKGW' ,(/ 'ilev' /), 'A','M2/S', &
       'Effective Kzz due to diffusion by gravity waves')

  if (history_waccm) then
     call add_default('EKGW', 1, ' ')
  end if

  call addfld ('UTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', &
       'Total U tendency due to gravity wave drag')
  call addfld ('VTGW_TOTAL',    (/ 'lev' /), 'A','m s-2', &
       'Total V tendency due to gravity wave drag')
  call register_vector_field('UTGW_TOTAL', 'VTGW_TOTAL')

  ! Total temperature tendency output.
  call addfld ('TTGW', (/ 'lev' /), 'A', 'K s-1',  &
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


end subroutine gw_drag_cam_init

!==========================================================================

subroutine gw_drag_cam_rdg_diag_init(use_gw_rdg_beta,use_gw_rdg_gamma)
  use cam_history, only: addfld, add_default, register_vector_field

  logical, intent(in) ::  use_gw_rdg_beta
  logical, intent(in) ::  use_gw_rdg_gamma

  character(len=1) :: cn
  integer          :: i
  logical :: history_waccm
  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  ! Used to decide whether temperature tendencies should be output.
  call phys_getopts( history_waccm_out = history_waccm)

  if (use_gw_rdg_beta) then

     call addfld ('WBR_HT1',     horiz_only,  'I','m', &
          'Wave breaking height for DSW')
     call addfld ('TLB_HT1',     horiz_only,  'I','m', &
          'Form drag layer height')
     call addfld ('BWV_HT1',     horiz_only,  'I','m', &
          'Bottom of freely-propagating OGW regime')
     call addfld ('TAUDSW1',     horiz_only,  'I','Nm-2', &
          'DSW enhanced drag')
     call addfld ('TAUORO1',     horiz_only,  'I','Nm-2', &
          'lower BC on propagating wave stress')
     call addfld ('UBMSRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level on-ridge wind')
     call addfld ('USRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level Zonal wind')
     call addfld ('VSRC1',     horiz_only,  'I','ms-1', &
          'below-peak-level Meridional wind')
     call addfld ('NSRC1',     horiz_only,  'I','s-1', &
          'below-peak-level stratification')
     call addfld ('MXDIS1',     horiz_only,  'I','m', &
          'Ridge/obstacle height')
     call addfld ('ANGLL1',     horiz_only,  'I','degrees', &
          'orientation clockwise w/resp north-south')
     call addfld ('ANIXY1',     horiz_only,  'I','1', &
          'Ridge quality')
     call addfld ('HWDTH1',     horiz_only,  'I','km', &
          'Ridge width')
     call addfld ('CLNGT1',     horiz_only,  'I','km', &
          'Ridge length')
     call addfld ('GBXAR1',     horiz_only,  'I','km+2', &
          'grid box area')

     call addfld ('Fr1_DIAG',     horiz_only,  'I','1', &
          'Critical Froude number for linear waves')
     call addfld ('Fr2_DIAG',     horiz_only,  'I','1', &
          'Critical Froude number for blocked flow')
     call addfld ('Frx_DIAG',     horiz_only,  'I','1', &
          'Obstacle Froude Number')

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

     call addfld('TAUM1_DIAG' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('TAU1RDGBETAM' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld('UBM1BETA',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'On-ridge wind profile           ' )
     call addfld('UBT1RDGBETA' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGBETAY' , (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGBETAX' , (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call register_vector_field('TAU'//cn//'RDGBETAX','TAU'//cn//'RDGBETAY')
        call addfld('UT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm s-1', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGBETA',    (/ 'lev' /),  'I', 'm s-1', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGBETA','VT'//cn//'RDGBETA')
     end do

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

     call addfld ('TAU1RDGGAMMAM' , (/ 'ilev' /) , 'I'  ,'N m-2' , &
          'Ridge based momentum flux profile')
     call addfld ('UBM1GAMMA',  (/ 'lev' /) , 'A'  ,'s-1' ,  &
          'On-ridge wind profile           ' )
     call addfld ('UBT1RDGGAMMA' , (/ 'lev' /) , 'I'  ,'m s-1' , &
          'On-ridge wind tendency from ridge 1     ')

     do i = 1, 6
        write(cn, '(i1)') i
        call addfld('TAU'//cn//'RDGGAMMAY', (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('TAU'//cn//'RDGGAMMAX', (/ 'ilev' /), 'I', 'N m-2', &
          'Ridge based momentum flux profile')
        call addfld('UT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm s-1', &
          'U wind tendency from ridge '//cn)
        call addfld('VT'//cn//'RDGGAMMA' , (/ 'lev' /),  'I', 'm s-1', &
          'V wind tendency from ridge '//cn)
        call register_vector_field('UT'//cn//'RDGGAMMA','VT'//cn//'RDGGAMMA')
     end do

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
  use gw_drag,         only: gw_drag_run
  use gw_rdg,          only: gw_rdg_run
  use gravity_wave_drag_interstitials, only: gravity_wave_drag_prepare_profiles_timestep_finalize

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

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8) :: taucd(state%ncol,pver+1,4)

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
  real(r8) :: sgharr(state%ncol)

  ! Frontogenesis
  real(r8), pointer :: frontgf(:,:)
  real(r8), pointer :: frontga(:,:)

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

  character(len=64)               :: scheme_name
  character(len=512)              :: errmsg
  integer                         :: errflg

  real(r8)          :: ttend_sh_arr(state%ncol,pver)

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

  ! Totals that accumulate over different sources.
  egwdffi_tot = 0._r8
  flx_heat = 0._r8

  if (use_gw_front .or. use_gw_front_igw) then
     ! Get frontogenesis physics buffer fields set by dynamics.
     call pbuf_get_field(pbuf, frontgf_idx, frontgf)
     call pbuf_get_field(pbuf, frontga_idx, frontga)
  end if

  if (use_gw_movmtn_pbl) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)

     !   New couplings from CLUBB
     call pbuf_get_field(pbuf, ttend_clubb_idx, ttend_clubb)
     call pbuf_get_field(pbuf, thlp2_clubb_gw_idx, thlp2_clubb_gw)
     call pbuf_get_field(pbuf, wpthlp_clubb_gw_idx, wpthlp_clubb_gw)
     call pbuf_get_field(pbuf, upwp_clubb_gw_idx, upwp_clubb_gw)
     call pbuf_get_field(pbuf, vpwp_clubb_gw_idx, vpwp_clubb_gw)
     call pbuf_get_field(pbuf, vort4gw_idx, vort4gw)
  end if
  if (use_gw_convect_dp) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
  end if
  if (use_gw_convect_sh) then
     ! Set up heating
     call pbuf_get_field(pbuf, ttend_sh_idx, ttend_sh)
  end if
  if (use_gw_oro) then
     call pbuf_get_field(pbuf, sgh_idx, sgh)
     sgharr=sgh
  else
     sgharr=0._r8
  end if
!!!!jt  There was a problem passing an unassociated pointer (ttend_sh_arr) it was temporarily replaced with real array.
!!!!jt  Only associated when running shallow convective gravity waves.  Fix this

  rhoi(:,:) = 0._r8
  nm(:,:) = 0._r8
  ni(:,:) = 0._r8

  ! Call the CCPPized subroutine to compute necessary profiles for gravity wave drag parameterizations.
  call gravity_wave_drag_prepare_profiles_run( &
    ncol   = ncol, &
    pver   = pver, &
    cpair  = cpair, &
    rair   = rair, &
    gravit = gravit, &
    pint   = state1%pint(:ncol,:pver+1), &
    t      = state1%t(:ncol,:pver), &
    ! below output
    p      = p, &
    rhoi   = rhoi(:ncol,:pver+1), &
    nm     = nm(:ncol,:pver), &
    ni     = ni(:ncol,:pver+1), &
    errmsg = errmsg, &
    errflg = errflg)

  ! Call the CCPPized subroutine
  call gw_drag_run( &
       ncol             = ncol, &
       pcnst            = pcnst, &
       pver             = pver, &
       dt               = dt, &
       cpair            = cpair, &
       cpairv           = cpairv(:ncol,:,lchnk), &
       pi               = pi, &
       vramp            = vramp, &
       frontgf          = frontgf(:ncol,:pver), &
       frontga          = frontga(:ncol,:pver), &
       pint             = state1%pint(:ncol,:), &
       piln             = state1%lnpint(:ncol,:), &
       pdel             = state1%pdel(:ncol,:), &
       pdeldry          = state1%pdeldry(:ncol,:), &
       zm               = state1%zm(:ncol,:), &
       zi               = state1%zi(:ncol,:), &
       lat              = state1%lat(:ncol), &
       landfrac         = cam_in%landfrac(:ncol), &
       dse              = state1%s(:ncol,:), &
       state_t          = state1%t(:ncol,:), &
       state_u          = state1%u(:ncol,:), &
       state_v          = state1%v(:ncol,:), &
       state_q          = state1%q(:ncol,:,:), &
       vorticity        = vort4gw(:ncol,:pver), &
       sgh              = sgharr(:ncol), &
       p                = p, &
       rhoi             = rhoi(:ncol,:pver+1), &
       nm               = nm(:ncol,:pver), &
       ni               = ni(:ncol,:pver+1), &
       kvtt             = kvtt(:ncol,:pver+1), &
       ttend_dp         = ttend_dp(:ncol,:pver), &
       ttend_sh         = ttend_sh_arr(:ncol,:pver), &
       ttend_clubb      = ttend_clubb(:ncol,:pver), &
       thlp2_clubb_gw   = thlp2_clubb_gw(:ncol,:pver+1), &
       wpthlp_clubb_gw  = wpthlp_clubb_gw(:ncol,:pver+1), &
       upwp_clubb_gw    = upwp_clubb_gw(:ncol,:pver+1), &
       vpwp_clubb_gw    = vpwp_clubb_gw(:ncol,:pver+1), &
       s_tend           = ptend%s(:ncol,:), &
       q_tend           = ptend%q(:ncol,:,:), &
       u_tend           = ptend%u(:ncol,:), &
       v_tend           = ptend%v(:ncol,:), &
       scheme_name      = scheme_name, &
       nbot_molec       = nbot_molec, &
       egwdffi_tot      = egwdffi_tot(:ncol,:pver+1), &
       flx_heat         = flx_heat(:ncol), &
       errmsg           = errmsg, &
       errflg           = errflg)

  if(errflg /= 0) then
    call endrun("gw_drag_run: " // errmsg)
  endif

  if (use_gw_rdg_beta .or. use_gw_rdg_gamma) then
     ! Save state at top of routine
     ! Useful for unit testing checks
      call outfld('UEGW', state1%u ,  ncol, lchnk)
      call outfld('VEGW', state1%v ,  ncol, lchnk)
      call outfld('TEGW', state1%t ,  ncol, lchnk)
      call outfld('ZEGW', state1%zi , ncol, lchnk)
      call outfld('ZMGW', state1%zm , ncol, lchnk)

      call gw_rdg_run( &
        ncol                    = ncol, &
        pver                    = pver, &
        pcnst                   = pcnst, &
        dt                      = dt, &
        pi                      = pi, &
        use_gw_rdg_beta         = use_gw_rdg_beta, &
        use_gw_rdg_gamma        = use_gw_rdg_gamma, &
        vramp                   = vramp, &
        p                       = p, &
        n_rdg_beta              = n_rdg_beta, &
        n_rdg_gamma             = n_rdg_gamma, &
        u                       = state1%u(:ncol,:), &
        v                       = state1%v(:ncol,:), &
        t                       = state1%t(:ncol,:), &
        q                       = state1%q(:ncol,:,:), &
        dse                     = state1%s(:ncol,:), &
        piln                    = state1%lnpint(:ncol,:), &
        zm                      = state1%zm(:ncol,:), &
        zi                      = state1%zi(:ncol,:), &
        nm                      = nm(:ncol,:), &
        ni                      = ni(:ncol,:), &
        rhoi                    = rhoi(:ncol,:), &
        kvtt                    = kvtt(:ncol,:), &
        effgw_rdg_resid         = effgw_rdg_resid, &
        use_gw_rdg_resid        = use_gw_rdg_resid, &
        effgw_rdg_beta          = effgw_rdg_beta, &
        effgw_rdg_beta_max      = effgw_rdg_beta_max, &
        effgw_rdg_gamma         = effgw_rdg_gamma, &
        effgw_rdg_gamma_max     = effgw_rdg_gamma_max, &
        rdg_beta_cd_llb         = rdg_beta_cd_llb, &
        trpd_leewv_rdg_beta     = trpd_leewv_rdg_beta, &
        rdg_gamma_cd_llb        = rdg_gamma_cd_llb, &
        trpd_leewv_rdg_gamma    = trpd_leewv_rdg_gamma, &
        ! Input data for beta/gamma waves - stored for all chunks
        gbxar                   = rdg_gbxar(:ncol,lchnk), &
        isovar                  = rdg_isovar(:ncol,lchnk), &
        isowgt                  = rdg_isowgt(:ncol,lchnk), &
        hwdth                   = rdg_hwdth(:ncol,:,lchnk), &
        clngt                   = rdg_clngt(:ncol,:,lchnk), &
        mxdis                   = rdg_mxdis(:ncol,:,lchnk), &
        anixy                   = rdg_anixy(:ncol,:,lchnk), &
        angll                   = rdg_angll(:ncol,:,lchnk), &
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
        errmsg                  = errmsg, &
        errflg                  = errflg)
  end if

  ! Call the CCPPized subroutine to clean up
  call gravity_wave_drag_prepare_profiles_timestep_finalize( &
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
