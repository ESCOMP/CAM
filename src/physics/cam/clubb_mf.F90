module clubb_mf

! =============================================================================== !
! Mass-flux module for use with CLUBB                                             !
! Together (CLUBB+MF) they comprise a eddy-diffusivity mass-flux approach (EDMF)  !
! =============================================================================== !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun
  use time_manager,  only: is_first_step, get_nstep
  use spmd_utils,    only: iam
  use physconst,     only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                           cpwv, cpliq, rh2o, zvir, pi

  implicit none
  private
  save

  public :: integrate_mf, &
            clubb_mf_readnl, &
            do_clubb_mf, &
            do_clubb_mf_diag, &
            clubb_mf_nup, &
            do_clubb_mf_rad, &
            clubb_mf_Lopt, &
            clubb_mf_ddalph, &
            clubb_mf_up_ndt, &
            clubb_mf_cp_ndt, &
            do_clubb_mf_cmt

  !
  ! Lopt 0 = fixed L0
  !      1 = tke_clubb L0
  !      2 = wpthlp_clubb L0
  !      3 = test plume L0
  !      4 = lel
  !      4 = cape
  !      6 = ztopm1
  !      7 = rel.hum. at 500 hPa
  !      8 = column int. rel.hum.
  integer  :: clubb_mf_Lopt    = 0
  real(r8) :: clubb_mf_a0      = 0._r8
  real(r8) :: clubb_mf_b0      = 0._r8
  real(r8) :: clubb_mf_L0      = 0._r8
  real(r8) :: clubb_mf_ent0    = 0._r8
  real(r8) :: clubb_mf_alphturb= 0._r8
  real(r8) :: clubb_mf_max_L0  = 0._r8
  real(r8) :: clubb_mf_fdd     = 0._r8
  real(r8) :: clubb_mf_ddalph  = 0._r8  
  real(r8) :: clubb_mf_ddbeta  = 0._r8
  real(r8) :: clubb_mf_pwfac   = 0._r8
  real(r8) :: clubb_mf_ddexp   = 0._r8
  integer  :: clubb_mf_up_ndt  = 1
  integer  :: clubb_mf_cp_ndt  = 1
  integer, protected :: clubb_mf_nup     = 0
  logical, protected :: do_clubb_mf = .false.
  logical, protected :: do_clubb_mf_diag = .false.
  logical, protected :: do_clubb_mf_rad = .false.
  logical, protected :: do_clubb_mf_coldpool = .false.
  logical, protected :: do_clubb_mf_ustar = .false.
  logical, protected :: do_clubb_mf_mixd = .false.
  logical, protected :: do_clubb_mf_precip = .false.
  logical, protected :: do_clubb_mf_rhtke = .false.
  logical, protected :: do_clubb_mf_cmt = .false.
  logical :: tht_tweaks = .true.
  integer :: mf_num_cin = 5

  contains

  subroutine clubb_mf_readnl(nlfile)

  ! =============================================================================== !
  ! MF namelists                                                                    !
  ! =============================================================================== !

    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8, mpi_integer, mpi_logical

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    character(len=*), parameter :: sub = 'clubb_mf_readnl'

    integer :: iunit, read_status, ierr


    namelist /clubb_mf_nl/ clubb_mf_Lopt, clubb_mf_a0, clubb_mf_b0, clubb_mf_L0, clubb_mf_ent0, clubb_mf_alphturb, &
                           clubb_mf_nup, clubb_mf_max_L0, do_clubb_mf, do_clubb_mf_diag, do_clubb_mf_precip, do_clubb_mf_rad, &
                           clubb_mf_fdd, do_clubb_mf_coldpool, clubb_mf_ddalph, clubb_mf_ddbeta, clubb_mf_pwfac, do_clubb_mf_ustar, &
                           clubb_mf_ddexp, do_clubb_mf_mixd, clubb_mf_up_ndt, clubb_mf_cp_ndt, do_clubb_mf_rhtke, do_clubb_mf_cmt

    if (masterproc) then
      open( newunit=iunit, file=trim(nlfile), status='old' )
      call find_group_name(iunit, 'clubb_mf_nl', status=read_status)
      if (read_status == 0) then
         read(iunit, clubb_mf_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun('clubb_mf_readnl: ERROR reading namelist')
         end if
      end if
      close(iunit)
    end if

    call mpi_bcast(clubb_mf_Lopt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_Lopt")
    call mpi_bcast(clubb_mf_a0,   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_a0")
    call mpi_bcast(clubb_mf_b0,   1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_b0")
    call mpi_bcast(clubb_mf_L0,   1, mpi_real8,   mstrid, mpicom, ierr) 
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_L0")
    call mpi_bcast(clubb_mf_ent0, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ent0")
    call mpi_bcast(clubb_mf_alphturb,1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_alphturb")
    call mpi_bcast(clubb_mf_nup,  1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_nup")
    call mpi_bcast(clubb_mf_max_L0,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_max_L0")
    call mpi_bcast(do_clubb_mf,      1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf")
    call mpi_bcast(do_clubb_mf_diag, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_diag")
    call mpi_bcast(do_clubb_mf_precip, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_precip")
    call mpi_bcast(do_clubb_mf_rad, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_rad")
    call mpi_bcast(clubb_mf_fdd,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_fdd")
    call mpi_bcast(do_clubb_mf_coldpool, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_coldpool")
    call mpi_bcast(clubb_mf_ddalph,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddalph")
    call mpi_bcast(clubb_mf_ddbeta,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddbeta")
    call mpi_bcast(clubb_mf_pwfac,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_pwfac")
    call mpi_bcast(clubb_mf_up_ndt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_up_ndt")
    call mpi_bcast(clubb_mf_cp_ndt, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_cp_ndt")
    call mpi_bcast(do_clubb_mf_ustar, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_ustar")
    call mpi_bcast(clubb_mf_ddexp,  1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ddexp")
    call mpi_bcast(do_clubb_mf_mixd, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_mixd")
    call mpi_bcast(do_clubb_mf_rhtke, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_rhtke")
    call mpi_bcast(do_clubb_mf_cmt, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_clubb_mf_cmt")

    if ((.not. do_clubb_mf) .and. do_clubb_mf_diag ) then
       call endrun('clubb_mf_readnl: Error - cannot turn on do_clubb_mf_diag without also turning on do_clubb_mf')
    end if
    

  end subroutine clubb_mf_readnl

  subroutine integrate_mf( nz,                                                      & ! input
                           rho_zm,  dzm,     zm,      p_zm,      iexner_zm,         & ! input
                           rho_zt,  dzt,     zt,      p_zt,      iexner_zt,         & ! input
                           u,       v,       thl,     qt,        thv,               & ! input
                                             th,      qv,        qc,                & ! input
                                             thl_zm,  qt_zm,     thv_zm,            & ! input
                                             th_zm,   qv_zm,     qc_zm,             & ! input
                           ustar,      ths,  wthl,    wqt,       pblh,              & ! input
                           wpthlp_env, tke,  tpert,  ztopm1,     rhinv,             & ! input
                           mcape,      ddcp, cbm1,                                  & ! output
                           upa,     dna,                                            & ! output
                           upw,     dnw,                                            & ! output
                           upmf,                                                    & ! output
                           upqt,    dnqt,                                           & ! output
                           upthl,   dnthl,                                          & ! output
                           upthv,   dnthv,                                          & ! output
                           upth,    dnth,                                           & ! output
                           upqc,    dnqc,                                           & ! output
                           upbuoy,                                                  & ! output
                           upent,                                                   & ! output
                           updet,                                                   & ! output
                           dry_a,   moist_a,                                        & ! output
                           dry_w,   moist_w,                                        & ! output
                           dry_qt,  moist_qt,                                       & ! output
                           dry_thl, moist_thl,                                      & ! output
                           dry_u,   moist_u,                                        & ! output
                           dry_v,   moist_v,                                        & ! output
                                    moist_qc,                                       & ! output
                           ae,                                                      &
                           ac,      aup,     adn,                                   &
                           aw,      awup,    awdn,                                  & 
                           aww,     awwup,   awwdn,                                 &
                           awthlup, awqtup,  awuup, awvup,                          & ! output
                           awthldn, awqtdn,  awudn, awvdn,                          & ! output
                           awthl,   awqt,                                           & ! output
                           awu,     awv,                                            & ! output
                           thlflxup,qtflxup, uflxup, vflxup,                        & ! output
                           thlflxdn,qtflxdn, uflxdn, vflxdn,                        & ! output
                           thlflx,  qtflx,   uflx,   vflx,                          & ! output - variables needed for solver
                           thvflx,                                                  & ! output
                           sqtup,   sthlup,                                         & ! output
                           sqtdn,   sthldn,                                         & ! output
                           sqt,     sthl,                                           & ! output - variables needed for solver
                           precc,                                                   & ! output
                           ztop,    dynamic_L0 )

  ! ================================================================================= !
  ! Mass-flux algorithm                                                               !
  !                                                                                   !
  ! Provides rtm and thl fluxes due to mass flux ensemble,                            !
  ! which are fed into the mixed explicit/implicit clubb solver as explicit terms     !
  !                                                                                   !
  ! Mass flux variables are computed on edges (i.e. momentum grid):                   !
  ! upa,upw,upqt,...                                                                  !
  ! dry_a,moist_a,dry_w,moist_w, ...                                                  !
  !                                                                                   ! 
  ! In CLUBB (unlike CAM) nlevs of momentum grid = nlevs of thermodynamic grid,       !
  ! due to a subsurface thermodynamic layer. To avoid confusion, below the variables  !  
  ! are grouped by the grid they are on.                                              !
  !                                                                                   !
  ! *note that state on the lowest thermo level is equal to state on the lowest       !
  ! momentum level due to state_zt(1) = state_zt(2), and lowest momentum level        !
  ! is a weighted combination of the lowest two thermodynamic levels.                 !
  !                                                                                   !
  ! ---------------------------------Authors----------------------------------------  !
  ! Marcin Kurowski, JPL                                                              !
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM       !
  ! Additional modifications by Adam Herrington, NCAR                                 !
  ! ================================================================================= !

     use wv_saturation,      only : qsat

     integer,  intent(in)                :: nz
     real(r8), dimension(nz), intent(in) :: u,      v,            & ! thermodynamic grid
                                            thl,    thv,          & ! thermodynamic grid
                                            th,     qv,           & ! thermodynamic grid
                                            qt,     qc,           & ! thermodynamic grid
                                            p_zt,   iexner_zt,    & ! thermodynamic grid
                                            dzt,    rho_zt,       & ! thermodynamic grid
                                            zt,                   & ! thermodynamic grid
                                            thl_zm, thv_zm,       & ! momentum grid
                                            th_zm,  qv_zm,        &
                                            qt_zm,  qc_zm,        & ! momentum grid
                                            p_zm,   iexner_zm,    & ! momentum grid
                                            dzm,    rho_zm,       & ! momentum grid
                                            zm,                   & ! momentum grid
                                            tke,    wpthlp_env      ! momentum grid

     real(r8), intent(in)                :: wthl,wqt
     real(r8), intent(in)                :: pblh,tpert
     real(r8), intent(in)                :: rhinv
     real(r8), intent(in)                :: ths,ustar
     real(r8), intent(inout)             :: ztopm1,ddcp,cbm1

     real(r8),dimension(nz,clubb_mf_nup), intent(out) :: upa,     & ! momentum grid
                                                         upw,     & ! momentum grid
                                                         upmf,    & ! momentum grid
                                                         upqt,    & ! momentum grid
                                                         upthl,   & ! momentum grid
                                                         upthv,   & ! momentum grid
                                                         upth,    & ! momentum grid
                                                         upqc,    & ! momentum grid
                                                         upbuoy,  & ! momentum grid
                                                         upent,   & ! momentum grid
                                                         updet
     !
     real(r8),dimension(nz,clubb_mf_nup), intent(out) :: dna,     & ! momentum grid
                                                         dnw,     & ! momentum grid
                                                         dnqt,    & ! momentum grid
                                                         dnthl,   & ! momentum grid
                                                         dnthv,   & ! momentum grid
                                                         dnth,    & ! momentum grid
                                                         dnqc
     !
     real(r8),dimension(nz), intent(out) :: dry_a,   moist_a,     & ! momentum grid
                                            dry_w,   moist_w,     & ! momentum grid
                                            dry_qt,  moist_qt,    & ! momentum grid
                                            dry_thl, moist_thl,   & ! momentum grid
                                            dry_u,   moist_u,     & ! momentum grid
                                            dry_v,   moist_v,     & ! momentum grid
                                                     moist_qc       ! momentum grid
     !
     real(r8),dimension(nz), intent(out) :: ae,                                &
                                            ac,      aup,     adn,             &
                                            aw,      awup,    awdn,            &
                                            aww,     awwup,  awwdn,            &
                                            awthlup, awqtup, awuup, awvup,     & ! momentum grid
                                            awthldn, awqtdn, awudn, awvdn,     & ! momentum grid
                                            awthl,   awqt,                     & ! momentum grid
                                            awu,     awv,                      & ! momentum grid
                                            thlflxup,qtflxup, uflxup, vflxup,  & ! momentum grid
                                            thlflxdn,qtflxdn, uflxdn, vflxdn,  & ! momentum grid
                                            thlflx,  qtflx,   uflx,   vflx,    & ! momentum grid
                                            thvflx,                            &
                                            sqtup,   sthlup,                   & ! thermodynamic grid
                                            sqtdn,   sthldn,                   & ! thermodynamic grid
                                            sqt,     sthl,                     & ! thermodynamic grid 
                                            precc

     real(r8), intent(out)               :: ztop,    dynamic_L0,  &
                                            mcape   
     ! =============================================================================== !
     ! INTERNAL VARIABLES
     !
     ! sums over all plumes
     real(r8), dimension(nz)              :: moist_th,   dry_th,       & 
                                             thl_env,    qt_env,       & 
                                             thv_env,                  &
                                             thvflxup,   thvflxdn,     &
                                             awthvup,    awthvdn
     !
     ! updraft properties
     real(r8), dimension(nz,clubb_mf_nup) :: upqv,     upqs,           & ! momentum grid
                                             upql,     upqi,           & ! momentum grid
                                             upu,      upv,            & ! momentum grid 
                                             uplmix,   upauto            ! momentum grid
     !
     ! downdraft properties
     real(r8), dimension(nz,clubb_mf_nup) ::           dnqs,           & ! momentum grid
                                             dnql,     dnqi,           & ! momentum grid
                                             dnu,      dnv,            & ! momentum grid 
                                             dnlmix                      ! momentum grid
     !
     ! microphyiscs terms
     real(r8), dimension(nz,clubb_mf_nup) :: supqt,    supthl,         & ! thermodynamic grid 
                                             sdnqt,    sdnthl,         & ! thermodynamic grid
                                             uprr,     dnrr                        
     !
     ! entrainment profiles
     real(r8), dimension(nz,clubb_mf_nup) :: entf,     mix               ! thermodynamic grid
     integer,  dimension(nz,clubb_mf_nup) :: enti                        ! thermodynamic grid
     ! 
     ! other variables
     integer                              :: k,i,kstart,ddtop,kcb
     integer,  dimension(clubb_mf_nup)    :: ddbot,kcbarr
     real(r8), dimension(clubb_mf_nup)    :: zcb
     real(r8)                             :: zcb_unset,       cpfac,   &
                                             wthv,   ddint,   iddcp,   &
                                             wstar,  qstar,   thvstar, & 
                                             sigmaw, sigmaqt, sigmathv,&
                                             convh,  wmin,    wmax,    & 
                                             wlv,    wtv,     wp,      & 
                                             B,                        & ! thermodynamic grid
                                             entexp, entexpu, entw,    & ! thermodynamic grid
                                             Mn,                       & ! momentum grid
                                             eturb,  det,     lmixt,   & ! thermodynamic grid
                                             qtovqs, sevap,   taum1,   & ! thermodynamic grid
                                             sqtint, sthlint, alphint, &
                                             qtmp,   betathl, betaqt,  & ! thermodynamic grid        
                                             thln,   thvn,    thn,     & ! momentum grid
                                             qtn,    qsn,              & ! momentum grid
                                             qcn,    qln,     qin,     & ! momentum grid
                                             un,     vn,      wn2,     & ! momentum grid
                                             wn,                       & ! momentum grid
                                             lmixn,   srfarea,         & ! momentum grid
                                             srfwqtu, srfwthvu,        &
                                             facqtu,  facthvu,         &
                                             zsub,    wcb,    rh_L0

!     !
!     ! cape variables
!     real(r8), dimension(nz)                :: t_zt
!     real(r8), dimension(nz-1)              :: tp,       qstp
!     !real(r8), dimension(nz-1,clubb_mf_nup) :: dmpdz
!     !real(r8), dimension(clubb_mf_nup)      :: tl,                     &
!     !                                          cape,     cin
!     !integer,  dimension(clubb_mf_nup)      :: lcl,      lel
!     real(r8), dimension(nz-1,1)            :: dmpdz
!     real(r8), dimension(1)                 :: tl,                     &
!                                               cape,     cin
!     integer,  dimension(1)                 :: lcl,      lel
!     real(r8)                               :: landfrac
!     integer                                :: kpbl,     msg,          &
!                                               lon,      mx
     !
     ! limit convective area
     logical                                :: limarea = .false.
     real(r8),parameter                     :: amax = 0.6_r8
     !
     ! buoyancy sorting variables
     logical                                :: bsort = .false.
     real(r8),parameter                     :: rle = 0.1_r8
     integer                                :: niter_xc = 1
     integer                                :: kk,      status,  iter_xc
     real(r8)                               :: tlm,     excessm, qsm,     &   
                                               tln,     excessn, es,      &
                                               xc,      xsat,    x_en,    &
                                               x_cu,    xs1,     xs2,     &
                                               aquad,   bquad,   cquad,   &
                                               thlxsat, thvxsat, qtxsat,  &
                                               thv_x0,  thv_x1,  cridis,  &
                                               thln0,   qtn0,    wn0,     &
                                               entn,    detn,    mfn,     &
                                               ee2,     ud2
                                               
     !
     ! parameters defining initial conditions for updrafts
     real(r8),parameter                   :: pwmin = 1.5_r8,           &
                                             pwmax = 3._r8

     !
     ! alpha relates star qunataties to stddev after Suselj etal 2019
     real(r8),parameter                   :: alphw   = 0.572_r8,       &
                                             alphqt  = 2.890_r8,       &     
                                             alphthv = 2.890_r8
     !
     ! w' covariance after Suselj etal 2019
     real(r8),parameter                   :: cwqt  = 0.32_r8,          &
                                             cwthv = 0.58_r8
     !
     ! virtual mass coefficients for w-eqn after Suselj etal 2019
     real(r8),parameter                   :: wa = 1.0_r8,              &
                                             wb = 1.5_r8
     !
     ! min values to avoid singularities
     real(r8),parameter                   :: wstarmin = 1.e-3_r8,      &
                                             pblhmin  = 100._r8
     !
     ! evaporation efficiency after Suselj etal 2019
     real(r8),parameter                   :: ke = 2.5e-4_r8
     !
     ! height here downdrafts feel the surface
     real(r8),parameter                   :: z00dn = 1.e3_r8, &
                                             tinynum = 1.e-7_r8
     !
     ! to fix entrainmnet rate
     logical                              :: fixent = .false.
     !
     ! fixed entrainment rate
     real(r8),parameter                   :: fixent_ent = 2.e-4_r8
     !
     ! Arakawa and Schubert detrainment limiter
     logical                              :: do_aspd = .false.
     !
     ! Lower limit on entrainment length scale
     real(r8),parameter                   :: min_L0 = 0.5_r8
     !
     ! limiter for tke enahnced fractional entrainment
     ! (only used when do_aspd = .true.)
     real(r8),parameter                   :: max_eturb = 10._r8
     !
     ! to condensate or not to condensate
     logical                              :: do_condensation = .true.
     !
     ! use implicit method for plume updraft velocity
     logical                              :: do_implicit = .false.
     !
     ! to scale surface fluxes
     logical                              :: scalesrf = .false. 
     !
     ! minimum downdraft speed
     real(r8),parameter                   :: mindnw = 1.E-2_r8
     !
     ! limiter on cold pool effects
     real(r8),parameter                   :: max_cpfac = 5._r8

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! INITIALIZE OUTPUT VARIABLES
     dry_a     = 0._r8
     moist_a   = 0._r8
     dry_w     = 0._r8
     moist_w   = 0._r8
     dry_qt    = 0._r8
     moist_qt  = 0._r8
     dry_thl   = 0._r8
     moist_thl = 0._r8
     dry_u     = 0._r8
     moist_u   = 0._r8
     dry_v     = 0._r8
     moist_v   = 0._r8
     moist_qc  = 0._r8

     ! this is the environmental area - by default 1.
     ae        = 1._r8
     ac        = 0._r8
     aup       = 0._r8
     adn       = 0._r8
     aw        = 0._r8 
     awup      = 0._r8
     awdn      = 0._r8
     aww       = 0._r8
     awwup     = 0._r8
     awwdn     = 0._r8
     awuup     = 0._r8
     awvup     = 0._r8
     awudn     = 0._r8
     awvdn     = 0._r8
     awthvup   = 0._r8
     awthvdn   = 0._r8
     awthlup   = 0._r8
     awqtup    = 0._r8
     awthldn   = 0._r8
     awqtdn    = 0._r8
     awthl     = 0._r8
     awqt      = 0._r8
     awu       = 0._r8
     awv       = 0._r8
     thvflxup  = 0._r8
     thlflxup  = 0._r8
     qtflxup   = 0._r8
     thvflxdn  = 0._r8
     thlflxdn  = 0._r8
     qtflxdn   = 0._r8
     thvflx    = 0._r8
     thlflx    = 0._r8
     uflxup    = 0._r8
     vflxup    = 0._r8
     uflxdn    = 0._r8
     vflxdn    = 0._r8
     uflx      = 0._r8
     vflx      = 0._r8
     qtflx     = 0._r8
     sqtup     = 0._r8
     sthlup    = 0._r8
     sqtdn     = 0._r8
     sthldn    = 0._r8
     sqt       = 0._r8
     sthl      = 0._r8
     precc     = 0._r8
    
     mix       = 0._r8
     entf      = 0._r8
     enti      = 0
     det       = 0._r8

     ! START MAIN COMPUTATION
     upw   = 0._r8
     upth  = 0._r8
     upthl = 0._r8
     upthv = 0._r8
     upqt  = 0._r8
     upa   = 0._r8
     upmf  = 0._r8
     upu   = 0._r8
     upv   = 0._r8
     upqc  = 0._r8
     upth  = 0._r8
     upql  = 0._r8
     upqi  = 0._r8
     upqv  = 0._r8
     upqs  = 0._r8
     upbuoy= 0._r8
     uplmix= 0._r8
     uprr  = 0._r8
     supqt = 0._r8
     supthl= 0._r8
     upent = 0._r8
     updet = 0._r8
     upauto= 0._r8

     dnw   = 0._r8
     dna   = 0._r8
     dnu   = 0._r8
     dnv   = 0._r8
     dnqt  = 0._r8
     dnthl = 0._r8
     dnthv = 0._r8
     dnrr  = 0._r8
     dnth  = 0._r8
     dnqc  = 0._r8
     dnql  = 0._r8
     dnqi  = 0._r8
     dnqs  = 0._r8
     dnlmix= 0._r8
     sdnqt = 0._r8
     sdnthl= 0._r8

     dynamic_L0 = 0._r8
     ztop = 0._r8
     ddbot= 0

     if (bsort) then
       niter_xc = 3
       limarea = .true.
     end if

     ! unique identifier
     zcb_unset = 9999999._r8
     zcb       = zcb_unset

     ! surface buoyancy flux
     wthv = wthl+zvir*ths*wqt

     ! if surface buoyancy is positive then do mass-flux
     if ( wthv > 0._r8 ) then

       if (do_clubb_mf_mixd) then
         convh = max(cbm1,pblhmin)
       else
         convh = max(pblh,pblhmin)
       end if

       ! --------------------------------------------------------- !
       ! Initialize using Deardorff convective velocity scale      ! 
       ! --------------------------------------------------------- !

       wstar = max( wstarmin, (gravit/thv(1)*wthv*convh)**(1._r8/3._r8) )

       ! --------------------------------------------------------- !
       ! Compute cold pool feedback parameter                      ! 
       ! --------------------------------------------------------- !

       cpfac = 1._r8
       if (do_clubb_mf_coldpool) cpfac = min( (max(ddcp/wstar,1._r8))**clubb_mf_ddbeta, max_cpfac ) 
 
       ! --------------------------------------------------------- !
       ! Construct tri-variate PDF at the surface from wstar       ! 
       ! and initialize plume thv, qt, w                           !
       ! --------------------------------------------------------- !

       if (do_clubb_mf_ustar) then
         qstar   = wqt / max(wstarmin,ustar)
         thvstar = wthv / max(wstarmin,ustar)
       else
         qstar   = wqt / wstar
         thvstar = wthv / wstar
       end if
!+++ARH
       sigmaw   = alphw * wstar !* cpfac
       sigmaqt  = alphqt * abs(qstar) !* cpfac
       sigmathv = alphthv * abs(thvstar) !* cpfac

       wmin = sigmaw * pwmin
       wmax = sigmaw * pwmax

       do i=1,clubb_mf_nup
         wlv = wmin + (wmax-wmin) / (real(clubb_mf_nup,r8)) * (real(i-1, r8))
         wtv = wmin + (wmax-wmin) / (real(clubb_mf_nup,r8)) * real(i,r8)

         upw(1,i) = 0.5_r8 * (wlv+wtv)
         upa(1,i) = 0.5_r8 * erf( wtv/(sqrt(2._r8)*sigmaw) ) &
                    - 0.5_r8 * erf( wlv/(sqrt(2._r8)*sigmaw) )

         upmf(1,i)= rho_zm(1)*upa(1,i)*upw(1,i)

         upu(1,i) = u(1)
         upv(1,i) = v(1)

         upqt(1,i)  = cwqt * upw(1,i) * sigmaqt/sigmaw
         upthv(1,i) = cwthv * upw(1,i) * sigmathv/sigmaw
       enddo

       facqtu=1._r8
       facthvu=1._r8
       if (scalesrf) then
         ! scale surface fluxes
         ! (req'd for conservation if not running with zero-flux b.c.'s)
         srfwqtu = 0._r8
         srfwthvu = 0._r8
         srfarea = 0._r8
         do i=1,clubb_mf_nup
             srfwqtu=srfwqtu+upqt(1,i)*upw(1,i)*upa(1,i)
             srfwthvu=srfwthvu+upthv(1,i)*upw(1,i)*upa(1,i)
             srfarea = srfarea+upa(1,i)
         end do
         facqtu=srfarea*wqt/srfwqtu
         facthvu=srfarea*wthv/srfwthvu
       end if

       do i=1,clubb_mf_nup

         betaqt = (qt(4)-qt(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
         betathl = (thv(4)-thv(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))

         upqt(1,i)= qt(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))+facqtu*upqt(1,i)
         upthv(1,i)= thv(2)-betathl*0.5_r8*(dzt(2)+dzt(1))+facthvu*upthv(1,i)

         upthl(1,i) = upthv(1,i) / (1._r8+zvir*upqt(1,i))
         upth(1,i)  = upthl(1,i)

         ! get cloud, lowest momentum level 
         if (do_condensation) then
           call condensation_mf(upqt(1,i), upthl(1,i), p_zm(1), iexner_zm(1), &
                                thvn, qcn, thn, qln, qin, qsn, lmixn)
           upthv(1,i) = thvn
           upqc(1,i)  = qcn
           upql(1,i)  = qln
           upqi(1,i)  = qin
           upqs(1,i)  = qsn
           upth(1,i)  = thn
           if (qcn > 0._r8) zcb(i) = zm(1)
         else
           ! assume no cldliq
           upqc(1,i)  = 0._r8
         end if
       end do

       ! --------------------------------------------------------- !
       ! Calculate ztop and dynamic_L based on value of namelist   ! 
       ! --------------------------------------------------------- !
       call get_Lscale (nz, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                        wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                        do_condensation, qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1, dynamic_L0, ztop, mcape)

       ! cold pool feedback on the entrainmnet length scale
       dynamic_L0 = dynamic_L0 * cpfac

       ! limit max/min
       dynamic_L0 = max(min_L0,dynamic_L0)
       dynamic_L0 = min(clubb_mf_max_L0,dynamic_L0)

       ! --------------------------------------------------------- !
       ! Stochastic entrainmnet calculation                        ! 
       ! From Suselj et al 2019, after Romps and Kuang 2010        !
       ! (ideally we wouldn't fill the entire arrray w/ the RNG,   !
       ! but the RNG doesn't work properly when it operates on     !
       ! the entire array. I'm not sure why this is happening.)    !
       ! --------------------------------------------------------- !
       do k=1,nz-1
         ! get entrainment coefficient, dz/L0
         entf(k,:) = dzt(k) / dynamic_L0
       end do

       ! get poisson, P(dz/L0)
       call poisson( nz, clubb_mf_nup, entf, enti, u(2:5))

       ! --------------------------------------------------------- !
       ! Main upward sweep to compute updraft properties           ! 
       !                                                           !
       ! --------------------------------------------------------- !
       do i=1,clubb_mf_nup
         do k=1,nz-1

           ! get microphysics, autoconversion
           if (do_clubb_mf_precip .and. upqc(k,i) > 0._r8) then
             call precip_mf(upqs(k,i),upqt(k,i),upw(k,i),dzt(k+1),zm(k+1)-zcb(i),supqt(k+1,i))
             supthl(k+1,i) = -1._r8*lmixn*supqt(k+1,i)*iexner_zt(k+1)/cpair
           else
             supqt(k+1,i)  = 0._r8
             supthl(k+1,i) = 0._r8
           end if

           ! compute mixing rate
           if (fixent) then
             mix(k+1,i) = fixent_ent
           else
             ! get entrainment, ent=ent0/dz*P(dz/L0)
             mix(k+1,i) = real( enti(k+1,i))*clubb_mf_ent0/dzt(k+1)
           end if

           do iter_xc = 1, niter_xc

             if (bsort) then
               if (iter_xc==1) then
                 qtn  = upqt(k,i)
                 thln = upthl(k,i)
                 wn   = upw(k,i)
               else
                 qtn  = 0.5_r8*(qtn + qtn0)
                 thln = 0.5_r8*(thln + thln0)
                 wn = 0.5_r8*(wn + wn0)
               end if

               ! save this iteration
               qtn0  = qtn
               thln0 = thln
               wn0 = wn

               ! --------------------------------------------------------- !
               ! Compute excess water to derive neutral mixing fraction    ! 
               ! after Bretherton et al 2014                               !
               ! --------------------------------------------------------- !

               ! qexcess of the envrionment
               tlm = thl_zm(k+1)/iexner_zm(k+1)
               call qsat(tlm,p_zm(k+1),es,qsm)
               excessm = qt_zm(k+1) - qsm

               ! qexcess in plume
               tln  = thln/iexner_zm(k+1)
               call qsat(tln,p_zm(k+1),es,qsn)
               excessn = qtn - qsn

               call condensation_mf(qtn, thln, p_zm(k+1), iexner_zm(k+1), &
                                    thvn, qcn, thn, qln, qin, qsn, lmixn)

               ! critical stopping distance
               cridis = rle*ztopm1

               ! ----------------------------------------------------------------- !
               ! Case 1 : When both cumulus and env. are unsaturated or saturated. !
               ! ----------------------------------------------------------------- !
               if (excessm*excessn > 0._r8) then
                 xc = min(1._r8,max(0._r8,1._r8-2._r8*wa*gravit*cridis/wn**2._r8*(1._r8-thvn/thv_zm(k+1))))
                 aquad = 0._r8
                 bquad = 0._r8
                 cquad = 0._r8
               else
               ! -------------------------------------------------- !
               ! Case 2 : When either cumulus or env. is saturated. !
               ! -------------------------------------------------- !
                 xsat    = excessn / ( excessn - excessm );
                 thlxsat = thln + xsat * ( thl_zm(k+1) - thln );
                 qtxsat  = qtn  + xsat * ( qt_zm(k+1) - qtn );
                 call condensation_mf(qtxsat, thlxsat, p_zm(k+1), iexner_zm(k+1), &
                                      thvxsat, qcn, thn, qln, qin, qsn, lmixn)
                 ! -------------------------------------------------- !
                 ! kk=1 : Cumulus Segment, kk=2 : Environment Segment !
                 ! -------------------------------------------------- ! 
                 do kk = 1, 2
                   if( kk .eq. 1 ) then
                     thv_x0 = thvn
                     thv_x1 = ( 1._r8 - 1._r8/xsat ) * thvn + ( 1._r8/xsat ) * thvxsat
                   else
                     thv_x1 = thv_zm(k+1)
                     thv_x0 = ( xsat / ( xsat - 1._r8 ) ) * thv_zm(k+1) + ( 1._r8/( 1._r8 - xsat ) ) * thvxsat
                   endif
                   aquad =  wn**2
                   bquad =  2._r8*wa*gravit*cridis*(thv_x1 - thv_x0)/thv_zm(k+1) - 2._r8*wn**2
                   cquad =  2._r8*wa*gravit*cridis*(thv_x0 - thv_zm(k+1))/thv_zm(k+1)  + wn**2
                   if( kk .eq. 1 ) then
                     if( ( bquad**2-4._r8*aquad*cquad ) .ge. 0._r8 ) then
                       call roots(aquad,bquad,cquad,xs1,xs2,status)
                       x_cu = min(1._r8,max(0._r8,min(xsat,min(xs1,xs2))))
                     else
                       x_cu = xsat
                     endif
                   else
                     if( ( bquad**2-4._r8*aquad*cquad) .ge. 0._r8 ) then
                       call roots(aquad,bquad,cquad,xs1,xs2,status)
                       x_en = min(1._r8,max(0._r8,max(xsat,min(xs1,xs2))))
                     else
                       x_en = 1._r8
                     endif
                   endif
                 enddo
                 if( x_cu .eq. xsat ) then
                   xc = max(x_cu, x_en)
                 else
                   xc = x_cu
                 endif
               endif

               ee2 = xc**2
               ud2 = 1._r8 - 2._r8*xc + xc**2

               ! detrainment rate
               detn  = mix(k+1,i) * ud2

             else !no bsort
               ee2 = 1._r8
               ud2 = 1._r8
             end if

             ! entrainment rate
             entn = mix(k+1,i) * ee2

             ! --------------------------------------------------------- !
             ! TKE enhanced entrainment                                  ! 
             ! switches off when dynamic_L0 > max_L0                     !
             ! --------------------------------------------------------- !
             eturb = (1._r8 + clubb_mf_alphturb*sqrt(tke(k))/upw(k,i))
             if (do_clubb_mf_rhtke) then
               rh_L0 = 50._r8*(rhinv**3._r8)
               if (rh_L0 >= 733.34_r8) eturb = 1._r8
             else
               if (dynamic_L0 >= clubb_mf_max_L0) eturb = 1._r8
             end if
             entn = entn * eturb

             ! integrate updraft
             entexp  = exp(-entn*eturb*dzt(k+1))
             entexpu = exp(-entn*dzt(k+1)/3._r8)

             qtn  = qt(k+1) *(1._r8-entexp ) + upqt (k,i)*entexp + supqt(k+1,i)
             thln = thl(k+1)*(1._r8-entexp ) + upthl(k,i)*entexp + supthl(k+1,i)           
             un   = u(k+1)  *(1._r8-entexpu) + upu  (k,i)*entexpu
             vn   = v(k+1)  *(1._r8-entexpu) + upv  (k,i)*entexpu

             ! convert source terms to a tendency (convert from S*dz/w to S)
             supqt(k+1,i) = supqt(k+1,i)*upw(k,i)/dzt(k+1)
             upauto(k+1,i) = supqt(k+1,i)
             supthl(k+1,i) = supthl(k+1,i)*upw(k,i)/dzt(k+1)

             ! get cloud, momentum levels
             if (do_condensation) then
               call condensation_mf(qtn, thln, p_zm(k+1), iexner_zm(k+1), &
                                    thvn, qcn, thn, qln, qin, qsn, lmixn)
               if (zcb(i).eq.zcb_unset .and. qcn > 0._r8) zcb(i) = zm(k+1)
             else
               thvn = thln*(1._r8+zvir*qtn)
             end if

             ! get buoyancy
             B=gravit*(0.5_r8*(thvn + upthv(k,i))/thv(k+1)-1._r8)

             if (do_implicit) then
               wp = clubb_mf_alphturb*wb*entn*sqrt(0.5_r8*(tke(k+1)+tke(k)))*dzt(k+1)
               wn = (-wp + sqrt(wp**2._r8 + (1._r8 + 2._r8*wb*entn*dzt(k+1))* &
                     (upw(k,i)**2._r8 + 2._r8*wa*B*dzt(k+1))) )/(1._r8 + 2._r8*wb*entn*dzt(k+1))
             else
               ! get wn2
               wp = wb*entn*eturb
               if (wp==0._r8) then
                 wn2 = upw(k,i)**2._r8+2._r8*wa*B*dzt(k+1)
               else
                 entw = exp(-2._r8*wp*dzt(k+1))
                 wn2 = entw*upw(k,i)**2._r8+(1._r8-entw)*wa*B/wp
               end if
               wn = sqrt(max(wn2, 0._r8))
             end if

           end do !iter_xc

           if (wn>0._r8) then

             upthv(k+1,i) = thvn
             upthl(k+1,i) = thln
             upqt(k+1,i)  = qtn
             upqc(k+1,i)  = qcn
             upqs(k+1,i)  = qsn
             upu(k+1,i)   = un
             upv(k+1,i)   = vn
             upql(k+1,i)  = qln
             upqi(k+1,i)  = qin
             upqv(k+1,i)  = qtn - qcn
             uplmix(k+1,i)= lmixn
             upth(k+1,i)  = thn

             if (bsort) then
               mfn = upmf(k,i)*exp( dzt(k+1)*( entn - detn ))
               upa(k+1,i) = mfn/(wn*rho_zm(k+1))
             else
               upa(k+1,i) = upa(k,i)
               mfn = rho_zm(k+1)*upa(k+1,i)*wn
               detn = entn - (mfn - rho_zm(k)*upa(k,i)*upw(k,i)) &
                             /(rho_zm(k)*upa(k,i)*upw(k,i)*dzt(k+1))
             end if

             upbuoy(k+1,i)= B
             upw(k+1,i)   = wn
             upmf(k+1,i)  = mfn
             upent(k+1,i) = entn
             updet(k+1,i) = detn

           else
             ! zero out plumes that terminate at k<3
             if (k<3) then
               supqt(:,i) = 0._r8
               upauto(:,i)= 0._r8
               supthl(:,i)= 0._r8

               upa(:,i)   = 0._r8
               upbuoy(:,i)= 0._r8
               upw(:,i)   = 0._r8
               upmf(:,i)  = 0._r8
               upent(:,i) = 0._r8
               updet(:,i) = 0._r8

               upthv(:,i) = 0._r8
               upthl(:,i) = 0._r8
               upqt(:,i)  = 0._r8
               upqc(:,i)  = 0._r8
               upqs(:,i)  = 0._r8
               upu(:,i)   = 0._r8
               upv(:,i)   = 0._r8
               upql(:,i)  = 0._r8
               upqi(:,i)  = 0._r8
               upqv(:,i)  = 0._r8
               uplmix(:,i)= 0._r8
               upth(:,i)  = 0._r8
             end if
             ! exit updraft integration
             exit
             !
           end if
         enddo
       enddo

       ! --------------------------------------------------------- !
       ! downward sweep for rain evaporation, snow melting         !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip) then
         do i=1,clubb_mf_nup
           do k=nz,2,-1
             ! get rain evaporation
             if ((upqs(k,i) + upqs(k-1,i)).le.0._r8) then
               qtovqs = 0._r8
             else
               qtovqs = (upqt(k,i) + upqt(k-1,i))/(upqs(k,i) + upqs(k-1,i))
             end if
             qtovqs = min(1._r8,qtovqs)
             sevap = ke*(1._r8 - qtovqs)*sqrt(max(uprr(k,i),0._r8))

             ! limit evaporation to available precip
             sevap = min(sevap,( uprr(k,i)/(rho_zt(k)*dzt(k)) - supqt(k,i)*(1._r8-clubb_mf_fdd) ))

             ! get rain rate
             uprr(k-1,i) = uprr(k,i) &
                         - rho_zt(k)*dzt(k)*( supqt(k,i)*(1._r8-clubb_mf_fdd) + sevap )

             ! update source terms
             lmixt = 0.5_r8*(uplmix(k,i)+uplmix(k-1,i))
             supqt(k,i) = supqt(k,i) + sevap
             supthl(k,i) = supthl(k,i) - lmixt*sevap*iexner_zt(k)/cpair
           end do
         end do
       end if

       ! --------------------------------------------------------- !
       ! begin computing downdrafts                                !
       ! --------------------------------------------------------- !
       if (do_clubb_mf_precip .and. clubb_mf_fdd > 0._r8) then       

         do i=1,clubb_mf_nup

           ! find cloud base
           do k = 1,nz
             if (upqc(k,i) > 0._r8) then
               ddbot(i) = k
               exit
             end if
           end do
 
           ! find cloud top
           ddtop = 0
           do k = 1,nz
             if (uprr(k,i) > 0._r8) ddtop = k
           end do

           if (ddtop /= 0) then
             ! initilaize downdrafts

             ! Kay initializes using negative of the updraft velocity
             ! this causes anomalouly large downdrafts at the initializaiton level
             ! I am intializing with zero velocity as that is more physically defensible
             dnw(ddtop,i)   = -1._r8*mindnw !upw(ddtop,i) ! 0._r8
             dna(ddtop,i)   = upa(ddtop,i)
             dnu(ddtop,i)   = 0.5_r8*(u(ddtop)+u(ddtop+1)) 
             dnv(ddtop,i)   = 0.5_r8*(v(ddtop)+v(ddtop+1))
             dnqt(ddtop,i)  = qt_zm(ddtop)
 
             ! no cloud in downdrafts, set to cloud free thl
             dnthl(ddtop,i) = thl_zm(ddtop)
             dnthv(ddtop,i) = thv_zm(ddtop) ! includes condensate loading (!)

             ! get rain generated in the updraft, appropriate it to the downdraft
             dnrr(ddtop,i)  = -1._r8*dzt(ddtop)*rho_zt(ddtop)*upauto(ddtop,i)*clubb_mf_fdd

             if (fixent) then
               entn = fixent_ent
             else
               ! use deterministic mean entrainment
               entn = clubb_mf_ent0/dynamic_L0
             end if

             ! downdraft qsat
             call qsat(dnthl(ddtop,i)/iexner_zm(ddtop),p_zm(ddtop),es,dnqs(ddtop,i))

             do k = ddtop-1,1,-1

               ! assume fixed area with height
               dna(k,i) = dna(k+1,i)

               ! get rain evaporation in integrated form
               taum1 = ke*sqrt(dnrr(k+1,i))/dnqs(k+1,i)
               alphint = exp(dzt(k+1)*taum1/dnw(k+1,i))
               sqtint = max( (dnqs(k+1,i) - dnqt(k+1,i))*(1._r8 - alphint) ,0._r8)

               ! limit to available rain
               sqtint = min( sqtint, -1._r8*dnrr(k+1,i) / (rho_zt(k+1)*dzt(k+1)*dnw(k+1,i)) )
               sthlint = -1._r8*latvap*sqtint*iexner_zt(k+1)/cpair

               ! get rain evaporation in tendency form
               sdnqt(k,i) = max( (dnqs(k+1,i) - dnqt(k+1,i))*taum1, 0._r8 )
               sdnthl(k,i) = -1._r8*latvap*sdnqt(k,i)*iexner_zt(k+1)/cpair

               ! compute rain rate (rain above - evaporation + appropriate updraft rain)
               dnrr(k,i) = max( dnrr(k+1,i) &
                                - rho_zt(k+1)*dzt(k+1)*(sdnqt(k,i) + upauto(k+1,i)*clubb_mf_fdd) , 0._r8 ) 

               ! include eturb?
               entexp  = exp(-1._r8*entn*eturb*dzt(k+1))
               entexpu = exp(-1._r8*entn*dzt(k+1)/3._r8)

               ! integrate downward
               dnu(k,i)   = u(k+1)  *(1._r8-entexpu) + dnu  (k+1,i)*entexpu
               dnv(k,i)   = v(k+1)  *(1._r8-entexpu) + dnv  (k+1,i)*entexpu
               dnqt(k,i)  = qt(k+1) *(1._r8-entexp ) + dnqt (k+1,i)*entexp + sqtint
               dnthl(k,i) = thl(k+1)*(1._r8-entexp ) + dnthl(k+1,i)*entexp + sthlint

               ! get qsat
               call qsat(dnthl(k,i)/iexner_zm(k),p_zm(k),es,dnqs(k,i))

               ! no supersaturation in downdrafts             
               if (dnqt(k,i) > dnqs(k,i)) then
                 ! set qt to saturation vapor pressure
                 dnqt(k,i) = dnqs(k,i)

                 ! find evaporation that gives saturation vapor pressure
                 sqtint = dnqt(k,i) - (qt(k+1) *(1._r8-entexp ) + dnqt (k+1,i)*entexp)
 
                 ! limit to available rain
                 sqtint = min( sqtint, -1._r8*dnrr(k+1,i) / (rho_zt(k+1)*dzt(k+1)*dnw(k+1,i)) )
                 sthlint = -1._r8*latvap*sqtint*iexner_zt(k+1)/cpair

                 ! find new evap tendency
                 if ((alphint - 1._r8) /= 0._r8) then
                   qtmp = dnqs(k+1,i) + sqtint/(alphint - 1._r8)
                   sdnqt(k,i) = max( (dnqs(k+1,i) - qtmp)*taum1, 0._r8 )
                 else
                   sdnqt(k,i) = 0._r8
                 end if
                 sdnthl(k,i) = -1._r8*latvap*sdnqt(k,i)*iexner_zt(k+1)/cpair

                 ! re-compute thl with new evaporation rate
                 dnthl(k,i) = thl(k+1)*(1._r8-entexp ) + dnthl(k+1,i)*entexp + sthlint

                 ! adjust rain
                 dnrr(k,i) = max( dnrr(k+1,i) &
                                  - rho_zt(k+1)*dzt(k+1)*(sdnqt(k,i) + upauto(k+1,i)*clubb_mf_fdd) , 0._r8 )
               end if

               ! get virtual temperature
               dnthv(k,i) = dnthl(k,i)*(1._r8+zvir*dnqt(k,i))
     
               if (k > ddbot(i)) then
                 ! get virtual temperature
                 dnthv(k,i) = dnthl(k,i)*(1._r8+zvir*dnqt(k,i))

                 ! get buoyancy
                 ! (midpoint k is surrounded by interface k and k-1,
                 ! and therefore we can't compute B at the midpoint properly)
                 B = gravit*(dnthv(k,i)/thv(k)-1._r8)

                 ! get wn2
                 wp = wb*entn*eturb  &
                      + clubb_mf_pwfac/( 2._r8*zm(k)+tinynum ) * max( 1._r8 - exp( zm(k)/z00dn-1._r8), 0._r8 )
                 if (wp==0._r8) then
                   wn2 = dnw(k+1,i)**2._r8-2._r8*wa*B*dzt(k+1)
                 else
                   entw = exp(-2._r8*wp*dzt(k+1))
                   wn2 = entw*dnw(k+1,i)**2._r8-(1._r8-entw)*wa*B/wp
                 end if
                 wn2 = max(wn2,mindnw**2._r8)
                 dnw(k,i) = -1._r8*sqrt(wn2)

                 ! enforce net positive mass flux at cloud base
                 !if (k == (ddbot(i)+1)) then
                 !  if (sqrt(wn2) > upw(k,i)) dnw(k,i) = -1._r8*upw(k,i)
                 !end if

               else
                 zsub = zm(ddbot(i)+1)
                 wcb  = dnw(ddbot(i)+1,i)
                 dnw(k,i) = wcb - (wcb/(zsub**clubb_mf_ddexp))*(zsub - zm(k))**clubb_mf_ddexp
                 dnw(k,i) = min(dnw(k,i),-1._r8*mindnw)
               end if

             end do!k

           end if

         end do!i

!+++ARH this should be changed to only zero out above the downdraft (dnw<-mindw)
!+++ARH also this should zero out dna as well
         ! zero out downdraft fluxes for dnw == -mindnw
         do i=1,clubb_mf_nup
           do k=1,nz
             if ( dnw(k,i) == -1._r8*mindnw ) then
               dnw(k,i) = 0._r8
               dna(k,i) = 0._r8
             end if
           end do
         end do

       end if
       ! end computing downdrafts

       ! --------------------------------------------------------- !
       ! AS.pd limiter                                             ! 
       ! --------------------------------------------------------- !
       if (do_aspd) then
         do k=1,nz-1
           do i=1,clubb_mf_nup
             if (upw(k+1,i)>0._r8) then
               ! diagnose detrainment
               Mn = rho_zm(k)*upa(k,i)*upw(k,i)
               det = upent(k+1,i)*eturb - (rho_zm(k+1)*upa(k+1,i)*upw(k+1,i) - Mn) &
                             /(Mn*dzt(k+1))
               if (det < 0._r8) then
                 ! diagnose area to eliminate detrainment and conserve mass
                 Mn = rho_zm(k)*upa(k,i)*upw(k,i)*exp(upent(k+1,i)*eturb*dzt(k+1))
                 upa(k+1,i) = Mn/(rho_zm(k+1)*upw(k+1,i))
               end if
               !
             end if
           end do
         end do
       end if

       ! --------------------------------------------------------- !
       ! integrate for total convective area                       ! 
       ! --------------------------------------------------------- !
       do k=1,nz-1
         !
         do i=1,clubb_mf_nup
           aup(k) = aup(k) + upa(k,i)
           adn(k) = adn(k) + dna(k,i)
         end do
         ac(k) = aup(k) + adn(k)
         !
         if (limarea .and. ac(k) > amax) then
           upa(k,:) = upa(k,:)*amax/ac(k)
           ac(k) = amax
         end if
         ae(k) = ae(k) - ac(k)
         !
       end do

       ! --------------------------------------------------------- !
       ! updraft properties for output                             ! 
       ! --------------------------------------------------------- !
       do k=1,nz

         ! first sum over all i-updrafts
         do i=1,clubb_mf_nup
           if (upqc(k,i)>0._r8) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0._r8 ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0._r8
           dry_qt(k)  = 0._r8
           dry_thl(k) = 0._r8
           dry_u(k)   = 0._r8
           dry_v(k)   = 0._r8
         endif

         if ( moist_a(k) > 0._r8 ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0._r8
           moist_qt(k)  = 0._r8
           moist_thl(k) = 0._r8
           moist_u(k)   = 0._r8
           moist_v(k)   = 0._r8
           moist_qc(k)  = 0._r8
         endif

       enddo

       ! --------------------------------------------------------- !
       ! get ensemble mean                                         ! 
       ! --------------------------------------------------------- !
       do k=1,nz
         do i=1,clubb_mf_nup

           awup(k) = awup(k) + upa(k,i)*upw(k,i)
           awdn(k) = awdn(k) + dna(k,i)*dnw(k,i)

           awwup(k) = awwup(k) + upa(k,i)*upw(k,i)*upw(k,i)
           awwdn(k) = awwdn(k) + dna(k,i)*dnw(k,i)*dnw(k,i)

           awuup(k) = awuup(k) + upa(k,i)*upw(k,i)*upu(k,i)
           awudn(k) = awudn(k) + dna(k,i)*dnw(k,i)*dnu(k,i)
           awvup(k) = awvup(k) + upa(k,i)*upw(k,i)*upv(k,i)
           awvdn(k) = awvdn(k) + dna(k,i)*dnw(k,i)*dnv(k,i)

           awthvdn(k)= awthvdn(k)+ dna(k,i)*dnw(k,i)*dnthv(k,i)
           awthldn(k)= awthldn(k)+ dna(k,i)*dnw(k,i)*dnthl(k,i)
           awqtdn(k) = awqtdn(k) + dna(k,i)*dnw(k,i)*dnqt(k,i)

           awthvup(k)= awthvup(k)+ upa(k,i)*upw(k,i)*upthv(k,i)
           awthlup(k)= awthlup(k)+ upa(k,i)*upw(k,i)*upthl(k,i)
           awqtup(k) = awqtup(k) + upa(k,i)*upw(k,i)*upqt(k,i) 

           if (k > 1) then
             sqtup(k)  = sqtup(k)  + upa(k-1,i)*supqt(k,i)  
             sthlup(k) = sthlup(k) + upa(k-1,i)*supthl(k,i) 

             sqtdn(k)  = sqtdn(k)  + dna(k,i)*sdnqt(k,i)
             sthldn(k) = sthldn(k) + dna(k,i)*sdnthl(k,i)
           end if

         enddo

         aw (k) = awup(k)+ awdn(k)
         aww(k) = awwup(k)+ awwdn(k)
         awu(k) = awuup(k)+ awudn(k)
         awv(k) = awvup(k)+ awvdn(k)
         sqt(k) = sqtup(k) + sqtdn(k)
         sthl(k)= sthlup(k) + sthldn(k)

       enddo

       ! --------------------------------------------------------- !
       ! ztopm1 calculation                                        ! 
       ! --------------------------------------------------------- !
       do k=1,nz
         ! retrun if no convection at k=2
         if (k == 2 .and. ac(k) == 0._r8) then
           sqt(k) = 0_r8
           sthl(k) = 0._r8
           ztopm1 = zm(1)
           ddcp = 0._r8
           return
         end if
         ! height of the plume ensemble
         if (ac(k) > 0._r8) ztopm1 = zm(k)
       end do

       ! --------------------------------------------------------- !
       ! cloud base / mixing depth calculation                                        ! 
       ! --------------------------------------------------------- !
       cbm1 = 0._r8
       do i=1,clubb_mf_nup
         kcbarr(i) = 0
         do k=1,nz
           if (upqc(k,i) > 0._r8) then
             kcbarr(i) = k
             exit
           end if
         end do

         ! find height of dry plumes
         if (kcbarr(i) == 0) then
           do k=1,nz
             if (upw(k,i) <= 0._r8) then
               kcbarr(i) = k
               exit
             end if
           end do
         end if

         cbm1 = cbm1 + zm(kcbarr(i))

       end do
       cbm1 = cbm1/REAL(clubb_mf_nup)

       ! --------------------------------------------------------- !
       ! bulk downdraft velocity for coldpool parameterization     ! 
       ! --------------------------------------------------------- !
!+++ARH
!       ! reset ddcp
!       ddcp = 0._r8
!       do i=1,clubb_mf_nup
!         ! find cloud base
!         kcb = 0
!         do k=1,nz
!           if (upqc(k,i) > 0._r8) then
!             kcb = k
!             exit
!           end if
!         end do
! 
!         ! reset iddcp
!         iddcp = 0._r8
!         if (kcb == 0) then
!           continue
!         else if (kcb == 1) then
!           iddcp = iddcp + dna(k,i)*dnw(k,i)
!           continue
!         else
!           ddint = 0._r8
!           do k=1,kcb-1
!             ddint = ddint + dna(k,i)*dnw(k,i)*dzt(k+1)
!           end do
!           iddcp = iddcp + -1._r8*ddint/zm(kcb)
!         end if
!         ddcp = ddcp + iddcp 
!         !
!       end do
!

       ddcp = 0._r8
       if (do_clubb_mf_coldpool .and. clubb_mf_fdd > 0._r8) then
         ! use single level for cold pool param.
         ! reset ddcp
         do i=1,clubb_mf_nup
           if (ddbot(i) == 0) then
             continue
           else
             ddcp = ddcp + (-1._r8)*dna(ddbot(i)+1,i)*dnw(ddbot(i)+1,i)
           end if
         end do
       end if
!---ARH

       ! --------------------------------------------------------- !
       ! downward sweep to get ensemble mean precip                ! 
       ! --------------------------------------------------------- !
       do k = nz,2,-1
         precc(k-1) = precc(k) - rho_zt(k)*dzt(k)*sqt(k)
       end do

       ! --------------------------------------------------------- !
       ! get turbulent fluxes                                      ! 
       ! --------------------------------------------------------- !
       thv_env = thv
       thl_env = thl
       qt_env  = qt

       betathl = (thl_env(4)-thl_env(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))
       betaqt = (qt_env(4)-qt_env(2))/(0.5_r8*(dzt(4)+2._r8*dzt(3)+dzt(2)))

       thl_env(1) = thl_env(2)-betathl*0.5_r8*(dzt(2)+dzt(1))
       qt_env(1) = qt_env(2)-betaqt*0.5_r8*(dzt(2)+dzt(1))
       if (qt_env(1) < 0._r8) qt_env(1) = 0._r8

       kstart = 2
       if (scalesrf) then
         kstart = 1
       end if

       do k=kstart,nz-1
         thvflxup(k)= awthvup(k) - awup(k)*thv_env(k+1)
         thlflxup(k)= awthlup(k) - awup(k)*thl_env(k+1)
         qtflxup (k)= awqtup (k) - awup(k)*qt_env (k+1)

         uflxup  (k)= awuup(k) - awup(k)*u(k+1)
         vflxup  (k)= awvup(k) - awup(k)*v(k+1)

         ! if no downdrafts, should be zero since awdn should be zero
         thvflxdn(k)= awthvdn(k) - awdn(k)*thv_env(k)
         thlflxdn(k)= awthldn(k) - awdn(k)*thl_env(k)
         qtflxdn (k)= awqtdn (k) - awdn(k)*qt_env (k)

         uflxdn  (k)= awudn(k) - awdn(k)*u(k)
         vflxdn  (k)= awvdn(k) - awdn(k)*v(k)

         thvflx(k)  = thvflxup(k) + thvflxdn(k)
         thlflx(k)  = thlflxup(k) + thlflxdn(k)
         qtflx (k)  = qtflxup (k) + qtflxdn (k)
        
         uflx(k)    = uflxup(k) + uflxdn(k)
         vflx(k)    = vflxup(k) + vflxdn(k) 
       enddo
       !
     else
       ddcp = 0._r8
       ztopm1 = zm(1)
     end if  ! ( wthv > 0.0 )

  end subroutine integrate_mf

  subroutine get_Lscale(nz, zm, tke, wpthlp_env, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, th, &
                        wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, wa, wb,  &
                        do_condensation, qv, p_zt, zt, tpert, pblh, convh, rhinv, ztopm1, dynamic_L0, ztop, mcape) 
  ! --------------------------------------------------------- !
  ! Calculate ztop and dynamic_L based on value of namelist   ! 
  ! --------------------------------------------------------- !
     integer,  intent(in)                :: nz
     real(r8), dimension(nz), intent(in) :: thl,    thv,          &
                                            th,                   &
                                            qt,     qv,           &
                                            p_zt,   iexner_zt,    &
                                            dzt,    zt,           &
                                            p_zm,   iexner_zm,    &
                                            zm,                   &
                                            tke,    wpthlp_env     
     !
     real(r8), intent(in) ::                wmax,   wmin,       tpert, &
                                            sigmaw, sigmaqt, sigmathv, &
                                            cwqt,   cwthv,  zcb_unset, &
                                            wa,     wb,        ztopm1, &
                                            pblh,   convh,     rhinv
     !
     logical, intent(in) ::                 do_condensation                
     !
     real(r8), intent(out) ::               dynamic_L0, ztop, mcape 
     !
     ! local variables
     real(r8), dimension(nz)                :: t_zt
     real(r8), dimension(nz-1)              :: tp,       qstp
     !real(r8), dimension(nz-1,clubb_mf_nup) :: dmpdz
     !real(r8), dimension(clubb_mf_nup)      :: tl,                     &
     !                                          cape,     cin
     !integer,  dimension(clubb_mf_nup)      :: lcl,      lel
     real(r8), dimension(nz-1,1)            :: dmpdz
     real(r8), dimension(1)                 :: tl,                     &
                                               cape,     cin
     integer,  dimension(1)                 :: lcl,      lel
     real(r8)                               :: landfrac
     integer                                :: kpbl,     msg,          &
                                               lon,      mx,           &
                                               k

     ! intialize local variables
     cape      = 0._r8
     mcape     = 0._r8
     dmpdz     = 0._r8

     if (clubb_mf_Lopt==0) then
       !Constant L0
       dynamic_L0 = clubb_mf_L0
       ztop = clubb_mf_L0
     else if (clubb_mf_Lopt==1) then
       !TKE
       do k=nz-2,2,-1
         if (zm(k) < 20000 .and. tke(k) - tke(k+1) > 1e-5) then
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==2) then
       !Heat flux
       do k=nz-2,2,-1
         !if (zm(k) < 20000 .and. abs(abs(wpthlp_env(k))-abs(wpthlp_env(k-1))) > 1e-3) then
         if (zm(k) < 20000 .and. abs(abs(wpthlp_env(k))-abs(wpthlp_env(k-1))) > 1e-4) then
           ztop = zm(k)
           exit
         endif
       enddo
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==3) then
       !Test plume
       call oneplume( nz, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                      wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                      wa, wb, tke, do_condensation, do_clubb_mf_precip, ztop )

       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
       !ztop = ztop - 1600._r8
       !if (ztop < 1._r8) then
       !  dynamic_L0 = clubb_mf_a0
       !else
       !  dynamic_L0 = min(35._r8,clubb_mf_a0*(ztop**clubb_mf_b0))
       !end if
     else if (clubb_mf_Lopt==4 .or. clubb_mf_Lopt==5) then
       !dilute cape calculation
       !dmpdz = -1._r8*ent_zt(2:nz,:)
       dmpdz(:,:) = -1.E-3_r8
       t_zt = th/iexner_zt
       landfrac = 1._r8

       do k=2,nz
         if (zt(k-1) <= pblh) then
           kpbl = k
         end if
       end do

       do k=1,nz
         if (p_zt(k) > 40.e2_r8) then
           msg = k
         end if
       end do
       !call buoyan_dilute(nz-1       ,clubb_mf_nup ,dmpdz , &
       call buoyan_dilute(nz-1       ,1          ,dmpdz , &
                          qv(2:nz)   ,t_zt(2:nz) ,p_zt(2:nz)*0.01_r8 ,zt(2:nz) ,p_zm*0.01_r8 , &
                          tp         ,qstp       ,tl         ,cape     ,cin  , &
                          kpbl-1     ,lcl        ,lel        ,lon      ,mx   , &
                          msg-1      ,tpert      ,landfrac )

       !do i=1,clubb_mf_nup
       !  mcape = mcape + cape(i)
       !end do
       !mcape = mcape/REAL(clubb_mf_nup)
       mcape = max(cape(1),25._r8)

       if (clubb_mf_Lopt==4) then
         ztop = max(zt(lel(1)+1),convh)
       else if (clubb_mf_Lopt==5) then
         ztop = mcape
       end if
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)

     else if (clubb_mf_Lopt==6) then
       ! grab ztop from max height of ensemble in prior time-step(s)
       ztop = ztopm1
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
       !if (masterproc) write(iam+110,*) 'mf_ztop, dynamic_L0 ', ztop, dynamic_L0
     else if (clubb_mf_Lopt==7 .or. clubb_mf_Lopt==8) then
       ztop = rhinv
       dynamic_L0 = clubb_mf_a0*(ztop**clubb_mf_b0)
     end if

  end subroutine get_Lscale 

  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi, qs, lmix )
  ! =============================================================================== !
  ! zero or one condensation for edmf: calculates thv and qc                        !
  ! =============================================================================== !
     use physconst,          only: cpair, zvir, h2otrip
     use wv_saturation,      only : qsat

     real(r8),intent(in) :: qt,thl,p,iex
     real(r8),intent(out):: thv,qc,th,ql,qi,qs,lmix

     !local variables
     integer  :: niter,i
     real(r8) :: diff,t,qstmp,qcold,es,wf
     logical  :: noice = .true.

     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5_r8

     qc=0._r8
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter

       if (noice) then
         wf = 1._r8
       else    
         wf = get_watf(t)
       end if
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qstmp)
       qcold = qc
       qc = max(0.5_r8*qc+0.5_r8*(qt-qstmp),0._r8)
       if (abs(qc-qcold)<diff) exit
     enddo

     if (noice) then
       wf = 1._r8
     else
       wf = get_watf(t)
     end if
     t = thl/iex+get_alhl(wf)/cpair*qc

     call qsat(t,p,es,qs)
     qc = max(qt-qs,0._r8)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1._r8+zvir*(qt-qc)-qc)
     lmix = get_alhl(wf)
     th = t*iex
     qi = qc*(1._r8-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(r8)            :: t,get_watf,tc
       real(r8), parameter :: &
                              tmax=-10._r8, &
                              tmin=-40._r8

       tc=t-h2otrip

       if (tc>tmax) then
         get_watf=1._r8
       else if (tc<tmin) then
         get_watf=0._r8
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice
       real(r8) :: get_alhl,wf

       get_alhl = wf*latvap+(1._r8-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf

  subroutine precip_mf(qs,qt,w,dz,dzcld,Supqt)
  !**********************************************************************
  ! Precipitation microphysics
  ! By Adam Herrington, after Kay Suselj
  !**********************************************************************

       real(r8),intent(in)  :: qs,qt,w,dz,dzcld
       real(r8),intent(out) :: Supqt
       ! 
       ! local vars
       real(r8)            :: tauwgt, tau,       & ! time-scale vars
                              qstar                ! excess cloud liquid                   

       real(r8),parameter  :: tau0  = 15._r8,    & ! base time-scale
                              zmin  = 300._r8,   & ! small cloud thick
                              zmax  = 3000._r8,  & ! large cloud thick
                              qcmin = 0.00125_r8   ! supersat threshold 

       qstar = qs+qcmin
       
       if (qt > qstar) then
         ! get precip efficiency
         tauwgt = (dzcld-zmin)/(zmax-zmin)
         tauwgt = min(max(tauwgt,0._r8),1._r8)
         tau    = tauwgt/tau0
 
         ! get source for updraft
         Supqt = (qstar-qt)*(1._r8 - exp(-1._r8*tau*dz/w))
       else
         Supqt = 0._r8
       end if

  end subroutine precip_mf

  subroutine poisson(nz,nup,lambda,poi,state)
  !**********************************************************************
  ! Set a unique (but reproduceble) seed for the kiss RNG
  ! Call Poisson deviate
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       integer,                     intent(in)  :: nz,nup
       real(r8), dimension(4),      intent(in)  :: state
       real(r8), dimension(nz,nup), intent(in)  :: lambda
       integer,  dimension(nz,nup), intent(out) :: poi
       integer,  dimension(1,4)                 :: tmpseed
       integer                                  :: i,j
       type(ShrKissRandGen)                     :: kiss_gen

       ! Compute seed
       tmpseed(1,1) = int((state(1) - int(state(1))) * 1000000000._r8)
       tmpseed(1,2) = int((state(2) - int(state(2))) * 1000000000._r8)
       tmpseed(1,3) = int((state(3) - int(state(3))) * 1000000000._r8)
       tmpseed(1,4) = int((state(4) - int(state(4))) * 1000000000._r8)

       ! Set seed
       kiss_gen = ShrKissRandGen(tmpseed)

       do i=1,nz
         do j=1,nup
           call hybridRNG(kiss_gen,lambda(i,j),poi(i,j))
         enddo
       enddo

  end subroutine poisson

  subroutine hybridRNG(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Interface for the two poisson rng subroutines
  ! chooses the appropriate subroutine based on the value of lambda
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(r8),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       if (lambda < 10._r8) then
          call knuth(kiss_gen,lambda,kout)
       else
          call hormann(kiss_gen,lambda,kout)
       end if

  end subroutine hybridRNG

  subroutine knuth(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth 
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************
   use shr_RandNum_mod, only: ShrKissRandGen

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(r8),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       ! Local variables
       real(r8), dimension(1,1) :: tmpuni
       real(r8)                 :: puni, explam
       integer                  :: k

       k = 0
       explam = exp(-1._r8*lambda)
       puni = 1._r8
       do while (puni > explam)
         k = k + 1
         call kiss_gen%random(tmpuni)
         puni = puni*tmpuni(1,1)
       end do
       kout = k - 1

  end subroutine knuth

  subroutine hormann(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson
  ! Implements Poisson Transformed Rejection with Squeeze (PTRS) 
  ! from W. Hormann Insurance: Mathematics and Economics 12, 39-45 (1993) 
  ! By Jake Reschke
  !**********************************************************************
  use shr_RandNum_mod, only: ShrKissRandGen

      type(ShrKissRandGen), intent(inout) :: kiss_gen
      real(r8),             intent(in)    :: lambda
      integer,              intent(out)   :: kout

      ! Local variables
      real(r8), dimension(1,1) :: U,V
      real(r8)                 :: a,b,vr,alphinv,us,loggam
      integer                  :: k,i

      b = 0.931_r8 + 2.53_r8*sqrt(lambda)
      a = -0.059_r8 + 0.02483_r8*b
      vr = 0.9277_r8 - 3.6224_r8/(b - 2._r8)
      alphinv = 1.1239_r8 + 1.1328_r8/(b - 3.4_r8)

      do
         call kiss_gen%random(U)
         call kiss_gen%random(V)
         U(1,1) = U(1,1) - 0.5_r8
         us = 0.5_r8 - abs(U(1,1))
         k = floor( (2._r8*a/us + b)*U(1,1) + lambda + 0.43_r8 )
         if (us >= 0.07_r8 .and.  V(1,1) <= vr) then
            kout = k
            exit
         end if
         if (k <= 0 .or. (us < 0.013_r8 .and. V(1,1) > us)) then
            cycle
         end if
         ! compute log(k!). If k >=10 use stirling's approximation
         if (k < 10) then
            loggam = 0._r8
            do i = 1, k
               loggam = loggam + log(1._r8*i)
            end do
         else
            loggam = log(sqrt(2._r8*pi)) + (k + 0.5_r8)*log(1._r8*k) - k + (1._r8/12._r8 - 1._r8/(360._r8*k*k))/k
         end if
         if (log( V(1,1)*alphinv/(a/(us*us) + b) ) <= -1._r8*lambda + k*log(lambda) - loggam) then
            kout = k
            exit
         end if
      end do

  end subroutine hormann

  subroutine roots(a,b,c,r1,r2,status)
  ! --------------------------------------------------------- !
  ! Subroutine to solve the second order polynomial equation. !
  ! after uwshcu.F90                                          !
  ! --------------------------------------------------------- !
    real(r8), intent(in)  :: a
    real(r8), intent(in)  :: b
    real(r8), intent(in)  :: c
    real(r8), intent(out) :: r1
    real(r8), intent(out) :: r2
    integer , intent(out) :: status
    real(r8)              :: q

    status = 0

    if( a .eq. 0._r8 ) then                            ! Form b*x + c = 0
        if( b .eq. 0._r8 ) then                        ! Failure: c = 0
            status = 1
        else                                           ! b*x + c = 0
            r1 = -c/b
        endif
        r2 = r1
    else
        if( b .eq. 0._r8 ) then                        ! Form a*x**2 + c = 0
            if( a*c .gt. 0._r8 ) then                  ! Failure: x**2 = -c/a < 0
                status = 2
            else                                       ! x**2 = -c/a 
                r1 = sqrt(-c/a)
            endif
            r2 = -r1
       else                                            ! Form a*x**2 + b*x + c = 0
            if( (b**2 - 4._r8*a*c) .lt. 0._r8 ) then   ! Failure, no real roots
                 status = 3
            else
                 q  = -0.5_r8*(b + sign(1.0_r8,b)*sqrt(b**2 - 4._r8*a*c))
                 r1 =  q/a
                 r2 =  c/q
            endif
       endif
    endif

    return

  end subroutine roots

  subroutine oneplume( nz, zm, dzt, iexner_zm, iexner_zt, p_zm, qt, thv, thl, &
                       wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, cwthv, zcb_unset, &
                       wa, wb, tke, do_condensation, do_precip, plumeheight )
  !**********************************************************************
  ! Calculate a single plume with fixed entrainment
  ! to be used for a dynamic mixing length calculation
  ! By Rachel Storer
  !**********************************************************************
    use physconst,          only: cpair, gravit, zvir

    integer,  intent(in)                :: nz
    real(r8), dimension(nz), intent(in) :: zm, dzt, iexner_zm, iexner_zt,  &
                                           p_zm, qt, thv, thl, tke
    real(r8), intent(in)                :: wmax, wmin, sigmaw, sigmaqt, sigmathv, cwqt, &
                                           cwthv, zcb_unset, wa, wb
    logical, intent(in)               :: do_condensation, do_precip

    real(r8), intent(inout)             :: plumeheight

    !local variables
    integer                     :: k
    real(r8)                    :: thvn, qtn, thln, qcn, thn, qln, qin, qsn, lmixn, zcb, B, wn2, pentexp, pturb, pentw, wp
    real(r8), dimension(nz)     :: upw, upa, upqt, upthv, upthl, upth, upqs, &
                                   upqc, upql, upqi, supqt, supthl
    !
    ! fractional entrainment rate
    real(r8), parameter         :: pent = 1.E-3_r8
    !
    ! use tke enhanced entrainment
    logical                     :: do_tptke = .false.

    zcb = zcb_unset

    upw(1) = 0.5_r8 * wmax
    upa(1) = 0.5_r8 * erf( wmax/(sqrt(2._r8)*sigmaw) )

    upqt(1)  = cwqt * upw(1) * sigmaqt/sigmaw
    upthv(1) = cwthv * upw(1) * sigmathv/sigmaw

    upqt(1) = qt(1)+upqt(1)
    upthv(1) = thv(1)+upthv(1)
    upthl(1) = upthv(1) / (1._r8+zvir*upqt(1))
    upth(1)  = upthl(1)

    ! get cloud, lowest momentum level
    if (do_condensation) then
      call condensation_mf(upqt(1), upthl(1), p_zm(1), iexner_zm(1), &
                           thvn, qcn, thn, qln, qin, qsn, lmixn)
      upthv(1) = thvn
      upqc(1)  = qcn
      upql(1)  = qln
      upqi(1)  = qin
      upqs(1)  = qsn
      upth(1)  = thn
      if (qcn > 0._r8) zcb = zm(1)
    else
      ! assume no cldliq
      upqc(1)  = 0._r8
    end if

    do k=1,nz-1
      ! get microphysics, autoconversion
      if (do_precip .and. upqc(k) > 0._r8) then
        call precip_mf(upqs(k),upqt(k),upw(k),dzt(k+1),zm(k+1)-zcb,supqt(k+1))
        supthl(k+1) = -1._r8*lmixn*supqt(k+1)*iexner_zt(k+1)/cpair
      else
        supqt(k+1)  = 0._r8
        supthl(k+1) = 0._r8
      end if
      ! integrate updraft
      if (do_tptke) then
        pturb = (1._r8 + clubb_mf_alphturb*sqrt(tke(k))/upw(k))
      else
        pturb = 1._r8
      end if
      pentexp  = exp(-pent*pturb*dzt(k+1))
      qtn  = qt(k+1) *(1._r8-pentexp ) + upqt (k)*pentexp + supqt(k+1)
      thln = thl(k+1)*(1._r8-pentexp ) + upthl(k)*pentexp + supthl(k+1)

      ! convert source terms to a tendency
      supqt(k+1) = supqt(k+1)*upw(k)/dzt(k+1)
      supthl(k+1) = supthl(k+1)*upw(k)/dzt(k+1)

      ! get cloud, momentum levels
      if (do_condensation) then
        call condensation_mf(qtn, thln, p_zm(k+1), iexner_zm(k+1), &
                             thvn, qcn, thn, qln, qin, qsn, lmixn)
        if (zcb.eq.zcb_unset .and. qcn > 0._r8) zcb = zm(k+1)
      else
        thvn = thln*(1._r8+zvir*qtn)
      end if
      ! get buoyancy
      B=gravit*(0.5_r8*(thvn + upthv(k))/thv(k+1)-1._r8)

      ! get wn^2
      wp = wb*pent*pturb
      if (wp==0._r8) then
         wn2 = upw(k)**2._r8+2._r8*wa*B*dzt(k+1)
      else
         pentw = exp(-2._r8*wp*dzt(k+1))
         wn2 = pentw*upw(k)**2._r8+(1._r8-pentw)*wa*B/wp
      end if

      if (wn2>0._r8) then
        upw(k+1)   = sqrt(wn2)
        upthv(k+1) = thvn
        upthl(k+1) = thln
        upqt(k+1)  = qtn
        upqc(k+1)  = qcn
        upqs(k+1)  = qsn
        upa(k+1)   = upa(k)
        upql(k+1)  = qln
        upqi(k+1)  = qin
        upth(k+1)  = thn
      else
        plumeheight = zm(k)
        exit
      end if
    enddo

  end subroutine oneplume

subroutine buoyan_dilute(  nz      ,nup     ,dmpdz   , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,cape    ,cin     , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  msg     ,tpert   ,landfrac )                  
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: nz            ! vertical grid
   integer, intent(in) :: nup           ! number of plumes

!+tht
   !real(r8), intent(in), dimension(nz,nup) :: dmpdz ! Parcel fractional mass entrainment rate (/m) 3D
   real(r8), intent(in) :: dmpdz(nz,nup)
   !real(r8), intent(inout) :: dmpdz(nz,nup)
!-tht

   real(r8), intent(in) :: q(nz)        ! spec. humidity
   real(r8), intent(in) :: t(nz)        ! temperature
   real(r8), intent(in) :: p(nz)        ! pressure
   real(r8), intent(in) :: z(nz)        ! height
   real(r8), intent(in) :: pf(nz+1)     ! pressure at interfaces
   integer,  intent(in) :: pblt         ! index of pbl depth
   real(r8), intent(in) :: tpert        ! perturbation temperature by pbl processes
   real(r8), intent(in) :: landfrac

!
! output arguments
!
   real(r8), intent(out) :: tp(nz,nup)       ! parcel temperature
   real(r8), intent(out) :: qstp(nz,nup)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(r8), intent(out) :: tl(nup)          ! parcel temperature at lcl
   real(r8), intent(out) :: cape(nup)        ! convective aval. pot. energy.

   real(r8), intent(out) :: cin (nup)        !+tht: CIN

   integer, intent(out)  :: lcl(nup)                          !
   integer, intent(out)  :: lel(nup)                          !
   integer, intent(out)  :: lon                               ! level of onset of deep convection
   integer, intent(out)  :: mx                                ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
   integer lelten(nup,mf_num_cin)
   real(r8) capeten(nup,mf_num_cin)     ! provisional value of cape
   real(r8) cinten(nup,mf_num_cin)     !+tht provisional value of CIN
   real(r8) tv(nz)       
   real(r8) tpv(nz,nup)      
   real(r8) buoy(nz,nup)
   real(r8) pl(nup)

   real(r8) a1
   real(r8) a2
   real(r8) estp
   real(r8) plexp
   real(r8) hmax
   real(r8) hmn
   real(r8) y

   logical plge600(nup)
   integer knt(nup)

   real(r8) e

   integer i
   integer k
   integer msg
   integer n

   real(r8), parameter :: tiedke_add = 0.5_r8
!
!-----------------------------------------------------------------------
!
   do n = 1,mf_num_cin
      do i = 1,nup
         lelten(i,n)  = 1
         capeten(i,n) = 0._r8
         cinten (i,n) = 0._r8
      end do
   end do
!
   lon = 1
   mx   = lon
   hmax = 0._r8

   do i = 1,nup
      knt(i) = 0
      lel(i) = 1
      cape(i) = 0._r8
      tp(:nz,i) = t(:nz)
      qstp(:nz,i) = q(:nz)
   end do

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   if (tht_tweaks) then 
!+tht use system constants
    tv  (:nz) = t(:nz) *(1._r8+q(:nz)/epsilo)/ (1._r8+q(:nz)) !+tht
   else
    tv  (:nz) = t(:nz) *(1._r8+1.608_r8*q(:nz))/ (1._r8+q(:nz))
   endif
!-tht
   do i = 1,nup
     tpv (:nz,i) = tv(:nz)
   end do
   buoy(:nz,:) = 0._r8

!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
   do k = 1,msg-1
!+tht: use total mse -- moist thermo
      !hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
       hmn =(cpair+q(k)*cpliq)*t(k)/(1._r8+q(k)) + (1._r8+q(k)/epsilo)/(1._r8+q(k))*gravit*z(k) &
              +(latvap-(cpliq-cpwv)*(t(k)-tmelt))*q(k)
!-tht
       if (k <= pblt .and. k >= lon .and. hmn > hmax) then
          hmax = hmn
          mx = k
       end if
   end do

! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

   do i = 1,nup ! Initialise LCL variables.
      lcl(i) = mx
      tl(i) = t(mx)
      pl(i) = p(mx)
   end do

!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!

!+tht: add geop.height in argument to allow enthalpy mixing
   call parcel_dilute(nz, nup, msg, mx, p, z, t, q, &
   tpert, tp, tpv, qstp, pl, tl, lcl, &
   landfrac, dmpdz)
!-tht


! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,nup
      plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
   end do

!
! Main buoyancy calculation.
!
   do k = 1,msg-1
      do i=1,nup
         if (k >= mx .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
          if (tht_tweaks) then 
            tv(k) = t(k)* (1._r8+q(k)/epsilo)/ (1._r8+q(k))     !+tht
          else
            tv(k) = t(k)* (1._r8+1.608_r8*q(k))/ (1._r8+q(k)) !orig
          endif
! +0.5K or not? (arbitrary at this point - introduce in parcel_dilute instead? tht)
            buoy(k,i) = tpv(k,i) - tv(k) + tiedke_add  ! +0.5K or not?
         else
            qstp(k,i) = q(k)
            tp(k,i)   = t(k)            
            tpv(k,i)  = tv(k)
         endif
      end do
   end do

!-------------------------------------------------------------------------------
! beginning from one below top (first level p>40hPa, msg) check for at most 
! num_cin levels of neutral buoyancy (LELten) and compute CAPEten between LCL 
! and those (tht)
   do k = msg-2,1,-1
      do i = 1,nup
         if (k > lcl(i) .and. plge600(i)) then
            if (buoy(k-1,i) > 0._r8 .and. buoy(k,i) <= 0._r8) then
               knt(i) = min(mf_num_cin,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do

! calculate convective available potential energy (cape).
   do n = 1,mf_num_cin
      do k = msg-1,1,-1
         do i = 1,nup
            if (plge600(i) .and. k >= mx .and. k < lelten(i,n)) then
               !capeten(i,n) = capeten(i,n) + rair*buoy(k,i)*log(pf(k-1)/pf(k))
               capeten(i,n) = capeten(i,n) + rair*buoy(k,i)*log(pf(k)/pf(k+1))
!+tht also compute total CIN
               !cinten (i,n) = cinten (i,n) - rair*min(buoy(k,i),0._r8)*log(pf(k-1)/pf(k))
               cinten (i,n) = cinten (i,n) - rair*min(buoy(k,i),0._r8)*log(pf(k)/pf(k+1))
!-tht
            end if
         end do
      end do
   end do

!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,mf_num_cin
      do i = 1,nup
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            cin (i) = cinten (i,n) !+tht CIN
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,nup
      cape(i) = max(cape(i), 0._r8)
   end do
!
   return
end subroutine buoyan_dilute

!+tht
 subroutine parcel_dilute (nz, nup, msg, klaunch, p, z, t, q, &
  tpert, tp, tpv, qstp, pl, tl, lcl, &
  landfrac, dmpdz)
!-tht

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: nz
integer, intent(in) :: nup
integer, intent(in) :: msg
integer, intent(in) :: klaunch

real(r8), intent(in)                :: tpert ! PBL temperature perturbation.
real(r8), intent(in)                :: landfrac
real(r8), intent(in), dimension(nz) :: p
!+tht
real(r8), intent(in), dimension(nz) :: z
!-tht
real(r8), intent(in), dimension(nz) :: t
real(r8), intent(in), dimension(nz) :: q

real(r8), intent(inout), dimension(nz,nup) :: tp    ! Parcel temp.
real(r8), intent(inout), dimension(nz,nup) :: qstp  ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(nup)    :: tl         ! Actual temp of LCL.
real(r8), intent(inout), dimension(nup)    :: pl          ! Actual pressure of LCL. 
integer,  intent(inout), dimension(nup)    :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(nz,nup)   :: tpv   ! Define tpv within this routine.

!+tht
!real(r8), dimension(pcols)      :: dmpdz ! Parcel fractional mass entrainment rate (/m) 2D
 real(r8), dimension(nz,nup) :: dmpdz ! Parcel fractional mass entrainment rate (/m) 3D
!-tht

!--------------------

! Have to be careful as s is also dry static energy.
!+tht
! in the mods below, s is used both as enthalpy (moist s.e.) and entropy
!-tht

! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(nz,nup)        ! Tempertaure of the entraining parcel.
real(r8) qtmix(nz,nup)       ! Total water of the entraining parcel.
real(r8) qsmix(nz,nup)       ! Saturated mixing ratio at the tmix.
real(r8) smix(nz,nup)        ! Entropy of the entraining parcel.
real(r8) xsh2o(nz,nup)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(nz,nup)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(nz,nup)   ! Entropy change sue to freezing of precip.
real(r8) dmpdz2d(nz,nup)     ! variable detrainment rate

!+tht
real(r8) zl(nup) ! lcl
!-tht

real(r8) mp(nup)    ! Parcel mass flux.
real(r8) qtp(nup)   ! Parcel total water.
real(r8) sp(nup)    ! Parcel entropy.

real(r8) sp0(nup)    ! Parcel launch entropy.
real(r8) qtp0(nup)   ! Parcel launch total water.
real(r8) mp0(nup)    ! Parcel launch relative mass flux.

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
!real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
!+tht
real(r8) zenv
!-tht
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.
real(r8) dmpdz_lnd, dmpdz_mask

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat   ! Number of iterations for condensation/freezing loop.
integer i,k,ii      ! Loop counters.

real(r8) est
!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.

!+tht should not be necessary but for bit-reproducibility it turns out it is
 !if (.not.tht_tweaks) then 
 ! dmpdz    =-1.e-3_r8    ! Entrainment rate. (-ve for /m)
 ! dmpdz_lnd=-1.e-3_r8 ! idem, on land
 !endif
!-tht

!dmpdpc   = 3.e-2_r8   ! In cloud entrainment rate (/mb).

 lwmax    = 1.e-3_r8    ! Need to put formula in for this.
 tscool   = 0.0_r8   ! Temp at which water loading freezes in the cloud.
!+tht
!lwmax    = 1.e10_r8   ! tht: don't precipitate 
!tscool   =-10._r8     ! tht: allow even just mild supercooling?!
!-tht

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8
!+tht
zenv = 0._r8
!-tht

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

zl(:)=0._r8

! **** Begin loops ****

do k = 1,msg-1
   do i=1,nup 

! Initialize parcel values at launch level.

      if (k == klaunch) then 
         qtp0(i) = q(k)   ! Parcel launch total water (assuming subsaturated) - OK????.

!+tht: formulate dilution on enthalpy not on entropy
         if (tht_tweaks) then
          sp0(i)  = enthalpy(t(k),p(k),qtp0(i),z(k))  ! Parcel launch enthalpy.
         else
          sp0(i)  = entropy (t(k),p(k),qtp0(i))         ! Parcel launch entropy.
         endif
!-tht
         mp0(i)  = 1._r8       ! Parcel launch relative mass (=1 for dmpdp=0 i.e. undilute). 
         smix(k,i)  = sp0(i)
         qtmix(k,i) = qtp0(i)
!+tht: since the function to invert for T is *identical* with sp0(i)=entropy(t), unless there is
! a coding error (likely, given the mess) the result must be t(i,k) (verified 21/2/2014)
         if (tht_tweaks) then
          tmix(k,i) = t(k)
          call qsat_hPa(tmix(k,i),p(k), est, qsmix(k,i))
         else
          tfguess = t(k)
          rcall = 1
          call ientropy (rcall,smix(k,i),p(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
         endif
!-tht
      end if

! Entraining levels
      
      if (k > klaunch) then 

! Set environmental values for this level.                 
         
         dp = (p(k)-p(k-1)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5_r8*(q(k)+q(k-1))         ! Total water of environment.
         tenv  = 0.5_r8*(t(k)+t(k-1)) 
         penv  = 0.5_r8*(p(k)+p(k-1))
!+tht
         zenv  = 0.5_r8*(z(k)+z(k-1))
!-tht

!+tht: base plume dilution on enthalpy not on entropy
         if (tht_tweaks) then
          senv  = enthalpy(tenv,penv,qtenv,zenv) ! Enthalpy of environment.   
         else
          senv  = entropy (tenv,penv,qtenv)      ! Entropy  of environment.   
         endif
!-tht

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*gravit)/(rair*tenv) ! in mb/m since  p in mb.
         dzdp = 1._r8/dpdz                  ! in m/mb
!+tht
! NB: land fudge makes no sense to me - make dmpdz_lnd=dmpdz (as per default code, hard-wired to 1e-3)
        !dmpdp = dmpdz*dzdp
        !dmpdp = dmpdz(i)*dzdp              ! /mb Fractional entrainment 2D
         dmpdp = dmpdz(k,i)*dzdp            ! /mb Fractional entrainment 3D
!-tht

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv 
         qtp(i) = qtp(i) - dmpdp*dp*qtenv 
         mp(i)  = mp(i)  - dmpdp*dp
            
! Entrain s and qt to next level.

         smix(k,i)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(k,i) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

         tfguess = tmix(k-1,i)
         rcall = 2
!+tht
         if (tht_tweaks) then
          call ienthalpy(rcall,smix(k,i),p(k),z(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)   
         else
          call ientropy (rcall,smix(k,i),p(k),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)   
         endif
!-tht

!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.
         if (qsmix(k,i) <= qtmix(k,i) .and. qsmix(k-1,i) > qtmix(k-1,i)) then
            lcl(i) = k
            qxsk   = qtmix(k,i) - qsmix(k,i)
            qxskp1 = qtmix(k-1,i) - qsmix(k-1,i)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(k-1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
!+tht
            zl(i)  = z(k-1) - qxskp1/dqxsdp *dzdp
!-tht
            dsdp   = (smix(k,i)  - smix(k-1,i))/dp
            dqtdp  = (qtmix(k,i) - qtmix(k-1,i))/dp
            slcl   = smix(k-1,i)  +  dsdp* (pl(i)-p(k-1))  
            qtlcl  = qtmix(k-1,i) +  dqtdp*(pl(i)-p(k-1))

            tfguess = tmix(k,i)
            rcall = 3
!+tht
         if (tht_tweaks) then
            call ienthalpy(rcall,slcl,pl(i),zl(i),qtlcl,tl(i),qslcl,tfguess)
         else
            call ientropy (rcall,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)
         endif
!-tht

!            write(iulog,*)' '
!            write(iulog,*)' p',p(i,k+1),pl(i),p(i,lcl(i))
!            write(iulog,*)' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
!            write(iulog,*)' s',smix(i,k+1),slcl,smix(i,lcl(i))
!            write(iulog,*)'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
!            write(iulog,*)'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

         endif
!         
      end if !  k < klaunch

 
   end do ! Levels loop
end do ! Columns loop


!   if ( masterproc ) then
!     do k = 1,msg-1
!         do i = 1,nup
!            write(iulog,*) "after, k, nup, dmpdz ", k, i, dmpdz(k,i)
!         end do
!     end do
!   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

xsh2o = 0._r8
ds_xsh2o = 0._r8
ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do k = 1, msg-1
   do i=1,nup    
      
! Initialize variables at k=klaunch
      
      if (k == klaunch) then
            
! Set parcel values at launch level assume no liquid water.            

         tp(k,i)    = tmix(k,i)
         qstp(k,i)  = q(k) 
         if (tht_tweaks) then
           tpv(k,i)   =  (tp(k,i) + tpert) * (1._r8+qstp(k,i)/epsilo) / (1._r8+qstp(k,i)) !+tht OK with mx ratio
         else
           tpv(k,i)   =  (tp(k,i) + tpert) * (1._r8+1.608_r8*qstp(k,i)) / (1._r8+qstp(k,i))
         endif
         
      end if

      if (k > klaunch) then

         if (tht_tweaks) then           
           smix(k,i)=entropy(tmix(k,i),p(k),qtmix(k,i)) !+tht make sure to use entropy here
         endif
   
!----
! Initiate loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.
! Iterate nit_lheat times for s,qt changes.
         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
            xsh2o(k,i) = max (0._r8, qtmix(k,i) - qsmix(k,i) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)
            ds_xsh2o(k,i) = ds_xsh2o(k-1,i) - cpliq * log (tmix(k,i)/tmelt) * max(0._r8,(xsh2o(k,i)-xsh2o(k-1,i)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
! 
            if (tmix(k,i) <= tmelt+tscool .and. ds_freeze(k-1,i) == 0._r8) then ! One off freezing of condensate. 
               ds_freeze(k,i) = (latice/tmix(k,i)) * max(0._r8,qtmix(k,i)-qsmix(k,i)-xsh2o(k,i)) ! Gain of LH
            end if
            
            if (tmix(k,i) <= tmelt+tscool .and. ds_freeze(k-1,i) /= 0._r8) then ! Continual freezing of additional condensate.
               ds_freeze(k,i) = ds_freeze(k-1,i)+(latice/tmix(k,i)) * max(0._r8,(qsmix(k-1,i)-qsmix(k,i)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).
            new_s = smix(k,i) + ds_xsh2o(k,i) + ds_freeze(k,i) 

! Adjust liquid water and accordingly to xsh2o.
            new_q = qtmix(k,i) - xsh2o(k,i)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(k,i)
            rcall =4
            call ientropy (rcall,new_s, p(k), new_q, tmix(k,i), qsmix(k,i), tfguess)
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(k,i)    = tmix(k,i)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)
         if (new_q > qsmix(k,i)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(k,i) = qsmix(k,i)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(k,i) = new_q
         end if

         if (tht_tweaks) then
           tpv(k,i) = (tp(k,i)+tpert)* (1._r8+qstp(k,i)/epsilo) / (1._r8+ new_q) !+tht
         else
           tpv(k,i) = (tp(k,i)+tpert)* (1._r8+1.608_r8*qstp(k,i)) / (1._r8+ new_q) 
         endif

      end if ! k > klaunch
      
   end do ! Loop for columns
   
end do  ! Loop for vertical levels.


return
end subroutine parcel_dilute

!-----------------------------------------------------------------------------------------
real(r8) function entropy(TK,p,qtot)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(r8), intent(in) :: p,qtot,TK
     real(r8) :: qv,qst,e,est,L
     real(r8), parameter :: pref = 1000._r8

L = latvap - (cpliq - cpwv)*(TK-tmelt)         ! T IN CENTIGRADE

call qsat_hPa(TK, p, est, qst)

qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
e = qv*p / (epsilo +qv)

entropy = (cpair + qtot*cpliq)*log( TK/tmelt) - rair*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qst)

end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
SUBROUTINE ientropy (rcall,s,p,qt,T,qst,Tfg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

  integer, intent(in) :: rcall
  real(r8), intent(in)  :: s, p, Tfg, qt
  real(r8), intent(out) :: qst, T
  real(r8) :: est 
  real(r8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
  integer :: i

  logical :: converged

  ! Max number of iteration loops.
  integer, parameter :: LOOPMAX = 100
  real(r8), parameter :: EPS = 3.e-8_r8

  converged = .false.

  ! Invert the entropy equation -- use Brent's method
  ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

  T = Tfg                  ! Better first guess based on Tprofile from conv.

  a = Tfg-10    !low bracket
  b = Tfg+10    !high bracket

  fa = entropy(a, p, qt) - s
  fb = entropy(b, p, qt) - s

  c=b
  fc=fb
  tol=0.001_r8

  converge: do i=0, LOOPMAX
     if ((fb > 0.0_r8 .and. fc > 0.0_r8) .or. &
          (fb < 0.0_r8 .and. fc < 0.0_r8)) then
        c=a
        fc=fa
        d=b-a
        ebr=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if

     tol1=2.0_r8*EPS*abs(b)+0.5_r8*tol
     xm=0.5_r8*(c-b)
     converged = (abs(xm) <= tol1 .or. fb == 0.0_r8)
     if (converged) exit converge

     if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr=fb/fa
        if (a == c) then
           pbr=2.0_r8*xm*sbr
           qbr=1.0_r8-sbr
        else
           qbr=fa/fc
           rbr=fb/fc
           pbr=sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
           qbr=(qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
        end if
        if (pbr > 0.0_r8) qbr=-qbr
        pbr=abs(pbr)
        if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
           ebr=d
           d=pbr/qbr
        else
           d=xm
           ebr=d
        end if
     else
        d=xm
        ebr=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

     fb = entropy(b, p, qt) - s

  end do converge

  T = b
  call qsat_hPa(T, p, est, qst)

  if (.not. converged) then
     call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

end SUBROUTINE ientropy

! Wrapper for qsat_water that does translation between Pa and hPa
! qsat_water uses Pa internally, so get it right, need to pass in Pa.
! Afterward, set es back to hPa.
subroutine qsat_hPa(t, p, es, qm)
  use wv_saturation, only: qsat_water

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature (K)
  real(r8), intent(in) :: p    ! Pressure (hPa)
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure (hPa)
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

  call qsat_water(t, p*100._r8, es, qm)

  es = es*0.01_r8

end subroutine qsat_hPa

!-----------------------------------------------------------------------------------------
real(r8) function enthalpy(TK,p,qtot,z)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
!
     real(r8), intent(in) :: p,qtot,TK,z
     real(r8) :: qv,qst,e,est,L

L = latvap - (cpliq - cpwv)*(TK-tmelt)

call qsat_hPa(TK, p, est, qst)
qv = min(qtot,qst)                         ! Partition qtot into vapor part only.

!enthalpy = (cpres + qtot*cpliq)*(TK-tfreez) + L*qv + (1._r8+qtot)*grav*z
 enthalpy = (cpair + qtot*cpliq)* TK         + L*qv + (1._r8+qtot)*gravit*z
 
return
end FUNCTION enthalpy

!-----------------------------------------------------------------------------------------
 SUBROUTINE ienthalpy (rcall,s,p,z,qt,T,qst,Tfg) !identical with iENTROPY, only function calls swapped
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

  integer, intent(in) :: rcall
  real(r8), intent(in)  :: s, p, z, Tfg, qt
  real(r8), intent(out) :: qst, T
  real(r8) :: est
  real(r8) :: a,b,c,d,ebr,fa,fb,fc,pbr,qbr,rbr,sbr,tol1,xm,tol
  integer :: i

  logical :: converged

  ! Max number of iteration loops.
  integer, parameter :: LOOPMAX = 100
  real(r8), parameter :: EPS = 3.e-8_r8

  converged = .false.

  ! Invert the entropy equation -- use Brent's method
  ! Brent, R. P. Ch. 3-4 in Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973.

  T = Tfg                  ! Better first guess based on Tprofile from conv.

  a = Tfg-10    !low bracket
  b = Tfg+10    !high bracket

  fa = enthalpy(a, p, qt,z) - s
  fb = enthalpy(b, p, qt,z) - s

  c=b
  fc=fb
  tol=0.001_r8

  converge: do i=0, LOOPMAX
     if ((fb > 0.0_r8 .and. fc > 0.0_r8) .or. &
          (fb < 0.0_r8 .and. fc < 0.0_r8)) then
        c=a
        fc=fa
        d=b-a
        ebr=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if

     tol1=2.0_r8*EPS*abs(b)+0.5_r8*tol
     xm=0.5_r8*(c-b)
     converged = (abs(xm) <= tol1 .or. fb == 0.0_r8)
     if (converged) exit converge

     if (abs(ebr) >= tol1 .and. abs(fa) > abs(fb)) then
        sbr=fb/fa
        if (a == c) then
           pbr=2.0_r8*xm*sbr
           qbr=1.0_r8-sbr
        else
           qbr=fa/fc
           rbr=fb/fc
           pbr=sbr*(2.0_r8*xm*qbr*(qbr-rbr)-(b-a)*(rbr-1.0_r8))
           qbr=(qbr-1.0_r8)*(rbr-1.0_r8)*(sbr-1.0_r8)
        end if
        if (pbr > 0.0_r8) qbr=-qbr
        pbr=abs(pbr)
        if (2.0_r8*pbr  <  min(3.0_r8*xm*qbr-abs(tol1*qbr),abs(ebr*qbr))) then
           ebr=d
           d=pbr/qbr
        else
           d=xm
           ebr=d
        end if
     else
        d=xm
        ebr=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )

     fb = enthalpy(b, p, qt,z) - s

  end do converge

  T = b
  call qsat_hPa(T, p, est, qst)

  if (.not. converged) then
     call endrun('**** ZM_CONV IENTHALPY: Tmix did not converge ****')
  end if

100 format (A,I1,I4,I4,7(A,F6.2))

 end SUBROUTINE ienthalpy

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module clubb_mf
