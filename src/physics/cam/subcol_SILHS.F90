module subcol_SILHS
   !---------------------------------------------------------------------------
   ! Purpose:
   !
   ! Implement a subcolumn scheme based on the Subgrid Importance Latin Hypercube 
   ! Sampling (SILHS) functionality of the CLUBB moist turbulence parameterization.
   !
   !---------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8, r4=>shr_kind_r4, i4=>shr_kind_i4
   use physics_types,    only: physics_state, physics_tend, physics_ptend
   use ppgrid,           only: pcols, psubcols, pver, pverp
   use constituents,     only: pcnst, cnst_get_ind
   use cam_abortutils,   only: endrun
   use cam_logfile,      only: iulog
   use cam_history,      only: addfld, add_default, outfld, horiz_only
#ifdef CLUBB_SGS
#ifdef SILHS
   use clubb_intr,       only: pdf_params_chnk
   use clubb_api_module, only: &
        hmp2_ip_on_hmm2_ip_slope_type, &
        hmp2_ip_on_hmm2_ip_intrcpt_type

   use silhs_api_module, only: &
        silhs_config_flags_type
#endif
#endif
   use physconst,     only: cpair, gravit, latvap, latice, rair

   implicit none
   private
   save

   public :: subcol_register_SILHS  ! 
   public :: subcol_init_SILHS      ! Initialize 
   public :: subcol_gen_SILHS       ! Generate subcolumn fields by calling SILHS 
   public :: subcol_readnl_SILHS    ! SILHS namelist reader
   public :: subcol_ptend_avg_SILHS
   public :: subcol_SILHS_var_covar_driver
   public :: subcol_SILHS_fill_holes_conserv
   public :: subcol_SILHS_hydromet_conc_tend_lim
   private :: fill_holes_sedimentation
   private :: fill_holes_same_phase_vert
#ifdef SILHS
   private :: Abs_Temp_profile
   private :: StaticEng_profile
   ! Calc subcol mean ! Calc subcol variance
   private :: meansc
   private :: stdsc
#endif

   !-----
   ! Private module vars
   !-----

   ! constituent indicies
   integer :: &
      ixq      = 0, &
      ixcldliq = 0, &
      ixnumliq = 0, &
      ixcldice = 0, &
      ixnumice = 0, &
      ixrain   = 0, &
      ixnumrain= 0, &
      ixsnow   = 0, &
      ixnumsnow= 0
   
   ! Pbuf indicies
   integer :: thlm_idx, rcm_idx, rtm_idx, ice_supersat_idx, &
              alst_idx, cld_idx, qrain_idx, qsnow_idx, &
              nrain_idx, nsnow_idx, ztodt_idx, tke_idx, kvh_idx, &
              prec_pcw_idx, snow_pcw_idx, prec_str_idx, snow_str_idx, &
              qcsedten_idx, qrsedten_idx, qisedten_idx, qssedten_idx, &
              vtrmc_idx, umr_idx, vtrmi_idx, ums_idx, qcsevap_idx, qisevap_idx

   logical :: subcol_SILHS_weight  ! if set, sets up weights for averaging subcolumns for SILHS
   integer :: subcol_SILHS_numsubcol ! number of subcolumns for this run
   logical :: docldfracscaling = .false. ! Weight tendencies by cloud fraction

   character(len=256) :: subcol_SILHS_corr_file_path
   character(len=16)  :: subcol_SILHS_corr_file_name

   logical :: subcol_SILHS_q_to_micro, &
              subcol_SILHS_n_to_micro, &
              subcol_SILHS_use_clear_col, &
              subcol_SILHS_meanice, &
              subcol_SILHS_constrainmn

   logical :: subcol_SILHS_var_covar_src

   real(r8) :: subcol_SILHS_ncnp2_on_ncnm2

   ! There may or may not be a better place to put this.
   real(r8), parameter :: p0_clubb = 100000._r8


!   real(r8) :: subcol_SILHS_c6rt, subcol_SILHS_c7, subcol_SILHS_c8, subcol_SILHS_c11, &
!               subcol_SILHS_c11b, subcol_SILHS_gamma_coef, &
!               subcol_SILHS_mult_coef, subcol_SILHS_mu

   real(r8) :: ztodt  ! model timestep
#ifdef CLUBB_SGS
#ifdef SILHS
    type(hmp2_ip_on_hmm2_ip_slope_type) :: subcol_SILHS_hmp2_ip_on_hmm2_ip_slope    
    type(hmp2_ip_on_hmm2_ip_intrcpt_type) :: subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt

    type(silhs_config_flags_type) :: silhs_config_flags
#endif
#endif

contains

   subroutine subcol_register_SILHS()

      !--------------------------------
      ! Register fields needed by SILHS in the physics buffer
      ! Currently, most fields needed by SILHS but calculated by CLUBB are registered
      ! by clubb in clubb_intr.F90.
      ! 
      !--------------------------------
#ifdef CLUBB_SGS
#ifdef SILHS
#endif
#endif
   end subroutine subcol_register_SILHS


   subroutine subcol_readnl_SILHS(nlfile)
#ifdef CLUBB_SGS
#ifdef SILHS
      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: masterproc, masterprocid, mpicom
      use spmd_utils,      only: mpi_integer, mpi_logical, mpi_character, mpir8
      use clubb_api_module,only: core_rknd
#endif
#endif
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
#ifdef CLUBB_SGS
#ifdef SILHS
      namelist /subcol_SILHS_nl/ subcol_SILHS_weight, &
                                 subcol_SILHS_numsubcol, &
                                 subcol_SILHS_corr_file_path, &
                                 subcol_SILHS_corr_file_name, &
                                 subcol_SILHS_q_to_micro, &
                                 subcol_SILHS_n_to_micro, &
                                 subcol_SILHS_ncnp2_on_ncnm2, &
                                 subcol_SILHS_hmp2_ip_on_hmm2_ip_slope, &
                                 subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt, &
                                 subcol_SILHS_meanice, &
                                 subcol_SILHS_use_clear_col, &
                                 subcol_SILHS_constrainmn, &
                                 subcol_SILHS_var_covar_src
!                                 subcol_SILHS_c6rt, subcol_SILHS_c7, &
!                                 subcol_SILHS_c8, subcol_SILHS_c11, subcol_SILHS_c11b, &
!                                 subcol_SILHS_gamma_coef, subcol_SILHS_mult_coef, subcol_SILHS_mu

      !-----------------------------------------------------------------------------
      ! Set defaults

      ! Eric Raut changed a default.
      subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%Ni = 0.0_core_rknd
      subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%Ni = 0.5_core_rknd

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'subcol_SILHS_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, subcol_SILHS_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun('subcol_readnl_SILHS: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpi_bcast(subcol_SILHS_weight,    1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_numsubcol, 1, mpi_integer, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_path, len(subcol_SILHS_corr_file_path), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_corr_file_name, len(subcol_SILHS_corr_file_name), &
                     mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_use_clear_col, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_constrainmn, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_meanice, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_q_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_n_to_micro, 1, mpi_logical, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_var_covar_src,1,mpi_logical,masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_ncnp2_on_ncnm2, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%rr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%Nr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%ri, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%Ni, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%rs, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_slope%Ns, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%rr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%Nr, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%ri, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%Ni, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%rs, 1, mpir8, masterprocid, mpicom, ierr)
      call mpi_bcast(subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt%Ns, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c6rt, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c7, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c8, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_c11b, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_gamma_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mult_coef, 1, mpir8, masterprocid, mpicom, ierr)
!      call mpi_bcast(subcol_SILHS_mu, 1, mpir8, masterprocid, mpicom, ierr)

! SPMD
#endif
! SILHS
#endif
! CLUBB_SGS
#endif
   end subroutine subcol_readnl_SILHS


   subroutine subcol_init_SILHS(pbuf2d)

      !--------------------------------
      ! Read in parameters and initialize SILHS PDF fields.
      ! Set up indexes into Pbuf fields.
      ! Register history outputs.
      !--------------------------------

      use physics_buffer,          only: physics_buffer_desc, pbuf_get_field, &
                                         dtype_r8, pbuf_get_index
      use units,                   only: getunit, freeunit 
#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,        only: core_rknd, &
                                         pdf_dim, &
                                         setup_corr_varnce_array_api, &
                                         init_pdf_hydromet_arrays_api, &
                                         Ncnp2_on_Ncnm2, &
                                         set_clubb_debug_level_api

      use silhs_api_module,        only: set_default_silhs_config_flags_api, &
                                         initialize_silhs_config_flags_type_api, &
                                         print_silhs_config_flags_api

      use spmd_utils,              only: iam

      use clubb_intr,              only: init_clubb_config_flags, &
                                         clubb_config_flags

#endif
#endif

      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CLUBB_SGS
#ifdef SILHS

      integer :: iunit = 501 ! Default value, will get iunit from CAM 
      !character(len=*), parameter :: default_corr_case = "arm_97"
      character(len=*), parameter :: &
            cloud_file_ext  = "_corr_array_cloud.in", & ! File extensions for corr files
            below_file_ext  = "_corr_array_below.in"
      character(len=256) :: corr_file_path_cloud, corr_file_path_below

      ! To set up CLUBB hydromet indices
      integer :: &
          hydromet_dim, & ! Number of enabled hydrometeors
          iirr,         & ! Hydrometeor array index for rain water mixing ratio, rr
          iirs,         & ! Hydrometeor array index for snow mixing ratio, rs
          iiri,         & ! Hydrometeor array index for ice mixing ratio, ri
          iirg,         & ! Hydrometeor array index for graupel mixing ratio, rg
          iiNr,         & ! Hydrometeor array index for rain drop concentration, Nr
          iiNs,         & ! Hydrometeor array index for snow concentration, Ns
          iiNi,         & ! Hydrometeor array index for ice concentration, Ni
          iiNg            ! Hydrometeor array index for graupel concentration, Ng

      integer :: &
          cluster_allocation_strategy

      logical :: &
          l_lh_importance_sampling, &
          l_Lscale_vert_avg, &
          l_lh_straight_mc, &
          l_lh_clustered_sampling, &
          l_rcm_in_cloud_k_lh_start, &
          l_random_k_lh_start, &
          l_max_overlap_in_cloud, &
          l_lh_instant_var_covar_src, &
          l_lh_limit_weights, &
          l_lh_var_frac, &
          l_lh_normalize_weights


      ! Set CLUBB's debug level
      ! This is called in module clubb_intr; no need to do it here.
!      call set_clubb_debug_level_api( 0 )

      !-------------------------------
      ! CLUBB-SILHS Parameters (global module variables)
      !-------------------------------
      call set_default_silhs_config_flags_api( cluster_allocation_strategy, &
                                               l_lh_importance_sampling, &
                                               l_Lscale_vert_avg, &
                                               l_lh_straight_mc, &
                                               l_lh_clustered_sampling, &
                                               l_rcm_in_cloud_k_lh_start, &
                                               l_random_k_lh_start, &
                                               l_max_overlap_in_cloud, &
                                               l_lh_instant_var_covar_src, &
                                               l_lh_limit_weights, &
                                               l_lh_var_frac, &
                                               l_lh_normalize_weights )

      call init_clubb_config_flags( clubb_config_flags ) ! In/Out
      clubb_config_flags%l_fix_w_chi_eta_correlations = .true.
      l_lh_importance_sampling = .true.
      clubb_config_flags%l_diagnose_correlations = .false.
      clubb_config_flags%l_calc_w_corr = .false.
!      l_prescribed_avg_deltaz = .false.
      clubb_config_flags%l_use_cloud_cover = .false.
      clubb_config_flags%l_const_Nc_in_cloud = .true.

      call initialize_silhs_config_flags_type_api( cluster_allocation_strategy, &
                                                   l_lh_importance_sampling, &
                                                   l_Lscale_vert_avg, &
                                                   l_lh_straight_mc, &
                                                   l_lh_clustered_sampling, &
                                                   l_rcm_in_cloud_k_lh_start, &
                                                   l_random_k_lh_start, &
                                                   l_max_overlap_in_cloud, &
                                                   l_lh_instant_var_covar_src, &
                                                   l_lh_limit_weights, &
                                                   l_lh_var_frac, &
                                                   l_lh_normalize_weights, &
                                                   silhs_config_flags )

      ! Print the SILHS configurable flags
      call print_silhs_config_flags_api( iulog, silhs_config_flags ) ! Intent(in)

      ! Values from the namelist
      docldfracscaling = subcol_SILHS_use_clear_col

      ! Namelist "tuning" or set correlations
      ! KTC Todo: Move these to a tuning "in" file or into the namelist
      ! JHTODO: we might want these on CLUBB's API and ultimatively on a namelist for tuning
!      C6rt = subcol_SILHS_c6rt
!      C7 = subcol_SILHS_c7                                      ! to all ice clouds
!      C8 = subcol_SILHS_c8
!      C11 = subcol_SILHS_c11
!      C11b = subcol_SILHS_c11b
!      gamma_coef = subcol_SILHS_gamma_coef
!      mult_coef = subcol_SILHS_mult_coef
!      mu = subcol_SILHS_mu

      !call set_clubb_debug_level( 0 )  !#KTCtodo: Add a namelist variable to set debug level
     

      ! Get constituent indices
      call cnst_get_ind('Q', ixq)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('NUMLIQ', ixnumliq)
      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('NUMICE', ixnumice)
      call cnst_get_ind('RAINQM', ixrain, abort=.false.)
      call cnst_get_ind('NUMRAI', ixnumrain, abort=.false.)
      call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)
      call cnst_get_ind('NUMSNO', ixnumsnow, abort=.false.)

      ! Get physics buffer indexes
      thlm_idx = pbuf_get_index('THLM')
      rcm_idx = pbuf_get_index('RCM')
      rtm_idx = pbuf_get_index('RTM')
      cld_idx = pbuf_get_index('CLD')
      alst_idx = pbuf_get_index('ALST')  ! SILHS expects clubb's cloud_frac liq stratus fraction
      ztodt_idx = pbuf_get_index('ZTODT')
      ice_supersat_idx = pbuf_get_index('ISS_FRAC')
      tke_idx = pbuf_get_index('tke')
      kvh_idx = pbuf_get_index('kvh')
      prec_pcw_idx = pbuf_get_index('PREC_PCW')
      snow_pcw_idx = pbuf_get_index('SNOW_PCW')
      prec_str_idx = pbuf_get_index('PREC_STR')
      snow_str_idx = pbuf_get_index('SNOW_STR')
      qcsedten_idx = pbuf_get_index('QCSEDTEN')
      qrsedten_idx = pbuf_get_index('QRSEDTEN')
      qisedten_idx = pbuf_get_index('QISEDTEN')
      qssedten_idx = pbuf_get_index('QSSEDTEN')
      vtrmc_idx = pbuf_get_index('VTRMC')
      umr_idx = pbuf_get_index('UMR')
      vtrmi_idx = pbuf_get_index('VTRMI')
      ums_idx = pbuf_get_index('UMS')
      qcsevap_idx = pbuf_get_index('QCSEVAP')
      qisevap_idx = pbuf_get_index('QISEVAP')
      qrain_idx = pbuf_get_index('QRAIN')
      qsnow_idx = pbuf_get_index('QSNOW')
      nrain_idx = pbuf_get_index('NRAIN')
      nsnow_idx = pbuf_get_index('NSNOW')
     
      !-------------------------------
      ! Set up SILHS hydrometeors #KTCtodo: move microphys specification to config time,
      !        Steve wants to set up a microphysics query so I can ask the microphysics
      !        scheme which hydrometeors to use. For the future.
      !-------------------------------
      iirr = 1
      iirs = 3
      iiri  = 5
      iirg = -1

      iiNr    = 2
      iiNs = 4
      iiNi    = 6
      iiNg = -1

      hydromet_dim = 6

 
      ! Set up pdf indices, hydromet indicies, hydromet arrays, and hydromet variance ratios
      call init_pdf_hydromet_arrays_api( 1.0_core_rknd, 1.0_core_rknd,  & ! intent(in)
                                         hydromet_dim,                  & ! intent(in)
                                         iirr, iiri, iirs, iirg,        & ! intent(in)
                                         iiNr, iiNi, iiNs, iiNg,        & ! intent(in)
                                         subcol_SILHS_hmp2_ip_on_hmm2_ip_slope,      & ! optional(in)
                                         subcol_SILHS_hmp2_ip_on_hmm2_ip_intrcpt )     ! optional(in)

      Ncnp2_on_Ncnm2 = subcol_SILHS_ncnp2_on_ncnm2

      !-------------------------------
      ! Set up hydrometeors and correlation arrays for SILHS
      !-------------------------------
      corr_file_path_cloud = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//cloud_file_ext
      corr_file_path_below = trim( subcol_SILHS_corr_file_path )//trim( subcol_SILHS_corr_file_name )//below_file_ext

      iunit = getunit()


      call setup_corr_varnce_array_api( corr_file_path_cloud, corr_file_path_below, &
                                        iunit, &
                                        clubb_config_flags%l_fix_w_chi_eta_correlations )
      call freeunit(iunit) 

      !-------------------------------
      ! Register output fields from SILHS
      ! #KTCtodo: Remove these from the default output list
      !-------------------------------
      call addfld('SILHS_NCLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Cloud Number Concentration', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm^-3', &
           'Subcolumn Number Concentration of Rain from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_OMEGA_SCOL', (/'psubcols', 'ilev    '/), 'I', 'Pa/s', &
           'Subcolumn vertical pressure velocity', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RCM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_NICLD_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Cloud Ice Number Conc from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RRAIN_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg', &
           'Subcolumn Precipitating Liquid Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_RT_SCOL', (/'psubcols', 'ilev    '/), 'I', 'kg/kg ', &
           'Subcolumn Total Water from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_THLM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'K', &
           'Subcolumn liquid water pot temperature', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WEIGHT_SCOL', (/'psubcols'/), 'I', 'frac', &
           'Weights for each subcolumn', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('SILHS_WM_SCOL', (/'psubcols', 'ilev    '/), 'I', 'm/s', &
           'Subcolumn vertical velocity from SILHS', flag_xyfill=.true., fill_value=1.e30_r8)

      call addfld('NR_IN_LH', (/ 'lev' /), 'I', 'm^-3', &
                  'Num Rain Conc as input to SILHS')
     call addfld('RTM_CLUBB', (/ 'ilev' /), 'I', 'kg/kg', &
                  'Input total water mixing ratio')
     call addfld('THLM_CLUBB', (/ 'ilev' /), 'I', 'K', &
                  'Input liquid water potential temperature')
     call addfld('SILHS_QC_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud water mixing ratio')
     call addfld('SILHS_QI_IN', (/ 'lev' /), 'I', 'kg/kg', &
                  'Input cloud ice mixing ratio')
     call addfld('SILHS_NC_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud water number concentration')
     call addfld('SILHS_NI_IN', (/ 'lev' /), 'I', '#/kg', &
                  'Input cloud ice number concentration')
     call addfld('AKM_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Exact Kessler autoconversion')
     call addfld('AKM_LH_CLUBB', (/ 'ilev' /), 'I', '(kg/kg)/s', &
                  'Monte Carlo estimate of Kessler autoconversion')
     call addfld('INVS_EXNER', (/ 'lev' /), 'I', 'none', &
                  'inverse EXNER function from state in subcol_SILHS')
     call addfld('SILHS_ZTODT', horiz_only, 'I', 's', & 
                  'Length of Physics timestep (for debugging)')
     if ( subcol_SILHS_constrainmn ) then
        call addfld('SILHS_MSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                    'Mean Cloud Ice across subcolumns')
        call addfld('SILHS_STDSC_CLDICE', (/ 'lev' /), 'A', 'kg/kg', &
                    'Standard deviation of Ice across subcolumns')
        if ( ixsnow > 0 ) then
           call addfld('SILHS_MSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                       'Mean Cloud Liquid across subcolumns')
           call addfld('SILHS_STDSC_CLDLIQ', (/ 'lev' /), 'A', 'kg/kg', &
                       'Standard deviation of Liquid across subcolumns')
           call addfld('SILHS_MSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                       'Mean water vapor across subcolumns')
           call addfld('SILHS_STDSC_Q', (/ 'lev' /), 'A', 'kg/kg', &
                       'Standard deviation of water vapor across subcolumns')
        endif ! ixsnow > 0
     endif ! subcol_SILHS_constrainmn
     call addfld('SILHS_EFF_CLDFRAC', (/ 'lev' /), 'A', 'frac', &
                  'Calculated cloud fraction from subcolumn liq or ice') 

     call addfld('SILHS_CLUBB_PRECIP_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Precipitation fraction from CLUBB (set_up_pdf_params_incl_hydromet)')
     call addfld('SILHS_CLUBB_ICE_SS_FRAC', (/ 'lev' /), 'A', 'frac', &
                  'Ice supersaturation fraction from CLUBB')

     call addfld ('QVHFTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Water vapor mixing ratio tendency from hole filling')
     call addfld ('QCHFTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud water mixing ratio tendency from hole filling')
     call addfld ('QRHFTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Rain water mixing ratio tendency from hole filling')
     call addfld ('QIHFTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud ice mixing ratio tendency from hole filling')
     call addfld ('QSHFTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Snow mixing ratio tendency from hole filling')
     call addfld ('THFTEN', (/ 'lev' /), 'A', 'K/s', 'Temperature tendency from hole filling')

#endif
#endif
   end subroutine subcol_init_SILHS
   
   subroutine subcol_gen_SILHS(state, tend, state_sc, tend_sc, pbuf)
      !-------------------------------
      ! This is where the subcolumns are created, and the call to
      !      generate_silhs_sample_mod_api
      !    goes out. Variables needed to make this call are pulled from the 
      !    pbuf, from module data, and calculated based on the CAM state.
      !-------------------------------

      use physics_buffer,         only : physics_buffer_desc, pbuf_get_index, &
                                         pbuf_get_field
      use ppgrid,                 only : pver, pverp, pcols
      use ref_pres,               only : top_lev => trop_cloud_top_lev
      use time_manager,           only : get_nstep
      use subcol_utils,           only : subcol_set_subcols, subcol_set_weight
      use phys_control,           only : phys_getopts
      use spmd_utils,             only : masterproc
      use shr_const_mod,          only : SHR_CONST_PI, SHR_CONST_RHOFW

#ifdef CLUBB_SGS
#ifdef SILHS
      use clubb_api_module,       only : hydromet_dim, &

                                         setup_pdf_parameters_api, &

                                         l_stats_samp, &

                                         hydromet_pdf_parameter, &

                                         zm2zt_api, setup_grid_heights_api, gr, &

                                         iirr, iiNr, iirs, iiri, &
                                         iirg, iiNs, &
                                         iiNi, iiNg, &

                                         core_rknd, &

                                         w_tol_sqd, zero_threshold, &
                                         em_min, cloud_frac_min, & ! rc_tol, &

                                         pdf_dim, &
                                         corr_array_n_cloud, &
                                         corr_array_n_below, &
                                         iiPDF_chi, iiPDF_rr, &
                                         iiPDF_w, iiPDF_Nr, &
                                         iiPDF_ri, iiPDF_Ni, &
                                         iiPDF_Ncn, iiPDF_rs, iiPDF_Ns, &

                                         genrand_intg, genrand_init_api, &

                                         nparams, ic_K, &
                                         read_parameters_api
   
      use silhs_api_module, only :       generate_silhs_sample_api, & ! Ncn_to_Nc, &
                                         clip_transform_silhs_output_api, &
                                         est_kessler_microphys_api

      use clubb_intr, only:              clubb_config_flags
#endif
#endif
      
      ! CAM data structures
      type(physics_state), intent(inout) :: state
      type(physics_tend),  intent(inout) :: tend
      type(physics_state), intent(inout) :: state_sc        ! sub-column state
      type(physics_tend),  intent(inout) :: tend_sc         ! sub-column tend
      type(physics_buffer_desc), pointer :: pbuf(:)

#ifdef CLUBB_SGS
#ifdef SILHS
      !----------------
      ! Local variables
      !----------------
      logical, parameter :: &
                 l_implemented = .true.   ! Implemented in a host model
      logical, parameter :: rx_Nc = .false. ! Use NC calculated based on grid mean effective radius
      integer, parameter :: &
                 grid_type = 3            ! The 3 grid centered on momentum levels
      real(r8), parameter :: cldmin = 0.001_r8 ! To use when cld frac = 0.0 to be consistant with micro_mg
      real(r8), parameter :: min_num_conc = 1.0e-12_r8
      real(r8), parameter :: qsmall = 1.0e-18_r8  ! Microphysics cut-off for cloud

      integer :: i, j, k, ngrdcol, ncol, lchnk, stncol
      integer :: begin_height, end_height ! Output from setup_grid call
      real(r8) :: sfc_elevation  ! Surface elevation
      real(r8), dimension(pverp-top_lev+1) :: zt_g, zi_g ! Thermo & Momentum grids for clubb
      real(r8), dimension(pverp) :: scfrac     ! cloud fraction based on sc distributions
      real(r8) :: msc, std, maxcldfrac, maxsccldfrac
      real(r8) :: scale = 1.0_r8

      real(r8), dimension(nparams) :: clubb_params ! Adjustable CLUBB parameters

      real(r8) :: c_K ! CLUBB parameter c_K (for eddy diffusivity)

      integer( kind = genrand_intg ) :: &
        lh_seed    ! Seed used in random number generator that will be different
                   ! for each column, yet reproducible for a restart run

      !----------------
      ! Required for set_up_pdf_params_incl_hydromet
      !----------------
      real(r8), dimension(pverp-top_lev+1) :: cld_frac_in  ! Cloud fraction
      type(hydromet_pdf_parameter), dimension(pverp-top_lev+1) :: &
                                    hydromet_pdf_params  ! Hydrometeor PDF parameters
      real(r8), dimension(:,:,:), allocatable :: &       ! Correlation matrix for pdf components
                                    corr_array_1, corr_array_2 
      real(r8), dimension(:,:), allocatable :: &
                                    mu_x_1, mu_x_2, &    ! Mean array for PDF components
                                    sigma_x_1, sigma_x_2 ! Std dev arr for PDF components
      real(r8), dimension(:,:,:), allocatable :: &       ! Transposed corr cholesky mtx
                                    corr_cholesky_mtx_1, corr_cholesky_mtx_2
      real(r8), dimension(pverp-top_lev+1) :: Nc_in_cloud
      real(r8), dimension(pverp-top_lev+1) :: ice_supersat_frac_in
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: hydrometp2


      !----------------
      ! Input to generate_silhs_sample
      !----------------
      integer :: iter                            ! CLUBB iteration 
      integer :: num_subcols                     ! Number of subcolumns
      integer, dimension(pcols) :: numsubcol_arr ! To set up the state struct
      integer, parameter :: sequence_length = 1  ! Number of timesteps btn subcol calls
      real(r8), dimension(pverp-top_lev+1) :: rho_ds_zt    ! Dry static density (kg/m^3) on thermo levs
      real(r8), dimension(pver)  :: dz_g         ! thickness of layer
      real(r8), dimension(pverp-top_lev+1) :: delta_zm     ! Difference in u wind altitudes
      real(r8), dimension(pverp-top_lev+1) :: invs_dzm     ! 1/delta_zm
      real(r8), dimension(pverp-top_lev+1) :: rcm_in       ! Cld water mixing ratio on CLUBB levs
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: hydromet  ! Hydrometeor species
      real(r8), dimension(pverp-top_lev+1,hydromet_dim) :: wphydrometp  ! Hydrometeor flux
      real(r8), dimension(pverp-top_lev+1)              :: Ncm ! Mean cloud droplet concentration, <N_c>

      real(r8), dimension(pverp-top_lev+1) :: tke       ! TKE
      real(r8), dimension(pverp-top_lev+1) :: khzm      ! Eddy diffusivity coef
      real(r8), dimension(pverp-top_lev+1) :: Lscale_zm ! CLUBB's length scale on momentum (zm) levels
      real(r8), dimension(pverp-top_lev+1) :: Lscale    ! CLUBB's length scale

      logical, parameter :: &  
         l_calc_weights_all_levs = .false. ! .false. if all time steps use the same
                                          !   weights at all vertical grid levels 
      logical :: & 
        l_calc_weights_all_levs_itime, & ! .true. if we calculate sample weights separately at all 
                                         !    grid levels at the current time step   
        l_rad_itime                      ! .true. if we calculate radiation at the current time step  
      
      !---------------
      !Output from generate_silhs_sample
      !--------------
      real(r8), allocatable, dimension(:,:,:) :: X_nl_all_levs ! Sample transformed to normal-lognormal
      real(r8), allocatable, dimension(:,:)   :: lh_sample_point_weights ! Subcolumn weights
      integer,  allocatable, dimension(:,:)    :: X_mixt_comp_all_levs ! Which Mixture Component

      real(r8), allocatable, dimension(:,:) :: rc_all_points ! Calculate RCM from LH output
      real(r8), allocatable, dimension(:,:) :: rain_all_pts  ! Calculate Rain from LH output
      real(r8), allocatable, dimension(:,:) :: nrain_all_pts ! Calculate Rain Conc from LH
      real(r8), allocatable, dimension(:,:) :: snow_all_pts  ! Calculate Snow from LH output
      real(r8), allocatable, dimension(:,:) :: nsnow_all_pts ! Calculate Snow Conc from LH
      real(r8), allocatable, dimension(:,:) :: w_all_points  ! Calculate W from LH output
      ! real(r8), allocatable, dimension(:,:) :: RVM_lh_out    ! Vapor mixing ratio sent away
      real(r8), allocatable, dimension(:,:) :: ice_all_pts   ! Calculate Cld Ice from LH output
      real(r8), allocatable, dimension(:,:) :: nice_all_pts  ! Calculate Num cld ice from LH
      real(r8), allocatable, dimension(:,:) :: nclw_all_pts  ! Calculate Num cld wat from LH

      !----------------
      ! Output from clip_transform_silhs_output_api
      !----------------
      real( kind = core_rknd ), dimension(:,:), allocatable :: &
        lh_rt_clipped,  & ! rt generated from silhs sample points
        lh_thl_clipped, & ! thl generated from silhs sample points
        lh_rc_clipped,  & ! rc generated from silhs sample points
        lh_rv_clipped,  & ! rv generated from silhs sample points
        lh_Nc_clipped     ! Nc generated from silhs sample points

      logical, parameter :: &
        l_use_Ncn_to_Nc = .true.  ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                                  ! Ncn_to_Nc might cause problems with the MG microphysics 
                                  ! since the changes made here (Nc-tendency) are not fed into 
                                  ! the microphysics
        

      !----------------
      ! Output to history
      !----------------
      ! V. Larson note: These variables are on the zt (full) levels: why do they
      ! have dimension pverp?  The pverp level corresponds to the CLUBB
      ! below-ground level.
      ! The variables in this paragraph are oriented like CAM variables (k=1 is
      ! the model top).
      ! They are flipped versions of CLUBB variables, for the entire chunk.
      real(r8), dimension(pcols*psubcols, pverp) :: RT_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: THL_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: OMEGA_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: WM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RVM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RCM_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NCLW_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: ICE_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NICE_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: RAIN_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NRAIN_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: SNOW_lh_out
      real(r8), dimension(pcols*psubcols, pverp) :: NSNOW_lh_out

      real(r8), dimension(state_sc%psetcols) :: weights ! Subcol weights

      real(r8), dimension(pcols, pver) :: meansc_ice
      real(r8), dimension(pcols, pver) :: stdsc_ice

      real(r8), dimension(pcols, pver) :: meansc_liq
      real(r8), dimension(pcols, pver) :: stdsc_liq

      real(r8), dimension(pcols, pver) :: meansc_vap
      real(r8), dimension(pcols, pver) :: stdsc_vap
      real(r8), dimension(pcols, pver) :: grmn_eff_rad
      real(r8), dimension(pcols, pver) :: eff_cldfrac
      real(r8), dimension(pcols, pver) :: precip_frac_out

      real(r8) :: tmp_mean, diff_mean, rcubed

      !----------------
      ! Output from Est_Kessler_microphys
      !----------------
      real(r8), dimension(pverp-top_lev+1) :: lh_Akm     ! Monte Carlo estimate of Kessler Autoconversion
      real(r8), dimension(pverp-top_lev+1) :: AKm        ! Exact Kessler autoconversion
      real(r8), dimension(pverp-top_lev+1) :: AKstd      ! Exact Stdev of gba Kessler
      real(r8), dimension(pverp-top_lev+1) :: AKstd_cld  ! Exact w/in cloud stdev of gba Kessler
      real(r8), dimension(pverp-top_lev+1) :: AKm_rcm    ! Exact local gba Kessler auto based on rcm
      real(r8), dimension(pverp-top_lev+1) :: AKm_rcc    ! Exact local gba Kessler based on w/in cloud rc
      real(r8), dimension(pverp-top_lev+1) :: lh_rcm_avg ! LH estimate of grid box avg liquid water
      real(r8), dimension(pcols,pverp) :: lh_AKm_out, AKm_out

      !----------------
      ! Needed to update State
      !----------------
      real(r8), dimension(pver)  :: Temp_prof  ! Subcolumn LWPT converted to Abs Temp
      real(r8), dimension(pver)  :: SE_prof    ! Static Energy calculated from Abs Temp
      real(r8), dimension(pver)  :: No_cloud = 0.0_r8     ! Clear air condensate profile
      real(r8), dimension(pcols, pver)  :: invs_exner  ! inverse exner sent to conversion codw
                                                       ! pcols for output to history
      real(r8) :: eff_rad_coef = 1.0_r8/(4.0_r8/3.0_r8*SHR_CONST_RHOFW*SHR_CONST_PI)
      real(r8), dimension(pver) :: eff_rad_prof ! r^3 as calculated from grid mean MR & NC
     
      !----------------
      ! Pointers
      !----------------
      real(r8), pointer, dimension(:) :: ztodt_ptr
      real(r8), pointer, dimension(:,:) :: thlm      ! Mean temperature
      real(r8), pointer, dimension(:,:) :: ice_supersat_frac ! ice cloud fraction
      real(r8), pointer, dimension(:,:) :: rcm       ! CLUBB cld water mr
      real(r8), pointer, dimension(:,:) :: rtm       ! mean moisture mixing ratio
      real(r8), pointer, dimension(:,:) :: cld       ! CAM cloud fraction
      real(r8), pointer, dimension(:,:) :: alst      ! CLUBB liq cloud fraction
      real(r8), pointer, dimension(:,:) :: qrain     ! micro_mg rain from previous step
      real(r8), pointer, dimension(:,:) :: qsnow     
      real(r8), pointer, dimension(:,:) :: nrain     ! micro_mg rain num conc 
      real(r8), pointer, dimension(:,:) :: nsnow

      real(r8), pointer, dimension(:,:) :: tke_in    ! TKE
      real(r8), pointer, dimension(:,:) :: khzm_in   ! Eddy diffusivity coef

      if (.not. allocated(state_sc%lat)) then
         call endrun('subcol_gen error: state_sc must be allocated before calling subcol_gen')
      end if

      ! Determine num of columns and which chunk we're working on and what timestep
      ngrdcol = state%ngrdcol
      ncol = state%ncol
      lchnk = state%lchnk
      iter = get_nstep() ! #KTCtodo: The model iteration is passed into SILHS without taking
                         !           substepping into account. I may need to change this in 
                         !           the future. Also, why does SILHS need an iter, but CLUBB
                         !           does not?
                         ! #ERDBG:   The model iteration number is not used in SILHS unless
                         !           sequence_length > 1, but nobody runs with that option.
      !----------------
      ! Establish associations between pointers and physics buffer fields
      !----------------
      call pbuf_get_field(pbuf, thlm_idx, thlm)
      call pbuf_get_field(pbuf, ztodt_idx, ztodt_ptr)
      call pbuf_get_field(pbuf, ice_supersat_idx, ice_supersat_frac)
      call pbuf_get_field(pbuf, rcm_idx, rcm)
      call pbuf_get_field(pbuf, rtm_idx, rtm)
      call pbuf_get_field(pbuf, alst_idx, alst)
      call pbuf_get_field(pbuf, cld_idx, cld)
      call pbuf_get_field(pbuf, qrain_idx, qrain)
      call pbuf_get_field(pbuf, qsnow_idx, qsnow)
      call pbuf_get_field(pbuf, nrain_idx, nrain)
      call pbuf_get_field(pbuf, nsnow_idx, nsnow)
      call pbuf_get_field(pbuf, tke_idx, tke_in)
      call pbuf_get_field(pbuf, kvh_idx, khzm_in)

      ! Read the clubb parameters in order to extract c_K.
      call read_parameters_api( -99, "", clubb_params )

      ! Pull c_K from clubb parameters.
      c_K = clubb_params(ic_K)

      !----------------
      ! Copy state and populate numbers and values of sub-columns
      !----------------
      ztodt = ztodt_ptr(1)
      numsubcol_arr(:) = 0  ! Start over each chunk
      numsubcol_arr(:ngrdcol) = subcol_SILHS_numsubcol ! Only set for valid grid columns
      call subcol_set_subcols(state, tend, numsubcol_arr, state_sc, tend_sc)

      ! The number of vertical grid levels used in CLUBB is pverp, which is originally
      ! set in the call to setup_clubb_core_api from subroutine clubb_ini_cam.  This
      ! is stored in CLUBB in the object gr%nz.  This isn't changed in CLUBB.
      ! However, when SILHS is used, SILHS only uses pverp - top_lev + 1 vertical grid
      ! levels and also uses the gr%nz object.  The value of gr%nz needs to be reset
      ! for SILHS here and then set again for CLUBB in subroutine clubb_tend_cam.
      gr%nz = pverp - top_lev + 1

      !----------------
      ! Loop over all the active grid columns in the chunk
      !----------------
      do i = 1, ngrdcol
      
         ! JHDBG: Big suspicion about that code
         ! V. Larson: I don't know what happens to arrays allocated with size
         ! num_subcols if num_subcols varies with the grid column.
         num_subcols = numsubcol_arr(i)
         stncol = 0         ! Each grid column needs to know how many subcolumns have gone by
         do k = 1, i-1
            ! stncol = stncol + numsubcol_arr(i-1)
            ! Eric Raut replaced i-1 with k in line immediately above.
            stncol = stncol + numsubcol_arr(k)
         enddo

         ! Setup the CLUBB vertical grid object. This must be done for each
         ! column as the z-distance between hybrid pressure levels can 
         ! change easily.
         sfc_elevation = state%zi(i,pverp)
         ! Define the CLUBB momentum grid (in height, units of m)
         do k = 1, pverp-top_lev+1
            zi_g(k) = state%zi(i,pverp-k+1)-sfc_elevation
         enddo
         ! Define the CLUBB thermodynamic grid (in units of m)
         do k = 1, pver-top_lev+1
            zt_g(k+1) = state%zm(i,pver-k+1)-state%zi(i,pverp)
         enddo
         ! Thermodynamic ghost point is below surface
         zt_g(1) = -1._r8*zt_g(2)
         ! Calculate the distance between grid levels on the host model grid,
         ! using host model grid indices.
         do k = top_lev, pver
            dz_g(k) = state%zi(i,k)-state%zi(i,k+1)
         enddo
         ! allocate grid object
         call setup_grid_heights_api( l_implemented, grid_type, &
                                      zi_g(2), zi_g(1), zi_g(1:pverp-top_lev+1), &
                                      zt_g(1:pverp-top_lev+1) )

         ! Inverse delta_zm is required for the 3-level L-scale averaging
         do k = 1, pver-top_lev+1
            delta_zm(k+1) = state%zi(i,pverp-k)-state%zi(i,pverp-k+1)
            invs_dzm(k+1) = 1.0_r8/delta_zm(k+1)
         enddo
         ! Handle CLUBB sub-sfc ghost point as done in clubb grid_class.F90
         delta_zm(1) = delta_zm(2) 
         invs_dzm(1) = invs_dzm(2)

         ! Compute dry static density on CLUBB vertical grid
         do k = 1, pver-top_lev+1
            rho_ds_zt(k+1) = (1._r8/gravit)*state%pdel(i,pver-k+1)/dz_g(pver-k+1)
         enddo
         ! CLUBB ghost point under the surface
         rho_ds_zt(1) = rho_ds_zt(2)

         ! Set up hydromet array, flipped from CAM vert grid to CLUBB
         do k = 1, pver-top_lev+1
            if ( iirr > 0 ) then
              ! If ixrain and family are greater than zero, then MG2 is
              ! being used, and rain and snow are part of state. Otherwise,
              ! diagnostic rain and snow from MG1 are used in hydromet.
               if (ixrain > 0) then
                  hydromet(k+1,iirr) = state%q(i,pver-k+1,ixrain)
               else
                  hydromet(k+1,iirr) = qrain(i,pver-k+1)
               endif
            endif
            if ( iiNr > 0 ) then
               if (ixnumrain > 0) then
                  hydromet(k+1,iiNr) = state%q(i,pver-k+1,ixnumrain)
               else
                  hydromet(k+1,iiNr) = nrain(i,pver-k+1)
               endif
            endif
            if ( iirs > 0 ) then
               if (ixsnow > 0) then
                  hydromet(k+1,iirs) = state%q(i,pver-k+1,ixsnow)
               else
                  hydromet(k+1,iirs) = qsnow(i,pver-k+1)
               endif
            endif
            if ( iiNs > 0 ) then
               if (ixnumsnow > 0) then
                  hydromet(k+1,iiNs) = state%q(i,pver-k+1,ixnumsnow)
               else
                  hydromet(k+1,iiNs) = nsnow(i,pver-k+1)
               endif
            endif
            if ( iiri > 0 ) then
               hydromet(k+1,iiri) = state%q(i,pver-k+1,ixcldice)
            endif
            if ( iiNi > 0 ) then
               hydromet(k+1,iiNi) = state%q(i,pver-k+1,ixnumice)
            endif
     
            Ncm(k+1) = state%q(i,pver-k+1,ixnumliq)

         enddo

         do k = 1, hydromet_dim ! ghost point below the surface
            hydromet(1,k) = hydromet(2,k)                  
         enddo

         Ncm(1) = Ncm(2)

         do k = top_lev, pver
            ! Calculate effective radius cubed, CAM-grid oriented for use in subcolumns
            eff_rad_prof(k) = eff_rad_coef*state%q(i,k,ixcldliq)/state%q(i,k,ixnumliq)
            ! Test a fixed effective radius
            ! eff_rad_prof(k) = 5.12e-16_r8 ! 8 microns
         enddo

         ! Allocate arrays for set_up_pdf_params_incl_hydromet
         allocate( corr_array_1(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( corr_array_2(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( mu_x_1(pdf_dim, pverp-top_lev+1) )
         allocate( mu_x_2(pdf_dim, pverp-top_lev+1) )
         allocate( sigma_x_1(pdf_dim, pverp-top_lev+1) )
         allocate( sigma_x_2(pdf_dim, pverp-top_lev+1) )
         allocate( corr_cholesky_mtx_1(pdf_dim, pdf_dim, pverp-top_lev+1) )
         allocate( corr_cholesky_mtx_2(pdf_dim, pdf_dim, pverp-top_lev+1) )
         ! Allocate arrays for SILHS output
         allocate( lh_sample_point_weights(pverp-top_lev+1,num_subcols) )
         allocate( X_mixt_comp_all_levs(pverp-top_lev+1,num_subcols) )
         allocate( X_nl_all_levs(pverp-top_lev+1,num_subcols,pdf_dim) )
         allocate( lh_rt_clipped(pverp-top_lev+1,num_subcols) )
         allocate( lh_thl_clipped(pverp-top_lev+1,num_subcols) )
         allocate( lh_rc_clipped(pverp-top_lev+1,num_subcols) )
         allocate( lh_rv_clipped(pverp-top_lev+1,num_subcols) )
         allocate( lh_Nc_clipped(pverp-top_lev+1,num_subcols) )
         ! Allocate arrays for output to either history files or for updating state_sc
         allocate( rc_all_points(pverp-top_lev+1, num_subcols) )
         allocate( rain_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nrain_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( snow_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nsnow_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( w_all_points(pverp-top_lev+1, num_subcols) )
         ! allocate( RVM_lh_out(num_subcols, pverp) )  ! This one used only to update state
         allocate( ice_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nice_all_pts(pverp-top_lev+1, num_subcols) )
         allocate( nclw_all_pts(pverp-top_lev+1, num_subcols) )
         
         ! Convert from CAM vertical grid to CLUBB
         do k = 1, pverp-top_lev+1 
            rcm_in(k)  = rcm(i,pverp-k+1)
            ice_supersat_frac_in(k) = ice_supersat_frac(i,pverp-k+1)
         enddo
         do k = 1, pver-top_lev+1
            cld_frac_in(k+1) = alst(i,pver-k+1)
         enddo
         cld_frac_in(1) = cld_frac_in(2) ! Ghost pt below surface
         ! Calculate a clubb-specific exner function
         ! (This is grid mean, as pressure levels do not change in 
         !  the subcolumn state)
         invs_exner(i,:) = ((state%pmid(i,:)/p0_clubb)**(rair/cpair))

         ! Call setup_pdf_parameters to get the CLUBB PDF ready for SILHS
         ! Compute Num concentration of cloud nuclei
         Nc_in_cloud = Ncm / max( cld_frac_in, cloud_frac_min )

         ! The variable wphydrometp is only used when l_calc_w_corr is enabled.
         ! The l_calc_w_corr flag is turned off by default, so wphydrometp will
         ! simply be set to 0 to simplify matters.
         wphydrometp = 0.0_r8
     
         ! make the call
         call setup_pdf_parameters_api( pverp-top_lev+1, pdf_dim, ztodt, &                 ! In
                                        Nc_in_cloud, rcm_in, cld_frac_in, &                ! In
                                        ice_supersat_frac_in, hydromet, wphydrometp, &     ! In
                                        corr_array_n_cloud, corr_array_n_below, &          ! In
                                        pdf_params_chnk(i,lchnk), l_stats_samp, &          ! In
                                        clubb_config_flags%l_use_precip_frac, &            ! In
                                        clubb_config_flags%l_predict_upwp_vpwp, &          ! In
                                        clubb_config_flags%l_diagnose_correlations, &      ! In
                                        clubb_config_flags%l_calc_w_corr, &                ! In
                                        clubb_config_flags%l_const_Nc_in_cloud, &          ! In
                                        clubb_config_flags%l_fix_w_chi_eta_correlations, & ! In
                                        hydrometp2, &                                      ! Out
                                        mu_x_1, mu_x_2, &                                  ! Out
                                        sigma_x_1, sigma_x_2, &                            ! Out
                                        corr_array_1, corr_array_2, &                      ! Out
                                        corr_cholesky_mtx_1, corr_cholesky_mtx_2, &        ! Out
                                        hydromet_pdf_params )                              ! Out

         ! Calculate radiation only once in a while
         ! l_rad_itime = (mod( itime, floor(dt_rad/dt_main) ) == 0 .or. itime == 1)  

         ! Calculate sample weights separately at all grid levels when
         ! radiation is not called  
         ! l_calc_weights_all_levs_itime = l_calc_weights_all_levs .and. .not.
         ! l_rad_itime  
         l_calc_weights_all_levs_itime = .false. ! subcol_utils cannot compute weighted avgs
                                                 !   when the weights vary with height.   
                                                 !   Don't set to true until this is fixed!!

         ! In order for Lscale to be used properly, it needs to be passed out of
         ! advance_clubb_core, saved to the pbuf, and then pulled out of the
         ! pbuf for use here.  The profile of Lscale is passed into subroutine
         ! generate_silhs_sample_api for use in calculating the vertical
         ! correlation coefficient.  Rather than output Lscale directly, its
         ! value can be calculated from other fields that are already output to
         ! pbuf.  The equation relating Lscale to eddy diffusivity is:
         !
         ! Kh = c_K * Lscale * sqrt( TKE ).
         !
         ! Both Kh and TKE are written to the pbuf, and c_K is easily extracted
         ! from CLUBB's tunable parameters.  The equation for Lscale is:
         !
         ! Lscale = Kh / ( c_K * sqrt( TKE ) ).
         !
         ! Since Kh and TKE are output on momentum (interface) grid levels, the
         ! resulting calculation of Lscale is also found on momentum levels.  It
         ! needs to be interpolated back to thermodynamic (midpoint) grid levels
         ! for further use.
         do k = 1, pverp-top_lev+1
            khzm(k) = khzm_in(i,pverp-k+1)
            tke(k)  = tke_in(i,pverp-k+1)
         enddo
         Lscale_zm = khzm / ( c_K * sqrt( max( tke, em_min ) ) )

         ! Interpolate Lscale_zm back to thermodynamic grid levels.
         Lscale = max( zm2zt_api( Lscale_zm ), 0.01_r8 )

         ! Set the seed to the random number generator based on a quantity that
         ! will be reproducible for restarts.
         lh_seed = int( 1.0e4_r8 * rtm(i,pver), kind = genrand_intg )
         call genrand_init_api( put=lh_seed )

         ! Let's generate some subcolumns!!!!!
         call generate_silhs_sample_api &
              ( iter, pdf_dim, num_subcols, sequence_length, pverp-top_lev+1, & ! In
                l_calc_weights_all_levs_itime, &                   ! In 
                pdf_params_chnk(i,lchnk), delta_zm, rcm_in, Lscale, & ! In
                rho_ds_zt, mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, & ! In 
                corr_cholesky_mtx_1, corr_cholesky_mtx_2, &        ! In
                hydromet_pdf_params, silhs_config_flags, &         ! In
                clubb_config_flags%l_uv_nudge, &                   ! In
                clubb_config_flags%l_tke_aniso, &                  ! In
                clubb_config_flags%l_standard_term_ta, &           ! In
                clubb_config_flags%l_single_C2_Skw, &              ! In
                X_nl_all_levs, X_mixt_comp_all_levs, &             ! Out
                lh_sample_point_weights)                           ! Out

         ! Extract clipped variables from subcolumns
         call clip_transform_silhs_output_api( pverp-top_lev+1, num_subcols, &   ! In
                                               pdf_dim, hydromet_dim, & ! In
                                               X_mixt_comp_all_levs, & ! In
                                               X_nl_all_levs, &        ! In
                                               pdf_params_chnk(i,lchnk), & ! In
                                               l_use_Ncn_to_Nc, & ! In
                                               lh_rt_clipped, lh_thl_clipped, & ! Out
                                               lh_rc_clipped, lh_rv_clipped, & ! Out
                                               lh_Nc_clipped ) ! Out

         ! Test subcolumns by comparing to an estimate of kessler autoconversion
         call est_kessler_microphys_api &
              ( pverp-top_lev+1, num_subcols, pdf_dim, X_nl_all_levs, &
                pdf_params_chnk(i,lchnk), &
                rcm_in, cld_frac_in, X_mixt_comp_all_levs, lh_sample_point_weights, &
                silhs_config_flags%l_lh_importance_sampling, &
                lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, lh_rcm_avg)

         ! Calc column liquid water for output (rcm)
         rc_all_points = lh_rc_clipped(:,:)

         if ( iiPDF_rr > 0 ) then
             ! Calc subcolumn precipitating liq water for output (rrm)
             rain_all_pts = real( X_nl_all_levs(:,:,iiPDF_rr), kind=r8 )
         end if

         if ( iiPDF_Nr > 0 ) then
             ! Calc subcolumn number rain conc for output (nrainm)
             nrain_all_pts = real( X_nl_all_levs(:,:,iiPDF_Nr), kind=r8 )
         end if

         if ( iiPDF_rs > 0 ) then
             ! Calc subcolumn precipitating snow      for output (rsm)
             snow_all_pts = real( X_nl_all_levs(:,:,iiPDF_rs), kind=r8 )
         end if

         if ( iiPDF_Ns > 0 ) then
             ! Calc subcolumn precipitating snow conc for output (Nsm)
             nsnow_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ns), kind=r8 )
         end if

         if ( iiPDF_ri > 0 ) then
             ! Calc subcolumn cloud ice mixing ratio
             ice_all_pts = real( X_nl_all_levs(:,:,iiPDF_ri), kind=r8)
         end if

         if ( iiPDF_Ni > 0 ) then
             ! Calc subcolumn cloud ice number
             nice_all_pts = real( X_nl_all_levs(:,:,iiPDF_Ni), kind=r8)
         end if

         ! Calc subcolumn vert velocity for output (wm)
         w_all_points = real( X_nl_all_levs(:,:,iiPDF_w), kind=r8 )
         ! Calc cloud liq water number conc 
         nclw_all_pts = lh_Nc_clipped(:,:)
         ! Calc mean liquid water potential temp for clear air
         !call THL_profile(pver, state%t(i,:), invs_exner(i,:), No_cloud, Temp_prof)

         ! Calc effective cloud fraction for testing
         eff_cldfrac(:,:) = 0.0_r8
         do k = top_lev, pver
            do j=1, num_subcols

               if ( ( rc_all_points(pverp-k+1,j) .gt. qsmall ) &
                      .or. ( ice_all_pts(pverp-k+1,j) .gt. qsmall ) ) then
                  eff_cldfrac(i,k) = eff_cldfrac(i,k)+lh_sample_point_weights(pverp-k+1,j)
               endif
            enddo 

            eff_cldfrac(i,k) = eff_cldfrac(i,k)/real(num_subcols, kind=r8)
         enddo

         ! Pack precip_frac for output
         do k = 2, pverp-top_lev+1
           precip_frac_out(i,pver-k+2) = hydromet_pdf_params(k)%precip_frac
         enddo

         ! Pack up weights for output
         do j = 1, num_subcols      
            if (subcol_SILHS_weight) then 
               weights(stncol+j) = lh_sample_point_weights(2,j) ! Using grid level 2 always won't work 
                                                                !   if weights vary with height.
            else
               weights(stncol+j) = 1._r8
            endif
         enddo

         ! Convert from CLUBB vertical grid to CAM grid for history output and
         ! Updating state variables
         do k = top_lev, pverp
            do j = 1, num_subcols
               RT_lh_out(    stncol+j,k ) = lh_rt_clipped(pverp-k+1,j)
               RCM_lh_out(   stncol+j,k ) = rc_all_points(pverp-k+1,j)
               NCLW_lh_out(  stncol+j,k ) = nclw_all_pts(pverp-k+1,j)
               ICE_lh_out(   stncol+j,k ) = ice_all_pts(pverp-k+1,j)
               NICE_lh_out(  stncol+j,k ) = nice_all_pts(pverp-k+1,j)
!               RVM_lh_out(j,k) = RT_lh_out(stncol+j,k)-RCM_lh_out(stncol+j,k)-ICE_lh_out(stncol+j,k)
               RVM_lh_out(   stncol+j,k ) = lh_rv_clipped(pverp-k+1,j)
               THL_lh_out(   stncol+j,k ) = lh_thl_clipped(pverp-k+1,j)
               RAIN_lh_out(  stncol+j,k ) = rain_all_pts(pverp-k+1,j)
               NRAIN_lh_out( stncol+j,k ) = nrain_all_pts(pverp-k+1,j)
               SNOW_lh_out(  stncol+j,k ) = snow_all_pts(pverp-k+1,j)
               NSNOW_lh_out( stncol+j,k ) = nsnow_all_pts(pverp-k+1,j)
               WM_lh_out(    stncol+j,k ) = w_all_points(pverp-k+1,j)
               OMEGA_lh_out( stncol+j,k ) = -1._r8*WM_lh_out(stncol+j,k)*rho_ds_zt(pverp-k+1)*gravit
               AKm_out(i,k) = AKm(pverp-k+1)
               lh_AKm_out(i,k) = lh_AKm(pverp-k+1)
            enddo
         enddo

         ! Constrain the sample distribution of cloud water and ice to the same mean
         ! as the grid to prevent negative condensate errors
         if(subcol_SILHS_constrainmn) then
            call subcol_constrainmn( num_subcols, ICE_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldice), meansc_ice(i,:), stdsc_ice(i,:) )
            if ( ixrain > 0 ) &
            call subcol_constrainmn( num_subcols, RAIN_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixrain) )
            if ( ixsnow > 0 ) &
            call subcol_constrainmn( num_subcols, SNOW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixsnow) )
            call subcol_constrainmn( num_subcols, RCM_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixcldliq), meansc_liq(i,:), stdsc_liq(i,:) )
            call subcol_constrainmn( num_subcols, RVM_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixq), meansc_vap(i,:), stdsc_vap(i,:) )
            call subcol_constrainmn( num_subcols, NICE_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumice) )
            if ( ixnumrain > 0 ) &
            call subcol_constrainmn( num_subcols, NRAIN_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumrain) )
            if ( ixnumsnow > 0 ) &
            call subcol_constrainmn( num_subcols, NSNOW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumsnow) )
            call subcol_constrainmn( num_subcols, NCLW_lh_out(stncol+1:stncol+num_subcols,:), &
                                     weights(stncol+1:stncol+num_subcols), &
                                     state%q(i,:,ixnumliq) )
            do k = top_lev, pver
               ! Look for exceptionally large values of condensate
               if(ANY(ICE_lh_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)) then
                  ! Clip the large values
                  where(ICE_lh_out(stncol+1:stncol+num_subcols,k) .gt. 0.01_r8)
                     ICE_lh_out(stncol+1:stncol+num_subcols,k) = 0.01_r8
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = 1.5e+7_r8
                  end where
                  ! Recalculate the weighted subcolumn mean
                  tmp_mean = meansc( ICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  ! Add the difference to each subcolumn
                  ICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                     ICE_lh_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  ! Recalculate the weight subcolumn mean for ice num conc
                  tmp_mean = meansc( NICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  ! Calculate the difference between the weighted mean and grid mean
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
                  ! Add the difference to each subcolumn
                  if(diff_mean.gt.0.0_r8) then
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                         NICE_lh_out(stncol+1:stncol+num_subcols,k)+diff_mean
                  else ! just use the grid mean in each subcolumn
                     NICE_lh_out(stncol+1:stncol+num_subcols,k) = &
                         state%q(i,k,ixnumice)
                  end if
                  ! Test adjusted means for debugging
                  tmp_mean = meansc( ICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixcldice)-tmp_mean
                  tmp_mean = meansc( NICE_lh_out( stncol+1:stncol+num_subcols, k ), & 
                                         weights(stncol+1:stncol+num_subcols), &
                                         real(num_subcols,r8) )
                  diff_mean = state%q(i,k,ixnumice)-tmp_mean
               endif
            enddo ! k = top_lev, pver
         endif ! subcol_silhs_constrainm

         ! Code to update the state variables for interactive runs
         ! Set state variables
         do j = 1, numsubcol_arr(i)

            call Abs_Temp_profile( pver-top_lev+1, THL_lh_out(stncol+j,top_lev:pver), &
                                   invs_exner(i,top_lev:pver), RCM_lh_out(stncol+j,top_lev:pver), &
                                   Temp_prof(top_lev:pver) )
            state_sc%t(stncol+j,top_lev:pver) = Temp_prof(top_lev:pver)
            call StaticEng_profile( pver-top_lev+1, Temp_prof(top_lev:pver), &
                                    state%zm(i,top_lev:pver), state%phis(i), &
                                    SE_prof(top_lev:pver) )
            state_sc%s(stncol+j,top_lev:pver) = SE_prof(top_lev:pver)

            ! Vertical Velocity is not part of the energy conservation checks, but
            ! we need to be careful here, because the SILHS output VV is noisy.
            state_sc%omega(stncol+j,top_lev:pver) = OMEGA_lh_out(stncol+j,top_lev:pver)
            state_sc%q(stncol+j,top_lev:pver,ixq) = RVM_lh_out(stncol+j,top_lev:pver) 

            if( rx_Nc ) then
               call endrun('subcol_gen_SILHS: rx_Nc not enabled')
            endif


            if (subcol_SILHS_meanice) then
               call endrun('subcol_gen_SILHS: subcol_SILHS_meanice = T not currently available')
                state_sc%q(stncol+j,top_lev:pver,ixcldice) = state%q(i,top_lev:pver,ixcldice)
                state_sc%q(stncol+j,top_lev:pver,ixnumice) = state%q(i,top_lev:pver,ixnumice)
                state_sc%q(stncol+j,top_lev:pver,ixcldliq) = RCM_lh_out(stncol+j,top_lev:pver)
                state_sc%q(stncol+j,top_lev:pver,ixnumliq) = NCLW_lh_out(stncol+j,top_lev:pver)
            else
               if (subcol_SILHS_q_to_micro) then ! Send SILHS predicted constituents to microp
                   state_sc%q(stncol+j,top_lev:pver,ixcldliq) = RCM_lh_out(stncol+j,top_lev:pver)
                   state_sc%q(stncol+j,top_lev:pver,ixcldice) = ICE_lh_out(stncol+j,top_lev:pver)
                   if (ixrain > 0) &
                      state_sc%q(stncol+j,top_lev:pver,ixrain) = RAIN_lh_out(stncol+j,top_lev:pver)
                   if (ixsnow > 0) &
                      state_sc%q(stncol+j,top_lev:pver,ixsnow) = SNOW_lh_out(stncol+j,top_lev:pver)
               else            
                  state_sc%q(stncol+j,top_lev:pver,ixcldliq) = state%q(i,top_lev:pver,ixcldliq)
                  state_sc%q(stncol+j,top_lev:pver,ixcldice) = state%q(i,top_lev:pver,ixcldice)
                  if (ixrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixrain) = state%q(i,top_lev:pver,ixrain)
                  if (ixsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixsnow) = state%q(i,top_lev:pver,ixsnow)
               endif
               if (subcol_SILHS_n_to_micro) then ! Send SILHS predicted number conc to microp
                  state_sc%q(stncol+j,top_lev:pver,ixnumice) = NICE_lh_out(stncol+j,top_lev:pver)
                  state_sc%q(stncol+j,top_lev:pver,ixnumliq) = NCLW_lh_out(stncol+j,top_lev:pver)
                  if (ixnumrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumrain) = NRAIN_lh_out(stncol+j,top_lev:pver)
                  if (ixnumsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = NSNOW_lh_out(stncol+j,top_lev:pver)
               else            
                  state_sc%q(stncol+j,top_lev:pver,ixnumliq) = state%q(i,top_lev:pver,ixnumliq)
                  state_sc%q(stncol+j,top_lev:pver,ixnumice) = state%q(i,top_lev:pver,ixnumice)
                  if (ixnumrain > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumrain) = state%q(i,top_lev:pver,ixnumrain)
                  if (ixnumsnow > 0) &
                     state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = state%q(i,top_lev:pver,ixnumsnow)
               endif
            endif ! meanice

            ! Change liq and ice (and rain and snow) num conc zeros to min values (1e-12)
            where (state_sc%q(stncol+j,top_lev:pver,ixnumliq) .lt. min_num_conc) 
               state_sc%q(stncol+j,top_lev:pver,ixnumliq) = min_num_conc
            end where
            where (state_sc%q(stncol+j,top_lev:pver,ixnumice) .lt. min_num_conc)
               state_sc%q(stncol+j,top_lev:pver,ixnumice) = min_num_conc
            end where
            if (ixnumrain > 0) then
               where(state_sc%q(stncol+j,top_lev:pver,ixnumrain) .lt. min_num_conc)
                  state_sc%q(stncol+j,top_lev:pver,ixnumrain) = min_num_conc
               end where
            endif
            if (ixnumsnow > 0) then
               where(state_sc%q(stncol+j,top_lev:pver,ixnumsnow) .lt. min_num_conc)
                  state_sc%q(stncol+j,top_lev:pver,ixnumsnow) = min_num_conc
               end where
            endif
               
         enddo

         ! Only use weights if namelist variable turned on
         if (subcol_SILHS_weight) call subcol_set_weight(state_sc%lchnk, weights)

     
         ! Deallocate the dynamic arrays used
         deallocate( lh_sample_point_weights, X_mixt_comp_all_levs, &
                     X_nl_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
                     lh_rv_clipped, lh_Nc_clipped, &
                     corr_array_1, corr_array_2, mu_x_1, mu_x_2, sigma_x_1, &
                     sigma_x_2, corr_cholesky_mtx_1, corr_cholesky_mtx_2 )
         ! deallocate( RVM_lh_out ) 
         deallocate( rc_all_points, rain_all_pts, nrain_all_pts, snow_all_pts, nsnow_all_pts, ice_all_pts, &
                     nice_all_pts, nclw_all_pts, w_all_points )
      enddo ! ngrdcol

      call outfld( 'SILHS_THLM_SCOL', THL_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RT_SCOL', RT_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_OMEGA_SCOL', OMEGA_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WM_SCOL', WM_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RCM_SCOL', RCM_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RICLD_SCOL', ICE_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NICLD_SCOL', NICE_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NCLD_SCOL', NCLW_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_RRAIN_SCOL', RAIN_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_NRAIN_SCOL', NRAIN_lh_out, pcols*psubcols, lchnk )
      call outfld( 'SILHS_WEIGHT_SCOL', weights, pcols*psubcols, lchnk )
      call outfld( 'NR_IN_LH', nrain, pcols, lchnk )
      call outfld( 'RTM_CLUBB', rtm, pcols, lchnk )
      call outfld( 'THLM_CLUBB', thlm, pcols, lchnk )
      call outfld( 'SILHS_QC_IN', state%q(:,:,ixcldliq), pcols, lchnk )
      call outfld( 'SILHS_QI_IN', state%q(:,:,ixcldice), pcols, lchnk )
      call outfld( 'SILHS_NC_IN', state%q(:,:,ixnumliq), pcols, lchnk )
      call outfld( 'SILHS_NI_IN', state%q(:,:,ixnumice), pcols, lchnk )
      call outfld( 'AKM_CLUBB', AKm_out, pcols, lchnk )
      call outfld( 'AKM_LH_CLUBB', lh_AKm_out, pcols, lchnk )
      call outfld( 'INVS_EXNER', invs_exner, pcols, lchnk )
      call outfld( 'SILHS_ZTODT', ztodt_ptr, pcols, lchnk )
      if ( subcol_SILHS_constrainmn ) then
         call outfld( 'SILHS_MSC_CLDICE', meansc_ice, pcols, lchnk )
         call outfld( 'SILHS_STDSC_CLDICE', stdsc_ice, pcols, lchnk )
         if ( ixsnow > 0 ) then
            call outfld( 'SILHS_MSC_CLDLIQ', meansc_liq, pcols, lchnk )
            call outfld( 'SILHS_STDSC_CLDLIQ', stdsc_liq, pcols, lchnk )
            call outfld( 'SILHS_MSC_Q', meansc_vap, pcols, lchnk )
            call outfld( 'SILHS_STDSC_Q', stdsc_vap, pcols, lchnk )
         endif ! ixsnow > 0
      endif ! subcol_SILHS_constrainmn
      call outfld( 'SILHS_EFF_CLDFRAC', eff_cldfrac, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_PRECIP_FRAC', precip_frac_out, pcols, lchnk )
      call outfld( 'SILHS_CLUBB_ICE_SS_FRAC', ice_supersat_frac, pcols, lchnk )

#endif
#endif
   end subroutine subcol_gen_SILHS

   subroutine subcol_ptend_avg_SILHS(ptend_sc, ngrdcol, lchnk, ptend)
      use physics_buffer,   only: physics_buffer_desc
      use subcol_utils,     only: subcol_ptend_get_firstsubcol, subcol_ptend_avg_shr, &
                                  subcol_get_weight, subcol_get_filter, &
                                  is_filter_set, is_weight_set

      !-----------------------------------
      ! Average the subcolumns dimension (pcols*psubcols) to the grid dimension (pcols)
      !-----------------------------------

      type(physics_ptend), intent(in)             :: ptend_sc        ! intent in
      integer,  intent(in)                        :: ngrdcol       ! # grid cols
      integer,  intent(in)                        :: lchnk         ! chunk index
      type(physics_ptend), intent(inout)          :: ptend
      ! Because we can't get a state passed in here, we might have to use values from the 
      ! subcolumn generation. This would make any conservation checks invalid if this
      ! function is called after another parameterization... hmm.

       call subcol_ptend_avg_shr(ptend_sc, ngrdcol, lchnk, ptend, is_filter_set(), is_weight_set())

   end subroutine subcol_ptend_avg_SILHS

   subroutine subcol_SILHS_var_covar_driver &
              ( ztodt, state_sc, ptend_sc, &
                pbuf )

     ! This subroutine calculates microphysical effects on five variances and
     ! covariances: rtp2, thlp2, wprtp, wpthlp, and rtpthlp.
     !
     ! This code is experimental!!

     use physics_buffer,          only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
#ifdef CLUBB_SGS
#ifdef SILHS
     use ref_pres,                only: top_lev => trop_cloud_top_lev
     use subcol_utils,            only: subcol_get_weight
     use subcol_pack_mod,         only: subcol_unpack, subcol_get_nsubcol
     use clubb_api_module,        only: T_in_K2thlm_api
     use silhs_api_module,        only: lh_microphys_var_covar_driver_api
#endif
#endif

     implicit none

     ! Parameters
     !  This fill value is set to catch errors; it should not be read.
     real(r8), parameter                   :: fillvalue = -999._r8

     ! Input Variables
     real(r8), intent(in)                  :: ztodt        ! model time increment
     type(physics_state), intent(in)       :: state_sc     ! state for sub-columns
     type(physics_ptend), intent(in)       :: ptend_sc     ! ptend for sub-columns

     ! Pointers
     type(physics_buffer_desc), pointer    :: pbuf(:)

#ifdef CLUBB_SGS
#ifdef SILHS
     ! Local Variables
     integer :: lchnk, ngrdcol, igrdcol, isubcol, ns, k
     integer, dimension(pcols) :: nsubcol
     real(r8), dimension(pcols*psubcols)       :: weights_packed
     real(r8), dimension(pcols,psubcols)       :: weights
     real(r8), dimension(pcols,psubcols,pverp) :: rc_all, rv_all, rt_all, w_all, thl_all
     real(r8), dimension(pcols,psubcols,pver ) :: s_all, t_all, zm_all, omega_all, pmid_all
     real(r8), dimension(pcols,psubcols)       :: phis_all
     real(r8), dimension(pcols,psubcols,pver ) :: stend, ttend
     real(r8), dimension(pcols,psubcols,pverp) :: thltend, qctend, qvtend

     real(r8), dimension(pcols,psubcols,pver)  :: dz_g, pdel_all, rho
     real(r8), dimension(pcols,psubcols,pverp) :: zi_all
 
     real(r8), dimension(pcols,psubcols,pver ) :: exner

     ! Inputs to lh_microphys_var_covar_driver
     real(r8), dimension(pcols,pverp,psubcols) :: rt_all_clubb, thl_all_clubb, w_all_clubb, &
                                                  qctend_clubb, qvtend_clubb, thltend_clubb
     real(r8), dimension(pcols,pverp-top_lev+1,psubcols) :: height_depndt_weights

     ! Outputs from lh_microphys_var_covar_driver
     real(r8), dimension(:,:), pointer :: rtp2_mc_zt, thlp2_mc_zt, wprtp_mc_zt, &
                                          wpthlp_mc_zt, rtpthlp_mc_zt

     ! pbuf indices
     integer :: &
       rtp2_mc_zt_idx,    &
       thlp2_mc_zt_idx,   &
       wprtp_mc_zt_idx,   &
       wpthlp_mc_zt_idx,  &
       rtpthlp_mc_zt_idx

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. subcol_SILHS_var_covar_src ) return

     lchnk = state_sc%lchnk
     ngrdcol  = state_sc%ngrdcol

     ! Obtain indices
     rtp2_mc_zt_idx = pbuf_get_index('rtp2_mc_zt')
     thlp2_mc_zt_idx = pbuf_get_index('thlp2_mc_zt')
     wprtp_mc_zt_idx = pbuf_get_index('wprtp_mc_zt')
     wpthlp_mc_zt_idx = pbuf_get_index('wpthlp_mc_zt')
     rtpthlp_mc_zt_idx = pbuf_get_index('rtpthlp_mc_zt')

     ! Obtain pbuf fields for output
     call pbuf_get_field(pbuf, rtp2_mc_zt_idx, rtp2_mc_zt)
     call pbuf_get_field(pbuf, thlp2_mc_zt_idx, thlp2_mc_zt)
     call pbuf_get_field(pbuf, wprtp_mc_zt_idx, wprtp_mc_zt)
     call pbuf_get_field(pbuf, wpthlp_mc_zt_idx, wpthlp_mc_zt)
     call pbuf_get_field(pbuf, rtpthlp_mc_zt_idx, rtpthlp_mc_zt)

     ! Unpack needed tendencies from subcolumn ptends
     call subcol_unpack(lchnk, ptend_sc%s(:,:), stend, fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixcldliq), qctend(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, ptend_sc%q(:,:,ixq), qvtend(:,:,1:pver), fillvalue)

     ! Unpack sample point values from subcolumn states
     call subcol_unpack(lchnk, state_sc%q(:,:,ixcldliq), rc_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%q(:,:,ixq),      rv_all(:,:,1:pver), fillvalue)
     call subcol_unpack(lchnk, state_sc%omega(:,:),      omega_all (:,:,:),  fillvalue)
     call subcol_unpack(lchnk, state_sc%s(:,:),          s_all,              fillvalue)
     call subcol_unpack(lchnk, state_sc%zm,              zm_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%phis,            phis_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%zi,              zi_all,             fillvalue)
     call subcol_unpack(lchnk, state_sc%pdel,            pdel_all,           fillvalue)
     call subcol_unpack(lchnk, state_sc%pmid,            pmid_all,           fillvalue)

     ! Initialize fields to fillvalue.
     rt_all  = fillvalue
     thl_all = fillvalue
     w_all   = fillvalue
     qctend  = fillvalue
     qvtend  = fillvalue
     thltend = fillvalue

     ! How many subcolumns in each column?
     call subcol_get_nsubcol(lchnk, nsubcol)

     do igrdcol = 1, ngrdcol
        do isubcol = 1, nsubcol(igrdcol)

           rt_all(igrdcol,isubcol,top_lev:pver) = rc_all(igrdcol,isubcol,top_lev:pver) &
                                                  + rv_all(igrdcol,isubcol,top_lev:pver)

           ! Compute dry static density on CLUBB vertical grid
           do k = top_lev, pver
              dz_g(igrdcol,isubcol,k) = zi_all(igrdcol,isubcol,k) - zi_all(igrdcol,isubcol,k+1) ! thickness
              rho(igrdcol,isubcol,k) = (1._r8/gravit)*pdel_all(igrdcol,isubcol,k)/dz_g(igrdcol,isubcol,k)
           enddo

           ! Compute w from omega
           w_all(igrdcol,isubcol,top_lev:pver) = -omega_all(igrdcol,isubcol,top_lev:pver) &
                                                  / ( rho(igrdcol,isubcol,top_lev:pver) * gravit )

           ! Convert stend and s_all to ttend and t_all
           !  Note 1: With subcolumns, cpair is truly a constant (I think).
           !  Note 2: For tendencies, the extra terns zm and phis should
           !          not be included in the calculation.
           ttend(igrdcol,isubcol,top_lev:pver) = stend(igrdcol,isubcol,top_lev:pver) / cpair

           do k = top_lev, pver
              t_all(igrdcol,isubcol,k) = ( s_all(igrdcol,isubcol,k) &
                                           - gravit * zm_all(igrdcol,isubcol,k) &
                                           - phis_all(igrdcol,isubcol) ) / cpair
           enddo ! k = 1, pver

           ! This formula is taken from earlier in this file.
           exner(igrdcol,isubcol,top_lev:pver) &
           = ( pmid_all(igrdcol,isubcol,top_lev:pver) / p0_clubb )**(rair/cpair)

           ! Note: all tendencies or all means should be used in the call to
           !       T_in_K2thlm_api (with the exception of exner)
           do k = top_lev, pver
              thltend(igrdcol,isubcol,k) &
              = T_in_K2thlm_api( ttend(igrdcol,isubcol,k), exner(igrdcol,isubcol,k), &
                                 qctend(igrdcol,isubcol,k) )
              thl_all(igrdcol,isubcol,k) &
              = T_in_K2thlm_api( t_all(igrdcol,isubcol,k), exner(igrdcol,isubcol,k), &
                                 rc_all(igrdcol,isubcol,k) )
           enddo ! k = 1, pver

           ! Add ghost points
           rt_all (igrdcol,isubcol,pverp) = rt_all (igrdcol,isubcol,pver)
           thl_all(igrdcol,isubcol,pverp) = thl_all(igrdcol,isubcol,pver)
           w_all  (igrdcol,isubcol,pverp) = w_all  (igrdcol,isubcol,pver)
           qctend (igrdcol,isubcol,pverp) = qctend (igrdcol,isubcol,pver)
           qvtend (igrdcol,isubcol,pverp) = qvtend (igrdcol,isubcol,pver)
           thltend(igrdcol,isubcol,pverp) = thltend(igrdcol,isubcol,pver)

           ! Flip inputs to CLUBB's grid. Note the dimension ordering change.
           rt_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( rt_all(igrdcol,isubcol,1:pverp) )
           thl_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( thl_all(igrdcol,isubcol,1:pverp) )
           w_all_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( w_all(igrdcol,isubcol,1:pverp) )
           qctend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( qctend(igrdcol,isubcol,1:pverp) )
           qvtend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( qvtend(igrdcol,isubcol,1:pverp) )
           thltend_clubb(igrdcol,1:pverp,isubcol) = clubb_flip_grid( thltend(igrdcol,isubcol,1:pverp) )

        enddo ! isubcol = 1, nsubcol(igrdcol)
     enddo ! igrdcol = 1, ngrdcol

     ! Obtain weights
     call subcol_get_weight(lchnk, weights_packed)
     call subcol_unpack(lchnk, weights_packed, weights, fillvalue)

     ! Call lh_microphys_var_covar_driver for each column
     do igrdcol=1, ngrdcol
       ns = nsubcol(igrdcol)

       ! This code assumes that the weights are height independent.
       ! It will have to change once the weights vary with altitude!
       ! I'm not sure whether the grid will need to be flipped.
       do k = 1, pverp-top_lev+1
          height_depndt_weights(igrdcol,k,1:ns) = weights(igrdcol,1:ns)
       end do

       ! Make the call!!!!!
       call lh_microphys_var_covar_driver_api &
            ( pverp-top_lev+1, ns, ztodt, height_depndt_weights(igrdcol,1:pverp-top_lev+1,1:ns), &
              pdf_params_chnk(igrdcol,lchnk), &
              rt_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), thl_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              w_all_clubb(igrdcol,1:pverp-top_lev+1,1:ns), qctend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              qvtend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), thltend_clubb(igrdcol,1:pverp-top_lev+1,1:ns), &
              silhs_config_flags%l_lh_instant_var_covar_src, &
              rtp2_mc_zt(igrdcol,1:pverp-top_lev+1), thlp2_mc_zt(igrdcol,1:pverp-top_lev+1), &
              wprtp_mc_zt(igrdcol,1:pverp-top_lev+1), wpthlp_mc_zt(igrdcol,1:pverp-top_lev+1), &
              rtpthlp_mc_zt(igrdcol,1:pverp-top_lev+1) )

       ! The *_mc_zt microphysics tendencies are passed out of SILHS and back
       ! to CLUBB without being used at all in the rest of the host model code.
       ! The arrays aren't flipped for the *_mc_zt microphysics tendencies, and
       ! they don't need to be.

       ! CLUBB used pverp vertical levels, but SILHS only uses
       ! pverp - top_lev + 1 vertical levels.
       ! Fill the upper levels with 0s when necessary.
       if ( pverp > pverp-top_lev+1 ) then
          rtp2_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          thlp2_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          wprtp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          wpthlp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
          rtpthlp_mc_zt(igrdcol,pverp-top_lev+2:pverp) = 0.0_r8
       endif ! pverp > pverp-top_lev+1

     enddo ! igrdcol = 1, ngrdcol
#endif
#endif

     return
   end subroutine subcol_SILHS_var_covar_driver
#ifdef SILHS
   real(r8) function meansc(arr_in, w_in, ns) result(val)
      real(r8), intent(in) :: ns                         ! Length of Array
      real(r8), dimension(int(ns)), intent(in) :: arr_in      ! Input array
      real(r8), dimension(int(ns)), intent(in) :: w_in        ! Weights
      real(r8) :: acc  ! accumulator
      integer :: i
      acc = 0
      val = 0
      do i=1,ns
         acc = acc + arr_in(i)*w_in(i)
      enddo
      val = acc/ns
   end function

   real(r8) function stdsc(arr_in, w_in, mn_in, ns) result(val)
      real(r8), intent(in) :: ns  ! Number of elements (subcolumns)
      real(r8), dimension(int(ns)), intent(in) :: arr_in, w_in  !Input array and weights
      real(r8), intent(in) :: mn_in   ! The mean of arr_in
      real(r8) :: accvar, var
      integer :: i
      accvar = 0
      do i=1,ns
         accvar = accvar + ((arr_in(i)-mn_in)**2)*w_in(i)
      enddo
      var = accvar/ns
      val = sqrt(var)
   end function

   subroutine Abs_Temp_profile(nz, LWPT_prof, ex_prof, rcm_prof, ABST_prof)

      use clubb_api_module,              only : thlm2T_in_K_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: LWPT_prof  ! Temp prof in LWPT
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: ABST_prof  ! Abs Temp prof
      integer :: i
 
      do i=1,nz
         ABST_prof(i) = thlm2T_in_K_api(LWPT_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine THL_profile(nz, ABST_prof, ex_prof, rcm_prof, THL_prof)

      use clubb_api_module,              only : T_in_K2thlm_api

      integer,                 intent(in)  :: nz         ! Num vert levels
      real(r8), dimension(nz), intent(in)  :: ABST_prof  ! Abs Temp prof
      real(r8), dimension(nz), intent(in)  :: ex_prof    ! Profile of Exner func
      real(r8), dimension(nz), intent(in)  :: rcm_prof   ! Profile of Cld Wat MR
      real(r8), dimension(nz), intent(out) :: THL_prof  ! LWPT prof
      integer :: i
 
      do i=1,nz
         THL_prof(i) = T_in_K2thlm_api(ABST_prof(i), ex_prof(i), rcm_prof(i))
      enddo
      
   end subroutine

   subroutine StaticEng_profile(nz, ABST_prof, zm_prof, zsfc, s_prof)
      integer,                 intent(in) :: nz
      real(r8), dimension(nz), intent(in) :: ABST_prof
      real(r8), dimension(nz), intent(in) :: zm_prof
      real(r8),                intent(in) :: zsfc
      real(r8), dimension(nz), intent(out) :: s_prof
      integer :: i

      do i=1,nz
         s_prof(i) = cpair*(ABST_prof(i)) + gravit*zm_prof(i)+zsfc
      enddo

   end subroutine

   subroutine subcol_constrainmn( num_subcols, samples, weights, grid_mean, mean_sc, std_sc )

      ! Input/Output Variables
      integer, intent(in) :: num_subcols
      real(r8), dimension(num_subcols, pverp), intent(inout) :: samples
      real(r8), dimension(num_subcols), intent(in) :: weights
      real(r8), dimension(pverp), intent(in) :: grid_mean
      real(r8), dimension(pver), intent(out), optional :: mean_sc, std_sc

      ! Local Variables
      real(r8) :: meansc_loc, adj_rat
      integer :: k
   !------------------------------------------------------------------
      !----- Begin Code -----
      do k=1, pver
         meansc_loc = meansc( samples(:,k), weights(:), real(num_subcols, r8) )

         if (present(mean_sc)) &
            mean_sc(k) = meansc_loc
         if (present(std_sc)) &
            std_sc(k) = stdsc( samples(:,k), weights(:), meansc_loc, &
                               real(num_subcols, r8) )

         if ( meansc_loc > 0.0_r8 ) then
            adj_rat = grid_mean(k)/meansc_loc
         else 
            ! If the mean is zero, then zero out all subcolumns to avoid
            ! negative samples
            adj_rat = 0.0_r8
         end if
         samples(:,k) = samples(:,k) * adj_rat
      end do
   end subroutine subcol_constrainmn

   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !
   function clubb_flip_grid ( profile ) result( profile_flipped )

     ! Description:
     !   Swaps the elements in profile so they are in reverse order. CAM and
     !   CLUBB's grids are flipped with respect to each other.
     !
     !   Usage:
     !     clubb_var = clubb_flip_grid( cam_var )
     !     cam_var   = clubb_flip_grid( clubb_var )

     implicit none

     ! Input Variable
     real(r8), dimension(pverp), intent(in) :: profile

     ! Output Variable
     real(r8), dimension(pverp) :: profile_flipped

     ! Local Variable
     integer :: k

     do k=1, pverp
       profile_flipped(k) = profile(pverp-k+1)
     end do ! k=1, pverp

     return
   end function clubb_flip_grid
   ! =============================================================================== !
   !                                                                                 !
   ! =============================================================================== !
#endif
   !============================================================================
   subroutine subcol_SILHS_fill_holes_conserv( state, dt, ptend, pbuf )

     ! The William F. Buckley Jr. Conservative Hole Filler.

     ! Description:
     ! Stops holes from forming in a hydrometeor mixing ratio by reducing the
     ! microphysics tendency of that hydrometeor mixing ratio which would
     ! otherwise cause that hydrometeor mixing ratio to have a negative value
     ! once the microphysics tendency is applied.  This code is used to prevent
     ! holes in water mass, not number concentration.
     !
     ! This subroutine is called after microphysics has completed and after
     ! microphysics fields from subcolumns have been averaged back to grid
     ! columns, but before the grid-column microphysics tendencies have been
     ! applied in physics_update.  This code is meant for use with the SILHS
     ! subcolumn approach.  This code needs to be applied to grid columns, not
     ! subcolumns.
     !
     ! This code adjusts the tendencies (ptend) before they are used to update
     ! the grid mean fields (state variables).
     !
     ! The column-integrated total water needs to be conserved during
     ! microphysics.  The conserved amount includes the amount of water that
     ! precipitated to the ground from sedimentation during microphysics.
     ! The conservation equation for each grid column is:
     !
     ! SUM(k=top_lev:pver) ( rv_start(k) + rc_start(k) + rr_start(k)
     !                       + ri_start(k) + rs_start(k) ) * pdel(k) / g
     ! = SUM(k=top_lev:pver) ( rv(k) + rc(k) + rr(k) + ri(k) + rs(k) )
     !                       * pdel(k) / g
     !   + prect * dt * 1000;
     !
     ! where rv_start, rc_start, rr_start, ri_start, and rs_start are water
     ! vapor, cloud water, rain water, cloud ice, and snow mixing ratios before
     ! microphysics is called; rv, rc, rr, ri, and rs are water vapor, cloud
     ! water, rain water, cloud ice, and snow mixing ratios after being updated
     ! by microphysics; pdel is the pressure difference between vertical levels,
     ! g is gravity, and prect * dt * 1000 is the total amount of water (from
     ! all precipitating hydrometeors) that sedimented to the ground during
     ! microphysics (dt is the timestep used for microphysics).  The units of
     ! column-integrated total water are kg (water) / m^2.
     !
     ! All the updated hydrometeor fields are related to the hydrometeor fields
     ! at the start by:
     !
     ! rv(k) = rv_start(k) + rv_tend(k) * dt;
     ! rc(k) = rc_start(k) + rc_tend(k) * dt;
     ! rr(k) = rr_start(k) + rr_tend(k) * dt;
     ! ri(k) = ri_start(k) + ri_tend(k) * dt; and
     ! rs(k) = rs_start(k) + rs_tend(k) * dt;
     !
     ! where rv_tend, rc_tend, rr_tend, ri_tend, and rs_tend are water vapor,
     ! cloud water, rain water, cloud ice, and snow mixing ratio tendencies
     ! from microphysics, which includes the sum of microphysics process rates
     ! and sedimentation.  When these equations are applied to the equation
     ! for column-integrated total water, that equation becomes:
     !
     ! SUM(k=top_lev:pver) ( rv_tend(k) + rc_tend(k) + rr_tend(k)
     !                       + ri_tend(k) + rs_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! As stated above, the hydrometeor tendencies are the sum of tendencies
     ! from microphysics process rates and tendencies from sedimentation:
     !
     ! rv_tend(k) = rv_mc_tend(k);
     ! rc_tend(k) = rc_mc_tend(k) + rc_sed_tend(k);
     ! rr_tend(k) = rr_mc_tend(k) + rr_sed_tend(k);
     ! ri_tend(k) = ri_mc_tend(k) + ri_sed_tend(k); and
     ! rs_tend(k) = rs_mc_tend(k) + rs_sed_tend(k);
     !
     ! where rv_mc_tend, rc_mc_tend, rr_mc_tend, ri_mc_tend, and rs_mc_tend are
     ! the tendencies of water vapor, cloud water, rain water, cloud ice, and
     ! snow from microphysics process rates, and rc_sed_tend, rr_sed_tend,
     ! ri_sed_tend, and rs_sed_tend are the tendencies of cloud water,
     ! rain water, cloud ice, and snow from sedimentation.  When these equations
     ! are applied to the equation for column-integrated total water, that
     ! equation becomes:
     !
     ! SUM(k=top_lev:pver) ( rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     !                       + ri_mc_tend(k) + rs_mc_tend(k) )
     !                     * dt * pdel(k) / g
     ! + SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) + ri_sed_tend(k)
     !                         + rs_sed_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! At any vertical level, the tendencies from microphysics process rates
     ! (mc_tend variables) must balance:
     !
     ! rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     ! + ri_mc_tend(k) + rs_mc_tend(k) = 0; for all k from top_lev to pver.
     !
     ! The column-integrated total water equation can be applied to
     ! sedimentation:
     !
     ! SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) + ri_sed_tend(k)
     !                       + rs_sed_tend(k) ) * dt * pdel(k) / g
     ! + prect * dt * 1000 = 0.
     !
     ! The total precipitation rate, prect, can be split into liquid
     ! precipitation rate, precl, and frozen precipitation rate, preci:
     !
     ! prect = precl + preci.
     !
     ! The microphysics code outputs prect and preci, so precl can be calculated
     ! by precl = prect - preci.  The column-integrated total water equation can
     ! be split into:
     !
     ! SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
     !                     * dt * pdel(k) / g
     ! + precl * dt * 1000 = 0; and
     !
     ! SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
     !                     * dt * pdel(k) / g
     ! + preci * dt * 1000 = 0.
     !
     ! Overall, the conservation methods used in this subroutine are:
     !
     ! 1) When adjusting the tendencies from microphysics process rates,
     !    conserve:
     !
     !    rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     !    + ri_mc_tend(k) + rs_mc_tend(k) = 0; for all k from top_lev to pver.
     !
     ! 2) When adjusting the tendencies from microphysics process rates, adjust
     !    dry static energy appropriately.  The change in dry static energy
     !    is necessary because of phase changes.  This "puts back" the extra dry
     !    static energy that was "taken out" when an excessive phase-changing
     !    process rate was produced by microphysics.
     !
     ! 3) When adjusting the hydrometeor tendency from sedimentation of a
     !    liquid hydrometeor (cloud water or rain water), conserve:
     !
     !    SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + precl * dt * 1000 = 0.
     !
     ! 4) When adjusting the hydrometeor tendency from sedimentation of a
     !    frozen hydrometeor (cloud ice or snow), conserve:
     !
     !    SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + preci * dt * 1000 = 0.
     !
     ! The conservative hole filler works as follows.  The total microphysics
     ! tendency for each hydrometeor is provided in ptend.  This is the sum of
     ! the microphysics process rate tendency and sedimentation tendency for
     ! each hydrometeor.  The sedimentation tendency is provided in pbuf.  The
     ! sedimentation tendency is subtracted off the total microphysics tendency
     ! to produce the microphysics process rate tendency for each hydrometeor.
     ! The microphysics process rate tendency is adjusted when necessary so that
     ! holes in the hydrometeor are not produced by microphysics process rates.
     ! When a hydrometeor's negative microphysics process rate tendency needs to
     ! be made smaller in magnitude to avoid a hole, all hydrometeor tendencies
     ! that are positive at that grid level are also decreased proportionately
     ! to maintain a balance.  Dry static energy tendency is also adjusted
     ! appropriately when necessary.  After this, the vertical integral of each
     ! hydrometeor species is greater than or equal to 0.
     !
     ! The sedimentation tendency is then added back onto the new microphysics
     ! process rate tendency to produce a new total microphysics tendency for
     ! each hydrometeor.  Since the sedimentation tendency was based on the old
     ! value of hydrometeor, before the hole-filling adjustment, it is possible
     ! that the new total microphysics tendency may produce holes.  When this
     ! happens, sedimentation hole filling fills holes in the vertical profile
     ! of each hydrometeor.  Holes are filled using mass from other vertical
     ! levels for the same hydrometeor (or from a same-phase hydrometeor when
     ! necessary).  Since the vertical integral of sedimentation tendency
     ! (including surface precipitation rate) is 0, the vertical integral of the
     ! hydrometeor must be greater than or equal to 0, which means that all
     ! holes can be filled.  The result is that all holes in any hydrometeor
     ! mixing ratio are filled completely and conservatively.  The value of
     ! ptend is updated appropriately so that it can be applied later in
     ! physics_update.

     !----------------------------------------------------------------------

     use physics_buffer, only: &
         physics_buffer_desc, &
         pbuf_get_field

     use ppgrid, only: &
         pcols

     use constituents, only: &
         qmin

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     type(physics_state), intent(in) :: state     ! Physics state variables
     real(r8), intent(in) :: dt                   ! Time step duration

     ! Input/Output Variables
     type(physics_ptend),  intent(inout) :: ptend  ! Parameterization tendencies
     type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer

     ! Local Variables
     real(r8), dimension(pcols,pver) :: &
       rv_start, & ! Water vapor mixing ratio at start of microphysics  [kg/kg]
       rc_start, & ! Cloud water mixing ratio at start of microphysics  [kg/kg]
       rr_start, & ! Rain water mixing ratio at start of microphysics   [kg/kg]
       ri_start, & ! Cloud ice mixing ratio at start of microphysics    [kg/kg]
       rs_start    ! Snow mixing ratio at start of microphysics         [kg/kg]

     real(r8), dimension(pcols,pver) :: &
       rv_tend, & ! Water vapor mixing ratio tendency  [kg/kg/s]
       rc_tend, & ! Cloud water mixing ratio tendency  [kg/kg/s]
       rr_tend, & ! Rain water mixing ratio tendency   [kg/kg/s]
       ri_tend, & ! Cloud ice mixing ratio tendency    [kg/kg/s]
       rs_tend, & ! Snow mixing ratio tendency         [kg/kg/s]
       stend      ! Dry static energy tendency         [J/kg/s]

     real(r8), dimension(:), pointer :: &
       prect,    & ! Total microphysics precipitation rate (surface)      [m/s]
       preci,    & ! Ice-phase microphysics precipitation rate (surface)  [m/s]
       prec_str, & ! Total surface precipitation rate from stratoform     [m/s]
       snow_str    ! Snow surface precipitation rate from stratoform      [m/s]

     real(r8), dimension(:,:), pointer :: &
       rc_sed_tend, & ! Mean cloud water sedimentation tendency        [kg/kg/s]
       rr_sed_tend, & ! Mean rain water sedimentation tendency         [kg/kg/s]
       ri_sed_tend, & ! Mean cloud ice sedimentation tendency          [kg/kg/s]
       rs_sed_tend, & ! Mean snow sedimentation tendency               [kg/kg/s]
       vtrmc,       & ! Mean cloud water sedimentation velocity        [m/s]
       umr,         & ! Mean rain water sedimentation velocity         [m/s]
       vtrmi,       & ! Mean cloud ice sedimentation velocity          [m/s]
       ums,         & ! Mean snow sedimentation velocity               [m/s]
       rc_sed_evap, & ! Mean evap of cloud water during sedimentation  [kg/kg/s]
       ri_sed_subl    ! Mean subl of cloud ice during sedimentation    [kg/kg/s]

     real(r8), dimension(pcols,pver) :: &
       rv_mc_tend, & ! Water vapor mixing ratio microphysics tendency  [kg/kg/s]
       rc_mc_tend, & ! Cloud water mixing ratio microphysics tendency  [kg/kg/s]
       rr_mc_tend, & ! Rain water mixing ratio microphysics tendency   [kg/kg/s]
       ri_mc_tend, & ! Cloud ice mixing ratio microphysics tendency    [kg/kg/s]
       rs_mc_tend    ! Snow mixing ratio microphysics tendency         [kg/kg/s]

     real(r8) :: &
       rv_curr, & ! Current water vapor mixing ratio    [kg/kg]
       rc_curr, & ! Current cloud water mixing ratio    [kg/kg]
       rr_curr, & ! Current rain water mixing ratio     [kg/kg]
       ri_curr, & ! Current cloud ice mixing ratio      [kg/kg]
       rs_curr    ! Current snow mixing ratio           [kg/kg]

     logical :: &
       l_pos_rv_mc_tend, & ! Flag for positive water vapor mixing ratio mc tend.
       l_pos_rc_mc_tend, & ! Flag for positive cloud water mixing ratio mc tend.
       l_pos_rr_mc_tend, & ! Flag for positive rain water mixing ratio mc tend.
       l_pos_ri_mc_tend, & ! Flag for positive cloud ice mixing ratio mc tend.
       l_pos_rs_mc_tend    ! Flag for positive snow mixing ratio mc tend.

     real(r8) :: &
       mc_tend_max_mag,     & ! Max. allowable mag. of (neg.) mc tend [kg/kg/s]
       mc_tend_correction,  & ! Amnt. correction necessary to mc tend [kg/kg/s]
       total_mc_positive,   & ! Total of all positive mc tendencies   [kg/kg/s]
       mc_correction_ratio    ! Ratio: mc_tend_correction/total_mc_positive [-]

     real(r8), dimension(pcols) :: &
       precl    ! Liquid-phase precipitation rate (surface)        [m/s]

     ! Budgeting terms for hole filling.
     ! These variables are for use in stats output.
     real(r8), dimension(pcols,pver) :: &
       rv_hf_tend, & ! Water vapor mixing ratio hole-filling tendency  [kg/kg/s]
       rc_hf_tend, & ! Cloud water mixing ratio hole-filling tendency  [kg/kg/s]
       rr_hf_tend, & ! Rain water mixing ratio hole-filling tendency   [kg/kg/s]
       ri_hf_tend, & ! Cloud ice mixing ratio hole-filling tendency    [kg/kg/s]
       rs_hf_tend, & ! Snow mixing ratio hole-filling tendency         [kg/kg/s]
       s_hf_tend     ! Dry static energy hole-filling tendency         [J/kg/s]

     integer :: ncol  ! Number of grid columns

     integer :: icol, k  ! Loop indices

     ! Flag to perform hole filling after the original sedimentation tendency
     ! is added back on to the new microphysics process tendency.  This calls
     ! the sedimentation hole filler.
     logical, parameter :: &
       l_sed_hole_fill = .true.

     logical, parameter :: &
       l_check_conservation = .true. ! Flag to perform water conservation check

     ! Vertically-integrated grand total water (rv+rc+rr+ri+rs)    [kg/m^2]
     real(r8), dimension(pcols) :: &
       grand_total_water_column_start,  & ! Column integral at start
       grand_total_water_column_finish    ! Column integral at finish

     ! Vertically-integrated total water energy    [J/m^2]
     real(r8), dimension(pcols) :: &
       total_energy_column_start,  & ! Column integral at start
       total_energy_column_finish    ! Column integral at finish

     real(r8), dimension(pcols) :: &
       tot_water_rel_err,  & ! Relative error: vert-integrated grand total water
       tot_energy_rel_err    ! Relative error: vert-integrated total energy

     real(r8), parameter :: &
       err_thresh = 1.0e-14_r8  ! Threshold of relative error


     ! Get the number of grid columns.
     ncol = state%ncol

     ! Get fields from the pbuf.
     call pbuf_get_field(pbuf, prec_pcw_idx, prect)
     call pbuf_get_field(pbuf, snow_pcw_idx, preci)
     call pbuf_get_field(pbuf, prec_str_idx, prec_str)
     call pbuf_get_field(pbuf, snow_str_idx, snow_str)
     call pbuf_get_field(pbuf, qcsedten_idx, rc_sed_tend)
     call pbuf_get_field(pbuf, qrsedten_idx, rr_sed_tend)
     call pbuf_get_field(pbuf, qisedten_idx, ri_sed_tend)
     call pbuf_get_field(pbuf, qssedten_idx, rs_sed_tend)
     call pbuf_get_field(pbuf, vtrmc_idx, vtrmc)
     call pbuf_get_field(pbuf, umr_idx, umr)
     call pbuf_get_field(pbuf, vtrmi_idx, vtrmi)
     call pbuf_get_field(pbuf, ums_idx, ums)
     call pbuf_get_field(pbuf, qcsevap_idx, rc_sed_evap)
     call pbuf_get_field(pbuf, qisevap_idx, ri_sed_subl)

     ! Calculate liquid precipitation rate (precl) from the total precipitation
     ! rate (prect) and the frozen preciptation rate (preci).  This should never
     ! be negative, but just to be safe, threshold at 0.
     precl(:ncol) = max( prect(:ncol) - preci(:ncol), 0.0_r8 )

     ! Perform total water and total energy conservation checks.
     if ( l_check_conservation ) then

        ! Calculate total water in each column.
        ! This calculation is the vertically-integrated grand total water (where
        ! grand total water is the sum of water vapor, cloud water, rain water,
        ! cloud ice, and snow, as well as the amount of water that precipitated
        ! to the surface) in each grid column after microphysics, but at the
        ! start of hole filling.
        do icol = 1, ncol
           grand_total_water_column_start(icol) = 0.0_r8
           do k = top_lev, pver
              grand_total_water_column_start(icol) &
              = grand_total_water_column_start(icol) &
                + ( state%q(icol,k,1) + ptend%q(icol,k,1) * dt &
                    + state%q(icol,k,ixcldliq) &
                    + ptend%q(icol,k,ixcldliq) * dt &
                    + state%q(icol,k,ixcldice) &
                    + ptend%q(icol,k,ixcldice) * dt ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 grand_total_water_column_start(icol) &
                 = grand_total_water_column_start(icol) &
                   + ( state%q(icol,k,ixrain) + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
              if ( ixsnow > 0 ) then
                 grand_total_water_column_start(icol) &
                 = grand_total_water_column_start(icol) &
                   + ( state%q(icol,k,ixsnow) + ptend%q(icol,k,ixsnow) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
           grand_total_water_column_start(icol) &
           = grand_total_water_column_start(icol) &
             + prect(icol) * dt * 1000.0_r8
        enddo ! icol = 1, ncol

        ! Calculate total energy in each column.
        ! This calculation is the vertically-integrated total energy in each
        ! grid column after microphysics, but at the start of hole filling.
        ! Since the microphysics and hole filling code does not directly change
        ! kinetic energy, 0.5 * ( u^2 + v^2 ), it can be skipped as part of the
        ! energy conservation check.
        do icol = 1, ncol
           total_energy_column_start(icol) = 0.0_r8
           do k = top_lev, pver
              total_energy_column_start(icol) &
              = total_energy_column_start(icol) &
                + ( state%s(icol,k) + ptend%s(icol,k) * dt &
                    + ( latvap + latice ) &
                      * ( state%q(icol,k,1) + ptend%q(icol,k,1) * dt ) &
                    + latice * ( state%q(icol,k,ixcldliq) &
                                 + ptend%q(icol,k,ixcldliq) * dt ) ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 total_energy_column_start(icol) &
                 = total_energy_column_start(icol) &
                   + latice * ( state%q(icol,k,ixrain) &
                                + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
           total_energy_column_start(icol) &
           = total_energy_column_start(icol) &
             + latice * precl(icol) * dt * 1000.0_r8
        enddo ! icol = 1, ncol

     endif ! l_check_conservation

     ! The fields within state haven't been updated yet, since this is before
     ! the call to physics_update.
     rv_start = state%q(:,:,1)
     rc_start = state%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_start = state%q(:,:,ixrain)
     endif
     ri_start = state%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_start = state%q(:,:,ixsnow)
     endif

     ! Unpack the current total tendencies for hydrometeor mixing ratio fields.
     rv_tend = ptend%q(:,:,1)
     rc_tend = ptend%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_tend = ptend%q(:,:,ixrain)
     endif
     ri_tend = ptend%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_tend = ptend%q(:,:,ixsnow)
     endif

     ! Unpack the current tendency for dry static energy.
     stend = ptend%s

     ! The total hydrometeor tendencies are the sum of microphysics process
     ! rates and sedimentation rates.  Calculate the microphysics process
     ! tendencies by subtracting the sedimentation tendencies from the overall
     ! tendencies.
     ! The sedimentation tendencies for cloud water (rc_sed_tend) and cloud ice
     ! (ri_sed_tend) include the evaporation of cloud water during sedimentation
     ! and the sublimation of cloud ice during sedimentation, respectively.  The
     ! true sedimentation of cloud water is the sum of rc_sed_tend and
     ! rc_sed_evap, and the true sedimentation of cloud ice is the sum of
     ! ri_sed_tend and ri_sed_subl.  Subtract off only the true sedimentation
     ! rates, as evaporation and sublimation need to be included in the
     ! microphysics process rates.
     rv_mc_tend(:ncol,:) = rv_tend(:ncol,:)
     rc_mc_tend(:ncol,:) = rc_tend(:ncol,:) - ( rc_sed_tend(:ncol,:) + rc_sed_evap(:ncol,:) )
     if ( ixrain > 0 ) then
        rr_mc_tend(:ncol,:) = rr_tend(:ncol,:) - rr_sed_tend(:ncol,:)
     endif
     ri_mc_tend(:ncol,:) = ri_tend(:ncol,:) - ( ri_sed_tend(:ncol,:) + ri_sed_subl(:ncol,:) )
     if ( ixsnow > 0 ) then
        rs_mc_tend(:ncol,:) = rs_tend(:ncol,:) - rs_sed_tend(:ncol,:)
     endif

     ! This section adjusts microphysics process rate tendencies so that the
     ! resulting values of all hydrometeor mixing ratios are greater than or
     ! equal to qmin after this section is complete.  Once sedimentation is
     ! added back on after this section, some of the hydrometeor mixing ratios
     ! may become less than qmin again.
     !
     ! This section, which again is concerned only with adjusting microphysics
     ! process rates, makes use of the following two principles:
     !
     ! 1) When adjusting the tendencies from microphysics process rates,
     !    conserve:
     !
     !    rv_mc_tend(k) + rc_mc_tend(k) + rr_mc_tend(k)
     !    + ri_mc_tend(k) + rs_mc_tend(k) = 0; for all k from top_lev to pver.
     !
     ! 2) When adjusting the tendencies from microphysics process rates, adjust
     !    dry static energy appropriately.  The change in dry static energy
     !    is necessary because of phase changes.  This "puts back" the extra dry
     !    static energy that was "taken out" when an excessive phase-changing
     !    process rate was produced by microphysics.

     ! Loop over all columns, performing any tendency adjustments one column
     ! at a time.
     do icol = 1, ncol

        ! Loop over all vertical levels, performing any microphysics process
        ! tendency adjustments one level at a time.
        do k = top_lev, pver

           ! Find which hydrometeors have positive microphysics process
           ! tendencies at this level.
           if ( rv_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_rv_mc_tend = .true.
           else
              l_pos_rv_mc_tend = .false.
           endif
           if ( rc_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_rc_mc_tend = .true.
           else
              l_pos_rc_mc_tend = .false.
           endif
           if ( ixrain > 0 ) then
              if ( rr_mc_tend(icol,k) >= 0.0_r8 ) then
                 l_pos_rr_mc_tend = .true.
              else
                 l_pos_rr_mc_tend = .false.
              endif
           endif
           if ( ri_mc_tend(icol,k) >= 0.0_r8 ) then
              l_pos_ri_mc_tend = .true.
           else
              l_pos_ri_mc_tend = .false.
           endif
           if ( ixsnow > 0 ) then
              if ( rs_mc_tend(icol,k) >= 0.0_r8 ) then
                 l_pos_rs_mc_tend = .true.
              else
                 l_pos_rs_mc_tend = .false.
              endif
           endif

           !!! Check for holes in water vapor mixing ratio
           if ( .not. l_pos_rv_mc_tend ) then

              ! Calculate the water vapor mixing ratio as it would be with the
              ! current microphysics process tendency.
              rv_curr = rv_start(icol,k) + rv_mc_tend(icol,k) * dt

              if ( rv_curr < qmin(1) ) then

                 ! Microphysics processes are causing a hole in water vapor
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) water
                 ! vapor microphysics process tendency.
                 mc_tend_max_mag = ( qmin(1) - rv_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the water vapor mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rv_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principle, this should never be greater than 1 outside of
                 ! numerical round-off errors.  This is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30_r8 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative water vapor mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of water vapor).
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latvap * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latvap * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to water vapor cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - ( latvap + latice ) &
                        * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to water vapor cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - ( latvap + latice ) &
                        * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new water vapor mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rv_mc_tend(icol,k) &
                 = rv_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rv_curr < qmin(1)

           endif ! .not. l_pos_rv_mc_tend

           !!! Check for holes in cloud water mixing ratio
           if ( .not. l_pos_rc_mc_tend ) then

              ! Calculate the cloud water mixing ratio as it would be with the
              ! current microphysics process tendency.
              rc_curr = rc_start(icol,k) + rc_mc_tend(icol,k) * dt

              if ( rc_curr < qmin(ixcldliq) ) then

                 ! Microphysics processes are causing a hole in cloud water
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) cloud
                 ! water microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixcldliq) - rc_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the cloud water mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rc_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principle, this should never be greater than 1 outside of
                 ! numerical round-off errors.  This is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30_r8 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative cloud water mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of cloud water).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to cloud water heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latvap * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to cloud water does not change
                    ! dry static energy.
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to cloud water cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to cloud water cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new cloud water mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rc_mc_tend(icol,k) &
                 = rc_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rc_curr < qmin(ixcldliq)

           endif ! .not. l_pos_rc_mc_tend

           !!! Check for holes in rain water mixing ratio
           if ( ixrain > 0 .and. ( .not. l_pos_rr_mc_tend ) ) then

              ! Calculate the rain water mixing ratio as it would be with the
              ! current microphysics process tendency.
              rr_curr = rr_start(icol,k) + rr_mc_tend(icol,k) * dt

              if ( rr_curr < qmin(ixrain) ) then

                 ! Microphysics processes are causing a hole in rain water
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) rain
                 ! water microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixrain) - rr_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the rain water mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rr_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principle, this should never be greater than 1 outside of
                 ! numerical round-off errors.  This is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30_r8 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative rain water mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of rain water).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to rain water heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latvap * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to rain water does not change
                    ! dry static energy.
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to rain water cools and reduces
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * ri_mc_tend(icol,k)
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to rain water cools and reduces dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      - latice * mc_correction_ratio * rs_mc_tend(icol,k)
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new rain water mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 rr_mc_tend(icol,k) &
                 = rr_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rr_curr < qmin(ixrain)

           endif ! ixrain > 0 .and. ( .not. l_pos_rr_mc_tend )

           !!! Check for holes in cloud ice mixing ratio
           if ( .not. l_pos_ri_mc_tend ) then

              ! Calculate the cloud ice mixing ratio as it would be with the
              ! current microphysics process tendency.
              ri_curr = ri_start(icol,k) + ri_mc_tend(icol,k) * dt

              if ( ri_curr < qmin(ixcldice) ) then

                 ! Microphysics processes are causing a hole in cloud ice
                 ! mixing ratio.

                 ! Calculate the maximum allowable magnitude of (negative) cloud
                 ! ice microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixcldice) - ri_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the cloud ice mixing ratio microphysics process
                 ! tendency.  This number is positive.
                 mc_tend_correction = mc_tend_max_mag - ri_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    total_mc_positive = total_mc_positive + rs_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principle, this should never be greater than 1 outside of
                 ! numerical round-off errors.  This is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30_r8 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative cloud ice mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of cloud ice).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + ( latvap + latice ) &
                        * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to cloud ice heats and increases
                    ! dry static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixsnow > 0 .and. l_pos_rs_mc_tend ) then
                    ! Changing snow to cloud ice does not change dry static
                    ! energy.
                    ! Update snow mixing ratio microphysics tendency.
                    rs_mc_tend(icol,k) &
                    = rs_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new cloud ice mixing ratio microphysics
                 ! process tendency.  This should be equal to the maximum
                 ! magnitude (negative) amount allowed, mc_tend_max_mag.
                 ri_mc_tend(icol,k) &
                 = ri_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! ri_curr < qmin(ixcldice)

           endif ! .not. l_pos_ri_mc_tend

           !!! Check for holes in snow mixing ratio
           if ( ixsnow > 0 .and. ( .not. l_pos_rs_mc_tend ) ) then

              ! Calculate the snow mixing ratio as it would be with the
              ! current microphysics process tendency.
              rs_curr = rs_start(icol,k) + rs_mc_tend(icol,k) * dt

              if ( rs_curr < qmin(ixsnow) ) then

                 ! Microphysics processes are causing a hole in snow mixing
                 ! ratio.

                 ! Calculate the maximum allowable magnitude of (negative) snow
                 ! microphysics process tendency.
                 mc_tend_max_mag = ( qmin(ixsnow) - rs_start(icol,k) ) / dt

                 ! Calculate the amount of the correction that needs to be made
                 ! to the snow mixing ratio microphysics process tendency.
                 ! This number is positive.
                 mc_tend_correction = mc_tend_max_mag - rs_mc_tend(icol,k)

                 ! Calculate the total amount of positive microphysics process
                 ! tendencies for all hydrometeor mixing ratios.
                 total_mc_positive = 0.0_r8
                 if ( l_pos_rv_mc_tend ) then
                    total_mc_positive = total_mc_positive + rv_mc_tend(icol,k)
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    total_mc_positive = total_mc_positive + rc_mc_tend(icol,k)
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    total_mc_positive = total_mc_positive + rr_mc_tend(icol,k)
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    total_mc_positive = total_mc_positive + ri_mc_tend(icol,k)
                 endif

                 ! Calculate the correction ratio.
                 ! In principle, this should never be greater than 1 outside of
                 ! numerical round-off errors.  This is limited at 1 to be safe.
                 mc_correction_ratio &
                 = min( mc_tend_correction &
                        / max( total_mc_positive, 1.0e-30_r8 ), 1.0_r8 )

                 ! Adjust (decrease) the tendencies of all positive hydrometeor
                 ! mixing ratio tendencies to balance the adjustment (increase)
                 ! to the excessively negative snow mixing ratio.
                 ! Transfer dry static energy appropriately (in response to the
                 ! excessive depletion of snow).
                 if ( l_pos_rv_mc_tend ) then
                    ! Changing water vapor to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + ( latvap + latice ) &
                        * mc_correction_ratio * rv_mc_tend(icol,k)
                    ! Update water vapor mixing ratio microphysics tendency.
                    rv_mc_tend(icol,k) &
                    = rv_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_rc_mc_tend ) then
                    ! Changing cloud water to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rc_mc_tend(icol,k)
                    ! Update cloud water mixing ratio microphysics tendency.
                    rc_mc_tend(icol,k) &
                    = rc_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( ixrain > 0 .and. l_pos_rr_mc_tend ) then
                    ! Changing rain water to snow heats and increases dry
                    ! static energy.
                    stend(icol,k) &
                    = stend(icol,k) &
                      + latice * mc_correction_ratio * rr_mc_tend(icol,k)
                    ! Update rain water mixing ratio microphysics tendency.
                    rr_mc_tend(icol,k) &
                    = rr_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif
                 if ( l_pos_ri_mc_tend ) then
                    ! Changing cloud ice to snow does not change dry static
                    ! energy.
                    ! Update cloud ice mixing ratio microphysics tendency.
                    ri_mc_tend(icol,k) &
                    = ri_mc_tend(icol,k) * ( 1.0_r8 - mc_correction_ratio )
                 endif

                 ! Calculate the new snow mixing ratio microphysics process
                 ! tendency.  This should be equal to the maximum magnitude
                 ! (negative) amount allowed, mc_tend_max_mag.
                 rs_mc_tend(icol,k) &
                 = rs_mc_tend(icol,k) &
                   + mc_correction_ratio * total_mc_positive

              endif ! rs_curr < qmin(ixsnow)

           endif ! ixsnow > 0 .and. ( .not. l_pos_rs_mc_tend )

        enddo ! k = top_lev, pver

     enddo ! icol = 1, ncol

     ! Calculate the new overall tendencies by adding the sedimentation
     ! tendencies back onto the new microphysics process tendencies.
     ! For cloud water and cloud ice, the sedimentation tendencies that are
     ! added back on are the true sedimentation tendencies.  For cloud water,
     ! this is the sum of rc_sed_tend and rc_sed_evap, and for cloud ice, this
     ! is the sum of ri_sed_tend and ri_sed_subl.
     rv_tend(:ncol,:) = rv_mc_tend(:ncol,:)
     rc_tend(:ncol,:) = rc_mc_tend(:ncol,:) + ( rc_sed_tend(:ncol,:) + rc_sed_evap(:ncol,:) )
     if ( ixrain > 0 ) then
        rr_tend(:ncol,:) = rr_mc_tend(:ncol,:) + rr_sed_tend(:ncol,:)
     endif
     ri_tend(:ncol,:) = ri_mc_tend(:ncol,:) + ( ri_sed_tend(:ncol,:) + ri_sed_subl(:ncol,:) )
     if ( ixsnow > 0 ) then
        rs_tend(:ncol,:) = rs_mc_tend(:ncol,:) + rs_sed_tend(:ncol,:)
     endif

     ! Now that the original sedimentation tendency has been added to the
     ! new microphysics process tendency, the new total microphysics tendency
     ! can still cause holes to form.  After the microphysics process rates were
     ! adjusted, the values of the hydrometeor fields were greater than or equal
     ! to 0 at all grid levels, which means their vertical integrals were also
     ! greater than or equal to 0.  Sedimentation by itself has a vertical
     ! integral of 0 (including the amount that sedimented to the surface).
     ! This means that after the microphysics process rates have been adjusted
     ! and sedimentation has been added back on, the resulting hydrometeor
     ! fields all still have vertical integrals that are greater than or equal
     ! to 0.  Holes that develop at any particular grid level can be filled.
     ! These holes can be filled conservatively using the sedimentation hole
     ! filler.
     if ( l_sed_hole_fill ) then

        ! This section makes use of the following principle:
        !
        ! 3) When adjusting the hydrometeor tendency from sedimentation of a
        !    liquid hydrometeor (cloud water or rain water), conserve:
        !
        !    SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
        !                        * dt * pdel(k) / g
        !    + precl * dt * 1000 = 0.

        ! Call the sedimentation hole filler for rain water mixing ratio.
        ! This can update rr_tend and precl.
        if ( ixrain > 0 ) then
           call fill_holes_sedimentation( dt, ncol, rr_start, state%pdel, &
                                          umr, state%zi, qmin(ixrain), &
                                          rr_tend, precl )
        endif ! ixrain > 0

        ! Call the sedimentation hole filler for cloud water mixing ratio.
        ! This can update rc_tend and precl.
        call fill_holes_sedimentation( dt, ncol, rc_start, state%pdel, &
                                       vtrmc, state%zi, qmin(ixcldliq), &
                                       rc_tend, precl )

        ! Occasionally, a situation can occur where filling a hole in rain can
        ! deplete all the surface liquid-phase precipitation (precl), resulting
        ! in not enough water mass in the vertical profile of cloud water to
        ! fill a hole in cloud water.  When this happens, there must be liquid
        ! water found in the vertical profile of rain, so pull the water from
        ! rain to fill any remaining holes in cloud water.
        if ( ixrain > 0 ) then
           call fill_holes_same_phase_vert( dt, ncol, rc_start, rr_start, &
                                            state%pdel, qmin(ixcldliq), &
                                            qmin(ixrain), &
                                            rc_tend, rr_tend )
        endif ! ixrain > 0

        ! This section makes use of the following principle:
        !
        ! 4) When adjusting the hydrometeor tendency from sedimentation of a
        !    frozen hydrometeor (cloud ice or snow), conserve:
        !
        !    SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
        !                        * dt * pdel(k) / g
        !    + preci * dt * 1000 = 0.

        ! Call the sedimentation hole filler for snow mixing ratio.
        ! This can update rs_tend and preci.
        if ( ixsnow > 0 ) then
           call fill_holes_sedimentation( dt, ncol, rs_start, state%pdel, &
                                          ums, state%zi, qmin(ixsnow), &
                                          rs_tend, preci )
        endif ! ixsnow > 0

        ! Call the sedimentation hole filler for cloud ice mixing ratio.
        ! This can update ri_tend and preci.
        call fill_holes_sedimentation( dt, ncol, ri_start, state%pdel, &
                                       vtrmi, state%zi, qmin(ixcldice), &
                                       ri_tend, preci )

        ! Occasionally, a situation can occur where filling a hole in snow can
        ! deplete all the surface ice-phase precipitation (preci), resulting
        ! in not enough water mass in the vertical profile of cloud ice to
        ! fill a hole in cloud ice.  When this happens, there must be ice-phase
        ! water found in the vertical profile of snow, so pull the water from
        ! snow to fill any remaining holes in cloud ice.
        if ( ixsnow > 0 ) then
           call fill_holes_same_phase_vert( dt, ncol, ri_start, rs_start, &
                                            state%pdel, qmin(ixcldice), &
                                            qmin(ixsnow), &
                                            ri_tend, rs_tend )
        endif  ! ixsnow > 0

        ! Update the total precipitation rate (prect) from the updated liquid
        ! precipitation rate (precl) and the updated frozen preciptation rate
        ! (preci).
        prect(:ncol) = precl(:ncol) + preci(:ncol)

        ! The MG code sets prec_str equal to prect (prec_pcw) and snow_str equal
        ! to preci (snow_pcw).  The prec_str and snow_str variables are used
        ! in the calculations for energy and water conservation.  Since prect
        ! and preci are adjusted here, when necessary, prec_str and snow_str
        ! also need to be adjusted.
        prec_str(:ncol) = prect(:ncol)
        snow_str(:ncol) = preci(:ncol)

     endif ! l_sed_hole_fill

     ! The updated total microphysics tendencies after hole filling have not
     ! been used to update ptend yet, so record the budget terms for hole
     ! filling first.
     rv_hf_tend = rv_tend - ptend%q(:,:,1)
     rc_hf_tend = rc_tend - ptend%q(:,:,ixcldliq)
     if ( ixrain > 0 ) then
        rr_hf_tend = rr_tend - ptend%q(:,:,ixrain)
     endif ! ixrain > 0
     ri_hf_tend = ri_tend - ptend%q(:,:,ixcldice)
     if ( ixsnow > 0 ) then
        rs_hf_tend = rs_tend - ptend%q(:,:,ixsnow)
     endif ! ixsnow > 0

     ! The updated dry static energy tendency after hole filling has not been
     ! used to update ptend yet, so record the budget term for hole filling
     ! first.
     s_hf_tend = stend - ptend%s

     ! Pack the current total tendencies for hydrometeor mixing ratio fields.
     ptend%q(:,:,1) = rv_tend
     ptend%q(:,:,ixcldliq) = rc_tend
     if ( ixrain > 0 ) then
        ptend%q(:,:,ixrain) = rr_tend
     endif
     ptend%q(:,:,ixcldice) = ri_tend
     if ( ixsnow > 0 ) then
        ptend%q(:,:,ixsnow) = rs_tend
     endif

     ! Pack the current tendency for dry static energy.
     ptend%s = stend

     ! Output stats for hole filling tendencies.
     call outfld( 'QVHFTEN', rv_hf_tend, pcols, state%lchnk )
     call outfld( 'QCHFTEN', rc_hf_tend, pcols, state%lchnk )
     call outfld( 'QRHFTEN', rr_hf_tend, pcols, state%lchnk )
     call outfld( 'QIHFTEN', ri_hf_tend, pcols, state%lchnk )
     call outfld( 'QSHFTEN', rs_hf_tend, pcols, state%lchnk )
     call outfld( 'THFTEN', s_hf_tend / cpair, pcols, state%lchnk )

     ! Perform total water and total energy conservation checks.
     if ( l_check_conservation ) then

        ! Calculate total water in each grid column.
        ! This calculation is the vertically-integrated grand total water
        ! in each grid column updated for all microphysics and hole filling.
        ! This includes the amount that precipitated to the surface.
        do icol = 1, ncol
           grand_total_water_column_finish(icol) = 0.0_r8
           do k = top_lev, pver
              grand_total_water_column_finish(icol) &
              = grand_total_water_column_finish(icol) &
                + ( state%q(icol,k,1) &
                    + ptend%q(icol,k,1) * dt &
                    + state%q(icol,k,ixcldliq) &
                    + ptend%q(icol,k,ixcldliq) * dt &
                    + state%q(icol,k,ixcldice) &
                    + ptend%q(icol,k,ixcldice) * dt ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 grand_total_water_column_finish(icol) &
                 = grand_total_water_column_finish(icol) &
                   + ( state%q(icol,k,ixrain) + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
              if ( ixsnow > 0 ) then
                 grand_total_water_column_finish(icol) &
                 = grand_total_water_column_finish(icol) &
                   + ( state%q(icol,k,ixsnow) + ptend%q(icol,k,ixsnow) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
           grand_total_water_column_finish(icol) &
           = grand_total_water_column_finish(icol) &
             + prect(icol) * dt * 1000.0_r8
        enddo ! icol = 1, ncol

        ! Calculate total energy in each column.
        ! This calculation is the vertically-integrated total energy in each
        ! grid column updated for all microphysics and hole filling.  This
        ! includes the amount that precipitated to the surface.  Since, the
        ! microphysics code does not directly change kinetic energy,
        ! 0.5 * ( u^2 + v^2 ), it can be skipped as part of the energy
        ! conservation check.
        do icol = 1, ncol
           total_energy_column_finish(icol) = 0.0_r8
           do k = top_lev, pver
              total_energy_column_finish(icol) &
              = total_energy_column_finish(icol) &
                + ( state%s(icol,k) + ptend%s(icol,k) * dt &
                    + ( latvap + latice ) &
                      * ( state%q(icol,k,1) + ptend%q(icol,k,1) * dt ) &
                    + latice * ( state%q(icol,k,ixcldliq) &
                                 + ptend%q(icol,k,ixcldliq) * dt ) ) &
                  * state%pdel(icol,k) / gravit
              if ( ixrain > 0 ) then
                 total_energy_column_finish(icol) &
                 = total_energy_column_finish(icol) &
                   + latice * ( state%q(icol,k,ixrain) &
                                + ptend%q(icol,k,ixrain) * dt ) &
                     * state%pdel(icol,k) / gravit
              endif
           enddo ! k = top_lev, pver
           total_energy_column_finish(icol) &
           = total_energy_column_finish(icol) &
             + latice * precl(icol) * dt * 1000.0_r8
        enddo ! icol = 1, ncol

        ! Calculate the total relative error in each grid column.
        do icol = 1, ncol

           tot_water_rel_err(icol) &
           = abs( ( grand_total_water_column_finish(icol) &
                    - grand_total_water_column_start(icol) ) ) &
             / min( grand_total_water_column_finish(icol), &
                    grand_total_water_column_start(icol) )

           tot_energy_rel_err(icol) &
           = abs( ( total_energy_column_finish(icol) &
                    - total_energy_column_start(icol) ) ) &
             / min( total_energy_column_finish(icol), &
                    total_energy_column_start(icol) )

        enddo ! icol = 1, ncol

        ! Print an error message if any total water relative error is found to
        ! be greater than the threshold.
        if ( any( tot_water_rel_err(:ncol) >= err_thresh ) ) then
           write(iulog,*) "Water conservation error reported in hole filling"
           do icol = 1, ncol
              if ( tot_water_rel_err(icol) >= err_thresh ) then
                 write(iulog,*) "Column = ", icol, &
                          "Relative error = ", tot_water_rel_err(icol), &
                          "Column-integrated grand total water at start = ", &
                          grand_total_water_column_start(icol), &
                          "Column-integrated grand total water at finish = ", &
                          grand_total_water_column_finish(icol)
              endif ! tot_water_rel_err(icol) >= err_thresh
           enddo ! icol = 1, ncol
        endif ! any( tot_water_rel_err >= err_thresh )

        ! Print an error message if any total energy relative error is found to
        ! be greater than the threshold.
        if ( any( tot_energy_rel_err(:ncol) >= err_thresh ) ) then
           write(iulog,*) "Energy conservation error reported in hole filling"
           do icol = 1, ncol
              if ( tot_energy_rel_err(icol) >= err_thresh ) then
                 write(iulog,*) "Column = ", icol, &
                          "Relative error = ", tot_energy_rel_err(icol), &
                          "Column-integrated total energy at start = ", &
                          total_energy_column_start(icol), &
                          "Column-integrated total energy at finish = ", &
                          total_energy_column_finish(icol)
              endif ! tot_energy_rel_err(icol) >= err_thresh
           enddo ! icol = 1, ncol
        endif ! any( tot_energy_rel_err >= err_thresh )

     endif ! l_check_conservation


     return

   end subroutine subcol_SILHS_fill_holes_conserv

   !============================================================================
   subroutine fill_holes_sedimentation( dt, ncol, hm_start, pdel, &
                                        fallspeed_m_per_s, zi, qmin_hm, &
                                        hm_tend, prec )

     ! Description:
     ! After hydrometeor tendencies from microphysics processes were adjusted
     ! so that holes don't form in a hydrometeor field from microphysics
     ! processes, the sedimentation tendency was added back on to produce an
     ! updated total microphysics tendency.  The first-order "up-gravity"
     ! sedimentation method that was originally used is positive definite.
     ! However, since the microphysics process tendencies were altered so that
     ! holes don't form, it is possible that adding the old sedimentation
     ! tendencies back onto the new microphysics process tendencies could
     ! produce new total microphysics tendencies that cause holes to form.
     !
     ! In this subroutine, holes in a hydrometeor field are checked for after
     ! the updated microphysics tendency is applied.  If any are found, they are
     ! filled from positive hydrometeor mass found at grid levels below where
     ! the hole is found.  The levels that are used to fill are within range
     ! based on fallspeed of the hydrometeor.  If the level that contains the
     ! hole is within proximity to the surface, then the water that sedimented
     ! to the surface can be used to fill the hole, as well.
     !
     ! If there isn't enough total hydrometeor mass within the fall range to
     ! fill the hole, then positive hydrometeor mass from levels below the fall
     ! range is to be added to the total available mass to fill the hole.  Mass
     ! is added one level at a time until enough mass is found to fill the hole
     ! or until the surface is reached and the surface precipitation is added to
     ! the total available fill mass.
     !
     ! After this, if there still isn't enough available mass to fill the hole,
     ! then positive hydrometeor mass is added from all levels above the hole to
     ! the total mass that is available to fill the hole.

     !----------------------------------------------------------------------

     use ppgrid, only: &
         pcols

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     real(r8), intent(in) :: dt                   ! Time step duration

     integer, intent(in) :: ncol                  ! Number of grid columns

     real(r8), dimension(pcols,pver), intent(in) :: &
       hm_start, & ! Hydrometeor mixing ratio at start of microphysics  [kg/kg]
       pdel        ! Pressure difference between grid levels            [Pa]

     real(r8), dimension(pcols,pver), intent(in) :: &
       fallspeed_m_per_s    ! Hydrometeor mixing ratio fall speed     [m/s]

     real(r8), dimension(pcols,pverp), intent(in) :: &
       zi    ! Height of momentum (interface) grid levels    [m]

     real(r8), intent(in) :: &
       qmin_hm    ! Minimum threshold value of hydrometeor mixing ratio  [kg/kg]

     ! Input/Output Variables
     real(r8), dimension(pcols,pver), intent(inout) :: &
       hm_tend    ! Hydrometeor mixing ratio tendency  [kg/kg/s]

     real(r8), dimension(pcols), intent(inout) :: &
       prec    ! Precipitation rate (surface)        [m/s]

     ! Local Variables
     real(r8), dimension(pver) :: &
       hm_update, & ! Hydrometeor mixing ratio; start of sed. hole fill  [kg/kg]
       hm_curr      ! Current value of hydrometeor mixing ratio          [kg/kg]

     real(r8) :: &
       total_hole,          & ! Total mass of hole in hydrometeor       [kg/m^2]
       total_fill_mass,     & ! Total mass available to fill hole       [kg/m^2]
       hole_fillmass_ratio, & ! Ratio: total_hole / total_fill_mass     [-]
       fallspeed_Pa_per_s,  & ! Hydrometeor mixing ratio fall speed     [Pa/s]
       total_fall_Pa,       & ! Pressure "distance" hydrometeor fell    [Pa]
       sum_pdel               ! Sum of pdel over levels                 [Pa]

     logical, dimension(pver) :: &
       l_pos_hm  ! Flag for a hydrometeor having a positive (>= qmin_hm) value

     ! Flag for whether surface precipitation mass needs to be included in
     ! the total_fill_mass for hole filling.
     logical :: l_reached_surface

     ! Flag for whether hydrometeor mass from levels above the hole needs to be
     ! included in the total_fill_mass for hole filling.
     logical :: l_fill_from_above

     integer :: icol  ! Grid column index

     integer :: k, idx  ! Vertical grid level indices

     ! Index of the lowest vertical grid level that needs to be included in the
     ! total_fill_mass for hole filling.
     integer :: lowest_level_idx


     ! Loop over all columns, performing any adjustments one column at a time.
     do icol = 1, ncol

        ! Calculate the updated value of the hydrometeor field based on the
        ! updated microphysics tendency.  Since the original sedimentation
        ! tendency has been added to the updated microphysics process tendency
        ! to produce the updated total microphysics tendency (hm_tend), the
        ! updated value of the hydrometeor field (hm_update) could be negative.
        hm_update = hm_start(icol,:) + hm_tend(icol,:) * dt
        hm_curr = hm_update

        ! Check for any holes in the vertical profile
        if ( any( hm_curr(top_lev:pver) < qmin_hm ) ) then

           ! At least one hole is found in this hydrometeor species in this
           ! grid column.  The holes must be filled conservatively.

           ! Check which levels have values of the hydrometeor that are at or
           ! above the minimum threshold value.
           do k = top_lev, pver
              if ( hm_curr(k) >= qmin_hm ) then
                l_pos_hm(k) = .true.
              else ! hm_curr < qmin_hm
                l_pos_hm(k) = .false.
              endif ! hm_curr >= qmin_hm
           enddo ! k = top_lev, pver

           do k = pver, top_lev, -1

              if ( .not. l_pos_hm(k) ) then

                 ! A hole is found in the hydrometeor at this grid level.

                 ! Calculate the total hydrometeor mass of the hole that needs
                 ! to be filled.
                 ! The value of the hydrometeor mixing ratio is negative, but
                 ! the value of total_hole is positive.
                 total_hole = ( qmin_hm - hm_curr(k) ) * pdel(icol,k) / gravit

                 ! Calculate the total hydrometeor mass available from below
                 ! to fill the hole.
                 if ( k == pver ) then

                    ! A hole is found at the lowermost level.
                    ! The only place the hydrometeor could have sedimented
                    ! to is the surface, so fill from only the surface.
                    l_reached_surface = .true.

                    ! Calculate the available amount of hydrometeor mass to
                    ! fill the hole.
                    total_fill_mass = prec(icol) * dt * 1000.0_r8

                 else ! top_lev <= k < pver

                    ! Calculate the hydrometeor fallspeed in Pa/s.
                    ! In MG2, the equation for this is given by:
                    !
                    ! fallspeed([Pa/s]) = g * rho * fallspeed([m/s]).
                    !
                    ! The value of rho is typically calculated from the
                    ! hydrostatic approximation:
                    !
                    ! rho = - ( 1 / g ) * dp/dz.
                    !
                    ! The equation for fallspeed in Pa/s becomes:
                    !
                    ! fallspeed([Pa/s]) = - dp/dz * fallspeed([m/s]).
                    fallspeed_Pa_per_s &
                    = fallspeed_m_per_s(icol,k) &
                      * pdel(icol,k) / ( zi(icol,k) - zi(icol,k+1) )

                    ! Calculate the fall "distance" in Pa.
                    total_fall_Pa = fallspeed_Pa_per_s * dt

                    ! Find the index of the vertical level that the hydrometeor
                    ! sedimented to in one timestep.  It must sediment at least
                    ! one level.
                    sum_pdel = 0.0_r8
                    idx = k + 1
                    do
                       ! Update the total pressure difference between the
                       ! level of origin and the current level.
                       sum_pdel = sum_pdel + pdel(icol,idx)
                       if ( sum_pdel >= total_fall_Pa ) then
                          ! The total pressure difference between the level of
                          ! origin and the current level exceeds the total
                          ! hydrometeor fall "distance" (in Pa).
                          lowest_level_idx = idx
                          l_reached_surface = .false.
                          exit
                       else ! sum_pdel < total_fall_Pa
                          ! The total hydrometeor fall "distance" (in Pa)
                          ! exceeds the total pressure difference between the
                          ! level of origin and the current level.
                          if ( idx == pver ) then
                             ! The lowest level of the model has been reached.
                             ! The hydrometeor sedimented to the surface.
                             lowest_level_idx = pver
                             l_reached_surface = .true.
                             exit
                          else ! idx < pver
                             ! Increment idx and keep going.
                             idx = idx + 1
                          endif ! idx == pver
                       endif ! sum_pdel >= total_fall_Pa
                    enddo

                    ! Calculate the available amount of hydrometeor mass to
                    ! fill the hole.
                    total_fill_mass = 0.0_r8
                    if ( l_reached_surface ) then
                       ! The hydrometeor sedimented to the surface, so
                       ! automatically loop down to pver and include the
                       ! surface mass.
                       do idx = k+1, pver, 1
                          if ( l_pos_hm(idx) ) then
                             total_fill_mass &
                             = total_fill_mass &
                               + ( hm_curr(idx) - qmin_hm ) &
                                 * pdel(icol,idx) / gravit
                          endif ! l_pos_hm(idx)
                       enddo ! idx = k+1, pver, 1
                       ! Contribution to total fill mass from the surface.
                       total_fill_mass &
                       = total_fill_mass + prec(icol) * dt * 1000.0_r8
                    else ! .not. l_reached_surface
                       ! The hydrometeor sedimented to lowest_level_idx.
                       idx = k + 1
                       do
                          if ( l_pos_hm(idx) ) then
                             total_fill_mass &
                             = total_fill_mass &
                               + ( hm_curr(idx) - qmin_hm ) &
                                 * pdel(icol,idx) / gravit
                          endif ! l_pos_hm(idx)
                          if ( idx >= lowest_level_idx ) then
                             ! Check if enough mass has been gathered in
                             ! total_fill_mass to fill the hole.
                             if ( total_fill_mass >= total_hole ) then
                                ! There has been enough total_fill_mass
                                ! gathered to completely fill the hole.
                                lowest_level_idx = idx
                                exit
                             else ! total_fill_mass < total_hole
                                ! Even though lowest_level_idx has been reached,
                                ! more total_fill_mass needs to be added in
                                ! order to completely fill the hole, so keep
                                ! going.
                                if ( idx == pver ) then
                                   ! The lowest vertical level has already been
                                   ! reached, so go to the surface.
                                   lowest_level_idx = pver
                                   l_reached_surface = .true.
                                   ! Contribution to total fill mass from the
                                   ! surface.
                                   total_fill_mass &
                                   = total_fill_mass &
                                     + prec(icol) * dt * 1000.0_r8
                                   exit
                                else ! idx < pver
                                   ! Haven't reached pver yet, so increment
                                   ! and keep going.
                                   idx = idx + 1
                                endif ! idx == pver
                             endif ! total_fill_mass >= total_hole
                          else ! idx < lowest_level_idx
                             ! Haven't reached lowest_level_idx yet, so
                             ! increment and keep going.
                             idx = idx + 1
                          endif ! idx >= lowest_level_idx
                       enddo
                    endif ! l_reached_surface

                 endif ! k == pver

                 ! If mass has been added all the way down to the surface and
                 ! there's still not enough mass to fill the hole, then fill the
                 ! hole pulling mass from above.
                 if ( total_fill_mass >= total_hole ) then
                    l_fill_from_above = .false.
                 else ! total_fill_mass < total_hole
                    l_fill_from_above = .true.
                    do idx = top_lev, k-1, 1
                       if ( l_pos_hm(idx) ) then
                          total_fill_mass &
                          = total_fill_mass &
                            + ( hm_curr(idx) - qmin_hm ) &
                              * pdel(icol,idx) / gravit
                       endif ! l_pos_hm(idx)
                    enddo ! idx = top_lev, k-1, 1
                 endif ! total_fill_mass >= total_hole

                 ! Calculate the ratio of total hole to total fill mass.  This
                 ! should not exceed 1 except as a result of numerical round-off
                 ! errors.  Use thresholding to be safe.
                 hole_fillmass_ratio &
                 = min( total_hole / max( total_fill_mass, 1.0e-30_r8 ), &
                        1.0_r8 )

                 if ( k < pver ) then
                    ! Modify (reduce) the amount of the hydrometeor at levels
                    ! that were used to fill the hole.
                    do idx = k+1, lowest_level_idx
                       if ( l_pos_hm(idx) ) then
                          ! Since pdel at a grid level does not change and
                          ! gravit is constant, the only variable that needs to
                          ! be modified proportionately is hm_curr.
                          hm_curr(idx) &
                          = qmin_hm &
                            + ( hm_curr(idx) - qmin_hm ) &
                              * ( 1.0_r8 - hole_fillmass_ratio )
                       endif ! l_pos_hm(idx)
                    enddo ! idx = k+1, lowest_level_idx
                 endif ! k < pver

                 if ( l_reached_surface ) then
                    ! Modify (reduce) the amount of surface precipitation in
                    ! order to fill the hole.  Since dt and 1000 are constants,
                    ! the only variable that needs to be modified
                    ! proportionately is prec.
                    prec(icol) = prec(icol) * ( 1.0_r8 - hole_fillmass_ratio )
                 endif ! l_reached_surface

                 if ( l_fill_from_above ) then
                    ! Modify (reduce) the amount of the hydrometeor at levels
                    ! that were used to fill the hole.
                    do idx = top_lev, k-1
                       if ( l_pos_hm(idx) ) then
                          ! Since pdel at a grid level does not change and
                          ! gravit is constant, the only variable that needs to
                          ! be modified proportionately is hm_curr.
                          hm_curr(idx) &
                          = qmin_hm &
                            + ( hm_curr(idx) - qmin_hm ) &
                              * ( 1.0_r8 - hole_fillmass_ratio )
                       endif ! l_pos_hm(idx)
                    enddo ! idx = top_lev, k-1
                 endif ! l_fill_from_above

                 ! Update the value of the hydrometeor at the level where the
                 ! hole was found.  Mathematically, as long as the available
                 ! mass was able to fill the entire hole, the new value of the
                 ! hydrometeor mixing ratio (hm_curr) should be qmin_hm.
                 hm_curr(k) &
                 = hm_curr(k) &
                   + hole_fillmass_ratio * total_fill_mass &
                     * gravit / pdel(icol,k)

              endif ! .not. l_pos_hm(k)

           enddo ! k = pver, top_lev, -1

        endif ! any( hm_curr(top_lev:pver) < qmin_hm )

        ! Update the value of total microphysics tendency after hole filling.
        hm_tend(icol,:) = hm_tend(icol,:) + ( hm_curr - hm_update ) / dt

     enddo ! icol = 1, ncol


     return

   end subroutine fill_holes_sedimentation

   !============================================================================
   subroutine fill_holes_same_phase_vert( dt, ncol, hm_start, hm_start_filler, &
                                          pdel, qmin_hm, qmin_hm_filler, &
                                          hm_tend, hm_tend_filler )

     ! Description:
     ! Fills remaining holes in a hydrometeor with mass from the the vertical
     ! profile of another hydrometeor of the same phase.  Remaining holes in
     ! cloud water are filled with rain water and remaining holes in snow are
     ! filled with cloud ice.
     !
     ! This subroutine, combined with subroutine fill_holes_sedimentation, fill
     ! holes making use of the following principles:
     !
     ! 3) When adjusting the hydrometeor tendency from sedimentation of a
     !    liquid hydrometeor (cloud water or rain water), conserve:
     !
     !    SUM(k=top_lev:pver) ( rc_sed_tend(k) + rr_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + precl * dt * 1000 = 0.
     !
     ! 4) When adjusting the hydrometeor tendency from sedimentation of a
     !    frozen hydrometeor (cloud ice or snow), conserve:
     !
     !    SUM(k=top_lev:pver) ( ri_sed_tend(k) + rs_sed_tend(k) )
     !                        * dt * pdel(k) / g
     !    + preci * dt * 1000 = 0.
     !
     ! These two equations (one for liquid-phase hydrometeors and one for
     ! ice-phase hydrometeors) could be further split into one equation for
     ! each hydrometeor if there was prec output for each hydrometeor.  However,
     ! there's only prec output for ice-phase precipitation rate and total
     ! precipitation rate (liquid preciptation rate is total rate minus
     ! ice-phase rate).
     !
     ! Since only liquid-phase precipitation rate (precl) and ice-phase
     ! precipitation rate (preci) are available, and there are two hydrometeors
     ! in each category, one hydrometeor from each category must fill before
     ! the other hydrometeor from its category and get priority access to precl
     ! or preci.  Since a vast majority of liquid precipitation comes from rain
     ! rather than sedimenting cloud water, rain is filled before cloud water
     ! and gets priority access to precl.  Likewise, since a vast majority of
     ! frozen precipitation comes from snow rather than sedimenting cloud ice,
     ! snow is filled before cloud ice and gets priority access to preci.
     !
     ! The order of sedimentation hole filling is as follows.  First, a level
     ! with a hole in it is identified.  The fall distance for the hydrometeor
     ! that originated at a level is calculated.  Total mass to fill the hole is
     ! calculated from all levels within the fall range that have positive
     ! values of the hydrometeor.  The amount that precipitated to the surface
     ! is also included if the hydrometeor fell that far.  If that isn't enough
     ! mass to fill the hole, then levels that are lower in the profile are
     ! included (if the hydrometeor has a positive value) until enough mass is
     ! found to fill the hole or until the surface is reached.  If there isn't
     ! enough mass found in all levels below the hole, including the amount that
     ! precipitated to the ground, to fill the hole, then the hydrometeor mass
     ! from all levels above the hole (again, where a positive value of the
     ! hydrometeor is found) are included in the total available mass to fill
     ! the hole.
     !
     ! Occasionally, a situation can occur where both hydrometeors in a category
     ! contributed to surface precipitation rate, and filling a hole in rain
     ! (or snow) can deplete all the surface precl (or preci), resulting in not
     ! enough water mass in the vertical profile (including the surface) of
     ! cloud water (or cloud ice) to fill a hole in cloud water (or cloud ice).
     ! When this happens, there must still be liquid water (or frozen water)
     ! found in the vertical profile of rain (or snow), so pull the water from
     ! rain (or snow) to fill any remaining holes in cloud water (or cloud ice).

     !----------------------------------------------------------------------

     use ppgrid, only: &
         pcols

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     real(r8), intent(in) :: dt                   ! Time step duration

     integer, intent(in) :: ncol                  ! Number of grid columns

     real(r8), dimension(pcols,pver), intent(in) :: &
       hm_start,        & ! Hydrometeor mixing ratio (microphys start)   [kg/kg]
       hm_start_filler, & ! Filler hydromet mix ratio (microphys start)  [kg/kg]
       pdel               ! Pressure difference between grid levels      [Pa]

     real(r8), intent(in) :: &
       qmin_hm,        & ! Minimum threshold hydrometeor mixing ratio  [kg/kg]
       qmin_hm_filler    ! Min threshold filler hydromet mixing ratio  [kg/kg]

     ! Input/Output Variables
     real(r8), dimension(pcols,pver), intent(inout) :: &
       hm_tend,        & ! Hydrometeor mixing ratio tendency         [kg/kg/s]
       hm_tend_filler    ! Filler hydrometeor mixing ratio tendency  [kg/kg/s]

     ! Local Variables
     real(r8), dimension(pver) :: &
       hm_update,        & ! Hydrometeor mixing ratio; start           [kg/kg]
       hm_update_filler, & ! Filler Hydrometeor mixing ratio; start    [kg/kg]
       hm_curr,          & ! Current hydrometeor mixing ratio          [kg/kg]
       hm_curr_filler      ! Current filler hydrometeor mixing ratio   [kg/kg]

     real(r8) :: &
       total_hole,          & ! Total mass of hole in hydrometeor       [kg/m^2]
       total_fill_mass,     & ! Total mass available to fill hole       [kg/m^2]
       hole_fillmass_ratio    ! Ratio: total_hole / total_fill_mass     [-]

     logical, dimension(pver) :: &
       l_pos_hm,        & ! Flag: hydrometeor has positive (>= qmin_hm) value
       l_pos_hm_filler    ! Flag: filler hydrometeor has positive value

     integer :: icol  ! Grid column index

     integer :: k, idx  ! Vertical grid level indices


     ! Loop over all columns, performing any adjustments one column at a time.
     do icol = 1, ncol

        ! Calculate the updated value of the hydrometeor field based on the
        ! updated microphysics tendency.
        hm_update = hm_start(icol,:) + hm_tend(icol,:) * dt
        hm_curr = hm_update

        ! Calculate the updated value of the filler hydrometeor field based on
        ! the updated microphysics tendency.
        hm_update_filler = hm_start_filler(icol,:) + hm_tend_filler(icol,:) * dt
        hm_curr_filler = hm_update_filler

        ! Check for any holes in the vertical profile
        if ( any( hm_curr(top_lev:pver) < qmin_hm ) ) then

           ! At least one hole is found in this hydrometeor species in this
           ! grid column.  The holes must be filled conservatively.

           ! Check which levels have values of the hydrometeor that are at or
           ! above the minimum threshold value.
           do k = top_lev, pver
              ! Check for the hydrometeor that might need to be filled.
              if ( hm_curr(k) >= qmin_hm ) then
                l_pos_hm(k) = .true.
              else ! hm_curr < qmin_hm
                l_pos_hm(k) = .false.
              endif ! hm_curr >= qmin_hm
              ! Check for the filler hydrometeor, as some levels might have
              ! numerical round-off level, small negative values.
              if ( hm_curr_filler(k) >= qmin_hm_filler ) then
                l_pos_hm_filler(k) = .true.
              else ! hm_curr_filler < qmin_hm_filler
                l_pos_hm_filler(k) = .false.
              endif ! hm_curr_filler >= qmin_hm_filler
           enddo ! k = top_lev, pver

           do k = top_lev, pver

              if ( .not. l_pos_hm(k) ) then

                 ! A hole is found in the hydrometeor at this grid level.

                 ! Calculate the total hydrometeor mass of the hole that needs
                 ! to be filled.
                 ! The value of the hydrometeor mixing ratio is negative, but
                 ! the value of total_hole is positive.
                 total_hole = ( qmin_hm - hm_curr(k) ) * pdel(icol,k) / gravit

                 ! Calculate the total hydrometeor mass available from the
                 ! filler hydrometeor to fill the hole.
                 total_fill_mass = 0.0_r8
                 do idx = top_lev, pver, 1
                    if ( l_pos_hm_filler(idx) ) then
                        total_fill_mass &
                        = total_fill_mass &
                          + ( hm_curr_filler(idx) - qmin_hm_filler ) &
                            * pdel(icol,idx) / gravit
                     endif ! l_pos_hm_filler(idx)
                 enddo ! idx = top_lev, pver, 1

                 ! Calculate the ratio of total hole to total fill mass.  This
                 ! should not exceed 1 except as a result of numerical round-off
                 ! errors.  Use thresholding to be safe.
                 hole_fillmass_ratio &
                 = min( total_hole / max( total_fill_mass, 1.0e-30_r8 ), &
                        1.0_r8 )

                 ! Modify (reduce) the amount of the filler hydrometeor.
                 do idx = top_lev, pver
                    if ( l_pos_hm_filler(idx) ) then
                       ! Since pdel at a grid level does not change and gravit
                       ! is constant, the only variable that needs to be
                       ! modified proportionately is hm_curr_filler.
                       hm_curr_filler(idx) &
                       = qmin_hm_filler &
                         + ( hm_curr_filler(idx) - qmin_hm_filler ) &
                           * ( 1.0_r8 - hole_fillmass_ratio )
                    endif ! l_pos_hm_filler(idx)
                 enddo ! idx = top_lev, pver

                 ! Update the value of the hydrometeor at the level where the
                 ! hole was found.  Mathematically, as long as the available
                 ! mass was able to fill the entire hole, the new value of the
                 ! hydrometeor mixing ratio (hm_curr) should be qmin_hm.
                 hm_curr(k) &
                 = hm_curr(k) &
                   + hole_fillmass_ratio * total_fill_mass &
                     * gravit / pdel(icol,k)

              endif ! .not. l_pos_hm(k)

           enddo ! k = top_lev, pver

        endif ! any( hm_curr(top_lev:pver) < qmin_hm )

        ! Update the value of total microphysics tendency after hole filling.
        hm_tend(icol,:) = hm_tend(icol,:) + ( hm_curr - hm_update ) / dt

        ! Update the value of total microphysics tendency after hole filling for
        ! the filler hydrometeor.
        hm_tend_filler(icol,:) &
        = hm_tend_filler(icol,:) + ( hm_curr_filler - hm_update_filler ) / dt

     enddo ! icol = 1, ncol


     return

   end subroutine fill_holes_same_phase_vert

   !============================================================================
   subroutine subcol_SILHS_hydromet_conc_tend_lim( state, dt, ptend )

     ! Description:
     ! Limits the values of mean hydrometeor concentrations so that the mean
     ! drop size for the hydrometeor type remains reasonable and does not become
     ! too large.

     !----------------------------------------------------------------------

     use shr_const_mod, only: &
         shr_const_pi, &
         shr_const_rhofw

     use constituents, only: &
         qmin

     use ref_pres, only: &
         top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     type(physics_state), intent(in) :: state     ! Physics state variables
     real(r8), intent(in) :: dt                   ! Time step duration

     ! Input/Output Variable
     type(physics_ptend),  intent(inout) :: ptend  ! Parameterization tendencies

     ! Local Variables
     real( r8 ) :: &
       rcm_update, & ! New value of mean cloud water mixing ratio    [kg/kg]
       rrm_update, & ! New value of mean rain water mixing ratio     [kg/kg]
       rim_update, & ! New value of mean ice mixing ratio            [kg/kg]
       rsm_update    ! New value of mean snow mixing ratio           [kg/kg]

     real( r8 ) :: &
       Nc_tend_min, & ! Minimum value of cloud droplet conc. tendency [num/kg/s]
       Nr_tend_min, & ! Minimum value of rain drop conc. tendency     [num/kg/s]
       Ni_tend_min, & ! Minimum value of ice conc. tendency           [num/kg/s]
       Ns_tend_min    ! Minimum value of snow conc. tendency          [num/kg/s]

     real( r8 ), parameter :: &
       four_thirds = 4.0_r8/3.0_r8,    & ! 4/3
       rho_ice     = 917.0_r8,         & ! Density of ice            [kg/m^3]
       rho_lw      = shr_const_rhofw,  & ! Density of liquid water   [kg/m^3]
       pi          = shr_const_pi        ! Pi

     real( r8 ), parameter :: &
       mvr_cloud_max = 1.6E-5_r8, & ! Max. avg. mean vol. rad. cloud   [m]
       mvr_rain_max  = 5.0E-3_r8, & ! Max. avg. mean vol. rad. rain    [m]
       mvr_ice_max   = 1.3E-4_r8, & ! Max. avg. mean vol. rad. ice     [m]
       mvr_snow_max  = 1.0E-2_r8    ! Max. avg. mean vol. rad. snow    [m]

     ! Calculate the coefficient for the minimum mean cloud droplet
     ! concentration, where <Nc>|_min = Ncm_min_coef * <rc> and has units of
     ! 1/kg.
     real( r8 ), parameter :: &
       Ncm_min_coef = 1.0_r8 / ( four_thirds * pi * rho_lw * mvr_cloud_max**3 )

     ! Calculate the coefficient for the minimum mean rain drop concentration,
     ! where <Nr>|_min = Nrm_min_coef * <rr> and has units of 1/kg.
     real( r8 ), parameter :: &
       Nrm_min_coef = 1.0_r8 / ( four_thirds * pi * rho_lw * mvr_rain_max**3 )

     ! Calculate the coefficient for the minimum mean ice crystal concentration,
     ! where <Ni>|_min = Nim_min_coef * <ri> and has units of 1/kg.
     real( r8 ), parameter :: &
       Nim_min_coef = 1.0_r8 / ( four_thirds * pi * rho_ice * mvr_ice_max**3 )

     ! Calculate the coefficient for the minimum mean snow flake concentration,
     ! where <Ns>|_min = Nsm_min_coef * <rs> and has units of 1/kg.
     real( r8 ), parameter :: &
       Nsm_min_coef = 1.0_r8 / ( four_thirds * pi * rho_ice * mvr_snow_max**3 )

     integer :: ncol  ! Number of grid columns

     integer :: icol    ! Column loop index

     integer :: k    ! Vertical level loop index


     ! Get the number of grid columns.
     ncol = state%ncol

     ! Loop over all grid columns.
     do icol = 1, ncol

        ! Loop over all vertical levels from top_lev to pver.
        do k = top_lev, pver

           ! Cloud droplet concentration
           if ( ixcldliq > 0 .and. ixnumliq > 0 ) then

              ! Calculate the value of cloud water mixing ratio after the
              ! update.
              rcm_update &
              = max( state%q(icol,k,ixcldliq) + ptend%q(icol,k,ixcldliq) * dt, &
                     qmin(ixcldliq) )

              ! Calculate the limiting cloud droplet concentration tendency so
              ! that cloud maintains a reasonable (not too big) mean volume
              ! radius.
              Nc_tend_min &
              = ( Ncm_min_coef * rcm_update - state%q(icol,k,ixnumliq) ) / dt

              ! The cloud droplet concentration tendency needs to be the greater
              ! of the current Nc_tend and Nc_tend_min.
              ptend%q(icol,k,ixnumliq) &
              = max( ptend%q(icol,k,ixnumliq), Nc_tend_min )

           endif ! ixcldliq > 0 .and. ixnumliq > 0

           ! Rain drop concentration
           if ( ixrain > 0 .and. ixnumrain > 0 ) then

              ! Calculate the value of rain water mixing ratio after the update.
              rrm_update &
              = max( state%q(icol,k,ixrain) + ptend%q(icol,k,ixrain) * dt, &
                     qmin(ixrain) )

              ! Calculate the limiting rain drop concentration tendency so that
              ! rain maintains a reasonable (not too big) mean volume radius.
              Nr_tend_min &
              = ( Nrm_min_coef * rrm_update - state%q(icol,k,ixnumrain) ) / dt

              ! The rain drop concentration tendency needs to be the greater of
              ! the current Nr_tend and Nr_tend_min.
              ptend%q(icol,k,ixnumrain) &
              = max( ptend%q(icol,k,ixnumrain), Nr_tend_min )

           endif ! ixrain > 0 .and. ixnumrain > 0

           ! Ice crystal concentration
           if ( ixcldice > 0 .and. ixnumice > 0 ) then

              ! Calculate the value of ice mixing ratio after the update.
              rim_update &
              = max( state%q(icol,k,ixcldice) + ptend%q(icol,k,ixcldice) * dt, &
                     qmin(ixcldice) )

              ! Calculate the limiting ice crystal concentration tendency so
              ! that ice maintains a reasonable (not too big) mean volume
              ! radius.
              Ni_tend_min &
              = ( Nim_min_coef * rim_update - state%q(icol,k,ixnumice) ) / dt

              ! The ice crystal concentration tendency needs to be the greater
              ! of the current Ni_tend and Ni_tend_min.
              ptend%q(icol,k,ixnumice) &
              = max( ptend%q(icol,k,ixnumice), Ni_tend_min )

           endif ! ixcldice > 0 .and. ixnumice > 0

           ! Snow flake concentration
           if ( ixsnow > 0 .and. ixnumsnow > 0 ) then

              ! Calculate the value of snow mixing ratio after the update.
              rsm_update &
              = max( state%q(icol,k,ixsnow) + ptend%q(icol,k,ixsnow) * dt, &
                     qmin(ixsnow) )

              ! Calculate the limiting snow flake concentration tendency so that
              ! snow maintains a reasonable (not too big) mean volume radius.
              Ns_tend_min &
              = ( Nsm_min_coef * rsm_update - state%q(icol,k,ixnumsnow) ) / dt

              ! The snow flake concentration tendency needs to be the greater of
              ! the current Ns_tend and Ns_tend_min.
              ptend%q(icol,k,ixnumsnow) &
              = max( ptend%q(icol,k,ixnumsnow), Ns_tend_min )

           endif ! ixsnow > 0 .and. ixnumsnow > 0

        enddo ! k = top_lev, pver

     enddo ! icol = 1, ncol


     return

   end subroutine subcol_SILHS_hydromet_conc_tend_lim

   !============================================================================
   
end module subcol_SILHS 
