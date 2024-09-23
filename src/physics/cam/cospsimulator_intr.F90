module cospsimulator_intr
  ! ######################################################################################
  ! Purpose: CAM interface to
  !         Name:         CFMIP Observational Simulator Package Version 2 (COSP2)
  !         What:         Simulate ISCCP/CloudSat/CALIPSO/MISR/MODIS cloud products from 
  !                       GCM inputs
  !         Version:      v2.1.4 (August 2019)
  !         Authors:      Dustin Swales (dustin.swales@noaa.gov)
  !
  ! Modifications:
  !
  ! ######################################################################################
  use shr_kind_mod,         only: r8 => shr_kind_r8
  use spmd_utils,           only: masterproc
  use ppgrid,               only: pcols, pver, pverp, begchunk, endchunk
  use ref_pres,             only: ktop => trop_cloud_top_lev
  use perf_mod,             only: t_startf, t_stopf
  use cam_abortutils,       only: endrun, handle_allocate_error
  use phys_control,         only: cam_physpkg_is
  use cam_logfile,          only: iulog
#ifdef USE_COSP
  use quickbeam,            only: radar_cfg
  use mod_quickbeam_optics, only: size_distribution
  use mod_cosp,             only: cosp_outputs, cosp_optical_inputs, cosp_column_inputs
  use mod_cosp_config,      only: pres_binCenters, pres_binEdges, tau_binCenters,      &
       tau_binEdges, cloudsat_binCenters, cloudsat_binEdges, calipso_binCenters,       &
       calipso_binEdges, misr_histHgtCenters, misr_histHgtEdges,  PARASOL_SZA,         &
       R_UNDEF, PARASOL_NREFL, LIDAR_NCAT,SR_BINS, N_HYDRO, RTTOV_MAX_CHANNELS,        &
       numMISRHgtBins, CLOUDSAT_DBZE_BINS, LIDAR_NTEMP, calipso_histBsct,              &
       numMODISTauBins, numMODISPresBins, numMODISReffIceBins, numMODISReffLiqBins,    &
       numISCCPTauBins, numISCCPPresBins, numMISRTauBins, reffICE_binEdges,            &
       reffICE_binCenters, reffLIQ_binEdges, reffLIQ_binCenters, LIDAR_NTYPE,          &
       nCloudsatPrecipClass, &
       nsza_cosp         => PARASOL_NREFL,       &
       nprs_cosp         => npres,               &
       ntau_cosp         => ntau,                &
       ntau_cosp_modis   => ntau,                &
       nsr_cosp          => SR_BINS,             &
       nhtmisr_cosp      => numMISRHgtBins,      &
       nhydro            => N_HYDRO, &
       cloudsat_preclvl
    use mod_cosp_stats,       only: cosp_change_vertical_grid
#endif
  implicit none
  private
  save
   
  ! Public functions/subroutines
  public :: &
       cospsimulator_intr_readnl,  &
       cospsimulator_intr_register,&
       cospsimulator_intr_init,    &
       cospsimulator_intr_run

  ! ######################################################################################
  ! Public declarations
  ! ######################################################################################
  ! Whether to do COSP calcs and I/O, default is false. If docosp is specified in 
  ! the atm_in namelist, this value is overwritten and cosp is run
  logical, public, protected :: docosp = .false.  

  ! Frequency at which cosp is called, every cosp_nradsteps radiation timestep
  integer, public, protected :: cosp_nradsteps = 1
  
#ifdef USE_COSP

  ! ######################################################################################  
  ! Local declarations
  ! ######################################################################################
  integer ::  &
       nlay,        &     ! Number of CAM layers used by COSP.
       nlayp,       &     ! Number of CAM layer interfaces used by COSP.
       nscol_cosp,  &     ! Number of subcolumns, allow namelist input to set.
       nht_cosp           ! Number of height for COSP radar and calipso simulator outputs.  
                          !  *set to 40 if csat_vgrid=.true., else set to Nlr*

  
  ! ######################################################################################
  ! Bin-boundaries for mixed dimensions. Calculated in cospsetupvales OR in cosp_config.F90
  ! ######################################################################################
  real(r8), target :: prsmid_cosp(nprs_cosp)            ! pressure midpoints of COSP ISCCP output
  real(r8), target :: prslim_cosp(2,nprs_cosp)
  real(r8), target :: taumid_cosp(ntau_cosp)            ! optical depth midpoints of COSP ISCCP output
  real(r8), target :: taulim_cosp(2,ntau_cosp)
  real(r8), target :: srmid_cosp(nsr_cosp)              ! sr midpoints of COSP lidar output     
  real(r8), target :: srlim_cosp(2,nsr_cosp)
  real(r8), target :: sza_cosp(nsza_cosp)
  real(r8), target :: dbzemid_cosp(CLOUDSAT_DBZE_BINS)          ! dbze midpoints of COSP radar output
  real(r8), target :: dbzelim_cosp(2,CLOUDSAT_DBZE_BINS)
  real(r8), target :: htmisrmid_cosp(nhtmisr_cosp)      ! htmisr midpoints of COSP misr simulator output
  real(r8), target :: htmisrlim_cosp(2,nhtmisr_cosp)
  real(r8), target :: taumid_cosp_modis(ntau_cosp_modis)! optical depth midpoints of COSP MODIS output
  real(r8), target :: taulim_cosp_modis(2,ntau_cosp_modis)
  real(r8), target :: reffICE_binEdges_cosp(2,numMODISReffIceBins)
  real(r8), target :: reffLIQ_binEdges_cosp(2,numMODISReffLiqBins)
  real(r8), target :: reffICE_binCenters_cosp(numMODISReffIceBins)
  real(r8), target :: reffLIQ_binCenters_cosp(numMODISReffLiqBins)

  integer  :: prstau_cosp(nprs_cosp*ntau_cosp)             ! ISCCP mixed output dimension index
  integer  :: prstau_cosp_modis(nprs_cosp*ntau_cosp_modis) ! MODIS mixed output dimension index
  integer  :: htmisrtau_cosp(nhtmisr_cosp*ntau_cosp)       ! MISR mixed output dimension index
  real(r8) :: prstau_prsmid_cosp(nprs_cosp*ntau_cosp)
  real(r8) :: prstau_taumid_cosp(nprs_cosp*ntau_cosp)
  real(r8) :: prstau_prsmid_cosp_modis(nprs_cosp*ntau_cosp_modis)
  real(r8) :: prstau_taumid_cosp_modis(nprs_cosp*ntau_cosp_modis)
  real(r8) :: htmisrtau_htmisrmid_cosp(nhtmisr_cosp*ntau_cosp)
  real(r8) :: htmisrtau_taumid_cosp(nhtmisr_cosp*ntau_cosp)
  real(r8),allocatable         :: htmlmid_cosp(:)          ! Model level height midpoints for output (nlay)
  real(r8),allocatable, public :: htdbze_dbzemid_cosp(:)   ! (nht_cosp*CLOUDSAT_DBZE_BINS)
  real(r8),allocatable, target :: htlim_cosp(:,:)          ! height limits for COSP outputs (nht_cosp+1)
  real(r8),allocatable, target :: htmid_cosp(:)            ! height midpoints of COSP radar/lidar output (nht_cosp)
  real(r8),allocatable         :: htlim_cosp_1d(:)         ! height limits for COSP outputs (nht_cosp+1)
  real(r8),allocatable         :: htdbze_htmid_cosp(:)     ! (nht_cosp*CLOUDSAT_DBZE_BINS)
  real(r8),allocatable         :: htsr_htmid_cosp(:)       ! (nht_cosp*nsr_cosp)
  real(r8),allocatable         :: htsr_srmid_cosp(:)       ! (nht_cosp*nsr_cosp)
  real(r8),allocatable         :: htmlscol_htmlmid_cosp(:) ! (nlay*nscol_cosp)
  real(r8),allocatable         :: htmlscol_scol_cosp(:)    ! (nlay*nscol_cosp) 
  integer, allocatable, target :: scol_cosp(:)             ! sub-column number (nscol_cosp)
  integer, allocatable         :: htdbze_cosp(:)           ! radar CFAD mixed output dimension index (nht_cosp*CLOUDSAT_DBZE_BINS)
  integer, allocatable         :: htsr_cosp(:)             ! lidar CFAD mixed output dimension index (nht_cosp*nsr_cosp)
  integer, allocatable         :: htmlscol_cosp(:)         ! html-subcolumn mixed output dimension index (nlay*nscol_cosp)

  ! ######################################################################################
  ! Default CAM namelist settings
  ! ######################################################################################
  logical :: cosp_amwg             = .false.
  logical :: cosp_lite             = .false.
  logical :: cosp_passive          = .false.
  logical :: cosp_active           = .false.
  logical :: cosp_isccp            = .false.
  logical :: cosp_lradar_sim       = .false.
  logical :: cosp_llidar_sim       = .false.
  logical :: cosp_lisccp_sim       = .false.
  logical :: cosp_lmisr_sim        = .false.
  logical :: cosp_lmodis_sim       = .false.
  logical :: cosp_histfile_aux     = .false.
  logical :: cosp_lfrac_out        = .false.
  logical :: cosp_runall           = .false.
  integer :: cosp_ncolumns         = 50
  integer :: cosp_histfile_num     =  1
  integer :: cosp_histfile_aux_num = -1
  
  ! COSP
  logical :: lradar_sim       = .false.
  logical :: llidar_sim       = .false.
  logical :: lparasol_sim     = .false.
  logical :: lgrLidar532      = .false.
  logical :: latlid           = .false.
  logical :: lisccp_sim       = .false.
  logical :: lmisr_sim        = .false.
  logical :: lmodis_sim       = .false.
  logical :: lrttov_sim       = .false.
  logical :: lfrac_out        = .false.

  ! ######################################################################################  
  ! COSP parameters
  ! ######################################################################################
  integer, parameter :: Npoints_it = 10000       ! Max # gridpoints to be processed in one iteration (10,000)
  integer :: ncolumns = 50                       ! Number of subcolumns in SCOPS (50)
  integer :: nlr = 40                            ! Number of levels in statistical outputs 
                                                 ! (only used if USE_VGRID=.true.)  (40)
  logical :: use_vgrid = .true.                  ! Use fixed vertical grid for outputs? 
                                                 ! (if .true. then define # of levels with nlr)  (.true.)
  logical :: csat_vgrid = .true.                 ! CloudSat vertical grid? 

  ! Variables for COSP input related to radar simulator
  real(r8) :: radar_freq = 94.0_r8               ! CloudSat radar frequency (GHz) (94.0)
  integer :: surface_radar = 0                   ! surface=1, spaceborne=0 (0)
  integer :: use_gas_abs = 1                     ! include gaseous absorption? yes=1,no=0 (1)
  integer :: do_ray = 0                          ! calculate/output Rayleigh refl=1, not=0 (0)
  real(r8) :: k2 = -1                            ! |K|^2, -1=use frequency dependent default (-1)

  ! Variables for COSP input related to lidar simulator
  integer, parameter :: Nprmts_max_hydro = 12    ! Max # params for hydrometeor size distributions (12)
  integer, parameter :: Naero = 1                ! Number of aerosol species (Not used) (1)
  integer, parameter :: Nprmts_max_aero = 1      ! Max # params for aerosol size distributions (not used) (1)
  integer :: lidar_ice_type = 0                  ! Ice particle shape in lidar calculations
                                                 ! (0=ice-spheres ; 1=ice-non-spherical) (0)
  integer, parameter :: overlap = 3              ! overlap type: 1=max, 2=rand, 3=max/rand (3)

  ! Variables for COSP input related to ISCCP simulator
  integer :: isccp_topheight = 1                 ! 1 = adjust top height using both a computed infrared
                                                 ! brightness temperature and the visible
                                                 ! optical depth to adjust cloud top pressure.
                                                 ! Note that this calculation is most appropriate to compare
                                                 ! to ISCCP data during sunlit hours.
                                                 ! 2 = do not adjust top height, that is cloud top pressure
                                                 ! is the actual cloud top pressure in the model
                                                 ! 3 = adjust top height using only the computed infrared
                                                 ! brightness temperature. Note that this calculation is most
                                                 ! appropriate to compare to ISCCP IR only algortihm (i.e.
                                                 ! you can compare to nighttime ISCCP data with this option) (1)
  integer :: isccp_topheight_direction = 2       ! direction for finding atmosphere pressure level with
                                                 ! interpolated temperature equal to the radiance
                                                 ! determined cloud-top temperature
                                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                                 ! with interpolated temperature
                                                 ! equal to the radiance determined cloud-top temperature
                                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                                 ! 1 = default setting in COSP v1.1, matches all versions of
                                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                                 ! 2 = default setting in COSP v1.3. default since V4.0 of ISCCP simulator

  ! ######################################################################################
  ! Other variables
  ! ######################################################################################  
  logical,allocatable :: first_run_cosp(:)      !.true. if run_cosp has been populated (allocatable->begchunk:endchunk)
  logical,allocatable :: run_cosp(:,:)          !.true. if cosp should be run by column and
                                                !       chunk (allocatable->1:pcols,begchunk:endchunk)
  ! pbuf indices
  integer :: cld_idx, concld_idx, lsreffrain_idx, lsreffsnow_idx, cvreffliq_idx
  integer :: cvreffice_idx
  integer :: gb_totcldliqmr_idx, gb_totcldicemr_idx
  integer :: dpflxprc_idx
  integer :: dpflxsnw_idx, shflxprc_idx, shflxsnw_idx, lsflxprc_idx, lsflxsnw_idx
  integer :: rei_idx, rel_idx
  
  ! ######################################################################################
  ! Declarations specific to COSP2
  ! ######################################################################################
  type(radar_cfg)              :: rcfg_cloudsat ! Radar configuration (Cloudsat)
  type(radar_cfg), allocatable :: rcfg_cs(:)    ! chunked version of rcfg_cloudsat
  type(size_distribution)              :: sd       ! Size distribution used by radar simulator
  type(size_distribution), allocatable :: sd_cs(:) ! chunked version of sd
  character(len=64)         :: cloudsat_micro_scheme = 'MMF_v3.5_single_moment'
  
  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel
  
  ! Stratiform and convective clouds in frac_out (scops output).
  integer, parameter :: &
       I_LSC = 1, & ! Large-scale clouds
       I_CVC = 2    ! Convective clouds    
  
  ! Microphysical settings for the precipitation flux to mixing ratio conversion
  real(r8),parameter,dimension(nhydro) :: &
                 !  LSL     LSI         LSR         LSS       CVL     CVI        CVR          CVS          LSG
       N_ax    = (/-1._r8, -1._r8,     8.e6_r8,     3.e6_r8, -1._r8, -1._r8,     8.e6_r8,     3.e6_r8,     4.e6_r8/),&
       N_bx    = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       alpha_x = (/-1._r8, -1._r8,      0.0_r8,      0.0_r8, -1._r8, -1._r8,      0.0_r8,      0.0_r8,      0.0_r8/),&
       c_x     = (/-1._r8, -1._r8,    842.0_r8,     4.84_r8, -1._r8, -1._r8,    842.0_r8,     4.84_r8,     94.5_r8/),&
       d_x     = (/-1._r8, -1._r8,      0.8_r8,     0.25_r8, -1._r8, -1._r8,      0.8_r8,     0.25_r8,      0.5_r8/),&
       g_x     = (/-1._r8, -1._r8,      0.5_r8,      0.5_r8, -1._r8, -1._r8,      0.5_r8,      0.5_r8,      0.5_r8/),&
       a_x     = (/-1._r8, -1._r8,    524.0_r8,    52.36_r8, -1._r8, -1._r8,    524.0_r8,    52.36_r8,   209.44_r8/),&
       b_x     = (/-1._r8, -1._r8,      3.0_r8,      3.0_r8, -1._r8, -1._r8,      3.0_r8,      3.0_r8,      3.0_r8/),&
       gamma_1 = (/-1._r8, -1._r8, 17.83725_r8, 8.284701_r8, -1._r8, -1._r8, 17.83725_r8, 8.284701_r8, 11.63230_r8/),&
       gamma_2 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/),&
       gamma_3 = (/-1._r8, -1._r8,      2.0_r8,      2.0_r8, -1._r8, -1._r8,      2.0_r8,      2.0_r8,      2.0_r8/),&
       gamma_4 = (/-1._r8, -1._r8,      6.0_r8,      6.0_r8, -1._r8, -1._r8,      6.0_r8,      6.0_r8,      6.0_r8/)       
#endif

CONTAINS

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_readnl
  ! ######################################################################################
  subroutine cospsimulator_intr_readnl(nlfile)
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,    only: mpicom, mpilog, mpiint
#endif

    character(len=*), intent(in) :: nlfile  ! file containing namelist input  (nlfile=atm_in)

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cospsimulator_intr_readnl'

#ifdef USE_COSP
    namelist /cospsimulator_nl/ docosp, cosp_ncolumns, cosp_nradsteps,           &
       cosp_amwg, cosp_lite, cosp_passive, cosp_active, cosp_isccp, cosp_runall, &
       cosp_lfrac_out, cosp_lradar_sim, cosp_llidar_sim, cosp_lisccp_sim,        &
       cosp_lmisr_sim, cosp_lmodis_sim,                                          &
       cosp_histfile_num, cosp_histfile_aux, cosp_histfile_aux_num
    
    !! read in the namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'cospsimulator_nl', status=ierr)   
       if (ierr == 0) then
          read(unitn, cospsimulator_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
    
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(docosp,               1,  mpilog, 0, mpicom)
    call mpibcast(cosp_amwg,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lite,            1,  mpilog, 0, mpicom)
    call mpibcast(cosp_passive,         1,  mpilog, 0, mpicom)
    call mpibcast(cosp_active,          1,  mpilog, 0, mpicom)
    call mpibcast(cosp_isccp,           1,  mpilog, 0, mpicom)
    call mpibcast(cosp_runall,          1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lfrac_out,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lradar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_llidar_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lisccp_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmisr_sim,       1,  mpilog, 0, mpicom)
    call mpibcast(cosp_lmodis_sim,      1,  mpilog, 0, mpicom)
    call mpibcast(cosp_ncolumns,        1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_num,    1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux_num,1,  mpiint, 0, mpicom)
    call mpibcast(cosp_histfile_aux,    1,  mpilog, 0, mpicom)
    call mpibcast(cosp_nradsteps,       1,  mpiint, 0, mpicom)
#endif   
    
    if (cosp_lfrac_out) then
       lfrac_out = .true.
    end if
    if (cosp_lradar_sim) then
       lradar_sim = .true.
    end if
    if (cosp_llidar_sim) then
       llidar_sim = .true.
       lparasol_sim = .true.
    end if
    if (cosp_lisccp_sim) then
       lisccp_sim = .true.
    end if
    if (cosp_lmisr_sim) then
       lmisr_sim = .true.
    end if
    if (cosp_lmodis_sim) then
       lmodis_sim = .true.
    end if
    
    if (cosp_histfile_aux .and. cosp_histfile_aux_num == -1) then
       cosp_histfile_aux_num = cosp_histfile_num
    end if
    
    if (cosp_lite) then
       llidar_sim = .true.
       lparasol_sim = .true.
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_passive) then
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_active) then
       lradar_sim = .true.
       llidar_sim = .true.
       lparasol_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_isccp) then
       lisccp_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    if (cosp_runall) then
       lradar_sim = .true.
       llidar_sim = .true.
       lparasol_sim = .true.
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       lfrac_out = .true.
    end if
    
    !! if no simulators are turned on at all and docosp is, set cosp_amwg = .true.
    if((docosp) .and. (.not.lradar_sim) .and. (.not.llidar_sim) .and. (.not.lisccp_sim) .and. &
         (.not.lmisr_sim) .and. (.not.lmodis_sim)) then
       cosp_amwg = .true.
    end if
    if (cosp_amwg) then
       lradar_sim = .true.
       llidar_sim = .true.
       lparasol_sim = .true.
       lisccp_sim = .true.
       lmisr_sim = .true.
       lmodis_sim = .true.
       cosp_ncolumns = 10
       cosp_nradsteps = 3
    end if
    
    ! Set number of sub-columns, from namelist
    ncolumns   = cosp_ncolumns
    nscol_cosp = cosp_ncolumns
        
    if (masterproc) then
       if (docosp) then 
          write(iulog,*)'COSP configuration:'
          write(iulog,*)'  Number of COSP subcolumns                = ', cosp_ncolumns
          write(iulog,*)'  COSP frequency in radiation steps        = ', cosp_nradsteps
          write(iulog,*)'  Enable radar simulator                   = ', lradar_sim
          write(iulog,*)'  Enable calipso simulator                 = ', llidar_sim
          write(iulog,*)'  Enable ISCCP simulator                   = ', lisccp_sim
          write(iulog,*)'  Enable MISR simulator                    = ', lmisr_sim
          write(iulog,*)'  Enable MODIS simulator                   = ', lmodis_sim
          write(iulog,*)'  RADAR_SIM microphysics scheme            = ', trim(cloudsat_micro_scheme)
          write(iulog,*)'  Write COSP output to history file        = ', cosp_histfile_num
          write(iulog,*)'  Write COSP input fields                  = ', cosp_histfile_aux
          write(iulog,*)'  Write COSP input fields to history file  = ', cosp_histfile_aux_num
          write(iulog,*)'  Write COSP subcolumn fields              = ', lfrac_out
       else
          write(iulog,*)'COSP not enabled'
       end if
    end if
#endif
  end subroutine cospsimulator_intr_readnl

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_register
  ! ######################################################################################
  subroutine cospsimulator_intr_register()

    ! The coordinate variables used for COSP output are defined here.  This
    ! needs to be done before the call to read_restart_history in order for
    ! restarts to work.

    use cam_history_support, only: add_hist_coord
    !---------------------------------------------------------------------------
    
#ifdef USE_COSP
    ! Set number of levels used by COSP to the number of levels used by
    ! CAM's cloud macro/microphysics parameterizations.
    nlay = pver - ktop + 1
    nlayp = nlay + 1

    ! Set COSP coordinate arrays
    call setcosp2values()

    ! Define coordinate variables for COSP outputs.
    if (lisccp_sim .or. lmodis_sim) then
       call add_hist_coord('cosp_prs', nprs_cosp, 'COSP Mean ISCCP pressure',  &
            'hPa', prsmid_cosp, bounds_name='cosp_prs_bnds', bounds=prslim_cosp)
    end if
    
    if (lisccp_sim .or. lmisr_sim) then
       call add_hist_coord('cosp_tau', ntau_cosp,                              &
            'COSP Mean ISCCP optical depth', '1', taumid_cosp,                 &
            bounds_name='cosp_tau_bnds', bounds=taulim_cosp)
    end if
    
    if (lisccp_sim .or. llidar_sim .or. lradar_sim .or. lmisr_sim) then
       call add_hist_coord('cosp_scol', nscol_cosp, 'COSP subcolumn',          &
            values=scol_cosp)
    end if
    
    if (llidar_sim .or. lradar_sim) then
       call add_hist_coord('cosp_ht', nht_cosp,                                &
            'COSP Mean Height for calipso and radar simulator outputs', 'm',   &
            htmid_cosp, bounds_name='cosp_ht_bnds', bounds=htlim_cosp,         &
            vertical_coord=.true.)
    end if
    
    if (llidar_sim) then
       call add_hist_coord('cosp_sr', nsr_cosp,                                &
            'COSP Mean Scattering Ratio for calipso simulator CFAD output', '1', &
            srmid_cosp, bounds_name='cosp_sr_bnds', bounds=srlim_cosp)
    end if
    
    if (llidar_sim) then
       call add_hist_coord('cosp_sza', nsza_cosp, 'COSP Parasol SZA',          &
            'degrees', sza_cosp)
    end if
    
    if (lradar_sim) then
       call add_hist_coord('cosp_dbze', CLOUDSAT_DBZE_BINS,                    &
            'COSP Mean dBZe for radar simulator CFAD output', 'dBZ',           &
            dbzemid_cosp, bounds_name='cosp_dbze_bnds', bounds=dbzelim_cosp)
    end if
    
    if (lmisr_sim) then
       call add_hist_coord('cosp_htmisr', nhtmisr_cosp, 'COSP MISR height', &
            'km', htmisrmid_cosp,                                           &
            bounds_name='cosp_htmisr_bnds', bounds=htmisrlim_cosp)
    end if
    
    if (lmodis_sim) then
       call add_hist_coord('cosp_tau_modis', ntau_cosp_modis,                  &
            'COSP Mean MODIS optical depth', '1', taumid_cosp_modis,           &
            bounds_name='cosp_tau_modis_bnds', bounds=taulim_cosp_modis)
       call add_hist_coord('cosp_reffice',numMODISReffIceBins,                 &
            'COSP Mean MODIS effective radius (ice)', 'microns', reffICE_binCenters_cosp, &
            bounds_name='cosp_reffice_bnds',bounds=reffICE_binEdges_cosp)
       call add_hist_coord('cosp_reffliq',numMODISReffLiqBins,                 &
            'COSP Mean MODIS effective radius (liquid)', 'microns', reffLIQ_binCenters_cosp, &
            bounds_name='cosp_reffliq_bnds',bounds=reffLIQ_binEdges_cosp)      
    end if
    
#endif
  end subroutine cospsimulator_intr_register
  
  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_init
  ! ######################################################################################
  subroutine cospsimulator_intr_init()

#ifdef USE_COSP     

    use cam_history,         only: addfld, add_default, horiz_only
    use physics_buffer,      only: pbuf_get_index

    integer :: i, ierr, istat
    character(len=*), parameter :: sub = 'cospsimulator_intr_init'
    !---------------------------------------------------------------------------
    
    ! The COSP init method (setcosp2values) was run from cospsimulator_intr_register in order to add
    ! the history coordinate variables earlier as needed for the restart time sequencing.

    ! ISCCP OUTPUTS
    if (lisccp_sim) then
       call addfld('FISCCP1_COSP', (/'cosp_tau','cosp_prs'/), 'A', 'percent', &
            'Grid-box fraction covered by each ISCCP D level cloud type', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_ISCCP', horiz_only, 'A', 'percent', &
            'Total Cloud Fraction Calculated by the ISCCP Simulator ',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MEANCLDALB_ISCCP', horiz_only, 'A', '1', &
            'Mean cloud albedo*CLDTOT_ISCCP', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MEANPTOP_ISCCP', horiz_only, 'A', 'Pa', &
            'Mean cloud top pressure*CLDTOT_ISCCP',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MEANTAU_ISCCP', horiz_only, 'A', '1', &
            'Mean optical thickness*CLDTOT_ISCCP',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MEANTB_ISCCP', horiz_only, 'A', 'K', &
            'Mean Infrared Tb from ISCCP simulator',flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MEANTBCLR_ISCCP', horiz_only, 'A', 'K', &
            'Mean Clear-sky Infrared Tb from ISCCP simulator', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAU_ISCCP', (/'cosp_scol'/), 'I', '1', &
            'Optical Depth in each Subcolumn', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDPTOP_ISCCP', (/'cosp_scol'/), 'I', 'Pa', &
            'Cloud Top Pressure in each Subcolumn', flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default('FISCCP1_COSP',cosp_histfile_num,' ')
       call add_default('CLDTOT_ISCCP',cosp_histfile_num,' ')
       call add_default('MEANCLDALB_ISCCP',cosp_histfile_num,' ')
       call add_default('MEANPTOP_ISCCP',cosp_histfile_num,' ')
       call add_default('MEANTAU_ISCCP',cosp_histfile_num,' ')
       call add_default('MEANTB_ISCCP',cosp_histfile_num,' ')
       call add_default('MEANTBCLR_ISCCP',cosp_histfile_num,' ')
      
    end if

    ! CALIPSO SIMULATOR OUTPUTS
    if (llidar_sim) then
       call addfld('CLDLOW_CAL', horiz_only, 'A', 'percent', &
            'Calipso Low-level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL', horiz_only, 'A', 'percent', &
            'Calipso Mid-level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL', horiz_only, 'A', 'percent', &
            'Calipso High-level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL', horiz_only, 'A', 'percent', &
            'Calipso Total Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL', (/'cosp_ht'/), 'A', 'percent', &
            'Calipso Cloud Fraction (532 nm)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('RFL_PARASOL', (/'cosp_sza'/), 'A', 'fraction', &
            'PARASOL-like mono-directional reflectance ', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CFAD_SR532_CAL', (/'cosp_sr','cosp_ht'/), 'A', 'fraction', &
            'Calipso Scattering Ratio CFAD (532 nm)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('MOL532_CAL', (/'trop_pref'/), 'A', 'm-1 sr-1', &
            'Calipso Molecular Backscatter (532 nm) ', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('ATB532_CAL', (/'cosp_scol','trop_pref'/), 'I', 'no_unit_log10(x)', &
            'Calipso Attenuated Total Backscatter (532 nm) in each Subcolumn', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_LIQ', (/'cosp_ht'/), 'A', 'percent', &
            'Calipso Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_ICE', (/'cosp_ht'/), 'A', 'percent', &
            'Calipso Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_UN', (/'cosp_ht'/), 'A', 'percent',  &
            'Calipso Undefined-Phase Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMP', (/'cosp_ht'/), 'A', 'K', &
            'Calipso Cloud Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPLIQ', (/'cosp_ht'/), 'A', 'K', &
            'Calipso Liquid Cloud Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPICE', (/'cosp_ht'/), 'A', 'K', &
            'Calipso Ice Cloud Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_TMPUN', (/'cosp_ht'/), 'A', 'K', &
            'Calipso Undefined-Phase Cloud Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_ICE', horiz_only, 'A', 'percent', &
            'Calipso Total Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_LIQ', horiz_only, 'A', 'percent', &
            'Calipso Total Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CAL_UN', horiz_only, 'A', 'percent', &
            'Calipso Total Undefined-Phase Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_ICE', horiz_only, 'A', 'percent', &
            'Calipso High-level Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_LIQ', horiz_only, 'A', 'percent', &
            'Calipso High-level Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDHGH_CAL_UN', horiz_only, 'A', 'percent', &
            'Calipso High-level Undefined-Phase Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_ICE', horiz_only, 'A', 'percent', &
            'Calipso Mid-level Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_LIQ', horiz_only, 'A', 'percent', &
            'Calipso Mid-level Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDMED_CAL_UN', horiz_only, 'A', 'percent', &
            'Calipso Mid-level Undefined-Phase Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_ICE', horiz_only, 'A', 'percent', &
            'Calipso Low-level Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_LIQ', horiz_only, 'A', 'percent', &
            'Calipso Low-level Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDLOW_CAL_UN', horiz_only, 'A', 'percent', &
            'Calipso Low-level Undefined-Phase Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
    
       call add_default('CLDLOW_CAL',cosp_histfile_num,' ')
       call add_default('CLDMED_CAL',cosp_histfile_num,' ')
       call add_default('CLDHGH_CAL',cosp_histfile_num,' ')
       call add_default('CLDTOT_CAL',cosp_histfile_num,' ')
       call add_default('CLD_CAL',cosp_histfile_num,' ')
       call add_default('RFL_PARASOL',cosp_histfile_num,' ')
       call add_default('CFAD_SR532_CAL',cosp_histfile_num,' ')
       call add_default('CLD_CAL_LIQ',cosp_histfile_num,' ')
       call add_default('CLD_CAL_ICE',cosp_histfile_num,' ')
       call add_default('CLD_CAL_UN',cosp_histfile_num,' ')
       call add_default('CLDTOT_CAL_ICE',cosp_histfile_num,' ')
       call add_default('CLDTOT_CAL_LIQ',cosp_histfile_num,' ')
       call add_default('CLDTOT_CAL_UN',cosp_histfile_num,' ')
       call add_default('CLDHGH_CAL_ICE',cosp_histfile_num,' ')
       call add_default('CLDHGH_CAL_LIQ',cosp_histfile_num,' ')
       call add_default('CLDHGH_CAL_UN',cosp_histfile_num,' ')
       call add_default('CLDMED_CAL_ICE',cosp_histfile_num,' ')
       call add_default('CLDMED_CAL_LIQ',cosp_histfile_num,' ')
       call add_default('CLDMED_CAL_UN',cosp_histfile_num,' ')
       call add_default('CLDLOW_CAL_ICE',cosp_histfile_num,' ')
       call add_default('CLDLOW_CAL_LIQ',cosp_histfile_num,' ')
       call add_default('CLDLOW_CAL_UN',cosp_histfile_num,' ')

       if ((.not.cosp_amwg) .and. (.not.cosp_lite) .and. (.not.cosp_passive) .and. (.not.cosp_active) &
            .and. (.not.cosp_isccp)) then
          call add_default('MOL532_CAL',cosp_histfile_num,' ')
       end if
    end if

    ! RADAR SIMULATOR OUTPUTS
    allocate(sd_cs(begchunk:endchunk), rcfg_cs(begchunk:endchunk), stat=istat)
    call handle_allocate_error(istat, sub, 'sd_cs,rcfg_cs')
    if (lradar_sim) then

       do i = begchunk, endchunk
          sd_cs(i)   = sd
          rcfg_cs(i) = rcfg_cloudsat
       end do

       call addfld('CFAD_DBZE94_CS',(/'cosp_dbze','cosp_ht  '/), 'A', 'fraction', &
            'Radar Reflectivity Factor CFAD (94 GHz)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLD_CAL_NOTCS', (/'cosp_ht'/), 'A', 'percent', &
            'Cloud occurrence seen by CALIPSO but not CloudSat ', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CALCS', horiz_only, 'A', 'percent', &
            'Calipso and Radar Total Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CS', horiz_only, 'A', 'percent', &
            'Radar total cloud amount', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLDTOT_CS2', horiz_only, 'A', 'percent', &
            'Radar total cloud amount without the data for the first kilometer above surface ', &
            flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('DBZE_CS', (/'cosp_scol','trop_pref'/), 'I', 'dBZe', &
            'Radar dBZe (94 GHz) in each Subcolumn', flag_xyfill=.true., fill_value=R_UNDEF)

       ! Cloudsat near-sfc precipitation diagnostics
       call addfld('CS_NOPRECIP',  horiz_only, 'A', '1', &
            'CloudSat No Rain Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINPOSS',  horiz_only, 'A', '1', &
            'Cloudsat Rain Possible Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINPROB',  horiz_only, 'A', '1', &
            'CloudSat Rain Probable Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINCERT',  horiz_only, 'A', '1', &
            'CloudSat Rain Certain Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_SNOWPOSS',  horiz_only, 'A', '1', &
            'CloudSat Snow Possible Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_SNOWCERT',  horiz_only, 'A', '1', &
            'CloudSat Snow Certain Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_MIXPOSS',   horiz_only, 'A', '1', &
            'CloudSat Mixed Possible Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_MIXCERT',   horiz_only, 'A', '1', &
            'CloudSat Mixed Certain Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_RAINHARD',  horiz_only, 'A', '1', &
            'CloudSat Heavy Rain Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_UN',        horiz_only, 'A', '1', &
            'CloudSat Unclassified Precipitation Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CS_PIA',       horiz_only, 'A', 'dBZ', &
            'CloudSat Radar Path Integrated Attenuation', flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default('CFAD_DBZE94_CS',cosp_histfile_num,' ')
       call add_default('CLD_CAL_NOTCS', cosp_histfile_num,' ')
       call add_default('CLDTOT_CALCS',  cosp_histfile_num,' ')
       call add_default('CLDTOT_CS',     cosp_histfile_num,' ')
       call add_default('CLDTOT_CS2',    cosp_histfile_num,' ')
       call add_default('CS_NOPRECIP',   cosp_histfile_num,' ')
       call add_default('CS_RAINPOSS',   cosp_histfile_num,' ')
       call add_default('CS_RAINPROB',   cosp_histfile_num,' ')
       call add_default('CS_RAINCERT',   cosp_histfile_num,' ')
       call add_default('CS_SNOWPOSS',   cosp_histfile_num,' ')
       call add_default('CS_SNOWCERT',   cosp_histfile_num,' ')
       call add_default('CS_MIXPOSS',    cosp_histfile_num,' ')
       call add_default('CS_MIXCERT',    cosp_histfile_num,' ')
       call add_default('CS_RAINHARD',   cosp_histfile_num,' ')
       call add_default('CS_UN',         cosp_histfile_num,' ')
       call add_default('CS_PIA',        cosp_histfile_num,' ')
    end if
    
    ! MISR SIMULATOR OUTPUTS
    if (lmisr_sim) then
       call addfld('CLD_MISR', (/'cosp_tau   ','cosp_htmisr'/), 'A', 'percent', &
            'Cloud Fraction from MISR Simulator', flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default('CLD_MISR',cosp_histfile_num,' ')
    end if

    ! MODIS OUTPUT
    if (lmodis_sim) then
       call addfld('CLTMODIS', horiz_only, 'A', '%', &
            'MODIS Total Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLWMODIS', horiz_only, 'A', '%', &
            'MODIS Liquid Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLIMODIS', horiz_only, 'A', '%', &
            'MODIS Ice Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLHMODIS', horiz_only, 'A', '%', &
            'MODIS High Level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLMMODIS', horiz_only, 'A', '%', &
            'MODIS Mid Level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLLMODIS', horiz_only, 'A', '%', &
            'MODIS Low Level Cloud Fraction', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUTMODIS', horiz_only, 'A', '1', &
            'MODIS Total Cloud Optical Thickness*CLTMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUWMODIS', horiz_only, 'A', '1', &
            'MODIS Liquid Cloud Optical Thickness*CLWMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUIMODIS', horiz_only, 'A', '1', &
            'MODIS Ice Cloud Optical Thickness*CLIMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUTLOGMODIS', horiz_only, 'A', '1', &
            'MODIS Total Cloud Optical Thickness (Log10 Mean)*CLTMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUWLOGMODIS', horiz_only, 'A', '1', &
            'MODIS Liquid Cloud Optical Thickness (Log10 Mean)*CLWMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('TAUILOGMODIS', horiz_only, 'A', '1', &
            'MODIS Ice Cloud Optical Thickness (Log10 Mean)*CLIMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('REFFCLWMODIS', horiz_only, 'A', 'm', &
            'MODIS Liquid Cloud Particle Size*CLWMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('REFFCLIMODIS', horiz_only, 'A', 'm', &
            'MODIS Ice Cloud Particle Size*CLIMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('PCTMODIS', horiz_only, 'A', 'Pa', &
            'MODIS Cloud Top Pressure*CLTMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('LWPMODIS', horiz_only, 'A', 'kg m-2', &
            'MODIS Cloud Liquid Water Path*CLWMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('IWPMODIS', horiz_only, 'A', 'kg m-2', &
            'MODIS Cloud Ice Water Path*CLIMODIS', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLMODIS', (/'cosp_tau_modis','cosp_prs      '/), 'A', '%', &
            'MODIS Cloud Area Fraction (tau-pressure histogram)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLRIMODIS', (/'cosp_tau_modis','cosp_reffice  '/), 'A', '%', &
            'MODIS Cloud Area Fraction (tau-reffice histogram)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld('CLRLMODIS', (/'cosp_tau_modis','cosp_reffliq  '/), 'A', '%', &
            'MODIS Cloud Area Fraction (tau-reffliq histogram)', flag_xyfill=.true., fill_value=R_UNDEF)
       
       call add_default('CLTMODIS',cosp_histfile_num,' ')
       call add_default('CLWMODIS',cosp_histfile_num,' ')
       call add_default('CLIMODIS',cosp_histfile_num,' ')
       call add_default('CLHMODIS',cosp_histfile_num,' ')
       call add_default('CLMMODIS',cosp_histfile_num,' ')
       call add_default('CLLMODIS',cosp_histfile_num,' ')
       call add_default('TAUTMODIS',cosp_histfile_num,' ')
       call add_default('TAUWMODIS',cosp_histfile_num,' ')
       call add_default('TAUIMODIS',cosp_histfile_num,' ')
       call add_default('TAUTLOGMODIS',cosp_histfile_num,' ')
       call add_default('TAUWLOGMODIS',cosp_histfile_num,' ')
       call add_default('TAUILOGMODIS',cosp_histfile_num,' ')
       call add_default('REFFCLWMODIS',cosp_histfile_num,' ')
       call add_default('REFFCLIMODIS',cosp_histfile_num,' ')
       call add_default('PCTMODIS',cosp_histfile_num,' ')
       call add_default('LWPMODIS',cosp_histfile_num,' ')
       call add_default('IWPMODIS',cosp_histfile_num,' ')
       call add_default('CLMODIS',cosp_histfile_num,' ')
       call add_default('CLRIMODIS',cosp_histfile_num,' ')
       call add_default('CLRLMODIS',cosp_histfile_num,' ')
    end if
    
    ! SUB-COLUMN OUTPUT
    if (lfrac_out) then
       call addfld('SCOPS_OUT', (/'cosp_scol','trop_pref'/), 'I', '0=nocld,1=strcld,2=cnvcld', &
            'SCOPS Subcolumn output', flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default('SCOPS_OUT',cosp_histfile_num,' ')

       if (lisccp_sim) then
          call add_default('TAU_ISCCP',cosp_histfile_num,' ')
          call add_default('CLDPTOP_ISCCP',cosp_histfile_num,' ')
       end if

       if (llidar_sim) then
          call add_default('ATB532_CAL',cosp_histfile_num,' ')
       end if

       if (lradar_sim) then
          call add_default('DBZE_CS',cosp_histfile_num,' ')
       end if
    end if
    
    !! ADDFLD, ADD_DEFAULT, OUTFLD CALLS FOR COSP OUTPUTS IF RUNNING COSP OFF-LINE
    if (cosp_histfile_aux) then
       call addfld ('PS_COSP',         horiz_only,            'I','Pa', &
          'COSP Surface Pressure', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('TS_COSP',         horiz_only,            'I','K',  &
          'COSP Skin Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('P_COSP',          (/            'trop_pref'/), 'I','Pa', &
          'COSP Pressure (layer midpoint)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('PH_COSP',         (/            'trop_prefi'/), 'I','Pa', &
          'COSP Pressure (layer interface)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_COSP',       (/            'trop_pref'/), 'I','m',  &
          'COSP Height (layer midpoint)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('ZLEV_HALF_COSP',  (/            'trop_prefi'/), 'I','m',  &
          'COSP Height (layer interface)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('T_COSP',          (/            'trop_pref'/), 'I','K',  &
          'COSP Temperature', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('Q_COSP',          (/            'trop_pref'/), 'I','percent', &
          'COSP Specific Humidity', flag_xyfill=.true., fill_value=R_UNDEF)

       call addfld ('TAU_067',         (/'cosp_scol','trop_pref'/), 'I','1', &
          'Subcolumn 0.67micron optical depth', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('EMISS_11',        (/'cosp_scol','trop_pref'/), 'I','1', &
          'Subcolumn 11micron emissivity', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_fracliq',   (/'cosp_scol','trop_pref'/), 'I','1', &
          'Fraction of tau from liquid water', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_asym',      (/'cosp_scol','trop_pref'/), 'I','1', &
          'Asymmetry parameter (MODIS)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('MODIS_ssa',       (/'cosp_scol','trop_pref'/), 'I','1', &
          'Single-scattering albedo (MODIS)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_z_vol',        (/'cosp_scol','trop_pref'/), 'I','1', &
          'Effective reflectivity factor (CLOUDSAT)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_kr_vol',       (/'cosp_scol','trop_pref'/), 'I','1', &
          'Attenuation coefficient (hydro) (CLOUDSAT)', flag_xyfill=.true., fill_value=R_UNDEF)
       call addfld ('CS_g_vol',        (/'cosp_scol','trop_pref'/), 'I','1', &
          'Attenuation coefficient (gases) (CLOUDSAT)', flag_xyfill=.true., fill_value=R_UNDEF)

       call add_default('PS_COSP',         cosp_histfile_aux_num,' ')
       call add_default('TS_COSP',         cosp_histfile_aux_num,' ')
       call add_default('P_COSP',          cosp_histfile_aux_num,' ')
       call add_default('PH_COSP',         cosp_histfile_aux_num,' ')
       call add_default('ZLEV_COSP',       cosp_histfile_aux_num,' ')
       call add_default('ZLEV_HALF_COSP',  cosp_histfile_aux_num,' ')
       call add_default('T_COSP',          cosp_histfile_aux_num,' ')
       call add_default('Q_COSP',          cosp_histfile_aux_num,' ')
       call add_default('TAU_067',         cosp_histfile_aux_num,' ')
       call add_default('EMISS_11',        cosp_histfile_aux_num,' ')
       call add_default('MODIS_fracliq',   cosp_histfile_aux_num,' ')
       call add_default('MODIS_asym',      cosp_histfile_aux_num,' ')
       call add_default('MODIS_ssa',       cosp_histfile_aux_num,' ')
       call add_default('CS_z_vol',        cosp_histfile_aux_num,' ')
       call add_default('CS_kr_vol',       cosp_histfile_aux_num,' ')
       call add_default('CS_g_vol',        cosp_histfile_aux_num,' ')
    end if
    
    rei_idx        = pbuf_get_index('REI')
    rel_idx        = pbuf_get_index('REL')
    cld_idx        = pbuf_get_index('CLD')
    concld_idx     = pbuf_get_index('CONCLD')
    lsreffrain_idx = pbuf_get_index('LS_REFFRAIN')
    lsreffsnow_idx = pbuf_get_index('LS_REFFSNOW')
    cvreffliq_idx  = pbuf_get_index('CV_REFFLIQ')
    cvreffice_idx  = pbuf_get_index('CV_REFFICE')
    gb_totcldliqmr_idx = pbuf_get_index('GB_TOTCLDLIQMR') ! grid box total cloud liquid water mr (kg/kg)
    gb_totcldicemr_idx = pbuf_get_index('GB_TOTCLDICEMR') ! grid box total cloud ice water mr (kg/kg)
    dpflxprc_idx   = pbuf_get_index('DP_FLXPRC')
    dpflxsnw_idx   = pbuf_get_index('DP_FLXSNW')
    shflxprc_idx   = pbuf_get_index('SH_FLXPRC', errcode=ierr)
    shflxsnw_idx   = pbuf_get_index('SH_FLXSNW', errcode=ierr)
    lsflxprc_idx   = pbuf_get_index('LS_FLXPRC')
    lsflxsnw_idx   = pbuf_get_index('LS_FLXSNW')
    
    allocate(first_run_cosp(begchunk:endchunk), run_cosp(1:pcols,begchunk:endchunk), &
       stat=istat)
    call handle_allocate_error(istat, sub, '*run_cosp')
    first_run_cosp(begchunk:endchunk)=.true.
    run_cosp(1:pcols,begchunk:endchunk)=.false.

#endif    
  end subroutine cospsimulator_intr_init

  ! ######################################################################################
  ! SUBROUTINE setcosp2values
  ! ######################################################################################
#ifdef USE_COSP
  subroutine setcosp2values()
    use mod_cosp,             only: cosp_init 
    use mod_cosp_config,      only: vgrid_zl, vgrid_zu, vgrid_z
    use mod_quickbeam_optics, only: hydro_class_init, quickbeam_optics_init
    
    ! Local
    logical :: ldouble=.false.
    logical :: lsingle=.true. ! Default is to use single moment
    integer :: k
    integer :: istat
    character(len=*), parameter :: sub = 'setcosp2values'
    !--------------------------------------------------------------------------------------

    prsmid_cosp  = pres_binCenters
    prslim_cosp  = pres_binEdges
    taumid_cosp  = tau_binCenters
    taulim_cosp  = tau_binEdges
    srmid_cosp   = calipso_binCenters
    srlim_cosp   = calipso_binEdges
    sza_cosp     = parasol_sza
    dbzemid_cosp = cloudsat_binCenters
    dbzelim_cosp = cloudsat_binEdges
    htmisrmid_cosp = misr_histHgtCenters
    htmisrlim_cosp = misr_histHgtEdges
    taumid_cosp_modis = tau_binCenters
    taulim_cosp_modis = tau_binEdges
    reffICE_binCenters_cosp = reffICE_binCenters
    reffICE_binEdges_cosp   = reffICE_binEdges
    reffLIQ_binCenters_cosp = reffLIQ_binCenters
    reffLIQ_binEdges_cosp   = reffLIQ_binEdges
                                  
    ! Initialize the distributional parameters for hydrometeors in radar simulator. In COSPv1.4, this was declared in
    ! cosp_defs.f.
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true. 
       lsingle = .false.
    endif
    call hydro_class_init(lsingle,ldouble,sd)
    call quickbeam_optics_init()

    ! DS2017: The setting up of the vertical grid for regridding the CALIPSO and Cloudsat products is 
    !         now done in cosp_init, but these fields are stored in cosp_config.F90.
    !         Additionally all static fields used by the individual simulators are set up by calls
    !         to _init functions in cosp_init.
    ! DS2019: Add logicals, default=.false., for new Lidar simuldators (Earthcare (atlid) and ground-based
    !         lidar at 532nm)
    call COSP_INIT(Lisccp_sim, Lmodis_sim, Lmisr_sim, Lradar_sim, Llidar_sim, LgrLidar532, &
         Latlid, Lparasol_sim, Lrttov_sim, radar_freq, k2, use_gas_abs, do_ray,              &
         isccp_topheight, isccp_topheight_direction, surface_radar, rcfg_cloudsat,           &
         use_vgrid, csat_vgrid, Nlr, nlay, cloudsat_micro_scheme)

    if (use_vgrid) then      !! using fixed vertical grid
       if (csat_vgrid) then
          nht_cosp = 40
       else
          nht_cosp = Nlr
       endif
    endif
    
    ! DJS2017: In COSP2, most of the bin boundaries, centers, and edges are declared in src/cosp_config.F90.
    !          Above I just assign them accordingly in the USE statement. Other bin bounds needed by CAM 
    !          are calculated here.

    allocate( &
       htmlmid_cosp(nlay), &
       htdbze_dbzemid_cosp(nht_cosp*CLOUDSAT_DBZE_BINS), &
       htlim_cosp(2,nht_cosp), &
       htmid_cosp(nht_cosp), &
       htlim_cosp_1d(nht_cosp+1), &
       htdbze_htmid_cosp(nht_cosp*CLOUDSAT_DBZE_BINS), &
       htsr_htmid_cosp(nht_cosp*nsr_cosp), &
       htsr_srmid_cosp(nht_cosp*nsr_cosp), &
       htmlscol_htmlmid_cosp(nlay*nscol_cosp), &
       htmlscol_scol_cosp(nlay*nscol_cosp), &
       scol_cosp(nscol_cosp), &
       htdbze_cosp(nht_cosp*CLOUDSAT_DBZE_BINS), &
       htsr_cosp(nht_cosp*nsr_cosp), &
       htmlscol_cosp(nlay*nscol_cosp), stat=istat)
    call handle_allocate_error(istat, sub, 'htmlmid_cosp,..,htmlscol_cosp')
    
    ! DJS2017: Just pull from cosp_config
    if (use_vgrid) then
       htlim_cosp_1d(1)            = vgrid_zu(1)
       htlim_cosp_1d(2:nht_cosp+1) = vgrid_zl
    endif
    htmid_cosp      = vgrid_z
    htlim_cosp(1,:) = vgrid_zu
    htlim_cosp(2,:) = vgrid_zl

    scol_cosp(:) = (/(k,k=1,nscol_cosp)/)
    
    !  Just using an index here, model height is a prognostic variable
    htmlmid_cosp(:) = (/(k,k=1,nlay)/)
    
    ! assign mixed dimensions an integer index for cam_history.F90
    do k=1,nprs_cosp*ntau_cosp
       prstau_cosp(k) = k
    end do
    do k=1,nprs_cosp*ntau_cosp_modis
       prstau_cosp_modis(k) = k
    end do
    do k=1,nht_cosp*CLOUDSAT_DBZE_BINS
       htdbze_cosp(k) = k
    end do
    do k=1,nht_cosp*nsr_cosp
       htsr_cosp(k) = k
    end do
    do k=1,nlay*nscol_cosp
       htmlscol_cosp(k) = k
    end do
    do k=1,nhtmisr_cosp*ntau_cosp
       htmisrtau_cosp(k) = k
    end do
    
    ! next, assign collapsed reference vectors for cam_history.F90
    ! convention for saving output = prs1,tau1 ... prs1,tau7 ; prs2,tau1 ... prs2,tau7 etc.
    ! actual output is specified in cospsimulator_intr_init.
    do k=1,nprs_cosp
       prstau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
       prstau_prsmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=prsmid_cosp(k)
       prstau_taumid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=taumid_cosp_modis(1:ntau_cosp_modis)
       prstau_prsmid_cosp_modis(ntau_cosp_modis*(k-1)+1:k*ntau_cosp_modis)=prsmid_cosp(k)
    enddo
    
    do k=1,nht_cosp
       htdbze_dbzemid_cosp(CLOUDSAT_DBZE_BINS*(k-1)+1:k*CLOUDSAT_DBZE_BINS)=dbzemid_cosp(1:CLOUDSAT_DBZE_BINS)
       htdbze_htmid_cosp(CLOUDSAT_DBZE_BINS*(k-1)+1:k*CLOUDSAT_DBZE_BINS)=htmid_cosp(k)
    enddo
    
    do k=1,nht_cosp
       htsr_srmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=srmid_cosp(1:nsr_cosp)
       htsr_htmid_cosp(nsr_cosp*(k-1)+1:k*nsr_cosp)=htmid_cosp(k)
    enddo
    
    do k=1,nlay
       htmlscol_scol_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=scol_cosp(1:nscol_cosp)
       htmlscol_htmlmid_cosp(nscol_cosp*(k-1)+1:k*nscol_cosp)=htmlmid_cosp(k)
    enddo
    
    do k=1,nhtmisr_cosp
       htmisrtau_taumid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=taumid_cosp(1:ntau_cosp)
       htmisrtau_htmisrmid_cosp(ntau_cosp*(k-1)+1:k*ntau_cosp)=htmisrmid_cosp(k)
    enddo
    
  end subroutine setcosp2values
#endif    

  ! ######################################################################################
  ! SUBROUTINE cospsimulator_intr_run
  ! ######################################################################################
  subroutine cospsimulator_intr_run(state, pbuf, cam_in, emis, coszrs,     &
                                    cld_swtau_in, snow_tau_in, snow_emis_in)    

    use physics_types,        only: physics_state
    use physics_buffer,       only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use camsrfexch,           only: cam_in_t
    use constituents,         only: cnst_get_ind
    use rad_constituents,     only: rad_cnst_get_gas
    use interpolate_data,     only: lininterp_init,lininterp,lininterp_finish,interp_type
    use physconst,            only: pi, inverse_gravit => rga
    use cam_history,          only: outfld,hist_fld_col_active 
    use cam_history_support,  only: max_fieldname_len

#ifdef USE_COSP
    use mod_cosp_config,      only: R_UNDEF,parasol_nrefl, Nlvgrid
    use mod_cosp,             only: cosp_simulator
    use mod_quickbeam_optics, only: size_distribution
#endif

    ! ######################################################################################
    ! Inputs
    ! ######################################################################################
    type(physics_state), intent(in),target  :: state
    type(physics_buffer_desc),      pointer :: pbuf(:)
    type(cam_in_t),      intent(in)         :: cam_in
    real(r8), intent(in) :: emis(pcols,pver)                  ! cloud longwave emissivity
    real(r8), intent(in) :: coszrs(pcols)                     ! cosine solar zenith angle (to tell if day or night)
    real(r8), intent(in),optional :: cld_swtau_in(pcols,pver) ! RRTM cld_swtau_in, read in using this variable
    real(r8), intent(in),optional :: snow_tau_in(pcols,pver)  ! RRTM grid-box mean SW snow optical depth, used for CAM5 simulations 
    real(r8), intent(in),optional :: snow_emis_in(pcols,pver) ! RRTM grid-box mean LW snow optical depth, used for CAM5 simulations 

#ifdef USE_COSP
    ! ######################################################################################
    ! Local variables
    ! ######################################################################################
    integer :: lchnk                             ! chunk identifier
    integer :: ncol                              ! number of active atmospheric columns
    integer :: i, k, kk
    integer :: itim_old
    integer :: ip, it
    integer :: ipt
    integer :: ih, ihd, ihs, ihsc, ihm, ihmt, ihml
    integer :: isc
    integer :: is
    integer :: id

    real(r8), parameter :: rad2deg = 180._r8/pi
    
    ! Microphysics variables
    integer :: ixcldliq                                   ! cloud liquid amount index for state%q
    integer :: ixcldice                                   ! cloud ice amount index
    
    ! COSP-related local vars
    type(cosp_outputs)        :: cospOUT                  ! COSP simulator outputs
    type(cosp_optical_inputs) :: cospIN                   ! COSP optical (or derived?) fields needed by simulators
    type(cosp_column_inputs)  :: cospstateIN              ! COSP model fields needed by simulators
    
    ! COSP input variables that depend on CAM
    integer :: Npoints                                    ! Number of gridpoints COSP will process
    real(r8), parameter :: emsfc_lw = 0.99_r8             ! longwave emissivity of surface at 10.5 microns 
    
    ! Local vars related to calculations to go from CAM input to COSP input
    ! cosp convective value includes both deep and shallow convection
    real(r8), allocatable :: &
       zmid(:,:), &           ! layer midpoint height asl (m)
       zint(:,:), &           ! layer interface height asl (m)
       surf_hgt(:), &         ! surface height (m)
       landmask(:), &         ! landmask (0 or 1)
       mr_ccliq(:,:), &       ! mixing_ratio_convective_cloud_liquid (kg/kg)
       mr_ccice(:,:), &       ! mixing_ratio_convective_cloud_ice (kg/kg)
       mr_lsliq(:,:), &       ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
       mr_lsice(:,:), &       ! mixing_ratio_large_scale_cloud_ice (kg/kg)
       rain_cv(:,:), &        ! interface flux_convective_cloud_rain (kg m^-2 s^-1)
       snow_cv(:,:), &        ! interface flux_convective_cloud_snow (kg m^-2 s^-1)
       rain_cv_interp(:,:), & ! midpoint flux_convective_cloud_rain (kg m^-2 s^-1)
       snow_cv_interp(:,:), & ! midpoint flux_convective_cloud_snow (kg m^-2 s^-1)
       rain_ls_interp(:,:), & ! midpoint ls rain flux (kg m^-2 s^-1)
       snow_ls_interp(:,:), & ! midpoint ls snow flux
       grpl_ls_interp(:,:), & ! midpoint ls grp flux, set to 0
       reff_cosp(:,:,:), &    ! effective radius for cosp input
       dtau_s(:,:), &         ! Optical depth of stratiform cloud at 0.67 um
       dtau_c(:,:), &         ! Optical depth of convective cloud at 0.67 um
       dtau_s_snow(:,:), &    ! Grid-box mean Optical depth of stratiform snow at 0.67 um
       dem_s(:,:), &          ! Longwave emis of stratiform cloud at 10.5 um
       dem_c(:,:), &          ! Longwave emis of convective cloud at 10.5 um
       dem_s_snow(:,:)        ! Grid-box mean Optical depth of stratiform snow at 10.5 um

    integer :: cam_sunlit(pcols) ! cam_sunlit - Sunlit flag(1-sunlit/0-dark).
    integer :: nSunLit           ! Number of sunlit (not sunlit) scenes.
    
    ! ######################################################################################
    ! Simulator output info
    ! ######################################################################################
    integer, parameter :: nf_radar=17                    ! number of radar outputs
    integer, parameter :: nf_calipso=28                  ! number of calipso outputs
    integer, parameter :: nf_isccp=9                     ! number of isccp outputs
    integer, parameter :: nf_misr=1                      ! number of misr outputs
    integer, parameter :: nf_modis=20                    ! number of modis outputs
    
    ! Cloudsat outputs
    character(len=max_fieldname_len),dimension(nf_radar),parameter ::          &
         fname_radar = (/'CFAD_DBZE94_CS', 'CLD_CAL_NOTCS ', 'DBZE_CS       ', &
                         'CLDTOT_CALCS  ', 'CLDTOT_CS     ', 'CLDTOT_CS2    ', &
                         'CS_NOPRECIP   ', 'CS_RAINPOSS   ', 'CS_RAINPROB   ', &
                         'CS_RAINCERT   ', 'CS_SNOWPOSS   ', 'CS_SNOWCERT   ', &
                         'CS_MIXPOSS    ', 'CS_MIXCERT    ', 'CS_RAINHARD   ', &
                         'CS_UN         ', 'CS_PIA        '/)

    ! CALIPSO outputs
    character(len=max_fieldname_len),dimension(nf_calipso),parameter :: &
         fname_calipso=(/'CLDLOW_CAL     ','CLDMED_CAL     ','CLDHGH_CAL     ','CLDTOT_CAL     ','CLD_CAL        ',&
                         'RFL_PARASOL    ','CFAD_SR532_CAL ','ATB532_CAL     ','MOL532_CAL     ','CLD_CAL_LIQ    ',&
                         'CLD_CAL_ICE    ','CLD_CAL_UN     ','CLD_CAL_TMP    ','CLD_CAL_TMPLIQ ','CLD_CAL_TMPICE ',&
                         'CLD_CAL_TMPUN  ','CLDTOT_CAL_ICE ','CLDTOT_CAL_LIQ ','CLDTOT_CAL_UN  ','CLDHGH_CAL_ICE ',&
                         'CLDHGH_CAL_LIQ ','CLDHGH_CAL_UN  ','CLDMED_CAL_ICE ','CLDMED_CAL_LIQ ','CLDMED_CAL_UN  ',&
                         'CLDLOW_CAL_ICE ','CLDLOW_CAL_LIQ ','CLDLOW_CAL_UN  '/)
    ! ISCCP outputs
    character(len=max_fieldname_len),dimension(nf_isccp),parameter :: &
         fname_isccp=(/'FISCCP1_COSP    ','CLDTOT_ISCCP    ','MEANCLDALB_ISCCP',&
                       'MEANPTOP_ISCCP  ','TAU_ISCCP       ','CLDPTOP_ISCCP   ','MEANTAU_ISCCP   ',&
                       'MEANTB_ISCCP    ','MEANTBCLR_ISCCP '/)
    ! MISR outputs 
    character(len=max_fieldname_len),dimension(nf_misr),parameter :: &
         fname_misr=(/'CLD_MISR '/)
    ! MODIS outputs
    character(len=max_fieldname_len),dimension(nf_modis) :: &
         fname_modis=(/'CLTMODIS    ','CLWMODIS    ','CLIMODIS    ','CLHMODIS    ','CLMMODIS    ',&
                       'CLLMODIS    ','TAUTMODIS   ','TAUWMODIS   ','TAUIMODIS   ','TAUTLOGMODIS',&
                       'TAUWLOGMODIS','TAUILOGMODIS','REFFCLWMODIS','REFFCLIMODIS',&
                       'PCTMODIS    ','LWPMODIS    ','IWPMODIS    ','CLMODIS     ','CLRIMODIS   ',&
                       'CLRLMODIS   '/)
    
    logical :: run_radar(nf_radar,pcols)                 ! logical telling you if you should run radar simulator
    logical :: run_calipso(nf_calipso,pcols)             ! logical telling you if you should run calipso simulator
    logical :: run_isccp(nf_isccp,pcols)                 ! logical telling you if you should run isccp simulator
    logical :: run_misr(nf_misr,pcols)                   ! logical telling you if you should run misr simulator
    logical :: run_modis(nf_modis,pcols)                 ! logical telling you if you should run modis simulator
    
    ! CAM pointers to get variables from radiation interface (get from rad_cnst_get_gas)
    real(r8), pointer, dimension(:,:) :: q               ! specific humidity (kg/kg)
    real(r8), pointer, dimension(:,:) :: o3              ! Mass mixing ratio 03
    
    ! CAM pointers to get variables from the physics buffer
    real(r8), pointer, dimension(:,:) :: cld             ! cloud fraction, tca - total_cloud_amount (0-1)
    real(r8), pointer, dimension(:,:) :: concld          ! concld fraction, cca - convective_cloud_amount (0-1)
    real(r8), pointer, dimension(:,:) :: rel             ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei             ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffrain     ! rain effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: ls_reffsnow     ! snow effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffliq      ! convective cld liq effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: cv_reffice      ! convective cld ice effective drop size (microns)
    
    !! precip flux pointers
    real(r8), target, dimension(pcols,pverp) :: zero_ifc ! zero array for interface fields not in the pbuf
    real(r8), pointer, dimension(:,:) :: dp_flxprc       ! deep interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: dp_flxsnw       ! deep interface gbm flux_convective_cloud_snow (kg m^-2 s^-1) 
    real(r8), pointer, dimension(:,:) :: sh_flxprc       ! shallow interface gbm flux_convective_cloud_rain+snow (kg m^-2 s^-1) 
    real(r8), pointer, dimension(:,:) :: sh_flxsnw       ! shallow interface gbm flux_convective_cloud_snow (kg m^-2 s^-1)
    real(r8), pointer, dimension(:,:) :: ls_flxprc       ! stratiform interface gbm flux_cloud_rain+snow (kg m^-2 s^-1) 
    real(r8), pointer, dimension(:,:) :: ls_flxsnw       ! stratiform interface gbm flux_cloud_snow (kg m^-2 s^-1)
    
    !! grid box total cloud mixing ratio (large-scale + convective)
    real(r8), pointer, dimension(:,:) :: totg_liq       ! gbm total cloud liquid water (kg/kg)
    real(r8), pointer, dimension(:,:) :: totg_ice       ! gbm total cloud ice water (kg/kg)
    
    ! Output CAM variables
    ! Multiple "mdims" are collapsed because CAM history buffers only support one mdim.
    ! MIXED DIMS: ntau_cosp*nprs_cosp, CLOUDSAT_DBZE_BINS*nht_cosp, nsr_cosp*nht_cosp, nscol_cosp*nlay,
    !             ntau_cosp*nhtmisr_cosp
    real(r8) :: clisccp2(pcols,ntau_cosp,nprs_cosp)
    real(r8) :: cfad_dbze94(pcols,CLOUDSAT_DBZE_BINS,nht_cosp)
    real(r8) :: cfad_lidarsr532(pcols,nsr_cosp,nht_cosp)
    real(r8) :: dbze94(pcols,nscol_cosp,nlay)
    real(r8) :: atb532(pcols,nscol_cosp,nlay)
    real(r8) :: clMISR(pcols,ntau_cosp,nhtmisr_cosp)
    real(r8) :: frac_out(pcols,nscol_cosp,nlay)
    real(r8) :: cldtot_isccp(pcols)
    real(r8) :: meancldalb_isccp(pcols)
    real(r8) :: meanptop_isccp(pcols)
    real(r8) :: cldlow_cal(pcols)
    real(r8) :: cldmed_cal(pcols)
    real(r8) :: cldhgh_cal(pcols)
    real(r8) :: cldtot_cal(pcols)
    real(r8) :: cldtot_cal_ice(pcols)
    real(r8) :: cldtot_cal_liq(pcols)
    real(r8) :: cldtot_cal_un(pcols)
    real(r8) :: cldhgh_cal_ice(pcols)
    real(r8) :: cldhgh_cal_liq(pcols)
    real(r8) :: cldhgh_cal_un(pcols)
    real(r8) :: cldmed_cal_ice(pcols)
    real(r8) :: cldmed_cal_liq(pcols)
    real(r8) :: cldmed_cal_un(pcols)
    real(r8) :: cldlow_cal_ice(pcols)
    real(r8) :: cldlow_cal_liq(pcols)
    real(r8) :: cldlow_cal_un(pcols)
    real(r8) :: cld_cal(pcols,nht_cosp)
    real(r8) :: cld_cal_liq(pcols,nht_cosp)
    real(r8) :: cld_cal_ice(pcols,nht_cosp)
    real(r8) :: cld_cal_un(pcols,nht_cosp)
    real(r8) :: cld_cal_tmp(pcols,nht_cosp)
    real(r8) :: cld_cal_tmpliq(pcols,nht_cosp)
    real(r8) :: cld_cal_tmpice(pcols,nht_cosp)
    real(r8) :: cld_cal_tmpun(pcols,nht_cosp)
    real(r8) :: cfad_dbze94_cs(pcols,nht_cosp*CLOUDSAT_DBZE_BINS)
    real(r8) :: cfad_sr532_cal(pcols,nht_cosp*nsr_cosp)
    real(r8) :: tau_isccp(pcols,nscol_cosp)
    real(r8) :: cldptop_isccp(pcols,nscol_cosp)
    real(r8) :: meantau_isccp(pcols)
    real(r8) :: meantb_isccp(pcols)
    real(r8) :: meantbclr_isccp(pcols)
    real(r8) :: dbze_cs(pcols,nlay*nscol_cosp)
    real(r8) :: cldtot_calcs(pcols)
    real(r8) :: cldtot_cs(pcols)
    real(r8) :: cldtot_cs2(pcols)
    real(r8) :: ptcloudsatflag0(pcols)
    real(r8) :: ptcloudsatflag1(pcols)
    real(r8) :: ptcloudsatflag2(pcols)
    real(r8) :: ptcloudsatflag3(pcols)
    real(r8) :: ptcloudsatflag4(pcols)
    real(r8) :: ptcloudsatflag5(pcols)
    real(r8) :: ptcloudsatflag6(pcols)
    real(r8) :: ptcloudsatflag7(pcols)
    real(r8) :: ptcloudsatflag8(pcols)
    real(r8) :: ptcloudsatflag9(pcols)
    real(r8) :: cloudsatpia(pcols)
    real(r8) :: cld_cal_notcs(pcols,nht_cosp)
    real(r8) :: atb532_cal(pcols,nlay*nscol_cosp)
    real(r8) :: mol532_cal(pcols,nlay)
    real(r8) :: cld_misr(pcols,nhtmisr_cosp*ntau_cosp)
    real(r8) :: refl_parasol(pcols,nsza_cosp)
    real(r8) :: scops_out(pcols,nlay*nscol_cosp)
    real(r8) :: cltmodis(pcols)
    real(r8) :: clwmodis(pcols)
    real(r8) :: climodis(pcols)
    real(r8) :: clhmodis(pcols)
    real(r8) :: clmmodis(pcols)
    real(r8) :: cllmodis(pcols)
    real(r8) :: tautmodis(pcols)
    real(r8) :: tauwmodis(pcols)
    real(r8) :: tauimodis(pcols)
    real(r8) :: tautlogmodis(pcols)
    real(r8) :: tauwlogmodis(pcols)
    real(r8) :: tauilogmodis(pcols)
    real(r8) :: reffclwmodis(pcols)
    real(r8) :: reffclimodis(pcols)
    real(r8) :: pctmodis(pcols)
    real(r8) :: lwpmodis(pcols)
    real(r8) :: iwpmodis(pcols)
    real(r8) :: clmodis_cam(pcols,ntau_cosp_modis*nprs_cosp)
    real(r8) :: clmodis(pcols,ntau_cosp_modis,nprs_cosp)
    real(r8) :: clrimodis_cam(pcols,ntau_cosp*numMODISReffIceBins)
    real(r8) :: clrimodis(pcols,ntau_cosp,numMODISReffIceBins)
    real(r8) :: clrlmodis_cam(pcols,ntau_cosp*numMODISReffLiqBins)
    real(r8) :: clrlmodis(pcols,ntau_cosp,numMODISReffLiqBins)
    real(r8), dimension(pcols,nlay*nscol_cosp) :: &
       tau067_out, emis11_out, fracliq_out, asym34_out, ssa34_out

    type(interp_type)  :: interp_wgts
    integer, parameter :: extrap_method = 1   ! sets extrapolation method to boundary value (1)
    
    ! COSPv2 stuff
    character(len=256),dimension(100) :: cosp_status
    integer :: nerror

    integer :: istat
    character(len=*), parameter :: sub = 'cospsimulator_intr_run'
    !--------------------------------------------------------------------------------------

    call t_startf("init_and_stuff")
    ! ######################################################################################
    ! Initialization
    ! ######################################################################################

    lchnk = state%lchnk    ! chunk ID
    ncol  = state%ncol     ! number of columns in the chunk
    Npoints = ncol         ! number of COSP gridpoints
    
    zero_ifc = 0._r8

    ! Initialize CAM variables as R_UNDEF, important for history files because it will exclude these from averages
    ! initialize over all pcols, not just ncol.  missing values needed in chunks where ncol<pcols
    clisccp2(1:pcols,1:ntau_cosp,1:nprs_cosp)     = R_UNDEF
    cfad_dbze94(1:pcols,1:CLOUDSAT_DBZE_BINS,1:nht_cosp)  = R_UNDEF
    cfad_lidarsr532(1:pcols,1:nsr_cosp,1:nht_cosp)= R_UNDEF
    dbze94(1:pcols,1:nscol_cosp,1:nlay)           = R_UNDEF
    atb532(1:pcols,1:nscol_cosp,1:nlay)           = R_UNDEF
    clMISR(1:pcols,ntau_cosp,1:nhtmisr_cosp)      = R_UNDEF
    frac_out(1:pcols,1:nscol_cosp,1:nlay)         = R_UNDEF
    
    ! (all CAM output variables. including collapsed variables)
    cldtot_isccp(1:pcols)                            = R_UNDEF
    meancldalb_isccp(1:pcols)                        = R_UNDEF
    meanptop_isccp(1:pcols)                          = R_UNDEF
    cldlow_cal(1:pcols)                              = R_UNDEF
    cldmed_cal(1:pcols)                              = R_UNDEF
    cldhgh_cal(1:pcols)                              = R_UNDEF
    cldtot_cal(1:pcols)                              = R_UNDEF
    cldtot_cal_ice(1:pcols)                          = R_UNDEF !+cosp1.4
    cldtot_cal_liq(1:pcols)                          = R_UNDEF
    cldtot_cal_un(1:pcols)                           = R_UNDEF
    cldhgh_cal_ice(1:pcols)                          = R_UNDEF
    cldhgh_cal_liq(1:pcols)                          = R_UNDEF
    cldhgh_cal_un(1:pcols)                           = R_UNDEF
    cldmed_cal_ice(1:pcols)                          = R_UNDEF
    cldmed_cal_liq(1:pcols)                          = R_UNDEF
    cldmed_cal_un(1:pcols)                           = R_UNDEF
    cldlow_cal_liq(1:pcols)                          = R_UNDEF
    cldlow_cal_ice(1:pcols)                          = R_UNDEF
    cldlow_cal_un(1:pcols)                           = R_UNDEF !+cosp1.4
    cld_cal(1:pcols,1:nht_cosp)                      = R_UNDEF
    cld_cal_liq(1:pcols,1:nht_cosp)                  = R_UNDEF !+cosp1.4
    cld_cal_ice(1:pcols,1:nht_cosp)                  = R_UNDEF
    cld_cal_un(1:pcols,1:nht_cosp)                   = R_UNDEF
    cld_cal_tmp(1:pcols,1:nht_cosp)                  = R_UNDEF
    cld_cal_tmpliq(1:pcols,1:nht_cosp)               = R_UNDEF
    cld_cal_tmpice(1:pcols,1:nht_cosp)               = R_UNDEF
    cld_cal_tmpun(1:pcols,1:nht_cosp)                = R_UNDEF
    cfad_dbze94_cs(1:pcols,1:nht_cosp*CLOUDSAT_DBZE_BINS)    = R_UNDEF
    cfad_sr532_cal(1:pcols,1:nht_cosp*nsr_cosp)      = R_UNDEF
    tau_isccp(1:pcols,1:nscol_cosp)                  = R_UNDEF
    cldptop_isccp(1:pcols,1:nscol_cosp)              = R_UNDEF
    meantau_isccp(1:pcols)                           = R_UNDEF
    meantb_isccp(1:pcols)                            = R_UNDEF
    meantbclr_isccp(1:pcols)                         = R_UNDEF     
    dbze_cs(1:pcols,1:nlay*nscol_cosp)               = R_UNDEF
    ptcloudsatflag0(1:pcols)                         = R_UNDEF 
    ptcloudsatflag1(1:pcols)                         = R_UNDEF 
    ptcloudsatflag2(1:pcols)                         = R_UNDEF 
    ptcloudsatflag3(1:pcols)                         = R_UNDEF 
    ptcloudsatflag4(1:pcols)                         = R_UNDEF 
    ptcloudsatflag5(1:pcols)                         = R_UNDEF 
    ptcloudsatflag6(1:pcols)                         = R_UNDEF 
    ptcloudsatflag7(1:pcols)                         = R_UNDEF 
    ptcloudsatflag8(1:pcols)                         = R_UNDEF 
    ptcloudsatflag9(1:pcols)                         = R_UNDEF 
    cloudsatpia(1:pcols)                             = R_UNDEF 
    cldtot_calcs(1:pcols)                            = R_UNDEF
    cldtot_cs(1:pcols)                               = R_UNDEF
    cldtot_cs2(1:pcols)                              = R_UNDEF
    cld_cal_notcs(1:pcols,1:nht_cosp)                = R_UNDEF
    atb532_cal(1:pcols,1:nlay*nscol_cosp)            = R_UNDEF
    mol532_cal(1:pcols,1:nlay)                       = R_UNDEF
    cld_misr(1:pcols,1:nhtmisr_cosp*ntau_cosp)       = R_UNDEF
    refl_parasol(1:pcols,1:nsza_cosp)                = R_UNDEF
    scops_out(1:pcols,1:nlay*nscol_cosp)             = R_UNDEF
    cltmodis(1:pcols)                                = R_UNDEF
    clwmodis(1:pcols)                                = R_UNDEF
    climodis(1:pcols)                                = R_UNDEF
    clhmodis(1:pcols)                                = R_UNDEF
    clmmodis(1:pcols)                                = R_UNDEF
    cllmodis(1:pcols)                                = R_UNDEF
    tautmodis(1:pcols)                               = R_UNDEF
    tauwmodis(1:pcols)                               = R_UNDEF
    tauimodis(1:pcols)                               = R_UNDEF
    tautlogmodis(1:pcols)                            = R_UNDEF
    tauwlogmodis(1:pcols)                            = R_UNDEF
    tauilogmodis(1:pcols)                            = R_UNDEF
    reffclwmodis(1:pcols)                            = R_UNDEF
    reffclimodis(1:pcols)                            = R_UNDEF
    pctmodis(1:pcols)                                = R_UNDEF
    lwpmodis(1:pcols)                                = R_UNDEF
    iwpmodis(1:pcols)                                = R_UNDEF
    clmodis_cam(1:pcols,1:ntau_cosp_modis*nprs_cosp) = R_UNDEF
    clmodis(1:pcols,1:ntau_cosp_modis,1:nprs_cosp)   = R_UNDEF
    clrimodis_cam(1:pcols,1:ntau_cosp_modis*numMODISReffIceBins) = R_UNDEF ! +cosp2
    clrimodis(1:pcols,1:ntau_cosp_modis,1:numMODISReffIceBins)   = R_UNDEF ! +cosp2
    clrlmodis_cam(1:pcols,1:ntau_cosp_modis*numMODISReffLiqBins) = R_UNDEF ! +cosp2
    clrlmodis(1:pcols,1:ntau_cosp_modis,1:numMODISReffLiqBins)   = R_UNDEF ! +cosp2
    tau067_out(1:pcols,1:nlay*nscol_cosp)            = R_UNDEF ! +cosp2
    emis11_out(1:pcols,1:nlay*nscol_cosp)            = R_UNDEF ! +cosp2
    asym34_out(1:pcols,1:nlay*nscol_cosp)            = R_UNDEF ! +cosp2
    ssa34_out(1:pcols,1:nlay*nscol_cosp)             = R_UNDEF ! +cosp2
    fracLiq_out(1:pcols,1:nlay*nscol_cosp)           = R_UNDEF ! +cosp2

    ! ######################################################################################
    ! DECIDE WHICH COLUMNS YOU ARE GOING TO RUN COSP ON....
    ! ######################################################################################
    
    !! run_cosp is set for each column in each chunk in the first timestep of the run
    !! hist_fld_col_active in cam_history.F90 is used to decide if you need to run cosp.
    if (first_run_cosp(lchnk)) then
       !! initalize to run logicals as false
       run_cosp(1:ncol,lchnk)=.false.
       run_radar(1:nf_radar,1:ncol)=.false.
       run_calipso(1:nf_calipso,1:ncol)=.false.
       run_isccp(1:nf_isccp,1:ncol)=.false.
       run_misr(1:nf_misr,1:ncol)=.false.
       run_modis(1:nf_modis,1:ncol)=.false.
       
       if (lradar_sim) then
          do i=1,nf_radar
             run_radar(i,1:pcols)=hist_fld_col_active(fname_radar(i),lchnk,pcols)
          end do
       end if
       if (llidar_sim) then
          do i=1,nf_calipso
             run_calipso(i,1:pcols)=hist_fld_col_active(fname_calipso(i),lchnk,pcols)
          end do
       end if
       if (lisccp_sim) then
          do i=1,nf_isccp
             run_isccp(i,1:pcols)=hist_fld_col_active(fname_isccp(i),lchnk,pcols)
          end do
       end if
       if (lmisr_sim) then
          do i=1,nf_misr
             run_misr(i,1:pcols)=hist_fld_col_active(fname_misr(i),lchnk,pcols)
          end do
       end if
       if (lmodis_sim) then
          do i=1,nf_modis
             run_modis(i,1:pcols)=hist_fld_col_active(fname_modis(i),lchnk,pcols)
          end do
       end if
       
       do i=1,ncol
          if ((any(run_radar(:,i))) .or. (any(run_calipso(:,i))) .or. (any(run_isccp(:,i))) &
               .or. (any(run_misr(:,i))) .or. (any(run_modis(:,i)))) then
             run_cosp(i,lchnk)=.true.
          end if
       end do
       
       first_run_cosp(lchnk)=.false.
    endif
    
    ! ######################################################################################
    ! GET CAM GEOPHYSICAL VARIABLES NEEDED FOR COSP INPUT
    ! ######################################################################################
    ! state variables (prognostic variables, see physics_types.F90)
    ! state%lat   ! lat (radians) 
    ! state%lon   ! lon (radians) 
    ! state%t     ! temperature (K)
    ! state%ps    ! surface pressure (Pa)
    ! state%pint  ! p - p_in_full_levels (Pa)
    ! state%pmid  ! ph - p_in_half_levels (Pa)
    ! state%zm    ! geopotential height above surface at midpoints (m), pver
    ! state%zi    ! geopotential height above surface at interfaces (m), pverp
    ! state%phis  ! surface geopotential (m2/s2)
    ! NOTE: The state variables state%q(:,:,ixcldliq)/state%q(:,:,ixcldice) are grid-box
    ! quantities for the stratiform clouds only.  stratiform water * stratiform cloud fraction
    ! state%q(:,:,ixcldliq) ! for CAM4: cldliq = stratiform incld water content * total cloud fraction
    ! state%q(:,:,ixcldice) ! for CAM4: cldice = stratiform incld ice content * total cloud fraction
    
    ! advected constiutent indices
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    
    ! radiative constituents (prognostic or data)
    call rad_cnst_get_gas(0,'H2O', state, pbuf,  q)                     
    call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3)
    
    ! fields from physics buffer
    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, concld_idx, concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, rel_idx, rel  )
    call pbuf_get_field(pbuf, rei_idx, rei)
    call pbuf_get_field(pbuf, lsreffrain_idx, ls_reffrain  )
    call pbuf_get_field(pbuf, lsreffsnow_idx, ls_reffsnow  )
    call pbuf_get_field(pbuf, cvreffliq_idx,  cv_reffliq   )
    call pbuf_get_field(pbuf, cvreffice_idx,  cv_reffice   )
    
    !! grid box total cloud mixing ratios
    call pbuf_get_field(pbuf, gb_totcldliqmr_idx, totg_liq)
    call pbuf_get_field(pbuf, gb_totcldicemr_idx, totg_ice)
    
    !! precipitation fluxes
    call pbuf_get_field(pbuf, dpflxprc_idx, dp_flxprc  )
    call pbuf_get_field(pbuf, dpflxsnw_idx, dp_flxsnw  )
    if (shflxprc_idx > 0) then
       call pbuf_get_field(pbuf, shflxprc_idx, sh_flxprc  )
    else
       sh_flxprc => zero_ifc
    end if
    if (shflxsnw_idx > 0) then
       call pbuf_get_field(pbuf, shflxsnw_idx, sh_flxsnw  )
    else
       sh_flxsnw => zero_ifc
    end if
    call pbuf_get_field(pbuf, lsflxprc_idx, ls_flxprc  )
    call pbuf_get_field(pbuf, lsflxsnw_idx, ls_flxsnw  )
   
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! CALCULATE COSP INPUT VARIABLES FROM CAM VARIABLES, done for all columns within chunk
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! These arrays are dimensioned to only include active columns (ncol), and the number
    ! of layers (nlay) and layer interfaces (nlayp) operated on by COSP.
    allocate( &
       zmid(ncol,nlay), &
       zint(ncol,nlayp), &
       surf_hgt(ncol),   &
       landmask(ncol),   &
       mr_ccliq(ncol,nlay), &
       mr_ccice(ncol,nlay), &
       mr_lsliq(ncol,nlay), &
       mr_lsice(ncol,nlay), &
       rain_cv(ncol,nlayp), &
       snow_cv(ncol,nlayp), &
       rain_cv_interp(ncol,nlay), &
       snow_cv_interp(ncol,nlay), &
       rain_ls_interp(ncol,nlay), &
       snow_ls_interp(ncol,nlay), &
       grpl_ls_interp(ncol,nlay), &
       reff_cosp(ncol,nlay,nhydro), &
       dtau_s(ncol,nlay), &
       dtau_c(ncol,nlay), &
       dtau_s_snow(ncol,nlay), &
       dem_s(ncol,nlay), &
       dem_c(ncol,nlay), &
       dem_s_snow(ncol,nlay), stat=istat)
    call handle_allocate_error(istat, sub, 'zmid,..,dem_s_snow')
    
    ! add surface height (surface geopotential/gravity) to convert CAM heights based on
    ! geopotential above surface into height above sea level
    surf_hgt = state%phis(:ncol)*inverse_gravit
    do k = 1, nlay
       zmid(:,k) = state%zm(:ncol,ktop+k-1) + surf_hgt
       zint(:,k) = state%zi(:ncol,ktop+k-1) + surf_hgt
    end do
    zint(:,nlayp) = surf_hgt

    landmask = 0._r8
    do i = 1, ncol
       if (cam_in%landfrac(i) > 0.01_r8) landmask(i)= 1
    end do
    
    ! Add together deep and shallow convection precipitation fluxes.
    ! Note: sh_flxprc and dp_flxprc variables are rain+snow
    rain_cv = (sh_flxprc(:ncol,ktop:pverp) - sh_flxsnw(:ncol,ktop:pverp)) + &
              (dp_flxprc(:ncol,ktop:pverp) - dp_flxsnw(:ncol,ktop:pverp))
    snow_cv = sh_flxsnw(:ncol,ktop:pverp) + dp_flxsnw(:ncol,ktop:pverp)
    
    ! interpolate interface precip fluxes to mid points
    do i = 1, ncol
       ! find weights
       call lininterp_init(state%zi(i,ktop:pverp), nlayp, state%zm(i,ktop:pver), nlay, &
                           extrap_method, interp_wgts)
       ! interpolate  lininterp(arrin, nin, arrout, nout, interp_wgts)
       call lininterp(rain_cv(i,:), nlayp, rain_cv_interp(i,:), nlay, interp_wgts)
       call lininterp(snow_cv(i,:), nlayp, snow_cv_interp(i,:), nlay, interp_wgts)
       call lininterp(ls_flxprc(i,ktop:pverp), nlayp, rain_ls_interp(i,:), nlay, interp_wgts)
       call lininterp(ls_flxsnw(i,ktop:pverp), nlayp, snow_ls_interp(i,:), nlay, interp_wgts)
       call lininterp_finish(interp_wgts)
       !! ls_flxprc is for rain+snow, find rain_ls_interp by subtracting off snow_ls_interp
       rain_ls_interp(i,:) = rain_ls_interp(i,:) - snow_ls_interp(i,:)
    end do

    !! Make sure interpolated values are not less than 0
    do k = 1, nlay
       do i = 1, ncol
          if (rain_ls_interp(i,k) < 0._r8) then
             rain_ls_interp(i,k) = 0._r8
          end if
          if (snow_ls_interp(i,k) < 0._r8) then
             snow_ls_interp(i,k) = 0._r8
          end if
          if (rain_cv_interp(i,k) < 0._r8) then
             rain_cv_interp(i,k) = 0._r8
          end if
          if (snow_cv_interp(i,k) < 0._r8) then
             snow_cv_interp(i,k) = 0._r8
          end if
       end do
    end do
    
    grpl_ls_interp = 0._r8
    
    ! subroutine subsample_and_optics provides separate arguments to pass
    ! the large scale and convective cloud condensate.  Below the grid box
    ! total cloud water mixing ratios are passed in the arrays for the
    ! large scale contributions and the arrays for the convective
    ! contributions are set to zero.  This is consistent with the treatment
    ! of cloud water by the radiation code.
    mr_ccliq = 0._r8
    mr_ccice = 0._r8
    mr_lsliq = 0._r8
    mr_lsice = 0._r8
    do k = 1, nlay
       kk = ktop + k -1
       do i = 1, ncol
          if (cld(i,k) > 0._r8) then
             mr_lsliq(i,k) = totg_liq(i,kk)
             mr_lsice(i,k) = totg_ice(i,kk)
          end if
       end do
    end do
    
    !! The specification of reff_cosp now follows e-mail discussion with Yuying in January 2011.
    !! The values from the physics buffer are in microns... convert to meters for COSP.
    reff_cosp(:,:,I_LSCLIQ) = rel(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_LSCICE) = rei(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_LSRAIN) = ls_reffrain(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_LSSNOW) = ls_reffsnow(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_CVCLIQ) = cv_reffliq(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_CVCICE) = cv_reffice(:ncol,ktop:pver)*1.e-6_r8
    reff_cosp(:,:,I_CVRAIN) = ls_reffrain(:ncol,ktop:pver)*1.e-6_r8  !! same as stratiform per Andrew
    reff_cosp(:,:,I_CVSNOW) = ls_reffsnow(:ncol,ktop:pver)*1.e-6_r8  !! same as stratiform per Andrew
    reff_cosp(:,:,I_LSGRPL) = 0._r8                                  !! using radar default reff
 
    ! assign optical depths and emissivities
    ! CAM4 assumes same radiative properties for stratiform and convective clouds, 
    !   (see ISCCP_CLOUD_TYPES subroutine call in cloudsimulator.F90)
    !   Assume CAM5 is doing the same thing based on the ISCCP simulator calls within RRTM's radiation.F90
    ! COSP wants in-cloud values.  CAM5 values cld_swtau are in-cloud.
    ! snow_tau_in and snow_emis_in are passed without modification to COSP
    dtau_s      = cld_swtau_in(:ncol,ktop:pver)
    dtau_c      = cld_swtau_in(:ncol,ktop:pver)
    dtau_s_snow = snow_tau_in(:ncol,ktop:pver)
    dem_s       = emis(:ncol,ktop:pver)
    dem_c       = emis(:ncol,ktop:pver)
    dem_s_snow  = snow_emis_in(:ncol,ktop:pver)

    ! ######################################################################################
    ! Compute sunlit flag. If cosp_runall=.true., then run on all points.
    ! ######################################################################################
    cam_sunlit(:) = 0
    if (cosp_runall) then
       cam_sunlit(:) = 1
       nSunLit   = ncol
    else
       nSunLit   = 0
       do i=1,ncol
          if ((coszrs(i) > 0.0_r8) .and. (run_cosp(i,lchnk))) then
             cam_sunlit(i) = 1
             nSunLit   = nSunLit+1
          endif
       enddo
    endif
    call t_stopf("init_and_stuff")

    ! ######################################################################################
    ! Construct COSP output derived type.
    ! ######################################################################################
    call t_startf("construct_cosp_outputs")
    call construct_cosp_outputs(ncol, nscol_cosp, nlay, Nlvgrid, cospOUT)
    call t_stopf("construct_cosp_outputs")
    
    ! ######################################################################################
    ! Construct and populate COSP input types
    ! ######################################################################################
    ! Model state
    call t_startf("construct_cospstateIN")

    call construct_cospstateIN(ncol, nlay, 0, cospstateIN)      

    ! convert to degrees.  Lat in range [-90,..,90], Lon in range [0,..,360]
    cospstateIN%lat             = state%lat(:ncol)*rad2deg
    cospstateIN%lon             = state%lon(:ncol)*rad2deg
    cospstateIN%at              = state%t(:ncol,ktop:pver)
    cospstateIN%qv              = q(:ncol,ktop:pver)
    cospstateIN%o3              = o3(:ncol,ktop:pver)
    cospstateIN%sunlit          = cam_sunlit(:ncol)
    cospstateIN%skt             = cam_in%ts(:ncol)
    cospstateIN%land            = landmask
    cospstateIN%pfull           = state%pmid(:ncol,ktop:pver)
    cospstateIN%phalf           = state%pint(:ncol,ktop:pverp)
    cospstateIN%hgt_matrix      = zmid
    cospstateIN%hgt_matrix_half = zint
    cospstateIN%surfelev        = surf_hgt
    call t_stopf("construct_cospstateIN")

    ! Optical inputs
    call t_startf("construct_cospIN")
    call construct_cospIN(ncol, nscol_cosp, nlay, cospIN)
    cospIN%emsfc_lw = emsfc_lw
    if (lradar_sim) cospIN%rcfg_cloudsat = rcfg_cs(lchnk)
    call t_stopf("construct_cospIN")

    call t_startf("subsample_and_optics")
    ! The arrays passed here contain only active columns and the limited vertical
    ! domain operated on by COSP.  Unsubscripted array arguments have already been
    ! allocated to the correct size.  Arrays the size of a CAM chunk (pcol,pver)
    ! need to pass the correct section (:ncol,ktop:pver).
    call subsample_and_optics( &
       ncol, nlay, nscol_cosp, nhydro, overlap, &
       lidar_ice_type, sd_cs(lchnk), &
       cld(:ncol,ktop:pver), concld(:ncol,ktop:pver), &
       rain_ls_interp, snow_ls_interp, grpl_ls_interp, rain_cv_interp, &
       snow_cv_interp, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, &
       reff_cosp, dtau_c, dtau_s ,dem_c, dem_s, dtau_s_snow, &
       dem_s_snow, state%ps(:ncol), cospstateIN, cospIN)
    call t_stopf("subsample_and_optics")
    
    ! ######################################################################################
    ! Call COSP
    ! ######################################################################################
    call t_startf("cosp_simulator")
    cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx=1, stop_idx=ncol,debug=.false.)

    ! Check status flags
    nerror = 0
    do i = 1, ubound(cosp_status, 1)
       if (len_trim(cosp_status(i)) > 0) then
          write(iulog,*) "cosp_simulator: ERROR: "//trim(cosp_status(i))
          nerror = nerror + 1
       end if
    end do
    if (nerror > 0) then
       call endrun('cospsimulator_intr_run: error return from cosp_simulator')
    end if
    call t_stopf("cosp_simulator")
  
    ! ######################################################################################
    ! Write COSP inputs to output file for offline use.
    ! ######################################################################################
    call t_startf("cosp_histfile_aux")
    if (cosp_histfile_aux) then
       ! 1D outputs
       call outfld('PS_COSP',        state%ps(1:ncol),             ncol,lchnk)
       call outfld('TS_COSP',        cospstateIN%skt,              ncol,lchnk)
       
       ! 2D outputs
       call outfld('P_COSP',         cospstateIN%pfull,            ncol,lchnk)
       call outfld('PH_COSP',        cospstateIN%phalf,            ncol,lchnk)
       call outfld('ZLEV_COSP',      cospstateIN%hgt_matrix,       ncol,lchnk)
       call outfld('ZLEV_HALF_COSP', cospstateIN%hgt_matrix_half,  ncol,lchnk)
       call outfld('T_COSP',         cospstateIN%at,               ncol,lchnk)
       call outfld('Q_COSP',         cospstateIN%qv,               ncol,lchnk)

       ! 3D outputs, but first compress to 2D
       do i=1,ncol
          do ihml=1,nlay
             do isc=1,nscol_cosp
                ihsc = (ihml-1)*nscol_cosp+isc                 
                tau067_out(i,ihsc)  = cospIN%tau_067(i,isc,ihml)
                emis11_out(i,ihsc)  = cospIN%emiss_11(i,isc,ihml)
                ssa34_out(i,ihsc)   = cospIN%ss_alb(i,isc,ihml)
                asym34_out(i,ihsc)  = cospIN%asym(i,isc,ihml)
                fracLiq_out(i,ihsc) = cospIN%fracLiq(i,isc,ihml)
             end do
          end do
       end do
       call outfld('TAU_067',      tau067_out, pcols,lchnk)
       call outfld('EMISS_11',     emis11_out, pcols,lchnk)
       call outfld('MODIS_asym',   asym34_out, pcols,lchnk)
       call outfld('MODIS_ssa',    ssa34_out,  pcols,lchnk)
       call outfld('MODIS_fracliq',fracLiq_out,pcols,lchnk)
    end if
    call t_stopf("cosp_histfile_aux")

    ! ######################################################################################
    ! Set dark-scenes to fill value. Only done for passive simulators and when cosp_runall=F
    ! ######################################################################################
    call t_startf("sunlit_passive")
    if (.not. cosp_runall) then
       ! ISCCP simulator
       if (lisccp_sim) then
          ! 1D
          where(cam_sunlit(1:ncol) .eq. 0)
             cospOUT%isccp_totalcldarea(1:ncol)  = R_UNDEF
             cospOUT%isccp_meanptop(1:ncol)      = R_UNDEF
             cospOUT%isccp_meantaucld(1:ncol)    = R_UNDEF
             cospOUT%isccp_meanalbedocld(1:ncol) = R_UNDEF
             cospOUT%isccp_meantb(1:ncol)        = R_UNDEF
             cospOUT%isccp_meantbclr(1:ncol)     = R_UNDEF
          end where
          ! 2D
          do i=1,nscol_cosp
             where (cam_sunlit(1:ncol) .eq. 0)
                cospOUT%isccp_boxtau(1:ncol,i)  = R_UNDEF
                cospOUT%isccp_boxptop(1:ncol,i) = R_UNDEF
             end where
          enddo
          ! 3D
          do i=1,nprs_cosp
             do k=1,ntau_cosp
                where(cam_sunlit(1:ncol) .eq. 0)
                   cospOUT%isccp_fq(1:ncol,k,i) = R_UNDEF
                end where
             end do
          end do
       endif

       ! MISR simulator
       if (lmisr_sim) then
          do i=1,nhtmisr_cosp
             do k=1,ntau_cosp
                where(cam_sunlit(1:ncol) .eq. 0)
                   cospOUT%misr_fq(1:ncol,k,i) = R_UNDEF
                end where
             end do
          end do
       end if

       ! MODIS simulator
       if (lmodis_sim) then
          ! 1D
          where(cam_sunlit(1:ncol) .eq. 0)
             cospOUT%modis_Cloud_Fraction_Total_Mean(1:ncol)       = R_UNDEF
             cospOUT%modis_Cloud_Fraction_Water_Mean(1:ncol)       = R_UNDEF
             cospOUT%modis_Cloud_Fraction_Ice_Mean(1:ncol)         = R_UNDEF
             cospOUT%modis_Cloud_Fraction_High_Mean(1:ncol)        = R_UNDEF
             cospOUT%modis_Cloud_Fraction_Mid_Mean(1:ncol)         = R_UNDEF
             cospOUT%modis_Cloud_Fraction_Low_Mean(1:ncol)         = R_UNDEF
             cospOUT%modis_Optical_Thickness_Total_Mean(1:ncol)    = R_UNDEF
             cospOUT%modis_Optical_Thickness_Water_Mean(1:ncol)    = R_UNDEF
             cospOUT%modis_Optical_Thickness_Ice_Mean(1:ncol)      = R_UNDEF
             cospOUT%modis_Optical_Thickness_Total_LogMean(1:ncol) = R_UNDEF
             cospOUT%modis_Optical_Thickness_Water_LogMean(1:ncol) = R_UNDEF
             cospOUT%modis_Optical_Thickness_Ice_LogMean(1:ncol)   = R_UNDEF
             cospOUT%modis_Cloud_Particle_Size_Water_Mean(1:ncol)  = R_UNDEF
             cospOUT%modis_Cloud_Particle_Size_Ice_Mean(1:ncol)    = R_UNDEF
             cospOUT%modis_Cloud_Top_Pressure_Total_Mean(1:ncol)   = R_UNDEF
             cospOUT%modis_Liquid_Water_Path_Mean(1:ncol)          = R_UNDEF
             cospOUT%modis_Ice_Water_Path_Mean(1:ncol)             = R_UNDEF
          endwhere
          ! 3D
          do i=1,ntau_cosp_modis
             do k=1,nprs_cosp
                where(cam_sunlit(1:ncol) .eq. 0)
                   cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(1:ncol,i,k) = R_UNDEF 
                end where
             enddo
             do k=1,numMODISReffIceBins
                where(cam_sunlit(1:ncol) .eq. 0)
                   cospOUT%modis_Optical_Thickness_vs_ReffICE(1:ncol,i,k) = R_UNDEF
                end where
             end do
             do k=1,numMODISReffLiqBins
                where(cam_sunlit(1:ncol) .eq. 0)
                   cospOUT%modis_Optical_Thickness_vs_ReffLIQ(1:ncol,i,k) = R_UNDEF
                end where
             enddo
          enddo
       end if
    end if
    call t_stopf("sunlit_passive")

    ! ######################################################################################
    ! Copy COSP outputs to CAM fields.
    ! ######################################################################################
    call t_startf("output_copying")
    if (allocated(cospIN%frac_out)) &
         frac_out(1:ncol,1:nscol_cosp,1:nlay) = cospIN%frac_out
    
    ! Cloudsat
    if (lradar_sim) then 
       cfad_dbze94(1:ncol,1:CLOUDSAT_DBZE_BINS,1:nht_cosp) = cospOUT%cloudsat_cfad_ze
       dbze94(1:ncol,1:nscol_cosp,1:nlay)    = cospOUT%cloudsat_Ze_tot
       cldtot_cs(1:ncol)  = 0._r8
       cldtot_cs2(1:ncol) = 0._r8
       ! *NOTE* These two fields are joint-simulator products, but in CAM they are controlled
       !        by the radar simulator control.
       cldtot_calcs(1:ncol) = cospOUT%radar_lidar_tcc
       cld_cal_notcs(1:ncol,1:nht_cosp) = cospOUT%lidar_only_freq_cloud

       ! Cloudsat near-surface precipitation diagnostics
       ptcloudsatflag0(1:ncol) = cospOUT%cloudsat_precip_cover(:,1)
       ptcloudsatflag1(1:ncol) = cospOUT%cloudsat_precip_cover(:,2)
       ptcloudsatflag2(1:ncol) = cospOUT%cloudsat_precip_cover(:,3)
       ptcloudsatflag3(1:ncol) = cospOUT%cloudsat_precip_cover(:,4)
       ptcloudsatflag4(1:ncol) = cospOUT%cloudsat_precip_cover(:,5)
       ptcloudsatflag5(1:ncol) = cospOUT%cloudsat_precip_cover(:,6)
       ptcloudsatflag6(1:ncol) = cospOUT%cloudsat_precip_cover(:,7)
       ptcloudsatflag7(1:ncol) = cospOUT%cloudsat_precip_cover(:,8)
       ptcloudsatflag8(1:ncol) = cospOUT%cloudsat_precip_cover(:,9)
       ptcloudsatflag9(1:ncol) = cospOUT%cloudsat_precip_cover(:,10)
       cloudsatpia(1:ncol)     = cospOUT%cloudsat_pia
       
    endif
    
    ! CALIPSO
    if (llidar_sim) then
       cldlow_cal(1:ncol)                = cospOUT%calipso_cldlayer(:,1)
       cldmed_cal(1:ncol)                = cospOUT%calipso_cldlayer(:,2)
       cldhgh_cal(1:ncol)                = cospOUT%calipso_cldlayer(:,3)
       cldtot_cal(1:ncol)                = cospOUT%calipso_cldlayer(:,4)
       cldlow_cal_ice(1:ncol)            = cospOUT%calipso_cldlayerphase(:,1,1)
       cldmed_cal_ice(1:ncol)            = cospOUT%calipso_cldlayerphase(:,2,1)
       cldhgh_cal_ice(1:ncol)            = cospOUT%calipso_cldlayerphase(:,3,1)
       cldtot_cal_ice(1:ncol)            = cospOUT%calipso_cldlayerphase(:,4,1)
       cldlow_cal_liq(1:ncol)            = cospOUT%calipso_cldlayerphase(:,1,2)
       cldmed_cal_liq(1:ncol)            = cospOUT%calipso_cldlayerphase(:,2,2)
       cldhgh_cal_liq(1:ncol)            = cospOUT%calipso_cldlayerphase(:,3,2)
       cldtot_cal_liq(1:ncol)            = cospOUT%calipso_cldlayerphase(:,4,2)
       cldlow_cal_un(1:ncol)             = cospOUT%calipso_cldlayerphase(:,1,3)
       cldmed_cal_un(1:ncol)             = cospOUT%calipso_cldlayerphase(:,2,3)
       cldhgh_cal_un(1:ncol)             = cospOUT%calipso_cldlayerphase(:,3,3)
       cldtot_cal_un(1:ncol)             = cospOUT%calipso_cldlayerphase(:,4,3)
       cld_cal_ice(1:ncol,1:nht_cosp)    = cospOUT%calipso_lidarcldphase(:,:,1)
       cld_cal_liq(1:ncol,1:nht_cosp)    = cospOUT%calipso_lidarcldphase(:,:,2)
       cld_cal_un(1:ncol,1:nht_cosp)     = cospOUT%calipso_lidarcldphase(:,:,3)
       cld_cal_tmp(1:ncol,1:nht_cosp)    = cospOUT%calipso_lidarcldtmp(:,:,1)
       cld_cal_tmpliq(1:ncol,1:nht_cosp) = cospOUT%calipso_lidarcldtmp(:,:,2)
       cld_cal_tmpice(1:ncol,1:nht_cosp) = cospOUT%calipso_lidarcldtmp(:,:,3)
       cld_cal_tmpun(1:ncol,1:nht_cosp)  = cospOUT%calipso_lidarcldtmp(:,:,4)
       cld_cal(1:ncol,1:nht_cosp)        = cospOUT%calipso_lidarcld(:,1:nht_cosp)
       mol532_cal(1:ncol,1:nlay)         = cospOUT%calipso_beta_mol
       atb532(1:ncol,1:nscol_cosp,1:nlay)= cospOUT%calipso_beta_tot
       cfad_lidarsr532(1:ncol,1:nsr_cosp,1:nht_cosp) = cospOUT%calipso_cfad_sr(:,:,:)
       refl_parasol(1:ncol,1:nsza_cosp)  = cospOUT%parasolGrid_refl
    endif
    
    ! ISCCP
    if (lisccp_sim) then
       clisccp2(1:ncol,1:ntau_cosp,1:nprs_cosp) = cospOUT%isccp_fq
       tau_isccp(1:ncol,1:nscol_cosp)           = cospOUT%isccp_boxtau
       cldptop_isccp(1:ncol,1:nscol_cosp)       = cospOUT%isccp_boxptop
       cldtot_isccp(1:ncol)                     = cospOUT%isccp_totalcldarea
       meanptop_isccp(1:ncol)                   = cospOUT%isccp_meanptop
       meantau_isccp(1:ncol)                    = cospOUT%isccp_meantaucld
       meancldalb_isccp(1:ncol)                 = cospOUT%isccp_meanalbedocld
       meantb_isccp(1:ncol)                     = cospOUT%isccp_meantb
       meantbclr_isccp(1:ncol)                  = cospOUT%isccp_meantbclr
    endif
    
    ! MISR
    if (lmisr_sim) then
       clMISR(1:ncol,1:ntau_cosp,1:nhtmisr_cosp) = cospOUT%misr_fq
    endif
    
    ! MODIS
    if (lmodis_sim) then
       cltmodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_Total_Mean
       clwmodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_Water_Mean
       climodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_Ice_Mean
       clhmodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_High_Mean
       clmmodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_Mid_Mean
       cllmodis(1:ncol)     = cospOUT%modis_Cloud_Fraction_Low_Mean
       tautmodis(1:ncol)    = cospOUT%modis_Optical_Thickness_Total_Mean
       tauwmodis(1:ncol)    = cospOUT%modis_Optical_Thickness_Water_Mean
       tauimodis(1:ncol)    = cospOUT%modis_Optical_Thickness_Ice_Mean
       tautlogmodis(1:ncol) = cospOUT%modis_Optical_Thickness_Total_LogMean
       tauwlogmodis(1:ncol) = cospOUT%modis_Optical_Thickness_Water_LogMean
       tauilogmodis(1:ncol) = cospOUT%modis_Optical_Thickness_Ice_LogMean
       reffclwmodis(1:ncol) = cospOUT%modis_Cloud_Particle_Size_Water_Mean
       reffclimodis(1:ncol) = cospOUT%modis_Cloud_Particle_Size_Ice_Mean
       pctmodis(1:ncol)     = cospOUT%modis_Cloud_Top_Pressure_Total_Mean
       lwpmodis(1:ncol)     = cospOUT%modis_Liquid_Water_Path_Mean
       iwpmodis(1:ncol)     = cospOUT%modis_Ice_Water_Path_Mean
       clmodis(1:ncol,1:ntau_cosp_modis,1:nprs_cosp)  = cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure 
       clrimodis(1:ncol,1:ntau_cosp_modis,1:numMODISReffIceBins) = cospOUT%modis_Optical_Thickness_vs_ReffICE
       clrlmodis(1:ncol,1:ntau_cosp_modis,1:numMODISReffLiqBins) = cospOUT%modis_Optical_Thickness_vs_ReffLIQ
    endif
    
    ! Use COSP output to populate CAM collapsed output variables
    do i=1,ncol
       if (lradar_sim) then
          do ih=1,nht_cosp
             do id=1,CLOUDSAT_DBZE_BINS
                ihd=(ih-1)*CLOUDSAT_DBZE_BINS+id                     
                cfad_dbze94_cs(i,ihd) = cfad_dbze94(i,id,ih)
             end do
          end do
          do ihml=1,nlay
             do isc=1,nscol_cosp
                ihsc=(ihml-1)*nscol_cosp+isc                 
                dbze_cs(i,ihsc) = dbze94(i,isc,ihml)
             end do
          end do
       endif
       
       if (llidar_sim) then
          do ih=1,nht_cosp
             do is=1,nsr_cosp
                ihs=(ih-1)*nsr_cosp+is                       
                cfad_sr532_cal(i,ihs) = cfad_lidarsr532(i,is,ih)
             end do
          end do
          do ihml=1,nlay
             do isc=1,nscol_cosp
                ihsc=(ihml-1)*nscol_cosp+isc                 
                atb532_cal(i,ihsc) = atb532(i,isc,ihml)
             end do
          end do
       endif
       
       if (lmisr_sim) then
          do ihm=1,nhtmisr_cosp
             do it=1,ntau_cosp
                ihmt=(ihm-1)*ntau_cosp+it                    
                cld_misr(i,ihmt) = clMISR(i,it,ihm) 
             end do
          end do
       endif
       
       if (lmodis_sim) then
          do ip=1,nprs_cosp
             do it=1,ntau_cosp_modis
                ipt=(ip-1)*ntau_cosp_modis+it
                clmodis_cam(i,ipt) = clmodis(i,it,ip)
             end do
          end do
          do ip=1,numMODISReffIceBins
             do it=1,ntau_cosp_modis
                ipt=(ip-1)*ntau_cosp_modis+it
                clrimodis_cam(i,ipt) = clrimodis(i,it,ip)
             end do
          end do
          do ip=1,numMODISReffLiqBins
             do it=1,ntau_cosp_modis
                ipt=(ip-1)*ntau_cosp_modis+it
                clrlmodis_cam(i,ipt) = clrlmodis(i,it,ip)
             end do
          end do
       endif
       
       ! Subcolums 
       do ihml=1,nlay
          do isc=1,nscol_cosp
             ihsc=(ihml-1)*nscol_cosp+isc                 
             scops_out(i,ihsc) = frac_out(i,isc,ihml)
          end do
       end do   
    end do
    call t_stopf("output_copying")

    ! ######################################################################################
    ! Clean up
    ! ######################################################################################
    call t_startf("destroy_cospIN")
    call destroy_cospIN(cospIN)
    call t_stopf("destroy_cospIN")
    call t_startf("destroy_cospstateIN")
    call destroy_cospstateIN(cospstateIN)
    call t_stopf("destroy_cospstateIN")
    call t_startf("destroy_cospOUT")
    call destroy_cosp_outputs(cospOUT) 
    call t_stopf("destroy_cospOUT")
    
    ! ######################################################################################
    ! OUTPUT
    ! ######################################################################################
    call t_startf("writing_output")
    ! ISCCP OUTPUTS
    if (lisccp_sim) then
       call outfld('FISCCP1_COSP',clisccp2,     pcols,lchnk)
       call outfld('CLDTOT_ISCCP',cldtot_isccp, pcols,lchnk)
       !! weight meancldalb_isccp by the cloud fraction
       !! where there is no isccp cloud fraction, set meancldalb_isccp = R_UNDEF
       !! weight meanptop_isccp  by the cloud fraction
       !! where there is no isccp cloud fraction, set meanptop_isccp = R_UNDEF
       !! weight meantau_isccp by the cloud fraction
       !! where there is no isccp cloud fraction, set meantau_isccp = R_UNDEF
       where (cldtot_isccp(:ncol) .eq. R_UNDEF)
          meancldalb_isccp(:ncol) = R_UNDEF
          meanptop_isccp(:ncol)   = R_UNDEF
          meantau_isccp(:ncol)    = R_UNDEF
       elsewhere
          meancldalb_isccp(:ncol) = meancldalb_isccp(:ncol)*cldtot_isccp(:ncol)
          meanptop_isccp(:ncol)   = meanptop_isccp(:ncol)*cldtot_isccp(:ncol)
          meantau_isccp(:ncol)    = meantau_isccp(:ncol)*cldtot_isccp(:ncol)
       end where
       call outfld('MEANCLDALB_ISCCP',meancldalb_isccp,pcols,lchnk)
       call outfld('MEANPTOP_ISCCP',  meanptop_isccp,  pcols,lchnk)
       call outfld('MEANTAU_ISCCP',   meantau_isccp,   pcols,lchnk)
       call outfld('MEANTB_ISCCP',    meantb_isccp,    pcols,lchnk)
       call outfld('MEANTBCLR_ISCCP', meantbclr_isccp, pcols,lchnk)
    end if
    
    ! CALIPSO SIMULATOR OUTPUTS
    if (llidar_sim) then
       call outfld('CLDLOW_CAL',    cldlow_cal,     pcols,lchnk)
       call outfld('CLDMED_CAL',    cldmed_cal,     pcols,lchnk)
       call outfld('CLDHGH_CAL',    cldhgh_cal,     pcols,lchnk)
       call outfld('CLDTOT_CAL',    cldtot_cal,     pcols,lchnk)
       call outfld('CLDTOT_CAL_ICE',cldtot_cal_ice, pcols,lchnk) !+1.4
       call outfld('CLDTOT_CAL_LIQ',cldtot_cal_liq, pcols,lchnk)
       call outfld('CLDTOT_CAL_UN', cldtot_cal_un,  pcols,lchnk)
       call outfld('CLDHGH_CAL_ICE',cldhgh_cal_ice, pcols,lchnk)
       call outfld('CLDHGH_CAL_LIQ',cldhgh_cal_liq, pcols,lchnk)
       call outfld('CLDHGH_CAL_UN', cldhgh_cal_un,  pcols,lchnk)
       call outfld('CLDMED_CAL_ICE',cldmed_cal_ice, pcols,lchnk)
       call outfld('CLDMED_CAL_LIQ',cldmed_cal_liq, pcols,lchnk)
       call outfld('CLDMED_CAL_UN', cldmed_cal_un,  pcols,lchnk)
       call outfld('CLDLOW_CAL_ICE',cldlow_cal_ice, pcols,lchnk)
       call outfld('CLDLOW_CAL_LIQ',cldlow_cal_liq, pcols,lchnk)
       call outfld('CLDLOW_CAL_UN', cldlow_cal_un,  pcols,lchnk) !+1.4
       where (cld_cal(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air).  
          !! I'm not sure why COSP produces a mix of R_UNDEF and realvalue in the nht_cosp dimension.
          cld_cal(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL',        cld_cal,       pcols,lchnk)  !! fails check_accum if 'A'
       call outfld('MOL532_CAL',     mol532_cal,    pcols,lchnk)
       
       where (cfad_sr532_cal(:ncol,:nht_cosp*nsr_cosp) .eq. R_UNDEF)
          !! fails check_accum if this is set... with ht_cosp set relative to sea level, mix of R_UNDEF and realvalue
          !!            cfad_sr532_cal(:ncol,:nht_cosp*nsr_cosp) = R_UNDEF
          cfad_sr532_cal(:ncol,:nht_cosp*nsr_cosp) = 0.0_r8
       end where
       call outfld('CFAD_SR532_CAL',cfad_sr532_cal    ,pcols,lchnk)
       
       where (refl_parasol(:ncol,:nsza_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air).  
          refl_parasol(:ncol,:nsza_cosp) = 0
       end where
       call outfld('RFL_PARASOL',refl_parasol   ,pcols,lchnk) !!
       
       where (cld_cal_liq(:ncol,:nht_cosp) .eq. R_UNDEF) !+cosp1.4
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_liq(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_LIQ',cld_cal_liq    ,pcols,lchnk)  !!
       
       where (cld_cal_ice(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_ice(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_ICE',cld_cal_ice    ,pcols,lchnk)  !!
       
       where (cld_cal_un(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_un(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_UN',cld_cal_un    ,pcols,lchnk)  !!
       
       where (cld_cal_tmp(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_tmp(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_TMP',cld_cal_tmp    ,pcols,lchnk)  !!
       
       where (cld_cal_tmpliq(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_tmpliq(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_TMPLIQ',cld_cal_tmpliq    ,pcols,lchnk)  !!
       
       where (cld_cal_tmpice(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_tmpice(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_TMPICE',cld_cal_tmpice    ,pcols,lchnk)  !!
       
       where (cld_cal_tmpun(:ncol,:nht_cosp) .eq. R_UNDEF)
          !! setting missing values to 0 (clear air), likely below sea level
          cld_cal_tmpun(:ncol,:nht_cosp) = 0.0_r8
       end where
       call outfld('CLD_CAL_TMPUN',cld_cal_tmpun    ,pcols,lchnk)  !!  !+cosp1.4 

    end if
    
    ! RADAR SIMULATOR OUTPUTS
    if (lradar_sim) then
       where (cfad_dbze94_cs(:ncol,:nht_cosp*CLOUDSAT_DBZE_BINS) .eq. R_UNDEF)
          !! fails check_accum if this is set... with ht_cosp set relative to sea level, mix of R_UNDEF and realvalue 
          !           cfad_dbze94_cs(:ncol,:nht_cosp*CLOUDSAT_DBZE_BINS) = R_UNDEF
          cfad_dbze94_cs(:ncol,:nht_cosp*CLOUDSAT_DBZE_BINS) = 0.0_r8
       end where
       call outfld('CFAD_DBZE94_CS',cfad_dbze94_cs,   pcols, lchnk)
       call outfld('CLDTOT_CALCS',  cldtot_calcs,     pcols, lchnk)
       call outfld('CLDTOT_CS',     cldtot_cs,        pcols, lchnk)
       call outfld('CLDTOT_CS2',    cldtot_cs2,       pcols, lchnk)
       call outfld('CLD_CAL_NOTCS', cld_cal_notcs,    pcols, lchnk)
       call outfld('CS_NOPRECIP',   ptcloudsatflag0,  pcols, lchnk)
       call outfld('CS_RAINPOSS',   ptcloudsatflag1,  pcols, lchnk)
       call outfld('CS_RAINPROB',   ptcloudsatflag2,  pcols, lchnk)
       call outfld('CS_RAINCERT',   ptcloudsatflag3,  pcols, lchnk)
       call outfld('CS_SNOWPOSS',   ptcloudsatflag4,  pcols, lchnk)
       call outfld('CS_SNOWCERT',   ptcloudsatflag5,  pcols, lchnk)
       call outfld('CS_MIXPOSS',    ptcloudsatflag6,  pcols, lchnk)
       call outfld('CS_MIXCERT',    ptcloudsatflag7,  pcols, lchnk)
       call outfld('CS_RAINHARD',   ptcloudsatflag8,  pcols, lchnk)
       call outfld('CS_UN',         ptcloudsatflag9,  pcols, lchnk)
       call outfld('CS_PIA',        cloudsatpia,      pcols, lchnk)
    end if
    
    ! MISR SIMULATOR OUTPUTS
    if (lmisr_sim) then
       call outfld('CLD_MISR',cld_misr    ,pcols,lchnk)
    end if
    
    ! MODIS SIMULATOR OUTPUTS
    if (lmodis_sim) then
       call outfld('CLTMODIS',cltmodis    ,pcols,lchnk)
       call outfld('CLWMODIS',clwmodis    ,pcols,lchnk)
       call outfld('CLIMODIS',climodis    ,pcols,lchnk)
       call outfld('CLHMODIS',clhmodis    ,pcols,lchnk)
       call outfld('CLMMODIS',clmmodis    ,pcols,lchnk)
       call outfld('CLLMODIS',cllmodis    ,pcols,lchnk)
       
       !! where there is no cloud fraction or no retrieval, set to R_UNDEF, 
       !! otherwise weight retrieval by cloud fraction
       where ((cltmodis(:ncol) .eq. R_UNDEF) .or. (tautmodis(:ncol) .eq. R_UNDEF))
          tautmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction cltmodis
          tautmodis(:ncol) = tautmodis(:ncol)*cltmodis(:ncol)
       end where
       call outfld('TAUTMODIS',tautmodis    ,pcols,lchnk)
       
       where ((tauwmodis(:ncol) .eq. R_UNDEF) .or. (clwmodis(:ncol) .eq. R_UNDEF))
          tauwmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction clwmodis
          tauwmodis(:ncol) = tauwmodis(:ncol)*clwmodis(:ncol)
       end where
       call outfld('TAUWMODIS',tauwmodis    ,pcols,lchnk)
       
       where ((tauimodis(:ncol) .eq. R_UNDEF) .or. (climodis(:ncol) .eq. R_UNDEF))
          tauimodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction climodis
          tauimodis(:ncol) = tauimodis(:ncol)*climodis(:ncol)
       end where
       call outfld('TAUIMODIS',tauimodis    ,pcols,lchnk)
       
       where ((tautlogmodis(:ncol)  .eq. R_UNDEF) .or. (cltmodis(:ncol) .eq. R_UNDEF))
          tautlogmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction cltmodis
          tautlogmodis(:ncol) = tautlogmodis(:ncol)*cltmodis(:ncol)
       end where
       call outfld('TAUTLOGMODIS',tautlogmodis    ,pcols,lchnk)
       
       where ((tauwlogmodis(:ncol)  .eq. R_UNDEF) .or. (clwmodis(:ncol) .eq. R_UNDEF))
          tauwlogmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction clwmodis
          tauwlogmodis(:ncol) = tauwlogmodis(:ncol)*clwmodis(:ncol)
       end where
       call outfld('TAUWLOGMODIS',tauwlogmodis    ,pcols,lchnk)
       
       where ((tauilogmodis(:ncol)  .eq. R_UNDEF) .or. (climodis(:ncol) .eq. R_UNDEF)) 
          tauilogmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction climodis
          tauilogmodis(:ncol) = tauilogmodis(:ncol)*climodis(:ncol)
       end where
       call outfld('TAUILOGMODIS',tauilogmodis    ,pcols,lchnk)
       
       where ((reffclwmodis(:ncol)  .eq. R_UNDEF) .or. (clwmodis(:ncol) .eq. R_UNDEF)) 
          reffclwmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction clwmodis
          reffclwmodis(:ncol) = reffclwmodis(:ncol)*clwmodis(:ncol)
       end where
       call outfld('REFFCLWMODIS',reffclwmodis    ,pcols,lchnk)
       
       where ((reffclimodis(:ncol)  .eq. R_UNDEF) .or. (climodis(:ncol) .eq. R_UNDEF))
          reffclimodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction climodis
          reffclimodis(:ncol) = reffclimodis(:ncol)*climodis(:ncol)
       end where
       call outfld('REFFCLIMODIS',reffclimodis    ,pcols,lchnk)
       
       where ((pctmodis(:ncol)  .eq. R_UNDEF) .or. ( cltmodis(:ncol) .eq. R_UNDEF))
          pctmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction cltmodis
          pctmodis(:ncol) = pctmodis(:ncol)*cltmodis(:ncol)
       end where
       call outfld('PCTMODIS',pctmodis    ,pcols,lchnk)
       
       where ((lwpmodis(:ncol)  .eq. R_UNDEF) .or. (clwmodis(:ncol) .eq. R_UNDEF))
          lwpmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction clwmodis
          lwpmodis(:ncol) = lwpmodis(:ncol)*clwmodis(:ncol)
       end where
       call outfld('LWPMODIS',lwpmodis    ,pcols,lchnk)
       
       where ((iwpmodis(:ncol)  .eq. R_UNDEF) .or. (climodis(:ncol) .eq. R_UNDEF))
          iwpmodis(:ncol) = R_UNDEF
       elsewhere
          !! weight by the cloud fraction climodis
          iwpmodis(:ncol) = iwpmodis(:ncol)*climodis(:ncol)
       end where
       call outfld('IWPMODIS',iwpmodis    ,pcols,lchnk)
       
       call outfld('CLMODIS',clmodis_cam  ,pcols,lchnk) 
       call outfld('CLRIMODIS',clrimodis_cam  ,pcols,lchnk) 
       call outfld('CLRLMODIS',clrlmodis_cam  ,pcols,lchnk) 
    end if
    
    ! SUB-COLUMN OUTPUT
    if (lfrac_out) then
       call outfld('SCOPS_OUT',scops_out   ,pcols,lchnk)!!!-1.00000E+30 !! fails check_accum if 'A'
       if (lisccp_sim) then
          call outfld('TAU_ISCCP',    tau_isccp,    pcols,lchnk) !! fails check_accum if 'A'
          call outfld('CLDPTOP_ISCCP',cldptop_isccp,pcols,lchnk) !! fails check_accum if 'A'
       end if
       if (llidar_sim) then
          call outfld('ATB532_CAL',atb532_cal,pcols,lchnk) !! fails check_accum if 'A'
       end if
       if (lradar_sim) then
          call outfld('DBZE_CS',dbze_cs,pcols,lchnk) !! fails check_accum if 'A'
       end if
    end if
    call t_stopf("writing_output")
#endif
  end subroutine cospsimulator_intr_run

#ifdef USE_COSP
  ! ######################################################################################
  ! SUBROUTINE subsample_and_optics
  ! ######################################################################################
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro,overlap,            &
                                  lidar_ice_type, sd, tca, cca,                          &
                                  fl_lsrainIN, fl_lssnowIN, fl_lsgrplIN, fl_ccrainIN,    &
                                  fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,   &
                                  reffIN, dtau_c, dtau_s, dem_c, dem_s, dtau_s_snow,     &
                                  dem_s_snow, sfcP, cospstateIN, cospIN)
    ! Dependencies
    use cosp_kinds,           only: wp
    use mod_rng,              only: rng_state, init_rng
    use mod_cosp_config,      only: R_UNDEF
    use mod_scops,            only: scops
    use mod_prec_scops,       only: prec_scops
    use mod_cosp_utils,       only: cosp_precip_mxratio
    use mod_quickbeam_optics, only: quickbeam_optics, gases
    use cosp_optics,          only: cosp_simulator_optics,lidar_optics,modis_optics,    &
                                    modis_optics_partition
    use mod_cosp_config,      only: Nlvgrid, vgrid_zl, vgrid_zu
    use mod_cosp_stats,       only: cosp_change_vertical_grid
    ! Inputs
    integer,intent(in) :: &
         nPoints,      & ! Number of gridpoints
         nLevels,      & ! Number of vertical levels
         nColumns,     & ! Number of subcolumns
         nHydro,       & ! Number pf hydrometeor types
         overlap,      & ! Overlap assumption (1/2/3)
         lidar_ice_type  ! Ice type assumption used by lidar optics
    real(wp),intent(in),dimension(nPoints,nLevels) :: &
         tca,          & ! Total cloud amount (0-1)
         cca,          & ! Convective cloud amount (0-1)
         mr_lsliq,     & ! Mixing ratio (kg/kg)
         mr_lsice,     & ! Mixing ratio (kg/kg)
         mr_ccliq,     & ! Mixing ratio (kg/kg)
         mr_ccice,     & ! Mixing ratio (kg/kg)
         dtau_c,       & ! 0.67-micron optical depth (convective)
         dtau_s,       & ! 0.67-micron optical depth (stratiform)
         dem_c,        & ! 11-micron emissivity (convective)
         dem_s,        & ! 11-micron emissivity (stratiform)
         fl_lsrainIN,  & ! Precipitation flux
         fl_lssnowIN,  & ! Precipitation flux
         fl_lsgrplIN,  & ! Precipitation flux
         fl_ccrainIN,  & ! Precipitation flux
         fl_ccsnowIN     ! Precipitation flux
    real(wp),intent(inout),dimension(nPoints,nLevels) :: &    
         dtau_s_snow,  & ! 0.67-micron optical depth (snow)
         dem_s_snow      ! 11-micron emissivity (snow)
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: &
         reffIN          !
    real(wp),intent(in),dimension(nPoints) :: &
         sfcP            ! Surface pressure 
    type(size_distribution),intent(inout) :: &
         sd
    
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN
    
    ! Local variables
    integer :: i, j, k, istat
    real(wp),dimension(nPoints,nLevels)      :: column_frac_out,column_prec_out,         &
                                                fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain, &
                                                fl_ccsnow
    real(wp),dimension(nPoints,nLevels,nHydro) :: ReffTemp                                                
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable         :: seed
    real(wp),dimension(:,:),allocatable      :: ls_p_rate,cv_p_rate,frac_ls,frac_cv,     &
                                                prec_ls,prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable  :: frac_prec,&
                                                 MODIS_cloudWater,MODIS_cloudIce,        &
                                                 MODIS_watersize,MODIS_iceSize,          &
                                                 MODIS_snowSize,MODIS_cloudSnow,         &
                                                 MODIS_opticalThicknessLiq,              &
                                                 MODIS_opticalThicknessSnow,             &
                                                 MODIS_opticalThicknessIce,              &
                                                 fracPrecipIce, fracPrecipIce_statGrid
    real(wp),dimension(:,:,:,:),allocatable   :: mr_hydro,Reff,Np

    character(len=*), parameter :: sub = 'subsample_and_optics'
    !--------------------------------------------------------------------------------------
             
    call t_startf("scops")
    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumn generation
       allocate(rngs(nPoints), seed(nPoints), stat=istat)
       call handle_allocate_error(istat, sub, 'rngs, seed')
       seed = int(sfcP)
       if (Npoints .gt. 1) seed=(sfcP-int(sfcP))*1000000 
       call init_rng(rngs, seed)
   
       ! Call scops
       call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
       deallocate(seed,rngs)
       
       ! Sum up precipitation rates.
       allocate(ls_p_rate(nPoints,nLevels), cv_p_rate(nPoints,Nlevels), stat=istat)
       call handle_allocate_error(istat, sub, 'ls_p_rate, cv_p_rate')
       ls_p_rate(:,1:nLevels) = fl_lsrainIN + fl_lssnowIN + fl_lsgrplIN
       cv_p_rate(:,1:nLevels) = fl_ccrainIN + fl_ccsnowIN
       
       ! Call PREC_SCOPS
       allocate(frac_prec(nPoints,nColumns,nLevels), stat=istat)
       call handle_allocate_error(istat, sub, 'frac_prec')
       call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
       deallocate(ls_p_rate,cv_p_rate)
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute precipitation fraction in each gridbox
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels), &
                frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels), stat=istat)
       call handle_allocate_error(istat, sub, 'frac_ls,..,prec_cv')

       ! Initialize
       frac_ls(1:nPoints,1:nLevels) = 0._wp
       prec_ls(1:nPoints,1:nLevels) = 0._wp
       frac_cv(1:nPoints,1:nLevels) = 0._wp
       prec_cv(1:nPoints,1:nLevels) = 0._wp
       do j=1,nPoints
          do k=1,nLevels
             do i=1,nColumns
                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 1)         prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 2)         prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)         prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)         prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns

             ! Adjust grid-box mean snow properties to local properties
             ! Convert longwave optical depth to longwave emissivity
             if (prec_ls(j,k) .ne. 0._r8 .and. dtau_s_snow(j,k) .gt. 0._r8) then
                dtau_s_snow(j,k) = dtau_s_snow(j,k)/prec_ls(j,k) 
             end if
             if (prec_ls(j,k) .ne. 0._r8 .and. dem_s_snow(j,k) .gt. 0._r8) then
                dem_s_snow(j,k) = dem_s_snow(j,k)/prec_ls(j,k)
                dem_s_snow(j,k) = 1._r8 - exp ( -1._r8*dem_s_snow(j,k))
             end if !!+JEK
          enddo
       enddo
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
       ! and precipitation
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro), &
                Reff(nPoints,nColumns,nLevels,nHydro),     &
                Np(nPoints,nColumns,nLevels,nHydro), stat=istat)
       call handle_allocate_error(istat, sub, 'mr_hydro,Reff,Np')

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       
       do k=1,nColumns
          ! Subcolumn clouds
          column_frac_out = cospIN%frac_out(:,k,:)
               
          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq
             mr_hydro(:,k,:,I_LSCICE) = mr_lsice
             Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,:,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = ReffIN(:,:,I_LSCICE)
          ! CONV clouds   
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq
             mr_hydro(:,k,:,I_CVCICE) = mr_ccice
             Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,:,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = ReffIN(:,:,I_CVCICE)
          end where
          
          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,k,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,k,:,I_LSRAIN) = ReffIN(:,:,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = ReffIN(:,:,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = ReffIN(:,:,I_LSGRPL)
          ! CONV precipitation   
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = ReffIN(:,:,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = ReffIN(:,:,I_CVSNOW)
          end where
       enddo

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
       ! the fraction-based values
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       do k=1,nLevels
          do j=1,nPoints
             ! Clouds
             if (frac_ls(j,k) .ne. 0._r8) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             endif
             if (frac_cv(j,k) .ne. 0._r8) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             endif
             
             ! Precipitation
             if (prec_ls(j,k) .ne. 0._r8) then
                fl_lsrain(j,k) = fl_lsrainIN(j,k)/prec_ls(j,k)
                fl_lssnow(j,k) = fl_lssnowIN(j,k)/prec_ls(j,k)
                fl_lsgrpl(j,k) = fl_lsgrplIN(j,k)/prec_ls(j,k)
             endif
             if (prec_cv(j,k) .ne. 0._r8) then
                fl_ccrain(j,k) = fl_ccrainIN(j,k)/prec_cv(j,k)
                fl_ccsnow(j,k) = fl_ccsnowIN(j,k)/prec_cv(j,k)
             endif
          enddo
       enddo
             
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! LS rain
       call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
            cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
            alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
            a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
            gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
            mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
       ! LS snow
       call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
            cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
            alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
            a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
            gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
            mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
       ! CV rain
       call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
            cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
            alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
            a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
            gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
            mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
       ! CV snow
       call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
            cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
            alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
            a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
            gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
            mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
       ! LS groupel.
       call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
            cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
            alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
            a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
            gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
            mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))

    else
       cospIN%frac_out(:,:,:) = 1  
       allocate(mr_hydro(nPoints, 1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),      &
                Np(nPoints,1,nLevels,nHydro), stat=istat)
       call handle_allocate_error(istat, sub, 'mr_hydro,Reff,Np')
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice
       Reff(:,1,:,:)            = ReffIN
    endif
    call t_stopf("scops")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("cloudsat_optics")
    if (lradar_sim) then
       ! Compute gaseous absorption (assume identical for each subcolun)
       allocate(g_vol(nPoints,nLevels), stat=istat)
       call handle_allocate_error(istat, sub, 'g_vol')
       g_vol(:,:)=0._wp
       do i = 1, nPoints
          do j = 1, nLevels
             if (cospIN%rcfg_cloudsat%use_gas_abs == 1 .or. &
                (cospIN%rcfg_cloudsat%use_gas_abs == 2 .and. j == 1)) then
                g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),    &
                                   cospstateIN%qv(i,j), cospIN%rcfg_cloudsat%freq)
             endif
             cospIN%g_vol_cloudsat(i,:,j) = g_vol(i,j)
          end do
       end do

       ! Loop over all subcolumns
       allocate(fracPrecipIce(nPoints,nColumns,nLevels), stat=istat)
       call handle_allocate_error(istat, sub, 'fracPrecipIce')
       fracPrecipIce(:,:,:) = 0._wp
       do k=1,nColumns
          call quickbeam_optics(sd, cospIN%rcfg_cloudsat, nPoints, nLevels, R_UNDEF, &
               mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,      &
               Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,                &
               cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),                 &
               cospIN%kr_vol_cloudsat(1:nPoints,k,:))

        ! At each model level, what fraction of the precipitation is frozen?
          where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
             fracPrecipIce(:,k,:) = (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + &
                  mr_hydro(:,k,:,I_LSGRPL)) / &
                  (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                  mr_hydro(:,k,:,I_LSRAIN)  + mr_hydro(:,k,:,I_CVRAIN))
          elsewhere
             fracPrecipIce(:,k,:) = 0._wp
          endwhere
       enddo

       ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
       allocate(fracPrecipIce_statGrid(nPoints,nColumns,Nlvgrid), stat=istat)
       call handle_allocate_error(istat, sub, 'fracPrecipIce_statGrid')
       fracPrecipIce_statGrid(:,:,:) = 0._wp
       call cosp_change_vertical_grid(Npoints, Ncolumns, Nlevels, cospstateIN%hgt_matrix(:,Nlevels:1:-1), &
            cospstateIN%hgt_matrix_half(:,Nlevels:1:-1), fracPrecipIce(:,:,Nlevels:1:-1), Nlvgrid,  &
            vgrid_zl(Nlvgrid:1:-1),  vgrid_zu(Nlvgrid:1:-1), fracPrecipIce_statGrid(:,:,Nlvgrid:1:-1))

       ! For near-surface diagnostics, we only need the frozen fraction at one layer.
       cospIN%fracPrecipIce(:,:) = fracPrecipIce_statGrid(:,:,cloudsat_preclvl)
       
    endif
    call t_stopf("cloudsat_optics")
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CALIPSO Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("calipso_optics")
    if (Llidar_sim) then
       ReffTemp = ReffIN
       call lidar_optics(nPoints,nColumns,nLevels,5,lidar_ice_type,                      &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSCLIQ),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSCICE),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_CVCLIQ),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_CVCICE),              &
                         mr_hydro(1:nPoints,1:nColumns,1:nLevels,I_LSSNOW),              &
                         ReffTemp(1:nPoints,1:nLevels,I_LSCLIQ),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_LSCICE),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_CVCLIQ),                         &
                         ReffTemp(1:nPoints,1:nLevels,I_CVCICE),                         & 
                         ReffTemp(1:nPoints,1:nLevels,I_LSSNOW),                         &
                         cospstateIN%pfull(1:nPoints,1:nLevels),                         &
                         cospstateIN%phalf(1:nPoints,1:nLevels+1),                       &
                         cospstateIN%at(1:nPoints,1:nLevels),                            &
                         cospIN%beta_mol_calipso(1:nPoints,1:nLevels),                   &
                         cospIN%betatot_calipso(1:nPoints,1:nColumns,1:nLevels),         &
                         cospIN%tau_mol_calipso(1:nPoints,1:nLevels),                    &
                         cospIN%tautot_calipso(1:nPoints,1:nColumns,1:nLevels),          &
                         cospIN%tautot_S_liq(1:nPoints,1:nColumns),                      &
                         cospIN%tautot_S_ice(1:nPoints,1:nColumns),                      &
                         cospIN%betatot_ice_calipso(1:nPoints,1:nColumns,1:nLevels),     &
                         cospIN%betatot_liq_calipso(1:nPoints,1:nColumns,1:nLevels),     &
                         cospIN%tautot_ice_calipso(1:nPoints,1:nColumns,1:nLevels),      &
                         cospIN%tautot_liq_calipso(1:nPoints,1:nColumns,1:nLevels)) 
    endif
    call t_stopf("calipso_optics")

    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Compute optical fields for passive simulators (i.e. only sunlit points)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity (needed by the ISCCP simulator)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("11micron_emissivity")
    if (Lisccp_sim) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,  &
            cospIN%emiss_11)
       ! Add in contributions from radiative snow 
       do j=1,nColumns
          where(frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3)
             cospIN%emiss_11(:,j,:) = 1._wp - (1- cospIN%emiss_11(:,j,:))*(1-dem_s_snow)
          endwhere
       enddo
    endif
    call t_stopf("11micron_emissivity")
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth (needed by ISCCP, MISR and MODIS simulators)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("067tau")
    if (Lisccp_sim .or. Lmisr_sim .or. Lmodis_sim) then
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,&
            cospIN%tau_067)
       
       ! Add in contributions from snow 
       do j=1,nColumns
          where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
               Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_s_snow .gt. 0._r8)
             cospIN%tau_067(:,j,:)  = cospIN%tau_067(:,j,:)+dtau_s_snow
          endwhere
       enddo
    endif
    call t_stopf("067tau")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call t_startf("modis_optics")
    if (lmodis_sim) then
       allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                              &
                MODIS_cloudIce(nPoints,nColumns,nLevels),                                &
                MODIS_cloudSnow(nPoints,nColumns,nLevels),                               &
                MODIS_waterSize(nPoints,nColumns,nLevels),                               &
                MODIS_iceSize(nPoints,nColumns,nLevels),                                 &
                MODIS_snowSize(nPoints,nColumns,nLevels),                                &
                MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                     &
                MODIS_opticalThicknessIce(nPoints,nColumns,nLevels),                     &
                MODIS_opticalThicknessSnow(nPoints,nColumns,nLevels), stat=istat)
       call handle_allocate_error(istat, sub, 'MODIS_*')

       ! Cloud water
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
       ! Cloud ice
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)  
       ! Cloud water droplet size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
       ! Cloud ice crystal size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,              &
            Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)
       
       ! Cloud snow and size	
       MODIS_snowSize(:,:,:)  = Reff(:,:,:,I_LSSNOW)
       do j=1,nColumns
          where((frac_prec(:,j,:) .eq. 1 .or. frac_prec(:,j,:) .eq. 3) .and. &
               Reff(:,j,:,I_LSSNOW) .gt. 0._r8 .and. dtau_s_snow .gt. 0._r8)
             MODIS_cloudSnow(:,j,:) = mr_hydro(:,j,:,I_LSSNOW)
             MODIS_snowSize(:,j,:)  = Reff(:,j,:,I_LSSNOW)
          elsewhere
             MODIS_snowSize(:,j,:)  = 0._wp
             MODIS_cloudSnow(:,j,:) = 0._wp
          endwhere
       enddo
       
       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,     &
            MODIS_cloudIce, MODIS_cloudSnow, MODIS_waterSize, MODIS_iceSize,         &
            MODIS_snowSize, cospIN%tau_067, MODIS_opticalThicknessLiq,               &
            MODIS_opticalThicknessIce, MODIS_opticalThicknessSnow)                            
       
       ! Compute asymmetry parameter and single scattering albedo 
       call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,      &
            MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                     &
            MODIS_iceSize*1.0e6_wp, MODIS_opticalThicknessSnow,                      &
            MODIS_snowSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)

    endif ! MODIS simulator optics
    call t_stopf("modis_optics")

  end subroutine subsample_and_optics
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(npoints,ncolumns,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y

    ! local
    integer :: istat
    character(len=*), parameter :: sub = 'construct_cospIN'
    !--------------------------------------------------------------------------------------
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    
    allocate(y%tau_067(            npoints, ncolumns, nlevels),&
             y%emiss_11(           npoints, ncolumns, nlevels),&
             y%frac_out(           npoints, ncolumns, nlevels),&
             y%betatot_calipso(    npoints, ncolumns, nlevels),&
             y%betatot_ice_calipso(npoints, ncolumns, nlevels),&
             y%fracLiq(            npoints, ncolumns, nlevels),&
             y%betatot_liq_calipso(npoints, ncolumns, nlevels),&
             y%tautot_calipso(     npoints, ncolumns, nlevels),&
             y%tautot_ice_calipso( npoints, ncolumns, nlevels),&
             y%tautot_liq_calipso( npoints, ncolumns, nlevels),&
             y%z_vol_cloudsat(     npoints, ncolumns, nlevels),&
             y%kr_vol_cloudsat(    npoints, ncolumns, nlevels),&
             y%g_vol_cloudsat(     npoints, ncolumns, nlevels),&
             y%asym(               npoints, ncolumns, nlevels),&
             y%ss_alb(             npoints, ncolumns, nlevels),&
             y%beta_mol_calipso(   npoints,           nlevels),&
             y%tau_mol_calipso(    npoints,           nlevels),&
             y%tautot_S_ice(       npoints, ncolumns         ),&
             y%tautot_S_liq(       npoints, ncolumns)         ,&
             y%fracPrecipIce(npoints,   ncolumns), stat=istat)
    call handle_allocate_error(istat, sub, 'tau_067,..,fracPrecipIce')

  end subroutine construct_cospIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y
    
    ! local
    integer :: istat
    character(len=*), parameter :: sub = 'construct_cospstateIN'
    !--------------------------------------------------------------------------------------

    allocate( &
       y%sunlit(npoints), &
       y%at(npoints,nlevels), &
       y%pfull(npoints,nlevels), &
       y%phalf(npoints,nlevels+1), &
       y%qv(npoints,nlevels), &
       y%hgt_matrix(npoints,nlevels), &
       y%hgt_matrix_half(npoints,nlevels+1), &
       y%land(npoints), &
       y%skt(npoints), &
       y%surfelev(nPoints), &
       y%emis_sfc(nchan), &
       y%u_sfc(npoints), &
       y%v_sfc(npoints), &
       y%seaice(npoints), &
       y%lat(npoints), &
       y%lon(nPoints), &
       y%o3(npoints,nlevels), &
       y%tca(nPoints,nLevels), &
       y%cloudIce(nPoints,nLevels), &
       y%cloudLiq(nPoints,nLevels), &
       y%fl_rain(nPoints,nLevels), &
       y%fl_snow(nPoints,nLevels), stat=istat)
    call handle_allocate_error(istat, sub, 'sunlit,..,fl_snow')

  end subroutine construct_cospstateIN
  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################  
  subroutine construct_cosp_outputs(Npoints,Ncolumns,Nlevels,Nlvgrid,x)
    ! Inputs
    integer,intent(in) :: &
         Npoints,         & ! Number of sampled points
         Ncolumns,        & ! Number of subgrid columns
         Nlevels,         & ! Number of model levels
         Nlvgrid            ! Number of levels in L3 stats computation
    
    ! Outputs
    type(cosp_outputs),intent(out) :: &
         x           ! COSP output structure  

    ! local
    integer :: istat
    character(len=*), parameter :: sub = 'construct_cosp_outputs'
    !--------------------------------------------------------------------------------------
  
     ! ISCCP simulator outputs
    if (lisccp_sim) then
       allocate( &
          x%isccp_boxtau(Npoints,Ncolumns), &
          x%isccp_boxptop(Npoints,Ncolumns), &
          x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins), &
          x%isccp_totalcldarea(Npoints), &
          x%isccp_meanptop(Npoints), &
          x%isccp_meantaucld(Npoints), &
          x%isccp_meantb(Npoints), &
          x%isccp_meantbclr(Npoints), &
          x%isccp_meanalbedocld(Npoints), stat=istat)
       call handle_allocate_error(istat, sub, 'isccp_*')
    endif

    ! MISR simulator
    if (lmisr_sim) then 
       allocate( &
          x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins), &
          ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
          !        they are still computed. Should probably have a logical to control these
          !        outputs.
          x%misr_dist_model_layertops(Npoints,numMISRHgtBins), &
          x%misr_meanztop(Npoints), &
          x%misr_cldarea(Npoints), stat=istat)
       call handle_allocate_error(istat, sub, 'misr_*')
    endif
    
    ! MODIS simulator
    if (lmodis_sim) then
       allocate( &
          x%modis_Cloud_Fraction_Total_Mean(Npoints), &
          x%modis_Cloud_Fraction_Water_Mean(Npoints), &
          x%modis_Cloud_Fraction_Ice_Mean(Npoints), &
          x%modis_Cloud_Fraction_High_Mean(Npoints), &
          x%modis_Cloud_Fraction_Mid_Mean(Npoints), &
          x%modis_Cloud_Fraction_Low_Mean(Npoints), &
          x%modis_Optical_Thickness_Total_Mean(Npoints), &
          x%modis_Optical_Thickness_Water_Mean(Npoints), &
          x%modis_Optical_Thickness_Ice_Mean(Npoints), &
          x%modis_Optical_Thickness_Total_LogMean(Npoints), &
          x%modis_Optical_Thickness_Water_LogMean(Npoints), &
          x%modis_Optical_Thickness_Ice_LogMean(Npoints), &
          x%modis_Cloud_Particle_Size_Water_Mean(Npoints), &
          x%modis_Cloud_Particle_Size_Ice_Mean(Npoints), &
          x%modis_Cloud_Top_Pressure_Total_Mean(Npoints), &
          x%modis_Liquid_Water_Path_Mean(Npoints), &
          x%modis_Ice_Water_Path_Mean(Npoints), &
          x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins), &
          x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins), &   
          x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins), &
          stat=istat)
       call handle_allocate_error(istat, sub, 'modis_*')
    endif
    
    ! CALIPSO simulator
    if (llidar_sim) then
       allocate( &
          x%calipso_beta_mol(Npoints,Nlevels), &
          x%calipso_beta_tot(Npoints,Ncolumns,Nlevels), &
          x%calipso_srbval(SR_BINS+1), &
          x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid), &
          x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels), &  
          x%calipso_lidarcld(Npoints,Nlvgrid), &
          x%calipso_cldlayer(Npoints,LIDAR_NCAT), &        
          x%calipso_lidarcldphase(Npoints,Nlvgrid,6), &
          x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5), &
          x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6), &     
          x%calipso_tau_tot(Npoints,Ncolumns,Nlevels), &       
          x%calipso_temp_tot(Npoints,Nlevels), stat=istat)               
       call handle_allocate_error(istat, sub, 'calipso_*')
    endif 
      
    ! PARASOL
    if (lparasol_sim) then
       allocate( &
          x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL), &
          x%parasolGrid_refl(Npoints,PARASOL_NREFL), stat=istat)
       call handle_allocate_error(istat, sub, 'parasol*')
    endif

    ! Cloudsat simulator
    if (lradar_sim) then
       allocate( &
          x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels), &
          x%cloudsat_cfad_ze(Npoints,CLOUDSAT_DBZE_BINS,Nlvgrid), &
          x%lidar_only_freq_cloud(Npoints,Nlvgrid), &
          x%radar_lidar_tcc(Npoints), &
          x%cloudsat_precip_cover(Npoints,nCloudsatPrecipClass), &
          x%cloudsat_pia(Npoints), stat=istat)
       call handle_allocate_error(istat, sub, 'cloudsat*')
    endif

  end subroutine construct_cosp_outputs

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y

    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%fracPrecipIce))       deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%surfelev))        deallocate(y%surfelev)
    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)    
    
  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol))          then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)     
     endif
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)     
     endif
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)    
        nullify(y%calipso_beta_tot)     
     endif
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot) 
        nullify(y%calipso_tau_tot)     
     endif
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)     
     endif
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)     
     endif
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)     
     endif
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)     
     endif
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     endif
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     endif
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)     
     endif
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl) 
        nullify(y%parasolGrid_refl)     
     endif
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot) 
        nullify(y%cloudsat_Ze_tot)  
     endif
     if (associated(y%cloudsat_precip_cover)) then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia)) then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)     
     endif
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc) 
        nullify(y%radar_lidar_tcc)  
     endif
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)     
     endif
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea) 
        nullify(y%isccp_totalcldarea)  
     endif
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb) 
        nullify(y%isccp_meantb)     
     endif
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)  
     endif
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)     
     endif
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld) 
        nullify(y%isccp_meantaucld)       
     endif
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)     
     endif
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)       
     endif
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)     
     endif
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)       
     endif
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq) 
        nullify(y%misr_fq)     
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)       
     endif
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)     
     endif
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)      
     endif
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)     
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
        nullify(y%modis_Cloud_Fraction_Total_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
        nullify(y%modis_Cloud_Fraction_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
        nullify(y%modis_Cloud_Fraction_Water_Mean)           
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
        nullify(y%modis_Cloud_Fraction_High_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
        nullify(y%modis_Cloud_Fraction_Mid_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
        nullify(y%modis_Cloud_Fraction_Low_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
        nullify(y%modis_Optical_Thickness_Total_Mean)  
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
        nullify(y%modis_Optical_Thickness_Water_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
        nullify(y%modis_Optical_Thickness_Ice_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
        nullify(y%modis_Optical_Thickness_Total_LogMean)    
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
        nullify(y%modis_Optical_Thickness_Water_LogMean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
        nullify(y%modis_Optical_Thickness_Ice_LogMean)     
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)     
        nullify(y%modis_Liquid_Water_Path_Mean)     
     endif
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)       
        nullify(y%modis_Ice_Water_Path_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
     if (associated(y%calipso_cldtype)) then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)
     endif
     if (associated(y%calipso_cldtypetemp)) then
        deallocate(y%calipso_cldtypetemp) 
        nullify(y%calipso_cldtypetemp) 
     endif
     if (associated(y%calipso_cldtypemeanz)) then
        deallocate(y%calipso_cldtypemeanz) 
        nullify(y%calipso_cldtypemeanz) 
     endif
     if (associated(y%calipso_cldtypemeanzse)) then
        deallocate(y%calipso_cldtypemeanzse) 
        nullify(y%calipso_cldtypemeanzse) 
     endif
     if (associated(y%calipso_cldthinemis)) then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif
     if (associated(y%calipso_lidarcldtype)) then
        deallocate(y%calipso_lidarcldtype)
        nullify(y%calipso_lidarcldtype)
     endif
        
   end subroutine destroy_cosp_outputs
#endif

!#######################################################################
end module cospsimulator_intr
