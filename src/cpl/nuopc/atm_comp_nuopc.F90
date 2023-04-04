module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CAM
  !----------------------------------------------------------------------------

   use ESMF                , only : operator(<=), operator(>), operator(==), operator(+)
   use ESMF                , only : ESMF_MethodRemove
   use ESMF                , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_State, ESMF_StateGet
   use ESMF                , only : ESMF_Grid, ESMF_GridCreateNoPeriDimUfrm, ESMF_Field, ESMF_FieldGet
   use ESMF                , only : ESMF_DistGrid, ESMF_DistGridCreate
   use ESMF                , only : ESMF_Mesh, ESMF_MeshCreate, ESMF_MeshGet, ESMF_FILEFORMAT_ESMFMESH
   use ESMF                , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockGetNextTime, ESMF_ClockAdvance
   use ESMF                , only : ESMF_Time, ESMF_TimeGet
   use ESMF                , only : ESMF_Alarm, ESMF_ClockGetAlarm, ESMF_AlarmRingerOff, ESMF_AlarmIsRinging
   use ESMF                , only : ESMF_ClockGetAlarmList, ESMF_ALARMLIST_ALL, ESMF_AlarmSet
   use ESMF                , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
   use ESMF                , only : ESMF_CalKind_Flag, ESMF_MAXSTR, ESMF_KIND_I8
   use ESMF                , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
   use ESMF                , only : ESMF_GridCompSetEntryPoint
   use ESMF                , only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
   use ESMF                , only : ESMF_LOGMSG_INFO, ESMF_LOGERR_PASSTHRU
   use ESMF                , only : ESMF_LogWrite, ESMF_LogSetError, ESMF_LogFoundError
   use ESMF                , only : ESMF_SUCCESS, ESMF_METHOD_INITIALIZE, ESMF_FAILURE, ESMF_RC_NOT_VALID
   use ESMF                , only : ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
   use NUOPC               , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
   use NUOPC               , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
   use NUOPC               , only : NUOPC_CompAttributeGet, NUOPC_Advertise
   use NUOPC               , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
   use NUOPC_Model         , only : model_routine_SS           => SetServices
   use NUOPC_Model         , only : SetVM
   use NUOPC_Model         , only : model_label_Advance        => label_Advance
   use NUOPC_Model         , only : model_label_DataInitialize => label_DataInitialize
   use NUOPC_Model         , only : model_label_SetRunClock    => label_SetRunClock
   use NUOPC_Model         , only : model_label_Finalize       => label_Finalize
   use NUOPC_Model         , only : NUOPC_ModelGet
   use shr_kind_mod        , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
   use shr_sys_mod         , only : shr_sys_abort
   use shr_log_mod         , only : shr_log_getlogunit, shr_log_setlogunit
   use shr_cal_mod         , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
   use shr_const_mod       , only : shr_const_pi
   use shr_orb_mod         , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
   use cam_instance        , only : cam_instance_init, inst_suffix, inst_index
   use cam_comp            , only : cam_init, cam_run1, cam_run2, cam_run3, cam_run4, cam_final
   use camsrfexch          , only : cam_out_t, cam_in_t
   use radiation           , only : nextsw_cday
   use cam_logfile         , only : iulog
   use spmd_utils          , only : spmdinit, masterproc, iam, mpicom
   use time_manager        , only : get_curr_calday, advance_timestep, get_curr_date, get_nstep, get_step_size
   use atm_import_export   , only : read_surface_fields_namelists, advertise_fields, realize_fields
   use atm_import_export   , only : import_fields, export_fields
   use nuopc_shr_methods   , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
   use nuopc_shr_methods   , only : set_component_logging, get_component_instance, log_clock_advance
   use perf_mod            , only : t_startf, t_stopf
   use ppgrid              , only : pcols, begchunk, endchunk
   use dyn_grid            , only : get_horiz_grid_dim_d
   use phys_grid           , only : get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
   use phys_grid           , only : ngcols=>num_global_phys_cols
   use cam_control_mod     , only : cam_ctrl_set_orbit
   use cam_pio_utils       , only : cam_pio_createfile, cam_pio_openfile, cam_pio_closefile, pio_subsystem
   use cam_initfiles       , only : cam_initfiles_get_caseid, cam_initfiles_get_restdir
   use cam_history_support , only : fillvalue
   use filenames           , only : interpret_filename_spec
   use pio                 , only : file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim, PIO_MAX_NAME
   use pio                 , only : pio_closefile, pio_put_att, pio_enddef, pio_nowrite
   use pio                 , only : pio_inq_dimid, pio_inq_varid, pio_inquire_dimension, pio_def_var
   use pio                 , only : pio_initdecomp, pio_freedecomp
   use pio                 , only : pio_read_darray, pio_write_darray
   use pio                 , only : pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling
   use pio                 , only : pio_def_var, pio_get_var, pio_put_var, PIO_INT
   use ioFileMod
   !$use omp_lib           , only : omp_set_num_threads

  implicit none
  private ! except

  public :: SetServices
  public :: SetVM

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize
  private :: cam_read_srfrest
  private :: cam_write_srfrest
  private :: cam_orbital_init
  private :: cam_orbital_update
  private :: cam_set_mesh_for_single_column
  private :: cam_pio_checkerr

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: flds_scalar_index_nextsw_cday = 0
  integer                      :: nthrds
  integer         , parameter  :: dbug_flag = 0
  type(cam_in_t)  , pointer    :: cam_in(:)
  type(cam_out_t) , pointer    :: cam_out(:)
  integer         , pointer    :: dof(:)              ! global index space decomposition
  character(len=256)           :: rsfilename_spec_cam ! Filename specifier for restart surface file
  character(*)    ,parameter   :: modName =  "(atm_comp_nuopc)"
  character(*)    ,parameter   :: u_FILE_u = &
       __FILE__

  logical                      :: dart_mode = .false.
  logical                      :: mediator_present

  character(len=CL)            :: orb_mode            ! attribute - orbital mode
  integer                      :: orb_iyear           ! attribute - orbital year
  integer                      :: orb_iyear_align     ! attribute - associated with model year
  real(R8)                     :: orb_obliq           ! attribute - obliquity in degrees
  real(R8)                     :: orb_mvelp           ! attribute - moving vernal equinox longitude
  real(R8)                     :: orb_eccen           ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  real(R8) , parameter         :: grid_tol = 1.e-2_r8 ! tolerance for calculated lat/lon vs read in

  type(ESMF_Mesh)  :: model_mesh     ! model_mesh
  type(ESMF_Clock) :: model_clock    ! model_clock

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS

    ! the NUOPC gcomp component will register the generic methods

    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD version v03

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetServices

  !===============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! intput/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer           :: n
    integer           :: localpet
    character(len=CL) :: cvalue
    character(len=CL) :: logmsg
    logical           :: isPresent, isSet
    integer           :: shrlogunit          ! original log unit
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VmGet(vm, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, localpet==0, iulog, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_log_setLogUnit (iulog)

    !----------------------------------------------------------------------------
    ! advertise import/export fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    ! read mediator fields namelists
    call read_surface_fields_namelists()

    call NUOPC_CompAttributeGet(gcomp, name="mediator_present", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read (cvalue,*) mediator_present
       if (mediator_present) then
          call advertise_fields(gcomp, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    else
       call shr_sys_abort(subname//'Need to set attribute mediator_present')
    endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if
  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    use ESMF, only : ESMF_VMGet

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: currTime                          ! Current time
    type(ESMF_Time)         :: startTime                         ! Start time
    type(ESMF_Time)         :: stopTime                          ! Stop time
    type(ESMF_Time)         :: refTime                           ! Ref time
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_CalKind_Flag) :: esmf_caltype                      ! esmf calendar type
    type(ESMF_DistGrid)     :: distGrid
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(R8), pointer       :: ownedElemCoords(:)
    real(r8), pointer       :: lat(:), latMesh(:)
    real(r8), pointer       :: lon(:), lonMesh(:)
    real(r8)                :: lats(pcols)                       ! array of chunk latitudes
    real(r8)                :: lons(pcols)                       ! array of chunk longitude
    integer                 :: hdim1_d, hdim2_d                  ! dims of rect horizontal grid data (If 1D data struct, hdim2_d==1)
    integer                 :: ncols                             ! number of local columns
    integer                 :: start_ymd                         ! Start date (YYYYMMDD)
    integer                 :: start_tod                         ! Start time of day (sec)
    integer                 :: curr_ymd                          ! Start date (YYYYMMDD)
    integer                 :: curr_tod                          ! Start time of day (sec)
    integer                 :: stop_ymd                          ! Stop date (YYYYMMDD)
    integer                 :: stop_tod                          ! Stop time of day (sec)
    integer                 :: ref_ymd                           ! Reference date (YYYYMMDD)
    integer                 :: ref_tod                           ! Reference time of day (sec)
    character(len=cs)       :: calendar                          ! Calendar type
    integer                 :: dtime                             ! time step increment (sec)
    integer                 :: atm_cpl_dt                        ! driver atm coupling time step
    integer                 :: nstep                             ! CAM nstep
    real(r8)                :: caldayp1                          ! CAM calendar day for for next cam time step
    integer                 :: yy,mm,dd                          ! Temporaries for time query
    logical                 :: perpetual_run                     ! If in perpetual mode or not
    integer                 :: perpetual_ymd                     ! Perpetual date (YYYYMMDD)
    character(CL)           :: cvalue
    character(ESMF_MAXSTR)  :: convCIM, purpComp
    integer                 :: lsize                             ! local size ofarrays
    integer                 :: n,c,g,i,j                         ! indices
    character(len=cs)       :: start_type                        ! infodata start type
    character(len=cl)       :: caseid                            ! case ID
    character(len=cl)       :: ctitle                            ! case title
    character(len=cl)       :: model_doi_url                     ! DOI for CESM model run
    logical                 :: aqua_planet                       ! Flag to run model in "aqua planet" mode
    logical                 :: brnch_retain_casename             ! true => branch run has same caseid as run being branched from
    logical                 :: single_column = .false.
    character(len=cl)       :: single_column_lnd_domainfile
    real(r8)                :: scol_lon
    real(r8)                :: scol_lat
    real(r8)                :: scol_spval
    real(r8)                :: eccen
    real(r8)                :: obliqr
    real(r8)                :: lambm0
    real(r8)                :: mvelpp
    !character(len=cl)      :: atm_resume_all_inst(num_inst_atm) ! atm resume file
    integer                 :: lbnum
    character(CS)           :: inst_name
    integer                 :: inst_index
    character(CS)           :: inst_suffix
    integer                 :: lmpicom
    logical                 :: isPresent, isSet
    character(len=512)      :: diro
    character(len=512)      :: logfile
    integer                 :: compid                            ! component id
    integer                 :: localPet, localPeCount
    logical                 :: initial_run                       ! startup mode which only requires a minimal initial file
    logical                 :: restart_run                       ! continue a previous run; requires a restart file
    logical                 :: branch_run                        ! branch from a previous run; requires a restart file
    character(len=CL)       :: tempc1,tempc2
    integer                 :: shrlogunit          ! original log unit
    real(r8)        , parameter :: radtodeg = 180.0_r8/shr_const_pi
    integer         , parameter :: aqua_perpetual_ymd = 321
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    call shr_log_setLogUnit (iulog)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, localpet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, pet=localPet, peCount=localPeCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if(localPeCount == 1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthrds
    else
       nthrds = localPeCount
    endif

!$  call omp_set_num_threads(nthrds)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------
    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    inst_name = 'ATM'//inst_suffix
    ! Set filename specifier for restart surface file
    ! (%c=caseid, $y=year, $m=month, $d=day, $s=seconds in day)
    rsfilename_spec_cam = '%c.cam' // trim(inst_suffix) // '.rs.%y-%m-%d-%s.nc'

    !----------------------------------------------------------------------------
    ! initialize cam mpi (needed for masterproc below)
    !----------------------------------------------------------------------------

    call spmdinit(lmpicom)

    !----------------------
    ! Initialize cam - needed in realize phase to get grid information
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    read(cvalue,*) compid
    call cam_instance_init(compid, inst_name, inst_index, inst_suffix)

    !----------------------
    ! Initialize cam - needed in realize phase to get grid information
    !----------------------

    if (masterproc) then
       write(iulog,format) "CAM atm model initialization"
    end if

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif

    !----------------------
    ! Obtain and load orbital values
    !----------------------

    call cam_orbital_init(gcomp, iulog, masterproc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call cam_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)

    !----------------------
    ! Obtain attributes
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle=caseid

    ! starting info
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_type
    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    ! single column input
    call NUOPC_CompAttributeGet(gcomp, name='scol_lon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lon
    call NUOPC_CompAttributeGet(gcomp, name='scol_lat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_lat
    call NUOPC_CompAttributeGet(gcomp, name='scol_spval', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scol_spval

    ! For single column mode in cam need to have a valid single_column_lnd_domainfile for the mask
    call NUOPC_CompAttributeGet(gcomp, name='single_column_lnd_domainfile', value=single_column_lnd_domainfile, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (scol_lon > scol_spval .and. scol_lat > scol_spval) then
       if (trim(single_column_lnd_domainfile) /= 'UNSET') then
          single_column = .true.
       else
          call shr_sys_abort('single_column_lnd_domainfile cannot be null for single column mode')
       end if
    else
       single_column = .false.
    end if

    ! aqua planet input
    call NUOPC_CompAttributeGet(gcomp, name='aqua_planet', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aqua_planet

    ! perpetual input
    call NUOPC_CompAttributeGet(gcomp, name='perpetual', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) perpetual_run
    call NUOPC_CompAttributeGet(gcomp, name='perpetual_ymd', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) perpetual_ymd

    ! TODO: query the config attributes for the number of instances - ASSUMES multi-driver

    ! TODO: must obtain model_doi_url from gcomp - for now hardwire to 'not_set'
    model_doi_url = 'not_set'

    ! Initialize CAM, allocate cam_in and cam_out and determine
    ! atm decomposition (needed to initialize gsmap)
    ! for an initial run, cam_in and cam_out are allocated in cam_init
    ! for a restart/branch run, cam_in and cam_out are allocated in restart
    !
    !TODO: the following strings must not be hard-wired - must have module variables
    ! like seq_infodata_start_type_type - maybe another entry in flds_mod?

    initial_run = .false.
    restart_run = .false.
    branch_run  = .false.
    if (trim(start_type) == trim('startup')) then
       initial_run = .true.
    else if (trim(start_type) == trim('continue') ) then
       restart_run = .true.
    else if (trim(start_type) == trim('branch')) then
       branch_run = .true.
    else
       call shr_sys_abort( subname//' ERROR: unknown start_type' )
    end if

    ! DART always starts up as an initial run.
    call NUOPC_CompAttributeGet(gcomp, name='data_assimilation_atm', value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) dart_mode
    end if
    if (dart_mode) then
       initial_run = .true.
       restart_run = .false.
       branch_run  = .false.
    end if

    ! Get properties from clock
    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    ! Initialize module orbital values and update orbital
    call cam_orbital_init(gcomp, iulog, masterproc, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call cam_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize CAM
    if (aqua_planet) then
       perpetual_run = .true.
       perpetual_ymd = aqua_perpetual_ymd
    end if

    call cam_init( &
         caseid=caseid, ctitle=ctitle, model_doi_url=model_doi_url, &
         initial_run_in=initial_run, restart_run_in=restart_run, &
         branch_run_in=branch_run, post_assim_in=dart_mode, &
         calendar=calendar, brnch_retain_casename=brnch_retain_casename, aqua_planet=aqua_planet, &
         single_column=single_column, scmlat=scol_lat, scmlon=scol_lon, &
         eccen=eccen, obliqr=obliqr, lambm0=lambm0, mvelpp=mvelpp,  &
         perpetual_run=perpetual_run, perpetual_ymd=perpetual_ymd, &
         dtime=dtime, start_ymd=start_ymd, start_tod=start_tod, ref_ymd=ref_ymd, ref_tod=ref_tod, &
         stop_ymd=stop_ymd, stop_tod=stop_tod, curr_ymd=curr_ymd, curr_tod=curr_tod, &
         cam_out=cam_out,  cam_in=cam_in)

    if (mediator_present) then

       if (single_column) then

          call cam_set_mesh_for_single_column(scol_lon, scol_lat, model_mesh, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(dof(1))
          dof(1) = 1

       else

          ! generate the dof
          lsize = 0
          do c = begchunk, endchunk
             do i = 1, get_ncols_p(c)
                lsize = lsize + 1
             end do
          end do
          allocate(dof(lsize))
          n = 0
          do c = begchunk, endchunk
             do i = 1, get_ncols_p(c)
                n = n+1
                dof(n) = get_gcol_p(c,i)
             end do
          end do

          ! create distGrid from global index array
          DistGrid = ESMF_DistGridCreate(arbSeqIndexList=dof, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! read in the mesh
          call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          model_mesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
               elementDistgrid=Distgrid, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (masterproc) then
             write(iulog,*)'mesh file for cam domain is ',trim(cvalue)
          end if

          ! obtain mesh lats and lons
          call ESMF_MeshGet(model_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (numOwnedElements /= lsize) then
             write(tempc1,'(i10)') numOwnedElements
             write(tempc2,'(i10)') lsize
             call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements "// trim(tempc1) // &
                  " not equal to local size "// trim(tempc2), ESMF_LOGMSG_INFO, rc=rc)
             rc = ESMF_FAILURE
             return
          end if
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          allocate(lonMesh(lsize), latMesh(lsize))
          call ESMF_MeshGet(model_mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             lonMesh(n) = ownedElemCoords(2*n-1)
             latMesh(n) = ownedElemCoords(2*n)
          end do

          ! obtain internally generated cam lats and lons
          allocate(lon(lsize)); lon(:) = 0._r8
          allocate(lat(lsize)); lat(:) = 0._r8
          n=0
          do c = begchunk, endchunk
             ncols = get_ncols_p(c)
             ! latitudes and longitudes returned in radians
             call get_rlat_all_p(c, ncols, lats)
             call get_rlon_all_p(c, ncols, lons)
             do i=1,ncols
                n = n+1
                lat(n) = lats(i)*radtodeg
                lon(n) = lons(i)*radtodeg
             end do
          end do

          ! error check differences between internally generated lons and those read in
          do n = 1,lsize
             if (abs(lonMesh(n) - lon(n)) > grid_tol .and. abs(lonMesh(n) - lon(n)) /= 360._r8) then
                write(6,100)n,lon(n),lonMesh(n), abs(lonMesh(n)-lon(n))
100             format('ERROR: CAM n, lonmesh(n), lon(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
                call shr_sys_abort()
             end if
             if (abs(latMesh(n) - lat(n)) > grid_tol) then
                write(6,100)n,lat(n),latMesh(n), abs(latMesh(n)-lat(n))
101             format('ERROR: CAM n, latmesh(n), lat(n), diff_lat = ',i6,2(f21.13,3x),d21.5)
                call shr_sys_abort()
             end if
          end do

          ! deallocate memory
          deallocate(ownedElemCoords)
          deallocate(lon, lonMesh)
          deallocate(lat, latMesh)

       end if ! end of if single_column

       ! realize the actively coupled fields
       call realize_fields(gcomp, model_mesh, flds_scalar_name, flds_scalar_num, single_column, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Create model_clock as a module variable - needed for generating streams
       model_clock = clock

       ! Create cam export array and set the state scalars
       call export_fields( gcomp, model_mesh, model_clock, cam_out, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
       call State_SetScalar(dble(hdim1_d), flds_scalar_index_nx, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call State_SetScalar(dble(hdim2_d), flds_scalar_index_ny, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! diagnostics
       if (dbug_flag > 1) then
          call State_diagnose(exportState,subname//':ES',rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

    end if ! end of mediator_present if-block

    call shr_log_setLogUnit (shrlogunit)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if

  end subroutine InitializeRealize

  !===============================================================================
  subroutine DataInitialize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)                   :: clock
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Time)                    :: currTime      ! Current time
    type(ESMF_TimeInterval)            :: timeStep
    type(ESMF_Field)                   :: field
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    character(ESMF_MAXSTR)             :: fieldName
    integer                            :: n, fieldCount
    integer                            :: shrlogunit    ! original log unit
    integer(ESMF_KIND_I8)              :: stepno        ! time step
    integer                            :: atm_cpl_dt    ! driver atm coupling time step
    logical                            :: importDone    ! true => import data is valid
    logical                            :: atCorrectTime ! true => field is at correct time
    character(CL)                      :: cvalue
    character(len=*),parameter         :: subname=trim(modName)//':(DataInitialize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    end if

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_DataInitialize:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call log_clock_advance(clock, 'CAM', iulog, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------------------------------
    if (mediator_present) then
    !---------------------------------------------------------------

       ! Determine if all the import state has been initialized
       ! And if not initialized, then return

       call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       importDone = .true.
       do n=1, fieldCount
          call ESMF_StateGet(importState, itemName=fieldNameList(n), field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          atCorrectTime = NUOPC_IsAtTime(field, currTime, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (.not. atCorrectTime) then
             call ESMF_LogWrite("CAM - Initialize-Data-Dependency NOT YET SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             importDone = .false.
             exit  ! break out of the loop when first not satisfied found
          end if
       end do
       deallocate(fieldNameList)

       ! *** Import state has not been initialized - RETURN ****

       if (.not. importDone) then
          ! Simply return if the import has not been initialized
          call ESMF_LogWrite("CAM - Initialize-Data-Dependency Returning to mediator without doing tphysbc", &
               ESMF_LOGMSG_INFO, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          RETURN
       end if

       ! *** Import state has been initialized - continue with tphysbc ***

       call ESMF_LogWrite("CAM - Initialize-Data-Dependency doing tphysbc", ESMF_LOGMSG_INFO, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! get the current step number and coupling interval
       call ESMF_ClockGet( clock, TimeStep=timeStep, advanceCount=stepno, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet( timeStep, s=atm_cpl_dt, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! For initial run, unpack the import state, run cam radiation/clouds and return
       ! For restart run, read the import state from the restart and run radiation/clouds and return

       ! Note - cam_run1 is called only for the purposes of finishing the
       ! flux averaged calculation to compute cam-out
       ! Note - cam_run1 is called on restart only to have cam internal state consistent with the
       ! cam_out state sent to the coupler

       if (stepno == 0) then
          call import_fields( gcomp, cam_in, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call cam_run1 ( cam_in, cam_out )
          call export_fields( gcomp, model_mesh, model_clock, cam_out, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call cam_read_srfrest( gcomp, clock, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call import_fields( gcomp, cam_in, restart_init=.true., rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call cam_run1 ( cam_in, cam_out )
          call export_fields( gcomp, model_mesh, model_clock, cam_out, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Compute time of next radiation computation
       call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! diagnostics
       if (dbug_flag > 1) then
          call State_diagnose(exportState,subname//':ES',rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! CAM data is now fully initialized

       call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       do n=1, fieldCount
          call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end do
       deallocate(fieldNameList)

       ! check whether all Fields in the exportState are "Updated"
       if (NUOPC_IsUpdated(exportState)) then
          call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_LogWrite("CAM - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    !---------------------------------------------------------------
    else  ! mediator is not present
    !---------------------------------------------------------------

       call cam_run1 ( cam_in, cam_out )

       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    ! End redirection of share output to cam log
    call shr_log_setLogUnit (shrlogunit)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_DataInitialize:end::',lbnum)
    endif
#endif

    if (dbug_flag > 5) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)
    end if

  end subroutine DataInitialize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    use ESMF, only : ESMF_GridCompGet, esmf_vmget, esmf_vm
    ! Run CAM

    ! Input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: time
    type(ESMF_Time)         :: currTime    ! Current time
    type(ESMF_Time)         :: nextTime    ! Next timestep time
    type(ESMF_TimeInterval) :: timeStep    ! Clock, time-step
    type(ESMF_State)        :: importState
    type(ESMF_State)        :: exportState
    character(CL)           :: cvalue
    integer                 :: shrlogunit  ! original log unit
    character(CL)           :: case_name   ! case name
    real(r8)                :: eccen
    real(r8)                :: obliqr
    real(r8)                :: lambm0
    real(r8)                :: mvelpp
    logical                 :: dosend      ! true => send data back to driver
    integer                 :: dtime       ! time step increment (sec)
    integer                 :: ymd_sync    ! Sync ymd
    integer                 :: yr_sync     ! Sync current year
    integer                 :: mon_sync    ! Sync current month
    integer                 :: day_sync    ! Sync current day
    integer                 :: tod_sync    ! Sync current time of day (sec)
    integer                 :: ymd         ! CAM current date (YYYYMMDD)
    integer                 :: yr          ! CAM current year
    integer                 :: mon         ! CAM current month
    integer                 :: day         ! CAM current day
    integer                 :: tod         ! CAM current time of day (sec)
    logical                 :: rstwr       ! .true. ==> write restart file before returning
    logical                 :: nlend       ! Flag signaling last time-step
    integer                 :: lbnum
    integer                 :: localPet, localPeCount
    logical                 :: first_time = .true.
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

!$  call omp_set_num_threads(nthrds)

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call log_clock_advance(clock, 'CAM', iulog, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !--------------------------------
    ! Determine current time
    !--------------------------------

    call ESMF_ClockGet( clock, currTime=currTime)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Update and load orbital parameters
    !----------------------

    if (trim(orb_mode) == trim(orb_variable_year) .or. first_time) then
       call cam_orbital_update(clock, iulog, masterproc, eccen, obliqr, lambm0, mvelpp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)
    end if
    first_time = .false.

    !--------------------------------
    ! Run cam
    !--------------------------------

    ! Unpack import state
    if (mediator_present) then
       call t_startf ('CAM_import')
       call State_diagnose(importState, string=subname//':IS', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fields( gcomp, cam_in, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf  ('CAM_import')
    end if

    dosend = .false.
    do while (.not. dosend)

       ! TODO: This is currently hard-wired - is there a better way for nuopc?
       ! Need to not return when nstep = 0 and return when nstep = 1
       ! Note that the model clock is updated at the end of the time step not at the beginning
       if (get_nstep() > 0) then
          dosend = .true.
       end if

       ! Determine if time to write restart

       call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          rstwr = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          rstwr = .false.
       endif

       ! Determine if time to stop

       call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          nlend = .true.
       else
          nlend = .false.
       endif

       ! Run CAM (run2, run3, run4)

       call t_startf ('CAM_run2')
       call cam_run2( cam_out, cam_in )
       call t_stopf  ('CAM_run2')

       call t_startf ('CAM_run3')
       call cam_run3( cam_out )
       call t_stopf  ('CAM_run3')

       call t_startf ('CAM_run4')
       call cam_run4( cam_out, cam_in, rstwr, nlend, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf  ('CAM_run4')

       ! Advance cam time step

       call t_startf ('CAM_adv_timestep')
       call advance_timestep()
       call t_stopf  ('CAM_adv_timestep')

       ! Run cam radiation/clouds (run1)

       call t_startf ('CAM_run1')
       call cam_run1 ( cam_in, cam_out )
       call t_stopf  ('CAM_run1')

    end do

    if (mediator_present) then
       ! Set export fields
       call t_startf ('CAM_export')
       call export_fields( gcomp, model_mesh, model_clock, cam_out, rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf ('CAM_export')

       ! Set the coupling scalars
       ! Return time of next radiation calculation - albedos will need to be
       ! calculated by each surface model at this time
       call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
            flds_scalar_name, flds_scalar_num, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! diagnostics
       if (dbug_flag > 1) then
          call State_diagnose(exportState, string=subname//':ES',rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (masterproc) then
             call log_clock_advance(clock, 'CAM', iulog, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       endif

       ! Write merged surface data restart file if appropriate
       if (rstwr) then
          call cam_write_srfrest( gcomp, &
               yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    else

       ! if there is no mediator, then write the clock info to a driver restart file
       if (rstwr) then
          call cam_write_clockrest( clock, yr_sync, mon_sync, day_sync, tod_sync, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    end if

    ! Check for consistency of internal cam clock with master sync clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! CAM is actually 2 coupling intervals (or physics time steps) ahead of the driver clock
    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) )then
       write(iulog,*)' cam ymd=',ymd     ,'  cam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( subname//': CAM clock is not in sync with master Sync Clock' )
    end if

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_ModelAdvance:end::',lbnum)
    endif
#endif

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_log_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelSetRunClock

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: shrlogunit            ! original log unit
    character(*), parameter :: F00   = "('(atm_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(atm_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)

    call cam_final( cam_out, cam_in )

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'CAM: end of main integration loop'
       write(iulog,F91)
    end if

    call shr_log_setLogUnit (shrlogunit)

  end subroutine ModelFinalize

  !===============================================================================
  subroutine cam_orbital_init(gcomp, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp) , intent(in)    :: gcomp
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(cam_orbital_init)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
             write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       orb_obliq = SHR_ORB_UNDEF_REAL
       orb_eccen = SHR_ORB_UNDEF_REAL
       orb_mvelp = SHR_ORB_UNDEF_REAL
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       orb_iyear = SHR_ORB_UNDEF_INT
       orb_iyear_align = SHR_ORB_UNDEF_INT
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
           orb_obliq == SHR_ORB_UNDEF_REAL .or. &
           orb_mvelp == SHR_ORB_UNDEF_REAL) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif

  end subroutine cam_orbital_init

  !===============================================================================
  subroutine cam_orbital_update(clock, logunit,  mastertask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: mastertask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical, save     :: logprint = .true.
    character(len=*) , parameter :: subname = "(cam_orbital_update)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
    else
       orb_year = orb_iyear
    end if
    if(.not. (logprint .and. mastertask)) then
       logprint = .false.
    endif

    eccen = orb_eccen

    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, logprint)
    logprint = .false.
    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine cam_orbital_update

  !===============================================================================
  subroutine cam_read_srfrest( gcomp, clock, rc )

    ! input/output variables
    type(ESMF_GridComp)             :: gcomp
    type(ESMF_Clock), intent(inout) :: clock
    integer         , intent(out)   :: rc

    ! local variables
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Field)                   :: lfield
    integer                            :: lrank
    integer                            :: rcode        ! return error code
    integer                            :: nf,n
    type(file_desc_t)                  :: file
    type(io_desc_t)                    :: iodesc
    integer                            :: fieldCount
    character(ESMF_MAXSTR),allocatable :: fieldNameList(:)
    type(var_desc_t)                   :: varid
    real(r8), pointer                  :: fldptr(:)
    real(r8), pointer                  :: tmpptr(:)
    real(r8), pointer                  :: fldptr2d(:,:)
    type(ESMF_Time)                    :: currTime      ! time at previous interval
    integer                            :: yr_spec       ! Current year
    integer                            :: mon_spec      ! Current month
    integer                            :: day_spec      ! Current day
    integer                            :: sec_spec      ! Current time of day (sec)
    character(len=256)                 :: fname_srf_cam ! surface restart filename
    character(len=256)                 :: pname_srf_cam ! surface restart full pathname
    character(len=PIO_MAX_NAME)        :: varname
    integer                            :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fieldds
    integer                            :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fieldds
    integer                            :: lsize
    character(len=8)                   :: cvalue
    integer                            :: nloop
    character(len=4)                   :: prefix
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! ------------------------------
    ! Get surface restart dataset
    ! ------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_spec, mm=mon_spec, dd=day_spec, s=sec_spec, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=cam_initfiles_get_caseid(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(cam_initfiles_get_restdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)

    ! ------------------------------
    ! Open restart file
    ! ------------------------------

    call cam_pio_openfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)
    call pio_seterrorhandling(File, pio_bcast_error)

    ! ------------------------------
    ! Read in import and export fields
    ! ------------------------------

    do nloop = 1,2

       if (nloop == 1) then
          prefix = 'x2a_' ! import fields
          call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldnameList(fieldCount))
          call ESMF_StateGet(importState, itemNameList=fieldnameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          prefix = 'a2x_' ! export fields
          call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldnameList(fieldCount))
          call ESMF_StateGet(exportState, itemNameList=fieldnameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Loop over fields in import or export state
       do nf = 1,fieldCount

          if (trim(fieldnameList(nf)) == flds_scalar_name) CYCLE

          ! Determine dimension of field
          if (nloop == 1) then
             call ESMF_StateGet(importState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_StateGet(exportState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (lrank == 1) then

             call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             varname = trim(prefix)//trim(fieldnameList(nf))
             rcode = pio_inq_varid(File,trim(varname) ,varid)
             if (rcode == pio_noerr) then
                call pio_read_darray(File, varid, iodesc, fldptr, rcode)
             else
                if (masterproc) then
                   write(iulog,*)'cam_read_srfrest warning: field ',trim(varname),' is not on restart file'
                   write(iulog,*)'for backwards compatibility will set it to 0'
                end if
                fldptr(:) = 0._r8
             end if

          else if (lrank == 2) then

             ! There is an output variable for each element of the undistributed dimension
             call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             if (gridToFieldMap(1) == 1) then
                lsize = size(fldptr2d, dim=1)
             else if (gridToFieldMap(1) == 2) then
                lsize = size(fldptr2d, dim=2)
             end if

             allocate(tmpptr(lsize))
             do n = 1,ungriddedUBound(1)
                write(cvalue,'(i0)') n
                varname = trim(prefix)//trim(fieldnameList(nf))//trim(cvalue)
                rcode = pio_inq_varid(File,trim(varname) ,varid)

                if (rcode == pio_noerr) then
                   call pio_read_darray(File, varid, iodesc, tmpptr, rcode)
                else
                   if (masterproc) then
                      write(iulog,*)'cam_read_srfrest warning: field ',trim(varname),' is not on restart file'
                      write(iulog,*)'for backwards compatibility will set it to 0'
                   end if
                   tmpptr(:) = 0._r8
                end if
                if (gridToFieldMap(1) == 1) then
                   fldptr2d(:,n) = tmpptr(:)
                else if (gridToFieldMap(1) == 2) then
                   fldptr2d(n,:) = tmpptr(:)
                end if
             end do
             deallocate(tmpptr)

          end if ! end lrank if block
       end do
       deallocate(fieldnameList)
    end do

    ! ------------------------------
    ! Close file
    ! ------------------------------

    call pio_seterrorhandling(File, pio_internal_error)
    call pio_freedecomp(File, iodesc)
    call cam_pio_closefile(File)

  end subroutine cam_read_srfrest

  !===========================================================================================
  subroutine cam_write_srfrest( gcomp, yr_spec, mon_spec, day_spec, sec_spec, rc )

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    integer , intent(in)  :: yr_spec  ! Simulation year
    integer , intent(in)  :: mon_spec ! Simulation month
    integer , intent(in)  :: day_spec ! Simulation day
    integer , intent(in)  :: sec_spec ! Seconds into current simulation day
    integer , intent(out) :: rc       ! error code

    ! Local variables
    type(ESMF_State)                   :: importState, exportState
    type(ESMF_Field)                   :: lField
    integer                            :: lrank
    integer                            :: rcode        ! return error code
    integer                            :: dimid(1), nf, n
    type(file_desc_t)                  :: file
    type(io_desc_t)                    :: iodesc
    integer                            :: fieldCount
    character(ESMF_MAXSTR),allocatable :: fieldnameList(:)
    type(var_desc_t)                   :: varid
    real(r8), pointer                  :: fldptr1d(:)
    real(r8), pointer                  :: fldptr2d(:,:)
    character(len=PIO_MAX_NAME)        :: varname
    character(len=256)                 :: fname_srf_cam      ! surface restart filename
    character(len=8)                   :: cvalue
    integer                            :: nloop
    character(len=4)                   :: prefix
    integer                            :: ungriddedUBound(1) ! currently the size must equal 1 for rank 2 fieldds
    integer                            :: gridToFieldMap(1)  ! currently the size must equal 1 for rank 2 fieldds
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! ----------------------
    ! Get import and export states
    ! ----------------------

    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ----------------------
    ! Open surface restart dataset
    ! ----------------------

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )

    call cam_pio_createfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)

    ! ----------------------
    ! Define dimensions
    ! ----------------------

    rcode = pio_def_dim(File, 'x2a_nx', ngcols, dimid(1))
    rcode = pio_def_dim(File, 'a2x_nx', ngcols, dimid(1))

    ! ----------------------
    ! Define import and export variable ids
    ! ----------------------

    do nloop = 1,2

       if (nloop == 1) then
          prefix = 'x2a_' ! import fields
          call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          prefix = 'a2x_' ! export fields
          call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       do nf = 1,fieldCount

          if (trim(fieldNameList(nf)) == flds_scalar_name) CYCLE

          if (nloop == 1) then
             call ESMF_StateGet(importState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_StateGet(exportState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (lrank == 1) then

             varname = trim(prefix)//trim(fieldNameList(nf))
             rcode = pio_def_var(File,trim(varname), PIO_DOUBLE, dimid, varid)
             rcode = pio_put_att(File, varid, "_fillvalue", fillvalue)

          else if (lrank == 2) then

             ! Determine the size of the ungridded dimension and the
             ! index where the undistributed dimension is located
             call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, gridToFieldMap=gridToFieldMap, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             ! Output for each ungriddedUbound index
             do n = 1,ungriddedUBound(1)
                write(cvalue,'(i0)') n
                varname = trim(prefix)//trim(fieldNameList(nf))//trim(cvalue)
                rcode = pio_def_var(File,trim(varname), PIO_DOUBLE, dimid, varid)
                rcode = pio_put_att(File, varid, "_fillvalue", fillvalue)
             end do

          end if ! end if-block over rank size

       end do ! end loop over import or export fieldsfields
       deallocate(fieldNameList)
    end do

    ! ----------------------
    ! End definition phase
    ! ----------------------

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    ! ----------------------
    ! Write the restart data for the import fields and export fields
    ! ----------------------

    do nloop = 1,2

       if (nloop == 1) then
          prefix = 'x2a_' ! import fields
          call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          prefix = 'a2x_' ! export fields
          call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(fieldNameList(fieldCount))
          call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       do nf = 1,fieldCount

          if (trim(fieldNameList(nf)) == flds_scalar_name) CYCLE

          if (nloop == 1) then
             call ESMF_StateGet(importState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call ESMF_StateGet(exportState, itemName=trim(fieldnameList(nf)), field=lfield, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call ESMF_FieldGet(lfield, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (lrank == 1) then

             varname = trim(prefix)//trim(fieldNameList(nf))
             rcode = pio_inq_varid(File, trim(varname), varid)
             call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call pio_write_darray(File, varid, iodesc, fldptr1d, rcode)

          else if (lrank == 2) then

             ! There is an output variable for each element of the undistributed dimension
             call ESMF_FieldGet(lfield, ungriddedUBound=ungriddedUBound, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             do n = 1,ungriddedUBound(1)
                write(cvalue,'(i0)') n
                varname = trim(prefix)//trim(fieldNameList(nf))//trim(cvalue)
                rcode = pio_inq_varid(File, trim(varname), varid)
                if (gridToFieldMap(1) == 1) then
                   call pio_write_darray(File, varid, iodesc, fldptr2d(:,n), rcode, fillval=fillvalue)
                else if (gridToFieldMap(1) == 2) then
                   call pio_write_darray(File, varid, iodesc, fldptr2d(n,:), rcode, fillval=fillvalue)
                end if
             end do

          end if
       end do ! end loop over import or export fields
       deallocate(fieldNameList)

    end do ! end of nloop

    ! ----------------------
    ! close the file
    ! ----------------------

    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(File)

  end subroutine cam_write_srfrest

  !===============================================================================
  subroutine cam_write_clockrest( clock, yr_spec, mon_spec, day_spec, sec_spec, rc )

    ! When there is no mediator, the driver needs to have restart information to start up
    ! This routine writes this out and the driver reads it back in on a restart run

    ! Arguments
    type(ESMF_Clock) , intent(in)  :: clock
    integer          , intent(in)  :: yr_spec  ! Simulation year
    integer          , intent(in)  :: mon_spec ! Simulation month
    integer          , intent(in)  :: day_spec ! Simulation day
    integer          , intent(in)  :: sec_spec ! Seconds into current simulation day
    integer          , intent(out) :: rc       ! error code

    ! Local variables
    type(ESMF_Time)   :: startTime
    type(ESMF_Time)   :: currTime
    type(ESMF_Time)   :: nextTime
    integer           :: unitn
    type(file_desc_t) :: File
    integer           :: start_ymd
    integer           :: start_tod
    integer           :: curr_ymd
    integer           :: curr_tod
    integer           :: yy,mm,dd  ! Temporaries for time query
    type(var_desc_t)  :: varid_start_ymd
    type(var_desc_t)  :: varid_start_tod
    type(var_desc_t)  :: varid_curr_ymd
    type(var_desc_t)  :: varid_curr_tod
    integer           :: rcode
    character(ESMF_MAXSTR)  :: restart_pfile
    character(ESMF_MAXSTR)  :: restart_file
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get properties from clock
    call ESMF_ClockGet( clock, startTime=startTime, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( nextTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
   !call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    ! Open clock info restart dataset
    restart_file = interpret_filename_spec( '%c.cpl.r.%y-%m-%d-%s.nc', &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )

    if (masterproc) then
       write(iulog,*) " In this configuration, there is no mediator"
       write(iulog,*) " Normally, the mediator restart file provides the restart time info"
       write(iulog,*) " In this case, CAM will create the rpointer.cpl and cpl restart file"
       write(iulog,*) " containing this information"
       write(iulog,*) " writing rpointer file for driver clock info, rpointer.cpl"
       write(iulog,*) " writing restart clock info for driver= "//trim(restart_file)
       open(newunit=unitn, file='rpointer.cpl', form='FORMATTED')
       write(unitn,'(a)') trim(restart_file)
       close(unitn)
    endif

    call cam_pio_createfile(File, trim(restart_file), 0)
    rcode = pio_def_var(File, 'start_ymd', PIO_INT, varid_start_ymd)
    rcode = pio_def_var(File, 'start_tod', PIO_INT, varid_start_tod)
    rcode = pio_def_var(File, 'curr_ymd' , PIO_INT, varid_curr_ymd)
    rcode = pio_def_var(File, 'curr_tod' , PIO_INT, varid_curr_tod)
    rcode = pio_enddef(File)
    rcode = pio_put_var(File, varid_start_ymd, start_ymd)
    rcode = pio_put_var(File, varid_start_tod, start_tod)
    rcode = pio_put_var(File, varid_curr_ymd, curr_ymd)
    rcode = pio_put_var(File, varid_curr_tod, curr_tod)
    call cam_pio_closefile(File)

  end subroutine cam_write_clockrest

  !===============================================================================
  subroutine cam_set_mesh_for_single_column(scol_lon, scol_lat, mesh, rc)

    ! Generate a mesh for single column
    use netcdf

    ! input/output variables
    real(r8)        , intent(in)  :: scol_lon
    real(r8)        , intent(in)  :: scol_lat
    type(ESMF_Mesh) , intent(out) :: mesh
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_Grid)        :: lgrid
    integer                :: maxIndex(2)
    real(r8)               :: mincornerCoord(2)
    real(r8)               :: maxcornerCoord(2)
    character(len=*), parameter :: subname= ' (lnd_set_mesh_for_single_column) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Use center and come up with arbitrary area delta lon and lat = .1 degree
    maxIndex(1)       = 1                ! number of lons
    maxIndex(2)       = 1                ! number of lats
    mincornerCoord(1) = scol_lon - .1_r8 ! min lon
    mincornerCoord(2) = scol_lat - .1_r8 ! min lat
    maxcornerCoord(1) = scol_lon + .1_r8 ! max lon
    maxcornerCoord(2) = scol_lat + .1_r8 ! max lat

    ! create the ESMF grid
    lgrid = ESMF_GridCreateNoPeriDimUfrm (maxindex=maxindex, &
         mincornercoord=mincornercoord, maxcornercoord= maxcornercoord, &
         staggerloclist=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create the mesh from the lgrid
    mesh = ESMF_MeshCreate(lgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine cam_set_mesh_for_single_column

  !===============================================================================
  subroutine cam_pio_checkerr(ierror, description)
    use pio, only : PIO_NOERR
    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    if (ierror /= PIO_NOERR) then
       write (*,'(6a)') 'ERROR ', trim(description)
       call shr_sys_abort()
    endif
  end subroutine cam_pio_checkerr

end module atm_comp_nuopc
