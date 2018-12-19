module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CAM
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use med_constants_mod     , only : IN, R8, I8, CXX, CS,CL
  use shr_sys_mod           , only : shr_sys_abort
  use shr_log_mod           , only : shr_log_Unit
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_time_mod    , only : shr_nuopc_time_AlarmInit

  use atm_import_export     , only : advertise_fields, realize_fields
  use atm_import_export     , only : import_fields, export_fields
  use cam_instance          , only : cam_instance_init, inst_suffix, inst_index
  use cam_comp              , only : cam_init, cam_run1, cam_run2, cam_run3, cam_run4, cam_final
  use radiation             , only : radiation_nextsw_cday
  use camsrfexch            , only : cam_out_t, cam_in_t
  use cam_logfile           , only : iulog
  use spmd_utils            , only : spmdinit, masterproc, iam, mpicom
  use time_manager          , only : get_curr_calday, advance_timestep, get_curr_date, get_nstep, get_step_size
  use ioFileMod
  use perf_mod

  !TODO: add this to gcomp for the atm
  use seq_comm_mct, only : num_inst_atm

  implicit none
  private ! except

  public :: SetServices

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

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type fld_list_type
    character(len=128) :: stdname
  end type fld_list_type

  integer         , parameter :: dbug_flag = 6
  type(cam_in_t)  , pointer   :: cam_in(:)
  type(cam_out_t) , pointer   :: cam_out(:)
  integer         , pointer   :: dof(:)              ! global index space decomposition
  integer                     :: shrlogunit          ! original log unit
  integer                     :: shrloglev           ! original log level
  character(len=256)          :: rsfilename_spec_cam ! Filename specifier for restart surface file
  character(*)    ,parameter  :: modName =  "(atm_comp_nuopc)"
  character(*)    ,parameter  :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: dbrc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods

    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD version v03

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    use shr_nuopc_utils_mod, only : shr_nuopc_set_component_logging
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer :: dbrc
    integer :: n
    integer :: localpet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_VMGetCurrent(vm, rc=rc)
    call ESMF_VmGet(vm, localPet=localPet, rc=rc)
    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_nuopc_set_component_logging(gcomp, localpet==0, iulog, shrlogunit, shrloglev)

    call advertise_fields(gcomp, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    use ppgrid          , only : pcols, begchunk, endchunk
    use phys_grid       , only : get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
    use dyn_grid        , only : get_horiz_grid_dim_d
    use cam_control_mod , only : dart_mode, cam_ctrl_set_orbit
    use shr_const_mod   , only : shr_const_pi

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Time)         :: currTime                          ! Current time
    type(ESMF_Time)         :: startTime                         ! Start time
    type(ESMF_Time)         :: stopTime                          ! Stop time
    type(ESMF_Time)         :: refTime                           ! Ref time
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Calendar)     :: esmf_calendar                     ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype                      ! esmf calendar type
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh, EmeshTemp
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
    logical                 :: single_column
    real(r8)                :: scmlat
    real(r8)                :: scmlon
    real(r8)                :: eccen
    real(r8)                :: obliqr
    real(r8)                :: lambm0
    real(r8)                :: mvelpp
    logical                 :: dart_mode_in
    character(len=cl)       :: atm_resume_all_inst(num_inst_atm) ! atm resume file
    integer                 :: lbnum
    character(CS)           :: inst_name
    integer                 :: inst_index
    character(CS)           :: inst_suffix
    integer                 :: lmpicom
    type(ESMF_VM)           :: vm
    logical                 :: isPresent
    character(len=512)      :: diro
    character(len=512)      :: logfile
    integer                 :: compid                            ! component id
    logical                 :: initial_run                       ! startup mode which only requires a minimal initial file
    logical                 :: restart_run                       ! continue a previous run; requires a restart file
    logical                 :: branch_run                        ! branch from a previous run; requires a restart file
    integer                 :: dbrc
    character(len=CL)       :: tempc1,tempc2
    real(r8)        , parameter :: radtodeg = 180.0_r8/shr_const_pi
    integer         , parameter :: aqua_perpetual_ymd = 321
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_setLogUnit (iulog)
    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------
    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       cvalue = inst_suffix(2:)
       read(cvalue, *) inst_index
    else
       inst_suffix = ''
       inst_index = 1
    end if
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
    ! Obtain orbital values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid

    ! load orbital parameters
    call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)

    !TODO: ctitle does not appear in the esm_AddAttributes in esm.F90
    ! just hard-wire from now - is this even needed?
    ! call NUOPC_CompAttributeGet(gcomp, name='ctitle', value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) ctitle
    ctitle='UNSET'

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) single_column

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_type

    call NUOPC_CompAttributeGet(gcomp, name='aqua_planet', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aqua_planet

    call NUOPC_CompAttributeGet(gcomp, name='perpetual', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) perpetual_run

    call NUOPC_CompAttributeGet(gcomp, name='perpetual_ymd', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) perpetual_ymd

    ! TODO: must obtain model_doi_url from gcomp - for now hardwire to 'not_set'
    model_doi_url = 'not_set'

    ! TODO: atm_resume_all_inst must be added to gcomp or infodata
    atm_resume_all_inst(:) = 'FALSE'

    dart_mode_in = .false.
    if (trim(atm_resume_all_inst(MIN(num_inst_atm,inst_index))) == 'TRUE') dart_mode_in = .true.

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

    ! Get properties from clock
    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    ! Initialize CAM
    if (aqua_planet) then
       perpetual_run = .true.
       perpetual_ymd = aqua_perpetual_ymd
    end if

    call cam_init( &
         caseid=caseid, &
         ctitle=ctitle, &
         model_doi_url=model_doi_url, &
         dart_mode=dart_mode_in, &
         initial_run_in=initial_run, &
         restart_run_in=restart_run, &
         branch_run_in=branch_run, &
         calendar=calendar, &
         brnch_retain_casename=brnch_retain_casename, &
         aqua_planet=aqua_planet, &
         single_column=single_column, &
         scmlat=scmlat, &
         scmlon=scmlon, &
         eccen=eccen, &
         obliqr=obliqr, &
         lambm0=lambm0, &
         mvelpp=mvelpp,  &
         perpetual_run=perpetual_run, &
         perpetual_ymd=perpetual_ymd, &
         dtime=dtime, &
         start_ymd=start_ymd, &
         start_tod=start_tod, &
         ref_ymd=ref_ymd, &
         ref_tod=ref_tod, &
         stop_ymd=stop_ymd, &
         stop_tod=stop_tod, &
         curr_ymd=curr_ymd, &
         curr_tod=curr_tod, &
         cam_out=cam_out, &
         cam_in=cam_in)

    !--------------------------------
    ! generate the mesh
    !--------------------------------

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    EMeshTemp = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       write(iulog,*)'mesh file for cam domain is ',trim(cvalue)
    end if

    ! recreate the mesh using the above distGrid
    EMesh = ESMF_MeshCreate(EMeshTemp, elementDistgrid=Distgrid, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! obtain mesh lats and lons
    call ESMF_MeshGet(Emesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
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
    call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
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
       if (abs(lonMesh(n) - lon(n)) > 1.e-12_r8) then
          write(6,100)n,lon(n),lonMesh(n), abs(lonMesh(n)-lon(n))
100       format('ERROR: CAM n, lonmesh(n), lon(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
       end if
       if (abs(latMesh(n) - lat(n)) > 1.e-12_r8) then
          write(6,100)n,lat(n),latMesh(n), abs(latMesh(n)-lat(n))
101       format('ERROR: CAM n, latmesh(n), lat(n), diff_lat = ',i6,2(f21.13,3x),d21.5)
       end if
    end do

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)

    !--------------------------------
    ! realize the actively coupled fields
    !--------------------------------

    call realize_fields(gcomp,  Emesh, rc)

    !--------------------------------
    ! Create cam export array and set the state scalars
    !--------------------------------

    call export_fields( gcomp, cam_out, rc=rc  )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
    call shr_nuopc_methods_State_SetScalar(dble(hdim1_d), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(hdim2_d), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug_flag > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "CAM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Community Atmosphere Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", "Community Atmosphere Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Atmosphere", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", TBD, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

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
    integer                            :: shrloglev     ! original log level
    integer(ESMF_KIND_I8)              :: stepno        ! time step
    integer                            :: dtime         ! time step increment (sec)
    integer                            :: atm_cpl_dt    ! driver atm coupling time step
    integer                            :: nstep         ! CAM nstep
    real(r8)                           :: caldayp1      ! CAM calendar day for for next cam time step
    real(r8)                           :: nextsw_cday   ! calendar of next atm shortwave
    logical                            :: importDone    ! true => import data is valid
    logical                            :: atCorrectTime ! true => field is at correct time
    character(CL)                      :: cvalue
    integer                            :: dbrc
    character(len=*),parameter         :: subname=trim(modName)//':(DataInitialize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (iulog)

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    endif

    !--------------------------------
    ! Determine if all the import state has been initialized
    ! And if not initialized, then return
    !--------------------------------

    call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    importDone = .true.
    do n=1, fieldCount
       call ESMF_StateGet(importState, itemName=fieldNameList(n), field=field, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       atCorrectTime = NUOPC_IsAtTime(field, currTime, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (.not. atCorrectTime) then
          call ESMF_LogWrite("CAM - Initialize-Data-Dependency NOT YET SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          importDone = .false.
          exit  ! break out of the loop when first not satisfied found
       end if
    end do
    deallocate(fieldNameList)

    !--------------------------------
    ! Import state has not been initialized - RETURN
    !--------------------------------

    if (.not. importDone) then
       ! Simply return if the import has not been initialized
       call ESMF_LogWrite("CAM - Initialize-Data-Dependency Returning to mediator without doing tphysbc", &
            ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       RETURN
    end if

    !--------------------------------
    ! Import state has been initialized - continue with tphysbc
    !--------------------------------

    call ESMF_LogWrite("CAM - Initialize-Data-Dependency doing tphysbc", ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! get the current step number and coupling interval
    !--------------------------------

    call ESMF_ClockGet( clock, TimeStep=timeStep, advanceCount=stepno, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet( timeStep, s=atm_cpl_dt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! For initial run, unpack the import state, run cam radiation/clouds and return
    ! For restart run, read the import state from the restart and run radiation/clouds and return
    !--------------------------------

    ! Note - cam_run1 is called only for the purposes of finishing the
    ! flux averaged calculation to compute cam-out
    ! Note - cam_run1 is called on restart only to have cam internal state consistent with the
    ! cam_out state sent to the coupler

    if (stepno == 0) then
       call import_fields( gcomp, cam_in, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call cam_run1 ( cam_in, cam_out )

       call export_fields( gcomp, cam_out, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call cam_read_srfrest( gcomp, clock, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fields( gcomp, cam_in, restart_init=.true., rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call cam_run1 ( cam_in, cam_out )

       call export_fields( gcomp, cam_out, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Compute time of next radiation computation, like in run method for exact restart
    dtime = get_step_size()
    nstep = get_nstep()
    if (nstep < 1 .or. dtime < atm_cpl_dt) then
       nextsw_cday = radiation_nextsw_cday()
    else if (dtime == atm_cpl_dt) then
       caldayp1 = get_curr_calday(offset=int(dtime))
       nextsw_cday = radiation_nextsw_cday()
       if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
    else
       call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
    end if

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug_flag > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !--------------------------------
    ! CAM data is now fully initialized
    !--------------------------------

    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n=1, fieldCount
       call ESMF_StateGet(exportState, itemName=fieldNameList(n), field=field, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_SetAttribute(field, name="Updated", value="true", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do
    deallocate(fieldNameList)

    ! check whether all Fields in the exportState are "Updated"
    if (NUOPC_IsUpdated(exportState)) then
       call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

       call ESMF_LogWrite("CAM - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! End redirection of share output to cam log
    !--------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_DataInitialize:end::',lbnum)
    endif
#endif

  end subroutine DataInitialize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    use time_manager, only : advance_timestep
    use time_manager, only : get_curr_calday, get_curr_date, get_nstep, get_step_size

    ! !DESCRIPTION:
    ! Run CAM

    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! !LOCAL VARIABLES:
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
    integer                 :: shrloglev   ! original log level
    character(CL)           :: case_name   ! case name
    real(r8)                :: eccen
    real(r8)                :: obliqr
    real(r8)                :: lambm0
    real(r8)                :: mvelpp
    logical                 :: dosend      ! true => send data back to driver
    integer                 :: dtime       ! time step increment (sec)
    integer                 :: atm_cpl_dt  ! driver atm coupling time step
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
    real(r8)                :: caldayp1    ! CAM calendar day for for next cam time step
    real(r8)                :: nextsw_cday ! calendar of next atm shortwave
    logical                 :: rstwr       ! .true. ==> write restart file before returning
    logical                 :: nlend       ! Flag signaling last time-step
    integer                 :: lbnum
    integer                 :: dbrc
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (iulog)

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !--------------------------------
    ! Determine current time
    !--------------------------------

    call ESMF_ClockGet( clock, currTime=currTime)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_GetScalar(importState, flds_scalar_index_nextsw_cday, nextsw_cday, &
         flds_scalar_name, flds_scalar_num, rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------
    ! Obtain orbital values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call t_startf ('CAM_import')
    call shr_nuopc_methods_State_diagnose(importState, string=subname//':IS', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call import_fields( gcomp, cam_in, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf  ('CAM_import')

    !--------------------------------
    ! Run cam
    !--------------------------------

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
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          rstwr = .true.
          call ESMF_AlarmRingerOff( alarm, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          rstwr = .false.
       endif

       ! Determine if time to stop

       call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

       ! Map output from cam to nuopc state fields

       call t_startf ('CAM_export')
       call export_fields( gcomp, cam_out, rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call t_stopf ('CAM_export')

    end do

    !--------------------------------
    ! Set the coupling scalars
    !--------------------------------

    ! Return time of next radiation calculation - albedos will need to be
    ! calculated by each surface model at this time

    call ESMF_ClockGet( clock, TimeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=atm_cpl_dt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    dtime = get_step_size()
    if (dtime < atm_cpl_dt) then
       nextsw_cday = radiation_nextsw_cday()
    else if (dtime == atm_cpl_dt) then
       caldayp1 = get_curr_calday(offset=int(dtime))
       nextsw_cday = radiation_nextsw_cday()
       if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
    else
       call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
    end if

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Write merged surface data restart file if appropriate
    !--------------------------------

    if (rstwr) then
       call cam_write_srfrest( gcomp, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Check for consistency of internal cam clock with master sync clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! CAM is actually 2 coupling intervals (or physics time steps) ahead of the driver clock
    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-2*dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) )then
       write(iulog,*)' cam ymd=',ymd     ,'  cam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( subname//': CAM clock is not in sync with master Sync Clock' )
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug_flag > 1) then
       call shr_nuopc_methods_State_diagnose(exportState, string=subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing CAM from: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','atm_comp_nuopc_ModelAdvance:end::',lbnum)
    endif
#endif

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    use shr_nuopc_time_mod, only : shr_nuopc_time_set_component_stop_alarm
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
    character(len=128)       :: name
    integer                  :: alarmcount
    integer                  :: dbrc
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO, rc=dbrc)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call shr_nuopc_time_set_component_stop_alarm(gcomp, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer :: shrlogunit            ! original log unit
    integer :: shrloglev             ! original log level
    integer :: dbrc
    character(*), parameter :: F00   = "('(atm_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(atm_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'CAM: end of main integration loop'
       write(iulog,F91)
    end if

    call shr_file_setLogUnit (shrlogunit)

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine cam_read_srfrest( gcomp, clock, rc )

    use cam_pio_utils     , only : cam_pio_openfile, cam_pio_closefile, pio_subsystem
    use cam_initfiles     , only : cam_initfiles_get_caseid, cam_initfiles_get_restdir
    use atm_import_export , only : state_getfldptr
    use filenames         , only : interpret_filename_spec
    use phys_grid         , only : ngcols
    use pio               , only : file_desc_t, io_desc_t, var_desc_t
    use pio               , only : file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim
    use pio               , only : pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp
    use pio               , only : pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid
    use pio               , only : pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling

    ! Arguments
    type(ESMF_GridComp)             :: gcomp
    type(ESMF_Clock), intent(inout) :: clock
    integer         , intent(out)   :: rc

    ! Local variables
    type(ESMF_State)                   :: importState, exportState
    integer                            :: rcode        ! return error code
    integer                            :: nf
    type(file_desc_t)                  :: file
    type(io_desc_t)                    :: iodesc
    integer                            :: fieldCount_import
    integer                            :: fieldCount_export
    character(ESMF_MAXSTR),allocatable :: fieldNameList_import(:)
    character(ESMF_MAXSTR),allocatable :: fieldNameList_export(:)
    type(var_desc_t)                   :: varid
    real(r8), pointer                  :: fldptr(:)
    logical                            :: exists
    type(ESMF_Time)                    :: currTime      ! time at previous interval
    integer                            :: yr_spec       ! Current year
    integer                            :: mon_spec      ! Current month
    integer                            :: day_spec      ! Current day
    integer                            :: sec_spec      ! Current time of day (sec)
    character(len=256)                 :: fname_srf_cam ! surface restart filename
    character(len=256)                 :: pname_srf_cam ! surface restart full pathname
    character(len=*),parameter  :: subname=trim(modName)//':(cam_read_srfrest)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine and surface restart dataset

    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yr_spec, mm=mon_spec, dd=day_spec, s=sec_spec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=cam_initfiles_get_caseid(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(cam_initfiles_get_restdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)

    ! Open restart file

    call cam_pio_openfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)

    ! Read in import state data

    call ESMF_StateGet(importState, itemCount=fieldCount_import, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldNameList_import(fieldCount_import))
    call ESMF_StateGet(importState, itemNameList=fieldNameList_import, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do nf=1,fieldCount_import
       call pio_seterrorhandling(File, pio_bcast_error)
       if (trim(fieldNameList_import(nf)) /= flds_scalar_name) then
          call state_getfldptr(importState, trim(fieldNameList_import(nf)), fldptr, exists, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       rcode = pio_inq_varid(File,'x2a_'//trim(fieldNameList_import(nf)) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, fldptr, rcode)
       else
	  if (masterproc) then
             write(iulog,*)'cam_read_srfrest warning: field ',trim(fieldNameList_import(nf)),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          fldptr(:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do
    deallocate(fieldNameList_import)

    ! Read in export state data

    call ESMF_StateGet(exportState, itemCount=fieldCount_export, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldNameList_export(fieldCount_export))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList_export, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do nf=1,fieldCount_export
       call pio_seterrorhandling(File, pio_bcast_error)
       if (trim(fieldNameList_export(nf)) /= flds_scalar_name) then
          call state_getfldptr(exportState, trim(fieldNameList_export(nf)), fldptr, exists, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       rcode = pio_inq_varid(File,'x2a_'//trim(fieldNameList_export(nf)) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, fldptr, rcode)
       else
	  if (masterproc) then
             write(iulog,*)'cam_read_srfrest warning: field ',trim(fieldNameList_export(nf)),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          fldptr(:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do
    deallocate(fieldNameList_export)

    ! Close file
    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(File)

  end subroutine cam_read_srfrest

  !===========================================================================================

  subroutine cam_write_srfrest( gcomp, yr_spec, mon_spec, day_spec, sec_spec, rc )

    use cam_pio_utils       , only : cam_pio_createfile, cam_pio_closefile, pio_subsystem
    use cam_history_support , only : fillvalue
    use filenames           , only : interpret_filename_spec
    use atm_import_export   , only : state_getfldptr
    use phys_grid           , only : ngcols
    use pio                 , only : file_desc_t, io_desc_t, var_desc_t
    use pio                 , only : file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim
    use pio                 , only : pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp
    use pio                 , only : pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid
    use pio                 , only : pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer , intent(in) :: yr_spec  ! Simulation year
    integer , intent(in) :: mon_spec ! Simulation month
    integer , intent(in) :: day_spec ! Simulation day
    integer , intent(in) :: sec_spec ! Seconds into current simulation day
    integer , intent(out):: rc       ! error code
    !
    ! Local variables
    type(ESMF_State)                   :: importState, exportState
    integer                            :: rcode        ! return error code
    integer                            :: dimid(1), nf, n
    type(file_desc_t)                  :: file
    type(io_desc_t)                    :: iodesc
    integer                            :: fieldCount_import
    integer                            :: fieldCount_export
    character(ESMF_MAXSTR),allocatable :: fieldNameList_import(:)
    character(ESMF_MAXSTR),allocatable :: fieldNameList_export(:)
    type(var_desc_t), pointer          :: varid_import(:)
    type(var_desc_t), pointer          :: varid_export(:)
    real(r8), pointer                  :: fldptr(:)
    logical                            :: exists
    character(len=256)                 :: fname_srf_cam       ! surface restart filename
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine and open surface restart dataset
    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )

    call cam_pio_createfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)

    ! Get import and export states
    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Define import and export state variable ids
    call ESMF_StateGet(importState, itemCount=fieldCount_import, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldNameList_import(fieldCount_import))
    call ESMF_StateGet(importState, itemNameList=fieldNameList_import, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(exportState, itemCount=fieldCount_export, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(fieldNameList_export(fieldCount_export))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList_export, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    allocate(varid_import(fieldCount_import))
    allocate(varid_export(fieldCount_export))

    rcode = pio_def_dim(File, 'x2a_nx', ngcols, dimid(1))
    do nf = 1,fieldCount_import
       if (trim(fieldNameList_import(nf)) /= flds_scalar_name) then
          rcode = pio_def_var(File,'x2a_'//trim(fieldNameList_import(nf)), PIO_DOUBLE, dimid, varid_import(nf))
          rcode = pio_put_att(File, varid_import(nf), "_fillvalue", fillvalue)
       end if
    end do
    rcode = pio_def_dim(File, 'a2x_nx', ngcols, dimid(1))
    do nf = 1,fieldCount_export
       if (trim(fieldNameList_export(nf)) /= flds_scalar_name) then
          rcode = pio_def_var(File,'a2x_'//trim(fieldNameList_export(nf)), PIO_DOUBLE, dimid, varid_export(nf))
          rcode = pio_put_att(File, varid_export(nf), "_fillvalue", fillvalue)
       end if
    end do

    ! End definition phase
    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    ! Write the restart data for the import and export states
    do nf = 1,fieldCount_import
       if (trim(fieldNameList_import(nf)) /= flds_scalar_name) then
          call state_getfldptr(importstate, trim(fieldNameList_import(nf)), fldptr, exists, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (exists) then
             call pio_write_darray(File, varid_import(nf), iodesc, fldptr, rcode)
          end if
       end if
    end do
    do nf = 1,fieldCount_export
       if (trim(fieldNameList_export(nf)) /= flds_scalar_name) then
          call state_getfldptr(exportstate, trim(fieldNameList_export(nf)), fldptr, exists, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (exists) then
             call pio_write_darray(File, varid_export(nf), iodesc, fldptr, rcode)
          end if
       end if
    end do

    deallocate(varid_import)
    deallocate(varid_export)

    ! close the file
    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(File)

  end subroutine cam_write_srfrest

end module atm_comp_nuopc
