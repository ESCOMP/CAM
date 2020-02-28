module cam_mpas_subdriver

!-------------------------------------------------------------------------------
!
! Handles the initialization of MPAS infrastructure in several phases
!
! This module mimics the functionality provided by stand-alone MPAS-Atmosphere's
! sub-driver (i.e., mpas_subdriver.F), but with minor modifications to code to
! accomodate the fact that MPAS-A is being driven by CAM.
!
!-------------------------------------------------------------------------------


    use mpas_derived_types, only : core_type, dm_info, domain_type, MPAS_Clock_type

    public :: cam_mpas_init_phase1, &
              cam_mpas_init_phase2, &
              cam_mpas_init_phase3, &
              cam_mpas_init_phase4, &
              cam_mpas_get_global_dims, &
              cam_mpas_get_global_coords, &
              cam_mpas_get_global_blocks, &
              cam_mpas_read_static, &
              cam_mpas_compute_unit_vectors, &
              cam_mpas_update_halo, &
              cam_mpas_cell_to_edge_winds, &
              cam_mpas_run, &
              cam_mpas_finalize
    public :: corelist, domain_ptr

    private


    type (core_type), pointer :: corelist => null()
    type (domain_type), pointer :: domain_ptr => null()
    type (MPAS_Clock_type), pointer :: clock => null()

    !
    ! This interface should be compatible with CAM's endrun routine
    !
    abstract interface
       subroutine halt_model(mesg, ierr)
           use shr_kind_mod, only : shr_kind_in
           character(len=*), intent(in), optional :: mesg
           integer(kind=shr_kind_in), intent(in), optional :: ierr
       end subroutine halt_model
    end interface


contains


    !-----------------------------------------------------------------------
    !  routine cam_mpas_init_phase1
    !
    !> \brief Tracks mpas_init up to the point of reading namelist
    !> \author Michael Duda
    !> \date   19 April 2019
    !> \details
    !>  This routine follows the stand-alone MPAS subdriver up to, but not
    !>  including, the point where namelists are read.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_init_phase1(mpicom, endrun, logUnits)

       use mpas_domain_routines, only : mpas_allocate_domain
       use mpas_framework, only : mpas_framework_init_phase1
       use atm_core_interface, only : atm_setup_core, atm_setup_domain
       use mpas_pool_routines, only : mpas_pool_add_config

       implicit none

       ! Dummy argument
       integer, intent(in) :: mpicom
       procedure(halt_model) :: endrun
       integer, dimension(2), intent(in) :: logUnits

       ! Local variables
       integer :: ierr

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase1'


       allocate(corelist)
       nullify(corelist % next)

       allocate(corelist % domainlist)
       nullify(corelist % domainlist % next)

       domain_ptr => corelist % domainlist
       domain_ptr % core => corelist

       call mpas_allocate_domain(domain_ptr)


       !
       ! Initialize MPAS infrastructure (principally, the mpas_dmpar module)
       !
       call mpas_framework_init_phase1(domain_ptr % dminfo, mpi_comm=mpicom)

       call atm_setup_core(corelist)
       call atm_setup_domain(domain_ptr)


       ! Set up the log manager as early as possible so we can use it for any errors/messages during subsequent init steps
       ! We need:
       ! 1) domain_ptr to be allocated,
       ! 2) dmpar_init complete to access dminfo,
       ! 3) *_setup_core to assign the setup_log function pointer
       ierr = domain_ptr % core % setup_log(domain_ptr % logInfo, domain_ptr, unitNumbers=logUnits)
       if ( ierr /= 0 ) then
          call endrun(subname//': FATAL: Log setup failed for MPAS-A dycore')
       end if

    end subroutine cam_mpas_init_phase1


    !-----------------------------------------------------------------------
    !  routine cam_mpas_init_phase2
    !
    !> \brief Tracks mpas_init after namelists have been read
    !> \author Michael Duda
    !> \date   19 April 2019
    !> \details
    !>  This routine follows the stand-alone MPAS subdriver from the point
    !>  where we call the second phase of MPAS framework initialization up
    !>  to the check on the existence of the streams.<core> file.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_init_phase2(pio_subsystem, endrun, cam_calendar)

       use mpas_log, only : mpas_log_write
       use mpas_kind_types, only : ShortStrKIND
       use mpas_derived_types, only : MPAS_LOG_ERR
       use pio_types, only : iosystem_desc_t

       use mpas_framework, only : mpas_framework_init_phase2

       implicit none

       type (iosystem_desc_t), pointer :: pio_subsystem
       procedure(halt_model) :: endrun
       character(len=*), intent(in) :: cam_calendar

       integer :: ierr
       logical :: streamsExists

       character(len=ShortStrKIND) :: mpas_calendar

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase2'


       !
       ! Translate between CAM calendar names and MPAS calendar names
       !
       select case(trim(cam_calendar))
          case ('noleap')
             mpas_calendar = 'gregorian_noleap'
          case ('gregorian')
             mpas_calendar = 'gregorian'
          case default
             call endrun(subname//': FATAL: Unrecognized calendar type '''//trim(cam_calendar)//'''')
       end select

       ! 4) Continue with normal procedure from MPAS subdriver
       call mpas_framework_init_phase2(domain_ptr, io_system=pio_subsystem, calendar=trim(mpas_calendar))

       ierr = domain_ptr % core % define_packages(domain_ptr % packages)
       if ( ierr /= 0 ) then
          call endrun(subname//': FATAL: Package definition failed for core '//trim(domain_ptr % core % coreName))
       end if

       ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % iocontext)
       if ( ierr /= 0 ) then
          call endrun(subname//': FATAL: Package setup failed for core '//trim(domain_ptr % core % coreName))
       end if

       ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
       if ( ierr /= 0 ) then
          call endrun(subname//': FATAL: Decomposition setup failed for core '//trim(domain_ptr % core % coreName))
       end if

       ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
       if ( ierr /= 0 ) then
          call endrun(subname//': FATAL: Clock setup failed for core '//trim(domain_ptr % core % coreName))
       end if

#ifdef MPAS_USE_STREAMS
       call mpas_log_write('Reading streams configuration from file '//trim(domain_ptr % streams_filename))
       inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

       if ( .not. streamsExists ) then
          call endrun(subname//': FATAL: Streams file '//trim(domain_ptr % streams_filename)//' does not exist.')
       end if
#endif

       ! At this point, we should be ready to set up decompositions, build halos, allocate blocks, etc. in dyn_grid_init

    end subroutine cam_mpas_init_phase2


    !-----------------------------------------------------------------------
    !  routine cam_mpas_init_phase3
    !
    !> \brief Finish MPAS initialization
    !> \author Michael Duda
    !> \date   19 April 2019
    !> \details
    !>  This routine completes the initialization of the MPAS infrastructure and
    !>  the MPAS core.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_init_phase3(fh_ini, endrun)

       use mpas_log, only : mpas_log_write
       use mpas_derived_types, only : MPAS_LOG_ERR
       use pio, only : file_desc_t
       use iso_c_binding, only : c_int, c_char, c_ptr, c_loc

       use mpas_derived_types, only : MPAS_Time_type, MPAS_TimeInterval_type
       use mpas_derived_types, only : MPAS_IO_PNETCDF, MPAS_IO_PNETCDF5, MPAS_IO_NETCDF, MPAS_IO_NETCDF4
       use mpas_derived_types, only : MPAS_START_TIME
       use mpas_derived_types, only : MPAS_STREAM_MGR_NOERR
       use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time, mpas_expand_string, mpas_set_time, &
                                    mpas_set_timeInterval
       use mpas_stream_manager, only : MPAS_stream_mgr_init, mpas_build_stream_filename, MPAS_stream_mgr_validate_streams
       use mpas_kind_types, only : StrKIND
       use mpas_c_interfacing, only : mpas_c_to_f_string, mpas_f_to_c_string
       use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2

       implicit none

       type (file_desc_t), pointer :: fh_ini
       procedure(halt_model) :: endrun

       integer :: ierr
       character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
       integer(kind=c_int) :: c_comm
       integer(kind=c_int) :: c_ierr
       type (c_ptr) :: mgr_p
       character(len=StrKIND) :: mesh_stream
       character(len=StrKIND) :: mesh_filename
       character(len=StrKIND) :: mesh_filename_temp
       character(len=StrKIND) :: ref_time_temp
       character(len=StrKIND) :: filename_interval_temp
       character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
       character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
       character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
       character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
       character(kind=c_char), dimension(StrKIND+1) :: c_iotype
       type (MPAS_Time_type) :: start_time
       type (MPAS_Time_type) :: ref_time
       type (MPAS_TimeInterval_type) :: filename_interval
       character(len=StrKIND) :: start_timestamp
       character(len=StrKIND) :: iotype
       logical :: streamsExists
       integer :: mesh_iotype
       integer :: blockID
       character(len=StrKIND) :: timeStamp

       interface
          subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
             use iso_c_binding, only : c_char, c_ptr, c_int
             character(kind=c_char), dimension(*), intent(in) :: xmlname
             type (c_ptr), intent(inout) :: mgr_p
             integer(kind=c_int), intent(inout) :: comm
             integer(kind=c_int), intent(out) :: ierr
          end subroutine xml_stream_parser

          subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
             use iso_c_binding, only : c_char, c_int
             character(kind=c_char), dimension(*), intent(in) :: xmlname
             character(kind=c_char), dimension(*), intent(in) :: streamname
             integer(kind=c_int), intent(inout) :: comm
             character(kind=c_char), dimension(*), intent(out) :: filename
             character(kind=c_char), dimension(*), intent(out) :: ref_time
             character(kind=c_char), dimension(*), intent(out) :: filename_interval
             character(kind=c_char), dimension(*), intent(out) :: io_type
             integer(kind=c_int), intent(out) :: ierr
          end subroutine xml_stream_get_attributes
       end interface

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase3'

   
#ifdef MPAS_USE_STREAMS
       !
       ! Using information from the namelist, a graph.info file, and a file containing
       !    mesh fields, build halos and allocate blocks in the domain
       !
       ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
       if ( ierr /= 0 ) then
          call endrun('Failed to find mesh stream for core '//trim(domain_ptr % core % coreName))
       end if

       call mpas_f_to_c_string(domain_ptr % streams_filename, c_filename)
       call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
       c_comm = domain_ptr % dminfo % comm
       call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                      c_mesh_filename_temp, c_ref_time_temp, &
                                      c_filename_interval_temp, c_iotype, c_ierr)
       if (c_ierr /= 0) then
          call endrun('stream xml get attribute failed: '//trim(domain_ptr % streams_filename))
       end if
       call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
       call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
       call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
       call mpas_c_to_f_string(c_iotype, iotype)

       if (trim(iotype) == 'pnetcdf') then
          mesh_iotype = MPAS_IO_PNETCDF
       else if (trim(iotype) == 'pnetcdf,cdf5') then
          mesh_iotype = MPAS_IO_PNETCDF5
       else if (trim(iotype) == 'netcdf') then
          mesh_iotype = MPAS_IO_NETCDF
       else if (trim(iotype) == 'netcdf4') then
          mesh_iotype = MPAS_IO_NETCDF4
       else
          mesh_iotype = MPAS_IO_PNETCDF
       end if

       start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
       if ( trim(ref_time_temp) == 'initial_time' ) then
           call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
       end if

       blockID = -1
       if ( trim(filename_interval_temp) == 'none' ) then
           call mpas_expand_string(ref_time_temp, blockID, mesh_filename_temp, mesh_filename)
       else
           call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
           call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
           call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, blockID, mesh_filename, ierr)
       end if

       call mpas_log_write(' ** Attempting to bootstrap MPAS framework using stream: ' // trim(mesh_stream))
#else
       mesh_filename = 'external mesh file'
#endif

! Use the call below if we intend to bootstrap from an input stream; also define MPAS_USE_STREAMS
!       call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

! Use the call below if we intend to supply the bootstrapping process with an external PIO file_desc_t, fh_ini
       mesh_iotype = MPAS_IO_NETCDF  ! Not actually used
       call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype, pio_file_desc=fh_ini)

#ifdef MPAS_USE_STREAMS
       !
       ! Set up run-time streams
       !
       call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % ioContext, domain_ptr % clock, &
                                 domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

       call add_stream_attributes(domain_ptr)

       ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
       if ( ierr /= 0 ) then
          call endrun('Immutable streams setup failed for core '//trim(domain_ptr % core % coreName))
       end if

       mgr_p = c_loc(domain_ptr % streamManager)
       call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
       if (c_ierr /= 0) then
          call endrun('xml stream parser failed: '//trim(domain_ptr % streams_filename))
       end if

       !
       ! Validate streams after set-up
       !
       call mpas_log_write(' ** Validating streams')
       call MPAS_stream_mgr_validate_streams(domain_ptr % streamManager, ierr = ierr)
       if ( ierr /= MPAS_STREAM_MGR_NOERR ) then
          call endrun('ERROR: Validation of streams failed for core ' // trim(domain_ptr % core % coreName))
       end if
#endif

       !
       ! Finalize the setup of blocks and fields
       !
! Use the call below if we intend to bootstrap from an input stream; also define MPAS_USE_STREAMS
!       call mpas_bootstrap_framework_phase2(domain_ptr)

! Use the call below if we intend to supply the bootstrapping process with an external PIO file_desc_t, fh_ini
       call mpas_bootstrap_framework_phase2(domain_ptr, pio_file_desc=fh_ini)

    end subroutine cam_mpas_init_phase3


    !-----------------------------------------------------------------------
    !  routine cam_mpas_init_phase4
    !
    !> \brief Finish MPAS initialization
    !> \author Michael Duda
    !> \date   29 February 2020
    !> \details
    !>  This routine completes the initialization of the MPAS core, essentially
    !>  following what is done in mpas_atm_core.F::atm_core_init(), but without
    !>  any calls to the MPAS-A diagnostics framework or the MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_init_phase4(endrun)

       use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time, MPAS_START_TIME
       use mpas_kind_types, only : StrKIND, RKIND
       use mpas_atm_dimensions, only : mpas_atm_set_dims
       use mpas_atm_threading, only : mpas_atm_threading_init
       use mpas_derived_types, only : mpas_pool_type, field2DReal, MPAS_Time_type
       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_config, &
                                      mpas_pool_get_field, mpas_pool_get_array, mpas_pool_initialize_time_levels
       use atm_core, only : atm_mpas_init_block, core_clock => clock
       use mpas_dmpar, only : mpas_dmpar_exch_halo_field

       implicit none

       procedure(halt_model) :: endrun

       real (kind=RKIND), pointer :: dt

       character(len=StrKIND) :: timeStamp
       integer :: i
       logical, pointer :: config_do_restart

       type (mpas_pool_type), pointer :: state
       type (mpas_pool_type), pointer :: mesh
       type (mpas_pool_type), pointer :: diag
       type (field2DReal), pointer :: u_field, pv_edge_field, ru_field, rw_field
       character (len=StrKIND), pointer :: xtime
       character (len=StrKIND), pointer :: initial_time1, initial_time2
       type (MPAS_Time_Type) :: startTime

       integer, pointer :: nVertLevels, maxEdges, maxEdges2, num_scalars

       integer :: ierr
       character(len=StrKIND) :: startTimeStamp


       !
       ! Setup threading
       !
       call mpas_atm_threading_init(domain_ptr % blocklist, ierr)
       if ( ierr /= 0 ) then
          call endrun('Threading setup failed for core '//trim(domain_ptr % core % coreName))
       end if


       !
       ! Set up inner dimensions used by arrays in optimized dynamics routines
       !
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_get_dimension(state, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(state, 'maxEdges', maxEdges)
       call mpas_pool_get_dimension(state, 'maxEdges2', maxEdges2)
       call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
       call mpas_atm_set_dims(nVertLevels, maxEdges, maxEdges2, num_scalars)

       !
       ! Set "local" clock to point to the clock contained in the domain type
       !
       clock => domain_ptr % clock
       core_clock => domain_ptr % clock


       call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_do_restart', config_do_restart)
       call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)


       if (.not. config_do_restart) then
           call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
           call mpas_pool_initialize_time_levels(state)
       end if


       !
       ! Set startTimeStamp based on the start time of the simulation clock
       !
       startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
       call mpas_get_time(startTime, dateTimeString=startTimeStamp) 


       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_get_field(state, 'u', u_field, 1)
       call mpas_dmpar_exch_halo_field(u_field)

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)

       call atm_mpas_init_block(domain_ptr % dminfo, domain_ptr % streamManager, domain_ptr % blocklist, mesh, dt)

       call mpas_pool_get_array(state, 'xtime', xtime, 1)
       xtime = startTimeStamp

       ! Initialize initial_time in second time level. We need to do this because initial state
       ! is read into time level 1, and if we write output from the set of state arrays that
       ! represent the original time level 2, the initial_time field will be invalid.
       call mpas_pool_get_array(state, 'initial_time', initial_time1, 1)
       call mpas_pool_get_array(state, 'initial_time', initial_time2, 2)
       initial_time2 = initial_time1

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)
       call mpas_pool_get_field(diag, 'pv_edge', pv_edge_field)
       call mpas_dmpar_exch_halo_field(pv_edge_field)

       call mpas_pool_get_field(diag, 'ru', ru_field)
       call mpas_dmpar_exch_halo_field(ru_field)

       call mpas_pool_get_field(diag, 'rw', rw_field)
       call mpas_dmpar_exch_halo_field(rw_field)

    end subroutine cam_mpas_init_phase4


    !-----------------------------------------------------------------------
    !  routine cam_mpas_get_global_dims
    !
    !> \brief  Returns global mesh dimensions
    !> \author Michael Duda
    !> \date   22 August 2019
    !> \details
    !>  This routine returns on all tasks the number of global cells, edges,
    !>  vertices, maxEdges, vertical layers, and the maximum number of cells owned by any task.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_get_global_dims(nCellsGlobal, nEdgesGlobal, nVerticesGlobal, maxEdges, nVertLevels, maxNCells)

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension
       use mpas_derived_types, only : mpas_pool_type
       use mpas_dmpar, only : mpas_dmpar_sum_int, mpas_dmpar_max_int

       implicit none

       integer, intent(out) :: nCellsGlobal
       integer, intent(out) :: nEdgesGlobal
       integer, intent(out) :: nVerticesGlobal
       integer, intent(out) :: maxEdges
       integer, intent(out) :: nVertLevels
       integer, intent(out) :: maxNCells

       integer, pointer :: nCellsSolve
       integer, pointer :: nEdgesSolve
       integer, pointer :: nVerticesSolve
       integer, pointer :: maxEdgesLocal
       integer, pointer :: nVertLevelsLocal

       type (mpas_pool_type), pointer :: meshPool


       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
       call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdgesLocal)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsLocal)

       call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)
       call mpas_dmpar_sum_int(domain_ptr % dminfo, nEdgesSolve, nEdgesGlobal)
       call mpas_dmpar_sum_int(domain_ptr % dminfo, nVerticesSolve, nVerticesGlobal)

       maxEdges = maxEdgesLocal
       nVertLevels = nVertLevelsLocal

       call mpas_dmpar_max_int(domain_ptr % dminfo, nCellsSolve, maxNCells)

    end subroutine cam_mpas_get_global_dims


    !-----------------------------------------------------------------------
    !  routine cam_mpas_get_global_coords
    !
    !> \brief  Returns global coordinate arrays
    !> \author Michael Duda
    !> \date   22 August 2019
    !> \details
    !>  This routine returns on all tasks arrays of latitude, longitude, and cell
    !>  area for all (global) cells.
    !>
    !>  It is assumed that latCellGlobal, lonCellGlobal, and areaCellGlobal have
    !>  been allocated by the caller with a size equal to the global number of
    !>  cells in the mesh.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_get_global_coords(latCellGlobal, lonCellGlobal, areaCellGlobal)

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_derived_types, only : mpas_pool_type
       use mpas_kind_types, only : RKIND
       use mpas_dmpar, only : mpas_dmpar_sum_int, mpas_dmpar_max_int, mpas_dmpar_max_real_array

       implicit none

       real (kind=RKIND), dimension(:), intent(out) :: latCellGlobal
       real (kind=RKIND), dimension(:), intent(out) :: lonCellGlobal
       real (kind=RKIND), dimension(:), intent(out) :: areaCellGlobal

       integer :: iCell

       integer, pointer :: nCellsSolve
       integer, dimension(:), pointer :: indexToCellID

       type (mpas_pool_type), pointer :: meshPool
       integer :: nCellsGlobal

       real (kind=RKIND), dimension(:), pointer :: latCell
       real (kind=RKIND), dimension(:), pointer :: lonCell
       real (kind=RKIND), dimension(:), pointer :: areaCell
       real (kind=RKIND), dimension(:), pointer :: temp


       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

       call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

       ! check: size(latCellGlobal) ?= nCellsGlobal

       allocate(temp(nCellsGlobal))
       
       !
       ! latCellGlobal
       !
       temp(:) = -huge(temp(0))
       do iCell=1,nCellsSolve
           temp(indexToCellID(iCell)) = latCell(iCell)
       end do

       call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, latCellGlobal)

       !
       ! lonCellGlobal
       !
       temp(:) = -huge(temp(0))
       do iCell=1,nCellsSolve
           temp(indexToCellID(iCell)) = lonCell(iCell)
       end do

       call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, lonCellGlobal)

       !
       ! areaCellGlobal
       !
       temp(:) = -huge(temp(0))
       do iCell=1,nCellsSolve
           temp(indexToCellID(iCell)) = areaCell(iCell)
       end do

       call mpas_dmpar_max_real_array(domain_ptr % dminfo, nCellsGlobal, temp, areaCellGlobal)

       deallocate(temp)

    end subroutine cam_mpas_get_global_coords


    !-----------------------------------------------------------------------
    !  routine cam_mpas_get_global_blocks
    !
    !> \brief  Returns global block indexing arrays
    !> \author Michael Duda
    !> \date   22 August 2019
    !> \details
    !>  Returns arrays with information about the number of columns in each global block,
    !>  which column indices are in each block, and which global block contains each
    !>  global column.
    !>
    !>  It is assumed that nCellsPerBlock, indexToCellIDBlock, indexToBlockID, and
    !>  localCellIDBlock have been allocated by the caller with dimensions:
    !>
    !>  nCellsPerBlock(num_blocks_global)
    !>  indexToCellIDBlock(maxNCells, num_blocks_global)
    !>  indexToBlockID(nCellsGlobal)
    !>  localCellIDBlock(nCellsGlobal)
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_get_global_blocks(nCellsPerBlock, indexToCellIDBlock, indexToBlockID, localCellIDBlock)

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_derived_types, only : mpas_pool_type
       use mpas_dmpar, only : mpas_dmpar_max_int_array

       implicit none

       integer, dimension(:), intent(out) :: nCellsPerBlock
       integer, dimension(:,:), intent(out) :: indexToCellIDBlock
       integer, dimension(:), intent(out) :: indexToBlockID
       integer, dimension(:), intent(out) :: localCellIDBlock

       integer :: iCell
       integer :: owningBlock, localCellID
       type (mpas_pool_type), pointer :: meshPool
       integer, pointer :: nCellsSolve
       integer, dimension(:), pointer :: indexToCellID
       integer, dimension(:), pointer :: temp1d


       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

       !
       ! nCellsPerBlock
       !
       allocate(temp1d(size(nCellsPerBlock)))
       temp1d(:) = 0
       temp1d(domain_ptr % dminfo % my_proc_id + 1) = nCellsSolve

       call mpas_dmpar_max_int_array(domain_ptr % dminfo, size(temp1d), temp1d, nCellsPerBlock)

       deallocate(temp1d)

       !
       ! indexToBlockID
       !
       allocate(temp1d(size(indexToBlockID)))
       temp1d(:) = -1
       do iCell=1,nCellsSolve
          temp1d(indexToCellID(iCell)) = domain_ptr % dminfo % my_proc_id + 1   ! 1-based block indices?
       end do

       call mpas_dmpar_max_int_array(domain_ptr % dminfo, size(temp1d), temp1d, indexToBlockID)

       deallocate(temp1d)

       !
       ! localCellIDBlock
       !
       allocate(temp1d(size(localCellIDBlock)))
       temp1d(:) = 0
       do iCell = 1, nCellsSolve
          temp1d(indexToCellID(iCell)) = iCell
       end do

       call mpas_dmpar_max_int_array(domain_ptr % dminfo, size(temp1d), temp1d, localCellIDBlock)

       deallocate(temp1d)

       !
       ! indexToCellIDBlock
       !
       indexToCellIDBlock(:,:) = 0
       do iCell = 1, size(localCellIDBlock) ! nCellsGlobal
          owningBlock = indexToBlockID(iCell) + 0            ! 1-based block indices?
          localCellID = localCellIDBlock(iCell)
          indexToCellIDBlock(localCellID, owningBlock) = iCell
       end do

    end subroutine cam_mpas_get_global_blocks


    !-----------------------------------------------------------------------
    !  routine cam_mpas_read_static
    !
    !> \brief  Reads time-invariant ("static") fields from an MPAS-A mesh file
    !> \author Michael Duda
    !> \date   6 January 2020
    !> \details
    !>  This routine takes as input an opened PIO file descriptor and a routine
    !>  to call if catastrophic errors are encountered. An MPAS stream is constructed
    !>  from this file descriptor, and most of the fields that exist in MPAS's
    !>  "mesh" pool are read from this stream.
    !>  Upon successful completion, valid mesh fields may be accessed from the mesh
    !>  pool.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_read_static(fh_ini, endrun)

       use pio, only : file_desc_t

       use mpas_io_streams, only : MPAS_createStream, MPAS_closeStream, MPAS_streamAddField, MPAS_readStream
       use mpas_derived_types, only : MPAS_IO_READ, MPAS_IO_NETCDF, MPAS_Stream_type, MPAS_pool_type, &
                                      field0DReal, field1DReal, field2DReal, field3DReal, field1DInteger, field2DInteger
       use mpas_pool_routines, only : MPAS_pool_get_subpool, MPAS_pool_get_field, MPAS_pool_create_pool, MPAS_pool_destroy_pool, &
                                      MPAS_pool_add_config
       use mpas_dmpar, only : MPAS_dmpar_exch_halo_field
       use mpas_stream_manager, only : postread_reindex

       implicit none

       type (file_desc_t), pointer :: fh_ini
       procedure(halt_model) :: endrun

       integer :: ierr
       type (MPAS_pool_type), pointer :: meshPool
       type (MPAS_pool_type), pointer :: reindexPool
       type (field1DReal), pointer :: latCell, lonCell, xCell, yCell, zCell
       type (field1DReal), pointer :: latEdge, lonEdge, xEdge, yEdge, zEdge
       type (field1DReal), pointer :: latVertex, lonVertex, xVertex, yVertex, zVertex
       type (field1DInteger), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
       type (field1DReal), pointer :: areaCell, areaTriangle, dcEdge, dvEdge, angleEdge
       type (field2DReal), pointer :: kiteAreasOnVertex, weightsOnEdge
       type (field1DReal), pointer :: meshDensity
       type (field1DInteger), pointer :: nEdgesOnCell, nEdgesOnEdge
       type (field2DInteger), pointer :: cellsOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, verticesOnCell, &
                                         verticesOnEdge, edgesOnVertex, cellsOnVertex
       type (field0DReal), pointer :: cf1, cf2, cf3
       type (field1DReal), pointer :: rdzw, dzu, rdzu, fzm, fzp
       type (field2DReal), pointer :: zgrid, zxu, zz
       type (field3DReal), pointer :: zb, zb3, deriv_two, cellTangentPlane, coeffs_reconstruct

       type (field2DReal), pointer :: edgeNormalVectors, localVerticalUnitVectors, defc_a, defc_b

       type (MPAS_Stream_type) :: mesh_stream


       call MPAS_createStream(mesh_stream, domain_ptr % ioContext, 'not_used', MPAS_IO_NETCDF, MPAS_IO_READ, &
                              pio_file_desc=fh_ini, ierr=ierr)

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)

       call mpas_pool_get_field(meshPool, 'latCell', latCell)
       call mpas_pool_get_field(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_field(meshPool, 'xCell', xCell)
       call mpas_pool_get_field(meshPool, 'yCell', yCell)
       call mpas_pool_get_field(meshPool, 'zCell', zCell)

       call mpas_pool_get_field(meshPool, 'latEdge', latEdge)
       call mpas_pool_get_field(meshPool, 'lonEdge', lonEdge)
       call mpas_pool_get_field(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_field(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_field(meshPool, 'zEdge', zEdge)

       call mpas_pool_get_field(meshPool, 'latVertex', latVertex)
       call mpas_pool_get_field(meshPool, 'lonVertex', lonVertex)
       call mpas_pool_get_field(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_field(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_field(meshPool, 'zVertex', zVertex)

       call mpas_pool_get_field(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_field(meshPool, 'indexToEdgeID', indexToEdgeID)
       call mpas_pool_get_field(meshPool, 'indexToVertexID', indexToVertexID)

       call mpas_pool_get_field(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_field(meshPool, 'areaTriangle', areaTriangle)
       call mpas_pool_get_field(meshPool, 'dcEdge', dcEdge)
       call mpas_pool_get_field(meshPool, 'dvEdge', dvEdge)
       call mpas_pool_get_field(meshPool, 'angleEdge', angleEdge)
       call mpas_pool_get_field(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
       call mpas_pool_get_field(meshPool, 'weightsOnEdge', weightsOnEdge)

       call mpas_pool_get_field(meshPool, 'meshDensity', meshDensity)

       call mpas_pool_get_field(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_field(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)

       call mpas_pool_get_field(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_field(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_field(meshPool, 'edgesOnEdge', edgesOnEdge)
       call mpas_pool_get_field(meshPool, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_field(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_field(meshPool, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_field(meshPool, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_field(meshPool, 'cellsOnVertex', cellsOnVertex)

       call mpas_pool_get_field(meshPool, 'cf1', cf1)
       call mpas_pool_get_field(meshPool, 'cf2', cf2)
       call mpas_pool_get_field(meshPool, 'cf3', cf3)

       call mpas_pool_get_field(meshPool, 'rdzw', rdzw)
       call mpas_pool_get_field(meshPool, 'dzu', dzu)
       call mpas_pool_get_field(meshPool, 'rdzu', rdzu)
       call mpas_pool_get_field(meshPool, 'fzm', fzm)
       call mpas_pool_get_field(meshPool, 'fzp', fzp)

       call mpas_pool_get_field(meshPool, 'zgrid', zgrid)
       call mpas_pool_get_field(meshPool, 'zxu', zxu)
       call mpas_pool_get_field(meshPool, 'zz', zz)
       call mpas_pool_get_field(meshPool, 'zb', zb)
       call mpas_pool_get_field(meshPool, 'zb3', zb3)

       call mpas_pool_get_field(meshPool, 'deriv_two', deriv_two)
       call mpas_pool_get_field(meshPool, 'cellTangentPlane', cellTangentPlane)
       call mpas_pool_get_field(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

       call mpas_pool_get_field(meshPool, 'edgeNormalVectors', edgeNormalVectors)
       call mpas_pool_get_field(meshPool, 'localVerticalUnitVectors', localVerticalUnitVectors)
       call mpas_pool_get_field(meshPool, 'defc_a', defc_a)
       call mpas_pool_get_field(meshPool, 'defc_b', defc_b)

       call MPAS_streamAddField(mesh_stream, latCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, lonCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, xCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, yCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zCell, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, latEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, lonEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, xEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, yEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zEdge, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, latVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, lonVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, xVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, yVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zVertex, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, indexToCellID, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, indexToEdgeID, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, indexToVertexID, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, areaCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, areaTriangle, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, dcEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, dvEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, angleEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, kiteAreasOnVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, weightsOnEdge, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, meshDensity, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, nEdgesOnCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, nEdgesOnEdge, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, cellsOnEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, edgesOnCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, edgesOnEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, cellsOnCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, verticesOnCell, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, verticesOnEdge, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, edgesOnVertex, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, cellsOnVertex, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, cf1, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, cf2, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, cf3, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, rdzw, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, dzu, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, rdzu, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, fzm, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, fzp, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, zgrid, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zxu, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zz, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zb, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, zb3, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, deriv_two, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, cellTangentPlane, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, coeffs_reconstruct, ierr=ierr)

       call MPAS_streamAddField(mesh_stream, edgeNormalVectors, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, localVerticalUnitVectors, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, defc_a, ierr=ierr)
       call MPAS_streamAddField(mesh_stream, defc_b, ierr=ierr)

       call MPAS_readStream(mesh_stream, 0, ierr=ierr)

       call MPAS_closeStream(mesh_stream, ierr=ierr)

       !
       ! Perform halo updates for all decomposed fields (i.e., fields with
       ! an outermost dimension of nCells, nVertices, or nEdges)
       !
       call MPAS_dmpar_exch_halo_field(latCell)
       call MPAS_dmpar_exch_halo_field(lonCell)
       call MPAS_dmpar_exch_halo_field(xCell)
       call MPAS_dmpar_exch_halo_field(yCell)
       call MPAS_dmpar_exch_halo_field(zCell)

       call MPAS_dmpar_exch_halo_field(latEdge)
       call MPAS_dmpar_exch_halo_field(lonEdge)
       call MPAS_dmpar_exch_halo_field(xEdge)
       call MPAS_dmpar_exch_halo_field(yEdge)
       call MPAS_dmpar_exch_halo_field(zEdge)

       call MPAS_dmpar_exch_halo_field(latVertex)
       call MPAS_dmpar_exch_halo_field(lonVertex)
       call MPAS_dmpar_exch_halo_field(xVertex)
       call MPAS_dmpar_exch_halo_field(yVertex)
       call MPAS_dmpar_exch_halo_field(zVertex)

       call MPAS_dmpar_exch_halo_field(indexToCellID)
       call MPAS_dmpar_exch_halo_field(indexToEdgeID)
       call MPAS_dmpar_exch_halo_field(indexToVertexID)

       call MPAS_dmpar_exch_halo_field(areaCell)
       call MPAS_dmpar_exch_halo_field(areaTriangle)
       call MPAS_dmpar_exch_halo_field(dcEdge)
       call MPAS_dmpar_exch_halo_field(dvEdge)
       call MPAS_dmpar_exch_halo_field(angleEdge)
       call MPAS_dmpar_exch_halo_field(kiteAreasOnVertex)
       call MPAS_dmpar_exch_halo_field(weightsOnEdge)

       call MPAS_dmpar_exch_halo_field(meshDensity)

       call MPAS_dmpar_exch_halo_field(nEdgesOnCell)
       call MPAS_dmpar_exch_halo_field(nEdgesOnEdge)

       call MPAS_dmpar_exch_halo_field(cellsOnEdge)
       call MPAS_dmpar_exch_halo_field(edgesOnCell)
       call MPAS_dmpar_exch_halo_field(edgesOnEdge)
       call MPAS_dmpar_exch_halo_field(cellsOnCell)
       call MPAS_dmpar_exch_halo_field(verticesOnCell)
       call MPAS_dmpar_exch_halo_field(verticesOnEdge)
       call MPAS_dmpar_exch_halo_field(edgesOnVertex)
       call MPAS_dmpar_exch_halo_field(cellsOnVertex)

       call MPAS_dmpar_exch_halo_field(zgrid)
       call MPAS_dmpar_exch_halo_field(zxu)
       call MPAS_dmpar_exch_halo_field(zz)
       call MPAS_dmpar_exch_halo_field(zb)
       call MPAS_dmpar_exch_halo_field(zb3)

       call MPAS_dmpar_exch_halo_field(deriv_two)
       call MPAS_dmpar_exch_halo_field(cellTangentPlane)
       call MPAS_dmpar_exch_halo_field(coeffs_reconstruct)

       call MPAS_dmpar_exch_halo_field(edgeNormalVectors)
       call MPAS_dmpar_exch_halo_field(localVerticalUnitVectors)
       call MPAS_dmpar_exch_halo_field(defc_a)
       call MPAS_dmpar_exch_halo_field(defc_b)

       !
       ! Re-index from global index space to local index space
       !
       call MPAS_pool_create_pool(reindexPool)

       call MPAS_pool_add_config(reindexPool, 'cellsOnEdge', 1)
       call MPAS_pool_add_config(reindexPool, 'edgesOnCell', 1)
       call MPAS_pool_add_config(reindexPool, 'edgesOnEdge', 1)
       call MPAS_pool_add_config(reindexPool, 'cellsOnCell', 1)
       call MPAS_pool_add_config(reindexPool, 'verticesOnCell', 1)
       call MPAS_pool_add_config(reindexPool, 'verticesOnEdge', 1)
       call MPAS_pool_add_config(reindexPool, 'edgesOnVertex', 1)
       call MPAS_pool_add_config(reindexPool, 'cellsOnVertex', 1)

       call postread_reindex(meshPool, reindexPool)

       call MPAS_pool_destroy_pool(reindexPool)

    end subroutine cam_mpas_read_static


    !-----------------------------------------------------------------------
    !  routine cam_mpas_compute_unit_vectors
    !
    !> \brief  Computes local unit north, east, and edge-normal vectors
    !> \author Michael Duda
    !> \date   15 January 2020
    !> \details
    !>  This routine computes the local unit north and east vectors at all cell
    !>  centers, storing the resulting fields in the mesh pool as 'north' and
    !>  'east'. It also computes the edge-normal unit vectors by calling
    !>  the mpas_initialize_vectors routine. Before this routine is called,
    !>  the mesh pool must contain 'latCell' and 'lonCell' fields that are valid
    !>  for all cells (not just solve cells), plus any fields that are required
    !>  by the mpas_initialize_vectors routine.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_compute_unit_vectors()

       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_derived_types, only : mpas_pool_type
       use mpas_kind_types, only : RKIND
       use mpas_vector_operations, only : mpas_initialize_vectors

       implicit none

       type (mpas_pool_type), pointer :: meshPool
       real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real(kind=RKIND), dimension(:,:), pointer :: east, north
       integer, pointer :: nCells
       integer :: iCell

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'east', east)
       call mpas_pool_get_array(meshPool, 'north', north)

       do iCell = 1, nCells

          east(1,iCell) = -sin(lonCell(iCell))
          east(2,iCell) =  cos(lonCell(iCell))
          east(3,iCell) =  0.0

          ! Normalize
          east(1:3,iCell) = east(1:3,iCell) / sqrt(sum(east(1:3,iCell) * east(1:3,iCell)))

          north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
          north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
          north(3,iCell) =  cos(latCell(iCell))

          ! Normalize
          north(1:3,iCell) = north(1:3,iCell) / sqrt(sum(north(1:3,iCell) * north(1:3,iCell)))

       end do

       call mpas_initialize_vectors(meshPool)

    end subroutine cam_mpas_compute_unit_vectors


    !-----------------------------------------------------------------------
    !  routine cam_mpas_update_halo
    !
    !> \brief  Updates the halo of the named field
    !> \author Michael Duda
    !> \date   16 January 2020
    !> \details
    !>  Given the name of a field that is defined in the MPAS Registry.xml file,
    !>  this routine updates the halo for that field.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_update_halo(fieldName)

       use mpas_derived_types, only : field1DReal, field2DReal, field3DReal, field4DReal, field5DReal, &
                                      field1DInteger, field2DInteger, field3DInteger, &
                                      mpas_pool_field_info_type, MPAS_POOL_REAL, MPAS_POOL_INTEGER
       use mpas_pool_routines, only : MPAS_pool_get_field_info, MPAS_pool_get_field
       use mpas_dmpar, only : MPAS_dmpar_exch_halo_field

       implicit none

       character(len=*), intent(in) :: fieldName

       type (mpas_pool_field_info_type) :: fieldInfo
       type (field1DReal), pointer :: field_real1d
       type (field2DReal), pointer :: field_real2d
       type (field3DReal), pointer :: field_real3d
       type (field4DReal), pointer :: field_real4d
       type (field5DReal), pointer :: field_real5d
       type (field1DInteger), pointer :: field_int1d
       type (field2DInteger), pointer :: field_int2d
       type (field3DInteger), pointer :: field_int3d


       call MPAS_pool_get_field_info(domain_ptr % blocklist % allFields, trim(fieldName), fieldInfo)

       if (fieldInfo % fieldType == MPAS_POOL_REAL) then
           if (fieldInfo % nDims == 1) then
               nullify(field_real1d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_real1d)
               if (associated(field_real1d)) then
                   call MPAS_dmpar_exch_halo_field(field_real1d)
               end if
           else if (fieldInfo % nDims == 2) then
               nullify(field_real2d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_real2d)
               if (associated(field_real2d)) then
                   call MPAS_dmpar_exch_halo_field(field_real2d)
               end if
           else if (fieldInfo % nDims == 3) then
               nullify(field_real3d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_real3d)
               if (associated(field_real3d)) then
                   call MPAS_dmpar_exch_halo_field(field_real3d)
               end if
           else if (fieldInfo % nDims == 4) then
               nullify(field_real4d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_real4d)
               if (associated(field_real4d)) then
                   call MPAS_dmpar_exch_halo_field(field_real4d)
               end if
           else if (fieldInfo % nDims == 5) then
               nullify(field_real5d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_real5d)
               if (associated(field_real5d)) then
                   call MPAS_dmpar_exch_halo_field(field_real5d)
               end if
           else
               ! Error...
           end if
       else if (fieldInfo % fieldType == MPAS_POOL_INTEGER) then
           if (fieldInfo % nDims == 1) then
               nullify(field_int1d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_int1d)
               if (associated(field_int1d)) then
                   call MPAS_dmpar_exch_halo_field(field_int1d)
               end if
           else if (fieldInfo % nDims == 2) then
               nullify(field_int2d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_int2d)
               if (associated(field_int2d)) then
                   call MPAS_dmpar_exch_halo_field(field_int2d)
               end if
           else if (fieldInfo % nDims == 3) then
               nullify(field_int3d)
               call MPAS_pool_get_field(domain_ptr % blocklist % allFields, trim(fieldName), field_int3d)
               if (associated(field_int3d)) then
                   call MPAS_dmpar_exch_halo_field(field_int3d)
               end if
           else
               ! Error...
           end if
       else
           ! Error...
       end if

    end subroutine cam_mpas_update_halo


    !-----------------------------------------------------------------------
    !  routine cam_mpas_cell_to_edge_winds
    !
    !> \brief  Projects cell-centered winds to the normal component of velocity on edges
    !> \author Michael Duda
    !> \date   16 January 2020
    !> \details
    !>  Given zonal and meridional winds at cell centers, unit vectors in the east
    !>  and north directions at cell centers, and unit vectors in the normal
    !>  direction at edges, this routine projects the cell-centered winds onto
    !>  the normal vectors.
    !>
    !>  Prior to calling this routine, the halos for the zonal and meridional
    !>  components of cell-centered winds should be updated. It is also critical
    !>  that the east, north, uZonal, and uMerid field are all allocated with
    !>  a "garbage" element; this is handled automatically for fields allocated
    !>  by the MPAS infrastructure.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_cell_to_edge_winds(nEdges, uZonal, uMerid, east, north, edgeNormalVectors, &
                                           cellsOnEdge, uNormal)

       use mpas_kind_types, only : RKIND

       implicit none

       integer, intent(in) :: nEdges
       real(kind=RKIND), dimension(:,:), intent(in) :: uZonal, uMerid
       real(kind=RKIND), dimension(:,:), intent(in) :: east, north, edgeNormalVectors
       integer, dimension(:,:), intent(in) :: cellsOnEdge
       real(kind=RKIND), dimension(:,:), intent(out) :: uNormal

       integer :: iEdge, cell1, cell2


       do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          uNormal(:,iEdge) =  uZonal(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                            + uMerid(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)   &
                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell1)   &
                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell1))  &
                            + uZonal(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                                    +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                                    +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                            + uMerid(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)   &
                                                    +  edgeNormalVectors(2,iEdge) * north(2,cell2)   &
                                                    +  edgeNormalVectors(3,iEdge) * north(3,cell2))
       end do

    end subroutine cam_mpas_cell_to_edge_winds


    !-----------------------------------------------------------------------
    !  routine cam_mpas_run
    !
    !> \brief  Integrate dynamical state for the specified length of time
    !> \author Michael Duda
    !> \date   29 February 2020
    !> \details
    !>  This routine calls the dynamical solver in a loop, with each iteration
    !>  of the loop stepping the dynamical state forward by one dynamics
    !>  time step and stopping after the state has been advanced by the time
    !>  interval specified by the integrationLength argument.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_run(integrationLength)

       use atm_core, only : atm_do_timestep
       use mpas_derived_types, only : MPAS_Time_type, MPAS_TimeInterval_type, mpas_pool_type
       use mpas_kind_types, only : StrKIND, RKIND
       use mpas_log, only : mpas_log_write
       use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_subpool, mpas_pool_shift_time_levels
       use mpas_timekeeping, only : mpas_advance_clock, mpas_get_clock_time, mpas_get_time, MPAS_NOW, &
                                    operator(.lt.), operator(+)
       use mpas_timer, only : mpas_timer_start, mpas_timer_stop

       implicit none

       integer :: ierr

       type (MPAS_TimeInterval_type) :: integrationLength

       real (kind=RKIND), pointer :: dt
       type (MPAS_Time_Type) :: currTime
       type (MPAS_Time_type) :: runUntilTime
       character(len=StrKIND) :: timeStamp
       type (mpas_pool_type), pointer :: state

       integer, save :: itimestep = 1

       ! Eventually, dt should be domain specific
       call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)

       ! During integration, time level 1 stores the model state at the beginning of the
       !   time step, and time level 2 stores the state advanced dt in time by timestep(...)

       currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
       runUntilTime = currTime + integrationLength

       do while (currTime < runUntilTime)
          call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
          call mpas_log_write('Dynamics timestep beginning at '//trim(timeStamp))

          call mpas_timer_start('time integration')
          call atm_do_timestep(domain_ptr, dt, itimestep)
          call mpas_timer_stop('time integration')   

          ! Move time level 2 fields back into time level 1 for next time step
          call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
          call mpas_pool_shift_time_levels(state)
         
          ! Advance clock before writing output
          itimestep = itimestep + 1
          call mpas_advance_clock(clock)
          currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
       end do

    end subroutine cam_mpas_run


    !-----------------------------------------------------------------------
    !  routine cam_mpas_finalize
    !
    !> \brief  Finalize the MPAS core and infrastructure
    !> \author Michael Duda
    !> \date   29 February 2020
    !> \details
    !>  This routine finalizes the MPAS-A dycore and any infrastructure that
    !>  was set-up during the simulation. The work here mirrors that done in
    !>  mpas_atm_core.F::atm_core_finalize(), except there is no need to finalize
    !>  the MPAS-A diagnostics framework or stand-alone physics modules.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_finalize()

       use mpas_decomp, only : mpas_decomp_destroy_decomp_list
       use mpas_timekeeping, only : mpas_destroy_clock
       use mpas_atm_threading, only : mpas_atm_threading_finalize

       implicit none

       integer :: ierr

       call mpas_destroy_clock(clock, ierr)
       call mpas_decomp_destroy_decomp_list(domain_ptr % decompositions)
       call mpas_atm_threading_finalize(domain_ptr % blocklist)

    end subroutine cam_mpas_finalize


    !-----------------------------------------------------------------------
    !  routine add_stream_attributes
    !
    !> \brief  Adds default attributes to all MPAS streams
    !> \author Michael Duda
    !> \date   14 May 2019
    !> \details
    !>  ...
    !
    !-----------------------------------------------------------------------
    subroutine add_stream_attributes(domain)

       use mpas_stream_manager, only : MPAS_stream_mgr_add_att
       use mpas_derived_types, only : MPAS_Pool_iterator_type
       use mpas_derived_types, only : MPAS_POOL_CONFIG, MPAS_POOL_REAL, MPAS_POOL_INTEGER, MPAS_POOL_CHARACTER, MPAS_POOL_LOGICAL
       use mpas_kind_types, only : RKIND, StrKIND
       use mpas_pool_routines, only : mpas_pool_begin_iteration, mpas_pool_get_next_member, mpas_pool_get_config

       implicit none

       type (domain_type), intent(inout) :: domain

       type (MPAS_Pool_iterator_type) :: itr
       integer, pointer :: intAtt
       logical, pointer :: logAtt
       character (len=StrKIND), pointer :: charAtt
       real (kind=RKIND), pointer :: realAtt
       character (len=StrKIND) :: histAtt

       integer :: local_ierr


       if (domain % dminfo % nProcs < 10) then
           write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       else if (domain % dminfo % nProcs < 100) then
           write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       else if (domain % dminfo % nProcs < 1000) then
           write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       else if (domain % dminfo % nProcs < 10000) then
           write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       else if (domain % dminfo % nProcs < 100000) then
           write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       else
           write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
       end if

       call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'source', domain % core % source)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

       call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'is_periodic', domain % is_periodic)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'x_period', domain % x_period)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'y_period', domain % y_period)
       ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
       ! add a shortened version.
!      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain %  parent_id)
       call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)

       call mpas_pool_begin_iteration(domain % configs)
       do while (mpas_pool_get_next_member(domain % configs, itr))

          if ( itr % memberType == MPAS_POOL_CONFIG) then

             if ( itr % dataType == MPAS_POOL_REAL ) then
                call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
                call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
             else if ( itr % dataType == MPAS_POOL_INTEGER ) then
                call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
                call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
             else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
                call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
                call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
             else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
                call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
                if (logAtt) then
                   call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
                else
                   call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
                end if
             end if

           end if
       end do

    end subroutine add_stream_attributes


end module cam_mpas_subdriver
