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

    use cam_abortutils, only: endrun
    use mpas_derived_types, only : core_type, dm_info, domain_type, MPAS_Clock_type

    implicit none

    public :: cam_mpas_init_phase1, &
              cam_mpas_init_phase2, &
              cam_mpas_init_phase3, &
              cam_mpas_init_phase4, &
              cam_mpas_define_scalars, &
              cam_mpas_get_global_dims, &
              cam_mpas_get_global_coords, &
              cam_mpas_get_global_blocks, &
              cam_mpas_read_static, &
              cam_mpas_setup_restart, &
              cam_mpas_read_restart, &
              cam_mpas_write_restart, &
              cam_mpas_compute_unit_vectors, &
              cam_mpas_update_halo, &
              cam_mpas_cell_to_edge_winds, &
              cam_mpas_run, &
              cam_mpas_finalize, &
              cam_mpas_debug_stream, &
              cam_mpas_global_sum_real

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
    subroutine cam_mpas_init_phase1(mpicom, endrun, logUnits, realkind)

       use mpas_domain_routines, only : mpas_allocate_domain
       use mpas_framework, only : mpas_framework_init_phase1
       use atm_core_interface, only : atm_setup_core, atm_setup_domain
       use mpas_pool_routines, only : mpas_pool_add_config
       use mpas_kind_types, only : RKIND

       ! Dummy argument
       integer, intent(in) :: mpicom
       procedure(halt_model) :: endrun
       integer, dimension(2), intent(in) :: logUnits
       integer, intent(in) :: realkind

       ! Local variables
       integer :: ierr

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase1'


       allocate(corelist, stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate corelist array')
       nullify(corelist % next)

       allocate(corelist % domainlist, stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate corelist%domainlist%next array')
       nullify(corelist % domainlist % next)

       domain_ptr => corelist % domainlist
       domain_ptr % core => corelist

       call mpas_allocate_domain(domain_ptr)
       domain_ptr % domainID = 0


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

       ! CAM does not yet allow running the dycore at a different precision than
       ! the physics package.  Check that the real kinds are the same.
       if (realkind /= RKIND) then
          call endrun(subname//': FATAL: CAM and MPAS real kinds do not match')
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
       use pio_types, only : iosystem_desc_t

       use mpas_framework, only : mpas_framework_init_phase2
       use mpas_timer, only : mpas_timer_start

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

       call mpas_timer_start('total time')

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
    !>  the MPAS core, including the allocation of all fields managed by MPAS.
    !>  The num_scalars argument should be set to CAM's value for PCNST,
    !>  the number of constituents.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_init_phase3(fh_ini, num_scalars, endrun)

       use mpas_log, only : mpas_log_write
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
       use mpas_pool_routines, only : mpas_pool_add_config

       type (file_desc_t), intent(inout) :: fh_ini
       integer, intent(in) :: num_scalars
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

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase3'


       mesh_filename = 'external mesh file'

       !
       ! Adding a config named 'cam_pcnst' with the number of constituents will indicate to
       ! MPAS-A setup code that it is operating as a CAM dycore, and that it is necessary to
       ! allocate scalars separately from other Registry-defined fields
       !
       call mpas_pool_add_config(domain_ptr % configs, 'cam_pcnst', num_scalars)

       mesh_iotype = MPAS_IO_NETCDF  ! Not actually used
       call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype, pio_file_desc=fh_ini)

       !
       ! Finalize the setup of blocks and fields
       !
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
       use atm_time_integration, only : mpas_atm_dynamics_init

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

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_init_phase4'

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
       if ( ierr /= 0 ) then
          call endrun(subname//': failed to get MPAS_START_TIME')
       end if
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

       !
       ! Prepare the dynamics for integration
       !
       call mpas_atm_dynamics_init(domain_ptr)

    end subroutine cam_mpas_init_phase4


    !-----------------------------------------------------------------------
    !  routine cam_mpas_define_scalars
    !
    !> \brief  Define the names of constituents at run-time
    !> \author Michael Duda
    !> \date   21 May 2020
    !> \details
    !>  Given an array of constituent names, which must have size equal to the number
    !>  of scalars that were set in the call to cam_mpas_init_phase3, and given
    !>  a function to identify which scalars are moisture species, this routine defines
    !>  scalar constituents for the MPAS-A dycore.
    !>  Because the MPAS-A dycore expects all moisture constituents to appear in
    !>  a contiguous range of constituent indices, this routine may in general need
    !>  to reorder the constituents; to allow for mapping of indices between CAM
    !>  physics and the MPAS-A dycore, this routine returns index mapping arrays
    !>  mpas_from_cam_cnst and cam_from_mpas_cnst.
    !>  
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_define_scalars(block, mpas_from_cam_cnst, cam_from_mpas_cnst, ierr)

       use mpas_derived_types, only : block_type

       use mpas_derived_types, only : mpas_pool_type, field3dReal
       use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field, &
                                      mpas_pool_get_dimension, mpas_pool_add_dimension
       use mpas_attlist, only : mpas_add_att
       use mpas_log, only : mpas_log_write
       use mpas_derived_types, only : MPAS_LOG_ERR

       use constituents, only: cnst_name, cnst_is_a_water_species

       ! Arguments
       type (block_type), pointer :: block
       integer, dimension(:), pointer :: mpas_from_cam_cnst, cam_from_mpas_cnst
       integer, intent(out) :: ierr

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_define_scalars'

       integer :: i, j, timeLevs
       integer, pointer :: num_scalars
       integer :: num_moist
       integer :: idx_passive
       type (mpas_pool_type), pointer :: statePool
       type (mpas_pool_type), pointer :: tendPool
       type (field3dReal), pointer :: scalarsField
       character(len=128) :: tempstr
       character :: moisture_char


       ierr = 0

       !
       ! Define scalars
       !
       nullify(statePool)
       call mpas_pool_get_subpool(block % structs, 'state', statePool)

       if (.not. associated(statePool)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''state'' pool was not found.', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       nullify(num_scalars)
       call mpas_pool_get_dimension(statePool, 'num_scalars', num_scalars)

       !
       ! The num_scalars dimension should have been defined by atm_core_interface::atm_allocate_scalars, and
       ! if this dimension does not exist, something has gone wrong
       !
       if (.not. associated(num_scalars)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''num_scalars'' dimension does not exist in the ''state'' pool.', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       !
       ! If at runtime there are not num_scalars names in the array of constituent names provided by CAM,
       ! something has gone wrong
       !
       if (size(cnst_name) /= num_scalars) then
          call mpas_log_write(trim(subname)//': ERROR: The number of constituent names is not equal to the num_scalars dimension', &
                              messageType=MPAS_LOG_ERR)
          call mpas_log_write('size(cnst_name) = $i, num_scalars = $i', intArgs=[size(cnst_name), num_scalars], &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       !
       ! In CAM, the first scalar (if there are any) is always Q (specific humidity); if this is not
       ! the case, something has gone wrong
       !
       if (size(cnst_name) > 0) then
          if (trim(cnst_name(1)) /= 'Q') then
             call mpas_log_write(trim(subname)//': ERROR: The first constituent is not Q', messageType=MPAS_LOG_ERR)
             ierr = 1
             return
          end if
       end if

       !
       ! Determine which of the constituents are moisture species
       !
       allocate(mpas_from_cam_cnst(num_scalars), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate mpas_from_cam_cnst array')
       mpas_from_cam_cnst(:) = 0
       num_moist = 0
       do i = 1, size(cnst_name)
          if (cnst_is_a_water_species(cnst_name(i))) then
                num_moist = num_moist + 1
                mpas_from_cam_cnst(num_moist) = i
          end if
       end do

       !
       ! If CAM has no scalars, let the only scalar in MPAS be 'qv' (a moisture species)
       !
       if (num_scalars == 1 .and. size(cnst_name) == 0) then
          num_moist = 1
       end if

       !
       ! Assign non-moisture constituents to mpas_from_cam_cnst(num_moist+1:size(cnst_name))
       !
       idx_passive = num_moist + 1
       do i = 1, size(cnst_name)

          ! If CAM constituent i is not already mapped as a moist constituent
          if (.not. cnst_is_a_water_species(cnst_name(i))) then
                mpas_from_cam_cnst(idx_passive) = i
                idx_passive = idx_passive + 1
          end if
       end do

       !
       ! Create inverse map, cam_from_mpas_cnst
       !
       allocate(cam_from_mpas_cnst(num_scalars), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate cam_from_mpas_cnst array')
       cam_from_mpas_cnst(:) = 0

       do i = 1, size(cnst_name)
          cam_from_mpas_cnst(mpas_from_cam_cnst(i)) = i
       end do

       timeLevs = 2

       do i = 1, timeLevs
          nullify(scalarsField)
          call mpas_pool_get_field(statePool, 'scalars', scalarsField, timeLevel=i)

          if (.not. associated(scalarsField)) then
             call mpas_log_write(trim(subname)//': ERROR: The ''scalars'' field was not found in the ''state'' pool', &
                                 messageType=MPAS_LOG_ERR)
             ierr = 1
             return
          end if

          if (i == 1) call mpas_pool_add_dimension(statePool, 'index_qv', 1)
          scalarsField % constituentNames(1) = 'qv'
          call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg kg^{-1}')
          call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')

          do j = 2, size(cnst_name)
             scalarsField % constituentNames(j) = trim(cnst_name(mpas_from_cam_cnst(j)))
          end do

       end do

       call mpas_pool_add_dimension(statePool, 'moist_start', 1)
       call mpas_pool_add_dimension(statePool, 'moist_end', num_moist)

       !
       ! Print a tabular summary of the mapping between constituent indices
       !
       call mpas_log_write('')
       call mpas_log_write('  i MPAS constituent mpas_from_cam_cnst(i)       i CAM constituent  cam_from_mpas_cnst(i)')
       call mpas_log_write('------------------------------------------     ------------------------------------------')
       do i = 1, min(num_scalars, size(cnst_name))
          if (i <= num_moist) then
             moisture_char = '*'
          else
             moisture_char = ' '
          end if
          write(tempstr, '(i3,1x,a16,1x,i18,8x,i3,1x,a16,1x,i18)') i, trim(scalarsField % constituentNames(i))//moisture_char, &
                                                                   mpas_from_cam_cnst(i), &
                                                                   i, trim(cnst_name(i)), &
                                                                   cam_from_mpas_cnst(i)
          call mpas_log_write(trim(tempstr))
       end do
       call mpas_log_write('------------------------------------------     ------------------------------------------')
       call mpas_log_write('* = constituent used as a moisture species in MPAS-A dycore')
       call mpas_log_write('')


       !
       ! Define scalars_tend
       !
       nullify(tendPool)
       call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

       if (.not. associated(tendPool)) then
          call mpas_log_write(trim(subname)//': ERROR: The ''tend'' pool was not found.', &
                              messageType=MPAS_LOG_ERR)
          ierr = 1
          return
       end if

       timeLevs = 1

       do i = 1, timeLevs
          nullify(scalarsField)
          call mpas_pool_get_field(tendPool, 'scalars_tend', scalarsField, timeLevel=i)

          if (.not. associated(scalarsField)) then
             call mpas_log_write(trim(subname)//': ERROR: The ''scalars_tend'' field was not found in the ''tend'' pool', &
                                 messageType=MPAS_LOG_ERR)
             ierr = 1
             return
          end if

          if (i == 1) call mpas_pool_add_dimension(tendPool, 'index_qv', 1)
          scalarsField % constituentNames(1) = 'tend_qv'
          call mpas_add_att(scalarsField % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
          call mpas_add_att(scalarsField % attLists(1) % attList, 'long_name', 'Tendency of water vapor mixing ratio')

          do j = 2, size(cnst_name)
             scalarsField % constituentNames(j) = 'tend_'//trim(cnst_name(mpas_from_cam_cnst(j)))
          end do
       end do

       call mpas_pool_add_dimension(tendPool, 'moist_start', 1)
       call mpas_pool_add_dimension(tendPool, 'moist_end', num_moist)

    end subroutine cam_mpas_define_scalars


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

       real (kind=RKIND), dimension(:), intent(out) :: latCellGlobal
       real (kind=RKIND), dimension(:), intent(out) :: lonCellGlobal
       real (kind=RKIND), dimension(:), intent(out) :: areaCellGlobal

       integer :: iCell

       integer, pointer :: nCellsSolve
       integer, dimension(:), pointer :: indexToCellID

       type (mpas_pool_type), pointer :: meshPool
       integer :: nCellsGlobal,ierr

       real (kind=RKIND), dimension(:), pointer :: latCell
       real (kind=RKIND), dimension(:), pointer :: lonCell
       real (kind=RKIND), dimension(:), pointer :: areaCell
       real (kind=RKIND), dimension(:), pointer :: temp

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_get_global_coords'


       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

       call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

       ! check: size(latCellGlobal) ?= nCellsGlobal

       allocate(temp(nCellsGlobal), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate temp array')
       
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
       use string_utils, only: int2str

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
       integer :: ierr
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_get_global_blocks'

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

       !
       ! nCellsPerBlock
       !
       allocate(temp1d(size(nCellsPerBlock)), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate temp1d array at line:'//int2str(__LINE__))
       temp1d(:) = 0
       temp1d(domain_ptr % dminfo % my_proc_id + 1) = nCellsSolve

       call mpas_dmpar_max_int_array(domain_ptr % dminfo, size(temp1d), temp1d, nCellsPerBlock)

       deallocate(temp1d)

       !
       ! indexToBlockID
       !
       allocate(temp1d(size(indexToBlockID)), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate temp1d array at line:'//int2str(__LINE__))
       temp1d(:) = -1
       do iCell=1,nCellsSolve
          temp1d(indexToCellID(iCell)) = domain_ptr % dminfo % my_proc_id + 1   ! 1-based block indices?
       end do

       call mpas_dmpar_max_int_array(domain_ptr % dminfo, size(temp1d), temp1d, indexToBlockID)

       deallocate(temp1d)

       !
       ! localCellIDBlock
       !
       allocate(temp1d(size(localCellIDBlock)), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate temp1d array at line:'//int2str(__LINE__))
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
          owningBlock = indexToBlockID(iCell)
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

       use mpas_kind_types, only : StrKIND
       use mpas_io_streams, only : MPAS_createStream, MPAS_closeStream, MPAS_streamAddField, MPAS_readStream
       use mpas_derived_types, only : MPAS_IO_READ, MPAS_IO_NETCDF, MPAS_Stream_type, MPAS_pool_type, &
                                      field0DReal, field1DReal, field2DReal, field3DReal, field1DInteger, field2DInteger, &
                                      MPAS_STREAM_NOERR
       use mpas_pool_routines, only : MPAS_pool_get_subpool, MPAS_pool_get_field, MPAS_pool_create_pool, MPAS_pool_destroy_pool, &
                                      MPAS_pool_add_config
       use mpas_dmpar, only : MPAS_dmpar_exch_halo_field
       use mpas_stream_manager, only : postread_reindex

       ! Arguments
       type (file_desc_t), pointer :: fh_ini
       procedure(halt_model) :: endrun

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_read_static'

       character(len=StrKIND) :: errString

       integer :: ierr
       integer :: ierr_total
       type (MPAS_pool_type), pointer :: meshPool
       type (MPAS_pool_type), pointer :: reindexPool
       type (field1DReal), pointer :: latCell, lonCell, xCell, yCell, zCell
       type (field1DReal), pointer :: latEdge, lonEdge, xEdge, yEdge, zEdge
       type (field1DReal), pointer :: latVertex, lonVertex, xVertex, yVertex, zVertex
       type (field1DInteger), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
       type (field1DReal), pointer :: fEdge, fVertex
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
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to create static input stream.')
       end if

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

       call mpas_pool_get_field(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_field(meshPool, 'fVertex', fVertex)

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

       ierr_total = 0

       call MPAS_streamAddField(mesh_stream, latCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, lonCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, xCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, yCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, latEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, lonEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, xEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, yEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, latVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, lonVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, xVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, yVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, indexToCellID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, indexToEdgeID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, indexToVertexID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, fEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, fVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, areaCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, areaTriangle, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, dcEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, dvEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, angleEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, kiteAreasOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, weightsOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, meshDensity, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, nEdgesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, nEdgesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, cellsOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, edgesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, edgesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, cellsOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, verticesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, verticesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, edgesOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, cellsOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, cf1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, cf2, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, cf3, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, rdzw, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, dzu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, rdzu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, fzm, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, fzp, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, zgrid, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zxu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zz, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zb, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, zb3, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, deriv_two, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, cellTangentPlane, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, coeffs_reconstruct, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(mesh_stream, edgeNormalVectors, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, localVerticalUnitVectors, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, defc_a, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(mesh_stream, defc_b, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       if (ierr_total > 0) then
           write(errString, '(a,i0,a)') subname//': FATAL: Failed to add ', ierr_total, ' fields to static input stream.'
           call endrun(trim(errString))
       end if

       call MPAS_readStream(mesh_stream, 1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to read static input stream.')
       end if

       call MPAS_closeStream(mesh_stream, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to close static input stream.')
       end if

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

       call MPAS_dmpar_exch_halo_field(fEdge)
       call MPAS_dmpar_exch_halo_field(fVertex)

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
    !  routine cam_mpas_setup_restart
    !
    !> \brief  Set up a restart stream, but do not read or write the stream
    !> \author Michael Duda
    !> \date   21 July 2020
    !> \details
    !>  This routine prepares an MPAS stream with fields needed to restart
    !>  the MPAS-A dynamics. The stream will read or write from an existing
    !>  PIO file descriptor, and whether the stream is set up for reading
    !>  or writing depends on the direction argument, which must be set to
    !>  either MPAS_IO_READ or MPAS_IO_WRITE.
    !>
    !>  This routine does not actually read or write the stream. A subsequent
    !>  call to either cam_mpas_read_restart or cam_mpas_write_restart must
    !>  be made to do this.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_setup_restart(fh_rst, restart_stream, direction, endrun)

       use pio, only : file_desc_t

       use mpas_kind_types, only : StrKIND
       use mpas_io_streams, only : MPAS_createStream, MPAS_streamAddField, MPAS_writeStreamAtt
       use mpas_derived_types, only : MPAS_IO_NETCDF, MPAS_Stream_type, MPAS_pool_type, &
                                      field0DReal, field1DReal, field2DReal, field3DReal, &
                                      field1DInteger, field2DInteger, field0DChar, &
                                      MPAS_IO_WRITE, MPAS_STREAM_NOERR
       use mpas_pool_routines, only : MPAS_pool_get_field

       ! Arguments
       type (file_desc_t), intent(inout) :: fh_rst
       type (MPAS_Stream_type), intent(inout) :: restart_stream
       integer, intent(in) :: direction
       procedure(halt_model) :: endrun

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_setup_restart'

       character(len=StrKIND) :: errString

       integer :: ierr
       integer :: ierr_total
       type (MPAS_pool_type), pointer :: allFields
       type (field1DReal), pointer :: latCell, lonCell, xCell, yCell, zCell
       type (field1DReal), pointer :: latEdge, lonEdge, xEdge, yEdge, zEdge
       type (field1DReal), pointer :: latVertex, lonVertex, xVertex, yVertex, zVertex
       type (field1DInteger), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
       type (field1DReal), pointer :: fEdge, fVertex
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

       type (field0DChar), pointer :: initial_time
       type (field0DChar), pointer :: xtime
       type (field2DReal), pointer :: u
       type (field2DReal), pointer :: w
       type (field2DReal), pointer :: rho_zz
       type (field2DReal), pointer :: theta_m
       type (field3DReal), pointer :: scalars

       type (field1DReal), pointer :: meshScalingDel2
       type (field1DReal), pointer :: meshScalingDel4
       type (field2DReal), pointer :: dss
       type (field2DReal), pointer :: east
       type (field2DReal), pointer :: north
       type (field2DReal), pointer :: pressure_p
       type (field2DReal), pointer :: rho
       type (field2DReal), pointer :: theta
       type (field2DReal), pointer :: relhum
       type (field2DReal), pointer :: uReconstructZonal
       type (field2DReal), pointer :: uReconstructMeridional
       type (field2DReal), pointer :: circulation
       type (field2DReal), pointer :: exner
       type (field2DReal), pointer :: exner_base
       type (field2DReal), pointer :: rtheta_base
       type (field2DReal), pointer :: pressure_base
       type (field2DReal), pointer :: rho_base
       type (field2DReal), pointer :: theta_base
       type (field2DReal), pointer :: ru
       type (field2DReal), pointer :: ru_p
       type (field2DReal), pointer :: rw
       type (field2DReal), pointer :: rw_p
       type (field2DReal), pointer :: rtheta_p
       type (field2DReal), pointer :: rho_p
       type (field1DReal), pointer :: surface_pressure
       type (field2DReal), pointer :: t_init

       type (field1DReal), pointer :: u_init
       type (field1DReal), pointer :: qv_init


       call MPAS_createStream(restart_stream, domain_ptr % ioContext, 'not_used', MPAS_IO_NETCDF, &
                              direction, pio_file_desc=fh_rst, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to create restart stream.')
       end if

       allFields => domain_ptr % blocklist % allFields

       call mpas_pool_get_field(allFields, 'latCell', latCell)
       call mpas_pool_get_field(allFields, 'lonCell', lonCell)
       call mpas_pool_get_field(allFields, 'xCell', xCell)
       call mpas_pool_get_field(allFields, 'yCell', yCell)
       call mpas_pool_get_field(allFields, 'zCell', zCell)

       call mpas_pool_get_field(allFields, 'latEdge', latEdge)
       call mpas_pool_get_field(allFields, 'lonEdge', lonEdge)
       call mpas_pool_get_field(allFields, 'xEdge', xEdge)
       call mpas_pool_get_field(allFields, 'yEdge', yEdge)
       call mpas_pool_get_field(allFields, 'zEdge', zEdge)

       call mpas_pool_get_field(allFields, 'latVertex', latVertex)
       call mpas_pool_get_field(allFields, 'lonVertex', lonVertex)
       call mpas_pool_get_field(allFields, 'xVertex', xVertex)
       call mpas_pool_get_field(allFields, 'yVertex', yVertex)
       call mpas_pool_get_field(allFields, 'zVertex', zVertex)

       call mpas_pool_get_field(allFields, 'indexToCellID', indexToCellID)
       call mpas_pool_get_field(allFields, 'indexToEdgeID', indexToEdgeID)
       call mpas_pool_get_field(allFields, 'indexToVertexID', indexToVertexID)

       call mpas_pool_get_field(allFields, 'fEdge', fEdge)
       call mpas_pool_get_field(allFields, 'fVertex', fVertex)

       call mpas_pool_get_field(allFields, 'areaCell', areaCell)
       call mpas_pool_get_field(allFields, 'areaTriangle', areaTriangle)
       call mpas_pool_get_field(allFields, 'dcEdge', dcEdge)
       call mpas_pool_get_field(allFields, 'dvEdge', dvEdge)
       call mpas_pool_get_field(allFields, 'angleEdge', angleEdge)
       call mpas_pool_get_field(allFields, 'kiteAreasOnVertex', kiteAreasOnVertex)
       call mpas_pool_get_field(allFields, 'weightsOnEdge', weightsOnEdge)

       call mpas_pool_get_field(allFields, 'meshDensity', meshDensity)

       call mpas_pool_get_field(allFields, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_field(allFields, 'nEdgesOnEdge', nEdgesOnEdge)

       call mpas_pool_get_field(allFields, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_field(allFields, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_field(allFields, 'edgesOnEdge', edgesOnEdge)
       call mpas_pool_get_field(allFields, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_field(allFields, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_field(allFields, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_field(allFields, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_field(allFields, 'cellsOnVertex', cellsOnVertex)

       call mpas_pool_get_field(allFields, 'cf1', cf1)
       call mpas_pool_get_field(allFields, 'cf2', cf2)
       call mpas_pool_get_field(allFields, 'cf3', cf3)

       call mpas_pool_get_field(allFields, 'rdzw', rdzw)
       call mpas_pool_get_field(allFields, 'dzu', dzu)
       call mpas_pool_get_field(allFields, 'rdzu', rdzu)
       call mpas_pool_get_field(allFields, 'fzm', fzm)
       call mpas_pool_get_field(allFields, 'fzp', fzp)

       call mpas_pool_get_field(allFields, 'zgrid', zgrid)
       call mpas_pool_get_field(allFields, 'zxu', zxu)
       call mpas_pool_get_field(allFields, 'zz', zz)
       call mpas_pool_get_field(allFields, 'zb', zb)
       call mpas_pool_get_field(allFields, 'zb3', zb3)

       call mpas_pool_get_field(allFields, 'deriv_two', deriv_two)
       call mpas_pool_get_field(allFields, 'cellTangentPlane', cellTangentPlane)
       call mpas_pool_get_field(allFields, 'coeffs_reconstruct', coeffs_reconstruct)

       call mpas_pool_get_field(allFields, 'edgeNormalVectors', edgeNormalVectors)
       call mpas_pool_get_field(allFields, 'localVerticalUnitVectors', localVerticalUnitVectors)
       call mpas_pool_get_field(allFields, 'defc_a', defc_a)
       call mpas_pool_get_field(allFields, 'defc_b', defc_b)

       call mpas_pool_get_field(allFields, 'initial_time', initial_time, timeLevel=1)
       call mpas_pool_get_field(allFields, 'xtime', xtime, timeLevel=1)
       call mpas_pool_get_field(allFields, 'u', u, timeLevel=1)
       call mpas_pool_get_field(allFields, 'w', w, timeLevel=1)
       call mpas_pool_get_field(allFields, 'rho_zz', rho_zz, timeLevel=1)
       call mpas_pool_get_field(allFields, 'theta_m', theta_m, timeLevel=1)
       call mpas_pool_get_field(allFields, 'scalars', scalars, timeLevel=1)

       call mpas_pool_get_field(allFields, 'meshScalingDel2', meshScalingDel2)
       call mpas_pool_get_field(allFields, 'meshScalingDel4', meshScalingDel4)
       call mpas_pool_get_field(allFields, 'dss', dss)
       call mpas_pool_get_field(allFields, 'east', east)
       call mpas_pool_get_field(allFields, 'north', north)
       call mpas_pool_get_field(allFields, 'pressure_p', pressure_p)
       call mpas_pool_get_field(allFields, 'rho', rho)
       call mpas_pool_get_field(allFields, 'theta', theta)
       call mpas_pool_get_field(allFields, 'relhum', relhum)
       call mpas_pool_get_field(allFields, 'uReconstructZonal', uReconstructZonal)
       call mpas_pool_get_field(allFields, 'uReconstructMeridional', uReconstructMeridional)
       call mpas_pool_get_field(allFields, 'circulation', circulation)
       call mpas_pool_get_field(allFields, 'exner', exner)
       call mpas_pool_get_field(allFields, 'exner_base', exner_base)
       call mpas_pool_get_field(allFields, 'rtheta_base', rtheta_base)
       call mpas_pool_get_field(allFields, 'pressure_base', pressure_base)
       call mpas_pool_get_field(allFields, 'rho_base', rho_base)
       call mpas_pool_get_field(allFields, 'theta_base', theta_base)
       call mpas_pool_get_field(allFields, 'ru', ru)
       call mpas_pool_get_field(allFields, 'ru_p', ru_p)
       call mpas_pool_get_field(allFields, 'rw', rw)
       call mpas_pool_get_field(allFields, 'rw_p', rw_p)
       call mpas_pool_get_field(allFields, 'rtheta_p', rtheta_p)
       call mpas_pool_get_field(allFields, 'rho_p', rho_p)
       call mpas_pool_get_field(allFields, 'surface_pressure', surface_pressure)
       call mpas_pool_get_field(allFields, 't_init', t_init)

       call mpas_pool_get_field(allFields, 'u_init', u_init)
       call mpas_pool_get_field(allFields, 'qv_init', qv_init)

       ierr_total = 0

       call MPAS_streamAddField(restart_stream, latCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, lonCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, xCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, yCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, latEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, lonEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, xEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, yEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, latVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, lonVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, xVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, yVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, indexToCellID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, indexToEdgeID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, indexToVertexID, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, fEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, fVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, areaCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, areaTriangle, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, dcEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, dvEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, angleEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, kiteAreasOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, weightsOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, meshDensity, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, nEdgesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, nEdgesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, cellsOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, edgesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, edgesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, cellsOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, verticesOnCell, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, verticesOnEdge, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, edgesOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, cellsOnVertex, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, cf1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, cf2, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, cf3, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, rdzw, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, dzu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rdzu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, fzm, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, fzp, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, zgrid, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zxu, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zz, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zb, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, zb3, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, deriv_two, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, cellTangentPlane, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, coeffs_reconstruct, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, edgeNormalVectors, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, localVerticalUnitVectors, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, defc_a, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, defc_b, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, initial_time, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, xtime, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, u, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, w, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rho_zz, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, theta_m, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, scalars, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, meshScalingDel2, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, meshScalingDel4, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, dss, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, east, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, north, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, pressure_p, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rho, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, theta, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, relhum, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, uReconstructZonal, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, uReconstructMeridional, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, circulation, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, exner, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, exner_base, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rtheta_base, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, pressure_base, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rho_base, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, theta_base, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, ru, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, ru_p, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rw, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rw_p, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rtheta_p, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, rho_p, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, surface_pressure, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, t_init, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       call MPAS_streamAddField(restart_stream, u_init, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1
       call MPAS_streamAddField(restart_stream, qv_init, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) ierr_total = ierr_total + 1

       if (ierr_total > 0) then
           write(errString, '(a,i0,a)') subname//': FATAL: Failed to add ', ierr_total, ' fields to restart stream.'
           call endrun(trim(errString))
       end if

       if (direction == MPAS_IO_WRITE) then
          !
          ! Add global attributes to the stream
          !
          if (domain_ptr % on_a_sphere) then
             call MPAS_writeStreamAtt(restart_stream, 'on_a_sphere', 'YES')
          else
             call MPAS_writeStreamAtt(restart_stream, 'on_a_sphere', 'NO')
          end if
          call MPAS_writeStreamAtt(restart_stream, 'sphere_radius', domain_ptr % sphere_radius)
          if (domain_ptr % is_periodic) then
             call MPAS_writeStreamAtt(restart_stream, 'is_periodic', 'YES')
          else
             call MPAS_writeStreamAtt(restart_stream, 'is_periodic', 'NO')
          end if
          call MPAS_writeStreamAtt(restart_stream, 'x_period', domain_ptr % x_period)
          call MPAS_writeStreamAtt(restart_stream, 'y_period', domain_ptr % y_period)
          call MPAS_writeStreamAtt(restart_stream, 'parent_id', domain_ptr %  parent_id)
          call MPAS_writeStreamAtt(restart_stream, 'mesh_spec', domain_ptr % mesh_spec)
       end if

    end subroutine cam_mpas_setup_restart


    !-----------------------------------------------------------------------
    !  routine cam_mpas_read_restart
    !
    !> \brief  Reads a restart stream that was previously set up
    !> \author Michael Duda
    !> \date   22 July 2020
    !> \details
    !>  From a restart stream previously set up with a call to
    !>  cam_mpas_setup_restart, read the stream, update halos for all fields
    !>  that were read, and re-index mesh indexing fields (cellsOnCell,
    !>  edgesOnCell, etc.) from global to local index space.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_read_restart(restart_stream, endrun)

       use pio, only : file_desc_t

       use mpas_io_streams, only : MPAS_readStream, MPAS_closeStream
       use mpas_derived_types, only : MPAS_Stream_type, MPAS_pool_type, MPAS_STREAM_NOERR
       use mpas_pool_routines, only : MPAS_pool_create_pool, MPAS_pool_destroy_pool, MPAS_pool_add_config
       use mpas_stream_manager, only : postread_reindex

       ! Arguments
       type (MPAS_Stream_type), intent(inout) :: restart_stream
       procedure(halt_model) :: endrun

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_read_restart'

       integer :: ierr
       type (MPAS_pool_type), pointer :: reindexPool

       call MPAS_readStream(restart_stream, 1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to read restart stream.')
       end if

       call MPAS_closeStream(restart_stream, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to close restart stream.')
       end if

       !
       ! Perform halo updates for all decomposed fields (i.e., fields with
       ! an outermost dimension of nCells, nVertices, or nEdges)
       !
       call cam_mpas_update_halo('latCell', endrun)
       call cam_mpas_update_halo('lonCell', endrun)
       call cam_mpas_update_halo('xCell', endrun)
       call cam_mpas_update_halo('yCell', endrun)
       call cam_mpas_update_halo('zCell', endrun)

       call cam_mpas_update_halo('latEdge', endrun)
       call cam_mpas_update_halo('lonEdge', endrun)
       call cam_mpas_update_halo('xEdge', endrun)
       call cam_mpas_update_halo('yEdge', endrun)
       call cam_mpas_update_halo('zEdge', endrun)

       call cam_mpas_update_halo('latVertex', endrun)
       call cam_mpas_update_halo('lonVertex', endrun)
       call cam_mpas_update_halo('xVertex', endrun)
       call cam_mpas_update_halo('yVertex', endrun)
       call cam_mpas_update_halo('zVertex', endrun)

       call cam_mpas_update_halo('indexToCellID', endrun)
       call cam_mpas_update_halo('indexToEdgeID', endrun)
       call cam_mpas_update_halo('indexToVertexID', endrun)

       call cam_mpas_update_halo('fEdge', endrun)
       call cam_mpas_update_halo('fVertex', endrun)

       call cam_mpas_update_halo('areaCell', endrun)
       call cam_mpas_update_halo('areaTriangle', endrun)
       call cam_mpas_update_halo('dcEdge', endrun)
       call cam_mpas_update_halo('dvEdge', endrun)
       call cam_mpas_update_halo('angleEdge', endrun)
       call cam_mpas_update_halo('kiteAreasOnVertex', endrun)
       call cam_mpas_update_halo('weightsOnEdge', endrun)

       call cam_mpas_update_halo('meshDensity', endrun)

       call cam_mpas_update_halo('nEdgesOnCell', endrun)
       call cam_mpas_update_halo('nEdgesOnEdge', endrun)

       call cam_mpas_update_halo('cellsOnEdge', endrun)
       call cam_mpas_update_halo('edgesOnCell', endrun)
       call cam_mpas_update_halo('edgesOnEdge', endrun)
       call cam_mpas_update_halo('cellsOnCell', endrun)
       call cam_mpas_update_halo('verticesOnCell', endrun)
       call cam_mpas_update_halo('verticesOnEdge', endrun)
       call cam_mpas_update_halo('edgesOnVertex', endrun)
       call cam_mpas_update_halo('cellsOnVertex', endrun)

       call cam_mpas_update_halo('zgrid', endrun)
       call cam_mpas_update_halo('zxu', endrun)
       call cam_mpas_update_halo('zz', endrun)
       call cam_mpas_update_halo('zb', endrun)
       call cam_mpas_update_halo('zb3', endrun)

       call cam_mpas_update_halo('deriv_two', endrun)
       call cam_mpas_update_halo('cellTangentPlane', endrun)
       call cam_mpas_update_halo('coeffs_reconstruct', endrun)

       call cam_mpas_update_halo('edgeNormalVectors', endrun)
       call cam_mpas_update_halo('localVerticalUnitVectors', endrun)
       call cam_mpas_update_halo('defc_a', endrun)
       call cam_mpas_update_halo('defc_b', endrun)

       call cam_mpas_update_halo('u', endrun)
       call cam_mpas_update_halo('w', endrun)
       call cam_mpas_update_halo('rho_zz', endrun)
       call cam_mpas_update_halo('theta_m', endrun)
       call cam_mpas_update_halo('scalars', endrun)

       call cam_mpas_update_halo('meshScalingDel2', endrun)
       call cam_mpas_update_halo('meshScalingDel4', endrun)
       call cam_mpas_update_halo('dss', endrun)
       call cam_mpas_update_halo('east', endrun)
       call cam_mpas_update_halo('north', endrun)
       call cam_mpas_update_halo('pressure_p', endrun)
       call cam_mpas_update_halo('rho', endrun)
       call cam_mpas_update_halo('theta', endrun)
       call cam_mpas_update_halo('relhum', endrun)
       call cam_mpas_update_halo('uReconstructZonal', endrun)
       call cam_mpas_update_halo('uReconstructMeridional', endrun)
       call cam_mpas_update_halo('circulation', endrun)
       call cam_mpas_update_halo('exner', endrun)
       call cam_mpas_update_halo('exner_base', endrun)
       call cam_mpas_update_halo('rtheta_base', endrun)
       call cam_mpas_update_halo('pressure_base', endrun)
       call cam_mpas_update_halo('rho_base', endrun)
       call cam_mpas_update_halo('theta_base', endrun)
       call cam_mpas_update_halo('ru', endrun)
       call cam_mpas_update_halo('ru_p', endrun)
       call cam_mpas_update_halo('rw', endrun)
       call cam_mpas_update_halo('rw_p', endrun)
       call cam_mpas_update_halo('rtheta_p', endrun)
       call cam_mpas_update_halo('rho_p', endrun)
       call cam_mpas_update_halo('surface_pressure', endrun)
       call cam_mpas_update_halo('t_init', endrun)

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

       call postread_reindex(domain_ptr % blocklist % allFields, reindexPool)

       call MPAS_pool_destroy_pool(reindexPool)

    end subroutine cam_mpas_read_restart


    !-----------------------------------------------------------------------
    !  routine cam_mpas_write_restart
    !
    !> \brief  Writes a restart stream that was previously set up
    !> \author Michael Duda
    !> \date   22 July 2020
    !> \details
    !>  From a restart stream previously set up with a call to
    !>  cam_mpas_setup_restart, re-index mesh indexing fields (cellsOnCell,
    !>  edgesOnCell, etc.) from local to global index space, and write
    !>  the stream.
    !
    !-----------------------------------------------------------------------
    subroutine cam_mpas_write_restart(restart_stream, endrun)

       use pio, only : file_desc_t

       use mpas_io_streams, only : MPAS_writeStream, MPAS_closeStream
       use mpas_derived_types, only : MPAS_Stream_type, MPAS_pool_type, MPAS_STREAM_NOERR
       use mpas_pool_routines, only : MPAS_pool_create_pool, MPAS_pool_destroy_pool, MPAS_pool_add_config
       use mpas_stream_manager, only : prewrite_reindex, postwrite_reindex

       ! Arguments
       type (MPAS_Stream_type), intent(inout) :: restart_stream
       procedure(halt_model) :: endrun

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_write_restart'

       integer :: ierr
       type (MPAS_pool_type), pointer :: reindexPool

       !
       ! Re-index from local index space to global index space
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

       call prewrite_reindex(domain_ptr % blocklist % allFields, reindexPool)

       call MPAS_writeStream(restart_stream, 1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to write restart stream.')
       end if

       call postwrite_reindex(domain_ptr % blocklist % allFields, reindexPool)

       call MPAS_pool_destroy_pool(reindexPool)

       call MPAS_closeStream(restart_stream, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
           call endrun(subname//': FATAL: Failed to close restart stream.')
       end if

    end subroutine cam_mpas_write_restart


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

       type (mpas_pool_type), pointer :: meshPool
       real(kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real(kind=RKIND), dimension(:,:), pointer :: east, north
       integer, pointer :: nCells
       integer :: iCell

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_compute_unit_vectors'

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'east', east)
       call mpas_pool_get_array(meshPool, 'north', north)

       do iCell = 1, nCells

          east(1,iCell) = -sin(lonCell(iCell))
          east(2,iCell) =  cos(lonCell(iCell))
          east(3,iCell) =  0.0_RKIND

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
    subroutine cam_mpas_update_halo(fieldName, endrun)

       use mpas_derived_types, only : field1DReal, field2DReal, field3DReal, field4DReal, field5DReal, &
                                      field1DInteger, field2DInteger, field3DInteger, &
                                      mpas_pool_field_info_type, MPAS_POOL_REAL, MPAS_POOL_INTEGER
       use mpas_pool_routines, only : MPAS_pool_get_field_info, MPAS_pool_get_field
       use mpas_dmpar, only : MPAS_dmpar_exch_halo_field
       use mpas_kind_types, only : StrKIND

       ! Arguments
       character(len=*), intent(in) :: fieldName
       procedure(halt_model) :: endrun

       ! Local variables
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_update_halo'

       character(len=StrKIND) :: errString

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
               write(errString, '(a,i0,a)') subname//': FATAL: Unhandled dimensionality ', &
                                            fieldInfo % nDims, ' for real-valued field'
               call endrun(trim(errString))
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
               write(errString, '(a,i0,a)') subname//': FATAL: Unhandled dimensionality ', &
                                            fieldInfo % nDims, ' for integer-valued field'
               call endrun(trim(errString))
           end if
       else
           write(errString, '(a,i0,a)') subname//': FATAL: Unhandled field type ', fieldInfo % fieldType
           call endrun(trim(errString))
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

       integer, intent(in) :: nEdges
       real(kind=RKIND), dimension(:,:), intent(in) :: uZonal, uMerid
       real(kind=RKIND), dimension(:,:), intent(in) :: east, north, edgeNormalVectors
       integer, dimension(:,:), intent(in) :: cellsOnEdge
       real(kind=RKIND), dimension(:,:), intent(out) :: uNormal

       integer :: iEdge, cell1, cell2

       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_cell_to_edge_winds'

       do iEdge = 1, nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          uNormal(:,iEdge) =  uZonal(:,cell1)*0.5_RKIND*(edgeNormalVectors(1,iEdge)*east(1,cell1)   &
                                                       + edgeNormalVectors(2,iEdge)*east(2,cell1)   &
                                                       + edgeNormalVectors(3,iEdge)*east(3,cell1))  &
                            + uMerid(:,cell1)*0.5_RKIND*(edgeNormalVectors(1,iEdge)*north(1,cell1)  &
                                                       + edgeNormalVectors(2,iEdge)*north(2,cell1)  &
                                                       + edgeNormalVectors(3,iEdge)*north(3,cell1)) &
                            + uZonal(:,cell2)*0.5_RKIND*(edgeNormalVectors(1,iEdge)*east(1,cell2)   &
                                                       + edgeNormalVectors(2,iEdge)*east(2,cell2)   &
                                                       + edgeNormalVectors(3,iEdge)*east(3,cell2))  &
                            + uMerid(:,cell2)*0.5_RKIND*(edgeNormalVectors(1,iEdge)*north(1,cell2)  &
                                                       + edgeNormalVectors(2,iEdge)*north(2,cell2)  &
                                                       + edgeNormalVectors(3,iEdge)*north(3,cell2))
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
       use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_dimension, mpas_pool_get_array, &
                                      mpas_pool_get_subpool, mpas_pool_shift_time_levels
       use mpas_timekeeping, only : mpas_advance_clock, mpas_get_clock_time, mpas_get_time, MPAS_NOW, &
                                    operator(.lt.), operator(+)
#ifdef __NVCOMPILER
       !
       ! Some versions of the nvfortran compiler complain about the illegal use
       ! of an operator on a derived type if the following specific
       ! implementations of the < and + operators are not explicitly imported
       !
       use mpas_timekeeping, only : lt_t_t, add_t_ti
#endif
       use mpas_timer, only : mpas_timer_start, mpas_timer_stop
       use mpas_constants, only : Rv_over_Rd => rvord

       ! Arguments
       type (MPAS_TimeInterval_type), intent(in) :: integrationLength

       ! Local variables
       integer :: ierr

       real (kind=RKIND), pointer :: dt
       type (MPAS_Time_Type) :: currTime
       type (MPAS_Time_type) :: runUntilTime
       character(len=StrKIND) :: timeStamp
       type (mpas_pool_type), pointer :: state, diag, mesh

       integer, pointer :: index_qv
       integer, pointer :: nCellsSolve
       real(kind=RKIND), dimension(:,:), pointer :: theta_m, rho_zz, zz, theta, rho
       real(kind=RKIND), dimension(:,:,:), pointer :: scalars

       integer, save :: itimestep = 1
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_run'

       ! Eventually, dt should be domain specific
       call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_dt', dt)

       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state)
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag)
       call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh)

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

       !
       ! Compute diagnostic fields from the final prognostic state
       !
       call mpas_pool_get_dimension(state, 'index_qv', index_qv)
       call mpas_pool_get_dimension(state, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(state, 'theta_m', theta_m, timeLevel=1)
       call mpas_pool_get_array(state, 'rho_zz', rho_zz, timeLevel=1)
       call mpas_pool_get_array(state, 'scalars', scalars, timeLevel=1)
       call mpas_pool_get_array(mesh, 'zz', zz)
       call mpas_pool_get_array(diag, 'theta', theta)
       call mpas_pool_get_array(diag, 'rho', rho)

       rho(:,1:nCellsSolve) = rho_zz(:,1:nCellsSolve) * zz(:,1:nCellsSolve)
       theta(:,1:nCellsSolve) = theta_m(:,1:nCellsSolve) / (1.0_RKIND + Rv_over_Rd * scalars(index_qv,:,1:nCellsSolve))

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
       use mpas_timer, only : mpas_timer_write_header, mpas_timer_write, mpas_timer_finalize
       use mpas_log, only : mpas_log_finalize
       use mpas_timer, only : mpas_timer_stop
       use mpas_framework, only : mpas_framework_finalize
       use atm_time_integration, only : mpas_atm_dynamics_finalize

       ! Local variables
       integer :: ierr
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_finalize'


       !
       ! Finalize the dynamics
       !
       call mpas_atm_dynamics_finalize(domain_ptr)

       call mpas_destroy_clock(clock, ierr)
       call mpas_decomp_destroy_decomp_list(domain_ptr % decompositions)
       call mpas_atm_threading_finalize(domain_ptr % blocklist)

       call mpas_timer_stop('total time')

       call mpas_timer_write_header()
       call mpas_timer_write()
       call mpas_timer_finalize(domain_ptr)

       !
       ! Finalize infrastructure
       !

       ! Print out log stats and close log file
       !   (Do this after timer stats are printed and stream mgr finalized,
       !    but before framework is finalized because domain is destroyed there.)
       call mpas_log_finalize(ierr)

       call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr)

       deallocate(corelist % domainlist)
       deallocate(corelist)

    end subroutine cam_mpas_finalize


    subroutine cam_mpas_debug_stream(domain, filename, timeLevel)

       use mpas_io_streams, only : MPAS_createStream, MPAS_closeStream, MPAS_streamAddField, MPAS_writeStream
       use mpas_derived_types, only : MPAS_IO_WRITE, MPAS_IO_NETCDF, MPAS_STREAM_NOERR, MPAS_Stream_type, MPAS_pool_type, &
                                      field0DReal, field1DReal, field2DReal, field3DReal, field4DReal, field5DReal, &
                                      field1DInteger, field2DInteger, field3DInteger
       use mpas_pool_routines, only : MPAS_pool_get_subpool, MPAS_pool_get_field, MPAS_pool_create_pool, MPAS_pool_destroy_pool, &
                                      MPAS_pool_add_config

       use mpas_derived_types, only : MPAS_Pool_iterator_type, MPAS_POOL_FIELD, MPAS_POOL_REAL, MPAS_POOL_INTEGER
       use mpas_pool_routines, only : mpas_pool_begin_iteration, mpas_pool_get_next_member, mpas_pool_get_config

       type (domain_type), intent(inout) :: domain
       character(len=*), intent(in) :: filename
       integer, intent(in), optional :: timeLevel

       type (MPAS_Pool_iterator_type) :: itr

       integer :: ierr
       type (MPAS_pool_type), pointer :: allFields

       type (field0DReal), pointer :: field_real0d
       type (field1DReal), pointer :: field_real1d
       type (field2DReal), pointer :: field_real2d
       type (field3DReal), pointer :: field_real3d
       type (field4DReal), pointer :: field_real4d
       type (field5DReal), pointer :: field_real5d
       type (field1DInteger), pointer :: field_int1d
       type (field2DInteger), pointer :: field_int2d
       type (field3DInteger), pointer :: field_int3d

       type (MPAS_Stream_type) :: stream
       character(len=*), parameter :: subname = 'cam_mpas_subdriver::cam_mpas_debug_stream'


       call MPAS_createStream(stream, domain % ioContext, trim(filename), MPAS_IO_NETCDF, MPAS_IO_WRITE, &
                              clobberFiles=.true., clobberRecords=.true., truncateFiles=.true., ierr=ierr)

       allFields => domain % blocklist % allFields

       call mpas_pool_begin_iteration(allFields)
       do while (mpas_pool_get_next_member(allFields, itr))

          if (index(trim(itr % memberName), 'OwnedIndices') /= 0) then
             cycle
          end if

          if ( itr % memberType == MPAS_POOL_FIELD) then

             if (itr % dataType == MPAS_POOL_REAL) then
                if (itr % nDims == 0) then
                   nullify(field_real0d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real0d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real0d)
                   end if
                   if (associated(field_real0d)) then
                      call MPAS_streamAddField(stream, field_real0d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 1) then
                   nullify(field_real1d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real1d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real1d)
                   end if
                   if (associated(field_real1d)) then
                      call MPAS_streamAddField(stream, field_real1d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 2) then
                   nullify(field_real2d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real2d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real2d)
                   end if
                   if (associated(field_real2d)) then
                      call MPAS_streamAddField(stream, field_real2d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 3) then
                   nullify(field_real3d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real3d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real3d)
                   end if
                   if (associated(field_real3d)) then
                      call MPAS_streamAddField(stream, field_real3d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 4) then
                   nullify(field_real4d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real4d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real4d)
                   end if
                   if (associated(field_real4d)) then
                      call MPAS_streamAddField(stream, field_real4d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 5) then
                   nullify(field_real5d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real5d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_real5d)
                   end if
                   if (associated(field_real5d)) then
                      call MPAS_streamAddField(stream, field_real5d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                end if
             else if (itr % dataType == MPAS_POOL_INTEGER) then
                if (itr % nDims == 1) then
                   nullify(field_int1d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int1d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int1d)
                   end if
                   if (associated(field_int1d)) then
                      call MPAS_streamAddField(stream, field_int1d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 2) then
                   nullify(field_int2d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int2d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int2d)
                   end if
                   if (associated(field_int2d)) then
                      call MPAS_streamAddField(stream, field_int2d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                else if (itr % nDims == 3) then
                   nullify(field_int3d)
                   if (itr % nTimeLevels > 1) then
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int3d, timeLevel=timeLevel)
                   else
                      call mpas_pool_get_field(allFields, trim(itr % memberName), field_int3d)
                   end if
                   if (associated(field_int3d)) then
                      call MPAS_streamAddField(stream, field_int3d, ierr=ierr)
                      if (ierr /= MPAS_STREAM_NOERR) then
                         write(0,*) '*** Failed to add field '//trim(itr % memberName)
                      end if
                   else
                      write(0,*) '*** Failed to get field '//trim(itr % memberName)
                   end if
                end if
             end if

           end if
       end do

       call MPAS_writeStream(stream, 1, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
          write(0,*) '*** Error writing stream ', ierr
       end if

       call MPAS_closeStream(stream, ierr=ierr)
       if (ierr /= MPAS_STREAM_NOERR) then
          write(0,*) '*** Error closing stream ', ierr
       end if

    end subroutine cam_mpas_debug_stream


    !-----------------------------------------------------------------------
    !  routine cam_mpas_global_sum_real
    !
    !> \brief  Compute the global sum of real array
    !> \author Miles Curry
    !> \date   25 February 2021
    !> \details
    !>  This routine computes a global sum of a real array across all tasks
    !> in a communicator and returns that sum to all tasks.
    !>
    !
    !-----------------------------------------------------------------------
    function cam_mpas_global_sum_real(rarray) result(global_sum)

       use mpas_kind_types, only: RKIND
       use mpas_dmpar, only: mpas_dmpar_sum_real

       ! Input variables
       real (RKIND), dimension(:), intent(in) :: rarray
       real (RKIND) :: global_sum

       real (RKIND) :: local_sum

       local_sum = sum(rarray)
       call mpas_dmpar_sum_real(domain_ptr % dminfo, local_sum, global_sum)

    end function cam_mpas_global_sum_real


end module cam_mpas_subdriver
