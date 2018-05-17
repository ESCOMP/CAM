!================================================================================
! offline unit driver utility module
!
!================================================================================
module offline_driver

  use shr_kind_mod,     only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog
  use drv_input_data,   only: drv_input_data_t
  use drv_input_data,   only: drv_input_data_open, drv_input_data_close
  use tracer_data,      only: incr_filename

  implicit none
  private
  save

  public :: offline_driver_init
  public :: offline_driver_run
  public :: offline_driver_dorun
  public :: offline_driver_readnl
  public :: offline_driver_reg

  integer :: recno = 1
  logical, protected :: offline_driver_dorun = .false.

  character(len=cl) :: current_file = ' '
  character(len=cl) :: next_file = ' '
  character(len=cl) :: offline_driver_fileslist = ' '

  type(drv_input_data_t) :: curr_indata

contains

!================================================================================
!================================================================================
  subroutine offline_driver_reg()
    use unit_driver, only: unit_driver_reg
    call unit_driver_reg()
  end subroutine offline_driver_reg

!================================================================================
!================================================================================
  subroutine offline_driver_run( phys_state, pbuf2d, cam_out, cam_in )

    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use camsrfexch,     only: cam_out_t, cam_in_t     
    use physics_buffer, only: physics_buffer_desc
    use time_manager,   only: get_curr_date
    use unit_driver,    only: unit_driver_run

    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    integer :: yr, mon, day
    integer :: curr_model_date, curr_model_tod
    logical :: active_step

    if (.not.offline_driver_dorun) return
    
    ! check model date/time against input data date/time
    call get_curr_date(yr, mon, day, curr_model_tod)
    curr_model_date = yr*10000 + mon*100 + day
    if ( recno <= curr_indata%ntimes ) then
       active_step = curr_model_date==curr_indata%dates(recno) .and. curr_model_tod==curr_indata%secs(recno)
    else
       active_step = .false.
    endif

    if (active_step) then

       call unit_driver_run(curr_indata, phys_state, pbuf2d, cam_out, cam_in, recno)

       recno = recno+1
       
       if ( recno > curr_indata%ntimes ) then
          call drv_input_data_close(curr_indata)
          current_file = next_file
          if ( current_file/='NOT_FOUND' ) then
             call drv_input_data_open( current_file, curr_indata )
             recno = 1
             next_file = incr_filename( current_file, filenames_list=offline_driver_fileslist, abort=.false.)
          endif
       endif

    endif

  end subroutine offline_driver_run

!================================================================================
!================================================================================
  subroutine offline_driver_init()
    use unit_driver,    only: unit_driver_init
    use drv_input_data, only: drv_input_data_freq
    use shr_const_mod,  only: SHR_CONST_CDAY
    use infnan,         only: nan, assignment(=)

    type(drv_input_data_t) :: next_indata

    if (.not.offline_driver_dorun) return

    call drv_input_data_open( current_file, curr_indata )

    drv_input_data_freq = nan

    if (curr_indata%ntimes>1) then
       drv_input_data_freq = (curr_indata%times(2) - curr_indata%times(1))*SHR_CONST_CDAY ! seconds
    else
       if ( next_file/='NOT_FOUND' ) then
          call drv_input_data_open( next_file, next_indata )
          drv_input_data_freq = (next_indata%times(1) - curr_indata%times(1))*SHR_CONST_CDAY ! seconds
          call drv_input_data_close(next_indata)
       endif
    endif

    call unit_driver_init()

  endsubroutine offline_driver_init

!=================================================================================
!=================================================================================
  subroutine offline_driver_readnl( nlfile )

    use namelist_utils,    only: find_group_name
    use units,             only: getunit, freeunit
    use cam_abortutils,    only: endrun
    use mpishorthand

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    character(len=cl) :: offline_driver_infile = ' '

    namelist /offline_driver_nl/ offline_driver_infile, offline_driver_fileslist

    if (masterproc) then

       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'offline_driver_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, offline_driver_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('offline_driver_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)

    end if

#ifdef SPMD
    call mpibcast (offline_driver_infile,     len(offline_driver_infile),    mpichar, 0, mpicom)
    call mpibcast (offline_driver_fileslist,  len(offline_driver_fileslist), mpichar, 0, mpicom)
#endif

    current_file = 'NOT_FOUND'
    next_file = 'NOT_FOUND'
    offline_driver_dorun = .false.
    
    if ( len_trim(offline_driver_infile) > 0 ) then
       current_file = trim(offline_driver_infile)
    elseif ( len_trim(offline_driver_fileslist) > 0 ) then
       current_file = incr_filename( offline_driver_infile, filenames_list=offline_driver_fileslist )
    else 
       offline_driver_dorun = .false.
       return
    endif

    if ( trim(current_file)/='NOT_FOUND' .and. len_trim(current_file) > 0 ) then
       offline_driver_dorun = .true.
       if ( len_trim(offline_driver_fileslist) > 0 ) then
          next_file = incr_filename( current_file, filenames_list=offline_driver_fileslist, abort=.false.)
       endif
    endif

  endsubroutine offline_driver_readnl

end module offline_driver
