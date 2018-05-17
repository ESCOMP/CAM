!================================================================================
! utility module for driver input data
!================================================================================
module drv_input_data

  use shr_kind_mod,     only: r8=>SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,      only: iulog
  use pio,              only: file_desc_t
  use time_manager,     only: get_step_size

  implicit none
  private
  save

  public :: drv_input_data_open
  public :: drv_input_data_read
  public :: drv_input_data_close
  public :: drv_input_data_freq
  public :: drv_input_data_t
  public :: drv_input_data_get

  public :: drv_input_4d_t
  public :: drv_input_3d_t
  public :: drv_input_2d_t
  public :: drv_input_2di_t

  interface drv_input_data_get
    module procedure get_data3d
    module procedure get_data2d
    module procedure get_idata2d
  end interface

  real(r8) :: drv_input_data_freq != nan
  
  type drv_input_data_t
     integer :: ntimes
     integer, allocatable :: dates(:)
     integer, allocatable :: secs(:)
     real(r8), allocatable :: times(:)
     type(file_desc_t) :: piofile
  endtype drv_input_data_t

  type drv_input_4d_t
     real(r8), pointer :: array(:,:,:)
  endtype drv_input_4d_t
  type drv_input_3d_t
     real(r8), pointer :: array(:,:)
  endtype drv_input_3d_t
  type drv_input_2d_t
     real(r8), pointer :: array(:)
  endtype drv_input_2d_t
  type drv_input_2di_t
     integer, pointer :: array(:)
  endtype drv_input_2di_t

  character(len=4) :: lonname = ' '
  character(len=4) :: latname = ' '

  interface drv_input_data_read
    module procedure drv_input_data_read_2d
    module procedure drv_input_data_read_3d
  end interface

contains

!=================================================================================
!=================================================================================
  subroutine drv_input_data_open( infile, indata )

    use cam_pio_utils, only: cam_pio_openfile
    use pio,           only: PIO_NOCLOBBER, pio_inq_dimid, pio_inq_dimlen
    use pio,           only: pio_inq_varid, pio_get_var
    use pio,           only: pio_seterrorhandling, PIO_INTERNAL_ERROR, PIO_BCAST_ERROR, PIO_NOERR
    use dyn_grid,      only: get_horiz_grid_dim_d

    implicit none

    character(len=*), intent(in) :: infile
    type(drv_input_data_t), intent(out) :: indata

    integer :: id, ierr
    integer :: hdim1_d,hdim2_d, nlons
    integer :: dtime
    integer :: data_dtime
    character(len=*), parameter :: sub = 'drv_input_data_open: '

    dtime = get_step_size()

    ! open file and get fileid
    !
    call cam_pio_openfile( indata%piofile, infile, PIO_NOCLOBBER)

    if(masterproc) write(iulog,*) sub // 'opened: ',trim(infile)

    !
    ! check horizontal grid ...
    !
    call pio_seterrorhandling( indata%piofile, PIO_BCAST_ERROR)
    lonname = 'ncol'
    latname =  ' '
    ierr = pio_inq_dimid( indata%piofile, lonname, id )
    if (ierr/=PIO_NOERR) then
       lonname = 'lon'
       latname = 'lat'
    endif

    ierr = pio_inq_dimid( indata%piofile, lonname, id )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find dimid for lonname')
    ierr = pio_inq_dimlen( indata%piofile, id, nlons )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find dimlen for lonname')

    call get_horiz_grid_dim_d(hdim1_d,hdim2_d)

    if (hdim1_d /= nlons) then
      call endrun('drv_input_data_open: input file has incorrect horizontal resolution')
    endif

    !
    ! get time/date info ...
    !
    ierr = pio_inq_dimid( indata%piofile, 'time', id )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find dimid for time')
    ierr = pio_inq_dimlen( indata%piofile, id, indata%ntimes )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find dimlen for time')
 
    allocate( indata%dates(indata%ntimes), indata%secs(indata%ntimes), indata%times(indata%ntimes) )

    ierr = pio_inq_varid( indata%piofile, 'date',  id  )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find varid for date')
    ierr = pio_get_var( indata%piofile, id, indata%dates )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to get values for date')

    ierr = pio_inq_varid( indata%piofile, 'datesec',  id  )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find varid for datesec')
    ierr = pio_get_var( indata%piofile, id, indata%secs )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to get values for datesec')

    ierr = pio_inq_varid( indata%piofile, 'time',  id  )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find varid for time')
    ierr = pio_get_var( indata%piofile, id, indata%times )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to get values for time')

    ierr = pio_inq_varid( indata%piofile, 'mdt',  id  )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to find varid for mdt')
    ierr = pio_get_var( indata%piofile, id, data_dtime )
    if (ierr/=PIO_NOERR) call endrun(sub//'failed to get value for mdt')

    call pio_seterrorhandling( indata%piofile, PIO_INTERNAL_ERROR)

    if ( .not. (data_dtime == dtime)) then
       write( iulog, * )  sub//'data mdt does not match dtime... use dtime = ', data_dtime
       call endrun(sub//'data mdt does not match dtime.')
    endif

  end subroutine drv_input_data_open

!================================================================================================
!================================================================================================
  subroutine drv_input_data_close(indata)
    use pio, only: pio_closefile
    implicit none

    type(drv_input_data_t), intent(inout) :: indata

    deallocate( indata%dates, indata%secs, indata%times )

    call pio_closefile( indata%piofile )

  end subroutine drv_input_data_close

  !=================================================================================
  !=================================================================================
  function drv_input_data_read_2d( indata, fldname, recno, abort ) result(field_array)
    use ncdio_atm,        only: infld

    implicit none

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*), intent(in) :: fldname
    integer,          intent(in) :: recno
    logical, optional,intent(in) :: abort

    logical  :: found, abort_run
    real(r8) :: field_array(pcols,begchunk:endchunk)

    abort_run = .false.
    if (present(abort)) then
       abort_run = abort
    endif

    call infld( fldname, indata%piofile, trim(lonname), trim(latname), 1,pcols, begchunk,endchunk, &
                field_array, found, gridname='physgrid',timelevel=recno)

    if (.not.found) then
       if ( abort_run ) then
          call endrun('drv_input_data_read_2d: did not find '// trim(fldname))
       else
          if (masterproc) write( iulog, * )  'drv_input_data_read_2d: ' // trim(fldname) // ' set to zero '
          field_array = 0._r8
       endif
    endif

  endfunction drv_input_data_read_2d

  !=================================================================================
  !=================================================================================
  function drv_input_data_read_3d( indata, fldname, vertname, vertsize, recno, abort ) result(field_array)
    use ncdio_atm,        only: infld
    implicit none

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*), intent(in) :: fldname
    character(len=*), intent(in) :: vertname
    integer,          intent(in) :: vertsize
    integer,          intent(in) :: recno
    logical, optional,intent(in) :: abort

    logical  :: found, abort_run
    real(r8) :: field_array(pcols,vertsize,begchunk:endchunk)

    real(r8), allocatable :: tmp_array(:,:,:)
    
    abort_run = .false.
    if (present(abort)) then
       abort_run = abort
    endif

    call infld( fldname, indata%piofile, lonname, vertname, latname, 1,pcols, 1,vertsize, begchunk,endchunk, &
                field_array, found, gridname='physgrid',timelevel=recno)

    if (.not.found) then
       if ( abort_run ) then
          call endrun('drv_input_data_read_3d: did not find '// trim(fldname))
       else
          if (masterproc) write( iulog, * )  'drv_input_data_read_3d: ' // trim(fldname) // ' set to zero '
          field_array = 0._r8
       endif
    endif

  endfunction drv_input_data_read_3d

  !================================================================================================
  !================================================================================================
  subroutine get_data3d(indata, infld_name, lev_name, nlev, recno, chunk_ptrs)

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*),       intent(in)    :: infld_name
    character(len=*),       intent(in)    :: lev_name
    integer,                intent(in)    :: nlev
    integer,                intent(in)    :: recno
    type(drv_input_3d_t),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:,:)

    integer :: c, ncol

    allocate( data (pcols, nlev,  begchunk:endchunk) )

    data = drv_input_data_read( indata, infld_name, lev_name, nlev, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:,:) = data(:,:,c)
    enddo

    deallocate( data )

  end subroutine get_data3d

  !================================================================================================
  !================================================================================================
  subroutine get_data2d(indata, infld_name, recno, chunk_ptrs)

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*),       intent(in)    :: infld_name
    integer,                intent(in)    :: recno
    type(drv_input_2d_t),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:)

    integer :: c, ncol

    allocate( data (pcols,  begchunk:endchunk) )

    data = drv_input_data_read( indata, infld_name, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:) = data(:,c)
    enddo

    deallocate( data )

  end subroutine get_data2d

  !================================================================================================
  !================================================================================================
  subroutine get_idata2d(indata, infld_name, recno, chunk_ptrs)

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*),       intent(in)    :: infld_name
    integer,                intent(in)    :: recno
    type(drv_input_2di_t),  intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8), allocatable :: data (:,:)

    integer :: c, ncol

    allocate( data (pcols,  begchunk:endchunk) )

    data = drv_input_data_read(indata,  infld_name, recno )
    do c=begchunk,endchunk
       chunk_ptrs(c)%array(:) = int(data(:,c))
    enddo

    deallocate( data )

  end subroutine get_idata2d

end module drv_input_data
