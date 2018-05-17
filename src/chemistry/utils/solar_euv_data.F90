!-----------------------------------------------------------------------
! Solar EUV irradiance data
!-----------------------------------------------------------------------
module solar_euv_data
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use cam_pio_utils,    only: cam_pio_openfile
  use cam_logfile,      only: iulog
  use pio,              only: pio_get_var, pio_inq_varid, pio_inq_dimid, pio_inq_dimlen, &
                              file_desc_t
  use input_data_utils, only : time_coordinate

  implicit none

  save
  private
  public :: solar_euv_init
  public :: solar_euv_advance
  public :: solar_euv_data_etf
  public :: solar_euv_data_active

  real(r8), target, allocatable :: solar_euv_data_etf(:)
  logical, protected :: solar_euv_data_active = .false.

  integer :: nbins
  real(r8), allocatable :: irradi(:,:)

  type(file_desc_t) :: file_id
  integer :: ssi_vid

  logical :: initialized = .false.

  real(r8), allocatable :: dellam(:)
  real(r8), allocatable :: lambda(:)
  real(r8), allocatable :: we(:)

  integer, parameter :: nrecords = 2
  logical, parameter :: debug = .false.

  type(time_coordinate) :: time_coord

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_euv_init(filepath, fixed, fixed_ymd, fixed_tod)

    use ioFileMod, only : getfil

    ! arguments
    character(len=*), intent(in) :: filepath
    logical, intent(in) :: fixed
    integer, intent(in) :: fixed_ymd
    integer, intent(in) :: fixed_tod

    ! local variables
    integer :: astat, dimid, vid
    character(len=256) :: filen   

    integer :: ierr

    solar_euv_data_active = (filepath.ne.'NONE') 
    if ( .not.solar_euv_data_active ) return

    call time_coord%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, fixed_tod=fixed_tod )

    call getfil( filepath, filen, 0 )
    call cam_pio_openfile( file_id, filen, 0 )

    if(masterproc)  write(iulog,*)'solar_euv_data_init: data file = ',trim(filen)

    ierr = pio_inq_varid( file_id, 'ssi', ssi_vid )
    ierr = pio_inq_dimid( file_id, 'bin', dimid )
    ierr = pio_inq_dimlen( file_id, dimid, nbins )
    
    allocate(irradi(nbins,nrecords), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_euv_data_init: failed to allocate irradi; error = ',astat
       call endrun('solar_data_init')
    end if

    allocate(lambda(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_euv_data_init: failed to allocate lambda; error = ',astat
       call endrun('solar_euv_data_init')
    end if
    allocate(dellam(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_euv_data_init: failed to allocate dellam; error = ',astat
       call endrun('solar_euv_data_init')
    end if
    allocate(solar_euv_data_etf(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_euv_data_init: failed to allocate solar_euv_data_etf; error = ',astat
       call endrun('solar_euv_data_init')
    end if

    ierr = pio_inq_varid( file_id, 'wavelength', vid )
    ierr = pio_get_var( file_id, vid, lambda )
    ierr = pio_inq_varid( file_id, 'band_width', vid  )
    ierr = pio_get_var( file_id, vid, dellam )
    
    allocate(we(nbins+1), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_euv_data_init: failed to allocate we; error = ',astat
       call endrun('solar_euv_data_init')
    end if

    we(:nbins)  = lambda(:nbins) - 0.5_r8*dellam(:nbins)
    we(nbins+1) = lambda(nbins)  + 0.5_r8*dellam(nbins)

    deallocate(lambda)
    deallocate(dellam)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call solar_euv_advance()
    initialized = .true.

  end subroutine solar_euv_init

!-----------------------------------------------------------------------
! Reads in the ETF data for the current date.  
!-----------------------------------------------------------------------
  subroutine solar_euv_advance()

    integer  :: index
    logical  :: read_data
    integer  :: ierr
    integer  :: offset(2), count(2)
    real(r8) :: delt

    if (.not.solar_euv_data_active) return

    index = -1

    read_data = time_coord%read_more() .or. .not.initialized
    call time_coord%advance()

    if ( read_data ) then

       index = time_coord%indxs(1)

       ! get the surrounding time slices
       offset = (/ 1, index /)
       count =  (/ nbins, nrecords /)

       ierr = pio_get_var( file_id, ssi_vid, offset, count, irradi )
    endif

    delt = time_coord%wghts(2)

    solar_euv_data_etf(:) = irradi(:,1) + delt*( irradi(:,2) - irradi(:,1) )

  end subroutine solar_euv_advance

end module solar_euv_data
