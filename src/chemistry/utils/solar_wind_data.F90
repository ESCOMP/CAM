!-------------------------------------------------------------------------------
! solar wind data -- IMF components, wind velocity and density
!-------------------------------------------------------------------------------
module solar_wind_data

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use input_data_utils, only : time_coordinate
  use infnan,           only : nan, assignment(=)

  implicit none

  private
  save

 ! public interface

  public :: solar_wind_init
  public :: solar_wind_advance

  logical, public :: solar_wind_on = .false.

 ! time-interpolated quantities

  real(r8), public, protected :: solar_wind_byimf
  real(r8), public, protected :: solar_wind_bzimf
  real(r8), public, protected :: solar_wind_swvel
  real(r8), public, protected :: solar_wind_swden

 ! private data

  real(r8), allocatable :: byimf_in(:)
  real(r8), allocatable :: bzimf_in(:)
  real(r8), allocatable :: swvel_in(:)
  real(r8), allocatable :: swden_in(:)

  type(time_coordinate) :: time_coord

contains

  subroutine solar_wind_init(filepath, fixed, fixed_ymd, fixed_tod)
    !---------------------------------------------------------------
    !	... initialize solar parmaters
    !---------------------------------------------------------------

    use ioFileMod
    use error_messages, only: alloc_err
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, &
                              pio_inq_varid, pio_closefile, pio_nowrite

    !---------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------
    character(len=*), intent(in) :: filepath
    logical, intent(in) :: fixed
    integer, intent(in) :: fixed_ymd
    integer, intent(in) :: fixed_tod

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    type(file_desc_t)  :: ncid
    type(var_desc_t)  :: varid
    integer  :: astat
    character(len=shr_kind_cl) :: locfn
    integer :: ierr

    solar_wind_byimf = nan
    solar_wind_bzimf = nan
    solar_wind_swvel = nan
    solar_wind_swden = nan

    solar_wind_on = (trim(filepath).ne.'NONE' .and. len_trim(filepath)>0)

    if (.not.solar_wind_on) return

    !-----------------------------------------------------------------------
    !	... readin the solar parms dataset
    !-----------------------------------------------------------------------

    call getfil(filepath,  locfn, 0)
    call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)

    call time_coord%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, &
                                fixed_tod=fixed_tod, force_time_interp=.true. )

    !---------------------------------------------------------------
    !	... allocate and read solar parms
    !---------------------------------------------------------------
    allocate( byimf_in(time_coord%ntimes), bzimf_in(time_coord%ntimes), &
              swvel_in(time_coord%ntimes), swden_in(time_coord%ntimes), &
              stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'solar_wind_init', 'byimf_in ... swden_in ', &
                       time_coord%ntimes )
    end if

    ierr = pio_inq_varid( ncid, 'by', varid )
    ierr = pio_get_var( ncid, varid, byimf_in )
    ierr = pio_inq_varid( ncid, 'bz', varid )
    ierr = pio_get_var( ncid, varid, bzimf_in )
    ierr = pio_inq_varid( ncid, 'swvel', varid )
    ierr = pio_get_var( ncid, varid, swvel_in )
    ierr = pio_inq_varid( ncid, 'swden', varid )
    ierr = pio_get_var( ncid, varid, swden_in )

    call pio_closefile( ncid )

end subroutine solar_wind_init

subroutine solar_wind_advance
  !---------------------------------------------------------------
  ! time interpolate space wx indices
  !---------------------------------------------------------------

  if (solar_wind_on) then
     call time_coord%advance()
    !  time interpolate
     solar_wind_byimf = time_coord%wghts(1)*byimf_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*byimf_in(time_coord%indxs(2))
     solar_wind_bzimf = time_coord%wghts(1)*bzimf_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*bzimf_in(time_coord%indxs(2))
     solar_wind_swvel = time_coord%wghts(1)*swvel_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*swvel_in(time_coord%indxs(2))
     solar_wind_swden = time_coord%wghts(1)*swden_in(time_coord%indxs(1)) &
                      + time_coord%wghts(2)*swden_in(time_coord%indxs(2))
  endif

end subroutine solar_wind_advance

end module solar_wind_data
