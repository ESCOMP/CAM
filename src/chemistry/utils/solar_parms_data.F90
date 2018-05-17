!-------------------------------------------------------------------------------
! solar variability parameters -- space weather indices
!-------------------------------------------------------------------------------
module solar_parms_data 

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use input_data_utils, only : time_coordinate
  use infnan,           only : nan, assignment(=)

  implicit none

  private
  save

 ! public interface

  public :: solar_parms_init
  public :: solar_parms_advance

  logical, public :: solar_parms_on = .false.

 ! time-interpolated quantities

  real(r8), public, protected :: solar_parms_f107
  real(r8), public, protected :: solar_parms_f107a
  real(r8), public, protected :: solar_parms_f107p ! previous day
  real(r8), public, protected :: solar_parms_kp
  real(r8), public, protected :: solar_parms_ap

 ! private data

  real(r8), allocatable :: f107_in(:)
  real(r8), allocatable :: f107a_in(:)
  real(r8), allocatable :: kp_in(:)
  real(r8), allocatable :: ap_in(:)

  type(time_coordinate) :: time_coord_curr ! for current model time interpolation
  type(time_coordinate) :: time_coord_prev ! for previous day time interpolation

contains

  subroutine solar_parms_init(filepath, fixed, fixed_ymd, fixed_tod)
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

    solar_parms_f107  = nan
    solar_parms_f107a = nan
    solar_parms_f107p = nan
    solar_parms_kp    = nan
    solar_parms_ap    = nan

    solar_parms_on = (filepath.ne.'NONE')

    if (.not.solar_parms_on) return

    !-----------------------------------------------------------------------
    !	... readin the solar parms dataset
    !-----------------------------------------------------------------------

    call getfil(filepath,  locfn, 0)
    call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)

    call time_coord_prev%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, fixed_tod=fixed_tod, &
                                     force_time_interp=.true., try_dates=.true., delta_days=-1._r8 )

    call time_coord_curr%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, fixed_tod=fixed_tod, &
                                     force_time_interp=.true., try_dates=.true. )

    !---------------------------------------------------------------
    !	... allocate and read solar parms
    !---------------------------------------------------------------
    allocate( f107_in(time_coord_curr%ntimes), f107a_in(time_coord_curr%ntimes), &
              kp_in(time_coord_curr%ntimes), ap_in(time_coord_curr%ntimes), stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'solar_parms_init', 'f107_in ... ap_in ', time_coord_curr%ntimes )
    end if
    ierr = pio_inq_varid( ncid, 'f107', varid )
    ierr = pio_get_var( ncid, varid, f107_in )
    ierr = pio_inq_varid( ncid, 'f107a', varid )
    ierr = pio_get_var( ncid, varid, f107a_in )
    ierr = pio_inq_varid( ncid, 'kp', varid )
    ierr = pio_get_var( ncid, varid, kp_in )
    ierr = pio_inq_varid( ncid, 'ap', varid )
    ierr = pio_get_var( ncid, varid, ap_in )

    call pio_closefile( ncid )

end subroutine solar_parms_init

subroutine solar_parms_advance
  !---------------------------------------------------------------
  ! time interpolate space wx indices
  !---------------------------------------------------------------

  integer  :: ndx1, ndx2
  real(r8) :: wgt1, wgt2

  if (solar_parms_on) then
     call time_coord_curr%advance()
     ndx1=time_coord_curr%indxs(1)
     ndx2=time_coord_curr%indxs(2)
     wgt1=time_coord_curr%wghts(1)
     wgt2=time_coord_curr%wghts(2)
     solar_parms_f107  = wgt1*f107_in(ndx1) + wgt2*f107_in(ndx2) 
     solar_parms_f107a = wgt1*f107a_in(ndx1) + wgt2*f107a_in(ndx2)
     solar_parms_kp    = wgt1*kp_in(ndx1) + wgt2*kp_in(ndx2)
     solar_parms_ap    = wgt1*ap_in(ndx1) + wgt2*ap_in(ndx2)

     call time_coord_prev%advance()
     ndx1=time_coord_prev%indxs(1)
     ndx2=time_coord_prev%indxs(2)
     wgt1=time_coord_prev%wghts(1)
     wgt2=time_coord_prev%wghts(2)
     solar_parms_f107p = wgt1*f107_in(ndx1) + wgt2*f107_in(ndx2)
  endif

end subroutine solar_parms_advance

end module solar_parms_data
