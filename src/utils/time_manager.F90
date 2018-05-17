module time_manager

! Provide CAM specific time management.  This is a wrapper layer for the ESMF
! time manager utility.

use ESMF

use shr_kind_mod,     only: r8 => shr_kind_r8, SHR_KIND_CS
use shr_cal_mod,      only: shr_cal_noleap, shr_cal_gregorian
use spmd_utils,       only: masterproc
use string_utils,     only: to_upper
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

! Public methods

public ::&
   timemgr_init,             &! time manager initialization
   advance_timestep,         &! increment the clocks current time
   get_step_size,            &! return step size in seconds
   get_nstep,                &! return timestep number
   get_curr_date,            &! return date components at end of current timestep
   get_prev_date,            &! return date components at beginning of current timestep
   get_start_date,           &! return components of the start date
   get_ref_date,             &! return components of the reference date
   get_perp_date,            &! return components of the perpetual date, and current time of day
   get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
   get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
   get_curr_calday,          &! return calendar day at end of current timestep
   get_calday,               &! return calendar day from input date
   is_first_step,            &! return true on first step of initial run
   is_first_restart_step,    &! return true on first step of restart or branch run
   timemgr_is_caltype,       &! return true if incoming calendar type string matches actual calendar type in use
   timemgr_get_calendar_cf,  &! return cf standard for calendar type
   is_end_curr_day,          &! return true on last timestep in current day
   is_end_curr_month,        &! return true on last timestep in current month
   is_last_step,             &! return true on last timestep
   is_perpetual,             &! return true if perpetual calendar is in use
   timemgr_datediff,         &! calculate difference between two time instants
   timemgr_time_ge,          &! check if time2 is later than or equal to time1
   timemgr_time_inc,         &! increment time instant by a given interval
   timemgr_set_date_time,    &! set the current date and time
   set_time_float_from_date, &! returns a float representation of time given  yr, mon, day, sec
   set_date_from_time_float   ! returns yr, mon, day, sec given time float


! Private module data

integer, parameter :: uninit_int = -999999999

integer :: dtime = uninit_int               ! timestep in seconds

character(len=32) :: calendar               ! Calendar type
logical :: tm_first_restart_step = .false.  ! true for first step of a restart or branch run
logical :: tm_perp_calendar = .false.       ! true when using perpetual calendar
integer :: cal_type = uninit_int            ! calendar type

! The target attribute for tm_cal is needed (at least by NAG) because there are
! pointers to this object inside ESMF_Time objects.
type(ESMF_Calendar), target :: tm_cal        ! calendar
type(ESMF_Clock)            :: tm_clock      ! Model clock   
type(ESMF_Time)             :: tm_perp_date  ! perpetual date

!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_init( &
   dtime_in, calendar_in, start_ymd, start_tod, ref_ymd, &
   ref_tod, stop_ymd, stop_tod, curr_ymd, curr_tod,      &
   perpetual_run, perpetual_ymd, initial_run)

   ! Initialize the time manager.

   ! Arguments
   integer,          intent(in) :: dtime_in       ! Coupling period (sec)
   character(len=*), intent(IN) :: calendar_in    ! Calendar type
   integer,          intent(IN) :: start_ymd      ! Start date (YYYYMMDD)
   integer,          intent(IN) :: start_tod      ! Start time of day (sec)
   integer,          intent(IN) :: ref_ymd        ! Reference date (YYYYMMDD)
   integer,          intent(IN) :: ref_tod        ! Reference time of day (sec)
   integer,          intent(IN) :: stop_ymd       ! Stop date (YYYYMMDD)
   integer,          intent(IN) :: stop_tod       ! Stop time of day (sec)
   integer,          intent(IN) :: curr_ymd       ! current date (YYYYMMDD)
   integer,          intent(IN) :: curr_tod       ! current time of day (sec)
   logical,          intent(IN) :: perpetual_run  ! If in perpetual mode or not
   integer,          intent(IN) :: perpetual_ymd  ! Perpetual date (YYYYMMDD)
   logical,          intent(in) :: initial_run    ! true => initial (or startup) run

   ! Local variables
   character(len=*), parameter :: sub = 'timemgr_init'
   integer :: rc                            ! return code
   type(ESMF_Time) :: start_date            ! start date for run
   type(ESMF_Time) :: stop_date             ! stop date for run
   type(ESMF_Time) :: curr_date             ! temporary date used in logic
   type(ESMF_Time) :: ref_date              ! reference date for time coordinate
   !----------------------------------------------------------------------------------------

   dtime = dtime_in

   calendar = trim(calendar_in)
   call init_calendar()

   ! create ESMF time instant objects
   start_date = TimeSetymd(start_ymd, start_tod, "start_date")
   stop_date  = TimeSetymd(stop_ymd, stop_tod, "stop_date")
   ref_date   = TimeSetymd(ref_ymd, ref_tod, "ref_date")
   curr_date  = TimeSetymd(curr_ymd, curr_tod, "curr_date")

   ! In order that an initial run start with nstep=0 (which is the criteria
   ! for the is_first_step() function to return true), the current time
   ! must be the same as the start time.  To allow the driver to be doing a
   ! continue run and at the same time force CAM into an initial run mode
   ! this code forces the start time to equal the current time, even though
   ! that may not be true for the data received from the driver.
   if (initial_run) then
      start_date = curr_date
   end if

   ! Initialize ESMF clock
   call initialize_clock(start_date, ref_date, curr_date, stop_date)

   if (.not. initial_run) then

      ! Advance the timestep.  Data from the restart file corresponds to the
      ! last timestep of the previous run.
      call advance_timestep()
      tm_first_restart_step = .true.
   end if

   ! Initialize date used for perpetual calendar day calculation.
   if (perpetual_run) then
      tm_perp_calendar = .true.
      tm_perp_date = TimeSetymd(perpetual_ymd, 0, "tm_perp_date")
   end if

   ! Print configuration summary to log file (stdout).
   if (masterproc) then
      call timemgr_print()
   end if

end subroutine timemgr_init

!=========================================================================================

subroutine initialize_clock( start_date, ref_date, curr_date, stop_date )

! Create an ESMF clock based on the stepsize, start_date, stop_date, and ref_date
! Then advance the clock to the curr_date.

   ! Input variables
   type(ESMF_Time), intent(inout) :: start_date  ! start date for run
   type(ESMF_Time), intent(inout) :: ref_date    ! reference date for time coordinate
   type(ESMF_Time), intent(inout) :: curr_date   ! current date (equal to start_date)
   type(ESMF_Time), intent(inout) :: stop_date   ! stop date for run

   ! Local variables
   character(len=*), parameter :: sub = 'initialize_clock'
   type(ESMF_TimeInterval) :: step_size       ! timestep size
   type(ESMF_Time) :: current     ! current date (from clock)
   integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
   integer :: rc                  ! return code

   if ( mod(86400,dtime) /= 0 ) then
      call endrun (sub//': timestep must divide evenly into 1 day')
   end if

   call ESMF_TimeIntervalSet(step_size, s=dtime, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

   ! check for valid input
   if ( stop_date <= start_date ) then
      write(iulog,*)sub, ': stop date must be specified later than start date: '
      call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
      call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
      call endrun
   end if
   if ( curr_date >= stop_date ) then
      write(iulog,*)sub, ': stop date must be specified later than current date: '
      call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
      call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
      write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
      call endrun
   end if

   ! Initialize the clock
   tm_clock = ESMF_ClockCreate(name="CAM Time-manager clock", timeStep=step_size, startTime=start_date, &
                               stopTime=stop_date, refTime=ref_date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockSetup')

   ! Advance clock to the current time (in case of a branch or restart)
   call ESMF_ClockGet(tm_clock, currTime=current, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   do while( curr_date > current )
      call ESMF_ClockAdvance( tm_clock, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
      call ESMF_ClockGet(tm_clock, currTime=current )
      call chkrc(rc, sub//': error return from ESMF_ClockGet')
   end do

end subroutine initialize_clock

!=========================================================================================

function TimeSetymd( ymd, tod, desc )
!
! Set the time by an integer as YYYYMMDD and integer seconds in the day
!
   integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
   integer, intent(in) :: tod            ! Time of day in seconds
   character(len=*), intent(in) :: desc  ! Description of time to set

   type(ESMF_Time) :: TimeSetymd    ! Return value

   character(len=*), parameter :: sub = 'TimeSetymd'
   integer :: yr, mon, day          ! Year, month, day as integers
   integer :: rc                    ! return code

   if ( (ymd < 0) .or. (tod < 0) .or. (tod > 24*3600) )then
      write(iulog,*) sub//': error yymmdd is a negative number or time-of-day out of bounds', &
                 ymd, tod
      call endrun
   end if
   yr  = ymd / 10000
   mon = (ymd - yr*10000) / 100
   day =  ymd - yr*10000 - mon*100
   call ESMF_TimeSet( TimeSetymd, yy=yr, mm=mon, dd=day, s=tod, &
                      calendar=tm_cal, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeSet: setting '//trim(desc))
end function TimeSetymd

!=========================================================================================

subroutine set_time_float_from_date( time, year, month, day, sec )
!
! Set the time as a float given year, month, day, sec
!
  implicit none

  real(r8),intent(out):: time
  integer, intent(in) :: year
  integer, intent(in) :: month
  integer, intent(in) :: day
  integer, intent(in) :: sec

  integer :: rc                    ! return code
  character(len=*), parameter :: sub = 'set_time_float_from_date'
  type(ESMF_Time) :: date, ref_date
  type(ESMF_TimeInterval) :: diff
  integer :: useday

  call ESMF_TimeSet( date, yy=year, mm=month, dd=day, s=sec, calendar=tm_cal, rc=rc)
  !
  ! If the subroutine returned error, check if it is Feb 29 of a non-leap year
  ! (legitimately used by the time-interpolation routines in tracer_data.F90)
  ! in which case, substitute Feb 28 for the day
  !
  if ( rc .ne. ESMF_SUCCESS ) then
     if ( ( month .eq. 2 ) .and. ( day .eq. 29 ) ) then ! assume the failure is because it is leap day
        useday = 28
        call ESMF_TimeSet( date, yy=year, mm=month, dd=useday, s=sec, calendar=tm_cal, rc=rc)
     else  ! legitimate error, let the model quit
        call chkrc(rc, sub//': error return from ESMF_TimeSet for set_time_float_from_date')        
     endif 
  endif  

  call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet for set_time_float_from_date')

  diff = date - ref_date

  call ESMF_TimeIntervalGet( diff, d_r8=time, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet for set_time_float_from_date')

endsubroutine set_time_float_from_date

!=========================================================================================

subroutine set_date_from_time_float( time, year, month, day, sec )
!
! Set year, month, day, sec given the time as a float
!
  implicit none

  real(r8),intent(in)  :: time
  integer, intent(out) :: year
  integer, intent(out) :: month
  integer, intent(out) :: day
  integer, intent(out) :: sec

  integer :: rc                    ! return code
  character(len=*), parameter :: sub = 'set_date_from_time_float'
  type(ESMF_Time) :: date, ref_date
  type(ESMF_TimeInterval) :: diff

  call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet for rset_date_from_time_float')

  call ESMF_TimeIntervalSet( diff, d_r8=time, rc=rc)

  call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet for set_date_from_time_float')

  date = ref_date + diff

  call ESMF_TimeGet( date, yy=year, mm=month, dd=day, s=sec, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeGet for set_date_from_time_float')

endsubroutine set_date_from_time_float

!=========================================================================================

integer function TimeGetymd( date, tod )
!
! Get the date and time of day in ymd from ESMF Time.
!
  type(ESMF_Time), intent(inout) :: date    ! Input date to convert to ymd
  integer, intent(out), optional :: tod  ! Time of day in seconds

  character(len=*), parameter :: sub = 'TimeGetymd'
  integer :: yr, mon, day
  integer :: rc                          ! return code

  call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeGet')
  TimeGetymd = yr*10000 + mon*100 + day
  if ( present( tod ) )then
     call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
     call chkrc(rc, sub//': error return from ESMF_TimeGet')
  end if
  if ( yr < 0 )then
     write(iulog,*) sub//': error year is less than zero', yr
     call endrun
  end if
end function TimeGetymd

!=========================================================================================

subroutine timemgr_set_date_time( new_ymd, new_tod )

   implicit none

   integer,          intent(IN) :: new_ymd  ! date (YYYYMMDD)
   integer,          intent(IN) :: new_tod  ! time of day (sec)

   type(ESMF_Time) :: new_time   ! new time obj
   type(ESMF_Time) :: clk_time   ! clock's time obj
   integer :: rc                 ! return code

   character(len=*), parameter :: sub = 'timemgr_set_date_time'

   new_time = TimeSetymd( new_ymd, new_tod, "curr_date" )

! Advance clock to the new time

   call ESMF_ClockGet(tm_clock, currTime=clk_time, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   do while( new_time > clk_time )
      call ESMF_ClockAdvance( tm_clock, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
      call ESMF_ClockGet(tm_clock, currTime=clk_time )
      call chkrc(rc, sub//': error return from ESMF_ClockGet')
   end do
! Print configuration summary to log file (stdout).

   if (masterproc) then
      call timemgr_print()
   end if

end subroutine timemgr_set_date_time

!=========================================================================================

subroutine init_calendar( )
!
! Initialize calendar
!
! Local variables
   character(len=*), parameter :: sub = 'init_calendar'
   type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
   character(len=len(calendar)) :: caltmp
   integer :: rc                              ! return code

   caltmp = to_upper(trim(calendar) )
   if ( trim(caltmp) == trim(shr_cal_noleap) ) then
      cal_type = ESMF_CALKIND_NOLEAP
   else if ( trim(caltmp) == trim(shr_cal_gregorian) ) then
      cal_type = ESMF_CALKIND_GREGORIAN
   else
      write(iulog,*)sub,': unrecognized calendar specified: ',calendar
      call endrun
   end if
   tm_cal = ESMF_CalendarCreate( name=caltmp, calkindflag=cal_type, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_CalendarSet')
end subroutine init_calendar

!=========================================================================================

subroutine timemgr_print()

! Local variables
   character(len=*), parameter :: sub = 'timemgr_print'
   integer :: rc
   integer :: yr, mon, day
   integer ::&
      nstep     = uninit_int,  &! current step number
      step_sec  = uninit_int,  &! timestep size seconds
      start_yr  = uninit_int,  &! start year
      start_mon = uninit_int,  &! start month
      start_day = uninit_int,  &! start day of month
      start_tod = uninit_int,  &! start time of day
      stop_yr   = uninit_int,  &! stop year
      stop_mon  = uninit_int,  &! stop month
      stop_day  = uninit_int,  &! stop day of month
      stop_tod  = uninit_int,  &! stop time of day
      ref_yr    = uninit_int,  &! reference year
      ref_mon   = uninit_int,  &! reference month
      ref_day   = uninit_int,  &! reference day of month
      ref_tod   = uninit_int,  &! reference time of day
      curr_yr   = uninit_int,  &! current year
      curr_mon  = uninit_int,  &! current month
      curr_day  = uninit_int,  &! current day of month
      curr_tod  = uninit_int    ! current time of day
   integer(ESMF_KIND_I8) :: step_no
   type(ESMF_Time) :: start_date! start date for run
   type(ESMF_Time) :: stop_date ! stop date for run
   type(ESMF_Time) :: curr_date ! current date for run
   type(ESMF_Time) :: ref_date  ! reference date
   type(ESMF_TimeInterval) :: step ! Time-step
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, &
                       refTime=ref_date, stopTime=stop_date, timeStep=step, &
                       advanceCount=step_no, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   nstep = step_no

   write(iulog,*)' ********** Time Manager Configuration **********'

   call ESMF_TimeIntervalGet( step, s=step_sec, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

   call ESMF_TimeGet( start_date, yy=start_yr, mm=start_mon, dd=start_day, &
                      s=start_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( stop_date, yy=stop_yr, mm=stop_mon, dd=stop_day, &
                      s=stop_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( ref_date, yy=ref_yr, mm=ref_mon, dd=ref_day, s=ref_tod, &
                      rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeGet( curr_date, yy=curr_yr, mm=curr_mon, dd=curr_day, &
                      s=curr_tod, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

   write(iulog,*)'  Calendar type:            ',trim(calendar)
   write(iulog,*)'  Timestep size (seconds):  ', step_sec
   write(iulog,*)'  Start date (yr mon day tod):     ', start_yr, start_mon, &
                start_day, start_tod
   write(iulog,*)'  Stop date (yr mon day tod):      ', stop_yr, stop_mon, &
                stop_day, stop_tod
   write(iulog,*)'  Reference date (yr mon day tod): ', ref_yr, ref_mon, &
                ref_day, ref_tod
   write(iulog,*)'  Current step number:      ', nstep
   write(iulog,*)'  Current date (yr mon day tod):   ', curr_yr, curr_mon, &
                curr_day, curr_tod

   if ( tm_perp_calendar ) then
      call ESMF_TimeGet( tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_TimeGet')
      write(iulog,*)'  Use perpetual diurnal cycle date (yr mon day): ', &
                yr, mon, day
   end if

   write(iulog,*)' ************************************************'

end subroutine timemgr_print
!=========================================================================================

subroutine advance_timestep()

! Increment the timestep number.

! Local variables
   character(len=*), parameter :: sub = 'advance_timestep'
   integer :: rc
!-----------------------------------------------------------------------------------------

   call ESMF_ClockAdvance( tm_clock, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockAdvance')

! Set first step flag off.

   tm_first_restart_step = .false.

end subroutine advance_timestep
!=========================================================================================

integer function get_step_size()

! Return the step size in seconds.

! Local variables
   character(len=*), parameter :: sub = 'get_step_size'
   type(ESMF_TimeInterval) :: step_size       ! timestep size
   integer :: rc
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, timeStep=step_size, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')

end function get_step_size
!=========================================================================================

integer function get_nstep()

! Return the timestep number.

! Local variables
   character(len=*), parameter :: sub = 'get_nstep'
   integer :: rc
   integer(ESMF_KIND_I8) :: step_no
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   get_nstep = step_no

end function get_nstep
!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

! Return date components valid at end of current timestep with an optional
! offset (positive or negative) in seconds.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

! Local variables
   character(len=*), parameter :: sub = 'get_curr_date'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: off
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_curr_date
!=========================================================================================

subroutine get_perp_date(yr, mon, day, tod, offset)

! Return time of day valid at end of current timestep and the components
! of the perpetual date (with an optional offset (positive or negative) in seconds.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

! Local variables
   character(len=*), parameter :: sub = 'get_perp_date'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: DelTime
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   ! Get time of day add it to perpetual date
   ! Get year, month, day so that seconds are time-of-day rather than since start time
   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   call ESMF_TimeIntervalSet(DelTime, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
   date = tm_perp_date + DelTime
   if ( present(offset) )then
      call ESMF_TimeIntervalSet(DelTime, s=offset, rc=rc)
      call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
      date = date + DelTime
   end if
   ! Get time of day from the result
   ! Get year, month, day so that seconds are time-of-day rather than since start time
   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)

   ! Get the date from the fixed perpetual date (in case it overflows to next day)
   call ESMF_TimeGet(tm_perp_date, yy=yr, mm=mon, dd=day, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')


end subroutine get_perp_date
!=========================================================================================

subroutine get_prev_date(yr, mon, day, tod)

! Return date components valid at beginning of current timestep.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_prev_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_prev_date
!=========================================================================================

subroutine get_start_date(yr, mon, day, tod)

! Return date components valid at beginning of initial run.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_start_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_start_date
!=========================================================================================

subroutine get_ref_date(yr, mon, day, tod)

! Return date components of the reference date.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'get_ref_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, refTime=date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_ref_date
!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'get_curr_time'
   integer :: rc
   type(ESMF_Time) :: cdate, rdate
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=cdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_ClockGet( tm_clock, refTime=rdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   diff = cdate - rdate

   call ESMF_TimeIntervalGet(diff, d=days, s=seconds, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

end subroutine get_curr_time
!=========================================================================================

subroutine get_prev_time(days, seconds)

! Return time components valid at beg of current timestep.
! prev time is the time interval between the prev date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'get_prev_time'
   integer :: rc
   type(ESMF_Time) :: date, ref_date
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet for prevTime')
   call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet for refTime')
   diff = date - ref_date
   call ESMF_TimeIntervalGet( diff, d=days, s=seconds, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeintervalGet')

end subroutine get_prev_time
!=========================================================================================

function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
   real(r8) :: get_curr_calday

! Local variables
   character(len=*), parameter :: sub = 'get_curr_calday'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: off, diurnal
   integer :: year, month, day, tod
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if

   if ( tm_perp_calendar ) then
!     Get current time-of-day from clock
      call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
      call chkrc(rc, sub//': error return from ESMF_TimeGet')
!     Get date from perpetual date add time-of-day to it
      call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
      call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
      date = tm_perp_date + diurnal
!!!!  write(iulog,*) ' tod = ', tod
!!!!  call ESMF_TimePrint( date, "string" )
   end if

   call ESMF_TimeGet( date, dayOfYear_r8=get_curr_calday, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

!
!                  WARNING: Gregorian calendar fakes day 366
!
! The zenith angle calculation is only capable of using a 365-day calendar.
! If a Gregorian calendar is being used, the last day of a leap year (day 366)
! is sent to the model as a repetition of the previous day (day 365). 
! This is done by decrementing calday by 1 immediately below.
! bundy, July 2008
!
   if (( get_curr_calday > 366.0_r8 ) .and. ( get_curr_calday <= 367.0_r8 ) &
        .and. (timemgr_is_caltype(trim(shr_cal_gregorian)))) then
      get_curr_calday = get_curr_calday - 1.0_r8
   endif

   if ( (get_curr_calday < 1.0_r8) .or. (get_curr_calday > 366.0_r8) )then
      write(iulog,*) 'atm '//sub//' calday = ', get_curr_calday
      if ( present(offset) ) write(iulog,*) 'offset = ', offset
      call endrun( sub//': error get_curr_calday out of bounds' )
   end if

end function get_curr_calday
!=========================================================================================

function get_calday(ymd, tod)

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, intent(in) :: &
      ymd,   &! date in yearmmdd format
      tod     ! time of day (seconds past 0Z)

! Return value
   real(r8) :: get_calday

! Local variables
   character(len=*), parameter :: sub = 'get_calday'
   integer :: rc                 ! return code
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   date = TimeSetymd( ymd, tod, "get_calday" )
   call ESMF_TimeGet( date, dayOfYear_r8=get_calday, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

!
!                  WARNING: Gregorian calendar fakes day 366
!
! The zenith angle calculation is only capable of using a 365-day calendar.
! If a Gregorian calendar is being used, the last day of a leap year (day 366)
! is sent to the model as a repetition of the previous day (day 365). 
! This is done by decrementing calday by 1 immediately below.
! bundy, July 2008
!
   if (( get_calday > 366.0_r8 ) .and. ( get_calday <= 367.0_r8 ) &
        .and. (timemgr_is_caltype(trim(shr_cal_gregorian)))) then
      get_calday = get_calday - 1.0_r8
   endif

   if ( (get_calday < 1.0_r8) .or. (get_calday > 366.0_r8) )then
      write(iulog,*) 'atm '//sub//' calday = ', get_calday
      call endrun( sub//': error calday out of range' )
   end if

end function get_calday
!=========================================================================================

character(len=SHR_KIND_CS) function timemgr_get_calendar_cf()

! Return cf standard for calendar type

! Local variables
   character(len=*), parameter  :: sub = 'timemgr_get_calendar_cf'
   character(len=len(calendar)) :: caltmp
!-----------------------------------------------------------------------------------------

   caltmp = to_upper(trim(calendar) )
   if ( trim(caltmp) == trim(shr_cal_noleap) ) then
      timemgr_get_calendar_cf = 'noleap'
   else if ( trim(caltmp) == trim(shr_cal_gregorian) ) then
      timemgr_get_calendar_cf = 'gregorian'
   else
      write(iulog,*)sub,': unrecognized calendar specified: ',calendar
      call endrun
   end if

end function timemgr_get_calendar_cf
!=========================================================================================
 
function timemgr_is_caltype( cal_in )

! Return true if incoming calendar type string matches actual calendar type in use

   character(len=*), intent(in) :: cal_in

! Return value
   logical :: timemgr_is_caltype

!-----------------------------------------------------------------------------------------

   timemgr_is_caltype = ( to_upper(trim(calendar)) == to_upper(trim(cal_in)) )

end function timemgr_is_caltype
!=========================================================================================
 
function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

! Return value
   logical :: is_end_curr_day

! Local variables
   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
!-----------------------------------------------------------------------------------------

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_day = (tod == 0)

end function is_end_curr_day
!=========================================================================================

logical function is_end_curr_month()

! Return true if current timestep is last timestep in current month.

! Local variables
   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)
!-----------------------------------------------------------------------------------------

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_month = (day == 1  .and.  tod == 0)

end function is_end_curr_month
!=========================================================================================

logical function is_first_step()

! Return true on first step of initial run only.

! Local variables
   character(len=*), parameter :: sub = 'is_first_step'
   integer :: rc
   integer :: nstep
   integer(ESMF_KIND_I8) :: step_no
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   nstep = step_no
   is_first_step = (nstep == 0)

end function is_first_step
!=========================================================================================

logical function is_first_restart_step()

! Return true on first step of restart run only.

!-----------------------------------------------------------------------------------------

   is_first_restart_step = tm_first_restart_step

end function is_first_restart_step
!=========================================================================================

logical function is_last_step()

! Return true on last timestep.

! Local variables
   character(len=*), parameter :: sub = 'is_last_step'
   type(ESMF_Time) :: stop_date
   type(ESMF_Time) :: curr_date
   type(ESMF_TimeInterval) :: time_step
   integer :: rc
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, stopTime=stop_date, &
                       currTime=curr_date, TimeStep=time_step, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   if ( curr_date+time_step > stop_date ) then
      is_last_step = .true.
   else
      is_last_step = .false.
   end if

end function is_last_step
!=========================================================================================

logical function is_perpetual()

! Return true on last timestep.

!-----------------------------------------------------------------------------------------

   is_perpetual = tm_perp_calendar

end function is_perpetual

!=========================================================================================

subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)

! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.

! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1,    &! time of day relative to date1 (seconds past 0Z)
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days

! Local variables
   character(len=*), parameter :: sub = 'timemgr_datediff'
   integer :: rc   ! return code

   type(ESMF_Time) :: date1
   type(ESMF_Time) :: date2
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   date1 = TimeSetymd( ymd1, tod1, "date1" )
   date2 = TimeSetymd( ymd2, tod2, "date2" )
   diff = date2 - date1
   call ESMF_TimeIntervalGet( diff, d_r8=days, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

end subroutine timemgr_datediff

!=========================================================================================

subroutine timemgr_time_ge(ymd1, tod1, ymd2, tod2, time2_ge_time1)

! time2_ge_time1 is set to true if (ymd2,tod2) is later than or equal to (ymd1,tod1)

! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1,    &! time of day relative to date1 (seconds past 0Z)
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   logical :: time2_ge_time1

! Local variables
   character(len=*), parameter :: sub = 'timemgr_time_ge'
   integer :: rc   ! return code

   type(ESMF_Time) :: time1, time2
!-----------------------------------------------------------------------------------------

   time1 = TimeSetymd( ymd1, tod1, "date1" )
   time2 = TimeSetymd( ymd2, tod2, "date2" )
   time2_ge_time1 = (time2 >= time1)

end subroutine timemgr_time_ge

!=========================================================================================

subroutine timemgr_time_inc(ymd1, tod1, ymd2, tod2, inc_s, inc_h, inc_d)

! Increment the time instant (ymd1,tod1) by an interval and return the resulting
! time instant (ymd2,tod2).

   ! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1      ! time of day relative to date1 (seconds past 0Z)

   integer, intent(out) ::&
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   integer, intent(in), optional ::&
      inc_s,   &! number of seconds in interval
      inc_h,   &! number of hours in interval
      inc_d     ! number of days in interval

   ! Local variables
   character(len=*), parameter :: sub = 'timemgr_time_inc'
   integer :: rc   ! return code

   type(ESMF_Time) :: date1
   type(ESMF_Time) :: date2
   type(ESMF_TimeInterval) :: t_interval
   integer :: year, month, day
!-----------------------------------------------------------------------------------------

   ! set esmf time object
   date1 = TimeSetymd( ymd1, tod1, "date1" )

   ! set esmf time interval object
   if (present(inc_s)) then
      call ESMF_TimeIntervalSet(t_interval, s=inc_s, rc=rc)
   else if (present(inc_h)) then
      call ESMF_TimeIntervalSet(t_interval, h=inc_h, rc=rc)
   else if (present(inc_d)) then
      call ESMF_TimeIntervalSet(t_interval, d=inc_d, rc=rc)
   else
      call endrun(sub//': one of the args inc_s, inc_h, or inc_d must be set')
   end if
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')

   ! increment the time instant
   date2 = date1 + t_interval

   ! extract the time components
   call ESMF_TimeGet(date2, yy=year, mm=month, dd=day, s=tod2, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')
   ymd2 = year*10000 + month*100 + day

end subroutine timemgr_time_inc

!=========================================================================================

subroutine chkrc(rc, mes)
   integer, intent(in)          :: rc   ! return code from time management library
   character(len=*), intent(in) :: mes  ! error message
   if ( rc == ESMF_SUCCESS ) return
   write(iulog,*) mes
   call endrun ('CHKRC')
end subroutine chkrc

end module time_manager
