module cam_control_mod
!------------------------------------------------------------------------------------------------
!
! High level control variables.  Information received from the driver/coupler is
! stored here.
!
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
use spmd_utils,       only: masterproc
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
public
save

! Public Routines:
!
!   cam_ctrl_init
!   cam_ctrl_set_orbit
!   cam_ctrl_set_physics_type

character(len=cl), protected :: caseid  ! case ID
character(len=cl), protected :: ctitle  ! case title

logical, protected :: initial_run  ! startup mode which only requires a minimal initial file
logical, protected :: restart_run  ! continue a previous run; requires a restart file
logical, protected :: branch_run   ! branch from a previous run; requires a restart file
logical, protected :: post_assim    ! We are resuming after a pause

logical, protected :: adiabatic         ! true => no physics
logical, protected :: ideal_phys        ! true => run Held-Suarez (1994) physics
logical, protected :: kessler_phys      ! true => run Kessler physics
logical, protected :: tj2016_phys       ! true => run tj2016 physics
logical, protected :: simple_phys       ! true => adiabatic or ideal_phys or kessler_phys
                                        !         or tj2016
logical, protected :: aqua_planet       ! Flag to run model in "aqua planet" mode
logical, protected :: moist_physics     ! true => moist physics enabled, i.e.,
                                        ! (.not. ideal_phys) .and. (.not. adiabatic)

logical, protected :: brnch_retain_casename ! true => branch run may use same caseid as
                                            !         the run being branched from

real(r8), protected :: eccen       ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)
real(r8), protected :: obliqr      ! Earth's obliquity in radians
real(r8), protected :: lambm0      ! Mean longitude of perihelion at the
                                   ! vernal equinox (radians)
real(r8), protected :: mvelpp      ! Earth's moving vernal equinox longitude
                                   ! of perihelion plus pi (radians)

!================================================================================================
contains
!================================================================================================

subroutine cam_ctrl_init( &
   caseid_in, ctitle_in, &
   initial_run_in, restart_run_in, branch_run_in, post_assim_in, &
   aqua_planet_in, brnch_retain_casename_in)

   character(len=cl), intent(in) :: caseid_in            ! case ID
   character(len=cl), intent(in) :: ctitle_in            ! case title
   logical,           intent(in) :: initial_run_in       ! true => inital run
   logical,           intent(in) :: restart_run_in       ! true => restart run
   logical,           intent(in) :: branch_run_in        ! true => branch run
   logical,           intent(in) :: post_assim_in        ! true => resume mode
   logical,           intent(in) :: aqua_planet_in       ! Flag to run model in "aqua planet" mode
   logical,           intent(in) :: brnch_retain_casename_in ! Flag to allow a branch to use the same
                                                             ! caseid as the run being branched from.

   integer :: unitn, ierr

   character(len=*), parameter :: sub='cam_ctrl_init'
   character(len=128) :: errmsg
   !---------------------------------------------------------------------------------------------

   caseid = caseid_in
   ctitle = ctitle_in

   initial_run = initial_run_in
   restart_run = restart_run_in
   branch_run  = branch_run_in
   post_assim  = post_assim_in

   aqua_planet = aqua_planet_in

   brnch_retain_casename = brnch_retain_casename_in

   if (masterproc) then
      write(iulog,*)' '
      write(iulog,*)' ------------------------------------------'
      write(iulog,*)' *********** CAM LOG OUTPUT ***************'
      write(iulog,*)' ------------------------------------------'
      if (restart_run) then
         write(iulog,*) '  Restart of an earlier run'
      else if (branch_run) then
         write(iulog,*) '  Branch of an earlier run'
      else if (post_assim) then
         write(iulog,*) '  DART run using CAM initial mode'
      else
         write(iulog,*) '         Initial run'
      end if
      write(iulog,*) ' ********** CASE = ',trim(caseid),' **********'
      write(iulog,'(1x,a)') ctitle


      if (aqua_planet) write(iulog,*) 'Run model in "AQUA_PLANET" mode'

   end if

end subroutine cam_ctrl_init

!--------------------------------------------------------------------------------------------------

subroutine cam_ctrl_set_orbit(eccen_in, obliqr_in, lambm0_in, mvelpp_in)

   real(r8), intent(in) :: eccen_in
   real(r8), intent(in) :: obliqr_in
   real(r8), intent(in) :: lambm0_in
   real(r8), intent(in) :: mvelpp_in

   eccen  = eccen_in
   obliqr = obliqr_in
   lambm0 = lambm0_in
   mvelpp = mvelpp_in

end subroutine cam_ctrl_set_orbit

!--------------------------------------------------------------------------------------------------

subroutine cam_ctrl_set_physics_type(phys_package)
  ! Dummy argument
  character(len=*), intent(in) :: phys_package
  ! Local variable
  character(len=*), parameter :: subname = 'cam_ctrl_set_physics_type'

  adiabatic = trim(phys_package) == 'adiabatic'
  ideal_phys = trim(phys_package) == 'held_suarez'
  kessler_phys = trim(phys_package) == 'kessler'
  tj2016_phys = trim(phys_package) == 'tj2016'

  simple_phys = adiabatic .or. ideal_phys .or. kessler_phys .or. tj2016_phys

  moist_physics = .not. (adiabatic .or. ideal_phys)

  if ((.not. moist_physics) .and. aqua_planet) then
    call endrun (subname//': FATAL: AQUA_PLANET not compatible with dry physics package, ('//trim(phys_package)//')')
  end if

  if (masterproc) then
    if (adiabatic) then
      write(iulog,*) 'Run model ADIABATICALLY (i.e. no physics)'
      write(iulog,*) '  Global energy fixer is on for non-Eulerian dycores.'
    else if (ideal_phys) then
      write(iulog,*) 'Run model with Held-Suarez physics forcing'
    else if (kessler_phys) then
      write(iulog,*) 'Run model with Kessler warm-rain physics forcing'
    else if (tj2016_phys) then
      write(iulog,*) 'Run model with Thatcher-Jablonowski (2016) physics forcing (moist Held-Suarez)'
    end if
  end if

end subroutine cam_ctrl_set_physics_type

end module cam_control_mod
