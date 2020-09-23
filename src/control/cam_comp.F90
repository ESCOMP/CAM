module cam_comp
!-----------------------------------------------------------------------
!
! Community Atmosphere Model (CAM) component interfaces.
!
! This interface layer is CAM specific, i.e., it deals entirely with CAM
! specific data structures.  It is the layer above this, either atm_comp_mct
! or atm_comp_esmf, which translates between CAM and either MCT or ESMF
! data structures in order to interface with the driver/coupler.
!
!-----------------------------------------------------------------------

use shr_kind_mod,      only: r8 => SHR_KIND_R8, cl=>SHR_KIND_CL, cs=>SHR_KIND_CS
use shr_sys_mod,       only: shr_sys_flush

use spmd_utils,        only: masterproc, mpicom
use cam_control_mod,   only: cam_ctrl_init, cam_ctrl_set_orbit
use runtime_opts,      only: read_namelist
use time_manager,      only: timemgr_init, get_step_size, &
                             get_nstep, is_first_step, is_first_restart_step

use camsrfexch,        only: cam_out_t, cam_in_t
use ppgrid,            only: begchunk, endchunk
use physics_types,     only: physics_state, physics_tend
use dyn_comp,          only: dyn_import_t, dyn_export_t

use physics_buffer,    only: physics_buffer_desc
use offline_driver,    only: offline_driver_init, offline_driver_dorun, offline_driver_run

use perf_mod
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun

implicit none
private

public cam_init      ! First phase of CAM initialization
public cam_run1      ! CAM run method phase 1
public cam_run2      ! CAM run method phase 2
public cam_run3      ! CAM run method phase 3
public cam_run4      ! CAM run method phase 4
public cam_final     ! CAM Finalization

type(dyn_import_t) :: dyn_in   ! Dynamics import container
type(dyn_export_t) :: dyn_out  ! Dynamics export container

type(physics_state),       pointer :: phys_state(:) => null()
type(physics_tend ),       pointer :: phys_tend(:) => null()
type(physics_buffer_desc), pointer :: pbuf2d(:,:) => null()

real(r8) :: dtime_phys         ! Time step for physics tendencies.  Set by call to
                               ! stepon_run1, then passed to the phys_run*

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

subroutine cam_init(                                             &
   caseid, ctitle, model_doi_url,                                &
   initial_run_in, restart_run_in, branch_run_in, post_assim_in, &
   calendar, brnch_retain_casename, aqua_planet,                 &
   single_column, scmlat, scmlon,                                &
   eccen, obliqr, lambm0, mvelpp,                                &
   perpetual_run, perpetual_ymd,                                 &
   dtime, start_ymd, start_tod, ref_ymd, ref_tod,                &
   stop_ymd, stop_tod, curr_ymd, curr_tod,                       &
   cam_out, cam_in)

   !-----------------------------------------------------------------------
   !
   ! CAM component initialization.
   !
   !-----------------------------------------------------------------------

   use history_defaults, only: bldfld
   use cam_initfiles,    only: cam_initfiles_open
   use dyn_grid,         only: dyn_grid_init
   use phys_grid,        only: phys_grid_init
   use physpkg,          only: phys_register, phys_init
   use chem_surfvals,    only: chem_surfvals_init
   use dyn_comp,         only: dyn_init
   use cam_restart,      only: cam_read_restart
   use stepon,           only: stepon_init
   use ionosphere_interface, only: ionosphere_init
   use camsrfexch,       only: hub2atm_alloc, atm2hub_alloc
   use cam_history,      only: intht
   use history_scam,     only: scm_intht
   use cam_pio_utils,    only: init_pio_subsystem
   use cam_instance,     only: inst_suffix
   use cam_snapshot,     only: cam_snapshot_deactivate
   use physconst,        only: composition_init
#if (defined BFB_CAM_SCAM_IOP)
   use history_defaults, only: initialize_iop_history
#endif


   ! Arguments
   character(len=cl), intent(in) :: caseid                ! case ID
   character(len=cl), intent(in) :: ctitle                ! case title
   character(len=cl), intent(in) :: model_doi_url         ! CESM model DOI
   logical,           intent(in) :: initial_run_in        ! true => inital run
   logical,           intent(in) :: restart_run_in        ! true => restart run
   logical,           intent(in) :: branch_run_in         ! true => branch run
   logical,           intent(in) :: post_assim_in         ! true => resume mode
   character(len=cs), intent(in) :: calendar              ! Calendar type
   logical,           intent(in) :: brnch_retain_casename ! Flag to allow a branch to use the same
                                                          ! caseid as the run being branched from.
   logical,           intent(in) :: aqua_planet           ! Flag to run model in "aqua planet" mode

   logical,           intent(in) :: single_column
   real(r8),          intent(in) :: scmlat
   real(r8),          intent(in) :: scmlon

   real(r8),          intent(in) :: eccen
   real(r8),          intent(in) :: obliqr
   real(r8),          intent(in) :: lambm0
   real(r8),          intent(in) :: mvelpp

   logical,           intent(in) :: perpetual_run    ! true => perpetual mode enabled
   integer,           intent(in) :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
   integer,           intent(in) :: dtime                 ! model timestep (sec)

   integer,           intent(in) :: start_ymd             ! Start date (YYYYMMDD)
   integer,           intent(in) :: start_tod             ! Start time of day (sec)
   integer,           intent(in) :: curr_ymd              ! Start date (YYYYMMDD)
   integer,           intent(in) :: curr_tod              ! Start time of day (sec)
   integer,           intent(in) :: stop_ymd              ! Stop date (YYYYMMDD)
   integer,           intent(in) :: stop_tod              ! Stop time of day (sec)
   integer,           intent(in) :: ref_ymd               ! Reference date (YYYYMMDD)
   integer,           intent(in) :: ref_tod               ! Reference time of day (sec)

   type(cam_out_t),   pointer    :: cam_out(:)       ! Output from CAM to surface
   type(cam_in_t) ,   pointer    :: cam_in(:)        ! Merged input state to CAM

   ! Local variables
   character(len=cs) :: filein      ! Input namelist filename
   !-----------------------------------------------------------------------

   call init_pio_subsystem()

   ! Initializations using data passed from coupler.
   call cam_ctrl_init(               &
      caseid_in=caseid,              &
      ctitle_in=ctitle,              &
      initial_run_in=initial_run_in, &
      restart_run_in=restart_run_in, &
      branch_run_in=branch_run_in,   &
      post_assim_in=post_assim_in,   &
      aqua_planet_in=aqua_planet,    &
      brnch_retain_casename_in=brnch_retain_casename)

   call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)

   call timemgr_init( &
      dtime, calendar, start_ymd, start_tod, ref_ymd,  &
      ref_tod, stop_ymd, stop_tod, curr_ymd, curr_tod, &
      perpetual_run, perpetual_ymd, initial_run_in)

   ! Read CAM namelists.
   filein = "atm_in" // trim(inst_suffix)
   call read_namelist(filein, single_column, scmlat, scmlon)

   ! Open initial or restart file, and topo file if specified.
   call cam_initfiles_open()

   ! Initialize grids and dynamics grid decomposition
   call dyn_grid_init()

   ! Initialize physics grid decomposition
   call phys_grid_init()

   ! Register advected tracers and physics buffer fields
   call phys_register ()

   ! Initialize ghg surface values before default initial distributions
   ! are set in dyn_init
   call chem_surfvals_init()

   call composition_init()
   ! initialize ionosphere
   call ionosphere_init()

   if (initial_run_in) then

      call dyn_init(dyn_in, dyn_out)

      ! Allocate and setup surface exchange data
      call atm2hub_alloc(cam_out)
      call hub2atm_alloc(cam_in)

   else

      call cam_read_restart(cam_in, cam_out, dyn_in, dyn_out, pbuf2d, stop_ymd, stop_tod)

#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
   end if

   call phys_init( phys_state, phys_tend, pbuf2d, cam_in, cam_out )

   call bldfld ()       ! master field list (if branch, only does hash tables)

   call stepon_init(dyn_in, dyn_out)

   call offline_driver_init()

   if (single_column) call scm_intht()
   call intht(model_doi_url)

   call cam_snapshot_deactivate()

end subroutine cam_init

!
!-----------------------------------------------------------------------
!
subroutine cam_run1(cam_in, cam_out)
!-----------------------------------------------------------------------
!
! Purpose:   First phase of atmosphere model run method.
!            Runs first phase of dynamics and first phase of
!            physics (before surface model updates).
!
!-----------------------------------------------------------------------

   use physpkg,          only: phys_run1
   use stepon,           only: stepon_run1
   use ionosphere_interface,only: ionosphere_run1

   type(cam_in_t)  :: cam_in(begchunk:endchunk)
   type(cam_out_t) :: cam_out(begchunk:endchunk)

   !-----------------------------------------------------------------------

   if (offline_driver_dorun) return

   !----------------------------------------------------------
   ! First phase of dynamics (at least couple from dynamics to physics)
   ! Return time-step for physics from dynamics.
   !----------------------------------------------------------
   call t_barrierf ('sync_stepon_run1', mpicom)
   call t_startf ('stepon_run1')
   call stepon_run1( dtime_phys, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out )
   call t_stopf  ('stepon_run1')

   !----------------------------------------------------------
   ! first phase of ionosphere -- write to IC file if needed
   !----------------------------------------------------------
   call ionosphere_run1(pbuf2d)

   !
   !----------------------------------------------------------
   ! PHYS_RUN Call the Physics package
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_phys_run1', mpicom)
   call t_startf ('phys_run1')
   call phys_run1(phys_state, dtime_phys, phys_tend, pbuf2d,  cam_in, cam_out)
   call t_stopf  ('phys_run1')

end subroutine cam_run1

!
!-----------------------------------------------------------------------
!

subroutine cam_run2( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:   Second phase of atmosphere model run method.
!            Run the second phase physics, run methods that
!            require the surface model updates.  And run the
!            second phase of dynamics that at least couples
!            between physics to dynamics.
!
!-----------------------------------------------------------------------

   use physpkg,          only: phys_run2
   use stepon,           only: stepon_run2
   use ionosphere_interface, only: ionosphere_run2

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
   type(cam_in_t),  intent(inout) :: cam_in(begchunk:endchunk)

   if (offline_driver_dorun) then
      call offline_driver_run( phys_state, pbuf2d, cam_out, cam_in )
      return
   endif

   !
   ! Second phase of physics (after surface model update)
   !
   call t_barrierf ('sync_phys_run2', mpicom)
   call t_startf ('phys_run2')
   call phys_run2(phys_state, dtime_phys, phys_tend, pbuf2d,  cam_out, cam_in )
   call t_stopf  ('phys_run2')

   !
   ! Second phase of dynamics (at least couple from physics to dynamics)
   !
   call t_barrierf ('sync_stepon_run2', mpicom)
   call t_startf ('stepon_run2')
   call stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
   call t_stopf  ('stepon_run2')

   !
   ! Ion transport
   !
   call t_startf('ionosphere_run2')
   call ionosphere_run2( phys_state, dyn_in, pbuf2d )
   call t_stopf ('ionosphere_run2')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run2_memusage')
      call t_stopf  ('cam_run2_memusage')
   end if
end subroutine cam_run2

!
!-----------------------------------------------------------------------
!

subroutine cam_run3( cam_out )
!-----------------------------------------------------------------------
!
! Purpose:  Third phase of atmosphere model run method. This consists
!           of the third phase of the dynamics. For some dycores
!           this will be the actual dynamics run, for others the
!           dynamics happens before physics in phase 1.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_run3

   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!-----------------------------------------------------------------------

   if (offline_driver_dorun) return

   !
   ! Third phase of dynamics
   !
   call t_barrierf ('sync_stepon_run3', mpicom)
   call t_startf ('stepon_run3')
   call stepon_run3( dtime_phys, cam_out, phys_state, dyn_in, dyn_out )

   call t_stopf  ('stepon_run3')

   if (is_first_step() .or. is_first_restart_step()) then
      call t_startf ('cam_run3_memusage')
      call t_stopf  ('cam_run3_memusage')
   end if
end subroutine cam_run3

!
!-----------------------------------------------------------------------
!

subroutine cam_run4( cam_out, cam_in, rstwr, nlend, &
                     yr_spec, mon_spec, day_spec, sec_spec )

!-----------------------------------------------------------------------
!
! Purpose:  Final phase of atmosphere model run method. This consists
!           of all the restart output, history writes, and other
!           file output.
!
!-----------------------------------------------------------------------
   use cam_history,      only: wshist, wrapup
   use cam_restart,      only: cam_write_restart
   use qneg_module,      only: qneg_print_summary
   use time_manager,     only: is_last_step

   type(cam_out_t), intent(inout)        :: cam_out(begchunk:endchunk)
   type(cam_in_t) , intent(inout)        :: cam_in(begchunk:endchunk)
   logical            , intent(in)           :: rstwr           ! true => write restart file
   logical            , intent(in)           :: nlend           ! true => this is final timestep
   integer            , intent(in), optional :: yr_spec         ! Simulation year
   integer            , intent(in), optional :: mon_spec        ! Simulation month
   integer            , intent(in), optional :: day_spec        ! Simulation day
   integer            , intent(in), optional :: sec_spec        ! Seconds into current simulation day

   !----------------------------------------------------------
   ! History and restart logic: Write and/or dispose history tapes if required
   !----------------------------------------------------------
   !
   call t_barrierf ('sync_wshist', mpicom)
   call t_startf ('wshist')
   call wshist ()
   call t_stopf  ('wshist')

   !
   ! Write restart files
   !
   if (rstwr) then
      call t_startf ('cam_write_restart')
      if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
         call cam_write_restart(cam_in, cam_out, dyn_out, pbuf2d, &
              yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
      else
         call cam_write_restart(cam_in, cam_out, dyn_out, pbuf2d )
      end if
      call t_stopf  ('cam_write_restart')
   end if

   call t_startf ('cam_run4_wrapup')
   call wrapup(rstwr, nlend)
   call t_stopf  ('cam_run4_wrapup')

   call qneg_print_summary(is_last_step())

   call shr_sys_flush(iulog)

end subroutine cam_run4

!
!-----------------------------------------------------------------------
!

subroutine cam_final( cam_out, cam_in )
!-----------------------------------------------------------------------
!
! Purpose:  CAM finalization.
!
!-----------------------------------------------------------------------
   use stepon,           only: stepon_final
   use physpkg,          only: phys_final
   use cam_initfiles,    only: cam_initfiles_close
   use camsrfexch,       only: atm2hub_deallocate, hub2atm_deallocate
   use ionosphere_interface, only: ionosphere_final
   use cam_control_mod,  only: initial_run

   !
   ! Arguments
   !
   type(cam_out_t), pointer :: cam_out(:) ! Output from CAM to surface
   type(cam_in_t),  pointer :: cam_in(:)   ! Input from merged surface to CAM

   ! Local variables
   integer :: nstep           ! Current timestep number.
   !-----------------------------------------------------------------------

   call phys_final( phys_state, phys_tend , pbuf2d)
   call stepon_final(dyn_in, dyn_out)
   call ionosphere_final()

   if (initial_run) then
      call cam_initfiles_close()
   end if

   call hub2atm_deallocate(cam_in)
   call atm2hub_deallocate(cam_out)

   ! This flush attempts to ensure that asynchronous diagnostic prints from all
   ! processes do not get mixed up with the "END OF MODEL RUN" message printed
   ! by masterproc below.  The test-model script searches for this message in the
   ! output log to figure out if CAM completed successfully.
   call shr_sys_flush( 0 )       ! Flush all output to standard error
   call shr_sys_flush( iulog )   ! Flush all output to the CAM log file

   if (masterproc) then
      nstep = get_nstep()
      write(iulog,9300) nstep-1,nstep
9300  format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
                ' partially done to provide convectively adjusted and ', &
                'time filtered values for history tape.')
      write(iulog,*)' '
      write(iulog,*)'******* END OF MODEL RUN *******'
   end if

end subroutine cam_final

!-----------------------------------------------------------------------

end module cam_comp
