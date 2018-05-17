module cam_restart

! Coordinate reading and writing of restart files.

use shr_kind_mod,     only: cl=>shr_kind_cl
use spmd_utils,       only: masterproc
use cam_control_mod,  only: restart_run, caseid
use ioFileMod,        only: opnfil
use camsrfexch,       only: cam_in_t, cam_out_t     
use dyn_comp,         only: dyn_import_t, dyn_export_t
use physics_buffer,   only: physics_buffer_desc
use units,            only: getunit, freeunit
use pio,              only: file_desc_t, pio_global, pio_enddef, &
                            pio_put_att, pio_closefile

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use ionosphere_interface, only: ionosphere_init_restart, ionosphere_write_restart, ionosphere_read_restart

implicit none
private
save

public :: &
   cam_write_restart,  &  ! Driver for writing restart files
   cam_read_restart       ! Driver for reading restart files

!=========================================================================================
contains
!=========================================================================================

subroutine cam_read_restart(cam_in, cam_out, dyn_in, dyn_out, pbuf2d, &
                            stop_ymd, stop_tod)

   use cam_initfiles,    only: initial_file_get_id
   use restart_dynamics, only: read_restart_dynamics
   use restart_physics,  only: read_restart_physics
   use camsrfexch,       only: atm2hub_alloc, hub2atm_alloc
   use cam_history,      only: read_restart_history
   use cam_pio_utils,    only: clean_iodesc_list

   ! Arguments
   type(cam_in_t),            pointer       :: cam_in(:)
   type(cam_out_t),           pointer       :: cam_out(:)
   type(dyn_import_t),        intent(inout) :: dyn_in
   type(dyn_export_t),        intent(inout) :: dyn_out
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)
   integer,                   intent(in)    :: stop_ymd       ! Stop date (YYYYMMDD)
   integer,                   intent(in)    :: stop_tod       ! Stop time of day (sec)

   ! Local workspace
   type(file_desc_t), pointer :: fh_ini

   character(len=*), parameter :: sub = 'cam_read_restart'
   !---------------------------------------------------------------------------
  
   ! get filehandle pointer to primary restart file
   fh_ini => initial_file_get_id()

   call read_restart_dynamics(fh_ini, dyn_in, dyn_out)   
   call ionosphere_read_restart(fh_ini)

   call hub2atm_alloc(cam_in)
   call atm2hub_alloc(cam_out)

   call read_restart_physics(fh_ini, cam_in, cam_out, pbuf2d)

   if (restart_run) then
      call read_restart_history (fh_ini)
   end if

   call clean_iodesc_list()

end subroutine cam_read_restart

!=========================================================================================

subroutine cam_write_restart(cam_in, cam_out, dyn_out, pbuf2d, &
                             yr_spec, mon_spec, day_spec, sec_spec )

   use filenames,        only: interpret_filename_spec
   use cam_pio_utils,    only: cam_pio_createfile
   use restart_dynamics, only: write_restart_dynamics, init_restart_dynamics
   use restart_physics,  only: write_restart_physics, init_restart_physics
   use cam_history,      only: write_restart_history, init_restart_history
   use cam_instance,     only: inst_suffix

   ! Arguments
   type(cam_in_t),          intent(in) :: cam_in(:)
   type(cam_out_t),         intent(in) :: cam_out(:)
   type(dyn_export_t),      intent(in) :: dyn_out
   type(physics_buffer_desc), pointer  :: pbuf2d(:,:)
   integer,       optional, intent(in) :: yr_spec         ! Simulation year
   integer,       optional, intent(in) :: mon_spec        ! Simulation month
   integer,       optional, intent(in) :: day_spec        ! Simulation day
   integer,       optional, intent(in) :: sec_spec        ! Seconds into current simulation day

   ! Local workspace
   character(len=cl) :: rfilename_spec ! filename specifier for primary restart file
   character(len=cl) :: fname  ! Restart filename
   type(file_desc_t) :: fh
   integer           :: ierr
   !-----------------------------------------------------------------------

   ! Set template for primary restart filename based on instance suffix
   ! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = number)
   rfilename_spec = '%c.cam' // trim(inst_suffix) //'.r.%y-%m-%d-%s.nc'

   if (present(yr_spec).and.present(mon_spec).and.present(day_spec).and.present(sec_spec)) then
      fname = interpret_filename_spec( rfilename_spec, &
              yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
   else
      fname = interpret_filename_spec( rfilename_spec )
   end if

   call cam_pio_createfile(fh, trim(fname), 0)

   call init_restart_dynamics(fh, dyn_out)
   call ionosphere_init_restart(fh)
   call init_restart_physics(fh, pbuf2d)
   call init_restart_history(fh)

   ierr = pio_put_att(fh, pio_global, 'caseid', caseid)

   ierr = pio_enddef(fh)

   !-----------------------------------------------------------------------
   ! Dynamics, physics, History
   !-----------------------------------------------------------------------

   call write_restart_dynamics(fh, dyn_out)
   call ionosphere_write_restart(fh)
   call write_restart_physics(fh, cam_in, cam_out, pbuf2d)

   if (present(yr_spec).and.present(mon_spec).and.&
      present(day_spec).and.present(sec_spec)) then
      call write_restart_history(fh, yr_spec=yr_spec, mon_spec=mon_spec, &
                                 day_spec=day_spec, sec_spec= sec_spec )
   else
      call write_restart_history(fh)
   end if

   ! Close the primary restart file
   call pio_closefile(fh)
      
   ! Update the restart pointer file
   call write_rest_pfile(fname)

end subroutine cam_write_restart

!========================================================================================

subroutine write_rest_pfile(restart_file)

   ! Write the restart pointer file

   use cam_initfiles, only: rest_pfile

   character(len=*), intent(in) :: restart_file

   integer :: nsds, ierr
   character(len=*), parameter :: sub='write_rest_pfile'
   !---------------------------------------------------------------------------
   
   if (masterproc) then

      nsds = getunit()
      call opnfil(rest_pfile, nsds, 'f')
      rewind nsds
      write(nsds, '(a)', iostat=ierr) trim(restart_file)
      if (ierr /= 0) then
         call endrun(sub//': ERROR: writing rpointer file')
      end if
      close(nsds)
      call freeunit(nsds)

      write(iulog,*)'(WRITE_REST_PFILE): successfully wrote local restart pointer file ',&
         trim(rest_pfile)
      write(iulog,'("---------------------------------------")')
   end if

end subroutine write_rest_pfile

!========================================================================================

end module cam_restart
