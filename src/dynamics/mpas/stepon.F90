#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module stepon

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: iam, masterproc, mpicom

use constituents,   only: pcnst, cnst_name, cnst_longname

use ppgrid,         only: begchunk, endchunk
use physics_types,  only: physics_state, physics_tend
use physics_buffer, only: physics_buffer_desc

use dyn_comp,       only: dyn_import_t, dyn_export_t, dyn_run, dyn_final

use dp_coupling,    only: d_p_coupling, p_d_coupling

use camsrfexch,     only: cam_out_t     

use time_manager,   only: get_step_size
use perf_mod,       only: t_startf, t_stopf, t_barrierf
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog
 
implicit none
private
save

public :: &
   stepon_init, &
   stepon_run1, &
   stepon_run2, &
   stepon_run3, &
   stepon_final


!=========================================================================================
contains
!=========================================================================================

subroutine stepon_init(dyn_in, dyn_out)

   use cam_history,    only: addfld, add_default, horiz_only
   use cam_history,    only: register_vector_field

   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

   integer :: m

   character(len=*), parameter :: subname = 'stepon::stepon_init'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   ! Forcing from physics
   call addfld('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term', gridname='physgrid')
   call addfld('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term', gridname='physgrid')
   call register_vector_field('FU', 'FV')

   call addfld('U&IC',   (/ 'lev' /),  'I', 'm/s', 'Zonal wind', gridname='physgrid' )
   call addfld('V&IC',   (/ 'lev' /),  'I', 'm/s', 'Meridional wind', gridname='physgrid' )
   call add_default ('U&IC       ',0, 'I')
   call add_default ('V&IC       ',0, 'I')

   call addfld('PS&IC', horiz_only,  'I', 'Pa', 'Surface pressure',gridname='physgrid')
   call addfld('T&IC',  (/ 'lev' /), 'I', 'K',  'Temperature', gridname='physgrid')
   call add_default ('PS&IC      ',0, 'I')
   call add_default ('T&IC       ',0, 'I')

   do m = 1,pcnst
      call addfld(trim(cnst_name(m))//'&IC', (/ 'lev' /), 'I', 'kg/kg', cnst_longname(m), gridname='physgrid')
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
   end do

end subroutine stepon_init


!-----------------------------------------------------------------------
!  routine stepon_run1
!
!> \brief First step in dynamics integration sequence
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine stepon_run1(dtime_out, phys_state, phys_tend, &
                       pbuf2d, dyn_in, dyn_out)

   real(r8),            intent(out)    :: dtime_out   ! Time-step
   type(physics_state), intent(inout)  :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout)  :: phys_tend(begchunk:endchunk)
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)
   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_run1'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   dtime_out = get_step_size()
   
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   ! Move data into phys_state structure.
   call d_p_coupling (phys_state, phys_tend, pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')
   
end subroutine stepon_run1


!-----------------------------------------------------------------------
!  routine stepon_run2
!
!> \brief Second step in dynamics integration sequence
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)

   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_run2'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)
 
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   ! copy from phys structures -> dynamics structures
   call p_d_coupling(phys_state, phys_tend, dyn_in)
   call t_stopf('p_d_coupling')

end subroutine stepon_run2


!-----------------------------------------------------------------------
!  routine stepon_run3
!
!> \brief Third step in dynamics integration sequence
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)

   ! arguments
   real(r8),            intent(in)    :: dtime
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_run3'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf('dyn_run')
   call dyn_run(dyn_in, dyn_out)	
   call t_stopf('dyn_run')

end subroutine stepon_run3


!-----------------------------------------------------------------------
!  routine stepon_final
!
!> \brief Finalizes dynamics timestepping
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine stepon_final(dyn_in, dyn_out)

   type(dyn_import_t), intent(inout) :: dyn_in
   type(dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_final'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   call dyn_final(dyn_in, dyn_out)

end subroutine stepon_final

end module stepon
