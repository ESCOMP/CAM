#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module stepon

use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: iam, masterproc, mpicom

use pmgrid,         only: plev, plevp
use constituents,   only: pcnst, cnst_name, cnst_longname

use ppgrid,         only: begchunk, endchunk
use physics_types,  only: physics_state, physics_tend
use physics_buffer, only: physics_buffer_desc

use dyn_comp,       only: dyn_import_t, dyn_export_t, dyn_run, dyn_final, &
                          swap_time_level_ptrs

use dp_coupling,    only: d_p_coupling, p_d_coupling

use camsrfexch,     only: cam_out_t     

use cam_history,    only: addfld, outfld

use time_manager,   only: get_step_size, get_nstep
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
   !----------------------------------------------------------------------------


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   ! dycore state variables on Spherical Centroidal Voronoi grid
   call addfld ('u',     (/ 'lev' /),  'A', 'm/s', 'normal velocity at edges', gridname='mpas_edge')
   call addfld ('w',     (/ 'ilev' /), 'A', 'm/s', 'vertical velocity', gridname='mpas_cell')
   call addfld ('theta', (/ 'lev' /),  'A', 'K',   'potential temperature', gridname='mpas_cell')
   call addfld ('rho',   (/ 'lev' /),  'A', 'kg/m^3', 'dry air density', gridname='mpas_cell')
   call addfld ('qv',    (/ 'lev' /),  'A', 'kg/kg', 'water vapor dry mmr', gridname='mpas_cell')
   call addfld ('uReconstructZonal', (/ 'lev' /),  'A', 'm/s', &
                'zonal velocity at cell centers', gridname='mpas_cell')
   call addfld ('uReconstructMeridional', (/ 'lev' /),  'A', 'm/s', &
                'meridional velocity at cell centers', gridname='mpas_cell')

   ! physics forcings on dycore grids
   call addfld ('ru_tend',     (/ 'lev' /),  'A', 'kg/m^2/s', &
                'physics tendency of normal horizontal momentum', gridname='mpas_edge')
   call addfld ('rtheta_tend', (/ 'lev' /),  'A', 'kg K/m^3/s', &
                'physics tendency of rho*theta/zz', gridname='mpas_cell')
   call addfld ('rho_tend', (/ 'lev' /),  'A', 'kg/m^3/s', &
                'physics tendency of dry air density', gridname='mpas_cell')

end subroutine stepon_init

!=========================================================================================

subroutine stepon_run1(dtime_out, phys_state, phys_tend, &
                       pbuf2d, dyn_in, dyn_out)

   real(r8),            intent(out)    :: dtime_out   ! Time-step
   type(physics_state), intent(inout)  :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout)  :: phys_tend(begchunk:endchunk)
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)
   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   ! local variables
   integer :: nstep

   character(len=*), parameter :: subname = 'stepon::stepon_run1'
   !----------------------------------------------------------------------------


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   nstep     = get_nstep()
   dtime_out = get_step_size()

   ! This call writes the inputs to the physics package on the dynamics
   ! grids.  It will output the initial fields the first time it is called
   ! since dyn_in is copied to dyn_out for the initialization sequence.  On
   ! subsequent calls dyn_out will contain the dycore output.
   call write_dynvar(dyn_out)
   
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   ! Move data into phys_state structure.
   call d_p_coupling (phys_state, phys_tend, pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')
   
   ! Update pointers for prognostic fields if necessary.  Note that this shift
   ! should not take place the first time stepon_run1 is called which is during
   ! the CAM initialization sequence before the dycore is called.
   if (nstep > 0 .and. swap_time_level_ptrs) call shift_time_levels(dyn_in, dyn_out)

end subroutine stepon_run1

!=========================================================================================

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)

   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_run2'
   !----------------------------------------------------------------------------


   MPAS_DEBUG_WRITE(1, 'begin '//subname)
 
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   ! copy from phys structures -> dynamics structures
   call p_d_coupling(phys_state, phys_tend, dyn_in)
   call t_stopf('p_d_coupling')

   ! write physics forcings input to dycore
   call write_forcings(dyn_in)

end subroutine stepon_run2

!=========================================================================================

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)

   ! arguments
   real(r8),            intent(in)    :: dtime
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_run3'
   !----------------------------------------------------------------------------


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf('dyn_run')
   call dyn_run(dyn_in, dyn_out)	
   call t_stopf('dyn_run')

end subroutine stepon_run3

!=========================================================================================

subroutine stepon_final(dyn_in, dyn_out)

   type(dyn_import_t), intent(inout) :: dyn_in
   type(dyn_export_t), intent(inout) :: dyn_out

   character(len=*), parameter :: subname = 'stepon::stepon_final'
   !----------------------------------------------------------------------------

   call dyn_final(dyn_in, dyn_out)

end subroutine stepon_final

!=========================================================================================
! Private
!=========================================================================================

subroutine write_dynvar(dyn_out)

   ! agruments
   type(dyn_export_t), intent(in) :: dyn_out

   ! local variables
   integer :: i, k, kk
   integer :: nCellsSolve, nEdgesSolve
   integer :: qv_idx
   real(r8), allocatable :: arr2d(:,:)

   character(len=*), parameter :: subname = 'stepon::write_dynvar'
   !----------------------------------------------------------------------------


   nCellsSolve = dyn_out%nCellsSolve
   nEdgesSolve = dyn_out%nEdgesSolve
   qv_idx      = dyn_out%index_qv

   allocate(arr2d(nEdgesSolve,plev))
   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nEdgesSolve
         arr2d(i,k) = dyn_out%uperp(kk,i)
      end do
   end do
   call outfld('u', arr2d, nEdgesSolve, 1)
   deallocate(arr2d)

   allocate(arr2d(nCellsSolve,plevp))
   do k = 1, plevp
      kk = plevp - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%w(kk,i)
      end do
   end do
   call outfld('w', arr2d, nCellsSolve, 1)
   deallocate(arr2d)

   allocate(arr2d(nCellsSolve,plev))

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%theta(kk,i)
      end do
   end do
   call outfld('theta', arr2d, nCellsSolve, 1)

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%rho(kk,i)
      end do
   end do
   call outfld('rho', arr2d, nCellsSolve, 1)

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%tracers(qv_idx,kk,i)
      end do
   end do
   call outfld('qv', arr2d, nCellsSolve, 1)

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%ux(kk,i)
      end do
   end do
   call outfld('uReconstructZonal', arr2d, nCellsSolve, 1)

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_out%uy(kk,i)
      end do
   end do
   call outfld('uReconstructMeridional', arr2d, nCellsSolve, 1)

   deallocate(arr2d)

end subroutine write_dynvar

!=========================================================================================

subroutine write_forcings(dyn_in)

   ! agruments
   type(dyn_import_t), intent(in) :: dyn_in

   ! local variables
   integer :: i, k, kk
   integer :: nCellsSolve, nEdgesSolve
   real(r8), allocatable :: arr2d(:,:)

   character(len=*), parameter :: subname = 'stepon::write_forcings'
   !----------------------------------------------------------------------------


   nCellsSolve = dyn_in%nCellsSolve
   nEdgesSolve = dyn_in%nEdgesSolve

   allocate(arr2d(nEdgesSolve,plev))
   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nEdgesSolve
         arr2d(i,k) = dyn_in%ru_tend(kk,i)
      end do
   end do
   call outfld('ru_tend', arr2d, nEdgesSolve, 1)
   deallocate(arr2d)

   allocate(arr2d(nCellsSolve,plev))

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_in%rtheta_tend(kk,i)
      end do
   end do
   call outfld('rtheta_tend', arr2d, nCellsSolve, 1)

   do k = 1, plev
      kk = plev - k + 1
      do i = 1, nCellsSolve
         arr2d(i,k) = dyn_in%rho_tend(kk,i)
      end do
   end do
   call outfld('rho_tend', arr2d, nCellsSolve, 1)

   deallocate(arr2d)

end subroutine write_forcings

!========================================================================================

subroutine shift_time_levels(dyn_in, dyn_out)

   ! The MPAS dycore swaps the pool time indices after each timestep
   ! (mpas_dt).  If an odd number of these shifts occur during the CAM
   ! timestep (i.e., the dynamics/physics coupling interval), then CAM
   ! needs a corresponding update to the pointers in the dyn_in and dyn_out
   ! objects.

   ! arguments
   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   ! local variables
   real(r8), dimension(:,:),   pointer :: ptr2d
   real(r8), dimension(:,:,:), pointer :: ptr3d
   !--------------------------------------------------------------------------------------

   ptr2d             => dyn_out % uperp
   dyn_out % uperp   => dyn_in % uperp
   dyn_in % uperp    => ptr2d

   ptr2d             => dyn_out % w
   dyn_out % w       => dyn_in % w
   dyn_in % w        => ptr2d

   ptr2d             => dyn_out % theta_m
   dyn_out % theta_m => dyn_in % theta_m
   dyn_in % theta_m  => ptr2d

   ptr2d             => dyn_out % rho_zz
   dyn_out % rho_zz  => dyn_in % rho_zz
   dyn_in % rho_zz   => ptr2d

   ptr3d             => dyn_out % tracers
   dyn_out % tracers => dyn_in % tracers
   dyn_in % tracers  => ptr3d

end subroutine shift_time_levels

end module stepon
