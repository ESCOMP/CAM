module stepon

use cam_abortutils, only: endrun
use shr_kind_mod,   only: r8 => shr_kind_r8
use spmd_utils,     only: mpicom

use pmgrid,         only: plev, plevp

use ppgrid,         only: begchunk, endchunk
use physics_types,  only: physics_state, physics_tend
use physics_buffer, only: physics_buffer_desc

use dyn_comp,       only: dyn_import_t, dyn_export_t, dyn_run, dyn_final

use dp_coupling,    only: d_p_coupling, p_d_coupling

use camsrfexch,     only: cam_out_t

use cam_history,    only: addfld, outfld, hist_fld_active, write_inithist

use time_manager,   only: get_step_size, get_nstep, is_first_step, is_first_restart_step

use perf_mod,       only: t_startf, t_stopf, t_barrierf

use aerosol_properties_mod, only: aerosol_properties
use aerosol_state_mod,      only: aerosol_state
use microp_aero,            only: aerosol_state_object, aerosol_properties_object

implicit none
private
save

public :: &
   stepon_init, &
   stepon_run1, &
   stepon_run2, &
   stepon_run3, &
   stepon_final

class(aerosol_properties), pointer :: aero_props_obj => null()
logical :: aerosols_transported = .false.

!=========================================================================================
contains
!=========================================================================================

subroutine stepon_init(dyn_in, dyn_out)

   ! arguments
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   !----------------------------------------------------------------------------

   ! dycore state variables on MPAS grids
   call addfld ('u',     (/ 'lev' /),  'A', 'm/s', 'normal velocity at edges', gridname='mpas_edge')
   call addfld ('w',     (/ 'ilev' /), 'A', 'm/s', 'vertical velocity', gridname='mpas_cell')
   call addfld ('theta', (/ 'lev' /),  'A', 'K',   'potential temperature', gridname='mpas_cell')
   call addfld ('rho',   (/ 'lev' /),  'A', 'kg/m^3', 'dry air density', gridname='mpas_cell')
   call addfld ('qv',    (/ 'lev' /),  'A', 'kg/kg', 'water vapor dry mmr', gridname='mpas_cell')
   call addfld ('uReconstructZonal', (/ 'lev' /),  'A', 'm/s', &
                'zonal velocity at cell centers', gridname='mpas_cell')
   call addfld ('uReconstructMeridional', (/ 'lev' /),  'A', 'm/s', &
                'meridional velocity at cell centers', gridname='mpas_cell')
   call addfld ('divergence', (/ 'lev' /), 'A', '1/s', &
                'Horizontal velocity divergence at cell center', gridname='mpas_cell')
   call addfld ('vorticity', (/ 'lev' /), 'A', '1/s', &
                'Relative vorticity at vertices', gridname='mpas_vertex')

   ! physics forcings on MPAS grids
   call addfld ('ru_tend',     (/ 'lev' /),  'A', 'kg/m^2/s', &
                'physics tendency of normal horizontal momentum', gridname='mpas_edge')
   call addfld ('rtheta_tend', (/ 'lev' /),  'A', 'kg K/m^3/s', &
                'physics tendency of rho*theta/zz', gridname='mpas_cell')
   call addfld ('rho_tend', (/ 'lev' /),  'A', 'kg/m^3/s', &
                'physics tendency of dry air density', gridname='mpas_cell')

   ! get aerosol properties
   aero_props_obj => aerosol_properties_object()

   if (associated(aero_props_obj)) then
      ! determine if there are transported aerosol contistuents
      aerosols_transported = aero_props_obj%number_transported()>0
   end if

end subroutine stepon_init

!=========================================================================================

subroutine stepon_run1(dtime_out, phys_state, phys_tend, &
                       pbuf2d, dyn_in, dyn_out)

   ! arguments
   real(r8),            intent(out)    :: dtime_out   ! Time-step
   type(physics_state), intent(inout)  :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout)  :: phys_tend(begchunk:endchunk)
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)
   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   ! local variables
   integer :: nstep

   integer :: c
   class(aerosol_state), pointer :: aero_state_obj
   nullify(aero_state_obj)

   !----------------------------------------------------------------------------

   nstep     = get_nstep()
   dtime_out = get_step_size()

   ! This call writes the dycore output (on the dynamics grids) to the
   ! history file.  Note that when nstep=0, these fields will be the result
   ! of the dynamics initialization (done in dyn_init) since the dycore
   ! does not run and dyn_in points to the same memory as dyn_out.
   call write_dynvar(dyn_out)

   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   ! Move data into phys_state structure.
   call d_p_coupling (phys_state, phys_tend, pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')

   !----------------------------------------------------------
   ! update aerosol state object from CAM physics state constituents
   !----------------------------------------------------------
   if (aerosols_transported) then

      do c = begchunk,endchunk
         aero_state_obj => aerosol_state_object(c)
         ! pass number mass or number mixing ratios of aerosol constituents
         ! to aerosol state object
         call aero_state_obj%set_transported(phys_state(c)%q)
      end do

   end if

end subroutine stepon_run1

!=========================================================================================

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)

   ! arguments
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out
   !----------------------------------------------------------------------------

   integer :: c
   class(aerosol_state), pointer :: aero_state_obj

   !----------------------------------------------------------
   ! update physics state with aerosol constituents
   !----------------------------------------------------------
   nullify(aero_state_obj)

   if (aerosols_transported) then
      do c = begchunk,endchunk
         aero_state_obj => aerosol_state_object(c)
         ! get mass or number mixing ratios of aerosol constituents
         call aero_state_obj%get_transported(phys_state(c)%q)
      end do
   end if

   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   ! copy from phys structures -> dynamics structures
   call p_d_coupling(phys_state, phys_tend, dyn_in)
   call t_stopf('p_d_coupling')

   ! write physics forcings which are input to dycore
   call write_forcings(dyn_in)

end subroutine stepon_run2

!=========================================================================================

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)

   ! arguments
   real(r8),            intent(in)    :: dtime
   type(cam_out_t),     intent(inout) :: cam_out(:) ! CAM export to surface models
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in
   type (dyn_export_t), intent(inout) :: dyn_out
   !----------------------------------------------------------------------------

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf('dyn_run')
   call dyn_run(dyn_in, dyn_out)
   call t_stopf('dyn_run')

end subroutine stepon_run3

!=========================================================================================

subroutine stepon_final(dyn_in, dyn_out)

   ! arguments
   type(dyn_import_t), intent(inout) :: dyn_in
   type(dyn_export_t), intent(inout) :: dyn_out
   !----------------------------------------------------------------------------

   call dyn_final(dyn_in, dyn_out)

end subroutine stepon_final

!=========================================================================================
! Private
!=========================================================================================

subroutine write_dynvar(dyn_out)

   ! Output from the internal MPAS data structures to CAM history files.
   ! Make call to MPAS to write an initial file when requested.

   use string_utils, only: int2str

   ! agruments
   type(dyn_export_t), intent(in) :: dyn_out

   ! local variables
   integer :: i, k, kk
   integer :: nCellsSolve, nEdgesSolve, nVerticesSolve
   integer :: qv_idx
   real(r8), allocatable :: arr2d(:,:)
   integer :: ierr
   character(len=*), parameter :: subname = 'stepon::write_dynvar'
   !----------------------------------------------------------------------------

   nCellsSolve    = dyn_out%nCellsSolve
   nEdgesSolve    = dyn_out%nEdgesSolve
   nVerticesSolve = dyn_out%nVerticesSolve
   qv_idx         = dyn_out%index_qv

   if (hist_fld_active('u')) then
      allocate(arr2d(nEdgesSolve,plev), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nEdgesSolve
            arr2d(i,k) = dyn_out%uperp(kk,i)
         end do
      end do
      call outfld('u', arr2d, nEdgesSolve, 1)
      deallocate(arr2d)
   end if

   if (hist_fld_active('w')) then
      allocate(arr2d(nCellsSolve,plevp), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))
      do k = 1, plevp
         kk = plevp - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%w(kk,i)
         end do
      end do
      call outfld('w', arr2d, nCellsSolve, 1)
      deallocate(arr2d)
   end if

   allocate(arr2d(nCellsSolve,plev), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))

   if (hist_fld_active('theta')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%theta(kk,i)
         end do
      end do
      call outfld('theta', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('rho')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%rho(kk,i)
         end do
      end do
      call outfld('rho', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('qv')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%tracers(qv_idx,kk,i)
         end do
      end do
      call outfld('qv', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('uReconstructZonal')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%ux(kk,i)
         end do
      end do
      call outfld('uReconstructZonal', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('uReconstructMeridional')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%uy(kk,i)
         end do
      end do
      call outfld('uReconstructMeridional', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('divergence')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_out%divergence(kk,i)
         end do
      end do
      call outfld('divergence', arr2d, nCellsSolve, 1)
   end if

   deallocate(arr2d)

   if (hist_fld_active('vorticity')) then
      allocate(arr2d(nVerticesSolve,plev), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nVerticesSolve
            arr2d(i,k) = dyn_out%vorticity(kk,i)
         end do
      end do
      call outfld('vorticity', arr2d, nVerticesSolve, 1)
      deallocate(arr2d)
   end if

   if (write_inithist()) then
      call write_initial_file()
   end if

end subroutine write_dynvar

!=========================================================================================

subroutine write_forcings(dyn_in)

   ! Output from the internal MPAS data structures to CAM history files.

   use string_utils, only: int2str

   ! agruments
   type(dyn_import_t), intent(in) :: dyn_in

   ! local variables
   integer :: i, k, kk
   integer :: nCellsSolve, nEdgesSolve
   real(r8), allocatable :: arr2d(:,:)
   integer :: ierr
   character(len=*), parameter :: subname = 'dyn_grid::write_forcings'

   !----------------------------------------------------------------------------


   nCellsSolve = dyn_in%nCellsSolve
   nEdgesSolve = dyn_in%nEdgesSolve

   if (hist_fld_active('ru_tend')) then
      allocate(arr2d(nEdgesSolve,plev), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nEdgesSolve
            arr2d(i,k) = dyn_in%ru_tend(kk,i)
         end do
      end do
      call outfld('ru_tend', arr2d, nEdgesSolve, 1)
      deallocate(arr2d)
   end if

   allocate(arr2d(nCellsSolve,plev), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate arr2d array at line:'//int2str(__LINE__))

   if (hist_fld_active('rtheta_tend')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_in%rtheta_tend(kk,i)
         end do
      end do
      call outfld('rtheta_tend', arr2d, nCellsSolve, 1)
   end if

   if (hist_fld_active('rho_tend')) then
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            arr2d(i,k) = dyn_in%rho_tend(kk,i)
         end do
      end do
      call outfld('rho_tend', arr2d, nCellsSolve, 1)
   end if

   deallocate(arr2d)

end subroutine write_forcings

!========================================================================================

subroutine write_initial_file()

   ! Make use of the MPAS functionality for writting a restart file to
   ! write an initial file.

   use shr_kind_mod,     only: cl=>shr_kind_cl
   use cam_instance,     only: inst_suffix
   use time_manager,     only: get_curr_date, get_stop_date, timemgr_datediff
   use filenames,        only: interpret_filename_spec
   use pio,              only: file_desc_t,  pio_enddef, &
                               pio_seterrorhandling, PIO_BCAST_ERROR
   use cam_pio_utils,    only: cam_pio_createfile, cam_pio_closefile
   use cam_abortutils,   only: endrun

   use mpas_derived_types, only: MPAS_Stream_type, MPAS_IO_WRITE
   use cam_mpas_subdriver, only: cam_mpas_setup_restart, cam_mpas_write_restart

   ! Local variables
   integer           :: yr, mon, day, tod1, tod2, ymd1, ymd2
   real(r8)          :: days
   character(len=cl) :: filename_spec ! filename specifier
   character(len=cl) :: fname         ! initial filename
   type(file_desc_t) :: fh
   integer           :: ierr

   type (MPAS_Stream_type) :: initial_stream
   !----------------------------------------------------------------------------

   ! Check whether the current time is during the final partial timestep taken by
   ! CAM.  Don't write the initial file during that time.  This avoids the problem
   ! of having an initial file written with a timestamp that is after the stop date.
   call get_curr_date(yr, mon, day, tod1)
   ymd1 = 10000*yr + 100*mon + day
   call get_stop_date(yr, mon, day, tod2)
   ymd2 = 10000*yr + 100*mon + day
   ! (ymd2,tod2) - (ymd1,tod1)
   call timemgr_datediff(ymd1, tod1, ymd2, tod2, days)
   if (days < 0._r8) return

   ! Set filename template for initial file based on instance suffix
   ! (%c = caseid, %y = year, %m = month, %d = day, %s = seconds in day)
   filename_spec = '%c.cam' // trim(inst_suffix) //'.i.%y-%m-%d-%s.nc'

   fname = interpret_filename_spec( filename_spec )

   call cam_pio_createfile(fh, trim(fname), 0)

   call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

   call cam_mpas_setup_restart(fh, initial_stream, MPAS_IO_WRITE, endrun)

   ierr = pio_enddef(fh)

   call cam_mpas_write_restart(initial_stream, endrun)

   call cam_pio_closefile(fh)

end subroutine write_initial_file

!========================================================================================

end module stepon
