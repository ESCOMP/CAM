module restart_dynamics

! Writing and reading grid and dynamics state information to/from restart files is
! delegated to MPAS utility code.  This module provides the CAM interfaces for the
! restart functionality.  CAM just provides MPAS with the PIO filehandle to the
! restart file.

use dyn_comp,           only: dyn_import_t, dyn_export_t, dyn_init
use pio,                only: file_desc_t

use cam_abortutils,     only: endrun

use mpas_derived_types, only: MPAS_Stream_type, MPAS_IO_WRITE, MPAS_IO_READ
use cam_mpas_subdriver, only: domain_ptr, cam_mpas_setup_restart, cam_mpas_write_restart, &
                              cam_mpas_read_restart, cam_mpas_define_scalars

implicit none
private
save

public :: &
   init_restart_dynamics,  &
   write_restart_dynamics, &
   read_restart_dynamics

! The restart_stream is set up in init_restart_dynamics and used later in
! write_restart_dynamics and read_restart_dynamics.
type (MPAS_Stream_type) :: restart_stream

!=========================================================================================
contains
!=========================================================================================

subroutine init_restart_dynamics(file, dyn_out)

   ! arguments
   type(file_desc_t),  target     :: File
   type(dyn_export_t), intent(in) :: dyn_out
   !----------------------------------------------------------------------------

   call cam_mpas_setup_restart(file, restart_stream, MPAS_IO_WRITE, endrun)

end subroutine init_restart_dynamics

!=========================================================================================

subroutine write_restart_dynamics(File, dyn_out)

   ! arguments
   type(File_desc_t), intent(inout) :: File
   type(dyn_export_t), intent(in)  :: dyn_out
   !----------------------------------------------------------------------------

   call cam_mpas_write_restart(restart_stream, endrun)

end subroutine write_restart_dynamics

!=========================================================================================

subroutine read_restart_dynamics(File, dyn_in, dyn_out)

   ! arguments
   type(File_desc_t),  intent(inout) :: File
   type(dyn_import_t), intent(out)   :: dyn_in
   type(dyn_export_t), intent(out)   :: dyn_out

   ! local variables
   character(len=*), parameter :: subname = 'restart_dynamics::read_restart_dynamics'
   integer :: ierr
   !----------------------------------------------------------------------------

   ! Before setting up the restart stream, names for each scalar constitutent must be defined
   call cam_mpas_define_scalars(domain_ptr % blocklist, dyn_in % mpas_from_cam_cnst, &
                                dyn_out % cam_from_mpas_cnst, ierr)
   if (ierr /= 0) then
      call endrun(subname//': Set-up of constituents for MPAS-A dycore failed.')
   end if

   call cam_mpas_setup_restart(file, restart_stream, MPAS_IO_READ, endrun)
   call cam_mpas_read_restart(restart_stream, endrun)

   ! Finish initializing the dynamics
   call dyn_init(dyn_in, dyn_out)

end subroutine read_restart_dynamics

end module restart_dynamics
