#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module restart_dynamics

use shr_kind_mod,   only: r8 => shr_kind_r8, r4=> shr_kind_r4
use spmd_utils,     only: iam

use pmgrid,         only: plev, plevp
use constituents,   only: pcnst

use dyn_grid,       only: get_horiz_grid_dim_d

use dyn_comp,       only: dyn_import_t, dyn_export_t, dyn_init

use pio,            only: pio_global, pio_double, pio_unlimited, &
                          pio_seterrorhandling, pio_bcast_error, pio_noerr, &
                          file_desc_t, io_desc_t, var_desc_t, &
                          pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                          pio_def_dim, pio_def_var, &
                          pio_put_att, pio_put_var, pio_write_darray, &
                          pio_get_att, pio_read_darray

use time_manager,   only: get_step_size, get_curr_time

use cam_logfile,    only: iulog
use shr_sys_mod,    only: shr_sys_flush
use perf_mod

use mpas_derived_types, only : MPAS_Stream_type

!    use mpas_cam_interface, only : MPAS_init_restart
!    use mpas_cam_interface, only : MPAS_write_restart

implicit none
private
save

public :: &
   init_restart_dynamics,  &
   write_restart_dynamics, &
   read_restart_dynamics

type(var_desc_t) :: PHISdesc, PSdesc, timedesc, Tdesc, Tracerdesc, UVPERPdesc
type(var_desc_t) :: UXdesc, UYdesc, OMEGAdesc
integer :: ncol_dimid, nlev_dimid, nlevp_dimid
real(r8), parameter ::  SECONDS_IN_DAY          =  86400_r8
logical, save :: initialized=.false.

! The restart_stream is set up in init_restart_dynamics and used later in write_restart_dynamics
type (MPAS_Stream_type) :: restart_stream


!=========================================================================================
contains
!=========================================================================================


!-----------------------------------------------------------------------
!  routine init_restart_dynamics
!
!> \brief Prepares PIO filehandle for later writing of dynamics restart state
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine init_restart_dynamics(file, dyn_out)

   use cam_mpas_subdriver, only : cam_mpas_setup_restart
   use mpas_derived_types, only : MPAS_IO_WRITE
   use cam_abortutils,     only : endrun

   type(file_desc_t), target :: File
   type(dyn_export_t), intent(in)  :: dyn_out

   integer :: ncols
   integer :: pcnst_dimid, time_dimid
   integer :: ierr

   character(len=*), parameter :: subname = 'restart_dynamics::init_restart_dynamics'


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   call cam_mpas_setup_restart(file, restart_stream, MPAS_IO_WRITE, endrun)

end subroutine init_restart_dynamics


!-----------------------------------------------------------------------
!  routine write_restart_dynamics
!
!> \brief Writes dynamics restart state to a PIO filehandle
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine write_restart_dynamics(File, dyn_out)

   use cam_mpas_subdriver, only : cam_mpas_write_restart
   use cam_abortutils,     only : endrun

   ! ARGUMENTS
   type(File_desc_t), intent(inout) :: File     ! Unit number
   type(dyn_export_t), intent(in)  :: dyn_out

   ! LOCALS
   type(io_desc_t)       :: iodesc2d, iodesc3d, iodesctr, iodescu, iodescw  ! I/O descriptors
   real(r8)              :: time                ! current time 
   integer               :: ierr                ! error flag
   real(kind=r8),pointer :: var3d(:), var2d(:), vartr(:), varu(:)
   integer :: i, j, k, icnt

   character(len=*), parameter :: subname = 'restart_dynamics::write_restart_dynamics'


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   call cam_mpas_write_restart(restart_stream, endrun)

end subroutine write_restart_dynamics


!-----------------------------------------------------------------------
!  routine read_restart_dynamics
!
!> \brief Reads dynamics restart state from PIO filehandle
!> \details
!>  More details go here...
!
!-----------------------------------------------------------------------
subroutine read_restart_dynamics(File, dyn_in, dyn_out)

   use cam_mpas_subdriver, only : domain_ptr, cam_mpas_setup_restart, cam_mpas_read_restart, &
                                  cam_mpas_define_scalars
   use mpas_derived_types, only : MPAS_IO_READ
   use cam_abortutils,     only : endrun

   ! ARGUMENTS
   type(File_desc_t),  intent(inout) :: File
   type(dyn_import_t), intent(out)   :: dyn_in
   type(dyn_export_t), intent(out)   :: dyn_out

   ! LOCALS
   type(io_desc_t)       :: iodesc2d, iodesc3d, iodesctr, iodescu, iodescw  ! I/O descriptors
   real(kind=r8),pointer :: var3d(:), var2d(:), vartr(:), varu(:)
   integer               :: i, j, icnt, ierr, k
   integer               :: fplev, fnCells, fnEdges, fnVertices, fnumcols, fpcnst

   character(len=*), parameter :: subname = 'restart_dynamics::read_restart_dynamics'


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   !
   ! Before setting up the restart stream, names for each scalar constitutent must be defined
   !
   if (cam_mpas_define_scalars(domain_ptr % blocklist, dyn_in % mpas_from_cam_cnst, &
                               dyn_out % cam_from_mpas_cnst) /= 0) then
      call endrun(subname//': Set-up of constituents for MPAS-A dycore failed.')
   end if

   call cam_mpas_setup_restart(file, restart_stream, MPAS_IO_READ, endrun)
   call cam_mpas_read_restart(restart_stream, endrun)

   ! Initialize the dynamics
   call dyn_init(dyn_in, dyn_out)

end subroutine read_restart_dynamics

end module restart_dynamics
