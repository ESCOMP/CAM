module rayleigh_friction_cam

!---------------------------------------------------------------------------------
! This contains the residual code required to read namelists for the
! Rayliegh Friction scheme. All of the functional code (the init and run subroutines)
! has been moved to ncar_ccpp code.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pver
use spmd_utils,        only: masterproc
use phys_control,      only: use_simple_phys
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun

implicit none
private
save
  
! Public interfaces
public :: rayleigh_friction_readnl  ! read namelist
! Rayleigh friction namelist parameters for use in physpkg                                                                                            
integer, public  :: rf_nl_k0 = 2           ! vertical level at which rayleigh friction term is centered                                               
real(r8), public :: rf_nl_krange = 0._r8   ! range of rayleigh friction profile                                                                       
                                           ! if 0, range is set to satisfy x=2 (see below)                                                            
real(r8), public :: rf_nl_tau0 = 0._r8     ! approximate value of decay time at model top (days)                                                      
                                           ! if 0., no rayleigh friction is applied                                                                   

!===============================================================================
contains
!===============================================================================

subroutine rayleigh_friction_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, rayk0
   real (r8) :: raykrange, raytau0
   character(len=*), parameter :: sub = 'rayleigh_friction_readnl'

   namelist /rayleigh_friction_nl/ rayk0, raykrange, raytau0
   !-----------------------------------------------------------------------------

   if (use_simple_phys) return

   ! Initialize with default values
   rayk0 = rf_nl_k0
   raykrange = rf_nl_krange
   raytau0 = rf_nl_tau0

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'rayleigh_friction_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, rayleigh_friction_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': FATAL: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(rayk0, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rayk0")
   call mpi_bcast(raykrange, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: raykrange")
   call mpi_bcast(raytau0, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: raytau0")

   ! Set module variables
   rf_nl_tau0 = raytau0
   rf_nl_krange = raykrange
   rf_nl_k0 = rayk0
   
   if (masterproc) then
      if (raytau0 > 0._r8) then
         write (iulog,*) 'Rayleigh friction options: '
         write (iulog,*) '  rayk0     = ', rf_nl_k0
         write (iulog,*) '  raykrange = ', rf_nl_krange
         write (iulog,*) '  raytau0   = ', rf_nl_tau0
      else
         write (iulog,*) 'Rayleigh friction not enabled.'
      end if
   end if

end subroutine rayleigh_friction_readnl

!=========================================================================================

end module rayleigh_friction_cam
