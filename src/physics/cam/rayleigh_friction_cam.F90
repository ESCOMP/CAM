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
use rayleigh_friction, only: rayleigh_friction_init, rayleigh_friction_run

implicit none
private
save
  
! Public interfaces
public :: &
     rayleigh_friction_readnl,  & ! read namelist
     rayleigh_friction_initO,   & ! Initialization
     rayleigh_friction_tendO      ! Computation of tendencies

! Namelist variables
integer  :: rayk0 = 2           ! vertical level at which rayleigh friction term is centered
real(r8) :: raykrange = 0._r8   ! range of rayleigh friction profile 
                                ! if 0, range is set to satisfy x=2 (see below)
real(r8) :: raytau0 = 0._r8     ! approximate value of decay time at model top (days)
                                ! if 0., no rayleigh friction is applied
! Local
real (r8) :: krange         ! range of rayleigh friction profile 
real (r8) :: tau0           ! approximate value of decay time at model top
real (r8) :: otau0          ! inverse of tau0
real (r8) :: otau(pver)     ! inverse decay time versus vertical level

!===============================================================================
contains
!===============================================================================

subroutine rayleigh_friction_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'rayleigh_friction_readnl'

   namelist /rayleigh_friction_nl/ rayk0, raykrange, raytau0
   !-----------------------------------------------------------------------------

   if (use_simple_phys) return

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

   if (masterproc) then
      if (raytau0 > 0._r8) then
         write (iulog,*) 'Rayleigh friction options: '
         write (iulog,*) '  rayk0     = ', rayk0
         write (iulog,*) '  raykrange = ', raykrange
         write (iulog,*) '  raytau0   = ', raytau0
      else
         write (iulog,*) 'Rayleigh friction not enabled.'
      end if
   end if

end subroutine rayleigh_friction_readnl

!=========================================================================================

subroutine rayleigh_friction_initO()

   !---------------------------Local storage-------------------------------
   character(len=512) errmsg
   integer errflg

   call rayleigh_friction_init(pver, raytau0, raykrange, rayk0, masterproc, iulog, errmsg, errflg)
   if (errflg /= 0) call endrun(errmsg)

end subroutine rayleigh_friction_initO
  
!=========================================================================================

subroutine rayleigh_friction_tendO(                                     &
     ztodt    ,state    ,ptend    )

   !-----------------------------------------------------------------------
   ! compute tendencies for rayleigh friction
   !-----------------------------------------------------------------------
   use physics_types, only: physics_state, physics_ptend, physics_ptend_init

   !------------------------------Arguments--------------------------------
   real(r8),            intent(in) :: ztodt      ! physics timestep
   type(physics_state), intent(in) :: state      ! physics state variables
    
   type(physics_ptend), intent(out):: ptend      ! individual parameterization tendencies
  !---------------------------Local storage-------------------------------
   character(len=512) errmsg
   integer errflg
   integer ncol                                ! number of atmospheric columns
   real(r8) rztodt                             ! 1./ztodt

   call physics_ptend_init(ptend, state%psetcols, 'rayleigh friction', ls=.true., lu=.true., lv=.true.)

   if (otau0 .eq. 0._r8) return

   ncol  = state%ncol

   call rayleigh_friction_run(ncol, pver, ztodt, state%u, state%v, ptend%u, ptend%v, ptend%s, errmsg, errflg)
   if (errflg /= 0) call endrun(errmsg)

end subroutine rayleigh_friction_tendO

end module rayleigh_friction_cam
