module dadadj_cam

! CAM interfaces for the dry adiabatic adjustment parameterization

use shr_kind_mod,    only: r8=>shr_kind_r8, cs=>shr_kind_cs, cm=>shr_kind_cm
use ppgrid,          only: pcols, pver, pverp
use constituents,    only: pcnst
use air_composition, only: cappav, cpairv
use physconst,       only: pi
use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
use phys_control,    only: use_simple_phys
use cam_abortutils,  only: endrun
use cam_logfile,     only: iulog
use error_messages,  only: handle_errmsg

use spmd_utils,      only: masterproc, masterprocid, mpicom, mpi_integer
use namelist_utils,  only: find_group_name
use units,           only: getunit, freeunit

use dadadj,          only: dadadj_init, dadadj_run

implicit none
private
save

public :: &
   dadadj_readnl, &
   dadadj_cam_init, &
   dadadj_tend

! Namelist variables
integer :: dadadj_nlvdry = 3  ! number of layers from top of model to apply the adjustment
integer :: dadadj_niter = 15  ! number of iterations for convergence

!===============================================================================
contains
!===============================================================================

subroutine dadadj_readnl(filein)

   character(len=cs), intent(in) :: filein ! Input namelist filename

   namelist /dadadj_nl/ dadadj_nlvdry, dadadj_niter

   integer                            :: unitn, ierr
   integer                            :: errflg             ! CCPP physics scheme error flag
   character(len=512)                 :: errmsg             ! CCPP physics scheme error message
   character(len=*), parameter        :: sub='dadadj_readnl'
   !------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(filein), status='old')
      call find_group_name(unitn, 'dadadj_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, dadadj_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun( sub//':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(dadadj_nlvdry, 1, mpi_integer, masterprocid, mpicom)
   call mpibcast(dadadj_niter, 1, mpi_integer, masterprocid, mpicom)
#endif

   call dadadj_init(dadadj_nlvdry, dadadj_niter, pver, errmsg, errflg)
   if (errflg /=0) then
      call endrun('dadadj_readnl: Error returned from dadadj_init: '//trim(errmsg))
   end if

   if (masterproc .and. .not. use_simple_phys) then
      write(iulog,*)'Dry adiabatic adjustment applied to top N layers; N=', &
           dadadj_nlvdry
      write(iulog,*)'Dry adiabatic adjustment number of iterations for convergence =', &
           dadadj_niter
   end if

end subroutine dadadj_readnl


!===============================================================================

subroutine dadadj_cam_init()
    use cam_history,   only: addfld

    call addfld('DADADJ_PD', (/ 'lev' /), 'A', 'probability', 'dry adiabatic adjustment probability')

end subroutine dadadj_cam_init


!===============================================================================

subroutine dadadj_tend(dt, state, ptend)
   use cam_history,   only: outfld

   real(r8),                  intent(in)  :: dt         ! Time step [s]
   type(physics_state),       intent(in)  :: state      ! Physics state variables
   type(physics_ptend),       intent(out) :: ptend      ! parameterization tendencies

   character(len=512)                     :: errstring  ! Error string
   character(len=512)                     :: errmsg     ! CCPP physics scheme error message
   character(len=64)                      :: scheme_name! CCPP physics scheme name (not used in CAM)
   integer                                :: icol_err
   integer                                :: lchnk
   integer                                :: ncol
   integer                                :: errflg     ! CCPP physics scheme error flag
   logical                                :: lq(pcnst)
   real(r8)                               :: dadpdf(pcols, pver)

   !------------------------------------------------------------------
   ncol  = state%ncol
   lchnk = state%lchnk
   lq(:) = .FALSE.
   lq(1) = .TRUE.
   call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)

   !REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   dadpdf = 0._r8
   ptend%s = 0._r8
   ptend%q = 0._r8
   !REMOVECAM_END

   ! dadadj_run returns t tend, we are passing the ptend%s array to receive the t tendency and will convert it to s
   ! before it is returned to CAM..
   call dadadj_run( &
        ncol, pver, dt, state%pmid(:ncol,:), state%pint(:ncol,:), state%pdel(:ncol,:), &
        state%t(:ncol,:), state%q(:ncol,:,1), cappav(:ncol,:,lchnk), cpairv(:ncol,:,lchnk), ptend%s(:ncol,:), &
        ptend%q(:ncol,:,1), dadpdf(:ncol,:), scheme_name, errmsg, errflg)

   ! error exit
   if (errflg /= 0) then
      ! If this is a Convergence error then output lat lon of problem column using column index (errflg)
      if(index('Convergence', errmsg) /= 0)then
         write(errstring, *) trim(adjustl(errmsg)),' lat:',state%lat(errflg)*180._r8/pi,' lon:', &
              state%lon(errflg)*180._r8/pi
      else
         errstring=trim(errmsg)
      end if
      call endrun('Error dadadj_tend:'//trim(errstring))
   end if

   call outfld('DADADJ_PD',  dadpdf(:ncol,:),  ncol, lchnk)

end subroutine dadadj_tend

!===============================================================================
end module dadadj_cam
