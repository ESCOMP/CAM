module cldfrc2m

! two-moment cloud fraction calculations - CAM interface
! Subroutine implementations are now in compute_cloud_fraction_two_moment module.
! this module only exists to provide public constants for the rest of CAM

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use compute_cloud_fraction_two_moment, only: CAMstfrac

implicit none
private
save

public ::           &
   cldfrc2m_readnl, &
   cldfrc2m_init,   &
   CAMstfrac,       &
   rhmini_const,    &
   rhmaxi_const,    &
   rhminis_const,   &
   rhmaxis_const,   &
   rhminl_const,    &
   rhminl_adj_land_const, &
   rhminh_const

! Namelist variables
real(r8) :: cldfrc2m_rhmini            ! Minimum rh for ice cloud fraction > 0.
real(r8) :: cldfrc2m_rhmaxi
real(r8) :: cldfrc2m_rhminis           ! Minimum rh for ice cloud fraction > 0 in the stratsophere.
real(r8) :: cldfrc2m_rhmaxis
real(r8) :: cldfrc2m_qist_min          ! Minimum in-stratus IWC constraint [ kg/kg ]
real(r8) :: cldfrc2m_qist_max          ! Maximum in-stratus IWC constraint [ kg/kg ]
logical  :: cldfrc2m_do_subgrid_growth = .false.
logical  :: cldfrc2m_do_avg_aist_algs = .false.
! -------------------------- !
! Parameters for Ice Stratus !
! -------------------------- !
real(r8),  protected :: rhmini_const                 ! Minimum rh for ice cloud fraction > 0.
real(r8),  protected :: rhmaxi_const
real(r8),  protected :: rhminis_const                ! Minimum rh for ice cloud fraction > 0.
real(r8),  protected :: rhmaxis_const

! ----------------------------- !
! Parameters for Liquid Stratus !
! ----------------------------- !

real(r8),  protected :: rhminl_const              ! Critical RH for low-level liquid stratus clouds
real(r8),  protected :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
real(r8),  protected :: rhminh_const              ! Critical RH for high-level liquid stratus clouds

! Internal variables (used only by init)
real(r8)             :: premit                    ! Top    height for mid-level liquid stratus fraction
real(r8)             :: premib                    ! Bottom height for mid-level liquid stratus fraction
integer              :: iceopt                    ! option for ice cloud closure
real(r8)             :: icecrit                   ! Critical RH for ice clouds in Wilson & Ballard closure

!================================================================================================
contains
!================================================================================================

subroutine cldfrc2m_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, masterprocid, mpi_logical, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldfrc2m_readnl'

   namelist /cldfrc2m_nl/ cldfrc2m_rhmini, cldfrc2m_rhmaxi, cldfrc2m_rhminis, cldfrc2m_rhmaxis, cldfrc2m_do_subgrid_growth, &
                          cldfrc2m_qist_min, cldfrc2m_qist_max, cldfrc2m_do_avg_aist_algs
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldfrc2m_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldfrc2m_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      rhmini_const  = cldfrc2m_rhmini
      rhmaxi_const  = cldfrc2m_rhmaxi
      rhminis_const = cldfrc2m_rhminis
      rhmaxis_const = cldfrc2m_rhmaxis
   end if

   ! Broadcast namelist variables
   call mpi_bcast(rhmini_const,               1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(rhmaxi_const,               1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(rhminis_const,              1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(rhmaxis_const,              1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(cldfrc2m_qist_min,          1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(cldfrc2m_qist_max,          1, mpi_real8,  masterprocid, mpicom, ierr)
   call mpi_bcast(cldfrc2m_do_subgrid_growth, 1, mpi_logical,masterprocid, mpicom, ierr)
   call mpi_bcast(cldfrc2m_do_avg_aist_algs,  1, mpi_logical,masterprocid, mpicom, ierr)

end subroutine cldfrc2m_readnl

!================================================================================================

subroutine cldfrc2m_init()

   use cloud_fraction, only: cldfrc_getparams
   use physconst, only: rair
   use compute_cloud_fraction_two_moment, only: compute_cloud_fraction_two_moment_init

   character(len=512) :: errmsg
   integer :: errflg

   call cldfrc_getparams(rhminl_out=rhminl_const, rhminl_adj_land_out=rhminl_adj_land_const,  &
                         rhminh_out=rhminh_const, premit_out=premit, premib_out=premib, &
                         iceopt_out=iceopt, icecrit_out=icecrit)

   call compute_cloud_fraction_two_moment_init( &
       masterproc, iulog, &
       premit, premib, iceopt, icecrit, &
       cldfrc2m_qist_min, cldfrc2m_qist_max, cldfrc2m_do_subgrid_growth, cldfrc2m_do_avg_aist_algs, &
       rair, &
       errmsg, errflg)

   if( masterproc ) then
      write(iulog,*) 'cldfrc2m parameters:'
      write(iulog,*) '  rhminl            = ', rhminl_const
      write(iulog,*) '  rhminl_adj_land   = ', rhminl_adj_land_const
      write(iulog,*) '  rhminh            = ', rhminh_const
      write(iulog,*) '  premit            = ', premit
      write(iulog,*) '  premib            = ', premib
      write(iulog,*) '  iceopt            = ', iceopt
      write(iulog,*) '  icecrit           = ', icecrit
      write(iulog,*) '  rhmini            = ', rhmini_const
      write(iulog,*) '  rhmaxi            = ', rhmaxi_const
      write(iulog,*) '  rhminis           = ', rhminis_const
      write(iulog,*) '  rhmaxis           = ', rhmaxis_const
      write(iulog,*) '  do_subgrid_growth = ', cldfrc2m_do_subgrid_growth
      write(iulog,*) '  do_avg_aist_algs  = ', cldfrc2m_do_avg_aist_algs
      write(iulog,*) '  cldfrc2m_qist_min = ', cldfrc2m_qist_min
      write(iulog,*) '  cldfrc2m_qist_max = ', cldfrc2m_qist_max
   end if

end subroutine cldfrc2m_init

!================================================================================================

end module cldfrc2m
