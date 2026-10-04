!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum --- Calculate x-sum bit-wise consistently
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum(grid, a, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8
      use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded
      use cam_logfile,   only : iulog
      use FVperf_module, only : FVstartclock, FVstopclock

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      ! input vector to be summed
      real (r8), intent(in) :: a(grid%ifirstxy:grid%ilastxy,ltot)

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)               ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !CALLED FROM:
!     te_map
!
! !REVISION HISTORY:
!
!     AAM 00.11.01 : Created
!     WS  03.10.22 : pmgrid removed (now spmd_dyn)
!     WS  04.10.04 : added grid as an argument; removed spmd_dyn
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     PW  08.06.25 : added fixed point reproducible sum
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8) :: rel_diff(2,ltot)

      integer :: im,lim

      logical :: write_warning
      logical :: tol_exceeded

      im  = grid%im
      lim = (grid%ilastxy-grid%ifirstxy) + 1

      call FVstartclock(grid,'xsum_reprosum')
      call shr_reprosum_calc(a, sum, lim, lim, ltot, gbl_count=im, &
                     commid=grid%commxy_x, rel_diff=rel_diff)
      call FVstopclock(grid,'xsum_reprosum')

      ! Warn if the nonreproducible floating point check sum differs from the
      ! integer vector sum by more than reprosum_diffmax. The integer vector sum is
      ! exact, so the result is kept regardless.
      write_warning = .false.
      if (grid%myidxy_x == 0) write_warning = .true.
      tol_exceeded = shr_reprosum_tolExceeded('par_xsum', ltot, write_warning, &
           iulog, rel_diff)

      return
!EOC
      end subroutine par_xsum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: par_xsum_r4 --- Calculate x-sum bit-wise consistently (real4)
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine par_xsum_r4(grid, a, ltot, sum)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      use dynamics_vars, only : T_FVDYCORE_GRID
      use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
      use shr_reprosum_mod, only : shr_reprosum_calc, shr_reprosum_tolExceeded
      use cam_logfile,   only : iulog
      use FVperf_module, only : FVstartclock, FVstopclock

      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      integer, intent(in) :: ltot       ! number of quantities to be summed
      real (r4) a(grid%ifirstxy:grid%ilastxy,ltot)    ! input vector to be summed

! !OUTPUT PARAMETERS:
      real (r8) sum(ltot)               ! sum of all vector entries

! !DESCRIPTION:
!     This subroutine calculates the sum of "a" in a reproducible
!     (sequentialized) fashion which should give bit-wise identical
!     results irrespective of the number of MPI processes.
!
! !REVISION HISTORY:
!
!     WS  05.04.08 : Created from par_xsum
!     WS  05.05.25 : removed ifirst, ilast, im as arguments (in grid)
!     WS  06.06.28 : Fixed bug in sequential version
!     PW  08.06.25 : added fixed point reproducible sum
!
!EOP
!---------------------------------------------------------------------
!BOC
 
! !Local
      real(r8) :: a8(grid%ifirstxy:grid%ilastxy,ltot)
      real(r8) :: rel_diff(2,ltot)

      integer :: im,lim

      logical :: write_warning
      logical :: tol_exceeded

      im  = grid%im
      lim = (grid%ilastxy-grid%ifirstxy) + 1

      call FVstartclock(grid,'xsum_r4_reprosum')
      a8(:,:) = a(:,:)
      call shr_reprosum_calc(a8, sum, lim, lim, ltot, gbl_count=im, &
                     commid=grid%commxy_x, rel_diff=rel_diff)
      call FVstopclock(grid,'xsum_r4_reprosum')

      ! Warn if the nonreproducible floating point check sum differs from the
      ! integer vector sum by more than reprosum_diffmax. The integer vector sum is
      ! exact, so the result is kept regardless.
      write_warning = .false.
      if (grid%myidxy_x == 0) write_warning = .true.
      tol_exceeded = shr_reprosum_tolExceeded('par_xsum_r4', ltot, write_warning, &
           iulog, rel_diff)

      return
!EOC
      end subroutine par_xsum_r4
!-----------------------------------------------------------------------
