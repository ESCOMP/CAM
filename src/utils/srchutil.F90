module srchutil

implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: Module containing Fortran equivalents to Cray library functions
!          NOTE: some aspects of this code may not meet the CCM coding standard
! 
! Author: Lifted from Cray manuals
! 
!-----------------------------------------------------------------------
#if (! defined UNICOSMP )

CONTAINS

!===============================================================================

   subroutine whenieq (n, array, inc, target, index, nval)
!----------------------------------------------------------------------- 
! 
! Purpose: Determine indices of "array" which equal "target"
! 
!-----------------------------------------------------------------------

!
! Arguments
!
      integer, intent(in) :: array(*)    ! array to be searched
      integer, intent(in) :: target      ! value to compare against
      integer, intent(in) :: inc         ! increment to move through array

      integer, intent(out) :: nval       ! number of values meeting criteria
      integer, intent(out) :: index(*)   ! output array of indices
!
! Local workspace
!
      integer :: i
      integer :: n
      integer :: ina

      ina=1
      nval=0
      if (inc .lt. 0) ina=(-inc)*(n-1)+1
      do i=1,n
         if(array(ina) .eq. target) then
           nval=nval+1
           index(nval)=i
         end if
         ina=ina+inc
      enddo
      return
   end subroutine whenieq

#endif

end module srchutil
