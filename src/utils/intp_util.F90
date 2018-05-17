module intp_util

implicit none

private

public :: findplb

contains

!#######################################################################

subroutine findplb( x, nx, xval, index )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8

   integer, intent(in) ::   nx         ! size of x
   real(r8), intent(in) ::  x(nx)      ! strictly increasing array
   real(r8), intent(in) ::  xval       ! value to be searched for in x
   
   integer, intent(out) ::  index

   ! Local variables:
   integer i
   !-----------------------------------------------------------------------

   if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
      index = nx
      return
   end if

   do i = 2, nx
      if ( xval .lt. x(i) ) then
         index = i-1
         return
      end if
   end do

end subroutine findplb

end module intp_util
