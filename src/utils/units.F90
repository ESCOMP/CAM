module units

use shr_file_mod,   only: shr_file_getUnit, shr_file_freeUnit

implicit none
private

public :: getunit, freeunit

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

integer function getunit(iu)

   ! return an available unit number for i/o

   integer, intent(in), optional :: iu   ! desired unit number

   getunit = shr_file_getUnit(iu)

end function getunit

!-------------------------------------------------------------------------------

subroutine freeunit(iu)

   ! release the unit

   integer, intent(in) :: iu       ! unit number to be freed

   call shr_file_freeUnit(iu)

end subroutine freeunit

end module units
