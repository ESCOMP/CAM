module radiation

! stub module

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl

implicit none
private
save

public :: &
   radiation_readnl,         &
   radiation_do

real(r8), public, protected :: nextsw_cday = -1._r8 ! future radiation calday for surface models

!========================================================================================
contains
!========================================================================================

subroutine radiation_readnl(nlfile)

   ! this stub can be called, but does nothing

   character(len=*), intent(in) :: nlfile

end subroutine radiation_readnl

!========================================================================================

function radiation_do(op, timestep)

   ! Returns true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value
   !---------------------------------------------------------------------------

   radiation_do = .false.

end function radiation_do

!========================================================================================

end module radiation

