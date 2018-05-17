module dycore

implicit none
private
save

public :: dycore_is

!=========================================================================================
contains
!=========================================================================================

logical function dycore_is (name)

   ! Determine the dynamical core in use.

   character(len=*) :: name
   !-----------------------------------------------------------------------

   if (name == 'lr' .or. name == 'LR') then
      dycore_is = .true.
   else
      dycore_is = .false.
   end if

end function dycore_is

!=========================================================================================

end module dycore
