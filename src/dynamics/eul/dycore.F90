module dycore

implicit none
private

public :: dycore_is

!=========================================================================================
CONTAINS
!=========================================================================================

logical function dycore_is(name)

   character(len=*), intent(in) :: name
      
   if (name == 'eul' .or. name == 'EUL') then
      dycore_is = .true.
   else
      dycore_is = .false.
   end if
      
end function dycore_is

!=========================================================================================

end module dycore


