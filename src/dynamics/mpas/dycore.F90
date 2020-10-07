module dycore

implicit none

public :: dycore_is

!=========================================================================================
contains
!=========================================================================================

logical function dycore_is(name)

   ! Identifies that the MPAS dycore is being used.
   ! Identifies that the MPAS dycore uses an 'unstructured' grid.

   character(len=*), intent(in) :: name

   dycore_is = .false.

   if (name == 'unstructured' .or. &
       name == 'UNSTRUCTURED' .or. &
       name == 'mpas' .or. &
       name == 'MPAS') then

      dycore_is = .true.

   end if

end function dycore_is

!=========================================================================================

end module dycore
