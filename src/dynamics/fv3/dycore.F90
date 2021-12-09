module dycore

    implicit none
    private

    public :: dycore_is

!=======================================================================
contains
!=======================================================================

logical function dycore_is(name)

    character(len=*) :: name
  
    dycore_is = .false.
    if (name == 'unstructured' .or. name == 'UNSTRUCTURED' .or. name == 'fv3' .or. name == 'FV3') then
        dycore_is = .true.
    end if
  
    return
end function dycore_is

end module dycore
