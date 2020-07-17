module dycore

    implicit none

    public :: dycore_is


contains


    !-----------------------------------------------------------------------
    !  routine dycore_is
    !
    !> \brief Indicates whether the dycore is the named entity
    !> \details
    !>  Returns true if the dycore is the named entity or possesses the named
    !>  quality, and false otherwise.
    !
    !-----------------------------------------------------------------------
    logical function dycore_is(name)

        character(len=*), intent(in) :: name

        dycore_is = .false.

        if (name == 'unstructured' .or. &
            name == 'UNSTRUCTURED' .or. &
            name == 'mpas' .or. &
            name == 'MPAS') then

            dycore_is = .true.

        end if

    end function dycore_is

end module dycore
