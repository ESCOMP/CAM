! Utility module for aerosol species types
module aerosol_spec_utils

  implicit none
  private
  public :: spec_type_in_list

contains

  ! returns TRUE if species type is found in type_list
  logical function spec_type_in_list(spec_type, type_list)
    character(len=*), intent(in) :: spec_type
    character(len=*), intent(in) :: type_list(:)
    integer :: i
    spec_type_in_list = .false.
    do i = 1, size(type_list)
       if (len_trim(type_list(i)) == 0) cycle
       if (trim(spec_type) == trim(type_list(i))) then
          spec_type_in_list = .true.
          return
       end if
    end do
  end function spec_type_in_list

end module aerosol_spec_utils
