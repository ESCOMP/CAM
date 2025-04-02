module atmos_phys_string_utils
    ! String utils

    implicit none
    private

    public :: to_lower
    public :: to_upper

contains

    pure function to_lower(input_string) result(lowercase_string)
       character(len=*), intent(in) :: input_string
       character(len=*)             :: lowercase_string
       ! Local variables

       integer :: i                ! Index
       integer :: aseq             ! ascii collating sequence
       integer :: upper_to_lower   ! integer to convert case
       character(len=1) :: ctmp    ! Character temporary
    !-----------------------------------------------------------------------
       upper_to_lower = iachar("a") - iachar("A")

       do i = 1, len(input_string)
          ctmp = input_string(i:i)
          aseq = iachar(ctmp)
          if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
               ctmp = achar(aseq + upper_to_lower)
          lowercase_string(i:i) = ctmp
       end do

    end function to_lower

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

    pure function to_upper(input_string) result(uppercase_string)
       character(len=*), intent(in) :: input_string
       character(len=*)             :: uppercase_string

       integer :: i                ! Index
       integer :: aseq             ! ascii collating sequence
       integer :: lower_to_upper   ! integer to convert case
       character(len=1) :: ctmp    ! Character temporary
    !-----------------------------------------------------------------------
       lower_to_upper = iachar("A") - iachar("a")

       do i = 1, len(input_string)
          ctmp = input_string(i:i)
          aseq = iachar(ctmp)
          if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
               ctmp = achar(aseq + lower_to_upper)
          uppercase_string(i:i) = ctmp
       end do

    end function to_upper

end module atmos_phys_string_utils
