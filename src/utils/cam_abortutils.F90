module cam_abortutils

   use shr_kind_mod, only: SHR_KIND_CL
   use shr_sys_mod,  only: endrun => shr_sys_abort

   implicit none
   private
   save

   public :: endrun
   public :: handle_allocate_error

CONTAINS

   subroutine handle_allocate_error(retval, subname, fieldname, errmsg)
      ! if <retval> is not zero, generate an error message and abort
      ! Dummy arguments
      integer,          intent(in) :: retval
      character(len=*), intent(in) :: subname
      character(len=*), intent(in) :: fieldname
      character(len=*), optional, intent(in) :: errmsg

      ! Local variable
      character(len=SHR_KIND_CL)   :: abort_msg

      if (retval /= 0) then
         write(abort_msg, '(4a,i0)') trim(subname), ' error allocating ',      &
              trim(fieldname), ', error = ', retval
         if (present(errmsg)) then
            abort_msg = trim(abort_msg) // new_line('a') // &
                        "Allocation failed with: " // trim(errmsg)
         end if
         call endrun(abort_msg)
      end if
   end subroutine handle_allocate_error

end module cam_abortutils
