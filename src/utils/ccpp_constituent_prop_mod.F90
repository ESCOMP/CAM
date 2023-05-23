! This module is a placeholder for the CCPP generated module of the same name
module ccpp_constituent_prop_mod

   implicit none

   !Define stub version of constituent properties mod
   type, public :: ccpp_constituent_prop_ptr_t
   contains
      procedure :: standard_name => ccp_get_standard_name
   end type ccpp_constituent_prop_ptr_t

contains

   subroutine ccp_get_standard_name(this, std_name, errcode, errmsg)
   ! Return this constituent's standard name

   ! Dummy arguments
   class(ccpp_constituent_prop_ptr_t),   intent(in)  :: this
   character(len=*),                     intent(out) :: std_name
   integer,          optional,           intent(out) :: errcode
   character(len=*), optional,           intent(out) :: errmsg

   std_name = 'Not Used!'

   if(present(errcode)) then
      errcode = 0
   end if
   if(present(errmsg)) then
      errmsg = 'Still Not Used!'
   end if

   end subroutine ccp_get_standard_name

end module ccpp_constituent_prop_mod
