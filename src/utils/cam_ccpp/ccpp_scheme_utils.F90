module ccpp_scheme_utils

   ! Module of utilities available to CCPP schemes; CAM stubs to enable CCPPized schemes to build

   implicit none
   private

   !! Public interfaces
   public :: ccpp_constituent_index          ! Lookup index constituent by name
   public :: ccpp_constituent_indices        ! Lookup indices of consitutents by name

contains
   subroutine ccpp_constituent_index(standard_name, const_index, errcode, errmsg)
      ! Dummy arguments
      character(len=*),           intent(in)  :: standard_name
      integer,                    intent(out) :: const_index
      integer,          optional, intent(out) :: errcode
      character(len=*), optional, intent(out) :: errmsg

      ! Local variable
      character(len=*), parameter :: subname = 'ccpp_constituent_index'

      ! STUB DOES NOTHING

   end subroutine ccpp_constituent_index

   subroutine ccpp_constituent_indices(standard_names, const_inds, errcode, errmsg)
      ! Dummy arguments
      character(len=*),           intent(in)  :: standard_names(:)
      integer,                    intent(out) :: const_inds(:)
      integer,          optional, intent(out) :: errcode
      character(len=*), optional, intent(out) :: errmsg

      ! Local variables
      integer                     :: indx
      character(len=*), parameter :: subname = 'ccpp_constituent_indices'

      ! STUB DOES NOTHING

   end subroutine ccpp_constituent_indices

end module ccpp_scheme_utils
