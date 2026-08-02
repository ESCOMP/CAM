module ccpp_scheme_utils

   ! Module of utilities available to CCPP schemes.
   ! CAM implementation using ccpp_const_props from ccpp_constituent_prop_mod.

   use ccpp_constituent_prop_mod, only: ccpp_const_props, int_unassigned

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

      ! Local variables
      character(len=*), parameter :: subname = 'ccpp_constituent_index'
      character(len=256)          :: t_const_name
      integer                     :: idx, ierr
      character(len=512)          :: local_errmsg

      const_index = int_unassigned

      if (.not. allocated(ccpp_const_props)) then
         if (present(errcode)) errcode = 1
         if (present(errmsg)) errmsg = subname//': constituent properties not initialized'
         return
      end if

      do idx = lbound(ccpp_const_props, 1), ubound(ccpp_const_props, 1)
         call ccpp_const_props(idx)%standard_name(t_const_name, ierr, &
              local_errmsg)
         if (ierr /= 0) then
            if (present(errcode)) errcode = ierr
            if (present(errmsg)) errmsg = trim(local_errmsg)
            return
         end if
         if (trim(t_const_name) == trim(standard_name)) then
            const_index = idx
            if (present(errcode)) errcode = 0
            if (present(errmsg)) errmsg = ''
            return
         end if
      end do

      if (present(errcode)) errcode = 1
      if (present(errmsg)) errmsg = trim(subname) // ': constituent "' // &
           trim(standard_name) // '" not found'

   end subroutine ccpp_constituent_index

   subroutine ccpp_constituent_indices(standard_names, const_inds, errcode, errmsg)
      ! Dummy arguments
      character(len=*),           intent(in)  :: standard_names(:)
      integer,                    intent(out) :: const_inds(:)
      integer,          optional, intent(out) :: errcode
      character(len=*), optional, intent(out) :: errmsg

      ! Local variables
      integer                     :: indx
      integer                     :: local_errcode
      character(len=512)          :: local_errmsg
      character(len=*), parameter :: subname = 'ccpp_constituent_indices'

      const_inds = int_unassigned

      if (.not. allocated(ccpp_const_props)) then
         if (present(errcode)) errcode = 1
         if (present(errmsg)) errmsg = subname//': constituent properties not initialized'
         return
      end if

      if (size(const_inds) < size(standard_names)) then
         if (present(errcode)) errcode = 1
         if (present(errmsg)) then
            write(errmsg, '(3a)') subname, ': const_inds array too small. ', &
                 'Must be >= size of standard_names'
         end if
         return
      end if

      do indx = 1, size(standard_names)
         call ccpp_constituent_index(standard_names(indx), &
              const_inds(indx), local_errcode, local_errmsg)
         if (local_errcode /= 0) then
            if (present(errcode)) errcode = local_errcode
            if (present(errmsg)) errmsg = trim(local_errmsg)
            return
         end if
      end do

      if (present(errcode)) errcode = 0
      if (present(errmsg)) errmsg = ''

   end subroutine ccpp_constituent_indices

end module ccpp_scheme_utils
