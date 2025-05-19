module pio_reader
   use pio,            only: file_desc_t
   use ccpp_io_reader, only: abstract_netcdf_reader_t
   use shr_kind_mod,   only: cl=>shr_kind_cl

   implicit none
   private

   public :: pio_reader_t

   !Error code parameters
   integer, parameter :: pio_inq_dim_id_err      = 1
   integer, parameter :: pio_inq_dim_len_err     = 2
   integer, parameter :: pio_inq_var_id_err      = 3
   integer, parameter :: pio_inq_var_info_err    = 4
   integer, parameter :: bad_var_rank_err        = 5
   integer, parameter :: too_high_rank_err       = 6
   integer, parameter :: pio_get_var_err         = 7
   integer, parameter :: not_char_type_err       = 8
   integer, parameter :: file_not_open_err       = 9
   integer, parameter :: pio_get_msg_err         = 10

   type :: file_handle_t
      logical            :: is_file_open = .false.  !Is NetCDF file currently open?
      type(file_desc_t)  :: pio_fh                  !PIO File handle type
      character(len=cl)  :: file_path = ''          !Local path to NetCDF file
   end type

   type, extends(abstract_netcdf_reader_t) :: pio_reader_t
      private
      type(file_handle_t) :: sima_pio_fh !PIO File handle type
   contains
      procedure :: open_file       => open_netcdf_file
      procedure :: close_file      => close_netcdf_file
      procedure :: get_var_int_0d  => get_netcdf_var_int_0d
      procedure :: get_var_int_1d  => get_netcdf_var_int_1d
      procedure :: get_var_int_2d  => get_netcdf_var_int_2d
      procedure :: get_var_int_3d  => get_netcdf_var_int_3d
      procedure :: get_var_int_4d  => get_netcdf_var_int_4d
      procedure :: get_var_int_5d  => get_netcdf_var_int_5d
      procedure :: get_var_real_0d => get_netcdf_var_real_0d
      procedure :: get_var_real_1d => get_netcdf_var_real_1d
      procedure :: get_var_real_2d => get_netcdf_var_real_2d
      procedure :: get_var_real_3d => get_netcdf_var_real_3d
      procedure :: get_var_real_4d => get_netcdf_var_real_4d
      procedure :: get_var_real_5d => get_netcdf_var_real_5d
      procedure :: get_var_char_0d => get_netcdf_var_char_0d
      procedure :: get_var_char_1d => get_netcdf_var_char_1d
      procedure :: get_var_char_2d => get_netcdf_var_char_2d
      procedure :: get_var_char_3d => get_netcdf_var_char_3d
      procedure :: get_var_char_4d => get_netcdf_var_char_4d
      procedure :: get_var_char_5d => get_netcdf_var_char_5d
   end type pio_reader_t

contains

   subroutine open_netcdf_file(this, file_path, errmsg, errcode)
      use ioFileMod,        only: getfil
      use cam_pio_utils,    only: cam_pio_openfile
      use pio,              only: PIO_NOWRITE

      class(pio_reader_t), intent(inout)  :: this
      character(len=*),    intent(in)  :: file_path
      integer,             intent(out) :: errcode
      character(len=*),    intent(out) :: errmsg

      character(len=cl)  :: local_file_path  !NetCDF file path on local file system

      if(this%sima_pio_fh%is_file_open) then
         errcode = 1
         errmsg = "Trying to reuse pio_reader already used for: '"//this%sima_pio_fh%file_path
         return
      end if

      if(file_path == 'UNSET_PATH') then
         errcode = 1
         errmsg = "Found UNSET_PATH trying to open file"
         return
      end if

      call getfil(file_path, local_file_path)
      call cam_pio_openfile(this%sima_pio_fh%pio_fh, local_file_path, PIO_NOWRITE)

      !Set file handle metadata
      this%sima_pio_fh%file_path    = local_file_path
      this%sima_pio_fh%is_file_open = .true.

      !File was successfully opened
      errcode = 0
      errmsg = ''
   end subroutine open_netcdf_file

   subroutine close_netcdf_file(this, errmsg, errcode)
      use pio, only: pio_closefile

      class(pio_reader_t), intent(inout)  :: this
      integer,             intent(out) :: errcode
      character(len=*),    intent(out) :: errmsg

      if(.not.this%sima_pio_fh%is_file_open) then
         !File is already closed, so just exit
         !quietly.
         errcode = 0
         errmsg = ''
         return
      end if

      !Close NetCDF File:
      call pio_closefile(this%sima_pio_fh%pio_fh)

      !Inidcate that file handle array id is no longer in use:
      this%sima_pio_fh%is_file_open = .false.
      this%sima_pio_fh%file_path = ''

      !File was successfully closed
      errcode = 0
      errmsg = ''
   end subroutine close_netcdf_file

   ! ------------------------------------------------------------------
   ! Integer interfaces
   ! ------------------------------------------------------------------

   subroutine get_netcdf_var_int_0d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var     !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 0) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      var = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var)

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_0d

   subroutine get_netcdf_var_int_1d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var(:)  !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 1) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension size
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:) = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var(:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_1d

   subroutine get_netcdf_var_int_2d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var(:,:) !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 2) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:) = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_2d

   subroutine get_netcdf_var_int_3d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var(:,:,:) !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 3) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:) = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_3d

   subroutine get_netcdf_var_int_4d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var(:,:,:,:) !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 4) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3), dim_sizes(4)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:) = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_4d

   subroutine get_netcdf_var_int_5d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t), intent(in)  :: this
      character(len=*),    intent(in)  :: varname
      integer, pointer,    intent(out) :: var(:,:,:,:,:) !Integer variable that file data will be read to.
      integer,             intent(out) :: errcode !Error code
      character(len=*),    intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 5) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3), dim_sizes(4), dim_sizes(5)), &
               stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:,:) = huge(1)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_int_5d

   ! ------------------------------------------------------------------
   ! Real interfaces
   ! ------------------------------------------------------------------

   subroutine get_netcdf_var_real_0d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var     !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 0) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      var = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var)

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_0d

   subroutine get_netcdf_var_real_1d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var(:) !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 1) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:) = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var(:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_1d

   subroutine get_netcdf_var_real_2d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var(:,:) !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 2) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:) = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_2d

   subroutine get_netcdf_var_real_3d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var(:,:,:) !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 3) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:) = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_3d

   subroutine get_netcdf_var_real_4d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var(:,:,:,:) !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 4) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3), dim_sizes(4)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:) = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_4d

   subroutine get_netcdf_var_real_5d(this, varname, var, errmsg, errcode)
      use ccpp_kinds, only: kind_phys
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR

      class(pio_reader_t),           intent(in)  :: this
      character(len=*),              intent(in)  :: varname
      real(kind_phys), pointer,      intent(out) :: var(:,:,:,:,:) !Real variable that file data will be read to.
      integer,                       intent(out) :: errcode !Error code
      character(len=*),              intent(out) :: errmsg  !Error message

      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes

      integer :: i !loop control variable
      !----------------------

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:
      errcode = 0
      if(ndims /= 5) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data:
      allocate(var(dim_sizes(1), dim_sizes(2), dim_sizes(3), dim_sizes(4), dim_sizes(5)), &
               stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:,:) = huge(1._kind_phys)
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_real_5d

   ! ------------------------------------------------------------------
   ! Character interfaces
   ! ------------------------------------------------------------------

   subroutine get_netcdf_var_char_0d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var     !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 1) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var, stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var)

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_0d

   subroutine get_netcdf_var_char_1d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var(:) !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 2) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var(dim_sizes(2)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:) = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var(:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_1d

   subroutine get_netcdf_var_char_2d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var(:,:) !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 3) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var(dim_sizes(2), dim_sizes(3)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:) = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_2d

   subroutine get_netcdf_var_char_3d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var(:,:,:) !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 4) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var(dim_sizes(2), dim_sizes(3), dim_sizes(4)), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:) = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_3d

   subroutine get_netcdf_var_char_4d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var(:,:,:,:) !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 5) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var(dim_sizes(2), dim_sizes(3), dim_sizes(4), dim_sizes(5)), &
               stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:) = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_4d

   subroutine get_netcdf_var_char_5d(this, varname, var, errmsg, errcode)
      use pio,        only: pio_inq_varid
      use pio,        only: pio_inq_dimlen
      use pio,        only: pio_inquire_variable
      use pio,        only: pio_seterrorhandling
      use pio,        only: pio_get_var
      use pio,        only: PIO_NOERR
      use pio,        only: PIO_BCAST_ERROR
      use pio_types,  only: PIO_CHAR

      class(pio_reader_t),       intent(in)  :: this
      character(len=*),          intent(in)  :: varname
      character(len=:), pointer, intent(out) :: var(:,:,:,:,:) !Character variable that file data will be read to.
      integer,                   intent(out) :: errcode !Error code
      character(len=*),          intent(out) :: errmsg  !Error message

      !Local variables:
      type(file_desc_t)    :: pio_file_handle !File handle type used by PIO
      character(len=cl)    :: file_path       !Path to NetCDF file
      integer              :: err_handling    !PIO error handling code
      integer              :: var_id          !NetCDF variable ID
      integer              :: nc_type         !NetCDF variable type
      integer              :: ndims           !Number of variable dimensions on NetCDF file
      integer, allocatable :: dim_ids(:)      !Variable dimension IDs
      integer, allocatable :: dim_sizes(:)    !Variable dimension sizes
      integer              :: i !loop control variable

      !Check if file is open:
      if(.not.this%sima_pio_fh%is_file_open) then
         !File isn't actually open, so throw an error
         errcode = file_not_open_err
         errmsg = "File '"//this%sima_pio_fh%file_path//"' is not open, need to call 'open_file' first."
         return
      end if

      !Extract open file information:
      pio_file_handle = this%sima_pio_fh%pio_fh
      file_path       = this%sima_pio_fh%file_path

      !Force PIO to send an error code instead of dying:
      call pio_seterrorhandling(pio_file_handle, PIO_BCAST_ERROR, oldmethod=err_handling)

      !Look for variable on file:
      errcode = pio_inq_varid(pio_file_handle, varname, var_id)
      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_id_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable type and number of variable dimensions on file:
      errcode = pio_inquire_variable(pio_file_handle, var_id, xtype=nc_type, ndims=ndims)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that variable is a character array
      !(as we cannot currently handle string-type variables):
      if(nc_type /= PIO_CHAR) then
         errcode = not_char_type_err
         errmsg = "NetCDF Variable '"//varname//"' is not a character array.  File can be found here: "//file_path
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Check that the variable rank as specified by the caller
      !matches what is found on the NetCDF file:

      !NOTE:  NetCDF supports both character arrays and string-type
      !data depending on the NetCDF version, so the dimensions
      !might be one larger than the actual array size if it
      !includes the character length as a dimension as well.
      !Ideally the actual type would be checked and handled
      !differently, but for now just assume a character array
      !and check for ndims = rank+1
      errcode = 0
      if(ndims /= 6) then
         errcode = bad_var_rank_err
         errmsg  = "Variable '"//varname//"' isn't declared with the correct number of dimensions"
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension sizes:
      !Allocate NetCDF variable dimension ID array:
      allocate(dim_ids(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get variable dimension IDs:
      errcode = pio_inquire_variable(pio_file_handle, var_id, dimids=dim_ids)
      if(errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_inq_var_info_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Allocate NetCDF variable dimension sizes array:
      allocate(dim_sizes(ndims), stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Get dimension sizes:
      do i = 1, ndims
         errcode = pio_inq_dimlen(pio_file_handle, dim_ids(i), dim_sizes(i))
         if(errcode /= PIO_NOERR) then
            !Extract error message from PIO:
            call get_pio_errmsg(pio_inq_dim_len_err, errcode, errmsg)

            !Reset PIO back to original error handling method:
            call pio_seterrorhandling(pio_file_handle, err_handling)
            return
         end if
      end do

      !Now attempt to allocate and initialize variable, and
      !read-in the NetCDF data.  Note that the first dimenstion
      !is the length of the character array, so need to start
      !the dim_sizes allocation count at index two:
      allocate(character(dim_sizes(1)) :: var(dim_sizes(2), dim_sizes(3), dim_sizes(4), dim_sizes(5), dim_sizes(6)), &
               stat=errcode, errmsg=errmsg)
      if(errcode /= 0) then
         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if
      var(:,:,:,:,:) = 'UNSET'
      errcode = pio_get_var(pio_file_handle, var_id, var(:,:,:,:,:))

      if (errcode /= PIO_NOERR) then
         !Extract error message from PIO:
         call get_pio_errmsg(pio_get_var_err, errcode, errmsg)

         !Reset PIO back to original error handling method:
         call pio_seterrorhandling(pio_file_handle, err_handling)
         return
      end if

      !Reset PIO back to original error handling method:
      call pio_seterrorhandling(pio_file_handle, err_handling)

      !Variable was successfully read, so properly set the error
      !code and message:
      errcode = 0
      errmsg = ''
   end subroutine get_netcdf_var_char_5d

   subroutine get_pio_errmsg(caller_errcode, errcode, errmsg)
      !Set error message based off PIO error code,
      !and then reset PIO error code to caller-specified
      !error code.

      !Note that if an internal error occurs when
      !attempting to grab the error message both
      !the error code and error message will be
      !updated.

      use pio, only: pio_strerror
      use pio, only: PIO_NOERR

      !Input/output arguments:
      integer,          intent(in)    :: caller_errcode !New error code caller wants.
      integer,          intent(inout) :: errcode        !Error code
      character(len=*), intent(inout) :: errmsg         !Error message

      !Local variables:
      integer :: strerr !Error code returned if pio_strerror fails

      strerr = pio_strerror(errcode, errmsg)
      if(strerr /= PIO_NOERR) then
         write(errmsg, *) "Failed to get error message for PIO code: ", errcode
         errcode = pio_get_msg_err
      else
         errcode = caller_errcode
      end if
   end subroutine get_pio_errmsg

end module pio_reader
