module interp_mod

! This is a stub module.  Online interpolation is not currently available.

use shr_kind_mod,        only: r8=>shr_kind_r8

use cam_history_support, only: interp_info_t

use pio,                 only: file_desc_t, var_desc_t

use cam_logfile,         only: iulog
use cam_abortutils,      only: endrun

implicit none
private
save

public :: setup_history_interpolation
public :: set_interp_hfile
public :: write_interpolated

interface write_interpolated
   module procedure write_interpolated_scalar
   module procedure write_interpolated_vector
end interface write_interpolated

integer, parameter :: nlat=0, nlon=0

!=========================================================================================
contains
!=========================================================================================

subroutine setup_history_interpolation(interp_ok, mtapes, interp_output, interp_info)

   logical,             intent(inout) :: interp_ok
   integer,             intent(in)    :: mtapes
   logical,             intent(in)    :: interp_output(:)
   type(interp_info_t), intent(inout) :: interp_info(:)
   !----------------------------------------------------------------------------

   interp_ok = .false.

   write (iulog,*) 'INFO - setup_history_interpolation is a no-op'

end subroutine setup_history_interpolation

!=========================================================================================

subroutine set_interp_hfile(hfilenum, interp_info)

   ! arguments
   integer,             intent(in)    :: hfilenum
   type(interp_info_t), intent(inout) :: interp_info(:)
   !----------------------------------------------------------------------------


end subroutine set_interp_hfile

!=========================================================================================

subroutine write_interpolated_scalar(File, varid, fld, numlev, data_type, decomp_type)

   type(file_desc_t), intent(inout) :: File
   type(var_desc_t), intent(inout) :: varid
   real(r8), intent(in) :: fld(:,:,:)
   integer, intent(in) :: numlev, data_type, decomp_type
   !----------------------------------------------------------------------------


   call endrun('FATAL - write_interpolated_scalar is a stub, you shouldnt get here')

end subroutine write_interpolated_scalar

!=========================================================================================

subroutine write_interpolated_vector(File, varidu, varidv, fldu, fldv, &
   numlev, data_type, decomp_type) 

   type(file_desc_t), intent(inout) :: File
   type(var_desc_t), intent(inout) :: varidu, varidv
   real(r8), intent(in) :: fldu(:,:,:), fldv(:,:,:)
   integer, intent(in) :: numlev, data_type, decomp_type
   !----------------------------------------------------------------------------

   call endrun('FATAL - write_interpolated_vector is a stub, you shouldnt get here')

end subroutine write_interpolated_vector

!=========================================================================================

end module interp_mod
