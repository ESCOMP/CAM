!------------------------------------------------------------------------------
! ESMF error handler
!------------------------------------------------------------------------------
module esmf_check_error_mod
  use shr_kind_mod, only: cl=>SHR_KIND_CL
  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun
  use ESMF, only: ESMF_SUCCESS

  implicit none

  private
  public :: check_esmf_error

contains

  subroutine check_esmf_error( rc, errmsg )

    integer, intent(in) :: rc
    character(len=*), intent(in) :: errmsg

    character(len=cl) :: errstr

    if (rc /= ESMF_SUCCESS) then
       write(errstr,'(a,i6)') 'esmf_zonal_mod::'//trim(errmsg)//' -- ESMF ERROR code: ',rc
       if (masterproc) write(iulog,*) trim(errstr)
       call endrun(trim(errstr))
    end if

  end subroutine check_esmf_error

end module esmf_check_error_mod
