module spcam_drivers

! stub module

use shr_kind_mod,     only: r8 => shr_kind_r8
use physics_types,    only: physics_state, physics_tend
use physics_buffer,   only: physics_buffer_desc
use camsrfexch,       only: cam_out_t, cam_in_t
use cam_abortutils,   only: endrun

implicit none
private
save

public :: tphysbc_spcam, spcam_register, spcam_init

!========================================================================================
contains
!========================================================================================

subroutine tphysbc_spcam (ztodt, state,   &
       tend,    pbuf,                     &
       cam_out, cam_in )

   real(r8),            intent(in)    :: ztodt
   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),     intent(inout) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in
   !---------------------------------------------------------------------------

   call endrun('tphysbc_spcam: ERROR: this is a stub')

end subroutine tphysbc_spcam

!========================================================================================

subroutine spcam_register()

end subroutine spcam_register

!========================================================================================

subroutine spcam_init(pbuf2d)
   
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

end subroutine spcam_init

!========================================================================================

end module spcam_drivers

