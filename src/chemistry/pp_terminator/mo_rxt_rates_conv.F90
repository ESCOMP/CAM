module mo_rxt_rates_conv
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  private
  public :: set_rates
contains
   subroutine set_rates( rxt_rates, sol, ncol )
      real(r8), intent(inout) :: rxt_rates(:,:,:)
      real(r8), intent(in) :: sol(:,:,:)
      integer, intent(in) :: ncol
      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,     2)                                                ! rate_const*CL2
      rxt_rates(:ncol,:,     2) = rxt_rates(:ncol,:,     2)*sol(:ncol,:,     1)*sol(:ncol,:,     1)                            ! rate_const*CL*CL
  end subroutine set_rates
end module mo_rxt_rates_conv
