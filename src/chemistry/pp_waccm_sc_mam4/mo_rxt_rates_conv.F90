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
      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,    10)                                                ! rate_const*H2O2
      rxt_rates(:ncol,:,     2) = rxt_rates(:ncol,:,     2)*sol(:ncol,:,    26)                                                ! rate_const*soa_a1
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,    27)                                                ! rate_const*soa_a2
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,    30)                                                ! rate_const*H2O
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,    10)                                                ! rate_const*OH*H2O2
                                                                                                                               ! rate_const
      rxt_rates(:ncol,:,     7) = rxt_rates(:ncol,:,     7)*sol(:ncol,:,    12)                                                ! rate_const*N2O
      rxt_rates(:ncol,:,     8) = rxt_rates(:ncol,:,     8)*sol(:ncol,:,     3)                                                ! rate_const*CFC11
      rxt_rates(:ncol,:,     9) = rxt_rates(:ncol,:,     9)*sol(:ncol,:,     4)                                                ! rate_const*CFC12
      rxt_rates(:ncol,:,    10) = rxt_rates(:ncol,:,    10)*sol(:ncol,:,     5)                                                ! rate_const*CH4
      rxt_rates(:ncol,:,    11) = rxt_rates(:ncol,:,    11)*sol(:ncol,:,     6)                                                ! rate_const*NO3*DMS
      rxt_rates(:ncol,:,    12) = rxt_rates(:ncol,:,    12)*sol(:ncol,:,     6)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,    13) = rxt_rates(:ncol,:,    13)*sol(:ncol,:,    22)                                                ! rate_const*OH*M*SO2
      rxt_rates(:ncol,:,    14) = rxt_rates(:ncol,:,    14)*sol(:ncol,:,     6)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,    15) = rxt_rates(:ncol,:,    15)*sol(:ncol,:,    28)                                                ! rate_const*SOAE
  end subroutine set_rates
end module mo_rxt_rates_conv
