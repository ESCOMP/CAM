module epotential_params
  use shr_kind_mod, only: r8 => shr_kind_r8

  logical,  public :: epot_active = .false.
  real(r8), public :: epot_crit_colats(2) = -huge(1._r8)

end module epotential_params

