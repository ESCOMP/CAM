module aurora_params

  use shr_kind_mod  ,only: r8 => shr_kind_r8

  implicit none

  ! global variables shared between column physics and waccmx ionosphere
  ! -- this assumes these share the same MPI tasks --

  logical  :: aurora_params_set = .false.
  real(r8) :: hpower = -huge(1.0_r8)
  real(r8) :: plevel = -huge(1.0_r8)
  real(r8) :: ctpoten = -huge(1.0_r8)
  real(r8) :: theta0(2) = -huge(1.0_r8)
  real(r8) :: offa(2) = -huge(1.0_r8)
  real(r8) :: dskofa(2) = -huge(1.0_r8)
  real(r8) :: phid(2) = -huge(1.0_r8)
  real(r8) :: rrad(2) = -huge(1.0_r8)
  real(r8) :: offc(2) = -huge(1.0_r8)
  real(r8) :: dskofc(2) = -huge(1.0_r8)
  real(r8) :: phin(2) = -huge(1.0_r8)

  logical  :: amie_period = .false. ! true during a period of prescribed high-latitude electric potential

end module aurora_params
