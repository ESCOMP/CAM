module vdiff_lu_solver

! This module provides a function returning the matrix decomposition for
! an implicit finite volume solver for vertical diffusion. It accepts
! diffusion coefficients, time/grid spacing, and boundary condition
! objects, and returns a TriDiagDecomp object that can be used to diffuse
! an array for one time step with the "left_div" method.

use coords_1d, only: Coords1D
use linear_1d_operators, only: TriDiagOp, operator(+), TriDiagDecomp

implicit none
private
save

! Public interfaces
public :: vd_lu_decomp
public :: fin_vol_lu_decomp

! 8-byte real.
integer, parameter :: r8 = selected_real_kind(12)

contains

! ========================================================================!

! Designed to solve the equation:
! dq/dt = c1 q'' + c2 q' + c q

function vd_lu_decomp(dt, dp, coef_q,  coef_q_d, coef_q_d2, upper_bndry, &
     lower_bndry) result(decomp)

  use linear_1d_operators, only: &
       identity_operator, &
       diagonal_operator, &
       first_derivative, &
       second_derivative, &
       BoundaryType

  ! ---------------------- !
  ! Input-Output Arguments !
  ! ---------------------- !

  ! Time step.
  real(r8), intent(in) :: dt
  ! Grid spacing (deltas).
  real(r8), USE_CONTIGUOUS intent(in) :: dp(:,:)

  ! Coefficients for q, q', and q''.
  real(r8), USE_CONTIGUOUS intent(in), optional :: coef_q(:,:), &
       coef_q_d(:,:), coef_q_d2(:,:)

  ! Boundary conditions (optional, default to 0 flux through boundary).
  class(BoundaryType), target, intent(in), optional :: &
       upper_bndry, lower_bndry

  ! Output decomposition.
  type(TriDiagDecomp) :: decomp

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  ! Operator objects.
  type(TriDiagOp) :: add_term
  type(TriDiagOp) :: net_operator

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  if (present(coef_q)) then
     net_operator = diagonal_operator(1._r8 - dt*coef_q)
  else
     net_operator = identity_operator(size(dp, 1), size(dp, 2) + 1)
  end if

  if (present(coef_q_d)) then
     add_term = first_derivative(dp, upper_bndry, lower_bndry)
     call add_term%lmult_as_diag(-dt*coef_q_d)
     call net_operator%add(add_term)
  end if

  if (present(coef_q_d2)) then
     add_term = second_derivative(dp, upper_bndry, lower_bndry)
     call add_term%lmult_as_diag(-dt*coef_q_d2)
     call net_operator%add(add_term)
  end if

  decomp = TriDiagDecomp(net_operator)

  call net_operator%finalize()
  call add_term%finalize()

end function vd_lu_decomp

! ========================================================================!

! Designed to solve the equation:
!
! w * dq/dt = d/dp (D q' - v q) + c q
!
! where q is a grid-cell average, and p is the vertical coordinate
! (presumably pressure).
!
! In this function, coef_q_weight == w, coef_q_diff == D,
! coef_q_adv == v, and coef_q == c. All these are optional; omitting a
! coefficient is equivalent to setting the entire array to 0.
!
! coef_q_diff and coef_q_adv are defined at the level interfaces, while
! coef_q and coef_q_weight are grid-cell averages.

function fin_vol_lu_decomp(dt, p, coef_q, coef_q_diff, coef_q_adv, &
     coef_q_weight, upper_bndry, lower_bndry, graft_decomp) result(decomp)

  use linear_1d_operators, only: &
       zero_operator, &
       diagonal_operator, &
       diffusion_operator, &
       advection_operator, &
       BoundaryType

  ! ---------------------- !
  ! Input-Output Arguments !
  ! ---------------------- !

  ! Time step.
  real(r8), intent(in) :: dt
  ! Grid spacings.
  type(Coords1D), intent(in) :: p

  ! Coefficients for diffusion and advection.
  !
  ! The sizes must be consistent among all the coefficients that are
  ! actually present, i.e. coef_q_diff and coef_q_adv should be one level
  ! bigger than coef_q and coef_q_weight, and have the same column number.
  real(r8), USE_CONTIGUOUS intent(in), optional :: coef_q(:,:), &
       coef_q_diff(:,:), coef_q_adv(:,:), coef_q_weight(:,:)

  ! Boundary conditions (optional, default to 0 flux through boundary).
  class(BoundaryType), target, intent(in), optional :: &
       upper_bndry, lower_bndry

  ! Decomposition to graft onto. If this is provided, you can pass in
  ! smaller coefficients.
  type(TriDiagDecomp), intent(in), optional :: graft_decomp

  ! Output decomposition.
  type(TriDiagDecomp) :: decomp

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  ! Operator objects.
  type(TriDiagOp) :: add_term
  type(TriDiagOp) :: net_operator

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  ! A diffusion term is probably present, so start with that. Otherwise
  ! start with an operator of all 0s.

  if (present(coef_q_diff)) then
     net_operator = diffusion_operator(p, coef_q_diff, &
          upper_bndry, lower_bndry)
  else
     net_operator = zero_operator(p%n, p%d)
  end if

  ! Constant term (damping).
  if (present(coef_q)) then
     add_term = diagonal_operator(coef_q)
     call net_operator%add(add_term)
  end if

  ! Effective advection.
  if (present(coef_q_adv)) then
     add_term = advection_operator(p, coef_q_adv, &
          upper_bndry, lower_bndry)
     call net_operator%add(add_term)
  end if

  ! We want I-dt*(w^-1)*A for a single time step, implicit method, where
  ! A is the right-hand-side operator (i.e. what net_operator is now).
  if (present(coef_q_weight)) then
     call net_operator%lmult_as_diag(-dt/coef_q_weight)
  else
     call net_operator%lmult_as_diag(-dt)
  end if
  call net_operator%add_to_diag(1._r8)

  ! Decompose, grafting on an optional input decomp. The graft is a way to
  ! avoid re-calculating the ending (bottom) levels when the coefficients
  ! have only changed at the beginning (top), e.g. for different
  ! constituents in the molecular diffusion.
  decomp = TriDiagDecomp(net_operator, graft_decomp=graft_decomp)

  ! Ensure local objects are deallocated.
  call net_operator%finalize()
  call add_term%finalize()

end function fin_vol_lu_decomp

end module vdiff_lu_solver
