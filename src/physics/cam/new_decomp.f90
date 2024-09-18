module new_decomp

implicit none
private
save

public :: fin_vol_solve
public :: vd_lu_solve

contains

function vd_lu_solve(                                                     &
   pcols , pver   , ncol       ,u, fixed_ubc  , cnst_fixed_ubflx, mw, &
   kv    , kq_scal, mw_facm    , dpidz_sq   , p,                          &
   interface_boundary, molec_boundary,                                    &
   tint  , ztodt  , nbot_molec ,                                          &
   lchnk , t, alphath, waccmx_mode, ubc_mmr, ubc_flux) &
       result(du)

   use coords_1d,           only: Coords1D
   use linear_1d_operators, only: BoundaryType,  &
                                  TriDiagDecomp, &
                                  BoundaryData,  &
                                  BoundaryFlux
   use vdiff_lu_solver,     only: fin_vol_lu_decomp
   use shr_kind_mod,        only: r8 => shr_kind_r8
   use air_composition,     only: mbarv
   use physconst,           only: mwdry, gravit

   integer,  intent(in)    :: pcols
   integer,  intent(in)    :: pver
   integer,  intent(in)    :: ncol                  ! Number of atmospheric columns
   integer,  intent(in)    :: u(pcols, pver)
   integer,  intent(in)    :: nbot_molec

   logical,  intent(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
   logical,  intent(in)    :: cnst_fixed_ubflx
   real(r8), intent(in)    :: kv(pcols,pver+1)      ! Eddy diffusivity
   real(r8), intent(in)    :: kq_scal(pcols,pver+1) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
   real(r8), intent(in)    :: mw                    ! Molecular weight for this constituent
   real(r8), intent(in)    :: mw_facm(pcols,pver+1) ! composition dependent sqrt(1/M_q + 1/M_d) for this constituent
   real(r8), intent(in)    :: dpidz_sq(ncol,pver+1) ! (g*rho)**2 (square of vertical derivative of pint)
   type(Coords1D), intent(in) :: p                  ! Pressure coordinates
   type(BoundaryType), intent(in) :: interface_boundary ! Boundary on grid edge.
   type(BoundaryType), intent(in) :: molec_boundary ! Boundary at edge of molec_diff region.
   real(r8), intent(in)    :: tint(pcols,pver+1)    ! Interface temperature [ K ]
   real(r8), intent(in)    :: ztodt                 ! 2 delta-t [ s ]

   integer,  intent(in)    :: lchnk		    ! Chunk number
   real(r8), intent(in)    :: t(pcols,pver)	    ! temperature
   logical,  intent(in), optional :: waccmx_mode  ! running waccmx
   real(r8), intent(in), optional :: alphath
   real(r8), intent(in), optional :: ubc_mmr(pcols)  ! Upper boundary mixing ratios [ kg/kg ]
   real(r8), intent(in), optional :: ubc_flux(pcols)! Upper boundary flux [ kg/s/m^2 ]

   real(r8) :: du(pcols, pver)

   ! LU decomposition information for solver.
   type(TriDiagDecomp) :: decomp

   ! --------------- !
   ! Local Variables !
   ! --------------- !

   ! Level index.
   integer :: k

   ! Molecular diffusivity for constituent.
   real(r8)                :: kmq(ncol,nbot_molec+1)

   ! Term for drift due to molecular separation: (m_i/m - 1) / p
   real(r8)                :: mw_term(ncol,nbot_molec+1)

   ! Diffusion coefficient.
   real(r8)                :: diff_coef(ncol,nbot_molec+1)
   ! Advection velocity.
   real(r8)                :: advect_v(ncol,nbot_molec+1)

   ! 1/mbar * d(mbar)/dp
   real(r8)                :: gradm(ncol,nbot_molec+1)

   ! alphaTh/T * dT/dp, for now alphaTh is non-zero only for H.
   real(r8)                :: gradt(ncol,nbot_molec+1)

   ! mbarv at interface
   real(r8)                :: mbarvi(ncol)

   ! ----------------------- !
   ! Main Computation Begins !
   ! ----------------------- !

   ! --------------------------------------------------------------------- !
   ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
   ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
   ! a combination of ca and cc; they are not required by the solver.      !
   !---------------------------------------------------------------------- !

   kmq  = 0._r8
   mw_term = 0._r8
   gradm = 0._r8
   gradt = 0._r8

   ! Compute difference between scale heights of constituent and dry air

   if ( waccmx_mode ) then

      ! Top level first.
      k = 1
      mbarvi = .75_r8*mbarv(:ncol,k,lchnk)+0.5_r8*mbarv(:ncol,k+1,lchnk) &
           -.25_r8*mbarv(:ncol,k+2,lchnk)
      mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
      gradm(:,k) = (mbarv(:ncol,k,lchnk)-mbarvi)/ &
           (p%mid(:,k)-p%ifc(:,k))/ &
           (mbarv(:ncol,k,lchnk)+mbarvi)*2._r8

      if (alphath /= 0._r8) then
         gradt(:,k) = alphath*(t(:ncol,k)-tint(:ncol,k))/ &
              (p%mid(:ncol,k)-p%ifc(:ncol,k))/ &
              (t(:ncol,k)+tint(:ncol,k))*2._r8
      end if

      ! Interior of molecular diffusion region.
      do k = 2, nbot_molec
         mbarvi = 0.5_r8 * (mbarv(:ncol,k-1,lchnk)+mbarv(:ncol,k,lchnk))
         mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
         gradm(:,k) = (mbarv(:ncol,k,lchnk)-mbarv(:ncol,k-1,lchnk)) * &
              p%rdst(:,k-1)/mbarvi
      enddo

      if (alphath /= 0._r8) then
         do k = 2, nbot_molec
            gradt(:,k) = alphath*(t(:ncol,k)-t(:ncol,k-1)) &
                 *p%rdst(:,k-1)/tint(:ncol,k)
         end do
      end if

      ! Leave nbot_molec+1 terms as zero, because molecular diffusion is
      ! small at the lower boundary.

   else

      do k = 1, nbot_molec
         mw_term(:,k) = (mw/mwdry - 1._r8) / p%ifc(:ncol,k)
      enddo

   endif

   !-------------------- !
   ! Molecular diffusion !
   !-------------------- !

   ! Start with non-molecular portion of diffusion.

   ! Molecular diffusion coefficient.
   do k = 1, nbot_molec
      kmq(:,k)  = kq_scal(:ncol,k) * mw_facm(:ncol,k)
   end do

   diff_coef = kv(:ncol,:nbot_molec+1) + kmq

   ! "Drift" terms.
   advect_v = kmq*mw_term
   if ( waccmx_mode ) then
      advect_v = advect_v - kmq*gradt - &
           (kv(:ncol,:nbot_molec+1) + kmq)*gradm
   end if

   ! Convert from z to pressure representation.
   diff_coef = dpidz_sq(:,:nbot_molec+1) * diff_coef
   advect_v = dpidz_sq(:,:nbot_molec+1) * advect_v

   if( fixed_ubc ) then
      decomp = fin_vol_lu_decomp(ztodt, p, &
           coef_q_diff=diff_coef, coef_q_adv=advect_v, &
           upper_bndry=interface_boundary, &
           lower_bndry=molec_boundary)
   else
      decomp = fin_vol_lu_decomp(ztodt, p, &
           coef_q_diff=diff_coef, coef_q_adv=advect_v, &
           lower_bndry=molec_boundary)
   end if

   du = u
   if (cnst_fixed_ubflx) then
      call decomp%left_div(du(:ncol,:), &
                      l_cond=BoundaryFlux( &
                      -gravit*ubc_flux(:ncol), ztodt, &
                      p%del(:,1)))
   else
      call decomp%left_div(du(:ncol,:), &
            l_cond=BoundaryData(ubc_mmr(:ncol)))
   end if
   du = du - u

end function vd_lu_solve

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

function fin_vol_solve(dt, p, u, ncols, pver, coef_q, coef_q_diff, coef_q_adv, &
   coef_q_weight, upper_bndry, lower_bndry, l_cond, r_cond)  result(du)

   use linear_1d_operators, only: &
     zero_operator,               &
     diagonal_operator,           &
     diffusion_operator,          &
     advection_operator,          &
     BoundaryType,                &
     TriDiagDecomp,               &
     TriDiagOp,                   &
     BoundaryCond,                &
     operator(+)
   use shr_kind_mod,        only: r8 => shr_kind_r8
   use coords_1d,           only: Coords1D

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !

   ! Time step.
   real(r8), intent(in) :: dt
   ! Grid spacings.
   type(Coords1D), intent(in) :: p

   ! Matrix to decomp from.
   ! real(r8), intent(in) :: u(ncols,pver)
   integer,  intent(in)    :: ncols
   integer,  intent(in)    :: pver
   real(r8), intent(in) :: u(ncols,pver)

   ! Coefficients for diffusion and advection.
   !
   ! The sizes must be consistent among all the coefficients that are
   ! actually present, i.e. coef_q_diff and coef_q_adv should be one level
   ! bigger than coef_q and coef_q_weight, and have the same column number.
   real(r8), contiguous, intent(in), optional :: coef_q(:,:), &
      coef_q_diff(:,:), coef_q_adv(:,:), coef_q_weight(:,:)

   ! Boundary conditions (optional, default to 0 flux through boundary).
   class(BoundaryType), target, intent(in), optional :: &
      upper_bndry, lower_bndry

   ! Objects representing boundary conditions.
   class(BoundaryCond), intent(in), optional :: l_cond, r_cond

   real(r8) :: du(ncols,pver)

   ! decomposition.
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

   ! Decompose
   decomp = TriDiagDecomp(net_operator)
   du = u
   call decomp%left_div(du(:ncols, :), l_cond=l_cond)
   du = du - u
   ! Ensure local objects are deallocated.
   call net_operator%finalize()
   call add_term%finalize()
   ! du = u
   ! call decomp%left_div(du(:ncols, :), l_cond, r_cond)
   ! call decomp%left_div(du(:ncols, :), l_cond=l_cond)
   call decomp%finalize()
   ! du = u - du

end function fin_vol_solve

end module new_decomp
