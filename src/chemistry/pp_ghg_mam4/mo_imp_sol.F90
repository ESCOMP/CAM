module mo_imp_sol
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : clscnt4, gas_pcnst, clsmap
  use cam_logfile, only : iulog
  implicit none
  private
  public :: imp_slv_inti, imp_sol
  save
  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter :: itermax = 11
  integer, parameter :: cut_limit = 5
  real(r8), parameter :: small = 1.e-40_r8
  real(r8) :: epsilon(clscnt4)
  logical :: factor(itermax)
contains
  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls, only : get_spc_ndx
    implicit none
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: m, ox_ndx, o3a_ndx
    real(r8) :: eps(gas_pcnst)
    factor(:) = .true.
    eps(:) = rel_err
    ox_ndx = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    m = get_spc_ndx( 'NO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'N2O5' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'OH' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    o3a_ndx = get_spc_ndx( 'O3A' )
    if( o3a_ndx > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do
  end subroutine imp_slv_inti
  subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
                      ncol,nlev, lchnk, prod_out, loss_out )
    !-----------------------------------------------------------------------
    ! ... imp_sol advances the volumetric mixing ratio
    ! forward one time step via the fully implicit euler scheme.
    ! this source is meant for small l1 cache machines such as
    ! the intel pentium and itanium cpus
    !-----------------------------------------------------------------------
    use chem_mods, only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname, only : solsym
    use mo_lin_matrix, only : linmat
    use mo_nln_matrix, only : nlnmat
    use mo_lu_factor, only : lu_fac
    use mo_lu_solve, only : lu_slv
    use mo_prod_loss, only : imp_prod_loss
    use mo_indprd, only : indprd
    use time_manager, only : get_nstep
    use perf_mod, only : t_startf, t_stopf
    implicit none
    !-----------------------------------------------------------------------
    ! ... dummy args
    !-----------------------------------------------------------------------
    integer, intent(in) :: ncol ! columns in chunck
    integer, intent(in) :: nlev
    integer, intent(in) :: lchnk ! chunk id
    real(r8), intent(in) :: delt ! time step (s)
    real(r8), intent(in) :: reaction_rates(ncol,nlev,max(1,rxntot)) ! rxt rates (1/cm^3/s)
    real(r8), intent(in) :: extfrc(ncol,nlev,max(1,extcnt)) ! external in-situ forcing (1/cm^3/s)
    real(r8), intent(in) :: het_rates(ncol,nlev,max(1,gas_pcnst)) ! washout rates (1/s)
    real(r8), intent(inout) :: base_sol(ncol,nlev,gas_pcnst) ! species mixing ratios (vmr)
    real(r8), intent(out) :: prod_out(ncol,nlev,max(1,clscnt4))
    real(r8), intent(out) :: loss_out(ncol,nlev,max(1,clscnt4))
    !-----------------------------------------------------------------------
    ! ... local variables
    !-----------------------------------------------------------------------
    integer :: nr_iter, &
         lev, &
         i, &
         j, &
         k, l, &
         m
    integer :: fail_cnt, cut_cnt, stp_con_cnt
    integer :: nstep
    real(r8) :: interval_done, dt, dti
    real(r8) :: max_delta(max(1,clscnt4))
    real(r8) :: sys_jac(max(1,nzcnt))
    real(r8) :: lin_jac(max(1,nzcnt))
    real(r8), dimension(max(1,clscnt4)) :: &
         solution, &
         forcing, &
         iter_invariant, &
         prod, &
         loss
    real(r8) :: lrxt(max(1,rxntot))
    real(r8) :: lsol(max(1,gas_pcnst))
    real(r8) :: lhet(max(1,gas_pcnst))
    real(r8), dimension(ncol,nlev,max(1,clscnt4)) :: &
         ind_prd
    logical :: convergence
    logical :: frc_mask, iter_conv
    logical :: converged(max(1,clscnt4))
    solution(:) = 0._r8
    !-----------------------------------------------------------------------
    ! ... class independent forcing
    !-----------------------------------------------------------------------
    if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
       call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
            reaction_rates, ncol )
    else
       do m = 1,max(1,clscnt4)
          ind_prd(:,:,m) = 0._r8
       end do
    end if
    level_loop : do lev = 1,nlev
       column_loop : do i = 1,ncol
          !-----------------------------------------------------------------------
          ! ... transfer from base to local work arrays
          !-----------------------------------------------------------------------
          do m = 1,rxntot
             lrxt(m) = reaction_rates(i,lev,m)
          end do
          if( gas_pcnst > 0 ) then
             do m = 1,gas_pcnst
                lhet(m) = het_rates(i,lev,m)
             end do
          end if
          !-----------------------------------------------------------------------
          ! ... time step loop
          !-----------------------------------------------------------------------
          dt = delt
          cut_cnt = 0
          fail_cnt = 0
          stp_con_cnt = 0
          interval_done = 0._r8
          time_step_loop : do
             dti = 1._r8 / dt
             !-----------------------------------------------------------------------
             ! ... transfer from base to local work arrays
             !-----------------------------------------------------------------------
             do m = 1,gas_pcnst
                lsol(m) = base_sol(i,lev,m)
             end do
             !-----------------------------------------------------------------------
             ! ... transfer from base to class array
             !-----------------------------------------------------------------------
             do k = 1,clscnt4
                j = clsmap(k,4)
                m = permute(k,4)
                solution(m) = lsol(j)
             end do
             !-----------------------------------------------------------------------
             ! ... set the iteration invariant part of the function f(y)
             !-----------------------------------------------------------------------
             if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
                do m = 1,clscnt4
                   iter_invariant(m) = dti * solution(m) + ind_prd(i,lev,m)
                end do
             else
                do m = 1,clscnt4
                   iter_invariant(m) = dti * solution(m)
                end do
             end if
             !-----------------------------------------------------------------------
             ! ... the linear component
             !-----------------------------------------------------------------------
             if( cls_rxt_cnt(2,4) > 0 ) then
                call t_startf( 'lin_mat' )
                call linmat( lin_jac, lsol, lrxt, lhet )
                call t_stopf( 'lin_mat' )
             end if
             !=======================================================================
             ! the newton-raphson iteration for f(y) = 0
             !=======================================================================
             iter_loop : do nr_iter = 1,itermax
                !-----------------------------------------------------------------------
                ! ... the non-linear component
                !-----------------------------------------------------------------------
                if( factor(nr_iter) ) then
                   call t_startf( 'nln_mat' )
                   call nlnmat( sys_jac, lsol, lrxt, lin_jac, dti )
                   call t_stopf( 'nln_mat' )
                   !-----------------------------------------------------------------------
                   ! ... factor the "system" matrix
                   !-----------------------------------------------------------------------
                   call t_startf( 'lu_fac' )
                   call lu_fac( sys_jac )
                   call t_stopf( 'lu_fac' )
                end if
                !-----------------------------------------------------------------------
                ! ... form f(y)
                !-----------------------------------------------------------------------
                call t_startf( 'prod_loss' )
                call imp_prod_loss( prod, loss, lsol, lrxt, lhet )
                call t_stopf( 'prod_loss' )
                do m = 1,clscnt4
                   forcing(m) = solution(m)*dti - (iter_invariant(m) + prod(m) - loss(m))
                end do
                !-----------------------------------------------------------------------
                ! ... solve for the mixing ratio at t(n+1)
                !-----------------------------------------------------------------------
                call t_startf( 'lu_slv' )
                call lu_slv( sys_jac, forcing )
                call t_stopf( 'lu_slv' )
                do m = 1,clscnt4
                   solution(m) = solution(m) + forcing(m)
                end do
                !-----------------------------------------------------------------------
                ! ... convergence measures
                !-----------------------------------------------------------------------
                if( nr_iter > 1 ) then
                   do k = 1,clscnt4
                      m = permute(k,4)
                      if( abs(solution(m)) > 1.e-20_r8 ) then
                         max_delta(k) = abs( forcing(m)/solution(m) )
                      else
                         max_delta(k) = 0._r8
                      end if
                   end do
                end if
                !-----------------------------------------------------------------------
                ! ... limit iterate
                !-----------------------------------------------------------------------
                where( solution(:) < 0._r8 )
                   solution(:) = 0._r8
                endwhere
                !-----------------------------------------------------------------------
                ! ... transfer latest solution back to work array
                !-----------------------------------------------------------------------
                do k = 1,clscnt4
                   j = clsmap(k,4)
                   m = permute(k,4)
                   lsol(j) = solution(m)
                end do
                !-----------------------------------------------------------------------
                ! ... check for convergence
                !-----------------------------------------------------------------------
                converged(:) = .true.
                if( nr_iter > 1 ) then
                   do k = 1,clscnt4
                      m = permute(k,4)
                      frc_mask = abs( forcing(m) ) > small
                      if( frc_mask ) then
                         converged(k) = abs(forcing(m)) <= epsilon(k)*abs(solution(m))
                      else
                         converged(k) = .true.
                      end if
                   end do
                   convergence = all( converged(:) )
                   if( convergence ) then
                      exit
                   end if
                end if
             end do iter_loop
             !-----------------------------------------------------------------------
             ! ... check for newton-raphson convergence
             !-----------------------------------------------------------------------
             if( .not. convergence ) then
                !-----------------------------------------------------------------------
                ! ... non-convergence
                !-----------------------------------------------------------------------
                fail_cnt = fail_cnt + 1
                nstep = get_nstep()
                write(iulog,'('' imp_sol: Time step '',1p,e21.13,'' failed to converge @ (lchnk,lev,col,nstep) = '',4i6)') &
                     dt,lchnk,lev,i,nstep
                stp_con_cnt = 0
                if( cut_cnt < cut_limit ) then
                   cut_cnt = cut_cnt + 1
                   if( cut_cnt < cut_limit ) then
                      dt = .5_r8 * dt
                   else
                      dt = .1_r8 * dt
                   end if
                   cycle time_step_loop
                else
                   write(iulog,'('' imp_sol: Failed to converge @ (lchnk,lev,col,nstep,dt,time) = '',4i6,1p,2e21.13)') &
                        lchnk,lev,i,nstep,dt,interval_done+dt
                   do m = 1,clscnt4
                      if( .not. converged(m) ) then
                         write(iulog,'(1x,a8,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
                      end if
                   end do
                end if
             end if
             !-----------------------------------------------------------------------
             ! ... check for interval done
             !-----------------------------------------------------------------------
             interval_done = interval_done + dt
             if( abs( delt - interval_done ) <= .0001_r8 ) then
                if( fail_cnt > 0 ) then
                   write(iulog,*) 'imp_sol : @ (lchnk,lev,col) = ',lchnk,lev,i,' failed ',fail_cnt,' times'
                end if
                exit time_step_loop
             else
                !-----------------------------------------------------------------------
                ! ... transfer latest solution back to base array
                !-----------------------------------------------------------------------
                if( convergence ) then
                   stp_con_cnt = stp_con_cnt + 1
                end if
                do m = 1,gas_pcnst
                   base_sol(i,lev,m) = lsol(m)
                end do
                if( stp_con_cnt >= 2 ) then
                   dt = 2._r8*dt
                   stp_con_cnt = 0
                end if
                dt = min( dt,delt-interval_done )
                ! write(iulog,'('' imp_sol: New time step '',1p,e21.13)') dt
             end if
          end do time_step_loop
          !-----------------------------------------------------------------------
          ! ... Transfer latest solution back to base array
          !-----------------------------------------------------------------------
          cls_loop: do k = 1,clscnt4
             j = clsmap(k,4)
             m = permute(k,4)
             base_sol(i,lev,j) = solution(m)
             ! output diagnostics
             prod_out(i,lev,k) = prod(k) + ind_prd(i,lev,k)
             loss_out(i,lev,k) = loss(k)
          end do cls_loop
       end do column_loop
    end do level_loop
  end subroutine imp_sol
end module mo_imp_sol
