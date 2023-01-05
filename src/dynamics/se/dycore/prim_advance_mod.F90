module prim_advance_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8
  use edgetype_mod,   only: EdgeBuffer_t
  use perf_mod,       only: t_startf, t_stopf, t_adj_detailf !, t_barrierf _EXTERNAL
  use cam_abortutils, only: endrun
  use parallel_mod,   only: parallel_t, HME_BNDRY_P2P!,HME_BNDRY_A2A
  use thread_mod ,    only: horz_num_threads, vert_num_threads, omp_set_nested

  implicit none
  private
  save

  public :: prim_advance_exp, prim_advance_init, applyCAMforcing, calc_tot_energy_dynamics, compute_omega

  type (EdgeBuffer_t) :: edge3,edgeOmega,edgeSponge
  real (kind=r8), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init(par, elem)
    use edge_mod,       only: initEdgeBuffer
    use element_mod,    only: element_t
    use dimensions_mod, only: nlev,ksponge_end
    use control_mod,    only: qsplit

    type (parallel_t)                       :: par
    type (element_t), target, intent(inout) :: elem(:)
    integer                                 :: i

    call initEdgeBuffer(par,edge3   ,elem,4*nlev   ,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
    call initEdgeBuffer(par,edgeSponge,elem,4*ksponge_end,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)
    call initEdgeBuffer(par,edgeOmega ,elem,nlev         ,bndry_type=HME_BNDRY_P2P, nthreads=horz_num_threads)

    if(.not. allocated(ur_weights)) allocate(ur_weights(qsplit))
    ur_weights(:)=0.0_r8

    if(mod(qsplit,2).NE.0)then
      ur_weights(1)=1.0_r8/qsplit
      do i=3,qsplit,2
        ur_weights(i)=2.0_r8/qsplit
      enddo
    else
      do i=2,qsplit,2
        ur_weights(i)=2.0_r8/qsplit
      enddo
    endif
  end subroutine prim_advance_init

  subroutine prim_advance_exp(elem, fvm, deriv, hvcoord, hybrid,dt, tl,  nets, nete)
    use control_mod,       only: tstep_type, qsplit
    use derivative_mod,    only: derivative_t
    use dimensions_mod,    only: np, nlev
    use element_mod,       only: element_t
    use hybvcoord_mod,     only: hvcoord_t
    use hybrid_mod,        only: hybrid_t
    use time_mod,          only: TimeLevel_t,  timelevel_qdp, tevolve
    use dimensions_mod,    only: lcp_moist
    use fvm_control_volume_mod, only: fvm_struct
    use control_mod,       only: raytau0
    use physconst,         only: get_cp, thermodynamic_active_species_num
    use physconst,         only: get_kappa_dry, dry_air_species_num
    use physconst,         only: thermodynamic_active_species_idx_dycore
    use physconst,         only: cpair, rair
    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type(fvm_struct)     , intent(in) :: fvm(:)
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=r8), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! Local
    real (kind=r8) :: dt_vis, eta_ave_w
    real (kind=r8) :: dp(np,np)
    integer        :: ie,nm1,n0,np1,k,qn0,m_cnst, nq
    real (kind=r8) :: inv_cp_full(np,np,nlev,nets:nete)
    real (kind=r8) :: qwater(np,np,nlev,thermodynamic_active_species_num,nets:nete)
    integer        :: qidx(thermodynamic_active_species_num)
    real (kind=r8) :: kappa(np,np,nlev,nets:nete)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1

    call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    !
    !   tstep_type=1  RK2-SSP 3 stage (as used by tracers)           CFL=.58
    !                    optimal in terms of SSP CFL, but not        CFLSSP=2
    !                    optimal in terms of CFL
    !                    typically requires qsplit=3
    !                    but if windspeed > 340m/s, could use this
    !                    with qsplit=1
    !   tstep_type=2  classic RK3                                    CFL=1.73 (sqrt(3))
    !
    !   tstep_type=3  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
    !                 should we replace by standard RK4 (CFL=sqrt(8))?
    !                 (K&G 1st order method has CFL=3)
    !   tstep_type=4  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
    !                 From Paul Ullrich.  3rd order for nonlinear terms also
    !                 K&G method is only 3rd order for linear
    !                 optimal: for windspeeds ~120m/s,gravity: 340m/2
    !                 run with qsplit=1
    !                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
    !

    if (dry_air_species_num > 0) &
      call endrun('ERROR: SE dycore not ready for species dependent thermodynamics - ABORT')

    call omp_set_nested(.true.)

    ! default weights for computing mean dynamics fluxes
    eta_ave_w = 1_r8/qsplit

    ! ==================================
    ! Take timestep
    ! ==================================
    do nq=1,thermodynamic_active_species_num
      qidx(nq) = nq
    end do
    do ie=nets,nete
      do nq=1,thermodynamic_active_species_num
        m_cnst = thermodynamic_active_species_idx_dycore(nq)
        !
        ! make sure Q is updated
        !
        qwater(:,:,:,nq,ie)      = elem(ie)%state%Qdp(:,:,:,m_cnst,qn0)/elem(ie)%state%dp3d(:,:,:,n0)
      end do
    end do
    !
    ! compute Cp and kappa=Rdry/cpdry here and not in RK-stages since Q stays constant => Cp and kappa also stays constant
    !
    if (lcp_moist) then
      do ie=nets,nete
        call get_cp(1,np,1,np,1,nlev,thermodynamic_active_species_num,qwater(:,:,:,:,ie),&
             .true.,inv_cp_full(:,:,:,ie),active_species_idx_dycore=qidx)
      end do
    else
      do ie=nets,nete
        inv_cp_full(:,:,:,ie) = 1.0_r8/cpair
      end do
    end if
    do ie=nets,nete
      call get_kappa_dry(1,np,1,np,1,nlev,nlev,thermodynamic_active_species_num,qwater(:,:,:,:,ie),qidx,kappa(:,:,:,ie))
    end do


    dt_vis = dt

    if (raytau0>0) call rayleigh_friction(elem,n0,nets,nete,dt)
    if (tstep_type==1) then
      ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
      ! not optimal for regular CFL
      ! u1 = u0 + dt/2 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater,qidx,kappa)
      ! u2 = u1 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,np1,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater,qidx,kappa)
      ! u3 = u2 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,np1,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/3,inv_cp_full,qwater,qidx,kappa)

      ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
      do ie=nets,nete
        elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
             + 2*elem(ie)%state%v(:,:,:,:,np1)/3
        elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
             + 2*elem(ie)%state%T(:,:,:,np1)/3
        elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
             + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
      enddo
    else if (tstep_type==2) then
      ! classic RK3  CFL=sqrt(3)
      ! u1 = u0 + dt/3 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! u2 = u0 + dt/2 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! u3 = u0 + dt RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w,inv_cp_full,qwater,qidx,kappa)
    else if (tstep_type==3) then
      ! KG 4th order 4 stage:   CFL=sqrt(8)
      ! low storage version of classic RK4
      ! u1 = u0 + dt/4 RHS(u0)
      call compute_and_apply_rhs(np1,n0,n0,dt/4,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! u2 = u0 + dt/3 RHS(u1)
      call compute_and_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! u3 = u0 + dt/2 RHS(u2)
      call compute_and_apply_rhs(np1,n0,np1,dt/2,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! u4 = u0 + dt RHS(u3)
      call compute_and_apply_rhs(np1,n0,np1,dt,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w,inv_cp_full,qwater,qidx,kappa)
    else if (tstep_type==4) then
      !
      ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
      ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
      ! rhs: t=t
      call compute_and_apply_rhs(nm1,n0,n0,dt/5,elem,hvcoord,hybrid,&
           deriv,nets,nete,eta_ave_w/4,inv_cp_full,qwater,qidx,kappa)
      !
      ! u2 = u0 + dt/5 RHS(u1); rhs: t=t+dt/5
      !
      call compute_and_apply_rhs(np1,n0,nm1,dt/5,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      !
      ! u3 = u0 + dt/3 RHS(u2); rhs: t=t+2*dt/5
      !
      call compute_and_apply_rhs(np1,n0,np1,dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      !
      ! u4 = u0 + 2dt/3 RHS(u3); rhs: t=t+2*dt/5+dt/3
      !
      call compute_and_apply_rhs(np1,n0,np1,2*dt/3,elem,hvcoord,hybrid,&
           deriv,nets,nete,0.0_r8,inv_cp_full,qwater,qidx,kappa)
      ! compute (5*u1/4 - u0/4) in timelevel nm1:
      do ie=nets,nete
        elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
             - elem(ie)%state%v(:,:,:,:,n0) ) /4
        elem(ie)%state%T(:,:,:,nm1)= (5*elem(ie)%state%T(:,:,:,nm1) &
             - elem(ie)%state%T(:,:,:,n0) )/4
        elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
             - elem(ie)%state%dp3d(:,:,:,n0) )/4
      enddo
      ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
      !
      ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
      !
      call compute_and_apply_rhs(np1,nm1,np1,3*dt/4,elem,hvcoord,hybrid,&
           deriv,nets,nete,3*eta_ave_w/4,inv_cp_full,qwater,qidx,kappa)
      ! final method is the same as:
      ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
    else
      call endrun('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================

    call t_startf('advance_hypervis')

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent

    ! forward-in-time, hypervis applied to dp3d
    call advance_hypervis_dp(edge3,elem,fvm,hybrid,deriv,np1,qn0,nets,nete,dt_vis,eta_ave_w,&
         inv_cp_full,hvcoord)

    call t_stopf('advance_hypervis')
    !
    ! update psdry
    !
    do ie=nets,nete
      elem(ie)%state%psdry(:,:) = hvcoord%hyai(1)*hvcoord%ps0
      do k=1,nlev
        elem(ie)%state%psdry(:,:) = elem(ie)%state%psdry(:,:)+elem(ie)%state%dp3d(:,:,k,np1)
      end do
    end do
    tevolve=tevolve+dt

    call omp_set_nested(.false.)

    call t_stopf('prim_advance_exp')
  end subroutine prim_advance_exp


  subroutine applyCAMforcing(elem,fvm,np1,np1_qdp,dt_dribble,dt_phys,nets,nete,nsubstep)
    use dimensions_mod,         only: np, nc, nlev, qsize, ntrac
    use element_mod,            only: element_t
    use control_mod,            only: ftype, ftype_conserve
    use fvm_control_volume_mod, only: fvm_struct
    use physconst,              only: get_dp, thermodynamic_active_species_idx_dycore
    type (element_t)     , intent(inout) :: elem(:)
    type(fvm_struct)     , intent(inout) :: fvm(:)
    real (kind=r8), intent(in) :: dt_dribble, dt_phys
    integer,  intent(in) :: np1,nets,nete,np1_qdp,nsubstep

    ! local
    integer :: i,j,k,ie,q
    real (kind=r8) :: v1,dt_local, dt_local_tracer,tmp
    real (kind=r8) :: dt_local_tracer_fvm
    real (kind=r8) :: ftmp(np,np,nlev,qsize,nets:nete) !diagnostics
    real (kind=r8) :: pdel(np,np,nlev)
    real (kind=r8), allocatable :: ftmp_fvm(:,:,:,:,:) !diagnostics


    if (ntrac>0) allocate(ftmp_fvm(nc,nc,nlev,ntrac,nets:nete))

    if (ftype==0) then
      !
      ! "Dribble" tendencies: divide total adjustment with nsplit and
      !                       add adjustments to state after each
      !                       vertical remap
      !
      dt_local            = dt_dribble
      dt_local_tracer     = dt_dribble
      dt_local_tracer_fvm = dt_dribble
    else if (ftype==1) then
      !
      ! CAM-FV-stype forcing, i.e. equivalent to updating state once in the
      ! beginning of dynamics
      !
      dt_local            = dt_phys
      dt_local_tracer     = dt_phys
      dt_local_tracer_fvm = dt_phys
      if (nsubstep.ne.1) then
        !
        ! do nothing
        !
        dt_local            = 0.0_r8
        dt_local_tracer     = 0.0_r8
        dt_local_tracer_fvm = 0.0_r8
      end if
    else if (ftype==2) then
      !
      ! do state-update for tracers and "dribbling" forcing for u,v,T
      !
      dt_local            = dt_dribble
      if (ntrac>0) then
        dt_local_tracer     = dt_dribble
        dt_local_tracer_fvm = dt_phys
        if (nsubstep.ne.1) then
          dt_local_tracer_fvm = 0.0_r8
        end if
      else
        dt_local_tracer     = dt_phys
        dt_local_tracer_fvm = dt_phys
        if (nsubstep.ne.1) then
          dt_local_tracer     = 0.0_r8
          dt_local_tracer_fvm = 0.0_r8
        end if
      end if
    end if

    do ie=nets,nete
      !
      ! tracers
      !
      if (qsize>0.and.dt_local_tracer>0) then
#if (defined COLUMN_OPENMP)
    !$omp parallel do num_threads(tracer_num_threads) private(q,k,i,j,v1)
#endif
        do q=1,qsize
          do k=1,nlev
            do j=1,np
              do i=1,np
                !
                ! FQ holds q-tendency: (qnew-qold)/dt_physics
                !
                v1 = dt_local_tracer*elem(ie)%derived%FQ(i,j,k,q)
                if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                  if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                    v1=0  ! Q already negative, dont make it more so
                  else
                    v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                  endif
                endif
                elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                ftmp(i,j,k,q,ie) = dt_local_tracer*&
                     elem(ie)%derived%FQ(i,j,k,q)-v1 !Only used for diagnostics!
              enddo
            enddo
          enddo
        enddo
      else
        ftmp(:,:,:,:,ie) = 0.0_r8
      end if
      if (ntrac>0.and.dt_local_tracer_fvm>0) then
        !
        ! Repeat for the fvm tracers: fc holds tendency (fc_new-fc_old)/dt_physics
        !
        do q = 1, ntrac
          do k = 1, nlev
            do j = 1, nc
              do i = 1, nc
                tmp = dt_local_tracer_fvm*fvm(ie)%fc(i,j,k,q)/fvm(ie)%dp_fvm(i,j,k)
                v1 = tmp
                if (fvm(ie)%c(i,j,k,q) + v1 < 0 .and. v1<0) then
                  if (fvm(ie)%c(i,j,k,q) < 0 ) then
                    v1 = 0  ! C already negative, dont make it more so
                  else
                    v1 = -fvm(ie)%c(i,j,k,q)
                  end if
                end if
                fvm(ie)%c(i,j,k,q) = fvm(ie)%c(i,j,k,q)+ v1
                ftmp_fvm(i,j,k,q,ie) = tmp-v1 !Only used for diagnostics!
              end do
            end do
          end do
        end do
      else
        if (ntrac>0) ftmp_fvm(:,:,:,:,ie) = 0.0_r8
      end if


      if (ftype_conserve==1) then
        call get_dp(1,np,1,np,1,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp),2, &
            thermodynamic_active_species_idx_dycore,elem(ie)%state%dp3d(:,:,:,np1),pdel)
        do k=1,nlev
          do j=1,np
            do i = 1,np
              pdel(i,j,k)=elem(ie)%derived%FDP(i,j,k)/pdel(i,j,k)

              elem(ie)%state%T(i,j,k,np1) = elem(ie)%state%T(i,j,k,np1) + &
                   dt_local*elem(ie)%derived%FT(i,j,k)*pdel(i,j,k)
              !
              ! momentum conserving: dp*u
              !
              elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%state%v(i,j,1,k,np1) + &
                   dt_local*elem(ie)%derived%FM(i,j,1,k)*pdel(i,j,k)!elem(ie)%state%dp3d(i,j,k,np1)
              elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%state%v(i,j,2,k,np1) + &
                   dt_local*elem(ie)%derived%FM(i,j,2,k)*pdel(i,j,k)!/elem(ie)%state%dp3d(i,j,k,np1)
            end do
          end do
        end do
      else
        elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
             dt_local*elem(ie)%derived%FT(:,:,:)
        elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
             dt_local*elem(ie)%derived%FM(:,:,:,:)
      end if
    end do
    if (ntrac>0) then
      call output_qdp_var_dynamics(ftmp_fvm(:,:,:,:,:),nc,ntrac,nets,nete,'PDC')
    else
      call output_qdp_var_dynamics(ftmp(:,:,:,:,:),np,qsize,nets,nete,'PDC')
    end if
    if (ftype==1.and.nsubstep==1) call calc_tot_energy_dynamics(elem,fvm,nets,nete,np1,np1_qdp,'p2d')
    if (ntrac>0) deallocate(ftmp_fvm)
  end subroutine applyCAMforcing


  subroutine advance_hypervis_dp(edge3,elem,fvm,hybrid,deriv,nt,qn0,nets,nete,dt2,eta_ave_w,inv_cp_full,hvcoord)
    !
    !  take one timestep of:
    !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
    !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
    !
    !
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    !
    !
    use physconst,      only: gravit, cappa, cpair, tref, lapse_rate, get_dp_ref
    use dimensions_mod, only: np, nlev, nc, ntrac, npsq, qsize
    use dimensions_mod, only: hypervis_dynamic_ref_state,ksponge_end
    use dimensions_mod, only: nu_scale_top,nu_lev,kmvis_ref,kmcnd_ref,rho_ref,km_sponge_factor
    use dimensions_mod, only: kmvisi_ref,kmcndi_ref,rhoi_ref
    use control_mod,    only: nu, nu_s, hypervis_subcycle,hypervis_subcycle_sponge, nu_p, nu_top
    use control_mod,    only: molecular_diff
    use hybrid_mod,     only: hybrid_t!, get_loop_ranges
    use element_mod,    only: element_t
    use derivative_mod, only: derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vlaplace_sphere_wk_mol
    use derivative_mod, only: subcell_Laplace_fluxes, subcell_dss_fluxes
    use edge_mod,       only: edgevpack, edgevunpack, edgeDGVunpack
    use edgetype_mod,   only: EdgeBuffer_t, EdgeDescriptor_t
    use bndry_mod,      only: bndry_exchange
    use viscosity_mod,  only: biharmonic_wk_dp3d
    use hybvcoord_mod,  only: hvcoord_t
    use fvm_control_volume_mod, only: fvm_struct
    use physconst,       only: thermodynamic_active_species_idx_dycore
    use physconst,       only: get_molecular_diff_coef,get_rho_dry
    use cam_history,     only: outfld, hist_fld_active

    type (hybrid_t)    , intent(in)   :: hybrid
    type (element_t)   , intent(inout), target :: elem(:)
    type(fvm_struct)   , intent(in)   :: fvm(:)
    type (EdgeBuffer_t), intent(inout):: edge3
    type (derivative_t), intent(in  ) :: deriv
    integer            , intent(in)   :: nets,nete, nt, qn0
    real (kind=r8)     , intent(in)   :: inv_cp_full(np,np,nlev,nets:nete)
    type (hvcoord_t)   , intent(in)   :: hvcoord
    real (kind=r8) :: eta_ave_w  ! weighting for mean flux terms
    real (kind=r8) :: dt2
    ! local
    integer :: k,kptr,i,j,ie,ic
    integer :: kbeg, kend, kblk
    real (kind=r8), dimension(np,np,2,nlev,nets:nete)      :: vtens
    real (kind=r8), dimension(np,np,nlev,nets:nete)        :: ttens, dptens
    real (kind=r8), dimension(np,np,nlev,nets:nete)        :: dp3d_ref, T_ref
    real (kind=r8), dimension(np,np,nets:nete)             :: ps_ref
    real (kind=r8), dimension(0:np+1,0:np+1,nlev)          :: corners
    real (kind=r8), dimension(2,2,2)                       :: cflux
    real (kind=r8)                                         :: temp      (np,np,nlev)
    real (kind=r8)                                         :: tempflux(nc,nc,4)
    real (kind=r8), dimension(nc,nc,4,nlev,nets:nete)      :: dpflux
    type (EdgeDescriptor_t)                                :: desc

    real (kind=r8), dimension(np,np)            :: lap_t,lap_dp
    real (kind=r8), dimension(np,np)            :: tmp, tmp2
    real (kind=r8), dimension(np,np,ksponge_end,nets:nete):: kmvis,kmcnd,rho_dry
    real (kind=r8), dimension(np,np,ksponge_end+1):: kmvisi,kmcndi
    real (kind=r8), dimension(np,np,ksponge_end+1):: pint,rhoi_dry
    real (kind=r8), dimension(np,np,ksponge_end  ):: pmid
    real (kind=r8), dimension(np,np,nlev)       :: tmp_kmvis,tmp_kmcnd
    real (kind=r8), dimension(np,np,2)          :: lap_v
    real (kind=r8)                              :: v1,v2,v1new,v2new,dt,heating,T0,T1,dp0
    real (kind=r8)                              :: laplace_fluxes(nc,nc,4)
    real (kind=r8)                              :: rhypervis_subcycle
    real (kind=r8)                              :: nu_ratio1, ptop, inv_rho
    real (kind=r8), dimension(ksponge_end)      :: dtemp,du,dv
    real (kind=r8)                              :: nu_temp, nu_dp, nu_velo

    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;

    ptop = hvcoord%hyai(1)*hvcoord%ps0

    if (hypervis_dynamic_ref_state) then
      !
      ! use dynamic reference pressure (P. Callaghan)
      !
      call calc_dp3d_reference(elem,edge3,hybrid,nets,nete,nt,hvcoord,dp3d_ref)
      do ie=nets,nete
        ps_ref(:,:,ie) = ptop + sum(elem(ie)%state%dp3d(:,:,:,nt),3)
      end do
    else
      !
      ! use static reference pressure (hydrostatic balance incl. effect of topography)
      !
      do ie=nets,nete
        call get_dp_ref(hvcoord%hyai, hvcoord%hybi, hvcoord%ps0,1,np,1,np,1,nlev,&
             elem(ie)%state%phis(:,:),dp3d_ref(:,:,:,ie),ps_ref(:,:,ie))
      end do
    endif
    !
    ! reference temperature profile (Simmons and Jiabin, 1991, QJRMS, Section 2a)
    !
    !  Tref = T0+T1*Exner
    !  T1 = .0065*Tref*Cp/g ! = ~191
    !  T0 = Tref-T1         ! = ~97
    !
    T1 = lapse_rate*Tref*cpair/gravit
    T0 = Tref-T1
    do ie=nets,nete
      do k=1,nlev
        dp3d_ref(:,:,k,ie) = ((hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 + &
                              (hvcoord%hybi(k+1)-hvcoord%hybi(k))*ps_ref(:,:,ie))
        dp0                = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
        !
        ! pel@ucar.edu: resolved noise issue over Antarctica
        ! dp3d_ref is the reference-level thickness that includes topography
        ! dp0 is the reference thickness without topography
        !
        dp3d_ref(:,:,k,ie) = dp3d_ref(:,:,k,ie) - dp0
        tmp                = hvcoord%hyam(k)*hvcoord%ps0+hvcoord%hybm(k)*ps_ref(:,:,ie)
        tmp2               = (tmp/hvcoord%ps0)**cappa
        T_ref(:,:,k,ie)    = (T0+T1*tmp2)
      end do
    end do

    kbeg=1; kend=nlev

    kblk = kend - kbeg + 1

    dt=dt2/hypervis_subcycle

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  hyper viscosity
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do ic=1,hypervis_subcycle
      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,'dBH')

      rhypervis_subcycle=1.0_r8/real(hypervis_subcycle,kind=r8)
      call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend,&
           dp3d_ref,T_ref)

      do ie=nets,nete
        ! compute mean flux
        if (nu_p>0) then
          do k=kbeg,kend
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                elem(ie)%derived%dpdiss_ave(i,j,k)=elem(ie)%derived%dpdiss_ave(i,j,k)+&
                     rhypervis_subcycle*eta_ave_w*(elem(ie)%state%dp3d(i,j,k,nt)-dp3d_ref(i,j,k,ie))
                elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+&
                     rhypervis_subcycle*eta_ave_w*dptens(i,j,k,ie)
              enddo
            enddo
          enddo
        endif
        !$omp parallel do num_threads(vert_num_threads) private(lap_t,lap_dp,lap_v,laplace_fluxes,nu_ratio1,i,j,k)
        do k=kbeg,kend
          ! advace in time.
          ! note: DSS commutes with time stepping, so we can time advance and then DSS.
          ! note: weak operators alreayd have mass matrix "included"

          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              ttens(i,j,k,ie)   = -nu_s*ttens(i,j,k,ie)
              dptens(i,j,k,ie)  = -nu_p*dptens(i,j,k,ie)
              vtens(i,j,1,k,ie) = -nu_lev(k)*vtens(i,j,1,k,ie)
              vtens(i,j,2,k,ie) = -nu_lev(k)*vtens(i,j,2,k,ie)
            enddo
          enddo

          if (ntrac>0) then
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,nc
              do i=1,nc
                !
                ! del4 mass flux for CSLAM
                !
                elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - &
                     rhypervis_subcycle*eta_ave_w*nu_p*dpflux(i,j,:,k,ie)
              enddo
            enddo
          endif

          ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              elem(ie)%state%dp3d(i,j,k,nt) = elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%spheremp(i,j)&
                   + dt*dptens(i,j,k,ie)
            enddo
          enddo

        enddo

        kptr = kbeg - 1
        call edgeVpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)

        kptr = kbeg - 1 + nlev
        call edgeVpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)

        kptr = kbeg - 1 + 2*nlev
        call edgeVpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

        kptr = kbeg - 1 + 3*nlev
        call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)
      enddo

      call bndry_exchange(hybrid,edge3,location='advance_hypervis_dp2')

      do ie=nets,nete

        kptr = kbeg - 1
        call edgeVunpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)

        kptr = kbeg - 1 + nlev
        call edgeVunpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)

        kptr = kbeg - 1 + 2*nlev
        call edgeVunpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

        if (ntrac>0) then
          do k=kbeg,kend
            temp(:,:,k) = elem(ie)%state%dp3d(:,:,k,nt) / elem(ie)%spheremp  ! STATE before DSS
            corners(0:np+1,0:np+1,k) = 0.0_r8
            corners(1:np  ,1:np  ,k) = elem(ie)%state%dp3d(1:np,1:np,k,nt) ! fill in interior data of STATE*mass
          enddo
        endif
        kptr = kbeg - 1 + 3*nlev
        call edgeVunpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)

        if (ntrac>0) then
          desc = elem(ie)%desc

          kptr = kbeg - 1 + 3*nlev
          call edgeDGVunpack(edge3,corners(:,:,kbeg:kend),kblk,kptr,ie)
          do k=kbeg,kend
            corners(:,:,k) = corners(:,:,k)/dt !note: array size is 0:np+1
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                temp(i,j,k) =  elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,nt) - temp(i,j,k)
                temp(i,j,k) =  temp(i,j,k)/dt
              enddo
            enddo

            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)

            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)

            call subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux,tempflux)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 rhypervis_subcycle*eta_ave_w*tempflux
          end do
        endif

        ! apply inverse mass matrix, accumulate tendencies
        !$omp parallel do num_threads(vert_num_threads) private(k,i,j)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              vtens(i,j,1,k,ie)=dt*vtens(i,j,1,k,ie)*elem(ie)%rspheremp(i,j)
              vtens(i,j,2,k,ie)=dt*vtens(i,j,2,k,ie)*elem(ie)%rspheremp(i,j)
              ttens(i,j,k,ie)=dt*ttens(i,j,k,ie)*elem(ie)%rspheremp(i,j)
              elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%rspheremp(i,j)
            enddo
          enddo
        enddo

        !$omp parallel do num_threads(vert_num_threads) private(k,i,j)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              ! update v first (gives better results than updating v after heating)
              elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                   vtens(i,j,:,k,ie)
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   +ttens(i,j,k,ie)
            enddo
          enddo
        enddo
      end do

      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,'dCH')
      do ie=nets,nete
        !$omp parallel do num_threads(vert_num_threads), private(k,i,j,v1,v2,heating)
        do k=kbeg,kend
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              v1new=elem(ie)%state%v(i,j,1,k,nt)
              v2new=elem(ie)%state%v(i,j,2,k,nt)
              v1   =elem(ie)%state%v(i,j,1,k,nt)- vtens(i,j,1,k,ie)
              v2   =elem(ie)%state%v(i,j,2,k,nt)- vtens(i,j,2,k,ie)
              heating = 0.5_r8*(v1new*v1new+v2new*v2new-(v1*v1+v2*v2))

              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   -heating*inv_cp_full(i,j,k,ie)
            enddo
          enddo
        enddo
      enddo
      call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,'dAH')
    end do

    !
    !***************************************************************
    !
    ! sponge layer damping
    !
    !***************************************************************
    !
    !
    ! vertical diffusion
    !
    call t_startf('vertical_molec_diff')
    if (molecular_diff>1) then
      do ie=nets,nete
        call get_rho_dry(1,np,1,np,ksponge_end,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,qn0),  &
             elem(ie)%state%T(:,:,:,nt),ptop,elem(ie)%state%dp3d(:,:,:,nt),&
             .true.,rhoi_dry=rhoi_dry(:,:,:),                           &
             active_species_idx_dycore=thermodynamic_active_species_idx_dycore,&
             pint_out=pint,pmid_out=pmid)
        !
        ! constant coefficients
        !
        do k=1,ksponge_end+1
           kmvisi(:,:,k) = kmvisi_ref(k)*rhoi_dry(:,:,k)
           kmcndi(:,:,k) = kmcndi_ref(k)*rhoi_dry(:,:,k)
        end do
        !
        ! do vertical diffusion
        !
        do j=1,np
          do i=1,np
            call solve_diffusion(dt2,np,nlev,i,j,ksponge_end,pmid,pint,kmcndi(:,:,:)/cpair,elem(ie)%state%T(:,:,:,nt),&
                 0,dtemp)
            call solve_diffusion(dt2,np,nlev,i,j,ksponge_end,pmid,pint,kmvisi(:,:,:),elem(ie)%state%v(:,:,1,:,nt),1,du)
            call solve_diffusion(dt2,np,nlev,i,j,ksponge_end,pmid,pint,kmvisi(:,:,:),elem(ie)%state%v(:,:,2,:,nt),1,dv)
            do k=1,ksponge_end
              v1    = elem(ie)%state%v(i,j,1,k,nt)
              v2    = elem(ie)%state%v(i,j,2,k,nt)
              v1new = v1 + du(k)
              v2new = v2 + dv(k)
              !
              ! frictional heating
              !
              heating = 0.5_r8*((v1new*v1new+v2new*v2new) - (v1*v1+v2*v2))
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   -heating*inv_cp_full(i,j,k,ie)+dtemp(k)
              elem(ie)%state%v(i,j,1,k,nt)=v1new
              elem(ie)%state%v(i,j,2,k,nt)=v2new
            end do
          end do
        end do
      end do
    end if
    call t_stopf('vertical_molec_diff')
    call t_startf('sponge_diff')
    !
    ! compute coefficients for horizontal diffusion
    !
    if (molecular_diff>0) then
      do ie=nets,nete
        call get_rho_dry(1,np,1,np,ksponge_end,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,qn0),  &
             elem(ie)%state%T(:,:,:,nt),ptop,elem(ie)%state%dp3d(:,:,:,nt),&
             .true.,rho_dry=rho_dry(:,:,:,ie),                                              &
             active_species_idx_dycore=thermodynamic_active_species_idx_dycore)
      end do

      if (molecular_diff==1) then
        do ie=nets,nete
          !
          ! compute molecular diffusion and thermal conductivity coefficients at mid-levels
          !
          call get_molecular_diff_coef(1,np,1,np,ksponge_end,nlev,&
               elem(ie)%state%T(:,:,:,nt),0,km_sponge_factor(1:ksponge_end),kmvis(:,:,:,ie),kmcnd(:,:,:,ie),qsize,&
               elem(ie)%state%Qdp(:,:,:,1:qsize,qn0),fact=1.0_r8/elem(ie)%state%dp3d(:,:,1:ksponge_end,nt),&
               active_species_idx_dycore=thermodynamic_active_species_idx_dycore)
        end do
      else
        !
        ! constant coefficients
        !
        do ie=nets,nete
          do k=1,ksponge_end
            kmvis  (:,:,k,ie) = kmvis_ref(k)
            kmcnd  (:,:,k,ie) = kmcnd_ref(k)
          end do
        end do
      end if
      !
      ! diagnostics
      !
      if (hist_fld_active('nu_kmvis')) then
        do ie=nets,nete
          tmp_kmvis = 0.0_r8
          do k=1,ksponge_end
            tmp_kmvis(:,:,k) = kmvis(:,:,k,ie)/rho_dry(:,:,k,ie)
          end do
          call outfld('nu_kmvis',RESHAPE(tmp_kmvis(:,:,:), (/npsq,nlev/)), npsq, ie)
        end do
      end if
      if (hist_fld_active('nu_kmcnd')) then
        do ie=nets,nete
          tmp_kmcnd = 0.0_r8
          do k=1,ksponge_end
            tmp_kmcnd(:,:,k) = kmcnd(:,:,k,ie)*inv_cp_full(:,:,k,ie)/rho_dry(:,:,k,ie)
          end do
          call outfld('nu_kmcnd',RESHAPE(tmp_kmcnd(:,:,:), (/npsq,nlev/)), npsq, ie)
        end do
      end if
      if (hist_fld_active('nu_kmcnd_dp')) then
        do ie=nets,nete
          tmp_kmcnd = 0.0_r8
          do k=1,ksponge_end
            tmp_kmcnd(:,:,k) = kmcnd(:,:,k,ie)/(cpair*rho_ref(k))
          end do
          call outfld('nu_kmcnd_dp',RESHAPE(tmp_kmcnd(:,:,:), (/npsq,nlev/)), npsq, ie)
        end do
      end if

      !
      ! scale by reference value
      !
      do ie=nets,nete
        do k=1,ksponge_end
          kmcnd(:,:,k,ie) = kmcnd(:,:,k,ie)/kmcnd_ref(k)
          kmvis(:,:,k,ie) = kmvis(:,:,k,ie)/kmvis_ref(k)
        end do
      end do
    end if
    !
    ! Horizontal Laplacian diffusion
    !
    dt=dt2/hypervis_subcycle_sponge
    call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,'dBS')
    kblk = ksponge_end
    do ic=1,hypervis_subcycle_sponge
      rhypervis_subcycle=1.0_r8/real(hypervis_subcycle_sponge,kind=r8)
      do ie=nets,nete
        do k=1,ksponge_end
          if (nu_top>0.or.molecular_diff>1) then
            !**************************************************************
            !
            ! traditional sponge formulation (constant coefficients)
            !
            !**************************************************************
            call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_t,var_coef=.false.)
            call laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),lap_dp,var_coef=.false.)
            nu_ratio1=1.0_r8
            call vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),.true.,lap_v, var_coef=.false.,&
                 nu_ratio=nu_ratio1)

            nu_dp   = nu_scale_top(k)*nu_top
            nu_temp = nu_scale_top(k)*nu_top
            nu_velo = nu_scale_top(k)*nu_top
            if (molecular_diff>1) then
              nu_dp   = nu_dp   + kmcnd_ref(k)/(cpair*rho_ref(k))
            end if

            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                ttens(i,j,k,ie)   = nu_temp*lap_t(i,j)
                dptens(i,j,k,ie)  = nu_dp  *lap_dp(i,j)
                vtens(i,j,1,k,ie) = nu_velo*lap_v(i,j,1)
                vtens(i,j,2,k,ie) = nu_velo*lap_v(i,j,2)
              enddo
            enddo
          end if
          if (molecular_diff>0) then
            !************************************************************************
            !
            ! sponge formulation using molecular diffusion and thermal conductivity
            !
            !************************************************************************
            call vlaplace_sphere_wk_mol(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),.false.,kmvis(:,:,k,ie),lap_v)
            call laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),lap_t ,var_coef=.false.,mol_nu=kmcnd(:,:,k,ie))

            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                inv_rho = 1.0_r8/rho_dry(i,j,k,ie)
                ttens(i,j,k,ie)   = ttens(i,j,k,ie)  + kmcnd_ref(k)*inv_cp_full(i,j,k,ie)*inv_rho*lap_t(i,j)
                vtens(i,j,1,k,ie) = vtens(i,j,1,k,ie)+ kmvis_ref(k)*inv_rho*lap_v(i,j,1)
                vtens(i,j,2,k,ie) = vtens(i,j,2,k,ie)+ kmvis_ref(k)*inv_rho*lap_v(i,j,2)
              end do
            end do
          end if

          if (ntrac>0.and.nu_dp>0) then
            !
            ! mass flux for CSLAM due to sponge layer diffusion on dp
            !
            call subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc,laplace_fluxes)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 rhypervis_subcycle*eta_ave_w*nu_dp*laplace_fluxes
          endif

          ! NOTE: we will DSS all tendencies, EXCEPT for dp3d, where we DSS the new state
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              elem(ie)%state%dp3d(i,j,k,nt) = elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%spheremp(i,j)&
                   + dt*dptens(i,j,k,ie)
            enddo
          enddo

        enddo


        kptr = 0
        call edgeVpack(edgeSponge,ttens(:,:,1:ksponge_end,ie),kblk,kptr,ie)

        kptr = ksponge_end
        call edgeVpack(edgeSponge,vtens(:,:,1,1:ksponge_end,ie),kblk,kptr,ie)

        kptr = 2*ksponge_end
        call edgeVpack(edgeSponge,vtens(:,:,2,1:ksponge_end,ie),kblk,kptr,ie)

        kptr = 3*ksponge_end
        call edgeVpack(edgeSponge,elem(ie)%state%dp3d(:,:,1:ksponge_end,nt),kblk,kptr,ie)
      enddo

      call bndry_exchange(hybrid,edgeSponge,location='advance_hypervis_sponge')

      do ie=nets,nete

        kptr = 0
        call edgeVunpack(edgeSponge,ttens(:,:,1:ksponge_end,ie),kblk,kptr,ie)

        kptr = ksponge_end
        call edgeVunpack(edgeSponge,vtens(:,:,1,1:ksponge_end,ie),kblk,kptr,ie)

        kptr = 2*ksponge_end
        call edgeVunpack(edgeSponge,vtens(:,:,2,1:ksponge_end,ie),kblk,kptr,ie)

        if (ntrac>0.and.nu_dp>0.0_r8) then
          do k=1,ksponge_end
            temp(:,:,k) = elem(ie)%state%dp3d(:,:,k,nt) / elem(ie)%spheremp  ! STATE before DSS
            corners(0:np+1,0:np+1,k) = 0.0_r8
            corners(1:np  ,1:np  ,k) = elem(ie)%state%dp3d(1:np,1:np,k,nt) ! fill in interior data of STATE*mass
          enddo
        endif
        kptr = 3*ksponge_end
        call edgeVunpack(edgeSponge,elem(ie)%state%dp3d(:,:,1:ksponge_end,nt),kblk,kptr,ie)

        if (ntrac>0.and.nu_dp>0.0_r8) then
          desc = elem(ie)%desc

          kptr = 3*ksponge_end
          call edgeDGVunpack(edgeSponge,corners(:,:,1:ksponge_end),kblk,kptr,ie)
          do k=1,ksponge_end
            corners(:,:,k) = corners(:,:,k)/dt !note: array size is 0:np+1
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                temp(i,j,k) =  elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,nt) - temp(i,j,k)
                temp(i,j,k) =  temp(i,j,k)/dt
              enddo
            enddo

            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)

            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)

            call subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux,tempflux)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 rhypervis_subcycle*eta_ave_w*tempflux
          end do
        endif

        ! apply inverse mass matrix, accumulate tendencies
        !$omp parallel do num_threads(vert_num_threads) private(k,i,j)
        do k=1,ksponge_end
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              vtens(i,j,1,k,ie)=dt*vtens(i,j,1,k,ie)*elem(ie)%rspheremp(i,j)
              vtens(i,j,2,k,ie)=dt*vtens(i,j,2,k,ie)*elem(ie)%rspheremp(i,j)
              ttens(i,j,k,ie)=dt*ttens(i,j,k,ie)*elem(ie)%rspheremp(i,j)
              elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%rspheremp(i,j)
            enddo
          enddo
        enddo
        !$omp parallel do num_threads(vert_num_threads) private(k,i,j,v1,v2,v1new,v2new)
        do k=1,ksponge_end
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              ! update v first (gives better results than updating v after heating)
              elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                   vtens(i,j,:,k,ie)
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   +ttens(i,j,k,ie)

              v1new=elem(ie)%state%v(i,j,1,k,nt)
              v2new=elem(ie)%state%v(i,j,2,k,nt)
              v1   =elem(ie)%state%v(i,j,1,k,nt)- vtens(i,j,1,k,ie)
              v2   =elem(ie)%state%v(i,j,2,k,nt)- vtens(i,j,2,k,ie)
              !
              ! frictional heating
              !
              heating = 0.5_r8*(v1new*v1new+v2new*v2new-(v1*v1+v2*v2))
              elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                   -heating*inv_cp_full(i,j,k,ie)
            enddo
          enddo
        enddo
      end do
    end do
    call t_stopf('sponge_diff')
    call calc_tot_energy_dynamics(elem,fvm,nets,nete,nt,qn0,'dAS')
  end subroutine advance_hypervis_dp



   subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,elem,hvcoord,hybrid,&
        deriv,nets,nete,eta_ave_w,inv_cp_full,qwater,qidx,kappa)
     ! ===================================
     ! compute the RHS, accumulate into u(np1) and apply DSS
     !
     !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
     !
     ! This subroutine is normally called to compute a leapfrog timestep
     ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
     ! accomodated.  For example, setting nm1=np1=n0 this routine will
     ! take a forward euler step, overwriting the input with the output.
     !
     ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
     !
     ! Combining the RHS and DSS pack operation in one routine
     ! allows us to fuse these two loops for more cache reuse
     !
     ! Combining the dt advance and DSS unpack operation in one routine
     ! allows us to fuse these two loops for more cache reuse
     !
     ! ===================================
     use dimensions_mod,  only: np, nc, nlev, ntrac, ksponge_end
     use hybrid_mod,      only: hybrid_t
     use element_mod,     only: element_t
     use derivative_mod,  only: derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
     use derivative_mod,  only: subcell_div_fluxes, subcell_dss_fluxes
     use edge_mod,        only: edgevpack, edgevunpack, edgeDGVunpack
     use edgetype_mod,    only: edgedescriptor_t
     use bndry_mod,       only: bndry_exchange
     use hybvcoord_mod,   only: hvcoord_t
     use physconst,       only: epsilo, get_gz_given_dp_Tv_Rdry
     use physconst,       only: thermodynamic_active_species_num, get_virtual_temp, get_cp_dry
     use physconst,       only: thermodynamic_active_species_idx_dycore,get_R_dry
     use physconst,       only: dry_air_species_num,get_exner
     use time_mod, only : tevolve

     implicit none
     integer,        intent(in) :: np1,nm1,n0,nets,nete
     real (kind=r8), intent(in) :: dt2

     type (hvcoord_t)     , intent(in) :: hvcoord
     type (hybrid_t)      , intent(in) :: hybrid
     type (element_t)     , intent(inout), target :: elem(:)
     type (derivative_t)  , intent(in) :: deriv
     real (kind=r8)       , intent(in) :: inv_cp_full(np,np,nlev,nets:nete)
     real (kind=r8)       , intent(in) :: qwater(np,np,nlev,thermodynamic_active_species_num,nets:nete)
     integer              , intent(in) :: qidx(thermodynamic_active_species_num)
     real (kind=r8)       , intent(in) :: kappa(np,np,nlev,nets:nete)

     real (kind=r8) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

    ! local
     real (kind=r8), dimension(np,np,nlev)                         :: phi
     real (kind=r8), dimension(np,np,nlev)                         :: omega_full
     real (kind=r8), dimension(np,np,nlev)                         :: divdp_dry
     real (kind=r8), dimension(np,np,nlev)                         :: divdp_full
     real (kind=r8), dimension(np,np,2)                            :: vtemp
     real (kind=r8), dimension(np,np,2)                            :: grad_kappa_term
     real (kind=r8), dimension(np,np,2,nlev)                       :: vdp_dry
     real (kind=r8), dimension(np,np,2,nlev)                       :: vdp_full
     real (kind=r8), dimension(np,np,nlev)                         :: vgrad_p_full
     real (kind=r8), dimension(np,np,2     )                       :: v            !
     real (kind=r8), dimension(np,np)                              :: vgrad_T      ! v.grad(T)
     real (kind=r8), dimension(np,np)                              :: Ephi         ! kinetic energy + PHI term
     real (kind=r8), dimension(np,np,2,nlev)                       :: grad_p_full
     real (kind=r8), dimension(np,np,2,nlev)                       :: grad_p_m_pmet! gradient(p - p_met)
     real (kind=r8), dimension(np,np,nlev)                         :: vort         ! vorticity
     real (kind=r8), dimension(np,np,nlev)                         :: p_dry        ! pressure dry
     real (kind=r8), dimension(np,np,nlev)                         :: dp_dry       ! delta pressure dry
     real (kind=r8), dimension(np,np,nlev)                         :: R_dry, cp_dry!
     real (kind=r8), dimension(np,np,nlev)                         :: p_full       ! pressure
     real (kind=r8), dimension(np,np,nlev)                         :: dp_full
     real (kind=r8), dimension(np,np)                              :: exner
     real (kind=r8), dimension(0:np+1,0:np+1,nlev)                 :: corners
     real (kind=r8), dimension(2,2,2)                              :: cflux
     real (kind=r8), dimension(np,np)                              :: suml
     real (kind=r8) :: vtens1(np,np,nlev),vtens2(np,np,nlev),ttens(np,np,nlev)
     real (kind=r8) :: stashdp3d (np,np,nlev),tempdp3d(np,np), tempflux(nc,nc,4)
     real (kind=r8) :: ckk, term, T_v(np,np,nlev)
     real (kind=r8), dimension(np,np,2) :: grad_exner
     real (kind=r8), dimension(np,np)   :: theta_v

     type (EdgeDescriptor_t):: desc

     real (kind=r8) :: sum_water(np,np,nlev), density_inv(np,np)
     real (kind=r8) :: E,v1,v2,glnps1,glnps2
     integer        :: i,j,k,kptr,ie
     real (kind=r8) :: u_m_umet, v_m_vmet, t_m_tmet, ptop

!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
     call t_adj_detailf(+1)
     call t_startf('compute_and_apply_rhs')
     ptop = hvcoord%hyai(1)*hvcoord%ps0
     do ie=nets,nete
       !
       ! compute virtual temperature and sum_water
       !
       call get_virtual_temp(1,np,1,np,1,nlev,thermodynamic_active_species_num,qwater(:,:,:,:,ie),&
            t_v(:,:,:),temp=elem(ie)%state%T(:,:,:,n0),sum_q =sum_water(:,:,:),&
            active_species_idx_dycore=qidx)
       call get_R_dry(1,np,1,np,1,nlev,1,nlev,thermodynamic_active_species_num,&
            qwater(:,:,:,:,ie),qidx,R_dry)
       call get_cp_dry(1,np,1,np,1,nlev,1,nlev,thermodynamic_active_species_num,&
            qwater(:,:,:,:,ie),qidx,cp_dry)

       do k=1,nlev
         dp_dry(:,:,k)  = elem(ie)%state%dp3d(:,:,k,n0)
         dp_full(:,:,k) = sum_water(:,:,k)*dp_dry(:,:,k)
       end do
       call get_gz_given_dp_Tv_Rdry(1,np,1,np,nlev,dp_full,T_v,R_dry,elem(ie)%state%phis,ptop,phi,pmid=p_full)
       do k=1,nlev
         ! vertically lagrangian code: we advect dp3d instead of ps
         ! we also need grad(p) at all levels (not just grad(ps))
         !p(k)= hyam(k)*ps0 + hybm(k)*ps
         !    = .5_r8*(hyai(k+1)+hyai(k))*ps0 + .5_r8*(hybi(k+1)+hybi(k))*ps
         !    = .5_r8*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
         !
         ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
         !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2

         call gradient_sphere(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full(:,:,:,k))
         ! ==============================
         ! compute vgrad_lnps - for omega_full
         ! ==============================
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)
             vgrad_p_full(i,j,k) = (v1*grad_p_full(i,j,1,k) + v2*grad_p_full(i,j,2,k))
             vdp_dry(i,j,1,k) = v1*dp_dry(i,j,k)
             vdp_dry(i,j,2,k) = v2*dp_dry(i,j,k)
             vdp_full(i,j,1,k) = v1*dp_full(i,j,k)
             vdp_full(i,j,2,k) = v2*dp_full(i,j,k)
           end do
         end do
         ! ================================
         ! Accumulate mean Vel_rho flux in vn0
         ! ================================
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             elem(ie)%derived%vn0(i,j,1,k)=elem(ie)%derived%vn0(i,j,1,k)+eta_ave_w*vdp_dry(i,j,1,k)
             elem(ie)%derived%vn0(i,j,2,k)=elem(ie)%derived%vn0(i,j,2,k)+eta_ave_w*vdp_dry(i,j,2,k)
           enddo
         enddo
         !divdp_dry(:,:,k)
         ! =========================================
         !
         ! Compute relative vorticity and divergence
         !
         ! =========================================
         call divergence_sphere(vdp_dry(:,:,:,k),deriv,elem(ie),divdp_dry(:,:,k))
         call divergence_sphere(vdp_full(:,:,:,k),deriv,elem(ie),divdp_full(:,:,k))
         call vorticity_sphere(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie),vort(:,:,k))
       enddo

       ! ====================================================
       ! Compute omega_full
       ! ====================================================
       ckk         = 0.5_r8
       suml(:,:  ) = 0
#if (defined COLUMN_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(k,j,i,ckk,term)
#endif
       do k=1,nlev
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
         do j=1,np   !   Loop inversion (AAM)
           do i=1,np
             term         = -divdp_full(i,j,k)

             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)

             omega_full(i,j,k) = suml(i,j) + ckk*term+vgrad_p_full(i,j,k)
             suml(i,j)    = suml(i,j)   + term
           end do
         end do
       end do
#if (defined COLUMN_OPENMP)
     !$omp parallel do num_threads(vert_num_threads) private(k)
#endif
       do k=1,nlev  !  Loop index added (AAM)
         elem(ie)%derived%omega(:,:,k) = &
              elem(ie)%derived%omega(:,:,k) + eta_ave_w*omega_full(:,:,k)
       enddo
       ! ==============================================
       ! Compute phi + kinetic energy term: 10*nv*nv Flops
       ! ==============================================
#if (defined COLUMN_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(k,i,j,v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2)
#endif
       vertloop: do k=1,nlev
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             E = 0.5_r8*( v1*v1 + v2*v2 )
             Ephi(i,j)=E+phi(i,j,k)
           end do
         end do
         ! ================================================
         ! compute gradp term (ps/p)*(dp/dps)*T
         ! ================================================
         call gradient_sphere(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp)
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
         end do


         ! vtemp = grad ( E + PHI )
         ! vtemp = gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv)
         call gradient_sphere(Ephi(:,:),deriv,elem(ie)%Dinv,vtemp)
         density_inv(:,:) = R_dry(:,:,k)*T_v(:,:,k)/p_full(:,:,k)

         if (dry_air_species_num==0) then        
           exner(:,:)=(p_full(:,:,k)/hvcoord%ps0)**kappa(:,:,k,ie)
           theta_v(:,:)=T_v(:,:,k)/exner(:,:)
           call gradient_sphere(exner(:,:),deriv,elem(ie)%Dinv,grad_exner)

           grad_exner(:,:,1) = cp_dry(:,:,k)*theta_v(:,:)*grad_exner(:,:,1)
           grad_exner(:,:,2) = cp_dry(:,:,k)*theta_v(:,:)*grad_exner(:,:,2)
         else
           exner(:,:)=(p_full(:,:,k)/hvcoord%ps0)**kappa(:,:,k,ie)
           theta_v(:,:)=T_v(:,:,k)/exner(:,:)
           call gradient_sphere(exner(:,:),deriv,elem(ie)%Dinv,grad_exner)

           call gradient_sphere(kappa(:,:,k,ie),deriv,elem(ie)%Dinv,grad_kappa_term)
           suml = exner(:,:)*LOG(p_full(:,:,k)/hvcoord%ps0)
           grad_kappa_term(:,:,1)=-suml*grad_kappa_term(:,:,1)
           grad_kappa_term(:,:,2)=-suml*grad_kappa_term(:,:,2)

           grad_exner(:,:,1) = cp_dry(:,:,k)*theta_v(:,:)*(grad_exner(:,:,1)+grad_kappa_term(:,:,1))
           grad_exner(:,:,2) = cp_dry(:,:,k)*theta_v(:,:)*(grad_exner(:,:,2)+grad_kappa_term(:,:,2))
         end if

         do j=1,np
           do i=1,np
             glnps1 = grad_exner(i,j,1)
             glnps2 = grad_exner(i,j,2)
             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)

             vtens1(i,j,k) =   &
                  + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,1) - glnps1

             vtens2(i,j,k) =   &
                  - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,2) - glnps2
             ttens(i,j,k)  =  - vgrad_T(i,j) + &
                  density_inv(i,j)*omega_full(i,j,k)*inv_cp_full(i,j,k,ie)
           end do
         end do

       end do vertloop

       ! =========================================================
       ! local element timestep, store in np1.
       ! note that we allow np1=n0 or nm1
       ! apply mass matrix
       ! =========================================================
#if (defined COLUMN_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(k)
#endif
       do k=1,nlev
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,1,k,nm1) + dt2*vtens1(i,j,k) )
             elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,2,k,nm1) + dt2*vtens2(i,j,k) )
             elem(ie)%state%T(i,j,k,np1) = elem(ie)%spheremp(i,j)*(elem(ie)%state%T(i,j,k,nm1) + dt2*ttens(i,j,k))
             elem(ie)%state%dp3d(i,j,k,np1) = &
                  elem(ie)%spheremp(i,j) * (elem(ie)%state%dp3d(i,j,k,nm1) - &
                  dt2 * (divdp_dry(i,j,k)))
           enddo
         enddo


         if (ntrac>0.and.eta_ave_w.ne.0._r8) then
           !OMP_COLLAPSE_SIMD
           !DIR_VECTOR_ALIGNED
           do j=1,np
             do i=1,np
               v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp_dry(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp_dry(i,j,2,k)
               v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp_dry(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp_dry(i,j,2,k)
             enddo
           enddo
           call subcell_div_fluxes(v, np, nc, elem(ie)%metdet,tempflux)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - eta_ave_w*tempflux
         end if
       enddo
       ! =========================================================
       !
       ! Pack
       !
       ! =========================================================
       kptr=0
       call edgeVpack(edge3, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)

       kptr=nlev
       call edgeVpack(edge3, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)

       kptr=kptr+2*nlev
       call edgeVpack(edge3, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     end do

     ! =============================================================
     ! Insert communications here: for shared memory, just a single
     ! sync is required
     ! =============================================================
     call bndry_exchange(hybrid,edge3,location='edge3')
     do ie=nets,nete
       ! ===========================================================
       ! Unpack the edges for vgrad_T and v tendencies...
       ! ===========================================================
       kptr=0
       call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie)

       kptr=nlev
       call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

       if (ntrac>0.and.eta_ave_w.ne.0._r8) then
         do k=1,nlev
           stashdp3d(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)/elem(ie)%spheremp(:,:)
         end do
       endif

       corners = 0.0_r8
       corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)
       kptr=kptr+2*nlev
       call edgeVunpack(edge3, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)

       if  (ntrac>0.and.eta_ave_w.ne.0._r8) then
         desc = elem(ie)%desc

         call edgeDGVunpack(edge3, corners, nlev, kptr, ie)

         corners = corners/dt2

         do k=1,nlev
           tempdp3d = elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
           tempdp3d = tempdp3d - stashdp3d(:,:,k)
           tempdp3d = tempdp3d/dt2

           call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)

           cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
           cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
           cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
           cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)

           call subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux,tempflux)
           elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + eta_ave_w*tempflux
         end do
       end if

       ! ====================================================
       ! Scale tendencies by inverse mass matrix
       ! ====================================================

#if (defined COLUMN_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(k)
#endif
       do k=1,nlev
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
         do i=1,np
             elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,np1)
             elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
             elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
         enddo
         enddo
       end do

       ! vertically lagrangian: complete dp3d timestep:
       do k=1,nlev
         elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
       enddo
     end do

     call t_stopf('compute_and_apply_rhs')
     call t_adj_detailf(-1)
   end subroutine compute_and_apply_rhs


   !
   ! corner fluxes for CSLAM
   !
   subroutine distribute_flux_at_corners(cflux, corners, getmapP)
     use dimensions_mod, only : np, max_corner_elem
     use control_mod,    only : swest

     real(r8), intent(out) :: cflux(2,2,2)
     real(r8), intent(in)  :: corners(0:np+1,0:np+1)
     integer,  intent(in)  :: getmapP(:)

     cflux = 0.0_r8
     if (getmapP(swest+0*max_corner_elem) /= -1) then
       cflux(1,1,1) =                (corners(0,1) - corners(1,1))
       cflux(1,1,1) = cflux(1,1,1) + (corners(0,0) - corners(1,1)) / 2.0_r8
       cflux(1,1,1) = cflux(1,1,1) + (corners(0,1) - corners(1,0)) / 2.0_r8

       cflux(1,1,2) =                (corners(1,0) - corners(1,1))
       cflux(1,1,2) = cflux(1,1,2) + (corners(0,0) - corners(1,1)) / 2.0_r8
       cflux(1,1,2) = cflux(1,1,2) + (corners(1,0) - corners(0,1)) / 2.0_r8
     else
       cflux(1,1,1) =                (corners(0,1) - corners(1,1))
       cflux(1,1,2) =                (corners(1,0) - corners(1,1))
     endif

     if (getmapP(swest+1*max_corner_elem) /= -1) then
       cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))
       cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,0) - corners(np,1)) / 2.0_r8
       cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,1) - corners(np,0)) / 2.0_r8

       cflux(2,1,2) =                (corners(np  ,0) - corners(np,  1))
       cflux(2,1,2) = cflux(2,1,2) + (corners(np+1,0) - corners(np,  1)) / 2.0_r8
       cflux(2,1,2) = cflux(2,1,2) + (corners(np  ,0) - corners(np+1,1)) / 2.0_r8
     else
       cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))
       cflux(2,1,2) =                (corners(np  ,0) - corners(np,1))
     endif

     if (getmapP(swest+2*max_corner_elem) /= -1) then
       cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))
       cflux(1,2,1) = cflux(1,2,1) + (corners(0,np+1) - corners(1,np  )) / 2.0_r8
       cflux(1,2,1) = cflux(1,2,1) + (corners(0,np  ) - corners(1,np+1)) / 2.0_r8

       cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))
       cflux(1,2,2) = cflux(1,2,2) + (corners(0,np+1) - corners(1,np  )) / 2.0_r8
       cflux(1,2,2) = cflux(1,2,2) + (corners(1,np+1) - corners(0,np  )) / 2.0_r8
     else
       cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))
       cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))
     endif

     if (getmapP(swest+3*max_corner_elem) /= -1) then
       cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))
       cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np+1) - corners(np,np  )) / 2.0_r8
       cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np  ) - corners(np,np+1)) / 2.0_r8

       cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))
       cflux(2,2,2) = cflux(2,2,2) + (corners(np+1,np+1) - corners(np,np  )) / 2.0_r8
       cflux(2,2,2) = cflux(2,2,2) + (corners(np  ,np+1) - corners(np+1,np)) / 2.0_r8
     else
       cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))
       cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))
     endif
   end subroutine distribute_flux_at_corners

  subroutine calc_tot_energy_dynamics(elem,fvm,nets,nete,tl,tl_qdp,outfld_name_suffix)
    use dimensions_mod,         only: npsq,nlev,np,lcp_moist,nc,ntrac,qsize
    use physconst,              only: gravit, cpair, rearth,omega
    use element_mod,            only: element_t
    use cam_history,            only: outfld, hist_fld_active
    use constituents,           only: cnst_get_ind
    use string_utils,           only: strlist_get_ind
    use hycoef,                 only: hyai, ps0
    use fvm_control_volume_mod, only: fvm_struct
    use physconst,              only: get_dp, get_cp
    use physconst,              only: thermodynamic_active_species_idx_dycore
    use dimensions_mod,         only: cnst_name_gll
    !------------------------------Arguments--------------------------------

    type (element_t) , intent(in) :: elem(:)
    type(fvm_struct) , intent(in) :: fvm(:)
    integer          , intent(in) :: tl, tl_qdp,nets,nete
    character*(*)    , intent(in) :: outfld_name_suffix ! suffix for "outfld" names

    !---------------------------Local storage-------------------------------

    real(kind=r8) :: se(npsq)                          ! Dry Static energy (J/m2)
    real(kind=r8) :: ke(npsq)                          ! kinetic energy    (J/m2)

    real(kind=r8) :: cdp_fvm(nc,nc,nlev)
    real(kind=r8) :: se_tmp
    real(kind=r8) :: ke_tmp
    real(kind=r8) :: ps(np,np)
    real(kind=r8) :: pdel(np,np,nlev)
    !
    ! global axial angular momentum (AAM) can be separated into one part (mr) associatedwith the relative motion
    ! of the atmosphere with respect to the planets surface (also known as wind AAM) and another part (mo)
    ! associated with the angular velocity OMEGA (2*pi/d, where d is the length of the day) of the planet
    ! (also known as mass AAM)
    !
    real(kind=r8) :: mr(npsq)  ! wind AAM
    real(kind=r8) :: mo(npsq)  ! mass AAM
    real(kind=r8) :: mr_cnst, mo_cnst, cos_lat, mr_tmp, mo_tmp
    real(kind=r8) :: cp(np,np,nlev)

    integer :: ie,i,j,k
    integer :: ixwv,ixcldice, ixcldliq, ixtt ! CLDICE, CLDLIQ and test tracer indices
    character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6

    !-----------------------------------------------------------------------

    name_out1 = 'SE_'   //trim(outfld_name_suffix)
    name_out2 = 'KE_'   //trim(outfld_name_suffix)
    name_out3 = 'WV_'   //trim(outfld_name_suffix)
    name_out4 = 'WL_'   //trim(outfld_name_suffix)
    name_out5 = 'WI_'   //trim(outfld_name_suffix)
    name_out6 = 'TT_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4).or.hist_fld_active(name_out5).or.hist_fld_active(name_out6)) then

      if (ntrac>0) then
        ixwv = 1
        call cnst_get_ind('CLDLIQ' , ixcldliq, abort=.false.)
        call cnst_get_ind('CLDICE' , ixcldice, abort=.false.)
      else
        !
        ! when using CSLAM the condensates on the GLL grid may be located in a different index than in physics
        !
        ixwv = -1
        call strlist_get_ind(cnst_name_gll, 'CLDLIQ' , ixcldliq, abort=.false.)
        call strlist_get_ind(cnst_name_gll, 'CLDICE' , ixcldice, abort=.false.)
      end if
      call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)
      !
      ! Compute frozen static energy in 3 parts:  KE, SE, and energy associated with vapor and liquid
      !
      do ie=nets,nete
        se    = 0.0_r8
        ke    = 0.0_r8
        call get_dp(1,np,1,np,1,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,tl_qdp),2,thermodynamic_active_species_idx_dycore,&
             elem(ie)%state%dp3d(:,:,:,tl),pdel,ps=ps,ptop=hyai(1)*ps0)
        call get_cp(1,np,1,np,1,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,tl_qdp),&
             .false.,cp,dp_dry=elem(ie)%state%dp3d(:,:,:,tl),&
             active_species_idx_dycore=thermodynamic_active_species_idx_dycore)
        do k = 1, nlev
          do j=1,np
            do i = 1, np
              !
              ! kinetic energy
              !
              ke_tmp   = 0.5_r8*(elem(ie)%state%v(i,j,1,k,tl)**2+ elem(ie)%state%v(i,j,2,k,tl)**2)*pdel(i,j,k)/gravit
              if (lcp_moist) then
                se_tmp = cp(i,j,k)*elem(ie)%state%T(i,j,k,tl)*pdel(i,j,k)/gravit
              else
                !
                ! using CAM physics definition of internal energy
                !
                se_tmp   = cpair*elem(ie)%state%T(i,j,k,tl)*pdel(i,j,k)/gravit
              end if
              se   (i+(j-1)*np) = se   (i+(j-1)*np) + se_tmp
              ke   (i+(j-1)*np) = ke   (i+(j-1)*np) + ke_tmp
            end do
          end do
        end do

        do j=1,np
          do i = 1, np
            se(i+(j-1)*np) = se(i+(j-1)*np) + elem(ie)%state%phis(i,j)*ps(i,j)/gravit
          end do
        end do
        !
        ! Output energy diagnostics on GLL grid
        !
        call outfld(name_out1  ,se       ,npsq,ie)
        call outfld(name_out2  ,ke       ,npsq,ie)
        !
        ! mass variables are output on CSLAM grid if using CSLAM else GLL grid
        !
        if (ntrac>0) then
          if (ixwv>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixwv)*fvm(ie)%dp_fvm(1:nc,1:nc,:)
            call util_function(cdp_fvm,nc,nlev,name_out3,ie)
          end if
          if (ixcldliq>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixcldliq)*fvm(ie)%dp_fvm(1:nc,1:nc,:)
            call util_function(cdp_fvm,nc,nlev,name_out4,ie)
          end if
          if (ixcldice>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixcldice)*fvm(ie)%dp_fvm(1:nc,1:nc,:)
            call util_function(cdp_fvm,nc,nlev,name_out5,ie)
          end if
          if (ixtt>0) then
            cdp_fvm = fvm(ie)%c(1:nc,1:nc,:,ixtt)*fvm(ie)%dp_fvm(1:nc,1:nc,:)
            call util_function(cdp_fvm,nc,nlev,name_out6,ie)
          end if
        else
          call util_function(elem(ie)%state%qdp(:,:,:,1       ,tl_qdp),np,nlev,name_out3,ie)
          if (ixcldliq>0) call util_function(elem(ie)%state%qdp(:,:,:,ixcldliq,tl_qdp),np,nlev,name_out4,ie)
          if (ixcldice>0) call util_function(elem(ie)%state%qdp(:,:,:,ixcldice,tl_qdp),np,nlev,name_out5,ie)
          if (ixtt>0    ) call util_function(elem(ie)%state%qdp(:,:,:,ixtt    ,tl_qdp),np,nlev,name_out6,ie)
        end if
      end do
    end if
    !
    ! Axial angular momentum diagnostics
    !
    ! Code follows
    !
    ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
    ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
    ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140,
    ! doi:10.1002/2013MS000268
    !
    ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
    ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)
    !
    name_out1 = 'MR_'   //trim(outfld_name_suffix)
    name_out2 = 'MO_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2)) then
      call strlist_get_ind(cnst_name_gll, 'CLDLIQ', ixcldliq, abort=.false.)
      call strlist_get_ind(cnst_name_gll, 'CLDICE', ixcldice, abort=.false.)
      mr_cnst = rearth**3/gravit
      mo_cnst = omega*rearth**4/gravit
      do ie=nets,nete
        mr    = 0.0_r8
        mo    = 0.0_r8
        call get_dp(1,np,1,np,1,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,tl_qdp),2,thermodynamic_active_species_idx_dycore,&
             elem(ie)%state%dp3d(:,:,:,tl),pdel,ps=ps,ptop=hyai(1)*ps0)
        do k = 1, nlev
          do j=1,np
            do i = 1, np
              cos_lat = cos(elem(ie)%spherep(i,j)%lat)
              mr_tmp   = mr_cnst*elem(ie)%state%v(i,j,1,k,tl)*pdel(i,j,k)*cos_lat
              mo_tmp   = mo_cnst*pdel(i,j,k)*cos_lat**2

              mr   (i+(j-1)*np) = mr   (i+(j-1)*np) + mr_tmp
              mo   (i+(j-1)*np) = mo   (i+(j-1)*np) + mo_tmp
            end do
          end do
        end do
        call outfld(name_out1  ,mr       ,npsq,ie)
        call outfld(name_out2  ,mo       ,npsq,ie)
      end do
    end if


  end subroutine calc_tot_energy_dynamics

  subroutine output_qdp_var_dynamics(qdp,nx,num_trac,nets,nete,outfld_name)
    use dimensions_mod, only: nlev,ntrac
    use cam_history   , only: outfld, hist_fld_active
    use constituents  , only: cnst_get_ind
    !------------------------------Arguments--------------------------------

    integer      ,intent(in) :: nx,num_trac,nets,nete
    real(kind=r8) :: qdp(nx,nx,nlev,num_trac,nets:nete)
    character*(*),intent(in) :: outfld_name

    !---------------------------Local storage-------------------------------

    integer :: ie
    integer :: ixcldice, ixcldliq, ixtt
    character(len=16) :: name_out1,name_out2,name_out3,name_out4
    !-----------------------------------------------------------------------

    name_out1 = 'WV_'   //trim(outfld_name)
    name_out2 = 'WI_'   //trim(outfld_name)
    name_out3 = 'WL_'   //trim(outfld_name)
    name_out4 = 'TT_'   //trim(outfld_name)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4)) then

      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
      call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)

      do ie=nets,nete
        call util_function(qdp(:,:,:,1,ie),nx,nlev,name_out1,ie)
        if (ixcldice>0) call util_function(qdp(:,:,:,ixcldice,ie),nx,nlev,name_out2,ie)
        if (ixcldliq>0) call util_function(qdp(:,:,:,ixcldliq,ie),nx,nlev,name_out3,ie)
        if (ixtt>0    ) call util_function(qdp(:,:,:,ixtt    ,ie),nx,nlev,name_out4,ie)
      end do
    end if
  end subroutine output_qdp_var_dynamics

  !
  ! column integrate mass-variable and outfld
  !
  subroutine util_function(f_in,nx,nz,name_out,ie)
    use physconst,   only: gravit
    use cam_history, only: outfld, hist_fld_active
    integer,           intent(in) :: nx,nz,ie
    real(kind=r8),     intent(in) :: f_in(nx,nx,nz)
    character(len=16), intent(in) :: name_out
    real(kind=r8)       :: f_out(nx*nx)
    integer             :: i,j,k
    real(kind=r8)       :: inv_g
    if (hist_fld_active(name_out)) then
      f_out = 0.0_r8
      inv_g = 1.0_r8/gravit
      do k = 1, nz
        do j = 1, nx
          do i = 1, nx
            f_out(i+(j-1)*nx) = f_out(i+(j-1)*nx) + f_in(i,j,k)
          end do
        end do
      end do
      f_out = f_out*inv_g
      call outfld(name_out,f_out,nx*nx,ie)
    end if
  end subroutine util_function

   subroutine compute_omega(hybrid,n0,qn0,elem,deriv,nets,nete,dt,hvcoord)
     use control_mod,    only : nu_p, hypervis_subcycle
     use dimensions_mod, only : np, nlev, qsize
     use hybrid_mod,     only : hybrid_t
     use element_mod,    only : element_t
     use derivative_mod, only : divergence_sphere, derivative_t,gradient_sphere
     use hybvcoord_mod,  only : hvcoord_t
     use edge_mod,       only : edgevpack, edgevunpack
     use bndry_mod,      only : bndry_exchange
     use viscosity_mod,  only: biharmonic_wk_omega
     use physconst,      only: thermodynamic_active_species_num, get_dp
     use physconst,      only: thermodynamic_active_species_idx_dycore
     implicit none
     type (hybrid_t)      , intent(in)            :: hybrid
     type (element_t)     , intent(inout), target :: elem(:)
     type (derivative_t)  , intent(in)            :: deriv
     integer              , intent(in)            :: nets,nete,n0,qn0
     real (kind=r8)       , intent(in)            :: dt
     type (hvcoord_t)     , intent(in)            :: hvcoord

     integer        :: i,j,k,ie,kptr,ic
     real (kind=r8) :: ckk, suml(np,np), v1, v2, term
     real (kind=r8) :: dp_full(np,np,nlev)
     real (kind=r8) :: p_full(np,np,nlev),grad_p_full(np,np,2),vgrad_p_full(np,np,nlev)
     real (kind=r8) :: divdp_full(np,np,nlev),vdp_full(np,np,2)
     real(kind=r8)  :: Otens(np,np  ,nlev,nets:nete), dt_hyper, sum_water(np,np,nlev)

     logical, parameter  :: del4omega = .true.

     do ie=nets,nete
        call get_dp(1,np,1,np,1,nlev,qsize,elem(ie)%state%Qdp(:,:,:,1:qsize,qn0),2,&
           thermodynamic_active_species_idx_dycore,elem(ie)%state%dp3d(:,:,:,n0),dp_full)
        do k=1,nlev
           if (k==1) then
              p_full(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0 + dp_full(:,:,k)/2
           else
              p_full(:,:,k)=p_full(:,:,k-1) + dp_full(:,:,k-1)/2 + dp_full(:,:,k)/2
           endif
           call gradient_sphere(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full (:,:,:))
         do j=1,np
           do i=1,np
             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)
             vdp_full(i,j,1) = dp_full(i,j,k)*v1
             vdp_full(i,j,2) = dp_full(i,j,k)*v2
             vgrad_p_full(i,j,k) = (v1*grad_p_full(i,j,1) + v2*grad_p_full(i,j,2))
           end do
         end do
         call divergence_sphere(vdp_full(:,:,:),deriv,elem(ie),divdp_full(:,:,k))
       end do
       ckk       = 0.5_r8
       suml(:,:  ) = 0
       do k=1,nlev
         do j=1,np   !   Loop inversion (AAM)
           do i=1,np
             term         = -divdp_full(i,j,k)

             v1 = elem(ie)%state%v(i,j,1,k,n0)
             v2 = elem(ie)%state%v(i,j,2,k,n0)

             elem(ie)%derived%omega(i,j,k) = suml(i,j) + ckk*term+vgrad_p_full(i,j,k)

             suml(i,j)    = suml(i,j)   + term
           end do
         end do
       end do
     end do
     do ie=nets,nete
       do k=1,nlev
         elem(ie)%derived%omega(:,:,k) = elem(ie)%spheremp(:,:)*elem(ie)%derived%omega(:,:,k)
       end do
       kptr=0
       call edgeVpack(edgeOmega, elem(ie)%derived%omega(:,:,:),nlev,kptr, ie)
     end do
     call bndry_exchange(hybrid,edgeOmega,location='compute_omega #1')
     do ie=nets,nete
       kptr=0
       call edgeVunpack(edgeOmega, elem(ie)%derived%omega(:,:,:),nlev,kptr, ie)
       do k=1,nlev
         elem(ie)%derived%omega(:,:,k) = elem(ie)%rspheremp(:,:)*elem(ie)%derived%omega(:,:,k)
       end do
     end do

     if  (del4omega) then
       dt_hyper=dt/hypervis_subcycle
       do ic=1,hypervis_subcycle
         do ie=nets,nete
           Otens(:,:,:,ie) = elem(ie)%derived%omega(:,:,:)
         end do
         call biharmonic_wk_omega(elem,Otens,deriv,edgeOmega,hybrid,nets,nete,1,nlev)
         do ie=nets,nete
           do k=1,nlev
             Otens(:,:,k,ie) = -dt_hyper*nu_p*Otens(:,:,k,ie)
           end do
           kptr=0
           call edgeVpack(edgeOmega,Otens(:,:,:,ie) ,nlev,kptr, ie)
         end do
         call bndry_exchange(hybrid,edgeOmega,location='compute_omega #2')
         do ie=nets,nete
           kptr=0
           call edgeVunpack(edgeOmega, Otens(:,:,:,ie),nlev,kptr, ie)
         end do
         do ie=nets,nete
           do k=1,nlev
             elem(ie)%derived%omega(:,:,k) =elem(ie)%derived%omega(:,:,k)+&
                  elem(ie)%rspheremp(:,:)*Otens(:,:,k,ie)
           end do
         end do
       end do
     end if
     !call FreeEdgeBuffer(edgeOmega)
   end subroutine compute_omega


  subroutine calc_dp3d_reference(elem,edge3,hybrid,nets,nete,nt,hvcoord,dp3d_ref)
    !
    ! calc_dp3d_reference: When the del^4 horizontal damping is applied to dp3d
    !                      the values are implicitly affected by natural variations
    !                      due to surface topography.
    !
    !                    To account for these physicaly correct variations, use
    !                    the current state values to compute appropriate
    !                    reference values for the current (lagrangian) ETA-surfaces.
    !                    Damping should then be applied to values relative to
    !                    this reference.
    !=======================================================================
    use hybvcoord_mod  ,only: hvcoord_t
    use physconst      ,only: rair,cappa
    use element_mod,    only: element_t
    use dimensions_mod, only: np,nlev
    use hybrid_mod,     only: hybrid_t
    use edge_mod,       only: edgevpack, edgevunpack
    use bndry_mod,      only: bndry_exchange
    !
    ! Passed variables
    !-------------------
    type(element_t   ),target,intent(inout):: elem(:)
    type(EdgeBuffer_t)       ,intent(inout):: edge3
    type(hybrid_t    )       ,intent(in   ):: hybrid
    integer                  ,intent(in   ):: nets,nete
    integer                  ,intent(in   ):: nt
    type(hvcoord_t   )       ,intent(in   ):: hvcoord
    real(kind=r8)            ,intent(out  ):: dp3d_ref(np,np,nlev,nets:nete)
    !
    ! Local Values
    !--------------
    real(kind=r8):: Phis_avg(np,np,     nets:nete)
    real(kind=r8):: Phi_avg (np,np,nlev,nets:nete)
    real(kind=r8):: RT_avg  (np,np,nlev,nets:nete)
    real(kind=r8):: P_val   (np,np,nlev)
    real(kind=r8):: Ps_val  (np,np)
    real(kind=r8):: Phi_val (np,np,nlev)
    real(kind=r8):: Phi_ival(np,np)
    real(kind=r8):: I_Phi   (np,np,nlev+1)
    real(kind=r8):: Alpha   (np,np,nlev  )
    real(kind=r8):: I_P     (np,np,nlev+1)
    real(kind=r8):: DP_avg  (np,np,nlev)
    real(kind=r8):: P_avg   (np,np,nlev)
    real(kind=r8):: Ps_avg  (np,np)
    real(kind=r8):: Ps_ref  (np,np)
    real(kind=r8):: RT_lapse(np,np)
    real(kind=r8):: dlt_Ps  (np,np)
    real(kind=r8):: dPhi    (np,np,nlev)
    real(kind=r8):: dPhis   (np,np)
    real(kind=r8):: E_Awgt,E_phis,E_phi(nlev),E_T(nlev),Lapse0,Expon0
    integer      :: ie,ii,jj,kk,kptr

    ! Loop over elements
    !--------------------
    do ie=nets,nete

      ! Calculate Pressure values from dp3dp
      !--------------------------------------
      P_val(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0 + elem(ie)%state%dp3d(:,:,1,nt)*0.5_r8
      do kk=2,nlev
        P_val(:,:,kk) =               P_val(:,:,kk-1)           &
                      + elem(ie)%state%dp3d(:,:,kk-1,nt)*0.5_r8 &
                      + elem(ie)%state%dp3d(:,:,kk  ,nt)*0.5_r8
      end do
      Ps_val(:,:) = P_val(:,:,nlev) + elem(ie)%state%dp3d(:,:,nlev,nt)*0.5_r8

      ! Calculate (dry) geopotential values
      !--------------------------------------
      dPhi    (:,:,:)    = 0.5_r8*(rair*elem(ie)%state%T   (:,:,:,nt) &
                                      *elem(ie)%state%dp3d(:,:,:,nt) &
                                                    /P_val(:,:,:)    )
      Phi_val (:,:,nlev) = elem(ie)%state%phis(:,:) + dPhi(:,:,nlev)
      Phi_ival(:,:)      = elem(ie)%state%phis(:,:) + dPhi(:,:,nlev)*2._r8
      do kk=(nlev-1),1,-1
        Phi_val (:,:,kk) = Phi_ival(:,:)    + dPhi(:,:,kk)
        Phi_ival(:,:)    = Phi_val (:,:,kk) + dPhi(:,:,kk)
      end do

      ! Calculate Element averages
      !----------------------------
      E_Awgt   = 0.0_r8
      E_phis   = 0.0_r8
      E_phi(:) = 0._r8
      E_T  (:) = 0._r8
      do jj=1,np
      do ii=1,np
        E_Awgt    = E_Awgt    + elem(ie)%spheremp(ii,jj)
        E_phis    = E_phis    + elem(ie)%spheremp(ii,jj)*elem(ie)%state%phis(ii,jj)
        E_phi (:) = E_phi (:) + elem(ie)%spheremp(ii,jj)*Phi_val(ii,jj,:)
        E_T   (:) = E_T   (:) + elem(ie)%spheremp(ii,jj)*elem(ie)%state%T(ii,jj,:,nt)
      end do
      end do

      Phis_avg(:,:,ie) = E_phis/E_Awgt
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = E_phi(kk)     /E_Awgt
        RT_avg (:,:,kk,ie) = E_T  (kk)*rair/E_Awgt
      end do
    end do ! ie=nets,nete

    ! Boundary Exchange of average values
    !-------------------------------------
    do ie=nets,nete
      Phis_avg(:,:,ie) = elem(ie)%spheremp(:,:)*Phis_avg(:,:,ie)
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = elem(ie)%spheremp(:,:)*Phi_avg(:,:,kk,ie)
        RT_avg (:,:,kk,ie) = elem(ie)%spheremp(:,:)*RT_avg (:,:,kk,ie)
      end do
      kptr = 0
      call edgeVpack(edge3,Phi_avg(:,:,:,ie),nlev,kptr,ie)
      kptr = nlev
      call edgeVpack(edge3,RT_avg (:,:,:,ie),nlev,kptr,ie)
      kptr = 2*nlev
      call edgeVpack(edge3,Phis_avg (:,:,ie),1   ,kptr,ie)
    end do ! ie=nets,nete

    call bndry_exchange(hybrid,edge3,location='calc_dp3d_reference')

    do ie=nets,nete
      kptr = 0
      call edgeVunpack(edge3,Phi_avg(:,:,:,ie),nlev,kptr,ie)
      kptr = nlev
      call edgeVunpack(edge3,RT_avg (:,:,:,ie),nlev,kptr,ie)
      kptr = 2*nlev
      call edgeVunpack(edge3,Phis_avg (:,:,ie),1   ,kptr,ie)
      Phis_avg(:,:,ie) = elem(ie)%rspheremp(:,:)*Phis_avg(:,:,ie)
      do kk=1,nlev
        Phi_avg(:,:,kk,ie) = elem(ie)%rspheremp(:,:)*Phi_avg(:,:,kk,ie)
        RT_avg (:,:,kk,ie) = elem(ie)%rspheremp(:,:)*RT_avg (:,:,kk,ie)
      end do
    end do ! ie=nets,nete

    ! Loop over elements
    !--------------------
    do ie=nets,nete

      ! Fill elements with uniformly varying average values
      !-----------------------------------------------------
      call fill_element(Phis_avg(1,1,ie))
      do kk=1,nlev
        call fill_element(Phi_avg(1,1,kk,ie))
        call fill_element(RT_avg (1,1,kk,ie))
      end do

      ! Integrate upward to compute Alpha == (dp3d/P)
      !----------------------------------------------
      I_Phi(:,:,nlev+1) = Phis_avg(:,:,ie)
      do kk=nlev,1,-1
        I_Phi(:,:,kk) = 2._r8* Phi_avg(:,:,kk,ie) - I_Phi(:,:,kk+1)
        Alpha(:,:,kk) = 2._r8*(Phi_avg(:,:,kk,ie) - I_Phi(:,:,kk+1))/RT_avg(:,:,kk,ie)
      end do

      ! Integrate downward to compute corresponding average pressure values
      !---------------------------------------------------------------------
      I_P(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0
      do kk=1,nlev
        DP_avg(:,:,kk  ) = I_P(:,:,kk)*(2._r8 * Alpha(:,:,kk))/(2._r8 - Alpha(:,:,kk))
        P_avg (:,:,kk  ) = I_P(:,:,kk)*(2._r8                )/(2._r8 - Alpha(:,:,kk))
        I_P   (:,:,kk+1) = I_P(:,:,kk)*(2._r8 + Alpha(:,:,kk))/(2._r8 - Alpha(:,:,kk))
      end do
      Ps_avg(:,:) = I_P(:,:,nlev+1)

      ! Determine an appropriate d<T>/d<PHI> lapse rate near the surface
      ! OPTIONALLY: Use dry adiabatic lapse rate or environmental lapse rate.
      !-----------------------------------------------------------------------
      if(.FALSE.) then
        ! DRY ADIABATIC laspe rate
        !------------------------------
        RT_lapse(:,:) = -cappa
      else
        ! ENVIRONMENTAL (empirical) laspe rate
        !--------------------------------------
        RT_lapse(:,:) =  (RT_avg (:,:,nlev-1,ie)-RT_avg (:,:,nlev,ie)) &
                        /(Phi_avg(:,:,nlev-1,ie)-Phi_avg(:,:,nlev,ie))
      endif

      ! Calcualte reference surface pressure
      !--------------------------------------
      dPhis(:,:) = elem(ie)%state%phis(:,:)-Phis_avg(:,:,ie)
      do jj=1,np
      do ii=1,np
        if (abs(RT_lapse(ii,jj)) .gt. 1.e-3_r8) then
          Lapse0 = RT_lapse(ii,jj)/RT_avg(ii,jj,nlev,ie)
          Expon0 = (-1._r8/RT_lapse(ii,jj))
          Ps_ref(ii,jj) = Ps_avg(ii,jj)*((1._r8 + Lapse0*dPhis(ii,jj))**Expon0)
        else
          Ps_ref(ii,jj) = Ps_avg(ii,jj)*exp(-dPhis(ii,jj)/RT_avg(ii,jj,nlev,ie))
        endif
      end do
      end do

      ! Calculate reference dp3d values
      !---------------------------------
      dlt_Ps(:,:) = Ps_ref(:,:) - Ps_avg(:,:)
      do kk=1,nlev
        dp3d_ref(:,:,kk,ie) = DP_avg(:,:,kk) + (hvcoord%hybi(kk+1)            &
                                               -hvcoord%hybi(kk  ))*dlt_Ps(:,:)
      end do

    end do ! ie=nets,nete

    ! End Routine
    !------------
    return
  end subroutine calc_dp3d_reference
  !=============================================================================


  !=============================================================================
  subroutine fill_element(Eval)
    !
    ! fill_element_bilin: Fill in element gridpoints using local bi-linear
    !                     interpolation of nearby average values.
    !
    !                     NOTE: This routine is hard coded for NP=4, if a
    !                           different value of NP is used... bad things
    !                           will happen.
    !=======================================================================
    use dimensions_mod,only: np
    !
    ! Passed variables
    !-------------------
    real(kind=r8),intent(inout):: Eval(np,np)
    !
    ! Local Values
    !--------------
    real(kind=r8):: X0
    real(kind=r8):: S1,S2,S3,S4
    real(kind=r8):: C1,C2,C3,C4
    real(kind=r8):: E1,E2,E3,E4,E0

    X0 = sqrt(1._r8/5._r8)

    ! Set the "known" values Eval
    !----------------------------
    S1 = (Eval(1 ,2 )+Eval(1 ,3 ))/2._r8
    S2 = (Eval(2 ,np)+Eval(3 ,np))/2._r8
    S3 = (Eval(np,2 )+Eval(np,3 ))/2._r8
    S4 = (Eval(2 ,1 )+Eval(3 ,1 ))/2._r8
    C1 = Eval(1 ,1 )
    C2 = Eval(1 ,np)
    C3 = Eval(np,np)
    C4 = Eval(np,1 )

    ! E0 OPTION: Element Center value:
    !---------------------------------
    IF(.FALSE.) THEN
      ! Use ELEMENT AVERAGE value contained in (2,2)
      !----------------------------------------------
      E0 = Eval(2,2)
    ELSE
      ! Use AVG OF SIDE VALUES after boundary exchange of E0 (smooting option)
      !-----------------------------------------------------------------------
      E0 = (S1 + S2 + S3 + S4)/4._r8
    ENDIF

    ! Calc interior values along center axes
    !----------------------------------------
    E1 = E0 + X0*(S1-E0)
    E2 = E0 + X0*(S2-E0)
    E3 = E0 + X0*(S3-E0)
    E4 = E0 + X0*(S4-E0)

    ! Calculate Side Gridpoint Values for Eval
    !------------------------------------------
    Eval(1 ,2 ) = S1 + X0*(C1-S1)
    Eval(1 ,3 ) = S1 + X0*(C2-S1)
    Eval(2 ,np) = S2 + X0*(C2-S2)
    Eval(3 ,np) = S2 + X0*(C3-S2)
    Eval(np,2 ) = S3 + X0*(C4-S3)
    Eval(np,3 ) = S3 + X0*(C3-S3)
    Eval(2 ,1 ) = S4 + X0*(C1-S4)
    Eval(3 ,1 ) = S4 + X0*(C4-S4)

    ! Calculate interior values
    !---------------------------
    Eval(2 ,2 ) = E1 + X0*(Eval(2 ,1 )-E1)
    Eval(2 ,3 ) = E1 + X0*(Eval(2 ,np)-E1)
    Eval(3 ,2 ) = E3 + X0*(Eval(3 ,1 )-E3)
    Eval(3 ,3 ) = E3 + X0*(Eval(3 ,np)-E3)

    ! End Routine
    !------------
    return
  end subroutine fill_element

  subroutine rayleigh_friction(elem,nt,nets,nete,dt)
    use dimensions_mod, only: nlev, otau
    use hybrid_mod,     only: hybrid_t!, get_loop_ranges
    use element_mod,    only: element_t

    type (element_t)   , intent(inout), target :: elem(:)
    integer            , intent(in)   :: nets,nete, nt
    real(r8)                          :: dt

    real(r8) :: c1, c2
    integer  :: k,ie

    do ie=nets,nete
      do k=1,nlev
        c2 = 1._r8 / (1._r8 + otau(k)*dt)
        c1 = -otau(k) * c2 * dt
        elem(ie)%state%v(:,:,:,k,nt) = elem(ie)%state%v(:,:,:,k,nt)+c1 * elem(ie)%state%v(:,:,:,k,nt)
!         ptend%s(:ncol,k) = c3 * (state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
      enddo
    end do
  end subroutine rayleigh_friction



  subroutine solve_diffusion(dt,nx,nlev,i,j,nlay,pmid,pint,km,fld,boundary_condition,dfld)
    use physconst,      only: gravit
    real(kind=r8), intent(in)    :: dt
    integer      , intent(in)    :: nlay, nlev,nx, i, j
    real(kind=r8), intent(in)    :: pmid(nx,nx,nlay),pint(nx,nx,nlay+1),km(nx,nx,nlay+1)
    real(kind=r8), intent(in)    :: fld(nx,nx,nlev)
    real(kind=r8), intent(out)   :: dfld(nlay)
    integer :: boundary_condition
    !
    real(kind=r8), dimension(nlay) :: current_guess,next_iterate
    real(kind=r8)                  :: alp, alm, value_level0
    integer                        :: k,iter, niterations=4

    ! Make the guess for the next time step equal to the initial value
    current_guess(:)= fld(i,j,1:nlay)
    do iter = 1, niterations
      ! two formulations of the upper boundary condition
      !next_iterate(1) = (initial_value(1) + alp * current_guess(i+1) + alm * current_guess(1)) /(1. + alp + alm) ! top BC, u'=0
      if (boundary_condition==0) then
        next_iterate(1) = fld(i,j,1) ! u doesn't get prognosed by diffusion at top
      else if (boundary_condition==1) then
        value_level0 = 0.75_r8*fld(i,j,1) ! value above sponge
        k=1
        alp = dt*(km(i,j,k+1)*gravit*gravit/(pmid(i,j,k)-pmid(i,j,k+1)))/(pint(i,j,k)-pint(i,j,k+1))
        alm = dt*(km(i,j,k  )*gravit*gravit/(0.5_r8*(pmid(i,j,1)-pmid(i,j,2))))/(pint(i,j,k)-pint(i,j,k+1))
        next_iterate(k) = (fld(i,j,k) + alp * current_guess(k+1) + alm * value_level0)/(1._r8 + alp + alm)
      else
        !
        ! set fld'=0 at model top
        !
        k=1
        alp = dt*(km(i,j,k+1)*gravit*gravit/(pmid(i,j,k)-pmid(i,j,k+1)))/(pint(i,j,k)-pint(i,j,k+1))
        alm = dt*(km(i,j,k  )*gravit*gravit/(0.5_r8*(pmid(i,j,1)-pmid(i,j,2))))/(pint(i,j,k)-pint(i,j,k+1))
        next_iterate(k) = (fld(i,j,1) + alp * current_guess(2) + alm * current_guess(1))/(1._r8 + alp + alm)
      end if
      do k = 2, nlay-1
        alp = dt*(km(i,j,k+1)*gravit*gravit/(pmid(i,j,k  )-pmid(i,j,k+1)))/(pint(i,j,k)-pint(i,j,k+1))
        alm = dt*(km(i,j,k  )*gravit*gravit/(pmid(i,j,k-1)-pmid(i,j,k  )))/(pint(i,j,k)-pint(i,j,k+1))
        next_iterate(k) = (fld(i,j,k) + alp * current_guess(k+1) + alm * current_guess(k-1))/(1._r8 + alp + alm)
      end do
      next_iterate(nlay) = (fld(i,j,nlay) + alp * fld(i,j,nlay) + alm * current_guess(nlay-1))/(1._r8 + alp + alm) ! bottom BC

      ! before the next iterate, make the current guess equal to the values of the last iteration
      current_guess(:) = next_iterate(:)
    end do
    dfld(:) = next_iterate(:) - fld(i,j,1:nlay)

  end subroutine solve_diffusion


end module prim_advance_mod
