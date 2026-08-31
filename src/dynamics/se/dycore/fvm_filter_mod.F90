#define FVM_TIMERS .FALSE.
module fvm_filter_mod
  !
  ! CSLAM-grid del4 tracer filter: an optional, mass-conservative fourth-order
  ! hyperdiffusion of water vapor on the CSLAM finite-volume grid, applied as
  ! two successive flux-form Laplacians.
  ! The Laplacian at cell i is the two-point flux approximation
  !   Sum_faces (l_f/d_f)*(Q_j - Q_i) / A_i,
  ! with spherical face arc lengths l_f, center-to-center great-circle
  ! distances d_f, and cell areas A_i; the second pass is layer-mass (dp)
  ! weighted, so tracer mass is conserved exactly and fluxes are antisymmetric
  ! across element and panel boundaries.  The coefficient is
  ! nu4 = se_cslam_q_filter_nu_fac * nu_p
  ! subcycled to satisfy the checkerboard stability
  ! bound dt*nu4*(8/A_min)^2 < 2 (cf. Andrews et al. 2025, arXiv:2505.05624).
  ! Two optional refinements:
  !  (i) a non-orthogonality (cross-diffusion) correction adding the
  !      tangential-gradient contribution the two-point flux omits on the
  !      skewed gnomonic mesh -- the discrete analogue of the g12 metric term
  !      in analytic cubed-sphere Laplacians (Ullrich et al. 2010, JCP);
  ! (ii) Zalesak (1979, JCP) flux-corrected limiting of the del4 fluxes, with
  !      a common coefficient on both sides of each face, so the filter is
  !      shape-preserving (each cell stays within its pre-filter 5-point
  !      min/max) while remaining conservative.
  !
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use dimensions_mod,         only: nc, nlev, cslam_q_filter
  use element_mod,            only: element_t
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod,             only: hybrid_t
  use edgetype_mod,           only: EdgeBuffer_t
  implicit none
  private
  save

  type (EdgeBuffer_t), public :: ghostBufQfilter

  public :: apply_cslam_q_filter_del4   ! the filter (called from run_consistent_se_cslam)
  public :: cslam_q_filter_geom_init    ! precompute face weights + alloc buffer (called from dyn_grid)

contains

  subroutine apply_cslam_q_filter_del4(fvm, hybrid, nets, nete, kmin, kmax, dt_fvm, limiter, xdiff)
    use edge_mod       , only: ghostpack, ghostunpack
    use bndry_mod      , only: ghost_exchange
    use physconst      , only: rearth
    use control_mod    , only: cslam_q_filter_nsub, nu_p
    use dimensions_mod , only: cslam_q_filter_nu_fac
    implicit none
    type (fvm_struct)   , intent(inout) :: fvm(:)
    type (hybrid_t)     , intent(in)    :: hybrid
    integer             , intent(in)    :: nets, nete, kmin, kmax
    real (kind=r8)      , intent(in)    :: dt_fvm
    logical             , intent(in)    :: limiter  ! Zalesak FCT limiter on face fluxes
    logical             , intent(in)    :: xdiff    ! cross-diffusion (non-orthogonality) correction

    integer       , parameter :: ixwv_fvm = 1
    real (kind=r8)            :: rearth4_inv        ! = 1 / rearth^4
    real (kind=r8)            :: nu4k               ! nu_p / rearth^4 (constant in the vertical)

    real (kind=r8), allocatable :: L(:,:,:,:)      ! (0:nc+1, 0:nc+1, kmin:kmax, nets:nete); ring-1 halo only
    !
    ! Limiter work arrays (allocated only when limiter):
    !   Fx/Fy : raw del4 fluxes at e-w/n-s faces (positive toward increasing index)
    !   Rp/Rm : Zalesak ratios for incoming/outgoing flux vs the pre-filter
    !           5-point [min,max]; ghost-exchanged so both sides of a face
    !           use the same pair (keeps fluxes antisymmetric).
    !
    real (kind=r8), allocatable :: Fx(:,:,:,:)     ! (0:nc, 1:nc, kmin:kmax, nets:nete)
    real (kind=r8), allocatable :: Fy(:,:,:,:)     ! (1:nc, 0:nc, kmin:kmax, nets:nete)
    real (kind=r8), allocatable :: Rp(:,:,:,:)     ! (0:nc+1, 0:nc+1, kmin:kmax, nets:nete)
    real (kind=r8), allocatable :: Rm(:,:,:,:)     ! (0:nc+1, 0:nc+1, kmin:kmax, nets:nete)
    real (kind=r8) :: dqw, dqe, dqs, dqn, pp, pm, qmx, qmn, q0
    real (kind=r8) :: c_w, c_e, c_s, c_n
    !
    ! Face-metric weights, copied per element from fvm%qfilter_*
    ! (precomputed at init; see compute_face_weights):
    !   w_ew(i,j) : face between cells (i,j) and (i+1,j), i=0..nc
    !   w_ns(i,j) : face between cells (i,j) and (i,j+1), j=0..nc
    !
    real (kind=r8) :: w_ew(0:nc  , 1:nc  )
    real (kind=r8) :: w_ns(1:nc  , 0:nc  )
    !
    ! xdiff coefficients: e-w face flux = xw_ew*(Q(i+1,j)-Q(i,j)) + xt_ew*T,
    ! T = (Q(i,j+1)+Q(i+1,j+1)) - (Q(i,j-1)+Q(i+1,j-1)); n-s analogous.
    ! Gx/Gy hold per-face fluxes (positive toward increasing index).
    !
    real (kind=r8) :: xw_ew(0:nc  , 1:nc  ), xt_ew(0:nc  , 1:nc  )
    real (kind=r8) :: xw_ns(1:nc  , 0:nc  ), xt_ns(1:nc  , 0:nc  )
    real (kind=r8) :: Gx(0:nc, 1:nc), Gy(1:nc, 0:nc)
    real (kind=r8) :: dp_fw, dp_fe, dp_fs, dp_fn
    real (kind=r8) :: flux_w, flux_e, flux_s, flux_n
    real (kind=r8) :: q_new(nc,nc), inv_dp_area, inv_area
    real (kind=r8) :: dt_sub
    integer        :: ie, k, q, i, j, kblk, kptr, isub

    q           = ixwv_fvm
    kblk        = kmax - kmin + 1
    rearth4_inv = 1.0_r8 / (rearth*rearth*rearth*rearth)
    dt_sub      = dt_fvm / real(cslam_q_filter_nsub, r8)

    allocate(L(0:nc+1, 0:nc+1, kmin:kmax, nets:nete))
    L = 0.0_r8
    if (limiter) then
       allocate(Fx(0:nc  , 1:nc  , kmin:kmax, nets:nete))
       allocate(Fy(1:nc  , 0:nc  , kmin:kmax, nets:nete))
       allocate(Rp(0:nc+1, 0:nc+1, kmin:kmax, nets:nete))
       allocate(Rm(0:nc+1, 0:nc+1, kmin:kmax, nets:nete))
       Rp = 1.0_r8
       Rm = 1.0_r8
    end if

    !
    ! Subcycle loop (cslam_q_filter_nsub auto-set in print_cfl); both halo
    ! exchanges must be inside the loop.
    !
    do isub = 1, cslam_q_filter_nsub

    !----- 1. Halo exchange of dp_fvm and Q_wv -----
    ! ghostBufQfilter is 1-deep: both stencil passes only read halo ring 1.
    do ie = nets, nete
       kptr = kmin - 1
       call ghostpack(ghostBufQfilter, fvm(ie)%dp_fvm(0:nc+1, 0:nc+1, kmin:kmax), kblk, kptr, ie)
       kptr = kmin - 1 + nlev
       call ghostpack(ghostBufQfilter, fvm(ie)%c(0:nc+1, 0:nc+1, kmin:kmax, q), kblk, kptr, ie)
    end do
    call ghost_exchange(hybrid, ghostBufQfilter, location='cslam_q_filter_del4_passQ')
    do ie = nets, nete
       kptr = kmin - 1
       call ghostunpack(ghostBufQfilter, fvm(ie)%dp_fvm(0:nc+1, 0:nc+1, kmin:kmax), kblk, kptr, ie)
       kptr = kmin - 1 + nlev
       call ghostunpack(ghostBufQfilter, fvm(ie)%c(0:nc+1, 0:nc+1, kmin:kmax, q), kblk, kptr, ie)
    end do

    !----- 2. Pass 1: discrete Laplacian L on own cells (1:nc, 1:nc), with metric weights -----
    do ie = nets, nete
       if (xdiff) then
          xw_ew = fvm(ie)%qfilter_xw_ew ; xt_ew = fvm(ie)%qfilter_xt_ew
          xw_ns = fvm(ie)%qfilter_xw_ns ; xt_ns = fvm(ie)%qfilter_xt_ns
       else
          w_ew  = fvm(ie)%qfilter_w_ew  ; w_ns  = fvm(ie)%qfilter_w_ns
       end if
       do k = kmin, kmax
          if (xdiff) then
             ! corrected face fluxes, then divergence per cell
             do j = 1, nc
                do i = 0, nc
                   Gx(i, j) = xw_ew(i, j) * (fvm(ie)%c(i+1, j, k, q) - fvm(ie)%c(i, j, k, q)) + &
                              xt_ew(i, j) * ( (fvm(ie)%c(i  , j+1, k, q) + fvm(ie)%c(i+1, j+1, k, q))   &
                                            - (fvm(ie)%c(i  , j-1, k, q) + fvm(ie)%c(i+1, j-1, k, q)) )
                end do
             end do
             do j = 0, nc
                do i = 1, nc
                   Gy(i, j) = xw_ns(i, j) * (fvm(ie)%c(i, j+1, k, q) - fvm(ie)%c(i, j, k, q)) + &
                              xt_ns(i, j) * ( (fvm(ie)%c(i+1, j  , k, q) + fvm(ie)%c(i+1, j+1, k, q))   &
                                            - (fvm(ie)%c(i-1, j  , k, q) + fvm(ie)%c(i-1, j+1, k, q)) )
                end do
             end do
             do j = 1, nc
                do i = 1, nc
                   inv_area = 1.0_r8 / fvm(ie)%area_sphere(i, j)
                   L(i, j, k, ie) = (Gx(i, j) - Gx(i-1, j) + Gy(i, j) - Gy(i, j-1)) * inv_area
                end do
             end do
          else
             do j = 1, nc
                do i = 1, nc
                   inv_area = 1.0_r8 / fvm(ie)%area_sphere(i, j)
                   L(i, j, k, ie) = (                                                       &
                        w_ew(i-1, j  ) * (fvm(ie)%c(i-1, j  , k, q) - fvm(ie)%c(i, j, k, q)) + &
                        w_ew(i  , j  ) * (fvm(ie)%c(i+1, j  , k, q) - fvm(ie)%c(i, j, k, q)) + &
                        w_ns(i  , j-1) * (fvm(ie)%c(i  , j-1, k, q) - fvm(ie)%c(i, j, k, q)) + &
                        w_ns(i  , j  ) * (fvm(ie)%c(i  , j+1, k, q) - fvm(ie)%c(i, j, k, q))   &
                        ) * inv_area
                end do
             end do
          end if
       end do
    end do

    !----- 3. Halo exchange of L -----
    do ie = nets, nete
       kptr = kmin - 1
       call ghostpack(ghostBufQfilter, L(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
    end do
    call ghost_exchange(hybrid, ghostBufQfilter, location='cslam_q_filter_del4_passL')
    do ie = nets, nete
       kptr = kmin - 1
       call ghostunpack(ghostBufQfilter, L(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
    end do

    !----- 4. Pass 2: flux-form del2 of L, dp-weighted with metric, MINUS sign for del4 -----
    if (.not. limiter) then
    !
    ! Unlimited del4
    !
    do ie = nets, nete
       if (xdiff) then
          xw_ew = fvm(ie)%qfilter_xw_ew ; xt_ew = fvm(ie)%qfilter_xt_ew
          xw_ns = fvm(ie)%qfilter_xw_ns ; xt_ns = fvm(ie)%qfilter_xt_ns
       else
          w_ew  = fvm(ie)%qfilter_w_ew  ; w_ns  = fvm(ie)%qfilter_w_ns
       end if

       do k = kmin, kmax
          nu4k = cslam_q_filter_nu_fac * nu_p * rearth4_inv   ! constant in the vertical
          if (xdiff) then
             ! corrected dp-weighted face fluxes of L, then divergence per cell
             do j = 1, nc
                do i = 0, nc
                   Gx(i, j) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i+1, j, k)) * ( &
                        xw_ew(i, j) * (L(i+1, j, k, ie) - L(i, j, k, ie)) + &
                        xt_ew(i, j) * ( (L(i  , j+1, k, ie) + L(i+1, j+1, k, ie))   &
                                      - (L(i  , j-1, k, ie) + L(i+1, j-1, k, ie)) ) )
                end do
             end do
             do j = 0, nc
                do i = 1, nc
                   Gy(i, j) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i, j+1, k)) * ( &
                        xw_ns(i, j) * (L(i, j+1, k, ie) - L(i, j, k, ie)) + &
                        xt_ns(i, j) * ( (L(i+1, j  , k, ie) + L(i+1, j+1, k, ie))   &
                                      - (L(i-1, j  , k, ie) + L(i-1, j+1, k, ie)) ) )
                end do
             end do
             do j = 1, nc
                do i = 1, nc
                   inv_dp_area = 1.0_r8 / (fvm(ie)%dp_fvm(i, j, k) * fvm(ie)%area_sphere(i, j))
                   q_new(i, j) = fvm(ie)%c(i, j, k, q) - &
                        dt_sub * (Gx(i, j) - Gx(i-1, j) + Gy(i, j) - Gy(i, j-1)) * inv_dp_area
                end do
             end do
          else
          do j = 1, nc
             do i = 1, nc
                dp_fw = 0.5_r8 * (fvm(ie)%dp_fvm(i  , j  , k) + fvm(ie)%dp_fvm(i-1, j  , k))
                dp_fe = 0.5_r8 * (fvm(ie)%dp_fvm(i  , j  , k) + fvm(ie)%dp_fvm(i+1, j  , k))
                dp_fs = 0.5_r8 * (fvm(ie)%dp_fvm(i  , j  , k) + fvm(ie)%dp_fvm(i  , j-1, k))
                dp_fn = 0.5_r8 * (fvm(ie)%dp_fvm(i  , j  , k) + fvm(ie)%dp_fvm(i  , j+1, k))

                flux_w = nu4k * dp_fw * w_ew(i-1, j  ) * (L(i-1, j  , k, ie) - L(i, j, k, ie))
                flux_e = nu4k * dp_fe * w_ew(i  , j  ) * (L(i+1, j  , k, ie) - L(i, j, k, ie))
                flux_s = nu4k * dp_fs * w_ns(i  , j-1) * (L(i  , j-1, k, ie) - L(i, j, k, ie))
                flux_n = nu4k * dp_fn * w_ns(i  , j  ) * (L(i  , j+1, k, ie) - L(i, j, k, ie))

                inv_dp_area = 1.0_r8 / (fvm(ie)%dp_fvm(i, j, k) * fvm(ie)%area_sphere(i, j))
                q_new(i, j) = fvm(ie)%c(i, j, k, q) - &
                     dt_sub * (flux_w + flux_e + flux_s + flux_n) * inv_dp_area
             end do
          end do
          end if ! xdiff
          fvm(ie)%c(1:nc, 1:nc, k, q) = q_new(1:nc, 1:nc)
       end do
    end do

    else
    !
    ! Zalesak FCT del4: each face flux scaled by C = min(Rm_giver,
    ! Rp_receiver) in [0,1] so no cell leaves its pre-filter 5-point
    ! [min,max]; same C on both sides of a face -> conservative.
    !
    !--- 4a. Raw face fluxes and per-cell Zalesak ratios ---
    do ie = nets, nete
       if (xdiff) then
          xw_ew = fvm(ie)%qfilter_xw_ew ; xt_ew = fvm(ie)%qfilter_xt_ew
          xw_ns = fvm(ie)%qfilter_xw_ns ; xt_ns = fvm(ie)%qfilter_xt_ns
       else
          w_ew  = fvm(ie)%qfilter_w_ew  ; w_ns  = fvm(ie)%qfilter_w_ns
       end if

       do k = kmin, kmax
          nu4k = cslam_q_filter_nu_fac * nu_p * rearth4_inv   ! constant in the vertical
          if (xdiff) then
             ! corrected raw fluxes; the FCT limiter below applies unchanged
             ! (it only ever reduces fluxes, so bounds/conservation hold)
             do j = 1, nc
                do i = 0, nc
                   Fx(i, j, k, ie) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i+1, j, k)) * ( &
                        xw_ew(i, j) * (L(i+1, j, k, ie) - L(i, j, k, ie)) + &
                        xt_ew(i, j) * ( (L(i  , j+1, k, ie) + L(i+1, j+1, k, ie))   &
                                      - (L(i  , j-1, k, ie) + L(i+1, j-1, k, ie)) ) )
                end do
             end do
             do j = 0, nc
                do i = 1, nc
                   Fy(i, j, k, ie) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i, j+1, k)) * ( &
                        xw_ns(i, j) * (L(i, j+1, k, ie) - L(i, j, k, ie)) + &
                        xt_ns(i, j) * ( (L(i+1, j  , k, ie) + L(i+1, j+1, k, ie))   &
                                      - (L(i-1, j  , k, ie) + L(i-1, j+1, k, ie)) ) )
                end do
             end do
          else
          do j = 1, nc
             do i = 0, nc
                Fx(i, j, k, ie) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i+1, j, k)) * &
                     w_ew(i, j) * (L(i+1, j, k, ie) - L(i, j, k, ie))
             end do
          end do
          do j = 0, nc
             do i = 1, nc
                Fy(i, j, k, ie) = nu4k * 0.5_r8 * (fvm(ie)%dp_fvm(i, j, k) + fvm(ie)%dp_fvm(i, j+1, k)) * &
                     w_ns(i, j) * (L(i, j+1, k, ie) - L(i, j, k, ie))
             end do
          end do
          end if ! xdiff
          do j = 1, nc
             do i = 1, nc
                inv_dp_area = 1.0_r8 / (fvm(ie)%dp_fvm(i, j, k) * fvm(ie)%area_sphere(i, j))
                ! signed would-be update contribution of each face to cell (i,j)
                dqw =  dt_sub * Fx(i-1, j  , k, ie) * inv_dp_area
                dqe = -dt_sub * Fx(i  , j  , k, ie) * inv_dp_area
                dqs =  dt_sub * Fy(i  , j-1, k, ie) * inv_dp_area
                dqn = -dt_sub * Fy(i  , j  , k, ie) * inv_dp_area
                pp  = max(0.0_r8, dqw) + max(0.0_r8, dqe) + max(0.0_r8, dqs) + max(0.0_r8, dqn)
                pm  = max(0.0_r8,-dqw) + max(0.0_r8,-dqe) + max(0.0_r8,-dqs) + max(0.0_r8,-dqn)
                q0  = fvm(ie)%c(i, j, k, q)
                qmx = max(q0, fvm(ie)%c(i-1, j, k, q), fvm(ie)%c(i+1, j, k, q), &
                              fvm(ie)%c(i, j-1, k, q), fvm(ie)%c(i, j+1, k, q))
                qmn = min(q0, fvm(ie)%c(i-1, j, k, q), fvm(ie)%c(i+1, j, k, q), &
                              fvm(ie)%c(i, j-1, k, q), fvm(ie)%c(i, j+1, k, q))
                Rp(i, j, k, ie) = 1.0_r8
                Rm(i, j, k, ie) = 1.0_r8
                if (pp > 0.0_r8) Rp(i, j, k, ie) = min(1.0_r8, (qmx - q0) / pp)
                if (pm > 0.0_r8) Rm(i, j, k, ie) = min(1.0_r8, (q0 - qmn) / pm)
             end do
          end do
       end do
    end do

    !--- 4b. Exchange Zalesak ratios (ring 1) ---
    do ie = nets, nete
       kptr = kmin - 1
       call ghostpack(ghostBufQfilter, Rp(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
       kptr = kmin - 1 + nlev
       call ghostpack(ghostBufQfilter, Rm(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
    end do
    call ghost_exchange(hybrid, ghostBufQfilter, location='cslam_q_filter_del4_passR')
    do ie = nets, nete
       kptr = kmin - 1
       call ghostunpack(ghostBufQfilter, Rp(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
       kptr = kmin - 1 + nlev
       call ghostunpack(ghostBufQfilter, Rm(0:nc+1, 0:nc+1, kmin:kmax, ie), kblk, kptr, ie)
    end do

    !--- 4c. Limited update ---
    do ie = nets, nete
       do k = kmin, kmax
          do j = 1, nc
             do i = 1, nc
                ! per-face coefficient: min(Rm of giving cell, Rp of receiving cell)
                if (Fx(i-1, j, k, ie) >= 0.0_r8) then      ! west face: (i-1,j) -> (i,j)
                   c_w = min(Rm(i-1, j, k, ie), Rp(i, j, k, ie))
                else                                       ! (i,j) -> (i-1,j)
                   c_w = min(Rp(i-1, j, k, ie), Rm(i, j, k, ie))
                end if
                if (Fx(i  , j, k, ie) >= 0.0_r8) then      ! east face: (i,j) -> (i+1,j)
                   c_e = min(Rm(i, j, k, ie), Rp(i+1, j, k, ie))
                else
                   c_e = min(Rp(i, j, k, ie), Rm(i+1, j, k, ie))
                end if
                if (Fy(i, j-1, k, ie) >= 0.0_r8) then      ! south face: (i,j-1) -> (i,j)
                   c_s = min(Rm(i, j-1, k, ie), Rp(i, j, k, ie))
                else
                   c_s = min(Rp(i, j-1, k, ie), Rm(i, j, k, ie))
                end if
                if (Fy(i, j  , k, ie) >= 0.0_r8) then      ! north face: (i,j) -> (i,j+1)
                   c_n = min(Rm(i, j, k, ie), Rp(i, j+1, k, ie))
                else
                   c_n = min(Rp(i, j, k, ie), Rm(i, j+1, k, ie))
                end if

                inv_dp_area = 1.0_r8 / (fvm(ie)%dp_fvm(i, j, k) * fvm(ie)%area_sphere(i, j))
                q_new(i, j) = fvm(ie)%c(i, j, k, q) - dt_sub * ( &
                     c_w * (-Fx(i-1, j  , k, ie)) + c_e * Fx(i, j, k, ie) + &
                     c_s * (-Fy(i  , j-1, k, ie)) + c_n * Fy(i, j, k, ie) ) * inv_dp_area
             end do
          end do
          fvm(ie)%c(1:nc, 1:nc, k, q) = q_new(1:nc, 1:nc)
       end do
    end do
    end if ! limiter

    end do ! isub subcycle loop

    deallocate(L)
    if (limiter) deallocate(Fx, Fy, Rp, Rm)
  end subroutine apply_cslam_q_filter_del4

  !
  ! One-time precompute of the static face-metric weights (w and xw/xt)
  ! into fvm%qfilter_*.  Call after fvm_init3 (vtx_cart/flux_orient/ifct
  ! halos final).  No-op when cslam_q_filter = .false.
  !
  subroutine cslam_q_filter_geom_init(elem, fvm, hybrid, nets, nete)
    use edge_mod  , only: initghostbuffer
    use thread_mod, only: horz_num_threads, vert_num_threads
    implicit none
    type (element_t) , intent(in)    :: elem(:)
    type (fvm_struct), intent(inout) :: fvm(:)
    type (hybrid_t)  , intent(in)    :: hybrid
    integer          , intent(in)    :: nets, nete

    real (kind=r8) :: c3d_own(3, 1:nc, 1:nc)   ! byproduct, not stored
    integer        :: ie

    !
    ! 1-deep halo buffer for the CSLAM Q filter (its stencils only read
    ! halo ring 1); 2*nlev layers hold dp_fvm and Q in one exchange.
    !
    call initghostbuffer(hybrid%par,ghostBufQfilter,elem,2*nlev,1,nc,nthreads=vert_num_threads*horz_num_threads)
    do ie = nets, nete
       call compute_face_weights(fvm(ie), elem(ie)%FaceNum, c3d_own,   &
            fvm(ie)%qfilter_w_ew , fvm(ie)%qfilter_w_ns ,              &
            fvm(ie)%qfilter_xw_ew, fvm(ie)%qfilter_xt_ew,              &
            fvm(ie)%qfilter_xw_ns, fvm(ie)%qfilter_xt_ns)
    end do
  end subroutine cslam_q_filter_geom_init

  !
  ! Per-face metric weights for one element:
  !   w = face arc-length / center-to-center arc distance.
  ! Cell centers = mean-of-4-vertices in tan-plane, mapped to 3D (exact
  ! great-circle midpoints under gnomonic projection).  Halo-cell centers
  ! are evaluated in the OWNER's tan-plane frame (panel from flux_orient),
  ! so both elements sharing a face agree to ULP -> antisymmetric fluxes
  ! across element and panel boundaries.  c3d_own returned as byproduct.
  !
  subroutine compute_face_weights(f, face_num, c3d_own, w_ew, w_ns, &
                                  xw_ew, xt_ew, xw_ns, xt_ns)
    implicit none
    type (fvm_struct), intent(in)  :: f
    integer          , intent(in)  :: face_num
    real (kind=r8)   , intent(out) :: c3d_own(3, 1:nc, 1:nc)
    real (kind=r8)   , intent(out) :: w_ew(0:nc, 1:nc)
    real (kind=r8)   , intent(out) :: w_ns(1:nc, 0:nc)
    !
    ! Optional xdiff coefficients.  The two-point (TPFA) flux w*(Q_hi - Q_lo)
    ! misses the part of grad(Q) not along the center-to-center direction d.
    ! xdiff adds it back: with a tangential difference T along t (diagonal
    ! ring-1 pairs), the dual basis of {d,t} about the face anchor m gives
    ! the full tangent-plane gradient, hence two-point plus transverse terms
    !   F = face_len*(g.n) = xw*(Q_hi - Q_lo) + xt*T
    !   xw = face_len*(e1.n)/dist,  xt = 0.5*face_len*(e2.n)/arc_t
    !   e1 = (m x t)/den,  e2 = -(m x d)/den,  den = m.(t x d)
    ! Exact for tangent-plane-linear Q; reduces to TPFA (xw = w, xt = 0) for
    ! d normal to the face, and falls back to it where the tangential stencil
    ! hits a missing wedge cell (f%ifct == 0).  Both sides of a face use the
    ! same owner-frame points, and take the fallback together, so F' = -F to
    ! ULP (conservation).
    !
    real (kind=r8)   , intent(out), optional :: xw_ew(0:nc, 1:nc)
    real (kind=r8)   , intent(out), optional :: xt_ew(0:nc, 1:nc)
    real (kind=r8)   , intent(out), optional :: xw_ns(1:nc, 0:nc)
    real (kind=r8)   , intent(out), optional :: xt_ns(1:nc, 0:nc)

    real (kind=r8) :: c3d(3, 0:nc+1, 0:nc+1)
    real (kind=r8) :: v1(3), v2(3)
    real (kind=r8) :: cx, cy, face_len, dist
    integer        :: i, j, ipanel
    logical        :: do_x
    ! xdiff work variables
    real (kind=r8) :: vmid(3), vd(3), vt(3), vp(3), vm(3), vn(3)
    real (kind=r8) :: e1(3), e2(3), den, arc_t, rnorm
    ! den is a sine of the angle between the d and t directions; below this
    ! the dual basis is ill-conditioned (never approached on a sane mesh)
    real (kind=r8), parameter :: den_min = 0.2_r8

    do_x = present(xw_ew)

    !-- Cell centers in 3D: own cells + face-adjacent ring halo cells --
    ! (+ diagonal ring corners when the xdiff tangential stencil needs them;
    !  skip non-existent wedge cells at cube-vertex elements, f%ifct == 0,
    !  whose vtx_cart/flux_orient halo entries are placeholders)
    c3d = 0.0_r8
    do j = 0, nc+1
       do i = 0, nc+1
          if ((i==0 .or. i==nc+1) .and. (j==0 .or. j==nc+1)) then
             if (.not. do_x) cycle
             if (f%ifct(i,j) == 0) cycle
          end if
          cx = 0.25_r8 * ( f%vtx_cart(1,1,i,j) + f%vtx_cart(2,1,i,j) &
                         + f%vtx_cart(3,1,i,j) + f%vtx_cart(4,1,i,j) )
          cy = 0.25_r8 * ( f%vtx_cart(1,2,i,j) + f%vtx_cart(2,2,i,j) &
                         + f%vtx_cart(3,2,i,j) + f%vtx_cart(4,2,i,j) )
          ipanel = NINT(f%flux_orient(1,i,j))
          call tan_to_3d(cx, cy, ipanel, c3d(:, i, j))
       end do
    end do
    c3d_own(:,:,:) = c3d(:, 1:nc, 1:nc)

    !-- East-west face weights: index i = face between cell i and cell i+1 --
    do j = 1, nc
       do i = 0, nc
          if (i == 0) then
             ! West-boundary face of own cell (1, j): vertices 1 & 4 of (1, j)
             call tan_to_3d(f%vtx_cart(1,1,1,j), f%vtx_cart(1,2,1,j), face_num, v1)
             call tan_to_3d(f%vtx_cart(4,1,1,j), f%vtx_cart(4,2,1,j), face_num, v2)
          else
             ! East face of own cell (i, j): vertices 2 & 3 of (i, j)
             call tan_to_3d(f%vtx_cart(2,1,i,j), f%vtx_cart(2,2,i,j), face_num, v1)
             call tan_to_3d(f%vtx_cart(3,1,i,j), f%vtx_cart(3,2,i,j), face_num, v2)
          end if
          face_len   = arc_length(v1, v2)
          dist       = arc_length(c3d(:, i, j), c3d(:, i+1, j))
          w_ew(i, j) = face_len / dist

          
          if (do_x) then
             xw_ew(i, j) = w_ew(i, j)
             xt_ew(i, j) = 0.0_r8
             if ( f%ifct(i,j+1) /= 0 .and. f%ifct(i+1,j+1) /= 0 .and. &
                  f%ifct(i,j-1) /= 0 .and. f%ifct(i+1,j-1) /= 0 ) then
                ! tangent-plane anchor: normalized mean of the two face cells
                vmid  = c3d(:, i, j) + c3d(:, i+1, j)
                vmid  = vmid / sqrt(vmid(1)**2 + vmid(2)**2 + vmid(3)**2)
                ! primary direction d: center-to-center, projected to tangent plane
                vd    = c3d(:, i+1, j) - c3d(:, i, j)
                vd    = vd - (vd(1)*vmid(1) + vd(2)*vmid(2) + vd(3)*vmid(3)) * vmid
                vd    = vd / sqrt(vd(1)**2 + vd(2)**2 + vd(3)**2)
                ! tangential direction t: -j-side pair mean -> +j-side pair mean
                vp    = c3d(:, i, j+1) + c3d(:, i+1, j+1)
                vp    = vp / sqrt(vp(1)**2 + vp(2)**2 + vp(3)**2)
                vm    = c3d(:, i, j-1) + c3d(:, i+1, j-1)
                vm    = vm / sqrt(vm(1)**2 + vm(2)**2 + vm(3)**2)
                arc_t = arc_length(vp, vm)
                vt    = vp - vm
                vt    = vt - (vt(1)*vmid(1) + vt(2)*vmid(2) + vt(3)*vmid(3)) * vmid
                vt    = vt / sqrt(vt(1)**2 + vt(2)**2 + vt(3)**2)
                ! unit face normal, oriented with d ((v2-v1) x m is tangent by construction)
                vn    = cross3(v2 - v1, vmid)
                vn    = vn / sqrt(vn(1)**2 + vn(2)**2 + vn(3)**2)
                if (vn(1)*vd(1) + vn(2)*vd(2) + vn(3)*vd(3) < 0.0_r8) vn = -vn
                den   = vmid(1)*(vt(2)*vd(3) - vt(3)*vd(2)) &
                      + vmid(2)*(vt(3)*vd(1) - vt(1)*vd(3)) &
                      + vmid(3)*(vt(1)*vd(2) - vt(2)*vd(1))
                if (abs(den) > den_min .and. arc_t > 0.0_r8) then
                   e1 = cross3(vmid, vt) / den
                   e2 = -cross3(vmid, vd) / den
                   xw_ew(i, j) = face_len * (e1(1)*vn(1) + e1(2)*vn(2) + e1(3)*vn(3)) / dist
                   xt_ew(i, j) = 0.5_r8 * face_len * (e2(1)*vn(1) + e2(2)*vn(2) + e2(3)*vn(3)) / arc_t
                end if
             end if
          end if
       end do
    end do

    !-- North-south face weights: index j = face between cell j and cell j+1 --
    do j = 0, nc
       do i = 1, nc
          if (j == 0) then
             ! South-boundary face of own cell (i, 1): vertices 1 & 2 of (i, 1)
             call tan_to_3d(f%vtx_cart(1,1,i,1), f%vtx_cart(1,2,i,1), face_num, v1)
             call tan_to_3d(f%vtx_cart(2,1,i,1), f%vtx_cart(2,2,i,1), face_num, v2)
          else
             ! North face of own cell (i, j): vertices 4 & 3 of (i, j)
             call tan_to_3d(f%vtx_cart(4,1,i,j), f%vtx_cart(4,2,i,j), face_num, v1)
             call tan_to_3d(f%vtx_cart(3,1,i,j), f%vtx_cart(3,2,i,j), face_num, v2)
          end if
          face_len   = arc_length(v1, v2)
          dist       = arc_length(c3d(:, i, j), c3d(:, i, j+1))
          w_ns(i, j) = face_len / dist

          if (do_x) then
             xw_ns(i, j) = w_ns(i, j)
             xt_ns(i, j) = 0.0_r8
             if ( f%ifct(i+1,j) /= 0 .and. f%ifct(i+1,j+1) /= 0 .and. &
                  f%ifct(i-1,j) /= 0 .and. f%ifct(i-1,j+1) /= 0 ) then
                vmid  = c3d(:, i, j) + c3d(:, i, j+1)
                vmid  = vmid / sqrt(vmid(1)**2 + vmid(2)**2 + vmid(3)**2)
                vd    = c3d(:, i, j+1) - c3d(:, i, j)
                vd    = vd - (vd(1)*vmid(1) + vd(2)*vmid(2) + vd(3)*vmid(3)) * vmid
                vd    = vd / sqrt(vd(1)**2 + vd(2)**2 + vd(3)**2)
                ! tangential direction t: -i-side pair mean -> +i-side pair mean
                vp    = c3d(:, i+1, j) + c3d(:, i+1, j+1)
                vp    = vp / sqrt(vp(1)**2 + vp(2)**2 + vp(3)**2)
                vm    = c3d(:, i-1, j) + c3d(:, i-1, j+1)
                vm    = vm / sqrt(vm(1)**2 + vm(2)**2 + vm(3)**2)
                arc_t = arc_length(vp, vm)
                vt    = vp - vm
                vt    = vt - (vt(1)*vmid(1) + vt(2)*vmid(2) + vt(3)*vmid(3)) * vmid
                vt    = vt / sqrt(vt(1)**2 + vt(2)**2 + vt(3)**2)
                vn    = cross3(v2 - v1, vmid)
                vn    = vn / sqrt(vn(1)**2 + vn(2)**2 + vn(3)**2)
                if (vn(1)*vd(1) + vn(2)*vd(2) + vn(3)*vd(3) < 0.0_r8) vn = -vn
                den   = vmid(1)*(vt(2)*vd(3) - vt(3)*vd(2)) &
                      + vmid(2)*(vt(3)*vd(1) - vt(1)*vd(3)) &
                      + vmid(3)*(vt(1)*vd(2) - vt(2)*vd(1))
                if (abs(den) > den_min .and. arc_t > 0.0_r8) then
                   e1 = cross3(vmid, vt) / den
                   e2 = -cross3(vmid, vd) / den
                   xw_ns(i, j) = face_len * (e1(1)*vn(1) + e1(2)*vn(2) + e1(3)*vn(3)) / dist
                   xt_ns(i, j) = 0.5_r8 * face_len * (e2(1)*vn(1) + e2(2)*vn(2) + e2(3)*vn(3)) / arc_t
                end if
             end if
          end if
       end do
    end do

  end subroutine compute_face_weights

  function cross3(a, b) result(c)
    implicit none
    real(r8), intent(in) :: a(3), b(3)
    real(r8) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  !
  ! (X, Y) on face `face_no` of the unit-edge-length cube (i.e., X = tan(alpha),
  ! Y = tan(beta)) => 3D unit vector on the sphere.  Derived from
  ! coordinate_systems_mod::unit_face_based_cube_to_unit_sphere.
  !
  subroutine tan_to_3d(X, Y, face_no, p)
    implicit none
    real(r8), intent(in)  :: X, Y
    integer , intent(in)  :: face_no
    real(r8), intent(out) :: p(3)
    real(r8) :: r
    r = 1.0_r8 / sqrt(1.0_r8 + X*X + Y*Y)
    select case (face_no)
    case (1); p(1) =  r    ; p(2) =  X*r  ; p(3) =  Y*r
    case (2); p(1) = -X*r  ; p(2) =  r    ; p(3) =  Y*r
    case (3); p(1) = -r    ; p(2) = -X*r  ; p(3) =  Y*r
    case (4); p(1) =  X*r  ; p(2) = -r    ; p(3) =  Y*r
    case (5); p(1) =  Y*r  ; p(2) =  X*r  ; p(3) = -r
    case (6); p(1) = -Y*r  ; p(2) =  X*r  ; p(3) =  r
    case default; p(:) = 0.0_r8
    end select
  end subroutine tan_to_3d

  !
  ! Great-circle arc-length between two unit vectors on the sphere.
  ! Uses the 2*asin(|a-b|/2) form for numerical stability at small arcs
  ! (acos of near-unity dot product loses precision).
  !
  function arc_length(a, b) result(arc)
    implicit none
    real(r8), intent(in) :: a(3), b(3)
    real(r8) :: arc, d
    d = sqrt( (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2 )
    arc = 2.0_r8 * asin( min(1.0_r8, 0.5_r8 * d) )
  end function arc_length

end module fvm_filter_mod
