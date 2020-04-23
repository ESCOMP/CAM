!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'fv_dynamics' is the top-level routine for the dynamical core.
!>@details It executes on the dt_atmos (or p_split) loop.

module fv_dynamics_mod

! Modules Included:
! <table>
!   <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>boundary_mod</td>
!     <td>nested_grid_BC_apply_intT</td>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>grav, pi=>pi_8, radius, hlv, rdgas, omega, rvgas, cp_vapor</td>
!   </tr>
!   <tr>
!     <td>diag_manager_mod</td>
!     <td>send_data</td>
!   </tr>
!   <tr>
!     <td>dyn_core_mod</td>
!     <td>dyn_core, del2_cubed, init_ijk_mem</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>fv_time, prt_mxm, range_check, prt_minmax</td>
!   </tr>
!   <tr>
!     <td>fv_fill_mod</td>
!     <td>fill2D</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>cubed_to_latlon, c2l_ord2, g_sum</td>
!   </tr>
!   <tr>
!     <td>fv_mapz_mod</td>
!     <td>compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is_master,group_halo_update_type, start_group_halo_update, complete_group_halo_update</td>
!   </tr>
!   <tr>
!     <td>fv_nesting_mod</td>
!     <td>setup_nested_grid_BCs</td>
!   </tr>
!   <tr>
!     <td>fv_nwp_mod</td>
!     <td>do_adiabatic_init</td>
!   </tr>
!   <tr>
!     <td>fv_restart_mod</td>
!     <td>d2a_setup, d2c_setup</td>
!   </tr>
!   <tr>
!     <td>fv_sg_mod</td>
!     <td>neg_adj3</td>
!   </tr>
!   <tr>
!     <td>fv_sg_mod</td>
!     <td>neg_adj2</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>fv_tracer2d_mod</td>
!     <td>tracer_2d, tracer_2d_1L, tracer_2d_nested</td>
!   </tr>
!   <tr>
!     <td>mpp_mod/td>
!     <td>mpp_pe</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod/td>
!     <td>DGRID_NE, CGRID_NE, mpp_update_domains, domain2D</td>
!   </tr>
!   <tr>
!     <td>fv_sg_mod</td>
!     <td>neg_adj3</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_index</td>
!   </tr>
! </table>

   use constants_mod,       only: grav, pi=>pi_8, radius, hlv, rdgas, omega, rvgas, cp_vapor
   use dyn_core_mod,        only: dyn_core, del2_cubed, init_ijk_mem
   use fv_mapz_mod,         only: compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp
   use fv_tracer2d_mod,     only: tracer_2d, tracer_2d_1L, tracer_2d_nested
   use fv_grid_utils_mod,   only: cubed_to_latlon, c2l_ord2, g_sum
   use fv_fill_mod,         only: fill2D
   use fv_mp_mod,           only: is_master
   use fv_mp_mod,           only: group_halo_update_type
   use fv_mp_mod,           only: start_group_halo_update, complete_group_halo_update
   use fv_timing_mod,       only: timing_on, timing_off
   use diag_manager_mod,    only: send_data
   use fv_diagnostics_mod,  only: fv_time, prt_mxm, range_check, prt_minmax
   use mpp_domains_mod,     only: DGRID_NE, CGRID_NE, mpp_update_domains, domain2D
   use mpp_mod,             only: mpp_pe
   use field_manager_mod,   only: MODEL_ATMOS
   use tracer_manager_mod,  only: get_tracer_index
   use fv_sg_mod,           only: neg_adj3, neg_adj2
   use fv_nesting_mod,      only: setup_nested_grid_BCs
   use fv_regional_mod,     only: regional_boundary_update, set_regional_BCs
   use fv_regional_mod,     only: dump_field, H_STAGGER, U_STAGGER, V_STAGGER
   use fv_regional_mod,     only: a_step, p_step, k_step
   use fv_regional_mod,     only: current_time_in_seconds
   use boundary_mod,        only: nested_grid_BC_apply_intT
   use fv_arrays_mod,       only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_nest_type, fv_diag_type, fv_grid_bounds_type
   use fv_nwp_nudge_mod,    only: do_adiabatic_init

implicit none
   logical :: RF_initialized = .false.
   logical :: pt_initialized = .false.
   logical :: bad_range = .false.
   real, allocatable ::  rf(:)
   integer :: kmax=1
   integer :: k_rf = 0

   real :: agrav
#ifdef HIWPP
   real, allocatable:: u00(:,:,:), v00(:,:,:)
#endif
private
public :: fv_dynamics,calc_tot_energy_dynamics1

contains

!-----------------------------------------------------------------------
!     fv_dynamics :: FV dynamical core driver
!-----------------------------------------------------------------------
 
  subroutine fv_dynamics(npx, npy, npz, nq_tot,  ng, bdt, consv_te, fill,             &
                        reproduce_sum, kappa, cp_air, zvir, ptop, ks, ncnst, n_split, &
                        q_split, u, v, w, delz, hydrostatic, pt, delp, q,             &
                        ps, pe, pk, peln, pkz, phis, q_con, omga, ua, va, uc, vc,     &
                        ak, bk, mfx, mfy, cx, cy, ze0, hybrid_z,                      &
                        gridstruct, flagstruct, neststruct, idiag, bd,                &
                        parent_grid, domain, diss_est,                                &
#if ( defined CALC_ENERGY )
                        qsize,qsize_condensate_loading,qsize_condensate_loading_idx,qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
                        se,ke,wv,wl,wi,wr,wsno,wg,tt,mo,mr,gravit, cpair, rearth,omega_cam,fv3_lcp_moist,fv3_lcv_moist,    &
#endif
                        time_total)

#ifdef CCPP
    use mpp_mod,   only: FATAL, mpp_error
    use CCPP_data, only: CCPP_interstitial
#endif

    real, intent(IN) :: bdt  !< Large time-step
    real, intent(IN) :: consv_te
    real, intent(IN) :: kappa, cp_air
    real, intent(IN) :: zvir, ptop
    real, intent(IN), optional :: time_total

    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: nq_tot             !< transported tracers
    integer, intent(IN) :: ng
    integer, intent(IN) :: ks
    integer, intent(IN) :: ncnst
    integer, intent(IN) :: n_split        !< small-step horizontal dynamics
    integer, intent(IN) :: q_split        !< tracer
    logical, intent(IN) :: fill
    logical, intent(IN) :: reproduce_sum
    logical, intent(IN) :: hydrostatic
    logical, intent(IN) :: hybrid_z       !< Using hybrid_z for remapping

    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u !< D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v !< D grid meridional wind (m/s)
    real, intent(inout) :: w(   bd%isd:  ,bd%jsd:  ,1:) !<  W (m/s)
    real, intent(inout) :: pt(  bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< temperature (K)
    real, intent(inout) :: delp(bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz)  !< pressure thickness (pascal)
    real, intent(inout) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) !< specific humidity and constituents
    real, intent(inout) :: delz(bd%isd:,bd%jsd:,1:)   !< delta-height (m); non-hydrostatic only
    real, intent(inout) ::  ze0(bd%is:, bd%js: ,1:) !< height at edges (m); non-hydrostatic
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed,npz) :: diss_est! diffusion estimate for SKEB
#if ( defined CALC_ENERGY )
   real, intent(IN) :: gravit,cpair,rearth,omega_cam
   real, intent(OUT), dimension(bd%is:bd%ie  ,bd%js:bd%je  ,5)  ::  se,ke,wv,wl,wi,wr,wsno,wg,tt,mo,mr
   integer, intent(IN)  :: qsize
   integer, intent(IN)  :: qsize_condensate_loading
   integer, intent(IN)  :: qsize_condensate_loading_idx(qsize_condensate_loading),qsize_tracer_idx_cam2dyn(qsize)
   real, intent(IN)     :: qsize_condensate_loading_cp(qsize_condensate_loading),qsize_condensate_loading_cv(qsize_condensate_loading)
   logical, intent(IN) :: fv3_lcp_moist,fv3_lcv_moist
#endif
! ze0 no longer used

!-----------------------------------------------------------------------
! Auxilliary pressure arrays:    
! The 5 vars below can be re-computed from delp and ptop.
!-----------------------------------------------------------------------
! dyn_aux:
    real, intent(inout) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           !< Surface pressure (pascal)
    real, intent(inout) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  !< edge pressure (pascal)
    real, intent(inout) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          !< pe**kappa
    real, intent(inout) :: peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           !< ln(pe)
    real, intent(inout) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             !< finite-volume mean pk
    real, intent(inout):: q_con(bd%isd:, bd%jsd:, 1:)
    
!-----------------------------------------------------------------------
! Others:
!-----------------------------------------------------------------------
    real, intent(inout) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)       !< Surface geopotential (g*Z_surf)
    real, intent(inout) :: omga(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   !< Vertical pressure velocity (pa/s)
    real, intent(inout) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) !< (uc,vc) mostly used as the C grid winds
    real, intent(inout) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz)

    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz):: ua, va
    real, intent(in),    dimension(npz+1):: ak, bk

! Accumulated Mass flux arrays: the "Flux Capacitor"
    real, intent(inout) ::  mfx(bd%is:bd%ie+1, bd%js:bd%je,   npz)
    real, intent(inout) ::  mfy(bd%is:bd%ie  , bd%js:bd%je+1, npz)
! Accumulated Courant number arrays
    real, intent(inout) ::  cx(bd%is:bd%ie+1, bd%jsd:bd%jed, npz)
    real, intent(inout) ::  cy(bd%isd:bd%ied ,bd%js:bd%je+1, npz)

    type(fv_grid_type),  intent(inout), target :: gridstruct
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_nest_type),  intent(INOUT) :: neststruct
    type(domain2d),      intent(INOUT) :: domain
    type(fv_atmos_type), pointer, intent(IN) :: parent_grid
    type(fv_diag_type),  intent(IN)    :: idiag

! Local Arrays
      real :: temp(bd%isd:bd%ied ,bd%jsd:bd%jed,npz)  !< temperature (K)
      real :: ws(bd%is:bd%ie,bd%js:bd%je)
      real :: te_2d(bd%is:bd%ie,bd%js:bd%je)
      real ::   teq(bd%is:bd%ie,bd%js:bd%je)
      real :: ps2(bd%isd:bd%ied,bd%jsd:bd%jed)
      real :: m_fac(bd%is:bd%ie,bd%js:bd%je)
      real :: pfull(npz)
      real, dimension(bd%is:bd%ie):: cvm
      real, allocatable :: dp1(:,:,:), dtdt_m(:,:,:), cappa(:,:,:)
      real:: akap, rdg, ph1, ph2, mdt, gam, amdt, u0
      real:: recip_k_split,reg_bc_update_time
      integer :: kord_tracer(ncnst)
      integer :: i,j,k, n, iq, n_map, nq, nwat, k_split
      integer :: sphum, liq_wat = -999, ice_wat = -999      ! GFDL physics
      integer :: rainwat = -999, snowwat = -999, graupel = -999, cld_amt = -999
      integer :: theta_d = -999
      logical used, last_step, do_omega
      integer, parameter :: max_packs=12
      type(group_halo_update_type), save :: i_pack(max_packs)
      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed
      real    :: dt2

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed


!     cv_air =  cp_air - rdgas
      agrav   = 1. / grav
        dt2   = 0.5*bdt
      k_split = flagstruct%k_split
      recip_k_split=1./real(k_split)
      nwat    = flagstruct%nwat
      nq      = nq_tot - flagstruct%dnats
      rdg     = -rdgas * agrav
      allocate ( dp1(isd:ied, jsd:jed, 1:npz) )
      call init_ijk_mem(isd,ied, jsd,jed, npz, dp1, 0.)

#if ( defined CALC_ENERGY )
       se=0.
       ke=0.
       wv=0.
       wl=0.
       wi=0.
       wr=0.
       wsno=0.
       wg=0.
       tt=0.
       mo=0.
       mr=0.
#endif
#ifdef MOIST_CAPPA
      allocate ( cappa(isd:ied,jsd:jed,npz) )
      call init_ijk_mem(isd,ied, jsd,jed, npz, cappa, 0.)
#else
      allocate ( cappa(isd:isd,jsd:jsd,1) )
      cappa = 0.
#endif
      !We call this BEFORE converting pt to virtual potential temperature, 
      !since we interpolate on (regular) temperature rather than theta.
      if (gridstruct%nested .or. ANY(neststruct%child_grids)) then
                                           call timing_on('NEST_BCs')
         call setup_nested_grid_BCs(npx, npy, npz, zvir, ncnst, &
              u, v, w, pt, delp, delz, q, uc, vc, pkz, &
              neststruct%nested, flagstruct%inline_q, flagstruct%make_nh, ng, &
              gridstruct, flagstruct, neststruct, &
              neststruct%nest_timestep, neststruct%tracer_nest_timestep, &
              domain, bd, nwat)

#ifndef SW_DYNAMICS
         if (gridstruct%nested) then
          !Correct halo values have now been set up for BCs; we can go ahead and apply them too...
            call nested_grid_BC_apply_intT(pt, &
                 0, 0, npx, npy, npz, bd, 1., 1., &
                 neststruct%pt_BC, bctype=neststruct%nestbctype  )
#ifdef USE_COND
            call nested_grid_BC_apply_intT(q_con, &
                 0, 0, npx, npy, npz, bd, 1., 1., &
                 neststruct%q_con_BC, bctype=neststruct%nestbctype  )            
#ifdef MOIST_CAPPA
            call nested_grid_BC_apply_intT(cappa, &
                 0, 0, npx, npy, npz, bd, 1., 1., &
                 neststruct%cappa_BC, bctype=neststruct%nestbctype  )            
#endif
#endif
         endif
#endif
                                           call timing_off('NEST_BCs')
      endif

    ! For the regional domain set values valid the beginning of the
    ! current large timestep at the boundary points of the pertinent
    ! prognostic arrays.

      if (flagstruct%regional) then
        call timing_on('Regional_BCs')

        reg_bc_update_time=current_time_in_seconds
        call set_regional_BCs          & !<-- Insert values into the boundary region valid for the start of this large timestep.
              (delp,delz,w,pt                                     &
#ifdef USE_COND
               ,q_con                                             &
#endif
#ifdef MOIST_CAPPA
               ,cappa                                             &
#endif
               ,q,u,v,uc,vc, bd, npz, reg_bc_update_time )

        call timing_off('Regional_BCs')
      endif

      if ( flagstruct%no_dycore ) then
         if ( nwat.eq.2 .and. (.not.hydrostatic) ) then
            sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
         endif
         goto 911
      endif

      if ( nwat==0 ) then
             sphum = 1
           cld_amt = -1   ! to cause trouble if (mis)used
      else
             sphum = get_tracer_index (MODEL_ATMOS, 'sphum')
           liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
           ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
           rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
           snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
           graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
           cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
      endif

      theta_d = get_tracer_index (MODEL_ATMOS, 'theta_d')

#ifdef SW_DYNAMICS
      akap  = 1.
      pfull(1) = 0.5*flagstruct%p_ref
#else
      akap  = kappa

!$OMP parallel do default(none) shared(npz,ak,bk,flagstruct,pfull) &
!$OMP                          private(ph1, ph2,k)
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
         pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo

    if ( hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,nwat,q,q_con,sphum,liq_wat, &
!$OMP      rainwat,ice_wat,snowwat,graupel) private(cvm,i,j,k)
      do k=1,npz
         do j=js,je
#ifdef USE_COND
             call moist_cp(is,ie,isd,ied,jsd,jed, npz, j, k, nwat, sphum, liq_wat, rainwat,    &
                           ice_wat, snowwat, graupel, q, q_con(is:ie,j,k), cvm)
#endif
            do i=is,ie
               dp1(i,j,k) = zvir*q(i,j,k,sphum)
            enddo
         enddo
      enddo
    else
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,zvir,q,q_con,sphum,liq_wat, &
!$OMP                                  rainwat,ice_wat,snowwat,graupel,pkz,flagstruct, & 
!$OMP                                  cappa,kappa,rdg,delp,pt,delz,nwat)              &
!$OMP                          private(cvm,i,j,k)
       do k=1,npz
         if ( flagstruct%moist_phys ) then
           do j=js,je
#ifdef MOIST_CAPPA
             call moist_cv(is,ie,isd,ied,jsd,jed, npz, j, k, nwat, sphum, liq_wat, rainwat,    &
                           ice_wat, snowwat, graupel, q, q_con(is:ie,j,k), cvm)
#endif
             do i=is,ie
                dp1(i,j,k) = zvir*q(i,j,k,sphum)
#ifdef MOIST_CAPPA
               cappa(i,j,k) = rdgas/(rdgas + cvm(i)/(1.+dp1(i,j,k)))
               pkz(i,j,k) = exp(cappa(i,j,k)*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                            (1.+dp1(i,j,k))*(1.-q_con(i,j,k))/delz(i,j,k)) )
#else
               pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
                            (1.+dp1(i,j,k))/delz(i,j,k)) )
! Using dry pressure for the definition of the virtual potential temperature
!              pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*pt(i,j,k)*    &
!                                      (1.-q(i,j,k,sphum))/delz(i,j,k)) )
#endif
             enddo
           enddo
         else
           do j=js,je
             do i=is,ie
               dp1(i,j,k) = 0.
               pkz(i,j,k) = exp(kappa*log(rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k)))
             enddo
           enddo
         endif
       enddo
    endif
      if ( flagstruct%fv_debug ) then
#ifdef MOIST_CAPPA
         call prt_mxm('cappa', cappa, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
#endif
         call prt_mxm('PS',        ps, is, ie, js, je, ng,   1, 0.01, gridstruct%area_64, domain)
         call prt_mxm('T_dyn_b',   pt, is, ie, js, je, ng, npz, 1.,   gridstruct%area_64, domain)
         if ( .not. hydrostatic) call prt_mxm('delz',    delz, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
         call prt_mxm('delp_b ', delp, is, ie, js, je, ng, npz, 0.01, gridstruct%area_64, domain)
         call prt_mxm('pk_b',    pk, is, ie, js, je, 0, npz+1, 1.,gridstruct%area_64, domain)
         call prt_mxm('pkz_b',   pkz,is, ie, js, je, 0, npz,   1.,gridstruct%area_64, domain)
      endif

!---------------------
! Compute Total Energy
!---------------------
      if ( consv_te > 0.  .and. (.not.do_adiabatic_init) ) then
           call compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, npz,        &
                                     u, v, w, delz, pt, delp, q, dp1, pe, peln, phis, &
                                     gridstruct%rsin2, gridstruct%cosa_s, &
                                     zvir, cp_air, rdgas, hlv, te_2d, ua, va, teq,        &
                                     flagstruct%moist_phys, nwat, sphum, liq_wat, rainwat,   &
                                     ice_wat, snowwat, graupel, hydrostatic, idiag%id_te)
           if( idiag%id_te>0 ) then
               used = send_data(idiag%id_te, teq, fv_time)
!              te_den=1.E-9*g_sum(teq, is, ie, js, je, ng, area, 0)/(grav*4.*pi*radius**2)
!              if(is_master())  write(*,*) 'Total Energy Density (Giga J/m**2)=',te_den
           endif
      endif

      if( (flagstruct%consv_am.or.idiag%id_amdt>0) .and. (.not.do_adiabatic_init) ) then
          call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                           ptop, ua, va, u, v, delp, teq, ps2, m_fac)
      endif

      if( .not.flagstruct%RF_fast .and. flagstruct%tau > 0. ) then
        if ( gridstruct%grid_type<4 ) then
!         if ( flagstruct%RF_fast ) then
!            call Ray_fast(abs(dt), npx, npy, npz, pfull, flagstruct%tau, u, v, w,  &
!                          dp_ref, ptop, hydrostatic, flagstruct%rf_cutoff, bd)
!         else
             call Rayleigh_Super(abs(bdt), npx, npy, npz, ks, pfull, phis, flagstruct%tau, u, v, w, pt,  &
                  ua, va, delz, gridstruct%agrid, cp_air, rdgas, ptop, hydrostatic,    &
                 (.not. (neststruct%nested .or. flagstruct%regional)), flagstruct%rf_cutoff, gridstruct, domain, bd)
!         endif
        else
             call Rayleigh_Friction(abs(bdt), npx, npy, npz, ks, pfull, flagstruct%tau, u, v, w, pt,  &
                  ua, va, delz, cp_air, rdgas, ptop, hydrostatic, .true., flagstruct%rf_cutoff, gridstruct, domain, bd)
        endif
      endif

#if ( defined CALC_ENERGY )
      call calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
           ptop,delp,u,v,ua,va,q,pt,phis, 'dAT', &
           qsize,qsize_condensate_loading, &
           qsize_condensate_loading_idx, &
           qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
           se(is:ie,js:je,4), ke(is:ie,js:je,4), wv(is:ie,js:je,4),wl(is:ie,js:je,4), &
           wi(is:ie,js:je,4), wr(is:ie,js:je,4),wsno(is:ie,js:je,4),wg(is:ie,js:je,4), &
           tt(is:ie,js:je,4),mo(is:ie,js:je,4),mr(is:ie,js:je,4), &
           gravit, cpair, rearth,omega_cam,1,.true.,fv3_lcp_moist,fv3_lcv_moist)
#endif
      

      
#endif

#ifndef SW_DYNAMICS
! Convert pt to virtual potential temperature on the first timestep
  if ( flagstruct%adiabatic .and. flagstruct%kord_tm>0 ) then
     if ( .not.pt_initialized )then
!$OMP parallel do default(none) shared(theta_d,is,ie,js,je,npz,pt,pkz,q)
       do k=1,npz
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
             enddo
          enddo
          if ( theta_d>0 ) then
             do j=js,je
                do i=is,ie
                   q(i,j,k,theta_d) = pt(i,j,k)
                enddo
             enddo
          endif
       enddo
       pt_initialized = .true.
     endif
  else
!$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,dp1,pkz,q_con)
  do k=1,npz
     do j=js,je
        do i=is,ie
#ifdef USE_COND
           pt(i,j,k) = pt(i,j,k)*(1.+dp1(i,j,k))*(1.-q_con(i,j,k))/pkz(i,j,k)
#else
           pt(i,j,k) = pt(i,j,k)*(1.+dp1(i,j,k))/pkz(i,j,k)
#endif
        enddo
     enddo
  enddo
  endif
#endif

  last_step = .false.
  mdt = bdt / real(k_split)

  allocate ( dtdt_m(is:ie,js:je,npz) )
  if ( idiag%id_mdt > 0 .and. (.not. do_adiabatic_init) ) then
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m)
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = 0.
             enddo
          enddo
       enddo
  endif

                                                  call timing_on('FV_DYN_LOOP')
  do n_map=1, k_split   ! first level of time-split
      k_step = n_map

      if ( n_map==k_split ) last_step = .true.

#ifdef USE_COND
      !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp,q_con)
#else
      !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp)
#endif
      do k=1,npz
         do j=js,je
            do i=is,ie
#ifdef USE_COND
               temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/((1.+zvir*q(i,j,k,sphum))*(1.-q_con(i,j,k)))
#else
               temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/(1.+zvir*q(i,j,k,sphum))
#endif
            enddo
         enddo
      enddo

#if ( defined CALC_ENERGY )
      call calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
           ptop,delp,u,v,ua,va,q,temp,phis, 'dAF', &
           qsize,qsize_condensate_loading, &
           qsize_condensate_loading_idx, &
           qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
           se(is:ie,js:je,1), ke(is:ie,js:je,1), wv(is:ie,js:je,1),wl(is:ie,js:je,1), &
           wi(is:ie,js:je,1), wr(is:ie,js:je,1),wsno(is:ie,js:je,1),wg(is:ie,js:je,1), &
           tt(is:ie,js:je,1),mo(is:ie,js:je,1),mr(is:ie,js:je,1), &
           gravit, cpair, rearth,omega_cam,k_split,last_step,fv3_lcp_moist,fv3_lcv_moist)
#endif

                                           call timing_on('COMM_TOTAL')
#ifdef USE_COND
      call start_group_halo_update(i_pack(11), q_con, domain)
#ifdef MOIST_CAPPA
      call start_group_halo_update(i_pack(12), cappa, domain)
#endif
#endif
      call start_group_halo_update(i_pack(1), delp, domain, complete=.false.)
      call start_group_halo_update(i_pack(1), pt,   domain, complete=.true.)
#ifndef ROT3
      call start_group_halo_update(i_pack(8), u, v, domain, gridtype=DGRID_NE)
#endif
                                           call timing_off('COMM_TOTAL')
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,npz,dp1,delp)
      do k=1,npz
         do j=jsd,jed
            do i=isd,ied
               dp1(i,j,k) = delp(i,j,k)
            enddo
         enddo
      enddo

      if ( n_map==k_split ) last_step = .true.

#ifdef USE_COND
                                           call timing_on('COMM_TOTAL')
     call complete_group_halo_update(i_pack(11), domain)
#ifdef MOIST_CAPPA
     call complete_group_halo_update(i_pack(12), domain)
#endif
                                           call timing_off('COMM_TOTAL')
#endif

                                           call timing_on('DYN_CORE')
      call dyn_core(npx, npy, npz, ng, sphum, nq, mdt, n_split, zvir, cp_air, akap, cappa, grav, hydrostatic, &
                    u, v, w, delz, pt, q, delp, pe, pk, phis, ws, omga, ptop, pfull, ua, va,           & 
                    uc, vc, mfx, mfy, cx, cy, pkz, peln, q_con, ak, bk, ks, &
                    gridstruct, flagstruct, neststruct, idiag, bd, &
                    domain, n_map==1, i_pack, last_step, diss_est,time_total)
                                           call timing_off('DYN_CORE')
#ifdef USE_COND
      !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp,q_con)
#else
      !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp)
#endif
      do k=1,npz
         do j=js,je
            do i=is,ie
#ifdef USE_COND
               temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/((1.+zvir*q(i,j,k,sphum))*(1.-q_con(i,j,k)))
#else
               temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/(1.+zvir*q(i,j,k,sphum))
#endif
            enddo
         enddo
      enddo
      
#if ( defined CALC_ENERGY )
      call calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
           ptop,delp,u,v,ua,va,q,temp,phis, 'dAD', &
           qsize,qsize_condensate_loading, &
           qsize_condensate_loading_idx, &
           qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
           se(is:ie,js:je,2), ke(is:ie,js:je,2), wv(is:ie,js:je,2),wl(is:ie,js:je,2), &
           wi(is:ie,js:je,2), wr(is:ie,js:je,2),wsno(is:ie,js:je,2), wg(is:ie,js:je,2), &
           tt(is:ie,js:je,2),mo(is:ie,js:je,2),mr(is:ie,js:je,2), &
           gravit, cpair, rearth,omega_cam,k_split,last_step,fv3_lcp_moist,fv3_lcv_moist)
#endif      
#ifdef SW_DYNAMICS
!!$OMP parallel do default(none) shared(is,ie,js,je,ps,delp,agrav)
      do j=js,je
         do i=is,ie
            ps(i,j) = delp(i,j,1) * agrav
         enddo
      enddo
#else
      if( .not. flagstruct%inline_q .and. nq /= 0 ) then    
!--------------------------------------------------------
! Perform large-time-step scalar transport using the accumulated CFL and
! mass fluxes
                                              call timing_on('tracer_2d')
       !!! CLEANUP: merge these two calls?
       if (gridstruct%nested .or. flagstruct%regional) then
         call tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%nord_tr, flagstruct%trdm2, &
                        k_split, neststruct, parent_grid, flagstruct%lim_fac,flagstruct%regional)
       else
         if ( flagstruct%z_tracer ) then
         call tracer_2d_1L(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac,flagstruct%regional)
         else
         call tracer_2d(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz, nq,    &
                        flagstruct%hord_tr, q_split, mdt, idiag%id_divg, i_pack(10), &
                        flagstruct%nord_tr, flagstruct%trdm2, flagstruct%lim_fac,flagstruct%regional)
         endif
       endif
                                             call timing_off('tracer_2d')

#ifdef FILL2D
     if ( flagstruct%hord_tr<8 .and. flagstruct%moist_phys ) then
                                                  call timing_on('Fill2D')
       if ( liq_wat > 0 )  &
        call fill2D(is, ie, js, je, ng, npz, q(isd,jsd,1,liq_wat), delp, gridstruct%area, domain, neststruct%nested, gridstruct%regional, npx, npy)
       if ( rainwat > 0 )  &
        call fill2D(is, ie, js, je, ng, npz, q(isd,jsd,1,rainwat), delp, gridstruct%area, domain, neststruct%nested, gridstruct%regional, npx, npy)
       if ( ice_wat > 0  )  &
        call fill2D(is, ie, js, je, ng, npz, q(isd,jsd,1,ice_wat), delp, gridstruct%area, domain, neststruct%nested, gridstruct%regional, npx, npy)
       if ( snowwat > 0 )  &
        call fill2D(is, ie, js, je, ng, npz, q(isd,jsd,1,snowwat), delp, gridstruct%area, domain, neststruct%nested, gridstruct%regional, npx, npy)
       if ( graupel > 0 )  &
        call fill2D(is, ie, js, je, ng, npz, q(isd,jsd,1,graupel), delp, gridstruct%area, domain, neststruct%nested, gridstruct%regional, npx, npy)
                                                  call timing_off('Fill2D')
     endif
#endif

         if( last_step .and. idiag%id_divg>0 ) then
             used = send_data(idiag%id_divg, dp1, fv_time) 
             if(flagstruct%fv_debug) call prt_mxm('divg',  dp1, is, ie, js, je, 0, npz, 1.,gridstruct%area_64, domain)
         endif
      endif

      if ( npz > 4 ) then
!------------------------------------------------------------------------
! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.
!------------------------------------------------------------------------

         do iq=1,nq
                                kord_tracer(iq) = flagstruct%kord_tr
            if ( iq==cld_amt )  kord_tracer(iq) = 9      ! monotonic
         enddo

         do_omega = hydrostatic .and. last_step
                                                  call timing_on('Remapping')
#ifdef AVEC_TIMERS
                                                  call avec_timer_start(6)
#endif

         call Lagrangian_to_Eulerian(last_step, consv_te, ps, pe, delp,                 &
                     pkz, pk, mdt, bdt, npz, is,ie,js,je, isd,ied,jsd,jed,              &
                     nq, nwat, sphum, q_con, u,  v, w, delz, pt, q, phis,               &
                     zvir, cp_air, akap, cappa, flagstruct%kord_mt, flagstruct%kord_wz, &
                     kord_tracer, flagstruct%kord_tm, peln, te_2d,                      &
                     ng, ua, va, omga, dp1, ws, fill, reproduce_sum,                    &
                     idiag%id_mdt>0, dtdt_m, ptop, ak, bk, pfull, gridstruct, domain,   &
                     flagstruct%do_sat_adj, hydrostatic, hybrid_z, do_omega,            &
                     flagstruct%adiabatic, do_adiabatic_init)

#if ( defined CALC_ENERGY )
         if (last_step) then
            call calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
                 ptop,delp,u,v,ua,va,q,pt,phis, 'dAR', &
                 qsize,qsize_condensate_loading, &
                 qsize_condensate_loading_idx, &
                 qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
                 se(is:ie,js:je,3), ke(is:ie,js:je,3), wv(is:ie,js:je,3),wl(is:ie,js:je,3), &
                 wi(is:ie,js:je,3), wr(is:ie,js:je,3),wsno(is:ie,js:je,3),wg(is:ie,js:je,3), &
                 tt(is:ie,js:je,3),mo(is:ie,js:je,3),mr(is:ie,js:je,3), &
                 gravit, cpair, rearth,omega_cam,k_split,last_step,fv3_lcp_moist,fv3_lcv_moist)

         else
#ifdef USE_COND
            !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp,q_con)
#else
            !$OMP parallel do default(none) shared(is,ie,js,je,npz,pt,q,pkz,zvir,sphum,temp)
#endif
            do k=1,npz
               do j=js,je
                  do i=is,ie
#ifdef USE_COND
                     temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/((1.+zvir*q(i,j,k,sphum))*(1.-q_con(i,j,k)))
#else
                     temp(i,j,k)=pt(i,j,k)*pkz(i,j,k)/(1.+zvir*q(i,j,k,sphum))
#endif
                  enddo
               enddo
             enddo
             call calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
                 ptop,delp,u,v,ua,va,q,temp,phis, 'dAR', &
                 qsize,qsize_condensate_loading, &
                 qsize_condensate_loading_idx, &
                 qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
                 se(is:ie,js:je,3), ke(is:ie,js:je,3),wv(is:ie,js:je,3),wl(is:ie,js:je,3), &
                 wi(is:ie,js:je,3),wr(is:ie,js:je,3),wsno(is:ie,js:je,3),wg(is:ie,js:je,3), &
                 tt(is:ie,js:je,3),mo(is:ie,js:je,3),mr(is:ie,js:je,3), &
                 gravit, cpair, rearth,omega_cam,k_split,last_step,fv3_lcp_moist,fv3_lcv_moist)
            
         end if
#endif
#ifdef AVEC_TIMERS
                                                  call avec_timer_stop(6)
#endif
                                                  call timing_off('Remapping')
#ifdef MOIST_CAPPA
         if ( neststruct%nested .and. .not. last_step) then
            call nested_grid_BC_apply_intT(cappa, &
                 0, 0, npx, npy, npz, bd, real(n_map+1), real(k_split), &
                 neststruct%cappa_BC, bctype=neststruct%nestbctype  )
         endif

         if ( flagstruct%regional .and. .not. last_step) then
            reg_bc_update_time=current_time_in_seconds+(n_map+1)*mdt
            call regional_boundary_update(cappa, 'cappa', &
                                          isd, ied, jsd, jed, npz, &
                                          is,  ie,  js,  je,       &
                                          isd, ied, jsd, jed,      &
                                          reg_bc_update_time )
         endif
#endif

         if( last_step )  then
            if( .not. hydrostatic ) then
!$OMP parallel do default(none) shared(is,ie,js,je,npz,omga,delp,delz,w)
               do k=1,npz
                  do j=js,je
                     do i=is,ie
                        omga(i,j,k) = delp(i,j,k)/delz(i,j,k)*w(i,j,k)
                     enddo
                  enddo
               enddo
            endif
!--------------------------
! Filter omega for physics:
!--------------------------
            if(flagstruct%nf_omega>0)    &
            call del2_cubed(omga, 0.18*gridstruct%da_min, gridstruct, domain, npx, npy, npz, flagstruct%nf_omega, bd)
         endif
      end if
#endif
  enddo    ! n_map loop
                                                  call timing_off('FV_DYN_LOOP')
  if ( idiag%id_mdt > 0 .and. (.not.do_adiabatic_init) ) then
! Output temperature tendency due to inline moist physics:
!$OMP parallel do default(none) shared(is,ie,js,je,npz,dtdt_m,bdt)
       do k=1,npz
          do j=js,je
             do i=is,ie
                dtdt_m(i,j,k) = dtdt_m(i,j,k) *( 86400.0 / bdt)
             enddo
          enddo
       enddo
!      call prt_mxm('Fast DTDT (deg/Day)', dtdt_m, is, ie, js, je, 0, npz, 1., gridstruct%area_64, domain)
       used = send_data(idiag%id_mdt, dtdt_m, fv_time)
  endif
  deallocate ( dtdt_m )

  if( nwat == 6 ) then
     if (cld_amt > 0) then
      call neg_adj3(is, ie, js, je, ng, npz,        &
                    flagstruct%hydrostatic,         &
                    peln, delz,                     &
                    pt, delp, q(isd,jsd,1,sphum),   &
                              q(isd,jsd,1,liq_wat), &
                              q(isd,jsd,1,rainwat), &
                              q(isd,jsd,1,ice_wat), &
                              q(isd,jsd,1,snowwat), &
                              q(isd,jsd,1,graupel), &
                              q(isd,jsd,1,cld_amt), flagstruct%check_negative)
     else
        call neg_adj3(is, ie, js, je, ng, npz,        &
                      flagstruct%hydrostatic,         &
                      peln, delz,                     &
                      pt, delp, q(isd,jsd,1,sphum),   &
                                q(isd,jsd,1,liq_wat), &
                                q(isd,jsd,1,rainwat), &
                                q(isd,jsd,1,ice_wat), &
                                q(isd,jsd,1,snowwat), &
                                q(isd,jsd,1,graupel), check_negative=flagstruct%check_negative)
     endif
     if ( flagstruct%fv_debug ) then
       call prt_mxm('T_dyn_a3',    pt, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
       call prt_mxm('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('rainwat_dyn', q(isd,jsd,1,rainwat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('graupel_dyn', q(isd,jsd,1,graupel), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
     endif
  endif

  if( nwat == 5 ) then
     if (cld_amt > 0) then
      call neg_adj2(is, ie, js, je, ng, npz,        &
                    flagstruct%hydrostatic,         &
                    peln, delz,                     &
                    pt, delp, q(isd,jsd,1,sphum),   &
                              q(isd,jsd,1,liq_wat), &
                              q(isd,jsd,1,rainwat), &
                              q(isd,jsd,1,ice_wat), &
                              q(isd,jsd,1,snowwat), &
                              q(isd,jsd,1,cld_amt), flagstruct%check_negative)
     else
        call neg_adj2(is, ie, js, je, ng, npz,        &
                      flagstruct%hydrostatic,         &
                      peln, delz,                     &
                      pt, delp, q(isd,jsd,1,sphum),   &
                                q(isd,jsd,1,liq_wat), &
                                q(isd,jsd,1,rainwat), &
                                q(isd,jsd,1,ice_wat), &
                                q(isd,jsd,1,snowwat), &
                                check_negative=flagstruct%check_negative)
     endif
     if ( flagstruct%fv_debug ) then
       call prt_mxm('T_dyn_a3',    pt, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
       call prt_mxm('SPHUM_dyn',   q(isd,jsd,1,sphum  ), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('liq_wat_dyn', q(isd,jsd,1,liq_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('rainwat_dyn', q(isd,jsd,1,rainwat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('ice_wat_dyn', q(isd,jsd,1,ice_wat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
       call prt_mxm('snowwat_dyn', q(isd,jsd,1,snowwat), is, ie, js, je, ng, npz, 1.,gridstruct%area_64, domain)
     endif
  endif

  if( (flagstruct%consv_am.or.idiag%id_amdt>0.or.idiag%id_aam>0) .and. (.not.do_adiabatic_init)  ) then
      call compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                       ptop, ua, va, u, v, delp, te_2d, ps, m_fac)
      if( idiag%id_aam>0 ) then
          used = send_data(idiag%id_aam, te_2d, fv_time)
          if ( prt_minmax ) then
             gam = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0) 
             if( is_master() ) write(6,*) 'Total AAM =', gam
          endif
      endif
  endif

  if( (flagstruct%consv_am.or.idiag%id_amdt>0) .and. (.not.do_adiabatic_init)  ) then
!$OMP parallel do default(none) shared(is,ie,js,je,te_2d,teq,dt2,ps2,ps,idiag) 
      do j=js,je
         do i=is,ie
! Note: the mountain torque computation contains also numerical error
! The numerical error is mostly from the zonal gradient of the terrain (zxg)
            te_2d(i,j) = te_2d(i,j)-teq(i,j) + dt2*(ps2(i,j)+ps(i,j))*idiag%zxg(i,j)
         enddo
      enddo
      if( idiag%id_amdt>0 ) used = send_data(idiag%id_amdt, te_2d/bdt, fv_time)

      if ( flagstruct%consv_am .or. prt_minmax ) then
         amdt = g_sum( domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.) 
         u0 = -radius*amdt/g_sum( domain, m_fac, is, ie, js, je, ng, gridstruct%area_64, 0,reproduce=.true.)
         if(is_master() .and. prt_minmax)         &
         write(6,*) 'Dynamic AM tendency (Hadleys)=', amdt/(bdt*1.e18), 'del-u (per day)=', u0*86400./bdt
      endif

    if( flagstruct%consv_am ) then
!$OMP parallel do default(none) shared(is,ie,js,je,m_fac,u0,gridstruct)
      do j=js,je
         do i=is,ie
            m_fac(i,j) = u0*cos(gridstruct%agrid(i,j,2))
         enddo
      enddo
!$OMP parallel do default(none) shared(is,ie,js,je,npz,hydrostatic,pt,m_fac,ua,cp_air, &
!$OMP                                  u,u0,gridstruct,v )
      do k=1,npz
      do j=js,je+1
         do i=is,ie
            u(i,j,k) = u(i,j,k) + u0*gridstruct%l2c_u(i,j)
         enddo
      enddo
      do j=js,je
         do i=is,ie+1
            v(i,j,k) = v(i,j,k) + u0*gridstruct%l2c_v(i,j)
         enddo
      enddo
      enddo
    endif   !  consv_am
  endif

911  call cubed_to_latlon(u, v, ua, va, gridstruct, &
          npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)

  deallocate(dp1)
  deallocate(cappa)

  if ( flagstruct%fv_debug ) then
     call prt_mxm('UA', ua, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
     call prt_mxm('VA', va, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
     call prt_mxm('TA', pt, is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
     if (.not. hydrostatic) call prt_mxm('W ', w,  is, ie, js, je, ng, npz, 1., gridstruct%area_64, domain)
  endif

  if ( flagstruct%range_warn ) then
       call range_check('UA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -280., 280., bad_range)
       call range_check('VA_dyn', ua, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -280., 280., bad_range)
       call range_check('TA_dyn', pt, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         150., 335., bad_range)
       if ( .not. hydrostatic ) &
            call range_check('W_dyn', w, is, ie, js, je, ng, npz, gridstruct%agrid,   &
                         -50., 100., bad_range)
  endif

  end subroutine fv_dynamics

#ifdef USE_RF_FAST
  subroutine Rayleigh_fast(dt, npx, npy, npz, pfull, tau, u, v, w,  &
                          ks, dp, ptop, hydrostatic, rf_cutoff, bd)
! Simple "inline" version of the Rayleigh friction
    real, intent(in):: dt
    real, intent(in):: tau              !< time scale (days)
    real, intent(in):: ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pfull
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) !< D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) !< D grid meridional wind (m/s)
    real, intent(inout):: w(bd%isd:      ,bd%jsd:      ,1: ) !< cell center vertical wind (m/s)
    real, intent(in):: dp(npz)
!
    real(kind=R_GRID):: rff(npz)
    real, parameter:: sday = 86400. !< seconds per day
    real, dimension(bd%is:bd%ie+1):: dmv
    real, dimension(bd%is:bd%ie):: dmu
    real:: tau0, dm
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

     if ( .not. RF_initialized ) then
          tau0 = tau * sday
          allocate( rf(npz) )
          rf(:) = 1.

          if( is_master() ) write(6,*) 'Fast Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pfull(k) < rf_cutoff ) then
                  rff(k) = dt/tau0*sin(0.5*pi*log(rf_cutoff/pfull(k))/log(rf_cutoff/ptop))**2
! Re-FACTOR rf
                  if( is_master() ) write(6,*) k, 0.01*pfull(k), dt/(rff(k)*sday)
                  kmax = k
                  rff(k) = 1.d0 / (1.0d0+rff(k))
                   rf(k) = rff(k)
             else
                  exit
             endif
          enddo
          dm = 0.
          do k=1,ks
             if ( pfull(k) < 100.E2 ) then
                  dm = dm + dp(k)
                  k_rf = k
             else
                  exit
             endif
          enddo
          if( is_master() ) write(6,*) 'k_rf=', k_rf, 0.01*pfull(k_rf), 'dm=', dm
          RF_initialized = .true.
     endif

!$OMP parallel do default(none) shared(k_rf,is,ie,js,je,kmax,pfull,rf_cutoff,w,rf,dp,u,v,hydrostatic) &
!$OMP          private(dm, dmu, dmv)
     do j=js,je+1

        dm = 0.
        do k=1, k_rf
           dm = dm + dp(k)
        enddo

        dmu(:) = 0.
        dmv(:) = 0.
        do k=1,kmax
           do i=is,ie
                dmu(i) = dmu(i) + (1.-rf(k))*dp(k)*u(i,j,k)
              u(i,j,k) = rf(k)*u(i,j,k)
           enddo
           if ( j/=je+1 ) then
              do i=is,ie+1
                   dmv(i) = dmv(i) + (1.-rf(k))*dp(k)*v(i,j,k)
                 v(i,j,k) = rf(k)*v(i,j,k)
              enddo
              if ( .not. hydrostatic ) then
                do i=is,ie
                   w(i,j,k) = rf(k)*w(i,j,k)
                enddo
              endif
           endif
        enddo

        do i=is,ie
           dmu(i) = dmu(i) / dm
        enddo
        if ( j/=je+1 ) then
           do i=is,ie+1
              dmv(i) = dmv(i) / dm
           enddo
        endif

        do k=1, k_rf
           do i=is,ie
              u(i,j,k) = u(i,j,k) + dmu(i)
           enddo
           if ( j/=je+1 ) then
              do i=is,ie+1
                 v(i,j,k) = v(i,j,k) + dmv(i)
              enddo
           endif
        enddo

     enddo

 end subroutine Rayleigh_fast
#endif



 subroutine Rayleigh_Super(dt, npx, npy, npz, ks, pm, phis, tau, u, v, w, pt,  &
                           ua, va, delz, agrid, cp, rg, ptop, hydrostatic,     &
                           conserve, rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              !< time scale (days)
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) !< D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) !< D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:      ,bd%jsd:      ,1: ) !< cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) !< temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:    ,bd%jsd:      ,1: ) !< delta-height (m); non-hydrostatic only
    real,   intent(in) :: agrid(bd%isd:bd%ied,  bd%jsd:bd%jed,2)
    real, intent(in) :: phis(bd%isd:bd%ied,bd%jsd:bd%jed)     !< Surface geopotential (g*Z_surf)
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
!
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: u0   = 60.   !< scaling velocity
    real, parameter:: sday = 86400.
    real rcv, tau0
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

    rcv = 1. / (cp - rg)

     if ( .not. RF_initialized ) then
#ifdef HIWPP
          allocate ( u00(is:ie,  js:je+1,npz) )
          allocate ( v00(is:ie+1,js:je  ,npz) )
!$OMP parallel do default(none) shared(is,ie,js,je,npz,u00,u,v00,v)
          do k=1,npz
             do j=js,je+1
                do i=is,ie
                   u00(i,j,k) = u(i,j,k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v00(i,j,k) = v(i,j,k)
                enddo
             enddo
          enddo
#endif
#ifdef SMALL_EARTH
          tau0 = tau
#else
          tau0 = tau * sday
#endif
          allocate( rf(npz) )
          rf(:) = 0.

          do k=1, ks+1
             if( is_master() ) write(6,*) k, 0.01*pm(k)
          enddo
          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/tau0*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
     endif

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd, gridstruct%nested)

    allocate( u2f(isd:ied,jsd:jed,kmax) )

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,hydrostatic,ua,va,agrid, &
!$OMP                                  u2f,rf,w)
    do k=1,kmax
       if ( pm(k) < rf_cutoff ) then
          u2f(:,:,k) = 1. / (1.+rf(k))
       else
          u2f(:,:,k) = 1.
       endif
    enddo
                                        call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,pm,rf_cutoff,w,rf,u,v, &
#ifdef HIWPP
!$OMP                                  u00,v00, &
#endif
!$OMP                                  conserve,hydrostatic,pt,ua,va,u2f,cp,rg,ptop,rcv)
     do k=1,kmax
        if ( pm(k) < rf_cutoff ) then
#ifdef HIWPP
           if (.not. hydrostatic) then
              do j=js,je
                 do i=is,ie
                    w(i,j,k) = w(i,j,k)/(1.+rf(k))
                 enddo
              enddo
           endif
             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = (u(i,j,k)+rf(k)*u00(i,j,k))/(1.+rf(k))
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = (v(i,j,k)+rf(k)*v00(i,j,k))/(1.+rf(k))
                enddo
             enddo
#else
! Add heat so as to conserve TE
          if ( conserve ) then
             if ( hydrostatic ) then
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2)*(1.-u2f(i,j,k)**2)/(cp-rg*ptop/pm(k))
                  enddo
               enddo
             else
               do j=js,je
                  do i=is,ie
                     pt(i,j,k) = pt(i,j,k) + 0.5*(ua(i,j,k)**2+va(i,j,k)**2+w(i,j,k)**2)*(1.-u2f(i,j,k)**2)*rcv
                  enddo
               enddo
             endif
          endif

             do j=js,je+1
                do i=is,ie
                   u(i,j,k) = 0.5*(u2f(i,j-1,k)+u2f(i,j,k))*u(i,j,k)
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   v(i,j,k) = 0.5*(u2f(i-1,j,k)+u2f(i,j,k))*v(i,j,k)
                enddo
             enddo
          if ( .not. hydrostatic ) then
             do j=js,je
                do i=is,ie
                   w(i,j,k) = u2f(i,j,k)*w(i,j,k)
                enddo
             enddo
          endif
#endif
        endif
     enddo

     deallocate ( u2f )

 end subroutine Rayleigh_Super


 subroutine Rayleigh_Friction(dt, npx, npy, npz, ks, pm, tau, u, v, w, pt,  &
                              ua, va, delz, cp, rg, ptop, hydrostatic, conserve, &
                              rf_cutoff, gridstruct, domain, bd)
    real, intent(in):: dt
    real, intent(in):: tau              !< time scale (days)
    real, intent(in):: cp, rg, ptop, rf_cutoff
    real, intent(in),  dimension(npz):: pm
    integer, intent(in):: npx, npy, npz, ks
    logical, intent(in):: hydrostatic
    logical, intent(in):: conserve
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout):: u(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) !< D grid zonal wind (m/s)
    real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz) !< D grid meridional wind (m/s)
    real, intent(inout)::  w(bd%isd:      ,bd%jsd:      ,1: ) !< cell center vertical wind (m/s)
    real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz) !< temp
    real, intent(inout):: ua(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: va(bd%isd:bd%ied,bd%jsd:bd%jed,npz) ! 
    real, intent(inout):: delz(bd%isd:    ,bd%jsd:      ,1: ) !< delta-height (m); non-hydrostatic only
    type(fv_grid_type), intent(IN) :: gridstruct
    type(domain2d), intent(INOUT) :: domain
! local:
    real, allocatable ::  u2f(:,:,:)
    real, parameter:: sday = 86400.  !< seconds per day
    real, parameter:: u000 = 4900.   !< scaling velocity  **2
    real  rcv
    integer i, j, k

    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    
    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed


    rcv = 1. / (cp - rg)

    if ( .not. RF_initialized ) then
          allocate( rf(npz) )
          if( is_master() ) write(6,*) 'Rayleigh friction E-folding time (days):'
          do k=1, npz
             if ( pm(k) < rf_cutoff ) then
                  rf(k) = dt/(tau*sday)*sin(0.5*pi*log(rf_cutoff/pm(k))/log(rf_cutoff/ptop))**2
                  if( is_master() ) write(6,*) k, 0.01*pm(k), dt/(rf(k)*sday)
                  kmax = k
             else
                  exit
             endif
          enddo
          RF_initialized = .true.
    endif

    allocate( u2f(isd:ied,jsd:jed,kmax) )

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd, gridstruct%nested)

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,u2f,hydrostatic,ua,va,w)
    do k=1,kmax
        if ( hydrostatic ) then
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2
              enddo
           enddo
        else
           do j=js,je
              do i=is,ie
                 u2f(i,j,k) = ua(i,j,k)**2 + va(i,j,k)**2 + w(i,j,k)**2
              enddo
           enddo
        endif
    enddo

                                        call timing_on('COMM_TOTAL')
    call mpp_update_domains(u2f, domain)
                                        call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,kmax,conserve,hydrostatic,pt,u2f,cp,rg, &
!$OMP                                  ptop,pm,rf,delz,rcv,u,v,w)
     do k=1,kmax

        if ( conserve ) then
           if ( hydrostatic ) then
             do j=js,je
                do i=is,ie
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k)/(cp-rg*ptop/pm(k))      &
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                enddo
             enddo
           else
             do j=js,je
                do i=is,ie
                   delz(i,j,k) = delz(i,j,k) / pt(i,j,k)
                   pt(i,j,k) = pt(i,j,k) + 0.5*u2f(i,j,k) * rcv      &
                             * ( 1. - 1./(1.+rf(k)*sqrt(u2f(i,j,k)/u000))**2 )
                   delz(i,j,k) = delz(i,j,k) * pt(i,j,k)
                enddo
             enddo
           endif
        endif

        do j=js-1,je+1
           do i=is-1,ie+1
              u2f(i,j,k) = rf(k)*sqrt(u2f(i,j,k)/u000)
           enddo
        enddo

        do j=js,je+1
           do i=is,ie
              u(i,j,k) = u(i,j,k) / (1.+0.5*(u2f(i,j-1,k)+u2f(i,j,k)))
           enddo
        enddo
        do j=js,je
           do i=is,ie+1
              v(i,j,k) = v(i,j,k) / (1.+0.5*(u2f(i-1,j,k)+u2f(i,j,k)))
           enddo
        enddo

        if ( .not. hydrostatic ) then
              do j=js,je
                 do i=is,ie
                    w(i,j,k) = w(i,j,k) / (1.+u2f(i,j,k))
                 enddo
              enddo
        endif

     enddo

     deallocate ( u2f )

 end subroutine Rayleigh_Friction

!>@brief The subroutine 'compute_aam' computes vertically (mass) integrated Atmospheric Angular Momentum.
 subroutine compute_aam(npz, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd,   &
                        ptop, ua, va, u, v, delp, aam, ps, m_fac)
! Compute vertically (mass) integrated Atmospheric Angular Momentum
    integer, intent(in):: npz
    integer, intent(in):: is,  ie,  js,  je
    integer, intent(in):: isd, ied, jsd, jed
    real, intent(in):: ptop
    real, intent(inout):: u(isd:ied  ,jsd:jed+1,npz) !< D grid zonal wind (m/s)
    real, intent(inout):: v(isd:ied+1,jsd:jed,npz) !< D grid meridional wind (m/s)
    real, intent(inout):: delp(isd:ied,jsd:jed,npz)
    real, intent(inout), dimension(isd:ied,jsd:jed, npz):: ua, va
    real, intent(out):: aam(is:ie,js:je)
    real, intent(out):: m_fac(is:ie,js:je)
    real, intent(out):: ps(isd:ied,jsd:jed)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN) :: gridstruct
! local:
    real, dimension(is:ie):: r1, r2, dm
    integer i, j, k

    call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd, gridstruct%nested)
    
!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,aam,m_fac,ps,ptop,delp,agrav,ua) &
!$OMP                          private(r1, r2, dm)
  do j=js,je
     do i=is,ie
        r1(i) = radius*cos(gridstruct%agrid(i,j,2))
        r2(i) = r1(i)*r1(i)
        aam(i,j) = 0.
        m_fac(i,j) = 0.
        ps(i,j) = ptop
     enddo
     do k=1,npz
        do i=is,ie
           dm(i) = delp(i,j,k)
           ps(i,j) = ps(i,j) + dm(i)
           dm(i) = dm(i)*agrav
           aam(i,j) = aam(i,j) + (r2(i)*omega + r1(i)*ua(i,j,k)) * dm(i)
           m_fac(i,j) = m_fac(i,j) + dm(i)*r2(i)
        enddo
     enddo
  enddo

 end subroutine compute_aam

!===============================================================================
  subroutine calc_tot_energy_dynamics1(gridstruct, flagstruct,domain,bd,npx,npy,npz,ncnst,ps, &
                                       ptop,delp,u,v,ua_in,va_in,q,pt,phis,outfld_name_suffix, &
                                       qsize,qsize_condensate_loading,qsize_condensate_loading_idx,qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
                                       se, ke,wv,wl,wi,wr,ws,wg,tt,mo,mr,gravit, cpair, rearth,omega_cam,ncnt,last_step,fv3_lcp_moist,fv3_lcv_moist)
    use constants_mod,          only: grav, radius, cp_air, omega
    use field_manager_mod,      only: MODEL_ATMOS
    use fv_mp_mod,              only: ng
    use fv_grid_utils_mod,      only:  c2l_ord2,g_sum
    !------------------------------Arguments--------------------------------
    
    type(fv_grid_type),  intent(IN), target                          :: gridstruct
    type(fv_flags_type), intent(IN)                                  :: flagstruct
    type(domain2d), intent(inout)                                    :: domain
    type(fv_grid_bounds_type), intent(IN)                            :: bd
    integer, intent(in)                                              :: npx,npy,npz,ncnt
    logical, intent(in)                                              :: last_step
    logical, intent(in)                                              :: fv3_lcp_moist
    logical, intent(in)                                              :: fv3_lcv_moist
    integer, intent(in)                                              :: ncnst
    real, intent(in)                                                 :: gravit,cpair,rearth,omega_cam,ptop
    real, intent(in), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed)        :: phis,ps
    real, intent(in), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)   :: pt,delp
    real, intent(inout), dimension(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,npz) :: u !< D grid zonal wind (m/s)
    real, intent(inout), dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,npz) :: v !< D grid meridional wind (m/s)
    real, intent(in) :: q(   bd%isd:bd%ied  ,bd%jsd:bd%jed  ,npz, ncnst) !< specific humidity and constituents
    real, intent(inout), dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)   :: ua_in, va_in
    real, intent(inout), dimension(bd%is:bd%ie ,bd%js:bd%je)           :: se,ke,wv,wl,wi,wr,ws,wg,tt,mo,mr
    character*(*)    , intent(in)                                    :: outfld_name_suffix ! suffix for "outfld" names
    integer, intent(IN)                                              :: qsize
    integer, intent(IN)                                              :: qsize_condensate_loading
    integer, intent(IN),dimension(qsize_condensate_loading)          :: qsize_condensate_loading_idx
    integer, intent(IN),dimension(qsize)                             :: qsize_tracer_idx_cam2dyn
    real,    intent(IN),dimension(qsize_condensate_loading)          :: qsize_condensate_loading_cp
    real,    intent(IN),dimension(qsize_condensate_loading)          :: qsize_condensate_loading_cv
    real, dimension(bd%isd:bd%ied ,bd%jsd:bd%jed )                   :: areasqrad
    real, dimension(bd%is:bd%ie ,bd%js:bd%je)                        :: senew,kenew,wvnew,wlnew,winew, &
                                                                        wrnew,wsnew,wgnew,ttnew,monew,mrnew
    real, dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)  :: dp
    real, dimension(bd%isd:bd%ied ,bd%jsd:bd%jed)       :: ps_local
    real, dimension(bd%isd:bd%ied ,bd%jsd:bd%jed ,npz)   :: ua, va
    !---------------------------Local storage-------------------------------
    
    real :: se_tmp
    real :: ke_tmp
    real :: wv_tmp
    real :: wl_tmp
    real :: wi_tmp
    real :: wr_tmp
    real :: ws_tmp
    real :: wg_tmp
    real :: tt_tmp

    real :: pdel
    real :: rncnt
    !
    ! global axial angular momentum (AAM) can be separated into one part (mr) associatedwith the relative motion 
    ! of the atmosphere with respect to the planets surface (also known as wind AAM) and another part (mo) 
    ! associated with the angular velocity OMEGA (2*pi/d, where d is the length of the day) of the planet 
    ! (also known as mass AAM)
    !
    real :: mr_cnst, mo_cnst, cos_lat, mr_tmp, mo_tmp

    real :: se_glob, ke_glob, wv_glob, wl_glob, wi_glob, wr_glob, ws_glob, wg_glob, tt_glob, mr_glob, mo_glob

    integer :: i,j,k,nq,m_cnst,n,idim,m_cnst_ffsl
    integer :: ixcldice, ixcldliq, ixtt,ixcldliq_ffsl,ixcldice_ffsl,ixtt_ffsl ! CLDICE, CLDLIQ and test tracer indices
    integer :: ixrain, ixsnow, ixgraupel,ixrain_ffsl, ixsnow_ffsl, ixgraupel_ffsl
!!$    character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6

    integer :: lchnk,ncol
    integer :: is,ie,js,je,isd,ied,jsd,jed
!!jt    logical :: fv3_lcp_moist
    logical :: printglobals = .false.
    !-----------------------------------------------------------------------

!!jt    fv3_lcp_moist                = .false.
    dp = delp
    
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    se_glob = 0.
    ke_glob = 0.
    wv_glob = 0.
    wl_glob = 0.
    wi_glob = 0.
    wr_glob = 0.
    ws_glob = 0.
    wg_glob = 0.
    tt_glob = 0.
    mr_glob = 0.
    mo_glob = 0.

    senew=0.
    kenew=0.
    wvnew=0.
    wlnew=0.
    winew=0.
    wrnew=0.
    wsnew=0.
    wgnew=0.
    ttnew=0.
    mrnew=0.
    monew=0.
    ua=0.
    va=0.
    call cubed_to_latlon(u, v, ua, va, gridstruct, &
         npx, npy, npz, 1, gridstruct%grid_type, domain, gridstruct%nested, flagstruct%c2l_ord, bd)
    
!!$    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
!!$         hist_fld_active(name_out4).or.hist_fld_active(name_out5).or.hist_fld_active(name_out6).or..true.) then
       if (qsize_condensate_loading>1) then
          ixcldliq_ffsl  = get_tracer_index (MODEL_ATMOS, 'liq_wat' )
          ixcldice_ffsl  = get_tracer_index (MODEL_ATMOS, 'ice_wat' )
          ixrain_ffsl    = get_tracer_index (MODEL_ATMOS, 'rainwat' )
          ixsnow_ffsl    = get_tracer_index (MODEL_ATMOS, 'snowwat' )
          ixgraupel_ffsl = get_tracer_index (MODEL_ATMOS, 'graupel' )
       else
          ixcldliq_ffsl = -1
          ixcldice_ffsl = -1
          ixrain_ffsl = -1
          ixsnow_ffsl = -1 
          ixgraupel_ffsl = -1
       end if
       
       ixtt_ffsl = -1
       
       
       !
       ! Compute frozen static energy in 3 parts:  KE, SE, and energy associated with vapor and liquid
       !
       do k = 1, npz
          n=0
          do j=js,je
             do i = is, ie
               n=n+1
               !
               ! make energy consistent with CAM physics (only water vapor and dry air in pressure)
               !
               if ((.not.fv3_lcp_moist).and.(.not.fv3_lcv_moist).and.qsize_condensate_loading>1) then
                 do nq=2,qsize_condensate_loading
                   m_cnst_ffsl=qsize_condensate_loading_idx(nq)
!jt                   m_cnst=qsize_condensate_loading_idx(nq)
!jt                   m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
                   dp(i,j,k) = dp(i,j,k) - delp(i,j,k)*q(i,j,k,m_cnst_ffsl)
                 end do
               end if
                !
                ! kinetic energy
                !
                ke_tmp   = 0.5*(va(i,j,k)**2+ ua(i,j,k)**2)*dp(i,j,k)/gravit
                if (fv3_lcp_moist.or.fv3_lcv_moist) then
                   !
                   ! Internal energy formula including all condensates and corresponding heat capacities
                   !
                   ! Start with energy of dry air and add energy of condensates
                   dp(i,j,k) = delp(i,j,k)             
                   do nq=1,qsize_condensate_loading                 
                      m_cnst_ffsl=qsize_condensate_loading_idx(nq)
                      dp(i,j,k) = dp(i,j,k)-delp(i,j,k)*q(i,j,k,m_cnst_ffsl)
                   end do
                   se_tmp = cpair*dp(i,j,k)
                   do nq=1,qsize_condensate_loading
                      m_cnst_ffsl=qsize_condensate_loading_idx(nq)
                      if (fv3_lcv_moist) then
                         se_tmp = se_tmp+qsize_condensate_loading_cv(nq)*q(i,j,k,m_cnst_ffsl)*delp(i,j,k)
                      end if
                      if (fv3_lcp_moist) then
                         se_tmp = se_tmp+qsize_condensate_loading_cp(nq)*q(i,j,k,m_cnst_ffsl)*delp(i,j,k)
                   end if
                   end do
                   se_tmp = se_tmp*pt(i,j,k)/gravit
                   ! reset dp to delp to use for ps_local calc below
                   dp(i,j,k) = delp(i,j,k)             
                else
                   !
                   ! using CAM physics definition of internal energy
                   !
                  se_tmp   = cpair*pt(i,j,k)*dp(i,j,k)/gravit
                end if
                wv_tmp   =  q(i,j,k,1)*delp(i,j,k)/gravit
                
                senew(i,j) = senew(i,j) + se_tmp
                kenew(i,j) = kenew(i,j) + ke_tmp
                wvnew(i,j) = wvnew(i,j) + wv_tmp

             end do
          end do
       end do
       
       do j=js,je
          do i = is,ie
             ps_local(i,j) =   ptop+sum(dp(i,j,:))
          end do
       end do

       do j=js,je
          do i = is,ie
            senew(i,j) = senew(i,j) + phis(i,j)*ps_local(i,j)/gravit
          end do
       end do
       
       ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
       
       if (ixcldliq_ffsl > 1) then
!!$          ixcldliq_ffsl = qsize_tracer_idx_cam2dyn(ixcldliq)
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   wl_tmp   = q(i,j,k,ixcldliq_ffsl)*delp(i,j,k)/gravit
                   wlnew   (i,j) = wlnew(i,j) + wl_tmp
                end do
             end do
          end do
       end if
       
       if (ixcldice_ffsl > 1) then
!!$          ixcldice_ffsl = qsize_tracer_idx_cam2dyn(ixcldice)
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   wi_tmp   = q(i,j,k,ixcldice_ffsl)*delp(i,j,k)/gravit
                   winew(i,j)    = winew(i,j) + wi_tmp
                end do
             end do
          end do
       end if

       if (ixrain_ffsl > 1) then
!!$          ixrain_ffsl = qsize_tracer_idx_cam2dyn(ixrain)
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   wr_tmp   = q(i,j,k,ixrain_ffsl)*delp(i,j,k)/gravit
                   wrnew   (i,j) = wrnew(i,j) + wr_tmp
                end do
             end do
          end do
       end if
       
       if (ixsnow_ffsl > 1) then
!!$          ixsnow_ffsl = qsize_tracer_idx_cam2dyn(ixsnow)
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   ws_tmp   = q(i,j,k,ixsnow_ffsl)*delp(i,j,k)/gravit
                   wsnew(i,j)    = wsnew(i,j) + wi_tmp
                end do
             end do
          end do
       end if

       if (ixgraupel_ffsl > 1) then
!!$          ixgraupel_ffsl = qsize_tracer_idx_cam2dyn(ixgraupel)
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   wg_tmp   = q(i,j,k,ixgraupel_ffsl)*delp(i,j,k)/gravit
                   wgnew(i,j)    = wgnew(i,j) + wg_tmp
                end do
             end do
          end do
       end if

       
       if (ixtt_ffsl > 1) then
          do k = 1, npz
             do j = js, je
                do i = is, ie
                   tt_tmp   = q(i,j,k,ixtt_ffsl)*delp(i,j,k)/gravit
                   ttnew   (i,j) = ttnew(i,j) + tt_tmp
                end do
             end do
          end do
       end if
    
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
      mr_cnst = rearth**3/gravit
      mo_cnst = omega_cam*rearth**4/gravit
      do k = 1, npz
         do j=js,je
            do i = is,ie
               pdel = delp(i,j,k)
               cos_lat = cos(gridstruct%agrid_64(i,j,2))
               mr_tmp   = mr_cnst*ua(i,j,k)*pdel*cos_lat
               mo_tmp   = mo_cnst*pdel*cos_lat**2
               
               mrnew(i,j) = mrnew(i,j) + mr_tmp
               monew(i,j) = monew(i,j) + mo_tmp
            end do
         end do
      end do

      se(is:ie,js:je)=se(is:ie,js:je)+senew(is:ie,js:je)
      ke(is:ie,js:je)=ke(is:ie,js:je)+kenew(is:ie,js:je)
      wv(is:ie,js:je)=wv(is:ie,js:je)+wvnew(is:ie,js:je)
      wl(is:ie,js:je)=wl(is:ie,js:je)+wlnew(is:ie,js:je)
      wi(is:ie,js:je)=wi(is:ie,js:je)+winew(is:ie,js:je)
      wr(is:ie,js:je)=wr(is:ie,js:je)+wrnew(is:ie,js:je)
      ws(is:ie,js:je)=ws(is:ie,js:je)+wsnew(is:ie,js:je)
      wg(is:ie,js:je)=wg(is:ie,js:je)+wgnew(is:ie,js:je)
      mr(is:ie,js:je)=mr(is:ie,js:je)+mrnew(is:ie,js:je)
      mo(is:ie,js:je)=mo(is:ie,js:je)+monew(is:ie,js:je)

      if (last_step) then 
         rncnt=ncnt
         se(is:ie,js:je)=se(is:ie,js:je)/rncnt
         ke(is:ie,js:je)=ke(is:ie,js:je)/rncnt
         wv(is:ie,js:je)=wv(is:ie,js:je)/rncnt
         wl(is:ie,js:je)=wl(is:ie,js:je)/rncnt
         wi(is:ie,js:je)=wi(is:ie,js:je)/rncnt
         wr(is:ie,js:je)=wr(is:ie,js:je)/rncnt
         ws(is:ie,js:je)=ws(is:ie,js:je)/rncnt
         wg(is:ie,js:je)=wg(is:ie,js:je)/rncnt
         mr(is:ie,js:je)=mr(is:ie,js:je)/rncnt
         mo(is:ie,js:je)=mo(is:ie,js:je)/rncnt


         if (printglobals) then
         areasqrad(isd:ied,jsd:jed) = gridstruct%area_64(isd:ied,jsd:jed)/(rearth*rearth)

         se_glob=g_sum(domain, se(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         ke_glob=g_sum(domain, ke(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         wv_glob=g_sum(domain, wv(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         wl_glob=g_sum(domain, wl(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         wi_glob=g_sum(domain, wi(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         wr_glob=g_sum(domain, wr(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         ws_glob=g_sum(domain, ws(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         wg_glob=g_sum(domain, wg(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         mr_glob=g_sum(domain, mr(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         mo_glob=g_sum(domain, mo(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         if (ixtt_ffsl > 1) tt_glob=g_sum(domain, tt(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
      
         mr_glob=g_sum(domain, mr(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)
         mo_glob=g_sum(domain, mo(is:ie,js:je), is, ie, js, je, ng, gridstruct%area_64, 1)

      
         if (mpp_pe().eq.0) then
            if (fv3_lcp_moist) then
               write(6, '(a,e25.17)') 'global wet static energy se_'//trim(outfld_name_suffix)//')            = ',se_glob
            else
               write(6, '(a,e25.17)') 'global dry static energy se_'//trim(outfld_name_suffix)//')            = ',se_glob
            end if
            write(6, '(a,e25.17)') 'global kinetic energy ke_'//trim(outfld_name_suffix)//')               = ',ke_glob
            write(6, '(a,e25.17)') 'global total energy se_plus_ke_'//trim(outfld_name_suffix)//')         = ',(ke_glob+se_glob)
            write(6, '(a,e25.17)') 'global column integrated vapor wv_'//trim(outfld_name_suffix)//'       = ',wv_glob
            write(6, '(a,e25.17)') 'global column integrated liquid wl_'//trim(outfld_name_suffix)//'      = ',wl_glob
            write(6, '(a,e25.17)') 'global column integrated ice wi_'//trim(outfld_name_suffix)//'         = ',wi_glob
            write(6, '(a,e25.17)') 'global column integrated liquid rain wr_'//trim(outfld_name_suffix)//'       = ',wr_glob
            write(6, '(a,e25.17)') 'global column integrated liquid snow ws_'//trim(outfld_name_suffix)//'      = ',ws_glob
            write(6, '(a,e25.17)') 'global column integrated graupel wg_'//trim(outfld_name_suffix)//'         = ',wg_glob
            write(6, '(a,e25.17)') 'global column integrated wind AAM mr_'//trim(outfld_name_suffix)//'         = ',mr_glob
            write(6, '(a,e25.17)') 'global column integrated mass AAM mo_'//trim(outfld_name_suffix)//'         = ',mo_glob
            if (ixtt_ffsl > 1)      write(6, '(a,e25.17)') 'global column integrated test tracer tt_'//trim(outfld_name_suffix)//' = ',tt_glob
         end if
         end if
      end if
    end subroutine calc_tot_energy_dynamics1


end module fv_dynamics_mod
