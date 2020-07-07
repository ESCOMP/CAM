module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use cam_abortutils,    only: endrun
use cam_logfile,       only: iulog
use constituents,      only: pcnst
use dimensions_mod,    only: npx,npy,nlev, &
                             cnst_name_ffsl, cnst_longname_ffsl,fv3_lcp_moist,fv3_lcv_moist, &
                             qsize_tracer_idx_cam2dyn,fv3_scale_ttend
use dyn_comp,          only: dyn_export_t, dyn_import_t
use dyn_grid,          only: get_gcol_block_d,mytile
use fv_grid_utils_mod, only: g_sum
use hycoef,            only: hyam, hybm, hyai, hybi, ps0
use mpp_domains_mod,   only: mpp_update_domains, domain2D, DGRID_NE
use perf_mod,          only: t_startf, t_stopf, t_barrierf
use physconst,         only: cpair, gravit, rair, zvir, cappa, rairv
use phys_grid,         only: get_ncols_p, get_gcol_all_p, block_to_chunk_send_pters, &
                             transpose_block_to_chunk, block_to_chunk_recv_pters, &
                             chunk_to_block_send_pters, transpose_chunk_to_block, &
                             chunk_to_block_recv_pters
use physics_types,     only: physics_state, physics_tend
use ppgrid,            only: begchunk, endchunk, pcols, pver, pverp
use shr_kind_mod,      only: r8=>shr_kind_r8, i8 => shr_kind_i8
use spmd_dyn,          only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
use spmd_utils,        only: mpicom, iam, npes,masterproc

implicit none
private
public :: d_p_coupling, p_d_coupling

!=======================================================================
contains
!=======================================================================

subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)

  ! Convert the dynamics output state into the physics input state.
  ! Note that all pressures and tracer mixing ratios coming from the FV3 dycore are based on
  ! wet air mass.


  use cam_abortutils,     only: endrun
  use fv_arrays_mod,      only: fv_atmos_type
  use fv_grid_utils_mod,  only: cubed_to_latlon
  use physics_buffer,     only: physics_buffer_desc

  ! arguments
  type (dyn_export_t),  intent(inout)                               :: dyn_out    ! dynamics export
  type (physics_buffer_desc), pointer                               :: pbuf2d(:,:)
  type (physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type (physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

  ! LOCAL VARIABLES

  integer :: ib                     ! indices over elements
  integer :: ioff
  integer :: lchnk, icol, ilyr      ! indices over chunks, columns, layers
  integer :: m, m_ffsl, n, i, j, k

  integer :: cpter(pcols,              0:pver)    ! offsets into chunk buffer for unpacking data

  integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
  integer :: tsize                 ! amount of data per grid point passed to physics
  type (fv_atmos_type),  pointer :: Atm(:)

  integer                                   :: is,ie,js,je
  integer                                   :: ncols

  ! LOCAL Allocatables
  integer, allocatable,  dimension(:,:)     :: bpter    !((ie-is+1)*(je-js+1),0:pver) ! packing data block buffer offset
  real(r8), allocatable, dimension(:)       :: bbuffer, cbuffer ! transpose buffers
  real(r8), allocatable, dimension(:,:)     :: phis_tmp !((ie-is+1)*(je-js+1),     1) ! temporary array to hold phis
  real(r8), allocatable, dimension(:,:)     :: ps_tmp   !((ie-is+1)*(je-js+1),     1) ! temporary array to hold ps
  real(r8), allocatable, dimension(:,:,:)   :: T_tmp    !((ie-is+1)*(je-js+1),pver,1) ! temporary array to hold T
  real(r8), allocatable, dimension(:,:,:)   :: omega_tmp!((ie-is+1)*(je-js+1),pver,1) ! temporary array to hold omega
  real(r8), allocatable, dimension(:,:,:)   :: pdel_tmp !((ie-is+1)*(je-js+1),pver,1) ! temporary array to hold pdel
  real(r8), allocatable, dimension(:,:,:)   :: u_tmp !((ie-is+1)*(je-js+1),pver,1) ! temp array to hold u
  real(r8), allocatable, dimension(:,:,:)   :: v_tmp !((ie-is+1)*(je-js+1),pver,1) ! temp array to hold v
  real(r8), allocatable, dimension(:,:,:,:) :: q_tmp !((ie-is+1)*(je-js+1),pver,pcnst,1) ! temp to hold advected constituents

  !-----------------------------------------------------------------------

  Atm=>dyn_out%atm

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  ! Allocate temporary arrays to hold data for physics decomposition
  allocate(ps_tmp   ((ie-is+1)*(je-js+1),           1))
  allocate(phis_tmp ((ie-is+1)*(je-js+1),           1))
  allocate(T_tmp    ((ie-is+1)*(je-js+1),pver,      1))
  allocate(u_tmp    ((ie-is+1)*(je-js+1),pver,      1))
  allocate(v_tmp    ((ie-is+1)*(je-js+1),pver,      1))
  allocate(omega_tmp((ie-is+1)*(je-js+1),pver,      1))
  allocate(pdel_tmp ((ie-is+1)*(je-js+1),pver,      1))
  allocate(Q_tmp    ((ie-is+1)*(je-js+1),pver,pcnst, 1))

  ps_tmp   = 0._r8
  phis_tmp = 0._r8
  T_tmp    = 0._r8
  u_tmp    = 0._r8
  v_tmp    = 0._r8
  omega_tmp= 0._r8
  pdel_tmp = 0._r8
  Q_tmp    = 0._r8

  n = 1
  do j = js, je
     do i = is, ie
        ps_tmp  (n, 1) = Atm(mytile)%ps  (i, j)
        phis_tmp(n, 1) = Atm(mytile)%phis(i, j)
        do k = 1, pver
           T_tmp    (n, k, 1) = Atm(mytile)%pt  (i, j, k)
           u_tmp    (n, k, 1) = Atm(mytile)%ua (i, j, k)
           v_tmp    (n, k, 1) = Atm(mytile)%va (i, j, k)
           omega_tmp(n, k, 1) = Atm(mytile)%omga(i, j, k)
           pdel_tmp (n, k, 1) = Atm(mytile)%delp(i, j, k)
           !
           ! The fv3 constituent array may be in a different order than the cam array, remap here.
           !
           do m = 1, pcnst
              m_ffsl=qsize_tracer_idx_cam2dyn(m)
              Q_tmp(n, k, m, 1) = Atm(mytile)%q(i, j, k, m_ffsl)
           end do
        end do
        n = n + 1
     end do
  end do

  call t_startf('dpcopy')
  if (local_dp_map) then

     !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ib, ioff, ilyr, m)
     do lchnk = begchunk, endchunk
        ncols = get_ncols_p(lchnk)
        call get_gcol_all_p(lchnk, pcols, pgcols)
        do icol = 1, ncols
           call get_gcol_block_d(pgcols(icol), 1, idmb1, idmb2, idmb3)
           ib   = idmb3(1)
           ioff = idmb2(1)
           phys_state(lchnk)%ps(icol)   = ps_tmp  (ioff,ib)
           phys_state(lchnk)%phis(icol) = phis_tmp(ioff,ib)
           do ilyr = 1, pver
              phys_state(lchnk)%t    (icol,ilyr) = T_tmp    (ioff,ilyr,ib)
              phys_state(lchnk)%u    (icol,ilyr) = u_tmp    (ioff,ilyr,ib)
              phys_state(lchnk)%v    (icol,ilyr) = v_tmp    (ioff,ilyr,ib)
              phys_state(lchnk)%omega(icol,ilyr) = omega_tmp(ioff,ilyr,ib)
              phys_state(lchnk)%pdel(icol,ilyr)  = pdel_tmp (ioff,ilyr,ib)
              do m = 1, pcnst
                 phys_state(lchnk)%q(icol,ilyr,m) = Q_tmp(ioff,ilyr,m,ib)
              end do
           end do
        end do

     end do


  else  ! .not. local_dp_map

     tsize = 5 + pcnst
     ib = 1

     allocate(bbuffer(tsize*block_buf_nrecs))
     allocate(cbuffer(tsize*chunk_buf_nrecs))
     allocate(bpter((ie-is+1)*(je-js+1),0:pver))

     if (iam < npes) then
        call block_to_chunk_send_pters(iam+1, (ie-is+1)*(je-js+1), pver+1, tsize, bpter)
        do icol = 1, (ie-is+1)*(je-js+1)
           bbuffer(bpter(icol,0)+2:bpter(icol,0)+tsize-1) = 0.0_r8
           bbuffer(bpter(icol,0))   = ps_tmp  (icol,ib)
           bbuffer(bpter(icol,0)+1) = phis_tmp(icol,ib)
           do ilyr = 1, pver
              bbuffer(bpter(icol,ilyr))   = T_tmp(icol,ilyr,ib)
              bbuffer(bpter(icol,ilyr)+1) = u_tmp(icol,ilyr,ib)
              bbuffer(bpter(icol,ilyr)+2) = v_tmp(icol,ilyr,ib)
              bbuffer(bpter(icol,ilyr)+3) = omega_tmp(icol,ilyr,ib)
              bbuffer(bpter(icol,ilyr)+4) = pdel_tmp (icol,ilyr,ib)
              do m = 1, pcnst
                 bbuffer(bpter(icol,ilyr)+tsize-pcnst-1+m) = Q_tmp(icol,ilyr,m,ib)
              end do
           end do
        end do
     else
        bbuffer(:) = 0._r8
     end if

     call t_barrierf ('sync_blk_to_chk', mpicom)
     call t_startf ('block_to_chunk')
     call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
     call t_stopf  ('block_to_chunk')

     do lchnk = begchunk,endchunk
        ncols = phys_state(lchnk)%ncol
        call block_to_chunk_recv_pters(lchnk, pcols, pver+1, tsize, cpter)
        do icol = 1, ncols
           phys_state(lchnk)%ps   (icol) = cbuffer(cpter(icol,0))
           phys_state(lchnk)%phis (icol) = cbuffer(cpter(icol,0)+1)
           do ilyr = 1, pver
              phys_state(lchnk)%t     (icol,ilyr)   = cbuffer(cpter(icol,ilyr))
              phys_state(lchnk)%u     (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+1)
              phys_state(lchnk)%v     (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+2)
              phys_state(lchnk)%omega (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+3)
              phys_state(lchnk)%pdel  (icol,ilyr)   = cbuffer(cpter(icol,ilyr)+4)
              do m = 1, pcnst
                 phys_state(lchnk)%q (icol,ilyr,m) = cbuffer(cpter(icol,ilyr)+tsize-pcnst-1+m)
              end do
           end do
        end do
     end do

     deallocate( bbuffer )
     deallocate( cbuffer )
     deallocate( bpter )

  end if

  deallocate(ps_tmp   )
  deallocate(phis_tmp )
  deallocate(T_tmp    )
  deallocate(u_tmp    )
  deallocate(v_tmp    )
  deallocate(omega_tmp)
  deallocate(pdel_tmp )
  deallocate(Q_tmp    )

  call t_stopf('dpcopy')

  ! derive the physics state from the dynamics state converting to proper vapor loading
  ! and setting dry mixing ratio variables based on cnst_type - no need to call wet_to_dry
  ! since derived_phys_dry takes care of that.

  call t_startf('derived_phys_dry')
  call derived_phys_dry(phys_state, phys_tend, pbuf2d)
  call t_stopf('derived_phys_dry')

end subroutine d_p_coupling

!=======================================================================

subroutine p_d_coupling(phys_state, phys_tend, dyn_in)

  ! Convert the physics output state into the dynamics input state.

  use cam_history,            only: outfld
  use constants_mod,          only: cp_air, kappa
  use dyn_comp,               only: calc_tot_energy_dynamics
  use fms_mod,                only: set_domain
  use fv_arrays_mod,          only: fv_atmos_type
  use fv_grid_utils_mod,      only: cubed_to_latlon
  use physconst,              only: thermodynamic_active_species_num,thermodynamic_active_species_idx_dycore
  use physconst,              only: thermodynamic_active_species_cp,thermodynamic_active_species_cv,dry_air_species_num
  use physics_types,          only: set_state_pdry
  use time_manager,           only: get_step_size

  ! arguments
  type (physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type (physics_tend),  intent(inout), dimension(begchunk:endchunk) :: phys_tend
  type (dyn_import_t),  intent(inout) :: dyn_in

  ! LOCAL VARIABLES

  integer :: cpter(pcols,0:pver)    ! offsets into chunk buffer for unpacking data
  integer :: ib                     ! indices over elements
  integer :: idim
  integer :: ioff
  integer :: is,isd,ie,ied,js,jsd,je,jed
  integer :: lchnk, icol, ilyr      ! indices over chunks, columns, layers
  integer :: m, n, i, j, k,m_ffsl,nq
  integer :: ncols
  integer :: pgcols(pcols), idmb1(1), idmb2(1), idmb3(1)
  integer :: tsize                 ! amount of data per grid point passed to physics
  integer :: num_wet_species       ! total number of wet species (first tracers in FV3 tracer array)

  integer, allocatable, dimension(:,:) :: bpter   !((ie-is+1)*(je-js+1),0:pver)    ! packing data block buffer offsets
  real(r8),  allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

  real (r8)                            :: dt
  real (r8)                            :: fv3_totwatermass, fv3_airmass
  real (r8)                            :: qall,cpfv3
  real (r8)                            :: tracermass(pcnst)

  type (fv_atmos_type),  pointer :: Atm(:)

  real(r8),  allocatable, dimension(:,:,:)   :: delpdry       ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: pdel_tmp      ! temporary to hold
  real(r8),  allocatable, dimension(:,:,:)   :: pdeldry_tmp   ! temporary to hold
  real(r8),  allocatable, dimension(:,:,:)   :: t_dt          ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: t_dt_tmp      ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: t_tendadj     ! temporary array to temperature tendency adjustment
  real(r8),  allocatable, dimension(:,:,:)   :: u_dt          ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: u_dt_tmp      ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: u_tmp         ! temporary array to hold u and v
  real(r8),  allocatable, dimension(:,:,:)   :: v_dt          ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: v_dt_tmp      ! temporary to hold tendencies
  real(r8),  allocatable, dimension(:,:,:)   :: v_tmp         ! temporary array to hold u and v
  real(r8),  allocatable, dimension(:,:,:,:) :: q_tmp         ! temporary to hold

  !-----------------------------------------------------------------------

  Atm=>dyn_in%atm

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed

  call set_domain ( Atm(mytile)%domain )

  allocate(delpdry(isd:ied,jsd:jed,nlev))
  allocate(t_dt_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(u_dt_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(v_dt_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(pdel_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(pdeldry_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(U_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(V_tmp((ie-is+1)*(je-js+1),pver,1))
  allocate(Q_tmp((ie-is+1)*(je-js+1),pver,pcnst,1))
  allocate(u_dt(isd:ied,jsd:jed,nlev))
  allocate(v_dt(isd:ied,jsd:jed,nlev))
  allocate(t_dt(is:ie,js:je,nlev))
  allocate(t_tendadj(is:ie,js:je,nlev))

  Atm=>dyn_in%atm

  if (local_dp_map) then
!$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ib, ioff, ilyr, m)
     do lchnk = begchunk, endchunk
        ncols = get_ncols_p(lchnk)
        call get_gcol_all_p(lchnk, pcols, pgcols)
        call set_state_pdry(phys_state(lchnk)) ! First get dry pressure to use for this timestep
        do icol = 1, ncols
           call get_gcol_block_d(pgcols(icol), 1, idmb1, idmb2, idmb3)
           ib   = idmb3(1)
           ioff = idmb2(1)
           do ilyr = 1, pver
              t_dt_tmp(ioff,ilyr,ib) = phys_tend(lchnk)%dtdt(icol,ilyr)
              u_tmp(ioff,ilyr,ib) = phys_state(lchnk)%u(icol,ilyr)
              v_tmp(ioff,ilyr,ib) = phys_state(lchnk)%v(icol,ilyr)
              u_dt_tmp(ioff,ilyr,ib) = phys_tend(lchnk)%dudt(icol,ilyr)
              v_dt_tmp(ioff,ilyr,ib) = phys_tend(lchnk)%dvdt(icol,ilyr)
              pdel_tmp(ioff,ilyr,ib) = phys_state(lchnk)%pdel(icol,ilyr)
              pdeldry_tmp(ioff,ilyr,ib) = phys_state(lchnk)%pdeldry(icol,ilyr)
              do m=1, pcnst
                 Q_tmp(ioff,ilyr,m,ib) = phys_state(lchnk)%q(icol,ilyr,m)
              end do
           end do
        end do
     end do

  else

     tsize = 7 + pcnst
     ib = 1

     allocate(bbuffer(tsize*block_buf_nrecs))
     allocate(cbuffer(tsize*chunk_buf_nrecs))
     allocate(bpter((ie-is+1)*(je-js+1),0:pver))    ! offsets into block buffer for packing data

!$omp parallel do private (lchnk, ncols, cpter, i, icol, ilyr, m)
     do lchnk = begchunk, endchunk

        call set_state_pdry(phys_state(lchnk))	 ! First get dry pressure to use for this timestep
        ncols = get_ncols_p(lchnk)

        call chunk_to_block_send_pters(lchnk, pcols, pver+1, tsize, cpter)

        do i=1,ncols
           cbuffer(cpter(i,0):cpter(i,0)+6+pcnst) = 0.0_r8
        end do

        do icol = 1, ncols

           do ilyr = 1, pver
              cbuffer(cpter(icol,ilyr))   = phys_tend(lchnk)%dtdt(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+1) = phys_state(lchnk)%u(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+2) = phys_state(lchnk)%v(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+3) = phys_tend(lchnk)%dudt(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+4) = phys_tend(lchnk)%dvdt(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+5) = phys_state(lchnk)%pdel(icol,ilyr)
              cbuffer(cpter(icol,ilyr)+6) = phys_state(lchnk)%pdeldry(icol,ilyr)
              do m = 1, pcnst
                 cbuffer(cpter(icol,ilyr)+6+m) = phys_state(lchnk)%q(icol,ilyr,m)
              end do
           end do

        end do

     end do

     call t_barrierf('sync_chk_to_blk', mpicom)
     call t_startf ('chunk_to_block')
     call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
     call t_stopf  ('chunk_to_block')

     if (iam < npes) then

        call chunk_to_block_recv_pters(iam+1, (ie-is+1)*(je-js+1), pver+1, tsize, bpter)
        do icol = 1, (ie-is+1)*(je-js+1)
           do ilyr = 1, pver
              t_dt_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr))
              u_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+1)
              v_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+2)
              u_dt_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+3)
              v_dt_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+4)
              pdel_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+5)
              pdeldry_tmp(icol,ilyr,ib) = bbuffer(bpter(icol,ilyr)+6)
              do m = 1, pcnst
                 Q_tmp(icol,ilyr,m,ib) = bbuffer(bpter(icol,ilyr)+6+m)
              end do
           end do
        end do

     end if

     deallocate(bbuffer)
     deallocate(cbuffer)
     deallocate(bpter)

  end if

  dt = get_step_size()

  idim=ie-is+1

! pt_dt is adjusted below.
  n = 1
  do j = js, je
     do i = is, ie
        do k = 1, pver
           t_dt(i, j, k) = t_dt_tmp    (n, k, 1)
           u_dt(i, j, k) = u_dt_tmp    (n, k, 1)
           v_dt(i, j, k) = v_dt_tmp    (n, k, 1)
           Atm(mytile)%ua(i, j, k) =  Atm(mytile)%ua(i, j, k) + u_dt(i, j, k)*dt
           Atm(mytile)%va(i, j, k) =  Atm(mytile)%va(i, j, k) + v_dt(i, j, k)*dt
           Atm(mytile)%delp(i, j, k) = pdel_tmp    (n, k, 1)
           delpdry(i, j, k) = pdeldry_tmp    (n, k, 1)
           do m = 1, pcnst
              ! dynamics tracers may be in a different order from cam tracer array
              m_ffsl=qsize_tracer_idx_cam2dyn(m)
              Atm(mytile)%q(i, j, k, m_ffsl) = Q_tmp(n, k, m, 1)
           end do
        end do
        n = n + 1
     end do
   end do

  ! Update delp and mixing ratios to account for the difference between CAM and FV3 total air mass
  ! CAM total air mass (pdel)  = (dry + vapor)
  ! FV3 total air mass (delp at beg of phys * mix ratio) =
  ! drymass + (vapor + condensate [liq_wat,ice_wat,rainwat,snowwat,graupel])*mix ratio
  ! FV3 tracer mixing ratios = tracer mass / FV3 total air mass
  ! convert the (dry+vap) mixing ratios to be based off of FV3 condensate loaded airmass (dry+vap+cond). When
  ! d_p_coupling/derive_phys_dry is called the mixing ratios are again parsed out into wet and
  ! dry for physics.
  num_wet_species=thermodynamic_active_species_num-dry_air_species_num
   ! recalculate ps based on new delp
  Atm(mytile)%ps(:,:)=hyai(1)*ps0
  do k=1,pver
     do j = js,je
        do i = is,ie
           do m = 1,pcnst
              tracermass(m)=Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m)
           end do
           fv3_totwatermass=sum(tracermass(thermodynamic_active_species_idx_dycore(1:num_wet_species)))
           fv3_airmass =  delpdry(i,j,k) + fv3_totwatermass
           Atm(mytile)%delp(i,j,k) = fv3_airmass
           Atm(mytile)%q(i,j,k,1:pcnst) = tracermass(1:pcnst)/fv3_airmass
           Atm(mytile)%ps(i,j)=Atm(mytile)%ps(i,j)+Atm(mytile)%delp(i, j, k)
        end do
     end do
  end do

   ! update dynamics temperature from physics tendency
   ! if using fv3_lcv_moist adjust temperature tendency to conserve energy across phys/dynamics
   ! interface accounting for differences in the moist/wet assumptions

   do k = 1, pver
      do j = js, je
         do i = is, ie
           if (fv3_scale_ttend) then
              qall=0._r8
              cpfv3=0._r8
              do nq=1,thermodynamic_active_species_num
                 m_ffsl = thermodynamic_active_species_idx_dycore(nq)
                 qall=qall+Atm(mytile)%q(i,j,k,m_ffsl)
                 if (fv3_lcp_moist) cpfv3 = cpfv3+thermodynamic_active_species_cp(nq)*Atm(mytile)%q(i,j,k,m_ffsl)
                 if (fv3_lcv_moist) cpfv3 = cpfv3+thermodynamic_active_species_cv(nq)*Atm(mytile)%q(i,j,k,m_ffsl)
              end do
              cpfv3=(1._r8-qall)*cp_air+cpfv3
              ! scale factor for t_dt so temperature tendency derived from CAM moist air (dry+vap - constant pressure)
              ! can be applied to FV3 wet air (dry+vap+cond - constant volume)

              t_tendadj(i,j,k)=cp_air/cpfv3

              if (.not.Atm(mytile)%flagstruct%hydrostatic) then
                 ! update to nonhydrostatic variable delz to account for phys temperature adjustment.
                 Atm(mytile)%delz(i, j, k) = Atm(mytile)%delz(i,j,k)/Atm(mytile)%pt(i, j, k)
                 Atm(mytile)%pt (i, j, k) = Atm(mytile)%pt (i, j, k) + t_dt(i, j, k)*dt*t_tendadj(i,j,k)
                 Atm(mytile)%delz(i, j, k) = Atm(mytile)%delz(i,j,k)*Atm(mytile)%pt (i, j, k)
              else
                 Atm(mytile)%pt (i, j, k) = Atm(mytile)%pt (i, j, k) + t_dt(i, j, k)*dt*t_tendadj(i,j,k)
              end if
           else
              Atm(mytile)%pt (i, j, k) = Atm(mytile)%pt (i, j, k) + t_dt(i, j, k)*dt
           end if
        end do
     end do
  end do

  !$omp parallel do private(i, j)
  do j=js,je
     do i=is,ie
        Atm(mytile)%pe(i,1,j)   = Atm(mytile)%ptop
        Atm(mytile)%pk(i,j,1)   = Atm(mytile)%ptop ** kappa
        Atm(mytile)%peln(i,1,j) = log(Atm(mytile)%ptop )
     enddo
  enddo

!$omp parallel do private(i,j,k)
  do j=js,je
     do k=1,pver
        do i=is,ie
           Atm(mytile)%pe(i,k+1,j)   = Atm(mytile)%pe(i,k,j) + Atm(mytile)%delp(i,j,k)
        enddo
     enddo
  enddo

!$omp parallel do private(i,j,k)
  do j=js,je
     do k=1,pver
        do i=is,ie
           Atm(mytile)%pk(i,j,k+1)= Atm(mytile)%pe(i,k+1,j) ** kappa
           Atm(mytile)%peln(i,k+1,j) = log(Atm(mytile)%pe(i,k+1,j))
           Atm(mytile)%pkz(i,j,k) = (Atm(mytile)%pk(i,j,k+1)-Atm(mytile)%pk(i,j,k))/ &
                                    (kappa*(Atm(mytile)%peln(i,k+1,j)-Atm(mytile)%peln(i,k,j)))
        enddo
     enddo
  enddo

  do j = js, je
     call outfld('FU', RESHAPE(u_dt(is:ie, j, :),(/idim,pver/)), idim, j)
     call outfld('FV', RESHAPE(v_dt(is:ie, j, :),(/idim,pver/)), idim, j)
     call outfld('FT', RESHAPE(t_dt(is:ie, j, :),(/idim,pver/)), idim, j)
  end do

  call calc_tot_energy_dynamics(dyn_in%atm,'dAP')


  !set the D-Grid winds from the physics A-grid winds/tendencies.
  if ( Atm(mytile)%flagstruct%dwind_2d ) then
     call endrun('dwind_2d update is not implemented')
  else
     call atend2dstate3d(  u_dt, v_dt, Atm(mytile)%u ,Atm(mytile)%v, is,  ie,  js,  je, &
                           isd, ied, jsd, jed, npx,npy, nlev, Atm(mytile)%gridstruct, Atm(mytile)%domain, dt)
  endif

  ! Again we are rederiving the A winds from the Dwinds to give our energy dynamics a consistent wind.
  call cubed_to_latlon(Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%gridstruct, &
                       npx, npy, nlev, 1, Atm(mytile)%gridstruct%grid_type, Atm(mytile)%domain, &
                       Atm(mytile)%gridstruct%nested, Atm(mytile)%flagstruct%c2l_ord, Atm(mytile)%bd)

  !$omp parallel do private(i, j)
  do j=js,je
     do i=is,ie
        Atm(mytile)%u_srf=Atm(mytile)%ua(i,j,pver)
        Atm(mytile)%v_srf=Atm(mytile)%va(i,j,pver)
     enddo
  enddo

  ! update halo regions
  call mpp_update_domains( Atm(mytile)%delp, Atm(mytile)%domain )
  call mpp_update_domains( Atm(mytile)%ps, Atm(mytile)%domain )
  call mpp_update_domains( Atm(mytile)%phis, Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%ps,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,   Atm(mytile)%domain, gridtype=DGRID_NE, complete=.true. )
  call mpp_update_domains( atm(mytile)%pt,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%q,    Atm(mytile)%domain )

  deallocate(delpdry)
  deallocate(t_dt_tmp)
  deallocate(u_dt_tmp)
  deallocate(v_dt_tmp)
  deallocate(pdel_tmp)
  deallocate(pdeldry_tmp)
  deallocate(U_tmp)
  deallocate(V_tmp)
  deallocate(Q_tmp)
  deallocate(u_dt)
  deallocate(v_dt)
  deallocate(t_dt)
  deallocate(t_tendadj)

end subroutine p_d_coupling

!=======================================================================

subroutine derived_phys_dry(phys_state, phys_tend, pbuf2d)

  use check_energy,   only: check_energy_timestep_init
  use constituents,   only: qmin
  use geopotential,   only: geopotential_t
  use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
  use physics_types,  only: set_wet_to_dry
  use physconst,      only: thermodynamic_active_species_num,thermodynamic_active_species_idx_dycore
  use physconst,      only: thermodynamic_active_species_idx,dry_air_species_num
  use ppgrid,         only: pver
  use qneg_module,    only: qneg3
  use shr_vmath_mod,  only: shr_vmath_log

  ! arguments
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
  type(physics_buffer_desc),      pointer     :: pbuf2d(:,:)

  ! local variables

  integer                         :: num_wet_species       ! total number of wet species (first tracers in FV3 tracer array)
  integer                         :: lchnk
  integer                         :: m, i, k, ncol

  real(r8)                        :: cam_totwatermass, cam_airmass
  real(r8), dimension(pcnst)      :: tracermass
  real(r8), dimension(pcols,pver) :: zvirv    ! Local zvir array pointer

  !----------------------------------------------------------------------------

  type(physics_buffer_desc), pointer :: pbuf_chnk(:)

  !
  ! Evaluate derived quantities
  !
  ! At this point the phys_state has been filled in from dynamics, rearranging tracers to match CAM tracer order.
  ! pdel is consistent with tracer array.
  ! All tracer mixing rations at this point are calculated using dry+vap+condensates - we need to convert
  ! to cam physics wet mixing ration based off of dry+vap.
  ! Following this loop call wet_to_dry to convert CAM's dry constituents to their dry mixing ratio.

!!! omp parallel do private (lchnk, ncol, k, i, zvirv, pbuf_chnk,m,cam_airmass,cam_totwatermass)
  num_wet_species=thermodynamic_active_species_num-dry_air_species_num
  do lchnk = begchunk,endchunk
     ncol = get_ncols_p(lchnk)
     do k=1,pver
        do i=1,ncol
           phys_state(lchnk)%pdeldry(i,k) = &
                phys_state(lchnk)%pdel(i,k) * &
                (1._r8-sum(phys_state(lchnk)%q(i,k,thermodynamic_active_species_idx(1:num_wet_species))))
           do m = 1,pcnst
              tracermass(m)=phys_state(lchnk)%pdel(i,k)*phys_state(lchnk)%q(i,k,m)
           end do
           cam_totwatermass=tracermass(1)
           cam_airmass =  phys_state(lchnk)%pdeldry(i,k) + cam_totwatermass
           phys_state(lchnk)%pdel(i,k) = cam_airmass
           phys_state(lchnk)%q(i,k,1:pcnst) = tracermass(1:pcnst)/cam_airmass
        end do
     end do

! Physics state now has CAM pdel (dry+vap) and pdeldry and all constituents are dry+vap
! Convert dry type constituents from moist to dry mixing ratio
!
     call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.

!
! Derive the rest of the pressure variables using pdel and pdeldry
!

     do i = 1, ncol
        phys_state(lchnk)%psdry(i) = hyai(1)*ps0 + sum(phys_state(lchnk)%pdeldry(i,:))
     end do

     do i = 1, ncol
        phys_state(lchnk)%pintdry(i,1) = hyai(1)*ps0
     end do
     call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,1), &
          phys_state(lchnk)%lnpintdry(1:ncol,1),ncol)
     do k = 1, pver
        do i = 1, ncol
           phys_state(lchnk)%pintdry(i,k+1) = phys_state(lchnk)%pintdry(i,k) + &
                phys_state(lchnk)%pdeldry(i,k)
        end do
        call shr_vmath_log(phys_state(lchnk)%pintdry(1:ncol,k+1),&
             phys_state(lchnk)%lnpintdry(1:ncol,k+1),ncol)
     end do

     do k=1,pver
        do i=1,ncol
           phys_state(lchnk)%rpdeldry(i,k) = 1._r8/phys_state(lchnk)%pdeldry(i,k)
           phys_state(lchnk)%pmiddry (i,k) = 0.5_r8*(phys_state(lchnk)%pintdry(i,k+1) + &
                phys_state(lchnk)%pintdry(i,k))
        end do
        call shr_vmath_log(phys_state(lchnk)%pmiddry(1:ncol,k), &
             phys_state(lchnk)%lnpmiddry(1:ncol,k),ncol)
     end do

     ! initialize moist pressure variables

     do i=1,ncol
        phys_state(lchnk)%ps(i)     = phys_state(lchnk)%pintdry(i,1)
        phys_state(lchnk)%pint(i,1) = phys_state(lchnk)%pintdry(i,1)
     end do
     do k = 1, pver
        do i=1,ncol
           phys_state(lchnk)%pint(i,k+1) =  phys_state(lchnk)%pint(i,k)+phys_state(lchnk)%pdel(i,k)
           phys_state(lchnk)%pmid(i,k)   = (phys_state(lchnk)%pint(i,k+1)+phys_state(lchnk)%pint(i,k))/2._r8
           phys_state(lchnk)%ps  (i)     =  phys_state(lchnk)%ps(i) + phys_state(lchnk)%pdel(i,k)
        end do
        call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
        call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
     end do
     call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

     do k = 1, pver
        do i = 1, ncol
           phys_state(lchnk)%rpdel(i,k)  = 1._r8/phys_state(lchnk)%pdel(i,k)
           phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
                / phys_state(lchnk)%pmid(i,k))**cappa
        end do
     end do

     ! fill zvirv 2D variables to be compatible with geopotential_t interface
     zvirv(:,:) = zvir

     ! Compute initial geopotential heights - based on full pressure
     call geopotential_t (phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
          phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
          phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk),  gravit,  zvirv       , &
          phys_state(lchnk)%zi    , phys_state(lchnk)%zm      , ncol                )

     ! Compute initial dry static energy, include surface geopotential
     do k = 1, pver
        do i = 1, ncol
           phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
        end do
     end do
     ! Ensure tracers are all positive
     call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
          1, pcnst, qmin  ,phys_state(lchnk)%q)

     ! Compute energy and water integrals of input state
     pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
     call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

  end do  ! lchnk

end subroutine derived_phys_dry

subroutine atend2dstate3d(u_dt, v_dt, u, v, is,  ie,  js,  je, isd, ied, jsd, jed, npx,npy, nlev, gridstruct, domain, dt)
!----------------------------------------------------------------------------
! This routine adds the a-grid wind tendencies returned by the physics to the d-state
! wind being sent to the dynamics.
!----------------------------------------------------------------------------

  use fv_arrays_mod,      only: fv_grid_type
  use mpp_domains_mod,    only: mpp_update_domains,  DGRID_NE

  ! arguments
  integer, intent(in)  :: npx,npy, nlev
  integer, intent(in)  :: is,  ie,  js,  je,&
                          isd, ied, jsd, jed
  real(r8), intent(in) :: dt
  real(r8), intent(inout), dimension(isd:ied,jsd:jed,nlev)     :: u_dt, v_dt
  real(r8), intent(inout), dimension(isd:ied,  jsd:jed+1,nlev) :: u
  real(r8), intent(inout), dimension(isd:ied+1,jsd:jed  ,nlev) :: v
  type(domain2d), intent(inout) :: domain
  type(fv_grid_type), intent(in), target :: gridstruct

  ! local:

  integer i, j, k, im2, jm2
  real(r8) dt5
  real(r8), dimension(is-1:ie+1,js:je+1,3)   :: ue    ! 3D winds at edges
  real(r8), dimension(is-1:ie+1,js-1:je+1,3) :: v3
  real(r8), dimension(is:ie+1,js-1:je+1,  3) :: ve    ! 3D winds at edges
  real(r8), dimension(is:ie)                 :: ut1, ut2, ut3
  real(r8), dimension(js:je)                 :: vt1, vt2, vt3
  real(r8), pointer, dimension(:)            :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n
  real(r8), pointer, dimension(:,:,:)        :: vlon, vlat
  real(r8), pointer, dimension(:,:,:,:)      :: es, ew

  !----------------------------------------------------------------------------

  es   => gridstruct%es
  ew   => gridstruct%ew
  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  edge_vect_w => gridstruct%edge_vect_w
  edge_vect_e => gridstruct%edge_vect_e
  edge_vect_s => gridstruct%edge_vect_s
  edge_vect_n => gridstruct%edge_vect_n

  call mpp_update_domains(u_dt, domain, complete=.false.)
  call mpp_update_domains(v_dt, domain, complete=.true.)

  dt5 = 0.5_r8 * dt
  im2 = (npx-1)/2
  jm2 = (npy-1)/2

!$OMP parallel do default(none) shared(is,ie,js,je,nlev,gridstruct,u,dt5,u_dt,v,v_dt,  &
!$OMP                                  vlon,vlat,jm2,edge_vect_w,npx,edge_vect_e,im2, &
!$OMP                                  edge_vect_s,npy,edge_vect_n,es,ew)             &
!$OMP                          private(ut1, ut2, ut3, vt1, vt2, vt3, ue, ve, v3)
  do k=1, nlev

     ! Compute 3D wind/tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = u_dt(i,j,k)*vlon(i,j,1) + v_dt(i,j,k)*vlat(i,j,1)
             v3(i,j,2) = u_dt(i,j,k)*vlon(i,j,2) + v_dt(i,j,k)*vlat(i,j,2)
             v3(i,j,3) = u_dt(i,j,k)*vlon(i,j,3) + v_dt(i,j,k)*vlat(i,j,3)
          enddo
       enddo

       ! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

       ! --- E_W edges (for v-wind):
       if (.not. gridstruct%nested) then
          if ( is==1) then
             i = 1
             do j=js,je
                if ( j>jm2 ) then
                   vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1._r8-edge_vect_w(j))*ve(i,j,1)
                   vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1._r8-edge_vect_w(j))*ve(i,j,2)
                   vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1._r8-edge_vect_w(j))*ve(i,j,3)
                else
                   vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1._r8-edge_vect_w(j))*ve(i,j,1)
                   vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1._r8-edge_vect_w(j))*ve(i,j,2)
                   vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1._r8-edge_vect_w(j))*ve(i,j,3)
                endif
             enddo
             do j=js,je
                ve(i,j,1) = vt1(j)
                ve(i,j,2) = vt2(j)
                ve(i,j,3) = vt3(j)
             enddo
          endif

          if ( (ie+1)==npx ) then
             i = npx
             do j=js,je
                if ( j>jm2 ) then
                   vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1._r8-edge_vect_e(j))*ve(i,j,1)
                   vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1._r8-edge_vect_e(j))*ve(i,j,2)
                   vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1._r8-edge_vect_e(j))*ve(i,j,3)
                else
                   vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1._r8-edge_vect_e(j))*ve(i,j,1)
                   vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1._r8-edge_vect_e(j))*ve(i,j,2)
                   vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1._r8-edge_vect_e(j))*ve(i,j,3)
                endif
             enddo
             do j=js,je
                ve(i,j,1) = vt1(j)
                ve(i,j,2) = vt2(j)
                ve(i,j,3) = vt3(j)
             enddo
          endif
          ! N-S edges (for u-wind):
          if ( js==1) then
             j = 1
             do i=is,ie
                if ( i>im2 ) then
                   ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1._r8-edge_vect_s(i))*ue(i,j,1)
                   ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1._r8-edge_vect_s(i))*ue(i,j,2)
                   ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1._r8-edge_vect_s(i))*ue(i,j,3)
                else
                   ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1._r8-edge_vect_s(i))*ue(i,j,1)
                   ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1._r8-edge_vect_s(i))*ue(i,j,2)
                   ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1._r8-edge_vect_s(i))*ue(i,j,3)
                endif
             enddo
             do i=is,ie
                ue(i,j,1) = ut1(i)
                ue(i,j,2) = ut2(i)
                ue(i,j,3) = ut3(i)
             enddo
          endif
          if ( (je+1)==npy ) then
             j = npy
             do i=is,ie
                if ( i>im2 ) then
                   ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1._r8-edge_vect_n(i))*ue(i,j,1)
                   ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1._r8-edge_vect_n(i))*ue(i,j,2)
                   ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1._r8-edge_vect_n(i))*ue(i,j,3)
                else
                   ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1._r8-edge_vect_n(i))*ue(i,j,1)
                   ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1._r8-edge_vect_n(i))*ue(i,j,2)
                   ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1._r8-edge_vect_n(i))*ue(i,j,3)
                endif
             enddo
             do i=is,ie
                ue(i,j,1) = ut1(i)
                ue(i,j,2) = ut2(i)
                ue(i,j,3) = ut3(i)
             enddo
          endif

       endif ! .not. nested

       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) + dt5*( ue(i,j,1)*es(1,i,j,1) +  &
                  ue(i,j,2)*es(2,i,j,1) +  &
                  ue(i,j,3)*es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) + dt5*( ve(i,j,1)*ew(1,i,j,2) +  &
                  ve(i,j,2)*ew(2,i,j,2) +  &
                  ve(i,j,3)*ew(3,i,j,2) )
          enddo
       enddo
    enddo         ! k-loop

    call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)

end subroutine atend2dstate3d


subroutine fv3_tracer_diags(atm)

  ! Dry/Wet surface pressure diagnostics

  use constituents,          only: pcnst
  use dimensions_mod,        only: nlev,cnst_name_ffsl
  use dyn_grid,              only: mytile
  use fv_arrays_mod,         only: fv_atmos_type
  use physconst,             only: thermodynamic_active_species_num,thermodynamic_active_species_idx_dycore, &
                                   dry_air_species_num

  ! arguments
  type (fv_atmos_type), intent(in),  pointer :: Atm(:)

  ! Locals
  integer                      :: i, j ,k, m,is,ie,js,je
  integer                      :: num_wet_species       ! total number of wet species
  integer                      :: kstrat,ng
  real(r8)                     :: global_ps,global_dryps
  real(r8)                     :: qm_strat
  real(r8)                     :: qtot(pcnst), psum
  real(r8), allocatable, dimension(:,:,:)  :: delpdry, psq
  real(r8), allocatable, dimension(:,:)    :: psdry, q_strat

  !----------------------------------------------------------------------------

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  ng = Atm(mytile)%ng

  allocate(delpdry(is:ie,js:je,nlev))
  allocate(psdry(is:ie,js:je))
  allocate(psq(is:ie,js:je,pcnst))
  allocate(q_strat(is:ie,js:je))
  num_wet_species=thermodynamic_active_species_num-dry_air_species_num
  do k=1,nlev
     do j = js, je
        do i = is, ie
           delpdry(i,j,k) = Atm(mytile)%delp(i,j,k) * &
                            (1.0_r8-sum(Atm(mytile)%q(i,j,k,thermodynamic_active_species_idx_dycore(1:num_wet_species))))
        end do
     end do
  end do
  !
  ! get psdry
  !
  do j = js, je
     do i = is, ie
        psdry(i,j) = hyai(1)*ps0 + sum(delpdry(i,j,:))
     end do
  end do

  global_ps = g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, &
              Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
  global_dryps = g_sum(Atm(mytile)%domain, psdry(is:ie,js:je), is, ie, js, je, &
                 Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
!-------------------
! Vertical mass sum for all tracers
!-------------------
  psq(:,:,:) = 0._r8
  do m=1,pcnst
     call z_sum(Atm,is,ie,js,je,nlev,Atm(mytile)%q(is:ie,js:je,1:nlev,m),psq(is:ie,js:je,m))
  end do
! Mean water vapor in the "stratosphere" (75 mb and above):
  qm_strat = 0._r8
  if ( Atm(mytile)%idiag%phalf(2)< 75._r8 ) then
     kstrat = 1
     do k=2,nlev
        if ( Atm(mytile)%idiag%phalf(k+1) > 75._r8 ) exit
        kstrat = k
     enddo
     call z_sum(Atm,is,ie,js,je, kstrat, Atm(mytile)%q(is:ie,js:je,1:kstrat,1 ), q_strat,psum)
     qm_strat = g_sum(Atm(mytile)%domain, q_strat(is:ie,js:je), is, ie, js, je, &
                Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1) * 1.e6_r8 / psum
  endif

  !-------------------
  ! Get global mean mass for all tracers
  !-------------------
  do m=1,pcnst
     qtot(m) = g_sum(Atm(mytile)%domain, psq(is,js,m), is, ie, js, je, &
               Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)/gravit
  enddo

  if (masterproc) then
     write(iulog,*)'Total Surface Pressure (mb)                  = ',global_ps/100.0_r8,"hPa"
     write(iulog,*)'Mean Dry Surface Pressure (mb)               = ',global_dryps/100.0_r8,"hPa"
     write(iulog,*)'Mean specific humidity (mg/kg) above 75 mb   = ',qm_strat
     do m=1,pcnst
        write(iulog,*)' Total '//cnst_name_ffsl(m)//' (kg/m**2)     = ',qtot(m)
     enddo
  end if


  deallocate(delpdry)
  deallocate(psdry)
  deallocate(psq)
  deallocate(q_strat)
end subroutine fv3_tracer_diags


subroutine z_sum(atm,is,ie,js,je,km,q,msum,gpsum)

  ! vertical integral

  use fv_arrays_mod,   only: fv_atmos_type

  ! arguments

  type (fv_atmos_type), intent(in),  pointer        :: Atm(:)
  integer, intent(in)                               :: is, ie, js, je
  integer, intent(in)                               :: km
  real(r8), intent(in), dimension(is:ie, js:je, km) :: q
  real(r8), intent(out), dimension(is:ie,js:je)     :: msum
  real(r8), intent(out), optional                   :: gpsum

  ! LOCAL VARIABLES
  integer :: i,j,k
  real(r8), dimension(is:ie,js:je)           :: psum
  !----------------------------------------------------------------------------
  msum=0._r8
  psum=0._r8
  do j=js,je
     do i=is,ie
        msum(i,j) = Atm(mytile)%delp(i,j,1)*q(i,j,1)
        psum(i,j) = Atm(mytile)%delp(i,j,1)
     enddo
     do k=2,km
        do i=is,ie
           msum(i,j) = msum(i,j) + Atm(mytile)%delp(i,j,k)*q(i,j,k)
           psum(i,j) = psum(i,j) + Atm(mytile)%delp(i,j,k)
        enddo
     enddo
  enddo
  if (present(gpsum)) then
     gpsum = g_sum(Atm(mytile)%domain, psum, is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
  end if
end subroutine z_sum

end module dp_coupling
