module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use pmgrid,         only: plev
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_type
use physconst,      only: gravit, cappa, rh2o, zvir
use air_composition,only: cpairv, rairv

use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
use spmd_utils,     only: mpicom, iam, masterproc

use dyn_comp,       only: dyn_export_t, dyn_import_t

use physics_types,  only: physics_state, physics_tend, physics_cnst_limit
use phys_grid,      only: get_dyn_col_p, get_chunk_info_p, get_ncols_p, get_gcol_all_p
use phys_grid,      only: columns_on_task

use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field

use cam_logfile,    only: iulog
use perf_mod,       only: t_startf, t_stopf, t_barrierf
use cam_abortutils, only: endrun
use air_composition,only: thermodynamic_active_species_num,thermodynamic_active_species_idx,thermodynamic_active_species_idx_dycore
implicit none
private
save
logical :: compute_energy_diags=.false.
integer :: index_qv_phys = -1

public :: &
   d_p_coupling, &
   p_d_coupling

!=========================================================================================
contains
!=========================================================================================

subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)

   ! Convert the dynamics output state into the physics input state.
   ! Note that all pressures and tracer mixing ratios coming from the dycore are based on
   ! dry air mass.
   use cam_history,    only : hist_fld_active
   use mpas_constants, only : Rv_over_Rd => rvord

   ! arguments
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_export_t),        intent(inout) :: dyn_out
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! local variables

   ! Variables from dynamics export container
   integer :: nCellsSolve
   integer :: index_qv
   integer, dimension(:), pointer :: cam_from_mpas_cnst

   real(r8), pointer :: exner(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8), pointer :: w(:,:)
   real(r8), pointer :: theta_m(:,:)
   real(r8), pointer :: tracers(:,:,:)


   integer :: lchnk, icol, icol_p, k, kk      ! indices over chunks, columns, physics columns and layers
   integer :: i, m, ncols, blockid
   integer :: block_index
   integer, dimension(:), pointer :: block_offset

   integer :: pgcols(pcols)
   integer :: tsize                    ! amount of data per grid point passed to physics
   integer, allocatable :: bpter(:,:)  ! offsets into block buffer for packing data
   integer, allocatable :: cpter(:,:)  ! offsets into chunk buffer for unpacking data

   real(r8), allocatable:: pmid(:,:)      !mid-level hydrostatic pressure consistent with MPAS discrete state
   real(r8), allocatable:: pintdry(:,:)   !interface hydrostatic pressure consistent with MPAS discrete state
   real(r8), allocatable:: pmiddry(:,:)   !mid-level hydrostatic dry pressure consistent with MPAS discrete state

   integer :: ierr
   character(len=*), parameter :: subname = 'd_p_coupling'
   !----------------------------------------------------------------------------

   compute_energy_diags=&
        (hist_fld_active('SE_dBF').or.hist_fld_active('SE_dAP').or.hist_fld_active('SE_dAM').or.&
         hist_fld_active('KE_dBF').or.hist_fld_active('KE_dAP').or.hist_fld_active('KE_dAM').or.&
         hist_fld_active('WV_dBF').or.hist_fld_active('WV_dAP').or.hist_fld_active('WV_dAM'))

   nCellsSolve = dyn_out % nCellsSolve
   index_qv    = dyn_out % index_qv
   cam_from_mpas_cnst => dyn_out % cam_from_mpas_cnst

   zint     => dyn_out % zint
   zz       => dyn_out % zz
   rho_zz   => dyn_out % rho_zz
   ux       => dyn_out % ux
   uy       => dyn_out % uy
   w        => dyn_out % w
   theta_m  => dyn_out % theta_m
   exner    => dyn_out % exner
   tracers  => dyn_out % tracers

   if (compute_energy_diags) then
     call tot_energy(nCellsSolve, plev,size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), &
          rho_zz(:,1:nCellsSolve), theta_m(:,1:nCellsSolve), tracers(:,:,1:nCellsSolve),&
          ux(:,1:nCellsSolve),uy(:,1:nCellsSolve),'dBF')
   end if
   !
   ! diagnose pintdry, pmiddry, pmid
   !
   allocate(pmid(plev, nCellsSolve),      stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate pmid array')

   allocate(pmiddry(plev, nCellsSolve),   stat=ierr)!note: .neq. dyn_out % pmiddry since it is non-hydrostatic
   if( ierr /= 0 ) call endrun(subname//':failed to allocate pmiddry array')

   allocate(pintdry(plev+1, nCellsSolve), stat=ierr)!note: .neq. dyn_out % pintdry since it is non-hydrostatic
   if( ierr /= 0 ) call endrun(subname//':failed to allocate pintdry array')

   call hydrostatic_pressure( &
        nCellsSolve, plev, zz, zint, rho_zz, theta_m, tracers(index_qv,:,:),&
        pmiddry, pintdry, pmid)

   call t_startf('dpcopy')

   ncols = columns_on_task

   allocate(block_offset(0), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate block_offset array')

   do icol = 1, ncols
      call get_dyn_col_p(icol, block_index, block_offset) ! Get the dynamic column (block_index)
      call get_chunk_info_p(icol, lchnk, icol_p)          ! Get the matching physics column (icol_p)
      i = block_index

      phys_state(lchnk)%psdry(icol_p) = pintdry(1,i)
      phys_state(lchnk)%phis(icol_p) = zint(1,i) * gravit

      do k = 1, pver                                  ! vertical index in physics chunk
         kk = pver - k + 1                            ! vertical index in dynamics block

         phys_state(lchnk)%t(icol_p,k)       = theta_m(kk,i) / (1.0_r8 + &
                                             Rv_over_Rd * tracers(index_qv,kk,i)) * exner(kk,i)
         phys_state(lchnk)%u(icol_p,k)       = ux(kk,i)
         phys_state(lchnk)%v(icol_p,k)       = uy(kk,i)
         phys_state(lchnk)%omega(icol_p,k)   = -rho_zz(kk,i)*zz(kk,i)*gravit*0.5_r8*(w(kk,i)+w(kk+1,i))   ! omega
         phys_state(lchnk)%pmiddry(icol_p,k) = pmiddry(kk,i)
         phys_state(lchnk)%pmid(icol_p,k)    = pmid(kk,i)
      end do

      do k = 1, pverp
         kk = pverp - k + 1
         phys_state(lchnk)%pintdry(icol_p,k) = pintdry(kk,i)
      end do

      do m = 1, pcnst
         do k = 1, pver
            kk = pver - k + 1
            phys_state(lchnk)%q(icol_p,k,m) = tracers(cam_from_mpas_cnst(m),kk,i)
         end do
      end do
   end do

   call t_stopf('dpcopy')

   call t_startf('derived_phys')
   call derived_phys(phys_state, phys_tend, pbuf2d)
   call t_stopf('derived_phys')

   deallocate(pmid,pintdry,pmiddry)

end subroutine d_p_coupling

!=========================================================================================

subroutine p_d_coupling(phys_state, phys_tend, dyn_in)

   ! Convert the physics output state and tendencies into the dynamics
   ! input state.  Begin by redistributing the output of the physics package
   ! to the block data structure.  Then derive the tendencies required by
   ! MPAS.

   use cam_mpas_subdriver, only : domain_ptr
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_field, mpas_pool_get_array
   use mpas_field_routines, only : mpas_allocate_scratch_field, mpas_deallocate_scratch_field
   use mpas_derived_types, only : mpas_pool_type, field2DReal
   use time_manager, only : get_step_size

   ! Arguments
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ), intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_import_t),  intent(inout) :: dyn_in

   ! Local variables
   integer :: lchnk, icol, icol_p, k, kk      ! indices over chunks, columns, layers
   integer :: i, m, ncols, blockid
   integer :: block_index
   integer, dimension(:), pointer :: block_offset

   real(r8) :: factor
   real(r8) :: dt_phys

   ! Variables from dynamics import container
   integer :: nCellsSolve
   integer :: nCells
   integer :: nEdgesSolve
   integer :: index_qv
   integer, dimension(:), pointer :: mpas_from_cam_cnst

   real(r8), pointer :: tracers(:,:,:)

   ! CAM physics output redistributed to blocks.
   real(r8), allocatable :: t_tend(:,:)
   real(r8), allocatable :: q_tend(:,:,:)
   real(r8), pointer :: u_tend(:,:)
   real(r8), pointer :: v_tend(:,:)

   integer :: idx_phys, idx_dycore

   integer :: pgcols(pcols)
   integer :: tsize                    ! amount of data per grid point passed to dynamics
   integer, allocatable :: bpter(:,:)  ! offsets into block buffer for unpacking data
   integer, allocatable :: cpter(:,:)  ! offsets into chunk buffer for packing data

   type (mpas_pool_type), pointer :: tend_physics
   type (field2DReal), pointer :: tend_uzonal, tend_umerid

   integer :: ierr
   character(len=*), parameter :: subname = 'dp_coupling::p_d_coupling'
   !----------------------------------------------------------------------------

   nCellsSolve = dyn_in % nCellsSolve
   nCells      = dyn_in % nCells
   index_qv    = dyn_in % index_qv
   mpas_from_cam_cnst => dyn_in % mpas_from_cam_cnst

   tracers => dyn_in % tracers

   allocate( t_tend(pver,nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate t_tend array')
   allocate( q_tend(thermodynamic_active_species_num,pver,nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate q_tend array')

   nullify(tend_physics)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics)

   nullify(tend_uzonal)
   nullify(tend_umerid)
   call mpas_pool_get_field(tend_physics, 'tend_uzonal', tend_uzonal)
   call mpas_pool_get_field(tend_physics, 'tend_umerid', tend_umerid)
   call mpas_allocate_scratch_field(tend_uzonal)
   call mpas_allocate_scratch_field(tend_umerid)
   call mpas_pool_get_array(tend_physics, 'tend_uzonal', u_tend)
   call mpas_pool_get_array(tend_physics, 'tend_umerid', v_tend)

   ! Physics coupling interval, used to compute tendency of qv
   dt_phys = get_step_size()

   call t_startf('pd_copy')

   ncols = columns_on_task ! should this be nCellsSolve?

   allocate(block_offset(0), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate block_offset array')
   do icol = 1, ncols                                   ! column index in physics chunk
      ! Get dynamics block
      call get_chunk_info_p(icol, lchnk, icol_p)          ! Get the matching physics column (icol_p)
      call get_dyn_col_p(icol, block_index, block_offset) ! Get the dynamic column (block_index)
      i = block_index

      do k = 1, pver                              ! vertical index in physics chunk
         kk = pver - k + 1                        ! vertical index in dynamics block

         t_tend(kk,i) = phys_tend(lchnk)%dtdt(icol_p,k)
         u_tend(kk,i) = phys_tend(lchnk)%dudt(icol_p,k)
         v_tend(kk,i) = phys_tend(lchnk)%dvdt(icol_p,k)

         ! convert wet mixing ratios to dry
         factor = phys_state(lchnk)%pdel(icol_p,k)/phys_state(lchnk)%pdeldry(icol_p,k)
         !
         ! compute tendencies for thermodynamic active species
         !
         do m=1,thermodynamic_active_species_num
           idx_phys   = thermodynamic_active_species_idx(m)
           idx_dycore = thermodynamic_active_species_idx_dycore(m)
           if (idx_dycore==index_qv) index_qv_phys = m
           if (cnst_type(idx_phys) == 'wet') then
             q_tend(m,kk,i) = (phys_state(lchnk)%q(icol_p,k,idx_phys)*factor - tracers(idx_dycore,kk,i)) / dt_phys
           else
             q_tend(m,kk,i) = (phys_state(lchnk)%q(icol_p,k,idx_phys) - tracers(idx_dycore,kk,i)) / dt_phys
           end if
         end do

         do m = 1, pcnst
           if (cnst_type(mpas_from_cam_cnst(m)) == 'wet') then
             tracers(m,kk,i) = phys_state(lchnk)%q(icol_p,k,mpas_from_cam_cnst(m))*factor
           else
             tracers(m,kk,i) = phys_state(lchnk)%q(icol_p,k,mpas_from_cam_cnst(m))
           end if
         end do
      end do
   end do

   call t_stopf('pd_copy')

   call t_startf('derived_tend')
   call derived_tend(nCellsSolve, nCells, t_tend, u_tend, v_tend, q_tend, dyn_in)
   call t_stopf('derived_tend')

   call mpas_deallocate_scratch_field(tend_uzonal)
   call mpas_deallocate_scratch_field(tend_umerid)

end subroutine p_d_coupling

!=========================================================================================

subroutine derived_phys(phys_state, phys_tend, pbuf2d)

   ! Compute fields in the physics state object which are diagnosed from the
   ! MPAS prognostic fields.

   use geopotential,    only: geopotential_t
   use check_energy,    only: check_energy_timestep_init
   use shr_vmath_mod,   only: shr_vmath_log
   use phys_control,    only: waccmx_is
   use cam_thermo,      only: cam_thermo_update
   use air_composition, only: rairv
   use qneg_module,     only: qneg3
   use shr_const_mod,   only: shr_const_rwv
   use constituents,    only: qmin
   ! Arguments
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! Local variables

   integer :: i, k, lchnk, m, ncol

   real(r8) :: factor(pcols,pver)
   real(r8) :: zvirv(pcols,pver)

   real(r8), parameter :: pref = 1.e5_r8 ! reference pressure (Pa)

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   character(len=*), parameter :: subname = 'dp_coupling::derived_phys'

   !$omp parallel do private (lchnk, ncol, k, factor)
   do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)

      ! The dry pressure profiles are derived using hydrostatic formulas
      ! and the dry air mass from MPAS.

      ! Derived variables for dry pressure profiles:

      do k = 1, pver

         phys_state(lchnk)%pdeldry(:ncol,k) = phys_state(lchnk)%pintdry(:ncol,k+1) - &
                                              phys_state(lchnk)%pintdry(:ncol,k)

         phys_state(lchnk)%rpdeldry(:ncol,k) = 1._r8 / phys_state(lchnk)%pdeldry(:ncol,k)

         call shr_vmath_log(phys_state(lchnk)%pintdry(:ncol,k), &
                            phys_state(lchnk)%lnpintdry(:ncol,k), ncol)

         call shr_vmath_log(phys_state(lchnk)%pmiddry(:ncol,k), &
                            phys_state(lchnk)%lnpmiddry(:ncol,k), ncol)

      end do

      call shr_vmath_log(phys_state(lchnk)%pintdry(:ncol,pverp), &
                         phys_state(lchnk)%lnpintdry(:ncol,pverp), ncol)


      ! Add in the water vapor mass to compute the moist pressure profiles
      ! used by CAM's physics packages.
      ! **N.B.** The input water vapor mixing ratio in phys_state is based on dry air.  It
      !          gets converted to a wet basis later.

      do k = 1, pver
         ! To be consistent with total energy formula in physic's check_energy module only
         ! include water vapor in moist pdel.
         factor(:ncol,k) = 1._r8 + phys_state(lchnk)%q(:ncol,k,1)
         phys_state(lchnk)%pdel(:ncol,k)  = phys_state(lchnk)%pdeldry(:ncol,k)*factor(:ncol,k)
         phys_state(lchnk)%rpdel(:ncol,k) = 1._r8 / phys_state(lchnk)%pdel(:ncol,k)
      end do

      ! Assume no water vapor above top of model.
      phys_state(lchnk)%pint(:ncol,1) = phys_state(lchnk)%pintdry(:ncol,1)
      do k = 1, pver
         phys_state(lchnk)%pint(:ncol,k+1) = phys_state(lchnk)%pint(:ncol,k) + &
                                             phys_state(lchnk)%pdel(:ncol,k)
         call shr_vmath_log(phys_state(lchnk)%pint(:ncol,k), &
                            phys_state(lchnk)%lnpint(:ncol,k), ncol)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(:ncol,pverp), &
                         phys_state(lchnk)%lnpint(:ncol,pverp), ncol)

      phys_state(lchnk)%ps(:ncol) = phys_state(lchnk)%pint(:ncol,pverp)

      do k = 1, pver
         call shr_vmath_log(phys_state(lchnk)%pmid(:ncol,k), &
                            phys_state(lchnk)%lnpmid(:ncol,k), ncol)
      end do

      do k = 1, pver
         phys_state(lchnk)%exner(:ncol,k) = (pref / phys_state(lchnk)%pmid(:ncol,k))**cappa
      end do

      ! Tracers from MPAS are in dry mixing ratio units.  CAM's physics package expects constituents
      ! which have been declared to be type 'wet' when they are registered to be represented by mixing
      ! ratios based on moist air mass (dry air + water vapor).  Do appropriate conversion here.
      factor(:ncol,:) = 1._r8/factor(:ncol,:)
      do m = 1,pcnst
         if (cnst_type(m) == 'wet') then
            phys_state(lchnk)%q(:ncol,:,m) = factor(:ncol,:)*phys_state(lchnk)%q(:ncol,:,m)
         end if
      end do


      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
        !------------------------------------------------------------
        ! Apply limiters to mixing ratios of major species
        !------------------------------------------------------------
        call physics_cnst_limit( phys_state(lchnk) )
        !-----------------------------------------------------------------------------
        ! Call cam_thermo_update to compute cpairv, rairv, mbarv, and cappav as
        ! constituent dependent variables.
        ! Compute molecular viscosity(kmvis) and conductivity(kmcnd).
        ! Fill local zvirv variable; calculated for WACCM-X.
        !-----------------------------------------------------------------------------
        call cam_thermo_update(phys_state(lchnk)%q, phys_state(lchnk)%t, lchnk, ncol)
        zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) -1._r8
      else
        zvirv(:,:) = zvir
      endif

      ! Compute geopotential height above surface - based on full pressure
      ! Note that phys_state%zi(:,plev+1) = 0 whereas zint in MPAS is surface height
      !
      call geopotential_t( &
         phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid,   phys_state(lchnk)%pint,          &
         phys_state(lchnk)%pmid,   phys_state(lchnk)%pdel,     phys_state(lchnk)%rpdel,         &
         phys_state(lchnk)%t,      phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk), gravit, zvirv, &
         phys_state(lchnk)%zi,     phys_state(lchnk)%zm,       ncol)

      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
         phys_state(lchnk)%s(:ncol,k) = cpairv(:ncol,k,lchnk)*phys_state(lchnk)%t(:ncol,k) &
            + gravit*phys_state(lchnk)%zm(:ncol,k) + phys_state(lchnk)%phis(:ncol)
      end do

      ! Ensure tracers are all positive
      call qneg3('D_P_COUPLING',lchnk  ,ncol    ,pcols   ,pver    , &
           1, pcnst, qmin  ,phys_state(lchnk)%q)


      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

   end do

end subroutine derived_phys

!=========================================================================================

subroutine derived_tend(nCellsSolve, nCells, t_tend, u_tend, v_tend, q_tend, dyn_in)

   ! Derive the physics tendencies required by MPAS from the tendencies produced by
   ! CAM's physics package.
   use mpas_constants, only: p0,cv,rgas,cp
   use cam_mpas_subdriver, only : cam_mpas_cell_to_edge_winds, cam_mpas_update_halo
   use mpas_constants,     only : Rv_over_Rd => rvord
   use time_manager,       only : get_step_size

   ! Arguments
   integer,             intent(in)    :: nCellsSolve
   integer,             intent(in)    :: nCells
   real(r8),            intent(in)    :: t_tend(pver,nCellsSolve)  ! physics dtdt
   real(r8),            intent(in)    :: q_tend(thermodynamic_active_species_num,pver,nCellsSolve) ! physics dqvdt
   real(r8),            intent(inout) :: u_tend(pver,nCells+1)     ! physics dudt
   real(r8),            intent(inout) :: v_tend(pver,nCells+1)     ! physics dvdt
   type(dyn_import_t),  intent(inout) :: dyn_in


   ! Local variables
   real(r8) :: dtime
   ! variables from dynamics import container
   integer :: nEdges
   real(r8), pointer :: ru_tend(:,:)
   real(r8), pointer :: rtheta_tend(:,:)
   real(r8), pointer :: rho_tend(:,:)

   real(r8), pointer :: normal(:,:)
   real(r8), pointer :: east(:,:)
   real(r8), pointer :: north(:,:)
   integer, pointer :: cellsOnEdge(:,:)

   real(r8), pointer :: theta(:,:)
   real(r8), pointer :: exner(:,:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: tracers(:,:,:)

   integer :: index_qv,m,idx_dycore
   real(r8) :: rhok,thetavk,thetak,pk,exnerk,tempk,tempk_new,exnerk_new,thetak_new,thetak_m_new,rhodk,tknew,thetaknew
   !
   ! variables for energy diagnostics
   !
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: theta_m(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8)          :: theta_m_new(pver,nCellsSolve) !modified potential temperature after various physics updates
   real(r8)          :: rtheta_param(pver,nCellsSolve)!tendency from temperature change only (for diagnostics)
   real(r8)          :: qk (thermodynamic_active_species_num,pver,nCellsSolve) !water species before physics (diagnostics)
   real(r8)          :: qwv(pver,nCellsSolve)                                  !water vapor before physics
   real(r8)          :: facnew, facold
   real(r8), allocatable :: tracers_old(:,:,:)

   integer  :: iCell,k

   character(len=*), parameter :: subname = 'dp_coupling:derived_tend'
   !----------------------------------------------------------------------------

   nEdges = dyn_in % nEdges
   ru_tend     => dyn_in % ru_tend
   rtheta_tend => dyn_in % rtheta_tend
   rho_tend    => dyn_in % rho_tend

   east        => dyn_in % east
   north       => dyn_in % north
   normal      => dyn_in % normal
   cellsOnEdge => dyn_in % cellsOnEdge

   theta       => dyn_in % theta
   exner       => dyn_in % exner
   rho_zz      => dyn_in % rho_zz
   tracers     => dyn_in % tracers
   index_qv    =  dyn_in % index_qv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Momentum tendency
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Couple u and v tendencies with rho_zz
   !
   u_tend(:,:) = u_tend(:,:) * rho_zz(:,:)
   v_tend(:,:) = v_tend(:,:) * rho_zz(:,:)

   !
   ! Update halos for u_tend and v_tend
   !
   call cam_mpas_update_halo('tend_uzonal', endrun)   ! dyn_in % u_tend
   call cam_mpas_update_halo('tend_umerid', endrun)   ! dyn_in % v_tend

   !
   ! Project u and v tendencies to edge normal tendency
   !
   call cam_mpas_cell_to_edge_winds(nEdges, u_tend, v_tend, east, north, normal, &
                                    cellsOnEdge, ru_tend)

   !
   ! Update halo for edge normal tendency
   !
   call cam_mpas_update_halo('tend_ru_physics', endrun)

   dtime = get_step_size()

   zz       => dyn_in % zz
   theta_m  => dyn_in % theta_m
   zint     => dyn_in % zint
   ux       => dyn_in % ux
   uy       => dyn_in % uy
   !
   ! Compute q not updated by physics
   !
   qwv = tracers(index_qv,:,1:nCellsSolve)-dtime*q_tend(index_qv_phys,:,1:nCellsSolve)

   do iCell = 1, nCellsSolve
     do k = 1, pver
       rhodk     = zz(k,iCell) * rho_zz(k,iCell)
       facold    = 1.0_r8 + Rv_over_Rd *qwv(k,iCell)
       thetak    = theta_m(k,iCell)/facold

       exnerk    = (rgas*rhodk*theta_m(k,iCell)/p0)**(rgas/cv)
       tknew     = exnerk*thetak+(cp/cv)*dtime*t_tend(k,icell)


       thetaknew = (tknew**(cv/cp))*((rgas*rhodk*facold)/p0)**(-rgas/cp)
       !
       ! calculate theta_m tendency due to parameterizations (but no water adjustment)
       !
       rtheta_param(k,iCell) = (thetaknew-thetak)/dtime
       rtheta_param(k,iCell) = rtheta_param(k,iCell)*(1.0_r8 + Rv_over_Rd *qwv(k,iCell)) !convert to thetam
       rtheta_param(k,iCell) = rtheta_param(k,iCell)*rho_zz(k,iCell)
       !
       ! include water change in theta_m
       !
       facnew               = 1.0_r8 + Rv_over_Rd *tracers(index_qv,k,iCell)
       thetaknew            = (tknew**(cv/cp))*((rgas*rhodk*facnew)/p0)**(-rgas/cp)
       rtheta_tend(k,iCell) = (thetaknew*facnew-thetak*facold)/dtime
       rtheta_tend(k,iCell) = rtheta_tend(k,iCell) * rho_zz(k,iCell)
     end do
   end do

   if (compute_energy_diags) then
     !
     ! compute energy based on parameterization increment (excl. water change)
     !
     theta_m_new = theta_m(:,1:nCellsSolve)+dtime*rtheta_param(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve)
     !
     ! temporarily save thermodynamic active species (n+1)
     !
     do m=1,thermodynamic_active_species_num
       idx_dycore                         = thermodynamic_active_species_idx_dycore(m)
       qk(m,:,: )                         = tracers(idx_dycore,:,1:nCellsSolve)
       tracers(idx_dycore,:,1:nCellsSolve)= qk(m,:,: )-dtime*q_tend(m,:,1:nCellsSolve)
     end do

     call tot_energy( &
          nCellsSolve, plev, size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), rho_zz(:,1:nCellsSolve), &
          theta_m_new,  tracers(:,:,1:nCellsSolve),   &
          ux(:,1:nCellsSolve)+dtime*u_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),       &
          uy(:,1:nCellsSolve)+dtime*v_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),'dAP')
     ! revert
     do m=1,thermodynamic_active_species_num
       idx_dycore                         = thermodynamic_active_species_idx_dycore(m)
       tracers(idx_dycore,:,1:nCellsSolve)= qk(m,:,: )
     end do
     !
     ! compute energy incl. water change
     !
     theta_m_new = theta_m(:,1:nCellsSolve)+dtime*rtheta_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve)
     call tot_energy( &
          nCellsSolve, plev, size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), &
          rho_zz(:,1:nCellsSolve), theta_m_new, tracers(:,:,1:nCellsSolve),    &
          ux(:,1:nCellsSolve)+dtime*u_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),       &
          uy(:,1:nCellsSolve)+dtime*v_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),'dAM')
   end if
   !
   ! Update halo for rtheta_m tendency
   !
   call cam_mpas_update_halo('tend_rtheta_physics', endrun)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Density tendency
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   rho_tend = 0.0_r8
end subroutine derived_tend

!=========================================================================================
subroutine hydrostatic_pressure(nCells, nVertLevels, zz, zgrid, rho_zz, theta_m, q, pmiddry, pintdry,pmid)

   ! Compute dry hydrostatic pressure at layer interfaces and midpoints
   !
   ! Given arrays of zz, zgrid, rho_zz, and theta_m from the MPAS-A prognostic
   ! state, compute dry hydrostatic pressure at layer interfaces and midpoints.
   ! The vertical dimension for 3-d arrays is innermost, and k=1 represents
   ! the lowest layer or level in the fields.
   !
  use mpas_constants, only : cp, rgas, cv, gravity, p0, Rv_over_Rd => rvord
  use physconst,      only:  rair, cpair

   ! Arguments
   integer, intent(in) :: nCells
   integer, intent(in) :: nVertLevels
   real(r8), dimension(nVertLevels, nCells),   intent(in) :: zz      ! d(zeta)/dz [-]
   real(r8), dimension(nVertLevels+1, nCells), intent(in) :: zgrid   ! geometric heights of layer interfaces [m]
   real(r8), dimension(nVertLevels, nCells),   intent(in) :: rho_zz  ! dry density / zz [kg m^-3]
   real(r8), dimension(nVertLevels, nCells),   intent(in) :: theta_m ! modified potential temperature
   real(r8), dimension(nVertLevels, nCells),   intent(in) :: q       ! water vapor dry mixing ratio
   real(r8), dimension(nVertLevels, nCells),   intent(out):: pmiddry ! layer midpoint dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels+1, nCells), intent(out):: pintdry ! layer interface dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels, nCells),   intent(out):: pmid    ! layer midpoint hydrostatic pressure [Pa]

   ! Local variables
   integer :: iCell, k
   real(r8), dimension(nVertLevels)   :: dz    ! Geometric layer thickness in column
   real(r8), dimension(nVertLevels+1) :: pint  ! hydrostatic pressure at interface
   real(r8) :: pi, t
   real(r8) :: pk,rhok,rhodryk,theta,thetavk,kap1,kap2

   !
   ! For each column, integrate downward from model top to compute dry hydrostatic pressure at layer
   ! midpoints and interfaces. The pressure averaged to layer midpoints should be consistent with
   ! the ideal gas law using the rho_zz and theta values prognosed by MPAS at layer midpoints.
   !
   kap1 = p0**(-rgas/cp)           ! pre-compute constants
   kap2 = cp/cv                    ! pre-compute constants
   do iCell = 1, nCells

      dz(:) = zgrid(2:nVertLevels+1,iCell) - zgrid(1:nVertLevels,iCell)

      k = nVertLevels
      rhok    = (1.0_r8+q(k,iCell))*zz(k,iCell) * rho_zz(k,iCell) !full CAM physics density
      thetavk = theta_m(k,iCell)/ (1.0_r8 + q(k,iCell))           !convert modified theta to virtual theta
      pk      = (rhok*rgas*thetavk*kap1)**kap2                    !mid-level top pressure
      !
      ! model top pressure consistently diagnosed using the assumption that the mid level
      ! is at height z(nVertLevels-1)+0.5*dz
      !
      pintdry(nVertLevels+1,iCell) = pk-0.5_r8*dz(nVertLevels)*rhok*gravity  !hydrostatic
      pint   (nVertLevels+1)       = pintdry(nVertLevels+1,iCell)
      do k = nVertLevels, 1, -1
        !
        ! compute hydrostatic dry interface pressure so that (pintdry(k+1)-pintdry(k))/g is pseudo density
        !
        rhodryk = zz(k,iCell) * rho_zz(k,iCell)
        rhok    = (1.0_r8+q(k,iCell))*rhodryk
        pintdry(k,iCell) = pintdry(k+1,iCell) + gravity * rhodryk * dz(k)
        pint   (k)       = pint   (k+1)       + gravity * rhok    * dz(k)
      end do

      do k = nVertLevels, 1, -1
        !hydrostatic mid-level pressure - MPAS full pressure is (rhok*rgas*thetavk*kap1)**kap2 
        pmid   (k,iCell) = 0.5_r8*(pint(k+1)+pint(k))       
        !hydrostatic dry mid-level dry pressure - 
        !MPAS non-hydrostatic dry pressure is pmiddry(k,iCell) = (rhodryk*rgas*theta*kap1)**kap2
        pmiddry(k,iCell) = 0.5_r8*(pintdry(k+1,iCell)+pintdry(k,iCell))  
      end do
    end do
end subroutine hydrostatic_pressure


subroutine tot_energy(nCells, nVertLevels, qsize, index_qv, zz, zgrid, rho_zz, theta_m, q, ux,uy,outfld_name_suffix)
  use physconst,         only: rair, cpair, gravit,cappa!=R/cp (dry air)
  use mpas_constants,    only: p0,cv,rv,rgas,cp
  use cam_history,       only: outfld, hist_fld_active
  use mpas_constants,    only: Rv_over_Rd => rvord
  use air_composition,   only: thermodynamic_active_species_ice_idx_dycore,thermodynamic_active_species_liq_idx_dycore
  use air_composition,   only: thermodynamic_active_species_ice_num,thermodynamic_active_species_liq_num
  ! Arguments
  integer, intent(in) :: nCells
  integer, intent(in) :: nVertLevels
  integer, intent(in) :: qsize
  integer, intent(in) :: index_qv
  real(r8), dimension(nVertLevels, nCells),       intent(in) :: zz      ! d(zeta)/dz [-]
  real(r8), dimension(nVertLevels+1, nCells),     intent(in) :: zgrid   ! geometric heights of layer interfaces [m]
  real(r8), dimension(nVertLevels, nCells),       intent(in) :: rho_zz  ! dry density / zz [kg m^-3]
  real(r8), dimension(nVertLevels, nCells),       intent(in) :: theta_m ! modified potential temperature
  real(r8), dimension(qsize,nVertLevels, nCells), intent(in) :: q       ! tracer array
  real(r8), dimension(nVertLevels, nCells),       intent(in) :: ux      ! A-grid zonal velocity component
  real(r8), dimension(nVertLevels, nCells),       intent(in) :: uy      ! A-grid meridional velocity component
  character*(*),                                  intent(in) :: outfld_name_suffix ! suffix for "outfld" names

  ! Local variables
  integer :: iCell, k, idx
  real(r8) :: rho_dz,zcell,temperature,theta,pk,ptop,exner
  real(r8), dimension(nVertLevels, nCells) :: rhod, dz
  real(r8), dimension(nCells)              :: kinetic_energy,potential_energy,internal_energy,water_vapor,water_liq,water_ice

  real(r8), dimension(nCells) :: liq !total column integrated liquid
  real(r8), dimension(nCells) :: ice !total column integrated ice

  character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5

  name_out1 = 'SE_'   //trim(outfld_name_suffix)
  name_out2 = 'KE_'   //trim(outfld_name_suffix)
  name_out3 = 'WV_'   //trim(outfld_name_suffix)
  name_out4 = 'WL_'   //trim(outfld_name_suffix)
  name_out5 = 'WI_'   //trim(outfld_name_suffix)

  if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
       hist_fld_active(name_out4).or.hist_fld_active(name_out5)) then

    kinetic_energy   = 0.0_r8
    potential_energy = 0.0_r8
    internal_energy  = 0.0_r8
    water_vapor      = 0.0_r8

    do iCell = 1, nCells
      do k = 1, nVertLevels
        dz(k,iCell)   = zgrid(k+1,iCell) - zgrid(k,iCell)
        zcell         = 0.5_r8*(zgrid(k,iCell)+zgrid(k+1,iCell))
        rhod(k,iCell) = zz(k,iCell) * rho_zz(k,iCell)
        rho_dz        = (1.0_r8+q(index_qv,k,iCell))*rhod(k,iCell)*dz(k,iCell)
        theta         = theta_m(k,iCell)/(1.0_r8 + Rv_over_Rd *q(index_qv,k,iCell))!convert theta_m to theta

        exner         = (rgas*rhod(k,iCell)*theta_m(k,iCell)/p0)**(rgas/cv)
        temperature   = exner*theta

        water_vapor(iCell)      = water_vapor(iCell) + rhod(k,iCell)*q(index_qv,k,iCell)*dz(k,iCell)
        kinetic_energy(iCell)   = kinetic_energy(iCell)  + &
                                  0.5_r8*(ux(k,iCell)**2._r8+uy(k,iCell)**2._r8)*rho_dz
        potential_energy(iCell) = potential_energy(iCell)+ rho_dz*gravit*zcell
        internal_energy(iCell)  = internal_energy(iCell) + rho_dz*cv*temperature
      end do
      internal_energy(iCell)  = internal_energy(iCell) + potential_energy(iCell) !static energy
    end do
    call outfld(name_out1,internal_energy,ncells,1)
    call outfld(name_out2,kinetic_energy ,ncells,1)
    call outfld(name_out3,water_vapor    ,ncells,1)
    !
    ! vertical integral of total liquid water
    !
    if (hist_fld_active(name_out4)) then
      liq = 0._r8
      do idx = 1,thermodynamic_active_species_liq_num
        do iCell = 1, nCells
          do k = 1, nVertLevels
            liq(iCell) = liq(iCell) + &
                 q(thermodynamic_active_species_liq_idx_dycore(idx),k,iCell)*rhod(k,iCell)*dz(k,iCell)
          end do
        end do
      end do
      call outfld(name_out4,liq,ncells,1)
    end if
    !
    ! vertical integral of total frozen (ice) water
    !
    if (hist_fld_active(name_out5)) then
      ice = 0._r8
      do idx = 1,thermodynamic_active_species_ice_num
        do iCell = 1, nCells
          do k = 1, nVertLevels
            ice(iCell) = ice(iCell) + &
                 q(thermodynamic_active_species_ice_idx_dycore(idx),k,iCell)*rhod(k,iCell)*dz(k,iCell)
          end do
        end do
      end do
      call outfld(name_out5,ice,ncells,1)
    end if
  end if
 end subroutine tot_energy

end module dp_coupling
