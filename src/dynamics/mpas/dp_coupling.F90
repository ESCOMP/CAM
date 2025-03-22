module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use pmgrid,         only: plev
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_type
use physconst,      only: gravit, cappa, zvir
use air_composition,only: cpairv
use air_composition,only: dry_air_species_num
use dyn_comp,       only: dyn_export_t, dyn_import_t
use physics_types,  only: physics_state, physics_tend, physics_cnst_limit
use phys_grid,      only: get_dyn_col_p, get_chunk_info_p, get_ncols_p
use phys_grid,      only: columns_on_task

use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field

use perf_mod,       only: t_startf, t_stopf
use cam_abortutils, only: endrun
use air_composition,only: thermodynamic_active_species_num,thermodynamic_active_species_idx, &
                          thermodynamic_active_species_idx_dycore

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
   use cam_mpas_subdriver, only: cam_mpas_update_halo

   ! Convert the dynamics output state into the physics input state.
   ! Note that all pressures and tracer mixing ratios coming from the dycore are based on
   ! dry air mass.
   use cam_history,    only: hist_fld_active
   use dyn_comp,       only: frontgf_idx, frontga_idx
   use mpas_constants, only: Rv_over_Rd => rvord
   use phys_control,   only: use_gw_front, use_gw_front_igw
   use cam_budget,     only : thermo_budget_history

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

   !
   ! mesh information and coefficients needed for
   ! frontogenesis function calculation
   !
   real(r8), pointer :: defc_a(:,:)
   real(r8), pointer :: defc_b(:,:)
   real(r8), pointer :: cell_gradient_coef_x(:,:)
   real(r8), pointer :: cell_gradient_coef_y(:,:)
   real(r8), pointer :: edgesOnCell_sign(:,:)
   real(r8), pointer :: dvEdge(:)
   real(r8), pointer :: areaCell(:)

   integer, pointer :: cellsOnEdge(:,:)
   integer, pointer :: edgesOnCell(:,:)
   integer, pointer :: nEdgesOnCell(:)

   real(r8), pointer :: uperp(:,:)
   real(r8), pointer :: utangential(:,:)

   !
   ! local storage for frontogenesis function and angle
   !
   real(r8), pointer :: frontogenesisFunction(:,:)
   real(r8), pointer :: frontogenesisAngle(:,:)
   real(r8), pointer :: pbuf_frontgf(:,:)
   real(r8), pointer :: pbuf_frontga(:,:)
   real(r8), allocatable :: frontgf_phys(:,:,:)
   real(r8), allocatable :: frontga_phys(:,:,:)

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   integer :: lchnk, icol, icol_p, k, kk      ! indices over chunks, columns, physics columns and layers
   integer :: i, m, ncols
   integer :: block_index
   integer, dimension(:), pointer :: block_offset

   real(r8), allocatable:: pmid(:,:)      !mid-level hydrostatic pressure consistent with MPAS discrete state
   real(r8), allocatable:: pintdry(:,:)   !interface hydrostatic pressure consistent with MPAS discrete state
   real(r8), allocatable:: pmiddry(:,:)   !mid-level hydrostatic dry pressure consistent with MPAS discrete state

   integer :: ierr
   character(len=*), parameter :: subname = 'd_p_coupling'
   !----------------------------------------------------------------------------

   compute_energy_diags=thermo_budget_history

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
     call tot_energy_dyn(nCellsSolve, plev,size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), &
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
        nCellsSolve, plev, size(tracers, 1), index_qv, zz, zint, rho_zz, theta_m, exner, tracers,&
        pmiddry, pintdry, pmid)

   if (use_gw_front .or. use_gw_front_igw) then
      call cam_mpas_update_halo('scalars', endrun)   ! scalars is the name of tracers in the MPAS state pool
      nullify(pbuf_chnk)
      nullify(pbuf_frontgf)
      nullify(pbuf_frontga)
      !
      ! compute frontogenesis function and angle for gravity wave scheme
      !
      defc_a => dyn_out % defc_a
      defc_b => dyn_out % defc_b
      cell_gradient_coef_x => dyn_out % cell_gradient_coef_x
      cell_gradient_coef_y => dyn_out % cell_gradient_coef_y
      edgesOnCell_sign => dyn_out % edgesOnCell_sign
      dvEdge => dyn_out % dvEdge
      areaCell => dyn_out % areaCell
      cellsOnEdge => dyn_out % cellsOnEdge
      edgesOnCell => dyn_out % edgesOnCell
      nEdgesOnCell => dyn_out % nEdgesOnCell
      uperp => dyn_out % uperp
      utangential => dyn_out % utangential

      allocate(frontogenesisFunction(plev, nCellsSolve), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate frontogenesisFunction array')
      allocate(frontogenesisAngle(plev, nCellsSolve), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate frontogenesisAngle array')

      allocate(frontgf_phys(pcols, pver, begchunk:endchunk), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate frontgf_phys array')
      allocate(frontga_phys(pcols, pver, begchunk:endchunk), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//':failed to allocate frontga_phys array')


      call calc_frontogenesis( frontogenesisFunction, frontogenesisAngle,  &
                               theta_m, tracers(index_qv,:,:),             &
                               uperp, utangential, defc_a, defc_b,         &
                               cell_gradient_coef_x, cell_gradient_coef_y, &
                               areaCell, dvEdge, cellsOnEdge, edgesOnCell, &
                               nEdgesOnCell, edgesOnCell_sign,             &
                               plev, nCellsSolve )

   end if

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

         if (use_gw_front .or. use_gw_front_igw) then
            frontgf_phys(icol_p, k, lchnk) = frontogenesisFunction(kk, i)
            frontga_phys(icol_p, k, lchnk) = frontogenesisAngle(kk, i)
         end if
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

   if (use_gw_front .or. use_gw_front_igw) then

      !$omp parallel do private (lchnk, ncols, icol, k, pbuf_chnk, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)
         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
         call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
         call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
         do k = 1, pver
            do icol = 1, ncols
               pbuf_frontgf(icol, k) = frontgf_phys(icol, k, lchnk)
               pbuf_frontga(icol, k) = frontga_phys(icol, k, lchnk)
            end do
         end do
      end do
      deallocate(frontgf_phys)
      deallocate(frontga_phys)
      deallocate(frontogenesisFunction)
      deallocate(frontogenesisAngle)
   end if

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
   integer :: i, m, ncols
   integer :: block_index
   integer, dimension(:), pointer :: block_offset

   real(r8) :: factor
   real(r8) :: dt_phys

   ! Variables from dynamics import container
   integer :: nCellsSolve
   integer :: nCells
   integer :: index_qv
   integer, dimension(:), pointer :: mpas_from_cam_cnst

   real(r8), pointer :: tracers(:,:,:)

   ! CAM physics output redistributed to blocks.
   real(r8), allocatable :: t_tend(:,:)
   real(r8), allocatable :: q_tend(:,:,:)
   real(r8), pointer :: u_tend(:,:)
   real(r8), pointer :: v_tend(:,:)

   integer :: idx_phys, idx_dycore

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
   allocate( q_tend(thermodynamic_active_species_num-dry_air_species_num,pver,nCellsSolve), stat=ierr)
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
         do m=dry_air_species_num + 1,thermodynamic_active_species_num
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
   use cam_thermo,      only: cam_thermo_dry_air_update, cam_thermo_water_update
   use air_composition, only: rairv, dry_air_species_num
   use qneg_module,     only: qneg3
   use shr_const_mod,   only: shr_const_rwv
   use constituents,    only: qmin
   use dyn_tests_utils, only: vcoord=>vc_height
   ! Arguments
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! Local variables

   integer :: k, lchnk, m, ncol, m_cnst

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
         factor(:ncol,k) = 1.0_r8
         do m_cnst=dry_air_species_num + 1,thermodynamic_active_species_num
           m = thermodynamic_active_species_idx(m_cnst)
           ! at this point all q's are dry
           factor(:ncol,k) = factor(:ncol,k)+phys_state(lchnk)%q(:ncol,k,m)
         end do
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




      if (dry_air_species_num>0) then
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
        call cam_thermo_dry_air_update(phys_state(lchnk)%q, phys_state(lchnk)%t, lchnk, ncol)
        zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) -1._r8
      else
        zvirv(:,:) = zvir
      endif
      !
      ! update cp_dycore in module air_composition.
      ! (note: at this point q is dry)
      !
      call cam_thermo_water_update(phys_state(lchnk)%q(1:ncol,:,:), lchnk, ncol, vcoord)
      ! Tracers from MPAS are in dry mixing ratio units.  CAM's physics package expects constituents
      ! which have been declared to be type 'wet' when they are registered to be represented by mixing
      ! ratios based on moist air mass (dry air + water vapor).  Do appropriate conversion here.
      factor(:ncol,:) = 1._r8/factor(:ncol,:)
      do m = 1,pcnst
         if (cnst_type(m) == 'wet') then
            phys_state(lchnk)%q(:ncol,:,m) = factor(:ncol,:)*phys_state(lchnk)%q(:ncol,:,m)
         end if
      end do

      ! Compute geopotential height above surface - based on full pressure
      ! Note that phys_state%zi(:,plev+1) = 0 whereas zint in MPAS is surface height
      !
      call geopotential_t( &
         phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid,   phys_state(lchnk)%pint,          &
         phys_state(lchnk)%pmid,   phys_state(lchnk)%pdel,     phys_state(lchnk)%rpdel,         &
         phys_state(lchnk)%t,      phys_state(lchnk)%q(:,:,:), rairv(:,:,lchnk), gravit, zvirv, &
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
   use air_composition,    only: get_R
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

   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: tracers(:,:,:)

   integer :: index_qv,m,idx_dycore
   real(r8) :: thetak,exnerk,rhodk,tknew,thetaknew
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
   real(r8)          :: Rold(nCellsSolve,pver)
   real(r8)          :: Rnew(nCellsSolve,pver)
   real(r8)          :: qk    (thermodynamic_active_species_num,pver,nCellsSolve) !water species before physics (diagnostics)
   real(r8)          :: qktmp (nCellsSolve,pver,thermodynamic_active_species_num)
   integer           :: idx_thermo (thermodynamic_active_species_num)
   real(r8)          :: qwv(pver,nCellsSolve)                                  !water vapor before physics
   real(r8)          :: facnew, facold

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

   if (compute_energy_diags) then
     !
     ! Rnew and Rold are only needed for diagnostics purposes
     !
     do m=dry_air_species_num+1,thermodynamic_active_species_num
       idx_thermo(m) = m
       idx_dycore                   = thermodynamic_active_species_idx_dycore(m)
       do iCell = 1, nCellsSolve
         do k = 1, pver
           qktmp(iCell,k,m)         = tracers(idx_dycore,k,iCell)
         end do
       end do
     end do
     call get_R(qktmp,idx_thermo,Rnew)
     Rnew = Rnew*cv/Rgas

     do m=dry_air_species_num+1,thermodynamic_active_species_num
       idx_dycore    = thermodynamic_active_species_idx_dycore(m)
       do iCell = 1, nCellsSolve
         do k = 1, pver
           qktmp(iCell,k,m)  = tracers(idx_dycore,k,iCell)-dtime*q_tend(m,k,iCell)
         end do
       end do
     end do
     call get_R(qktmp,idx_thermo,Rold)
     Rold=Rold*cv/Rgas
   else
     Rnew = 0.0_r8
     Rold = 1.0_r8
   end if
   !
   ! Compute q not updated by physics
   !
   qwv = tracers(index_qv,:,1:nCellsSolve)-dtime*q_tend(index_qv_phys,:,1:nCellsSolve)
   !
   ! for energy diagnostics compute state with physics tendency (no water change) first
   ! and then add water changes (parameterizations + dme_adjust)
   !
   do iCell = 1, nCellsSolve
     do k = 1, pver
       rhodk     = zz(k,iCell) * rho_zz(k,iCell)
       facold    = 1.0_r8 + Rv_over_Rd *qwv(k,iCell)
       thetak    = theta_m(k,iCell)/facold
       exnerk    = (rgas*rhodk*theta_m(k,iCell)/p0)**(rgas/cv)
       !
       ! for compute_energy_diags only
       !
       tknew     = exnerk*thetak+(cp/Rold(iCell,k))*(Rnew(iCell,k)/cp)*dtime*t_tend(k,icell)!for diags only
       thetaknew = (tknew**(cv/cp))*((rgas*rhodk*facold)/p0)**(-rgas/cp)                    !for diags only
       !
       ! calculate theta_m tendency due to parameterizations (but no water adjustment)
       ! (for diagnostics only)
       !
       rtheta_param(k,iCell) = (thetaknew-thetak)/dtime                                     !for diags only
       rtheta_param(k,iCell) = rtheta_param(k,iCell)*(1.0_r8 + Rv_over_Rd *qwv(k,iCell))    !for diags only
       !convert to thetam
       rtheta_param(k,iCell) = rtheta_param(k,iCell)*rho_zz(k,iCell)                        !for diags only
       !
       ! include water change in theta_m
       !
       facnew               = 1.0_r8 + Rv_over_Rd *tracers(index_qv,k,iCell)
       tknew                = exnerk*thetak+dtime*t_tend(k,icell)
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
     do m=dry_air_species_num+1,thermodynamic_active_species_num
       idx_dycore                         = thermodynamic_active_species_idx_dycore(m)
       qk(m,:,: )                         = tracers(idx_dycore,:,1:nCellsSolve)
       tracers(idx_dycore,:,1:nCellsSolve)= qk(m,:,: )-dtime*q_tend(m,:,1:nCellsSolve)
     end do

     call tot_energy_dyn( &
          nCellsSolve, plev, size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), rho_zz(:,1:nCellsSolve), &
          theta_m_new,  tracers(:,:,1:nCellsSolve),   &
          ux(:,1:nCellsSolve)+dtime*u_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),       &
          uy(:,1:nCellsSolve)+dtime*v_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),'dAP')
     ! revert
     do m=dry_air_species_num+1,thermodynamic_active_species_num
       idx_dycore                         = thermodynamic_active_species_idx_dycore(m)
       tracers(idx_dycore,:,1:nCellsSolve)= qk(m,:,: )
     end do
     !
     ! compute energy incl. water change
     !
     theta_m_new = theta_m(:,1:nCellsSolve)+dtime*rtheta_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve)
     call tot_energy_dyn( &
          nCellsSolve, plev, size(tracers, 1), index_qv, zz(:,1:nCellsSolve), zint(:,1:nCellsSolve), &
          rho_zz(:,1:nCellsSolve), theta_m_new, tracers(:,:,1:nCellsSolve),    &
          ux(:,1:nCellsSolve)+dtime*u_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),       &
          uy(:,1:nCellsSolve)+dtime*v_tend(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve),'dAM')
   end if
   !
   ! compute energy based on parameterization increment (excl. water change)
   !
   theta_m_new = theta_m(:,1:nCellsSolve)+dtime*rtheta_param(:,1:nCellsSolve)/rho_zz(:,1:nCellsSolve)

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
subroutine hydrostatic_pressure(nCells, nVertLevels, qsize, index_qv, zz, zgrid, rho_zz, theta_m, &
     exner, q, pmiddry, pintdry,pmid)
   ! Compute dry hydrostatic pressure at layer interfaces and midpoints
   !
   ! Given arrays of zz, zgrid, rho_zz, and theta_m from the MPAS-A prognostic
   ! state, compute dry hydrostatic pressure at layer interfaces and midpoints.
   ! The vertical dimension for 3-d arrays is innermost, and k=1 represents
   ! the lowest layer or level in the fields.
   !
   use mpas_constants, only: cp, rgas, cv, gravity, p0, Rv_over_Rd => rvord

   ! Arguments
   integer, intent(in) :: nCells
   integer, intent(in) :: nVertLevels
   integer, intent(in) :: qsize
   integer, intent(in) :: index_qv
   real(r8), dimension(nVertLevels, nCells),       intent(in) :: zz      ! d(zeta)/dz [-]
   real(r8), dimension(nVertLevels+1, nCells),     intent(in) :: zgrid   ! geometric heights of layer interfaces [m]
   real(r8), dimension(nVertLevels, nCells),       intent(in) :: rho_zz  ! dry density / zz [kg m^-3]
   real(r8), dimension(nVertLevels, nCells),       intent(in) :: theta_m ! modified potential temperature
   real(r8), dimension(nVertLevels, nCells),       intent(in) :: exner   ! Exner function
   real(r8), dimension(qsize,nVertLevels, nCells), intent(in) :: q       ! water vapor dry mixing ratio
   real(r8), dimension(nVertLevels, nCells),       intent(out):: pmiddry ! layer midpoint dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels+1, nCells),     intent(out):: pintdry ! layer interface dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels, nCells),       intent(out):: pmid    ! layer midpoint hydrostatic pressure [Pa]

   ! Local variables
   integer :: iCell, k, idx
   real(r8), dimension(nVertLevels)          :: dz       ! Geometric layer thickness in column
   real(r8), dimension(nVertLevels)          :: dp,dpdry ! Pressure thickness
   real(r8), dimension(nVertLevels+1,nCells) :: pint  ! hydrostatic pressure at interface
   real(r8) :: sum_water
   real(r8) :: pk,rhok,rhodryk,thetavk,kap1,kap2,tvk,tk
   !
   ! For each column, integrate downward from model top to compute dry hydrostatic pressure at layer
   ! midpoints and interfaces. The pressure averaged to layer midpoints should be consistent with
   ! the ideal gas law using the rho_zz and theta values prognosed by MPAS at layer midpoints.
   !
   do iCell = 1, nCells
      dz(:) = zgrid(2:nVertLevels+1,iCell) - zgrid(1:nVertLevels,iCell)
      do k = nVertLevels, 1, -1
        rhodryk  = zz(k,iCell)* rho_zz(k,iCell) !full CAM physics density
        rhok = 1.0_r8
        do idx=dry_air_species_num+1,thermodynamic_active_species_num
          rhok = rhok+q(thermodynamic_active_species_idx_dycore(idx),k,iCell)
        end do
        rhok     = rhok*rhodryk
        dp(k)    = gravit*dz(k)*rhok
        dpdry(k) = gravit*dz(k)*rhodryk
      end do

      k = nVertLevels
      sum_water = 1.0_r8
      do idx=dry_air_species_num+1,thermodynamic_active_species_num
        sum_water = sum_water+q(thermodynamic_active_species_idx_dycore(idx),k,iCell)
      end do
      rhok     = sum_water*zz(k,iCell) * rho_zz(k,iCell)
      thetavk  = theta_m(k,iCell)/sum_water
      tvk      = thetavk*exner(k,iCell)
      pk       = dp(k)*rgas*tvk/(gravit*dz(k))
      !
      ! model top pressure consistently diagnosed using the assumption that the mid level
      ! is at height z(nVertLevels-1)+0.5*dz
      !
      pintdry(nVertLevels+1,iCell) = pk-0.5_r8*dz(nVertLevels)*rhok*gravity  !hydrostatic
      pint   (nVertLevels+1,iCell) = pintdry(nVertLevels+1,iCell)
      do k = nVertLevels, 1, -1
        !
        ! compute hydrostatic dry interface pressure so that (pintdry(k+1)-pintdry(k))/g is pseudo density
        !
        sum_water = 1.0_r8
        do idx=dry_air_species_num+1,thermodynamic_active_species_num
          sum_water = sum_water+q(thermodynamic_active_species_idx_dycore(idx),k,iCell)
        end do
        thetavk = theta_m(k,iCell)/sum_water!convert modified theta to virtual theta
        tvk     = thetavk*exner(k,iCell)
        tk      = tvk*sum_water/(1.0_r8+Rv_over_Rd*q(index_qv,k,iCell))
        pint   (k,iCell) = pint   (k+1,iCell)+dp(k)
        pintdry(k,iCell) = pintdry(k+1,iCell)+dpdry(k)
        pmid(k,iCell)    = dp(k)   *rgas*tvk/(gravit*dz(k))
        pmiddry(k,iCell) = dpdry(k)*rgas*tk /(gravit*dz(k))
      end do
    end do
end subroutine hydrostatic_pressure

subroutine tot_energy_dyn(nCells, nVertLevels, qsize, index_qv, zz, zgrid, rho_zz, theta_m, q, ux,uy,outfld_name_suffix)
  use physconst,         only: rair, gravit
  use mpas_constants,    only: p0,cv,rv,rgas,cp
  use cam_history,       only: outfld, hist_fld_active
  use mpas_constants,    only: Rv_over_Rd => rvord
  use air_composition,   only: thermodynamic_active_species_ice_idx_dycore,thermodynamic_active_species_liq_idx_dycore
  use air_composition,   only: thermodynamic_active_species_ice_num,thermodynamic_active_species_liq_num
  use air_composition,   only: dry_air_species_num, thermodynamic_active_species_R
  use cam_thermo,        only: wvidx,wlidx,wiidx,seidx,poidx,keidx,teidx,thermo_budget_num_vars
  use cam_thermo,        only: get_hydrostatic_energy,thermo_budget_vars
  use dyn_tests_utils,   only: vcoord=>vc_height
  use cam_history_support,    only: max_fieldname_len

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
  integer :: iCell, k, idx, idx_tmp
  integer :: i
  real(r8) :: rho_dz,theta,pk,ptop,exner,dz,rhod
  real(r8), dimension(nCells,nVertLevels)       :: temperature, pdeldry, cp_or_cv, zcell, u, v
  real(r8), dimension(nCells)                   :: phis
  real(r8), dimension(nCells,nVertLevels,qsize) :: tracers
  real(r8), dimension(nCells)                   :: kinetic_energy,potential_energy,internal_energy,water_vapor

  real(r8), dimension(nCells) :: liq !total column integrated liquid
  real(r8), dimension(nCells) :: ice !total column integrated ice
  real(r8) :: sum_species

  character(len=max_fieldname_len) :: name_out(thermo_budget_num_vars)


  do i=1,thermo_budget_num_vars
     name_out(i)=trim(thermo_budget_vars(i))//'_'//trim(outfld_name_suffix)
  end do

  kinetic_energy   = 0.0_r8
  potential_energy = 0.0_r8
  internal_energy  = 0.0_r8
  water_vapor      = 0.0_r8
  tracers          = 0.0_r8

  do iCell = 1, nCells
     do k = 1, nVertLevels
        dz              = zgrid(k+1,iCell) - zgrid(k,iCell)
        zcell(iCell,k)  = 0.5_r8*(zgrid(k,iCell)+zgrid(k+1,iCell))-zgrid(1,iCell)
        rhod            = zz(k,iCell) * rho_zz(k,iCell)
        theta           = theta_m(k,iCell)/(1.0_r8 + Rv_over_Rd *q(index_qv,k,iCell))!convert theta_m to theta
        exner           = (rgas*rhod*theta_m(k,iCell)/p0)**(rgas/cv)

        temperature(iCell,k)   = exner*theta
        pdeldry(iCell,k)       = gravit*rhod*dz
        !
        ! internal energy coefficient for MPAS
        ! (equation 92 in Eldred et al. 2023; https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.4353)
        !
        cp_or_cv(iCell,k)      = rair
        sum_species            = 1.0_r8
        do idx=dry_air_species_num + 1,thermodynamic_active_species_num
          idx_tmp = thermodynamic_active_species_idx_dycore(idx)
          cp_or_cv(iCell,k) = cp_or_cv(iCell,k)+thermodynamic_active_species_R(idx)*q(idx_tmp,k,iCell)
          sum_species       = sum_species+q(idx_tmp,k,iCell)
        end do
        cp_or_cv(iCell,k)      = cv*cp_or_cv(iCell,k)/(sum_species*rair)
        u(iCell,k)             = ux(k,iCell)
        v(iCell,k)             = uy(k,iCell)
        phis(iCell)            = zgrid(1,iCell)*gravit
        do idx=dry_air_species_num+1,thermodynamic_active_species_num
           idx_tmp = thermodynamic_active_species_idx_dycore(idx)
           tracers(iCell,k,idx_tmp) = q(idx_tmp,k,iCell)
        end do
     end do
  enddo
  call get_hydrostatic_energy(tracers, .false., pdeldry, cp_or_cv, u, v, temperature, &
       vcoord=vcoord, phis = phis, z_mid=zcell, dycore_idx=.true.,                    &
       se=internal_energy, po=potential_energy, ke=kinetic_energy,                    &
       wv=water_vapor    , liq=liq             , ice=ice)

  call outfld(name_out(seidx),internal_energy ,ncells,1)
  call outfld(name_out(poidx),potential_energy,ncells,1)
  call outfld(name_out(keidx),kinetic_energy  ,ncells,1)
  call outfld(name_out(wvidx),water_vapor     ,ncells,1)
  call outfld(name_out(wlidx),liq             ,ncells,1)
  call outfld(name_out(wiidx),ice             ,ncells,1)
  call outfld(name_out(teidx),potential_energy+internal_energy+kinetic_energy,ncells,1)

end subroutine tot_energy_dyn

 subroutine calc_frontogenesis( frontogenesisFunction, frontogenesisAngle,             &
      theta_m, qv, u,v, defc_a, defc_b, cell_gradient_coef_x, cell_gradient_coef_y, &
      areaCell, dvEdge, cellsOnEdge, edgesOnCell, nEdgesOnCell, edgesOnCell_sign,   &
      nVertLevels, nCellsSolve )

   use mpas_constants, only: rvord

   ! inputs

   integer, intent(in) :: nVertLevels, nCellsSolve
   real(r8), dimension(:,:), intent(in) :: theta_m, qv
   real(r8), dimension(:,:), intent(in) :: u, v
   real(r8), dimension(:,:), intent(in) :: defc_a
   real(r8), dimension(:,:), intent(in) :: defc_b
   real(r8), dimension(:,:), intent(in) :: cell_gradient_coef_x
   real(r8), dimension(:,:), intent(in) :: cell_gradient_coef_y
   real(r8), dimension(:,:), intent(in) :: edgesOnCell_sign
   real(r8), dimension(:), intent(in) :: dvEdge
   real(r8), dimension(:), intent(in) :: areaCell
   integer, dimension(:,:), intent(in) :: cellsOnEdge
   integer, dimension(:,:), intent(in) :: edgesOnCell
   integer, dimension(:), intent(in) :: nEdgesOnCell

   ! outputs

   real(r8), dimension(:,:), intent(out) :: frontogenesisFunction(:,:)
   real(r8), dimension(:,:), intent(out) :: frontogenesisAngle(:,:)

   ! local storage

   integer :: iCell, iEdge, k, cell1, cell2
   real(r8), dimension(nVertLevels) :: d_diag, d_off_diag, divh, theta_x, theta_y
   real(r8) :: edge_sign, thetaEdge

   !
   ! for each column, compute frontogenesis function and del(theta) angle
   !

   do iCell = 1,nCellsSolve

      d_diag(1:nVertLevels) = 0.0_r8
      d_off_diag(1:nVertLevels) = 0.0_r8
      divh(1:nVertLevels) = 0.0_r8
      theta_x(1:nVertLevels) = 0.0_r8
      theta_y(1:nVertLevels) = 0.0_r8

      !
      ! Integrate over edges to compute cell-averaged divergence, deformation,
      ! d(theta)/dx, and d(theta)/dy.  (x,y) are aligned with (lon,lat) at the
      ! cell center in the 2D tangent-plane approximation used here.  This alignment
      ! is set in the initialization routine for the coefficients
      ! defc_a, defc_b, cell_gradient_coef_x and cell_gradient_coef_y that is
      ! part of the MPAS mesh initialization.  The horizontal divergence is calculated
      ! as it is in the MPAS solver, i.e. on the sphere as opposed to on the tangent plane.
      !
      do iEdge=1,nEdgesOnCell(iCell)

         edge_sign = edgesOnCell_sign(iEdge,iCell) * dvEdge(edgesOnCell(iEdge,iCell)) / areaCell(iCell)
         cell1 = cellsOnEdge(1,edgesOnCell(iEdge,iCell))
         cell2 = cellsOnEdge(2,edgesOnCell(iEdge,iCell))

         do k=1,nVertLevels

            d_diag(k)     = d_diag(k)     + defc_a(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                          - defc_b(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
            d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                          + defc_a(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
            divh(k) = divh(k) + edge_sign * u(k,EdgesOnCell(iEdge,iCell))
            thetaEdge = 0.5_r8*( theta_m(k,cell1)/(1.0_r8 + rvord*qv(k,cell1))  &
                                +theta_m(k,cell2)/(1.0_r8 + rvord*qv(k,cell2)) )
            theta_x(k) = theta_x(k) + cell_gradient_coef_x(iEdge,iCell)*thetaEdge
            theta_y(k) = theta_y(k) + cell_gradient_coef_y(iEdge,iCell)*thetaEdge

         end do

      end do

      !
      ! compute the frontogenesis function:
      !  1/2  |del(theta)/dt)| = 1/2 (
      !                - Div * |del(theta)|^2
      !                - E  (d(theta)/dx)^2
      !                - 2F (d(theta)/dx)*(d(theta)/dy)
      !                + E  (d(theta)/dy)  )
      ! where
      !        Div = u_x + v_y (horizontal velocity divergence)
      !        E = u_x - v_y (stretching deformation)
      !        F = v_x + u_y (shearing deformation)
      !
      do k=1, nVertLevels

         frontogenesisFunction(k,iCell) = 0.5_r8*(         &
              -divh(k)*(theta_x(k)**2 + theta_y(k)**2)  &
              -d_diag(k)*theta_x(k)**2                  &
              -2.0_r8*d_off_diag(k)*theta_x(k)*theta_y(k)  &
              +d_diag(k)*theta_y(k)**2                )
         frontogenesisAngle(k,iCell) = atan2(theta_y(k),theta_x(k))

      end do

   end do

 end subroutine calc_frontogenesis

end module dp_coupling
