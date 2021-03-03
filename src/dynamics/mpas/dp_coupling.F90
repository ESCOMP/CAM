module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use pmgrid,         only: plev
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_type
use physconst,      only: gravit, cpairv, cappa, rairv, rh2o, zvir

use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
use spmd_utils,     only: mpicom, iam, masterproc

use dyn_grid,       only: get_gcol_block_d
use dyn_comp,       only: dyn_export_t, dyn_import_t

use physics_types,  only: physics_state, physics_tend
use phys_grid,      only: get_ncols_p, get_gcol_all_p, block_to_chunk_send_pters, &
                          transpose_block_to_chunk, block_to_chunk_recv_pters,    &
                          chunk_to_block_send_pters, transpose_chunk_to_block,    &
                          chunk_to_block_recv_pters

use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field

use cam_logfile,    only: iulog
use perf_mod,       only: t_startf, t_stopf, t_barrierf
use cam_abortutils, only: endrun

implicit none
private
save

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

   real(r8), pointer :: pmiddry(:,:)
   real(r8), pointer :: pintdry(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8), pointer :: w(:,:)
   real(r8), pointer :: theta_m(:,:)
   real(r8), pointer :: exner(:,:)
   real(r8), pointer :: tracers(:,:,:)


   integer :: lchnk, icol, k, kk      ! indices over chunks, columns, layers
   integer :: i, m, ncols, blockid
   integer :: blk(1), bcid(1)

   integer :: pgcols(pcols)
   integer :: tsize                    ! amount of data per grid point passed to physics
   integer, allocatable :: bpter(:,:)  ! offsets into block buffer for packing data
   integer, allocatable :: cpter(:,:)  ! offsets into chunk buffer for unpacking data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   character(len=*), parameter :: subname = 'd_p_coupling'
   !----------------------------------------------------------------------------

   nCellsSolve = dyn_out % nCellsSolve
   index_qv    = dyn_out % index_qv
   cam_from_mpas_cnst => dyn_out % cam_from_mpas_cnst

   pmiddry  => dyn_out % pmiddry
   pintdry  => dyn_out % pintdry
   zint     => dyn_out % zint
   zz       => dyn_out % zz
   rho_zz   => dyn_out % rho_zz
   ux       => dyn_out % ux
   uy       => dyn_out % uy
   w        => dyn_out % w
   theta_m  => dyn_out % theta_m
   exner    => dyn_out % exner
   tracers  => dyn_out % tracers

   ! diagnose pintdry, pmiddry
   call dry_hydrostatic_pressure( &
      nCellsSolve, plev, zz, zint, rho_zz, theta_m, pmiddry, pintdry)

   call t_startf('dpcopy')

   if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, icol, i, k, kk, m, pgcols, blk, bcid)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)                         ! number of columns in this chunk
         call get_gcol_all_p(lchnk, pcols, pgcols)          ! global column indices in chunk

         do icol = 1, ncols                                   ! column index in physics chunk
            call get_gcol_block_d(pgcols(icol), 1, blk, bcid) ! column index in dynamics block
            i = bcid(1)

            phys_state(lchnk)%psdry(icol) = pintdry(1,i)
            phys_state(lchnk)%phis(icol) = zint(1,i) * gravit

            do k = 1, pver                                  ! vertical index in physics chunk
               kk = pver - k + 1                            ! vertical index in dynamics block

               phys_state(lchnk)%t(icol,k)       = theta_m(kk,i) / (1.0_r8 + &
                                                   Rv_over_Rd * tracers(index_qv,kk,i)) * exner(kk,i)
               phys_state(lchnk)%u(icol,k)       = ux(kk,i)
               phys_state(lchnk)%v(icol,k)       = uy(kk,i)
               phys_state(lchnk)%omega(icol,k)   = -rho_zz(kk,i)*zz(kk,i)*gravit*0.5_r8*(w(kk,i)+w(kk+1,i))   ! omega
               phys_state(lchnk)%pmiddry(icol,k) = pmiddry(kk,i)
            end do

            do k = 1, pverp
               kk = pverp - k + 1
               phys_state(lchnk)%pintdry(icol,k) = pintdry(kk,i)
            end do

            do m = 1, pcnst
               do k = 1, pver
                  kk = pver - k + 1
                  phys_state(lchnk)%q(icol,k,m) = tracers(cam_from_mpas_cnst(m),kk,i)
               end do
            end do
         end do
      end do

   else  ! .not. local_dp_map

      tsize = 6 + pcnst
      allocate(bbuffer(tsize*block_buf_nrecs))    ! block buffer
      bbuffer = 0.0_r8
      allocate(cbuffer(tsize*chunk_buf_nrecs))    ! chunk buffer
      cbuffer = 0.0_r8

      allocate( bpter(nCellsSolve,0:pver) )
      allocate( cpter(pcols,0:pver) )

      blockid = iam + 1   ! global block index
      call block_to_chunk_send_pters(blockid, nCellsSolve, pverp, tsize, bpter)

      do i = 1, nCellsSolve                           ! column index in block

         bbuffer(bpter(i,0))   = pintdry(1,i)         ! psdry
         bbuffer(bpter(i,0)+1) = zint(1,i) * gravit   ! phis

         do k = 1, pver
            bbuffer(bpter(i,k))   = theta_m(k,i) / (1.0_r8 + &
                                       Rv_over_Rd * tracers(index_qv,k,i)) * exner(k,i)
            bbuffer(bpter(i,k)+1) = ux(k,i)
            bbuffer(bpter(i,k)+2) = uy(k,i)
            bbuffer(bpter(i,k)+3) = -rho_zz(k,i) * zz(k,i) * gravit * 0.5_r8 * (w(k,i) + w(k+1,i))   ! omega
            bbuffer(bpter(i,k)+4) = pmiddry(k,i)
            do m=1,pcnst
               bbuffer(bpter(i,k)+4+m) = tracers(cam_from_mpas_cnst(m),k,i)
            end do
         end do

         do k = 1, pverp
            bbuffer(bpter(i,k-1)+5+pcnst) = pintdry(k,i)
         end do
      end do

      call t_barrierf ('sync_blk_to_chk', mpicom)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')

      !$omp parallel do private (lchnk, ncols, icol, k, kk, m, cpter)
      do lchnk = begchunk, endchunk
         ncols = phys_state(lchnk)%ncol

         call block_to_chunk_recv_pters(lchnk, pcols, pverp, tsize, cpter)

         do icol = 1, ncols

            phys_state(lchnk)%psdry(icol)    = cbuffer(cpter(icol,0))
            phys_state(lchnk)%phis(icol)     = cbuffer(cpter(icol,0)+1)

            ! do the vertical reorder here when assigning to phys_state
            do k = 1, pver
               kk = pver - k + 1
               phys_state(lchnk)%t      (icol,kk)  = cbuffer(cpter(icol,k))
               phys_state(lchnk)%u      (icol,kk)  = cbuffer(cpter(icol,k)+1)
               phys_state(lchnk)%v      (icol,kk)  = cbuffer(cpter(icol,k)+2)
               phys_state(lchnk)%omega  (icol,kk)  = cbuffer(cpter(icol,k)+3)
               phys_state(lchnk)%pmiddry(icol,kk)  = cbuffer(cpter(icol,k)+4)
               do m = 1, pcnst
                  phys_state(lchnk)%q  (icol,kk,m) = cbuffer(cpter(icol,k)+4+m)
               end do
            end do

            do k = 0, pver
               kk = pverp - k
               phys_state(lchnk)%pintdry(icol,kk) = cbuffer(cpter(icol,k)+5+pcnst)
            end do
         end do
      end do

      deallocate( bbuffer, bpter )
      deallocate( cbuffer, cpter )

   end if
   call t_stopf('dpcopy')

   call t_startf('derived_phys')
   call derived_phys(phys_state, phys_tend, pbuf2d)
   call t_stopf('derived_phys')

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
   integer :: lchnk, icol, k, kk      ! indices over chunks, columns, layers
   integer :: i, m, ncols, blockid
   integer :: blk(1), bcid(1)

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
   real(r8), allocatable :: qv_tend(:,:)
   real(r8), pointer :: u_tend(:,:)
   real(r8), pointer :: v_tend(:,:)


   integer :: pgcols(pcols)
   integer :: tsize                    ! amount of data per grid point passed to dynamics
   integer, allocatable :: bpter(:,:)  ! offsets into block buffer for unpacking data
   integer, allocatable :: cpter(:,:)  ! offsets into chunk buffer for packing data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   type (mpas_pool_type), pointer :: tend_physics
   type (field2DReal), pointer :: tend_uzonal, tend_umerid

   character(len=*), parameter :: subname = 'dp_coupling::p_d_coupling'
   !----------------------------------------------------------------------------

   nCellsSolve = dyn_in % nCellsSolve
   nCells      = dyn_in % nCells
   index_qv    = dyn_in % index_qv
   mpas_from_cam_cnst => dyn_in % mpas_from_cam_cnst

   tracers => dyn_in % tracers

   allocate( t_tend(pver,nCellsSolve) )
   allocate( qv_tend(pver,nCellsSolve) )

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
   if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, icol, i, k, kk, m, pgcols, blk, bcid)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)                     ! number of columns in this chunk
         call get_gcol_all_p(lchnk, pcols, pgcols)      ! global column indices

         do icol = 1, ncols                                   ! column index in physics chunk
            call get_gcol_block_d(pgcols(icol), 1, blk, bcid) ! column index in dynamics block
            i = bcid(1)
            
            do k = 1, pver                              ! vertical index in physics chunk
               kk = pver - k + 1                        ! vertical index in dynamics block

               t_tend(kk,i) = phys_tend(lchnk)%dtdt(icol,k)
               u_tend(kk,i) = phys_tend(lchnk)%dudt(icol,k)
               v_tend(kk,i) = phys_tend(lchnk)%dvdt(icol,k)

               ! convert wet mixing ratios to dry
               factor = phys_state(lchnk)%pdel(icol,k)/phys_state(lchnk)%pdeldry(icol,k)
               do m = 1, pcnst
                  if (cnst_type(mpas_from_cam_cnst(m)) == 'wet') then
                     if (m == index_qv) then
                        qv_tend(kk,i) = (phys_state(lchnk)%q(icol,k,mpas_from_cam_cnst(m))*factor - tracers(index_qv,kk,i)) / dt_phys
                     end if
                     tracers(m,kk,i) = phys_state(lchnk)%q(icol,k,mpas_from_cam_cnst(m))*factor
                  else
                     if (m == index_qv) then
                        qv_tend(kk,i) = (phys_state(lchnk)%q(icol,k,mpas_from_cam_cnst(m)) - tracers(index_qv,kk,i)) / dt_phys
                     end if
                     tracers(m,kk,i) = phys_state(lchnk)%q(icol,k,mpas_from_cam_cnst(m))
                  end if
               end do

            end do
         end do
      end do

   else

      tsize = 3 + pcnst
      allocate( bbuffer(tsize*block_buf_nrecs) )
      bbuffer = 0.0_r8
      allocate( cbuffer(tsize*chunk_buf_nrecs) )
      cbuffer = 0.0_r8

      allocate( bpter(nCellsSolve,0:pver) )
      allocate( cpter(pcols,0:pver) )

      !$omp parallel do private (lchnk, ncols, icol, k, m, cpter)
      do lchnk = begchunk, endchunk
         ncols = get_ncols_p(lchnk)

         call chunk_to_block_send_pters(lchnk, pcols, pverp, tsize, cpter)

         do icol = 1, ncols

            do k = 1, pver
               cbuffer(cpter(icol,k))   = phys_tend(lchnk)%dtdt(icol,k)
               cbuffer(cpter(icol,k)+1) = phys_tend(lchnk)%dudt(icol,k)
               cbuffer(cpter(icol,k)+2) = phys_tend(lchnk)%dvdt(icol,k)

               ! convert wet mixing ratios to dry
               factor = phys_state(lchnk)%pdel(icol,k)/phys_state(lchnk)%pdeldry(icol,k)
               do m = 1, pcnst
                  if (cnst_type(m) == 'wet') then
                     cbuffer(cpter(icol,k)+2+m) = phys_state(lchnk)%q(icol,k,m)*factor
                  else
                     cbuffer(cpter(icol,k)+2+m) = phys_state(lchnk)%q(icol,k,m)
                  end if
               end do

            end do
         end do
      end do

      call t_barrierf('sync_chk_to_blk', mpicom)
      call t_startf ('chunk_to_block')
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      call t_stopf  ('chunk_to_block')

      blockid = iam + 1   !  global block index
   
      call chunk_to_block_recv_pters(blockid, nCellsSolve, pverp, tsize, bpter)

      do i = 1, nCellsSolve                       ! index in dynamics block

         ! flip vertical index here
         do k = 1, pver                        ! vertical index in physics chunk
            kk = pver - k + 1                  ! vertical index in dynamics block

            t_tend(kk,i) = bbuffer(bpter(i,k))
            u_tend(kk,i) = bbuffer(bpter(i,k)+1)
            v_tend(kk,i) = bbuffer(bpter(i,k)+2)

            do m = 1, pcnst
               if (m == index_qv) then
                  qv_tend(kk,i) = (bbuffer(bpter(i,k)+2+mpas_from_cam_cnst(m)) - tracers(index_qv,kk,i)) / dt_phys
               end if
               tracers(m,kk,i) = bbuffer(bpter(i,k)+2+mpas_from_cam_cnst(m))
            end do

         end do
      end do

      deallocate( bbuffer, bpter )
      deallocate( cbuffer, cpter )

   end if
   call t_stopf('pd_copy')

   call t_startf('derived_tend')
   call derived_tend(nCellsSolve, nCells, t_tend, u_tend, v_tend, qv_tend, dyn_in)
   call t_stopf('derived_tend')

   call mpas_deallocate_scratch_field(tend_uzonal)
   call mpas_deallocate_scratch_field(tend_umerid)

end subroutine p_d_coupling

!=========================================================================================

subroutine derived_phys(phys_state, phys_tend, pbuf2d)

   ! Compute fields in the physics state object which are diagnosed from the
   ! MPAS prognostic fields.

   use geopotential,  only: geopotential_t
   use check_energy,  only: check_energy_timestep_init
   use shr_vmath_mod, only: shr_vmath_log

   ! Arguments
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! Local variables

   integer :: k, lchnk, m, ncol

   real(r8) :: factor(pcols,pver)
   real(r8) :: zvirv(pcols,pver)

   real(r8), parameter :: pref = 1.e5_r8 ! reference pressure (Pa)

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   character(len=*), parameter :: subname = 'dp_coupling::derived_phys'
   !----------------------------------------------------------------------------

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
         phys_state(lchnk)%pmid(:ncol,k) = (phys_state(lchnk)%pint(:ncol,k+1) + &
                                            phys_state(lchnk)%pint(:ncol,k)) / 2._r8

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

      ! fill zvirv 2D variables to be compatible with geopotential_t interface
      zvirv(:,:) = zvir

      ! Compute geopotential height above surface - based on full pressure
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

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

   end do

end subroutine derived_phys

!=========================================================================================

subroutine derived_tend(nCellsSolve, nCells, t_tend, u_tend, v_tend, qv_tend, dyn_in)

   ! Derive the physics tendencies required by MPAS from the tendencies produced by
   ! CAM's physics package.

   use cam_mpas_subdriver, only : cam_mpas_cell_to_edge_winds, cam_mpas_update_halo
   use mpas_constants, only : Rv_over_Rd => rvord

   ! Arguments
   integer,             intent(in)    :: nCellsSolve
   integer,             intent(in)    :: nCells
   real(r8),            intent(in)    :: t_tend(pver,nCellsSolve)  ! physics dtdt
   real(r8),            intent(in)    :: qv_tend(pver,nCellsSolve) ! physics dqvdt
   real(r8),            intent(inout) :: u_tend(pver,nCells+1)     ! physics dudt
   real(r8),            intent(inout) :: v_tend(pver,nCells+1)     ! physics dvdt
   type(dyn_import_t),  intent(inout) :: dyn_in


   ! Local variables

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

   integer :: index_qv

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


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Temperature tendency
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

   !
   ! Convert temperature tendency to potential temperature tendency
   !
   rtheta_tend(:,1:nCellsSolve) = t_tend(:,1:nCellsSolve) / exner(:,1:nCellsSolve)

   !
   ! Couple theta tendency with rho_zz
   !
   rtheta_tend(:,1:nCellsSolve) = rtheta_tend(:,1:nCellsSolve) * rho_zz(:,1:nCellsSolve)

   !
   ! Modify with moisture terms
   !
   rtheta_tend(:,1:nCellsSolve) = rtheta_tend(:,1:nCellsSolve) * (1.0_r8 + Rv_over_Rd * tracers(index_qv,:,1:nCellsSolve))
   rtheta_tend(:,1:nCellsSolve) = rtheta_tend(:,1:nCellsSolve) + Rv_over_Rd * theta(:,1:nCellsSolve) * qv_tend(:,1:nCellsSolve)

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

subroutine dry_hydrostatic_pressure(nCells, nVertLevels, zz, zgrid, rho_zz, theta_m, pmiddry, pintdry)

   ! Compute dry hydrostatic pressure at layer interfaces and midpoints
   !
   ! Given arrays of zz, zgrid, rho_zz, and theta_m from the MPAS-A prognostic
   ! state, compute dry hydrostatic pressure at layer interfaces and midpoints.
   ! The vertical dimension for 3-d arrays is innermost, and k=1 represents
   ! the lowest layer or level in the fields.
   !
   ! IMPORTANT NOTE: At present, this routine is probably not correct when there
   !                 is moisture in the atmosphere.

   use mpas_constants, only : cp, rgas, cv, gravity, p0

   ! Arguments
   integer, intent(in) :: nCells
   integer, intent(in) :: nVertLevels
   real(r8), dimension(nVertLevels, nCells), intent(in) :: zz         ! d(zeta)/dz [-]
   real(r8), dimension(nVertLevels+1, nCells), intent(in) :: zgrid    ! geometric heights of layer interfaces [m]
   real(r8), dimension(nVertLevels, nCells), intent(in) :: rho_zz     ! dry density / zz [kg m^-3]
   real(r8), dimension(nVertLevels, nCells), intent(in) :: theta_m    ! potential temperature * (1 + Rv/Rd * qv)
   real(r8), dimension(nVertLevels, nCells), intent(out) :: pmiddry   ! layer midpoint dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels+1, nCells), intent(out) :: pintdry ! layer interface dry hydrostatic pressure [Pa]

   ! Local variables
   integer :: iCell, k
   real(r8), dimension(nCells) :: ptop_int   ! Extrapolated pressure at top of the model
   real(r8), dimension(nCells) :: ptop_mid   ! Full non-hydrostatic pressure at top layer midpoint
   real(r8), dimension(nCells) :: ttop_mid   ! Temperature at top layer midpoint
   real(r8), dimension(nVertLevels) :: dz    ! Geometric layer thickness in column
   real(r8) :: pi, t

   !
   ! Compute full non-hydrostatic pressure and temperature at top layer midpoint
   !
   ptop_mid(:) = p0 * (rgas * rho_zz(nVertLevels,:) * zz(nVertLevels,:) * theta_m(nVertLevels,:) / p0)**(cp/cv)
   ttop_mid(:) = theta_m(nVertLevels,:) * &
                 (zz(nVertLevels,:) * rgas * rho_zz(nVertLevels,:) * theta_m(nVertLevels,:) / p0)**(rgas/(cp-rgas))


   !
   ! Extrapolate upward from top layer midpoint to top of the model
   ! The formula used here results from combination of the hypsometric equation with the equation
   ! for the layer mid-point pressure (i.e., (pint_top + pint_bot)/2 = pmid)
   !
   ! TODO: Should temperature here be virtual temperature?
   !
   ptop_int(:) = 2.0_r8 * ptop_mid(:) &
                 / (1.0_r8 + exp( (zgrid(nVertLevels+1,:) - zgrid(nVertLevels,:)) * gravity / rgas / ttop_mid(:)))


   !
   ! For each column, integrate downward from model top to compute dry hydrostatic pressure at layer
   ! midpoints and interfaces. The pressure averaged to layer midpoints should be consistent with
   ! the ideal gas law using the rho_zz and theta_m values prognosed by MPAS at layer midpoints.
   !
   ! TODO: Should temperature here be virtual temperature?
   ! TODO: Is it problematic that the computed temperature is consistent with the non-hydrostatic pressure?
   !
   do iCell = 1, nCells

      dz(:) = zgrid(2:nVertLevels+1,iCell) - zgrid(1:nVertLevels,iCell)

      pintdry(nVertLevels+1,iCell) = ptop_int(iCell)
      do k = nVertLevels, 1, -1
         pintdry(k,iCell) = pintdry(k+1,iCell) + gravity * zz(k,iCell) * rho_zz(k,iCell) * dz(k)
         pmiddry(k,iCell) = 0.5_r8 * (pintdry(k+1,iCell) + pintdry(k,iCell))
      end do
   end do

end subroutine dry_hydrostatic_pressure

!=========================================================================================

end module dp_coupling
