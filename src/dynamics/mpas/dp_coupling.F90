#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

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

use dyn_grid,       only: max_col_per_block=>maxNCells, get_block_gcol_cnt_d, &
                          col_indices_in_block, global_to_local_cell=>local_col_index

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


integer, parameter :: nblocks_per_pe = 1

real(r8), parameter :: pref = 1.e5_r8 ! reference pressure (Pa)

!=========================================================================================
contains
!=========================================================================================

subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)

   ! Convert the dynamics output state into the physics input state.
   ! Note that all pressures and tracer mixing ratios coming from the dycore are based on
   ! dry air mass.

   ! arguments
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_export_t),        intent(inout) :: dyn_out
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! local variables

   ! Variables from dynamics export container
   integer :: nCellsSolve
   integer :: index_qv

   real(r8), pointer :: pmiddry(:,:)
   real(r8), pointer :: pintdry(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8), pointer :: w(:,:)
   real(r8), pointer :: theta_m(:,:)
   real(r8), pointer :: tracers(:,:,:)


   integer :: lchnk, icol, k, kk      ! indices over chunks, columns, layers
   integer :: ncols, ig, nblk, nb, m, i

   integer :: pgcols(pcols)
   integer :: tsize                           ! amount of data per grid point passed to physics
   integer :: bpter(max_col_per_block,0:pver) ! offsets into block buffer for packing data
   integer :: cpter(pcols,0:pver)             ! offsets into chunk buffer for unpacking data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   character(len=*), parameter :: subname = 'd_p_coupling'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   nCellsSolve = dyn_out % nCellsSolve
   index_qv    = dyn_out % index_qv

   pmiddry  => dyn_out % pmiddry
   pintdry  => dyn_out % pintdry
   zint     => dyn_out % zint
   zz       => dyn_out % zz
   rho_zz   => dyn_out % rho_zz
   ux       => dyn_out % ux
   uy       => dyn_out % uy
   w        => dyn_out % w
   theta_m  => dyn_out % theta_m
   tracers  => dyn_out % tracers

   ! diagnose pintdry, pmiddry
   call dry_hydrostatic_pressure( &
      nCellsSolve, plev, zz, zint, rho_zz, theta_m, pmiddry, pintdry)

   call t_startf('dpcopy')

   if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, icol, i, k, kk, m, pgcols)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)                         ! number of columns in this chunk
         call get_gcol_all_p(lchnk, pcols, pgcols)          ! global column indices in chunk

         do icol = 1, ncols                                 ! column index in physics chunk
            i = global_to_local_cell(pgcols(icol))          ! column index in dynamics block

            phys_state(lchnk)%psdry(icol) = pintdry(1,i)
            phys_state(lchnk)%phis(icol) = zint(1,i) * gravit

            do k = 1, pver                                  ! vertical index in physics chunk
               kk = pver - k + 1                            ! vertical index in dynamics block

               phys_state(lchnk)%t(icol,k)       = theta_m(kk,i)  ! convert to temperature in derived_phys
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
                  ! *** needs conversion of constituent indices
                  phys_state(lchnk)%q(icol,k,m) = tracers(m,kk,i)
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

      !$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, m, bpter)
      do nb = 1, nblocks_per_pe
         nblk = iam * nblocks_per_pe + nb   ! global block index
         ncols = get_block_gcol_cnt_d(nblk) ! number of columns in this block

         call block_to_chunk_send_pters(nblk, max_col_per_block, pverp, tsize, bpter)

         do icol = 1, ncols                        ! column index in physics chunk
            ig = col_indices_in_block(icol,nblk)   ! global column index
            i = global_to_local_cell(ig)           ! column index in dynamics block

            bbuffer(bpter(icol,0))   = pintdry(1,i)         ! psdry
            bbuffer(bpter(icol,0)+1) = zint(1,i) * gravit   ! phis

            do k = 1, pver
               bbuffer(bpter(icol,k))   = theta_m(k,i) ! convert to temperature in derived_phys
               bbuffer(bpter(icol,k)+1) = ux(k,i)
               bbuffer(bpter(icol,k)+2) = uy(k,i)
               bbuffer(bpter(icol,k)+3) = -rho_zz(k,i) * zz(k,i) * gravit * 0.5_r8 * (w(k,i) + w(k+1,i))   ! omega
               bbuffer(bpter(icol,k)+4) = pmiddry(k,i)
               do m=1,pcnst
                  ! *** needs conversion of constituent indices
                  bbuffer(bpter(icol,k)+4+m) = tracers(m,k,i)
               end do
            end do

            do k = 1, pverp
               bbuffer(bpter(icol,k-1)+5+pcnst) = pintdry(k,i)
            end do
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

      deallocate( bbuffer )
      deallocate( cbuffer )

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

   ! Arguments
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ), intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_import_t),  intent(inout) :: dyn_in

   ! Local variables
   integer :: lchnk, icol, k, kk      ! indices over chunks, columns, layers
   integer :: i, ig, m, nb, nblk, ncols

   real(r8) :: factor

   ! Variables from dynamics import container
   integer :: nCellsSolve
   integer :: index_qv

   real(r8), pointer :: tracers(:,:,:)

   ! CAM physics output redistributed to blocks.
   real(r8), allocatable :: t_tend(:,:)
   real(r8), allocatable :: u_tend(:,:)
   real(r8), allocatable :: v_tend(:,:)


   integer :: pgcols(pcols)
   integer :: tsize                           ! amount of data per grid point passed to dynamics
   integer :: cpter(pcols,0:pver)             ! offsets into chunk buffer for packing data
   integer :: bpter(max_col_per_block,0:pver) ! offsets into block buffer for unpacking data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   character(len=*), parameter :: subname = 'dp_coupling::p_d_coupling'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   nCellsSolve = dyn_in % nCellsSolve
   index_qv    = dyn_in % index_qv

   tracers => dyn_in % tracers

   allocate( &
      t_tend(pver,nCellsSolve), &
      u_tend(pver,nCellsSolve), &
      v_tend(pver,nCellsSolve)  )

   call t_startf('pd_copy')
   if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, icol, i, k, m, pgcols)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)                     ! number of columns in this chunk
         call get_gcol_all_p(lchnk, pcols, pgcols)      ! global column indices

         do icol = 1, ncols                             ! column index in physics chunk
            i = global_to_local_cell(pgcols(icol))      ! column index in dynamics block

            do k = 1, pver                              ! vertical index in physics chunk
               kk = pver - k + 1                        ! vertical index in dynamics block

               t_tend(kk,i) = phys_tend(lchnk)%dtdt(icol,k)
               u_tend(kk,i) = phys_tend(lchnk)%dudt(icol,k)
               v_tend(kk,i) = phys_tend(lchnk)%dvdt(icol,k)

               ! convert wet mixing ratios to dry
               factor = phys_state(lchnk)%pdel(icol,k)/phys_state(lchnk)%pdeldry(icol,k)
               do m = 1, pcnst
                  ! *** will need conversion between CAM and MPAS tracer indices
                  if (cnst_type(m) == 'wet') then
                     tracers(m,kk,i) = phys_state(lchnk)%q(icol,k,m)*factor
                  else
                     tracers(m,kk,i) = phys_state(lchnk)%q(icol,k,m)
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

      !$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, kk, m, bpter)
      do nb = 1, nblocks_per_pe
         nblk = iam * nblocks_per_pe + nb   !  global block index
         ncols = get_block_gcol_cnt_d(nblk) !  number of columns in this block

         call chunk_to_block_recv_pters(nblk, max_col_per_block, pverp, tsize, bpter)

         do icol = 1, ncols                       ! column index in physics chunk
            ig = col_indices_in_block(icol,nblk)  ! global column index
            i = global_to_local_cell(ig)          ! column index in dynamics block

            ! flip vertical index here
            do k = 1, pver                        ! vertical index in physics chunk
               kk = pver - k + 1                  ! vertical index in dynamics block

               t_tend(kk,i) = bbuffer(bpter(icol,k))
               u_tend(kk,i) = bbuffer(bpter(icol,k)+1)
               v_tend(kk,i) = bbuffer(bpter(icol,k)+2)

               do m = 1, pcnst
                  tracers(m,kk,i) = bbuffer(bpter(icol,k)+2+m)
               end do

            end do
         end do
      end do

      deallocate( bbuffer )
      deallocate( cbuffer )

   end if
   call t_stopf('pd_copy')

   call t_startf('derived_tend')
   call derived_tend(nCellsSolve, t_tend, u_tend, v_tend, dyn_in)
   call t_stopf('derived_tend')


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

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   character(len=*), parameter :: subname = 'dp_coupling::derived_phys'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(1, 'begin '//subname)

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
         ! include water vapor in in moist pdel.  
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

      ! convert the MPAS modified moist potential temperature to temperature
      do k = 1, pver
         phys_state(lchnk)%t(:ncol,k) = phys_state(lchnk)%t(:ncol,k) &
            / (1.0_r8 + (rh2o/rairv(:ncol,k,lchnk))*phys_state(lchnk)%q(:ncol,k,1)) &
            * (phys_state(lchnk)%pmiddry(:ncol,k)/pref)**(rairv(:ncol,k,lchnk)/cpairv(:ncol,k,lchnk))
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
#if FIX_TOTE
         ! general formula:  E = CV_air T + phis + gravit*zi )
         ! hydrostatic case: integrate zi term by parts, use CP=CV+R to get:
         ! E = CP_air T + phis   (Holton Section 8.3)
         ! to use this, update geopotential.F90, and other not-yet-found physics routines:
         ! (check boundary layer code, others which have gravit and zi() or zm()
         phys_state(lchnk)%s(:ncol,k) = cpairv(:ncol,k,lchnk)*phys_state(lchnk)%t(:ncol,k) &
            + phys_state(lchnk)%phis(:ncol)
#else
         phys_state(lchnk)%s(:ncol,k) = cpairv(:ncol,k,lchnk)*phys_state(lchnk)%t(:ncol,k) &
            + gravit*phys_state(lchnk)%zm(:ncol,k) + phys_state(lchnk)%phis(:ncol)
#endif
      end do

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)

   end do

end subroutine derived_phys

!=========================================================================================

subroutine derived_tend(nCellsSolve, t_tend, u_tend, v_tend, dyn_in)

   ! Derive the physics tendencies required by MPAS from the tendencies produced by
   ! CAM's physics package.

   ! Arguments
   integer,             intent(in)    :: nCellsSolve
   real(r8),            intent(in)    :: t_tend(pver,nCellsSolve) ! physics dtdt
   real(r8),            intent(in)    :: u_tend(pver,nCellsSolve) ! physics dudt
   real(r8),            intent(in)    :: v_tend(pver,nCellsSolve) ! physics dvdt
   type(dyn_import_t),  intent(inout) :: dyn_in

   ! Local variables


   ! variables from dynamics import container
   integer :: nEdgesSolve
   real(r8), pointer :: ru_tend(:,:)
   real(r8), pointer :: rtheta_tend(:,:)
   real(r8), pointer :: rho_tend(:,:)


   character(len=*), parameter :: subname = 'dp_coupling:derived_tend'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   nEdgesSolve = dyn_in % nEdgesSolve
   ru_tend     => dyn_in % ru_tend
   rtheta_tend => dyn_in % rtheta_tend
   rho_tend    => dyn_in % rho_tend

   ! insert code to derive MPAS dycore tendencies from CAM physics tendencies.
   ! set to zero for now...

   ru_tend     = 0._r8
   rtheta_tend = 0._r8
   rho_tend    = 0._r8

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

   ! Arguments
   integer, intent(in) :: nCells
   integer, intent(in) :: nVertLevels
   real(r8), dimension(nVertLevels, nCells), intent(in) :: zz         ! d(zeta)/dz [-]
   real(r8), dimension(nVertLevels+1, nCells), intent(in) :: zgrid    ! geometric heights of layer interfaces [m]
   real(r8), dimension(nVertLevels, nCells), intent(in) :: rho_zz     ! dry density / zz [kg m^-3]
   real(r8), dimension(nVertLevels, nCells), intent(in) :: theta_m    ! potential temperature * (1 + Rv/Rd * qv)
   real(r8), dimension(nVertLevels, nCells), intent(out) :: pmiddry   ! layer midpoint dry hydrostatic pressure [Pa]
   real(r8), dimension(nVertLevels+1, nCells), intent(out) :: pintdry ! layer interface dry hydrostatic pressure [Pa]

   ! Constants (which should probably come from elsewhere?)
   real(r8), parameter :: cp = 1004.5_r8
   real(r8), parameter :: rgas = 287.0_r8
   real(r8), parameter :: cv = 717.5_r8
   real(r8), parameter :: p0 = 1.0e5_r8
   real(r8), parameter :: g = 9.806_r8

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
   ptop_int(:) = 2.0 * ptop_mid(:) / (1.0 + exp( (zgrid(nVertLevels+1,:) - zgrid(nVertLevels,:)) * g / rgas / ttop_mid(:)))


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
         pintdry(k,iCell) = pintdry(k+1,iCell) + g * zz(k,iCell) * rho_zz(k,iCell) * dz(k)
         pmiddry(k,iCell) = 0.5 * (pintdry(k+1,iCell) + pintdry(k,iCell))
      end do
   end do

end subroutine dry_hydrostatic_pressure

!=========================================================================================

end module dp_coupling
