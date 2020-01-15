#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module dp_coupling

!-------------------------------------------------------------------------------
! dynamics - physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use pmgrid,         only: plev
use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp
use constituents,   only: pcnst, cnst_name

use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
use spmd_utils,     only: mpicom, iam

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

!  use mpas_cam_interface,  only: mpas_to_cam, cam_to_mpas

implicit none
private
save

public :: &
   d_p_coupling, &
   p_d_coupling


integer, parameter :: nblocks_per_pe = 1

!=========================================================================================
contains
!=========================================================================================


!-----------------------------------------------------------------------
!  routine d_p_coupling
!
!> \brief Dynamics to physics coupling
!> \details
!>  Handles coupling of current dynamics state to the physics by filling
!>  in arrays in the phys_state based on the contents of the dynamics export
!>  state in dyn_out.
!
!-----------------------------------------------------------------------
subroutine d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)

   use dyn_grid, only : nCellsSolve, nVertLevelsSolve

   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_export_t),        intent(inout) :: dyn_out
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   ! LOCAL VARIABLES

   ! Variables from dynamics export container
   real(r8), pointer :: pmid(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: fzm(:)
   real(r8), pointer :: fzp(:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8), pointer :: w(:,:)
   real(r8), pointer :: theta_m(:,:)
   real(r8), pointer :: tracers(:,:,:)

   integer :: lchnk, icol, k      ! indices over chunks, columns, layers
   integer :: ncols, ig, nblk, nb, m, i
   integer :: index_qv

   integer :: pgcols(pcols)
   integer :: tsize                 ! amount of data per grid point passed to physics
   integer :: bpter(max_col_per_block,0:pver)    ! offsets into block buffer for packing data
   integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   character(len=*), parameter :: subname = 'dp_coupling::d_p_coupling'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   pmid => dyn_out % pmid
   zint => dyn_out % zint
   zz => dyn_out % zz
   fzm => dyn_out % fzm
   fzp => dyn_out % fzp
   rho_zz => dyn_out % rho_zz
   ux => dyn_out % ux
   uy => dyn_out % uy
   w => dyn_out % w
   theta_m => dyn_out % theta_m
   tracers => dyn_out % tracers
   index_qv = dyn_out % index_qv

   call t_startf('dpcopy')
   if (local_dp_map) then

#ifdef USE_TOTALLY_UNTESTED_MPAS_CODE
      !$omp parallel do private (lchnk, ncols, icol, i, k, m, pgcols)
      do lchnk = begchunk, endchunk

         ncols = get_ncols_p(lchnk)                         ! number of columns in this chunk
         call get_gcol_all_p(lchnk, pcols, pgcols)          ! global column indices

         do icol = 1, ncols
            i = global_to_local_cell(pgcols(icol))        ! local (to process) column index

            phys_state(lchnk)%psdry(icol) = rho_zz(1,i) * zz(1,i) * 9.806 * 0.5 * (zint(2,i) - zint(1,i)) + pmid(1,i)   ! psd
            phys_state(lchnk)%phis(icol) = zint(1,i) * 9.806   ! phis

            do k = 1, pver
               phys_state(lchnk)%t(icol,k) = theta_m(k,i) / (1.0_r8 + 1.61_r8 * tracers(index_qv,k,i)) * (pmid(k,i) / 1.0e5)**(287.0_r8 / 1003.0_r8)   ! t
               phys_state(lchnk)%u(icol,k) = ux(k,i)
               phys_state(lchnk)%v(icol,k) = uy(k,i)
               phys_state(lchnk)%omega(icol,k) = -rho_zz(k,i) * zz(k,i) * 9.806 * 0.5 * (w(k,i) + w(k+1,i))   ! omega
               phys_state(lchnk)%pmiddry(icol,k)=pmid(k,i)   ! Wrong: these are full pressures at present
               phys_state(lchnk)%pmid(icol,k) = pmid(k,i)
               phys_state(lchnk)%zi(icol,k) = zint(k,i)
               phys_state(lchnk)%zm(icol,k) = 0.5_r8 * (zint(k,i) + zint(k+1,i))   ! zmid
            end do
            phys_state(lchnk)%zi(icol,pverp) = zint(pverp,i)

            phys_state(lchnk)%pint(icol,1) = rho_zz(1,i) * zz(1,i) * 9.806 * 0.5 * (zint(2,i) - zint(1,i)) + pmid(1,i)
            phys_state(lchnk)%pintdry(icol,1) = rho_zz(1,i) * zz(1,i) * 9.806 * 0.5 * (zint(2,i) - zint(1,i)) + pmid(1,i)
            do k=2,pver
               phys_state(lchnk)%pint(icol,k) = fzp(k) * pmid(k-1,i) + fzm(k) * pmid(k,i)
               phys_state(lchnk)%pintdry(icol,k) = fzp(k) * pmid(k-1,i) + fzm(k) * pmid(k,i)   ! Wrong: this is still full pressure
            end do
            phys_state(lchnk)%pint(icol,pverp) = rho_zz(pver,i) * zz(pver,i) * 9.806 * 0.5 * (zint(pver,i) - zint(pver+1,i)) + pmid(pver,i)
            phys_state(lchnk)%pintdry(icol,pverp) = rho_zz(pver,i) * zz(pver,i) * 9.806 * 0.5 * (zint(pver,i) - zint(pver+1,i)) + pmid(pver,i)

            do m=1,pcnst
               do k=1,pver
                  phys_state(lchnk)%q(icol,k,m) = tracers(m,k,i)
               end do
            end do
         end do
      end do
#endif


   else  ! .not. local_dp_map

      tsize = 8 + pcnst
      allocate(bbuffer(tsize*block_buf_nrecs))    ! block buffer
      bbuffer = 0.0_r8
      allocate(cbuffer(tsize*chunk_buf_nrecs))    ! chunk buffer
      cbuffer = 0.0_r8

      !$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, m, bpter)
      do nb = 1, nblocks_per_pe
         nblk = iam * nblocks_per_pe + nb   ! global block index
         ncols = get_block_gcol_cnt_d(nblk) ! number of columns in this block

         call block_to_chunk_send_pters(nblk,max_col_per_block,pver+1,tsize,bpter)

         do icol=1,ncols
            ig = col_indices_in_block(icol,nblk)   !  global column index
            i = global_to_local_cell(ig)           !  local (to process) column index

            bbuffer(bpter(icol,0))   = rho_zz(1,i) * zz(1,i) * 9.806 * 0.5 * (zint(2,i) - zint(1,i)) + pmid(1,i)   ! psd
            bbuffer(bpter(icol,0)+1) = zint(1,i) * 9.806   ! phis

            do k=1,pver
               bbuffer(bpter(icol,k))   = theta_m(k,i) / (1.0_r8 + 1.61_r8 * tracers(index_qv,k,i)) * (pmid(k,i) / 1.0e5)**(287.0_r8 / 1003.0_r8)   ! t
               bbuffer(bpter(icol,k)+1) = ux(k,i)
               bbuffer(bpter(icol,k)+2) = uy(k,i)
               bbuffer(bpter(icol,k)+3) = -rho_zz(k,i) * zz(k,i) * 9.806 * 0.5 * (w(k,i) + w(k+1,i))   ! omega
               bbuffer(bpter(icol,k)+4) = pmid(k,i)    ! Full pressure, not dry pressure
               bbuffer(bpter(icol,k)+5) = 0.5_r8 * (zint(k,i) + zint(k+1,i))   ! zmid

               do m=1,pcnst
                  bbuffer(bpter(icol,k)+5+m) = tracers(m,k,i)
               end do
            end do

            do k=0,pver
               bbuffer(bpter(icol,k)+7+pcnst) = zint(k+1,i)
            end do

            bbuffer(bpter(icol,0)+6+pcnst) = rho_zz(1,i) * zz(1,i) * 9.806 * 0.5 * (zint(2,i) - zint(1,i)) + pmid(1,i)
            do k=2,pver
               bbuffer(bpter(icol,k-1)+6+pcnst) = fzp(k) * pmid(k-1,i) + fzm(k) * pmid(k,i)    ! Full pressure, not dry pressure
            end do
            bbuffer(bpter(icol,pver)+6+pcnst) = rho_zz(pver,i) * zz(pver,i) * 9.806 * 0.5 * (zint(pver,i) - zint(pver+1,i)) + pmid(pver,i)

         end do

      end do

      call t_barrierf ('sync_blk_to_chk', mpicom)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')

      !$omp parallel do private (lchnk, ncols, icol, k, m, cpter)
      do lchnk = begchunk,endchunk
         ncols = phys_state(lchnk)%ncol

         call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

         do icol = 1, ncols

            phys_state(lchnk)%psdry(icol)    = cbuffer(cpter(icol,0))
            phys_state(lchnk)%phis(icol)     = cbuffer(cpter(icol,0)+1)

            do k=1,pver

               phys_state(lchnk)%t      (icol,k)  = cbuffer(cpter(icol,k))
               phys_state(lchnk)%u      (icol,k)  = cbuffer(cpter(icol,k)+1)
               phys_state(lchnk)%v      (icol,k)  = cbuffer(cpter(icol,k)+2)
               phys_state(lchnk)%omega  (icol,k)  = cbuffer(cpter(icol,k)+3)
               phys_state(lchnk)%pmiddry(icol,k)  = cbuffer(cpter(icol,k)+4)
               phys_state(lchnk)%zm     (icol,k)  = cbuffer(cpter(icol,k)+5)

               do m=1,pcnst
                  phys_state(lchnk)%q  (icol,k,m) = cbuffer(cpter(icol,k)+5+m)
               end do

            end do

            do k=0,pver
               phys_state(lchnk)%pintdry(icol,k+1) = cbuffer(cpter(icol,k)+6+pcnst)
               phys_state(lchnk)%zi     (icol,k+1) = cbuffer(cpter(icol,k)+7+pcnst)
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


!-----------------------------------------------------------------------
!  routine p_d_coupling
!
!> \brief Physics to dynamics coupling
!> \details
!>  Handles coupling of physics to the dynamics by filling in arrays in
!>  the dynamics import state based on the contents of updated state and
!>  tendencies from the physics.
!
!-----------------------------------------------------------------------
subroutine p_d_coupling(phys_state, phys_tend, dyn_in)

   ! arguments
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ), intent(inout) :: phys_tend(begchunk:endchunk)
   type(dyn_import_t),  intent(inout) :: dyn_in

   ! LOCAL VARIABLES
   integer :: ic , ncols                            ! index
   integer :: lchnk, icol, k      ! indices over chunks, columns, layers

   integer :: m, i, j, nb, nblk, ig
   integer :: pgcols(pcols)

   integer :: tsize                 ! amount of data per grid point passed to physics
   integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data
   integer :: bpter(max_col_per_block,0:pver)    ! offsets into block buffer for unpacking data

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers

   ! Variables from dynamics import container
   real(r8), pointer :: pmid(:,:)
   real(r8), pointer :: zint(:,:)
   real(r8), pointer :: zz(:,:)
   real(r8), pointer :: rho_zz(:,:)
   real(r8), pointer :: ux(:,:)
   real(r8), pointer :: uy(:,:)
   real(r8), pointer :: theta(:,:)
   real(r8), pointer :: w(:,:)
   real(r8), pointer :: tracers(:,:,:)
   real(r8), pointer :: ru_tend(:,:)
   real(r8), pointer :: rtheta_tend(:,:)
   real(r8), pointer :: rho_tend(:,:)

   character(len=*), parameter :: subname = 'dp_coupling::p_d_coupling'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   pmid => dyn_in % pmid
   zint => dyn_in % zint
   zz => dyn_in % zz
   rho_zz => dyn_in % rho_zz
   ux => dyn_in % ux
   uy => dyn_in % uy
   theta => dyn_in % theta
   w => dyn_in % w
   tracers => dyn_in % tracers
   ru_tend => dyn_in % ru_tend
   rtheta_tend => dyn_in % rtheta_tend
   rho_tend => dyn_in % rho_tend

   call t_startf('pd_copy')
   if (local_dp_map) then

      !$omp parallel do private (lchnk, ncols, icol, i, k, m, pgcols)
      do lchnk=begchunk,endchunk
         ncols=get_ncols_p(lchnk)                         ! number of columns in this chunk
         call get_gcol_all_p(lchnk,pcols,pgcols)          ! global column indices
         do icol=1,ncols
            i = global_to_local_cell(pgcols(icol))        ! local (to process) column index

            do k=1,pver
#ifdef USE_UNCOMPILABLE_MPAS_CODE
               temp_tend(k,i) = phys_tend(lchnk)%dtdt(icol,k)
               ux_tend(k,i)   = phys_tend(lchnk)%dudt(icol,k)
               uy_tend(k,i)   = phys_tend(lchnk)%dvdt(icol,k)
#endif
               do m=1,pcnst
                  tracers(m,k,i) = phys_state(lchnk)%q(icol,k,m)
               end do
            end do

         end do

      end do

   else

      tsize = 3 + pcnst

      allocate( bbuffer(tsize*block_buf_nrecs) )
      allocate( cbuffer(tsize*chunk_buf_nrecs) )

      !$omp parallel do private (lchnk, ncols, icol, k, m, cpter)
      do lchnk = begchunk,endchunk
         ncols = get_ncols_p(lchnk)

         call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

         do icol=1,ncols
            cbuffer(cpter(icol,0):cpter(icol,0)+2+pcnst) = 0.0_r8
         end do

         do icol=1,ncols

            do k=1,pver
               cbuffer   (cpter(icol,k))     = phys_tend(lchnk)%dtdt(icol,k)
               cbuffer   (cpter(icol,k)+1)   = phys_tend(lchnk)%dudt(icol,k)
               cbuffer   (cpter(icol,k)+2)   = phys_tend(lchnk)%dvdt(icol,k)

               do m=1,pcnst
                  cbuffer(cpter(icol,k)+2+m) = phys_state(lchnk)%q(icol,k,m)
               end do
            end do

         end do

      end do

      call t_barrierf('sync_chk_to_blk', mpicom)
      call t_startf ('chunk_to_block')
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      call t_stopf  ('chunk_to_block')

      !$omp parallel do private (nb, nblk, ncols, icol, ig, i, k, m, bpter)
      do nb = 1, nblocks_per_pe
         nblk = iam * nblocks_per_pe + nb   !  global block index
         ncols = get_block_gcol_cnt_d(nblk) !  number of columns in this block

         call chunk_to_block_recv_pters(nblk,max_col_per_block,pver+1,tsize,bpter)

         do icol=1,ncols
            ig = col_indices_in_block(icol,nblk)   !  global column index
            i = global_to_local_cell(ig)           !  local (to process) column index

            do k=1,pver

#ifdef USE_UNCOMPILABLE_MPAS_CODE
               temp_tend(k,i) = bbuffer(bpter(icol,k))
               ux_tend  (k,i) = bbuffer(bpter(icol,k)+1)
               uy_tend  (k,i) = bbuffer(bpter(icol,k)+2)
#endif

               do m=1,pcnst
                  tracers(m,k,i) = bbuffer(bpter(icol,k)+2+m)
               end do

            end do
            
         end do

      end do

      deallocate( bbuffer )
      deallocate( cbuffer )

   end if

   call t_stopf('pd_copy')

end subroutine p_d_coupling


!-----------------------------------------------------------------------
!  routine derived_phys
!
!> \brief Derives physics fields
!> \details
!>  Derives fields needed by the physics after some fields have been
!>  set by d_p_coupling and before the call to physics.
!
!-----------------------------------------------------------------------
subroutine derived_phys(phys_state, phys_tend, pbuf2d)

   use constituents,  only: qmin
   use physconst,     only: cpair, gravit, rair, zvir, cappa
   use spmd_utils,    only: masterproc
   use ppgrid,        only: pver
   use physics_types, only: set_state_pdry, set_wet_to_dry
   use check_energy,  only: check_energy_timestep_init
   use shr_vmath_mod, only: shr_vmath_log
   use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk

   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)

   integer :: lchnk
   real(r8) :: qbot                 ! bottom level q before change
   real(r8) :: qbotm1               ! bottom-1 level q before change
   real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
   real(r8) :: qmavl                ! available q at level pver-1

   real(r8) :: ke(pcols,begchunk:endchunk)   
   real(r8) :: se(pcols,begchunk:endchunk)   
   real(r8) :: ke_glob(1),se_glob(1)

   integer :: m, i, k, ncol
   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   real(r8):: p00

   character(len=*), parameter :: subname = 'dp_coupling::derived_phys'


   MPAS_DEBUG_WRITE(1, 'begin '//subname)

   p00   = 1.e5_r8

   ! Evaluate derived quantities

   !$omp parallel do private (lchnk, ncol, k, i)
   do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)

!       do k=1,pver
!          call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
!          call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
!       end do
!       call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)
!
!       do k=1,pver
!          do i=1,ncol
!             phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
!             phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
!             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) &
!                                             / phys_state(lchnk)%pmid(i,k))**cappa
!          end do
!       end do

      do k=1,pver
         do i=1,ncol
            phys_state(lchnk)%pdeldry(i,k) = phys_state(lchnk)%pintdry(i,k+1) - phys_state(lchnk)%pintdry(i,k)
         end do
      end do
      phys_state(lchnk)%rpdeldry(:ncol,:) = 1._r8/phys_state(lchnk)%pdeldry(:ncol,:)
      phys_state(lchnk)%lnpmiddry(:ncol,:) = log(phys_state(lchnk)%pmiddry(:ncol,:))
      phys_state(lchnk)%lnpintdry(:ncol,:) = log(phys_state(lchnk)%pintdry(:ncol,:))

      do k=1,pver
         phys_state(lchnk)%pdel(:ncol,k) = phys_state(lchnk)%pdeldry(:ncol,k)/(1._r8-phys_state(lchnk)%q(:ncol,k,1))

         !SHP-PS
         !phys_state(lchnk)%pdel(:ncol,k) = phys_state(lchnk)%pdeldry(:ncol,k)*(1._r8 + &
         !                                  phys_state(lchnk)%q(:ncol,k,1)/(rair/rh2o*(1._r8-phys_state(lchnk)%q(:ncol,k,1))))
      end do
       
      phys_state(lchnk)%ps(:ncol) = phys_state(lchnk)%pintdry(:ncol,1)
      phys_state(lchnk)%pint(:ncol,1) = phys_state(lchnk)%pintdry(:ncol,1)
      do k=1,pver
         phys_state(lchnk)%pint(:ncol,k+1) = phys_state(lchnk)%pint(:ncol,k)+phys_state(lchnk)%pdel(:ncol,k)
         phys_state(lchnk)%pmid(:ncol,k) = (phys_state(lchnk)%pint(:ncol,k+1)+phys_state(lchnk)%pint(:ncol,k))/2._r8
         phys_state(lchnk)%ps(:ncol) = phys_state(lchnk)%ps(:ncol) + phys_state(lchnk)%pdel(:ncol,k)
      end do

      do k=1,pver
         call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,k),phys_state(lchnk)%lnpint(1:ncol,k),ncol)
         call shr_vmath_log(phys_state(lchnk)%pmid(1:ncol,k),phys_state(lchnk)%lnpmid(1:ncol,k),ncol)
      end do
      call shr_vmath_log(phys_state(lchnk)%pint(1:ncol,pverp),phys_state(lchnk)%lnpint(1:ncol,pverp),ncol)

      do k=1,pver
         do i=1,ncol
            phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
            phys_state(lchnk)%rpdel(i,k) = 1._r8/phys_state(lchnk)%pdel(i,k)
!             phys_state(lchnk)%exner (i,k) = (phys_state(lchnk)%pint(i,pver+1) / phys_state(lchnk)%pmid(i,k))**cappa 
            phys_state(lchnk)%exner (i,k) = (p00 / phys_state(lchnk)%pmid(i,k))**cappa
         end do
      end do


      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
         do i=1,ncol
#if FIX_TOTE
            ! general formula:  E = CV_air T + phis + gravit*zi )
            ! hydrostatic case: integrate zi term by parts, use CP=CV+R to get:
            ! E = CP_air T + phis   (Holton Section 8.3)
            ! to use this, update geopotential.F90, and other not-yet-found physics routines:
            ! (check boundary layer code, others which have gravit and zi() or zm()
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                                      + phys_state(lchnk)%phis(i)
#else
            phys_state(lchnk)%s(i,k) = cpair*phys_state(lchnk)%t(i,k) &
                  + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
#endif
         end do
      end do

      ! Compute energy and water integrals of input state
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), phys_tend(lchnk), pbuf_chnk)


#if 0
      ke(:,lchnk) = 0._r8
      se(:,lchnk) = 0._r8
      do k = 1, pver
         do i = 1, ncol
            ke(i,lchnk) = ke(i,lchnk) + ( 0.5_r8*(phys_state(lchnk)%u(i,k)**2 + phys_state(lchnk)%v(i,k)**2)*phys_state(lchnk)%pdel(i,k) )/gravit
            se(i,lchnk) = se(i,lchnk) + phys_state(lchnk)%s(i,k         )*phys_state(lchnk)%pdel(i,k)/gravit
         end do
      end do
#endif

   end do

end subroutine derived_phys

end module dp_coupling
