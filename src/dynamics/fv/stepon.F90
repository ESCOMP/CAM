module stepon

!----------------------------------------------------------------------
! stepon provides the interface layer that allows the different dynamical
! cores to be called from different locations in the time loop.  It also
! provides a standard interface that is called from the higher level CAM   
! component run methods while leaving non-standardized dycore interface
! methods to be called from this layer.  Ideally only the run methods
! which allow flexibility in the dynamics/physics calling sequence should
! remain.  The init and finalize methods should be removed and their
! functionality incorporated in the dycore init and finalize.
!----------------------------------------------------------------------

use shr_kind_mod,       only: r8 => shr_kind_r8

use spmd_utils,         only: mpicom, iam, masterproc
use cam_control_mod,    only: initial_run, moist_physics, simple_phys
use ppgrid,             only: begchunk, endchunk
use physconst,          only: zvir, cappa

use physics_types,      only: physics_state, physics_tend

use dyn_comp,           only: dyn_import_t, dyn_export_t, initial_mr
use dynamics_vars,      only: t_fvdycore_state, t_fvdycore_grid
use dyn_internal_state, only: get_dyn_state, get_dyn_state_grid

use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun
use perf_mod,           only: t_startf, t_stopf, t_barrierf

implicit none
private
save

public :: &
   stepon_init,  &! Initialization
   stepon_run1,  &! run method phase 1
   stepon_run2,  &! run method phase 2
   stepon_run3,  &! run method phase 3
   stepon_final   ! Finalization

integer :: pdt       ! Physics time step
real(r8) :: dtime    ! Physics time step
real(r8) :: te0            ! Total energy before dynamics

! for fv_out
logical, parameter :: fv_monitor=.true.  ! Monitor Mean/Max/Min fields
                                         ! This is CPU-time comsuming;
                                         ! set it to false for production runs
real (r8) :: ptop

!=========================================================================================
contains
!=========================================================================================

subroutine stepon_init(dyn_in, dyn_out)

   use constituents, only: pcnst
   use time_manager, only: get_step_size
   use physconst,    only: physconst_calc_kappav, rair, cpair
   use inic_analytic,      only: analytic_ic_active

   type (dyn_import_t)   :: dyn_in             ! Dynamics import container
   type (dyn_export_t)   :: dyn_out            ! Dynamics export container

   ! local variables:
   type (t_fvdycore_grid), pointer :: grid

   integer :: im, km
   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
   integer :: i,k,j,m             ! longitude, level, latitude and tracer indices
   logical :: nlres = .false.  ! true => restart or branch run

   integer :: ks
   real (r8), pointer :: ak(:)
   real (r8), pointer :: bk(:)

   real(r8), allocatable :: delpdryxy(:,:,:)
   real(r8), allocatable :: cap3vi(:,:,:), cappa3v(:,:,:)
   !----------------------------------------------------------------------------

   if (.not. initial_run) nlres=.true.

   grid => get_dyn_state_grid()
   im      =  grid%im
   km      =  grid%km


   ifirstxy  =  grid%ifirstxy
   ilastxy  =  grid%ilastxy
   jfirstxy  =  grid%jfirstxy
   jlastxy  =  grid%jlastxy

   ks     =  grid%ks
   ptop   =  grid%ptop
   ak     => grid%ak
   bk     => grid%bk

   pdt = get_step_size()    ! Physics time step
   dtime = pdt

   do j = jfirstxy, jlastxy
      do i=ifirstxy, ilastxy
         dyn_in%pe(i,1,j) = ptop
      enddo
   enddo

   if ( nlres) then ! restart or branch run
      !
      ! read_restart_dynamics delivers phis, ps, u3s, v3s, delp, pt
      ! in XY decomposition

      !
      ! Do not recalculate delta pressure (delp) if this is a restart run.
      ! Re. SJ Lin: The variable "delp" (pressure thikness for a Lagrangian
      ! layer) must be in the restart file. This is because delp will be
      ! modified "after" the physics update (to account for changes in water
      ! vapor), and it can not be reproduced by surface pressure and the
      ! ETA coordinate's a's and b's.

!$omp parallel do private(i,j,k)
      do j = jfirstxy, jlastxy
        do k=1, km
          do i=ifirstxy, ilastxy
            dyn_in%pe(i,k+1,j) = dyn_in%pe(i,k,j) + dyn_in%delp(i,j,k)
          enddo
        enddo
      enddo
   else
 
      ! Initial run --> generate pe and delp from the surface pressure
 
!$omp parallel do private(i,j,k)
         do j = jfirstxy, jlastxy
            do k=1,km+1
               do i=ifirstxy, ilastxy
                  dyn_in%pe(i,k,j) = ak(k) + bk(k) * dyn_in%ps(i,j)
               enddo
            enddo
         enddo

!$omp parallel do private(i,j,k)
         do k = 1, km
            do j = jfirstxy, jlastxy
               do i= ifirstxy, ilastxy
                  dyn_in%delp(i,j,k) = dyn_in%pe(i,k+1,j) - dyn_in%pe(i,k,j)
               enddo
            enddo
         enddo
   endif

   !----------------------------------------------------------
   ! Check total dry air mass; set to 982.22 mb if initial run
   ! Print out diagnostic message if restart run
   !----------------------------------------------------------

   if (.not. simple_phys) then
      call dryairm( grid, .true., dyn_in%ps, dyn_in%tracer,  &
                    dyn_in%delp, dyn_in%pe, nlres )
   endif

   if (grid%iam < grid%npes_xy) then

      allocate( cappa3v(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( cap3vi(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      if (grid%high_alt) then
         call physconst_calc_kappav( ifirstxy,ilastxy,jfirstxy,jlastxy,1,km, grid%ntotq, dyn_in%tracer, cappa3v )

!$omp parallel do private(i,j,k)
         do k=2,km
            do j=jfirstxy,jlastxy
               do i=ifirstxy,ilastxy
                  cap3vi(i,j,k) = 0.5_r8*(cappa3v(i,j,k-1)+cappa3v(i,j,k))
               enddo
            enddo
         enddo
         cap3vi(:,:,1) = 1.5_r8 * cappa3v(:,:,1) - 0.5_r8 * cappa3v(:,:,2)
         cap3vi(:,:,km+1) = 1.5_r8 * cappa3v(:,:,km) - 0.5_r8 * cappa3v(:,:,km-1)
      else
         cappa3v = rair/cpair
         cap3vi = rair/cpair
      endif

      ! Initialize pk, edge pressure to the cappa power.  Do this with constituent dependent cappa

!$omp parallel do private(i,j,k)
      do k = 1, km+1
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pk(i,j,k) = dyn_in%pe(i,k,j)**cap3vi(i,j,k)
            enddo
         enddo
      enddo

      ! Generate pkz, the conversion factor betw pt and t3

      call pkez(1,      im,   km,       jfirstxy,  jlastxy,              &
                1,      km,   ifirstxy, ilastxy,    dyn_in%pe,    &
                dyn_in%pk, cappa3v,  ks, dyn_out%peln, dyn_out%pkz,  .false., grid%high_alt )

      deallocate( cappa3v, cap3vi )

   endif

   if (initial_run) then

      ! Compute pt for initial run: scaled virtual potential temperature
      ! defined as (virtual temp deg K)/pkz. pt will be written to restart (SJL)

!$omp parallel do private(i,j,k)
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pt(i,j,k) =  dyn_in%t3(i,j,k)*            &
                (1._r8 + zvir*dyn_in%tracer(i,j,k,1))    &
                /dyn_in%pkz(i,j,k) 
            enddo
         enddo
      enddo

      !----------------------------------------------------------------
      ! Convert mixing ratios initialized as dry to moist for dynamics
      !----------------------------------------------------------------

      ! on initial time step, dry mixing ratio advected constituents have been
      !    initialized to dry mixing ratios. dynpkg expects moist m.r. so convert here.

      ! first calculate delpdry. The set_pdel_state subroutine
      !   is called after the dynamics in d_p_coupling to set more variables.
      !   This is not in tracers.F90 because it is only used by LR dynamics.
      allocate (delpdryxy(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km))
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               delpdryxy(i,j,k) = dyn_in%delp(i,j,k)*          &
                    (1._r8 - dyn_in%tracer(i,j,k,1))
            enddo
         enddo
      enddo
      do m = 1,pcnst
         if (initial_mr(m) == 'dry') then
            do k=1, km
               do j = jfirstxy, jlastxy
                  do i = ifirstxy, ilastxy
                     dyn_in%tracer(i,j,k,m) =               &
                             dyn_in%tracer(i,j,k,m)*        &
                             delpdryxy(i,j,k)/dyn_in%delp(i,j,k)
                  end do
               end do
            end do
         end if
      end do
      deallocate (delpdryxy)
      
   end if

end subroutine stepon_init

!=========================================================================================

subroutine stepon_run1( dtime_out, phys_state, phys_tend, pbuf2d,        &
                        dyn_in, dyn_out )

   ! Phase 1 run of FV dynamics. Run the dynamics, and couple to physics.

   use dp_coupling,       only: d_p_coupling
   use dyn_comp,          only: dyn_run
   
   use physics_buffer,    only: physics_buffer_desc
   use advect_tend,       only: compute_adv_tends_xyz

   ! arguments
   real(r8),            intent(out)   :: dtime_out   ! Time-step
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(dyn_import_t)                 :: dyn_in  ! Dynamics import container
   type(dyn_export_t)                 :: dyn_out ! Dynamics export container

   type(T_FVDYCORE_STATE), pointer :: dyn_state

   integer  :: rc 

   dtime_out = dtime
   dyn_state => get_dyn_state()

   ! Dump state variables to IC file
   call t_barrierf('sync_diag_dynvar_ic', mpicom)
   call t_startf ('diag_dynvar_ic')
   call diag_dynvar_ic (dyn_state%grid, dyn_out%phis, dyn_out%ps,             &
                        dyn_out%t3, dyn_out%u3s, dyn_out%v3s, dyn_out%tracer  )
   call t_stopf  ('diag_dynvar_ic')

   call t_startf ('comp_adv_tends1')
   call compute_adv_tends_xyz(dyn_state%grid, dyn_in%tracer )
   call t_stopf  ('comp_adv_tends1')
   !
   !--------------------------------------------------------------------------
   ! Perform finite-volume dynamics -- this dynamical core contains some 
   ! yet to be published algorithms. Its use in the CAM is
   ! for software development purposes only. 
   ! Please contact S.-J. Lin (Shian-Jiann.Lin@noaa.gov)
   ! if you plan to use this mudule for scientific purposes. Contact S.-J. Lin
   ! or Will Sawyer (sawyer@gmao.gsfc.nasa.gov) if you plan to modify the
   ! software.
   !--------------------------------------------------------------------------

   !----------------------------------------------------------
   ! For 2-D decomposition, phisxy is input to dynpkg, and the other
   ! xy variables are output. Some are computed through direct
   ! transposes, and others are derived.
   !----------------------------------------------------------
   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(ptop,      pdt,     te0,         &
                dyn_state, dyn_in,  dyn_out,  rc )
   if ( rc /= 0 ) then
     write(iulog,*) "STEPON_RUN: dyn_run returned bad error code", rc
     write(iulog,*) "Quitting."
     call endrun
   endif 
   call t_stopf  ('dyn_run')

   call t_startf ('comp_adv_tends2')
   call compute_adv_tends_xyz(dyn_state%grid, dyn_out%tracer )
   call t_stopf  ('comp_adv_tends2')

   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling(dyn_state%grid, phys_state, phys_tend,  pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')

!EOC
end subroutine stepon_run1

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run2 -- second phase run method
!
! !INTERFACE:
subroutine stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
! !USES:
   use dp_coupling,      only: p_d_coupling
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

   type (T_FVDYCORE_GRID), pointer :: grid
!
! !DESCRIPTION:
!
! Second phase run method. Couple from physics to dynamics.
!
!EOP
!-----------------------------------------------------------------------
!BOC

!-----------------------------------------------------------------------

   !----------------------------------------------------------
   ! Update dynamics variables using phys_state & phys_tend.
   ! 2-D decomposition: Compute ptxy and q3xy; for ideal
   !   physics, scale ptxy by (old) pkzxy; then transpose to yz variables
   ! 1-D decomposition: Compute dudt, dvdt, pt and q3; for ideal physics,
   !   scale pt by old pkz.
   ! Call uv3s_update to update u3s and v3s from dudt and dvdt.
   ! Call p_d_adjust to update pt, q3, pe, delp, ps, piln, pkz and pk.
   ! For adiabatic case, transpose to yz variables.
   !----------------------------------------------------------
   grid => get_dyn_state_grid()

   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf ('p_d_coupling')
   call p_d_coupling(grid, phys_state, phys_tend, &
                     dyn_in, dtime, zvir, cappa, ptop)
   call t_stopf  ('p_d_coupling')

!EOC
end subroutine stepon_run2

!-----------------------------------------------------------------------

subroutine stepon_run3(dtime, cam_out, phys_state,             &
                        dyn_in, dyn_out )
! !USES:
   use time_manager,     only: get_curr_date
   use fv_prints,        only: fv_out
   use camsrfexch,       only: cam_out_t    
!
! !INPUT PARAMETERS:
!
   type(physics_state), intent(in):: phys_state(begchunk:endchunk)
   real(r8), intent(in) :: dtime            ! Time-step
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
!
! !INPUT/OUTPUT PARAMETERS:
!
   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
!
! !DESCRIPTION:
!
!	Final run phase of dynamics. Some printout and time index updates.
!
! !HISTORY:
!   2005.09.16  Kluzek     Creation
!   2006.04.13  Sawyer     Removed shift_time_indices (not needed in FV)
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

   type(t_fvdycore_state), pointer :: state
   type(t_fvdycore_grid),  pointer :: grid
   integer :: ncdate            ! current date in integer format [yyyymmdd]
   integer :: ncsec             ! time of day relative to current date [seconds]
   integer :: yr, mon, day      ! year, month, day components of a date
   integer :: ncsecp
   integer :: freq_diag

   !----------------------------------------------------------
   ! Monitor max/min/mean of selected fields
   !
   !  SEE BELOW  ****  SEE BELOW  ****  SEE BELOW
   
   ! Beware that fv_out uses both dynamics and physics instantiations.
   ! However, I think that they are used independently, so that the
   ! answers are correct. Still, this violates the notion that the
   ! physics state is no longer active after p_d_coupling.
   !----------------------------------------------------------
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
   ncsecp = ncsec + pdt      !  step complete, but nstep not incremented yet

   state => get_dyn_state()
   freq_diag = state%check_dt

   if (fv_monitor .and. mod(ncsecp, freq_diag) == 0) then
      grid => state%grid

      call t_barrierf('sync_fv_out', mpicom)
      call t_startf('fv_out')
      call fv_out(grid, dyn_out%pk, dyn_out%pt,         &
                  ptop, dyn_out%ps, dyn_out%tracer,     &
                  dyn_out%delp, dyn_out%pe, cam_out,    &
                   phys_state, ncdate, ncsecp, moist_physics)
      call t_stopf('fv_out')
   endif

!EOC
end subroutine stepon_run3

!-----------------------------------------------------------------------

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)

! !PARAMETERS:
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC

!!! Not yet ready for the call to dyn_final
!!! call dyn_final( RESTART_FILE, dyn_state, dyn_in, dyn_out )
!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------

end module stepon
