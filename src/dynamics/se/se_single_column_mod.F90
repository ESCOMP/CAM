module se_single_column_mod
!--------------------------------------------------------
! 
! Module for the SE single column model

use shr_kind_mod,           only: r8=>shr_kind_r8
use element_mod,            only: element_t
use scamMod,                only: have_t, have_q, have_u, have_v, have_ps, have_numliq, &
                                  have_cldliq, have_numice, have_cldice, have_omega, use_camiop, &
                                  tobs, qobs,have_numliq, numliqobs, cldliqobs, numiceobs, cldiceobs, &
                                  wfld, psobs,uobs,vobs,tobs,divt,divQ,divT3d,divq3d,precobs,lhflxobs, &
                                  shflxobs, tground, have_ps, have_tg, have_lhflx, have_shflx, have_t, &
                                  have_omega, have_cldliq, have_divt, have_divq, have_divt3d, have_divq3d, &
                                  use_3dfrc
use constituents,           only: cnst_get_ind, pcnst
use dimensions_mod,         only: nelemd, np, nlev
use time_manager,           only: get_nstep, is_first_step, get_step_size, is_first_restart_step
use ppgrid,                 only: begchunk
use time_mod,               only: timelevel_qdp
use cam_history,            only: outfld

implicit none

private
save

public scm_setinitial
public scm_setfield
public apply_SC_forcing
public iop_broadcast

integer                      :: tl_f, tl_fqdp

!=========================================================================
contains
!=========================================================================

subroutine scm_setinitial(elem)

  use constituents, only: qmin
  use dyn_grid,     only: TimeLevel
  use control_mod,  only: qsplit

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, cix, ie, thelev
  integer inumliq, inumice, icldliq, icldice
  integer              :: tl_f, tl_fqdp

  tl_f = timelevel%n0
  call TimeLevel_Qdp(timelevel, qsplit, tl_fqdp)

  if (.not. use_camiop .and. get_nstep() .eq. 0) then
    call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
    call cnst_get_ind('NUMICE', inumice, abort=.false.)
    call cnst_get_ind('CLDLIQ', icldliq)
    call cnst_get_ind('CLDICE', icldice)

    do ie=1,nelemd
      do j=1,np
        do i=1,np

          ! Find level where tobs is no longer zero
          thelev=1
          do k=1, NLEV
            if (tobs(k) .ne. 0) then
              thelev=k
              go to 1000
            endif
          enddo

1000 continue
           
          if (get_nstep() .le. 1) then
            do k=1,thelev-1
              tobs(k)=elem(ie)%state%T(i,j,k,tl_f)
              qobs(k)=elem(ie)%state%qdp(i,j,k,1,tl_fqdp)/elem(ie)%state%dp3d(i,j,k,tl_f)
            enddo
          else
            tobs(:)=elem(ie)%state%T(i,j,:,tl_f)
            qobs(:)=elem(ie)%state%qdp(i,j,:,1,tl_fqdp)/elem(ie)%state%dp3d(i,j,:,tl_f)
          endif

          if (get_nstep() .eq. 0) then
            do cix = 1, pcnst
!jt               if (scm_zero_non_iop_tracers) elem(ie)%state%qdp(i,j,:,cix,tl_qdp_np0) = qmin(cix)*elem(ie)%state%dp3d(i,j,:,tl_qdp_np0)
               elem(ie)%state%qdp(i,j,:,cix,tl_fqdp) = qmin(cix)*elem(ie)%state%dp3d(i,j,:,tl_f)
            end do
            do k=thelev, NLEV
              if (have_t) elem(ie)%state%T(i,j,k,tl_f)=tobs(k)
              if (have_q) elem(ie)%state%qdp(i,j,k,1,tl_fqdp)=qobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)
!jt              if (have_q) elem(ie)%state%qdp(i,j,k,1,tl_f)=qobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)

           enddo

            do k=1,NLEV
              if (have_ps) elem(ie)%state%psdry(i,j) = psobs
              if (have_u) elem(ie)%state%v(i,j,1,k,tl_f) = uobs(k)
              if (have_v) elem(ie)%state%v(i,j,2,k,tl_f) = vobs(k)
              if (have_numliq) elem(ie)%state%qdp(i,j,k,inumliq,tl_fqdp) = numliqobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)
              if (have_cldliq) elem(ie)%state%qdp(i,j,k,icldliq,tl_fqdp) = cldliqobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)
              if (have_numice) elem(ie)%state%qdp(i,j,k,inumice,tl_fqdp) = numiceobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)
              if (have_cldice) elem(ie)%state%qdp(i,j,k,icldice,tl_fqdp) = cldiceobs(k)*elem(ie)%state%dp3d(i,j,k,tl_f)
              if (have_omega) elem(ie)%derived%omega(i,j,k) = wfld(k)
            enddo

          endif

        enddo
      enddo
    enddo
  endif

end subroutine scm_setinitial

subroutine scm_setfield(elem,iop_update_phase1)

!---------------------------------------------------------
! Purpose: Update various fields based on available data
!   provided by IOP file
!----------------------------------------------------------

  implicit none

  logical, intent(in) :: iop_update_phase1
  type(element_t), intent(inout) :: elem(:)

  integer i, j, k, ie

  do ie=1,nelemd
    if (have_ps .and. use_camiop .and. .not. iop_update_phase1) elem(ie)%state%psdry(:,:) = psobs
    if (have_ps .and. .not. use_camiop) elem(ie)%state%psdry(:,:) = psobs
    do i=1, NLEV
      if (have_omega .and. iop_update_phase1) elem(ie)%derived%omega(:,:,i)=wfld(i)  !     set t to tobs at first
    end do
  end do

end subroutine scm_setfield

subroutine apply_SC_forcing(elem,hvcoord,tl,n,t_before_advance,nets,nete)
! 
    use scamMod,        only: single_column, use_3dfrc
    use dimensions_mod, only : np, nlev, npsq,qsize_d
    
    use hybvcoord_mod,  only : hvcoord_t
    use element_mod,    only : element_t
    use physconst,      only: rair
    use time_mod
    use time_manager,   only: get_nstep
    use shr_const_mod,  only: SHR_CONST_PI
    use control_mod,    only: qsplit
    use apply_iop_forcing_mod, only:advance_iop_forcing, advance_iop_nudging

    integer :: n,nets,nete
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (TimeLevel_t), intent(in)       :: tl
    logical :: t_before_advance

    integer :: tl_qdp_np0,tl_qdp_np1
    integer :: ie,k,i,j,t,ii,jj,m
    real (r8), dimension(nlev)  :: p
    real (r8) ::dt

    integer ::nelemd_todo, np_todo
    logical ::scm_multcols = .false.
    logical ::iop_nudge_tq = .false.
    real (r8), dimension(nlev,pcnst) :: stateQ_in, q_update, q_phys_frc
    real (r8), dimension(nlev) :: t_phys_frc, t_update, u_update, v_update
    real (r8), dimension(nlev) :: t_in, u_in, v_in
    real (r8), dimension(nlev) :: relaxt, relaxq
    real (r8), dimension(nlev) :: tdiff_dyn, qdiff_dyn
    real (r8), dimension(npsq,nlev) :: tdiff_out, qdiff_out
    real (r8) :: dpscm(nlev)

!----------------------------------------------------------------------- 

    tl_f = tl%n0

    call TimeLevel_Qdp(tl, qsplit, tl_fqdp)

    ! For SCM only one column is considered
    ie = 35
    ii=3
    jj=4

    dt = get_step_size()

    ! Set initial profiles for current column
    do m=1,pcnst
       stateQ_in(:nlev,m) =             elem(ie)%state%Qdp(ii,jj,:nlev,m,tl_fqdp)/elem(ie)%state%dp3d(ii,jj,:nlev,tl_f)
    end do
    t_in(:nlev) = elem(ie)%state%T(ii,jj,:nlev,tl_f)
    u_in(:nlev) = elem(ie)%state%v(ii,jj,1,:nlev,tl_f)
    v_in(:nlev) = elem(ie)%state%v(ii,jj,2,:nlev,tl_f)
    
!!$    if (.not. use_3dfrc ) then
!!$       t_phys_frc(:) = 0.0_r8
!!$    else
    t_phys_frc(:)   = elem(ie)%derived%fT(ii,jj,:)
    q_phys_frc(:,:) = elem(ie)%derived%fQ(ii,jj,:,:)/dt
!!$    endif

    ! Call the main subroutine to update t, q, u, and v according to
    !  large scale forcing as specified in IOP file.
    call advance_iop_forcing(dt,elem(ie)%state%psdry(ii,jj),& ! In
         u_in,v_in,t_in,stateQ_in,t_phys_frc, q_phys_frc, hvcoord, &            ! In
         u_update,v_update,t_update,q_update)                      ! Out
    
    ! Nudge to observations if desired, for T & Q only if in SCM mode
    if (iop_nudge_tq ) then
       call advance_iop_nudging(dt,elem(ie)%state%psdry(ii,jj),& ! In
            t_update,q_update(:,1), hvcoord, &                   ! Inn
            t_update,q_update(:,1),relaxt,relaxq)                ! Out
    endif

    if (use_3dfrc) then    ! vertical remap of dynamics not run need to update state%dp3d using new psdry
       do k=1,nlev
          elem(ie)%state%dp3d(ii,jj,k,tl_f) = (hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 + (hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem(ie)%state%psdry(ii,jj)
       end do
    end if

    ! Update qdp using new dp3d
    do m=1,pcnst
       ! Update the Qdp array
       elem(ie)%state%Qdp(ii,jj,:nlev,m,tl_fqdp) = &
            q_update(:nlev,m) * elem(ie)%state%dp3d(ii,jj,:nlev,tl_f)
    enddo

    ! Update prognostic variables to the current values
    elem(ie)%state%T(ii,jj,:,tl_f) = t_update(:)
    elem(ie)%state%v(ii,jj,1,:,tl_f) = u_update(:)
    elem(ie)%state%v(ii,jj,2,:,tl_f) = v_update(:)

    ! Evaluate the differences in state information from observed
    !  (done for diganostic purposes only)
    do k = 1, nlev
       tdiff_dyn(k) = t_update(k)   - tobs(k)
       qdiff_dyn(k) = q_update(k,1) - qobs(k)
    end do

    ! Add various diganostic outfld calls
    call outfld('TDIFF',tdiff_dyn,1,begchunk)
    call outfld('QDIFF',qdiff_dyn,1,begchunk)
    call outfld('TOBS',tobs,1,begchunk)
    call outfld('QOBS',qobs,1,begchunk)
    call outfld('DIVQ',divq,1,begchunk)
    call outfld('DIVT',divt,1,begchunk)
    call outfld('DIVQ3D',divq3d,1,begchunk)
    call outfld('DIVT3D',divt3d,1,begchunk)
    call outfld('PRECOBS',precobs,1,begchunk)
    call outfld('LHFLXOBS',lhflxobs,1,begchunk)
    call outfld('SHFLXOBS',shflxobs,1,begchunk)

    call outfld('TRELAX',relaxt,1,begchunk)
    call outfld('QRELAX',relaxq,1,begchunk)


    end subroutine apply_SC_forcing
!=========================================================================
    subroutine iop_broadcast()
      
      !---------------------------------------------------------
      ! Purpose: When running DP-CRM, broadcast relevant logical 
      !   flags and data to all processors
      !----------------------------------------------------------
      
      use spmd_utils,   only: mpi_logical, mpi_real8, masterproc, iam, mpicom, mstrid=>masterprocid
      use dimensions_mod,         only: nlev

      integer :: ierr
#ifdef SPMD  
      
      call mpi_bcast(have_ps,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_tg,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_lhflx,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_shflx,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_t,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_q,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_u,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_v,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_omega,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_cldliq,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_divt,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_divq,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_divt3d,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(have_divq3d,1,mpi_logical,mstrid,mpicom,ierr)
      call mpi_bcast(use_3dfrc,1,mpi_logical,mstrid,mpicom,ierr)
      
      call mpi_bcast(psobs,1,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(tground,1,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(lhflxobs,1,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(shflxobs,1,mpi_real8,mstrid,mpicom,ierr)
      
      call mpi_bcast(tobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(qobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(uobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(vobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(cldliqobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(wfld,nlev,mpi_real8,mstrid,mpicom,ierr) 
      
      call mpi_bcast(divt,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(divq,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(divt3d,nlev,mpi_real8,mstrid,mpicom,ierr)
      call mpi_bcast(divq3d,nlev,mpi_real8,mstrid,mpicom,ierr)
      
#endif
      
    end subroutine iop_broadcast
    
 end module se_single_column_mod
