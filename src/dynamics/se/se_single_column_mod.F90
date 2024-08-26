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
                                  use_3dfrc,scmlat,scmlon
use constituents,           only: cnst_get_ind, pcnst
use dimensions_mod,         only: nelemd, np, nlev, qsize
use time_manager,           only: get_nstep, is_first_step, get_step_size, is_first_restart_step
use ppgrid,                 only: begchunk
use se_dyn_time_mod,        only: timelevel_qdp
use cam_history,            only: outfld

implicit none

private
save

public scm_setinitial
public scm_setfield
public apply_SC_forcing
public iop_broadcast
public scm_dyn_grid_indicies

integer, public :: indx_scm, ie_scm, i_scm, j_scm

integer         :: tl_f, tl_fqdp, thelev

!=========================================================================
contains
!=========================================================================

subroutine scm_setinitial(elem)

  use dyn_grid,     only: TimeLevel
  use control_mod,  only: qsplit

  implicit none

  type(element_t), intent(inout) :: elem(:)

  integer              :: k
  integer              :: inumliq, inumice, icldliq, icldice

  call scm_dyn_grid_indicies(elem,scmlat,scmlon,ie_scm,i_scm,j_scm,indx_scm)

  tl_f = timelevel%n0
  call TimeLevel_Qdp(timelevel, qsplit, tl_fqdp)

  if (.not. use_camiop .and. get_nstep()  ==  0) then
    call cnst_get_ind('NUMLIQ', inumliq, abort=.false.)
    call cnst_get_ind('NUMICE', inumice, abort=.false.)
    call cnst_get_ind('CLDLIQ', icldliq)
    call cnst_get_ind('CLDICE', icldice)

    ! Find level where tobs is no longer zero
    thelev=minloc(abs(tobs), 1, mask=abs(tobs) > 0)

    if (get_nstep()  <=  1) then
       do k=1,thelev-1
          tobs(k)=elem(ie_scm)%state%T(i_scm,j_scm,k,tl_f)
          qobs(k)=elem(ie_scm)%state%qdp(i_scm,j_scm,k,1,tl_fqdp)/elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
       enddo
    else
       tobs(:)=elem(ie_scm)%state%T(i_scm,j_scm,:,tl_f)
       qobs(:)=elem(ie_scm)%state%qdp(i_scm,j_scm,:,1,tl_fqdp)/elem(ie_scm)%state%dp3d(i_scm,j_scm,:,tl_f)
    endif

    if (get_nstep()  ==  0) then
       do k=thelev, NLEV
          if (have_t) elem(ie_scm)%state%T(i_scm,j_scm,k,tl_f)=tobs(k)
          if (have_q) elem(ie_scm)%state%qdp(i_scm,j_scm,k,1,tl_fqdp)=qobs(k)*elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
       enddo

       do k=1,NLEV
          if (have_ps) elem(ie_scm)%state%psdry(i_scm,j_scm) = psobs
          if (have_u) elem(ie_scm)%state%v(i_scm,j_scm,1,k,tl_f) = uobs(k)
          if (have_v) elem(ie_scm)%state%v(i_scm,j_scm,2,k,tl_f) = vobs(k)
          if (have_numliq) elem(ie_scm)%state%qdp(i_scm,j_scm,k,inumliq,tl_fqdp) = &
               numliqobs(k)*elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
          if (have_cldliq) elem(ie_scm)%state%qdp(i_scm,j_scm,k,icldliq,tl_fqdp) = &
               cldliqobs(k)*elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
          if (have_numice) elem(ie_scm)%state%qdp(i_scm,j_scm,k,inumice,tl_fqdp) = &
               numiceobs(k)*elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
          if (have_cldice) elem(ie_scm)%state%qdp(i_scm,j_scm,k,icldice,tl_fqdp) = &
               cldiceobs(k)*elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
          if (have_omega) elem(ie_scm)%derived%omega(i_scm,j_scm,k) = wfld(k)
       enddo

    endif

 endif

end subroutine scm_setinitial

subroutine scm_setfield(elem,iop_update_phase1)

!---------------------------------------------------------
! Purpose: Update various fields based on available data
!   provided by IOP file
!----------------------------------------------------------

  use control_mod,  only: qsplit
  use dyn_grid,     only: TimeLevel

  implicit none

  logical, intent(in) :: iop_update_phase1
  type(element_t), intent(inout) :: elem(:)

  integer              :: k
  integer              :: tl_f, tl_fqdp

  tl_f = timelevel%n0
  call TimeLevel_Qdp(timelevel, qsplit, tl_fqdp)

  if (have_ps .and. use_camiop .and. .not. iop_update_phase1) elem(ie_scm)%state%psdry(:,:) = psobs
  if (have_ps .and. .not. use_camiop) elem(ie_scm)%state%psdry(:,:) = psobs
  do k=1, NLEV
     if (have_omega .and. iop_update_phase1) elem(ie_scm)%derived%omega(:,:,k)=wfld(k)  !     set t to tobs at first
     if (k < thelev) then
        tobs(k) = elem(ie_scm)%state%T(i_scm,j_scm,k,tl_f)
        qobs(k) = elem(ie_scm)%state%qdp(i_scm,j_scm,k,1,tl_fqdp)/elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f)
        uobs(k) = elem(ie_scm)%state%v(i_scm,j_scm,1,k,tl_f)
        vobs(k) = elem(ie_scm)%state%v(i_scm,j_scm,2,k,tl_f)
     end if
  end do

end subroutine scm_setfield

subroutine apply_SC_forcing(elem,hvcoord,tl,n,t_before_advance)
!
    use scamMod,        only: single_column, use_3dfrc
    use hybvcoord_mod,  only: hvcoord_t
    use se_dyn_time_mod,only: TimeLevel_t
    use control_mod,    only: qsplit
    use apply_iop_forcing_mod, only:advance_iop_forcing, advance_iop_nudging

    type (element_t), intent(inout), target :: elem(:)
    type (hvcoord_t), intent(in)            :: hvcoord
    type (TimeLevel_t), intent(in)          :: tl
    logical, intent(in)                     :: t_before_advance
    integer, intent(in)                     :: n

    integer                                 :: k, m
    real (r8)                               :: dt
    logical                                 :: iop_nudge_tq = .false.
    real (r8), dimension(nlev,pcnst)        :: stateQ_in, q_update, q_phys_frc
    real (r8), dimension(nlev)              :: t_phys_frc, t_update, u_update, v_update
    real (r8), dimension(nlev)              :: t_in, u_in, v_in
    real (r8), dimension(nlev)              :: relaxt, relaxq
    real (r8), dimension(nlev)              :: tdiff_dyn, qdiff_dyn

!-----------------------------------------------------------------------

    tl_f = tl%n0

    call TimeLevel_Qdp(tl, qsplit, tl_fqdp)

    dt = get_step_size()

    ! Set initial profiles for current column
    do m=1,pcnst
       stateQ_in(:nlev,m) = elem(ie_scm)%state%Qdp(i_scm,j_scm,:nlev,m,tl_fqdp)/elem(ie_scm)%state%dp3d(i_scm,j_scm,:nlev,tl_f)
    end do
    t_in(:nlev) = elem(ie_scm)%state%T(i_scm,j_scm,:nlev,tl_f)
    u_in(:nlev) = elem(ie_scm)%state%v(i_scm,j_scm,1,:nlev,tl_f)
    v_in(:nlev) = elem(ie_scm)%state%v(i_scm,j_scm,2,:nlev,tl_f)

    t_phys_frc(:)   = elem(ie_scm)%derived%fT(i_scm,j_scm,:)
    q_phys_frc(:,:qsize) = elem(ie_scm)%derived%fQ(i_scm,j_scm,:,:qsize)/dt

    ! Call the main subroutine to update t, q, u, and v according to
    !  large scale forcing as specified in IOP file.
    call advance_iop_forcing(dt,elem(ie_scm)%state%psdry(i_scm,j_scm),& ! In
         u_in,v_in,t_in,stateQ_in,t_phys_frc, q_phys_frc, hvcoord, &            ! In
         u_update,v_update,t_update,q_update)                      ! Out

    ! Nudge to observations if desired, for T & Q only if in SCM mode
    if (iop_nudge_tq ) then
       call advance_iop_nudging(dt,elem(ie_scm)%state%psdry(i_scm,j_scm),& ! In
            t_update,q_update,u_update,v_update, hvcoord, &                ! Inout
            relaxt,relaxq)                                                 ! Out
    endif

    if (use_3dfrc) then    ! vertical remap of dynamics not run need to update state%dp3d using new psdry
       do k=1,nlev
          elem(ie_scm)%state%dp3d(i_scm,j_scm,k,tl_f) = (hvcoord%hyai(k+1)-hvcoord%hyai(k))*hvcoord%ps0 + &
               (hvcoord%hybi(k+1)-hvcoord%hybi(k))*elem(ie_scm)%state%psdry(i_scm,j_scm)
       end do
    end if

    ! Update qdp using new dp3d
    do m=1,pcnst
       ! Update the Qdp array
       elem(ie_scm)%state%Qdp(i_scm,j_scm,:nlev,m,tl_fqdp) = &
            q_update(:nlev,m) * elem(ie_scm)%state%dp3d(i_scm,j_scm,:nlev,tl_f)
    enddo

    ! Update prognostic variables to the current values
    elem(ie_scm)%state%T(i_scm,j_scm,:,tl_f) = t_update(:)
    elem(ie_scm)%state%v(i_scm,j_scm,1,:,tl_f) = u_update(:)
    elem(ie_scm)%state%v(i_scm,j_scm,2,:,tl_f) = v_update(:)

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
      ! Purpose: Broadcast relevant logical
      !   flags and data to all processors
      !----------------------------------------------------------

      use spmd_utils,             only: mpi_logical, mpi_real8, masterproc, iam, mpicom, mstrid=>masterprocid
      use cam_abortutils,         only: endrun

      integer :: ierr
      character(len=*), parameter :: sub = 'radiation_readnl'

#ifdef SPMD
      call mpi_bcast(have_ps,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_ps")
      call mpi_bcast(have_tg,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_tg")
      call mpi_bcast(have_lhflx,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_lhflx")
      call mpi_bcast(have_shflx,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_shflx")
      call mpi_bcast(have_t,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_t")
      call mpi_bcast(have_q,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_q")
      call mpi_bcast(have_u,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_u")
      call mpi_bcast(have_v,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_v")
      call mpi_bcast(have_omega,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_omega")
      call mpi_bcast(have_cldliq,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_cldliq")
      call mpi_bcast(have_divt,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_divt")
      call mpi_bcast(have_divq,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_divq")
      call mpi_bcast(have_divt3d,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_divt3d")
      call mpi_bcast(have_divq3d,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: have_divq3d")
      call mpi_bcast(use_3dfrc,1,mpi_logical,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_3dfrc")

      call mpi_bcast(psobs,1,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: psobs")
      call mpi_bcast(tground,1,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: tground")
      call mpi_bcast(lhflxobs,1,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: lhflxobs")
      call mpi_bcast(shflxobs,1,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: shflxobs")

      call mpi_bcast(tobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: tobs")
      call mpi_bcast(qobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: qobs")
      call mpi_bcast(uobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: uobs")
      call mpi_bcast(vobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: vobs")
      call mpi_bcast(cldliqobs,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: cldliqobs")
      call mpi_bcast(wfld,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: wfld")

      call mpi_bcast(divt,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: divt")
      call mpi_bcast(divq,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: divq")
      call mpi_bcast(divt3d,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: divt3d")
      call mpi_bcast(divq3d,nlev,mpi_real8,mstrid,mpicom,ierr)
      if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: divq3d")

#endif

    end subroutine iop_broadcast

!=========================================================================
    subroutine scm_dyn_grid_indicies(elem,scmlat,scmlon,ie_scm,i_scm,j_scm,indx_scm)

      !---------------------------------------------------------
      ! Purpose: Determine closest column index in the IOP file
      !   based on the input scm latitude and longitude
      !----------------------------------------------------------

      use shr_const_mod,          only: SHR_CONST_PI
      use cam_abortutils,         only: endrun

      type(element_t), intent(in) :: elem(:)
      real (r8),       intent(in) :: scmlat,scmlon
      integer,         intent(out) :: ie_scm, j_scm, i_scm, indx_scm

      integer :: i, j, indx, ie
      real(r8) :: scmposlon, minpoint, testlat, testlon, testval
      integer :: ierr
      real(r8), parameter :: rad2deg = 180.0_r8 / SHR_CONST_PI
      character(len=*), parameter :: sub = 'scm_dyn_grid_indicies'

      ie_scm=0
      i_scm=0
      j_scm=0
      indx_scm=0
      minpoint = 1000
      scmposlon = mod(scmlon + 360._r8,360._r8)
      do ie=1, nelemd
         indx=1
         do j=1, np
            do i=1, np
               testlat=elem(ie)%spherep(i,j)%lat * rad2deg
               testlon=elem(ie)%spherep(i,j)%lon * rad2deg
               if (testlon  <  0._r8) testlon=testlon+360._r8
               testval=abs(scmlat-testlat)+abs(scmposlon-testlon)
               if (testval  <  minpoint) then
                  ie_scm=ie
                  indx_scm=indx
                  i_scm=i
                  j_scm=j
                  minpoint=testval
                  if (minpoint  <  1.e-7_r8) minpoint=0._r8
               endif
               indx=indx+1
            enddo
         enddo
      enddo

      if (ie_scm == 0 .or. i_scm == 0 .or. j_scm == 0 .or. indx_scm == 0) then
         call endrun(sub//':FATAL: Could not find closest SCM point on input datafile')
      endif

    end subroutine scm_dyn_grid_indicies

 end module se_single_column_mod
