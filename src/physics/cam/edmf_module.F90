module edmf_module

! =============================================================================== !
! Mass-flux module for use with CLUBB                                             !
! Together (CLUBB+MF) they comprise a eddy-diffusivity mass-flux approach (EDMF)  !
! =============================================================================== !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: integrate_mf, &
            clubb_mf_readnl

  real(r8) :: clubb_mf_L0   = 0._r8
  real(r8) :: clubb_mf_ent0 = 0._r8
  real(r8) :: clubb_mf_wa   = 0._r8
  real(r8) :: clubb_mf_wb   = 0._r8
  integer  :: clubb_mf_nup  = 0

  contains

  subroutine clubb_mf_readnl(nlfile)

  ! =============================================================================== !
  ! MF namelists                                                                    !
  ! =============================================================================== !

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use clubb_api_module,only: l_stats, l_output_rad_files
    use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical, mpi_real8, mpi_integer

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    character(len=*), parameter :: sub = 'clubb_mf_readnl'

    integer :: iunit, read_status, ierr

    namelist /clubb_mf_nl/ clubb_mf_L0, clubb_mf_ent0, clubb_mf_wa, clubb_mf_wb, clubb_mf_nup

    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' )

      call find_group_name(iunit, 'clubb_mf_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_mf_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_mf_readnl:  error reading namelist')
         end if
      end if

      close(unit=iunit)
      call freeunit(iunit)
    end if

    call mpi_bcast(clubb_mf_L0, 1, mpi_real8,   mstrid, mpicom, ierr) 
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_L0")
    call mpi_bcast(clubb_mf_ent0, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_ent0")
    call mpi_bcast(clubb_mf_wa, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_wa")
    call mpi_bcast(clubb_mf_wb, 1, mpi_real8,   mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_wb")
    call mpi_bcast(clubb_mf_nup, 1, mpi_integer,mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: clubb_mf_nup")

  end subroutine clubb_mf_readnl

  subroutine integrate_mf( nz,      dzt,     zm,      p_zm,      iexner_zm,         & ! input
                                                      p_zt,      iexner_zt,         & ! input
                           u,       v,       thl,     qt,        thv,               & ! input
                                             thl_zm,  qt_zm,                        & ! input
                                             wthl,    wqt,       pblh,              & ! input
                           dry_a,   moist_a,                                        & ! output - plume diagnostics
                           dry_w,   moist_w,                                        & ! output - plume diagnostics
                           dry_qt,  moist_qt,                                       & ! output - plume diagnostics
                           dry_thl, moist_thl,                                      & ! output - plume diagnostics
                           dry_u,   moist_u,                                        & ! output - plume diagnostics
                           dry_v,   moist_v,                                        & ! output - plume diagnostics
                                    moist_qc,                                       & ! output - plume diagnostics
                           ae,      aw,                                             & ! output - diagnosed fluxes BEFORE mean field update
                           awthl,   awqt,                                           & ! output - diagnosed fluxes BEFORE mean field update
                           awql,    awqi,                                           & ! output - diagnosed fluxes BEFORE mean field update
                           awu,     awv,                                            & ! output - diagnosed fluxes BEFORE mean field update
                           thlflx,  qtflx )                                           ! output - variables needed for solver

  ! ================================================================================= !
  ! Mass-flux algorithm                                                               !
  !                                                                                   !
  ! Provides rtm and thl fluxes due to mass flux ensemble,                            !
  ! which are fed into the mixed explicit/implicit clubb solver as explicit terms     !
  !                                                                                   !
  ! Variables needed for solver                                                       !
  ! ae = sum_i (1-a_i)                                                                !
  ! aw3 = sum (a_i w_i)                                                               ! 
  ! aws3 = sum (a_i w_i*s_i); s=thl*cp                                                !
  ! aws3,awqv3,awql3,awqi3,awu3,awv3 similar as above except for different variables  !
  !                                                                                   !
  ! Mass flux variables are computed on edges (i.e. momentum grid):                   !
  ! upa,upw,upqt,...                                                                  !
  ! dry_a,moist_a,dry_w,moist_w, ...                                                  !
  !                                                                                   ! 
  ! In CLUBB (unlike CAM) nlevs of momentum grid = nlevs of thermodynamic grid,       !
  ! due to a subsurface thermodynamic layer. To avoid confusion, below the variables  !  
  ! are grouped by the grid they are on.                                              !
  !                                                                                   !
  ! *note that state on the lowest thermo level is equal to state on the lowest       !
  ! momentum level due to state_zt(1) = state_zt(2), and lowest momentum level        !
  ! is a weighted combination of the lowest two thermodynamic levels.                 !
  !                                                                                   !
  ! ---------------------------------Authors----------------------------------------  !
  ! Marcin Kurowski, JPL                                                              !
  ! Modified heavily by Mikael Witte, UCLA/JPL for implementation in CESM2/E3SM       !
  ! Additional modifications by Adam Herrington, NCAR                                 !
  ! ================================================================================= !

     use physconst,          only: rair, cpair, gravit, latvap, latice, zvir

     real(r8), dimension(nz), intent(in) :: u,      v,            & ! thermodynamic grid
                                            thl,    thv,          & ! thermodynamic grid
                                            qt,                   & ! thermodynamic grid
                                            dzt,                  & ! thermodynamic grid
                                            p_zt,   iexner_zt,    & ! thermodynamic grid
                                            thl_zm,               & ! momentum grid
                                            qt_zm,                & ! momentum grid
                                            zm,                   & ! momentum grid
                                            p_zm,   iexner_zm       ! momentum grid

     integer,  intent(in)                :: nz
     real(r8), intent(in)                :: wthl,wqt
     real(r8), intent(inout)             :: pblh

     real(r8),dimension(nz), intent(out) :: dry_a,   moist_a,     & ! momentum grid
                                            dry_w,   moist_w,     & ! momentum grid
                                            dry_qt,  moist_qt,    & ! momentum grid
                                            dry_thl, moist_thl,   & ! momentum grid
                                            dry_u,   moist_u,     & ! momentum grid
                                            dry_v,   moist_v,     & ! momentum grid
                                                     moist_qc,    & ! momentum grid
                                            ae,      aw,          & ! momentum grid
                                            awthl,   awqt,        & ! momentum grid
                                            awql,    awqi,        & ! momentum grid
                                            awu,     awv,         & ! momentum grid
                                            thlflx,  qtflx          ! momentum grid

     ! =============================================================================== !
     ! INTERNAL VARIABLES
     !
     ! sums over all plumes
     real(r8), dimension(nz)              :: moist_th, dry_th,         & ! momentum grid
                                             awqv,     awth              ! momentum grid
     !
     ! updraft properties
     real(r8), dimension(nz,clubb_mf_nup) :: upw,      upa,            & ! momentum grid
                                             upqt,     upqc,           & ! momentum grid
                                             upqv,     upqs,           & ! momentum grid
                                             upql,     upqi,           & ! momentum grid
                                             upth,     upthv,          & ! momentum grid
                                                       upthl,          & ! momentum grid
                                             upu,      upv               ! momentum grid 
     !
     ! entrainment profiles
     real(r8), dimension(nz,clubb_mf_nup) :: ent,      entf              ! thermodynamic grid
     integer,  dimension(nz,clubb_mf_nup) :: enti                        ! thermodynamic grid
     ! 
     ! other variables
     integer                              :: k,i,ih
     real(r8), dimension(clubb_mf_nup)    :: zcb
     real(r8)                             :: zcb_unset,                &
                                             wthv,                     &
                                             wstar,  qstar,   thlstar, & 
                                             sigmaw, sigmaqt, sigmathl,&
                                                     wmin,    wmax,    & 
                                             wlv,    wtv,     wp,      & 
                                             B,                        & ! thermodynamic grid
                                             entexp, entexpu, entw,    & ! thermodynamic grid
                                             thln,   thvn,    thn,     & ! momentum grid
                                             qtn,    qsn,              & ! momentum grid
                                             qcn,    qln,     qin,     & ! momentum grid
                                             un,     vn,      wn2,     & ! momentum grid
                                             lmixn,                    & ! momentum grid
                                             supqt,  supthl              ! thermodynamic grid
     !
     ! parameters defining initial conditions for updrafts
     real(r8),parameter                   :: pwmin = 1.5_r8,           &
                                             pwmax = 3._r8

     !
     ! alpha, z-scores after Suselj etal 2019
     real(r8),parameter                   :: alphw   = 0.572_r8,       &
                                             alphqt  = 2.890_r8,       &     
                                             alphthl = 2.890_r8
     !
     ! w' covariance after Suselj etal 2019
     real(r8),parameter                   :: cwqt  = 0.32_r8,          &
                                             cwthl = 0.58_r8
     !
     ! min values to avoid singularities
     real(r8),parameter                   :: wstarmin = 1.e-3_r8,      &
                                             pblhmin  = 100._r8
     !
     ! to condensate or not to condensate
     logical                              :: do_condensation = .true.
     !
     ! to precip or not to precip
     logical                              :: do_precip = .false.
     !
     ! to debug flag
     logical                              :: debug = .false.

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! INITIALIZE OUTPUT VARIABLES
     ! set updraft properties to zero
     dry_a     = 0._r8
     moist_a   = 0._r8
     dry_w     = 0._r8
     moist_w   = 0._r8
     dry_qt    = 0._r8
     moist_qt  = 0._r8
     dry_thl   = 0._r8
     moist_thl = 0._r8
     dry_u     = 0._r8
     moist_u   = 0._r8
     dry_v     = 0._r8
     moist_v   = 0._r8
     moist_qc  = 0._r8
     ! outputs - variables needed for solver
     aw        = 0._r8
     awthl     = 0._r8
     awqt      = 0._r8
     awqv      = 0._r8
     awql      = 0._r8
     awqi      = 0._r8
     awu       = 0._r8
     awv       = 0._r8
     thlflx    = 0._r8
     qtflx     = 0._r8

     ! this is the environmental area - by default 1.
     ae = 1._r8

     ! START MAIN COMPUTATION
     upw   = 0._r8
     upthl = 0._r8
     upthv = 0._r8
     upqt  = 0._r8
     upa   = 0._r8
     upu   = 0._r8
     upv   = 0._r8
     upqc  = 0._r8
     ent   = 0._r8
     upth  = 0._r8
     upql  = 0._r8
     upqi  = 0._r8
     upqv  = 0._r8

     ! unique identifier
     zcb_unset = 9999999._r8
     zcb       = zcb_unset

     pblh = max(pblh,pblhmin)
     wthv = wthl+zvir*thv(1)*wqt

     ! if surface buoyancy is positive then do mass-flux
     if ( wthv > 0._r8 ) then

       ! get entrainment coefficient, dz/L0
       do i=1,clubb_mf_nup
         do k=2,nz
           entf(k,i) = dzt(k) / clubb_mf_L0
         enddo
       enddo

       ! get poisson, P(dz/L0)
       if (debug) then
         enti(:,:) = 4
       else
         call poisson( nz, clubb_mf_nup, entf, enti, thl(nz))
       end if

       ! get entrainment, ent=ent0/dz*P(dz/L0)
       do i=1,clubb_mf_nup
         do k=2,nz
           ent(k,i) = real( enti(k,i))*clubb_mf_ent0/dzt(k)
         enddo
       enddo

       ! get surface conditions
       wstar   = max( wstarmin, (gravit/thv(1)*wthv*pblh)**(1./3.) )
       qstar   = wqt / wstar
       thlstar = wthl / wstar

       sigmaw   = alphw * wstar
       sigmaqt  = alphqt * abs(qstar)
       sigmathl = alphthl * abs(thlstar)

       wmin = sigmaw * pwmin
       wmax = sigmaw * pwmax

       do i=1,clubb_mf_nup

         wlv = wmin + (wmax-wmin) / (real(clubb_mf_nup)) * (real(i)-1.)
         wtv = wmin + (wmax-wmin) / (real(clubb_mf_nup)) * real(i)

         upw(1,i) = 0.5_r8 * (wlv+wtv)
         upa(1,i) = 0.5_r8 * erf( wtv/(sqrt(2.)*sigmaw) ) &
                    - 0.5_r8 * erf( wlv/(sqrt(2.)*sigmaw) )

         upu(1,i) = u(1)
         upv(1,i) = v(1)

         upqt(1,i)  = qt(1)  + cwqt * upw(1,i) * sigmaqt/sigmaw
         upthl(1,i) = thl(1) + cwthl * upw(1,i) * sigmathl/sigmaw

         ! get cloud, lowest momentum level 
         if (do_condensation) then
           call condensation_mf(upqt(1,i), upthl(1,i), p_zm(1), iexner_zm(1), &
                                thvn, qcn, thn, qln, qin, qsn, lmixn)
           upthv(1,i) = thvn
           upqc(1,i)  = qcn
           upql(1,i)  = qln
           upqi(1,i)  = qin
           upqs(1,i)  = qsn
           if (qcn.gt.0._r8) zcb(i) = zm(1)
         else
           ! assume no cldliq
           upqc(1,i)  = 0._r8
           upthv(1,i) = upthl(1,i)*(1._r8+zvir*upqt(1,i))
           ! assume sigmathl = sigmath
           !!upthv(1,i) = thv(1) + 0.58_r8 * upw(1,i) * sigmathl/sigmaw
         end if

       enddo

       ! get updraft properties
       do i=1,clubb_mf_nup
         do k=1,nz-1

           ! get microphysics, autoconversion
           if (do_precip .and. upqc(k,i).gt.0._r8) then
             call precip_mf(upqs(k,i),upqt(k,i),upw(k,i),dzt(k+1),zm(k+1)-zcb(i),supqt)

             supthl = -1._r8*lmixn*supqt*iexner_zt(k+1)/cpair
           else
             supqt  = 0._r8
             supthl = 0._r8
           end if

           ! integrate updraft
           entexp  = exp(-ent(k+1,i)*dzt(k+1))
           entexpu = exp(-ent(k+1,i)*dzt(k+1)/3._r8)
           
           qtn  = qt(k+1) *(1.-entexp ) + upqt (k,i)*entexp + supqt
           thln = thl(k+1)*(1.-entexp ) + upthl(k,i)*entexp + supthl
           un   = u(k+1)  *(1.-entexpu) + upu  (k,i)*entexpu
           vn   = v(k+1)  *(1.-entexpu) + upv  (k,i)*entexpu

           ! get cloud, momentum levels
           if (do_condensation) then
             call condensation_mf(qtn, thln, p_zm(k+1), iexner_zm(k+1), &
                                  thvn, qcn, thn, qln, qin, qsn, lmixn)
             if (zcb(i).eq.zcb_unset .and. qcn.gt.0._r8) zcb(i) = zm(k+1)
           else
             thvn = thln*(1._r8+zvir*qtn)
           end if

           ! get buoyancy
           B=gravit*(0.5_r8*(thvn + upthv(k,i))/thv(k+1)-1.)
           if (debug) then
             if ( masterproc ) then
               write(iulog,*) "B(k,i), k, i ", B, k, i
             end if
           end if

           ! get wn^2
           wp = clubb_mf_wb*ent(k+1,i)
           if (wp==0._r8) then
             wn2 = upw(k,i)**2._r8+2._r8*clubb_mf_wa*B*dzt(k+1)
           else
             entw = exp(-2._r8*wp*dzt(k+1))
             wn2 = entw*upw(k,i)**2._r8+clubb_mf_wa*B/(clubb_mf_wb*ent(k+1,i))*(1._r8-entw)
           end if

           if (wn2>0._r8) then
             upw(k+1,i)   = sqrt(wn2)
             upthv(k+1,i) = thvn
             upthl(k+1,i) = thln
             upqt(k+1,i)  = qtn
             upqc(k+1,i)  = qcn
             upqs(k+1,i)  = qsn
             upu(k+1,i)   = un
             upv(k+1,i)   = vn
             upa(k+1,i)   = upa(k,i)
             upql(k+1,i)  = qln
             upqi(k+1,i)  = qin
             upqv(k+1,i)  = qtn - qcn
           else
             exit
           end if
         enddo
       enddo

       ! writing updraft properties for output
       ! all variables, except Areas are now multipled by the area
       do k=1,nz

         ! first sum over all i-updrafts
         do i=1,clubb_mf_nup
           if (upqc(k,i)>0._r8) then
             moist_a(k)   = moist_a(k)   + upa(k,i)
             moist_w(k)   = moist_w(k)   + upa(k,i)*upw(k,i)
             moist_qt(k)  = moist_qt(k)  + upa(k,i)*upqt(k,i)
             moist_thl(k) = moist_thl(k) + upa(k,i)*upthl(k,i)
             moist_u(k)   = moist_u(k)   + upa(k,i)*upu(k,i)
             moist_v(k)   = moist_v(k)   + upa(k,i)*upv(k,i)
             moist_qc(k)  = moist_qc(k)  + upa(k,i)*upqc(k,i)
           else
             dry_a(k)     = dry_a(k)     + upa(k,i)
             dry_w(k)     = dry_w(k)     + upa(k,i)*upw(k,i)
             dry_qt(k)    = dry_qt(k)    + upa(k,i)*upqt(k,i)
             dry_thl(k)   = dry_thl(k)   + upa(k,i)*upthl(k,i)
             dry_u(k)     = dry_u(k)     + upa(k,i)*upu(k,i)
             dry_v(k)     = dry_v(k)     + upa(k,i)*upv(k,i)
           endif
         enddo

         if ( dry_a(k) > 0._r8 ) then
           dry_w(k)   = dry_w(k)   / dry_a(k)
           dry_qt(k)  = dry_qt(k)  / dry_a(k)
           dry_thl(k) = dry_thl(k) / dry_a(k)
           dry_u(k)   = dry_u(k)   / dry_a(k)
           dry_v(k)   = dry_v(k)   / dry_a(k)
         else
           dry_w(k)   = 0._r8
           dry_qt(k)  = 0._r8
           dry_thl(k) = 0._r8
           dry_u(k)   = 0._r8
           dry_v(k)   = 0._r8
         endif

         if ( moist_a(k) > 0._r8 ) then
           moist_w(k)   = moist_w(k)   / moist_a(k)
           moist_qt(k)  = moist_qt(k)  / moist_a(k)
           moist_thl(k) = moist_thl(k) / moist_a(k)
           moist_u(k)   = moist_u(k)   / moist_a(k)
           moist_v(k)   = moist_v(k)   / moist_a(k)
           moist_qc(k)  = moist_qc(k)  / moist_a(k)
         else
           moist_w(k)   = 0._r8
           moist_qt(k)  = 0._r8
           moist_thl(k) = 0._r8
           moist_u(k)   = 0._r8
           moist_v(k)   = 0._r8
           moist_qc(k)  = 0._r8
         endif

       enddo

       do k=1,nz
         do i=1,clubb_mf_nup
           ae  (k) = ae  (k) - upa(k,i)
           aw  (k) = aw  (k) + upa(k,i)*upw(k,i)
           awu (k) = awu (k) + upa(k,i)*upw(k,i)*upu(k,i)
           awv (k) = awv (k) + upa(k,i)*upw(k,i)*upv(k,i)
           awthl(k)= awthl(k)+ upa(k,i)*upw(k,i)*upthl(k,i) 
           awqt(k) = awqt(k) + upa(k,i)*upw(k,i)*upqt(k,i)
           awqv(k) = awqv(k) + upa(k,i)*upw(k,i)*upqv(k,i)
           awql(k) = awql(k) + upa(k,i)*upw(k,i)*upql(k,i)
           awqi(k) = awqi(k) + upa(k,i)*upw(k,i)*upqi(k,i)
         enddo
       enddo

       do k=1,nz
         thlflx(k)= awthl(k) - aw(k)*thl_zm(k)
         qtflx(k)= awqt(k) - aw(k)*qt_zm(k)
       enddo

     end if  ! ( wthv > 0.0 )

  end subroutine integrate_mf

  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi, qs, lmix )
  ! =============================================================================== !
  ! zero or one condensation for edmf: calculates thv and qc                        !
  ! =============================================================================== !
     use physconst,          only: cpair, zvir
     use wv_saturation,      only : qsat

     real(r8),intent(in) :: qt,thl,p,iex
     real(r8),intent(out):: thv,qc,th,ql,qi,qs,lmix

     !local variables
     integer  :: niter,i
     real(r8) :: diff,t,qstmp,qcold,es,wf

     ! max number of iterations
     niter=50
     ! minimum difference
     diff=2.e-5_r8

     qc=0._r8
     t=thl/iex

     !by definition:
     ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
     ! Thl = Th - L/cp*ql/Exner              (2)
     !so:
     ! Th  = Thl + L/cp*ql/Exner             (3)
     ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
     !     = Thl*Exner + L/cp*ql
     do i=1,niter
       wf = get_watf(t)
       t = thl/iex+get_alhl(wf)/cpair*qc   !as in (4)

       ! qsat, p is in pascal (check!)
       call qsat(t,p,es,qstmp)
       qcold = qc
       qc = max(0.5_r8*qc+0.5_r8*(qt-qstmp),0._r8)
       if (abs(qc-qcold)<diff) exit
     enddo

     wf = get_watf(t)
     t = thl/iex+get_alhl(wf)/cpair*qc

     call qsat(t,p,es,qs)
     qc = max(qt-qs,0._r8)
     thv = (thl+get_alhl(wf)/cpair*iex*qc)*(1._r8+zvir*(qt-qc)-qc)
     lmix = get_alhl(wf)
     th = t*iex
     qi = qc*(1._r8-wf)
     ql = qc*wf

     contains

     function get_watf(t)
       real(r8)            :: t,get_watf,tc
       real(r8), parameter :: &
                              tmax=-10._r8, &
                              tmin=-40._r8

       tc=t-273.16_r8

       if (tc>tmax) then
         get_watf=1._r8
       else if (tc<tmin) then
         get_watf=0._r8
       else
         get_watf=(tc-tmin)/(tmax-tmin);
       end if

     end function get_watf


     function get_alhl(wf)
     !latent heat of the mixture based on water fraction
       use physconst,        only : latvap , latice
       real(r8) :: get_alhl,wf

       get_alhl = wf*latvap+(1._r8-wf)*(latvap+latice)

     end function get_alhl

  end subroutine condensation_mf

  subroutine precip_mf(qs,qt,w,dz,dzcld,Supqt)
  !**********************************************************************
  ! Precipitation microphysics
  ! By Adam Herrington, after Kay Suselj
  !**********************************************************************

       real(r8),intent(in)  :: qs,qt,w,dz,dzcld
       real(r8),intent(out) :: Supqt
       ! 
       ! local vars
       real(r8) :: tauwgt, tau,       & ! time-scale vars
                   qstar                ! excess cloud liquid                   

       real(r8) :: tau0  = 15._r8,    & ! base time-scale
                   zmin  = 300._r8,   & ! small cloud thick
                   zmax  = 3000._r8,  & ! large cloud thick
                   qcmin = 0.00125_r8   ! supersat threshold 

       qstar = qs+qcmin
       
       if (qt.gt.qstar) then
         ! get precip efficiency
         tauwgt = (dzcld-zmin)/(zmax-zmin)
         tauwgt = min(max(tauwgt,0._r8),1._r8)
         tau    = tauwgt/tau0
 
         ! get source for updraft
         Supqt = (qstar-qt)*(1._r8 - exp(-1._r8*tau/w))
       else
         Supqt = 0._r8
       end if

  end subroutine precip_mf

  subroutine poisson(nz,nup,lambda,poi,state)

       integer, intent(in)                     :: nz,nup
       real(r8), intent(in)                    :: state
       real(r8), dimension(nz,nup), intent(in) :: lambda
       integer, dimension(nz,nup), intent(out) :: poi
       integer                                 :: i,j

       call set_seed_from_state(state)

       do i=1,nz
         do j=1,nup
           call knuth(lambda(i,j),poi(i,j))
         enddo
       enddo

  end subroutine poisson

  subroutine set_seed_from_state(state)
  !**********************************************************************
  ! Set a unique (but reproduceble) seed using state
  ! By Adam Herrington
  !**********************************************************************

       real(r8),intent(in)  :: state
       integer, allocatable :: seed(:)
       integer               :: i,n,tmpseed

       call random_seed(size = n)

       if (allocated(seed)) deallocate(seed)
       allocate(seed(n))

       tmpseed = int((state - int(state)) * 1000000000._r8)
       do i=1,n
         seed(i) = tmpseed
       end do

       call random_seed(put=seed)
       deallocate(seed)

  end subroutine set_seed_from_state

  subroutine knuth(lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth 
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************

       real(r8), intent(in) :: lambda
       integer, intent(out) :: kout

       !Local variables
       real(r8) :: puni, tmpuni, explam
       integer  :: k

       k = 0
       explam = exp(-1._r8*lambda)
       puni = 1._r8
       do while (puni.gt.explam)
         k = k + 1
         call random_number(tmpuni)
         puni = puni*tmpuni
       end do
       kout = k - 1

  end subroutine knuth

end module edmf_module
