module prim_state_mod
  use shr_kind_mod,     only: r8=>shr_kind_r8
  use cam_logfile,      only: iulog
  use dimensions_mod,   only: nlev, np, nc, qsize_d, ntrac_d
  use parallel_mod,     only: ordered
  use hybrid_mod,       only: hybrid_t
  use time_mod,         only: timelevel_t, TimeLevel_Qdp, time_at
  use control_mod,      only: qsplit, statediag_numtrac
  use global_norms_mod, only: global_integrals_general
  use element_mod,      only: element_t
  use reduction_mod,    only: parallelmax,parallelmin
  use fvm_control_volume_mod, only: fvm_struct

  implicit none
  private

  public :: prim_printstate, adjust_nsplit

CONTAINS

  subroutine prim_printstate(elem, tl,hybrid,nets,nete, fvm, omega_cn)
    use dimensions_mod,         only: ntrac
    use constituents,           only: cnst_name
    use physconst,              only: thermodynamic_active_species_idx_dycore, dry_air_species_num
    use physconst,              only: thermodynamic_active_species_num,thermodynamic_active_species_idx
    use cam_control_mod,        only: initial_run
    use time_mod,               only: tstep
    use control_mod,            only: rsplit, qsplit
    use perf_mod,       only: t_startf, t_stopf
    type (element_t),             intent(inout) :: elem(:)
    type (TimeLevel_t), target,   intent(in)    :: tl
    type (hybrid_t),              intent(in)    :: hybrid
    integer,                      intent(in)    :: nets,nete
    type(fvm_struct),             intent(inout) :: fvm(:)        
    real (kind=r8), optional,     intent(in)    :: omega_cn(2,nets:nete)
    ! Local variables...
    integer            :: k,ie,m_cnst
    integer, parameter :: type=ORDERED

    integer, parameter :: vmax=11+2*MAX(qsize_d,ntrac_d)

    character(len=10) :: varname(vmax)

    real (kind=r8), dimension(nets:nete,vmax) :: min_local,max_local
    real (kind=r8), dimension(vmax)           :: min_p,max_p,mass,mass_chg
    real (kind=r8), dimension(np,np,nets:nete):: moist_ps
    real (kind=r8), dimension(nc,nc,nets:nete):: moist_ps_fvm

    real (kind=r8) :: tmp_gll(np,np,vmax,nets:nete),tmp_mass(vmax)!
    real (kind=r8) :: tmp_fvm(nc,nc,vmax,nets:nete)
    real (kind=r8) :: tmp_q(np,np,nlev)
    integer        :: n0, n0_qdp, q, nm, nm2
    real(kind=r8)  :: da_gll(np,np,nets:nete),da_fvm(nc,nc,nets:nete)

    !dynamics variables in n0 are at time =  'time': time=tl%nstep*tstep
    if (hybrid%masterthread) then
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    ! dynamics timelevels
    n0=tl%n0
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    ! moist surface pressure
    if (ntrac>0) then
      do ie=nets,nete
        moist_ps_fvm(:,:,ie)=SUM(fvm(ie)%dp_fvm(1:nc,1:nc,:),DIM=3)
        do q=dry_air_species_num+1,thermodynamic_active_species_num
          m_cnst = thermodynamic_active_species_idx(q)
          do k=1,nlev
            moist_ps_fvm(:,:,ie) = moist_ps_fvm(:,:,ie)+&
                 fvm(ie)%dp_fvm(1:nc,1:nc,k)*fvm(ie)%c(1:nc,1:nc,k,m_cnst)
          end do
        end do
      enddo
    end if
    do ie=nets,nete
      moist_ps(:,:,ie)=elem(ie)%state%psdry(:,:)
      do q=dry_air_species_num+1,thermodynamic_active_species_num
        m_cnst = thermodynamic_active_species_idx_dycore(q)
        do k=1,nlev
          moist_ps(:,:,ie) = moist_ps(:,:,ie)+&
               elem(ie)%state%Qdp(:,:,k,m_cnst,n0_qdp)
        end do
      end do
    enddo
    ! weights/areas for global integrals
    do ie=nets,nete
      da_gll(:,:,ie) = elem(ie)%mp(:,:)*elem(ie)%metdet(:,:)
    enddo
    if (ntrac>0) then
      do ie=nets,nete
        da_fvm(:,:,ie) = fvm(ie)%area_sphere(:,:)
      enddo
    end if
    !
    !*********************************************
    !
    ! min/max of u,v,T,PS,OMEGA
    !
    !*********************************************
    !
    varname(1)           = 'U         '
    varname(2)           = 'V         '
    varname(3)           = 'T         '
    varname(4)           = 'OMEGA     '
    varname(5)           = 'OMEGA CN  '
    if (ntrac>0) then
      varname(6)         = 'PSDRY(fvm)'
      varname(7)         = 'PS(fvm)   '
      varname(8)         = 'PSDRY(gll)'
      varname(9)         = 'PS(gll)   '
      nm  = 9                  !number of vars before tracers
      nm2 = nm+statediag_numtrac!number of vars after tracers
    else
      varname(6)         = 'PSDRY     '
      varname(7)         = 'PS        '
      nm  = 7                  !number of vars before tracers
      nm2 = nm+statediag_numtrac!number of vars after tracers
    end if

    do ie=nets,nete      
      min_local(ie,1)  = MINVAL(elem(ie)%state%v(:,:,1,:,n0))
      max_local(ie,1)  = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
      min_local(ie,2)  = MINVAL(elem(ie)%state%v(:,:,2,:,n0))
      max_local(ie,2)  = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))
      min_local(ie,3)  = MINVAL(elem(ie)%state%T(:,:,:,n0))
      max_local(ie,3)  = MAXVAL(elem(ie)%state%T(:,:,:,n0))
      min_local(ie,4)  = MINVAL(elem(ie)%derived%Omega(:,:,:))
      max_local(ie,4)  = MAXVAL(elem(ie)%derived%Omega(:,:,:))
      if (present(omega_cn)) then
        min_local(ie,5)  = omega_cn(1,ie)
        max_local(ie,5)  = omega_cn(2,ie)
      else
        min_local(ie,5)  = 0.0_r8
        max_local(ie,5)  = 0.0_r8
      end if
      if (ntrac>0) then
        min_local(ie,6) = MINVAL(SUM(fvm(ie)%dp_fvm(1:nc,1:nc,:),DIM=3))
        max_local(ie,6) = MAXVAL(SUM(fvm(ie)%dp_fvm(1:nc,1:nc,:),DIM=3))
        min_local(ie,7) = MINVAL(moist_ps_fvm(:,:,ie))
        max_local(ie,7) = MINVAL(moist_ps_fvm(:,:,ie))
        min_local(ie,8)  = MINVAL(elem(ie)%state%psdry(:,:))
        max_local(ie,8)  = MAXVAL(elem(ie)%state%psdry(:,:))
        min_local(ie,9)  = MINVAL(moist_ps(:,:,ie))
        max_local(ie,9)  = MAXVAL(moist_ps(:,:,ie))      
        do q=1,statediag_numtrac
          varname(nm+q)         = TRIM(cnst_name(q))
          min_local(ie,nm+q) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,q))
          max_local(ie,nm+q) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,q))
        end do
      else
        min_local(ie,6)  = MINVAL(elem(ie)%state%psdry(:,:))
        max_local(ie,6)  = MAXVAL(elem(ie)%state%psdry(:,:))
        min_local(ie,7)  = MINVAL(moist_ps(:,:,ie))
        max_local(ie,7)  = MAXVAL(moist_ps(:,:,ie))        
        do q=1,statediag_numtrac
          varname(nm+q)         = TRIM(cnst_name(q))
          tmp_q = elem(ie)%state%Qdp(:,:,:,q,n0_qdp)/elem(ie)%state%dp3d(:,:,:,n0)
          min_local(ie,nm+q) = MINVAL(tmp_q)
          max_local(ie,nm+q) = MAXVAL(tmp_q)
        end do
      end if
      !
      ! forcing diagnostics
      !
      varname(nm2+1)         = 'FT     '
      varname(nm2+2)         = 'FM     '
      min_local(ie,nm2+1)  = MINVAL(elem(ie)%derived%FT(:,:,:))
      max_local(ie,nm2+1)  = MAXVAL(elem(ie)%derived%FT(:,:,:))
      min_local(ie,nm2+2)  = MINVAL(elem(ie)%derived%FM(:,:,:,:))
      max_local(ie,nm2+2)  = MAXVAL(elem(ie)%derived%FM(:,:,:,:))
      if (ntrac>0) then
        do q=1,statediag_numtrac
          varname(nm2+2+q)         = TRIM('F'//TRIM(cnst_name(q)))
          min_local(ie,nm2+2+q) = MINVAL(fvm(ie)%fc(1:nc,1:nc,:,q))
          max_local(ie,nm2+2+q) = MAXVAL(fvm(ie)%fc(1:nc,1:nc,:,q))
        end do
      else
        do q=1,statediag_numtrac
          varname(nm2+2+q)         = TRIM('F'//TRIM(cnst_name(q)))
          tmp_q = elem(ie)%derived%FQ(:,:,:,q)
          min_local(ie,nm2+2+q) = MINVAL(tmp_q)
          max_local(ie,nm2+2+q) = MAXVAL(tmp_q)
        end do
      end if

    end do
    !JMD This is a Thread Safe Reduction
    do k = 1, nm2+2+statediag_numtrac
       if (k==1) call t_startf('parallelMin')
      min_p(k) = ParallelMin(min_local(:,k),hybrid)
      if (k==1) call t_stopf('parallelMin')
      max_p(k) = ParallelMax(max_local(:,k),hybrid)
    end do
    !
    !*********************************************
    !
    ! Mass diagnostics
    !
    !*********************************************
    !
    ! tracers
    !
    mass = -1.0_r8
    if (ntrac>0) then
      do ie=nets,nete
        do q=1,statediag_numtrac
          tmp_fvm(:,:,q,ie) = SUM(fvm(ie)%c(1:nc,1:nc,:,q)*fvm(ie)%dp_fvm(1:nc,1:nc,:),DIM=3)
        end do
        q=statediag_numtrac+1
        tmp_fvm(:,:,q,ie) = SUM(fvm(ie)%dp_fvm(1:nc,1:nc,:),DIM=3)
        q=statediag_numtrac+2
        tmp_fvm(:,:,q,ie) = moist_ps_fvm(:,:,ie)
      end do
      call global_integrals_general(tmp_fvm(:,:,1:statediag_numtrac+2,nets:nete),hybrid,nc,da_fvm,statediag_numtrac+2,&
           nets,nete,tmp_mass(1:statediag_numtrac+2))
      mass(nm+1:nm+statediag_numtrac)=tmp_mass(1:statediag_numtrac)*0.01_r8
      mass(6:7)=tmp_mass(statediag_numtrac+1:statediag_numtrac+2)*0.01_r8
      do ie=nets,nete
        tmp_gll(:,:,1,ie)=elem(ie)%state%psdry(:,:)
        tmp_gll(:,:,2,ie)=moist_ps(:,:,ie)
      end do
      call global_integrals_general(tmp_gll(:,:,1:2,nets:nete),hybrid,np,da_gll,2,&
           nets,nete,tmp_mass(1:2))
      mass(8:9)=tmp_mass(1:2)*0.01_r8
    else
      do ie=nets,nete
        do q=1,statediag_numtrac
          tmp_gll(:,:,q,ie)=sum(elem(ie)%state%Qdp(:,:,:,q,n0_qdp),DIM=3)
        end do
        q=statediag_numtrac+1
        tmp_gll(:,:,q,ie)=elem(ie)%state%psdry(:,:)
        q=statediag_numtrac+2
        tmp_gll(:,:,q,ie)=moist_ps(:,:,ie)
      end do
      call global_integrals_general(tmp_gll(:,:,1:statediag_numtrac+2,nets:nete),hybrid,np,da_gll,statediag_numtrac+2,&
           nets,nete,tmp_mass(1:statediag_numtrac+2))
      mass(nm+1:nm+statediag_numtrac)=tmp_mass(1:statediag_numtrac)*0.01_r8
      mass(6:7)=tmp_mass(statediag_numtrac+1:statediag_numtrac+2)*0.01_r8
    end if
    !
    ! compute relative mass change
    !
    if (tl%nstep==0.or..not. initial_run) then
      mass_chg(:) = 0.0_R8
      elem(nets)%derived%mass(nm+1:nm+statediag_numtrac)   = mass(nm+1:nm+statediag_numtrac)
      if (ntrac>0) then
        elem(nets)%derived%mass(6:9)   = mass(6:9)
      else
        elem(nets)%derived%mass(6:7)   = mass(6:7)
      end if
    else
      mass_chg(:) = 0.0_r8
      do q=1,nm2!statediag_numtrac
        if (mass(q).ne.-1.0_r8) then
          if (ABS(elem(nets)%derived%mass(q))<1.0e-12_r8) then
            mass_chg(q) =mass(q) - elem(nets)%derived%mass(q)
          else
            mass_chg(q) =(mass(q) - elem(nets)%derived%mass(q))/elem(nets)%derived%mass(q)
          end if
        end if
      end do
    end if
    !
    ! write diagnostics to log file
    !
    if(hybrid%masterthread) then
      write(iulog,*)   '  '
      write(iulog,*)   'STATE DIAGNOSTICS'
      write(iulog,*)   '  '
      write(iulog,101) ' ','MIN','MAX','AVE (hPa)','REL. MASS. CHANGE'
      do k=1,nm+statediag_numtrac
        if (mass(k)==-1.0_r8) then
          write(iulog,100) varname(k),min_p(k),max_p(k)
        else
          write(iulog,100) varname(k),min_p(k),max_p(k),mass(k),mass_chg(k)
        end if
      end do
      !
      ! forcing diagnostics
      !
      write(iulog,*)   '  '
      write(iulog,*)   'FORCING DIAGNOSTICS'
      write(iulog,*)   '  '
      write(iulog,101) ' ','MIN','MAX'
      do k=nm2+1,nm2+2+statediag_numtrac
        write(iulog,100) varname(k),min_p(k),max_p(k)
      end do
    end if
    
100 format (A12,4(E23.15))
101 format (A12,A23,A23,A23,A23)

#ifdef waccm_debug
    call prim_printstate_cslam_gamma(elem, tl,hybrid,nets,nete, fvm)
#endif
    call prim_printstate_U(elem, tl,hybrid,nets,nete, fvm) 
  end subroutine prim_printstate


#ifdef waccm_debug
  subroutine prim_printstate_cslam_gamma(elem, tl,hybrid,nets,nete, fvm)
    type (element_t),             intent(inout) :: elem(:)
    type(fvm_struct),             intent(inout) :: fvm(:)
    type (TimeLevel_t), target,   intent(in)    :: tl
    type (hybrid_t),              intent(in)    :: hybrid
    integer,                      intent(in)    :: nets,nete

    ! Local variables...
    integer            :: k,ie

    real (kind=r8), dimension(nets:nete,nlev) :: max_local
    real (kind=r8), dimension(nlev)           :: max_p
    integer        :: n0, n0_qdp, q, nm, nm2

    !dt=tstep*qsplit
    !dt = tstep*qsplit*rsplit  ! vertical REMAP timestep
    !dynamics variables in n0 are at time =  'time': time=tl%nstep*tstep
    if (hybrid%masterthread) then
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    ! dynamics timelevels
    n0=tl%n0
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)

    do ie=nets,nete
      do k=1,nlev
         max_local(ie,k)  = MAXVAL(fvm(ie)%CSLAM_gamma(:,:,k,1))
      end do
    end do
    !JMD This is a Thread Safe Reduction
    do k = 1, nlev
      max_p(k) = Parallelmax(max_local(:,k),hybrid)
    end do
    if (hybrid%masterthread) then
       write(iulog,*)   '  '
       write(iulog,*)   'Gamma max'
       write(iulog,*)   '  '
       do k=1,nlev
          write(iulog,*) 'k,gamma= ',k,max_p(k)
       end do
    end if
  end subroutine prim_printstate_cslam_gamma
#endif

  subroutine adjust_nsplit(elem, tl,hybrid,nets,nete, fvm, omega_cn)
    use dimensions_mod,         only: ksponge_end
    use dimensions_mod,         only: fvm_supercycling, fvm_supercycling_jet
    use time_mod,               only: tstep
    use control_mod,            only: rsplit, qsplit
    use perf_mod,               only: t_startf, t_stopf
    use time_mod,               only: nsplit, nsplit_baseline,rsplit_baseline
    use control_mod,            only: qsplit, rsplit
    use time_manager,           only: get_step_size
    use cam_abortutils,         only: endrun
    use control_mod,    only: nu_top
    !
    type (element_t),             intent(inout) :: elem(:)    
    type (TimeLevel_t), target,   intent(in)    :: tl
    type (hybrid_t),              intent(in)    :: hybrid
    integer,                      intent(in)    :: nets,nete
    type(fvm_struct),             intent(inout) :: fvm(:)        
    real (kind=r8),               intent(in)    :: omega_cn(2,nets:nete)
    ! Local variables...
    integer            :: k,ie
    real (kind=r8), dimension(1) :: min_o
    real (kind=r8), dimension(1) :: max_o
    real (kind=r8)               :: dtime
    character(len=128)           :: errmsg
    real (kind=r8)               :: threshold=0.90_r8
    real (kind=r8)               :: max_abs_omega_cn(nets:nete)
    real (kind=r8)               :: min_abs_omega_cn(nets:nete)
    !
    ! The threshold values for when to double nsplit are empirical.
    ! In FW2000climo runs the Courant numbers are large in the sponge
    !
    ! The model was found to be stable if regular del4 is increased
    ! in the sponge and nu_top is increased (when nsplit doubles)
    !
    !
    do ie=nets,nete
      max_abs_omega_cn(ie) = MAXVAL(ABS(omega_cn(:,ie)))
    end do

    !JMD This is a Thread Safe Reduction
    do k = 1,1
      max_o(k) = ParallelMax(max_abs_omega_cn(:),hybrid)
!      min_o(k) = ParallelMin(min_abs_omega_cn(:),hybrid)
    end do
    if (max_o(1)>threshold.and.nsplit==nsplit_baseline) then
      !
      ! change vertical remap time-step
      !
       nsplit=2*nsplit_baseline
       fvm_supercycling     = rsplit
       fvm_supercycling_jet = rsplit
       nu_top=2.0_r8*nu_top       
      !
      ! write diagnostics to log file
      !
       if(hybrid%masterthread) then
          !dynamics variables in n0 are at time =  'time': time=tl%nstep*tstep
          !dt=tstep*qsplit
          !    dt_remap = tstep*qsplit*rsplit  ! vertical REMAP timestep
          !
          write(iulog,*)   'adj. nsplit: doubling nsplit; t=',Time_at(tl%nstep)/(24*3600)," [day]; max OMEGA",max_o(1)
       end if
       dtime = get_step_size()
       tstep = dtime / real(nsplit*qsplit*rsplit, r8)
       
    else if (nsplit.ne.nsplit_baseline.and.max_o(1)<0.4_r8*threshold) then
      !
      ! should nsplit be reduced again?
      !
       nsplit=nsplit_baseline
       rsplit=rsplit_baseline
       fvm_supercycling     = rsplit
       fvm_supercycling_jet = rsplit
       nu_top=nu_top/2.0_r8
       
!       nu_div_scale_top(:) = 1.0_r8
       
       dtime = get_step_size()
       tstep = dtime / real(nsplit*qsplit*rsplit, r8)
       if(hybrid%masterthread) then
         write(iulog,*)   'adj. nsplit: reset nsplit   ; t=',Time_at(tl%nstep)/(24*3600)," [day]; max OMEGA",max_o(1)
       end if
    end if
  end subroutine adjust_nsplit

  subroutine prim_printstate_U(elem, tl,hybrid,nets,nete, fvm)
    type (element_t),             intent(inout) :: elem(:)
    type(fvm_struct),             intent(inout) :: fvm(:)
    type (TimeLevel_t), target,   intent(in)    :: tl
    type (hybrid_t),              intent(in)    :: hybrid
    integer,                      intent(in)    :: nets,nete

    ! Local variables...
    integer            :: k,ie

    real (kind=r8), dimension(nets:nete,nlev) :: max_local
    real (kind=r8), dimension(nets:nete,nlev) :: min_local    
    real (kind=r8), dimension(nlev)           :: max_p
    real (kind=r8), dimension(nlev)           :: min_p
    integer        :: n0, n0_qdp, q, nm, nm2

    !dt=tstep*qsplit
    !dt = tstep*qsplit*rsplit  ! vertical REMAP timestep
    !dynamics variables in n0 are at time =  'time': time=tl%nstep*tstep
    if (hybrid%masterthread) then
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    ! dynamics timelevels
    n0=tl%n0
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)

    do ie=nets,nete
      do k=1,nlev
        max_local(ie,k)  = MAXVAL(elem(ie)%state%v(:,:,:,k,n0))
        min_local(ie,k)  = MINVAL(elem(ie)%state%v(:,:,:,k,n0))
      end do
    end do
    !JMD This is a Thread Safe Reduction
    do k = 1, nlev
      max_p(k) = Parallelmax(max_local(:,k),hybrid)
      min_p(k) = Parallelmin(min_local(:,k),hybrid)      
    end do
    if (hybrid%masterthread) then
       write(iulog,*)   '  '
       write(iulog,*)   'min/max of wind components in each layer'
       write(iulog,*)   '  '
       do k=1,nlev
          write(iulog,*) 'k,V (min max)= ',k,min_p(k),max_p(k)
       end do
    end if
  end subroutine prim_printstate_U
end module prim_state_mod
