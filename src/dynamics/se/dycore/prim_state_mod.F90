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

  public :: prim_printstate

CONTAINS

  subroutine prim_printstate(elem, tl,hybrid,nets,nete, fvm)
    use fvm_control_volume_mod, only: n0_fvm
    use dimensions_mod,         only: ntrac
    use constituents,           only: cnst_name
    use dimensions_mod,         only: qsize_condensate_loading,qsize_condensate_loading_idx_gll
    use cam_control_mod,        only: initial_run

    type (element_t),             intent(inout) :: elem(:)
    type(fvm_struct),             intent(inout) :: fvm(:)
    type (TimeLevel_t), target,   intent(in)    :: tl
    type (hybrid_t),              intent(in)    :: hybrid
    integer,                      intent(in)    :: nets,nete

    ! Local variables...
    integer            :: k,ie,m_cnst
    integer, parameter :: type=ORDERED

    integer, parameter :: vmax=8+2*MAX(qsize_d,ntrac_d)

    character(len=7) :: varname(vmax)

    real (kind=r8), dimension(nets:nete,vmax) :: min_local,max_local
    real (kind=r8), dimension(vmax)           :: min_p,max_p,mass,mass_chg
    real (kind=r8), dimension(np,np,nets:nete):: moist_ps

    real (kind=r8) :: tmp_gll(np,np,vmax,nets:nete),tmp_mass(vmax)!
    real (kind=r8) :: tmp_fvm(nc,nc,vmax,nets:nete)
    real (kind=r8) :: tmp_q(np,np,nlev)
    integer        :: n0, n0_qdp, q, nm, nm2
    real(kind=r8)  :: da_gll(np,np,nets:nete),da_fvm(nc,nc,nets:nete)

    !dt=tstep*qsplit
    !dt = tstep*qsplit*rsplit  ! vertical REMAP timestep
    !dynamics variables in n0 are at time =  'time': time=tl%nstep*tstep
    if (hybrid%masterthread) then
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    ! dynamics timelevels
    n0=tl%n0
    call TimeLevel_Qdp( tl, qsplit, n0_qdp)
    ! moist surface pressure
    do ie=nets,nete
      moist_ps(:,:,ie)=elem(ie)%state%psdry(:,:,n0)
      do q=1,qsize_condensate_loading
        m_cnst = qsize_condensate_loading_idx_gll(q)
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
    varname(1)         = 'U      '
    varname(2)         = 'V      '
    varname(3)         = 'T      '
    varname(4)         = 'OMEGA  '
    varname(5)         = 'PSDRY  '
    varname(6)         = 'PS     '
    nm  = 6                   ! number of vars before tracers
    nm2 = 6+statediag_numtrac ! number of state diagnostics
    do ie=nets,nete      
      min_local(ie,1)  = MINVAL(elem(ie)%state%v(:,:,1,:,n0))
      max_local(ie,1)  = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
      min_local(ie,2)  = MINVAL(elem(ie)%state%v(:,:,2,:,n0))
      max_local(ie,2)  = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))
      min_local(ie,3)  = MINVAL(elem(ie)%state%T(:,:,:,n0))
      max_local(ie,3)  = MAXVAL(elem(ie)%state%T(:,:,:,n0))
      min_local(ie,4)  = MINVAL(elem(ie)%derived%Omega(:,:,:))
      max_local(ie,4)  = MAXVAL(elem(ie)%derived%Omega(:,:,:))
      min_local(ie,5)  = MINVAL(elem(ie)%state%psdry(:,:,n0))
      max_local(ie,5)  = MAXVAL(elem(ie)%state%psdry(:,:,n0))
      min_local(ie,6)  = MINVAL(moist_ps(:,:,ie))
      max_local(ie,6)  = MAXVAL(moist_ps(:,:,ie))
      
      if (ntrac>0) then
        do q=1,statediag_numtrac
          varname(nm+q)         = TRIM(cnst_name(q))
          min_local(ie,nm+q) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0_fvm))
          max_local(ie,nm+q) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0_fvm))
        end do
      else
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
      min_p(k) = ParallelMin(min_local(:,k),hybrid)
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
    mass(1:nm) = -1.0_r8
    if (ntrac>0) then
      do ie=nets,nete
        do q=1,statediag_numtrac
          tmp_fvm(:,:,q,ie) = SUM(fvm(ie)%c(1:nc,1:nc,:,q,n0_fvm)*fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm),DIM=3)
        end do
      end do
      call global_integrals_general(tmp_fvm(:,:,1:statediag_numtrac,nets:nete),hybrid,nc,da_fvm,statediag_numtrac,&
           nets,nete,tmp_mass(1:statediag_numtrac))
    else
      do ie=nets,nete
        do q=1,statediag_numtrac
          tmp_gll(:,:,q,ie)=sum(elem(ie)%state%Qdp(:,:,:,q,n0_qdp),DIM=3)
        end do
      end do
      call global_integrals_general(tmp_gll(:,:,1:statediag_numtrac,nets:nete),hybrid,np,da_gll,statediag_numtrac,&
           nets,nete,tmp_mass(1:statediag_numtrac))
    end if
    !
    ! convert to weight in hPa
    !
    mass(nm+1:nm+statediag_numtrac)=tmp_mass(1:statediag_numtrac)*0.01_r8
    !
    ! compute dry and moist average PS 
    !
    do ie=nets,nete
      tmp_gll(:,:,1,ie)=elem(ie)%state%psdry(:,:,n0)
      tmp_gll(:,:,2,ie)=moist_ps(:,:,ie)
    enddo
    call global_integrals_general(tmp_gll(:,:,1:2,nets:nete),hybrid,np,da_gll,2,&
         nets,nete,tmp_mass(1:2))
    !
    ! convert to hPa
    !
    mass(5) = tmp_mass(1)*0.01_r8
    mass(6) = tmp_mass(2)*0.01_r8
    !
    ! compute relative mass change
    !
    if (tl%nstep==0.or..not. initial_run) then
      mass_chg(:) = 0.0_R8
      elem(nets)%derived%mass(1:statediag_numtrac)   = mass(nm+1:nm+statediag_numtrac)
      elem(nets)%derived%mass(statediag_numtrac+1)   = mass(5)
      elem(nets)%derived%mass(statediag_numtrac+2)   = mass(6)
    else
      mass_chg(:) = 0.0_r8
      do q=1,statediag_numtrac
        if (ABS(elem(nets)%derived%mass(q))<1.0e-12_r8) then
          mass_chg(nm+q) =mass(nm+q) - elem(nets)%derived%mass(q)
        else
          mass_chg(nm+q) =(mass(nm+q) - elem(nets)%derived%mass(q))/elem(nets)%derived%mass(q)
        end if
      end do
      mass_chg(5) =(mass(5) - elem(nets)%derived%mass(statediag_numtrac+1))/&
           elem(nets)%derived%mass(statediag_numtrac+1)
      mass_chg(6) =(mass(6) - elem(nets)%derived%mass(statediag_numtrac+2))/&
           elem(nets)%derived%mass(statediag_numtrac+2)
    end if
    !
    ! write diagnostics to log file
    !
    if(hybrid%masterthread) then
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
    
100 format (A8,4(E23.15))
101 format (A8,A23,A23,A23,A23)

  end subroutine prim_printstate


end module prim_state_mod
