subroutine stepout(nstatsteps)

use crmx_vars
!use rad, only: qrad
use crmx_sgs, only: tk, sgs_print
use crmx_crmtracers
use crmx_microphysics, only: micro_print
use crmx_params
implicit none	
	
integer i,j,k,ic,jc,nstatsteps
integer n
real div, divmax, divmin
real rdx, rdy, rdz, coef
integer im,jm,km
real wmax, qnmax(1), qnmax1(1)
real(kind=selected_real_kind(12)) buffer(6), buffer1(6)
real(kind=selected_real_kind(12)) qi0(nzm)

#ifdef CLUBB_CRM
real(8) buffer_e(7), buffer1_e(7)
#endif



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Print stuff out:

!call t_startf ('print_out')

if(masterproc) print *,'NSTEP = ',nstep,'    NCYCLE=',ncycle

if(mod(nstep,nprint).eq.0) then
	

 divmin=1.e20
 divmax=-1.e20
	 
 rdx = 1./dx
 rdy = 1./dy

 wmax=0.
 do k=1,nzm
  coef = rho(k)*adz(k)*dz
  rdz = 1./coef
  if(ny.ne.1) then
   do j=1,ny-1*YES3D
    jc = j+1*YES3D
    do i=1,nx-1
     ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx + (v(i,jc,k)-v(i,j,k))*rdy + &
		  (w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
   end do
  else
    j = 1
    do i=1,nx-1
    ic = i+1
     div = (u(ic,j,k)-u(i,j,k))*rdx +(w(i,j,k+1)*rhow(k+1)-w(i,j,k)*rhow(k))*rdz
     divmax = max(divmax,div)
     divmin = min(divmin,div)
     if(w(i,j,k).gt.wmax) then
	wmax=w(i,j,k)
	im=i
	jm=j
	km=k
     endif
    end do
  endif
 end do

 if(dompi) then
   buffer(1) = total_water_before
   buffer(2) = total_water_after
   buffer(3) = total_water_evap
   buffer(4) = total_water_prec
   buffer(5) = total_water_ls
#ifdef CLUBB_CRM
   buffer(6) = total_water_clubb

   buffer_e(1) = total_energy_before
   buffer_e(2) = total_energy_after
   buffer_e(3) = total_energy_evap
   buffer_e(4) = total_energy_prec
   buffer_e(5) = total_energy_ls
   buffer_e(6) = total_energy_clubb
   buffer_e(7) = total_energy_rad
#endif
   call task_sum_real8(buffer, buffer1,6)
   total_water_before = buffer1(1)
   total_water_after = buffer1(2)
   total_water_evap = buffer1(3)
   total_water_prec = buffer1(4)
   total_water_ls = buffer1(5)
#ifdef CLUBB_CRM
   total_water_clubb = buffer1(6)

   call task_sum_real8(buffer_e, buffer1_e,7)
   total_energy_before = buffer1_e(1)
   total_energy_after = buffer1_e(2)
   total_energy_evap = buffer1_e(3)
   total_energy_prec = buffer1_e(4)
   total_energy_ls = buffer1_e(5)
   total_energy_clubb = buffer1_e(6)
   total_energy_rad = buffer1_e(7)
#endif
 end if

!print*,rank,minval(u(1:nx,1:ny,:)),maxval(u(1:nx,1:ny,:))
!print*,rank,'min:',minloc(u(1:nx,1:ny,:))
!print*,rank,'max:',maxloc(u(1:nx,1:ny,:))

!if(masterproc) then

!print*,'--->',tk(27,1,1)
!print*,'tk->:'
!write(6,'(16f7.2)')((tk(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'p->:'
!write(6,'(16f7.2)')((p(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'u->:'
!write(6,'(16f7.2)')((u(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'v->:'
!write(6,'(16f7.2)')((v(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'w->:'
!write(6,'(16f7.2)')((w(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'qcl:'
!write(6,'(16f7.2)')((qcl(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'qpl:'
!write(6,'(16f7.2)')((qpl(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'qrad:'
!write(6,'(16f7.2)')((qrad(i,13,k)*3600.,i=1,16),k=30,1,-1)
!print*,'qv:'
!write(6,'(16f7.2)')((qv(i,13,k)*1000.,i=1,16),k=30,1,-1)
!print*,'tabs:'
!write(6,'(16f7.2)')((tabs(i,13,k),i=1,16),k=30,1,-1)
!
!end if

!--------------------------------------------------------
 if(masterproc) then
	
    print*,'DAY = ',day	
    write(6,*) 'NSTEP=',nstep
    write(6,*) 'div:',divmax,divmin
    if(.not.dodynamicocean) write(6,*) 'SST=',tabs_s 
    write(6,*) 'surface pressure=',pres0

 endif

 call fminmax_print('u:',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm)
 call fminmax_print('v:',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm-5)
 call fminmax_print('w:',w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz)
 call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
 call fminmax_print('t:',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tabs:',tabs,1,nx,1,ny,nzm)
 call fminmax_print('qv:',qv,1,nx,1,ny,nzm)
 if(dosgs) call sgs_print()
#ifdef CLUBB_CRM
 if(docloud.or.doclubb) then
#else
 if(docloud) then
#endif /*CLUBB_CRM*/
   call fminmax_print('qcl:',qcl,1,nx,1,ny,nzm)
   call fminmax_print('qci:',qci,1,nx,1,ny,nzm)
   call micro_print()
 end if
 if(doprecip) then
   call fminmax_print('qpl:',qpl,1,nx,1,ny,nzm)
   call fminmax_print('qpi:',qpi,1,nx,1,ny,nzm)
 end if
! if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',qrad*86400.,1,nx,1,ny,nzm)
 if(dotracers) then
   do k=1,ntracers
     call fminmax_print(trim(tracername(k))//':',tracer(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
   end do
 end if
 call fminmax_print('shf:',fluxbt*cp*rhow(1),1,nx,1,ny,1)
 call fminmax_print('lhf:',fluxbq*lcond*rhow(1),1,nx,1,ny,1)
 call fminmax_print('uw:',fluxbu,1,nx,1,ny,1)
 call fminmax_print('vw:',fluxbv,1,nx,1,ny,1)
 call fminmax_print('sst:',sstxy,0,nx,1-YES3D,ny,1)

end if ! (mod(nstep,nprint).eq.0)

!call t_stopf ('print_out')

end
