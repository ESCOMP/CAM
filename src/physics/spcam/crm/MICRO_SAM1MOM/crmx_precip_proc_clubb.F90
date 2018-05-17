#define CLDFRAC
#ifdef CLDFRAC 
subroutine precip_proc_clubb

#ifdef CLUBB_CRM
use crmx_vars
use crmx_microphysics
use crmx_micro_params
use crmx_params
use crmx_vars, only: CF3D

implicit none

integer i,j,k
real autor, autos, accrr, accris, accrcs, accrig, accrcg
real dq, omn, omp, omg, qsatt
real pows1, pows2, powg1, powg2, powr1, powr2, tmp
real qii, qcc, qrr, qss, qgg

real cld3d(nx, ny, nzm), cldmax(nx, ny, nzm)
real cld3d_temp(nx, ny, nzm)
real cloud_frac_thresh
real qclr
real dqpsrc, dqpevp

powr1 = (3 + b_rain) / 4.
powr2 = (5 + b_rain) / 8.
pows1 = (3 + b_snow) / 4.
pows2 = (5 + b_snow) / 8.
powg1 = (3 + b_grau) / 4.
powg2 = (5 + b_grau) / 8.
      
!call t_startf ('precip_proc_clubb')

! Get cloud fraction of non-precipitating condensate 
! and precipitating condensate
cloud_frac_thresh = 0.005
do j=1, ny  
 do i=1, nx
  do k=nzm, 1, -1
    cld3d(i, j, k) = CF3D(i,j,k)
    cld3d_temp(i, j, k) = min(0.999, max(CF3D(i,j,k), cloud_frac_thresh)) 
  end do
  cldmax(i,j,nzm)=cld3d_temp(i,j,nzm)

  do k=nzm-1, 1, -1
  ! if precipitating condensate is smaller than threshold, set cldmax
  ! to cloud fraction at current level
    if(qp(i, j, k+1).ge.qp_threshold) then
       cldmax(i,j,k) = max(cldmax(i,j,k+1), cld3d_temp(i,j,k))
    else
       cldmax(i,j,k) = cld3d_temp(i,j,k)
    end if
    
!    if(cld3d(i,j,k).le.cloud_frac_thresh .and. qp(i,j,k).gt.qp_threshold) then
!       if(cldmax(i,j,k).lt.0.1) then
!         cldmax(i,j,k) = 0.50
!       end if
!    end if
  end do 
! test: assume precipitating hydrometer fill the whole grid box
!  cldmax(i,j,:) = 0.999

 end do
end do

 
do k=1,nzm
 qpsrc(k)=0.
 qpevp(k)=0.
 do j=1,ny
  do i=1,nx	  
   dqpsrc = 0.0
   dqpevp = 0.0
	  
!-------     Autoconversion/accretion 

   if(qn(i,j,k)+qp(i,j,k).gt.0.) then


         omn = max(0.,min(1.,(tabs(i,j,k)-tbgmin)*a_bg))
         omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
         omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))

!	 if(qn(i,j,k).gt.0.) then
         if(cld3d(i,j,k).gt.0.) then  ! the generation of precipitating condensate
     
           qcc = qn(i,j,k) * omn /cld3d_temp(i,j,k)
           qii = qn(i,j,k) * (1.-omn)/cld3d_temp(i,j,k)

           if(qcc .gt. qcw0) then
            autor = alphaelq
           else
            autor = 0.
           endif 

           if(qii .gt. qci0) then
            autos = betaelq*coefice(k)
           else
            autos = 0.
           endif 

           accrr = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp / cldmax(i,j,k)
             accrr = accrrc(k) * qrr ** powr1
           end if
           accrcs = 0.
           accris = 0. 
           if(omp.lt.0.999.and.omg.lt.0.999) then
             qss = qp(i,j,k) * (1.-omp)*(1.-omg) / cldmax(i,j,k)
             tmp = qss ** pows1
             accrcs = accrsc(k) * tmp
             accris = accrsi(k) * tmp 
           end if
           accrcg = 0.
           accrig = 0. 
           if(omp.lt.0.999.and.omg.gt.0.001) then
             qgg = qp(i,j,k) * (1.-omp)*omg / cldmax(i,j,k)
             tmp = qgg ** powg1
             accrcg = accrgc(k) * tmp
             accrig = accrgi(k) * tmp 
           endif
           qcc = (qcc+dtn*autor*qcw0)/(1.+dtn*(accrr+accrcs+accrcg+autor))
           qii = (qii+dtn*autos*qci0)/(1.+dtn*(accris+accrig+autos))
           dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+ &
             (accris+accrig)*qii + (accrcs+accrcg)*qcc + autos*(qii-qci0))

           dq = dq * cld3d(i,j,k)  ! convert fro the in-cloud value to grid-mean value

           dq = min(dq,qn(i,j,k))
!           qp(i,j,k) = qp(i,j,k) + dq
!           q(i,j,k) = q(i,j,k) - dq
!           qn(i,j,k) = qn(i,j,k) - dq
           dqpsrc = dq
	   qpsrc(k) = qpsrc(k) + dq

         end if

         !elseif(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then
         ! Evaporation is only allowed when cldmax exceeds cld3d_temp
!         if(qp(i,j,k).gt.qp_threshold.and.cldmax(i,j,k).gt.cld3d_temp(i,j,k)) then 
        if(qp(i,j,k).gt.qp_threshold.and.qn(i,j,k).eq.0.) then

           qsatt = 0.
           if(omn.gt.0.001) qsatt = qsatt + omn*qsatw_crm(tabs(i,j,k),pres(k))
           if(omn.lt.0.999) qsatt = qsatt + (1.-omn)*qsati_crm(tabs(i,j,k),pres(k))
           dq = 0.
           if(omp.gt.0.001) then
             qrr = qp(i,j,k) * omp /cldmax(i,j,k)
             dq = dq + evapr1(k)*sqrt(qrr) + evapr2(k)*qrr**powr2 
           end if
           if(omp.lt.0.999.and.omg.lt.0.999) then
             qss = qp(i,j,k) * (1.-omp)*(1.-omg) / cldmax(i,j,k)
             dq = dq + evaps1(k)*sqrt(qss) + evaps2(k)*qss**pows2 
           end if
           if(omp.lt.0.999.and.omg.gt.0.001) then
             qgg = qp(i,j,k) * (1.-omp)*omg /cldmax(i,j,k)
             dq = dq + evapg1(k)*sqrt(qgg) + evapg2(k)*qgg**powg2
           end if

!           dq = dq * dtn * (q(i,j,k) /qsatt-1.) 
           qclr = max(0., (q(i,j,k)-qn(i,j,k)-qsatt * cld3d(i,j,k)))/max(0.001, (1-cld3d(i,j,k)))
           qclr = min(qclr, qsatt)
           dq = dq * dtn * (qclr/qsatt-1.)
           dq = dq * (cldmax(i,j,k) - cld3d_temp(i,j,k))  ! convert this to the grid-mean value 

           dq = max(-0.5*qp(i,j,k),dq) 
!           qp(i,j,k) = qp(i,j,k) + dq
!           q(i,j,k) = q(i,j,k) - dq
           dqpevp = dq
	   qpevp(k) = qpevp(k) + dq

         end if
	
         if(qp(i,j,k).le.qp_threshold .and. cld3d(i,j,k).le.0) then
!           q(i,j,k) = q(i,j,k) + qp(i,j,k)
           dqpevp = dqpevp - qp(i,j,k)
	   qpevp(k) = qpevp(k) - qp(i,j,k)
!           qp(i,j,k) = 0.
         endif

    endif

    qp(i,j,k) = qp(i,j,k) + dqpsrc + dqpevp
    q(i,j,k) = q(i,j,k) - dqpsrc - dqpevp
    qn(i,j,k) = qn(i,j,k) - dqpsrc

    dq = qp(i,j,k)
    qp(i,j,k)=max(0.,qp(i,j,k))
    q(i,j,k) = q(i,j,k) + (dq-qp(i,j,k))

  end do
 enddo
enddo

!call t_stopf ('precip_proc_clubb')

#endif /*CLUBB_CRM*/
end subroutine precip_proc_clubb
#endif

