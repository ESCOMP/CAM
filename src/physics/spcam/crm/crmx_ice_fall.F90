
subroutine ice_fall()


! Sedimentation of ice:

use crmx_vars
use crmx_microphysics, only: micro_field, index_cloud_ice
!use micro_params
use crmx_params

implicit none

integer i,j,k, kb, kc, kmax, kmin, ici
real coef,dqi,lat_heat,vt_ice
real omnu, omnc, omnd, qiu, qic, qid, tmp_theta, tmp_phi
real fz(nx,ny,nz)

kmax=0
kmin=nzm+1

do k = 1,nzm
 do j = 1, ny
  do i = 1, nx
      if(qcl(i,j,k)+qci(i,j,k).gt.0..and. tabs(i,j,k).lt.273.15) then
        kmin = min(kmin,k)
        kmax = max(kmax,k)
      end if
  end do
 end do
end do

do k = 1,nzm
   qifall(k) = 0.
   tlatqi(k) = 0.
end do

if(index_cloud_ice.eq.-1) return

!call t_startf ('ice_fall')

fz = 0.

! Compute cloud ice flux (using flux limited advection scheme, as in
! chapter 6 of Finite Volume Methods for Hyperbolic Problems by R.J.
! LeVeque, Cambridge University Press, 2002). 
do k = max(1,kmin-1),kmax
   ! Set up indices for x-y planes above and below current plane.
   kc = min(nzm,k+1)
   kb = max(1,k-1)
   ! CFL number based on grid spacing interpolated to interface i,j,k-1/2
   coef = dtn/(0.5*(adz(kb)+adz(k))*dz)
   do j = 1,ny
      do i = 1,nx
         ! Compute cloud ice density in this cell and the ones above/below.
         ! Since cloud ice is falling, the above cell is u (upwind),
         ! this cell is c (center) and the one below is d (downwind). 

         qiu = rho(kc)*qci(i,j,kc)
         qic = rho(k) *qci(i,j,k) 
         qid = rho(kb)*qci(i,j,kb) 

         ! Ice sedimentation velocity depends on ice content. The fiting is
         ! based on the data by Heymsfield (JAS,2003). -Marat
         vt_ice = min(0.4,8.66*(max(0.,qic)+1.e-10)**0.24)   ! Heymsfield (JAS, 2003, p.2607)

         ! Use MC flux limiter in computation of flux correction.
         ! (MC = monotonized centered difference).
!         if (qic.eq.qid) then
         if (abs(qic-qid).lt.1.0e-25) then  ! when qic, and qid is very small, qic_qid can still be zero
                                            ! even if qic is not equal to qid. so add a fix here +++mhwang
            tmp_phi = 0.
         else
            tmp_theta = (qiu-qic)/(qic-qid)
            tmp_phi = max(0.,min(0.5*(1.+tmp_theta),2.,2.*tmp_theta))
         end if

         ! Compute limited flux.
         ! Since falling cloud ice is a 1D advection problem, this
         ! flux-limited advection scheme is monotonic.
         fz(i,j,k) = -vt_ice*(qic - 0.5*(1.-coef*vt_ice)*tmp_phi*(qic-qid))
      end do
   end do
end do
fz(:,:,nz) = 0.

ici = index_cloud_ice

do k=max(1,kmin-2),kmax
   coef=dtn/(dz*adz(k)*rho(k))
   do j=1,ny
      do i=1,nx
         ! The cloud ice increment is the difference of the fluxes.
         dqi=coef*(fz(i,j,k)-fz(i,j,k+1))
         ! Add this increment to both non-precipitating and total water.
         micro_field(i,j,k,ici)  = micro_field(i,j,k,ici)  + dqi
         ! Include this effect in the total moisture budget.
         qifall(k) = qifall(k) + dqi

         ! The latent heat flux induced by the falling cloud ice enters
         ! the liquid-ice static energy budget in the same way as the
         ! precipitation.  Note: use latent heat of sublimation. 
         lat_heat  = (fac_cond+fac_fus)*dqi
         ! Add divergence of latent heat flux to liquid-ice static energy.
         t(i,j,k)  = t(i,j,k)  - lat_heat
         ! Add divergence to liquid-ice static energy budget.
         tlatqi(k) = tlatqi(k) - lat_heat
      end do
   end do
end do

coef=dtn/dz
do j=1,ny
   do i=1,nx
      dqi=-coef*fz(i,j,1)
      precsfc(i,j) = precsfc(i,j)+dqi
      precssfc(i,j) = precssfc(i,j)+dqi
   end do
end do

!call t_stopf ('ice_fall')

end subroutine ice_fall

