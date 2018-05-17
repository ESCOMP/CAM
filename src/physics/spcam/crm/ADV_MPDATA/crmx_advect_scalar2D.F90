
subroutine advect_scalar2D (f, u, w, rho, rhow, flux)
 	
!     positively definite monotonic advection with non-oscillatory option

use crmx_grid
use crmx_params, only: dowallx
implicit none


real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
real rho(nzm)
real rhow(nz)
real flux(nz)
	
real mx (0:nxp1,1,nzm)
real mn (0:nxp1,1,nzm)
real uuu(-1:nxp3,1,nzm)
real www(-1:nxp2,1,nz)

real eps, dd
integer i,j,k,ic,ib,kc,kb
logical nonos
real iadz(nzm),irho(nzm),irhow(nzm)

real x1, x2, a, b, a1, a2, y
real andiff,across,pp,pn
andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
across(x1,a1,a2)=0.03125*a1*a2*x1
pp(y)= max(0.,y)
pn(y)=-min(0.,y)
	
nonos = .true.
eps = 1.e-10

j=1

www(:,:,nz)=0.

if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
       do i=dimx1_u,1
         u(i,j,k) = 0.
       end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
       do i=nx+1,dimx2_u
         u(i,j,k) = 0.
       end do
    end do
  end if

end if

!-----------------------------------------
	 	 
if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
    mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
   end do
 end do

end if  ! nonos

do k=1,nzm
  kb=max(1,k-1)
  do i=-1,nxp3
    uuu(i,j,k)=max(0.,u(i,j,k))*f(i-1,j,k)+min(0.,u(i,j,k))*f(i,j,k)
  end do
  do i=-1,nxp2
    www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
  end do
  flux(k) = 0.
  do i=1,nx
    flux(k) = flux(k) + www(i,j,k)	
  end do
end do

do k=1,nzm
  irho(k) = 1./rho(k)
  iadz(k) = 1./adz(k)
   do i=-1,nxp2
      f(i,j,k) = f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) & 
                        + (www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)            
   end do
end do 


do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  dd=2./(kc-kb)/adz(k)
  irhow(k)=1./(rhow(k)*adz(k))
  do i=0,nxp2
   ib=i-1
   uuu(i,j,k)=andiff(f(ib,j,k),f(i,j,k),u(i,j,k),irho(k)) &
      - across(dd*(f(ib,j,kc)+f(i,j,kc)-f(ib,j,kb)-f(i,j,kb)), &
              u(i,j,k), w(ib,j,k)+w(ib,j,kc)+w(i,j,k)+w(i,j,kc)) *irho(k)
  end do
          

  do i=0,nxp1
   ib=i-1
   ic=i+1
   www(i,j,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),irhow(k)) &
      -across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), &
        w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) *irho(k)
  end do
end do
www(:,:,1) = 0.
!---------- non-osscilatory option ---------------

if(nonos) then

 do k=1,nzm
   kc=min(nzm,k+1)
   kb=max(1,k-1)
   do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mx(i,j,k))
    mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mn(i,j,k))
   end do
 end do

 do k=1,nzm
   kc=min(nzm,k+1)
   do i=0,nxp1
    ic=i+1
     mx(i,j,k)=rho(k)*(mx(i,j,k)-f(i,j,k))/(pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+&
               iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)	
     mn(i,j,k)=rho(k)*(f(i,j,k)-mn(i,j,k))/(pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+&
               iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)	
   end do
 end do

 do k=1,nzm
  kb=max(1,k-1)
   do i=1,nxp1
    ib=i-1
    uuu(i,j,k)= pp(uuu(i,j,k))*min(1.,mx(i,j,k), mn(ib,j,k)) &
              - pn(uuu(i,j,k))*min(1.,mx(ib,j,k),mn(i,j,k))
   end do
   do i=1,nx
    www(i,j,k)= pp(www(i,j,k))*min(1.,mx(i,j,k), mn(i,j,kb)) &
              - pn(www(i,j,k))*min(1.,mx(i,j,kb),mn(i,j,k))
    flux(k) = flux(k) + www(i,j,k)	
   end do
 end do


endif ! nonos


 do k=1,nzm
  kc=k+1
   do i=1,nx
 ! MK: added fix for very small negative values (relative to positive values) 
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values is
 !     most likely truncation error.
    f(i,j,k)= max(0., f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) &
                     +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k))
   end do
 end do 

end subroutine advect_scalar2D


