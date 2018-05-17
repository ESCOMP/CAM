
subroutine advect2_mom_xy		
	
!        momentum tendency due to 2nd-order-central horizontal advection

use crmx_vars

implicit none
	
real fu(0:nx,1-YES3D:ny,nzm) 
real fv(0:nx,1-YES3D:ny,nzm)
real fw(0:nx,1-YES3D:ny,nzm)
real dx25, dy25, irho

integer i, j, k, kc, kcu, ic, jb, ib, jc

dx25 = 0.25 / dx
dy25 = 0.25 / dy


if(RUN3D) then

do k = 1,nzm	
 kc= k+1
 kcu =min(kc, nzm)
 irho = 1./(rhow(kc)*adzw(kc))
	 
 do j = 1, ny	
  jb = j-1
  do i = 0, nx 
   ic = i+1			
   fu(i,j,k)=dx25*(u(ic,j,k)+u(i,j,k))*(u(i,j,k)+u(ic,j,k))
   fv(i,j,k)=dx25*(u(ic,j,k)+u(ic,jb,k))*(v(i,j,k)+v(ic,j,k))
   fw(i,j,k)=dx25*(u(ic,j,k)*rho(k)*adz(k)+ &
             u(ic,j,kcu)*rho(kcu)*adz(kcu))*(w(i,j,kc)+w(ic,j,kc))
  end do 
  do i = 1, nx	  
   ib = i-1
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - (fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - (fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na) = dwdt(i,j,kc,na)-irho*(fw(i,j,k)-fw(ib,j,k))
  end do 
 end do 

 do j = 0, ny 
  jc = j+1	
  do i = 1, nx
   ib = i-1
   fu(i,j,k)=dy25*(v(i,jc,k)+v(ib,jc,k))*(u(i,j,k)+u(i,jc,k))
   fv(i,j,k)=dy25*(v(i,jc,k)+v(i,j,k))*(v(i,j,k)+v(i,jc,k))
   fw(i,j,k)=dy25*(v(i,jc,k)*rho(k)*adz(k)+ &
             v(i,jc,kcu)*rho(kcu)*adz(kcu))*(w(i,j,kc)+w(i,jc,kc))
  end do
 end do 
 do j = 1,ny	
  jb = j-1
  do i = 1, nx
   dudt(i,j,k,na) = dudt(i,j,k,na) - (fu(i,j,k) - fu(i,jb,k))
   dvdt(i,j,k,na) = dvdt(i,j,k,na) - (fv(i,j,k) - fv(i,jb,k))
   dwdt(i,j,kc,na)= dwdt(i,j,kc,na)-irho*(fw(i,j,k)-fw(i,jb,k))
  end do
 end do 
 
end do ! k


else

j=1

do k = 1,nzm	
 kc= k+1
 kcu =min(kc, nzm)
 irho = 1./(rhow(kc)*adzw(kc))	 

  do i = 0, nx 
   ic = i+1			
   fu(i,j,k)=dx25*(u(ic,j,k)+u(i,j,k))*(u(i,j,k)+u(ic,j,k))
   fv(i,j,k)=dx25*(u(ic,j,k)+u(i,j,k))*(v(i,j,k)+v(ic,j,k))
   fw(i,j,k)=dx25*(u(ic,j,k)*rho(k)*adz(k)+ &
             u(ic,j,kcu)*rho(kcu)*adz(kcu))*(w(i,j,kc)+w(ic,j,kc))
  end do 
  do i = 1, nx	  
   ib = i-1
    dudt(i,j,k,na)  = dudt(i,j,k,na)  - (fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)  = dvdt(i,j,k,na)  - (fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na) = dwdt(i,j,kc,na)-irho*(fw(i,j,k)-fw(ib,j,k))
  end do 

end do ! k

endif

end subroutine advect2_mom_xy

