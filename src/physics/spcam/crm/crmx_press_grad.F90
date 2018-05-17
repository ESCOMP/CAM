
subroutine press_grad
	
!       pressure term of the momentum equations

use crmx_vars
use crmx_params, only: dowallx, dowally
implicit none
	
real *8 rdx,rdy,rdz
integer i,j,k,kb,jb,ib

rdx=1./dx
rdy=1./dy

do k=1,nzm
 kb=max(1,k-1)
 rdz = 1./(dz*adzw(k))
 do j=1,ny
  jb=j-YES3D
  do i=1,nx
   ib=i-1 
   dudt(i,j,k,na)=dudt(i,j,k,na)-(p(i,j,k)-p(ib,j,k))*rdx	
   dvdt(i,j,k,na)=dvdt(i,j,k,na)-(p(i,j,k)-p(i,jb,k))*rdy	
   dwdt(i,j,k,na)=dwdt(i,j,k,na)-(p(i,j,k)-p(i,j,kb))*rdz	
  end do ! i
 end do ! j	
end do ! k

do k=1,nzm
 do j=1-YES3D,ny !bloss: 0,n* fixes computation of dp/d* in stats.               
  do i=0,nx                                                                      
    p(i,j,k)=p(i,j,k)*rho(k)  ! convert p'/rho to p'
  end do
 end do 
end do  

if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then

    do k=1,nzm
     do j=1,ny
      dudt(1,j,k,na) = 0.
     end do
    end do

end if

if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

    do k=1,nzm
     do i=1,nx
      dvdt(i,1,k,na) = 0.
     end do
    end do

end if

if(dompi) then
   call task_bound_duvdt()
else
   call bound_duvdt()	   
endif

call task_barrier()
	
end subroutine press_grad



