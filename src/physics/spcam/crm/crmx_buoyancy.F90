
subroutine buoyancy()

use crmx_vars
use crmx_params
implicit none
	
integer i,j,k,kb
real betu, betd

if(docolumn) return

do k=2,nzm	
 kb=k-1
 betu=adz(kb)/(adz(k)+adz(kb))
 betd=adz(k)/(adz(k)+adz(kb))
 do j=1,ny
  do i=1,nx

   dwdt(i,j,k,na)=dwdt(i,j,k,na) +  &
      bet(k)*betu* &
     ( tabs0(k)*(epsv*(qv(i,j,k)-qv0(k))-(qcl(i,j,k)+qci(i,j,k)-qn0(k)+qpl(i,j,k)+qpi(i,j,k)-qp0(k))) &
       +(tabs(i,j,k)-tabs0(k))*(1.+epsv*qv0(k)-qn0(k)-qp0(k)) ) &
    + bet(kb)*betd* &
     ( tabs0(kb)*(epsv*(qv(i,j,kb)-qv0(kb))-(qcl(i,j,kb)+qci(i,j,kb)-qn0(kb)+qpl(i,j,kb)+qpi(i,j,kb)-qp0(kb))) &
       +(tabs(i,j,kb)-tabs0(kb))*(1.+epsv*qv0(kb)-qn0(kb)-qp0(kb)) )  

  end do ! i
 end do ! j
end do ! k

end subroutine buoyancy


