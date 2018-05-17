subroutine uvw
	
! update the velocity field 

use crmx_vars
use crmx_params
implicit none
	
u(1:nx,1:ny,1:nzm) = dudt(1:nx,1:ny,1:nzm,nc)
v(1:nx,1:ny,1:nzm) = dvdt(1:nx,1:ny,1:nzm,nc)
w(1:nx,1:ny,1:nzm) = dwdt(1:nx,1:ny,1:nzm,nc)

end subroutine uvw
