
subroutine advect_scalar (f,fadv,flux,f2leadv,f2legrad,fwleadv,doit)
 	
!     positively definite monotonic advection with non-oscillatory option

use crmx_grid
use crmx_vars, only: u, v, w, rho, rhow
use crmx_params, only: docolumn

implicit none

real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real flux(nz), fadv(nz)
real f2leadv(nz),f2legrad(nz),fwleadv(nz)
logical doit

real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
integer i,j,k

if(docolumn) then
  flux = 0.
  return
end if

!call t_startf ('advect_scalars')
	
 df(:,:,:) = f(:,:,:)

if(RUN3D) then
  call advect_scalar3D(f, u, v, w, rho, rhow, flux)
else
  call advect_scalar2D(f, u, w, rho, rhow, flux)	  
endif

  do k=1,nzm
    fadv(k)=0.
    do j=1,ny
     do i=1,nx
      fadv(k)=fadv(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

!call t_stopf ('advect_scalars')

end subroutine advect_scalar

