subroutine diffuse_mom

!  Interface to the diffusion routines

use crmx_vars
implicit none
integer i,j,k

!call t_startf ('diffuse_mom')

if(RUN3D) then
   call diffuse_mom3D()
else
   call diffuse_mom2D()
endif

!call t_stopf ('diffuse_mom')

end subroutine diffuse_mom

