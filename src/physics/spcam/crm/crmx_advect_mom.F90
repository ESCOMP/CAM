subroutine advect_mom

use crmx_vars
use crmx_params, only: docolumn

implicit none
integer i,j,k

if(docolumn) return

!call t_startf ('advect_mom')

call advect2_mom_xy()
call advect2_mom_z()

!call t_stopf ('advect_mom')

end subroutine advect_mom

