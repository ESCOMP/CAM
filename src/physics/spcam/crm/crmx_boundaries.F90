
subroutine boundaries(flag)

use crmx_grid, only: dompi

        	
implicit none
integer flag

!call t_startf ('boundaries')

if(dompi) then
  call task_boundaries(flag)
else
  call periodic(flag)
end if

!call t_stopf ('boundaries')

end subroutine boundaries
