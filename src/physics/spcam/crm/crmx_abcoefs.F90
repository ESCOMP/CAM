
subroutine abcoefs

!      coefficients for the Adams-Bashforth scheme

use crmx_grid 

implicit none

real alpha, beta
	
if(nstep.ge.3.and.nadams.eq.3.or.nrestart.eq.2) then
  alpha = dt3(nb) / dt3(na)
  beta = dt3(nc) / dt3(na)
  ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta)
  bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha)
  at = 1. - bt - ct
else if(nstep.ge.2) then
  at = 3./2.
  bt = -1./2.
  ct = 0.
else
  at = 1.
  bt = 0.
  ct = 0.
end if

end subroutine abcoefs
