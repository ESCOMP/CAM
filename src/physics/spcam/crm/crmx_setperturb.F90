
subroutine setperturb(iseed)

!  Random noise
!  This surboutine has been updated for SPCAM5 (Minghuai.Wang@pnnl.gov, April, 2012). 
!  Now the random generator is seeded based on the global column id, which gets rid
!  of the dependence of the SPCAM reulst on pcols. 

use crmx_vars
use crmx_sgs, only: setperturb_sgs

implicit none

integer, intent(in) :: iseed

integer i,j,k
real rrr,ranf_
integer, allocatable :: rndm_seed(:)
integer :: rndm_seed_sz
real :: t02(nzm)
real :: tke02(nzm)

!call ranset_(30*rank)
call random_seed(size=rndm_seed_sz)
allocate(rndm_seed(rndm_seed_sz))

rndm_seed = iseed
call random_seed(put=rndm_seed)

call setperturb_sgs(0)  ! set sgs fields

t02 = 0.0
tke02 = 0.0
do k=1,nzm
 do j=1,ny
  do i=1,nx
    rrr=1.-2.*ranf_()

    if(k.le.5) then
      t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
    endif
    t02(k) = t02(k) + t(i,j,k)/(nx*ny)
  end do
 end do

! energy conservation +++mhwang (2012-06)
 do j=1, ny
  do i=1, nx
    if(k.le.5) then
      t(i,j,k) = t(i,j,k) * t0(k)/t02(k)
    end if 
  end do
 end do
end do

deallocate(rndm_seed)

end

