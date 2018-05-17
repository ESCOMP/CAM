subroutine bound_exchange(f,dimx1,dimx2,dimy1,dimy2,dimz,i_1, i_2, j_1, j_2, id)
	
! periodic boundary exchange


use crmx_grid
implicit none
	
integer dimx1, dimx2, dimy1, dimy2, dimz
integer i_1, i_2, j_1, j_2
real f(dimx1:dimx2, dimy1:dimy2, dimz)
integer id   ! id of the sent field (dummy variable)
	
real buffer((nx+ny)*3*nz)	! buffer for sending data
	
integer i, j, k, n
integer i1, i2, j1, j2
	
i1 = i_1 - 1
i2 = i_2 - 1
j1 = j_1 - 1
j2 = j_2 - 1

!----------------------------------------------------------------------
!  Send buffers to neighbors
!----------------------------------------------------------------------


	if(RUN3D) then

! "North" -> "South":	

	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

! "North-East" -> "South-West":	
	
	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

! "South-East" -> "North-West":

	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

! "South" -> "North":

	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=1,nx
	           n = n+1
	           buffer(n) = f(i,j,k) 
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=1,nx
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

! "South-West" -> "North-East":
	  
	     n=0
	     do k=1,dimz
	       do j=1,1+j2
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k) 
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do


! To "North-West" -> "South-East":
	  	  
	     n=0
	     do k=1,dimz
	       do j=ny-j1,ny
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=-j1,0
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do
	     

	endif

!  "East" -> "West":
 	  
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=nx-i1,nx
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=-i1,0
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do

! "West" -> "East":

	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=1,1+i2
	           n = n+1
	           buffer(n) = f(i,j,k)
	         end do
	       end do
	     end do
	     n=0
	     do k=1,dimz
	       do j=1,ny
	         do i=nxp1,nxp1+i2
	           n = n+1
	           f(i,j,k) = buffer(n)
	         end do
	       end do
	     end do


end subroutine bound_exchange
	     
	     
