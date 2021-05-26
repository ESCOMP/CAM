module edyn_mudcom
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: cor2
  public :: factri
  public :: factrp
  public :: swk2
  public :: trsfc2
  public :: prolon2
  public :: res2
  public :: sgfa
  public :: sgsl
  public :: transp

!-----------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------
!
!     file mudcom.f
!  .                                                             .
!  .                      MUDPACK version 4.0                    .
!
! ... author and specialist
!
!          John C. Adams (National Center for Atmospheric Research) (retired)

! ... For MUDPACK information, visit the website:
!     (https://www2.cisl.ucar.edu/resources/legacy/mudpack)
!
! ... purpose
!
!     mudcom.f is a common subroutines file containing subroutines
!     called by some or all of the real two- and three-dimensional
!     mudpack solvers.  mudcom.f must be loaded with any real mudpack
!     solver.
!
!     cb mud2cr1: call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
!
!-----------------------------------------------------------------------
      subroutine swk2(nfx,nfy,phif,rhsf,phi,rhs)
!
!     set phif,rhsf input in arrays which include
!     virtual boundaries for phi (for all 2-d real codes)
!
      integer nfx,nfy,i,j
      real(r8) :: phif(nfx,nfy),rhsf(nfx,nfy)
      real(r8) :: phi(0:nfx+1,0:nfy+1),rhs(nfx,nfy)
      do j=1,nfy
        do i=1,nfx
          phi(i,j) = phif(i,j)
          rhs(i,j) = rhsf(i,j)
        end do
      end do
!
!     set virtual boundaries in phi to zero
!
      do j=0,nfy+1
        phi(0,j) = 0.0_r8
        phi(nfx+1,j) = 0.0_r8
      end do
      do i=0,nfx+1
        phi(i,0) = 0.0_r8
        phi(i,nfy+1) = 0.0_r8
      end do
      return
      end subroutine swk2
!-----------------------------------------------------------------------
      subroutine trsfc2(nx,ny,phi,rhs,ncx,ncy,phic,rhsc)
!
!     transfer fine grid to coarse grid
!
      integer nx,ny,ncx,ncy,i,j,ic,jc
      real(r8) :: phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real(r8) :: phic(0:ncx+1,0:ncy+1),rhsc(ncx,ncy)
!
!     set virtual boundaries in phic to zero
!
      do jc=0,ncy+1
        phic(0,jc) = 0.0_r8
        phic(ncx+1,jc) = 0.0_r8
      end do
      do ic=0,ncx+1
        phic(ic,0) = 0.0_r8
        phic(ic,ncy+1) = 0.0_r8
      end do
      if (ncx.lt.nx .and. ncy.lt.ny) then
!
!     coarsening in both x and y
!
        do jc=1,ncy
          j = jc+jc-1
          do ic=1,ncx
            i = ic+ic-1
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      else if (ncx.lt.nx .and. ncy.eq.ny) then
!
!     coarsening in x only
!
        do jc=1,ncy
          j = jc
          do ic=1,ncx
            i = ic+ic-1
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      else
!
!     coarsening in y only
!
        do jc=1,ncy
          j = jc+jc-1
          do ic=1,ncx
            i = ic
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      end if
      return
      end subroutine trsfc2
!-----------------------------------------------------------------------
      subroutine res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)

      integer nx,ny,ncx,ncy,nxa,nxb,nyc,nyd
      integer i,j,ic,jc,im1,ip1,jm1,jp1,ix,jy
!
!     restrict fine grid residual in resf to coarse grid in rhsc
!     using full weighting for all real 2d codes
!
      real(r8) :: resf(nx,ny),rhsc(ncx,ncy)
!
!     set x,y coarsening integer subscript scales
!
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
!
!     restrict on interior
!
      if (ncy.lt.ny .and. ncx.lt.nx) then
!
!     coarsening in both directions
!
        do jc=2,ncy-1
          j = jc+jc-1
          do ic=2,ncx-1
            i = ic+ic-1
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+ &
                           resf(i+1,j+1)+2._r8*(resf(i-1,j)+resf(i+1,j)+ &
                           resf(i,j-1)+resf(i,j+1))+4._r8*resf(i,j))*.0625_r8
          end do
        end do
      else if (ncy.eq.ny) then
!
!     no coarsening in y but coarsening in x
!
        do jc=2,ncy-1
          j = jc
          do ic=2,ncx-1
            i = ic+ic-1
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+ &
                           resf(i+1,j+1)+2._r8*(resf(i-1,j)+resf(i+1,j)+ &
                           resf(i,j-1)+resf(i,j+1))+4._r8*resf(i,j))*.0625_r8
          end do
        end do
      else
!
!     no coarsening in x but coarsening in y
!
        do jc=2,ncy-1
          j = jc+jc-1
          do ic=2,ncx-1
            i = ic
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+ &
                           resf(i+1,j+1)+2._r8*(resf(i-1,j)+resf(i+1,j)+ &
                           resf(i,j-1)+resf(i,j+1))+4._r8*resf(i,j))*.0625_r8
          end do
        end do
      end if
!
!     set residual on boundaries
!
      do jc=1,ncy,ncy-1
!
!     y=yc,yd boundaries
!
        j = jc+jy*(jc-1)
        jm1 = max0(j-1,2)
        jp1 = min0(j+1,ny-1)
        if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
        if (j.eq.ny .and. nyc.eq.0) jp1 = 2
!
!     y=yc,yd and x=xa,xb cornors
!
        do ic=1,ncx,ncx-1
          i = ic+ix*(ic-1)
          im1 = max0(i-1,2)
          ip1 = min0(i+1,nx-1)
          if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
          if (i.eq.nx .and. nxa.eq.0) ip1 = 2
          rhsc(ic,jc) = (resf(im1,jm1)+resf(ip1,jm1)+resf(im1,jp1)+ &
                         resf(ip1,jp1)+2._r8*(resf(im1,j)+resf(ip1,j)+ &
                         resf(i,jm1)+resf(i,jp1))+4._r8*resf(i,j))*.0625_r8
        end do
!
!     set y=yc,yd interior edges
!
        do ic=2,ncx-1
          i = ic+ix*(ic-1)
          rhsc(ic,jc) = (resf(i-1,jm1)+resf(i+1,jm1)+resf(i-1,jp1)+ &
                         resf(i+1,jp1)+2._r8*(resf(i-1,j)+resf(i+1,j)+ &
                         resf(i,jm1)+resf(i,jp1))+4._r8*resf(i,j))*.0625_r8
        end do
      end do
!
!     set x=xa,xb interior edges
!
      do ic=1,ncx,ncx-1
        i = ic+ix*(ic-1)
        im1 = max0(i-1,2)
        ip1 = min0(i+1,nx-1)
        if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
        if (i.eq.nx .and. nxa.eq.0) ip1 = 2
        do jc=2,ncy-1
          j = jc+jy*(jc-1)
          rhsc(ic,jc) = (resf(im1,j-1)+resf(ip1,j-1)+resf(im1,j+1)+ &
                         resf(ip1,j+1)+2._r8*(resf(im1,j)+resf(ip1,j)+ &
                         resf(i,j-1)+resf(i,j+1))+4._r8*resf(i,j))*.0625_r8
        end do
      end do
!
!     set coarse grid residual zero on specified boundaries
!
      if (nxa.eq.1) then
        do jc=1,ncy
          rhsc(1,jc) = 0.0_r8
        end do
      end if
      if (nxb.eq.1) then
        do jc=1,ncy
          rhsc(ncx,jc) = 0.0_r8
        end do
      end if
      if (nyc.eq.1) then
        do ic=1,ncx
          rhsc(ic,1) = 0.0_r8
        end do
      end if
      if (nyd.eq.1) then
        do ic=1,ncx
          rhsc(ic,ncy) = 0.0_r8
        end do
      end if
      return
      end subroutine res2
!-----------------------------------------------------------------------
!
!     prolon2 modified from rgrd2u 11/20/97
!
      subroutine prolon2(ncx,ncy,p,nx,ny,q,nxa,nxb,nyc,nyd,intpol)

      integer ncx,ncy,nx,ny,intpol,nxa,nxb,nyc,nyd
      real(r8) :: p(0:ncx+1,0:ncy+1),q(0:nx+1,0:ny+1)
      integer i,j,jc,ist,ifn,jst,jfn,joddst,joddfn
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      joddst = 1
      joddfn = ny
      if (nxa.eq.1) then
        ist = 2
      end if
      if (nxb.eq.1) then
        ifn = nx-1
      end if
      if (nyc.eq.1) then
        jst = 2
        joddst = 3
      end if
      if (nyd.eq.1) then
        jfn = ny-1
        joddfn = ny-2
      end if
      if (intpol.eq.1 .or. ncy.lt.4) then
!
!     linearly interpolate in y
!
        if (ncy .lt. ny) then
!
!     ncy grid is an every other point subset of ny grid
!     set odd j lines interpolating in x and then set even
!     j lines by averaging odd j lines
!
          do j=joddst,joddfn,2
            jc = j/2+1
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
          do j=2,jfn,2
            do i=ist,ifn
              q(i,j) = 0.5_r8*(q(i,j-1)+q(i,j+1))
            end do
          end do
!
!     set periodic virtual boundaries if necessary
!
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        else
!
!     ncy grid is equals ny grid so interpolate in x only
!
          do j=jst,jfn
            jc = j
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
!
!     set periodic virtual boundaries if necessary
!
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        end if
      else
!
!     cubically interpolate in y
!
        if (ncy .lt. ny) then
!
!     set every other point of ny grid by interpolating in x
!
          do j=joddst,joddfn,2
            jc = j/2+1
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
!
!     set deep interior of ny grid using values just
!     generated and symmetric cubic interpolation in y
!
          do j=4,ny-3,2
            do i=ist,ifn
            q(i,j)=(-q(i,j-3)+9._r8*(q(i,j-1)+q(i,j+1))-q(i,j+3))*.0625_r8
            end do
          end do
!
!     interpolate from q at j=2 and j=ny-1
!
          if (nyc.ne.0) then
!
!     asymmetric formula near nonperiodic y boundaries
!
            do i=ist,ifn
              q(i,2)=(5._r8*q(i,1)+15._r8*q(i,3)-5._r8*q(i,5)+q(i,7))*.0625_r8
              q(i,ny-1)=(5._r8*q(i,ny)+15._r8*q(i,ny-2)-5._r8*q(i,ny-4)+ &
                          q(i,ny-6))*.0625_r8
            end do
          else
!
!     periodicity in y alows symmetric formula near bndys
!
            do i=ist,ifn
              q(i,2) = (-q(i,ny-2)+9._r8*(q(i,1)+q(i,3))-q(i,5))*.0625_r8
              q(i,ny-1)=(-q(i,ny-4)+9._r8*(q(i,ny-2)+q(i,ny))-q(i,3))*.0625_r8
              q(i,ny+1) = q(i,2)
              q(i,0) = q(i,ny-1)
            end do
          end if
          return
        else
!
!     ncy grid is equals ny grid so interpolate in x only
!
          do j=jst,jfn
            jc = j
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
!
!     set periodic virtual boundaries if necessary
!
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        end if
      end if
      end subroutine prolon2
!-----------------------------------------------------------------------
!
!     11/20/97  modification of rgrd1u.f for mudpack
!
      subroutine prolon1(ncx,p,nx,q,nxa,nxb,intpol)

      integer intpol,nxa,nxb,ncx,nx,i,ic,ist,ifn,ioddst,ioddfn
      real(r8) :: p(0:ncx+1),q(0:nx+1)
      ist = 1
      ioddst = 1
      ifn = nx
      ioddfn = nx
      if (nxa.eq.1) then
        ist = 2
        ioddst = 3
      end if
      if (nxb.eq.1) then
        ifn = nx-1
        ioddfn = nx-2
      end if
      if (intpol.eq.1 .or. ncx.lt.4) then
!
!     linear interpolation in x
!
        if (ncx .lt. nx) then
!
!     every other point of nx grid is ncx grid
!
          do i=ioddst,ioddfn,2
            ic = (i+1)/2
            q(i) = p(ic)
          end do
          do i=2,ifn,2
            q(i) = 0.5_r8*(q(i-1)+q(i+1))
          end do
        else
!
!     nx grid equals ncx grid
!
          do i=ist,ifn
            q(i) = p(i)
          end do
        end if
!
!     set virtual end points if periodic
!
        if (nxa.eq.0) then
          q(0) = q(nx-1)
          q(nx+1) = q(2)
        end if
        return
      else
!
!     cubic interpolation in x
!
        if (ncx .lt. nx) then
          do i=ioddst,ioddfn,2
            ic = (i+1)/2
            q(i) = p(ic)
          end do
!
!      set deep interior with symmetric formula
!
          do i=4,nx-3,2
            q(i)=(-q(i-3)+9._r8*(q(i-1)+q(i+1))-q(i+3))*.0625_r8
          end do
!
!     interpolate from q at i=2 and i=nx-1
!
          if (nxa.ne.0) then
!
!     asymmetric formula near nonperiodic bndys
!
            q(2)=(5._r8*q(1)+15._r8*q(3)-5._r8*q(5)+q(7))*.0625_r8
            q(nx-1)=(5._r8*q(nx)+15._r8*q(nx-2)-5._r8*q(nx-4)+q(nx-6))*.0625_r8
          else
!
!     periodicity in x alows symmetric formula near bndys
!
            q(2) = (-q(nx-2)+9._r8*(q(1)+q(3))-q(5))*.0625_r8
            q(nx-1) = (-q(nx-4)+9._r8*(q(nx-2)+q(nx))-q(3))*.0625_r8
            q(nx+1) = q(2)
            q(0) = q(nx-1)
          end if
          return
        else
!
!     ncx grid equals nx grid
!
          do i=ist,ifn
            q(i) = p(i)
          end do
          if (nxa.eq.0) then
            q(0) = q(nx-1)
            q(nx+1) = q(2)
          end if
          return
        end if
      end if
      end subroutine prolon1
!-----------------------------------------------------------------------
      subroutine cor2(nx,ny,phif,ncx,ncy,phic,nxa,nxb,nyc,nyd,intpol,phcor)
!
!     add coarse grid correction in phic to fine grid approximation
!     in phif using linear or cubic interpolation
!
      integer i,j,nx,ny,ncx,ncy,nxa,nxb,nyc,nyd,intpol,ist,ifn,jst,jfn
      real(r8) :: phif(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real(r8) :: phcor(0:nx+1,0:ny+1)
      do j=0,ny+1
        do i=0,nx+1
          phcor(i,j) = 0.0_r8
        end do
      end do
!
!     lift correction in phic to fine grid in phcor
!
      call prolon2(ncx,ncy,phic,nx,ny,phcor,nxa,nxb,nyc,nyd,intpol)
!
!     add correction in phcor to phif on nonspecified boundaries
!
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      do j=jst,jfn
        do i=ist,ifn
          phif(i,j) = phif(i,j) + phcor(i,j)
        end do
      end do
!
!     add periodic points if necessary
!
      if (nyc.eq.0) then
        do i=ist,ifn
          phif(i,0) = phif(i,ny-1)
          phif(i,ny+1) = phif(i,2)
        end do
      end if
      if (nxa.eq.0) then
        do j=jst,jfn
          phif(0,j) = phif(nx-1,j)
          phif(nx+1,j) = phif(2,j)
        end do
      end if
      end subroutine cor2
!-----------------------------------------------------------------------
!
!     factri and factrip are:
!     subroutines called by any real mudpack solver which uses line
!     relaxation(s) within multigrid iteration.  these subroutines do
!     a vectorized factorization of m simultaneous tridiagonal systems
!     of order n arising from nonperiodic or periodic discretizations
!
      subroutine factri(m,n,a,b,c)
!
!     factor the m simultaneous tridiagonal systems of order n
!
      integer m,n,i,j
      real(r8) :: a(n,m),b(n,m),c(n,m)
      do i=2,n
        do j=1,m
          a(i-1,j) = a(i-1,j)/b(i-1,j)
          b(i,j) = b(i,j)-a(i-1,j)*c(i-1,j)
       end do
      end do
      return
      end subroutine factri
!-----------------------------------------------------------------------
      subroutine factrp(m,n,a,b,c,d,e,sum)
!
!     factor the m simultaneous "tridiagonal" systems of order n
!     from discretized periodic system (leave out periodic n point)
!     (so sweeps below only go from i=1,2,...,n-1) n > 3 is necessary
!
      integer m,n,i,j
      real(r8) :: a(n,m),b(n,m),c(n,m),d(n,m),e(n,m),sum(m)
      do j=1,m
        d(1,j) = a(1,j)
      end do
      do i=2,n-2
        do j=1,m
          a(i,j) = a(i,j)/b(i-1,j)
          b(i,j) = b(i,j)-a(i,j)*c(i-1,j)
          d(i,j) = -a(i,j)*d(i-1,j)
       end do
      end do
!
!     correct computation of last d element
!
      do j=1,m
        d(n-2,j) = c(n-2,j)+d(n-2,j)
      end do
      do j=1,m
        e(1,j) = c(n-1,j)/b(1,j)
      end do
      do i=2,n-3
        do j=1,m
          e(i,j) = -e(i-1,j)*c(i-1,j)/b(i,j)
        end do
      end do
      do j=1,m
        e(n-2,j) = (a(n-1,j)-e(n-3,j)*c(n-3,j))/b(n-2,j)
      end do
!
!     compute  inner product (e,d) for each j in sum(j)
!
      do j=1,m
        sum(j) = 0._r8
      end do
      do i=1,n-2
        do j=1,m
          sum(j) = sum(j)+e(i,j)*d(i,j)
        end do
      end do
!
!     set last diagonal element
!
      do j=1,m
        b(n-1,j) = b(n-1,j)-sum(j)
      end do
      return
      end subroutine factrp
!-----------------------------------------------------------------------
      subroutine transp(n,amat)
!
!     transpose n by n real matrix
!
      integer n,i,j
      real(r8) :: amat(n,n),temp
      do i=1,n-1
        do j=i+1,n
          temp = amat(i,j)
          amat(i,j) = amat(j,i)
          amat(j,i) = temp
        end do
      end do
      return
      end subroutine transp
!-----------------------------------------------------------------------
      subroutine sgfa (a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      real(r8) :: a(lda,*)
      real(r8) :: t
      integer :: j,k,kp1,l,nm1
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
         l = isfmax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
         if (a(l,k) .eq. 0.0e0_r8) go to 40
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
            t = -1.0e0_r8/a(k,k)
            call sscl(n-k,t,a(k+1,k),1)
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call sxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0_r8) info = n
      return
      end subroutine sgfa
!-----------------------------------------------------------------------
      subroutine sgsl (a,lda,n,ipvt,b,job)

      integer lda,n,ipvt(*),job
      real(r8) :: a(lda,*),b(*)
      real(r8) :: t
      integer k,kb,l,nm1
      nm1 = n - 1
      if (job .ne. 0) go to 50
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call sxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call sxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
         do 60 k = 1, n
            t = sdt(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdt(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end subroutine sgsl
!-----------------------------------------------------------------------
      function sdt(n,sx,incx,sy,incy) result(sdtx)

      real(r8), intent(in) :: sx(*),sy(*)
      integer,  intent(in) :: n, incx, incy

      integer :: i,ix,iy,m,mp1
      real(r8) :: sdtx
      real(r8) :: stemp

      stemp = 0.0e0_r8
      sdtx = 0.0e0_r8
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdtx = stemp
      return
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
         sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdtx = stemp
      return
      end function sdt
!-----------------------------------------------------------------------
      integer function isfmax(n,sx,incx)

      real(r8) :: sx(*),smax
      integer i,incx,ix,n
      isfmax = 0
      if( n .lt. 1 ) return
      isfmax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isfmax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isfmax = i
         smax = abs(sx(i))
   30 continue
      return
      end function isfmax
!-----------------------------------------------------------------------
      subroutine sxpy(n,sa,sx,incx,sy,incy)

      real(r8) :: sx(*),sy(*),sa
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      if (sa .eq. 0.0_r8) return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end subroutine sxpy
!-----------------------------------------------------------------------
      subroutine sscl(n,sa,sx,incx)

      real(r8) :: sa,sx(*)
      integer i,incx,m,mp1,n,nincx
      if(n.le.0)return
      if(incx.eq.1)go to 20
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
   10 continue
      return
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end subroutine sscl
!-----------------------------------------------------------------------
end module edyn_mudcom
