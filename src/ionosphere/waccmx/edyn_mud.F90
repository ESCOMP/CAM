!-----------------------------------------------------------------------
      subroutine mud(pe,jntl,isolve)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      use cam_abortutils   ,only: endrun
      use edyn_solve,only: nc,ncee,cee
!
      implicit none
      integer,intent(in) :: isolve
      integer jntl
!
!     set grid size params
!
      integer,parameter :: iixp = 5 , jjyq = 3, iiex = 5, jjey = 5
      integer,parameter :: nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1
!
!     estimate work space for point relaxation (see mud2cr.d)
!
      integer,parameter :: llwork=(7*(nnx+2)*(nny+2)+76*nnx*nny)/3
      real(r8) :: phi(nnx,nny),rhs(nnx,nny),work(llwork)
      real(r8) :: time0,time1
!
!     put integer and floating point argument names in contiguous
!     storage for labelling in vectors iprm,fprm
!
! btf 1/21/14: dimension iprm(17) to match iprm in edyn_muh2cr.F90
!     integer iprm(16),mgopt(4)
      integer iprm(17),mgopt(4)
      real(r8) :: fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                    iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                    iguess,maxcy,method,nwork,lwrkqd,itero
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real(r8) :: PE(NNX,1)
      integer maxcya
      DATA MAXCYA/150/
      integer mm,nn,jj,jjj
      real(r8) :: pi
!
!     set input integer arguments
!
      MM = NNX
      NN = NNY
      PI = 4._r8*ATAN(1._r8)
!
!     SET INPUT INTEGER PARAMETERS
!
      INTL = JNTL
!
!     set boundary condition flags
!
      nxa = 0
      nxb = 0
      nyc = 2
      nyd = 1
!
!     set grid sizes from parameter statements
!
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey
      nx = nnx
      ny = nny
!
!     set multigrid arguments (w(2,1) cycling with fully weighted
!     residual restriction and cubic prolongation)
!
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
!
!     set for one cycle
!
      maxcy = maxcya
!
!     set no initial guess forcing full multigrid cycling
!
      iguess = 0
!
!     set work space length approximation from parameter statement
!
      nwork = llwork
!
!     set line z relaxation
!
      method = 3
!
!     set end points of solution rectangle in (x,y) space
!
      xa = -pi
      xb =  pi
      yc = 0.0_r8
      yd = 0.5_r8*pi
!
!     set error control flag
!
      tolmax = 0.01_r8
!
!     set right hand side in rhs
!     initialize phi to zero
!
      if (isolve >= 0) then ! called from dynamo
        do i=1,nx
          do j=1,ny
            RHS(I,J) = CEE(I+(J-1)*NX+9*NX*NY)
            phi(i,j) = 0.0_r8
          end do
        end do
!
!     set specified boundaries in phi
!
        DO I=1,NX
          PHI(I,NY) = RHS(I,NY)/CEE(I+(NY-1)*NX+8*NX*NY)
        END DO
!
!     set specified boundaries in phi
!
      endif ! isolve
!
!     intialization call
!
      call mud2cr(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
      if (ierror.gt.0) call endrun('mud call init mud2cr')
!
!     attempt solution
!
      intl = 1
      call mud2cr(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
      if (ierror.gt.0) call endrun('mud call solve mud2cr')
!
!     COPY PHI TO PE
!
      DO J = 1,NY
        JJ = NY+J-1
        JJJ = NY+1-J
        DO I = 1,NX
          PE(I,JJ) = PHI(I,J)
          PE(I,JJJ) = PHI(I,J)
        END DO
      END DO
!     ITRANS = 0
!     CALL EZCNTR(PE(1,JMX0),IMX0,JMX0)
!     ITRANS = 1
!     CALL SET(.05,.95,.05,.95,-1.,1.,-1.,1.,1)
!     CALL CONREC(PE(1,JMX0),IMX0,IMX0,JMX0,0.,0.,0.,1,0,-1430B)
!     CALL FRAME
!     ITRANS = 0
!     CALL EZCNTR(PE(1,JMX0),IMX0,JMX0)
!     ITRANS = 1
      end subroutine mud
!-----------------------------------------------------------------------
!
!     file mud2cr.f  (version 4.0 modified for Cicley 2/99)
!                                                             .
!  .                      MUDPACK version 4.0                    .
!  .                                                             .
! ... author and specialist
!
!          John C. Adams (National Center for Atmospheric Research) (retired)
!
! ... For MUDPACK information, visit the website:
!     (https://www2.cisl.ucar.edu/resources/legacy/mudpack)
!
! ... purpose 
!
!     mud2cr attempts to produce a second order finite difference
!     approximation to the two dimensional nonseparable elliptic
!     partial differential equation with cross derivative
!
!       cxx(x,y)*pxx + cxy(x,y)*pxy + cyy(x,y)*pyy +
!
!       cx(x,y)*px + cy(x,y)*py + ce(x,y)*pe(x,y) = r(x,y)
!
! ... documentation
!
!     see the documentation on above website for a complete discussion
!     of how to use subroutine mud2cr.  
!
! ... required MUDPACK files
!
!     mudcom.f
!
!
!
      subroutine mud2cr(iparm,fparm,work,rhs,phi,mgopt, &
                        ierror,isolve)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      implicit none
      integer,intent(in) :: isolve
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      real(r8) :: fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      dimension iparm(17),fparm(6),mgopt(4)
      real(r8) :: work(*),phi(*),rhs(*)
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50), &
        nxk(50),nyk(50),isx,jsy
      data int / 0 /
      save int
!
      ierror = 1
      intl = iparm(1)    ! set and check intl on all calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
        int = 1
        if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
!
!     set  arguments internally
!     these will not be rechecked if intl=1!
!
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      ixp = iparm(6)
      jyq = iparm(7)
      iex = iparm(8)
      jey = iparm(9)
      ngrid = max0(iex,jey)
      nfx = iparm(10)
      nfy = iparm(11)
      iguess = iparm(12)
      maxcy = iparm(13)
      method = iparm(14)
      nwork = iparm(15)
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
!       set defaults
        kcycle = 2
        iprer = 2
        ipost = 1
        intpol = 3
      else
        iprer = mgopt(2)
        ipost = mgopt(3)
        intpol = mgopt(4)
      end if
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      tolmax = fparm(5)
      if (intl .eq. 0) then  ! intialization call
!
!     check input arguments
!
        ierror = 2   ! check boundary condition flags
        if (max0(nxa,nxb,nyc,nyd).gt.2) return
        if (min0(nxa,nxb,nyc,nyd).lt.0) return
        if (nxa.eq.0.and.nxb.ne.0) return
        if (nxa.ne.0.and.nxb.eq.0) return
        if (nyc.eq.0.and.nyd.ne.0) return
        if (nyc.ne.0.and.nyd.eq.0) return
        ierror = 3   ! check grid sizes
        if (ixp.lt.2) return
        if (jyq.lt.2) return
        ierror = 4
        ngrid = max0(iex,jey)
        if (iex.lt.1) return
        if (jey.lt.1) return
        if (ngrid.gt.50) return
        ierror = 5
        if (nfx.ne.ixp*2**(iex-1)+1) return
        if (nfy.ne.jyq*2**(jey-1)+1) return
        ierror = 6
        if (iguess*(iguess-1).ne.0) return
        ierror = 7
        if (maxcy.lt.1) return
        ierror = 8
        if (method.lt.0 .or. method.gt.3) return
        ierror = 9
!       compute and test minimum work space
        isx = 0
        if (method.eq.1 .or. method.eq.3) then
          if (nxa.ne.0) isx = 3
          if (nxa.eq.0) isx = 5
        end if
        jsy = 0
        if (method.eq.2 .or. method.eq.3) then
          if (nyc.ne.0) jsy = 3
          if (nyc.eq.0) jsy = 5
        end if
        kps = 1
        do k=1,ngrid
!       set subgrid sizes
          nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
          nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
          nx = nxk(k)
          ny = nyk(k)
          kps = kps+(nx+2)*(ny+2)+nx*ny*(10+isx+jsy)
        end do
        iparm(16) = kps+(nfx+2)*(nfy+2)   ! exact minimum work space
        lwork = iparm(16)
        if (lwork .gt. nwork) return
        ierror = 10   ! check solution region
        if (xb.le.xa .or. yd.le.yc) return
        ierror = 11
        if (tolmax .lt. 0.0_r8) return
        ierror = 12   ! multigrid parameters
        if (kcycle.lt.0) return
        if (min0(iprer,ipost).lt.1) return
        if ((intpol-1)*(intpol-3).ne.0) return
        if (max0(kcycle,iprer,ipost).gt.2) then
          ierror = -5   ! inefficient multigrid cycling
        end if
        if (ierror .gt. 0) ierror = 0   ! no fatal errors
!
!     set work space pointers and discretize pde at each grid level
!
        iw = 1
        do kb=1,ngrid
          k = ngrid-kb+1
          nx = nxk(k)
          ny = nyk(k)
          kpbgn(k) = iw
          kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
          ktxbgn(k) = kcbgn(k)+10*nx*ny
          ktybgn(k) = ktxbgn(k)+isx*nx*ny
          iw = ktybgn(k)+jsy*nx*ny
          ic = kcbgn(k)
          itx = ktxbgn(k)
          ity = ktybgn(k)
          klevel = k
          call dismd2cr(nx,ny,work(ic),work(itx),work(ity), &
                        work,ierror,isolve)
          end do
        return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call mud2cr1(nx,ny,rhs,phi,work)
      iparm(17) = itero
      if (tolmax.gt.0.0_r8) then   ! check for convergence
        fparm(6) = relmax
        if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end subroutine mud2cr
!-----------------------------------------------------------------------
      subroutine mud2cr1(nx,ny,rhsf,phif,wk)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      implicit none
      integer nx,ny
      real(r8) :: phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      real(r8) :: xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),&
      nxk(50),nyk(50),isx,jsy
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+9*nx*ny
!
!     set phif,rhsf in wk and adjust right hand side
!
      call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
!
!     no initial guess at finest grid level!
!
        do kb=2,ngrid
          k = ngrid-kb+1
          nx = nxk(k+1)
          ny = nyk(k+1)
          ip = kpbgn(k+1)
          ir = kcbgn(k+1)+9*nx*ny
          ncx = nxk(k)
          ncy = nyk(k)
          ipc = kpbgn(k)
          icc = kcbgn(k)
          irc = icc+9*ncx*ncy
!
!     transfer down to all grid levels
!
          call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,&
                      wk(ipc),wk(irc))
        end do
!
!     adjust right hand side at all grid levels in case
!     rhs or specified b.c. in phi or gbdy changed
!
        do k=1,ngrid
          nx = nxk(k)
          ny = nyk(k)
          ip = kpbgn(k)
          ic = kcbgn(k)
          call adjmd2cr(nx,ny,wk(ip),wk(ic))
        end do
!
!     execute one full multigrid cycle
!
        do k=1,ngrid-1
          kcur = k
          call kcymd2cr(wk)
          nx = nxk(k+1)
          ny = nyk(k+1)
          ip = kpbgn(k+1)
          ipc = kpbgn(k)
          ncx = nxk(k)
          ncy = nyk(k)
!
!     lift or prolong approximation from k to k+1
!
          call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,&
                       nyc,nyd,intpol)
        end do
      else
!
!     adjust rhs at finest grid level only
!
        nx = nxk(ngrid)
        ny = nyk(ngrid)
        ip = kpbgn(ngrid)
        ic = kcbgn(ngrid)
        call adjmd2cr(nx,ny,wk(ip),wk(ic))
      end if
!
!     execute maxcy more multigrid k cycles from finest level
!
      kcur = ngrid
      do iter=1,maxcy
        itero = iter
        call kcymd2cr(wk)
        if (tolmax.gt.0.0_r8) then
!
!      error control
!
          relmax = 0.0_r8
          phmax = 0.0_r8
          do j=1,nfy
            jj = j*(nfx+2)
            do i=1,nfx
              ij = jj+i+1
!             phmax = amax1(phmax,abs(wk(ij)))
!             relmax = amax1(relmax,abs(wk(ij)-phif(i,j)))
              phmax = max(phmax,abs(wk(ij)))
              relmax = max(relmax,abs(wk(ij)-phif(i,j)))
              phif(i,j) = wk(ij)
            end do
          end do
!
!     set maximum relative difference and check for convergence
!
          if (phmax.gt.0.0_r8) relmax = relmax/phmax
          if (relmax.le.tolmax) return
        end if
      end do
!
!     set final interate after maxcy cycles in phif
!
      do j=1,nfy
        jj = j*(nfx+2)
        do i=1,nfx
          ij = jj+i+1
          phif(i,j) = wk(ij)
        end do
      end do
      return
      end subroutine mud2cr1
!-----------------------------------------------------------------------
      subroutine kcymd2cr(wk)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     execute multigrid k cycle from kcur grid level
!     kcycle=1 for v cycles, kcycle=2 for w cycles
!
      implicit none
      real(r8) :: wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),&
        nxk(50),nyk(50),isx,jsy
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
!
!     prerelax at current finest grid level
!
      do l=1,iprer
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
!
!     restrict residual to kcur-1 level
!
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+9*ncx*ncy
      call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
!
!    set counter for grid levels to zero
!
      do l = 1,kcur
        kount(l) = 0
      end do
!
!    set new grid level and continue k-cycling
!
      klevel = kcur-1
      nrel = iprer
!
!   kcycle control point
!
   10 continue
!
!      post relax when kcur revisited
!
      if (klevel .eq. kcur) go to 5
!
!   count hit at current level
!
      kount(klevel) = kount(klevel)+1
!
!   relax at current level
!
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
!
!     kcycle complete at klevel
!
        ipc = ip
        ip = kpbgn(klevel+1)
        ncx = nxk(klevel)
        ncy = nyk(klevel)
        nx = nxk(klevel+1)
        ny = nyk(klevel+1)
!
!    inject correction to finer grid
!
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,&
                  intpol,wk(kps))
!
!    reset counter to zero
!
        kount(klevel) = 0
!
!     ascend to next higher level and set to postrelax there
!
        klevel = klevel+1
        nrel = ipost
        go to 10
      else
        if (klevel .gt. 1) then
!
!    kcycle not complete so descend unless at coarsest grid
!
          ipc = kpbgn(klevel-1)
          ncx = nxk(klevel-1)
          ncy = nyk(klevel-1)
          irc = kcbgn(klevel-1)+9*ncx*ncy
          call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),&
                      wk(kps))
!
!     prerelax at next coarser level
!
          klevel = klevel-1
          nrel = iprer
          go to 10
        else
!
!    postrelax at coarsest level
!
          do l=1,ipost
            call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
          end do
          ipc = ip
          ip = kpbgn(2)
          ncx = nxk(1)
          ncy = nyk(1)
          nx = nxk(2)
          ny = nyk(2)
!
!    inject correction to level 2
!
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,&
                  intpol,wk(kps))
!
!     set to postrelax at level 2
!
          nrel = ipost
          klevel = 2
          go to 10
        end if
      end if
    5 continue
!
!     post relax at current finest grid level
!
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ic = kcbgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
        call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end subroutine kcymd2cr
!-----------------------------------------------------------------------
      subroutine dismd2cr(nx,ny,cf,tx,ty,wk,ier,isolve)
      use edyn_solve,only:    nc,ncee,cee,ceee
      use shr_kind_mod ,only: r8 => shr_kind_r8
      use cam_abortutils   ,only: endrun
!
!     discretize elliptic pde for mud2cr, set nonfatal errors
!
      implicit none
      integer,intent(in) :: isolve
      integer nx,ny,i,j,l,im1,jm1,ier,nnx,nny
      real(r8) :: cf(nx,ny,10),tx(nx,ny,*),ty(ny,nx,*)
      real(r8) :: wk(*)
      integer        intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps

      real(r8) ::           xa,xb,yc,yd,tolmax,relmax
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
!
!     CHECK FOR CONSISTENCYT WRT KLEVEL
!
      NNX = ixp*2**(KLEVEL-1)+1
      NNY = jyq*2**(KLEVEL-1)+1
      IF(NNX.NE.NX.OR.NNY.NE.NY)THEN
        call endrun('dismd2cr in mud')
      ENDIF
      if (isolve >= 0) then
        call ceee(cee(nc(6-klevel)),nx,ny,cf)
      endif
!
!     set coefficient for specified boundaries
!
      if (nxa.eq.1) then
        i = 1
        do j=1,ny
          do l=1,9
            cf(i,j,l) = 0.0_r8
          end do
          cf(i,j,9) = 1.0_r8
        end do
      end if
      if (nxb.eq.1) then
        i = nx
        do j=1,ny
          do l=1,9
            cf(i,j,l) = 0.0_r8
          end do
          cf(i,j,9) = 1.0_r8
        end do
      end if
      if (nyc.eq.1) then
        j = 1
        do i=1,nx
          do l=1,9
            cf(i,j,l) = 0.0_r8
          end do
          cf(i,j,9) = 1.0_r8
        end do
      end if
      if (nyd.eq.1) then
        j = ny
        do i=1,nx
          do l=1,9
            cf(i,j,l) = 0.0_r8
          end do
          cf(i,j,9) = 1.0_r8
        end do
      end if
!
!     set and factor tridiagonal matrices for line relaxation(s) if flagged
!
      if (method.eq.1.or.method.eq.3) then
        if (nxa.ne.0) then
!
!    nonperiodic x line relaxation
!
          do i=1,nx
            im1 = max0(i-1,1)
            do j=1,ny
              tx(im1,j,1) = cf(i,j,5)
              tx(i,j,2) = cf(i,j,9)
              tx(i,j,3) = cf(i,j,1)
            end do
          end do
          call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
        else
!
!     periodic x line relaxation
!
          if (nx .gt. 3) then
!
!     set and factor iff nx > 3
!
            do i=1,nx-1
              do j=1,ny
                tx(i,j,1) = cf(i,j,5)
                tx(i,j,2) = cf(i,j,9)
                tx(i,j,3) = cf(i,j,1)
              end do
            end do
            call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),&
                        tx(1,1,5),wk(kps))
          end if
        end if
      end if

      if (method.eq.2.or.method.eq.3) then
        if (nyc.ne.0) then
!
!     nonperiodic y line relaxation
!
          do j=1,ny
            jm1 = max0(j-1,1)
            do i=1,nx
              ty(jm1,i,1) = cf(i,j,7)
              ty(j,i,2) = cf(i,j,9)
              ty(j,i,3) = cf(i,j,3)
            end do
          end do
          call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
        else
!
!      periodic y line relaxation
!
          if (ny .gt. 3) then
!
!     set and factor iff ny > 3
!
            do j=1,ny-1
              do i=1,nx
                ty(j,i,1) = cf(i,j,7)
                ty(j,i,2) = cf(i,j,9)
                ty(j,i,3) = cf(i,j,3)
              end do
            end do
            call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),&
                        ty(1,1,5),wk(kps))
          end if
        end if
      end if
      return
      end subroutine dismd2cr
!-----------------------------------------------------------------------
      subroutine adjmd2cr(nx,ny,phi,cf)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     adjust righthand side in cf(i,j,10) for boundary conditions
!
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j
      real(r8) :: cf(nx,ny,10),phi(0:nx+1,0:ny+1)
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
!
!     set specified boundaries in rhs from phi
!
      if (nxa.eq.1) then
        i = 1
        do j=1,ny
          cf(i,j,10) = phi(i,j)
        end do
      end if
      if (nxb.eq.1) then
        i = nx
        do j=1,ny
          cf(i,j,10) = phi(i,j)
        end do
      end if
      if (nyc.eq.1) then
        j = 1
        do i=1,nx
          cf(i,j,10) = phi(i,j)
        end do
      end if
      if (nyd.eq.1) then
        j = ny
        do i=1,nx
          cf(i,j,10) = phi(i,j)
        end do
      end if
      return
      end subroutine adjmd2cr
!-----------------------------------------------------------------------
      subroutine resmd2cr(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     restrict residual from fine to coarse mesh using fully weighted
!     residual restriction
!
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real(r8) :: rhsc(ncx,ncy),resf(nx,ny)
      real(r8) :: phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real(r8) :: cof(nx,ny,10)
!
!     set phic zero
!
      do jc=0,ncy+1
        do ic=0,ncx+1
          phic(ic,jc) = 0.0_r8
        end do
      end do
!
!     compute residual on fine mesh in resf
!
      do j=1,ny
        do i=1,nx
          resf(i,j) = cof(i,j,10)-(              &
                      cof(i,j,1)*phi(i+1,j)+     &
                      cof(i,j,2)*phi(i+1,j+1)+   &
                      cof(i,j,3)*phi(i,j+1)+     &
                      cof(i,j,4)*phi(i-1,j+1)+   &
                      cof(i,j,5)*phi(i-1,j)+     &
                      cof(i,j,6)*phi(i-1,j-1)+   &
                      cof(i,j,7)*phi(i,j-1)+     &
                      cof(i,j,8)*phi(i+1,j-1)+   &
                      cof(i,j,9)*phi(i,j))
        end do
      end do
!
!     restrict resf to coarse mesh in rhsc
!
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end subroutine resmd2cr
!-----------------------------------------------------------------------
      subroutine relmd2cr(nx,ny,phi,cof,tx,ty,sum)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     relaxation for mud2
!
      implicit none
      integer nx,ny
      real(r8) :: phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
        call relmd2crp(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
        call slxmd2cr(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
        call slymd2cr(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
        call slxmd2cr(nx,ny,phi,cof,tx,sum)
        call slymd2cr(nx,ny,phi,cof,ty,sum)
      end if
      return
      end subroutine relmd2cr
!-----------------------------------------------------------------------
      subroutine relmd2crp(nx,ny,phi,cof)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     gauss-seidel four color point relaxation
!
      implicit none
      integer nx,ny,i,j,lcolor,i1,i2,i3,i4,it
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real(r8) :: phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      i1 = 1
      i2 = 4
      i3 = 3
      i4 = 2
!
!     sweep four colored grid points
!
      do lcolor=1,4
!$OMP PARALLEL DO SHARED(i1,cof,phi,nx,ny) PRIVATE(i,j)
        do j=1,ny,4
          do i=i1,nx,4
              phi(i,j) = (cof(i,j,10) - (           &
                          cof(i,j,1)*phi(i+1,j)   + &
                          cof(i,j,2)*phi(i+1,j+1) + &
                          cof(i,j,3)*phi(i,j+1)   + &
                          cof(i,j,4)*phi(i-1,j+1) + &
                          cof(i,j,5)*phi(i-1,j)   + &
                          cof(i,j,6)*phi(i-1,j-1) + &
                          cof(i,j,7)*phi(i,j-1)   + &
                          cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
          end do
        end do
!$OMP PARALLEL DO SHARED(i2,cof,phi,nx,ny) PRIVATE(i,j)
        do j=2,ny,4
          do i=i2,nx,4
              phi(i,j) = (cof(i,j,10) - (           &
                          cof(i,j,1)*phi(i+1,j)   + &
                          cof(i,j,2)*phi(i+1,j+1) + &
                          cof(i,j,3)*phi(i,j+1)   + &
                          cof(i,j,4)*phi(i-1,j+1) + &
                          cof(i,j,5)*phi(i-1,j)   + &
                          cof(i,j,6)*phi(i-1,j-1) + &
                          cof(i,j,7)*phi(i,j-1)   + &
                          cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
          end do
        end do
!$OMP PARALLEL DO SHARED(i3,cof,phi,nx,ny) PRIVATE(i,j)
        do j=3,ny,4
          do i=i3,nx,4
              phi(i,j) = (cof(i,j,10) - (           &
                          cof(i,j,1)*phi(i+1,j)   + &
                          cof(i,j,2)*phi(i+1,j+1) + &
                          cof(i,j,3)*phi(i,j+1)   + &
                          cof(i,j,4)*phi(i-1,j+1) + &
                          cof(i,j,5)*phi(i-1,j)   + &
                          cof(i,j,6)*phi(i-1,j-1) + &
                          cof(i,j,7)*phi(i,j-1)   + &
                          cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
          end do
        end do
!$OMP PARALLEL DO SHARED(i4,cof,phi,nx,ny) PRIVATE(i,j)
        do j=4,ny,4
          do i=i4,nx,4
              phi(i,j) = (cof(i,j,10) - (           &
                          cof(i,j,1)*phi(i+1,j)   + &
                          cof(i,j,2)*phi(i+1,j+1) + &
                          cof(i,j,3)*phi(i,j+1)   + &
                          cof(i,j,4)*phi(i-1,j+1) + &
                          cof(i,j,5)*phi(i-1,j)   + &
                          cof(i,j,6)*phi(i-1,j-1) + &
                          cof(i,j,7)*phi(i,j-1)   + &
                          cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
          end do
        end do
!
!     set periodic virtual boundaries as necessary
!
        if (nxa.eq.0) then
          do j=1,ny
            phi(0,j) = phi(nx-1,j)
            phi(nx+1,j) = phi(2,j)
          end do
        end if
        if (nyc.eq.0) then
          do i=1,nx
            phi(i,0) = phi(i,ny-1)
            phi(i,ny+1) = phi(i,2)
          end do
        end if
!
!    permute (i1,i2,i3,i4) for next color
!
        it = i4
        i4 = i3
        i3 = i2
        i2 = i1
        i1 = it
      end do
      return
      end subroutine relmd2crp
!-----------------------------------------------------------------------
      subroutine slxmd2cr(nx,ny,phi,cof,tx,sum)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     line relaxation in the x direction (periodic or nonperiodic)
!
      implicit none

      integer nx,ny,i,ib,j,ii
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,&
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,&
                   kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,&
                     iguess, maxcy,method,nwork,lwork,itero,ngrid,&
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real(r8) :: phi(0:nx+1,0:ny+1),cof(nx,ny,10),tx(nx,ny,*),sum(ny)
      real(r8) :: starttime,endtime
!
!     replace line x with point gauss-seidel if
!     x direction is periodic and nx = 3 (coarsest)
!
      if (nxa .eq. 0 .and. nx .eq. 3) then
        call relmd2crp(nx,ny,phi,cof)
        return
      end if
!
!     set periodic y virtual boundary if necessary
!
      if (nyc.eq.0) then
        do i=1,nx
          phi(i,0) = phi(i,ny-1)
          phi(i,ny+1) = phi(i,2)
        end do
      end if

      if (nxa.ne.0) then
!
!     x direction not periodic, sweep odd j lines
!
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
        do j=1,ny,2
          do i=1,nx
            phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,3)*phi(i,j+1)+   &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,7)*phi(i,j-1)+   &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
!
!     forward sweep
!
          do i=2,nx
            phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
          end do
!
!     backward sweep
!
          phi(nx,j) = phi(nx,j)/tx(nx,j,2)
          do ib=2,nx
            i = nx-ib+1
            phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
          end do
        end do
!
!     sweep even j lines forward and back
!
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
        do j=2,ny,2
          do i=1,nx
            phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,3)*phi(i,j+1)+   &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,7)*phi(i,j-1)+   &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
          do i=2,nx
            phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
          end do
          phi(nx,j) = phi(nx,j)/tx(nx,j,2)
          do ib=2,nx
            i = nx-ib+1
            phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
          end do
        end do
      else
!
!     x direction periodic
!
        do j=1,ny
          sum(j) = 0.0_r8
          phi(0,j) = phi(nx-1,j)
          phi(nx+1,j) = phi(2,j)
        end do
!
!      sweep odd lines forward and back
!
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
        do j=1,ny,2
          do i=1,nx-1
            phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,3)*phi(i,j+1)+   &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,7)*phi(i,j-1)+   &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
!
!     forward sweep
!
          do i=2,nx-2
            phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
          end do
          do i=1,nx-2
            sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
          end do
          phi(nx-1,j) = phi(nx-1,j)-sum(j)
!
!     backward sweep
!
          phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
          phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/ &
                         tx(nx-2,j,2)
          do ib=4,nx
            i = nx-ib+1
            phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)* &
                       phi(nx-1,j))/tx(i,j,2)
          end do
        end do
!
!     set periodic and virtual points for j odd
!
        do j=1,ny,2
          phi(nx,j) = phi(1,j)
          phi(0,j) = phi(nx-1,j)
          phi(nx+1,j) = phi(2,j)
        end do
!
!     sweep even j lines
!
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
        do j=2,ny,2
          do i=1,nx-1
            phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,3)*phi(i,j+1)+   &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,7)*phi(i,j-1)+   &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
!
!     forward sweep
!
          do i=2,nx-2
            phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
          end do
          do i=1,nx-2
            sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
          end do
          phi(nx-1,j) = phi(nx-1,j)-sum(j)
!
!     backward sweep
!
          phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
          phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/ &
                         tx(nx-2,j,2)
          do ib=4,nx
            i = nx-ib+1
            phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)* &
                       phi(nx-1,j))/tx(i,j,2)
          end do
        end do
!
!     set periodic and virtual points for j even
!
        do j=2,ny,2
          phi(nx,j) = phi(1,j)
          phi(0,j) = phi(nx-1,j)
          phi(nx+1,j) = phi(2,j)
        end do
      end if
!
!     set periodic y virtual boundaries if necessary
!
      if (nyc.eq.0) then
        do i=1,nx
          phi(i,0) = phi(i,ny-1)
          phi(i,ny+1) = phi(i,2)
        end do
      end if
      return
      end subroutine slxmd2cr
!-----------------------------------------------------------------------
      subroutine slymd2cr(nx,ny,phi,cof,ty,sum)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      implicit none

      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real(r8) :: phi(0:nx+1,0:ny+1),cof(nx,ny,10),ty(ny,nx,*),sum(nx)
      real(r8) :: starttime,endtime
!
!     replace line y with point gauss-seidel if
!     y direction is periodic and ny = 3
!
      if (nyc .eq. 0 .and. ny .eq. 3) then
        call relmd2crp(nx,ny,phi,cof)
        return
      end if
!
!      set periodic and virtual x boundaries if necessary
!
      if (nxa.eq.0) then
        do j=1,ny
          phi(0,j) = phi(nx-1,j)
          phi(nx,j) = phi(1,j)
          phi(nx+1,j) = phi(2,j)
        end do
      end if

      if (nyc.ne.0) then
!
!     y direction not periodic
!
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
        do i=1,nx,2
          do j=1,ny
            phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+   &
                                    cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,5)*phi(i-1,j)+   &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
!
!     forward sweep thru odd x lines
!
          do j=2,ny
            phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
          end do
!
!      backward sweep
!
          phi(i,ny) = phi(i,ny)/ty(ny,i,2)
          do jb=2,ny
            j = ny-jb+1
            phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
          end do
        end do
!
!     forward sweep even x lines
!
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
        do i=2,nx,2
          do j=1,ny
            phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+   &
                                    cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,5)*phi(i-1,j)+   &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
          do j=2,ny
            phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
          end do
!
!      backward sweep
!
          phi(i,ny) = phi(i,ny)/ty(ny,i,2)
          do jb=2,ny
            j = ny-jb+1
            phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
          end do
        end do
      else
!
!     y direction periodic
!
        do i=1,nx
          sum(i) = 0.0_r8
          phi(i,0) = phi(i,ny-1)
          phi(i,ny) = phi(i,1)
          phi(i,ny+1) = phi(i,2)
        end do
!
!     forward sweep odd x lines
!
!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
        do i=1,nx,2
          do j=1,ny-1
            phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+   &
                                    cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,5)*phi(i-1,j)+   &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
          do j=2,ny-2
            phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
          end do
          do j=1,ny-2
            sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
          end do
          phi(i,ny-1) = phi(i,ny-1)-sum(i)
!
!     backward sweep
!
          phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
          phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/ &
                         ty(ny-2,i,2)
          do jb=4,ny
            j = ny-jb+1
            phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)* &
                        phi(i,ny-1))/ty(j,i,2)
          end do
        end do
!
!       set odd periodic and virtual y boundaries
!
        do i=1,nx,2
          phi(i,0) = phi(i,ny-1)
          phi(i,ny) = phi(i,1)
          phi(i,ny+1) = phi(i,2)
        end do
!
!     forward sweep even x lines
!
!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
        do i=2,nx,2
          do j=1,ny-1
            phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+   &
                                    cof(i,j,2)*phi(i+1,j+1)+ &
                                    cof(i,j,4)*phi(i-1,j+1)+ &
                                    cof(i,j,5)*phi(i-1,j)+   &
                                    cof(i,j,6)*phi(i-1,j-1)+ &
                                    cof(i,j,8)*phi(i+1,j-1))
          end do
          do j=2,ny-2
            phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
          end do
          do j=1,ny-2
            sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
          end do
          phi(i,ny-1) = phi(i,ny-1)-sum(i)
!
!     backward sweep
!
          phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
          phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/ &
                         ty(ny-2,i,2)
          do jb=4,ny
            j = ny-jb+1
            phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)* &
                        phi(i,ny-1))/ty(j,i,2)
          end do
        end do
!
!       set even periodic and virtual y boundaries
!
        do i=2,nx,2
          phi(i,0) = phi(i,ny-1)
          phi(i,ny) = phi(i,1)
          phi(i,ny+1) = phi(i,2)
        end do
      end if
!
!      set periodic and virtual x boundaries if necessary
!
      if (nxa.eq.0) then
        do j=1,ny
          phi(0,j) = phi(nx-1,j)
          phi(nx+1,j) = phi(2,j)
        end do
      end if

      return
      end subroutine slymd2cr
!-----------------------------------------------------------------------
