!-----------------------------------------------------------------------
      subroutine mudmod(pe,phi_out,jntl,isolve,ier)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      use cam_abortutils   ,only: endrun
      use edyn_solve   ,only: cee
      use cam_logfile  ,only: iulog

      implicit none

      integer jntl,ier  ! output: not converged ier < 0
      integer,intent(in) :: isolve
!
!     set grid size params
!
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
      parameter (iixp = 5 , jjyq = 3, iiex = 5, jjey = 5 )
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
!
!     estimate work space for point relaxation (see mud2cr.d)
!
      parameter (llwork=(7*(nnx+2)*(nny+2)+76*nnx*nny)/3 )
      real(r8) :: phi(nnx,nny),rhs(nnx,nny),work(llwork)
      real(r8) :: phi_out(0:nnx+1,0:nny+1)
      real(r8) :: time0,time1
!
!     put integer and floating point argument names in contiguous
!     storage for labelling in vectors iprm,fprm
!
      integer iprm(17),mgopt(4)
      real(r8) :: fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,&
                    iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny, &
                    iguess,maxcy,method,nwork,lwrkqd,itero
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real(r8) :: PE(NNX,*)
      integer maxcya
      DATA MAXCYA/50/
      integer mm,nn,jj,jjj,ij
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
      mgopt(2) = 3
      mgopt(3) = 2
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
      tolmax = 0.001_r8
!
!     set right hand side in rhs
!     initialize phi to zero
!
      do i=1,nx
        do j=1,ny
          rhs(i,j) = cee(i+(j-1)*nx+9*nx*ny)
          phi(i,j) = 0.0_r8
        end do
      end do
!
!     set specified boundaries in phi
!
      do i=1,nx
        phi(i,ny) = rhs(i,ny)/cee(i+(ny-1)*nx+8*nx*ny)
      end do
      
!     write(iulog,100)
! 100 format(//' mud2cr test ')
!     write (iulog,101) (iprm(i),i=1,15)
! 101 format(/,' integer input arguments ',/,
!    |  ' intl =  ',i2,/,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,
!    |  ' nyd =   ',i2,/,' ixp = ',i2,' jyq = ',i2,' iex = ',i2,
!    |  ' jey =   ',i2,/,' nx =  ',i3,' ny =  ',i3,' iguess = ',i2,
!    |  ' maxcy = ',i3,/,' method = ',i2, ' work space estimate = ',i7)
!     write (iulog,102) (mgopt(i),i=1,4)
! 102 format(/' multigrid option arguments ',
!    |  /,' kcycle = ',i2,
!    |  /,' iprer = ',i2,
!    |  /,' ipost = ',i2
!    |  /,' intpol = ',i2)
!     write(iulog,103) xa,xb,yc,yd,tolmax
! 103 format(/' floating point input parameters ',
!    |  /,' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
!    |  /,' tolerance (error control) =   ',e10.3)
!     write(iulog,"('fprm(1-5) (xa,xb,yc,yd,tolmax=',6f8.3)") fprm(1:5)
!
!     intialization call
!
!     write(iulog,104) intl
! 104 format(/' discretization call to mud2cr', ' intl = ', i2)

      call mud2cm(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)

!     write (iulog,200) ierror,iprm(16)
! 200 format(' ierror = ',i2, ' minimum work space = ',i7)
!     if (ierror.gt.0) call exit(0)
!
!     attempt solution
!
      intl = 1
!     write(iulog,106) intl,method,iguess
! 106 format(/' approximation call to mud2cr',
!    +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      
      call mud2cm(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
      ier = ierror ! ier < 0 not converged
      if(ier < 0 )  goto 108
      
!     write (iulog,107) ierror
! 107 format(' ierror = ',i2)
      if (ierror.gt.0) call endrun('mudmod call mud2cm')
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

! am 8/10 for calculating residual: convert work array (solution) into array
! sized as coefficient stencil (c0, cofum) including values at index 0, nmlon0+1
! and nmlat0+1     

      do j=0,ny+1
        jj = j*(nx+2)
        do i=0,nx+1
          ij = jj+i+1
          phi_out(i,j) = work(ij)
        end do
      end do
      
  108 continue   
      end subroutine mudmod
!-------------------------------------------------------------------
      subroutine mud2cm(iparm,fparm,work,rhs,phi,mgopt,ierror,isolve)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      use cam_logfile  ,only: iulog
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
      call mud2c1m(nx,ny,rhs,phi,work)
      iparm(17) = itero
      if (tolmax.gt.0.0_r8) then   ! check for convergence
        fparm(6) = relmax
        if (relmax.gt.tolmax) then
        
!          ierror = -1   ! flag convergenc failure
           write(iulog,*) "no convergence with mudmod"
!         
           iguess = 1 
           iparm(12)= iguess                               
           call mud2cr1(nx,ny,rhs,phi,work) !  solve with modified stencils
           
           fparm(6) = relmax
           if (relmax.gt.tolmax) then
             write(iulog,*) "no convergence with mud"
             ierror = -1   ! flag convergenc failure
           end if
           
        end if
      end if
      
      return
      end subroutine mud2cm
!------------------------------------------------------------------------      
      subroutine mud2c1m(nx,ny,rhsf,phif,wk)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      implicit none
      integer nx,ny
      real(r8) :: phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      real(r8) :: xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      integer iw,itx,ity,ierror
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50), &
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
          call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,wk(ipc),wk(irc))
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
          call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,nyc,nyd,intpol)
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
        call kcym2cm(wk)
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
      end subroutine mud2c1m

!------------------------------------------------------------------------
      subroutine kcym2cm(wk)
      use shr_kind_mod ,only: r8 => shr_kind_r8
      use edyn_solve,only: cofum
!
!     execute multigrid k cycle from kcur grid level
!     kcycle=1 for v cycles, kcycle=2 for w cycles
!
      implicit none
      real(r8) :: wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50), &
        nxk(50),nyk(50),isx,jsy
      integer kount(50)
!     real(r8) :: ::  cofum
!     common/mudmd/cofum(1) 
      
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
!     call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))

      call resm2cm(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic), &
                   wk(kps),cofum)
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
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd, &
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
          call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic), &
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
        call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd, &
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
      end subroutine kcym2cm
!----------------------------------------------------------------------      
      subroutine resm2cm(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf,cofum)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     restrict residual from fine to coarse mesh using fully weighted
!     residual restriction
!
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real(r8) :: rhsc(ncx,ncy),resf(nx,ny)
      real(r8) :: phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real(r8) :: cof(nx,ny,10),cofum(nx,ny,9)
      real(r8) :: l2norm
!
!     set phic zero
!
      do jc=0,ncy+1
        do ic=0,ncx+1
          phic(ic,jc) = 0.0_r8
        end do
      end do
      
      call bnd2cm(nx,ny,cofum)
!
!     compute residual on fine mesh in resf
!
      l2norm = 0._r8    
!$OMP PARALLEL DO SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
        do i=1,nx
          resf(i,j) = cof(i,j,10)-(              &
                      cofum(i,j,1)*phi(i+1,j)+   &
                      cofum(i,j,2)*phi(i+1,j+1)+ &
                      cofum(i,j,3)*phi(i,j+1)+   &
                      cofum(i,j,4)*phi(i-1,j+1)+ &
                      cofum(i,j,5)*phi(i-1,j)+   &
                      cofum(i,j,6)*phi(i-1,j-1)+ &
                      cofum(i,j,7)*phi(i,j-1)+   &
                      cofum(i,j,8)*phi(i+1,j-1)+ &
                      cofum(i,j,9)*phi(i,j))
          
              l2norm = l2norm + resf(i,j)*resf(i,j)
        end do
      end do
!
!     restrict resf to coarse mesh in rhsc
!
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end subroutine resm2cm

!-----------------------------------------------------------------------
      subroutine bnd2cm(nx,ny,cf)
      use shr_kind_mod ,only: r8 => shr_kind_r8
!
!     set stencil & boundary condition for finest stencil
!
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess, &
                   maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur, &
                   kcycle,iprer,ipost,intpol,kps
      real(r8) :: xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j,kbdy,l,im1,jm1,ier,jc,nnx,nny
      real(r8) :: cf(nx,ny,*)
      real(r8) :: dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,alfmax,cemax
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy, &
                     iguess, maxcy,method,nwork,lwork,itero,ngrid, &
                     klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
 
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
      return
      end subroutine bnd2cm
!-----------------------------------------------------------------------
