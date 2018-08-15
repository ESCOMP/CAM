module getapex
!
! Calculate quantities needed to transform scalar fields between geographic 
! and geomagnetic coordinate systems. 
!
  use shr_kind_mod     ,only : r8 => shr_kind_r8
  use cam_logfile      ,only: iulog
  use cam_abortutils   ,only: endrun
  use edyn_geogrid     ,only: nlon,nlonp1,ylatg,ylong,dlong,&
                              jspole,jnpole
  use edyn_maggrid     ,only: nmlon,nmlonp1,nmlat,ylatm,ylonm,dlonm

  implicit none
  save

  private

  public :: get_apex
  public :: magfield, bx, by, bz, bmod2, bmod, xb, yb, zb, be3arr, dddarr, dvec
  public :: alatm, alonm, gdlondeg, gdlatdeg
  public :: rjac

  integer ::             &
    ig(nmlonp1,nmlat),   & ! geog lon grid containing each geomag point
    jg(nmlonp1,nmlat)      ! geog lat grid containing each geomag point

  real(r8) ::              &
    wt(4,nmlonp1,nmlat)      ! interpolation weights for geo2mag
 
  real(r8),dimension(nmlonp1,nmlat) :: & ! geo lat,lon coords on mag grid
                           gdlatdeg,   & ! geographic latitude of each magnetic grid point (deg)
                           gdlondeg      ! geographic longitude of each magnetic grid point (deg)
!
! Variables on geographic grid needed by other modules must
! be allocated dynamically to be grid-independent (sub alloc_apex):
!
  integer,allocatable :: & ! (nlonp1,jspole:jnpole))
    im(:,:),             & ! geomag lon grid containing each geog point
    jm(:,:)                ! geomag lat grid containing each geog point

  real(r8),allocatable ::  & ! (nlonp1,jspole:jnpole)
    dim(:,:),              & ! fraction in lon for grid interp
    djm(:,:)                 ! fraction in lat for grid interp

  real(r8),allocatable ::  & ! (nlonp1,jspole:jnpole,3,2)
    dvec(:,:,:,:)            ! vectors from apxmall

  real(r8),allocatable ::  & ! (nlonp1,jspole:jnpole)
    dddarr(:,:),           & ! from apxmall
    be3arr(:,:)              ! from apxmall

  real(r8),allocatable ::  & ! (nlonp1,jspole:jnpole)
    alatm(:,:),            & ! geomagnetic latitude at each geographic grid point (radians)
    alonm(:,:),            & ! geomagnetic longitude at each geographic grid point (radians)
    xb(:,:),               & ! northward component of magnetic field
    yb(:,:),               & ! eastward component of magnetic field
    zb(:,:),               & ! downward component of magnetic field (gauss)
    bmod(:,:)                ! magnitude of magnetic field (gauss)
!
! rjac: scaled derivatives of geomagnetic coords wrt geographic coordinates.
! rjac(1,1) = cos(thetas)/cos(theta)*d(lamdas)/d(lamda)
! rjac(1,2) = cos(thetas)*d(lamdas)/d(theta)
! rjac(2,1) = 1./cos(theta)*d(thetas)/d(lamda)
! rjac(2,2) = d(thetas)/d(theta)
! where (lamda,theta) are geographic coordinates
!       (lamdas,thetas) are geomagnetic coordinates
!
    real(r8),allocatable :: &
      rjac(:,:,:,:) ! (nlon+1,jspole:jnpole,2,2)
!
! Parameters defined by sub magfield (allocated in alloc_magfield):
!
    real(r8),allocatable,dimension(:,:) :: & ! (0:nlon+1,jspole-1:jnpole+1)
      bx,by,bz,bmod2

  contains
!-----------------------------------------------------------------------
  subroutine get_apex( )
!
! This is called once per run from main.
!
  use edyn_params,only: re_dyn,h0,hs,dtr,rtd
  use apex,       only: apex_mall,apex_q2g
  use edyn_geogrid,only: glat_edyn_geo => glat, glon_edyn_geo => glon

!
! Local:
    integer :: i,j,ier,jjm,jjg
    integer,parameter :: nalt=2
    real(r8) :: real8

    real(r8) :: rekm,h0km,alt,hr,ror03,glat,glon,&
      xlonmi,qdlon,qdlat,gdlon,gdlat,xlongi,frki,frkj

!
! Non-scalar arguments returned by apxmall:
    real(r8) ::           &
      b(3),bhat(3),       &
      d1(3),d2(3),d3(3),  &
      e1(3),e2(3),e3(3),  &
      f1(2),f2(2)
    real(r8) :: bmag,alon,xlatm,vmp,w,d,be3,si,sim,xlatqd,f

!
! Allocate arrays that are needed by other modules:
    call alloc_apex
    call alloc_magfield

    rekm = re_dyn*1.e-5_r8  ! earth radius (km)
    h0km = h0*1.e-5_r8
    alt  = hs*1.e-5_r8  ! modified apex reference altitude (km)
    hr   = alt
    ror03= ((rekm + alt)/(rekm + h0km))**3
!
! Loop over 2d geographic grid:
!
    do j=jspole,jnpole
      glat = glat_edyn_geo(j)
      do i=1,nlonp1
         if (i.eq.nlonp1) then
            glon = glon_edyn_geo(1)
         else
            glon = glon_edyn_geo(i)
         endif

        call apex_mall (                           & 
          glat,glon,alt,hr,                        & !Inputs
          b,bhat,bmag,si,                          & !Mag Fld
          alon,                                    & !Apx Lon
          xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
          xlatqd,f,f1,f2 , ier)                      !Qsi-Dpl

        if (ier /= 0) call endrun('get_apex: apxmall error')

        alatm(i,j)    =  xlatm*dtr
        alonm(i,j)    =  alon *dtr
        xb   (i,j)    =  b(2)*1.e-5_r8 ! nT -> gauss
        yb   (i,j)    =  b(1)*1.e-5_r8 ! nT -> gauss
        zb   (i,j)    = -b(3)*1.e-5_r8 ! nT -> gauss
        bmod (i,j)    =  bmag*1.e-5_r8 ! nT -> gauss

        rjac (i,j,1,1) =  f2(2)
        rjac (i,j,1,2) = -f2(1)
        rjac (i,j,2,1) = -f1(2)
        rjac (i,j,2,2) =  f1(1)
!
! Set up parameters for magnetic to geographic interpolation.
!
        xlonmi = (alonm(i,j) - ylonm(1))/dlonm
        real8 = dble(nmlon)
        if (xlonmi < 0._r8) xlonmi = xlonmi + real8
        im(i,j) = xlonmi
        real8 = dble(im(i,j))
        dim(i,j) = xlonmi - real8
        im(i,j) = im(i,j) + 1
        if (im(i,j) >= nmlonp1) im(i,j) = im(i,j) - nmlon
        alatm(i,j) = min(alatm(i,j),ylatm(nmlat))
        do jjm=2,nmlat
          if (alatm(i,j) > ylatm(jjm)) cycle
          jm(i,j) = jjm - 1
          djm(i,j) = (alatm(i,j) - ylatm(jm(i,j)))/  &
                     (ylatm(jjm) - ylatm(jm(i,j)))
          exit
        enddo
        if (j /= jspole .and. j /= jnpole) then
          dvec(i,j,1,1) = d1(1)
          dvec(i,j,2,1) = d1(2)
          dvec(i,j,3,1) = d1(3)
          dvec(i,j,1,2) = d2(1)
          dvec(i,j,2,2) = d2(2)
          dvec(i,j,3,2) = d2(3)
          dddarr(i,j)   = d
!
! Scale be3 from 130 km to a reference height of 90 km.
          be3arr(i,j)   = be3*ror03
        endif
      enddo ! i=1,nlonp1
    enddo ! j=jspole,jnpole
!
! Set up parameters for geographic to magnetic interpolation
    do i=1,nmlonp1
      qdlon = ylonm(i)*rtd
      do j=1,nmlat
        qdlat = ylatm(j)*rtd
!
! Convert from Quasi-Dipole to geographic coordinates.
! gdlat,gdlon are returned by apxq2g.
!
        call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ier)
        if (ier /= 0) then
          write(iulog,"(i3,i3,i3)") '>>> Error from apex_q2g: ier=',ier, &
            ' i=',i,' j=',j
          call endrun('get_apex: apex_q2g ier')
        endif
        gdlat = gdlat*dtr
        gdlon = gdlon*dtr
        xlongi = (gdlon - ylong(1))/dlong
        real8 = dble(nlon)
        if (xlongi < 0._r8) xlongi = xlongi + real8
        ig(i,j) = xlongi
        real8 = dble(ig(i,j))
        frki = xlongi - real8
        ig(i,j) = ig(i,j) + 1
        if (ig(i,j) >= nlonp1) ig(i,j) = ig(i,j) - nlon
        gdlat = min(gdlat,ylatg(jnpole))
        do jjg=1,jnpole
          if (gdlat > ylatg(jjg)) cycle
          jg(i,j) = jjg - 1
          frkj = (gdlat - ylatg(jg(i,j)))/(ylatg(jjg) - ylatg(jg(i,j)))
!
! 99/2/25b Add one to JG to account for the fact that AG in geo2mag has
!  a second (J) index starting at 1, while the second index of the
!  array in the calling arguments begins at 0.
!
          jg(i,j) = jg(i,j) + 1
          exit
        enddo
        wt(1,i,j) = (1._r8 - frki)*(1._r8 - frkj)
        wt(2,i,j) =          frki *(1._r8 - frkj)
        wt(3,i,j) =          frki *frkj
        wt(4,i,j) = (1._r8 - frki)*frkj
!
! gdlatdeg,gdlondeg will be coordY,coordX of the mag grid for ESMF 
! regridding (see edyn_esmf.F)
!
        gdlatdeg(i,j) = gdlat*rtd
        gdlondeg(i,j) = gdlon*rtd
      enddo ! j=1,nmlat
    enddo ! i=1,nmlonp1
  end subroutine get_apex
!-----------------------------------------------------------------------
  subroutine magfield
!
! Calculate magnetic field parameters (bx,by,bz)
! (see also TIEGCM magfield.F)
! This is called once per run and when crossing year boundary from edyn_init, after get_apex.
! All arrays are on the global domain, all processors execute.
!
! Local:
    integer :: i,j
!
! QUESTION: in TIEGCM, dipmin is resolution dependent - how do we
!   handle this for different resolutions in WACCM?
!
!   real(r8),parameter :: dipmin=0.17 ! set for 5.0-deg TIEGCM (also known as sin10)
    real(r8),parameter :: dipmin=0.24_r8 ! set for 2.5-deg TIEGCM (also known as sin10)
    real(r8) :: cos10

    cos10 = sqrt(1._r8-dipmin**2)
    do j=jspole,jnpole ! 1,nlat
      do i=1,nlon
        bx(i,j) = yb(i,j)/bmod(i,j)
        by(i,j) = xb(i,j)/bmod(i,j)
        bz(i,j) = -zb(i,j)/bmod(i,j)
        bmod2(i,j) = bmod(i,j)
!
! Set minimum dip to 10 degrees
        if (abs(bz(i,j))-dipmin < 0._r8) then
          bx(i,j) = bx(i,j)*(cos10/sqrt(1._r8-bz(i,j)**2))
          by(i,j) = by(i,j)*(cos10/sqrt(1._r8-bz(i,j)**2))
          bz(i,j) = sign(dipmin,bz(i,j))
        endif
      enddo ! i=1,nlon
    enddo ! j=jspole,jnpole

!
! Values at jspole-1:
    j=jspole-1 ! j=0
    do i=1,nlon
      bx(i,j)   =    -bx(1+mod(i-1+nlon/2,nlon),jspole)
      by(i,j)   =    -by(1+mod(i-1+nlon/2,nlon),jspole)
      bz(i,j)   =     bz(1+mod(i-1+nlon/2,nlon),jspole)
      bmod2(i,j) = bmod2(1+mod(i-1+nlon/2,nlon),jspole)
    enddo
!
! Values at jnpole+1:
    j=jnpole+1 ! j=nlat+1
    do i=1,nlon
      bx(i,j) =      -bx(1+mod(i-1+nlon/2,nlon),jnpole)
      by(i,j) =      -by(1+mod(i-1+nlon/2,nlon),jnpole)
      bz(i,j) =       bz(1+mod(i-1+nlon/2,nlon),jnpole)
      bmod2(i,j) = bmod2(1+mod(i-1+nlon/2,nlon),jnpole)
    enddo
!
! Periodic points:
! FIX: not sure about this, but
!   I am following tiegcm, but with a single point on each end instead of 2
!
    do j=jspole-1,jnpole+1
      bx(nlonp1,j) = bx(1,j)
      by(nlonp1,j) = by(1,j)
      bz(nlonp1,j) = bz(1,j)
      bmod2(nlonp1,j) = bmod2(1,j)

      bx(0,j) = bx(nlon,j)
      by(0,j) = by(nlon,j)
      bz(0,j) = bz(nlon,j)
      bmod2(0,j) = bmod2(nlon,j)
    enddo

  end subroutine magfield
!-----------------------------------------------------------------------
  subroutine alloc_magfield

!------------------------------------------------------------------------------------------
! Do allocations, checking if previously allocated in case of year boundary crossing
!------------------------------------------------------------------------------------------
    if (.not.allocated(bx)) allocate(bx(0:nlonp1,jspole-1:jnpole+1))
    if (.not.allocated(by)) allocate(by(0:nlonp1,jspole-1:jnpole+1))
    if (.not.allocated(bz)) allocate(bz(0:nlonp1,jspole-1:jnpole+1))
    if (.not.allocated(bmod2)) allocate(bmod2(0:nlonp1,jspole-1:jnpole+1))

  end subroutine alloc_magfield
!-----------------------------------------------------------------------
 
  subroutine alloc_apex

!------------------------------------------------------------------------------------------
! Do allocations, checking if previously allocated in case of year boundary crossing
!------------------------------------------------------------------------------------------
    if (.not.allocated(im)) allocate(im (nlonp1,jspole:jnpole))
    if (.not.allocated(jm)) allocate(jm (nlonp1,jspole:jnpole))
    if (.not.allocated(dim)) allocate(dim(nlonp1,jspole:jnpole))
    if (.not.allocated(djm)) allocate(djm(nlonp1,jspole:jnpole))

    if (.not.allocated(xb)) allocate(xb   (nlonp1,jspole:jnpole))
    if (.not.allocated(yb)) allocate(yb   (nlonp1,jspole:jnpole))
    if (.not.allocated(zb)) allocate(zb   (nlonp1,jspole:jnpole))
    if (.not.allocated(bmod)) allocate(bmod (nlonp1,jspole:jnpole))
    if (.not.allocated(alatm)) allocate(alatm(nlonp1,jspole:jnpole))
    if (.not.allocated(alonm))allocate(alonm(nlonp1,jspole:jnpole))

    if (.not.allocated(dvec)) allocate(dvec  (nlonp1,jspole:jnpole,3,2))
    if (.not.allocated(dddarr)) allocate(dddarr(nlonp1,jspole:jnpole))
    if (.not.allocated(be3arr)) allocate(be3arr(nlonp1,jspole:jnpole))

    if (.not.allocated(rjac)) allocate(rjac(nlon+1,jspole:jnpole,2,2))

  end subroutine alloc_apex
!-----------------------------------------------------------------------
end module getapex
