!-----------------------------------------------------------------------
! extracted from components/cam/src/ionosphere/waccmx/heelis.F90
! to be used by waccm physics
!-----------------------------------------------------------------------
module heelis_mod

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use aurora_params, only: offc, dskofc, phin, phid, theta0
  use aurora_params, only: ctpoten, hpower, plevel, dskofa, rrad
  
  implicit none
  private

  public :: heelis_update
  public :: heelis_flwv32

! private data
! Auroral parameters (taken from aurora.F of timegcm):
! (dimension of 2 is for south,north hemispheres)
!
  real(r8) ::    &
    psim(2)     ,& ! night convection entrance in MLT converted to radians (f(By))
    psie(2)     ,& !
    pcen(2)     ,& !
    phidp0(2)   ,& !
    phidm0(2)   ,& !
    phinp0(2)   ,& !
    phinm0(2)   ,& !
    rr1(2)         !

contains
!-----------------------------------------------------------------------
  subroutine heelis_update()
    use physconst, only: pi
    use mag_parms, only: get_mag_parms
!
! This is called at every timestep because ctpoten may change with time.
! Time-dependent ctpoten (kV) is read from TIMEGCM and WACCM input files, 
! unless it was provided as a constant by the user via namelist (namelist.F90).
!
    real(r8) :: rcp, rhp, arad
    real(r8) :: byloc      ! local By; now is just a hook, and set to 0. 
    integer,  parameter :: isouth = 1
    integer,  parameter :: inorth = 2
    real(r8), parameter :: h2deg = 15._r8   ! hour to degree
    real(r8), parameter :: dtr = pi/180._r8

    call get_mag_parms( hpower = hpower, ctpoten = ctpoten )

    byloc  = 0._r8

    offc(:)   =  1._r8*dtr
    dskofc(:) =  0._r8
    phin(:)   =  180._r8*dtr

    phid(isouth) = (9.39_r8 + 0.21_r8*byloc - 12._r8) * h2deg * dtr   ! In keeping with TIE-GCM2.0, phid also changed in mo_aurora.F90
    phid(inorth) = (9.39_r8 - 0.21_r8*byloc - 12._r8) * h2deg * dtr
    phin(isouth) = (23.50_r8 + 0.15_r8*byloc - 12._r8) * h2deg * dtr
    phin(inorth) = (23.50_r8 - 0.15_r8*byloc - 12._r8) * h2deg * dtr
    psim(:) =  0.44_r8 * ctpoten * 1000._r8
    psie(:) = -0.56_r8 * ctpoten * 1000._r8
    pcen(isouth) = (-0.168_r8 + 0.027_r8*byloc) * ctpoten * 1000._r8
    pcen(inorth) = (-0.168_r8 - 0.027_r8*byloc) * ctpoten * 1000._r8

    phidp0(:) =  90._r8*dtr
    phidm0(:) =  90._r8*dtr
    phinp0(:) =  90._r8*dtr
    phinm0(:) =  90._r8*dtr
    rr1(:)    = -2.6_r8
    theta0(:) = (-3.80_r8+8.48_r8*(ctpoten**0.1875_r8))*dtr

    dskofa(:) = 0._r8

    if( hpower >= 1.0_r8 ) then
       plevel = 2.09_r8*log( hpower )
    else
       plevel = 0._r8
    end if
    rhp          = 14.20_r8 + 0.96_r8*plevel
    rcp          = -0.43_r8 + 9.69_r8 * (ctpoten**.1875_r8)
    arad         = max( rcp,rhp )
    rrad(:)      = arad*dtr
    
  end subroutine heelis_update

!-----------------------------------------------------------------------
  subroutine heelis_flwv32(dlat,dlon,ratio,pi_in,iflag,nmlon,poten)
!     
! Calculate heelis potential at current magnetic latitude mlat.
!     
! 
    ! Args:
    real(r8),                    intent(in)    :: pi_in
    integer,                     intent(in)    :: nmlon
    integer,                     intent(inout) :: iflag(nmlon)
    real(r8),dimension(nmlon),   intent(in)    :: dlat,dlon,ratio
    real(r8),dimension(nmlon+1), intent(out)   :: poten
!
! Local:
    integer :: i,n,ihem
    real(r8),parameter :: eps=1.e-6_r8
    real(r8) :: &
      pi2,pih,sinthr1,psi(8),phirc,sinth0, &
      ofdc,cosofc(2),sinofc(2),aslonc(2),  &
      phdpmx(2),phnpmx(2),phnmmx(2),phdmmx(2)
    real(r8),dimension(nmlon) :: sinlat,coslat,sinlon,coslon,alon, &
      colat,wk1,wk2,wk3,phifun,phifn2
    integer :: ifn(nmlon)
    real(r8) :: phi(nmlon,8)
!
    pi2 = 2.0_r8*pi_in
    pih = 0.5_r8*pi_in

    do n=1,2
!
      ofdc = sqrt(offc(n)**2+dskofc(n)**2)
      cosofc(n) = cos(ofdc)
      sinofc(n) = sin(ofdc)
      aslonc(n) = asin(dskofc(n)/ofdc)
!
      if (phin(n) < phid(n)) phin(n) = phin(n)+pi2  ! modifies aurora phin
      phdpmx(n) = .5_r8*min(pi_in,(phin(n)-phid(n)))
      phnpmx(n) = .5_r8*min(pi_in,(phid(n)-phin(n)+pi2))
      phnmmx(n) = phdpmx(n)
      phdmmx(n) = phnpmx(n)
    enddo ! n=1,2
!
! Set ihem=1,2 for South,North hemisphere: 
!
    ihem = int( dlat(max0(1,nmlon/2))*2._r8/3.1416_r8 + 2._r8 )
    sinth0 = sin(theta0(ihem))
!
! Average amie results show r1=-2.6 for 11.3 degrees
!   (0.1972 rad) beyond theta0.
!
    sinthr1 = sin(theta0(ihem)+0.1972_r8)
    psi(1) = psie(ihem)
    psi(3) = psim(ihem)
    do n=2,4,2
      psi(n) = psi(n-1)
    enddo ! n=2,4,2
    do n=1,4
      psi(n+4) = psi(n)
    enddo ! n=1,4
!
! Transform to auroral circle coordinates:
!
    do i=1,nmlon
      sinlat(i) = sin(abs(dlat(i)))
      coslat(i) = cos(dlat(i))
      sinlon(i) = sin(dlon(i)+aslonc(ihem))
      coslon(i) = cos(dlon(i)+aslonc(ihem))
      colat(i) = cosofc(ihem)*sinlat(i)-sinofc(ihem)*coslat(i)* &
        coslon(i) 
      alon(i) = mod(atan2(sinlon(i)*coslat(i),sinlat(i)* &
        sinofc(ihem)+cosofc(ihem)*coslat(i)*coslon(i))- &
        aslonc(ihem)+3._r8*pi_in,pi2)-pi_in
      colat(i) = acos(colat(i))*sqrt(ratio(i))
!
! Boundaries for longitudinal function:
!
      wk1(i) = ((colat(i)-theta0(ihem))/theta0(ihem))**2
      phi(i,4)=phid(ihem)+eps-min(phidm0(ihem)+wk1(i)* &
        (pih-phidm0(ihem)),phdmmx(ihem))
      phi(i,5)=phid(ihem)-eps+min(phidp0(ihem)+wk1(I)* &
        (pih-phidp0(ihem)),phdpmx(ihem))
      phi(i,6)=phin(ihem)+eps-min(phinm0(ihem)+wk1(i)* &
        (pih-phinm0(ihem)),phnmmx(ihem))
      phi(i,7)=phin(ihem)-eps+min(phinp0(ihem)+wk1(i)* &
        (pih-phinp0(ihem)),phnpmx(ihem))
      phi(i,1)=phi(i,5)-pi2
      phi(i,2)=phi(i,6)-pi2
      phi(i,3)=phi(i,7)-pi2
      phi(i,8)=phi(i,4)+pi2
      phifun(i)=0._r8
      phifn2(i) = 0._r8
      if (colat(i)-theta0(ihem) >= 0._r8) then
        ifn(i) = 3
      else
        ifn(i) = 2
      endif
      if (iflag(i) == 1) iflag(i) = ifn(i)
!
! Add ring current rotation to potential (phirc)
!
      phirc = 0._r8
      wk2(i) = mod(alon(i)+phirc+2._r8*pi2+pi_in,pi2)-pi_in
      wk3(i) = mod(alon(i)+phirc+3._r8*pi2,pi2)-pi_in
    enddo ! i=1,nmlon
!
! Longitudinal variation:
!
    do n=1,7
      do i=1,nmlon
        phifun(i)=phifun(i)+.25_r8*(psi(n)+psi(n+1)+(psi(n)-       &
          psi(n+1))*cos(mod(pi_in*(wk2(i)-phi(i,n))/(phi(i,n+1)-  &
          phi(i,n)),pi2)))*(1._r8-sign(1._r8,(wk2(i)-phi(i,n))*    &
          (wk2(i)-phi(i,n+1))))
        phifn2(i)=phifn2(i)+.25_r8*(psi(n)+psi(n+1)+(psi(n)-       &
          psi(n+1))*cos(mod(pi_in*(wk3(i)-phi(i,n))/(phi(i,n+1)-  &
          phi(i,n)),pi2)))*(1._r8-sign(1._r8,(wk3(i)-phi(i,n))*    &
          (wk3(i)-phi(i,n+1))))
      enddo
    enddo
!
! Evaluate total potential:
!

    do i=1,nmlon
      if (iflag(i)==2) then
        poten(i) = (2._r8*(pcen(ihem)-phifun(i))+(phifun(i)-phifn2(i))* &
          0.75_r8)*(colat(i)/theta0(ihem))**3 +                         &
          (1.5_r8*(phifun(i)+phifn2(i))-3._r8*pcen(ihem))*(colat(i)/    &
          theta0(ihem))**2 + 0.75_r8*(phifun(i)-phifn2(i))*(colat(i)/   &
          theta0(ihem)) + pcen(ihem)
      else
        poten(i) = phifun(i)*(max(sin(colat(i)),sinth0)/sinth0)**rr1(ihem)* &
          exp(7._r8*(1._r8-max(sin(colat(i)),sinthr1)/sinthr1))
      endif
    enddo

  end subroutine heelis_flwv32
!-----------------------------------------------------------------------
end module heelis_mod
