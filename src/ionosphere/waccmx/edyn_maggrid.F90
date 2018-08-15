module edyn_maggrid
  use shr_kind_mod,   only : r8 => shr_kind_r8            ! 8-byte reals
  use cam_logfile,    only: iulog
  implicit none
  save

!
! Global geomagnetic grid:
!
  integer, parameter ::    &
    nmlat = 97,           & ! number of mag latitudes
    nmlath = (nmlat+1)/2, & ! index of magnetic equator
    nmlon = 80,           & ! number of mag longitudes
    nmlonp1 = nmlon+1       ! number of longitudes plus periodic point
!
! Mag grid coordinates:
!
  real(r8) ::       &
    ylatm(nmlat),   & ! magnetic latitudes (radians)
    ylonm(nmlonp1), & ! magnetic longitudes (radians)
    gmlat(nmlat),   & ! magnetic latitudes (degrees)
    gmlon(nmlonp1), & ! magnetic longitudes (degrees)
    dlonm,dlatm
!
! Level coordinates will be same as geographic levels:
!
  integer :: nmlev    ! number of levels (same as nlev in geographic)

  real(r8) ::       &
    rcos0s(nmlat),  & ! cos(theta0)/cos(thetas)
    dt0dts(nmlat),  & ! d(theta0)/d(thetas)
    dt1dts(nmlat)     ! dt0dts/abs(sinim) (non-zero at equator)

  real(r8) :: table(91,2)

  logical :: debug=.false. ! set true for prints to stdout at each call

  contains
!-----------------------------------------------------------------------
  subroutine set_maggrid
    use edyn_params ,only: pi,pi_dyn,rtd,r0
    use edyn_geogrid,only: nlev
!
! Local:
    integer :: i,j,n
    real(r8) :: tanths2,dtheta,real8
    real(r8) ::      &
      tanth0(nmlat), &
      tanths(nmlat), &
      theta0(nmlat), &
      hamh0(nmlat)

    real(r8),parameter :: &
      e=1.e-6_r8,    &
      r1=1.06e7_r8,  &
      alfa=1.668_r8

    real(r8) :: table2(91,3:5)

    real8 = dble(nmlat-1)
    dlatm = pi_dyn/real8
    real8 = dble(nmlon)
    dlonm = 2._r8*pi_dyn/real8
!
! ylatm is equally spaced in theta0, but holds corresponding value of thetas.
!
    do j=1,nmlat
      real8 = dble(j-1)
      theta0(j) = -pi_dyn/2._r8+real8*dlatm ! note use of pi_dyn
    enddo ! j=1,nmlat
    do j=2,nmlat-1
      tanth0(j) = abs(tan(theta0(j)))
      hamh0(j) = r1*tanth0(j)+r0*tanth0(j)**(2._r8+2._r8*alfa)/ &
        (1._r8+tanth0(j)**2)**alfa
      tanths(j) = sqrt(hamh0(j)/r0)
      ylatm(j) = sign(atan(tanths(j)),theta0(j))
      rcos0s(j) = sqrt((1._r8+tanths(j)**2)/(1._r8+tanth0(j)**2))
!
! Timegcm has an alternate calculation for dt1dts and dt0dts if dynamo
! is not called.
!
      tanths2  = tanths(j)**2
      dt1dts(j) = &
        (r0*sqrt(1._r8+4._r8*tanths2)*(1._r8+tanths2))/                   &
        (r1*(1._r8+tanth0(j)**2)+2._r8*r0*tanth0(j)**(2._r8*alfa+1._r8)*  &
        (1._r8+alfa+tanth0(j)**2)/(1._r8+tanth0(j)**2)**alfa)
      dt0dts(j) = dt1dts(j)*2._r8*tanths(j)/sqrt(1._r8+4._r8*tanths2)
    enddo ! j=2,nmlat-1
!
! Magnetic poles:
!
    ylatm(1) = theta0(1)
    ylatm(nmlat)  = theta0(nmlat)
    rcos0s(1)     = 1._r8
    rcos0s(nmlat) = 1._r8
    dt0dts(1)     = 1._r8
    dt0dts(nmlat) = 1._r8
!
! Magnetic longitudes:
!
    do i=1,nmlonp1
      real8 = dble(i-1)
      ylonm(i) = -pi+real8*dlonm
!     ylonm(i) = real8*dlonm
    enddo ! i=1,nmlonp1
!
! Define mag grid in degrees, and mag levels:
!
    gmlat(:) = ylatm(:)*rtd
    gmlon(:) = ylonm(:)*rtd
!
! Magnetic levels are same as midpoint geographic levels:
!
    nmlev = nlev

!
! Calculate table:
!
    table(1,1) = 0._r8
    table(1,2) = 0._r8
    dtheta = pi/180._r8
    do i=2,91
      table(i,1) = table(i-1,1)+dtheta
    enddo
    do i=2,90
      table2(i,4) = tan(table(i,1))
      table(i,2) = table(i,1)
    enddo ! i=2,90
    table(91,2) = table(91,1)
    do n=1,7
      do i=2,90
        table2(i,3) = table(i,2)
        table(i,2) = tan(table2(i,3))
        table2(i,5) = sqrt(r1/r0*table(i,2)+table(i,2)**(2._r8*(1._r8+alfa))/ &
          (1._r8+table(i,2)**2)**alfa)
        table(i,2) = table2(i,3)-(table2(i,5)-table2(i,4))*2._r8*      &
          table2(i,5)/(r1/r0*(1._r8+table(i,2)**2)+2._r8*table(i,2)**  &
          (2._r8*alfa+1._r8)*(1._r8+alfa+table(i,2)**2)/               &
          (1._r8+table(i,2)**2)**alfa)
      enddo ! i=2,90
    enddo ! n=1,7
    
    if (debug) then
      write(iulog,"('set_maggrid: table= ',/,(6e12.4))") table
      write(iulog,"('set_maggrid: table2=',/,(6e12.4))") table2
    endif
  
  end subroutine set_maggrid
!-----------------------------------------------------------------------
end module edyn_maggrid
