module edyn_maggrid
   use shr_kind_mod,   only : r8 => shr_kind_r8            ! 8-byte reals
   use cam_logfile,    only: iulog
   use edyn_params,    only: finit

   implicit none

   !
   ! Global geomagnetic grid:
   !
   integer, protected ::       &
        nmlat, &   ! number of mag latitudes
        nmlath, &  ! index of magnetic equator
        nmlon, &   ! number of mag longitudes
        nmlonp1    ! number of longitudes plus periodic point

   !
   ! geomagnetic grid resolution parameters:
   !
   integer, protected :: res_nlev
   integer, protected :: res_ngrid

   !
   ! Mag grid coordinates:
   !
   real(r8), allocatable, protected :: &
        ylatm(:),   & ! magnetic latitudes (radians)
        ylonm(:),   & ! magnetic longitudes (radians)
        gmlat(:),   & ! magnetic latitudes (degrees)
        gmlon(:)      ! magnetic longitudes (degrees)
   real(r8), protected :: dlonm,dlatm
   !
   ! Level coordinates will be same as geographic levels:
   !
   integer, protected :: nmlev ! number of levels (same as nlev in geographic)

   real(r8), allocatable, protected :: &
        rcos0s(:),    & ! cos(theta0)/cos(thetas)
        dt0dts(:),    & ! d(theta0)/d(thetas)
        dt1dts(:)       ! dt0dts/abs(sinim) (non-zero at equator)


   real(r8), protected :: table(91,2) = finit

   logical, private :: debug = .false. ! set true for prints to stdout at each call

 contains

   !-----------------------------------------------------------------------
   subroutine alloc_maggrid( mag_nlon, mag_nlat, mag_nlev, mag_ngrid )

     integer, intent(in) :: mag_nlon, mag_nlat, mag_nlev, mag_ngrid

     res_nlev = mag_nlev
     res_ngrid = mag_ngrid

     nmlat   = mag_nlat    ! number of mag latitudes
     nmlath  = (nmlat+1)/2 ! index of magnetic equator
     nmlon   = mag_nlon    ! number of mag longitudes
     nmlonp1 = nmlon+1     ! number of longitudes plus periodic point

     allocate(ylatm(nmlat))
     allocate(ylonm(nmlonp1))
     allocate(gmlat(nmlat))
     allocate(gmlon(nmlonp1))
     allocate(rcos0s(nmlat))
     allocate(dt0dts(nmlat))
     allocate(dt1dts(nmlat))

   end subroutine alloc_maggrid

   !-----------------------------------------------------------------------
   subroutine set_maggrid()
      use edyn_params, only: pi, pi_dyn, rtd, r0
      use edyn_mpi,    only: nlev => nlev_geo
      !
      ! Local:
      integer :: i, j, n
      real(r8) :: tanths2, dtheta, real8
      real(r8) :: tanth0(nmlat)
      real(r8) :: tanths(nmlat)
      real(r8) :: theta0(nmlat)
      real(r8) :: hamh0(nmlat)

      real(r8), parameter :: e    = 1.e-6_r8
      real(r8), parameter :: r1   = 1.06e7_r8
      real(r8), parameter :: alfa = 1.668_r8

      real(r8) :: table2(91, 3:5)

      real8 = real(nmlat-1, r8)
      dlatm = pi_dyn / real8
      real8 = real(nmlon, r8)
      dlonm = 2._r8 * pi_dyn / real8
      !
      ! ylatm is equally spaced in theta0, but holds the corresponding
      !   value of thetas.
      !
      do j = 1, nmlat
         real8 = real(j-1, r8)
         theta0(j) = -pi_dyn/2._r8+real8*dlatm ! note use of pi_dyn
      end do ! j=1,nmlat
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
         dt1dts(j) =                                                          &
              (r0*sqrt(1._r8+4._r8*tanths2)*(1._r8+tanths2))/                 &
              (r1*(1._r8+tanth0(j)**2)+2._r8*r0*tanth0(j)**(2._r8*alfa+1._r8)* &
              (1._r8+alfa+tanth0(j)**2)/(1._r8+tanth0(j)**2)**alfa)
         dt0dts(j) = dt1dts(j)*2._r8*tanths(j)/sqrt(1._r8+4._r8*tanths2)
      end do ! j=2,nmlat-1
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
         real8 = real(i-1, r8)
         ylonm(i) = -pi+real8*dlonm
         !     ylonm(i) = real8*dlonm
      end do ! i=1,nmlonp1
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
      dtheta = pi / 180._r8
      do i = 2, 91
         table(i,1) = table(i-1,1)+dtheta
      end do
      do i=2,90
         table2(i,4) = tan(table(i,1))
         table(i,2) = table(i,1)
      end do ! i=2,90
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
         end do ! i=2,90
      end do ! n=1,7

      if (debug) then
         write(iulog,"('set_maggrid: table= ',/,(6e12.4))") table
         write(iulog,"('set_maggrid: table2=',/,(6e12.4))") table2
      end if

   end subroutine set_maggrid
   !-----------------------------------------------------------------------
end module edyn_maggrid
