module filter_module
  use shr_kind_mod,only: r8 => shr_kind_r8
  use cam_logfile, only: iulog
  use cam_abortutils,only: endrun
  use edyn_geogrid,only: nlon, nlat

  implicit none
  private
  public :: ntrigs, ifax
  public :: filter_init
  public :: filter1, filter2, ringfilter
  public :: kut1, kut2
!
! Coefficients and factors for fft. Sub setfft is called once per run from edyn_init
!
  integer :: ntrigs ! = 3*nlon/2+1
  real(r8), allocatable :: trigs(:)
  integer :: ifax(13)
!--------------------------------------------------------------------------
!
! For filter1:
!
! This is used by TIEGCM for basic filtering (t,u,v, et.al.),
! when nlat=72 (2.5 deg res):
!
!      integer,parameter :: kut(nlat) =
!    |   (/1  ,1  ,2  ,2  ,4  ,4  ,8  ,8  ,10 ,10 ,12 ,12,
!    |     15 ,15 ,18 ,18 ,22 ,22 ,26 ,26 ,30 ,30 ,32 ,32,
!    |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
!    |     34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34,
!    |     32 ,32 ,30 ,30 ,26 ,26 ,22 ,22 ,18 ,18 ,15 ,15,
!    |     12 ,12 ,10 ,10 ,8  ,8  ,4  ,4  ,2  ,2  ,1  ,1/)

  integer, allocatable, protected :: kut1(:)
  integer, allocatable, protected :: kut2(:)
  integer, allocatable, protected :: nn(:)

  integer, parameter :: nxlat = 96
  integer,parameter :: kut1_x(nxlat) =                   &
    (/0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,2  ,2  ,2  ,2 , &
      3  ,3  ,3  ,3  ,4  ,4  ,4  ,4  ,6  ,6  ,6  ,6 , &
      8  ,8  ,8  ,8  ,10 ,10 ,10 ,10 ,12 ,12 ,12 ,12, &
     16  ,16 ,18 ,18 ,20 ,20 ,22 ,22 ,24 ,24 ,26 ,26, &
     26  ,26 ,24 ,24 ,22 ,22 ,20 ,20 ,18 ,18 ,16 ,16, &
     12  ,12 ,12 ,12 ,10 ,10 ,10 ,10 ,8  ,8  ,8  ,8 , &
      6  ,6  ,6  ,6  ,4  ,4  ,4  ,4  ,3  ,3  ,3  ,3 , &
      2  ,2  ,2  ,2  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  /)
!--------------------------------------------------------------------------
!
! For filter2:
!
! This is used by TIEGCM for O+ filtering when nlat=72 (2.5 deg res):
!
!      kut2=(/0, 0, 1, 2, 4, 4, 6, 6, 8, 8,10,10,12,12,15,15,18,18,
!    |      20,20,20,20,18,18,15,12, 8, 8, 4, 4, 4, 4, 2, 2, 1, 1,
!    |       1, 1, 2, 2, 4, 4, 4, 4, 8, 8,12,15,18,18,20,20,20,20,
!    |      18,18,15,15,12,12,10,10, 8, 8, 6, 6, 4, 4, 2, 1, 0, 0/) ! 2.5 deg
!
!     nn=(/90,90,40,40,22,22,14,14,10,10, 8, 8, 6, 6, 4, 4, 2, 2,
!    |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
!    |      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
!    |      2, 2, 4, 4, 6, 6, 8, 8,10,10,14,14,22,22,40,40,90,90/) ! 2.5 deg
!
! At 1.9 deg resolution, nlat==96
!
  integer, parameter :: kut2_x(nxlat) = &
    (/  0,   0,   0,   0,   1,   2,   3,   4,   5,   5,   6,   7, &
        8,   8,   9,  10,  11,  11,  13,  14,  15,  16,  17,  18, &
       19,  19,  20,  20,  19,  19,  18,  17,  15,  13,  10,   8, &
        7,   5,   4,   4,   4,   3,   3,   2,   1,   1,   1,   1, &
        1,   1,   1,   1,   2,   3,   3,   4,   4,   4,   5,   7, &
        8,  10,  13,  15,  17,  18,  19,  19,  20,  20,  19,  19, &
       18,  17,  16,  15,  14,  13,  11,  11,  10,   9,   8,   8, &
        7,   6,   5,   5,   4,   3,   2,   1,   0,   0,   0,   0 /)

  integer, parameter :: nn_x(nxlat) = &
    (/255, 171, 104,  60,  42,  32,  26,  20,  17,  15,  12,  11, &
        9,   9,   8,   7,   6,   6,   5,   4,   3,   3,   2,   1, &
        1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
        1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
        1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
        1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, &
        1,   2,   3,   3,   4,   5,   6,   6,   7,   8,   9,   9, &
       11,  12,  15,  17,  20,  26,  32,  42,  60, 104, 171, 255 /)

  ! for ring filter
  integer :: nlat_filter
  integer, allocatable :: chunk_array(:)
       
  contains
!-----------------------------------------------------------------------
  subroutine filter_init(use_ringfilter)
     use interpolate_data, only : lininterp

     logical, intent(in) :: use_ringfilter

     real(r8) :: xlats(nxlat), lats(nlat), kut1_out(nlat),  kut2_out(nlat)
     real(r8) :: nn_out(nlat)
     integer :: i
     character(len=128) :: errmsg

     if (use_ringfilter) then

        select case (nlon)
        case (72)
           ! 5 deg
           nlat_filter = 6
           allocate(chunk_array(nlat_filter))
           chunk_array(:) = &
                (/9,18,36,36,72,72/)
        case (144)
           ! 2.5deg
           nlat_filter = 16
           allocate(chunk_array(nlat_filter))
           chunk_array(:) = &
                (/9,9,9,18,18,18,36,36,36,36,36,72,72,72,72,72/)
        case (288)
           ! 1.25 deg 
           nlat_filter = 26
           allocate(chunk_array(nlat_filter))
           chunk_array(:) = &
                (/9,9,9,18,18,18,36,36,36,36,36,36,72,72,72,72,72,72,144,144,144,144, &
                144,144,144,144/)    
        case (576)
           ! 0.625 deg
           nlat_filter = 55
           allocate(chunk_array(nlat_filter))
           chunk_array(:) = &
                (/9,9,9,9,9,18,18,18,18,18,36,36,36,36,36, &
                72,72,72,72,72,72,72,72,72,72, &
                144,144,144,144,144,144,144,144,144,144, &
                288,288,288,288,288,288,288,288,288,288/)
        case (1152)
           ! 0.3125 deg
           nlat_filter = 111
           allocate(chunk_array(nlat_filter))
           chunk_array(:) = &
                (/9,9,9,9,9,9,9,9,9, &
                18,18,18,18,18,18,18,18,18, &
                36,36,36,36,36,36,36,36,36, &
                72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72,72, &
                144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144,144, &
                288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288,288, &
                576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576,576/)
        case default
           write (errmsg,'(a,i6)') 'ringfilter -- resolution is not recognized: number of global longitudes = ',nlon
           call endrun(errmsg)
        end select
        
     else
        
        ntrigs = 3*nlon/2+1
        allocate(trigs(ntrigs),kut1(nlat),kut2(nlat),nn(nlat))

        call set99(trigs,ifax,nlon) ! initialize fft for O+ polar filtering

        do i = 1,nlat
           lats(i) = -90._r8 + (i-1)*(180._r8/(nlat-1))
        enddo

        do i = 1,nxlat
           xlats(i) = -90._r8 + (i-1)*(180._r8/(nxlat-1))
        end do

        call lininterp( dble(kut1_x), xlats, nxlat, kut1_out, lats, nlat )
        call lininterp( dble(kut2_x), xlats, nxlat, kut2_out, lats, nlat )

        kut1 = int( kut1_out )
        kut2 = int( kut2_out )

        call lininterp( dble(nn_x), xlats, nxlat, nn_out, lats, nlat )
        nn = int( nn_out )

     endif

  end subroutine filter_init
!-----------------------------------------------------------------------
  subroutine filter1(f,lev0,lev1,lat)
!
! Remove longitudinal waves of prognostic variables with global fft.
! Remove wave numbers greater than kut(nlat).  This is called after 
! mp_gatherlons, and only by tasks with mytidi==0. On entry, task must 
! have global longitude data defined (mp_gatherlons).
!
! Args:
  integer,intent(in) :: lev0,lev1,lat
  real(r8),intent(inout) :: f(nlon,lev0:lev1)
!
! Local:
  integer :: n1,n2,k,i,nlevs
  real(r8) :: fx(nlon+2,lev1-lev0+1), wfft(nlon+1,lev1-lev0+1)
 
  nlevs = lev1-lev0+1
  n1 = 2*kut1(lat)+3 ! nyquist freq (?)
  n2 = nlon+2
  if (n1 > n2) then
    write(iulog,"('filter1: lat=',i2,' kutj=',i2,' n1,2=',2i3,' n1 > n2')") &
      lat,kut1(lat),n1,n2
    return
  endif
!
! Load fx from f for the fft:
  fx(:,:) = 0._r8
  do k=lev0,lev1
    do i=1,nlon
      fx(i,k) = f(i,k)
    enddo
  enddo

  call fft991(fx, wfft, trigs, ifax,1,nlon+2,nlon,nlevs,-1) 
!
! Remove wave numbers greater than kut(lat)
  do k = 1,nlevs
    do i=n1,n2
      fx(i,k) = 0.0_r8
    enddo
  enddo
!
! Inverse transform fourier back to gridpoint:
!
  call fft991(fx, wfft, trigs, ifax,1,nlon+2,nlon,nlevs,1)
!
! Redefine f from fx:
  do k=lev0,lev1
    do i=1,nlon
      f(i,k) = fx(i,k)
    enddo
  enddo
  end subroutine filter1
!-----------------------------------------------------------------------
  subroutine filter2(f,lev0,lev1,lat)
    use edyn_geogrid,only : dlamda
!
! Remove longitudinal waves of prognostic variables with global fft.
! Remove wave numbers greater than kut2(nlat).  This is called after 
! mp_gatherlons, and only by tasks with mytidi==0. On entry, task must 
! have global longitude data defined (mp_gatherlons).
!
! Args:
    integer,intent(in) :: lev0,lev1,lat
    real(r8),intent(inout) :: f(nlon,lev0:lev1)
!
! Local:
    integer :: n1,k,i,nlevs
    real(r8) :: fx(nlon+2,lev1-lev0+1), wfft(nlon+1,lev1-lev0+1)
    real(r8) :: smoothfunc,coslon
!
    nlevs = lev1-lev0+1
!
! Load local fx from inout f subdomain for the fft:
!
    fx(:,:) = 0._r8
    do k=lev0,lev1
      do i=1,nlon
        fx(i,k) = f(i,k)
      enddo
    enddo

    call fft991(fx, wfft, trigs, ifax,1,nlon+2,nlon,nlevs,-1)
!
! Wenbin's comments from TIEGCM:
! Change filters so that it does not over filtering at high latitudes, it will be the
! same as filter for low wavenumer, but wrapping up smoothly for large wavenumers, not a
! sharp transition, so there is still filtering effect in the lower latitudes
!                                                       Wenbin Wang  06/11/13
    n1=2*(kut2(lat)-1)+3
!
! Multiply by smoothing function:
! Test coslon to avoid underflow in smoothfunc at the poles
!

    do k=lev0,lev1
       do i=n1,nlon
          coslon = cos(((i-n1)/2._r8)*dlamda/2._r8)
          if ((coslon<0.1_r8) .or. (coslon<0.9_r8 .and. nn(lat)>50)) then
             fx(i,k) = 0._r8
          else
             smoothfunc = coslon**(2*nn(lat))
             fx(i,k) = fx(i,k)*smoothfunc
          endif
       enddo ! i=1,nlon
    enddo ! k=lev0,lev1
!
! Inverse transform fourier back to gridpoint:
!
    call fft991(fx, wfft, trigs, ifax,1,nlon+2,nlon,nlevs,1)
!
! Redefine f from fx:
    do k=lev0,lev1
      do i=1,nlon
        f(i,k) = fx(i,k)
      enddo
    enddo
  end subroutine filter2
!-----------------------------------------------------------------------

  subroutine ringfilter(lev0,lev1,lat,f)
!
! Ringfilter for the second order of FFT
! keep first and second order of fourier series, and filter orders
! coded by Dang, 2017
! Args:
    integer,intent(in) :: lat,lev0,lev1
    real(r8),intent(inout) :: f(nlon,lev0:lev1)
!
! Local:
    real(r8) :: fx(nlon),average(200),f_out(nlon)
    real(r8) :: w(nlon),wm1(nlon),a0,a1,b1,theta,dtheta
    real(r8) :: fL,fR,fm2,fm1,ff,fp1,fp2,a,b,c
    integer :: i,k,points,n,chunk,j_index,m

    near_pole: if(lat .LE. nlat_filter .OR. lat .GE. (nlat-nlat_filter+1) ) then

      dtheta=2._r8*3.14159_r8/real(nlon)

      do k=lev0,lev1

! Load field data into w
! Fourier expansion: f(x)=a0+a1*cos(x)+b1*sin(x)+others
        a1=0._r8
        b1=0._r8
        do i=1,nlon
           w(i) = f(i,k)
           theta=dtheta*i
           a1=a1+w(i)*cos(theta)
           b1=b1+w(i)*sin(theta)
        enddo
        a1=a1*2._r8/real(nlon)
        b1=b1*2._r8/real(nlon)
        a0=sum(w)/real(nlon)

! Chunk numbers in this latitude
        if(lat .LE. nlat_filter) chunk=chunk_array(lat)
        if(lat .GE. (nlat-nlat_filter+1)) chunk=chunk_array(nlat-lat+1)

! w(i)=wm1(i)+fx(i), then filter fx(i)
        do i=1,nlon
          theta=dtheta*i
          wm1(i)=a0+a1*cos(theta)+b1*sin(theta)
!          fx(i)=w(i)-wm1(i)      ! This option is to apply the ring filter to perturbations with wavenumber larger than 1
          fx(i)=w(i)             ! This is to apply the ring filter to the full field 
        enddo
       
! Start the ring average filtering

! Grid points in each chunk
        points=nlon/chunk
        n=points

! Calculate the average value in each chunk
        do i=1,chunk    ! i is the chunk number in each ring
          average(i)=sum(fx((i-1)*points+1:i*points))/real(points)
        enddo
        
! Then do the linear interpolation between each fL, fR
        do i=1,chunk  ! i is the chunk number in each ring

! Calculate f,fL,fR 
          if(i == 1) then
             
             fm2 = average(chunk-1)
             fm1 = average(chunk)
             ff  = average(i)
             fp1 = average(i+1)
             fp2 = average(i+2)

          else if(i == 2) then

             fm2 = average(chunk)
             fm1 = average(i-1)
             ff  = average(i)
             fp1 = average(i+1)
             fp2 = average(i+2)

          else if(i == chunk-1) then

             fm2 = average(i-2)
             fm1 = average(i-1)
             ff  = average(i)
             fp1 = average(i+1)
             fp2 = average(1)

          else if(i == chunk) then

             fm2 = average(i-2)
             fm1 = average(i-1)
             ff  = average(i)
             fp1 = average(1)
             fp2 = average(2)

          else

             fm2 = average(i-2)
             fm1 = average(i-1)
             ff  = average(i)
             fp1 = average(i+1)
             fp2 = average(i+2)

          endif

          fL = (-fm2+7._r8*fm1+7._r8*ff-fp1)/12._r8
          fR = (-fm1+7._r8*ff+7._r8*fp1-fp2)/12._r8

          a = 3._r8*(fL + fR - 2._r8*ff)
          b = 2._r8*(3._r8*ff - fR - 2._r8*fL)
          c = fL

! Calculate the filtered data at j_index
          do m=1,n
            j_index=m+(i-1)*points
            f_out(j_index)=(a/3.0_r8)*(3*m*m-3*m+1)/(n*n) &
             + 0.5_r8*b*(2*m-1)/n + c
          enddo

        enddo ! i=1,chunk

        fx(:)=f_out(:)

! Save filtered field:
        do i=1,nlon
          f(i,k) = fx(i)                ! This corresponds to the option applying ring filter to the full field.
        enddo ! i=1,nlon

      enddo ! k=lev0,lev1

    endif near_pole
  end subroutine ringfilter
!-----------------------------------------------------------------------
end module filter_module
