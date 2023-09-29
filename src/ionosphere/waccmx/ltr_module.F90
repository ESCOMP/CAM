module ltr_module
  !
  ! Module used to read data from the LFM/LTR outputs (POT,mean energy,
  !   and energy flux).
  !

  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc
  use edyn_maggrid,   only: nmlat, nmlonp1
  use edyn_maggrid,   only: ylonm     ! magnetic latitudes (nmlat) (radians)
  use edyn_maggrid,   only: ylatm     ! magnetic longtitudes (nmlonp1) (radians)
  use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile
  use pio,            only: pio_inq_dimid, pio_inquire_dimension
  use pio,            only: pio_inquire, pio_inq_varid
  use pio,            only: file_desc_t, pio_noerr, pio_nowrite, pio_get_var
  use utils_mod,      only: check_ncerr, check_alloc
  use edyn_mpi,       only: ntask, mytid
  use edyn_params,    only: pi, dtr, rtd
  use input_data_utils, only: time_coordinate

  implicit none

  private
  public :: init_ltr
  public :: getltr

  ! Grid dimension sizes for LTR input data file:
  integer :: lonp1,latp1
  !
  ! Define fields for LTR input data file:
  ! electric potential in Volt
  ! mean energy in KeV
  ! energy flux in W/m^2
  ! Time interpolated LTR outputs with suffix _ltr
  !
  real(r8),allocatable,dimension(:,:,:) :: & ! (lonp1,latp1,ntimes)
       pot_input, ekv_input, efx_input
  real(r8),allocatable,dimension(:,:) :: & ! (lonp1,latp1)
       pot_ltr, ekv_ltr, efx_ltr
  integer, allocatable,dimension(:) :: & ! (ntimes)
       year,month,day,jday
  real(r8), allocatable,dimension(:) :: & ! (ntimes)
       hpi_input, pcp_input, ltr_ut
  real(r8) :: hpi_ltr, pcp_ltr
  !
  type(file_desc_t) :: ncid

  character(len=cl), allocatable :: ltr_files(:)
  integer :: num_files, file_ndx

  type(time_coordinate) :: time_coord

contains

  !-----------------------------------------------------------------------
  subroutine init_ltr(ltr_list)

    character(len=*),intent(in) :: ltr_list(:)

    integer :: n, nfiles

    nfiles = size(ltr_list)
    num_files = 0

    count_files: do n = 1,nfiles
       if (len_trim(ltr_list(n))<1 .or. &
           trim(ltr_list(n))=='NONE') then
          exit count_files
       else
          num_files = num_files + 1
       end if
    end do count_files

    allocate(ltr_files(num_files))
    ltr_files(:num_files) = ltr_list(:num_files)
    file_ndx = 1
    call open_files()

  end subroutine init_ltr

  !-----------------------------------------------------------------------
  subroutine rdltr(ltrfile)
    !
    ! Read LTR data
    !
    character(len=*), intent(in) :: ltrfile
    ! Local:
    integer                      :: istat, ntimes, ndims, nvars, ngatts, ier
    integer                      :: idunlim
    integer                      :: id_lon, id_lat, id_time
    integer                      :: idv_year, idv_mon, idv_day, idv_jday
    integer                      :: idv_ut, idv_hpi, idv_pcp
    character(len=*), parameter  :: subname = 'rdltr'

    !
    !
    if (masterproc) then
       write(iulog, "(/, 72('-'))")
       write(iulog, "(a, ': read LTR data:')") subname
    end if
    !
    ! Open netcdf file:
    call cam_pio_openfile(ncid, ltrfile, pio_nowrite)
    !
    ! Get LTR grid dimension:
    istat = pio_inq_dimid(ncid, 'lon', id_lon)
    istat = pio_inquire_dimension(ncid, id_lon, len=lonp1)
    call check_ncerr(istat, subname, 'LTR longitude dimension')

    istat = pio_inq_dimid(ncid, 'lat', id_lat)
    istat = pio_inquire_dimension(ncid, id_lat, len=latp1)
    call check_ncerr(istat, subname, 'LTR latitude dimension')

    call time_coord%initialize( ltrfile, set_weights=.false. )
    !
    ! Get time dimension:
    istat = pio_inquire(ncid, unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid, id_time, len=ntimes)
    call check_ncerr(istat, subname, 'LTR time dimension')
    !
    ! Search for requested LTR output fields
    istat = pio_inquire(ncid,ndims,nvars,ngatts,idunlim)
    !
    ! Get 1-D LTR fields (ntimes)
    if (.not. allocated(year)) then
       allocate(year(ntimes), stat=ier)
       call check_alloc(ier, subname, 'year', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'year', idv_year)
    call check_ncerr(istat, subname, 'LTR year id')
    istat = pio_get_var(ncid, idv_year, year)
    call check_ncerr(istat, subname, 'LTR year')

    if (.not. allocated(month)) then
       allocate(month(ntimes), stat=ier)
       call check_alloc(ier, subname, 'month', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'month', idv_mon)
    call check_ncerr(istat, subname, 'LTR month id')
    istat = pio_get_var(ncid, idv_mon, month)
    call check_ncerr(istat, subname, 'LTR month')
    if (.not. allocated(day)) then
       allocate(day(ntimes), stat=ier)
       call check_alloc(ier, subname, 'day', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'day', idv_day)
    call check_ncerr(istat, subname, 'LTR day id')
    istat = pio_get_var(ncid, idv_day, day)
    call check_ncerr(istat, subname, 'LTR day')

    if (.not. allocated(jday)) then
       allocate(jday(ntimes), stat=ier)
       call check_alloc(ier, subname, 'jday', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid, 'jday', idv_jday)
    call check_ncerr(istat, subname, 'LTR jday id')
    istat = pio_get_var(ncid, idv_jday, jday)
    call check_ncerr(istat, subname, 'LTR jday')
    !
    ! Allocate 1-d fields:
    if (.not. allocated(ltr_ut)) then
       allocate(ltr_ut(ntimes), stat=ier)
       call check_alloc(ier, subname, 'ltr_ut', ntimes=ntimes)
    end if
    if (.not. allocated(hpi_input)) then
       allocate(hpi_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'hpi_input', ntimes=ntimes)
    end if
    if (.not. allocated(pcp_input)) then
       allocate(pcp_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'pcp_input', ntimes=ntimes)
    end if
    !
    ! Get ut
    istat = pio_inq_varid(ncid, 'ut', idv_ut)
    call check_ncerr(istat, subname, 'LTR ut id')
    istat = pio_get_var(ncid, idv_ut, ltr_ut)
    call check_ncerr(istat, subname, 'LTR ut')
    !
    ! Get HPI
    istat = pio_inq_varid(ncid, 'hpiN', idv_hpi)
    call check_ncerr(istat, subname, 'LTR hpi id')
    istat = pio_get_var(ncid, idv_hpi, hpi_input)
    call check_ncerr(istat, subname, 'LTR hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid, 'pcpN', idv_pcp)
    call check_ncerr(istat, subname, 'LTR pcp id')
    istat = pio_get_var(ncid, idv_pcp, pcp_input)
    call check_ncerr(istat, subname, 'LTR pcp')
    !
    ! Allocate 2-d fields:
    if (.not. allocated(pot_ltr)) then
       allocate(pot_ltr(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'pot_ltr', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(ekv_ltr)) then
       allocate(ekv_ltr(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'ekv_ltr', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(efx_ltr)) then
       allocate(efx_ltr(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'efx_ltr', lonp1=lonp1, latp1=latp1)
    end if
    !
    ! Allocate 3-d fields:
    if (.not. allocated(pot_input)) then
       allocate(pot_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'pot_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(ekv_input)) then
       allocate(ekv_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'ekv_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(efx_input)) then
       allocate(efx_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'efx_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
  end subroutine rdltr

  !-----------------------------------------------------------------------
  subroutine update_3d_fields( ncid, offset, kount, pot_3d,ekv_3d,efx_3d )

    type(file_desc_t), intent(in) :: ncid
    integer, intent(in)  :: offset(:)
    integer, intent(in)  :: kount(:)
    real(r8),intent(out) :: pot_3d(:,:,:)
    real(r8),intent(out) :: ekv_3d(:,:,:)
    real(r8),intent(out) :: efx_3d(:,:,:)


    integer :: istat
    integer :: idv_pot,idv_ekv, idv_efx
    character(len=*), parameter :: subname = 'update_3d_fields'

    !
    ! Get 3-D fields (lon,lat,ntimes)
    !
    ! electric potential
    istat = pio_inq_varid(ncid, 'pot', idv_pot)
    call check_ncerr(istat, subname, 'LTR pot id')
    istat = pio_get_var(ncid, idv_pot, offset, kount, pot_3d)
    call check_ncerr(istat, subname, 'LTR pot')
    !
    ! mean energy
    istat = pio_inq_varid(ncid, 'ekv', idv_ekv)
    call check_ncerr(istat, subname, 'LTR ekv id')
    istat = pio_get_var(ncid, idv_ekv, offset, kount, ekv_3d)
    call check_ncerr(istat, subname, 'LTR ekv')
    !
    ! energy flux
    istat = pio_inq_varid(ncid, 'efx', idv_efx)
    call check_ncerr(istat, subname, 'LTR efx id')
    istat = pio_get_var(ncid, idv_efx, offset, kount, efx_3d)
    call check_ncerr(istat, subname, 'LTR efx')

  end subroutine update_3d_fields

  !-----------------------------------------------------------------------
  subroutine getltr(iyear, imo, iday, iutsec, sunlon, iprint,  &
                     iltr, phihm, ltr_efxm, ltr_kevm)
    use cam_history_support, only: fillvalue
    use rgrd_mod, only: rgrd2
    !
    !    Read LTR outputs from ltr_ncfile file, returning electric potential,
    !    auroral mean energy and energy flux at current date and time,
    !    and the data is linearly interpolated to the model time
    !
    !
    !    Args:

    integer,  intent(in)    :: iyear
    integer,  intent(in)    :: imo
    integer,  intent(in)    :: iday
    real(r8), intent(in)    :: sunlon
    integer,  intent(in)    :: iutsec
    integer,  intent(in)    :: iprint
    integer,  intent(out)   :: iltr
    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: ltr_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: ltr_kevm(nmlonp1,nmlat) ! on geomag grid
    !
    !
    !     Local:
    real(r8)                    :: potm(lonp1,latp1)
    real(r8)                    :: efxm(lonp1,latp1), ekvm(lonp1,latp1)
    real(r8)                    :: alat(latp1), alon(lonp1)
    real(r8)                    :: alatm(latp1), alonm(lonp1)
    integer                     :: ier, lw, liw, intpol(2)
    integer,  allocatable       :: iw(:)
    real(r8), allocatable       :: w(:)
    integer                     :: i, j, ithmx
    integer                     :: nn, iset2, iset1, m, mp1, n
    real(r8)                    :: f1, f2
    real(r8)                    :: del, xmlt, dmlat, dlatm, dlonm, dmltm, rot
    integer                     :: offset(3), kount(3)
    character(len=*), parameter :: subname = 'getltr'

    phihm = fillvalue
    ltr_efxm = fillvalue
    ltr_kevm = fillvalue

    if (iprint > 0 .and. masterproc) then
       write(iulog,"(/,72('-'))")
       write(iulog,"(a,':')") subname
       write(iulog,"(a,i4,', iday = ',i3,', iutsec = ',i10)")                &
            'Initial requested iyear= ', iyear, iday, iutsec
    end if

    nn = size(ltr_ut)

    !
    !     Check times:
    !
    iltr = 1 - time_coord%times_check()

    check_loop: do while( iltr/=1 )

       if (masterproc) write(iulog,*) 'file_ndx = ',file_ndx

       if (iltr==2) then
          if (masterproc) then
             write(iulog, "(a,': Model date prior to LTR first date:',3I5)")   &
                  subname, year(1), month(1), day(1)
          end if
          return
       endif

       if (iltr==0) then
          if (masterproc) then
             write(iulog, "(a,': Model date beyond the LTR last Data:',3I5)")  &
                  subname, year(nn), month(nn), day(nn)
          end if

          if (file_ndx<num_files) then
             file_ndx = file_ndx+1
             call close_files()
             call open_files()
             iltr = 1 - time_coord%times_check()
          else
             return
          end if

       endif
    end do check_loop

    !
    !     get LTR data
    pot_ltr(:,:) = 0._r8
    ekv_ltr(:,:) = 0._r8
    efx_ltr(:,:) = 0._r8
    hpi_ltr = 0._r8
    pcp_ltr = 0._r8

    call time_coord%advance()

    iset1 = time_coord%indxs(1)
    iset2 = time_coord%indxs(2)

    f1 = time_coord%wghts(1)
    f2 = time_coord%wghts(2)

    hpi_ltr = (f1*hpi_input(iset1) + f2*hpi_input(iset2))
    pcp_ltr = (f1*pcp_input(iset1) + f2*pcp_input(iset2))

    offset = (/1,1,iset1/)
    kount = (/lonp1,latp1,2/)
    call update_3d_fields( ncid, offset, kount, pot_input,ekv_input,efx_input )
    pot_ltr(:,:) = (f1*pot_input(:,:,1) + f2*pot_input(:,:,2))
    ekv_ltr(:,:) = (f1*ekv_input(:,:,1) + f2*ekv_input(:,:,2))
    efx_ltr(:,:) = (f1*efx_input(:,:,1) + f2*efx_input(:,:,2))

    active_task: if ( mytid<ntask ) then

       !     mlongitude starts from 180 degree
       rot = sunlon*rtd
       if(rot < 0) then
          rot = rot + 360._r8    !  0 to 360 degrees
       end if
       rot = rot / 15._r8        !  convert from degree to hrs

       dmltm = 24._r8 / real(lonp1, kind=r8)

       do i = 1, lonp1
          xmlt = (real(i-1, kind=r8) * dmltm) - rot + 24._r8
          xmlt = MOD(xmlt, 24._r8)
          m = int(xmlt/dmltm + 1.001_r8)
          mp1 = m + 1
          if (mp1 > lonp1) mp1 = 2
          del = xmlt - (m-1)*dmltm
          !     Put in LTR arrays from south pole to north pole
          do j=1,latp1
             potm(i,j) = (1._r8-del)*pot_ltr(m,j) + &
                  del*pot_ltr(mp1,j)
             ekvm(i,j) = (1._r8-del)*ekv_ltr(m,j) + &
                  del*ekv_ltr(mp1,j)
             if (ekvm(i,j) == 0._r8) ekvm(i,j)=1._r8
             efxm(i,j) = (1._r8-del)*efx_ltr(m,j) + &
                  del*efx_ltr(mp1,j)
          end do

       end do

       !     Set up coeffs to go between EPOTM(IMXMP,JMNH) and TIEPOT(IMAXM,JMAXMH)

       !     ****     SET GRID SPACING DLATM, DLONG, DLONM
       !     DMLAT=lat spacing in degrees of LTR apex grid
       dmlat = 180._r8 / real(latp1-1, kind=r8)
       dlatm = dmlat * dtr
       dlonm = 2._r8 * pi / real(lonp1, kind=r8)
       dmltm = 24._r8 / real(lonp1, kind=r8)
       !     ****
       !     ****     SET ARRAY YLATM (LATITUDE VALUES FOR GEOMAGNETIC GRID
       !     ****
       alatm(1) = -pi / 2._r8
       alat(1) = -90._r8
       alatm(latp1) = pi / 2._r8
       alat(latp1) = 90._r8
       ithmx = (latp1+1)/2
       do i = 2, ithmx
          alat(i) = alat(i-1)+dlatm*rtd
          alat(latp1+1-i) = alat(latp1+2-i)-dlatm*rtd
          alatm(i) = alatm(i-1)+dlatm
          alatm(latp1+1-i) = alatm(latp1+2-i)-dlatm
       end do
       alon(1) = -pi*rtd
       alonm(1) = -pi
       do i = 2, lonp1
          alon(i) = alon(i-1) + dlonm*rtd
          alonm(i) = alonm(i-1) + dlonm
       end do

       !     ylatm and ylonm are arrays of latitudes and longitudes of the
       !     distorted magnetic grids in radian - from consdyn.h
       !     Convert from apex magnetic grid to distorted magnetic grid
       !
       !     Allocate workspace for regrid routine rgrd_mod:
       lw = nmlonp1+nmlat+2*nmlonp1
       if (.not. allocated(w)) then
          allocate(w(lw), stat=ier)
          call check_alloc(ier, 'getltr', 'w', lw=lw)
       end if
       liw = nmlonp1 + nmlat
       if (.not. allocated(iw)) then
          allocate(iw(liw), stat=ier)
          call check_alloc(ier, 'getltr', 'iw', lw=liw)
       end if
       intpol(:) = 1             ! linear (not cubic) interp in both dimensions
       if (alatm(1) > ylatm(1)) then
          alatm(1) = ylatm(1)
       end if
       if (alatm(latp1) < ylatm(nmlat)) then
          alatm(latp1) = ylatm(nmlat)
       end if
       if (alonm(1) > ylonm(1)) then
          alonm(1) = ylonm(1)
       end if
       if (alonm(lonp1) < ylonm(nmlonp1)) then
          alonm(lonp1) = ylonm(nmlonp1)
       end if

       !     ylatm from -pi/2 to pi/2, and ylonm from -pi to pi
       call rgrd2(lonp1, latp1, alonm, alatm, potm, nmlonp1, nmlat,  &
            ylonm, ylatm, phihm, intpol, w, lw, iw, liw, ier)
       call rgrd2(lonp1, latp1, alonm, alatm, ekvm, nmlonp1, nmlat,  &
            ylonm, ylatm, ltr_kevm, intpol, w, lw, iw, liw, ier)
       call rgrd2(lonp1, latp1, alonm, alatm, efxm, nmlonp1, nmlat,  &
            ylonm, ylatm, ltr_efxm, intpol, w, lw, iw, liw, ier)

       if (iprint > 0 .and. masterproc) then
          write(iulog, *) subname, ': Max, min ltr_efxm = ', &
               maxval(ltr_efxm), minval(ltr_efxm)
          write(iulog, "('getltr: LTR data interpolated to date and time')")
          write(iulog,"('getltr: iyear,imo,iday,iutsec = ',3i6,i10)")  &
               iyear,imo,iday,iutsec
          write(iulog,"('getltr: LTR iset1 f1,f2,year,mon,day,ut = ',  &
               i6,2F9.5,3I6,f10.4)")  &
               iset1,f1,f2,year(iset1),month(iset1),day(iset1),ltr_ut(iset1)
          write(iulog,*)'getltr: max,min phihm= ', maxval(phihm),minval(phihm)
       end if

    end if active_task

  end subroutine getltr
  !-------------------------------------------------------------------

  subroutine close_files

    deallocate( year,month,day )
    deallocate( hpi_input, pcp_input, ltr_ut )

    call cam_pio_closefile(ncid)

  end subroutine close_files
  !-----------------------------------------------------------------------
  subroutine open_files()

    call rdltr(ltr_files(file_ndx))

  end subroutine open_files

end module ltr_module
