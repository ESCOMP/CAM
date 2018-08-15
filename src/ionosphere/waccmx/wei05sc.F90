module wei05sc
!
! The Weimer model of high-latitude potential created by Daniel Weimer and 
! if extracted, distributed, or used for any purpose other than as implemented 
! in the NCAR TIEGCM and CESM/WACCM models, please contact Dan Weimer for 
! further information and discussion.
!
! 2005 Version of the electric and magnetic potential (FAC) models
! by Dan Weimer.  Uses Spherical Cap Harmonic Analysis (SCHA) functions.
! Model description is in:
! Weimer, D. R., Predicting Surface Geomagnetic Variations Using Ionospheric
!    Electrodynamic Models, J. Geophys. Res., 110, A12307, doi:10.1029/
!    2005JA011270, 2005.
! Some information about the model (such as outer boundary calculation)
! is also in the earlier paper:
! Weimer, D. R. (2005), Improved ionospheric electrodynamic models and
!    application to calculating Joule heating rates,  J. Geophys. Res., 110,
!    A05306, doi:10.1029/2004JA010884.
!
! For information about the SCHA, see the paper:
! Haines, G. V., Spherical cap harmonic analysis, J. Geophys. Res., 90, B3,
!  2583, 1985.  (Note that this is in JGR-B, "Solid Earth", rather than JGR-A)
!
! April, 2008:
! This f90 module of the Electric Potential model was translated
!   from the original IDL by Ben Foster (NCAR, foster@ucar.edu)
! Netcdf data file wei05sc.nc was written from original IDL save files
!   W05scBndy.xdr, W05scEpot.xdr, W05scBpot.xdr, and SCHAtable.dat
!
! September, 2015 btf:
! Modified for free-format fortran, and for CESM/WACCM (r8, etc).
!
  use shr_kind_mod  ,only: r8 => shr_kind_r8
  use shr_kind_mod  ,only: shr_kind_cl
  use spmd_utils    ,only: masterproc
#ifdef WACCMX_IONOS
  use cam_logfile   ,only: iulog
  use cam_abortutils,only: endrun  
  use time_manager  ,only: get_curr_date
  use edyn_maggrid  ,only: nmlat,nmlon,nmlonp1
#endif

  use edyn_maggrid,only: &
    ylonm,    & ! magnetic latitudes (nmlat) (radians)
    ylatm       ! magnetic longtitudes (nmlonp1) (radians)
  use edyn_solve,only: &
    nmlat0,   & ! (nmlat+1)/2
    phihm       ! output: high-latitude potential (nmlonp1,nmlat)

  use physconst, only: pi
  use aurora_params, only: aurora_params_set, hpower, ctpoten, theta0
  use aurora_params, only: offa, dskofa, dskofc, phid, rrad, offc, phin
  implicit none
  private

#ifdef WACCMX_IONOS
!
! Coefficients read from netcdf data file wei05sc.nc:
!
  integer,parameter :: &
    na=6, nb=7, nex=2, n1_scha=19, n2_scha=7, n3_scha=68, &
    csize=28, n_schfits=15, n_alschfits=18
  integer :: maxk_scha, maxm_scha, maxl_pot, maxm_pot
  real(r8) :: bndya(na), bndyb(nb), ex_bndy(nex), ex_epot(nex),ex_bpot(nex)
  real(r8) :: th0s(n3_scha), allnkm(n1_scha,n2_scha,n3_scha)
  integer :: ab(csize), ls(csize), ms(csize)
  real(r8) :: epot_alschfits(n_alschfits,csize), bpot_alschfits(n_alschfits,csize)
  real(r8) :: bpot_schfits(n_schfits,csize),epot_schfits(n_schfits,csize)
!
! Intermediate calculations:
!
  integer,parameter :: mxtablesize=500
  real(r8) :: rad2deg,deg2rad           ! set by setmodel
  real(r8) :: bndyfitr                  ! calculated by setboundary
  real(r8) :: esphc(csize),bsphc(csize) ! calculated by setmodel
  real(r8) :: tmat(3,3) !,ttmat(3,3)      ! from setboundary
  real(r8) :: plmtable(mxtablesize,csize),colattable(mxtablesize)
  real(r8) :: nlms(csize)
  real(r8) :: wei05sc_fac(nmlonp1,nmlat)  ! field-aligned current output

! 05/08 bae:  Have ctpoten from both hemispheres from Weimer
  real(r8) :: weictpoten(2),phimin,phimax

  real(r8) :: real8,real8a ! for type conversion to 8-byte real

!
! Several items in the public list are for efield.F90 (chemistry/mozart)
! (dpie_coupling calls the weimer05 driver, but efield calls the individual
!  routines, not the driver)
!
  public :: weimer05
  public :: weimer05_init

#endif

  real(r8), parameter :: r2d = 180._r8/pi  ! radians to degrees
  real(r8), parameter :: d2r = pi/180._r8  ! degrees to radians

  contains

!-----------------------------------------------------------------------
  subroutine weimer05_init(wei05_ncfile)
    use infnan, only: nan, assignment(=)

    character(len=*),intent(in) :: wei05_ncfile

    hpower = nan
    ctpoten = nan
    phin = nan
    phid = nan
    theta0 = nan
    offa = nan
    dskofa = nan
    rrad = nan
    offc = nan
    dskofc = nan
    
    bndya = nan
    bndyb = nan
    ex_bndy = nan
    ex_bpot = nan
    th0s = nan
    allnkm = nan
    bpot_schfits = nan
    bpot_alschfits = nan

    if (wei05_ncfile.ne.'NONE') then
       call read_wei05_ncfile(wei05_ncfile)
       aurora_params_set = .true.
    endif

  end subroutine weimer05_init

!-----------------------------------------------------------------------
  subroutine weimer05(by,bz_in,swvel,swden,sunlons)
!
! 9/16/15 btf: Driver to call Weimer 2005 model for waccm[x].
!

  implicit none
!
! Args:
  real(r8),intent(in) :: bz_in,by,swvel,swden
  real(r8),intent(in) :: sunlons(:)

#ifdef WACCMX_IONOS
!
! Local:

  real(r8) :: angl,angle,bt
  integer :: i,j
  real(r8) :: rmlt,mlat,tilt,htilt,hem,ut,secs
  real(r8),parameter :: fill=0._r8
  integer :: iyear,imon,iday,isecs
  logical :: debug = .false.
  real(r8) :: bz

  bz = bz_in

  hpower = hp_from_bz_swvel(bz,swvel)
!
! Get current date and time:
!
  call get_curr_date(iyear,imon,iday,isecs)
!
! Get sun's location (longitude at all latitudes):
!
  real8 = dble(isecs)
  secs = real8

!
! At least one of by,bz must be non-zero: 
  if (by==0._r8.and.bz==0._r8) then
     if (masterproc) then
        write(iulog,"(/,'>>> WARNING: by and bz cannot both be zero',&
                        ' when calling the Weimer model: am setting bz=0.01')")
     endif
     bz = 0.01_r8
  endif
!
  bt = sqrt(by**2+bz**2)
  angl = atan2(by,bz)*r2d
!
! Convert from day-of-year to month,day and get tilt from date and ut:
!
  ut = secs/3600._r8    ! decimal hours
!
! Given year and day-of-year, cvt2md returns month and day of month.
! We do not need this, since get_curr_date returns month and day of month.
! call cvt2md(iulog,iyear,idoy,imon,iday)  ! given iyear,idoy, return imo,ida
!
  if (debug) write(iulog,"('weimer05: iyear,imon,iday=',3i5,' ut=',f8.2)") &
    iyear,imon,iday,ut
  tilt = get_tilt(iyear,imon,iday,ut)
  if (debug) write(iulog,"('weimer05: tilt=',e12.4)") tilt

  phihm = 0._r8 ! whole-array init (nmlonp1,nmlat)
!
! Call Weimer model for southern hemisphere electric potential:
!
  hem = -1._r8
  htilt = hem * tilt
  angle = hem * angl
  if (debug) write(iulog,"('weimer05 call setmodel for SH potential')")
  call setmodel(angle,bt,htilt,swvel,swden,'epot')
  if (debug) write(iulog,"('weimer05 after setmodel for SH potential')")
  do j=1,nmlat0 ! Spole to equator
    do i=1,nmlon
!
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
!
      rmlt = (ylonm(i)-sunlons(1)) * r2d / 15._r8 + 12._r8
      mlat = abs(ylatm(j))*r2d
!
! Obtain electric potential and convert from kV to V
!
      call epotval(mlat,rmlt,fill,phihm(i,j))
      phihm(i,j) = phihm(i,j)*1000._r8
    enddo ! i=1,nmlon
  enddo ! j=1,nmlat0
  if (debug) write(iulog,"('weimer05: SH phihm min,max=',2es12.4)") &
    minval(phihm(1:nmlon,1:nmlat0)),maxval(phihm(1:nmlon,1:nmlat0))
!
! Re-calculate SH values of offa, dskofa, arad, and phid and phin from
!    Weimer 2005 setboundary values of offc, dskofc, and theta0
!
  call wei05loc (1, by, hpower, sunlons)
!
! Call Weimer model for southern hemisphere fac:
!
  if (debug) write(iulog,"('weimer05 call setmodel for SH fac')")
  call setmodel(angle,bt,htilt,swvel,swden,'bpot')
  if (debug) write(iulog,"('weimer05 after setmodel for SH fac')")
  do j=1,nmlat0
    do i=1,nmlon
      rmlt = (ylonm(i)-sunlons(1)) * r2d / 15._r8 + 12._r8
      mlat = abs(ylatm(j))*r2d
      call mpfac(mlat,rmlt,fill,wei05sc_fac(i,j))
    enddo ! i=1,nmlon
  enddo ! j=1,nmlat0
!
! Call Weimer model for northern hemisphere epot:
!
  hem = 1._r8
  htilt = hem * tilt
  angle = hem * angl
  if (debug) write(iulog,"('weimer05 call setmodel for NH potential')")
  call setmodel(angle,bt,htilt,swvel,swden,'epot')
  if (debug) write(iulog,"('weimer05 after setmodel for NH potential')")
  do j=nmlat0+1,nmlat
    do i=1,nmlon
!
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc) in rad
      rmlt = (ylonm(i)-sunlons(1)) * r2d / 15._r8 + 12._r8
      mlat = abs(ylatm(j))*r2d
!
! Obtain electric potential and convert from kV to V
      call epotval(mlat,rmlt,fill,phihm(i,j))
      phihm(i,j) = phihm(i,j)*1000._r8
    enddo ! i=1,nmlon
  enddo ! j=1,nmlat0+1,nmlat
  if (debug) write(iulog,"('weimer05: NH phihm min,max=',2es12.4)") &
    minval(phihm(1:nmlon,nmlat0+1:nmlat)),maxval(phihm(1:nmlon,nmlat0+1:nmlat))
!
! Re-calculate NH values of offa, dskofa, arad, and Heelis phid and phin from
!   Weimer 2005 setboundary values of offc, dskofc, and theta0
!
 call wei05loc (2, by, hpower, sunlons)
!
! Call Weimer model for northern hemisphere fac:
  if (debug) write(iulog,"('weimer05 call setmodel for NH fac')")
  call setmodel(angle,bt,htilt,swvel,swden,'bpot')
  if (debug) write(iulog,"('weimer05 after setmodel for NH fac')")
  do j=nmlat0+1,nmlat
    do i=1,nmlon
      rmlt = (ylonm(i)-sunlons(1)) * r2d / 15._r8 + 12._r8
      mlat = abs(ylatm(j))*r2d
      call mpfac(mlat,rmlt,fill,wei05sc_fac(i,j))
    enddo ! i=1,nmlon
  enddo ! j=1,nmlat0
!
! Periodic points:
  do j=1,nmlat
    phihm(nmlonp1,j) = phihm(1,j)
    wei05sc_fac(nmlonp1,j) = wei05sc_fac(1,j)
  enddo ! j=1,nmlat
!
! Calculate ctpoten for each hemisphere:
! South:
!
  phimax = -1.e36_r8
  phimin =  1.e36_r8
  do j=1,nmlat0 ! SH
    do i=1,nmlon
      if (phihm(i,j) > phimax) phimax = phihm(i,j)
      if (phihm(i,j) < phimin) phimin = phihm(i,j)
    enddo
  enddo 
  weictpoten(1) = 0.001_r8 * (phimax - phimin)
!
! North:
!
  phimax = -1.e36_r8
  phimin =  1.e36_r8
  do j=nmlat0+1,nmlat ! NH
    do i=1,nmlon
      if (phihm(i,j) > phimax) phimax = phihm(i,j)
      if (phihm(i,j) < phimin) phimin = phihm(i,j)
    enddo
  enddo 
  weictpoten(2) = 0.001_r8 * (phimax - phimin)
!
! average of the SH and NH in ctpoten
  ctpoten = 0.5_r8*(weictpoten(1)+weictpoten(2))

  if (masterproc) then
    write(iulog,"('weimer05: ctpoten=',f8.2,' phihm min,max=',2es12.4)") ctpoten,minval(phihm),maxval(phihm)
  endif
!

#endif
  end subroutine weimer05
!-----------------------------------------------------------------------
  subroutine read_wei05_ncfile(file)

    use ioFileMod,     only: getfil
    use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile
    use pio, only: file_desc_t, pio_nowrite, pio_inq_dimid, pio_inquire_dimension, &
                   pio_inq_varid, pio_get_var
!
! Read coefficients and other data from netcdf data file.
!
  implicit none
!
! Arg:
  character(len=*),intent(in) :: file
#ifdef WACCMX_IONOS
!
! Local:
  integer :: istat
  integer :: rd_na,rd_nb,rd_nex,rd_n1_scha,rd_n2_scha,rd_n3_scha,&
    rd_csize,rd_n_schfits,rd_n_alschfits
  integer :: id
  character(len=shr_kind_cl) :: filen
  type(file_desc_t) :: ncid
!
! Open netcdf file for reading:
!
  call getfil( file, filen, 0 )
  call cam_pio_openfile(ncid, filen, PIO_NOWRITE)
  
  write(iulog,"('wei05sc: opened netcdf data file',a)") trim(filen)
!
! Read and check dimensions:
!
! na=6
  istat = pio_inq_dimid(ncid,'na',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_na)
  if (rd_na /= na) then
     write(iulog,"(/,'>>> wei05sc: rd_na /= na: rd_na=',i4,' na=',i4)") rd_na,na
     call endrun('wei05sc: rd_na /= na')
  endif
!
! nb=7
!
  istat = pio_inq_dimid(ncid,'nb',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_nb)
  if (rd_na /= na) then
    write(iulog,"(/,'>>> wei05sc: rd_nb /= nb: rd_nb=',i4,' nb=',i4)") rd_nb,nb
    call endrun('wei05sc: rd_nb /= nb: rd_nb')
  endif
!
! nex=2
!
  istat = pio_inq_dimid(ncid,'nex',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_nex)
  if (rd_nex /= nex) then
    write(iulog,"(/,'>>> wei05sc: rd_nex /= nex: rd_nex=',i4,' nex=',i4)") &
      rd_nex,nex
    call endrun('wei05sc')
  endif
!
! n1_scha=19
!
  istat = pio_inq_dimid(ncid,'n1_scha',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_n1_scha)
  if (rd_n1_scha /= n1_scha) then
    write(iulog,"(/,'>>> wei05sc: rd_n1_scha /= n1_scha: rd_n1_scha=',i4,' n1_scha=',i4)") &
      rd_n1_scha,n1_scha
    call endrun('wei05sc')
  endif
!
! n2_scha=7
!
  istat = pio_inq_dimid(ncid,'n2_scha',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_n2_scha)
  if (rd_n2_scha /= n2_scha) then
    write(iulog,"(/,'>>> wei05sc: rd_n2_scha /= n2_scha: rd_n2_scha=',i4,' n2_scha=',i4)") &
      rd_n2_scha,n2_scha
    call endrun('wei05sc')
  endif
!
! n3_scha=68
!
  istat = pio_inq_dimid(ncid,'n3_scha',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_n3_scha)
  if (rd_n3_scha /= n3_scha) then
    write(6,"(/,'>>> wei05sc: rd_n3_scha /= n3_scha: rd_n3_scha=',i4,' n3_scha=',i4)") &
      rd_n3_scha,n3_scha
    call endrun('wei05sc')
  endif
!
! csize=28
!
  istat = pio_inq_dimid(ncid,'csize',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_csize)
  if (rd_csize /= csize) then
    write(iulog,"(/,'>>> wei05sc: rd_csize /= csize: rd_csize=',i4,' csize=',i4)") &
      rd_csize,csize
    call endrun('wei05sc')
  endif
!
! n_schfits=15
!
  istat = pio_inq_dimid(ncid,'n_schfits',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_n_schfits)
  if (rd_n_schfits /= n_schfits) then
    write(iulog,"(/,'>>> wei05sc: rd_n_schfits /= n_schfits: rd_n_schfits=',i4,' n_schfits=',i4)") &
      rd_n_schfits,n_schfits
    call endrun('wei05sc')
  endif
!
! n_alschfits=18
!
  istat = pio_inq_dimid(ncid,'n_alschfits',id)
  istat = pio_inquire_dimension(ncid,id,len=rd_n_alschfits)
  if (rd_n_alschfits /= n_alschfits) then
    write(iulog,"(/,'>>> wei05sc: rd_n_alschfits /= n_alschfits: rd_n_alschfits=',i4,' n_alschfits=',i4)") & 
      rd_n_alschfits,n_alschfits
    call endrun('wei05sc')
  endif
!
! integer :: maxk_scha, maxm_scha, maxl_pot, maxm_pot
! maxk_scha = 18 ;
! maxm_scha = 6 ;
! maxl_pot = 12 ;
! maxm_pot = 2 ;
!
  istat = pio_inq_dimid(ncid,"maxk_scha",id) 
  istat = pio_inquire_dimension(ncid,id,len=maxk_scha)
  istat = pio_inq_dimid(ncid,"maxm_scha",id) 
  istat = pio_inquire_dimension(ncid,id,len=maxm_scha)
  istat = pio_inq_dimid(ncid,"maxl_pot",id) 
  istat = pio_inquire_dimension(ncid,id,len=maxl_pot)
  istat = pio_inq_dimid(ncid,"maxm_pot",id) 
  istat = pio_inquire_dimension(ncid,id,len=maxm_pot)

! write(iulog,"('wei05sc: maxk_scha=',i3,' maxm_scha=',i3)") &
!   maxk_scha,maxm_scha
! write(iulog,"('wei05sc: maxl_pot=',i3,' maxm_pot=',i3)") &
!   maxl_pot,maxm_pot
!
! Read variables:
!
! double bndya(na):
  istat = pio_inq_varid(ncid,'bndya',id)
  istat = pio_get_var(ncid,id,bndya)
! write(iulog,"('wei05sc: bndya=',/,(8f8.3))") bndya
!
! double bndyb(nb):
  istat = pio_inq_varid(ncid,'bndyb',id)
  istat = pio_get_var(ncid,id,bndyb)
! write(iulog,"('wei05sc: bndyb=',/,(8f8.3))") bndyb
!
! double ex_bndy(nex):
  istat = pio_inq_varid(ncid,'ex_bndy',id)
  istat = pio_get_var(ncid,id,ex_bndy)
! write(iulog,"('wei05sc: ex_bndy=',/,(8f8.3))") ex_bndy
!
! double th0s(n3_scha):
  istat = pio_inq_varid(ncid,'th0s',id)
  istat = pio_get_var(ncid,id,th0s)
! write(iulog,"('wei05sc: th0s=',/,(8f8.3))") th0s
!
! double allnkm(n1_scha,n2_scha,n3_scha):
  istat = pio_inq_varid(ncid,'allnkm',id)
  istat = pio_get_var(ncid,id,allnkm)
! write(iulog,"('wei05sc: allnkm min,max=',2e12.4)") minval(allnkm),maxval(allnkm)
!
! int ab(csize):
  istat = pio_inq_varid(ncid,'ab',id)
  istat = pio_get_var(ncid,id,ab)
! write(iulog,"('wei05sc: ab=',/,(10i4))") ab
!
! int ls(csize):
  istat = pio_inq_varid(ncid,'ls',id)
  istat = pio_get_var(ncid,id,ls)
! write(iulog,"('wei05sc: ls=',/,(10i4))") ls
!
! int ms(csize):
  istat = pio_inq_varid(ncid,'ms',id)
  istat = pio_get_var(ncid,id,ms)
! write(iulog,"('wei05sc: ms=',/,(10i4))") ms
!
! double ex_epot(nex):
  istat = pio_inq_varid(ncid,'ex_epot',id)
  istat = pio_get_var(ncid,id,ex_epot)
! write(iulog,"('wei05sc: ex_epot=',/,(8f8.3))") ex_epot
!
! double ex_bpot(nex):
  istat = pio_inq_varid(ncid,'ex_bpot',id)
  istat = pio_get_var(ncid,id,ex_bpot)
! write(iulog,"('wei05sc: ex_bpot=',/,(8f8.3))") ex_bpot
!
! double epot_schfits(csize,n_schfits):
  istat = pio_inq_varid(ncid,'epot_schfits',id)
  istat = pio_get_var(ncid,id,epot_schfits)
! write(iulog,"('wei05sc: epot_schfits min,max=',2e12.4)") &
!   minval(epot_schfits),maxval(epot_schfits)
!
! double bpot_schfits(csize,n_schfits):
  istat = pio_inq_varid(ncid,'bpot_schfits',id)
  istat = pio_get_var(ncid,id,bpot_schfits)
! write(iulog,"('wei05sc: bpot_schfits min,max=',2e12.4)") &
!   minval(bpot_schfits),maxval(bpot_schfits)
!
! double epot_alschfits(csize,n_alschfits):
  istat = pio_inq_varid(ncid,'epot_alschfits',id)
  istat = pio_get_var(ncid,id,epot_alschfits)
! write(iulog,"('wei05sc: epot_alschfits min,max=',2e12.4)") &
!   minval(epot_alschfits),maxval(epot_alschfits)
!
! double bpot_alschfits(csize,n_alschfits):
  istat = pio_inq_varid(ncid,'bpot_alschfits',id)
  istat = pio_get_var(ncid,id,bpot_alschfits)
! write(iulog,"('wei05sc: bpot_alschfits min,max=',2e12.4)") &
!   minval(bpot_alschfits),maxval(bpot_alschfits)
!
! Close file:
  call cam_pio_closefile(ncid)
  if(masterproc) write(iulog,"('wei05sc: completed read of file ',a)") trim(file)
#endif
  end subroutine read_wei05_ncfile
#ifdef WACCMX_IONOS
!-----------------------------------------------------------------------
  subroutine setmodel(angle,bt,tilt,swvel,swden,model)
!
! Calculate the complete set of the models' SCHA coeficients,
!   given an aribitrary IMF angle (degrees from northward toward +Y),
!   given byimf, bzimf, solar wind velocity (km/sec), and density.
!
  implicit none
!
! Args:
  real(r8),intent(in) :: angle,bt,tilt,swvel,swden
  character(len=*),intent(in) :: model
!
! Local:
  integer :: i,j
  real(r8) :: pi,stilt,stilt2,sw,swp,swe,c0,rang,cosa,sina,cos2a,sin2a
  real(r8) :: a(n_schfits)
!
  if (trim(model) /= 'epot'.and.trim(model) /= 'bpot') then
    write(iulog,"('>>> model=',a)") trim(model)
    write(iulog,"('>>> setmodel: model must be either','''epot'' or ''bpot''')")
    call endrun('setmodel')
  endif
!
  pi = 4._r8*atan(1._r8)
  rad2deg = 180._r8/pi
  deg2rad = pi/180._r8
!
! write(iulog,"('setmodel call setboundary: model=',a,' swvel=',e12.4)") &
!   model, swvel

  call setboundary(angle,bt,swvel,swden)
!
  stilt = sin(tilt*deg2rad)
  stilt2 = stilt**2
  sw = bt*swvel/1000._r8
  if (trim(model) == 'epot') then
    swe = (1._r8-exp(-sw*ex_epot(2)))*sw**ex_epot(1)
  else
    swe = (1._r8-exp(-sw*ex_bpot(2)))*sw**ex_bpot(1)
  endif
  c0 = 1._r8
  swp = swvel**2 * swden*1.6726e-6_r8
  rang = angle*deg2rad
  cosa = cos(rang)
  sina = sin(rang)
  cos2a = cos(2._r8*rang)
  sin2a = sin(2._r8*rang)
  if (bt < 1._r8) then ! remove angle dependency for IMF under 1 nT
    cosa = -1._r8+bt*(cosa+1._r8)
    cos2a = 1._r8+bt*(cos2a-1._r8)
    sina = bt*sina
    sin2a = bt*sin2a
  endif
  a = (/c0      , swe       , stilt      , stilt2     , swp, & 
    swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa, &        
    swe*sina, stilt*sina, stilt2*sina, swp*sina, &        
    swe*cos2a,swe*sin2a/)
  if (trim(model) == 'epot') then
    esphc(:) = 0._r8
    do j=1,csize
      do i=1,n_schfits
        esphc(j) = esphc(j)+epot_schfits(i,j)*a(i)
      enddo
    enddo
!   write(iulog,"('setmodel: esphc=',/,(6e12.4))") esphc
  else
    bsphc(:) = 0._r8
    do j=1,csize
      do i=1,n_schfits
        bsphc(j) = bsphc(j)+bpot_schfits(i,j)*a(i)
      enddo
    enddo
!   write(iulog,"('setmodel: bsphc=',/,(6e12.4))") bsphc
  endif
  end subroutine setmodel

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine wei05loc (ih, byimf, power, sunlons)
! ih=1,2 for SH,NH called from weimer05
!
! (dimension 2 is for south, north hemispheres)
!  Calculate offa, dskofa, rrad, phid, and phin from Weimer 2005 offc, dskofc, theta0
!   Use Fig 8 of Heelis et al. [JGR, 85, 3315-3324, 1980]
!     This shows:  arad = 18.7 deg, crad = 16.7 deg (so arad = crad + 2 deg)
!           offa = offc = 3 deg (so offa = offc)
!           dskofc = 2 deg, dskofa = -0.5 deg  (so dskofa = dskofc - 2.5 deg)
!   Parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
!   (In aurora_cons, phid=0., phin=180.*rtd)
!     (For zero By, should be phid=21.39MLT*15*rtd, phin=11.5*15*rtd)
!  05/08:  But formulae for ra-rc using IMF CP between 5-7 deg (not 2) so use
!          difference of ra(max IMF CP or HP)-rc(IMF CP) as in aurora.F
! These are the dimensions and descriptions (corrected phid,n) from aurora.F:
!      theta0(2), ! convection reversal boundary in radians
!      offa(2),   ! offset of oval towards 0 MLT relative to magnetic pole (rad)
!      dskofa(2), ! offset of oval in radians towards 18 MLT (f(By))
!      phid(2),   ! dayside convection entrance in MLT-12 converted to radians (f(By))
!                      phid is the MLT-12 location of the cusp on the dayside
!      phin(2),   ! night convection entrance in MLT-12 converted to radians (f(By))
!      rrad(2),   ! radius of auroral circle in radians
!      offc(2),   ! offset of convection towards 0 MLT relative to mag pole (rad)
!      dskofc(2)  ! offset of convection in radians towards 18 MLT (f(By))
! sunlons(nlat): sun's longitude in dipole coordinates (see sub sunloc)
!

!
! Args:
      integer,intent(in) :: ih
      real(r8),intent(in) :: byimf
      real(r8),intent(in) :: power
      real(r8),intent(in) :: sunlons(:)
!
! Local:
      real(r8) :: rccp,racp,rahp,ramx,diffrac,plevel,tmltmin,tmltmax
      real(r8) :: offcdegp(2)
      integer :: i,j,j1,j2
      real(r8) :: vnx(2,2),hem,mltd,mltn
      integer :: inx(2,2)
      real(r8) :: offcdeg,dskof,arad,crad
      real(r8) :: byloc
        
! Limit size of byimf in phin and phid calculations (as in aurora.F) 
!  NOTE:  This byloc is assymetric in hemisphere, which is probably not correct
      byloc = byimf
      if (byloc .gt. 7._r8) byloc = 7._r8
      if (byloc .lt. -11._r8) byloc = -11._r8
!
!  ih=1 is SH, ih=2 is NH
	if (ih .eq. 1) then
	  j1 = 1
	  j2 = nmlat0
	  hem = -1._r8
	else
	  j1 = nmlat0 + 1
	  j2 = nmlat
	  hem = 1._r8
	endif
! Print out un-revised values:
!       write (6,"(1x,'Original convection/oval params (hem,By,off,dsk',
!    |    ',rad,phid,n=',10f9.4)") hem,byimf,offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.
!  Find min/max
	vnx(ih,1) = 0._r8
	vnx(ih,2) = 0._r8
	do j=j1,j2
	  do i=1,nmlonp1-1
	    if (phihm(i,j) .gt. vnx(ih,2)) then
	      vnx(ih,2) = phihm(i,j)
	      inx(ih,2) = i
	    endif
	    if (phihm(i,j) .lt. vnx(ih,1)) then
	      vnx(ih,1) = phihm(i,j)
	      inx(ih,1) = i
	    endif
	  enddo  !  i=1,nmlonp1-1
	enddo  !  j=j1,j2
! 05/08: Calculate weictpoten in kV from Weimer model min/max in V
	weictpoten(ih) = 0.001_r8 * (vnx(ih,2) - vnx(ih,1))
	tmltmin = (ylonm(inx(ih,1))-sunlons(1)) * r2d/15._r8 + 12._r8
	if (tmltmin .gt. 24._r8) tmltmin = tmltmin - 24._r8
	tmltmax = (ylonm(inx(ih,2))-sunlons(1)) * r2d/15._r8 + 12._r8
	if (tmltmax .gt. 24._r8) tmltmax = tmltmax - 24._r8
!       write (6,"('ih Bz By Hp ctpoten,wei min/max potV,lat,mlt=',i2,
!    |    5f8.2,2x,e12.4,2f8.2,2x,e12.4,2f8.2))") ih,bzimf,byimf,power,
!    |    ctpoten,weictpoten(ih),
!    |    vnx(ih,1),ylatm(jnx(ih,1))*rtd,tmltmin,
!    |    vnx(ih,2),ylatm(jnx(ih,2))*rtd,tmltmax
! 05/08: From aurora_cons, calculate convection and aurora radii using IMF convection
!   and power (plevel);  racp (DMSP/NOAA) - rccp (AMIE) = 5.32 (Bz>0) to 6.62 (Bz<0) deg
!  Heelis et al [1980, JGR, 85, pp 3315-3324] Fig 8: ra=rc+2deg, and is 2.5 deg to dusk
	rccp = -3.80_r8+8.48_r8*(weictpoten(ih)**0.1875_r8)
	racp = -0.43_r8+9.69_r8*(weictpoten(ih)**0.1875_r8)
	plevel = 0._r8
	if (power >=1.00_r8) plevel = 2.09_r8*log(power)
	rahp = 14.20_r8 + 0.96_r8*plevel
	ramx = max(racp,rahp)
	diffrac = ramx - rccp

!  Set default values
!  Use parameterization defaults for phid (phid(MLT)=9.39 +/- 0.21By - 12)
!                             and phin (phin(MLT)=23.50 +/- 0.15By - 12)
	mltd = 9.39_r8 - hem*0.21_r8*byloc
	mltn = 23.50_r8 - hem*0.15_r8*byloc
	phid(ih) = (mltd-12._r8) * 15._r8 *d2r
	phin(ih) = (mltn-12._r8) * 15._r8 *d2r
! 05/18/08:  Note that phid,phin are only for Heelis and are irrelevant for Weimer
!       write (6,"(1x,'mltd mltn phid,n =',4f8.2)")
!    |   mltd,mltn,phid(ih)*rtd/15.,phin(ih)*rtd/15.
!  Use default constant value of offcdegp from setboundary in Weimer 2005
	offcdeg = 4.2_r8
	offcdegp(ih) = offcdeg
	offc(ih) = offcdegp(ih) *d2r
	offa(ih) = offcdegp(ih) *d2r
!       write (6,"(1x,'offcdeg,rad =',2e12.4)") offcdeg,offc(ih)
	dskof = 0._r8
        dskofc(ih) = dskof *d2r
!  oval offset is 2.5 deg towards dawn (more neg dskof)
	dskofa(ih) = (dskof-2.5_r8) *d2r
!       write (6,"(1x,'dskof,c,a=',3f8.2)")
!    |    dskof,dskofc(ih)*rtd,dskofa(ih)*rtd
! Set crad from bndyfitr/2 of setboundary of Weimer 2005
	crad = bndyfitr/2._r8
!      write (6,"(1x,'wei05loc: ih,bz,y,crad =',i2,3f8.2)") 
!    |    ih,bzimf,byimf,crad
!  Fig 8 Heelis et al [1980]: ra=rc+2deg, and shifted 2.5 deg to dusk
	arad = crad + 2._r8
! 05/08:  Make ra=rc+diffrac(=ramx-rccp) - same difference as in aurora.F
! Choose to have arad=crad(Weimer) + diffrac(same diff as in aurora.F)
	arad = crad + diffrac
! 08/08: OR make ra=ramx=max(racp,rahp) so diffrac=arad-crad
!	diffrac2 = ramx - crad
! Choose to have arad=ramx (same as in aurora.F as determined by P/CP)
!       arad = ramx
	theta0(ih) = crad *d2r
	rrad(ih) = arad *d2r
!       write (6,"(1x,'radius: crad,rccp,racp,rahp diffa-c',
!    |   '(aurF,ramx-Weic) ramx,Weic+d,arad deg=',9f8.2)") crad,rccp,
!    |   racp,rahp,diffrac,diffrac2,ramx,crad+diffrac,arad

! Print out revised values (revised 05/08):
!       write (6,"(1x,'Revised convection/oval params (off,dsk,',
!    |    'rad,phid,n=',8f9.4)")offc(ih)*rtd,offa(ih)*rtd,
!    |    dskofc(ih)*rtd,dskofa(ih)*rtd,theta0(ih)*rtd,rrad(ih)*rtd,
!    |    phid(ih)*rtd/15.+12.,phin(ih)*rtd/15.+12.

      end subroutine wei05loc

!-----------------------------------------------------------------------
! for now this is here ... might need to move to a gen util module      
!-----------------------------------------------------------------------
    function hp_from_bz_swvel(bz,swvel) result(hp)
!
! Calculate hemispheric power from bz, swvel:
! Emery, et.al., (2008, in press, JGR)
! 6/3/08: Enforce minimum hp of 4.0 before *fac.
! 6/6/08: Reset minimum hp from 4.0 to 2.5 before *fac,
!         as per Emery email of 6/5/08.
!
      real(r8),intent(in) :: bz,swvel  ! in
      real(r8) :: hp                   ! out

      real(r8), parameter :: fac = 2.0_r8
!
      if (bz < 0._r8) then
        hp = 6.0_r8 + 3.3_r8*abs(bz) + (0.05_r8 + 0.003_r8*abs(bz))* (min(swvel,700._r8)-300._r8)
      else
        hp = 5.0_r8 + 0.05_r8 * (min(swvel,700._r8)-300._r8)
      endif
      hp = max(2.5_r8,hp)*fac

    end function hp_from_bz_swvel

!-----------------------------------------------------------------------
  subroutine setboundary(angle,bt,swvel,swden)
!
! Sets the coefficients that define the low-latitude boundary model,
!   given the IMF and solar wind values.
!
  implicit none
!
! Args:
  real(r8),intent(in) :: angle,bt,swvel,swden
!
! Local:
  integer :: i
  real(r8) :: swp,xc,theta,ct,st,cosa,btx,x(na),c(na)
!
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
  xc = 4.2_r8
  theta = xc*(deg2rad)
  ct = cos(theta)
  st = sin(theta)
!
  tmat(1,:) = (/ ct, 0._r8, st/) 
  tmat(2,:) = (/ 0._r8, 1._r8, 0._r8/) 
  tmat(3,:) = (/-st, 0._r8, ct/)
!
!  ttmat(1,:) = (/ct, 0._r8,-st/)
!  ttmat(2,:) = (/ 0._r8,1._r8, 0._r8/)
!  ttmat(3,:) = (/st, 0._r8, ct/)
!
  swp = swden*swvel**2*1.6726e-6_r8 ! pressure
  cosa = cos(angle*deg2rad)
  btx = 1._r8-exp(-bt*ex_bndy(1))
  if (bt > 1._r8) then
    btx = btx*bt**ex_bndy(2)
  else
    cosa = 1._r8+bt*(cosa-1._r8) ! remove angle dependency for IMF under 1 nT
  endif
  x = (/1._r8, cosa, btx, btx*cosa, swvel, swp/)
  c = bndya
  bndyfitr = 0._r8
  do i=1,na
    bndyfitr = bndyfitr+x(i)*c(i)

!   write(iulog,"('setboundry: i=',i3,' bndyfitr=',e12.4)") i,bndyfitr

  enddo
  end subroutine setboundary
!-----------------------------------------------------------------------
  subroutine epotval(lat,mlt,fill,epot)
!
! Return the Potential (in kV) at given combination of def. latitude 
!   (lat) and MLT, in geomagnetic apex coordinates (practically identical 
!   to AACGM).  
! If the location is outside of the model's low-latitude boundary, then 
!   the value "fill" is returned.
!
  implicit none
!
! Args:
  real(r8),intent(in) :: lat,mlt,fill
  real(r8),intent(out) :: epot
!
! Local:
  integer :: inside,j,m,skip
  real(r8) :: z,phir,plm,colat,nlm
  real(r8) :: phim(2),cospm(2),sinpm(2)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
  call checkinputs(lat,mlt,inside,phir,colat)
  if (inside == 0) then
    epot = fill
    return
  endif
!
! IDL code: 
! phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
! Here, maxm == maxm_pot == 2, and phir is a scalar. The above IDL statement then 
!   becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where phim will be 
!   dimensioned [1,2]
!
  phim(1) = phir
  phim(2) = phir*2._r8
  cospm(:) = cos(phim(:))
  sinpm(:) = sin(phim(:))
!
  z = 0._r8
  skip = 0 ! Added by B.Foster, 4/23/14
  do j=1,csize
    if (skip == 1) then
      skip = 0
      cycle
    endif
    m = ms(j)
    if (ab(j)==1) then
      plm = scplm(j,colat,nlm) ! scplm function is in this module
      skip = 0
      if (m == 0) then
        z = z+plm*esphc(j)
      else
        z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
        skip = 1
      endif
    endif ! ab(j)
  enddo
  epot = z 
  end subroutine epotval
!-----------------------------------------------------------------------
  subroutine mpfac(lat,mlt,fill,fac)
  implicit none
!
! Args:
  real(r8),intent(in) :: lat,mlt,fill
  real(r8),intent(out) :: fac
!
! Local:
  integer :: j,m,inside,skip
  real(r8) :: phim(2),cospm(2),sinpm(2),cfactor
  real(r8) :: re,z,phir,plm,colat,nlm,pi
!
  re = 6371.2_r8 + 110._r8 ! km radius (allow default ht=110)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
  call checkinputs(lat,mlt,inside,phir,colat)
  if (inside == 0) then
    fac = fill
    return
  endif
!
  phim(1) = phir
  phim(2) = phir*2._r8
  cospm(:) = cos(phim(:))
  sinpm(:) = sin(phim(:))
!
  z = 0._r8
  skip = 0 ! Added by B.Foster, 4/23/14
  jloop: do j=1,csize
    if (skip == 1) then
      skip = 0
      cycle
    endif
    if (ls(j) >= 11) exit jloop
    m = ms(j)
    if (ab(j) == 1) then
      plm = scplm(j,colat,nlm) ! colat and nlm are returned (both reals)
      plm = plm*(nlm*(nlm+1._r8))
!
! bsphc was calculated in setmodel (when setmodel called with 'bpot')
      if (m==0) then
        z = z-plm*bsphc(j)
      else
        z = z-(plm*(bsphc(j)*cospm(m)+bsphc(j+1)*sinpm(m)))
        skip = 1
      endif
    endif
  enddo jloop ! j=1,csize
  pi = 4._r8*atan(1._r8)
  cfactor = -1.e5_r8/(4._r8*pi*re**2) ! convert to uA/m2
  z = z*cfactor
  fac = z
! write(iulog,"('mpfac: lat=',f8.3,' mlt=',f8.3,' fac=',1pe12.4)") lat,mlt,fac
  end subroutine mpfac
!-----------------------------------------------------------------------
  real(r8) function scplm(index,colat,nlm)
!
! Return Spherical Cap Harmonic Associated Legendre values, given colat 
!   values and index i into array of L and M values.
!
  implicit none
!
! Args:
  integer,intent(in) :: index
  real(r8),intent(in) :: colat
  real(r8),intent(out) :: nlm
!
! Local:
  integer :: i,j,l,m,skip
  real(r8) :: th0,out(1),colata(1)
  real(r8) :: cth(mxtablesize)
  real(r8),save :: prevth0=1.e36_r8
  integer,save :: tablesize
!
  scplm = 0._r8
  skip = 0 ! Added by B.Foster, 4/23/14
  th0 = bndyfitr
  if (prevth0 /= th0) then
    tablesize = 3*nint(th0)
    if (tablesize > mxtablesize) then 
      write(iulog,"('>>> tablesize > mxtablesize: tablesize=',i8,' mxtablesize=',i8,' th0=',e12.4)") &
        tablesize,mxtablesize,th0
      call endrun('tablesize')
    endif
    do i=1,tablesize
      real8 = dble(i-1)
      real8a = dble(tablesize-1)
      colattable(i) = real8*(th0/real8a)
      cth(i) = cos(colattable(i)*deg2rad)
    enddo
    prevth0 = th0
    nlms = 0._r8 ! whole array init 
    do j=1,csize
      if (skip == 1) then
        skip = 0
        cycle
      endif
      l = ls(j)
      m = ms(j)
      nlms(j) = nkmlookup(l,m,th0) ! nkmlookup in this module

! real :: plmtable(mxtablesize,csize)
      call pm_n(m,nlms(j),cth,plmtable(1:tablesize,j),tablesize)
      skip = 0
      if (m /= 0 .and. ab(j) > 0) then
        plmtable(1,j+1) = plmtable(1,j)
        nlms(j+1) = nlms(j)
        skip = 1
      endif
    enddo ! j=1,csize
  endif ! prevth0
  nlm = nlms(index)
  colata(1) = colat
  call interpol_quad(plmtable(1:tablesize,index), &
    colattable(1:tablesize),colata,out)
  scplm = out(1)
  end function scplm
!-----------------------------------------------------------------------
  subroutine pm_n(m,r,cth,plmtable,tablesize)
!
! Another SCHA function, returns the SCHA version of the associated 
! Legendre Polynomial, Pmn
!
  implicit none
!
! Args:
  integer,intent(in) :: m,tablesize
  real(r8),intent(in) :: r
  real(r8),intent(in) :: cth(tablesize)
  real(r8),intent(out) :: plmtable(tablesize)
!
! Local:
  integer :: i,k
  real(r8) :: rm,rk,div,ans,xn
  real(r8),dimension(tablesize) :: a,x,tmp,table
!
  if (m == 0) then 
    a = 1._r8 ! whole array op
  else
    do i=1,tablesize
      a(i) = sqrt(1._r8-cth(i)**2)**m
    enddo
  endif
  xn = r*(r+1._r8)
  x(:) = (1._r8-cth(:))/2._r8
  table = a ! whole array init
  k = 1
  pmn_loop: do         ! repeat-until loop in idl code
    do i=1,tablesize
      real8 = dble(m)
      rm = real8
      real8 = dble(k)
      rk = real8
      a(i) = a(i)*(x(i)*((rk+rm-1._r8)*(rk+rm)-xn)/(rk*(rk+rm)))
      table(i) = table(i)+a(i) ! "result" in idl code
    enddo
    k = k+1
    do i=1,tablesize
      div = abs(table(i))
      if (div <= 1.e-6_r8) div = 1.e-6_r8
      tmp(i) = abs(a(i)) / div
    enddo
    if (maxval(tmp) < 1.e-6_r8) exit pmn_loop
  enddo pmn_loop
  ans = km_n(m,r)

  plmtable(:) = table(:)*ans
  end subroutine pm_n
!-----------------------------------------------------------------------
  real(r8) function km_n(m,rn)
!
! A normalization function used by the SCHA routines.  See Haines.
!
  implicit none
!
! Args:
  integer,intent(in) :: m
  real(r8),intent(in) :: rn
!
! Local:
  real(r8) :: rm
!
  if (m == 0) then 
    km_n = 1._r8
    return
  endif
  real8 = dble(m)
  rm = real8
  km_n = sqrt(2._r8*exp(lngamma(rn+rm+1._r8)-lngamma(rn-rm+1._r8))) / &
    (2._r8**m*factorial(m))
  end function km_n
!-----------------------------------------------------------------------
  real(r8) function nkmlookup(k,m,th0)
!
! Given the size of a spherical cap, defined by the polar cap angle, th0, 
!   and also the values of integers k and m, returns the value of n, a 
!   real number (see Haines).
! It uses interpolation from a lookup table that had been precomputed, 
!   in order to reduce the computation time.
!
  implicit none
!
! Args:
  integer,intent(in) :: k,m
  real(r8),intent(in) :: th0
!
! Local:
  integer :: kk,mm
  real(r8) :: th0a(1),out(1)

  if (th0 == 90._r8) then
    real8 = dble(k)
    nkmlookup = real8
    return
  endif
  th0a(1) = th0
  kk = k+1
  mm = m+1
  if (kk > maxk_scha) then
    call interpol_quad(allnkm(maxk_scha,mm,:),th0s,th0a,out)
  endif
  if (mm > maxm_scha) then
    call interpol_quad(allnkm(kk,maxm_scha,:),th0s,th0a,out)
  endif
  if (th0 < th0s(1)) then
    write(iulog,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,' th0s(1)=',e12.4)") &
      th0,th0s(1)
  endif
  call interpol_quad(allnkm(kk,mm,:),th0s,th0a,out)
  nkmlookup = out(1)
  end function nkmlookup
!-----------------------------------------------------------------------
  subroutine checkinputs(lat,mlt,inside,phir,colat)
  implicit none
!
! Args:
  real(r8),intent(in) :: lat,mlt
  integer,intent(out) :: inside
  real(r8),intent(out) :: phir,colat
!
! Local:
  real(r8) :: lon,tlat,tlon,radii
!
  lon = mlt*15._r8
  call dorotation(lat,lon,tlat,tlon)
  radii = 90._r8-tlat
  inside = 0
  if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
  phir = tlon*deg2rad
  colat = radii
  end subroutine checkinputs
!-----------------------------------------------------------------------
  subroutine dorotation(latin,lonin,latout,lonout)
!
! Uses transformation matrices tmat and ttmat, to convert between
!   the given geomagnetic latatud/longitude, and the coordinate 
!   system that is used within the model,that is offset from the pole.
!
! Rotate Lat/Lon spherical coordinates with the transformation given 
!   by saved matrix. The coordinates are assumed to be on a sphere of 
!   Radius=1. Uses cartesian coordinates as an intermediate step.
!
  implicit none
!
! Args:
  real(r8),intent(in) :: latin,lonin
  real(r8),intent(out) :: latout,lonout
!
! Local:
  real(r8) :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
  integer :: i
!
  latr = latin*deg2rad
  lonr = lonin*deg2rad
  stc = sin(latr)
  ctc = cos(latr)
  sf = sin(lonr)
  cf = cos(lonr)
  a = ctc*cf
  b = ctc*sf
!
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
!   Pos is then (1,3)
!
  do i=1,3
    pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
  enddo
  latout = asin(pos(3))*rad2deg
  lonout = atan2(pos(2),pos(1))*rad2deg
  end subroutine dorotation
!-----------------------------------------------------------------------
  subroutine interpol_quad(v,x,u,p)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
  implicit none
!
! Args:
  real(r8),intent(in) :: v(:),x(:),u(:)
  real(r8),intent(out) :: p(:)
!
! Local:
  integer :: nv,nx,nu,i,ix
  real(r8) :: x0,x1,x2
!
  nv = size(v)
  nx = size(x)
  nu = size(u)
  if (nx /= nv) then
    p(:) = 0._r8
    return
  endif
  do i=1,nu
    ix = value_locate(x,u(i))
! 01/14 bae: interpol_quad in wei05sc.F is called when inside=1 or radii<bndryfit 
!  for Weimer 2005. The fix to ix<=1 and ix>=nx assures epot is non-zero near 
!  the pole (85.8mlat,0MLT) and the boundary (bndryfit).
    if (ix <=1) ix = 2        ! bug fix by bae 01/28/14
    if (ix >=nx) ix = nx-1            ! bug fix by bae 01/29/14
!   if (ix <= 1.or.ix >= nx) then ! bug fix by btf 12/23/09
!      p(i) = 0._r8
!      cycle                       ! bug fix by btf 12/23/09
!    endif
    x1 = x(ix)
    x0 = x(ix-1)
    x2 = x(ix+1)
    p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) + & 
           v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) + & 
           v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
  enddo
  end subroutine interpol_quad
!-----------------------------------------------------------------------
  integer function value_locate(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val >= vec(i+1)
! Input vec must be monotonically increasing
!
  implicit none
!
! Args:
  real(r8),intent(in) :: vec(:),val
!
! Local:
  integer :: n,i
!
  value_locate = 0
  n = size(vec)
  if (val < vec(1)) return
  if (val > vec(n)) then
    value_locate = n
    return
  endif
  do i=1,n-1
    if (val >= vec(i) .and. val <= vec(i+1)) then
      value_locate = i
      return
    endif
  enddo
  end function value_locate
!-----------------------------------------------------------------------
  real(r8) function lngamma(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
  implicit none
  real(r8),intent(in) :: xx
  real(r8) :: x,y,tmp,ser
  real(r8) :: cof(6) = (/76.18009172947146_r8, -86.50532032941677_r8,   &
    24.01409824083091_r8, -1.231739572450155_r8, 0.1208650973866179e-2_r8, &
    -0.5395239384953e-5_r8/)
  integer :: j
!
  y = xx
  x = xx
  tmp = x+5.5_r8
  tmp = tmp-(x+0.5_r8)*log(tmp)
  ser = 1.000000000190015_r8
  do j=1,5
    y = y+1
    ser = ser+cof(j)/y
  enddo
  lngamma = -tmp+log(2.5066282746310005_r8*ser/x)
  end function lngamma
!-----------------------------------------------------------------------
  real(r8) function factorial(n)
  implicit none
  integer,intent(in) :: n
  integer :: m
  if (n <= 0) then
    factorial = 0._r8
    return
  endif
  if (n == 1) then
    factorial = 1._r8
    return
  endif
  real8 = dble(n)
  factorial = real8
  do m = n-1,1,-1
    real8 = dble(m)
    factorial = factorial * real8
  enddo
  end function factorial
!-----------------------------------------------------------------------
!*********************** Copyright 1996,2001 Dan Weimer/MRC ***********************
! COORDINATE TRANSFORMATION UTILITIES

!NCAR      Feb 01:  Changed TRANS to GET_TILT s.t. the dipole tilt angle is
!          returned.

        real(r8) FUNCTION GET_TILT (YEAR,MONTH,DAY,HOUR)
!       SUBROUTINE TRANS(YEAR,MONTH,DAY,HOUR,IDBUG)
        implicit none
        real(r8) :: B3,B32,B33
        integer :: IYR,JD,MJD,I,J,K
!NCAR

        INTEGER YEAR,MONTH,DAY,IDBUG
        real(r8) :: HOUR
!         
!      THIS SUBROUTINE DERIVES THE ROTATION MATRICES AM(I,J,K) FOR 11
!      TRANSFORMATIONS, IDENTIFIED BY K.
!          K=1 TRANSFORMS GSE to GEO
!          K=2     "      GEO to MAG
!          K=3     "      GSE to MAG
!          K=4     "      GSE to GSM
!          K=5     "      GEO to GSM
!          K=6     "      GSM to MAG
!          K=7     "      GSE to GEI
!          K=8     "      GEI to GEO
!          K=9     "      GSM to SM 
!          K=10    "      GEO to SM 
!          K=11    "      MAG to SM 
!
!      IF IDBUG IS NOT 0, THEN OUTPUTS DIAGNOSTIC INFORMATION TO
!      FILE UNIT=IDBUG
!       
        INTEGER GSEGEO,GEOGSE,GEOMAG,MAGGEO
        INTEGER GSEMAG,MAGGSE,GSEGSM,GSMGSE
        INTEGER GEOGSM,GSMGEO,GSMMAG,MAGGSM
        INTEGER GSEGEI,GEIGSE,GEIGEO,GEOGEI
        INTEGER GSMSM,SMGSM,GEOSM,SMGEO,MAGSM,SMMAG

        PARAMETER (GSEGEO= 1,GEOGSE=-1,GEOMAG= 2,MAGGEO=-2)
        PARAMETER (GSEMAG= 3,MAGGSE=-3,GSEGSM= 4,GSMGSE=-4)
        PARAMETER (GEOGSM= 5,GSMGEO=-5,GSMMAG= 6,MAGGSM=-6)
        PARAMETER (GSEGEI= 7,GEIGSE=-7,GEIGEO= 8,GEOGEI=-8)
        PARAMETER (GSMSM = 9,SMGSM =-9,GEOSM =10,SMGEO=-10)
        PARAMETER (MAGSM =11,SMMAG =-11)
!
!      The formal names of the coordinate systems are:
!       GSE - Geocentric Solar Ecliptic
!       GEO - Geographic
!       MAG - Geomagnetic
!       GSM - Geocentric Solar Magnetospheric
!       SM  - Solar Magnetic
!       
!      THE ARRAY CX(I) ENCODES VARIOUS ANGLES, STORED IN DEGREES
!      ST(I) AND CT(I) ARE SINES & COSINES.       
!
!      Program author:  D. R. Weimer
!
!      Some of this code has been copied from subroutines which had been
!      obtained from D. Stern, NASA/GSFC.  Other formulas are from "Space 
!      Physics Coordinate Transformations: A User Guide" by M. Hapgood (1991).
!
!      The formulas for the calculation of Greenwich mean sidereal time (GMST)
!      and the sun's location are from "Almanac for Computers 1990",
!      U.S. Naval Observatory.
!
        real(r8) :: UT,T0,GMSTD,GMSTH,ECLIP,MA,LAMD,SUNLON

!NCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
!          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
!          are ADJUST and JULDAY
!NCAR      Nov 02: Commons MFIELD and TRANSDAT now only in TRANS (GET_TILT)
!                  So eliminate them as commons.  For Fortran 90, eliminate
!                  the DATA statement for assignments (not block_data)
!       COMMON/MFIELD/EPOCH,TH0,PH0,DIPOLE
!       COMMON/TRANSDAT/CX(9),ST(6),CT(6),AM(3,3,11)
!         
        real(r8) TH0,PH0 !,DIPOLE
        real(r8) CX(9),ST(6),CT(6),AM(3,3,11)
!         
!  TH0 = geog co-lat of NH magnetic pole
!  PH0 = geog longitude of NH magnetic pole
!  DIPOLE = magnitude of the B field in gauss at the equator

        TH0 = 11.19_r8
        PH0 = -70.76_r8
!        DIPOLE = 0.30574_r8
!NCAR

!NCAR      Feb 01:  Prevent debug printing to a file
        IDBUG = 0
!NCAR

        IF(YEAR.LT.1900)THEN
          IYR=1900+YEAR
        ELSE
          IYR=YEAR
        ENDIF
        UT=HOUR
        JD=JULDAY(MONTH,DAY,IYR)
        MJD=JD-2400001
        real8 = dble(MJD)
        T0=(real8-51544.5_r8)/36525.0_r8
        GMSTD=100.4606184_r8 + 36000.770_r8*T0 + 3.87933E-4_r8*T0*T0 + &
              15.0410686_r8*UT
        CALL ADJUST(GMSTD)
        GMSTH=GMSTD*24._r8/360._r8
        ECLIP=23.439_r8 - 0.013_r8*T0
        MA=357.528_r8 + 35999.050_r8*T0 + 0.041066678_r8*UT
        CALL ADJUST(MA)
        LAMD=280.460_r8 + 36000.772_r8*T0 + 0.041068642_r8*UT
        CALL ADJUST(LAMD)
        SUNLON=LAMD + (1.915_r8-0.0048_r8*T0)*SIND(MA) + 0.020_r8*SIND(2._r8*MA)
        CALL ADJUST(SUNLON)
        IF(IDBUG.NE.0)THEN
          WRITE(IDBUG,*) YEAR,MONTH,DAY,HOUR
          WRITE(IDBUG,*) 'MJD=',MJD
          WRITE(IDBUG,*) 'T0=',T0
          WRITE(IDBUG,*) 'GMSTH=',GMSTH
          WRITE(IDBUG,*) 'ECLIPTIC OBLIQUITY=',ECLIP
          WRITE(IDBUG,*) 'MEAN ANOMALY=',MA
          WRITE(IDBUG,*) 'MEAN LONGITUDE=',LAMD
          WRITE(IDBUG,*) 'TRUE LONGITUDE=',SUNLON
        ENDIF

        CX(1)= GMSTD
        CX(2) = ECLIP
        CX(3) = SUNLON
        CX(4) = TH0
        CX(5) = PH0
! Derived later:
!       CX(6) = Dipole tilt angle  
!       CX(7) = Angle between sun and magnetic pole
!       CX(8) = Subsolar point latitude
!       CX(9) = Subsolar point longitude

        DO I=1,5
          ST(I) = SIND(CX(I))
          CT(I) = COSD(CX(I))
        ENDDO
!         
      AM(1,1,GSEGEI) = CT(3)
      AM(1,2,GSEGEI) = -ST(3)
      AM(1,3,GSEGEI) = 0._r8
      AM(2,1,GSEGEI) = ST(3)*CT(2)
      AM(2,2,GSEGEI) = CT(3)*CT(2)
      AM(2,3,GSEGEI) = -ST(2)
      AM(3,1,GSEGEI) = ST(3)*ST(2)
      AM(3,2,GSEGEI) = CT(3)*ST(2)
      AM(3,3,GSEGEI) = CT(2)      
!         
      AM(1,1,GEIGEO) = CT(1)      
      AM(1,2,GEIGEO) = ST(1)      
      AM(1,3,GEIGEO) = 0._r8         
      AM(2,1,GEIGEO) = -ST(1)     
      AM(2,2,GEIGEO) = CT(1)      
      AM(2,3,GEIGEO) = 0._r8         
      AM(3,1,GEIGEO) = 0._r8         
      AM(3,2,GEIGEO) = 0._r8         
      AM(3,3,GEIGEO) = 1._r8         
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSEGEO) = AM(I,1,GEIGEO)*AM(1,J,GSEGEI) + &
          AM(I,2,GEIGEO)*AM(2,J,GSEGEI) +  AM(I,3,GEIGEO)*AM(3,J,GSEGEI)
      ENDDO
      ENDDO
!         
      AM(1,1,GEOMAG) = CT(4)*CT(5) 
      AM(1,2,GEOMAG) = CT(4)*ST(5) 
      AM(1,3,GEOMAG) =-ST(4)       
      AM(2,1,GEOMAG) =-ST(5)       
      AM(2,2,GEOMAG) = CT(5)       
      AM(2,3,GEOMAG) = 0._r8
      AM(3,1,GEOMAG) = ST(4)*CT(5) 
      AM(3,2,GEOMAG) = ST(4)*ST(5) 
      AM(3,3,GEOMAG) = CT(4)       
!         
      DO I=1,3   
      DO J=1,3   
       AM(I,J,GSEMAG) = AM(I,1,GEOMAG)*AM(1,J,GSEGEO) + &
         AM(I,2,GEOMAG)*AM(2,J,GSEGEO) +  AM(I,3,GEOMAG)*AM(3,J,GSEGEO)
      ENDDO
      ENDDO
!         
      B32 = AM(3,2,GSEMAG)         
      B33 = AM(3,3,GSEMAG)         
      B3  = SQRT(B32*B32+B33*B33)       
      IF (B33.LE.0._r8) B3 = -B3    
!         
      AM(2,2,GSEGSM) = B33/B3      
      AM(3,3,GSEGSM) = AM(2,2,GSEGSM)   
      AM(3,2,GSEGSM) = B32/B3      
      AM(2,3,GSEGSM) =-AM(3,2,GSEGSM)   
      AM(1,1,GSEGSM) = 1._r8
      AM(1,2,GSEGSM) = 0._r8
      AM(1,3,GSEGSM) = 0._r8
      AM(2,1,GSEGSM) = 0._r8
      AM(3,1,GSEGSM) = 0._r8
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOGSM) = AM(I,1,GSEGSM)*AM(J,1,GSEGEO) + &
          AM(I,2,GSEGSM)*AM(J,2,GSEGEO) + AM(I,3,GSEGSM)*AM(J,3,GSEGEO)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GSMMAG) = AM(I,1,GEOMAG)*AM(J,1,GEOGSM) + &
         AM(I,2,GEOMAG)*AM(J,2,GEOGSM) + AM(I,3,GEOMAG)*AM(J,3,GEOGSM)
      ENDDO
      ENDDO
!
        ST(6) = AM(3,1,GSEMAG)       
        CT(6) = SQRT(1._r8-ST(6)*ST(6))      
        CX(6) = ASIND(ST(6))     

        AM(1,1,GSMSM) = CT(6)
        AM(1,2,GSMSM) = 0._r8
        AM(1,3,GSMSM) = -ST(6)
        AM(2,1,GSMSM) = 0._r8
        AM(2,2,GSMSM) = 1._r8
        AM(2,3,GSMSM) = 0._r8
        AM(3,1,GSMSM) = ST(6)
        AM(3,2,GSMSM) = 0._r8
        AM(3,3,GSMSM) = CT(6)
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,GEOSM) = AM(I,1,GSMSM)*AM(1,J,GEOGSM) + &
          AM(I,2,GSMSM)*AM(2,J,GEOGSM) +  AM(I,3,GSMSM)*AM(3,J,GEOGSM)
      ENDDO
      ENDDO
!         
      DO I=1,3   
      DO J=1,3   
        AM(I,J,MAGSM) = AM(I,1,GSMSM)*AM(J,1,GSMMAG) + &
         AM(I,2,GSMSM)*AM(J,2,GSMMAG) + AM(I,3,GSMSM)*AM(J,3,GSMMAG)
      ENDDO
      ENDDO
!
      CX(7)=ATAN2D( AM(2,1,11) , AM(1,1,11) )
      CX(8)=ASIND( AM(3,1,1) )
      CX(9)=ATAN2D( AM(2,1,1) , AM(1,1,1) )

      IF(IDBUG.NE.0)THEN
          WRITE(IDBUG,*) 'Dipole tilt angle=',CX(6)
          WRITE(IDBUG,*) 'Angle between sun and magnetic pole=',CX(7)
          WRITE(IDBUG,*) 'Subsolar point latitude=',CX(8)
          WRITE(IDBUG,*) 'Subsolar point longitude=',CX(9)

        DO K=1,11
         WRITE(IDBUG,1001) K
         DO I=1,3
           WRITE(IDBUG,1002) (AM(I,J,K),J=1,3)
         ENDDO
        ENDDO
 1001   FORMAT(' ROTATION MATRIX ',I2)
 1002   FORMAT(3F9.5)
      ENDIF

!NCAR      Mar 96: return the dipole tilt from this function call.
      GET_TILT = CX(6)
!NCAR

      RETURN
      end function get_tilt
!-----------------------------------------------------------------------
!NCAR      Feb 01:  Eliminate unused routines from translib.for: ROTATE,
!          ROTATEV, FROMCART, TOCART, MLT, MAGLONG, SUNLOC.  Remaining
!          are ADJUST and JULDAY
!NCAR
        SUBROUTINE ADJUST(ANGLE)
        implicit none
        real(r8) :: angle
!       ADJUST AN ANGLE IN DEGREES TO BE IN RANGE OF 0 TO 360.
 10     CONTINUE
        IF(ANGLE.LT.0._r8)THEN
          ANGLE=ANGLE+360._r8
          GOTO 10
        ENDIF
 20     CONTINUE
        IF(ANGLE.GE.360._r8)THEN
          ANGLE=ANGLE-360._r8
          GOTO 20
        ENDIF
        end subroutine adjust
!-----------------------------------------------------------------------
      integer FUNCTION JULDAY(MM,ID,IYYY)
      implicit none
      integer :: igreg, iyyy, mm, id, jy, jm, ja
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY.EQ.0) call endrun('There is no Year Zero.')
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25_r8*JY)+INT(30.6001_r8*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY).GE.IGREG) THEN
        JA=INT(0.01_r8*JY)
        JULDAY=JULDAY+2-JA+INT(0.25_r8*JA)
      ENDIF
      end function julday
!-----------------------------------------------------------------------
      SUBROUTINE CVT2MD(iulog,IYEAR,NDA,MON,DAY)
!          This sub converts NDA, the day number of the year, IYEAR,
!          into the appropriate month and day of month (integers)
      implicit none
      integer :: iulog,iyear,nda,mon,miss,numd,i
      INTEGER DAY
      INTEGER LMON(12)
      PARAMETER (MISS=-32767)
      SAVE LMON
      DATA LMON/31,28,31,30,31,30,31,31,30,31,30,31/
 
      LMON(2)=28
      IF(MOD(IYEAR,4) .EQ. 0)LMON(2)=29
 
      NUMD=0
      DO 100 I=1,12
      IF(NDA.GT.NUMD .AND. NDA.LE.NUMD+LMON(I))GO TO 200
      NUMD=NUMD+LMON(I)
  100 CONTINUE
      WRITE(iulog,'('' CVT2MD:  Unable to convert year & day of year'', &
                    I5,'','',I5,''to month & day of month'')')IYEAR,NDA
      MON = MISS
      DAY = MISS
      RETURN
  200 MON=I
      DAY=NDA-NUMD
      end subroutine cvt2md
!-----------------------------------------------------------------------
!
!NCAR      Routines added to work around non-ANSI trig functions which
!          input degrees instead of radians:  SIND, COSD, ASIND, ATAN2D

      FUNCTION SIND (DEG)
      implicit none
      real(r8) :: sind,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                  R2D = 57.2957795130823208767981548147_r8)
      SIND = SIN (DEG * D2R)
      end function sind
!-----------------------------------------------------------------------
      FUNCTION COSD (DEG)
      implicit none
      real(r8) :: cosd,d2r,r2d,deg
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                  R2D = 57.2957795130823208767981548147_r8)

      COSD = COS (DEG * D2R)
      end function cosd
!-----------------------------------------------------------------------
      FUNCTION ASIND (RNUM)
      implicit none
      real(r8) :: asind,d2r,r2d,rnum
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                  R2D = 57.2957795130823208767981548147_r8)
      ASIND = R2D * ASIN (RNUM)
      end function asind
!-----------------------------------------------------------------------
      FUNCTION ATAN2D (RNUM1,RNUM2)
      implicit none
      real(r8) :: atan2d,d2r,r2d,rnum1,rnum2
      PARAMETER ( D2R =  0.0174532925199432957692369076847_r8 , &
                  R2D = 57.2957795130823208767981548147_r8)
      ATAN2D = R2D * ATAN2 (RNUM1,RNUM2)
      end function atan2d
#endif
!-----------------------------------------------------------------------
end module wei05sc
