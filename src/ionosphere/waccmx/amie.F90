module amie_module
  !
  ! Module used to read data from the AMIE outputs (POT,mean energy,
  !   and energy flux).
  !

  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_logfile    ,only: iulog
  use spmd_utils     ,only: masterproc
  use cam_abortutils ,only: endrun
  use edyn_maggrid,   only: nmlat,nmlonp1
  use edyn_mpi       ,only: mlon0,mlon1,mlat0,mlat1, &
                            lon0,lon1,lat0,lat1
#ifdef WACCMX_EDYN_ESMF
  use edyn_params    ,only: finit
  use edyn_maggrid   ,only:  &
                            ylonm,  &     ! magnetic latitudes (nmlat) (radians)
                            ylatm         ! magnetic longtitudes (nmlonp1) (radians)
  use edyn_esmf      ,only: mag_efx,mag_kev,geo_efx,geo_kev
  use esmf           ,only: ESMF_FIELD        ! ESMF library module
  use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile
  use pio, only: pio_inq_dimid, pio_inquire_dimension, pio_inquire, pio_inq_varid
  use pio, only: file_desc_t, pio_noerr, pio_nowrite, pio_get_var
#endif
  implicit none

  private
  public :: init_amie, getamie
#ifdef WACCMX_EDYN_ESMF

  ! Define parameters for AMIE input data file:
  integer, parameter ::  &
       mxgdays = 10,     & ! maximum number of days of AMIE data
       mxtimes = 5881,   & ! maximum number of times of AMIE data per day
       ithtrns = 30,     & ! corresponding to trans lat 40-deg
       ithmx = 55,       & ! maximum number of latitudes of AMIE data
       jmxm = 2*ithmx-1, & ! maximum number of global latitudes
       lonmx = 36          ! maximum number of longitudes of AMIE data
  integer :: lonp1,latp1
  !     integer,dimension(mxtimes) :: year,month,day,jday
  ! Define AMIE output fields
  real(r8) :: &
       tiepot(nmlonp1,nmlat),tieekv(nmlonp1,nmlat),  &
       tieefx(nmlonp1,nmlat)
  ! defined output AMIE fields in TGCM geographic grid
  !     real,dimension(nlonp4,nlat) ::
  !    |  potg_sech, ekvg_sech, efxg_sech
  !     real,dimension(nmlonp1,-2:nlevp1) :: tiepot_sech
  !
  ! Define fields for AMIE input data file:
  ! electric potential in Volt
  ! mean energy in KeV
  ! energy flux in W/m^2
  ! amie_cusplat_nh(sh) and amie_cuspmlt_nh(sh) are
  !   AMIE cusp latitude and MLT in NH and SH
  ! amie_hpi_nh(sh) are AMIE hemi-integrated power
  ! amie_pcp_nh(sh) are AMIE polar-cap potential drop
  ! Saved AMIE outputs with suffix _amie
  !
  real(r8),allocatable,dimension(:,:,:),save :: & ! (lonp1,latp1,ntimes)
       amie_pot_nh, amie_pot_sh, amie_ekv_nh, amie_ekv_sh, & 
       amie_efx_nh, amie_efx_sh
  real(r8),allocatable,dimension(:,:),save ::  &  ! (lonp1,latp1)
       pot_nh_amie,pot_sh_amie, ekv_nh_amie,ekv_sh_amie, &
       efx_nh_amie,efx_sh_amie
  integer, allocatable,dimension(:),save :: & ! (ntimes)
       year,month,day,jday
  real(r8), allocatable,dimension(:),save :: & ! (ntimes)
       amie_cusplat_nh, amie_cuspmlt_nh, amie_hpi_nh, &
       amie_pcp_nh, amie_nh_ut, &
       amie_cusplat_sh, amie_cuspmlt_sh, amie_hpi_sh, &
       amie_pcp_sh, amie_sh_ut
  real(r8) :: &
       cusplat_nh_amie, cuspmlt_nh_amie, cusplat_sh_amie, &
       cuspmlt_sh_amie, hpi_sh_amie, hpi_nh_amie, pcp_sh_amie, &
       pcp_nh_amie
  !
#endif

contains
  !-----------------------------------------------------------------------
  subroutine init_amie(amienh,amiesh)
    !
    ! Called from tgcm.F
    ! (this is not in init.F to avoid circular dependencies)
    !
    character(len=*),intent(in) :: amienh, amiesh

#ifdef WACCMX_EDYN_ESMF
    ! read north hemisphere file:
    if (len_trim(amienh) > 0) then
       if (masterproc) write(iulog,"('Reading AMIENH file ',a)") trim(amienh)
       call rdamie_nh(amienh)
    end if
    !
    ! Read south hemisphere file:
    if (len_trim(amiesh) > 0) then
       if (masterproc) write(iulog,"('Reading AMIESH file ',a)") trim(amiesh)
       call rdamie_sh(amiesh)
    end if
#else
    call endrun('Cannot use AMIE without electro-dynamo active.')
#endif
  end subroutine init_amie
#ifdef WACCMX_EDYN_ESMF
  !-----------------------------------------------------------------------
  subroutine rdamie_nh(amienh)
    !
    ! Read AMIE data for the northern hemisphere from amienh
    !
    ! Local:

    character(len=*),intent(in) :: amienh
    integer :: istat,ntimes,ndims,nvars,ngatts,idunlim,ier
    integer :: id_lon,id_lat,id_time, &
         idv_year,idv_mon,idv_day,idv_jday, &
         idv_ut,idv_pot,idv_ekv, &
         idv_efx,idv_cusplat,idv_cuspmlt,idv_hpi,idv_pcp
    type(file_desc_t) :: ncid
    !
    if (masterproc) write(iulog,"(/,72('-'))")
    if (masterproc) write(iulog,"('RDAMIE_NH: read AMIE data for northern hemisphere:')")
    !
    ! Open netcdf file:
    call cam_pio_openfile(ncid, amienh, pio_nowrite)
    !
    ! Get AMIE grid dimension:
    istat = pio_inq_dimid(ncid,'lon',id_lon)
    istat = pio_inquire_dimension(ncid,id_lon,len=lonp1)
    if (istat /= pio_noerr) call rpt_ncerr(istat, 'rdamie_nh: Error getting AMIE longitude dimension')

    istat = pio_inq_dimid(ncid,'lat',id_lat)
    istat = pio_inquire_dimension(ncid,id_lat,len=latp1)
    if (istat /= pio_noerr) call rpt_ncerr(istat, 'rdamie_nh: Error getting AMIE latitude dimension')
    !     write(iulog,"('lonp1=',i3,' latp1=',i3)") lonp1,latp1
    !
    ! Get time dimension:
    istat = pio_inquire(ncid,unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid,id_time,len=ntimes)
    !
    ! Search for requested AMIE output fields
    istat = pio_inquire(ncid,ndims,nvars,ngatts,idunlim)
    !
    ! Get 1-D AMIE fields (ntimes)
    if (.not. allocated(year)) allocate(year(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'year',idv_year)
    istat = pio_get_var(ncid,idv_year,year)
    !     write(iulog,*)'rdamie_nh: year=', year(1:10)
    if (.not. allocated(month)) allocate(month(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'month',idv_mon)
    istat = pio_get_var(ncid,idv_mon,month)
    if (.not. allocated(day)) allocate(day(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'day',idv_day)
    istat = pio_get_var(ncid,idv_day,day)
    !     write(iulog,*)'rdamie_nh: day=', day(1:10)
    if (.not. allocated(jday)) allocate(jday(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'jday',idv_jday)
    istat = pio_get_var(ncid,idv_jday,jday)
    !
    ! Allocate 1-d fields:
    if (.not. allocated(amie_nh_ut)) &
         allocate(amie_nh_ut(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_nh_ut: ntimes=',i3)")ntimes
    if (.not. allocated(amie_cusplat_nh))  &
         allocate(amie_cusplat_nh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_cusplat_nh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_cuspmlt_nh))  &
         allocate(amie_cuspmlt_nh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_cuspmlt_nh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_hpi_nh)) &
         allocate(amie_hpi_nh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_hpi_nh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_pcp_nh)) &
         allocate(amie_pcp_nh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_pcp_nh: ntimes=',i3)")ntimes
    !
    ! Get ut
    istat = pio_inq_varid(ncid,'ut',idv_ut)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE UT id')
    istat = pio_get_var(ncid,idv_ut,amie_nh_ut)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable ut')
    !
    ! Get HPI
    istat = pio_inq_varid(ncid,'hpi',idv_hpi)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE hpi id')
    istat = pio_get_var(ncid,idv_hpi,amie_hpi_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid,'pcp',idv_pcp)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE pcp id')
    istat = pio_get_var(ncid,idv_pcp,amie_pcp_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid,'cusplat',idv_cusplat)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE cusplat id')
    istat = pio_get_var(ncid,idv_cusplat,amie_cusplat_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid,'cuspmlt',idv_cuspmlt)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE cusplat id')
    istat = pio_get_var(ncid,idv_cuspmlt,amie_cuspmlt_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable cuspmlt')
    !
    ! Allocate 2-d fields:
    if (.not. allocated(pot_nh_amie)) &
         allocate(pot_nh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' pot_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    if (.not. allocated(ekv_nh_amie)) &
         allocate(ekv_nh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' ekv_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    if (.not. allocated(efx_nh_amie)) &
         allocate(efx_nh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' efx_nh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    !
    ! Allocate 3-d fields:
    if (.not. allocated(amie_pot_nh))  &
         allocate(amie_pot_nh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) WRITE(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_pot_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    if (.not. allocated(amie_ekv_nh))  &
         allocate(amie_ekv_nh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_ekv_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    if (.not. allocated(amie_efx_nh))  &
         allocate(amie_efx_nh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_nh: error allocating',  &
         ' amie_efx_nh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    !
    ! Get 3-D AMIE fields (lon,lat,ntimes)
    !
    ! AMIE electric potential
    istat = pio_inq_varid(ncid,'pot',idv_pot)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE electric potential id')
    istat = pio_get_var(ncid,idv_pot,amie_pot_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable pot')
    !
    ! AMIE mean energy
    istat = pio_inq_varid(ncid,'ekv',idv_ekv)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE mean energy id')
    istat = pio_get_var(ncid,idv_ekv,amie_ekv_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable ekv')
    !
    ! AMIE energy flux
    istat = pio_inq_varid(ncid,'efx',idv_efx)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE energy flux id')
    istat = pio_get_var(ncid,idv_efx,amie_efx_nh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_nh: Error getting NH AMIE variable efx')
    !
    ! Close the file:
    call cam_pio_closefile(ncid)
    if (masterproc)  &
         write(iulog,"('Completed read from NH AMIE data file ',a)") trim(amienh)
    if (masterproc) write(iulog,"(72('-'),/)")
  end subroutine rdamie_nh
  !-----------------------------------------------------------------------
  subroutine rdamie_sh(amiesh)
    !
    ! Read AMIE data for the northern hemisphere from amiesh
    !
    ! Local:

    character(len=*),intent(in) :: amiesh
    integer :: istat,ntimes,ndims,nvars,ngatts,idunlim,ier
    integer :: id_lon,id_lat,id_time, &
         idv_year,idv_mon,idv_day,idv_jday, &
         idv_ut,idv_pot,idv_ekv, &
         idv_efx,idv_cusplat,idv_cuspmlt,idv_hpi,idv_pcp
    type(file_desc_t) :: ncid
    !
    if (masterproc) write(iulog,"(/,72('-'))")
    if (masterproc) write(iulog,"('RDAMIE_SH: read AMIE data for northern hemisphere:')")
    !
    ! Open netcdf file:
    call cam_pio_openfile(ncid, amiesh, pio_nowrite)
    !
    ! Get AMIE grid dimension:
    istat = pio_inq_dimid(ncid,'lon',id_lon)
    istat = pio_inquire_dimension(ncid,id_lon,len=lonp1)
    if (istat /= pio_noerr) call rpt_ncerr(istat, 'rdamie_sh: Error getting AMIE longitude dimension')

    istat = pio_inq_dimid(ncid,'lat',id_lat)
    istat = pio_inquire_dimension(ncid,id_lat,len=latp1)
    if (istat /= pio_noerr) call rpt_ncerr(istat, 'rdamie_sh: Error getting AMIE latitude dimension')
    !     write(iulog,"('lonp1=',i3,' latp1=',i3)") lonp1,latp1
    !
    ! Get time dimension:
    istat = pio_inquire(ncid,unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid,id_time,len=ntimes)
    !
    ! Search for requested AMIE output fields
    istat = pio_inquire(ncid,ndims,nvars,ngatts,idunlim)
    !
    ! Get 1-D AMIE fields (ntimes)
    if (.not. allocated(year)) allocate(year(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'year',idv_year)
    istat = pio_get_var(ncid,idv_year,year)
    !     write(iulog,*)'rdamie_sh: year=', year(1:10)
    if (.not. allocated(month)) allocate(month(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'month',idv_mon)
    istat = pio_get_var(ncid,idv_mon,month)
    if (.not. allocated(day)) allocate(day(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'day',idv_day)
    istat = pio_get_var(ncid,idv_day,day)
    !     write(iulog,*)'rdamie_sh: day=', day(1:10)
    if (.not. allocated(jday)) allocate(jday(ntimes),stat=ier)
    istat = pio_inq_varid(ncid,'jday',idv_jday)
    istat = pio_get_var(ncid,idv_jday,jday)
    !
    ! Allocate 1-d fields:
    if (.not. allocated(amie_sh_ut)) &
         allocate(amie_sh_ut(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_sh_ut: ntimes=',i3)")ntimes
    if (.not. allocated(amie_cusplat_sh))  &
         allocate(amie_cusplat_sh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_cusplat_sh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_cuspmlt_sh))  &
         allocate(amie_cuspmlt_sh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_cuspmlt_sh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_hpi_sh)) &
         allocate(amie_hpi_sh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_hpi_sh: ntimes=',i3)")ntimes
    if (.not. allocated(amie_pcp_sh)) &
         allocate(amie_pcp_sh(ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_pcp_sh: ntimes=',i3)")ntimes
    !
    ! Get ut
    istat = pio_inq_varid(ncid,'ut',idv_ut)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE UT id')
    istat = pio_get_var(ncid,idv_ut,amie_sh_ut)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable ut')
    !
    ! Get HPI
    istat = pio_inq_varid(ncid,'hpi',idv_hpi)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE hpi id')
    istat = pio_get_var(ncid,idv_hpi,amie_hpi_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid,'pcp',idv_pcp)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE pcp id')
    istat = pio_get_var(ncid,idv_pcp,amie_pcp_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid,'cusplat',idv_cusplat)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE cusplat id')
    istat = pio_get_var(ncid,idv_cusplat,amie_cusplat_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid,'cuspmlt',idv_cuspmlt)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE cusplat id')
    istat = pio_get_var(ncid,idv_cuspmlt,amie_cuspmlt_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable cuspmlt')
    !
    ! Allocate 2-d fields:
    if (.not. allocated(pot_sh_amie)) &
         allocate(pot_sh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' pot_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    if (.not. allocated(ekv_sh_amie)) &
         allocate(ekv_sh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' ekv_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    if (.not. allocated(efx_sh_amie)) &
         allocate(efx_sh_amie(lonp1,latp1),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' efx_sh_amie: lonp1=',i3,' latp1=',i3)")lonp1,latp1
    !
    ! Allocate 3-d fields:
    if (.not. allocated(amie_pot_sh))  &
         allocate(amie_pot_sh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_pot_sh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    if (.not. allocated(amie_ekv_sh))  &
         allocate(amie_ekv_sh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_ekv_sh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    if (.not. allocated(amie_efx_sh))  &
         allocate(amie_efx_sh(lonp1,latp1,ntimes),stat=ier)
    if (ier /= 0) write(iulog,"('>>> rdamie_sh: error allocating',  &
         ' amie_efx_sh: lonp1=',i3,' latp1=',i3,' ntimes=',i3)") &
         lonp1,latp1,ntimes
    !
    ! Get 3-D AMIE fields (lon,lat,ntimes)
    !
    ! AMIE electric potential
    istat = pio_inq_varid(ncid,'pot',idv_pot)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE electric potential id')
    istat = pio_get_var(ncid,idv_pot,amie_pot_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable pot')
    !
    ! AMIE mean energy
    istat = pio_inq_varid(ncid,'ekv',idv_ekv)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE mean energy id')
    istat = pio_get_var(ncid,idv_ekv,amie_ekv_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable ekv')
    !
    ! AMIE energy flux
    istat = pio_inq_varid(ncid,'efx',idv_efx)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE energy flux id')
    istat = pio_get_var(ncid,idv_efx,amie_efx_sh)
    if (istat /= pio_noerr) call rpt_ncerr(istat,  &
         'rdamie_sh: Error getting SH AMIE variable efx')
    !
    ! Close the file:
    call cam_pio_closefile(ncid)
    if (masterproc)  &
         write(iulog,"('Completed read from SH AMIE data file ',a)") trim(amiesh)
    if (masterproc) write(iulog,"(72('-'),/)")
  end subroutine rdamie_sh
#endif
  !-----------------------------------------------------------------------
  subroutine getamie(iyear,imo,iday,iutsec,sunlon,amie_ibkg,iprint,  &
                     iamie,phihm,amie_efxm,amie_kevm,crad, efxg,kevg)
    use cam_history_support, only: fillvalue
    use rgrd_mod, only: rgrd2
    !
    !     Read AMIE outputs from amie_ncfile file, returning electric potential,
    !     auroral mean energy and energy flux at current date and time,
    !     and the data is linearly interpolated to the model time
    !     gl - 12/07/2002
    !
    !
    !     Args:

    integer,  intent(in)    :: iyear
    integer,  intent(in)    :: imo
    integer,  intent(in)    :: iday
    real(r8), intent(in)    :: sunlon
    integer,  intent(in)    :: iutsec
    integer,  intent(in)    :: amie_ibkg
    integer,  intent(in)    :: iprint
    integer,  intent(inout) :: iamie
    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: amie_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: amie_kevm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: crad(2)
    real(r8), intent(out)   :: efxg(lon0:lon1,lat0:lat1) ! on geographic grid
    real(r8), intent(out)   :: kevg(lon0:lon1,lat0:lat1) ! on geographic grid
#ifdef WACCMX_EDYN_ESMF

    !
    !     Local:
    real(r8) :: potm(lonp1,jmxm),efxm(lonp1,jmxm),ekvm(lonp1,jmxm),  &
                alat(jmxm),alon(lonp1),alatm(jmxm),alonm(lonp1)
    integer :: ier,lw,liw,intpol(2)
    integer, allocatable :: iw(:)
    real(r8),allocatable :: w(:)
    integer :: i,j
    integer :: nn, iset, iset1, m, mp1, n
    integer :: iboxcar
    real(r8) :: model_ut, denoma, f1, f2
    real(r8) :: del,xmlt,dmlat,dlatm,dlonm,dmltm,rot,dtr,rtd
    integer :: idate,bdate,edate
    real(r8) :: pi

    pi = 4._r8*atan(1._r8)
    dtr = pi/180._r8          ! degrees to radians
    rtd = 180._r8/pi
    !     

    phihm = fillvalue
    amie_efxm = fillvalue
    amie_kevm = fillvalue
    efxg = fillvalue
    kevg = fillvalue
    crad = fillvalue

    !
    if (iprint > 0 .and. masterproc) then
       write(iulog,"(/,72('-'))")
       write(iulog,"('GETAMIE:')")
       write(iulog,"('Initial requested iyear=',i4,  &
            ' iday=',i3,' iutsec=', i10)") iyear,iday,iutsec
    end if

    !
    !     Check times:
    !
    nn = size(amie_sh_ut)
    bdate = year(1)*10000+month(1)*100+day(1)
    edate = year(nn)*10000+month(nn)*100+day(nn)
    idate = iyear*10000+imo*100+iday

    if (idate<bdate) then
       if (masterproc) write(iulog, "('getamie: Model date prior to AMIE first date:',3I5)") &
            year(1),month(1),day(1)
       iamie = 2
       return
    endif
    if (idate>edate) then
       if (masterproc) write(iulog, "('getamie: Model date beyond the AMIE last Data:',3I5)") &
            year(nn),month(nn),day(nn)
       iamie = 0
       return
    endif

    if (iamie/=1) return
    
    model_ut = dble(iutsec)/3600._r8

    !
    !     interpolate AMIE data to modeltime iutsec
    !     amie_ibkg = 0  use real UT AMIE data
    !     = 1  use the first AMIE volumne as the background
    !     = 2  use the 24-hr average AMIE volumne as the background
    pot_sh_amie(:,:) = 0._r8
    ekv_sh_amie(:,:) = 0._r8
    efx_sh_amie(:,:) = 0._r8
    cusplat_sh_amie = 0._r8
    cuspmlt_sh_amie = 0._r8
    hpi_sh_amie = 0._r8
    pcp_sh_amie = 0._r8
    !

    iboxcar = 0

    if (amie_ibkg == 0) then

       iset = nn
       iset1 = nn
       do i=1,nn
          !     if (amie_sh_ut(i) < model_ut) iset = i
          if (amie_sh_ut(i) < model_ut+(iday-day(i))*24._r8) iset = i
       end do
       !     write(iulog,"('getamie: AMIE SH Data nn,iset,day1,day2=',4i5)")
       !     |    nn,iset,jday(1),jday(nn)
       iset1 = iset + 1
       if (iset == nn) iset1 = iset

       denoma = amie_sh_ut(iset1) - amie_sh_ut(iset)
       if (denoma > 1._r8) then
          write(iulog, "('getamie: Finding a gap in the AMIE Data set:',  &
               'modelday, amieday =',2I5)") iday,day(n)
          iamie = 2
          return
       end if
       if (denoma == 0._r8) then
          f1 = 1._r8
          f2 = 0._r8
       else
          !     f1 = (amie_sh_ut(iset1) - model_ut)/denoma
          !     f2 = (model_ut - amie_sh_ut(iset))/denoma
          f1 = (amie_sh_ut(iset1) - (model_ut+(iday- &
               day(iset1))*24._r8))/denoma
          f2 = (model_ut+(iday-day(iset1))*24._r8 - &
               amie_sh_ut(iset))/denoma
       end if
       !     write(iulog,"('getamie: AMIE SH Data n,iset,modeltime,f1,f2 =',
       !     |    4i5,2f5.2)")n,iset,iday,day(iset1),f1,f2
       !     write(iulog,"('getamie: AMIE SH Data model_day,model_ut,amie_day,',
       !     |    'amie_ut,f1,f2,iset,iset1 =',i4,f7.1,i4,f7.1,2f5.2,2i3)")
       !     |    iday,model_ut,day(iset),amie_sh_ut(iset),f1,f2,
       !     |    iset,iset1
       cusplat_sh_amie = (f1*amie_cusplat_sh(iset1) + &
            f2*amie_cusplat_sh(iset))
       cuspmlt_sh_amie = (f1*amie_cuspmlt_sh(iset1) + &
            f2*amie_cuspmlt_sh(iset))
       hpi_sh_amie = (f1*amie_hpi_sh(iset1) + f2*amie_hpi_sh(iset))
       pcp_sh_amie = (f1*amie_pcp_sh(iset1) + f2*amie_pcp_sh(iset))
       if (iboxcar == 0) then
          pot_sh_amie(:,:) = (f1*amie_pot_sh(:,:,iset1) + &
               f2*amie_pot_sh(:,:,iset))
          ekv_sh_amie(:,:) = (f1*amie_ekv_sh(:,:,iset1) + &
               f2*amie_ekv_sh(:,:,iset))
          efx_sh_amie(:,:) = (f1*amie_efx_sh(:,:,iset1) + &
               f2*amie_efx_sh(:,:,iset))
       else
          call boxcar_ave(amie_pot_sh,pot_sh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
          call boxcar_ave(amie_efx_sh,efx_sh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
          call boxcar_ave(amie_ekv_sh,ekv_sh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
       end if
    else
       if (amie_ibkg == 1) then
          pot_sh_amie(:,:) = amie_pot_sh(:,:,1)
          ekv_sh_amie(:,:) = amie_ekv_sh(:,:,1)
          efx_sh_amie(:,:) = amie_efx_sh(:,:,1)
          cusplat_sh_amie = amie_cusplat_sh(1)
          cuspmlt_sh_amie = amie_cuspmlt_sh(1)
          hpi_sh_amie = amie_hpi_sh(1)
          pcp_sh_amie = amie_pcp_sh(1)
       else if (amie_ibkg == 3) then
          pot_sh_amie(:,:) = amie_pot_sh(:,:,241)
          ekv_sh_amie(:,:) = amie_ekv_sh(:,:,241)
          efx_sh_amie(:,:) = amie_efx_sh(:,:,241)
          cusplat_sh_amie = amie_cusplat_sh(241)
          cuspmlt_sh_amie = amie_cuspmlt_sh(241)
          hpi_sh_amie = amie_hpi_sh(241)
          pcp_sh_amie = amie_pcp_sh(241)
       else
          do i=1,nn
             pot_sh_amie(:,:) = pot_sh_amie(:,:) + amie_pot_sh(:,:,1)
             ekv_sh_amie(:,:) = ekv_sh_amie(:,:) + amie_ekv_sh(:,:,1)
             efx_sh_amie(:,:) = efx_sh_amie(:,:) + amie_efx_sh(:,:,1)
             cusplat_sh_amie = cusplat_sh_amie + amie_cusplat_sh(1)
             cuspmlt_sh_amie = cuspmlt_sh_amie + amie_cuspmlt_sh(1)
             hpi_sh_amie = hpi_sh_amie + amie_hpi_sh(1)
             pcp_sh_amie = pcp_sh_amie + amie_pcp_sh(1)
          end do
          pot_sh_amie(:,:) = pot_sh_amie(:,:)/nn
          ekv_sh_amie(:,:) = ekv_sh_amie(:,:)/nn
          efx_sh_amie(:,:) = efx_sh_amie(:,:)/nn
          cusplat_sh_amie = cusplat_sh_amie/nn
          cuspmlt_sh_amie = cuspmlt_sh_amie/nn
          hpi_sh_amie = hpi_sh_amie/nn
          pcp_sh_amie = pcp_sh_amie/nn
       end if
    end if

    !
    !     get NH AMIE data
    pot_nh_amie(:,:) = 0._r8
    ekv_nh_amie(:,:) = 0._r8
    efx_nh_amie(:,:) = 0._r8
    cusplat_nh_amie = 0._r8
    cuspmlt_nh_amie = 0._r8
    hpi_nh_amie = 0._r8
    pcp_nh_amie = 0._r8

    iboxcar = 0
    !     write(iulog,"('getamie: Interpolate AMIE NH Data nn=',i3)")nn
    if (amie_ibkg == 0) then
       iset = 0
       iset1 = nn
       do i=1,nn
          if (amie_nh_ut(i) < model_ut+(iday-day(i))*24._r8) iset = i
       end do
       iset1 = iset + 1
       if (iset == 0) iset = 1
       if (iset == nn) iset1 = iset

       denoma = amie_nh_ut(iset1) - amie_nh_ut(iset)
       if (denoma > 1._r8) then
          write(iulog, "('getamie: Finding a gap in the AMIE Data set:',  &
               'modelday, amieday =',2I5)") iday,day(n)
          iamie = 2
          return
       end if
       if (denoma == 0._r8) then
          f1 = 1._r8
          f2 = 0._r8
       else
          !     f1 = (amie_nh_ut(iset1) - model_ut)/denoma
          !     f2 = (model_ut - amie_nh_ut(iset))/denoma
          f1 = (amie_nh_ut(iset1) - (model_ut+(iday- &
               day(iset1))*24._r8))/denoma
          f2 = (model_ut+(iday-day(iset1))*24._r8 - &
               amie_nh_ut(iset))/denoma
       end if
       !     write(iulog,"('getamie: AMIE NH Data model_day,model_ut,amie_day,',
       !     |    'amie_ut,f1,f2,iset,iset1 =',i4,f7.1,i4,f7.1,2f5.2,2i3)")
       !     |    iday,model_ut,day(iset),amie_nh_ut(iset),f1,f2,
       !     |    iset,iset1
       !
       cusplat_nh_amie = (f1*amie_cusplat_nh(iset1) + &
            f2*amie_cusplat_nh(iset))
       cuspmlt_nh_amie = (f1*amie_cuspmlt_nh(iset1) + &
            f2*amie_cuspmlt_nh(iset))
       hpi_nh_amie = (f1*amie_hpi_nh(iset1) + f2*amie_hpi_nh(iset))
       pcp_nh_amie = (f1*amie_pcp_nh(iset1) + f2*amie_pcp_nh(iset))
       if (iboxcar == 0) then
          pot_nh_amie(:,:) = (f1*amie_pot_nh(:,:,iset1) + &
               f2*amie_pot_nh(:,:,iset))
          ekv_nh_amie(:,:) = (f1*amie_ekv_nh(:,:,iset1) + &
               f2*amie_ekv_nh(:,:,iset))
          efx_nh_amie(:,:) = (f1*amie_efx_nh(:,:,iset1) + &
               f2*amie_efx_nh(:,:,iset))
          !     write(iulog,"('ekv_nh_amie min, max = ',2e12.4)")
          !     |       minval(ekv_nh_amie),maxval(ekv_nh_amie)
       else
          call boxcar_ave(amie_pot_nh,pot_nh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
          !     call fminmax(amie_pot_nh(:,:,iset),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('AMIE pot max,min = ',2f8.0)")fmax,fmin
          !     call fminmax(pot_nh_amie(:,:),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('boxcar_ave AMIE pot max,min= ',2f8.0)")fmax,fmin
          call boxcar_ave(amie_efx_nh,efx_nh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
          !     call fminmax(amie_efx_nh(:,:,iset),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('AMIE efx max,min = ',2f8.0)")fmax,fmin
          !     call fminmax(efx_nh_amie(:,:),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('boxcar_ave AMIE efx max,min= ',2f8.0)")fmax,fmin
          call boxcar_ave(amie_ekv_nh,ekv_nh_amie,lonp1,latp1, &
               nn,iset,iboxcar)
          !     call fminmax(amie_ekv_nh(:,:,iset),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('AMIE ekv max,min = ',2f8.0)")fmax,fmin
          !     call fminmax(ekv_nh_amie(:,:),lonp1*latp1,fmin,fmax)
          !     write(iulog,"('boxcar_ave AMIE ekv max,min= ',2f8.0)")fmax,fmin
       end if
    else
       if (amie_ibkg == 1) then
          pot_nh_amie(:,:) = amie_pot_nh(:,:,1)
          ekv_nh_amie(:,:) = amie_ekv_nh(:,:,1)
          efx_nh_amie(:,:) = amie_efx_nh(:,:,1)
          cusplat_nh_amie = amie_cusplat_nh(1)
          cuspmlt_nh_amie = amie_cuspmlt_nh(1)
          hpi_nh_amie = amie_hpi_nh(1)
          pcp_nh_amie = amie_pcp_nh(1)
       else if (amie_ibkg == 3) then
          pot_nh_amie(:,:) = amie_pot_nh(:,:,241)
          ekv_nh_amie(:,:) = amie_ekv_nh(:,:,241)
          efx_nh_amie(:,:) = amie_efx_nh(:,:,241)
          cusplat_nh_amie = amie_cusplat_nh(241)
          cuspmlt_nh_amie = amie_cuspmlt_nh(241)
          hpi_nh_amie = amie_hpi_nh(241)
          pcp_nh_amie = amie_pcp_nh(241)
       else
          do i=1,nn
             pot_nh_amie(:,:) = pot_nh_amie(:,:) + amie_pot_nh(:,:,1)
             ekv_nh_amie(:,:) = ekv_nh_amie(:,:) + amie_ekv_nh(:,:,1)
             efx_nh_amie(:,:) = efx_nh_amie(:,:) + amie_efx_nh(:,:,1)
             cusplat_nh_amie = cusplat_nh_amie + amie_cusplat_nh(1)
             cuspmlt_nh_amie = cuspmlt_nh_amie + amie_cuspmlt_nh(1)
             hpi_nh_amie = hpi_nh_amie + amie_hpi_nh(1)
             pcp_nh_amie = pcp_nh_amie + amie_pcp_nh(1)
          end do
          pot_nh_amie(:,:) = pot_nh_amie(:,:)/nn
          ekv_nh_amie(:,:) = ekv_nh_amie(:,:)/nn
          efx_nh_amie(:,:) = efx_nh_amie(:,:)/nn
          cusplat_nh_amie = cusplat_nh_amie/nn
          cuspmlt_nh_amie = cuspmlt_nh_amie/nn
          hpi_nh_amie = hpi_nh_amie/nn
          pcp_nh_amie = pcp_nh_amie/nn
       end if
    end if
    !
    !     The OLTMAX latitude also defines the co-latitude theta0, which in
    !     turn determines crit1(+2.5deg) and crit2(-12.5deg) which are used
    !     in TIE-GCM as the boundaries of the polar cap and the region of
    !     influence of the high-lat potential versus the low-lat dynamo potential
    !     Define this latitude to be between 70 and 77.5 degrees
    !
    !     if (cusplat_sh_amie > 65.0) then
    !     cusplat_sh_amie = 65.0
    !     cuspmlt_sh_amie = 11.
    !     endif
    if (cusplat_sh_amie > 75.0_r8) then
       cusplat_sh_amie = 75.0_r8
       cuspmlt_sh_amie = 11._r8
    end if
    if (cusplat_sh_amie < 60.0_r8) then
       cusplat_sh_amie = 60.0_r8
       cuspmlt_sh_amie = 11._r8
    end if
    if (cusplat_nh_amie > 75.0_r8) then
       cusplat_nh_amie = 75.0_r8
       cuspmlt_nh_amie = 11._r8
    end if
    if (cusplat_nh_amie < 60.0_r8) then
       cusplat_nh_amie = 60.0_r8
       cuspmlt_nh_amie = 11._r8
    end if
    !     cusplat_nh_amie = amin1(65.0,cusplat_nh_amie)
    if (cuspmlt_sh_amie > 12.5_r8) cuspmlt_sh_amie = 12.5_r8
    if (cuspmlt_sh_amie < 11.0_r8) cuspmlt_sh_amie = 11.0_r8
    if (cuspmlt_nh_amie > 12.5_r8) cuspmlt_nh_amie = 12.5_r8
    if (cuspmlt_nh_amie < 11.0_r8) cuspmlt_nh_amie = 11.0_r8
    crad(1) = (90._r8-cusplat_sh_amie)*pi/180._r8
    crad(2) = (90._r8-cusplat_nh_amie)*pi/180._r8

    !     mlongitude starts from 180 degree
    rot = sunlon*rtd
    if(rot.lt.0) rot = rot + 360._r8    !  0 to 360 degrees
    rot = rot/15._r8                    !  convert from degree to hrs

    dmltm = 24._r8/dble(lonmx)
    do i=1,lonp1
       xmlt = dble(i-1)*dmltm - rot + 24._r8
       xmlt = MOD(xmlt,24._r8)
       m = int(xmlt/dmltm + 1.01_r8)
       mp1 = m + 1
       if (mp1 > lonp1) mp1 = 2
       del = xmlt - (m-1)*dmltm
       !     Initialize arrays around equator
       do j=latp1+1,ithmx
          potm(i,j) = 0._r8
          potm(i,jmxm+1-j) = 0._r8
          ekvm(i,j) = (1._r8-del)*ekv_sh_amie(m,latp1) + &
               del*ekv_sh_amie(mp1,latp1)
          ekvm(i,jmxm+1-j) = (1._r8-del)*ekv_nh_amie(m,latp1) +  &
               del*ekv_nh_amie(mp1,latp1)
          efxm(i,j) = 0._r8
          efxm(i,jmxm+1-j) = 0._r8
       end do
       !     Put in AMIE arrays from pole to latp1
       do j=1,latp1
          potm(i,j) = (1._r8-del)*pot_sh_amie(m,j) + &
               del*pot_sh_amie(mp1,j)
          potm(i,jmxm+1-j) = (1._r8-del)*pot_nh_amie(m,j) + &
               del*pot_nh_amie(mp1,j)
          ekvm(i,j) = (1._r8-del)*ekv_sh_amie(m,j) + &
               del*ekv_sh_amie(mp1,j)
          ekvm(i,jmxm+1-j) = (1._r8-del)*ekv_nh_amie(m,j) + &
               del*ekv_nh_amie(mp1,j)
          efxm(i,j) = (1._r8-del)*efx_sh_amie(m,j) + &
               del*efx_sh_amie(mp1,j)
          efxm(i,jmxm+1-j) = (1._r8-del)*efx_nh_amie(m,j) + &
               del*efx_nh_amie(mp1,j)
       end do

    end do

    !     Set up coeffs to go between EPOTM(IMXMP,JMNH) and TIEPOT(IMAXM,JMAXMH)

    !     ****     SET GRID SPACING DLATM, DLONG, DLONM
    !     DMLAT=lat spacing in degrees of AMIE apex grid
    dtr = pi/180._r8
    dmlat = 180._r8 / dble(jmxm-1)
    dlatm = dmlat*dtr
    dlonm = 2._r8*pi/dble(lonmx)
    dmltm = 24._r8/dble(lonmx)
    !     ****
    !     ****     SET ARRAY YLATM (LATITUDE VALUES FOR GEOMAGNETIC GRID
    !     ****
    alatm(1) = -pi/2._r8
    alat(1) = -90._r8
    alatm(jmxm) = pi/2._r8
    alat(jmxm) = 90._r8
    do i = 2,ithmx
       alat(i) = alat(i-1)+dlatm/dtr
       alat(jmxm+1-i) = alat(jmxm+2-i)-dlatm/dtr
       alatm(i) = alatm(i-1)+dlatm
       alatm(jmxm+1-i) = alatm(jmxm+2-i)-dlatm
    end do
    alon(1) = -pi/dtr
    alonm(1) = -pi
    do i=2,lonp1
       alon(i) = alon(i-1) + dlonm/dtr
       alonm(i) = alonm(i-1) + dlonm
    end do

    !     ylatm and ylonm are arrays of latitudes and longitudes of the
    !     distored magnetic grids in radian - from consdyn.h
    !     Convert from apex magnetic grid to distorted magnetic grid
    !
    !     Allocate workspace for regrid routine rgrd2.f:
    lw = nmlonp1+nmlat+2*nmlonp1
    if (.not. allocated(w)) allocate(w(lw),stat=ier)
    IF (ier /= 0) WRITE(iulog,"('>>> horizontal_interp: error allocating',  &
         ' w(lw): lw=',i6,' ier=',i4)") lw,ier
    liw = nmlonp1 + nmlat
    if (.not. allocated(iw)) allocate(iw(liw),stat=ier)
    if (ier /= 0) write(iulog,"('>>> horzontal_interp: error allocating',  &
         ' iw(liw): liw=',i6,' ier=',i4)") liw,ier
    intpol(:) = 1             ! linear (not cubic) interp in both dimensions
    if (alatm(1) > ylatm(1)) alatm(1) = ylatm(1)
    if (alatm(jmxm) < ylatm(nmlat)) alatm(jmxm) = ylatm(nmlat)
    if (alonm(1) > ylonm(1)) alonm(1) = ylonm(1)
    if (alonm(lonp1) < ylonm(nmlonp1)) alonm(lonp1) = ylonm(nmlonp1)
    !     write(iulog,"('  AMIE: ylatm =',/,(6e12.4))") ylatm
    !     write(iulog,"('  AMIE: ylonm =',/,(6e12.4))") ylonm
    !     write(iulog,"('  AMIE: potm(1,:) =',/,(6e12.4))") potm(1,:)
    !     ylatm from -pi/2 to pi/2, and ylonm from -pi to pi
    call rgrd2(lonp1,jmxm,alonm,alatm,potm,nmlonp1,nmlat,  &
         ylonm,ylatm,tiepot,intpol,w,lw,iw,liw,ier)
    call rgrd2(lonp1,jmxm,alonm,alatm,ekvm,nmlonp1,nmlat,  &
         ylonm,ylatm,tieekv,intpol,w,lw,iw,liw,ier)
    call rgrd2(lonp1,jmxm,alonm,alatm,efxm,nmlonp1,nmlat,  &
         ylonm,ylatm,tieefx,intpol,w,lw,iw,liw,ier)
    !     write(iulog,"('  AMIE: tiepot(1,:) =',/,(6e12.4))") tiepot(1,:)
    phihm(:,:) = tiepot(:,:)
    amie_efxm(:,:) = tieefx(:,:)
    amie_kevm(:,:) = tieekv(:,:)

    !     Convert from WACCM-X distorted magnetic grid to geographic one
    !     call mag2geo(tiepot(1,1),potg(1,0),im(1,0),jm(1,0),
    !     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
    !     call mag2geo(tieekv(1,1),ekvg(1,0),im(1,0),jm(1,0),
    !     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)
    !     call mag2geo(tieefx(1,1),efxg(1,0),im(1,0),jm(1,0),
    !     |    dim(1,0),djm(1,0),nlonp1,nmlonp1,nlon,nlat+2,nmlon,nmlat)

    call mag2geo_2d(amie_efxm(mlon0:mlon1,mlat0:mlat1),  &
         efxg, mag_efx,geo_efx,'MEFXAMIE')
    call mag2geo_2d(amie_kevm(mlon0:mlon1,mlat0:mlat1),  &
         kevg, mag_kev,geo_kev,'MKEVAMIE')

    !     call mag2geo_2d(amie_kevm,amie_kevg,mag_kev,geo_kev,'KEVM')
    if (iprint > 0 .and. masterproc) write(iulog,*) 'Max,min amie_efxm = ',  &
         maxval(amie_efxm),minval(amie_efxm)
    if (iprint > 0 .and. masterproc) write(iulog,*)  &
         'Max,min efxg = ',maxval(efxg),minval(efxg)
    !     ****
    !     ****     INSERT PERIODIC POINTS
    !     ****
    !     DO j = 1,nlat
    !     ekvg(nlonp1,j) = ekvg(1,j)
    !     efxg(nlonp1,j) = efxg(1,j)
    !     potg(nlonp1,j) = potg(1,j)
    !     ENDDO
    !
    if (iprint > 0 .and. masterproc) then
       write(iulog, "('getamie: AMIE data interpolated to date and time')")
       write(iulog,"('getamie: iyear,imo,iday,iutsec = ',3i6,i10)")  &
            iyear,imo,iday,iutsec
       write(iulog,"('getamie: AMIE iset f1,f2,year,mon,day,ut = ',  &
            i6,2F9.5,3I6,f10.4)")  &
            iset,f1,f2,year(iset),month(iset),day(iset),amie_nh_ut(iset)
       write(iulog,*)'getamie: max,min phihm= ', maxval(phihm),minval(phihm)
       !     write(iulog,*)'getamie: max,min phihm,amie_efx,amie_kev = ',
       !     |    maxval(phihm),minval(tiepot),maxval(amie_efx),
       !     |    minval(amie_efx),maxval(amie_kev),minval(amie_kev)
    end if
#else
    call endrun('Cannot use AMIE without electro-dynamo active.')
#endif
  end subroutine getamie
#ifdef WACCMX_EDYN_ESMF
  !-------------------------------------------------------------------
  subroutine boxcar_ave(x,y,lon,lat,mtime,itime,ibox)
    !
    ! perform boxcar average
    !
    ! Args:
    integer,  intent(in)  :: lon
    integer,  intent(in)  :: lat
    integer,  intent(in)  :: mtime
    integer,  intent(in)  :: itime
    integer,  intent(in)  :: ibox
    real(r8), intent(in)  :: x(lon,lat,mtime)
    real(r8), intent(out) :: y(lon,lat)
    
    ! Local:
    integer :: i, iset, iset1
    !
    iset = itime - ibox/2
    if (iset < 1) iset = 1
    iset1 = iset + ibox
    if (iset1 > mtime) then
       iset1 = mtime
       iset = iset1 - ibox
    end if
    !     write(iulog,"('boxcar_ave: mtime,itime,ibox',3i5)")
    !    |  mtime,itime,ibox
    !
    y(:,:) = 0._r8
    do i=iset,iset1
       y(:,:) = y(:,:) + x(:,:,i)
    end do
    if (ibox > 0) y(:,:) = y(:,:)/ibox
    !
  end subroutine boxcar_ave
  !-----------------------------------------------------------------------
  subroutine mag2geo(am,ag,im,jm,dim,djm,lg,lm,nlong,nlatg)
    !
    ! Args:
    integer,  intent(in)  :: lg
    integer,  intent(in)  :: lm
    real(r8), intent(in)  :: am(lm,*)
    real(r8), intent(out) :: ag(lg,*)
    integer,  intent(in)  :: im(lg,*)
    integer,  intent(in)  :: jm(lg,*)
    real(r8), intent(in)  :: dim(lg,*)
    real(r8), intent(in)  :: djm(lg,*)
    integer,  intent(in)  :: nlong
    integer,  intent(in)  :: nlatg
    !
    ! Local:
    integer :: ig,jg
    !
    do jg=1,nlatg
       do ig=1,nlong
          ag(ig,jg) =  &
               am(im(ig,jg)  ,jm(ig,jg))  *(1._r8-dim(ig,jg))*(1._r8-djm(ig,jg))+  &
               am(im(ig,jg)+1,jm(ig,jg))  *    dim(ig,jg) *(1._r8-djm(ig,jg))+  &
               am(im(ig,jg)  ,jm(ig,jg)+1)*(1._r8-dim(ig,jg))*djm(ig,jg)+  &
               am(im(ig,jg)+1,jm(ig,jg)+1)*    dim(ig,jg) *djm(ig,jg)
       end do ! ig=1,nlong
    end do ! jg=1,nlatg
  end subroutine mag2geo
  !-----------------------------------------------------------------------
  subroutine mag2geo_2d(fmag,fgeo,ESMF_mag,ESMF_geo,fname)
    !
    ! Convert field on geomagnetic grid fmag to geographic grid in fgeo.
    !
    use edyn_esmf,only: edyn_esmf_set2d_mag,edyn_esmf_regrid,  &
                        edyn_esmf_get_2dfield
    !
    ! Args:
    real(r8),         intent(in)    :: fmag(mlon0:mlon1,mlat0:mlat1)
    real(r8),         intent(out)   :: fgeo(lon0:lon1,lat0:lat1)
    type(ESMF_Field), intent(inout) :: ESMF_mag, ESMF_geo
    character(len=*), intent(in)    :: fname
    !
    ! Local:
    integer :: j
    character (len=8) :: fnames(1)
    type(ESMF_Field) :: magfields(1)
    real(r8),pointer,dimension(:,:) :: fptr

    fgeo = finit
    fnames(1) = fname
    magfields(1) = ESMF_mag
    !
    ! Put fmag into ESMF mag field on mag source grid:
    call edyn_esmf_set2d_mag(magfields,fnames,fmag,1, &
         mlon0,mlon1,mlat0,mlat1)
    !
    ! Regrid to geographic destination grid, defining ESMF_geo:
    call edyn_esmf_regrid(ESMF_mag,ESMF_geo,'mag2geo',2)
    !
    ! Put regridded geo field into pointer:
    call edyn_esmf_get_2dfield(ESMF_geo,fptr,fname)
    !      write(iulog,*) 'mag2geo: Max,min fptr = ',maxval(fptr),minval(fptr)
    !
    ! Transfer from pointer to output arg:
    do j=lat0,lat1
       fgeo(:,j) = fptr(:,j)
    end do
    !      write(iulog,*) 'mag2geo: max,min fmag = ',maxval(fmag),minval(fmag)
    !      write(iulog,*) 'mag2geo: max,min fgeo = ',maxval(fgeo),minval(fgeo)
  end subroutine mag2geo_2d
  !-----------------------------------------------------------------------
  subroutine rpt_ncerr(istat,msg)
    !
    ! Handle a netcdf lib error:
    !
    integer,         intent(in) :: istat
    character(len=*),intent(in) :: msg
    !
    write(iulog,"(/72('-'))")
    write(iulog,"('>>> Error from netcdf library:')")
    write(iulog,"(a)") trim(msg)

    write(iulog,"('istat=',i5)") istat
    write(iulog,"(72('-')/)")
    return
  end subroutine rpt_ncerr

#endif

end module amie_module
