module amie_module
  !
  ! Module used to read data from the AMIE outputs (POT,mean energy,
  !   and energy flux).
  !

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc
  use edyn_maggrid,   only: nmlat, nmlonp1
#ifdef WACCMX_EDYN_ESMF
  use edyn_maggrid,   only: ylonm     ! magnetic latitudes (nmlat) (radians)
  use edyn_maggrid,   only: ylatm     ! magnetic longtitudes (nmlonp1) (radians)
  use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile
  use pio,            only: pio_inq_dimid, pio_inquire_dimension
  use pio,            only: pio_inquire, pio_inq_varid
  use pio,            only: file_desc_t, pio_noerr, pio_nowrite, pio_get_var
  use utils_mod,      only: check_ncerr, check_alloc, boxcar_ave
  use edyn_mpi,       only: ntask, mytid
#else
  use cam_abortutils, only: endrun
#endif
  implicit none

  private
  public :: init_amie
  public :: getamie
#ifdef WACCMX_EDYN_ESMF

  ! Define parameters for AMIE input data file:
  integer, parameter ::  &
       ithmx = 55,       & ! maximum number of latitudes of AMIE data
       jmxm = 2*ithmx-1, & ! maximum number of global latitudes
       lonmx = 36          ! maximum number of longitudes of AMIE data
  integer :: lonp1,latp1
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
  real(r8), allocatable, dimension(:,:,:), save :: & ! (lonp1,latp1,ntimes)
       amie_pot_nh, amie_pot_sh,                   &
       amie_ekv_nh, amie_ekv_sh,                   &
       amie_efx_nh, amie_efx_sh
  real(r8), allocatable, dimension(:,:), save ::           &  ! (lonp1,latp1)
       pot_nh_amie, pot_sh_amie, ekv_nh_amie, ekv_sh_amie, &
       efx_nh_amie, efx_sh_amie
  integer,  allocatable, dimension(:), save ::                 & ! (ntimes)
       year, month, day, jday
  real(r8), allocatable, dimension(:), save ::                 & ! (ntimes)
       amie_cusplat_nh, amie_cuspmlt_nh, amie_hpi_nh,          &
       amie_pcp_nh, amie_nh_ut,                                &
       amie_cusplat_sh, amie_cuspmlt_sh, amie_hpi_sh,          &
       amie_pcp_sh, amie_sh_ut
  real(r8) ::                                                  &
       cusplat_nh_amie, cuspmlt_nh_amie, cusplat_sh_amie,      &
       cuspmlt_sh_amie, hpi_sh_amie, hpi_nh_amie, pcp_sh_amie, &
       pcp_nh_amie
  !
  type(file_desc_t) :: ncid_nh
  type(file_desc_t) :: ncid_sh

  character(len=256), allocatable :: amienh_files(:)
  character(len=256), allocatable :: amiesh_files(:)
  integer :: num_files, file_ndx
#endif

contains
  !-----------------------------------------------------------------------
  subroutine init_amie(amienh_list,amiesh_list)

    character(len=*),intent(in) :: amienh_list(:)
    character(len=*),intent(in) :: amiesh_list(:)

#ifdef WACCMX_EDYN_ESMF
    integer :: n, nfiles

    nfiles = min( size(amienh_list), size(amiesh_list) )
    num_files = 0

    count_files: do n = 1,nfiles
       if (len_trim(amienh_list(n))<1 .or. len_trim(amiesh_list(n))<1 .or. &
           trim(amienh_list(n))=='NONE' .or. trim(amiesh_list(n))=='NONE') then
          exit count_files
       else
          num_files = num_files + 1
       end if
    end do count_files
       
    allocate(amienh_files(num_files), amiesh_files(num_files))
    amienh_files(:num_files) = amienh_list(:num_files)
    amiesh_files(:num_files) = amiesh_list(:num_files)
    file_ndx = 1
    call open_files()
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

    ! Dummy argument
    character(len=*), intent(in) :: amienh
    ! Local variables:
    integer                      :: istat, ntimes, ndims, nvars, ngatts
    integer                      :: idunlim, ier
    integer                      :: id_lon, id_lat, id_time
    integer                      :: idv_year, idv_mon, idv_day, idv_jday
    integer                      :: idv_ut, idv_cusplat, idv_cuspmlt
    integer                      :: idv_hpi, idv_pcp
    character(len=*), parameter  :: subname = 'rdamie_nh'
    !
    if (masterproc) then
       write(iulog, "(/,72('-'))")
       write(iulog, "(a,': read AMIE data for northern hemisphere:')") subname
    end if
    !
    ! Open netcdf file:
    call cam_pio_openfile(ncid_nh, amienh, pio_nowrite)
    !
    ! Get AMIE grid dimension:
    istat = pio_inq_dimid(ncid_nh, 'lon', id_lon)
    istat = pio_inquire_dimension(ncid_nh, id_lon, len=lonp1)
    call check_ncerr(istat, subname, 'AMIE longitude dimension')

    istat = pio_inq_dimid(ncid_nh, 'lat', id_lat)
    istat = pio_inquire_dimension(ncid_nh, id_lat, len=latp1)
    call check_ncerr(istat, subname, 'AMIE latitude dimension')
    !     write(iulog, "('lonp1=', i3, ' latp1=', i3)") lonp1, latp1
    !
    ! Get time dimension:
    istat = pio_inquire(ncid_nh, unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid_nh, id_time, len=ntimes)
    !
    ! Search for requested AMIE output fields
    istat = pio_inquire(ncid_nh, ndims, nvars, ngatts, idunlim)
    !
    ! Get 1-D AMIE fields (ntimes)
    if (.not. allocated(year)) then
       allocate(year(ntimes), stat=ier)
       call check_alloc(ier, subname, 'year', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_nh, 'year', idv_year)
    call check_ncerr(istat, subname, 'AMIE year id')
    istat = pio_get_var(ncid_nh, idv_year, year)
    call check_ncerr(istat, subname, 'AMIE year')
    !     write(iulog, *)'rdamie_nh: year=', year(1:10)
    if (.not. allocated(month)) then
       allocate(month(ntimes), stat=ier)
       call check_alloc(ier, subname, 'month', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_nh, 'month', idv_mon)
    call check_ncerr(istat, subname, 'AMIE month id')
    istat = pio_get_var(ncid_nh, idv_mon, month)
    call check_ncerr(istat, subname, 'AMIE month')
    if (.not. allocated(day)) then
       allocate(day(ntimes), stat=ier)
       call check_alloc(ier, subname, 'day', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_nh, 'day', idv_day)
    call check_ncerr(istat, subname, 'AMIE day id')
    istat = pio_get_var(ncid_nh, idv_day, day)
    call check_ncerr(istat, subname, 'AMIE day')
    !     write(iulog, *)'rdamie_nh: day=', day(1:10)
    if (.not. allocated(jday)) then
       allocate(jday(ntimes), stat=ier)
       call check_alloc(ier, subname, 'jday', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_nh, 'jday', idv_jday)
    call check_ncerr(istat, subname, 'AMIE jday id')
    istat = pio_get_var(ncid_nh, idv_jday, jday)
    call check_ncerr(istat, subname, 'AMIE jday')
    !
    ! Allocate 1-d fields:
    if (.not. allocated(amie_nh_ut)) then
       allocate(amie_nh_ut(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_nh_ut', ntimes=ntimes)
    end if
    if (.not. allocated(amie_cusplat_nh)) then
       allocate(amie_cusplat_nh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_cusplat_nh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_cuspmlt_nh)) then
       allocate(amie_cuspmlt_nh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_cuspmlt_nh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_hpi_nh)) then
       allocate(amie_hpi_nh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_hpi_nh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_pcp_nh)) then
       allocate(amie_pcp_nh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_pcp_nh', ntimes=ntimes)
    end if
    !
    ! Get ut
    istat = pio_inq_varid(ncid_nh, 'ut', idv_ut)
    call check_ncerr(istat, subname, 'AMIE ut id')
    istat = pio_get_var(ncid_nh, idv_ut, amie_nh_ut)
    call check_ncerr(istat, subname, 'AMIE ut')
    !
    ! Get HPI
    istat = pio_inq_varid(ncid_nh, 'hpi', idv_hpi)
    call check_ncerr(istat, subname, 'AMIE hpi id')
    istat = pio_get_var(ncid_nh, idv_hpi, amie_hpi_nh)
    call check_ncerr(istat, subname, 'AMIE hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid_nh, 'pcp', idv_pcp)
    call check_ncerr(istat, subname, 'AMIE pcp id')
    istat = pio_get_var(ncid_nh, idv_pcp, amie_pcp_nh)
    call check_ncerr(istat, subname, 'AMIE pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid_nh, 'cusplat', idv_cusplat)
    call check_ncerr(istat, subname, 'AMIE cusplat id')
    istat = pio_get_var(ncid_nh, idv_cusplat, amie_cusplat_nh)
    call check_ncerr(istat, subname, 'AMIE cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid_nh, 'cuspmlt', idv_cuspmlt)
    call check_ncerr(istat, subname, 'AMIE cuspmlt id')
    istat = pio_get_var(ncid_nh, idv_cuspmlt, amie_cuspmlt_nh)
    call check_ncerr(istat, subname, 'AMIE cuspmlt')
    !
    ! Allocate 2-d fields:
    if (.not. allocated(pot_nh_amie)) then
       allocate(pot_nh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'pot_nh_amie', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(ekv_nh_amie)) then
       allocate(ekv_nh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'ekv_nh_amie', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(efx_nh_amie)) then
       allocate(efx_nh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'efx_nh_amie', lonp1=lonp1, latp1=latp1)
    end if
    !
    ! Allocate 3-d fields:
    if (.not. allocated(amie_pot_nh)) then
       allocate(amie_pot_nh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_pot_nh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(amie_ekv_nh)) then
       allocate(amie_ekv_nh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_ekv_nh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(amie_efx_nh)) then
       allocate(amie_efx_nh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_efx_nh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
  end subroutine rdamie_nh

  !-----------------------------------------------------------------------
  subroutine rdamie_sh(amiesh)
    !
    ! Read AMIE data for the southern hemisphere from amiesh
    !
    character(len=*), intent(in) :: amiesh
    ! Local:
    integer                      :: istat, ntimes, ndims, nvars, ngatts, ier
    integer                      :: idunlim
    integer                      :: id_lon, id_lat, id_time
    integer                      :: idv_year, idv_mon, idv_day, idv_jday
    integer                      :: idv_ut
    integer                      :: idv_cusplat, idv_cuspmlt, idv_hpi, idv_pcp
    character(len=*), parameter  :: subname = 'rdamie_sh'
    !
    if (masterproc) then
       write(iulog, "(/, 72('-'))")
       write(iulog, "(a, ': read AMIE data for southern hemisphere:')") subname
    end if
    !
    ! Open netcdf file:
    call cam_pio_openfile(ncid_sh, amiesh, pio_nowrite)
    !
    ! Get AMIE grid dimension:
    istat = pio_inq_dimid(ncid_sh, 'lon', id_lon)
    istat = pio_inquire_dimension(ncid_sh, id_lon, len=lonp1)
    call check_ncerr(istat, subname, 'AMIE longitude dimension')

    istat = pio_inq_dimid(ncid_sh, 'lat', id_lat)
    istat = pio_inquire_dimension(ncid_sh, id_lat, len=latp1)
    call check_ncerr(istat, subname, 'AMIE latitude dimension')
    !     write(iulog, "('lonp1=', i3, ' latp1=', i3)") lonp1, latp1
    !
    ! Get time dimension:
    istat = pio_inquire(ncid_sh, unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid_sh, id_time, len=ntimes)
    !
    ! Search for requested AMIE output fields
    istat = pio_inquire(ncid_sh, ndims, nvars, ngatts, idunlim)
    !
    ! Get 1-D AMIE fields (ntimes)
    if (.not. allocated(year)) then
       allocate(year(ntimes), stat=ier)
       call check_alloc(ier, subname, 'year', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_sh, 'year', idv_year)
    call check_ncerr(istat, subname, 'AMIE year id')
    istat = pio_get_var(ncid_sh, idv_year, year)
    call check_ncerr(istat, subname, 'AMIE year')
    !     write(iulog,*)'rdamie_sh: year=', year(1:10)
    if (.not. allocated(month)) then
       allocate(month(ntimes), stat=ier)
       call check_alloc(ier, subname, 'month', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_sh, 'month', idv_mon)
    call check_ncerr(istat, subname, 'AMIE month id')
    istat = pio_get_var(ncid_sh, idv_mon, month)
    call check_ncerr(istat, subname, 'AMIE month')
    if (.not. allocated(day)) then
       allocate(day(ntimes), stat=ier)
       call check_alloc(ier, subname, 'day', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_sh, 'day', idv_day)
    call check_ncerr(istat, subname, 'AMIE day id')
    istat = pio_get_var(ncid_sh, idv_day, day)
    call check_ncerr(istat, subname, 'AMIE day')
    !     write(iulog,*)'rdamie_sh: day=', day(1:10)
    if (.not. allocated(jday)) then
       allocate(jday(ntimes), stat=ier)
       call check_alloc(ier, subname, 'jday', ntimes=ntimes)
    end if
    istat = pio_inq_varid(ncid_sh, 'jday', idv_jday)
    call check_ncerr(istat, subname, 'AMIE jday id')
    istat = pio_get_var(ncid_sh, idv_jday, jday)
    call check_ncerr(istat, subname, 'AMIE jday')
    !
    ! Allocate 1-d fields:
    if (.not. allocated(amie_sh_ut)) then
       allocate(amie_sh_ut(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_sh_ut', ntimes=ntimes)
    end if
    if (.not. allocated(amie_cusplat_sh)) then
       allocate(amie_cusplat_sh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_cusplat_sh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_cuspmlt_sh)) then
       allocate(amie_cuspmlt_sh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_cuspmlt_sh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_hpi_sh)) then
       allocate(amie_hpi_sh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_hpi_sh', ntimes=ntimes)
    end if
    if (.not. allocated(amie_pcp_sh)) then
       allocate(amie_pcp_sh(ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_pcp_sh', ntimes=ntimes)
    end if
    !
    ! Get ut
    istat = pio_inq_varid(ncid_sh, 'ut', idv_ut)
    call check_ncerr(istat, subname, 'AMIE ut id')
    istat = pio_get_var(ncid_sh, idv_ut, amie_sh_ut)
    call check_ncerr(istat, subname, 'AMIE ut')
    !
    ! Get HPI
    istat = pio_inq_varid(ncid_sh, 'hpi', idv_hpi)
    call check_ncerr(istat, subname, 'AMIE hpi id')
    istat = pio_get_var(ncid_sh, idv_hpi, amie_hpi_sh)
    call check_ncerr(istat, subname, 'AMIE hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid_sh, 'pcp', idv_pcp)
    call check_ncerr(istat, subname, 'AMIE pcp id')
    istat = pio_get_var(ncid_sh, idv_pcp, amie_pcp_sh)
    call check_ncerr(istat, subname, 'AMIE pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid_sh, 'cusplat', idv_cusplat)
    call check_ncerr(istat, subname, 'AMIE cusplat id')
    istat = pio_get_var(ncid_sh, idv_cusplat, amie_cusplat_sh)
    call check_ncerr(istat, subname, 'AMIE cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid_sh, 'cuspmlt', idv_cuspmlt)
    call check_ncerr(istat, subname, 'AMIE cuspmlt id')
    istat = pio_get_var(ncid_sh, idv_cuspmlt, amie_cuspmlt_sh)
    call check_ncerr(istat, subname, 'AMIE cuspmlt')
    !
    ! Allocate 2-d fields:
    if (.not. allocated(pot_sh_amie)) then
       allocate(pot_sh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'pot_sh_amie', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(ekv_sh_amie)) then
       allocate(ekv_sh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'ekv_sh_amie', lonp1=lonp1, latp1=latp1)
    end if
    if (.not. allocated(efx_sh_amie)) then
       allocate(efx_sh_amie(lonp1, latp1), stat=ier)
       call check_alloc(ier, subname, 'efx_sh_amie', lonp1=lonp1, latp1=latp1)
    end if
    !
    ! Allocate 3-d fields:
    if (.not. allocated(amie_pot_sh)) then
       allocate(amie_pot_sh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_pot_sh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(amie_ekv_sh)) then
       allocate(amie_ekv_sh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_ekv_sh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(amie_efx_sh)) then
       allocate(amie_efx_sh(lonp1, latp1, ntimes), stat=ier)
       call check_alloc(ier, subname, 'amie_efx_sh', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
  end subroutine rdamie_sh

  !-----------------------------------------------------------------------
  subroutine update_3d_fields( ncid, offset, kount, pot_3d,ekv_3d,efx_3d )

    type(file_desc_t), intent(in) :: ncid
    integer, intent(in)  :: offset(:)
    integer, intent(in)  :: kount(:)
    real(r8),intent(out) :: pot_3d(:,:,:)
    real(r8),intent(out) :: ekv_3d(:,:,:)
    real(r8),intent(out) :: efx_3d(:,:,:)

    
    integer :: istat
    integer :: idv_pot, idv_ekv, idv_efx
    character(len=*), parameter :: subname = 'update_3d_fields'

    !
    ! Get 3-D fields (lon,lat,ntimes)
    !
    ! electric potential
    istat = pio_inq_varid(ncid, 'pot', idv_pot)
    call check_ncerr(istat, subname, 'AMIE pot id')
    istat = pio_get_var(ncid, idv_pot, offset, kount, pot_3d)
    call check_ncerr(istat, subname, 'AMIE pot')
    !
    ! mean energy
    istat = pio_inq_varid(ncid, 'ekv', idv_ekv)
    call check_ncerr(istat, subname, 'AMIE ekv id')
    istat = pio_get_var(ncid, idv_ekv, offset, kount, ekv_3d)
    call check_ncerr(istat, subname, 'AMIE ekv')
    !
    ! energy flux
    istat = pio_inq_varid(ncid, 'efx', idv_efx)
    call check_ncerr(istat, subname, 'AMIE efx id')
    istat = pio_get_var(ncid, idv_efx, offset, kount, efx_3d)
    call check_ncerr(istat, subname, 'AMIE efx')

  end subroutine update_3d_fields
#endif
  !-----------------------------------------------------------------------
  subroutine getamie(iyear, imo, iday, iutsec, sunlon, iprint,  &
                     iamie, phihm, amie_efxm, amie_kevm, crad)
    use cam_history_support, only: fillvalue
    use rgrd_mod,            only: rgrd2

    !
    !    Read AMIE outputs from amie_ncfile file, returning electric potential,
    !    auroral mean energy and energy flux at current date and time,
    !    and the data is linearly interpolated to the model time
    !    gl - 12/07/2002
    !
    !
    !    Args:

    integer,  intent(in)    :: iyear
    integer,  intent(in)    :: imo
    integer,  intent(in)    :: iday
    real(r8), intent(in)    :: sunlon
    integer,  intent(in)    :: iutsec
    integer,  intent(in)    :: iprint
    integer,  intent(out)   :: iamie
    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: amie_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: amie_kevm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: crad(2)
#ifdef WACCMX_EDYN_ESMF
    !     
    !
    !     Local:
    real(r8)                    :: potm(lonp1,jmxm)
    real(r8)                    :: efxm(lonp1,jmxm), ekvm(lonp1,jmxm)
    real(r8)                    :: alat(jmxm), alon(lonp1)
    real(r8)                    :: alatm(jmxm), alonm(lonp1)
    integer                     :: ier, lw, liw, intpol(2)
    integer,  allocatable       :: iw(:)
    real(r8), allocatable       :: w(:)
    integer                     :: i, j
    integer                     :: nn, iset, iset1, m, mp1, n
    integer                     :: iboxcar
    integer                     :: idate, bdate, edate
    real(r8)                    :: model_ut, denoma, f1, f2
    real(r8)                    :: del, xmlt, dmlat, dlatm, dlonm, dmltm, rot
    real(r8)                    :: pi, dtr, rtd
    integer                     :: offset(3), kount(3)
    character(len=*), parameter :: subname = 'getamie'

    pi = 4._r8 * atan(1._r8)
    dtr = pi / 180._r8          ! degrees to radians
    rtd = 180._r8 / pi          ! radians to degrees
    !

    phihm = fillvalue
    amie_efxm = fillvalue
    amie_kevm = fillvalue
    crad = fillvalue

    !
    if (iprint > 0 .and. masterproc) then
       write(iulog,"(/,72('-'))")
       write(iulog,"(a,':')") subname
       write(iulog,"(a,i4,', iday = ',i3,', iutsec = ',i10)")                &
            'Initial requested iyear= ', iyear, iday, iutsec
    end if

    nn = size(amie_sh_ut)
    bdate = year(1)*10000+month(1)*100+day(1)
    edate = year(nn)*10000+month(nn)*100+day(nn)
    idate = iyear*10000+imo*100+iday
    !
    !     Check times:
    !
    iamie=-1
    check_loop: do while( iamie/=1 )
       if (masterproc) write(iulog,*) 'file_ndx = ',file_ndx

       iamie = 1

       if (idate<bdate) then
          if (masterproc) then
             write(iulog, "(a,': Model date prior to AMIE first date:',3I5)")   &
                  subname, year(1), month(1), day(1)
          end if
          iamie = 2
          return
       endif
       if (idate>edate) then
          if (masterproc) then
             write(iulog, "(a,': Model date beyond the AMIE last Data:',3I5)")  &
                  subname, year(nn), month(nn), day(nn)
          end if
          iamie = 0

          if (file_ndx<num_files) then
             file_ndx = file_ndx+1
             call close_files()
             call open_files()
          else
             return
          end if

       endif
    end do check_loop

    model_ut = real(iutsec, kind=r8) / 3600._r8

!     get SH AMIE data
    pot_sh_amie(:,:) = 0._r8
    ekv_sh_amie(:,:) = 0._r8
    efx_sh_amie(:,:) = 0._r8
    cusplat_sh_amie = 0._r8
    cuspmlt_sh_amie = 0._r8
    hpi_sh_amie = 0._r8
    pcp_sh_amie = 0._r8
    !

    iboxcar = 0
    iset = 0
    iset1 = nn
    do i=1,nn
       if (amie_sh_ut(i) < model_ut+(iday-day(i))*24._r8) iset = i
    end do
    !     write(iulog,"('getamie: AMIE SH Data nn,iset,day1,day2=',4i5)")
    !     |    nn,iset,jday(1),jday(nn)
    if (iset == 0) iset = 1
    if (iset == nn) iset = nn-1
    iset1 = iset + 1

    denoma = amie_sh_ut(iset1) - amie_sh_ut(iset)
    if (denoma > 1._r8) then
       write(iulog, "(a,2(a,i0))")                                        &
            subname, ': Finding a gap in the AMIE Data set: modelday = ', &
            iday, ', amieday = ', day(n)
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

    offset = (/1,1,iset/)
    kount = (/lonp1,latp1,2/)

    call update_3d_fields( ncid_sh, offset, kount, amie_pot_sh,amie_ekv_sh,amie_efx_sh )
    if (iboxcar == 0) then       
       pot_sh_amie(:,:) = (f1*amie_pot_sh(:,:,2) + &
            f2*amie_pot_sh(:,:,1))
       ekv_sh_amie(:,:) = (f1*amie_ekv_sh(:,:,2) + &
            f2*amie_ekv_sh(:,:,1))
       efx_sh_amie(:,:) = (f1*amie_efx_sh(:,:,2) + &
            f2*amie_efx_sh(:,:,1))
    else
       call boxcar_ave(amie_pot_sh,pot_sh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
       call boxcar_ave(amie_efx_sh,efx_sh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
       call boxcar_ave(amie_ekv_sh,ekv_sh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
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
    iset = 0
    iset1 = nn

    do i=1,nn
       if (amie_nh_ut(i) < model_ut+(iday-day(i))*24._r8) iset = i
    end do
    if (iset == 0) iset = 1
    if (iset == nn) iset = nn-1
    iset1 = iset + 1

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

    offset = (/1,1,iset/)
    kount = (/lonp1,latp1,2/)

    call update_3d_fields( ncid_nh, offset, kount, amie_pot_nh,amie_ekv_nh,amie_efx_nh )

    if (iboxcar == 0) then
       pot_nh_amie(:,:) = (f1*amie_pot_nh(:,:,2) + &
            f2*amie_pot_nh(:,:,1))
       ekv_nh_amie(:,:) = (f1*amie_ekv_nh(:,:,2) + &
            f2*amie_ekv_nh(:,:,1))
       efx_nh_amie(:,:) = (f1*amie_efx_nh(:,:,2) + &
            f2*amie_efx_nh(:,:,1))
       !     write(iulog,"('ekv_nh_amie min, max = ',2e12.4)")
       !     |       minval(ekv_nh_amie),maxval(ekv_nh_amie)
    else
       call boxcar_ave(amie_pot_nh,pot_nh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
       call boxcar_ave(amie_efx_nh,efx_nh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
       call boxcar_ave(amie_ekv_nh,ekv_nh_amie,lonp1,latp1, &
            nn,iset,iboxcar)
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

    active_task: if ( mytid<ntask ) then

       !     mlongitude starts from 180 degree
       rot = sunlon*rtd
       if(rot < 0) then
          rot = rot + 360._r8    !  0 to 360 degrees
       end if
       rot = rot / 15._r8        !  convert from degree to hrs

       dmltm = 24._r8 / real(lonmx, kind=r8)
       do i = 1, lonp1
          xmlt = (real(i-1, kind=r8) * dmltm) - rot + 24._r8
          xmlt = MOD(xmlt, 24._r8)
          m = int(xmlt/dmltm + 1.01_r8)
          mp1 = m + 1
          if (mp1 > lonp1) mp1 = 2
          del = xmlt - (m-1)*dmltm
          !     Initialize arrays around equator
          do j = latp1+1, ithmx
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
          do j = 1, latp1
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

       !     ****     SET GRID SPACING DLATM, DLONM
       !     DMLAT=lat spacing in degrees of AMIE apex grid
       dmlat = 180._r8 / real(jmxm-1, kind=r8)
       dlatm = dmlat * dtr
       dlonm = 2._r8 * pi / real(lonmx, kind=r8)
       dmltm = 24._r8 / real(lonmx, kind=r8)
       !     ****
       !     ****     SET ARRAY YLATM (LATITUDE VALUES FOR GEOMAGNETIC GRID
       !     ****
       alatm(1) = -pi / 2._r8
       alat(1) = -90._r8
       alatm(jmxm) = pi / 2._r8
       alat(jmxm) = 90._r8
       do i = 2, ithmx
          alat(i) = alat(i-1)+dlatm/dtr
          alat(jmxm+1-i) = alat(jmxm+2-i)-dlatm/dtr
          alatm(i) = alatm(i-1)+dlatm
          alatm(jmxm+1-i) = alatm(jmxm+2-i)-dlatm
       end do
       alon(1) = -pi/dtr
       alonm(1) = -pi
       do i = 2, lonp1
          alon(i) = alon(i-1) + dlonm/dtr
          alonm(i) = alonm(i-1) + dlonm
       end do

       !     ylatm and ylonm are arrays of latitudes and longitudes of the
       !     distored magnetic grids in radian - from consdyn.h
       !     Convert from apex magnetic grid to distorted magnetic grid
       !
       !     Allocate workspace for regrid routine rgrd2.f:
       lw = nmlonp1+nmlat+2*nmlonp1
       if (.not. allocated(w)) then
          allocate(w(lw), stat=ier)
          call check_alloc(ier, 'getamie', 'w', lw=lw)
       end if
       liw = nmlonp1 + nmlat
       if (.not. allocated(iw)) then
          allocate(iw(liw), stat=ier)
          call check_alloc(ier, 'getamie', 'iw', lw=liw)
       end if
       intpol(:) = 1             ! linear (not cubic) interp in both dimensions
       if (alatm(1) > ylatm(1)) then
          alatm(1) = ylatm(1)
       end if
       if (alatm(jmxm) < ylatm(nmlat)) then
          alatm(jmxm) = ylatm(nmlat)
       end if
       if (alonm(1) > ylonm(1)) then
          alonm(1) = ylonm(1)
       end if
       if (alonm(lonp1) < ylonm(nmlonp1)) then
          alonm(lonp1) = ylonm(nmlonp1)
       end if
       !     write(iulog,"('  AMIE: ylatm =',/,(6e12.4))") ylatm
       !     write(iulog,"('  AMIE: ylonm =',/,(6e12.4))") ylonm
       !     write(iulog,"('  AMIE: potm(1,:) =',/,(6e12.4))") potm(1,:)
       !     ylatm from -pi/2 to pi/2, and ylonm from -pi to pi
       call rgrd2(lonp1, jmxm, alonm, alatm, potm, nmlonp1, nmlat,  &
            ylonm, ylatm, phihm, intpol, w, lw, iw, liw, ier)
       call rgrd2(lonp1, jmxm, alonm, alatm, ekvm, nmlonp1, nmlat,  &
            ylonm, ylatm, amie_kevm, intpol, w, lw, iw, liw, ier)
       call rgrd2(lonp1, jmxm, alonm, alatm, efxm, nmlonp1, nmlat,  &
            ylonm, ylatm, amie_efxm, intpol, w, lw, iw, liw, ier)

       if (iprint > 0 .and. masterproc) then
          write(iulog, *) subname, ': Max, min amie_efxm = ', &
               maxval(amie_efxm), minval(amie_efxm)
       end if
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
          write(iulog, "(a,': AMIE data interpolated to date and time')") subname
          write(iulog,"(a,': iyear,imo,iday,iutsec = ',3i6,i10)") subname,       &
               iyear, imo, iday, iutsec
          write(iulog,"(2a,i6,2F9.5,3I6,f10.4)")                                 &
               subname, ': AMIE iset f1,f2,year,mon,day,ut = ', iset,            &
               f1, f2, year(iset), month(iset), day(iset), amie_nh_ut(iset)
          write(iulog,*) subname, ': max,min phihm= ', maxval(phihm), minval(phihm)
          !     write(iulog,*)'getamie: max,min phihm,amie_efx,amie_kev = ',
          !     |    maxval(phihm),minval(tiepot),maxval(amie_efx),
          !     |    minval(amie_efx),maxval(amie_kev),minval(amie_kev)
       end if
    end if active_task
#else
    call endrun('Cannot use AMIE without electro-dynamo active.')
#endif
  end subroutine getamie
  
#ifdef WACCMX_EDYN_ESMF
  !-----------------------------------------------------------------------
  subroutine close_files

    deallocate( year,month,day )
    deallocate( amie_cusplat_nh, amie_cuspmlt_nh, amie_hpi_nh, &
                amie_pcp_nh, amie_nh_ut, &
                amie_cusplat_sh, amie_cuspmlt_sh, amie_hpi_sh, &
                amie_pcp_sh, amie_sh_ut )

    call cam_pio_closefile(ncid_nh)
    call cam_pio_closefile(ncid_sh)

    
  end subroutine close_files
  !-----------------------------------------------------------------------
  subroutine open_files()

    call rdamie_nh(amienh_files(file_ndx))
    call rdamie_sh(amiesh_files(file_ndx))

  end subroutine open_files
#endif

end module amie_module
