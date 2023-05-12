module amie_module
  !
  ! Module used to read data from the AMIE outputs (POT,mean energy,
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
  use utils_mod,      only: check_ncerr, check_alloc, boxcar_ave
  use edyn_mpi,       only: ntask, mytid
  use edyn_params,    only: pi, dtr, rtd
  use input_data_utils, only: time_coordinate

  implicit none

  private
  public :: init_amie
  public :: getamie

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
  ! cusplat_nh_input(sh) and cuspmlt_nh_input(sh) are
  !   AMIE cusp latitude and MLT in NH and SH
  ! hpi_nh(sh) are AMIE hemi-integrated power
  ! pcp_nh(sh) are AMIE polar-cap potential drop
  ! Time interpolated AMIE outputs with suffix _amie
  !
  real(r8), allocatable, dimension(:,:,:), save :: & ! (lonp1,latp1,ntimes)
       pot_nh_input, pot_sh_input,                   &
       ekv_nh_input, ekv_sh_input,                   &
       efx_nh_input, efx_sh_input
  real(r8), allocatable, dimension(:,:), save ::           &  ! (lonp1,latp1)
       pot_nh_amie, pot_sh_amie, ekv_nh_amie, ekv_sh_amie, &
       efx_nh_amie, efx_sh_amie
  integer,  allocatable, dimension(:), save ::                 & ! (ntimes)
       year, month, day, jday
  real(r8), allocatable, dimension(:), save ::                 & ! (ntimes)
       cusplat_nh_input, cuspmlt_nh_input, hpi_nh_input,          &
       pcp_nh_input, amie_nh_ut,                                &
       cusplat_sh_input, cuspmlt_sh_input, hpi_sh_input,          &
       pcp_sh_input, amie_sh_ut
  real(r8) ::                                                  &
       cusplat_nh_amie, cuspmlt_nh_amie, cusplat_sh_amie,      &
       cuspmlt_sh_amie, hpi_sh_amie, hpi_nh_amie, pcp_sh_amie, &
       pcp_nh_amie
  !
  type(file_desc_t) :: ncid_nh
  type(file_desc_t) :: ncid_sh

  character(len=cl), allocatable :: amienh_files(:)
  character(len=cl), allocatable :: amiesh_files(:)
  integer :: num_files, file_ndx

  type(time_coordinate) :: time_coord_nh
  type(time_coordinate) :: time_coord_sh

contains

  !-----------------------------------------------------------------------
  subroutine init_amie(amienh_list,amiesh_list)

    character(len=*),intent(in) :: amienh_list(:)
    character(len=*),intent(in) :: amiesh_list(:)

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

  end subroutine init_amie

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

    call time_coord_nh%initialize( amienh, set_weights=.false. )

    !
    ! Get time dimension:
    istat = pio_inquire(ncid_nh, unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid_nh, id_time, len=ntimes)
    call check_ncerr(istat, subname, 'AMIE time dimension')
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
    if (.not. allocated(cusplat_nh_input)) then
       allocate(cusplat_nh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'cusplat_nh_input', ntimes=ntimes)
    end if
    if (.not. allocated(cuspmlt_nh_input)) then
       allocate(cuspmlt_nh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'cuspmlt_nh_input', ntimes=ntimes)
    end if
    if (.not. allocated(hpi_nh_input)) then
       allocate(hpi_nh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'hpi_nh_input', ntimes=ntimes)
    end if
    if (.not. allocated(pcp_nh_input)) then
       allocate(pcp_nh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'pcp_nh_input', ntimes=ntimes)
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
    istat = pio_get_var(ncid_nh, idv_hpi, hpi_nh_input)
    call check_ncerr(istat, subname, 'AMIE hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid_nh, 'pcp', idv_pcp)
    call check_ncerr(istat, subname, 'AMIE pcp id')
    istat = pio_get_var(ncid_nh, idv_pcp, pcp_nh_input)
    call check_ncerr(istat, subname, 'AMIE pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid_nh, 'cusplat', idv_cusplat)
    call check_ncerr(istat, subname, 'AMIE cusplat id')
    istat = pio_get_var(ncid_nh, idv_cusplat, cusplat_nh_input)
    call check_ncerr(istat, subname, 'AMIE cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid_nh, 'cuspmlt', idv_cuspmlt)
    call check_ncerr(istat, subname, 'AMIE cuspmlt id')
    istat = pio_get_var(ncid_nh, idv_cuspmlt, cuspmlt_nh_input)
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
    if (.not. allocated(pot_nh_input)) then
       allocate(pot_nh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'pot_nh_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(ekv_nh_input)) then
       allocate(ekv_nh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'ekv_nh_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(efx_nh_input)) then
       allocate(efx_nh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'efx_nh_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
  end subroutine rdamie_nh

  !-----------------------------------------------------------------------
  subroutine rdamie_sh(amiesh)
    !
    ! Read AMIE data for the southern hemisphere from amiesh
    !

    ! Dummy argument
    character(len=*), intent(in) :: amiesh
    ! Local variables:
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

    call time_coord_sh%initialize( amiesh, set_weights=.false. )

    !
    ! Get time dimension:
    istat = pio_inquire(ncid_sh, unlimiteddimid=id_time)
    istat = pio_inquire_dimension(ncid_sh, id_time, len=ntimes)
    call check_ncerr(istat, subname, 'AMIE time dimension')
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
    if (.not. allocated(cusplat_sh_input)) then
       allocate(cusplat_sh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'cusplat_sh_input', ntimes=ntimes)
    end if
    if (.not. allocated(cuspmlt_sh_input)) then
       allocate(cuspmlt_sh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'cuspmlt_sh_input', ntimes=ntimes)
    end if
    if (.not. allocated(hpi_sh_input)) then
       allocate(hpi_sh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'hpi_sh_input', ntimes=ntimes)
    end if
    if (.not. allocated(pcp_sh_input)) then
       allocate(pcp_sh_input(ntimes), stat=ier)
       call check_alloc(ier, subname, 'pcp_sh_input', ntimes=ntimes)
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
    istat = pio_get_var(ncid_sh, idv_hpi, hpi_sh_input)
    call check_ncerr(istat, subname, 'AMIE hpi')
    !
    ! Get PCP
    istat = pio_inq_varid(ncid_sh, 'pcp', idv_pcp)
    call check_ncerr(istat, subname, 'AMIE pcp id')
    istat = pio_get_var(ncid_sh, idv_pcp, pcp_sh_input)
    call check_ncerr(istat, subname, 'AMIE pcp')
    !
    ! Get cusplat
    istat = pio_inq_varid(ncid_sh, 'cusplat', idv_cusplat)
    call check_ncerr(istat, subname, 'AMIE cusplat id')
    istat = pio_get_var(ncid_sh, idv_cusplat, cusplat_sh_input)
    call check_ncerr(istat, subname, 'AMIE cusplat')
    !
    ! Get cuspmlt
    istat = pio_inq_varid(ncid_sh, 'cuspmlt', idv_cuspmlt)
    call check_ncerr(istat, subname, 'AMIE cuspmlt id')
    istat = pio_get_var(ncid_sh, idv_cuspmlt, cuspmlt_sh_input)
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
    if (.not. allocated(pot_sh_input)) then
       allocate(pot_sh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'pot_sh_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(ekv_sh_input)) then
       allocate(ekv_sh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'ekv_sh_input', &
            lonp1=lonp1, latp1=latp1, ntimes=ntimes)
    end if
    if (.not. allocated(efx_sh_input)) then
       allocate(efx_sh_input(lonp1, latp1, 2), stat=ier)
       call check_alloc(ier, subname, 'efx_sh_input', &
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
    integer                     :: nn, iset1, iset2, m, mp1, n
    integer                     :: iboxcar
    real(r8)                    :: f1, f2
    real(r8)                    :: del, xmlt, dmlat, dlatm, dlonm, dmltm, rot
    integer                     :: offset(3), kount(3)
    character(len=*), parameter :: subname = 'getamie'

    phihm = fillvalue
    amie_efxm = fillvalue
    amie_kevm = fillvalue
    crad = fillvalue

    if (iprint > 0 .and. masterproc) then
       write(iulog,"(/,72('-'))")
       write(iulog,"(a,':')") subname
       write(iulog,"(a,i4,', iday = ',i3,', iutsec = ',i10)")                &
            'Initial requested iyear= ', iyear, iday, iutsec
    end if

    nn = size(amie_sh_ut)
    !
    !     Check times:
    !

    iamie = 1 - time_coord_sh%times_check()
    check_loop: do while( iamie/=1 )
       if (iamie==2) then
          if (masterproc) then
             write(iulog, "(a,': Model date prior to AMIE first date:',3I5)")   &
                  subname, year(1), month(1), day(1)
          end if
          return
       end if

       if (iamie==0) then
          if (masterproc) then
             write(iulog, "(a,': Model date beyond the AMIE last Data:',3I5)")  &
                  subname, year(nn), month(nn), day(nn)
          end if

          if (file_ndx<num_files) then
             file_ndx = file_ndx+1
             call close_files()
             call open_files()
             iamie = 1 - time_coord_sh%times_check()
          else
             return
          end if
       end if
    end do check_loop

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

    call time_coord_sh%advance()

    iset1 = time_coord_sh%indxs(1)
    iset2 = time_coord_sh%indxs(2)

    f1 = time_coord_sh%wghts(1)
    f2 = time_coord_sh%wghts(2)

    cusplat_sh_amie = (f1*cusplat_sh_input(iset1) + &
                       f2*cusplat_sh_input(iset2))
    cuspmlt_sh_amie = (f1*cuspmlt_sh_input(iset1) + &
                       f2*cuspmlt_sh_input(iset2))
    hpi_sh_amie = (f1*hpi_sh_input(iset1) + f2*hpi_sh_input(iset2))
    pcp_sh_amie = (f1*pcp_sh_input(iset1) + f2*pcp_sh_input(iset2))

    offset = (/1,1,iset1/)
    kount = (/lonp1,latp1,2/)

    call update_3d_fields( ncid_sh, offset, kount, pot_sh_input,ekv_sh_input,efx_sh_input )
    if (iboxcar == 0) then
       pot_sh_amie(:,:) = (f1*pot_sh_input(:,:,1) + &
                           f2*pot_sh_input(:,:,2))
       ekv_sh_amie(:,:) = (f1*ekv_sh_input(:,:,1) + &
                           f2*ekv_sh_input(:,:,2))
       efx_sh_amie(:,:) = (f1*efx_sh_input(:,:,1) + &
                           f2*efx_sh_input(:,:,2))
    else
       call boxcar_ave(pot_sh_input,pot_sh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
       call boxcar_ave(efx_sh_input,efx_sh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
       call boxcar_ave(ekv_sh_input,ekv_sh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
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

    call time_coord_nh%advance()

    iset1 = time_coord_nh%indxs(1)
    iset2 = time_coord_nh%indxs(2)

    f1 = time_coord_nh%wghts(1)
    f2 = time_coord_nh%wghts(2)

    cusplat_nh_amie = (f1*cusplat_nh_input(iset1) + &
                       f2*cusplat_nh_input(iset2))
    cuspmlt_nh_amie = (f1*cuspmlt_nh_input(iset1) + &
                       f2*cuspmlt_nh_input(iset2))
    hpi_nh_amie = (f1*hpi_nh_input(iset1) + f2*hpi_nh_input(iset2))
    pcp_nh_amie = (f1*pcp_nh_input(iset1) + f2*pcp_nh_input(iset2))

    offset = (/1,1,iset1/)
    kount = (/lonp1,latp1,2/)

    call update_3d_fields( ncid_nh, offset, kount, pot_nh_input,ekv_nh_input,efx_nh_input )

    if (iboxcar == 0) then
       pot_nh_amie(:,:) = (f1*pot_nh_input(:,:,1) + &
                           f2*pot_nh_input(:,:,2))
       ekv_nh_amie(:,:) = (f1*ekv_nh_input(:,:,1) + &
                           f2*ekv_nh_input(:,:,2))
       efx_nh_amie(:,:) = (f1*efx_nh_input(:,:,1) + &
                           f2*efx_nh_input(:,:,2))
    else
       call boxcar_ave(pot_nh_input,pot_nh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
       call boxcar_ave(efx_nh_input,efx_nh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
       call boxcar_ave(ekv_nh_input,ekv_nh_amie,lonp1,latp1, &
            nn,iset1,iboxcar)
    end if
    !
    !     The OLTMAX latitude also defines the co-latitude theta0, which in
    !     turn determines crit1(+2.5deg) and crit2(-12.5deg) which are used
    !     in TIE-GCM as the boundaries of the polar cap and the region of
    !     influence of the high-lat potential versus the low-lat dynamo potential
    !     Define this latitude to be between 70 and 77.5 degrees
    !
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
          alat(i) = alat(i-1)+dlatm*rtd
          alat(jmxm+1-i) = alat(jmxm+2-i)-dlatm*rtd
          alatm(i) = alatm(i-1)+dlatm
          alatm(jmxm+1-i) = alatm(jmxm+2-i)-dlatm
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
          write(iulog, "(a,': AMIE data interpolated to date and time')") subname
          write(iulog,"(a,': iyear,imo,iday,iutsec = ',3i6,i10)") subname,       &
               iyear, imo, iday, iutsec
          write(iulog,"(2a,i6,2F9.5,3I6,f10.4)")                                 &
               subname, ': AMIE iset1 f1,f2,year,mon,day,ut = ', iset1,          &
               f1, f2, year(iset1), month(iset1), day(iset1), amie_nh_ut(iset1)
          write(iulog,*) subname, ': max,min phihm= ', maxval(phihm), minval(phihm)
       end if
    end if active_task

  end subroutine getamie

  !-----------------------------------------------------------------------
  subroutine close_files

    deallocate( year,month,day )
    deallocate( cusplat_nh_input, cuspmlt_nh_input, hpi_nh_input, &
                pcp_nh_input, amie_nh_ut, &
                cusplat_sh_input, cuspmlt_sh_input, hpi_sh_input, &
                pcp_sh_input, amie_sh_ut )

    call cam_pio_closefile(ncid_nh)
    call cam_pio_closefile(ncid_sh)


  end subroutine close_files
  !-----------------------------------------------------------------------
  subroutine open_files()

    call rdamie_nh(amienh_files(file_ndx))
    call rdamie_sh(amiesh_files(file_ndx))

  end subroutine open_files

end module amie_module
