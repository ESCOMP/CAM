module mo_drydep

  !---------------------------------------------------------------------
  !       ... Dry deposition
  !---------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use chem_mods,        only : gas_pcnst
  use pmgrid,           only : plev
  use spmd_utils,       only : masterproc
  use ppgrid,           only : pcols, begchunk, endchunk
  use mo_tracname,      only : solsym
  use cam_abortutils,   only : endrun
  use ioFileMod,        only : getfil
  use pio
  use cam_pio_utils,    only : cam_pio_openfile, cam_pio_closefile
  use cam_logfile,      only : iulog
  use dyn_grid,         only : get_dyn_grid_parm, get_horiz_grid_d
  use scamMod,          only : single_column

  use shr_drydep_mod,   only : nddvels =>  n_drydep, drydep_list, mapping
  use physconst,        only : karman

  use infnan,                only : nan, assignment(=)

  implicit none

  save

  interface drydep_inti
     module procedure dvel_inti_xactive
  end interface

  interface drydep
     module procedure drydep_fromlnd
  end interface

  private

  public :: drydep_inti, drydep, has_drydep
  public :: drydep_update
  public :: n_land_type, fraction_landuse, drydep_srf_file

  integer :: pan_ndx, mpan_ndx, o3_ndx, ch4_ndx, co_ndx, h2_ndx, ch3cooh_ndx
  integer :: sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx

  integer :: so2_ndx, ch3cn_ndx, hcn_ndx, hcooh_ndx

  integer :: o3a_ndx,xpan_ndx,xmpan_ndx

  integer :: cohc_ndx=-1, come_ndx=-1
  integer, parameter :: NTAGS = 50
  integer :: cotag_ndx(NTAGS)
  integer :: tag_cnt

  real(r8), parameter    :: small_value = 1.e-36_r8
  real(r8), parameter    :: large_value = 1.e36_r8
  real(r8), parameter    :: diffm       = 1.789e-5_r8
  real(r8), parameter    :: diffk       = 1.461e-5_r8
  real(r8), parameter    :: difft       = 2.060e-5_r8
  real(r8), parameter    :: vonkar      = karman
  real(r8), parameter    :: ric         = 0.2_r8
  real(r8), parameter    :: r           = 287.04_r8
  real(r8), parameter    :: cp          = 1004._r8
  real(r8), parameter    :: grav        = 9.81_r8
  real(r8), parameter    :: p00         = 100000._r8
  real(r8), parameter    :: wh2o        = 18.0153_r8
  real(r8), parameter    :: ph          = 1.e-5_r8
  real(r8), parameter    :: ph_inv      = 1._r8/ph
  real(r8), parameter    :: rovcp = r/cp

  logical, public :: has_dvel(gas_pcnst) = .false.
  integer         :: map_dvel(gas_pcnst) = 0

  real(r8), protected, allocatable  :: fraction_landuse(:,:,:)
  real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
  real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer
  integer, parameter :: n_land_type = 11

  integer, allocatable :: spc_ndx(:) ! nddvels
  real(r8), public :: crb

  type lnd_dvel_type
     real(r8), pointer :: dvel(:,:)   ! deposition velocity over land (cm/s)
  end type lnd_dvel_type

  type(lnd_dvel_type), allocatable :: lnd(:)
  character(len=SHR_KIND_CL) :: drydep_srf_file

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  subroutine dvel_inti_fromlnd
    use mo_chem_utls,         only : get_spc_ndx
    use cam_abortutils,       only : endrun

    integer :: ispc

    allocate(spc_ndx(nddvels))
    allocate( lnd(begchunk:endchunk) )

    do ispc = 1,nddvels

       spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
       if (spc_ndx(ispc) < 1) then
          write(*,*) 'drydep_inti: '//trim(drydep_list(ispc))//' is not included in species set'
          call endrun('drydep_init: invalid dry deposition species')
       endif

    enddo

    crb = (difft/diffm)**(2._r8/3._r8) !.666666_r8

  endsubroutine dvel_inti_fromlnd

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_update( state, cam_in )
    use physics_types,   only : physics_state
    use camsrfexch,      only : cam_in_t

    type(physics_state), intent(in) :: state           ! Physics state variables
    type(cam_in_t),  intent(in) :: cam_in

    if (nddvels<1) return

    lnd(state%lchnk)%dvel => cam_in%depvel

  end subroutine drydep_update

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_fromlnd( ocnfrac, icefrac, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvelocity, dflx, mmr, &
                             tv, ncol, lchnk )

    !-------------------------------------------------------------------------------------
    ! combines the deposition velocities provided by the land model with deposition
    ! velocities over ocean and sea ice
    !-------------------------------------------------------------------------------------

    use ppgrid,         only : pcols
    use chem_mods,      only : gas_pcnst

#if (defined OFFLINE_DYN)
    use metdata, only: get_met_fields
#endif

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------

    real(r8), intent(in)      :: icefrac(pcols)
    real(r8), intent(in)      :: ocnfrac(pcols)
    integer,  intent(in)      :: ncol
    integer,  intent(in)      :: lchnk                    ! chunk number
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvelocity(ncol,gas_pcnst)    ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: ocnice_dvel(ncol,gas_pcnst)
    real(r8) :: ocnice_dflx(pcols,gas_pcnst)

    real(r8), dimension(ncol) :: term    ! work array
    integer  :: ispec
    real(r8)  :: lndfrac(pcols)
#if (defined OFFLINE_DYN)
    real(r8)  :: met_ocnfrac(pcols)
    real(r8)  :: met_icefrac(pcols)
#endif
    integer :: i

    lndfrac(:ncol) = 1._r8 - ocnfrac(:ncol) - icefrac(:ncol)

    where( lndfrac(:ncol) < 0._r8 )
       lndfrac(:ncol) = 0._r8
    endwhere

#if (defined OFFLINE_DYN)
    call get_met_fields(lndfrac, met_ocnfrac, met_icefrac, lchnk, ncol)
#endif

    !-------------------------------------------------------------------------------------
    !   ... initialize
    !-------------------------------------------------------------------------------------
    dvelocity(:,:) = 0._r8

    !-------------------------------------------------------------------------------------
    !   ... compute the dep velocities over ocean and sea ice
    !       land type 7 is used for ocean
    !       land type 8 is used for sea ice
    !-------------------------------------------------------------------------------------
    call drydep_xactive( sfc_temp, pressure_sfc,  &
                         wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                         snow, solar_flux, ocnice_dvel, ocnice_dflx, mmr, &
                         tv, ncol, lchnk, &
#if (defined OFFLINE_DYN)
                         ocnfrc=met_ocnfrac,icefrc=met_icefrac, beglandtype=7, endlandtype=8 )
#else
                         ocnfrc=ocnfrac,icefrc=icefrac, beglandtype=7, endlandtype=8 )
#endif
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... merge the land component with the non-land component
       !            ocn and ice already have fractions factored in
       !-------------------------------------------------------------------------------------
       dvelocity(:ncol,spc_ndx(ispec)) = lnd(lchnk)%dvel(:ncol,ispec)*lndfrac(:ncol) &
                                       + ocnice_dvel(:ncol,spc_ndx(ispec))
    enddo

    !-------------------------------------------------------------------------------------
    !        ... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx>0 ) then
       dvelocity(:ncol,mpan_ndx) = dvelocity(:ncol,mpan_ndx)/3._r8
    endif
    if( xmpan_ndx>0 ) then
       dvelocity(:ncol,xmpan_ndx) = dvelocity(:ncol,xmpan_ndx)/3._r8
    endif
    if( hcn_ndx>0 ) then
       dvelocity(:ncol,hcn_ndx) = ocnice_dvel(:ncol,hcn_ndx) ! should be zero over land
    endif
    if( ch3cn_ndx>0 ) then
       dvelocity(:ncol,ch3cn_ndx) = ocnice_dvel(:ncol,ch3cn_ndx) ! should be zero over land
    endif

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0 .and. ch3cooh_ndx > 0 ) then
       if( has_dvel(hcooh_ndx) ) then
          dvelocity(:ncol,hcooh_ndx) = dvelocity(:ncol,ch3cooh_ndx)
       end if
    end if

    !-------------------------------------------------------------------------------------
    !        ... assign CO tags to CO
    ! put this kludge in for now ...
    !  -- should be able to set all these via the table mapping in shr_drydep_mod
    !-------------------------------------------------------------------------------------
    if( cohc_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cohc_ndx)
    endif
    if( come_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,come_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,come_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,come_ndx)
    endif

    if ( co_ndx>0 ) then
       do i=1,tag_cnt
          dvelocity(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx)
          dflx(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cotag_ndx(i))
       enddo
    endif

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... compute the deposition flux
       !-------------------------------------------------------------------------------------
       dflx(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec)) * term(:ncol) * mmr(:ncol,plev,spc_ndx(ispec))
    end do

  end subroutine drydep_fromlnd

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine dvel_inti_xactive( depvel_lnd_file )
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------
    use dycore,        only : dycore_is
    use mo_chem_utls,  only : get_spc_ndx
    use phys_control,  only : phys_getopts

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    character(len=*), intent(in) :: depvel_lnd_file

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i
    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: dimid
    integer :: m
    integer :: astat
    integer :: plon, plat
    integer :: ierr, ndx

    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg_edge(:)

    character(len=32) :: test_name
    character(len=4) :: tag_name
    type(file_desc_t) :: piofile
    type(var_desc_t) :: vid

    character(len=shr_kind_cl) :: locfn
    logical :: prog_modal_aero

    ! determine if modal aerosols are active so that fraction_landuse array is initialized for modal aerosal dry dep
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    call dvel_inti_fromlnd()

    if( masterproc ) then
       write(iulog,*) 'drydep_inti: following species have dry deposition'
       do i=1,nddvels
          if( len_trim(drydep_list(i)) > 0 ) then
             write(iulog,*) 'drydep_inti: '//trim(drydep_list(i))//' is requested to have dry dep'
          endif
       enddo
       write(iulog,*) 'drydep_inti:'
    endif

    !-------------------------------------------------------------------------------------
    ! 	... get species indices
    !-------------------------------------------------------------------------------------
    xpan_ndx      = get_spc_ndx( 'XPAN' )
    xmpan_ndx     = get_spc_ndx( 'XMPAN' )
    o3a_ndx       = get_spc_ndx( 'O3A' )

    ch4_ndx      = get_spc_ndx( 'CH4' )
    h2_ndx       = get_spc_ndx( 'H2' )
    co_ndx       = get_spc_ndx( 'CO' )
    pan_ndx      = get_spc_ndx( 'PAN' )
    mpan_ndx     = get_spc_ndx( 'MPAN' )
    o3_ndx       = get_spc_ndx( 'OX' )
    if( o3_ndx < 0 ) then
       o3_ndx  = get_spc_ndx( 'O3' )
    end if
    so2_ndx     = get_spc_ndx( 'SO2' )
    ch3cooh_ndx = get_spc_ndx( 'CH3COOH')

    sogm_ndx   = get_spc_ndx( 'SOGM' )
    sogi_ndx   = get_spc_ndx( 'SOGI' )
    sogt_ndx   = get_spc_ndx( 'SOGT' )
    sogb_ndx   = get_spc_ndx( 'SOGB' )
    sogx_ndx   = get_spc_ndx( 'SOGX' )

    hcn_ndx     = get_spc_ndx( 'HCN')
    ch3cn_ndx   = get_spc_ndx( 'CH3CN')

    cohc_ndx     = get_spc_ndx( 'COhc' )
    come_ndx     = get_spc_ndx( 'COme' )

    tag_cnt=0
    cotag_ndx(:)=-1
    do i = 1,NTAGS
       write(tag_name,'(a2,i2.2)') 'CO',i
       ndx = get_spc_ndx(tag_name)
       if (ndx>0) then
          tag_cnt = tag_cnt+1
          cotag_ndx(tag_cnt) = ndx
       endif
    enddo

    do i=1,nddvels
       if ( mapping(i) > 0 ) then
          test_name = drydep_list(i)
          m = get_spc_ndx( test_name )
          has_dvel(m) = .true.
          map_dvel(m) = i
       endif
    enddo

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !---------------------------------------------------------------------------
    ! 	... allocate module variables
    !---------------------------------------------------------------------------
    allocate( dep_ra(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
       call endrun('dvel_inti: failed to allocate dep_ra')
    end if
    allocate( dep_rb(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
       call endrun('dvel_inti: failed to allocate dep_rb')
    end if

    if (.not.prog_modal_aero) then
       return
    endif

    allocate( fraction_landuse(pcols,n_land_type, begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun('dvel_inti: failed to allocate fraction_landuse')
    end if
    fraction_landuse = nan

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

    if(dycore_is('UNSTRUCTURED') ) then
       call get_landuse_and_soilw_from_file()
    else
       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------
       call getfil (depvel_lnd_file, locfn, 0)
       call cam_pio_openfile (piofile, trim(locfn), PIO_NOWRITE)
       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       ierr = pio_inq_dimid( piofile, 'lon', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlon_veg )
       ierr = pio_inq_dimid( piofile, 'lat', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, nlat_veg )
       ierr = pio_inq_dimid( piofile, 'pft', dimid )
       ierr = pio_inq_dimlen( piofile, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation_map; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation_map')
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation_map; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation_map')
       end if
       allocate( lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegetation lon, lat arrays; error = ',astat
          call endrun('dvel_inti: failed to allocate vegetation lon, lat arrays')
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       ierr = pio_inq_varid( piofile, 'PCT_PFT', vid )
       ierr = pio_get_var( piofile, vid, vegetation_map )

       ierr = pio_inq_varid( piofile, 'LANDMASK', vid )
       ierr = pio_get_var( piofile, vid, landmask )

       ierr = pio_inq_varid( piofile, 'PCT_URBAN', vid )
       ierr = pio_get_var( piofile, vid, urban )

       ierr = pio_inq_varid( piofile, 'PCT_LAKE', vid )
       ierr = pio_get_var( piofile, vid, lake )

       ierr = pio_inq_varid( piofile, 'PCT_WETLAND', vid )
       ierr = pio_get_var( piofile, vid, wetland )

       call cam_pio_closefile( piofile )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
       urban(:,:)            = .01_r8 * urban(:,:)
#ifdef DEBUG
       if(masterproc) then
          write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
          write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
          write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
       end if
#endif
       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------
       call interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
            lon_veg_edge, landmask, urban, lake, &
            wetland, vegetation_map )

       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
    endif  ! Unstructured grid

  end subroutine dvel_inti_xactive

  !-------------------------------------------------------------------------------------
  subroutine get_landuse_and_soilw_from_file()
    use ncdio_atm, only : infld

    logical :: readvar

    type(file_desc_t) :: piofile
    character(len=shr_kind_cl) :: locfn
    logical :: lexist

    if (len_trim(drydep_srf_file) == 0) then
       write(iulog,*)'**************************************'
       write(iulog,*)' get_landuse_and_soilw_from_file: INFO:'
       write(iulog,*)' drydep_srf_file not set:'
       write(iulog,*)' setting fraction_landuse to zero'
       write(iulog,*)'**************************************'
       fraction_landuse = 0._r8
       return
    end if

    call getfil (drydep_srf_file, locfn, 1, lexist)
    if(lexist) then
       call cam_pio_openfile(piofile, locfn, PIO_NOWRITE)

       call infld('fraction_landuse', piofile, 'ncol','class',1,pcols,1,n_land_type, begchunk,endchunk, &
            fraction_landuse, readvar, gridname='physgrid')
       if (.not. readvar) then
          write(iulog,*)'**************************************'
          write(iulog,*)'get_landuse_and_soilw_from_file: INFO:'
          write(iulog,*)' fraction_landuse not read from file: '
          write(iulog,*)' ', trim(locfn)
          write(iulog,*)' setting all values to zero'
          write(iulog,*)'**************************************'
          fraction_landuse = 0._r8
       end if

       call cam_pio_closefile(piofile)
    else
       call endrun('Unstructured grids require drydep_srf_file ')
    end if


  end subroutine get_landuse_and_soilw_from_file

  !-------------------------------------------------------------------------------------
  subroutine interp_map( plon, plat, nlon_veg, nlat_veg, npft_veg, lat_veg_edge, &
                         lon_veg_edge, landmask, urban, lake, &
                         wetland, vegetation_map )

    use mo_constants, only : r2d
    use scamMod, only : latiop,loniop,scmlat,scmlon,scm_cambfb_mode
    use shr_scam_mod  , only: shr_scam_getCloseLatLon  ! Standardized system subroutines
    use cam_initfiles, only: initial_file_get_id
    use dycore, only : dycore_is
    use phys_grid,     only : get_rlat_all_p, get_rlon_all_p, get_ncols_p

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer,  intent(in)         :: plon, plat, nlon_veg, nlat_veg, npft_veg
    real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
    real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
    real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8) :: closelat,closelon
    integer :: latidx,lonidx

    integer, parameter           :: veg_ext = 20
    type(file_desc_t), pointer   :: piofile
    integer                      :: i, j, ii, jj, i_ndx, n
    integer, dimension(plon+1)   :: ind_lon
    integer, dimension(plat+1)  :: ind_lat
    real(r8)                         :: total_land
    real(r8), dimension(plon+1)      :: lon_edge
    real(r8), dimension(plat+1)     :: lat_edge
    real(r8)                         :: lat1, lon1
    real(r8)                         :: x1, x2, y1, y2, dx, dy
    real(r8)                         :: area, total_area
    real(r8), dimension(npft_veg+3)  :: fraction
    real(r8), dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
    integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext

    real(r8), allocatable :: lam(:), phi(:)

    logical, parameter :: has_npole = .true.
    integer :: ploniop,platiop
    real(r8) :: tmp_frac_lu(plon,n_land_type,plat)

    real(r8):: rlats(pcols), rlons(pcols)
    integer :: lchnk, ncol, icol
    logical :: found

    if(dycore_is('UNSTRUCTURED') ) then
       call endrun('mo_drydep::interp_map called for UNSTRUCTURED grid')
    endif

    allocate(lam(plon), phi(plat))
    call get_horiz_grid_d(plat, clat_d_out=phi)
    call get_horiz_grid_d(plon, clon_d_out=lam)

    if (single_column) then
       if (scm_cambfb_mode) then
          piofile => initial_file_get_id()
          call shr_scam_getCloseLatLon(piofile,scmlat,scmlon,closelat,closelon,latidx,lonidx)
          ploniop=size(loniop)
          platiop=size(latiop)
       else
          latidx=1
          lonidx=1
          ploniop=1
          platiop=1
       end if

       lon_edge(1) = loniop(lonidx) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d

       if (lonidx.lt.ploniop) then
          lon_edge(2) = loniop(lonidx+1) * r2d - .5_r8*(loniop(2) - loniop(1)) * r2d
       else
          lon_edge(2) = lon_edge(1) + (loniop(2) - loniop(1)) * r2d
       end if

       lat_edge(1) = latiop(latidx) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d

       if (latidx.lt.platiop) then
          lat_edge(2) = latiop(latidx+1) * r2d - .5_r8*(latiop(2) - latiop(1)) * r2d
       else
          lat_edge(2) = lat_edge(1) + (latiop(2) - latiop(1)) * r2d
       end if
    else
       do i = 1,plon
          lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
       end do
       lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
       if( .not. has_npole ) then
          do j = 1,plat+1
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
       else
          do j = 1,plat
             lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
          end do
          lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
       end if
    end if
    do j = 1,plat+1
       lat_edge(j) = min( lat_edge(j), 90._r8 )
       lat_edge(j) = max( lat_edge(j),-90._r8 )
    end do

    !-------------------------------------------------------------------------------------
    ! wrap around the longitudes
    !-------------------------------------------------------------------------------------
    do i = -veg_ext,0
       lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
       mapping_ext     (i) =              nlon_veg+i
    end do
    do i = 1,nlon_veg
       lon_veg_edge_ext(i) = lon_veg_edge(i)
       mapping_ext     (i) =              i
    end do
    do i = nlon_veg+1,nlon_veg+veg_ext
       lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
       mapping_ext     (i) =              i-nlon_veg
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : lon_edge ',lon_edge
    write(iulog,*) 'interp_map : lat_edge ',lat_edge
    write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
#endif
    do j = 1,plon+1
       lon1 = lon_edge(j)
       do i = -veg_ext,nlon_veg+veg_ext
          dx = lon_veg_edge_ext(i  ) - lon1
          dy = lon_veg_edge_ext(i+1) - lon1
          if( dx*dy <= 0._r8 ) then
             ind_lon(j) = i
             exit
          end if
       end do
    end do

    do j = 1,plat+1
       lat1 = lat_edge(j)
       do i = 1,nlat_veg
          dx = lat_veg_edge(i  ) - lat1
          dy = lat_veg_edge(i+1) - lat1
          if( dx*dy <= 0._r8 ) then
             ind_lat(j) = i
             exit
          end if
       end do
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : ind_lon ',ind_lon
    write(iulog,*) 'interp_map : ind_lat ',ind_lat
#endif
    lat_loop : do j = 1,plat
       lon_loop : do i = 1,plon
          total_area       = 0._r8
          fraction         = 0._r8
          do jj = ind_lat(j),ind_lat(j+1)
             y1 = max( lat_edge(j),lat_veg_edge(jj) )
             y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) )
             dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
             do ii =ind_lon(i),ind_lon(i+1)
                i_ndx = mapping_ext(ii)
                x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) )
                dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                area = dx * dy
                total_area = total_area + area
                !-----------------------------------------------------------------
                ! 	... special case for ocean grid point
                !-----------------------------------------------------------------
                if( nint(landmask(i_ndx,jj)) == 0 ) then
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area
                else
                   do n = 1,npft_veg
                      fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                   end do
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                   fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                   fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                   !-----------------------------------------------------------------
                   ! 	... check if land accounts for the whole area.
                   !           If not, the remaining area is in the ocean
                   !-----------------------------------------------------------------
                   total_land = sum(vegetation_map(i_ndx,jj,:)) &
                              + urban  (i_ndx,jj) &
                              + lake   (i_ndx,jj) &
                              + wetland(i_ndx,jj)
                   if( total_land < 1._r8 ) then
                      fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                   end if
                end if
             end do
          end do
          !-------------------------------------------------------------------------------------
          ! 	... divide by total area of grid box
          !-------------------------------------------------------------------------------------
          fraction(:) = fraction(:)/total_area
          !-------------------------------------------------------------------------------------
          ! 	... make sure we don't have too much or too little
          !-------------------------------------------------------------------------------------
          if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
             fraction(:) = fraction(:)/sum(fraction)
          end if
          !-------------------------------------------------------------------------------------
          ! 	... map to Wesely land classification
          !-------------------------------------------------------------------------------------
          tmp_frac_lu(i, 1, j) =     fraction(20)
          tmp_frac_lu(i, 2, j) = sum(fraction(16:17))
          tmp_frac_lu(i, 3, j) = sum(fraction(13:15))
          tmp_frac_lu(i, 4, j) = sum(fraction( 5: 9))
          tmp_frac_lu(i, 5, j) = sum(fraction( 2: 4))
          tmp_frac_lu(i, 6, j) =     fraction(19)
          tmp_frac_lu(i, 7, j) =     fraction(18)
          tmp_frac_lu(i, 8, j) =     fraction( 1)
          tmp_frac_lu(i, 9, j) = 0._r8
          tmp_frac_lu(i,10, j) = 0._r8
          tmp_frac_lu(i,11, j) = sum(fraction(10:12))
       end do lon_loop
    end do lat_loop

    do lchnk = begchunk, endchunk
       ncol = get_ncols_p(lchnk)
       call get_rlat_all_p(lchnk, ncol, rlats(:ncol))
       call get_rlon_all_p(lchnk, ncol, rlons(:ncol))
       do icol= 1,ncol
          found=.false.
          find_col: do j = 1,plat
             do i = 1,plon
                if (rlats(icol)==phi(j) .and. rlons(icol)==lam(i)) then
                   found=.true.
                   exit find_col
                endif
             enddo
          enddo find_col

          if (.not.found) call endrun('mo_drydep::interp_map not able find physics column coordinate')
          fraction_landuse(icol,1:n_land_type,lchnk) =  tmp_frac_lu(i,1:n_land_type,j)

       end do

       !-------------------------------------------------------------------------------------
       ! 	... make sure there are no out of range values
       !-------------------------------------------------------------------------------------
       where (fraction_landuse(:ncol,:n_land_type,lchnk) < 0._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 0._r8
       where (fraction_landuse(:ncol,:n_land_type,lchnk) > 1._r8) fraction_landuse(:ncol,:n_land_type,lchnk) = 1._r8
    end do

  end subroutine interp_map

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_xactive( sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvel, dflx, mmr, &
                             tv, ncol, lchnk, &
                             ocnfrc, icefrc, beglandtype, endlandtype )
    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

    use shr_drydep_mod, only: z0, rgso, rgss, ri, rclo, rcls, rlu, rac
    use shr_drydep_mod, only: shr_drydep_setHCoeff, foxd, drat
    use physconst,      only: tmelt

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)

    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
    real(r8), intent(in)      :: mmr(pcols,plev,gas_pcnst)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvel(ncol,gas_pcnst)        ! deposition velocity (cm/s)
    real(r8), intent(inout)   :: dflx(pcols,gas_pcnst)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   lchnk                   ! chunk number

    integer, intent(in), optional     ::  beglandtype
    integer, intent(in), optional     ::  endlandtype

    real(r8), intent(in), optional      :: ocnfrc(pcols)
    real(r8), intent(in), optional      :: icefrc(pcols)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8), parameter :: scaling_to_cm_per_s = 100._r8
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx

    real(r8) :: slope = 0._r8
    real(r8) :: z0water ! revised z0 over water
    real(r8) :: p       ! pressure at midpoint first layer
    real(r8) :: pg      ! surface pressure
    real(r8) :: es      ! saturation vapor pressure
    real(r8) :: ws      ! saturation mixing ratio
    real(r8) :: hvar    ! constant to compute xmol
    real(r8) :: h       ! constant to compute xmol
    real(r8) :: psih    ! stability correction factor
    real(r8) :: rs      ! constant for calculating rsmx
    real(r8) :: rmx     ! resistance by vegetation
    real(r8) :: zovl    ! ratio of z to  m-o length
    real(r8) :: cvarb   ! cvar averaged over landtypes
    real(r8) :: bb      ! b averaged over landtypes
    real(r8) :: ustarb  ! ustar averaged over landtypes
    real(r8) :: tc(ncol)  ! temperature in celsius
    real(r8) :: cts(ncol) ! correction to rlu rcl and rgs for frost

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,nddvels) :: heff

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(ncol,n_land_type)
    real(r8), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
    real(r8), dimension(ncol) :: thg     ! ground virtual potential temperature
    real(r8), dimension(ncol) :: z       ! height of lowest level
    real(r8), dimension(ncol) :: va      ! magnitude of v on cross points
    real(r8), dimension(ncol) :: ribn    ! richardson number
    real(r8), dimension(ncol) :: qs      ! saturation specific humidity
    real(r8), dimension(ncol) :: crs     ! multiplier to calculate crs
    real(r8), dimension(ncol) :: rdc     ! part of lower canopy resistance
    real(r8), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
    real(r8), dimension(ncol) :: z0b     ! average roughness length over grid
    real(r8), dimension(ncol) :: wrk     ! work array
    real(r8), dimension(ncol) :: term    ! work array
    real(r8), dimension(ncol) :: resc    ! work array
    real(r8), dimension(ncol) :: lnd_frc ! work array
    logical,  dimension(ncol) :: unstable
    logical,  dimension(ncol) :: has_rain
    logical,  dimension(ncol) :: has_dew

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type) :: rds   ! resistance for deposition of sulfate
    real(r8), dimension(ncol,n_land_type) :: b     ! buoyancy parameter for unstable conditions
    real(r8), dimension(ncol,n_land_type) :: cvar  ! height parameter
    real(r8), dimension(ncol,n_land_type) :: ustar ! friction velocity
    real(r8), dimension(ncol,n_land_type) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rclx  ! lower canopy resistance
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
    real(r8), dimension(ncol,n_land_type,gas_pcnst) :: rgsx  ! ground resistance
    real(r8) :: vds
    logical  :: fr_lnduse(ncol,n_land_type)           ! wrking array
    real(r8) :: dewm                                  ! multiplier for rs when dew occurs

    real(r8) :: lcl_frc_landuse(ncol,n_land_type)

    integer :: beglt, endlt

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(r8) :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)

    if (present( beglandtype)) then
      beglt = beglandtype
    else
      beglt = 1
    endif
    if (present( endlandtype)) then
      endlt = endlandtype
    else
      endlt = n_land_type
    endif

    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,gas_pcnst
       dvel(:,m) = 0._r8
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
    call shr_drydep_setHCoeff( ncol, sfc_temp, heff )

    do lt = 1,n_land_type
       dep_ra (:,lt,lchnk)   = 0._r8
       dep_rb (:,lt,lchnk)   = 0._r8
       rds(:,lt)   = 0._r8
    end do

    !-------------------------------------------------------------------------------------
    ! season index only for ocn and sea ice
    !-------------------------------------------------------------------------------------
    index_season = 4
    !-------------------------------------------------------------------------------------
    ! special case for snow covered terrain
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       if( snow(i) > .01_r8 ) then
          index_season(i,:) = 4
       end if
    end do
    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_r8,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._r8
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_r8*es/(pg - es)
       qs(i) = ws/(1._r8 + ws)
       has_dew(i) = .false.
       if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
       end if
       if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
       end if
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       tc(i) = sfc_temp(i) - tmelt
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc(i)*(40._r8 - tc(i))))
       else
          crs(i) = large_value
       end if
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._r8*(1._r8 + 1000._r8/(solar_flux(i) + 10._r8))/(1._r8 + 1000._r8*slope)
    end do col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working arrays
    !-------------------------------------------------------------------------------------
    lcl_frc_landuse(:,:) = 0._r8

    if ( present(ocnfrc) .and. present(icefrc) ) then
       do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,7) = ocnfrc(i)
          lcl_frc_landuse(i,8) = icefrc(i)
       enddo
    endif
    do lt = 1,n_land_type
       do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._r8
       enddo
    end do

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._r8
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
       else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
       end if
       uustar(i) = va(i)*ustarb
    end do

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( unstable(i) ) then
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*b(i,lt))) ) )
             else
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
             end if
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = 7
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
             b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8+ 7.4_r8*b(i,lt))) ) )
          else
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log(z(i)/z0water))
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
          end if
       end if
    end do

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then                      ! unstable
                h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
             else
                h = hvar/((1._r8+4.7_r8*ribn(i))**2)
             end if
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._r8 ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._r8,zovl )
                psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
                vds  = 2.e-3_r8*ustar(i,lt) * (1._r8 + (300/(-xmol(i,lt)))**0.666_r8)
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._r8,zovl )
                psih = -5._r8 * zovl
                vds  = 2.e-3_r8*ustar(i,lt)
             end if
             dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
             rds(i,lt) = 1._r8/vds
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    species_loop1 :  do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglt,endlt
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._r8
                   else
                      rmx = 1._r8/(heff(i,m)/3000._r8 + 100._r8*foxd(m))
                   end if
                   cts(i) = 1000._r8*exp( - tc(i) - 4._r8 )                 ! correction for frost
                   rgsx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                         rgsx(i,lt,ispec) = large_value
                      end if
                   end if
                   if( lt == 7 ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      if ( has_dew(i) .or. has_rain(i) ) then
                         dewm = 3._r8
                      else
                         dewm = 1._r8
                      end if
                      rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                      !-------------------------------------------------------------------------------------
                      ! jfl : special case for PAN
                      !-------------------------------------------------------------------------------------
                      if( ispec == pan_ndx .or. ispec == xpan_ndx ) then
                         dv_pan =  c0_pan(lt) * (1._r8 - exp( -k_pan(lt)*(dewm*rs*drat(m))*1.e-2_r8 ))
                         if( dv_pan > 0._r8 .and. sndx /= 4 ) then
                            rsmx(i,lt,ispec) = ( 1._r8/dv_pan )
                         end if
                      end if
                      rclx(i,lt,ispec) = cts(i) + 1._r8/((heff(i,m)/(1.e5_r8*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_r8*heff(i,m) + foxd(m))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1

    do lt = beglt,endlt
       if( lt /= 7 ) then
          do i = 1,ncol
             if( fr_lnduse(i,lt) ) then
                sndx = index_season(i,lt)
                !-------------------------------------------------------------------------------------
                ! 	... no effect if sfc_temp < O C
                !-------------------------------------------------------------------------------------
                if( sfc_temp(i) > tmelt ) then
                   if( has_dew(i) ) then
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                   if( has_rain(i) ) then
                      ! rlux(i,lt,o3_ndx) = 1./(1.e-3 + (1./(3.*rlu(sndx,lt))))
                      rlux_o3(i,lt)     = 3000._r8*rlu(sndx,lt)/(1000._r8 + 3._r8*rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                end if

                if ( o3_ndx > 0 ) then
                   rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
                end if
                if ( o3a_ndx > 0 ) then
                   rclx(i,lt,o3a_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3a_ndx) = cts(i) + rlux(i,lt,o3a_ndx)
                end if

             end if
          end do
       end if
    end do

    species_loop2 : do ispec = 1,gas_pcnst
       m = map_dvel(ispec)
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( has_dew(i) ) then
                               rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                                    + 1.e-7_r8*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                            end if
                         end if

                      end if
                   end do
                end if
             end do
          else if( ispec == so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( qs(i) <= spec_hum(i) ) then
                               rlux(i,lt,ispec) = 100._r8
                            end if
                            if( has_rain(i) ) then
                               !                               rlux(i,lt,ispec) = 1./(2.e-4 + (1./(3.*rlu(index_season(i,lt),lt))))
                               rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                            end if
                         end if
                         rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                         rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                      end if
                   end do
                end if
             end do
             do i = 1,ncol
                if( fr_lnduse(i,1) .and. (has_dew(i) .or. has_rain(i)) ) then
                   rlux(i,1,ispec) = 50._r8
                end if
             end do
          end if
       end if
    end do species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))
    species_loop3 : do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          wrk(:) = 0._r8
          lt_loop: do lt = beglt,endlt
             do i = 1,ncol
                if (fr_lnduse(i,lt)) then
                   resc(i) = 1._r8/( 1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                                   + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                                   + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

                   resc(i) = max( 10._r8,resc(i) )

                   lnd_frc(i) = lcl_frc_landuse(i,lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
             select case( solsym(ispec) )
             case( 'SO2' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      ! assume no surface resistance for SO2 over water`
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk))
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:))
                   endwhere
                end if

                !  JFL - increase in dry deposition of SO2 to improve bias over US/Europe
                wrk(:) = wrk(:) * 2._r8

             case( 'SO4' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
                endwhere
             case( 'NH4', 'NH4NO3', 'XNH4NO3' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + 0.5_r8*rds(:,lt))
                endwhere

             !-------------------------------------------------------------------------------------
             !  ... special case for Pb (for consistency with offline code)
             !-------------------------------------------------------------------------------------
             case( 'Pb' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_r8
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if

             !-------------------------------------------------------------------------------------
             !  ... special case for carbon aerosols
             !-------------------------------------------------------------------------------------
             case( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB','SOAX' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.10e-2_r8
                endwhere

             !-------------------------------------------------------------------------------------
             ! deposition over ocean for HCN, CH3CN
             !    velocity estimated from aircraft measurements (E.Apel, INTEX-B)
             !-------------------------------------------------------------------------------------
             case( 'HCN','CH3CN' )
                if( lt == 7 ) then ! over ocean only
                   where( fr_lnduse(:ncol,lt) .and. snow(:ncol) < 0.01_r8  )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
                   endwhere
                end if
             case default
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
                endwhere
             end select
          end do lt_loop
          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
       end if

    end do species_loop3

    if ( beglt > 1 ) return

    !-------------------------------------------------------------------------------------
    ! 	... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx > 0 ) then
       if( has_dvel(mpan_ndx) ) then
          dvel(:ncol,mpan_ndx) = dvel(:ncol,mpan_ndx)/3._r8
          dflx(:ncol,mpan_ndx) = term(:ncol) * dvel(:ncol,mpan_ndx) * mmr(:ncol,plev,mpan_ndx)
       end if
    end if
    if( xmpan_ndx > 0 ) then
       if( has_dvel(xmpan_ndx) ) then
          dvel(:ncol,xmpan_ndx) = dvel(:ncol,xmpan_ndx)/3._r8
          dflx(:ncol,xmpan_ndx) = term(:ncol) * dvel(:ncol,xmpan_ndx) * mmr(:ncol,plev,xmpan_ndx)
       end if
    end if

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0) then
       if( has_dvel(hcooh_ndx) ) then
          dvel(:ncol,hcooh_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,hcooh_ndx) = term(:ncol) * dvel(:ncol,hcooh_ndx) * mmr(:ncol,plev,hcooh_ndx)
       end if
    end if
!
! SOG species
!
    if( sogm_ndx > 0) then
       if( has_dvel(sogm_ndx) ) then
          dvel(:ncol,sogm_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogm_ndx) = term(:ncol) * dvel(:ncol,sogm_ndx) * mmr(:ncol,plev,sogm_ndx)
       end if
    end if
    if( sogi_ndx > 0) then
       if( has_dvel(sogi_ndx) ) then
          dvel(:ncol,sogi_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogi_ndx) = term(:ncol) * dvel(:ncol,sogi_ndx) * mmr(:ncol,plev,sogi_ndx)
       end if
    end if
    if( sogt_ndx > 0) then
       if( has_dvel(sogt_ndx) ) then
          dvel(:ncol,sogt_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogt_ndx) = term(:ncol) * dvel(:ncol,sogt_ndx) * mmr(:ncol,plev,sogt_ndx)
       end if
    end if
    if( sogb_ndx > 0) then
       if( has_dvel(sogb_ndx) ) then
          dvel(:ncol,sogb_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogb_ndx) = term(:ncol) * dvel(:ncol,sogb_ndx) * mmr(:ncol,plev,sogb_ndx)
       end if
    end if
    if( sogx_ndx > 0) then
       if( has_dvel(sogx_ndx) ) then
          dvel(:ncol,sogx_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogx_ndx) = term(:ncol) * dvel(:ncol,sogx_ndx) * mmr(:ncol,plev,sogx_ndx)
       end if
    end if

  end subroutine drydep_xactive

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  function has_drydep( name )

    character(len=*), intent(in) :: name

    logical :: has_drydep
    integer :: i

    has_drydep = .false.

    do i=1,nddvels
       if ( trim(name) == trim(drydep_list(i)) ) then
         has_drydep = .true.
         exit
       endif
    enddo

  endfunction has_drydep

end module mo_drydep
