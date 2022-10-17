!----------------------------------------------------------------------------------
! circulation diagnostics -- terms of the Transformed Eulerian Mean (TEM) equation
!
!----------------------------------------------------------------------------------
module phys_grid_ctem
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: begchunk, endchunk, pcols, pver, pverp
  use ref_pres,      only: pref_edge
  use interpolate_data, only: vertinterp
  use physics_types, only: physics_state
  use cam_history,   only: addfld, outfld
  use zonal_mean_mod,only: ZonalAverage_t, ZonalMean_t
  use physconst,     only: pi
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun
  use namelist_utils,only: find_group_name
  use spmd_utils,    only: masterproc, mpi_integer, masterprocid, mpicom

  implicit none

  private
  public :: phys_grid_ctem_readnl
  public :: phys_grid_ctem_reg
  public :: phys_grid_ctem_init
  public :: phys_grid_ctem_diags

  type(ZonalMean_t) :: ZMobj
  type(ZonalAverage_t) :: ZAobj

  integer :: nzalat = -huge(1)
  integer :: nzmbas = -huge(1)

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine phys_grid_ctem_readnl(nlfile)
    character(len=*), intent(in) :: nlfile
    integer :: ierr, unitn

    character(len=*), parameter :: prefix = 'phys_grid_ctem_readnl: '
    integer :: phys_grid_ctem_zm_nbas
    integer :: phys_grid_ctem_za_nlat

    namelist /phys_grid_ctem_opts/ phys_grid_ctem_zm_nbas, phys_grid_ctem_za_nlat

    ! Read in namelist values
    !------------------------
    if(masterproc) then
       open(newunit=unitn, file=trim(nlfile), status='old')
       call find_group_name(unitn, 'phys_grid_ctem_opts', status=ierr)
       if(ierr == 0) then
          read(unitn,phys_grid_ctem_opts,iostat=ierr)
          if(ierr /= 0) then
             call endrun(prefix//'ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call MPI_bcast(phys_grid_ctem_zm_nbas, 1, mpi_integer, masterprocid, mpicom, ierr)
    call MPI_bcast(phys_grid_ctem_za_nlat, 1, mpi_integer, masterprocid, mpicom, ierr)

    if (masterproc) then
       write(iulog,*) 'phys_grid_ctem_readnl... phys_grid_ctem_zm_nbas: ',phys_grid_ctem_zm_nbas
       write(iulog,*) 'phys_grid_ctem_readnl... phys_grid_ctem_za_nlat: ',phys_grid_ctem_za_nlat
    endif

    nzalat = phys_grid_ctem_za_nlat
    nzmbas = phys_grid_ctem_zm_nbas

  end subroutine phys_grid_ctem_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine phys_grid_ctem_reg

    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register

    type(horiz_coord_t), pointer :: zalon_coord
    type(horiz_coord_t), pointer :: zalat_coord
    integer(iMap),       pointer :: grid_map(:,:)

    real(r8) :: zalats(nzalat)
    real(r8) :: area(nzalat)
    real(r8) :: zalons(1)
    real(r8) :: dlatrad, dlatdeg, lat1, lat2
    real(r8) :: total_area
    real(r8) :: total_wght
    integer :: j

    real(r8), parameter :: latdeg0 = -90._r8
    real(r8), parameter :: latrad0 = -pi*0.5_r8
    real(r8), parameter :: fourpi = pi*4._r8

    integer, parameter :: ctem_zavg_phys_decomp = 201 ! Must be unique within CAM

    nullify(zalat_coord)
    nullify(zalon_coord)
    nullify(grid_map)

    zalons(1) = 0._r8

    dlatrad = pi/real(nzalat,kind=r8)
    dlatdeg = 180._r8/real(nzalat,kind=r8)
    total_area = 0._r8
    total_wght = 0._r8

    do j = 1,nzalat
       zalats(j) = latdeg0 + (real(j,kind=r8)-0.5_r8)*dlatdeg
       lat1 = latrad0 + real(j-1,kind=r8)*dlatrad
       lat2 = latrad0 + real(j  ,kind=r8)*dlatrad
       area(j) = 2._r8*pi*(sin(lat2)-sin(lat1))
       total_area = total_area + area(j)
       total_wght = total_wght + 0.5_r8*(sin(lat2)-sin(lat1))
    end do

    if ( abs(1._r8-total_wght)>1.e-12_r8 .or. abs(fourpi-total_area)>1.e-12_r8 ) then
       call endrun('zmean_phys_fields_reg: problem with area/wght calc')
    end if

    call ZAobj%init(zalats,area,nzalat,GEN_GAUSSLATS=.false.)
    call ZMobj%init(nzmbas)

    ! Zonal average grid

    zalat_coord => horiz_coord_create('zalat', '', nzalat, 'latitude',                &
         'degrees_north', 1, nzalat, zalats)
    zalon_coord => horiz_coord_create('zalon', '', 1, 'longitude',                &
         'degrees_east', 1, 1, zalons)

    ! grid decomposition map
    allocate(grid_map(4,nzalat))

    do j = 1,nzalat
       grid_map(1,j) = 1
       grid_map(2,j) = j
       if (masterproc) then
          grid_map(3,j) = 1
          grid_map(4,j) = j
       else
          grid_map(3,j) = 0
          grid_map(4,j) = 0
       end if
    end do

    ! register the zonal average grid
    call cam_grid_register('ctem_zavg_phys', ctem_zavg_phys_decomp, zalat_coord, &
         zalon_coord, grid_map, unstruct=.false., zonal_grid=.true.)

  end subroutine phys_grid_ctem_reg

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine phys_grid_ctem_init

    call addfld ('VTHzaphys',(/'ilev'/), 'A', 'MK/S', 'Meridional Heat Flux:', gridname='ctem_zavg_phys')
    call addfld ('WTHzaphys',(/'ilev'/), 'A', 'MK/S', 'Vertical Heat Flux:', gridname='ctem_zavg_phys')
    call addfld ('UVzaphys', (/'ilev'/), 'A', 'M2/S2','Meridional Flux of Zonal Momentum', gridname='ctem_zavg_phys')
    call addfld ('UWzaphys', (/'ilev'/), 'A', 'M2/S2','Vertical Flux of Zonal Momentum', gridname='ctem_zavg_phys')

    call addfld ('VTHzm3d',(/'ilev' /), 'A','MK/S', 'Meridional Heat Flux: 3D zon. mean', gridname='physgrid' )
    call addfld ('WTHzm3d',(/'ilev' /), 'A','MK/S', 'Vertical Heat Flux: 3D zon. mean', gridname='physgrid' )
    call addfld ('UVzm3d', (/'ilev' /), 'A','M2/S2','Meridional Flux of Zonal Momentum: 3D zon. mean', gridname='physgrid' )
    call addfld ('UWzm3d', (/'ilev' /), 'A','M2/S2','Vertical Flux of Zonal Momentum: 3D zon. mean', gridname='physgrid' )

    call addfld ('Uzm3d',  (/'ilev' /), 'A','M/S',  'Zonal-Mean zonal wind - defined on ilev', gridname='physgrid')
    call addfld ('Vzm3d',  (/'ilev' /), 'A','M/S',  'Zonal-Mean meridional wind - defined on ilev', gridname='physgrid' )
    call addfld ('Wzm3d',  (/'ilev' /), 'A','M/S',  'Zonal-Mean vertical wind - defined on ilev', gridname='physgrid' )
    call addfld ('THzm3d', (/'ilev' /), 'A',  'K',  'Zonal-Mean potential temp - defined on ilev', gridname='physgrid' )
    call addfld ('THphys', (/'ilev' /), 'A',  'K',  'Zonal-Mean potential temp - defined on ilev', gridname='physgrid' )

  end subroutine phys_grid_ctem_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine phys_grid_ctem_diags(phys_state)
    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8) :: ui(pcols,pverp,begchunk:endchunk)
    real(r8) :: vi(pcols,pverp,begchunk:endchunk)
    real(r8) :: wi(pcols,pverp,begchunk:endchunk)

    real(r8) :: uzm(pcols,pverp,begchunk:endchunk)
    real(r8) :: vzm(pcols,pverp,begchunk:endchunk)
    real(r8) :: wzm(pcols,pverp,begchunk:endchunk)

    real(r8) :: ud(pcols,pverp,begchunk:endchunk)
    real(r8) :: vd(pcols,pverp,begchunk:endchunk)
    real(r8) :: wd(pcols,pverp,begchunk:endchunk)
    real(r8) :: thd(pcols,pverp,begchunk:endchunk)

    real(r8) :: uvp(pcols,pverp,begchunk:endchunk)
    real(r8) :: uwp(pcols,pverp,begchunk:endchunk)
    real(r8) :: vthp(pcols,pverp,begchunk:endchunk)
    real(r8) :: wthp(pcols,pverp,begchunk:endchunk)

    real(r8) :: uv(pcols,pverp,begchunk:endchunk)
    real(r8) :: uw(pcols,pverp,begchunk:endchunk)
    real(r8) :: vth(pcols,pverp,begchunk:endchunk)
    real(r8) :: wth(pcols,pverp,begchunk:endchunk)

    integer  :: lchnk, ncol, j, k
    real(r8) :: fld_tmp(pcols,pverp)

    real(r8) :: theta(pcols,pver,begchunk:endchunk) ! potential temperature
    real(r8) :: thi(pcols,pverp,begchunk:endchunk)
    real(r8) :: thzm(pcols,pverp,begchunk:endchunk)

    real(r8) :: w(pcols,pver,begchunk:endchunk)

    real(r8) :: uvza(nzalat,pverp)
    real(r8) :: uwza(nzalat,pverp)
    real(r8) :: vthza(nzalat,pverp)
    real(r8) :: wthza(nzalat,pverp)

    real(r8), parameter :: hscale = 7000._r8          ! pressure scale height

    ui(:,:,:) = 0._r8
    vi(:,:,:) = 0._r8
    wi(:,:,:) = 0._r8
    thi(:,:,:) = 0._r8

    uzm(:,:,:) = 0._r8
    vzm(:,:,:) = 0._r8
    wzm(:,:,:) = 0._r8
    thzm(:,:,:) = 0._r8

    ud(:,:,:) = 0._r8
    vd(:,:,:) = 0._r8
    uvp(:,:,:) = 0._r8

    do lchnk = begchunk,endchunk

       ncol = phys_state(lchnk)%ncol

       theta(:ncol,:,lchnk) = phys_state(lchnk)%t(:ncol,:) * phys_state(lchnk)%exner(:ncol,:)
       w(:ncol,:,lchnk)  = - hscale *  phys_state(lchnk)%omega(:ncol,:) / phys_state(lchnk)%pmid(:ncol,:)

       do k = 1,pverp
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), phys_state(lchnk)%u(:,:), ui(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), phys_state(lchnk)%v(:,:), vi(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), theta(:,:,lchnk), thi(:,k,lchnk) )
          call vertinterp( ncol, pcols, pver, phys_state(lchnk)%pmid(:,:), pref_edge(k), w(:,:,lchnk), wi(:,k,lchnk) )
       end do

    end do

    ! these need to be evaluated on the physics grid (3D)
    ! to be used in the deviations calculation below
    uzm(:,:,:) = zmean_fld(ui(:,:,:))
    vzm(:,:,:) = zmean_fld(vi(:,:,:))
    wzm(:,:,:) = zmean_fld(wi(:,:,:))
    thzm(:,:,:) = zmean_fld(thi(:,:,:))

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol

       fld_tmp(:ncol,:) = thi(:ncol,:,lchnk)
       call outfld( 'THphys', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = thzm(:ncol,:,lchnk)
       call outfld( 'THzm3d', fld_tmp(:ncol,:), ncol, lchnk)

       fld_tmp(:ncol,:) = uzm(:ncol,:,lchnk)
       call outfld( 'Uzm3d', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = vzm(:ncol,:,lchnk)
       call outfld( 'Vzm3d', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = wzm(:ncol,:,lchnk)
       call outfld( 'Wzm3d', fld_tmp(:ncol,:), ncol, lchnk)
    end do


    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do k = 1,pverp
          ! zonal deviations
          thd(:ncol,k,lchnk) = thi(:ncol,k,lchnk) - thzm(:ncol,k,lchnk)
          ud(:ncol,k,lchnk) = ui(:ncol,k,lchnk) - uzm(:ncol,k,lchnk)
          vd(:ncol,k,lchnk) = vi(:ncol,k,lchnk) - vzm(:ncol,k,lchnk)
          wd(:ncol,k,lchnk) = wi(:ncol,k,lchnk) - wzm(:ncol,k,lchnk)
          ! fluxes
          uvp(:ncol,k,lchnk) = ud(:ncol,k,lchnk) * vd(:ncol,k,lchnk)
          uwp(:ncol,k,lchnk) = ud(:ncol,k,lchnk) * wd(:ncol,k,lchnk)
          vthp(:ncol,k,lchnk) = vd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
          wthp(:ncol,k,lchnk) = wd(:ncol,k,lchnk) * thd(:ncol,k,lchnk)
       end do
    end do

    ! evaluate and output these on the zonal-average grid
    call ZAobj%binAvg(uvp, uvza)
    call ZAobj%binAvg(uwp, uwza)
    call ZAobj%binAvg(vthp, vthza)
    call ZAobj%binAvg(wthp, wthza)

    do j = 1,nzalat
       call outfld('UVzaphys',uvza(j,:),1,j)
       call outfld('UWzaphys',uwza(j,:),1,j)
       call outfld('VTHzaphys',vthza(j,:),1,j)
       call outfld('WTHzaphys',wthza(j,:),1,j)
    end do

    ! not needed --- only for sanity checks
    uv(:,:,:) = zmean_fld(uvp(:,:,:))
    uw(:,:,:) = zmean_fld(uwp(:,:,:))
    vth(:,:,:) = zmean_fld(vthp(:,:,:))
    wth(:,:,:) = zmean_fld(wthp(:,:,:))

    do lchnk = begchunk, endchunk
       ncol = phys_state(lchnk)%ncol
       fld_tmp(:ncol,:) = uv(:ncol,:,lchnk)
       call outfld( 'UVzm3d', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = uw(:ncol,:,lchnk)
       call outfld( 'UWzm3d', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = vth(:ncol,:,lchnk)
       call outfld( 'VTHzm3d', fld_tmp(:ncol,:), ncol, lchnk)
       fld_tmp(:ncol,:) = wth(:ncol,:,lchnk)
       call outfld( 'WTHzm3d', fld_tmp(:ncol,:), ncol, lchnk)
    end do

  contains

    !------------------------------------------------------------------------------
    ! utility function for evaluating 3D zonal mean fields
    !------------------------------------------------------------------------------
    function zmean_fld( fld ) result(fldzm)

      real(r8), intent(in) :: fld(pcols,pverp,begchunk:endchunk)

      real(r8) :: fldzm(pcols,pverp,begchunk:endchunk)

      real(r8) :: Zonal_Bamp3d(nzmbas,pverp)

      call ZMobj%calc_amps(fld,Zonal_Bamp3d)
      call ZMobj%eval_grid(Zonal_Bamp3d,fldzm)

    end function zmean_fld

  end subroutine phys_grid_ctem_diags

end module phys_grid_ctem
