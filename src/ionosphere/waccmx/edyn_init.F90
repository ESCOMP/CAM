module edyn_init
!
! Initialize edynamo
!
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use spmd_utils,     only: masterproc
   use edynamo,        only: alloc_edyn, debug_hist

   implicit none

   private
   public :: edynamo_init

contains
!-----------------------------------------------------------------------
  subroutine edynamo_init(mpicomm,ionos_debug_hist)

      !
      ! One-time initialization, called from ionosphere_init
      ! before dyn_init and phys_init
      !
      use edyn_maggrid, only: set_maggrid, gmlat, nmlonp1, nmlat, nmlath, nmlev
      use edyn_mpi,     only: mp_exchange_tasks
      use edyn_mpi,     only: mp_distribute_mag
      use edyn_grid_comp, only: edyn_grid_comp_init
      use edyn_solve, only: edyn_solve_init

      !
      ! Args:
      integer, intent(in) :: mpicomm
      logical, intent(in) :: ionos_debug_hist

      debug_hist = ionos_debug_hist

      if (masterproc) then
         write(iulog,"('Enter edynamo_init:')")
      endif

      call set_maggrid ()   ! set parameter-based global magnetic grid

      call edyn_solve_init

      call mp_distribute_mag(nmlonp1, nmlat, nmlath, nmlev)

      call register_grids()
      call mp_exchange_tasks(mpicomm, 0, gmlat) ! single arg is iprint

      call alloc_edyn()      ! allocate dynamo arrays
      call edyn_grid_comp_init(mpicomm)

      call add_fields()      ! add fields to WACCM history master list

   end subroutine edynamo_init

   !-----------------------------------------------------------------------
   subroutine add_fields
      use cam_history,  only: addfld, horiz_only, add_default
      use phys_control, only: phys_getopts !Method used to get flag for waccmx ionosphere output variables

      logical :: history_waccmx

      call addfld ('PED_MAG'   ,(/ 'lev' /), 'I', 'S/m  ','Pedersen Conductivity'            ,gridname='gmag_grid')
      call addfld ('HAL_MAG'   ,(/ 'lev' /), 'I', 'S/m  ','Hall Conductivity'                ,gridname='gmag_grid')
      call addfld ('PHIHM'     , horiz_only, 'I', 'VOLTS','High Latitude Electric Potential' ,gridname='gmag_grid')
      call addfld ('PHIM2D'    , horiz_only, 'I', 'VOLTS','PHIM2D: Electric Potential'       ,gridname='gmag_grid')
      call addfld ('ED1'       , horiz_only, 'I', 'V/m  ','ED1: Eastward Electric Field'     ,gridname='gmag_grid')
      call addfld ('ED2'       , horiz_only, 'I', 'V/m  ','ED2: Equatorward Electric Field'  ,gridname='gmag_grid')
      call addfld ('PHIM3D'    ,(/ 'lev' /), 'I', 'VOLTS','PHIM3D: 3d Electric Potential'    ,gridname='gmag_grid')
      call addfld ('ED13D'     ,(/ 'lev' /), 'I', 'V/m  ','ED13D: Eastward Electric Field'   ,gridname='gmag_grid')
      call addfld ('ED23D'     ,(/ 'lev' /), 'I', 'V/m  ','ED23D: Equatorward Electric Field',gridname='gmag_grid')
      call addfld ('ZPOT_MAG'  ,(/ 'lev' /), 'I', 'cm   ','Geopotential on mag grid (h0 min)',gridname='gmag_grid')

      call addfld ('HALL_CONDUCTANCE',horiz_only, 'I', 'S','Hall Conductance',    gridname='gmag_grid')
      call addfld ('PED_CONDUCTANCE', horiz_only, 'I', 'S','Pedersen Conductance',gridname='gmag_grid')

      call addfld ('POTEN',(/ 'lev' /), 'I', 'Volts','POTEN: Electric Potential', gridname='geo_grid')

      call addfld ('EX',(/'lev'/),'I','V/m','Geographic zonal component of Electric Field',  gridname='geo_grid')
      call addfld ('EY',(/'lev'/),'I','V/m','Geographic meridional component of Electric Field', gridname='geo_grid')
      call addfld ('EZ',(/'lev'/),'I','V/m','Vertical component of Electric Field', gridname='geo_grid')

      call addfld ('BMOD', horiz_only, 'I', 'gauss','magnitude of magnetic field',gridname='geo_grid')
      call addfld ('XB',   horiz_only, 'I', 'gauss','northward component of magnetic field',gridname='geo_grid')
      call addfld ('YB',   horiz_only, 'I', 'gauss','eastward component of magnetic field',gridname='geo_grid')
      call addfld ('ZB',   horiz_only, 'I', 'gauss','downward component of magnetic field',gridname='geo_grid')

      if (debug_hist) then
         call addfld ('EPHI3D'    ,(/ 'lev' /), 'I', ' ','EPHI3D'    ,gridname='gmag_grid')
         call addfld ('ELAM3D'    ,(/ 'lev' /), 'I', ' ','ELAM3D'    ,gridname='gmag_grid')
         call addfld ('EMZ3D'     ,(/ 'lev' /), 'I', ' ','EMZ3D'     ,gridname='gmag_grid')

         call addfld ('ADOTV1_MAG',(/ 'lev' /), 'I', '     ','ADOTV1 on mag grid'               ,gridname='gmag_grid')
         call addfld ('ADOTV2_MAG',(/ 'lev' /), 'I', '     ','ADOTV2 on mag grid'               ,gridname='gmag_grid')
         call addfld ('adota1_mag_a', horiz_only, 'I', ' ','adota1_mag_a',gridname='gmag_grid')
         call addfld ('ZIGM11_a', horiz_only, 'I', ' ','ZIGM11_a',gridname='gmag_grid')
         call addfld ('EDYN_ZIGM11_0', horiz_only, 'I', ' ','EDYN_ZIGM11_0',gridname='gmag_grid')
         call addfld ('EDYN_ZIGM11', horiz_only, 'I', ' ','EDYN_ZIGM11',gridname='gmag_grid')
         call addfld ('EDYN_ZIGM22', horiz_only, 'I', ' ','EDYN_ZIGM22',gridname='gmag_grid')
         call addfld ('EDYN_ZIGMC' , horiz_only, 'I', ' ','EDYN_ZIGMC' ,gridname='gmag_grid')
         call addfld ('EDYN_ZIGM2' , horiz_only, 'I', ' ','EDYN_ZIGM2' ,gridname='gmag_grid')
         call addfld ('EDYN_RIM1'  , horiz_only, 'I', ' ','EDYN_RIM1'  ,gridname='gmag_grid')
         call addfld ('EDYN_RIM2'  , horiz_only, 'I', ' ','EDYN_RIM2'  ,gridname='gmag_grid')

         ! rjac: scaled derivatives of geomagnetic coords wrt geographic coordinates.
         call addfld ('RJAC11',(/'lev'/), 'I', '1','cos(thetas)/cos(theta)*d(lamdas)/d(lamda)'  ,gridname='geo_grid')
         call addfld ('RJAC12',(/'lev'/), 'I', '1','cos(thetas)*d(lamdas)/d(theta)'  ,gridname='geo_grid')
         call addfld ('RJAC21',(/'lev'/), 'I', '1','1./cos(theta)*d(thetas)/d(lamda)'  ,gridname='geo_grid')
         call addfld ('RJAC22',(/'lev'/), 'I', '1','d(thetas)/d(theta)'  ,gridname='geo_grid')
      endif
      !-------------------------------------------------------------------------------
      !  Set default values for ionosphere history variables
      !-------------------------------------------------------------------------------
      call phys_getopts(history_waccmx_out=history_waccmx)

      if (history_waccmx) then
         call add_default ('PED_CONDUCTANCE', 1, ' ')
         call add_default ('HALL_CONDUCTANCE' , 1, ' ')
      end if

   end subroutine add_fields
   !-----------------------------------------------------------------------

   subroutine register_grids()

      use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
      use cam_grid_support, only: cam_grid_register
      use edyn_mpi,         only: mlat0, mlat1, mlon0, omlon1, ntask, mytid
      use edyn_mpi,         only:  lat0,  lat1,  lon0,   lon1
      use edyn_maggrid,     only: gmlat, gmlon, nmlat, nmlon
      use edyn_geogrid,     only:  glat,  glon,  nlat,  nlon

      integer, parameter :: mag_decomp = 111 ! Must be unique within CAM
      integer, parameter :: geo_decomp = 112 ! Must be unique within CAM

      type(horiz_coord_t), pointer :: lat_coord => null()
      type(horiz_coord_t), pointer :: lon_coord => null()
      integer(iMap),       pointer :: grid_map(:,:) => null()
      integer(iMap),       pointer :: coord_map(:) => null()
      integer                      :: i, j, ind

      if (mytid>=ntask) then

         if (mlon0/=1) then
            call endrun('register_grids: mlon0 needs to be 1 on inactive PEs')
         end if
         if (omlon1/=0) then
            call endrun('register_grids: omlon1 needs to be 0 on inactive PEs')
         end if
         if (mlat0/=1) then
            call endrun('register_grids: mlat0 needs to be 1 on inactive PEs')
         end if
         if (mlat1/=0) then
            call endrun('register_grids: mlat1 needs to be 0 on inactive PEs')
         end if

         if (lon0/=1) then
            call endrun('register_grids: lon0 needs to be 1 on inactive PEs')
         end if
         if (lon1/=0) then
            call endrun('register_grids: lon1 needs to be 0 on inactive PEs')
         end if
         if (lat0/=1) then
            call endrun('register_grids: lat0 needs to be 1 on inactive PEs')
         end if
         if (lat1/=0) then
            call endrun('register_grids: lat1 needs to be 0 on inactive PEs')
         end if

      endif

      allocate(grid_map(4, ((omlon1 - mlon0 + 1) * (mlat1 - mlat0 + 1))))
      ind = 0
      do i = mlat0, mlat1
         do j = mlon0, omlon1
            ind = ind + 1
            grid_map(1, ind) = j
            grid_map(2, ind) = i
            grid_map(3, ind) = j
            grid_map(4, ind) = i
         end do
      end do

      allocate(coord_map(mlat1 - mlat0 + 1))
      if (mlon0==1) then
         coord_map = (/ (i, i = mlat0, mlat1) /)
      else
         coord_map = 0
      end if
      lat_coord => horiz_coord_create('mlat', '', nmlat, 'latitude',           &
           'degrees_north', mlat0, mlat1, gmlat(mlat0:mlat1),                  &
           map=coord_map)
      nullify(coord_map)

      allocate(coord_map(omlon1 - mlon0 + 1))
      if (mlat0==1) then
         coord_map = (/ (i, i = mlon0, omlon1) /)
      else
         coord_map = 0
      end if
      lon_coord => horiz_coord_create('mlon', '', nmlon, 'longitude',          &
           'degrees_east', mlon0, omlon1, gmlon(mlon0:omlon1),                 &
           map=coord_map)
      deallocate(coord_map)
      nullify(coord_map)

      call cam_grid_register('gmag_grid', mag_decomp, lat_coord, lon_coord,    &
           grid_map, unstruct=.false.)
      nullify(grid_map)


      ! for the Oplus geo grid
      allocate(grid_map(4, ((lon1 - lon0 + 1) * (lat1 - lat0 + 1))))
      ind = 0
      do i = lat0, lat1
         do j = lon0, lon1
            ind = ind + 1
            grid_map(1, ind) = j
            grid_map(2, ind) = i
            grid_map(3, ind) = j
            grid_map(4, ind) = i
         end do
      end do

      allocate(coord_map(lat1 - lat0 + 1))
      if (lon0==1) then
         coord_map = (/ (i, i = lat0, lat1) /)
      else
         coord_map = 0
      end if
      lat_coord => horiz_coord_create('glat', '', nlat, 'latitude',           &
           'degrees_north', lat0, lat1, glat(lat0:lat1),                  &
           map=coord_map)
      nullify(coord_map)

      allocate(coord_map(lon1 - lon0 + 1))
      if (lat0==1) then
         coord_map = (/ (i, i = lon0, lon1) /)
      else
         coord_map = 0
      end if

      lon_coord => horiz_coord_create('glon', '', nlon, 'longitude',          &
           'degrees_east', lon0, lon1, glon(lon0:lon1),                 &
           map=coord_map)
      deallocate(coord_map)
      nullify(coord_map)

      call cam_grid_register('geo_grid', geo_decomp, lat_coord, lon_coord,    &
           grid_map, unstruct=.false.)
      nullify(grid_map)

   end subroutine register_grids

!-----------------------------------------------------------------------
end module edyn_init
