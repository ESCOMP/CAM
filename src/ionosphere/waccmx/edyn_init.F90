module edyn_init
!
! Initialize edynamo
!
   use shr_kind_mod,   only: r8 => shr_kind_r8 ! 8-byte reals
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use spmd_utils,     only: masterproc

   implicit none

   private
   public :: edynamo_init, lonshift_global

contains
!-----------------------------------------------------------------------
  subroutine edynamo_init(mpicomm)

      !
      ! One-time initialization, called from ionosphere_init
      ! before dyn_init and phys_init
      !
      use edyn_maggrid, only: set_maggrid, gmlat, nmlonp1, nmlat, nmlath, nmlev
      use edyn_mpi,     only: mp_exchange_tasks
      use edyn_mpi,     only: mp_distribute_mag
      use edynamo,      only: alloc_edyn
      use edyn_grid_comp, only: edyn_grid_comp_init
      !
      ! Args:
      integer, intent(in) :: mpicomm

      if (masterproc) then
         write(iulog,"('Enter edynamo_init:')")
      endif

      call set_maggrid ()   ! set parameter-based global magnetic grid

      call mp_distribute_mag(nmlonp1, nmlat, nmlath, nmlev)

      call register_grids()
      call mp_exchange_tasks(mpicomm, 0, gmlat) ! single arg is iprint

      call alloc_edyn()      ! allocate dynamo arrays
      call edyn_grid_comp_init(mpicomm)

      call add_fields()      ! add fields to WACCM history master list

   end subroutine edynamo_init
!-----------------------------------------------------------------------
   subroutine lonshift_global(f,nlon,lonseq,iscoord)
!
! Shift longitude vector f(nlon) forward 180 degrees according to input
! string lonseq. Input f can be either arbitrary field values or
! the coordinate array itself. Shift f in the 'lonseq' manner, as follows:
!
! If lonseq='-180to180', then shift from 0->360 to -180->+180
! If lonseq='zeroto360', then shift from -180->+180 to 0->360
!
! WARNING: This routine works with WACCM-X history files, where nlon=144, 72, or 80
!          It has not been tested with other models or resolutions.
!          (e.g., there is no test for center point, its assumed to be nlon/2)
!
! Args:
      integer,intent(in) :: nlon
      real(r8),intent(inout) :: f(nlon)
      character(len=*),intent(in) :: lonseq
      logical,intent(in) :: iscoord ! if true, f is a coordinate, otherwise it is data
!
! Local:
      character(len=80) :: msg
      integer :: ihalf,i

      if (lonseq /= '-180to180'.and.lonseq /= 'zeroto360') then
         write(msg,"('shift_lon: bad lonseq=',a,' must be either ''-180to180'' or ''zeroto360''')") &
              lonseq
         call endrun
      endif

      ihalf = nlon/2
      if (lonseq == '-180to180') then ! shift to -180 -> +180
         f = cshift(f,ihalf)           ! cshift is circular shift intrinsic
         if (iscoord) then
            do i=1,ihalf
               f(i) = f(i)-360._r8
            enddo
         endif
      else                           ! shift to 0 -> 360
         f = cshift(f,ihalf)          ! cshift is circular shift intrinsic
         if (iscoord) then
            do i=ihalf+1,nlon
               f(i) = f(i)+360._r8
            enddo
         endif
      endif
   end subroutine lonshift_global
!-----------------------------------------------------------------------
   subroutine add_fields
      use cam_history,  only: addfld, horiz_only, add_default
      use phys_control, only: phys_getopts !Method used to get flag for waccmx ionosphere output variables

      logical :: history_waccmx

      ! Geomagnetic fields are in waccm format, in CGS units):
      call addfld ('PED_MAG'   ,(/ 'lev' /), 'I', 'S/m  ','Pedersen Conductivity'            ,gridname='gmag_grid')
      call addfld ('HAL_MAG'   ,(/ 'lev' /), 'I', 'S/m  ','Hall Conductivity'                ,gridname='gmag_grid')
      call addfld ('PHIM2D'    , horiz_only, 'I', 'VOLTS','PHIM2D: Electric Potential'       ,gridname='gmag_grid')
      call addfld ('ED1'       , horiz_only, 'I', 'V/m  ','ED1: Eastward Electric Field'     ,gridname='gmag_grid')
      call addfld ('ED2'       , horiz_only, 'I', 'V/m  ','ED2: Equatorward Electric Field'  ,gridname='gmag_grid')
      call addfld ('PHIM3D'    ,(/ 'lev' /), 'I', 'VOLTS','PHIM3D: 3d Electric Potential'    ,gridname='gmag_grid')

      call addfld ('EPHI3D'    ,(/ 'lev' /), 'I', ' ','EPHI3D'    ,gridname='gmag_grid')
      call addfld ('ELAM3D'    ,(/ 'lev' /), 'I', ' ','ELAM3D'    ,gridname='gmag_grid')
      call addfld ('EMZ3D'     ,(/ 'lev' /), 'I', ' ','EMZ3D'     ,gridname='gmag_grid')

      call addfld ('ED13D'     ,(/ 'lev' /), 'I', 'V/m  ','ED13D: Eastward Electric Field'   ,gridname='gmag_grid')
      call addfld ('ED23D'     ,(/ 'lev' /), 'I', 'V/m  ','ED23D: Equatorward Electric Field',gridname='gmag_grid')
      call addfld ('ZPOT_MAG'  ,(/ 'lev' /), 'I', 'cm   ','Geopotential on mag grid (h0 min)',gridname='gmag_grid')
      call addfld ('ADOTV1_MAG',(/ 'lev' /), 'I', '     ','ADOTV1 on mag grid'               ,gridname='gmag_grid')
      call addfld ('ADOTV2_MAG',(/ 'lev' /), 'I', '     ','ADOTV2 on mag grid'               ,gridname='gmag_grid')
      !
      call addfld ('prescr_phihm' , horiz_only, 'I','VOLTS','Prescribed Electric Potential-mag grid' ,gridname='gmag_grid')
      call addfld ('prescr_efxm'  , horiz_only, 'I','mW/m2','Prescribed energy flux on mag grid'     ,gridname='gmag_grid')
      call addfld ('prescr_kevm'  , horiz_only, 'I','keV  ','Prescribed mean energy on mag grid'     ,gridname='gmag_grid')

      !
      ! Dynamo inputs from sub dynamo_input (edynamo.F90):

      call addfld ('EDYN_ADOTV1 ',(/ 'lev' /), 'I', '        ','EDYN_ADOTV1' , gridname='geo_grid')
      call addfld ('EDYN_ADOTV2 ',(/ 'lev' /), 'I', '        ','EDYN_ADOTV2' , gridname='geo_grid')
      !
      ! 2d dynamo input fields on geo grid (edynamo.F90):
      call addfld ('EDYN_SINI   ', horiz_only , 'I', '        ','EDYN_SINI'   , gridname='geo_grid')
      call addfld ('EDYN_ADOTA1 ', horiz_only , 'I', '        ','EDYN_ADOTA1' , gridname='geo_grid')
      call addfld ('EDYN_ADOTA2 ', horiz_only , 'I', '        ','EDYN_ADOTA2' , gridname='geo_grid')
      call addfld ('EDYN_A1DTA2 ', horiz_only , 'I', '        ','EDYN_A1DTA2' , gridname='geo_grid')
      call addfld ('EDYN_BE3    ', horiz_only , 'I', '        ','EDYN_BE3'    , gridname='geo_grid')


      call addfld ('ADOTA1_MAG', horiz_only , 'I', ' ','ADOTA1 in geo-mag coords' , gridname='gmag_grid')
      call addfld ('SINI_MAG',   horiz_only , 'I', ' ','sini in geo-mag coords' , gridname='gmag_grid')

      call addfld ('adota1_mag_a', horiz_only, 'I', ' ','EDYN_ZIGM11',gridname='gmag_grid')
      call addfld ('ZIGM11_a', horiz_only, 'I', ' ','EDYN_ZIGM11',gridname='gmag_grid')
      call addfld ('EDYN_ZIGM11_0', horiz_only, 'I', ' ','EDYN_ZIGM11',gridname='gmag_grid')
      call addfld ('EDYN_ZIGM11', horiz_only, 'I', ' ','EDYN_ZIGM11',gridname='gmag_grid')
      call addfld ('EDYN_ZIGM11_PED', horiz_only, 'I', 'S','Pedersen Conductance',gridname='gmag_grid')
      call addfld ('EDYN_ZIGM22', horiz_only, 'I', ' ','EDYN_ZIGM22',gridname='gmag_grid')
      call addfld ('EDYN_ZIGMC' , horiz_only, 'I', ' ','EDYN_ZIGMC' ,gridname='gmag_grid')
      call addfld ('EDYN_ZIGM2' , horiz_only, 'I', ' ','EDYN_ZIGM2' ,gridname='gmag_grid')
      call addfld ('EDYN_ZIGM2_HAL' , horiz_only, 'I', 'S','Hall Conductance' ,gridname='gmag_grid')
      call addfld ('EDYN_RIM1'  , horiz_only, 'I', ' ','EDYN_RIM1'  ,gridname='gmag_grid')
      call addfld ('EDYN_RIM2'  , horiz_only, 'I', ' ','EDYN_RIM2'  ,gridname='gmag_grid')

      call addfld ('POTEN'  ,(/ 'lev' /), 'I', 'Volts','POTEN: Electric Potential',&
           gridname='geo_grid')
      call addfld ('EX'     ,(/ 'lev' /), 'I', 'V/m'  ,'EX: Zonal component of Electric Field',&
           gridname='geo_grid')
      call addfld ('EY'     ,(/ 'lev' /), 'I', 'V/m'  ,'EY: Meridional component of Electric Field',&
           gridname='geo_grid')
      call addfld ('EZ'     ,(/ 'lev' /), 'I', 'V/m'  ,'EZ: Vertical component of Electric Field',&
           gridname='geo_grid')

      call addfld ('BMOD'       , horiz_only, 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('XB'       ,   horiz_only, 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('YB'       ,   horiz_only, 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('ZB'       ,   horiz_only, 'I', '  ',' '  ,gridname='geo_grid')

      call addfld ('RJAC11'       ,(/'lev'/), 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('RJAC12'       ,(/'lev'/), 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('RJAC21'       ,(/'lev'/), 'I', '  ',' '  ,gridname='geo_grid')
      call addfld ('RJAC22'       ,(/'lev'/), 'I', '  ',' '  ,gridname='geo_grid')

      call addfld ('OPLUS',   (/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')
      call addfld ('OPtm1i',   (/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')
      call addfld ('OPtm1o',   (/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)',    gridname='geo_grid')

      call addfld ('TN_phys',(/ 'lev' /), 'I', 'deg K   ','T'   , gridname='physgrid')
      call addfld ('TN_geo', (/ 'lev' /), 'I', 'deg K   ','T'   , gridname='geo_grid')

      call addfld ('Te_phys',(/ 'lev' /), 'I', 'deg K   ','T'   , gridname='physgrid')
      call addfld ('Te_geo', (/ 'lev' /), 'I', 'deg K   ','T'   , gridname='geo_grid')

      call addfld ('Ti_phys',(/ 'lev' /), 'I', 'deg K   ','T'   , gridname='physgrid')
      call addfld ('Ti_geo', (/ 'lev' /), 'I', 'deg K   ','T'   , gridname='geo_grid')

      call addfld ('ZPOT_phys',(/ 'lev' /), 'I', 'cm','T'   , gridname='physgrid')
      call addfld ('ZPOT_geo', (/ 'lev' /), 'I', 'cm','T'   , gridname='geo_grid')

      call addfld ('U_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('U_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('V_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('V_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('W_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('W_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('Ui_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('Ui_geo',(/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('Vi_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('Vi_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('Wi_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('Wi_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('OMEGA_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('OMEGA_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('O2_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('O2_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('N2_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('N2_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('O_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('O_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('Op_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('Op_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')
      call addfld ('Op_geo2', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('Optm1_phys',(/ 'lev' /), 'I', 'cm-3',' '   , gridname='physgrid')
      call addfld ('Optm1_geo', (/ 'lev' /), 'I', 'cm-3',' '   , gridname='geo_grid')

      call addfld ('MBAR_phys',(/ 'lev' /), 'I', ' ',' '   , gridname='physgrid')
      call addfld ('MBAR_geo', (/ 'lev' /), 'I', ' ',' '   , gridname='geo_grid')

      call addfld ('PED_phys',(/ 'lev' /), 'I', 'S/m','Pedersen Conductivity'  , gridname='physgrid')
      call addfld ('HAL_phys',(/ 'lev' /), 'I', 'S/m','Hall Conductivity'   , gridname='physgrid')

      !-------------------------------------------------------------------------------
      !  Set default values for ionosphere history variables
      !-------------------------------------------------------------------------------
      call phys_getopts(history_waccmx_out=history_waccmx)

      if (history_waccmx) then
         call add_default ('EDYN_ZIGM11_PED'         , 1, ' ')
         call add_default ('EDYN_ZIGM2_HAL'          , 1, ' ')
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

      character(len=300) :: dbgstr

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
      coord_map = (/ (i, i = mlat0, mlat1) /)
      lat_coord => horiz_coord_create('mlat', '', nmlat, 'latitude',           &
           'degrees_north', mlat0, mlat1, gmlat(mlat0:mlat1),                  &
           map=coord_map)
      nullify(coord_map)

      allocate(coord_map(omlon1 - mlon0 + 1))
      coord_map = (/ (i, i = mlon0, omlon1) /)
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
      coord_map = (/ (i, i = lat0, lat1) /)
      lat_coord => horiz_coord_create('glat', '', nlat, 'latitude',           &
           'degrees_north', lat0, lat1, glat(lat0:lat1),                  &
           map=coord_map)
      nullify(coord_map)

      allocate(coord_map(lon1 - lon0 + 1))
      coord_map = (/ (i, i = lon0, lon1) /)
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
