   module edyn_init
!
! Initialize edynamo 
!
  use shr_kind_mod   ,only: r8 => shr_kind_r8 ! 8-byte reals
  use shr_const_mod,  only: pi => shr_const_pi  
  use cam_logfile    ,only: iulog
  use cam_abortutils ,only: endrun
  use spmd_utils,     only: masterproc
  use infnan,         only: nan, assignment(=)

  use edyn_geogrid   ,only: nlon,nlat,nlev,nilev,glon,glat,zlev,zilev,&
                            nlonp1,nlonp2,nlatp1,jspole,jnpole,dlatg,dlong,&
                            ylatg,ylong,dphi,dlamda,cs,expz,zp
  use edyn_params    ,only: kbotdyn, pbotdyn

  implicit none

  private
  public :: edynamo_init, lonshift_global

  logical :: debug=.false. ! set true for prints to stdout at each call

  contains
!-----------------------------------------------------------------------
  subroutine edynamo_init(mpicomm, nlon_in,nlat_in,nlev_in, lonndx0,lonndx1,latndx0,latndx1,levndx0,levndx1, ntaski,ntaskj, &
                          glon_in, glat_in, pres_in, pres_edge_in )
!
! One-time initialization, called from inital.F90 after dyn_init and initcom.
!
    use edyn_maggrid  ,only: set_maggrid
    use edyn_mpi      ,only: mp_init,mp_distribute_geo,mp_distribute_mag,&
                             mp_exchange_tasks
#ifdef WACCMX_EDYN_ESMF
    use edynamo       ,only: alloc_edyn
    use edyn_esmf     ,only: edyn_esmf_init  ! initialize ESMF
#endif
!
! Args:
    integer, intent(in) :: mpicomm
    integer, intent(in) :: nlon_in,nlat_in,nlev_in
    integer, intent(in) :: lonndx0,lonndx1,latndx0,latndx1,levndx0,levndx1, ntaski,ntaskj
    real(r8),intent(in) :: glon_in(:), glat_in(:)
    real(r8),intent(in) :: pres_in(:), pres_edge_in(:)

    if (masterproc) then
       write(iulog,"('Enter edynamo_init:')")
    endif

    call mp_init(mpicomm) ! get ntask,mytid
    call set_geogrid(nlon_in,nlat_in,nlev_in, glon_in, glat_in, pres_in, pres_edge_in) ! set global geographic grid 
    call set_maggrid ()   ! set parameter-based global magnetic grid

    call mp_distribute_geo(lonndx0,lonndx1,latndx0,latndx1,levndx0,levndx1, ntaski,ntaskj)
    call mp_distribute_mag
    call register_maggrid
    call mp_exchange_tasks(0) ! single arg is iprint

#ifdef WACCMX_EDYN_ESMF
    call alloc_edyn      ! allocate dynamo arrays
    call edyn_esmf_init(mpicomm)       ! initialize ESMF
#endif

    call add_fields      ! add fields to WACCM history master list

  end subroutine edynamo_init
!-----------------------------------------------------------------------
  subroutine set_geogrid( nlon_in,nlat_in,nlev_in, glon_in, glat_in, pres_in, pres_edge_in )

  ! Args
    integer, intent(in) :: nlon_in,nlat_in,nlev_in
    real(r8),intent(in) :: glon_in(:), glat_in(:)
    real(r8),intent(in) :: pres_in(:), pres_edge_in(:)
!
! Local:
    integer :: i,j,js, k
    real(r8) :: real8,phi
    real(r8),parameter :: eps = 1.e-6_r8

    real(r8) :: pmid(nlev_in)

    nlon = nlon_in
    nlat = nlat_in
    nlev = nlev_in

    nilev = nlev+1

    nlonp1 = nlon+1
    nlonp2 = nlon+2
    nlatp1 = nlat+1

    jspole = 1
    jnpole = nlat

    dphi   = pi/dble(nlat)
    dlamda = 2._r8*pi/dble(nlon)

!
! Allocate coordinate variables:
!
    allocate(glon(nlon))
    glon(:nlon) = glon_in(:nlon)

    allocate(glat(nlat))
    glat(:nlat) = glat_in(:nlat)

    allocate(zlev(nlev))
    allocate(zilev(nilev))
!
! zp and expz are not set until oplus is called from dpie_coupling.
    allocate(zp(nlev))      ! log pressure (as in TIEGCM)
    allocate(expz(nlev))    ! exp(-zp)
    zp = nan
    expz = nan
!
!
    call lonshift_global(glon,nlon,'-180to180',.true.) ! shift to +/-180
!
! Hybrid-sigma levels from ref_pres module:
!
    zlev(:nlev)  = pres_in(:)  ! midpoints vertical coord (top down)
    zilev(:nilev) = pres_edge_in(:nilev)  ! interfaces vertical coord

    ! do bottom up search for kbotdyn
    pmid(:nlev) = zlev(nlev:1:-1)
    kloop: do k=1,nlev
       if ( pmid(k) <= pbotdyn) then
          kbotdyn = k
          exit kloop
       end if
    enddo kloop
    if ( kbotdyn < 1 ) then
       call endrun('set_geogrid: kbotdyn is not set')
    endif
    if (debug) then
       write(iulog,"('set_geogrid: kbotdyn=',i4,' pmid(kbotdyn)=',es12.4)") kbotdyn,pmid(kbotdyn)
    endif

!
! Set horizontal geographic grid in radians (for apex code):
!
    allocate(ylatg(nlat))   ! waccm grid includes poles
    allocate(ylong(nlonp1)) ! single periodic point
    real8 = dble(nlat) ; dlatg = pi/real8
    real8 = dble(nlon) ; dlong = 2._r8*pi/real8
    ylatg(1)    = -pi/2._r8+eps ! south pole
    ylatg(nlat) =  pi/2._r8-eps ! north pole
    do j=2,nlat-1
      real8 = dble(j-1)
      ylatg(j) = -0.5_r8*(pi-dlatg)+real8*dlatg
    enddo
    do i=1,nlonp1
      real8 = dble(i-1)
      ylong(i) = -pi+real8*dlong
    enddo
!
! Calculate cosine of latitude
!
    allocate(cs(0:nlat+1))
    js = -(nlat/2)
    do j=1,nlat
      phi = (j+js-.5_r8)*dphi
      cs(j) = cos(phi)
    enddo
    cs(0) = -cs(1)
    cs(nlat+1) = -cs(nlat)

  end subroutine set_geogrid
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
    call addfld ('ZMAG'      ,(/ 'lev' /), 'I', 'cm   ','ZMAG: Geopotential'               ,gridname='gmag_grid')
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
    call addfld ('amie_phihm' , horiz_only, 'I','VOLTS','AMIE Electric Potential-mag grid' ,gridname='gmag_grid')
    call addfld ('amie_efxm'  , horiz_only, 'I','mW/m2','AMIE energy flux on mag grid'     ,gridname='gmag_grid')
    call addfld ('amie_kevm'  , horiz_only, 'I','keV  ','AMIE mean energy on mag grid'     ,gridname='gmag_grid')
    call addfld ('amie_efxg'  , horiz_only, 'I','mW/m2','AMIE energy flux on geo grid'     ,gridname='fv_centers')
    call addfld ('amie_kevg'  , horiz_only, 'I','keV  ','AMIE mean energy on geo grid'     ,gridname='fv_centers')

!
! Dynamo inputs from sub dynamo_input (edynamo.F90):
    call addfld ('EDYN_TN   ',(/ 'lev' /), 'I', 'deg K   ','EDYN_TN'   , gridname='fv_centers') 
    call addfld ('EDYN_UN   ',(/ 'lev' /), 'I', 'cm/s    ','EDYN_UN'   , gridname='fv_centers')
    call addfld ('EDYN_VN   ',(/ 'lev' /), 'I', 'cm/s    ','EDYN_VN'   , gridname='fv_centers')
    call addfld ('EDYN_OMG  ',(/ 'lev' /), 'I', 's-1     ','EDYN_OMG'  , gridname='fv_centers')
    call addfld ('EDYN_Z    ',(/ 'lev' /), 'I', 'cm      ','EDYN_ZHT'  , gridname='fv_centers')
    call addfld ('EDYN_BARM ',(/ 'lev' /), 'I', '        ','EDYN_MBAR' , gridname='fv_centers')
    call addfld ('EDYN_PED  ',(/ 'lev' /), 'I', 'S/m     ','EDYN_PED'  , gridname='fv_centers')
    call addfld ('EDYN_HALL ',(/ 'lev' /), 'I', 'S/m     ','EDYN_HALL' , gridname='fv_centers')
 
!    call addfld ('EDYN_SCHT   ',(/ 'lev' /), 'I', '        ','EDYN_SCHT  ' , gridname='fv_centers')
    call addfld ('EDYN_WN     ',(/ 'lev' /), 'I', 'm/s     ','EDYN_WN    ' , gridname='fv_centers')
    call addfld ('EDYN_ADOTV1 ',(/ 'lev' /), 'I', '        ','EDYN_ADOTV1' , gridname='fv_centers')
    call addfld ('EDYN_ADOTV2 ',(/ 'lev' /), 'I', '        ','EDYN_ADOTV2' , gridname='fv_centers')
!
! 2d dynamo input fields on geo grid (edynamo.F90):
    call addfld ('EDYN_SINI   ', horiz_only , 'I', '        ','EDYN_SINI'   , gridname='fv_centers')
    call addfld ('EDYN_ADOTA1 ', horiz_only , 'I', '        ','EDYN_ADOTA1' , gridname='fv_centers')
    call addfld ('EDYN_ADOTA2 ', horiz_only , 'I', '        ','EDYN_ADOTA2' , gridname='fv_centers')
    call addfld ('EDYN_A1DTA2 ', horiz_only , 'I', '        ','EDYN_A1DTA2' , gridname='fv_centers')
    call addfld ('EDYN_BE3    ', horiz_only , 'I', '        ','EDYN_BE3'    , gridname='fv_centers')
 

    call addfld ('ADOTA1', horiz_only , 'I', ' ','ADOTA1' , gridname='fv_centers')
    call addfld ('ADOTA1_MAG', horiz_only , 'I', ' ','ADOTA1 in geo-mag coords' , gridname='fv_centers')

! 3d ion drifts and 2d conductances at end of dpie_coupling 
! (from either edynamo or time3d):
!
!       call addfld ('TIME3D_ZIGM11',horiz_only,'I',' ','TIME3D_ZIGM11',gridname='gmag_grid)
!       call addfld ('TIME3D_ZIGM22',horiz_only,'I',' ','TIME3D_ZIGM22',gridname='gmag_grid)
!       call addfld ('TIME3D_ZIGMC' ,horiz_only,'I',' ','TIME3D_ZIGMC' ,gridname='gmag_grid)
!       call addfld ('TIME3D_ZIGM2' ,horiz_only,'I',' ','TIME3D_ZIGM2' ,gridname='gmag_grid)
!       call addfld ('TIME3D_RIM1'  ,horiz_only,'I',' ','TIME3D_RIM1'  ,gridname='gmag_grid)
!       call addfld ('TIME3D_RIM2'  ,horiz_only,'I',' ','TIME3D_RIM2'  ,gridname='gmag_grid)

!       call addfld ('TIME3D_UI',(/ 'lev' /),'I',' ','TIME3D_UI')
!       call addfld ('TIME3D_VI',(/ 'lev' /),'I',' ','TIME3D_VI')
!       call addfld ('TIME3D_WI',(/ 'lev' /),'I',' ','TIME3D_WI')

!       call addfld ('T3D_OP_2WACCM',(/ 'lev' /),'I',' ','T3D_OP_2WACCM')
!       call addfld ('DPIE_OP',(/ 'lev' /),'I',' ','DPIE_OP') ! this is also below

    call addfld ('QEP',(/ 'lev' /), 'I', 'm^3/s'   ,'Photo-Electron Production', gridname='fv_centers')
    call addfld ('QOP',(/ 'lev' /), 'I', 'm^3/s'   ,'O+ Production Rate'       , gridname='fv_centers')
    call addfld ('OpO2',(/ 'lev' /), 'I', 'cm^3/s' ,'Op+O2 Loss Rate'          , gridname='fv_centers')
    call addfld ('OpN2',(/ 'lev' /), 'I', 'cm^3/s' ,'Op+N2 Loss Rate'          , gridname='fv_centers')
    call addfld ('LOP',(/ 'lev' /), 'I', 'cm^3/s'  ,'O+ Loss Rate'             , gridname='fv_centers')
    call addfld ('SIGMA_PED' ,(/ 'lev' /), 'I', ' ','Pederson Conductivity'    , gridname='fv_centers')
    call addfld ('SIGMA_HALL',(/ 'lev' /), 'I', ' ','Hall Conductivity'        , gridname='fv_centers')

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

    call addfld ('EDYN_UI',(/ 'lev' /), 'I', 'cm/s','EDYN_UI', gridname='fv_centers')
    call addfld ('EDYN_VI',(/ 'lev' /), 'I', 'cm/s','EDYN_VI', gridname='fv_centers')
    call addfld ('EDYN_WI',(/ 'lev' /), 'I', 'cm/s','EDYN_WI', gridname='fv_centers')

    call addfld ('POTEN'  ,(/ 'lev' /), 'I', 'Volts','POTEN: Electric Potential',&
         gridname='fv_centers')
    call addfld ('EX'     ,(/ 'lev' /), 'I', 'V/m'  ,'EX: Zonal component of Electric Field',&
         gridname='fv_centers')
    call addfld ('EY'     ,(/ 'lev' /), 'I', 'V/m'  ,'EY: Meridional component of Electric Field',&
         gridname='fv_centers')
    call addfld ('EZ'     ,(/ 'lev' /), 'I', 'V/m'  ,'EZ: Vertical component of Electric Field',&
         gridname='fv_centers')

    call addfld ('ZEDYN360   ' ,(/ 'lev' /), 'I', 'm     ','Geopotential 0 to 360 lon grid', gridname='fv_centers')
    call addfld ('ZEDYN180   ',(/ 'lev' /), 'I', 'm     ','Geopotential -180 to 180 lon grid', gridname='fv_centers')

    call addfld ('BMOD'       , horiz_only, 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('XB'       ,   horiz_only, 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('YB'       ,   horiz_only, 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('ZB'       ,   horiz_only, 'I', '  ',' '  ,gridname='fv_centers')

    call addfld ('RJAC11'       ,(/'lev'/), 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('RJAC12'       ,(/'lev'/), 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('RJAC21'       ,(/'lev'/), 'I', '  ',' '  ,gridname='fv_centers')
    call addfld ('RJAC22'       ,(/'lev'/), 'I', '  ',' '  ,gridname='fv_centers')

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

  subroutine register_maggrid

  use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register
  use edyn_mpi,         only: mlat0,mlat1,mlon0,omlon1
  use edyn_maggrid,     only: gmlat, gmlon, nmlat, nmlon
  integer, parameter :: mag_decomp         = 111 !arbitrary value

  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  integer(iMap),       pointer :: grid_map(:,:)
  integer(iMap),       pointer :: coord_map(:)
  integer                      :: i,j,ind

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
  lat_coord => horiz_coord_create('mlat', '', nmlat, 'latitude',                &
       'degrees_north', mlat0, mlat1, gmlat(mlat0:mlat1),                      &
       map=coord_map)
  nullify(coord_map)

    allocate(coord_map(omlon1 - mlon0 + 1))
    coord_map = (/ (i, i = mlon0, omlon1) /)
    lon_coord => horiz_coord_create('mlon', '', nmlon, 'longitude',             &
         'degrees_east', mlon0, omlon1, gmlon(mlon0:omlon1),                   &
         map=coord_map)
    deallocate(coord_map)
    nullify(coord_map)

  call cam_grid_register('gmag_grid', mag_decomp, lat_coord, lon_coord,      &
       grid_map, unstruct=.false.)
  nullify(grid_map)

  end subroutine register_maggrid

!-----------------------------------------------------------------------
end module edyn_init
