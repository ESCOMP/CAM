module edyn_esmf
#ifdef WACCMX_EDYN_ESMF

  use esmf           ,only: ESMF_Grid, ESMF_Field, ESMF_RouteHandle, & ! ESMF library module
                            ESMF_SUCCESS, ESMF_KIND_R8, ESMF_KIND_I4, &
		            ESMF_FieldGet, ESMF_STAGGERLOC_CENTER, ESMF_FieldRegridStore, &
		            ESMF_REGRIDMETHOD_BILINEAR, ESMF_POLEMETHOD_ALLAVG, ESMF_FieldSMMStore, &
                            ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL, ESMF_GridAddCoord, ESMF_GridGetCoord, &
                            ESMF_TYPEKIND_R8, ESMF_FieldCreate, ESMF_Array, ESMF_ArraySpec, ESMF_DistGrid, &
                            ESMF_GridGet, ESMF_ArraySpecSet, ESMF_ArrayCreate, ESMF_FieldGet, ESMF_FieldSMM, &
                            ESMF_TERMORDER_SRCSEQ
  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_logfile    ,only: iulog
  use cam_abortutils ,only: endrun
  use edyn_mpi       ,only: ntask,ntaski,ntaskj,tasks,lon0,lon1,lat0,lat1,&
                            nmagtaski,nmagtaskj,mlon0,mlon1,mlat0,mlat1
  use getapex        ,only: gdlatdeg,gdlondeg
  use edyn_geogrid   ,only: nlon,nlat,nlev,glon,glat,jspole,jnpole  ! dynamically allocated geo grid
  use edyn_maggrid   ,only: nmlev,gmlat,gmlon

#endif

  implicit none
  save
  private

  public :: edyn_esmf_update

#ifdef WACCMX_EDYN_ESMF

  public :: nf_3dgeo,f_3dgeo
  public :: edyn_esmf_update_flag
  public :: edyn_esmf_init, edyn_esmf_final, edyn_esmf_update_step, edyn_esmf_regrid
  public :: edyn_esmf_get_2dfield, edyn_esmf_set2d_geo, edyn_esmf_get_3dfield, edyn_esmf_set3d_mag, edyn_esmf_set3d_geo 
  public :: edyn_esmf_set2d_mag
  
  public :: mag_be3, mag_adota1,mag_adota2,mag_a1dta2,mag_sini,mag_adotv2,mag_adotv1,mag_scht
  public :: mag_efx, mag_kev
  public :: mag_zpot,mag_hal,mag_ped, mag_phi3d
  public :: geo_be3,geo_adotv2,geo_a1dta2,geo_adota2,geo_adota1,geo_adotv1,geo_sini,geo_scht,geo_zpot
  public :: geo_efx, geo_kev
  public :: geo_hal, geo_ped, mag_des_grid, geo_src_grid, geo_phi3d, geo_emz3d, geo_elam3d, geo_ephi3d
  public :: mag_emz3d, mag_elam3d, mag_ephi3d

  type(ESMF_Grid) :: &
    geo_src_grid, mag_src_grid, & ! source grids (will not have periodic pts)
    geo_des_grid, mag_des_grid    ! destination grids (will have periodic pts)
!
! 3d (i,j,k) ESMF Fields on geographic subdomains:
!
  type(ESMF_Field) :: & ! 3d ESMF fields on geographic grid
    geo_ped,        & ! pedersen conductivity
    geo_hal,        & ! hall conductivity
    geo_zpot,       & ! geopotential height (cm)
    geo_scht,       & ! scale height (cm)
    geo_adotv1,     & ! ue1 (m/s)
    geo_adotv2        ! ue2 (m/s)
  integer,parameter :: nf_3dgeo=6            ! number of 3d fields on geographic grid
  type(ESMF_Field) :: f_3dgeo(nf_3dgeo) ! fields on 3d geo grid (could be bundled?)
!
! 2d (i,j) ESMF fields on geographic subdomains:
!
  type(ESMF_Field) :: & ! 2d ESMF fields on geographic grid
    geo_sini,       & ! sin(I_m)
    geo_adota1,     & ! d(1)**2/D
    geo_adota2,     & ! d(2)**2/D
    geo_a1dta2,     & ! (d(1) dot d(2)) /D
    geo_be3,        & ! mag field strength (T)
    geo_efx, geo_kev  ! amie fields
!
! 3d (i,j,k) ESMF fields regridded to magnetic subdomains:
!
  type(ESMF_Field) :: & ! 3d ESMF fields on geomagnetic grid
    mag_ped,        & ! pedersen conductivity
    mag_hal,        & ! hall conductivity
    mag_zpot,       & ! geopotential height (cm)
    mag_scht,       & ! scale height (cm)
    mag_adotv1,     & ! ue1 (m/s)
    mag_adotv2        ! ue2 (m/s)
!
! 2d (i,j) ESMF fields on magnetic subdomains:
!
  type(ESMF_Field) :: & ! 2d fields on geomagnetic grid
    mag_sini,       & ! sin(I_m)
    mag_adota1,     & ! d(1)**2/D
    mag_adota2,     & ! d(2)**2/D
    mag_a1dta2,     & ! (d(1) dot d(2)) /D
    mag_be3,        & ! mag field strength (T)
    mag_efx, mag_kev  ! amie fields
!
! 3d electric potential and electric field for mag to geo regridding:
!
  type(ESMF_Field) :: mag_phi3d,mag_ephi3d,mag_elam3d,mag_emz3d
  type(ESMF_Field) :: geo_phi3d,geo_ephi3d,geo_elam3d,geo_emz3d

  type(ESMF_RouteHandle) :: & ! ESMF route handles for regridding
    routehandle_geo2mag,    & ! for geo to mag regrid 
    routehandle_mag2geo,    & ! for mag to geo regrid
    routehandle_geo2mag_2d, & ! for 2d geo to mag
    routehandle_mag2geo_2d    ! for 2d mag to geo for AMIE fields
!
  real(r8),allocatable :: unitv(:)
!
  private routehandle_geo2mag, routehandle_mag2geo,&
    routehandle_geo2mag_2d

  logical, protected :: edyn_esmf_update_step = .true.
  logical :: debug=.false. ! set true for prints to stdout at each call
#endif

  contains
#ifdef WACCMX_EDYN_ESMF
!-----------------------------------------------------------------------
  subroutine edyn_esmf_init( mpi_comm )

    integer, intent(in) :: mpi_comm

  end subroutine edyn_esmf_init

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine edyn_esmf_final

  end subroutine edyn_esmf_final

#endif

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine edyn_esmf_update
    use getapex, only: get_apex,magfield, alonm
    use mo_apex, only: geomag_year_updated

#ifdef WACCMX_EDYN_ESMF
! Create ESMF grids for geographic and magnetic, and create ESMF fields
!   as necessary on both grids. Define the 2d coordinates for each grid,
!   and save an ESMF routehandles for geo2mag and mag2geo regridding.
!
! Local:
    integer :: rc ! return code for ESMF calls
    real(ESMF_KIND_R8),pointer :: fptr(:,:,:)
    integer :: lbnd_destgeo(3),ubnd_destgeo(3) ! 3d bounds of destination geo grid
    integer :: lbnd_destmag(3),ubnd_destmag(3) ! 3d bounds of destination mag grid
    integer :: lbnd_srcgeo(3),ubnd_srcgeo(3)   ! 3d bounds of source geo grid
    integer :: lbnd_srcmag(3),ubnd_srcmag(3)   ! 3d bounds of source mag grid
    integer(ESMF_KIND_I4),pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),pointer :: factorList(:)
    integer :: smm_srctermproc,  smm_pipelinedep
#endif

    if (.not.geomag_year_updated .and. allocated(alonm)) return
!
! Get apex coordinates.
!
    call get_apex( )     ! get apex coordinates
    call magfield        ! calculate magnetic field parameters

#ifdef WACCMX_EDYN_ESMF

    smm_srctermproc = 0
    smm_pipelinedep = 16
!
! Set unit vector (this routine called once per run unless crossing year boundary):
! Handle year boundary by checking if field is allocated
!
    if (.not.allocated(unitv)) allocate(unitv(nlon))
    unitv(:) = 1._r8    
!
! Make magnetic and geographic grids for geo2mag regridding:
!
    call create_geo_grid(geo_src_grid,'src')  ! geo source grid
    call create_mag_grid(mag_des_grid,'des')  ! mag destination grid
! 
! Make grids for mag2geo regridding:
!
    call create_mag_grid(mag_src_grid,'src')
    call create_geo_grid(geo_des_grid,'des')
! 
! Create empty fields on geographic grid that will be transformed to
!   the magnetic grid and passed as input to the dynamo. This does not
!   assign any values.
!
! 3d fields on source geo grid (these exclude periodic points):
!
    call edyn_esmf_create_geofield(geo_ped,geo_src_grid,   'PED     ',nlev)
    call edyn_esmf_create_geofield(geo_hal ,geo_src_grid,  'HAL     ',nlev)
    call edyn_esmf_create_geofield(geo_zpot,geo_src_grid,  'ZPOT    ',nlev)
    call edyn_esmf_create_geofield(geo_scht,geo_src_grid,  'SCHT    ',nlev)
    call edyn_esmf_create_geofield(geo_adotv1,geo_src_grid,'ADOTV1  ',nlev)
    call edyn_esmf_create_geofield(geo_adotv2,geo_src_grid,'ADOTV2  ',nlev)
!
! Get 3d bounds of source geo field:
!
    call ESMF_FieldGet(geo_ped,localDe=0,farrayPtr=fptr, &
      computationalLBound=lbnd_srcgeo,                   &
      computationalUBound=ubnd_srcgeo,rc=rc)

    if (debug) then
      write(iulog,"('Bounds of source geo field: lbnd_srcgeo=',3i4,' ubnd_srcgeo=',3i4,' glon=',2f9.3)") &
        lbnd_srcgeo,ubnd_srcgeo
    endif
!
! 2d fields on source geo grid (these exclude periodic points):
!
    call edyn_esmf_create_geofield(geo_sini  ,geo_src_grid,'SINI    ',0)
    call edyn_esmf_create_geofield(geo_adota1,geo_src_grid,'ADOTA1  ',0)
    call edyn_esmf_create_geofield(geo_adota2,geo_src_grid,'ADOTA2  ',0)
    call edyn_esmf_create_geofield(geo_a1dta2,geo_src_grid,'A1DTA2  ',0) 
    call edyn_esmf_create_geofield(geo_be3   ,geo_src_grid,'BE3     ',0) 
!
! 3d fields on destination mag grid (will include periodic point):
!
    call edyn_esmf_create_magfield(mag_ped ,mag_des_grid,  'PED     ',nmlev)
    call edyn_esmf_create_magfield(mag_hal ,mag_des_grid,  'HAL     ',nmlev)
    call edyn_esmf_create_magfield(mag_zpot,mag_des_grid,  'ZPOT    ',nmlev)
    call edyn_esmf_create_magfield(mag_scht,mag_des_grid,  'SCHT    ',nmlev)
    call edyn_esmf_create_magfield(mag_adotv1,mag_des_grid,'ADOTV1  ',nmlev)
    call edyn_esmf_create_magfield(mag_adotv2,mag_des_grid,'ADOTV2  ',nmlev)
!
! Get 3d bounds of destination mag field:
!
    call ESMF_FieldGet(mag_ped,localDe=0,farrayPtr=fptr, &
      computationalLBound=lbnd_destmag,                  &
      computationalUBound=ubnd_destmag,rc=rc)

    if (debug) then
      write(iulog,"('Bounds of destination mag field: lbnd_destmag=',3i4,' ubnd_destmag=',3i4,' gmlon=',2f9.3)") &
        lbnd_destmag,ubnd_destmag
      write(iulog,"('esmf_init: lon bnd_destmag =',2i4,' gmlon=',2f9.3)") &
        lbnd_destmag(1),ubnd_destmag(1),gmlon(lbnd_destmag(1)),gmlon(ubnd_destmag(1))
      write(iulog,"('esmf_init: lat bnd_destmag =',2i4,' gmlat=',2f9.3)") &
        lbnd_destmag(2),ubnd_destmag(2),gmlat(lbnd_destmag(2)),gmlat(ubnd_destmag(2))
    endif
!
! 2d fields on destination mag grid (will include periodic point):
!
    call edyn_esmf_create_magfield(mag_sini  ,mag_des_grid,'SINI    ',0)
    call edyn_esmf_create_magfield(mag_adota1,mag_des_grid,'ADOTA1  ',0)
    call edyn_esmf_create_magfield(mag_adota2,mag_des_grid,'ADOTA2  ',0)
    call edyn_esmf_create_magfield(mag_a1dta2,mag_des_grid,'A1DTA2  ',0)
    call edyn_esmf_create_magfield(mag_be3   ,mag_des_grid,'BE3     ',0)
!
! 3d fields on source mag grid for mag2geo:
!
    call edyn_esmf_create_magfield(mag_phi3d ,mag_src_grid,'PHIM3D  ',nmlev)
    call edyn_esmf_create_magfield(mag_ephi3d,mag_src_grid,'EPHI3D  ',nmlev)
    call edyn_esmf_create_magfield(mag_elam3d,mag_src_grid,'ELAM3D  ',nmlev)
    call edyn_esmf_create_magfield(mag_emz3d ,mag_src_grid,'EMZ3D   ',nmlev)
    call edyn_esmf_create_magfield(mag_efx   ,mag_src_grid,'MEFXAMIE',0)
    call edyn_esmf_create_magfield(mag_kev   ,mag_src_grid,'MKEVAMIE',0)
!
! 3d fields on destination geo grid for mag2geo:
!
    call edyn_esmf_create_geofield(geo_phi3d ,geo_des_grid,'PHIG3D  ',nlev)
    call edyn_esmf_create_geofield(geo_ephi3d,geo_des_grid,'EPHI3D  ',nlev)
    call edyn_esmf_create_geofield(geo_elam3d,geo_des_grid,'ELAM3D  ',nlev)
    call edyn_esmf_create_geofield(geo_emz3d ,geo_des_grid,'EMZ3D   ',nlev)
    call edyn_esmf_create_geofield(geo_efx   ,geo_des_grid,'GEFXAMIE',0)
    call edyn_esmf_create_geofield(geo_kev   ,geo_des_grid,'GKEVAMIE',0)
!
! Get 3d bounds of source mag field:
    call ESMF_FieldGet(mag_phi3d,localDe=0,farrayPtr=fptr,&
      computationalLBound=lbnd_srcmag,                    &
      computationalUBound=ubnd_srcmag,rc=rc)

    if (debug) then
      write(iulog,"('esmf_init: lon bnd_srcmag =',2i4,' gmlon=',2f9.3)") &
        lbnd_srcmag(1),ubnd_srcmag(1)
      write(iulog,"('esmf_init: lat bnd_srcmag =',2i4,' gmlat=',2f9.3)") &
        lbnd_srcmag(2),ubnd_srcmag(2)
    endif
!
! Get 3d bounds of destination geo field:
!
    call ESMF_FieldGet(geo_phi3d,localDe=0,farrayPtr=fptr,&
      computationalLBound=lbnd_destgeo,                   &
      computationalUBound=ubnd_destgeo,rc=rc)

    if (debug) then
      write(iulog,"('esmf_init: lon bnd_destgeo=',2i4,' glon=',2f9.3)") &
        lbnd_destgeo(1),ubnd_destgeo(1)
      write(iulog,"('esmf_init: lat bnd_destgeo=',2i4,' glat=',2f9.3)") &
        lbnd_destgeo(2),ubnd_destgeo(2)
    endif
!
! Save route handles for grid transformations in both directions
! geo2mag and mag2geo. FieldRegridStore needs to be called only 
! once for each transformation before the timestep loop (src and 
! dest fields are still required, so just use ped here). Once inside
! the timestep loop, the same routehandle can be used for all fields
! that are regridded in the given direction.
!
! These calls will leave *.vtk info files in execdir:
!   call ESMF_GridWriteVTK(geo_src_grid, &
!     staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoGrid",rc=rc)
!   call ESMF_GridWriteVTK(mag_des_grid, &
!     staggerloc=ESMF_STAGGERLOC_CENTER, filename="magGrid",rc=rc)
!
! Save route handle and get esmf indices and weights for geo2mag:
!
    call ESMF_FieldRegridStore(srcField=geo_ped,dstField=mag_ped,       &
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                          &
      polemethod=ESMF_POLEMETHOD_ALLAVG,                                &
      routeHandle=routehandle_geo2mag,factorIndexList=factorIndexList,  &
      factorList=factorList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(a,a,i4)") '>>> edyn_esmf_update: error return from ', &
        'ESMF_FieldRegridStore for 3d geo2mag: rc=',rc
      call endrun('edyn_esmf_update: ESMF_FieldRegridStore ped')
    endif
!
! Store route handle for geo2mag 3d fields.
!
    call ESMF_FieldSMMStore(geo_ped,mag_ped,routehandle_geo2mag, &
      factorList,factorIndexList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> edyn_esmf_update: error return from ESMF_FieldSMMStore for ',&
        '3d geo2mag: rc=',rc
      call endrun('edyn_esmf_update: ESMF_FieldSMMStore for 3d geo2mag ped')
    endif
!    
! Store route handle geo2mag 2d fields:
!    
    call ESMF_FieldSMMStore(geo_sini,mag_sini,routehandle_geo2mag_2d, &
      factorList,factorIndexList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> edyn_esmf_update: error return from ESMF_FieldSMMStore',&
        ' for 2d geo2mag: rc=',rc
      call endrun('edyn_esmf_update: ESMF_FieldSMMStore for 2d geo2mag sini')
    endif
!
! Save route handle and get esmf indices and weights for mag2geo:
! (this overwrites factorIndexList and factorList from geo2mag call above)
!
! These calls will leave *.vtk info files in execdir:
!   call ESMF_GridWriteVTK(mag_src_grid, &
!     staggerloc=ESMF_STAGGERLOC_CENTER, filename="magSrcGrid",rc=rc)
!   call ESMF_GridWriteVTK(geo_des_grid, &
!     staggerloc=ESMF_STAGGERLOC_CENTER, filename="geoDesGrid",rc=rc)

! Save route handle and get esmf indices and weights for mag2geo:
!
    call ESMF_FieldRegridStore(srcField=mag_phi3d,dstField=geo_phi3d, &
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                        &
      polemethod=ESMF_POLEMETHOD_ALLAVG,                              &
      routeHandle=routehandle_mag2geo,factorIndexList=factorIndexList,&
      factorList=factorList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> edyn_esmf_update: error return from ',&
        'ESMF_FieldRegridStore for 3d mag2geo: rc=',rc
      call endrun('edyn_esmf_update: ESMF_FieldRegridStore for 3d mag2geo phi3d')
    endif
!
! mag2geo 3d fields:
!
    call ESMF_FieldSMMStore(mag_phi3d,geo_phi3d,routehandle_mag2geo,&
      factorList,factorIndexList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> edyn_esmf_update: error return from ESMF_FieldSMMStore ',&
        'for 3d mag2geo: rc=',rc
      call endrun('edyn_esmf_update: ESMF_FieldSMMStore for 3d geo2mag phi3d')
    endif

! amie fields
    call ESMF_FieldRegridStore(srcField=mag_efx,dstField=geo_efx, &
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                        &
      polemethod=ESMF_POLEMETHOD_ALLAVG,                              &
      routeHandle=routehandle_mag2geo_2d,factorIndexList=factorIndexList,&
      factorList=factorList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(6,"(2a,i4)") '>>> esmf_init: error return from ',&
        'ESMF_FieldRegridStore for 2d mag2geo_2d: rc=',rc
      call endrun
    endif
!
! mag2geo 2d fields:
!
    call ESMF_FieldSMMStore(mag_efx,geo_efx,routehandle_mag2geo_2d,&
      factorList,factorIndexList,srcTermProcessing=smm_srctermproc,pipelineDepth=smm_pipelinedep,rc=rc)


    edyn_esmf_update_step = .true.
#endif
  end subroutine edyn_esmf_update

#ifdef WACCMX_EDYN_ESMF
!-----------------------------------------------------------------------
  real(r8) function select_wt_mag2geo(n,dimx,djmy)
    integer,intent(in) :: n
    real(r8),intent(in) :: dimx,djmy

    select_wt_mag2geo = 0._r8
    select case (n)
      case(1)
        select_wt_mag2geo = (1._r8-dimx)*(1._r8-djmy)
      case(2)
        select_wt_mag2geo = dimx*(1._r8-djmy)
      case(3)
        select_wt_mag2geo = dimx*djmy
      case(4)
        select_wt_mag2geo = (1._r8-dimx)*djmy
    end select
  end function select_wt_mag2geo
!-----------------------------------------------------------------------
  subroutine create_mag_grid(grid_out,srcdes)
!
! Create ESMF geomagnetic grid, w/ lon,lat coordinates.
! This is called from esmf_init during model initialization.
!
! Args:
    type(ESMF_Grid),intent(out) :: grid_out
    character(len=*),intent(in) :: srcdes
!
! Local:
    integer :: i,j,n,rc
    real(ESMF_KIND_R8),pointer :: coordX(:,:),coordY(:,:)
    integer :: lbnd(2),ubnd(2)
    integer :: nmlons_task(ntaski) ! number of lons per task 
    integer :: nmlats_task(ntaskj) ! number of lats per task
!
! We are creating either a source grid or a destination grid:
!
    if (srcdes /= 'src' .and. srcdes /= 'des') then
      write(iulog,"(a)") '>>> create_mag_grid: srcdes = ''',srcdes, &
        ''' but must be either ''src'' or ''des'''
      call endrun('create_mag_grid: srcdes')
    endif
!
! nmlons_task(nmagtaski) = number of mag lons per task in lon dim
!
    do i=1,nmagtaski
      loop: do n=0,ntask-1
        if (tasks(n)%magtidi==i-1) then
          nmlons_task(i) = tasks(n)%nmaglons
          exit loop
        endif
      enddo loop
    enddo
!
! Exclude periodic points (1 point fewer for mpi tasks at east end)
! for source grids (this overwrites above for eastern-most tasks):
!
    if (srcdes == 'src') then
      do n=0,ntask-1
        if (tasks(n)%magtidi==nmagtaski-1) then  ! east edge of proc matrix
          nmlons_task(tasks(n)%magtidi+1) = tasks(n)%nmaglons-1
        endif
      enddo
    endif
!
! nmlats_task(nmagtaskj) = number of mag lats per task in lat dim
!
    do j=1,nmagtaskj
      loop1: do n=0,ntask-1
        if (tasks(n)%magtidj==j-1) then
          nmlats_task(j) = tasks(n)%nmaglats
          exit loop1
        endif
      enddo loop1
    enddo
!
! Create curvilinear magnetic grid (both coords depend
! on both dimensions, i.e., lon(i,j),lat(i,j)):
!
    grid_out = ESMF_GridCreate1PeriDim(               &
      countsPerDEDim1=nmlons_task, coordDep1=(/1,2/), &
      countsPerDEDim2=nmlats_task, coordDep2=(/1,2/), &
      indexflag=ESMF_INDEX_GLOBAL,rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> create_mag_grid: error return from ',&
        'ESMF_GridCreateShapeTile: rc=',rc
      call endrun('create_mag_grid: ESMF_GridCreate1PeriDim')
    endif
!
! Allocate coordinates:
!
    call ESMF_GridAddCoord(grid_out,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)

    if (rc /=ESMF_SUCCESS) then
      write(iulog,"(2a,i4)") '>>> create_mag_grid: error return from ',&
        'ESMF_GridAddCoord: rc=',rc
      call endrun('create_mag_grid: ESMF_GridAddCoord mag_grid')
    endif
!
! Get pointer and set mag grid longitude coordinates:
!
    call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0,  &
      computationalLBound=lbnd, computationalUBound=ubnd,    &
      farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(i4)") '>>> create_mag_grid: error return from ',  &
        'ESMF_GridGetCoord for longitude coords: rc=',rc
      call endrun('create_mag_grid: ESMF_GridGetCoord mag grid longitude')
    endif

    do j=lbnd(2),ubnd(2)
      do i=lbnd(1),ubnd(1)
        coordX(i,j) = gdlondeg(i,j)
      enddo
    enddo
! 
! Get pointer and set mag grid latitude coordinates:
!
    call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0, &
      computationalLBound=lbnd, computationalUBound=ubnd,   &
      farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(i4)") '>>> create_mag_grid: error return from ',&
        'ESMF_GridGetCoord for latitude coords: rc=',rc
      call endrun('create_mag_grid: ESMF_GridGetCoord latitude')
    endif

    do j=lbnd(2),ubnd(2)
      do i=lbnd(1),ubnd(1)
        coordY(i,j) = gdlatdeg(i,j)
      enddo
    enddo

    if (debug) then
      write(iulog,"(4a,2i4,a,2i4,a,2i4,a,2i4)") 'Created ESMF ',srcdes,' mag grid: ',  &
        ' lbnd,ubnd_lon=',lbnd(1),ubnd(1),' mlon0,1=',mlon0,mlon1, &
        ' lbnd,ubnd_lat=',lbnd(2),ubnd(2),' mlat0,1=',mlat0,mlat1
    endif

  end subroutine create_mag_grid
!-----------------------------------------------------------------------
  subroutine create_geo_grid(grid_out,srcdes)
!
! Args:
    type(ESMF_Grid),intent(out) :: grid_out
    character(len=*),intent(in) :: srcdes
!
! Local:
    integer :: i,j,n,rc
    integer :: lbnd_lat,ubnd_lat,lbnd_lon,ubnd_lon,lbnd(1),ubnd(1)
    real(ESMF_KIND_R8),pointer :: coordX(:),coordY(:)
    integer :: nlons_task(ntaski) ! number of lons per task 
    integer :: nlats_task(ntaskj) ! number of lats per task
    logical :: has_poles
!
! We are creating either a source grid or a destination grid:
!
    if (srcdes /= 'src' .and. srcdes /= 'des') then
      write(iulog,"(a)") '>>> create_geo_grid: srcdes = ''',srcdes, &
        ''' but must be either ''src'' or ''des'''
      call endrun('create_geo_grid: srcdes')
    endif
!
! nlons_task(ntaski) = number of geo lons per task.
!
    do i=1,ntaski
      loop: do n=0,ntask-1
        if (tasks(n)%mytidi==i-1) then
          nlons_task(i) = tasks(n)%nlons
          exit loop
        endif
      enddo loop
    enddo
! 
! Exclude periodic points (2 points fewer for procs at each end)
! for source grids only (east and west edges of task table).
! (TIMEGCM only)
!
!   if (srcdes == 'src'.and.trim(model_name)=='TIMEGCM') then
!     do n=0,ntask-1
! east or west edge of task table:
!       if (tasks(n)%mytidi==ntaski-1.or.tasks(n)%mytidi==0) &
!         nlons_task(tasks(n)%mytidi+1) = tasks(n)%nlons-2
!     enddo
!   endif
!
! nlats_task(ntaskj) = number of geo lats per task.
!
    do j=1,ntaskj
      loop1: do n=0,ntask-1
        if (tasks(n)%mytidj==j-1) then
          nlats_task(j) = tasks(n)%nlats
          exit loop1
        endif
      enddo loop1
    enddo
!
! Check to see if global glat(nlat) has poles (WACCM does, TIMEGCM does not):
    has_poles = .false.
    do j=1,nlat
      if (abs(glat(j))==90._r8) has_poles = .true. 
    enddo

    if (debug) write(iulog,"('create_geo_grid: srcdes=',a,' has_poles=',l1)") srcdes,has_poles
!
! If making destination grid and glat does not have poles, add extra points 
! at north and south edges of task table:
!
    if (.not.has_poles.and.srcdes=='des') then ! probably TIMEGCM
      do n=0,ntask-1
! north or south edge of task table: add 1 lat for pole
        if (tasks(n)%mytidj==ntaskj-1.or.tasks(n)%mytidj==0) &
          nlats_task(tasks(n)%mytidj+1) = tasks(n)%nlats+1
      enddo
!
! Create 2d geographic destination grid (minimum lat index is 0 to include poles):
      grid_out = ESMF_GridCreate1PeriDim(            &
        countsPerDEDim1=nlons_task, coordDep1=(/1/), &
        countsPerDEDim2=nlats_task, coordDep2=(/2/), &
        indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,0/),rc=rc)

    elseif (has_poles) then  ! geo source grid does not have poles
!
! Create 2d geographic source grid (without poles)
      grid_out = ESMF_GridCreate1PeriDim(            &
        countsPerDEDim1=nlons_task, coordDep1=(/1/), &
        countsPerDEDim2=nlats_task, coordDep2=(/2/), &
        indexflag=ESMF_INDEX_GLOBAL,minIndex=(/1,1/),rc=rc)
    else
      write(iulog,*) 'No capability for ESMF to handle source grid without poles'
      call endrun('create_geo_grid: No ESMF capability for source grid without poles')    
    endif

    if (rc /=ESMF_SUCCESS) then
      write(iulog,"(/,2a,i4)") '>>> create_geo_grid: error return from ',&
        'ESMF_GridCreate1PeriDim: rc=',rc
      call endrun('create_geo_grid: ESMF_GridCreate1PeriDim')
    endif
!
! Allocate coordinates:
!
    call ESMF_GridAddCoord(grid_out,staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)

    if (rc /=ESMF_SUCCESS) then
      write(iulog,"(/,a)") '>>> create_geo_grid: error return from ESMF_GridAddCoord'
      call endrun('create_geo_grid: ESMF_GridAddCoord')
    endif
!
! Get pointer and set geo grid longitude coordinates:
!
    call ESMF_GridGetCoord(grid_out, coordDim=1, localDE=0, &
      computationalLBound=lbnd, computationalUBound=ubnd,   &
      farrayPtr=coordX, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if (rc /=ESMF_SUCCESS) then
      write(iulog,"(/,2a)") '>>> create_geo_grid: error return from ',&
        'ESMF_GridGetCoord for longitude coords'
      call endrun('create_geo_grid: ESMF_GridGetCoord longitude')
    endif
!
! Note glon was shifted to +/-180 by sub set_geogrid (edyn_init.F90)
!
    lbnd_lon = lbnd(1) ; ubnd_lon = ubnd(1)
    do i=lbnd_lon,ubnd_lon
       coordX(i) = glon(i)          ! 1 -> 72
    enddo
!
! Get pointer and set geo grid latitude coordinates, including poles:
!
    call ESMF_GridGetCoord(grid_out, coordDim=2, localDE=0, &
      computationalLBound=lbnd, computationalUBound=ubnd,   &
      farrayPtr=coordY, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)

    if (rc /=ESMF_SUCCESS) then
      write(iulog,"(/,2a)") '>>> create_geo_grid: error return from ',&
        'ESMF_GridGetCoord for latitude coords'
      call endrun('create_geo_grid: ESMF_GridGetCoord latitude')
    endif

    lbnd_lat = lbnd(1) ; ubnd_lat = ubnd(1)

    if (.not.has_poles.and.srcdes=='des') then ! geo destination grid has poles
      do j=lbnd_lat,ubnd_lat
        if (j==jspole) then
          coordY(j) = -90._r8
        elseif (j==jnpole) then
          coordY(j) = +90._r8
        else
          coordY(j) = glat(j)
        endif
      enddo
    elseif (has_poles) then
      do j=lbnd_lat,ubnd_lat
        coordY(j) = glat(j)
      enddo
    else
      write(iulog,*) 'No capability for ESMF to handle source grid without poles'
      call endrun('create_geo_grid: No ESMF capability for source grid without poles')    
    endif

    if (debug) then
      write(iulog,"(4a,2i4,a,2i4,a,2i4,a,2i4)") 'Created ESMF ',srcdes,' geo grid: ',  &
        ' lbnd,ubnd_lon=',lbnd_lon,ubnd_lon,' lon0,1=',lon0,lon1, &
        ' lbnd,ubnd_lat=',lbnd_lat,ubnd_lat,' lat0,1=',lat0,lat1

      write(iulog,"('coordX for ',a,' geo grid = ',/,(8f10.4))") srcdes,coordX
      write(iulog,"('coordY for ',a,' geo grid = ',/,(8f10.4))") srcdes,coordY
    endif

  end subroutine create_geo_grid
!-----------------------------------------------------------------------
  subroutine edyn_esmf_create_geofield(field,grid,name,nlev)
!
! Create ESMF field (2d or 3d) on geo grid (will exclude periodic points)
! If nlev == 0, field is 2d (i,j), otherwise field is 3d,
!   and 3rd dimension is ungridded
!
! Args:
    integer,intent(in) :: nlev ! if nlev == 0, field is 2d (i,j)
    type(ESMF_Grid),intent(in) :: grid
    character(len=*),intent(in) :: name
    type(ESMF_Field),intent(out) :: field
!
! Local:
    integer :: rc
    type(ESMF_ArraySpec) :: arrayspec
!
! Create 3d field (i,j,k), with non-distributed vertical dimension:
    if (nlev > 0) then
      call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_geofield: ESMF_ArraySpecSet 3d field')
      field = ESMF_FieldCreate(grid, arrayspec,ungriddedLBound=(/1/), &
        ungriddedUBound=(/nlev/),staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_geofield: ESMF_FieldCreate 3d field')
!
! Create 2d field (i,j):
    else                ! create 2d field
      call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8,rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_geofield: ESMF_ArraySpecSet 2d field')
      field = ESMF_FieldCreate(grid, arrayspec,&
        staggerloc=ESMF_STAGGERLOC_CENTER,rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_geofield: ESMF_FieldCreate 2d field')
    endif
  end subroutine edyn_esmf_create_geofield
!-----------------------------------------------------------------------
  subroutine edyn_esmf_create_magfield(field,grid,name,nlev)
!
! Create ESMF field (2d or 3d) on mag grid. This will include the
!   mag periodic point, which will be zero after regridding.
! If nlev == 0, field is 2d (i,j), otherwise field is 3d,
!   and 3rd dimension is ungridded
!
! Args:
    integer,intent(in) :: nlev ! if nlev == 0, field is 2d (i,j)
    type(ESMF_Grid),intent(in) :: grid
    character(len=*),intent(in) :: name
    type(ESMF_Field),intent(out) :: field
!
! Local:
    integer :: rc
    type(ESMF_ArraySpec) :: arrayspec
    type(ESMF_Array) :: array3d,array2d
    type(ESMF_DistGrid) :: distgrid
!
! Get necessary information from the mag grid:
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER,&
      distgrid=distgrid,rc=rc)
    if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_magfield: ESMF_GridGet')
!
! Create 3d mag field (i,j,k), with non-distributed vertical dimension:
! (add periodic point in longitude with computationalEdgeUWidth)
!
    if (nlev > 0) then
      call ESMF_ArraySpecSet(arrayspec,3,ESMF_TYPEKIND_R8,rc=rc)
      if (rc /= ESMF_SUCCESS)call endrun('edyn_esmf_create_magfield: ESMF_ArraySpecSet 3d field')

      array3d = ESMF_ArrayCreate(arrayspec=arrayspec,      &
        distgrid=distgrid,computationalEdgeUWidth=(/1,0/), &
        undistLBound=(/1/),undistUBound=(/nlev/),          &
        indexflag=ESMF_INDEX_GLOBAL,rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_magfield: ESMF_ArrayCreate 3d field')

      field = ESMF_FieldCreate(grid, array3d,              &
        ungriddedLBound=(/1/), ungriddedUBound=(/nlev/),   &
        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_magfield: ESMF_FieldCreate 3d field')
!
! Create 2d mag field (i,j):
! (add periodic point in longitude with computationalEdgeUWidth)
!
    else                ! create 2d field
      call ESMF_ArraySpecSet(arrayspec,2,ESMF_TYPEKIND_R8,rc=rc)
      if (rc /= ESMF_SUCCESS)call endrun('edyn_esmf_create_magfield: ESMF_ArraySpecSet 2d field')

      array2d = ESMF_ArrayCreate(arrayspec=arrayspec,      &
        distgrid=distgrid,computationalEdgeUWidth=(/1,0/), &
        indexflag=ESMF_INDEX_GLOBAL,rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_magfield: ESMF_ArrayCreate 2d field')
      field = ESMF_FieldCreate(grid, array2d,              &
        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (rc /= ESMF_SUCCESS) call endrun('edyn_esmf_create_magfield: ESMF_FieldCreate 2d field')
    endif
  end subroutine edyn_esmf_create_magfield

!-----------------------------------------------------------------------
  subroutine edyn_esmf_set3d_geo(fields,fnames,f,nf,ilev0,ilev1,&
    ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 3d ESMF field on geographic source grid, prior to
!   geographic to magnetic grid transformation.
! Periodic points are excluded, geographic poles are at j==jspole and jnpole
! Note dimension order changes from input (k,i,j) to output (i,j,k).
!
! Args:
    integer,intent(in) :: nf
    type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on geo grid
    character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains (including periodic points)
! (note esmf source field excludes periodic points)
!
    integer,intent(in) :: ilev0,ilev1,ilon0,ilon1,ilat0,ilat1
    real(r8),intent(in) :: f(ilev0:ilev1,ilon0:ilon1,ilat0:ilat1,nf)
!
! Local:
    integer :: i,ii,j,k,rc,n,istat
    integer,parameter :: mxf=8 ! for call by dynamo_inputs
    integer :: lbnd(3),ubnd(3) ! lower,upper bounds of 3d field
!
! fptr is esmf pointer (i,j,k) to 3d field, set by this subroutine
    real(ESMF_KIND_R8),pointer :: fptr(:,:,:)
    real(r8),allocatable :: ftmp(:,:,:,:) ! esmf bounds, plus nf

    if (nf > mxf) then
      write(iulog,"('>>> esmf_set3d_geo: nf cannot be greater than mxf: nf=',i4,' mxf=',i4)") &
        nf,mxf
      call endrun('edyn_esmf_set3d_geo: nf > mxf')
    endif
    !
    ! This routine is called every timestep from dynamo_inputs for 8 fields,
    ! and is called once per run for a single field from geo2mag_3d  
    ! (called by define_phim3d, which is called from rdsource.F).
    !

!
! Get array bounds:
    call ESMF_FieldGet(fields(1),localDe=0,farrayPtr=fptr, &
         computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
    if (rc /= ESMF_SUCCESS) then
       write(iulog,"('>>> esmf_set3d_geo: error from ESMF_FieldGet: rc=',i4)") rc
       call endrun('edyn_esmf_set3d_geo: ESMF_FieldGet field 1')
    endif
!
! Do the allocation:
    allocate(ftmp(lbnd(1):ubnd(1),lbnd(2):ubnd(2),lbnd(3):ubnd(3),mxf),stat=istat)
    if (istat /= 0) then
       write(iulog,"('>>> esmf_set3d_geo: error allocating ftmp')")
       call endrun('edyn_esmf_set3d_geo: allocating ftmp')
    endif
!
! Fields loop:
    do n=1,nf
      ftmp(:,:,:,n) = 0._r8
!
! Set interior latitudes (ftmp(i,j,k,n) <- f(k,i,j,n))
! ftmp excludes periodic points.
!
      do j=lbnd(2),ubnd(2)      ! lat
        if (j /= jspole .and. j /= jnpole) then ! interior latitudes (not poles)
          do i=lbnd(1),ubnd(1)    ! lon
            ii = i
            do k=lbnd(3),ubnd(3)  ! lev
              ftmp(i,j,k,n) = f(k,ii,j,n)
            enddo ! lev
          enddo   ! lon
        endif     ! poles or interior
      enddo       ! lat
    enddo         ! n=1,nf
!
! Get and set pointer to the field:
    do n=1,nf
      call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr, &
        computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
      if (rc /= ESMF_SUCCESS) then
        write(iulog,"(a,i4)") '>>> esmf_set3d_geo: error from ESMF_FieldGet: rc=',rc
        call endrun('edyn_esmf_set3d_geo: ESMF_FieldGet field')
      endif
      fptr(:,:,:) = ftmp(:,:,:,n)
    enddo ! n=1,nf

    deallocate(ftmp)

  end subroutine edyn_esmf_set3d_geo
!-----------------------------------------------------------------------
  subroutine edyn_esmf_set2d_geo(field,grid,fname,f,ilon0,ilon1,ilat0,ilat1)
!   
! Set values of a 2d ESMF field on geographic source grid, prior to
!   geographic to magnetic grid transformation. (Essentially the same
!   as esmf_set3d_geo, except for 2d fields instead of 3d)
! Periodic points are excluded, geographic poles are at j==jspole and jnpole
!     
! Args: 
    type(ESMF_Field) ,intent(in) :: field
    type(ESMF_Grid)  ,intent(in) :: grid
    character(len=*) ,intent(in) :: fname  ! field name
    integer          ,intent(in) :: ilon0,ilon1,ilat0,ilat1
    real(r8)         ,intent(in) :: f(ilon0:ilon1,ilat0:ilat1)
!
! Local:
    integer :: i,ii,j,rc
    real(ESMF_KIND_R8),pointer :: fptr(:,:) ! i,j
    integer :: lbnd(2),ubnd(2)
!
! Get pointer to the field:
    call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr,&
      computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(a,i4)") '>>> esmf_set2d_geo: error from ESMF_FieldGet: rc=',rc
      call endrun('edyn_esmf_set2d_geo: ESMF_FieldGet')
    endif
!
    fptr(:,:) = 0._r8 ! init
!
! Set interior latitudes (excluding poles):
    do j=lbnd(2),ubnd(2)
      if (j /= jspole .and. j /= jnpole) then
        do i=lbnd(1),ubnd(1)
          ii = i
          fptr(i,j) = f(ii,j)
        enddo
      endif       ! interior latitudes only
    enddo

    if (debug) &
      write(iulog,"('esmf_set2d_geo field ',a,': lon bnds=',2i4, &
        ' lat bnds=',2i4,' 2d mnmx=',2e12.4)") &
        fname,lbnd(1),ubnd(1),lbnd(2),ubnd(2), &
        minval(fptr(:,:)),maxval(fptr(:,:))    

  end subroutine edyn_esmf_set2d_geo
!-----------------------------------------------------------------------
  subroutine edyn_esmf_set3d_mag(fields,fnames,f,nf,ilev0,ilev1,ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 3d ESMF field on magnetic grid, prior to magnetic to
! geographic grid transformation.
!
! Args:
  integer,intent(in) :: nf
  type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on mag grid
  character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains:
!
  integer,intent(in) :: ilev0,ilev1,ilon0,ilon1,ilat0,ilat1
  real(r8),intent(in) :: f(ilon0:ilon1,ilat0:ilat1,ilev0:ilev1,nf)
!
! Local:
  integer :: i,j,k,rc,n
  integer :: lbnd(3),ubnd(3) ! lower,upper bounds of 3d field
!
! fptr is esmf pointer (i,j,k) to 3d field, set by this subroutine
  real(ESMF_KIND_R8),pointer :: fptr(:,:,:)
!
! Fields loop:
  do n=1,nf
    call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr,&
      computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(a,i4)") '>>> esmf_set3d_mag: error from ESMF_FieldGet: rc=',rc
      call endrun('edyn_esmf_set3d_mag: ESMF_FieldGet')
    endif
!
    fptr(:,:,:) = 0._r8
!
! Set ESMF pointer:
!
    do j=lbnd(2),ubnd(2)      ! lat
      do i=lbnd(1),ubnd(1)    ! lon
        do k=lbnd(3),ubnd(3)  ! lev
          fptr(i,j,k) = f(i,j,k,n)
        enddo ! mlev
      enddo   ! mlon
    enddo     ! mlat
  enddo       ! n=1,nf
  end subroutine edyn_esmf_set3d_mag
!-----------------------------------------------------------------------
!
  subroutine edyn_esmf_set2d_mag(fields,fnames,f,nf,ilon0,ilon1,ilat0,ilat1)
!
! Set values of a 2d ESMF field on magnetic grid, prior to magnetic to
! geographic grid transformation.
!
! Args:
  integer,intent(in) :: nf
  type(ESMF_Field) ,intent(in) :: fields(nf) ! esmf fields on mag grid
  character(len=*) ,intent(in) :: fnames(nf) ! field names
!
! f is input data on model subdomains:
!
  integer,intent(in) :: ilon0,ilon1,ilat0,ilat1
  real(r8),intent(in) :: f(ilon0:ilon1,ilat0:ilat1,nf)
!
! Local:
  integer :: i,j,rc,n
  integer :: lbnd(2),ubnd(2) ! lower,upper bounds of 2d field
!
! fptr is esmf pointer (i,j,k) to 2d field, set by this subroutine
  real(ESMF_KIND_R8),pointer :: fptr(:,:)
!
! Fields loop:
  do n=1,nf
    call ESMF_FieldGet(fields(n),localDe=0,farrayPtr=fptr,&
      computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
    if (rc /= ESMF_SUCCESS) then
      write(iulog,"(a,i4)") '>>> esmf_set2d_mag: error from ESMF_FieldGet: rc=',rc
      call endrun('edyn_esmf_set2d_mag: ESMF_FieldGet')
    endif
!
    fptr(:,:) = 0._r8
!
! Set ESMF pointer:
!
    do j=lbnd(2),ubnd(2)      ! lat
      do i=lbnd(1),ubnd(1)    ! lon
        fptr(i,j) = f(i,j,n)
      enddo   ! mlon
    enddo     ! mlat
  enddo       ! n=1,nf
!
  end subroutine edyn_esmf_set2d_mag
!-----------------------------------------------------------------------
  subroutine edyn_esmf_get_3dfield(field, fptr, name)
!
! Get pointer to 3d esmf field (i,j,k):
!
! Args:
    type(ESMF_field),intent(in) :: field
    real(r8),pointer,dimension(:,:,:),intent(out) :: fptr
    character(len=*),intent(in) :: name
!
! Local:
    integer :: rc,lbnd(3),ubnd(3)
    character(len=80) :: errmsg

    call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr, &
      computationalLBound=lbnd,computationalUBound=ubnd,rc=rc)
    if (rc /= ESMF_SUCCESS) then
      write(errmsg,"('esmf_get_field 3d field ',a)") trim(name)
      call endrun('edyn_esmf_get_3dfield: ESMF_FieldGet')
    endif
  end subroutine edyn_esmf_get_3dfield
!-----------------------------------------------------------------------
  subroutine edyn_esmf_get_2dfield(field, fptr, name)
!
! Get pointer to 2d esmf field (i,j):
!
! Args:
    type(ESMF_field),intent(in) :: field
    real(r8),pointer,dimension(:,:),intent(out) :: fptr
    character(len=*),intent(in) :: name
!
! Local:
    integer :: rc
    character(len=80) :: errmsg

    call ESMF_FieldGet(field,localDe=0,farrayPtr=fptr,rc=rc)
    if (rc /= ESMF_SUCCESS) then
      write(errmsg,"('edyn_esmf_get_2dfield ',a)") trim(name)
      call endrun('edyn_esmf_get_2dfield: ESMF_FieldGet')
    endif

  end subroutine edyn_esmf_get_2dfield
!-----------------------------------------------------------------------
  subroutine edyn_esmf_regrid(srcfield,dstfield,direction,ndim)
!
! Args:
    integer :: ndim
    type(ESMF_Field),intent(inout) :: srcfield,dstfield
    character(len=*),intent(in) :: direction
!
! Local:
    integer :: rc
    type(ESMF_RouteHandle) :: routehandle
!
! Direction is either geo2mag or mag2geo.
! Use corresponding route handle (module data)
!
    select case(trim(direction))
      case ('geo2mag')
        routehandle = routehandle_geo2mag
        if (ndim==2) then
!
! Do sparse matrix multiply for 2d geo2mag.
!
          routehandle = routehandle_geo2mag_2d
          call ESMF_FieldSMM(srcfield,dstfield,routehandle,termorderflag=ESMF_TERMORDER_SRCSEQ,rc=rc)

          if (rc /= ESMF_SUCCESS) then
              write(iulog,"(/,4a,i4)") '>>> edyn_esmf_regrid: error return from ',&
                'ESMF_FieldSMM for 2d ',trim(direction),': rc=',rc
              call endrun('edyn_esmf_regrid: ESMF_FieldSMM 2d')
          endif
        else ! 3d geo2mag
!
! Do sparse matrix multiply for 3d geo2mag. 
!
          routehandle = routehandle_geo2mag
          call ESMF_FieldSMM(srcfield,dstfield,routehandle,termorderflag=ESMF_TERMORDER_SRCSEQ,rc=rc)
          if (rc /= ESMF_SUCCESS) then
            write(iulog,"(/,4a,i4)") '>>> edyn_esmf_regrid: error return from ',&
              'ESMF_FieldSMM for 3d ',trim(direction),': rc=',rc
            call endrun('edyn_esmf_regrid: ESMF_FieldSMM 3d')
          endif
        endif
!
! Do sparse matrix multiply for 3d mag2geo.
! btf 6/18/14: mag2geo is not working due to error return rc=51 from the 
!   below call. Calls to mag2geo_3d at end of sub pefield (edynamo.F90)
!   are commented out (mag2geo_3d calls this routine with direction='mag2geo').
!
      case ('mag2geo')
         if (ndim==2) then
            routehandle = routehandle_mag2geo_2d
         else
            routehandle = routehandle_mag2geo
         endif
         call ESMF_FieldSMM(srcfield,dstfield,routehandle,termorderflag=ESMF_TERMORDER_SRCSEQ,checkflag=.true.,rc=rc)
         if (rc /= ESMF_SUCCESS) then
            write(iulog,"(/,4a,i4)") '>>> edyn_esmf_regrid: error return from ',&
                 'ESMF_FieldSMM for 3d ',trim(direction),': rc=',rc
            call endrun('edyn_esmf_regrid: ESMF_FieldSMM magtogeo')
         endif
      case default
        write(iulog,"('>>> edyn_esmf_regrid: bad direction=',a)") trim(direction)
        call endrun
    end select
  end subroutine edyn_esmf_regrid
!-----------------------------------------------------------------------

  subroutine edyn_esmf_update_flag( flag )
    logical, intent(in) :: flag
    edyn_esmf_update_step=flag
  end subroutine edyn_esmf_update_flag
  
#endif
end module edyn_esmf
