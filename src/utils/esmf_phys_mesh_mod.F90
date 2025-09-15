!-------------------------------------------------------------------------------
! Encapsulates the CAM physics grid mesh
!-------------------------------------------------------------------------------
module esmf_phys_mesh_mod
  use shr_kind_mod,   only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use ESMF,           only: ESMF_DistGrid, ESMF_DistGridCreate, ESMF_MeshCreate
  use ESMF,           only: ESMF_FILEFORMAT_ESMFMESH,ESMF_MeshGet,ESMF_Mesh, ESMF_SUCCESS
  use ESMF,           only: ESMF_MeshDestroy, ESMF_DistGridDestroy
  use esmf_check_error_mod, only: check_esmf_error

  implicit none

  private

  public :: esmf_phys_mesh_init
  public :: esmf_phys_mesh_destroy
  public :: physics_grid_mesh

  ! phys_mesh: Local copy of physics grid
  type(ESMF_Mesh), protected :: physics_grid_mesh

  ! dist_grid_2d: DistGrid for 2D fields
  type(ESMF_DistGrid) :: dist_grid_2d

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine esmf_phys_mesh_init()
    use phys_control, only: phys_getopts
    use phys_grid,    only: get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
    use ppgrid,       only: pcols, begchunk, endchunk
    use shr_const_mod,only: shr_const_pi

    ! Local variables
    integer                               :: ncols
    integer                               :: chnk, col, dindex
    integer,                allocatable   :: decomp(:)
    character(len=cl)                     :: grid_file
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(r8), pointer       :: ownedElemCoords(:)
    real(r8), pointer       :: lat(:), latMesh(:)
    real(r8), pointer       :: lon(:), lonMesh(:)
    real(r8)                :: lats(pcols)                       ! array of chunk latitudes
    real(r8)                :: lons(pcols)                       ! array of chunk longitude
    character(len=cs)       :: tempc1,tempc2
    character(len=300)      :: errstr

    integer :: i, c, n, total_cols, rc

    real(r8), parameter :: abstol =  1.e-3_r8
    real(r8), parameter :: radtodeg = 180.0_r8/shr_const_pi
    character(len=*), parameter :: subname = 'esmf_phys_mesh_init: '

    ! Find the physics grid file
    call phys_getopts(physics_grid_out=grid_file)

    ! Compute the local decomp
    total_cols = 0
    do chnk = begchunk, endchunk
       total_cols = total_cols + get_ncols_p(chnk)
    end do
    allocate(decomp(total_cols), stat=rc)
    if (rc/=0) then
       call endrun(subname//'not able to allocate decomp')
    end if

    dindex = 0
    do chnk = begchunk, endchunk
       ncols = get_ncols_p(chnk)
       do col = 1, ncols
          dindex = dindex + 1
          decomp(dindex) = get_gcol_p(chnk, col)
       end do
    end do

    ! Create a DistGrid based on the physics decomp
    dist_grid_2d = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_DistGridCreate')

    ! Create an ESMF_mesh for the physics decomposition
    physics_grid_mesh = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH,  &
                                        elementDistgrid=dist_grid_2d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_MeshCreate')

    ! Check that the mesh coordinates are consistent with the model physics column coordinates

    ! obtain mesh lats and lons
    call ESMF_MeshGet(physics_grid_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_MeshGet')

    if (numOwnedElements /= total_cols) then
       write(tempc1,'(i10)') numOwnedElements
       write(tempc2,'(i10)') total_cols
       call endrun(subname//"ERROR numOwnedElements "// &
            trim(tempc1) //" not equal to local size "// trim(tempc2))
    end if

    allocate(ownedElemCoords(spatialDim*numOwnedElements), stat=rc)
    if (rc/=0) then
       call endrun(subname//'not able to allocate ownedElemCoords')
    end if

    allocate(lonMesh(total_cols), stat=rc)
    if (rc/=0) then
       call endrun(subname//'not able to allocate lonMesh')
    end if

    allocate(latMesh(total_cols), stat=rc)
    if (rc/=0) then
       call endrun(subname//'not able to allocate latMesh')
    end if

    call ESMF_MeshGet(physics_grid_mesh, ownedElemCoords=ownedElemCoords)
    call check_esmf_error(rc, subname//'ESMF_MeshGet')

    do n = 1,total_cols
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    ! obtain internally generated cam lats and lons
    allocate(lon(total_cols), stat=rc);
    if (rc/=0) then
       call endrun(subname//'not able to allocate lon')
    end if

    lon(:) = 0._r8

    allocate(lat(total_cols), stat=rc);
    if (rc/=0) then
       call endrun(subname//'not able to allocate lat')
    end if

    lat(:) = 0._r8

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       ! latitudes and longitudes returned in radians
       call get_rlat_all_p(c, ncols, lats)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          lat(n) = lats(i)*radtodeg
          lon(n) = lons(i)*radtodeg
       end do
    end do

    errstr = ''
    ! error check differences between internally generated lons and those read in
    do n = 1,total_cols
       if (abs(lonMesh(n) - lon(n)) > abstol) then
          if ( (abs(lonMesh(n)-lon(n)) > 360._r8+abstol) .or. (abs(lonMesh(n)-lon(n)) < 360._r8-abstol) ) then
             write(errstr,100) n,lon(n),lonMesh(n), abs(lonMesh(n)-lon(n))
             write(iulog,*) trim(errstr)
          endif
       end if
       if (abs(latMesh(n) - lat(n)) > abstol) then
          ! poles in the 4x5 SCRIP file seem to be off by 1 degree
          if (.not.( (abs(lat(n))>88.0_r8) .and. (abs(latMesh(n))>88.0_r8) )) then
             write(errstr,101) n,lat(n),latMesh(n), abs(latMesh(n)-lat(n))
             write(iulog,*) trim(errstr)
          endif
       end if
    end do

    if ( len_trim(errstr) > 0 ) then
       call endrun(subname//'physics mesh coords do not match model coords')
    end if

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)
    deallocate(decomp)

100 format('esmf_phys_mesh_init: coord mismatch... n, lon(n), lonmesh(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
101 format('esmf_phys_mesh_init: coord mismatch... n, lat(n), latmesh(n), diff_lat = ',i6,2(f21.13,3x),d21.5)

  end subroutine esmf_phys_mesh_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine esmf_phys_mesh_destroy()

    integer :: rc
    character(len=*), parameter :: subname = 'esmf_phys_mesh_destroy: '

    call ESMF_MeshDestroy(physics_grid_mesh, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_MeshDestroy phys_mesh')

    call ESMF_DistGridDestroy(dist_grid_2d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_DistGridDestroy dist_grid_2d')

  end subroutine esmf_phys_mesh_destroy

end module esmf_phys_mesh_mod
