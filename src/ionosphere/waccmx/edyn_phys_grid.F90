!-------------------------------------------------------------------------------
! Initializes the CAM physics grid mesh
!-------------------------------------------------------------------------------
module edyn_phys_grid
  use shr_kind_mod,   only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none

  private

  public :: edyn_phys_grid_init

contains

  subroutine edyn_phys_grid_init()
    use ESMF,         only: ESMF_DistGrid, ESMF_DistGridCreate, ESMF_MeshCreate
    use ESMF,         only: ESMF_FILEFORMAT_ESMFMESH,ESMF_MeshGet,ESMF_Mesh
    use phys_control, only: phys_getopts
    use phys_grid,    only: get_ncols_p, get_gcol_p, get_rlon_all_p, get_rlat_all_p
    use ppgrid,       only: begchunk, endchunk
    use edyn_esmf,    only: edyn_esmf_chkerr, edyn_esmf_update_phys_mesh
    use shr_const_mod,only: shr_const_pi
    use ppgrid,       only: pcols
    use error_messages,only: alloc_err

    ! Local variables
    integer                               :: ncols
    integer                               :: chnk, col, dindex
    integer,                allocatable   :: decomp(:)
    character(len=cl)                     :: grid_file
    character(len=*),           parameter :: subname = 'edyn_gcomp_init'
    real(r8)        ,           parameter :: radtodeg = 180.0_r8/shr_const_pi
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(r8), pointer       :: ownedElemCoords(:)
    real(r8), pointer       :: lat(:), latMesh(:)
    real(r8), pointer       :: lon(:), lonMesh(:)
    real(r8)                :: lats(pcols)                       ! array of chunk latitudes
    real(r8)                :: lons(pcols)                       ! array of chunk longitude
    integer :: i, c, n
    character(len=cs)       :: tempc1,tempc2
    character(len=300)      :: errstr

    ! dist_grid_2d: DistGrid for 2D fields
    type(ESMF_DistGrid) :: dist_grid_2d

    ! phys_mesh: Local copy of physics grid
    type(ESMF_Mesh) :: phys_mesh

    real(r8), parameter :: abstol =  1.e-6_r8
    integer :: total_cols, rc

    ! Find the physics grid file
    call phys_getopts(physics_grid_out=grid_file)
    ! Compute the local decomp
    total_cols = 0
    do chnk = begchunk, endchunk
       total_cols = total_cols + get_ncols_p(chnk)
    end do
    allocate(decomp(total_cols), stat=rc)
    call alloc_err(rc,subname,'decomp',total_cols)

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
    call edyn_esmf_chkerr(subname, 'ESMF_DistGridCreate phys decomp', rc)

    ! Create an ESMF_mesh for the physics decomposition
    phys_mesh = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH,  &
         elementDistgrid=dist_grid_2d, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_MeshCreateFromFile', rc)

    call edyn_esmf_update_phys_mesh(phys_mesh)

    ! Check that the mesh coordinates are consistent with the model physics column coordinates

    ! obtain mesh lats and lons
    call ESMF_MeshGet(phys_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    call edyn_esmf_chkerr(subname, 'ESMF_MeshGet', rc)

    if (numOwnedElements /= total_cols) then
       write(tempc1,'(i10)') numOwnedElements
       write(tempc2,'(i10)') total_cols
       call endrun(trim(subname)//": ERROR numOwnedElements "// &
            trim(tempc1) //" not equal to local size "// trim(tempc2))
    end if

    allocate(ownedElemCoords(spatialDim*numOwnedElements), stat=rc)
    call alloc_err(rc,subname,'ownedElemCoords',spatialDim*numOwnedElements)

    allocate(lonMesh(total_cols), stat=rc)
    call alloc_err(rc,subname,'lonMesh',total_cols)

    allocate(latMesh(total_cols), stat=rc)
    call alloc_err(rc,subname,'latMesh',total_cols)

    call ESMF_MeshGet(phys_mesh, ownedElemCoords=ownedElemCoords)

    do n = 1,total_cols
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    ! obtain internally generated cam lats and lons
    allocate(lon(total_cols), stat=rc);
    call alloc_err(rc,subname,'lon',total_cols)

    lon(:) = 0._r8

    allocate(lat(total_cols), stat=rc);
    call alloc_err(rc,subname,'lat',total_cols)

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
       call endrun(subname//': physics mesh coords do not match model coords')
    end if

    ! deallocate memory
    deallocate(ownedElemCoords)
    deallocate(lon, lonMesh)
    deallocate(lat, latMesh)
    deallocate(decomp)

100 format('edyn_gcomp_init: coord mismatch... n, lon(n), lonmesh(n), diff_lon = ',i6,2(f21.13,3x),d21.5)
101 format('edyn_gcomp_init: coord mismatch... n, lat(n), latmesh(n), diff_lat = ',i6,2(f21.13,3x),d21.5)

  end subroutine edyn_phys_grid_init


end module edyn_phys_grid
