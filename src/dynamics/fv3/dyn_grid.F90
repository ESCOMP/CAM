module dyn_grid
!-------------------------------------------------------------------------------
! Define FV3 computational grids on the dynamics decomposition.
!
! The grid used by the FV3 dynamics is called the FSSL grid and is a
! gnomonic cubed sphere consisting of 6 tiled faces. Each tile consists
! of an array of cells whose coordinates are great circles. The grid
! nomenclature (C96, C384, etc.) describes the number of cells along
! the top and side of a tile face (square). All prognostic variables
! are 3-D cell-mean values (cell center), except for the horizontal winds,
! which are 2-D face-mean values located on the cell walls (D-Grid winds).
! Each tile can be decomposed into a number of subdomains (consisting of
! one or more cells) which correspond to "blocks" in the physics/dynamics
! coupler terminology. The namelist variable "layout" consists of 2 integers
! and determines the size/shape of the blocks by dividing the tile into a
! number of horizonal and vertical sections. The total number of blocks in
! the global domain is therefore layout(1)*layout(2)*ntiles. The decomposition
! and communication infrastructure is provided by the GFDL FMS library.
!
! Module responsibilities:
!
! . Provide the physics/dynamics coupler (in module phys_grid) with data for the
!   physics grid on the dynamics decomposition.
!
! . Create CAM grid objects that are used by the I/O functionality to read
!   data from an unstructured grid format to the dynamics data structures, and
!   to write from the dynamics data structures to unstructured grid format.  The
!   global column ordering for the unstructured grid is determined by the FV3 dycore.
!
!-------------------------------------------------------------------------------

    use cam_abortutils,   only: endrun
    use cam_grid_support, only: iMap
    use cam_logfile,      only: iulog
    use dimensions_mod,   only: npx, npy, ntiles
    use fms_mod,          only: fms_init, write_version_number
    use fv_arrays_mod,    only: fv_atmos_type
    use fv_control_mod,   only: ngrids,fv_init
    use fv_mp_mod,        only: mp_bcst
    use mpp_mod,          only: mpp_pe, mpp_root_pe
    use physconst,        only: rearth,pi
    use shr_kind_mod,     only: r8 => shr_kind_r8
    use spmd_utils,       only: mpicom, masterproc

    implicit none
    private
    save

    ! The FV3 dynamics grids
    integer, parameter :: dyn_decomp   = 101
    integer, parameter :: dyn_decomp_ew = 102
    integer, parameter :: dyn_decomp_ns = 103
    integer, parameter :: dyn_decomp_hist   = 104
    integer, parameter :: dyn_decomp_hist_ew = 105
    integer, parameter :: dyn_decomp_hist_ns = 106

    integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore

    integer :: mytile = 1
    integer :: p_split  = 1
    integer, allocatable :: pelist(:)

    real(r8), parameter :: rad2deg = 180._r8/pi

    logical, allocatable :: grids_on_this_pe(:)
    type(fv_atmos_type), allocatable, target :: Atm(:)


public ::             &
   dyn_decomp,        &
   p_split,           &
   grids_on_this_pe,  &
   ptimelevels

!-----------------------------------------------------------------------
! Calculate Global Index

integer, allocatable, target, dimension(:,:) :: mygindex
integer, allocatable, target, dimension(:,:) :: mylindex
integer, allocatable, target, dimension(:,:) :: myblkidx
real(r8), allocatable, target, dimension(:,:,:) ::   locidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   blkidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   gindex_g

real(r8), allocatable :: block_extents_g(:,:)

integer :: uniqpts_glob = 0     ! number of dynamics columns
integer :: uniqpts_glob_ew = 0     ! number of dynamics columns for D grid ew
integer :: uniqpts_glob_ns = 0     ! number of dynamics columns for D grid ns

real(r8), pointer, dimension(:,:,:) :: grid_ew, grid_ns

public :: mygindex
public :: mylindex
!-----------------------------------------------------------------------
public :: &
   dyn_grid_init,            &
   get_block_bounds_d,       & ! get first and last indices in global block ordering
   get_block_gcol_d,         & ! get column indices for given block
   get_block_gcol_cnt_d,     & ! get number of columns in given block
   get_block_lvl_cnt_d,      & ! get number of vertical levels in column
   get_block_levels_d,       & ! get vertical levels in column
   get_block_owner_d,        & ! get process "owning" given block
   get_gcol_block_d,         & ! get global block indices and local columns
                               ! index for given global column index
   get_gcol_block_cnt_d,     & ! get number of blocks containing data
                               ! from a given global column index
   get_horiz_grid_dim_d,     &
   get_horiz_grid_d,         & ! get horizontal grid coordinates
   get_dyn_grid_parm,        &
   get_dyn_grid_parm_real1d, &
   dyn_grid_get_elem_coords, & ! get coordinates of a specified block element
   dyn_grid_get_colndx,      & ! get element block/column and MPI process indices
                               ! corresponding to a specified global column index
   physgrid_copy_attributes_d

public Atm, mytile

!=======================================================================
contains
!=======================================================================

subroutine dyn_grid_init()

   ! Initialize FV grid, decomposition

   use block_control_mod,  only: block_control_type, define_blocks_packed
   use cam_initfiles,      only: initial_file_get_id
   use constants_mod,      only: constants_init
   use fv_mp_mod,          only: switch_current_Atm,mp_gather, mp_bcst
   use hycoef,             only: hycoef_init, hyai, hybi, hypi, hypm, nprlev
   use mpp_mod,            only: mpp_init, mpp_npes, mpp_get_current_pelist,mpp_gather
   use pmgrid,             only: plev
   use ref_pres,           only: ref_pres_init
   use time_manager,       only: get_step_size
   use pio,                only: file_desc_t

   ! Local variables

   type(file_desc_t),      pointer :: fh_ini

   character(len=*), parameter :: sub='dyn_grid_init'
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

   real(r8) :: dt_atmos_real = 0._r8

   integer :: i, j, k, tile
   integer :: is,ie,js,je,n,nx,ny
   character(len=128) :: errmsg

   !-----------------------------------------------------------------------
   !  from couple_main initialize atm structure - initializes fv3 grid
   !-----------------------------------------------------------------------

   call fms_init(mpicom)
   call mpp_init()
   call constants_init

!-----------------------------------------------------------------------
! initialize atmospheric model -----

   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)

!---- compute physics/atmos time step in seconds ----

   dt_atmos_real = get_step_size()

!----- initialize FV dynamical core -----

   call fv_init( Atm, dt_atmos_real, grids_on_this_pe, p_split)  ! allocates Atm components

   do n=1,ngrids
      if (grids_on_this_pe(n)) mytile = n
   enddo

!----- write version and namelist to log file -----
   call write_version_number ( version, tagname )

   call switch_current_Atm(Atm(mytile))

!!  set up dimensions_mod convenience variables.

   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je
   npx   = Atm(mytile)%flagstruct%npx
   npy   = Atm(mytile)%flagstruct%npy
   ntiles = Atm(mytile)%gridstruct%ntiles_g
   tile = Atm(mytile)%tile

   if (Atm(mytile)%flagstruct%npz /= plev) then
      write(errmsg,*) 'FV3 dycore levels (npz),',Atm(mytile)%flagstruct%npz,' do not match model levels (plev)',plev
      call endrun(sub//':'//errmsg)
   end if

   ! Get file handle for initial file
   fh_ini => initial_file_get_id()

   ! Initialize hybrid coordinate arrays
   call hycoef_init(fh_ini)

   ! Initialize reference pressures
   call ref_pres_init(hypi, hypm, nprlev)

   ! Hybrid coordinate info for FV grid object
   Atm(mytile)%ks = plev
   do k = 1, plev+1
      Atm(mytile)%ak(k) = hyai(k) * 1.e5_r8
      Atm(mytile)%bk(k) = hybi(k)
      if ( Atm(mytile)%bk(k) == 0._r8)  Atm(mytile)%ks = k-1
   end do
   Atm(mytile)%ptop = Atm(mytile)%ak(1)

   ! Define the CAM grids
   call define_cam_grids(Atm)

   ! Define block index arrays that are part of dyn_in and
   ! global array for mapping columns to block decompositions

   allocate(mygindex(is:ie,js:je))
   allocate(mylindex(is:ie,js:je))

   nx=npx-1
   ny=npy-1

   n = 1
   do j = js, je
      do i = is, ie
         mygindex(i,j)=((j-1)*(npx-1)+i)+((npx-1)*(npy-1)*(tile-1))
         mylindex(i,j)=n
         n = n + 1
      end do
   end do

   ! create globalID index on block decomp
   allocate(gindex_g(nx,ny,ntiles))
   if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.(gindex_g)'
   gindex_g(is:ie,js:je,tile)=mygindex(is:ie,js:je)
   call mp_gather(gindex_g, is, ie, js, je, nx, ny, ntiles)
   call mp_bcst(gindex_g, nx, ny, ntiles)

   ! create global blockID index on block decomp
   if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.(blkidx_g)'
   allocate(blkidx_g(nx,ny,ntiles))
   blkidx_g(is:ie,js:je,tile)= mpp_pe() + 1
   call mp_gather(blkidx_g, is, ie, js, je, nx ,ny, ntiles)
   call mp_bcst(blkidx_g, nx, ny, ntiles)

   ! create global block index on block decomp
   if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.(locidx_g)'
   allocate(locidx_g(nx,ny,ntiles))
   locidx_g(is:ie,js:je,tile)= mylindex(is:ie,js:je)
   call mp_gather(locidx_g, is, ie, js, je, nx ,ny, ntiles)
   call mp_bcst(locidx_g, nx, ny, ntiles)

end subroutine dyn_grid_init

!=======================================================================

subroutine get_block_bounds_d(block_first, block_last)

    ! Return first and last indices used in global block ordering

    use spmd_utils,    only : npes

    ! arguments
    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks
    !----------------------------------------------------------------------------

    block_first = 1
    block_last = npes

end subroutine get_block_bounds_d

!=======================================================================

subroutine get_block_gcol_d(blockid, size, cdex)

    ! Return number of dynamics columns in indicated block

    use fv_mp_mod,          only: mp_bcst
    use mpp_mod,            only: mpp_npes, mpp_gather

    ! arguments
    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size
    integer, intent(out):: cdex(size)   ! global column indices

   ! Local variables
    integer, parameter :: be_arrlen = 5

    real(r8),allocatable :: rtmp(:)
    real(r8) :: block_extents(be_arrlen)
    integer, allocatable :: be_size(:)
    integer :: i, j, n,is,ie,js,je,tile,npes
    !----------------------------------------------------------------------------
    !--- get block extents for each task/pe

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    if (.not. allocated(block_extents_g)) then
       npes=mpp_npes()
       allocate(block_extents_g(be_arrlen,npes))
       allocate(rtmp(be_arrlen*npes))
       allocate(be_size(npes))
       be_size(:)=be_arrlen
       block_extents(1)=is
       block_extents(2)=ie
       block_extents(3)=js
       block_extents(4)=je
       block_extents(5)=Atm(mytile)%tile

       call mpp_gather(block_extents,be_arrlen,rtmp,be_size)
       call mp_bcst(rtmp,be_arrlen*npes)
       block_extents_g=reshape(rtmp,(/be_arrlen,npes/))

       deallocate(rtmp)
       deallocate(be_size)
    end if

    is=block_extents_g(1,blockid)
    ie=block_extents_g(2,blockid)
    js=block_extents_g(3,blockid)
    je=block_extents_g(4,blockid)
    tile=block_extents_g(5,blockid)

    if (size .ne. (ie - is + 1) * (je - js + 1)) then
       call endrun ('get_block_gcol_d: block sizes are not consistent.')
    end if
    ! the following algorithm for cdex calculates global ids for a block
    ! given the tile,and i,j column locations on tile.
    n=1
    do j = js, je
        do i = is, ie
            cdex(n)= ((j-1)*(npx-1)+i)+((npx-1)*(npy-1)*(tile-1))
            n=n+1
         end do
    end do

end subroutine get_block_gcol_d

!=======================================================================

integer function get_block_gcol_cnt_d(blockid)

    ! Return number of dynamics columns in indicated block

    ! arguments
    integer, intent(in) :: blockid
    !----------------------------------------------------------------------------

    get_block_gcol_cnt_d=count(blkidx_g == blockid)

end function get_block_gcol_cnt_d

!=======================================================================

integer function get_block_lvl_cnt_d(blockid, bcid)

   ! Return number of levels in indicated column. If column
   ! includes surface fields, then it is defined to also
   ! include level 0.

    use pmgrid, only: plevp

    ! arguments
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid     ! column index within block
    !----------------------------------------------------------------------------

    get_block_lvl_cnt_d = plevp

end function get_block_lvl_cnt_d

!=======================================================================

subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

    use pmgrid, only: plev

    ! Return level indices in indicated column. If column
    ! includes surface fields, then it is defined to also
    ! include level 0.

    ! arguments
    integer, intent(in)  :: blockid        ! global block id
    integer, intent(in)  :: bcid           ! column index within block
    integer, intent(in)  :: lvlsiz         ! dimension of levels array
    integer, intent(out) :: levels(lvlsiz) ! levels indices for block

    ! local variables
    integer :: k
    character(len=128) :: errmsg
    !---------------------------------------------------------------------------

    if (lvlsiz < plev + 1) then
       write(errmsg,*) 'levels array not large enough (', lvlsiz,' < ',plev + 1,')'
       call endrun('GET_BLOCK_LEVELS_D: '//trim(errmsg))
    else
        do k = 0, plev
            levels(k+1) = k
        enddo
        do k = plev + 2, lvlsiz
            levels(k) = -1
        enddo
    end if

end subroutine get_block_levels_d

!=======================================================================

integer function get_block_owner_d(blockid)

    ! Return id of processor that "owns" the indicated block

    ! arguments
    integer, intent(in) :: blockid  ! global block id

    get_block_owner_d = blockid - 1

end function get_block_owner_d

!=======================================================================

subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)

    ! Return global block index and local column index for given global column index.
    !
    ! The FV3 dycore assigns each global column to a singe element.  So cnt is assumed
    ! to be 1.

    use dimensions_mod,     only: npx, npy
    use fv_mp_mod,          only: mp_gather, mp_bcst

    ! arguments
    integer, intent(in)  :: gcol         ! global column index
    integer, intent(in)  :: cnt          ! size of blockid and bcid arrays
    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)

    ! local variables
    integer :: tot
    integer :: ijk(3)
    !----------------------------------------------------------------------------

    if (cnt /= 1) then
       call endrun ('get_gcol_block_d: cnt is not equal to 1:.')
    end if
    tot=(npx-1)*(npy-1)*6
    if (gcol < 1.or.gcol > tot) then
       call endrun ('get_gcol_block_d: global column number is out of bounds')
    else

       ijk=maxloc(blkidx_g,mask=gindex_g == gcol)
       blockid(1) = blkidx_g(ijk(1),ijk(2),ijk(3))

       ijk=maxloc(locidx_g,mask=gindex_g == gcol)
       bcid(1) = locidx_g(ijk(1),ijk(2),ijk(3))
    end if

    if (present(localblockid)) then
        localblockid(cnt) = 1
    end if

end subroutine get_gcol_block_d

!=======================================================================

integer function get_gcol_block_cnt_d(gcol)

    ! Return number of blocks containg data for the vertical column with the
    ! given global column index.

    ! For FV3 dycore each column is contained in a single block, so this routine
    ! always returns 1.

    ! arguments
    integer, intent(in) :: gcol     ! global column index
    !----------------------------------------------------------------------------

    get_gcol_block_cnt_d = 1

end function get_gcol_block_cnt_d

!=======================================================================

subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

    ! Return global arrays of latitude and longitude (in radians), column
    ! surface area (in radians squared) and surface integration weights for
    ! global column indices that will be passed to/from physics

    ! arguments
    integer, intent(in)             :: nxy           ! array sizes
    real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes
    real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes
    real(r8), intent(out), optional :: area_d_out(:) ! column surface area
    real(r8), intent(out), optional :: wght_d_out(:) ! column integration
    real(r8), intent(out), optional :: lat_d_out(:)  ! column degree latitudes
    real(r8), intent(out), optional :: lon_d_out(:)  ! column degree longitudes

    ! local variables
    character(len=*), parameter :: sub = 'get_horiz_grid_d'
    real(r8), allocatable       :: tmparr(:,:)
    real(r8), pointer           :: area(:,:)
    real(r8), pointer           :: agrid(:,:,:)
    integer                     :: is,ie,js,je
    !----------------------------------------------------------------------------

    area  => Atm(mytile)%gridstruct%area_64
    agrid => Atm(mytile)%gridstruct%agrid_64
    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    if (present(clon_d_out)) then
       if (size(clon_d_out) /= nxy) call endrun(sub//': bad clon_d_out array size')
       call create_global(is,ie,js,je,agrid(is:ie,js:je,1), clon_d_out)
    end if
    if (present(clat_d_out)) then
       if (size(clat_d_out) /= nxy) call endrun(sub//': bad clat_d_out array size')
       call create_global(is,ie,js,je,agrid(is:ie,js:je,2), clat_d_out)
    end if
    if (present(area_d_out).or.present(wght_d_out)) then
       allocate(tmparr(is:ie,js:je))
       tmparr(is:ie,js:je) = area  (is:ie,js:je) / (rearth * rearth)
       if (present(area_d_out)) then
          if (size(area_d_out) /= nxy) call endrun(sub//': bad area_d_out array size')
          call create_global(is,ie,js,je,tmparr, area_d_out)
       end if
       if (present(wght_d_out)) then
          if (size(wght_d_out) /= nxy) call endrun(sub//': bad wght_d_out array size')
          call create_global(is,ie,js,je,tmparr, wght_d_out)
       end if
       deallocate(tmparr)
    end if
    if (present(lon_d_out)) then
       if (size(lon_d_out) /= nxy) call endrun(sub//': bad clon_d_out array size')
       call create_global(is,ie,js,je,agrid(is:ie,js:je,1), lon_d_out)
       lon_d_out=lon_d_out*rad2deg
    end if
    if (present(lat_d_out)) then
       if (size(lat_d_out) /= nxy) call endrun(sub//': bad clat_d_out array size')
       call create_global(is,ie,js,je,agrid(is:ie,js:je,2), lat_d_out)
       lat_d_out=lat_d_out*rad2deg
    end if

  end subroutine get_horiz_grid_d

!=======================================================================

subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)

    ! Returns declared horizontal dimensions of computational grid.
    ! For non-lon/lat grids, declare grid to be one-dimensional,

    use dimensions_mod, only: npx,npy,ntiles

    ! arguments
    integer, intent(out) :: hdim1_d                 ! first horizontal dimension
    integer, intent(out), optional :: hdim2_d       ! second horizontal dimension
    !-----------------------------------------------------------------------

    hdim1_d = (npx-1)*(npy-1)*ntiles
    if (present(hdim2_d)) hdim2_d = 1

end subroutine get_horiz_grid_dim_d

!=======================================================================

subroutine define_cam_grids(Atm)

  ! Create grid objects on the dynamics decomposition for grids used by
  ! the dycore.  The decomposed grid object contains data for the elements
  ! in each task and information to map that data to the global grid.
  !
  ! Notes on dynamic memory management:
  !
  ! . Coordinate values and the map passed to the horiz_coord_create
  !   method are copied to the object.  The memory may be deallocated
  !   after the object is created.
  !
  ! . The area values passed to cam_grid_attribute_register are only pointed
  !   to by the attribute object, so that memory cannot be deallocated.  But the
  !   map is copied.
  !
  ! . The grid_map passed to cam_grid_register is just pointed to.
  !   Cannot be deallocated.

  use cam_grid_support,  only: horiz_coord_t, horiz_coord_create
  use cam_grid_support,  only: cam_grid_register, cam_grid_attribute_register
  use fv_grid_utils_mod, only: mid_pt_sphere
  use mpp_mod,           only: mpp_pe
  use physconst,         only: rearth

  ! arguments
  type(fv_atmos_type), target, intent(in) :: Atm(:)

  ! local variables
  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord

  integer(iMap), pointer :: grid_map(:,:)

  integer, allocatable, target, dimension(:,:) :: mygid, mygid_ew,mygid_ns
  integer                                      :: mybindex
  integer :: i, j, mapind,is,ie,js,je,isd,ied,jsd,jed,tile
  real(r8), pointer, dimension(:,:,:) :: agrid
  real(r8), pointer, dimension(:,:,:) :: grid
  real(r8), pointer, dimension(:,:)   :: area
  real(r8), pointer :: area_ffsl(:)                   !fv3 cell centered grid area in sq radians
  real(r8), pointer :: pelon_deg(:)
  real(r8), pointer :: pelat_deg(:)
  real(r8), pointer :: pelon_deg_ew(:)
  real(r8), pointer :: pelat_deg_ew(:)
  real(r8), pointer :: pelon_deg_ns(:)
  real(r8), pointer :: pelat_deg_ns(:)
  real(r8)               :: lonrad,latrad
  integer(iMap), pointer :: pemap(:)
  integer(iMap), pointer :: pemap_ew(:)
  integer(iMap), pointer :: pemap_ns(:)
  integer                :: iend, jend

  !-----------------------------------------------------------------------

  area  => Atm(mytile)%gridstruct%area_64
  agrid => Atm(mytile)%gridstruct%agrid_64
  grid => Atm(mytile)%gridstruct%grid_64
  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed
  tile = Atm(mytile)%tile

  allocate(area_ffsl((ie-is+1)*(je-js+1)))
  allocate(grid_ew(isd:ied+1,jsd:jed,2))
  allocate(grid_ns(isd:ied,jsd:jed+1,2))
  allocate(pelon_deg((ie-is+1)*(je-js+1)))
  allocate(pelon_deg_ns((ie-is+1)*(je-js+2)))
  allocate(pelon_deg_ew((ie-is+2)*(je-js+1)))
  allocate(pelat_deg((ie-is+1)*(je-js+1)))
  allocate(pelat_deg_ew((ie-is+2)*(je-js+1)))
  allocate(pelat_deg_ns((ie-is+1)*(je-js+2)))
  allocate(pemap((ie-is+1)*(je-js+1)))
  allocate(pemap_ew((ie-is+2)*(je-js+1)))
  allocate(pemap_ns((ie-is+1)*(je-js+2)))

  do j=jsd,jed
     do i=isd,ied+1
        call mid_pt_sphere(grid(i,  j,1:2), grid(i,  j+1,1:2), grid_ew(i,j,:))
     end do
  end do

  do j=jsd,jed+1
     do i=isd,ied
        call mid_pt_sphere(grid(i,j  ,1:2), grid(i+1,j  ,1:2), grid_ns(i,j,:))
     end do
  end do

  allocate(mygid(is:ie,js:je))
  allocate(mygid_ew(is:ie+1,js:je))
  allocate(mygid_ns(is:ie,js:je+1))

  mygid=0

  mybindex = mpp_pe() + 1

  do j = js, je
     do i = is, ie
        mygid(i,j)=((j-1)*(npx-1)+i)+((npx-1)*(npy-1)*(tile-1))
     end do
  end do

  !  calculate local portion of global NS index array
  !  unique global indexing bottom left to top right of each tile consecutively. Dups reported as 0
  !  North tile edges of 2,4,6 are duplicates of south edge of 3,5,1 and are reported as 0 in mygid array
  mygid_ns=0
  if (je+1 == npy) then
     jend = je+mod(tile,2)
  else
     jend = je+1
  end if
  do j = js, jend
     do i = is, ie
        mygid_ns(i,j)=(i-1)*(npy-(mod(tile-1,2))) + j + (int((tile-1)/2)*(npx-1)*(npy-1)) + (int(tile/2)*(npx-1)*(npy))
     end do
  end do
  ! appropriate tile boundaries already 0'd  need to
  ! zero inner tile je+1 boundaries (These are also repeated points between tasks in ns direction))
  if (je+1 /= npy) mygid_ns(is:ie,je+1)=0

  !  calculate local portion of global EW index array
  !  unique global indexing bottom left to top right of each tile consecutively. Dups reported as 0
  !  East tile edges of 1,3,5 are duplicates of west edge of 2,4,6 and are reported as 0 in mygid array
  mygid_ew=0
  if (ie+1 == npx) then
     iend=ie+mod(tile-1,2)
  else
     iend=ie+1
  end if
  do j = js, je
     do i = is, iend
        mygid_ew(i,j)=(j-1)*(npx-(mod(tile,2))) + i + (int(tile/2)*(npx-1)*(npy-1)) + (int((tile-1)/2)*(npx)*(npy-1))
     end do
  end do

  ! appropriate east tile boundaries already 0'd from above need to
  ! zero inner tile ie+1 boundaries on appropriate processors
  ! (These are also repeated points between tasks in ew direction)
  if (ie+1 /= npx) mygid_ew(ie+1,js:je)=0

  !-----------------------
  ! Create FFSL grid object
  !-----------------------

  ! Calculate the mapping between FFSL points and file order (tile1 thru tile6)
  mapind = 1
  do j = js, je
     do i = is, ie
        pelon_deg(mapind) = agrid(i,j,1) * rad2deg
        pelat_deg(mapind) = agrid(i,j,2) * rad2deg
        area_ffsl(mapind) = area(i,j)/(rearth*rearth)
        pemap(mapind)     = mygid(i,j)
        mapind = mapind + 1
     end do
  end do

  mapind = 1
  do j = js, je
     do i = is, ie+1
        lonrad=grid_ew(i,j,1)
        latrad=grid_ew(i,j,2)
        pelon_deg_ew(mapind) = lonrad * rad2deg
        pelat_deg_ew(mapind) = latrad * rad2deg
        pemap_ew(mapind)     = mygid_ew(i,j)
        mapind = mapind + 1
     end do
  end do

  mapind = 1
  do j = js, je+1
     do i = is, ie
        lonrad=grid_ns(i,j,1)
        latrad=grid_ns(i,j,2)
        pelon_deg_ns(mapind) = lonrad * rad2deg
        pelat_deg_ns(mapind) = latrad * rad2deg
        pemap_ns(mapind)     = mygid_ns(i,j)
        mapind = mapind + 1
     end do
  end do

  allocate(grid_map(3, (ie-is+1)*(je-js+1)))
  grid_map = 0
  mapind = 1
  do j = js, je
     do i = is, ie
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap(mapind)
        mapind = mapind + 1
     end do
  end do

  ! output local and global uniq points
  uniqpts_glob=(npx-1)*(npy-1)*6

  lat_coord => horiz_coord_create('lat_d', 'ncol_d', uniqpts_glob, 'latitude',      &
       'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
  lon_coord => horiz_coord_create('lon_d', 'ncol_d', uniqpts_glob, 'longitude',     &
       'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

  ! register dynamic A-grid, src_in(/1,2/) allows ilev,jlev,nlev ordering for restart IO
  call cam_grid_register('FFSL', dyn_decomp, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.,src_in=(/1,2/))
  call cam_grid_attribute_register('FFSL', 'cell', '', 1)
  call cam_grid_attribute_register('FFSL', 'area_d', 'FFSL grid areas', &
         'ncol_d', area_ffsl, map=pemap)

  ! register grid for writing dynamics A-Grid fields in history files 
  call cam_grid_register('FFSLHIST', dyn_decomp_hist, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSLHIST', 'cell', '', 1)
  call cam_grid_attribute_register('FFSLHIST', 'area_d', 'FFSLHIST grid areas', &
         'ncol_d', area_ffsl, map=pemap)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  ! lat_coord and lon_coord belong to grid so can't be deleted. It can be nullified
  nullify(lat_coord)
  nullify(lon_coord)
  ! area_ffsl cannot be deallocated as the attribute object is just pointing
  ! to that memory.  It can be nullified since the attribute object has
  ! the reference.
  nullify(area_ffsl)


  ! global EW uniq points
  uniqpts_glob_ew=((2*npx)-1)*(npy-1)*3

  lat_coord => horiz_coord_create('lat_d_ew', 'ncol_d_ew', uniqpts_glob_ew, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ew), pelat_deg_ew, map=pemap_ew)
  lon_coord => horiz_coord_create('lon_d_ew', 'ncol_d_ew',  uniqpts_glob_ew, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ew), pelon_deg_ew, map=pemap_ew)

  allocate(grid_map(3, (ie-is+2)*(je-js+1)))
  grid_map = 0
  mapind = 1
  do j = js, je
     do i = is, ie+1
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_ew(mapind)
        mapind = mapind + 1
     end do
  end do

  ! register dynamic D-grid, src_in(/1,2/) allows ilev,jlev,nlev ordering for restart IO
  call cam_grid_register('FFSL_EW', dyn_decomp_ew, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.,src_in=(/1,2/))
  call cam_grid_attribute_register('FFSL_EW', 'cell', '', 1)

  ! register grid for writing dynamics D-Grid fields in history files
  call cam_grid_register('FFSLHIST_EW', dyn_decomp_hist_ew, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSLHIST_EW', 'cell', '', 1)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  ! lat_coord and lon_coord belong to grid so can't be deleted. It can be nullified
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid


  ! output local and global uniq points
  uniqpts_glob_ns=((2*npy)-1)*(npx-1)*3

  lat_coord => horiz_coord_create('lat_d_ns', 'ncol_d_ns',  uniqpts_glob_ns, 'latitude',      &
       'degrees_north', 1, size(pelat_deg_ns), pelat_deg_ns, map=pemap_ns)
  lon_coord => horiz_coord_create('lon_d_ns', 'ncol_d_ns',  uniqpts_glob_ns, 'longitude',     &
       'degrees_east', 1, size(pelon_deg_ns), pelon_deg_ns, map=pemap_ns)

  allocate(grid_map(3, (ie-is+1)*(je-js+2)))
  grid_map = 0
  mapind = 1
  do j = js, je+1
     do i = is, ie
        grid_map(1, mapind) = i
        grid_map(2, mapind) = j
        grid_map(3, mapind) = pemap_ns(mapind)
        mapind = mapind + 1
     end do
  end do

  ! register dynamic D-grid, src_in(/1,2/) allows ilev,jlev,nlev ordering for restart IO
  call cam_grid_register('FFSL_NS', dyn_decomp_ns, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.,src_in=(/1,2/))
  call cam_grid_attribute_register('FFSL_NS', 'cell', '', 1)

  ! register grid for writing dynamics D-Grid fields in history files
  call cam_grid_register('FFSLHIST_NS', dyn_decomp_hist_ns, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSLHIST_NS', 'cell', '', 1)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  ! lat_coord and lon_coord belong to grid so can't be deleted. It can be nullified
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid

  deallocate(pelon_deg)
  deallocate(pelat_deg)
  deallocate(pelon_deg_ns)
  deallocate(pelat_deg_ns)
  deallocate(pelon_deg_ew)
  deallocate(pelat_deg_ew)
  deallocate(pemap)
  deallocate(pemap_ew)
  deallocate(pemap_ns)
  deallocate(mygid)
  deallocate(mygid_ew)
  deallocate(mygid_ns)

end subroutine define_cam_grids

!=========================================================================================

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)

  ! create list of attributes for the physics grid that should be copied
  ! from the corresponding grid object on the dynamics decomposition

  use cam_grid_support, only: max_hcoordname_len

  ! arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)
  !-----------------------------------------------------------------------

  gridname = 'FFSL'
  allocate(grid_attribute_names(1))
  ! For standard CAM-FV3, we need to copy the area attribute.
  ! For physgrid, the physics grid will create area
  grid_attribute_names(1) = 'cell'

end subroutine physgrid_copy_attributes_d

!=======================================================================

integer function get_dyn_grid_parm(name) result(ival)

   ! This function is in the process of being deprecated, but is still needed
   ! as a dummy interface to satisfy external references from some chemistry routines.

    use pmgrid,  only: plon, plev, plat, plevp

    character(len=*), intent(in) :: name
    integer is,ie,js,je

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    if (name == 'plat') then
        ival = plat
    else if (name == 'plon') then
        ival = (je-js+1)*(ie-is+1)
    else if (name == 'plev') then
        ival = plev
    else if (name == 'plevp') then
        ival = plevp
    else
        call endrun('get_dyn_grid_parm: undefined name: '//adjustl(trim(name)))
    end if

end function get_dyn_grid_parm

!=======================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

    ! This routine is not used for FV3, but still needed as a dummy interface to satisfy
    ! references from mo_synoz.F90 and phys_gmean.F90

    ! arguments
    character(len=*), intent(in) :: name
    real(r8), pointer :: rval(:)
    !----------------------------------------------------------------------------

    if(name == 'w') then
       call endrun('get_dyn_grid_parm_real1d: w not defined')
    else if(name == 'clat') then
       call endrun('get_dyn_grid_parm_real1d: clat not supported, use get_horiz_grid_d')
    else if(name == 'latdeg') then
       call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
    else
       nullify(rval)
    end if

end function get_dyn_grid_parm_real1d

!=========================================================================================

subroutine dyn_grid_get_colndx( igcol, ncols, owners, indx, jndx)
  use spmd_utils,       only: iam

  ! For each global column index return the owning task.  If the column is owned
  ! by this task, then also return the MPI process indicies for that column


  ! arguments
  integer, intent(in)  :: ncols
  integer, intent(in)  :: igcol(ncols)
  integer, intent(out) :: owners(ncols)
  integer, intent(out) :: indx(ncols)
  integer, intent(out) :: jndx(ncols)

  ! local variables
  integer  :: i,is,ie,js,je
  integer  :: blockid(1), bcid(1), lclblockid(1), ind(2)
  !----------------------------------------------------------------------------

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  do i = 1,ncols

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam == owners(i) ) then
        if (minval(abs(bcid(1)-mylindex)) == 0) then
           ind = minloc(abs(bcid(1)-mylindex))
           indx(i) = is+ind(1)-1
           jndx(i) = js+ind(2)-1
        end if
     else
        indx(i) = -1
        jndx(i) = -1
     endif

  end do

end subroutine dyn_grid_get_colndx

!=======================================================================

subroutine dyn_grid_get_elem_coords(ie, rlon, rlat, cdex)

   ! Returns coordinates of a specified block element of the dyn grid
   !

    ! arguments
    integer, intent(in) :: ie ! block element index
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
    integer, optional, intent(out) :: cdex(:) ! global column index
    !----------------------------------------------------------------------------

    call endrun('dyn_grid_get_elem_coords: currently not avaliable.')

end subroutine dyn_grid_get_elem_coords

!=========================================================================================

subroutine create_global(is,ie,js,je,arr_d, global_out)

  ! Gather global array of columns for the physics grid,
  ! reorder to global column order, then broadcast it to all tasks.

  use fv_mp_mod,         only: mp_gather, mp_bcst

  ! arguments
  integer,  intent(in)  :: is, ie, js, je
  real(r8), intent(in)  :: arr_d(is:ie,js:je)    ! input array
  real(r8), intent(out) :: global_out(:) ! global output in block order

  ! local variables
  integer          :: i, j, k
  integer          :: tile
  real(r8), allocatable           :: globid(:,:,:)
  real(r8), allocatable           :: globarr_tmp(:,:,:)
  !----------------------------------------------------------------------------

  tile = Atm(mytile)%tile

  if (.not. allocated(globarr_tmp)) then
     if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.(globarr_tmp)'
     allocate(globarr_tmp(npx-1, npy-1, ntiles))
  end if

  globarr_tmp(is:ie,js:je,tile)=arr_d(is:ie,js:je)
  call mp_gather(globarr_tmp, is, ie, js, je, npx-1, npy-1, ntiles)
  if (masterproc) then
     do k = 1, ntiles
        do j = 1, npy-1
           do i = 1, npx-1
              global_out(gindex_g(i,j,k)) = globarr_tmp(i,j,k)
           end do
        end do
     end do
  end if
  call mp_bcst(global_out, (npx-1)*(npy-1)*ntiles)
  deallocate(globarr_tmp)

end subroutine create_global

end module dyn_grid
