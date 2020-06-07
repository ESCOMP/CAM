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
    use dimensions_mod,   only: npx, npy
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

    integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore

    integer :: mytile = 1
    integer :: p_split  = 1
    integer, allocatable :: pelist(:)

    real(r8), parameter :: rad2deg = 180._r8/pi

    logical, allocatable :: grids_on_this_pe(:)
    type(fv_atmos_type), allocatable, target :: Atm(:)
    
    real(r8), allocatable :: block_extents_g(:,:)
    real(r8), allocatable :: tile_decomp(:,:)

public ::             &
   dyn_decomp,        &
   p_split,           &
   grids_on_this_pe,  &
   ptimelevels

!-----------------------------------------------------------------------
! Calculate Global Index

integer, allocatable, target, dimension(:,:) :: mygindex
integer, allocatable, target, dimension(:,:) :: mylindex
real(r8), allocatable, target, dimension(:,:,:) ::   locidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   blkidx_g
real(r8), allocatable, target, dimension(:,:,:) ::   gindex_g
real(r8), allocatable, target, dimension(:,:) ::   gblidx_g

integer, public :: uniqpts_glob = 0     ! number of dynamics columns
integer, public :: uniqpts_glob_ew = 0     ! number of dynamics columns for D grid ew
integer, public :: uniqpts_glob_ns = 0     ! number of dynamics columns for D grid ns

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
   use fv_mp_mod,          only: switch_current_Atm
   use hycoef,             only: hycoef_init, hyai, hybi, hypi, hypm, nprlev
   use mpp_mod,            only: mpp_init, mpp_npes, mpp_get_current_pelist,mpp_gather
   use pmgrid,             only: plev
   use ref_pres,           only: ref_pres_init
   use time_manager,       only: get_step_size
   use pio,                only: file_desc_t

   ! Local variables

   type(file_desc_t),      pointer :: fh_ini
   type (block_control_type), target   :: Atm_block

   character(len=*), parameter :: sub='dyn_grid_init'
   integer, parameter :: be_arrlen = 7
   character(len=128) :: version = '$Id$'
   character(len=128) :: tagname = '$Name$'

   real(r8) :: dt_atmos_real = 0.
   real(r8),allocatable :: rtmp(:)
   real(r8) :: block_extents(be_arrlen)
   integer, allocatable :: be_size(:)

   integer :: i, k, lat
   integer :: is,ie,js,je,npes,n

   !-----------------------------------------------------------------------
   !  from couple_main initialize atm structure - initializes fv3 grid
   !-----------------------------------------------------------------------

   !!********From Coupler_main and coupler_init.F90
   call fms_init(mpicom)
   call mpp_init()
   
   call fms_init
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

   if (Atm(mytile)%flagstruct%npz /= plev) call endrun('dyn_grid_init: FV3 dycore levels (npz) does not match model levels (plev)')

!-----------------------------------------------------------------------
!--- get block extents for each task/pe
!-----------------------------------------------------------------------
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
   block_extents(6)=mpp_pe() + 1
   block_extents(7)=((js-1)*(npx-1)+is)+((npx-1)*(npy-1)*(Atm(mytile)%tile-1))

   call mpp_gather(block_extents,be_arrlen,rtmp,be_size)
   call mp_bcst(rtmp,be_arrlen*npes)
   block_extents_g=reshape(rtmp,(/be_arrlen,npes/))

   deallocate(rtmp)
   deallocate(be_size)

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

end subroutine dyn_grid_init

!=======================================================================

subroutine get_block_bounds_d(block_first, block_last)

  use spmd_utils,    only : npes

    implicit none

    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks

    block_first = 1
    block_last = npes

end subroutine get_block_bounds_d

!=======================================================================

subroutine get_block_gcol_d(blockid, size, cdex)

    implicit none

    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size
    integer, intent(out):: cdex(size)   ! global column indices

    integer :: i, j, n,is,ie,js,je,tile

    is=block_extents_g(1,blockid)
    ie=block_extents_g(2,blockid)
    js=block_extents_g(3,blockid)
    je=block_extents_g(4,blockid)
    tile=block_extents_g(5,blockid)

    if (size .ne. (ie - is + 1) * (je - js + 1)) then
        call endrun ('get_block_gcol_d: block sizes are not consistent.')
    end if

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

    implicit none

    integer, intent(in) :: blockid

    get_block_gcol_cnt_d=count(blkidx_g.eq.blockid)

end function get_block_gcol_cnt_d

!=======================================================================

integer function get_block_lvl_cnt_d(blockid, bcid)

    use pmgrid, only: plevp

    implicit none

    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid     ! column index within block

    get_block_lvl_cnt_d = plevp

end function get_block_lvl_cnt_d

!=======================================================================

subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

    use pmgrid, only: plev

    implicit none

    integer, intent(in)  :: blockid        ! global block id
    integer, intent(in)  :: bcid           ! column index within block
    integer, intent(in)  :: lvlsiz         ! dimension of levels array
    integer, intent(out) :: levels(lvlsiz) ! levels indices for block

    integer :: k

    if (lvlsiz < plev + 1) then
        write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (',lvlsiz,' < ',plev + 1,')'
        call endrun
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

    implicit none

    integer, intent(in) :: blockid  ! global block id

    get_block_owner_d = blockid - 1

end function get_block_owner_d

!=======================================================================

subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)
 
    use dimensions_mod,     only: npx, npy
    use spmd_utils,        only: iam
    implicit none
    integer, intent(in)  :: gcol         ! global column index
    integer, intent(in)  :: cnt          ! size of blockid and bcid arrays
    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)

    integer :: tot

    if (cnt .ne. 1) then
       call endrun ('get_gcol_block_d: cnt is not equal to 1:.')
    end if
    tot=(npx-1)*(npy-1)*6
    if (gcol.lt.1.or.gcol.gt.tot) then
       call endrun ('get_gcol_block_d: global column number is out of bounds')
    else
       blockid(1) = gblidx_g(gcol,1)
       bcid(1) = gblidx_g(gcol,2)
    end if

    if (present(localblockid)) then
        localblockid(cnt) = 1
    end if

end subroutine get_gcol_block_d

!=======================================================================

integer function get_gcol_block_cnt_d(gcol)

    implicit none

    integer, intent(in) :: gcol     ! global column index

    get_gcol_block_cnt_d = 1

end function get_gcol_block_cnt_d

!=======================================================================

subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

    use fv_mp_mod,         only: mp_gather, mp_bcst
    implicit none

    integer, intent(in)   :: nxy                       ! array sizes
    real(r8), intent(out), optional :: clat_d_out(nxy) ! column latitudes
    real(r8), intent(out), optional :: clon_d_out(nxy) ! column longitudes
    real(r8), intent(out), optional :: area_d_out(nxy) ! column surface area
    real(r8), intent(out), optional :: wght_d_out(nxy) ! column integration
    real(r8), intent(out), optional :: lat_d_out(nxy)  ! column degree latitudes
    real(r8), intent(out), optional :: lon_d_out(nxy)  ! column degree longitudes

    integer :: i, j, k

    real(r8) :: clon_d_mpi(nxy)
    real(r8) :: clat_d_mpi(nxy)
    real(r8) :: area_d_mpi(nxy)
    real(r8) :: wght_d_mpi(nxy)

    real(r8), allocatable :: clon_d(:,:,:)
    real(r8), allocatable :: clat_d(:,:,:)
    real(r8), allocatable :: area_d(:,:,:)
    real(r8), allocatable :: indx_d(:,:,:)

    real(r8), pointer, dimension(:,:,:) :: agrid
    real(r8), pointer, dimension(:,:)   :: area
    integer,  pointer                   :: ntiles_g,tile
    integer                             :: is,ie,js,je,ind
!
    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    area  => Atm(mytile)%gridstruct%area_64
    agrid => Atm(mytile)%gridstruct%agrid_64
    ntiles_g => Atm(mytile)%gridstruct%ntiles_g
    tile=>  Atm(mytile)%tile

    clon_d_mpi = -999.0_r8
    clat_d_mpi = -999.0_r8
    area_d_mpi = -999.0_r8
    wght_d_mpi = -999.0_r8

    allocate(clon_d(npx-1, npy-1, ntiles_g))
    allocate(clat_d(npx-1, npy-1, ntiles_g))
    allocate(area_d(npx-1, npy-1, ntiles_g))
    allocate(indx_d(npx-1, npy-1, ntiles_g))

    do j = js, je
        do i = is, ie
            clon_d(i,j,tile) = agrid (i,j,1)
            clat_d(i,j,tile) = agrid (i,j,2)
            area_d(i,j,tile) = area  (i,j) / (rearth * rearth)
            indx_d(i,j,tile) = mygindex  (i,j)
        end do
    end do

    call mp_gather(indx_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
    if (masterproc) then
        if (maxval(indx_d) .gt. nxy .or. minval(indx_d) .lt. 1) then
            write(iulog,*)'GET_HORIZ_GRID_D: indx_d is out of bound. nxy = ',nxy
            call endrun
        end if
    end if

    if (present(clat_d_out)) then
        call mp_gather(clat_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       clat_d_mpi(ind) = clat_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clat_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clat_d_mpi .eq. -999.0_r8)) then
            clat_d_out = clat_d_mpi
        else
            call endrun('clat_d_mpi is not full filled')
        end if
    end if

    if (present(clon_d_out)) then
        call mp_gather(clon_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       clon_d_mpi(ind) = clon_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clon_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clon_d_mpi .eq. -999.0_r8)) then
            clon_d_out = clon_d_mpi
        else
            call endrun('clon_d_mpi is not full filled')
        end if
    end if

    if (present(lat_d_out)) then
        call mp_gather(clat_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       clat_d_mpi(ind) = clat_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clat_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clat_d_mpi .eq. -999.0_r8)) then
            lat_d_out = clat_d_mpi * rad2deg
        else
            call endrun('clat_d_mpi is not full filled')
        end if
    end if

    if (present(lon_d_out)) then
        call mp_gather(clon_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       clon_d_mpi(ind) = clon_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(clon_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(clon_d_mpi .eq. -999.0_r8)) then
            lon_d_out = clon_d_mpi * rad2deg
        else
            call endrun('clon_d_mpi is not full filled')
        end if
    end if

    if (present(area_d_out)) then
        call mp_gather(area_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       area_d_mpi(ind) = area_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(area_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(area_d_mpi .eq. -999.0_r8)) then
            area_d_out = area_d_mpi
        else
            call endrun('area_d_mpi is not full filled')
        end if
    end if

    if (present(wght_d_out)) then
        call mp_gather(area_d, is, ie, js, je, npx-1, npy-1, ntiles_g)
        if (masterproc) then
            do k = 1, ntiles_g
                do j = 1, npy-1
                    do i = 1, npx-1
                       ind=indx_d(i,j,k)
                       wght_d_mpi(ind) = area_d(i,j,k)
                    end do
                end do
            end do
        end if
        call mp_bcst(wght_d_mpi, (npx-1)*(npy-1)*ntiles_g)
        if (.not. any(wght_d_mpi .eq. -999.0_r8)) then
            wght_d_out = wght_d_mpi
        else
            call endrun('wght_d_mpi is not full filled')
        end if
    end if

    deallocate(clon_d)
    deallocate(clat_d)
    deallocate(area_d)
    deallocate(indx_d)

end subroutine get_horiz_grid_d

!=======================================================================

subroutine get_horiz_grid_d1(nxy, clat_d_out, clon_d_out, area_d_out, wght_d_out, lat_d_out, lon_d_out)

  implicit none

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

  area  => Atm(mytile)%gridstruct%area_64
  agrid => Atm(mytile)%gridstruct%agrid_64
  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  if (present(clon_d_out)) then
     if (size(clon_d_out) /= nxy) call endrun(sub//': bad clon_d_out array size')
     call create_global(agrid(is:ie,js:je,1), clon_d_out)
  end if
  if (present(clat_d_out)) then
     if (size(clat_d_out) /= nxy) call endrun(sub//': bad clat_d_out array size')
     call create_global(agrid(is:ie,js:je,2), clat_d_out)
  end if
  if (present(area_d_out).or.present(wght_d_out)) then
     allocate(tmparr(is:ie,js:je))
     tmparr(is:ie,js:je) = area  (is:ie,js:je) / (rearth * rearth)
     if (present(area_d_out)) then
        if (size(area_d_out) /= nxy) call endrun(sub//': bad area_d_out array size')
        call create_global(tmparr, area_d_out)
     end if
     if (present(wght_d_out)) then
        if (size(wght_d_out) /= nxy) call endrun(sub//': bad wght_d_out array size')
        call create_global(tmparr, wght_d_out)
     end if
     deallocate(tmparr)
  end if
  if (present(lon_d_out)) then
     if (size(lon_d_out) /= nxy) call endrun(sub//': bad clon_d_out array size')
     call create_global(agrid(is:ie,js:je,1), lon_d_out)
     lon_d_out=lon_d_out*rad2deg
  end if
  if (present(lat_d_out)) then
     if (size(lat_d_out) /= nxy) call endrun(sub//': bad clat_d_out array size')
     call create_global(agrid(is:ie,js:je,2), lat_d_out)
     lat_d_out=lat_d_out*rad2deg
  end if

end subroutine get_horiz_grid_d1

!=======================================================================

subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)
    use spmd_utils,    only : npes
    implicit none

    integer, intent(out) :: hdim1_d                 ! first horizontal dimension
    integer, intent(out), optional :: hdim2_d       ! second horizontal dimension
    integer is,ie,js,je

    hdim1_d = count(gindex_g.gt.0)
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
  use fv_mp_mod,         only: mp_gather, mp_bcst
  use spmd_utils,        only: iam
  use physconst,         only: rearth
  implicit none
  
  type(fv_atmos_type), target, intent(in) :: Atm(:)


  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  
  integer(iMap), pointer :: grid_map(:,:)

  integer, allocatable, target, dimension(:,:) :: mygindex_ew,mygindex_ns
  integer                                      :: mybindex
  integer :: n, i, j, mapind,is,ie,js,je,isd,ied,jsd,jed,tile,nregions
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
  integer :: nx,ny,ind

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
  nregions = Atm(mytile)%gridstruct%ntiles_g

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

  allocate(locidx_g(npx-1,npy-1,nregions))
  allocate(blkidx_g(npx-1,npy-1,nregions))
  allocate(gindex_g(npx-1,npy-1,nregions))
  allocate(gblidx_g( (npx-1)*(npy-1)*nregions,2))
  allocate(mygindex(is:ie,js:je))
  allocate(mygindex_ew(is:ie+1,js:je))
  allocate(mygindex_ns(is:ie,js:je+1))
  allocate(mylindex(is:ie,js:je))

  !  calculate local portion of global A-Grid index array
  !  unique global indexing bottom left to top right of each tile consecutively. Dups reported as 0
  mygindex=0
  mylindex=0
 
  mybindex = mpp_pe() + 1

  mapind = 1
  do j = js, je
     do i = is, ie
        mygindex(i,j)=((j-1)*(npx-1)+i)+((npx-1)*(npy-1)*(tile-1))
        mylindex(i,j)=mapind
        mapind = mapind + 1
     end do
  end do

  ! output local and global uniq points
  uniqpts_glob=(npx-1)*(npy-1)*6

  !  calculate local portion of global NS index array
  !  unique global indexing bottom left to top right of each tile consecutively. Dups reported as 0
  !  North tile edges of 2,4,6 are duplicates of south edge of 3,5,1 and are reported as 0 in mygindex array
  mygindex_ns=0
  if (je+1.eq.npy) then
     do j = js, je+mod(tile,2)
        do i = is, ie
           mygindex_ns(i,j)=(i-1)*(npy-(mod(tile-1,2))) + j + (int((tile-1)/2)*(npx-1)*(npy-1)) + (int(tile/2)*(npx-1)*(npy))
        end do
     end do
  else
     do j = js, je+1
        do i = is, ie
           mygindex_ns(i,j)=(i-1)*(npy-(mod(tile-1,2))) + j + (int((tile-1)/2)*(npx-1)*(npy-1)) + (int(tile/2)*(npx-1)*(npy))
        end do
     end do
  end if
  ! appropriate tile boundaries already 0'd  need to
  ! zero inner tile je+1 boundaries (These are also repeated points between tasks in ns direction))
  if (je+1.ne.npy) mygindex_ns(is:ie,je+1)=0

  ! output local and global uniq points
  uniqpts_glob_ns=((2*npy)-1)*(npx-1)*3


  !  calculate local portion of global EW index array
  !  unique global indexing bottom left to top right of each tile consecutively. Dups reported as 0
  !  East tile edges of 1,3,5 are duplicates of west edge of 2,4,6 and are reported as 0 in mygindex array
  mygindex_ew=0
  do j = js, je
     if (ie+1.eq.npx) then
        do i = is, ie+mod(tile-1,2)
           mygindex_ew(i,j)=(j-1)*(npx-(mod(tile,2))) + i + (int(tile/2)*(npx-1)*(npy-1)) + (int((tile-1)/2)*(npx)*(npy-1))
        end do
     else
        do i = is, ie+1
           mygindex_ew(i,j)=(j-1)*(npx-(mod(tile,2))) + i + (int(tile/2)*(npx-1)*(npy-1)) + (int((tile-1)/2)*(npx)*(npy-1))
        end do
     end if
  end do

  ! appropriate east tile boundaries already 0'd from above need to
  ! zero inner tile ie+1 boundaries on appropriate processors
  ! (These are also repeated points between tasks in ew direction)
  if (ie+1.ne.npx) mygindex_ew(ie+1,js:je)=0

  ! global EW uniq points
  uniqpts_glob_ew=((2*npx)-1)*(npy-1)*3

  blkidx_g(is:ie,js:je,tile)=mybindex
  locidx_g(is:ie,js:je,tile)=mylindex(is:ie,js:je)
  gindex_g(is:ie,js:je,tile)=mygindex(is:ie,js:je)

  nx=npx-1
  ny=npy-1
  call mp_gather(locidx_g, is, ie, js, je, nx ,ny, nregions)
  call mp_gather(blkidx_g, is, ie, js, je, nx ,ny, nregions)
  call mp_gather(gindex_g, is, ie, js, je, nx, ny, nregions)
  call mp_bcst(locidx_g, nx, ny, nregions)
  call mp_bcst(blkidx_g, nx, ny, nregions)
  call mp_bcst(gindex_g, nx, ny, nregions)

  do n = 1, nregions
     do j=1,nx
        do i=1,ny
           ind=gindex_g(i,j,n)
           gblidx_g(ind,1)=blkidx_g(i,j,n)
           gblidx_g(ind,2)=locidx_g(i,j,n)
        end do
     end do
  end do

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
        pemap(mapind)     = mygindex(i,j)
        mapind = mapind + 1
     end do
  end do

  pelon_deg_ew = 0
  pelat_deg_ew = 0
  mapind = 1
  do j = js, je
     do i = is, ie+1
        lonrad=grid_ew(i,j,1)
        latrad=grid_ew(i,j,2)
        pelon_deg_ew(mapind) = lonrad * rad2deg
        pelat_deg_ew(mapind) = latrad * rad2deg
        pemap_ew(mapind)     = mygindex_ew(i,j)
        mapind = mapind + 1
     end do
  end do
  pelon_deg_ns = 0
  pelat_deg_ns = 0
  mapind = 1
  do j = js, je+1
     do i = is, ie
        lonrad=grid_ns(i,j,1)
        latrad=grid_ns(i,j,2)
        pelon_deg_ns(mapind) = lonrad * rad2deg
        pelat_deg_ns(mapind) = latrad * rad2deg
        pemap_ns(mapind)     = mygindex_ns(i,j)
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

  lat_coord => horiz_coord_create('lat_d', 'ncol_d', uniqpts_glob, 'latitude',      &
       'degrees_north', 1, size(pelat_deg), pelat_deg, map=pemap)
  lon_coord => horiz_coord_create('lon_d', 'ncol_d', uniqpts_glob, 'longitude',     &
       'degrees_east', 1, size(pelon_deg), pelon_deg, map=pemap)

  call cam_grid_register('FFSL', dyn_decomp, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL', 'cell', '', 1)
  call cam_grid_attribute_register('FFSL', 'area_d', 'FFSL grid areas', &
         'ncol_d', area_ffsl, map=pemap)         
  nullify(grid_map)
  nullify(lat_coord)
  nullify(lon_coord)
  nullify(area_ffsl)

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

  call cam_grid_register('FFSL_EW', dyn_decomp_ew, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_EW', 'cell', '', 1)
  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
  nullify(lat_coord)         ! Belongs to grid
  nullify(lon_coord)         ! Belongs to grid

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

  call cam_grid_register('FFSL_NS', dyn_decomp_ns, lat_coord, lon_coord,          &
       grid_map, block_indexed=.false., unstruct=.true.)
  call cam_grid_attribute_register('FFSL_NS', 'cell', '', 1)

  ! grid_map cannot be deallocated as the cam_filemap_t object just points
  ! to it.  It can be nullified.
  nullify(grid_map)
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
  deallocate(mygindex_ew)
  deallocate(mygindex_ns)

end subroutine define_cam_grids

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
  use cam_grid_support, only: max_hcoordname_len

  ! Dummy arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

  gridname = 'FFSL'
  allocate(grid_attribute_names(1))
  ! For standard CAM-FV3, we need to copy the area attribute.
  ! For physgrid, the physics grid will create area
  grid_attribute_names(1) = 'cell'

end subroutine physgrid_copy_attributes_d

!=======================================================================

integer function get_dyn_grid_parm(name) result(ival)

    use pmgrid,  only: plon, plev, plat, plevp

    implicit none

    character(len=*), intent(in) :: name
    integer is,ie,js,je

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je

    if(name.eq.'ne') then
        ival = -1
    else if (name.eq.'np') then
        ival = -1
    else if (name .eq. 'plat') then
        ival = plat
    else if (name .eq. 'plon') then
        ival = (je-js+1)*(ie-is+1)
    else if (name .eq. 'plev') then
        ival = plev
    else if (name .eq. 'plevp') then
        ival = plevp
    else if (name .eq. 'beglat') then
        ival = 0
    else if (name .eq. 'endlat') then
        ival = 0
    else if (name .eq. 'beglon') then
        ival = 0
    else if (name .eq. 'endlon') then
        ival = 0
    else if (name .eq. 'beglonxy') then
        ival = 1
    else if (name .eq. 'endlonxy') then
        ival = (ie - is + 1) * (je - js + 1)  ! number of column per block
    else if (name .eq. 'beglatxy') then
        ival = 1
    else if (name .eq. 'endlatxy') then
        ival = 1                              ! number of block per MPI task
    else if (name .eq. 'splon') then
        ival = -1
    else
        call endrun('get_dyn_grid_parm: undefined name: '//adjustl(trim(name)))
    end if

end function get_dyn_grid_parm

!=======================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

    implicit none

    character(len=*), intent(in) :: name
    real(r8), pointer :: rval(:)

    if(name.eq.'w') then
       call endrun('get_dyn_grid_parm_real1d: w not defined')
    else if(name.eq.'clat') then
       call endrun('get_dyn_grid_parm_real1d: clat not supported, use get_horiz_grid_d')
    else if(name.eq.'latdeg') then
       call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
    else
       nullify(rval)
    end if

end function get_dyn_grid_parm_real1d

!#######################################################################
subroutine dyn_grid_get_colndx( igcol, ncols, owners, indx, jndx)
  use spmd_utils,       only: iam
   ! For each global column index return the owning task.  If the column is owned
   ! by this task, then also return the MPI process indicies for that column

  integer, intent(in)  :: ncols
  integer, intent(in)  :: igcol(ncols)
  integer, intent(out) :: owners(ncols)
  integer, intent(out) :: indx(ncols)
  integer, intent(out) :: jndx(ncols)

  integer  :: i,is,ie,js,je
  integer  :: blockid(1), bcid(1), lclblockid(1), ind(2)

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  do i = 1,ncols

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        if (minval(abs(bcid(1)-mylindex)) .eq. 0) then
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

    implicit none

    integer, intent(in) :: ie ! block element index
    real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
    real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
    integer, optional, intent(out) :: cdex(:) ! global column index

    call endrun('dyn_grid_get_elem_coords: currently not avaliable.')

end subroutine dyn_grid_get_elem_coords

!=========================================================================================

subroutine create_global(arr_d, global_out)

  use fv_mp_mod,         only: mp_gather, mp_bcst
  implicit none
  
  real(r8), intent(in)  :: arr_d(:,:)    ! input array
  real(r8), intent(out) :: global_out(:) ! global output in block order
  
  ! local variables
  integer :: i, j, k, is, ie, js, je
  integer, pointer :: ntiles_g
  integer, pointer :: tile
  real(r8), allocatable     :: globid(:,:,:)
  real(r8), allocatable     :: globarr_by_tile(:,:,:)
  
  ntiles_g => Atm(mytile)%gridstruct%ntiles_g
  tile=>  Atm(mytile)%tile
  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  
  if (.not. allocated(globid)) then
     if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.'
     allocate(globid(npx-1, npy-1, ntiles_g))
     globid(is:ie,js:je,tile) = mygindex(is:ie,js:je)
     call mp_gather(globid, is, ie, js, je, npx-1, npy-1, ntiles_g)
   end if

  if (.not. allocated(globarr_by_tile)) then
     if (masterproc) write(iulog, *) 'INFO: Non-scalable action: Allocating global blocks in FV3 dycore.'
     allocate(globarr_by_tile(npx-1, npy-1, ntiles_g))
  end if

  globarr_by_tile(is:ie,js:je,tile)=arr_d(is:ie,js:je)
  call mp_gather(globarr_by_tile, is, ie, js, je, npx-1, npy-1, ntiles_g)
  if (masterproc) then
     do k = 1, ntiles_g
        do j = 1, npy-1
           do i = 1, npx-1
              global_out(globid(i,j,k)) = globarr_by_tile(i,j,k)
           end do
        end do
     end do
  end if
  call mp_bcst(global_out, (npx-1)*(npy-1)*ntiles_g)
end subroutine create_global

end module dyn_grid
