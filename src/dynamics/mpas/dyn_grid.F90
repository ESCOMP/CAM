#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module dyn_grid

!-------------------------------------------------------------------------------
!
! Define MPAS computational grids on the dynamics decomposition.
!
!
! Module responsibilities:
!
! . Provide the physics/dynamics coupler (in module phys_grid) with data for the
!   physics grid on the dynamics decomposition.
!
! . Create CAM grid objects that are used by the I/O functionality to read
!   data from an unstructured grid format to the dynamics data structures, and
!   to write from the dynamics data structures to unstructured grid format.  The
!   global column ordering for the unstructured grid is determined by the dycore.
!
!-------------------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8
use spmd_utils,        only: iam, masterproc, mpicom, npes

use pmgrid,            only: plev, plevp
use physconst,         only: pi, rearth

use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun

use pio,               only: file_desc_t, var_desc_t, &
                             pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                             pio_double, pio_def_dim, pio_def_var, &
                             pio_put_var, pio_get_var, &
                             pio_seterrorhandling, PIO_BCAST_ERROR, PIO_NOERR

use cam_mpas_subdriver, only : domain_ptr

implicit none
private
save

integer, parameter :: dyn_decomp  = 101 ! cell center grid
integer, parameter :: ptimelevels = 2

integer ::      &
   maxNCells         ! maximum number of cells for any task

public :: &
   dyn_decomp, &
   ptimelevels, &
   maxNCells,   &
   dyn_grid_init, &
   get_block_bounds_d, &
   get_block_gcol_d, &
   get_block_gcol_cnt_d, &
   get_block_lvl_cnt_d, &
   get_block_levels_d, &
   get_block_owner_d, &
   get_gcol_block_d, &
   get_gcol_block_cnt_d, &
   get_horiz_grid_dim_d, &
   get_horiz_grid_d, &
   get_dyn_grid_parm, &
   get_dyn_grid_parm_real1d, &
   dyn_grid_get_elem_coords, &
   dyn_grid_get_colndx, &
   physgrid_copy_attributes_d

real(r8), parameter :: rad2deg=180.0_r8/pi ! convert radians to degrees

! vertical reference heights (m)
real(r8) :: zw(plevp), zw_mid(plev)

integer ::      &
   nCells,      &    ! global number of cells/columns
   nEdges,      &    ! global number of edges
   nVertices,   &    ! global number of vertices
   maxEdges,    &    ! maximum number of edges per cell
   nVertLevels       ! number of vertical layers (yes, layers and not layer interfaces...)


! global grid data

integer, public, allocatable :: col_indices_in_block(:,:)  !  global column indices (used in dp_coupling)
integer,         allocatable :: num_col_per_block(:)
integer,         allocatable :: global_blockid(:)
integer, public, allocatable :: local_col_index(:)  !  local to block

real(r8), dimension(:), pointer :: lonCell               ! global cell longitudes
real(r8), dimension(:), pointer :: latCell               ! global cell latitudes
real(r8), dimension(:), pointer :: areaCell              ! global cell areas

integer, pointer, public :: nCellsSolve, &
                            nEdgesSolve, &
                            nVerticesSolve, &
                            nVertLevelsSolve

!=========================================================================================
contains
!=========================================================================================


!-----------------------------------------------------------------------
!  routine dyn_grid_init
!
!> \brief Initialize dynamics grid
!> \details
!>  Prepares module variables describing the dynamics grid on the dynamics
!>  decomposition for later use by interface routines get_gcol_block_d, etc.
!
!-----------------------------------------------------------------------
subroutine dyn_grid_init()

   use ref_pres,            only: std_atm_pres, ref_pres_init
   use time_manager,        only: get_step_size

   use cam_initfiles,       only: initial_file_get_id

   use cam_history_support, only: add_vert_coord

   use cam_mpas_subdriver,  only: cam_mpas_init_phase3, cam_mpas_read_geometry,         &
                                  cam_mpas_get_global_dims, cam_mpas_get_global_coords, &
                                  cam_mpas_get_global_blocks
        
   type(file_desc_t), pointer :: fh_ini

   integer  :: k
   integer  :: num_pr_lev       ! number of top levels using pure pressure representation
   real(r8) :: pref_edge(plevp) ! reference pressure at layer edges (Pa)
   real(r8) :: pref_mid(plev)   ! reference pressure at layer midpoints (Pa)

   character(len=*), parameter :: subname = 'dyn_grid::dyn_grid_init'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(0, 'begin '//subname)


   ! Get filehandle for initial file
   fh_ini => initial_file_get_id()

   call cam_mpas_init_phase3(fh_ini, endrun)

   ! Read reference heights.
   call ref_height_read(fh_ini)

   ! Compute reference pressures from reference heights.
   call std_atm_pres(zw, pref_edge)
   pref_mid = (pref_edge(1:plev) + pref_edge(2:plevp)) * 0.5_r8

   num_pr_lev = 0
   call ref_pres_init(pref_edge, pref_mid, num_pr_lev)

   ! Vertical coordinates for output streams
   call add_vert_coord('lev', plev,                                         &
         'zeta level at vertical midpoints', 'm', zw_mid)
   call add_vert_coord('ilev', plevp,                                       &
         'zeta level at vertical interfaces', 'm', zw)

   if (masterproc) then
      write(iulog,'(a)')' Reference Layer Locations: '
      write(iulog,'(a)')' index      height (m)              pressure (hPa) '
      do k= 1, plev
         write(iulog,9830) k, zw(k), pref_edge(k)
         write(iulog,9840)    zw_mid(k), pref_mid(k)
      end do
      write(iulog,9830) plevp, zw(plevp), pref_edge(plevp)
   end if

   ! Read distributed latitude, longitude, and area of mesh cells,
   ! with the resulting fields stored in an MPAS pool
   call cam_mpas_read_geometry(fh_ini, endrun)

   ! Query global grid dimensions from MPAS
   call cam_mpas_get_global_dims(nCells, nEdges, nVertices, maxEdges, nVertLevels, maxNCells)

   ! Temporary global arrays needed by phys_grid_init
   allocate(lonCell(nCells))
   allocate(latCell(nCells))
   allocate(areaCell(nCells))
   call cam_mpas_get_global_coords(latCell, lonCell, areaCell)
   
   allocate(num_col_per_block(npes))
   allocate(col_indices_in_block(maxNCells,npes))
   allocate(global_blockid(nCells))
   allocate(local_col_index(nCells))
   call cam_mpas_get_global_blocks(num_col_per_block, col_indices_in_block, global_blockID, local_col_index)
   
   ! Define the dynamics and physics grids on the dynamics decompostion.
   ! Physics grid on the physics decomposition is defined in phys_grid_init.
   call define_cam_grids()
   
9830 format(1x, i3, f15.4, 9x, f15.4)
9840 format(1x, 3x, 12x, f15.4, 9x, f15.4)

end subroutine dyn_grid_init


!-----------------------------------------------------------------------
!  routine get_block_bounds_d
!
!> \brief Return first and last indices used in global block ordering
!> \details
!>  Return first and last indices used in global block ordering.
!>
!>  Are these the blocks owned by the calling MPI task, or all blocks
!>  across all MPI tasks?
!>
!>  Is this a 1-based ordering? Why would the first block ever not be 1?
!
!-----------------------------------------------------------------------
subroutine get_block_bounds_d( block_first, block_last )

   integer, intent(out) :: block_first  ! first global index used for blocks
   integer, intent(out) :: block_last   ! last global index used for blocks

   character(len=*), parameter :: subname = 'dyn_grid::get_block_bounds_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   block_first = 1
   block_last = npes

end subroutine get_block_bounds_d


!-----------------------------------------------------------------------
!  routine get_block_gcol_cnt_d
!
!> \brief Return number of dynamics columns in a block
!> \details
!>  Returns the number of dynamics columns in the block with the specified
!>  global block ID.
!>
!>  Will the blockid values be only among those block IDs owned by the calling
!>  MPI task, or can the blockid be the ID of a block owned by any MPI task?
!
!-----------------------------------------------------------------------
integer function get_block_gcol_cnt_d( blockid )

   integer, intent(in) :: blockid

   character(len=*), parameter :: subname = 'dyn_grid::get_block_gcol_cnt_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   get_block_gcol_cnt_d = num_col_per_block(blockid)

end function get_block_gcol_cnt_d


!-----------------------------------------------------------------------
!  routine get_block_gcol_d
!
!> \brief Return list of dynamics column indices in a block
!> \details
!>  Return list of global dynamics column indices in the block with
!>  the specified global block ID.
!>
!>  Will the blockid values be among only those block IDs owned by the calling
!>  MPI task, or can the blockid be the ID of a block owned by any MPI task?
!
!-----------------------------------------------------------------------
subroutine get_block_gcol_d(blockid, asize, cdex)

   integer, intent(in) :: blockid      ! global block id
   integer, intent(in) :: asize        ! array size

   integer, intent(out):: cdex(asize)  ! global column indices

   integer :: icol

   character(len=*), parameter :: subname = 'dyn_grid::get_block_gcol_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   do icol = 1, num_col_per_block(blockid)
      cdex(icol) = col_indices_in_block(icol, blockid)
   end do
   do icol = num_col_per_block(blockid)+1, asize
      cdex(icol) = 0
   end do

end subroutine get_block_gcol_d
   
   
!-----------------------------------------------------------------------
!  routine get_block_lvl_cnt_d
!
!> \brief Return number of levels in a column
!> \details
!>  Returns the number of levels in the specified column of the specified block.
!>  If column includes surface fields, then it is defined to also
!>  include level 0.
!>
!>  Will the blockid values be among only those block IDs owned by the calling
!>  MPI task, or can the blockid be the ID of a block owned by any MPI task?
!>
!>  Is bcid a global column ID, or a column ID in the local index space of
!>  the block?
!
!-----------------------------------------------------------------------
integer function get_block_lvl_cnt_d(blockid, bcid)

   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid     ! column index within block

   character(len=*), parameter :: subname = 'dyn_grid::get_block_lvl_cnt_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   ! At present, all blocks have the same number of levels
   get_block_lvl_cnt_d = plevp

end function get_block_lvl_cnt_d


!-----------------------------------------------------------------------
!  routine get_block_levels_d
!
!> \brief Return level indices in a column
!> \details
!>  Returns the level indices in the column of the specified global block.
!>
!>  Will the blockid values be among only those block IDs owned by the calling
!>  MPI task, or can the blockid be the ID of a block owned by any MPI task?
!>
!>  Is bcid a global column ID, or a column ID in the local index space of
!>  the block?
!
!-----------------------------------------------------------------------
subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

   integer, intent(in) :: blockid  ! global block id
   integer, intent(in) :: bcid    ! column index within block
   integer, intent(in) :: lvlsiz   ! dimension of levels array

   integer, intent(out) :: levels(lvlsiz) ! levels indices for block

   integer :: k
   character(len=128) :: errmsg

   character(len=*), parameter :: subname = 'dyn_grid::get_block_levels_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   if ( lvlsiz < plev + 1 ) then
      write(errmsg,*) ': levels array not large enough (', lvlsiz,' < ',plev + 1,')'
      call endrun( subname // trim(errmsg) )
   else
      do k = 0, plev
         levels(k+1) = k
      end do
      do k = plev+2, lvlsiz
         levels(k) = -1
      end do
   end if

end subroutine get_block_levels_d


!-----------------------------------------------------------------------
!  routine get_gcol_block_cnt_d
!
!> \brief Return number of blocks containing data for a column
!> \details
!>  Return number of blocks containing data for the vertical column
!>  with the specified global column index.
!>
!>  Under which conditions would the returned value potentially differ from 1?
!
!-----------------------------------------------------------------------
integer function get_gcol_block_cnt_d(gcol)

   integer, intent(in) :: gcol     ! global column index

   character(len=*), parameter :: subname = 'dyn_grid::get_gcol_block_cnt_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   ! For MPAS, each column is contained by exactly one block
   get_gcol_block_cnt_d = 1

end function get_gcol_block_cnt_d


!-----------------------------------------------------------------------
!  routine get_gcol_block_d
!
!> \brief Return global block index and local column index for a global column
!> \details
!>  Return global block index and local column index for a global column index.
!>
!>  Can this routine be called for global columns that are not owned by
!>  the calling task?
!
!-----------------------------------------------------------------------
subroutine get_gcol_block_d(gcol, cnt, blockid, bcid, localblockid)

   integer, intent(in) :: gcol     ! global column index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
   integer, intent(out), optional :: localblockid(cnt)

   integer :: j

   character(len=*), parameter :: subname = 'dyn_grid::get_gcol_block_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   if ( cnt < 1 ) then
      call endrun( subname // ':: arrays not large enough' )
   end if

   blockid(1) = global_blockid(gcol)
   bcid(1) = local_col_index(gcol)

   do j=2,cnt
      blockid(j) = -1
      bcid(j)    = -1
   end do

end subroutine get_gcol_block_d


!-----------------------------------------------------------------------
!  routine get_block_owner_d
!
!> \brief Return ID of task that owns a block
!> \details
!>  Returns the ID of the task that owns the indicated global block.
!>
!>  Should we assume that task IDs are 0-based (as in MPI)?
!
!-----------------------------------------------------------------------
integer function get_block_owner_d(blockid)

   integer, intent(in) :: blockid  ! global block id

   character(len=*), parameter :: subname = 'dyn_grid::get_block_owner_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   get_block_owner_d = (blockid - 1)

end function get_block_owner_d


!-----------------------------------------------------------------------
!  routine get_horiz_grid_dim_d
!
!> \brief Return declared horizontal dimensions of computational grid
!> \details
!>  Return declared horizontal dimensions of computational grid.
!>  For non-lon/lat grids, declare grid to be one-dimensional,
!>  i.e., (ncols x 1).
!>
!>  Is this the global dimension, or the task-local dimension?
!>  I.e., is the number of cells the total number of cells in the mesh
!>  or just the number of cells in blocks owned by this MPI task?
!
!-----------------------------------------------------------------------
subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)

   integer, intent(out) :: hdim1_d             ! first horizontal dimension
   integer, intent(out), optional :: hdim2_d   ! second horizontal dimension

   character(len=*), parameter :: subname = 'dyn_grid::get_horiz_grid_dim_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   hdim1_d = nCells

   if( present(hdim2_d) ) hdim2_d = 1

end subroutine get_horiz_grid_dim_d


!-----------------------------------------------------------------------
!  routine get_horiz_grid_d
!
!> \brief Returns lat, lon, area, and interp weight for a global column
!> \details
!>  Return latitude and longitude (in radians), column surface
!>  area (in radians squared) and surface integration weights
!>  for global column indices that will be passed to/from physics.
!>
!>  Can this routine be passed global column IDs for columns that are not
!>  owned by the calling MPI task?
!>
!>  What do the interpolation weights represent, or how are they used?
!
!-----------------------------------------------------------------------
subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, &
       wght_d_out, lat_d_out, lon_d_out)

   integer, intent(in) :: nxy                     ! array sizes

   real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes
   real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes
   real(r8), intent(out), target, optional :: area_d_out(:)
   real(r8), intent(out), target, optional :: wght_d_out(:) !  weight
   real(r8), intent(out), optional :: lat_d_out(:)  ! column degree latitudes
   real(r8), intent(out), optional :: lon_d_out(:)  ! column degree longitudes

   character(len=*), parameter :: subname = 'dyn_grid::get_horiz_grid_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   if ( nxy /= nCells ) then
      call endrun( subname // ':: incorrect number of cells' )
   end if

   if ( present( clat_d_out ) ) then
      clat_d_out(:) = latCell(:)
   end if

   if ( present( clon_d_out ) ) then
      clon_d_out(:) = lonCell(:)
   end if

   if ( present( area_d_out ) ) then
      area_d_out(:) = areaCell(:) / (6371229.0_r8**2.0_r8)
   end if

   if ( present( wght_d_out ) ) then
      wght_d_out(:) = areaCell(:) / (6371229.0_r8**2.0_r8)
   end if

   if ( present( lat_d_out ) ) then
      lat_d_out(:) = latCell(:) * rad2deg
   end if

   if ( present( lon_d_out ) ) then
      lon_d_out(:) = lonCell(:) * rad2deg
   end if

end subroutine get_horiz_grid_d


!-----------------------------------------------------------------------
!  routine physgrid_copy_attributes_d
!
!> \brief Create list of attributes for physics grid to copy from dynamics grid
!> \details
!>  Create list of attributes for the physics grid that should be copied
!>  from the corresponding grid object on the dynamics decomposition
!>
!>  The purpose of this routine is not entirely clear.
!
!-----------------------------------------------------------------------
subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)

   use cam_grid_support, only: max_hcoordname_len

   character(len=max_hcoordname_len),          intent(out) :: gridname
   character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

   character(len=*), parameter :: subname = 'dyn_grid::physgrid_copy_attributes_d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   gridname = 'mpas_cell'
   allocate(grid_attribute_names(1))
   grid_attribute_names(1) = 'area'

end subroutine physgrid_copy_attributes_d


!-----------------------------------------------------------------------
!  routine get_dyn_grid_parm_real1d
!
!> \brief Not used for unstructured grids
!> \details
!>  This routine is not used for unstructured grids, but still needed as a
!>  dummy interface to satisfy references from mo_synoz.F90 and phys_gmean.F90
!>
!>  If this routine is unused but called, do we need to ensure, e.g., that
!>  rval is nullified before returning?
!
!-----------------------------------------------------------------------
function get_dyn_grid_parm_real1d(name) result(rval)

   character(len=*), intent(in) :: name
   real(r8), pointer :: rval(:)

   character(len=*), parameter :: subname = 'dyn_grid::get_dyn_grid_parm_real1d'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

!   if(name.eq.'w') then
!      call endrun('get_dyn_grid_parm_real1d: w not defined')
!   else if(name.eq.'clat') then
!      call endrun('get_dyn_grid_parm_real1d: clat not supported, use get_horiz_grid_d')
!   else if(name.eq.'latdeg') then
!      call endrun('get_dyn_grid_parm_real1d: latdeg not defined')
!   else
!      nullify(rval)
!   end if

end function get_dyn_grid_parm_real1d


!-----------------------------------------------------------------------
!  routine get_dyn_grid_parm
!
!> \brief Deprecated
!> \details
!>  This function is in the process of being deprecated, but is still needed
!>  as a dummy interface to satisfy external references from some chemistry routines.
!> 
!>  Until this routine is deleted, what values must be returned?
!
!-----------------------------------------------------------------------
integer function get_dyn_grid_parm(name) result(ival)

   character(len=*), intent(in) :: name

   character(len=*), parameter :: subname = 'dyn_grid::get_dyn_grid_parm'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

!   if (name == 'plat') then
!      ival = 1
!   else if (name == 'plon') then
!      ival = nCells
!   else if(name == 'plev') then
!      ival = plev
!   else	
!      ival = -1
!   end if

end function get_dyn_grid_parm


!-----------------------------------------------------------------------
!  routine dyn_grid_get_colndx
!
!> \brief Not sure?
!> \details
!>  The purpose of this routine is unclear. Does it essentially just call
!>  get_gcol_block_d and get_block_owner_d for an array of global columns?
!>
!>  As with get_gcol_block_d, can this routine be called for global columns
!>  that are not owned by the calling task?
!>
!>  When is this routine needed? Previous implementations suggested that it may
!>  never be called for MPAS.
!
!-----------------------------------------------------------------------
subroutine dyn_grid_get_colndx(igcol, ncols, owners, col, lbk )

   integer, intent(in)  :: ncols
   integer, intent(in)  :: igcol(ncols)
   integer, intent(out) :: owners(ncols)
   integer, intent(out) :: col(ncols)
   integer, intent(out) :: lbk(ncols)

   integer  :: i
   integer :: blockid(1), bcid(1), lclblockid(1)

   character(len=*), parameter :: subname = 'dyn_grid::dyn_grid_get_colndx'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

!     do i = 1,ncols
!  
!       call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
!       owners(i) = get_block_owner_d(blockid(1))
!  
!       if ( iam==owners(i) ) then
!          lbk(i) = lclblockid(1)
!          col(i) = bcid(1)
!       else
!          lbk(i) = -1
!          col(i) = -1
!       end if
!  
!    end do

   call endrun('dyn_grid_get_colndx not supported for mpas dycore')

end subroutine dyn_grid_get_colndx


!-----------------------------------------------------------------------
!  routine dyn_grid_get_elem_coords
!
!> \brief Return coordinates of the columns of a block of the dynamics grid
!> \details
!>  Returns the latitude and longitude coordinates, as well as global IDs,
!>  for the columns in a block.
!>
!>  Is the block index a global block index?
!>  Can the block index be that of a block that is not owned by the calling
!>  MPI task?
!>
!>  When is this routine needed? Previous implementations suggested that it may
!>  never be called for MPAS.
!
!-----------------------------------------------------------------------
subroutine dyn_grid_get_elem_coords(ie, rlon, rlat, cdex )

   integer, intent(in) :: ie ! block index

   real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the block
   real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the block
   integer, optional, intent(out) :: cdex(:) ! global column index

   character(len=*), parameter :: subname = 'dyn_grid::dyn_grid_get_elem_coords'


!   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   call endrun('dyn_grid_get_elem_coords not supported for mpas dycore')

end subroutine dyn_grid_get_elem_coords


!=========================================================================================
! Private routines.
!=========================================================================================


subroutine ref_height_read(File)

   ! This code is used both for initial and restart reading.

   type(file_desc_t), intent(inout) :: File

   integer :: ilev, ilev_dimid, ierr
   integer :: pio_errtype
   integer :: start(1), cnt(1)

   type(var_desc_t) :: zw_desc
   real(r8) :: zw_in(plevp)

   character(len=*), parameter :: routine = 'ref_height_read'
   !----------------------------------------------------------------------------

   ! Set PIO to return error codes.
   call pio_seterrorhandling(file, PIO_BCAST_ERROR, pio_errtype)

   ierr = PIO_Inq_DimID(File, 'nVertLevelsP1', ilev_dimid)
   if (ierr == PIO_NOERR) then
      ierr = PIO_Inq_dimlen(File, ilev_dimid, ilev)
      if (ierr == PIO_NOERR) then
         if (plevp /= ilev) then
            write(iulog,*) routine//': ERROR: number file interface levels does not match model. lev (file, model):',ilev, plevp
            call endrun(routine//': ERROR: file nVertLevelsP1 does not match model.')
         end if
      end if
   else
      write(iulog,*) routine//': ERROR: interface level dimension nVertLevelsP1 not found'
      call endrun(routine//': ERROR: interface level dimension nVertLevelsP1 not found')
   end if

   ierr = pio_inq_varid(File, 'zw', zw_desc)
   if (ierr == PIO_NOERR) then
      start=(/1/)
      cnt=(/plevp/)
      ierr = pio_get_var(File, zw_desc, start, cnt, zw_in)
      if (ierr /= PIO_NOERR) then
         write(iulog,*) routine//': ERROR: failed to read data from zw'
         call endrun(routine//': ERROR: failed to read data from zw')
      end if
   else
      write(iulog,*) routine//': ERROR: height data in zw not found'
      call endrun(routine//': ERROR: height data in zw not found')
   end if

   ! reverse level ordering for CAM which uses top to bottom order
   zw = zw_in(plevp:1:-1)

   ! compute layer midpoints
   zw_mid = (zw(1:plev) + zw(2:plevp)) / 2._r8

   ! Put the error handling back the way it was
   call pio_seterrorhandling(file, pio_errtype)

end subroutine ref_height_read

!-----------------------------------------------------------------------
!  routine define_cam_grids
!
!> \brief Define the dynamics and physics grids on the dynamics decomposition
!> \details
!>  Defines the dynamics and physics grids on the dynamics decompostion.
!>  The physics grid on the physics decomposition is defined in phys_grid_init.
!
!-----------------------------------------------------------------------
subroutine define_cam_grids()

   use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
   use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
 
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type

   ! Local variables
   integer :: i, j

   type(horiz_coord_t), pointer     :: lat_coord
   type(horiz_coord_t), pointer     :: lon_coord
   integer(iMap),       allocatable :: coord_map(:)
   integer(iMap),       pointer     :: grid_map(:,:)

   type(mpas_pool_type),   pointer :: meshPool
   integer,  dimension(:), pointer :: indexToCellID
   real(r8), dimension(:), pointer :: latCell
   real(r8), dimension(:), pointer :: lonCell
   real(r8), dimension(:), pointer :: areaCell

   character(len=*), parameter :: subname = 'dyn_grid::define_cam_grids'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(0, 'begin '//subname)
 
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsSolve)   ! MPAS always solves over the full column
   call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

   allocate(coord_map(nCellsSolve))
   coord_map = indexToCellID(1:nCellsSolve)

   lat_coord => horiz_coord_create('lat', 'ncol', nCells, 'latitude',      &
          'degrees_north', 1, nCellsSolve, latCell(1:nCellsSolve)*rad2deg, map=coord_map)
   lon_coord => horiz_coord_create('lon', 'ncol', nCells, 'longitude',     &
          'degrees_east', 1, nCellsSolve, lonCell(1:nCellsSolve)*rad2deg, map=coord_map)
 
   ! Map for cell centers grid
   allocate(grid_map(3, nCellsSolve))
   do i = 1, nCellsSolve
      grid_map(1, i) = i
      grid_map(2, i) = 0
      grid_map(3, i) = coord_map(i)
   end do

   ! cell center grid
   call cam_grid_register('mpas_cell', dyn_decomp, lat_coord, lon_coord,     &
          grid_map, block_indexed=.false., unstruct=.true., src_in=[1,0])
   call cam_grid_attribute_register('mpas_cell', 'area', 'cell areas',  &
          'ncol', areaCell/6371229.0_r8**2.0_r8, coord_map)    ! MGD Should we pass areaCell(1:nCellsSolve) instead?

   ! grid_map cannot be deallocated as the cam_filemap_t object just points
   ! to it.  It can be nullified.
   nullify(grid_map) ! Map belongs to grid now
 
end subroutine define_cam_grids

end module dyn_grid
