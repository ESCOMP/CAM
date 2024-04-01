module dyn_grid

!-------------------------------------------------------------------------------
!
! Define MPAS computational grids on the dynamics decomposition.
!
! Module responsibilities:
!
! . Provide the physics/dynamics coupler (in module phys_grid) with data for the
!   physics grid (cell centers) on the dynamics decomposition.
!
! . Create CAM grid objects that are used by the I/O functionality to read
!   data from an unstructured grid format to the dynamics data structures, and
!   to write from the dynamics data structures to unstructured grid format.  The
!   global column ordering for the unstructured grid is determined by the dycore.
!
! The MPAS grid is decomposed into "blocks" which contain the cells that are solved
! plus a set of halo cells.  The dycore assigns one block per task.
!
!-------------------------------------------------------------------------------


use cam_abortutils,      only: endrun
use cam_logfile,         only: iulog
use cam_mpas_subdriver,  only: domain_ptr, cam_mpas_init_phase3, cam_mpas_get_global_dims, &
                               cam_mpas_read_static, cam_mpas_compute_unit_vectors
use mpas_derived_types,  only: mpas_pool_type
use mpas_pool_routines,  only: mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
use physconst,           only: pi
use physics_column_type, only: physics_column_t
use pio,                 only: file_desc_t, pio_global, pio_get_att
use pmgrid,              only: plev, plevp
use shr_kind_mod,        only: r8 => shr_kind_r8
use spmd_utils,          only: iam, masterproc, mpicom, npes
use string_utils,        only: int2str

implicit none
private
save

integer, parameter :: dyn_decomp    = 101 ! cell center grid (this parameter is public to provide a dycore
                                          ! independent way to identify the physics grid on the dynamics
                                          ! decomposition)
integer, parameter :: cam_cell_decomp = 104 ! same grid decomp as dyn_decomp, but the grid definition
                                            ! uses ncol, lat, lon
integer, parameter :: edge_decomp   = 102 ! edge node grid
integer, parameter :: vertex_decomp = 103 ! vertex node grid
integer, parameter :: ptimelevels = 2

public :: &
   dyn_decomp, &
   ptimelevels, &
   dyn_grid_init, &
   get_dyn_grid_info, &
   get_horiz_grid_dim_d, &
   get_horiz_grid_d, &
   get_dyn_grid_parm, &
   get_dyn_grid_parm_real1d, &
   dyn_grid_get_elem_coords, &
   dyn_grid_get_colndx, &
   physgrid_copy_attributes_d

! vertical reference heights (m) in CAM top to bottom order.
! These arrays are targets of the real_values pointers in the hist_coords
! objects in the cam_history_support module.
real(r8), target :: zw(plevp), zw_mid(plev)

integer ::      &
   maxNCells,   &    ! maximum number of cells for any task (nCellsSolve <= maxNCells)
   maxEdges,    &    ! maximum number of edges per cell
   nVertLevels       ! number of vertical layers (midpoints)

integer, pointer :: &
   nCellsSolve,     & ! number of cells that a task solves
   nEdgesSolve,     & ! number of edges (velocity) that a task solves
   nVerticesSolve,  & ! number of vertices (vorticity) that a task solves
   nVertLevelsSolve

real(r8), parameter :: rad2deg=180.0_r8/pi ! convert radians to degrees

! sphere_radius is a global attribute in the MPAS initial file.  It is needed to
! normalize the cell areas to a unit sphere.
real(r8) :: sphere_radius

! global grid data

integer ::      &
   nCells_g,    &    ! global number of cells/columns
   nEdges_g,    &    ! global number of edges
   nVertices_g       ! global number of vertices

integer, allocatable :: col_indices_in_block(:,:)  ! global column indices in each block
integer, allocatable :: num_col_per_block(:)       ! number of columns in each block
integer, allocatable :: global_blockid(:)          ! block id for each global column
integer, allocatable :: local_col_index(:)         ! local column index (in block) for each global column

!=========================================================================================
contains
!=========================================================================================

subroutine dyn_grid_init()

   ! Initialize grids on the dynamics decomposition and create associated
   ! grid objects for use by I/O utilities.  The current physics/dynamics
   ! coupling code requires constructing global fields for the cell center
   ! grid which is used by the physics parameterizations.

   use ref_pres,            only: ref_pres_init
   use std_atm_profile,     only: std_atm_pres
   use time_manager,        only: get_step_size

   use cam_initfiles,       only: initial_file_get_id

   use cam_history_support, only: add_vert_coord

   use constituents,        only: pcnst

   type(file_desc_t), pointer :: fh_ini

   integer  :: k, ierr
   integer  :: num_pr_lev       ! number of top levels using pure pressure representation
   real(r8) :: pref_edge(plevp) ! reference pressure at layer edges (Pa)
   real(r8) :: pref_mid(plev)   ! reference pressure at layer midpoints (Pa)

   character(len=*), parameter :: subname = 'dyn_grid::dyn_grid_init'
   !----------------------------------------------------------------------------

   ! Get filehandle for initial file
   fh_ini => initial_file_get_id()

   ! MPAS-A always requires at least one scalar (qv).  CAM has the same requirement
   ! and it is enforced by the configure script which sets the cpp macrop PCNST.
   call cam_mpas_init_phase3(fh_ini, pcnst)

   ! Read or compute all time-invariant fields for the MPAS-A dycore
   ! Time-invariant fields are stored in the MPAS mesh pool.  This call
   ! also sets the module data zw and zw_mid.
   call setup_time_invariant(fh_ini)

   ! Read the global sphere_radius attribute.  This is needed to normalize the cell areas.
   ierr = pio_get_att(fh_ini, pio_global, 'sphere_radius', sphere_radius)

   ! Compute reference pressures from reference heights.
   call std_atm_pres(zw, pref_edge)
   pref_mid = (pref_edge(1:plev) + pref_edge(2:plevp)) * 0.5_r8

   num_pr_lev = 0
   call ref_pres_init(pref_edge, pref_mid, num_pr_lev)

   ! Vertical coordinates for output streams
   call add_vert_coord('lev', plev,                       &
         'zeta level at vertical midpoints', 'm', zw_mid)
   call add_vert_coord('ilev', plevp,                     &
         'zeta level at vertical interfaces', 'm', zw)

   if (masterproc) then
      write(iulog,'(a)')' Reference Layer Locations: '
      write(iulog,'(a)')' index      height (m)              pressure (hPa) '
      do k= 1, plev
         write(iulog,9830) k, zw(k), pref_edge(k)/100._r8
         write(iulog,9840)    zw_mid(k), pref_mid(k)/100._r8
      end do
      write(iulog,9830) plevp, zw(plevp), pref_edge(plevp)/100._r8

9830  format(1x, i3, f15.4, 9x, f15.4)
9840  format(1x, 3x, 12x, f15.4, 9x, f15.4)
   end if

   ! Query global grid dimensions from MPAS
   call cam_mpas_get_global_dims(nCells_g, nEdges_g, nVertices_g, maxEdges, nVertLevels, maxNCells)

   ! Define the dynamics grids on the dynamics decompostion.  The cell
   ! centered grid is used by the physics parameterizations.  The physics
   ! decomposition of the cell centered grid is defined in phys_grid_init.
   call define_cam_grids()
   
end subroutine dyn_grid_init

!=========================================================================================

subroutine get_dyn_grid_info(hdim1_d, hdim2_d, num_levels, index_model_top_layer, index_surface_layer, unstructured, dyn_columns)
   !------------------------------------------------------------
   !
   ! get_dyn_grid_info return dynamics grid column information
   !
   ! Return dynamic grid columns information and other dycore information. After
   ! this function is called, dyn_columns will be allocated to size nCellsSolve
   ! and each entry will represent a dynamics column (cell) of the MPAS dynamics.
   !
   !------------------------------------------------------------

   use shr_const_mod, only: SHR_CONST_PI ! TODO: Is this the correct PI constant to use?

   ! Input variables
   integer, intent(out) :: hdim1_d ! Global Longitudes or global grid size (nCells_g)
   integer, intent(out) :: hdim2_d ! Latitudes or 1 for unstructured grids
   integer, intent(out) :: num_levels ! Number of levels
   integer, intent(out) :: index_model_top_layer
   integer, intent(out) :: index_surface_layer
   logical, intent(out) :: unstructured
   type (physics_column_t), allocatable, intent(out):: dyn_columns(:)

   ! Local variables
   type(mpas_pool_type),   pointer :: meshPool
   integer, pointer :: nCellsSolve ! Cells owned by this task, excluding halo cells
   integer, pointer :: nVertLevels ! number of vertical layers (midpoints)
   integer,  dimension(:), pointer :: indexToCellID ! global indices of cell centers
   real(r8), dimension(:), pointer :: latCell   ! cell center latitude (radians)
   real(r8), dimension(:), pointer :: lonCell   ! cell center longitudes (radians)
   real(r8), dimension(:), pointer :: areaCell  ! cell areas in m^2
   integer :: iCell
   integer :: ierr
   integer :: my_proc_id
   character(len=*), parameter :: subname = 'get_dyn_grid_info'

   ! TODO: Is it possible we can guarantee that local_dyn_columns will never be allocated before this function?
   if (allocated(dyn_columns)) then
       call endrun(subname//': dyn_columns must be unallocated')
   end if

   ! Retrieve MPAS grid dimensions and variables from the mesh pool
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

   call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

   allocate(dyn_columns(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate dyn_columns array')

   ! Fill out subroutine arguments
   hdim1_d = nCells_g
   hdim2_d = 1
   num_levels = nVertLevels
   index_model_top_layer = nVertLevels
   index_surface_layer = 1
   unstructured = .True.

   my_proc_id = domain_ptr % dminfo % my_proc_id

   !
   ! Fill out physics_t_column information, one member per cell (column)
   !
   do iCell = 1, nCellsSolve
       ! Column information
       dyn_columns(iCell) % lat_rad = latCell(iCell)
       dyn_columns(iCell) % lon_rad = lonCell(iCell)
       dyn_columns(iCell) % lat_deg = latCell(iCell) * rad2deg
       dyn_columns(iCell) % lon_deg = lonCell(iCell) * rad2deg
       ! Normalize cell areas and cell weights to a unit sphere
       dyn_columns(iCell) % area = areaCell(iCell) / (sphere_radius**2)
       dyn_columns(iCell) % weight = areaCell(iCell) / (sphere_radius**2)

       ! File information
       dyn_columns(iCell) % global_col_num = indexToCellID(iCell)

       ! Dynamics decomposition
       dyn_columns(iCell) % dyn_task = my_proc_id
       dyn_columns(iCell) % local_dyn_block = iCell
       dyn_columns(iCell) % global_dyn_block = indexToCellID(iCell)

       ! dyn_block_index is not used, but it needs to be allocate to a 0 size
       allocate(dyn_columns(iCell) % dyn_block_index(0), stat=ierr)
       if( ierr /= 0 ) call endrun(subname//':failed to allocate dyn_columns%dyn_block_index array')

   end do

end subroutine get_dyn_grid_info

!=========================================================================================

subroutine get_horiz_grid_dim_d(hdim1_d, hdim2_d)

   ! Return declared horizontal dimensions of global grid.
   ! For non-lon/lat grids, declare grid to be one-dimensional,
   ! i.e., (ngcols,1) where ngcols is total number of columns in grid.

   integer, intent(out) :: hdim1_d             ! first horizontal dimension
   integer, intent(out), optional :: hdim2_d   ! second horizontal dimension
   !----------------------------------------------------------------------------

   hdim1_d = nCells_g

   if( present(hdim2_d) ) hdim2_d = 1

end subroutine get_horiz_grid_dim_d

!=========================================================================================

subroutine get_horiz_grid_d(nxy, clat_d_out, clon_d_out, area_d_out, &
       wght_d_out, lat_d_out, lon_d_out)

   ! Return global arrays of latitude and longitude (in radians), column
   ! surface area (in radians squared) and surface integration weights for
   ! columns in physics grid (cell centers)

   integer, intent(in) :: nxy                     ! array sizes

   real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes (radians)
   real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes (radians)
   real(r8), intent(out), target, optional :: area_d_out(:) ! sum to 4*pi (radians^2)
   real(r8), intent(out), target, optional :: wght_d_out(:) ! normalized to sum to 4*pi
   real(r8), intent(out), optional :: lat_d_out(:)  ! column latitudes (degrees)
   real(r8), intent(out), optional :: lon_d_out(:)  ! column longitudes (degrees)

   character(len=*), parameter :: subname = 'dyn_grid::get_horiz_grid_d'
   !----------------------------------------------------------------------------

   call endrun(subname//': NOT SUPPORTED WITH WEAK SCALING FIX')

end subroutine get_horiz_grid_d

!=========================================================================================

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)

   ! Create list of attributes for the physics grid that should be copied
   ! from the corresponding grid object on the dynamics decomposition

   use cam_grid_support, only: max_hcoordname_len

   character(len=max_hcoordname_len),          intent(out) :: gridname
   character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)
   integer :: ierr
   character(len=*), parameter :: subname = 'dyn_grid::physgrid_copy_attributes_d'
   !----------------------------------------------------------------------------


   ! Do not let the physics grid copy the mpas_cell "area" attribute because
   ! it is using a different dimension name.
   gridname = 'mpas_cell'
   allocate(grid_attribute_names(0), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate grid_attribute_names array')


end subroutine physgrid_copy_attributes_d

!=========================================================================================

function get_dyn_grid_parm_real1d(name) result(rval)

   ! This routine is not used for unstructured grids, but still needed as a
   ! dummy interface to satisfy references (for linking executable) from mo_synoz.F90
   ! and phys_gmean.F90.

   character(len=*), intent(in) :: name
   real(r8), pointer :: rval(:)

   character(len=*), parameter :: subname = 'dyn_grid::get_dyn_grid_parm_real1d'
   !----------------------------------------------------------------------------

   if (name .eq. 'w') then
      call endrun(subname//': w not defined')
   else if( name .eq. 'clat') then
      call endrun(subname//': clat not supported, use get_horiz_grid_d')
   else if( name .eq. 'latdeg') then
      call endrun(subname//': latdeg not defined')
   else
      nullify(rval)
   end if

end function get_dyn_grid_parm_real1d

!=========================================================================================

integer function get_dyn_grid_parm(name) result(ival)

   ! This function is in the process of being deprecated, but is still needed
   ! as a dummy interface to satisfy external references from some chemistry routines.

   character(len=*), intent(in) :: name
   !----------------------------------------------------------------------------

   if (name == 'plat') then
      ival = 1
   else if (name == 'plon') then
      ival = nCells_g
   else if(name == 'plev') then
      ival = plev
   else	
      ival = -1
   end if

end function get_dyn_grid_parm

!=========================================================================================

subroutine dyn_grid_get_colndx(igcol, ncols, owners, col, lbk )

   ! For each global column index return the owning task.  If the column is owned
   ! by this task, then also return the local block number and column index in that
   ! block.

   integer, intent(in)  :: ncols
   integer, intent(in)  :: igcol(ncols)
   integer, intent(out) :: owners(ncols)
   integer, intent(out) :: col(ncols)
   integer, intent(out) :: lbk(ncols)

   integer  :: i
   integer :: blockid(1), bcid(1)
   !----------------------------------------------------------------------------

   call endrun('dyn_grid_get_colndx: not implemented for unstructured grids')

end subroutine dyn_grid_get_colndx

!=========================================================================================

subroutine dyn_grid_get_elem_coords(ie, rlon, rlat, cdex )

   ! Returns the latitude and longitude coordinates, as well as global IDs,
   ! for the columns in a block.

   integer, intent(in) :: ie ! block index

   real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the block
   real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the block
   integer, optional, intent(out) :: cdex(:) ! global column index

   character(len=*), parameter :: subname = 'dyn_grid::dyn_grid_get_elem_coords'
   !----------------------------------------------------------------------------

   ! This routine is called for history output when local time averaging is requested
   ! for a field on a dynamics decomposition.  The code in hbuf_accum_addlcltime appears
   ! to also assume that the field is on the physics grid since there is no argument
   ! passed to specify which dynamics grid the coordinates are for.
   
   call endrun(subname//': not implemented for the MPAS grids')

end subroutine dyn_grid_get_elem_coords

!=========================================================================================
! Private routines.
!=========================================================================================

subroutine setup_time_invariant(fh_ini)

   ! Initialize all time-invariant fields needed by the MPAS-Atmosphere dycore,
   ! by reading these fields from the initial file.

   use mpas_rbf_interpolation,     only : mpas_rbf_interp_initialize
   use mpas_vector_reconstruction, only : mpas_init_reconstruct
   use string_utils,               only: int2str

   ! Arguments
   type(file_desc_t), pointer :: fh_ini

   ! Local variables
   type(mpas_pool_type),   pointer :: meshPool
   real(r8), pointer     :: rdzw(:)
   real(r8), allocatable :: dzw(:)

   integer :: k, kk
   integer :: ierr

   character(len=*), parameter :: subname = 'dyn_grid::setup_time_invariant'
   !----------------------------------------------------------------------------

   ! Read time-invariant fields
   call cam_mpas_read_static(fh_ini, endrun)

   ! Compute unit vectors giving the local north and east directions as well as
   ! the unit normal vector for edges
   call cam_mpas_compute_unit_vectors()

   ! Access dimensions that are made public via this module
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)
   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
   call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
   call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsSolve) ! MPAS always solves over the full column

   ! check that number of vertical layers matches MPAS grid data
   if (plev /= nVertLevelsSolve) then
      write(iulog,*) subname//': ERROR: number of levels in IC file does not match plev: file, plev=', &
                     nVertLevelsSolve, plev
      call endrun(subname//': ERROR: number of levels in IC file ('//int2str(nVertLevelsSolve)// &
                           ') does not match plev ('//int2str(nVertLevelsSolve)//').')
   end if

   ! Initialize fields needed for reconstruction of cell-centered winds from edge-normal winds
   ! Note: This same pair of calls happens a second time later in the initialization of
   !       the MPAS-A dycore (in atm_mpas_init_block), but the redundant calls do no harm
   call mpas_rbf_interp_initialize(meshPool)
   call mpas_init_reconstruct(meshPool)

   ! Compute the zeta coordinate at layer interfaces and midpoints.  Store
   ! in arrays using CAM vertical index order (top to bottom of atm) for use
   ! in CAM coordinate objects.
   call mpas_pool_get_array(meshPool, 'rdzw', rdzw)

   allocate(dzw(plev), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate dzw array')

   dzw = 1._r8 / rdzw
   zw(plev+1) = 0._r8
   do k = plev, 1, -1
      kk = plev - k + 1
      zw(k) = zw(k+1) + dzw(kk)
      zw_mid(k) = 0.5_r8 * (zw(k+1) + zw(k))
   end do

   deallocate(dzw)

end subroutine setup_time_invariant

!=========================================================================================

subroutine define_cam_grids()

   ! Define the dynamics grids on the dynamics decompostion.  The 'physics'
   ! grid contains the same nodes as the dynamics cell center grid, but is
   ! on the physics decomposition and is defined in phys_grid_init.
   !
   ! Note that there are two versions of cell center grid defined here.
   ! The 'mpas_cell' grid uses 'nCells' rather than 'ncol' as the dimension
   ! name and 'latCell', 'lonCell' rather than 'lat' and 'lon' as the
   ! coordinate names.  This allows us to read the same initial file that
   ! is used by the standalone MPAS-A model.  The second cell center grid
   ! is called 'cam_cell' and uses the standard CAM names: ncol, lat, and
   ! lon.  This grid allows us to read the PHIS field from the CAM topo
   ! file.  There is just a single version of the grids to read data on the
   ! cell edge and vertex locations.  These are used to read data from the
   ! initial file and to write data from the dynamics decomposition to the
   ! CAM history file.

   use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
   use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
   use shr_const_mod,    only: PI => SHR_CONST_PI
 
   ! Local variables
   integer :: i, j

   type(horiz_coord_t), pointer     :: lat_coord
   type(horiz_coord_t), pointer     :: lon_coord
   integer(iMap),       allocatable :: gidx(:)        ! global indices
   integer(iMap),       pointer     :: grid_map(:,:)

   type(mpas_pool_type),   pointer :: meshPool

   integer,  dimension(:), pointer :: indexToCellID ! global indices of cell centers
   real(r8), dimension(:), pointer :: latCell   ! cell center latitude (radians)
   real(r8), dimension(:), pointer :: lonCell   ! cell center longitude (radians)
   real(r8), dimension(:), pointer :: areaCell  ! cell areas in m^2
   real(r8), dimension(:), pointer :: areaWeight! normalized cell areas weights

   integer,  dimension(:), pointer :: indexToEdgeID ! global indices of edge nodes
   real(r8), dimension(:), pointer :: latEdge   ! edge node latitude (radians)
   real(r8), dimension(:), pointer :: lonEdge   ! edge node longitude (radians)

   integer,  dimension(:), pointer :: indexToVertexID ! global indices of vertex nodes
   real(r8), dimension(:), pointer :: latVertex ! vertex node latitude (radians)
   real(r8), dimension(:), pointer :: lonVertex ! vertex node longitude (radians)
   integer :: ierr
   character(len=*), parameter :: subname = 'dyn_grid::define_cam_grids'
   integer              :: hdim1_d ! Global Longitudes or global grid size (nCells_g)
   integer              :: hdim2_d ! Latitudes or 1 for unstructured grids
   integer              :: num_levels ! Number of levels
   integer              :: index_model_top_layer
   integer              :: index_surface_layer
   logical              :: unstructured
   type (physics_column_t), allocatable  :: dyn_cols(:)
   !----------------------------------------------------------------------------

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', meshPool)

   !-------------------------------------------------------------!
   ! Construct coordinate and grid objects for cell center grid. !
   !-------------------------------------------------------------!

   call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
   call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

   allocate(gidx(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate gidx array at line:'//int2str(__LINE__))

   gidx = indexToCellID(1:nCellsSolve)

   lat_coord => horiz_coord_create('latCell', 'nCells', nCells_g, 'latitude',      &
          'degrees_north', 1, nCellsSolve, latCell(1:nCellsSolve)*rad2deg, map=gidx)
   lon_coord => horiz_coord_create('lonCell', 'nCells', nCells_g, 'longitude',     &
          'degrees_east', 1, nCellsSolve, lonCell(1:nCellsSolve)*rad2deg, map=gidx)
 
   allocate(areaWeight(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate area_weight :'//int2str(__LINE__))
   call get_dyn_grid_info(hdim1_d, hdim2_d, num_levels, index_model_top_layer, index_surface_layer, unstructured, dyn_cols)


   ! Map for cell centers grid
   allocate(grid_map(3, nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate grid_map array at line:'//int2str(__LINE__))

   do i = 1, nCellsSolve
      grid_map(1, i) = i
      grid_map(2, i) = 1
      grid_map(3, i) = gidx(i)
      areaWeight(i) = dyn_cols(i)%weight/(4.0_r8*PI)
   end do

   ! cell center grid for I/O using MPAS names
   call cam_grid_register('mpas_cell', dyn_decomp, lat_coord, lon_coord,     &
          grid_map, block_indexed=.false., unstruct=.true.)
   call cam_grid_attribute_register('mpas_cell', 'area_cell', 'mpas cell areas', &
         'nCells', areaCell, map=gidx)
   call cam_grid_attribute_register('mpas_cell', 'area_weight_mpas', 'mpas area weight', &
         'nCells', areaWeight, map=gidx)

   nullify(areaWeight) ! areaWeight belongs to grid now
   nullify(areaCell) ! areaCell belongs to grid now

   ! create new coordinates and grid using CAM names
   lat_coord => horiz_coord_create('lat', 'ncol', nCells_g, 'latitude',      &
          'degrees_north', 1, nCellsSolve, latCell(1:nCellsSolve)*rad2deg, map=gidx)
   lon_coord => horiz_coord_create('lon', 'ncol', nCells_g, 'longitude',     &
          'degrees_east', 1, nCellsSolve, lonCell(1:nCellsSolve)*rad2deg, map=gidx)
   call cam_grid_register('cam_cell', cam_cell_decomp, lat_coord, lon_coord, &
          grid_map, block_indexed=.false., unstruct=.true.)

   ! gidx can be deallocated.  Values are copied into the coordinate and attribute objects.
   deallocate(gidx)

   deallocate(dyn_cols)

   ! grid_map memory cannot be deallocated.  The cam_filemap_t object just points
   ! to it.  Pointer can be disassociated.
   nullify(grid_map) ! Map belongs to grid now

   ! pointers to coordinate objects can be nullified.  Memory is now pointed to by the
   ! grid object.
   nullify(lat_coord)
   nullify(lon_coord)

   !-----------------------------------------------------------!
   ! Construct coordinate and grid objects for edge node grid. !
   !-----------------------------------------------------------!

   call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
   call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
   call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)

   allocate(gidx(nEdgesSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate gidx array at line:'//int2str(__LINE__))

   gidx = indexToEdgeID(1:nEdgesSolve)

   lat_coord => horiz_coord_create('latEdge', 'nEdges', nEdges_g, 'latitude',      &
          'degrees_north', 1, nEdgesSolve, latEdge(1:nEdgesSolve)*rad2deg, map=gidx)
   lon_coord => horiz_coord_create('lonEdge', 'nEdges', nEdges_g, 'longitude',     &
          'degrees_east', 1, nEdgesSolve, lonEdge(1:nEdgesSolve)*rad2deg, map=gidx)
 
   ! Map for edge node grid
   allocate(grid_map(3, nEdgesSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate grid_map array at line:'//int2str(__LINE__))

   do i = 1, nEdgesSolve
      grid_map(1, i) = i
      grid_map(2, i) = 1
      grid_map(3, i) = gidx(i)
   end do

   ! Edge node grid object
   call cam_grid_register('mpas_edge', edge_decomp, lat_coord, lon_coord,     &
          grid_map, block_indexed=.false., unstruct=.true.)

   deallocate(gidx)
   nullify(grid_map)
   nullify(lat_coord)
   nullify(lon_coord)

   !-------------------------------------------------------------!
   ! Construct coordinate and grid objects for vertex node grid. !
   !-------------------------------------------------------------!

   call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
   call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
   call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

   allocate(gidx(nVerticesSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate gidx array at line:'//int2str(__LINE__))

   gidx = indexToVertexID(1:nVerticesSolve)

   lat_coord => horiz_coord_create('latVertex', 'nVertices', nVertices_g, 'latitude',      &
          'degrees_north', 1, nVerticesSolve, latVertex(1:nVerticesSolve)*rad2deg, map=gidx)
   lon_coord => horiz_coord_create('lonVertex', 'nVertices', nVertices_g, 'longitude',     &
          'degrees_east', 1, nVerticesSolve, lonVertex(1:nVerticesSolve)*rad2deg, map=gidx)
 
   ! Map for vertex node grid
   allocate(grid_map(3, nVerticesSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//':failed to allocate grid_map array at line:'//int2str(__LINE__))

   do i = 1, nVerticesSolve
      grid_map(1, i) = i
      grid_map(2, i) = 1
      grid_map(3, i) = gidx(i)
   end do

   ! Vertex node grid object
   call cam_grid_register('mpas_vertex', vertex_decomp, lat_coord, lon_coord,     &
          grid_map, block_indexed=.false., unstruct=.true.)

   deallocate(gidx)
   nullify(grid_map)
   nullify(lat_coord)
   nullify(lon_coord)
   
end subroutine define_cam_grids

end module dyn_grid
