module phys_grid

!------------------------------------------------------------------------------
!
! The phys_grid module represents the CAM physics decomposition.
!
!  phys_grid_init receives the physics column info (area, weight, centers)
!                 from the dycore.
!                 The routine then creates the physics decomposition which
!                 is the arrangement of columns across the atmosphere model's
!                 MPI tasks as well as the arrangement into groups to
!                 facilitate efficient threading.
!                 The routine then creates a grid object to allow for data
!                 to be read into and written from this decomposition.
! The phys_grid module also provides interfaces for retrieving information
! about the decomposition
!
! Note: This current implementation does not perform load balancing,
!       physics columns ae always on the same task as the corresponding
!       column received from the dycore.
!
!------------------------------------------------------------------------------
   use shr_kind_mod,        only: r8 => shr_kind_r8
   use ppgrid,              only: begchunk, endchunk
   use physics_column_type, only: physics_column_t
   use perf_mod,            only: t_adj_detailf, t_startf, t_stopf

   implicit none
   private
   save

!!XXgoldyXX: v This needs to be removed to complete the weak scaling transition.
   public :: SCATTER_FIELD_TO_CHUNK
!!XXgoldyXX: ^ This needs to be removed to complete the weak scaling transition.

   ! Physics grid management
   public :: phys_grid_init     ! initialize the physics grid
   public :: phys_grid_readnl   ! Read the phys_grid_nl namelist
   public :: phys_grid_initialized
   ! Local task interfaces
   public :: get_nlcols_p       ! Number of local columns
   public :: get_area_p         ! area of a physics column in radians squared
   public :: get_wght_p         ! weight of a physics column in radians squared
   public :: get_rlat_p         ! latitude of a physics column in radians
   public :: get_rlon_p         ! longitude of a physics column in radians
   public :: get_rlat_all_p     ! latitudes of physics cols in chunk (radians)
   public :: get_rlon_all_p     ! longitudes of physics cols in chunk (radians)
   public :: get_lat_p          ! latitude of a physics column in degrees
   public :: get_lon_p          ! longitude of a physics column in degrees
   public :: get_lat_all_p      ! latitudes of physics cols in chunk (degrees)
   public :: get_lon_all_p      ! longitudes of physics cols in chunk (degrees)
   public :: get_area_all_p     ! areas of physics cols in chunk
   public :: get_wght_all_p     ! weights of physics cols in chunk
   public :: get_ncols_p        ! number of columns in a chunk
   public :: get_gcol_p         ! global column index of a physics column
   public :: get_gcol_all_p     ! global col index of all phys cols in a chunk
   public :: get_dyn_col_p      ! dynamics local blk number and blk offset(s)
   public :: get_chunk_info_p   ! chunk index and col # of a physics column
   public :: get_grid_dims      ! return grid dimensions
   ! Physics-dynamics coupling
   public :: phys_decomp_to_dyn ! Transfer physics data to dynamics decomp
   public :: dyn_decomp_to_phys ! Transfer dynamics data to physics decomp

   ! The identifier for the physics grid
   integer, parameter, public          :: phys_decomp = 100

   !! PUBLIC TYPES

   ! Physics chunking (thread blocking) data
   ! Note that chunks cover local data
   type, public :: chunk
      integer, private :: ncols          =  1  ! # of grid columns in this chunk
      integer, private :: chunk_index    = -1  ! Local index of this chunk
      integer, private, allocatable :: phys_cols(:) ! phys column indices
   end type chunk

   !! PRIVATE DATA

   ! dynamics field grid information
   ! hdim1_d and hdim2_d are dimensions of rectangular horizontal grid
   ! data structure, If 1D data structure, then hdim2_d == 1.
   integer                             :: hdim1_d, hdim2_d

   ! Physics decomposition information
   type(physics_column_t), allocatable :: phys_columns(:)

   type(chunk), private, pointer :: chunks(:) => NULL() ! (begchunk:endchunk)

   logical                       :: phys_grid_set = .false.

   logical                       :: calc_memory_increase = .false.

   interface get_dyn_col_p
      module procedure :: get_dyn_col_p_chunk
      module procedure :: get_dyn_col_p_index
   end interface get_dyn_col_p

   ! Private interfaces
   private :: chunk_info_to_index_p

!!XXgoldyXX: v temporary interface to allow old code to compile
   interface get_lat_all_p
      module procedure :: get_lat_all_p_r8 ! The new version
      module procedure :: get_lat_all_p_int ! calls endun
   end interface get_lat_all_p

   interface get_lon_all_p
      module procedure :: get_lon_all_p_r8 ! The new version
      module procedure :: get_lon_all_p_int ! calls endun
   end interface get_lon_all_p
!!XXgoldyXX: ^ temporary interface to allow old code to compile


   integer,          protected, public :: pver = 0
   integer,          protected, public :: pverp = 0
   integer,          protected, public :: num_global_phys_cols = 0
   integer,          protected, public :: columns_on_task = 0
   integer,          protected, public :: index_top_layer = 0
   integer,          protected, public :: index_bottom_layer = 0
   integer,          protected, public :: index_top_interface = 1
   integer,          protected, public :: index_bottom_interface = 0

!==============================================================================
CONTAINS
!==============================================================================

   subroutine phys_grid_readnl(nlfile)
      use cam_abortutils, only: endrun
      use namelist_utils,  only: find_group_name
      use cam_logfile,     only: iulog
      use spmd_utils,      only: mpicom, mstrid=>masterprocid, masterproc
      use spmd_utils,      only: mpi_integer
      use ppgrid,          only: pcols

      character(len=*), intent(in) :: nlfile

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: sub = 'phys_grid_readnl'

      integer :: phys_alltoall = -HUGE(1)
      integer :: phys_loadbalance = -HUGE(1)
      integer :: phys_twin_algorithm = -HUGE(1)
      integer :: phys_chnk_per_thd = -HUGE(1)

      namelist /phys_grid_nl/ phys_alltoall, phys_loadbalance,                &
           phys_twin_algorithm, phys_chnk_per_thd
      !------------------------------------------------------------------------

      ! Read namelist
      if (masterproc) then
         open(newunit=unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'phys_grid_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, phys_grid_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(sub//': FATAL: reading namelist')
            end if
         end if
         close(unitn)
      end if

      call mpi_bcast(phys_alltoall, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpi_bcast(phys_loadbalance, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpi_bcast(phys_twin_algorithm, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpi_bcast(phys_chnk_per_thd, 1, mpi_integer, mstrid, mpicom, ierr)

      if (masterproc) then
         write(iulog,*) 'PHYS_GRID options:'
         write(iulog,*) '  Using PCOLS         =', pcols
         write(iulog,*) '  phys_loadbalance    = (not used)'
         write(iulog,*) '  phys_twin_algorithm = (not used)'
         write(iulog,*) '  phys_alltoall       = (not used)'
         write(iulog,*) '  chunks_per_thread   = (not used)'
      end if

   end subroutine phys_grid_readnl

   !========================================================================

   subroutine phys_grid_init()
      use mpi,              only: MPI_INTEGER, MPI_REAL8, MPI_MIN, MPI_MAX
      use shr_mem_mod,      only: shr_mem_getusage
      use cam_abortutils,   only: endrun
      use cam_logfile,      only: iulog
      use spmd_utils,       only: npes, mpicom, masterprocid, masterproc, iam
      use ppgrid,           only: pcols
      use dyn_grid,         only: get_dyn_grid_info, physgrid_copy_attributes_d
      use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register
      use cam_grid_support, only: iMap, hclen => max_hcoordname_len
      use cam_grid_support, only: horiz_coord_t, horiz_coord_create
      use cam_grid_support, only: cam_grid_attribute_copy, cam_grid_attr_exists
      use shr_const_mod,    only: PI => SHR_CONST_PI

      ! Local variables
      integer                             :: index
      integer                             :: col_index, phys_col
      integer                             :: ichnk, icol, ncol, gcol
      integer                             :: num_chunks
      type(physics_column_t), allocatable :: dyn_columns(:) ! Dyn decomp
      ! Maps and values for physics grid
      real(r8),               pointer     :: lonvals(:)
      real(r8),               pointer     :: latvals(:)
      real(r8)                            :: lonmin, latmin
      integer(iMap),          pointer     :: grid_map(:,:)
      integer(iMap),          allocatable :: coord_map(:)
      type(horiz_coord_t),    pointer     :: lat_coord
      type(horiz_coord_t),    pointer     :: lon_coord
      real(r8),               pointer     :: area_d(:)
      real(r8),               pointer     :: areawt_d(:)
      real(r8)                            :: mem_hw_beg, mem_hw_end
      real(r8)                            :: mem_beg, mem_end
      logical                             :: unstructured
      real(r8)                            :: temp ! For MPI
      integer                             :: ierr ! For MPI
      character(len=hclen),   pointer     :: copy_attributes(:)
      character(len=hclen)                :: copy_gridname
      character(len=*),       parameter   :: subname = 'phys_grid_init: '
      real(r8),               parameter   :: rarea_sphere = 1.0_r8 / (4.0_r8*PI)

      nullify(lonvals)
      nullify(latvals)
      nullify(grid_map)
      nullify(lat_coord)
      nullify(lon_coord)
      nullify(area_d)
      nullify(areawt_d)
      nullify(copy_attributes)

      if (calc_memory_increase) then
         call shr_mem_getusage(mem_hw_beg, mem_beg)
      end if

      call t_adj_detailf(-2)
      call t_startf("phys_grid_init")

      ! Gather info from the dycore
      call get_dyn_grid_info(hdim1_d, hdim2_d, pver, index_top_layer,         &
           index_bottom_layer, unstructured, dyn_columns)
      ! hdim1_d * hdim2_d is the total number of columns
      num_global_phys_cols = hdim1_d * hdim2_d
      pverp = pver + 1
      !!XXgoldyXX: Can we enforce interface numbering separate from dycore?
      !!XXgoldyXX: This will work for both CAM and WRF/MPAS physics
      !!XXgoldyXX: This only has a 50% chance of working on a single level model
      if (index_top_layer < index_bottom_layer) then
         index_top_interface = index_top_layer
         index_bottom_interface = index_bottom_layer + 1
      else
         index_bottom_interface = index_bottom_layer
         index_top_interface = index_top_layer + 1
      end if

      ! Set up the physics decomposition
      columns_on_task = size(dyn_columns)
      if (allocated(phys_columns)) then
         deallocate(phys_columns)
      end if
      allocate(phys_columns(columns_on_task))
      if (columns_on_task > 0) then
         col_index = columns_on_task
         num_chunks = col_index / pcols
         if ((num_chunks * pcols) < col_index) then
            num_chunks = num_chunks + 1
         end if
         begchunk = 1
         endchunk = begchunk + num_chunks - 1
      else
         ! We do not support tasks with no physics columns
         call endrun(subname//'No columns on task, use fewer tasks')
      end if
      allocate(chunks(begchunk:endchunk))
      col_index = 0
      ! Simple chunk assignment
      do index = begchunk, endchunk
         chunks(index)%ncols = MIN(pcols, (columns_on_task - col_index))
         chunks(index)%chunk_index = index
         allocate(chunks(index)%phys_cols(chunks(index)%ncols))
         do phys_col = 1, chunks(index)%ncols
            col_index = col_index + 1
            ! Copy information supplied by the dycore
            phys_columns(col_index) = dyn_columns(col_index)
            ! Fill in physics decomp info
            phys_columns(col_index)%phys_task = iam
            phys_columns(col_index)%local_phys_chunk = index
            phys_columns(col_index)%phys_chunk_index = phys_col
            chunks(index)%phys_cols(phys_col) = col_index
         end do
      end do

      deallocate(dyn_columns)

      ! Add physics-package grid to set of CAM grids
      ! physgrid always uses 'lat' and 'lon' as coordinate names; If dynamics
      !    grid is different, it will use different coordinate names

      ! First, create a map for the physics grid
      ! It's structure will depend on whether or not the physics grid is
      ! unstructured
      if (unstructured) then
         allocate(grid_map(3, pcols * (endchunk - begchunk + 1)))
      else
         allocate(grid_map(4, pcols * (endchunk - begchunk + 1)))
      end if
      grid_map = 0_iMap
      allocate(latvals(size(grid_map, 2)))
      allocate(lonvals(size(grid_map, 2)))

      lonmin = 1000.0_r8 ! Out of longitude range
      latmin = 1000.0_r8 ! Out of latitude range
      index = 0
      do ichnk = begchunk, endchunk
         ncol = chunks(ichnk)%ncols ! Too soon to call get_ncols_p
         do icol = 1, pcols
            index = index + 1
            if (icol <= ncol) then
               col_index = chunks(ichnk)%phys_cols(icol)
               latvals(index) = phys_columns(col_index)%lat_deg
               if (latvals(index) < latmin) then
                  latmin = latvals(index)
               end if
               lonvals(index) = phys_columns(col_index)%lon_deg
               if (lonvals(index) < lonmin) then
                  lonmin = lonvals(index)
               end if
            else
               col_index = -1
               latvals(index) = 1000.0_r8
               lonvals(index) = 1000.0_r8
            end if
            grid_map(1, index) = int(icol, iMap)
            grid_map(2, index) = int(ichnk, iMap)
            if (icol <= ncol) then
               if (unstructured) then
                  gcol = phys_columns(col_index)%global_col_num
                  if (gcol > 0) then
                     grid_map(3, index) = int(gcol, iMap)
                  end if ! else entry remains 0
               else
                  ! lon
                  gcol = phys_columns(col_index)%coord_indices(1)
                  if (gcol > 0) then
                     grid_map(3, index) = int(gcol, iMap)
                  end if ! else entry remains 0
                  ! lat
                  gcol = phys_columns(col_index)%coord_indices(2)
                  if (gcol > 0) then
                     grid_map(4, index) = gcol
                  end if ! else entry remains 0
               end if
            end if ! Else entry remains 0
         end do
      end do

      ! Note that if the dycore is using the same points as the physics grid,
      !      it will have already set up 'lat' and 'lon' axes for
      !      the physics grid
      !      However, these will be in the dynamics decomposition

      if (unstructured) then
         lon_coord => horiz_coord_create('lon', 'ncol', num_global_phys_cols, &
              'longitude', 'degrees_east', 1, size(lonvals), lonvals,         &
              map=grid_map(3,:))
         lat_coord => horiz_coord_create('lat', 'ncol', num_global_phys_cols, &
              'latitude', 'degrees_north', 1, size(latvals), latvals,         &
              map=grid_map(3,:))
      else
         allocate(coord_map(size(grid_map, 2)))
         ! We need a global minimum longitude and latitude
         if (npes > 1) then
            temp = lonmin
            call MPI_allreduce(temp, lonmin, 1, MPI_INTEGER, MPI_MIN,         &
                 mpicom, ierr)
            temp = latmin
            call MPI_allreduce(temp, latmin, 1, MPI_INTEGER, MPI_MIN,         &
                 mpicom, ierr)
            ! Create lon coord map which only writes from one of each unique lon
            where(latvals == latmin)
               coord_map(:) = grid_map(3, :)
            elsewhere
               coord_map(:) = 0_iMap
            end where
            lon_coord => horiz_coord_create('lon', 'lon', hdim1_d,            &
                 'longitude', 'degrees_east', 1, size(lonvals), lonvals,      &
                 map=coord_map)

            ! Create lat coord map which only writes from one of each unique lat
            where(lonvals == lonmin)
               coord_map(:) = grid_map(4, :)
            elsewhere
               coord_map(:) = 0_iMap
            end where
            lat_coord => horiz_coord_create('lat', 'lat', hdim2_d,            &
                 'latitude', 'degrees_north', 1, size(latvals), latvals,      &
                 map=coord_map)
            deallocate(coord_map)
         end if
      end if
      call cam_grid_register('physgrid', phys_decomp, lat_coord, lon_coord,   &
           grid_map, unstruct=unstructured, block_indexed=.true.)
      ! Copy required attributes from the dynamics array
      nullify(copy_attributes)
      call physgrid_copy_attributes_d(copy_gridname, copy_attributes)
      do index = 1, size(copy_attributes)
         call cam_grid_attribute_copy(copy_gridname, 'physgrid',              &
              copy_attributes(index))
      end do

      if (.not. cam_grid_attr_exists('physgrid', 'area')) then
         ! Physgrid always needs an area attribute.
         if (unstructured) then
            ! If we did not inherit one from the dycore (i.e., physics and
            ! dynamics are on different grids), create that attribute here
            ! (Note, a separate physics grid is only supported for
            !  unstructured grids).
            allocate(area_d(size(grid_map, 2)))
            do col_index = 1, columns_on_task
               area_d(col_index) = phys_columns(col_index)%area
            end do
            call cam_grid_attribute_register('physgrid', 'area',              &
                 'physics column areas', 'ncol', area_d, map=grid_map(3,:))
            nullify(area_d) ! Belongs to attribute now

            allocate(areawt_d(size(grid_map, 2)))
            do col_index = 1, columns_on_task
               areawt_d(col_index) = phys_columns(col_index)%weight*rarea_sphere
            end do
            call cam_grid_attribute_register('physgrid', 'areawt',              &
                 'physics column area weight', 'ncol', areawt_d, map=grid_map(3,:))
            nullify(areawt_d) ! Belongs to attribute now
         else
            call endrun(subname//"No 'area' attribute from dycore")
         end if
      end if
      ! Cleanup pointers (they belong to the grid now)
      nullify(grid_map)
      deallocate(latvals)
      nullify(latvals)
      deallocate(lonvals)
      nullify(lonvals)
      ! Cleanup, we are responsible for copy attributes
      if (associated(copy_attributes)) then
         deallocate(copy_attributes)
         nullify(copy_attributes)
      end if

      ! Set flag indicating physics grid is now set
      phys_grid_set = .true.

      call t_stopf("phys_grid_init")
      call t_adj_detailf(+2)

      if (calc_memory_increase) then
         call shr_mem_getusage(mem_hw_end, mem_end)
         temp = mem_end - mem_beg
         call MPI_reduce(temp, mem_end, 1, MPI_REAL8, MPI_MAX, masterprocid,  &
              mpicom, ierr)
         if (masterproc) then
            write(iulog, *) 'phys_grid_init: Increase in memory usage = ',    &
                 mem_end, ' (MB)'
         end if
         temp = mem_hw_end - mem_hw_beg
         call MPI_reduce(temp, mem_hw_end, 1, MPI_REAL8, MPI_MAX,             &
              masterprocid, mpicom, ierr)
         if (masterproc) then
            write(iulog, *) subname, 'Increase in memory highwater = ',       &
                 mem_end, ' (MB)'
         end if
      end if

   end subroutine phys_grid_init

   !========================================================================

   integer function chunk_info_to_index_p(lcid, col, subname_in)
      use cam_logfile,    only: iulog
      use cam_abortutils, only: endrun
      ! Return the physics column index indicated by
      ! <lcid> (chunk) and <col> (column).

      ! Dummy arguments
      integer,                    intent(in) :: lcid ! local chunk id
      integer,                    intent(in) :: col  ! Column index
      character(len=*), optional, intent(in) :: subname_in
      ! Local variables
      character(len=128)          :: errmsg
      character(len=*), parameter :: subname = 'chunk_info_to_index_p: '

      if (.not. phys_grid_initialized()) then
         if (present(subname_in)) then
            call endrun(trim(subname_in)//'physics grid not initialized')
         else
            call endrun(subname//'physics grid not initialized')
         end if
      else if ((lcid < begchunk) .or. (lcid > endchunk)) then
         if (present(subname_in)) then
            write(errmsg, '(a,3(a,i0))') trim(subname_in), 'lcid (', lcid,    &
                 ') out of range (', begchunk, ' to ', endchunk
         else
            write(errmsg, '(a,3(a,i0))') subname, 'lcid (', lcid,             &
                 ') out of range (', begchunk, ' to ', endchunk
         end if
         write(iulog, *) trim(errmsg)
         call endrun(trim(errmsg))
      else if ((col < 1) .or. (col > get_ncols_p(lcid))) then
         if (present(subname_in)) then
            write(errmsg, '(a,2(a,i0))') trim(subname_in), 'col (', col,      &
                 ') out of range (1 to ', get_ncols_p(lcid)
         else
            write(errmsg, '(a,2(a,i0))') subname, 'col (', col,               &
                 ') out of range (1 to ', get_ncols_p(lcid)
         end if
         write(iulog, *) trim(errmsg)
         call endrun(trim(errmsg))
      end if
      chunk_info_to_index_p = chunks(lcid)%phys_cols(col)
   end function chunk_info_to_index_p

   !========================================================================

   logical function phys_grid_initialized()
      ! Return .true. if the physics grid is initialized, otherwise .false.
      phys_grid_initialized = phys_grid_set
   end function phys_grid_initialized

   !========================================================================

   integer function get_nlcols_p()
      get_nlcols_p = columns_on_task
   end function get_nlcols_p

   !========================================================================

   real(r8) function get_rlat_p(lcid, col)
      !-----------------------------------------------------------------------
      !
      ! get_rlat_p: latitude of a physics column in radians
      !
      !-----------------------------------------------------------------------

      ! Dummy argument
      integer, intent(in) :: lcid
      integer, intent(in) :: col
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_rlat_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_rlat_p = phys_columns(index)%lat_rad

   end function get_rlat_p

   !========================================================================

   real(r8) function get_rlon_p(lcid, col)
      !-----------------------------------------------------------------------
      !
      ! get_rlon_p: longitude of a physics column in radians
      !
      !-----------------------------------------------------------------------

      ! Dummy argument
      integer, intent(in) :: lcid
      integer, intent(in) :: col
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_rlon_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_rlon_p = phys_columns(index)%lon_rad

   end function get_rlon_p

   !========================================================================

   subroutine get_rlat_all_p(lcid, rlatdim, rlats)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_rlat_all_p: Return all latitudes (in radians) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid           ! local chunk id
      integer,  intent(in)  :: rlatdim        ! declared size of output array
      real(r8), intent(out) :: rlats(rlatdim) ! array of latitudes

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_rlat_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), rlatdim)
         phys_ind = chunks(lcid)%phys_cols(index)
         rlats(index) = phys_columns(phys_ind)%lat_rad
      end do

   end subroutine get_rlat_all_p

   !========================================================================

   subroutine get_rlon_all_p(lcid, rlondim, rlons)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_rlon_all_p:: Return all longitudes (in radians) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid           ! local chunk id
      integer,  intent(in)  :: rlondim        ! declared size of output array
      real(r8), intent(out) :: rlons(rlondim) ! array of longitudes

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_rlon_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), rlondim)
         phys_ind = chunks(lcid)%phys_cols(index)
         rlons(index) = phys_columns(phys_ind)%lon_rad
      end do

   end subroutine get_rlon_all_p

   !========================================================================

   real(r8) function get_lat_p(lcid, col)
      !-----------------------------------------------------------------------
      !
      ! get_lat_p: latitude of a physics column in degrees
      !
      !-----------------------------------------------------------------------

      ! Dummy argument
      integer, intent(in) :: lcid
      integer, intent(in) :: col
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_lat_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_lat_p = phys_columns(index)%lat_deg

   end function get_lat_p

   !========================================================================

   real(r8) function get_lon_p(lcid, col)
      !-----------------------------------------------------------------------
      !
      ! get_lon_p: longitude of a physics column in degrees
      !
      !-----------------------------------------------------------------------

      ! Dummy argument
      integer, intent(in) :: lcid
      integer, intent(in) :: col
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_lon_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_lon_p = phys_columns(index)%lon_deg

   end function get_lon_p

   !========================================================================

   subroutine get_lat_all_p_r8(lcid, latdim, lats)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_lat_all_p: Return all latitudes (in degrees) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid         ! local chunk id
      integer,  intent(in)  :: latdim       ! declared size of output array
      real(r8), intent(out) :: lats(latdim) ! array of latitudes

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_lat_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), latdim)
         phys_ind = chunks(lcid)%phys_cols(index)
         lats(index) = phys_columns(phys_ind)%lat_deg
      end do

   end subroutine get_lat_all_p_r8

   !========================================================================

   subroutine get_lon_all_p_r8(lcid, londim, lons)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_lon_all_p:: Return all longitudes (in degrees) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid          ! local chunk id
      integer,  intent(in)  :: londim        ! declared size of output array
      real(r8), intent(out) :: lons(londim)  ! array of longitudes

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_lon_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), londim)
         phys_ind = chunks(lcid)%phys_cols(index)
         lons(index) = phys_columns(phys_ind)%lon_deg
      end do

   end subroutine get_lon_all_p_r8

   !========================================================================

   subroutine get_area_all_p(lcid, areadim, areas)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_area_all_p: Return all areas for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid           ! local chunk id
      integer,  intent(in)  :: areadim        ! declared size of output array
      real(r8), intent(out) :: areas(areadim) ! array of areas

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_area_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), areadim)
         phys_ind = chunks(lcid)%phys_cols(index)
         areas(index) = phys_columns(phys_ind)%area
      end do

   end subroutine get_area_all_p

   !========================================================================

   subroutine get_wght_all_p(lcid, wghtdim, wghts)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_wght_all_p: Return all weights for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid           ! local chunk id
      integer,  intent(in)  :: wghtdim        ! declared size of output array
      real(r8), intent(out) :: wghts(wghtdim) ! array of weights

      ! Local variables
      integer                     :: index ! loop index
      integer                     :: phys_ind
      character(len=*), parameter :: subname = 'get_wght_all_p: '

      !-----------------------------------------------------------------------
      if ((lcid < begchunk) .or. (lcid > endchunk)) then
         call endrun(subname//'chunk index out of range')
      end if
      do index = 1, MIN(get_ncols_p(lcid, subname_in=subname), wghtdim)
         phys_ind = chunks(lcid)%phys_cols(index)
         wghts(index) = phys_columns(phys_ind)%weight
      end do

   end subroutine get_wght_all_p

   !========================================================================

   integer function get_ncols_p(lcid, subname_in)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_ncols_p: Return number of columns in chunk given the local chunk id.
      !
      !-----------------------------------------------------------------------
      ! Dummy arguments
      integer, intent(in)  :: lcid      ! local chunk id
      character(len=*), optional, intent(in) :: subname_in

      if (.not. phys_grid_initialized()) then
         if (present(subname_in)) then
            call endrun(trim(subname_in)//'physics grid not initialized')
         else
            call endrun('get_ncols_p: physics grid not initialized')
         end if
      else
         get_ncols_p = chunks(lcid)%ncols
      end if

   end function get_ncols_p

   !========================================================================

   real(r8) function get_area_p(lcid, col)
      ! area of a physics column in radians squared

      ! Dummy arguments
      integer, intent(in) :: lcid ! Chunk number
      integer, intent(in) :: col  ! <lcid> column
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_area_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_area_p = phys_columns(index)%area

   end function get_area_p

   !========================================================================

   real(r8) function get_wght_p(lcid, col)
      ! weight of a physics column in radians squared

      ! Dummy arguments
      integer, intent(in) :: lcid ! Chunk number
      integer, intent(in) :: col  ! <lcid> column
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_wght_p'

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_wght_p = phys_columns(index)%weight

   end function get_wght_p

   !========================================================================

   integer function get_gcol_p(lcid, col)
      ! global column index of a physics column

      ! Dummy arguments
      integer, intent(in)  :: lcid          ! local chunk id
      integer, intent(in)  :: col           ! column index
      ! Local variables
      integer                     :: index
      character(len=*), parameter :: subname = 'get_gcol_p: '

      index = chunk_info_to_index_p(lcid, col, subname_in=subname)
      get_gcol_p = phys_columns(index)%global_col_num

   end function get_gcol_p

   !========================================================================

   subroutine get_dyn_col_p_chunk(lcid, col, blk_num, blk_ind, caller)
      use cam_abortutils, only: endrun
      ! Return the dynamics local block number and block offset(s) for
      ! the physics column indicated by <lcid> (chunk) and <col> (column).

      ! Dummy arguments
      integer, intent(in)  :: lcid          ! local chunk id
      integer, intent(in)  :: col           ! Column index
      integer, intent(out) :: blk_num       ! Local dynamics block index
      integer, intent(out) :: blk_ind(:)    ! Local dynamics block offset(s)
      character(len=*), optional, intent(in) :: caller ! Calling routine
      ! Local variables
      integer                     :: index
      integer                     :: off_size
      character(len=*), parameter :: subname = 'get_dyn_col_p_chunk: '

      index = chunk_info_to_index_p(lcid, col)
      off_size = SIZE(phys_columns(index)%dyn_block_index, 1)
      if (SIZE(blk_ind, 1) < off_size) then
         if (present(caller)) then
            call endrun(trim(caller)//': blk_ind too small')
         else
            call endrun(subname//'blk_ind too small')
         end if
      end if
      blk_num = phys_columns(index)%local_dyn_block
      blk_ind(1:off_size) = phys_columns(index)%dyn_block_index(1:off_size)
      if (SIZE(blk_ind, 1) > off_size) then
         blk_ind(off_size+1:) = -1
      end if

   end subroutine get_dyn_col_p_chunk

   !========================================================================

   subroutine get_dyn_col_p_index(index, blk_num, blk_ind)
      use cam_logfile,    only: iulog
      use cam_abortutils, only: endrun
      ! Return the dynamics local block number and block offset(s) for
      ! the physics column indicated by <index>.

      ! Dummy arguments
      integer, intent(in)  :: index         ! index of local physics column
      integer, intent(out) :: blk_num       ! Local dynamics block index
      integer, intent(out) :: blk_ind(:)    ! Local dynamics block offset(s)
      ! Local variables
      integer                     :: off_size
      character(len=128)          :: errmsg
      character(len=*), parameter :: subname = 'get_dyn_col_p_index: '

      if (.not. phys_grid_initialized()) then
         call endrun(subname//'physics grid not initialized')
      else if ((index < 1) .or. (index > columns_on_task)) then
         write(errmsg, '(a,2(a,i0))') subname, 'index (', index,              &
              ') out of range (1 to ', columns_on_task
         write(iulog, *) trim(errmsg)
         call endrun(trim(errmsg))
      else
         off_size = SIZE(phys_columns(index)%dyn_block_index, 1)
         if (SIZE(blk_ind, 1) < off_size) then
            call endrun(subname//'blk_ind too small')
         end if
         blk_num = phys_columns(index)%local_dyn_block
         blk_ind(1:off_size) = phys_columns(index)%dyn_block_index(1:off_size)
         if (SIZE(blk_ind, 1) > off_size) then
            blk_ind(off_size+1:) = -1
         end if
      end if

   end subroutine get_dyn_col_p_index

   !========================================================================

   subroutine get_gcol_all_p(lcid, gdim, gcols)
      use cam_logfile,    only: iulog
      use cam_abortutils, only: endrun
      use spmd_utils,     only: masterproc
      ! collect global column indices of all physics columns in a chunk

      ! Dummy arguments
      integer, intent(in)  :: lcid          ! local chunk id
      integer, intent(in)  :: gdim          ! gcols dimension
      integer, intent(out) :: gcols(:)      ! global column indices
      ! Local variables
      integer                     :: ncol, col_ind
      character(len=128)          :: errmsg
      character(len=*), parameter :: subname = 'get_gcol_all_p: '

      if (.not. phys_grid_initialized()) then
         call endrun(subname//'physics grid not initialized')
      else if ((lcid < begchunk) .or. (lcid > endchunk)) then
         write(errmsg, '(a,3(a,i0))') subname, 'lcid (', lcid,                &
              ') out of range (', begchunk, ' to ', endchunk
         write(iulog, *) trim(errmsg)
         call endrun(trim(errmsg))
      else
         ncol = chunks(lcid)%ncols
         if (gdim < ncol) then
            if (masterproc) then
               write(iulog, '(2a,2(i0,a))') subname, 'WARNING: gdim (', gdim, &
                    ') < ncol (', ncol,'), not all indices will be filled.'
            end if
            gcols(gdim+1:ncol) = -1
         end if
         do col_ind = 1, MIN(ncol, gdim)
            gcols(col_ind) = get_gcol_p(lcid, col_ind)
         end do
      end if

   end subroutine get_gcol_all_p

   !========================================================================

   subroutine get_chunk_info_p(index, lchnk, icol)
      use cam_logfile,    only: iulog
      use cam_abortutils, only: endrun
      ! local chunk index and column number of a physics column

      ! Dummy arguments
      integer, intent(in)  :: index
      integer, intent(out) :: lchnk
      integer, intent(out) :: icol
      ! Local variables
      character(len=128)          :: errmsg
      character(len=*), parameter :: subname = 'get_chunk_info_p: '

      if (.not. phys_grid_initialized()) then
         call endrun(subname//': physics grid not initialized')
      else if ((index < 1) .or. (index > columns_on_task)) then
         write(errmsg, '(a,2(a,i0))') subname, 'index (', index,              &
              ') out of range (1 to ', columns_on_task
         write(iulog, *) errmsg
         call endrun(errmsg)
      else
         lchnk = phys_columns(index)%local_phys_chunk
         icol = phys_columns(index)%phys_chunk_index
      end if

   end subroutine get_chunk_info_p

   !========================================================================

   subroutine get_grid_dims(hdim1_d_out, hdim2_d_out)
      use cam_abortutils, only: endrun
      ! retrieve dynamics field grid information
      ! hdim1_d and hdim2_d are dimensions of rectangular horizontal grid
      ! data structure, If 1D data structure, then hdim2_d == 1.
      integer, intent(out) :: hdim1_d_out
      integer, intent(out) :: hdim2_d_out

      if (.not. phys_grid_initialized()) then
         call endrun('get_grid_dims: physics grid not initialized')
      end if
      hdim1_d_out = hdim1_d
      hdim2_d_out = hdim2_d

   end subroutine get_grid_dims

   !========================================================================

   ! Note: This routine is a stub for future load-balancing
   subroutine phys_decomp_to_dyn()
      !-----------------------------------------------------------------------
      !
      ! phys_decomp_to_dyn: Transfer physics data to dynamics decomp
      !
      !-----------------------------------------------------------------------
   end subroutine phys_decomp_to_dyn

   !========================================================================

   ! Note: This routine is a stub for future load-balancing
   subroutine dyn_decomp_to_phys()
      !-----------------------------------------------------------------------
      !
      ! dyn_decomp_to_phys: Transfer dynamics data to physics decomp
      !
      !-----------------------------------------------------------------------

   end subroutine dyn_decomp_to_phys

   !========================================================================

   subroutine dump_grid_map(grid_map)
      use spmd_utils,       only: iam, npes, mpicom
      use cam_grid_support, only: iMap

      integer(iMap), pointer :: grid_map(:,:)

      integer                :: num_cols
      integer                :: penum, icol
      logical                :: unstruct
      integer                :: file
      integer                :: ierr

      unstruct = SIZE(grid_map, 1) == 3
      num_cols = SIZE(grid_map, 2)
      if (iam == 0) then
         open(newunit=file, file='physgrid_map.csv', status='replace')
         if (unstruct) then
            write(file, *) '"iam","col","block","map pos"'
         else
            write(file, *) '"iam","col","block","lon","lat"'
         end if
         close(unit=file)
      end if
      do penum = 0, npes - 1
         if (iam == penum) then
            open(newunit=file, file='physgrid_map.csv', status='old',         &
                 action='readwrite', position='append')
            do icol = 1, num_cols
               if (unstruct) then
                  write(file, '(3(i0,","),i0)') iam, int(grid_map(1,icol)),   &
                       int(grid_map(2,icol)), int(grid_map(3,icol))
               else
                  write(file, '(4(i0,","),i0)') iam, int(grid_map(1,icol)),   &
                       int(grid_map(2,icol)), int(grid_map(3,icol)),          &
                       int(grid_map(4,icol))
               end if
            end do
            close(unit=file)
         end if
         call MPI_barrier(mpicom, ierr)
      end do
   end subroutine dump_grid_map

!=============================================================================
!==
!!!!!! DUMMY INTERFACEs TO TEST WEAK SCALING INFRASTRUCTURE, SHOULD GO AWAY
!==
!=============================================================================

   subroutine scatter_field_to_chunk(fdim,mdim,ldim, &
                                     hdim1d,globalfield,localchunks)
      use cam_abortutils, only: endrun
      use ppgrid, only: pcols
      !-----------------------------------------------------------------------
      !
      ! Purpose: DUMMY FOR WEAK SCALING TESTS
      !
      !------------------------------Arguments--------------------------------
      integer, intent(in) :: fdim      ! declared length of first dimension
      integer, intent(in) :: mdim      ! declared length of middle dimension
      integer, intent(in) :: ldim      ! declared length of last dimension
      integer, intent(in) :: hdim1d    ! declared first horizontal index
      real(r8), intent(in) :: globalfield(fdim,hdim1d,mdim,hdim2_d,ldim)
      real(r8), intent(out):: localchunks(fdim,pcols,mdim, &
           begchunk:endchunk,ldim)

      call endrun('scatter_field_to_chunk: NOT SUPPORTED WITH WEAK SCALING')
   end subroutine scatter_field_to_chunk

   !========================================================================

   subroutine get_lat_all_p_int(lcid, latdim, lats)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_lat_all_p: Return all latitudes (in degrees) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid         ! local chunk id
      integer,  intent(in)  :: latdim       ! declared size of output array
      integer, intent(out) :: lats(latdim) ! array of latitudes

      call endrun('get_lat_all_p: deprecated interface')

   end subroutine get_lat_all_p_int

   !========================================================================

   subroutine get_lon_all_p_int(lcid, londim, lons)
      use cam_abortutils, only: endrun
      !-----------------------------------------------------------------------
      !
      ! get_lon_all_p:: Return all longitudes (in degrees) for chunk, <lcid>
      !
      !-----------------------------------------------------------------------
      ! Dummy Arguments
      integer,  intent(in)  :: lcid          ! local chunk id
      integer,  intent(in)  :: londim        ! declared size of output array
      integer, intent(out) :: lons(londim)  ! array of longitudes

      call endrun('get_lon_all_p: deprecated interface')

   end subroutine get_lon_all_p_int

   !========================================================================

end module phys_grid
