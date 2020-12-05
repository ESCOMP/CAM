module physics_column_type

   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   private
   save

type, public :: physics_column_t
   ! A type to hold all grid and task information for a single physics column
   ! Column information
   real(r8)             :: lat_rad = -HUGE(1.0_r8) ! Latitude in radians
   real(r8)             :: lon_rad = -HUGE(1.0_r8) ! Longitude in radians
   real(r8)             :: lat_deg = -HUGE(1.0_r8) ! Latitude in degrees
   real(r8)             :: lon_deg = -HUGE(1.0_r8) ! Longitude in degrees
   real(r8)             :: area = -1.0_r8          ! Column area
   real(r8)             :: weight = -1.0_r8        ! Column integration weight
   ! File decomposition
   integer              :: global_col_num = -1     ! Location on data file
   integer              :: coord_indices(2) = -1   ! Global lon/lat (if used)
   ! Dynamics decomposition
   integer              :: dyn_task = -1           ! Dynamics MPI task
   integer              :: local_dyn_block = -1    ! Block number for this task
   integer              :: global_dyn_block = -1   ! Global dyn block number
   !    If there is more than one block index, they are in the same order
   !    as in the dynamics block structure
   integer, allocatable :: dyn_block_index(:)      ! Index(cies) into block
   ! Physics decomposition
   integer              :: phys_task = -1          ! Physics MPI task
   integer              :: local_phys_chunk = -1   ! Local phys 'block' number
   integer              :: phys_chunk_index = -1   ! Index into physics chunk
end type physics_column_t


!==============================================================================
CONTAINS
!==============================================================================

end module physics_column_type
