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
      real(r8)             :: weight = -1.0_r8        ! Col integration weight
      ! File decomposition
      integer              :: global_col_num = -1     ! Location on data file
      integer              :: coord_indices(2) = -1   ! Global lon/lat (if used)
      ! Dynamics decomposition
      integer              :: dyn_task = -1           ! Dynamics MPI task
      integer              :: local_dyn_block = -1    ! Block num for this task
      integer              :: global_dyn_block = -1   ! Global dyn block number
      !    If there is more than one block index, they are in the same order
      !    as in the dynamics block structure
      integer, allocatable :: dyn_block_index(:)      ! Index(cies) into block
      ! Physics decomposition
      integer              :: phys_task = -1          ! Physics MPI task
      integer              :: local_phys_chunk = -1   ! Local phys 'block' num
      integer              :: phys_chunk_index = -1   ! Index into physics chunk
   contains
      procedure :: copyColumn
      generic :: assignment(=) => copyColumn
   end type physics_column_t

!==============================================================================
CONTAINS
!==============================================================================

   subroutine copyColumn(outCol, inCol)
      ! Dummy arguments
      class(physics_column_t), intent(inout) :: outCol
      type(physics_column_t),  intent(in)    :: inCol
      ! Local variables
      integer                                :: nind ! # dynamics indices

      outCol%lat_rad            = inCol%lat_rad
      outCol%lon_rad            = inCol%lon_rad
      outCol%lat_deg            = inCol%lat_deg
      outCol%lon_deg            = inCol%lon_deg
      outCol%area               = inCol%area
      outCol%weight             = inCol%weight
      outCol%global_col_num     = inCol%global_col_num
      outCol%coord_indices(:)   = inCol%coord_indices(2)
      outCol%dyn_task           = inCol%dyn_task
      outCol%local_dyn_block    = inCol%local_dyn_block
      outCol%global_dyn_block   = inCol%global_dyn_block
      nind                      = SIZE(inCol%dyn_block_index)
      allocate(outCol%dyn_block_index(nind))
      outCol%dyn_block_index(:) = inCol%dyn_block_index(:)
      outCol%phys_task          = inCol%phys_task
      outCol%local_phys_chunk   = inCol%local_phys_chunk
      outCol%phys_chunk_index   = inCol%phys_chunk_index
   end subroutine copyColumn

end module physics_column_type
