module edgetype_mod 

  use shr_kind_mod,           only: r8=>shr_kind_r8, i8=>shr_kind_i8
  use coordinate_systems_mod, only : cartesian3D_t
  use gbarriertype_mod,       only : gbarrier_t

  implicit none 
  private 
  save 

  integer, public :: initedgebuffer_callid = 0

  type, public :: rotation_t
    integer                 :: nbr                ! nbr direction: north south east west
    integer                 :: reverse            ! 0 = do not reverse order
    ! 1 = reverse order
    real (kind=r8), pointer :: R(:,:,:) => null() !  rotation matrix
  end type rotation_t

  type, public :: EdgeDescriptor_t
     integer                      :: use_rotation
     integer                      :: padding
     integer,             pointer :: putmapP(:) => null()
     integer,             pointer :: getmapP(:) => null()
     integer,             pointer :: putmapP_ghost(:) => null()
     integer,             pointer :: getmapP_ghost(:) => null()
     integer,             pointer :: putmapS(:) => null()
     integer,             pointer :: getmapS(:) => null()
     integer,             pointer :: globalID(:) => null()
     integer,             pointer :: loc2buf(:) => null()
     type(cartesian3D_t), pointer :: neigh_corners(:,:) => null()
     integer                      :: actual_neigh_edges
     logical,             pointer :: reverse(:) => null()
     type (rotation_t),   pointer :: rot(:) => null() !  Identifies list of edges
     !  that must be rotated, and how
  end type EdgeDescriptor_t

  type, public :: EdgeBuffer_t
     real (kind=r8), allocatable :: buf(:)
     real (kind=r8), allocatable :: receive(:)
     integer,        allocatable :: putmap(:,:)
     integer,        allocatable :: getmap(:,:)
     logical,        allocatable :: reverse(:,:)
     integer,        allocatable :: moveLength(:)
     integer,        allocatable :: movePtr(:)
     integer,        allocatable :: rcountsFull(:)
     integer,        allocatable :: scountsFull(:)
     integer,        allocatable :: sdisplsFull(:)
     integer,        allocatable :: rdisplsFull(:)
     integer,        allocatable :: rcountsInter(:)
     integer,        allocatable :: scountsInter(:)
     integer,        allocatable :: sdisplsInter(:)
     integer,        allocatable :: rdisplsInter(:)
     integer,        allocatable :: rcountsIntra(:)
     integer,        allocatable :: scountsIntra(:)
     integer,        allocatable :: sdisplsIntra(:)
     integer,        allocatable :: rdisplsIntra(:)
     integer,        allocatable :: getDisplsFull(:)
     integer,        allocatable :: putDisplsFull(:)
     integer,        allocatable :: Rrequest(:),Srequest(:)
     integer,        allocatable :: status(:,:)
     type (gbarrier_t) :: gbarrier
     integer                     :: nlyr    ! Number of layers
     integer                     :: nbuf    ! total size of message passing buffer, includes vertical levels
     integer                     :: ndepth  ! Depth of halo 
     integer                     :: npoints ! length of edge
     integer                     :: lb,ub   ! lower and upper bound of arrays
     integer                     :: nInter, nIntra
     integer                     :: id
     integer                     :: bndry_type
     integer                     :: tag
     integer                     :: win
     integer(kind=i8)            :: winsize 
  end type EdgeBuffer_t

  type, public :: LongEdgeBuffer_t
     integer          :: nlyr
     integer          :: nbuf
     integer, allocatable :: buf(:,:)
     integer, allocatable :: receive(:,:)
  end type LongEdgeBuffer_t

  type, public :: GhostBuffer3D_t
     real (kind=r8), dimension(:,:,:,:), allocatable :: buf
     real (kind=r8), dimension(:,:,:,:), allocatable :: receive
     integer :: nlyr ! Number of layers
     integer :: nhc  ! Number of layers of ghost cells
     integer :: np   ! Number of points in a cell
     integer :: nbuf ! size of the horizontal dimension of the buffers.
     integer :: elem_size ! size of 2D array (first two dimensions of buf())
  end type GhostBuffer3D_t

 
end module edgetype_mod
