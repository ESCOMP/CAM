module cam_shmem_mod
!-------------------------------------------------------------------------------
! Per-node MPI-3 shared-memory allocation helper for large, read-only lookup
! tables that are identical on every MPI rank (e.g. photolysis cross-section /
! radiative-source-function tables, aerosol-optics property tables).  Instead of
! every rank holding its own copy, one physical copy is allocated per
! shared-memory node and mapped (read-only) into every rank on that node.  This
! frees (ranks_per_node - 1) copies per node.
!
! Usage pattern (collective over the global communicator):
!   1. call cam_shmem_alloc_rX_Nd(ptr, win, dims...)   ! all ranks
!   2. call cam_shmem_fence(win)                        ! open epoch
!   3. node leaders fill the table (read from file / broadcast among leaders)
!   4. call cam_shmem_fence(win)                        ! publish; all ranks may now read
! After step 4 the data is static and may be read for the rest of the run with no
! further synchronization.  Free it with cam_shmem_free(ptr, win) at finalize
! (collective over the node communicator); otherwise MPI_Finalize reclaims it.
!
! This module uses the F90 'mpi' module (not mpif.h).  MPI-3.0 makes the TYPE(C_PTR)
! overloads of MPI_WIN_ALLOCATE_SHARED / MPI_WIN_SHARED_QUERY optional in the 'mpi'
! module, and some implementations (cray-mpich) do not expose those two routines
! there at all, so explicit interfaces for them are declared below instead.
!-------------------------------------------------------------------------------

   use shr_kind_mod,   only: r4 => shr_kind_r4, r8 => shr_kind_r8
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use spmd_utils,     only: masterproc

#ifdef SPMD
   use spmd_utils,     only: mpicom
   use mpi,            only: MPI_ADDRESS_KIND, MPI_COMM_NULL, MPI_COMM_TYPE_SHARED,   &
                             MPI_INFO_NULL, MPI_SUCCESS, MPI_UNDEFINED, MPI_WIN_NULL, &
                             mpi_comm_rank, mpi_comm_size, mpi_comm_split,            &
                             mpi_comm_split_type, mpi_win_fence, mpi_win_free,        &
                             MPI_INTEGER, MPI_MIN
   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_associated
#endif

   implicit none
   private

   public :: cam_shmem_alloc_r4_4d   ! allocate node-shared real(r4) rank-4 table
   public :: cam_shmem_alloc_r4_5d   ! allocate node-shared real(r4) rank-5 table
   public :: cam_shmem_alloc_r8_2d   ! allocate node-shared real(r8) rank-2 table
   public :: cam_shmem_alloc_r8_3d   ! allocate node-shared real(r8) rank-3 table
   public :: cam_shmem_alloc_r8_4d   ! allocate node-shared real(r8) rank-4 table
   public :: cam_shmem_alloc_r8_5d   ! allocate node-shared real(r8) rank-5 table
   public :: cam_shmem_free          ! free a node-shared table (MPI_Win_free)
   public :: cam_shmem_fence         ! synchronize a window (publish writes)
   public :: cam_shmem_is_leader     ! .true. on the leader (rank 0) of this node
   public :: cam_shmem_leader_comm   ! communicator containing only node leaders
   public :: cam_shmem_npes_per_node ! number of ranks sharing this node
   public :: cam_shmem_init          ! force one-time node-comm setup at an early collective point
   public :: cam_shmem_readnl        ! read the cam_shmem_nl namelist

   ! Generic finalizer: free whatever node-shared table the pointer aliases.  Safe
   ! to call blindly on an unallocated table (win == -1 / null pointer -> no-op).
   interface cam_shmem_free
      module procedure cam_shmem_free_r4_4d
      module procedure cam_shmem_free_r4_5d
      module procedure cam_shmem_free_r8_2d
      module procedure cam_shmem_free_r8_3d
      module procedure cam_shmem_free_r8_4d
      module procedure cam_shmem_free_r8_5d
   end interface cam_shmem_free

   integer, parameter :: bytes_r4 = 4
   integer, parameter :: bytes_r8 = 8

#ifdef SPMD
   ! MPI_WIN_ALLOCATE_SHARED and MPI_WIN_SHARED_QUERY take a TYPE(C_PTR) baseptr in
   ! MPI-3, but that overload is not exposed by every implementation's F90 'mpi'
   ! module (cray-mpich does not expose these two routines at all), so they cannot be
   ! named in the 'only' list above.  Declare explicit interfaces here rather than
   ! relying on implicit ones, so the argument types are still checked.  Both forms
   ! resolve to the same mpi_win_allocate_shared_ / mpi_win_shared_query_ bindings,
   ! so this is equivalent on implementations that do expose them.
   interface
      subroutine mpi_win_allocate_shared(nbytes, disp_unit, info, comm, baseptr, win, ierror)
         import :: MPI_ADDRESS_KIND, c_ptr
         integer(kind=MPI_ADDRESS_KIND), intent(in)  :: nbytes
         integer,                        intent(in)  :: disp_unit
         integer,                        intent(in)  :: info
         integer,                        intent(in)  :: comm
         type(c_ptr),                    intent(out) :: baseptr
         integer,                        intent(out) :: win
         integer,                        intent(out) :: ierror
      end subroutine mpi_win_allocate_shared

      subroutine mpi_win_shared_query(win, rank, segsize, disp_unit, baseptr, ierror)
         import :: MPI_ADDRESS_KIND, c_ptr
         integer,                        intent(in)  :: win
         integer,                        intent(in)  :: rank
         integer(kind=MPI_ADDRESS_KIND), intent(out) :: segsize
         integer,                        intent(out) :: disp_unit
         type(c_ptr),                    intent(out) :: baseptr
         integer,                        intent(out) :: ierror
      end subroutine mpi_win_shared_query
   end interface

   logical, save :: initialized = .false.
   integer, save :: node_comm   = MPI_COMM_NULL  ! ranks sharing a node
   integer, save :: leader_comm = MPI_COMM_NULL  ! one rank per node (the leaders)
   integer, save :: node_rank   = 0
   integer, save :: node_size   = 1
   logical, save :: is_leader   = .true.
#else
   ! Non-SPMD build: single task, nothing is shared.
   integer, save :: node_size   = 1
   logical, save :: is_leader   = .true.
#endif

   ! Namelist switch (cam_shmem_nl): .true. disables node-shared memory, so every
   ! rank holds private copies.  Answers must not change; the SHM system test
   ! compares a run with this set against a default run.
   logical, save :: cam_shmem_force_private = .false.

contains

!===============================================================================

#ifdef SPMD
   subroutine init_comms()
      ! Lazily build the node-local and node-leader communicators.  Collective
      ! over mpicom; safe to call from every shared-memory allocation request.
      integer :: ierr, color, wrank, nsize0
      logical :: shared_ok

      if (initialized) return

      call mpi_comm_split_type(mpicom, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &
                               node_comm, ierr)
      call mpi_comm_rank(node_comm, node_rank, ierr)
      call mpi_comm_size(node_comm, node_size, ierr)
      is_leader = (node_rank == 0)
      nsize0 = node_size

      ! Verify that MPI-3 shared windows really work before relying on them: some
      ! MPI setups (e.g. MVAPICH2 under CESM on izumi) build node_comm correctly but
      ! return windows that are not shared-memory backed, so peers get no address
      ! from MPI_Win_shared_query.  If any rank fails the probe, or the namelist
      ! switch asks for it, every rank becomes its own "node": node_size 1, every
      ! rank a leader, leader_comm = all ranks.  Callers then hold private copies
      ! and their leader broadcasts reach every rank, so answers are unchanged.
      shared_ok = probe_shared_windows()
      if (cam_shmem_force_private .or. .not. shared_ok) then
         call mpi_comm_free(node_comm, ierr)
         call mpi_comm_rank(mpicom, wrank, ierr)
         call mpi_comm_split(mpicom, wrank, 0, node_comm, ierr)
         node_rank = 0
         node_size = 1
         is_leader = .true.
      end if
      if (masterproc) then
         if (shared_ok) then
            write(iulog,'(a,i0,a)') 'cam_shmem_mod: MPI-3 shared-memory probe PASSED (', nsize0, &
                 ' ranks on the first node)'
         else
            write(iulog,'(a)') 'cam_shmem_mod: MPI-3 shared-memory probe FAILED (peers got no ' // &
                 'address from MPI_Win_shared_query)'
         end if
         if (cam_shmem_force_private) then
            write(iulog,'(a)') 'cam_shmem_mod: cam_shmem_force_private = .true.: using per-rank ' // &
                 'private memory'
         else if (shared_ok) then
            write(iulog,'(a)') 'cam_shmem_mod: node-shared tables enabled'
         else
            write(iulog,'(a)') 'cam_shmem_mod: falling back to per-rank private memory'
         end if
      end if

      ! Communicator of node leaders only (used to distribute data between nodes).
      ! masterproc (global rank 0) is a node leader and is rank 0 here.
      if (is_leader) then
         color = 0
      else
         color = MPI_UNDEFINED
      end if
      call mpi_comm_split(mpicom, color, 0, leader_comm, ierr)

      initialized = .true.
   end subroutine init_comms

!===============================================================================

   logical function probe_shared_windows()
      ! Collective over mpicom: .true. on every rank iff, on every node, a small
      ! window allocated by the node leader is visible to all its peers.
      integer(kind=MPI_ADDRESS_KIND) :: wsize, qsize
      type(c_ptr) :: base
      integer :: win, ierr, qdisp, lok, gok

      wsize = 0_MPI_ADDRESS_KIND
      if (is_leader) wsize = 64_MPI_ADDRESS_KIND
      lok = 1
      call mpi_win_allocate_shared(wsize, 8, MPI_INFO_NULL, node_comm, base, win, ierr)
      if (ierr /= MPI_SUCCESS) then
         lok = 0
      else
         if (.not. is_leader) then
            call mpi_win_shared_query(win, 0, qsize, qdisp, base, ierr)
            if (ierr /= MPI_SUCCESS .or. qsize < 64_MPI_ADDRESS_KIND) lok = 0
         end if
         if (.not. c_associated(base)) lok = 0
         call mpi_win_free(win, ierr)
      end if
      call mpi_allreduce(lok, gok, 1, MPI_INTEGER, MPI_MIN, mpicom, ierr)
      probe_shared_windows = (gok == 1)
   end function probe_shared_windows

!===============================================================================

   subroutine shmem_alloc_bytes(nbytes, disp_unit, win, baseptr)
      ! Allocate a node-shared window of nbytes (only the node leader requests
      ! storage; peers map the leader's segment) and return the base C pointer of
      ! the leader's contiguous segment.  Type/shape-agnostic core shared by all
      ! of the typed cam_shmem_alloc_* wrappers.
      integer(kind=MPI_ADDRESS_KIND), intent(in)  :: nbytes
      integer,                        intent(in)  :: disp_unit
      integer,                        intent(out) :: win
      type(c_ptr),                    intent(out) :: baseptr

      integer(kind=MPI_ADDRESS_KIND) :: winsize, qsize
      integer :: ierr, qdisp

      call init_comms()

      if (is_leader) then
         winsize = nbytes
      else
         winsize = 0_MPI_ADDRESS_KIND
      end if

      call mpi_win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, node_comm, &
                                   baseptr, win, ierr)
      if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_allocate_shared failed')

      ! Non-leaders learn the address of the leader's (rank 0) contiguous segment.
      if (.not. is_leader) then
         call mpi_win_shared_query(win, 0, qsize, qdisp, baseptr, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_shared_query failed')
      end if
   end subroutine shmem_alloc_bytes
#endif

!===============================================================================

   subroutine cam_shmem_readnl(nlfile)
      ! Read namelist group cam_shmem_nl (cam_shmem_force_private).  Must be called
      ! before cam_shmem_init / the first cam_shmem_alloc_* call.
      use namelist_utils, only: find_group_name
#ifdef SPMD
      use spmd_utils,     only: masterprocid, mpi_logical
#endif
      character(len=*), intent(in) :: nlfile
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'cam_shmem_readnl'
      namelist /cam_shmem_nl/ cam_shmem_force_private

      if (masterproc) then
         open(newunit=unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'cam_shmem_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, cam_shmem_nl, iostat=ierr)
            if (ierr /= 0) call endrun(subname//': ERROR reading namelist cam_shmem_nl')
         end if
         close(unitn)
      end if
#ifdef SPMD
      call mpi_bcast(cam_shmem_force_private, 1, mpi_logical, masterprocid, mpicom, ierr)
#endif
      if (masterproc) write(iulog,*) subname//': cam_shmem_force_private = ', cam_shmem_force_private
   end subroutine cam_shmem_readnl

!===============================================================================

   subroutine cam_shmem_init()
      ! Force the one-time node-local / node-leader communicator setup NOW, at a
      ! controlled early all-ranks collective point, instead of lazily on the first
      ! cam_shmem_alloc_* during physics init.  At large rank counts the first
      ! MPI_Comm_split_type(MPI_COMM_TYPE_SHARED) can be expensive; doing it early
      ! (rather than mid-physprop) keeps it predictable.  Idempotent; non-SPMD no-op.
#ifdef SPMD
      call init_comms()
#endif
   end subroutine cam_shmem_init

!===============================================================================

   subroutine cam_shmem_alloc_r4_4d(ptr, win, n1, n2, n3, n4)
      real(r4), pointer, intent(out) :: ptr(:,:,:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2, n3, n4
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(n3,MPI_ADDRESS_KIND)*int(n4,MPI_ADDRESS_KIND) &
              *int(bytes_r4,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r4, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2, n3, n4])
#else
      win = -1
      allocate(ptr(n1,n2,n3,n4))
#endif
   end subroutine cam_shmem_alloc_r4_4d

!===============================================================================

   subroutine cam_shmem_alloc_r4_5d(ptr, win, n1, n2, n3, n4, n5)
      real(r4), pointer, intent(out) :: ptr(:,:,:,:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2, n3, n4, n5
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(n3,MPI_ADDRESS_KIND)*int(n4,MPI_ADDRESS_KIND) &
              *int(n5,MPI_ADDRESS_KIND)*int(bytes_r4,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r4, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2, n3, n4, n5])
#else
      win = -1
      allocate(ptr(n1,n2,n3,n4,n5))
#endif
   end subroutine cam_shmem_alloc_r4_5d

!===============================================================================

   subroutine cam_shmem_alloc_r8_2d(ptr, win, n1, n2)
      real(r8), pointer, intent(out) :: ptr(:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(bytes_r8,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r8, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2])
#else
      win = -1
      allocate(ptr(n1,n2))
#endif
   end subroutine cam_shmem_alloc_r8_2d

!===============================================================================

   subroutine cam_shmem_alloc_r8_3d(ptr, win, n1, n2, n3)
      real(r8), pointer, intent(out) :: ptr(:,:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2, n3
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(n3,MPI_ADDRESS_KIND)*int(bytes_r8,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r8, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2, n3])
#else
      win = -1
      allocate(ptr(n1,n2,n3))
#endif
   end subroutine cam_shmem_alloc_r8_3d

!===============================================================================

   subroutine cam_shmem_alloc_r8_4d(ptr, win, n1, n2, n3, n4)
      real(r8), pointer, intent(out) :: ptr(:,:,:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2, n3, n4
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(n3,MPI_ADDRESS_KIND)*int(n4,MPI_ADDRESS_KIND) &
              *int(bytes_r8,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r8, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2, n3, n4])
#else
      win = -1
      allocate(ptr(n1,n2,n3,n4))
#endif
   end subroutine cam_shmem_alloc_r8_4d

!===============================================================================

   subroutine cam_shmem_alloc_r8_5d(ptr, win, n1, n2, n3, n4, n5)
      real(r8), pointer, intent(out) :: ptr(:,:,:,:,:)
      integer,           intent(out) :: win
      integer,           intent(in)  :: n1, n2, n3, n4, n5
#ifdef SPMD
      integer(kind=MPI_ADDRESS_KIND) :: nbytes
      type(c_ptr) :: baseptr
      nbytes = int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
              *int(n3,MPI_ADDRESS_KIND)*int(n4,MPI_ADDRESS_KIND) &
              *int(n5,MPI_ADDRESS_KIND)*int(bytes_r8,MPI_ADDRESS_KIND)
      call shmem_alloc_bytes(nbytes, bytes_r8, win, baseptr)
      call c_f_pointer(baseptr, ptr, [n1, n2, n3, n4, n5])
#else
      win = -1
      allocate(ptr(n1,n2,n3,n4,n5))
#endif
   end subroutine cam_shmem_alloc_r8_5d

!===============================================================================

   subroutine cam_shmem_fence(win)
      ! Collective over the node communicator; synchronizes the window so that
      ! stores by the leader become visible to all ranks on the node.
      integer, intent(in) :: win
#ifdef SPMD
      integer :: ierr
      call mpi_win_fence(0, win, ierr)
      if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_fence failed')
#endif
   end subroutine cam_shmem_fence

!===============================================================================
! cam_shmem_free generic procedures.  Each frees the node-shared window (in SPMD)
! and disassociates the pointer, or deallocates it (non-SPMD).  Collective over
! the node communicator in SPMD; a no-op when win == -1 (table never shared).
!===============================================================================

   subroutine cam_shmem_free_r4_4d(ptr, win)
      real(r4), pointer      :: ptr(:,:,:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r4_4d

   subroutine cam_shmem_free_r4_5d(ptr, win)
      real(r4), pointer      :: ptr(:,:,:,:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r4_5d

   subroutine cam_shmem_free_r8_2d(ptr, win)
      real(r8), pointer      :: ptr(:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r8_2d

   subroutine cam_shmem_free_r8_3d(ptr, win)
      real(r8), pointer      :: ptr(:,:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r8_3d

   subroutine cam_shmem_free_r8_4d(ptr, win)
      real(r8), pointer      :: ptr(:,:,:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r8_4d

   subroutine cam_shmem_free_r8_5d(ptr, win)
      real(r8), pointer      :: ptr(:,:,:,:,:)
      integer, intent(inout) :: win
#ifdef SPMD
      integer :: ierr
      if (win /= -1 .and. win /= MPI_WIN_NULL) then
         call mpi_win_free(win, ierr)
         if (ierr /= MPI_SUCCESS) call endrun('cam_shmem_mod: MPI_Win_free failed')
      end if
      if (associated(ptr)) nullify(ptr)
#else
      if (associated(ptr)) deallocate(ptr)
#endif
      win = -1
   end subroutine cam_shmem_free_r8_5d

!===============================================================================

   logical function cam_shmem_is_leader()
#ifdef SPMD
      call init_comms()
#endif
      cam_shmem_is_leader = is_leader
   end function cam_shmem_is_leader

!===============================================================================

   integer function cam_shmem_leader_comm()
#ifdef SPMD
      call init_comms()
      cam_shmem_leader_comm = leader_comm
#else
      cam_shmem_leader_comm = -1
#endif
   end function cam_shmem_leader_comm

!===============================================================================

   integer function cam_shmem_npes_per_node()
#ifdef SPMD
      call init_comms()
#endif
      cam_shmem_npes_per_node = node_size
   end function cam_shmem_npes_per_node

!===============================================================================

end module cam_shmem_mod
