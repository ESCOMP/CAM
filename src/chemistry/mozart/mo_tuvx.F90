module mo_tuvx
  !----------------------------------------------------------------------
  !     ... wrapper for TUV-x photolysis rate constant calculator
  !----------------------------------------------------------------------

#ifdef _MPI
  use mpi
#endif
  use tuvx_core,        only : core_t

  implicit none

  private

  public :: tuvx_init
  public :: tuvx_finalize

  ! TUV-x calculator for each OMP thread
  type :: tuvx_ptr
    type(core_t), pointer :: core_ => null( )
  end type tuvx_ptr
  type(tuvx_ptr), allocatable :: tuvx_ptrs(:)

  ! TODO where should this path be stored?
  character(len=*), parameter :: tuvx_config_path = "tuvx_config.json"
  ! TODO how to know what MPI communicator to use?
#ifdef _MPI
  integer, parameter :: tuvx_comm = MPI_COMM_WORLD
#else
  integer, parameter :: tuvx_comm = 0
#endif

!================================================================================================
contains
!================================================================================================

  subroutine tuvx_init( )
!-----------------------------------------------------------------------
!
! Purpose: initialize TUV-x for photolysis calculations
!
!-----------------------------------------------------------------------

    use cam_logfile,    only : iulog
    use musica_string,  only : string_t, to_char
    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(core_t), pointer :: core
    character, allocatable :: buffer(:)
    type(string_t) :: config_path
    integer :: pack_size, pos, i_core, i_err
    character(len=255) :: cwd

    call getcwd(cwd)
    config_path = tuvx_config_path
    if( is_main_task ) then
      write(iulog,*) "Initializing TUV-x on MPI Task "//trim( to_char( main_task ) ) &
                     //" and OpenMP thread "//trim( to_char( thread_id( ) ) ) &
                     //" for "//trim( to_char( max_threads( ) ) )//" threads."
      write(iulog,*) "TUV-x working dir: '"//trim(cwd) &
                      //"' config path: '"//config_path//"'"
    end if

#if 0
    ! construct a core on the primary process and pack it onto an MPI buffer
    if( is_main_task ) then
      core => core_t( config_path )
      pack_size = core%pack_size( tuvx_comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call core%mpi_pack( buffer, pos, tuvx_comm )
      deallocate( core )
    end if

    ! broadcast the core data to all MPI processes
#ifdef _MPI
    call mpi_bcast( pack_size, 1, MPI_INTEGER, main_task, tuvx_comm, i_err )
    if( .not. is_main_task ) allocate( buffer( pack_size ) )
    call mpi_bcast( buffer, pack_size, MPI_CHARACTER, main_task, tuvx_comm, i_err )
#endif

    ! unpack the core for each OMP thread on every MPI process
    allocate( tuvx_ptrs( max_threads( ) ) )
    do i_core = 1, size( tuvx_ptrs )
    associate( tuvx => tuvx_ptrs( i_core ) )
      allocate( tuvx%core_ )
      pos = 0
      call tuvx%core_%mpi_unpack( buffer, pos, tuvx_comm )
    end associate
    end do
#endif

  end subroutine tuvx_init

!================================================================================================

  subroutine tuvx_finalize( )
!-----------------------------------------------------------------------
!
! Purpose: clean up memory associated with TUV-x calculators
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: i_core

    if( allocated( tuvx_ptrs ) ) then
      do i_core = 1, size( tuvx_ptrs )
        if( associated( tuvx_ptrs( i_core )%core_ ) ) then
          deallocate( tuvx_ptrs( i_core )%core_ )
        end if
      end do
    end if

  end subroutine tuvx_finalize

!================================================================================================

  integer function thread_id( )
!-----------------------------------------------------------------------
!
! Purpose: returns the id of the current OpenMP thread, or 1 if not
!          using OpenMP
!
!-----------------------------------------------------------------------
#ifdef _OPENMP
    use omp_lib,        only : omp_get_thread_num
#endif

#ifdef _OPENMP
    thread_id = omp_get_thread_num( )
#else
    thread_id = 1
#endif

  end function thread_id

!================================================================================================

  integer function max_threads( )
!-----------------------------------------------------------------------
!
! Purpose: returns the number of threads available for calculations at
!          runtime
!
!-----------------------------------------------------------------------
#ifdef _OPENMP
    use omp_lib,        only : omp_get_max_threads
#endif

#ifdef _OPENMP
    max_threads = omp_get_max_threads( )
#else
    max_threads = 1
#endif

  end function max_threads

!================================================================================================

end module mo_tuvx
