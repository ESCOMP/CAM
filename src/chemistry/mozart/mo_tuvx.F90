module mo_tuvx
!----------------------------------------------------------------------
!     ... wrapper for TUV-x photolysis rate constant calculator
!----------------------------------------------------------------------

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

  ! TODO where should this file and other TUV-x data be stored?
  ! TODO how should this path be set and communicated to this wrapper?
  character(len=*), parameter :: tuvx_config_path = "tuvx_config.json"

!================================================================================================
contains
!================================================================================================

  subroutine tuvx_init( )
!-----------------------------------------------------------------------
!
! Purpose: initialize TUV-x for photolysis calculations
!
!-----------------------------------------------------------------------

#ifdef HAVE_MPI
    use mpi
#endif
    use musica_assert,  only : assert_msg
    use musica_string,  only : string_t, to_char
    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc, &
                               mpicom

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(core_t), pointer :: core
    character, allocatable :: buffer(:)
    type(string_t) :: config_path
    integer :: pack_size, pos, i_core, i_err

    config_path = tuvx_config_path
    if( is_main_task ) call log_initialization( )

#ifndef HAVE_MPI
    call assert_msg( 113937299, is_main_task, "Multiple tasks present without " &
                     //"MPI support enabled for TUV-x")
#endif

    ! construct a core on the primary process and pack it onto an MPI buffer
    if( is_main_task ) then
      core => core_t( config_path )
      pack_size = core%pack_size( mpicom )
      allocate( buffer( pack_size ) )
      pos = 0
      call core%mpi_pack( buffer, pos, mpicom )
      deallocate( core )
    end if

    ! broadcast the core data to all MPI processes
#ifdef HAVE_MPI
    call mpi_bcast( pack_size, 1, MPI_INTEGER, main_task, mpicom, i_err )
    if( .not. is_main_task ) allocate( buffer( pack_size ) )
    call mpi_bcast( buffer, pack_size, MPI_CHARACTER, main_task, mpicom, i_err )
#endif

    ! unpack the core for each OMP thread on every MPI process
    allocate( tuvx_ptrs( max_threads( ) ) )
    do i_core = 1, size( tuvx_ptrs )
    associate( tuvx => tuvx_ptrs( i_core ) )
      allocate( tuvx%core_ )
      pos = 0
      call tuvx%core_%mpi_unpack( buffer, pos, mpicom )
    end associate
    end do

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

  subroutine log_initialization( )
!-----------------------------------------------------------------------
!
! Purpose: prints initialization conditions to the log file
!
!-----------------------------------------------------------------------

    use cam_logfile,    only : iulog
    use musica_string,  only : to_char
    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc

    if( is_main_task ) then
      write(iulog,*) "Initializing TUV-x"
#ifdef HAVE_MPI
      write(iulog,*) "  - with MPI support on task "//trim( to_char( main_task ) )
#else
      write(iulog,*) "  - without MPI support"
#endif
#ifdef _OPENMP
      write(iulog,*) "  - with OpenMP support for "// &
                     trim( to_char( max_threads( ) ) )//" threads, on thread" &
                     //trim( to_char( thread_id( ) ) )
#else
      write(iulog,*) "  - without OpenMP support"
#endif
      write(iulog,*) "  - with configuration file: '"//tuvx_config_path//"'"
    end if

  end subroutine log_initialization
!================================================================================================

end module mo_tuvx
