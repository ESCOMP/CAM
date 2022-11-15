module mo_tuvx
!----------------------------------------------------------------------
!     ... wrapper for TUV-x photolysis rate constant calculator
!----------------------------------------------------------------------

  use musica_string,    only : string_t
  use shr_kind_mod,     only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use tuvx_core,        only : core_t
  use cam_logfile,      only : iulog
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun

  implicit none

  private

  public :: tuvx_readnl
  public :: tuvx_init
  public :: tuvx_get_photo_rates
  public :: tuvx_finalize
  public :: tuvx_active

  ! TUV-x calculator for each OMP thread
  type :: tuvx_ptr
    type(core_t), pointer :: core_ => null( )  ! TUV-x calculator
    integer               :: n_photo_rates_    ! number of photo reactions in TUV-x
    real(r8), allocatable :: photo_rates_(:,:) ! photolysis rate constants
                                               !   (vertical level, reaction) [s-1]
  end type tuvx_ptr
  type(tuvx_ptr), allocatable :: tuvx_ptrs(:)

  ! namelist options
  character(len=cl) :: tuvx_config_path = 'NONE'  ! absolute path to TUVX configuration file
  logical, protected :: tuvx_active = .false.

!================================================================================================
contains
!================================================================================================

  !-----------------------------------------------------------------------
  ! read namelist options
  !-----------------------------------------------------------------------
  subroutine tuvx_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_logical

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'tuvx_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /tuvx_opts/ tuvx_active, tuvx_config_path

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'tuvx_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, tuvx_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(tuvx_config_path, len(tuvx_config_path), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(tuvx_active,      1,                     mpi_logical,   masterprocid, mpicom, ierr)

    if (tuvx_active .and. tuvx_config_path == 'NONE') then
       call endrun(subname // ' : must set tuvx_config_path when TUV-X is active')
    end if

    if (masterproc) then
       write(iulog,*) 'tuvx_readnl: tuvx_config_path = ', trim(tuvx_config_path)
       write(iulog,*) 'tuvx_readnl: tuvx_active = ', tuvx_active
    end if

  end subroutine tuvx_readnl
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  subroutine tuvx_init( )
!-----------------------------------------------------------------------
!
! Purpose: initialize TUV-x for photolysis calculations
!
!-----------------------------------------------------------------------

    use cam_logfile,    only : iulog
    use musica_assert,  only : assert_msg
    use musica_mpi,     only : musica_mpi_rank
    use musica_string,  only : string_t, to_char
    use tuvx_grid,      only : grid_t
    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc, &
                               mpicom, mpi_character, mpi_integer, mpi_success

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(core_t), pointer :: core
    character, allocatable :: buffer(:)
    type(string_t) :: config_path
    class(grid_t), pointer :: height
    integer :: pack_size, pos, i_core, i_err, i_thread

    if (.not.tuvx_active) return

    config_path = trim(tuvx_config_path)

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
    call mpi_bcast( pack_size, 1, MPI_INTEGER, main_task, mpicom, i_err )
    if( i_err /= MPI_SUCCESS ) then
      write(iulog,*) "TUV-x MPI int bcast error"
      call mpi_abort( mpicom, 1, i_err )
    end if
    if( .not. is_main_task ) allocate( buffer( pack_size ) )
    call mpi_bcast( buffer, pack_size, MPI_CHARACTER, main_task, mpicom, i_err )
    if( i_err /= MPI_SUCCESS ) then
      write(iulog,*) "TUV-x MPI char array bcast error"
      call mpi_abort( mpicom, 1, i_err )
    end if

    ! unpack the core for each OMP thread on every MPI process
    allocate( tuvx_ptrs( max_threads( ) ) )
    do i_core = 1, size( tuvx_ptrs )
    associate( tuvx => tuvx_ptrs( i_core ) )
      allocate( tuvx%core_ )
      pos = 0
      call tuvx%core_%mpi_unpack( buffer, pos, mpicom )
    end associate
    end do

    ! Set up map between CAM and TUV-x photolysis reactions for each thread
    do i_thread = 1, max_threads( )
    associate( tuvx => tuvx_ptrs( i_thread ) )
      tuvx%n_photo_rates_ = tuvx%core_%number_of_photolysis_reactions( )
      height => tuvx%core_%get_grid( "height", "km" )
      ! Temporary for development
      allocate( tuvx%photo_rates_( height%ncells_ + 1, tuvx%n_photo_rates_ ) )
    end associate
    end do

  end subroutine tuvx_init

!================================================================================================

  subroutine tuvx_get_photo_rates( ncol )
!-----------------------------------------------------------------------
!
! Purpose: calculate and return photolysis rate constants
!
!-----------------------------------------------------------------------

    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc, &
                               mpicom
!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------

    integer, intent(in) :: ncol ! Number of colums to calculated photolysis for

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: i_col                                ! column index

    if (.not.tuvx_active) return

    associate( tuvx => tuvx_ptrs( thread_id( ) ) )
      do i_col = 1, ncol
        ! Temporary fix SZA for development
        call tuvx%core_%run( 45.0_r8, photolysis_rate_constants = tuvx%photo_rates_ )
      end do
    end associate

  end subroutine tuvx_get_photo_rates

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
      associate( tuvx => tuvx_ptrs( i_core ) )
        if( associated( tuvx%core_ ) ) deallocate( tuvx%core_ )
      end associate
      end do
    end if

  end subroutine tuvx_finalize

!================================================================================================

  integer function thread_id( )
!-----------------------------------------------------------------------
!
! Purpose: returns the id of the current OpenMP thread, or 1 if not
!          using OpenMP (1 <= id <= max_threads())
!
!-----------------------------------------------------------------------
#ifdef _OPENMP
    use omp_lib,        only : omp_get_thread_num
#endif

#ifdef _OPENMP
    thread_id = 1 + omp_get_thread_num( )
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
