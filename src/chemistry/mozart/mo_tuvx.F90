module mo_tuvx
!----------------------------------------------------------------------
!     ... wrapper for TUV-x photolysis rate constant calculator
!----------------------------------------------------------------------

  use musica_string,          only : string_t
  use shr_kind_mod,           only : r8 => shr_kind_r8
  use tuvx_core,              only : core_t
  use tuvx_grid_from_host,    only : grid_updater_t
  use tuvx_profile_from_host, only : profile_updater_t

  implicit none

  private

  public :: tuvx_init
  public :: tuvx_get_photo_rates
  public :: tuvx_finalize

  integer, parameter :: NUM_WAVELENGTHS = 157 ! TEMPORARY FOR DEVELOPMENT

  ! Inidices for grid updaters
  integer, parameter :: NUM_GRIDS = 2             ! number of grids that CAM will update at runtime
  integer, parameter :: GRID_INDEX_HEIGHT     = 1 ! Height grid index
  integer, parameter :: GRID_INDEX_WAVELENGTH = 2 ! Wavelength grid index

  ! Indices for profile updaters
  integer, parameter :: NUM_PROFILES = 8               ! number of profiles that CAM will update at runtime
  integer, parameter :: PROFILE_INDEX_TEMPERATURE  = 1 ! Temperature profile index
  integer, parameter :: PROFILE_INDEX_ALBEDO       = 2 ! Surface albedo profile index
  integer, parameter :: PROFILE_INDEX_ET_FLUX      = 3 ! Extraterrestrial flux profile index
  integer, parameter :: PROFILE_INDEX_AIR          = 4 ! Air density profile index
  integer, parameter :: PROFILE_INDEX_O3           = 5 ! Ozone profile index
  integer, parameter :: PROFILE_INDEX_O2           = 6 ! Molecular oxygen profile index
  integer, parameter :: PROFILE_INDEX_SO2          = 7 ! Sulfur dioxide profile index
  integer, parameter :: PROFILE_INDEX_NO2          = 8 ! Nitrogen dioxide profile index

  ! TODO how should this path be set and communicated to this wrapper?
  character(len=*), parameter :: tuvx_config_path = "tuvx_config.json"

  ! TUV-x calculator for each OMP thread
  type :: tuvx_ptr
    type(core_t), pointer   :: core_ => null( )        ! TUV-x calculator
    integer                 :: n_photo_rates_          ! number of photo reactions in TUV-x
    real(r8), allocatable   :: photo_rates_(:,:)       ! photolysis rate constants
                                                       !   (vertical level, reaction) [s-1]
    type(grid_updater_t)    :: grids_(NUM_GRIDS)       ! grid updaters
    type(profile_updater_t) :: profiles_(NUM_PROFILES) ! profile updaters
    real(r8), allocatable   :: heights_(:)             ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable   :: height_values_(:)       ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable   :: height_mid_values_(:)   ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable   :: wavelength_values_(:)   ! TEMPORARY FOR DEVELOPMENT
  end type tuvx_ptr
  type(tuvx_ptr), allocatable :: tuvx_ptrs(:)

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
    use cam_logfile,            only : iulog
    use musica_assert,          only : assert_msg
    use musica_mpi,             only : musica_mpi_rank
    use musica_string,          only : string_t, to_char
    use tuvx_grid,              only : grid_t
    use tuvx_grid_warehouse,    only : grid_warehouse_t
    use tuvx_profile_warehouse, only : profile_warehouse_t
    use spmd_utils,             only : main_task => masterprocid, &
                                       is_main_task => masterproc, &
                                       mpicom

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(core_t), pointer :: core
    character, allocatable :: buffer(:)
    type(string_t) :: config_path
    class(grid_t), pointer :: height
    class(grid_warehouse_t), pointer :: cam_grids
    class(profile_warehouse_t), pointer :: cam_profiles
    integer :: pack_size, pos, i_core, i_err

    config_path = tuvx_config_path
    if( is_main_task ) call log_initialization( )

#ifndef HAVE_MPI
    call assert_msg( 113937299, is_main_task, "Multiple tasks present without " &
                     //"MPI support enabled for TUV-x")
#endif

    ! Create the set of TUV-x grids and profiles that CAM will update at runtime
    cam_grids => get_cam_grids( )
    cam_profiles => get_cam_profiles( )

    ! construct a core on the primary process and pack it onto an MPI buffer
    if( is_main_task ) then
      core => core_t( config_path, cam_grids, cam_profiles )
      pack_size = core%pack_size( mpicom )
      allocate( buffer( pack_size ) )
      pos = 0
      call core%mpi_pack( buffer, pos, mpicom )
      deallocate( core )
    end if

    ! broadcast the core data to all MPI processes
#ifdef HAVE_MPI
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
#endif

    ! unpack the core for each OMP thread on every MPI process
    allocate( tuvx_ptrs( max_threads( ) ) )
    do i_core = 1, size( tuvx_ptrs )
    associate( tuvx => tuvx_ptrs( i_core ) )
      allocate( tuvx%core_ )
      pos = 0
      call tuvx%core_%mpi_unpack( buffer, pos, mpicom )

      ! Set up map between CAM and TUV-x photolysis reactions for each thread
      call create_updaters( tuvx, cam_grids, cam_profiles )

      ! TEMPORARY FOR DEVELOPMENT
      tuvx%n_photo_rates_ = tuvx%core_%number_of_photolysis_reactions( )
      height => tuvx%core_%get_grid( "height", "km" )
      allocate( tuvx%photo_rates_( height%ncells_ + 1, tuvx%n_photo_rates_ ) )
      deallocate( height )

    end associate
    end do

    deallocate( cam_grids    )
    deallocate( cam_profiles )

  end subroutine tuvx_init

!================================================================================================

  subroutine tuvx_get_photo_rates( ncol )
!-----------------------------------------------------------------------
!
! Purpose: calculate and return photolysis rate constants
!
!-----------------------------------------------------------------------

    use cam_logfile,    only : iulog
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

    associate( tuvx => tuvx_ptrs( thread_id( ) ) )
      do i_col = 1, ncol

        ! set conditions for this column in TUV-x
        call set_heights( tuvx, i_col )
        call set_temperatures( tuvx, i_col )
        call set_surface_albedo( tuvx, i_col )
        call set_et_flux( tuvx, i_col )
        call set_radiator_profiles( tuvx, i_col )

        ! Calculate photolysis rate constants for this column
        ! TEMPORARY FOR DEVELOPMENT - fix SZA
        call tuvx%core_%run( solar_zenith_angle = 45.0_r8, &
                             photolysis_rate_constants = tuvx%photo_rates_ )

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
!================================================================================================
!
! Support functions
!
!================================================================================================
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
                     trim( to_char( max_threads( ) ) )//" threads, on thread " &
                     //trim( to_char( thread_id( ) ) )
#else
      write(iulog,*) "  - without OpenMP support"
#endif
      write(iulog,*) "  - with configuration file: '"//tuvx_config_path//"'"
    end if

  end subroutine log_initialization

!================================================================================================

  function get_cam_grids( ) result( grids )
!-----------------------------------------------------------------------
!
! Purpose: creates and loads a grid warehouse with grids that CAM will
!          update at runtime
!
!-----------------------------------------------------------------------

    use ppgrid,              only : pver ! number of vertical levels
    use tuvx_grid_from_host, only : grid_from_host_t
    use tuvx_grid_warehouse, only : grid_warehouse_t

    class(grid_warehouse_t), pointer :: grids ! collection of grids to be updated by CAM

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(grid_from_host_t), pointer :: host_grid
    type(grid_updater_t)             :: updater
    real(r8) :: wavelengths(121) ! TEMPORARY FOR DEVELOPMENT
    integer :: i_wavelength

    grids => grid_warehouse_t( )

    ! Height grid will be ... \todo figure out how height grid should translate
    ! to CAM vertical grid
    host_grid => grid_from_host_t( "height", "km", pver )
    call grids%add( host_grid )
    deallocate( host_grid )

    ! Wavelength grid wil be ... /todo figure out where to get wavelength grid
    ! from (wavelengths must be set prior to construction of the TUV-x core)
    do i_wavelength = 1, size( wavelengths )
      wavelengths( i_wavelength ) = 199.0_r8 + i_wavelength
    end do
    host_grid => grid_from_host_t( "wavelength", "nm", 120 )
    updater = grid_updater_t( host_grid )
    call updater%update( edges = wavelengths )
    call grids%add( host_grid )
    deallocate( host_grid )

  end function get_cam_grids

!================================================================================================

  function get_cam_profiles( ) result( profiles )
!-----------------------------------------------------------------------
!
! Purpose: creates and loads a profile warehouse with profiles that CAM
!          will update at runtime
!
!-----------------------------------------------------------------------

    use ppgrid,                 only : pver ! number of vertical levels
    use tuvx_profile_from_host, only : profile_from_host_t
    use tuvx_profile_warehouse, only : profile_warehouse_t

    class(profile_warehouse_t), pointer :: profiles ! collection of profiles to be updated by CAM

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(profile_from_host_t), pointer :: host_profile

    profiles => profile_warehouse_t( )

    ! Temperature profile on height grid
    host_profile => profile_from_host_t( "temperature", "K", pver )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Surface albedo on wavelength grid
    host_profile => profile_from_host_t( "surface albedo", "none", NUM_WAVELENGTHS )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Extraterrestrial flux on wavelength grid
    host_profile => profile_from_host_t( "extraterrestrial flux", "photon cm-2 s-1", &
                                         NUM_WAVELENGTHS )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Air profile
    host_profile => profile_from_host_t( "air", "molecule cm-3", pver )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! O3 profile
    ! TODO optionally include if available
    host_profile => profile_from_host_t( "O3", "molecule cm-3", pver )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! O2 profile
    host_profile => profile_from_host_t( "O2", "molecule cm-3", pver )
    call profiles%add( host_profile )
    deallocate( host_profile )

  end function get_cam_profiles

!================================================================================================

  subroutine create_updaters( this, grids, profiles )
!-----------------------------------------------------------------------
!
! Purpose: creates updaters for each grid and profile that CAM will use
!          to update TUV-x at each timestep
!
!-----------------------------------------------------------------------

    use musica_assert,          only : assert
    use tuvx_grid_from_host,    only : grid_updater_t
    use tuvx_grid,              only : grid_t
    use tuvx_grid_warehouse,    only : grid_warehouse_t
    use tuvx_profile,           only : profile_t
    use tuvx_profile_from_host, only : profile_updater_t
    use tuvx_profile_warehouse, only : profile_warehouse_t

    class(tuvx_ptr),            intent(inout) :: this
    class(grid_warehouse_t),    intent(in)    :: grids
    class(profile_warehouse_t), intent(in)    :: profiles

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(grid_t),    pointer :: host_grid
    class(profile_t), pointer :: host_profile
    logical                   :: found

    ! Grid updaters

    host_grid => grids%get_grid( "height", "km" )
    this%grids_( GRID_INDEX_HEIGHT ) = this%core_%get_updater( host_grid, found )
    call assert( 213798815, found )
    allocate( this%heights_(           host_grid%size( ) + 1 ) ) ! TEMPORARY FOR DEVELOPMENT
    allocate( this%height_values_(     host_grid%size( ) + 1 ) ) ! TEMPORARY FOR DEVELOPMENT
    allocate( this%height_mid_values_( host_grid%size( )     ) ) ! TEMPORARY FOR DEVELOPMENT
    deallocate( host_grid )

    ! wavelength grid cannot be updated at runtime
    allocate( this%wavelength_values_( NUM_WAVELENGTHS + 1 ) ) ! TEMPORARY FOR DEVELOPMENT

    ! Profile updaters

    host_profile => profiles%get_profile( "temperature", "K" )
    this%profiles_( PROFILE_INDEX_TEMPERATURE ) = this%core_%get_updater( host_profile, found )
    call assert( 418735162, found )
    deallocate( host_profile )

    host_profile => profiles%get_profile( "surface albedo", "none" )
    this%profiles_( PROFILE_INDEX_ALBEDO ) = this%core_%get_updater( host_profile, found )
    call assert( 720785186, found )
    deallocate( host_profile )

    host_profile => profiles%get_profile( "extraterrestrial flux", "photon cm-2 s-1" )
    this%profiles_( PROFILE_INDEX_ET_FLUX ) = this%core_%get_updater( host_profile, found )
    call assert( 550628282, found )
    deallocate( host_profile )

    host_profile => profiles%get_profile( "air", "molecule cm-3" )
    this%profiles_( PROFILE_INDEX_AIR ) = this%core_%get_updater( host_profile, found )
    call assert( 380471378, found )
    deallocate( host_profile )

    host_profile => profiles%get_profile( "O3", "molecule cm-3" )
    this%profiles_( PROFILE_INDEX_O3 ) = this%core_%get_updater( host_profile, found )
    call assert( 210314474, found )
    deallocate( host_profile )

    host_profile => profiles%get_profile( "O2", "molecule cm-3" )
    this%profiles_( PROFILE_INDEX_O2 ) = this%core_%get_updater( host_profile, found )
    call assert( 105165970, found )
    deallocate( host_profile )

  end subroutine create_updaters

!================================================================================================

  subroutine set_heights( this, i_col )
!-----------------------------------------------------------------------
!
! Purpose: sets the height values in TUV-x for the given column
!
! TODO: Describe how CAM vertical profile is mapped to TUV-x heights
!
!-----------------------------------------------------------------------

    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT
    integer :: i_level
    do i_level = 1, size( this%heights_ )
      this%heights_( i_level ) = i_level * 1.0_r8 - 1.0_r8
    end do
    call this%grids_( GRID_INDEX_HEIGHT )%update( &
        edges = this%heights_(:) )

  end subroutine set_heights

!================================================================================================

  subroutine set_temperatures( this, i_col )
!-----------------------------------------------------------------------
!
! Purpose: sets the temperatures in TUV-x for the given column
!
! TODO: Describe how CAM temperature profile is mapped to TUV-x heights
!
!-----------------------------------------------------------------------

    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT
    this%height_values_(:) = 298.15_r8
    call this%profiles_( PROFILE_INDEX_TEMPERATURE )%update( &
        edge_values = this%height_values_(:) )

  end subroutine set_temperatures

!================================================================================================

  subroutine set_surface_albedo( this, i_col )
!-----------------------------------------------------------------------
!
! Purpose: sets the surface albedo in TUV-x for the given column
!
! TODO: Describe how CAM surface albedo profile is mapped to TUV-x wavelengths
!
!-----------------------------------------------------------------------

    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT
    this%wavelength_values_(:) = 0.1_r8
    call this%profiles_( PROFILE_INDEX_ALBEDO )%update( &
        edge_values = this%wavelength_values_(:) )

  end subroutine set_surface_albedo

!================================================================================================

  subroutine set_et_flux( this, i_col )
!-----------------------------------------------------------------------
!
! Purpose: sets the extraterrestrial flux in TUV-x for the given column
!
! TODO: Describe how CAM extraterrestrial flux profile is mapped to TUV-x wavelengths
!
!-----------------------------------------------------------------------

    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT
    this%wavelength_values_(:) = 1000.0_r8
    call this%profiles_( PROFILE_INDEX_ET_FLUX )%update( &
        edge_values = this%wavelength_values_(:) )

  end subroutine set_et_flux

!================================================================================================

  subroutine set_radiator_profiles( this, i_col )
!-----------------------------------------------------------------------
!
! Purpose: sets the profiles of optically active atmospheric constituents
!          in TUV-x for the given column
!
! TODO: Describe how CAM profiles are mapped to TUV-x heights
!
!-----------------------------------------------------------------------

    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT - air
    this%height_values_(:) = 2.54e19_r8
    this%height_mid_values_(:) = 2.54e21_r8
    call this%profiles_( PROFILE_INDEX_AIR )%update( &
        mid_point_values = this%height_values_(1:size(this%height_values_)-1), &
        edge_values = this%height_values_(:), &
        layer_densities = this%height_mid_values_(:) )

    ! TEMPORARY FOR DEVELOPMENT - O2
    this%height_values_(:) = 1.0e17_r8
    this%height_mid_values_(:) = 1.0e19_r8
    call this%profiles_( PROFILE_INDEX_O2 )%update( &
        mid_point_values = this%height_values_(1:size(this%height_values_)-1), &
        edge_values = this%height_values_(:), &
        layer_densities = this%height_mid_values_(:) )

    ! TEMPORARY FOR DEVELOPMENT - O3
    this%height_values_(:) = 1.0e13_r8
    this%height_mid_values_(:) = 1.0e15_r8
    call this%profiles_( PROFILE_INDEX_O3 )%update( &
        mid_point_values = this%height_values_(1:size(this%height_values_)-1), &
        edge_values = this%height_values_(:), &
        layer_densities = this%height_mid_values_(:) )

  end subroutine set_radiator_profiles

!================================================================================================

end module mo_tuvx
