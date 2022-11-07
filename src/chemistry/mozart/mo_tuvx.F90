module mo_tuvx
!----------------------------------------------------------------------
!     ... wrapper for TUV-x photolysis rate constant calculator
!----------------------------------------------------------------------

  use musica_string,           only : string_t
  use ppgrid,                  only : pver                ! number of vertical layers
  use shr_kind_mod,            only : r8 => shr_kind_r8
  use tuvx_core,               only : core_t
  use tuvx_grid_from_host,     only : grid_updater_t
  use tuvx_profile_from_host,  only : profile_updater_t
  use tuvx_radiator_from_host, only : radiator_updater_t

  implicit none

  private

  public :: tuvx_init
  public :: tuvx_get_photo_rates
  public :: tuvx_finalize

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

  ! Indices for radiator updaters
  integer, parameter :: NUM_RADIATORS = 1          ! number of radiators that CAM will update at runtime
  integer, parameter :: RADIATOR_INDEX_AEROSOL = 1 ! Aerosol radiator index

  ! Information needed to access CAM species state data
  logical :: is_fixed_O2 = .false. ! indicates whether O2 concentrations are fixed
  logical :: is_fixed_O3 = .false. ! indicates whether O3 concentrations are fixed
  integer :: index_O2 = 0 ! index for O2 in concentration array
  integer :: index_O3 = 0 ! index for O3 in concentration array

  ! TODO how should this path be set and communicated to this wrapper?
  character(len=*), parameter :: tuvx_config_path = "tuvx_config.json"

  ! TUV-x calculator for each OMP thread
  type :: tuvx_ptr
    type(core_t), pointer    :: core_ => null( )          ! TUV-x calculator
    integer                  :: n_photo_rates_            ! number of photo reactions in TUV-x
    real(r8), allocatable    :: photo_rates_(:,:)         ! photolysis rate constants
                                                          !   (vertical level, reaction) [s-1]
    type(grid_updater_t)     :: grids_(NUM_GRIDS)         ! grid updaters
    type(profile_updater_t)  :: profiles_(NUM_PROFILES)   ! profile updaters
    type(radiator_updater_t) :: radiators_(NUM_RADIATORS) ! radiator updaters
    real(r8)                 :: height_delta_(pver)       ! change in height in each
                                                          !   vertical layer (km)
    real(r8), allocatable    :: heights_(:)               ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable    :: height_values_(:)         ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable    :: height_mid_values_(:)     ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable    :: wavelength_values_(:)     ! TEMPORARY FOR DEVELOPMENT
    real(r8), allocatable    :: optics_values_(:,:)       ! TEMPORARY FOR DEVELOPMENT
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
    use cam_logfile,             only : iulog
    use mo_chem_utls,            only : get_spc_ndx, get_inv_ndx
    use musica_assert,           only : assert_msg
    use musica_mpi,              only : musica_mpi_rank
    use musica_string,           only : string_t, to_char
    use tuvx_grid,               only : grid_t
    use tuvx_grid_warehouse,     only : grid_warehouse_t
    use tuvx_profile_warehouse,  only : profile_warehouse_t
    use tuvx_radiator_warehouse, only : radiator_warehouse_t
    use spmd_utils,              only : main_task => masterprocid, &
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
    class(radiator_warehouse_t), pointer :: cam_radiators
    integer :: pack_size, pos, i_core, i_err

    config_path = tuvx_config_path
    if( is_main_task ) call log_initialization( )

#ifndef HAVE_MPI
    call assert_msg( 113937299, is_main_task, "Multiple tasks present without " &
                     //"MPI support enabled for TUV-x")
#endif

    ! Create the set of TUV-x grids and profiles that CAM will update at runtime
    cam_grids => get_cam_grids( )
    cam_profiles => get_cam_profiles( cam_grids )
    cam_radiators => get_cam_radiators( cam_grids )

    ! construct a core on the primary process and pack it onto an MPI buffer
    if( is_main_task ) then
      core => core_t( config_path, cam_grids, cam_profiles, cam_radiators )
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
      call create_updaters( tuvx, cam_grids, cam_profiles, cam_radiators )

      ! TEMPORARY FOR DEVELOPMENT
      tuvx%n_photo_rates_ = tuvx%core_%number_of_photolysis_reactions( )
      height => tuvx%core_%get_grid( "height", "km" )
      allocate( tuvx%photo_rates_( height%ncells_ + 1, tuvx%n_photo_rates_ ) )
      deallocate( height )

    end associate
    end do

    deallocate( cam_grids     )
    deallocate( cam_profiles  )
    deallocate( cam_radiators )

    ! Get index info for CAM species concentrations
    index_O2 = get_inv_ndx( 'O2' )
    is_fixed_O2 = index_O2 > 0
    if( .not. is_fixed_O2 ) index_O2 = get_spc_ndx( 'O2' )
    index_O3 = get_inv_ndx( 'O3' )
    is_fixed_O3 = index_O3 > 0
    if( .not. is_fixed_O3 ) index_O3 = get_spc_ndx( 'O3' )


  end subroutine tuvx_init

!================================================================================================

  subroutine tuvx_get_photo_rates( ncol, height_mid, height_int, temperature_mid, &
      surface_temperature, fixed_species_conc, species_vmr, surface_albedo )
!-----------------------------------------------------------------------
!
! Purpose: calculate and return photolysis rate constants
!
!-----------------------------------------------------------------------

    use cam_logfile,    only : iulog        ! log info output unit
    use chem_mods,      only : gas_pcnst, & ! number of non-fixed species
                               nfs          ! number of fixed species
    use ppgrid,         only : pcols        ! maximum number of columns
    use spmd_utils,     only : main_task => masterprocid, &
                               is_main_task => masterproc, &
                               mpicom
!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    integer,  intent(in) :: ncol                        ! Number of colums to calculated photolysis for
    real(r8), intent(in) :: height_mid(pcols,pver)      ! height at mid-points (km)
    real(r8), intent(in) :: height_int(pcols,pver+1)    ! height at interfaces (km)
    real(r8), intent(in) :: temperature_mid(pcols,pver) ! midpoint temperature (K)
    real(r8), intent(in) :: surface_temperature(pcols)  ! surface temperature (K)
    real(r8), intent(in) :: fixed_species_conc(ncol,pver,max(1,nfs)) ! fixed species densities
                                                                     !   (molecule cm-3)
    real(r8), intent(in) :: species_vmr(ncol,pver,max(1,gas_pcnst))  ! species volume mixing
                                                                     !   ratios (mol mol-1)
    real(r8), intent(in) :: surface_albedo(pcols)       ! surface albedo (unitless)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: i_col                                ! column index

    associate( tuvx => tuvx_ptrs( thread_id( ) ) )
      do i_col = 1, ncol

        ! update grid heights
        call set_heights( tuvx, i_col, ncol, height_mid, height_int )

        ! set conditions for this column in TUV-x
        call set_temperatures( tuvx, i_col, temperature_mid, surface_temperature )
        call set_surface_albedo( tuvx, i_col, surface_albedo )
        call set_et_flux( tuvx, i_col )
        call set_radiator_profiles( tuvx, i_col, ncol, fixed_species_conc, species_vmr )

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

    use cam_logfile,    only : iulog ! log info output unit
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

    use tuvx_grid_from_host, only : grid_from_host_t
    use tuvx_grid_warehouse, only : grid_warehouse_t

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
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

  function get_cam_profiles( grids ) result( profiles )
!-----------------------------------------------------------------------
!
! Purpose: creates and loads a profile warehouse with profiles that CAM
!          will update at runtime
!
!-----------------------------------------------------------------------

    use tuvx_grid,              only : grid_t
    use tuvx_grid_warehouse,    only : grid_warehouse_t
    use tuvx_profile_from_host, only : profile_from_host_t
    use tuvx_profile_warehouse, only : profile_warehouse_t

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(grid_warehouse_t),    intent(in) :: grids    ! CAM grids used in TUV-x
    class(profile_warehouse_t), pointer    :: profiles ! collection of profiles to be updated by CAM

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(profile_from_host_t), pointer :: host_profile
    class(grid_t),              pointer :: height, wavelength

    profiles   => profile_warehouse_t( )
    height     => grids%get_grid(     "height", "km" )
    wavelength => grids%get_grid( "wavelength", "nm" )

    ! Temperature profile on height grid
    host_profile => profile_from_host_t( "temperature", "K", height%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Surface albedo on wavelength grid
    host_profile => profile_from_host_t( "surface albedo", "none", wavelength%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Extraterrestrial flux on wavelength grid
    host_profile => profile_from_host_t( "extraterrestrial flux", "photon cm-2 s-1", &
                                         wavelength%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! Air profile
    host_profile => profile_from_host_t( "air", "molecule cm-3", height%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! O3 profile
    ! TODO optionally include if available
    host_profile => profile_from_host_t( "O3", "molecule cm-3", height%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    ! O2 profile
    host_profile => profile_from_host_t( "O2", "molecule cm-3", height%size( ) )
    call profiles%add( host_profile )
    deallocate( host_profile )

    deallocate( height )
    deallocate( wavelength )

  end function get_cam_profiles

!================================================================================================

  function get_cam_radiators( grids ) result( radiators )
!-----------------------------------------------------------------------
!
! Purpose: creates and loads a radiator warehouse with radiators that CAM
!          will update at runtime
!
!-----------------------------------------------------------------------

    use tuvx_grid,               only : grid_t
    use tuvx_grid_warehouse,     only : grid_warehouse_t
    use tuvx_radiator_from_host, only : radiator_from_host_t
    use tuvx_radiator_warehouse, only : radiator_warehouse_t

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    type(grid_warehouse_t),      intent(in) :: grids     ! CAM grids used in TUV-x
    class(radiator_warehouse_t), pointer    :: radiators ! collection of radiators to be updated by CAM

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(radiator_from_host_t), pointer :: host_radiator
    class(grid_t), pointer :: height, wavelength

    radiators  => radiator_warehouse_t( )
    height     => grids%get_grid(     "height", "km" )
    wavelength => grids%get_grid( "wavelength", "nm" )

    ! Aerosol radiator
    host_radiator => radiator_from_host_t( "aerosol", height, wavelength )
    call radiators%add( host_radiator )
    deallocate( host_radiator )

    deallocate( height )
    deallocate( wavelength )

  end function get_cam_radiators

!================================================================================================

  subroutine create_updaters( this, grids, profiles, radiators )
!-----------------------------------------------------------------------
!
! Purpose: creates updaters for each grid and profile that CAM will use
!          to update TUV-x at each timestep
!
!-----------------------------------------------------------------------

    use musica_assert,           only : assert
    use tuvx_grid,               only : grid_t
    use tuvx_grid_from_host,     only : grid_updater_t
    use tuvx_grid_warehouse,     only : grid_warehouse_t
    use tuvx_profile,            only : profile_t
    use tuvx_profile_from_host,  only : profile_updater_t
    use tuvx_profile_warehouse,  only : profile_warehouse_t
    use tuvx_radiator,           only : radiator_t
    use tuvx_radiator_from_host, only : radiator_updater_t
    use tuvx_radiator_warehouse, only : radiator_warehouse_t

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr),             intent(inout) :: this
    class(grid_warehouse_t),     intent(in)    :: grids
    class(profile_warehouse_t),  intent(in)    :: profiles
    class(radiator_warehouse_t), intent(in)    :: radiators

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    class(grid_t),     pointer :: height, wavelength
    class(profile_t),  pointer :: host_profile
    class(radiator_t), pointer :: host_radiator
    logical                    :: found

    ! Grid updaters

    height => grids%get_grid( "height", "km" )
    this%grids_( GRID_INDEX_HEIGHT ) = this%core_%get_updater( height, found )
    call assert( 213798815, found )
    allocate( this%heights_(           height%size( ) + 1 ) ) ! TEMPORARY FOR DEVELOPMENT
    allocate( this%height_values_(     height%size( ) + 1 ) ) ! TEMPORARY FOR DEVELOPMENT
    allocate( this%height_mid_values_( height%size( )     ) ) ! TEMPORARY FOR DEVELOPMENT

    ! wavelength grid cannot be updated at runtime
    wavelength => grids%get_grid( "wavelength", "nm" )
    allocate( this%wavelength_values_( wavelength%size( ) + 1 ) ) ! TEMPORARY FOR DEVELOPMENT

    ! optical property working array
    allocate( this%optics_values_( height%size( ), wavelength%size( ) ) )

    deallocate( height )
    deallocate( wavelength )

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

    ! radiator updaters

    host_radiator => radiators%get_radiator( "aerosol" )
    this%radiators_( RADIATOR_INDEX_AEROSOL ) = this%core_%get_updater( host_radiator, found )
    call assert( 675200430, found )
    nullify( host_radiator )

  end subroutine create_updaters

!================================================================================================

  subroutine set_heights( this, i_col, ncol, height_mid, height_int)
!-----------------------------------------------------------------------
!
! Purpose: sets the height values in TUV-x for the given column
!
! NOTE: This function must be called before updating any profile data.
!
!  CAM to TUV-x height grid mapping
!
!  TUV-x heights are "bottom-up" and require atmospheric constituent
!  concentrations at interfaces. Therefore, CAM mid-points are used as
!  TUV-x grid interfaces, with an additional layer introduced between
!  the surface and the lowest CAM mid-point.
!
!  ---- (interface)  ===== (mid-point)
!
!        CAM                                  TUV-x
! ------(top)------ i_int = 1                              (exo values)
! ================= i_mid = 1           -------(top)------ i_int = pver + 1
! ----------------- i_int = 2           ================== i_mid = pver
!                                       ------------------ i_int = pver
!        ||
!        ||                                     ||
!                                               ||
! ----------------- i_int = pver
! ================= i_imd = pver        ------------------ i_int = 2
!                                       ================== i_mid = 1
! -----(ground)---- i_int = pver+1      -----(ground)----- i_int = 1
!
!-----------------------------------------------------------------------

    use ppgrid,         only : pcols ! maximum number of columns

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr), intent(inout) :: this                     ! TUV-x calculator
    integer,         intent(in)    :: i_col                    ! column to set conditions for
    integer,         intent(in)    :: ncol                     ! number of colums to calculated photolysis for
    real(r8),        intent(in)    :: height_mid(pcols,pver)   ! height above the surface at mid-points (km)
    real(r8),        intent(in)    :: height_int(pcols,pver+1) ! height above the surface at interfaces (km)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer :: i_level
    real(r8) :: edges(pver+1)
    real(r8) :: mid_points(pver)

    edges(1) = 0.0_r8
    edges(2:pver+1) = height_mid(i_col,pver:1:-1)
    mid_points(1) = height_mid(i_col,pver) * 0.5_r8
    mid_points(2:pver) = height_int(i_col,pver:2:-1)
    call this%grids_( GRID_INDEX_HEIGHT )%update( edges = edges, mid_points = mid_points )
    this%height_delta_(1:pver) = edges(2:pver+1) - edges(1:pver)

  end subroutine set_heights

!================================================================================================

  subroutine set_temperatures( this, i_col, temperature_mid, surface_temperature )
!-----------------------------------------------------------------------
!
! Purpose: sets the temperatures in TUV-x for the given column
!
! See description of `set_heights` for CAM <-> TUV-x vertical grid
! mapping.
!
!-----------------------------------------------------------------------

    use ppgrid,         only : pcols ! maximum number of columns

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr), intent(inout) :: this                        ! TUV-x calculator
    integer,         intent(in)    :: i_col                       ! column to set conditions for
    real(r8),        intent(in)    :: temperature_mid(pcols,pver) ! midpoint temperature (K)
    real(r8),        intent(in)    :: surface_temperature(pcols)  ! surface temperature (K)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    real(r8) :: edges(pver+1)

    edges(1) = surface_temperature(i_col)
    edges(2:pver+1) = temperature_mid(i_col,pver:1:-1)
    call this%profiles_( PROFILE_INDEX_TEMPERATURE )%update( edge_values = edges )

  end subroutine set_temperatures

!================================================================================================

  subroutine set_surface_albedo( this, i_col, surface_albedo )
!-----------------------------------------------------------------------
!
! Purpose: sets the surface albedo in TUV-x for the given column
!
! CAM uses a single value for surface albedo at all wavelengths
!
!-----------------------------------------------------------------------

    use ppgrid,         only : pcols ! maximum number of columns

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr), intent(inout) :: this                  ! TUV-x calculator
    integer,         intent(in)    :: i_col                 ! column to set conditions for
    real(r8),        intent(in)    :: surface_albedo(pcols) ! surface albedo (unitless)

    this%wavelength_values_(:) = surface_albedo(i_col)
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

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! Column to set conditions for

    ! TEMPORARY FOR DEVELOPMENT
    this%wavelength_values_(:) = 1000.0_r8
    call this%profiles_( PROFILE_INDEX_ET_FLUX )%update( &
        edge_values = this%wavelength_values_(:) )

  end subroutine set_et_flux

!================================================================================================

  subroutine set_radiator_profiles( this, i_col, ncol, fixed_species_conc, species_vmr )
!-----------------------------------------------------------------------
!
! Purpose: sets the profiles of optically active atmospheric constituents
!          in TUV-x for the given column
!
! TODO: Describe how CAM profiles are mapped to TUV-x heights
!
!-----------------------------------------------------------------------

    use chem_mods, only : gas_pcnst, & ! number of non-fixed species
                          nfs,       & ! number of fixed species
                          indexm       ! index for air density in fixed species array

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
    integer,         intent(in)    :: i_col ! column to set conditions for
    integer,         intent(in)    :: ncol  ! number of columns
    real(r8),        intent(in)    :: fixed_species_conc(ncol,pver,max(1,nfs)) ! fixed species densities
                                                                               !   (molecule cm-3)
    real(r8),        intent(in)    :: species_vmr(ncol,pver,max(1,gas_pcnst))  ! species volume mixing
                                                                               !   ratios (mol mol-1)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    real(r8) :: edges(pver+1), densities(pver)
    real(r8) :: km2cm = 1.0e5 ! conversion from km to cm

    ! air
    edges(1) = fixed_species_conc(i_col,pver,indexm)
    edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,indexm)
    densities(1:pver) = this%height_delta_(1:pver) * km2cm * &
                        sqrt(edges(1:pver)) + sqrt(edges(2:pver+1))
    call this%profiles_( PROFILE_INDEX_AIR )%update( &
        edge_values = edges, layer_densities = densities )

    ! O2
    if( is_fixed_O2 ) then
      edges(1) = fixed_species_conc(i_col,pver,index_O2)
      edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,index_O2)
    else if( index_O2 > 0 ) then
      edges(1) = species_vmr(i_col,pver,index_O2) * &
                 fixed_species_conc(i_col,pver,indexm)
      edges(2:pver+1) = species_vmr(i_col,pver:1:-1,index_O2) * &
                        fixed_species_conc(i_col,pver:1:-1,indexm)
    else
      edges(:) = 0.0_r8
    end if
    densities(1:pver) = this%height_delta_(1:pver) * km2cm * &
                        sqrt(edges(1:pver)) + sqrt(edges(2:pver+1))
    call this%profiles_( PROFILE_INDEX_O2 )%update( &
        edge_values = edges, layer_densities = densities )

    ! O3
    if( is_fixed_O3 ) then
      edges(1) = fixed_species_conc(i_col,pver,index_O3)
      edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,index_O3)
    else if( index_O3 > 0 ) then
      edges(1) = species_vmr(i_col,pver,index_O3) * &
                 fixed_species_conc(i_col,pver,indexm)
      edges(2:pver+1) = species_vmr(i_col,pver:1:-1,index_O3) * &
                        fixed_species_conc(i_col,pver:1:-1,indexm)
    else
      edges(:) = 0.0_r8
    end if
    densities(1:pver) = this%height_delta_(1:pver) * km2cm * &
                        sqrt(edges(1:pver)) + sqrt(edges(2:pver+1))
    call this%profiles_( PROFILE_INDEX_O3 )%update( &
        edge_values = edges, layer_densities = densities )

    ! TEMPORARY FOR DEVELOPMENT - aerosols
    this%optics_values_(:,:) = 0.01_r8
    call this%radiators_( RADIATOR_INDEX_AEROSOL )%update( &
        optical_depths = this%optics_values_ )

  end subroutine set_radiator_profiles

!================================================================================================

end module mo_tuvx
