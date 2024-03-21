!----------------------------------------------------------------------
! Wrapper for TUV-x photolysis rate constant calculator
!----------------------------------------------------------------------
module mo_tuvx

   use musica_map,              only : map_t
   use musica_string,           only : string_t
   use ppgrid,                  only : pver,pverp                ! number of vertical layers
   use shr_kind_mod,            only : r8 => shr_kind_r8, cl=>shr_kind_cl
   use tuvx_core,               only : core_t
   use tuvx_grid_from_host,     only : grid_updater_t
   use tuvx_profile_from_host,  only : profile_updater_t
   use tuvx_radiator_from_host, only : radiator_updater_t

   implicit none

   private

   public :: tuvx_readnl
   public :: tuvx_register
   public :: tuvx_init
   public :: tuvx_timestep_init
   public :: tuvx_get_photo_rates
   public :: tuvx_finalize
   public :: tuvx_active

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
   integer, parameter :: NUM_RADIATORS = 2          ! number of radiators that CAM will update at runtime
   integer, parameter :: RADIATOR_INDEX_AEROSOL = 1 ! Aerosol radiator index
   integer, parameter :: RADIATOR_INDEX_CLOUDS  = 2 ! Cloud radiator index

   ! Definition of the MS93 wavelength grid  TODO add description of this
   integer,       parameter :: NUM_BINS_MS93 = 4
   real(kind=r8), parameter :: WAVELENGTH_EDGES_MS93(NUM_BINS_MS93+1) = &
      (/ 181.6_r8, 183.1_r8, 184.6_r8, 190.2_r8, 192.5_r8 /)

   ! Heating rate indices
   integer :: number_of_heating_rates = 0 ! number of heating rates in TUV-x
   integer :: index_cpe_jo2_a = -1 ! index for jo2_a in heating rate array
   integer :: index_cpe_jo2_b = -1 ! index for jo2_b in heating rate array
   integer :: index_cpe_jo3_a = -1 ! index for jo3_a in heating rate array
   integer :: index_cpe_jo3_b = -1 ! index for jo3_b in heating rate array
   integer :: cpe_jo2_a_pbuf_index = -1 ! index in physics buffer for jo2_a heating rate
   integer :: cpe_jo2_b_pbuf_index = -1 ! index in physics buffer for jo2_b heating rate
   integer :: cpe_jo3_a_pbuf_index = -1 ! index in physics buffer for jo3_a heating rate
   integer :: cpe_jo3_b_pbuf_index = -1 ! index in physics buffer for jo3_b heating rate

   ! Information needed to access CAM species state data
   logical :: is_fixed_N2 = .false. ! indicates whether N2 concentrations are fixed
   logical :: is_fixed_O  = .false. ! indicates whether O concentrations are fixed
   logical :: is_fixed_O2 = .false. ! indicates whether O2 concentrations are fixed
   logical :: is_fixed_O3 = .false. ! indicates whether O3 concentrations are fixed
   logical :: is_fixed_NO = .false. ! indicates whether NO concentrations are fixed
   integer :: index_N2 = 0 ! index for N2 in concentration array
   integer :: index_O  = 0 ! index for O in concentration array
   integer :: index_O2 = 0 ! index for O2 in concentration array
   integer :: index_O3 = 0 ! index for O3 in concentration array
   integer :: index_NO = 0 ! index for NO in concentration array

   ! Information needed to access aerosol and cloud optical properties
   logical :: do_aerosol = .false. ! indicates whether aerosol optical properties
   !   are available and should be used in radiative
   !   transfer calculations
   logical :: do_clouds  = .false. ! indicates whether cloud optical properties
   !   should be calculated and used in radiative
   !   transfer calculations

   ! Information needed to set extended-UV photo rates
   logical :: do_euv = .false.              ! Indicates whether to calculate
   !   extended-UV photo rates
   integer :: ion_rates_pbuf_index = 0      ! Index in physics buffer for
   !   ionization rates

   ! Information needed to do special NO photolysis rate calculation
   logical :: do_jno     = .false. ! Indicates whether to calculate jno
   integer :: jno_index  = 0       ! Index in tuvx_ptr::photo_rates_ array for jno

   ! Cutoff solar zenith angle for doing photolysis rate calculations [degrees]
   real(r8) :: max_sza = 0.0_r8

   ! TODO how should these paths be set and communicated to this wrapper?
   character(len=*), parameter :: wavelength_config_path = &
      "data/grids/wavelength/cam.csv"
   logical, parameter :: enable_diagnostics = .true.

   ! TUV-x calculator for each OMP thread
   type :: tuvx_ptr
      type(core_t), pointer    :: core_ => null( )          ! TUV-x calculator
      integer                  :: n_photo_rates_ = 0        ! number of photo reactions in TUV-x
      integer                  :: n_euv_rates_ = 0          ! number of extreme-UV rates
      integer                  :: n_special_rates_ = 0      ! number of special photo rates
      real(r8), allocatable    :: photo_rates_(:,:,:)       ! photolysis rate constants
      !   (column, vertical level, reaction) [s-1]
      type(map_t)              :: photo_rate_map_           ! map between TUV-x and CAM
      !   photo rate constant arrays
      type(grid_updater_t)     :: grids_(NUM_GRIDS)         ! grid updaters
      type(profile_updater_t)  :: profiles_(NUM_PROFILES)   ! profile updaters
      type(radiator_updater_t) :: radiators_(NUM_RADIATORS) ! radiator updaters
      real(r8)                 :: height_delta_(pver+1)     ! change in height in each
      !   vertical layer (km)
      real(r8), allocatable    :: wavelength_edges_(:)      ! TUV-x wavelength bin edges (nm)
      real(r8), allocatable    :: wavelength_values_(:)     ! Working array for interface values
      !   on the TUV-x wavelength grid
      real(r8), allocatable    :: wavelength_mid_values_(:) ! Working array for mid-point values
      !   on the TUV-x wavelength grid
      real(r8), allocatable    :: optical_depth_(:,:,:)            ! (column, vertical level, wavelength) [unitless]
      real(r8), allocatable    :: single_scattering_albedo_(:,:,:) ! (column, vertical level, wavelength) [unitless]
      real(r8), allocatable    :: asymmetry_factor_(:,:,:)         ! (column, vertical level, wavelength) [unitless]
      real(r8)                 :: et_flux_ms93_(NUM_BINS_MS93)     ! extraterrestrial flux on the MS93 grid
      !   [photon cm-2 nm-1 s-1]
   end type tuvx_ptr
   type(tuvx_ptr), allocatable :: tuvx_ptrs(:)

   ! Diagnostic photolysis rate constant output
   type :: diagnostic_t
      character(len=:), allocatable :: name_ ! Name of the output field
      integer :: index_                      ! index of the photolysis rate constant from TUV-x
   end type diagnostic_t
   type(diagnostic_t), allocatable :: diagnostics(:)

   ! namelist options
   character(len=cl) :: tuvx_config_path = 'NONE'  ! absolute path to TUVX configuration file
   logical, protected :: tuvx_active = .false.

!================================================================================================
contains
!================================================================================================

   !-----------------------------------------------------------------------
   ! registers fields in the physics buffer
   !-----------------------------------------------------------------------
   subroutine tuvx_register( )

      use mo_jeuv,        only : nIonRates
      use physics_buffer, only : pbuf_add_field, dtype_r8
      use ppgrid,         only : pcols ! maximum number of columns

      if( .not. tuvx_active ) return

      ! add photo-ionization rates to physics buffer for WACCMX Ionosphere module
      call pbuf_add_field( 'IonRates', 'physpkg', dtype_r8, (/ pcols, pver, nIonRates /), &
         ion_rates_pbuf_index ) ! Ionization rates for O+, O2+, N+, N2+, NO+

      call pbuf_add_field( 'CPE_jO2a', 'global', dtype_r8, (/ pcols, pver /), cpe_jo2_a_pbuf_index )
      call pbuf_add_field( 'CPE_jO2b', 'global', dtype_r8, (/ pcols, pver /), cpe_jo2_b_pbuf_index )
      call pbuf_add_field( 'CPE_jO3a', 'global', dtype_r8, (/ pcols, pver /), cpe_jo3_a_pbuf_index )
      call pbuf_add_field( 'CPE_jO3b', 'global', dtype_r8, (/ pcols, pver /), cpe_jo3_b_pbuf_index )

   end subroutine tuvx_register

!================================================================================================

   !-----------------------------------------------------------------------
   ! read namelist options
   !-----------------------------------------------------------------------
   subroutine tuvx_readnl(nlfile)

#ifdef HAVE_MPI
      use mpi
#endif
      use cam_abortutils, only : endrun
      use cam_logfile,    only : iulog ! log file output unit
      use namelist_utils, only : find_group_name
      use spmd_utils,     only : mpicom, is_main_task => masterproc, &
         main_task => masterprocid

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
      if (is_main_task) then
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
#ifdef HAVE_MPI
      call mpi_bcast(tuvx_config_path, len(tuvx_config_path), mpi_character, main_task, mpicom, ierr)
      call mpi_bcast(tuvx_active,      1,                     mpi_logical,   main_task, mpicom, ierr)
#endif

      if (tuvx_active .and. tuvx_config_path == 'NONE') then
         call endrun(subname // ' : must set tuvx_config_path when TUV-X is active')
      end if

      if (is_main_task) then
         write(iulog,*) 'tuvx_readnl: tuvx_config_path = ', trim(tuvx_config_path)
         write(iulog,*) 'tuvx_readnl: tuvx_active = ', tuvx_active
      end if

   end subroutine tuvx_readnl

!================================================================================================

   !-----------------------------------------------------------------------
   ! Initializes TUV-x for photolysis calculations
   !-----------------------------------------------------------------------
   subroutine tuvx_init( photon_file, electron_file, max_solar_zenith_angle, pbuf2d )

#ifdef HAVE_MPI
      use mpi
#endif
      use cam_history,             only : addfld
      use cam_logfile,             only : iulog ! log file output unit
      use infnan,                  only : nan, assignment(=)
      use mo_chem_utls,            only : get_spc_ndx, get_inv_ndx
      use mo_jeuv,                 only : neuv ! number of extreme-UV rates
      use musica_assert,           only : assert_msg, die_msg
      use musica_config,           only : config_t
      use musica_mpi,              only : musica_mpi_rank, &
         musica_mpi_pack_size, &
         musica_mpi_pack, &
         musica_mpi_unpack
      use musica_string,           only : string_t, to_char
      use physics_buffer,          only : physics_buffer_desc
      use physics_buffer,          only : pbuf_set_field
      use ppgrid,                  only : pcols ! maximum number of columns
      use shr_const_mod,           only : pi => shr_const_pi
      use solar_irrad_data,        only : has_spectrum
      use spmd_utils,              only : main_task => masterprocid, &
         is_main_task => masterproc, &
         mpicom
      use tuvx_grid,               only : grid_t
      use tuvx_grid_warehouse,     only : grid_warehouse_t
      use tuvx_profile_warehouse,  only : profile_warehouse_t
      use tuvx_radiator_warehouse, only : radiator_warehouse_t
      use time_manager,            only : is_first_step

      character(len=*), intent(in) :: photon_file   ! photon file used in extended-UV module setup
      character(len=*), intent(in) :: electron_file ! electron file used in extended-UV module setup
      real(r8),         intent(in) :: max_solar_zenith_angle ! cutoff solar zenith angle for
      !    photo rate calculations [degrees]
      type(physics_buffer_desc), pointer :: pbuf2d(:,:) ! Physics buffer

      character(len=*), parameter :: my_name = "TUV-x wrapper initialization"
      class(core_t), pointer :: core
      character, allocatable :: buffer(:)
      type(string_t) :: config_path
      type(config_t) :: tuvx_config, cam_config, map_config
      type(map_t) :: map
      class(grid_t), pointer :: height
      class(grid_t), pointer :: wavelength
      class(grid_warehouse_t), pointer :: cam_grids
      class(profile_warehouse_t), pointer :: cam_profiles
      class(radiator_warehouse_t), pointer :: cam_radiators
      integer :: pack_size, pos, i_core, i_err
      logical :: disable_aerosols, disable_clouds
      type(string_t) :: required_keys(1), optional_keys(2)
      logical, save :: is_initialized = .false.

      type(string_t), allocatable :: labels(:)
      character(len=16) :: label
      integer :: i
      real(r8) :: nanval

      if( .not. tuvx_active ) return
      if( is_initialized ) return
      is_initialized = .true.

      nanval=nan

      call pbuf_set_field( pbuf2d, ion_rates_pbuf_index, nanval )

      if( is_first_step( ) ) then
        call pbuf_set_field( pbuf2d, cpe_jo2_a_pbuf_index, 0.0_r8 )
        call pbuf_set_field( pbuf2d, cpe_jo2_b_pbuf_index, 0.0_r8 )
        call pbuf_set_field( pbuf2d, cpe_jo3_a_pbuf_index, 0.0_r8 )
        call pbuf_set_field( pbuf2d, cpe_jo3_b_pbuf_index, 0.0_r8 )
      end if

      if( is_main_task ) write(iulog,*) "Beginning TUV-x Initialization"

      config_path = trim(tuvx_config_path)

      ! ===============================
      ! CAM TUV-x configuration options
      ! ===============================
      required_keys(1) = "aliasing"
      optional_keys(1) = "disable aerosols"
      optional_keys(2) = "disable clouds"

#ifndef HAVE_MPI
      call assert_msg( 113937299, is_main_task, "Multiple tasks present without " &
         //"MPI support enabled for TUV-x" )
#endif

      ! ===============================================================
      ! set the maximum solar zenith angle to calculate photo rates for
      ! ===============================================================
      max_sza = max_solar_zenith_angle
      if( max_sza <= 0.0_r8 .or. max_sza > 180.0_r8 ) then
         call die_msg( 723815691, "TUV-x max solar zenith angle must be between 0 and 180 degress" )
      end if

      ! =================================
      ! initialize the extended-UV module
      ! =================================
      call initialize_euv( photon_file, electron_file, do_euv )

      ! ==========================================================================
      ! create the set of TUV-x grids and profiles that CAM will update at runtime
      ! ==========================================================================
      cam_grids => get_cam_grids( wavelength_config_path )
      cam_profiles => get_cam_profiles( cam_grids )
      cam_radiators => get_cam_radiators( cam_grids )

      ! ==================================================================
      ! construct a core and a map between TUV-x and CAM photolysis arrays
      ! on the primary process and pack them onto an MPI buffer
      ! ==================================================================
      if( is_main_task ) then
         call tuvx_config%from_file( config_path%to_char( ) )
         call tuvx_config%get( "__CAM options", cam_config, my_name )
         call assert_msg( 973680295, &
            cam_config%validate( required_keys, optional_keys ), &
            "Bad configuration for CAM TUV-x options." )
         call cam_config%get( "disable aerosols", disable_aerosols, my_name, &
            default = .false. )
         call cam_config%get( "disable clouds", disable_clouds, my_name, &
            default = .false. )
         call cam_config%get( "aliasing", map_config, my_name )
         core => core_t( config_path, cam_grids, cam_profiles, cam_radiators )
         call set_photo_rate_map( core, map_config, do_euv, do_jno, jno_index, map )
         pack_size = core%pack_size( mpicom ) + &
            map%pack_size( mpicom ) + &
            musica_mpi_pack_size( do_jno, mpicom ) + &
            musica_mpi_pack_size( jno_index, mpicom ) + &
            musica_mpi_pack_size( disable_aerosols, mpicom ) + &
            musica_mpi_pack_size( disable_clouds, mpicom )
         allocate( buffer( pack_size ) )
         pos = 0
         call core%mpi_pack( buffer, pos, mpicom )
         call map%mpi_pack(  buffer, pos, mpicom )
         call musica_mpi_pack( buffer, pos, do_jno,           mpicom )
         call musica_mpi_pack( buffer, pos, jno_index,        mpicom )
         call musica_mpi_pack( buffer, pos, disable_aerosols, mpicom )
         call musica_mpi_pack( buffer, pos, disable_clouds,   mpicom )
         deallocate( core )
      end if

#ifdef HAVE_MPI
      ! ====================================================
      ! broadcast the core and map data to all MPI processes
      ! ====================================================
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

      ! ================================================================
      ! unpack the core and map for each OMP thread on every MPI process
      ! ================================================================
      allocate( tuvx_ptrs( max_threads( ) ) )
      do i_core = 1, size( tuvx_ptrs )
         associate( tuvx => tuvx_ptrs( i_core ) )
            allocate( tuvx%core_ )
            pos = 0
            call tuvx%core_%mpi_unpack( buffer, pos, mpicom )
            call tuvx%photo_rate_map_%mpi_unpack( buffer, pos, mpicom )
            call musica_mpi_unpack( buffer, pos, do_jno,           mpicom )
            call musica_mpi_unpack( buffer, pos, jno_index,        mpicom )
            call musica_mpi_unpack( buffer, pos, disable_aerosols, mpicom )
            call musica_mpi_unpack( buffer, pos, disable_clouds,   mpicom )

            ! ===================================================================
            ! Set up connections between CAM and TUV-x input data for each thread
            ! ===================================================================
            call create_updaters( tuvx, cam_grids, cam_profiles, cam_radiators, &
               disable_aerosols, disable_clouds )

            ! ===============================================================
            ! Create a working array for calculated photolysis rate constants
            ! ===============================================================
            tuvx%n_photo_rates_ = tuvx%core_%number_of_photolysis_reactions( )
            if( do_euv ) tuvx%n_euv_rates_ = neuv
            if( do_jno ) tuvx%n_special_rates_ = tuvx%n_special_rates_ + 1
            height => tuvx%core_%get_grid( "height", "km" )
            allocate( tuvx%photo_rates_( pcols, height%ncells_ + 1, &
               tuvx%n_photo_rates_ + tuvx%n_euv_rates_ + &
               tuvx%n_special_rates_ ) )
            deallocate( height )

         end associate
      end do
      deallocate( cam_grids     )
      deallocate( cam_profiles  )
      deallocate( cam_radiators )

      ! =============================================
      ! Get index info for CAM species concentrations
      ! =============================================
      index_N2  = get_inv_ndx( 'N2' )
      is_fixed_N2 = index_N2 > 0
      if( .not. is_fixed_N2 ) index_N2 = get_spc_ndx( 'N2' )
      index_O  = get_inv_ndx( 'O' )
      is_fixed_O = index_O > 0
      if( .not. is_fixed_O ) index_O = get_spc_ndx( 'O' )
      index_O2 = get_inv_ndx( 'O2' )
      is_fixed_O2 = index_O2 > 0
      if( .not. is_fixed_O2 ) index_O2 = get_spc_ndx( 'O2' )
      index_O3 = get_inv_ndx( 'O3' )
      is_fixed_O3 = index_O3 > 0
      if( .not. is_fixed_O3 ) index_O3 = get_spc_ndx( 'O3' )
      index_NO = get_inv_ndx( 'NO' )
      is_fixed_NO = index_NO > 0
      if( .not. is_fixed_NO ) index_NO = get_spc_ndx( 'NO' )

      ! ====================================================
      ! make sure extraterrestrial flux values are available
      ! ====================================================
      call assert_msg( 170693514, has_spectrum, &
         "Solar irradiance spectrum needed for TUV-x" )

      ! ============================================
      ! set up diagnostic output of photolysis rates
      ! ============================================
      call initialize_diagnostics( tuvx_ptrs( 1 ) )

      ! ===============================
      ! set up map to CAM heating rates
      ! ===============================
      labels = tuvx_ptrs(1)%core_%heating_rate_labels( )
      do i = 1, size( labels )
         label = trim( labels( i )%to_char( ) )
         select case( label )
         case( 'jo2_a' )
            index_cpe_jo2_a = i
            number_of_heating_rates = number_of_heating_rates + 1
            call addfld('CPE_jO2a',(/ 'lev' /), 'A', 'joules sec-1', &
                        trim(label)//' chemical potential energy')
         case( 'jo2_b' )
            index_cpe_jo2_b = i
            number_of_heating_rates = number_of_heating_rates + 1
            call addfld('CPE_jO2b',(/ 'lev' /), 'A', 'joules sec-1', &
                        trim(label)//' chemical potential energy')
         case( 'jo3_a' )
            index_cpe_jo3_a = i
            number_of_heating_rates = number_of_heating_rates + 1
            call addfld('CPE_jO3a',(/ 'lev' /), 'A', 'joules sec-1', &
                        trim(label)//' chemical potential energy')
         case( 'jo3_b' )
            index_cpe_jo3_b = i
            number_of_heating_rates = number_of_heating_rates + 1
            call addfld('CPE_jO3b',(/ 'lev' /), 'A', 'joules sec-1', &
                        trim(label)//' chemical potential energy')
         end select
      end do
      call assert_msg( 398372957, &
                       number_of_heating_rates == size( labels ), &
                       "TUV-x heating rate mismatch. Expected "// &
                       trim( to_char( size( labels ) )// &
                       " rates, but only matched "// &
                       trim( to_char( number_of_heating_rates ) )//"." ) )

      if( is_main_task ) call log_initialization( labels )

   end subroutine tuvx_init

!================================================================================================

   !-----------------------------------------------------------------------
   ! Updates TUV-x profiles that depend on time but not space
   !-----------------------------------------------------------------------
   subroutine tuvx_timestep_init( )

      integer :: i_thread

      if( .not. tuvx_active ) return

      do i_thread = 1, size( tuvx_ptrs )
         associate( tuvx => tuvx_ptrs( i_thread ) )
            call set_et_flux( tuvx )
         end associate
      end do

   end subroutine tuvx_timestep_init

!================================================================================================

   !-----------------------------------------------------------------------
   ! Calculates and returns photolysis rate constants
   !-----------------------------------------------------------------------
   subroutine tuvx_get_photo_rates( state, pbuf, ncol, lchnk, height_mid, &
      height_int, temperature_mid, surface_temperature, fixed_species_conc, &
      species_vmr, exo_column_conc, surface_albedo, solar_zenith_angle, &
      earth_sun_distance, pressure_delta, cloud_fraction, liquid_water_content, &
      photolysis_rates )

      use cam_history,      only : outfld
      use cam_logfile,      only : iulog        ! log info output unit
      use chem_mods,        only : phtcnt,    & ! number of photolysis reactions
         gas_pcnst, & ! number of non-fixed species
         nfs,       & ! number of fixed species
         nabscol      ! number of absorbing species (radiators)
      use physics_types,    only : physics_state
      use physics_buffer,   only : physics_buffer_desc
      use physics_buffer,   only : pbuf_get_field
      use ppgrid,           only : pcols        ! maximum number of columns
      use shr_const_mod,    only : pi => shr_const_pi
      use spmd_utils,       only : main_task => masterprocid, &
         is_main_task => masterproc, &
         mpicom

      type(physics_state),       target,  intent(in)    :: state
      type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
      integer,  intent(in)    :: ncol                        ! number of active columns on this thread
      integer,  intent(in)    :: lchnk                       ! identifier for this thread
      real(r8), intent(in)    :: height_mid(ncol,pver)       ! height at mid-points (km)
      real(r8), intent(in)    :: height_int(ncol,pver+1)     ! height at interfaces (km)
      real(r8), intent(in)    :: temperature_mid(pcols,pver) ! midpoint temperature (K)
      real(r8), intent(in)    :: surface_temperature(pcols)  ! surface temperature (K)
      real(r8), intent(in)    :: fixed_species_conc(ncol,pver,max(1,nfs))    ! fixed species densities
      !   (molecule cm-3)
      real(r8), intent(in)    :: species_vmr(ncol,pver,max(1,gas_pcnst))     ! species volume mixing
      !   ratios (mol mol-1)
      real(r8), intent(in)    :: exo_column_conc(ncol,0:pver,max(1,nabscol)) ! layer column densities
      !   (molecule cm-2)
      real(r8), intent(in)    :: surface_albedo(pcols)       ! surface albedo (unitless)
      real(r8), intent(in)    :: solar_zenith_angle(ncol)    ! solar zenith angle (radians)
      real(r8), intent(in)    :: earth_sun_distance          ! Earth-Sun distance (AU)
      real(r8), intent(in)    :: pressure_delta(pcols,pver)  ! pressure delta about midpoints (Pa)
      real(r8), intent(in)    :: cloud_fraction(ncol,pver)   ! cloud fraction (unitless)
      real(r8), intent(in)    :: liquid_water_content(ncol,pver)    ! liquid water content (kg/kg)
      real(r8), intent(inout) :: photolysis_rates(ncol,pver,phtcnt) ! photolysis rate
      !   constants (1/s)

      integer  :: i_col   ! column index
      integer  :: i_level ! vertical level index
      real(r8) :: sza     ! solar zenith angle [degrees]
      real(r8) :: cpe_rates(ncol,pverp+1,number_of_heating_rates) ! heating rates from TUV-x
      real(r8), pointer :: cpe_jo2_a(:,:) ! heating rate for jo2_a in physics buffer
      real(r8), pointer :: cpe_jo2_b(:,:) ! heating rate for jo2_b in physics buffer
      real(r8), pointer :: cpe_jo3_a(:,:) ! heating rate for jo3_a in physics buffer
      real(r8), pointer :: cpe_jo3_b(:,:) ! heating rate for jo3_b in physics buffer

      if( .not. tuvx_active ) return

      call pbuf_get_field(pbuf, cpe_jo2_a_pbuf_index, cpe_jo2_a)
      call pbuf_get_field(pbuf, cpe_jo2_b_pbuf_index, cpe_jo2_b)
      call pbuf_get_field(pbuf, cpe_jo3_a_pbuf_index, cpe_jo3_a)
      call pbuf_get_field(pbuf, cpe_jo3_b_pbuf_index, cpe_jo3_b)

      cpe_rates(:,:,:) = 0.0_r8
      cpe_jo2_a(:,:) = 0.0_r8
      cpe_jo2_b(:,:) = 0.0_r8
      cpe_jo3_a(:,:) = 0.0_r8
      cpe_jo3_b(:,:) = 0.0_r8

      associate( tuvx => tuvx_ptrs( thread_id( ) ) )

         tuvx%photo_rates_(:,:,:) = 0.0_r8

         ! ==============================================
         ! set aerosol optical properties for all columns
         ! ==============================================
         call get_aerosol_optical_properties( tuvx, state, pbuf )

         do i_col = 1, ncol

            ! ===================================
            ! skip columns in near total darkness
            ! ===================================
            sza = solar_zenith_angle(i_col) * 180.0_r8 / pi
            if( sza < 0.0_r8 .or. sza > max_sza ) cycle

            ! ===================
            ! update grid heights
            ! ===================
            call set_heights( tuvx, i_col, ncol, height_mid, height_int )

            ! =======================================
            ! set conditions for this column in TUV-x
            ! =======================================
            call set_temperatures( tuvx, i_col, temperature_mid, surface_temperature )
            call set_surface_albedo( tuvx, i_col, surface_albedo )
            call set_radiator_profiles( tuvx, i_col, ncol, fixed_species_conc, &
               species_vmr, exo_column_conc, &
               pressure_delta(1:ncol,:), cloud_fraction, &
               liquid_water_content )

            ! ===================================================
            ! Calculate photolysis rate constants for this column
            ! ===================================================
            call tuvx%core_%run( solar_zenith_angle = sza, &
               earth_sun_distance = earth_sun_distance, &
               photolysis_rate_constants = &
               tuvx%photo_rates_(i_col,:,1:tuvx%n_photo_rates_), &
               heating_rates = cpe_rates(i_col,:,:) )

            ! ==============================
            ! Calculate the extreme-UV rates
            ! ==============================
            if( do_euv ) then
               associate( euv_begin => tuvx%n_photo_rates_ + 1, &
                  euv_end   => tuvx%n_photo_rates_ + tuvx%n_euv_rates_ )
                  call calculate_euv_rates( sza, &
                     fixed_species_conc(i_col,:,:), &
                     species_vmr(i_col,:,:), &
                     height_mid(i_col,:), &
                     height_int(i_col,:), &
                     tuvx%photo_rates_(i_col,2:pver+1,euv_begin:euv_end) )
               end associate
            end if

            ! =============================
            ! Calculate special photo rates
            ! =============================
            if( do_jno ) then
               call calculate_jno( sza, &
                  tuvx%et_flux_ms93_, &
                  fixed_species_conc(i_col,:,:), &
                  species_vmr(i_col,:,:), &
                  height_int(i_col,:), &
                  tuvx%photo_rates_(i_col,2:pver+1,jno_index) )
            end if
         end do

         ! =====================
         ! Filter negative rates
         ! =====================
         tuvx%photo_rates_(:,:,:) = max( 0.0_r8, tuvx%photo_rates_(:,:,:) )

         ! ============================================
         ! Return the photolysis rates on the CAM grids
         ! ============================================
         do i_col = 1, ncol
            do i_level = 1, pver
               call tuvx%photo_rate_map_%apply( tuvx%photo_rates_(i_col,pver-i_level+2,:), &
                  photolysis_rates(i_col,i_level,:) )
            end do
         end do

         call output_diagnostics( tuvx, ncol, lchnk )

      end associate

      if (index_cpe_jo2_a>0) then
         do i_level = 1, pver
            cpe_jo2_a(:ncol,i_level) = &
                0.5_r8 * ( cpe_rates(:ncol,pver-i_level+2,index_cpe_jo2_a) &
                           + cpe_rates(:ncol,pver-i_level+1,index_cpe_jo2_a) )
         end do
         call outfld('CPE_jO2a', cpe_jo2_a(:ncol,:), ncol, lchnk )
      end if
      if (index_cpe_jo2_b>0) then
         do i_level = 1, pver
            cpe_jo2_b(:ncol,i_level) = &
                0.5_r8 * ( cpe_rates(:ncol,pver-i_level+2,index_cpe_jo2_b) &
                           + cpe_rates(:ncol,pver-i_level+1,index_cpe_jo2_b) )
         end do
         call outfld('CPE_jO2b', cpe_jo2_b(:ncol,:), ncol, lchnk )
      end if

      if (index_cpe_jo3_a>0) then
         do i_level = 1, pver
            cpe_jo3_a(:ncol,i_level) = &
                0.5_r8 * ( cpe_rates(:ncol,pver-i_level+2,index_cpe_jo3_a) &
                           + cpe_rates(:ncol,pver-i_level+1,index_cpe_jo3_a) )
         end do
         call outfld('CPE_jO3a', cpe_jo3_a(:ncol,:), ncol, lchnk )
      end if
      if (index_cpe_jo3_b>0) then
         do i_level = 1, pver
            cpe_jo3_b(:ncol,i_level) = &
                0.5_r8 * ( cpe_rates(:ncol,pver-i_level+2,index_cpe_jo3_b) &
                           + cpe_rates(:ncol,pver-i_level+1,index_cpe_jo3_b) )
         end do
         call outfld('CPE_jO3b', cpe_jo3_b(:ncol,:), ncol, lchnk )
      end if

   end subroutine tuvx_get_photo_rates

!================================================================================================

   !-----------------------------------------------------------------------
   ! Cleans up memory associated with TUV-x calculators
   !-----------------------------------------------------------------------
   subroutine tuvx_finalize( )

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

   !-----------------------------------------------------------------------
   ! Returns the id of the current OpenMP thread, or 1 if not
   !   using OpenMP (1 <= id <= max_threads())
   !-----------------------------------------------------------------------
   integer function thread_id( )
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

   !-----------------------------------------------------------------------
   ! Returns the number of threads available for calculations at
   !   runtime
   !-----------------------------------------------------------------------
   integer function max_threads( )
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

   !-----------------------------------------------------------------------
   ! Prints initialization conditions to the log file
   !-----------------------------------------------------------------------
   subroutine log_initialization( heating_rate_labels )

      use cam_logfile,    only : iulog ! log info output unit
      use musica_string,  only : to_char
      use spmd_utils,     only : main_task => masterprocid, &
         is_main_task => masterproc

      type(string_t), intent(in) :: heating_rate_labels(:) ! heating rate labels

      integer :: i

      if( is_main_task ) then
         write(iulog,*) "Initialized TUV-x"
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
         write(iulog,*) "  - with configuration file: '"//trim( tuvx_config_path )//"'"
         if( do_aerosol ) then
            write(iulog,*) "  - with on-line aerosols"
         else
            write(iulog,*) "  - without on-line aerosols"
         end if
         if( do_clouds ) then
            write(iulog,*) "  - with on-line clouds"
         else
            write(iulog,*) "  - without on-line clouds"
         end if
         if( index_N2 > 0 ) write(iulog,*) "  - including N2"
         if( index_O  > 0 ) write(iulog,*) "  - including O"
         if( index_O2 > 0 ) write(iulog,*) "  - including O2"
         if( index_O3 > 0 ) write(iulog,*) "  - including O3"
         if( index_NO > 0 ) write(iulog,*) "  - including NO"
         if( do_euv ) write(iulog,*) "  - doing Extreme-UV calculations"
         if( do_jno ) write(iulog,*) "  - including special jno rate calculation"
         write(iulog,*) "  - max solar zenith angle [degrees]:", max_sza
         if( size( heating_rate_labels ) > 0 ) then
            write(iulog,*) "  - with heating rates:"
            do i = 1, size( heating_rate_labels )
               write(iulog,*) "    - "//trim( heating_rate_labels( i )%to_char( ) )
            end do
         end if
      end if

   end subroutine log_initialization

!================================================================================================

   !-----------------------------------------------------------------------
   ! initializes the external extreme-UV module
   !-----------------------------------------------------------------------
   subroutine initialize_euv( photon_file, electron_file, do_euv )

      use chem_mods,    only : phtcnt    ! number of CAM-Chem photolysis reactions
      use mo_jeuv,      only : jeuv_init ! extreme-UV initialization

      character(len=*),     intent(in)  :: photon_file   ! photon file used in extended-UV module
      !   setup
      character(len=*),     intent(in)  :: electron_file ! electron file used in extended-UV
      !   module setup
      logical,              intent(out) :: do_euv        ! indicates whether extreme-UV
      !   calculations are needed

      integer, allocatable :: euv_index_map(:)

      allocate( euv_index_map( phtcnt ) )
      euv_index_map(:) = 0
      call jeuv_init( photon_file, electron_file, euv_index_map )
      do_euv = any( euv_index_map(:) > 0 )

   end subroutine initialize_euv

!================================================================================================

   !-----------------------------------------------------------------------
   ! Registers fields for diagnostic output
   !-----------------------------------------------------------------------
   subroutine initialize_diagnostics( this )

      use cam_history,   only : addfld
      use musica_assert, only : assert
      use musica_string, only : string_t

      type(tuvx_ptr), intent(in) :: this

      type(string_t), allocatable :: labels(:), all_labels(:)
      integer :: i_label

      if( .not. enable_diagnostics ) then
         allocate( diagnostics( 0 ) )
         return
      end if

      ! ==========================================================
      ! add output for specific photolysis reaction rate constants
      ! ==========================================================
      labels = this%core_%photolysis_reaction_labels( )
      allocate( all_labels( size( labels ) + this%n_special_rates_ ) )
      all_labels( 1 : size( labels ) ) = labels(:)
      i_label = size( labels ) + 1
      if( do_jno ) then
         all_labels( i_label ) = "jno"
         i_label = i_label + 1
      end if
      call assert( 522515214, i_label == size( all_labels ) + 1 )
      allocate( diagnostics( size( all_labels ) ) )
      do i_label = 1, size( all_labels )
         diagnostics( i_label )%name_  = trim( all_labels( i_label )%to_char( ) )
         diagnostics( i_label )%index_ = i_label
         call addfld( "tuvx_"//diagnostics( i_label )%name_, (/ 'lev' /), 'A', 'sec-1', &
            'photolysis rate constant' )
      end do

   end subroutine initialize_diagnostics

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets up a map between the TUV-x and CAM photolysis rate arrays
   !-----------------------------------------------------------------------
   subroutine set_photo_rate_map( core, config, do_euv, do_jno, jno_index, map )

      use cam_logfile,   only : iulog       ! log info output unit
      use chem_mods,     only : phtcnt, &   ! number of photolysis reactions
         rxt_tag_lst ! labels for all chemical reactions
      ! NOTE photolysis reactions are
      ! expected to appear first
      use mo_jeuv,       only : neuv        ! number of extreme-UV rates
      use musica_config, only : config_t
      use musica_string, only : string_t

      type(core_t),   intent(in)    :: core      ! TUV-x core
      type(config_t), intent(inout) :: config    ! CAM<->TUV-x map configuration
      logical,        intent(in)    :: do_euv    ! indicates whether to include
      !   extreme-UV rates in the mapping
      logical,        intent(out)   :: do_jno    ! indicates whether jno should be
      !   calculated
      integer,        intent(out)   :: jno_index ! index for jno in source photo
      !   rate array
      type(map_t),    intent(out)   :: map

      integer :: i_label, i_start, i_end
      type(string_t) :: str_label
      type(string_t), allocatable :: tuvx_labels(:), euv_labels(:), special_labels(:), &
         all_labels(:), cam_labels(:)

      ! ==================
      ! MOZART photo rates
      ! ==================
      allocate( cam_labels( phtcnt ) )
      do i_label = 1, phtcnt
         cam_labels( i_label ) = trim( rxt_tag_lst( i_label ) )
      end do

      ! =================
      ! TUV-x photo rates
      ! =================
      tuvx_labels = core%photolysis_reaction_labels( )

      ! ======================
      ! Extreme-UV photo rates
      ! ======================
      if( do_euv ) then
         allocate( euv_labels( neuv ) )
         do i_label = 1, neuv
            str_label = i_label
            euv_labels( i_label ) = "jeuv_"//str_label
         end do
      else
         allocate( euv_labels(0) )
      end if

      ! ===============================
      ! Special photo rate calculations
      ! ===============================
      do_jno = .false.
      jno_index = 0
      do i_label = 1, size( cam_labels )
         if( cam_labels( i_label ) == "jno" ) then
            do_jno = .true.
            exit
         end if
      end do
      if( do_jno ) then
         allocate( special_labels(1) )
         special_labels(1) = "jno"
      else
         allocate( special_labels(0) )
      end if

      ! ==========================
      ! Combine photo rate sources
      ! ==========================
      allocate( all_labels( size( tuvx_labels ) + size( euv_labels ) + &
         size( special_labels ) ) )
      i_end = 0
      if( size( tuvx_labels ) > 0 ) then
         i_start = i_end + 1
         i_end   = i_start + size( tuvx_labels ) - 1
         all_labels( i_start : i_end ) = tuvx_labels(:)
      end if
      if( size( euv_labels ) > 0 ) then
         i_start = i_end + 1
         i_end   = i_start + size( euv_labels ) - 1
         all_labels( i_start : i_end ) = euv_labels(:)
      end if
      if( size( special_labels ) > 0 ) then
         i_start = i_end + 1
         i_end   = i_start + size( special_labels ) - 1
         all_labels( i_start : i_end ) = special_labels(:)
         jno_index = i_start
      end if

      ! ==========
      ! Create map
      ! ==========
      map = map_t( config, all_labels, cam_labels )
      write(iulog,*)
      write(iulog,*) "TUV-x --> CAM-Chem photolysis rate constant map"
      call map%print( all_labels, cam_labels, iulog )

   end subroutine set_photo_rate_map

!================================================================================================

   !-----------------------------------------------------------------------
   ! Outputs diagnostic information for the current time step
   !-----------------------------------------------------------------------
   subroutine output_diagnostics( this, ncol, lchnk )

      use cam_history, only : outfld

      type(tuvx_ptr), intent(in) :: this
      integer,        intent(in) :: ncol  ! number of active columns on this thread
      integer,        intent(in) :: lchnk ! identifier for this thread

      integer :: i_diag

      if( .not. enable_diagnostics ) return

      do i_diag = 1, size( diagnostics )
         associate( diag => diagnostics( i_diag ) )
            call outfld( "tuvx_"//diag%name_, this%photo_rates_(:ncol,pver+1:2:-1,diag%index_), &
               ncol, lchnk )
         end associate
      end do

   end subroutine output_diagnostics

!================================================================================================

   !-----------------------------------------------------------------------
   ! Creates and loads a grid warehouse with grids that CAM will
   !   update at runtime
   !-----------------------------------------------------------------------
   function get_cam_grids( wavelength_path ) result( grids )

      use musica_assert,       only : assert_msg
      use musica_config,       only : config_t
      use tuvx_grid,           only : grid_t
      use tuvx_grid_factory,   only : grid_builder
      use tuvx_grid_from_host, only : grid_from_host_t
      use tuvx_grid_warehouse, only : grid_warehouse_t

      character(len=*),        intent(in) :: wavelength_path ! path to the wavelength data file
      class(grid_warehouse_t), pointer    :: grids ! collection of grids to be updated by CAM

      character(len=*), parameter :: my_name = "CAM grid creator"
      class(grid_t),    pointer   :: host_grid
      type(config_t)              :: config

      grids => grid_warehouse_t( )

      ! =========================
      ! heights above the surface
      ! =========================
      host_grid => grid_from_host_t( "height", "km", pver+1 )
      call grids%add( host_grid )
      deallocate( host_grid )

      ! ============================================
      ! wavelengths (will not be updated at runtime)
      ! ============================================
      call config%empty( )
      call config%add( "type", "from csv file", my_name )
      call config%add( "name", "wavelength", my_name )
      call config%add( "units", "nm", my_name )
      call config%add( "file path", wavelength_path, my_name )
      host_grid => grid_builder( config )
      call grids%add( host_grid )
      deallocate( host_grid )

   end function get_cam_grids

!================================================================================================

   !-----------------------------------------------------------------------
   ! Creates and loads a profile warehouse with profiles that CAM
   !   will update at runtime
   !-----------------------------------------------------------------------
   function get_cam_profiles( grids ) result( profiles )

      use tuvx_grid,              only : grid_t
      use tuvx_grid_warehouse,    only : grid_warehouse_t
      use tuvx_profile_from_host, only : profile_from_host_t
      use tuvx_profile_warehouse, only : profile_warehouse_t

      class(grid_warehouse_t),    intent(in) :: grids       ! CAM grids used in TUV-x
      class(profile_warehouse_t), pointer    :: profiles    ! collection of profiles to be updated by CAM

      class(profile_from_host_t), pointer :: host_profile
      class(grid_t),              pointer :: height, wavelength

      profiles   => profile_warehouse_t( )
      height     => grids%get_grid(     "height", "km" )
      wavelength => grids%get_grid( "wavelength", "nm" )

      ! ==================================
      ! Temperature profile on height grid
      ! ==================================
      host_profile => profile_from_host_t( "temperature", "K", height%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      ! =================================
      ! Surface albedo on wavelength grid
      ! =================================
      host_profile => profile_from_host_t( "surface albedo", "none", wavelength%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      ! ========================================
      ! Extraterrestrial flux on wavelength grid
      ! ========================================
      host_profile => profile_from_host_t( "extraterrestrial flux", "photon cm-2 s-1", &
         wavelength%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      ! ===========
      ! Air profile
      ! ===========
      host_profile => profile_from_host_t( "air", "molecule cm-3", height%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      ! ==========
      ! O3 profile
      ! ==========
      host_profile => profile_from_host_t( "O3", "molecule cm-3", height%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      ! ==========
      ! O2 profile
      ! ==========
      host_profile => profile_from_host_t( "O2", "molecule cm-3", height%size( ) )
      call profiles%add( host_profile )
      deallocate( host_profile )

      deallocate( height )
      deallocate( wavelength )

   end function get_cam_profiles

!================================================================================================

   !-----------------------------------------------------------------------
   ! Creates and loads a radiator warehouse with radiators that CAM will
   !    update at runtime
   !-----------------------------------------------------------------------
   function get_cam_radiators( grids ) result( radiators )

      use tuvx_grid,               only : grid_t
      use tuvx_grid_warehouse,     only : grid_warehouse_t
      use tuvx_radiator_from_host, only : radiator_from_host_t
      use tuvx_radiator_warehouse, only : radiator_warehouse_t

      type(grid_warehouse_t),      intent(in) :: grids     ! CAM grids used in TUV-x
      class(radiator_warehouse_t), pointer    :: radiators ! collection of radiators to be updated by CAM

      class(radiator_from_host_t), pointer :: host_radiator
      class(grid_t), pointer :: height, wavelength

      radiators  => radiator_warehouse_t( )
      height     => grids%get_grid(     "height", "km" )
      wavelength => grids%get_grid( "wavelength", "nm" )

      ! ================
      ! Aerosol radiator
      ! ================
      host_radiator => radiator_from_host_t( "aerosol", height, wavelength )
      call radiators%add( host_radiator )
      deallocate( host_radiator )

      ! ==============
      ! Cloud radiator
      ! ==============
      host_radiator => radiator_from_host_t( "clouds", height, wavelength )
      call radiators%add( host_radiator )
      deallocate( host_radiator )

      deallocate( height )
      deallocate( wavelength )

   end function get_cam_radiators

!================================================================================================

   !-----------------------------------------------------------------------
   ! Forms connections between CAM and TUV-x data structures
   !
   ! - creates updaters for each grid and profile that CAM will use
   !   to update TUV-x at each timestep.
   ! - allocates working arrays when needed for interpolation between
   !   grids
   ! - initializes CAM modules if necessary and sets parameters needed
   !   for runtime access of CAM data
   !
   !-----------------------------------------------------------------------
   subroutine create_updaters( this, grids, profiles, radiators, disable_aerosols, &
      disable_clouds )

      use ppgrid,                  only : pcols ! maximum number of columns
      use rad_constituents,        only : rad_cnst_get_info
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

      class(tuvx_ptr),             intent(inout) :: this
      class(grid_warehouse_t),     intent(in)    :: grids
      class(profile_warehouse_t),  intent(in)    :: profiles
      class(radiator_warehouse_t), intent(in)    :: radiators
      logical,                     intent(in)    :: disable_aerosols
      logical,                     intent(in)    :: disable_clouds

      class(grid_t),     pointer :: height, wavelength
      class(profile_t),  pointer :: host_profile
      class(radiator_t), pointer :: host_radiator
      integer                    :: n_modes
      logical                    :: found

      ! =============
      ! Grid updaters
      ! =============

      height => grids%get_grid( "height", "km" )
      this%grids_( GRID_INDEX_HEIGHT ) = this%core_%get_updater( height, found )
      call assert( 213798815, found )

      ! ============================================
      ! wavelength grid cannot be updated at runtime
      ! ============================================
      wavelength => grids%get_grid( "wavelength", "nm" )
      allocate( this%wavelength_edges_(      wavelength%size( ) + 1 ) )
      allocate( this%wavelength_values_(     wavelength%size( ) + 1 ) )
      allocate( this%wavelength_mid_values_( wavelength%size( )     ) )
      this%wavelength_edges_(:) = wavelength%edge_(:)

      ! ==============================
      ! optical property working array
      ! ==============================
      allocate( this%optical_depth_(            pcols, height%size( ), wavelength%size( ) ) )
      allocate( this%single_scattering_albedo_( pcols, height%size( ), wavelength%size( ) ) )
      allocate( this%asymmetry_factor_(         pcols, height%size( ), wavelength%size( ) ) )

      deallocate( height )
      deallocate( wavelength )

      ! ================
      ! Profile updaters
      ! ================

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

      ! =================
      ! radiator updaters
      ! =================

      ! ====================================================================
      ! determine if aerosol optical properties will be available, and if so
      ! intialize the aerosol optics module
      ! ====================================================================
      call rad_cnst_get_info( 0, nmodes = n_modes )
      if( n_modes > 0 .and. .not. do_aerosol .and. .not. disable_aerosols ) then
         do_aerosol = .true.
         ! TODO update to use new aerosol_optics class
         ! call modal_aer_opt_init( )
      else
         ! TODO are there default aerosol optical properties that should be used
         !      when an aerosol module is not available?
         this%optical_depth_(:,:,:)            = 0.0_r8
         this%single_scattering_albedo_(:,:,:) = 0.0_r8
         this%asymmetry_factor_(:,:,:)         = 0.0_r8
      end if
      host_radiator => radiators%get_radiator( "aerosol" )
      this%radiators_( RADIATOR_INDEX_AEROSOL ) = &
         this%core_%get_updater( host_radiator, found )
      call assert( 675200430, found )
      nullify( host_radiator )

      ! =====================================
      ! get an updater for the cloud radiator
      ! =====================================
      do_clouds = .not. disable_clouds
      host_radiator => radiators%get_radiator( "clouds" )
      this%radiators_( RADIATOR_INDEX_CLOUDS ) = &
         this%core_%get_updater( host_radiator, found )
      call assert( 993715720, found )
      nullify( host_radiator )

   end subroutine create_updaters

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets the height values in TUV-x for the given column
   !
   ! NOTE: This function must be called before updating any profile data.
   !
   !  CAM to TUV-x height grid mapping
   !
   !  TUV-x heights are "bottom-up" and require atmospheric constituent
   !  concentrations at interfaces. Therefore, CAM mid-points are used as
   !  TUV-x grid interfaces, with an additional layer introduced between
   !  the surface and the lowest CAM mid-point, and a layer at the
   !  top of the TUV-x grid to hold species densities above the top CAM
   !  mid-point.
   !
   !  ---- (interface)  ===== (mid-point)
   !
   !        CAM                                  TUV-x
   ! ------(top)------ i_int = 1           -------(top)------ i_int = pver + 2
   ! ************************ (exo values) *****************************
   !                                       ================== i_mid = pver + 1
   ! ================= i_mid = 1           ------------------ i_int = pver + 1
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
   subroutine set_heights( this, i_col, ncol, height_mid, height_int )

      use ppgrid,         only : pcols ! maximum number of columns

      class(tuvx_ptr), intent(inout) :: this                     ! TUV-x calculator
      integer,         intent(in)    :: i_col                    ! column to set conditions for
      integer,         intent(in)    :: ncol                     ! number of colums to calculated photolysis for
      real(r8),        intent(in)    :: height_mid(ncol,pver)    ! height above the surface at mid-points (km)
      real(r8),        intent(in)    :: height_int(ncol,pver+1)  ! height above the surface at interfaces (km)

      integer :: i_level
      real(r8) :: edges(pver+2)
      real(r8) :: mid_points(pver+1)

      edges(1) = height_int(i_col,pver+1)
      edges(2:pver+1) = height_mid(i_col,pver:1:-1)
      edges(pver+2) = height_int(i_col,1)
      mid_points(1) = ( height_mid(i_col,pver) - height_int(i_col,pver+1) ) * 0.5_r8 &
         + height_int(i_col,pver+1)
      mid_points(2:pver) = height_int(i_col,pver:2:-1)
      mid_points(pver+1) = 0.5_r8 * ( edges(pver+1) + edges(pver+2) )
      call this%grids_( GRID_INDEX_HEIGHT )%update( edges = edges, mid_points = mid_points )
      this%height_delta_(1:pver+1) = edges(2:pver+2) - edges(1:pver+1)

   end subroutine set_heights

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets the temperatures in TUV-x for the given column
   !
   ! See description of `set_heights` for CAM <-> TUV-x vertical grid
   ! mapping.
   !-----------------------------------------------------------------------
   subroutine set_temperatures( this, i_col, temperature_mid, surface_temperature )

      use ppgrid,         only : pcols ! maximum number of columns

      class(tuvx_ptr), intent(inout) :: this                        ! TUV-x calculator
      integer,         intent(in)    :: i_col                       ! column to set conditions for
      real(r8),        intent(in)    :: temperature_mid(pcols,pver) ! midpoint temperature (K)
      real(r8),        intent(in)    :: surface_temperature(pcols)  ! surface temperature (K)

      real(r8) :: edges(pver+2)

      edges(1) = surface_temperature(i_col)
      edges(2:pver+1) = temperature_mid(i_col,pver:1:-1)
      edges(pver+2) = temperature_mid(i_col,1) ! Use upper mid-point temperature for top edge
      call this%profiles_( PROFILE_INDEX_TEMPERATURE )%update( edge_values = edges )

   end subroutine set_temperatures

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets the surface albedo in TUV-x for the given column
   !
   ! CAM uses a single value for surface albedo at all wavelengths
   !-----------------------------------------------------------------------
   subroutine set_surface_albedo( this, i_col, surface_albedo )

      use ppgrid,         only : pcols ! maximum number of columns

      class(tuvx_ptr), intent(inout) :: this                  ! TUV-x calculator
      integer,         intent(in)    :: i_col                 ! column to set conditions for
      real(r8),        intent(in)    :: surface_albedo(pcols) ! surface albedo (unitless)

      this%wavelength_values_(:) = surface_albedo( i_col )
      call this%profiles_( PROFILE_INDEX_ALBEDO )%update( &
         edge_values = this%wavelength_values_(:) )

   end subroutine set_surface_albedo

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets the extraterrestrial flux in TUV-x for the given column
   !
   ! Extraterrestrial flux is read from data files and interpolated to the
   ! TUV-x wavelength grid. CAM ET Flux values are multiplied by the
   ! width of the wavelength bins to get the TUV-x units of photon cm-2 s-1
   !
   ! NOTE: TUV-x only uses mid-point values for ET Flux
   !-----------------------------------------------------------------------
   subroutine set_et_flux( this )

      use mo_util,          only : rebin
      use solar_irrad_data, only : nbins,   & ! number of wavelength bins
         we,      & ! wavelength bin edges
         sol_etf    ! extraterrestrial flux
      !   (photon cm-2 nm-1 s-1)

      class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator

      real(r8) :: et_flux_orig(nbins)
      integer  :: n_tuvx_bins, i_bin

      ! ===============================================
      ! regrid normalized flux to TUV-x wavelength grid
      !================================================
      et_flux_orig(:) = sol_etf(:)
      n_tuvx_bins = size(this%wavelength_mid_values_)
      call rebin( nbins, n_tuvx_bins, we, this%wavelength_edges_, et_flux_orig, &
         this%wavelength_mid_values_ )

      ! ========================================================
      ! convert normalized flux to flux on TUV-x wavelength grid
      ! ========================================================
      this%wavelength_mid_values_(:) = this%wavelength_mid_values_(:) * &
         ( this%wavelength_edges_(2:n_tuvx_bins+1) - &
         this%wavelength_edges_(1:n_tuvx_bins) )

      ! ====================================
      ! estimate unused edge values for flux
      ! ====================================
      this%wavelength_values_(1)  = this%wavelength_mid_values_(1) - &
         ( this%wavelength_mid_values_(2) - &
         this%wavelength_mid_values_(1) ) * 0.5_r8
      do i_bin = 2, n_tuvx_bins
         this%wavelength_values_(i_bin) = this%wavelength_mid_values_(i_bin-1) + &
            ( this%wavelength_mid_values_(i_bin) - &
            this%wavelength_mid_values_(i_bin-1) ) * 0.5_r8
      end do
      this%wavelength_values_(n_tuvx_bins+1) = &
         this%wavelength_mid_values_(n_tuvx_bins) + &
         ( this%wavelength_mid_values_(n_tuvx_bins) - &
         this%wavelength_mid_values_(n_tuvx_bins-1) ) * 0.5_r8

      ! ============================
      ! update TUV-x ET flux profile
      ! ============================
      call this%profiles_( PROFILE_INDEX_ET_FLUX )%update( &
         mid_point_values = this%wavelength_mid_values_, &
         edge_values      = this%wavelength_values_)

      ! ======================================================================
      ! rebin extraterrestrial flux to MS93 grid for use with jno calculations
      ! ======================================================================
      call rebin( nbins, NUM_BINS_MS93, we, WAVELENGTH_EDGES_MS93, et_flux_orig, &
         this%et_flux_ms93_ )

   end subroutine set_et_flux

!================================================================================================

   !-----------------------------------------------------------------------
   ! Sets the profiles of optically active atmospheric constituents
   !   in TUV-x for the given column
   !
   ! See `set_height` for a description of the CAM <-> TUV-x vertical grid
   ! mapping.
   !
   ! Above layer densities are calculated using a scale height for air
   ! and pre-calculated values for O2 and O3
   !-----------------------------------------------------------------------
   subroutine set_radiator_profiles( this, i_col, ncol, fixed_species_conc, species_vmr, &
      exo_column_conc, delta_pressure, cloud_fraction, &
      liquid_water_content )

      use chem_mods, only : gas_pcnst, & ! number of non-fixed species
         nfs,       & ! number of fixed species
         nabscol,   & ! number of absorbing species (radiators)
         indexm       ! index for air density in fixed species array

      class(tuvx_ptr), intent(inout) :: this  ! TUV-x calculator
      integer,         intent(in)    :: i_col ! column to set conditions for
      integer,         intent(in)    :: ncol  ! number of columns
      real(r8),        intent(in)    :: fixed_species_conc(ncol,pver,max(1,nfs))    ! fixed species densities
      !   (molecule cm-3)
      real(r8),        intent(in)    :: species_vmr(ncol,pver,max(1,gas_pcnst))     ! species volume mixing
      !   ratios (mol mol-1)
      real(r8),        intent(in)    :: exo_column_conc(ncol,0:pver,max(1,nabscol)) ! above column densities
      !   (molecule cm-2)
      real(r8),        intent(in)    :: delta_pressure(ncol,pver)       ! pressure delta about midpoints (Pa)
      real(r8),        intent(in)    :: cloud_fraction(ncol,pver)       ! cloud fraction (unitless)
      real(r8),        intent(in)    :: liquid_water_content(ncol,pver) ! liquid water content (kg/kg)

      integer  :: i_level
      real(r8) :: tmp(pver)
      real(r8) :: tau(pver+1, size(this%wavelength_mid_values_))
      real(r8) :: edges(pver+2), densities(pver+1)
      real(r8) :: exo_val
      real(r8), parameter :: rgrav = 1.0_r8 / 9.80616_r8 ! reciprocal of acceleration by gravity (s/m)
      real(r8), parameter :: km2cm = 1.0e5_r8 ! conversion from km to cm

      ! ===========
      ! air profile
      ! ===========
      edges(1) = fixed_species_conc(i_col,pver,indexm)
      edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,indexm)
      edges(pver+2) = fixed_species_conc(i_col,1,indexm) ! use upper mid-point value for top edge
      densities(1:pver+1) = this%height_delta_(1:pver+1) * km2cm * &
         sqrt(edges(1:pver+1)) * sqrt(edges(2:pver+2))
      call this%profiles_( PROFILE_INDEX_AIR )%update( &
         edge_values = edges, layer_densities = densities, &
         scale_height = 8.01_r8 ) ! scale height in [km]

      ! ==========
      ! O2 profile
      ! ==========
      if( is_fixed_O2 ) then
         edges(1) = fixed_species_conc(i_col,pver,index_O2)
         edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,index_O2)
         edges(pver+2) = fixed_species_conc(i_col,1,index_O2)
      else if( index_O2 > 0 ) then
         edges(1) = species_vmr(i_col,pver,index_O2) * &
            fixed_species_conc(i_col,pver,indexm)
         edges(2:pver+1) = species_vmr(i_col,pver:1:-1,index_O2) * &
            fixed_species_conc(i_col,pver:1:-1,indexm)
         edges(pver+2) = species_vmr(i_col,1,index_O2) * &
            fixed_species_conc(i_col,1,indexm)
      else
         edges(:) = 0.0_r8
      end if
      densities(1:pver+1) = this%height_delta_(1:pver+1) * km2cm * &
         sqrt(edges(1:pver+1)) * sqrt(edges(2:pver+2))
      call this%profiles_( PROFILE_INDEX_O2 )%update( &
         edge_values = edges, layer_densities = densities, &
         scale_height = 7.0_r8 )

      ! ==========
      ! O3 profile
      ! ==========
      if( is_fixed_O3 ) then
         edges(1) = fixed_species_conc(i_col,pver,index_O3)
         edges(2:pver+1) = fixed_species_conc(i_col,pver:1:-1,index_O3)
         edges(pver+2) = fixed_species_conc(i_col,1,index_O3)
      else if( index_O3 > 0 ) then
         edges(1) = species_vmr(i_col,pver,index_O3) * &
            fixed_species_conc(i_col,pver,indexm)
         edges(2:pver+1) = species_vmr(i_col,pver:1:-1,index_O3) * &
            fixed_species_conc(i_col,pver:1:-1,indexm)
         edges(pver+2) = species_vmr(i_col,1,index_O3) * &
            fixed_species_conc(i_col,1,indexm)
      else
         edges(:) = 0.0_r8
      end if
      if( nabscol >= 1 ) then
         densities(1) = 0.5_r8 * exo_column_conc(i_col,pver,1)
         densities(2:pver) = 0.5_r8 * ( exo_column_conc(i_col,pver-1:1:-1,1) &
            + exo_column_conc(i_col,pver:2:-1,1) )
         densities(pver+1) = exo_column_conc(i_col,0,1) &
            + 0.5_r8 * exo_column_conc(i_col,1,1)
         call this%profiles_( PROFILE_INDEX_O3 )%update( &
            edge_values = edges, layer_densities = densities, &
            exo_density = exo_column_conc(i_col,0,1) )
      else
         densities(1:pver+1) = this%height_delta_(1:pver+1) * km2cm * &
            ( edges(1:pver+1) + edges(2:pver+2) ) * 0.5_r8
         call this%profiles_( PROFILE_INDEX_O3 )%update( &
            edge_values = edges, layer_densities = densities, &
            scale_height = 7.0_r8 )
      end if

      ! ===============
      ! aerosol profile
      ! ===============
      if( do_aerosol ) then
         call this%radiators_( RADIATOR_INDEX_AEROSOL )%update( &
            optical_depths            = this%optical_depth_(i_col,:,:), &
            single_scattering_albedos = this%single_scattering_albedo_(i_col,:,:), &
            asymmetry_factors         = this%asymmetry_factor_(i_col,:,:) )
      end if

      ! =============
      ! cloud profile
      ! =============
      if( do_clouds ) then
         ! ===================================================
         ! estimate cloud optical depth as:
         !    liquid_water_path * 0.155 * cloud_fraction^(1.5)
         ! ===================================================
         associate( clouds => cloud_fraction(i_col,:) )
            where( clouds(:) /= 0.0_r8 )
               tmp(:) = ( rgrav * liquid_water_content(i_col,:) * delta_pressure(i_col,:) &
                  * 1.0e3_r8 / clouds(:) ) * 0.155_r8 * clouds(:)**1.5_r8
            elsewhere
               tmp(:) = 0.0_r8
            end where
         end associate
         do i_level = 1, pver
            tau(i_level,:) = tmp(pver-i_level+1)
         end do
         tau(pver+1,:) = 0.0_r8
         call this%radiators_( RADIATOR_INDEX_CLOUDS )%update( optical_depths = tau )
      end if

   end subroutine set_radiator_profiles

!================================================================================================

   !-----------------------------------------------------------------------
   ! Updates working arrays of aerosol optical properties for all
   !   columns from the aerosol package
   !-----------------------------------------------------------------------
   subroutine get_aerosol_optical_properties( this, state, pbuf )

      use aer_rad_props,    only : aer_rad_props_sw
      use mo_util,          only : rebin
      use physics_types,    only : physics_state
      use physics_buffer,   only : physics_buffer_desc
      use ppgrid,           only : pcols          ! maximum number of columns
      use radconstants,     only : nswbands, &    ! Number of CAM shortwave radiation bands
         get_sw_spectral_boundaries

      class(tuvx_ptr),                    intent(inout) :: this  ! TUV-x calculator
      type(physics_state),       target,  intent(in)    :: state
      type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)

      real(r8) :: wavelength_edges(nswbands+1)       ! CAM radiation wavelength grid edges [nm]
      real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
      real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
      real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
      real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau
      real(r8) :: low_bound(nswbands)  ! lower bound of CAM wavenumber bins
      real(r8) :: high_bound(nswbands) ! upper bound of CAM wavenumber bins
      integer  :: n_night              ! number of night columns
      integer  :: idx_night(pcols)     ! indices of night columns
      integer  :: n_tuvx_bins          ! number of TUV-x wavelength bins
      integer  :: i_col, i_level       ! column and level indices

      ! ============================================
      ! do nothing if no aerosol module is available
      ! ============================================
      if( .not. do_aerosol ) return

      ! TODO just assume all daylight columns for now
      !      can adjust later if necessary
      n_night = 0
      idx_night(:) = 0

      ! ===========================================================
      ! get aerosol optical properties on native CAM radiation grid
      ! ===========================================================
      call aer_rad_props_sw( 0, state, pbuf, n_night, idx_night, &
         aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f )

      ! =========================================================================
      ! Convert CAM wavenumber grid to wavelength grid and re-order optics arrays
      ! NOTE: CAM wavenumber grid is continuous and increasing, except that the
      !       last bin should be moved to the just before the first bin (!?!)
      ! =========================================================================
      call get_sw_spectral_boundaries( low_bound, high_bound, 'nm' )
      wavelength_edges(1:nswbands-1) = high_bound(nswbands-1:1:-1)
      wavelength_edges(nswbands    ) = high_bound(nswbands       )
      wavelength_edges(nswbands+1  ) = low_bound( nswbands       )
      call reorder_optics_array(     aer_tau )
      call reorder_optics_array(   aer_tau_w )
      call reorder_optics_array( aer_tau_w_g )

      ! =============================================================
      ! regrid optical properties to TUV-x wavelength and height grid
      ! =============================================================
      ! TODO is this the correct regridding scheme to use?
      n_tuvx_bins = size(this%wavelength_mid_values_)
      do i_col = 1, pcols
         do i_level = 1, pver
            call rebin(nswbands, n_tuvx_bins, wavelength_edges, this%wavelength_edges_, &
               aer_tau(i_col,pver-i_level,:), this%optical_depth_(i_col,i_level,:))
            call rebin(nswbands, n_tuvx_bins, wavelength_edges, this%wavelength_edges_, &
               aer_tau_w(i_col,pver-i_level,:), this%single_scattering_albedo_(i_col,i_level,:))
            call rebin(nswbands, n_tuvx_bins, wavelength_edges, this%wavelength_edges_, &
               aer_tau_w_g(i_col,pver-i_level,:), this%asymmetry_factor_(i_col,i_level,:))
         end do
         this%optical_depth_(i_col,pver+1,:) = &
            this%optical_depth_(i_col,pver,:)
         this%single_scattering_albedo_(i_col,pver+1,:) = &
            this%single_scattering_albedo_(i_col,pver,:)
         this%asymmetry_factor_(i_col,pver+1,:) = &
            this%asymmetry_factor_(i_col,pver,:)
      end do

      ! ================================================================
      ! back-calculate the single scattering albedo and asymmetry factor
      ! ================================================================
      associate( tau   => this%optical_depth_, &
         omega => this%single_scattering_albedo_, &
         g     => this%asymmetry_factor_ )
         where(omega > 0.0_r8)
            g = g / omega
         elsewhere
            g = 0.0_r8
         end where
         where(tau > 0.0_r8)
            omega = omega / tau
         elsewhere
            omega = 0.0_r8
         end where
      end associate

   end subroutine get_aerosol_optical_properties

!================================================================================================

   !-----------------------------------------------------------------------
   ! Reorders elements of an optical property array for conversion
   !   from wavenumber to wavelength grid and out-of-order final
   !   element in CAM wavenumber grid
   !-----------------------------------------------------------------------
   subroutine reorder_optics_array( optics_array )

      use ppgrid,           only : pcols    ! maximum number of columns
      use radconstants,     only : nswbands ! Number of CAM shortwave radiation bands

      real(r8), intent(inout) :: optics_array(pcols,0:pver,nswbands) ! optics array to reorder

      real(r8) :: working(pcols,0:pver,nswbands) ! working array

      working(:,:,:) = optics_array(:,:,:)
      optics_array(:,:,1:nswbands-1) = working(:,:,nswbands-1:1:-1)

   end subroutine reorder_optics_array

!================================================================================================

   !-----------------------------------------------------------------------
   ! Calculates extreme-UV ionization rates
   !
   ! NOTE This never includes an above-column layer
   !-----------------------------------------------------------------------
   subroutine calculate_euv_rates( solar_zenith_angle, fixed_species_conc, &
      species_vmr, height_mid, height_int, euv_rates )

      use chem_mods, only : gas_pcnst, & ! number of non-fixed species
         nfs, &       ! number of fixed species
         indexm       ! index for air density in fixed species array
      use mo_jeuv,   only : jeuv, neuv   ! number of extreme-UV rates
      use ref_pres,  only : ptop_ref     ! pressure at the top of the column (Pa)

      real(r8), intent(in)  :: solar_zenith_angle ! degrees
      real(r8), intent(in)  :: fixed_species_conc(pver,max(1,nfs))   ! fixed species densities
      !   (molecule cm-3)
      real(r8), intent(in)  :: species_vmr(pver,max(1,gas_pcnst))    ! species volume mixing
      !   ratios (mol mol-1)
      real(r8), intent(in)  :: height_mid(pver)     ! height at mid-points (km)
      real(r8), intent(in)  :: height_int(pver+1)   ! height at interfaces (km)
      real(r8), intent(out) :: euv_rates(pver,neuv) ! calculated extreme-UV rates

      real(r8) :: o_dens(pver), o2_dens(pver), n2_dens(pver), height_arg(pver)

      ! ==========
      ! N2 density
      ! ==========
      if( is_fixed_N2 ) then
         n2_dens(:) = fixed_species_conc(:pver,index_N2)
      else
         n2_dens(:) = species_vmr(:pver,index_N2) * fixed_species_conc(:pver,indexm)
      end if

      ! =========
      ! O density
      ! =========
      if( is_fixed_O ) then
         o_dens(:) = fixed_species_conc(:pver,index_O)
      else
         o_dens(:) = species_vmr(:pver,index_O) * fixed_species_conc(:pver,indexm)
      end if

      ! ==========
      ! O2 density
      ! ==========
      if( is_fixed_O2 ) then
         o2_dens(:) = fixed_species_conc(:pver,index_O2)
      else
         o2_dens(:) = species_vmr(:pver,index_O2) * fixed_species_conc(:pver,indexm)
      end if

      ! =======================
      ! special height argument
      ! =======================
      height_arg(:) = height_mid(:)

      call jeuv( pver, solar_zenith_angle, o_dens, o2_dens, n2_dens, height_arg, &
                 euv_rates(pver:1:-1,:) )

   end subroutine calculate_euv_rates

!================================================================================================

   !-----------------------------------------------------------------------
   ! Calculates NO photolysis rates
   !
   ! NOTE: Always includes an above-column layer
   !-----------------------------------------------------------------------
   subroutine calculate_jno( solar_zenith_angle, et_flux, fixed_species_conc, species_vmr, &
      height_int, jno )

      use chem_mods, only : gas_pcnst, & ! number of non-fixed species
         nfs, &       ! number of fixed species
         indexm       ! index for air density in fixed species array
      use mo_jshort, only : sphers, slant_col, calc_jno
      use ref_pres,  only : ptop_ref     ! pressure at the top of the column (Pa)

      real(r8), intent(in)  :: solar_zenith_angle                  ! degrees
      real(r8), intent(in)  :: et_flux(NUM_BINS_MS93)              ! extraterrestrial flux MS93 grid
      !   (photon cm-2 nm-1 s-1)
      real(r8), intent(in)  :: fixed_species_conc(pver,max(1,nfs)) ! fixed species densities
      !   (molecule cm-3)
      real(r8), intent(in)  :: species_vmr(pver,max(1,gas_pcnst))  ! species volume mixing
      !   ratios (mol mol-1)
      real(r8), intent(in)  :: height_int(pver+1)                  ! height at interfaces (km)
      real(r8), intent(out) :: jno(pver)                           ! calculated NO rate

      ! species column densities (molecule cm-3)
      real(kind=r8) :: n2_dens(pver+1), o2_dens(pver+1), o3_dens(pver+1), no_dens(pver+1)
      ! species slant column densities (molecule cm-2)
      real(kind=r8) :: o2_slant(pver+1), o3_slant(pver+1), no_slant(pver+1)
      ! working photo rate array
      real(kind=r8) :: work_jno(pver+1)
      ! parameters needed to calculate slant column densities
      ! (see sphers routine description for details)
      integer       :: nid(pver+1)
      real(kind=r8) :: dsdh(0:pver+1,pver+1)
      ! layer thickness (cm)
      real(kind=r8) :: delz(pver+1)
      ! conversion from km to cm
      real(kind=r8), parameter :: km2cm = 1.0e5_r8

      ! ==========
      ! N2 density
      ! ==========
      if( is_fixed_N2 ) then
         n2_dens(2:) = fixed_species_conc(:pver,index_N2)
      else
         n2_dens(2:) = species_vmr(:pver,index_N2) * fixed_species_conc(:pver,indexm)
      end if
      n2_dens(1) = n2_dens(2) * 0.9_r8

      ! ==========
      ! O2 density
      ! ==========
      if( is_fixed_O2 ) then
         o2_dens(2:) = fixed_species_conc(:pver,index_O2)
      else
         o2_dens(2:) = species_vmr(:pver,index_O2) * fixed_species_conc(:pver,indexm)
      end if
      o2_dens(1) = o2_dens(2) * 7.0_r8 / ( height_int(1) - height_int(2) )

      ! ==========
      ! O3 density
      ! ==========
      if( is_fixed_O3 ) then
         o3_dens(2:) = fixed_species_conc(:pver,index_O3)
      else
         o3_dens(2:) = species_vmr(:pver,index_O3) * fixed_species_conc(:pver,indexm)
      end if
      o3_dens(1) = o3_dens(2) * 7.0_r8 / ( height_int(1) - height_int(2) )

      ! ==========
      ! NO density
      ! ==========
      if( is_fixed_NO ) then
         no_dens(2:) = fixed_species_conc(:pver,index_NO)
      else
         no_dens(2:) = species_vmr(:pver,index_NO) * fixed_species_conc(:pver,indexm)
      end if
      no_dens(1) = no_dens(2) * 0.9_r8

      ! ================================
      ! calculate slant column densities
      ! ================================
      call sphers( pver+1, height_int, solar_zenith_angle, dsdh, nid )
      delz(1:pver) = km2cm * ( height_int(1:pver) - height_int(2:pver+1) )
      call slant_col( pver+1, delz, dsdh, nid, o2_dens, o2_slant )
      call slant_col( pver+1, delz, dsdh, nid, o3_dens, o3_slant )
      call slant_col( pver+1, delz, dsdh, nid, no_dens, no_slant )

      ! =========================================
      ! calculate the NO photolysis rate constant
      ! =========================================
      call calc_jno( pver+1, et_flux, n2_dens, o2_slant, o3_slant, no_slant, work_jno )
      jno(:) = work_jno(pver:1:-1)

   end subroutine calculate_jno

!================================================================================================

end module mo_tuvx
