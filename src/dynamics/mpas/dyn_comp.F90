module dyn_comp

! CAM component interfaces to the MPAS Dynamical Core

use shr_kind_mod,       only: r8=>shr_kind_r8
use spmd_utils,         only: masterproc, mpicom, npes
use physconst,          only: pi, gravit, rair, cpair

use pmgrid,             only: plev, plevp
use constituents,       only: pcnst, cnst_name, cnst_is_a_water_species, cnst_read_iv
use const_init,         only: cnst_init_default

use cam_control_mod,    only: initial_run
use cam_initfiles,      only: initial_file_get_id, topo_file_get_id

use cam_grid_support,   only: cam_grid_id, &
                              cam_grid_get_latvals, cam_grid_get_lonvals
use cam_map_utils,      only: iMap

use inic_analytic,      only: analytic_ic_active, dyn_set_inic_col
use dyn_tests_utils,    only: vcoord=>vc_height

use cam_history,        only: addfld, horiz_only
use string_utils,       only: date2yyyymmdd, sec2hms, int2str

use ncdio_atm,          only: infld
use pio,                only: file_desc_t
use cam_pio_utils,      only: clean_iodesc_list

use time_manager,       only: get_start_date, get_stop_date, get_run_duration, &
                              timemgr_get_calendar_cf, get_step_size

use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun

use mpas_timekeeping,   only : MPAS_TimeInterval_type
use cam_mpas_subdriver, only: cam_mpas_global_sum_real
use cam_budget,         only: cam_budget_em_snapshot, cam_budget_em_register


use phys_control,       only: use_gw_front, use_gw_front_igw

implicit none
private
save

public :: &
   dyn_import_t, &
   dyn_export_t, &
   dyn_readnl,   &
   dyn_register, &
   dyn_init,     &
   dyn_run,      &
   dyn_final

! Note that the fields in the import and export states are pointers into the MPAS dycore internal
! data structures.  These fields have the order of the vertical and horizontal dimensions swapped
! relative to the CAM convention, as well as having the vertical indices ordered from bottom to top
! of atm.  An exception is that the export state contains two fields, pmiddry and pintdry, not managed
! by the MPAS infrastructure.  These fields are only used by the physics package and are computed
! in the dp_coupling module.

type dyn_import_t
   !
   ! Number of cells, edges, vertices, and vertical layers in this block
   !
   integer                        :: nCells           ! Number of cells, including halo cells
   integer                        :: nEdges           ! Number of edges, including halo edges
   integer                        :: nVertices        ! Number of vertices, including halo vertices
   integer                        :: nVertLevels      ! Number of vertical layers

   integer                        :: nCellsSolve      ! Number of cells, excluding halo cells
   integer                        :: nEdgesSolve      ! Number of edges, excluding halo edges
   integer                        :: nVerticesSolve   ! Number of vertices, excluding halo vertices

   !
   ! State that is directly prognosed by the dycore
   !
   real(r8), dimension(:,:),   pointer :: uperp   ! Normal velocity at edges [m/s]  (nver,nedge)
   real(r8), dimension(:,:),   pointer :: w       ! Vertical velocity [m/s]        (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: theta_m ! Moist potential temperature [K]  (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_zz  ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), dimension(:,:,:), pointer :: tracers ! Tracers [kg/kg dry air]       (nq,nver,ncol)

   !
   ! Index map between MPAS tracers and CAM constituents
   !
   integer, dimension(:), pointer :: mpas_from_cam_cnst => null() ! indices into CAM constituent array

   !
   ! Base state variables
   !
   real(r8), dimension(:,:),   pointer :: rho_base    ! Base-state dry air density [kg/m^3]  (nver,ncol)
   real(r8), dimension(:,:),   pointer :: theta_base  ! Base-state potential temperature [K] (nver,ncol)

   !
   ! Indices of tracers
   !
   integer                             :: index_qv  ! Index in tracers array of water vapor
                                                    ! mixing ratio

   !
   ! Invariant -- the vertical coordinate in MPAS-A is a height coordinate
   !
   real(r8), dimension(:,:),   pointer :: zint    ! Geometric height [m]
                                                  ! at layer interfaces            (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: zz      ! Vertical coordinate metric [dimensionless]
                                                  ! at layer midpoints               (nver,ncol)
   real(r8), dimension(:),     pointer :: fzm     ! Interp weight from k layer midpoint to k layer
                                                  ! interface [dimensionless]             (nver)
   real(r8), dimension(:),     pointer :: fzp     ! Interp weight from k-1 layer midpoint to k
                                                  ! layer interface [dimensionless]       (nver)
   !
   ! Invariant -- cell area
   !
   real(r8), dimension(:),     pointer :: areaCell ! cell area (m^2)


   !
   ! Invariant -- needed to compute edge-normal velocities
   !
   real(r8), dimension(:,:),   pointer :: east    ! Cartesian components of unit east vector
                                                  ! at cell centers [dimensionless]     (3,ncol)
   real(r8), dimension(:,:),   pointer :: north   ! Cartesian components of unit north vector
                                                  ! at cell centers [dimensionless]     (3,ncol)
   real(r8), dimension(:,:),   pointer :: normal  ! Cartesian components of the vector normal
                                                  ! to an edge and tangential to the surface
                                                  ! of the sphere [dimensionless]       (3,ncol)
   integer, dimension(:,:), pointer :: cellsOnEdge ! Indices of cells separated by an edge (2,nedge)


   !
   ! State that may be directly derived from dycore prognostic state
   !
   real(r8), dimension(:,:),   pointer :: theta   ! Potential temperature [K]        (nver,ncol)
   real(r8), dimension(:,:),   pointer :: exner   ! Exner function [-]               (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ux      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), dimension(:,:),   pointer :: uy      ! Meridional veloc at center [m/s] (nver,ncol)

   !
   ! Tendencies from physics
   !
   real(r8), dimension(:,:),   pointer :: ru_tend ! Normal horizontal momentum tendency
                                                  ! from physics [kg/m^2/s]         (nver,nedge)
   real(r8), dimension(:,:),   pointer :: rtheta_tend ! Tendency of rho*theta/zz
                                                      ! from physics [kg K/m^3/s]    (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_tend ! Dry air density tendency
                                                   ! from physics [kg/m^3/s]         (nver,ncol)
end type dyn_import_t

type dyn_export_t
   !
   ! Number of cells, edges, vertices, and vertical layers in this block
   !
   integer                        :: nCells           ! Number of cells, including halo cells
   integer                        :: nEdges           ! Number of edges, including halo edges
   integer                        :: nVertices        ! Number of vertices, including halo vertices
   integer                        :: nVertLevels      ! Number of vertical layers

   integer                        :: nCellsSolve      ! Number of cells, excluding halo cells
   integer                        :: nEdgesSolve      ! Number of edges, excluding halo edges
   integer                        :: nVerticesSolve   ! Number of vertices, excluding halo vertices

   !
   ! State that is directly prognosed by the dycore
   !
   real(r8), dimension(:,:),   pointer :: uperp   ! Normal velocity at edges [m/s]  (nver,nedge)
   real(r8), dimension(:,:),   pointer :: w       ! Vertical velocity [m/s]        (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: theta_m ! Moist potential temperature [K]  (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_zz  ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), dimension(:,:,:), pointer :: tracers ! Tracers [kg/kg dry air]       (nq,nver,ncol)

   !
   ! Indices of tracers
   !
   integer                             :: index_qv  ! Index in tracers array of water vapor
                                                    ! mixing ratio

   !
   ! Index map between MPAS tracers and CAM constituents
   !
   integer, dimension(:), pointer :: cam_from_mpas_cnst => null() ! indices into MPAS tracers array

   !
   ! Invariant -- the vertical coordinate in MPAS-A is a height coordinate
   !
   real(r8), dimension(:,:),   pointer :: zint    ! Geometric height [m]
                                                  ! at layer interfaces            (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: zz      ! Vertical coordinate metric [dimensionless]
                                                  ! at layer midpoints               (nver,ncol)
   real(r8), dimension(:),     pointer :: fzm     ! Interp weight from k layer midpoint to k layer
                                                  ! interface [dimensionless]             (nver)
   real(r8), dimension(:),     pointer :: fzp     ! Interp weight from k-1 layer midpoint to k
   !
   ! Invariant -- needed for computing the frontogenesis function
   !
   real(r8), dimension(:,:),   pointer :: defc_a
   real(r8), dimension(:,:),   pointer :: defc_b
   real(r8), dimension(:,:),   pointer :: cell_gradient_coef_x
   real(r8), dimension(:,:),   pointer :: cell_gradient_coef_y
   real(r8), dimension(:,:),   pointer :: edgesOnCell_sign
   real(r8), dimension(:),     pointer :: dvEdge
   real(r8), dimension(:),     pointer :: areaCell ! cell area (m^2)

   integer, dimension(:,:), pointer :: edgesOnCell
   integer, dimension(:,:), pointer :: cellsOnEdge
   integer, dimension(:),   pointer :: nEdgesOnCell

   real(r8), dimension(:,:),     pointer :: utangential  ! velocity tangent to cell edge,
                                                         ! diagnosed by mpas

   !
   ! State that may be directly derived from dycore prognostic state
   !
   real(r8), dimension(:,:),   pointer :: theta   ! Potential temperature [K]        (nver,ncol)
   real(r8), dimension(:,:),   pointer :: exner   ! Exner function [-]               (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ux      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), dimension(:,:),   pointer :: uy      ! Meridional veloc at center [m/s] (nver,ncol)
   real(r8), dimension(:,:),   pointer :: pmiddry ! Dry hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nver,ncol)
   real(r8), dimension(:,:),   pointer :: pintdry ! Dry hydrostatic pressure [Pa]
                                                  ! at layer interfaces            (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: vorticity   ! Relative vertical vorticity [s^-1]
                                                      !                              (nver,nvtx)
   real(r8), dimension(:,:),   pointer :: divergence  ! Horizontal velocity divergence [s^-1]
                                                      !                              (nver,ncol)
end type dyn_export_t

! Frontogenesis indices
integer, public    :: frontgf_idx      = -1
integer, public    :: frontga_idx      = -1

real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8

! The global cell indices are used to seed the RNG which is used to apply
! random perturbations to the initial temperature field.  These global indices
! are just those for the local dynamics block.
integer, allocatable :: glob_ind(:)

type (MPAS_TimeInterval_type) :: integrationLength ! set to CAM's dynamics/physics coupling interval

!=========================================================================================
contains
!=========================================================================================

subroutine dyn_readnl(NLFileName)

   ! Read the dycore-relevant namelists from the input file.
   ! First must set up basic MPAS infrastructure to allow the MPAS-A dycore
   ! to save namelist options into MPAS-native datastructures called "pools".

   use units,         only: getunit
   use cam_pio_utils, only: pio_subsystem

   use cam_mpas_subdriver, only : domain_ptr, cam_mpas_init_phase1, cam_mpas_init_phase2
   use mpas_pool_routines, only : mpas_pool_add_config


   ! Dummy argument
   character(len=*), intent(in) :: NLFileName

   ! Local variables
   integer, dimension(2) :: logUnits   ! stdout and stderr for MPAS logging
   integer :: yr, mon, day, tod, ndate, nday, nsec
   character(len=*), parameter :: subname = 'dyn_comp:dyn_readnl'
   !----------------------------------------------------------------------------

   logUnits(1) = iulog
   logUnits(2) = getunit()

   call cam_mpas_init_phase1(mpicom, endrun, logUnits, r8)

   ! read namelist
   call cam_mpas_namelist_read(NLFileName, domain_ptr % configs)

   ! Set config_start_date, etc. (these will not appear in the dycore namelist)
   call get_start_date(yr, mon, day, tod)
   ndate = yr*10000 + mon*100 + day
   call mpas_pool_add_config(domain_ptr % configs, 'config_start_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

   call get_stop_date(yr, mon, day, tod)
   ndate = yr*10000 + mon*100 + day
   call mpas_pool_add_config(domain_ptr % configs, 'config_stop_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

   call get_run_duration(nday, nsec)
   call mpas_pool_add_config(domain_ptr % configs, 'config_run_duration', trim(int2str(nday))//'_'//sec2hms(nsec))

   ! Although the following namelist options are not expected to be used by CAM-MPAS, the MPAS-A dycore
   ! references these options, and they therefore must be defined in the configs pool
   call mpas_pool_add_config(domain_ptr % configs, 'config_restart_timestamp_name', 'restart_timestamp')
   call mpas_pool_add_config(domain_ptr % configs, 'config_IAU_option', 'off')
   call mpas_pool_add_config(domain_ptr % configs, 'config_do_DAcycling', .false.)
   call mpas_pool_add_config(domain_ptr % configs, 'config_halo_exch_method', 'mpas_halo')

   call cam_mpas_init_phase2(pio_subsystem, endrun, timemgr_get_calendar_cf())

end subroutine dyn_readnl

!=========================================================================================

subroutine dyn_register()

   ! Register fields that are computed by the dycore and passed to the physics via the
   ! physics buffer.

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver
   use phys_control,    only: use_gw_front, use_gw_front_igw
   !----------------------------------------------------------------------------

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

   if (use_gw_front .or. use_gw_front_igw) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), frontga_idx)
   end if

end subroutine dyn_register

!=========================================================================================

subroutine dyn_init(dyn_in, dyn_out)
   use air_composition,    only : thermodynamic_active_species_idx, thermodynamic_active_species_idx_dycore
   use air_composition,    only : thermodynamic_active_species_num
   use air_composition,    only : thermodynamic_active_species_liq_idx,thermodynamic_active_species_ice_idx
   use air_composition,    only : thermodynamic_active_species_liq_idx_dycore,thermodynamic_active_species_ice_idx_dycore
   use air_composition,    only : thermodynamic_active_species_liq_num, thermodynamic_active_species_ice_num
   use cam_mpas_subdriver, only : domain_ptr, cam_mpas_init_phase4
   use cam_mpas_subdriver, only : cam_mpas_define_scalars
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension, &
                                  mpas_pool_get_config
   use mpas_timekeeping,   only : MPAS_set_timeInterval
   use mpas_derived_types, only : mpas_pool_type
   use mpas_constants,     only : mpas_constants_compute_derived
   use dyn_tests_utils,    only : vc_dycore, vc_height, string_vc, vc_str_lgth
   use cam_budget,         only : thermo_budget_history

   ! arguments:
   type(dyn_import_t), intent(inout)  :: dyn_in
   type(dyn_export_t), intent(inout)  :: dyn_out

   ! Local variables:
   integer :: ierr

   type(mpas_pool_type), pointer :: mesh_pool
   type(mpas_pool_type), pointer :: state_pool
   type(mpas_pool_type), pointer :: diag_pool
   type(mpas_pool_type), pointer :: tend_physics_pool

   integer, pointer :: nCells
   integer, pointer :: nEdges
   integer, pointer :: nVertices
   integer, pointer :: nVertLevels
   integer, pointer :: nCellsSolve
   integer, pointer :: nEdgesSolve
   integer, pointer :: nVerticesSolve
   integer, pointer :: index_qv

   integer, pointer :: indexToCellID(:) ! global indices of cell centers of local block

   real(r8) :: dtime
   real(r8), pointer :: mpas_dt
   real(r8) :: dt_ratio
   character(len=128) :: errmsg

   character(len=*), parameter :: subname = 'dyn_comp::dyn_init'

   ! variables for initializing energy and axial angular momentum diagnostics
   integer, parameter                         :: num_stages = 6
   character (len = 8), dimension(num_stages) :: stage = (/"dBF     ","dAP     ","dAM     ","BD_dparm","BD_DMEA ","BD_phys "/)
   character (len = 55),dimension(num_stages) :: stage_txt = (/&
      " dynamics state before physics (d_p_coupling)       ",&
      " dynamics state with T,u,V increment but not q      ",&
      " dynamics state with full physics increment (incl.q)",&
      "dE/dt params+efix in dycore (dparam)(dAP-dBF)       ",&
      "dE/dt dry mass adjustment in dycore        (dAM-dAP)",&
      "dE/dt physics total in dycore (phys)       (dAM-dBF)" &
      /)

   integer :: istage, ivars, m
   character (len=108)         :: str1, str2, str3
   character (len=vc_str_lgth) :: vc_str
   !-------------------------------------------------------

   vc_dycore = vc_height
   if (masterproc) then
     call string_vc(vc_dycore,vc_str)
     write(iulog,*)'vertical coordinate dycore   : ',trim(vc_str)
   end if
   !----------------------------------------------------------------------------

   if (initial_run) then
      call cam_mpas_define_scalars(domain_ptr % blocklist, dyn_in % mpas_from_cam_cnst, &
                                   dyn_out % cam_from_mpas_cnst, ierr)
      if (ierr /= 0) then
         call endrun(subname//': Set-up of constituents for MPAS-A dycore failed.')
      end if
   end if

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics_pool)

   ! Let dynamics import state point to memory managed by MPAS-Atmosphere

   call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells)
   dyn_in % nCells = nCells

   call mpas_pool_get_dimension(mesh_pool, 'nEdges', nEdges)
   dyn_in % nEdges = nEdges

   call mpas_pool_get_dimension(mesh_pool, 'nVertices', nVertices)
   dyn_in % nVertices = nVertices

   call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels)
   dyn_in % nVertLevels = nVertLevels

   call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve)
   dyn_in % nCellsSolve = nCellsSolve

   call mpas_pool_get_dimension(mesh_pool, 'nEdgesSolve', nEdgesSolve)
   dyn_in % nEdgesSolve = nEdgesSolve

   call mpas_pool_get_dimension(mesh_pool, 'nVerticesSolve', nVerticesSolve)
   dyn_in % nVerticesSolve = nVerticesSolve

   ! In MPAS timeLevel=1 is the current state.  So the fields input to the dycore should
   ! be in timeLevel=1.

   call mpas_pool_get_array(state_pool, 'u',                      dyn_in % uperp,   timeLevel=1)
   call mpas_pool_get_array(state_pool, 'w',                      dyn_in % w,       timeLevel=1)
   call mpas_pool_get_array(state_pool, 'theta_m',                dyn_in % theta_m, timeLevel=1)
   call mpas_pool_get_array(state_pool, 'rho_zz',                 dyn_in % rho_zz,  timeLevel=1)
   call mpas_pool_get_array(state_pool, 'scalars',                dyn_in % tracers, timeLevel=1)

   call mpas_pool_get_array(diag_pool, 'rho_base',                dyn_in % rho_base)
   call mpas_pool_get_array(diag_pool, 'theta_base',              dyn_in % theta_base)

   call mpas_pool_get_dimension(state_pool, 'index_qv', index_qv)
   dyn_in % index_qv = index_qv

   call mpas_pool_get_array(mesh_pool,  'zgrid',                  dyn_in % zint)
   call mpas_pool_get_array(mesh_pool,  'zz',                     dyn_in % zz)
   call mpas_pool_get_array(mesh_pool,  'fzm',                    dyn_in % fzm)
   call mpas_pool_get_array(mesh_pool,  'fzp',                    dyn_in % fzp)
   call mpas_pool_get_array(mesh_pool,  'areaCell',               dyn_in % areaCell)

   call mpas_pool_get_array(mesh_pool,  'east',                   dyn_in % east)
   call mpas_pool_get_array(mesh_pool,  'north',                  dyn_in % north)
   call mpas_pool_get_array(mesh_pool,  'edgeNormalVectors',      dyn_in % normal)
   call mpas_pool_get_array(mesh_pool,  'cellsOnEdge',            dyn_in % cellsOnEdge)

   call mpas_pool_get_array(diag_pool,  'theta',                  dyn_in % theta)
   call mpas_pool_get_array(diag_pool,  'exner',                  dyn_in % exner)
   call mpas_pool_get_array(diag_pool,  'rho',                    dyn_in % rho)
   call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      dyn_in % ux)
   call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', dyn_in % uy)

   ! Let dynamics export state point to memory managed by MPAS-Atmosphere
   ! Exception: pmiddry and pintdry are not managed by the MPAS infrastructure

   dyn_out % nCells         = dyn_in % nCells
   dyn_out % nEdges         = dyn_in % nEdges
   dyn_out % nVertices      = dyn_in % nVertices
   dyn_out % nVertLevels    = dyn_in % nVertLevels
   dyn_out % nCellsSolve    = dyn_in % nCellsSolve
   dyn_out % nEdgesSolve    = dyn_in % nEdgesSolve
   dyn_out % nVerticesSolve = dyn_in % nVerticesSolve
   dyn_out % index_qv       = dyn_in % index_qv

   ! MPAS swaps pointers internally so that after a dycore timestep, the updated state is
   ! in timeLevel=1.  Thus we want dyn_out to also point to timeLevel=1.  Can just copy
   ! the pointers from dyn_in.

   dyn_out % uperp   => dyn_in % uperp
   dyn_out % w       => dyn_in % w
   dyn_out % theta_m => dyn_in % theta_m
   dyn_out % rho_zz  => dyn_in % rho_zz
   dyn_out % tracers => dyn_in % tracers

   ! These components don't have a time level index.
   dyn_out % zint  => dyn_in % zint
   dyn_out % zz    => dyn_in % zz
   dyn_out % fzm   => dyn_in % fzm
   dyn_out % fzp   => dyn_in % fzp

   dyn_out % theta => dyn_in % theta
   dyn_out % exner => dyn_in % exner
   dyn_out % rho   => dyn_in % rho
   dyn_out % ux    => dyn_in % ux
   dyn_out % uy    => dyn_in % uy

   ! for frontogenesis calc

   if (use_gw_front .or. use_gw_front_igw) then
      dyn_out % areaCell => dyn_in % areaCell
      dyn_out % cellsOnEdge => dyn_in % cellsOnEdge
      call mpas_pool_get_array(mesh_pool, 'defc_a',               dyn_out % defc_a)
      call mpas_pool_get_array(mesh_pool, 'defc_b',               dyn_out % defc_b)
      call mpas_pool_get_array(mesh_pool, 'cell_gradient_coef_x', dyn_out % cell_gradient_coef_x)
      call mpas_pool_get_array(mesh_pool, 'cell_gradient_coef_y', dyn_out % cell_gradient_coef_y)
      call mpas_pool_get_array(mesh_pool, 'edgesOnCell_sign',     dyn_out % edgesOnCell_sign)
      call mpas_pool_get_array(mesh_pool, 'dvEdge',               dyn_out % dvEdge)
      call mpas_pool_get_array(mesh_pool, 'edgesOnCell',          dyn_out % edgesOnCell)
      call mpas_pool_get_array(mesh_pool, 'nEdgesOnCell',         dyn_out % nEdgesOnCell)
      call mpas_pool_get_array(diag_pool, 'v',                    dyn_out % utangential)
   endif

   ! cam-required hydrostatic pressures

   allocate(dyn_out % pmiddry(nVertLevels,   nCells), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate dyn_out%pmiddry array')

   allocate(dyn_out % pintdry(nVertLevels+1, nCells), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate dyn_out%pintdry array')

   call mpas_pool_get_array(diag_pool, 'vorticity',  dyn_out % vorticity)
   call mpas_pool_get_array(diag_pool, 'divergence', dyn_out % divergence)

   call mpas_pool_get_array(mesh_pool, 'indexToCellID', indexToCellID)
   allocate(glob_ind(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate glob_ind array')

   glob_ind = indexToCellID(1:nCellsSolve)

   call mpas_constants_compute_derived()

   if (initial_run) then

      call read_inidat(dyn_in)
      call clean_iodesc_list()

   end if

   call cam_mpas_init_phase4(endrun)

   !
   ! Set pointers to tendency fields that are not allocated until the call to cam_mpas_init_phase4
   !
   call mpas_pool_get_array(tend_physics_pool, 'tend_ru_physics',     dyn_in % ru_tend)
   call mpas_pool_get_array(tend_physics_pool, 'tend_rtheta_physics', dyn_in % rtheta_tend)
   call mpas_pool_get_array(tend_physics_pool, 'tend_rho_physics',    dyn_in % rho_tend)

   ! Check that CAM's timestep, i.e., the dynamics/physics coupling interval, is an integer multiple
   ! of the MPAS timestep.

   ! Get CAM time step
   dtime = get_step_size()

   ! Get MPAS-A dycore time step
   call mpas_pool_get_config(domain_ptr % configs, 'config_dt', mpas_dt)

   ! Calculate time step ratio
   dt_ratio = dtime / mpas_dt

   ! Stop if the dycore time step does not evenly divide the CAM time step
   if (ceiling(dt_ratio) /= floor(dt_ratio)) then
      write(errmsg, '(a,f9.3,a,f9.3,a)') 'The ratio of the CAM timestep, ', dtime, &
                                         ' to the MPAS-A dycore timestep, ', mpas_dt, ' is not an integer'
      call endrun(subname//': '//trim(errmsg))
   end if

   ! dtime has no fractional part, but use nint to deal with any roundoff errors.
   ! Set the interval over which the dycore should integrate during each call to dyn_run.
   call MPAS_set_timeInterval(integrationLength, S=nint(dtime), S_n=0, S_d=1)

   !
   ! initialize history for MPAS energy budgets

   if (thermo_budget_history) then

      ! Define energy/mass snapshots using stage structure
      do istage = 1, num_stages
         call cam_budget_em_snapshot(TRIM(ADJUSTL(stage(istage))), 'dyn', longname=TRIM(ADJUSTL(stage_txt(istage))))
      end do
      !
      ! initialize MPAS energy budgets
      ! add budgets that are derived from stages
      !
      call cam_budget_em_register('dEdt_param_efix_in_dyn','dAP','dBF',pkgtype='dyn',optype='dif', &
                      longname="dE/dt parameterizations+efix in dycore (dparam)(dAP-dBF)")
      call cam_budget_em_register('dEdt_dme_adjust_in_dyn','dAM','dAP',pkgtype='dyn',optype='dif', &
                      longname="dE/dt dry mass adjustment in dycore (dAM-dAP)")
      call cam_budget_em_register('dEdt_phys_total_in_dyn','dAM','dBF',pkgtype='dyn',optype='dif', &
                      longname="dE/dt physics total in dycore (phys) (dAM-dBF)")
   end if

   !
   ! initialize CAM thermodynamic infrastructure
   !
   do m=1,thermodynamic_active_species_num
      thermodynamic_active_species_idx_dycore(m) = dyn_out % cam_from_mpas_cnst(thermodynamic_active_species_idx(m))
      if (masterproc) then
         write(iulog,'(a,2I4)') subname//": m,thermodynamic_active_species_idx_dycore: ", &
                                m,thermodynamic_active_species_idx_dycore(m)
      end if
   end do
   do m=1,thermodynamic_active_species_liq_num
      thermodynamic_active_species_liq_idx_dycore(m) = dyn_out % cam_from_mpas_cnst(thermodynamic_active_species_liq_idx(m))
      if (masterproc) then
         write(iulog,'(a,2I4)') subname//": m,thermodynamic_active_species_idx_liq_dycore: ", &
                                m,thermodynamic_active_species_liq_idx_dycore(m)
      end if
   end do
   do m=1,thermodynamic_active_species_ice_num
      thermodynamic_active_species_ice_idx_dycore(m) = dyn_out % cam_from_mpas_cnst(thermodynamic_active_species_ice_idx(m))
      if (masterproc) then
         write(iulog,'(a,2I4)') subname//": m,thermodynamic_active_species_idx_ice_dycore: ", &
                                m,thermodynamic_active_species_ice_idx_dycore(m)
      end if
   end do

 end subroutine dyn_init

!=========================================================================================

subroutine dyn_run(dyn_in, dyn_out)

   use cam_mpas_subdriver, only : cam_mpas_run
   use cam_mpas_subdriver, only : domain_ptr
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type

   ! Advances the dynamics state provided in dyn_in by one physics
   ! timestep to produce dynamics state held in dyn_out.

   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   ! Local variables
   type(mpas_pool_type), pointer :: state_pool
   character(len=*), parameter :: subname = 'dyn_comp:dyn_run'
   real(r8) :: dtime

   !----------------------------------------------------------------------------

   ! Call the MPAS-A dycore
   call cam_mpas_run(integrationLength)

   ! Update the dyn_in/dyn_out pointers to the current state of the prognostic fields.
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
   call mpas_pool_get_array(state_pool, 'u',       dyn_in % uperp,   timeLevel=1)
   call mpas_pool_get_array(state_pool, 'w',       dyn_in % w,       timeLevel=1)
   call mpas_pool_get_array(state_pool, 'theta_m', dyn_in % theta_m, timeLevel=1)
   call mpas_pool_get_array(state_pool, 'rho_zz',  dyn_in % rho_zz,  timeLevel=1)
   call mpas_pool_get_array(state_pool, 'scalars', dyn_in % tracers, timeLevel=1)
   dyn_out % uperp   => dyn_in % uperp
   dyn_out % w       => dyn_in % w
   dyn_out % theta_m => dyn_in % theta_m
   dyn_out % rho_zz  => dyn_in % rho_zz
   dyn_out % tracers => dyn_in % tracers

end subroutine dyn_run


subroutine dyn_final(dyn_in, dyn_out)

  use cam_mpas_subdriver, only : cam_mpas_finalize

   ! Deallocates the dynamics import and export states, and finalizes
   ! the MPAS dycore.

   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out
   character(len=*), parameter :: subname = 'dyn_comp:dyn_final'
   !----------------------------------------------------------------------------

   !
   ! Prevent any further access to MPAS-Atmosphere memory
   !
   dyn_in % nCells = 0
   dyn_in % nEdges = 0
   dyn_in % nVertices = 0
   dyn_in % nVertLevels = 0
   dyn_in % nCellsSolve = 0
   dyn_in % nEdgesSolve = 0
   dyn_in % nVerticesSolve = 0
   nullify(dyn_in % uperp)
   nullify(dyn_in % w)
   nullify(dyn_in % theta_m)
   nullify(dyn_in % rho_zz)
   nullify(dyn_in % tracers)
   deallocate(dyn_in % mpas_from_cam_cnst)
   nullify(dyn_in % rho_base)
   nullify(dyn_in % theta_base)
   dyn_in % index_qv = 0
   nullify(dyn_in % zint)
   nullify(dyn_in % zz)
   nullify(dyn_in % fzm)
   nullify(dyn_in % fzp)
   nullify(dyn_in % east)
   nullify(dyn_in % north)
   nullify(dyn_in % normal)
   nullify(dyn_in % cellsOnEdge)
   nullify(dyn_in % theta)
   nullify(dyn_in % exner)
   nullify(dyn_in % rho)
   nullify(dyn_in % ux)
   nullify(dyn_in % uy)
   nullify(dyn_in % ru_tend)
   nullify(dyn_in % rtheta_tend)
   nullify(dyn_in % rho_tend)

   !
   ! Prevent any further access to MPAS-Atmosphere memory
   !
   dyn_out % nCells = 0
   dyn_out % nEdges = 0
   dyn_out % nVertices = 0
   dyn_out % nVertLevels = 0
   dyn_out % nCellsSolve = 0
   dyn_out % nEdgesSolve = 0
   dyn_out % nVerticesSolve = 0
   nullify(dyn_out % uperp)
   nullify(dyn_out % w)
   nullify(dyn_out % theta_m)
   nullify(dyn_out % rho_zz)
   nullify(dyn_out % tracers)
   deallocate(dyn_out % cam_from_mpas_cnst)
   dyn_out % index_qv = 0
   nullify(dyn_out % zint)
   nullify(dyn_out % zz)
   nullify(dyn_out % fzm)
   nullify(dyn_out % fzp)
   nullify(dyn_out % theta)
   nullify(dyn_out % exner)
   nullify(dyn_out % rho)
   nullify(dyn_out % ux)
   nullify(dyn_out % uy)
   deallocate(dyn_out % pmiddry)
   deallocate(dyn_out % pintdry)
   nullify(dyn_out % vorticity)
   nullify(dyn_out % divergence)

   call cam_mpas_finalize()

end subroutine dyn_final

!=========================================================================================
! Private routines.
!=========================================================================================

subroutine read_inidat(dyn_in)

   ! Set initial conditions.  Either from analytic expressions or read from file.

   use cam_mpas_subdriver, only : domain_ptr, cam_mpas_update_halo, cam_mpas_cell_to_edge_winds
   use cam_initfiles, only : scale_dry_air_mass
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type
   use mpas_vector_reconstruction, only : mpas_reconstruct
   use mpas_constants, only : Rv_over_Rd => rvord
   use mpas_constants, only : rgas
   use mpas_constants, only : p0
   use mpas_constants, only : gravity
   use string_utils,   only : int2str
   use shr_kind_mod,   only : shr_kind_cx
   ! arguments
   type(dyn_import_t), target, intent(inout) :: dyn_in

   ! Local variables
   integer :: nCellsSolve, nEdgesSolve
   integer :: i, k, kk, m

   type(file_desc_t), pointer :: fh_ini
   type(file_desc_t), pointer :: fh_topo

   real(r8), allocatable :: latvals(:)
   real(r8), allocatable :: lonvals(:)
   real(r8), pointer     :: latvals_deg(:)
   real(r8), pointer     :: lonvals_deg(:)

   real(r8), pointer :: uperp(:,:)   ! Normal velocity at edges [m/s]  (nver,nedge)
   real(r8), pointer :: w(:,:)       ! Vertical velocity [m/s]        (nver+1,ncol)
   real(r8), pointer :: theta_m(:,:) ! Moist potential temperature [K]  (nver,ncol)
   real(r8), pointer :: rho_zz(:,:)  ! Dry density [kg/m^3]
                                     ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), pointer :: tracers(:,:,:) ! Tracers [kg/kg dry air]       (nq,nver,ncol)
   real(r8), pointer :: zint(:,:)    ! Geometric height [m]
                                     ! at layer interfaces            (nver+1,ncol)
   real(r8), pointer :: zz(:,:)      ! Vertical coordinate metric [dimensionless]
                                     ! at layer midpoints               (nver,ncol)
   real(r8), pointer :: theta(:,:)   ! Potential temperature [K]        (nver,ncol)
   real(r8), pointer :: rho(:,:)     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), pointer :: ux(:,:)      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), pointer :: uy(:,:)      ! Meridional veloc at center [m/s] (nver,ncol)
   real(r8), pointer :: theta_base(:,:)
   real(r8), pointer :: rho_base(:,:)

   integer :: ixqv
   integer, dimension(:), pointer :: mpas_from_cam_cnst

   integer,  allocatable :: m_ind(:)
   real(r8), allocatable :: &
      cam2d(:), cam3d(:,:), cam4d(:,:,:), zi(:,:) ! temp arrays using CAM data order
   real(r8), allocatable :: zsurf(:)

   ! temp arrays using MPAS data order
   real(r8), allocatable :: t(:,:)       ! temperature
   real(r8), allocatable :: pintdry(:,:) ! dry interface pressures
   real(r8), allocatable :: pmiddry(:,:) ! dry midpoint pressures
   real(r8), allocatable :: pmid(:,:)    ! midpoint pressures
   real(r8), allocatable :: mpas3d(:,:,:)

   real(r8), allocatable :: qv(:), tm(:)

   logical  :: readvar

   character(len=shr_kind_cx) :: str

   type(mpas_pool_type), pointer :: mesh_pool
   type(mpas_pool_type), pointer :: diag_pool

   real(r8), pointer :: uReconstructX(:,:)
   real(r8), pointer :: uReconstructY(:,:)
   real(r8), pointer :: uReconstructZ(:,:)

   integer :: mpas_idx, cam_idx, ierr
   character(len=32) :: trac_name

   character(len=*), parameter :: subname = 'dyn_comp:read_inidat'
   !--------------------------------------------------------------------------------------

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

   nCellsSolve = dyn_in % nCellsSolve
   nEdgesSolve = dyn_in % nEdgesSolve

   ixqv        = dyn_in % index_qv
   mpas_from_cam_cnst => dyn_in % mpas_from_cam_cnst

   uperp      => dyn_in % uperp
   w          => dyn_in % w
   theta_m    => dyn_in % theta_m
   rho_zz     => dyn_in % rho_zz
   tracers    => dyn_in % tracers

   zint       => dyn_in % zint
   zz         => dyn_in % zz
   theta      => dyn_in % theta
   rho        => dyn_in % rho
   ux         => dyn_in % ux
   uy         => dyn_in % uy
   rho_base   => dyn_in % rho_base
   theta_base => dyn_in % theta_base

   ! Check that number of advected tracers is consistent with MPAS.
   if (pcnst /= size(tracers, 1)) then
      write(iulog,*) subname//': number of tracers, pcnst:', size(tracers,1), pcnst
      call endrun(subname//': number of tracers /= pcnst')
   end if

   ! lat/lon needed in radians
   latvals_deg => cam_grid_get_latvals(cam_grid_id('mpas_cell'))
   lonvals_deg => cam_grid_get_lonvals(cam_grid_id('mpas_cell'))
   allocate(latvals(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate latvals array')

   allocate(lonvals(nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate lonvals array')
   latvals(:) = latvals_deg(:)*deg2rad
   lonvals(:) = lonvals_deg(:)*deg2rad

   ! Set ICs.  Either from analytic expressions or read from file.

   allocate( &
      ! temporary arrays using CAM indexing
      cam2d(nCellsSolve),              &
      cam3d(nCellsSolve,plev),         &
      cam4d(nCellsSolve,plev,pcnst),   &
      zi(nCellsSolve,plevp),           &
      ! temporary arrays using MPAS indexing
      t(plev,nCellsSolve),             &
      pintdry(plevp,nCellsSolve),      &
      pmiddry(plev,nCellsSolve),       &
      pmid(plev,nCellsSolve), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate tmp arrays using CAM and MPAS indexing')

   do k = 1, plevp
      kk = plevp - k + 1
      zi(:,kk) = zint(k,:nCellsSolve)
   end do

   ! If using a topo file check that PHIS is consistent with the surface z coordinate.
   if (associated(fh_topo)) then

      allocate(zsurf(nCellsSolve), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate zsurf array')

      call get_zsurf_from_topo(fh_topo, zsurf)

      do i = 1, nCellsSolve
        if (abs(zi(i,plevp) - zsurf(i)) > 0.001_r8) then
          write(str,*) 'zi= ', zi(i,plevp), ' zsurf= ', zsurf(i),' i= ',i
          write(iulog,*) subname//': ERROR: '//TRIM(str)
          call endrun(subname//': ERROR: PHIS not consistent with surface z coordinate; '//TRIM(str))
        end if
      end do

      deallocate(zsurf)
   else
      do i = 1, nCellsSolve
         if (abs(zi(i,plevp)) > 1.0E-12_r8) then
           write(str,*) 'zi= ', zi(i,plevp), ' but PHIS should be zero'
           write(iulog,*) subname//': ERROR: '//TRIM(str)
           call endrun(subname//': ERROR: PHIS not consistent with surface z coordinate; '//TRIM(str))
         end if
      end do
   end if

   if (analytic_ic_active()) then

      w(:,1:nCellsSolve) = 0.0_r8

      ! U, V cell center velocity components

      call dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, zint=zi, U=cam3d)
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            ux(kk,i) = cam3d(i,k)
         end do
      end do

      call dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, zint=zi, V=cam3d)
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            uy(kk,i) = cam3d(i,k)
         end do
      end do

      ! Compute uperp by projecting ux and uy from cell centers to edges
      call cam_mpas_update_halo('uReconstructZonal', endrun)       ! ux => uReconstructZonal
      call cam_mpas_update_halo('uReconstructMeridional', endrun)  ! uy => uReconstructMeridional
      call cam_mpas_cell_to_edge_winds(dyn_in % nEdges, ux, uy, dyn_in % east, dyn_in % north, &
                                       dyn_in % normal, dyn_in % cellsOnEdge, uperp)

      call cam_mpas_update_halo('u', endrun)         ! u is the name of uperp in the MPAS state pool

      ! Constituents

      allocate(m_ind(pcnst), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate m_ind array')
      do m = 1, pcnst
         m_ind(m) = m
      end do
      call dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, zint=zi, m_cnst=m_ind, Q=cam4d)
      do m = 1, pcnst  ! index into MPAS tracers array
         do k = 1, plev
            kk = plev - k + 1
            do i = 1, nCellsSolve
               tracers(m,kk,i) = cam4d(i,k,mpas_from_cam_cnst(m))
            end do
         end do
      end do
      deallocate(m_ind)

      ! Temperature

      call dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, zint=zi, T=cam3d)
      do k = 1, plev
         kk = plev - k + 1
         do i = 1, nCellsSolve
            t(kk,i) = cam3d(i,k)
         end do
      end do

      ! Pressures are needed to convert temperature to potential temperature.

      call dyn_set_inic_col(vcoord, latvals, lonvals, glob_ind, zint=zi, PS=cam2d)
      do i = 1, nCellsSolve
         pintdry(1,i) = cam2d(i)
      end do

      allocate(qv(plev), tm(plev), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate qv and tm arrays')

      do i = 1, nCellsSolve
         do k = 1, plev
            ! convert specific humidity to mixing ratio relative to dry air
            tracers(1,k,i) = tracers(1,k,i)/(1.0_r8 - tracers(1,k,i))
            qv(k) = tracers(1,k,i)
            ! convert temperature to tm (we are using dry air density and Rd in state eqn)
            tm(k) = t(k,i)*(1.0_r8+(Rv_over_Rd)*qv(k))
         end do

         ! integrate up from surface to first mid level.  This is full mid-level pressure (i.e. accounts for vapor).
         ! we are assuming that pintdry(1,i) is the full surface pressure here.
         pmid(1,i) = pintdry(1,i)/(1.0_r8+0.5_r8*(zint(2,i)-zint(1,i))*(1.0_r8+qv(1))*gravity/(rgas*tm(1)))

         ! integrate up the column
         do k=2,plev
            !  this is full mid-level pressure (i.e. accounts for vapor)
            pmid(k,i) = pmid(k-1,i)*(1.0_r8-0.5_r8*(zint(k  ,i)-zint(k-1,i))*gravity*(1.0_r8+qv(k-1))/(rgas*tm(k-1)))/ &
                                    (1.0_r8+0.5_r8*(zint(k+1,i)-zint(k  ,i))*gravity*(1.0_r8+qv(k  ))/(rgas*tm(k  )))
         end do

         do k=1,plev
            ! Note: this is theta and not theta_m
            theta(k,i) = tm(k) * ((p0/pmid(k,i))**(rair/cpair))/(1.0_r8+(Rv_over_Rd)*qv(k))
            ! Dry air density
            rho(k,i) = pmid(k,i) / (rgas * tm(k))
         end do

      end do

      deallocate(qv)
      deallocate(tm)

      rho_zz(:,1:nCellsSolve) = rho(:,1:nCellsSolve) / zz(:,1:nCellsSolve)

      ! Set theta_base and rho_base
      call set_base_state(dyn_in)

   else

      ! read uperp
      allocate( mpas3d(plev,nEdgesSolve,1), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate mpas3d array at line:'//int2str(__LINE__))

      call infld('u', fh_ini, 'lev', 'nEdges', 1, plev, 1, nEdgesSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_edge')
      if (readvar) then
         uperp(:,:nEdgesSolve) = mpas3d(:,:nEdgesSolve,1)
      else
         call endrun(subname//': failed to read u (uperp) from initial file')
      end if
      deallocate( mpas3d )

      call cam_mpas_update_halo('u', endrun)         ! u is the name of uperp in the MPAS state pool

      ! Reconstruct ux and uy from uperp.
      ! This is only needed because during CAM's initialization the physics package
      ! is called before the dycore advances a step.
      nullify(mesh_pool)
      nullify(diag_pool)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh_pool)
      call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag_pool)

      ! The uReconstruct{X,Y,Z} arguments to mpas_reconstruct are required, but these
      ! field already exist in the diag pool
      nullify(uReconstructX)
      nullify(uReconstructY)
      nullify(uReconstructZ)
      call mpas_pool_get_array(diag_pool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(diag_pool, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(diag_pool, 'uReconstructZ', uReconstructZ)

      call mpas_reconstruct(mesh_pool, uperp, &
         uReconstructX, uReconstructY, uReconstructZ, &
         ux, uy)

      ! read w
      allocate( mpas3d(plevp,nCellsSolve,1), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate mpas3d array at line:'//int2str(__LINE__))
      call infld('w', fh_ini, 'ilev', 'nCells', 1, plevp, 1, nCellsSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_cell')
      if (readvar) then
         w(:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
      else
         call endrun(subname//': failed to read w from initial file')
      end if
      deallocate( mpas3d )

      allocate( mpas3d(plev,nCellsSolve,1), stat=ierr)
      if( ierr /= 0 ) call endrun(subname//': failed to allocate mpas3d array at line:'//int2str(__LINE__))

      ! read theta
      call infld('theta', fh_ini, 'lev', 'nCells', 1, plev, 1, nCellsSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_cell')
      if (readvar) then
         theta(:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
      else
         call endrun(subname//': failed to read theta from initial file')
      end if

      ! read rho
      call infld('rho', fh_ini, 'lev', 'nCells', 1, plev, 1, nCellsSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_cell')
      if (readvar) then
         rho(:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
      else
         call endrun(subname//': failed to read rho from initial file')
      end if

      rho_zz(:,1:nCellsSolve) = rho(:,1:nCellsSolve) / zz(:,1:nCellsSolve)

      ! read theta_base
      call infld('theta_base', fh_ini, 'lev', 'nCells', 1, plev, 1, nCellsSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_cell')
      if (readvar) then
         theta_base(:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
      else
         call endrun(subname//': failed to read theta_base from initial file')
      end if

      ! read rho_base
      call infld('rho_base', fh_ini, 'lev', 'nCells', 1, plev, 1, nCellsSolve, 1, 1, &
                 mpas3d, readvar, gridname='mpas_cell')
      if (readvar) then
         rho_base(:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
      else
         call endrun(subname//': failed to read rho_base from initial file')
      end if

      deallocate( mpas3d )

   end if

   ! Finish initialization of tracer fields.
   !
   ! If analytic ICs are being used, we allow constituents in an initial
   ! file to overwrite mixing ratios set by the default constituent initialization
   ! except for the water species.

   allocate( mpas3d(plev,nCellsSolve,1), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate mpas3d array at line:'//int2str(__LINE__))

   do mpas_idx = 1, pcnst

      ! The names of the species in the MPAS initial file may be different from the
      ! names in the CAM constituent module.  Also the species order in the MPAS
      ! tracers array may be different from the order in the CAM constituent array.

      cam_idx = mpas_from_cam_cnst(mpas_idx)

      if (analytic_ic_active() .and. cnst_is_a_water_species(cnst_name(cam_idx))) cycle

      ! The name translation is hardcoded here temporarily...
      trac_name = cnst_name(cam_idx)
      if (mpas_idx == 1) trac_name = 'qv'


      readvar = .false.
      if (cnst_read_iv(cam_idx)) then

         ! read constituent from the initial file if present
         call infld(trac_name, fh_ini, 'lev', 'nCells', 1, plev, 1, nCellsSolve, 1, 1, &
                    mpas3d, readvar, gridname='mpas_cell')
         if (readvar) then
            tracers(mpas_idx,:,1:nCellsSolve) = mpas3d(:,:nCellsSolve,1)
         end if
      end if
      if (.not. readvar .and. .not. analytic_ic_active()) then
         ! default constituent initialization (this was already done if analytic ICs are active)
         call cnst_init_default(cam_idx, latvals, lonvals, cam3d)
         do k = 1, plev
            kk = plev - k + 1
            do i = 1, nCellsSolve
               tracers(mpas_idx,kk,i) = cam3d(i,k)
            end do
         end do

      end if
   end do

   deallocate( mpas3d )

   theta_m(:,1:nCellsSolve) = theta(:,1:nCellsSolve) * (1.0_r8 + Rv_over_Rd * tracers(ixqv,:,1:nCellsSolve))

   ! If scale_dry_air_mass > 0.0 then scale dry air mass to scale_dry_air_mass global average dry pressure
   if (scale_dry_air_mass > 0.0_r8) then
     call set_dry_mass(dyn_in, scale_dry_air_mass)
   end if


   ! Update halos for initial state fields
   ! halo for 'u' updated in both branches of conditional above
   call cam_mpas_update_halo('w', endrun)
   call cam_mpas_update_halo('scalars', endrun)   ! scalars is the name of tracers in the MPAS state pool
   call cam_mpas_update_halo('theta_m', endrun)
   call cam_mpas_update_halo('theta', endrun)
   call cam_mpas_update_halo('rho_zz', endrun)
   call cam_mpas_update_halo('rho', endrun)
   call cam_mpas_update_halo('rho_base', endrun)
   call cam_mpas_update_halo('theta_base', endrun)

   deallocate(cam2d, cam3d, cam4d, zi, t, pintdry, pmiddry, pmid)

end subroutine read_inidat

!========================================================================================

subroutine get_zsurf_from_topo(fh_topo, zsurf)

   ! Read PHIS from the topo file and convert it to a surface height field.

   ! Arguments
   type(file_desc_t), pointer :: fh_topo

   real(r8), intent(out) :: zsurf(:)

   ! Local variables
   integer :: zsurf_len
   real(r8), allocatable :: phis(:,:)
   logical :: readvar
   integer :: ierr
   character(len=*), parameter :: subname = 'dyn_comp:get_zsurf_from_topo'
   !--------------------------------------------------------------------------------------

   zsurf_len = size(zsurf)
   allocate(phis(zsurf_len,1), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate phis array')

   ! read theta
   call infld('PHIS', fh_topo, 'ncol', 1, zsurf_len, 1, 1, &
                 phis, readvar, gridname='cam_cell')
   if (readvar) then
      zsurf = phis(:,1) / gravit
   else
      call endrun(subname//': failed to read PHIS from topo file')
   end if

end subroutine get_zsurf_from_topo

!========================================================================================

subroutine set_base_state(dyn_in)

   use mpas_constants, only : gravity, cp, Rgas, p0

   ! Set base-state fields for dynamics assuming an isothermal atmosphere

   ! Arguments
   type(dyn_import_t), intent(inout) :: dyn_in

   ! Local variables
   real(r8), parameter :: t0b = 250.0_r8      ! Temperature [K]

   integer :: iCell, klev
   real(r8), dimension(:,:), pointer :: zint
   real(r8), dimension(:,:), pointer :: zz
   real(r8), dimension(:,:), pointer :: rho_base
   real(r8), dimension(:,:), pointer :: theta_base
   real(r8) :: zmid
   real(r8) :: pres
   real(r8) :: pres_kp1
   logical, parameter :: discrete_hydrostatic_base = .true.
   character(len=*), parameter :: subname = 'dyn_comp:get_zsurf_from_topo'
   !--------------------------------------------------------------------------------------

   zint       => dyn_in % zint
   zz         => dyn_in % zz
   rho_base   => dyn_in % rho_base
   theta_base => dyn_in % theta_base

   ! reference state with discrete MPAS hydrostatic balance

   if (discrete_hydrostatic_base) then

      do iCell = 1, dyn_in % nCellsSolve

         klev = dyn_in % nVertLevels
         ! reference pressure at the model top
         pres_kp1 = p0*exp(-gravity*zint(klev+1,iCell)/(Rgas*t0b))

         ! integrate down to first mid level, set referfence state
         pres = pres_kp1/(1.0_r8-0.5_r8*(zint(klev+1,iCell) - zint(klev,iCell))*gravity/(Rgas*t0b))
         theta_base(klev,iCell) = t0b / (pres / p0)**(Rgas/cp)
         rho_base(klev,iCell) = pres / ( Rgas * t0b * zz(klev,iCell))
         pres_kp1 = pres

         ! integrate down the column
         do klev = dyn_in % nVertLevels-1, 1, -1
            pres = pres_kp1*(1.0_r8+0.5_r8*(zint(klev+2,iCell)-zint(klev+1,iCell))*gravity/(rgas*t0b))/  &
                            (1.0_r8-0.5_r8*(zint(klev+1,iCell)-zint(klev  ,iCell))*gravity/(rgas*t0b))
            theta_base(klev,iCell) = t0b / (pres / p0)**(Rgas/cp)
            rho_base(klev,iCell) = pres / ( Rgas * t0b * zz(klev,iCell))
            pres_kp1 = pres
         end do
      end do

   else

      do iCell = 1, dyn_in % nCellsSolve
         do klev = 1, dyn_in % nVertLevels
            zmid = 0.5_r8 * (zint(klev,iCell) + zint(klev+1,iCell))   ! Layer midpoint geometric height
            pres = p0 * exp(-gravity * zmid / (Rgas * t0b))
            theta_base(klev,iCell) = t0b / (pres / p0)**(Rgas/cp)
            rho_base(klev,iCell) = pres / ( Rgas * t0b * zz(klev,iCell))
         end do
      end do

   end if

end subroutine set_base_state

!========================================================================================

subroutine cam_mpas_namelist_read(namelistFilename, configPool)

   ! Read MPAS-A dycore namelists and add the namelists to the MPAS configPool.
   !
   ! Only the CAM masterproc actually opens and reads from the specified file. Upon return,
   ! if no errors were encountered, all MPI ranks have valid namelists in their configPool.

   use spmd_utils,         only: mpicom, masterproc, masterprocid, &
                                 mpi_integer, mpi_real8, mpi_logical, mpi_character
   use namelist_utils,     only: find_group_name

   use mpas_derived_types, only: mpas_pool_type
   use mpas_kind_types,    only: StrKIND
   use mpas_pool_routines, only: mpas_pool_add_config

   ! Arguments
   character(len=*), intent(in) :: namelistFilename
   type (mpas_pool_type), intent(inout) :: configPool

   ! Local variables
   integer :: unitNumber

   integer :: ierr, ierr2
   integer :: mpi_ierr

   character (len=StrKIND) :: mpas_time_integration = 'SRK3'
   integer                 :: mpas_time_integration_order = 2
   real(r8)                :: mpas_dt = 720.0_r8
   logical                 :: mpas_split_dynamics_transport = .true.
   integer                 :: mpas_number_of_sub_steps = 2
   integer                 :: mpas_dynamics_split_steps = 3
   real(r8)                :: mpas_h_mom_eddy_visc2 = 0.0_r8
   real(r8)                :: mpas_h_mom_eddy_visc4 = 0.0_r8
   real(r8)                :: mpas_v_mom_eddy_visc2 = 0.0_r8
   real(r8)                :: mpas_h_theta_eddy_visc2 = 0.0_r8
   real(r8)                :: mpas_h_theta_eddy_visc4 = 0.0_r8
   real(r8)                :: mpas_v_theta_eddy_visc2 = 0.0_r8
   character (len=StrKIND) :: mpas_horiz_mixing = '2d_smagorinsky'
   real(r8)                :: mpas_len_disp = 120000.0_r8
   real(r8)                :: mpas_visc4_2dsmag = 0.05_r8
   real(r8)                :: mpas_del4u_div_factor = 10.0_r8
   integer                 :: mpas_w_adv_order = 3
   integer                 :: mpas_theta_adv_order = 3
   integer                 :: mpas_scalar_adv_order = 3
   integer                 :: mpas_u_vadv_order = 3
   integer                 :: mpas_w_vadv_order = 3
   integer                 :: mpas_theta_vadv_order = 3
   integer                 :: mpas_scalar_vadv_order = 3
   logical                 :: mpas_scalar_advection = .true.
   logical                 :: mpas_positive_definite = .false.
   logical                 :: mpas_monotonic = .true.
   real(r8)                :: mpas_coef_3rd_order = 0.25_r8
   real(r8)                :: mpas_smagorinsky_coef = 0.125_r8
   logical                 :: mpas_mix_full = .true.
   real(r8)                :: mpas_epssm = 0.1_r8
   real(r8)                :: mpas_smdiv = 0.1_r8
   real(r8)                :: mpas_apvm_upwinding = 0.5_r8
   logical                 :: mpas_h_ScaleWithMesh = .true.
   real(r8)                :: mpas_zd = 22000.0_r8
   real(r8)                :: mpas_xnutr = 0.2_r8
   real(r8)                :: mpas_cam_coef = 0.0_r8
   integer                 :: mpas_cam_damping_levels = 0
   logical                 :: mpas_rayleigh_damp_u = .true.
   real(r8)                :: mpas_rayleigh_damp_u_timescale_days = 5.0_r8
   integer                 :: mpas_number_rayleigh_damp_u_levels = 3
   logical                 :: mpas_apply_lbcs = .false.
   logical                 :: mpas_jedi_da = .false.
   character (len=StrKIND) :: mpas_block_decomp_file_prefix = 'x1.40962.graph.info.part.'
   logical                 :: mpas_do_restart = .false.
   logical                 :: mpas_print_global_minmax_vel = .true.
   logical                 :: mpas_print_detailed_minmax_vel = .false.
   logical                 :: mpas_print_global_minmax_sca = .false.

   namelist /nhyd_model/ &
           mpas_time_integration, &
           mpas_time_integration_order, &
           mpas_dt, &
           mpas_split_dynamics_transport, &
           mpas_number_of_sub_steps, &
           mpas_dynamics_split_steps, &
           mpas_h_mom_eddy_visc2, &
           mpas_h_mom_eddy_visc4, &
           mpas_v_mom_eddy_visc2, &
           mpas_h_theta_eddy_visc2, &
           mpas_h_theta_eddy_visc4, &
           mpas_v_theta_eddy_visc2, &
           mpas_horiz_mixing, &
           mpas_len_disp, &
           mpas_visc4_2dsmag, &
           mpas_del4u_div_factor, &
           mpas_w_adv_order, &
           mpas_theta_adv_order, &
           mpas_scalar_adv_order, &
           mpas_u_vadv_order, &
           mpas_w_vadv_order, &
           mpas_theta_vadv_order, &
           mpas_scalar_vadv_order, &
           mpas_scalar_advection, &
           mpas_positive_definite, &
           mpas_monotonic, &
           mpas_coef_3rd_order, &
           mpas_smagorinsky_coef, &
           mpas_mix_full, &
           mpas_epssm, &
           mpas_smdiv, &
           mpas_apvm_upwinding, &
           mpas_h_ScaleWithMesh

   namelist /damping/ &
           mpas_zd, &
           mpas_xnutr, &
           mpas_cam_coef, &
           mpas_cam_damping_levels, &
           mpas_rayleigh_damp_u, &
           mpas_rayleigh_damp_u_timescale_days, &
           mpas_number_rayleigh_damp_u_levels

   namelist /limited_area/ &
           mpas_apply_lbcs

   namelist /assimilation/ &
           mpas_jedi_da

   namelist /decomposition/ &
           mpas_block_decomp_file_prefix

   namelist /restart/ &
           mpas_do_restart

   namelist /printout/ &
           mpas_print_global_minmax_vel, &
           mpas_print_detailed_minmax_vel, &
           mpas_print_global_minmax_sca

   ! These configuration parameters must be set in the MPAS configPool, but can't
   ! be changed in CAM.
   integer                :: config_num_halos = 2
   integer                :: config_number_of_blocks = 0
   logical                :: config_explicit_proc_decomp = .false.
   character(len=StrKIND) :: config_proc_decomp_file_prefix = 'graph.info.part'

   character(len=*), parameter :: subname = 'dyn_comp::cam_mpas_namelist_read'
   !----------------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'Reading MPAS-A dycore namelist from ', trim(namelistFilename)
      open(newunit=unitNumber, file=trim(namelistFilename), status='old', form='formatted')
   end if

   ! Read namelist group &nhyd_model
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'nhyd_model', status=ierr)
      if (ierr == 0) then
         read(unitNumber, nhyd_model, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &nhyd_model')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &nhyd_model')
      end if
   end if

   call mpi_bcast(mpas_time_integration,       StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_time_integration_order,       1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_dt,                           1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_split_dynamics_transport,     1, mpi_logical,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_number_of_sub_steps,          1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_dynamics_split_steps,         1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_h_mom_eddy_visc2,             1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_h_mom_eddy_visc4,             1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_v_mom_eddy_visc2,             1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_h_theta_eddy_visc2,           1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_h_theta_eddy_visc4,           1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_v_theta_eddy_visc2,           1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_horiz_mixing,           StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_len_disp,                     1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_visc4_2dsmag,                 1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_del4u_div_factor,             1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_w_adv_order,                  1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_theta_adv_order,              1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_scalar_adv_order,             1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_u_vadv_order,                 1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_w_vadv_order,                 1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_theta_vadv_order,             1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_scalar_vadv_order,            1, mpi_integer,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_scalar_advection,             1, mpi_logical,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_positive_definite,            1, mpi_logical,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_monotonic,                    1, mpi_logical,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_coef_3rd_order,               1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_smagorinsky_coef,             1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_mix_full,                     1, mpi_logical,   masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_epssm,                        1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_smdiv,                        1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_apvm_upwinding,               1, mpi_real8,     masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_h_ScaleWithMesh,              1, mpi_logical,   masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_time_integration', mpas_time_integration)
   call mpas_pool_add_config(configPool, 'config_time_integration_order', mpas_time_integration_order)
   call mpas_pool_add_config(configPool, 'config_dt', mpas_dt)
   call mpas_pool_add_config(configPool, 'config_split_dynamics_transport', mpas_split_dynamics_transport)
   call mpas_pool_add_config(configPool, 'config_number_of_sub_steps', mpas_number_of_sub_steps)
   call mpas_pool_add_config(configPool, 'config_dynamics_split_steps', mpas_dynamics_split_steps)
   call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc2', mpas_h_mom_eddy_visc2)
   call mpas_pool_add_config(configPool, 'config_h_mom_eddy_visc4', mpas_h_mom_eddy_visc4)
   call mpas_pool_add_config(configPool, 'config_v_mom_eddy_visc2', mpas_v_mom_eddy_visc2)
   call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc2', mpas_h_theta_eddy_visc2)
   call mpas_pool_add_config(configPool, 'config_h_theta_eddy_visc4', mpas_h_theta_eddy_visc4)
   call mpas_pool_add_config(configPool, 'config_v_theta_eddy_visc2', mpas_v_theta_eddy_visc2)
   call mpas_pool_add_config(configPool, 'config_horiz_mixing', mpas_horiz_mixing)
   call mpas_pool_add_config(configPool, 'config_len_disp', mpas_len_disp)
   call mpas_pool_add_config(configPool, 'config_visc4_2dsmag', mpas_visc4_2dsmag)
   call mpas_pool_add_config(configPool, 'config_del4u_div_factor', mpas_del4u_div_factor)
   call mpas_pool_add_config(configPool, 'config_w_adv_order', mpas_w_adv_order)
   call mpas_pool_add_config(configPool, 'config_theta_adv_order', mpas_theta_adv_order)
   call mpas_pool_add_config(configPool, 'config_scalar_adv_order', mpas_scalar_adv_order)
   call mpas_pool_add_config(configPool, 'config_u_vadv_order', mpas_u_vadv_order)
   call mpas_pool_add_config(configPool, 'config_w_vadv_order', mpas_w_vadv_order)
   call mpas_pool_add_config(configPool, 'config_theta_vadv_order', mpas_theta_vadv_order)
   call mpas_pool_add_config(configPool, 'config_scalar_vadv_order', mpas_scalar_vadv_order)
   call mpas_pool_add_config(configPool, 'config_scalar_advection', mpas_scalar_advection)
   call mpas_pool_add_config(configPool, 'config_positive_definite', mpas_positive_definite)
   call mpas_pool_add_config(configPool, 'config_monotonic', mpas_monotonic)
   call mpas_pool_add_config(configPool, 'config_coef_3rd_order', mpas_coef_3rd_order)
   call mpas_pool_add_config(configPool, 'config_smagorinsky_coef', mpas_smagorinsky_coef)
   call mpas_pool_add_config(configPool, 'config_mix_full', mpas_mix_full)
   call mpas_pool_add_config(configPool, 'config_epssm', mpas_epssm)
   call mpas_pool_add_config(configPool, 'config_smdiv', mpas_smdiv)
   call mpas_pool_add_config(configPool, 'config_apvm_upwinding', mpas_apvm_upwinding)
   call mpas_pool_add_config(configPool, 'config_h_ScaleWithMesh', mpas_h_ScaleWithMesh)

   ! Read namelist group &damping
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'damping', status=ierr)
      if (ierr == 0) then
         read(unitNumber, damping, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &damping')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &damping')
      end if
   end if

   call mpi_bcast(mpas_zd,       1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_xnutr,    1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_cam_coef, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_cam_damping_levels,             1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_rayleigh_damp_u,                1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_rayleigh_damp_u_timescale_days, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_number_rayleigh_damp_u_levels,  1, mpi_integer, masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_zd', mpas_zd)
   call mpas_pool_add_config(configPool, 'config_xnutr', mpas_xnutr)
   call mpas_pool_add_config(configPool, 'config_mpas_cam_coef', mpas_cam_coef)
   call mpas_pool_add_config(configPool, 'config_number_cam_damping_levels', mpas_cam_damping_levels)
   call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u', mpas_rayleigh_damp_u)
   call mpas_pool_add_config(configPool, 'config_rayleigh_damp_u_timescale_days', mpas_rayleigh_damp_u_timescale_days)
   call mpas_pool_add_config(configPool, 'config_number_rayleigh_damp_u_levels', mpas_number_rayleigh_damp_u_levels)

   ! Read namelist group &limited_area
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'limited_area', status=ierr)
      if (ierr == 0) then
         read(unitNumber, limited_area, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &limited_area')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &limited_area')
      end if
   end if

   call mpi_bcast(mpas_apply_lbcs, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_apply_lbcs', mpas_apply_lbcs)

   ! Read namelist group &assimilation
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'assimilation', status=ierr)
      if (ierr == 0) then
         read(unitNumber, assimilation, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &assimilation')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &assimilation')
      end if
   end if

   call mpi_bcast(mpas_jedi_da, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_jedi_da', mpas_jedi_da)

   ! Read namelist group &decomposition if npes > 1
   if (masterproc .and. npes > 1) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'decomposition', status=ierr)
      if (ierr == 0) then
         read(unitNumber, decomposition, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &decomposition')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &decomposition. Required for multiprocessor execution.')
      end if
   end if

   call mpi_bcast(mpas_block_decomp_file_prefix, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix', mpas_block_decomp_file_prefix)

   ! Read namelist group &restart
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'restart', status=ierr)
      if (ierr == 0) then
         read(unitNumber, restart, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &restart')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &restart')
      end if
   end if

   call mpi_bcast(mpas_do_restart, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)

   ! Set mpas_do_restart based on information from the driver code.
   if (.not. initial_run) mpas_do_restart = .true.

   call mpas_pool_add_config(configPool, 'config_do_restart', mpas_do_restart)

   ! Read namelist group &printout
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'printout', status=ierr)
      if (ierr == 0) then
         read(unitNumber, printout, iostat=ierr2)
         if (ierr2 /= 0) then
            call endrun(subname // ':: Failed to read namelist group &printout')
         end if
      else
         call endrun(subname // ':: Failed to find namelist group &printout')
      end if
   end if

   call mpi_bcast(mpas_print_global_minmax_vel,   1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_print_detailed_minmax_vel, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   call mpi_bcast(mpas_print_global_minmax_sca,   1, mpi_logical, masterprocid, mpicom, mpi_ierr)

   call mpas_pool_add_config(configPool, 'config_print_global_minmax_vel', mpas_print_global_minmax_vel)
   call mpas_pool_add_config(configPool, 'config_print_detailed_minmax_vel', mpas_print_detailed_minmax_vel)
   call mpas_pool_add_config(configPool, 'config_print_global_minmax_sca', mpas_print_global_minmax_sca)

   if (masterproc) then
      close(unit=unitNumber)
   end if

   ! Set some configuration parameters that cannot be changed by CAM.
   call mpas_pool_add_config(configPool, 'config_num_halos', config_num_halos)
   call mpas_pool_add_config(configPool, 'config_number_of_blocks', config_number_of_blocks)
   call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', config_explicit_proc_decomp)
   call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', config_proc_decomp_file_prefix)


   if (masterproc) then
      write(iulog,*) 'MPAS-A dycore configuration:'
      write(iulog,*) '   mpas_time_integration = ', trim(mpas_time_integration)
      write(iulog,*) '   mpas_time_integration_order = ', mpas_time_integration_order
      write(iulog,*) '   mpas_dt = ', mpas_dt
      write(iulog,*) '   mpas_split_dynamics_transport = ', mpas_split_dynamics_transport
      write(iulog,*) '   mpas_number_of_sub_steps = ', mpas_number_of_sub_steps
      write(iulog,*) '   mpas_dynamics_split_steps = ', mpas_dynamics_split_steps
      write(iulog,*) '   mpas_h_mom_eddy_visc2 = ', mpas_h_mom_eddy_visc2
      write(iulog,*) '   mpas_h_mom_eddy_visc4 = ', mpas_h_mom_eddy_visc4
      write(iulog,*) '   mpas_v_mom_eddy_visc2 = ', mpas_v_mom_eddy_visc2
      write(iulog,*) '   mpas_h_theta_eddy_visc2 = ', mpas_h_theta_eddy_visc2
      write(iulog,*) '   mpas_h_theta_eddy_visc4 = ', mpas_h_theta_eddy_visc4
      write(iulog,*) '   mpas_v_theta_eddy_visc2 = ', mpas_v_theta_eddy_visc2
      write(iulog,*) '   mpas_horiz_mixing = ', trim(mpas_horiz_mixing)
      write(iulog,*) '   mpas_len_disp = ', mpas_len_disp
      write(iulog,*) '   mpas_visc4_2dsmag = ', mpas_visc4_2dsmag
      write(iulog,*) '   mpas_del4u_div_factor = ', mpas_del4u_div_factor
      write(iulog,*) '   mpas_w_adv_order = ', mpas_w_adv_order
      write(iulog,*) '   mpas_theta_adv_order = ', mpas_theta_adv_order
      write(iulog,*) '   mpas_scalar_adv_order = ', mpas_scalar_adv_order
      write(iulog,*) '   mpas_u_vadv_order = ', mpas_u_vadv_order
      write(iulog,*) '   mpas_w_vadv_order = ', mpas_w_vadv_order
      write(iulog,*) '   mpas_theta_vadv_order = ', mpas_theta_vadv_order
      write(iulog,*) '   mpas_scalar_vadv_order = ', mpas_scalar_vadv_order
      write(iulog,*) '   mpas_scalar_advection = ', mpas_scalar_advection
      write(iulog,*) '   mpas_positive_definite = ', mpas_positive_definite
      write(iulog,*) '   mpas_monotonic = ', mpas_monotonic
      write(iulog,*) '   mpas_coef_3rd_order = ', mpas_coef_3rd_order
      write(iulog,*) '   mpas_smagorinsky_coef = ', mpas_smagorinsky_coef
      write(iulog,*) '   mpas_mix_full = ', mpas_mix_full
      write(iulog,*) '   mpas_epssm = ', mpas_epssm
      write(iulog,*) '   mpas_smdiv = ', mpas_smdiv
      write(iulog,*) '   mpas_apvm_upwinding = ', mpas_apvm_upwinding
      write(iulog,*) '   mpas_h_ScaleWithMesh = ', mpas_h_ScaleWithMesh
      write(iulog,*) '   mpas_zd = ', mpas_zd
      write(iulog,*) '   mpas_xnutr = ', mpas_xnutr
      write(iulog,*) '   mpas_cam_coef = ', mpas_cam_coef
      write(iulog,*) '   mpas_cam_damping_levels = ', mpas_cam_damping_levels
      write(iulog,*) '   mpas_rayleigh_damp_u = ', mpas_rayleigh_damp_u
      write(iulog,*) '   mpas_rayleigh_damp_u_timescale_days = ', mpas_rayleigh_damp_u_timescale_days
      write(iulog,*) '   mpas_number_rayleigh_damp_u_levels = ', mpas_number_rayleigh_damp_u_levels
      write(iulog,*) '   mpas_apply_lbcs = ', mpas_apply_lbcs
      write(iulog,*) '   mpas_jedi_da = ', mpas_jedi_da
      write(iulog,*) '   mpas_block_decomp_file_prefix = ', trim(mpas_block_decomp_file_prefix)
      write(iulog,*) '   mpas_do_restart = ', mpas_do_restart
      write(iulog,*) '   mpas_print_global_minmax_vel = ', mpas_print_global_minmax_vel
      write(iulog,*) '   mpas_print_detailed_minmax_vel = ', mpas_print_detailed_minmax_vel
      write(iulog,*) '   mpas_print_global_minmax_sca = ', mpas_print_global_minmax_sca
   end if

end subroutine cam_mpas_namelist_read

!-----------------------------------------------------------------------
!  routine set_dry_mass
!
!> \brief Scale dry air mass
!> \author Bill Skamarock, Miles Curry
!> \date   25 April 2021
!> \details Given a target dry air mass surface pressure,
!> target_avg_dry_surface_pressure, scale the current dry air mass so
!> that the average dry surface pressure equals
!> target_avg_dry_surface_pressure. Water vapor is scaled for mass-
!> conservation; all other tracer mixing ratios are unaltered
!> (i.e. tracer mass is not conserved but gradients are during the
!> dry mass scaling process)
!
!-----------------------------------------------------------------------
subroutine set_dry_mass(dyn_in, target_avg_dry_surface_pressure)

   use mpas_constants, only : rgas, gravity, p0, Rv_over_Rd => rvord

   type(dyn_import_t), intent(in) :: dyn_in
   real(r8), intent(in) :: target_avg_dry_surface_pressure

   integer :: i, k
   integer :: nCellsSolve

   real(r8), pointer :: theta_m(:,:) ! Moist potential temperature [K]  (nver,ncol)
   real(r8), pointer :: zint(:,:)    ! Geometric height [m]
   real(r8), pointer :: areaCell(:)  ! cell area (m^2)
   real(r8), pointer :: theta(:,:)   ! Potential temperature [K]        (nver,ncol)
   real(r8), pointer :: rho(:,:)     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), pointer :: rho_zz(:,:)  ! Dry density [kg/m^3]
                                     ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), pointer :: tracers(:,:,:) ! Tracers [kg/kg dry air]       (nq,nver,ncol)
   real(r8), pointer :: zz(:,:)      ! Vertical coordinate metric [dimensionless]
                                     ! at layer midpoints               (nver,ncol)

   real(r8), allocatable :: preliminary_dry_surface_pressure(:), p_top(:), pm(:)
   real(r8) :: preliminary_avg_dry_surface_pressure, scaled_avg_dry_surface_pressure
   real(r8) :: scaling_ratio
   real(r8) :: sphere_surface_area

   integer :: ixqv,ierr

   character(len=*), parameter :: subname = 'dyn_comp:set_dry_mass'

   nCellsSolve = dyn_in % nCellsSolve
   ixqv        = dyn_in % index_qv
   theta_m    => dyn_in % theta_m
   theta      => dyn_in % theta
   zint       => dyn_in % zint
   areaCell   => dyn_in % areaCell
   rho        => dyn_in % rho
   rho_zz     => dyn_in % rho_zz
   zz         => dyn_in % zz
   tracers    => dyn_in % tracers

   allocate( p_top(nCellsSolve), preliminary_dry_surface_pressure(nCellsSolve), pm(plev), stat=ierr)
   if( ierr /= 0 ) call endrun(subname//': failed to allocate  arrays preliminary_dry_surface_pressure and pm')
   ! (1) calculate pressure at the lid
   do i=1, nCellsSolve
      p_top(i) = p0*(rgas*rho(plev,i)*theta_m(plev,i)/p0)**(cpair/(cpair-rgas))
      p_top(i) = p_top(i) - gravity*0.5_r8*(zint(plev+1,i)-zint(plev,i))*rho(plev,i)*(1.0_r8+tracers(ixqv,plev,i))
   end do

   ! (2) integrate dry mass in column
   do i=1, nCellsSolve
      preliminary_dry_surface_pressure(i) = 0.0_r8
      do k=1, plev
         preliminary_dry_surface_pressure(i) = preliminary_dry_surface_pressure(i) + gravity*(zint(k+1,i)-zint(k,i))*rho(k,i)
      end do
   end do

   ! (3) compute average global dry surface pressure
   preliminary_dry_surface_pressure(1:nCellsSolve) =  preliminary_dry_surface_pressure(1:nCellsSolve)*areaCell(1:nCellsSolve)
   sphere_surface_area = cam_mpas_global_sum_real(areaCell(1:nCellsSolve))
   preliminary_avg_dry_surface_pressure = cam_mpas_global_sum_real(preliminary_dry_surface_pressure(1:nCellsSolve)) &
                                                                   /sphere_surface_area

   if (masterproc) then
       write(iulog,*) '---------------------------- set_dry_mass ----------------------------'
       write(iulog,*) 'Initial dry globally average surface pressure = ', preliminary_avg_dry_surface_pressure/100._r8, 'hPa'
       write(iulog,*) 'target dry globally avg surface pressure = ', target_avg_dry_surface_pressure/100._r8, 'hPa'
   end if

   ! (4) scale dry air density
   scaling_ratio = target_avg_dry_surface_pressure / preliminary_avg_dry_surface_pressure
   rho(:,:) = rho(:,:)*scaling_ratio

   ! (4a) recompute dry mass after scaling
   do i = 1, nCellsSolve
       preliminary_dry_surface_pressure(i) = 0.0_r8
       do k = 1, plev
           preliminary_dry_surface_pressure(i) = preliminary_dry_surface_pressure(i) + gravity*(zint(k+1,i)-zint(k,i))*rho(k,i)
       end do
   end do
   preliminary_dry_surface_pressure(1:nCellsSolve) = preliminary_dry_surface_pressure(1:nCellsSolve)*areaCell(1:nCellsSolve)
   scaled_avg_dry_surface_pressure = cam_mpas_global_sum_real(preliminary_dry_surface_pressure(1:nCellsSolve)) &
                                                              / sphere_surface_area

   if (masterproc) then
      write(iulog,*) 'Average dry global surface pressure after scaling = ', scaled_avg_dry_surface_pressure/100._r8, 'hPa'
      write(iulog,*) 'Change in dry surface pressure = ', scaled_avg_dry_surface_pressure-preliminary_avg_dry_surface_pressure,'Pa'
   end if

   ! (5) reset qv to conserve mass
   tracers(ixqv,:,1:nCellsSolve) = tracers(ixqv,:,1:nCellsSolve)/scaling_ratio

   ! (6) integrate down the column to compute full pressure given the density and qv
   do i=1,nCellsSolve
      pm(plev) = p_top(i) + 0.5_r8*(zint(plev+1,i)-zint(plev,i))*gravity*rho(plev,i)*(1.0_r8+tracers(ixqv,plev,i))
      do k=plev-1,1,-1
         pm(k) = pm(k+1) + 0.5_r8*(zint(k+2,i)-zint(k+1,i))*gravity*rho(k+1,i)*(1.0_r8+tracers(ixqv,k+1,i)) &
                         + 0.5_r8*(zint(k+1,i)-zint(k  ,i))*gravity*rho(k  ,i)*(1.0_r8+tracers(ixqv,k  ,i))
      end do

   ! (7) compute theta_m from the state equation, compute rho_zz and theta while we are here

      do k=1,plev
         theta_m(k,i) = (pm(k)/p0)**((cpair-rgas)/cpair)*p0/rgas/rho(k,i)
         theta(k,i) = theta_m(k,i)/(1.0_r8 + Rv_over_Rd * tracers(ixqv,k,i))
         rho_zz(k,i) = rho(k,i)/zz(k,i)
      end do
   end do

  deallocate( p_top, preliminary_dry_surface_pressure, pm )

end subroutine set_dry_mass

end module dyn_comp
