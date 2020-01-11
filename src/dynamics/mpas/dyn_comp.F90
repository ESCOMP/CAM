#define MPAS_DEBUG_WRITE(print_task, x) if (iam == (print_task)) write(iulog,*) 'MPAS_DEBUG '//subname//' ', (x)

module dyn_comp

! CAM component interfaces to the MPAS Dynamical Core

use shr_kind_mod,       only: r8=>shr_kind_r8
use spmd_utils,         only: iam, masterproc, mpicom, npes
use physconst,          only: pi

use pmgrid,             only: plev, plevp
use constituents,       only: pcnst, cnst_name, cnst_read_iv

use cam_control_mod,    only: initial_run
use cam_initfiles,      only: initial_file_get_id, topo_file_get_id

use dyn_grid,           only: nCellsSolve, nVertLevelsSolve

use cam_grid_support,   only: cam_grid_id, cam_grid_get_gcid, &
                              cam_grid_dimensions, cam_grid_get_dim_names, &
                              cam_grid_get_latvals, cam_grid_get_lonvals,  &
                              max_hcoordname_len
use cam_map_utils,      only: iMap

use inic_analytic,      only: analytic_ic_active, analytic_ic_set_ic
use dyn_tests_utils,    only: vcoord=>vc_height

use cam_history,        only: addfld, add_default, horiz_only, register_vector_field, &
                              outfld, hist_fld_active
use cam_history_support, only: date2yyyymmdd, sec2hms, nday2str, &
                               max_fieldname_len
use cam_pio_utils,      only: clean_iodesc_list

use ncdio_atm,          only: infld
use pio,                only: file_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR, &
                              pio_inq_dimid, pio_inq_dimlen, PIO_NOERR

use time_manager,       only: get_start_date, get_stop_date, get_run_duration, &
                              timemgr_get_calendar_cf

use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun


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

type dyn_import_t
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
   ! Invariant -- the vertical coordinate in MPAS-A is a height coordinate
   !
   real(r8), dimension(:,:),   pointer :: zint    ! Geometric height [m]
                                                  ! at layer interfaces            (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: zz      ! Vertical coordinate metric [dimensionless]
                                                  ! at layer midpoints               (nver,ncol)

   !
   ! State that may be directly derived from dycore prognostic state
   !
   real(r8), dimension(:,:),   pointer :: theta   ! Potential temperature [K]        (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ux      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), dimension(:,:),   pointer :: uy      ! Meridional veloc at center [m/s] (nver,ncol)
   real(r8), dimension(:,:),   pointer :: pmid    ! Full non-hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nver,ncol)

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
   ! State that is directly prognosed by the dycore
   !
   real(r8), dimension(:,:),   pointer :: uperp   ! Normal velocity at edges [m/s]  (nver,nedge)
   real(r8), dimension(:,:),   pointer :: w       ! Vertical velocity [m/s]        (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: theta_m ! Moist potential temperature [K]  (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_zz  ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), dimension(:,:,:), pointer :: tracers ! Tracers [kg/kg dry air]       (nq,nver,ncol)

   !
   ! Invariant -- the vertical coordinate in MPAS-A is a height coordinate
   !
   real(r8), dimension(:,:),   pointer :: zint    ! Geometric height [m]
                                                  ! at layer interfaces            (nver+1,ncol)

   !
   ! State that may be directly derived from dycore prognostic state
   !
   real(r8), dimension(:,:),   pointer :: theta   ! Potential temperature [K]        (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ux      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), dimension(:,:),   pointer :: uy      ! Meridional veloc at center [m/s] (nver,ncol)
   real(r8), dimension(:,:),   pointer :: pmid    ! Full non-hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nver,ncol)
end type dyn_export_t

real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8


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
   integer :: ierr
   integer, dimension(2) :: logUnits   ! stdout and stderr for MPAS logging
   integer :: yr, mon, day, tod, ndate, nday, nsec
   character(len=10) :: date_str
   character(len=8)  :: tod_str

   character(len=*), parameter :: subname = 'dyn_comp::dyn_readnl'
   !----------------------------------------------------------------------------


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   logUnits(1) = iulog
   logUnits(2) = getunit()

   call cam_mpas_init_phase1(mpicom, endrun, logUnits)

   ! read namelist
   ierr = cam_mpas_namelist_read(NLFileName, domain_ptr % configs)
   if ( ierr /= 0 ) then
      call endrun(subname//': FATAL: Namelist setup failed for MPAS-A dycore')
   end if

   ! Set config_start_date, etc. (these will not appear in the dycore namelist)
   call get_start_date(yr, mon, day, tod)
   ndate = yr*10000 + mon*100 + day
   call mpas_pool_add_config(domain_ptr % configs, 'config_start_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

   call get_stop_date(yr, mon, day, tod)
   ndate = yr*10000 + mon*100 + day
   call mpas_pool_add_config(domain_ptr % configs, 'config_stop_time', date2yyyymmdd(ndate)//'_'//sec2hms(tod))

   call get_run_duration(nday, nsec)
   call mpas_pool_add_config(domain_ptr % configs, 'config_run_duration', trim(nday2str(nday))//'_'//sec2hms(nsec))

   call mpas_pool_add_config(domain_ptr % configs, 'config_restart_timestamp_name', 'restart_timestamp')
   call mpas_pool_add_config(domain_ptr % configs, 'config_IAU_option', 'off')

   call cam_mpas_init_phase2(pio_subsystem, endrun, timemgr_get_calendar_cf())

end subroutine dyn_readnl

!=========================================================================================

subroutine dyn_register()

   ! Register fields that are computed by the dycore and passed to the physics via the
   ! physics buffer.

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver

   character(len=*), parameter :: subname = 'dyn_comp::dyn_register'
   !----------------------------------------------------------------------------


end subroutine dyn_register

!=========================================================================================

subroutine dyn_init(dyn_in, dyn_out)

   use cam_mpas_subdriver, only : domain_ptr
   use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
   use mpas_derived_types, only : mpas_pool_type

   ! arguments:
   type(dyn_import_t), intent(out)  :: dyn_in
   type(dyn_export_t), intent(out)  :: dyn_out

   character(len=*), parameter :: subname = 'dyn_comp::dyn_init'

   type(mpas_pool_type), pointer :: mesh_pool
   type(mpas_pool_type), pointer :: state_pool
   type(mpas_pool_type), pointer :: diag_pool
   type(mpas_pool_type), pointer :: tend_physics_pool

   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_physics_pool)


   !
   ! Let dynamics import state point to memory managed by MPAS-Atmosphere
   !
   call mpas_pool_get_array(state_pool, 'u',                      dyn_in % uperp,   timeLevel=2)
   call mpas_pool_get_array(state_pool, 'w',                      dyn_in % w,       timeLevel=2)
   call mpas_pool_get_array(state_pool, 'theta_m',                dyn_in % theta_m, timeLevel=2)
   call mpas_pool_get_array(state_pool, 'rho_zz',                 dyn_in % rho_zz,  timeLevel=2)
   call mpas_pool_get_array(state_pool, 'scalars',                dyn_in % tracers, timeLevel=2)

   call mpas_pool_get_array(mesh_pool,  'zgrid',                  dyn_in % zint)
   call mpas_pool_get_array(mesh_pool,  'zz',                     dyn_in % zz)

   call mpas_pool_get_array(diag_pool,  'theta',                  dyn_in % theta)
   call mpas_pool_get_array(diag_pool,  'rho',                    dyn_in % rho)
   call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      dyn_in % ux)
   call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', dyn_in % uy)
   call mpas_pool_get_array(diag_pool,  'pressure',               dyn_in % pmid)

   call mpas_pool_get_array(tend_physics_pool, 'tend_ru_physics',     dyn_in % ru_tend)
   call mpas_pool_get_array(tend_physics_pool, 'tend_rtheta_physics', dyn_in % rtheta_tend)
   call mpas_pool_get_array(tend_physics_pool, 'tend_rho_physics',    dyn_in % rho_tend)


   !
   ! Let dynamics export state point to memory managed by MPAS-Atmosphere
   !
   call mpas_pool_get_array(state_pool, 'u',                      dyn_out % uperp,   timeLevel=1)
   call mpas_pool_get_array(state_pool, 'w',                      dyn_out % w,       timeLevel=1)
   call mpas_pool_get_array(state_pool, 'theta_m',                dyn_out % theta_m, timeLevel=1)
   call mpas_pool_get_array(state_pool, 'rho_zz',                 dyn_out % rho_zz,  timeLevel=1)
   call mpas_pool_get_array(state_pool, 'scalars',                dyn_out % tracers, timeLevel=1)

   call mpas_pool_get_array(mesh_pool,  'zgrid',                  dyn_out % zint)

   call mpas_pool_get_array(diag_pool,  'theta',                  dyn_out % theta)
   call mpas_pool_get_array(diag_pool,  'rho',                    dyn_out % rho)
   call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      dyn_out % ux)
   call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', dyn_out % uy)
   call mpas_pool_get_array(diag_pool,  'pressure',               dyn_out % pmid)

   if (initial_run) then
      call read_inidat(dyn_in)
      call clean_iodesc_list()
   end if

end subroutine dyn_init

!=========================================================================================

subroutine dyn_run(dyn_in, dyn_out)

   ! Advances the dynamics state provided in dyn_in by one physics
   ! timestep to produce dynamics state held in dyn_out.

   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   character(len=*), parameter :: subname = 'dyn_comp::dyn_run'
   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(0, 'begin '//subname)


end subroutine dyn_run

!=========================================================================================

subroutine dyn_final(dyn_in, dyn_out)

   ! Deallocates the dynamics import and export states, and finalizes
   ! the MPAS dycore.

   type (dyn_import_t), intent(inout)  :: dyn_in
   type (dyn_export_t), intent(inout)  :: dyn_out

   character(len=*), parameter :: subname = 'dyn_comp::dyn_final'
   !----------------------------------------------------------------------------

   !
   ! Prevent any further access to MPAS-Atmosphere memory
   !
   nullify(dyn_in % uperp)
   nullify(dyn_in % w)
   nullify(dyn_in % theta_m)
   nullify(dyn_in % rho_zz)
   nullify(dyn_in % tracers)
   nullify(dyn_in % zint)
   nullify(dyn_in % zz)
   nullify(dyn_in % theta)
   nullify(dyn_in % rho)
   nullify(dyn_in % ux)
   nullify(dyn_in % uy)
   nullify(dyn_in % pmid)
   nullify(dyn_in % ru_tend)
   nullify(dyn_in % rtheta_tend)
   nullify(dyn_in % rho_tend)

   !
   ! Prevent any further access to MPAS-Atmosphere memory
   !
   nullify(dyn_out % uperp)
   nullify(dyn_out % w)
   nullify(dyn_out % theta_m)
   nullify(dyn_out % rho_zz)
   nullify(dyn_out % tracers)
   nullify(dyn_out % zint)
   nullify(dyn_out % theta)
   nullify(dyn_out % rho)
   nullify(dyn_out % ux)
   nullify(dyn_out % uy)
   nullify(dyn_out % pmid)

end subroutine dyn_final

!=========================================================================================
! Private routines.
!=========================================================================================

subroutine read_inidat(dyn_in)

   use dyn_grid, only : nCellsSolve, nEdgesSolve, nVertLevelsSolve
   use cam_mpas_subdriver, only : domain_ptr

   ! Set initial conditions.  Either from analytic expressions or read from file.

   ! arguments
   type(dyn_import_t), target, intent(inout) :: dyn_in

   ! Local variables
   integer :: klev, icol

   type(file_desc_t), pointer :: fh_ini
   type(file_desc_t), pointer :: fh_topo

   real(r8), dimension(:,:),   pointer :: uperp   ! Normal velocity at edges [m/s]  (nver,nedge)
   real(r8), dimension(:,:),   pointer :: w       ! Vertical velocity [m/s]        (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: theta_m ! Moist potential temperature [K]  (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_zz  ! Dry density [kg/m^3]
                                                  ! divided by d(zeta)/dz            (nver,ncol)
   real(r8), dimension(:,:,:), pointer :: tracers ! Tracers [kg/kg dry air]       (nq,nver,ncol)
   real(r8), dimension(:,:),   pointer :: zint    ! Geometric height [m]
                                                  ! at layer interfaces            (nver+1,ncol)
   real(r8), dimension(:,:),   pointer :: zz      ! Vertical coordinate metric [dimensionless]
                                                  ! at layer midpoints               (nver,ncol)
   real(r8), dimension(:,:),   pointer :: theta   ! Potential temperature [K]        (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho     ! Dry density [kg/m^3]             (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ux      ! Zonal veloc at center [m/s]      (nver,ncol)
   real(r8), dimension(:,:),   pointer :: uy      ! Meridional veloc at center [m/s] (nver,ncol)
   real(r8), dimension(:,:),   pointer :: pmid    ! Full non-hydrostatic pressure [Pa]
                                                  ! at layer midpoints               (nver,ncol)
   real(r8), dimension(:,:),   pointer :: ru_tend ! Normal horizontal momentum tendency
                                                  ! from physics [kg/m^2/s]         (nver,nedge)
   real(r8), dimension(:,:),   pointer :: rtheta_tend ! Tendency of rho*theta/zz
                                                      ! from physics [kg K/m^3/s]    (nver,ncol)
   real(r8), dimension(:,:),   pointer :: rho_tend ! Dry air density tendency
                                                   ! from physics [kg/m^3/s]         (nver,ncol)

   real(r8) :: t

   real(r8), parameter, dimension(32) :: pmid_ref = [    370.0507, &
                772.0267, &
               1460.0401, &
               2456.1190, &
               3510.4428, &
               4394.9775, &
               5258.6650, &
               6260.5462, &
               7505.1291, &
               8937.0345, &
              10513.9623, &
              12368.9543, &
              14551.4143, &
              17119.0922, &
              20139.5649, &
              23688.5536, &
              27855.9967, &
              32754.2135, &
              38514.6948, &
              45288.6443, &
              53254.6026, &
              61865.3664, &
              70121.4348, &
              77405.1773, &
              83215.3208, &
              87125.5307, &
              89903.9705, &
              92494.3081, &
              94874.9828, &
              97027.0667, &
              98931.3386, &
             100572.0510 ]

   real(r8), parameter, dimension(33) :: pint_ref = [  228.8976, &
              511.2039, &
             1032.8496, &
             1887.2306, &
             3025.0074, &
             3995.8781, &
             4794.0770, &
             5723.2531, &
             6797.8393, &
             8212.4188, &
             9661.6502, &
            11366.2745, &
            13371.6342, &
            15731.1944, &
            18506.9901, &
            21772.1397, &
            25604.9675, &
            30107.0258, &
            35401.4011, &
            41627.9885, &
            48949.3000, &
            57559.9052, &
            66170.8275, &
            74072.0421, &
            80738.3125, &
            85692.3292, &
            88558.7323, &
            91249.2088, &
            93739.4075, &
            96010.5582, &
            98043.5752, &
            99819.1020, &
            101325.0000 ]

   real(r8), parameter, dimension(32) :: zmid_ref = [ 38524.2527, &
             33230.5418, &
             28850.1087, &
             25327.9309, &
             22900.7199, &
             21425.0307, &
             20281.2855, &
             19173.2056, &
             18028.1748, &
             16913.4272, &
             15882.8868, &
             14852.4477, &
             13821.9328, &
             12791.3613, &
             11760.8767, &
             10726.9401, &
              9674.2358, &
              8589.4962, &
              7470.6318, &
              6316.6734, &
              5126.5126, &
              3987.1676, &
              3008.7228, &
              2220.0977, &
              1632.9889, &
              1256.2718, &
               998.0675, &
               763.0978, &
               551.8000, &
               364.4534, &
               201.4641, &
                63.0587 ]

   real(r8), parameter, dimension(33) :: zint_ref = [ 41426.9190, &
             35621.5865, &
             30839.4971, &
             26860.7204, &
             23795.1414, &
             22006.2985, &
             20843.7630, &
             19718.8080, &
             18627.6033, &
             17428.7463, &
             16398.1081, &
             15367.6655, &
             14337.2299, &
             13306.6357, &
             12276.0870, &
             11245.6664, &
             10208.2138, &
              9140.2578, &
              8038.7346, &
              6902.5289, &
              5730.8179, &
              4522.2073, &
              3452.1278, &
              2565.3178, &
              1874.8775, &
              1391.1003, &
              1121.4432, &
               874.6918, &
               651.5038, &
               452.0961, &
               276.8108, &
               126.1175, &
                 0.0000 ]

   character(len=*), parameter :: subname = 'dyn_comp::read_inidat'


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   
   uperp => dyn_in % uperp
   w => dyn_in % w
   theta_m => dyn_in % theta_m
   rho_zz => dyn_in % rho_zz
   tracers => dyn_in % tracers
   zint => dyn_in % zint
   zz => dyn_in % zz
   theta => dyn_in % theta
   rho => dyn_in % rho
   ux => dyn_in % ux
   uy => dyn_in % uy
   pmid => dyn_in % pmid
   ru_tend => dyn_in % ru_tend
   rtheta_tend => dyn_in % rtheta_tend
   rho_tend => dyn_in % rho_tend

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()


   !
   ! Set placeholder state
   !

   ! No winds
   uperp(:,1:nEdgesSolve) = 0.0_r8
   w(:,1:nCellsSolve) = 0.0_r8
   ux(:,1:nCellsSolve) = 0.0_r8
   uy(:,1:nCellsSolve) = 0.0_r8

   ! No tracers
   tracers(:,:,1:nCellsSolve) = 0.0_r8

   ! Standard atmosphere (for MPAS, k=1 is the layer closest to the surface)
   do klev=1,nVertLevelsSolve+1
      zint(klev,1:nCellsSolve) = zint_ref(nVertLevelsSolve+1-klev+1)
   end do

   do klev=1,nVertLevelsSolve
      pmid(klev,1:nCellsSolve) = pmid_ref(nVertLevelsSolve-klev+1)
   end do

   do klev=1,nVertLevelsSolve
      t = (zint_ref(nVertLevelsSolve-klev+1) - zint_ref(nVertLevelsSolve-klev+1+1)) * 9.806/287.0 &
          / log(pint_ref(nVertLevelsSolve-klev+1+1) / pint_ref(nVertLevelsSolve-klev+1))
      theta(klev,1:nCellsSolve) = t * (1.0e5 / pmid_ref(nVertLevelsSolve-klev+1))**(287.0/1003.0)
      rho(klev,1:nCellsSolve) = pmid_ref(nVertLevelsSolve-klev+1) / (287.0 * t)
   end do

   theta_m(:,1:nCellsSolve) = theta(:,1:nCellsSolve)    ! With no moisture, theta_m := theta
   rho_zz(:,1:nCellsSolve) = rho(:,1:nCellsSolve) / zz(:,1:nCellsSolve)

   ! Initially, no tendencies from physics
   ru_tend(:,1:nEdgesSolve) = 0.0_r8
   rtheta_tend(:,1:nCellsSolve) = 0.0_r8
   rho_tend(:,1:nCellsSolve) = 0.0_r8

end subroutine read_inidat

!========================================================================================

!-----------------------------------------------------------------------
!  routine cam_mpas_namelist_read
!
!> \brief Reads MPAS-A dycore namelists and adds the namelists to the MPAS configPool
!> \details
!>  Given the name of a file containing namelists and an MPAS pool, reads the dycore
!>  namelists from that file and adds the namelist options to the pool.
!
!>  Only the CAM masterproc actually opens and reads from the specified file. Upon return,
!>  if no errors were encountered, all MPI ranks have valid namelists in their configPool.
!
!>  A value of zero is returned if no errors were encountered, and a non-zero value is returned
!>  if any errors were encountered in reading the namelist file.
!
!   WARNING: This routine was auto-generated based on the MPAS-Atmosphere Registry.xml file
!
!            Rather than editing this function directly, edit core_atmosphere/Registry.xml
!            and re-run the build_cam_namelist.py script to regenerate this function.
!
!-----------------------------------------------------------------------
function cam_mpas_namelist_read(namelistFilename, configPool) result(ierr)

   use units, only : getunit, freeunit
   use spmd_utils, only : mpicom, masterproc, masterprocid, &
                          mpi_integer, mpi_real8,  mpi_logical, mpi_character, mpi_success
   use shr_kind_mod, only : shr_kind_r8
   use cam_logfile, only : iulog
   use namelist_utils, only : find_group_name

   use mpas_derived_types, only : mpas_pool_type
   use mpas_kind_types, only : StrKIND
   use mpas_pool_routines, only : mpas_pool_add_config

   implicit none

   character(len=*), intent(in) :: namelistFilename
   type (mpas_pool_type), intent(inout) :: configPool

   integer :: ierr   ! Return value

   integer :: unitNumber

   integer :: mpi_ierr

   character (len=StrKIND) :: mpas_time_integration = 'SRK3'
   integer                 :: mpas_time_integration_order = 2
   real (kind=shr_kind_r8) :: mpas_dt = 720.0
   logical                 :: mpas_split_dynamics_transport = .true.
   integer                 :: mpas_number_of_sub_steps = 2
   integer                 :: mpas_dynamics_split_steps = 3
   real (kind=shr_kind_r8) :: mpas_h_mom_eddy_visc2 = 0.0
   real (kind=shr_kind_r8) :: mpas_h_mom_eddy_visc4 = 0.0
   real (kind=shr_kind_r8) :: mpas_v_mom_eddy_visc2 = 0.0
   real (kind=shr_kind_r8) :: mpas_h_theta_eddy_visc2 = 0.0
   real (kind=shr_kind_r8) :: mpas_h_theta_eddy_visc4 = 0.0
   real (kind=shr_kind_r8) :: mpas_v_theta_eddy_visc2 = 0.0
   character (len=StrKIND) :: mpas_horiz_mixing = '2d_smagorinsky'
   real (kind=shr_kind_r8) :: mpas_len_disp = 120000.0
   real (kind=shr_kind_r8) :: mpas_visc4_2dsmag = 0.05
   real (kind=shr_kind_r8) :: mpas_del4u_div_factor = 10.0
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
   real (kind=shr_kind_r8) :: mpas_coef_3rd_order = 0.25
   real (kind=shr_kind_r8) :: mpas_smagorinsky_coef = 0.125
   logical                 :: mpas_mix_full = .true.
   real (kind=shr_kind_r8) :: mpas_epssm = 0.1
   real (kind=shr_kind_r8) :: mpas_smdiv = 0.1
   real (kind=shr_kind_r8) :: mpas_apvm_upwinding = 0.5
   logical                 :: mpas_h_ScaleWithMesh = .true.
   integer                 :: mpas_num_halos = 2
   real (kind=shr_kind_r8) :: mpas_zd = 22000.0
   real (kind=shr_kind_r8) :: mpas_xnutr = 0.2
   character (len=StrKIND) :: mpas_block_decomp_file_prefix = 'x1.40962.graph.info.part.'
   integer                 :: mpas_number_of_blocks = 0
   logical                 :: mpas_explicit_proc_decomp = .false.
   character (len=StrKIND) :: mpas_proc_decomp_file_prefix = 'graph.info.part.'
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
           mpas_h_ScaleWithMesh, &
           mpas_num_halos

   namelist /damping/ &
           mpas_zd, &
           mpas_xnutr

   namelist /decomposition/ &
           mpas_block_decomp_file_prefix, &
           mpas_number_of_blocks, &
           mpas_explicit_proc_decomp, &
           mpas_proc_decomp_file_prefix

   namelist /restart/ &
           mpas_do_restart

   namelist /printout/ &
           mpas_print_global_minmax_vel, &
           mpas_print_detailed_minmax_vel, &
           mpas_print_global_minmax_sca

   if (masterproc) then
      write(iulog,*) 'Reading MPAS-A dycore namelist from ', trim(namelistFilename)
      unitNumber = getunit()
      open(unit=unitNumber, file=trim(namelistFilename), status='old', form='formatted')
   end if

   !
   ! Read namelist group &nhyd_model
   !
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'nhyd_model', status=ierr)
      if (ierr == 0) then
         read(unitNumber, nhyd_model, iostat=ierr)
      else
         close(unit=unitNumber)
         call freeunit(unitNumber)
      end if
   end if
   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (ierr /= 0) then
      if (masterproc) then
         write(iulog,*) 'Failed to read namelist group &nhyd_model'
      end if
      return
   end if
   call mpi_bcast(mpas_time_integration, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_time_integration'
      end if
      return
   end if
   call mpi_bcast(mpas_time_integration_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_time_integration_order'
      end if
      return
   end if
   call mpi_bcast(mpas_dt, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_dt'
      end if
      return
   end if
   call mpi_bcast(mpas_split_dynamics_transport, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_split_dynamics_transport'
      end if
      return
   end if
   call mpi_bcast(mpas_number_of_sub_steps, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_number_of_sub_steps'
      end if
      return
   end if
   call mpi_bcast(mpas_dynamics_split_steps, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_dynamics_split_steps'
      end if
      return
   end if
   call mpi_bcast(mpas_h_mom_eddy_visc2, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_h_mom_eddy_visc2'
      end if
      return
   end if
   call mpi_bcast(mpas_h_mom_eddy_visc4, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_h_mom_eddy_visc4'
      end if
      return
   end if
   call mpi_bcast(mpas_v_mom_eddy_visc2, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_v_mom_eddy_visc2'
      end if
      return
   end if
   call mpi_bcast(mpas_h_theta_eddy_visc2, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_h_theta_eddy_visc2'
      end if
      return
   end if
   call mpi_bcast(mpas_h_theta_eddy_visc4, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_h_theta_eddy_visc4'
      end if
      return
   end if
   call mpi_bcast(mpas_v_theta_eddy_visc2, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_v_theta_eddy_visc2'
      end if
      return
   end if
   call mpi_bcast(mpas_horiz_mixing, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_horiz_mixing'
      end if
      return
   end if
   call mpi_bcast(mpas_len_disp, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_len_disp'
      end if
      return
   end if
   call mpi_bcast(mpas_visc4_2dsmag, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_visc4_2dsmag'
      end if
      return
   end if
   call mpi_bcast(mpas_del4u_div_factor, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_del4u_div_factor'
      end if
      return
   end if
   call mpi_bcast(mpas_w_adv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_w_adv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_theta_adv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_theta_adv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_scalar_adv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_scalar_adv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_u_vadv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_u_vadv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_w_vadv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_w_vadv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_theta_vadv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_theta_vadv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_scalar_vadv_order, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_scalar_vadv_order'
      end if
      return
   end if
   call mpi_bcast(mpas_scalar_advection, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_scalar_advection'
      end if
      return
   end if
   call mpi_bcast(mpas_positive_definite, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_positive_definite'
      end if
      return
   end if
   call mpi_bcast(mpas_monotonic, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_monotonic'
      end if
      return
   end if
   call mpi_bcast(mpas_coef_3rd_order, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_coef_3rd_order'
      end if
      return
   end if
   call mpi_bcast(mpas_smagorinsky_coef, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_smagorinsky_coef'
      end if
      return
   end if
   call mpi_bcast(mpas_mix_full, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_mix_full'
      end if
      return
   end if
   call mpi_bcast(mpas_epssm, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_epssm'
      end if
      return
   end if
   call mpi_bcast(mpas_smdiv, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_smdiv'
      end if
      return
   end if
   call mpi_bcast(mpas_apvm_upwinding, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_apvm_upwinding'
      end if
      return
   end if
   call mpi_bcast(mpas_h_ScaleWithMesh, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_h_ScaleWithMesh'
      end if
      return
   end if
   call mpi_bcast(mpas_num_halos, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_num_halos'
      end if
      return
   end if

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
   call mpas_pool_add_config(configPool, 'config_num_halos', mpas_num_halos)

   !
   ! Read namelist group &damping
   !
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'damping', status=ierr)
      if (ierr == 0) then
         read(unitNumber, damping, iostat=ierr)
      else
         close(unit=unitNumber)
         call freeunit(unitNumber)
      end if
   end if
   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (ierr /= 0) then
      if (masterproc) then
         write(iulog,*) 'Failed to read namelist group &damping'
      end if
      return
   end if
   call mpi_bcast(mpas_zd, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_zd'
      end if
      return
   end if
   call mpi_bcast(mpas_xnutr, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_xnutr'
      end if
      return
   end if

   call mpas_pool_add_config(configPool, 'config_zd', mpas_zd)
   call mpas_pool_add_config(configPool, 'config_xnutr', mpas_xnutr)

   !
   ! Read namelist group &decomposition
   !
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'decomposition', status=ierr)
      if (ierr == 0) then
         read(unitNumber, decomposition, iostat=ierr)
      else
         close(unit=unitNumber)
         call freeunit(unitNumber)
      end if
   end if
   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (ierr /= 0) then
      if (masterproc) then
         write(iulog,*) 'Failed to read namelist group &decomposition'
      end if
      return
   end if
   call mpi_bcast(mpas_block_decomp_file_prefix, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_block_decomp_file_prefix'
      end if
      return
   end if
   call mpi_bcast(mpas_number_of_blocks, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_number_of_blocks'
      end if
      return
   end if
   call mpi_bcast(mpas_explicit_proc_decomp, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_explicit_proc_decomp'
      end if
      return
   end if
   call mpi_bcast(mpas_proc_decomp_file_prefix, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_proc_decomp_file_prefix'
      end if
      return
   end if

   call mpas_pool_add_config(configPool, 'config_block_decomp_file_prefix', mpas_block_decomp_file_prefix)
   call mpas_pool_add_config(configPool, 'config_number_of_blocks', mpas_number_of_blocks)
   call mpas_pool_add_config(configPool, 'config_explicit_proc_decomp', mpas_explicit_proc_decomp)
   call mpas_pool_add_config(configPool, 'config_proc_decomp_file_prefix', mpas_proc_decomp_file_prefix)

   !
   ! Read namelist group &restart
   !
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'restart', status=ierr)
      if (ierr == 0) then
         read(unitNumber, restart, iostat=ierr)
      else
         close(unit=unitNumber)
         call freeunit(unitNumber)
      end if
   end if
   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (ierr /= 0) then
      if (masterproc) then
         write(iulog,*) 'Failed to read namelist group &restart'
      end if
      return
   end if
   call mpi_bcast(mpas_do_restart, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_do_restart'
      end if
      return
   end if

   call mpas_pool_add_config(configPool, 'config_do_restart', mpas_do_restart)

   !
   ! Read namelist group &printout
   !
   if (masterproc) then
      rewind(unitNumber)
      call find_group_name(unitNumber, 'printout', status=ierr)
      if (ierr == 0) then
         read(unitNumber, printout, iostat=ierr)
      else
         close(unit=unitNumber)
         call freeunit(unitNumber)
      end if
   end if
   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)
   if (ierr /= 0) then
      if (masterproc) then
         write(iulog,*) 'Failed to read namelist group &printout'
      end if
      return
   end if
   call mpi_bcast(mpas_print_global_minmax_vel, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_print_global_minmax_vel'
      end if
      return
   end if
   call mpi_bcast(mpas_print_detailed_minmax_vel, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_print_detailed_minmax_vel'
      end if
      return
   end if
   call mpi_bcast(mpas_print_global_minmax_sca, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)
   if (mpi_ierr /= mpi_success) then
      if (masterproc) then
         write(iulog,*) 'MPI_Bcast failed for namelist option mpas_print_global_minmax_sca'
      end if
      return
   end if

   call mpas_pool_add_config(configPool, 'config_print_global_minmax_vel', mpas_print_global_minmax_vel)
   call mpas_pool_add_config(configPool, 'config_print_detailed_minmax_vel', mpas_print_detailed_minmax_vel)
   call mpas_pool_add_config(configPool, 'config_print_global_minmax_sca', mpas_print_global_minmax_sca)

   if (masterproc) then
      close(unit=unitNumber)
      call freeunit(unitNumber)
   end if

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
      write(iulog,*) '   mpas_num_halos = ', mpas_num_halos
      write(iulog,*) '   mpas_zd = ', mpas_zd
      write(iulog,*) '   mpas_xnutr = ', mpas_xnutr
      write(iulog,*) '   mpas_block_decomp_file_prefix = ', trim(mpas_block_decomp_file_prefix)
      write(iulog,*) '   mpas_number_of_blocks = ', mpas_number_of_blocks
      write(iulog,*) '   mpas_explicit_proc_decomp = ', mpas_explicit_proc_decomp
      write(iulog,*) '   mpas_proc_decomp_file_prefix = ', trim(mpas_proc_decomp_file_prefix)
      write(iulog,*) '   mpas_do_restart = ', mpas_do_restart
      write(iulog,*) '   mpas_print_global_minmax_vel = ', mpas_print_global_minmax_vel
      write(iulog,*) '   mpas_print_detailed_minmax_vel = ', mpas_print_detailed_minmax_vel
      write(iulog,*) '   mpas_print_global_minmax_sca = ', mpas_print_global_minmax_sca
   end if

end function cam_mpas_namelist_read

end module dyn_comp
