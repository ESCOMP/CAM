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
   real(r8), dimension(: ),     pointer     :: phis   ! Surface geopotential      (ncol)
   real(r8), dimension(: ),     pointer     :: psd    ! Dry surface pressure      (ncol)
   real(r8), dimension(:,:,: ), pointer     :: uvperp ! Normal velocity at edges  (ncol,nedge,nver)
   real(r8), dimension(:,:   ), pointer     :: ux     ! Lon veloc at center       (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: uy     ! Lat veloc at center       (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: t      ! Temperature               (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: omega  ! Omega                     (ncol,nver+1)
   real(r8), dimension(:,:,: ), pointer     :: tracer ! Tracers                   (ncol,nver,nq)
   real(r8), dimension(:,:   ), pointer     :: ux_tend! Lon veloc tend at center  (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: uy_tend! Lat veloc tend at center  (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: t_tend ! Temperature tendency      (ncol,nver)
end type dyn_import_t

type dyn_export_t
   real(r8), dimension(: ),     pointer     :: phis   ! Surface geopotential      (ncol)
   real(r8), dimension(: ),     pointer     :: psd    ! Dry surface pressure      (ncol)
   real(r8), dimension(:,: ),   pointer     :: pint   ! Dry pressure at layer interfaces (ncol,nver+1)
   real(r8), dimension(:,: ),   pointer     :: pmid   ! Dry pressure at layer mid-points (ncol,nver)
   real(r8), dimension(:,: ),   pointer     :: zint   ! Geopotential height 
                                                       !              at layer interfaces (ncol,nver+1)
   real(r8), dimension(:,: ),   pointer     :: zmid   ! Geopotential height 
                                                       !              at layer mid-points (ncol,nver)
   real(r8), dimension(:,:,: ), pointer     :: uvperp ! Normal velocity at edges  (ncol,nedge,nver)
   real(r8), dimension(:,:   ), pointer     :: ux     ! Lon veloc at center       (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: uy     ! Lat veloc at center       (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: t      ! Temperature               (ncol,nver)
   real(r8), dimension(:,:   ), pointer     :: omega  ! Omega                     (ncol,nver+1)
   real(r8), dimension(:,:,: ), pointer     :: tracer ! Tracers                   (ncol,nver,nq)
   real(r8), dimension(:,:   ), pointer     :: pressure! Pressure                 (ncol,nver)
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

   type(mpas_pool_type), pointer :: dyn_inPool
   type(mpas_pool_type), pointer :: dyn_outPool

   !----------------------------------------------------------------------------

   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'dyn_in', dyn_inPool)
   call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'dyn_out', dyn_outPool)


   !
   ! Let dynamics import state point to memory managed by MPAS-Atmosphere
   !
   call mpas_pool_get_array(dyn_inPool, 'phis', dyn_in % phis)
   call mpas_pool_get_array(dyn_inPool, 'psd', dyn_in % psd)
   call mpas_pool_get_array(dyn_inPool, 'ux', dyn_in % ux)
   call mpas_pool_get_array(dyn_inPool, 'uy', dyn_in % uy)
   call mpas_pool_get_array(dyn_inPool, 't', dyn_in % t)
   call mpas_pool_get_array(dyn_inPool, 'omega', dyn_in % omega)
   call mpas_pool_get_array(dyn_inPool, 'tracer', dyn_in % tracer)
   call mpas_pool_get_array(dyn_inPool, 'ux_tend', dyn_in % ux_tend)
   call mpas_pool_get_array(dyn_inPool, 'uy_tend', dyn_in % uy_tend)
   call mpas_pool_get_array(dyn_inPool, 't_tend', dyn_in % t_tend)

   !
   ! Let dynamics export state point to memory managed by MPAS-Atmosphere
   !
   call mpas_pool_get_array(dyn_outPool, 'phis', dyn_out % phis)
   call mpas_pool_get_array(dyn_outPool, 'psd', dyn_out % psd)
   call mpas_pool_get_array(dyn_outPool, 'pint', dyn_out % pint)
   call mpas_pool_get_array(dyn_outPool, 'pmid', dyn_out % pmid)
   call mpas_pool_get_array(dyn_outPool, 'zint', dyn_out % zint)
   call mpas_pool_get_array(dyn_outPool, 'zmid', dyn_out % zmid)
   call mpas_pool_get_array(dyn_outPool, 'ux', dyn_out % ux)
   call mpas_pool_get_array(dyn_outPool, 'uy', dyn_out % uy)
   call mpas_pool_get_array(dyn_outPool, 't', dyn_out % t)
   call mpas_pool_get_array(dyn_outPool, 'omega', dyn_out % omega)
   call mpas_pool_get_array(dyn_outPool, 'tracer', dyn_out % tracer)
   call mpas_pool_get_array(dyn_outPool, 'pressure', dyn_out % pressure)

   call read_phis(dyn_in)

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
   nullify(dyn_in % phis)
   nullify(dyn_in % psd)
   nullify(dyn_in % ux)
   nullify(dyn_in % uy)
   nullify(dyn_in % t)
   nullify(dyn_in % omega)
   nullify(dyn_in % tracer)
   nullify(dyn_in % ux_tend)
   nullify(dyn_in % uy_tend)
   nullify(dyn_in % t_tend)

   !
   ! Prevent any further access to MPAS-Atmosphere memory
   !
   nullify(dyn_out % phis)
   nullify(dyn_out % psd)
   nullify(dyn_out % pint)
   nullify(dyn_out % pmid)
   nullify(dyn_out % zint)
   nullify(dyn_out % zmid)
   nullify(dyn_out % ux)
   nullify(dyn_out % uy)
   nullify(dyn_out % t)
   nullify(dyn_out % omega)
   nullify(dyn_out % tracer)
   nullify(dyn_out % pressure)

end subroutine dyn_final

!=========================================================================================
! Private routines.
!=========================================================================================

subroutine read_inidat(dyn_in)

   ! Set initial conditions.  Either from analytic expressions or read from file.

   ! arguments
   type(dyn_import_t), target, intent(inout) :: dyn_in

   ! Local variables

   type(file_desc_t), pointer :: fh_ini
   type(file_desc_t), pointer :: fh_topo

   integer :: p, nCellsLocal, io_master_id, m, i, k
   integer, allocatable, dimension(:) :: proc_start_i
   real(r8), allocatable, dimension(:,:) :: tmp_in, tmp_out


   real(r8), pointer :: psd(:)        !  psd(numcols)               ! dry surface pressure
   real(r8), pointer :: phis(:)       !  phis(numcols)              ! surface geopotential
   real(r8), pointer :: pt(:,:)       !  pt(numcols,plev)           ! potential temperature
   real(r8), pointer :: ux(:,:)       !  ux(numcols,plev)           ! cell centered velocity
   real(r8), pointer :: uy(:,:)       !  uy(numcols,plev)           ! cell centered velocity
   real(r8), pointer :: uvperp(:,:,:) !  uvperp(numcols,maxEdges,plev) ! edge normal velocity
   real(r8), pointer :: tracer(:,:,:) !  tracer(numcols,plev,pcnst) ! tracers

   real(r8), allocatable :: pdry(:,:)
   real(r8), allocatable :: delpp(:,:)
   real(r8), allocatable :: ptot(:,:)
   real(r8), allocatable :: pcap(:,:)

   character(len=*), parameter :: subname = 'dyn_comp::read_inidat'


   MPAS_DEBUG_WRITE(0, 'begin '//subname)

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

! For CAM/MPAS simulations, PS contains the dry surface pressure
! For CAM with fixed dynamics, PS contains the total pressure

   !MGD IO -- call MPAS_read_initial(ncid_ini, ncid_topo?) 
   !MGD IO -- fill in dyn_in at this point?

end subroutine read_inidat

!========================================================================================

subroutine read_phis(dyn_in)

   ! Set PHIS according to the following rules.
   !
   ! 1) If a topo file is specified use it.  This option has highest precedence.
   ! 2) If not using topo file, but analytic_ic option is on, use analytic phis.
   ! 3) Set phis = 0.0.
   !
   ! If using the physics grid then the topo file will be on that grid since its
   ! contents are primarily for the physics parameterizations, and the values of
   ! PHIS should be consistent with the values of sub-grid variability (e.g., SGH)
   ! which are computed on the physics grid.  In this case phis on the physics grid
   ! will be interpolated to the GLL grid.


   ! Arguments
   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

   ! local variables
   type(file_desc_t), pointer       :: fh_topo

   real(r8), allocatable            :: phis_tmp(:,:)      ! (npsp,nelemd)
   real(r8), allocatable            :: phis_phys_tmp(:,:) ! (fv_nphys**2,nelemd)

   integer                          :: i, ie, indx, j, kptr
   integer                          :: ierr, pio_errtype

   character(len=max_fieldname_len) :: fieldname
   character(len=max_hcoordname_len):: grid_name
   integer                          :: dims(2)
   integer                          :: ncells
   integer                          :: ncol_did
   integer                          :: ncol_size

   integer(iMap), pointer           :: gcid(:)            ! map local column order to global order
   logical,  allocatable            :: pmask(:)           ! unique columns

   ! Variables for analytic initial conditions
   integer,  allocatable            :: glob_ind(:)
   logical,  allocatable            :: pmask_phys(:)
   real(r8), pointer                :: latvals_deg(:)
   real(r8), pointer                :: lonvals_deg(:)
   real(r8), allocatable            :: latvals(:)
   real(r8), allocatable            :: lonvals(:)
   real(r8), allocatable            :: latvals_phys(:)
   real(r8), allocatable            :: lonvals_phys(:)

   character(len=*), parameter      :: subname='read_phis'
   !----------------------------------------------------------------------------

   fh_topo => topo_file_get_id()


   ! Set mask to indicate which columns are active in GLL grid.
   nullify(gcid)
   call cam_grid_get_gcid(cam_grid_id('mpas_cell'), gcid)
   allocate(pmask(size(gcid)))
   pmask(:) = (gcid /= 0)
   deallocate(gcid)

   if (associated(fh_topo)) then

      ! Set PIO to return error flags.
      call pio_seterrorhandling(fh_topo, PIO_BCAST_ERROR, pio_errtype)

      ! Get number of global columns from the grid object and check that
      ! it matches the file data.
      call cam_grid_dimensions('mpas_cell', dims)
      ncells = dims(1)

      ! code to read phis...

      ! Put the error handling back the way it was
      call pio_seterrorhandling(fh_topo, pio_errtype)

   else if (analytic_ic_active()) then

      ! lat/lon needed in radians
      latvals_deg => cam_grid_get_latvals(cam_grid_id('mpas_cell'))
      lonvals_deg => cam_grid_get_lonvals(cam_grid_id('mpas_cell'))
      allocate(latvals(size(latvals_deg)))
      allocate(lonvals(size(lonvals_deg)))
      latvals(:) = latvals_deg(:)*deg2rad
      lonvals(:) = lonvals_deg(:)*deg2rad


!      call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind, &
!                              PHIS=phis_tmp, mask=pmask(:))



   end if

   deallocate(pmask)


!   deallocate(phis_tmp)


end subroutine read_phis

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
