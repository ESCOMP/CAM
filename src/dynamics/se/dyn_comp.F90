module dyn_comp

! CAM interfaces to the SE Dynamical Core

use shr_kind_mod,           only: r8=>shr_kind_r8, shr_kind_cl
use physconst,              only: pi
use spmd_utils,             only: iam, masterproc
use constituents,           only: pcnst, cnst_get_ind, cnst_name, cnst_longname, &
                                  cnst_read_iv, qmin, cnst_type, tottnam,        &
                                  cnst_is_a_water_species
use cam_control_mod,        only: initial_run, simple_phys
use cam_initfiles,          only: initial_file_get_id, topo_file_get_id, pertlim
use phys_control,           only: use_gw_front, use_gw_front_igw, waccmx_is
use dyn_grid,               only: ini_grid_name, timelevel, hvcoord, edgebuf

use cam_grid_support,       only: cam_grid_id, cam_grid_get_gcid, &
                                  cam_grid_dimensions, cam_grid_get_dim_names, &
                                  cam_grid_get_latvals, cam_grid_get_lonvals,  &
                                  max_hcoordname_len
use cam_map_utils,          only: iMap

use inic_analytic,          only: analytic_ic_active, analytic_ic_set_ic
use dyn_tests_utils,        only: vcoord=>vc_dry_pressure

use cam_history,            only: outfld, hist_fld_active, fieldname_len
use cam_history_support,    only: max_fieldname_len
use time_manager,           only: get_step_size

use ncdio_atm,              only: infld
use pio,                    only: file_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR, &
                                  pio_inq_dimid, pio_inq_dimlen, PIO_NOERR

use infnan,                 only: isnan
use cam_logfile,            only: iulog
use cam_abortutils,         only: endrun
use shr_sys_mod,            only: shr_sys_flush

use parallel_mod,           only: par
use hybrid_mod,             only: hybrid_t
use dimensions_mod,         only: nelemd, nlev, np, npsq, ntrac, nc, fv_nphys, &
                                  qsize
use element_mod,            only: element_t, elem_state_t
use fvm_control_volume_mod, only: fvm_struct
use time_mod,               only: nsplit
use edge_mod,               only: initEdgeBuffer, edgeVpack, edgeVunpack, FreeEdgeBuffer
use edgetype_mod,           only: EdgeBuffer_t
use bndry_mod,              only: bndry_exchange

implicit none
private
save

public ::          &
     dyn_import_t, &
     dyn_export_t, &
     dyn_readnl,   &
     dyn_register, &
     dyn_init,     &
     dyn_run,      &
     dyn_final

type dyn_import_t
  type (element_t),  pointer :: elem(:) => null()
  type (fvm_struct), pointer :: fvm(:) => null()
end type dyn_import_t

type dyn_export_t
  type (element_t),  pointer :: elem(:) => null()
  type (fvm_struct), pointer :: fvm(:) => null()
end type dyn_export_t

! Namelist
logical, public, protected :: write_restart_unstruct

! Frontogenesis indices
integer, public    :: frontgf_idx      = -1
integer, public    :: frontga_idx      = -1

! constituent indices for waccm-x dry air properties
integer, public, protected :: &
   ixo  = -1, &
   ixo2 = -1, &
   ixh  = -1, &
   ixh2 = -1

interface read_dyn_var
  module procedure read_dyn_field_2d
  module procedure read_dyn_field_3d
end interface read_dyn_var

real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8

!===============================================================================
contains
!===============================================================================

subroutine dyn_readnl(NLFileName)
   use physconst,      only: thermodynamic_active_species_num
   use namelist_utils, only: find_group_name
   use namelist_mod,   only: homme_set_defaults, homme_postprocess_namelist
   use units,          only: getunit, freeunit
   use spmd_utils,     only: masterproc, masterprocid, mpicom, npes
   use spmd_utils,     only: mpi_real8, mpi_integer, mpi_character, mpi_logical
   use dyn_grid,       only: se_write_grid_file, se_grid_filename, se_write_gll_corners
   use dp_mapping,     only: nphys_pts
   use native_mapping, only: native_mapping_readnl

   use control_mod,    only: hypervis_subcycle, hypervis_subcycle_sponge
   use control_mod,    only: hypervis_subcycle_q, statefreq, runtype
   use control_mod,    only: nu, nu_div, nu_p, nu_q, nu_top, qsplit, rsplit
   use control_mod,    only: vert_remap_uvTq_alg, vert_remap_tracer_alg
   use control_mod,    only: tstep_type, rk_stage_user
   use control_mod,    only: ftype, limiter_option, partmethod
   use control_mod,    only: topology, phys_dyn_cp, variable_nsplit
   use control_mod,    only: fine_ne, hypervis_power, hypervis_scaling
   use control_mod,    only: max_hypervis_courant, statediag_numtrac,refined_mesh
   use control_mod,    only: molecular_diff
   use dimensions_mod, only: ne, npart
   use dimensions_mod, only: lcp_moist
   use dimensions_mod, only: hypervis_dynamic_ref_state,large_Courant_incr
   use dimensions_mod, only: fvm_supercycling, fvm_supercycling_jet
   use dimensions_mod, only: kmin_jet, kmax_jet
   use params_mod,     only: SFCURVE
   use parallel_mod,   only: initmpi
   use thread_mod,     only: initomp, max_num_threads
   use thread_mod,     only: horz_num_threads, vert_num_threads, tracer_num_threads
   use physconst,      only: rearth
   ! Dummy argument
   character(len=*), intent(in) :: NLFileName

   ! Local variables
   integer                      :: unitn, ierr,k
   real(r8)                     :: uniform_res_hypervis_scaling,nu_fac

   ! SE Namelist variables
   integer                      :: se_fine_ne
   integer                      :: se_ftype
   integer                      :: se_statediag_numtrac
   integer                      :: se_fv_nphys
   real(r8)                     :: se_hypervis_power
   real(r8)                     :: se_hypervis_scaling
   integer                      :: se_hypervis_subcycle
   integer                      :: se_hypervis_subcycle_sponge
   integer                      :: se_hypervis_subcycle_q
   integer                      :: se_limiter_option
   real(r8)                     :: se_max_hypervis_courant
   character(len=SHR_KIND_CL)   :: se_mesh_file
   integer                      :: se_ne
   integer                      :: se_npes
   integer                      :: se_nsplit
   real(r8)                     :: se_nu
   real(r8)                     :: se_nu_div
   real(r8)                     :: se_nu_p
   real(r8)                     :: se_nu_top
   integer                      :: se_qsplit
   logical                      :: se_refined_mesh
   integer                      :: se_rsplit
   integer                      :: se_statefreq
   integer                      :: se_tstep_type
   character(len=32)            :: se_vert_remap_T
   character(len=32)            :: se_vert_remap_uvTq_alg
   character(len=32)            :: se_vert_remap_tracer_alg
   integer                      :: se_horz_num_threads
   integer                      :: se_vert_num_threads
   integer                      :: se_tracer_num_threads
   logical                      :: se_hypervis_dynamic_ref_state
   logical                      :: se_lcp_moist
   logical                      :: se_write_restart_unstruct
   logical                      :: se_large_Courant_incr
   integer                      :: se_fvm_supercycling
   integer                      :: se_fvm_supercycling_jet
   integer                      :: se_kmin_jet
   integer                      :: se_kmax_jet
   integer                      :: se_phys_dyn_cp
   real(r8)                     :: se_molecular_diff

   namelist /dyn_se_inparm/        &
      se_fine_ne,                  & ! For refined meshes
      se_ftype,                    & ! forcing type
      se_statediag_numtrac,        &
      se_fv_nphys,                 &
      se_hypervis_power,           &
      se_hypervis_scaling,         &
      se_hypervis_subcycle,        &
      se_hypervis_subcycle_sponge, &
      se_hypervis_subcycle_q,      &
      se_limiter_option,           &
      se_max_hypervis_courant,     &
      se_mesh_file,                & ! Refined mesh definition file
      se_ne,                       &
      se_npes,                     &
      se_nsplit,                   & ! # of dyn steps per physics timestep
      se_nu,                       &
      se_nu_div,                   &
      se_nu_p,                     &
      se_nu_top,                   &
      se_qsplit,                   &
      se_refined_mesh,             &
      se_rsplit,                   &
      se_statefreq,                & ! number of steps per printstate call
      se_tstep_type,               &
      se_vert_remap_T,             &
      se_vert_remap_uvTq_alg,      &
      se_vert_remap_tracer_alg,    &
      se_write_grid_file,          &
      se_grid_filename,            &
      se_write_gll_corners,        &
      se_horz_num_threads,         &
      se_vert_num_threads,         &
      se_tracer_num_threads,       &
      se_hypervis_dynamic_ref_state,&
      se_lcp_moist,                &
      se_write_restart_unstruct,   &
      se_large_Courant_incr,       &
      se_fvm_supercycling,         &
      se_fvm_supercycling_jet,     &
      se_kmin_jet,                 &
      se_kmax_jet,                 &
      se_phys_dyn_cp,              &
      se_molecular_diff

   !--------------------------------------------------------------------------

   ! defaults for variables not set by build-namelist
   se_fine_ne                  = -1
   se_hypervis_power           = 0
   se_hypervis_scaling         = 0
   se_max_hypervis_courant     = 1.0e99_r8
   se_mesh_file                = ''
   se_npes                     = npes
   se_write_restart_unstruct   = .false.

   ! Read the namelist (dyn_se_inparm)
   call MPI_barrier(mpicom, ierr)
   if (masterproc) then
      write(iulog, *) "dyn_readnl: reading dyn_se_inparm namelist..."
      unitn = getunit()
      open( unitn, file=trim(NLFileName), status='old' )
      call find_group_name(unitn, 'dyn_se_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, dyn_se_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun('dyn_readnl: ERROR reading dyn_se_inparm namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist values to all PEs
   call MPI_bcast(se_fine_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_ftype, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_statediag_numtrac, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_hypervis_power, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_hypervis_scaling, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_hypervis_subcycle, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_hypervis_subcycle_sponge, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_hypervis_subcycle_q, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_limiter_option, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_max_hypervis_courant, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_mesh_file, SHR_KIND_CL,  mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_ne, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_nsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_nu, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_nu_div, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_nu_p, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_nu_top, 1, mpi_real8, masterprocid, mpicom, ierr)
   call MPI_bcast(se_qsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_refined_mesh, 1, mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_rsplit, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_statefreq, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_tstep_type, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_vert_remap_T, 32, mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_vert_remap_uvTq_alg, 32, mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_vert_remap_tracer_alg, 32, mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_fv_nphys, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_write_grid_file, 16,  mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_grid_filename, shr_kind_cl, mpi_character, masterprocid, mpicom, ierr)
   call MPI_bcast(se_write_gll_corners, 1,  mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_horz_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)
   call MPI_bcast(se_vert_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)
   call MPI_bcast(se_tracer_num_threads, 1, MPI_integer, masterprocid, mpicom,ierr)
   call MPI_bcast(se_hypervis_dynamic_ref_state, 1, mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_lcp_moist, 1, mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_write_restart_unstruct, 1, mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_large_Courant_incr, 1, mpi_logical, masterprocid, mpicom, ierr)
   call MPI_bcast(se_fvm_supercycling, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_fvm_supercycling_jet, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_kmin_jet, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_kmax_jet, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_phys_dyn_cp, 1, mpi_integer, masterprocid, mpicom, ierr)
   call MPI_bcast(se_molecular_diff, 1, mpi_real8, masterprocid, mpicom, ierr)

   if (se_npes <= 0) then
      call endrun('dyn_readnl: ERROR: se_npes must be > 0')
   end if

   ! Initialize the SE structure that holds the MPI decomposition information
   par = initmpi(se_npes)
   call initomp()


   if (se_fvm_supercycling < 0) se_fvm_supercycling = se_rsplit
   if (se_fvm_supercycling_jet < 0) se_fvm_supercycling_jet = se_rsplit

   ! Go ahead and enforce ne = 0 for refined mesh runs
   if (se_refined_mesh) then
      se_ne = 0
   end if

   ! Set HOMME defaults
   call homme_set_defaults()
   ! Set HOMME variables not in CAM's namelist but with different CAM defaults
   partmethod               = SFCURVE
   npart                    = se_npes
   ! CAM requires forward-in-time, subcycled dynamics
   ! RK2 3 stage tracers, sign-preserving conservative
   rk_stage_user            = 3
   topology                 = "cube"
   ! Finally, set the HOMME variables which have different names
   fine_ne                  = se_fine_ne
   ftype                    = se_ftype
   statediag_numtrac        = MIN(se_statediag_numtrac,pcnst)
   hypervis_power           = se_hypervis_power
   hypervis_scaling         = se_hypervis_scaling
   hypervis_subcycle        = se_hypervis_subcycle
   if (hypervis_subcycle_sponge<0) then
     hypervis_subcycle_sponge = hypervis_subcycle
   else
     hypervis_subcycle_sponge = se_hypervis_subcycle_sponge
   end if
   hypervis_subcycle_q      = se_hypervis_subcycle_q
   limiter_option           = se_limiter_option
   max_hypervis_courant     = se_max_hypervis_courant
   refined_mesh             = se_refined_mesh
   ne                       = se_ne
   nsplit                   = se_nsplit
   nu                       = se_nu
   nu_div                   = se_nu_div
   nu_p                     = se_nu_p
   nu_q                     = se_nu_p !for tracer-wind consistency nu_q must me equal to nu_p
   nu_top                   = se_nu_top
   qsplit                   = se_qsplit
   rsplit                   = se_rsplit
   statefreq                = se_statefreq
   tstep_type               = se_tstep_type
   vert_remap_uvTq_alg      = set_vert_remap(se_vert_remap_T, se_vert_remap_uvTq_alg)
   vert_remap_tracer_alg    = set_vert_remap(se_vert_remap_T, se_vert_remap_tracer_alg)
   fv_nphys                 = se_fv_nphys
   hypervis_dynamic_ref_state = se_hypervis_dynamic_ref_state
   lcp_moist                = se_lcp_moist
   large_Courant_incr       = se_large_Courant_incr
   fvm_supercycling         = se_fvm_supercycling
   fvm_supercycling_jet     = se_fvm_supercycling_jet
   kmin_jet                 = se_kmin_jet
   kmax_jet                 = se_kmax_jet
   variable_nsplit          = .false.
   phys_dyn_cp              = se_phys_dyn_cp
   molecular_diff           = se_molecular_diff

   if (fv_nphys > 0) then
      ! Use finite volume physics grid and CSLAM for tracer advection
      nphys_pts = fv_nphys*fv_nphys
      qsize = thermodynamic_active_species_num ! number tracers advected by GLL
      ntrac = pcnst                    ! number tracers advected by CSLAM
   else
      ! Use GLL grid for physics and tracer advection
      nphys_pts = npsq
      qsize = pcnst
      ntrac = 0
   end if

   if (rsplit < 1) then
      call endrun('dyn_readnl: rsplit must be > 0')
   end if

   ! if restart or branch run
   if (.not. initial_run) then
      runtype = 1
   end if

   ! HOMME wants 'none' to indicate no mesh file
   if (len_trim(se_mesh_file) == 0) then
      se_mesh_file = 'none'
      if (se_refined_mesh) then
         call endrun('dyn_readnl ERROR: se_refined_mesh=.true. but no se_mesh_file')
      end if
   end if
   call homme_postprocess_namelist(se_mesh_file, par)

   ! Set threading numbers to reasonable values
   if ((se_horz_num_threads == 0) .and. (se_vert_num_threads == 0) .and. (se_tracer_num_threads == 0)) then
      ! The user has not set any threading values, choose defaults
      se_horz_num_threads = 1
      se_vert_num_threads = max_num_threads
      se_tracer_num_threads = se_vert_num_threads
   end if
   if (se_horz_num_threads < 1) then
      se_horz_num_threads = 1
   end if
   if (se_vert_num_threads < 1) then
      se_vert_num_threads = 1
   end if
   if (se_tracer_num_threads < 1) then
      se_tracer_num_threads = 1
   end if
   horz_num_threads = se_horz_num_threads
   vert_num_threads = se_vert_num_threads
   tracer_num_threads = se_tracer_num_threads

   write_restart_unstruct = se_write_restart_unstruct

   if (se_kmin_jet<0            ) kmin_jet             = 1
   if (se_kmax_jet<0            ) kmax_jet             = nlev

   if (masterproc) then
      write(iulog, '(a,i0)')   'dyn_readnl: se_ftype                    = ',ftype
      write(iulog, '(a,i0)')   'dyn_readnl: se_statediag_numtrac        = ',statediag_numtrac
      write(iulog, '(a,i0)')   'dyn_readnl: se_hypervis_subcycle        = ',se_hypervis_subcycle
      write(iulog, '(a,i0)')   'dyn_readnl: se_hypervis_subcycle_sponge = ',se_hypervis_subcycle_sponge
      write(iulog, '(a,i0)')   'dyn_readnl: se_hypervis_subcycle_q      = ',se_hypervis_subcycle_q
      write(iulog, '(a,l4)')   'dyn_readnl: se_large_Courant_incr       = ',se_large_Courant_incr
      write(iulog, '(a,i0)')   'dyn_readnl: se_limiter_option           = ',se_limiter_option
      if (.not. se_refined_mesh) then
         write(iulog, '(a,i0)')'dyn_readnl: se_ne                       = ',se_ne
      end if
      write(iulog, '(a,i0)')   'dyn_readnl: se_npes                     = ',se_npes
      write(iulog, '(a,i0)')   'dyn_readnl: se_nsplit                   = ',se_nsplit
      write(iulog, '(a,i0)')   'dyn_readnl: se_phys_dyn_cp              = ',se_phys_dyn_cp
      !
      ! se_nu<0 then coefficients are set automatically in module global_norms_mod
      !
      if (se_nu_div>0) &
           write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu                       = ',se_nu
      if (se_nu_div>0) &
           write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_div                   = ',se_nu_div
      if (se_nu_p>0) then
        write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_p                     = ',se_nu_p
        write(iulog, '(a)') 'Note that nu_q must be the same as nu_p for  mass / tracer inconsistency'
      end if
      write(iulog, '(a,e9.2)') 'dyn_readnl: se_nu_top                     = ',se_nu_top
      write(iulog, '(a,i0)')   'dyn_readnl: se_qsplit                     = ',se_qsplit
      write(iulog, '(a,i0)')   'dyn_readnl: se_rsplit                     = ',se_rsplit
      write(iulog, '(a,i0)')   'dyn_readnl: se_statefreq                  = ',se_statefreq
      write(iulog, '(a,i0)')   'dyn_readnl: se_tstep_type                 = ',se_tstep_type
      write(iulog, '(a,a)')    'dyn_readnl: se_vert_remap_T               = ',trim(se_vert_remap_T)
      write(iulog, '(a,a)')    'dyn_readnl: se_vert_remap_uvTq_alg        = ',trim(se_vert_remap_uvTq_alg)
      write(iulog, '(a,a)')    'dyn_readnl: se_vert_remap_tracer_alg      = ',trim(se_vert_remap_tracer_alg)
      write(iulog, '(a,l4)')   'dyn_readnl: se_hypervis_dynamic_ref_state = ',hypervis_dynamic_ref_state
      write(iulog, '(a,l4)')   'dyn_readnl: lcp_moist                     = ',lcp_moist
      write(iulog, '(a,i0)')   'dyn_readnl: se_fvm_supercycling           = ',fvm_supercycling
      write(iulog, '(a,i0)')   'dyn_readnl: se_fvm_supercycling_jet       = ',fvm_supercycling_jet
      write(iulog, '(a,i0)')   'dyn_readnl: se_kmin_jet                   = ',kmin_jet
      write(iulog, '(a,i0)')   'dyn_readnl: se_kmax_jet                   = ',kmax_jet
      if (se_refined_mesh) then
         write(iulog, '(a)') 'dyn_readnl: Refined mesh simulation'
         write(iulog, '(a)') 'dyn_readnl: se_mesh_file = ',trim(se_mesh_file)
         if (hypervis_power /= 0) then
           write(iulog, '(a)') 'Using scalar viscosity (Zarzycki et al 2014 JClim)'
           write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_power = ',se_hypervis_power, ', (tensor hyperviscosity)'
           write(iulog, '(a,e11.4)') 'dyn_readnl: se_max_hypervis_courant = ',se_max_hypervis_courant
         end if
         if (hypervis_scaling /= 0) then
           write(iulog, '(a)') 'Using tensor viscosity (Guba et al., 2014)'
           write(iulog, '(a,e11.4)') 'dyn_readnl: se_hypervis_scaling = ',se_hypervis_scaling
         end if
      end if

      if (fv_nphys > 0) then
         write(iulog, '(a)') 'dyn_readnl: physics will run on FVM points; advection by CSLAM'
         write(iulog,'(a,i0)') 'dyn_readnl: se_fv_nphys = ', fv_nphys
      else
         write(iulog, '(a)') 'dyn_readnl: physics will run on SE GLL points'
      end if
      write(iulog, '(a,i0)') 'dyn_readnl: se_horz_num_threads = ',horz_num_threads
      write(iulog, '(a,i0)') 'dyn_readnl: se_vert_num_threads = ',vert_num_threads
      write(iulog, '(a,i0)') 'dyn_readnl: se_tracer_num_threads = ',tracer_num_threads
      if (trim(se_write_grid_file) == 'SCRIP') then
         write(iulog,'(2a)') "dyn_readnl: write SCRIP grid file = ", trim(se_grid_filename)
      else
         write(iulog,'(a)') "dyn_readnl: do not write grid file"
      end if
      write(iulog,'(a,l1)') 'dyn_readnl: write gll corners to SEMapping.nc = ', &
                            se_write_gll_corners
      write(iulog,'(a,l1)') 'dyn_readnl: write restart data on unstructured grid = ', &
                            se_write_restart_unstruct

      write(iulog, '(a,e9.2)') 'dyn_readnl: se_molecular_diff  = ', molecular_diff
   end if

   call native_mapping_readnl(NLFileName)

   !---------------------------------------------------------------------------
   contains
   !---------------------------------------------------------------------------

      integer function set_vert_remap( remap_T, remap_alg )

         ! Convert namelist input strings to the internally used integers.

         character(len=*), intent(in) :: remap_T    ! scheme for remapping temperature
         character(len=*), intent(in) :: remap_alg  ! remapping algorithm

         ! check valid remap_T values:
         if (remap_T /= 'thermal_energy_over_P' .and. remap_T /= 'Tv_over_logP') then
            write(iulog,*)'set_vert_remap: invalid remap_T= ',trim(remap_T)
            call endrun('set_vert_remap: invalid remap_T')
         end if

         select case (remap_alg)
         case ('PPM_bc_mirror')
            set_vert_remap = 1
         case ('PPM_bc_PCoM')
            set_vert_remap = 2
         case ('PPM_bc_linear_extrapolation')
            set_vert_remap = 10
         case ('FV3_PPM')
            if (remap_T == 'thermal_energy_over_P') then
               set_vert_remap = -4
            else
               set_vert_remap = -40
            end if
         case ('FV3_CS')
            if (remap_T == 'thermal_energy_over_P') then
               set_vert_remap = -9
            else
               set_vert_remap = -90
            end if
         case ('FV3_CS_2dz_filter')
            if (remap_T == 'thermal_energy_over_P') then
               set_vert_remap = -10
            else
               set_vert_remap = -100
            end if
         case ('FV3_non_monotone_CS_2dz_filter')
            if (remap_T == 'thermal_energy_over_P') then
               set_vert_remap = -11
            else
               set_vert_remap = -110
            end if
         case default
            write(iulog,*)'set_vert_remap: invalid remap_alg= ',trim(remap_alg)
            call endrun('set_vert_remap: invalid remap_alg')
         end select

      end function set_vert_remap

end subroutine dyn_readnl

!=========================================================================================

subroutine dyn_register()

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

   if (use_gw_front .or. use_gw_front_igw) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/),       &
         frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/),       &
         frontga_idx)
   end if

end subroutine dyn_register

!=========================================================================================

subroutine dyn_init(dyn_in, dyn_out)
   use prim_advance_mod,   only: prim_advance_init
   use dyn_grid,           only: elem, fvm
   use cam_pio_utils,      only: clean_iodesc_list
   use physconst,          only: thermodynamic_active_species_num, thermodynamic_active_species_idx
   use physconst,          only: thermodynamic_active_species_idx_dycore, rair, cpair
   use cam_history,        only: addfld, add_default, horiz_only, register_vector_field
   use gravity_waves_sources, only: gws_init

   use thread_mod,         only: horz_num_threads
   use hybrid_mod,         only: get_loop_ranges, config_thread_region
   use dimensions_mod,     only: nu_scale_top, nu_lev, nu_div_lev
   use dimensions_mod,     only: ksponge_end, kmvis_ref, kmcnd_ref,rho_ref,km_sponge_factor
   use dimensions_mod,     only: cnst_name_gll, cnst_longname_gll
   use dimensions_mod,     only: irecons_tracer_lev, irecons_tracer, kord_tr, kord_tr_cslam
   use prim_driver_mod,    only: prim_init2
   use time_mod,           only: time_at
   use control_mod,        only: runtype, molecular_diff, nu_top
   use test_fvm_mapping,   only: test_mapping_addfld
   use phys_control,       only: phys_getopts
   use physconst,          only: get_molecular_diff_coef_reference
   use control_mod,        only: vert_remap_uvTq_alg, vert_remap_tracer_alg
   use std_atm_profile,    only: std_atm_height
   ! Dummy arguments:
   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out

   ! Local variables
   integer             :: ithr, nets, nete, ie, k, kmol_end
   real(r8), parameter :: Tinit = 300.0_r8
   real(r8)            :: press(1), ptop, tref,z(1)

   type(hybrid_t)      :: hybrid

   integer :: ixcldice, ixcldliq, ixrain, ixsnow, ixgraupel
   integer :: m_cnst, m

   ! variables for initializing energy and axial angular momentum diagnostics
   character (len = 3), dimension(12) :: stage = (/"dED","dAF","dBD","dAD","dAR","dBF","dBH","dCH","dAH",'dBS','dAS','p2d'/)
   character (len = 70),dimension(12) :: stage_txt = (/&
      " end of previous dynamics                           ",& !dED
      " from previous remapping or state passed to dynamics",& !dAF - state in beginning of nsplit loop
      " state after applying CAM forcing                   ",& !dBD - state after applyCAMforcing
      " before vertical remapping                          ",& !dAD - state before vertical remapping
      " after vertical remapping                           ",& !dAR - state at end of nsplit loop
      " state passed to parameterizations                  ",& !dBF
      " state before hypervis                              ",& !dBH
      " state after hypervis but before adding heating term",& !dCH
      " state after hypervis                               ",& !dAH
      " state before sponge layer diffusion                ",& !dBS - state before sponge del2
      " state after sponge layer diffusion                 ",& !dAS - state after sponge del2
      " phys2dyn mapping errors (requires ftype-1)         " & !p2d - for assessing phys2dyn mapping errors
      /)
   character (len = 2)  , dimension(8) :: vars  = (/"WV"  ,"WL"  ,"WI"  ,"SE"   ,"KE"   ,"MR"   ,"MO"   ,"TT"   /)
   !if ntrac>0 then tracers should be output on fvm grid but not energy (SE+KE) and AAM diags
   logical              , dimension(8) :: massv = (/.true.,.true.,.true.,.false.,.false.,.false.,.false.,.false./)
   character (len = 70) , dimension(8) :: vars_descriptor = (/&
      "Total column water vapor                ",&
      "Total column cloud water                ",&
      "Total column cloud ice                  ",&
      "Total column dry static energy          ",&
      "Total column kinetic energy             ",&
      "Total column wind axial angular momentum",&
      "Total column mass axial angular momentum",&
      "Total column test tracer                "/)
   character (len = 14), dimension(8)  :: &
      vars_unit = (/&
      "kg/m2        ","kg/m2        ","kg/m2        ","J/m2         ",&
      "J/m2         ","kg*m2/s*rad2 ","kg*m2/s*rad2 ","kg/m2        "/)

   integer :: istage, ivars
   character (len=108) :: str1, str2, str3

   logical :: history_budget      ! output tendencies and state variables for budgets
   integer :: budget_hfile_num

   character(len=*), parameter :: subname = 'dyn_init'

   real(r8) :: km_sponge_factor_local(nlev+1)
   !----------------------------------------------------------------------------

   ! Now allocate and set condenstate vars
   allocate(cnst_name_gll(qsize))     ! constituent names for gll tracers
   allocate(cnst_longname_gll(qsize)) ! long name of constituents for gll tracers

   allocate(kord_tr(qsize))
   kord_tr(:) = vert_remap_tracer_alg
   if (ntrac>0) then
     allocate(kord_tr_cslam(ntrac))
     kord_tr_cslam(:) = vert_remap_tracer_alg
   end if


   do m=1,qsize
     !
     ! The "_gll" index variables below are used to keep track of condensate-loading tracers
     ! since they are not necessarily indexed contiguously and not necessarily in the same
     ! order (physics is in charge of the order)
     !
     ! if running with CSLAM then the SE (gll) condensate-loading water tracers are always
     ! indexed contiguously (q,cldliq,cldice,rain,snow,graupel) - see above
     !
     ! CSLAM tracers are always indexed as in physics
     ! of no CSLAM then SE tracers are always indexed as in physics
     !
     if (ntrac>0) then
       !
       ! note that in this case qsize = thermodynamic_active_species_num
       !
       thermodynamic_active_species_idx_dycore(m) = m
       kord_tr_cslam(thermodynamic_active_species_idx(m)) = vert_remap_uvTq_alg
       kord_tr(m)                                 = vert_remap_uvTq_alg
       cnst_name_gll    (m)                       = cnst_name    (thermodynamic_active_species_idx(m))
       cnst_longname_gll(m)                       = cnst_longname(thermodynamic_active_species_idx(m))
     else
       !
       ! if not running with CSLAM then the condensate-loading water tracers are not necessarily
       ! indexed contiguously (are indexed as in physics)
       !
       if (m.le.thermodynamic_active_species_num) then
         thermodynamic_active_species_idx_dycore(m) = thermodynamic_active_species_idx(m)
         kord_tr(thermodynamic_active_species_idx_dycore(m)) = vert_remap_uvTq_alg
       end if
       cnst_name_gll    (m)                = cnst_name    (m)
       cnst_longname_gll(m)                = cnst_longname(m)

     end if
   end do



   !
   ! Initialize the import/export objects
   !
   if(iam < par%nprocs) then
     dyn_in%elem  => elem
     dyn_in%fvm   => fvm

     dyn_out%elem => elem
     dyn_out%fvm  => fvm
   else
     nullify(dyn_in%elem)
     nullify(dyn_in%fvm)
     nullify(dyn_out%elem)
     nullify(dyn_out%fvm)
   end if

   call set_phis(dyn_in)

   if (initial_run) then
     call read_inidat(dyn_in)
     call clean_iodesc_list()
   end if
   !
   ! initialize diffusion in dycore
   !
   kmol_end = 0
   if (molecular_diff>0) then
     !
     ! molecular diffusion and thermal conductivity reference values
     !
     if (masterproc) write(iulog,*) subname//": initialize molecular diffusion reference profiles"
     tref = 1000._r8     !mean value at model top for solar max
     km_sponge_factor = molecular_diff
     km_sponge_factor_local = molecular_diff
     !
     ! get rho, kmvis and kmcnd at mid-levels
     !
     call get_molecular_diff_coef_reference(1,nlev,tref,&
          (hvcoord%hyam(:)+hvcoord%hybm(:))*hvcoord%ps0,km_sponge_factor,&
          kmvis_ref,kmcnd_ref,rho_ref)

     do k=1,nlev
       ! only apply molecular viscosity where viscosity is > 1000 m/s^2
       if (MIN(kmvis_ref(k)/rho_ref(k),kmcnd_ref(k)/(cpair*rho_ref(k)))>1000.0_r8) then
         if (masterproc) then
           press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
           call std_atm_height(press,z)
           write(iulog,'(a,i3,3e11.4)') "k, p, z, km_sponge_factor                   :",k, &
                press, z,km_sponge_factor(k)
            write(iulog,'(a,2e11.4)') "kmvis_ref/rho_ref, kmcnd_ref/(cp*rho_ref): ", &
               kmvis_ref(k)/rho_ref(k),kmcnd_ref(k)/(cpair*rho_ref(k))
         end if
         kmol_end = k
       else
         kmvis_ref(k) = 1.0_r8
         kmcnd_ref(k) = 1.0_r8
       end if
     end do
   else
     ! -1.0E6 is an arbitrary unrealistic value.  But it is used in the calculation
     ! of a diagnostic quantity in global_norms_mod so can't be set to huge or nan.
     kmvis_ref(:) = -1.0E6_r8
     kmcnd_ref(:) = -1.0E6_r8
     rho_ref(:)   = -1.0E6_r8
   end if
   !
   irecons_tracer_lev(:) = irecons_tracer !use high-order CSLAM in all layers
   !
   ! compute scaling of traditional sponge layer damping (following cd_core.F90 in CAM-FV)
   !
   nu_scale_top(:) = 0.0_r8
   if (nu_top>0) then
     if (masterproc) write(iulog,*) subname//": sponge layer viscosity scaling factor"
     do k=1,nlev
       press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
       ptop  = hvcoord%hyai(1)*hvcoord%ps0
       nu_scale_top(k) = 8.0_r8*(1.0_r8+tanh(1.0_r8*log(ptop/press(1)))) ! tau will be maximum 8 at model top
       if (nu_scale_top(k).ge.0.15_r8) then
         ksponge_end = k
       else
         nu_scale_top(k) = 0.0_r8
       end if
     end do
   else
     ksponge_end = 0
   end if
   ksponge_end = MAX(MAX(ksponge_end,1),kmol_end)
   if (masterproc) then
     write(iulog,*) subname//": ksponge_end = ",ksponge_end
     if (nu_top>0) then
       do k=1,ksponge_end+1
         press = (hvcoord%hyam(k)+hvcoord%hybm(k))*hvcoord%ps0
         call std_atm_height(press,z)
         write(iulog,'(a,i3,4e11.4)') subname//": k, p, z, nu_scale_top, nu ",k,press,z,&
              nu_scale_top(k),nu_scale_top(k)*nu_top
       end do
     end if
   end if

   if (iam < par%nprocs) then
      call prim_advance_init(par,elem)
      !$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete)
      hybrid = config_thread_region(par,'horizontal')
      call get_loop_ranges(hybrid, ibeg=nets, iend=nete)
      call prim_init2(elem, fvm, hybrid, nets, nete, TimeLevel, hvcoord)
      !$OMP END PARALLEL

      if (use_gw_front .or. use_gw_front_igw) call gws_init(elem)
   end if  ! iam < par%nprocs

   call addfld ('nu_kmvis',   (/ 'lev' /), 'A', '', 'Molecular viscosity Laplacian coefficient'            , gridname='GLL')
   call addfld ('nu_kmcnd',   (/ 'lev' /), 'A', '', 'Thermal conductivity Laplacian coefficient'           , gridname='GLL')
   call addfld ('nu_kmcnd_dp',(/ 'lev' /), 'A', '', 'Thermal conductivity like Laplacian coefficient on dp', gridname='GLL')


   ! Forcing from physics on the GLL grid
   call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term on GLL grid',     gridname='GLL')
   call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term on GLL grid',gridname='GLL')
   call register_vector_field('FU', 'FV')
   call addfld ('FT',  (/ 'lev' /), 'A', 'K/s', 'Temperature forcing term on GLL grid',gridname='GLL')

   ! Tracer forcing on fvm (CSLAM) grid and internal CSLAM pressure fields
   if (ntrac>0) then
      do m = 1, ntrac
         call addfld (trim(cnst_name(m))//'_fvm',  (/ 'lev' /), 'I', 'kg/kg',   &
            trim(cnst_longname(m)), gridname='FVM')

         call addfld ('F'//trim(cnst_name(m))//'_fvm',  (/ 'lev' /), 'I', 'kg/kg/s',   &
            trim(cnst_longname(m))//' mixing ratio forcing term (q_new-q_old) on fvm grid', &
            gridname='FVM')
      end do

      call addfld ('dp_fvm' ,(/ 'lev' /), 'I', 'Pa','CSLAM Pressure level thickness', gridname='FVM')
      call addfld ('PSDRY_fvm',horiz_only, 'I','Pa','CSLAM dry surface pressure'    , gridname='FVM')
   end if

   do m_cnst = 1, qsize
     call addfld ('F'//trim(cnst_name_gll(m_cnst))//'_gll',  (/ 'lev' /), 'I', 'kg/kg/s',   &
          trim(cnst_longname(m_cnst))//' mixing ratio forcing term (q_new-q_old) on GLL grid', gridname='GLL')
   end do

   ! Energy diagnostics and axial angular momentum diagnostics
   call addfld ('ABS_dPSdt',  horiz_only, 'A', 'Pa/s', 'Absolute surface pressure tendency',gridname='GLL')

   if (ntrac>0) then
#ifdef waccm_debug
     call addfld ('CSLAM_gamma',  (/ 'lev' /), 'A', '', 'Courant number from CSLAM',     gridname='FVM')
#endif
     call addfld ('WV_PDC',   horiz_only, 'A', 'kg/m2','Total column water vapor lost in physics-dynamics coupling',gridname='FVM')
     call addfld ('WL_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud water lost in physics-dynamics coupling',gridname='FVM')
     call addfld ('WI_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud ice lost in physics-dynamics coupling'  ,gridname='FVM')
     call addfld ('TT_PDC',   horiz_only, 'A', 'kg/m2','Total column test tracer lost in physics-dynamics coupling',gridname='FVM')
   else
     call addfld ('WV_PDC',   horiz_only, 'A', 'kg/m2','Total column water vapor lost in physics-dynamics coupling',gridname='GLL')
     call addfld ('WL_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud water lost in physics-dynamics coupling',gridname='GLL')
     call addfld ('WI_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud ice lost in physics-dynamics coupling'  ,gridname='GLL')
     call addfld ('TT_PDC',   horiz_only, 'A', 'kg/m2','Total column test tracer lost in physics-dynamics coupling',gridname='GLL')
   end if

   do istage = 1,SIZE(stage)
      do ivars=1,SIZE(vars)
         write(str1,*) TRIM(ADJUSTL(vars(ivars))),TRIM(ADJUSTL("_")),TRIM(ADJUSTL(stage(istage)))
         write(str2,*) TRIM(ADJUSTL(vars_descriptor(ivars))),&
             TRIM(ADJUSTL(" ")),TRIM(ADJUSTL(stage_txt(istage)))
         write(str3,*) TRIM(ADJUSTL(vars_unit(ivars)))
         if (ntrac>0.and.massv(ivars)) then
           call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='FVM')
         else
           call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='GLL')
         end if
      end do
   end do

   !
   ! add dynamical core tracer tendency output
   !
   if (ntrac>0) then
     do m = 1, pcnst
       call addfld(tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert',  &
            gridname='FVM')
     end do
   else
     do m = 1, pcnst
       call addfld(tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert',  &
            gridname='GLL')
     end do
   end if
   call phys_getopts(history_budget_out=history_budget, history_budget_histfile_num_out=budget_hfile_num)
   if ( history_budget ) then
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      call add_default(tottnam(       1), budget_hfile_num, ' ')
      call add_default(tottnam(ixcldliq), budget_hfile_num, ' ')
      call add_default(tottnam(ixcldice), budget_hfile_num, ' ')
   end if

   ! constituent indices for waccm-x
   if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      call cnst_get_ind('O',  ixo)
      call cnst_get_ind('O2', ixo2)
      call cnst_get_ind('H',  ixh)
      call cnst_get_ind('H2', ixh2)
   end if

   call test_mapping_addfld
end subroutine dyn_init

!=========================================================================================

subroutine dyn_run(dyn_state)
   use physconst,        only: thermodynamic_active_species_num, dry_air_species_num
   use physconst,        only: thermodynamic_active_species_idx_dycore
   use prim_advance_mod, only: calc_tot_energy_dynamics
   use prim_driver_mod,  only: prim_run_subcycle
   use dimensions_mod,   only: cnst_name_gll
   use time_mod,         only: tstep, nsplit, timelevel_qdp
   use hybrid_mod,       only: config_thread_region, get_loop_ranges
   use control_mod,      only: qsplit, rsplit, ftype_conserve
   use thread_mod,       only: horz_num_threads
   use time_mod,         only: tevolve

   type(dyn_export_t), intent(inout) :: dyn_state

   type(hybrid_t) :: hybrid
   integer        :: tl_f
   integer        :: n
   integer        :: nets, nete, ithr
   integer        :: i, ie, j, k, m, nq, m_cnst
   integer        :: n0_qdp, nsplit_local
   logical        :: ldiag

   real(r8) :: ftmp(npsq,nlev,3)
   real(r8) :: dtime
   real(r8) :: rec2dt, pdel

   real(r8), allocatable, dimension(:,:,:) :: ps_before
   real(r8), allocatable, dimension(:,:,:) :: abs_ps_tend
   real (kind=r8)                          :: omega_cn(2,nelemd) !min and max of vertical Courant number
   !----------------------------------------------------------------------------

#ifdef debug_coupling
   return
#endif

   nsplit_local = nsplit
   tevolve = 0._r8

   if (iam >= par%nprocs) return

   ldiag = hist_fld_active('ABS_dPSdt')
   if (ldiag) then
      allocate(ps_before(np,np,nelemd))
      allocate(abs_ps_tend(np,np,nelemd))

   end if

   !$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n,ie,m,i,j,k,ftmp)
   hybrid = config_thread_region(par,'horizontal')
   call get_loop_ranges(hybrid, ibeg=nets, iend=nete)

   dtime = get_step_size()
   rec2dt = 1._r8/dtime

   tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics
   call TimeLevel_Qdp(TimeLevel, qsplit, n0_qdp)!get n0_qdp for diagnostics call

   ! output physics forcing
   if (hist_fld_active('FU') .or. hist_fld_active('FV') .or.hist_fld_active('FT')) then
      do ie = nets, nete
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  ftmp(i+(j-1)*np,k,1) = dyn_state%elem(ie)%derived%FM(i,j,1,k)
                  ftmp(i+(j-1)*np,k,2) = dyn_state%elem(ie)%derived%FM(i,j,2,k)
                  ftmp(i+(j-1)*np,k,3) = dyn_state%elem(ie)%derived%FT(i,j,k)
               end do
            end do
         end do

         call outfld('FU', ftmp(:,:,1), npsq, ie)
         call outfld('FV', ftmp(:,:,2), npsq, ie)
         call outfld('FT', ftmp(:,:,3), npsq, ie)
      end do
   end if

   do m = 1, qsize
     if (hist_fld_active('F'//trim(cnst_name_gll(m))//'_gll')) then
       do ie = nets, nete
         call outfld('F'//trim(cnst_name_gll(m))//'_gll',&
              RESHAPE(dyn_state%elem(ie)%derived%FQ(:,:,:,m), (/np*np,nlev/)), npsq, ie)
       end do
     end if
   end do



   ! convert elem(ie)%derived%fq to mass tendency
   do ie = nets, nete
      do m = 1, qsize
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  dyn_state%elem(ie)%derived%FQ(i,j,k,m) = dyn_state%elem(ie)%derived%FQ(i,j,k,m)* &
                     rec2dt*dyn_state%elem(ie)%state%dp3d(i,j,k,tl_f)
               end do
            end do
         end do
      end do
   end do


   if (ftype_conserve>0) then
     do ie = nets, nete
       do k=1,nlev
         do j=1,np
           do i = 1, np
             pdel     = dyn_state%elem(ie)%state%dp3d(i,j,k,tl_f)
             do nq=dry_air_species_num+1,thermodynamic_active_species_num
               m_cnst = thermodynamic_active_species_idx_dycore(nq)
               pdel = pdel + (dyn_state%elem(ie)%state%qdp(i,j,k,m_cnst,n0_qdp)+dyn_state%elem(ie)%derived%FQ(i,j,k,m_cnst)*dtime)
             end do
             dyn_state%elem(ie)%derived%FDP(i,j,k) = pdel
           end do
         end do
       end do
     end do
   end if


   if (ntrac > 0) then
      do ie = nets, nete
         do m = 1, ntrac
            do k = 1, nlev
               do j = 1, nc
                  do i = 1, nc
                     dyn_state%fvm(ie)%fc(i,j,k,m) = dyn_state%fvm(ie)%fc(i,j,k,m)* &
                        rec2dt!*dyn_state%fvm(ie)%dp_fvm(i,j,k)
                  end do
               end do
            end do
         end do
      end do
   end if

   if (ldiag) then
      abs_ps_tend(:,:,nets:nete) = 0.0_r8
   endif

   do n = 1, nsplit_local

      if (ldiag) then
         do ie = nets, nete
            ps_before(:,:,ie) = dyn_state%elem(ie)%state%psdry(:,:)
         end do
      end if

      ! forward-in-time RK, with subcycling
      call prim_run_subcycle(dyn_state%elem, dyn_state%fvm, hybrid, nets, nete, &
                             tstep, TimeLevel, hvcoord, n, omega_cn)

      if (ldiag) then
         do ie = nets, nete
            abs_ps_tend(:,:,ie) = abs_ps_tend(:,:,ie) +                                &
               ABS(ps_before(:,:,ie)-dyn_state%elem(ie)%state%psdry(:,:)) &
               /(tstep*qsplit*rsplit)
         end do
      end if

   end do

   if (ldiag) then
      do ie=nets,nete
         abs_ps_tend(:,:,ie)=abs_ps_tend(:,:,ie)/DBLE(nsplit)
         call outfld('ABS_dPSdt',RESHAPE(abs_ps_tend(:,:,ie),(/npsq/)),npsq,ie)
      end do
   end if


   call calc_tot_energy_dynamics(dyn_state%elem,dyn_state%fvm, nets, nete, TimeLevel%n0, n0_qdp,'dBF')
   !$OMP END PARALLEL

   if (ldiag) then
      deallocate(ps_before,abs_ps_tend)
   endif
   ! output vars on CSLAM fvm grid
   call write_dyn_vars(dyn_state)

end subroutine dyn_run

!===============================================================================

subroutine dyn_final(DYN_STATE, RESTART_FILE)

   type (elem_state_t), target     :: DYN_STATE
   character(LEN=*)   , intent(IN) :: RESTART_FILE

end subroutine dyn_final

!===============================================================================

subroutine read_inidat(dyn_in)
   use physconst,           only: thermodynamic_active_species_num, dry_air_species_num
   use shr_sys_mod,         only: shr_sys_flush
   use hycoef,              only: hyai, hybi, ps0
   use const_init,          only: cnst_init_default

   use element_mod,         only: timelevels
   use fvm_mapping,         only: dyn2fvm_mass_vars
   use control_mod,         only: runtype,initial_global_ave_dry_ps
   use prim_driver_mod,     only: prim_set_dry_mass
   use physconst,           only: thermodynamic_active_species_idx

   ! Arguments
   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

   ! Local variables

   integer(iMap), pointer           :: ldof(:) ! Basic (2D) grid dof

   type(file_desc_t), pointer       :: fh_ini, fh_topo

   type(element_t), pointer         :: elem(:)

   real(r8), allocatable            :: qtmp(:,:,:,:,:)    ! (np,np,nlev,nelemd,n)
   real(r8), allocatable            :: dbuf2(:,:)         ! (npsq,nelemd)
   real(r8), allocatable            :: dbuf3(:,:,:)       ! (npsq,nlev,nelemd)
   real(r8), allocatable            :: phis_tmp(:,:)      ! (npsp,nelemd)
   real(r8), allocatable            :: factor_array(:,:,:,:) ! (np,np,nlev,nelemd)
   logical,  allocatable            :: pmask(:)           ! (npsq*nelemd) unique grid vals

   character(len=max_hcoordname_len):: grid_name
   real(r8), allocatable            :: latvals(:),latvals_phys(:)
   real(r8), allocatable            :: lonvals(:),lonvals_phys(:)
   real(r8), pointer                :: latvals_deg(:)
   real(r8), pointer                :: lonvals_deg(:)

   integer                          :: ie, k, t
   character(len=max_fieldname_len) :: fieldname, fieldname2
   logical                          :: found
   logical                          :: inic_wet           ! true if initial condition is based on
                                                          ! wet pressure and water species
   integer                          :: kptr, m_cnst
   type(EdgeBuffer_t)               :: edge

   character(len=max_fieldname_len) :: dimname, varname
   integer                          :: ierr

   integer                          :: rndm_seed_sz
   integer, allocatable             :: rndm_seed(:)
   integer                          :: dims(2)
   integer                          :: pio_errtype
   real(r8)                         :: pertval
   integer                          :: i, j, indx, nq
   integer                          :: dyn_cols
   character(len=128)               :: errmsg
   character(len=*), parameter      :: subname='READ_INIDAT'

   ! fvm vars
   real(r8), allocatable            :: inv_dp_darea_fvm(:,:,:)
   real(r8)                         :: min_val, max_val

   real(r8)                         :: dp_tmp, pstmp(np,np)

   ! Variables for analytic initial conditions
   integer,  allocatable            :: glob_ind(:)
   integer,  allocatable            :: m_ind(:)
   real(r8), allocatable            :: dbuf4(:,:,:,:)
   !----------------------------------------------------------------------------

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

   if (iam < par%nprocs) then
      elem => dyn_in%elem
   else
      nullify(elem)
   end if

   allocate(qtmp(np,np,nlev,nelemd,pcnst))
   qtmp = 0._r8

   ! Set mask to indicate which columns are active
   nullify(ldof)
   call cam_grid_get_gcid(cam_grid_id(ini_grid_name), ldof)
   allocate(pmask(npsq*nelemd))
   pmask(:) = (ldof /= 0)

   ! lat/lon needed in radians
   latvals_deg => cam_grid_get_latvals(cam_grid_id(ini_grid_name))
   lonvals_deg => cam_grid_get_lonvals(cam_grid_id(ini_grid_name))
   allocate(latvals(np*np*nelemd))
   allocate(lonvals(np*np*nelemd))
   latvals(:) = latvals_deg(:)*deg2rad
   lonvals(:) = lonvals_deg(:)*deg2rad

   ! Set PIO to return error codes when reading data from IC file.
   call pio_seterrorhandling(fh_ini, PIO_BCAST_ERROR, pio_errtype)

   ! Get the number of columns in the global GLL grid.
   call cam_grid_dimensions(ini_grid_name, dims)
   dyn_cols = dims(1)

   ! Set ICs.  Either from analytic expressions or read from file.

   if (analytic_ic_active() .and. (iam < par%nprocs)) then

      ! PHIS has already been set by set_phis.  Get local copy for
      ! possible use in setting T and PS in the analytic IC code.
      allocate(phis_tmp(npsq,nelemd))
      do ie = 1, nelemd
         k = 1
         do j = 1, np
            do i = 1, np
               phis_tmp(k,ie) = elem(ie)%state%phis(i,j)
               k = k + 1
            end do
         end do
      end do

      inic_wet = .false.
      allocate(glob_ind(npsq * nelemd))
      j = 1
      do ie = 1, nelemd
         do i = 1, npsq
            ! Create a global(ish) column index
            glob_ind(j) = elem(ie)%GlobalId
            j = j + 1
         end do
      end do

      ! First, initialize all the variables, then assign
      allocate(dbuf4(npsq, nlev, nelemd, (qsize + 4)))
      dbuf4 = 0.0_r8
      allocate(m_ind(qsize))
      do m_cnst = 1, qsize
         m_ind(m_cnst) = m_cnst
      end do

      ! Init tracers on the GLL grid.  Note that analytic_ic_set_ic makes
      ! use of cnst_init_default for the tracers except water vapor.

      call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind,  &
         PS=dbuf4(:,1,:,(qsize+1)), U=dbuf4(:,:,:,(qsize+2)),      &
         V=dbuf4(:,:,:,(qsize+3)), T=dbuf4(:,:,:,(qsize+4)),       &
         Q=dbuf4(:,:,:,1:qsize), m_cnst=m_ind, mask=pmask(:),      &
         PHIS_IN=PHIS_tmp)
      deallocate(m_ind)
      deallocate(glob_ind)
      deallocate(phis_tmp)
      do ie = 1, nelemd
         indx = 1
         do j = 1, np
            do i = 1, np
               ! PS
               elem(ie)%state%psdry(i,j) = dbuf4(indx, 1, ie, (qsize+1))
               ! U
               elem(ie)%state%v(i,j,1,:,1) = dbuf4(indx, :, ie, (qsize+2))
               ! V
               elem(ie)%state%v(i,j,2,:,1) = dbuf4(indx, :, ie, (qsize+3))
               ! T
               elem(ie)%state%T(i,j,:,1) = dbuf4(indx, :, ie, (qsize+4))
               indx = indx + 1
            end do
         end do
      end do

      ! Tracers to be advected on GLL grid.
      ! Note that fvm tracers are initialized below.
      do m_cnst = 1, qsize
         do ie = 1, nelemd
            qtmp(:,:,:,ie,m_cnst) = 0.0_r8
            indx = 1
            do j = 1, np
               do i = 1, np
                  ! Set qtmp at the unique columns only
                  if (pmask(((ie - 1) * npsq) + indx)) then
                     qtmp(i,j,:,ie,m_cnst) = dbuf4(indx, :, ie, m_cnst)
                  end if
                  indx = indx + 1
               end do
            end do
         end do
      end do
      deallocate(dbuf4)


   else

      ! Read ICs from file.  Assume all fields in the initial file are on the GLL grid.

      allocate(dbuf2(npsq,nelemd))
      allocate(dbuf3(npsq,nlev,nelemd))

      ! Check that columns in IC file match grid definition.
      call check_file_layout(fh_ini, elem, dyn_cols, 'ncdata', .true., dimname)

      ! Read 2-D field

      fieldname  = 'PS'
      fieldname2 = 'PSDRY'
      if (dyn_field_exists(fh_ini, trim(fieldname), required=.false.)) then
         inic_wet = .true.
         call read_dyn_var(trim(fieldname), fh_ini, dimname, dbuf2)
      elseif (dyn_field_exists(fh_ini, trim(fieldname2), required=.false.)) then
         inic_wet = .false.
         call read_dyn_var(trim(fieldname2), fh_ini, dimname, dbuf2)
      else
         call endrun(trim(subname)//': PS or PSDRY must be on GLL grid')
      end if
#ifndef planet_mars
      if (iam < par%nprocs) then
         if (minval(dbuf2, mask=reshape(pmask, (/npsq,nelemd/))) < 10000._r8) then
            call endrun(trim(subname)//': Problem reading ps or psdry field -- bad values')
         end if
      end if
#endif
      do ie = 1, nelemd
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%psdry(i,j) = dbuf2(indx,ie) ! can be either wet or dry ps
               indx = indx + 1
            end do
         end do
      end do

      ! Read in 3-D fields

      if (dyn_field_exists(fh_ini, 'U')) then
         call read_dyn_var('U', fh_ini, dimname, dbuf3)
      else
         call endrun(trim(subname)//': U not found')
      end if
      do ie = 1, nelemd
         elem(ie)%state%v = 0.0_r8
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%v(i,j,1,:,1) = dbuf3(indx,:,ie)
               indx = indx + 1
            end do
         end do
      end do

      if (dyn_field_exists(fh_ini, 'V')) then
         call read_dyn_var('V', fh_ini, dimname, dbuf3)
      else
         call endrun(trim(subname)//': V not found')
      end if
      do ie = 1, nelemd
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%v(i,j,2,:,1) = dbuf3(indx,:,ie)
               indx = indx + 1
            end do
         end do
      end do

      if (dyn_field_exists(fh_ini, 'T')) then
         call read_dyn_var('T', fh_ini, dimname, dbuf3)
      else
         call endrun(trim(subname)//': T not found')
      end if
      do ie=1,nelemd
         elem(ie)%state%T = 0.0_r8
         indx = 1
         do j = 1, np
            do i = 1, np
               elem(ie)%state%T(i,j,:,1) = dbuf3(indx,:,ie)
               indx = indx + 1
            end do
         end do
      end do

      if (pertlim .ne. 0.0_r8) then
         if (masterproc) then
            write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
               'by +/- ', pertlim, ' to initial temperature field'
         end if

         call random_seed(size=rndm_seed_sz)
         allocate(rndm_seed(rndm_seed_sz))

         do ie = 1, nelemd
            ! seed random number generator based on element ID
            ! (possibly include a flag to allow clock-based random seeding)
            rndm_seed = elem(ie)%GlobalId
            call random_seed(put=rndm_seed)
            do i = 1, np
               do j = 1, np
                  do k = 1, nlev
                     call random_number(pertval)
                     pertval = 2.0_r8*pertlim*(0.5_r8 - pertval)
                     elem(ie)%state%T(i,j,k,1) = elem(ie)%state%T(i,j,k,1)*(1.0_r8 + pertval)
                  end do
               end do
            end do
         end do

         deallocate(rndm_seed)
      end if

      ! Cleanup
      deallocate(dbuf2)
      deallocate(dbuf3)

   end if ! analytic_ic_active

   ! Read in or cold-initialize all the tracer fields.
   ! Data is read in on the GLL grid.
   ! Both GLL and FVM tracer fields are initialized based on the
   ! dimension qsize or ntrac for GLL or FVM tracers respectively.
   ! Data is only read in on GLL so if FVM tracers are active,
   ! interpolation is performed.
   !
   ! If analytic ICs are being used, we allow constituents in an initial
   ! file to overwrite mixing ratios set by the default constituent initialization
   ! except for the water species.

   if (ntrac > qsize) then
      if (ntrac < pcnst) then
         write(errmsg, '(a,3(i0,a))') ': ntrac (',ntrac,') > qsize (',qsize, &
            ') but < pcnst (',pcnst,')'
         call endrun(trim(subname)//errmsg)
      end if
   else if (qsize < pcnst) then
      write(errmsg, '(a,2(i0,a))') ': qsize (',qsize,') < pcnst (',pcnst,')'
      call endrun(trim(subname)//errmsg)
   end if

   ! If using analytic ICs the initial file only needs the horizonal grid
   ! dimension checked in the case that the file contains constituent mixing
   ! ratios.
   do m_cnst = 1, pcnst
      if (cnst_read_iv(m_cnst) .and. .not. cnst_is_a_water_species(cnst_name(m_cnst))) then
         if (dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)), required=.false.)) then
            call check_file_layout(fh_ini, elem, dyn_cols, 'ncdata', .true., dimname)
            exit
         end if
      end if
   end do

   allocate(dbuf3(npsq,nlev,nelemd))

   do m_cnst = 1, pcnst

      if (analytic_ic_active() .and. cnst_is_a_water_species(cnst_name(m_cnst))) cycle

      found = .false.
      if (cnst_read_iv(m_cnst)) then
         found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)), required=.false.)
      end if

      if (found) then
         call read_dyn_var(trim(cnst_name(m_cnst)), fh_ini, dimname, dbuf3)
      else
         call cnst_init_default(m_cnst, latvals, lonvals, dbuf3, pmask)
      end if

      do ie = 1, nelemd
         ! Copy tracers defined on GLL grid into Eulerian array
         ! Make sure tracers have at least minimum value
         do k=1, nlev
            indx = 1
            do j = 1, np
               do i = 1, np
                  ! Set qtmp at the unique columns only: zero non-unique columns
                  if (pmask(((ie - 1) * npsq) + indx)) then
                     qtmp(i,j, k, ie, m_cnst) = max(qmin(m_cnst),dbuf3(indx,k,ie))
                  else
                     qtmp(i,j, k, ie, m_cnst) = 0.0_r8
                  end if
                  indx = indx + 1
               end do
            end do
         end do
      end do

   end do ! pcnst

   ! Cleanup
   deallocate(dbuf3)

   ! Put the error handling back the way it was
   call pio_seterrorhandling(fh_ini, pio_errtype)

   ! Cleanup
   deallocate(pmask)
   deallocate(latvals)
   deallocate(lonvals)

   if (associated(ldof)) then
      deallocate(ldof)
      nullify(ldof)
   end if

   ! once we've read or initialized all the fields we do a boundary exchange to
   ! update the redundent columns in the dynamics
   if(iam < par%nprocs) then
      call initEdgeBuffer(par, edge, elem, (3+pcnst)*nlev + 2 )
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVpack(edge, elem(ie)%state%psdry,1,kptr,ie)
      kptr = kptr + 1
      call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
      kptr = kptr + (2 * nlev)
      call edgeVpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
      kptr = kptr + nlev
      call edgeVpack(edge, qtmp(:,:,:,ie,:),nlev*pcnst,kptr,ie)
   end do
   if(iam < par%nprocs) then
      call bndry_exchange(par,edge,location='read_inidat')
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edge, elem(ie)%state%psdry,1,kptr,ie)
      kptr = kptr + 1
      call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,ie)
      kptr = kptr + (2 * nlev)
      call edgeVunpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,ie)
      kptr = kptr + nlev
      call edgeVunpack(edge, qtmp(:,:,:,ie,:),nlev*pcnst,kptr,ie)
   end do

   if (inic_wet) then
      !
      ! convert to dry
      !
      ! (this has to be done after edge-exchange since shared points between elements are only
      !  initialized in one element and not the other!)
      !
      if (par%masterproc) then
         write(iulog,*) 'Convert specific/wet mixing ratios to dry'
      end if

      allocate(factor_array(np,np,nlev,nelemd))
      !
      ! compute: factor_array = 1/(1-sum(q))
      !
      factor_array(:,:,:,:) = 1.0_r8
      do ie = 1, nelemd
         do k = dry_air_species_num+1, thermodynamic_active_species_num
            m_cnst = thermodynamic_active_species_idx(k)
            factor_array(:,:,:,ie) = factor_array(:,:,:,ie) - qtmp(:,:,:,ie,m_cnst)
         end do
      end do
      factor_array(:,:,:,:) = 1.0_r8/factor_array(:,:,:,:)

      do m_cnst = 1, pcnst
         if (cnst_type(m_cnst) == 'wet') then
            do ie = 1, nelemd
               do k = 1, nlev
                  do j = 1, np
                     do i = 1, np

                        ! convert wet mixing ratio to dry
                        qtmp(i,j,k,ie,m_cnst) = qtmp(i,j,k,ie,m_cnst) * factor_array(i,j,k,ie)

                        ! truncate negative values if they were not analytically specified
                        if (.not. analytic_ic_active()) then
                           qtmp(i,j,k,ie,m_cnst) = max(qmin(m_cnst), qtmp(i,j,k,ie,m_cnst))
                        end if
                     end do
                  end do
               end do
            end do
         end if
      end do

      ! initialize dp3d and qdp
      !
      ! compute: factor_array = 1/(1+sum(q))

      factor_array(:,:,:,:) = 1.0_r8
      do ie = 1, nelemd
         do k = dry_air_species_num+1, thermodynamic_active_species_num
            m_cnst = thermodynamic_active_species_idx(k)
            factor_array(:,:,:,ie) = factor_array(:,:,:,ie) + qtmp(:,:,:,ie,m_cnst)
         end do
      end do
      factor_array(:,:,:,:) = 1.0_r8/factor_array(:,:,:,:)
      do ie = 1, nelemd
         ! pstmp is the wet ps
         pstmp = elem(ie)%state%psdry(:,:)
         ! start accumulating the dry air pressure differences across each layer
         elem(ie)%state%psdry(:,:) = hyai(1)*ps0
         do k=1,nlev
            do j = 1,np
               do i = 1,np
                  dp_tmp = ((hyai(k+1) - hyai(k))*ps0) + &
                           ((hybi(k+1) - hybi(k))*pstmp(i,j))
                  if (.not. analytic_ic_active()) then

                     ! if analytic_ic then the surface pressure is already dry
                     ! (note that it is not correct to convert to moist pressure
                     ! in analytic_ic and not have the #ifndef statement here
                     ! since the dry levels are in a different location than
                     ! what is obtained from algorithm below)

                     ! convert dp_tmp to dry
                     dp_tmp = dp_tmp*factor_array(i,j,k,ie)
                  end if

                  elem(ie)%state%dp3d(i,j,k,:) = dp_tmp

                  ! compute dry surface pressure; note that at this point
                  !
                  ! dp3d .NE. (hyai(k+1) - hyai(k))*ps0 + (hybi(k+1) - hybi(k))*ps(i,j)

                  elem(ie)%state%psdry(i,j) = elem(ie)%state%psdry(i,j)+elem(ie)%state%dp3d(i,j,k,1)
               end do
            end do
         end do
      end do

      deallocate(factor_array)

   else

      ! initial condition is based on dry surface pressure and constituents
      !
      ! we only need to initialize state%dp3d

      do ie = 1, nelemd
         do k = 1, nlev
            do j = 1, np
               do i = 1, np
                  elem(ie)%state%dp3d(i,j,k,:) = (hyai(k+1) - hyai(k))*ps0 + &
                                                 (hybi(k+1) - hybi(k))*elem(ie)%state%psdry(i,j)
               end do
            end do
         end do
      end do
   end if

   ! scale PS to achieve prescribed dry mass following FV dycore (dryairm.F90)
#ifndef planet_mars
   if (runtype == 0) then
      initial_global_ave_dry_ps = 98288.0_r8
      if (.not. associated(fh_topo)) then
         initial_global_ave_dry_ps = 101325._r8 - 245._r8
      end if
      if (simple_phys) then
         initial_global_ave_dry_ps = 0                  !do not scale psdry
      end if
      if (iam < par%nprocs) then
        call prim_set_dry_mass(elem, hvcoord, initial_global_ave_dry_ps, qtmp)
      end if
   endif
#endif
   ! store Q values:
   !
   ! if CSLAM is NOT active then state%Qdp for all constituents
   ! if CSLAM active then we only advect water vapor and condensate
   ! loading tracers in state%qdp

   if (ntrac > 0) then
      do ie = 1, nelemd
         do nq = 1, thermodynamic_active_species_num
            m_cnst = thermodynamic_active_species_idx(nq)
            do k = 1, nlev
               do j = 1, np
                  do i = 1, np
                     elem(ie)%state%Qdp(i,j,k,nq,:) = &
                                 elem(ie)%state%dp3d(i,j,k,1)*qtmp(i,j,k,ie,m_cnst)
                  end do
               end do
            end do
         end do
      end do
   else
      do ie = 1, nelemd
         do m_cnst = 1, qsize
            do k = 1, nlev
               do j = 1, np
                  do i = 1, np
                     elem(ie)%state%Qdp(i,j,k,m_cnst,:)=&
                        elem(ie)%state%dp3d(i,j,k,1)*qtmp(i,j,k,ie,m_cnst)
                  end do
               end do
            end do
         end do
      end do
   end if

   ! interpolate fvm tracers and fvm pressure variables

   if (ntrac > 0) then
      if (par%masterproc) then
         write(iulog,*) 'Initializing dp_fvm from spectral element dp'
      end if

      do ie = 1, nelemd

         ! note that the area over fvm cells as computed from subcell_integration is up to 1.0E-6
         ! different than the areas (exact) computed by CSLAM
         !
         ! Map the constituents which are also to be transported by dycore
         call dyn2fvm_mass_vars(elem(ie)%state%dp3d(:,:,:,1),elem(ie)%state%psdry(:,:),&
            qtmp(:,:,:,ie,1:ntrac),&
            dyn_in%fvm(ie)%dp_fvm(1:nc,1:nc,:),dyn_in%fvm(ie)%psC(1:nc,1:nc),&
            dyn_in%fvm(ie)%c(1:nc,1:nc,:,1:ntrac),&
            ntrac,elem(ie)%metdet,dyn_in%fvm(ie)%inv_se_area_sphere(1:nc,1:nc))
      end do

      if(par%masterproc) then
         write(iulog,*) 'FVM tracers, FVM pressure variables and se_area_sphere initialized.'
      end if

   end if    ! (ntrac > 0)

   ! Cleanup
   deallocate(qtmp)

   do ie = 1, nelemd
      do t = 2, timelevels
         elem(ie)%state%v(:,:,:,:,t) = elem(ie)%state%v(:,:,:,:,1)
         elem(ie)%state%T(:,:,:,t)   = elem(ie)%state%T(:,:,:,1)
      end do
   end do

   if(iam < par%nprocs) then
      call FreeEdgeBuffer(edge)
   end if

end subroutine read_inidat


!========================================================================================

subroutine set_phis(dyn_in)

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

   type(element_t), pointer         :: elem(:)

   real(r8), allocatable            :: phis_tmp(:,:)      ! (npsp,nelemd)
   real(r8), allocatable            :: phis_phys_tmp(:,:) ! (fv_nphys**2,nelemd)

   integer                          :: i, ie, indx, j, kptr
   integer                          :: ierr, pio_errtype

   character(len=max_fieldname_len) :: fieldname
   character(len=max_hcoordname_len):: grid_name
   integer                          :: dims(2)
   integer                          :: dyn_cols
   integer                          :: ncol_did
   integer                          :: ncol_size

   integer(iMap), pointer           :: ldof(:)            ! Basic (2D) grid dof
   logical,  allocatable            :: pmask(:)           ! (npsq*nelemd) unique columns

   ! Variables for analytic initial conditions
   integer,  allocatable            :: glob_ind(:)
   logical,  allocatable            :: pmask_phys(:)
   real(r8), pointer                :: latvals_deg(:)
   real(r8), pointer                :: lonvals_deg(:)
   real(r8), allocatable            :: latvals(:)
   real(r8), allocatable            :: lonvals(:)
   real(r8), allocatable            :: latvals_phys(:)
   real(r8), allocatable            :: lonvals_phys(:)

   character(len=*), parameter      :: subname='set_phis'
   !----------------------------------------------------------------------------

   fh_topo => topo_file_get_id()

   if (iam < par%nprocs) then
      elem => dyn_in%elem
   else
      nullify(elem)
   end if

   allocate(phis_tmp(npsq,nelemd))
   phis_tmp = 0.0_r8

   if (fv_nphys > 0) then
      allocate(phis_phys_tmp(fv_nphys**2,nelemd))
      phis_phys_tmp = 0.0_r8
      do ie=1,nelemd
        elem(ie)%sub_elem_mass_flux=0.0_r8
#ifdef waccm_debug
        dyn_in%fvm(ie)%CSLAM_gamma = 0.0_r8
#endif
      end do
   end if

   ! Set mask to indicate which columns are active in GLL grid.
   nullify(ldof)
   call cam_grid_get_gcid(cam_grid_id('GLL'), ldof)
   allocate(pmask(npsq*nelemd))
   pmask(:) = (ldof /= 0)
   deallocate(ldof)

   if (associated(fh_topo)) then

      ! Set PIO to return error flags.
      call pio_seterrorhandling(fh_topo, PIO_BCAST_ERROR, pio_errtype)

      ! Set name of grid object which will be used to read data from file
      ! into internal data structure via PIO.
      if (fv_nphys == 0) then
         grid_name = 'GLL'
      else
         grid_name = 'physgrid_d'
      end if

      ! Get number of global columns from the grid object and check that
      ! it matches the file data.
      call cam_grid_dimensions(grid_name, dims)
      dyn_cols = dims(1)

      ! The dimension of the unstructured grid in the TOPO file is 'ncol'.
      ierr = pio_inq_dimid(fh_topo, 'ncol', ncol_did)
      if (ierr /= PIO_NOERR) then
         call endrun(subname//': dimension ncol not found in bnd_topo file')
      end if
      ierr = pio_inq_dimlen(fh_topo, ncol_did, ncol_size)
      if (ncol_size /= dyn_cols) then
         if (masterproc) then
            write(iulog,*) subname//': ncol_size=', ncol_size, ' : dyn_cols=', dyn_cols
         end if
         call endrun(subname//': ncol size in bnd_topo file does not match grid definition')
      end if

      fieldname = 'PHIS'
      if (dyn_field_exists(fh_topo, trim(fieldname))) then
         if (fv_nphys == 0) then
           call read_dyn_var(fieldname, fh_topo, 'ncol', phis_tmp)
         else
           call read_phys_field_2d(fieldname, fh_topo, 'ncol', phis_phys_tmp)
           call map_phis_from_physgrid_to_gll(dyn_in%fvm, elem, phis_phys_tmp, &
                phis_tmp, pmask)
         end if
      else
         call endrun(subname//': Could not find PHIS field on input datafile')
      end if

      ! Put the error handling back the way it was
      call pio_seterrorhandling(fh_topo, pio_errtype)

   else if (analytic_ic_active() .and. (iam < par%nprocs)) then

      ! lat/lon needed in radians
      latvals_deg => cam_grid_get_latvals(cam_grid_id('GLL'))
      lonvals_deg => cam_grid_get_lonvals(cam_grid_id('GLL'))
      allocate(latvals(np*np*nelemd))
      allocate(lonvals(np*np*nelemd))
      latvals(:) = latvals_deg(:)*deg2rad
      lonvals(:) = lonvals_deg(:)*deg2rad

      allocate(glob_ind(npsq*nelemd))
      j = 1
      do ie = 1, nelemd
         do i = 1, npsq
            ! Create a global(ish) column index
            glob_ind(j) = elem(ie)%GlobalId
            j = j + 1
         end do
      end do
      call analytic_ic_set_ic(vcoord, latvals, lonvals, glob_ind, &
                              PHIS_OUT=phis_tmp, mask=pmask(:))
      deallocate(glob_ind)

      if (fv_nphys > 0) then

         ! initialize PHIS on physgrid
         allocate(latvals_phys(fv_nphys*fv_nphys*nelemd))
         allocate(lonvals_phys(fv_nphys*fv_nphys*nelemd))
         indx = 1
         do ie = 1, nelemd
            do j = 1, fv_nphys
               do i = 1, fv_nphys
                  latvals_phys(indx) = dyn_in%fvm(ie)%center_cart_physgrid(i,j)%lat
                  lonvals_phys(indx) = dyn_in%fvm(ie)%center_cart_physgrid(i,j)%lon
                  indx = indx + 1
               end do
            end do
         end do

         allocate(pmask_phys(fv_nphys*fv_nphys*nelemd))
         pmask_phys(:) = .true.
         allocate(glob_ind(fv_nphys*fv_nphys*nelemd))

         j = 1
         do ie = 1, nelemd
            do i = 1, fv_nphys*fv_nphys
               ! Create a global(ish) column index
               glob_ind(j) = elem(ie)%GlobalId
               j = j + 1
            end do
         end do

         call analytic_ic_set_ic(vcoord, latvals_phys, lonvals_phys, glob_ind, &
                                 PHIS_OUT=phis_phys_tmp, mask=pmask_phys)

         deallocate(latvals_phys)
         deallocate(lonvals_phys)
         deallocate(pmask_phys)
         deallocate(glob_ind)
      end if

   end if

   deallocate(pmask)

   ! Set PHIS in element objects
   do ie = 1, nelemd
      elem(ie)%state%phis = 0.0_r8
      indx = 1
      do j = 1, np
         do i = 1, np
            elem(ie)%state%phis(i,j) = phis_tmp(indx, ie)
            indx = indx + 1
         end do
      end do
   end do
   if (fv_nphys > 0) then
      do ie = 1, nelemd
         dyn_in%fvm(ie)%phis_physgrid = RESHAPE(phis_phys_tmp(:,ie),(/fv_nphys,fv_nphys/))
      end do
   end if

   deallocate(phis_tmp)
   if (fv_nphys > 0) then
      deallocate(phis_phys_tmp)
   end if

   ! boundary exchange to update the redundent columns in the element objects
   do ie = 1, nelemd
      kptr = 0
      call edgeVpack(edgebuf, elem(ie)%state%phis, 1, kptr, ie)
   end do
   if(iam < par%nprocs) then
      call bndry_exchange(par, edgebuf, location=subname)
   end if
   do ie = 1, nelemd
      kptr = 0
      call edgeVunpack(edgebuf, elem(ie)%state%phis,1,kptr,ie)
   end do

end subroutine set_phis

!========================================================================================

subroutine check_file_layout(file, elem, dyn_cols, file_desc, dyn_ok, dimname)

   ! This routine is only called when data will be read from the initial file.  It is not
   ! called when the initial file is only supplying vertical coordinate info.

   type(file_desc_t), pointer       :: file
   type(element_t),   pointer       :: elem(:)
   integer,           intent(in)    :: dyn_cols
   character(len=*),  intent(in)    :: file_desc
   logical,           intent(in)    :: dyn_ok ! .true. iff ncol_d is okay
   character(len=*),  intent(out)   :: dimname

   integer                          :: ncol_did, ncol_size
   integer                          :: ierr
   integer                          :: ie, i, j
   integer                          :: grid_id
   integer                          :: indx
   real(r8)                         :: dbuf2(npsq, nelemd)
   logical                          :: found
   character(len=max_fieldname_len) :: dimname2, coordname

   character(len=*), parameter      :: subname = 'check_file_layout'
   !----------------------------------------------------------------------------

   ! Check that number of columns in IC file matches grid definition.

   call cam_grid_get_dim_names(cam_grid_id(ini_grid_name), dimname, dimname2)

   ierr = pio_inq_dimid(file, trim(dimname), ncol_did)
   if (ierr /= PIO_NOERR) then
      call endrun(subname//': ERROR: either ncol or ncol_d dimension not found in ' &
         //trim(file_desc)//' file')
   end if

   ierr = pio_inq_dimlen(file, ncol_did, ncol_size)
   if (ncol_size /= dyn_cols) then
      if (masterproc) then
         write(iulog, '(a,2(a,i0))') trim(subname), ': ncol_size=', ncol_size, &
             ' : dyn_cols=', dyn_cols
      end if
      call endrun(subname//': ERROR: dimension ncol size not same as in ncdata file')
   end if

   ! Set coordinate name associated with dimname.
   if (dimname == 'ncol') then
      coordname = 'lat'
   else
      coordname = 'lat_d'
   end if

   !! Check to make sure file is in correct order
   call read_dyn_var(coordname, file, dimname, dbuf2)
   found = .true.
   do ie = 1, nelemd
      indx = 1
      do j = 1, np
         do i = 1, np
            if (abs(dbuf2(indx,ie)) > 1.e-12_r8) then
               if (abs((elem(ie)%spherep(i,j)%lat*rad2deg - dbuf2(indx,ie)) / &
                    dbuf2(indx,ie)) > 1.0e-10_r8) then
                  write(iulog, '(2a,4(i0,a),f11.5,a,f11.5)')                  &
                       "ncdata file latitudes not in correct column order",   &
                       ' on task ', iam, ': elem(', ie, ')%spherep(', i,      &
                       ', ', j, ')%lat = ', elem(ie)%spherep(i,j)%lat,        &
                       ' /= ', dbuf2(indx, ie)*deg2rad
                  call shr_sys_flush(iulog)
                  found = .false.
               end if
            end if
            indx = indx + 1
         end do
      end do
   end do
   if (.not. found) then
      call endrun("ncdata file latitudes not in correct column order")
   end if

   if (dimname == 'ncol') then
      coordname = 'lon'
   else
      coordname = 'lon_d'
   end if

   call read_dyn_var(coordname, file, dimname, dbuf2)
   do ie = 1, nelemd
      indx = 1
      do j = 1, np
         do i = 1, np
            if (abs(dbuf2(indx,ie)) > 1.e-12_r8) then
               if (abs((elem(ie)%spherep(i,j)%lon*rad2deg - dbuf2(indx,ie)) / &
                    dbuf2(indx,ie)) > 1.0e-10_r8) then
                  write(iulog, '(2a,4(i0,a),f11.5,a,f11.5)')                  &
                       "ncdata file longitudes not in correct column order",  &
                       ' on task ', iam, ': elem(', ie, ')%spherep(', i,      &
                       ', ', j, ')%lon = ', elem(ie)%spherep(i,j)%lon,        &
                       ' /= ', dbuf2(indx, ie)*deg2rad
                  call shr_sys_flush(iulog)
                  found = .false.
               end if
            end if
            indx = indx + 1
         end do
      end do
   end do
   if (.not. found) then
      call endrun("ncdata file longitudes not in correct column order")
   end if
end subroutine check_file_layout

!========================================================================================

logical function dyn_field_exists(fh, fieldname, required)

   use pio,            only: var_desc_t, PIO_inq_varid
   use pio,            only: PIO_NOERR

   type(file_desc_t), intent(in) :: fh
   character(len=*),  intent(in) :: fieldname
   logical, optional, intent(in) :: required

   ! Local variables
   logical                  :: found
   logical                  :: field_required
   integer                  :: ret
   type(var_desc_t)         :: varid
   character(len=128)       :: errormsg
   !--------------------------------------------------------------------------

   if (present(required)) then
      field_required = required
   else
      field_required = .true.
   end if

   ret = PIO_inq_varid(fh, trim(fieldname), varid)
   found = (ret == PIO_NOERR)
   if (.not. found) then
      if (field_required) then
         write(errormsg, *) trim(fieldname),' was not present in the input file.'
         call endrun('DYN_FIELD_EXISTS: '//errormsg)
      end if
   end if

   dyn_field_exists = found

end function dyn_field_exists

!========================================================================================

subroutine read_dyn_field_2d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:, :)

   ! Local variables
   logical                  :: found
   real(r8)                 :: fillvalue
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), fh, dimname, 1, npsq, 1, nelemd, buffer,    &
        found, gridname=ini_grid_name, fillvalue=fillvalue)
   if(.not. found) then
      call endrun('READ_DYN_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   ! Set NaNs or fillvalue points to zero
   where (isnan(buffer) .or. (buffer==fillvalue)) buffer = 0.0_r8

end subroutine read_dyn_field_2d

!========================================================================================

subroutine read_dyn_field_3d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:,:,:)

   ! Local variables
   logical                  :: found
   real(r8)                 :: fillvalue
   !----------------------------------------------------------------------------

   buffer = 0.0_r8
   call infld(trim(fieldname), fh, dimname, 'lev',  1, npsq, 1, nlev,         &
        1, nelemd, buffer, found, gridname=ini_grid_name, fillvalue=fillvalue)
   if(.not. found) then
      call endrun('READ_DYN_FIELD_3D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

   ! This code allows use of compiler option to set uninitialized values
   ! to NaN.  In that case infld can return NaNs where the element GLL points
   ! are not "unique columns"
   ! Set NaNs or fillvalue points to zero
   where (isnan(buffer) .or. (buffer == fillvalue)) buffer = 0.0_r8

end subroutine read_dyn_field_3d

!========================================================================================

subroutine read_phys_field_2d(fieldname, fh, dimname, buffer)

   ! Dummy arguments
   character(len=*),  intent(in)    :: fieldname
   type(file_desc_t), intent(inout) :: fh
   character(len=*),  intent(in)    :: dimname
   real(r8),          intent(inout) :: buffer(:, :)

   ! Local variables
   logical                  :: found
   !----------------------------------------------------------------------------

   call infld(trim(fieldname), fh, dimname, 1, fv_nphys**2, 1, nelemd, buffer,    &
      found, gridname='physgrid_d')
   if(.not. found) then
      call endrun('READ_PHYS_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
   end if

end subroutine read_phys_field_2d

!========================================================================================

subroutine map_phis_from_physgrid_to_gll(fvm,elem,phis_phys_tmp,phis_tmp,pmask)

   use hybrid_mod,         only: get_loop_ranges, config_thread_region
   use dimensions_mod,     only: nhc_phys
   use fvm_mapping,        only: phys2dyn
   use thread_mod,         only: horz_num_threads

   type(element_t),   intent(inout) :: elem(:)
   type (fvm_struct), intent(in)    :: fvm(:)
   real(r8)         , intent(in)    :: phis_phys_tmp(fv_nphys**2,nelemd) !physgrid phis
   real(r8)         , intent(inout) :: phis_tmp(npsq,nelemd)             !gll phis
   logical          , intent(in)    :: pmask(npsq*nelemd)

   type(hybrid_t)                   :: hybrid
   integer                          :: nets, nete, ie,i,j,indx
   real(r8),            allocatable :: fld_phys(:,:,:,:,:),fld_gll(:,:,:,:,:)
   logical                          :: llimiter(1)
   !----------------------------------------------------------------------------

   !!$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,ie)
   !hybrid = config_thread_region(par,'horizontal')
   hybrid = config_thread_region(par,'serial')

   call get_loop_ranges(hybrid, ibeg=nets, iend=nete)

   allocate(fld_phys(1-nhc_phys:fv_nphys+nhc_phys,1-nhc_phys:fv_nphys+nhc_phys,1,1,nets:nete))
   allocate(fld_gll(np,np,1,1,nets:nete))
   fld_phys = 0.0_r8
   do ie = nets, nete
      fld_phys(1:fv_nphys,1:fv_nphys,1,1,ie) = RESHAPE(phis_phys_tmp(:,ie),(/fv_nphys,fv_nphys/))
   end do
   llimiter = .true.
   call phys2dyn(hybrid,elem,fld_phys,fld_gll,nets,nete,1,1,fvm,llimiter,halo_filled=.false.)
   do ie = nets,nete
      indx = 1
      do j = 1, np
         do i = 1, np
            if (pmask(((ie - 1) * npsq) + indx)) then
               phis_tmp(indx,ie) = fld_gll(i,j,1,1,ie)
            else
               phis_tmp(indx,ie) = 0.0_r8
            end if
            indx = indx + 1
         end do
      end do
   end do
   deallocate(fld_phys)
   deallocate(fld_gll)
   !!$OMP END PARALLEL
end subroutine map_phis_from_physgrid_to_gll

!========================================================================================

subroutine write_dyn_vars(dyn_out)

   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

   character(len=fieldname_len) :: tfname
   integer                      :: ie, m
   !----------------------------------------------------------------------------

   if (ntrac > 0) then
      do ie = 1, nelemd
         call outfld('dp_fvm', RESHAPE(dyn_out%fvm(ie)%dp_fvm(1:nc,1:nc,:), &
                                       (/nc*nc,nlev/)), nc*nc, ie)
         call outfld('PSDRY_fvm', RESHAPE(dyn_out%fvm(ie)%psc(1:nc,1:nc), &
                                          (/nc*nc/)), nc*nc, ie)
         do m = 1, ntrac
            tfname = trim(cnst_name(m))//'_fvm'
            call outfld(tfname, RESHAPE(dyn_out%fvm(ie)%c(1:nc,1:nc,:,m), &
                                        (/nc*nc,nlev/)), nc*nc, ie)

            tfname = 'F'//trim(cnst_name(m))//'_fvm'
            call outfld(tfname, RESHAPE(dyn_out%fvm(ie)%fc(1:nc,1:nc,:,m),&
                                        (/nc*nc,nlev/)), nc*nc, ie)
         end do
      end do
   end if

end subroutine write_dyn_vars

!=========================================================================================
end module dyn_comp
