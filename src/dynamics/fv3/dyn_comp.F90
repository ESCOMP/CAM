module dyn_comp
! CAM interfaces to the GFDL FV3 Dynamical Core

!-----------------------------------------------------------------------
! Five prognostic state variables for the fv3 dynamics
!-----------------------------------------------------------------------
! dyn_state:
! D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!     o--------u(i,j+1)----------o
!     |           |              |
!     |           |              |
!  v(i,j)------scalar(i,j)----v(i+1,j)
!     |           |              |
!     |           |              |
!     o--------u(i,j)------------o
!
! The C grid component is "diagnostic" in that it is predicted every time step
! from the D grid variables.
!----------------------------------------------------------------------
! hydrostatic state:
!----------------------------------------------------------------------
!      u     ! D grid zonal wind (m/s)
!      v     ! D grid meridional wind (m/s)
!      p     ! temperature (K)
!      delp  ! pressure thickness (pascal)
!      q     ! specific humidity and prognostic constituents
!      qdiag ! diagnostic tracers
!----------------------------------------------------------------------
! additional non-hydrostatic state:
!----------------------------------------------------------------------
!      w     ! cell center vertical wind (m/s)
!      delz  ! layer thickness (meters)
!      ze0   ! height at layer edges for remapping
!      q_con ! total condensates
!      
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    use fv_mp_mod,       only: mp_barrier
    use physics_types,   only: physics_state, physics_tend
    use dyn_grid,        only: mytile
    use constituents,    only: pcnst, cnst_name, cnst_longname, tottnam, cnst_get_ind

    use fv_arrays_mod,   only: fv_atmos_type, fv_grid_bounds_type
    use shr_kind_mod,    only: r8 => shr_kind_r8, r4 => shr_kind_r4, i8 => shr_kind_i8

    use shr_sys_mod,     only: shr_sys_flush
    use pio,             only: file_desc_t, pio_seterrorhandling, PIO_BCAST_ERROR, &
                               pio_inq_dimid, pio_inq_varid, pio_inq_dimlen, PIO_NOERR
    use constants_mod,   only: cp_air, kappa, rvgas, rdgas
    use fv_nesting_mod,  only: twoway_nesting
    use physconst,       only: gravit, cpair, rearth,omega
    use spmd_utils,      only: masterproc, masterprocid, mpicom, npes,iam
    use spmd_utils,      only: mpi_real8, mpi_integer, mpi_character, mpi_logical
    use infnan,          only: isnan
    use cam_logfile,     only: iulog
    use ppgrid,          only: pcols, pver, begchunk, endchunk

    use cam_abortutils,  only: endrun
    use mpp_domains_mod, only: mpp_update_domains, domain2D, DGRID_NE
    use mpp_mod,                only: mpp_set_current_pelist,mpp_pe,mpp_chksum,mpp_sum,mpp_npes
    use field_manager_mod,      only: MODEL_ATMOS
    use tracer_manager_mod,     only: get_tracer_index,tracer_manager_init
    use time_manager_mod, only: &
        date_to_string, increment_date, time_type, set_time, &
        set_date, set_calendar_type, &
        days_in_month,  operator(-),  operator(<=), &
        operator(+), operator (<), operator (>),   &
        operator (/=), operator (/), operator (==),&
        operator (*), THIRTY_DAY_MONTHS, JULIAN,   &
        NOLEAP, NO_CALENDAR, INVALID_CALENDAR
   use fms_io_mod,         only: set_domain, nullify_domain
   use dimensions_mod,     only: npx, npy, npz, ncnst, pnats, dnats,nq
   use dimensions_mod,     only: qsize_condensate_loading_idx,qsize_condensate_loading_cp,qsize_condensate_loading_cv, &
                                 qsize_condensate_loading_idx_gll, qsize_condensate_loading, &
                                 cnst_name_ffsl, cnst_longname_ffsl,qsize,fv3_lcp_moist,fv3_lcv_moist,qsize_tracer_idx_cam2dyn,fv3_scale_ttend
   use fv_grid_utils_mod, only: cubed_to_latlon, great_circle_dist, mid_pt_sphere,    &
                                   ptop_min, inner_prod, get_latlon_vector, get_unit_vect2, &
                                   g_sum, latlon2xyz, cart_to_latlon, make_eta_level, f_p, project_sphere_v

   use dyn_grid,           only:  grid_ew,grid_ns


    implicit none
    private
    save

    public ::                 &
         dyn_init,            &
         dyn_run,             &
         dyn_final,           &
         dyn_readnl,          &
         dyn_register,        &
         dyn_import_t,        &
         dyn_export_t

    public calc_tot_energy_dynamics
    public :: frontgf_idx, frontga_idx, uzm_idx

type dyn_import_t
  type (fv_atmos_type),  pointer :: Atm(:) => null()
end type dyn_import_t

type dyn_export_t
  type (fv_atmos_type),  pointer :: Atm(:) => null()
end type dyn_export_t

 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real(r8),                 pointer, dimension(:,:) :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real(r8),                 pointer, dimension(:,:) :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(r8), pointer, dimension(:,:) :: lon      => null() ! local longitude axis grid box centers in radians.
     real(r8), pointer, dimension(:,:) :: lat      => null() ! local latitude axis grid box centers in radians.
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     real(r8), pointer, dimension(:) :: ak
     real(r8), pointer, dimension(:) :: bk
     real(r8), pointer, dimension(:,:) :: dx
     real(r8), pointer, dimension(:,:) :: dy
     real(r8), pointer, dimension(:,:) :: area
 end type atmos_data_type

! The FV core is always called in its "full physics" mode.  We don't want
! the dycore to know what physics package is responsible for the forcing.
logical, parameter         :: convt = .true.
logical :: first_time = .true.
! Indices for fields that are computed in the dynamics and passed to the physics
! via the physics buffer
integer, protected :: frontgf_idx  = -1
integer, protected :: frontga_idx  = -1
integer, protected :: uzm_idx = -1

! Private interfaces
interface read_dyn_var
  module procedure read_dyn_field_2d
  module procedure read_dyn_field_3d
end interface read_dyn_var

real(r8), public, allocatable :: u_dt(:,:,:), v_dt(:,:,:), t_dt(:,:,:)

!These are convenience variables for local use only, and are set to values in Atm%
integer :: dt_atmos
integer(i8) :: checksum
real(r8) ::  zvir, dt_atmos_real
logical :: first_diag = .true.

integer :: ldof_size
integer :: grid_size
real(r8), parameter :: pi = 3.1415926535897931_r8
real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8

real(r8), allocatable,dimension(:,:,:)       :: se_dyn,ke_dyn,wv_dyn,wl_dyn,wi_dyn, &
                                                wr_dyn,ws_dyn,wg_dyn,tt_dyn,mo_dyn,mr_dyn
integer, public :: isav,jsav,icolsav,lchnksav
integer, public :: iamsav = -99

!=======================================================================
contains
!=======================================================================
subroutine dyn_readnl(nlfilename)

  ! Read dynamics namelist group from atm_in and write to fv3 input.nml file
  use units,           only: getunit, freeunit
  use namelist_utils,  only: find_group_name
  use constituents,    only: pcnst
  use cam_control_mod, only: initial_run

  ! args
  character(len=*), intent(in) :: nlfilename
  
  ! Local variables
  integer                      :: unitn,unito, ierr,i,ios

  ! FV3 Namelist variables
  integer                      :: fv3_qsize_condensate_loading, fv3_npes
  
  namelist /dyn_fv3_inparm/          &
       fv3_scale_ttend, &
       fv3_lcp_moist, &
       fv3_lcv_moist, &
       fv3_qsize_condensate_loading, &
       fv3_npes

  character(len = 256), dimension(5) :: group_names = (/"main_nml","fv_core_nml","surf_map_nml","test_case_nml","fms_nml"/)
  character(len=256) :: inrec  ! first 80 characters of input record
  character(len=256) :: inrec2 ! left adjusted input record

  !--------------------------------------------------------------------------
  
  ! defaults for variables not set by build-namelist
  fv3_npes                     = npes

  ! Read the namelist (dyn_se_inparm)
  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(NLFileName), status='old' )
     call find_group_name(unitn, 'dyn_fv3_inparm', status=ierr)
     if (ierr == 0) then
        read(unitn, dyn_fv3_inparm, iostat=ierr)
        if (ierr /= 0) then
           call endrun('dyn_readnl: ERROR reading dyn_fv3_inparm namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if
  if ((fv3_lcp_moist.eqv.fv3_lcv_moist) .and. (fv3_lcv_moist.eqv..true.)) call endrun('dyn_readnl: ERROR reading dyn_fv3_inparm namelist fv3_lcp_moist and fv3_lcv_moist can not both be true')

  ! Broadcast namelist values to all PEs
  call MPI_bcast(fv3_qsize_condensate_loading, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_lcv_moist, 1, mpi_logical, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_lcp_moist, 1, mpi_logical, masterprocid, mpicom, ierr)

  if (fv3_npes <= 0) then
     call endrun('dyn_readnl: ERROR: fv3_npes must be > 0')
  end if

  qsize_condensate_loading = fv3_qsize_condensate_loading
  qsize = pcnst

  ! Create the input.nml namelist needed by the fv3dycore.
  ! Read strings one at a time from the fv3 namelist groups, strip off the leading 'fv3_' from the variable names and write to input.nml.
  ! This could be replaced by also by writing to the internal namelist file

  if (masterproc) then

     ! write the file diag_table for the fv3 diag_manager  this is just a place holder until I get netcdf output working from fv3
     unito = getunit()
     ! overwrite file if it exists.
     open( unito, file='diag_table', status='replace' )
     write(unito,'(a)')'20161003.00Z.Cxx.64bit.non-mono'
     write(unito,'(a)')'2016 10 03 00 0 0'
     close(unito)
     call freeunit(unito)
     
     write(iulog,*) 'Creating fv3 input.nml file from atm_in fv3_xxx namelist parameters'
     ! Read the namelist (main_nml)
     ! open the file input.nml
     unito = getunit()
     ! overwrite file if it exists.
     open( unito, file='input.nml', status='replace' )

     unitn = getunit()
     open( unitn, file=trim(NLFileName), status='old' )
     
     do i=1,SIZE(group_names(:))
        rewind(unitn)        
        call find_group_name(unitn, trim(group_names(i)), status=ierr)
        
        if (ierr == 0) then ! we found it now copy each line until endding '/' to input.nml

           ! write group name to input.nml
           read(unitn, '(a)', iostat=ios, end=100) inrec
           if (ios .ne. 0) call endrun('ERROR: dyn_readnl - error reading fv3 namelist')
           write(unito,'(a)') trim(inrec)
           
           ios = 0
           do while (ios <= 0)
              
              read(unitn, '(a)', iostat=ios, end=100) inrec
              
              if (ios <= 0) then  ! ios < 0  indicates an end of record condition
                 
                 ! remove leading blanks and check for leading '/'
                 inrec2 = adjustl(inrec)
                 if (inrec2(1:4) == 'fv3_') then
                    inrec2(1:4) = '    '
                 end if
                 write(unito,'(a)') trim(inrec2)
                 if (inrec2(1:1) == '/') exit
              end if
           end do
        end if
     end do
     close(unitn)
     call freeunit(unitn)
     close(unito)
     call freeunit(unito)
  end if
  return
100 continue
  call endrun('ERROR: dyn_readnl: End of file encountered while reading fv3 namelist groups')

end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_register()

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use phys_control,    only: use_gw_front, use_gw_front_igw
   use qbo,             only: qbo_use_forcing

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

   if (use_gw_front .or. use_gw_front_igw) then
      call pbuf_add_field("FRONTGF", "global", dtype_r8, (/pcols,pver/), &
         frontgf_idx)
      call pbuf_add_field("FRONTGA", "global", dtype_r8, (/pcols,pver/), &
         frontga_idx)
   end if

   if (qbo_use_forcing) then
      call pbuf_add_field("UZM", "global", dtype_r8, (/pcols,pver/), &
         uzm_idx)
   end if

end subroutine dyn_register

!=============================================================================================

subroutine dyn_init(dyn_in, dyn_out)
  
  ! DESCRIPTION: Initialize the FV dynamical core
  
  ! Initialize FV dynamical core state variables
  
  use pmgrid,          only: plev
  use infnan,          only: inf, assignment(=)
  use dyn_grid,        only: Atm, grids_on_this_pe,p_split,define_cam_grids
  use hycoef,          only: hyai, hybi, ps0
  
  use cam_history,     only: addfld, add_default, horiz_only
  use cam_pio_utils,   only: clean_iodesc_list
  use phys_control,    only: phys_getopts
  use cam_control_mod, only: initial_run
  use cam_history,     only: register_vector_field
  use fv_mp_mod,       only: fill_corners, XDir, YDir, switch_current_Atm
  use fv_diagnostics_mod, only: fv_diag, fv_diag_init, fv_diag, fv_time, prt_maxmin
  use diag_manager_mod,only: diag_manager_init,get_base_date
  use fms_mod,         only: file_exist, &
                                 error_mesg, &
                                 uppercase, FATAL, NOTE
  use time_manager_mod,only: time_type, set_time
  use time_manager,    only: get_step_size,get_curr_date
  use physconst,       only: cpwv, cpliq, cpice

  ! arguments:
   type (dyn_import_t),     intent(out) :: dyn_in
   type (dyn_export_t),     intent(out) :: dyn_out
   
   type (time_type) :: Time_atmos, Time_init, Time_end, Time_start, &
                       Time_step_atmos, Time, &
                       Time_restart, Time_step_restart,Run_length

   character(len=*), parameter :: sub='dyn_init'
   real(r8)                    :: alpha

   real(r8), pointer, dimension(:,:)            :: fC,f0
   real(r8), pointer, dimension(:,:,:)          :: grid,agrid,delp
   logical, pointer :: cubed_sphere
   type(domain2d), pointer     :: domain
   integer                     :: i,j,m
   integer :: date(6),date_init(6),date_restart(6)
   character(len=32) :: timestamp   
   integer :: calendar_type = -99

   integer :: loading

   integer :: lonid
   integer :: latid
   integer :: morec            ! latitude dimension length from dataset
   
   integer :: yy             ! CAM current year
   integer :: mm             ! CAM current month
   integer :: dd             ! CAM current day
   integer :: tod             ! CAM current time of day (sec)
   integer :: num_atmos_calls, na

   character(len = 17) :: calendar = 'NOLEAP'
   character(len = 256) :: err_msg
   integer :: restart_interval(6) = 0
   integer :: months=0, days=0, hours=0, minutes=0, seconds=0
   integer ::   sphum, liq_wat, ice_wat, rainwat, snowwat,graupel
   integer :: ixcldice, ixcldliq, ixrain, ixsnow, ixgraupel
   character(len=*), parameter :: subname = 'dyn_init'

   ! variables for initializing energy and axial angular momentum diagnostics                                                                                                             
   character (len = 3), dimension(8) :: stage = (/"dED","dAP","dBD","dAT","dAF","dAD","dAR","dBF"/)
   character (len = 70),dimension(8) :: stage_txt = (/&
        " end of previous dynamics                           ",& !dED
        " after physics increment on A-grid                  ",& !dAP
        " state after applying CAM forcing                   ",& !dBD - state after applyCAMforcing
        " state after top of atmosphere damping (Rayleigh)   ",& !dAT 
        " from previous remapping or state passed to dynamics",& !dAF - state in beginning of ksplit loop                                                               
        " before vertical remapping                          ",& !dAD - state before vertical remapping                                                                                     
        " after vertical remapping                           ",& !dAR - state at end of nsplit loop                                                                                         
        " state passed to parameterizations                  " & !dBF                                                                                                                       
        /)
   character (len = 2)  , dimension(11) :: vars = (/"WV","WL","WI","WR","WS","WG","SE","KE","MR","MO","TT"/)
   character (len = 70), dimension(11)  :: vars_descriptor = (/&
      "Total column water vapor                ",&
      "Total column cloud water                ",&
      "Total column cloud ice                  ",&
      "Total column rain                       ",&
      "Total column snow                       ",&
      "Total column graupel                    ",&
      "Total column dry static energy          ",&
      "Total column kinetic energy             ",&
      "Total column wind axial angular momentum",&
      "Total column mass axial angular momentum",&
      "Total column test tracer                "/)
   character (len = 14), dimension(11)  :: &
      vars_unit = (/&
      "kg/m2        ","kg/m2        ","kg/m2        ",                &
      "kg/m2        ","kg/m2        ","kg/m2        ","J/m2         ",&
      "J/m2         ","kg*m2/s*rad2 ","kg*m2/s*rad2 ","kg/m2        "/)

   integer :: istage, ivars
   character (len=108) :: str1, str2, str3
   integer :: is,isc,isd,ie,iec,ied,js,jsc,jsd,je,jec,jed
   integer :: fv3idx,icnst_ffsl,k
   character(len=1024) :: fieldtable(pcnst)

   real, parameter:: cv_vap = 3.*rvgas        ! < 1384.5
   real, parameter:: cv_air =  cp_air - rdgas !< = rdgas * (7/2-1) = 2.5*rdgas=717.68

   !-----------------------------------------------------------------------

   ! Setup the condensate loading arrays and fv3/cam tracer mapping and
   ! finish initializing fv3 by allocating the tracer arrays in the fv3 atm structure 

   allocate(qsize_condensate_loading_idx(qsize_condensate_loading))
   allocate(qsize_condensate_loading_idx_gll(qsize_condensate_loading))
   allocate(qsize_tracer_idx_cam2dyn(qsize))
   qsize_tracer_idx_cam2dyn(:)=-1
   allocate(qsize_condensate_loading_cp(qsize_condensate_loading))
   allocate(qsize_condensate_loading_cv(qsize_condensate_loading))
   
   allocate(cnst_name_ffsl(qsize))     ! constituent names for ffsl tracers
   allocate(cnst_longname_ffsl(qsize)) ! long name of constituents for ffsl tracers


   ! set up the condensate loading array
   if (qsize_condensate_loading > 6) then
     call endrun(subname//': fv3_qsize_condensate_loading not setup for more than 6 forms of water')
   end if
   
   ! water vapor is always tracer 1
   icnst_ffsl = 1
   cnst_name_ffsl(icnst_ffsl)='sphum'
   cnst_longname_ffsl(1) = cnst_longname(1)
   qsize_condensate_loading_idx(1) = 1
   qsize_condensate_loading_idx_gll(1) = 1
   qsize_condensate_loading_cp(1) = cpwv
   qsize_condensate_loading_cv(1) = cv_vap
   qsize_tracer_idx_cam2dyn(1) = icnst_ffsl

   call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
   if (ixcldliq > 1) then
      icnst_ffsl= icnst_ffsl+1
      cnst_name_ffsl(icnst_ffsl)='liq_wat'
      cnst_longname_ffsl(icnst_ffsl) = cnst_longname(ixcldliq)
      qsize_condensate_loading_idx(icnst_ffsl) = icnst_ffsl
      qsize_condensate_loading_idx_gll(icnst_ffsl) = ixcldliq
      qsize_condensate_loading_cp(icnst_ffsl)  = cpliq
      qsize_condensate_loading_cv(icnst_ffsl)  = cpliq
      qsize_tracer_idx_cam2dyn(ixcldliq) = icnst_ffsl
   end if

   call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
   if (ixcldice  > 1) then
      icnst_ffsl= icnst_ffsl+1
      cnst_name_ffsl(icnst_ffsl)='ice_wat'
      cnst_longname_ffsl(icnst_ffsl) = cnst_longname(ixcldice)
      qsize_condensate_loading_idx(icnst_ffsl) = icnst_ffsl
      qsize_condensate_loading_idx_gll(icnst_ffsl) = ixcldice
      qsize_condensate_loading_cp(icnst_ffsl)  = cpice
      qsize_condensate_loading_cv(icnst_ffsl)  = cpice
      qsize_tracer_idx_cam2dyn(ixcldice) = icnst_ffsl
   end if

   call cnst_get_ind('RAINQM', ixrain, abort=.false.)
   if (ixrain > 1) then
      icnst_ffsl= icnst_ffsl+1
      cnst_name_ffsl(icnst_ffsl)='rainwat'
      cnst_longname_ffsl(icnst_ffsl) = cnst_longname(ixrain)
      qsize_condensate_loading_idx(icnst_ffsl) = icnst_ffsl
      qsize_condensate_loading_idx_gll(icnst_ffsl) = ixrain
      qsize_condensate_loading_cp(icnst_ffsl)  = cpliq
      qsize_condensate_loading_cv(icnst_ffsl)  = cpliq
      qsize_tracer_idx_cam2dyn(ixrain) = icnst_ffsl
   end if

   call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)
   if (ixsnow > 1) then
      icnst_ffsl= icnst_ffsl+1
      cnst_name_ffsl(icnst_ffsl)='snowwat'
      cnst_longname_ffsl(icnst_ffsl) = cnst_longname(ixsnow)
      qsize_condensate_loading_idx(icnst_ffsl) = icnst_ffsl
      qsize_condensate_loading_idx_gll(icnst_ffsl) = ixsnow
      qsize_condensate_loading_cp(icnst_ffsl)  = cpice
      qsize_condensate_loading_cv(icnst_ffsl)  = cpice
      qsize_tracer_idx_cam2dyn(ixsnow) = icnst_ffsl
   end if

   call cnst_get_ind('GRAUQM', ixgraupel, abort=.false.)
   if (ixgraupel > 1) then
      icnst_ffsl= icnst_ffsl+1
      cnst_name_ffsl(icnst_ffsl)='graupel'
      cnst_longname_ffsl(icnst_ffsl) = cnst_longname(ixgraupel)
      qsize_condensate_loading_idx(icnst_ffsl) = icnst_ffsl
      qsize_condensate_loading_idx_gll(icnst_ffsl) = ixgraupel
      qsize_condensate_loading_cp(icnst_ffsl)  = cpice
      qsize_condensate_loading_cv(icnst_ffsl)  = cpice
      qsize_tracer_idx_cam2dyn(ixgraupel) = icnst_ffsl
   end if

   if (icnst_ffsl /= qsize_condensate_loading) &
        call endrun(subname//': ERROR: qsize_condensate_loading not equal to the number of water constituents added to q array')

   !Now add all other CAM tracer after any of the condensates in the fv3 tracer array 
   do m=1,pcnst
      if (m.ne.1.and. &
           m.ne.ixcldliq.and. &
           m.ne.ixcldice.and. &
           m.ne.ixsnow.and. &
           m.ne.ixrain.and. &
           m.ne.ixgraupel) then
         icnst_ffsl=icnst_ffsl+1
         cnst_name_ffsl(icnst_ffsl)=cnst_name(m)
         cnst_longname_ffsl(m) = cnst_longname(m)
         qsize_tracer_idx_cam2dyn(m) = icnst_ffsl
      end if
   end do

   dt_atmos = get_step_size()
   dt_atmos_real=dt_atmos

   !---------Space in ATM structure for constituents was allocated in dyn_init.
   !---------now that cam has registered all tracers create entries in fms tracer_manager
   !---------we will build fms fieldtable internal file that can be read by tracermanager
   do i=1,pcnst
      write(fieldtable(i), '(a,a,a)') '"tracer" "atmos_mod" "'//trim(cnst_name_ffsl(i))//'" /'
   end do
   
   call tracer_manager_init(fieldtable)


   do m=1,pcnst  
      !  just check condensate loading tracers as they are mapped above
      if(qsize_tracer_idx_cam2dyn(m).le.qsize_condensate_loading) then
         fv3idx  = get_tracer_index (MODEL_ATMOS, cnst_name_ffsl(qsize_tracer_idx_cam2dyn(m)) )
         if (fv3idx.ne.qsize_tracer_idx_cam2dyn(m)) then
            write(6,*)'m,fv3idx,qsize_tracer_idx_cam2dyn=',m,fv3idx,qsize_tracer_idx_cam2dyn,cnst_name_ffsl
            call endrun(subname//': ERROR: CAM/FV3 Tracer mapping incorrect')
         end if
      end if
   end do

   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je
   isc = Atm(mytile)%bd%isc
   iec = Atm(mytile)%bd%iec
   jsc = Atm(mytile)%bd%jsc
   jec = Atm(mytile)%bd%jec
   isd = Atm(mytile)%bd%isd
   ied = Atm(mytile)%bd%ied
   jsd = Atm(mytile)%bd%jsd
   jed = Atm(mytile)%bd%jed

   ! Data initialization
   dyn_in%Atm  => Atm
   dyn_out%Atm => Atm

   allocate(u_dt(isd:ied,jsd:jed,npz))
   allocate(v_dt(isd:ied,jsd:jed,npz))
   allocate(t_dt(isd:ied,jsd:jed,npz))
   u_dt(:,:,:)=0._r8; v_dt(:,:,:)=0._r8; t_dt(:,:,:)=0._r8

   fC    => atm(mytile)%gridstruct%fC
   f0    => atm(mytile)%gridstruct%f0
   grid  => atm(mytile)%gridstruct%grid_64
   agrid => atm(mytile)%gridstruct%agrid_64
   domain=> Atm(mytile)%domain
   cubed_sphere => atm(mytile)%gridstruct%cubed_sphere
   delp =>  Atm(mytile)%delp

   f0(:,:) = inf
   fC(:,:) = inf
   alpha = 0._r8

   do j=jsd,jed+1
      do i=isd,ied+1
         fC(i,j) = 2._r8*omega*( -1._r8*cos(grid(i,j,1))*cos(grid(i,j,2))*sin(alpha) + &
              sin(grid(i,j,2))*cos(alpha) )
      enddo
   enddo
   do j=jsd,jed
      do i=isd,ied
         f0(i,j) = 2._r8*omega*( -1._r8*cos(agrid(i,j,1))*cos(agrid(i,j,2))*sin(alpha) + &
              sin(agrid(i,j,2))*cos(alpha) )
      enddo
   enddo
   call mpp_update_domains( f0, domain )
   if (cubed_sphere) call fill_corners(f0, npx, npy, YDir)
   
   delp(isd:is-1,jsd:js-1,1:npz)=0._r8
   delp(isd:is-1,je+1:jed,1:npz)=0._r8
   delp(ie+1:ied,jsd:js-1,1:npz)=0._r8
   delp(ie+1:ied,je+1:jed,1:npz)=0._r8
   
   if (initial_run) then
      
      ! Read in initial data
      call read_inidat(dyn_in)
      call clean_iodesc_list()
      
   end if

 call get_curr_date(yy, mm, dd, tod)
 
 ! use namelist value (either no restart or override flag on)
 date(1) = yy
 date(2) = mm
 date(3) = dd
 date(4) = int(tod / 3600)
 date(5) = int((tod - date(4) * 3600) / 60)
 date(6) = tod - date(4) * 3600 - date(5) * 60
 
 ! read date and calendar type from restart file 
 ! override calendar type with namelist value
 
select case (uppercase(trim(calendar)))
 case ('JULIAN')
    calendar_type = JULIAN
 case ('NOLEAP')
    calendar_type = NOLEAP
 case ('THIRTY_DAY')
    calendar_type = THIRTY_DAY_MONTHS
 case ('NO_CALENDAR')
    calendar_type = NO_CALENDAR
 end select
 
 call set_calendar_type(calendar_type, err_msg)
 if (err_msg .ne. '') then
    call error_mesg('subroutine dyn_grid_init', 'ERROR in dyn_grid_init: '//trim(err_msg), FATAL)
 end if
 
 ! initialize diagnostics manager  
 timestamp = date_to_string(set_date(date(1), date(2), date(3), date(4), date(5), date(6)))
 call diag_manager_init(TIME_INIT = date)
 
 ! always override initial/base date with diag_manager value
 ! get the base date in the diag_table from the diag_manager
 ! this base date is typically the starting date for the
 ! experiment and is subtracted from the current date

 date_init = date

 ! set initial and current time types
 ! set run length and compute ending time

 Time_init  = set_date(date_init(1), date_init(2), date_init(3), date_init(4), date_init(5), date_init(6))
 Time       = set_date(date     (1), date     (2), date     (3), date     (4), date     (5), date     (6))
 Time_start = Time

 ! compute the ending time
 Time_end = Time
 do m = 1, months
    Time_end = Time_end + set_time(0, days_in_month(Time_end))
 end do
 Time_end   = Time_end + set_time(hours * 3600 + minutes * 60 + seconds, days)
 Run_length = Time_end - Time

 ! get the time that last intermediate restart file was written out.
 date_restart = date
 
 if (all(restart_interval .eq. 0)) then
    Time_restart = increment_date(Time_end, 1, 0, 0, 0, 0, 0)   ! no intermediate restart                                                                                                                                         
 else
    Time_restart = set_date(date_restart(1), date_restart(2), date_restart(3), date_restart(4), date_restart(5), date_restart(6))
    Time_restart = increment_date(Time_restart, restart_interval(1), restart_interval(2), restart_interval(3), &
         restart_interval(4), restart_interval(5), restart_interval(6))
    if (Time_restart <= Time) &
         call error_mesg('subroutine dyn_grid_init', 'The first intermediate restart time is no larger than the start time', FATAL)
 end if

 ! compute the time steps
 ! determine number of iterations through the time integration loop
 ! must be evenly divisible

 Time_step_atmos = set_time(dt_atmos, 0)
 num_atmos_calls = Run_length / Time_step_atmos

 ! initial (base) time must not be greater than current time
 if (Time_init > Time) &
      call error_mesg('subroutine dyn_grid_init', 'initial time is greater than current time', FATAL)
 
 ! make sure run length is a multiple of atmos time step
 if (num_atmos_calls * Time_step_atmos /= Run_length)  &
      call error_mesg('subroutine dyn_grid_init', 'run length must be multiple of atmosphere time step', FATAL)

 Atm(mytile)%Time_init = Time_init

 ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
 ! This data is only needed for the COARSEST grid.
 call switch_current_Atm(Atm(mytile))
 
 call set_domain ( Atm(mytile)%domain )
 fv_time = Time

!----- initialize atmos_axes and fv_dynamics diagnostics
       !I've had trouble getting this to work with multiple grids at a time; worth revisiting?
   call fv_diag_init(Atm(mytile:mytile), Atm(mytile)%atmos_axes, Time, npx, npy, npz, Atm(mytile)%flagstruct%p_ref)

   ! Forcing from physics on the FFSL grid
   call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term on FFSL grid',     gridname='FFSL')
   call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term on FFSL grid',gridname='FFSL')
   call register_vector_field('FU', 'FV')
   call addfld ('FT',  (/ 'lev' /), 'A', 'K/s', 'Temperature forcing term on FFSL grid',gridname='FFSL')

   do m = 1, qsize
     call addfld ('F'//trim(cnst_name_ffsl(m))//'_ffsl',  (/ 'lev' /), 'I', 'kg/kg/s',   &
          trim(cnst_longname(m))//' mixing ratio forcing term (q_new-q_old) on FFSL grid', gridname='FFSL')
     call addfld(tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name_ffsl(m))//' horz + vert + fixer tendency ',  &
                  gridname='FFSL')
   end do

   ! Energy diagnostics and axial angular momentum diagnostics
   call addfld ('ABS_dPSdt',  horiz_only, 'A', 'Pa/s', 'Absolute surface pressure tendency',gridname='FFSL')

   call addfld ('WV_PDC',   horiz_only, 'A', 'kg/m2','Total column water vapor lost in physics-dynamics coupling',gridname='FFSL')
   call addfld ('WL_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud water lost in physics-dynamics coupling',gridname='FFSL')
   call addfld ('WI_PDC',   horiz_only, 'A', 'kg/m2','Total column cloud ice lost in physics-dynamics coupling'  ,gridname='FFSL')
   call addfld ('WR_PDC',   horiz_only, 'A', 'kg/m2','Total column rain water lost in physics-dynamics coupling',gridname='FFSL')
   call addfld ('WS_PDC',   horiz_only, 'A', 'kg/m2','Total column snow water lost in physics-dynamics coupling',gridname='FFSL')
   call addfld ('WG_PDC',   horiz_only, 'A', 'kg/m2','Total column graupel water lost in physics-dynamics coupling'  ,gridname='FFSL')
   call addfld ('TT_PDC',   horiz_only, 'A', 'kg/m2','Total column test tracer lost in physics-dynamics coupling'  ,gridname='FFSL')

   do istage = 1,SIZE(stage)
      do ivars=1,SIZE(vars)
         write(str1,*) TRIM(ADJUSTL(vars(ivars))),TRIM(ADJUSTL("_")),TRIM(ADJUSTL(stage(istage)))
         write(str2,*) TRIM(ADJUSTL(vars_descriptor(ivars))),&
             TRIM(ADJUSTL(" ")),TRIM(ADJUSTL(stage_txt(istage)))
         write(str3,*) TRIM(ADJUSTL(vars_unit(ivars)))
         call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)),gridname='FFSL')
      end do
   end do

   allocate(se_dyn(is:ie,js:je,5))
   allocate(ke_dyn(is:ie,js:je,5))
   allocate(wv_dyn(is:ie,js:je,5))
   allocate(wl_dyn(is:ie,js:je,5))
   allocate(wi_dyn(is:ie,js:je,5))
   allocate(wr_dyn(is:ie,js:je,5))
   allocate(ws_dyn(is:ie,js:je,5))
   allocate(wg_dyn(is:ie,js:je,5))
   allocate(tt_dyn(is:ie,js:je,5))
   allocate(mr_dyn(is:ie,js:je,5))
   allocate(mo_dyn(is:ie,js:je,5))

     if (iam.eq.960) write(6,*)'ua end dyn init ',atm(mytile)%ua(is:ie,js:je,1)

end subroutine dyn_init

!=======================================================================

subroutine dyn_run(dyn_state)

  ! DESCRIPTION: Driver for the NASA finite-volume dynamical core
  use time_manager,  only: get_step_size
  use fv_control_mod,         only: ngrids
  use tracer_manager_mod,     only: get_tracer_index,&
                                    NO_TRACER

  use fv_dynamics_mod,        only: fv_dynamics
  use fv_sg_mod,              only: fv_subgrid_z
  use dyn_grid,               only: p_split,grids_on_this_pe
  use dimensions_mod,         only: npz
  use cam_history,            only: outfld, hist_fld_active
  use fv_mp_mod,              only: switch_current_Atm

  implicit none
  
  type (dyn_export_t), intent(inout) :: dyn_state


  integer :: itrac, psc,idim
  integer :: k, w_diff, nt_dyn,j,i
  real(r8), allocatable :: testarr(:,:)
  type(fv_atmos_type), pointer         :: Atm(:)  
  integer :: is,isc,isd,ie,iec,ied,js,jsc,jsd,je,jec,jed
  
  !---- Call FV dynamics -----

  Atm => dyn_state%Atm

  !-----------------------------------------------------------------------

  call mpp_set_current_pelist(Atm(mytile)%pelist, no_sync=.TRUE.)

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isc = Atm(mytile)%bd%isc
  iec = Atm(mytile)%bd%iec
  jsc = Atm(mytile)%bd%jsc
  jec = Atm(mytile)%bd%jec
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed
  
  idim=ie-is+1

  nq = Atm(mytile)%ncnst - Atm(mytile)%flagstruct%pnats
  dt_atmos = get_step_size()
  dt_atmos_real=dt_atmos

  se_dyn = 0._r8
  ke_dyn = 0._r8
  wv_dyn = 0._r8
  wl_dyn = 0._r8
  wi_dyn = 0._r8
  wr_dyn = 0._r8
  ws_dyn = 0._r8
  wg_dyn = 0._r8
  tt_dyn = 0._r8
  mo_dyn = 0._r8
  mr_dyn = 0._r8

  zvir = rvgas/rdgas - 1._r8

  Atm(mytile)%parent_grid => Atm(mytile)

  do psc=1,abs(p_split)

     call fv_dynamics(npx, npy, npz, nq, Atm(mytile)%ng, dt_atmos_real/real(abs(p_split)),&
          Atm(mytile)%flagstruct%consv_te, Atm(mytile)%flagstruct%fill,  &
          Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
          Atm(mytile)%ptop, Atm(mytile)%ks, nq,                          &
          Atm(mytile)%flagstruct%n_split, Atm(mytile)%flagstruct%q_split,&
          Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%w, Atm(mytile)%delz,           &
          Atm(mytile)%flagstruct%hydrostatic,                       &
          Atm(mytile)%pt, Atm(mytile)%delp, Atm(mytile)%q, Atm(mytile)%ps,         &
          Atm(mytile)%pe, Atm(mytile)%pk, Atm(mytile)%peln,                   &
          Atm(mytile)%pkz, Atm(mytile)%phis, Atm(mytile)%q_con,               &
          Atm(mytile)%omga, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%uc,        &
          Atm(mytile)%vc, Atm(mytile)%ak, Atm(mytile)%bk, Atm(mytile)%mfx,         &
          Atm(mytile)%mfy, Atm(mytile)%cx, Atm(mytile)%cy, Atm(mytile)%ze0,        &
          Atm(mytile)%flagstruct%hybrid_z,                          &
          Atm(mytile)%gridstruct, Atm(mytile)%flagstruct,                &
          Atm(mytile)%neststruct, Atm(mytile)%idiag, Atm(mytile)%bd,          &
          Atm(mytile)%parent_grid, Atm(mytile)%domain, &
#if ( defined CALC_ENERGY )
          Atm(mytile)%diss_est, &
          qsize,qsize_condensate_loading,qsize_condensate_loading_idx,qsize_tracer_idx_cam2dyn,qsize_condensate_loading_cp,qsize_condensate_loading_cp, &
          se_dyn, ke_dyn, wv_dyn,wl_dyn,wi_dyn,wr_dyn,ws_dyn,wg_dyn,tt_dyn,mo_dyn,mr_dyn,gravit, cpair, rearth,omega,fv3_lcp_moist,fv3_lcv_moist)
#else
          Atm(mytile)%diss_est)
#endif

     if (ngrids > 1 .and. (psc < p_split .or. p_split < 0)) then
        call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir)
     endif
     
  end do !p_split
#if ( defined CALC_ENERGY )
  call write_dyn_var(se_dyn(is:ie,js:je,1),'SE_dAF',Atm(mytile)%bd)
  call write_dyn_var(ke_dyn(is:ie,js:je,1),'KE_dAF',Atm(mytile)%bd)
  call write_dyn_var(wv_dyn(is:ie,js:je,1),'WV_dAF',Atm(mytile)%bd)
  call write_dyn_var(wl_dyn(is:ie,js:je,1),'WL_dAF',Atm(mytile)%bd)
  call write_dyn_var(wi_dyn(is:ie,js:je,1),'WI_dAF',Atm(mytile)%bd)
  call write_dyn_var(wr_dyn(is:ie,js:je,1),'WR_dAF',Atm(mytile)%bd)
  call write_dyn_var(ws_dyn(is:ie,js:je,1),'WS_dAF',Atm(mytile)%bd)
  call write_dyn_var(wg_dyn(is:ie,js:je,1),'WG_dAF',Atm(mytile)%bd)
  call write_dyn_var(tt_dyn(is:ie,js:je,1),'TT_dAF',Atm(mytile)%bd)
  call write_dyn_var(mo_dyn(is:ie,js:je,1),'MO_dAF',Atm(mytile)%bd)
  call write_dyn_var(mr_dyn(is:ie,js:je,1),'MR_dAF',Atm(mytile)%bd)
  
  call write_dyn_var(se_dyn(is:ie,js:je,2),'SE_dAD',Atm(mytile)%bd)
  call write_dyn_var(ke_dyn(is:ie,js:je,2),'KE_dAD',Atm(mytile)%bd)
  call write_dyn_var(wv_dyn(is:ie,js:je,2),'WV_dAD',Atm(mytile)%bd)
  call write_dyn_var(wl_dyn(is:ie,js:je,2),'WL_dAD',Atm(mytile)%bd)
  call write_dyn_var(wi_dyn(is:ie,js:je,2),'WI_dAD',Atm(mytile)%bd)
  call write_dyn_var(wr_dyn(is:ie,js:je,2),'WR_dAD',Atm(mytile)%bd)
  call write_dyn_var(ws_dyn(is:ie,js:je,2),'WS_dAD',Atm(mytile)%bd)
  call write_dyn_var(wg_dyn(is:ie,js:je,2),'WG_dAD',Atm(mytile)%bd)
  call write_dyn_var(tt_dyn(is:ie,js:je,2),'TT_dAD',Atm(mytile)%bd)
  call write_dyn_var(mo_dyn(is:ie,js:je,2),'MO_dAD',Atm(mytile)%bd)
  call write_dyn_var(mr_dyn(is:ie,js:je,2),'MR_dAD',Atm(mytile)%bd)

  call write_dyn_var(se_dyn(is:ie,js:je,3),'SE_dAR',Atm(mytile)%bd)
  call write_dyn_var(ke_dyn(is:ie,js:je,3),'KE_dAR',Atm(mytile)%bd)
  call write_dyn_var(wv_dyn(is:ie,js:je,3),'WV_dAR',Atm(mytile)%bd)
  call write_dyn_var(wl_dyn(is:ie,js:je,3),'WL_dAR',Atm(mytile)%bd)
  call write_dyn_var(wi_dyn(is:ie,js:je,3),'WI_dAR',Atm(mytile)%bd)
  call write_dyn_var(wr_dyn(is:ie,js:je,3),'WR_dAR',Atm(mytile)%bd)
  call write_dyn_var(ws_dyn(is:ie,js:je,3),'WS_dAR',Atm(mytile)%bd)
  call write_dyn_var(wg_dyn(is:ie,js:je,3),'WG_dAR',Atm(mytile)%bd)
  call write_dyn_var(tt_dyn(is:ie,js:je,3),'TT_dAR',Atm(mytile)%bd)
  call write_dyn_var(mo_dyn(is:ie,js:je,3),'MO_dAR',Atm(mytile)%bd)
  call write_dyn_var(mr_dyn(is:ie,js:je,3),'MR_dAR',Atm(mytile)%bd)

  call write_dyn_var(se_dyn(is:ie,js:je,4),'SE_dAT',Atm(mytile)%bd)
  call write_dyn_var(ke_dyn(is:ie,js:je,4),'KE_dAT',Atm(mytile)%bd)
  call write_dyn_var(wv_dyn(is:ie,js:je,4),'WV_dAT',Atm(mytile)%bd)
  call write_dyn_var(wl_dyn(is:ie,js:je,4),'WL_dAT',Atm(mytile)%bd)
  call write_dyn_var(wi_dyn(is:ie,js:je,4),'WI_dAT',Atm(mytile)%bd)
  call write_dyn_var(wr_dyn(is:ie,js:je,4),'WR_dAT',Atm(mytile)%bd)
  call write_dyn_var(ws_dyn(is:ie,js:je,4),'WS_dAT',Atm(mytile)%bd)
  call write_dyn_var(wg_dyn(is:ie,js:je,4),'WG_dAT',Atm(mytile)%bd)
  call write_dyn_var(tt_dyn(is:ie,js:je,4),'TT_dAT',Atm(mytile)%bd)
  call write_dyn_var(mo_dyn(is:ie,js:je,4),'MO_dAT',Atm(mytile)%bd)
  call write_dyn_var(mr_dyn(is:ie,js:je,4),'MR_dAT',Atm(mytile)%bd)  
#endif

  !-----------------------------------------------------
  !--- COMPUTE SUBGRID Z
  !-----------------------------------------------------
  !--- zero out tendencies
  u_dt(:,:,:)   = 0._r8
  v_dt(:,:,:)   = 0._r8
  t_dt(:,:,:)   = 0._r8
  
  w_diff = get_tracer_index (MODEL_ATMOS, 'w_diff' )

  ! Perform grid-scale dry adjustment if fv_sg_adj > 0
  if ( Atm(mytile)%flagstruct%fv_sg_adj > 0 ) then
     nt_dyn = nq
     if ( w_diff /= NO_TRACER ) then
        nt_dyn = nq - 1
     endif
     call fv_subgrid_z(isd, ied, jsd, jed, isc, iec, jsc, jec, npz, &
          nt_dyn, dt_atmos_real, Atm(mytile)%flagstruct%fv_sg_adj,      &
          Atm(mytile)%flagstruct%nwat, Atm(mytile)%delp, Atm(mytile)%pe,     &
          Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%pt, Atm(mytile)%q,       &
          Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%flagstruct%hydrostatic,&
          Atm(mytile)%w, Atm(mytile)%delz, u_dt, v_dt, t_dt, Atm(mytile)%flagstruct%n_sponge)
  endif
  
#ifdef USE_Q_DT
  if ( .not. Atm(mytile)%flagstruct%hydrostatic .and. w_diff /= NO_TRACER ) then
!$OMP parallel do default (none) &
!$OMP              shared (isc, iec, jsc, jec, w_diff, n, Atm, q_dt) &
!$OMP             private (k)
     do k=1, npz
        Atm(mytile)%q(isc:iec,jsc:jec,k,w_diff) = Atm(mytile)%w(isc:iec,jsc:jec,k) + w0_big
        q_dt(:,:,k,w_diff) = 0._r8
     enddo
  endif
#endif
  
#if ( defined CALC_ENERGY )
  call calc_tot_energy_dynamics(atm, 1, 1, 1, 1, 'dBF')
#endif
  
end subroutine dyn_run

!=======================================================================

subroutine dyn_final(dyn_in, dyn_out, restart_file)
  use dyn_grid,           only: grids_on_this_pe
  use fms_mod,            only: fms_end
  use fv_diagnostics_mod, only: fv_time,fv_diag
  use fv_control_mod,     only: fv_end
  implicit none

  type (dyn_import_t),      intent(inout) :: dyn_in
  type (dyn_export_t),      intent(inout) :: dyn_out
  character(len=*),optional,intent(in)    :: restart_file
  type(fv_atmos_type), pointer         :: Atm(:)  
  
  !---- Call FV dynamics -----

  Atm => dyn_in%Atm

  call fv_end(Atm, grids_on_this_pe)

  deallocate( u_dt, v_dt, t_dt)

  ! fms finalization
  call fms_end()

end subroutine dyn_final

!=============================================================================================
! Private routines
!=============================================================================================

subroutine read_inidat(dyn_in)
  use inic_analytic,         only: analytic_ic_active, analytic_ic_set_ic
  use dyn_tests_utils,       only: vc_moist_pressure,vc_dry_pressure
  use pmgrid,                only: plev
  use constituents,          only: pcnst
  use pio,                   only: file_desc_t, pio_global, pio_double, pio_offset, &
                             pio_get_att, pio_inq_dimid, pio_inq_dimlen, pio_initdecomp, pio_inq_varid, &
                             pio_read_darray, pio_setframe, file_desc_t, io_desc_t, pio_double,pio_offset_kind,&
                             pio_seterrorhandling, pio_bcast_error,pio_get_local_array_size, pio_freedecomp, var_desc_t

  use ppgrid,                only: pver
  use ncdio_atm,             only: infld
  use cam_abortutils,        only: endrun
  use constituents,          only: pcnst, cnst_name, cnst_read_iv,qmin, cnst_type
  use const_init,            only: cnst_init_default
  use cam_initfiles,         only: initial_file_get_id, topo_file_get_id, pertlim
  use cam_grid_support,      only: cam_grid_id, cam_grid_get_gcid, &
       cam_grid_dimensions, cam_grid_get_decomp, &
       cam_grid_get_latvals, cam_grid_get_lonvals,  &
       iMap
  use cam_history_support,   only: max_fieldname_len
  use hycoef,                only: hyai, hybi, ps0
  use dyn_grid,              only: mygindex,mygindex_ew,mygindex_ns, &
                                   mylindex,mylindex_ew,mylindex_ns, &
                                   mygindexdups,mygindexdups_ew,mygindexdups_ns
  use cam_pio_utils,         only: cam_pio_handle_error
  use shr_infnan_mod,  only: shr_infnan_inf_type, assignment(=), &
                             shr_infnan_posinf, shr_infnan_neginf, &
                             shr_infnan_nan, &
                             shr_infnan_isnan, shr_infnan_isinf, &
                             shr_infnan_isposinf, shr_infnan_isneginf
  use shr_sys_mod, only: shr_sys_abort
  implicit none

  type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

  real(r8), allocatable :: tmp(:,:,:)

  type(io_desc_t) :: iodesc

  logical :: found

  character(len = 40) :: fieldname,fieldname2

  integer :: lsize

  integer :: i, j, k, m, n

  type(file_desc_t), pointer :: fh_topo
  type(file_desc_t)          :: fh_ini
  type(fv_atmos_type), pointer         :: Atm(:)  


  character(len=128)               :: errmsg
  character(len=*), parameter      :: subname='READ_INIDAT'
  real(r8), allocatable            :: phis_tmp(:,:) 
  integer(iMap), pointer           :: ldof(:),ldof_ew(:),ldof_ns(:) ! Basic (2D) grid dof
  logical,  allocatable            :: pmask(:)           ! (npsq*nelemd) unique grid vals

  ! Variables for analytic initial conditions
  integer,  allocatable            :: glob_inddups(:),glob_inddups_ew(:),glob_inddups_ns(:)
  integer,  allocatable            :: m_ind(:)
  integer                          :: vcoord
  integer                          :: pio_errtype
  integer                          :: ncol_did
  integer                          :: dims(2)
  real(r8), allocatable            :: dbuf2(:,:)         ! (pcol,nblk=1)
  real(r8), allocatable            :: dbuf3(:,:,:)       ! (pcol,plev,nblk=1)
  real(r8), allocatable            :: dbuf4(:,:,:,:)       ! (pcol,plev,nblk=1,pcnst)
  real(r8), pointer                :: latvals_deg(:),latvals_rad_ew(:),latvals_rad_ns(:)
  real(r8), pointer                :: lonvals_deg(:),lonvals_rad_ew(:),lonvals_rad_ns(:)
  real(r8), allocatable            :: latvals(:)
  real(r8), allocatable            :: lonvals(:)
  real(r8), allocatable            :: latvals_rad(:)
  real(r8), allocatable            :: lonvals_rad(:)
  integer                          :: rndm_seed_sz
  integer, allocatable             :: rndm_seed(:)
  logical                          :: inic_wet !initial condition is based on wet pressure and water species
  integer                          :: m_cnst,m_cnst_ffsl
  integer                          :: indx, nq
  integer                          :: ierr,ig,err_handling
  character(len=max_fieldname_len) :: dimname, varname
  integer                          :: ncol_size
  real(r8)                         :: pertval
  real(r8), allocatable            :: pstmp(:,:), psdry(:,:)
  real(r8), pointer, dimension(:,:,:)  :: agrid,grid
  real(r8), allocatable            :: gz(:,:,:),factor_mixd2dvc(:,:,:)
  real (r8)                        :: tracermass(pcnst),delpdry
  real (r8)                        :: fv3_totwatermass, fv3_airmass
  real (r8)                        :: initial_global_ave_dry_ps,reldif

  !-----------------------------------------------------------------------
  integer                :: is,ie,js,je,isd,ied,jsd,jed
  integer                :: blksize,blksize_ew,blksize_ns
  logical                :: fv2fv3_mixratio
  integer                :: sphum, liq_wat, ice_wat, rainwat, snowwat,graupel
  real(r8)               :: u1
  real(r8), dimension(2) :: pa
  real(r8), dimension(3) :: e1,ex,ey
  integer :: fnlev,nlev_dimid
  integer :: hdim_len, ncols_ns, ncols_ew, ncols
  
  integer :: npz_dimid
  integer :: ncol_dimid
  integer :: ncol_ns_dimid
  integer :: ncol_ew_dimid
  integer :: m_ffsl
  
  type(var_desc_t) :: omegadesc
  type(var_desc_t) :: delpdesc
  type(var_desc_t) :: udesc
  type(var_desc_t) :: vdesc
  type(var_desc_t) :: usdesc
  type(var_desc_t) :: vsdesc
  type(var_desc_t) :: tdesc
  type(var_desc_t) :: psdesc
  type(var_desc_t) :: phisdesc
  type(var_desc_t), allocatable :: qdesc(:)
  type(io_desc_t),pointer :: iodesc2d, iodesc3d,iodesc3d_ns,iodesc3d_ew,iodesc3d_ns_rst,iodesc3d_ew_rst
  integer :: array_lens_3d(3), array_lens_2d(2), array_lens_1d(1)
  integer :: file_lens_2d(2), file_lens_1d(1)
  integer :: grid_id,grid_id_ns,grid_id_ew,ilen,jlen,grid_id_ns_rst,grid_id_ew_rst
  integer :: grid_dimlens(2),grid_dimlens_ns(2),grid_dimlens_ew(2),grid_dimlens_ns_rst(2),grid_dimlens_ew_rst(2)
  real(r8), allocatable :: var3d(:,:,:), var3d_ew(:,:,:), var3d_ew_tmp(:,:,:), var3d_ew_rst(:,:,:), var3d_ns(:,:,:), var3d_ns_tmp(:,:,:), var3d_ns_rst(:,:,:), var2d(:,:)

  logical  nan_check,inf_check,inf_nan_gchecks
  integer  nan_count,inf_count,ii
  real(r8) :: tmparr(10000,pver)

  Atm => dyn_in%Atm

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed

  nullify(ldof)
  nullify(ldof_ew)
  nullify(ldof_ns)

!  fh_ini  => initial_file_get_id()
  fh_topo => topo_file_get_id()
  fh_ini  = initial_file_get_id()

  Atm => dyn_in%Atm
  grid => atm(mytile)%gridstruct%grid_64
  agrid => atm(mytile)%gridstruct%agrid_64

  ! The grid name is defined in dyn_grid::define_cam_grids.
  ! Get the number of columns in the global grid.
  call cam_grid_dimensions('FFSL', dims)
  grid_size = dims(1)

  ! Set mask to indicate which columns are active
  call cam_grid_get_gcid(cam_grid_id('FFSL'), ldof)
  ldof_size=(je-js+1)*(ie-is+1)
  allocate(phis_tmp(ldof_size,1))
  allocate(pmask(ldof_size))
  phis_tmp(:,:)=0._r8
  pmask(:) = .true.

  nullify(latvals_rad_ew)
  nullify(lonvals_rad_ew)
  nullify(latvals_rad_ns)
  nullify(lonvals_rad_ns)

  blksize=(ie-is+1)*(je-js+1)
  blksize_ew=(ie-is+2)*(je-js+1)
  blksize_ns=(ie-is+1)*(je-js+2)
  allocate(latvals_rad(blksize))
  allocate(lonvals_rad(blksize))
  allocate(latvals_rad_ew(blksize_ew))
  allocate(lonvals_rad_ew(blksize_ew))
  allocate(latvals_rad_ns(blksize_ns))
  allocate(lonvals_rad_ns(blksize_ns))
  allocate(glob_inddups(blksize))
  allocate(glob_inddups_ew(blksize_ew))
  allocate(glob_inddups_ns(blksize_ns))

  do j = js, je
     do i = is, ie
        n=mylindex(i,j)
        lonvals_rad(n) = agrid(i,j,1)
        latvals_rad(n) = agrid(i,j,2)
        glob_inddups(n) = mygindexdups(i,j)
     end do
  end do

  do j = js, je
     do i = is, ie+1
        n=mylindex_ew(i,j)
        lonvals_rad_ew(n) = grid_ew(i,j,1)
        latvals_rad_ew(n) = grid_ew(i,j,2)
        glob_inddups_ew(n) = mygindexdups_ew(i,j)
     end do
  end do

  do j = js, je+1
     do i = is, ie
        n=mylindex_ns(i,j)
        lonvals_rad_ns(n) = grid_ns(i,j,1)
        latvals_rad_ns(n) = grid_ns(i,j,2)
        glob_inddups_ns(n)    = mygindexdups_ns(i,j)
     end do
  end do
  ! Set ICs.  Either from analytic expressions or read from file.

  if (analytic_ic_active()) then
     vcoord = vc_moist_pressure
     inic_wet = .true.
     ! First, initialize all the variables, then assign
     allocate(dbuf2(blksize,1))
     allocate(dbuf3(blksize,plev,1))
     allocate(dbuf4(blksize,plev, 1,pcnst))
     dbuf2 = 0.0_r8
     dbuf3 = 0.0_r8
     dbuf4 = 0.0_r8

     allocate(m_ind(pcnst))
     do m_cnst = 1, pcnst
        m_ind(m_cnst) = m_cnst
     end do

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups,PS=dbuf2)
     do j = js, je
        do i = is, ie
           ! PS
           n=mylindex(i,j)
           atm(mytile)%ps(i,j) =   dbuf2(n, 1)
        end do
     end do
     deallocate(dbuf2)
     
     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups ,            &
          PHIS_OUT=phis_tmp(:,:))

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups,            &
          T=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! T
           n=mylindex(i,j)
           atm(mytile)%pt(i,j,:) = dbuf3(n, :, 1)
        end do
     end do


     dbuf3=0._r8
     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups,            &
          U=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! U a-grid
           n=mylindex(i,j)
           atm(mytile)%ua(i,j,:) = dbuf3(n, :, 1)
        end do
     end do

     dbuf3=0._r8
     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups,            &
          V=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! V a-grid
           n=mylindex(i,j)
           atm(mytile)%va(i,j,:) = dbuf3(n, :, 1)
        end do
     end do


     deallocate(dbuf3)

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_inddups,            &
          Q=dbuf4(:,:,:,1:pcnst), m_cnst=m_ind, mask=pmask(:))
     deallocate(m_ind)

     ! Tracers to be advected on FFSL grid.
     do m_cnst = 1, pcnst
        m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
        Atm(mytile)%q(:,:,:,m_cnst_ffsl) = 0.0_r8
        do j = js, je
           do i = is, ie
              indx=mylindex(i,j)
              Atm(mytile)%q(i,j,:,m_cnst_ffsl) = dbuf4(indx, :, 1, m_cnst)
           end do
        end do
     end do
     deallocate(dbuf4)

     allocate(dbuf3(blksize_ew,plev,1))

! Get just the U wind for the ew grid points as well as ns 
! This will become atm(mytile)%v
     call analytic_ic_set_ic(vcoord, latvals_rad_ew, lonvals_rad_ew, glob_inddups_ew,            &
          U=dbuf3(:,:,:))



     do j = js, je
        do i = is, ie+1
           ! calculate rotation (u1) for D-grid V coordinate from lat/lon A-grid representation
           call mid_pt_sphere(grid(i,j,1:2),grid(i,j+1,1:2),pa)
           call get_unit_vect2(grid(i,j,1:2),grid(i,j+1,1:2),e1)
           call get_latlon_vector(pa,ex,ey)
           u1 = inner_prod(e1,ex) !v components
           ! V
           n=mylindex_ew(i,j)
           atm(mytile)%v(i,j,:) = dbuf3(n, :, 1)*u1
        end do
     end do

     deallocate(dbuf3)
     allocate(dbuf3(blksize_ns,plev,1))
     
     call analytic_ic_set_ic(vcoord, latvals_rad_ns, lonvals_rad_ns,glob_inddups_ns,            &
          U=dbuf3(:,:,:))

     do j = js, je+1
        do i = is, ie
           ! calculate rotation (u1) for D-grid U coordinate from lat/lon representation
           call mid_pt_sphere(grid(i,j,1:2),grid(i+1,j,1:2),pa)
           call get_unit_vect2(grid(i,j,1:2),grid(i+1,j,1:2),e1)
           call get_latlon_vector(pa,ex,ey)
           u1 = inner_prod(e1,ex) !u components
           !jt don't need this           u2 = inner_prod(e1,ey)
           ! U
           n=mylindex_ns(i,j)
           atm(mytile)%u(i,j,:) = dbuf3(n, :, 1)*u1
        end do
     end do

     call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,   Atm(mytile)%domain, gridtype=DGRID_NE )
     
     call cubed_to_latlon(Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%gridstruct, &
          npx, npy, npz, 1, Atm(mytile)%gridstruct%grid_type, Atm(mytile)%domain, Atm(mytile)%gridstruct%nested, Atm(mytile)%flagstruct%c2l_ord, Atm(mytile)%bd)

     deallocate(dbuf3)
     deallocate(latvals_rad)
     deallocate(lonvals_rad)
     deallocate(latvals_rad_ew)
     deallocate(latvals_rad_ns)
     deallocate(lonvals_rad_ew)
     deallocate(lonvals_rad_ns)
     deallocate(glob_inddups)
     deallocate(glob_inddups_ew)
     deallocate(glob_inddups_ns)
     !-----------------------------------------------------------------------

  else
     ! Read ICs from file. 
     call pio_seterrorhandling(fh_ini, pio_bcast_error, err_handling)

     ierr = PIO_Inq_DimID(fh_ini, 'lev', nlev_dimid)
     ierr = PIO_Inq_dimlen(fh_ini, nlev_dimid, fnlev)
     if (npz /= fnlev) then
        write(iulog,*) 'Initial condition file nlev does not match model. nlev (file, namelist):', &
             fnlev, npz
        call endrun(subname//': Initial Condition File file nlev does not match model.')
     end if
     
     ! variable descriptors of required dynamics fields

     ierr = PIO_Inq_varid(fh_ini, 'U',     udesc)
     call cam_pio_handle_error(ierr, subname//': cannot find U')
     ierr = PIO_Inq_varid(fh_ini, 'V',     Vdesc)
     call cam_pio_handle_error(ierr, subname//': cannot find V')
     ! variable descriptors of required dynamics fields
     ierr = PIO_Inq_varid(fh_ini, 'US',     usdesc)
     call cam_pio_handle_error(ierr, subname//': cannot find US')
     ierr = PIO_Inq_varid(fh_ini, 'VS',     Vsdesc)
     call cam_pio_handle_error(ierr, subname//': cannot find VS')
     ierr = PIO_Inq_varid(fh_ini, 'T',     tdesc)
     call cam_pio_handle_error(ierr, subname//': cannot find T')
     ierr = PIO_Inq_varid(fh_ini, 'PS', psdesc)
     call cam_pio_handle_error(ierr, subname//': cannot find PS')

     fieldname  = 'PS'
     fieldname2 = 'PSDRY'
     if (dyn_field_exists(fh_ini, trim(fieldname), required=.false.)) then
        inic_wet = .true.
     elseif (dyn_field_exists(fh_ini, trim(fieldname2), required=.false.)) then
        inic_wet = .false.
     else
        call endrun(trim(subname)//': PS or PSDRY must be on ncdata')
     end if

     allocate(qdesc(pcnst))

     ! check whether the restart fields on the GLL grid contain unique columns
     ! or the fv3 task structure (ncol_ns = (ie-is+1)*(je-js+2)+npes columns)
     ! or the fv3 task structure (ncol_ew = (ie-is+2)*(je-js+1)+npes columns)
     
     ierr = PIO_Inq_DimID(fh_ini, 'ncol', ncol_dimid)
     call cam_pio_handle_error(ierr, subname//': cannot find ncol')
     ierr = PIO_Inq_dimlen(fh_ini, ncol_dimid, ncols)
     
     
     ierr = PIO_Inq_DimID(fh_ini, 'ncol_ns', ncol_ns_dimid)
     call cam_pio_handle_error(ierr, subname//': cannot find ncol_ns')
     ierr = PIO_Inq_dimlen(fh_ini, ncol_ns_dimid, ncols_ns)
     
     ierr = PIO_Inq_DimID(fh_ini, 'ncol_ew', ncol_ew_dimid)
     call cam_pio_handle_error(ierr, subname//': cannot find ncol_ew')
     ierr = PIO_Inq_dimlen(fh_ini, ncol_ew_dimid, ncols_ew)
     
     grid_id = cam_grid_id('FFSL')
     grid_id_ns = cam_grid_id('FFSL_NS')
     grid_id_ew = cam_grid_id('FFSL_EW')
     grid_id_ns_rst = cam_grid_id('FFSL_NS_RST')
     grid_id_ew_rst = cam_grid_id('FFSL_EW_RST')
     call cam_grid_dimensions(grid_id, grid_dimlens)
     call cam_grid_dimensions(grid_id_ew, grid_dimlens_ew)
     call cam_grid_dimensions(grid_id_ns, grid_dimlens_ns)
     call cam_grid_dimensions(grid_id_ew_rst, grid_dimlens_ew_rst)
     call cam_grid_dimensions(grid_id_ns_rst, grid_dimlens_ns_rst)
     if (masterproc) write(iulog,*)'reading grid dimensions',grid_id,grid_id_ew,grid_id_ns,grid_id_ew_rst,grid_id_ns_rst
     if (ncols_ns /= grid_dimlens_ns(1)) then
        write(iulog,*) 'Restart file ncol_ns does not match model. ncols_ns (file, model):',&
             ncols_ns, grid_dimlens_ns(1)
        call endrun(subname//': Restart file ncols_fvm does not match model.')
     end if
     
     ilen = ie-is+1
     jlen = je-js+1
     ! create map for distributed write of 2D fields
     array_lens_2d = (/ilen,jlen/)
     file_lens_1d  = (/grid_dimlens(1)/)
     call cam_grid_get_decomp(grid_id, array_lens_2d, file_lens_1d, pio_double, iodesc2d)
     
     ! create map for distributed write of 3D fields
     array_lens_3d = (/ilen,npz, jlen/)
     file_lens_2d  = (/grid_dimlens(1), npz/)
     call cam_grid_get_decomp(grid_id, array_lens_3d, file_lens_2d, pio_double, iodesc3d)
     
     ! create map for distributed write of 3D NS fields
     array_lens_3d = (/ilen, npz, jlen+1/)
     file_lens_2d  = (/grid_dimlens_ns(1), npz/)
     call cam_grid_get_decomp(grid_id_ns, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns)
     
     ! create map for distributed write of 3D NS RST fields (reading dups - map has gindex for dups only all others 0)
     array_lens_3d = (/ilen, npz, jlen+1/)
     file_lens_2d  = (/grid_dimlens_ns_rst(1), npz/)
     call cam_grid_get_decomp(grid_id_ns_rst, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ns_rst)
     
     ! create map for distributed write of 3D EW fields
     array_lens_3d = (/ilen+1, npz, jlen/)
     file_lens_2d  = (/grid_dimlens_ew(1), npz/)
     call cam_grid_get_decomp(grid_id_ew, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew)
    
     ! create map for distributed write of 3D EW RST fields  (reading dups - map has gindex for dups only all others 0)
     array_lens_3d = (/ilen+1, npz, jlen/)
     file_lens_2d  = (/grid_dimlens_ew_rst(1), npz/)
     call cam_grid_get_decomp(grid_id_ew_rst, array_lens_3d, file_lens_2d, pio_double, iodesc3d_ew_rst)
     
     allocate(var2d(is:ie,js:je))
     var2d = 0._r8
     ! PS
     call PIO_Read_Darray(fh_ini, psdesc, iodesc2d, var2d, ierr)
     atm(mytile)%ps(is:ie,js:je) = var2d
     checksum=mpp_chksum(atm(mytile)%ps(is:ie,js:je))
     if (masterproc) write(6,*)'iam ',iam,'read PS is:ie,js:je CHKSUM=',checksum

     ! PHIS
!jt     call PIO_Read_Darray(fh_ini, phisdesc, iodesc2d, var2d, ierr)
!jt     atm(mytile)%phis(is:ie,js:je) = var2d
!jt     deallocate(var2d)
     
     
     allocate(var3d(is:ie,npz,js:je))
     var3d = 0._r8
     
!jt     ! OMEGA
!jt     call PIO_Read_Darray(fh_ini, omegadesc, iodesc3d, var3d, ierr)
!jt     Atm(mytile)%omga(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))

!jt     ! DELP
!jt     call PIO_Read_Darray(fh_ini, delpdesc, iodesc3d, var3d, ierr)
!jt     atm(mytile)%delp(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
     
     ! T
     call PIO_Read_Darray(fh_ini, Tdesc, iodesc3d, var3d, ierr)
     atm(mytile)%pt(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
     checksum=mpp_chksum(atm(mytile)%pt(is:ie,js:je,1:npz))
     if (masterproc) write(6,*)'iam ',iam,' reading T',' CHKSUM=',checksum
     
     if (pertlim .ne. 0.0_r8) then
        if(masterproc) then
           write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
                'by +/- ', pertlim, ' to initial temperature field'
        end if

        call random_seed(size=rndm_seed_sz)
        allocate(rndm_seed(rndm_seed_sz))

        ! seed random number generator based on global index
        ! (possibly include a flag to allow clock-based random seeding)
        allocate(glob_inddups(blksize))

        do i=is,ie
           do j=js,je
              indx=mylindex(i,j)
              rndm_seed = glob_inddups(indx)
              call random_seed(put=rndm_seed)
              do k=1,plev
                 call random_number(pertval)
                 pertval = 2.0_r8*pertlim*(0.5_r8 - pertval)
                 atm(mytile)%pt(i,j,k) = atm(mytile)%pt(i,j,k)*(1.0_r8 + pertval)
              end do
           end do
        end do
        deallocate(rndm_seed)
        deallocate(glob_inddups)
     end if


     
     ! V
     call PIO_Read_Darray(fh_ini, Vdesc, iodesc3d, var3d, ierr)
     atm(mytile)%va(is:ie,js:je,1:npz)=RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
     checksum=mpp_chksum(atm(mytile)%va(is:ie,js:je,1:npz))
     if (masterproc) write(6,*)'iam ',iam,' reading VA',' CHKSUM=',checksum
     ii=0
     tmparr(:,:)=0.
     do j = js, je
        do i = is, ie
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%va(i,j,1:npz)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,26) real(nan_count,r8), real(inf_count,r8), iam, iam
26        format("SHR_REPROSUM_CALC: v Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif
     ! U
     call PIO_Read_Darray(fh_ini, Udesc, iodesc3d, var3d, ierr)
     atm(mytile)%ua(is:ie,js:je,1:npz)   =RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
     checksum=mpp_chksum(atm(mytile)%ua(is:ie,js:je,1:npz))
     if (masterproc) write(6,*)'iam ',iam,' reading UA',' CHKSUM=',checksum
     if (iam.eq.960) write(6,*)'read ua ',atm(mytile)%ua(is:ie,js:je,1)
     ii=0
     tmparr(:,:)=0.
     do j = js, je
        do i = is, ie
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%ua(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,27) real(nan_count,r8), real(inf_count,r8), iam, iam
27        format("SHR_REPROSUM_CALC: ua Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif
     ! Q 
        m=1
        ierr = PIO_Inq_varid(fh_ini, trim(cnst_name(m)), Qdesc(m))
        call cam_pio_handle_error(ierr, subname//': cannot find '//trim(cnst_name(m)))
        call PIO_Read_Darray(fh_ini, Qdesc(m), iodesc3d, var3d, ierr)
        atm(mytile)%q(is:ie,js:je,1:npz,m) = RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
        checksum=mpp_chksum(atm(mytile)%q(is:ie,js:je,1:npz,m))
        if (masterproc) write(6,*)'iam ',iam,' reading Q is:ie,js:je',' CHKSUM=',checksum
     
     ! Read in or cold-initialize all the tracer fields
     ! Copy tracers defined on unstructured grid onto distributed FFSL grid
     ! Make sure tracers have at least minimum value
        
     allocate(dbuf3(blksize,plev,1))
     do m_cnst = 2, pcnst
        m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
        found = .false.

        if(cnst_read_iv(m_cnst)) then
           found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)),            &
                required=.false.)
        end if

        if(found) then
           if (masterproc) write(6,*)'Found ',trim(cnst_name(m_cnst)),' constituent number',m_cnst
           var3d=0._r8
           ierr = PIO_Inq_varid(fh_ini, trim(cnst_name(m_cnst)), Qdesc(m_cnst))
           call cam_pio_handle_error(ierr, subname//': cannot find '//trim(cnst_name(m_cnst)))
           call PIO_Read_Darray(fh_ini, Qdesc(m_cnst), iodesc3d, var3d, ierr)
           atm(mytile)%q(is:ie,js:je,1:npz,m_cnst_ffsl) = RESHAPE(var3d,(/ilen,jlen,npz/),ORDER=(/1,3,2/))
           checksum=mpp_chksum(atm(mytile)%q(is:ie,js:je,1:npz,m_cnst_ffsl))
           if (masterproc) write(6,*)'iam ',iam,' reading Q cnst ffsl is:ie,js:je ',m_cnst_ffsl,' m_cnst ',m_cnst,' CHKSUM=',checksum
        else
           dbuf3=0._r8
           if (masterproc) write(6,*)'Missing ',trim(cnst_name(m_cnst)),' constituent number',m_cnst,size(latvals_rad),size(dbuf3)
           if (masterproc) write(6,*)'Initializing ',trim(cnst_name(m_cnst)),'fv3 constituent number ',m_cnst_ffsl,' to default'
           call cnst_init_default(m_cnst, latvals_rad, lonvals_rad, dbuf3, pmask)
           do k=1, plev
              indx = 1
              do j = js, je
                 do i = is, ie
                    indx=mylindex(i,j)
                    atm(mytile)%q(i,j, k, m_cnst_ffsl) = max(qmin(m_cnst),dbuf3(indx,k,1))
                 end do
              end do
           end do
           checksum=mpp_chksum(atm(mytile)%q(is:ie,js:je,1:npz,m_cnst_ffsl))
           if (masterproc) write(6,*)'iam ',iam,' initializing Q cnst is:ie,js:je ',m_cnst_ffsl,' CHKSUM=',checksum
        end if

     end do ! pcnst

     deallocate(dbuf3)
     deallocate(var3d)
     deallocate(qdesc)
     
!!$     !jt check - This is incorrect as the values on the initial condition file hav3 staggered winds in lat/lon coordinates these need to be
!!$     !jt check   to D-grid coordinates.  For now just interpolate from a grid already read in.
!!$
!!$     ! US - For ew and ns grids PIO 1 can not read a single file point to
!!$     !      multiple processors.  For now I am reading twice, once to get the
!!$     !      uniq points distributed and the second time to read those points
!!$     !      required on more than one processor. (ie those points at the 
!!$     !      n/s e/w DGRID boundaries between processors
!!$     
!!$     allocate(var3d_ns(is:ie, npz,js:je+1))
!!$     allocate(var3d_ns_tmp(is:ie, npz,js:je+1))
!!$     allocate(var3d_ns_rst(is:ie, npz, js:je+1))
!!$     var3d_ns = 0._r8
!!$     var3d_ns_tmp = 0._r8
!!$     var3d_ns_rst = 0._r8
!!$     call PIO_Read_Darray(fh_ini, USdesc, iodesc3d_ns, var3d_ns_tmp, ierr)
!!$     var3d_ns=var3d_ns_tmp
!!$
!!$     ! US hack to read in duplicate points on adjacent processor
!!$     ! should fill in zeros in the decomposed global array after first read
!!$     call PIO_Read_Darray(fh_ini, USdesc, iodesc3d_ns_rst, var3d_ns_rst, ierr)
!!$     where(var3d_ns_tmp.eq.0)
!!$        var3d_ns = var3d_ns_rst
!!$     end where
!!$
!!$     atm(mytile)%u(is:ie,js:je+1,1:npz) = RESHAPE(var3d_ns,(/ilen,jlen+1,npz/),ORDER=(/1,3,2/))
!!$     
!!$     deallocate(var3d_ns)
!!$     deallocate(var3d_ns_rst)
!!$     deallocate(var3d_ns_tmp)
!!$     
!!$     allocate(var3d_ew(is:ie+1, npz, js:je))
!!$     allocate(var3d_ew_tmp(is:ie+1, npz, js:je))
!!$     allocate(var3d_ew_rst(is:ie+1, npz, js:je))
!!$     var3d_ew = 0._r8
!!$     var3d_ew_tmp = 0._r8
!!$     var3d_ew_rst = 0._r8
!!$     
!!$     ! VS
!!$     call PIO_Read_Darray(fh_ini, VSdesc, iodesc3d_ew, var3d_ew_tmp, ierr)
!!$     var3d_ew=var3d_ew_tmp
!!$     
!!$     
!!$     ! VS hack to read in duplicate points on adjacent processor
!!$     ! should fill in zeros in the decomposed global array after first read
!!$     call PIO_Read_Darray(fh_ini, VSdesc, iodesc3d_ew_rst, var3d_ew_rst, ierr)
!!$     where(var3d_ew_tmp.eq.0)
!!$        var3d_ew = var3d_ew_rst
!!$     end where
!!$     atm(mytile)%v(is:ie+1,js:je,1:npz) = RESHAPE(var3d_ew,(/ilen+1,jlen,npz/),ORDER=(/1,3,2/))
!!$     
!!$     deallocate(var3d_ew)
!!$     deallocate(var3d_ew_tmp)
!!$     deallocate(var3d_ew_rst)
     if (iam.eq.960) write(6,*)'ua 1 ',atm(mytile)%ua(is:ie,js:je,1)

     call a2d3djt(atm(mytile)%ua, atm(mytile)%va, atm(mytile)%u, atm(mytile)%v, is,  ie,  js,  je, isd, ied, jsd, jed, npx,npy, npz, atm(mytile)%gridstruct, atm(mytile)%domain)
     if (iam.eq.960) write(6,*)'ua 2 ',atm(mytile)%ua(is:ie,js:je,1)

     ii=0
     tmparr(:,:)=0.
     do j = jsd, jed+1
        do i = isd, ied
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%u(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,29) real(nan_count,r8), real(inf_count,r8), iam, iam
29        format("SHR_REPROSUM_CALC: after a2d3 %u contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif

     ii=0
     tmparr(:,:)=0.
     do j = js, je
        do i = is, ie
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%ua(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,30) real(nan_count,r8), real(inf_count,r8), iam, iam
30        format("SHR_REPROSUM_CALC: after a2d3 %ua Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif

     ! recreating a set of A winds from D winds using cubed_to_latlon to be consistent with what is done in the energy diagnostics. 
     call cubed_to_latlon(Atm(mytile)%u, Atm(mytile)%v, Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%gridstruct, &
          npx, npy, npz, 1, Atm(mytile)%gridstruct%grid_type, Atm(mytile)%domain, Atm(mytile)%gridstruct%nested, Atm(mytile)%flagstruct%c2l_ord, Atm(mytile)%bd)
     ii=0
     tmparr(:,:)=0.
     do j = jsd, jed+1
        do i = isd, ied
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%u(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,31) real(nan_count,r8), real(inf_count,r8), iam, iam
31        format("SHR_REPROSUM_CALC: after cubed %u contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif

     ii=0
     tmparr(:,:)=0.
     do j = js, je
        do i = is, ie
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%ua(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,32) real(nan_count,r8), real(inf_count,r8), iam, iam
32        format("SHR_REPROSUM_CALC: after cubed %ua Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif


     ! Put the error handling back the way it was
     call pio_seterrorhandling(fh_ini, err_handling)

  end if ! analytic_ic_active



  ! If a topo file is specified use it.  This will overwrite the PHIS set by the
  ! analytic IC option.
  !
  ! If using the physics grid then the topo file will be on that grid since its
  ! contents are primarily for the physics parameterizations, and the values of
  ! PHIS should be consistent with the values of sub-grid variability (e.g., SGH)
  ! which are computed on the physics grid.
  if (associated(fh_topo)) then

     ! We need to be able to see the PIO return values
     call pio_seterrorhandling(fh_topo, PIO_BCAST_ERROR, pio_errtype)

     fieldname = 'PHIS'
     if (dyn_field_exists(fh_topo, trim(fieldname))) then
        call read_dyn_var(trim(fieldname), fh_topo, 'ncol', phis_tmp)
     else
        call endrun(trim(subname)//': ERROR: Could not find PHIS field on input datafile')
     end if

     ! Put the error handling back the way it was
     call pio_seterrorhandling(fh_topo, pio_errtype)
  end if

  ! Process phis_tmp
  atm(mytile)%phis = 0.0_r8
  do j = js, je
     do i = is, ie
        indx = mylindex(i,j)
        atm(mytile)%phis(i,j) = phis_tmp(indx,1)
     end do
  end do
  checksum=mpp_chksum(atm(mytile)%phis(is:ie,js:je))
  if (masterproc) write(6,*)'iam ',iam,'reading PHIS is:ie,js:je CHKSUM=',checksum

  !                                                                                                                                                                                     
  ! initialize delp and mixing ratios
  !                                                                                                                                                                                     

  if (inic_wet) then

!
!  Initial condition should be consistent with fv3 dynamics and wouldn't normally need any adjustment here but I am using an 
!  interpolated fv initial condition with mixing ratios based off of (dry mass + vapor) and needs to be (dry mass + vapor + condensates)
!
     fv2fv3_mixratio=.true.

     if (fv2fv3_mixratio) then 
        allocate(pstmp(isd:ied,jsd:jed))
        pstmp(:,:) = atm(mytile)%ps(:,:)
        atm(mytile)%ps(:,:)=hyai(1)*ps0
        do k=1,pver
           do j = js, je
              do i = is, ie
                 ! this delp is (dry+vap) using the moist ps read in.
                 Atm(mytile)%delp(i, j, k) = (((hyai(k+1) - hyai(k))*ps0)         +                &
                      ((hybi(k+1) - hybi(k))*pstmp(i,j)))
                 delpdry=Atm(mytile)%delp(i,j,k)*(1.0_r8-Atm(mytile)%q(i,j,k,1))
                 do m=1,pcnst
                    m_ffsl=qsize_tracer_idx_cam2dyn(m)
                    if (cnst_type(m) == 'wet') then
                       tracermass(m_ffsl)=Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl)
                    else
                       tracermass(m_ffsl)=delpdry*Atm(mytile)%q(i,j,k,m_ffsl)
                    end if
                 end do
                 fv3_totwatermass=sum(tracermass(qsize_condensate_loading_idx(1:qsize_condensate_loading)))
                 fv3_airmass =  delpdry + fv3_totwatermass
                 Atm(mytile)%delp(i,j,k) = fv3_airmass
                 Atm(mytile)%q(i,j,k,1:pcnst) = tracermass(1:pcnst)/fv3_airmass
                 Atm(mytile)%ps(i,j)=Atm(mytile)%ps(i,j)+Atm(mytile)%delp(i, j, k)
                 ! check new tracermass all should be wet and consistent with mass before conversion to fv3 mix ratio
                 do m=1,pcnst
                    m_ffsl=qsize_tracer_idx_cam2dyn(m)
                    if (tracermass(m_ffsl).ne.0) then
                       reldif=(Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl)-tracermass(m_ffsl))/tracermass(m_ffsl)
                       if (reldif.gt.1.0e-15_r8) &
                       write(6,*)'mass inconsistency new, old, relative error=',iam,cnst_name(m),Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl),tracermass(m_ffsl),reldif
                    end if
                 end do
              end do
           end do
        end do
        deallocate(pstmp)
     end if
  else
     allocate(pstmp(isd:ied,jsd:jed))
     pstmp(:,:) = atm(mytile)%ps(:,:)
     atm(mytile)%ps(:,:)=hyai(1)*ps0
     do k=1,pver
        do j = js, je
           do i = is, ie
              ! this delp is (dry+vap) using the moist ps read in.
              delpdry = (((hyai(k+1) - hyai(k))*ps0)         +                &
                   ((hybi(k+1) - hybi(k))*pstmp(i,j)))
              do m=1,pcnst
                 tracermass(m)=delpdry*Atm(mytile)%q(i,j,k,m)
              end do
              fv3_totwatermass=sum(tracermass(qsize_condensate_loading_idx(1:qsize_condensate_loading)))
              fv3_airmass =  delpdry + fv3_totwatermass
              Atm(mytile)%delp(i,j,k) = fv3_airmass
              Atm(mytile)%q(i,j,k,1:pcnst) = tracermass(1:pcnst)/fv3_airmass
              Atm(mytile)%ps(i,j)=Atm(mytile)%ps(i,j)+Atm(mytile)%delp(i, j, k)
              ! check new tracermass all should be wet and consistent with mass before conversion to fv3 mix ratio
              do m=1,pcnst
                 m_ffsl=qsize_tracer_idx_cam2dyn(m)
                 reldif=(Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl)-tracermass(m_ffsl))/tracermass(m_ffsl)
                 if (reldif.gt.1.0e-15_r8) &
                      write(6,*)'mass inconsistency new, old, relative error=',iam,cnst_name(m),Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl),tracermass(m_ffsl),reldif
              end do
           end do
        end do
     end do
     deallocate(pstmp)
  end if

  !
  ! scale PS to achieve prescribed dry mass following FV dycore (dryairm.F90)
  !
  initial_global_ave_dry_ps = 98288.0_r8
  if (.not. associated(fh_topo)) initial_global_ave_dry_ps = 101325._r8-245._r8
  if (analytic_ic_active()     ) initial_global_ave_dry_ps = 0                  !do not scale psdry
  call set_dry_mass(Atm, initial_global_ave_dry_ps)
  

  !$omp parallel do private(i, j)
  do j=js,je
     do i=is,ie
        Atm(mytile)%pe(i,1,j)   = Atm(mytile)%ptop
        Atm(mytile)%pk(i,j,1)   = Atm(mytile)%ptop ** kappa
        Atm(mytile)%peln(i,1,j) = log(Atm(mytile)%ptop )
     enddo
  enddo
  
! TODO: consider swapping loops for better OMP performance (vertical dependency)
!$omp parallel do private(i,j,k)
  do j=js,je
     do k=1,pver
        do i=is,ie
           Atm(mytile)%pe(i,k+1,j)   = Atm(mytile)%pe(i,k,j) + Atm(mytile)%delp(i,j,k)
        enddo
     enddo
  enddo
  
!$omp parallel do private(i,j,k)
  do j=js,je
     do k=1,pver
        do i=is,ie
           !$$jt todo check to make sure kappa is same as fv3 (ie kappa_moist or normal kappa)
           Atm(mytile)%pk(i,j,k+1)= Atm(mytile)%pe(i,k+1,j) ** kappa
           Atm(mytile)%peln(i,k+1,j) = log(Atm(mytile)%pe(i,k+1,j))
           Atm(mytile)%pkz(i,j,k) = (Atm(mytile)%pk(i,j,k+1)-Atm(mytile)%pk(i,j,k))/(kappa*(Atm(mytile)%peln(i,k+1,j)-Atm(mytile)%peln(i,k,j)))
        enddo
     enddo
  enddo
!!  Initialize non hydrostatic variables if needed
  if (.not. Atm(mytile)%flagstruct%hydrostatic) then
     do k=1,npz
        do j=js,je
           do i=is,ie
              Atm(mytile)%w ( i,j,k ) = 0.
              Atm(mytile)%delz ( i,j,k ) = -rdgas/gravit*Atm(mytile)%pt( i,j,k ) * ( Atm(mytile)%peln( i,k+1,j ) - Atm(mytile)%peln( i,k,j ) )
           enddo
        enddo
     enddo
  end if

  ! once we've read or initialized all the fields we call update_domains to
  ! update the redundent columns in the dynamics

  call mpp_update_domains( Atm(mytile)%phis, Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%ps,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,   Atm(mytile)%domain, gridtype=DGRID_NE, complete=.true. )
  call mpp_update_domains( atm(mytile)%pt,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%delp,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%q,    Atm(mytile)%domain )
     ii=0
     tmparr(:,:)=0.
     do j = jsd, jed
        do i = isd, ied
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%u(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,33) real(nan_count,r8), real(inf_count,r8), iam, iam
33        format("SHR_REPROSUM_CALC: end u Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif
     ii=0
     tmparr(:,:)=0.
     do j = js, je
        do i = is, ie
           ii=ii+1
           tmparr(ii,:)=atm(mytile)%ua(i,j,1:pver)
        end do
     end do
    nan_check = any(shr_infnan_isnan(tmparr(:,:)))
    inf_check = any(shr_infnan_isinf(tmparr(:,:)))
    nan_count = count(shr_infnan_isnan(tmparr(:,:)))
    inf_count = count(shr_infnan_isinf(tmparr(:,:)))
    if (nan_check.or.inf_check) then
       if ((nan_count > 0) .or. (inf_count > 0)) then
          write(iulog,34) real(nan_count,r8), real(inf_count,r8), iam, iam
34        format("SHR_REPROSUM_CALC: end ua Input contains ",e12.5, &
               " NaNs and ", e12.5, " INFs on process ", i7, i7)
          call shr_sys_abort("shr_reprosum_calc ERROR: NaNs or INFs in input")
       endif
    endif

  ! Cleanup
  deallocate(pmask)
  deallocate(phis_tmp)
  
  if (associated(ldof)) then
     deallocate(ldof)
     nullify(ldof)
  end if
     if (iam.eq.960) write(6,*)'ua 3 end of inidat ',atm(mytile)%ua(is:ie,js:je,1)
  
end subroutine read_inidat


!=======================================================================

subroutine get_dyn_decomp(is,ie,js,je,nlev, datatype, iodesc)

    use pio,           only: io_desc_t, pio_initdecomp
    use cam_pio_utils, only: pio_subsystem
    use dyn_grid,      only: get_horiz_grid_dim_d

    implicit none

    integer, intent(in) :: nlev, datatype
    type(io_desc_t), intent(out) :: iodesc
    integer, pointer :: ldof(:)
    integer :: dimlens(2), dimcnt
    integer :: is,ie,js,je

    dimcnt = 1
    call get_horiz_grid_dim_d(dimlens(1))
    if (nlev .gt. 1) then
        dimlens(2) = nlev
        dimcnt = dimcnt + 1
    end if

    ldof => get_ldof(is,ie,js,je,nlev)
    call pio_initdecomp(pio_subsystem, datatype, dimlens(1:dimcnt), ldof, iodesc)

    deallocate(ldof)

end subroutine get_dyn_decomp

!=======================================================================

function get_ldof(is,ie,js,je,nlev) result(ldof)

    use dyn_grid,          only: get_horiz_grid_dim_d,mygindex

    implicit none

    integer, intent(in) :: is,ie,js,je
    integer, intent(in) :: nlev
    integer, pointer :: ldof(:)

    integer :: lcnt, i, j, k, ig, offset, hdim

    call get_horiz_grid_dim_d(hdim)

    lcnt = nlev * (ie - is + 1) * (je - js + 1)
    allocate(ldof(lcnt))

    ig = 1
    ldof(:) = 0
    do k = 1, nlev
        do j = js, je
            do i = is, ie
                offset = mygindex(i, j)
                ldof(ig) = offset + (k - 1) * hdim
                ig = ig + 1
            end do
        end do
    end do

end function get_ldof

!=======================================================================
  subroutine calc_tot_energy_dynamics(atm,nets,nete,tl,tl_qdp,outfld_name_suffix)
    use physconst,              only: gravit, cpair, rearth,omega
    use cam_history,            only: outfld, hist_fld_active
    use constituents,           only: cnst_get_ind
    use hycoef,                 only: hyai, ps0
    use constituents,           only: pcnst
    use pmgrid,                 only: plev
    use fv_mp_mod,              only: ng
    !------------------------------Arguments--------------------------------
    
    type(fv_atmos_type), pointer, intent(in) :: Atm(:)  
    integer          , intent(in) :: tl, tl_qdp,nets,nete!, n_fvm
    character*(*)    , intent(in) :: outfld_name_suffix ! suffix for "outfld" names
    
    !---------------------------Local storage-------------------------------
    
    real(kind=r8), allocatable :: se(:,:)                          ! Dry Static energy (J/m2)
    real(kind=r8), allocatable :: ke(:,:)                          ! kinetic energy    (J/m2)
    real(kind=r8), allocatable :: wv(:,:),wl(:,:),wi(:,:), &       ! col integ vap/liq/ice/rain/sno/graup(kg/m2)
                                  wr(:,:),ws(:,:),wg(:,:)          ! column integrated vapor       (kg/m2)
    real(kind=r8), allocatable :: tt(:,:)                          ! column integrated test tracer (kg/m2)
    real(kind=r8), allocatable :: dp(:,:,:)
    real(kind=r8), allocatable :: ps_local(:,:) 
    real(kind=r8) :: se_tmp
    real(kind=r8) :: ke_tmp
    real(kind=r8) :: wv_tmp,wl_tmp,wi_tmp,wr_tmp,ws_tmp,wg_tmp
    real(kind=r8) :: tt_tmp

    real(kind=r8), allocatable :: areasqrad(:,:)
    !
    ! global axial angular momentum (AAM) can be separated into one part (mr) associatedwith the relative motion 
    ! of the atmosphere with respect to the planets surface (also known as wind AAM) and another part (mo) 
    ! associated with the angular velocity OMEGA (2*pi/d, where d is the length of the day) of the planet 
    ! (also known as mass AAM)
    !
    real(kind=r8), allocatable :: mr(:,:)  ! wind AAM                        
    real(kind=r8), allocatable :: mo(:,:)  ! mass AAM
    real(kind=r8) :: mr_cnst, mo_cnst, cos_lat, mr_tmp, mo_tmp

    real(kind=r8) :: se_glob, ke_glob, wv_glob, wl_glob, wi_glob, wr_glob, ws_glob, wg_glob, tt_glob, mr_glob, mo_glob

    integer :: i,j,k,nq,m_cnst,n,idim,m_cnst_ffsl
    integer :: ixcldice, ixcldliq, ixtt,ixcldliq_ffsl,ixcldice_ffsl ! CLDICE, CLDLIQ and test tracer indices
    integer :: ixrain, ixsnow, ixgraupel,ixrain_ffsl, ixsnow_ffsl, ixgraupel_ffsl
    character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6,name_out7,name_out8,name_out9

    integer :: is,ie,js,je,isd,ied,jsd,jed,lchnk,ncol
    logical :: printglobals = .true.
    !-----------------------------------------------------------------------

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je
    isd = Atm(mytile)%bd%isd
    ied = Atm(mytile)%bd%ied
    jsd = Atm(mytile)%bd%jsd
    jed = Atm(mytile)%bd%jed

    se_glob = 0._r8; 
    ke_glob = 0._r8; 
    wv_glob = 0._r8; 
    wl_glob = 0._r8; 
    wi_glob = 0._r8; 
    wr_glob = 0._r8; 
    ws_glob = 0._r8; 
    wg_glob = 0._r8; 
    tt_glob = 0._r8; 
    mr_glob = 0._r8; 
    mo_glob = 0._r8; 

    allocate(se(is:ie,js:je))
    allocate(ke(is:ie,js:je))
    allocate(wv(is:ie,js:je))
    allocate(wl(is:ie,js:je))
    allocate(wi(is:ie,js:je))
    allocate(wr(is:ie,js:je))
    allocate(ws(is:ie,js:je))
    allocate(wg(is:ie,js:je))
    allocate(tt(is:ie,js:je))
    allocate(mr(is:ie,js:je))
    allocate(mo(is:ie,js:je))
    allocate(dp(is:ie,js:je,npz))
    allocate(ps_local(is:ie,js:je))

    name_out1 = 'SE_'   //trim(outfld_name_suffix)
    name_out2 = 'KE_'   //trim(outfld_name_suffix)
    name_out3 = 'WV_'   //trim(outfld_name_suffix)
    name_out4 = 'WL_'   //trim(outfld_name_suffix)
    name_out5 = 'WI_'   //trim(outfld_name_suffix)
    name_out6 = 'WR_'   //trim(outfld_name_suffix)
    name_out7 = 'WS_'   //trim(outfld_name_suffix)
    name_out8 = 'WG_'   //trim(outfld_name_suffix)
    name_out9 = 'TT_'   //trim(outfld_name_suffix)


    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4).or.hist_fld_active(name_out5).or.hist_fld_active(name_out6).or..true.) then
       if (qsize_condensate_loading>1) then
          call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
          call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
          call cnst_get_ind('RAINQM', ixrain, abort=.false.)
          call cnst_get_ind('SNOWQM', ixsnow, abort=.false.)
          call cnst_get_ind('GRAUQM', ixgraupel, abort=.false.)
       else
          ixcldliq = -1
          ixcldice = -1
          ixrain = -1
          ixsnow = -1 
          ixgraupel = -1
       end if
       
       call cnst_get_ind('TT_LW', ixtt, abort=.false.)
       if (ixtt.le.0) ixtt = -1
          
       dp(is:ie,js:je,1:npz)=Atm(mytile)%delp(is:ie,js:je,1:npz)
       !
       ! Compute frozen static energy in 3 parts:  KE, SE, and energy associated with vapor and liquid
       !

       se    = 0.0_r8
       ke    = 0.0_r8
       wv    = 0.0_r8
       wl    = 0.0_r8
       wi    = 0.0_r8
       wr    = 0.0_r8
       ws    = 0.0_r8
       wg    = 0.0_r8
       tt    = 0.0_r8

       do k = 1, plev
          do j=js,je
             do i = is, ie
                !
                ! make energy consistent with CAM physics (only water vapor and dry air in pressure)
                !
                !
                if ((.not.fv3_lcp_moist).and.(.not.fv3_lcv_moist).and.qsize_condensate_loading>1) then
                   dp(i,j,k) = Atm(mytile)%delp(i,j,k)             
                   ! adjust dp to include just dry + vap to use below 
                   do nq=2,qsize_condensate_loading                 
                      m_cnst_ffsl=qsize_condensate_loading_idx(nq)
                      dp(i,j,k) = dp(i,j,k)-Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_cnst_ffsl)
                   end do
                end if
                !
                ! kinetic energy using cam dp only dry + vap
                !
                ke_tmp   = 0.5_r8*(Atm(mytile)%va(i,j,k)**2+ Atm(mytile)%ua(i,j,k)**2)*dp(i,j,k)/gravit
                if (abs(Atm(mytile)%ua(i,j,k))>200.0.or.abs(Atm(mytile)%va(i,j,k))>200.0) then
                   write(*,*) "yyy",Atm(mytile)%ua(i,j,k),Atm(mytile)%va(i,j,k),outfld_name_suffix
                end if

                if (fv3_lcp_moist.or.fv3_lcv_moist) then
                   !
                   ! Internal energy formula including all condensates and corresponding heat capacities
                   !
                   ! Start with energy of dry air and add energy of condensates
                   dp(i,j,k) = Atm(mytile)%delp(i,j,k)             
                   do nq=1,qsize_condensate_loading                 
                      m_cnst_ffsl=qsize_condensate_loading_idx(nq)
                      dp(i,j,k) = dp(i,j,k)-Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_cnst_ffsl)
                   end do
                   se_tmp = cpair*dp(i,j,k)
                   do nq=1,qsize_condensate_loading
                      m_cnst_ffsl=qsize_condensate_loading_idx(nq)
                      if (fv3_lcp_moist) then
                         se_tmp = se_tmp+qsize_condensate_loading_cp(nq)*Atm(mytile)%q(i,j,k,m_cnst_ffsl)*Atm(mytile)%delp(i,j,k)
                      end if
                      if (fv3_lcv_moist) then
                         se_tmp = se_tmp+qsize_condensate_loading_cv(nq)*Atm(mytile)%q(i,j,k,m_cnst_ffsl)*Atm(mytile)%delp(i,j,k)
                      end if
                   end do
                   se_tmp = se_tmp*Atm(mytile)%pt(i,j,k)/gravit
                   ! reset dp to delp to use below fv3_lcp_moist case
                   dp(i,j,k) = Atm(mytile)%delp(i,j,k)             
                else
                   !
                   ! using CAM physics definition of internal energy
                   !
                   se_tmp   = cpair*Atm(mytile)%pt(i,j,k)*dp(i,j,k)/gravit
                end if
                wv_tmp   =  Atm(mytile)%q(i,j,k,1)*Atm(mytile)%delp(i,j,k)/gravit
                
                se(i,j) = se(i,j) + se_tmp
                ke(i,j) = ke(i,j) + ke_tmp
                wv(i,j) = wv(i,j) + wv_tmp

             end do
          end do
       end do
       
       do j=js,je
          do i = is,ie
             ps_local(i,j) =   Atm(mytile)%ptop+sum(dp(i,j,:))
          end do
       end do

       do j=js,je
          do i = is,ie
            se(i,j) = se(i,j) + Atm(mytile)%phis(i,j)*ps_local(i,j)/gravit
          end do
       end do
       
       ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
       
       if (ixcldliq > 1) then
          ixcldliq_ffsl = qsize_tracer_idx_cam2dyn(ixcldliq)
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   wl_tmp   = Atm(mytile)%q(i,j,k,ixcldliq_ffsl)*Atm(mytile)%delp(i,j,k)/gravit
                   wl   (i,j) = wl(i,j) + wl_tmp
                end do
             end do
          end do
       end if
       
       if (ixcldice > 1) then
          ixcldice_ffsl = qsize_tracer_idx_cam2dyn(ixcldice)
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   wi_tmp   = Atm(mytile)%q(i,j,k,ixcldice_ffsl)*Atm(mytile)%delp(i,j,k)/gravit
                   wi(i,j)    = wi(i,j) + wi_tmp
                end do
             end do
          end do
       end if

       if (ixrain > 1) then
          ixrain_ffsl = qsize_tracer_idx_cam2dyn(ixrain)
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   wr_tmp   = Atm(mytile)%q(i,j,k,ixrain_ffsl)*Atm(mytile)%delp(i,j,k)/gravit
                   wr   (i,j) = wr(i,j) + wr_tmp
                end do
             end do
          end do
       end if
       
       if (ixsnow > 1) then
          ixsnow_ffsl = qsize_tracer_idx_cam2dyn(ixsnow)
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   ws_tmp   = Atm(mytile)%q(i,j,k,ixsnow_ffsl)*Atm(mytile)%delp(i,j,k)/gravit
                   ws(i,j)    = ws(i,j) + ws_tmp
                end do
             end do
          end do
       end if

       if (ixgraupel > 1) then
          ixgraupel_ffsl = qsize_tracer_idx_cam2dyn(ixgraupel)
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   wg_tmp   = Atm(mytile)%q(i,j,k,ixgraupel_ffsl)*Atm(mytile)%delp(i,j,k)/gravit
                   wg(i,j)    = wg(i,j) + wg_tmp
                end do
             end do
          end do
       end if

       
       if (ixtt > 1) then
          do k = 1, plev
             do j = js, je
                do i = is, ie
                   tt_tmp   = Atm(mytile)%q(i,j,k,ixtt)*Atm(mytile)%delp(i,j,k)/gravit
                   tt   (i,j) = tt(i,j) + tt_tmp
                end do
             end do
          end do
       end if
       idim=ie-is+1
       do j=js,je
          ! Output energy diagnostics
          call outfld(name_out1  ,se(:,j)       ,idim, j)
          call outfld(name_out2  ,ke(:,j)       ,idim, j)
          call outfld(name_out3  ,wv(:,j)       ,idim, j)
          call outfld(name_out4  ,wl(:,j)       ,idim, j)
          call outfld(name_out5  ,wi(:,j)       ,idim, j)
          call outfld(name_out6  ,wr(:,j)       ,idim, j)
          call outfld(name_out7  ,ws(:,j)       ,idim, j)
          call outfld(name_out8  ,wg(:,j)       ,idim, j)
          if (ixtt > 1) call outfld(name_out9  ,tt(:,j)       ,idim, j)
       end do

       if (printglobals) then
          se_glob=g_sum(Atm(mytile)%domain, se(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          ke_glob=g_sum(Atm(mytile)%domain, ke(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wv_glob=g_sum(Atm(mytile)%domain, wv(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wl_glob=g_sum(Atm(mytile)%domain, wl(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wi_glob=g_sum(Atm(mytile)%domain, wi(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wr_glob=g_sum(Atm(mytile)%domain, wr(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          ws_glob=g_sum(Atm(mytile)%domain, ws(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wg_glob=g_sum(Atm(mytile)%domain, wg(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          if (ixtt > 1) &
               tt_glob=g_sum(Atm(mytile)%domain, tt(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          if (masterproc) then
      
             write(iulog, *) 'Using default g_sum quick'
             if (fv3_lcp_moist) then
                write(iulog, '(a,e25.17)') 'global wet static energy se_'//trim(outfld_name_suffix)//')            = ',se_glob
             else
                write(iulog, '(a,e25.17)') 'global dry static energy se_'//trim(outfld_name_suffix)//')            = ',se_glob
             end if
             write(iulog, '(a,e25.17)') 'global kinetic energy ke_'//trim(outfld_name_suffix)//')               = ',ke_glob
             write(iulog, '(a,e25.17)') 'global total energy se_plus_ke_'//trim(outfld_name_suffix)//')         = ',(ke_glob+se_glob)
             write(iulog, '(a,e25.17)') 'global column integrated vapor wv_'//trim(outfld_name_suffix)//'       = ',wv_glob
             write(iulog, '(a,e25.17)') 'global column integrated liquid wl_'//trim(outfld_name_suffix)//'      = ',wl_glob
             write(iulog, '(a,e25.17)') 'global column integrated ice wi_'//trim(outfld_name_suffix)//'         = ',wi_glob
             write(iulog, '(a,e25.17)') 'global column integrated liquid rain wr_'//trim(outfld_name_suffix)//'       = ',wr_glob
             write(iulog, '(a,e25.17)') 'global column integrated liquid snow ws_'//trim(outfld_name_suffix)//'      = ',ws_glob
             write(iulog, '(a,e25.17)') 'global column integrated graupel wg_'//trim(outfld_name_suffix)//'         = ',wg_glob
             if (ixtt > 1)      write(iulog, '(a,e25.17)') 'global column integrated test tracer tt_'//trim(outfld_name_suffix)//' = ',tt_glob
          end if
       end if
    end if
    
    !
    ! Axial angular momentum diagnostics
    !
    ! Code follows 
    !
    ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
    ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
    ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140, 
    ! doi:10.1002/2013MS000268
    !
    ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
    ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)
    !
    name_out1 = 'MR_'   //trim(outfld_name_suffix)
    name_out2 = 'MO_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2)) then



      mr_cnst = rearth**3/gravit
      mo_cnst = omega*rearth**4/gravit
      mr    = 0.0_r8
      mo    = 0.0_r8
      do k = 1, plev
         do j=js,je
            do i = is,ie
               cos_lat = cos(Atm(mytile)%gridstruct%agrid_64(i,j,2))
               mr_tmp   = mr_cnst*Atm(mytile)%ua(i,j,k)*Atm(mytile)%delp(i,j,k)*cos_lat
               mo_tmp   = mo_cnst*Atm(mytile)%delp(i,j,k)*cos_lat**2
               
               mr   (i,j) = mr(i,j) + mr_tmp
               mo   (i,j) = mo(i,j) + mo_tmp
            end do
         end do
      end do
      do j=js,je
         call outfld(name_out1  ,mr(is:ie,j)       ,idim,j)
         call outfld(name_out2  ,mo(is:ie,j)       ,idim,j)
      end do

      if (printglobals) then
         mr_glob=g_sum(Atm(mytile)%domain, mr(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
         mo_glob=g_sum(Atm(mytile)%domain, mo(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
         if (masterproc) then
            write(iulog, '(a,e25.17)') 'global column integrated wind AAM name_out1_'//trim(outfld_name_suffix)//'         = ',mr_glob
            write(iulog, '(a,e25.17)') 'global column integrated mass AAM name_out1_'//trim(outfld_name_suffix)//'         = ',mo_glob
         end if
      end if
   end if

  deallocate(ps_local)
  deallocate(dp)
  deallocate(se)
  deallocate(ke)
  deallocate(wv)
  deallocate(wl)
  deallocate(wi)
  deallocate(wr)
  deallocate(ws)
  deallocate(wg)
  deallocate(tt)
  deallocate(mr)
  deallocate(mo)
  end subroutine calc_tot_energy_dynamics

!========================================================================================

logical function dyn_field_exists(fh, fieldname, required)

   use pio,            only: file_desc_t, var_desc_t, PIO_inq_varid
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
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld

    ! Dummy arguments
    character(len=*),  intent(in)    :: fieldname
    type(file_desc_t), intent(inout) :: fh
    character(len=*),  intent(in)    :: dimname
    real(r8),          intent(inout) :: buffer(:, :)

    ! Local variables
    logical                  :: found
    !--------------------------------------------------------------------------

    buffer = 0.0_r8
    call infld(trim(fieldname), fh, dimname, 1, ldof_size, 1, 1, buffer,    &
         found, gridname='FFSL')
    if(.not. found) then
      call endrun('READ_DYN_FIELD_2D: Could not find '//trim(fieldname)//' field on input datafile')
    end if

    ! This code allows use of compiler option to set uninitialized values
    ! to NaN.  In that case infld can return NaNs where the element FFSL points
    ! are not "unique columns"
    where (isnan(buffer)) buffer = 0.0_r8

  end subroutine read_dyn_field_2d

!========================================================================================

  subroutine read_dyn_field_3d(fieldname, fh, dimname, buffer)
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld

    ! Dummy arguments
    character(len=*),  intent(in)    :: fieldname
    type(file_desc_t), intent(inout) :: fh
    character(len=*),  intent(in)    :: dimname
    real(r8),          intent(inout) :: buffer(:,:,:)

    ! Local variables
    logical                  :: found
    !--------------------------------------------------------------------------

    buffer = 0.0_r8
    call infld(fieldname, fh, 'ncol', 'lev', 1, ldof_size, 1, pver,     &
               1, 1, buffer, found, gridname='FFSL')
    if(.not. found) then
      call endrun('READ_DYN_FIELD_3D: Could not find '//trim(fieldname)//' field on input datafile')
    end if

    ! This code allows use of compiler option to set uninitialized values
    ! to NaN.  In that case infld can return NaNs where the element FFSL points
    ! are not "unique columns"
    where (isnan(buffer)) buffer = 0.0_r8

  end subroutine read_dyn_field_3d

!=========================================================================================

subroutine write_dyn_var(field,outfld_name,bd)

  use cam_history,            only: outfld
  implicit none
  
  type(fv_grid_bounds_type), intent(IN) :: bd
  real(r8), intent(IN)          :: field(bd%is:bd%ie,bd%js:bd%je)
  character*(*)    , intent(IN) :: outfld_name ! suffix for "outfld" names

  ! local variables
  integer              :: idim, j
  
  !----------------------------------------------------------------------------
  idim=bd%ie-bd%is+1
  do j=bd%js,bd%je
     ! Output energy diagnostics
     call outfld(trim(outfld_name)  ,field(bd%is:bd%ie,j)       ,idim, j)
  end do
  
end subroutine write_dyn_var

!=========================================================================================

subroutine set_dry_mass(atm,fixed_global_ave_dry_ps)
  
  !----------------------------------------------------------------------------
  
  use constituents,          only: pcnst, cnst_name, cnst_read_iv,qmin, cnst_type
  use cam_logfile,           only: iulog
  use hycoef,                only: hyai, hybi, ps0
  use dimensions_mod,        only: npz
  use dyn_grid,              only: mytile

  type (fv_atmos_type), intent(in),  pointer :: Atm(:)
  real (kind=r8), intent(in)                 :: fixed_global_ave_dry_ps

  ! local
  real (kind=r8)               :: global_ave_ps_inic,global_ave_dryps_inic,global_ave_dryps_scaled,global_ave_ps_new,global_ave_dryps_new
  real (r8), allocatable       :: factor(:,:,:),delpwet(:,:,:),delpdry(:,:,:),newdelp(:,:,:),psdry(:,:),psdry_scaled(:,:),psdry_new(:,:)
  integer                      :: i, j ,k, m,is,ie,js,je,m_ffsl

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  allocate(factor(is:ie,js:je,npz))
  allocate(delpdry(is:ie,js:je,npz))
  allocate(delpwet(is:ie,js:je,npz))
  allocate(newdelp(is:ie,js:je,npz))
  allocate(psdry(is:ie,js:je))
  allocate(psdry_scaled(is:ie,js:je))
  allocate(psdry_new(is:ie,js:je))


  if (fixed_global_ave_dry_ps == 0) return;
  
  ! get_global_ave_surface_pressure - must use bitwise sum (reproducable) to get with different decompositions.
  !  global_ave_ps_inic=g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
  global_ave_ps_inic=g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
  
  do k=1,pver
     do j = js, je
        do i = is, ie
           delpdry(i,j,k)=Atm(mytile)%delp(i,j,k)*(1.0_r8-sum(Atm(mytile)%q(i,j,k,qsize_condensate_loading_idx(1:qsize_condensate_loading))))
           delpwet(i,j,k)=Atm(mytile)%delp(i,j,k)-delpdry(i,j,k)
        end do
     end do
  end do
  !
  ! get psdry and scale it
  !
  do j = js, je
     do i = is, ie
        psdry(i,j) = hyai(1)*ps0 + sum(delpdry(i,j,:))
     end do
  end do

!  global_ave_dryps_inic=g_sum(Atm(mytile)%domain, psdry(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
  global_ave_dryps_inic=g_sum(Atm(mytile)%domain, psdry(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  psdry_scaled = psdry*(fixed_global_ave_dry_ps/global_ave_dryps_inic)

  global_ave_dryps_scaled=g_sum(Atm(mytile)%domain, psdry_scaled(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  !use adjusted psdry to calculate new dp_dry throughout atmosphere
  do k=1,pver
     do j = js, je
        do i = is, ie
           delpdry(i,j,k)=(hyai(k+1)-hyai(k))*ps0+&
                (hybi(k+1)-hybi(k))*psdry_scaled(i,j)
           ! new dp is adjusted dp + total watermass
           newdelp(i,j,k)=(delpdry(i,j,k)+delpwet(i,j,k))
           ! factor to conserve mass once using the new dp
           factor(i,j,k)=Atm(mytile)%delp(i,j,k)/newdelp(i,j,k)
           Atm(mytile)%delp(i,j,k)=newdelp(i,j,k)
        end do
     end do
  end do
  !
  ! all tracers wet in fv3 so conserve initial condition mass of 'wet' tracers (following se prim_set_dry)
  !
  do m=1,pcnst
     do k=1,pver
        do j = js, je
           do i = is, ie
              Atm(mytile)%q(i,j,k,m)=Atm(mytile)%q(i,j,k,m)*factor(i,j,k)
              Atm(mytile)%q(i,j,k,m)=max(qmin(m),Atm(mytile)%q(i,j,k,m))
           end do
        end do
     end do
  end do

  do j = js, je
     do i = is, ie
        Atm(mytile)%ps(i,j)=hyai(1)*ps0+sum(Atm(mytile)%delp(i, j, :))
        psdry_new(i,j)=hyai(1)*ps0+sum(delpdry(i, j, :))
     end do
  end do
!  global_ave_ps_new=g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1)
  global_ave_ps_new=   g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
  global_ave_dryps_new=g_sum(Atm(mytile)%domain, psdry_new(is:ie,js:je), is, ie, js, je, Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  if (masterproc) then
     write (iulog,*) "------ info from set_dry_mass -----------------------------------------------------------"
     write (iulog,*) "Scaling dry surface pressure to global average of = ",&
          fixed_global_ave_dry_ps/100.0_r8,"hPa"
     write (iulog,*) "Average surface pressure in initial condition = ",global_ave_ps_inic/100.0_r8,"hPa"
     write (iulog,*) "Average dry surface pressure in initial condition = ",global_ave_dryps_inic/100.0_r8,"hPa"
     write (iulog,*) "Average surface pressure after scaling = ",global_ave_ps_new/100.0_r8,"hPa"
     write (iulog,*) "Average dry surface pressure after scaling = ",global_ave_dryps_new/100.0_r8,"hPa"
     write (iulog,*) "Change in surface pressure           = ",&
          global_ave_ps_new-global_ave_ps_inic,"Pa"
     write (iulog,*) "Change in dry surface pressure           = ",&
          global_ave_dryps_new-global_ave_dryps_inic,"Pa"
     write (iulog,*) "Mixing ratios have been scaled so that total mass of tracer is conserved"
     write (iulog,*) "Total precipitable water before scaling = ", (global_ave_ps_inic-global_ave_dryps_inic)/gravit, '(kg/m**2)'
     write (iulog,*) "Total precipitable water after  scaling = ", (global_ave_ps_new-global_ave_dryps_new)/gravit, '(kg/m**2)'
     write (iulog,*) "------ end info from set_dry_mass -------------------------------------------------------"
  endif
  
  deallocate(factor)
  deallocate(delpdry)
  deallocate(delpwet)
  deallocate(newdelp)
  deallocate(psdry)
  deallocate(psdry_scaled)
  deallocate(psdry_new)

end subroutine set_dry_mass
!=========================================================================================

subroutine a2d3djt(ua, va, u, v, is,  ie,  js,  je, isd, ied, jsd, jed, npx,npy, npz, gridstruct, domain)

  use mpp_domains_mod,    only: mpp_update_domains,  DGRID_NE
  use fv_arrays_mod,      only: fv_grid_type

  integer, intent(in):: is,  ie,  js,  je
  integer, intent(in):: isd, ied, jsd, jed
  integer, intent(IN) :: npx,npy, npz
  real(r8), intent(inout):: u(isd:ied,  jsd:jed+1,npz)
  real(r8), intent(inout):: v(isd:ied+1,jsd:jed  ,npz)
  real(r8), intent(inout), dimension(isd:ied,jsd:jed,npz):: ua, va
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain

  ! local:
  real(r8) v3(is-1:ie+1,js-1:je+1,3)
  real(r8) ue(is-1:ie+1,js:je+1,3)    ! 3D winds at edges
  real(r8) ve(is:ie+1,js-1:je+1,  3)    ! 3D winds at edges
  real(r8), dimension(is:ie):: ut1, ut2, ut3
  real(r8), dimension(js:je):: vt1, vt2, vt3
  integer i, j, k, m, im2, jm2

  real(r8), pointer, dimension(:,:,:) :: vlon, vlat
  real(r8), pointer, dimension(:,:,:,:) :: es, ew
  real(r8), pointer, dimension(:) :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n

  es   => gridstruct%es
  ew   => gridstruct%ew
  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  edge_vect_w => gridstruct%edge_vect_w
  edge_vect_e => gridstruct%edge_vect_e
  edge_vect_s => gridstruct%edge_vect_s
  edge_vect_n => gridstruct%edge_vect_n

  call mpp_update_domains(ua, domain, complete=.false.)
  call mpp_update_domains(va, domain, complete=.true.)

  im2 = (npx-1)/2
  jm2 = (npy-1)/2

!$OMP parallel do default(none) shared(is,ie,js,je,npz,gridstruct,u,ua,v,va,  &
!$OMP                                  vlon,vlat,jm2,edge_vect_w,npx,edge_vect_e,im2, &
!$OMP                                  edge_vect_s,npy,edge_vect_n,es,ew)             &
!$OMP                          private(i,j,k,ut1, ut2, ut3, vt1, vt2, vt3, ue, ve, v3)
  do k=1, npz

     ! Compute 3D wind/tendency on A grid
     do j=js-1,je+1
        do i=is-1,ie+1
           v3(i,j,1) = ua(i,j,k)*vlon(i,j,1) + va(i,j,k)*vlat(i,j,1)
           v3(i,j,2) = ua(i,j,k)*vlon(i,j,2) + va(i,j,k)*vlat(i,j,2)
           v3(i,j,3) = ua(i,j,k)*vlon(i,j,3) + va(i,j,k)*vlat(i,j,3)
        enddo
     enddo

     ! Interpolate to cell edges
     do j=js,je+1
        do i=is-1,ie+1
           ue(i,j,1) = 0.5_r8*(v3(i,j-1,1) + v3(i,j,1))
           ue(i,j,2) = 0.5_r8*(v3(i,j-1,2) + v3(i,j,2))
           ue(i,j,3) = 0.5_r8*(v3(i,j-1,3) + v3(i,j,3))
        enddo
     enddo

     do j=js-1,je+1
        do i=is,ie+1
           ve(i,j,1) = 0.5_r8*(v3(i-1,j,1) + v3(i,j,1))
           ve(i,j,2) = 0.5_r8*(v3(i-1,j,2) + v3(i,j,2))
           ve(i,j,3) = 0.5_r8*(v3(i-1,j,3) + v3(i,j,3))
        enddo
     enddo

     ! --- E_W edges (for v-wind):
     if (.not. gridstruct%nested) then
        if ( is==1) then
           i = 1
           do j=js,je
              if ( j>jm2 ) then
                 vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
                 vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
                 vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
              else
                 vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1.-edge_vect_w(j))*ve(i,j,1)
                 vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1.-edge_vect_w(j))*ve(i,j,2)
                 vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1.-edge_vect_w(j))*ve(i,j,3)
              endif
           enddo
           do j=js,je
              ve(i,j,1) = vt1(j)
              ve(i,j,2) = vt2(j)
              ve(i,j,3) = vt3(j)
           enddo
        endif

        if ( (ie+1)==npx ) then
           i = npx
           do j=js,je
              if ( j>jm2 ) then
                 vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
                 vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
                 vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
              else
                 vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1.-edge_vect_e(j))*ve(i,j,1)
                 vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1.-edge_vect_e(j))*ve(i,j,2)
                 vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1.-edge_vect_e(j))*ve(i,j,3)
              endif
           enddo
           do j=js,je
              ve(i,j,1) = vt1(j)
              ve(i,j,2) = vt2(j)
              ve(i,j,3) = vt3(j)
           enddo
        endif
        ! N-S edges (for u-wind):
        if ( js==1) then
           j = 1
           do i=is,ie
              if ( i>im2 ) then
                 ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
                 ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
                 ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
              else
                 ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1.-edge_vect_s(i))*ue(i,j,1)
                 ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1.-edge_vect_s(i))*ue(i,j,2)
                 ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1.-edge_vect_s(i))*ue(i,j,3)
              endif
           enddo
           do i=is,ie
              ue(i,j,1) = ut1(i)
              ue(i,j,2) = ut2(i)
              ue(i,j,3) = ut3(i)
           enddo
        endif
        if ( (je+1)==npy ) then
           j = npy
           do i=is,ie
              if ( i>im2 ) then
                 ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
                 ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
                 ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
              else
                 ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1.-edge_vect_n(i))*ue(i,j,1)
                 ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1.-edge_vect_n(i))*ue(i,j,2)
                 ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1.-edge_vect_n(i))*ue(i,j,3)
              endif
           enddo
           do i=is,ie
              ue(i,j,1) = ut1(i)
              ue(i,j,2) = ut2(i)
              ue(i,j,3) = ut3(i)
           enddo
        endif

     endif ! .not. nested

     do j=js,je+1
        do i=is,ie
           u(i,j,k) =  ue(i,j,1)*es(1,i,j,1) +  &
                ue(i,j,2)*es(2,i,j,1) +  &
                ue(i,j,3)*es(3,i,j,1)
        enddo
     enddo
     do j=js,je
        do i=is,ie+1
           v(i,j,k) =  ve(i,j,1)*ew(1,i,j,2) +  &
                ve(i,j,2)*ew(2,i,j,2) +  &
                ve(i,j,3)*ew(3,i,j,2)
        enddo
     enddo
  enddo         ! k-loop

  call mpp_update_domains(u, v, domain, gridtype=DGRID_NE)

end subroutine a2d3djt

end module dyn_comp


