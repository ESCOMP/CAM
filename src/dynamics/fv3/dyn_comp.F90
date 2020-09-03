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



    use cam_abortutils,  only: endrun
    use cam_logfile,     only: iulog
    use constants_mod,   only: cp_air, kappa, rvgas, rdgas
    use constituents,    only: pcnst, cnst_name, cnst_longname, tottnam
    use dimensions_mod,  only: npx, npy, nlev, &
                               cnst_name_ffsl,cnst_longname_ffsl, &
                               fv3_lcp_moist,fv3_lcv_moist,qsize_tracer_idx_cam2dyn,fv3_scale_ttend
    use dyn_grid,        only: mytile
    use field_manager_mod, only: MODEL_ATMOS
    use fms_io_mod,      only: set_domain, nullify_domain
    use fv_arrays_mod,   only: fv_atmos_type, fv_grid_bounds_type
    use fv_grid_utils_mod,only: cubed_to_latlon, g_sum
    use fv_nesting_mod,  only: twoway_nesting
    use infnan,          only: isnan
    use mpp_domains_mod, only: mpp_update_domains, domain2D, DGRID_NE
    use mpp_mod,         only: mpp_set_current_pelist,mpp_pe
    use physconst,       only: gravit, cpair, rearth, omega, pi
    use ppgrid,          only: pver
    use shr_kind_mod,    only: r8 => shr_kind_r8, r4 => shr_kind_r4, i8 => shr_kind_i8
    use spmd_utils,      only: masterproc, masterprocid, mpicom, npes,iam
    use spmd_utils,      only: mpi_integer, mpi_logical
    use tracer_manager_mod,     only: get_tracer_index

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

type dyn_import_t
  type (fv_atmos_type),  pointer :: Atm(:) => null()
  integer,               pointer :: mygindex(:,:) => null()
  integer,               pointer :: mylindex(:,:) => null()
end type dyn_import_t

type dyn_export_t
  type (fv_atmos_type),  pointer :: Atm(:) => null()
end type dyn_export_t

! Private interfaces
interface read_dyn_var
  module procedure read_dyn_field_2d
  module procedure read_dyn_field_3d
end interface read_dyn_var

real(r8), public, allocatable, dimension(:,:,:) :: u_dt, v_dt, t_dt

!These are convenience variables for local use only, and are set to values in Atm%
real(r8) ::  zvir, dt_atmos_real

integer :: ldof_size

real(r8), allocatable,dimension(:,:,:)       :: se_dyn,ke_dyn,wv_dyn,wl_dyn,wi_dyn, &
                                                wr_dyn,ws_dyn,wg_dyn,tt_dyn,mo_dyn,mr_dyn

real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: deg2rad = pi / 180.0_r8

!=======================================================================
contains
!=======================================================================
subroutine dyn_readnl(nlfilename)

  ! Read dynamics namelist group from atm_in and write to fv3 input.nml file
  use namelist_utils,  only: find_group_name
  use constituents,    only: pcnst

  ! args
  character(len=*), intent(in) :: nlfilename

  ! Local variables
  integer                      :: unitn,unito, ierr,i,ios

  ! FV3 Namelist variables
  integer                      :: fv3_npes

  ! fv_core namelist variables - these namelist variables defined in fv3 library without fv3_

  integer            :: fv3_consv_te, fv3_dnats, fv3_fv_sg_adj, fv3_grid_type, &
                        fv3_hord_dp, fv3_hord_mt, fv3_hord_tm, fv3_hord_tr, fv3_hord_vt, &
                        fv3_io_layout(2), fv3_k_split, fv3_kord_mt, fv3_kord_tm, fv3_kord_tr, &
                        fv3_kord_wz, fv3_layout(2), fv3_n_split, fv3_n_sponge, fv3_na_init, &
                        fv3_ncnst, fv3_nord, fv3_npx, fv3_npy, fv3_npz, fv3_ntiles, &
                        fv3_nwat, fv3_print_freq

  real(r8)           :: fv3_beta, fv3_d2_bg, fv3_d2_bg_k1, fv3_d2_bg_k2, fv3_d4_bg, &
                        fv3_d_con, fv3_d_ext, fv3_dddmp, fv3_delt_max, fv3_ke_bg, &
                        fv3_rf_cutoff, fv3_tau, fv3_vtdm4

  logical            :: fv3_adjust_dry_mass, fv3_consv_am, fv3_do_sat_adj, fv3_do_vort_damp, &
                        fv3_dwind_2d, fv3_fill, fv3_fv_debug, fv3_fv_diag, fv3_hydrostatic, &
                        fv3_make_nh, fv3_no_dycore, fv3_range_warn

  ! fms_nml namelist variables - these namelist variables defined in fv3 library without fv3_

  character(len=256) :: fv3_clock_grain
  integer            :: fv3_domains_stack_size
  integer            :: fv3_stack_size
  logical            :: fv3_print_memory_usage

  character(len=256) :: inrec  ! first 80 characters of input record
  character(len=256) :: inrec2 ! left adjusted input record

  character(len = 20), dimension(5) :: group_names = (/  &
  "main_nml            ", &
  "fv_core_nml         ", &
  "surf_map_nml        ", &
  "test_case_nml       ", &
  "fms_nml             "/)

  namelist /fms_nml/          &
       fv3_clock_grain, &
       fv3_domains_stack_size, &
       fv3_print_memory_usage, &
       fv3_stack_size

  namelist /dyn_fv3_inparm/          &
       fv3_scale_ttend, &
       fv3_lcp_moist, &
       fv3_lcv_moist, &
       fv3_npes

  namelist /fv_core_nml/          &
       fv3_adjust_dry_mass,fv3_beta,fv3_consv_am,fv3_consv_te,fv3_d2_bg, &
       fv3_d2_bg_k1,fv3_d2_bg_k2,fv3_d4_bg,fv3_d_con,fv3_d_ext,fv3_dddmp, &
       fv3_delt_max,fv3_dnats,fv3_do_sat_adj,fv3_do_vort_damp,fv3_dwind_2d, &
       fv3_fill,fv3_fv_debug,fv3_fv_diag,fv3_fv_sg_adj,fv3_grid_type, &
       fv3_hord_dp,fv3_hord_mt,fv3_hord_tm,fv3_hord_tr,fv3_hord_vt, &
       fv3_hydrostatic,fv3_io_layout,fv3_k_split,fv3_ke_bg,fv3_kord_mt, &
       fv3_kord_tm,fv3_kord_tr,fv3_kord_wz,fv3_layout,fv3_make_nh, &
       fv3_n_split,fv3_n_sponge,fv3_na_init,fv3_ncnst,fv3_no_dycore, &
       fv3_nord,fv3_npx,fv3_npy,fv3_npz,fv3_ntiles,fv3_nwat, &
       fv3_print_freq,fv3_range_warn,fv3_rf_cutoff,fv3_tau, &
       fv3_vtdm4
  !--------------------------------------------------------------------------

  ! defaults for namelist variables not set by build-namelist
  fv3_npes         = npes

  if (masterproc) then
  ! Read the namelist (dyn_fv3_inparm)
     open( newunit=unitn, file=trim(NLFileName), status='old' )
     call find_group_name(unitn, 'dyn_fv3_inparm', status=ierr)
     if (ierr == 0) then
        read(unitn, dyn_fv3_inparm, iostat=ierr)
        if (ierr /= 0) then
           call endrun('dyn_readnl: ERROR reading dyn_fv3_inparm namelist')
        end if
     end if
     close(unitn)
  ! Read the namelist (fms_nml)
     open( newunit=unitn, file=trim(NLFileName), status='old' )
     call find_group_name(unitn, 'fms_nml', status=ierr)
     if (ierr == 0) then
        read(unitn, fms_nml, iostat=ierr)
        if (ierr /= 0) then
           call endrun('dyn_readnl: ERROR reading fms_nml namelist')
        end if
     end if
     close(unitn)
  ! Read the namelist (fv_core_nml)
     open( newunit=unitn, file=trim(NLFileName), status='old' )
     call find_group_name(unitn, 'fv_core_nml', status=ierr)
     if (ierr == 0) then
        read(unitn, fv_core_nml, iostat=ierr)
        if (ierr /= 0) then
           call endrun('dyn_readnl: ERROR reading fv_core_nml namelist')
        end if
     end if
     close(unitn)
  end if

  ! Broadcast namelist values to all PEs
  call MPI_bcast(fv3_npes, 1, mpi_integer, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_scale_ttend, 1, mpi_logical, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_hydrostatic, 1, mpi_logical, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_lcv_moist, 1, mpi_logical, masterprocid, mpicom, ierr)
  call MPI_bcast(fv3_lcp_moist, 1, mpi_logical, masterprocid, mpicom, ierr)

  if ((fv3_lcp_moist.eqv.fv3_lcv_moist) .and. (fv3_lcv_moist.eqv..true.)) then
     call endrun('dyn_readnl: fv3_lcp_moist and fv3_lcv_moist can not both be true')
  endif

  if (fv3_npes <= 0) then
     call endrun('dyn_readnl: ERROR: fv3_npes must be > 0')
  end if

  ! Non-hydrostatic runs not currently supported
  if (.not.fv3_hydrostatic) &
       call endrun('dyn_readnl: ERROR FV3 Non-hydrostatic option is not supported, set namelist fv3_hydrostatic = .true.')

  !
  ! write fv3 dycore namelist options to log
  !
  if (masterproc) then
     write (iulog,*) 'FV3 dycore Options: '
     write (iulog,*) '  fv3_adjust_dry_mass       = ',fv3_adjust_dry_mass
     write (iulog,*) '  fv3_beta                  = ',fv3_beta
     write (iulog,*) '  fv3_clock_grain           = ',trim(fv3_clock_grain)
     write (iulog,*) '  fv3_consv_am              = ',fv3_consv_am
     write (iulog,*) '  fv3_consv_te              = ',fv3_consv_te
     write (iulog,*) '  fv3_d2_bg                 = ',fv3_d2_bg
     write (iulog,*) '  fv3_d2_bg_k1              = ',fv3_d2_bg_k1
     write (iulog,*) '  fv3_d2_bg_k2              = ',fv3_d2_bg_k2
     write (iulog,*) '  fv3_d4_bg                 = ',fv3_d4_bg
     write (iulog,*) '  fv3_d_con                 = ',fv3_d_con
     write (iulog,*) '  fv3_d_ext                 = ',fv3_d_ext
     write (iulog,*) '  fv3_dddmp                 = ',fv3_dddmp
     write (iulog,*) '  fv3_delt_max              = ',fv3_delt_max
     write (iulog,*) '  fv3_dnats                 = ',fv3_dnats
     write (iulog,*) '  fv3_do_sat_adj            = ',fv3_do_sat_adj
     write (iulog,*) '  fv3_do_vort_damp          = ',fv3_do_vort_damp
     write (iulog,*) '  fv3_dwind_2d              = ',fv3_dwind_2d
     write (iulog,*) '  fv3_fill                  = ',fv3_fill
     write (iulog,*) '  fv3_fv_debug              = ',fv3_fv_debug
     write (iulog,*) '  fv3_fv_diag               = ',fv3_fv_diag
     write (iulog,*) '  fv3_fv_sg_adj             = ',fv3_fv_sg_adj
     write (iulog,*) '  fv3_grid_type             = ',fv3_grid_type
     write (iulog,*) '  fv3_hord_dp               = ',fv3_hord_dp
     write (iulog,*) '  fv3_hord_mt               = ',fv3_hord_mt
     write (iulog,*) '  fv3_hord_tm               = ',fv3_hord_tm
     write (iulog,*) '  fv3_hord_tr               = ',fv3_hord_tr
     write (iulog,*) '  fv3_hord_vt               = ',fv3_hord_vt
     write (iulog,*) '  fv3_hydrostatic           = ',fv3_hydrostatic
     write (iulog,*) '  fv3_io_layout             = ',fv3_io_layout
     write (iulog,*) '  fv3_k_split               = ',fv3_k_split
     write (iulog,*) '  fv3_ke_bg                 = ',fv3_ke_bg
     write (iulog,*) '  fv3_kord_mt               = ',fv3_kord_mt
     write (iulog,*) '  fv3_kord_tm               = ',fv3_kord_tm
     write (iulog,*) '  fv3_kord_tr               = ',fv3_kord_tr
     write (iulog,*) '  fv3_kord_wz               = ',fv3_kord_wz
     write (iulog,*) '  fv3_layout                = ',fv3_layout
     write (iulog,*) '  fv3_lcp_moist             = ',fv3_lcp_moist
     write (iulog,*) '  fv3_lcv_moist             = ',fv3_lcv_moist
     write (iulog,*) '  fv3_make_nh               = ',fv3_make_nh
     write (iulog,*) '  fv3_n_split               = ',fv3_n_split
     write (iulog,*) '  fv3_n_sponge              = ',fv3_n_sponge
     write (iulog,*) '  fv3_na_init               = ',fv3_na_init
     write (iulog,*) '  fv3_ncnst                 = ',fv3_ncnst
     write (iulog,*) '  fv3_no_dycore             = ',fv3_no_dycore
     write (iulog,*) '  fv3_nord                  = ',fv3_nord
     write (iulog,*) '  fv3_npx                   = ',fv3_npx
     write (iulog,*) '  fv3_npy                   = ',fv3_npy
     write (iulog,*) '  fv3_npz                   = ',fv3_npz
     write (iulog,*) '  fv3_ntiles                = ',fv3_ntiles
     write (iulog,*) '  fv3_nwat                  = ',fv3_nwat
     write (iulog,*) '  fv3_print_freq            = ',fv3_print_freq
     write (iulog,*) '  fv3_domains_stack_size    = ',fv3_domains_stack_size
     write (iulog,*) '  fv3_range_warn            = ',fv3_range_warn
     write (iulog,*) '  fv3_rf_cutoff             = ',fv3_rf_cutoff
     write (iulog,*) '  fv3_scale_ttend           = ',fv3_scale_ttend
     write (iulog,*) '  fv3_stack_size            = ',fv3_stack_size
     write (iulog,*) '  fv3_tau                   = ',fv3_tau
     write (iulog,*) '  fv3_vtdm4                 = ',fv3_vtdm4
  end if

  ! Create the input.nml namelist needed by the fv3dycore.
  ! Read strings one at a time from the fv3 namelist groups,
  ! strip off the leading 'fv3_' from the variable names and write to input.nml.
  ! This could be replaced by also by writing to the internal namelist file

  if (masterproc) then

     write(iulog,*) 'Creating fv3 input.nml file from atm_in fv3_xxx namelist parameters'
     ! Read the namelist (main_nml)
     ! open the file input.nml
     ! overwrite file if it exists.
     open( newunit=unito, file='input.nml', status='replace' )

     open( newunit=unitn, file=trim(NLFileName), status='old' )

     do i=1,SIZE(group_names(:))
        rewind(unitn)
        call find_group_name(unitn, trim(group_names(i)), status=ierr)

        if (ierr == 0) then ! Found it. Copy each line to input.nml until '/' is encountered.

           ! write group name to input.nml
           read(unitn, '(a)', iostat=ios, end=100) inrec
           if (ios /= 0) call endrun('ERROR: dyn_readnl - error reading fv3 namelist')
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
     close(unito)
  end if
  return
100 continue
  call endrun('ERROR: dyn_readnl: End of file encountered while reading fv3 namelist groups')

end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_register()

   ! These fields are computed by the dycore and passed to the physics via the
   ! physics buffer.

end subroutine dyn_register

!=============================================================================================

subroutine dyn_init(dyn_in, dyn_out)

  ! DESCRIPTION: Initialize the FV dynamical core

  ! Initialize FV dynamical core state variables


  use cam_control_mod, only: initial_run
  use cam_history,     only: addfld, horiz_only
  use cam_history,     only: register_vector_field
  use cam_pio_utils,   only: clean_iodesc_list
  use dyn_grid,        only: Atm,mygindex,mylindex
  use fv_diagnostics_mod, only: fv_diag_init
  use fv_mp_mod,       only: fill_corners, YDir, switch_current_Atm
  use infnan,          only: inf, assignment(=)
  use physconst,       only: cpwv, cpliq, cpice
  use physconst,          only: thermodynamic_active_species_num, dry_air_species_num, thermodynamic_active_species_idx
  use physconst,          only: thermodynamic_active_species_idx_dycore, rair, cpair
  use tracer_manager_mod,     only: register_tracers

  ! arguments:
   type (dyn_import_t),     intent(out) :: dyn_in
   type (dyn_export_t),     intent(out) :: dyn_out

   ! Locals
   character(len=*), parameter :: sub='dyn_init'
   real(r8)                    :: alpha


   real(r8), pointer, dimension(:,:)            :: fC,f0   ! Coriolis parameters
   real(r8), pointer, dimension(:,:,:)          :: grid,agrid,delp
   logical, pointer :: cubed_sphere
   type(domain2d), pointer     :: domain
   integer                     :: i,j,m

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
   integer :: is,isd,ie,ied,js,jsd,je,jed
   integer :: fv3idx,idx

   integer        :: unito
   integer, parameter  :: ndiag = 5
   integer             :: ncnst, pnats, num_family, nt_prog
   character(len=128) :: errmsg
   logical            :: wet_thermo_species
   !-----------------------------------------------------------------------

   ! Setup the condensate loading arrays and fv3/cam tracer mapping and
   ! finish initializing fv3 by allocating the tracer arrays in the fv3 atm structure

   allocate(qsize_tracer_idx_cam2dyn(pcnst))
   qsize_tracer_idx_cam2dyn(:)=-1
   allocate(cnst_name_ffsl(pcnst))     ! constituent names for ffsl tracers
   allocate(cnst_longname_ffsl(pcnst)) ! long name of constituents for ffsl tracers


   ! set up the condensate loading array
   if (thermodynamic_active_species_num - dry_air_species_num > 6) then
     call endrun(subname//': fv3_thermodynamic_active_species_num is limited to 6 wet condensates')
   end if

   !For FV3 Q must be the first species in the fv3 tracer array followed by wet constituents
   idx=1
   do m=1,pcnst
      if ( trim(cnst_name(m)) == 'Q'.or.&
           trim(cnst_name(m)) == 'CLDLIQ'.or.&
           trim(cnst_name(m)) == 'CLDICE'.or.&
           trim(cnst_name(m)) == 'RAINQM'.or.&
           trim(cnst_name(m)) == 'SNOWQM'.or.&
           trim(cnst_name(m)) == 'GRAUQM') then
         idx=idx+1
         wet_thermo_species=any(thermodynamic_active_species_idx(dry_air_species_num+1:thermodynamic_active_species_num)==m)
         select case ( trim(cnst_name(m)) )
         case ( 'Q' )
            idx=idx-1
            cnst_name_ffsl(1)='sphum'
            cnst_longname_ffsl(1) = cnst_longname(m)
            qsize_tracer_idx_cam2dyn(m)  = 1
            if (wet_thermo_species) thermodynamic_active_species_idx_dycore(1)=1
         case ( 'CLDLIQ' )
            cnst_name_ffsl(idx)='liq_wat'
         case ( 'CLDICE' )
            cnst_name_ffsl(idx)='ice_wat'
         case ( 'RAINQM' )
            cnst_name_ffsl(idx)='rainwat'
         case ( 'SNOWQM' )
            cnst_name_ffsl(idx)='snowwat'
         case ( 'GRAUQM' )
            cnst_name_ffsl(idx)='graupel'
         end select

         if (trim(cnst_name(m))/='Q') then
            if (wet_thermo_species) thermodynamic_active_species_idx_dycore(idx)=idx
            cnst_longname_ffsl(idx) = cnst_longname(m)
            qsize_tracer_idx_cam2dyn(m)  = idx
         end if
      end if
   end do

   do m=1,pcnst
      if ( trim(cnst_name(m)) /= 'Q'.and.&
           trim(cnst_name(m)) /= 'CLDLIQ'.and.&
           trim(cnst_name(m)) /= 'CLDICE'.and.&
           trim(cnst_name(m)) /= 'RAINQM'.and.&
           trim(cnst_name(m)) /= 'SNOWQM'.and.&
           trim(cnst_name(m)) /= 'GRAUQM') then
         idx=idx+1
         cnst_name_ffsl(idx)=cnst_name(m)
         cnst_longname_ffsl(idx) = cnst_longname(m)
         qsize_tracer_idx_cam2dyn(m) = idx
      end if
   end do

   if (masterproc) then

      write(iulog,*) 'Creating field_table file to load tracer fields into fv3'
      ! overwrite file if it exists.
      open( newunit=unito, file='field_table', status='replace' )
      do i=1,pcnst
         write(unito, '(a,a,a)') '"tracer" "atmos_mod" "'//trim(cnst_name_ffsl(i))//'" /'
      end do
      close(unito)
   end if
   !---------must make sure the field_table file is written before reading across processors
   call mpibarrier (mpicom)
   call register_tracers (MODEL_ATMOS, ncnst, nt_prog, pnats, num_family)
   if (ncnst /= pcnst) then
      call endrun(subname//': ERROR: FMS tracer Manager has inconsistent tracer numbers')
   endif

   do m=1,pcnst
      !  just check condensate loading tracers as they are mapped above
      if(qsize_tracer_idx_cam2dyn(m) <= thermodynamic_active_species_num-dry_air_species_num) then
         fv3idx  = get_tracer_index (MODEL_ATMOS, cnst_name_ffsl(qsize_tracer_idx_cam2dyn(m)) )
         if (fv3idx /= qsize_tracer_idx_cam2dyn(m)) then
            write(errmsg,*) subname//': Physics index ',m,'and FV3 tracer index',fv3idx,' are inconsistent'
            call endrun(errmsg)
         end if
      end if
   end do

   is = Atm(mytile)%bd%is
   ie = Atm(mytile)%bd%ie
   js = Atm(mytile)%bd%js
   je = Atm(mytile)%bd%je
   isd = Atm(mytile)%bd%isd
   ied = Atm(mytile)%bd%ied
   jsd = Atm(mytile)%bd%jsd
   jed = Atm(mytile)%bd%jed

   ! Data initialization
   dyn_in%Atm  => Atm
   dyn_in%mygindex  => mygindex
   dyn_in%mylindex  => mylindex
   dyn_out%Atm => Atm

   allocate(u_dt(isd:ied,jsd:jed,nlev))
   allocate(v_dt(isd:ied,jsd:jed,nlev))
   allocate(t_dt(isd:ied,jsd:jed,nlev))
   u_dt(:,:,:) = 0._r8
   v_dt(:,:,:) = 0._r8
   t_dt(:,:,:) = 0._r8

   fC    => atm(mytile)%gridstruct%fC
   f0    => atm(mytile)%gridstruct%f0
   grid  => atm(mytile)%gridstruct%grid_64
   agrid => atm(mytile)%gridstruct%agrid_64
   domain=> Atm(mytile)%domain
   cubed_sphere => atm(mytile)%gridstruct%cubed_sphere
   delp =>  Atm(mytile)%delp

   ! initialize Coriolis parameters which are used in sw_core.
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

   delp(isd:is-1,jsd:js-1,1:nlev)=0._r8
   delp(isd:is-1,je+1:jed,1:nlev)=0._r8
   delp(ie+1:ied,jsd:js-1,1:nlev)=0._r8
   delp(ie+1:ied,je+1:jed,1:nlev)=0._r8

   if (initial_run) then

      ! Read in initial data
      call read_inidat(dyn_in)
      call clean_iodesc_list()

   end if

   call switch_current_Atm(Atm(mytile))
   call set_domain ( Atm(mytile)%domain )

   ! Forcing from physics on the FFSL grid
   call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term on FFSL grid',     gridname='FFSLHIST')
   call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term on FFSL grid',gridname='FFSLHIST')
   call register_vector_field('FU', 'FV')
   call addfld ('FT',  (/ 'lev' /), 'A', 'K/s', 'Temperature forcing term on FFSL grid',gridname='FFSLHIST')

   do m = 1, pcnst
     call addfld ('F'//trim(cnst_name_ffsl(m))//'_ffsl',  (/ 'lev' /), 'I', 'kg/kg/s',   &
          trim(cnst_longname(m))//' mixing ratio forcing term (q_new-q_old) on FFSL grid', gridname='FFSLHIST')
     call addfld(tottnam(m),(/ 'lev' /),'A','kg/kg/s', &
          trim(cnst_name_ffsl(m))//' horz + vert + fixer tendency ',  &
                  gridname='FFSLHIST')
   end do

   ! Energy diagnostics and axial angular momentum diagnostics
   do istage = 1,SIZE(stage)
      do ivars=1,SIZE(vars)
         write(str1,*) TRIM(ADJUSTL(vars(ivars))),TRIM(ADJUSTL("_")),TRIM(ADJUSTL(stage(istage)))
         write(str2,*) TRIM(ADJUSTL(vars_descriptor(ivars))),&
             TRIM(ADJUSTL(" ")),TRIM(ADJUSTL(stage_txt(istage)))
         write(str3,*) TRIM(ADJUSTL(vars_unit(ivars)))
         call addfld (TRIM(ADJUSTL(str1)),horiz_only,'A',TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)), &
              gridname='FFSLHIST')
      end do
   end do

   allocate(se_dyn(is:ie,js:je,ndiag))
   allocate(ke_dyn(is:ie,js:je,ndiag))
   allocate(wv_dyn(is:ie,js:je,ndiag))
   allocate(wl_dyn(is:ie,js:je,ndiag))
   allocate(wi_dyn(is:ie,js:je,ndiag))
   allocate(wr_dyn(is:ie,js:je,ndiag))
   allocate(ws_dyn(is:ie,js:je,ndiag))
   allocate(wg_dyn(is:ie,js:je,ndiag))
   allocate(tt_dyn(is:ie,js:je,ndiag))
   allocate(mr_dyn(is:ie,js:je,ndiag))
   allocate(mo_dyn(is:ie,js:je,ndiag))


end subroutine dyn_init

!=======================================================================

subroutine dyn_run(dyn_state)

  ! DESCRIPTION: Driver for the NASA finite-volume dynamical core


  use dimensions_mod,         only: nlev
  use dyn_grid,               only: p_split,grids_on_this_pe
  use fv_control_mod,         only: ngrids
  use fv_dynamics_mod,        only: fv_dynamics
  use fv_sg_mod,              only: fv_subgrid_z
  use physconst,              only: thermodynamic_active_species_num, thermodynamic_active_species_idx_dycore, &
                                    thermodynamic_active_species_cp,thermodynamic_active_species_cv,dry_air_species_num
  use time_manager,           only: get_step_size
  use tracer_manager_mod,     only: get_tracer_index, NO_TRACER

  ! Arguments
  type (dyn_export_t), intent(inout) :: dyn_state

  ! Locals
  integer :: psc,idim
  integer :: w_diff, nt_dyn
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

  dt_atmos_real=get_step_size()

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

     call fv_dynamics(npx, npy, nlev, pcnst, Atm(mytile)%ng, dt_atmos_real/real(abs(p_split), r8),&
          Atm(mytile)%flagstruct%consv_te, Atm(mytile)%flagstruct%fill,  &
          Atm(mytile)%flagstruct%reproduce_sum, kappa, cp_air, zvir,&
          Atm(mytile)%ptop, Atm(mytile)%ks, pcnst,                          &
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
          pcnst,thermodynamic_active_species_num,dry_air_species_num, &
          thermodynamic_active_species_idx_dycore, qsize_tracer_idx_cam2dyn, &
          thermodynamic_active_species_cp,thermodynamic_active_species_cv, se_dyn, ke_dyn, wv_dyn,wl_dyn, &
          wi_dyn,wr_dyn,ws_dyn,wg_dyn,tt_dyn,mo_dyn,mr_dyn,gravit,cpair,rearth,omega,fv3_lcp_moist,&
          fv3_lcv_moist)
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
     nt_dyn = pcnst
     if ( w_diff /= NO_TRACER ) then
        nt_dyn = pcnst - 1
     endif
     call fv_subgrid_z(isd, ied, jsd, jed, isc, iec, jsc, jec, nlev, &
          nt_dyn, dt_atmos_real, Atm(mytile)%flagstruct%fv_sg_adj,      &
          Atm(mytile)%flagstruct%nwat, Atm(mytile)%delp, Atm(mytile)%pe,     &
          Atm(mytile)%peln, Atm(mytile)%pkz, Atm(mytile)%pt, Atm(mytile)%q,       &
          Atm(mytile)%ua, Atm(mytile)%va, Atm(mytile)%flagstruct%hydrostatic,&
          Atm(mytile)%w, Atm(mytile)%delz, u_dt, v_dt, t_dt, Atm(mytile)%flagstruct%n_sponge)
  endif

#if ( defined CALC_ENERGY )
  call calc_tot_energy_dynamics(atm,'dBF')
#endif

end subroutine dyn_run

!=======================================================================

subroutine dyn_final(dyn_in, dyn_out, restart_file)

  ! Arguments
  type (dyn_import_t),      intent(inout) :: dyn_in
  type (dyn_export_t),      intent(inout) :: dyn_out
  character(len=*),optional,intent(in)    :: restart_file

  !----------------------------------------------------------------------------

  deallocate( u_dt, v_dt, t_dt)

end subroutine dyn_final

!=============================================================================================
! Private routines
!=============================================================================================

subroutine read_inidat(dyn_in)

  use cam_control_mod,       only: simple_phys
  use inic_analytic,         only: analytic_ic_active, analytic_ic_set_ic
  use dyn_tests_utils,       only: vc_moist_pressure,vc_dry_pressure
  use dimensions_mod,        only: nlev
  use constituents,          only: pcnst, cnst_is_a_water_species
  use physconst,             only: thermodynamic_active_species_num, dry_air_species_num, thermodynamic_active_species_idx_dycore
  use pio,                   only: file_desc_t, pio_seterrorhandling, pio_bcast_error
  use ppgrid,                only: pver
  use cam_abortutils,        only: endrun
  use constituents,          only: pcnst, cnst_name, cnst_read_iv,qmin, cnst_type
  use const_init,            only: cnst_init_default
  use cam_initfiles,         only: initial_file_get_id, topo_file_get_id, pertlim
  use cam_grid_support,      only: cam_grid_id, cam_grid_get_gcid, iMap, &
                                   cam_grid_get_latvals, cam_grid_get_lonvals
  use cam_history_support,   only: max_fieldname_len
  use hycoef,                only: hyai, hybi, ps0

  ! Arguments:
  type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

  ! Locals:
  logical :: found

  character(len = 40) :: fieldname,fieldname2

  integer :: i, j, k, m, n

  type(file_desc_t),     pointer :: fh_topo       => null()
  type(fv_atmos_type),   pointer :: Atm(:)        => null()
  integer,               pointer :: mylindex(:,:) => null()
  integer,               pointer :: mygindex(:,:) => null()
  type(file_desc_t)              :: fh_ini


  character(len=*), parameter      :: subname='READ_INIDAT'

  ! Variables for analytic initial conditions
  integer, allocatable, dimension(:)        :: glob_ind, m_ind,rndm_seed
  integer                                   :: is,ie,js,je,isd,ied,jsd,jed
  integer                                   :: blksize
  integer                                   :: indx
  integer                                   :: err_handling
  integer                                   :: m_cnst,m_cnst_ffsl
  integer                                   :: m_ffsl
  integer                                   :: ilen,jlen
  integer                                   :: num_wet_species! (wet species are first tracers in FV3 tracer array)
  integer                                   :: pio_errtype
  integer                                   :: rndm_seed_sz
  integer                                   :: vcoord
  real(r8), pointer, dimension(:)           :: latvals_deg(:)
  real(r8), pointer, dimension(:)           :: lonvals_deg(:)
  real(r8), allocatable, dimension(:)       :: latvals_rad, lonvals_rad
  real(r8), allocatable, dimension(:,:)     :: dbuf2
  real(r8), allocatable, dimension(:,:)     :: pstmp
  real(r8), allocatable, dimension(:,:)     :: phis_tmp, var2d
  real(r8), allocatable, dimension(:,:,:)   :: dbuf3, var3d
  real(r8), allocatable, dimension(:,:,:,:) :: dbuf4
  real(r8), pointer, dimension(:,:,:)       :: agrid,grid
  real(r8)                                  :: pertval
  real(r8)                                  :: tracermass(pcnst),delpdry
  real(r8)                                  :: fv3_totwatermass, fv3_airmass
  real(r8)                                  :: initial_global_ave_dry_ps,reldif
  logical                                   :: inic_wet !initial condition is based on wet pressure and water species

  !-----------------------------------------------------------------------

  Atm => dyn_in%Atm
  grid => Atm(mytile)%gridstruct%grid_64
  agrid => Atm(mytile)%gridstruct%agrid_64
  mylindex => dyn_in%mylindex
  mygindex => dyn_in%mygindex

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  isd = Atm(mytile)%bd%isd
  ied = Atm(mytile)%bd%ied
  jsd = Atm(mytile)%bd%jsd
  jed = Atm(mytile)%bd%jed

  fh_topo => topo_file_get_id()
  fh_ini  = initial_file_get_id()


  ! Set mask to indicate which columns are active
  ldof_size=(je-js+1)*(ie-is+1)
  allocate(phis_tmp(ldof_size,1))
  phis_tmp(:,:)=0._r8

  latvals_deg => cam_grid_get_latvals(cam_grid_id('FFSL'))
  lonvals_deg => cam_grid_get_lonvals(cam_grid_id('FFSL'))
  blksize=(ie-is+1)*(je-js+1)

  ! consistency check
  if (blksize /= SIZE(latvals_deg)) then
     call endrun(trim(subname)//': number of latitude values is inconsistent with dynamics block size.')
  end if

  allocate(latvals_rad(blksize))
  allocate(lonvals_rad(blksize))
  latvals_rad(:) = latvals_deg(:)*deg2rad
  lonvals_rad(:) = lonvals_deg(:)*deg2rad

  allocate(glob_ind(blksize))
  do j = js, je
     do i = is, ie
        n=mylindex(i,j)
        glob_ind(n) = mygindex(i,j)
     end do
  end do

  ! Set ICs.  Either from analytic expressions or read from file.

  if (analytic_ic_active()) then
     vcoord = vc_moist_pressure
     inic_wet = .true.
     ! First, initialize all the variables, then assign
     allocate(dbuf2(blksize,1))
     allocate(dbuf3(blksize,nlev,1))
     allocate(dbuf4(blksize,nlev, 1,pcnst))
     dbuf2 = 0.0_r8
     dbuf3 = 0.0_r8
     dbuf4 = 0.0_r8

     allocate(m_ind(pcnst))
     do m_cnst = 1, pcnst
        m_ind(m_cnst) = m_cnst
     end do

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,PS=dbuf2)
     do j = js, je
        do i = is, ie
           ! PS
           n=mylindex(i,j)
           atm(mytile)%ps(i,j) =   dbuf2(n, 1)
        end do
     end do

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind ,            &
          PHIS_OUT=phis_tmp(:,:))

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,            &
          T=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! T
           n=mylindex(i,j)
           atm(mytile)%pt(i,j,:) = dbuf3(n, :, 1)
        end do
     end do


     dbuf3=0._r8
     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,            &
          U=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! U a-grid
           n=mylindex(i,j)
           atm(mytile)%ua(i,j,:) = dbuf3(n, :, 1)
        end do
     end do

     dbuf3=0._r8
     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,            &
          V=dbuf3(:,:,:))

     do j = js, je
        do i = is, ie
           ! V a-grid
           n=mylindex(i,j)
           atm(mytile)%va(i,j,:) = dbuf3(n, :, 1)
        end do
     end do

     call analytic_ic_set_ic(vcoord, latvals_rad, lonvals_rad, glob_ind,            &
          Q=dbuf4(:,:,:,1:pcnst), m_cnst=m_ind)

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

     !-----------------------------------------------------------------------
     call a2d3djt(atm(mytile)%ua, atm(mytile)%va, atm(mytile)%u, atm(mytile)%v, is,  ie,  js,  je, &
                  isd, ied, jsd, jed, npx,npy, nlev, atm(mytile)%gridstruct, atm(mytile)%domain)

     deallocate(dbuf2)
     deallocate(dbuf3)
     deallocate(dbuf4)
     deallocate(m_ind)

  else
     ! Read ICs from file.

     allocate(dbuf3(blksize,nlev,1))
     allocate(var2d(is:ie,js:je))
     allocate(var3d(is:ie,js:je,nlev))

     call pio_seterrorhandling(fh_ini, pio_bcast_error, err_handling)
     ! PSDRY is unambiguous so use that field first if it exists and reset mixing ratios to
     ! wet for FV3. PS (inic_wet) is assumed to be DRY+All wet condensates but could also be
     ! DRY+Q (CAM physics)
     fieldname   = 'PSDRY'
     fieldname2  = 'PS'
     if (dyn_field_exists(fh_ini, trim(fieldname), required=.false.)) then
        inic_wet = .false.
        call read_dyn_var(trim(fieldname), fh_ini, 'ncol', var2d)
     elseif (dyn_field_exists(fh_ini, trim(fieldname2), required=.false.)) then
        inic_wet = .true.
        call read_dyn_var(trim(fieldname2), fh_ini, 'ncol', var2d)
     else
        call endrun(trim(subname)//': PS or PSDRY must be on ncdata')
     end if
     atm(mytile)%ps(is:ie,js:je) = var2d

     ilen = ie-is+1
     jlen = je-js+1

     ! T
     if (dyn_field_exists(fh_ini, 'T')) then
        call read_dyn_var('T', fh_ini, 'ncol_d', var3d)
        atm(mytile)%pt(is:ie,js:je,1:nlev)=var3d(is:ie,js:je,1:nlev)
     else
         call endrun(trim(subname)//': T not found')
     end if

     if (pertlim /= 0.0_r8) then
        if(masterproc) then
           write(iulog,*) trim(subname), ': Adding random perturbation bounded', &
                'by +/- ', pertlim, ' to initial temperature field'
        end if

        call random_seed(size=rndm_seed_sz)
        allocate(rndm_seed(rndm_seed_sz))

        do i=is,ie
           do j=js,je
              indx=mylindex(i,j)
              rndm_seed = glob_ind(indx)
              call random_seed(put=rndm_seed)
              do k=1,nlev
                 call random_number(pertval)
                 pertval = 2.0_r8*pertlim*(0.5_r8 - pertval)
                 atm(mytile)%pt(i,j,k) = atm(mytile)%pt(i,j,k)*(1.0_r8 + pertval)
              end do
           end do
        end do
        deallocate(rndm_seed)
     end if

     ! V
     if (dyn_field_exists(fh_ini, 'V')) then
        call read_dyn_var('V', fh_ini, 'ncol_d', var3d)
        atm(mytile)%va(is:ie,js:je,1:nlev)=var3d(is:ie,js:je,1:nlev)
     else
         call endrun(trim(subname)//': V not found')
     end if

     if (dyn_field_exists(fh_ini, 'U')) then
        call read_dyn_var('U', fh_ini, 'ncol_d', var3d)
        atm(mytile)%ua(is:ie,js:je,1:nlev)   =var3d(is:ie,js:je,1:nlev)
     else
         call endrun(trim(subname)//': U not found')
     end if

     m_cnst=1
     if (dyn_field_exists(fh_ini, 'Q')) then
        call read_dyn_var('Q', fh_ini, 'ncol_d', var3d)
        atm(mytile)%q(is:ie,js:je,1:nlev,m_cnst) = var3d(is:ie,js:je,1:nlev)
     else
         call endrun(trim(subname)//': Q not found')
     end if

     ! Read in or cold-initialize all the tracer fields
     ! Copy tracers defined on unstructured grid onto distributed FFSL grid
     ! Make sure tracers have at least minimum value

     do m_cnst = 2, pcnst
        m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
        found = .false.

        if(cnst_read_iv(m_cnst)) then
           found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)),            &
                required=.false.)
        end if

        if(found) then
           call read_dyn_var(trim(cnst_name(m_cnst)), fh_ini, 'ncol_d', var3d)
           atm(mytile)%q(is:ie,js:je,1:nlev,m_cnst_ffsl) =  var3d(is:ie,js:je,1:nlev)
        else
           dbuf3=0._r8
           if (masterproc) write(iulog,*)'Missing ',trim(cnst_name(m_cnst)),' constituent number', &
                m_cnst,size(latvals_rad),size(dbuf3)
           if (masterproc) write(iulog,*)'Initializing ',trim(cnst_name(m_cnst)),'fv3 constituent number ',&
                m_cnst_ffsl,' to default'
           call cnst_init_default(m_cnst, latvals_rad, lonvals_rad, dbuf3)
           do k=1, nlev
              indx = 1
              do j = js, je
                 do i = is, ie
                    indx=mylindex(i,j)
                    atm(mytile)%q(i,j, k, m_cnst_ffsl) = max(qmin(m_cnst),dbuf3(indx,k,1))
                 end do
              end do
           end do
        end if

     end do ! pcnst

     call a2d3djt(atm(mytile)%ua, atm(mytile)%va, atm(mytile)%u, atm(mytile)%v, is,  ie,  js,  je, &
                  isd, ied, jsd, jed, npx,npy, nlev, atm(mytile)%gridstruct, atm(mytile)%domain)

     ! Put the error handling back the way it was
     call pio_seterrorhandling(fh_ini, err_handling)

     deallocate(dbuf3)
     deallocate(var2d)
     deallocate(var3d)

  end if ! analytic_ic_active

  deallocate(latvals_rad)
  deallocate(lonvals_rad)
  deallocate(glob_ind)

  ! If analytic ICs are being used, we allow constituents in an initial
  ! file to overwrite mixing ratios set by the default constituent initialization
  ! except for the water species.

  call pio_seterrorhandling(fh_ini, pio_bcast_error, err_handling)
  allocate(var3d(is:ie,js:je,nlev))
  do m_cnst = 1, pcnst
     m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)

     if (analytic_ic_active() .and. cnst_is_a_water_species(cnst_name(m_cnst))) cycle

     found = .false.

     if(cnst_read_iv(m_cnst)) then
        found = dyn_field_exists(fh_ini, trim(cnst_name(m_cnst)),            &
             required=.false.)
     end if

     if(found) then
        call read_dyn_var(trim(cnst_name(m_cnst)), fh_ini, 'ncol_d', var3d)
        atm(mytile)%q(is:ie,js:je,1:nlev,m_cnst_ffsl) =  var3d(is:ie,js:je,1:nlev)
     end if
  end do
  deallocate(var3d)
  ! Put the error handling back the way it was
  call pio_seterrorhandling(fh_ini, err_handling)

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
  !
  ! initialize delp (and possibly mixing ratios) from IC fields.
  !
  if (inic_wet) then
     !
     ! /delp/mix ratios/ps consistent with fv3 airmass (dry+all wet tracers) assuming IC is CAM phys airmass (dry+q only)
     !
     allocate(pstmp(isd:ied,jsd:jed))
     pstmp(:,:) = atm(mytile)%ps(:,:)
     atm(mytile)%ps(:,:)=hyai(1)*ps0
     num_wet_species=thermodynamic_active_species_num-dry_air_species_num
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
              fv3_totwatermass=sum(tracermass(thermodynamic_active_species_idx_dycore(1:num_wet_species)))
              fv3_airmass =  delpdry + fv3_totwatermass
              Atm(mytile)%delp(i,j,k) = fv3_airmass
              Atm(mytile)%q(i,j,k,1:pcnst) = tracermass(1:pcnst)/fv3_airmass
              Atm(mytile)%ps(i,j)=Atm(mytile)%ps(i,j)+Atm(mytile)%delp(i, j, k)
           end do
        end do
     end do
     deallocate(pstmp)
  else
     !
     ! Make delp/mix ratios/ps consistent with fv3 airmass (dry+all wet constituents) assuming IC based off dry airmass
     !
     allocate(pstmp(isd:ied,jsd:jed))
     pstmp(:,:) = atm(mytile)%ps(:,:)
     atm(mytile)%ps(:,:)=hyai(1)*ps0
     num_wet_species=thermodynamic_active_species_num-dry_air_species_num
     do k=1,pver
        do j = js, je
           do i = is, ie
              ! this delp is assumed dry.
              delpdry = (((hyai(k+1) - hyai(k))*ps0)         +                &
                   ((hybi(k+1) - hybi(k))*pstmp(i,j)))
              do m=1,pcnst
                 tracermass(m)=delpdry*Atm(mytile)%q(i,j,k,m)
              end do
              fv3_totwatermass=sum(tracermass(thermodynamic_active_species_idx_dycore(1:num_wet_species)))
              fv3_airmass =  delpdry + fv3_totwatermass
              Atm(mytile)%delp(i,j,k) = fv3_airmass
              Atm(mytile)%q(i,j,k,1:pcnst) = tracermass(1:pcnst)/fv3_airmass
              Atm(mytile)%ps(i,j)=Atm(mytile)%ps(i,j)+Atm(mytile)%delp(i, j, k)
              ! check new tracermass
              do m=1,pcnst
                 m_ffsl=qsize_tracer_idx_cam2dyn(m)
                 reldif=(Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl)-tracermass(m_ffsl))/ &
                        tracermass(m_ffsl)
                 if (reldif > abs(1.0e-15_r8)) &
                      write(iulog,*)'mass inconsistency new, old, relative error=',iam,cnst_name(m), &
                      Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_ffsl),tracermass(m_ffsl),reldif
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
  if ( simple_phys ) initial_global_ave_dry_ps = 0                  !do not scale psdry
  call set_dry_mass(Atm, initial_global_ave_dry_ps)


  !$omp parallel do private(i, j)
  do j=js,je
     do i=is,ie
        Atm(mytile)%pe(i,1,j)   = Atm(mytile)%ptop
        Atm(mytile)%pk(i,j,1)   = Atm(mytile)%ptop ** kappa
        Atm(mytile)%peln(i,1,j) = log(Atm(mytile)%ptop )
     enddo
  enddo

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
           Atm(mytile)%pk(i,j,k+1)= Atm(mytile)%pe(i,k+1,j) ** kappa
           Atm(mytile)%peln(i,k+1,j) = log(Atm(mytile)%pe(i,k+1,j))
           Atm(mytile)%pkz(i,j,k) = (Atm(mytile)%pk(i,j,k+1)-Atm(mytile)%pk(i,j,k)) / &
           (kappa*(Atm(mytile)%peln(i,k+1,j)-Atm(mytile)%peln(i,k,j)))
        enddo
     enddo
  enddo
!!  Initialize non hydrostatic variables if needed
  if (.not. Atm(mytile)%flagstruct%hydrostatic) then
     do k=1,nlev
        do j=js,je
           do i=is,ie
              Atm(mytile)%w ( i,j,k ) = 0._r8
              Atm(mytile)%delz ( i,j,k ) = -rdgas/gravit*Atm(mytile)%pt( i,j,k ) * &
                   ( Atm(mytile)%peln( i,k+1,j ) - Atm(mytile)%peln( i,k,j ) )
           enddo
        enddo
     enddo
  end if

  ! once we've read or initialized all the fields we call update_domains to
  ! update the halo regions

  call mpp_update_domains( Atm(mytile)%phis, Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%ps,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%u,atm(mytile)%v,Atm(mytile)%domain,gridtype=DGRID_NE,complete=.true. )
  call mpp_update_domains( atm(mytile)%pt,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%delp,   Atm(mytile)%domain )
  call mpp_update_domains( atm(mytile)%q,    Atm(mytile)%domain )

  ! Cleanup
  deallocate(phis_tmp)

end subroutine read_inidat

!=======================================================================

  subroutine calc_tot_energy_dynamics(atm,suffix)
    use physconst,              only: gravit, cpair, rearth,omega
    use physconst,              only: thermodynamic_active_species_num,thermodynamic_active_species_idx_dycore
    use physconst,              only: thermodynamic_active_species_cp,thermodynamic_active_species_cv,dry_air_species_num
    use cam_history,            only: outfld, hist_fld_active
    use constituents,           only: cnst_get_ind
    use dimensions_mod,         only: nlev
    use fv_mp_mod,              only: ng
    !------------------------------Arguments--------------------------------

    type(fv_atmos_type), pointer, intent(in) :: Atm(:)
    character(len=*)    , intent(in) :: suffix ! suffix for "outfld" names

    !---------------------------Local storage-------------------------------

    real(kind=r8), allocatable, dimension(:,:) :: se,              &! Dry Static energy (J/m2)
                                                  ke,              &! kinetic energy    (J/m2)
                                                  ps_local          ! ps temp based on CAM or FV3 airmass
    real(kind=r8), allocatable, dimension(:,:) :: wv,wl,wi,wr,ws,wg ! col integ constiuents(kg/m2)
    real(kind=r8), allocatable, dimension(:,:) :: tt                ! column integrated test tracer (kg/m2)
    real(kind=r8), allocatable, dimension(:,:,:) :: dp,delpograv
    real(kind=r8) :: se_tmp, dpdry
    real(kind=r8) :: ke_tmp
    real(kind=r8) :: wv_tmp,wl_tmp,wi_tmp,wr_tmp,ws_tmp,wg_tmp
    real(kind=r8) :: tt_tmp

    !
    ! global axial angular momentum (AAM) can be separated into one part (mr)
    ! associated with the relative motion of the atmosphere with respect to the planet surface
    ! (also known as wind AAM) and another part (mo) associated with the angular velocity OMEGA
    ! (2*pi/d, where d is the length of the day) of the planet (also known as mass AAM)
    !
    real(kind=r8), allocatable, dimension(:,:) :: mr  ! wind AAM
    real(kind=r8), allocatable, dimension(:,:) :: mo  ! mass AAM
    real(kind=r8) :: mr_cnst, mo_cnst, cos_lat, mr_tmp, mo_tmp

    real(kind=r8) :: se_glob, ke_glob, wv_glob, wl_glob, wi_glob, &
                     wr_glob, ws_glob, wg_glob, tt_glob, mr_glob, mo_glob

    integer :: i,j,k,nq,idim,m_cnst_ffsl
    integer :: ixcldice, ixcldliq, ixtt,ixcldliq_ffsl,ixcldice_ffsl ! CLDICE, CLDLIQ and test tracer indices
    integer :: ixrain, ixsnow, ixgraupel,ixrain_ffsl, ixsnow_ffsl, ixgraupel_ffsl
    character(len=16) :: se_name,ke_name,wv_name,wl_name, &
                         wi_name,wr_name,ws_name,wg_name,tt_name,mo_name,mr_name

    integer :: is,ie,js,je,isd,ied,jsd,jed
    logical :: printglobals = .false.
    !-----------------------------------------------------------------------

    is = Atm(mytile)%bd%is
    ie = Atm(mytile)%bd%ie
    js = Atm(mytile)%bd%js
    je = Atm(mytile)%bd%je
    isd = Atm(mytile)%bd%isd
    ied = Atm(mytile)%bd%ied
    jsd = Atm(mytile)%bd%jsd
    jed = Atm(mytile)%bd%jed

    se_glob = 0._r8
    ke_glob = 0._r8
    wv_glob = 0._r8
    wl_glob = 0._r8
    wi_glob = 0._r8
    wr_glob = 0._r8
    ws_glob = 0._r8
    wg_glob = 0._r8
    tt_glob = 0._r8
    mr_glob = 0._r8
    mo_glob = 0._r8

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
    allocate(dp(is:ie,js:je,nlev))
    allocate(delpograv(is:ie,js:je,nlev))
    allocate(ps_local(is:ie,js:je))

    se_name = 'SE_'   //trim(suffix)
    ke_name = 'KE_'   //trim(suffix)
    wv_name = 'WV_'   //trim(suffix)
    wl_name = 'WL_'   //trim(suffix)
    wi_name = 'WI_'   //trim(suffix)
    wr_name = 'WR_'   //trim(suffix)
    ws_name = 'WS_'   //trim(suffix)
    wg_name = 'WG_'   //trim(suffix)
    tt_name = 'TT_'   //trim(suffix)


    if ( hist_fld_active(se_name).or.hist_fld_active(ke_name).or. &
         hist_fld_active(wv_name).or.hist_fld_active(wl_name).or. &
         hist_fld_active(wi_name).or.hist_fld_active(wr_name).or. &
         hist_fld_active(ws_name).or.hist_fld_active(wg_name).or. &
         hist_fld_active(tt_name)) then
       if (thermodynamic_active_species_num-dry_air_species_num > 1) then
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

       delpograv(is:ie,js:je,1:nlev) = Atm(mytile)%delp(is:ie,js:je,1:nlev)/gravit ! temporary

       !
       ! Calculate Energy, CAM or FV3 based on fv3_lcp_moist and fv3_lcv_moist
       !


       do k = 1, nlev
          do j=js,je
             do i = is, ie
                ! initialize dp with delp
                dp(i,j,k) = Atm(mytile)%delp(i,j,k)
                !
                ! if neither fv3_lcp_moist and fv3_lcv_moist is set then
                ! use cam definition of internal energy
                ! adjust dp to be consistent with CAM physics air mass (only water vapor and dry air in pressure)
                if ((.not.fv3_lcp_moist).and.(.not.fv3_lcv_moist)) then
                   if (thermodynamic_active_species_num-dry_air_species_num > 1) then
                      ! adjust dp to include just dry + vap to use below
                      do nq=2,thermodynamic_active_species_num-dry_air_species_num
                         m_cnst_ffsl=thermodynamic_active_species_idx_dycore(nq)
                         dp(i,j,k) = dp(i,j,k) - &
                              Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,m_cnst_ffsl)
                      end do
                   end if
                   se_tmp  = cpair*Atm(mytile)%pt(i,j,k)*dp(i,j,k)/gravit
                else
                   ! if either fv3_lcp_moist or fv3_lcv_moist is set then
                   ! use all condensates in calculation of energy and dp
                   ! Start with energy of dry air and add energy of condensates
                   dpdry = Atm(mytile)%delp(i,j,k)
                   do nq=1,thermodynamic_active_species_num-dry_air_species_num
                      m_cnst_ffsl=thermodynamic_active_species_idx_dycore(nq)
                      dpdry = dpdry - Atm(mytile)%delp(i,j,k)*Atm(mytile)%q(i,j,k,nq)
                   end do
                   se_tmp = cpair*dpdry
                   do nq=1,thermodynamic_active_species_num-dry_air_species_num
                      m_cnst_ffsl=thermodynamic_active_species_idx_dycore(nq)
                      if (fv3_lcp_moist) then
                         se_tmp = se_tmp + &
                              thermodynamic_active_species_cp(nq)*Atm(mytile)%q(i,j,k,m_cnst_ffsl) * &
                              Atm(mytile)%delp(i,j,k)
                      end if
                      if (fv3_lcv_moist) then
                         se_tmp = se_tmp + &
                              thermodynamic_active_species_cv(nq)*Atm(mytile)%q(i,j,k,m_cnst_ffsl) * &
                              Atm(mytile)%delp(i,j,k)
                      end if
                   end do
                   se_tmp = se_tmp*Atm(mytile)%pt(i,j,k)/gravit
                end if
                ke_tmp   = 0.5_r8*(Atm(mytile)%va(i,j,k)**2+ Atm(mytile)%ua(i,j,k)**2)*dp(i,j,k)/gravit
                wv_tmp   =  Atm(mytile)%q(i,j,k,1)*delpograv(i,j,k)

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
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   wl_tmp   = Atm(mytile)%q(i,j,k,ixcldliq_ffsl)*delpograv(i,j,k)
                   wl   (i,j) = wl(i,j) + wl_tmp
                end do
             end do
          end do
       end if

       if (ixcldice > 1) then
          ixcldice_ffsl = qsize_tracer_idx_cam2dyn(ixcldice)
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   wi_tmp   = Atm(mytile)%q(i,j,k,ixcldice_ffsl)*delpograv(i,j,k)
                   wi(i,j)    = wi(i,j) + wi_tmp
                end do
             end do
          end do
       end if

       if (ixrain > 1) then
          ixrain_ffsl = qsize_tracer_idx_cam2dyn(ixrain)
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   wr_tmp   = Atm(mytile)%q(i,j,k,ixrain_ffsl)*delpograv(i,j,k)
                   wr   (i,j) = wr(i,j) + wr_tmp
                end do
             end do
          end do
       end if

       if (ixsnow > 1) then
          ixsnow_ffsl = qsize_tracer_idx_cam2dyn(ixsnow)
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   ws_tmp   = Atm(mytile)%q(i,j,k,ixsnow_ffsl)*delpograv(i,j,k)
                   ws(i,j)    = ws(i,j) + ws_tmp
                end do
             end do
          end do
       end if

       if (ixgraupel > 1) then
          ixgraupel_ffsl = qsize_tracer_idx_cam2dyn(ixgraupel)
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   wg_tmp   = Atm(mytile)%q(i,j,k,ixgraupel_ffsl)*delpograv(i,j,k)
                   wg(i,j)    = wg(i,j) + wg_tmp
                end do
             end do
          end do
       end if


       if (ixtt > 1) then
          do k = 1, nlev
             do j = js, je
                do i = is, ie
                   tt_tmp   = Atm(mytile)%q(i,j,k,ixtt)*delpograv(i,j,k)
                   tt   (i,j) = tt(i,j) + tt_tmp
                end do
             end do
          end do
       end if
       idim=ie-is+1
       do j=js,je
          ! Output energy diagnostics
          call outfld(se_name  ,se(:,j)       ,idim, j)
          call outfld(ke_name  ,ke(:,j)       ,idim, j)
          call outfld(wv_name  ,wv(:,j)       ,idim, j)
          call outfld(wl_name  ,wl(:,j)       ,idim, j)
          call outfld(wi_name  ,wi(:,j)       ,idim, j)
          call outfld(wr_name  ,wr(:,j)       ,idim, j)
          call outfld(ws_name  ,ws(:,j)       ,idim, j)
          call outfld(wg_name  ,wg(:,j)       ,idim, j)
          if (ixtt > 1) call outfld(tt_name  ,tt(:,j)       ,idim, j)
       end do

       if (printglobals) then
          se_glob=g_sum(Atm(mytile)%domain, se(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          ke_glob=g_sum(Atm(mytile)%domain, ke(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wv_glob=g_sum(Atm(mytile)%domain, wv(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wl_glob=g_sum(Atm(mytile)%domain, wl(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wi_glob=g_sum(Atm(mytile)%domain, wi(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wr_glob=g_sum(Atm(mytile)%domain, wr(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          ws_glob=g_sum(Atm(mytile)%domain, ws(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          wg_glob=g_sum(Atm(mytile)%domain, wg(is:ie,js:je), is, ie, js, je, &
                  Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          if (ixtt > 1) &
               tt_glob=g_sum(Atm(mytile)%domain, tt(is:ie,js:je), is, ie, js, je, &
                       Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
          if (masterproc) then

             write(iulog, '(a,e25.17)') 'static energy se_'//trim(suffix)//')        = ',se_glob
             write(iulog, '(a,e25.17)') 'kinetic energy ke_'//trim(suffix)//')           = ',ke_glob
             write(iulog, '(a,e25.17)') 'total energy se_plus_ke_'//trim(suffix)//')     = ',(ke_glob+se_glob)
             write(iulog, '(a,e25.17)') 'integrated vapor wv_'//trim(suffix)//'   = ',wv_glob
             write(iulog, '(a,e25.17)') 'integrated liquid wl_'//trim(suffix)//'  = ',wl_glob
             write(iulog, '(a,e25.17)') 'integrated ice wi_'//trim(suffix)//'     = ',wi_glob
             write(iulog, '(a,e25.17)') 'integrated liquid rain wr_'//trim(suffix)//'   = ',wr_glob
             write(iulog, '(a,e25.17)') 'integrated liquid snow ws_'//trim(suffix)//'  = ',ws_glob
             write(iulog, '(a,e25.17)') 'integrated graupel wg_'//trim(suffix)//'     = ',wg_glob
             if (ixtt > 1) write(iulog, '(a,e25.17)') &
                  'global column integrated test tracer tt_'//trim(suffix)//' = ',tt_glob
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
    mr_name = 'MR_'   //trim(suffix)
    mo_name = 'MO_'   //trim(suffix)

    if ( hist_fld_active(mr_name).or.hist_fld_active(mo_name)) then



      mr_cnst = rearth**3/gravit
      mo_cnst = omega*rearth**4/gravit
      mr    = 0.0_r8
      mo    = 0.0_r8
      do k = 1, nlev
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
         call outfld(mr_name  ,mr(is:ie,j)       ,idim,j)
         call outfld(mo_name  ,mo(is:ie,j)       ,idim,j)
      end do

      if (printglobals) then
         mr_glob=g_sum(Atm(mytile)%domain, mr(is:ie,js:je), is, ie, js, je, &
                 Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
         mo_glob=g_sum(Atm(mytile)%domain, mo(is:ie,js:je), is, ie, js, je, &
                 Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
         if (masterproc) then
            write(iulog, '(a,e25.17)') 'integrated wind AAM '//trim(mr_name)//' = ',mr_glob
            write(iulog, '(a,e25.17)') 'integrated mass AAM '//trim(mo_name)//' = ',mo_glob
         end if
      end if
   end if

  deallocate(ps_local)
  deallocate(dp)
  deallocate(delpograv)
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

   ! Arguments
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
    call infld(fieldname, fh,dimname, 'lev', 1, ldof_size, 1, pver,     &
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

  ! Arguments
  type(fv_grid_bounds_type), intent(in) :: bd
  real(r8), intent(in)          :: field(bd%is:bd%ie,bd%js:bd%je)
  character(len=*)    , intent(in) :: outfld_name ! suffix for "outfld" names

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

  use constituents,          only: pcnst, qmin
  use cam_logfile,           only: iulog
  use hycoef,                only: hyai, hybi, ps0
  use dimensions_mod,        only: nlev
  use dyn_grid,              only: mytile
  use physconst,             only: thermodynamic_active_species_num,thermodynamic_active_species_idx_dycore,dry_air_species_num

  ! Arguments
  type (fv_atmos_type), intent(in),  pointer :: Atm(:)
  real (kind=r8), intent(in)                 :: fixed_global_ave_dry_ps

  ! local
  real (kind=r8)                             :: global_ave_ps_inic,global_ave_dryps_inic,global_ave_dryps_scaled, &
                                                global_ave_ps_new,global_ave_dryps_new
  real (r8), allocatable, dimension(:,:)     :: psdry, psdry_scaled, psdry_new
  real (r8), allocatable, dimension(:,:,:)   :: factor, delpwet, delpdry, newdelp
  integer                                    :: i, j ,k, m,is,ie,js,je
  integer                                    :: num_wet_species   ! first tracers in FV3 tracer array

  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je
  allocate(factor(is:ie,js:je,nlev))
  allocate(delpdry(is:ie,js:je,nlev))
  allocate(delpwet(is:ie,js:je,nlev))
  allocate(newdelp(is:ie,js:je,nlev))
  allocate(psdry(is:ie,js:je))
  allocate(psdry_scaled(is:ie,js:je))
  allocate(psdry_new(is:ie,js:je))


  if (fixed_global_ave_dry_ps == 0) return;

  ! get_global_ave_surface_pressure - must use bitwise sum (reproducable)
  global_ave_ps_inic=g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, &
                           Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  num_wet_species=thermodynamic_active_species_num-dry_air_species_num
  do k=1,pver
     do j = js, je
        do i = is, ie
           delpdry(i,j,k)=Atm(mytile)%delp(i,j,k) * (1.0_r8 - &
                          sum(Atm(mytile)%q(i,j,k,thermodynamic_active_species_idx_dycore(1:num_wet_species))))
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

  global_ave_dryps_inic=g_sum(Atm(mytile)%domain, psdry(is:ie,js:je), is, ie, js, je, &
                        Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  psdry_scaled = psdry*(fixed_global_ave_dry_ps/global_ave_dryps_inic)

  global_ave_dryps_scaled=g_sum(Atm(mytile)%domain, psdry_scaled(is:ie,js:je), is, ie, js, je, &
                          Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

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
  global_ave_ps_new=   g_sum(Atm(mytile)%domain, Atm(mytile)%ps(is:ie,js:je), is, ie, js, je, &
                             Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)
  global_ave_dryps_new=g_sum(Atm(mytile)%domain, psdry_new(is:ie,js:je), is, ie, js, je, &
                             Atm(mytile)%ng, Atm(mytile)%gridstruct%area_64, 1, .true.)

  if (masterproc) then
     write (iulog,*) "-------------------------- set_dry_mass---------------------------------------------"
     write (iulog,*) "Scaling dry surface pressure to global average of = ",&
          fixed_global_ave_dry_ps/100.0_r8,"hPa"
     write (iulog,*) "Average surface pressure in initial condition = ", &
          global_ave_ps_inic/100.0_r8,"hPa"
     write (iulog,*) "Average dry surface pressure in initial condition = ",&
          global_ave_dryps_inic/100.0_r8,"hPa"
     write (iulog,*) "Average surface pressure after scaling = ",global_ave_ps_new/100.0_r8,"hPa"
     write (iulog,*) "Average dry surface pressure after scaling = ",global_ave_dryps_new/100.0_r8,"hPa"
     write (iulog,*) "Change in surface pressure           = ",&
          global_ave_ps_new-global_ave_ps_inic,"Pa"
     write (iulog,*) "Change in dry surface pressure           = ",&
          global_ave_dryps_new-global_ave_dryps_inic,"Pa"
     write (iulog,*) "Mixing ratios have been scaled so that total mass of tracer is conserved"
     write (iulog,*) "Total precipitable water before scaling = ", &
          (global_ave_ps_inic-global_ave_dryps_inic)/gravit, '(kg/m**2)'
     write (iulog,*) "Total precipitable water after  scaling = ", &
          (global_ave_ps_new-global_ave_dryps_new)/gravit, '(kg/m**2)'
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

subroutine a2d3djt(ua, va, u, v, is,  ie,  js,  je, isd, ied, jsd, jed, npx,npy, nlev, gridstruct, domain)

! This routine interpolates cell centered a-grid winds to d-grid (cell edges)

  use mpp_domains_mod,    only: mpp_update_domains,  DGRID_NE
  use fv_arrays_mod,      only: fv_grid_type

  ! arguments
  integer, intent(in)                                          :: is,  ie,  js,  je
  integer, intent(in)                                          :: isd, ied, jsd, jed
  integer, intent(in)                                          :: npx,npy, nlev
  real(r8), intent(inout), dimension(isd:ied,  jsd:jed+1,nlev) :: u
  real(r8), intent(inout), dimension(isd:ied+1,jsd:jed  ,nlev) :: v
  real(r8), intent(inout), dimension(isd:ied,jsd:jed,nlev)     :: ua, va
  type(fv_grid_type), intent(in), target                       :: gridstruct
  type(domain2d), intent(inout)                                :: domain

  ! local:
  real(r8), dimension(is-1:ie+1,js-1:je+1,3) :: v3
  real(r8), dimension(is-1:ie+1,js:je+1,3)   :: ue    ! 3D winds at edges
  real(r8), dimension(is:ie+1,js-1:je+1,  3) :: ve    ! 3D winds at edges
  real(r8), dimension(is:ie)                 :: ut1, ut2, ut3
  real(r8), dimension(js:je)                 :: vt1, vt2, vt3
  integer                                    :: i, j, k, im2, jm2

  real(r8), pointer, dimension(:,:,:)        :: vlon, vlat
  real(r8), pointer, dimension(:,:,:,:)      :: es, ew
  real(r8), pointer, dimension(:)            :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n

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

!$OMP parallel do default(none) shared(is,ie,js,je,nlev,gridstruct,u,ua,v,va,  &
!$OMP                                  vlon,vlat,jm2,edge_vect_w,npx,edge_vect_e,im2, &
!$OMP                                  edge_vect_s,npy,edge_vect_n,es,ew)             &
!$OMP                          private(i,j,k,ut1, ut2, ut3, vt1, vt2, vt3, ue, ve, v3)
  do k=1, nlev

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
                 vt1(j) = edge_vect_w(j)*ve(i,j-1,1)+(1._r8-edge_vect_w(j))*ve(i,j,1)
                 vt2(j) = edge_vect_w(j)*ve(i,j-1,2)+(1._r8-edge_vect_w(j))*ve(i,j,2)
                 vt3(j) = edge_vect_w(j)*ve(i,j-1,3)+(1._r8-edge_vect_w(j))*ve(i,j,3)
              else
                 vt1(j) = edge_vect_w(j)*ve(i,j+1,1)+(1._r8-edge_vect_w(j))*ve(i,j,1)
                 vt2(j) = edge_vect_w(j)*ve(i,j+1,2)+(1._r8-edge_vect_w(j))*ve(i,j,2)
                 vt3(j) = edge_vect_w(j)*ve(i,j+1,3)+(1._r8-edge_vect_w(j))*ve(i,j,3)
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
                 vt1(j) = edge_vect_e(j)*ve(i,j-1,1)+(1._r8-edge_vect_e(j))*ve(i,j,1)
                 vt2(j) = edge_vect_e(j)*ve(i,j-1,2)+(1._r8-edge_vect_e(j))*ve(i,j,2)
                 vt3(j) = edge_vect_e(j)*ve(i,j-1,3)+(1._r8-edge_vect_e(j))*ve(i,j,3)
              else
                 vt1(j) = edge_vect_e(j)*ve(i,j+1,1)+(1._r8-edge_vect_e(j))*ve(i,j,1)
                 vt2(j) = edge_vect_e(j)*ve(i,j+1,2)+(1._r8-edge_vect_e(j))*ve(i,j,2)
                 vt3(j) = edge_vect_e(j)*ve(i,j+1,3)+(1._r8-edge_vect_e(j))*ve(i,j,3)
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
                 ut1(i) = edge_vect_s(i)*ue(i-1,j,1)+(1._r8-edge_vect_s(i))*ue(i,j,1)
                 ut2(i) = edge_vect_s(i)*ue(i-1,j,2)+(1._r8-edge_vect_s(i))*ue(i,j,2)
                 ut3(i) = edge_vect_s(i)*ue(i-1,j,3)+(1._r8-edge_vect_s(i))*ue(i,j,3)
              else
                 ut1(i) = edge_vect_s(i)*ue(i+1,j,1)+(1._r8-edge_vect_s(i))*ue(i,j,1)
                 ut2(i) = edge_vect_s(i)*ue(i+1,j,2)+(1._r8-edge_vect_s(i))*ue(i,j,2)
                 ut3(i) = edge_vect_s(i)*ue(i+1,j,3)+(1._r8-edge_vect_s(i))*ue(i,j,3)
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
                 ut1(i) = edge_vect_n(i)*ue(i-1,j,1)+(1._r8-edge_vect_n(i))*ue(i,j,1)
                 ut2(i) = edge_vect_n(i)*ue(i-1,j,2)+(1._r8-edge_vect_n(i))*ue(i,j,2)
                 ut3(i) = edge_vect_n(i)*ue(i-1,j,3)+(1._r8-edge_vect_n(i))*ue(i,j,3)
              else
                 ut1(i) = edge_vect_n(i)*ue(i+1,j,1)+(1._r8-edge_vect_n(i))*ue(i,j,1)
                 ut2(i) = edge_vect_n(i)*ue(i+1,j,2)+(1._r8-edge_vect_n(i))*ue(i,j,2)
                 ut3(i) = edge_vect_n(i)*ue(i+1,j,3)+(1._r8-edge_vect_n(i))*ue(i,j,3)
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
