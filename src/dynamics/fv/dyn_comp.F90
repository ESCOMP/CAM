module dyn_comp

!----------------------------------------------------------------------
!   This module contains the interfaces for the Finite-Volume
!   Dynamical Core used in the Community Atmospheric Model
!   (FVCAM). This component will hereafter be referred
!   to as the ``FVdycore'' gridded component.  FVdycore
!   consists of four sub-components,
!
!      cd_core:  The C/D-grid dycore component
!      te_map:   Vertical remapping algorithm
!      trac2d:   Tracer advection
!      benergy:  Energy balance
!
!  FVdycore maintains an internal state consisting of the
!  following fields:  control variables
!
!     U:    U winds on a D-grid (m/s)
!     V:    V winds on a D-grid (m/s)
!     PT:   Scaled Virtual Potential Temperature (T_v/PKZ)
!     PE:   Edge pressures
!     Q:    Tracers
!     PKZ:  Consistent mean for p^kappa
!
!  as well as a GRID and same additional run-specific variables
!  (dt, iord, jord, nsplit, nspltrac, nspltvrm)
!
! Note: PT is not updated if the flag CONVT is true.
!
! The internal state is updated each time FVdycore is called.
!
! REVISION HISTORY:
!
!   WS  05.06.10:  Adapted from FVdycore_GridCompMod
!   WS  05.09.20:  Renamed dyn_comp
!   WS  05.11.10:  Now using dyn_import/export_t containers
!   WS  06.03.01:  Removed tracertrans-related variables
!   WS  06.04.13:  dyn_state moved here from prognostics
!   CC  07.01.29:  Corrected calculation of OMGA
!   AM  07.10.31:  Supports overlap of trac2d and cd_core subcycles
!----------------------------------------------------------------------

use shr_kind_mod,       only: r8=>shr_kind_r8
use spmd_utils,         only: masterproc, iam
use physconst,          only: pi

use pmgrid,             only: plon, plat
use constituents,       only: pcnst, cnst_name, cnst_read_iv, qmin, cnst_type, &
                              cnst_is_a_water_species

use time_manager,       only: get_step_size

use dynamics_vars,      only: t_fvdycore_grid,            &
                              t_fvdycore_state, t_fvdycore_constants
use dyn_internal_state, only: get_dyn_state, get_dyn_state_grid

use dyn_grid,           only: get_horiz_grid_dim_d
use commap,             only: clat, clon, clat_staggered, londeg_st
use spmd_dyn,           only: spmd_readnl

use inic_analytic,      only: analytic_ic_active, analytic_ic_set_ic
use dyn_tests_utils,    only: vc_moist_pressure

use cam_control_mod,    only: initial_run, moist_physics
use phys_control,       only: phys_setopts

use cam_initfiles,      only: initial_file_get_id, topo_file_get_id, pertlim
use cam_pio_utils,      only: clean_iodesc_list
use ncdio_atm,          only: infld
use pio,                only: pio_seterrorhandling, pio_bcast_error, pio_noerr, &
                              file_desc_t, pio_inq_dimid, pio_inq_dimlen,       &
                              pio_inq_varid, pio_get_att

use perf_mod,           only: t_startf, t_stopf, t_barrierf
use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun

use par_vecsum_mod,     only: par_vecsum
use te_map_mod,         only: te_map

implicit none
private
save

public :: &
   dyn_readnl,   &
   dyn_register, &
   dyn_init,     &
   dyn_run,      &
   dyn_final,    &
   dyn_import_t, &
   dyn_export_t, &
   dyn_state,    &
   frontgf_idx,  &
   frontga_idx,  &
   uzm_idx,      &
   initial_mr

type (t_fvdycore_state), target :: dyn_state

type dyn_import_t
     real(r8), dimension(:,: ),    pointer :: phis   ! Surface geopotential
     real(r8), dimension(:,: ),    pointer :: ps     ! Surface pressure
     real(r8), dimension(:,:,:),   pointer :: u3s    ! U-winds (staggered)
     real(r8), dimension(:,:,:),   pointer :: v3s    ! V-winds (staggered)
     real(r8), dimension(:,:,:),   pointer :: pe     ! Pressure
     real(r8), dimension(:,:,:),   pointer :: pt     ! Potential temperature
     real(r8), dimension(:,:,:),   pointer :: t3     ! Temperatures
     real(r8), dimension(:,:,:),   pointer :: pk     ! Pressure to the kappa
     real(r8), dimension(:,:,:),   pointer :: pkz    ! Pressure to the kappa offset
     real(r8), dimension(:,:,:),   pointer :: delp   ! Delta pressure
     real(r8), dimension(:,:,:,:), pointer :: tracer ! Tracers
end type dyn_import_t

type dyn_export_t
     real(r8), dimension(:,: ),    pointer :: phis   ! Surface geopotential
     real(r8), dimension(:,: ),    pointer :: ps     ! Surface pressure
     real(r8), dimension(:,:,:),   pointer :: u3s    ! U-winds (staggered)
     real(r8), dimension(:,:,:),   pointer :: v3s    ! V-winds (staggered)
     real(r8), dimension(:,:,:),   pointer :: pe     ! Pressure
     real(r8), dimension(:,:,:),   pointer :: pt     ! Potential temperature
     real(r8), dimension(:,:,:),   pointer :: t3     ! Temperatures
     real(r8), dimension(:,:,:),   pointer :: pk     ! Pressure to the kappa
     real(r8), dimension(:,:,:),   pointer :: pkz    ! Pressure to the kappa offset
     real(r8), dimension(:,:,:),   pointer :: delp   ! Delta pressure
     real(r8), dimension(:,:,:,:), pointer :: tracer ! Tracers
     real(r8), dimension(:,:,:),   pointer :: peln   !
     real(r8), dimension(:,:,:),   pointer :: omga   ! Vertical velocity
     real(r8), dimension(:,:,:),   pointer :: mfx    ! Mass flux in X
     real(r8), dimension(:,:,:),   pointer :: mfy    ! Mass flux in Y
     real(r8), dimension(:,:,:),   pointer :: du3s   ! U-wind tend. from dycore (staggered)
     real(r8), dimension(:,:,:),   pointer :: dv3s   ! V-wind tend. from dycore (staggered)
     real(r8), dimension(:,:,:),   pointer :: dua3s  ! U-wind tend. from advection (stagg)
     real(r8), dimension(:,:,:),   pointer :: dva3s  ! V-wind tend. from advection (stagg)
     real(r8), dimension(:,:,:),   pointer :: duf3s  ! U-wind tend. from fixer (staggered)
end type dyn_export_t

! The FV core is always called in its "full physics" mode.  We don't want
! the dycore to know what physics package is responsible for the forcing.
logical, parameter         :: convt = .true.

! Indices for fields that are computed in the dynamics and passed to the physics
! via the physics buffer
integer, protected :: frontgf_idx  = -1
integer, protected :: frontga_idx  = -1
integer, protected :: uzm_idx = -1

character(len=3), protected :: initial_mr(pcnst) ! constituents initialized with wet or dry mr

logical :: readvar            ! inquiry flag:  true => variable exists on netCDF file

character(len=8)  :: fv_print_dpcoup_warn = "off"
public            :: fv_print_dpcoup_warn

!=============================================================================================
CONTAINS
!=============================================================================================

subroutine dyn_readnl(nlfilename)

   ! Read dynamics namelist group.
   use units,           only: getunit, freeunit
   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8, &
                              mpi_logical, mpi_character
   use ctem,            only: ctem_readnl
   use fill_module,     only: fill_readnl

   ! args
   character(len=*), intent(in) :: nlfilename

   ! Local variables
   integer :: ierr
   integer :: unitn
   character(len=*), parameter ::  sub="dyn_readnl"

   integer :: fv_nsplit      = 0       ! Lagrangian time splits
   integer :: fv_nspltrac    = 0       ! Tracer time splits
   integer :: fv_nspltvrm    = 0       ! Vertical re-mapping time splits
   integer :: fv_iord        = 4       ! scheme to be used in E-W direction
   integer :: fv_jord        = 4       ! scheme to be used in N-S direction
   integer :: fv_kord        = 4       ! scheme to be used for vertical mapping
                                       ! _ord = 1: first order upwind
                                       ! _ord = 2: 2nd order van Leer (Lin et al 1994)
                                       ! _ord = 3: standard PPM
                                       ! _ord = 4: enhanced PPM (default)
   logical :: fv_conserve    = .false. ! Flag indicating whether the dynamics is conservative
   integer :: fv_filtcw      = 0       ! flag for filtering c-grid winds
   integer :: fv_fft_flt     = 1       ! 0 => FFT/algebraic filter; 1 => FFT filter
   integer :: fv_div24del2flag = 2     ! 2 for 2nd order div damping, 4 for 4th order div damping,
                                       ! 42 for 4th order div damping plus 2nd order velocity damping
   real(r8):: fv_del2coef    = 3.e5_r8 ! strength of 2nd order velocity damping
   logical :: fv_high_altitude = .false. ! switch to apply variables appropriate for high-altitude physics

   logical :: fv_high_order_top=.false.! do not degrade calculation to 1st order near the model top

   logical :: fv_am_correction=.false. ! apply correction for angular momentum (AM) in SW eqns
   logical :: fv_am_geom_crrct=.false. ! apply correction for angular momentum (AM) in geometry
   logical :: fv_am_fixer     =.false. ! apply global fixer to conserve AM
   logical :: fv_am_fix_lbl   =.false. ! apply global AM fixer level by level
   logical :: fv_am_diag      =.false. ! turns on an AM diagnostic calculation written to log file

   namelist /dyn_fv_inparm/ fv_nsplit, fv_nspltrac, fv_nspltvrm, fv_iord, fv_jord, &
                            fv_kord, fv_conserve, fv_filtcw, fv_fft_flt,           &
                            fv_div24del2flag, fv_del2coef, fv_high_order_top,      &
                            fv_am_correction, fv_am_geom_crrct, &
                            fv_am_fixer, fv_am_fix_lbl, fv_am_diag, fv_high_altitude, &
                            fv_print_dpcoup_warn

   type(t_fvdycore_state), pointer :: dyn_state

   real(r8) :: dt
   !-----------------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'Read in dyn_fv_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for dyn_fv_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'dyn_fv_inparm', status=ierr)
      if (ierr == 0) then  ! found dyn_fv_inparm
         read(unitn, dyn_fv_inparm, iostat=ierr)  ! read the dyn_fv_inparm namelist group
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading dyn_fv_inparm')
         end if
      else
         call endrun(sub//': can''t find dyn_fv_inparm in file '//trim(nlfilename))
      end if
      close( unitn )
      call freeunit( unitn )
   endif

   call mpi_bcast(fv_nsplit, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_nsplit")

   call mpi_bcast(fv_nspltrac, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_nspltrac")

   call mpi_bcast(fv_nspltvrm, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_nspltvrm")

   call mpi_bcast(fv_iord, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_iord")

   call mpi_bcast(fv_jord, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_jord")

   call mpi_bcast(fv_kord, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_kord")

   call mpi_bcast(fv_conserve, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_conserve")

   call mpi_bcast(fv_filtcw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_filtcw")

   call mpi_bcast(fv_fft_flt, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_fft_flt")

   call mpi_bcast(fv_div24del2flag, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_div24del2flag")

   call mpi_bcast(fv_del2coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_del2coef")

   call mpi_bcast(fv_high_order_top, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_high_order_top")

   call mpi_bcast(fv_am_correction, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_am_correction")

   call mpi_bcast(fv_am_geom_crrct, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_am_geom_crrct")

   ! if fv_am_fix_lbl is true then fv_am_fixer must also be true.
   if (fv_am_fix_lbl .and. .not. fv_am_fixer) then
      fv_am_fixer = .true.
   end if

   call mpi_bcast(fv_am_fixer, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_am_fixer")

   call mpi_bcast(fv_am_fix_lbl, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_am_fix_lbl")

   call mpi_bcast(fv_am_diag, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_am_diag")

   call mpi_bcast(fv_high_altitude, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_high_altitude")

   call mpi_bcast(fv_print_dpcoup_warn, len(fv_print_dpcoup_warn), mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: fv_print_dpcoup_warn")

   ! Store namelist settings in fv state object
   dyn_state => get_dyn_state()

   dyn_state%grid%high_alt = fv_high_altitude

   ! Calculate nsplit if it was specified as 0
   if ( fv_nsplit <= 0 ) then
      dt = get_step_size()
      dyn_state%nsplit= init_nsplit(dt, plon, plat)
   else
      dyn_state%nsplit= fv_nsplit
   end if

   ! Calculate nspltrac if it was specified as 0
   if (fv_nspltrac <= 0) then
      dyn_state%nspltrac = max (1, dyn_state%nsplit/4)
   else
      dyn_state%nspltrac = fv_nspltrac
   end if

   ! Set nspltvrm to 1 if it was specified as 0
   if (fv_nspltvrm <= 0) then
      dyn_state%nspltvrm = 1
   else
      dyn_state%nspltvrm = fv_nspltvrm
   end if

   dyn_state%iord = fv_iord
   dyn_state%jord = fv_jord
   dyn_state%kord = fv_kord

   ! Calculation of orders for the C grid is fixed by D-grid IORD, JORD
   if( fv_iord <= 2 ) then
      dyn_state%icd =  1
   else
      dyn_state%icd = -2
   end if

   if( fv_jord <= 2 ) then
      dyn_state%jcd =  1
   else
      dyn_state%jcd =  -2
   end if

   dyn_state%consv         = fv_conserve
   dyn_state%filtcw        = fv_filtcw
   dyn_state%fft_flt       = fv_fft_flt
   dyn_state%div24del2flag = fv_div24del2flag
   dyn_state%del2coef      = fv_del2coef

   dyn_state%high_order_top= fv_high_order_top  
   dyn_state%am_correction = fv_am_correction
   dyn_state%am_geom_crrct = fv_am_geom_crrct
   dyn_state%am_fixer      = fv_am_fixer
   dyn_state%am_fix_lbl    = fv_am_fix_lbl
   dyn_state%am_diag       = fv_am_diag


   ! There is a mod for the AM correction in the vertical diffusion code.  Make use
   ! of the physics control module to communicate whether correction is to be applied there.
   call phys_setopts(fv_am_correction_in=fv_am_correction)

   if (masterproc) then
      write(iulog,*)'FV dycore configuration:'
      write(iulog,*)'  Lagrangian time splits (fv_nsplit)                = ', fv_nsplit
      write(iulog,*)'  Tracer time splits (fv_nslptrac)                  = ', fv_nspltrac
      write(iulog,*)'  Vertical re-mapping time splits (fv_nspltvrm)     = ', fv_nspltvrm
      write(iulog,*)'  Scheme in E-W direction (fv_iord)                 = ', fv_iord
      write(iulog,*)'  Scheme in N-S direction (fv_jord)                 = ', fv_jord
      write(iulog,*)'  Scheme for vertical mapping (fv_kord)             = ', fv_kord
      write(iulog,*)'  Conservative dynamics (fv_conserve)               = ', fv_conserve
      write(iulog,*)'  Filtering c-grid winds (fv_filcw)                 = ', fv_filtcw
      write(iulog,*)'  FFT filter (fv_fft_flt)                           = ', fv_fft_flt
      write(iulog,*)'  Divergence/velocity damping (fv_div24del2flag)    = ', fv_div24del2flag
      write(iulog,*)'  Coef for 2nd order velocity damping (fv_del2coef) = ', fv_del2coef
      write(iulog,*)'  High-order top                                     = ', fv_high_order_top
      write(iulog,*)'  Geometry & pressure corr. for AM (fv_am_geom_crrct) = ', fv_am_geom_crrct
      write(iulog,*)'  Angular momentum (AM) correction (fv_am_correction) = ', fv_am_correction
      write(iulog,*)'  Apply AM fixer (fv_am_fixer)                        = ', fv_am_fixer
      write(iulog,*)'  Level by level AM fixer (fv_am_fix_lbl)             = ', fv_am_fix_lbl
      write(iulog,*)'  Enable AM diagnostics (fv_am_diag)                  = ', fv_am_diag
      write(iulog,*)'  '
   end if

   call spmd_readnl(nlfilename)

   call ctem_readnl(nlfilename)
   call fill_readnl(nlfilename)

   !---------------------------------------------------------------------------
   contains
   !---------------------------------------------------------------------------

   integer function init_nsplit(dtime, im, jm)

      !-----------------------------------------------------------------------
      ! find proper value for nsplit if not specified
      !
      ! If nsplit=0 (module variable) then determine a good value
      ! for ns (used in dynpkg) based on resolution and the large-time-step
      ! (pdt). The user may have to set this manually if instability occurs.
      !
      ! REVISION HISTORY:
      !   00.10.19   Lin     Creation
      !   01.06.10   Sawyer  Modified for dynamics_init framework
      !   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
      !-----------------------------------------------------------------------

      ! arguments
      real (r8), intent(in) :: dtime      !  time step
      integer,   intent(in) :: im, jm     !  Global horizontal resolution

      ! LOCAL VARIABLES:
      real (r8)   pdt                       ! Time-step in seconds
      real (r8)   dim
      real (r8)   dim0                      ! base dimension
      real (r8)   dt0                       ! base time step
      real (r8)   ns0                       ! base nsplit for base dimension
      real (r8)   ns                        ! final value to be returned

      parameter ( dim0 = 191._r8  )
      parameter ( dt0  = 1800._r8 )
      parameter ( ns0  = 4._r8    )
      !-----------------------------------------------------------------------

      pdt = int(dtime)   ! dtime is a variable internal to this module
      dim = max ( im, 2*(jm-1) )
      ns  = int ( ns0*abs(pdt)*dim/(dt0*dim0) + 0.75_r8 )
      ns  = max ( 1._r8, ns )   ! for cases in which dt or dim is too small

      init_nsplit = ns

   end function init_nsplit
   !---------------------------------------------------------------------------

end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_register()

   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use ppgrid,          only: pcols, pver
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

   ! Initialize FV dynamical core state variables

   use physconst,       only: pi, omega, rearth, rair, cpair, zvir
   use infnan,          only: inf, assignment(=)

   use constituents,    only: pcnst, cnst_name, cnst_longname, tottnam, cnst_get_ind
   use cam_history,     only: addfld, add_default, horiz_only
   use phys_control,    only: phys_getopts

#if ( defined OFFLINE_DYN )
   use metdata,         only: metdata_dyn_init
#endif
   use ctem,            only: ctem_init
   use diag_module,     only: fv_diag_init

   ! arguments:
   type (dyn_import_t),     intent(out) :: dyn_in
   type (dyn_export_t),     intent(out) :: dyn_out

   ! Local variables
   type (t_fvdycore_state),     pointer :: dyn_state
   type (t_fvdycore_grid),      pointer :: grid
   type (t_fvdycore_constants), pointer :: constants

   real(r8) :: dt

   integer :: ifirstxy, ilastxy
   integer :: jfirstxy, jlastxy
   integer :: km
   integer :: ierr

   integer :: m, ixcldice, ixcldliq
   logical :: history_budget      ! output tendencies and state variables for budgets
   integer :: budget_hfile_num

   character(len=*), parameter :: sub='dyn_init'
   !----------------------------------------------------------------------------

   dyn_state => get_dyn_state()
   grid      => dyn_state%grid
   constants => dyn_state%constants

   if (grid%high_alt) then
      grid%ntotq = grid%ntotq + 1 ! advect Kappa
      grid%kthi = grid%kthi + 1
      grid%kthia(:) = grid%kthia(:) + 1
   endif

   ! Set constants
   constants%pi    = pi
   constants%omega = omega
   constants%ae    = rearth
   constants%rair  = rair
   constants%cp    = cpair
   constants%cappa = rair/cpair
   constants%zvir  = zvir

   dt = get_step_size()
   dyn_state%dt        = dt         ! Should this be part of state??

   dyn_state%check_dt  = 21600.0_r8 ! Check max and min every 6 hours.

   ! Create the dynamics import and export state objects
   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy
   km       = grid%km

   allocate(dyn_in%phis(  ifirstxy:ilastxy,jfirstxy:jlastxy),          &
            dyn_in%ps(    ifirstxy:ilastxy,jfirstxy:jlastxy),          &
            dyn_in%u3s(   ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%v3s(   ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%pe(    ifirstxy:ilastxy,km+1,jfirstxy:jlastxy),     &
            dyn_in%pt(    ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%t3(    ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%pk(    ifirstxy:ilastxy,jfirstxy:jlastxy,km+1),     &
            dyn_in%pkz(   ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%delp(  ifirstxy:ilastxy,jfirstxy:jlastxy,km),       &
            dyn_in%tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,km, grid%ntotq ), &
            stat=ierr)

   if ( ierr /= 0 ) then
      write(iulog,*) sub//': ERROR: allocating components of dyn_in.  ierr=', ierr
      call endrun(sub//': ERROR: allocating components of dyn_in')
   end if

   dyn_in%phis   = inf
   dyn_in%ps     = inf
   dyn_in%u3s    = inf
   dyn_in%v3s    = inf
   dyn_in%pe     = inf
   dyn_in%pt     = inf
   dyn_in%t3     = inf
   dyn_in%pk     = inf
   dyn_in%pkz    = inf
   dyn_in%delp   = inf
   dyn_in%tracer = inf

   ! Export object has all of these except phis
   dyn_out%phis   => dyn_in%phis
   dyn_out%ps     => dyn_in%ps
   dyn_out%u3s    => dyn_in%u3s
   dyn_out%v3s    => dyn_in%v3s
   dyn_out%pe     => dyn_in%pe
   dyn_out%pt     => dyn_in%pt
   dyn_out%t3     => dyn_in%t3
   dyn_out%pk     => dyn_in%pk
   dyn_out%pkz    => dyn_in%pkz
   dyn_out%delp   => dyn_in%delp
   dyn_out%tracer => dyn_in%tracer

   ! And several more which are not in the import container
   allocate(dyn_out%peln (ifirstxy:ilastxy,km+1,jfirstxy:jlastxy),&
            dyn_out%omga (ifirstxy:ilastxy,km,jfirstxy:jlastxy),  &
            dyn_out%mfx  (ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
            dyn_out%mfy  (ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
            stat=ierr)
   if ( ierr /= 0 ) then
      write(iulog,*) sub//': ERROR: allocating components of dyn_out.  ierr=', ierr
      call endrun(sub//': ERROR: allocating components of dyn_out')
   end if

   if (dyn_state%am_fixer .or. dyn_state%am_diag) then

      allocate( &
         dyn_out%duf3s(ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
         stat=ierr)
      if ( ierr /= 0 ) then
         write(iulog,*) sub//': ERROR: allocating duf3s components of dyn_out.  ierr=', ierr
         call endrun(sub//': ERROR: allocating duf3s components of dyn_out')
      end if
      dyn_out%duf3s= inf
   end if

   if (dyn_state%am_diag) then
      allocate( &
         dyn_out%du3s (ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
         dyn_out%dv3s (ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
         dyn_out%dua3s(ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
         dyn_out%dva3s(ifirstxy:ilastxy,jfirstxy:jlastxy,km),  &
         stat=ierr)
      if ( ierr /= 0 ) then
         write(iulog,*) sub//': ERROR: allocating du3s components of dyn_out.  ierr=', ierr
         call endrun(sub//': ERROR: allocating du3s components of dyn_out')
      end if
      dyn_out%du3s = inf
      dyn_out%dv3s = inf
      dyn_out%dua3s= inf
      dyn_out%dva3s= inf
   end if

   dyn_out%peln = inf
   dyn_out%omga = inf
   dyn_out%mfx  = inf
   dyn_out%mfy  = inf

#if ( defined OFFLINE_DYN )
   call metdata_dyn_init(grid)
#endif

   ! Setup circulation diagnostics
   call ctem_init()

   ! Diagnostics for AM
   if (dyn_state%am_diag) call fv_diag_init()

   call set_phis(dyn_in)

   if (initial_run) then

      ! Read in initial data
      call read_inidat(dyn_in)
      call clean_iodesc_list()

   end if

   ! History output

   call addfld('US',    (/ 'lev' /),'A','m/s','Zonal wind, staggered', gridname='fv_u_stagger')
   call addfld('VS',    (/ 'lev' /),'A','m/s','Meridional wind, staggered', gridname='fv_v_stagger')
   call addfld('US&IC', (/ 'lev' /),'I','m/s','Zonal wind, staggered', gridname='fv_u_stagger')
   call addfld('VS&IC', (/ 'lev' /),'I','m/s','Meridional wind, staggered', gridname='fv_v_stagger')
   call addfld('PS&IC', horiz_only, 'I','Pa', 'Surface pressure', gridname='fv_centers')
   call addfld('T&IC',  (/ 'lev' /),'I','K',  'Temperature', gridname='fv_centers')
   do m = 1, pcnst
      call addfld(trim(cnst_name(m))//'&IC',(/ 'lev' /),'I','kg/kg', cnst_longname(m), gridname='fv_centers')
   end do
   do m = 1, pcnst
      call addfld(tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency ',  &
                  gridname='fv_centers')
   end do

   call add_default('US&IC   ', 0, 'I')
   call add_default('VS&IC   ', 0, 'I')
   call add_default('PS&IC      ',0, 'I')
   call add_default('T&IC       ',0, 'I')
   do m = 1, pcnst
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
   end do

   call addfld('DUH',      (/ 'lev' /), 'A','K/s','U horizontal diffusive heating', gridname='fv_centers')
   call addfld('DVH',      (/ 'lev' /), 'A','K/s','V horizontal diffusive heating', gridname='fv_centers')
   call addfld('ENGYCORR', (/ 'lev' /), 'A','W/m2','Energy correction for over-all conservation', &
                                                   gridname='fv_centers')

   call addfld('FU',       (/ 'lev' /), 'A','m/s2','Zonal wind forcing term', gridname='fv_centers')
   call addfld('FV',       (/ 'lev' /), 'A','m/s2','Meridional wind forcing term', gridname='fv_centers')
   call addfld('FU_S',     (/ 'lev' /), 'A','m/s2','Zonal wind forcing term on staggered grid', &
                                                   gridname='fv_u_stagger')
   call addfld('FV_S',     (/ 'lev' /), 'A','m/s2','Meridional wind forcing term on staggered grid', &
                                                   gridname='fv_v_stagger')
   call addfld('TTEND',    (/ 'lev' /), 'A','K/s','Total T tendency (all processes)', gridname='fv_centers')
   call addfld('LPSTEN',   horiz_only,  'A','Pa/s','Surface pressure tendency',       gridname='fv_centers')
   call addfld('VAT',      (/ 'lev' /), 'A','K/s','Vertical advective tendency of T', gridname='fv_centers')
   call addfld('KTOOP',    (/ 'lev' /), 'A','K/s','(Kappa*T)*(omega/P)',              gridname='fv_centers')

   call phys_getopts(history_budget_out=history_budget, history_budget_histfile_num_out=budget_hfile_num)
   if ( history_budget ) then
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      call add_default(tottnam(       1), budget_hfile_num, ' ')
      call add_default(tottnam(ixcldliq), budget_hfile_num, ' ')
      call add_default(tottnam(ixcldice), budget_hfile_num, ' ')
      call add_default('TTEND   '       , budget_hfile_num, ' ')
   end if

end subroutine dyn_init

!=============================================================================================

subroutine dyn_run(ptop, ndt, te0, dyn_state, dyn_in, dyn_out, rc)

   ! Driver for the NASA finite-volume dynamical core
   !
   ! Developer: Shian-Jiann Lin, NASA/GSFC; email: lin@dao.gsfc.nasa.gov
   !
   ! Top view of D-grid prognostatic variables: u, v, and delp (and other scalars)
   !
   !               u(i,j+1)
   !                 |
   !      v(i,j)---delp(i,j)---v(i+1,j)
   !                 |
   !               u(i,j)
   !
   ! External routine required:
   !
   ! The user needs to supply a subroutine to set up
   ! Eulerian vertical coordinate" for remapping purpose.
   ! Currently this routine is named as set_eta()
   ! In principle any terrian following vertical
   ! coordinate can be used. The input to fvcore
   ! need not be on the same vertical coordinate
   ! as the output.
   ! If SPMD is defined the Pilgrim communication
   ! library developed by Will Sawyer will be needed.
   !
   ! Remarks:
   !
   ! Values at poles for both u and v need not be defined; but values for
   ! all other scalars needed to be defined at both poles (as polar cap mean
   ! quantities). Tracer advection is done "off-line" using the
   ! large time step. Consistency is maintained by using the time accumulated
   ! Courant numbers and horizontal mass fluxes for the FFSL algorithm.
   ! The input "pt" can be either dry potential temperature
   ! defined as T/pkz (adiabatic case) or virtual potential temperature
   ! defined as T*/pkz (full phys case). IF convt is true, pt is not updated.
   ! Instead, virtual temperature is ouput.
   ! ipt is updated if convt is false.
   ! The user may set the value of nx to optimize the SMP performance
   ! The optimal valuse of nx depends on the total number of available
   ! shared memory CPUs per node (NS). Assuming the maximm MPI
   ! decomposition is used in the y-direction, set nx=1 if the
   ! NS <=4; nx=4 if NS=16.
   !
   ! A 2D xy decomposition is used for handling the Lagrangian surface
   ! remapping, the ideal physics, and (optionally) the geopotential
   ! calculation.
   !
   ! The transpose from yz to xy decomposition takes place within dynpkg.
   ! The xy decomposed variables are then transposed directly to the
   ! physics decomposition within d_p_coupling.
   !
   ! The xy decomposed variables have names corresponding to the
   ! yz decomposed variables: simply append "xy". Thus, "uxy" is the
   ! xy decomposed version of "u".
   !
   ! This version supports overlap of trac2d and cd_core subcycles (Art Mirin, November 2007).
   ! This refers to the subcycles described by the "do n=1,n2" loop and has nothing to
   ! do with the "do it=1,nsplit" lower-level subcycling. Each trac2d call (n), other than the last,
   ! is overlapped with the subsequent cd_core 'series' (n+1). The controlling namelist variable
   ! is ct_overlap. The overlapping trac2d calls are carried out on the second set of
   ! npes_yz processes (npes_yz <= iam < 2*npes_yz). The tracer arrays are sent to the
   ! auxiliary processes prior to the do n=1,n2 loop. During each subcycle (other than the last),
   ! the dp0 array is sent prior to the cd_core series; arrays cx, cy, mfx, mfy are sent directly
   ! from cd_core during the last call in the series (it=nsplit). At the completion of the last
   ! auxiliary trac2d subcycle (n=n2-1), the updated tracer values are returned to the
   ! primary processes; the last tracer subcycle (n=n2) is carried out on the primary processes.
   ! Communication calls are nonblocking, with attempt to overlap computation to the extent
   ! possible. The CCSM mpi layer (wrap_mpi) is used. Tags with values greater than npes_xy
   ! are chosen to avoid possible interference between the messages sent from cd_core and
   ! the geopk-related transpose messages called from cd_core thereafter. The auxiliary
   ! processes must use values of jfirst, jlast, kfirst, klast corresponding to their primary
   ! process antecedents, whereas by design those values are (1,0,1,0), resp. (set in spmdinit_dyn).
   ! We therefore add auxiliary subdomain limits to the grid datatype: jfirstct, jlastct,
   ! kfirstct, klastct. For the primary processes, these are identical to the actual subdomain
   ! limits; for the secondary processes, these correspond to the subdomain limits of the
   ! antecedent primary process. These values are communicated to the auxiliary processes
   ! during initialization (spmd_vars_init). During the auxiliary calculations (and allocations)
   ! we temporarily set jfirst equal to jfirstct (etc.) and when done, restore to the original
   ! values. Other information needed by the auxiliary processes is obtained through the grid
   ! datatype.
   !
   ! This version supports tracer decomposition with trac2d (Art Mirin, January 2008).
   ! This option is mutually exclusive with ct_overlap. Variable "trac_decomp" is the size of the
   ! decomposition. The tracers are divided into trac_decomp groups, and the kth group is solved
   ! on the kth set of npes_yz processes. Much of the methodology is similar to that for ct_overlap.
   !
   ! REVISION HISTORY:
   !   SJL 99.04.13:  Initial SMP version delivered to Will Sawyer
   !   WS  99.10.03:  1D MPI completed and tested;
   !   WS  99.10.11:  Additional documentation
   !   WS  99.10.28:  benergy and te_map added; arrays pruned
   !   SJL 00.01.01:  SMP and MPI enhancements; documentation
   !   WS  00.07.13:  Changed PILGRIM API
   !   WS  00.08.28:  SPMD instead of MPI_ON
   !   AAM 00.08.10:  Add kfirst:klast
   !   WS  00.12.19:  phis now distr., LLNL2DModule initialized here
   !   WS  01.02.02:  bug fix: parsplit only called for FIRST time
   !   WS  01.04.09:  Added initialization of ghost regions
   !   WS  01.06.10:  Removed if(first) section; use module
   !   AAM 01.06.27:  Extract te_map call into separate routine
   !   AAM 01.07.13:  Get rid of dynpkg2; recombine te_map;
   !                  perform forward transposes for 2D decomposition
   !   WS  01.12.10:  Ghosted PT (changes benergy, cd_core, te_map, hswf)
   !   WS  03.08.05:  removed vars dcaf, rayf, ideal, call to hswf
   !                  (idealized physics is now in physics package)
   !   WS  03.08.13:  Removed ghost region from UXY
   !   WS  05.06.11:  Inserted into FVCAM_GridCompMod
   !   WS  06.03.03:  Added dyn_state as argument (for reentrancy)
   !   WS  06.06.28:  Using new version of benergy
   !   TT  16.12.11:  AM conservation options
   !   TT  17.30.01:  dynamic wind increments diagnostic
   !-----------------------------------------------------------------------


   use diag_module, only   : compute_vdot_gradp

#if defined( SPMD )
   use mpishorthand,   only: mpir8
   use mod_comm,       only: mp_sendirr, mp_recvirr, mp_send4d_ns, &
                             mp_send3d, mp_recv3d, &
                             mp_recv4d_ns, mp_sendtrirr, mp_recvtrirr
#endif
#if ( defined OFFLINE_DYN )
   use metdata, only: get_met_fields, advance_met, get_us_vs, met_rlx
   use pfixer,  only: adjust_press
#endif
   use metdata, only: met_fix_mass

   use shr_reprosum_mod, only: shr_reprosum_calc
   use physconst,        only: physconst_calc_kappav

#if defined( SPMD )
#include "mpif.h"
#endif

   ! arguments
   real(r8),            intent(in)  :: ptop      ! Pressure at model top (interface pres)
   integer,             intent(in)  :: ndt       ! the large time step in seconds
                                                 ! Also the mapping time step in this setup

   real(r8),            intent(out) :: te0       ! Total energy before dynamics
   type (T_FVDYCORE_STATE), target  :: dyn_state ! Internal state
   type (dyn_import_t)              :: dyn_in    ! Import container
   type (dyn_export_t)              :: dyn_out   ! Export container

   integer,             intent(out) :: rc        ! Return code

   integer, parameter  ::  DYN_RUN_SUCCESS           = 0
   integer, parameter  ::  DYN_RUN_FAILURE           = -1
   integer, parameter  ::  DYN_RUN_R4_NOT_SUPPORTED  = -10
   integer, parameter  ::  DYN_RUN_MUST_BE_2D_DECOMP = -20
   real(r8), parameter ::  D1_0                      = 1.0_r8

   ! Variables from the dynamics interface (import or export)

   real(r8), pointer :: phisxy(:,:)   ! surface geopotential (grav*zs)
   real(r8), pointer :: psxy(:,:)     ! Surface pressure (pa)
   real(r8), pointer :: t3xy(:,:,:)   ! temperature (K)
   real(r8), pointer :: ptxy(:,:,:)   ! scaled (virtual) potential temperature
   real(r8), pointer :: delpxy(:,:,:) ! Pressure thickness
   real(r8), pointer :: tracer(:,:,:,:) ! Tracers
   real(r8), pointer :: uxy(:,:,:)    ! u wind velocities, staggered grid
   real(r8), pointer :: vxy(:,:,:)    ! v wind velocities, staggered grid

   !--------------------------------------------------------------------------------------
   ! The arrays pexy, pkxy, pkzxy must be pre-computed as input to benergy().
   ! They are NOT needed if dyn_state%consv=.F.; updated on output (to be used
   ! by physdrv) Please refer to routine pkez on the algorithm for computing pkz
   ! from pe and pk
   !--------------------------------------------------------------------------------------

   real(r8), pointer :: pexy(:,:,:)   ! Pres at layer edges
   real(r8), pointer :: pkxy(:,:,:)   ! pe**cappa
   real(r8), pointer :: pkzxy(:,:,:)  ! finite-volume mean of pk

   ! Export state only variables
   real(r8), pointer :: pelnxy(:,:,:)   ! Natural logarithm of pe
   real(r8), pointer :: omgaxy(:,:,:)   ! vertical pressure velocity (pa/sec)
   real(r8), pointer :: mfxxy(:,:,:)    ! mass flux in X (Pa m^\2 / s)
   real(r8), pointer :: mfyxy(:,:,:)    ! mass flux in Y (Pa m^\2 / s)
   real(r8), pointer :: duxy(:,:,:)     ! u tot. tend. from dycore, staggered grid
   real(r8), pointer :: dvxy(:,:,:)     ! v tot. tend. from dycore, staggered grid
   real(r8), pointer :: ucxy(:,:,:)     ! u tend. from advection only, staggd grid
   real(r8), pointer :: vcxy(:,:,:)     ! v tend. from advection only, staggd grid
   real(r8), pointer :: dufix_xy(:,:,:) ! u tend. from AM fixer, staggered grid

   ! Other pointers (for convenience)
   type (T_FVDYCORE_GRID)      , pointer :: GRID      ! For convenience
   type (T_FVDYCORE_CONSTANTS) , pointer :: CONSTANTS ! For convenience

   !  YZ variables currently allocated on stack... should they be on the heap?

   real(r8) :: ps(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: phis(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: pe(dyn_state%grid%im,  &
                  dyn_state%grid%kfirst:dyn_state%grid%klast+1,&
                  dyn_state%grid%jfirst:dyn_state%grid%jlast)
   real(r8) :: delp(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                    dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: pk(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                  dyn_state%grid%kfirst:dyn_state%grid%klast+1)
   real(r8) :: pkz(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast, &
                   dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: u(dyn_state%grid%im,   &
                 dyn_state%grid%jfirst-dyn_state%grid%ng_d:dyn_state%grid%jlast+dyn_state%grid%ng_s,&
                 dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: v(dyn_state%grid%im,   &
                 dyn_state%grid%jfirst-dyn_state%grid%ng_s:dyn_state%grid%jlast+dyn_state%grid%ng_d,&
                 dyn_state%grid%kfirst:dyn_state%grid%klast)
   real(r8) :: pt(dyn_state%grid%im,  &
                  dyn_state%grid%jfirst-dyn_state%grid%ng_d:dyn_state%grid%jlast+dyn_state%grid%ng_d,&
                  dyn_state%grid%kfirst:dyn_state%grid%klast)

   real(r8) :: pi
   real(r8) :: om       ! angular velocity of earth's rotation
   real(r8) :: cp       ! heat capacity of air at constant pressure
   real(r8) :: ae       ! radius of the earth (m)

   real(r8) :: rair     ! Gas constant of the air
   real(r8) :: cappa    ! R/Cp
   real(r8) :: zvir     ! Virtual effect constant ( = rwv/rair-1 )

   logical :: consv     ! Energy conserved?

   integer :: im        ! dimension in east-west
   integer :: jm        ! dimension in North-South
   integer :: km        ! number of Lagrangian layers
   integer :: jfirst    ! starting latitude index for MPI
   integer :: jlast     ! ending latitude index for MPI
   integer :: kfirst    ! starting vertical index for MPI
   integer :: klast     ! ending vertical index for MPI
   integer :: ntotq     ! total # of tracers to be advected
   integer :: iord      ! parameter controlling monotonicity in E-W
                                   ! recommendation: iord=4
   integer :: jord      ! parameter controlling monotonicity in N-S
                                   ! recommendation: jord=4
   integer :: kord      ! parameter controlling monotonicity in mapping
                                   ! recommendation: kord=4
   integer :: icd       ! X algorithm order on C-grid
   integer :: jcd       ! Y algorithm order on C-grid
   integer :: ng_d      ! Ghosting width on D-grid
   integer :: ng_s      ! Ghosting width (staggered, for winds)
   integer :: ns        ! overall split
   integer :: div24del2flag
   real(r8):: del2coef

   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy  ! xy decomposition
   integer :: npr_z

   logical :: cd_penul

   real(r8), allocatable, target    :: q_internal(:,:,:,:)    ! Pointers to tracers
   integer i, j, k, iq          ! Loop indicies
   real(r8) umax                ! Maximum winds, m/s
   parameter (umax = 300.0_r8)

   integer    nx          ! # of split pieces in x-direction; for performance, the
#if defined( UNICOSMP )
   parameter (nx = 1)
#else
   parameter (nx = 4)     ! user may set nx=1 if there is NO shared memory multitasking
#endif
   integer :: ipe, it, iv
   integer :: nsplit, n, n2, nv
   integer :: mq

   ! Move the following 3D arrays to an initialization routine?
   real(r8), allocatable :: worka(:,:,:),workb(:,:,:),dp0(:,:,:),cx(:,:,:),cy(:,:,:)
   real(r8), allocatable :: mfx(:,:,:), mfy(:,:,:)
   real(r8), allocatable :: delpf(:,:,:), uc(:,:,:), vc(:,:,:)
   real(r8), allocatable :: dwz(:,:,:), pkc(:,:,:), wz(:,:,:)
   real(r8), allocatable :: dpt(:,:,:)
   real(r8), allocatable :: pkcc(:,:,:), wzc(:,:,:)

   ! The following variables are work arrays for xy=>yz transpose
   real(r8), allocatable :: pkkp(:,:,:), wzkp(:,:,:)

   ! The following variables are xy instantiations
   real(r8), allocatable :: tempxy(:,:,:), dp0xy(:,:,:), wzxy(:,:,:)

   ! psxy3 is dummy 3d variant of psxy
   real(r8), allocatable :: psxy3(:,:,:)

   ! phisxy3 is dummy 3d variant of phisxy
   real(r8), allocatable :: phisxy3(:,:,:)
   real(r8), pointer     :: q3xypt(:,:,:)
   real(r8), pointer     :: q3yzpt(:,:,:)
   real(r8)              :: tte(dyn_state%grid%jm)
   real(r8)              :: XXX(dyn_state%grid%km)

#if ( defined OFFLINE_DYN )
   real(r8), allocatable :: ps_obs(:,:)
   real(r8), allocatable :: ps_mod(:,:)
   real(r8), allocatable :: u_tmp(:,:,:)
   real(r8), allocatable :: v_tmp(:,:,:)
#endif

   logical :: fill

   real(r8) :: dt
   real(r8) :: bdt
   integer  :: filtcw
   integer  :: ct_overlap
   integer  :: trac_decomp

   integer modc_tracers, mlast

   ! cd_core / trac2d overlap and tracer decomposition data (AAM)
   integer :: commnyz                                         ! n*npes_yz communicator
   integer :: jfirstct, jlastct, kfirstct, klastct            ! primary subdomain limits
   integer :: jkstore(4)                                      ! storage for subdomain limits
   integer :: iamlocal                                        ! task number (global indexing)
   integer :: iremotea(dyn_state%grid%trac_decomp)            ! source/target; id array
   integer :: iremote                                         ! source/target; working id
   integer :: ndp0, ncx, ncy, nmfx, nmfy, ntrac               ! message sizes
   integer :: dp0tag, cxtag, cytag, mfxtag, mfytag, tractag   ! message tags
   integer :: cxtaga(dyn_state%grid%trac_decomp)              ! tag arrays for cd_core
   integer :: cytaga(dyn_state%grid%trac_decomp)              ! tag arrays for cd_core
   integer :: mfxtaga(dyn_state%grid%trac_decomp)             ! tag arrays for cd_core
   integer :: mfytaga(dyn_state%grid%trac_decomp)             ! tag arrays for cd_core
   logical :: ct_aux                                          ! true if auxiliary process
   logical :: s_trac                                          ! true for cd_core posting tracer-related sends
   integer, allocatable :: ctreq(:,:)                         ! used for nonblocking receive
   integer, allocatable :: ctstat(:,:,:)                      ! used for nonblocking receive
   integer, allocatable :: ctreqs(:,:)                        ! used for nonblocking send
   integer, allocatable :: ctstats(:,:,:)                     ! used for nonblocking send
   integer, allocatable :: cdcreqs(:,:)                       ! used for nonblocking send in cd_core
   integer, pointer :: ktloa(:)                               ! lower limit of tracer decomposition (global)
   integer, pointer :: kthia(:)                               ! upper limit of tracer decomposition (global)
   integer ktlo                                               ! lower limit of tracer decomposition (local)
   integer kthi                                               ! upper limit of tracer decomposition (local)
   integer kt, tagu, naux, kaux, ntg0

   logical :: print_subcycling = .true.
   logical :: c_dotrac, t_dotrac
   logical :: convt_local

   data fill  /.true./              ! perform a simple filling algorithm
                                    ! in case negatives were found

   ! C.-C. Chen, omega calculation
   real(r8) :: cx_om(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast, &
                     dyn_state%grid%kfirst:dyn_state%grid%klast)        ! Courant no. in X
   real(r8) :: cy_om(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast+1, &
                     dyn_state%grid%kfirst:dyn_state%grid%klast)        ! Courant no. in Y
   real(r8) :: pexy_om(dyn_state%grid%ifirstxy:dyn_state%grid%ilastxy,dyn_state%grid%km+1, &
                       dyn_state%grid%jfirstxy:dyn_state%grid%jlastxy)

   ! Non-constant air properties for high top models (waccmx).
   real(r8) :: cap3vi(dyn_state%grid%ifirstxy:dyn_state%grid%ilastxy,&
                      dyn_state%grid%jfirstxy:dyn_state%grid%jlastxy,dyn_state%grid%km+1)
   real(r8) :: cp3vc (dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                      dyn_state%grid%kfirst:dyn_state%grid%klast)         !C_p on yz
   real(r8) :: cap3vc(dyn_state%grid%im,dyn_state%grid%jfirst:dyn_state%grid%jlast,&
                      dyn_state%grid%kfirst:dyn_state%grid%klast)        !cappa on yz

   real(r8), dimension(dyn_state%grid%ifirstxy:dyn_state%grid%ilastxy,&
                       dyn_state%grid%jfirstxy:dyn_state%grid%jlastxy,dyn_state%grid%km) :: &
             cp3v,cap3v
   logical :: high_alt

   ! angular momentum (AM) conservation
   logical  :: am_correction         ! apply AM correction?
   logical  :: am_geom_crrct         ! apply AM geom. corr?
   logical  :: am_fixer              ! apply AM fixer?
   logical  :: am_fix_lbl            ! apply fixer separately on each shallow-water layer?
   logical  :: am_fix_taper=.false.  ! def. no tapering; modified if global fixer applied or high_order_top=.false.
   real(r8) :: tmpsum(1,2)
   real(r8) :: tmpresult(2)
   real(r8) :: am0, am1, me0

   real(r8) :: don(dyn_state%grid%jm,dyn_state%grid%km), & ! out of cd_core
               dod(dyn_state%grid%jm,dyn_state%grid%km)    ! out of cd_core
   real(r8) :: dons(dyn_state%grid%km), &                  ! sums over j
               dods(dyn_state%grid%km)

   real(r8), allocatable :: zpkck(:,:)
   real(r8) :: avgpk(dyn_state%grid%km)
   real(r8) :: taper(dyn_state%grid%km)
   real(r8) :: ptapk, xdlt2
   real(r8), parameter :: ptap =9000._r8
   real(r8), parameter :: dptap=1000._r8
   real(r8), parameter :: tiny=.1e-10_r8

   ! AM diagnostics
   logical  :: am_diag                  ! enable angular momentum diagnostic output
   logical  :: am_fix_out
   integer  :: kmtp                     ! range of levels (1:kmtp) where order is reduced
   real(r8) :: ame(dyn_state%grid%jm)
   real(r8) :: zpe(dyn_state%grid%jfirstxy:dyn_state%grid%jlastxy)
   real(r8) :: tmp
   real(r8) :: du_fix_g
   real(r8) :: du_fix(dyn_state%grid%km)
   real(r8) :: du_fix_s(dyn_state%grid%km)
   real(r8), allocatable :: du_fix_i(:,:,:)
   real(r8), allocatable :: du_k    (:,:)
   real(r8), allocatable :: du_north(:,:)
   real(r8), allocatable :: uc_s(:,:,:),vc_s(:,:,:)  ! workspace (accumulated uc,vc)
   real(r8), allocatable :: uc_i(:,:,:),vc_i(:,:,:)  ! workspace (transposed uc_s,vc_s)


   ! NOTE -- model behaviour with high_order_top=true is still under validation and may require
   !         some other form of enhanced damping in the top layer
   logical :: high_order_top

   !--------------------------------------------------------------------------------------
   kmtp=dyn_state%grid%km/8

   rc       =  DYN_RUN_FAILURE      ! Set initially to fail

   phisxy   => dyn_in%phis
   psxy     => dyn_in%ps
   uxy      => dyn_in%u3s
   vxy      => dyn_in%v3s
   t3xy     => dyn_in%t3
   ptxy     => dyn_in%pt
   delpxy   => dyn_in%delp
   tracer   => dyn_in%tracer
   pexy     => dyn_in%pe
   pkxy     => dyn_in%pk
   pkzxy    => dyn_in%pkz

   pelnxy   => dyn_out%peln
   omgaxy   => dyn_out%omga
   mfxxy    => dyn_out%mfx
   mfyxy    => dyn_out%mfy
   duxy     => dyn_out%du3s
   dvxy     => dyn_out%dv3s
   ucxy     => dyn_out%dua3s
   vcxy     => dyn_out%dva3s
   dufix_xy => dyn_out%duf3s

   grid => dyn_state%grid    ! For convenience
   constants => DYN_STATE%CONSTANTS

   ns   = dyn_state%nsplit   ! large split (will be subdivided later)
   n2   = dyn_state%nspltrac ! tracer split(will be subdivided later)
   nv   = dyn_state%nspltvrm ! vertical re-mapping split
   icd  = dyn_state%icd
   jcd  = dyn_state%jcd
   iord = dyn_state%iord
   jord = dyn_state%jord
   kord = dyn_state%kord
   div24del2flag = dyn_state%div24del2flag
   del2coef      = dyn_state%del2coef
   filtcw        = dyn_state%filtcw
   high_alt      = grid%high_alt

   consv      = dyn_state%consv
   high_order_top= dyn_state%high_order_top
   am_correction = dyn_state%am_correction
   am_geom_crrct = dyn_state%am_geom_crrct
   am_fixer   = dyn_state%am_fixer
   am_fix_lbl = dyn_state%am_fix_lbl
   am_diag    = dyn_state%am_diag

   pi   =  constants%pi
   om   =  constants%omega
   ae   =  constants%ae
   rair =  constants%rair
   cp   =  constants%cp
   cappa=  constants%cappa
   zvir =  constants%zvir

   im = grid%im
   jm = grid%jm
   km = grid%km

   ng_d  = grid%ng_d
   ng_s  = grid%ng_s

   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   jfirst   = grid%jfirst
   jlast    = grid%jlast
   kfirst   = grid%kfirst
   klast    = grid%klast

   ntotq    = grid%ntotq
   modc_tracers = grid%modc_tracers

   npr_z    = grid%npr_z

   ! cd_core/trac2d overlap and tracer decomposition
   ct_overlap  = grid%ct_overlap
   trac_decomp = grid%trac_decomp
   jfirstct    = grid%jfirstct
   jlastct     = grid%jlastct
   kfirstct    = grid%kfirstct
   klastct     = grid%klastct
   commnyz     = grid%commnyz
   iamlocal    = grid%iam

   ! kaux is an index describing the set of npes_yz processes; 0 for first set, 1 for second set, etc.
   kaux = iamlocal/grid%npes_yz

   ! ct_aux is true if current process is auxiliary, false otherwise
   ct_aux = ((ct_overlap .gt. 0 .and. kaux .eq. 1) .or.      &
             (trac_decomp .gt. 1 .and. kaux .ge. 1 .and. kaux .lt. trac_decomp))

   ! define message tags to exceed npes_xy so as not to interfere with geopotential transpose tags
   ! tags below correspond to communicated variables with ct_overlap and trac_decomp
   dp0tag  = grid%npes_xy + 5
   cxtag   = dp0tag + 1
   cytag   = dp0tag + 2
   mfxtag  = dp0tag + 3
   mfytag  = dp0tag + 4
   tractag = dp0tag + 5

   ! ntg0 is upper bound on number of needed tags beyond tracer tags for ct_overlap and trac_decomp
   ntg0 = 10

   ! set am_fix tapering parameters
   if (am_fixer.and..not.am_fix_lbl) then
      am_fix_taper = .true.   ! always apply tapering with global fixer
      ptapk        = ptap**constants%cappa
      xdlt2        = 2._r8/(log((ptap+.5_r8*dptap)/(ptap-.5_r8*dptap))*constants%cappa)
   end if

#if ( defined OFFLINE_DYN )

   ! advance the meteorology data
   call advance_met(grid)

   ! set the staggered winds (verticity winds) to offline meteorological data
   call get_us_vs( grid, u, v )
#endif

   if (high_alt) then
      call physconst_calc_kappav(ifirstxy,ilastxy,jfirstxy,jlastxy,1,km, grid%ntotq, tracer, cap3v, cpv=cp3v )
   else
      cp3v  = cp
      cp3vc = cp
      cap3v = cappa
      cap3vi= cappa
      cap3vc= cappa
   endif

   if ( km > 1 ) then         ! not shallow water equations

      if (consv) then

         if (grid%iam .lt. grid%npes_xy) then

            ! Tests indicate that t3 does not have consistent
            ! pole values, e.g. t3(:,1,k) are not all the same.
            ! Not clear why this is not the case: it may be that the pole
            ! values are not consistent on the restart file.  For the time being,
            ! perform a parallel sum over t3 and correct the pole values

            if ( jfirstxy == 1 ) then
               call par_xsum(grid, t3xy(:,1,:), km, XXX)
               do k = 1, km
                  do i = ifirstxy, ilastxy
                     t3xy(i,1,k) = XXX(k) / real(im,r8)
                  end do
               end do
            end if

            if ( jlastxy == jm ) then
               call par_xsum(grid, t3xy(:,jm,:), km, XXX)
               do k = 1, km
                  do i = ifirstxy, ilastxy
                     t3xy(i,jm,k) = XXX(k) / real(im,r8)
                  end do
               end do
            end if

            if (consv) then
               ! Compute globally integrated Total Energy (te0)
               call t_startf ('benergy')

               call benergy(grid, uxy, vxy, t3xy, delpxy,          &
                            tracer(:,:,:,1), pexy, pelnxy, phisxy, &
                            zvir, cp,  rair, tte, te0)

               call t_stopf('benergy')
            end if

         end if
      end if
   end if


   ! Allocate temporary work arrays
   ! Change later to use pointers for SMP performance???
   ! (prime candidates: uc, vc, delpf)

   call t_startf ('dyn_run_alloc')

   if (ct_aux) then
      ! Temporarily set subdomain limits in auxiliary process equal to those of antecedent
      ! to allow following arrays to have proper size
      ! (Normally, sizes of unneeded arrays for auxiliary processes will be deliberately small.)
      jkstore(1) = jfirst
      jkstore(2) = jlast
      jkstore(3) = kfirst
      jkstore(4) = klast
      jfirst = jfirstct
      jlast  = jlastct
      kfirst = kfirstct
      klast  = klastct
   endif

   allocate( worka(im,jfirst:     jlast,     kfirst:klast) )
   allocate( workb(im,jfirst:     jlast,     kfirst:klast) )
   allocate(   dp0(im,jfirst-1:   jlast,     kfirst:klast) )
   allocate(   mfx(im,jfirst:     jlast,     kfirst:klast) )
   allocate(   mfy(im,jfirst:     jlast+1,   kfirst:klast) )
   allocate(    cx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
   allocate(    cy(im,jfirst:     jlast+1,   kfirst:klast) )
   dp0(:,:,:) = 0._r8
   mfx(:,:,:) = 0._r8
   mfy(:,:,:) = 0._r8
   cx(:,:,:) = 0._r8
   cy(:,:,:) = 0._r8

   if (ct_aux) then
      ! Restore subdomain limits in auxiliary process
      jfirst = jkstore(1)
      jlast  = jkstore(2)
      kfirst = jkstore(3)
      klast  = jkstore(4)
   endif

   allocate( delpf(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
   allocate(    uc(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
   allocate(    vc(im,jfirst-2:   jlast+2,   kfirst:klast) )
   allocate(   dpt(im,jfirst-1:   jlast+1,   kfirst:klast) )
   allocate(   dwz(im,jfirst-1:    jlast,    kfirst:klast+1) )
   allocate(   pkc(im,jfirst-1:   jlast+1,   kfirst:klast+1) )
   allocate(    wz(im,jfirst-1:   jlast+1,   kfirst:klast+1) )
   allocate(  pkcc(im,jfirst  :   jlast  ,   kfirst:klast+1) )
   allocate(   wzc(im,jfirst  :   jlast  ,   kfirst:klast+1) )
   allocate(pkkp(im,jfirst:jlast,kfirst:klast+1))
   allocate(wzkp(im,jfirst:jlast,kfirst:klast+1))
   allocate(wzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1))
   allocate(tempxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
   allocate(dp0xy(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
   allocate(psxy3(ifirstxy:ilastxy,jfirstxy:jlastxy,npr_z))
   allocate(phisxy3(ifirstxy:ilastxy,jfirstxy:jlastxy,npr_z))

#if ( defined OFFLINE_DYN )
   allocate( ps_obs(im,jfirst:jlast) )
   allocate( ps_mod(im,jfirst:jlast) )
   allocate( u_tmp(im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
   allocate( v_tmp(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
#endif


   ! Allocation of tracers

   if (ct_aux) then
      ! Temporarily set subdomain limits in auxiliary process equal to those of antecedent
      ! to allow trac2d temporary storage to have proper size
      jfirst = jfirstct
      jlast  = jlastct
      kfirst = kfirstct
      klast  = klastct
   end if
   allocate ( q_internal(im, jfirst:jlast, kfirst:klast, ntotq) )

   ! Trac2d-related mpi quantities for ct_overlap and tracer decomposition
   allocate (ctreq(ntotq+ntg0,trac_decomp))
   allocate (ctreqs(ntotq+ntg0,trac_decomp))
   allocate (cdcreqs(trac_decomp,4))
   cdcreqs(:,:) = 0
#if defined(SPMD)
   allocate (ctstat(MPI_STATUS_SIZE,ntotq+ntg0,trac_decomp))
   allocate (ctstats(MPI_STATUS_SIZE,ntotq+ntg0,trac_decomp))
#endif

   ! Allocate the variables used in tapering
   if (am_fix_taper) then
      allocate(zpkck(dyn_state%grid%jm,dyn_state%grid%km))
   end if

   ! Allocate fields required for dycore diagnostic
   if (am_fixer .or. am_diag) then
      allocate(du_fix_i(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(du_k    (ifirstxy:ilastxy,jfirstxy:jlastxy+1))
      allocate(du_north(ifirstxy:ilastxy,km))
      allocate(uc_s(im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(vc_s(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(uc_i(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      allocate(vc_i(ifirstxy:ilastxy,jfirstxy:jlastxy,km))
      du_fix_i(:,:,:) = 0._r8
      uc_s (:,:,:)  = 0._r8
      vc_s (:,:,:)  = 0._r8
   end if

   ! Compute i.d.'s of remote processes for ct_overlap or trac_decomp
   naux = 0
   if ((ct_overlap .gt. 0 .and. kaux .lt. 2) .or.      &
      (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)) then

      ! Identify involved processes
      iremotea(:) = -1
      naux = max(1,trac_decomp-1)

      if (kaux .eq. 0) then

         ! Primary process - identify corresponding auxiliary process(es)
         do kt = 1, naux
            iremotea(kt) = iamlocal + kt*grid%npes_yz
            cxtaga(kt) = cxtag + (kt-1)*(ntotq+ntg0)
            cytaga(kt) = cytag + (kt-1)*(ntotq+ntg0)
            mfxtaga(kt) = mfxtag + (kt-1)*(ntotq+ntg0)
            mfytaga(kt) = mfytag + (kt-1)*(ntotq+ntg0)
         end do
      else

         ! Auxiliary process - identify corresponding primary process
         iremotea(1) = iamlocal - kaux*grid%npes_yz
      end if
      iremote = iremotea(1)
      ! Message sizes
      ndp0  = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      ncx   = im*(jlast-jfirst+2*ng_d+1)*(klast-kfirst+1)
      ncy   = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      nmfx  = im*(jlast-jfirst+1       )*(klast-kfirst+1)
      nmfy  = im*(jlast-jfirst+2       )*(klast-kfirst+1)
      ntrac = im*(jlast-jfirst+1       )*(klast-kfirst+1)
   end if

   if (ct_aux) then
      ! Restore subdomain limits in auxiliary process
      jfirst = jkstore(1)
      jlast  = jkstore(2)
      kfirst = jkstore(3)
      klast  = jkstore(4)
   end if

   ! Set tracer limits to be supplied to trac2d (needed even without tracer decomposition)
   ktloa => grid%ktloa
   kthia => grid%kthia
   ktlo = grid%ktlo
   kthi = grid%kthi

   call t_stopf  ('dyn_run_alloc')

   ! Determine splitting
   bdt = ndt

   ! Second/third level splitting (nsplit and n2 variables overloaded)
   n2     = (n2+nv   -1) / nv
   nsplit = (ns+n2*nv-1) / (n2*nv)
   dt     = bdt / real(nsplit*n2*nv,r8)

   if (print_subcycling) then
      print_subcycling = .false.
      if (masterproc) then
         write(iulog,*) 'FV subcycling - nv, n2, nsplit, dt = ', nv, n2, nsplit, dt
         if ( (nsplit*n2*nv /= dyn_state%nsplit) .or. (n2*nv /= dyn_state%nspltrac) ) then
            write(iulog,*) "ERROR:  Because of loop nesting, FV dycore can't use the specified namelist settings for subcycling"
            write(iulog,*) '  The original namelist settings were:'
            write(iulog,*) '  fv_nsplit   = ', dyn_state%nsplit
            write(iulog,*) '  fv_nspltrac = ', dyn_state%nspltrac
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '  fv_nspltvrm = ', dyn_state%nspltvrm
            write(iulog,*)
            write(iulog,*) '  fv_nsplit needs to be a multiple of fv_nspltrac'
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '    which in turn needs to be a multiple of fv_nspltvrm.'
            write(iulog,*) '  Suggested settings would be:'
            write(iulog,*) '  fv_nsplit   = ', nsplit*n2*nv
            write(iulog,*) '  fv_nspltrac = ', n2*nv
            if( dyn_state%nspltvrm /= 1 ) write(iulog,*) '  fv_nspltvrm = ', nv
            call endrun("Bad namelist settings for FV subcycling.")
         end if
      end if
   end if

   ! IF convt_local is false, pt is updated for the next iteration of the iv=1,nv loop
   ! On the last iteration, convt_local is set to convt
   convt_local = .false.

   ! initialise global non-conservation integrals
   am1=0._r8
   me0=1._r8

   if (am_fixer.or.am_diag) then
      du_fix_g        = 0._r8
      du_fix(:)       = 0._r8
      du_fix_s(:)     = 0._r8
      dufix_xy(:,:,:) = 0._r8
   end if

   if (am_diag) then
      ucxy = 0._r8
      vcxy = 0._r8

!$omp parallel do private(i,j,k)
      ! store old winds to get total increments
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               duxy(i,j,k)=uxy(i,j,k)
               dvxy(i,j,k)=vxy(i,j,k)
            enddo
         enddo
      enddo
   end if

   ! Begin vertical re-mapping sub-cycle loop
   do iv = 1, nv

      if (iv == nv) convt_local = convt

      ! Transpose XY arrays to YZ
      call t_barrierf('sync_xy_to_yz_1', grid%commdyn)
      call t_startf ('xy_to_yz')

      if (grid%iam .lt. grid%npes_xy) then

         if (grid%twod_decomp .eq. 1) then

#if defined( SPMD )


!$omp parallel do private(i,j,k)
            ! Embed psxy and phisxy in 3D array since transpose machinery cannot handle 2D arrays
            do k = 1, npr_z
               do j = jfirstxy, jlastxy
                  do i = ifirstxy, ilastxy
                     psxy3(i,j,k)   = psxy(i,j)
                     phisxy3(i,j,k) = phisxy(i,j)
                  end do
               end do
            end do

            if (grid%modc_onetwo .eq. 1) then
               call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                               modc=grid%modc_dynrun )
               call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                               modc=grid%modc_dynrun )

               call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, phisxy3, phis,                 &
                               modc=grid%modc_dynrun )
               call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, phisxy3, phis,                 &
                               modc=grid%modc_dynrun )
            else
               call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                               phisxy3, phis,                                             &
                               modc=grid%modc_dynrun )
               call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                               grid%xy2d_to_yz2d%RecvDesc, psxy3, ps,                     &
                               phisxy3, phis,                                             &
                               modc=grid%modc_dynrun )
            end if

            ! if OFFLINE_DYN is defined, u and v are filled at this point

#if defined( OFFLINE_DYN )
            call mp_sendirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                             grid%uxy_to_u%RecvDesc, uxy, u_tmp,                         &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                             grid%uxy_to_u%RecvDesc, uxy, u_tmp,                         &
                             modc=grid%modc_dynrun )

            call mp_sendirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                             grid%vxy_to_v%RecvDesc, vxy, v_tmp,                         &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                             grid%vxy_to_v%RecvDesc, vxy, v_tmp,                         &
                             modc=grid%modc_dynrun )

!$omp parallel do private(i,j,k)
            do k = kfirst, klast
               do j = jfirst, jlast
                  do i = 1, im
                     u(i,j,k) = (1._r8-met_rlx(k))*u_tmp(i,j,k) + met_rlx(k)*u(i,j,k)
                     v(i,j,k) = (1._r8-met_rlx(k))*v_tmp(i,j,k) + met_rlx(k)*v(i,j,k)
                  end do
               end do
            end do
#else
            call mp_sendirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                             grid%uxy_to_u%RecvDesc, uxy, u,                             &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%uxy_to_u%SendDesc,                        &
                             grid%uxy_to_u%RecvDesc, uxy, u,                             &
                             modc=grid%modc_dynrun )

            call mp_sendirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                             grid%vxy_to_v%RecvDesc, vxy, v,                             &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%vxy_to_v%SendDesc,                        &
                             grid%vxy_to_v%RecvDesc, vxy, v,                             &
                             modc=grid%modc_dynrun )
#endif

            call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                      &
                             grid%pexy_to_pe%RecvDesc, pexy, pe,                         &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                      &
                             grid%pexy_to_pe%RecvDesc, pexy, pe,                         &
                             modc=grid%modc_dynrun )

            call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                             grid%ijk_xy_to_yz%RecvDesc, delpxy, delp,                   &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                             grid%ijk_xy_to_yz%RecvDesc, delpxy, delp,                   &
                             modc=grid%modc_dynrun )

            call mp_sendirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,                     &
                             grid%pkxy_to_pkc%RecvDesc, pkxy, pk,                        &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%pkxy_to_pkc%SendDesc,                     &
                             grid%pkxy_to_pkc%RecvDesc, pkxy, pk,                        &
                             modc=grid%modc_dynrun )

            call mp_sendirr( grid%commxy, grid%ptxy_to_pt%SendDesc,                      &
                             grid%ptxy_to_pt%RecvDesc, ptxy, pt,                         &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%ptxy_to_pt%SendDesc,                      &
                             grid%ptxy_to_pt%RecvDesc, ptxy, pt,                         &
                             modc=grid%modc_dynrun )
            if (high_alt) then
               call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                                grid%ijk_xy_to_yz%RecvDesc, cp3v, cp3vc,                    &
                                modc=grid%modc_dynrun )
               call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                                grid%ijk_xy_to_yz%RecvDesc, cp3v, cp3vc,                    &
                                modc=grid%modc_dynrun )
               call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                                grid%ijk_xy_to_yz%RecvDesc, cap3v, cap3vc,                  &
                                modc=grid%modc_dynrun )
               call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                                grid%ijk_xy_to_yz%RecvDesc, cap3v, cap3vc,                  &
                                modc=grid%modc_dynrun )
            endif

            if (modc_tracers .eq. 0) then
               do mq = 1, ntotq
                  q3xypt => tracer(:,:,:,mq)
                  q3yzpt => q_internal(:,:,:,mq)
                  call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                                   grid%ijk_xy_to_yz%RecvDesc, q3xypt, q3yzpt,               &
                                   modc=grid%modc_dynrun )
                  call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                  &
                                   grid%ijk_xy_to_yz%RecvDesc, q3xypt, q3yzpt,               &
                                   modc=grid%modc_dynrun )
               end do
            else
               do mq = 1, ntotq, modc_tracers
                  mlast = min(mq+modc_tracers-1,ntotq)
                  call mp_sendtrirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                      &
                                     grid%ijk_xy_to_yz%RecvDesc, tracer, q_internal, mq, mlast, ntotq,    &
                                     grid%ifirstxy, grid%ilastxy, grid%jfirstxy, grid%jlastxy,     &
                                     1, grid%km,                                                   &
                                     1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, grid%klast, &
                                     modc=grid%modc_tracer )
                  call mp_recvtrirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                      &
                                     grid%ijk_xy_to_yz%RecvDesc, tracer, q_internal, mq, mlast, ntotq,    &
                                     grid%ifirstxy, grid%ilastxy, grid%jfirstxy, grid%jlastxy,     &
                                     1, grid%km,                                                   &
                                     1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, grid%klast, &
                                     modc=grid%modc_tracer )
               end do
            end if

#else
            write(iulog,*)'DYN_COMP:dyn_run -- SPMD must be defined for 2D decomp -- returning'
            rc = DYN_RUN_MUST_BE_2D_DECOMP
            return    ! Not possible to have 2D decomposition with SPMD undefined
#endif
         else ! if not twod_decomp

            do j = jfirst, jlast
               do i = 1, im
                  ps(i,j)   = psxy(i,j)
                  phis(i,j) = phisxy(i,j)
               end do
            end do

!$omp parallel do private(i,j,k)
            do j = jfirst, jlast
               do k = 1, km+1
                  do i = 1, im
                     pe(i,k,j) = pexy(i,k,j)
                  end do
               end do
            end do

!$omp parallel do private(i,j,k)
            do k = 1, km+1
               do j = jfirst, jlast
                  do i = 1, im
                     pk(i,j,k) = pkxy(i,j,k)
                  end do
               end do
            end do

!$omp parallel do private(i,j,k)
            do k = 1, km
               do j = jfirst, jlast
                  do i = 1, im
#if defined( OFFLINE_DYN )
                     u(i,j,k)    = (1._r8-met_rlx(k))*uxy(i,j,k) + met_rlx(k)*u(i,j,k)
                     v(i,j,k)    = (1._r8-met_rlx(k))*vxy(i,j,k) + met_rlx(k)*v(i,j,k)
#else
                     u(i,j,k)    = uxy(i,j,k)
                     v(i,j,k)    = vxy(i,j,k)
#endif
                     delp(i,j,k) = delpxy(i,j,k)
                     pt(i,j,k)   = ptxy(i,j,k)
                  end do
               end do
            end do
            if (high_alt) then
!$omp parallel do private(i,j,k)
               do k = 1, km
                  do j = jfirst, jlast
                     do i = 1, im
                        cp3vc(i,j,k) = cp3v(i,j,k)
                        cap3vc(i,j,k) = cap3v(i,j,k)
                     end do
                  end do
               end do
            endif

            do mq = 1, ntotq

               ! For now just copy in the contents of tracer; later, use pointers
               ! TODO:  q_internal(mq) => tracer(mq)    ! Make sure not to allocate q_internal in this case

               q_internal(1:im,jfirst:jlast,kfirst:klast,mq) = &
                   tracer(1:im,jfirst:jlast,kfirst:klast,mq)
            end do

         end if   ! (grid%twod_decomp .eq. 1)

      end if      ! (grid%iam .lt. grid%npes_xy)

#if defined(SPMD)

      ! Send tracers to auxiliary processes when overlapping
      if (ct_overlap .gt. 0 .and. n2 .gt. 1 .and. kaux .eq. 0) then
         do iq = 1, ntotq
            call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tractag+iq-1, commnyz, ctreqs(5+iq,1))
         end do
      end if

      ! Send tracers to auxiliary processes when decomposing
      if (trac_decomp .gt. 1 .and. kaux .eq. 0) then
         do kt = 2, trac_decomp
            do iq = ktloa(kt), kthia(kt)
               tagu = tractag+iq-1 + (kt-2)*(ntotq+ntg0)
               call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremotea(kt-1), tagu, commnyz, ctreqs(5+iq,kt-1))
            end do
         end do
      end if
#endif

      call t_stopf  ('xy_to_yz')

      omgaxy(:,:,:) = 0._r8

      if (am_fixer .or. am_diag) then
         du_fix_s (:)  = 0._r8
         uc_s (:,:,:)  = 0._r8
         vc_s (:,:,:)  = 0._r8
      endif

      ! Begin tracer sub-cycle loop
      do n = 1, n2

         if (ntotq > 0) then

            call t_barrierf('sync_small_ts_init', grid%commdyn)
            call t_startf('small_ts_init')

!$omp parallel do private(i, j, k)
            do k = kfirst, klast
               do j = jfirst, jlast
                  do i = 1, im
                     ! Save initial delp field before the small-time-step
                     ! Initialize the CFL number accumulators: (cx, cy)
                     ! Initialize total mass fluxes: (mfx, mfy)
                     dp0(i,j,k) = delp(i,j,k)
                     cx(i,j,k) = 0._r8
                     cy(i,j,k) = 0._r8
                     mfx(i,j,k) = 0._r8
                     mfy(i,j,k) = 0._r8
                  end do
               end do
            end do

#if defined( SPMD )
            if (grid%iam .lt. grid%npes_yz) then
               call mp_send4d_ns( grid%commyz, im, jm, km,                     &
                                  1, jfirst, jlast, kfirst, klast, 1, 0, dp0 )
               call mp_recv4d_ns( grid%commyz, im, jm, km,                     &
                                  1, jfirst, jlast, kfirst, klast, 1, 0, dp0 )
            end if
#endif

            call t_stopf  ('small_ts_init')

         end if  ! (ntotq > 0)

#if defined(SPMD)

         ! Send dp0 to auxiliary processes when overlapping or tracer decomposition
         if (kaux .eq. 0) then
            if (ct_overlap .gt. 0 .and. n .lt. n2) then
               call mpiisend(dp0, ndp0, mpir8, iremote, dp0tag, commnyz, ctreqs(1,1))
            end if

            if (trac_decomp .gt. 1) then
               do kt = 2, trac_decomp
                  tagu = dp0tag + (kt-2)*(ntotq+ntg0)
                  call mpiisend(dp0, ndp0, mpir8, iremotea(kt-1), tagu, commnyz, ctreqs(1,kt-1))
               end do
            end if
         end if
#endif

         ! Begin dynamics sub-cycle loop
         do it = 1, nsplit

            if (it == nsplit .and. n == n2) then
               ipe = 1                     ! end of cd_core; output pexy for te_map
            else if (it == 1 .and. n == 1) then
               ipe = -1                    ! start of cd_core
            else
               ipe = 0
            end if

            ! determine whether this is the second to last call to cd_core or not
            cd_penul = .false.
            if ( nsplit > 1 ) then
               if ( (n == n2) .and. (it == nsplit-1) ) cd_penul = .true.
            else if ( n2 > 1 ) then
               if ( n == n2-1 ) cd_penul = .true.
            end if

            if (cd_penul) then
               if (ipe == -1) then
                  ipe = -2   ! second to last is also the first
               else
                  ipe = 2
               end if
            end if

            ! s_trac is true if cd_core is to post sends for ct_overlap or trac_decomp
            ! such sends are posted during last inner cd_core subcycle
            s_trac = ((ct_overlap .gt. 0 .and. it .eq. nsplit .and. n .lt. n2) .or.     &
                      (trac_decomp .gt. 1 .and. it .eq. nsplit))


            if ((it == nsplit) .and. (n == n2) .and. (iv == nv)) then
!$omp parallel do private(j)
               do j = jfirstxy, jlastxy
                  pexy_om(ifirstxy:ilastxy,1:km+1,j) = pexy(ifirstxy:ilastxy,1:km+1,j)
               end do
            end if

            ! Call the Lagrangian dynamical core using small tme step

            call t_barrierf('sync_cd_core', grid%commdyn)
            call t_startf ('cd_core')

            if (grid%iam .lt. grid%npes_yz) then
               am_fix_out = am_fixer .or. am_diag
               call cd_core(grid,   nx,     u,   v,   pt,                    &
                            delp,   pe,     pk,  nsplit,  dt,                &
                            ptop,   umax,   pi, ae,                          &
                            cp3vc,  cap3vc, cp3v, cap3v,                     &
                            icd,    jcd, iord, jord,   ipe,                  &
                            div24del2flag, del2coef,                         &
                            om,     phis,     cx  ,  cy, mfx, mfy,           &
                            delpf, uc, vc, pkz, dpt, worka,                  &
                            dwz, pkc, wz,  phisxy, ptxy, pkxy,               &
                            pexy, pkcc, wzc, wzxy, delpxy,                   &
                            pkkp, wzkp, cx_om, cy_om, filtcw, s_trac,        &
                            naux, ncx, ncy, nmfx, nmfy, iremotea,            &
                            cxtaga, cytaga, mfxtaga, mfytaga, cdcreqs(1,1),  &
                            cdcreqs(1,2), cdcreqs(1,3), cdcreqs(1,4),        &
                            kmtp, am_correction, am_geom_crrct, am_fix_out,  &
                            dod, don, high_order_top)

               ctreqs(2,:) = cdcreqs(:,1)
               ctreqs(3,:) = cdcreqs(:,2)
               ctreqs(4,:) = cdcreqs(:,3)
               ctreqs(5,:) = cdcreqs(:,4)
            end if !  (grid%iam .lt. grid%npes_yz)

            call t_stopf  ('cd_core')

            ! AM fixer
            if (am_fixer.or.am_diag) then

               call t_barrierf('sync_lfix', grid%commdyn)
               call t_startf ('lfix')
               if (grid%iam .lt. grid%npes_yz) then

                  ! option for pressure tapering on AM fixer
                  if (am_fix_taper) then
                     zpkck(:,:)=0._r8
!$omp parallel do private(j, k)
                     do k=kfirst,klast
                        do j = jfirst, jlast
                           zpkck(j,k)=0.25_r8*sum(pkc(:,j,k))*grid%cose(j)
                        enddo
                     enddo
                     do k=kfirst,klast
                        call par_vecsum(jm, jfirst, jlast, zpkck(1:jm,k), me0, grid%comm_y, grid%npr_y)
                        avgpk(k)=me0/im/sum(grid%cose)
                        taper(k)=.5_r8*(1._r8+(1._r8-(ptapk/avgpk(k))**xdlt2)/(1._r8+(ptapk/avgpk(k))**xdlt2))
                     enddo
                  else
                     do k=kfirst,klast
                        taper(k)=1._r8
                     enddo
                  endif

                  ! always exclude fixer at top levels if top is not high order
                  if (.not.high_order_top) then
                     taper(1:kmtp)=0._r8
                  endif

                  do k = kfirst, klast
                     call par_vecsum(jm, jfirst, jlast, don(1:jm,k), am1, grid%comm_y, grid%npr_y)
                     dons(k) = am1
                  end do

                  do k = kfirst, klast
                     call par_vecsum(jm, jfirst, jlast, dod(1:jm,k), me0, grid%comm_y, grid%npr_y)
                     dods(k) = me0
                  end do

                  if (am_fix_lbl) then
!$omp parallel do private(i, j, k)
                     do k = kfirst, klast
                        do j = jfirst, jlast
                           do i = 1, im
                              u(i,j,k) = u(i,j,k) - dons(k)/dods(k)*grid%cose(j) * taper(k)
                           end do
                        end do
                     end do
                  endif

                  ! diagnose du_fix
                  if (am_fix_lbl) then ! output applied increment (tapered)
!$omp parallel do private(k)
                     do k = kfirst, klast
                        du_fix_s(k)=du_fix_s(k)-dons(k)/dods(k)*taper(k)
                     end do
                  elseif(am_diag) then ! output diagnosed increment (not tapered)
!$omp parallel do private(k)
                     do k = kfirst, klast
                        du_fix_s(k)=du_fix_s(k)-dons(k)/dods(k)
                     end do
                  endif

!$omp parallel do private(j, k)
                  do k=kfirst,klast
                     do j = jfirst, jlast
                        don(j,k)=don(j,k)*taper(k)
                        dod(j,k)=dod(j,k)*taper(k)
                     enddo
                  enddo
                  tmpsum(1,1) = SUM(don)
                  tmpsum(1,2) = SUM(dod)
                  call shr_reprosum_calc(tmpsum, tmpresult, 1, 1, 2, commid=grid%commyz)
                  am1 = tmpresult(1)
                  me0 = max(tmpresult(2),tiny)

                  if (am_fixer.and.(.not.am_fix_lbl)) then
!$omp parallel do private(i, j, k)
                     do k = kfirst, klast
                        do j = jfirst, jlast
                           do i = 1, im
                              u(i,j,k) = u(i,j,k) - am1/me0*grid%cose(j) *taper(k)
                           end do
                        end do
                     end do

!$omp parallel do private(k)
                     do k = kfirst, klast
                        du_fix_s(k)=du_fix_s(k)-am1/me0*taper(k)
                     end do
                  end if  ! (am_fix_lbl)

                  du_fix_g =du_fix_g -am1/me0
                  if (masterproc) then
                     if ((it == nsplit) .and. (n == n2) .and. (iv == nv)) then
                        write(iulog,'(1x,a21,1x,1x,e25.17)') "AM GLOBAL FIXER: ", du_fix_g
                     endif
                  endif
                  ! the following call is blocking, but probably cheaper than 3D transposition for du_fix
                  if ((it == nsplit) .and. (n == n2)) then
                     call par_vecsum(km, kfirst, klast, du_fix_s, tmp, grid%comm_z, grid%npr_z, return_sum_in=.true.)
                  endif
               end if     ! (grid%iam .lt. grid%npes_yz)
               call t_stopf  ('lfix')

!$omp parallel do private(i,j,k)
               do k=kfirst,klast
                  do j = jfirst, jlast
                     do i=1,im
                        uc_s(i,j,k)=uc_s(i,j,k)+uc(i,j,k)
                        vc_s(i,j,k)=vc_s(i,j,k)+vc(i,j,k)
                     enddo
                  enddo
               enddo

            end if    ! (am_fixer.or.am_diag)

            if ((it == nsplit) .and. (n == n2) .and. (iv == nv)) then
!$omp  parallel do     &
!$omp  default(shared) &
!$omp  private(i,j,k)
               do j = jfirstxy, jlastxy
                  do k = 1, km
                     do i = ifirstxy, ilastxy
                        omgaxy(i,k,j) = omgaxy(i,k,j) + 0.5_r8*(pexy(i,k,j) + pexy(i,k+1,j) - &
                                        pexy_om(i,k,j) - pexy_om(i,k+1,j))/dt
                     end do
                  end do
                  do k = 1, km+1
                     do i = ifirstxy, ilastxy
                        pexy_om(i,k,j) = 0.5_r8*(pexy_om(i,k,j) + pexy(i,k,j))
                     end do
                  end do
               end do

               !-----------------------------------------------------
               ! Add the v*grad(p) term to omega (dp/dt) for physics
               !-----------------------------------------------------
               call t_startf ('vdot_gradp')
               if (grid%iam .lt. grid%npes_xy) then
                  call compute_vdot_gradp( grid, dt, dt/dt, cx_om, cy_om, pexy_om, omgaxy )
               end if
               call t_stopf  ('vdot_gradp')

            end if

         end do  !  it = 1, nsplit - dynamics sub-cycle loop

         if (ntotq .ne. 0) then

#if ( defined OFFLINE_DYN )
            if (met_fix_mass) then
               ps_mod(:,:) = ps(:,:)
               ! get the observed PS interpolated to current substep
               call get_met_fields(grid, ps_obs, n2, n)

               ! adjust mass fluxes and edge pressures to be consistent with observed PS
               call adjust_press(grid, ps_mod, ps_obs, mfx, mfy, pexy)

               if (high_alt) then
!$omp parallel do private(i,j,k)
                  do k=2,km
                     do j=jfirstxy,jlastxy
                        do i=ifirstxy,ilastxy
                           cap3vi(i,j,k) = 0.5_r8*(cap3v(i,j,k-1)+cap3v(i,j,k))
                        enddo
                     enddo
                  enddo
                  cap3vi(:,:,1) = 1.5_r8 * cap3v(:,:,1) - 0.5_r8 * cap3v(:,:,2)
                  cap3vi(:,:,km+1) = 1.5_r8 * cap3v(:,:,km) - 0.5_r8 * cap3v(:,:,km-1)
              endif
!$omp parallel do private(i,j,k)
               ! make pkxy consistent with the adjusted pexy
               do i = ifirstxy, ilastxy
                  do j = jfirstxy, jlastxy
                     do k = 1, km+1
                        pkxy(i,j,k) = pexy(i,k,j)**cap3vi(i,j,k)
                     end do
                  end do
               end do

!$omp parallel do private(i,j,k)
               ! adjust courant numbers to be consistent with the adjusted mass fluxes
               do i = 1, im
                  do j = jfirst, jlast
                     do k = kfirst, klast
                        if (i .ne. 1) cx(i,j,k) = mfx(i,j,k)/(0.5_r8*(dp0(i-1,j,k)+dp0(i,j,k)))
                        if (i .eq. 1) cx(i,j,k) = mfx(i,j,k)/(0.5_r8*(dp0(1,j,k)+dp0(im,j,k)))
                     end do
                  end do
               end do

!$omp parallel do private(i,j,k)
               do i = 1, im
                  do j = jfirst, jlast
                     do k = kfirst, klast
                        if ((j .gt. 1) .and. (j .lt. jm)) cy(i,j,k) = &
                           mfy(i,j,k)/(0.5_r8*(dp0(i,j-1,k)+dp0(i,j,k)))/grid%cose(j)
                     end do
                  end do
               end do
            end if
#endif

            ! WS 2006-12-04 : this seems like the safest place to preprocess and
            !                 transpose the C-grid mass-flux and later the
            !                 Courant numbers for potential output

            ! Horizontal mass fluxes

            if (grid%iam .lt. grid%npes_xy) then

               if (grid%twod_decomp .eq. 1) then
#if defined( SPMD )
!$omp parallel do private(i,j,k)
                  do k = kfirst, klast
                     do j = jfirst, jlast
                        do i = 1, im
                           worka(i,j,k) = mfx(i,j,k)*(ae*grid%dp)*(grid%dl*ae*grid%cosp(j))/(ndt) ! Pa m^2/s
                           workb(i,j,k) = mfy(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                        end do
                     end do
                  end do
                  if (grid%modc_onetwo .eq. 1) then
                     call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                      modc=grid%modc_dynrun )
                     call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                      modc=grid%modc_dynrun )

                     call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, workb, mfyxy,                 &
                                      modc=grid%modc_dynrun )
                     call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, workb, mfyxy,                 &
                                      modc=grid%modc_dynrun )
                  else
                     call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                      workb, mfyxy,                                             &
                                      modc=grid%modc_dynrun )
                     call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                                      grid%ijk_yz_to_xy%RecvDesc, worka, mfxxy,                 &
                                      workb, mfyxy,                                             &
                                      modc=grid%modc_dynrun )
                  end if

#else
                  write(iulog,*)'DYN_COMP:dyn_run -- SPMD must be defined for 2D decomp -- returning'
                  rc = DYN_RUN_MUST_BE_2D_DECOMP
                  return    ! Not possible to have 2D decomposition with SPMD undefined
#endif
               else   ! if not twod_decomp   (1D or sequential)
!$omp parallel do private(i,j,k)
                  do k = kfirst, klast
                     do j = jfirst, jlast
                        do i = 1, im
                           mfxxy(i,j,k) = mfx(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                           mfyxy(i,j,k) = mfy(i,j,k)*(grid%dl*ae*grid%cosp(j))*(ae*grid%dp)/(ndt*grid%cose(j)) ! Pa m^2 / s
                        end do
                     end do
                  end do

               end if

            end if  !  (grid%iam .lt. grid%npes_xy)


            ! Perform large-tme-step scalar transport using the accumulated CFL and
            ! mass fluxes

            call t_barrierf('sync_trac2d', grid%commdyn)
            call t_startf ('trac2d')

            ! Overlap trac2d with subsequent cd_core set, or decompose over tracers

            if ((ct_overlap .gt. 0 .and. n .lt. n2 .and. kaux .lt. 2) .or.   &
               (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)) then

               if (kaux .eq. 0) then

                  ! Primary process

                  ! Send data to auxiliary yz decomposition
                  ! Communicate tracers on first subcycle only
                  ! Also post receive of new tracer values from aux processes

#if defined(SPMD)
                  if (n .eq. 1) then

                     ! Block on send of tracers to aux
                     if (ct_overlap .gt. 0) then
                        do iq = 1, ntotq
                           call mpiwait(ctreqs(5+iq,1), ctstats(1,5+iq,1))
                        end do
                     end if

                     if (trac_decomp .gt. 1) then
                        do kt = 2, trac_decomp
                           do iq = ktloa(kt), kthia(kt)
                              call mpiwait(ctreqs(5+iq,kt-1), ctstats(1,5+iq,kt-1))
                           enddo
                        enddo
                     endif

                     ! Post receive for updated tracers from aux
                     if (ct_overlap .gt. 0) then
                        do iq = 1, ntotq
                           call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremote,    &
                               tractag+iq-1, commnyz, ctreq(iq,1))
                        end do
                     end if

                     if (trac_decomp .gt. 1) then
                        do kt = 2, trac_decomp
                           do iq = ktloa(kt), kthia(kt)
                              tagu = tractag+iq-1 + (kt-2)*(ntotq+ntg0)
                              call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremotea(kt-1),   &
                                            tagu, commnyz, ctreq(iq,kt-1))
                           end do
                        end do
                     end if
                  end if  !  (n .eq. 1)

                  if (ct_overlap .gt. 0) then

                     ! Block on send of dp0 to aux
                     call mpiwait(ctreqs(1,1), ctstats(1,1,1))
                     ! Block on sends from cd_core to aux
                     call mpiwait(ctreqs(2,1), ctstats(1,2,1))
                     call mpiwait(ctreqs(3,1), ctstats(1,3,1))
                     call mpiwait(ctreqs(4,1), ctstats(1,4,1))
                     call mpiwait(ctreqs(5,1), ctstats(1,5,1))
                  endif

                  if (trac_decomp .gt. 1) then

                     do kt = 2, trac_decomp
                        ! Block on send of dp0 to aux
                        call mpiwait(ctreqs(1,kt-1), ctstats(1,1,kt-1))
                        ! Block on sends from cd_core to aux
                        call mpiwait(ctreqs(2,kt-1), ctstats(1,2,kt-1))
                        call mpiwait(ctreqs(3,kt-1), ctstats(1,3,kt-1))
                        call mpiwait(ctreqs(4,kt-1), ctstats(1,4,kt-1))
                        call mpiwait(ctreqs(5,kt-1), ctstats(1,5,kt-1))
                     end do
                  end if
#endif

               else

                  ! Auxiliary process

                  ! Temporarily set subdomain limits and process index in auxiliary process equal
                  ! to those of antecedent
                  jfirst = jfirstct
                  jlast  = jlastct
                  kfirst = kfirstct
                  klast  = klastct
                  grid%jfirst = jfirstct
                  grid%jlast  = jlastct
                  grid%kfirst = kfirstct
                  grid%klast  = klastct
                  ! Translate process index to frame of auxiliary yz decomposition for use with auxiliary
                  !    communication in trac2d
                  grid%iam = iremote

#if defined(SPMD)
                  ! Receive data from primary yz decomposition
                  ! Include tracers first subcycle only

                  if (n .eq. 1) then
                     do iq = ktlo, kthi
                        tagu = tractag+iq-1 + (kaux-1)*(ntotq+ntg0)
                        call mpiirecv(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tagu, commnyz, ctreq(5+iq,1))
                        call mpiwait(ctreq(5+iq,1), ctstat(1,5+iq,1))
                     end do
                  end if
                  tagu = dp0tag + (kaux-1)*(ntotq+ntg0)
                  call mpiirecv(dp0, ndp0, mpir8, iremote, tagu, commnyz, ctreq(1,1))
                  tagu = cxtag + (kaux-1)*(ntotq+ntg0)
                  call mpiirecv(cx, ncx, mpir8, iremote, tagu, commnyz, ctreq(2,1))
                  tagu = cytag + (kaux-1)*(ntotq+ntg0)
                  call mpiirecv(cy, ncy, mpir8, iremote, tagu, commnyz, ctreq(3,1))
                  tagu = mfxtag + (kaux-1)*(ntotq+ntg0)
                  call mpiirecv(mfx, nmfx, mpir8, iremote, tagu, commnyz, ctreq(4,1))
                  tagu = mfytag + (kaux-1)*(ntotq+ntg0)
                  call mpiirecv(mfy, nmfy, mpir8, iremote, tagu, commnyz, ctreq(5,1))
                  call mpiwait(ctreq(1,1), ctstat(1,1,1))
                  call mpiwait(ctreq(2,1), ctstat(1,2,1))
                  call mpiwait(ctreq(3,1), ctstat(1,3,1))
                  call mpiwait(ctreq(4,1), ctstat(1,4,1))
                  call mpiwait(ctreq(5,1), ctstat(1,5,1))
#endif

               end if  !  (kaux .eq. 0)

            else

#if defined(SPMD)
               ! Block on receive of updated tracers from aux (last subcycle)
               if (ct_overlap .gt. 0 .and. n .eq. n2 .and. n2 .gt. 1 .and. kaux .eq. 0) then
                  do iq = 1, ntotq
                     call mpiwait(ctreq(iq,1), ctstat(1,iq,1))
                  end do
               end if
#endif

            end if  !  (ct_overlap .gt. 0 .and. n .lt. n2 .and. kaux .lt. 2)
                    ! or (trac_decomp .gt. 1 .and. kaux .lt. trac_decomp)

            ! Call tracer advection

            c_dotrac = ct_overlap .gt. 0 .and.    &
                      ((n .lt. n2 .and. kaux .eq. 1) .or. (n .eq. n2 .and. kaux .eq. 0))
            t_dotrac = ct_overlap .eq. 0 .and. kaux .lt. trac_decomp
            high_alt1: if (high_alt) then
               !
               ! phl: overwrite last tracer with kappa
               !
               !$omp parallel do private(i,j,k)
               do k=grid%kfirst,grid%klast
                  do j=grid%jfirst,grid%jlast
                     do i=1,grid%im
                        q_internal(i,j,k,ntotq) = cap3vc(i,j,k)
                     end do
                  end do
               end do
            endif high_alt1
            if (c_dotrac .or. t_dotrac) then
               call trac2d( grid, dp0(:,jfirst:jlast,:),    q_internal,         &
                           cx,    cy,     mfx,    mfy,    iord,   jord,        &
                           fill,  ktlo,   kthi,   workb,  worka  )
            endif

#if defined(SPMD)
            ! Return data to primary yz decomposition
            ! For overlap, next-to-last subcycle only; for tracer decomp, last subcycle only
            if (ct_aux .and. ((ct_overlap .gt. 0 .and. n .eq. n2-1) .or.   &
               (trac_decomp .gt. 1 .and. n .eq. n2)))  then
               do iq = ktlo, kthi
                  tagu = tractag+iq-1 + (kaux-1)*(ntotq+ntg0)
                  call mpiisend(q_internal(:,:,:,iq), ntrac, mpir8, iremote, tagu, commnyz, ctreqs(5+iq,1))
                  call mpiwait(ctreqs(5+iq,1), ctstats(1,5+iq,1))
               end do
            end if
#endif

#if defined(SPMD)
            ! For tracer decomposition, block on receive of updated tracers from aux (last subcycle)
            if (trac_decomp .gt. 1 .and. n .eq. n2 .and. kaux .eq. 0) then
               do kt = 2, trac_decomp
                  do iq = ktloa(kt), kthia(kt)
                     call mpiwait(ctreq(iq,kt-1), ctstat(1,iq,kt-1))
                  end do
               end do
            end if
#endif

            ! Restore subdomain limits and process index in auxiliary process
            if (ct_aux) then
               jfirst = jkstore(1)
               jlast  = jkstore(2)
               kfirst = jkstore(3)
               klast  = jkstore(4)
               grid%jfirst = jkstore(1)
               grid%jlast  = jkstore(2)
               grid%kfirst = jkstore(3)
               grid%klast  = jkstore(4)
               grid%iam = iamlocal
            end if

            ! NOTE: for cd_core / trac2d overlap, tracer data is returned to primary processes
            !   prior to n=n2 call to trac2d

            call t_stopf('trac2d')

            trans_pexy: if (met_fix_mass .or. high_alt) then

               if (grid%twod_decomp .eq. 1) then
                 if (grid%iam .lt. grid%npes_yz) then
#if defined( SPMD )
                  call mp_sendirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                                  grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                                  modc=grid%modc_dynrun )
                  call mp_recvirr( grid%commxy, grid%pexy_to_pe%SendDesc,                    &
                                  grid%pexy_to_pe%RecvDesc, pexy, pe,                       &
                                  modc=grid%modc_dynrun )
#endif
                 endif
               else
!$omp parallel do private(i,j,k)
                  do j = jfirst, jlast
                     do k = kfirst, klast+1
                        do i = 1, im
                           pe(i,k,j) = pexy(i,k,j)
                        end do
                     end do
                  end do
               end if

#if ( defined OFFLINE_DYN )
               do j = jfirst,jlast
                  if (klast .eq. km) ps_mod(:,j) = pe(:,km+1,j)
               end do
#endif
            end if trans_pexy

            high_alt2: if (high_alt) then

               !+hi-waccm: perform potential temperature correction:
               !           1. Update kappa according to new major species
               !           2. calculate the difference between kappa from step 1 and kappa from advection
               !           3. calculate ln(p0/p) from the most recent p
               !           4. update pt, then transpose to ptxy.

               ! Since rairv is not defined on yz decomp, can retrieve it using cp3vc and cap3vc when needed
               ! Also will check if mbarv is needed somewhere in dynamics.
               ! These updates of cp3vc, cap3vc etc are currently not passed back to physics.
               ! This update is put here, after the transpose of pexy to pe, since we need pe (on yz decomp).

               call physconst_calc_kappav(1,im,jfirst,jlast,kfirst,klast, grid%ntotq, q_internal, cap3vc )

!$omp parallel do private(i,j,k)
               do k = kfirst,klast
                  do j = jfirst,jlast
                     do i = 1,im
                        pt(i,j,k) = pt(i,j,k) * (1._r8 - &
                             .5_r8*(log(pe(i,k,j))+log(pe(i,k+1,j))) * &
                             (cap3vc(i,j,k)-q_internal(i,j,k,ntotq)))
                     enddo
                  enddo
               enddo

            endif high_alt2
         end if  !          if (ntotq .ne. 0) then

      end do !   do n=1, n2 - tracer sub-cycle loop

      call t_barrierf('sync_yz_to_xy_1', grid%commdyn)

      if (grid%iam .lt. grid%npes_xy) then

         if (grid%twod_decomp .eq. 1) then

            ! Transpose ps, u, v, and tracer from yz to xy decomposition
            !
            ! Note: pt, pe and pk will have already been transposed through
            ! call to geopk in cd_core. geopk does not actually require
            ! secondary xy decomposition; direct 16-byte technique works just
            ! as well, perhaps better. However, transpose method is used on last
            ! call to avoid having to compute these three transposes now.

#if defined (SPMD)

            call t_startf ('yz_to_xy_psuv')

            ! Transpose ps
            ! Embed in 3D array since transpose machinery cannot handle 2D arrays

            call mp_sendirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                             grid%yz2d_to_xy2d%RecvDesc, ps, psxy3,                     &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                             grid%yz2d_to_xy2d%RecvDesc, ps, psxy3,                     &
                             modc=grid%modc_dynrun )

!$omp parallel do private(i,j)
            do j = jfirstxy, jlastxy
               do i = ifirstxy, ilastxy
                  psxy(i,j) = psxy3(i,j,1)
               end do
            end do

            ! Transpose u
            call mp_sendirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                             grid%u_to_uxy%RecvDesc, u, uxy,                            &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                             grid%u_to_uxy%RecvDesc, u, uxy,                            &
                             modc=grid%modc_dynrun )

            ! Transpose v
            call mp_sendirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                             grid%v_to_vxy%RecvDesc, v, vxy,                            &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                             grid%v_to_vxy%RecvDesc, v, vxy,                            &
                             modc=grid%modc_dynrun )


            if (am_fixer.or.am_diag) then

               ! Transpose uc_s
               call mp_sendirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                                grid%u_to_uxy%RecvDesc, uc_s, uc_i,                        &
                                modc=grid%modc_dynrun )
               call mp_recvirr( grid%commxy, grid%u_to_uxy%SendDesc,                       &
                                grid%u_to_uxy%RecvDesc, uc_s, uc_i,                        &
                                modc=grid%modc_dynrun )

               ! Transpose vc_s
               call mp_sendirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                                grid%v_to_vxy%RecvDesc, vc_s, vc_i,                        &
                                modc=grid%modc_dynrun )
               call mp_recvirr( grid%commxy, grid%v_to_vxy%SendDesc,                       &
                                grid%v_to_vxy%RecvDesc, vc_s, vc_i,                        &
                                modc=grid%modc_dynrun )
            end if

            call t_stopf  ('yz_to_xy_psuv')

            call t_startf ('yz_to_xy_q')

            if (modc_tracers .eq. 0) then
               do mq = 1, ntotq


                  ! Transpose
                  q3yzpt => q_internal(:,:,:,mq)
                  q3xypt => dyn_out%tracer(:,:,:,mq)
                  call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                                   grid%ijk_yz_to_xy%RecvDesc, q3yzpt, q3xypt,                 &
                                   modc=grid%modc_dynrun )
                  call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                                   grid%ijk_yz_to_xy%RecvDesc, q3yzpt, q3xypt,                 &
                                   modc=grid%modc_dynrun )
               end do
            else
               do mq = 1, ntotq, modc_tracers
                  mlast = min(mq+modc_tracers-1,ntotq)
                  call mp_sendtrirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                     &
                                     grid%ijk_yz_to_xy%RecvDesc, q_internal, dyn_out%tracer,      &
                                     mq, mlast, ntotq, 1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, &
                                     grid%klast, grid%ifirstxy, grid%ilastxy, grid%jfirstxy,      &
                                     grid%jlastxy, 1, grid%km,                                    &
                                     modc=grid%modc_tracer )
                  call mp_recvtrirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                     &
                                     grid%ijk_yz_to_xy%RecvDesc, q_internal, dyn_out%tracer,      &
                                     mq, mlast, ntotq, 1, grid%im, grid%jfirst, grid%jlast, grid%kfirst, &
                                     grid%klast, grid%ifirstxy, grid%ilastxy, grid%jfirstxy,      &
                                     grid%jlastxy, 1, grid%km,                                    &
                                     modc=grid%modc_tracer )
               end do
            end if

            call t_stopf  ('yz_to_xy_q')

            if (high_alt) then
               ! Transpose pt (because pt correction is done after cd_core)

               call t_startf ('yz_to_xy_pt')
               call mp_sendirr( grid%commxy, grid%pt_to_ptxy%SendDesc,                    &
                    grid%pt_to_ptxy%RecvDesc, pt, ptxy,                 &
                    modc=grid%modc_dynrun )
               call mp_recvirr( grid%commxy, grid%pt_to_ptxy%SendDesc,                    &
                    grid%pt_to_ptxy%RecvDesc, pt, ptxy,                 &
                    modc=grid%modc_dynrun )
               call t_stopf ('yz_to_xy_pt')
            endif

#endif

         else

            call t_startf ('yz_to_xy_psuv')

            do j = jfirst, jlast
               do i = 1, im
                  psxy(i,j) = ps(i,j)
               end do
            end do

!$omp parallel do private(i,j,k)
            do k = kfirst, klast
               do j = jfirst, jlast
                  do i = 1, im
                     uxy(i,j,k) = u(i,j,k)
                     vxy(i,j,k) = v(i,j,k)
                  end do
               end do
            end do

            if (am_fixer.or.am_diag) then
!$omp parallel do private(i,j,k)
               do k = kfirst, klast
                  do j = jfirst, jlast
                     do i = 1, im
                        uc_i(i,j,k)= uc_s(i,j,k)
                        vc_i(i,j,k)= vc_s(i,j,k)
                     end do
                  end do
               end do
            end if

            if (high_alt) then
!$omp parallel do private(i,j,k)
               do k = kfirst, klast
                  do j = jfirst, jlast
                     do i = 1, im
                        ptxy(i,j,k) = pt(i,j,k)
                     end do
                  end do
               end do
            end if

            call t_stopf  ('yz_to_xy_psuv')

            call t_startf ('yz_to_xy_q')

!$omp parallel do private(i,j,k,mq)
            do mq = 1, ntotq

               ! Temporary -- here the pointers will ultimately be set, not the contents copied
               do k = 1, km
                  do j = jfirst, jlast
                     do i = 1, im
                        dyn_out%tracer(i,j,k,mq) = q_internal(i,j,k,mq)
                     end do
                  end do
               end do
            end do

            call t_stopf  ('yz_to_xy_q')

         end if  !  (grid%twod_decomp .eq. 1)

      end if  !  (grid%iam .lt. grid%npes_xy)

      if ( km > 1 ) then           ! not shallow water equations

         ! Perform vertical remapping from Lagrangian control-volume to
         ! the Eulerian coordinate as specified by the routine set_eta.
         ! Note that this finite-volume dycore is otherwise independent of the vertical
         ! Eulerian coordinate.


         ! te_map requires uxy, vxy, psxy, pexy, pkxy, phisxy, q3xy, and ptxy

         call t_barrierf('sync_te_map', grid%commdyn)
         call t_startf ('te_map')

         if (grid%iam .lt. grid%npes_xy) then

            call te_map(grid,     consv,   convt_local, psxy, omgaxy,       &
                        pexy,     delpxy,  pkzxy,  pkxy,   ndt,             &
                        nx,       uxy,     vxy,    ptxy,   dyn_out%tracer,  &
                        phisxy,   cp3v,    cap3v,  kord,   pelnxy,          &
                        te0,      tempxy,  dp0xy,  mfxxy,  mfyxy,           &
                        uc_i,     vc_i,  du_fix_s, du_fix_i,                &
                        am_geom_crrct, (am_fixer.or.am_diag) )

            if (am_diag) then
!$omp parallel do private(i,j,k)
                do j=jfirstxy,jlastxy
                  do k=1,km
                     do i=ifirstxy,ilastxy
                        ucxy(i,j,k)=ucxy(i,j,k)+uc_i(i,j,k)
                        vcxy(i,j,k)=vcxy(i,j,k)+vc_i(i,j,k)
                     enddo
                  enddo
               enddo
            end if

            if (am_fixer .or. am_diag) then
!$omp parallel do private(i,j,k)
               do j=jfirstxy,jlastxy
                  do k=1,km
                     do i=ifirstxy,ilastxy
                        dufix_xy(i,j,k)=dufix_xy(i,j,k)+du_fix_i(i,j,k)*grid%cose(j)
                     enddo
                  enddo
               enddo
            endif

            if ( .not. convt_local ) then
!$omp parallel do private(i,j,k)
               do j=jfirstxy,jlastxy
                  do k=1,km
                     do i=ifirstxy,ilastxy
                        t3xy(i,j,k) = ptxy(i,j,k)*pkzxy(i,j,k)/ &
                             (D1_0+zvir*dyn_out%tracer(i,j,k,1))
                     end do
                  end do
               end do
            end if

         end if

         call t_stopf ('te_map')

      end if  !  ( km > 1 )

      ! te_map computes uxy, vxy, psxy, delpxy, pexy, pkxy, pkzxy,
      ! pelnxy, omgaxy, tracer, ptxy, mfxxy and mfyxy

   end do  !    do iv = 1, nv - vertical re-mapping sub-cycle loop

   ! get total wind increments from dynamics timestep
   if (am_diag) then
!$omp parallel do private(i,j,k)
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               duxy(i,j,k) = uxy(i,j,k) - duxy(i,j,k)
               dvxy(i,j,k) = vxy(i,j,k) - dvxy(i,j,k)
            enddo
         enddo
      enddo
   end if

   call t_startf ('dyn_run_dealloc')

   deallocate( worka )
   deallocate( workb )
   deallocate( dp0 )
   deallocate( mfx )
   deallocate( mfy )
   deallocate(  cx )
   deallocate(  cy )
   deallocate( delpf )
   deallocate( uc    )
   deallocate( vc    )
   deallocate( dpt   )
   deallocate( dwz   )
   deallocate( pkc   )
   deallocate(  wz   )
   deallocate( pkcc )
   deallocate( wzc )
   deallocate( pkkp )
   deallocate( wzkp )
   deallocate( wzxy )
   deallocate( tempxy )
   deallocate( dp0xy )
   deallocate( psxy3 )
   deallocate( phisxy3 )
   deallocate( q_internal )
   deallocate (ctreq)
   deallocate (ctreqs)
   deallocate (cdcreqs)
#if defined(SPMD)
   deallocate (ctstat)
   deallocate (ctstats)
#endif
#if ( defined OFFLINE_DYN )
   deallocate( ps_obs )
   deallocate( ps_mod )
   deallocate( u_tmp )
   deallocate( v_tmp )
#endif

   if (am_fix_taper) then
      deallocate(zpkck)
   end if
   if (am_fixer.or.am_diag) then
      deallocate(du_fix_i)
      deallocate(du_k)
      deallocate(du_north)
      deallocate(uc_s)
      deallocate(vc_s)
      deallocate(uc_i)
      deallocate(vc_i)
   end if

   call t_stopf  ('dyn_run_dealloc')

   rc = DYN_RUN_SUCCESS

end subroutine dyn_run

!=============================================================================================

subroutine dyn_final(restart_file, dyn_state, dyn_in, dyn_out)

   use dynamics_vars, only: dynamics_clean

   character(LEN=*)             , intent(IN   ) :: restart_file
   type (T_FVDYCORE_STATE), target              :: dyn_state
   type (dyn_import_t), intent(inout)           :: dyn_in
   type (dyn_export_t), intent(inout)           :: dyn_out
   !-----------------------------------------------------------------------

   call dynamics_clean    ( dyn_state%grid  )
   call dyn_free_interface( dyn_in, dyn_out )

   !=============================================================================================
   contains
   !=============================================================================================

   subroutine dyn_free_interface ( dyn_in, dyn_out )

      ! free the dynamics import and export

      ! arguments
      type (dyn_import_t), intent(inout) :: dyn_in
      type (dyn_export_t), intent(inout) :: dyn_out
      !-----------------------------------------------------------------------

      if ( associated(dyn_in%phis) ) deallocate( dyn_in%phis )
      if ( associated(dyn_in%ps) )   deallocate( dyn_in%ps )
      if ( associated(dyn_in%u3s) )  deallocate( dyn_in%u3s )
      if ( associated(dyn_in%v3s) )  deallocate( dyn_in%v3s )
      if ( associated(dyn_in%pe) )   deallocate( dyn_in%pe )
      if ( associated(dyn_in%pt) )   deallocate( dyn_in%pt )
      if ( associated(dyn_in%t3) )   deallocate( dyn_in%t3 )
      if ( associated(dyn_in%pk) )   deallocate( dyn_in%pk )
      if ( associated(dyn_in%pkz) )  deallocate( dyn_in%pkz )
      if ( associated(dyn_in%delp) ) deallocate( dyn_in%delp )
      if ( associated(dyn_in%tracer) ) deallocate( dyn_in%tracer)

      if ( associated(dyn_out%ps) )   nullify( dyn_out%ps )
      if ( associated(dyn_out%u3s) )  nullify( dyn_out%u3s )
      if ( associated(dyn_out%v3s) )  nullify( dyn_out%v3s )
      if ( associated(dyn_out%pe) )   nullify( dyn_out%pe )
      if ( associated(dyn_out%pt) )   nullify( dyn_out%pt )
      if ( associated(dyn_out%t3) )   nullify( dyn_out%t3 )
      if ( associated(dyn_out%pk) )   nullify( dyn_out%pk )
      if ( associated(dyn_out%pkz) )  nullify( dyn_out%pkz )
      if ( associated(dyn_out%delp) ) nullify( dyn_out%delp )
      if ( associated(dyn_out%tracer) ) nullify( dyn_out%tracer )

      if ( associated(dyn_out%omga) ) deallocate( dyn_out%omga )
      if ( associated(dyn_out%peln) ) deallocate( dyn_out%peln )
      if ( associated(dyn_out%mfx) )  deallocate( dyn_out%mfx )
      if ( associated(dyn_out%mfy) )  deallocate( dyn_out%mfy )

   end subroutine dyn_free_interface

end subroutine dyn_final

!=============================================================================================
! Private routines
!=============================================================================================

subroutine read_inidat(dyn_in)

  ! Read initial dataset

  ! Arguments
  type (dyn_import_t), target, intent(inout) :: dyn_in

  ! Local variables
  integer                         :: ifirstxy, ilastxy, jfirstxy, jlastxy, km
  integer                         :: m, ntotq

  character(len=16)               :: fieldname

  type(file_desc_t), pointer      :: fh_ini    ! PIO filehandle

  type (t_fvdycore_grid), pointer :: grid
  ! variables for analytic initial conditions
  integer,  allocatable           :: glob_ind(:)
  integer,  allocatable           :: m_cnst(:)
  real(r8), allocatable           :: clon_st(:)
  integer                         :: nglon, nglat
  integer                         :: i, j
  integer                         :: jf, gf, uf ! First indices for setting u3s
  integer                         :: ierr
  integer                         :: lonid
  integer                         :: latid
  integer                         :: mlon ! longitude dimension length from dataset
  integer                         :: mlat            ! latitude dimension length from dataset
  real(r8), parameter             :: deg2rad = pi/180._r8

  character(len=*), parameter     :: sub='read_inidat'
  !----------------------------------------------------------------------------

  fh_ini  => initial_file_get_id()

  grid     => get_dyn_state_grid()
  ifirstxy =  grid%ifirstxy
  ilastxy  =  grid%ilastxy
  jfirstxy =  grid%jfirstxy
  jlastxy  =  grid%jlastxy
  km       =  grid%km
  ntotq    =  grid%ntotq

  ! Set the array initial_mr assuming that constituents are initialized with mixing ratios
  ! that are consistent with their declared type in the constituents module.  This array
  ! may be modified below to provide backwards compatibility for reading old initial files
  ! that contain wet mixing ratios for all constituents regardless of how they were registered.
  do i = 1, pcnst
     initial_mr(i) = cnst_type(i)
  end do

  if (analytic_ic_active()) then
     readvar   = .false.
    if (jfirstxy == 1) then
      jf = 1
      uf = 2
      gf = (ilastxy - ifirstxy + 1) + 1 ! Skip the first block of longitudes
    else
      jf = jfirstxy-1
      uf = jfirstxy
      gf = 1
    end if
    allocate(glob_ind((ilastxy - ifirstxy + 1) * (jlastxy - jfirstxy + 1)))
    call get_horiz_grid_dim_d(nglon, nglat)
    m = 1
    do j = jfirstxy, jlastxy
      do i = ifirstxy, ilastxy
        ! Create a global column index
        glob_ind(m) = i + (j-1)*nglon
        m = m + 1
      end do
    end do
    allocate(m_cnst(ntotq))
    do i = 1, ntotq
      m_cnst(i) = i
    end do
    allocate(clon_st(ifirstxy:ilastxy))
    clon_st(ifirstxy:ilastxy) = londeg_st(ifirstxy:ilastxy,1) * deg2rad

    call analytic_ic_set_ic(vc_moist_pressure, clat_staggered(jf:jlastxy-1),  &
         clon(ifirstxy:ilastxy,1), glob_ind(gf:), U=dyn_in%u3s(:,uf:,:))
    call analytic_ic_set_ic(vc_moist_pressure, clat(jfirstxy:jlastxy),        &
         clon_st(ifirstxy:ilastxy), glob_ind, V=dyn_in%v3s)
    ! Note that analytic_ic_set_ic makes use of cnst_init_default for
    ! the tracers except water vapor.
    call analytic_ic_set_ic(vc_moist_pressure, clat(jfirstxy:jlastxy),        &
         clon(ifirstxy:ilastxy,1), glob_ind, T=dyn_in%t3, PS=dyn_in%ps,       &
         Q=dyn_in%tracer(:,:,:,1:ntotq), PHIS_IN=dyn_in%phis, m_cnst=m_cnst)
    do m = 1, ntotq
       call process_inidat(grid, dyn_in, 'CONSTS', m_cnst=m, fh_ini=fh_ini)
    end do
    deallocate(glob_ind)
    deallocate(m_cnst)
    deallocate(clon_st)

  else

    !-----------
    ! Check coord sizes
    !-----------
     ierr = pio_inq_dimid(fh_ini, 'lon' , lonid)
     ierr = pio_inq_dimid(fh_ini, 'lat' , latid)
     ierr = pio_inq_dimlen(fh_ini, lonid , mlon)
     ierr = pio_inq_dimlen(fh_ini, latid , mlat)
     if (mlon /= plon .or. mlat /= plat) then
        write(iulog,*) sub//': ERROR: model parameters do not match initial dataset parameters'
        write(iulog,*)'Model Parameters:    plon = ',plon,' plat = ',plat
        write(iulog,*)'Dataset Parameters:  dlon = ',mlon,' dlat = ',mlat
        call endrun(sub//': ERROR: model parameters do not match initial dataset parameters')
     end if

    !-----------
    ! 2-D fields
    !-----------

    fieldname = 'PS'
    call infld(fieldname, fh_ini, 'lon', 'lat', ifirstxy, ilastxy, jfirstxy, jlastxy, &
         dyn_in%ps, readvar, gridname='fv_centers')
    if (.not. readvar) call endrun(sub//': ERROR: PS not found')


    !-----------
    ! 3-D fields
    !-----------

    fieldname = 'US'
    call infld(fieldname, fh_ini, 'lon', 'slat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
         1, km, dyn_in%u3s, readvar, gridname='fv_u_stagger')
    if (.not. readvar) call endrun(sub//': ERROR: US not found')

    fieldname = 'VS'
    call infld(fieldname, fh_ini, 'slon', 'lat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
         1, km, dyn_in%v3s, readvar, gridname='fv_v_stagger')
    if (.not. readvar) call endrun(sub//': ERROR: VS not found')

    fieldname = 'T'
    call infld(fieldname, fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
         1, km, dyn_in%t3, readvar, gridname='fv_centers')
    if (.not. readvar) call endrun(sub//': ERROR: T not found')
  end if

  ! Constituents (read and process one at a time)
  !
  ! If analytic ICs are being used, we allow constituents in an initial
  ! file to overwrite mixing ratios set by the default constituent initialization
  ! except for the water species.
  !
  ! If using analytic ICs the initial file only needs the horizonal grid
  ! dimension checked in the case that the file contains constituent mixing
  ! ratios.
  if (analytic_ic_active()) then
     do m = 1, pcnst
        if (cnst_read_iv(m) .and. .not. cnst_is_a_water_species(cnst_name(m))) then
           if (dyn_field_exists(fh_ini, trim(cnst_name(m)), required=.false.)) then
              ierr = pio_inq_dimid(fh_ini, 'lon' , lonid)
              ierr = pio_inq_dimid(fh_ini, 'lat' , latid)
              ierr = pio_inq_dimlen(fh_ini, lonid , mlon)
              ierr = pio_inq_dimlen(fh_ini, latid , mlat)
              if (mlon /= plon .or. mlat /= plat) then
                 write(iulog,*) sub//': ERROR: model parameters do not match initial dataset parameters'
                 write(iulog,*)'Model Parameters:    plon = ',plon,' plat = ',plat
                 write(iulog,*)'Dataset Parameters:  dlon = ',mlon,' dlat = ',mlat
                 call endrun(sub//': ERROR: model parameters do not match initial dataset parameters')
              end if
              exit
           end if
        end if
     end do
  end if

  do m = 1, pcnst

    if (analytic_ic_active() .and. cnst_is_a_water_species(cnst_name(m))) cycle

    readvar   = .false.
    fieldname = cnst_name(m)
    if (cnst_read_iv(m) .and. dyn_field_exists(fh_ini, trim(cnst_name(m)), required=.false.)) then
      call infld(fieldname, fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
           1, km, dyn_in%tracer(:,:,:,m), readvar, gridname='fv_centers')
    end if
    call process_inidat(grid, dyn_in, 'CONSTS', m_cnst=m, fh_ini=fh_ini)
  end do

  ! Set u3s(:,1,:) to zero as it is used in interpolation routines
  if ((jfirstxy == 1) .and. (size(dyn_in%u3s) > 0)) then
    dyn_in%u3s(ifirstxy:ilastxy,jfirstxy,1:km) = 0.0_r8
  end if

  ! These always happen
  call process_inidat(grid, dyn_in, 'PS')
  call process_inidat(grid, dyn_in, 'T')

end subroutine read_inidat

!=========================================================================================

subroutine set_phis(dyn_in)

   ! Set PHIS according to the following rules.
   !
   ! 1) If a topo file is specified use it.  This option has highest precedence.
   ! 2) If not using topo file, but analytic_ic option is on, use analytic phis.
   ! 3) Set phis = 0.0.

   ! Arguments
   type (dyn_import_t), target, intent(inout) :: dyn_in   ! dynamics import

   ! local variables
   type(file_desc_t),      pointer :: fh_topo
   type (t_fvdycore_grid), pointer :: grid

   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
   integer :: ierr
   integer :: lonid
   integer :: latid
   integer :: mlon            ! longitude dimension length from dataset
   integer :: mlat            ! latitude dimension length from dataset
   integer :: nglon, nglat
   integer :: i, j, m

   integer, allocatable :: glob_ind(:)

   character(len=16)                :: fieldname
   character(len=*), parameter      :: sub='set_phis'
   !----------------------------------------------------------------------------

   fh_topo => topo_file_get_id()

   grid     => get_dyn_state_grid()
   ifirstxy =  grid%ifirstxy
   ilastxy  =  grid%ilastxy
   jfirstxy =  grid%jfirstxy
   jlastxy  =  grid%jlastxy

   if (associated(fh_topo)) then    
      !-----------
      ! Check coord sizes
      !-----------
      ierr = pio_inq_dimid(fh_topo, 'lon' , lonid)
      ierr = pio_inq_dimid(fh_topo, 'lat' , latid)
      ierr = pio_inq_dimlen(fh_topo, lonid , mlon)
      ierr = pio_inq_dimlen(fh_topo, latid , mlat)
      if (mlon /= plon .or. mlat /= plat) then
         write(iulog,*) sub//': ERROR: model parameters do not match topo dataset parameters'
         write(iulog,*)'Model Parameters:    plon = ',plon,' plat = ',plat
         write(iulog,*)'Dataset Parameters:  dlon = ',mlon,' dlat = ',mlat
         call endrun(sub//': ERROR: model parameters do not match topo dataset parameters')
      end if
    
      fieldname = 'PHIS'
      readvar   = .false.      
      call infld(fieldname, fh_topo, 'lon', 'lat', ifirstxy, ilastxy, jfirstxy, jlastxy, &
         dyn_in%phis, readvar, gridname='fv_centers')
      if (.not. readvar) call endrun(sub//': ERROR: PHIS not found')

   else if (analytic_ic_active()) then

      allocate(glob_ind((ilastxy - ifirstxy + 1) * (jlastxy - jfirstxy + 1)))
      call get_horiz_grid_dim_d(nglon, nglat)
      m = 1
      do j = jfirstxy, jlastxy
         do i = ifirstxy, ilastxy
            ! Create a global column index
            glob_ind(m) = i + (j-1)*nglon
            m = m + 1
         end do
      end do

      call analytic_ic_set_ic(vc_moist_pressure, clat(jfirstxy:jlastxy),      &
         clon(ifirstxy:ilastxy,1), glob_ind, PHIS_OUT=dyn_in%phis)

   else

      dyn_in%phis(:,:) = 0._r8

   end if

   call process_inidat(grid, dyn_in, 'PHIS')

end subroutine set_phis

!=========================================================================================

subroutine process_inidat(grid, dyn_in, fieldname, m_cnst, fh_ini)

   ! Post-process input fields
   use commap,              only: clat, clon
   use const_init,          only: cnst_init_default
   use inic_analytic,       only: analytic_ic_active

   ! arguments
   type(t_fvdycore_grid), target, intent(inout) :: grid        ! dynamics state grid
   type(dyn_import_t),    target, intent(inout) :: dyn_in      ! dynamics import
   character(len=*),              intent(in)    :: fieldname   ! field to be processed
   integer,             optional, intent(in)    :: m_cnst      ! constituent index
   type(file_desc_t),   optional, intent(inout) :: fh_ini

   ! Local variables
   integer :: i, j, k                     ! grid and constituent indices
   integer :: npes_xy
   integer :: im, jm, km
   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy

   integer :: nglon, nglat, rndm_seed_sz
   integer, allocatable :: rndm_seed(:)
   real(r8) :: pertval                       ! perturbation value

   real(r8), pointer :: phisxy(:,:), psxy(:,:), t3xy(:,:,:)
   real(r8), pointer :: tracer(:,:,:,:)

   real(r8) :: xsum(grid%km)               ! temp array for parallel sums

   integer :: err_handling
   integer :: varid                        ! variable id
   integer :: ret                          ! return values
   character(len=256) :: trunits           ! tracer untis
   character(len=3)   :: mixing_ratio

   character(len=*), parameter :: sub='process_inidat'
   !----------------------------------------------------------------------------

   psxy   => dyn_in%ps
   phisxy => dyn_in%phis
   t3xy   => dyn_in%t3

   npes_xy  = grid%npes_xy
   im       = grid%im
   jm       = grid%jm
   km       = grid%km
   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   select case (fieldname)

   case ('T')

      if (iam >= npes_xy) return

      ! Add random perturbation to temperature if requested
      if ((pertlim /= 0._r8) .and. (.not. analytic_ic_active())) then

         if (masterproc) then
            write(iulog,*) sub//':  Adding random perturbation bounded by +/-', &
                           pertlim,' to initial temperature field'
         end if

         call get_horiz_grid_dim_d(nglon, nglat)
         call random_seed(size=rndm_seed_sz)
         allocate(rndm_seed(rndm_seed_sz))

         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               ! seed random_number generator based on global column index
               rndm_seed = i + (j-1)*nglon
               call random_seed(put=rndm_seed)
               do k = 1, km
                  call random_number(pertval)
                  pertval = 2._r8*pertlim*(0.5_r8 - pertval)
                  t3xy(i,j,k) = t3xy(i,j,k)*(1._r8 + pertval)
               end do
            end do
         end do

         deallocate(rndm_seed)
      end if

      ! Average T at the poles.
      if (jfirstxy == 1) then
         call par_xsum(grid, t3xy(:,1,:), km, xsum)
         do k=1, km
            do i = ifirstxy, ilastxy
               t3xy(i,1,k) = xsum(k) / real(im,r8)
            end do
         end do
      end if
      if (jlastxy == jm) then
         call par_xsum(grid, t3xy(:,jm,:), km, xsum)
         do k = 1, km
            do i = ifirstxy, ilastxy
               t3xy(i,jm,k) = xsum(k) / real(im,r8)
            end do
         end do
      end if

   case ('CONSTS')

      if (.not. present(m_cnst)) then
         call endrun(sub//': ERROR:  m_cnst needs to be present in the'// &
                     ' argument list')
      end if

      tracer => dyn_in%tracer

      if (readvar) then

         if (.not. present(fh_ini)) then
            call endrun(sub//': ERROR:  fh_ini needs to be present in the'// &
                        ' argument list')
         end if

         ! Check that all tracer units are in mass mixing ratios
         ret = pio_inq_varid(fh_ini, cnst_name(m_cnst), varid)
         ret = pio_get_att (fh_ini, varid, 'units', trunits)
         if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
            call endrun(sub//': ERROR:  Units for tracer ' &
                  //trim(cnst_name(m_cnst))//' must be in KG/KG')
         end if

         ! Check for mixing_ratio attribute.  If present then use it to
         ! specify whether the initial file contains wet or dry values.  If
         ! not present then assume the mixing ratio is wet.  This is for
         ! backwards compatibility with old initial files that were written
         ! will all wet mixing ratios.

         ! We will handle errors for this routine
         call pio_seterrorhandling(fh_ini, pio_bcast_error, err_handling)

         ret = pio_get_att(fh_ini, varid, 'mixing_ratio', mixing_ratio)
         if (ret == pio_noerr) then
            initial_mr(m_cnst) = mixing_ratio
         else
            initial_mr(m_cnst) = 'wet'
         end if

         ! reset PIO to handle errors as before
         call pio_seterrorhandling(fh_ini, err_handling)


      else if (.not. analytic_ic_active()) then

         ! Constituents not read from initial file are initialized by the
         ! package that implements them.  Note that the analytic IC code calls
         ! cnst_init_default internally.

         if (iam >= npes_xy) return

         if (m_cnst == 1 .and. moist_physics) then
            call endrun(sub//': ERROR:  Q must be on Initial File')
         end if

         call cnst_init_default(m_cnst, clat(jfirstxy:jlastxy), clon(ifirstxy:ilastxy,1), tracer(:,:,:,m_cnst))
      end if

      if (.not. analytic_ic_active()) then
         do k = 1, km
            do j = jfirstxy, jlastxy
               do i = ifirstxy, ilastxy
                  tracer(i,j,k,m_cnst) = max(tracer(i,j,k,m_cnst), qmin(m_cnst))
               end do
            end do
         end do
      end if

      if (iam >= npes_xy) return

      ! Compute polar average
      if (jfirstxy == 1) then
         call par_xsum(grid, dyn_in%tracer(:,1,:,m_cnst), km, xsum)
         do k = 1, km
            do i = ifirstxy, ilastxy
               dyn_in%tracer(i,1,k,m_cnst) = xsum(k) / real(im,r8)
            end do
         end do
      end if
      if (jlastxy == jm) then
         call par_xsum(grid, dyn_in%tracer(:,jm,:,m_cnst), km, xsum)
         do k = 1, km
            do i = ifirstxy, ilastxy
               dyn_in%tracer(i,jm,k,m_cnst) = xsum(k) / real(im,r8)
            end do
         end do
      end if

   case ('PS')

      ! Average PS at the poles.
      if (jfirstxy == 1) then
         if (size(psxy,2) > 0) then
            call par_xsum(grid, psxy(:,1:1), 1, xsum(1:1))
            do i = ifirstxy, ilastxy
               psxy(i,1) = xsum(1) / real(im,r8)
            end do
         end if
      end if
      if (jlastxy == jm) then
         call par_xsum(grid, psxy(:,jm:jm), 1, xsum(1:1))
         do i = ifirstxy, ilastxy
            psxy(i,jm) = xsum(1) / real(im,r8)
         end do
      end if

   case ('PHIS')

      ! Average PHIS at the poles.
      if (jfirstxy == 1) then
         if (size(phisxy,2) > 0) then
            call par_xsum(grid, phisxy(:,1:1), 1, xsum(1:1))
            do i = ifirstxy, ilastxy
               phisxy(i,1) = xsum(1) / real(im,r8)
            end do
         end if
      end if
      if (jlastxy == jm) then
         call par_xsum(grid, phisxy(:,jm:jm), 1, xsum(1:1))
         do i = ifirstxy, ilastxy
            phisxy(i,jm) = xsum(1) / real(im,r8)
         end do
      end if

   end select

end subroutine process_inidat

!=========================================================================================

logical function dyn_field_exists(fh, fieldname, required)

   use pio,            only: var_desc_t, PIO_inq_varid
   use pio,            only: PIO_NOERR

   type(file_desc_t), intent(inout) :: fh ! needs to be inout because of pio_seterrorhandling
   character(len=*),  intent(in) :: fieldname
   logical, optional, intent(in) :: required

   ! Local variables
   logical                  :: found
   logical                  :: field_required
   integer                  :: ret
   integer                  :: pio_errtype
   type(var_desc_t)         :: varid
   character(len=128)       :: errormsg
   !--------------------------------------------------------------------------

   if (present(required)) then
      field_required = required
   else
      field_required = .true.
   end if

   ! Set PIO to return error codes when reading data from IC file.
   call pio_seterrorhandling(fh, pio_bcast_error, oldmethod=pio_errtype)

   ret = PIO_inq_varid(fh, trim(fieldname), varid)
   found = (ret == PIO_NOERR)
   if (.not. found) then
      if (field_required) then
         write(errormsg, *) trim(fieldname),' was not present in the input file.'
         call endrun('DYN_FIELD_EXISTS: '//errormsg)
      end if
   end if

   dyn_field_exists = found

   ! Put the error handling back the way it was
   call pio_seterrorhandling(fh, pio_errtype)

end function dyn_field_exists

!=========================================================================================

end module dyn_comp
