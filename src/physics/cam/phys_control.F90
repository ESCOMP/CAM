module phys_control
!-----------------------------------------------------------------------
! Purpose:
!
! Provides a control interface to CAM physics packages
!
! Revision history:
! 2006-05-01  D. B. Coleman,  Creation of module
! 2009-02-13  Eaton           Replace *_{default,set}opts methods with module namelist.
!                             Add vars to indicate physics version and chemistry type.
!-----------------------------------------------------------------------

use spmd_utils,     only: masterproc
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use shr_kind_mod,   only: r8 => shr_kind_r8

implicit none
private
save

public :: &
   phys_ctl_readnl,   &! read namelist from file
   phys_getopts,      &! generic query method
   phys_setopts,      &! generic set method
   phys_deepconv_pbl, &! return true if deep convection is allowed in the PBL
   phys_do_flux_avg,  &! return true to average surface fluxes
   cam_physpkg_is,    &! query for the name of the physics package
   cam_chempkg_is,    &! query for the name of the chemistry package
   waccmx_is

! Private module data

character(len=16), parameter :: unset_str = 'UNSET'
integer,           parameter :: unset_int = huge(1)

! Namelist variables:
character(len=16) :: cam_physpkg          = unset_str  ! CAM physics package
character(len=32) :: cam_chempkg          = unset_str  ! CAM chemistry package 
character(len=16) :: waccmx_opt           = unset_str  ! WACCMX run option [ionosphere | neutral | off
character(len=16) :: deep_scheme          = unset_str  ! deep convection package
character(len=16) :: shallow_scheme       = unset_str  ! shallow convection package
character(len=16) :: eddy_scheme          = unset_str  ! vertical diffusion package
character(len=16) :: microp_scheme        = unset_str  ! microphysics package
character(len=16) :: macrop_scheme        = unset_str  ! macrophysics package
character(len=16) :: radiation_scheme     = unset_str  ! radiation package
integer           :: srf_flux_avg         = unset_int  ! 1 => smooth surface fluxes, 0 otherwise

logical           :: use_subcol_microp    = .false.    ! if .true. then use sub-columns in microphysics

logical           :: atm_dep_flux         = .true.     ! true => deposition fluxes will be provided
                                                       ! to the coupler
logical           :: history_amwg         = .true.     ! output the variables used by the AMWG diag package
logical           :: history_vdiag        = .false.    ! output the variables used by the AMWG variability diag package
logical           :: history_aerosol      = .false.    ! output the MAM aerosol variables and tendencies
logical           :: history_aero_optics  = .false.    ! output the aerosol
logical           :: history_eddy         = .false.    ! output the eddy variables
logical           :: history_budget       = .false.    ! output tendencies and state variables for CAM4
                                                       ! temperature, water vapor, cloud ice and cloud
                                                       ! liquid budgets.
logical           :: convproc_do_aer      = .false.    ! switch for new convective scavenging treatment for modal aerosols

integer           :: history_budget_histfile_num = 1   ! output history file number for budget fields
logical           :: history_waccm        = .false.    ! output variables of interest for WACCM runs
logical           :: history_waccmx       = .false.    ! output variables of interest for WACCM-X runs
logical           :: history_chemistry    = .true.     ! output default chemistry-related variables
logical           :: history_carma        = .false.    ! output default CARMA-related variables
logical           :: history_clubb        = .true.     ! output default CLUBB-related variables
logical           :: history_cesm_forcing = .false.
logical           :: history_dust         = .false.
logical           :: history_scwaccm_forcing = .false.
logical           :: history_chemspecies_srf = .false.

logical           :: do_clubb_sgs
! Check validity of physics_state objects in physics_update.
logical           :: state_debug_checks   = .false.

! Macro/micro-physics co-substeps
integer           :: cld_macmic_num_steps = 1

logical           :: offline_driver       = .false.    ! true => offline driver is being used


logical, public, protected :: use_simple_phys = .false. ! true => simple physics configuration

logical :: use_spcam       ! true => use super parameterized CAM

logical :: prog_modal_aero ! determines whether prognostic modal aerosols are present in the run.

! Option to use heterogeneous freezing
logical, public, protected :: use_hetfrz_classnuc = .false.

! Which gravity wave sources are used?
logical, public, protected :: use_gw_oro = .true.         ! Orography.
logical, public, protected :: use_gw_front = .false.      ! Frontogenesis.
logical, public, protected :: use_gw_front_igw = .false.  ! Frontogenesis to inertial spectrum.
logical, public, protected :: use_gw_convect_dp = .false. ! Deep convection.
logical, public, protected :: use_gw_convect_sh = .false. ! Shallow convection.

! FV dycore angular momentum correction
logical, public, protected :: fv_am_correction = .false.

! CAM snapshot before/after file numbers and control
character(len=32) :: cam_take_snapshot_before = ''  ! Physics routine to take a snopshot "before"
character(len=32) :: cam_take_snapshot_after = ''   ! Physics routine to take a snopshot "after"
integer           :: cam_snapshot_before_num = -1   ! output history file number for CAM "before" snapshot
integer           :: cam_snapshot_after_num = -1    ! output history file number for CAM "after" snapshot

!======================================================================= 
contains
!======================================================================= 

subroutine phys_ctl_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpi_character, mpi_integer, mpi_logical, masterprocid, mpicom
   use cam_control_mod, only: cam_ctrl_set_physics_type

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_ctl_readnl'

   namelist /phys_ctl_nl/ cam_physpkg, use_simple_phys, cam_chempkg, waccmx_opt,  &
      deep_scheme, shallow_scheme, &
      eddy_scheme, microp_scheme,  macrop_scheme, radiation_scheme, srf_flux_avg, &
      use_subcol_microp, atm_dep_flux, history_amwg, history_vdiag, history_aerosol, history_aero_optics, &
      history_eddy, history_budget,  history_budget_histfile_num, history_waccm, &
      history_waccmx, history_chemistry, history_carma, history_clubb, history_dust, &
      history_cesm_forcing, history_scwaccm_forcing, history_chemspecies_srf, &
      do_clubb_sgs, state_debug_checks, use_hetfrz_classnuc, use_gw_oro, use_gw_front, &
      use_gw_front_igw, use_gw_convect_dp, use_gw_convect_sh, cld_macmic_num_steps, &
      offline_driver, convproc_do_aer, cam_snapshot_before_num, cam_snapshot_after_num, &
      cam_take_snapshot_before, cam_take_snapshot_after
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_ctl_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_ctl_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(deep_scheme,                 len(deep_scheme),      mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(cam_physpkg,                 len(cam_physpkg),      mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(use_simple_phys,             1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(cam_chempkg,                 len(cam_chempkg),      mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(waccmx_opt,                  len(waccmx_opt),       mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(shallow_scheme,              len(shallow_scheme),   mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(eddy_scheme,                 len(eddy_scheme),      mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(microp_scheme,               len(microp_scheme),    mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(radiation_scheme,            len(radiation_scheme), mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(macrop_scheme,               len(macrop_scheme),    mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(srf_flux_avg,                1,                     mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_subcol_microp,           1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(atm_dep_flux,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_amwg,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_vdiag,               1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_eddy,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_aerosol,             1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_aero_optics,         1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_budget,              1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_budget_histfile_num, 1,                     mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_waccm,               1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_waccmx,              1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_chemistry,           1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_carma,               1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_clubb,               1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_cesm_forcing,        1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_chemspecies_srf,     1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_dust,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(history_scwaccm_forcing,     1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(do_clubb_sgs,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(state_debug_checks,          1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_hetfrz_classnuc,         1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_gw_oro,                  1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_gw_front,                1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_gw_front_igw,            1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_gw_convect_dp,           1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(use_gw_convect_sh,           1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(cld_macmic_num_steps,        1,                     mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(offline_driver,              1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(convproc_do_aer,             1,                     mpi_logical,   masterprocid, mpicom, ierr)
   call mpi_bcast(cam_snapshot_before_num,     1,                     mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(cam_snapshot_after_num,      1,                     mpi_integer,   masterprocid, mpicom, ierr)
   call mpi_bcast(cam_take_snapshot_before,    len(cam_take_snapshot_before), mpi_character, masterprocid, mpicom, ierr)
   call mpi_bcast(cam_take_snapshot_after,     len(cam_take_snapshot_after),  mpi_character, masterprocid, mpicom, ierr)

   use_spcam       = (     cam_physpkg_is('spcam_sam1mom') &
                      .or. cam_physpkg_is('spcam_m2005'))

   call cam_ctrl_set_physics_type(cam_physpkg)

   ! Error checking:

   ! Check compatibility of eddy & shallow schemes
   if (( shallow_scheme .eq. 'UW' ) .and. ( eddy_scheme .ne. 'diag_TKE' )) then
      write(iulog,*)'Do you really want to run UW shallow scheme without diagnostic TKE eddy scheme? Quiting'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   if (( shallow_scheme .eq. 'Hack' ) .and. ( ( eddy_scheme .ne. 'HB' ) .and. ( eddy_scheme .ne. 'HBR' ))) then
      write(iulog,*)'Do you really want to run Hack shallow scheme with a non-standard eddy scheme? Quiting.'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   ! Check compatibility of PBL and Microphysics schemes
   if (( eddy_scheme .eq. 'diag_TKE' ) .and. ( microp_scheme .eq. 'RK' )) then
      write(iulog,*)'UW PBL is not compatible with RK microphysics.  Quiting'
      call endrun('PBL and Microphysics schemes incompatible')
   endif
   
   ! Add a check to make sure CLUBB and MG are used together
   if ( do_clubb_sgs .and. ( microp_scheme .ne. 'MG') .and. .not. use_spcam) then
      write(iulog,*)'CLUBB is only compatible with MG microphysics.  Quiting'
      call endrun('CLUBB and microphysics schemes incompatible')
   endif

   ! Check that eddy_scheme, macrop_scheme, shallow_scheme are all set to CLUBB_SGS if do_clubb_sgs is true
   if (do_clubb_sgs .and. .not. use_spcam) then
      if (eddy_scheme .ne. 'CLUBB_SGS' .or. macrop_scheme .ne. 'CLUBB_SGS' .or. shallow_scheme .ne. 'CLUBB_SGS') then
         write(iulog,*)'eddy_scheme, macrop_scheme and shallow_scheme must all be CLUBB_SGS.  Quiting'
         call endrun('CLUBB and eddy, macrop or shallow schemes incompatible')
      endif
   endif
      
   ! Macro/micro co-substepping support.
   if (cld_macmic_num_steps > 1) then
      if (microp_scheme /= "MG" .or. (macrop_scheme /= "park" .and. macrop_scheme /= "CLUBB_SGS")) then
         call endrun ("Setting cld_macmic_num_steps > 1 is only &
              &supported with Park or CLUBB macrophysics and MG microphysics.")
      end if
   end if

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   prog_modal_aero = index(cam_chempkg,'_mam')>0

end subroutine phys_ctl_readnl

!===============================================================================

logical function cam_physpkg_is(name)

   ! query for the name of the physics package

   character(len=*) :: name
   
   cam_physpkg_is = (trim(name) == trim(cam_physpkg))
end function cam_physpkg_is

!===============================================================================

logical function cam_chempkg_is(name)

   ! query for the name of the chemics package

   character(len=*) :: name
   
   cam_chempkg_is = (trim(name) == trim(cam_chempkg))
end function cam_chempkg_is

!===============================================================================

logical function waccmx_is(name)

   ! query for the name of the waccmx run option

   character(len=*) :: name
   
   waccmx_is = (trim(name) == trim(waccmx_opt))
end function waccmx_is

!===============================================================================

subroutine phys_getopts(deep_scheme_out, shallow_scheme_out, eddy_scheme_out, microp_scheme_out, &
                        radiation_scheme_out, use_subcol_microp_out, atm_dep_flux_out, &
                         history_amwg_out, history_vdiag_out, history_aerosol_out, history_aero_optics_out, history_eddy_out, &
                        history_budget_out, history_budget_histfile_num_out, &
                        history_waccm_out, history_waccmx_out, history_chemistry_out, &
                        history_carma_out, history_clubb_out, history_dust_out, &
                        history_cesm_forcing_out, history_scwaccm_forcing_out, history_chemspecies_srf_out, &
                        cam_chempkg_out, prog_modal_aero_out, macrop_scheme_out, &
                        do_clubb_sgs_out, use_spcam_out, state_debug_checks_out, cld_macmic_num_steps_out, &
                        offline_driver_out, convproc_do_aer_out, cam_snapshot_before_num_out, cam_snapshot_after_num_out,&
                        cam_take_snapshot_before_out, cam_take_snapshot_after_out)
!-----------------------------------------------------------------------
! Purpose: Return runtime settings
!          deep_scheme_out   : deep convection scheme
!          shallow_scheme_out: shallow convection scheme
!          eddy_scheme_out   : vertical diffusion scheme
!          microp_scheme_out : microphysics scheme
!          radiation_scheme_out : radiation_scheme
!	   SPCAM_microp_scheme_out : SPCAM microphysics scheme
!-----------------------------------------------------------------------

   character(len=16), intent(out), optional :: deep_scheme_out
   character(len=16), intent(out), optional :: shallow_scheme_out
   character(len=16), intent(out), optional :: eddy_scheme_out
   character(len=16), intent(out), optional :: microp_scheme_out
   character(len=16), intent(out), optional :: radiation_scheme_out
   character(len=16), intent(out), optional :: macrop_scheme_out
   logical,           intent(out), optional :: use_subcol_microp_out
   logical,           intent(out), optional :: use_spcam_out
   logical,           intent(out), optional :: atm_dep_flux_out
   logical,           intent(out), optional :: history_amwg_out
   logical,           intent(out), optional :: history_vdiag_out
   logical,           intent(out), optional :: history_eddy_out
   logical,           intent(out), optional :: history_aerosol_out
   logical,           intent(out), optional :: history_aero_optics_out
   logical,           intent(out), optional :: history_budget_out
   integer,           intent(out), optional :: history_budget_histfile_num_out
   logical,           intent(out), optional :: history_waccm_out
   logical,           intent(out), optional :: history_waccmx_out
   logical,           intent(out), optional :: history_chemistry_out
   logical,           intent(out), optional :: history_carma_out
   logical,           intent(out), optional :: history_clubb_out
   logical,           intent(out), optional :: history_cesm_forcing_out
   logical,           intent(out), optional :: history_chemspecies_srf_out
   logical,           intent(out), optional :: history_dust_out
   logical,           intent(out), optional :: history_scwaccm_forcing_out
   logical,           intent(out), optional :: do_clubb_sgs_out
   character(len=32), intent(out), optional :: cam_chempkg_out
   logical,           intent(out), optional :: prog_modal_aero_out
   logical,           intent(out), optional :: state_debug_checks_out
   integer,           intent(out), optional :: cld_macmic_num_steps_out
   logical,           intent(out), optional :: offline_driver_out
   logical,           intent(out), optional :: convproc_do_aer_out
   integer,           intent(out), optional :: cam_snapshot_before_num_out
   integer,           intent(out), optional :: cam_snapshot_after_num_out
   character(len=32), intent(out), optional :: cam_take_snapshot_before_out
   character(len=32), intent(out), optional :: cam_take_snapshot_after_out

   if ( present(deep_scheme_out         ) ) deep_scheme_out          = deep_scheme
   if ( present(shallow_scheme_out      ) ) shallow_scheme_out       = shallow_scheme
   if ( present(eddy_scheme_out         ) ) eddy_scheme_out          = eddy_scheme
   if ( present(microp_scheme_out       ) ) microp_scheme_out        = microp_scheme
   if ( present(radiation_scheme_out    ) ) radiation_scheme_out     = radiation_scheme
   if ( present(use_subcol_microp_out   ) ) use_subcol_microp_out    = use_subcol_microp
   if ( present(use_spcam_out           ) ) use_spcam_out            = use_spcam

   if ( present(macrop_scheme_out       ) ) macrop_scheme_out        = macrop_scheme
   if ( present(atm_dep_flux_out        ) ) atm_dep_flux_out         = atm_dep_flux
   if ( present(history_aerosol_out     ) ) history_aerosol_out      = history_aerosol
   if ( present(history_aero_optics_out ) ) history_aero_optics_out  = history_aero_optics
   if ( present(history_budget_out      ) ) history_budget_out       = history_budget
   if ( present(history_amwg_out        ) ) history_amwg_out         = history_amwg
   if ( present(history_vdiag_out       ) ) history_vdiag_out        = history_vdiag
   if ( present(history_eddy_out        ) ) history_eddy_out         = history_eddy
   if ( present(history_budget_histfile_num_out ) ) history_budget_histfile_num_out = history_budget_histfile_num
   if ( present(history_waccm_out       ) ) history_waccm_out        = history_waccm
   if ( present(history_waccmx_out      ) ) history_waccmx_out       = history_waccmx
   if ( present(history_chemistry_out   ) ) history_chemistry_out    = history_chemistry
   if ( present(history_cesm_forcing_out) ) history_cesm_forcing_out = history_cesm_forcing
   if ( present(history_chemspecies_srf_out) ) history_chemspecies_srf_out = history_chemspecies_srf
   if ( present(history_scwaccm_forcing_out) ) history_scwaccm_forcing_out = history_scwaccm_forcing
   if ( present(history_carma_out       ) ) history_carma_out        = history_carma
   if ( present(history_clubb_out       ) ) history_clubb_out        = history_clubb
   if ( present(history_dust_out        ) ) history_dust_out         = history_dust
   if ( present(do_clubb_sgs_out        ) ) do_clubb_sgs_out         = do_clubb_sgs
   if ( present(cam_chempkg_out         ) ) cam_chempkg_out          = cam_chempkg
   if ( present(prog_modal_aero_out     ) ) prog_modal_aero_out      = prog_modal_aero
   if ( present(state_debug_checks_out  ) ) state_debug_checks_out   = state_debug_checks
   if ( present(cld_macmic_num_steps_out) ) cld_macmic_num_steps_out = cld_macmic_num_steps
   if ( present(offline_driver_out      ) ) offline_driver_out       = offline_driver
   if ( present(convproc_do_aer_out     ) ) convproc_do_aer_out      = convproc_do_aer
   if ( present(cam_snapshot_before_num_out ) ) cam_snapshot_before_num_out = cam_snapshot_before_num
   if ( present(cam_snapshot_after_num_out  ) ) cam_snapshot_after_num_out  = cam_snapshot_after_num
   if ( present(cam_take_snapshot_before_out) ) cam_take_snapshot_before_out = cam_take_snapshot_before
   if ( present(cam_take_snapshot_after_out ) ) cam_take_snapshot_after_out  = cam_take_snapshot_after

end subroutine phys_getopts

!===============================================================================

subroutine phys_setopts(fv_am_correction_in)

   logical, intent(in), optional :: fv_am_correction_in

   if ( present(fv_am_correction_in) ) fv_am_correction = fv_am_correction_in

end subroutine phys_setopts

!===============================================================================

function phys_deepconv_pbl()

  logical phys_deepconv_pbl

   ! Don't allow deep convection in PBL if running UW PBL scheme
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq. 'UW' ) ) then
      phys_deepconv_pbl = .true.
   else
      phys_deepconv_pbl = .false.
   endif

   return

end function phys_deepconv_pbl

!===============================================================================

function phys_do_flux_avg()

   logical :: phys_do_flux_avg
   !----------------------------------------------------------------------

   phys_do_flux_avg = .false.
   if (srf_flux_avg == 1) phys_do_flux_avg = .true.

end function phys_do_flux_avg

!===============================================================================
end module phys_control
