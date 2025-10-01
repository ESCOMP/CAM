module vertical_diffusion

!----------------------------------------------------------------------------------------------------- !
! Module to compute vertical diffusion of momentum,  moisture, trace constituents                      !
! and static energy. Separate modules compute                                                          !
!   1. stresses associated with turbulent flow over orography                                          !
!      ( turbulent mountain stress )                                                                   !
!   2. eddy diffusivities, including nonlocal tranport terms                                           !
!   3. molecular diffusivities                                                                         !
! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
! differencing the diffused and initial states.                                                        !
!                                                                                                      !
! Calling sequence:                                                                                    !
!                                                                                                      !
!  vertical_diffusion_init      Initializes vertical diffustion constants and modules                  !
!        init_molec_diff        Initializes molecular diffusivity module                               !
!        init_eddy_diff         Initializes eddy diffusivity module (includes PBL)                     !
!        init_tms               Initializes turbulent mountain stress module                           !
!        init_vdiff             Initializes diffusion solver module                                    !
!  vertical_diffusion_ts_init   Time step initialization (only used for upper boundary condition)      !
!  vertical_diffusion_tend      Computes vertical diffusion tendencies                                 !
!        compute_tms            Computes turbulent mountain stresses                                   !
!        compute_eddy_diff      Computes eddy diffusivities and countergradient terms                  !
!        compute_vdiff          Solves vertical diffusion equations, including molecular diffusivities !
!                                                                                                      !
!----------------------------------------------------------------------------------------------------- !
! Some notes on refactoring changes made in 2015, which were not quite finished.                       !
!                                                                                                      !
!      - eddy_diff_tend should really only have state, pbuf, and cam_in as inputs. The process of      !
!        removing these arguments, and referring to pbuf fields instead, is not complete.              !
!                                                                                                      !
!      - compute_vdiff was intended to be split up into three components:                              !
!                                                                                                      !
!         1. Diffusion of winds and heat ("U", "V", and "S" in the fieldlist object).                  !
!                                                                                                      !
!         2. Turbulent diffusion of a single constituent                                               !
!                                                                                                      !
!         3. Molecular diffusion of a single constituent                                               !
!                                                                                                      !
!        This reorganization would allow the three resulting functions to each use a simpler interface !
!        than the current combined version, and possibly also remove the need to use the fieldlist     !
!        object at all.                                                                                !
!                                                                                                      !
!---------------------------Code history-------------------------------------------------------------- !
! J. Rosinski : Jun. 1992                                                                              !
! J. McCaa    : Sep. 2004                                                                              !
! S. Park     : Aug. 2006, Dec. 2008. Jan. 2010                                                        !
!----------------------------------------------------------------------------------------------------- !

use shr_kind_mod,     only : r8 => shr_kind_r8
use ppgrid,           only : pcols, pver, pverp
use constituents,     only : pcnst
use cam_abortutils,   only : endrun
use error_messages,   only : handle_errmsg
use physconst,        only :          &
     cpair  , &     ! Specific heat of dry air
     gravit , &     ! Acceleration due to gravity
     rair   , &     ! Gas constant for dry air
     zvir   , &     ! rh2o/rair - 1
     latvap , &     ! Latent heat of vaporization
     latice , &     ! Latent heat of fusion
     karman , &     ! von Karman constant
     mwdry  , &     ! Molecular weight of dry air
     avogad         ! Avogadro's number
use cam_history,      only : fieldname_len
use perf_mod
use cam_logfile,      only : iulog
use ref_pres,         only : do_molec_diff, nbot_molec
use phys_control,     only : phys_getopts
use time_manager,     only : is_first_step
implicit none
private
save

! ----------------- !
! Public interfaces !
! ----------------- !

public vd_readnl
public vd_register                                   ! Register multi-time-level variables with physics buffer
public vertical_diffusion_init                       ! Initialization
public vertical_diffusion_ts_init                    ! Time step initialization (only used for upper boundary condition)
public vertical_diffusion_tend                       ! Full vertical diffusion routine

! ------------ !
! Private data !
! ------------ !

character(len=16)    :: eddy_scheme                  ! Default set in phys_control.F90, use namelist to change
!     'HB'       = Holtslag and Boville (default)
!     'HBR'      = Holtslag and Boville and Rash
!     'diag_TKE' = Bretherton and Park ( UW Moist Turbulence Scheme )
!     'CLUBB_SGS'= CLUBB;
!       in this case, this module will only be used for:
!       1) applying non-water vapor constituent indices;
!       2) applying HB scheme above CLUBB (when do_hb_above_clubb = .true.)
!       many PBL diagnostics are suppressed in this module via the logical is_clubb_scheme = .true.
!
logical, parameter   :: wstarent = .true.            ! Use wstar (.true.) or TKE (.false.) entrainment closure
! ( when 'diag_TKE' scheme is selected )
logical              :: do_pseudocon_diff = .false.  ! If .true., do pseudo-conservative variables diffusion

character(len=16)    :: shallow_scheme               ! Shallow convection scheme

logical              :: do_diffusion_const_dry(pcnst)! Do vertical diffusion (as a dry mixing ratio) for this constituent?
logical              :: do_diffusion_const_wet(pcnst)! Do vertical diffusion (as a wet mixing ratio) for this constituent?
logical              :: do_molecular_diffusion_const(pcnst) ! Do molecular diffusion for this constituent?
                                                            ! In compute_vdiff, both do_diffusion_const and molecular_diffusion_const
                                                            ! have to be true for molec diff to happen

integer              :: tke_idx, kvh_idx, kvm_idx    ! TKE and eddy diffusivity indices for fields in the physics buffer
integer              :: kvt_idx                      ! Index for kinematic molecular conductivity
integer              :: tauresx_idx, tauresy_idx     ! Redisual stress for implicit surface stress

character(len=fieldname_len) :: vdiffnam(pcnst)      ! Names of vertical diffusion tendencies
integer              :: ixcldice, ixcldliq           ! Constituent indices for cloud liquid and ice water
integer              :: ixnumice, ixnumliq
integer              :: ixq

integer              :: pblh_idx, tpert_idx, qpert_idx

! pbuf fields for unicon
integer              :: qtl_flx_idx  = -1            ! for use in cloud macrophysics when UNICON is on
integer              :: qti_flx_idx  = -1            ! for use in cloud macrophysics when UNICON is on

! pbuf fields for tms
integer              :: ksrftms_idx  = -1
integer              :: tautmsx_idx  = -1
integer              :: tautmsy_idx  = -1

! pbuf fields for blj (Beljaars)
integer              :: dragblj_idx  = -1
integer              :: taubljx_idx  = -1
integer              :: taubljy_idx  = -1

! pbuf field for clubb top above which HB (Holtslag Boville) scheme may be enabled
integer              :: clubbtop_idx = -1

logical              :: diff_cnsrv_mass_check        ! do mass conservation check
logical              :: do_iss                       ! switch for implicit turbulent surface stress
logical              :: is_clubb_scheme = .false.
logical              :: waccmx_mode = .false.
logical              :: do_hb_above_clubb = .false.

contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
subroutine vd_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: masterproc, masterprocid, mpi_logical, mpicom
  use shr_log_mod,     only: errMsg => shr_log_errMsg
  use trb_mtn_stress_cam, only: trb_mtn_stress_readnl
  use beljaars_drag_cam, only: beljaars_drag_readnl
  use eddy_diff_cam,   only: eddy_diff_readnl

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'vd_readnl'

  namelist /vert_diff_nl/ diff_cnsrv_mass_check, do_iss
  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'vert_diff_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, vert_diff_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(diff_cnsrv_mass_check, 1, mpi_logical, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(do_iss,                1, mpi_logical, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")

  ! Get eddy_scheme setting from phys_control.
  call phys_getopts( eddy_scheme_out          =          eddy_scheme, &
       shallow_scheme_out       =       shallow_scheme )

  ! TMS reads its own namelist.
  call trb_mtn_stress_readnl(nlfile)

  ! Beljaars reads its own namelist.
  call beljaars_drag_readnl(nlfile)

  if (eddy_scheme == 'diag_TKE') call eddy_diff_readnl(nlfile)

end subroutine vd_readnl

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine vd_register()

  !------------------------------------------------ !
  ! Register physics buffer fields and constituents !
  !------------------------------------------------ !

  use physics_buffer,      only : pbuf_add_field, dtype_r8
  use trb_mtn_stress_cam,  only : trb_mtn_stress_register
  use beljaars_drag_cam,   only : beljaars_drag_register
  use eddy_diff_cam,       only : eddy_diff_register

  ! Add fields to physics buffer

  ! kvt is used by gw_drag.  only needs physpkg scope.
  call pbuf_add_field('kvt', 'physpkg', dtype_r8, (/pcols,pverp/), kvt_idx) ! molecular_kinematic_temperature_conductivity_at_interfaces


  if (eddy_scheme /= 'CLUBB_SGS') then
     call pbuf_add_field('kvh',      'global', dtype_r8, (/pcols, pverp/), kvh_idx) ! eddy_heat_diffusivity_at_interfaces
  end if

  call pbuf_add_field('kvm',      'global', dtype_r8, (/pcols, pverp/), kvm_idx ) ! eddy_momentum_diffusivity_at_interfaces
  call pbuf_add_field('pblh',     'global', dtype_r8, (/pcols/),        pblh_idx) ! atmosphere_boundary_layer_thickness
  call pbuf_add_field('tke',      'global', dtype_r8, (/pcols, pverp/), tke_idx) ! specific_turbulent_kinetic_energy_at_interfaces

  call pbuf_add_field('tauresx',  'global', dtype_r8, (/pcols/),        tauresx_idx) ! eastward_reserved_stress_at_surface_on_previous_timestep
  call pbuf_add_field('tauresy',  'global', dtype_r8, (/pcols/),        tauresy_idx) ! northward_reserved_stress_at_surface_on_previous_timestep

  call pbuf_add_field('tpert', 'global', dtype_r8, (/pcols/),                       tpert_idx) ! convective_temperature_perturbation_due_to_pbl_eddies
  call pbuf_add_field('qpert', 'global', dtype_r8, (/pcols/),                       qpert_idx) ! convective_water_vapor_wrt_moist_air_and_condensed_water_perturbation_due_to_pbl_eddies

  if (trim(shallow_scheme) == 'UNICON') then
     call pbuf_add_field('qtl_flx',  'global', dtype_r8, (/pcols, pverp/), qtl_flx_idx)
     call pbuf_add_field('qti_flx',  'global', dtype_r8, (/pcols, pverp/), qti_flx_idx)
  end if

  ! diag_TKE fields
  if (eddy_scheme == 'diag_TKE') then
     call eddy_diff_register()
  end if

  ! TMS fields
  call trb_mtn_stress_register()

  ! Beljaars fields
  call beljaars_drag_register()

end subroutine vd_register

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine vertical_diffusion_init(pbuf2d)

  !------------------------------------------------------------------!
  ! Initialization of time independent fields for vertical diffusion !
  ! Calls initialization routines for subsidiary modules             !
  !----------------------------------------------------------------- !

  use cam_history,       only : addfld, add_default, horiz_only
  use cam_history,       only : register_vector_field
  use eddy_diff_cam,     only : eddy_diff_init

  use holtslag_boville_diff, only: holtslag_boville_diff_init
  use vertical_diffusion_sponge_layer, only: vertical_diffusion_sponge_layer_init

  use holtslag_boville_diff_interstitials, only: hb_diff_set_vertical_diffusion_top_init

  use molec_diff,        only : init_molec_diff
  use diffusion_solver_cam, only : init_vdiff
  use constituents,      only : cnst_get_ind, cnst_name, cnst_get_molec_byind, cnst_ndropmixed
  use spmd_utils,        only : masterproc
  use ref_pres,          only : pref_mid
  use physics_buffer,    only : pbuf_set_field, pbuf_get_index, physics_buffer_desc
  use trb_mtn_stress_cam,only : trb_mtn_stress_init
  use beljaars_drag_cam, only : beljaars_drag_init
  use upper_bc,          only : ubc_init
  use phys_control,      only : waccmx_is, fv_am_correction
  use ref_pres,          only : ptop_ref

  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  character(128) :: errstring   ! Error status for init_vdiff
  integer        :: ntop_eddy   ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
  integer        :: nbot_eddy   ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
  integer        :: k           ! Vertical loop index

  logical :: history_amwg                 ! output the variables used by the AMWG diag package
  logical :: history_eddy                 ! output the eddy variables
  logical :: history_budget               ! Output tendencies and state variables for CAM4 T, qv, ql, qi
  integer :: history_budget_histfile_num  ! output history file number for budget fields
  logical :: history_waccm                ! output variables of interest for WACCM runs

  character(len=512)   :: errmsg
  integer              :: errflg

  !
  ! add sponge layer vertical diffusion
  !
  call vertical_diffusion_sponge_layer_init( &
    amIRoot  = masterproc, &
    iulog    = iulog, &
    ptop_ref = ptop_ref, &
    errmsg   = errmsg, &
    errflg   = errflg)

  ! Check to see if WACCM-X is on (currently we don't care whether the
  ! ionosphere is on or not, since this neutral diffusion code is the
  ! same either way).
  waccmx_mode = waccmx_is('ionosphere') .or. waccmx_is('neutral')

  ! ----------------------------------------------------------------- !
  ! Get indices of cloud liquid and ice within the constituents array !
  ! ----------------------------------------------------------------- !

  call cnst_get_ind( 'CLDLIQ', ixcldliq )
  call cnst_get_ind( 'CLDICE', ixcldice )
  ! These are optional; with the CAM4 microphysics, there are no number
  ! constituents.
  call cnst_get_ind( 'NUMLIQ', ixnumliq, abort=.false. )
  call cnst_get_ind( 'NUMICE', ixnumice, abort=.false. )

  call cnst_get_ind( 'Q', ixq ) ! water vapor index in const array.

  ! Initialize upper boundary condition module
  call ubc_init()

  ! ---------------------------------------------------------------------------------------- !
  ! Initialize molecular diffusivity module                                                  !
  ! Note that computing molecular diffusivities is a trivial expense, but constituent        !
  ! diffusivities depend on their molecular weights. Decomposing the diffusion matrix        !
  ! for each constituent is a needless expense unless the diffusivity is significant.        !
  ! ---------------------------------------------------------------------------------------- !

  !----------------------------------------------------------------------------------------
  ! Initialize molecular diffusion and get top and bottom molecular diffusion limits
  !----------------------------------------------------------------------------------------

  if( do_molec_diff ) then
     call init_molec_diff( r8, pcnst, mwdry, avogad, &
          errstring)

     call handle_errmsg(errstring, subname="init_molec_diff")

     call addfld( 'TTPXMLC', horiz_only, 'A', 'K/S', 'Top interf. temp. flux: molec. viscosity' )
     if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NBOT_MOLEC =', nbot_molec
  end if

  ! ---------------------------------- !
  ! Initialize eddy diffusivity module !
  ! ---------------------------------- !

  ! ntop_eddy must be 1 or <= nbot_molec
  ! Currently, it is always 1 except for WACCM-X.
  call hb_diff_set_vertical_diffusion_top_init(ntop_eddy=ntop_eddy, errmsg=errmsg, errflg=errflg)
  nbot_eddy  = pver

  if (masterproc) write(iulog, fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy

  call phys_getopts(do_hb_above_clubb_out=do_hb_above_clubb)

  select case ( eddy_scheme )
  case ( 'diag_TKE' )
     if( masterproc ) write(iulog,*) &
          'vertical_diffusion_init: eddy_diffusivity scheme: UW Moist Turbulence Scheme by Bretherton and Park'
     call eddy_diff_init(pbuf2d, ntop_eddy, nbot_eddy)
  case ( 'HB', 'HBR')
     if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'

     call holtslag_boville_diff_init( &
      amIRoot = masterproc, &
      iulog   = iulog, &
      pver    = pver, &
      pverp   = pverp, &
      karman  = karman, &
      pref_mid = pref_mid, &
      is_hbr_pbl_scheme = (eddy_scheme .eq. 'HBR'), &
      ntop_turb_in = ntop_eddy, &
      errmsg = errmsg, &
      errflg = errflg)

     if(errflg /= 0) then
        call endrun('holtslag_boville_diff_init: ' // errmsg)
     endif

     call addfld('HB_ri',      (/ 'lev' /),  'A', 'no',  'Richardson Number (HB Scheme), I' )
  case ( 'CLUBB_SGS' )
     is_clubb_scheme = .true.

     call holtslag_boville_diff_init( &
      amIRoot = masterproc, &
      iulog   = iulog, &
      pver    = pver, &
      pverp   = pverp, &
      karman  = karman, &
      pref_mid = pref_mid, &
      is_hbr_pbl_scheme = (eddy_scheme .eq. 'HBR'), &
      ntop_turb_in = ntop_eddy, &
      errmsg = errmsg, &
      errflg = errflg)

     if(errflg /= 0) then
        call endrun('holtslag_boville_diff_init: ' // errmsg)
     endif

     !
     ! run HB scheme where CLUBB is not active when running cam7 or cam6 physics
     ! else init_hb_diff is called just for diagnostic purposes
     !
     if (do_hb_above_clubb) then
       if( masterproc ) then
         write(iulog,*) 'vertical_diffusion_init: '
         write(iulog,*) 'eddy_diffusivity scheme where CLUBB is not active:  Holtslag and Boville'
       end if
       call addfld('HB_ri',      (/ 'lev' /),  'A', 'no',  'Richardson Number (HB Scheme), I' )
     end if
  end select

  ! ------------------------------------------- !
  ! Initialize turbulent mountain stress module !
  ! ------------------------------------------- !

  call trb_mtn_stress_init()

  ! ----------------------------------- !
  ! Initialize Beljaars SGO drag module !
  ! ----------------------------------- !

  call beljaars_drag_init()

  ! ---------------------------------- !
  ! Initialize diffusion solver module !
  ! ---------------------------------- !

  call init_vdiff(r8, iulog, rair, cpair, gravit, do_iss, fv_am_correction, errstring)
  call handle_errmsg(errstring, subname="init_vdiff")

  ! Set which fields will be diffused using dry or moist mixing ratios.
  ! All fields are diffused using moist mixing ratios by default.
  do_diffusion_const_dry(:) = .false.
  do_diffusion_const_wet(:) = .false.
  do_molecular_diffusion_const(:) = .false.

  const_loop: do k = 1, pcnst
    ! If constituent is treated in dropmixnuc (vertically mixed by ndrop activation process)
    ! then do not handle vertical diffusion in this module.
    if (cnst_ndropmixed(k)) then
      cycle const_loop
    endif

    ! Convert all constituents to wet before doing diffusion.
    do_diffusion_const_wet(k) = .true.

    ! Select constituents for molecular diffusion
    if (cnst_get_molec_byind(k) .eq. 'minor') then
       do_molecular_diffusion_const(k) = .true.
    endif
  end do const_loop

  ! ------------------------ !
  ! Diagnostic output fields !
  ! ------------------------ !

  do k = 1, pcnst
     vdiffnam(k) = 'VD'//cnst_name(k)
     if( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
     call addfld( vdiffnam(k), (/ 'lev' /), 'A', 'kg/kg/s', 'Vertical diffusion of '//cnst_name(k) )
  end do

  if (.not. is_clubb_scheme) then
     call addfld( 'PBLH'        , horiz_only    , 'A', 'm'      , 'PBL height'                                         )
     call addfld( 'QT'          , (/ 'lev' /)   , 'A', 'kg/kg'  , 'Total water mixing ratio'                           )
     call addfld( 'SL'          , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liquid water static energy'                         )
     call addfld( 'SLV'         , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liq wat virtual static energy'                      )
     call addfld( 'SLFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Liquid static energy flux'                          )
     call addfld( 'QTFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Total water flux'                                   )
     call addfld( 'TKE'         , (/ 'ilev' /)  , 'A', 'm2/s2'  , 'Turbulent Kinetic Energy'                           )
     call addfld( 'TPERT'       , horiz_only    , 'A', 'K'      , 'Perturbation temperature (eddies in PBL)'           )
     call addfld( 'QPERT'       , horiz_only    , 'A', 'kg/kg'  , 'Perturbation specific humidity (eddies in PBL)'     )

     call addfld( 'UFLX'        , (/ 'ilev' /)  , 'A', 'N/m2'   , 'Zonal momentum flux'                                )
     call addfld( 'VFLX'        , (/ 'ilev' /)  , 'A', 'N/m2'   , 'Meridional momentm flux'                            )
     call register_vector_field('UFLX', 'VFLX')

     ! For detailed analysis of UW moist turbulence scheme
     call addfld( 'qt_pre_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qt_prePBL'          )
     call addfld( 'sl_pre_PBL',   (/ 'lev' /)   , 'A', 'J/kg'   , 'sl_prePBL'          )
     call addfld( 'slv_pre_PBL',  (/ 'lev' /)   , 'A', 'J/kg'   , 'slv_prePBL'         )
     call addfld( 'u_pre_PBL',    (/ 'lev' /)   , 'A', 'm/s'    , 'u_prePBL'           )
     call addfld( 'v_pre_PBL',    (/ 'lev' /)   , 'A', 'm/s'    , 'v_prePBL'           )
     call addfld( 'qv_pre_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qv_prePBL'          )
     call addfld( 'ql_pre_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'ql_prePBL'          )
     call addfld( 'qi_pre_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qi_prePBL'          )
     call addfld( 't_pre_PBL',    (/ 'lev' /)   , 'A', 'K'      , 't_prePBL'           )
     call addfld( 'rh_pre_PBL',   (/ 'lev' /)   , 'A', '%'      , 'rh_prePBL'          )

     call addfld( 'qt_aft_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qt_afterPBL'        )
     call addfld( 'sl_aft_PBL',   (/ 'lev' /)   , 'A', 'J/kg'   , 'sl_afterPBL'        )
     call addfld( 'slv_aft_PBL',  (/ 'lev' /)   , 'A', 'J/kg'   , 'slv_afterPBL'       )
     call addfld( 'u_aft_PBL',    (/ 'lev' /)   , 'A', 'm/s'    , 'u_afterPBL'         )
     call addfld( 'v_aft_PBL',    (/ 'lev' /)   , 'A', 'm/s'    , 'v_afterPBL'         )
     call addfld( 'qv_aft_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qv_afterPBL'        )
     call addfld( 'ql_aft_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'ql_afterPBL'        )
     call addfld( 'qi_aft_PBL',   (/ 'lev' /)   , 'A', 'kg/kg'  , 'qi_afterPBL'        )
     call addfld( 't_aft_PBL',    (/ 'lev' /)   , 'A', 'K'      , 't_afterPBL'         )
     call addfld( 'rh_aft_PBL',   (/ 'lev' /)   , 'A', '%'      , 'rh_afterPBL'        )

     call addfld( 'slflx_PBL',    (/ 'ilev' /)  , 'A', 'J/m2/s' , 'sl flux by PBL'     )
     call addfld( 'qtflx_PBL',    (/ 'ilev' /)  , 'A', 'kg/m2/s', 'qt flux by PBL'     )
     call addfld( 'uflx_PBL',     (/ 'ilev' /)  , 'A', 'kg/m/s2', 'u flux by PBL'      )
     call addfld( 'vflx_PBL',     (/ 'ilev' /)  , 'A', 'kg/m/s2', 'v flux by PBL'      )

     call addfld( 'slflx_cg_PBL', (/ 'ilev' /)  , 'A', 'J/m2/s' , 'sl_cg flux by PBL'  )
     call addfld( 'qtflx_cg_PBL', (/ 'ilev' /)  , 'A', 'kg/m2/s', 'qt_cg flux by PBL'  )
     call addfld( 'uflx_cg_PBL',  (/ 'ilev' /)  , 'A', 'kg/m/s2', 'u_cg flux by PBL'   )
     call addfld( 'vflx_cg_PBL',  (/ 'ilev' /)  , 'A', 'kg/m/s2', 'v_cg flux by PBL'   )

     call addfld( 'qtten_PBL',    (/ 'lev' /)   , 'A', 'kg/kg/s', 'qt tendency by PBL' )
     call addfld( 'slten_PBL',    (/ 'lev' /)   , 'A', 'J/kg/s' , 'sl tendency by PBL' )
     call addfld( 'uten_PBL',     (/ 'lev' /)   , 'A', 'm/s2'   , 'u tendency by PBL'  )
     call addfld( 'vten_PBL',     (/ 'lev' /)   , 'A', 'm/s2'   , 'v tendency by PBL'  )
     call addfld( 'qvten_PBL',    (/ 'lev' /)   , 'A', 'kg/kg/s', 'qv tendency by PBL' )
     call addfld( 'qlten_PBL',    (/ 'lev' /)   , 'A', 'kg/kg/s', 'ql tendency by PBL' )
     call addfld( 'qiten_PBL',    (/ 'lev' /)   , 'A', 'kg/kg/s', 'qi tendency by PBL' )
     call addfld( 'tten_PBL',     (/ 'lev' /)   , 'A', 'K/s'    , 'T tendency by PBL'  )
     call addfld( 'rhten_PBL',    (/ 'lev' /)   , 'A', '%/s'    , 'RH tendency by PBL' )
  end if

  call addfld( 'USTAR'       , horiz_only    , 'A', 'm/s'    , 'Surface friction velocity'                          )
  call addfld( 'KVH'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion diffusivities (heat/moisture)'   )
  call addfld( 'KVM'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion diffusivities (momentum)'        )
  call addfld( 'KVT'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion kinematic molecular conductivity')
  call addfld( 'CGS'         , (/ 'ilev' /)  , 'A', 's/m2'   , 'Counter-gradient coeff on surface kinematic fluxes' )
  call addfld( 'DTVKE'       , (/ 'lev' /)   , 'A', 'K/s'    , 'dT/dt vertical diffusion KE dissipation'            )
  call addfld( 'DTV'         , (/ 'lev' /)   , 'A', 'K/s'    , 'T vertical diffusion'                               )
  call addfld( 'DUV'         , (/ 'lev' /)   , 'A', 'm/s2'   , 'U vertical diffusion'                               )
  call addfld( 'DVV'         , (/ 'lev' /)   , 'A', 'm/s2'   , 'V vertical diffusion'                               )

  call addfld ('ustar',horiz_only, 'A',     ' ',' ')
  call addfld ('obklen',horiz_only, 'A',    ' ',' ')

  ! ----------------------------
  ! determine default variables
  ! ----------------------------

  call phys_getopts( history_amwg_out = history_amwg, &
       history_eddy_out = history_eddy, &
       history_budget_out = history_budget, &
       history_budget_histfile_num_out = history_budget_histfile_num, &
       history_waccm_out = history_waccm)

  if (history_amwg) then
     call add_default(  vdiffnam(1), 1, ' ' )
     call add_default( 'DTV'       , 1, ' ' )
     if (.not. is_clubb_scheme) then
        call add_default( 'PBLH'      , 1, ' ' )
     end if
  endif

  if (history_eddy) then
     if (.not. is_clubb_scheme) then
        call add_default( 'UFLX    ', 1, ' ' )
        call add_default( 'VFLX    ', 1, ' ' )
     end if
  endif

  if( history_budget ) then
     call add_default( vdiffnam(ixcldliq), history_budget_histfile_num, ' ' )
     call add_default( vdiffnam(ixcldice), history_budget_histfile_num, ' ' )
     if( history_budget_histfile_num > 1 ) then
        call add_default(  vdiffnam(1), history_budget_histfile_num, ' ' )
        call add_default( 'DTV'       , history_budget_histfile_num, ' ' )
     end if
  end if

  if ( history_waccm ) then
     if (do_molec_diff) then
        call add_default ( 'TTPXMLC', 1, ' ' )
     end if
     call add_default( 'DUV'     , 1, ' ' )
     call add_default( 'DVV'     , 1, ' ' )
  end if
  ! ----------------------------


  ksrftms_idx = pbuf_get_index('ksrftms')
  tautmsx_idx = pbuf_get_index('tautmsx')
  tautmsy_idx = pbuf_get_index('tautmsy')

  dragblj_idx = pbuf_get_index('dragblj')
  taubljx_idx = pbuf_get_index('taubljx')
  taubljy_idx = pbuf_get_index('taubljy')

  if (eddy_scheme == 'CLUBB_SGS') then
     kvh_idx = pbuf_get_index('kvh')
  end if

  if (do_hb_above_clubb) then
     ! pbuf field denoting top of clubb
     clubbtop_idx = pbuf_get_index('clubbtop')
  end if

  ! Initialization of some pbuf fields
  if (is_first_step()) then
     ! Initialization of pbuf fields tke, kvh, kvm are done in phys_inidat
     call pbuf_set_field(pbuf2d, tauresx_idx,  0.0_r8)
     call pbuf_set_field(pbuf2d, tauresy_idx,  0.0_r8)
     if (trim(shallow_scheme) == 'UNICON') then
        call pbuf_set_field(pbuf2d, qtl_flx_idx,  0.0_r8)
        call pbuf_set_field(pbuf2d, qti_flx_idx,  0.0_r8)
     end if
  end if
end subroutine vertical_diffusion_init

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine vertical_diffusion_ts_init( pbuf2d, state )

  !-------------------------------------------------------------- !
  ! Timestep dependent setting,                                   !
  ! At present only invokes upper bc code                         !
  !-------------------------------------------------------------- !
  use upper_bc,       only : ubc_timestep_init
  use physics_types , only : physics_state
  use ppgrid        , only : begchunk, endchunk

  use physics_buffer, only : physics_buffer_desc

  type(physics_state), intent(in) :: state(begchunk:endchunk)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  call ubc_timestep_init( pbuf2d, state)

end subroutine vertical_diffusion_ts_init

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine vertical_diffusion_tend( &
     ztodt    , state    , cam_in,          &
     ustar    , obklen   , ptend    , &
     cldn     , pbuf)
  !---------------------------------------------------- !
  ! This is an interface routine for vertical diffusion !
  !---------------------------------------------------- !
  use physics_buffer,     only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
  use physics_types,      only : physics_state, physics_ptend, physics_ptend_init

  ! flag for FV correction
  use phys_control,       only: fv_am_correction
  ! flag for Beljaars
  use beljaars_drag_cam,    only: do_beljaars

  use camsrfexch,           only : cam_in_t
  use cam_history,          only : outfld

  use eddy_diff_cam,        only : eddy_diff_tend


  ! CCPP-ized HB scheme
  use holtslag_boville_diff, only: hb_pbl_independent_coefficients_run
  use holtslag_boville_diff, only: hb_pbl_dependent_coefficients_run
  use holtslag_boville_diff, only: hb_diff_exchange_coefficients_run

  ! CCPP-ized HB (free atmosphere) scheme
  use holtslag_boville_diff, only: hb_diff_free_atm_exchange_coefficients_run

  ! CCPP-ized sponge layer logic
  use vertical_diffusion_sponge_layer, only: vertical_diffusion_sponge_layer_run

  ! CCPP-ized vertical diffusion solver (for non-WACCM-X use)
  ! to replace compute_vdiff
  ! and interstitials that have been CCPP-ized
  use holtslag_boville_diff_interstitials, only: hb_diff_prepare_vertical_diffusion_inputs_run
  use holtslag_boville_diff_interstitials, only: hb_free_atm_diff_prepare_vertical_diffusion_inputs_run
  use diffusion_solver,     only: vertical_diffusion_interpolate_to_interfaces_run
  use diffusion_solver,     only: implicit_surface_stress_add_drag_coefficient_run
  use diffusion_stubs,      only: turbulent_mountain_stress_add_drag_coefficient_run
  use diffusion_solver,     only: vertical_diffusion_wind_damping_rate_run
  use diffusion_stubs,      only: beljaars_add_wind_damping_rate_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_horizontal_momentum_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_dry_static_energy_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_tracers_run
  use diffusion_solver,     only: vertical_diffusion_tendencies_run
  use ccpp_constituent_prop_mod, only: ccpp_const_props

  ! Diagnostic utility subroutine shared with CCPP-ized/SIMA
  use vertical_diffusion_diagnostic_utils, only: vertical_diffusion_diagnostic_profiles

  use diffusion_solver,     only: vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run
  use molec_diff,           only : compute_molec_diff, vd_lu_qdecomp
  use constituents,         only : qmincg, qmin, cnst_type

  ! Legacy version of diffusion solver supporting molecular diffusion for WACCM.
  use diffusion_solver_cam, only : compute_vdiff
  use air_composition,      only : cpairv, rairv !Needed for calculation of upward H flux
  use air_composition,      only : mbarv
  use time_manager,         only : get_nstep
  use constituents,         only : cnst_get_type_byind, cnst_name, &
                                   cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx, cnst_ndropmixed
  use physconst,            only : pi
  use atmos_phys_pbl_utils, only: calc_virtual_temperature, calc_ideal_gas_rrho, calc_friction_velocity,                             &
                                  calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, calc_kinematic_buoyancy_flux, &
                                  calc_obukhov_length
  use upper_bc,             only : ubc_get_vals, ubc_fixed_temp
  use upper_bc,             only : ubc_get_flxs
  use coords_1d,            only : Coords1D
  use phys_control,         only : cam_physpkg_is

  ! --------------- !
  ! Input Arguments !
  ! --------------- !

  type(physics_state), intent(inout) :: state                     ! Physics state variables
  type(cam_in_t),      intent(in)    :: cam_in                    ! Surface inputs

  real(r8),            intent(in)    :: ztodt                     ! 2 delta-t [ s ]
  real(r8),            intent(in)    :: cldn(pcols,pver)          ! New stratus fraction [ fraction ]

  ! ---------------------- !
  ! Input-Output Arguments !
  ! ---------------------- !

  type(physics_ptend), intent(out) :: ptend                       ! Individual parameterization tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)

  ! ---------------- !
  ! Output Arguments !
  ! ---------------- !

  real(r8),            intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
  real(r8),            intent(out)   :: obklen(pcols)             ! Obukhov length [ m ]

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  character(128) :: errstring                                     ! Error status for compute_vdiff

  integer  :: lchnk                                               ! Chunk identifier
  integer  :: ncol                                                ! Number of atmospheric columns
  integer  :: i, k, l, m                                          ! column, level, constituent indices

  real(r8) :: dtk(pcols,pver)                                     ! T tendency from KE dissipation
  real(r8), pointer   :: tke(:,:)                                 ! Turbulent kinetic energy [ m2/s2 ]

  real(r8), pointer   :: qtl_flx(:,:)                             ! overbar(w'qtl') where qtl = qv + ql
  real(r8), pointer   :: qti_flx(:,:)                             ! overbar(w'qti') where qti = qv + qi

  real(r8) :: cgs(pcols,pverp)                                    ! Counter-gradient star  [ cg/flux ]
  real(r8) :: cgh(pcols,pverp)                                    ! Counter-gradient term for heat
  real(r8) :: rztodt                                              ! 1./ztodt [ 1/s ]
  real(r8), pointer :: ksrftms(:)                                 ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
  real(r8), pointer :: tautmsx(:)                                 ! U component of turbulent mountain stress [ N/m2 ]
  real(r8), pointer :: tautmsy(:)                                 ! V component of turbulent mountain stress [ N/m2 ]
  real(r8) :: tautotx(pcols)                                      ! U component of total surface stress [ N/m2 ]
  real(r8) :: tautoty(pcols)                                      ! V component of total surface stress [ N/m2 ]

  real(r8), pointer :: dragblj(:,:)                               ! Beljaars SGO form drag profile [ 1/s ]
  real(r8), pointer :: taubljx(:)                                 ! U component of turbulent mountain stress [ N/m2 ]
  real(r8), pointer :: taubljy(:)                                 ! V component of turbulent mountain stress [ N/m2 ]

  real(r8), pointer :: kvh_in(:,:)                                ! kvh from previous timestep [ m2/s ]
  real(r8), pointer :: kvm_in(:,:)                                ! kvm from previous timestep [ m2/s ]
  real(r8), pointer :: kvt(:,:)                                   ! Molecular kinematic conductivity for temperature [  ]
  real(r8) :: kvq(pcols,pverp)                                    ! Eddy diffusivity for constituents [ m2/s ]
  real(r8) :: kvh(pcols,pverp)                                    ! Eddy diffusivity for heat [ m2/s ]
  real(r8) :: kvm(pcols,pverp)                                    ! Eddy diffusivity for momentum [ m2/s ]
  real(r8) :: kvm_temp(pcols,pverp)                               ! Dummy eddy diffusivity for momentum (unused) [ m2/s ]
  real(r8) :: dtk_temp(pcols,pverp)                               ! Unused output from second compute_vdiff call
  real(r8) :: tautmsx_temp(pcols)                                 ! Unused output from second compute_vdiff call
  real(r8) :: tautmsy_temp(pcols)                                 ! Unused output from second compute_vdiff call
  real(r8) :: topflx_temp(pcols)                                  ! Unused output from second compute_vdiff call
  real(r8) :: sprod(pcols,pverp)                                  ! Shear production of tke [ m2/s3 ]
  real(r8) :: sfi(pcols,pverp)                                    ! Saturation fraction at interfaces [ fraction ]
  real(r8) :: sl(pcols,pver)
  real(r8) :: qt(pcols,pver)
  real(r8) :: slv(pcols,pver)
  real(r8) :: sl_prePBL(pcols,pver)
  real(r8) :: qt_prePBL(pcols,pver)
  real(r8) :: slv_prePBL(pcols,pver)
  real(r8) :: slten(pcols,pver)
  real(r8) :: qtten(pcols,pver)
  real(r8) :: slflx(pcols,pverp)
  real(r8) :: qtflx(pcols,pverp)
  real(r8) :: uflx(pcols,pverp)
  real(r8) :: vflx(pcols,pverp)
  real(r8) :: slflx_cg(pcols,pverp)
  real(r8) :: qtflx_cg(pcols,pverp)
  real(r8) :: uflx_cg(pcols,pverp)
  real(r8) :: vflx_cg(pcols,pverp)
  real(r8) :: th(pcols,pver)                                      ! Potential temperature
  real(r8) :: topflx(pcols)                                       ! Molecular heat flux at top interface
  real(r8) :: rhoair

  real(r8) :: ri(pcols,pver)                                      ! richardson number (HB output)

  ! for obklen calculation outside HB
  real(r8) :: thvs(pcols)                                         ! Virtual potential temperature at surface
  real(r8) :: rrho(pcols)                                         ! Reciprocal of density at surface
  real(r8) :: khfs(pcols)                                         ! sfc kinematic heat flux [K m s-1]
  real(r8) :: kqfs(pcols)                                         ! sfc kinematic water vapor flux [kg kg-1 m s-1]
  real(r8) :: kbfs(pcols)                                         ! sfc kinematic buoyancy flux [m^2/s^3]

  real(r8) :: ftem_prePBL(pcols,pver)                             ! Saturation vapor pressure before PBL
  real(r8) :: ftem_aftPBL(pcols,pver)                             ! Saturation vapor pressure after PBL
  real(r8) :: t_aftPBL(pcols,pver)                                ! Temperature after PBL diffusion
  real(r8) :: tten(pcols,pver)                                    ! Temperature tendency by PBL diffusion
  real(r8) :: rhten(pcols,pver)                                   ! RH tendency by PBL diffusion
  real(r8) :: qv_aft_PBL(pcols,pver)                              ! qv after PBL diffusion
  real(r8) :: ql_aft_PBL(pcols,pver)                              ! ql after PBL diffusion
  real(r8) :: qi_aft_PBL(pcols,pver)                              ! qi after PBL diffusion
  real(r8) :: s_aft_PBL(pcols,pver)                               ! s after PBL diffusion
  real(r8) :: u_aft_PBL(pcols,pver)                               ! u after PBL diffusion
  real(r8) :: v_aft_PBL(pcols,pver)                               ! v after PBL diffusion
  real(r8) :: qv_pro(pcols,pver)
  real(r8) :: ql_pro(pcols,pver)
  real(r8) :: qi_pro(pcols,pver)
  real(r8) :: s_pro(pcols,pver)
  real(r8) :: t_pro(pcols,pver)
  real(r8), pointer :: tauresx(:)                                      ! Residual stress to be added in vdiff to correct
  real(r8), pointer :: tauresy(:)                                      ! for turb stress mismatch between sfc and atm accumulated.

  ! Interpolated interface values.
  real(r8) :: tint(pcols,pver+1)      ! Temperature [ K ]
  real(r8) :: rairi(pcols,pver+1)     ! Gas constant [ J/K/kg ]
  real(r8) :: rhoi(pcols,pver+1)      ! Density of air [ kg/m^3 ]
  real(r8) :: rhoi_dry(pcols,pver+1)  ! Density of air based on dry air pressure [ kg/m^3 ]
  real(r8) :: t_toai(pcols)           ! Temporary for temperature to use at interface above TOA [K]

  ! Upper boundary conditions
  real(r8) :: ubc_t(pcols)            ! Temperature [ K ]
  real(r8) :: ubc_mmr(pcols,pcnst)    ! Mixing ratios [ kg/kg ]
  real(r8) :: ubc_flux(pcols,pcnst)   ! Constituent upper boundary flux (kg/s/m^2)

  ! Pressure coordinates used by the solver.
  type(Coords1D) :: p
  type(Coords1D) :: p_dry

  real(r8), pointer :: tpert(:)
  real(r8), pointer :: qpert(:)
  real(r8), pointer :: pblh(:)

  real(r8) :: tmp1(pcols)                                         ! Temporary storage

  integer  :: nstep
  real(r8) :: sum1, sum2, sum3, pdelx
  real(r8) :: sflx

  ! Copy state so we can pass to intent(inout) routines that return
  ! new state instead of a tendency.
  real(r8) :: s_tmp(pcols,pver)
  real(r8) :: u_tmp(pcols,pver)
  real(r8) :: v_tmp(pcols,pver)
  real(r8) :: q_tmp(pcols,pver,pcnst)

  ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
  real(r8) :: kq_scal(pcols,pver+1)
  ! composition dependent mw_fac on interface level
  real(r8) :: mw_fac(pcols,pver+1,pcnst)

  ! Dry static energy top boundary condition.
  real(r8) :: dse_top(pcols)

  ! Copies of flux arrays used to zero out any parts that are applied
  ! elsewhere (e.g. by CLUBB).
  real(r8) :: taux(pcols)
  real(r8) :: tauy(pcols)
  real(r8) :: shflux(pcols)
  real(r8) :: cflux(pcols,pcnst)
  integer,  pointer :: clubbtop(:)   ! (pcols)
  real(r8) :: clubbtop_r(pcols)

  logical  :: lq(pcnst)

  ! Temporaries for CCPP-ized HB
  real(r8) :: s2(pcols,pver)  ! shear squared (HB output) [s-2]
  real(r8) :: thv(pcols,pver) ! virtual potential temperature [K]
  real(r8) :: wstar(pcols)    ! convective scale velocity [m s-1]
  real(r8) :: bge(pcols)      ! buoyancy gradient enhancement
  real(r8) :: q_wv_cflx(pcols)! water vapor surface upward flux for kinematic wv fluxes

  ! Temporaries for CCPP-ized diffusion solver
  real(r8) :: ksrf(pcols) ! total surface drag coefficient
  real(r8) :: tau_damp_rate(pcols, pver) ! wind damping rate
  real(r8) :: tautotx_ccpp(pcols)
  real(r8) :: tautoty_ccpp(pcols)
  real(r8) :: dpidz_sq(pcols, pverp) ! square of derivative of pressure with height
  logical  :: itaures
  character(len=64) :: scheme_name

  character(len=512)   :: errmsg
  integer              :: errflg

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  rztodt = 1._r8 / ztodt
  lchnk  = state%lchnk
  ncol   = state%ncol

  call pbuf_get_field(pbuf, tauresx_idx,  tauresx)
  call pbuf_get_field(pbuf, tauresy_idx,  tauresy)
  call pbuf_get_field(pbuf, tpert_idx,    tpert)
  call pbuf_get_field(pbuf, qpert_idx,    qpert)
  call pbuf_get_field(pbuf, pblh_idx,     pblh)

  ! Get upper boundary values
  call ubc_get_vals( state%lchnk, ncol, state%pint, state%zi, ubc_t, ubc_mmr )

  if(waccmx_mode) then
    call ubc_get_flxs( state%lchnk, ncol, state%pint, state%zi, state%t, state%q, state%phis, ubc_flux )
  endif

  ! For WACCM-X or fixed upper boundary condition temperature
  ! set temperature at TOA interface, otherwise use temperature at TOA
  if (waccmx_mode) then
     ! For WACCM-X, set ubc temperature to extrapolate from next two lower interface level temperatures
     ! the original formulation is:
     ! t_toai(:ncol) = 1.5_r8*tint(:ncol,2)-.5_r8*tint(:ncol,3)
     ! this appears to be:
     !               = tint(:ncol,2) + 0.5_r8*(tint(:ncol,2) - tint(:ncol,3))
     ! assuming that the extrapolated gradient is 1/2 of the temperature gradient between the lower interfaces
     ! because the interpolation will be done later in the CCPP-ized scheme, formulate this in terms of
     ! the temperature (at midpoints):
     t_toai(:ncol) = 1.5_r8*(state%t(:ncol,2)+state%t(:ncol,1))/2._r8-.5_r8*(state%t(:ncol,3)+state%t(:ncol,2))/2._r8
  else
    if(ubc_fixed_temp) then
      ! Fixed temperature at upper boundary condition
      t_toai(:ncol) = ubc_t(:ncol)
    else
      ! Default
      t_toai(:ncol) = state%t(:ncol,1)
    endif
  endif

  !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
  tint(:,:) = 0._r8
  rairi(:,:) = 0._r8
  rhoi(:,:) = 0._r8
  rhoi_dry(:,:) = 0._r8
  dpidz_sq(:,:) = 0._r8
  !REMOVECAM_END

  ! Interpolate t, rho (moist and dry), and set rairi
  ! to interfaces for use by vertical diffusion solver.
  ! The CCPP-ized subroutine does not support WACCM-X which relies on
  ! upper boundary condition values; these will be overwritten later.
  call vertical_diffusion_interpolate_to_interfaces_run( &
       ncol     = ncol, &
       pver     = pver, &
       pverp    = pverp, &
       gravit   = gravit, &
       rair     = rair, &
       rairv    = rairv(:ncol,:pver,state%lchnk), &
       ! only use constituent-dependent gas constant when in WACCM-X mode.
       flag_for_constituent_dependent_gas_constant = waccmx_mode, &
       t        = state%t(:ncol,:pver), &
       t_toai   = t_toai(:ncol), &
       pint     = state%pint(:ncol,:pverp), &
       pintdry  = state%pintdry(:ncol,:pverp), &
       ! below output
       ti       = tint(:ncol,:pverp), &
       rairi    = rairi(:ncol,:pverp), &
       rhoi     = rhoi(:ncol,:pverp), &
       rhoi_dry = rhoi_dry(:ncol,:pverp), &
       dpidz_sq = dpidz_sq(:ncol,:pverp), &
       errmsg   = errmsg, &
       errflg   = errflg)

  if(errflg /= 0) then
     call endrun('vertical_diffusion_interpolate_to_interfaces_run: ' // errmsg)
  endif

  ! Initialize total surface stresses
  ! these are used for HB diffusion scheme and later PBL diagnostics but
  ! not for the vertical diffusion solver, which uses surface stresses from the coupler
  ! or just zero (in the case of CLUBB)
  tautotx(:ncol) = cam_in%wsx(:ncol)
  tautoty(:ncol) = cam_in%wsy(:ncol)

  ! ---------------------------------------- !
  ! Add in turbulent mountain stress (CAM5)  !
  ! ---------------------------------------- !

  ! Consistent with the computation of 'normal' drag coefficient,
  ! the raw input (u,v) is used to compute 'ksrftms', not the provisionally-marched 'u,v'
  ! within the iteration loop of the PBL scheme.

  call pbuf_get_field(pbuf, ksrftms_idx, ksrftms)
  call pbuf_get_field(pbuf, tautmsx_idx, tautmsx)
  call pbuf_get_field(pbuf, tautmsy_idx, tautmsy)

  ! Add turbulent mountain stress to total surface stress for PBL scheme.
  tautotx(:ncol) = tautotx(:ncol) + tautmsx(:ncol)
  tautoty(:ncol) = tautoty(:ncol) + tautmsy(:ncol)

  ! ------------------------------------- !
  ! Add in Beljaars SGO form drag (CAM6+) !
  ! ------------------------------------- !

  call pbuf_get_field(pbuf, dragblj_idx, dragblj)
  call pbuf_get_field(pbuf, taubljx_idx, taubljx)
  call pbuf_get_field(pbuf, taubljy_idx, taubljy)

  ! Add Beljaars integrated drag to total surface stress for PBL scheme.
  tautotx(:ncol) = tautotx(:ncol) + taubljx(:ncol)
  tautoty(:ncol) = tautoty(:ncol) + taubljy(:ncol)

  ! -------------------------------------
  ! Preparation of HB/HB_free inputs to vertical diffusion
  ! -------------------------------------

  ! Use CCPPized interstitial schemes for setting
  ! the necessary inputs (taux,tauy,shflux,cflux) for vertical_diffusion_compute_run from the coupler
  ! as well as Coords1D (p) pressure coordinates used by the solver.

  !REMOVECAM: no longer needed when pcols no longer exists
  q_wv_cflx(:) = 0._r8
  !END REMOVECAM

  if(eddy_scheme == 'CLUBB_SGS') then
     ! If running CLUBB_SGS, use the hb_free_atm CCPPized interstitial.
     !
     ! In this case, vertical diffusion solver only applies constituent fluxes excluding water vapor
     ! (before CAM7); no fluxes are applied in CAM7.
     call hb_free_atm_diff_prepare_vertical_diffusion_inputs_run( &
          ncol              = ncol,                       &
          pverp             = pverp,                      &
          pcnst             = pcnst,                      &
          const_props       = ccpp_const_props,           &
          apply_nonwv_cflx  = (.not. cam_physpkg_is("cam7")), & ! does vertical diffusion apply ANY fluxes?
          cflx_from_coupler = cam_in%cflx(:ncol,:pcnst),  &
          pint              = state%pint(:ncol,:pverp),   &
          ! below output
          taux              = taux(:ncol),                & ! these are zero since handled by CLUBB.
          tauy              = tauy(:ncol),                & ! these are zero since handled by CLUBB.
          shflux            = shflux(:ncol),              & ! these are zero since handled by CLUBB.
          cflux             = cflux(:ncol,:pcnst),        & ! if apply_nonwv_cflx, contains non-wv. fluxes, otherwise 0
          itaures           = itaures,                    &
          p                 = p,                          &
          q_wv_cflx         = q_wv_cflx(:ncol),           & ! for use in HB for kinematic water vapor flux calc.
          errmsg            = errmsg,                     &
          errflg            = errflg)

     if(errflg /= 0) then
        call endrun('hb_free_atm_diff_prepare_vertical_diffusion_inputs_run: ' // errmsg)
     endif
  else
     call hb_diff_prepare_vertical_diffusion_inputs_run( &
          ncol               = ncol,                      &
          pverp              = pverp,                     &
          pcnst              = pcnst,                     &
          const_props        = ccpp_const_props,          &
          wsx_from_coupler   = cam_in%wsx(:ncol),         &
          wsy_from_coupler   = cam_in%wsy(:ncol),         &
          shf_from_coupler   = cam_in%shf(:ncol),         &
          cflx_from_coupler  = cam_in%cflx(:ncol,:pcnst), &
          pint               = state%pint(:ncol,:pverp),  &
          ! below output
          taux               = taux(:ncol),               &
          tauy               = tauy(:ncol),               &
          shflux             = shflux(:ncol),             &
          cflux              = cflux(:ncol,:pcnst),       &
          itaures            = itaures,                   &
          p                  = p,                         &
          q_wv_cflx          = q_wv_cflx(:ncol),          & ! for use in HB for kinematic water vapor flux calc.
          errmsg             = errmsg,                 &
          errflg             = errflg)

     if(errflg /= 0) then
        call endrun('hb_diff_prepare_vertical_diffusion_inputs_run: ' // errmsg)
     endif
  endif

  !----------------------------------------------------------------------- !
  !   Computation of eddy diffusivities - Select appropriate PBL scheme    !
  !----------------------------------------------------------------------- !
  call pbuf_get_field(pbuf, kvm_idx,  kvm_in)
  call pbuf_get_field(pbuf, kvh_idx,  kvh_in)
  call pbuf_get_field(pbuf, tke_idx,  tke)

  select case (eddy_scheme)
  case ( 'diag_TKE' )

     ! Get potential temperature.
     th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

     ! Set up pressure coordinates for solver calls.
     p = Coords1D(state%pint(:ncol,:))

     call eddy_diff_tend(state, pbuf, cam_in, &
          ztodt, p, tint, rhoi, cldn, wstarent, &
          kvm_in, kvh_in, ksrftms, dragblj, tauresx, tauresy, &
          rrho, ustar, pblh, kvm, kvh, kvq, cgh, cgs, tpert, qpert, &
          tke, sprod, sfi)

     ! The diag_TKE scheme does not calculate the Monin-Obukhov length, which is used in dry deposition calculations.
     ! Use the routines from pbl_utils to accomplish this. Assumes ustar and rrho have been set.
      thvs  (:ncol) = calc_virtual_temperature(th(:ncol,pver), state%q(:ncol,pver,1), zvir)

      khfs  (:ncol) = calc_kinematic_heat_flux(cam_in%shf(:ncol), rrho(:ncol), cpair)
      kqfs  (:ncol) = calc_kinematic_water_vapor_flux(cam_in%cflx(:ncol,1), rrho(:ncol))
      kbfs  (:ncol) = calc_kinematic_buoyancy_flux(khfs(:ncol), zvir, th(:ncol,pver), kqfs(:ncol))
      obklen(:ncol) = calc_obukhov_length(thvs(:ncol), ustar(:ncol), gravit, karman, kbfs(:ncol))

  case ( 'HB', 'HBR' )

     ! Modification : We may need to use 'taux' instead of 'tautotx' here, for
     !                consistency with the previous HB scheme.

     !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
     thv(:,:) = 0._r8
     ustar(:) = 0._r8
     khfs(:) = 0._r8
     kqfs(:) = 0._r8
     kbfs(:) = 0._r8
     obklen(:) = 0._r8
     ri(:,:) = 0._r8
     s2(:,:) = 0._r8
     !REMOVECAM_END
     ! call CCPP-ized HB scheme (piecewise)
     call hb_pbl_independent_coefficients_run( &
       ncol      = ncol,                     &
       pver      = pver,                     &
       zvir      = zvir,                     &
       rair      = rair,                     &
       cpair     = cpair,                    &
       gravit    = gravit,                   &
       karman    = karman,                   &
       exner     = state%exner(:ncol,:pver), &
       t         = state%t(:ncol,:pver),     &
       q_wv      = state%q(:ncol,:pver,ixq), &
       z         = state%zm(:ncol,:pver),    &
       pmid      = state%pmid(:ncol,:pver),  &
       u         = state%u(:ncol,:pver),     &
       v         = state%v(:ncol,:pver),     &
       taux      = tautotx(:ncol),           &
       tauy      = tautoty(:ncol),           &
       shflx     = cam_in%shf(:ncol),        &
       q_wv_flx  = q_wv_cflx(:ncol),         &
       ! Output variables
       thv       = thv(:ncol,:pver),         &
       ustar     = ustar(:ncol),             &
       khfs      = khfs(:ncol),              &
       kqfs      = kqfs(:ncol),              &
       kbfs      = kbfs(:ncol),              &
       obklen    = obklen(:ncol),            &
       s2        = s2(:ncol,:pver),          &
       ri        = ri(:ncol,:pver),          &
       errmsg    = errmsg,                   &
       errflg    = errflg)

     if(errflg /= 0) then
        call endrun('hb_pbl_independent_coefficients_run: ' // errmsg)
     endif

     !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
     pblh(:) = 0._r8
     wstar(:) = 0._r8
     bge(:) = 0._r8
     !REMOVECAM_END
     call hb_pbl_dependent_coefficients_run( &
       ncol      = ncol,                                      &
       pver      = pver,                                      &
       pverp     = pverp,                                     &
       gravit    = gravit,                                    &
       z         = state%zm(:ncol,:pver),                     &
       zi        = state%zi(:ncol,:pverp),                    &
       u         = state%u(:ncol,:pver),                      &
       v         = state%v(:ncol,:pver),                      &
       cldn      = cldn(:ncol,:pver),                         &
       ! Inputs from pbl_independent_coefficients
       thv       = thv(:ncol,:pver),                          &
       ustar     = ustar(:ncol),                              &
       kbfs      = kbfs(:ncol),                               &
       obklen    = obklen(:ncol),                             &
       ! Output variables
       pblh      = pblh(:ncol),                               &
       wstar     = wstar(:ncol),                              &
       bge       = bge(:ncol),                                &
       errmsg    = errmsg,                                    &
       errflg    = errflg)

     if(errflg /= 0) then
        call endrun('hb_pbl_dependent_coefficients_run: ' // errmsg)
     endif

     !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
     kvm(:,:) = 0._r8
     kvh(:,:) = 0._r8
     kvq(:,:) = 0._r8
     cgh(:,:) = 0._r8
     cgs(:,:) = 0._r8
     tpert(:) = 0._r8
     qpert(:) = 0._r8
     tke(:,:) = 0._r8
     !REMOVECAM_END
     call hb_diff_exchange_coefficients_run( &
       ncol      = ncol,                                      &
       pver      = pver,                                      &
       pverp     = pverp,                                     &
       karman    = karman,                                    &
       cpair     = cpair,                                     &
       z         = state%zm(:ncol,:pver),                     &
       is_hbr_pbl_scheme = (eddy_scheme .eq. 'HBR'),          &
       ! Input from hb_pbl_independent_coefficients
       kqfs      = kqfs(:ncol),                               &
       khfs      = khfs(:ncol),                               &
       kbfs      = kbfs(:ncol),                               &
       ustar     = ustar(:ncol),                              &
       obklen    = obklen(:ncol),                             &
       s2        = s2(:ncol,:pver),                           &
       ri        = ri(:ncol,:pver),                           &
       ! Input from hb_pbl_dependent_coefficients
       pblh      = pblh(:ncol),                               &
       wstar     = wstar(:ncol),                              &
       bge       = bge(:ncol),                                &
       ! Output variables
       kvm       = kvm(:ncol,:pverp),                         &
       kvh       = kvh(:ncol,:pverp),                         &
       kvq       = kvq(:ncol,:pverp),                         &
       cgh       = cgh(:ncol,:pverp),                         &
       cgs       = cgs(:ncol,:pverp),                         &
       tpert     = tpert(:ncol),                              &
       qpert     = qpert(:ncol),                              &
       tke       = tke(:ncol,:pverp),                         &
       errmsg    = errmsg,                                    &
       errflg    = errflg)

     if(errflg /= 0) then
        call endrun('hb_diff_exchange_coefficients_run: ' // errmsg)
     endif

     call outfld( 'HB_ri',          ri,         pcols,   lchnk )

  case ( 'CLUBB_SGS' )
    !
    ! run HB scheme where CLUBB is not active when running cam7
    !
    if (do_hb_above_clubb) then
      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      thv(:,:) = 0._r8
      ustar(:) = 0._r8
      khfs(:) = 0._r8
      kqfs(:) = 0._r8
      kbfs(:) = 0._r8
      obklen(:) = 0._r8
      s2(:,:) = 0._r8
      ri(:,:) = 0._r8
      !REMOVECAM_END
      ! call CCPP-ized HB scheme (piecewise)
      call hb_pbl_independent_coefficients_run( &
        ncol      = ncol,                     &
        pver      = pver,                     &
        zvir      = zvir,                     &
        rair      = rair,                     &
        cpair     = cpair,                    &
        gravit    = gravit,                   &
        karman    = karman,                   &
        exner     = state%exner(:ncol,:pver), &
        t         = state%t(:ncol,:pver),     &
        q_wv      = state%q(:ncol,:pver,1),   & ! NOTE: assumes wv at 1 (need to change to ixq?)
        z         = state%zm(:ncol,:pver),    &
        pmid      = state%pmid(:ncol,:pver),  &
        u         = state%u(:ncol,:pver),     &
        v         = state%v(:ncol,:pver),     &
        taux      = tautotx(:ncol),           &
        tauy      = tautoty(:ncol),           &
        shflx     = cam_in%shf(:ncol),        &
        q_wv_flx  = q_wv_cflx(:ncol),         &
        ! Output variables
        thv       = thv(:ncol,:pver),         &
        ustar     = ustar(:ncol),             &
        khfs      = khfs(:ncol),              &
        kqfs      = kqfs(:ncol),              &
        kbfs      = kbfs(:ncol),              &
        obklen    = obklen(:ncol),            &
        s2        = s2(:ncol,:pver),          &
        ri        = ri(:ncol,:pver),          &
        errmsg    = errmsg,                   &
        errflg    = errflg)

      if(errflg /= 0) then
         call endrun('hb_pbl_independent_coefficients_run: ' // errmsg)
      endif

      call pbuf_get_field(pbuf, clubbtop_idx, clubbtop)
      clubbtop_r = real(clubbtop, r8)

      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      kvm(:,:) = 0._r8
      kvh(:,:) = 0._r8
      kvq(:,:) = 0._r8
      cgh(:,:) = 0._r8
      cgs(:,:) = 0._r8
      !REMOVECAM_END
      call hb_diff_free_atm_exchange_coefficients_run( &
        ncol      = ncol,                                      &
        pver      = pver,                                      &
        pverp     = pverp,                                     &
        ! Input from hb_pbl_independent_coefficients
        s2        = s2(:ncol,:pver),                           &
        ri        = ri(:ncol,:pver),                           &
        ! Zero out HB below this level, where CLUBB is active:
        bottom_boundary = clubbtop_r(:ncol),                   &
        ! Output variables
        kvm       = kvm(:ncol,:pverp),                         &
        kvh       = kvh(:ncol,:pverp),                         &
        kvq       = kvq(:ncol,:pverp),                         &
        cgh       = cgh(:ncol,:pverp),                         &
        cgs       = cgs(:ncol,:pverp),                         &
        errmsg    = errmsg,                                    &
        errflg    = errflg)

      if(errflg /= 0) then
         call endrun('hb_diff_free_atm_exchange_coefficients_run: ' // errmsg)
      endif

      call outfld( 'HB_ri',          ri,         pcols,   lchnk )
    else
      ! CLUBB has only a bare-bones placeholder here. If using CLUBB, the
      ! PBL diffusion will happen before coupling, so vertical_diffusion
      ! is only handling other things, e.g. some boundary conditions, tms,
      ! and molecular diffusion.

      ! Get potential temperature.
      th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

      thvs  (:ncol) = calc_virtual_temperature(th(:ncol,pver), state%q(:ncol,pver,1), zvir)
      rrho  (:ncol) = calc_ideal_gas_rrho(rair, state%t(:ncol,pver), state%pmid(:ncol,pver))
      ustar (:ncol) = calc_friction_velocity(cam_in%wsx(:ncol), cam_in%wsy(:ncol), rrho(:ncol))
      khfs  (:ncol) = calc_kinematic_heat_flux(cam_in%shf(:ncol), rrho(:ncol), cpair)
      kqfs  (:ncol) = calc_kinematic_water_vapor_flux(cam_in%cflx(:ncol,1), rrho(:ncol))
      kbfs  (:ncol) = calc_kinematic_buoyancy_flux(khfs(:ncol), zvir, th(:ncol,pver), kqfs(:ncol))
      obklen(:ncol) = calc_obukhov_length(thvs(:ncol), ustar(:ncol), gravit, karman, kbfs(:ncol))

      ! These tendencies all applied elsewhere.
      kvm = 0._r8
      kvh = 0._r8
      kvq = 0._r8
      ! Not defined since PBL is not actually running here.
      cgh = 0._r8
      cgs = 0._r8
    end if
  end select

  ! Write diagnostic output after diffusion coefficients are calculated
  call outfld( 'ustar',   ustar(:), pcols, lchnk )
  call outfld( 'obklen', obklen(:), pcols, lchnk )

  !
  ! add sponge layer vertical diffusion
  !
  call vertical_diffusion_sponge_layer_run( &
    ncol   = ncol,   &
    pverp  = pverp,  &
    kvm    = kvm,    & ! in/out
    errmsg = errmsg, &
    errflg = errflg)

  if(errflg /= 0) then
     call endrun('vertical_diffusion_sponge_layer_run: ' // errmsg)
  endif

  ! kvh (in pbuf) is used by other physics parameterizations, and as an initial guess in compute_eddy_diff
  ! on the next timestep.  It is not updated by the compute_vdiff call below.
  call pbuf_set_field(pbuf, kvh_idx, kvh)

  ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
  ! The contributions for molecular diffusion made to kvm by the call to compute_vdiff below
  ! are not included in the pbuf as these are not needed in the initial guess by compute_eddy_diff.
  call pbuf_set_field(pbuf, kvm_idx, kvm)

  ! Get molecular_kinematic_temperature_conductivity_at_interfaces
  call pbuf_get_field(pbuf, kvt_idx, kvt)

  !------------------------------------ !
  !    Application of diffusivities     !
  !------------------------------------ !

  ! Set arrays from input state.
  q_tmp(:ncol,:,:) = state%q(:ncol,:,:)
  s_tmp(:ncol,:) = state%s(:ncol,:)
  u_tmp(:ncol,:) = state%u(:ncol,:)
  v_tmp(:ncol,:) = state%v(:ncol,:)

  ! --------------------------------------------------------------------------------- !
  ! Call the diffusivity solver and solve diffusion equation                          !
  ! The final two arguments are optional function references to                       !
  ! constituent-independent and constituent-dependent moleculuar diffusivity routines !
  ! --------------------------------------------------------------------------------- !

  ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
  !                separately print out as diagnostic output, because these are different from
  !                the explicit 'tautotx, tautoty' computed above.
  ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

  if(.not. do_molec_diff) then
    ! Dry static energy top boundary is zero if no molecular diffusion
    ! Note: this is an artifact to the way the vertical diffusion calculations are done;
    ! it does not represent physical reality.
    dse_top(:ncol) = 0._r8

    ksrf(:) = 0._r8
    tau_damp_rate(:,:) = 0._r8
    tautotx_ccpp(:) = 0._r8
    tautoty_ccpp(:) = 0._r8

    ! Calculate surface drag rate
    call implicit_surface_stress_add_drag_coefficient_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         do_iss          = do_iss,                       &
         taux            = taux(:ncol),                  &
         tauy            = tauy(:ncol),                  &
         u0              = state%u(:ncol,:pver),         &
         v0              = state%v(:ncol,:pver),         &
         ! below input/output:
         ksrf            = ksrf(:ncol),                  &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('implicit_surface_stress_add_drag_coefficient_run: ' // errmsg)
    endif

    ! Add TMS surface drag rate
    call turbulent_mountain_stress_add_drag_coefficient_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         ksrftms         = ksrftms(:ncol),               &
         ! below input/output:
         ksrf            = ksrf(:ncol),                  &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('turbulent_mountain_stress_add_drag_coefficient_run: ' // errmsg)
    endif

    ! Based on the drag coefficients, calculate wind damping rates
    call vertical_diffusion_wind_damping_rate_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         gravit          = gravit,                       &
         p               = p,                            & ! Coords1D, pressure coordinates [Pa]
         ksrf            = ksrf(:ncol),                  &
         ! below output:
         tau_damp_rate   = tau_damp_rate(:ncol,:pver),   &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('vertical_diffusion_wind_damping_rate_run: ' // errmsg)
    endif

    ! Add Beljaars wind damping rate
    call beljaars_add_wind_damping_rate_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         dragblj         = dragblj(:ncol,:pver),         &
         ! below input/output:
         tau_damp_rate   = tau_damp_rate(:ncol,:pver),   &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('beljaars_add_wind_damping_rate_run: ' // errmsg)
    endif

    ! If molecular diffusion is not done, use the CCPP-ized subroutine
    call vertical_diffusion_diffuse_horizontal_momentum_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         pverp           = pverp,                        &
         dt              = ztodt,                        &
         rair            = rair,                         &
         gravit          = gravit,                       &
         do_iss          = do_iss,                       &
         am_correction   = fv_am_correction,             &
         itaures         = itaures,                      &
         t               = state%t(:ncol,:pver),         &
         p               = p,                            & ! Coords1D, pressure coordinates [Pa]
         rhoi            = rhoi(:ncol,:pverp),           &
         taux            = taux(:ncol),                  &
         tauy            = tauy(:ncol),                  &
         tau_damp_rate   = tau_damp_rate(:ncol,:pver),   & ! tau damp rate from above
         kvm             = kvm(:ncol,:pverp),            &
         ksrftms         = ksrftms(:ncol),               &
         dragblj         = dragblj(:ncol,:pver),         &
         dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist
         u0              = state%u(:ncol,:pver),         &
         v0              = state%v(:ncol,:pver),         &
         dse0            = state%s(:ncol,:pver),         &
         ! input/output
         tauresx         = tauresx(:ncol),               &
         tauresy         = tauresy(:ncol),               &
         ! below output
         u1              = u_tmp(:ncol,:pver),           &
         v1              = v_tmp(:ncol,:pver),           &
         dse1            = s_tmp(:ncol,:pver),           &
         dtk             = dtk(:ncol,:),                 &
         tautmsx         = tautmsx(:ncol),               &
         tautmsy         = tautmsy(:ncol),               &
         ! arguments for Beljaars
         do_beljaars     = do_beljaars,                  &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('vertical_diffusion_diffuse_horizontal_momentum_run: ' // errmsg)
    endif

    ! Diffuse dry static energy
    call vertical_diffusion_diffuse_dry_static_energy_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         dt              = ztodt,                        &
         gravit          = gravit,                       &
         p               = p,                            & ! Coords1D, pressure coordinates [Pa]
         rhoi            = rhoi(:ncol,:pverp),           &
         shflx           = shflux(:ncol),                &
         dse_top         = dse_top(:ncol),               & ! zero
         kvh             = kvh(:ncol,:pverp),            &
         cgh             = cgh(:ncol,:pverp),            &
         dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist
         ! input/output
         dse             = s_tmp(:ncol,:pver),           &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('vertical_diffusion_diffuse_dry_static_energy_run: ' // errmsg)
    endif

    ! Diffuse tracers
    call vertical_diffusion_diffuse_tracers_run( &
         ncol            = ncol,                         &
         pver            = pver,                         &
         ncnst           = pcnst,                        &
         dt              = ztodt,                        &
         rair            = rair,                         &
         gravit          = gravit,                       &
         do_diffusion_const = do_diffusion_const_wet,    & ! moist constituents to diffuse
         p               = p,                            & ! Coords1D, pressure coordinates [Pa]
         t               = state%t(:ncol,:pver),         &
         rhoi            = rhoi(:ncol,:pverp),           &
         cflx            = cflux(:ncol,:pcnst),          &
         kvh             = kvh(:ncol,:pverp),            &
         kvq             = kvq(:ncol,:pverp),            &
         cgs             = cgs(:ncol,:pverp),            &
         qmincg          = qmincg(:pcnst),               &
         dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist
         ! upper boundary conditions from ubc module
         ubc_mmr         = ubc_mmr(:ncol,:pcnst),        &
         cnst_fixed_ubc  = cnst_fixed_ubc(:pcnst),       &
         q0              = state%q(:ncol,:pver,:pcnst),  &
         q1              = q_tmp(:ncol,:pver,:pcnst),    &
         errmsg          = errmsg,                       &
         errflg          = errflg)

    if(errflg /= 0) then
       call endrun('vertical_diffusion_diffuse_tracers_run: ' // errmsg)
    endif
  else
     ! Molecular diffusion is active, use old compute_vdiff

     ! Top boundary condition for dry static energy if molecular diffusion is active
     ! but not in WACCM-X mode
     if (.not. waccmx_mode) then
        !REMOVECAM: no longer need this after pcols no longer exists
        dse_top(:) = 0._r8
        !END REMOVECAM

        call vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run( &
             ncol = ncol, &
             gravit = gravit, &
             cpairv = cpairv(:ncol,:,state%lchnk), &
             zi = state%zi(:ncol,:), &
             tint = tint(:ncol,:), &
             dse_top = dse_top(:ncol), &
             errmsg = errmsg, &
             errflg = errflg)

       if(errflg /= 0) then
          call endrun('vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run: ' // errmsg)
       endif
     else
        dse_top(:ncol) = 0._r8
     end if

     if( any(do_diffusion_const_wet) ) then
        call compute_molec_diff(state%lchnk, pcols, pver, pcnst, ncol, &
                kvm, kvt, tint, rhoi, kq_scal, cnst_mw, &
                mw_fac, nbot_molec)

        call compute_vdiff( &
             ncol            = ncol,                                          &
             pver            = pver,                                          &
             pverp           = pverp,                                         &
             ncnst           = pcnst,                                         &
             ztodt           = ztodt,                                         &
             do_diffusion_u_v= .true.,                                        & ! horizontal winds and
             do_diffusion_s  = .true.,                                        & ! dry static energy are diffused
             do_diffusion_const = do_diffusion_const_wet,                     & ! together with moist constituents.
             do_molecular_diffusion_const = do_molecular_diffusion_const,     &
             itaures         = itaures,                                       &
             t               = state%t(:ncol,:pver),                          &
             tint            = tint(:ncol,:pverp),                            &
             p               = p,                                             & ! Coords1D, pressure coordinates [Pa]
             rhoi            = rhoi(:ncol,:pverp),                            &
             taux            = taux(:ncol),                                   &
             tauy            = tauy(:ncol),                                   &
             shflx           = shflux(:ncol),                                 &
             cflx            = cflux(:ncol,:pcnst),                           &
             dse_top         = dse_top(:ncol),                                &
             kvh             = kvh(:ncol,:pverp),                             &
             kvm             = kvm(:ncol,:pverp),                             &
             kvq             = kvq(:ncol,:pverp),                             &
             cgs             = cgs(:ncol,:pverp),                             &
             cgh             = cgh(:ncol,:pverp),                             &
             ksrftms         = ksrftms(:ncol),                                &
             dragblj         = dragblj(:ncol,:pver),                          &
             qmincg          = qmincg(:pcnst),                                &
             ! input/output
             u               = u_tmp(:ncol,:pver),                            &
             v               = v_tmp(:ncol,:pver),                            &
             q               = q_tmp(:ncol,:pver,:pcnst),                     &
             dse             = s_tmp(:ncol,:pver),                            &
             tauresx         = tauresx(:ncol),                                &
             tauresy         = tauresy(:ncol),                                &
             ! below output
             dtk             = dtk(:ncol,:),                                  &
             tautmsx         = tautmsx(:ncol),                                &
             tautmsy         = tautmsy(:ncol),                                &
             topflx          = topflx(:ncol),                                 &
             errmsg          = errstring,                                     &
             ! arguments for Beljaars
             do_beljaars     = do_beljaars,                                   &
             ! arguments for molecular diffusion only.
             do_molec_diff   = do_molec_diff,                                 &
             use_temperature_molec_diff = waccmx_mode,                        &
             cpairv          = cpairv(:ncol,:,state%lchnk),                   &
             rairv           = rairv(:ncol,:,state%lchnk),                    &
             mbarv           = mbarv(:ncol,:,state%lchnk),                    &
             vd_lu_qdecomp   = vd_lu_qdecomp,                                 &
             ubc_mmr         = ubc_mmr(:ncol,:pcnst),                         &
             ubc_flux        = ubc_flux(:ncol,:pcnst),                        &
             kvt             = kvt(:ncol,:pverp),                             &
             pmid            = state%pmid(:ncol,:pver),                       &
             cnst_mw         = cnst_mw(:pcnst),                               &
             cnst_fixed_ubc  = cnst_fixed_ubc(:pcnst),                        &
             cnst_fixed_ubflx= cnst_fixed_ubflx(:pcnst),                      &
             nbot_molec      = nbot_molec,                                    &
             kq_scal         = kq_scal(:ncol,:pverp),                         &
             mw_fac          = mw_fac(:ncol,:pverp,:pcnst))

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_wet call from vertical_diffusion.")
     end if

     if( any( do_diffusion_const_dry ) ) then
        ! kvm is unused in the output here (since it was assigned
        ! above), so we use a temp kvm for the inout argument, and
        ! ignore the value output by compute_molec_diff.
        kvm_temp = kvm
        call compute_molec_diff(state%lchnk, pcols, pver, pcnst, ncol, &
             kvm_temp, kvt, tint, rhoi_dry, kq_scal, cnst_mw, &
             mw_fac, nbot_molec)

        ! Set up dry pressure coordinates for solver call.
        p_dry = Coords1D(state%pintdry(:ncol,:))

        call compute_vdiff( &
             ncol            = ncol,                                          &
             pver            = pver,                                          &
             pverp           = pverp,                                         &
             ncnst           = pcnst,                                         &
             ztodt           = ztodt,                                         &
             do_diffusion_u_v= .false.,                                       &
             do_diffusion_s  = .false.,                                       &
             do_diffusion_const = do_diffusion_const_dry,                     &
             do_molecular_diffusion_const = do_molecular_diffusion_const,     &
             itaures         = itaures,                                       &
             t               = state%t(:ncol,:pver),                          &
             tint            = tint(:ncol,:pverp),                            &
             p               = p_dry,                                         & ! Coords1D, pressure coordinates [Pa]
             rhoi            = rhoi_dry(:ncol,:pverp),                        &
             taux            = taux(:ncol),                                   &
             tauy            = tauy(:ncol),                                   &
             shflx           = shflux(:ncol),                                 &
             cflx            = cflux(:ncol,:pcnst),                           &
             dse_top         = dse_top(:ncol),                                &
             kvh             = kvh(:ncol,:pverp),                             &
             kvm             = kvm(:ncol,:pverp),                             &
             kvq             = kvq(:ncol,:pverp),                             &
             cgs             = cgs(:ncol,:pverp),                             &
             cgh             = cgh(:ncol,:pverp),                             &
             ksrftms         = ksrftms(:ncol),                                &
             dragblj         = dragblj(:ncol,:pver),                          &
             qmincg          = qmincg(:pcnst),                                &
             ! input/output
             u               = u_tmp(:ncol,:pver),                            &
             v               = v_tmp(:ncol,:pver),                            &
             q               = q_tmp(:ncol,:pver,:pcnst),                     &
             dse             = s_tmp(:ncol,:pver),                            &
             tauresx         = tauresx(:ncol),                                &
             tauresy         = tauresy(:ncol),                                &
             ! below output
             dtk             = dtk_temp(:ncol,:),                             & ! unused dummy
             tautmsx         = tautmsx_temp(:ncol),                           & ! unused dummy
             tautmsy         = tautmsy_temp(:ncol),                           & ! unused dummy
             topflx          = topflx_temp(:ncol),                            & ! unused dummy
             errmsg          = errstring,                                     &
             ! arguments for Beljaars
             do_beljaars     = do_beljaars,                                   &
             ! arguments for molecular diffusion only.
             do_molec_diff   = do_molec_diff,                                 &
             use_temperature_molec_diff = waccmx_mode,                        &
             cpairv          = cpairv(:ncol,:,state%lchnk),                   &
             rairv           = rairv(:ncol,:,state%lchnk),                    &
             mbarv           = mbarv(:ncol,:,state%lchnk),                    &
             vd_lu_qdecomp   = vd_lu_qdecomp,                                 &
             ubc_mmr         = ubc_mmr(:ncol,:pcnst),                         &
             ubc_flux        = ubc_flux(:ncol,:pcnst),                        &
             kvt             = kvt(:ncol,:pverp),                             &
             pmid            = state%pmiddry(:ncol,:pver),                    &
             cnst_mw         = cnst_mw(:pcnst),                               &
             cnst_fixed_ubc  = cnst_fixed_ubc(:pcnst),                        &
             cnst_fixed_ubflx= cnst_fixed_ubflx(:pcnst),                      &
             nbot_molec      = nbot_molec,                                    &
             kq_scal         = kq_scal(:ncol,:pverp),                         &
             mw_fac          = mw_fac(:ncol,:pverp,:pcnst))

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_dry call from vertical_diffusion.")

     end if
   end if

  ! For species not diffused, i.e., treated in dropmixnuc (vertically mixed by ndrop activation process)
  ! Just add the explicit surface fluxes to the lowest layer.  **NOTE** This code assumes wet mmr.
  tmp1(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
  do l = 1, pcnst
     if (cnst_ndropmixed(l)) then
        q_tmp(:ncol,pver,l) = q_tmp(:ncol,pver,l) + tmp1(:ncol) * cflux(:ncol,l)
     end if
  end do

  ! --------------------------------------------------------------- !
  ! Convert the new profiles into vertical diffusion tendencies.    !
  ! Convert KE dissipative heat change into "temperature" tendency. !
  ! --------------------------------------------------------------- !
  ! All variables are modified by vertical diffusion

  lq(:) = .TRUE.
  call physics_ptend_init(ptend,state%psetcols, "vertical diffusion", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  call vertical_diffusion_tendencies_run( &
     ncol        = ncol, &
     pver        = pver, &
     pcnst       = pcnst, &
     const_props = ccpp_const_props, &
     dt          = ztodt, &
     pdel        = state%pdel(:ncol,:pver), &
     pdeldry     = state%pdeldry(:ncol,:pver), &
     u0          = state%u(:ncol,:pver), &
     v0          = state%v(:ncol,:pver), &
     s0          = state%s(:ncol,:pver), &
     q0          = state%q(:ncol,:pver,:pcnst), &
     u1          = u_tmp(:ncol,:pver), &
     v1          = v_tmp(:ncol,:pver), &
     s1          = s_tmp(:ncol,:pver), &
     q1          = q_tmp(:ncol,:pver,:pcnst), &
     ! below output
     tend_s      = ptend%s(:ncol,:pver), &
     tend_u      = ptend%u(:ncol,:pver), &
     tend_v      = ptend%v(:ncol,:pver), &
     tend_q      = ptend%q(:ncol,:pver,:pcnst), &
     scheme_name = scheme_name, &
     errmsg      = errmsg, &
     errflg      = errflg)

  if(errflg /= 0) then
     call endrun("vertical_diffusion_tendencies_run: " // errmsg)
  endif

  ! -------------------------------------------------------- !
  ! Diagnostics and output writing after applying PBL scheme !
  ! -------------------------------------------------------- !

  if (.not. is_clubb_scheme) then

   ! Unified computation of diagnostics shared with CCPP-ized scheme
   call vertical_diffusion_diagnostic_profiles( &
        ncol                 = ncol, &
        pver                 = pver, &
        pverp                = pverp, &
        dt                   = ztodt, &
        latvap               = latvap, &
        latice               = latice, &
        zvir                 = zvir, &
        cpair                = cpair, &
        gravit               = gravit, &
        rair                 = rair, &
        pint                 = state%pint(:ncol, :pverp), &
        pmid                 = state%pmid(:ncol, :pver), &
        zi                   = state%zi(:ncol, :pverp), &
        zm                   = state%zm(:ncol, :pver), &
        kvh                  = kvh(:ncol, :pverp), &
        kvm                  = kvm(:ncol, :pverp), &
        cgh                  = cgh(:ncol, :pverp), &
        cgs                  = cgs(:ncol, :pverp), &
        cam_in_shf           = cam_in%shf(:ncol), &
        cam_in_cflx_wv       = cam_in%cflx(:ncol, 1), &
        cam_in_cflx_cldliq   = cam_in%cflx(:ncol, ixcldliq), &
        cam_in_cflx_cldice   = cam_in%cflx(:ncol, ixcldice), &
        tautotx              = tautotx(:ncol), & ! these are the surface stresses to HB
        tautoty              = tautoty(:ncol), & ! these are the surface stresses to HB
        t0                   = state%t(:ncol, :pver), &
        q0_wv                = state%q(:ncol, :pver, ixq), &
        q0_cldliq            = state%q(:ncol, :pver, ixcldliq), &
        q0_cldice            = state%q(:ncol, :pver, ixcldice), &
        s0                   = state%s(:ncol, :pver), &
        u0                   = state%u(:ncol, :pver), &
        v0                   = state%v(:ncol, :pver), &
        s1                   = s_tmp(:ncol, :pver), &
        u1                   = u_tmp(:ncol, :pver), &
        v1                   = v_tmp(:ncol, :pver), &
        tend_s               = ptend%s(:ncol, :pver), &
        tend_u               = ptend%u(:ncol, :pver), &
        tend_v               = ptend%v(:ncol, :pver), &
        tend_q_wv            = ptend%q(:ncol, :pver, ixq), &
        tend_q_cldliq        = ptend%q(:ncol, :pver, ixcldliq), &
        tend_q_cldice        = ptend%q(:ncol, :pver, ixcldice), &
        qt_pre_PBL           = qt_prePBL(:ncol, :pver), &
        sl_pre_PBL           = sl_prePBL(:ncol, :pver), &
        slv_pre_PBL          = slv_prePBL(:ncol, :pver), &
        rh_pre_PBL           = ftem_prePBL(:ncol, :pver), &
        qt_aft_PBL           = qt(:ncol, :pver), &
        sl_aft_PBL           = sl(:ncol, :pver), &
        slv_aft_PBL          = slv(:ncol, :pver), &
        qv_aft_PBL           = qv_aft_PBL(:ncol, :pver), &
        ql_aft_PBL           = ql_aft_PBL(:ncol, :pver), &
        qi_aft_PBL           = qi_aft_PBL(:ncol, :pver), &
        t_aft_PBL            = t_aftPBL(:ncol, :pver), &
        rh_aft_PBL           = ftem_aftPBL(:ncol, :pver), &
        u_aft_PBL            = u_aft_PBL(:ncol, :pver), &
        v_aft_PBL            = v_aft_PBL(:ncol, :pver), &
        slflx                = slflx(:ncol, :pverp), &
        qtflx                = qtflx(:ncol, :pverp), &
        uflx                 = uflx(:ncol, :pverp), &
        vflx                 = vflx(:ncol, :pverp), &
        slflx_cg             = slflx_cg(:ncol, :pverp), &
        qtflx_cg             = qtflx_cg(:ncol, :pverp), &
        uflx_cg              = uflx_cg(:ncol, :pverp), &
        vflx_cg              = vflx_cg(:ncol, :pverp), &
        slten                = slten(:ncol, :pver), &
        qtten                = qtten(:ncol, :pver), &
        tten                 = tten(:ncol, :pver), &
        rhten                = rhten(:ncol, :pver))

     if (trim(shallow_scheme) == 'UNICON') then
        call pbuf_get_field(pbuf, qtl_flx_idx,  qtl_flx)
        call pbuf_get_field(pbuf, qti_flx_idx,  qti_flx)
        qtl_flx(:ncol,1) = 0._r8
        qti_flx(:ncol,1) = 0._r8
        do k = 2, pver
           do i = 1, ncol
              ! For use in the cloud macrophysics
              ! Note that density is not added here. Also, only consider local transport term.
              qtl_flx(i,k) = - kvh(i,k)*(q_tmp(i,k-1,1)-q_tmp(i,k,1)+q_tmp(i,k-1,ixcldliq)-q_tmp(i,k,ixcldliq))/&
                   (state%zm(i,k-1)-state%zm(i,k))
              qti_flx(i,k) = - kvh(i,k)*(q_tmp(i,k-1,1)-q_tmp(i,k,1)+q_tmp(i,k-1,ixcldice)-q_tmp(i,k,ixcldice))/&
                   (state%zm(i,k-1)-state%zm(i,k))
           end do
        end do
        do i = 1, ncol
           rhoair = state%pint(i,pverp)/(rair*((slv(i,pver)-gravit*state%zi(i,pverp))/cpair))
           qtl_flx(i,pverp) = cam_in%cflx(i,1)/rhoair
           qti_flx(i,pverp) = cam_in%cflx(i,1)/rhoair
        end do
     end if

  end if

  ! ------------------------------------------------------------ !
  ! In order to perform 'pseudo-conservative variable diffusion' !
  ! perform the following two stages:                            !
  !                                                              !
  ! I.  Re-set (1) 'qvten' by 'qtten', and 'qlten = qiten = 0'   !
  !            (2) 'sten'  by 'slten', and                       !
  !            (3) 'qlten = qiten = 0'                           !
  !                                                              !
  ! II. Apply 'positive_moisture'                                !
  !                                                              !
  ! ------------------------------------------------------------ !

  if( eddy_scheme .eq. 'diag_TKE' .and. do_pseudocon_diff ) then

     ptend%q(:ncol,:pver,1) = qtten(:ncol,:pver)
     ptend%s(:ncol,:pver)   = slten(:ncol,:pver)
     ptend%q(:ncol,:pver,ixcldliq) = 0._r8
     ptend%q(:ncol,:pver,ixcldice) = 0._r8
     if (ixnumliq > 0) ptend%q(:ncol,:pver,ixnumliq) = 0._r8
     if (ixnumice > 0) ptend%q(:ncol,:pver,ixnumice) = 0._r8

     do i = 1, ncol
        do k = 1, pver
           qv_pro(i,k) = state%q(i,k,1)        + ptend%q(i,k,1)             * ztodt
           ql_pro(i,k) = state%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)      * ztodt
           qi_pro(i,k) = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)      * ztodt
           s_pro(i,k)  = state%s(i,k)          + ptend%s(i,k)               * ztodt
           t_pro(i,k)  = state%t(i,k)          + (1._r8/cpair)*ptend%s(i,k) * ztodt
        end do
     end do

     call positive_moisture( &
        ncol   = ncol, &
        mkx    = pver, &
        dt     = ztodt, &
        cp     = cpair, &
        xlv    = latvap, &
        xls    = latvap+latice, &
        qvmin  = qmin(1), &
        qlmin  = qmin(ixcldliq), &
        qimin  = qmin(ixcldice), &
        dp     = state%pdel(:ncol,pver:1:-1), &
        qv     = qv_pro(:ncol,pver:1:-1), &
        ql     = ql_pro(:ncol,pver:1:-1), &
        qi     = qi_pro(:ncol,pver:1:-1), &
        t      = t_pro(:ncol,pver:1:-1), &
        s      = s_pro(:ncol,pver:1:-1), &
        qvten  = ptend%q(:ncol,pver:1:-1,1), &
        qlten  = ptend%q(:ncol,pver:1:-1,ixcldliq), &
        qiten  = ptend%q(:ncol,pver:1:-1,ixcldice), &
        sten   = ptend%s(:ncol,pver:1:-1))

  end if

  ! -------------------------------------------------------------- !
  ! mass conservation check.........
  ! -------------------------------------------------------------- !
  if (diff_cnsrv_mass_check) then

     ! Conservation check
     do m = 1, pcnst
        fixed_ubc: if ((.not.cnst_fixed_ubc(m)).and.(.not.cnst_fixed_ubflx(m))) then
           col_loop: do i = 1, ncol
              sum1 = 0._r8
              sum2 = 0._r8
              sum3 = 0._r8
              do k = 1, pver
                 if(cnst_get_type_byind(m).eq.'wet') then
                    pdelx = state%pdel(i,k)
                 else
                    pdelx = state%pdeldry(i,k)
                 endif
                 sum1 = sum1 + state%q(i,k,m)*pdelx/gravit                          ! total column
                 sum2 = sum2 +(state%q(i,k,m)+ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! total column after tendancy is applied
                 sum3 = sum3 +(               ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! rate of change in column
              enddo
              sum1 = sum1 + (cam_in%cflx(i,m) * ztodt) ! add in surface flux (kg/m2)
              sflx = (cam_in%cflx(i,m) * ztodt)
              if (sum1>1.e-36_r8) then
                 if( abs((sum2-sum1)/sum1) .gt. 1.e-12_r8  ) then
                    nstep = get_nstep()
                    write(iulog,'(a,a8,a,I4,2f8.3,5e25.16)') &
                         'MASSCHECK vert diff : nstep,lon,lat,mass1,mass2,sum3,sflx,rel-diff : ', &
                         trim(cnst_name(m)), ' : ', nstep, state%lon(i)*180._r8/pi, state%lat(i)*180._r8/pi, &
                         sum1, sum2, sum3, sflx, abs(sum2-sum1)/sum1
!xxx                    call endrun('vertical_diffusion_tend : mass not conserved' )
                 endif
              endif
           enddo col_loop
        endif fixed_ubc
     enddo
  endif

  ! ------------------------------------------- !
  ! Writing standard output variables           !
  ! ------------------------------------------- !

  if (.not. is_clubb_scheme) then

     ! Write profile output of before diffusion scheme
     call outfld( 'qt_pre_PBL   ', qt_prePBL,                 pcols, lchnk )
     call outfld( 'sl_pre_PBL   ', sl_prePBL,                 pcols, lchnk )
     call outfld( 'slv_pre_PBL  ', slv_prePBL,                pcols, lchnk )
     call outfld( 'u_pre_PBL    ', state%u,                   pcols, lchnk )
     call outfld( 'v_pre_PBL    ', state%v,                   pcols, lchnk )
     call outfld( 'qv_pre_PBL   ', state%q(:,:,1),            pcols, lchnk )
     call outfld( 'ql_pre_PBL   ', state%q(:,:,ixcldliq),     pcols, lchnk )
     call outfld( 'qi_pre_PBL   ', state%q(:,:,ixcldice),     pcols, lchnk )
     call outfld( 't_pre_PBL    ', state%t,                   pcols, lchnk )
     call outfld( 'rh_pre_PBL   ', ftem_prePBL,               pcols, lchnk )

     ! Standard output variables
     call outfld( 'QT'           , qt,                        pcols, lchnk )
     call outfld( 'SL'           , sl,                        pcols, lchnk )
     call outfld( 'SLV'          , slv,                       pcols, lchnk )
     call outfld( 'SLFLX'        , slflx,                     pcols, lchnk )
     call outfld( 'QTFLX'        , qtflx,                     pcols, lchnk )
     call outfld( 'UFLX'         , uflx,                      pcols, lchnk )
     call outfld( 'VFLX'         , vflx,                      pcols, lchnk )
     call outfld( 'TKE'          , tke,                       pcols, lchnk )

     call outfld( 'PBLH'         , pblh,                      pcols, lchnk )
     call outfld( 'TPERT'        , tpert,                     pcols, lchnk )
     call outfld( 'QPERT'        , qpert,                     pcols, lchnk )

     ! State variables after PBL scheme for detailed analysis
     call outfld( 'sl_aft_PBL'   , sl,                        pcols, lchnk )
     call outfld( 'qt_aft_PBL'   , qt,                        pcols, lchnk )
     call outfld( 'slv_aft_PBL'  , slv,                       pcols, lchnk )
     call outfld( 'u_aft_PBL'    , u_aft_PBL,                 pcols, lchnk )
     call outfld( 'v_aft_PBL'    , v_aft_PBL,                 pcols, lchnk )
     call outfld( 'qv_aft_PBL'   , qv_aft_PBL,                pcols, lchnk )
     call outfld( 'ql_aft_PBL'   , ql_aft_PBL,                pcols, lchnk )
     call outfld( 'qi_aft_PBL'   , qi_aft_PBL,                pcols, lchnk )
     call outfld( 't_aft_PBL '   , t_aftPBL,                  pcols, lchnk )
     call outfld( 'rh_aft_PBL'   , ftem_aftPBL,               pcols, lchnk )
     call outfld( 'slflx_PBL'    , slflx,                     pcols, lchnk )
     call outfld( 'qtflx_PBL'    , qtflx,                     pcols, lchnk )
     call outfld( 'uflx_PBL'     , uflx,                      pcols, lchnk )
     call outfld( 'vflx_PBL'     , vflx,                      pcols, lchnk )
     call outfld( 'slflx_cg_PBL' , slflx_cg,                  pcols, lchnk )
     call outfld( 'qtflx_cg_PBL' , qtflx_cg,                  pcols, lchnk )
     call outfld( 'uflx_cg_PBL'  , uflx_cg,                   pcols, lchnk )
     call outfld( 'vflx_cg_PBL'  , vflx_cg,                   pcols, lchnk )
     call outfld( 'slten_PBL'    , slten,                     pcols, lchnk )
     call outfld( 'qtten_PBL'    , qtten,                     pcols, lchnk )
     call outfld( 'uten_PBL'     , ptend%u,                   pcols, lchnk )
     call outfld( 'vten_PBL'     , ptend%v,                   pcols, lchnk )
     call outfld( 'qvten_PBL'    , ptend%q(:,:,1),            pcols, lchnk )
     call outfld( 'qlten_PBL'    , ptend%q(:,:,ixcldliq),     pcols, lchnk )
     call outfld( 'qiten_PBL'    , ptend%q(:,:,ixcldice),     pcols, lchnk )
     call outfld( 'tten_PBL'     , tten,                      pcols, lchnk )
     call outfld( 'rhten_PBL'    , rhten,                     pcols, lchnk )
  end if

  call outfld( 'USTAR'        , ustar,                     pcols, lchnk )
  call outfld( 'KVH'          , kvh,                       pcols, lchnk )
  call outfld( 'KVT'          , kvt,                       pcols, lchnk )
  call outfld( 'KVM'          , kvm,                       pcols, lchnk )
  call outfld( 'CGS'          , cgs,                       pcols, lchnk )
  dtk(:ncol,:) = dtk(:ncol,:) / cpair / ztodt      ! Normalize heating for history
  call outfld( 'DTVKE'        , dtk,                       pcols, lchnk )
  dtk(:ncol,:) = ptend%s(:ncol,:) / cpair          ! Normalize heating for history using dtk
  call outfld( 'DTV'          , dtk,                       pcols, lchnk )
  call outfld( 'DUV'          , ptend%u,                   pcols, lchnk )
  call outfld( 'DVV'          , ptend%v,                   pcols, lchnk )
  do m = 1, pcnst
     call outfld( vdiffnam(m) , ptend%q(1,1,m),            pcols, lchnk )
  end do
  if( do_molec_diff ) then
     call outfld( 'TTPXMLC'  , topflx,                    pcols, lchnk )
  end if

  call p%finalize()
  call p_dry%finalize()

end subroutine vertical_diffusion_tend

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine positive_moisture( cp, xlv, xls, ncol, mkx, dt, qvmin, qlmin, qimin, &
     dp, qv, ql, qi, t, s, qvten, qlten, qiten, sten )
  ! ------------------------------------------------------------------------------- !
  ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
  ! force them to be larger than minimum value by (1) condensating water vapor      !
  ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
  ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
  ! Update final state variables and tendencies associated with this correction.    !
  ! If any condensation happens, update (s,t) too.                                  !
  ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
  ! input tendencies.                                                               !
  ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          !
  ! ------------------------------------------------------------------------------- !
  integer,  intent(in)     :: ncol            ! Number of atmospheric columns [count]
  integer,  intent(in)     :: mkx             ! Number of vertical levels [count]
  real(r8), intent(in)     :: cp              ! Specific heat of dry air at constant pressure [J kg-1 K-1]
  real(r8), intent(in)     :: xlv             ! Latent heat of vaporization [J kg-1]
  real(r8), intent(in)     :: xls             ! Latent heat of sublimation [J kg-1]
  real(r8), intent(in)     :: dt              ! Time step [s]
  real(r8), intent(in)     :: qvmin           ! Minimum water vapor mixing ratio [kg kg-1]
  real(r8), intent(in)     :: qlmin           ! Minimum liquid water mixing ratio [kg kg-1]
  real(r8), intent(in)     :: qimin           ! Minimum ice water mixing ratio [kg kg-1]
  real(r8), intent(in)     :: dp(ncol,mkx)    ! Pressure thickness of layers [Pa]
  real(r8), intent(inout)  :: qv(ncol,mkx)    ! Water vapor mixing ratio [kg kg-1]
  real(r8), intent(inout)  :: ql(ncol,mkx)    ! Cloud liquid water mixing ratio [kg kg-1]
  real(r8), intent(inout)  :: qi(ncol,mkx)    ! Cloud ice water mixing ratio [kg kg-1]
  real(r8), intent(inout)  :: t(ncol,mkx)     ! Temperature [K]
  real(r8), intent(inout)  :: s(ncol,mkx)     ! Dry static energy [J kg-1]
  real(r8), intent(inout)  :: qvten(ncol,mkx) ! Water vapor tendency [kg kg-1 s-1]
  real(r8), intent(inout)  :: qlten(ncol,mkx) ! Liquid water tendency [kg kg-1 s-1]
  real(r8), intent(inout)  :: qiten(ncol,mkx) ! Ice water tendency [kg kg-1 s-1]
  real(r8), intent(inout)  :: sten(ncol,mkx)  ! Dry static energy tendency [J kg-1 s-1]

  integer :: i, k
  real(r8) :: dql, dqi, dqv, sum, aa, dum

  ! Modification : I should check whether this is exactly same as the one used in
  !                shallow convection and cloud macrophysics.

  do i = 1, ncol
     do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
        dql        = max(0._r8,1._r8*qlmin-ql(i,k))
        dqi        = max(0._r8,1._r8*qimin-qi(i,k))
        qlten(i,k) = qlten(i,k) +  dql/dt
        qiten(i,k) = qiten(i,k) +  dqi/dt
        qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
        sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
        ql(i,k)    = ql(i,k) +  dql
        qi(i,k)    = qi(i,k) +  dqi
        qv(i,k)    = qv(i,k) -  dql - dqi
        s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
        t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
        dqv        = max(0._r8,1._r8*qvmin-qv(i,k))
        qvten(i,k) = qvten(i,k) + dqv/dt
        qv(i,k)    = qv(i,k)    + dqv
        if( k .ne. 1 ) then
           qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
           qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
        endif
        qv(i,k) = max(qv(i,k),qvmin)
        ql(i,k) = max(ql(i,k),qlmin)
        qi(i,k) = max(qi(i,k),qimin)
     end do
     ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally
     ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
     ! preserves column moisture.
     if( dqv .gt. 1.e-20_r8 ) then
        sum = 0._r8
        do k = 1, mkx
           if( qv(i,k) .gt. 2._r8*qvmin ) sum = sum + qv(i,k)*dp(i,k)
        enddo
        aa = dqv*dp(i,1)/max(1.e-20_r8,sum)
        if( aa .lt. 0.5_r8 ) then
           do k = 1, mkx
              if( qv(i,k) .gt. 2._r8*qvmin ) then
                 dum        = aa*qv(i,k)
                 qv(i,k)    = qv(i,k) - dum
                 qvten(i,k) = qvten(i,k) - dum/dt
              endif
           enddo
        else
           write(iulog,*) 'Full positive_moisture is impossible in vertical_diffusion'
        endif
     endif
  end do

end subroutine positive_moisture

end module vertical_diffusion
