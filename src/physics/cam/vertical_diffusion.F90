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
!      - The conditionals controlled by "do_pbl_diags" are somewhat scattered. It might be better to   !
!        pull out these diagnostic calculations and outfld calls into separate functions.              !
!                                                                                                      !
!---------------------------Code history-------------------------------------------------------------- !
! J. Rosinski : Jun. 1992                                                                              !
! J. McCaa    : Sep. 2004                                                                              !
! S. Park     : Aug. 2006, Dec. 2008. Jan. 2010                                                        !
!----------------------------------------------------------------------------------------------------- !

use shr_kind_mod,     only : r8 => shr_kind_r8, i4=> shr_kind_i4
use ppgrid,           only : pcols, pver, pverp
use constituents,     only : pcnst
use diffusion_solver, only : vdiff_selector
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
logical, parameter   :: wstarent = .true.            ! Use wstar (.true.) or TKE (.false.) entrainment closure
! ( when 'diag_TKE' scheme is selected )
logical              :: do_pseudocon_diff = .false.  ! If .true., do pseudo-conservative variables diffusion

character(len=16)    :: shallow_scheme               ! Shallow convection scheme

type(vdiff_selector) :: fieldlist_wet                ! Logical switches for moist mixing ratio diffusion
type(vdiff_selector) :: fieldlist_dry                ! Logical switches for dry mixing ratio diffusion
type(vdiff_selector) :: fieldlist_molec              ! Logical switches for molecular diffusion
integer              :: tke_idx, kvh_idx, kvm_idx    ! TKE and eddy diffusivity indices for fields in the physics buffer
integer              :: kvt_idx                      ! Index for kinematic molecular conductivity
integer              :: turbtype_idx, smaw_idx       ! Turbulence type and instability functions
integer              :: tauresx_idx, tauresy_idx     ! Redisual stress for implicit surface stress

character(len=fieldname_len) :: vdiffnam(pcnst)      ! Names of vertical diffusion tendencies
integer              :: ixcldice, ixcldliq           ! Constituent indices for cloud liquid and ice water
integer              :: ixnumice, ixnumliq

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

logical              :: diff_cnsrv_mass_check        ! do mass conservation check
logical              :: do_iss                       ! switch for implicit turbulent surface stress
logical              :: prog_modal_aero = .false.    ! set true if prognostic modal aerosols are present
integer              :: pmam_ncnst = 0               ! number of prognostic modal aerosol constituents
integer, allocatable :: pmam_cnst_idx(:)             ! constituent indices of prognostic modal aerosols

logical              :: do_pbl_diags = .false.
logical              :: waccmx_mode = .false.

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

  if (eddy_scheme == 'diag_TKE' .or. eddy_scheme == 'SPCAM_m2005' ) call eddy_diff_readnl(nlfile)

end subroutine vd_readnl

! =============================================================================== !
!                                                                                 !
! =============================================================================== !

subroutine vd_register()

  !------------------------------------------------ !
  ! Register physics buffer fields and constituents !
  !------------------------------------------------ !

  use physics_buffer,      only : pbuf_add_field, dtype_r8, dtype_i4, pbuf_get_index
  use trb_mtn_stress_cam,  only : trb_mtn_stress_register
  use beljaars_drag_cam,   only : beljaars_drag_register
  use eddy_diff_cam,       only : eddy_diff_register

  integer :: err_idx, idx

  ! Add fields to physics buffer

  ! kvt is used by gw_drag.  only needs physpkg scope.
  call pbuf_add_field('kvt', 'physpkg', dtype_r8, (/pcols,pverp/), kvt_idx)

  idx = pbuf_get_index('kvh',errcode=err_idx)
  if (err_idx == -1) then
    call pbuf_add_field('kvh',      'global', dtype_r8, (/pcols, pverp/), kvh_idx)
  end if

  call pbuf_add_field('kvm',      'global', dtype_r8, (/pcols, pverp/), kvm_idx )
  call pbuf_add_field('pblh',     'global', dtype_r8, (/pcols/),        pblh_idx)
  call pbuf_add_field('tke',      'global', dtype_r8, (/pcols, pverp/), tke_idx)
  call pbuf_add_field('turbtype', 'global', dtype_i4, (/pcols, pverp/), turbtype_idx)
  call pbuf_add_field('smaw',     'global', dtype_r8, (/pcols, pverp/), smaw_idx)

  call pbuf_add_field('tauresx',  'global', dtype_r8, (/pcols/),        tauresx_idx)
  call pbuf_add_field('tauresy',  'global', dtype_r8, (/pcols/),        tauresy_idx)

  call pbuf_add_field('tpert', 'global', dtype_r8, (/pcols/),                       tpert_idx)
  call pbuf_add_field('qpert', 'global', dtype_r8, (/pcols,pcnst/),                 qpert_idx)

  if (trim(shallow_scheme) == 'UNICON') then
     call pbuf_add_field('qtl_flx',  'global', dtype_r8, (/pcols, pverp/), qtl_flx_idx)
     call pbuf_add_field('qti_flx',  'global', dtype_r8, (/pcols, pverp/), qti_flx_idx)
  end if

  ! diag_TKE fields
  if (eddy_scheme == 'diag_TKE' .or. eddy_scheme == 'SPCAM_m2005') then
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
  use hb_diff,           only : init_hb_diff
  use molec_diff,        only : init_molec_diff
  use diffusion_solver,  only : init_vdiff, new_fieldlist_vdiff, vdiff_select
  use constituents,      only : cnst_get_ind, cnst_get_type_byind, cnst_name, cnst_get_molec_byind
  use spmd_utils,        only : masterproc
  use ref_pres,          only : press_lim_idx, pref_mid
  use physics_buffer,    only : pbuf_set_field, pbuf_get_index, physics_buffer_desc
  use rad_constituents,  only : rad_cnst_get_info, rad_cnst_get_mode_num_idx, &
       rad_cnst_get_mam_mmr_idx
  use trb_mtn_stress_cam,only : trb_mtn_stress_init
  use beljaars_drag_cam, only : beljaars_drag_init
  use upper_bc,          only : ubc_init
  use phys_control,      only : waccmx_is, fv_am_correction, cam_physpkg_is

  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  character(128) :: errstring   ! Error status for init_vdiff
  integer        :: ntop_eddy   ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
  integer        :: nbot_eddy   ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
  integer        :: k           ! Vertical loop index

  real(r8), parameter :: ntop_eddy_pres = 1.e-7_r8 ! Pressure below which eddy diffusion is not done in WACCM-X. (Pa)

  integer :: im, l, m, nmodes, nspec

  logical :: history_amwg                 ! output the variables used by the AMWG diag package
  logical :: history_eddy                 ! output the eddy variables
  logical :: history_budget               ! Output tendencies and state variables for CAM4 T, qv, ql, qi
  integer :: history_budget_histfile_num  ! output history file number for budget fields
  logical :: history_waccm                ! output variables of interest for WACCM runs

  ! ----------------------------------------------------------------- !

  if (masterproc) then
     write(iulog,*)'Initializing vertical diffusion (vertical_diffusion_init)'
  end if

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

  ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
  call phys_getopts(prog_modal_aero_out=prog_modal_aero)
  if (prog_modal_aero) then

     ! Get the constituent indices of the number and mass mixing ratios of the modal
     ! aerosols.
     !
     ! N.B. - This implementation assumes that the prognostic modal aerosols are
     !        impacting the climate calculation (i.e., can get info from list 0).
     !

     ! First need total number of mam constituents
     call rad_cnst_get_info(0, nmodes=nmodes)
     do m = 1, nmodes
        call rad_cnst_get_info(0, m, nspec=nspec)
        pmam_ncnst = pmam_ncnst + 1 + nspec
     end do

     allocate(pmam_cnst_idx(pmam_ncnst))

     ! Get the constituent indicies
     im = 1
     do m = 1, nmodes
        call rad_cnst_get_mode_num_idx(m, pmam_cnst_idx(im))
        im = im + 1
        call rad_cnst_get_info(0, m, nspec=nspec)
        do l = 1, nspec
           call rad_cnst_get_mam_mmr_idx(m, l, pmam_cnst_idx(im))
           im = im + 1
        end do
     end do
  end if

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
  if ( waccmx_mode ) then
     ntop_eddy  = press_lim_idx(ntop_eddy_pres, top=.true.)
  else
     ntop_eddy = 1
  end if
  nbot_eddy  = pver

  if (masterproc) write(iulog, fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy

  select case ( eddy_scheme )
  case ( 'diag_TKE', 'SPCAM_m2005' )
     if( masterproc ) write(iulog,*) &
          'vertical_diffusion_init: eddy_diffusivity scheme: UW Moist Turbulence Scheme by Bretherton and Park'
     call eddy_diff_init(pbuf2d, ntop_eddy, nbot_eddy)
  case ( 'HB', 'HBR', 'SPCAM_sam1mom')
     if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'
     call init_hb_diff(gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, &
          karman, eddy_scheme)
     call addfld('HB_ri',      (/ 'lev' /),  'A', 'no',  'Richardson Number (HB Scheme), I' )
  case ( 'CLUBB_SGS' )
     do_pbl_diags = .true.
       call init_hb_diff(gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, &
            karman, eddy_scheme)
     !
     ! run HB scheme where CLUBB is not active when running cam_dev
     ! else init_hb_diff is called just for diagnostic purposes
     !
     if (cam_physpkg_is("cam_dev")) then
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

  ! Use fieldlist_wet to select the fields which will be diffused using moist mixing ratios ( all by default )
  ! Use fieldlist_dry to select the fields which will be diffused using dry   mixing ratios.

  fieldlist_wet = new_fieldlist_vdiff( pcnst)
  fieldlist_dry = new_fieldlist_vdiff( pcnst)
  fieldlist_molec = new_fieldlist_vdiff( pcnst)

  if( vdiff_select( fieldlist_wet, 'u' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'u' ) )
  if( vdiff_select( fieldlist_wet, 'v' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'v' ) )
  if( vdiff_select( fieldlist_wet, 's' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 's' ) )

  constit_loop: do k = 1, pcnst

     if (prog_modal_aero) then
        ! Do not diffuse droplet number - treated in dropmixnuc
        if (k == ixnumliq) cycle constit_loop
        ! Don't diffuse modal aerosol - treated in dropmixnuc
        do m = 1, pmam_ncnst
           if (k == pmam_cnst_idx(m)) cycle constit_loop
        enddo
     end if

     ! Convert all constituents to wet before doing diffusion.
     if( vdiff_select( fieldlist_wet, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'q', k ) )

     ! ----------------------------------------------- !
     ! Select constituents for molecular diffusion     !
     ! ----------------------------------------------- !
     if ( cnst_get_molec_byind(k) .eq. 'minor' ) then
        if( vdiff_select(fieldlist_molec,'q',k) .ne. '' ) call endrun( vdiff_select( fieldlist_molec,'q',k ) )
     endif

  end do constit_loop

  ! ------------------------ !
  ! Diagnostic output fields !
  ! ------------------------ !

  do k = 1, pcnst
     vdiffnam(k) = 'VD'//cnst_name(k)
     if( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
     call addfld( vdiffnam(k), (/ 'lev' /), 'A', 'kg/kg/s', 'Vertical diffusion of '//cnst_name(k) )
  end do

  if (.not. do_pbl_diags) then
     call addfld( 'PBLH'        , horiz_only    , 'A', 'm'      , 'PBL height'                                         )
     call addfld( 'QT'          , (/ 'lev' /)   , 'A', 'kg/kg'  , 'Total water mixing ratio'                           )
     call addfld( 'SL'          , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liquid water static energy'                         )
     call addfld( 'SLV'         , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liq wat virtual static energy'                      )
     call addfld( 'SLFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Liquid static energy flux'                          )
     call addfld( 'QTFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Total water flux'                                   )
     call addfld( 'TKE'         , (/ 'ilev' /)  , 'A', 'm2/s2'  , 'Turbulent Kinetic Energy'                           )
     call addfld( 'TPERT'       , horiz_only    , 'A', 'K'      , 'Perturbation temperature (eddies in PBL)'           )
     call addfld( 'QPERT'       , horiz_only    , 'A', 'kg/kg'  , 'Perturbation specific humidity (eddies in PBL)'     )

     call addfld( 'UFLX'        , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Zonal momentum flux'                                )
     call addfld( 'VFLX'        , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Meridional momentm flux'                            )
     call register_vector_field('UFLX', 'VFLX')
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

  ! ---------------------------------------------------------------------------- !
  ! Below ( with '_PBL') are for detailed analysis of UW Moist Turbulence Scheme !
  ! ---------------------------------------------------------------------------- !

  if (.not. do_pbl_diags) then

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
     if (.not. do_pbl_diags) then
        call add_default( 'PBLH'      , 1, ' ' )
     end if
  endif

  if (history_eddy) then
     call add_default( 'UFLX    ', 1, ' ' )
     call add_default( 'VFLX    ', 1, ' ' )
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

  ! Initialization of some pbuf fields
  if (is_first_step()) then
     ! Initialization of pbuf fields tke, kvh, kvm are done in phys_inidat
     call pbuf_set_field(pbuf2d, turbtype_idx, 0    )
     call pbuf_set_field(pbuf2d, smaw_idx,     0.0_r8)
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
  use physics_buffer,     only : pbuf_get_index
  use physics_types,      only : physics_state, physics_ptend, physics_ptend_init
  use physics_types,      only : set_dry_to_wet, set_wet_to_dry

  use camsrfexch,         only : cam_in_t
  use cam_history,        only : outfld

  use trb_mtn_stress_cam, only : trb_mtn_stress_tend
  use beljaars_drag_cam,  only : beljaars_drag_tend
  use eddy_diff_cam,      only : eddy_diff_tend
  use hb_diff,            only : compute_hb_diff
  use wv_saturation,      only : qsat
  use molec_diff,         only : compute_molec_diff, vd_lu_qdecomp
  use constituents,       only : qmincg, qmin, cnst_type
  use diffusion_solver,   only : compute_vdiff, any, operator(.not.)
  use air_composition,    only : cpairv, rairv !Needed for calculation of upward H flux
  use time_manager,       only : get_nstep
  use constituents,       only : cnst_get_type_byind, cnst_name, &
       cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx
  use physconst,          only : pi
  use pbl_utils,          only : virtem, calc_obklen, calc_ustar
  use upper_bc,           only : ubc_get_vals, ubc_fixed_temp
  use upper_bc,           only : ubc_get_flxs
  use coords_1d,          only : Coords1D
  use phys_control,       only : cam_physpkg_is

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
  integer(i4),pointer :: turbtype(:,:)                            ! Turbulent interface types [ no unit ]
  real(r8), pointer   :: smaw(:,:)                                ! Normalized Galperin instability function
  ! ( 0<= <=4.964 and 1 at neutral )

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
  real(r8) :: khfs(pcols)                                         ! sfc kinematic heat flux [mK/s]
  real(r8) :: kqfs(pcols)                                         ! sfc kinematic water vapor flux [m/s]
  real(r8) :: kbfs(pcols)                                         ! sfc kinematic buoyancy flux [m^2/s^3]

  real(r8) :: ftem(pcols,pver)                                    ! Saturation vapor pressure before PBL
  real(r8) :: ftem_prePBL(pcols,pver)                             ! Saturation vapor pressure before PBL
  real(r8) :: ftem_aftPBL(pcols,pver)                             ! Saturation vapor pressure after PBL
  real(r8) :: tem2(pcols,pver)                                    ! Saturation specific humidity and RH
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
  integer  :: clubbtop_idx
  integer,  pointer :: clubbtop(:)   ! (pcols)

  logical  :: lq(pcnst)

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  ! Assume 'wet' mixing ratios in diffusion code.
  call set_dry_to_wet(state)

  rztodt = 1._r8 / ztodt
  lchnk  = state%lchnk
  ncol   = state%ncol

  call pbuf_get_field(pbuf, tauresx_idx,  tauresx)
  call pbuf_get_field(pbuf, tauresy_idx,  tauresy)
  call pbuf_get_field(pbuf, tpert_idx,    tpert)
  call pbuf_get_field(pbuf, qpert_idx,    qpert)
  call pbuf_get_field(pbuf, pblh_idx,     pblh)
  call pbuf_get_field(pbuf, turbtype_idx, turbtype)

  ! Interpolate temperature to interfaces.
  do k = 2, pver
     do i = 1, ncol
        tint(i,k)  = 0.5_r8 * ( state%t(i,k) + state%t(i,k-1) )
     end do
  end do
  tint(:ncol,pver+1) = state%t(:ncol,pver)

  ! Get upper boundary values
  call ubc_get_vals( state%lchnk, ncol, state%pint, state%zi, ubc_t, ubc_mmr )

  if (waccmx_mode) then
     call ubc_get_flxs( state%lchnk, ncol, state%pint, state%zi, state%t, state%q, state%phis, ubc_flux )
     ! For WACCM-X, set ubc temperature to extrapolate from next two lower interface level temperatures
     tint(:ncol,1) = 1.5_r8*tint(:ncol,2)-.5_r8*tint(:ncol,3)
  else if(ubc_fixed_temp) then
     tint(:ncol,1) = ubc_t(:ncol)
  else
     tint(:ncol,1) = state%t(:ncol,1)
  end if

  ! Set up pressure coordinates for solver calls.
  p = Coords1D(state%pint(:ncol,:))
  p_dry = Coords1D(state%pintdry(:ncol,:))

  !------------------------------------------------------------------------
  !  Check to see if constituent dependent gas constant needed (WACCM-X)
  !------------------------------------------------------------------------
  if (waccmx_mode) then
     rairi(:ncol,1) = rairv(:ncol,1,lchnk)
     do k = 2, pver
        do i = 1, ncol
           rairi(i,k) = 0.5_r8 * (rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
        end do
     end do
     rairi(:ncol,pver+1) = rairv(:ncol,pver,lchnk)
  else
     rairi(:ncol,:pver+1) = rair
  endif

  ! Compute rho at interfaces.
  do k = 1, pver+1
     do i = 1, ncol
        rhoi(i,k)  = p%ifc(i,k) / (rairi(i,k)*tint(i,k))
     end do
  end do

  ! Compute rho_dry at interfaces.
  do k = 1, pver+1
     do i = 1, ncol
        rhoi_dry(i,k)  = p_dry%ifc(i,k) / (rairi(i,k)*tint(i,k))
     end do
  end do

  ! ---------------------------------------- !
  ! Computation of turbulent mountain stress !
  ! ---------------------------------------- !

  ! Consistent with the computation of 'normal' drag coefficient, we are using
  ! the raw input (u,v) to compute 'ksrftms', not the provisionally-marched 'u,v'
  ! within the iteration loop of the PBL scheme.

  call trb_mtn_stress_tend(state, pbuf, cam_in)

  call pbuf_get_field(pbuf, ksrftms_idx, ksrftms)
  call pbuf_get_field(pbuf, tautmsx_idx, tautmsx)
  call pbuf_get_field(pbuf, tautmsy_idx, tautmsy)

  tautotx(:ncol) = cam_in%wsx(:ncol) + tautmsx(:ncol)
  tautoty(:ncol) = cam_in%wsy(:ncol) + tautmsy(:ncol)

  ! ------------------------------------- !
  ! Computation of Beljaars SGO form drag !
  ! ------------------------------------- !

  call beljaars_drag_tend(state, pbuf, cam_in)

  call pbuf_get_field(pbuf, dragblj_idx, dragblj)
  call pbuf_get_field(pbuf, taubljx_idx, taubljx)
  call pbuf_get_field(pbuf, taubljy_idx, taubljy)

  ! Add Beljaars integrated drag

  tautotx(:ncol) = tautotx(:ncol) + taubljx(:ncol)
  tautoty(:ncol) = tautoty(:ncol) + taubljy(:ncol)

  !----------------------------------------------------------------------- !
  !   Computation of eddy diffusivities - Select appropriate PBL scheme    !
  !----------------------------------------------------------------------- !
  call pbuf_get_field(pbuf, kvm_idx,  kvm_in)
  call pbuf_get_field(pbuf, kvh_idx,  kvh_in)
  call pbuf_get_field(pbuf, smaw_idx, smaw)
  call pbuf_get_field(pbuf, tke_idx,  tke)

  ! Get potential temperature.
  th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

  select case (eddy_scheme)
  case ( 'diag_TKE', 'SPCAM_m2005' )

     call eddy_diff_tend(state, pbuf, cam_in, &
          ztodt, p, tint, rhoi, cldn, wstarent, &
          kvm_in, kvh_in, ksrftms, dragblj, tauresx, tauresy, &
          rrho, ustar, pblh, kvm, kvh, kvq, cgh, cgs, tpert, qpert, &
          tke, sprod, sfi, turbtype, smaw)

     ! The diag_TKE scheme does not calculate the Monin-Obukhov length, which is used in dry deposition calculations.
     ! Use the routines from pbl_utils to accomplish this. Assumes ustar and rrho have been set.
     call virtem(ncol, th(:ncol,pver),state%q(:ncol,pver,1), thvs(:ncol))
     call calc_obklen(ncol, th(:ncol,pver), thvs(:ncol), cam_in%cflx(:ncol,1), &
          cam_in%shf(:ncol), rrho(:ncol), ustar(:ncol), &
          khfs(:ncol),    kqfs(:ncol), kbfs(:ncol),   obklen(:ncol))


  case ( 'HB', 'HBR', 'SPCAM_sam1mom' )

     ! Modification : We may need to use 'taux' instead of 'tautotx' here, for
     !                consistency with the previous HB scheme.


     call compute_hb_diff( lchnk     , ncol     ,                                &
          th        , state%t  , state%q , state%zm , state%zi, &
          state%pmid, state%u  , state%v , tautotx  , tautoty , &
          cam_in%shf, cam_in%cflx(:,1), obklen  , ustar    , pblh    , &
          kvm       , kvh      , kvq     , cgh      , cgs     , &
          tpert     , qpert    , cldn    , cam_in%ocnfrac  , tke     , &
          ri        , &
          eddy_scheme )

     call outfld( 'HB_ri',          ri,         pcols,   lchnk )

  case ( 'CLUBB_SGS' )

    !
    ! run HB scheme where CLUBB is not active when running cam_dev
    !
    if (cam_physpkg_is("cam_dev")) then
      call compute_hb_diff( lchnk     , ncol     ,                                &
           th        , state%t  , state%q , state%zm , state%zi, &
           state%pmid, state%u  , state%v , tautotx  , tautoty , &
           cam_in%shf, cam_in%cflx(:,1), obklen  , ustar    , pblh    , &
           kvm       , kvh      , kvq     , cgh      , cgs     , &
           tpert     , qpert    , cldn    , cam_in%ocnfrac  , tke     , &
           ri        , &
           eddy_scheme )
      clubbtop_idx = pbuf_get_index('clubbtop')
      call pbuf_get_field(pbuf, clubbtop_idx, clubbtop)

      do i=1,ncol
        do k=clubbtop(i),pverp
          kvm(i,k) = 0.0_r8
          kvh(i,k) = 0.0_r8
          kvq(i,k) = 0.0_r8
          cgs(i,k) = 0.0_r8
          cgh(i,k) = 0.0_r8
        end do
      end do

      call outfld( 'HB_ri',          ri,         pcols,   lchnk )
    else
      ! CLUBB has only a bare-bones placeholder here. If using CLUBB, the
      ! PBL diffusion will happen before coupling, so vertical_diffusion
      ! is only handling other things, e.g. some boundary conditions, tms,
      ! and molecular diffusion.
      
      call virtem(ncol, th(:ncol,pver),state%q(:ncol,pver,1), thvs(:ncol))
      
      call calc_ustar( ncol, state%t(:ncol,pver), state%pmid(:ncol,pver), &
           cam_in%wsx(:ncol), cam_in%wsy(:ncol), rrho(:ncol), ustar(:ncol))
      ! Use actual qflux, not lhf/latvap as was done previously
      call calc_obklen( ncol, th(:ncol,pver), thvs(:ncol), cam_in%cflx(:ncol,1), &
           cam_in%shf(:ncol), rrho(:ncol), ustar(:ncol),  &
           khfs(:ncol), kqfs(:ncol), kbfs(:ncol), obklen(:ncol))
      ! These tendencies all applied elsewhere.
      kvm = 0._r8
      kvh = 0._r8
      kvq = 0._r8
      ! Not defined since PBL is not actually running here.
      cgh = 0._r8
      cgs = 0._r8
    end if
  end select

  call outfld( 'ustar',   ustar(:), pcols, lchnk )
  call outfld( 'obklen', obklen(:), pcols, lchnk )

  ! kvh (in pbuf) is used by other physics parameterizations, and as an initial guess in compute_eddy_diff
  ! on the next timestep.  It is not updated by the compute_vdiff call below.
  call pbuf_set_field(pbuf, kvh_idx, kvh)

  ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
  ! The contributions for molecular diffusion made to kvm by the call to compute_vdiff below
  ! are not included in the pbuf as these are not needed in the initial guess by compute_eddy_diff.
  call pbuf_set_field(pbuf, kvm_idx, kvm)

  !------------------------------------ !
  !    Application of diffusivities     !
  !------------------------------------ !

  ! Set arrays from input state.
  q_tmp(:ncol,:,:) = state%q(:ncol,:,:)
  s_tmp(:ncol,:) = state%s(:ncol,:)
  u_tmp(:ncol,:) = state%u(:ncol,:)
  v_tmp(:ncol,:) = state%v(:ncol,:)

  !------------------------------------------------------ !
  ! Write profile output before applying diffusion scheme !
  !------------------------------------------------------ !

  if (.not. do_pbl_diags) then
     sl_prePBL(:ncol,:pver)  = s_tmp(:ncol,:) -   latvap * q_tmp(:ncol,:,ixcldliq) &
          - ( latvap + latice) * q_tmp(:ncol,:,ixcldice)
     qt_prePBL(:ncol,:pver)  = q_tmp(:ncol,:,1) + q_tmp(:ncol,:,ixcldliq) &
          + q_tmp(:ncol,:,ixcldice)
     slv_prePBL(:ncol,:pver) = sl_prePBL(:ncol,:pver) * ( 1._r8 + zvir*qt_prePBL(:ncol,:pver) )

     do k = 1, pver
        call qsat(state%t(1:ncol,k), state%pmid(1:ncol,k), tem2(1:ncol,k), ftem(1:ncol,k), ncol)
     end do
     ftem_prePBL(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8

     call outfld( 'qt_pre_PBL   ', qt_prePBL,                 pcols, lchnk )
     call outfld( 'sl_pre_PBL   ', sl_prePBL,                 pcols, lchnk )
     call outfld( 'slv_pre_PBL  ', slv_prePBL,                pcols, lchnk )
     call outfld( 'u_pre_PBL    ', state%u,                   pcols, lchnk )
     call outfld( 'v_pre_PBL    ', state%v,                   pcols, lchnk )
     call outfld( 'qv_pre_PBL   ', state%q(:ncol,:,1),        pcols, lchnk )
     call outfld( 'ql_pre_PBL   ', state%q(:ncol,:,ixcldliq), pcols, lchnk )
     call outfld( 'qi_pre_PBL   ', state%q(:ncol,:,ixcldice), pcols, lchnk )
     call outfld( 't_pre_PBL    ', state%t,                   pcols, lchnk )
     call outfld( 'rh_pre_PBL   ', ftem_prePBL,               pcols, lchnk )

  end if

  ! --------------------------------------------------------------------------------- !
  ! Call the diffusivity solver and solve diffusion equation                          !
  ! The final two arguments are optional function references to                       !
  ! constituent-independent and constituent-dependent moleculuar diffusivity routines !
  ! --------------------------------------------------------------------------------- !

  ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
  !                separately print out as diagnostic output, because these are different from
  !                the explicit 'tautotx, tautoty' computed above.
  ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

  call pbuf_get_field(pbuf, kvt_idx, kvt)

  if (do_molec_diff .and. .not. waccmx_mode) then
     ! Top boundary condition for dry static energy
     dse_top(:ncol) = cpairv(:ncol,1,lchnk) * tint(:ncol,1) + &
          gravit * state%zi(:ncol,1)
  else
     dse_top(:ncol) = 0._r8
  end if

  select case (eddy_scheme)
  case ('CLUBB_SGS')
     ! CLUBB applies some fluxes itself, but we still want constituent
     ! fluxes applied here (except water vapor).
     taux = 0._r8
     tauy = 0._r8
     shflux = 0._r8
     cflux(:,1) = 0._r8
     if (cam_physpkg_is("cam_dev")) then
       ! surface fluxes applied in clubb emissions module
       cflux(:,2:) = 0._r8
     else
       cflux(:,2:) = cam_in%cflx(:,2:)
     end if
  case default
     taux = cam_in%wsx
     tauy = cam_in%wsy
     shflux = cam_in%shf
     cflux = cam_in%cflx
  end select

  if( any(fieldlist_wet) ) then

     if (do_molec_diff) then
        call compute_molec_diff(state%lchnk, pcols, pver, pcnst, ncol, &
             kvm, kvt, tint, rhoi, kq_scal, cnst_mw, &
             mw_fac, nbot_molec)
     end if

     call compute_vdiff( state%lchnk   ,                                                                     &
          pcols         , pver               , pcnst        , ncol          , tint          , &
          p    , state%t      , rhoi, ztodt         , taux          , &
          tauy          , shflux             , cflux        , &
          kvh           , kvm                , kvq          , cgs           , cgh           , &
          state%zi      , ksrftms            , dragblj      , &
          qmincg       , fieldlist_wet , fieldlist_molec,&
          u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
          tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
          tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), dse_top, &
          do_molec_diff, waccmx_mode, &
          vd_lu_qdecomp, &
          ubc_mmr, ubc_flux, kvt, state%pmid, &
          cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx, nbot_molec, &
          kq_scal, mw_fac)

     call handle_errmsg(errstring, subname="compute_vdiff", &
          extra_msg="Error in fieldlist_wet call from vertical_diffusion.")

  end if

  if( any( fieldlist_dry ) ) then

     if( do_molec_diff ) then
        ! kvm is unused in the output here (since it was assigned
        ! above), so we use a temp kvm for the inout argument, and
        ! ignore the value output by compute_molec_diff.
        kvm_temp = kvm
        call compute_molec_diff(state%lchnk, pcols, pver, pcnst, ncol, &
             kvm_temp, kvt, tint, rhoi_dry, kq_scal, cnst_mw, &
             mw_fac, nbot_molec)
     end if

     call compute_vdiff( state%lchnk   ,                                                                     &
          pcols         , pver               , pcnst        , ncol          , tint          , &
          p_dry , state%t      , rhoi_dry,  ztodt         , taux          , &
          tauy          , shflux             , cflux        , &
          kvh           , kvm                , kvq          , cgs           , cgh           , &
          state%zi      , ksrftms            , dragblj      , &
          qmincg       , fieldlist_dry , fieldlist_molec,&
          u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
          tautmsx_temp  , tautmsy_temp       , dtk_temp     , topflx_temp   , errstring     , &
          tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), dse_top, &
          do_molec_diff , waccmx_mode, &
          vd_lu_qdecomp, &
          ubc_mmr, ubc_flux, kvt, state%pmiddry, &
          cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx, nbot_molec, &
          kq_scal, mw_fac)

     call handle_errmsg(errstring, subname="compute_vdiff", &
          extra_msg="Error in fieldlist_dry call from vertical_diffusion.")

  end if

  if (prog_modal_aero) then

     ! Modal aerosol species not diffused, so just add the explicit surface fluxes to the
     ! lowest layer.  **NOTE** This code assumes wet mmr.

     tmp1(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
     do m = 1, pmam_ncnst
        l = pmam_cnst_idx(m)
        q_tmp(:ncol,pver,l) = q_tmp(:ncol,pver,l) + tmp1(:ncol) * cflux(:ncol,l)
     enddo
  end if

  ! -------------------------------------------------------- !
  ! Diagnostics and output writing after applying PBL scheme !
  ! -------------------------------------------------------- !

  if (.not. do_pbl_diags) then

     sl(:ncol,:pver)  = s_tmp(:ncol,:) -   latvap           * q_tmp(:ncol,:,ixcldliq) &
          - ( latvap + latice) * q_tmp(:ncol,:,ixcldice)
     qt(:ncol,:pver)  = q_tmp(:ncol,:,1) + q_tmp(:ncol,:,ixcldliq) &
          + q_tmp(:ncol,:,ixcldice)
     slv(:ncol,:pver) = sl(:ncol,:pver) * ( 1._r8 + zvir*qt(:ncol,:pver) )

     slflx(:ncol,1) = 0._r8
     qtflx(:ncol,1) = 0._r8
     uflx(:ncol,1)  = 0._r8
     vflx(:ncol,1)  = 0._r8

     slflx_cg(:ncol,1) = 0._r8
     qtflx_cg(:ncol,1) = 0._r8
     uflx_cg(:ncol,1)  = 0._r8
     vflx_cg(:ncol,1)  = 0._r8

     do k = 2, pver
        do i = 1, ncol
           rhoair     = state%pint(i,k) / ( rair * ( ( 0.5_r8*(slv(i,k)+slv(i,k-1)) - gravit*state%zi(i,k))/cpair ) )
           slflx(i,k) = kvh(i,k) * &
                ( - rhoair*(sl(i,k-1)-sl(i,k))/(state%zm(i,k-1)-state%zm(i,k)) &
                + cgh(i,k) )
           qtflx(i,k) = kvh(i,k) * &
                ( - rhoair*(qt(i,k-1)-qt(i,k))/(state%zm(i,k-1)-state%zm(i,k)) &
                + rhoair*(cam_in%cflx(i,1)+cam_in%cflx(i,ixcldliq)+cam_in%cflx(i,ixcldice))*cgs(i,k) )
           uflx(i,k)  = kvm(i,k) * &
                ( - rhoair*(u_tmp(i,k-1)-u_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
           vflx(i,k)  = kvm(i,k) * &
                ( - rhoair*(v_tmp(i,k-1)-v_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
           slflx_cg(i,k) = kvh(i,k) * cgh(i,k)
           qtflx_cg(i,k) = kvh(i,k) * rhoair * ( cam_in%cflx(i,1) + &
                cam_in%cflx(i,ixcldliq) + cam_in%cflx(i,ixcldice) ) * cgs(i,k)
           uflx_cg(i,k)  = 0._r8
           vflx_cg(i,k)  = 0._r8
        end do
     end do

     ! Modification : I should check whether slflx(:ncol,pverp) is correctly computed.
     !                Note also that 'tautotx' is explicit total stress, different from
     !                the ones that have been actually added into the atmosphere.

     slflx(:ncol,pverp) = cam_in%shf(:ncol)
     qtflx(:ncol,pverp) = cam_in%cflx(:ncol,1)
     uflx(:ncol,pverp)  = tautotx(:ncol)
     vflx(:ncol,pverp)  = tautoty(:ncol)

     slflx_cg(:ncol,pverp) = 0._r8
     qtflx_cg(:ncol,pverp) = 0._r8
     uflx_cg(:ncol,pverp)  = 0._r8
     vflx_cg(:ncol,pverp)  = 0._r8

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

  ! --------------------------------------------------------------- !
  ! Convert the new profiles into vertical diffusion tendencies.    !
  ! Convert KE dissipative heat change into "temperature" tendency. !
  ! --------------------------------------------------------------- !
  ! All variables are modified by vertical diffusion

  lq(:) = .TRUE.
  call physics_ptend_init(ptend,state%psetcols, "vertical diffusion", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  ptend%s(:ncol,:)       = ( s_tmp(:ncol,:) - state%s(:ncol,:) ) * rztodt
  ptend%u(:ncol,:)       = ( u_tmp(:ncol,:) - state%u(:ncol,:) ) * rztodt
  ptend%v(:ncol,:)       = ( v_tmp(:ncol,:) - state%v(:ncol,:) ) * rztodt
  ptend%q(:ncol,:pver,:) = ( q_tmp(:ncol,:pver,:) - state%q(:ncol,:pver,:) ) * rztodt

  ! Convert tendencies of dry constituents to dry basis.
  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        ptend%q(:ncol,:pver,m) = ptend%q(:ncol,:pver,m)*state%pdel(:ncol,:pver)/state%pdeldry(:ncol,:pver)
     endif
  end do
  ! convert wet mmr back to dry before conservation check
  call set_wet_to_dry(state)

  if (.not. do_pbl_diags) then
     slten(:ncol,:)         = ( sl(:ncol,:) - sl_prePBL(:ncol,:) ) * rztodt
     qtten(:ncol,:)         = ( qt(:ncol,:) - qt_prePBL(:ncol,:) ) * rztodt
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

    if( (eddy_scheme .eq. 'diag_TKE' .or. eddy_scheme .eq. 'SPCAM_m2005') .and. do_pseudocon_diff ) then

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

     call positive_moisture( cpair, latvap, latvap+latice, ncol, pver, ztodt, qmin(1), qmin(ixcldliq), qmin(ixcldice),    &
          state%pdel(:ncol,pver:1:-1), qv_pro(:ncol,pver:1:-1), ql_pro(:ncol,pver:1:-1), &
          qi_pro(:ncol,pver:1:-1), t_pro(:ncol,pver:1:-1), s_pro(:ncol,pver:1:-1),       &
          ptend%q(:ncol,pver:1:-1,1), ptend%q(:ncol,pver:1:-1,ixcldliq),                 &
          ptend%q(:ncol,pver:1:-1,ixcldice), ptend%s(:ncol,pver:1:-1) )

  end if

  ! ----------------------------------------------------------------- !
  ! Re-calculate diagnostic output variables after vertical diffusion !
  ! ----------------------------------------------------------------- !

  if (.not. do_pbl_diags) then

     qv_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,1)         + ptend%q(:ncol,:pver,1)        * ztodt
     ql_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,ixcldliq)  + ptend%q(:ncol,:pver,ixcldliq) * ztodt
     qi_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,ixcldice)  + ptend%q(:ncol,:pver,ixcldice) * ztodt
     s_aft_PBL(:ncol,:pver)   =   state%s(:ncol,:pver)           + ptend%s(:ncol,:pver)          * ztodt
     t_aftPBL(:ncol,:pver)    = ( s_aft_PBL(:ncol,:pver) - gravit*state%zm(:ncol,:pver) ) / cpair

     u_aft_PBL(:ncol,:pver)   =  state%u(:ncol,:pver)          + ptend%u(:ncol,:pver)            * ztodt
     v_aft_PBL(:ncol,:pver)   =  state%v(:ncol,:pver)          + ptend%v(:ncol,:pver)            * ztodt

     do k = 1, pver
        call qsat(t_aftPBL(1:ncol,k), state%pmid(1:ncol,k), tem2(1:ncol,k), ftem(1:ncol,k), ncol)
     end do
     ftem_aftPBL(:ncol,:pver) = qv_aft_PBL(:ncol,:pver) / ftem(:ncol,:pver) * 100._r8

     tten(:ncol,:pver)        = ( t_aftPBL(:ncol,:pver)    - state%t(:ncol,:pver) )              * rztodt
     rhten(:ncol,:pver)       = ( ftem_aftPBL(:ncol,:pver) - ftem_prePBL(:ncol,:pver) )          * rztodt

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

  ! -------------------------------------------------------------- !
  ! Writing state variables after PBL scheme for detailed analysis !
  ! -------------------------------------------------------------- !

  if (.not. do_pbl_diags) then

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
     call outfld( 'uten_PBL'     , ptend%u(:ncol,:),          pcols, lchnk )
     call outfld( 'vten_PBL'     , ptend%v(:ncol,:),          pcols, lchnk )
     call outfld( 'qvten_PBL'    , ptend%q(:ncol,:,1),        pcols, lchnk )
     call outfld( 'qlten_PBL'    , ptend%q(:ncol,:,ixcldliq), pcols, lchnk )
     call outfld( 'qiten_PBL'    , ptend%q(:ncol,:,ixcldice), pcols, lchnk )
     call outfld( 'tten_PBL'     , tten,                      pcols, lchnk )
     call outfld( 'rhten_PBL'    , rhten,                     pcols, lchnk )

  end if

  ! ------------------------------------------- !
  ! Writing the other standard output variables !
  ! ------------------------------------------- !

  if (.not. do_pbl_diags) then
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
  end if
  call outfld( 'USTAR'        , ustar,                     pcols, lchnk )
  call outfld( 'KVH'          , kvh,                       pcols, lchnk )
  call outfld( 'KVT'          , kvt,                       pcols, lchnk )
  call outfld( 'KVM'          , kvm,                       pcols, lchnk )
  call outfld( 'CGS'          , cgs,                       pcols, lchnk )
  dtk(:ncol,:) = dtk(:ncol,:) / cpair              ! Normalize heating for history
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
  implicit none
  integer,  intent(in)     :: ncol, mkx
  real(r8), intent(in)     :: cp, xlv, xls
  real(r8), intent(in)     :: dt, qvmin, qlmin, qimin
  real(r8), intent(in)     :: dp(ncol,mkx)
  real(r8), intent(inout)  :: qv(ncol,mkx), ql(ncol,mkx), qi(ncol,mkx), t(ncol,mkx), s(ncol,mkx)
  real(r8), intent(inout)  :: qvten(ncol,mkx), qlten(ncol,mkx), qiten(ncol,mkx), sten(ncol,mkx)
  integer   i, k
  real(r8)  dql, dqi, dqv, sum, aa, dum

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
  return

end subroutine positive_moisture

end module vertical_diffusion
