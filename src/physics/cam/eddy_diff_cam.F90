module eddy_diff_cam

use shr_kind_mod, only: i4 => shr_kind_i4, r8 => shr_kind_r8
use ppgrid, only: pcols, pver, pverp
use cam_logfile, only: iulog
use cam_abortutils, only: endrun
use physconst, only: gravit, cpair, rair, zvir, latvap, latice, karman
use diffusion_solver, only: vdiff_selector
use eddy_diff, only: ncvmax
use time_manager, only: is_first_step
use physics_buffer, only: physics_buffer_desc
use spmd_utils, only: masterproc
use phys_control, only: phys_getopts

implicit none
private

public :: eddy_diff_readnl
public :: eddy_diff_register
public :: eddy_diff_init
public :: eddy_diff_tend

! Is UNICON switched on (and thus interacting with eddy_diff via pbuf)?
logical :: unicon_is_on

! Number of iterations for solution
integer, parameter :: nturb = 5

! Logical switches for moist mixing ratio diffusion
type(vdiff_selector) :: fieldlist_wet
! Logical switches for molecular diffusion
! (Molecular diffusion is not done here.)
type(vdiff_selector) :: fieldlist_molec

integer :: ntop_eddy, nbot_eddy

! Cloud mass constituent indices
integer :: ixcldliq, ixcldice

! input pbuf field indices
integer :: qrl_idx   = -1
integer :: wsedl_idx = -1

! output pbuf field indices for UNICON
integer :: bprod_idx    = -1
integer :: ipbl_idx     = -1
integer :: kpblh_idx    = -1
integer :: wstarPBL_idx = -1
integer :: tkes_idx     = -1
integer :: went_idx     = -1

! Mixing lengths squared.
! Used for computing free air diffusivity.
real(r8) :: ml2(pver+1)

! Various namelist options to limit or tweak the effects of eddy diffusion.

! Pressure defining the bottom of the upper atmosphere for kvh scaling (Pa)
real(r8) :: kv_top_pressure = 0._r8
! Eddy diffusivity scale factor for upper atmosphere
real(r8) :: kv_top_scale = 1._r8
! Eddy diffusivity scale factor for the free troposphere
real(r8) :: kv_freetrop_scale = 1._r8

! The following all have to be set in all cases.
real(r8), parameter  :: unset_r8 = huge(1._r8)
! Maximum master length for diag_TKE
real(r8) :: eddy_lbulk_max = unset_r8
! Maximum dissipation length for diag_TKE
real(r8) :: eddy_leng_max = unset_r8
! Bottom pressure level (hPa) for eddy_leng_max
real(r8) :: eddy_max_bot_pressure = unset_r8
! Moist entrainment enhancement param
real(r8) :: eddy_moist_entrain_a2l = unset_r8

contains

subroutine eddy_diff_readnl(nlfile)
  use namelist_utils, only: find_group_name
  use units, only: getunit, freeunit
  use spmd_utils, only: masterprocid, mpi_real8, mpicom
  use shr_log_mod, only: errMsg => shr_log_errMsg

  ! filepath for file containing namelist input
  character(len=*), intent(in) :: nlfile

  ! file unit and error code
  integer :: unitn, ierr

  character(len=*), parameter :: subname = 'eddy_diff_readnl'

  namelist /eddy_diff_nl/ kv_top_pressure, kv_top_scale, &
       kv_freetrop_scale, eddy_lbulk_max, eddy_leng_max, &
       eddy_max_bot_pressure, eddy_moist_entrain_a2l

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'eddy_diff_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, eddy_diff_nl, iostat=ierr)
     end if
     if (ierr /= 0) then
        call endrun(subname // ':: ERROR reading namelist')
     end if
     close(unitn)
     call freeunit(unitn)
  end if

  call mpi_bcast(kv_top_pressure,   1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(kv_top_scale,      1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(kv_freetrop_scale, 1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")

  call mpi_bcast(eddy_lbulk_max,    1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(eddy_leng_max,     1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(eddy_max_bot_pressure,     1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")
  call mpi_bcast(eddy_moist_entrain_a2l,    1,   mpi_real8, masterprocid, mpicom, ierr)
  if (ierr /= 0) call endrun(errMsg(__FILE__, __LINE__)//" mpi_bcast error")

end subroutine eddy_diff_readnl

subroutine eddy_diff_register()
  use physics_buffer, only: pbuf_add_field, dtype_r8, dtype_i4

  character(len=16) :: shallow_scheme

  ! Check for UNICON and add relevant pbuf entries.
  call phys_getopts(shallow_scheme_out=shallow_scheme)

  unicon_is_on = (shallow_scheme == "UNICON")

  if (unicon_is_on) then
     call pbuf_add_field('bprod',    'global', dtype_r8, (/pcols,pverp/), bprod_idx)
     call pbuf_add_field('ipbl',     'global', dtype_i4, (/pcols/),       ipbl_idx)
     call pbuf_add_field('kpblh',    'global', dtype_i4, (/pcols/),       kpblh_idx)
     call pbuf_add_field('wstarPBL', 'global', dtype_r8, (/pcols/),       wstarPBL_idx)
     call pbuf_add_field('tkes',     'global', dtype_r8, (/pcols/),       tkes_idx)
     call pbuf_add_field('went',     'global', dtype_r8, (/pcols/),       went_idx)
  end if

end subroutine eddy_diff_register

subroutine eddy_diff_init(pbuf2d, ntop_eddy_in, nbot_eddy_in)

  use error_messages, only: handle_errmsg
  use cam_history, only: addfld, add_default, horiz_only
  use constituents, only: cnst_get_ind
  use ref_pres, only: pref_mid
  use diffusion_solver, only: new_fieldlist_vdiff, vdiff_select
  use eddy_diff, only: init_eddy_diff
  use physics_buffer, only: pbuf_set_field, pbuf_get_index

  type(physics_buffer_desc), pointer :: pbuf2d(:,:) ! Physics buffer
  integer,  intent(in) :: ntop_eddy_in ! Top interface level to which eddy vertical diffusivity is applied ( = 1 )
  integer,  intent(in) :: nbot_eddy_in ! Bottom interface level to which eddy vertical diffusivity is applied ( = pver )

  character(len=128) :: errstring

  real(r8) :: leng_max(pver)
  integer :: k

  logical :: history_amwg

  ntop_eddy = ntop_eddy_in
  nbot_eddy = nbot_eddy_in

  do k = 1, pver
     if (pref_mid(k) <= eddy_max_bot_pressure*1.e2_r8) then
        leng_max(k) = eddy_leng_max
     else
        leng_max(k) = 40.e3_r8
     end if
  end do

  if (masterproc) then
     write(iulog,*)'init_eddy_diff: nturb=',nturb
     write(iulog,*)'init_eddy_diff: eddy_leng_max=',eddy_leng_max,' lbulk_max=',eddy_lbulk_max
     do k = 1,pver
        write(iulog,*)'init_eddy_diff:',k,pref_mid(k),'leng_max=',leng_max(k)
     end do
  end if

  call init_eddy_diff(pver, gravit, cpair, rair, zvir, &
       latvap, latice, ntop_eddy, nbot_eddy, karman, &
       eddy_lbulk_max, leng_max, &
       eddy_moist_entrain_a2l, errstring)

  call handle_errmsg(errstring, subname="init_eddy_diff")

  ! Set the square of the mixing lengths.
  ml2(1:ntop_eddy) = 0._r8
  do k = ntop_eddy + 1, nbot_eddy
     ml2(k) = 30.0_r8**2
  end do
  ml2(nbot_eddy+1:pver+1) = 0._r8

  ! Get fieldlists to pass to diffusion solver.
  fieldlist_wet   = new_fieldlist_vdiff(1)
  fieldlist_molec = new_fieldlist_vdiff(1)

  call handle_errmsg(vdiff_select(fieldlist_wet,'s'), &
       subname="vdiff_select")
  call handle_errmsg(vdiff_select(fieldlist_wet,'q',1), &
       subname="vdiff_select")
  call handle_errmsg(vdiff_select(fieldlist_wet,'u'), &
       subname="vdiff_select")
  call handle_errmsg(vdiff_select(fieldlist_wet,'v'), &
       subname="vdiff_select")

  ! Cloud mass constituents
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)

  ! Input pbuf fields
  qrl_idx   = pbuf_get_index('QRL')
  wsedl_idx = pbuf_get_index('WSEDL')

  ! Initialize output pbuf fields
  if (is_first_step() .and. unicon_is_on) then
     call pbuf_set_field(pbuf2d, bprod_idx,    1.0e-5_r8)
     call pbuf_set_field(pbuf2d, ipbl_idx,     0    )
     call pbuf_set_field(pbuf2d, kpblh_idx,    1    )
     call pbuf_set_field(pbuf2d, wstarPBL_idx, 0.0_r8)
     call pbuf_set_field(pbuf2d, tkes_idx,     0.0_r8)
     call pbuf_set_field(pbuf2d, went_idx,     0.0_r8)
  end if

  ! Scheme-specific default output.
  call phys_getopts(history_amwg_out=history_amwg)

  call addfld('WGUSTD', horiz_only,  'A',          'm/s',  'wind gusts from turbulence'                            )
  if (history_amwg) then
     call add_default( 'WGUSTD  ', 1, ' ' )
  end if

  ! ------------------------------------------------------------------- !
  ! Writing outputs for detailed analysis of UW moist turbulence scheme !
  ! ------------------------------------------------------------------- !

  call addfld( 'BPROD',   ['ilev'],  'A',        'm2/s3',  'Buoyancy Production'                                   )
  call addfld( 'SFI',     ['ilev'],  'A',            '1',  'Interface-layer sat frac'                              )
  call addfld( 'SPROD',   ['ilev'],  'A',        'm2/s3',  'Shear Production'                                      )


  call addfld('UW_errorPBL',horiz_only,'A',         'm2/s',  'Error function of UW PBL')
  call addfld('UW_n2',        ['lev'], 'A',          's-2',  'Buoyancy Frequency, LI')
  call addfld('UW_s2',        ['lev'], 'A',          's-2',  'Shear Frequency, LI')
  call addfld('UW_ri',        ['lev'], 'A',            '1',  'Interface Richardson Number, I')
  call addfld('UW_sfuh',      ['lev'], 'A',            '1',  'Upper-Half Saturation Fraction, L')
  call addfld('UW_sflh',      ['lev'], 'A',            '1',  'Lower-Half Saturation Fraction, L')
  call addfld('UW_sfi',      ['ilev'], 'A',            '1',  'Interface Saturation Fraction, I')
  call addfld('UW_cldn',      ['lev'], 'A',            '1',  'Cloud Fraction, L')
  call addfld('UW_qrl',       ['lev'], 'A', 'gravity W/m2',  'LW cooling rate, L')
  call addfld('UW_ql',        ['lev'], 'A',        'kg/kg',  'ql(LWC), L')
  call addfld('UW_chu',      ['ilev'], 'A', 'gravity kg/J',  'Buoyancy Coefficient, chu, I')
  call addfld('UW_chs',      ['ilev'], 'A', 'gravity kg/J',  'Buoyancy Coefficient, chs, I')
  call addfld('UW_cmu',      ['ilev'], 'A','gravity/kg/kg',  'Buoyancy Coefficient, cmu, I')
  call addfld('UW_cms',      ['ilev'], 'A','gravity/kg/kg',  'Buoyancy Coefficient, cms, I')
  call addfld('UW_tke',      ['ilev'], 'A',        'm2/s2',  'TKE, I')
  call addfld('UW_wcap',     ['ilev'], 'A',        'm2/s2',  'Wcap, I')
  call addfld('UW_bprod',    ['ilev'], 'A',        'm2/s3',  'Buoyancy production, I')
  call addfld('UW_sprod',    ['ilev'], 'A',        'm2/s3',  'Shear production, I')
  call addfld('UW_kvh',      ['ilev'], 'A',         'm2/s',  'Eddy diffusivity of heat, I')
  call addfld('UW_kvm',      ['ilev'], 'A',         'm2/s',  'Eddy diffusivity of uv, I')
  call addfld('UW_pblh',   horiz_only, 'A',            'm',  'PBLH, 1')
  call addfld('UW_pblhp',  horiz_only, 'A',           'Pa',  'PBLH pressure, 1')
  call addfld('UW_tpert',  horiz_only, 'A',            'K',  'Convective T excess, 1')
  call addfld('UW_qpert',  horiz_only, 'A',        'kg/kg',  'Convective qt excess, I')
  call addfld('UW_wpert',  horiz_only, 'A',          'm/s',  'Convective W excess, I')
  call addfld('UW_ustar',  horiz_only, 'A',          'm/s',  'Surface Frictional Velocity, 1')
  call addfld('UW_tkes',   horiz_only, 'A',        'm2/s2',  'Surface TKE, 1')
  call addfld('UW_minpblh',horiz_only, 'A',            'm',  'Minimum PBLH, 1')
  call addfld('UW_turbtype', ['ilev'], 'A',            '1',  'Interface Turbulence Type, I')
  call addfld('UW_kbase_o',   ['lev'], 'A',            '1',  'Initial CL Base Exterbal Interface Index, CL')
  call addfld('UW_ktop_o',    ['lev'], 'A',            '1',  'Initial Top Exterbal Interface Index, CL')
  call addfld('UW_ncvfin_o',horiz_only,'A',            '1',  'Initial Total Number of CL regimes, CL')
  call addfld('UW_kbase_mg',  ['lev'], 'A',            '1',  'kbase after merging, CL')
  call addfld('UW_ktop_mg',   ['lev'], 'A',            '1',  'ktop after merging, CL')
  call addfld('UW_ncvfin_mg',horiz_only,'A',           '1',  'ncvfin after merging, CL')
  call addfld('UW_kbase_f',   ['lev'], 'A',            '1',  'Final kbase with SRCL, CL')
  call addfld('UW_ktop_f',    ['lev'], 'A',            '1',  'Final ktop with SRCL, CL')
  call addfld('UW_ncvfin_f',horiz_only,'A',            '1',  'Final ncvfin with SRCL, CL')
  call addfld('UW_wet',       ['lev'], 'A',          'm/s',  'Entrainment rate at CL top, CL')
  call addfld('UW_web',       ['lev'], 'A',          'm/s',  'Entrainment rate at CL base, CL')
  call addfld('UW_jtbu',      ['lev'], 'A',         'm/s2',  'Buoyancy jump across CL top, CL')
  call addfld('UW_jbbu',      ['lev'], 'A',         'm/s2',  'Buoyancy jump across CL base, CL')
  call addfld('UW_evhc',      ['lev'], 'A',            '1',  'Evaporative enhancement factor, CL')
  call addfld('UW_jt2slv',    ['lev'], 'A',         'J/kg',  'slv jump for evhc, CL')
  call addfld('UW_n2ht',      ['lev'], 'A',          's-2',  'n2 at just below CL top interface, CL')
  call addfld('UW_n2hb',      ['lev'], 'A',          's-2',  'n2 at just above CL base interface')
  call addfld('UW_lwp',       ['lev'], 'A',        'kg/m2',  'LWP in the CL top layer, CL')
  call addfld('UW_optdepth',  ['lev'], 'A',            '1',  'Optical depth of the CL top layer, CL')
  call addfld('UW_radfrac',   ['lev'], 'A',            '1',  'Fraction of radiative cooling confined in the CL top')
  call addfld('UW_radf',      ['lev'], 'A',        'm2/s3',  'Buoyancy production at the CL top by radf, I')
  call addfld('UW_wstar',     ['lev'], 'A',          'm/s',  'Convective velocity, Wstar, CL')
  call addfld('UW_wstar3fact',['lev'], 'A',            '1',  'Enhancement of wstar3 due to entrainment, CL')
  call addfld('UW_ebrk',      ['lev'], 'A',        'm2/s2',  'CL-averaged TKE, CL')
  call addfld('UW_wbrk',      ['lev'], 'A',        'm2/s2',  'CL-averaged W, CL')
  call addfld('UW_lbrk',      ['lev'], 'A',            'm',  'CL internal thickness, CL')
  call addfld('UW_ricl',      ['lev'], 'A',            '1',  'CL-averaged Ri, CL')
  call addfld('UW_ghcl',      ['lev'], 'A',            '1',  'CL-averaged gh, CL')
  call addfld('UW_shcl',      ['lev'], 'A',            '1',  'CL-averaged sh, CL')
  call addfld('UW_smcl',      ['lev'], 'A',            '1',  'CL-averaged sm, CL')
  call addfld('UW_gh',       ['ilev'], 'A',            '1',  'gh at all interfaces, I')
  call addfld('UW_sh',       ['ilev'], 'A',            '1',  'sh at all interfaces, I')
  call addfld('UW_sm',       ['ilev'], 'A',            '1',  'sm at all interfaces, I')
  call addfld('UW_ria',      ['ilev'], 'A',            '1',  'ri at all interfaces, I')
  call addfld('UW_leng',     ['ilev'], 'A',          'm/s',  'Turbulence length scale, I')
  ! For sedimentation-entrainment feedback analysis
  call addfld('UW_wsed',      ['lev'], 'A',          'm/s',  'Sedimentation velocity at CL top, CL')

end subroutine eddy_diff_init

subroutine eddy_diff_tend(state, pbuf, cam_in, &
     ztodt, p, tint, rhoi, cldn, wstarent, &
     kvm_in, kvh_in, ksrftms, dragblj,tauresx, tauresy, &
     rrho, ustar, pblh, kvm, kvh, kvq, cgh, cgs, tpert, qpert, &
     tke, sprod, sfi, turbtype, sm_aw)

  use physics_types, only: physics_state
  use camsrfexch, only: cam_in_t
  use coords_1d, only: Coords1D

  type(physics_state), intent(in) :: state
  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  type(cam_in_t), intent(in) :: cam_in
  real(r8), intent(in) :: ztodt
  type(Coords1D), intent(in) :: p
  real(r8), intent(in) :: tint(pcols,pver+1)
  real(r8), intent(in) :: rhoi(pcols,pver+1)
  real(r8), intent(in) :: cldn(pcols,pver)
  logical, intent(in) :: wstarent
  real(r8), intent(in) :: kvm_in(pcols,pver+1)
  real(r8), intent(in) :: kvh_in(pcols,pver+1)
  real(r8), intent(in) :: ksrftms(pcols)
  real(r8), intent(in) :: dragblj(pcols,pver)       ! Drag profile from Beljaars SGO form drag [ 1/s ]
  real(r8), intent(inout) :: tauresx(pcols)
  real(r8), intent(inout) :: tauresy(pcols)
  real(r8), intent(out) :: rrho(pcols)
  real(r8), intent(out) :: ustar(pcols)
  real(r8), intent(out) :: pblh(pcols)
  real(r8), intent(out) :: kvm(pcols,pver+1)
  real(r8), intent(out) :: kvh(pcols,pver+1)
  real(r8), intent(out) :: kvq(pcols,pver+1)
  real(r8), intent(out) :: cgh(pcols,pver+1)
  real(r8), intent(out) :: cgs(pcols,pver+1)
  real(r8), intent(out) :: tpert(pcols)
  real(r8), intent(out) :: qpert(pcols)
  real(r8), intent(out) :: tke(pcols,pver+1)
  real(r8), intent(out) :: sprod(pcols,pver+1)
  real(r8), intent(out) :: sfi(pcols,pver+1)
  integer(i4), intent(out) :: turbtype(pcols,pver+1)
  real(r8), intent(out) :: sm_aw(pcols,pver+1)

  integer :: i, k

  call compute_eddy_diff( pbuf, state%lchnk    ,                                     &
       pcols    , pver        , state%ncol       , state%t    , tint, state%q(:,:,1) , ztodt   , &
       state%q(:,:,ixcldliq)  , state%q(:,:,ixcldice)   ,                            &
       state%s  , p           , rhoi, cldn       , &
       state%zm , state%zi    , state%pmid , state%pint , state%u        , state%v , &
       cam_in%wsx, cam_in%wsy , cam_in%shf , cam_in%cflx(:,1) , wstarent           , &
       rrho     , ustar       , pblh       , kvm_in     , kvh_in         , kvm     , &
       kvh      , kvq         , cgh        ,                                         &
       cgs      , tpert       , qpert      , tke            ,                        &
       sprod    , sfi         ,                                                      &
       tauresx  , tauresy     , ksrftms    , dragblj , turbtype   , sm_aw )

  ! The diffusivities from diag_TKE can be much larger than from HB in the free
  ! troposphere and upper atmosphere. These seem to be larger than observations,
  ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
  ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
  ! or the upper atmosphere.
  !
  ! NOTE: Further investigation should be done as to why the diffusivities are
  ! larger in diag_TKE.
  if ((kv_freetrop_scale /= 1._r8) .or. ((kv_top_scale /= 1._r8) .and. (kv_top_pressure > 0._r8))) then
     do i = 1, state%ncol
        do k = 1, pverp
           ! Outside of the boundary layer?
           if (state%zi(i,k) > pblh(i)) then
              ! In the upper atmosphere?
              if (state%pint(i,k) <= kv_top_pressure) then
                 kvh(i,k) = kvh(i,k) * kv_top_scale
                 kvm(i,k) = kvm(i,k) * kv_top_scale
                 kvq(i,k) = kvq(i,k) * kv_top_scale
              else
                 kvh(i,k) = kvh(i,k) * kv_freetrop_scale
                 kvm(i,k) = kvm(i,k) * kv_freetrop_scale
                 kvq(i,k) = kvq(i,k) * kv_freetrop_scale
              end if
           else
              exit
           end if
        end do
     end do
  end if

end subroutine eddy_diff_tend

!=============================================================================== !
!                                                                                !
!=============================================================================== !

subroutine compute_eddy_diff( pbuf, lchnk  ,                                                      &
                              pcols  , pver   , ncol     , t       , tint, qv       , ztodt   ,   &
                              ql     , qi     , s        , p       , rhoi, cldn     ,             &
                              z      , zi     , pmid     , pi      , u        , v       ,         &
                              taux   , tauy   , shflx    , qflx    , wstarent ,           rrho  , &
                              ustar  , pblh   , kvm_in   , kvh_in  , kvm_out  , kvh_out , kvq   , &
                              cgh    , cgs    , tpert    , qpert   , tke     ,                    &
                              sprod  , sfi    ,                                                   &
                              tauresx, tauresy, ksrftms, dragblj, turbtype, sm_aw )

  !-------------------------------------------------------------------- !
  ! Purpose: Interface to compute eddy diffusivities.                   !
  !          Eddy diffusivities are calculated in a fully implicit way  !
  !          through iteration process.                                 !
  ! Author:  Sungsu Park. August. 2006.                                 !
  !                       May.    2008.                                 !
  !-------------------------------------------------------------------- !

  use diffusion_solver, only: compute_vdiff
  use cam_history,      only: outfld
  use phys_debug_util,  only: phys_debug_col
  use physconst,        only: cpairv
  use pbl_utils,        only: calc_ustar, austausch_atm
  use error_messages,   only: handle_errmsg
  use coords_1d,        only: Coords1D
  use wv_saturation,    only: qsat
  use eddy_diff,        only: trbintd, caleddy
  use physics_buffer,   only: pbuf_get_field

  ! --------------- !
  ! Input Variables !
  ! --------------- !

  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  integer,  intent(in)    :: lchnk
  integer,  intent(in)    :: pcols                     ! Number of atmospheric columns [ # ]
  integer,  intent(in)    :: pver                      ! Number of atmospheric layers  [ # ]
  integer,  intent(in)    :: ncol                      ! Number of atmospheric columns [ # ]
  logical,  intent(in)    :: wstarent                  ! .true. means use the 'wstar' entrainment closure.
  real(r8), intent(in)    :: ztodt                     ! Physics integration time step 2 delta-t [ s ]
  real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [ K ]
  real(r8), intent(in)    :: tint(pcols,pver+1)        ! Temperature defined on interfaces [ K ]
  real(r8), intent(in)    :: qv(pcols,pver)            ! Water vapor  specific humidity [ kg/kg ]
  real(r8), intent(in)    :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
  real(r8), intent(in)    :: qi(pcols,pver)            ! Ice specific humidity [ kg/kg ]
  real(r8), intent(in)    :: s(pcols,pver)             ! Dry static energy [ J/kg ]
  type(Coords1D), intent(in) :: p                      ! Pressure coordinates for solver [ Pa ]
  real(r8), intent(in)    :: rhoi(pcols,pver+1)        ! Density at interfaces [ kg/m3 ]
  real(r8), intent(in)    :: cldn(pcols,pver)          ! Stratiform cloud fraction [ fraction ]
  real(r8), intent(in)    :: z(pcols,pver)             ! Layer mid-point height above surface [ m ]
  real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface height above surface [ m ]
  real(r8), intent(in)    :: pmid(pcols,pver)          ! Layer mid-point pressure [ Pa ]
  real(r8), intent(in)    :: pi(pcols,pver+1)          ! Interface pressure [ Pa ]
  real(r8), intent(in)    :: u(pcols,pver)             ! Zonal velocity [ m/s ]
  real(r8), intent(in)    :: v(pcols,pver)             ! Meridional velocity [ m/s ]
  real(r8), intent(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
  real(r8), intent(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]
  real(r8), intent(in)    :: shflx(pcols)              ! Sensible heat flux at surface [ unit ? ]
  real(r8), intent(in)    :: qflx(pcols)               ! Water vapor flux at surface [ unit ? ]
  real(r8), intent(in)    :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep [ m2/s ]
  real(r8), intent(in)    :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep [ m2/s ]
  real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient of turbulent mountain stress [ unit ? ]
  real(r8), intent(in)    :: dragblj(pcols,pver)       ! Drag profile from Beljaars SGO form drag [ 1/s ]

  ! ---------------- !
  ! Output Variables !
  ! ---------------- !

  real(r8), intent(out)   :: kvm_out(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
  real(r8), intent(out)   :: kvh_out(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
  real(r8), intent(out)   :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ]
                                                       ! (note not having '_out')
  real(r8), intent(out)   :: rrho(pcols)               ! Reciprocal of density at the lowest layer
  real(r8), intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
  real(r8), intent(out)   :: pblh(pcols)               ! PBL top height [ m ]
  real(r8), intent(out)   :: cgh(pcols,pver+1)         ! Counter-gradient term for heat [ J/kg/m ]
  real(r8), intent(out)   :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
  real(r8), intent(out)   :: tpert(pcols)              ! Convective temperature excess [ K ]
  real(r8), intent(out)   :: qpert(pcols)              ! Convective humidity excess [ kg/kg ]
  real(r8), intent(out)   :: tke(pcols,pver+1)         ! Turbulent kinetic energy [ m2/s2 ]
  real(r8), intent(out)   :: sprod(pcols,pver+1)       ! Shear production [ m2/s3 ]
  real(r8), intent(out)   :: sfi(pcols,pver+1)         ! Interfacial layer saturation fraction [ fraction ]
  integer(i4), intent(out):: turbtype(pcols,pver+1)    ! Turbulence type identifier at all interfaces [ no unit ]
  real(r8), intent(out)   :: sm_aw(pcols,pver+1)       ! Normalized Galperin instability function for momentum [ no unit ]
                                                       ! This is 1 when neutral condition (Ri=0),
                                                       ! 4.964 for maximum unstable case, and 0 when Ri > Ricrit=0.19.

  ! ---------------------- !
  ! Input-Output Variables !
  ! ---------------------- !

  real(r8), intent(inout) :: tauresx(pcols)            ! Residual stress to be added in vdiff to correct for turb
  real(r8), intent(inout) :: tauresy(pcols)            ! Stress mismatch between sfc and atm accumulated in prior timesteps

  ! -------------- !
  ! pbuf Variables !
  ! -------------- !

  real(r8), pointer :: qrl(:,:)                        ! LW radiative cooling rate
  real(r8), pointer :: wsedl(:,:)                      ! Sedimentation velocity
                                                       ! of stratiform liquid cloud droplet [ m/s ]

  real(r8), pointer :: bprod(:,:)                      ! Buoyancy production of tke [ m2/s3 ]
  integer(i4), pointer :: ipbl(:)                      ! If 1, PBL is CL, while if 0, PBL is STL.
  integer(i4), pointer :: kpblh(:)                     ! Layer index containing PBL top within or at the base interface
  real(r8), pointer :: wstarPBL(:)                     ! Convective velocity within PBL [ m/s ]
  real(r8), pointer :: tkes(:)                         ! TKE at surface interface [ m2/s2 ]
  real(r8), pointer :: went(:)                         ! Entrainment rate at the PBL top interface [ m/s ]

  ! --------------- !
  ! Local Variables !
  ! --------------- !

  integer                    icol
  integer                    i, k, iturb, status

  character(2048)         :: warnstring                ! Warning(s) to print
  character(128)          :: errstring                 ! Error message

  real(r8)                :: kvf(pcols,pver+1)         ! Free atmospheric eddy diffusivity [ m2/s ]
  real(r8)                :: kvm(pcols,pver+1)         ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: kvm_preo(pcols,pver+1)    ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh_preo(pcols,pver+1)    ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: kvm_pre(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh_pre(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: errorPBL(pcols)           ! Error function showing whether PBL produced convergent solution or not.
                                                       ! [ unit ? ]
  real(r8)                :: s2(pcols,pver)            ! Shear squared, defined at interfaces except surface [ s-2 ]
  real(r8)                :: n2(pcols,pver)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
  real(r8)                :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
  real(r8)                :: pblhp(pcols)              ! PBL top pressure [ Pa ]
  real(r8)                :: minpblh(pcols)            ! Minimum PBL height based on surface stress

  real(r8)                :: qt(pcols,pver)            ! Total specific humidity [ kg/kg ]
  real(r8)                :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
  real(r8)                :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
  real(r8)                :: sl(pcols,pver)            ! Liquid water static energy [ J/kg ]
  real(r8)                :: slv(pcols,pver)           ! Liquid water virtual static energy [ J/kg ]
  real(r8)                :: slslope(pcols,pver)       ! Slope of 'sl' in each layer
  real(r8)                :: qtslope(pcols,pver)       ! Slope of 'qt' in each layer
  real(r8)                :: qvfd(pcols,pver)          ! Specific humidity for diffusion [ kg/kg ]
  real(r8)                :: tfd(pcols,pver)           ! Temperature for diffusion [ K ]
  real(r8)                :: slfd(pcols,pver)          ! Liquid static energy [ J/kg ]
  real(r8)                :: qtfd(pcols,pver)          ! Total specific humidity [ kg/kg ]
  real(r8)                :: qlfd(pcols,pver)          ! Liquid water specific humidity for diffusion [ kg/kg ]
  real(r8)                :: ufd(pcols,pver)           ! U-wind for diffusion [ m/s ]
  real(r8)                :: vfd(pcols,pver)           ! V-wind for diffusion [ m/s ]

  ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'

  real(r8)                :: chu(pcols,pver+1)         ! Heat buoyancy coef for dry states, defined at each interface, finally.
  real(r8)                :: chs(pcols,pver+1)         ! Heat buoyancy coef for sat states, defined at each interface, finally.
  real(r8)                :: cmu(pcols,pver+1)         ! Moisture buoyancy coef for dry states,
                                                       ! defined at each interface, finally.
  real(r8)                :: cms(pcols,pver+1)         ! Moisture buoyancy coef for sat states,
                                                       ! defined at each interface, finally.

  real(r8)                :: jnk1d(pcols)
  real(r8)                :: jnk2d(pcols,pver+1)
  real(r8)                :: zero(pcols)
  real(r8)                :: zero2d(pcols,pver+1)
  real(r8)                :: es                     ! Saturation vapor pressure
  real(r8)                :: qs                     ! Saturation specific humidity
  real(r8)                :: ep2, templ, temps

  ! ------------------------------- !
  ! Variables for diagnostic output !
  ! ------------------------------- !

  real(r8)                :: wpert(pcols)              ! Turbulent velocity excess [ m/s ]

  real(r8)                :: kbase_o(pcols,ncvmax)     ! Original external base interface index of CL from 'exacol'
  real(r8)                :: ktop_o(pcols,ncvmax)      ! Original external top  interface index of CL from 'exacol'
  real(r8)                :: ncvfin_o(pcols)           ! Original number of CLs from 'exacol'
  real(r8)                :: kbase_mg(pcols,ncvmax)    ! 'kbase' after extending-merging from 'zisocl'
  real(r8)                :: ktop_mg(pcols,ncvmax)     ! 'ktop' after extending-merging from 'zisocl'
  real(r8)                :: ncvfin_mg(pcols)          ! 'ncvfin' after extending-merging from 'zisocl'
  real(r8)                :: kbase_f(pcols,ncvmax)     ! Final 'kbase' after extending-merging & including SRCL
  real(r8)                :: ktop_f(pcols,ncvmax)      ! Final 'ktop' after extending-merging & including SRCL
  real(r8)                :: ncvfin_f(pcols)           ! Final 'ncvfin' after extending-merging & including SRCL
  real(r8)                :: wet(pcols,ncvmax)         ! Entrainment rate at the CL top  [ m/s ]
  real(r8)                :: web(pcols,ncvmax)         ! Entrainment rate at the CL base [ m/s ].
                                                       ! Set to zero if CL is based at surface.
  real(r8)                :: jtbu(pcols,ncvmax)        ! Buoyancy jump across the CL top  [ m/s2 ]
  real(r8)                :: jbbu(pcols,ncvmax)        ! Buoyancy jump across the CL base [ m/s2 ]
  real(r8)                :: evhc(pcols,ncvmax)        ! Evaporative enhancement factor at the CL top
  real(r8)                :: jt2slv(pcols,ncvmax)      ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
  real(r8)                :: n2ht(pcols,ncvmax)        ! n2 defined at the CL top  interface but using
                                                       ! sfuh(kt)   instead of sfi(kt) [ s-2 ]
  real(r8)                :: n2hb(pcols,ncvmax)        ! n2 defined at the CL base interface but using
                                                       ! sflh(kb-1) instead of sfi(kb) [ s-2 ]
  real(r8)                :: lwp(pcols,ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
  real(r8)                :: opt_depth(pcols,ncvmax)   ! Optical depth of the CL top layer
  real(r8)                :: radinvfrac(pcols,ncvmax)  ! Fraction of radiative cooling confined in the top portion of CL top layer
  real(r8)                :: radf(pcols,ncvmax)        ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
  real(r8)                :: wstar(pcols,ncvmax)       ! Convective velocity in each CL [ m/s ]
  real(r8)                :: wstar3fact(pcols,ncvmax)  ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
  real(r8)                :: ebrk(pcols,ncvmax)        ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
  real(r8)                :: wbrk(pcols,ncvmax)        ! Net mean normalized TKE (W) of CL,
                                                       ! 'ebrk/b1' including entrainment effect [ m2/s2 ]
  real(r8)                :: lbrk(pcols,ncvmax)        ! Energetic internal thickness of CL [m]
  real(r8)                :: ricl(pcols,ncvmax)        ! CL internal mean Richardson number
  real(r8)                :: ghcl(pcols,ncvmax)        ! Half of normalized buoyancy production of CL
  real(r8)                :: shcl(pcols,ncvmax)        ! Galperin instability function of heat-moisture of CL
  real(r8)                :: smcl(pcols,ncvmax)        ! Galperin instability function of mementum of CL
  real(r8)                :: ghi(pcols,pver+1)         ! Half of normalized buoyancy production at all interfaces
  real(r8)                :: shi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
  real(r8)                :: smi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
  real(r8)                :: rii(pcols,pver+1)         ! Interfacial Richardson number defined at all interfaces
  real(r8)                :: lengi(pcols,pver+1)       ! Turbulence length scale at all interfaces [ m ]
  real(r8)                :: wcap(pcols,pver+1)        ! Normalized TKE at all interfaces [ m2/s2 ]
  ! For sedimentation-entrainment feedback
  real(r8)                :: wsed(pcols,ncvmax)        ! Sedimentation velocity at the top of each CL [ m/s ]

  ! ---------- !
  ! Parameters !
  ! ---------- !

  logical,          parameter :: use_kvf        =  .false.      ! .true. (.false.) : initialize kvh/kvm =  kvf ( 0. )
  real(r8),         parameter :: lambda         =   0.5_r8      ! Under-relaxation factor ( 0 < lambda =< 1 )

  ! ---------- !
  ! Initialize !
  ! ---------- !

  zero(:)     = 0._r8
  zero2d(:,:) = 0._r8

  ! ---------------------------------------------- !
  ! Get LW radiative heating out of physics buffer !
  ! ---------------------------------------------- !
  call pbuf_get_field(pbuf, qrl_idx,   qrl)
  call pbuf_get_field(pbuf, wsedl_idx, wsedl)

  ! These fields are put into the pbuf for UNICON only.
  if (unicon_is_on) then
     call pbuf_get_field(pbuf, bprod_idx,    bprod)
     call pbuf_get_field(pbuf, ipbl_idx,     ipbl)
     call pbuf_get_field(pbuf, kpblh_idx,    kpblh)
     call pbuf_get_field(pbuf, wstarPBL_idx, wstarPBL)
     call pbuf_get_field(pbuf, tkes_idx,     tkes)
     call pbuf_get_field(pbuf, went_idx,     went)
  else
     allocate(bprod(pcols,pverp), ipbl(pcols), kpblh(pcols), wstarPBL(pcols), tkes(pcols), went(pcols))
  end if

  ! ----------------------- !
  ! Main Computation Begins !
  ! ----------------------- !

  ufd(:ncol,:)  = u(:ncol,:)
  vfd(:ncol,:)  = v(:ncol,:)
  tfd(:ncol,:)  = t(:ncol,:)
  qvfd(:ncol,:) = qv(:ncol,:)
  qlfd(:ncol,:) = ql(:ncol,:)

  do iturb = 1, nturb

     ! Total stress includes 'tms'.
     ! Here, in computing 'tms', we can use either iteratively changed 'ufd,vfd' or the
     ! initially given 'u,v' to the PBL scheme. Note that normal stress, 'taux, tauy'
     ! are not changed by iteration. In order to treat 'tms' in a fully implicit way,
     ! I am using updated wind, here.

     ! Compute ustar
     call calc_ustar( ncol, tfd(:ncol,pver), pmid(:ncol,pver), &
                      taux(:ncol) - ksrftms(:ncol) * ufd(:ncol,pver), & ! Zonal wind stress
                      tauy(:ncol) - ksrftms(:ncol) * vfd(:ncol,pver), & ! Meridional wind stress
                      rrho(:ncol), ustar(:ncol))
     minpblh(:ncol) = 100.0_r8 * ustar(:ncol)   ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'.

     ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)

     call trbintd( &
                   pcols    , pver    , ncol  , z       , ufd     , vfd     , tfd   , pmid    , &
                   s2       , n2      , ri    , zi      , pi      , cldn    , qtfd  , qvfd    , &
                   qlfd     , qi      , sfi   , sfuh    , sflh    , slfd    , slv   , slslope , &
                   qtslope  , chs     , chu   , cms     , cmu     )

     ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.
     ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables.

     if( iturb == 1 ) then
        qt(:ncol,:) = qtfd(:ncol,:)
        sl(:ncol,:) = slfd(:ncol,:)
     endif

     ! Get free atmosphere exchange coefficients. This 'kvf' is not used in UW moist PBL scheme
     if (use_kvf) then
        call austausch_atm(pcols, ncol, pver, ntop_eddy, nbot_eddy, &
             ml2, ri, s2, kvf )
     else
        kvf = 0._r8
     end if

     ! Initialize kvh/kvm to send to caleddy, depending on model timestep and iteration number
     ! This is necessary for 'wstar-based' entrainment closure.

     if( iturb == 1 ) then
        if( is_first_step() ) then
           ! First iteration of first model timestep : Use free tropospheric value or zero.
           kvh(:ncol,:) = kvf(:ncol,:)
           kvm(:ncol,:) = kvf(:ncol,:)
        else
           ! First iteration on any model timestep except the first : Use value from previous timestep
           kvh(:ncol,:) = kvh_in(:ncol,:)
           kvm(:ncol,:) = kvm_in(:ncol,:)
        endif
     else
        ! Not the first iteration : Use from previous iteration
        kvh(:ncol,:) = kvh_out(:ncol,:)
        kvm(:ncol,:) = kvm_out(:ncol,:)
     endif

     ! Calculate eddy diffusivity (kvh_out,kvm_out) and (tke,bprod,sprod) using
     ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
     ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
     ! caleddy after calculating (kvh_out,kvm_out)

     call caleddy( pcols     , pver      , ncol      ,                     &
                   slfd      , qtfd      , qlfd      , slv      ,ufd     , &
                   vfd       , pi        , z         , zi       ,          &
                   qflx      , shflx     , slslope   , qtslope  ,          &
                   chu       , chs       , cmu       , cms      ,sfuh    , &
                   sflh      , n2        , s2        , ri       ,rrho    , &
                   pblh      , ustar     ,                                 &
                   kvh       , kvm       , kvh_out   , kvm_out  ,          &
                   tpert     , qpert     , qrl       , kvf      , tke    , &
                   wstarent  , bprod     , sprod     , minpblh  , wpert  , &
                   tkes      , went      , turbtype  , sm_aw    ,          &
                   kbase_o   , ktop_o    , ncvfin_o  ,                     &
                   kbase_mg  , ktop_mg   , ncvfin_mg ,                     &
                   kbase_f   , ktop_f    , ncvfin_f  ,                     &
                   wet       , web       , jtbu      , jbbu     ,          &
                   evhc      , jt2slv    , n2ht      , n2hb     ,          &
                   lwp       , opt_depth , radinvfrac, radf     ,          &
                   wstar     , wstar3fact,                                 &
                   ebrk      , wbrk      , lbrk      , ricl     , ghcl   , &
                   shcl      , smcl      , ghi       , shi      , smi    , &
                   rii       , lengi     , wcap      , pblhp    , cldn   , &
                   ipbl      , kpblh     , wsedl     , wsed, &
                   warnstring, errstring)

     if (trim(warnstring) /= "") then
        write(iulog,*) "eddy_diff_cam: Messages from caleddy follow."
        write(iulog,*) warnstring
     end if

     call handle_errmsg(errstring, subname="caleddy")

     ! Calculate errorPBL to check whether PBL produced convergent solutions or not.

     if( iturb == nturb ) then
        do i = 1, ncol
           errorPBL(i) = 0._r8
           do k = 1, pver
              errorPBL(i) = errorPBL(i) + ( kvh(i,k) - kvh_out(i,k) )**2
           end do
           errorPBL(i) = sqrt(errorPBL(i)/pver)
        end do
     end if

     ! Eddy diffusivities which will be used for the initialization of (bprod,
     ! sprod) in 'caleddy' at the next iteration step.

     if( iturb > 1 .and. iturb < nturb ) then
        kvm_out(:ncol,:) = lambda * kvm_out(:ncol,:) + ( 1._r8 - lambda ) * kvm(:ncol,:)
        kvh_out(:ncol,:) = lambda * kvh_out(:ncol,:) + ( 1._r8 - lambda ) * kvh(:ncol,:)
     endif

     ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.

     cgh(:ncol,:) = 0._r8
     cgs(:ncol,:) = 0._r8

     if( iturb < nturb ) then

        ! Each time we diffuse the original state

        slfd(:ncol,:)  = sl(:ncol,:)
        qtfd(:ncol,:)  = qt(:ncol,:)
        ufd(:ncol,:)   = u(:ncol,:)
        vfd(:ncol,:)   = v(:ncol,:)

        ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
        ! In the below 'compute_vdiff', (slfd,qtfd,ufd,vfd) are 'inout' variables.

        call compute_vdiff( lchnk   ,                                                  &
                            pcols   , pver     , 1        , ncol         , tint, &
                            p        , t        , rhoi, ztodt        , taux      , &
                            tauy    , shflx    , qflx     , &
                            kvh_out , kvm_out  , kvh_out  , cgs          , cgh       , &
                            zi      , ksrftms  , dragblj  , & 
                            zero    , fieldlist_wet, fieldlist_molec, &
                            ufd     , vfd      , qtfd     , slfd         ,             &
                            jnk1d   , jnk1d    , jnk2d    , jnk1d        , errstring , &
                            tauresx , tauresy  , 0        , cpairv(:,:,lchnk), zero, &
                            .false., .false. )

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="compute_vdiff called from eddy_diff_cam")

        ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to
        ! use 'trbintd' at the next iteration.

        do k = 1, pver
           do i = 1, ncol
              ! ----------------------------------------------------- !
              ! Compute the condensate 'qlfd' in the updated profiles !
              ! ----------------------------------------------------- !
              ! Option.1 : Assume grid-mean condensate is homogeneously diffused by the moist turbulence scheme.
              !            This should be used if 'pseudodiff = .false.' in vertical_diffusion.F90.
              ! Modification : Need to be check whether below is correct in the presence of ice, qi.
              !                I should understand why the variation of ice, qi is neglected during diffusion.
              templ     = ( slfd(i,k) - gravit*z(i,k) ) / cpair
              call qsat( templ, pmid(i,k), es, qs)
              ep2       =  .622_r8
              temps     =   templ + ( qtfd(i,k) - qs ) / ( cpair / latvap + latvap * qs / ( rair * templ**2 ) )
              call qsat( temps, pmid(i,k), es, qs)
              qlfd(i,k) =   max( qtfd(i,k) - qi(i,k) - qs ,0._r8 )
              ! Option.2 : Assume condensate is not diffused by the moist turbulence scheme.
              !            This should bs used if 'pseudodiff = .true.'  in vertical_diffusion.F90.
              ! qlfd(i,k) = ql(i,k)
              ! ----------------------------- !
              ! Compute the other 'qvfd, tfd' !
              ! ----------------------------- !
              qvfd(i,k) = max( 0._r8, qtfd(i,k) - qi(i,k) - qlfd(i,k) )
              tfd(i,k)  = ( slfd(i,k) + latvap * qlfd(i,k) + (latvap+latice) * qi(i,k) - gravit*z(i,k)) / cpair
           end do
        end do
     endif

  end do  ! End of 'iturb' iteration

  kvq(:ncol,:) = kvh_out(:ncol,:)

  ! Compute 'wstar' within the PBL for use in the future convection scheme.

  do i = 1, ncol
     if(ipbl(i) == 1) then
        wstarPBL(i) = max( 0._r8, wstar(i,1) )
     else
        wstarPBL(i) = 0._r8
     endif
  end do

  ! --------------------------------------------------------------- !
  ! Writing for detailed diagnostic analysis of UW moist PBL scheme !
  ! --------------------------------------------------------------- !

  call outfld( 'WGUSTD' , wpert, pcols, lchnk )

  call outfld( 'BPROD   ', bprod, pcols, lchnk )
  call outfld( 'SPROD   ', sprod, pcols, lchnk )
  call outfld( 'SFI     ', sfi,   pcols, lchnk )

  call outfld( 'UW_errorPBL',    errorPBL,   pcols,   lchnk )

  call outfld( 'UW_n2',          n2,         pcols,   lchnk )
  call outfld( 'UW_s2',          s2,         pcols,   lchnk )
  call outfld( 'UW_ri',          ri,         pcols,   lchnk )

  call outfld( 'UW_sfuh',        sfuh,       pcols,   lchnk )
  call outfld( 'UW_sflh',        sflh,       pcols,   lchnk )
  call outfld( 'UW_sfi',         sfi,        pcols,   lchnk )

  call outfld( 'UW_cldn',        cldn,       pcols,   lchnk )
  call outfld( 'UW_qrl',         qrl,        pcols,   lchnk )
  call outfld( 'UW_ql',          qlfd,       pcols,   lchnk )

  call outfld( 'UW_chu',         chu,        pcols,   lchnk )
  call outfld( 'UW_chs',         chs,        pcols,   lchnk )
  call outfld( 'UW_cmu',         cmu,        pcols,   lchnk )
  call outfld( 'UW_cms',         cms,        pcols,   lchnk )

  call outfld( 'UW_tke',         tke,        pcols,   lchnk )
  call outfld( 'UW_wcap',        wcap,       pcols,   lchnk )
  call outfld( 'UW_bprod',       bprod,      pcols,   lchnk )
  call outfld( 'UW_sprod',       sprod,      pcols,   lchnk )

  call outfld( 'UW_kvh',         kvh_out,    pcols,   lchnk )
  call outfld( 'UW_kvm',         kvm_out,    pcols,   lchnk )

  call outfld( 'UW_pblh',        pblh,       pcols,   lchnk )
  call outfld( 'UW_pblhp',       pblhp,      pcols,   lchnk )
  call outfld( 'UW_tpert',       tpert,      pcols,   lchnk )
  call outfld( 'UW_qpert',       qpert,      pcols,   lchnk )
  call outfld( 'UW_wpert',       wpert,      pcols,   lchnk )

  call outfld( 'UW_ustar',       ustar,      pcols,   lchnk )
  call outfld( 'UW_tkes',        tkes,       pcols,   lchnk )
  call outfld( 'UW_minpblh',     minpblh,    pcols,   lchnk )

  call outfld( 'UW_turbtype',    real(turbtype,r8),   pcols,   lchnk )

  call outfld( 'UW_kbase_o',     kbase_o,    pcols,   lchnk )
  call outfld( 'UW_ktop_o',      ktop_o,     pcols,   lchnk )
  call outfld( 'UW_ncvfin_o',    ncvfin_o,   pcols,   lchnk )

  call outfld( 'UW_kbase_mg',    kbase_mg,   pcols,   lchnk )
  call outfld( 'UW_ktop_mg',     ktop_mg,    pcols,   lchnk )
  call outfld( 'UW_ncvfin_mg',   ncvfin_mg,  pcols,   lchnk )

  call outfld( 'UW_kbase_f',     kbase_f,    pcols,   lchnk )
  call outfld( 'UW_ktop_f',      ktop_f,     pcols,   lchnk )
  call outfld( 'UW_ncvfin_f',    ncvfin_f,   pcols,   lchnk )

  call outfld( 'UW_wet',         wet,        pcols,   lchnk )
  call outfld( 'UW_web',         web,        pcols,   lchnk )
  call outfld( 'UW_jtbu',        jtbu,       pcols,   lchnk )
  call outfld( 'UW_jbbu',        jbbu,       pcols,   lchnk )
  call outfld( 'UW_evhc',        evhc,       pcols,   lchnk )
  call outfld( 'UW_jt2slv',      jt2slv,     pcols,   lchnk )
  call outfld( 'UW_n2ht',        n2ht,       pcols,   lchnk )
  call outfld( 'UW_n2hb',        n2hb,       pcols,   lchnk )
  call outfld( 'UW_lwp',         lwp,        pcols,   lchnk )
  call outfld( 'UW_optdepth',    opt_depth,  pcols,   lchnk )
  call outfld( 'UW_radfrac',     radinvfrac, pcols,   lchnk )
  call outfld( 'UW_radf',        radf,       pcols,   lchnk )
  call outfld( 'UW_wstar',       wstar,      pcols,   lchnk )
  call outfld( 'UW_wstar3fact',  wstar3fact, pcols,   lchnk )
  call outfld( 'UW_ebrk',        ebrk,       pcols,   lchnk )
  call outfld( 'UW_wbrk',        wbrk,       pcols,   lchnk )
  call outfld( 'UW_lbrk',        lbrk,       pcols,   lchnk )
  call outfld( 'UW_ricl',        ricl,       pcols,   lchnk )
  call outfld( 'UW_ghcl',        ghcl,       pcols,   lchnk )
  call outfld( 'UW_shcl',        shcl,       pcols,   lchnk )
  call outfld( 'UW_smcl',        smcl,       pcols,   lchnk )

  call outfld( 'UW_gh',          ghi,        pcols,   lchnk )
  call outfld( 'UW_sh',          shi,        pcols,   lchnk )
  call outfld( 'UW_sm',          smi,        pcols,   lchnk )
  call outfld( 'UW_ria',         rii,        pcols,   lchnk )
  call outfld( 'UW_leng',        lengi,      pcols,   lchnk )

  call outfld( 'UW_wsed',        wsed,       pcols,   lchnk )

  if (.not. unicon_is_on) then
     deallocate(bprod, ipbl, kpblh, wstarPBL, tkes, went)
  end if

end subroutine compute_eddy_diff

end module eddy_diff_cam
