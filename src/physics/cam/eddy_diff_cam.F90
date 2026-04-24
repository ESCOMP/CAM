module eddy_diff_cam

use shr_kind_mod, only: i4 => shr_kind_i4, r8 => shr_kind_r8
use ppgrid, only: pcols, pver, pverp
use cam_logfile, only: iulog
use cam_abortutils, only: endrun
use physconst, only: gravit, cpair, rair, zvir, latvap, latice, karman
use time_manager, only: is_first_step
use physics_buffer, only: physics_buffer_desc
use spmd_utils, only: masterproc
use phys_control, only: phys_getopts

implicit none
private

public :: eddy_diff_readnl
public :: eddy_diff_init
public :: eddy_diff_tend

! Cloud mass constituent indices
integer :: ixcldliq, ixcldice

! input pbuf field indices
integer :: qrl_idx   = -1
integer :: wsedl_idx = -1

integer :: ncvmax

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

subroutine eddy_diff_init(ntop_eddy_in)

  use error_messages, only: handle_errmsg
  use cam_history, only: addfld, add_default, horiz_only
  use constituents, only: cnst_get_ind
  use ref_pres, only: pref_mid
  use physics_buffer, only: pbuf_get_index

  use bretherton_park_diff, only: bretherton_park_diff_init

  integer,  intent(in) :: ntop_eddy_in ! Top interface level to which eddy vertical diffusivity is applied ( = 1 )

  character(len=512)   :: errmsg
  integer              :: errflg

  logical :: history_amwg

  ! Call CCPPized subroutine:
  call bretherton_park_diff_init(masterproc, iulog, pver, &
   gravit, cpair, rair, zvir, latvap, latice, karman, &
   ntop_eddy_in, &
   pref_mid, &
   eddy_lbulk_max, eddy_leng_max, eddy_max_bot_pressure, eddy_moist_entrain_a2l, &
   ncvmax, errmsg, errflg)
  if(errflg /= 0) then
   call endrun('bretherton_park_diff_init: ' // errmsg)
  endif

  ! Cloud mass constituents
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)

  ! Input pbuf fields
  qrl_idx   = pbuf_get_index('QRL')
  wsedl_idx = pbuf_get_index('WSEDL')

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
     ztodt, do_iss, fv_am_correction, &
     p, tint, rhoi, dpidz_sq, cldn, wstarent, &
     kvm_in, kvh_in, ksrftms, dragblj,tauresx, tauresy, &
     rrho, ustar, pblh, kvm, kvh, kvq, cgh, cgs, tpert, qpert, &
     tke, sprod, sfi)

  use physics_types, only: physics_state
  use camsrfexch, only: cam_in_t
  use coords_1d, only: Coords1D
  use physics_buffer,       only: pbuf_get_field
  use cam_history,          only: outfld

  use constituents,              only: pcnst
  use ccpp_constituent_prop_mod, only: ccpp_const_props
  use beljaars_drag_cam,         only: do_beljaars

  ! CCPPized subroutines
  use bretherton_park_diff,      only: bretherton_park_diff_run
  use eddy_diffusivity_adjustment_above_pbl, only: eddy_diffusivity_adjustment_above_pbl_run

  type(physics_state), intent(in) :: state
  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  type(cam_in_t), intent(in) :: cam_in
  real(r8), intent(in) :: ztodt
  logical,         intent(in)       :: do_iss           ! Use implicit turbulent surface stress computation
  logical,         intent(in)       :: fv_am_correction    ! Do angular momentum conservation correction
  type(Coords1D), intent(in) :: p
  real(r8), intent(in) :: tint(pcols,pverp)
  real(r8), intent(in) :: rhoi(pcols,pverp)
  real(r8), intent(in) :: dpidz_sq(pcols,pverp)
  real(r8), intent(in) :: cldn(pcols,pver)
  logical,  intent(in) :: wstarent
  real(r8), intent(in) :: kvm_in(pcols,pverp)
  real(r8), intent(in) :: kvh_in(pcols,pverp)
  real(r8), intent(in) :: ksrftms(pcols)
  real(r8), intent(in) :: dragblj(pcols,pver)       ! Drag profile from Beljaars SGO form drag [ 1/s ]
  real(r8), intent(inout) :: tauresx(pcols)
  real(r8), intent(inout) :: tauresy(pcols)
  real(r8), intent(out) :: rrho(pcols)
  real(r8), intent(out) :: ustar(pcols)
  real(r8), intent(out) :: pblh(pcols)
  real(r8), intent(out) :: kvm(pcols,pverp)
  real(r8), intent(out) :: kvh(pcols,pverp)
  real(r8), intent(out) :: kvq(pcols,pverp)
  real(r8), intent(out) :: cgh(pcols,pverp)
  real(r8), intent(out) :: cgs(pcols,pverp)
  real(r8), intent(out) :: tpert(pcols)
  real(r8), intent(out) :: qpert(pcols)
  real(r8), intent(out) :: tke(pcols,pverp)
  real(r8), intent(out) :: sprod(pcols,pverp)
  real(r8), intent(out) :: sfi(pcols,pverp)

  ! pbuf fields
  real(r8), pointer :: qrl(:,:)                        ! LW radiative cooling rate [J Pa kg-1 s-1]
  real(r8), pointer :: wsedl(:,:)                      ! Sedimentation velocity of stratiform liquid cloud droplet [m s-1]

  integer :: i, k
  integer :: ncol, lchnk

  ! outputs from UW PBL scheme for history output
  real(r8) :: bprod(pcols,pverp)
  real(r8) :: s2(pcols,pver)            ! Shear squared, defined at interfaces except surface [ s-2 ]
  real(r8) :: n2(pcols,pver)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
  real(r8) :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
  real(r8) :: wpert(pcols)              ! Turbulent velocity excess [m s-1]
  real(r8) :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
  real(r8) :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
  real(r8) :: qlfd(pcols,pver)          ! Liquid water specific humidity for diffusion [ kg/kg ]
  ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'
  real(r8) :: chu(pcols,pverp)          ! Heat buoyancy coef for dry states, interfaces
  real(r8) :: chs(pcols,pverp)          ! Heat buoyancy coef for sat states, interfaces
  real(r8) :: cmu(pcols,pverp)          ! Moisture buoyancy coef for dry states, interfaces
  real(r8) :: cms(pcols,pverp)          ! Moisture buoyancy coef for sat states, interfaces
  real(r8) :: errorPBL(pcols)           ! Error function showing whether PBL produced convergent solution or not [m2 s-1]
  real(r8) :: pblhp(pcols)              ! PBL top pressure [Pa]
  real(r8) :: minpblh(pcols)            ! Minimum PBL height based on surface stress [m]
  real(r8) :: tkes(pcols)               ! Specific TKE at surface interface [ m2/s2 ]
  real(r8) :: wcap(pcols,pver+1)        ! Normalized TKE at all interfaces [ m2/s2 ]
  integer  :: turbtype(pcols,pverp)     ! Turbulence type identifier at all interfaces [ no unit ]
  real(r8) :: kbase_o(pcols,ncvmax)     ! Original external base interface index of CL from 'exacol'
  real(r8) :: ktop_o(pcols,ncvmax)      ! Original external top  interface index of CL from 'exacol'
  real(r8) :: ncvfin_o(pcols)           ! Original number of CLs from 'exacol'
  real(r8) :: kbase_mg(pcols,ncvmax)    ! 'kbase' after extending-merging from 'zisocl'
  real(r8) :: ktop_mg(pcols,ncvmax)     ! 'ktop' after extending-merging from 'zisocl'
  real(r8) :: ncvfin_mg(pcols)          ! 'ncvfin' after extending-merging from 'zisocl'
  real(r8) :: kbase_f(pcols,ncvmax)     ! Original external base interface index of CL from 'exacol'
  real(r8) :: ktop_f(pcols,ncvmax)      ! Original external top  interface index of CL from 'exacol'
  real(r8) :: ncvfin_f(pcols)           ! Original number of CLs from 'exacol'
  real(r8) :: wet(pcols,ncvmax)         ! Entrainment rate at the CL top, ncvmax  [m s-1]
  real(r8) :: web(pcols,ncvmax)         ! Entrainment rate at the CL base, ncvmax [m s-1] (Set to zero if CL is based at surface)
  real(r8) :: jtbu(pcols,ncvmax)        ! Buoyancy jump across the CL top, ncvmax  [m s-2]
  real(r8) :: jbbu(pcols,ncvmax)        ! Buoyancy jump across the CL base, ncvmax [m s-2]
  real(r8) :: evhc(pcols,ncvmax)        ! Evaporative enhancement factor at the CL top, ncvmax
  real(r8) :: jt2slv(pcols,ncvmax)      ! Jump of slv (liquid water virtual static energy) (across two layers)
                                        ! at CL top used only for evhc (evaporative enhancement factor at CL top), ncvmax [J kg-1]
  real(r8) :: n2ht(pcols,ncvmax)        ! n2 defined at the CL top  interface but using
                                        ! sfuh(kt)   instead of sfi(kt), ncvmax [s-2]
  real(r8) :: n2hb(pcols,ncvmax)        ! n2 defined at the CL base interface but using
                                        ! sflh(kb-1) instead of sfi(kb), ncvmax [s-2]
  real(r8) :: lwp(pcols,ncvmax)         ! LWP in the CL top layer, ncvmax [kg m-2]
  real(r8) :: opt_depth(pcols,ncvmax)   ! Optical depth of the CL top layer, ncvmax [1]
  real(r8) :: radinvfrac(pcols,ncvmax)  ! Fraction of radiative cooling confined in the top portion of CL top layer, ncvmax [fraction]
  real(r8) :: radf(pcols,ncvmax)        ! Buoyancy production at the CL top due to LW radiative cooling, ncvmax [m2 s-3]
  real(r8) :: wstar(pcols,ncvmax)       ! Convective velocity in each CL, ncvmax [m s-1]
  real(r8) :: wstar3fact(pcols,ncvmax)  ! Enhancement of 'wstar3' due to entrainment (inverse), ncvmax [1]
  real(r8) :: ebrk(pcols,ncvmax)        ! Net mean TKE of CL including entrainment effect, ncvmax [m2 s-2]
  real(r8) :: wbrk(pcols,ncvmax)        ! Net mean normalized TKE (W) of CL,
                                        ! 'ebrk/b1' including entrainment effect, ncvmax [m2 s-2]
  real(r8) :: lbrk(pcols,ncvmax)        ! Energetic internal thickness of CL, ncvmax [m]
  real(r8) :: ricl(pcols,ncvmax)        ! CL internal mean Richardson number, ncvmax [1]
  real(r8) :: ghcl(pcols,ncvmax)        ! Half of normalized buoyancy production of CL, ncvmax [1]
  real(r8) :: shcl(pcols,ncvmax)        ! Galperin instability function of heat-moisture of CL, ncvmax [1]
  real(r8) :: smcl(pcols,ncvmax)        ! Galperin instability function of mementum of CL, ncvmax [1]
  real(r8) :: ghi(pcols,pverp)          ! Half of normalized buoyancy production at all interfaces [1]
  real(r8) :: shi(pcols,pverp)          ! Galperin instability function of heat-moisture at all interfaces [1]
  real(r8) :: smi(pcols,pverp)          ! Galperin instability function of heat-moisture at all interfaces [1]
  real(r8) :: rii(pcols,pverp)          ! Interfacial Richardson number defined at all interfaces [1]
  real(r8) :: lengi(pcols,pverp)        ! Turbulence length scale at all interfaces [m]
  ! For sedimentation-entrainment feedback
  real(r8) :: wsed(pcols,ncvmax)        ! Sedimentation velocity at the top of each CL [ m/s ]

  character(len=512)   :: errmsg
  integer              :: errflg

  ncol = state%ncol
  lchnk = state%lchnk

  ! ---------------------------------------------- !
  ! Get LW radiative heating out of physics buffer !
  ! ---------------------------------------------- !
  call pbuf_get_field(pbuf, qrl_idx,   qrl)
  call pbuf_get_field(pbuf, wsedl_idx, wsedl)

  ! Update input values to run phase with values from previous timestep (pbuf)
  ! the pbuf field is not passed as inout directly here. This is because
  ! (from the original vertical_diffusion_tend comments:)
  !
  ! kvh (in pbuf) is used by other physics parameterizations,
  ! and as an initial guess in compute_eddy_diff on the next timestep.
  ! It is not updated by the diffusion solver call.
  !
  ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
  ! The contributions for molecular diffusion made to kvm by the call
  ! to the diffusion solver below are not included in the pbuf
  ! as these are not needed in the initial guess by compute_eddy_diff.
  !
  ! There is a pbuf_set_field call after the PBL scheme calls that updates
  ! kvm and kvh in pbuf from the pbuf fields.
  ! The entirety of vertical_diffusion_tend will be obsolete in CAM-SIMA,
  ! and thus the original logic is retained here without further refactoring.
  kvm(:ncol, :pverp) = kvm_in(:ncol, :pverp)
  kvh(:ncol, :pverp) = kvh_in(:ncol, :pverp)

  ! zero out output arrays to pcols
  s2         = 0._r8
  n2         = 0._r8
  ri         = 0._r8
  kvq        = 0._r8
  rrho       = 0._r8
  ustar      = 0._r8
  pblh       = 0._r8
  pblhp      = 0._r8
  minpblh    = 0._r8
  cgh        = 0._r8
  cgs        = 0._r8
  tpert      = 0._r8
  qpert      = 0._r8
  wpert      = 0._r8
  tke        = 0._r8
  tkes       = 0._r8
  wcap       = 0._r8
  wsed       = 0._r8
  turbtype   = 0._r8
  bprod      = 0._r8
  sprod      = 0._r8
  sfi        = 0._r8
  sfuh       = 0._r8
  sflh       = 0._r8
  qlfd       = 0._r8
  chu        = 0._r8
  chs        = 0._r8
  cmu        = 0._r8
  cms        = 0._r8
  kbase_o    = 0._r8
  ktop_o     = 0._r8
  ncvfin_o   = 0._r8
  kbase_mg   = 0._r8
  ktop_mg    = 0._r8
  ncvfin_mg  = 0._r8
  kbase_f    = 0._r8
  ktop_f     = 0._r8
  ncvfin_f   = 0._r8
  wet        = 0._r8
  web        = 0._r8
  jtbu       = 0._r8
  jbbu       = 0._r8
  evhc       = 0._r8
  jt2slv     = 0._r8
  n2ht       = 0._r8
  n2hb       = 0._r8
  lwp        = 0._r8
  opt_depth  = 0._r8
  radinvfrac = 0._r8
  radf       = 0._r8
  wstar      = 0._r8
  wstar3fact = 0._r8
  ebrk       = 0._r8
  wbrk       = 0._r8
  lbrk       = 0._r8
  ricl       = 0._r8
  ghcl       = 0._r8
  shcl       = 0._r8
  smcl       = 0._r8
  ghi        = 0._r8
  shi        = 0._r8
  smi        = 0._r8
  rii        = 0._r8
  lengi      = 0._r8
  errorPBL   = 0._r8

  ! TODO reorder arguments of the subroutine such that in, inout, out (in this order)
  ! Call CCPPized run phase subroutine
  call bretherton_park_diff_run( &
       ncol            = ncol,                          &
       pver            = pver,                          &
       pverp           = pverp,                         &
       pcnst           = pcnst,                         &
       ncvmax          = ncvmax,                        & ! max # of CLs.
       iulog           = iulog,                         &
       dt              = ztodt,                         &
       const_props     = ccpp_const_props,              &
       do_iss          = do_iss,                        &
       am_correction   = fv_am_correction,              &
       do_beljaars     = do_beljaars,                   &
       is_first_timestep= is_first_step(),              &
       gravit          = gravit,                        &
       cpair           = cpair,                         &
       rair            = rair,                          &
       latvap          = latvap,                        &
       latice          = latice,                        &
       t               = state%t(:ncol,:pver),          &
       tint            = tint(:ncol,:pverp),            &
       qv              = state%q(:ncol,:pver,1),        & ! assumes q_wv at 1
       ql              = state%q(:ncol,:pver,ixcldliq), &
       qi              = state%q(:ncol,:pver,ixcldice), &
       s               = state%s(:ncol,:pver),          &
       p               = p,                             &
       rhoi            = rhoi(:ncol,:pverp),            &
       dpidz_sq        = dpidz_sq(:ncol,:pverp),        &
       cldn            = cldn(:ncol,:pver),             &
       z               = state%zm(:ncol,:pver),         &
       zi              = state%zi(:ncol,:pverp),        &
       pmid            = state%pmid(:ncol,:pver),       &
       pint            = state%pint(:ncol,:pverp),      &
       u               = state%u(:ncol,:pver),          &
       v               = state%v(:ncol,:pver),          &
       taux            = cam_in%wsx(:ncol),             &
       tauy            = cam_in%wsy(:ncol),             &
       shflx           = cam_in%shf(:ncol),             &
       qflx            = cam_in%cflx(:ncol,:pcnst),     & ! will be subsetted to wv in run phase.
       wstarent        = wstarent,                      & ! use wstar entrainment? logical
       ksrftms         = ksrftms(:ncol),                &
       dragblj         = dragblj(:ncol,:pver),          &
       qrl             = qrl(:ncol,:pver),              &
       wsedl           = wsedl(:ncol,:pver),            &
       ! below input/output
       tauresx         = tauresx(:ncol),                &
       tauresy         = tauresy(:ncol),                &
       kvm             = kvm(:ncol,:pverp),             & ! in from prev timestep, out from curr timestep.
       kvh             = kvh(:ncol,:pverp),             & ! in from prev timestep, out from curr timestep.
       ! below output
       s2              = s2(:ncol,:pver),               &
       n2              = n2(:ncol,:pver),               &
       ri              = ri(:ncol,:pver),               &
       kvq             = kvq(:ncol,:pverp),             &
       rrho            = rrho(:ncol),                   &
       ustar           = ustar(:ncol),                  &
       pblh            = pblh(:ncol),                   &
       pblhp           = pblhp(:ncol),                  &
       minpblh         = minpblh(:ncol),                &
       cgh             = cgh(:ncol,:pverp),             &
       cgs             = cgs(:ncol,:pverp),             &
       tpert           = tpert(:ncol),                  &
       qpert           = qpert(:ncol),                  &
       wpert           = wpert(:ncol),                  &
       tke             = tke(:ncol,:pverp),             &
       tkes            = tkes(:ncol),                   &
       wcap            = wcap(:ncol,:pverp),            &
       wsed            = wsed(:ncol,:ncvmax),           & ! ncvmax = pver.
       turbtype        = turbtype(:ncol,:pverp),        &
       bprod           = bprod(:ncol,:pverp),           &
       sprod           = sprod(:ncol,:pverp),           &
       sfi             = sfi(:ncol,:pverp),             &
       sfuh            = sfuh(:ncol,:pver),             &
       sflh            = sflh(:ncol,:pver),             &
       qlfd            = qlfd(:ncol,:pver),             &
       chu             = chu(:ncol,:pverp),             &
       chs             = chs(:ncol,:pverp),             &
       cmu             = cmu(:ncol,:pverp),             &
       cms             = cms(:ncol,:pverp),             &
       kbase_o         = kbase_o(:ncol,:ncvmax),        &
       ktop_o          = ktop_o(:ncol,:ncvmax),         &
       ncvfin_o        = ncvfin_o(:ncol),               &
       kbase_mg        = kbase_mg(:ncol,:ncvmax),       &
       ktop_mg         = ktop_mg(:ncol,:ncvmax),        &
       ncvfin_mg       = ncvfin_mg(:ncol),              &
       kbase_f         = kbase_f(:ncol,:ncvmax),        &
       ktop_f          = ktop_f(:ncol,:ncvmax),         &
       ncvfin_f        = ncvfin_f(:ncol),               &
       wet             = wet(:ncol,:ncvmax),            &
       web             = web(:ncol,:ncvmax),            &
       jtbu            = jtbu(:ncol,:ncvmax),           &
       jbbu            = jbbu(:ncol,:ncvmax),           &
       evhc            = evhc(:ncol,:ncvmax),           &
       jt2slv          = jt2slv(:ncol,:ncvmax),         &
       n2ht            = n2ht(:ncol,:ncvmax),           &
       n2hb            = n2hb(:ncol,:ncvmax),           &
       lwp             = lwp(:ncol,:ncvmax),            &
       opt_depth       = opt_depth(:ncol,:ncvmax),      &
       radinvfrac      = radinvfrac(:ncol,:ncvmax),     &
       radf            = radf(:ncol,:ncvmax),           &
       wstar           = wstar(:ncol,:ncvmax),          &
       wstar3fact      = wstar3fact(:ncol,:ncvmax),     &
       ebrk            = ebrk(:ncol,:ncvmax),           &
       wbrk            = wbrk(:ncol,:ncvmax),           &
       lbrk            = lbrk(:ncol,:ncvmax),           &
       ricl            = ricl(:ncol,:ncvmax),           &
       ghcl            = ghcl(:ncol,:ncvmax),           &
       shcl            = shcl(:ncol,:ncvmax),           &
       smcl            = smcl(:ncol,:ncvmax),           &
       ghi             = ghi(:ncol,:pverp),             &
       shi             = shi(:ncol,:pverp),             &
       smi             = smi(:ncol,:pverp),             &
       rii             = rii(:ncol,:pverp),             &
       lengi           = lengi(:ncol,:pverp),           &
       errorPBL        = errorPBL(:ncol),               &
       errmsg          = errmsg,                        &
       errflg          = errflg)

  if(errflg /= 0) then
    call endrun('compute_eddy_diff: ' // errmsg)
  end if

  ! inputs into UW written out as debug:
  call outfld( 'UW_cldn',        cldn,       pcols,   lchnk )
  call outfld( 'UW_qrl',         qrl,        pcols,   lchnk )

  ! outputs from UW:
  call outfld( 'UW_errorPBL',    errorPBL,   pcols,   lchnk )

  call outfld( 'BPROD   ', bprod, pcols, lchnk )
  call outfld( 'UW_bprod',       bprod,      pcols,   lchnk )
  call outfld( 'SPROD   ', sprod, pcols, lchnk )
  call outfld( 'UW_sprod',       sprod,      pcols,   lchnk )

  call outfld( 'WGUSTD' , wpert, pcols, lchnk )
  call outfld( 'UW_wpert',       wpert,      pcols,   lchnk )

  call outfld( 'SFI     ', sfi,   pcols, lchnk )
  call outfld( 'UW_sfi',         sfi,        pcols,   lchnk )

  call outfld( 'UW_chu',         chu,        pcols,   lchnk )
  call outfld( 'UW_chs',         chs,        pcols,   lchnk )
  call outfld( 'UW_cmu',         cmu,        pcols,   lchnk )
  call outfld( 'UW_cms',         cms,        pcols,   lchnk )

  call outfld( 'UW_n2',          n2,         pcols,   lchnk )
  call outfld( 'UW_s2',          s2,         pcols,   lchnk )
  call outfld( 'UW_ri',          ri,         pcols,   lchnk )

  call outfld( 'UW_kvh',         kvh,        pcols,   lchnk )
  call outfld( 'UW_kvm',         kvm,        pcols,   lchnk )
  call outfld( 'UW_pblh',        pblh,       pcols,   lchnk )
  call outfld( 'UW_ustar',       ustar,      pcols,   lchnk )
  call outfld( 'UW_pblhp',       pblhp,      pcols,   lchnk )
  call outfld( 'UW_minpblh',     minpblh,    pcols,   lchnk )

  call outfld( 'UW_tpert',       tpert,      pcols,   lchnk )
  call outfld( 'UW_qpert',       qpert,      pcols,   lchnk )
  call outfld( 'UW_tke',         tke,        pcols,   lchnk )

  call outfld( 'UW_sfuh',        sfuh,       pcols,   lchnk )
  call outfld( 'UW_sflh',        sflh,       pcols,   lchnk )

  call outfld( 'UW_ql',          qlfd,       pcols,   lchnk )

  call outfld( 'UW_tkes',        tkes,       pcols,   lchnk )
  call outfld( 'UW_wcap',        wcap,       pcols,   lchnk )
  call outfld( 'UW_wsed',        wsed,       pcols,   lchnk )
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

  ! The diffusivities from diag_TKE can be much larger than from HB in the free
  ! troposphere and upper atmosphere. These seem to be larger than observations,
  ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
  ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
  ! or the upper atmosphere.
  !
  ! NOTE: Further investigation should be done as to why the diffusivities are
  ! larger in diag_TKE.
  call eddy_diffusivity_adjustment_above_pbl_run( &
       ncol = ncol, &
       pverp = pverp, &
       kv_top_pressure = kv_top_pressure, &
       kv_freetrop_scale = kv_freetrop_scale, &
       kv_top_scale = kv_top_scale, &
       zi = state%zi(:ncol,:pverp), &
       pint = state%pint(:ncol,:pverp), &
       pblh = pblh(:ncol), &
       ! below in/out
       kvh = kvh(:ncol,:pverp), &
       kvm = kvm(:ncol,:pverp), &
       kvq = kvq(:ncol,:pverp), &
       errmsg = errmsg, errflg = errflg)

  if(errflg /= 0) then
     call endrun('eddy_diffusivity_adjustment_above_pbl_run: ' // errmsg)
  end if

end subroutine eddy_diff_tend

end module eddy_diff_cam
