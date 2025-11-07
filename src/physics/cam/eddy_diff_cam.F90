module eddy_diff_cam

use shr_kind_mod, only: i4 => shr_kind_i4, r8 => shr_kind_r8
use ppgrid, only: pcols, pver, pverp
use cam_logfile, only: iulog
use cam_abortutils, only: endrun
use physconst, only: gravit, cpair, rair, zvir, latvap, latice, karman
use eddy_diff, only: ncvmax
use time_manager, only: is_first_step
use physics_buffer, only: physics_buffer_desc
use spmd_utils, only: masterproc
use phys_control, only: phys_getopts

implicit none
private

public :: eddy_diff_readnl
public :: eddy_diff_init
public :: eddy_diff_tend

! Number of iterations for solution
integer, parameter :: nturb = 5

! Logical switches for moist mixing ratio diffusion
! (molecular diffusion is not done here)
logical :: do_diffusion_const_wet(1)
logical :: do_molecular_diffusion_const(1)

integer :: ntop_eddy, nbot_eddy

! Cloud mass constituent indices
integer :: ixcldliq, ixcldice

! input pbuf field indices
integer :: qrl_idx   = -1
integer :: wsedl_idx = -1

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

subroutine eddy_diff_init(pbuf2d, ntop_eddy_in, nbot_eddy_in)

  use error_messages, only: handle_errmsg
  use cam_history, only: addfld, add_default, horiz_only
  use constituents, only: cnst_get_ind
  use ref_pres, only: pref_mid
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

 ! Only diffuse water vapor (constituent 1) and disable molecular diffusion
  do_diffusion_const_wet(:) = .false.
  do_molecular_diffusion_const(:) = .false.
  do_diffusion_const_wet(1) = .true.

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

  type(physics_state), intent(in) :: state
  type(physics_buffer_desc), pointer, intent(in) :: pbuf(:)
  type(cam_in_t), intent(in) :: cam_in
  real(r8), intent(in) :: ztodt
  logical,         intent(in)       :: do_iss           ! Use implicit turbulent surface stress computation
  logical,         intent(in)       :: fv_am_correction    ! Do angular momentum conservation correction
  type(Coords1D), intent(in) :: p
  real(r8), intent(in) :: tint(pcols,pver+1)
  real(r8), intent(in) :: rhoi(pcols,pver+1)
  real(r8), intent(in) :: dpidz_sq(pcols,pverp)
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

  ! pbuf fields
  real(r8), pointer :: qrl(:,:)                        ! LW radiative cooling rate [K s-1]
  real(r8), pointer :: wsedl(:,:)                      ! Sedimentation velocity of stratiform liquid cloud droplet [m s-1]

  integer :: i, k
  integer :: ncol
  character(len=512)   :: errmsg
  integer              :: errflg

  ncol = state%ncol

  ! ---------------------------------------------- !
  ! Get LW radiative heating out of physics buffer !
  ! ---------------------------------------------- !
  call pbuf_get_field(pbuf, qrl_idx,   qrl)
  call pbuf_get_field(pbuf, wsedl_idx, wsedl)

  ! TODO reorder arguments of the subroutine such that in, inout, out (in this order)
  ! Call CCPPized run phase subroutine
  call compute_eddy_diff( &
       lchnk           = state%lchnk,                                  & ! to remove
       pcols           = pcols,                                        & ! to remove
       ncol            = ncol,                                         &
       pver            = pver,                                         &
       ztodt           = ztodt,                                        &
       do_iss          = do_iss,                                       &
       am_correction   = fv_am_correction,                             &
       t               = state%t(:ncol,:pver),                         &
       tint            = tint(:ncol,:pverp),                           &
       qv              = state%q(:ncol,:pver,1),                       & ! assumes q_wv at 1
       ql              = state%q(:ncol,:pver,ixcldliq),                &
       qi              = state%q(:ncol,:pver,ixcldice),                &
       s               = state%s(:ncol,:pver),                         &
       p               = p,                                            &
       rhoi            = rhoi(:ncol,:pverp),                           &
       dpidz_sq        = dpidz_sq(:ncol,:pverp),                       &
       cldn            = cldn(:ncol,:pver),                            &
       z               = state%zm(:ncol,:pver),                        &
       zi              = state%zi(:ncol,:pverp),                       &
       pmid            = state%pmid(:ncol,:pver),                      &
       pi              = state%pint(:ncol,:pverp),                     &
       u               = state%u(:ncol,:pver),                         &
       v               = state%v(:ncol,:pver),                         &
       taux            = cam_in%wsx(:ncol),                            &
       tauy            = cam_in%wsy(:ncol),                            &
       shflx           = cam_in%shf(:ncol),                            &
       qflx            = cam_in%cflx(:ncol,:1),                        &
       wstarent        = wstarent,                                     & ! use wstar entrainment? logical
       rrho            = rrho(:ncol),                                  &
       ustar           = ustar(:ncol),                                 &
       pblh            = pblh(:ncol),                                  &
       kvm_in          = kvm_in(:ncol,:pverp),                         & ! TODO hplin why _in? diff from _out? from vdiff?
       kvh_in          = kvh_in(:ncol,:pverp),                         & ! TODO hplin why _in? diff from _out? from vdiff?
       qrl             = qrl(:ncol,:pver),                             & ! pbuf fields from above (in)
       wsedl           = wsedl(:ncol,:pver),                           & ! pbuf fields from above (in)
       kvm_out         = kvm(:ncol,:pverp),                            &
       kvh_out         = kvh(:ncol,:pverp),                            &
       kvq             = kvq(:ncol,:pverp),                            &
       cgh             = cgh(:ncol,:pverp),                            &
       cgs             = cgs(:ncol,:pverp),                            &
       tpert           = tpert(:ncol),                                 &
       qpert           = qpert(:ncol),                                 &
       tke             = tke(:ncol,:pverp),                            &
       sprod           = sprod(:ncol,:pverp),                          &
       sfi             = sfi(:ncol,:pverp),                            &
       tauresx         = tauresx(:ncol),                               &
       tauresy         = tauresy(:ncol),                               &
       ksrftms         = ksrftms(:ncol),                               &
       dragblj         = dragblj(:ncol,:pver),                         &
       errmsg          = errmsg,                                       &
       errflg          = errflg)

  if(errflg /= 0) then
    call endrun('compute_eddy_diff: ' // errmsg)
  end if

  ! TODO need to pass out all relevant history field arrays as intent(out) from compute_eddy_diff if not already
  ! so that:
  ! TODO need to add the history outflds here

  ! TODO this could be made into a CCPPized scheme to do the correction
  ! (and it could do the read of kv_freetrop_scale kv_top-scale kv_top_pressure namelist options if unused elsewhere here.)

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

! Interface to UW PBL scheme to compute eddy diffusivities.
! Eddy diffusivities are calculated in a fully implicit way through iterative process.
!
! Original author: Sungsu Park, August 2006, May 2008.
subroutine compute_eddy_diff( lchnk  ,                                                      &
                              pcols  , pver   , ncol     , &
                              do_iss, am_correction, &
                              t       , tint, qv       , ztodt   ,   &
                              ql     , qi     , s        , p       , rhoi, dpidz_sq, cldn     ,             &
                              z      , zi     , pmid     , pi      , u        , v       ,         &
                              taux   , tauy   , shflx    , qflx    , wstarent ,           rrho  , &
                              ustar  , pblh   , kvm_in   , kvh_in  , &
                              qrl, wsedl, &
                              kvm_out  , kvh_out , kvq   , &
                              cgh    , cgs    , tpert    , qpert   , tke     ,                    &
                              sprod  , sfi    ,                                                   &
                              tauresx, tauresy, ksrftms, dragblj, &
                              errmsg, errflg )

  use cam_history,          only: outfld
  use atmos_phys_pbl_utils, only: calc_eddy_flux_coefficient, calc_ideal_gas_rrho, calc_friction_velocity
  use error_messages,       only: handle_errmsg
  use coords_1d,            only: Coords1D
  use wv_saturation,        only: qsat

  ! Driver routines for UW PBL scheme.
  use eddy_diff,            only: trbintd, caleddy

  ! Driver routines for diffusion solver to be called iteratively within this run phase.
  use diffusion_solver,     only: implicit_surface_stress_add_drag_coefficient_run
  use diffusion_stubs,      only: turbulent_mountain_stress_add_drag_coefficient_run
  use diffusion_solver,     only: vertical_diffusion_wind_damping_rate_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_horizontal_momentum_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_dry_static_energy_run
  use diffusion_solver,     only: vertical_diffusion_diffuse_tracers_run

  use beljaars_drag_cam,    only: do_beljaars

  ! Input variables
  integer,  intent(in)    :: lchnk                       ! Chunk identifier [index]
  integer,  intent(in)    :: pcols                       ! Number of atmospheric columns [#]

  integer,        intent(in) :: pver                ! Number of atmospheric layers [#]
  integer,        intent(in) :: ncol                ! Number of atmospheric columns [#]
  logical,        intent(in) :: do_iss              ! Use implicit turbulent surface stress computation [flag]
  logical,        intent(in) :: am_correction       ! Do angular momentum conservation correction [flag]
  real(r8),       intent(in) :: t(:,:)              ! Temperature [K]
  real(r8),       intent(in) :: tint(:,:)           ! Temperature defined on interfaces [K]
  real(r8),       intent(in) :: qv(:,:)             ! Water vapor specific humidity [kg kg-1]
  real(r8),       intent(in) :: ztodt               ! Physics integration time step 2 delta-t [s]
  real(r8),       intent(in) :: ql(:,:)             ! Liquid water specific humidity [kg kg-1]
  real(r8),       intent(in) :: qi(:,:)             ! Ice specific humidity [kg kg-1]
  real(r8),       intent(in) :: s(:,:)              ! Dry static energy [J kg-1]
  type(Coords1D), intent(in) :: p                   ! Pressure coordinates for solver [Pa]
  real(r8),       intent(in) :: rhoi(:,:)           ! Density at interfaces [kg m-3]
  real(r8),       intent(in) :: dpidz_sq(:,:)       ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4], interfaces
  real(r8),       intent(in) :: cldn(:,:)           ! Stratiform cloud fraction [fraction]
  real(r8),       intent(in) :: z(:,:)              ! Layer mid-point height above surface [m]
  real(r8),       intent(in) :: zi(:,:)             ! Interface height above surface [m]
  real(r8),       intent(in) :: pmid(:,:)           ! Layer mid-point pressure [Pa]
  real(r8),       intent(in) :: pi(:,:)             ! Interface pressure [Pa]
  real(r8),       intent(in) :: u(:,:)              ! Zonal velocity [m s-1]
  real(r8),       intent(in) :: v(:,:)              ! Meridional velocity [m s-1]
  real(r8),       intent(in) :: taux(:)             ! Zonal wind stress at surface [N m-2]
  real(r8),       intent(in) :: tauy(:)             ! Meridional wind stress at surface [N m-2]
  real(r8),       intent(in) :: shflx(:)            ! Sensible heat flux at surface [W m-2]
  real(r8),       intent(in) :: qflx(:,:)           ! Water vapor flux at surface [kg m-2 s-1]
  logical,        intent(in) :: wstarent            ! .true. means use the 'wstar' entrainment closure [flag]
  real(r8),       intent(in) :: kvm_in(:,:)         ! kvm saved from last timestep [m2 s-1]
  real(r8),       intent(in) :: kvh_in(:,:)         ! kvh saved from last timestep [m2 s-1]
  real(r8),       intent(in) :: ksrftms(:)          ! Surface drag coefficient of turbulent mountain stress [kg m-2 s-1]
  real(r8),       intent(in) :: dragblj(:,:)        ! Drag profile from Beljaars SGO form drag [s-1]
  real(r8),       intent(in) :: qrl(:,:)            ! LW radiative cooling rate [K s-1]
  real(r8),       intent(in) :: wsedl(:,:)          ! Sedimentation velocity of stratiform liquid cloud droplet [m s-1]

  ! Input/output variables
  real(r8),    intent(inout) :: tauresx(:)          ! Residual stress to be added in vdiff to correct for turb [N m-2]
  real(r8),    intent(inout) :: tauresy(:)          ! Stress mismatch between sfc and atm accumulated in prior timesteps [N m-2]

  ! Output variables
  real(r8),      intent(out) :: kvm_out(:,:)        ! Eddy diffusivity for momentum [m2 s-1]
  real(r8),      intent(out) :: kvh_out(:,:)        ! Eddy diffusivity for heat [m2 s-1]
  real(r8),      intent(out) :: kvq(:,:)            ! Eddy diffusivity for constituents, moisture and tracers [m2 s-1]
  real(r8),      intent(out) :: rrho(:)             ! Reciprocal of density at the lowest layer [m3 kg-1]
  real(r8),      intent(out) :: ustar(:)            ! Surface friction velocity [m s-1]
  real(r8),      intent(out) :: pblh(:)             ! PBL top height [m]
  real(r8),      intent(out) :: cgh(:,:)            ! Counter-gradient term for heat [J kg-1 m-1]
  real(r8),      intent(out) :: cgs(:,:)            ! Counter-gradient star [cg flux-1]
  real(r8),      intent(out) :: tpert(:)            ! Convective temperature excess [K]
  real(r8),      intent(out) :: qpert(:)            ! Convective humidity excess [kg kg-1]
  real(r8),      intent(out) :: tke(:,:)            ! Turbulent kinetic energy [m2 s-2]
  real(r8),      intent(out) :: sprod(:,:)          ! Shear production [m2 s-3]
  real(r8),      intent(out) :: sfi(:,:)            ! Interfacial layer saturation fraction [fraction]
  character(len=512), intent(out) :: errmsg         ! Error message
  integer,            intent(out) :: errflg         ! Error flag

  ! --------------- !
  ! Local Variables !
  ! --------------- !
  integer                    icol
  integer                    i, k, iturb, status
  integer :: ipbl(ncol)     ! If 1, PBL is CL, while if 0, PBL is STL.
  integer :: kpblh(ncol)    ! Layer index containing PBL top within or at the base interface (NOT USED)

  character(2048)         :: warnstring                ! Warning(s) to print
  character(128)          :: errstring                 ! Error message

  real(r8) :: bprod(ncol,pverp) ! Buoyancy production of tke [ m2/s3 ]
  real(r8) :: tkes(ncol)        ! TKE at surface interface [ m2/s2 ]
  real(r8) :: went(ncol)        ! Entrainment rate at the PBL top interface [ m/s ] (NOT USED)

  real(r8)                :: kvf(ncol,pver+1)         ! Free atmospheric eddy diffusivity [ m2/s ]
  real(r8)                :: kvm(ncol,pver+1)         ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh(ncol,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: kvm_preo(ncol,pver+1)    ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh_preo(ncol,pver+1)    ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: kvm_pre(ncol,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
  real(r8)                :: kvh_pre(ncol,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
  real(r8)                :: errorPBL(ncol)           ! Error function showing whether PBL produced convergent solution or not.
                                                       ! [ unit ? ]
  real(r8)                :: s2(ncol,pver)            ! Shear squared, defined at interfaces except surface [ s-2 ]
  real(r8)                :: n2(ncol,pver)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
  real(r8)                :: ri(ncol,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
  real(r8)                :: pblhp(ncol)              ! PBL top pressure [ Pa ]
  real(r8)                :: minpblh(ncol)            ! Minimum PBL height based on surface stress

  real(r8)                :: qt(ncol,pver)            ! Total specific humidity [ kg/kg ]
  real(r8)                :: sfuh(ncol,pver)          ! Saturation fraction in upper half-layer [ fraction ]
  real(r8)                :: sflh(ncol,pver)          ! Saturation fraction in lower half-layer [ fraction ]
  real(r8)                :: sl(ncol,pver)            ! Liquid water static energy [ J/kg ]
  real(r8)                :: slv(ncol,pver)           ! Liquid water virtual static energy [ J/kg ]
  real(r8)                :: slslope(ncol,pver)       ! Slope of 'sl' in each layer
  real(r8)                :: qtslope(ncol,pver)       ! Slope of 'qt' in each layer
  real(r8)                :: qvfd(ncol,pver)          ! Specific humidity for diffusion [ kg/kg ]
  real(r8)                :: tfd(ncol,pver)           ! Temperature for diffusion [ K ]
  real(r8)                :: slfd(ncol,pver)          ! Liquid static energy [ J/kg ]
  real(r8)                :: qtfd(ncol,pver,1)        ! Total specific humidity [ kg/kg ]
  real(r8)                :: qlfd(ncol,pver)          ! Liquid water specific humidity for diffusion [ kg/kg ]
  real(r8)                :: ufd(ncol,pver)           ! U-wind for diffusion [ m/s ]
  real(r8)                :: vfd(ncol,pver)           ! V-wind for diffusion [ m/s ]

  ! Output arguments from iterative calls of diffusion solver
  real(r8)                :: ufd_out(ncol,pver)       ! U-wind for diffusion [ m/s ]
  real(r8)                :: vfd_out(ncol,pver)       ! V-wind for diffusion [ m/s ]
  real(r8)                :: slfd_out(ncol,pver)      ! Liquid static energy [ J/kg ]
  real(r8)                :: qtfd_out(ncol,pver,1)    ! Total specific humidity [ kg/kg ]

  ! Input arguments for CCPPized diffusion solver
  real(r8)                :: ksrf(ncol)
  real(r8)                :: tau_damp_rate(ncol, pver)
  real(r8)                :: ubc_mmr_dummy(ncol, 1)
  logical                 :: cnst_fixed_ubc(1)

  ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'

  real(r8)                :: chu(ncol,pver+1)         ! Heat buoyancy coef for dry states, defined at each interface, finally.
  real(r8)                :: chs(ncol,pver+1)         ! Heat buoyancy coef for sat states, defined at each interface, finally.
  real(r8)                :: cmu(ncol,pver+1)         ! Moisture buoyancy coef for dry states,
                                                       ! defined at each interface, finally.
  real(r8)                :: cms(ncol,pver+1)         ! Moisture buoyancy coef for sat states,
                                                       ! defined at each interface, finally.

  real(r8)                :: jnk1d(ncol)
  real(r8)                :: jnk2d(ncol,pver+1)
  real(r8)                :: zero(ncol)
  real(r8)                :: zero2d(ncol,pver+1)
  real(r8)                :: es                     ! Saturation vapor pressure
  real(r8)                :: qs                     ! Saturation specific humidity
  real(r8)                :: ep2, templ, temps

  ! ------------------------------- !
  ! Variables for diagnostic output !
  ! ------------------------------- !

  real(r8)                :: wpert(ncol)              ! Turbulent velocity excess [ m/s ]

  real(r8)                :: kbase_o(ncol,ncvmax)     ! Original external base interface index of CL from 'exacol'
  real(r8)                :: ktop_o(ncol,ncvmax)      ! Original external top  interface index of CL from 'exacol'
  real(r8)                :: ncvfin_o(ncol)           ! Original number of CLs from 'exacol'
  real(r8)                :: kbase_mg(ncol,ncvmax)    ! 'kbase' after extending-merging from 'zisocl'
  real(r8)                :: ktop_mg(ncol,ncvmax)     ! 'ktop' after extending-merging from 'zisocl'
  real(r8)                :: ncvfin_mg(ncol)          ! 'ncvfin' after extending-merging from 'zisocl'
  real(r8)                :: kbase_f(ncol,ncvmax)     ! Final 'kbase' after extending-merging & including SRCL
  real(r8)                :: ktop_f(ncol,ncvmax)      ! Final 'ktop' after extending-merging & including SRCL
  real(r8)                :: ncvfin_f(ncol)           ! Final 'ncvfin' after extending-merging & including SRCL
  real(r8)                :: wet(ncol,ncvmax)         ! Entrainment rate at the CL top  [ m/s ]
  real(r8)                :: web(ncol,ncvmax)         ! Entrainment rate at the CL base [ m/s ].
                                                       ! Set to zero if CL is based at surface.
  real(r8)                :: jtbu(ncol,ncvmax)        ! Buoyancy jump across the CL top  [ m/s2 ]
  real(r8)                :: jbbu(ncol,ncvmax)        ! Buoyancy jump across the CL base [ m/s2 ]
  real(r8)                :: evhc(ncol,ncvmax)        ! Evaporative enhancement factor at the CL top
  real(r8)                :: jt2slv(ncol,ncvmax)      ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
  real(r8)                :: n2ht(ncol,ncvmax)        ! n2 defined at the CL top  interface but using
                                                       ! sfuh(kt)   instead of sfi(kt) [ s-2 ]
  real(r8)                :: n2hb(ncol,ncvmax)        ! n2 defined at the CL base interface but using
                                                       ! sflh(kb-1) instead of sfi(kb) [ s-2 ]
  real(r8)                :: lwp(ncol,ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
  real(r8)                :: opt_depth(ncol,ncvmax)   ! Optical depth of the CL top layer
  real(r8)                :: radinvfrac(ncol,ncvmax)  ! Fraction of radiative cooling confined in the top portion of CL top layer
  real(r8)                :: radf(ncol,ncvmax)        ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
  real(r8)                :: wstar(ncol,ncvmax)       ! Convective velocity in each CL [ m/s ]
  real(r8)                :: wstar3fact(ncol,ncvmax)  ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
  real(r8)                :: ebrk(ncol,ncvmax)        ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
  real(r8)                :: wbrk(ncol,ncvmax)        ! Net mean normalized TKE (W) of CL,
                                                       ! 'ebrk/b1' including entrainment effect [ m2/s2 ]
  real(r8)                :: lbrk(ncol,ncvmax)        ! Energetic internal thickness of CL [m]
  real(r8)                :: ricl(ncol,ncvmax)        ! CL internal mean Richardson number
  real(r8)                :: ghcl(ncol,ncvmax)        ! Half of normalized buoyancy production of CL
  real(r8)                :: shcl(ncol,ncvmax)        ! Galperin instability function of heat-moisture of CL
  real(r8)                :: smcl(ncol,ncvmax)        ! Galperin instability function of mementum of CL
  real(r8)                :: ghi(ncol,pver+1)         ! Half of normalized buoyancy production at all interfaces
  real(r8)                :: shi(ncol,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
  real(r8)                :: smi(ncol,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
  real(r8)                :: rii(ncol,pver+1)         ! Interfacial Richardson number defined at all interfaces
  real(r8)                :: lengi(ncol,pver+1)       ! Turbulence length scale at all interfaces [ m ]
  real(r8)                :: wcap(ncol,pver+1)        ! Normalized TKE at all interfaces [ m2/s2 ]
  ! For sedimentation-entrainment feedback
  real(r8)                :: wsed(ncol,ncvmax)        ! Sedimentation velocity at the top of each CL [ m/s ]

  integer(i4)             :: turbtype(ncol,pver+1)    ! Turbulence type identifier at all interfaces [ no unit ]


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
     rrho(:ncol)   = calc_ideal_gas_rrho(rair, tfd(:ncol,pver), pmid(:ncol,pver))
     ustar(:ncol)  = calc_friction_velocity(taux(:ncol) - ksrftms(:ncol) * ufd(:ncol,pver), & ! Zonal wind stress
                                            tauy(:ncol) - ksrftms(:ncol) * vfd(:ncol,pver), & ! Meridional wind stress
                                            rrho(:ncol))

     minpblh(:ncol) = 100.0_r8 * ustar(:ncol)   ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'.

     ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)

     call trbintd( &
                   ncol    , pver    , ncol  , z       , ufd     , vfd     , tfd   , pmid    , &
                   s2       , n2      , ri    , zi      , pi      , cldn    , qtfd  , qvfd    , &
                   qlfd     , qi      , sfi   , sfuh    , sflh    , slfd    , slv   , slslope , &
                   qtslope  , chs     , chu   , cms     , cmu     )

     ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.
     ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables.

     if( iturb == 1 ) then
        qt(:ncol,:) = qtfd(:ncol,:,1)
        sl(:ncol,:) = slfd(:ncol,:)
     endif

     ! Get free atmosphere exchange coefficients. This 'kvf' is not used in UW moist PBL scheme
     if (use_kvf) then
        kvf(:ncol,:) = 0.0_r8
        do k = ntop_eddy, nbot_eddy-1
           do i = 1, ncol
              kvf(i,k+1) = calc_eddy_flux_coefficient(ml2(k), ri(i, k), s2(i, k))
           end do
        end do
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

     write(6, *) "b caleddy iturb = ", iturb

     call caleddy( ncol      , pver      , ncol      ,                     &
                   slfd      , qtfd      , qlfd      , slv      ,ufd     , &
                   vfd       , pi        , z         , zi       ,          &
                   qflx      , shflx     , slslope   , qtslope  ,          &
                   chu       , chs       , cmu       , cms      ,sfuh    , &
                   sflh      , n2        , s2        , ri       ,rrho    , &
                   pblh      , ustar     ,                                 &
                   kvh       , kvm       , kvh_out   , kvm_out  ,          &
                   tpert     , qpert     , qrl       , kvf      , tke    , &
                   wstarent  , bprod     , sprod     , minpblh  , wpert  , &
                   tkes      , went      , turbtype  ,                     &
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
        qtfd(:ncol,:,1)= qt(:ncol,:)
        ufd(:ncol,:)   = u(:ncol,:)
        vfd(:ncol,:)   = v(:ncol,:)

        ! TODO hplin: calculation of k and drag coef could probably be moved to outside the iterative loop as u/v is orig flx
        ! Calculate surface drag rate
        ksrf(:ncol) = 0._r8
        call implicit_surface_stress_add_drag_coefficient_run( &
             ncol            = ncol,                         &
             pver            = pver,                         &
             do_iss          = do_iss,                       &
             taux            = taux(:ncol),                  &
             tauy            = tauy(:ncol),                  &
             u0              = ufd(:ncol,:pver),             &
             v0              = vfd(:ncol,:pver),             &
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

        ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
        ! Call the CCPPized subroutines for the diffusion solver
        ! in iteration. This is not specified in the SDF but instead
        ! called internally because it is a iterative process.
        ! Notes:
        ! - there is no molecular diffusion used here.
        ! - in tracers, only water vapor is diffused (ncnst = 1)
        ufd_out(:,:) = 0._r8
        vfd_out(:,:) = 0._r8
        slfd_out(:,:) = 0._r8
        call vertical_diffusion_diffuse_horizontal_momentum_run( &
            ncol            = ncol,                         &
            pver            = pver,                         &
            pverp           = pverp,                        &
            dt              = ztodt,                        &
            rair            = rair,                         &
            gravit          = gravit,                       &
            do_iss          = do_iss,                       &
            am_correction   = am_correction,                &
            itaures         = .false.,                      &
            t               = t(:ncol,:pver),               &
            p               = p,                            & ! Coords1D, pressure coordinates [Pa]
            rhoi            = rhoi(:ncol,:pverp),           &
            taux            = taux(:ncol),                  &
            tauy            = tauy(:ncol),                  &
            tau_damp_rate   = tau_damp_rate(:ncol,:pver),   & ! tau damp rate from above
            kvm             = kvm_out(:ncol,:pverp),        &
            ksrftms         = ksrftms(:ncol),               &
            dragblj         = dragblj(:ncol,:pver),         &
            dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist
            u0              = ufd(:ncol,:pver),             &
            v0              = vfd(:ncol,:pver),             &
            dse0            = slfd(:ncol,:pver),            &
            ! input/output
            tauresx         = tauresx(:ncol),               &
            tauresy         = tauresy(:ncol),               &
            ! below output
            u1              = ufd_out(:ncol,:pver),         &
            v1              = vfd_out(:ncol,:pver),         &
            dse1            = slfd_out(:ncol,:pver),        &
            dtk             = jnk2d(:ncol,:),               & ! unused dummy.
            tautmsx         = jnk1d(:ncol),                 & ! unused dummy.
            tautmsy         = jnk1d(:ncol),                 & ! unused dummy.
            ! arguments for Beljaars
            do_beljaars     = do_beljaars,                  &
            errmsg          = errmsg,                       &
            errflg          = errflg)

        if(errflg /= 0) then
           call endrun('vertical_diffusion_diffuse_horizontal_momentum_run: ' // errmsg)
        endif

        ! Update u, v, dse with updated iterative values after diffusion solver
        ufd(:, :pver)  = ufd_out(:, :pver)
        vfd(:, :pver)  = vfd_out(:, :pver)
        slfd(:, :pver) = slfd_out(:, :pver)

        ! Diffuse dry static energy
        call vertical_diffusion_diffuse_dry_static_energy_run( &
             ncol            = ncol,                         &
             pver            = pver,                         &
             dt              = ztodt,                        &
             gravit          = gravit,                       &
             p               = p,                            & ! Coords1D, pressure coordinates [Pa]
             rhoi            = rhoi(:ncol,:pverp),           &
             shflx           = shflx(:ncol),                 &
             dse_top         = zero(:ncol),                  & ! = zero
             kvh             = kvh_out(:ncol,:pverp),        &
             cgh             = cgh(:ncol,:pverp),            &
             dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist
             ! input/output
             dse             = slfd(:ncol,:pver),            &
             errmsg          = errmsg,                       &
             errflg          = errflg)

        if(errflg /= 0) then
           call endrun('vertical_diffusion_diffuse_dry_static_energy_run: ' // errmsg)
        endif

        ! Diffuse tracers
        ubc_mmr_dummy(:ncol,:1) = 0._r8
        cnst_fixed_ubc(:1) = .false.

        qtfd_out(:,:,:) = 0._r8
        call vertical_diffusion_diffuse_tracers_run( &
             ncol            = ncol,                         &
             pver            = pver,                         &
             ncnst           = 1,                            & ! only water vapor is diffused here.
             dt              = ztodt,                        &
             rair            = rair,                         &
             gravit          = gravit,                       &
             do_diffusion_const = do_diffusion_const_wet,    & ! moist constituents to diffuse
             p               = p,                            & ! Coords1D, pressure coordinates [Pa]
             t               = t(:ncol,:pver),               &
             rhoi            = rhoi(:ncol,:pverp),           &
             cflx            = qflx(:ncol,:1),               & ! wv only. WARN: assumes wv at 1
             kvh             = kvh_out(:ncol,:pverp),        &
             kvq             = kvh_out(:ncol,:pverp),        & ! [sic] kvh_out is assigned to kvh, kvq. check
             cgs             = cgs(:ncol,:pverp),            &
             qmincg          = zero(:ncol),                  &
             dpidz_sq        = dpidz_sq(:ncol,:pverp),       & ! moist TODO
             ! upper boundary conditions from ubc module
             ubc_mmr         = ubc_mmr_dummy(:ncol,:1),      &
             cnst_fixed_ubc  = cnst_fixed_ubc(:1),           & ! = .false.
             q0              = qtfd(:ncol,:pver,:1),         &
             q1              = qtfd_out(:ncol,:pver,:1),     &
             errmsg          = errmsg,                       &
             errflg          = errflg)

        if(errflg /= 0) then
           call endrun('vertical_diffusion_diffuse_tracers_run: ' // errmsg)
        endif

        ! update q with after in iterative process.
        qtfd(:ncol, :pver, :1) = qtfd_out(:ncol, :pver, :1)

        ! TODO (hplin, 5/9/2025): after these are subset to ncol check if we
        ! need to initialize some outs to 0; compute_vdiff did not do this before

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
              temps     =   templ + ( qtfd(i,k,1) - qs ) / ( cpair / latvap + latvap * qs / ( rair * templ**2 ) )
              call qsat( temps, pmid(i,k), es, qs)
              qlfd(i,k) =   max( qtfd(i,k,1) - qi(i,k) - qs ,0._r8 )
              ! Option.2 : Assume condensate is not diffused by the moist turbulence scheme.
              !            This should bs used if 'pseudodiff = .true.'  in vertical_diffusion.F90.
              ! qlfd(i,k) = ql(i,k)
              ! ----------------------------- !
              ! Compute the other 'qvfd, tfd' !
              ! ----------------------------- !
              qvfd(i,k) = max( 0._r8, qtfd(i,k,1) - qi(i,k) - qlfd(i,k) )
              tfd(i,k)  = ( slfd(i,k) + latvap * qlfd(i,k) + (latvap+latice) * qi(i,k) - gravit*z(i,k)) / cpair
           end do
        end do
     endif

  end do  ! End of 'iturb' iteration

  kvq(:ncol,:) = kvh_out(:ncol,:)

  ! --------------------------------------------------------------- !
  ! Writing for detailed diagnostic analysis of UW moist PBL scheme !
  ! --------------------------------------------------------------- !

  ! call outfld( 'WGUSTD' , wpert, pcols, lchnk )

  ! call outfld( 'BPROD   ', bprod, pcols, lchnk )
  ! call outfld( 'SPROD   ', sprod, pcols, lchnk )
  ! call outfld( 'SFI     ', sfi,   pcols, lchnk )

  ! call outfld( 'UW_errorPBL',    errorPBL,   pcols,   lchnk )

  ! call outfld( 'UW_n2',          n2,         pcols,   lchnk )
  ! call outfld( 'UW_s2',          s2,         pcols,   lchnk )
  ! call outfld( 'UW_ri',          ri,         pcols,   lchnk )

  ! call outfld( 'UW_sfuh',        sfuh,       pcols,   lchnk )
  ! call outfld( 'UW_sflh',        sflh,       pcols,   lchnk )
  ! call outfld( 'UW_sfi',         sfi,        pcols,   lchnk )

  ! call outfld( 'UW_cldn',        cldn,       pcols,   lchnk )
  ! call outfld( 'UW_qrl',         qrl,        pcols,   lchnk )
  ! call outfld( 'UW_ql',          qlfd,       pcols,   lchnk )

  ! call outfld( 'UW_chu',         chu,        pcols,   lchnk )
  ! call outfld( 'UW_chs',         chs,        pcols,   lchnk )
  ! call outfld( 'UW_cmu',         cmu,        pcols,   lchnk )
  ! call outfld( 'UW_cms',         cms,        pcols,   lchnk )

  ! call outfld( 'UW_tke',         tke,        pcols,   lchnk )
  ! call outfld( 'UW_wcap',        wcap,       pcols,   lchnk )
  ! call outfld( 'UW_bprod',       bprod,      pcols,   lchnk )
  ! call outfld( 'UW_sprod',       sprod,      pcols,   lchnk )

  ! call outfld( 'UW_kvh',         kvh_out,    pcols,   lchnk )
  ! call outfld( 'UW_kvm',         kvm_out,    pcols,   lchnk )

  ! call outfld( 'UW_pblh',        pblh,       pcols,   lchnk )
  ! call outfld( 'UW_pblhp',       pblhp,      pcols,   lchnk )
  ! call outfld( 'UW_tpert',       tpert,      pcols,   lchnk )
  ! call outfld( 'UW_qpert',       qpert,      pcols,   lchnk )
  ! call outfld( 'UW_wpert',       wpert,      pcols,   lchnk )

  ! call outfld( 'UW_ustar',       ustar,      pcols,   lchnk )
  ! call outfld( 'UW_tkes',        tkes,       pcols,   lchnk )
  ! call outfld( 'UW_minpblh',     minpblh,    pcols,   lchnk )

  ! call outfld( 'UW_turbtype',    real(turbtype,r8),   pcols,   lchnk )

  ! call outfld( 'UW_kbase_o',     kbase_o,    pcols,   lchnk )
  ! call outfld( 'UW_ktop_o',      ktop_o,     pcols,   lchnk )
  ! call outfld( 'UW_ncvfin_o',    ncvfin_o,   pcols,   lchnk )

  ! call outfld( 'UW_kbase_mg',    kbase_mg,   pcols,   lchnk )
  ! call outfld( 'UW_ktop_mg',     ktop_mg,    pcols,   lchnk )
  ! call outfld( 'UW_ncvfin_mg',   ncvfin_mg,  pcols,   lchnk )

  ! call outfld( 'UW_kbase_f',     kbase_f,    pcols,   lchnk )
  ! call outfld( 'UW_ktop_f',      ktop_f,     pcols,   lchnk )
  ! call outfld( 'UW_ncvfin_f',    ncvfin_f,   pcols,   lchnk )

  ! call outfld( 'UW_wet',         wet,        pcols,   lchnk )
  ! call outfld( 'UW_web',         web,        pcols,   lchnk )
  ! call outfld( 'UW_jtbu',        jtbu,       pcols,   lchnk )
  ! call outfld( 'UW_jbbu',        jbbu,       pcols,   lchnk )
  ! call outfld( 'UW_evhc',        evhc,       pcols,   lchnk )
  ! call outfld( 'UW_jt2slv',      jt2slv,     pcols,   lchnk )
  ! call outfld( 'UW_n2ht',        n2ht,       pcols,   lchnk )
  ! call outfld( 'UW_n2hb',        n2hb,       pcols,   lchnk )
  ! call outfld( 'UW_lwp',         lwp,        pcols,   lchnk )
  ! call outfld( 'UW_optdepth',    opt_depth,  pcols,   lchnk )
  ! call outfld( 'UW_radfrac',     radinvfrac, pcols,   lchnk )
  ! call outfld( 'UW_radf',        radf,       pcols,   lchnk )
  ! call outfld( 'UW_wstar',       wstar,      pcols,   lchnk )
  ! call outfld( 'UW_wstar3fact',  wstar3fact, pcols,   lchnk )
  ! call outfld( 'UW_ebrk',        ebrk,       pcols,   lchnk )
  ! call outfld( 'UW_wbrk',        wbrk,       pcols,   lchnk )
  ! call outfld( 'UW_lbrk',        lbrk,       pcols,   lchnk )
  ! call outfld( 'UW_ricl',        ricl,       pcols,   lchnk )
  ! call outfld( 'UW_ghcl',        ghcl,       pcols,   lchnk )
  ! call outfld( 'UW_shcl',        shcl,       pcols,   lchnk )
  ! call outfld( 'UW_smcl',        smcl,       pcols,   lchnk )

  ! call outfld( 'UW_gh',          ghi,        pcols,   lchnk )
  ! call outfld( 'UW_sh',          shi,        pcols,   lchnk )
  ! call outfld( 'UW_sm',          smi,        pcols,   lchnk )
  ! call outfld( 'UW_ria',         rii,        pcols,   lchnk )
  ! call outfld( 'UW_leng',        lengi,      pcols,   lchnk )

  ! call outfld( 'UW_wsed',        wsed,       pcols,   lchnk )

end subroutine compute_eddy_diff

end module eddy_diff_cam
