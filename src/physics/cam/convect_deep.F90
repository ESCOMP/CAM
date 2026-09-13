
module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to several deep convection interfaces. Currently includes:
!    Zhang-McFarlane (default)
!    Kerry Emanuel
!
!
! Author: D.B. Coleman, Sep 2004
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp
   use cam_logfile,  only: iulog

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend,               &! return tendencies
      convect_deep_tend_2,             &! return tendencies
      deep_scheme_does_scav_trans             ! = .t. if scheme does scavenging and conv. transport

! Private module data
   character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change
! Physics buffer indices
   integer     ::  icwmrdp_idx      = 0
   integer     ::  rprddp_idx       = 0
   integer     ::  nevapr_dpcu_idx  = 0
   integer     ::  cldtop_idx       = 0
   integer     ::  cldbot_idx       = 0
   integer     ::  cld_idx          = 0
   integer     ::  fracis_idx       = 0

   integer     ::  pblh_idx        = 0
   integer     ::  tpert_idx       = 0
   integer     ::  prec_dp_idx     = 0
   integer     ::  snow_dp_idx     = 0

   integer     ::  ttend_dp_idx        = 0

   ! pbuf indices of the ZM gathered arrays, registered below when
   ! deep_scheme='CLUBB_MF' (the CLUBB-MF plume ensemble provides the deep
   ! transport); populated by clubb_tend_cam, used by clubb_mf_convtran2 below
   integer     ::  zm_mu_idx       = 0
   integer     ::  zm_eu_idx       = 0
   integer     ::  zm_du_idx       = 0
   integer     ::  zm_md_idx       = 0
   integer     ::  zm_ed_idx       = 0
   integer     ::  zm_dp_idx       = 0
   integer     ::  zm_dsubcld_idx  = 0
   integer     ::  zm_jt_idx       = 0
   integer     ::  zm_maxg_idx     = 0
   integer     ::  zm_ideep_idx    = 0
   ! plume-ensemble mean updraft speed (m/s, gathered like ZM_MU;
   ! populated by clubb_tend_cam, consumed by aero_convproc activation)
   integer     ::  mf_wup_idx      = 0

!=========================================================================================
  contains

!=========================================================================================
function deep_scheme_does_scav_trans()
!
! Function called by tphysbc to determine if it needs to do scavenging and convective transport
! or if those have been done by the deep convection scheme. Each scheme could have its own
! identical query function for a less-knowledgable interface but for now, we know that KE
! does scavenging & transport, and ZM doesn't
!

  logical deep_scheme_does_scav_trans

  deep_scheme_does_scav_trans = .false.

  if ( deep_scheme .eq. 'KE' ) deep_scheme_does_scav_trans = .true.

  return

end function deep_scheme_does_scav_trans

!=========================================================================================
subroutine convect_deep_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------


  use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4
  use zm_conv_intr, only: zm_conv_register
  use phys_control, only: phys_getopts, use_gw_convect_dp

  implicit none

  integer idx

  ! get deep_scheme setting from phys_control
  call phys_getopts(deep_scheme_out = deep_scheme)

  select case ( deep_scheme )
  case('ZM') !    Zhang-McFarlane (default)
     call zm_conv_register

  case('off') ! Off needs to setup the following fields
   call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
   call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
   call pbuf_add_field('SNOW_DP',    'physpkg',dtype_r8,(/pcols/),     snow_dp_idx)

  case('CLUBB_MF') ! The standard deep fields and the ZM gathered arrays are registered here
   call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
   call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
   call pbuf_add_field('SNOW_DP',    'physpkg',dtype_r8,(/pcols/),     snow_dp_idx)

   call pbuf_add_field('ZM_MU',      'physpkg', dtype_r8, (/pcols,pver/), zm_mu_idx)
   call pbuf_add_field('ZM_EU',      'physpkg', dtype_r8, (/pcols,pver/), zm_eu_idx)
   call pbuf_add_field('ZM_DU',      'physpkg', dtype_r8, (/pcols,pver/), zm_du_idx)
   call pbuf_add_field('ZM_MD',      'physpkg', dtype_r8, (/pcols,pver/), zm_md_idx)
   call pbuf_add_field('ZM_ED',      'physpkg', dtype_r8, (/pcols,pver/), zm_ed_idx)
   call pbuf_add_field('ZM_DP',      'physpkg', dtype_r8, (/pcols,pver/), zm_dp_idx)
   call pbuf_add_field('ZM_DSUBCLD', 'physpkg', dtype_r8, (/pcols/),      zm_dsubcld_idx)
   call pbuf_add_field('ZM_JT',      'physpkg', dtype_i4, (/pcols/),      zm_jt_idx)
   call pbuf_add_field('ZM_MAXG',    'physpkg', dtype_i4, (/pcols/),      zm_maxg_idx)
   call pbuf_add_field('ZM_IDEEP',   'physpkg', dtype_i4, (/pcols/),      zm_ideep_idx)
   call pbuf_add_field('MF_WUP',     'physpkg', dtype_r8, (/pcols,pver/), mf_wup_idx)

  end select

  ! If gravity waves from deep convection are on, output this field.
  if (use_gw_convect_dp .and. deep_scheme == 'ZM') then
     call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,(/pcols,pver/),ttend_dp_idx)
  end if

end subroutine convect_deep_register

!=========================================================================================



subroutine convect_deep_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld
  use pmgrid,         only: plevp
  use spmd_utils,     only: masterproc
  use zm_conv_intr,   only: zm_conv_init
  use cam_abortutils, only: endrun

  use physics_buffer, only: physics_buffer_desc, pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces

  select case ( deep_scheme )
  case('off')
     if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
  case('CLUBB_SGS')
     if (masterproc) write(iulog,*)'convect_deep: CLUBB_SGS selected'
  case('CLUBB_MF')
     if (masterproc) write(iulog,*) &
        'convect_deep: CLUBB-MF selected: the EDMF plume ensemble drives '// &
        'aero_convproc and convtran2 through the ZM pbuf arrays'
  case('ZM')
     if (masterproc) write(iulog,*)'convect_deep initializing Zhang-McFarlane convection'
     call zm_conv_init(pref_edge)
  case default
     if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme. May fail.'
  end select

  icwmrdp_idx     = pbuf_get_index('ICWMRDP')
  rprddp_idx      = pbuf_get_index('RPRDDP')
  nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
  prec_dp_idx     = pbuf_get_index('PREC_DP')
  snow_dp_idx     = pbuf_get_index('SNOW_DP')

  cldtop_idx = pbuf_get_index('CLDTOP')
  cldbot_idx = pbuf_get_index('CLDBOT')
  cld_idx    = pbuf_get_index('CLD')
  fracis_idx = pbuf_get_index('FRACIS')

  pblh_idx   = pbuf_get_index('pblh')
  tpert_idx  = pbuf_get_index('tpert')

  if (trim(deep_scheme) == 'CLUBB_MF') then
     if (masterproc) write(iulog,*) &
        'convect_deep: CLUBB-MF active: convtran2 will transport '// &
        'constituents using the MF plume ensemble mass fluxes'
     zm_mu_idx      = pbuf_get_index('ZM_MU')
     zm_eu_idx      = pbuf_get_index('ZM_EU')
     zm_du_idx      = pbuf_get_index('ZM_DU')
     zm_md_idx      = pbuf_get_index('ZM_MD')
     zm_ed_idx      = pbuf_get_index('ZM_ED')
     zm_dp_idx      = pbuf_get_index('ZM_DP')
     zm_dsubcld_idx = pbuf_get_index('ZM_DSUBCLD')
     zm_jt_idx      = pbuf_get_index('ZM_JT')
     zm_maxg_idx    = pbuf_get_index('ZM_MAXG')
     zm_ideep_idx   = pbuf_get_index('ZM_IDEEP')
  end if

  call addfld ('ICWMRDP', (/ 'lev' /), 'A', 'kg/kg', 'Deep Convection in-cloud water mixing ratio ' )

end subroutine convect_deep_init
!=========================================================================================
!subroutine convect_deep_tend(state, ptend, tdt, pbuf)

subroutine convect_deep_tend( &
     mcon    ,cme     ,          &
     zdu      , &
     rliq    , &
     ztodt   , &
     state   ,ptend   ,landfrac ,pbuf)


   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init

   use cam_history,    only: outfld
   use constituents,   only: pcnst
   use zm_conv_intr,   only: zm_conv_tend
   use cam_history,    only: outfld
   use physconst,      only: cpair
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field

! Arguments
   type(physics_state), intent(in ) :: state   ! Physics state variables
   type(physics_ptend), intent(out) :: ptend   ! individual parameterization tendencies


   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt               ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)     ! Land fraction


   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals

   real(r8), pointer :: prec(:)   ! total precipitation
   real(r8), pointer :: snow(:)   ! snow from ZM convection

   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql        ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd      ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp   ! Evaporation of deep convective precipitation

   real(r8), pointer :: pblh(:)                ! Planetary boundary layer height
   real(r8), pointer :: tpert(:)               ! Thermal temperature excess

   ! Temperature tendency from deep convection (pbuf pointer).
   real(r8), pointer, dimension(:,:) :: ttend_dp

   real(r8) zero(pcols, pver)

   integer i, k

   call pbuf_get_field(pbuf, cldtop_idx,  jctop )
   call pbuf_get_field(pbuf, cldbot_idx,  jcbot )
   call pbuf_get_field(pbuf, icwmrdp_idx, ql    )

  select case ( deep_scheme )
  case('off', 'CLUBB_SGS', 'CLUBB_MF')
    zero = 0
    mcon = 0
    cme = 0
    zdu = 0
    rliq = 0

    call physics_ptend_init(ptend, state%psetcols, 'convect_deep')

!
! Associate pointers with physics buffer fields
!

    call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1/),   kount=(/pcols,pver/) )
    call pbuf_get_field(pbuf, rprddp_idx,      rprd )
    call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
    call pbuf_get_field(pbuf, prec_dp_idx,     prec )
    call pbuf_get_field(pbuf, snow_dp_idx,     snow )

    prec=0
    snow=0

    jctop = pver
    jcbot = 1._r8
    cld = 0
    ql = 0
    rprd = 0
    ! with CLUBB-MF, fracis (insoluble fraction) must be 1.0 so that
    ! gases are transported by convtran2 in tphysac (wetdep later overwrites
    ! the aerosol entries); ql/rprd/evapcdp zeroed here are repopulated by
    ! clubb_tend_cam before any consumer reads them
    if (deep_scheme == 'CLUBB_MF') then
       fracis = 1._r8
    else
       fracis = 0
    end if
    evapcdp = 0

  case('ZM') !    1 ==> Zhang-McFarlane (default)
     call pbuf_get_field(pbuf, pblh_idx,  pblh)
     call pbuf_get_field(pbuf, tpert_idx, tpert)

     call zm_conv_tend( pblh    ,mcon    ,cme     , &
          tpert   ,zdu      , &
          rliq    , &
          ztodt   , &
          jctop, jcbot , &
          state   ,ptend   ,landfrac, pbuf)

  end select

  ! If we added temperature tendency to pbuf, set it now.

  if (ttend_dp_idx > 0) then
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
     ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
  end if

  call outfld( 'ICWMRDP ', ql  , pcols, state%lchnk )

end subroutine convect_deep_tend
!=========================================================================================


subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf)

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init

   use physics_buffer,  only: physics_buffer_desc
   use constituents, only: pcnst
   use zm_conv_intr, only: zm_conv_tend_2

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies

   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)


   if ( deep_scheme .eq. 'ZM' ) then  ! Zhang-McFarlane
      call zm_conv_tend_2( state,   ptend,  ztodt,  pbuf)
   else if ( deep_scheme .eq. 'CLUBB_MF' ) then
      call clubb_mf_convtran2( state, ptend, ztodt, pbuf )
   else
      call physics_ptend_init(ptend, state%psetcols, 'convect_deep')
   end if


end subroutine convect_deep_tend_2

!=========================================================================================

subroutine clubb_mf_convtran2( state, ptend, ztodt, pbuf)
! convective tracer transport driven by the CLUBB-MF plume
! ensemble.  Mirrors zm_conv_tend_2 (fully pbuf-driven), reading the ZM_*
! gathered arrays that clubb_tend_cam populated from the plume ensemble.

   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,   only: get_nstep
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use constituents,   only: pcnst, cnst_is_convtran2
   use ccpp_constituent_prop_mod, only: ccpp_const_props
   use zm_conv_convtran, only: zm_conv_convtran_run
   use cam_abortutils, only: endrun

! Arguments
   type(physics_state), intent(in )   :: state
   type(physics_ptend), intent(out)   :: ptend
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt

! Local variables
   integer :: i
   integer :: lengath          ! number of columns with active plumes
   integer :: nstep
   integer :: ncol

   logical :: lq(pcnst)
   real(r8), dimension(pcols,pver) :: dpdry

   ! physics buffer fields
   real(r8), pointer :: fracis(:,:,:)  ! fraction of transported species that are insoluble
   real(r8), pointer :: mu(:,:)
   real(r8), pointer :: eu(:,:)
   real(r8), pointer :: du(:,:)
   real(r8), pointer :: md(:,:)
   real(r8), pointer :: ed(:,:)
   real(r8), pointer :: dp(:,:)
   real(r8), pointer :: dsubcld(:)
   integer,  pointer :: jt(:)
   integer,  pointer :: maxg(:)
   integer,  pointer :: ideep(:)

   character(len=40)  :: scheme_name
   character(len=512) :: errmsg
   integer            :: errflg

   !-----------------------------------------------------------------------------------

   ! transport ONLY the convtran2 constituents (trace gases; aerosols are
   ! excluded via convproc_do_aer).  The convtran1 species are the PUMAS
   ! condensate/precip constituents (CLDLIQ/CLDICE/RAINQM/... registered with
   ! is_convtran1=.true.): transporting those with the plume mass fluxes
   ! double counts water already carried by the CLUBB-MF qt/thl fluxes and,
   ! without ZM's compensating heating/detrainment, destabilizes the dycore
   ! (NaN in w after ~2 weeks in testing).
   lq(1)  = .false.
   lq(2:) = cnst_is_convtran2(2:)
   call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

   call pbuf_get_field(pbuf, fracis_idx,     fracis)
   call pbuf_get_field(pbuf, zm_mu_idx,      mu)
   call pbuf_get_field(pbuf, zm_eu_idx,      eu)
   call pbuf_get_field(pbuf, zm_du_idx,      du)
   call pbuf_get_field(pbuf, zm_md_idx,      md)
   call pbuf_get_field(pbuf, zm_ed_idx,      ed)
   call pbuf_get_field(pbuf, zm_dp_idx,      dp)
   call pbuf_get_field(pbuf, zm_dsubcld_idx, dsubcld)
   call pbuf_get_field(pbuf, zm_jt_idx,      jt)
   call pbuf_get_field(pbuf, zm_maxg_idx,    maxg)
   call pbuf_get_field(pbuf, zm_ideep_idx,   ideep)

   ncol  = state%ncol
   nstep = get_nstep()

   lengath = count(ideep > 0)
   if (lengath > ncol) lengath = ncol

   if (any(ptend%lq(:)) .and. lengath > 0) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1, lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      ptend%q(:,:,:) = 0._r8

      call zm_conv_convtran_run (ncol, pver,          &
                  ptend%lq,state%q(:ncol,:,:), pcnst,  mu(:ncol,:), md(:ncol,:),   &
                  du(:ncol,:), eu(:ncol,:), ed(:ncol,:), dp(:ncol,:), dsubcld(:ncol),  &
                  jt(:ncol), maxg(:ncol), ideep(:ncol), 1, lengath,  &
                  nstep,   fracis(:ncol,:,:),  ptend%q(:ncol,:,:), dpdry(:ncol,:), ccpp_const_props, &
                  scheme_name, errmsg, errflg)

      if (errflg /= 0) then
         call endrun('clubb_mf_convtran2: From zm_conv_convtran_run: ' // errmsg)
      end if
   end if

end subroutine clubb_mf_convtran2


end module convect_deep
