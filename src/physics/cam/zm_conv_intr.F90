module zm_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair, epsilo, gravit, latvap, tmelt, rair
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv_evap,         only: zm_conv_evap_run
   use zm_convr,             only: zm_convr_init, zm_convr_run
   use zm_conv_convtran,     only: zm_conv_convtran_run
   use zm_conv_momtran,      only: zm_conv_momtran_run

   use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                               rad_cnst_get_aer_props, rad_cnst_get_mode_props !, &
   use ndrop_bam,        only: ndrop_bam_init
   use cam_abortutils,   only: endrun
   use physconst,        only: pi
   use spmd_utils,       only: masterproc
   use perf_mod
   use cam_logfile,  only: iulog
   use constituents, only: cnst_add

   implicit none
   private
   save

   ! Public methods

   public ::&
      zm_conv_register,           &! register fields in physics buffer
      zm_conv_readnl,             &! read namelist
      zm_conv_init,               &! initialize donner_deep module
      zm_conv_tend,               &! return tendencies
      zm_conv_tend_2               ! return tendencies

   public zmconv_ke, zmconv_ke_lnd,  zmconv_org  ! needed by convect_shallow

   integer ::& ! indices for fields in the physics buffer
      zm_mu_idx,      &
      zm_eu_idx,      &
      zm_du_idx,      &
      zm_md_idx,      &
      zm_ed_idx,      &
      zm_dp_idx,      &
      zm_dsubcld_idx, &
      zm_jt_idx,      &
      zm_maxg_idx,    &
      zm_ideep_idx,   &
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      ixorg,       &
      dp_cldice_idx, &
      dlfzm_idx,     &     ! detrained convective cloud water mixing ratio.
      difzm_idx,     &     ! detrained convective cloud ice mixing ratio.
      dnlfzm_idx,    &     ! detrained convective cloud water num concen.
      dnifzm_idx,    &     ! detrained convective cloud ice num concen.
      prec_dp_idx,   &
      snow_dp_idx,   &
      mconzm_idx           ! convective mass flux

   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: zmconv_c0_lnd = unset_r8
   real(r8) :: zmconv_c0_ocn = unset_r8
   real(r8) :: zmconv_ke     = unset_r8
   real(r8) :: zmconv_ke_lnd = unset_r8
   real(r8) :: zmconv_momcu  = unset_r8
   real(r8) :: zmconv_momcd  = unset_r8
   integer  :: zmconv_num_cin            ! Number of negative buoyancy regions that are allowed
                                         ! before the convection top and CAPE calculations are completed.
   logical  :: zmconv_org                ! Parameterization for sub-grid scale convective organization for the ZM deep
                                         ! convective scheme based on Mapes and Neale (2011)
   real(r8) :: zmconv_dmpdz = unset_r8        ! Parcel fractional mass entrainment rate
   real(r8) :: zmconv_tiedke_add = unset_r8   ! Convective parcel temperature perturbation
   real(r8) :: zmconv_capelmt = unset_r8      ! Triggering thereshold for ZM convection
   logical  :: zmconv_parcel_pbl = .false.             ! switch for parcel pbl calculation
   real(r8) :: zmconv_tau = unset_r8          ! Timescale for convection


!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0
   integer  ::    icwmrdp_idx      = 0
   integer  ::    rprddp_idx       = 0
   integer  ::    fracis_idx       = 0
   integer  ::    nevapr_dpcu_idx  = 0
   integer  ::    dgnum_idx        = 0

   integer :: nmodes
   integer :: nbulk

!=========================================================================================
contains
!=========================================================================================

subroutine zm_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8, dtype_i4

  implicit none

  integer idx

   call pbuf_add_field('ZM_MU', 'physpkg', dtype_r8, (/pcols,pver/), zm_mu_idx)
   call pbuf_add_field('ZM_EU', 'physpkg', dtype_r8, (/pcols,pver/), zm_eu_idx)
   call pbuf_add_field('ZM_DU', 'physpkg', dtype_r8, (/pcols,pver/), zm_du_idx)
   call pbuf_add_field('ZM_MD', 'physpkg', dtype_r8, (/pcols,pver/), zm_md_idx)
   call pbuf_add_field('ZM_ED', 'physpkg', dtype_r8, (/pcols,pver/), zm_ed_idx)

   ! wg layer thickness in mbs (between upper/lower interface).
   call pbuf_add_field('ZM_DP', 'physpkg', dtype_r8, (/pcols,pver/), zm_dp_idx)

   ! wg layer thickness in mbs between lcl and maxi.
   call pbuf_add_field('ZM_DSUBCLD', 'physpkg', dtype_r8, (/pcols/), zm_dsubcld_idx)

   ! wg top level index of deep cumulus convection.
   call pbuf_add_field('ZM_JT', 'physpkg', dtype_i4, (/pcols/), zm_jt_idx)

   ! wg gathered values of maxi.
   call pbuf_add_field('ZM_MAXG', 'physpkg', dtype_i4, (/pcols/), zm_maxg_idx)

   ! map gathered points to chunk index
   call pbuf_add_field('ZM_IDEEP', 'physpkg', dtype_i4, (/pcols/), zm_ideep_idx)

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx)

! Flux of snow from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx)

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)

   call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
   call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
   call pbuf_add_field('SNOW_DP',    'physpkg',dtype_r8,(/pcols/),     snow_dp_idx)

   ! detrained convective cloud water mixing ratio.
   call pbuf_add_field('DLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dlfzm_idx)
   ! detrained convective cloud ice mixing ratio.
   call pbuf_add_field('DIFZM', 'physpkg', dtype_r8, (/pcols,pver/), difzm_idx)
   ! convective mass fluxes
   call pbuf_add_field('CMFMC_DP', 'physpkg', dtype_r8, (/pcols,pverp/), mconzm_idx)

!CACNOTE - Is zm_org really a constituent or was it just a handy structure to use for an allocatable which persists in the run?
   if (zmconv_org) then
      call cnst_add('ZM_ORG',0._r8,0._r8,0._r8,ixorg,longname='organization parameter')
   endif

end subroutine zm_conv_register

!=========================================================================================

subroutine zm_conv_readnl(nlfile)

   use spmd_utils,      only: mpicom, masterproc, masterprocid, mpi_real8, mpi_integer, mpi_logical
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zm_conv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_num_cin, &
                        zmconv_ke, zmconv_ke_lnd, zmconv_org, &
                        zmconv_momcu, zmconv_momcd, &
                        zmconv_dmpdz, zmconv_tiedke_add, zmconv_capelmt, &
                        zmconv_parcel_pbl, zmconv_tau
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

   end if

   ! Broadcast namelist variables
   call mpi_bcast(zmconv_num_cin,           1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_num_cin")
   call mpi_bcast(zmconv_c0_lnd,            1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_c0_lnd")
   call mpi_bcast(zmconv_c0_ocn,            1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_c0_ocn")
   call mpi_bcast(zmconv_ke,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_ke")
   call mpi_bcast(zmconv_ke_lnd,            1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_ke_lnd")
   call mpi_bcast(zmconv_momcu,             1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_momcu")
   call mpi_bcast(zmconv_momcd,             1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_momcd")
   call mpi_bcast(zmconv_org,               1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_org")
   call mpi_bcast(zmconv_dmpdz,             1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_dmpdz")
   call mpi_bcast(zmconv_tiedke_add,        1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_tiedke_add")
   call mpi_bcast(zmconv_capelmt,           1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_capelmt")
   call mpi_bcast(zmconv_parcel_pbl,        1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_parcel_pbl")
   call mpi_bcast(zmconv_tau,               1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_tau")

end subroutine zm_conv_readnl

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld, add_default, horiz_only
  use ppgrid,         only: pcols, pver
  use zm_convr,       only: zm_convr_init
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces

  character(len=512) :: errmsg
  integer            :: errflg

  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

!
! Register fields with the output buffer
!

    if (zmconv_org) then
       call addfld ('ZM_ORG     ', (/ 'lev' /), 'A', '-       ','Organization parameter')
       call addfld ('ZM_ORG2D   ', (/ 'lev' /), 'A', '-       ','Organization parameter 2D')
    endif
    call addfld ('PRECZ',    horiz_only,   'A', 'm/s','total precipitation from ZM convection')
    call addfld ('ZMDT',     (/ 'lev' /),  'A', 'K/s','T tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDQ',     (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Zhang-McFarlane moist convection')
    call addfld ('ZMDICE',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud ice tendency - Zhang-McFarlane convection')
    call addfld ('ZMDLIQ',   (/ 'lev' /),  'A', 'kg/kg/s','Cloud liq tendency - Zhang-McFarlane convection')
    call addfld ('EVAPTZM',  (/ 'lev' /),  'A', 'K/s','T tendency - Evaporation/snow prod from Zhang convection')
    call addfld ('FZSNTZM',  (/ 'lev' /),  'A', 'K/s','T tendency - Rain to snow conversion from Zhang convection')
    call addfld ('EVSNTZM',  (/ 'lev' /),  'A', 'K/s','T tendency - Snow to rain prod from Zhang convection')
    call addfld ('EVAPQZM',  (/ 'lev' /),  'A', 'kg/kg/s','Q tendency - Evaporation from Zhang-McFarlane moist convection')

    call addfld ('ZMFLXPRC', (/ 'ilev' /), 'A', 'kg/m2/s','Flux of precipitation from ZM convection'       )
    call addfld ('ZMFLXSNW', (/ 'ilev' /), 'A', 'kg/m2/s','Flux of snow from ZM convection'                )
    call addfld ('ZMNTPRPD', (/ 'lev' /) , 'A', 'kg/kg/s','Net precipitation production from ZM convection')
    call addfld ('ZMNTSNPD', (/ 'lev' /) , 'A', 'kg/kg/s','Net snow production from ZM convection'         )
    call addfld ('ZMEIHEAT', (/ 'lev' /) , 'A', 'W/kg'   ,'Heating by ice and evaporation in ZM convection')

    call addfld ('CMFMC_DP', (/ 'ilev' /), 'A', 'kg/m2/s','Convection mass flux from ZM deep ')
    call addfld ('PRECCDZM', horiz_only,   'A', 'm/s','Convective precipitation rate from ZM deep')

    call addfld ('PCONVB',   horiz_only ,  'A', 'Pa'    ,'convection base pressure')
    call addfld ('PCONVT',   horiz_only ,  'A', 'Pa'    ,'convection top  pressure')

    call addfld ('CAPE',     horiz_only,   'A', 'J/kg', 'Convectively available potential energy')
    call addfld ('FREQZM',   horiz_only  , 'A', 'fraction', 'Fractional occurance of ZM convection')

    call addfld ('ZMMTT',    (/ 'lev' /),  'A', 'K/s', 'T tendency - ZM convective momentum transport')
    call addfld ('ZMMTU',    (/ 'lev' /),  'A', 'm/s2', 'U tendency - ZM convective momentum transport')
    call addfld ('ZMMTV',    (/ 'lev' /),  'A', 'm/s2', 'V tendency - ZM convective momentum transport')

    call addfld ('ZMMU',     (/ 'lev' /),  'A', 'kg/m2/s', 'ZM convection updraft mass flux')
    call addfld ('ZMMD',     (/ 'lev' /),  'A', 'kg/m2/s', 'ZM convection downdraft mass flux')

    call addfld ('ZMUPGU',   (/ 'lev' /),  'A', 'm/s2', 'zonal force from ZM updraft pressure gradient term')
    call addfld ('ZMUPGD',   (/ 'lev' /),  'A', 'm/s2', 'zonal force from ZM downdraft pressure gradient term')
    call addfld ('ZMVPGU',   (/ 'lev' /),  'A', 'm/s2', 'meridional force from ZM updraft pressure gradient term')
    call addfld ('ZMVPGD',   (/ 'lev' /),  'A', 'm/s2', 'merdional force from ZM downdraft pressure gradient term')

    call addfld ('ZMICUU',   (/ 'lev' /),  'A', 'm/s', 'ZM in-cloud U updrafts')
    call addfld ('ZMICUD',   (/ 'lev' /),  'A', 'm/s', 'ZM in-cloud U downdrafts')
    call addfld ('ZMICVU',   (/ 'lev' /),  'A', 'm/s', 'ZM in-cloud V updrafts')
    call addfld ('ZMICVD',   (/ 'lev' /),  'A', 'm/s', 'ZM in-cloud V downdrafts')

    call addfld ('DIFZM'   ,(/ 'lev' /), 'A','kg/kg/s ','Detrained ice water from ZM convection')
    call addfld ('DLFZM'   ,(/ 'lev' /), 'A','kg/kg/s ','Detrained liquid water from ZM convection')

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num)

    if (zmconv_org) then
       call add_default('ZM_ORG', 1, ' ')
       call add_default('ZM_ORG2D', 1, ' ')
    endif
    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')
       call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
    end if

!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if

    if (masterproc) then
       write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
    end if

    no_deep_pbl = phys_deepconv_pbl()
!CACNOTE - Need to check errflg and report errors
    call zm_convr_init(cpair, epsilo, gravit, latvap, tmelt, rair, &
                  limcnv,zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                  zmconv_momcu, zmconv_momcd, zmconv_num_cin, zmconv_org, &
                  no_deep_pbl, zmconv_tiedke_add, &
                  zmconv_capelmt, zmconv_dmpdz,zmconv_parcel_pbl, zmconv_tau, errmsg, errflg)

    cld_idx         = pbuf_get_index('CLD')
    fracis_idx      = pbuf_get_index('FRACIS')

end subroutine zm_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt)

subroutine zm_conv_tend(pblh    ,mcon    ,cme     , &
     tpert   ,pflx    ,zdu      , &
     rliq    ,rice    ,ztodt    , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac,  pbuf)


   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit, latice, latvap, tmelt, cpwv, cpliq, rh2o

   use phys_control,  only: cam_physpkg_is

   ! Arguments

   type(physics_state), intent(in),target   :: state          ! Physics state variables
   type(physics_ptend), intent(out)         :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer       :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: rice(pcols) ! reserved ice (not yet in cldice) for energy integrals


   ! Local variables
   character(len=512) :: errmsg
   integer            :: errflg

   integer :: i,k,l,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim_old                ! for physics buffer fields

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend),target :: ptend_loc     ! package tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice
   real(r8), pointer :: dlf(:,:)    ! detrained convective cloud water mixing ratio.
   real(r8), pointer :: dif(:,:)    ! detrained convective cloud ice mixing ratio.
   real(r8), pointer :: dnlf(:,:)   ! detrained convective cloud water num concen.
   real(r8), pointer :: dnif(:,:)   ! detrained convective cloud ice num concen.
   real(r8), pointer :: lambdadpcu(:,:) ! slope of cloud liquid size distr
   real(r8), pointer :: mudpcu(:,:)     ! width parameter of droplet size distr
   real(r8), pointer :: mconzm(:,:)     !convective mass fluxes

   real(r8), pointer :: mu(:,:)    ! (pcols,pver)
   real(r8), pointer :: eu(:,:)    ! (pcols,pver)
   real(r8), pointer :: du(:,:)    ! (pcols,pver)
   real(r8), pointer :: md(:,:)    ! (pcols,pver)
   real(r8), pointer :: ed(:,:)    ! (pcols,pver)
   real(r8), pointer :: dp(:,:)    ! (pcols,pver)
   real(r8), pointer :: dsubcld(:) ! (pcols)
   integer,  pointer :: jt(:)      ! (pcols)
   integer,  pointer :: maxg(:)    ! (pcols)
   integer,  pointer :: ideep(:)   ! (pcols)
   integer           :: lengath

   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.

   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   real(r8) :: mu_out(pcols,pver)
   real(r8) :: md_out(pcols,pver)

   ! used in momentum transport calculation
   real(r8) :: winds(pcols, pver, 2)
   real(r8) :: wind_tends(pcols, pver, 2)
   real(r8) :: pguall(pcols, pver, 2)
   real(r8) :: pgdall(pcols, pver, 2)
   real(r8) :: icwu(pcols,pver, 2)
   real(r8) :: icwd(pcols,pver, 2)
   real(r8) :: seten(pcols, pver)
   logical  :: l_windt(2)
   real(r8) :: tfinal1, tfinal2
   integer  :: ii

   real(r8),pointer :: zm_org2d(:,:)
   real(r8),allocatable :: orgt_alloc(:,:), org_alloc(:,:)

   real(r8) :: zm_org2d_ncol(state%ncol,pver)
   real(r8) :: orgt_ncol(state%ncol,pver), org_ncol(state%ncol,pver)

   logical  :: lq(pcnst)

   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)             ! copy state to local state1.

   lq(:) = .FALSE.
   lq(1) = .TRUE.
   if (zmconv_org) then
      lq(ixorg) = .TRUE.
   endif
   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr_run', ls=.true., lq=lq)! initialize local ptend type

!
! Associate pointers with physics buffer fields
!
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )

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

   call pbuf_get_field(pbuf, dlfzm_idx,  dlf)
   call pbuf_get_field(pbuf, difzm_idx,  dif)
   call pbuf_get_field(pbuf, mconzm_idx, mconzm)

   allocate(dnlf(pcols,pver), dnif(pcols,pver))

!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr_run')

   if (zmconv_org) then
      allocate(zm_org2d(pcols,pver))
      allocate(org_alloc(ncol,pver))
      allocate(orgt_alloc(ncol,pver))
      org_ncol(:ncol,:) = state%q(1:ncol,:,ixorg)
   endif

!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
   ptend_loc%q(:,:,1) = 0._r8
   ptend_loc%s(:,:) = 0._r8
   mcon(:,:) = 0._r8
   dlf(:,:) = 0._r8
   pflx(:,:) = 0._r8
   cme(:,:) = 0._r8
   cape(:) = 0._r8
   zdu(:,:) = 0._r8
   rprd(:,:) = 0._r8
   dif(:,:) = 0._r8
   dnlf(:,:) = 0._r8
   dnif(:,:) = 0._r8
   mu(:,:) = 0._r8
   eu(:,:) = 0._r8
   du(:,:) = 0._r8
   md(:,:) = 0._r8
   ed(:,:) = 0._r8
   dp(:,:) = 0._r8
   dsubcld(:) = 0._r8
   jctop(:) = 0._r8
   jcbot(:) = 0._r8
   prec(:) = 0._r8
   rliq(:) = 0._r8
   rice(:) = 0._r8
   ideep(:) = 0._r8
!REMOVECAM_END

!CACNOTE - Need to check errflg and report errors
   call zm_convr_run(ncol, pver, &
                    pverp, gravit, latice, cpwv, cpliq, rh2o,  &
                    state%t(:ncol,:), state%q(:ncol,:,1), prec(:ncol), jctop(:ncol), jcbot(:ncol), &
                    pblh(:ncol), state%zm(:ncol,:), state%phis, state%zi(:ncol,:), ptend_loc%q(:ncol,:,1), &
                    ptend_loc%s(:ncol,:), state%pmid(:ncol,:), state%pint(:ncol,:), state%pdel(:ncol,:), &
                    .5_r8*ztodt, mcon(:ncol,:), cme(:ncol,:), cape(:ncol),      &
                    tpert(:ncol), dlf(:ncol,:), pflx(:ncol,:), zdu(:ncol,:), rprd(:ncol,:), &
                    mu(:ncol,:), md(:ncol,:), du(:ncol,:), eu(:ncol,:), ed(:ncol,:),       &
                    dp(:ncol,:), dsubcld(:ncol), jt(:ncol), maxg(:ncol), ideep(:ncol),    &
                    ql(:ncol,:),  rliq(:ncol), landfrac(:ncol),                          &
                    org_ncol(:,:), orgt_ncol(:,:), zm_org2d_ncol(:,:),  &
                    dif(:ncol,:), dnlf(:ncol,:), dnif(:ncol,:),  &
                    rice(:ncol), errmsg, errflg)

   if (zmconv_org) then
      ptend_loc%q(:,:,ixorg)=orgt_ncol(:ncol,:)
      zm_org2d(:ncol,:) = zm_org2d_ncol(:ncol,:)
   endif

   lengath = count(ideep > 0)
   if (lengath > ncol) lengath = ncol  ! should not happen, but force it to not be larger than ncol for safety sake

   call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
!
! Output fractional occurance of ZM convection
!
   freqzm(:) = 0._r8
   do i = 1,lengath
      freqzm(ideep(i)) = 1.0_r8
   end do
   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   mcon(:ncol,:pverp) = mcon(:ncol,:pverp) * 100._r8/gravit
   mconzm(:ncol,:pverp) = mcon(:ncol,:pverp)

   call outfld('CMFMC_DP', mconzm, pcols, lchnk)

   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath
      do k=1,pver
         ii = ideep(i)
         mu_out(ii,k) = mu(i,k) * 100._r8/gravit
         md_out(ii,k) = md(i,k) * 100._r8/gravit
      end do
   end do

   call outfld('ZMMU', mu_out, pcols, lchnk)
   call outfld('ZMMD', md_out, pcols, lchnk)

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf ('zm_convr_run')

   call outfld('DIFZM'   ,dif            ,pcols, lchnk)
   call outfld('DLFZM'   ,dlf            ,pcols, lchnk)

   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   do i = 1,lengath
       if (maxg(i).gt.jt(i)) then
          pcont(ideep(i)) = state%pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
          pconb(ideep(i)) = state%pmid(ideep(i),maxg(i))! gathered array
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
    call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
    call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )

  call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc
  call physics_update(state1, ptend_loc, ztodt)

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  if (zmconv_org) then
     lq(ixorg) = .TRUE.
  endif
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap_run', ls=.true., lq=lq)

   call t_startf ('zm_conv_evap_run')
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8
!REMOVECAM - no longer need these when CAM is retired and pcols no longer exists
    flxprec(:,:) = 0._r8
    flxsnow(:,:) = 0._r8
    snow(:) = 0._r8
!REMOVECAM_END

    call zm_conv_evap_run(state1%ncol, pver, pverp, &
         gravit, latice, latvap, tmelt, &
         cpair, zmconv_ke, zmconv_ke_lnd, zmconv_org, &
         state1%t(:ncol,:),state1%pmid(:ncol,:),state1%pdel(:ncol,:),state1%q(:ncol,:pver,1), &
         landfrac(:ncol), &
         ptend_loc%s(:ncol,:), tend_s_snwprd(:ncol,:), tend_s_snwevmlt(:ncol,:), ptend_loc%q(:ncol,:pver,1), &
         rprd(:ncol,:), cld(:ncol,:), ztodt, &
         prec(:ncol), snow(:ncol), ntprprd(:ncol,:), ntsnprd(:ncol,:), flxprec(:ncol,:), flxsnow(:ncol,:) )

    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)

     if (zmconv_org) then
         ptend_loc%q(:ncol,:pver,ixorg) = min(1._r8,max(0._r8,(50._r8*1000._r8*1000._r8*abs(evapcdp(:ncol,:pver))) &
                                          -(state%q(:ncol,:pver,ixorg)/10800._r8)))
         ptend_loc%q(:ncol,:pver,ixorg) = (ptend_loc%q(:ncol,:pver,ixorg) - state%q(:ncol,:pver,ixorg))/ztodt
     endif

!
! Write out variables from zm_conv_evap_run
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call outfld('FZSNTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call outfld('EVSNTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
   call outfld('CMFMC_DP   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )


   call t_stopf ('zm_conv_evap_run')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc
  call physics_update(state1, ptend_loc, ztodt)


  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_momtran_run', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)

     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('zm_conv_momtran_run')

!REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
     wind_tends(:,:,:) = 0._r8
!REMOVECAM_END

     call zm_conv_momtran_run (ncol, pver, pverp,                    &
                   l_windt,winds(:ncol,:,:), 2,  mu(:ncol,:), md(:ncol,:),   &
                   zmconv_momcu, zmconv_momcd, &
                   du(:ncol,:), eu(:ncol,:), ed(:ncol,:), dp(:ncol,:), dsubcld(:ncol),  &
                   jt(:ncol), maxg(:ncol), ideep(:ncol), 1, lengath,  &
                   nstep,  wind_tends(:ncol,:,:), pguall(:ncol,:,:), pgdall(:ncol,:,:), &
                   icwu(:ncol,:,:), icwd(:ncol,:,:), ztodt, seten(:ncol,:) )
     call t_stopf ('zm_conv_momtran_run')

     ptend_loc%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
     ptend_loc%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
     ptend_loc%s(:ncol,:pver) = seten(:ncol,:pver)

     call physics_ptend_sum(ptend_loc,ptend_all, ncol)

     ! update physics state type state1 with ptend_loc
     call physics_update(state1, ptend_loc, ztodt)

     ftem(:ncol,:pver) = seten(:ncol,:pver)/cpair
     if (zmconv_org) then
        call outfld('ZM_ORG', state%q(:,:,ixorg), pcols, lchnk)
        call outfld('ZM_ORG2D', zm_org2d, pcols, lchnk)
     endif
     call outfld('ZMMTT', ftem             , pcols, lchnk)
     call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
     call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)

     ! Output apparent force from  pressure gradient
     call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
     call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
     call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
     call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)

     ! Output in-cloud winds
     call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
     call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
     call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
     call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)

   end if

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)


   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')

!REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   ptend_loc%q(:,:,:) = 0._r8
!REMOVECAM_END

   call zm_conv_convtran_run (ncol, pver,          &
                  ptend_loc%lq,state1%q(:ncol,:,:), pcnst,  mu(:ncol,:), md(:ncol,:),   &
                  du(:ncol,:), eu(:ncol,:), ed(:ncol,:), dp(:ncol,:), dsubcld(:ncol),  &
                  jt(:ncol), maxg(:ncol), ideep(:ncol), 1, lengath,  &
                  nstep,   fracis(:ncol,:,:),  ptend_loc%q(:ncol,:,:), fake_dpdry(:ncol,:), ztodt)
   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

   if (zmconv_org) then
      deallocate(zm_org2d)
   end if

   deallocate(dnlf, dnif)

end subroutine zm_conv_tend
!=========================================================================================


subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf)

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,   only: pcnst, cnst_is_convtran2

! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies

   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! Local variables
   integer :: i, lchnk, istat
   integer :: lengath          ! number of columns with deep convection
   integer :: nstep
   integer :: ncol

   real(r8), dimension(pcols,pver) :: dpdry

   ! physics buffer fields
   real(r8), pointer :: fracis(:,:,:)  ! fraction of transported species that are insoluble
   real(r8), pointer :: mu(:,:)    ! (pcols,pver)
   real(r8), pointer :: eu(:,:)    ! (pcols,pver)
   real(r8), pointer :: du(:,:)    ! (pcols,pver)
   real(r8), pointer :: md(:,:)    ! (pcols,pver)
   real(r8), pointer :: ed(:,:)    ! (pcols,pver)
   real(r8), pointer :: dp(:,:)    ! (pcols,pver)
   real(r8), pointer :: dsubcld(:) ! (pcols)
   integer,  pointer :: jt(:)      ! (pcols)
   integer,  pointer :: maxg(:)    ! (pcols)
   integer,  pointer :: ideep(:)   ! (pcols)
   !-----------------------------------------------------------------------------------


   call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=cnst_is_convtran2 )

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


   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   lengath = count(ideep > 0)
   if (lengath > ncol) lengath = ncol  ! should not happen, but force it to not be larger than ncol for safety sake

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1, lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      call t_startf ('convtran2')

!REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
   ptend%q(:,:,:) = 0._r8
!REMOVECAM_END

      call zm_conv_convtran_run (ncol, pver,          &
                  ptend%lq,state%q(:ncol,:,:), pcnst,  mu(:ncol,:), md(:ncol,:),   &
                  du(:ncol,:), eu(:ncol,:), ed(:ncol,:), dp(:ncol,:), dsubcld(:ncol),  &
                  jt(:ncol), maxg(:ncol), ideep(:ncol), 1, lengath,  &
                  nstep,   fracis(:ncol,:,:),  ptend%q(:ncol,:,:), dpdry(:ncol,:), ztodt)
      call t_stopf ('convtran2')
   end if

end subroutine zm_conv_tend_2

!=========================================================================================


end module zm_conv_intr
