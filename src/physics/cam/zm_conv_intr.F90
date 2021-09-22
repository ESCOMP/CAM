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
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap, zm_convr, convtran, momtran
   use zm_microphysics,  only: zm_aero_t, zm_conv_t
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

   public :: zmconv_microp

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
      snow_dp_idx

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
   logical  :: zmconv_microp = .false.             ! switch for microphysics


!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    
   integer  ::    dgnum_idx        = 0

   integer :: nmodes
   integer :: nbulk

   type(zm_aero_t), allocatable :: aero(:)   ! object contains information about the aerosols

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

   if (zmconv_microp) then
      ! Only add the number conc fields if the microphysics is active.

      ! detrained convective cloud water num concen.
      call pbuf_add_field('DNLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnlfzm_idx)
      ! detrained convective cloud ice num concen.
      call pbuf_add_field('DNIFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnifzm_idx)
   end if

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
                        zmconv_momcu, zmconv_momcd, zmconv_microp
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
   call mpi_bcast(zmconv_microp,            1, mpi_logical, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_readnl: FATAL: mpi_bcast: zmconv_microp")

end subroutine zm_conv_readnl

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: addfld, add_default, horiz_only
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

! Allocate the basic aero structure outside the zmconv_microp logical
! This allows the aero structure to be passed
! Note that all of the arrays inside this structure are conditionally allocated

  allocate(aero(begchunk:endchunk))

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
    
    call addfld ('CMFMCDZM', (/ 'ilev' /), 'A', 'kg/m2/s','Convection mass flux from ZM deep ')
    call addfld ('PRECCDZM', horiz_only,   'A', 'm/s','Convective precipitation rate from ZM deep')

    call addfld ('PCONVB',   horiz_only ,  'A', 'Pa'    ,'convection base pressure')
    call addfld ('PCONVT',   horiz_only ,  'A', 'Pa'    ,'convection top  pressure')

    call addfld ('CAPE',     horiz_only,   'A', 'J/kg', 'Convectively available potential energy')
    call addfld ('FREQZM',   horiz_only  , 'A', 'fraction', 'Fractional occurance of ZM convection') 

! RBN: Output variables for more detailed ZM analysis (+dynamica parcel and tau)

    call addfld ('TAUZM',   horiz_only,   'A', '/s      ', 'ZM deep convection timescale')  
    call addfld ('WINCLD', (/ 'lev' /),   'A', 'm/s     ', 'Deep convective in-cloud vertical velocity')
    call addfld ('KEPAR',  (/ 'lev' /),   'A', 'J/kg    ', 'Convective parcel kinetic energy')
    call addfld ('BUOY',   (/ 'lev' /),   'A', 'K       ', 'Buoyancy as temperature')

    call addfld ('MWINCLD', horiz_only,'A','m/s     ',    'Deep convective mean in-cloud vertical velocity')
    call addfld ('HMAX',    horiz_only,'A','unitless',    'Moist Static energy maximum')
    call addfld ('LCL',     horiz_only,'A','unitless',    'Lifting condensation model level index')
    call addfld ('LEL',     horiz_only,'A','unitless',    'Convective top negative buoyancy level index')
    call addfld ('KHMAX',   horiz_only,'A','unitless',    'Moist Static energy maximum level index') 
    call addfld ('TLCL',    horiz_only,'A','K       ',    'Temperature at the lifting condensation level')
    call addfld ('PLCL',    horiz_only,'A','K       ',    'Pressure at the lifting condensation level')
    
!!!!
    
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

    if (zmconv_microp) then
       call add_default ('DIFZM',    1, ' ')
       call add_default ('DLFZM',    1, ' ')
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
    call zm_convi(limcnv,zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                  zmconv_momcu, zmconv_momcd, zmconv_num_cin, zmconv_org, &
                  zmconv_microp, no_deep_pbl_in = no_deep_pbl)

    cld_idx         = pbuf_get_index('CLD')
    fracis_idx      = pbuf_get_index('FRACIS')

    if (zmconv_microp) call zm_conv_micro_init()

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

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit
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

   type(zm_conv_t)              :: conv

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
   real(r8),pointer :: orgt(:,:), org(:,:)

   logical  :: lq(pcnst)

   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   if (zmconv_microp) then
     allocate( &
       conv%qi(pcols,pver), &
       conv%qliq(pcols,pver), &
       conv%qice(pcols,pver), &
       conv%wu(pcols,pver), &
       conv%sprd(pcols,pver), &
       conv%qrain(pcols,pver), &
       conv%qsnow(pcols,pver), &
       conv%qnl(pcols,pver), &
       conv%qni(pcols,pver), &
       conv%qnr(pcols,pver), &
       conv%qns(pcols,pver), &
       conv%frz(pcols,pver), &
       conv%autolm(pcols,pver), &
       conv%accrlm(pcols,pver), &
       conv%bergnm(pcols,pver), &
       conv%fhtimm(pcols,pver), &
       conv%fhtctm(pcols,pver), &
       conv%fhmlm (pcols,pver), &
       conv%hmpim (pcols,pver), &
       conv%accslm(pcols,pver), &
       conv%dlfm  (pcols,pver), &
       conv%autoln(pcols,pver), &
       conv%accrln(pcols,pver), &
       conv%bergnn(pcols,pver), &
       conv%fhtimn(pcols,pver), &
       conv%fhtctn(pcols,pver), &
       conv%fhmln (pcols,pver), &
       conv%accsln(pcols,pver), &
       conv%activn(pcols,pver), &
       conv%dlfn  (pcols,pver), &
       conv%autoim(pcols,pver), &
       conv%accsim(pcols,pver), &
       conv%difm  (pcols,pver), &
       conv%nuclin(pcols,pver), &
       conv%autoin(pcols,pver), &
       conv%accsin(pcols,pver), &
       conv%hmpin (pcols,pver), &
       conv%difn  (pcols,pver), &
       conv%cmel  (pcols,pver), &
       conv%cmei  (pcols,pver), &
       conv%trspcm(pcols,pver), &
       conv%trspcn(pcols,pver), &
       conv%trspim(pcols,pver), &
       conv%trspin(pcols,pver), &
       conv%lambdadpcu(pcols,pver), &
       conv%mudpcu(pcols,pver), &
       conv%dcape(pcols) )
   end if

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
   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq)! initialize local ptend type

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

   if (zmconv_microp) then
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlf)
      call pbuf_get_field(pbuf, dnifzm_idx, dnif)
   else
       allocate(dnlf(pcols,pver), dnif(pcols,pver))
   end if

   if (zmconv_microp) then

      if (nmodes > 0) then

         ! Associate pointers with the modes and species that affect the climate
         ! (list 0)

         do m = 1, nmodes
            call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, aero(lchnk)%num_a(m)%val)
            call pbuf_get_field(pbuf, dgnum_idx, aero(lchnk)%dgnum(m)%val, start=(/1,1,m/), kount=(/pcols,pver,1/))

            do l = 1, aero(lchnk)%nspec(m)
               call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, aero(lchnk)%mmr_a(l,m)%val)
            end do
         end do

      else if (nbulk > 0) then

         ! Associate pointers with the bulk aerosols that affect the climate
         ! (list 0)

         do m = 1, nbulk
            call rad_cnst_get_aer_mmr(0, m, state, pbuf, aero(lchnk)%mmr_bulk(m)%val)
         end do

      end if
   end if

!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')

   if (zmconv_org) then
      allocate(zm_org2d(pcols,pver))
      org => state%q(:,:,ixorg) 
      orgt => ptend_loc%q(:,:,ixorg)
   endif

   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q(:,:,1),    state%omega,  prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    , state%pmid     ,state%pint    ,state%pdel     , &
                    .5_r8*ztodt    ,mcon    ,cme     , cape,      &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu,      md,      du,      eu,      ed,       &
                    dp,      dsubcld, jt,      maxg,    ideep,    &
                    ql,  rliq, landfrac,                          &
                    org, orgt, zm_org2d,  &
                    dif, dnlf, dnif,  conv, &
                    aero(lchnk), rice)

   lengath = count(ideep > 0)

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
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit

   call outfld('CMFMCDZM', mcon, pcols, lchnk)

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
   call t_stopf ('zm_convr')

   call outfld('DIFZM'   ,dif            ,pcols, lchnk)
   call outfld('DLFZM'   ,dlf            ,pcols, lchnk)

   if (zmconv_microp) call zm_conv_micro_outfld(conv, dnif, dnlf, lchnk, ncol)

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
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

   call t_startf ('zm_conv_evap')
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

    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         landfrac, &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow, conv%sprd)

    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)
    
     if (zmconv_org) then
         ptend_loc%q(:ncol,:pver,ixorg) = min(1._r8,max(0._r8,(50._r8*1000._r8*1000._r8*abs(evapcdp(:ncol,:pver))) &
                                          -(state%q(:ncol,:pver,ixorg)/10800._r8)))
         ptend_loc%q(:ncol,:pver,ixorg) = (ptend_loc%q(:ncol,:pver,ixorg) - state%q(:ncol,:pver,ixorg))/ztodt 
     endif    
    
!
! Write out variables from zm_conv_evap
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
   call outfld('CMFMCDZM   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )


   call t_stopf ('zm_conv_evap')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)


  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
   
     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('momtran')
     call momtran (lchnk, ncol,                                        &
                   l_windt,winds, 2,  mu, md,   &
                   du, eu, ed, dp, dsubcld,  &
                   jt, maxg, ideep, 1, lengath,  &
                   nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
     call t_stopf ('momtran')

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
   call convtran (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu, md,   &
                  du, eu, ed, dp, dsubcld,  &
                  jt,maxg, ideep, 1, lengath,  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry, ztodt)
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

   if (zmconv_microp) then
     deallocate( &
       conv%qi, &
       conv%qliq, &
       conv%qice, &
       conv%wu, &
       conv%sprd, &
       conv%qrain, &
       conv%qsnow, &
       conv%qnl, &
       conv%qni, &
       conv%qnr, &
       conv%qns, &
       conv%frz, &
       conv%autolm, &
       conv%accrlm, &
       conv%bergnm, &
       conv%fhtimm, &
       conv%fhtctm, &
       conv%fhmlm , &
       conv%hmpim , &
       conv%accslm, &
       conv%dlfm  , &
       conv%autoln, &
       conv%accrln, &
       conv%bergnn, &
       conv%fhtimn, &
       conv%fhtctn, &
       conv%fhmln , &
       conv%accsln, &
       conv%activn, &
       conv%dlfn  , &
       conv%autoim, &
       conv%accsim, &
       conv%difm  , &
       conv%nuclin, &
       conv%autoin, &
       conv%accsin, &
       conv%hmpin , &
       conv%difn  , &
       conv%cmel  , &
       conv%cmei  , &
       conv%trspcm, &
       conv%trspcn, &
       conv%trspim, &
       conv%trspin, &
       conv%lambdadpcu, &
       conv%mudpcu, &
       conv%dcape )

   else

      deallocate(dnlf, dnif)

   end if

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

   lengath = count(ideep > 0)

   lchnk = state%lchnk
   nstep = get_nstep()

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1, lengath
         dpdry(i,:) = state%pdeldry(ideep(i),:)/100._r8
      end do

      call t_startf ('convtran2')
      call convtran (lchnk,                                        &
                     ptend%lq,state%q, pcnst,  mu, md,   &
                     du, eu, ed, dp, dsubcld,  &
                     jt, maxg, ideep, 1, lengath,  &
                     nstep,   fracis,  ptend%q, dpdry, ztodt)
      call t_stopf ('convtran2')
   end if

end subroutine zm_conv_tend_2

!=========================================================================================

subroutine zm_conv_micro_init()

  use cam_history,    only: addfld, add_default, horiz_only
  use ppgrid,         only: pcols, pver
  use pmgrid,         only: plev,plevp
  use phys_control,   only: cam_physpkg_is
  use physics_buffer, only: pbuf_get_index
  use zm_microphysics, only: zm_mphyi

  implicit none

  integer :: i

  ! 
  ! Register fields with the output buffer
  !
    call addfld ('ICIMRDP', (/ 'lev' /), 'A','kg/kg',  'Deep Convection in-cloud ice mixing ratio ')
    call addfld ('CLDLIQZM',(/ 'lev' /), 'A','g/m3'    ,'Cloud liquid water - ZM convection')
    call addfld ('CLDICEZM',(/ 'lev' /), 'A','g/m3'    ,'Cloud ice water - ZM convection')
    call addfld ('CLIQSNUM',(/ 'lev' /), 'A','1'       ,'Cloud liquid water sample number - ZM convection')
    call addfld ('CICESNUM',(/ 'lev' /), 'A','1'       ,'Cloud ice water sample number - ZM convection')
    call addfld ('QRAINZM' ,(/ 'lev' /), 'A','g/m3'    ,'rain water - ZM convection')
    call addfld ('QSNOWZM' ,(/ 'lev' /), 'A','g/m3'    ,'snow - ZM convection')
    call addfld ('CRAINNUM',(/ 'lev' /), 'A','1'       ,'Cloud rain water sample number - ZM convection')
    call addfld ('CSNOWNUM',(/ 'lev' /), 'A','1'       ,'Cloud snow sample number - ZM convection')

    call addfld ('DNIFZM'  ,(/ 'lev' /), 'A','1/kg/s ' ,'Detrained ice water num concen from ZM convection')
    call addfld ('DNLFZM'  ,(/ 'lev' /), 'A','1/kg/s ' ,'Detrained liquid water num concen from ZM convection')
    call addfld ('WUZM'    ,(/ 'lev' /), 'A','m/s'     ,'vertical velocity - ZM convection')
    call addfld ('WUZMSNUM',(/ 'lev' /), 'A','1'       ,'vertical velocity sample number - ZM convection')

    call addfld ('QNLZM',(/ 'lev' /), 'A','1/m3'       ,'Cloud liquid water number concen - ZM convection')
    call addfld ('QNIZM',(/ 'lev' /), 'A','1/m3'       ,'Cloud ice number concen - ZM convection')
    call addfld ('QNRZM',(/ 'lev' /), 'A','1/m3'       ,'Cloud rain water number concen - ZM convection')
    call addfld ('QNSZM',(/ 'lev' /), 'A','1/m3'       ,'Cloud snow number concen - ZM convection')

    call addfld ('FRZZM',(/ 'lev' /), 'A','1/s'       ,'mass tendency due to freezing - ZM convection')

    call addfld ('AUTOL_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to autoconversion of droplets to rain')
    call addfld ('ACCRL_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to accretion of droplets by rain')
    call addfld ('BERGN_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to Bergeron process')
    call addfld ('FHTIM_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to immersion freezing')
    call addfld ('FHTCT_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to contact freezing')
    call addfld ('FHML_M'  ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to homogeneous freezing of droplet')
    call addfld ('HMPI_M'  ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to HM process')
    call addfld ('ACCSL_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to accretion of droplet by snow')
    call addfld ('DLF_M'   ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to detrainment of droplet')
    call addfld ('COND_M'  ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to condensation')

    call addfld ('AUTOL_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to autoconversion of droplets to rain')
    call addfld ('ACCRL_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to accretion of droplets by rain')
    call addfld ('BERGN_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to Bergeron process')
    call addfld ('FHTIM_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to immersion freezing')
    call addfld ('FHTCT_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to contact freezing')
    call addfld ('FHML_N'  ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to homogeneous freezing of droplet')
    call addfld ('ACCSL_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to accretion of droplet by snow')
    call addfld ('ACTIV_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to droplets activation')
    call addfld ('DLF_N'   ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to detrainment of droplet')

    call addfld ('AUTOI_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to autoconversion of ice to snow')
    call addfld ('ACCSI_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to accretion of ice by snow')
    call addfld ('DIF_M'   ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to detrainment of cloud ice')
    call addfld ('DEPOS_M' ,(/ 'lev' /), 'A','kg/kg/m' ,'mass tendency due to deposition')

    call addfld ('NUCLI_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to ice nucleation')
    call addfld ('AUTOI_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to autoconversion of ice to snow')
    call addfld ('ACCSI_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to accretion of ice by snow')
    call addfld ('HMPI_N'  ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to HM process')
    call addfld ('DIF_N'   ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency due to detrainment of cloud ice')

    call addfld ('TRSPC_M' ,(/ 'lev' /), 'A','kg/kg/m','mass tendency of droplets due to convective transport')
    call addfld ('TRSPC_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency of droplets due to convective transport')
    call addfld ('TRSPI_M' ,(/ 'lev' /), 'A','kg/kg/m','mass tendency of ice crystal due to convective transport')
    call addfld ('TRSPI_N' ,(/ 'lev' /), 'A','1/kg/m' ,'num tendency of ice crystal due to convective transport')


    call add_default ('CLDLIQZM', 1, ' ')
    call add_default ('CLDICEZM', 1, ' ')
    call add_default ('CLIQSNUM', 1, ' ')
    call add_default ('CICESNUM', 1, ' ')
    call add_default ('DNIFZM',   1, ' ')
    call add_default ('DNLFZM',   1, ' ')
    call add_default ('WUZM',     1, ' ')
    call add_default ('QRAINZM',  1, ' ')
    call add_default ('QSNOWZM',  1, ' ')
    call add_default ('CRAINNUM', 1, ' ')
    call add_default ('CSNOWNUM', 1, ' ')
    call add_default ('QNLZM',    1, ' ')
    call add_default ('QNIZM',    1, ' ')
    call add_default ('QNRZM',    1, ' ')
    call add_default ('QNSZM',    1, ' ')
    call add_default ('FRZZM',   1, ' ')

    ! Initialization for the microphysics

    call zm_mphyi()

    ! Initialize the aerosol object with data from the modes/species
    ! affecting climate,
    ! i.e., the list index is hardcoded to 0.

    call rad_cnst_get_info(0, nmodes=nmodes, naero=nbulk)


    do i = begchunk, endchunk
       call zm_aero_init(nmodes, nbulk, aero(i))
    end do

    if (nmodes > 0) then

       dgnum_idx = pbuf_get_index('DGNUM')

    else if (nbulk > 0 .and.  cam_physpkg_is('cam4')) then

       ! This call is needed to allow running the ZM microphysics with the
       ! cam4 physics package.
       call ndrop_bam_init()

    end if

 end subroutine zm_conv_micro_init


 subroutine zm_aero_init(nmodes, nbulk, aero)

  use pmgrid,         only: plev,plevp

    ! Initialize the zm_aero_t object for modal aerosols

    integer,         intent(in)  :: nmodes
    integer,         intent(in)  :: nbulk
    type(zm_aero_t), intent(out) :: aero

    integer :: iaer, l, m
    integer :: nspecmx   ! max number of species in a mode

    character(len=20), allocatable :: aername(:)
    character(len=32) :: str32
    character(len=*), parameter :: routine = 'zm_conv_init'

    real(r8) :: sigmag, dgnumlo, dgnumhi
    real(r8) :: alnsg
    !----------------------------------------------------------------------------------

    aero%nmodes = nmodes
    aero%nbulk  = nbulk

       if (nmodes > 0) then

          ! Initialize the modal aerosol information

          aero%scheme = 'modal'

          ! Get number of species in each mode, and find max.
          allocate(aero%nspec(aero%nmodes))
          nspecmx = 0
          do m = 1, aero%nmodes

             call rad_cnst_get_info(0, m, nspec=aero%nspec(m), mode_type=str32)

             nspecmx = max(nspecmx, aero%nspec(m))

             ! save mode index for specified mode types
             select case (trim(str32))
             case ('accum')
                aero%mode_accum_idx = m
             case ('aitken')
                aero%mode_aitken_idx = m
             case ('coarse')
                aero%mode_coarse_idx = m
             end select

          end do

          ! Check that required mode types were found
          if (aero%mode_accum_idx == -1 .or. aero%mode_aitken_idx == -1 .or. aero%mode_coarse_idx == -1) then
             write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
                aero%mode_accum_idx, aero%mode_aitken_idx, aero%mode_coarse_idx
             call endrun(routine//': ERROR required mode type not found')
          end if

          ! find indices for the dust and seasalt species in the coarse mode
          do l = 1, aero%nspec(aero%mode_coarse_idx)
             call rad_cnst_get_info(0, aero%mode_coarse_idx, l, spec_type=str32)
             select case (trim(str32))
             case ('dust')
                aero%coarse_dust_idx = l
             case ('seasalt')
                aero%coarse_nacl_idx = l
             end select
          end do
          ! Check that required modal specie types were found
          if (aero%coarse_dust_idx == -1 .or. aero%coarse_nacl_idx == -1) then
             write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
                aero%coarse_dust_idx, aero%coarse_nacl_idx
             call endrun(routine//': ERROR required mode-species type not found')
          end if

          allocate( &
             aero%num_a(nmodes), &
             aero%mmr_a(nspecmx,nmodes), &
             aero%numg_a(pcols,pver,nmodes), &
             aero%mmrg_a(pcols,pver,nspecmx,nmodes), &
             aero%voltonumblo(nmodes), &
             aero%voltonumbhi(nmodes), &
             aero%specdens(nspecmx,nmodes), &
             aero%spechygro(nspecmx,nmodes), &
             aero%dgnum(nmodes), &
             aero%dgnumg(pcols,pver,nmodes) )


          do m = 1, nmodes

             ! Properties of modes
             call rad_cnst_get_mode_props(0, m, &
                sigmag=sigmag, dgnumlo=dgnumlo, dgnumhi=dgnumhi)

             alnsg               = log(sigmag)
             aero%voltonumblo(m) = 1._r8 / ( (pi/6._r8)*(dgnumlo**3._r8)*exp(4.5_r8*alnsg**2._r8) )
             aero%voltonumbhi(m) = 1._r8 / ( (pi/6._r8)*(dgnumhi**3._r8)*exp(4.5_r8*alnsg**2._r8) )

             ! save sigmag of aitken mode
             if (m == aero%mode_aitken_idx) aero%sigmag_aitken = sigmag

             ! Properties of modal species
             do l = 1, aero%nspec(m)
                call rad_cnst_get_aer_props(0, m, l, density_aer=aero%specdens(l,m), &
                   hygro_aer=aero%spechygro(l,m))
             end do
          end do

       else if (nbulk > 0) then

          aero%scheme = 'bulk'

          ! Props needed for BAM number concentration calcs.
          allocate( &
             aername(nbulk),                   &
             aero%num_to_mass_aer(nbulk),      &
             aero%mmr_bulk(nbulk),             &
             aero%mmrg_bulk(pcols,plev,nbulk)  )

          do iaer = 1, aero%nbulk
             call rad_cnst_get_aer_props(0, iaer, &
                aername         = aername(iaer), &
                num_to_mass_aer = aero%num_to_mass_aer(iaer) )

             ! Look for sulfate aerosol in this list (Bulk aerosol only)
             if (trim(aername(iaer)) == 'SULFATE') aero%idxsul = iaer
             if (trim(aername(iaer)) == 'DUST1')   aero%idxdst1 = iaer
             if (trim(aername(iaer)) == 'DUST2')   aero%idxdst2 = iaer
             if (trim(aername(iaer)) == 'DUST3')   aero%idxdst3 = iaer
             if (trim(aername(iaer)) == 'DUST4')   aero%idxdst4 = iaer
             if (trim(aername(iaer)) == 'BCPHI')   aero%idxbcphi = iaer
          end do

       end if

    end subroutine zm_aero_init

   subroutine zm_conv_micro_outfld(conv, dnif, dnlf, lchnk, ncol)

   use cam_history,   only: outfld

   type(zm_conv_t),intent(in)  :: conv
   real(r8), intent(in) :: dnlf(:,:)   ! detrained convective cloud water num concen.
   real(r8), intent(in) :: dnif(:,:)   ! detrained convective cloud ice num concen.
   integer, intent(in)         :: lchnk
   integer, intent(in)         :: ncol

   integer :: i,k

   real(r8) :: cice_snum(pcols,pver)      ! convective cloud ice sample number.
   real(r8) :: cliq_snum(pcols,pver)      ! convective cloud liquid sample number.
   real(r8) :: crain_snum(pcols,pver)     ! convective rain water sample number.
   real(r8) :: csnow_snum(pcols,pver)     ! convective snow sample number.
   real(r8) :: wu_snum(pcols,pver)        ! vertical velocity sample number

   real(r8) :: qni_snum(pcols,pver)       ! convective cloud ice number sample number.
   real(r8) :: qnl_snum(pcols,pver)       ! convective cloud liquid number sample number.

       do k = 1,pver
          do i = 1,ncol
             if (conv%qice(i,k) .gt. 0.0_r8) then
                cice_snum(i,k) = 1.0_r8
             else
                cice_snum(i,k) = 0.0_r8
             end if
             if (conv%qliq(i,k) .gt. 0.0_r8) then
                cliq_snum(i,k) = 1.0_r8
             else
                cliq_snum(i,k) = 0.0_r8
             end if
             if (conv%qsnow(i,k) .gt. 0.0_r8) then
                csnow_snum(i,k) = 1.0_r8
             else
                csnow_snum(i,k) = 0.0_r8
             end if
             if (conv%qrain(i,k) .gt. 0.0_r8) then
                crain_snum(i,k) = 1.0_r8
             else
                crain_snum(i,k) = 0.0_r8
             end if

             if (conv%qnl(i,k) .gt. 0.0_r8) then
                qnl_snum(i,k) = 1.0_r8
             else
                qnl_snum(i,k) = 0.0_r8
             end if
             if (conv%qni(i,k) .gt. 0.0_r8) then
                qni_snum(i,k) = 1.0_r8
             else
                qni_snum(i,k) = 0.0_r8
             end if
             if (conv%wu(i,k) .gt. 0.0_r8) then
                wu_snum(i,k) = 1.0_r8
             else
                wu_snum(i,k) = 0.0_r8
             end if

          end do
       end do

       call outfld('ICIMRDP ',conv%qi        ,pcols, lchnk )
       call outfld('CLDLIQZM',conv%qliq      ,pcols, lchnk)
       call outfld('CLDICEZM',conv%qice      ,pcols, lchnk)
       call outfld('CLIQSNUM',cliq_snum      ,pcols, lchnk)
       call outfld('CICESNUM',cice_snum      ,pcols, lchnk)
       call outfld('QRAINZM' ,conv%qrain     ,pcols, lchnk)
       call outfld('QSNOWZM' ,conv%qsnow     ,pcols, lchnk)
       call outfld('CRAINNUM',crain_snum     ,pcols, lchnk)
       call outfld('CSNOWNUM',csnow_snum     ,pcols, lchnk)

       call outfld('WUZM'    ,conv%wu        ,pcols, lchnk)
       call outfld('WUZMSNUM',wu_snum        ,pcols, lchnk)
       call outfld('QNLZM'   ,conv%qnl       ,pcols, lchnk)
       call outfld('QNIZM'   ,conv%qni       ,pcols, lchnk)
       call outfld('QNRZM'   ,conv%qnr       ,pcols, lchnk)
       call outfld('QNSZM'   ,conv%qns       ,pcols, lchnk)
       call outfld('FRZZM'   ,conv%frz       ,pcols, lchnk)

       call outfld('AUTOL_M' ,conv%autolm    ,pcols, lchnk)
       call outfld('ACCRL_M' ,conv%accrlm    ,pcols, lchnk)
       call outfld('BERGN_M' ,conv%bergnm    ,pcols, lchnk)
       call outfld('FHTIM_M' ,conv%fhtimm    ,pcols, lchnk)
       call outfld('FHTCT_M' ,conv%fhtctm    ,pcols, lchnk)
       call outfld('FHML_M'  ,conv%fhmlm     ,pcols, lchnk)
       call outfld('HMPI_M'  ,conv%hmpim     ,pcols, lchnk)
       call outfld('ACCSL_M' ,conv%accslm    ,pcols, lchnk)
       call outfld('DLF_M'   ,conv%dlfm      ,pcols, lchnk)

       call outfld('AUTOL_N' ,conv%autoln    ,pcols, lchnk)
       call outfld('ACCRL_N' ,conv%accrln    ,pcols, lchnk)
       call outfld('BERGN_N' ,conv%bergnn    ,pcols, lchnk)
       call outfld('FHTIM_N' ,conv%fhtimn    ,pcols, lchnk)
       call outfld('FHTCT_N' ,conv%fhtctn    ,pcols, lchnk)
       call outfld('FHML_N'  ,conv%fhmln     ,pcols, lchnk)
       call outfld('ACCSL_N' ,conv%accsln    ,pcols, lchnk)
       call outfld('ACTIV_N' ,conv%activn    ,pcols, lchnk)
       call outfld('DLF_N'   ,conv%dlfn      ,pcols, lchnk)
       call outfld('AUTOI_M' ,conv%autoim    ,pcols, lchnk)
       call outfld('ACCSI_M' ,conv%accsim    ,pcols, lchnk)
       call outfld('DIF_M'   ,conv%difm      ,pcols, lchnk)
       call outfld('NUCLI_N' ,conv%nuclin    ,pcols, lchnk)
       call outfld('AUTOI_N' ,conv%autoin    ,pcols, lchnk)
       call outfld('ACCSI_N' ,conv%accsin    ,pcols, lchnk)
       call outfld('HMPI_N'  ,conv%hmpin     ,pcols, lchnk)
       call outfld('DIF_N'   ,conv%difn      ,pcols, lchnk)
       call outfld('COND_M'  ,conv%cmel      ,pcols, lchnk)
       call outfld('DEPOS_M' ,conv%cmei      ,pcols, lchnk)

       call outfld('TRSPC_M' ,conv%trspcm    ,pcols, lchnk)
       call outfld('TRSPC_N' ,conv%trspcn    ,pcols, lchnk)
       call outfld('TRSPI_M' ,conv%trspim    ,pcols, lchnk)
       call outfld('TRSPI_N' ,conv%trspin    ,pcols, lchnk)
       call outfld('DNIFZM'  ,dnif           ,pcols, lchnk)
       call outfld('DNLFZM'  ,dnlf           ,pcols, lchnk)

 end subroutine zm_conv_micro_outfld

end module zm_conv_intr
