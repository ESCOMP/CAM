module spcam_drivers


use camsrfexch,       only: cam_out_t, cam_in_t
use ppgrid,           only: pcols, pver
use camsrfexch ,      only: cam_export
use shr_kind_mod,     only: r8 => shr_kind_r8
#ifdef CRM
use crmdims,          only: crm_nx, crm_ny, crm_nz
#endif
use radiation,        only: rad_out_t
use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
use physics_types,    only: physics_state, physics_state_copy, physics_ptend
use pkg_cldoptics,    only: cldems, cldovrlap, cldefr
use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
use cam_history,      only: outfld
use cam_history_support, only : fillvalue

implicit none
save
private

type rad_avgdata_type_sam1mom
   real(r8), allocatable :: solin_m(:)         ! Solar incident flux
   real(r8), allocatable :: fsntoa_m(:)        ! Net solar flux at TOA
   real(r8), allocatable :: fsutoa_m(:)        ! upwelling solar flux at TOA
   real(r8), allocatable :: fsntoac_m(:)       ! Clear sky net solar flux at TOA
   real(r8), allocatable :: fsnirt_m(:)        ! Near-IR flux absorbed at toa
   real(r8), allocatable :: fsnrtc_m(:)        ! Clear sky near-IR flux absorbed at toa
   real(r8), allocatable :: fsnirtsq_m(:)      ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8), allocatable :: fsntc_m(:)         ! Clear sky total column abs solar flux
   real(r8), allocatable :: fsnsc_m(:)         ! Clear sky surface abs solar flux
   real(r8), allocatable :: fsdsc_m(:)         ! Clear sky surface downwelling solar flux
   real(r8), allocatable :: flut_m(:)          ! Upward flux at top of model
   real(r8), allocatable :: flutc_m(:)         ! Upward Clear Sky flux at top of model
   real(r8), allocatable :: flntc_m(:)         ! Clear sky lw flux at model top
   real(r8), allocatable :: flnsc_m(:)         ! Clear sky lw flux at srf (up-down)
   real(r8), allocatable :: fldsc_m(:)         ! Clear sky lw flux at srf (down)
   real(r8), allocatable :: flwds_m(:)          ! Down longwave flux at surface
   real(r8), allocatable :: fsns_m(:)          ! Surface solar absorbed flux
   real(r8), allocatable :: fsnr_m(:)
   real(r8), allocatable :: fsnt_m(:)          ! Net column abs solar flux at model top
   real(r8), allocatable :: flns_m(:)          ! Srf longwave cooling (up-down) flux
   real(r8), allocatable :: flnt_m(:)          ! Net outgoing lw flux at model top
   real(r8), allocatable :: flnr_m(:)
   real(r8), allocatable :: fsds_m(:)          ! Surface solar down flux
   real(r8), allocatable :: fln200_m(:)        ! net longwave flux interpolated to 200 mb
   real(r8), allocatable :: fln200c_m(:)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8), allocatable :: fsn200_m(:)        ! fns interpolated to 200 mb
   real(r8), allocatable :: fsn200c_m(:)       ! fcns interpolated to 200 mb
   real(r8), allocatable :: sols_m(:)          ! Solar downward visible direct  to surface
   real(r8), allocatable :: soll_m(:)          ! Solar downward near infrared direct  to surface
   real(r8), allocatable :: solsd_m(:)         ! Solar downward visible diffuse to surface
   real(r8), allocatable :: solld_m(:)         ! Solar downward near infrared diffuse to surface
   real(r8), allocatable :: qrs_m(:,:)
   real(r8), allocatable :: qrl_m(:,:)
   real(r8), allocatable :: qrsc_m(:,:)
   real(r8), allocatable :: qrlc_m(:,:)
   real(r8), allocatable :: rel_crm(:,:,:,:)
   real(r8), allocatable :: rei_crm(:,:,:,:)
   real(r8), allocatable :: qrl_crm(:,:,:,:)
   real(r8), allocatable :: qrs_crm(:,:,:,:)
   real(r8), allocatable :: fsdtoa_m(:)        ! Solar input = Flux Solar Downward Top of Atmosphere
   real(r8), allocatable :: flds_m(:)          ! Down longwave flux at surface

   real(r8), pointer :: t_rad (:,:,:,:) ! rad temperuture
   real(r8), pointer :: qv_rad(:,:,:,:) ! rad vapor
   real(r8), pointer :: qc_rad(:,:,:,:) ! rad cloud water
   real(r8), pointer :: qi_rad(:,:,:,:) ! rad cloud ice
   real(r8), pointer :: crm_qrad(:,:,:,:) ! rad heating

   real(r8), allocatable :: tot_cld_vistau_m(:,:)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: tot_icld_vistau_m(:,:)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: liq_icld_vistau_m(:,:)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: ice_icld_vistau_m(:,:)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: nct_tot_icld_vistau_m(:,:) ! the number of CRM columns that has in-cloud visible sw optical depth
   real(r8), allocatable :: nct_liq_icld_vistau_m(:,:) ! the number of CRM column that has liq in-cloud visible sw optical depth
   real(r8), allocatable :: nct_ice_icld_vistau_m(:,:) ! the number of CRM column that has ice in-cloud visible sw optical depth

   ! Just used in m2005 -- needed for compilation only
   real(r8), allocatable :: snow_icld_vistau_m(:,:) ! snow in-cloud visible sw optical depth for output on history files
   real(r8), allocatable :: nct_snow_icld_vistau_m(:,:) ! the number of CRM column that has snow in-cloud visible sw optical depth
   real(r8), allocatable ::  crm_aodvisz(:,:,:,:)   ! layer aerosol optical depth at 550nm at CRM grids
   real(r8), allocatable ::  crm_aodvis(:,:,:)   ! AOD at 550nm at CRM grids
   real(r8), allocatable ::  crm_aod400(:,:,:)   ! AOD at 400nm at CRM grids
   real(r8), allocatable ::  crm_aod700(:,:,:)   ! AOD at 700nm at CRM grids
   real(r8), allocatable ::  aod400(:)   ! AOD at 400nm at CRM grids
   real(r8), allocatable ::  aod700(:)   ! AOD at 700nm at CRM grids
   real(r8), allocatable :: cld_tau_crm(:,:,:,:)
   real(r8), allocatable :: crm_fsnt(:,:,:)   ! net shortwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_fsntc(:,:,:)   ! net clear-sky shortwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_fsns(:,:,:)   ! net shortwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_fsnsc(:,:,:)   ! net clear-sky shortwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_flnt(:,:,:)   ! net longwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_flntc(:,:,:)   ! net clear-sky longwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_flns(:,:,:)   ! net longwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_flnsc(:,:,:)   ! net clear-sky longwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_swcf(:,:,:)   ! shortwave cloud forcing at CRM grids
end type rad_avgdata_type_sam1mom

type rad_avgdata_type_m2005
   real(r8),allocatable :: solin_m(:)         ! Solar incident flux
   real(r8),allocatable :: fsntoa_m(:)        ! Net solar flux at TOA
   real(r8),allocatable :: fsutoa_m(:)        ! upwelling solar flux at TOA
   real(r8),allocatable :: fsntoac_m(:)       ! Clear sky net solar flux at TOA
   real(r8),allocatable :: fsnirt_m(:)        ! Near-IR flux absorbed at toa
   real(r8),allocatable :: fsnrtc_m(:)        ! Clear sky near-IR flux absorbed at toa
   real(r8),allocatable :: fsnirtsq_m(:)      ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8),allocatable :: fsntc_m(:)         ! Clear sky total column abs solar flux
   real(r8),allocatable :: fsnsc_m(:)         ! Clear sky surface abs solar flux
   real(r8),allocatable :: fsdsc_m(:)         ! Clear sky surface downwelling solar flux
   real(r8),allocatable :: flut_m(:)          ! Upward flux at top of model
   real(r8),allocatable :: flutc_m(:)         ! Upward Clear Sky flux at top of model
   real(r8),allocatable :: flntc_m(:)         ! Clear sky lw flux at model top
   real(r8),allocatable :: flnsc_m(:)         ! Clear sky lw flux at srf (up-down)
   real(r8),allocatable :: fldsc_m(:)         ! Clear sky lw flux at srf (down)
   real(r8),allocatable :: flwds_m(:)          ! Down longwave flux at surface
   real(r8),allocatable :: fsns_m(:)          ! Surface solar absorbed flux
   real(r8),allocatable :: fsnr_m(:)
   real(r8),allocatable :: fsnt_m(:)          ! Net column abs solar flux at model top
   real(r8),allocatable :: flns_m(:)          ! Srf longwave cooling (up-down) flux
   real(r8),allocatable :: flnt_m(:)          ! Net outgoing lw flux at model top
   real(r8),allocatable :: flnr_m(:)
   real(r8),allocatable :: fsds_m(:)          ! Surface solar down flux
   real(r8),allocatable :: fln200_m(:)        ! net longwave flux interpolated to 200 mb
   real(r8),allocatable :: fln200c_m(:)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8),allocatable :: fsn200_m(:)        ! fns interpolated to 200 mb
   real(r8),allocatable :: fsn200c_m(:)       ! fcns interpolated to 200 mb
   real(r8),allocatable :: sols_m(:)          ! Solar downward visible direct  to surface
   real(r8),allocatable :: soll_m(:)          ! Solar downward near infrared direct  to surface
   real(r8),allocatable :: solsd_m(:)         ! Solar downward visible diffuse to surface
   real(r8),allocatable :: solld_m(:)         ! Solar downward near infrared diffuse to surface
   real(r8),allocatable :: qrs_m(:,:)
   real(r8),allocatable :: qrl_m(:,:)
   real(r8),allocatable :: qrsc_m(:,:)
   real(r8),allocatable :: qrlc_m(:,:)
   real(r8),allocatable :: su_m(:,:,:)  ! shortwave spectral flux up
   real(r8),allocatable :: sd_m(:,:,:)  ! shortwave spectral flux down
   real(r8),allocatable :: lu_m(:,:,:)  ! longwave  spectral flux up
   real(r8),allocatable :: ld_m(:,:,:)  ! longwave  spectral flux down
   real(r8),pointer :: su(:,:,:)  ! shortwave spectral flux up
   real(r8),pointer :: sd(:,:,:)  ! shortwave spectral flux down
   real(r8),pointer :: lu(:,:,:)  ! longwave  spectral flux up
   real(r8),pointer :: ld(:,:,:)  ! longwave  spectral flux down
   real(r8), allocatable :: dei_crm(:,:,:,:)  ! cloud scale ice effective diameter for optics
   real(r8), allocatable :: mu_crm(:,:,:,:)   ! cloud scale gamma parameter for optics
   real(r8), allocatable :: lambdac_crm(:,:,:,:)  ! cloud scale slope of droplet distribution for optics
   real(r8), allocatable :: des_crm(:,:,:,:)  ! cloud scale snow crystal diameter (micro-meter)
   real(r8), allocatable :: rel_crm(:,:,:,:)
   real(r8), allocatable :: rei_crm(:,:,:,:)
   real(r8), allocatable :: cld_tau_crm(:,:,:,:)
   real(r8), allocatable :: qrl_crm(:,:,:,:)
   real(r8), allocatable :: qrs_crm(:,:,:,:)
   real(r8), allocatable :: crm_fsnt(:,:,:)   ! net shortwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_fsntc(:,:,:)   ! net clear-sky shortwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_fsns(:,:,:)   ! net shortwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_fsnsc(:,:,:)   ! net clear-sky shortwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_flnt(:,:,:)   ! net longwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_flntc(:,:,:)   ! net clear-sky longwave fluxes at TOA at CRM grids
   real(r8), allocatable :: crm_flns(:,:,:)   ! net longwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_flnsc(:,:,:)   ! net clear-sky longwave fluxes at surface at CRM grids
   real(r8), allocatable :: crm_swcf(:,:,:)   ! shortwave cloud forcing at CRM grids


   real(r8), allocatable ::  crm_aodvisz(:,:,:,:)   ! layer aerosol optical depth at 550nm at CRM grids
   real(r8), allocatable ::  crm_aodvis(:,:,:)   ! AOD at 550nm at CRM grids
   real(r8), allocatable ::  crm_aod400(:,:,:)   ! AOD at 400nm at CRM grids
   real(r8), allocatable ::  crm_aod700(:,:,:)   ! AOD at 700nm at CRM grids
   real(r8), allocatable ::  aod400(:)   ! AOD at 400nm at CRM grids
   real(r8), allocatable ::  aod700(:)   ! AOD at 700nm at CRM grids

   real(r8), pointer :: t_rad (:,:,:) ! rad temperuture
   real(r8), pointer :: qv_rad(:,:,:) ! rad vapor
   real(r8), pointer :: qc_rad(:,:,:) ! rad cloud water
   real(r8), pointer :: qi_rad(:,:,:) ! rad cloud ice
   real(r8), pointer :: crm_qrad(:,:,:) ! rad heating

   real(r8), allocatable :: tot_cld_vistau_m(:,:)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: tot_icld_vistau_m(:,:)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: liq_icld_vistau_m(:,:)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: ice_icld_vistau_m(:,:)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
   real(r8), allocatable :: nct_tot_icld_vistau_m(:,:) ! the number of CRM columns that has in-cloud visible sw optical depth
   real(r8), allocatable :: nct_liq_icld_vistau_m(:,:) ! the number of CRM column that has liq in-cloud visible sw optical depth
   real(r8), allocatable :: nct_ice_icld_vistau_m(:,:) ! the number of CRM column that has ice in-cloud visible sw optical depth

    ! These do not need N_DIAG dimension
    real(r8),allocatable :: snow_tau(:,:,:) ! snow extinction optical depth

    real(r8),allocatable :: snow_lw_abs (:,:,:)   ! snow absorption optics depth (LW)

   ! Just used in m2005
   real(r8),allocatable :: snow_icld_vistau_m(:,:) ! snow in-cloud visible sw optical depth for output on history files
   real(r8),allocatable :: nct_snow_icld_vistau_m(:,:) ! the number of CRM column that has snow in-cloud visible sw optical depth


end type rad_avgdata_type_m2005

public :: tphysbc_spcam, spcam_register, spcam_init

integer :: dei_idx          = -1
integer :: mu_idx           = -1
integer :: lambdac_idx      = -1
integer :: des_idx          = -1
integer :: dgnumwet_crm_idx = -1
integer :: qaerwat_crm_idx  = -1
integer :: rel_idx          = -1
integer :: rei_idx          = -1
integer :: landm_idx        = -1
integer :: iciwp_idx        = -1
integer :: iclwp_idx        = -1
integer :: icswp_idx        = -1
integer :: cld_idx          = -1
integer :: dgnumwet_idx     = -1
integer :: qaerwat_idx      = -1
integer :: crm_t_rad_idx    = -1
integer :: crm_qc_rad_idx   = -1
integer :: crm_qi_rad_idx   = -1
integer :: crm_qv_rad_idx   = -1
integer :: crm_qrad_idx     = -1
integer :: crm_cld_rad_idx  = -1
integer :: crm_nc_rad_idx   = -1
integer :: crm_ni_rad_idx   = -1
integer :: crm_qs_rad_idx   = -1
integer :: crm_ns_rad_idx   = -1
integer :: cicewp_idx       = -1
integer :: cliqwp_idx       = -1
integer :: cldemis_idx      = -1
integer :: cldtau_idx       = -1
integer :: pmxrgn_idx       = -1
integer :: nmxrgn_idx       = -1
integer :: qrs_idx          = -1
integer :: qrl_idx          = -1
integer :: fsns_idx         = -1
integer :: fsnt_idx         = -1
integer :: flns_idx         = -1
integer :: flnt_idx         = -1
integer :: fsds_idx         = -1
integer :: cldfsnow_idx     = -1

! Minghuai - todo -- CAC note
!    These values will be "averaged" as appropriate and stored back in the pbuf
!    They should no longer be "saved"  -- Probably will  want to put in rad_avgdata structure
!    Email from Minghaui - 10/10/14 said to put on todo list as he did not have
!    time to address it now
!    real(r8),allocatable :: cicewp(:,:)
!    real(r8),allocatable :: cliqwp(:,:)
!    real(r8),allocatable :: rel(:,:)
!    real(r8),allocatable :: rei(:,:)
!    real(r8),allocatable :: dei(:,:)
!    real(r8),allocatable :: mu(:,:)
!    real(r8),allocatable :: lambdac(:,:)
!    real(r8),allocatable :: des(:,:)
!    real(r8),allocatable :: cld(:,:)        ! cloud fraction
!    real(r8),allocatable :: cldfsnow(:,:)   ! cloud fraction of just "snow clouds- whatever they are"
!    real(r8),allocatable :: csnowp(:,:)
!    real(r8),allocatable :: dgnumwet(:,:,:) ! number mode diameter
!    real(r8),allocatable :: qaerwat(:,:,:)  ! aerosol water


integer           :: nmodes
logical           :: is_spcam_m2005, is_spcam_sam1mom
logical           :: prog_modal_aero

contains
subroutine tphysbc_spcam (ztodt, state,   &
       tend,    pbuf,                     &
       cam_out, cam_in )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Evaluate and apply physical processes that are calculated BEFORE
    ! coupling to land, sea, and ice models.
    !
    ! Processes currently included are:
    !
    !  o Resetting Negative Tracers to Positive
    !  o Global Mean Total Energy Fixer
    !  o Dry Adjustment
    !  o Asymmetric Turbulence Scheme : Deep Convection & Shallow Convection
    !  o Stratiform Macro-Microphysics
    !  o Wet Scavenging of Aerosol
    !  o Radiation
    !
    ! Method:
    !
    ! Each parameterization should be implemented with this sequence of calls:
    !  1)  Call physics interface
    !  2)  Check energy
    !  3)  Call physics_update
    ! See Interface to Column Physics and Chemistry Packages
    !   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
    !
    !-----------------------------------------------------------------------

    use physics_buffer,  only : pbuf_old_tim_idx, dyn_time_lvls
    use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, &
         physics_state_check
    use dadadj_cam,      only: dadadj_tend
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use time_manager,    only: get_nstep
    use check_energy,    only: check_energy_chng, check_energy_fix
    use check_energy,    only: check_tracers_data, check_tracers_init
    use dycore,          only: dycore_is
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
    use tropopause,      only: tropopause_output
    use cam_abortutils,  only: endrun
#ifdef CRM
    use crm_physics,     only: crm_physics_tend
#endif
    use phys_control,    only: phys_getopts
    use sslt_rebin,      only: sslt_rebin_adv
    use qneg_module,     only: qneg3

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in)    :: ztodt                         ! 2 delta t (model time increment)

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


#ifdef CRM
    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_loc

    integer :: nstep                          ! current timestep number

    real(r8) :: net_flx(pcols)

    real(r8) cldn(pcols,pver)


    integer lchnk                              ! chunk identifier
    integer ncol                               ! number of atmospheric columns

    integer  i                                 ! index
    integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.

    ! physics buffer fields to compute tendencies for stratiform package
    integer itim_old, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore

    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble


    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: flx_heat(pcols)
    type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes

    logical           :: state_debug_checks  ! Debug physics_state.


    type(rad_avgdata_type_sam1mom)   :: rad_avgdata_sam1mom
    type(rad_avgdata_type_m2005)     :: rad_avgdata_m2005
    type(rad_out_t)                  :: rd

    integer :: teout_idx, qini_idx, cldliqini_idx, cldiceini_idx
    integer :: ii, jj
    !-----------------------------------------------------------------------
    call t_startf('bc_init')
    zero = 0._r8

    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()

    teout_idx = pbuf_get_index('TEOUT')
    qini_idx  = pbuf_get_index('QINI')
    cldliqini_idx  = pbuf_get_index('CLDLIQINI')
    cldiceini_idx  = pbuf_get_index('CLDICEINI')

    call phys_getopts(state_debug_checks_out=state_debug_checks)

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim_old/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))

    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )
    fracis (:ncol,:,1:pcnst) = 1._r8

    ! Set physics tendencies to 0
    tend %dTdt(:ncol,:pver)  = 0._r8
    tend %dudt(:ncol,:pver)  = 0._r8
    tend %dvdt(:ncol,:pver)  = 0._r8

    call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )

    ! Validate state coming from the dynamics.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")
    !
    ! Dump out "before physics" state
    !
    call diag_state_b4_phys_write (state)

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    call t_stopf('bc_init')

    !===================================================
    ! Global mean total energy fixer
    !===================================================
    call t_startf('energy_fixer')

    if (dycore_is('LR') .or. dycore_is('SE'))  then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
       call outfld('EFIX', flx_heat, pcols,lchnk)
    end if
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)


    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! T tendency due to dynamics
    if( nstep > dyn_time_lvls-1 ) then
       dtcore(:ncol,:pver) = (state%t(:ncol,:pver) - dtcore(:ncol,:pver))/(ztodt)
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')

    call sslt_rebin_adv(pbuf, state)

    !
    !===================================================
    ! Dry adjustment
    ! This code block is not a good example of interfacing a parameterization
    !===================================================
    call t_startf('dry_adjustment')

    call dadadj_tend (ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')

    ! -------------------------------------------------------------------------------
    ! Call cloud resolving model
    ! -------------------------------------------------------------------------------

    call crm_physics_tend(ztodt, state, tend, ptend, pbuf, cam_in)
    call physics_update(state, ptend, ztodt, tend)

    !===================================================
    ! Moist physical parameteriztions complete:
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, pbuf)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    if (is_spcam_sam1mom) then
       call spcam_radiation_setup_sam1mom(cam_in, cldn, state, pbuf, rad_avgdata_sam1mom, state_loc)
    else if (is_spcam_m2005) then
       call spcam_radiation_setup_m2005(state, pbuf, rad_avgdata_m2005, state_loc)
    end if

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')

    if (is_spcam_sam1mom) then
       do jj=1,crm_ny
          do ii=1,crm_nx
             call spcam_radiation_col_setup_sam1mom(ii, jj, state_loc, pbuf, rad_avgdata_sam1mom)
             call radiation_tend(state_loc, ptend, pbuf, &
                  cam_out, cam_in, &
                  net_flx, rd)
             call spcam_radiation_col_finalize_sam1mom(state, ii, jj, pbuf, rd, cam_out, rad_avgdata_sam1mom)
          end do
       end do
       call spcam_radiation_finalize_sam1mom(cam_in, state, pbuf, rad_avgdata_sam1mom, cam_out, cldn, net_flx, ptend)

    else if(is_spcam_m2005) then
       do jj=1,crm_ny
          do ii=1,crm_nx
             call spcam_radiation_col_setup_m2005(ii, jj, ixcldice, ixcldliq, state_loc, pbuf, rad_avgdata_m2005)
              call radiation_tend(state_loc, ptend, pbuf, &
                   cam_out, cam_in, &
                   net_flx, rd)
             call spcam_radiation_col_finalize_m2005(state, ii, jj, pbuf, rd, cam_out, rad_avgdata_m2005)
          end do
       end do
       call spcam_radiation_finalize_m2005(cam_in, state, pbuf, rad_avgdata_m2005, cam_out, net_flx, ptend)
    end if

    ! Set net flux used by spectral dycores
    do i=1,ncol
       tend%flx_net(i) = net_flx(i)
    end do

    ! don't add radiative tendency to GCM temperature in case of superparameterization
    ! as it was added above as part of crm tendency.
    ptend%s = 0._r8

    call physics_update(state, ptend, ztodt, tend)

    call check_energy_chng(state, tend, "spradheat", nstep, ztodt, zero, zero, zero, zero)

    call t_stopf('radiation')

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

#endif
end subroutine tphysbc_spcam

!===============================================================================

subroutine spcam_register()
  use physics_buffer, only: pbuf_add_field, dtype_r8, dyn_time_lvls ! is dyn_time_lvls needed ???
  use phys_control,   only: cam_physpkg_is
#ifdef CRM
  use crm_physics,      only: crm_physics_register
  use crmx_vars,        only: naer, vaer, hgaer
  use crmx_grid
#ifdef MODAL_AERO
  use modal_aero_data,  only: ntot_amode

  allocate(naer(nzm, ntot_amode))     ! Aerosol number concentration [/m3]
  allocate(vaer(nzm, ntot_amode))     ! aerosol volume concentration [m3/m3]
  allocate(hgaer(nzm, ntot_amode))    ! hygroscopicity of aerosol mode
#endif


  call crm_physics_register()

#endif

  is_spcam_m2005   = cam_physpkg_is('spcam_m2005')
  is_spcam_sam1mom = cam_physpkg_is('spcam_sam1mom')

  if (is_spcam_m2005) then
     call pbuf_add_field('ICSWP',    'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)
     call pbuf_add_field('CLDFSNOW', 'physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx)
  endif

end subroutine spcam_register

!===============================================================================

subroutine spcam_init(pbuf2d)
   use physics_buffer,   only: pbuf_get_index
   use phys_control,     only: phys_getopts
#ifdef CRM
   use crm_physics,      only: crm_physics_init
#endif
   use rad_constituents, only: rad_cnst_get_info

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CRM

   call phys_getopts(prog_modal_aero_out     = prog_modal_aero)

   call rad_cnst_get_info(0, nmodes=nmodes)

   dei_idx         = pbuf_get_index('DEI')
   mu_idx          = pbuf_get_index('MU')
   lambdac_idx     = pbuf_get_index('LAMBDAC')
   des_idx         = pbuf_get_index('DES')
   rel_idx         = pbuf_get_index('REL')
   rei_idx         = pbuf_get_index('REI')
   landm_idx       = pbuf_get_index('LANDM')
   cld_idx         = pbuf_get_index('CLD')
   qrs_idx         = pbuf_get_index('QRS')
   qrl_idx         = pbuf_get_index('QRL')
   fsns_idx        = pbuf_get_index('FSNS')
   fsds_idx        = pbuf_get_index('FSDS')
   fsnt_idx        = pbuf_get_index('FSNT')
   flnt_idx        = pbuf_get_index('FLNT')
   flns_idx        = pbuf_get_index('FLNS')

   crm_t_rad_idx   = pbuf_get_index('CRM_T_RAD')
   crm_qc_rad_idx  = pbuf_get_index('CRM_QC_RAD')
   crm_qi_rad_idx  = pbuf_get_index('CRM_QI_RAD')
   crm_qv_rad_idx  = pbuf_get_index('CRM_QV_RAD')
   crm_qrad_idx    = pbuf_get_index('CRM_QRAD')
   crm_cld_rad_idx = pbuf_get_index('CRM_CLD_RAD')


   if (is_spcam_sam1mom) then
      cldemis_idx     = pbuf_get_index('CLDEMIS')
      cldtau_idx      = pbuf_get_index('CLDTAU')
      cicewp_idx      = pbuf_get_index('CICEWP')
      cliqwp_idx      = pbuf_get_index('CLIQWP')
      pmxrgn_idx      = pbuf_get_index('PMXRGN')
      nmxrgn_idx      = pbuf_get_index('NMXRGN')
   else if (is_spcam_m2005) then
      iciwp_idx       = pbuf_get_index('ICIWP')
      iclwp_idx       = pbuf_get_index('ICLWP')
      crm_nc_rad_idx  = pbuf_get_index('CRM_NC_RAD')
      crm_ni_rad_idx  = pbuf_get_index('CRM_NI_RAD')
      crm_qs_rad_idx  = pbuf_get_index('CRM_QS_RAD')
      crm_ns_rad_idx  = pbuf_get_index('CRM_NS_RAD')
   end if

   if (prog_modal_aero) then
      dgnumwet_idx     = pbuf_get_index('DGNUMWET')
      qaerwat_idx      = pbuf_get_index('QAERWAT')
      dgnumwet_crm_idx = pbuf_get_index('CRM_DGNUMWET')
      qaerwat_crm_idx  = pbuf_get_index('CRM_QAERWAT')
   end if

   ! Initialize the crm_physics layer
   call crm_physics_init(pbuf2d)

#endif
end subroutine spcam_init

!===============================================================================

subroutine spcam_radiation_setup_m2005(state, pbuf, rad_avgdata,  state_loc)

   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
   use physics_buffer,   only: pbuf_old_tim_idx

   type(physics_state),             intent(in)             :: state
   type(physics_buffer_desc),       intent(inout), pointer :: pbuf(:)

   type(rad_avgdata_type_m2005),    intent(out)            :: rad_avgdata
   type(physics_state),             intent(out)            :: state_loc

#ifdef m2005
   real(r8), pointer, dimension(:, :)  ::  cicewp
   real(r8), pointer, dimension(:, :)  ::  cliqwp
   real(r8), pointer, dimension(:, :)  ::  csnowp
   real(r8), pointer, dimension(:,:)   ::  rel      ! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:)   ::  rei      ! ice effective drop size (microns)
   real(r8), pointer, dimension(:,:)   ::  cld      ! cloud fraction
   real(r8), pointer, dimension(:,:)   ::  cldfsnow ! cloud fraction of just "snow clouds- whatever they are"
   real(r8), pointer, dimension(:, :)  ::  dei      ! ice effective diameter for optics (radiation)
   real(r8), pointer, dimension(:, :)  ::  mu       ! gamma parameter for optics (radiation)
   real(r8), pointer, dimension(:, :)  ::  lambdac  ! slope of droplet distribution for optics (radiation)
   real(r8), pointer, dimension(:, :)  ::  des      ! snow crystatl diameter for optics (mirometer, radiation)

   integer :: ncol                               ! number of atmospheric columns
   integer :: itim_old

   ncol  = state%ncol

   call physics_state_copy(state, state_loc)

   allocate(rad_avgdata%solin_m      (pcols))
   allocate(rad_avgdata%fsntoa_m     (pcols))
   allocate(rad_avgdata%fsutoa_m     (pcols))
   allocate(rad_avgdata%fsntoac_m    (pcols))
   allocate(rad_avgdata%fsnirt_m     (pcols))
   allocate(rad_avgdata%fsnrtc_m     (pcols))
   allocate(rad_avgdata%fsnirtsq_m   (pcols))
   allocate(rad_avgdata%fsntc_m      (pcols))
   allocate(rad_avgdata%fsnsc_m      (pcols))
   allocate(rad_avgdata%fsdsc_m      (pcols))
   allocate(rad_avgdata%flut_m       (pcols))
   allocate(rad_avgdata%flutc_m      (pcols))
   allocate(rad_avgdata%flntc_m      (pcols))
   allocate(rad_avgdata%flnsc_m      (pcols))
   allocate(rad_avgdata%fldsc_m      (pcols))
   allocate(rad_avgdata%flwds_m      (pcols))
   allocate(rad_avgdata%fsns_m       (pcols))
   allocate(rad_avgdata%fsnr_m       (pcols))
   allocate(rad_avgdata%fsnt_m       (pcols))
   allocate(rad_avgdata%flns_m       (pcols))
   allocate(rad_avgdata%flnt_m       (pcols))
   allocate(rad_avgdata%flnr_m       (pcols))
   allocate(rad_avgdata%fsds_m       (pcols))
   allocate(rad_avgdata%fln200_m     (pcols))
   allocate(rad_avgdata%fln200c_m    (pcols))
   allocate(rad_avgdata%fsn200_m     (pcols))
   allocate(rad_avgdata%fsn200c_m    (pcols))
   allocate(rad_avgdata%sols_m       (pcols))
   allocate(rad_avgdata%soll_m       (pcols))
   allocate(rad_avgdata%solsd_m      (pcols))
   allocate(rad_avgdata%solld_m      (pcols))
   allocate(rad_avgdata%qrs_m        (pcols,pver))
   allocate(rad_avgdata%qrl_m        (pcols,pver))
   allocate(rad_avgdata%qrsc_m       (pcols,pver))
   allocate(rad_avgdata%qrlc_m       (pcols,pver))
   allocate(rad_avgdata%rel_crm      (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%rei_crm      (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%cld_tau_crm  (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%qrl_crm      (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%qrs_crm      (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%crm_fsnt     (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_fsntc    (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_fsns     (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_fsnsc    (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_flnt     (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_flntc    (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_flns     (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_flnsc    (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_swcf     (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_aodvisz  (pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%crm_aodvis   (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_aod400   (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%crm_aod700   (pcols, crm_nx, crm_ny))
   allocate(rad_avgdata%aod400       (pcols))
   allocate(rad_avgdata%aod700       (pcols))

   allocate(rad_avgdata%tot_cld_vistau_m      (pcols,pver))
   allocate(rad_avgdata%tot_icld_vistau_m     (pcols,pver))
   allocate(rad_avgdata%liq_icld_vistau_m     (pcols,pver))
   allocate(rad_avgdata%ice_icld_vistau_m     (pcols,pver))
   allocate(rad_avgdata%nct_tot_icld_vistau_m (pcols,pver))
   allocate(rad_avgdata%nct_liq_icld_vistau_m (pcols,pver))
   allocate(rad_avgdata%nct_ice_icld_vistau_m (pcols,pver))
   allocate(rad_avgdata%snow_icld_vistau_m    (pcols,pver))
   allocate(rad_avgdata%nct_snow_icld_vistau_m(pcols,pver))

   allocate(rad_avgdata%dei_crm(pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%mu_crm(pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%lambdac_crm(pcols, crm_nx, crm_ny, crm_nz))
   allocate(rad_avgdata%des_crm(pcols, crm_nx, crm_ny, crm_nz))

   call pbuf_get_field(pbuf, iciwp_idx, cicewp)
   call pbuf_get_field(pbuf, iclwp_idx, cliqwp)
   call pbuf_get_field(pbuf, icswp_idx, csnowp)
   call pbuf_get_field(pbuf, rel_idx,    rel)
   call pbuf_get_field(pbuf, rei_idx,    rei)
   call pbuf_get_field(pbuf, dei_idx, dei)
   call pbuf_get_field(pbuf, mu_idx, mu)
   call pbuf_get_field(pbuf, lambdac_idx, lambdac)
   call pbuf_get_field(pbuf, des_idx, des)

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   endif

   ! Initialize the summation values

   rad_avgdata%solin_m     = 0._r8
   rad_avgdata%fsntoa_m    = 0._r8
   rad_avgdata%fsutoa_m    = 0._r8
   rad_avgdata%fsntoac_m   = 0._r8
   rad_avgdata%fsnirt_m    = 0._r8
   rad_avgdata%fsnrtc_m    = 0._r8
   rad_avgdata%fsnirtsq_m  = 0._r8
   rad_avgdata%fsntc_m     = 0._r8
   rad_avgdata%fsnsc_m     = 0._r8
   rad_avgdata%fsdsc_m     = 0._r8
   rad_avgdata%flut_m      = 0._r8
   rad_avgdata%flutc_m     = 0._r8
   rad_avgdata%flntc_m     = 0._r8
   rad_avgdata%flnsc_m     = 0._r8
   rad_avgdata%fldsc_m     = 0._r8
   rad_avgdata%flwds_m     = 0._r8
   rad_avgdata%fsns_m      = 0._r8
   rad_avgdata%fsnt_m      = 0._r8
   rad_avgdata%flns_m      = 0._r8
   rad_avgdata%flnt_m      = 0._r8
   rad_avgdata%flnr_m      = 0._r8
   rad_avgdata%fsds_m      = 0._r8
   rad_avgdata%fsnr_m      = 0._r8
   rad_avgdata%fln200_m    = 0._r8
   rad_avgdata%fln200c_m   = 0._r8
   rad_avgdata%fsn200_m    = 0._r8
   rad_avgdata%fsn200c_m   = 0._r8
   rad_avgdata%sols_m      = 0._r8
   rad_avgdata%soll_m      = 0._r8
   rad_avgdata%solsd_m     = 0._r8
   rad_avgdata%solld_m     = 0._r8
   rad_avgdata%qrs_m       = 0._r8
   rad_avgdata%qrl_m       = 0._r8
   rad_avgdata%qrsc_m      = 0._r8
   rad_avgdata%qrlc_m      = 0._r8
   rad_avgdata%qrs_crm     = 0._r8
   rad_avgdata%qrl_crm     = 0._r8
   rad_avgdata%cld_tau_crm = 0.0_r8
   rad_avgdata%crm_aodvisz =  0._r8
   rad_avgdata%crm_aodvis  =  0._r8

   rad_avgdata%crm_aod400  = 0._r8  ; rad_avgdata%crm_aod700 = 0._r8
   rad_avgdata%aod400      = 0._r8  ; rad_avgdata%aod700     = 0._r8
   rad_avgdata%crm_fsnt    = 0._r8  ; rad_avgdata%crm_fsntc  = 0._r8
   rad_avgdata%crm_fsns    = 0._r8  ; rad_avgdata%crm_fsnsc  = 0._r8
   rad_avgdata%crm_flnt    = 0._r8  ; rad_avgdata%crm_flntc  = 0._r8
   rad_avgdata%crm_flns    = 0._r8  ; rad_avgdata%crm_flnsc  = 0._r8
   rad_avgdata%crm_swcf    = 0._r8


   rad_avgdata%tot_cld_vistau_m   = 0._r8
   rad_avgdata%tot_icld_vistau_m  = 0._r8  ; rad_avgdata%nct_tot_icld_vistau_m  = 0._r8
   rad_avgdata%liq_icld_vistau_m  = 0._r8  ; rad_avgdata%nct_liq_icld_vistau_m  = 0._r8
   rad_avgdata%ice_icld_vistau_m  = 0._r8  ; rad_avgdata%nct_ice_icld_vistau_m  = 0._r8
   rad_avgdata%snow_icld_vistau_m = 0._r8  ; rad_avgdata%nct_snow_icld_vistau_m = 0._r8

   ! Initialize the pbuf values
   lambdac               = 0.0_r8
   des                   = 0.0_r8
   cicewp(1:ncol,1:pver) = 0.0_r8
   cliqwp(1:ncol,1:pver) = 0.0_r8
   csnowp(1:ncol,1:pver) = 0.0_r8
   cld                   = 0.0_r8
   cldfsnow              = 0.0_r8
   rel                   = 0.0_r8
   rei                   = 0.0_r8
   dei                   = 0.0_r8
   mu                    = 0.0_r8

#endif
end subroutine spcam_radiation_setup_m2005

!===============================================================================

subroutine spcam_radiation_col_setup_m2005(ii, jj, ixcldice, ixcldliq, state_loc, pbuf, rad_avgdata)

   use physics_buffer, only: pbuf_old_tim_idx
   use physconst,      only: gravit
#ifdef CRM
    use crm_physics,    only: m2005_effradius
#endif


   integer,                         intent(in)             :: ii,jj
   integer,                         intent(in)             :: ixcldice, ixcldliq ! constituent indices for cloud liq and ice water.

   type(physics_state),             intent(inout)          :: state_loc
   type(physics_buffer_desc),       intent(inout), pointer :: pbuf(:)
   type(rad_avgdata_type_m2005),    intent(inout)          :: rad_avgdata

#ifdef m2005
   real(r8),pointer :: nc_rad(:,:,:,:) ! rad cloud water droplet number (#/kg)
   real(r8),pointer :: ni_rad(:,:,:,:) ! rad cloud ice crystal nubmer (#/kg)
   real(r8),pointer :: qs_rad(:,:,:,:) ! rad cloud snow crystal mass (kg/kg)
   real(r8),pointer :: ns_rad(:,:,:,:) ! rad cloud snow crystal nubmer (#/kg)


   real(r8),pointer :: t_rad (:,:,:,:) ! rad temperuture
   real(r8),pointer :: qv_rad(:,:,:,:) ! rad vapor
   real(r8),pointer :: qc_rad(:,:,:,:) ! rad cloud water
   real(r8),pointer :: qi_rad(:,:,:,:) ! rad cloud ice
   real(r8),pointer :: crm_qrad(:,:,:,:) ! rad heating
   real(r8),pointer :: cld_rad(:,:,:,:) ! rad cloud fraction


   real(r8), pointer, dimension(:,:)   :: cicewp
   real(r8), pointer, dimension(:,:)   :: cliqwp
   real(r8), pointer, dimension(:,:)   :: csnowp
   real(r8), pointer, dimension(:,:)   :: rel      ! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:)   :: rei      ! ice effective drop size (microns)
   real(r8), pointer, dimension(:,:)   :: cld      ! cloud fraction
   real(r8), pointer, dimension(:,:)   :: cldfsnow ! cloud fraction of just "snow clouds- whatever they are"
   real(r8), pointer, dimension(:,:)   :: dei      ! ice effective diameter for optics (radiation)
   real(r8), pointer, dimension(:,:)   :: mu       ! gamma parameter for optics (radiation)
   real(r8), pointer, dimension(:,:)   :: lambdac  ! slope of droplet distribution for optics (radiation)
   real(r8), pointer, dimension(:,:)   :: des      ! snow crystatl diameter for optics (mirometer, radiation)
   real(r8), pointer, dimension(:,:,:) :: dgnumwet ! number mode diameter
   real(r8), pointer, dimension(:,:,:) :: qaerwat  ! aerosol water

   real(r8),pointer, dimension(:,:,:,:,:) :: qaerwat_crm  ! aerosol water
   real(r8),pointer, dimension(:,:,:,:,:) :: dgnumwet_crm ! wet mode dimaeter

   real(r8) :: qtot
   real(r8) :: effl     ! droplet effective radius [micrometer]
   real(r8) :: effi     ! ice crystal effective radius [micrometer]
   real(r8) :: effl_fn  ! effl for fixed number concentration of nlic = 1.e8

   real(r8) :: deffi    ! ice effective diameter for optics (radiation)
   real(r8) :: lamc     ! slope of droplet distribution for optics (radiation)
   real(r8) :: pgam     ! gamma parameter for optics (radiation)
   real(r8) :: dest     ! snow crystal effective diameters for optics (radiation) (micro-meter)


   integer :: itim_old
   integer :: m, k, i
   integer :: ncol                               ! number of atmospheric columns

   ncol  = state_loc%ncol

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, dei_idx, dei)
   call pbuf_get_field(pbuf, mu_idx, mu)
   call pbuf_get_field(pbuf, lambdac_idx, lambdac)
   call pbuf_get_field(pbuf, des_idx, des)
   if (prog_modal_aero) then
      call pbuf_get_field(pbuf, dgnumwet_crm_idx, dgnumwet_crm)
      call pbuf_get_field(pbuf, qaerwat_crm_idx,  qaerwat_crm)
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat)
   endif

   call pbuf_get_field(pbuf, rel_idx,    rel)
   call pbuf_get_field(pbuf, rei_idx,    rei)

   call pbuf_get_field(pbuf, crm_t_rad_idx,   t_rad)
   call pbuf_get_field(pbuf, crm_qc_rad_idx,  qc_rad)
   call pbuf_get_field(pbuf, crm_qi_rad_idx,  qi_rad)
   call pbuf_get_field(pbuf, crm_qv_rad_idx,  qv_rad)
   call pbuf_get_field(pbuf, crm_qrad_idx,    crm_qrad)
   call pbuf_get_field(pbuf, crm_cld_rad_idx, cld_rad)

   crm_qrad=0._r8


   call pbuf_get_field(pbuf, iciwp_idx, cicewp)
   call pbuf_get_field(pbuf, iclwp_idx, cliqwp)
   call pbuf_get_field(pbuf, icswp_idx, csnowp)

   call pbuf_get_field(pbuf, crm_nc_rad_idx, nc_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
   call pbuf_get_field(pbuf, crm_ni_rad_idx, ni_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
   call pbuf_get_field(pbuf, crm_qs_rad_idx, qs_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))
   call pbuf_get_field(pbuf, crm_ns_rad_idx, ns_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx, crm_ny, crm_nz/))

   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   endif

   do m=1,crm_nz
      k = pver-m+1
      do i=1,ncol

          qtot = qc_rad(i,ii,jj,m) + qi_rad(i,ii,jj,m)
          if(qtot.gt.1.e-9_r8) then
             cld(i,k)           = min(0.99_r8, cld_rad(i,ii,jj,m))

             ! In-cloud ice water path.
             cicewp(i,k)        = qi_rad(i,ii,jj,m)*state_loc%pdel(i,k)/gravit / max(0.01_r8,cld(i,k))
             ! In-cloud liquid water path.
             cliqwp(i,k)        = qc_rad(i,ii,jj,m)*state_loc%pdel(i,k)/gravit / max(0.01_r8,cld(i,k))
          else
             cld(i,k)           = 0._r8
             cicewp(i,k)        = 0._r8           ! In-cloud ice water path.
             cliqwp(i,k)        = 0._r8           ! In-cloud liquid water path.
          end if

          !
          ! snow water-related variables:
          ! snow water is an important component in m2005 microphysics, and is therefore taken
          ! account in the radiative calculation (snow water path is several times larger than ice water path in m2005 globally).
          !
          if( qs_rad(i, ii, jj, m).gt.1.0e-7_r8) then
             cldfsnow(i,k) = 0.99_r8
             csnowp(i,k)   = qs_rad(i,ii,jj,m)*state_loc%pdel(i,k)/gravit / max(0.001_r8,cldfsnow(i,k))
          else
             cldfsnow(i,k) = 0.0_r8
             csnowp(i,k)   = 0.0_r8
          end if


          ! update ice water, liquid water, water vapor, and temperature in state_loc
          state_loc%q(i,k,ixcldice) =  qi_rad(i,ii,jj,m)
          state_loc%q(i,k,ixcldliq) =  qc_rad(i,ii,jj,m)
          state_loc%q(i,k,1)        =  max(1.e-9_r8,qv_rad(i,ii,jj,m))
          state_loc%t(i,k)          =  t_rad(i, ii, jj, m)

          !    Using CRM scale aerosol water to calculate aerosol optical depth.
          !    Here we assume no aerosol water uptake at cloudy sky at CRM grids.
          !    This is not really phyisically correct. But if we assume 100% of relative humidity for
          !    aerosol water uptake, this will bias 'AODVIS' to be large, since 'AODVIS' is used
          !    to compare with observated clear sky AOD. In the future, AODVIS is needed to be calcualted
          !    from clear sky CRM AOD only. But before this is done, we will assume no water uptake at CCRM
          !    cloudy grids (The radiative effects of this assumption will be small, since in cloudy sky,
          !    aerosol effects is small anyway.
          !
          if (prog_modal_aero) then
             qaerwat(i, k, 1:nmodes)  = qaerwat_crm(i, ii, jj, m, 1:nmodes)
             dgnumwet(i, k, 1:nmodes) = dgnumwet_crm(i, ii, jj, m, 1:nmodes)
          endif
      end do ! i
   end do ! m


   ! update effective radius
   do m=1,crm_nz
      k = pver-m+1
      do i=1,ncol

         call m2005_effradius(qc_rad(i,ii,jj,m), nc_rad(i,ii,jj,m), qi_rad(i,ii,jj,m), &
                              ni_rad(i,ii,jj,m), qs_rad(i,ii,jj,m), ns_rad(i,ii,jj,m),  &
                              1.0_r8, state_loc%pmid(i,k), state_loc%t(i,k), effl, effi, effl_fn, deffi, lamc, pgam, dest)

         rel(i,k)     = effl
         rei(i,k)     = effi
         dei(i,k)     = deffi
         mu(i,k)      = pgam
         lambdac(i,k) = lamc
         des(i,k)     = dest

         rad_avgdata%dei_crm(i,ii,jj,m)     = dei(i,k)
         rad_avgdata%mu_crm(i,ii,jj,m)      = mu(i,k)
         rad_avgdata%lambdac_crm(i,ii,jj,m) = lambdac(i,k)
         rad_avgdata%des_crm(i,ii,jj,m)     = des(i,k)
         rad_avgdata%rel_crm(i,ii,jj,m)     = rel(i,k)
         rad_avgdata%rei_crm(i,ii,jj,m)     = rei(i,k)
      end do
   end do

#endif
end subroutine spcam_radiation_col_setup_m2005

!===============================================================================

subroutine spcam_radiation_finalize_m2005(cam_in, state, pbuf, rad_avgdata, cam_out, net_flx, ptend)

   use physconst,       only: cpair
   use rad_constituents,only: rad_cnst_out

   use physconst,       only: cappa
   use radiation_data,  only: rad_data_write
   use radheat,         only: radheat_tend
   use time_manager,    only: get_curr_calday
   use physics_buffer,  only: pbuf_old_tim_idx
   use radheat,         only: radheat_tend
   use orbit,           only: zenith

   type(cam_in_t),                  intent(in) :: cam_in
   type(physics_state),             intent(in) :: state


   type(physics_buffer_desc),       intent(inout), pointer      :: pbuf(:)
   type(rad_avgdata_type_m2005),    intent(inout) :: rad_avgdata
   type(cam_out_t),                 intent(inout) :: cam_out

   real(r8),                        intent(inout) :: net_flx(pcols)
   type(physics_ptend),             intent(out)   :: ptend            ! indivdual parameterization tendencies



#ifdef m2005

   real(r8), parameter :: factor_xy = 1._r8/dble(crm_nx*crm_ny)

   real(r8), pointer, dimension(:,:)   ::  cicewp
   real(r8), pointer, dimension(:,:)   ::  cliqwp
   real(r8), pointer, dimension(:,:)   ::  csnowp
   real(r8), pointer, dimension(:,:)   ::  rel      ! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:)   ::  rei      ! ice effective drop size (microns)
   real(r8), pointer, dimension(:,:)   ::  landm
   real(r8), pointer, dimension(:,:)   ::  cld      ! cloud fraction
   real(r8), pointer, dimension(:,:)   ::  cldfsnow ! cloud fraction of just "snow clouds- whatever they are"
   real(r8), pointer, dimension(:,:)   ::  dei      ! ice effective diameter for optics (radiation)
   real(r8), pointer, dimension(:,:)   ::  mu       ! gamma parameter for optics (radiation)
   real(r8), pointer, dimension(:,:)   ::  lambdac  ! slope of droplet distribution for optics (radiation)
   real(r8), pointer, dimension(:,:)   ::  des      ! snow crystatl diameter for optics (mirometer, radiation)
   real(r8), pointer, dimension(:,:,:) :: dgnumwet ! number mode diameter
   real(r8), pointer, dimension(:,:,:) :: qaerwat ! aerosol water
   real(r8), pointer, dimension(:,:,:,:) :: crm_qrad ! rad heating
   real(r8), pointer, dimension(:,:)   :: qrs
   real(r8), pointer, dimension(:,:)   :: qrl
   real(r8), pointer, dimension(:)     :: fsns      ! Surface solar absorbed flux
   real(r8), pointer, dimension(:)     :: fsnt      ! Net column abs solar flux at model top
   real(r8), pointer, dimension(:)     :: flns      ! Srf longwave cooling (up-down) flux
   real(r8), pointer, dimension(:)     :: flnt      ! Net outgoing lw flux at model top
   real(r8), pointer, dimension(:)     :: fsds      ! Surface solar down flux

   integer :: lchnk                              ! chunk identifier
   integer :: ncol                               ! number of atmospheric columns

   integer :: Nday                      ! Number of daylight columns
   integer :: Nnite                     ! Number of night columns
   integer :: itim_old
   integer :: i, k, m

   integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)
   real(r8) :: coszrs(pcols)                 ! Cosine solar zenith angle


   lchnk = state%lchnk
   ncol  = state%ncol

   calday = get_curr_calday()

   !
   ! Cosine solar zenith angle for current time step
   !
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call zenith (calday, clat, clon, coszrs, ncol)

   ! Gather night/day column indices.
   Nday = 0
   Nnite = 0
   do i = 1, ncol
      if ( coszrs(i) > 0.0_r8 ) then
         Nday = Nday + 1
      else
         Nnite = Nnite + 1
         IdxNite(Nnite) = i
      end if
   end do



   ! Shortwave

   ftem(:ncol,:pver) = rad_avgdata%qrs_m(:ncol,:pver)/cpair
   call outfld('QRS'//'    ',ftem  ,pcols,lchnk)
   ftem(:ncol,:pver) = rad_avgdata%qrsc_m(:ncol,:pver)/cpair
   call outfld('QRSC'//'    ',ftem  ,pcols,lchnk)
   call outfld('SOLIN'//'    ',rad_avgdata%solin_m(:) ,pcols,lchnk)
   call outfld('FSDS'//'    ',rad_avgdata%fsds_m(:)  ,pcols,lchnk)
   call outfld('FSNIRTOA'//'    ',rad_avgdata%fsnirt_m(:),pcols,lchnk)
   call outfld('FSNRTOAC'//'    ',rad_avgdata%fsnrtc_m(:),pcols,lchnk)
   call outfld('FSNRTOAS'//'    ',rad_avgdata%fsnirtsq_m(:),pcols,lchnk)
   call outfld('FSNT'//'    ',rad_avgdata%fsnt_m(:)  ,pcols,lchnk)
   call outfld('FSNS'//'    ',rad_avgdata%fsns_m(:)  ,pcols,lchnk)
   call outfld('FSNTC'//'    ',rad_avgdata%fsntc_m(:) ,pcols,lchnk)
   call outfld('FSNSC'//'    ',rad_avgdata%fsnsc_m(:) ,pcols,lchnk)
   call outfld('FSDSC'//'    ',rad_avgdata%fsdsc_m(:) ,pcols,lchnk)
   call outfld('FSNTOA'//'    ',rad_avgdata%fsntoa_m(:),pcols,lchnk)
   call outfld('FSUTOA'//'    ',rad_avgdata%fsutoa_m(:),pcols,lchnk)
   call outfld('FSNTOAC'//'    ',rad_avgdata%fsntoac_m(:),pcols,lchnk)
   call outfld('SOLS'//'    ',rad_avgdata%sols_m(:)  ,pcols,lchnk)
   call outfld('SOLL'//'    ',rad_avgdata%soll_m(:)  ,pcols,lchnk)
   call outfld('SOLSD'//'    ',rad_avgdata%solsd_m(:) ,pcols,lchnk)
   call outfld('SOLLD'//'    ',rad_avgdata%solld_m(:) ,pcols,lchnk)
   call outfld('FSN200'//'    ',rad_avgdata%fsn200_m(:),pcols,lchnk)
   call outfld('FSN200C'//'    ',rad_avgdata%fsn200c_m(:),pcols,lchnk)
   call outfld('SWCF'//'    ',rad_avgdata%fsntoa_m(:)-rad_avgdata%fsntoac_m(:)  ,pcols,lchnk)
   call outfld('FSNR'//'    ',rad_avgdata%fsnr_m(:) ,pcols,lchnk)

   do i = 1, nnite
      rad_avgdata%crm_aodvis(idxnite(i), :, :)     = fillvalue
      rad_avgdata%crm_aod400(idxnite(i), :, :)     = fillvalue
      rad_avgdata%crm_aod700(idxnite(i), :, :)     = fillvalue
      rad_avgdata%aod400(idxnite(i))               = fillvalue
      rad_avgdata%aod700(idxnite(i))               = fillvalue
      rad_avgdata%crm_aodvisz(idxnite(i), :, :, :) = fillvalue
      rad_avgdata%tot_cld_vistau_m(IdxNite(i),:)   = fillvalue
      rad_avgdata%tot_icld_vistau_m(IdxNite(i),:)  = fillvalue
      rad_avgdata%liq_icld_vistau_m(IdxNite(i),:)  = fillvalue
      rad_avgdata%ice_icld_vistau_m(IdxNite(i),:)  = fillvalue
      if (cldfsnow_idx > 0) then
         rad_avgdata%snow_icld_vistau_m(IdxNite(i),:) = fillvalue
      endif
   end do

   call outfld('CRM_FSNT',        rad_avgdata%crm_fsnt,          pcols, lchnk)
   call outfld('CRM_FSNTC',       rad_avgdata%crm_fsntc,         pcols, lchnk)
   call outfld('CRM_FSNS',        rad_avgdata%crm_fsns,          pcols, lchnk)
   call outfld('CRM_FSNSC',       rad_avgdata%crm_fsnsc,         pcols, lchnk)
   call outfld('CRM_AODVIS',      rad_avgdata%crm_aodvis,        pcols, lchnk)
   call outfld('CRM_AOD400',      rad_avgdata%crm_aod400,        pcols, lchnk)
   call outfld('CRM_AOD700',      rad_avgdata%crm_aod700,        pcols, lchnk)
   call outfld('AOD400',          rad_avgdata%aod400,            pcols, lchnk)
   call outfld('AOD700',          rad_avgdata%aod700,            pcols, lchnk)
   call outfld('CRM_AODVISZ',     rad_avgdata%crm_aodvisz,       pcols, lchnk)
   call outfld('TOT_CLD_VISTAU',  rad_avgdata%tot_cld_vistau_m,  pcols, lchnk)
   call outfld('TOT_ICLD_VISTAU', rad_avgdata%tot_icld_vistau_m, pcols, lchnk)
   call outfld('LIQ_ICLD_VISTAU', rad_avgdata%liq_icld_vistau_m, pcols, lchnk)
   call outfld('ICE_ICLD_VISTAU', rad_avgdata%ice_icld_vistau_m, pcols, lchnk)
   if (cldfsnow_idx > 0) then
      call outfld('SNOW_ICLD_VISTAU', rad_avgdata%snow_icld_vistau_m, pcols, lchnk)
   endif

   ! Longwave
   call outfld('QRL'//'    ',rad_avgdata%qrl_m (:ncol,:)/cpair,ncol,lchnk)
   call outfld('QRLC'//'    ',rad_avgdata%qrlc_m(:ncol,:)/cpair,ncol,lchnk)
   call outfld('FLNT'//'    ',rad_avgdata%flnt_m(:)  ,pcols,lchnk)
   call outfld('FLUT'//'    ',rad_avgdata%flut_m(:)  ,pcols,lchnk)
   call outfld('FLUTC'//'    ',rad_avgdata%flutc_m(:) ,pcols,lchnk)
   call outfld('FLNTC'//'    ',rad_avgdata%flntc_m(:) ,pcols,lchnk)
   call outfld('FLNS'//'    ',rad_avgdata%flns_m(:)  ,pcols,lchnk)

   call outfld('FLDSC'//'    ',rad_avgdata%fldsc_m(:) ,pcols,lchnk)
   call outfld('FLNSC'//'    ',rad_avgdata%flnsc_m(:) ,pcols,lchnk)
   call outfld('LWCF'//'    ',rad_avgdata%flutc_m(:)-rad_avgdata%flut_m(:)  ,pcols,lchnk)
   call outfld('FLN200'//'    ',rad_avgdata%fln200_m(:),pcols,lchnk)
   call outfld('FLN200C'//'    ',rad_avgdata%fln200c_m(:),pcols,lchnk)
   call outfld('FLDS'//'    ',rad_avgdata%flwds_m(:) ,pcols,lchnk)
   call outfld('FLNR'//'    ',rad_avgdata%flnr_m(:),pcols,lchnk)

   call outfld('CRM_FLNT',  rad_avgdata%crm_flnt,  pcols, lchnk)
   call outfld('CRM_FLNTC', rad_avgdata%crm_flntc, pcols, lchnk)
   call outfld('CRM_FLNS',  rad_avgdata%crm_flns,  pcols, lchnk)
   call outfld('CRM_FLNSC', rad_avgdata%crm_flnsc, pcols, lchnk)

   call outfld('CRM_REL',    rad_avgdata%rel_crm,     pcols, lchnk)
   call outfld('CRM_REI',    rad_avgdata%rei_crm,     pcols, lchnk)
   call outfld('CRM_MU',     rad_avgdata%mu_crm,      pcols, lchnk)
   call outfld('CRM_DEI',    rad_avgdata%dei_crm,     pcols, lchnk)
   call outfld('CRM_DES',    rad_avgdata%des_crm,     pcols, lchnk)
   call outfld('CRM_LAMBDA', rad_avgdata%lambdac_crm, pcols, lchnk)
   call outfld('CRM_TAU',    rad_avgdata%cld_tau_crm, pcols, lchnk)
   call outfld('CRM_QRL',    rad_avgdata%qrl_crm,     pcols, lchnk)
   call outfld('CRM_QRS',    rad_avgdata%qrs_crm,     pcols, lchnk)



   do i=1, ncol
      do k=1, pver
          rad_avgdata%tot_cld_vistau_m(i,k) = rad_avgdata%tot_icld_vistau_m(i,k) *  factor_xy
          if(rad_avgdata%nct_tot_icld_vistau_m(i,k).ge. 0.1_r8) then
             rad_avgdata%tot_icld_vistau_m(i,k)  = rad_avgdata%tot_icld_vistau_m(i,k)/rad_avgdata%nct_tot_icld_vistau_m(i,k)
          else
             rad_avgdata%tot_icld_vistau_m(i,k)  = 0.0_r8
          end if

          if(rad_avgdata%nct_liq_icld_vistau_m(i,k).ge. 0.1_r8) then
            rad_avgdata%liq_icld_vistau_m(i,k)  = rad_avgdata%liq_icld_vistau_m(i,k)/rad_avgdata%nct_liq_icld_vistau_m(i,k)
          else
             rad_avgdata%liq_icld_vistau_m(i,k)  = 0.0_r8
          end if

          if(rad_avgdata%nct_ice_icld_vistau_m(i,k).ge. 0.1_r8) then
            rad_avgdata%ice_icld_vistau_m(i,k)  = rad_avgdata%ice_icld_vistau_m(i,k)/rad_avgdata%nct_ice_icld_vistau_m(i,k)
          else
             rad_avgdata%ice_icld_vistau_m(i,k)  = 0.0_r8
          end if

          if(rad_avgdata%nct_snow_icld_vistau_m(i,k).ge. 0.1_r8) then
            rad_avgdata%snow_icld_vistau_m(i,k) = rad_avgdata%snow_icld_vistau_m(i,k)/rad_avgdata%nct_snow_icld_vistau_m(i,k)
          else
             rad_avgdata%snow_icld_vistau_m(i,k) = 0.0_r8
          end if

      end do
   end do

   ! Output aerosol mmr
   call rad_cnst_out(0, state, pbuf)


   ! restore to the non-spcam values

   call pbuf_get_field(pbuf, iciwp_idx,    cicewp)
   call pbuf_get_field(pbuf, iclwp_idx,    cliqwp)
   call pbuf_get_field(pbuf, icswp_idx,    csnowp)
   call pbuf_get_field(pbuf, rel_idx,      rel)
   call pbuf_get_field(pbuf, rei_idx,      rei)
   call pbuf_get_field(pbuf, landm_idx,    landm)
   call pbuf_get_field(pbuf, dei_idx,      dei)
   call pbuf_get_field(pbuf, mu_idx,       mu)
   call pbuf_get_field(pbuf, lambdac_idx,  lambdac)
   call pbuf_get_field(pbuf, des_idx,      des)
   call pbuf_get_field(pbuf, crm_qrad_idx, crm_qrad)

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   endif

   if (prog_modal_aero) then
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
   endif

   do m=1,crm_nz
      k = pver-m+1
      do i = 1,ncol
        ! for energy conservation
        crm_qrad(i,:,:,m) = (rad_avgdata%qrs_crm(i,:,:,m)+rad_avgdata%qrl_crm(i,:,:,m)) * state%pdel(i,k)
      end do
   end do

   ! output rad inputs and resulting heating rates
   call rad_data_write(  pbuf, state, cam_in, coszrs )

   ! Compute net radiative heating tendency
   call radheat_tend(state, pbuf,  ptend, rad_avgdata%qrl_m(:,:), rad_avgdata%qrs_m(:,:), rad_avgdata%fsns_m(:), &
                     rad_avgdata%fsnt_m(:), rad_avgdata%flns_m(:), rad_avgdata%flnt_m(:), cam_in%asdir, net_flx)

   ! Compute heating rate for dtheta/dt
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (rad_avgdata%qrs_m(i,k) + rad_avgdata%qrl_m(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )

   ! convert radiative heating rates to Q*dp for energy conservation
   call pbuf_get_field(pbuf, qrs_idx, qrs)
   call pbuf_get_field(pbuf, qrl_idx, qrl)
   do k =1 , pver
      do i = 1, ncol
         qrs(i,k) = rad_avgdata%qrs_m(i,k)*state%pdel(i,k)
         qrl(i,k) = rad_avgdata%qrl_m(i,k)*state%pdel(i,k)
      end do
   end do

   ! Output icall=0 (climate)
   cam_out%flwds(:ncol)  = rad_avgdata%flwds_m(:ncol)
   cam_out%netsw(:ncol)  = rad_avgdata%fsns_m(:ncol)
   cam_out%sols(:ncol)   = rad_avgdata%sols_m(:ncol)
   cam_out%soll(:ncol)   = rad_avgdata%soll_m(:ncol)
   cam_out%solsd(:ncol)  = rad_avgdata%solsd_m(:ncol)
   cam_out%solld(:ncol)  = rad_avgdata%solld_m(:ncol)


    call pbuf_get_field(pbuf, fsns_idx, fsns)
    call pbuf_get_field(pbuf, fsnt_idx, fsnt)
    call pbuf_get_field(pbuf, flns_idx, flns)
    call pbuf_get_field(pbuf, flnt_idx, flnt)
    call pbuf_get_field(pbuf, fsds_idx, fsds)
    fsns(:ncol) = rad_avgdata%fsns_m(:ncol)
    fsnt(:ncol) = rad_avgdata%fsnt_m(:ncol)
    flns(:ncol) = rad_avgdata%flns_m(:ncol)
    flnt(:ncol) = rad_avgdata%flnt_m(:ncol)
    fsds(:ncol) = rad_avgdata%fsds_m(:ncol)

   deallocate(rad_avgdata%solin_m)
   deallocate(rad_avgdata%fsntoa_m)
   deallocate(rad_avgdata%fsutoa_m)
   deallocate(rad_avgdata%fsntoac_m)
   deallocate(rad_avgdata%fsnirt_m)
   deallocate(rad_avgdata%fsnrtc_m)
   deallocate(rad_avgdata%fsnirtsq_m)
   deallocate(rad_avgdata%fsntc_m)
   deallocate(rad_avgdata%fsnsc_m)
   deallocate(rad_avgdata%fsdsc_m)
   deallocate(rad_avgdata%flut_m)
   deallocate(rad_avgdata%flutc_m)
   deallocate(rad_avgdata%flntc_m)
   deallocate(rad_avgdata%flnsc_m)
   deallocate(rad_avgdata%fldsc_m)
   deallocate(rad_avgdata%flwds_m)
   deallocate(rad_avgdata%fsns_m)
   deallocate(rad_avgdata%fsnr_m)
   deallocate(rad_avgdata%fsnt_m)
   deallocate(rad_avgdata%flns_m)
   deallocate(rad_avgdata%flnt_m)
   deallocate(rad_avgdata%flnr_m)
   deallocate(rad_avgdata%fsds_m)
   deallocate(rad_avgdata%fln200_m)
   deallocate(rad_avgdata%fln200c_m)
   deallocate(rad_avgdata%fsn200_m)
   deallocate(rad_avgdata%fsn200c_m)
   deallocate(rad_avgdata%sols_m)
   deallocate(rad_avgdata%soll_m)
   deallocate(rad_avgdata%solsd_m)
   deallocate(rad_avgdata%solld_m)
   deallocate(rad_avgdata%qrs_m)
   deallocate(rad_avgdata%qrl_m)
   deallocate(rad_avgdata%qrsc_m)
   deallocate(rad_avgdata%qrlc_m)
   deallocate(rad_avgdata%rel_crm)
   deallocate(rad_avgdata%rei_crm)
   deallocate(rad_avgdata%cld_tau_crm)
   deallocate(rad_avgdata%qrl_crm)
   deallocate(rad_avgdata%qrs_crm)
   deallocate(rad_avgdata%crm_fsnt)
   deallocate(rad_avgdata%crm_fsntc)
   deallocate(rad_avgdata%crm_fsns)
   deallocate(rad_avgdata%crm_fsnsc)
   deallocate(rad_avgdata%crm_flnt)
   deallocate(rad_avgdata%crm_flntc)
   deallocate(rad_avgdata%crm_flns)
   deallocate(rad_avgdata%crm_flnsc)
   deallocate(rad_avgdata%crm_swcf)
   deallocate(rad_avgdata%crm_aodvisz)
   deallocate(rad_avgdata%crm_aodvis)
   deallocate(rad_avgdata%crm_aod400)
   deallocate(rad_avgdata%crm_aod700)
   deallocate(rad_avgdata%aod400)
   deallocate(rad_avgdata%aod700)

   deallocate(rad_avgdata%tot_cld_vistau_m)
   deallocate(rad_avgdata%tot_icld_vistau_m)
   deallocate(rad_avgdata%liq_icld_vistau_m)
   deallocate(rad_avgdata%ice_icld_vistau_m)
   deallocate(rad_avgdata%nct_tot_icld_vistau_m)

   deallocate(rad_avgdata%nct_liq_icld_vistau_m)
   deallocate(rad_avgdata%nct_ice_icld_vistau_m)
   deallocate(rad_avgdata%snow_icld_vistau_m)
   deallocate(rad_avgdata%nct_snow_icld_vistau_m)

   deallocate(rad_avgdata%dei_crm)
   deallocate(rad_avgdata%mu_crm)
   deallocate(rad_avgdata%lambdac_crm)
   deallocate(rad_avgdata%des_crm)

#endif
end subroutine spcam_radiation_finalize_m2005

!===============================================================================

subroutine spcam_radiation_col_finalize_m2005(state, ii, jj, pbuf, rd, cam_out, rad_avgdata)

    use physconst,       only: cpair
    use physics_buffer,  only: pbuf_old_tim_idx
    use radiation,       only: radiation_do
    use cam_history,     only: hist_fld_active

    type(physics_state),               intent(in) :: state
    integer,                           intent(in) :: ii
    integer,                           intent(in) :: jj
    type(physics_buffer_desc), pointer            :: pbuf(:)
    type(rad_out_t),     intent(in)         :: rd
    type(cam_out_t),     intent(inout)      :: cam_out

    type(rad_avgdata_type_m2005),    intent(inout) :: rad_avgdata

#ifdef m2005

    real(r8), parameter :: cgs2mks = 1.e-3_r8
    real(r8), parameter :: factor_xy = 1._r8/dble(crm_nx*crm_ny)

    integer :: i, k, m
    integer :: ncol
    integer :: itim_old

    logical :: dosw, dolw

    real(r8), pointer, dimension(:,:) :: qrs, qrl, cld
    real(r8), pointer, dimension(:)   :: fsds, fsns, fsnt, flns, flnt

    ncol  = state%ncol

    dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
    dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,  cld,      start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, qrs_idx,  qrs)
    call pbuf_get_field(pbuf, qrl_idx,  qrl)

    call pbuf_get_field(pbuf, fsns_idx, fsns)
    call pbuf_get_field(pbuf, fsnt_idx, fsnt)
    call pbuf_get_field(pbuf, flns_idx, flns)
    call pbuf_get_field(pbuf, flnt_idx, flnt)
    call pbuf_get_field(pbuf, fsds_idx, fsds)

    ! convert radiative heating rates from Q*dp for energy conservation
    do k =1 , pver
       do i = 1, ncol
          qrs(i,k) = qrs(i,k)/state%pdel(i,k)
          qrl(i,k) = qrl(i,k)/state%pdel(i,k)
       end do
    end do

    do m=1,crm_nz
       k = pver-m+1
       do i=1,ncol
          rad_avgdata%cld_tau_crm(i,ii,jj,m)=  rd%cld_tau_cloudsim(i,k)
       end do ! i
    end do ! m

    if (dosw) then

       do i=1, ncol
          rad_avgdata%qrs_m(i,:pver)  = rad_avgdata%qrs_m(i,:pver)  + qrs(i,:pver)    *factor_xy
          rad_avgdata%fsds_m(i)       = rad_avgdata%fsds_m(i)       + fsds(i)         *factor_xy
          rad_avgdata%fsnt_m(i)       = rad_avgdata%fsnt_m(i)       + fsnt(i)         *factor_xy
          rad_avgdata%fsns_m(i)       = rad_avgdata%fsns_m(i)       + fsns(i)         *factor_xy
          rad_avgdata%qrsc_m(i,:pver) = rad_avgdata%qrsc_m(i,:pver) + rd%qrsc(i,:pver)   *factor_xy
          rad_avgdata%solin_m(i)      = rad_avgdata%solin_m(i)      + rd%solin(i)        *factor_xy
          rad_avgdata%fsnirt_m(i)     = rad_avgdata%fsnirt_m(i)     + rd%fsnirt(i)       *factor_xy
          rad_avgdata%fsnrtc_m(i)     = rad_avgdata%fsnrtc_m(i)     + rd%fsnrtc(i)       *factor_xy
          rad_avgdata%fsnirtsq_m(i)   = rad_avgdata%fsnirtsq_m(i)   + rd%fsnirtsq(i)     *factor_xy
          rad_avgdata%fsntc_m(i)      = rad_avgdata%fsntc_m(i)      + rd%fsntc(i)        *factor_xy
          rad_avgdata%fsnsc_m(i)      = rad_avgdata%fsnsc_m(i)      + rd%fsnsc(i)        *factor_xy
          rad_avgdata%fsdsc_m(i)      = rad_avgdata%fsdsc_m(i)      + rd%fsdsc(i)        *factor_xy
          rad_avgdata%fsntoa_m(i)     = rad_avgdata%fsntoa_m(i)     + rd%fsntoa(i)       *factor_xy
          rad_avgdata%fsutoa_m(i)     = rad_avgdata%fsutoa_m(i)     + rd%fsutoa(i)       *factor_xy
          rad_avgdata%fsntoac_m(i)    = rad_avgdata%fsntoac_m(i)    + rd%fsntoac(i)      *factor_xy
          rad_avgdata%sols_m(i)       = rad_avgdata%sols_m(i)       + cam_out%sols(i)         *factor_xy
          rad_avgdata%soll_m(i)       = rad_avgdata%soll_m(i)       + cam_out%soll(i)         *factor_xy
          rad_avgdata%solsd_m(i)      = rad_avgdata%solsd_m(i)      + cam_out%solsd(i)        *factor_xy
          rad_avgdata%solld_m(i)      = rad_avgdata%solld_m(i)      + cam_out%solld(i)        *factor_xy
          rad_avgdata%fsn200_m(i)     = rad_avgdata%fsn200_m(i)     + rd%fsn200(i)       *factor_xy
          rad_avgdata%fsn200c_m(i)    = rad_avgdata%fsn200c_m(i)    + rd%fsn200c(i)      *factor_xy
          if (hist_fld_active('FSNR')) then
             rad_avgdata%fsnr_m(i)       = rad_avgdata%fsnr_m(i)       + rd%fsnr(i)         *factor_xy
          end if
          rad_avgdata%crm_fsnt(i, ii, jj) = fsnt(i)
          rad_avgdata%crm_fsntc(i,ii,jj)  = rd%fsntc(i)
          rad_avgdata%crm_fsns(i, ii, jj) = fsns(i)
          rad_avgdata%crm_fsnsc(i,ii,jj)  = rd%fsnsc(i)
          rad_avgdata%crm_swcf(i,ii,jj)   = rd%fsntoa(i) - rd%fsntoac(i)
          rad_avgdata%crm_aodvis(i,ii,jj) = sum(rd%aer_tau550(i, :))
          rad_avgdata%crm_aod400(i,ii,jj) = sum(rd%aer_tau400(i, :))
          rad_avgdata%crm_aod700(i,ii,jj) = sum(rd%aer_tau700(i, :))
          rad_avgdata%aod400(i) = rad_avgdata%aod400(i)+rad_avgdata%crm_aod400(i,ii,jj) * factor_xy
          rad_avgdata%aod700(i) = rad_avgdata%aod700(i)+rad_avgdata%crm_aod700(i,ii,jj) * factor_xy
       end do
       do m=1,crm_nz
          k = pver-m+1
          rad_avgdata%qrs_crm(:ncol,ii,jj,m) = qrs(:ncol,k) / cpair
          rad_avgdata%crm_aodvisz(:ncol, ii, jj, m) = rd%aer_tau550(:ncol,k)
       end do

       do i=1, ncol
          do k=1, pver
             if(rd%tot_icld_vistau(i,k).gt.1.0e-10_r8) then
                 rad_avgdata%tot_icld_vistau_m(i,k)      = rad_avgdata%tot_icld_vistau_m(i,k) + &
                                                           rd%tot_icld_vistau(i,k)*cld(i,k)
                 rad_avgdata%nct_tot_icld_vistau_m(i,k)  = rad_avgdata%nct_tot_icld_vistau_m(i,k) + cld(i,k)
             end if
             if(rd%liq_icld_vistau(i,k).gt.1.0e-10_r8) then
                 rad_avgdata%liq_icld_vistau_m(i,k)     = rad_avgdata%liq_icld_vistau_m(i,k) + &
                                                                       rd%liq_icld_vistau(i,k)*cld(i,k)
                 rad_avgdata%nct_liq_icld_vistau_m(i,k) = rad_avgdata%nct_liq_icld_vistau_m(i,k) + cld(i,k)
             end if
             if(rd%ice_icld_vistau(i,k).gt.1.0e-10_r8) then
                 rad_avgdata%ice_icld_vistau_m(i,k)     = rad_avgdata%ice_icld_vistau_m(i,k) + &
                                                                       rd%ice_icld_vistau(i,k)*cld(i,k)
                 rad_avgdata%nct_ice_icld_vistau_m(i,k) = rad_avgdata%nct_ice_icld_vistau_m(i,k) + cld(i,k)
             end if
             if(rd%snow_icld_vistau(i,k).gt.1.0e-10_r8) then
                rad_avgdata%snow_icld_vistau_m(i,k) = rad_avgdata%snow_icld_vistau_m(i,k) + &
                                                                       rd%snow_icld_vistau(i,k)
                rad_avgdata%nct_snow_icld_vistau_m(i,k) = rad_avgdata%nct_snow_icld_vistau_m(i,k) + 1
             end if
          end do
       end do
    end if   ! dosw

    if (dolw) then

       do i=1, ncol
          rad_avgdata%qrl_m(i,:pver) = rad_avgdata%qrl_m(i,:pver) + qrl(i,:pver)*factor_xy
          rad_avgdata%qrlc_m(i,:pver) = rad_avgdata%qrlc_m(i,:pver) + rd%qrlc(i,:pver)*factor_xy
          rad_avgdata%flnt_m(i)      = rad_avgdata%flnt_m(i)      + flnt(i)     *factor_xy
          rad_avgdata%flut_m(i)  = rad_avgdata%flut_m(i)+rd%flut(i) *factor_xy
          rad_avgdata%flutc_m(i)  = rad_avgdata%flutc_m(i)+rd%flutc(i) *factor_xy
          rad_avgdata%flntc_m(i)  = rad_avgdata%flntc_m(i)+rd%flntc(i) *factor_xy
          rad_avgdata%flns_m(i)      = rad_avgdata%flns_m(i)      + flns(i)     *factor_xy
          rad_avgdata%flnsc_m(i)  = rad_avgdata%flnsc_m(i)+rd%flnsc(i) *factor_xy
          rad_avgdata%fldsc_m(i)  = rad_avgdata%fldsc_m(i)+rd%fldsc(i) *factor_xy
          rad_avgdata%flwds_m(i)  = rad_avgdata%flwds_m(i)+cam_out%flwds(i) *factor_xy
          rad_avgdata%fln200_m(i) = rad_avgdata%fln200_m(i)+rd%fln200(i) *factor_xy
          rad_avgdata%fln200c_m(i) = rad_avgdata%fln200c_m(i)+rd%fln200c(i) *factor_xy
          if (hist_fld_active('FLNR')) then
             rad_avgdata%flnr_m(i) = rad_avgdata%flnr_m(i)+rd%flnr(i) *factor_xy
          end if

          call pbuf_get_field(pbuf, fsns_idx, fsns)
          call pbuf_get_field(pbuf, fsnt_idx, fsnt)
          call pbuf_get_field(pbuf, flns_idx, flns)
          call pbuf_get_field(pbuf, flnt_idx, flnt)
          call pbuf_get_field(pbuf, fsds_idx, fsds)

          rad_avgdata%crm_flnt(i, ii, jj) = flnt(i)
          rad_avgdata%crm_flntc(i,ii,jj)  = rd%flntc(i)
          rad_avgdata%crm_flns(i, ii, jj) = flns(i)
          rad_avgdata%crm_flnsc(i,ii,jj)  = rd%flnsc(i)
          do m=1,crm_nz
             k = pver-m+1
                   rad_avgdata%qrl_crm(:ncol,ii,jj,m) = qrl(:ncol,k) / cpair
          end do

       end do

    end if  !dolw


#endif

end subroutine spcam_radiation_col_finalize_m2005

!===============================================================================

subroutine spcam_radiation_setup_sam1mom(cam_in, cldn, state, pbuf, rad_avgdata, state_loc)

    use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
    use physics_buffer,  only: pbuf_old_tim_idx

    type(cam_in_t),                  intent(in)             :: cam_in
    real(r8), dimension(:,:),        intent(out)            :: cldn
    type(physics_state),             intent(in)             :: state
    type(physics_buffer_desc),       intent(inout), pointer :: pbuf(:)

    type(rad_avgdata_type_sam1mom) :: rad_avgdata
    type(physics_state),             intent(inout)          :: state_loc

#ifdef sam1mom
    real(r8),pointer :: emis(:,:)              ! Cloud longwave emissivity
    real(r8),pointer :: cldtau(:,:)            ! Cloud longwave optical depth
    real(r8),pointer :: cicewp(:,:)            ! in-cloud cloud ice water path
    real(r8),pointer :: cliqwp(:,:)            ! in-cloud cloud liquid water path

    real(r8), pointer, dimension(:,:)      :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:)      :: rei     ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:)      :: cld
    real(r8), pointer, dimension(:)        :: landm   ! land fraction ramp


    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns
    integer :: itim_old

    ncol  = state%ncol
    lchnk = state%lchnk


    call physics_state_copy(state, state_loc)


    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ! Save the grid level cld values as cld will be overwritten with each crm-scale level value during radiation
    cldn = cld

    allocate(rad_avgdata%solin_m      (pcols))
    allocate(rad_avgdata%fsntoa_m     (pcols))
    allocate(rad_avgdata%fsutoa_m     (pcols))
    allocate(rad_avgdata%fsntoac_m    (pcols))
    allocate(rad_avgdata%fsnirt_m     (pcols))
    allocate(rad_avgdata%fsnrtc_m     (pcols))
    allocate(rad_avgdata%fsnirtsq_m   (pcols))
    allocate(rad_avgdata%fsntc_m      (pcols))
    allocate(rad_avgdata%fsnsc_m      (pcols))
    allocate(rad_avgdata%fsdsc_m      (pcols))
    allocate(rad_avgdata%flut_m       (pcols))
    allocate(rad_avgdata%flutc_m      (pcols))
    allocate(rad_avgdata%flntc_m      (pcols))
    allocate(rad_avgdata%flnsc_m      (pcols))
    allocate(rad_avgdata%fldsc_m      (pcols))
    allocate(rad_avgdata%flwds_m      (pcols))
    allocate(rad_avgdata%fsns_m       (pcols))
    allocate(rad_avgdata%fsnr_m       (pcols))
    allocate(rad_avgdata%fsnt_m       (pcols))
    allocate(rad_avgdata%flns_m       (pcols))
    allocate(rad_avgdata%flnt_m       (pcols))
    allocate(rad_avgdata%flnr_m       (pcols))
    allocate(rad_avgdata%fsds_m       (pcols))
    allocate(rad_avgdata%fln200_m     (pcols))
    allocate(rad_avgdata%fln200c_m    (pcols))
    allocate(rad_avgdata%fsn200_m     (pcols))
    allocate(rad_avgdata%fsn200c_m    (pcols))
    allocate(rad_avgdata%sols_m       (pcols))
    allocate(rad_avgdata%soll_m       (pcols))
    allocate(rad_avgdata%solsd_m      (pcols))
    allocate(rad_avgdata%solld_m      (pcols))
    allocate(rad_avgdata%qrs_m        (pcols,pver))
    allocate(rad_avgdata%qrl_m        (pcols,pver))
    allocate(rad_avgdata%qrsc_m       (pcols,pver))
    allocate(rad_avgdata%qrlc_m       (pcols,pver))
    allocate(rad_avgdata%rel_crm      (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%rei_crm      (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%cld_tau_crm  (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%qrl_crm      (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%qrs_crm      (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%crm_fsnt     (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_fsntc    (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_fsns     (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_fsnsc    (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_flnt     (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_flntc    (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_flns     (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_flnsc    (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_swcf     (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_aodvisz  (pcols, crm_nx, crm_ny, crm_nz))
    allocate(rad_avgdata%crm_aodvis   (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_aod400   (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%crm_aod700   (pcols, crm_nx, crm_ny))
    allocate(rad_avgdata%aod400       (pcols))
    allocate(rad_avgdata%aod700       (pcols))
    allocate(rad_avgdata%fsdtoa_m     (pcols))
    allocate(rad_avgdata%flds_m       (pcols))

    allocate(rad_avgdata%tot_cld_vistau_m    ( pcols,pver))
    allocate(rad_avgdata%tot_icld_vistau_m    (pcols,pver))
    allocate(rad_avgdata%liq_icld_vistau_m    (pcols,pver))
    allocate(rad_avgdata%ice_icld_vistau_m    (pcols,pver))
    allocate(rad_avgdata%nct_tot_icld_vistau_m(pcols,pver))
    allocate(rad_avgdata%nct_liq_icld_vistau_m(pcols,pver))
    allocate(rad_avgdata%nct_ice_icld_vistau_m(pcols,pver))

    call pbuf_get_field(pbuf, rel_idx,    rel)
    call pbuf_get_field(pbuf, rei_idx,    rei)
    call pbuf_get_field(pbuf, landm_idx,  landm)
    call pbuf_get_field(pbuf, crm_t_rad_idx,   rad_avgdata%t_rad)
    call pbuf_get_field(pbuf, crm_qc_rad_idx,  rad_avgdata%qc_rad)
    call pbuf_get_field(pbuf, crm_qi_rad_idx,  rad_avgdata%qi_rad)
    call pbuf_get_field(pbuf, crm_qv_rad_idx,  rad_avgdata%qv_rad)
    call pbuf_get_field(pbuf, crm_qrad_idx,    rad_avgdata%crm_qrad)


    ! pbuf cloud properties set in cloud_diagnostics
    call pbuf_get_field(pbuf, cicewp_idx, cicewp)
    call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
    call pbuf_get_field(pbuf, cldemis_idx, emis)
    call pbuf_get_field(pbuf, cldtau_idx, cldtau)


    rad_avgdata%solin_m    = 0._r8
    rad_avgdata%fsntoa_m   = 0._r8
    rad_avgdata%fsutoa_m   = 0._r8
    rad_avgdata%fsntoac_m  = 0._r8
    rad_avgdata%fsnirt_m   = 0._r8
    rad_avgdata%fsnrtc_m   = 0._r8
    rad_avgdata%fsnirtsq_m = 0._r8
    rad_avgdata%fsntc_m    = 0._r8
    rad_avgdata%fsdtoa_m   = 0._r8
    rad_avgdata%fsnsc_m    = 0._r8
    rad_avgdata%fsdsc_m    = 0._r8
    rad_avgdata%flut_m     = 0._r8
    rad_avgdata%flutc_m    = 0._r8
    rad_avgdata%flntc_m    = 0._r8
    rad_avgdata%flnsc_m    = 0._r8
    rad_avgdata%flds_m     = 0._r8
    rad_avgdata%fldsc_m    = 0._r8
    rad_avgdata%fsns_m     = 0._r8
    rad_avgdata%fsnt_m     = 0._r8
    rad_avgdata%flns_m     = 0._r8
    rad_avgdata%flnt_m     = 0._r8
    rad_avgdata%flnr_m     = 0._r8
    rad_avgdata%fsds_m     = 0._r8
    rad_avgdata%fsnr_m     = 0._r8
    rad_avgdata%fln200_m   = 0._r8
    rad_avgdata%fln200c_m  = 0._r8
    rad_avgdata%fsn200_m   = 0._r8
    rad_avgdata%fsn200c_m  = 0._r8
    rad_avgdata%sols_m     = 0._r8
    rad_avgdata%soll_m     = 0._r8
    rad_avgdata%solsd_m    = 0._r8
    rad_avgdata%solld_m    = 0._r8
    rad_avgdata%qrs_m      = 0._r8
    rad_avgdata%qrl_m      = 0._r8
    rad_avgdata%qrsc_m     = 0._r8
    rad_avgdata%qrlc_m     = 0._r8
    rad_avgdata%qrs_crm    = 0._r8
    rad_avgdata%qrl_crm    = 0._r8

    rad_avgdata%tot_cld_vistau_m =0._r8
    rad_avgdata%tot_icld_vistau_m=0._r8  ; rad_avgdata%nct_tot_icld_vistau_m=0._r8
    rad_avgdata%liq_icld_vistau_m=0._r8  ; rad_avgdata%nct_liq_icld_vistau_m=0._r8
    rad_avgdata%ice_icld_vistau_m=0._r8  ; rad_avgdata%nct_ice_icld_vistau_m=0._r8


    ! Compute effective sizes
    call cldefr(lchnk, ncol, cam_in%landfrac, state%t, rel, rei, state%ps, state%pmid, landm, cam_in%icefrac, cam_in%snowhland)

    cicewp(1:ncol,1:pver) = 0._r8
    cliqwp(1:ncol,1:pver) = 0._r8

#endif
end subroutine spcam_radiation_setup_sam1mom

!===============================================================================

subroutine spcam_radiation_col_setup_sam1mom(ii, jj, state_loc, pbuf, rad_avgdata)

    use physics_buffer, only: pbuf_old_tim_idx
    use physconst,      only: gravit

    integer,intent(in) :: ii,jj

    type(physics_state),       intent(inout)          :: state_loc
    type(physics_buffer_desc), intent(inout), pointer :: pbuf(:)
    type(rad_avgdata_type_sam1mom),    intent(inout)          :: rad_avgdata

#ifdef sam1mom

    real(r8),pointer :: emis(:,:)              ! Cloud longwave emissivity
    real(r8),pointer :: cldtau(:,:)            ! Cloud longwave optical depth
    real(r8),pointer :: cicewp(:,:)            ! in-cloud cloud ice water path
    real(r8),pointer :: cliqwp(:,:)            ! in-cloud cloud liquid water path

    real(r8), pointer, dimension(:,:)     :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:)     :: rei     ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:,:,:) :: cld_rad ! rad cloud fraction
    real(r8), pointer, dimension(:,:)     :: pmxrgn  ! Maximum values of pressure for each
						     !    maximally overlapped region.
						     !    0->pmxrgn(i,1) is range of pressure for
						     !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
						     !    2nd region, etc
    integer, pointer, dimension(:)  :: nmxrgn        ! pbuf pointer to Number of maximally overlapped regions

    real(r8)                          :: qtot
    real(r8), dimension(pcols,pver)   :: fice
    real(r8), dimension(pcols,pver)   :: tmp
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

    integer :: itim_old
    integer :: m, k, i

    integer :: lchnk                              ! chunk identifier
    integer :: ncol                               ! number of atmospheric columns

    lchnk = state_loc%lchnk
    ncol  = state_loc%ncol

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,         cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, pmxrgn_idx,      pmxrgn)
    call pbuf_get_field(pbuf, nmxrgn_idx,      nmxrgn)
    call pbuf_get_field(pbuf, rel_idx,         rel)
    call pbuf_get_field(pbuf, rei_idx,         rei)
    call pbuf_get_field(pbuf, crm_cld_rad_idx, cld_rad)

    call pbuf_get_field(pbuf, crm_t_rad_idx,  rad_avgdata%t_rad)
    call pbuf_get_field(pbuf, crm_qc_rad_idx, rad_avgdata%qc_rad)
    call pbuf_get_field(pbuf, crm_qi_rad_idx, rad_avgdata%qi_rad)
    call pbuf_get_field(pbuf, crm_qv_rad_idx, rad_avgdata%qv_rad)
    call pbuf_get_field(pbuf, crm_qrad_idx,   rad_avgdata%crm_qrad)


    ! pbuf cloud properties set in cloud_diagnostics
    call pbuf_get_field(pbuf, cicewp_idx, cicewp)
    call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
    call pbuf_get_field(pbuf, cldemis_idx, emis)
    call pbuf_get_field(pbuf, cldtau_idx, cldtau)

    fice(1:ncol,1:pver-crm_nz) = 0._r8

    do m=1,crm_nz
       k = pver-m+1
       do i=1,ncol

	  qtot = rad_avgdata%qc_rad(i,ii,jj,m) + rad_avgdata%qi_rad(i,ii,jj,m)
	  if(qtot.gt.1.e-9_r8) then
	    fice(i,k) = rad_avgdata%qi_rad(i,ii,jj,m)/qtot
	    ! In case CRM produces fractional cloudiness
	    cld(i,k) = min(0.99_r8, cld_rad(i,ii,jj,m))

	    cicewp(i,k) = rad_avgdata%qi_rad(i,ii,jj,m)*state_loc%pdel(i,k)/gravit*1000.0_r8 &
		    / max(0.01_r8,cld(i,k)) ! In-cloud ice water path.
	    cliqwp(i,k) = rad_avgdata%qc_rad(i,ii,jj,m)*state_loc%pdel(i,k)/gravit*1000.0_r8 &
		    / max(0.01_r8,cld(i,k)) ! In-cloud liquid water path.
	  else
	    fice(i,k)=0._r8
	    cld(i,k)=0._r8
	    cicewp(i,k) = 0._r8           ! In-cloud ice water path.
	    cliqwp(i,k) = 0._r8           ! In-cloud liquid water path.
	  end if
       end do ! i
    end do ! m

    !  Cloud emissivity.

    tmp(:ncol,:) = cicewp(:ncol,:) + cliqwp(:ncol,:)
    call cldems(lchnk, ncol, tmp, fice, rei, emis, cldtau)

    call cldovrlap(lchnk, ncol, state_loc%pint, cld, nmxrgn, pmxrgn)

    ! Setup the trad and qvrad variables (now in state)
    do m=1,crm_nz
       k = pver-m+1
       do i=1,ncol
          state_loc%q(i,k,1) = max(1.e-9_r8,rad_avgdata%qv_rad(i,ii,jj,m))
          state_loc%t(i,k)  = rad_avgdata%t_rad(i,ii,jj,m)
       end do
    end do


#endif
end subroutine spcam_radiation_col_setup_sam1mom

!===============================================================================

subroutine spcam_radiation_finalize_sam1mom(cam_in, state, pbuf, rad_avgdata, cam_out, cldn, net_flx, ptend)

    use physconst,       only: cpair
    use rad_constituents,only: rad_cnst_out

    use physconst,       only: cappa
    use radiation_data,  only: rad_data_write
    use radheat,         only: radheat_tend
    use time_manager,    only: get_curr_calday
    use physics_buffer,  only: pbuf_old_tim_idx
    use orbit,           only: zenith

    type(cam_in_t),                    intent(in) :: cam_in
    type(physics_state),               intent(in) :: state


    type(physics_buffer_desc),         intent(inout), pointer      :: pbuf(:)
    type(rad_avgdata_type_sam1mom),    intent(inout) :: rad_avgdata
    type(cam_out_t),                   intent(inout) :: cam_out
    real(r8), dimension(:,:),          intent(in)    :: cldn
    real(r8),                          intent(inout) :: net_flx(pcols)

    type(physics_ptend),               intent(out)   :: ptend            ! indivdual parameterization tendencies

#ifdef sam1mom

    integer  :: lchnk                    ! chunk identifier
    integer  :: ncol                     ! number of atmospheric columns
    integer  :: i, k, m
    real(r8) :: ftem(pcols,pver)         ! Temporary workspace for outfld variables

    real(r8), pointer, dimension(:,:) :: qrs, qrl, cld
    real(r8), pointer :: fsns(:)      ! Surface solar absorbed flux
    real(r8), pointer :: fsnt(:)      ! Net column abs solar flux at model top
    real(r8), pointer :: flns(:)      ! Srf longwave cooling (up-down) flux
    real(r8), pointer :: flnt(:)      ! Net outgoing lw flux at model top
    real(r8), pointer :: fsds(:)      ! Surface solar down flux



    real(r8) :: calday                   ! current calendar day
    real(r8) :: clat(pcols)              ! current latitudes(radians)
    real(r8) :: clon(pcols)              ! current longitudes(radians)
    real(r8) :: coszrs(pcols)            ! Cosine solar zenith angle
    real(r8) :: factor_xy

    integer :: Nday                      ! Number of daylight columns
    integer :: Nnite                     ! Number of night columns
    integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
    integer, dimension(pcols) :: IdxNite ! Indicies of night coumns
    integer :: itim_old

    lchnk = state%lchnk
    ncol  = state%ncol

    call pbuf_get_field(pbuf, qrs_idx, qrs)
    call pbuf_get_field(pbuf, qrl_idx, qrl)

    factor_xy = 1._r8/dble(crm_nx*crm_ny)

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    ! Reassign the grid level cld values since cld was overwritten with each crm-scale level value during radiation
    cld = cldn


    do m=1,crm_nz
       k = pver-m+1
       do i = 1,ncol
         ! for energy conservation
         rad_avgdata%crm_qrad(i,:,:,m) = (rad_avgdata%qrs_crm(i,:,:,m)+rad_avgdata%qrl_crm(i,:,:,m)) * state%pdel(i,k)
       end do
    end do

    call pbuf_get_field(pbuf, fsns_idx, fsns)
    call pbuf_get_field(pbuf, fsnt_idx, fsnt)
    call pbuf_get_field(pbuf, flns_idx, flns)
    call pbuf_get_field(pbuf, flnt_idx, flnt)
    call pbuf_get_field(pbuf, fsds_idx, fsds)

    fsns = rad_avgdata%fsns_m(:)
    fsnt = rad_avgdata%fsnt_m(:)
    flns = rad_avgdata%flns_m(:)
    flnt = rad_avgdata%flnt_m(:)
    fsds = rad_avgdata%fsds_m(:)

    calday = get_curr_calday()

    ! Cosine solar zenith angle for current time step
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)

    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
	  Nday = Nday + 1
	  IdxDay(Nday) = i
       else
	  Nnite = Nnite + 1
	  IdxNite(Nnite) = i
       end if
    end do

    cam_out%sols(:ncol) = rad_avgdata%sols_m(:ncol)
    cam_out%soll(:ncol) = rad_avgdata%soll_m(:ncol)
    cam_out%solsd(:ncol) = rad_avgdata%solsd_m(:ncol)
    cam_out%solld(:ncol) = rad_avgdata%solld_m(:ncol)

    call outfld('CRM_QRS ',rad_avgdata%qrs_crm,pcols,lchnk)
    call outfld('QRS     ',rad_avgdata%qrs_m(:,:)/cpair ,pcols,lchnk)
    call outfld('QRSC    ',rad_avgdata%qrsc_m/cpair,pcols,lchnk)
    call outfld('SOLIN   ',rad_avgdata%solin_m(:) ,pcols,lchnk)
    call outfld('FSDS    ',rad_avgdata%fsds_m(:)  ,pcols,lchnk)
    call outfld('FSNIRTOA',rad_avgdata%fsnirt_m(:),pcols,lchnk)
    call outfld('FSNRTOAC',rad_avgdata%fsnrtc_m(:),pcols,lchnk)
    call outfld('FSNRTOAS',rad_avgdata%fsnirtsq_m(:),pcols,lchnk)
    call outfld('FSNT    ',rad_avgdata%fsnt_m(:)  ,pcols,lchnk)
    call outfld('FSDTOA  ',rad_avgdata%fsdtoa_m(:),pcols,lchnk)
    call outfld('FSNS    ',rad_avgdata%fsns_m(:)  ,pcols,lchnk)
    call outfld('FSNTC   ',rad_avgdata%fsntc_m(:) ,pcols,lchnk)
    call outfld('FSNSC   ',rad_avgdata%fsnsc_m(:) ,pcols,lchnk)
    call outfld('FSDSC   ',rad_avgdata%fsdsc_m(:) ,pcols,lchnk)
    call outfld('FSNTOA  ',rad_avgdata%fsntoa_m(:),pcols,lchnk)
    call outfld('FSUTOA  ',rad_avgdata%fsutoa_m(:),pcols,lchnk)
    call outfld('FSNTOAC ',rad_avgdata%fsntoac_m(:),pcols,lchnk)
    call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
    call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
    call outfld('SOLSD   ',cam_out%solsd ,pcols,lchnk)
    call outfld('SOLLD   ',cam_out%solld ,pcols,lchnk)
    call outfld('FSN200  ',rad_avgdata%fsn200_m(:),pcols,lchnk)
    call outfld('FSN200C ',rad_avgdata%fsn200c_m(:),pcols,lchnk)
    call outfld('FSNR'    ,rad_avgdata%fsnr_m(:)  ,pcols,lchnk)
    call outfld('SWCF    ',rad_avgdata%fsntoa_m(:ncol)-rad_avgdata%fsntoac_m(:ncol)  ,pcols,lchnk)

    do i=1, Nday
       do k=1, pver
	  rad_avgdata%tot_cld_vistau_m(IdxDay(i),k) = rad_avgdata%tot_icld_vistau_m(IdxDay(i),k) *  factor_xy
	  if(rad_avgdata%nct_tot_icld_vistau_m(IdxDay(i),k).ge. 0.1_r8) then
	    rad_avgdata%tot_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%tot_icld_vistau_m(IdxDay(i),k)/&
						   rad_avgdata%nct_tot_icld_vistau_m(IdxDay(i),k)
	  else
	    rad_avgdata%tot_icld_vistau_m(IdxDay(i),k)  = 0.0_r8
	  end if
	  if(rad_avgdata%nct_liq_icld_vistau_m(IdxDay(i),k).ge. 0.1_r8) then
	     rad_avgdata%liq_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%liq_icld_vistau_m(IdxDay(i),k)/&
						    rad_avgdata%nct_liq_icld_vistau_m(IdxDay(i),k)
	  else
	      rad_avgdata%liq_icld_vistau_m(IdxDay(i),k)  = 0.0_r8
	  end if
	  if(rad_avgdata%nct_ice_icld_vistau_m(IdxDay(i),k).ge. 0.1_r8) then
	     rad_avgdata%ice_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%ice_icld_vistau_m(IdxDay(i),k)/&
						    rad_avgdata%nct_ice_icld_vistau_m(IdxDay(i),k)
	  else
	     rad_avgdata%ice_icld_vistau_m(IdxDay(i),k)  = 0.0_r8
	  end if
       end do
    end do

    ! add fillvalue for night columns
    do i = 1, Nnite
       rad_avgdata%tot_cld_vistau_m(IdxNite(i),:)   = fillvalue
       rad_avgdata%tot_icld_vistau_m(IdxNite(i),:)  = fillvalue
       rad_avgdata%liq_icld_vistau_m(IdxNite(i),:)  = fillvalue
       rad_avgdata%ice_icld_vistau_m(IdxNite(i),:)  = fillvalue
    end do

    call outfld ('TOT_CLD_VISTAU    ',rad_avgdata%tot_cld_vistau_m  ,pcols,lchnk)
    call outfld ('TOT_ICLD_VISTAU   ',rad_avgdata%tot_icld_vistau_m ,pcols,lchnk)
    call outfld ('LIQ_ICLD_VISTAU   ',rad_avgdata%liq_icld_vistau_m ,pcols,lchnk)
    call outfld ('ICE_ICLD_VISTAU   ',rad_avgdata%ice_icld_vistau_m ,pcols,lchnk)


    ! Longwave
    cam_out%flwds(:) = rad_avgdata%flds_m(:)
    call outfld('CRM_QRL ',rad_avgdata%qrl_crm,                     pcols, lchnk)
    call outfld('QRL     ',rad_avgdata%qrl_m(:ncol,:)/cpair,      ncol,  lchnk)
    call outfld('QRLC    ',rad_avgdata%qrlc_m(:ncol,:)/cpair,     ncol,  lchnk)
    call outfld('FLNT    ',rad_avgdata%flnt_m ,                     pcols, lchnk)
    call outfld('FLUT    ',rad_avgdata%flut_m,                      pcols, lchnk)
    call outfld('FLUTC   ',rad_avgdata%flutc_m,                     pcols, lchnk)
    call outfld('FLNTC   ',rad_avgdata%flntc_m,                     pcols, lchnk)
    call outfld('FLNS    ',rad_avgdata%flns_m,                      pcols, lchnk)
    call outfld('FLDS    ',rad_avgdata%flds_m,                      pcols, lchnk)
    call outfld('FLNSC   ',rad_avgdata%flnsc_m,                     pcols, lchnk)
    call outfld('FLDSC   ',rad_avgdata%fldsc_m,                     pcols, lchnk)
    call outfld('LWCF    ',rad_avgdata%flutc_m-rad_avgdata%flut_m,  pcols, lchnk)
    call outfld('FLN200  ',rad_avgdata%fln200_m,                    pcols, lchnk)
    call outfld('FLN200C ',rad_avgdata%fln200c_m,                   pcols, lchnk)
    call outfld('FLNR '   ,rad_avgdata%flnr_m,                      pcols, lchnk)

    ! Output aerosol mmr
    call rad_cnst_out(0, state, pbuf)

    ! output rad inputs and resulting heating rates
    call rad_data_write( pbuf, state, cam_in, coszrs )

    ! Compute net radiative heating tendency
    call radheat_tend(state, pbuf,  ptend, rad_avgdata%qrl_m, rad_avgdata%qrs_m, rad_avgdata%fsns_m, &
		      rad_avgdata%fsnt_m, rad_avgdata%flns_m, rad_avgdata%flnt_m, cam_in%asdir, net_flx)

    ! Compute heating rate for dtheta/dt
    do k=1,pver
       do i=1,ncol
	  ftem(i,k) = (rad_avgdata%qrs_m(i,k) + rad_avgdata%qrl_m(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
       end do
    end do
    call outfld('HR      ',ftem    ,pcols   ,lchnk   )

    do k =1 , pver
       do i = 1, ncol
          qrs(i,k) = rad_avgdata%qrs_m(i,k)*state%pdel(i,k)
          qrl(i,k) = rad_avgdata%qrl_m(i,k)*state%pdel(i,k)
       end do
    end do

    cam_out%netsw(:ncol) = rad_avgdata%fsns_m(:ncol)
    cam_out%flwds(:ncol) = rad_avgdata%flds_m(:ncol)

    deallocate(rad_avgdata%solin_m)
    deallocate(rad_avgdata%fsntoa_m)
    deallocate(rad_avgdata%fsutoa_m)
    deallocate(rad_avgdata%fsntoac_m)
    deallocate(rad_avgdata%fsnirt_m)
    deallocate(rad_avgdata%fsnrtc_m)
    deallocate(rad_avgdata%fsnirtsq_m)
    deallocate(rad_avgdata%fsntc_m)
    deallocate(rad_avgdata%fsnsc_m)
    deallocate(rad_avgdata%fsdsc_m)
    deallocate(rad_avgdata%flut_m)
    deallocate(rad_avgdata%flutc_m)
    deallocate(rad_avgdata%flntc_m)
    deallocate(rad_avgdata%flnsc_m)
    deallocate(rad_avgdata%fldsc_m)
    deallocate(rad_avgdata%flwds_m)
    deallocate(rad_avgdata%fsns_m)
    deallocate(rad_avgdata%fsnr_m)
    deallocate(rad_avgdata%fsnt_m)
    deallocate(rad_avgdata%flns_m)
    deallocate(rad_avgdata%flnt_m)
    deallocate(rad_avgdata%flnr_m)
    deallocate(rad_avgdata%fsds_m)
    deallocate(rad_avgdata%fln200_m)
    deallocate(rad_avgdata%fln200c_m)
    deallocate(rad_avgdata%fsn200_m)
    deallocate(rad_avgdata%fsn200c_m)
    deallocate(rad_avgdata%sols_m)
    deallocate(rad_avgdata%soll_m)
    deallocate(rad_avgdata%solsd_m)
    deallocate(rad_avgdata%solld_m)
    deallocate(rad_avgdata%qrs_m)
    deallocate(rad_avgdata%qrl_m)
    deallocate(rad_avgdata%qrsc_m)
    deallocate(rad_avgdata%qrlc_m)
    deallocate(rad_avgdata%rel_crm)
    deallocate(rad_avgdata%rei_crm)
    deallocate(rad_avgdata%cld_tau_crm)
    deallocate(rad_avgdata%qrl_crm)
    deallocate(rad_avgdata%qrs_crm)
    deallocate(rad_avgdata%crm_fsnt)
    deallocate(rad_avgdata%crm_fsntc)
    deallocate(rad_avgdata%crm_fsns)
    deallocate(rad_avgdata%crm_fsnsc)
    deallocate(rad_avgdata%crm_flnt)
    deallocate(rad_avgdata%crm_flntc)
    deallocate(rad_avgdata%crm_flns)
    deallocate(rad_avgdata%crm_flnsc)
    deallocate(rad_avgdata%crm_swcf)
    deallocate(rad_avgdata%crm_aodvisz)
    deallocate(rad_avgdata%crm_aodvis)
    deallocate(rad_avgdata%crm_aod400)
    deallocate(rad_avgdata%crm_aod700)
    deallocate(rad_avgdata%aod400)
    deallocate(rad_avgdata%aod700)
    deallocate(rad_avgdata%fsdtoa_m)
    deallocate(rad_avgdata%flds_m)

    deallocate(rad_avgdata%tot_cld_vistau_m)
    deallocate(rad_avgdata%tot_icld_vistau_m)
    deallocate(rad_avgdata%liq_icld_vistau_m)
    deallocate(rad_avgdata%ice_icld_vistau_m)
    deallocate(rad_avgdata%nct_tot_icld_vistau_m)
    deallocate(rad_avgdata%nct_liq_icld_vistau_m)
    deallocate(rad_avgdata%nct_ice_icld_vistau_m)
#endif

end subroutine spcam_radiation_finalize_sam1mom

subroutine spcam_radiation_col_finalize_sam1mom(state, ii, jj, pbuf, rd, cam_out, rad_avgdata)

    use physconst,       only: cpair
    use physics_buffer,  only: pbuf_old_tim_idx
    use orbit,           only: zenith
    use time_manager,    only: get_curr_calday
    use radiation,       only: radiation_do

    type(physics_state),               intent(in) :: state
    integer,                           intent(in) :: ii
    integer,                           intent(in) :: jj
    type(physics_buffer_desc), pointer            :: pbuf(:)
    type(rad_out_t),     intent(in)         :: rd
    type(cam_out_t),     intent(inout)      :: cam_out

    real(r8), parameter :: cgs2mks = 1.e-3_r8

    type(rad_avgdata_type_sam1mom),    intent(inout) :: rad_avgdata

#ifdef sam1mom

    real(r8), pointer :: fsns(:)      ! Surface solar absorbed flux
    real(r8), pointer :: fsnt(:)      ! Net column abs solar flux at model top
    real(r8), pointer :: flns(:)      ! Srf longwave cooling (up-down) flux
    real(r8), pointer :: flnt(:)      ! Net outgoing lw flux at model top
    real(r8), pointer :: fsds(:)      ! Surface solar down flux

    integer                           :: itim_old
    integer                           :: ncol
    integer                           :: i, m, k, lchnk


    logical :: dosw, dolw
    integer :: Nday                      ! Number of daylight columns
    integer :: Nnite                     ! Number of night columns
    integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns


    real(r8) :: calday                   ! current calendar day
    real(r8) :: clat(pcols)              ! current latitudes(radians)
    real(r8) :: clon(pcols)              ! current longitudes(radians)
    real(r8) :: coszrs(pcols)            ! Cosine solar zenith angle
    real(r8) :: factor_xy

    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: qrs
    real(r8), pointer, dimension(:,:) :: qrl

    ncol  = state%ncol
    lchnk = state%lchnk

    calday = get_curr_calday()

    ! Cosine solar zenith angle for current time step
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call zenith (calday, clat, clon, coszrs, ncol)

    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
          Nday = Nday + 1
          IdxDay(Nday) = i
       else
          Nnite = Nnite + 1
       end if
    end do

    dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
    dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?

    factor_xy = 1._r8/dble(crm_nx*crm_ny)

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, qrs_idx,qrs)
    call pbuf_get_field(pbuf, qrl_idx,qrl)
    call pbuf_get_field(pbuf, qrl_idx,qrl)

    ! convert radiative heating rates from Q*dp for energy conservation
    do k =1 , pver
       do i = 1, ncol
          qrs(i,k) = qrs(i,k)/state%pdel(i,k)
          qrl(i,k) = qrl(i,k)/state%pdel(i,k)
       end do
    end do

    if (dosw) then
       call pbuf_get_field(pbuf, fsds_idx, fsds)
       call pbuf_get_field(pbuf, fsns_idx, fsns)
       call pbuf_get_field(pbuf, fsnt_idx, fsnt)
       do i=1,ncol
          rad_avgdata%fsds_m(i)    = rad_avgdata%fsds_m(i)    +fsds(i) *factor_xy
          rad_avgdata%fsns_m(i)    = rad_avgdata%fsns_m(i)    +fsns(i) *factor_xy
          rad_avgdata%fsnt_m(i)    = rad_avgdata%fsnt_m(i)    +fsnt(i) *factor_xy

          rad_avgdata%solin_m(i)   = rad_avgdata%solin_m(i)   +rd%solin(i)*factor_xy
          rad_avgdata%fsnirt_m(i)  = rad_avgdata%fsnirt_m(i)  +rd%fsnirt(i)*factor_xy
          rad_avgdata%fsnrtc_m(i)  = rad_avgdata%fsnrtc_m(i)  +rd%fsnrtc(i)*factor_xy
          rad_avgdata%fsnirtsq_m(i)= rad_avgdata%fsnirtsq_m(i)+rd%fsnirtsq(i)*factor_xy
          rad_avgdata%fsdtoa_m(i)  = rad_avgdata%fsdtoa_m(i)  +rd%fsdtoa(i)*factor_xy
          rad_avgdata%fsntc_m(i)   = rad_avgdata%fsntc_m(i)   +rd%fsntc(i)*factor_xy
          rad_avgdata%fsnsc_m(i)   = rad_avgdata%fsnsc_m(i)   +rd%fsnsc(i)*factor_xy
          rad_avgdata%fsdsc_m(i)   = rad_avgdata%fsdsc_m(i)   +rd%fsdsc(i)*factor_xy
          rad_avgdata%fsntoa_m(i)  = rad_avgdata%fsntoa_m(i)  +rd%fsntoa(i)*factor_xy
          rad_avgdata%fsutoa_m(i)  = rad_avgdata%fsutoa_m(i)  +rd%fsutoa(i)*factor_xy
          rad_avgdata%fsntoac_m(i) = rad_avgdata%fsntoac_m(i) +rd%fsntoac(i)*factor_xy

          ! sols, soll, solsd, solld have unit of mks, so no conversion is needed
          rad_avgdata%sols_m(i)    = rad_avgdata%sols_m(i)    +cam_out%sols(i)    *factor_xy
          rad_avgdata%soll_m(i)    = rad_avgdata%soll_m(i)    +cam_out%soll(i)    *factor_xy
          rad_avgdata%solsd_m(i)   = rad_avgdata%solsd_m(i)   +cam_out%solsd(i)   *factor_xy
          rad_avgdata%solld_m(i)   = rad_avgdata%solld_m(i)   +cam_out%solld(i)   *factor_xy

          rad_avgdata%fsn200_m(i)  = rad_avgdata%fsn200_m(i)  +rd%fsn200(i)  *factor_xy
          rad_avgdata%fsn200c_m(i) = rad_avgdata%fsn200c_m(i) +rd%fsn200c(i) *factor_xy
          rad_avgdata%fsnr_m(i)    = rad_avgdata%fsnr_m(i)    +rd%fsnr(i)    *factor_xy
       end do
       rad_avgdata%qrs_m(:ncol,:pver)  = rad_avgdata%qrs_m(:ncol,:pver)  + qrs(:ncol,:pver) *factor_xy
       rad_avgdata%qrsc_m(:ncol,:pver) = rad_avgdata%qrsc_m(:ncol,:pver) + rd%qrsc(:ncol,:pver)*factor_xy
       do m=1,crm_nz
          k = pver-m+1
          rad_avgdata%qrs_crm(:ncol,ii,jj,m) = qrs(:ncol,k) / cpair
       end do
       do i=1, Nday
          do k=1, pver
             if((rd%liq_icld_vistau(IdxDay(i),k)+rd%ice_icld_vistau(IdxDay(i),k)).gt.1.0e-10_r8) then
                rad_avgdata%tot_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%tot_icld_vistau_m(IdxDay(i),k) +   &
                                      (rd%liq_icld_vistau(IdxDay(i),k)+rd%ice_icld_vistau(IdxDay(i),k)) * cld(i,k)
                rad_avgdata%nct_tot_icld_vistau_m(IdxDay(i),k) = rad_avgdata%nct_tot_icld_vistau_m(IdxDay(i),k) + cld(i,k)
             end if
             if(rd%liq_icld_vistau(IdxDay(i),k).gt.1.0e-10_r8) then
                rad_avgdata%liq_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%liq_icld_vistau_m(IdxDay(i),k) + &
                                                                 rd%liq_icld_vistau(IdxDay(i),k) * cld(i,k)
                rad_avgdata%nct_liq_icld_vistau_m(IdxDay(i),k) = rad_avgdata%nct_liq_icld_vistau_m(IdxDay(i),k) + cld(i,k)
             end if
             if(rd%ice_icld_vistau(IdxDay(i),k).gt.1.0e-10_r8) then
                rad_avgdata%ice_icld_vistau_m(IdxDay(i),k)  = rad_avgdata%ice_icld_vistau_m(IdxDay(i),k) + &
                                                                 rd%ice_icld_vistau(IdxDay(i),k) * cld(i,k)
                rad_avgdata%nct_ice_icld_vistau_m(IdxDay(i),k) = rad_avgdata%nct_ice_icld_vistau_m(IdxDay(i),k) + cld(i,k)
             end if
          end do
       end do
    end if ! dosw

    if (dolw) then
       call pbuf_get_field(pbuf, flns_idx, flns)
       call pbuf_get_field(pbuf, flnt_idx, flnt)
       do i=1,ncol
          rad_avgdata%flns_m(i)    = rad_avgdata%flns_m(i)    +flns(i)   *factor_xy
          rad_avgdata%flnt_m(i)    = rad_avgdata%flnt_m(i)    +flnt(i)   *factor_xy

          rad_avgdata%flut_m(i)    = rad_avgdata%flut_m(i)    +rd%flut(i)          *factor_xy
          rad_avgdata%flutc_m(i)   = rad_avgdata%flutc_m(i)   +rd%flutc(i)         *factor_xy
          rad_avgdata%flds_m(i)    = rad_avgdata%flds_m(i)    +cam_out%flwds(i) *factor_xy
          rad_avgdata%fldsc_m(i)   = rad_avgdata%fldsc_m(i)   +rd%fldsc(i)         *factor_xy
          rad_avgdata%flntc_m(i)   = rad_avgdata%flntc_m(i)   +rd%flntc(i)         *factor_xy
          rad_avgdata%fln200_m(i)  = rad_avgdata%fln200_m(i)  +rd%fln200(i)        *factor_xy
          rad_avgdata%fln200c_m(i) = rad_avgdata%fln200c_m(i) +rd%fln200c(i)       *factor_xy
          rad_avgdata%flnsc_m(i)   = rad_avgdata%flnsc_m(i)   +rd%flnsc(i)         *factor_xy
          rad_avgdata%flnr_m(i)    = rad_avgdata%flnr_m(i)    +rd%flnr(i)          *factor_xy
       end do
       rad_avgdata%qrl_m(:ncol,:pver)  = rad_avgdata%qrl_m(:ncol,:pver)  + qrl(:ncol,:pver)  *factor_xy
       rad_avgdata%qrlc_m(:ncol,:pver) = rad_avgdata%qrlc_m(:ncol,:pver) + rd%qrlc(:ncol,:pver) *factor_xy

       do m=1,crm_nz
          k = pver-m+1
          rad_avgdata%qrl_crm(:ncol,ii,jj,m) = qrl(:ncol,k) / cpair
       end do
    end if

    do m=1,crm_nz
       k = pver-m+1
       do i = 1,ncol
          ! for energy conservation
          rad_avgdata%crm_qrad(i,ii,jj,m) = (rad_avgdata%qrs_crm(i,ii,jj,m)+rad_avgdata%qrl_crm(i,ii,jj,m)) * state%pdel(i,k)
       end do
    end do

#endif
end subroutine spcam_radiation_col_finalize_sam1mom

end module spcam_drivers
