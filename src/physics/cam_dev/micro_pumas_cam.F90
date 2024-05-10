module micro_pumas_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use shr_kind_mod,   only: cl=>shr_kind_cl
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, pverp, psubcols
use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o, &
                          latvap, latice, mwh2o
use phys_control,   only: phys_getopts, use_hetfrz_classnuc


use physics_types,  only: physics_state, physics_ptend, &
                          physics_ptend_init, physics_state_copy, &
                          physics_update, physics_state_dealloc, &
                          physics_ptend_sum, physics_ptend_scale

use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dyn_time_lvls, &
                          pbuf_old_tim_idx, pbuf_get_index, dtype_r8, dtype_i4, &
                          pbuf_get_field, pbuf_set_field, col_type_subcol, &
                          pbuf_register_subcol
use constituents,   only: cnst_add, cnst_get_ind, &
                          cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst

use cldfrc2m,       only: rhmini=>rhmini_const

use cam_history,    only: addfld, add_default, outfld, horiz_only

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev

use micro_pumas_diags, only: proc_rates_type

use subcol_utils,   only: subcol_get_scheme

implicit none
private
save

public :: &
   micro_pumas_cam_readnl,          &
   micro_pumas_cam_register,        &
   micro_pumas_cam_init_cnst,       &
   micro_pumas_cam_implements_cnst, &
   micro_pumas_cam_init,            &
   micro_pumas_cam_tend,            &
   micro_mg_version,             &
   massless_droplet_destroyer

integer :: micro_mg_version     = 1      ! Version number for MG.
integer :: micro_mg_sub_version = 0      ! Second part of version number.

real(r8) :: micro_mg_dcs = -1._r8
real(r8), target, allocatable :: trop_levs(:)

logical :: microp_uniform       = .false.
logical :: micro_mg_adjust_cpt  = .false.

logical :: micro_do_massless_droplet_destroyer ! turn on/off destruction of massless droplets

character(len=16) :: micro_mg_precip_frac_method = 'max_overlap' ! type of precipitation fraction method

real(r8), parameter :: unset_r8 = huge(1.0_r8)

! Tunable namelist parameters (set in atm_in)
real(r8) :: micro_mg_berg_eff_factor   = unset_r8        ! berg efficiency factor
real(r8) :: micro_mg_accre_enhan_fact  = unset_r8        ! accretion enhancment factor
real(r8) :: micro_mg_autocon_fact      = unset_r8       ! autoconversion prefactor
real(r8) :: micro_mg_autocon_nd_exp    = unset_r8       ! autoconversion nd exponent
real(r8) :: micro_mg_autocon_lwp_exp   = unset_r8       ! autoconversion lwp exponent
real(r8) :: micro_mg_homog_size        = unset_r8     ! size of freezing homogeneous ice
real(r8) :: micro_mg_vtrmi_factor      = unset_r8        ! ice fall speed factor
real(r8) :: micro_mg_vtrms_factor      = unset_r8        ! snow fall speed factor
real(r8) :: micro_mg_effi_factor       = unset_r8        ! ice effective radius factor
real(r8) :: micro_mg_iaccr_factor      = unset_r8        ! ice accretion of cloud droplet
real(r8) :: micro_mg_max_nicons        = unset_r8  ! max allowed ice number concentration


logical, public :: do_cldliq ! Prognose cldliq flag
logical, public :: do_cldice ! Prognose cldice flag

integer :: num_steps ! Number of MG substeps

integer :: ncnst = 4       ! Number of constituents

! Namelist variables for option to specify constant cloud droplet/ice number
logical :: micro_mg_nccons = .false. ! set .true. to specify constant cloud droplet number
logical :: micro_mg_nicons = .false. ! set .true. to specify constant cloud ice number
logical :: micro_mg_ngcons = .false. ! set .true. to specify constant graupel/hail number
logical :: micro_mg_nrcons = .false. ! set .true. to specify constant rain number
logical :: micro_mg_nscons = .false. ! set .true. to specify constant snow number

! parameters for specified ice and droplet number concentration
! note: these are local in-cloud values, not grid-mean
real(r8) :: micro_mg_ncnst = 50.e6_r8 ! constant liquid droplet num concentration (m-3)
real(r8) :: micro_mg_ninst = 0.05e6_r8  ! ice num concentration when nicons=.true. (m-3)
real(r8) :: micro_mg_nrnst = 0.2e6_r8     ! rain  num concentration when nrcons=.true. (m-3)
real(r8) :: micro_mg_nsnst = 0.005e6_r8 ! snow num concentration when nscons=.true. (m-3)
real(r8) :: micro_mg_ngnst = 0.0005e6_r8 ! graupel/hail num concentration when ngcons=.true. (m-3)

logical, public ::   micro_mg_do_graupel
logical, public ::   micro_mg_do_hail

! switches for IFS like behavior
logical  ::  micro_mg_evap_sed_off = .false.      ! Turn off evaporation/sublimation based on cloud fraction for sedimenting condensate
logical  ::  micro_mg_icenuc_rh_off  = .false.    ! Remove RH conditional from ice nucleation
logical  ::  micro_mg_icenuc_use_meyers = .false. ! Meyers Ice Nucleation
logical  ::  micro_mg_evap_scl_ifs = .false.      ! Scale evaporation as IFS does
logical  ::  micro_mg_evap_rhthrsh_ifs = .false.  ! Evap RH threhold following IFS
logical  ::  micro_mg_rainfreeze_ifs = .false.    ! Rain freezing at 0C following IFS
logical  ::  micro_mg_ifs_sed = .false.           ! Snow sedimentation = 1 m/s following IFS
logical  ::  micro_mg_precip_fall_corr = .false.    ! Precip fall speed following IFS (does not go to zero)

logical  ::  micro_mg_implicit_fall = .false. !Implicit fall speed (sedimentation) for hydrometeors

logical  ::  micro_mg_accre_sees_auto = .false.    !Accretion sees autoconverted rain

character(len=10), parameter :: &      ! Constituent names
   cnst_names(10) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                     'RAINQM', 'SNOWQM','NUMRAI','NUMSNO','GRAUQM','NUMGRA'/)

integer :: &
   ixq = -1,           &! water vapor
   ixcldliq = -1,      &! cloud liquid amount index
   ixcldice = -1,      &! cloud ice amount index
   ixnumliq = -1,      &! cloud liquid number index
   ixnumice = -1,      &! cloud ice water index
   ixrain = -1,        &! rain index
   ixsnow = -1,        &! snow index
   ixnumrain = -1,     &! rain number index
   ixnumsnow = -1,     &! snow number index
   ixgraupel = -1,     &! graupel index
   ixnumgraupel = -1   ! graupel number index

! Physics buffer indices for fields registered by this module
integer :: &
   cldo_idx,           &
   qme_idx,            &
   prain_idx,          &
   nevapr_idx,         &
   wsedl_idx,          &
   rei_idx,            &
   sadice_idx,         &
   sadsnow_idx,        &
   rel_idx,            &
   dei_idx,            &
   mu_idx,             &
   prer_evap_idx,            &
   lambdac_idx,        &
   iciwpst_idx,        &
   iclwpst_idx,        &
   des_idx,            &
   icswp_idx,          &
   cldfsnow_idx,       &
   degrau_idx = -1,    &
   icgrauwp_idx = -1,  &
   cldfgrau_idx = -1,  &
   rate1_cw2pr_st_idx = -1, &
   ls_flxprc_idx,      &
   ls_flxsnw_idx,      &
   relvar_idx,         &
   cmeliq_idx,         &
   accre_enhan_idx

! Fields for UNICON
integer :: &
     am_evp_st_idx,      &! Evaporation area of stratiform precipitation
     evprain_st_idx,     &! Evaporation rate of stratiform rain [kg/kg/s]. >= 0.
     evpsnow_st_idx       ! Evaporation rate of stratiform snow [kg/kg/s]. >= 0.

! Fields needed as inputs to COSP
integer :: &
     ls_mrprc_idx,    ls_mrsnw_idx,    &
     ls_reffrain_idx, ls_reffsnow_idx, &
     cv_reffliq_idx,  cv_reffice_idx

! Fields needed by Park macrophysics
integer :: &
     cc_t_idx,  cc_qv_idx, &
     cc_ql_idx, cc_qi_idx, &
     cc_nl_idx, cc_ni_idx, &
     cc_qlst_idx

! Used to replace aspects of MG microphysics
! (e.g. by CARMA)
integer :: &
     tnd_qsnow_idx = -1, &
     tnd_nsnow_idx = -1, &
     re_ice_idx = -1

! Index fields for precipitation efficiency.
integer :: &
     acpr_idx = -1, &
     acgcme_idx = -1, &
     acnum_idx = -1

! Physics buffer indices for fields registered by other modules
integer :: &
   ast_idx = -1,            &
   cld_idx = -1,            &
   concld_idx = -1,         &
   qsatfac_idx = -1

! Pbuf fields needed for subcol_SILHS
integer :: &
     qrain_idx=-1, qsnow_idx=-1,    &
     nrain_idx=-1, nsnow_idx=-1,    &
     qcsedten_idx=-1, qrsedten_idx=-1, &
     qisedten_idx=-1, qssedten_idx=-1, &
     vtrmc_idx=-1, umr_idx=-1, &
     vtrmi_idx=-1, ums_idx=-1, &
     qcsevap_idx=-1, qisevap_idx=-1

integer :: &
   naai_idx = -1,           &
   naai_hom_idx = -1,       &
   npccn_idx = -1,          &
   rndst_idx = -1,          &
   nacon_idx = -1,          &
   prec_str_idx = -1,       &
   snow_str_idx = -1,       &
   prec_pcw_idx = -1,       &
   snow_pcw_idx = -1,       &
   prec_sed_idx = -1,       &
   snow_sed_idx = -1

! pbuf fields for heterogeneous freezing
integer :: &
   frzimm_idx = -1, &
   frzcnt_idx = -1, &
   frzdep_idx = -1

logical :: allow_sed_supersat                      ! allow supersaturated conditions after sedimentation loop
character(len=16) :: micro_mg_warm_rain= 'kk2000'  ! 'tau', 'emulated', 'sb2001' and ' kk2000'

integer :: bergso_idx = -1

!===============================================================================
contains
!===============================================================================

subroutine micro_pumas_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8, &
                             mpi_logical, mpi_character

  use stochastic_emulated_cam, only: stochastic_emulated_readnl
  use stochastic_tau_cam,      only: stochastic_tau_readnl

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  logical :: micro_mg_do_cldice = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
  logical :: micro_mg_do_cldliq = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq
  integer :: micro_mg_num_steps = 1      ! Number of substepping iterations done by MG (1.5 only for now).


  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'micro_pumas_cam_readnl'

  namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
       micro_mg_do_cldice, micro_mg_do_cldliq, micro_mg_num_steps, &
       microp_uniform, micro_mg_dcs, micro_mg_precip_frac_method, &
       micro_mg_berg_eff_factor, micro_mg_warm_rain, micro_mg_adjust_cpt, &
       micro_mg_do_hail, micro_mg_do_graupel, micro_mg_ngcons, micro_mg_ngnst, &
       micro_mg_vtrmi_factor, micro_mg_vtrms_factor, micro_mg_effi_factor, &
       micro_mg_iaccr_factor, micro_mg_max_nicons, micro_mg_accre_enhan_fact, &
       micro_mg_autocon_fact, micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, micro_mg_homog_size, &
       micro_mg_nccons, micro_mg_nicons, micro_mg_ncnst, micro_mg_ninst, &
       micro_mg_nrcons, micro_mg_nscons, micro_mg_nrnst, micro_mg_nsnst, &
       micro_do_massless_droplet_destroyer, &
       micro_mg_evap_sed_off, micro_mg_icenuc_rh_off, micro_mg_icenuc_use_meyers, &
       micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs, &
       micro_mg_rainfreeze_ifs, micro_mg_ifs_sed, micro_mg_precip_fall_corr, &
       micro_mg_accre_sees_auto, micro_mg_implicit_fall

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_mg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_mg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(sub // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! set local variables
     do_cldice = micro_mg_do_cldice
     do_cldliq = micro_mg_do_cldliq
     num_steps = micro_mg_num_steps

     ! Verify that version numbers are valid.
     select case (micro_mg_version)
     case (2)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 2.0
        case default
           call bad_version_endrun()
        end select
     case (3)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 3.0
        case default
           call bad_version_endrun()
        end select
     case default
        call bad_version_endrun()
     end select

     if (micro_mg_dcs < 0._r8) call endrun( "micro_pumas_cam_readnl: &
              &micro_mg_dcs has not been set to a valid value.")

     if (micro_mg_version < 3) then

        if(micro_mg_do_graupel .or. micro_mg_do_hail ) then
           call endrun ("micro_pumas_cam_readnl: Micro_mg_do_graupel and micro_mg_do_hail &
                &must be false for MG versions before MG3.")
        end if

     else ! micro_mg_version = 3 or greater

        if(micro_mg_do_graupel .and. micro_mg_do_hail ) then
           call endrun ("micro_pumas_cam_readnl: Only one of micro_mg_do_graupel or &
                &micro_mg_do_hail may be true at a time.")
        end if

     end if

  end if

  ! Broadcast namelist variables
  call mpi_bcast(micro_mg_version, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_version")

  call mpi_bcast(micro_mg_sub_version, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_sub_version")

  call mpi_bcast(do_cldice, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_cldice")

  call mpi_bcast(do_cldliq, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: do_cldliq")

  call mpi_bcast(num_steps, 1, mpi_integer, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: num_steps")

  call mpi_bcast(microp_uniform, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: microp_uniform")

  call mpi_bcast(micro_mg_dcs, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_dcs")

  call mpi_bcast(micro_mg_berg_eff_factor, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_berg_eff_factor")

  call mpi_bcast(micro_mg_accre_enhan_fact, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_accre_enhan_fact")

  call mpi_bcast(micro_mg_autocon_fact, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_autocon_fact")

  call mpi_bcast(micro_mg_autocon_nd_exp, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_autocon_nd_exp")

  call mpi_bcast(micro_mg_autocon_lwp_exp, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_autocon_lwp_exp")

  call mpi_bcast(micro_mg_homog_size, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_homog_size")

  call mpi_bcast(micro_mg_vtrmi_factor, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_vtrmi_factor")

  call mpi_bcast(micro_mg_vtrms_factor, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_vtrms_factor")

  call mpi_bcast(micro_mg_effi_factor, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_effi_factor")

  call mpi_bcast(micro_mg_iaccr_factor, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_iaccr_factor")

  call mpi_bcast(micro_mg_max_nicons, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_max_nicons")

  call mpi_bcast(micro_mg_precip_frac_method, 16, mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_precip_frac_method")

  call mpi_bcast(micro_mg_warm_rain, 16, mpi_character, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_warm_rain")

  call mpi_bcast(micro_mg_adjust_cpt, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_adjust_cpt")

  call mpi_bcast(micro_mg_nccons, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nccons")

  call mpi_bcast(micro_mg_nicons, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nicons")

  call mpi_bcast(micro_mg_nrcons, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nrcons")

  call mpi_bcast(micro_mg_nscons, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nscons")

  call mpi_bcast(micro_mg_ncnst, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_ncnst")

  call mpi_bcast(micro_mg_ninst, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_ninst")

  call mpi_bcast(micro_mg_nrnst, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nrnst")

  call mpi_bcast(micro_mg_nsnst, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_nsnst")

  call mpi_bcast(micro_mg_do_hail, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_do_hail")

  call mpi_bcast(micro_mg_do_graupel, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_do_graupel")

  call mpi_bcast(micro_mg_ngcons, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_ngcons")

  call mpi_bcast(micro_mg_ngnst, 1, mpi_real8, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_ngnst")

  call mpi_bcast(micro_do_massless_droplet_destroyer, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_do_massless_droplet_destroyer")

  call mpi_bcast(micro_mg_evap_sed_off, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_evap_sed_off")

  call mpi_bcast(micro_mg_icenuc_rh_off, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_icenuc_rh_off")

  call mpi_bcast(micro_mg_icenuc_use_meyers, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_icenuc_use_meyers")

  call mpi_bcast(micro_mg_evap_scl_ifs, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_evap_scl_ifs")

  call mpi_bcast(micro_mg_evap_rhthrsh_ifs, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_evap_rhthrsh_ifs")

  call mpi_bcast(micro_mg_rainfreeze_ifs, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_rainfreeze_ifs")

  call mpi_bcast(micro_mg_ifs_sed, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_ifs_sed")

  call mpi_bcast(micro_mg_precip_fall_corr, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_precip_fall_corr")

  call mpi_bcast(micro_mg_implicit_fall, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_implicit_fall")

  call mpi_bcast(micro_mg_accre_sees_auto, 1, mpi_logical, mstrid, mpicom, ierr)
  if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: micro_mg_accre_sees_auto")

  if(micro_mg_berg_eff_factor == unset_r8) call endrun(sub//": FATAL: micro_mg_berg_eff_factor is not set")
  if(micro_mg_accre_enhan_fact == unset_r8) call endrun(sub//": FATAL: micro_mg_accre_enhan_fact is not set")
  if(micro_mg_autocon_fact == unset_r8) call endrun(sub//": FATAL: micro_mg_autocon_fact is not set")
  if(micro_mg_autocon_nd_exp == unset_r8) call endrun(sub//": FATAL: micro_mg_autocon_nd_exp is not set")
  if(micro_mg_autocon_lwp_exp == unset_r8) call endrun(sub//": FATAL: micro_mg_autocon_lwp_exp is not set")
  if(micro_mg_homog_size == unset_r8) call endrun(sub//": FATAL: micro_mg_homog_size is not set")
  if(micro_mg_vtrmi_factor == unset_r8) call endrun(sub//": FATAL: micro_mg_vtrmi_factor is not set")
  if(micro_mg_vtrms_factor == unset_r8) call endrun(sub//": FATAL: micro_mg_vtrms_factor is not set")
  if(micro_mg_effi_factor == unset_r8) call endrun(sub//": FATAL: micro_mg_effi_factor is not set")
  if(micro_mg_iaccr_factor == unset_r8) call endrun(sub//": FATAL: micro_mg_iaccr_factor is not set")
  if(micro_mg_max_nicons == unset_r8) call endrun(sub//": FATAL: micro_mg_max_nicons is not set")

  if (masterproc) then

     write(iulog,*) 'MG microphysics namelist:'
     write(iulog,*) '  micro_mg_version            = ', micro_mg_version
     write(iulog,*) '  micro_mg_sub_version        = ', micro_mg_sub_version
     write(iulog,*) '  micro_mg_do_cldice          = ', do_cldice
     write(iulog,*) '  micro_mg_do_cldliq          = ', do_cldliq
     write(iulog,*) '  micro_mg_num_steps          = ', num_steps
     write(iulog,*) '  microp_uniform              = ', microp_uniform
     write(iulog,*) '  micro_mg_dcs                = ', micro_mg_dcs
     write(iulog,*) '  micro_mg_berg_eff_factor    = ', micro_mg_berg_eff_factor
     write(iulog,*) '  micro_mg_accre_enhan_fact   = ', micro_mg_accre_enhan_fact
     write(iulog,*) '  micro_mg_autocon_fact       = ' , micro_mg_autocon_fact
     write(iulog,*) '  micro_mg_autocon_nd_exp     = ' , micro_mg_autocon_nd_exp
     write(iulog,*) '  micro_mg_autocon_lwp_exp    = ' , micro_mg_autocon_lwp_exp
     write(iulog,*) '  micro_mg_homog_size         = ', micro_mg_homog_size
     write(iulog,*) '  micro_mg_vtrmi_factor       = ', micro_mg_vtrmi_factor
     write(iulog,*) '  micro_mg_vtrms_factor       = ', micro_mg_vtrms_factor
     write(iulog,*) '  micro_mg_effi_factor        = ', micro_mg_effi_factor
     write(iulog,*) '  micro_mg_iaccr_factor       = ', micro_mg_iaccr_factor
     write(iulog,*) '  micro_mg_max_nicons         = ', micro_mg_max_nicons
     write(iulog,*) '  micro_mg_precip_frac_method = ', micro_mg_precip_frac_method
     write(iulog,*) '  micro_mg_warm_rain          = ', micro_mg_warm_rain
     write(iulog,*) '  micro_mg_adjust_cpt         = ', micro_mg_adjust_cpt
     write(iulog,*) '  micro_mg_nccons             = ', micro_mg_nccons
     write(iulog,*) '  micro_mg_nicons             = ', micro_mg_nicons
     write(iulog,*) '  micro_mg_ncnst              = ', micro_mg_ncnst
     write(iulog,*) '  micro_mg_ninst              = ', micro_mg_ninst
     write(iulog,*) '  micro_mg_ngcons             = ', micro_mg_ngcons
     write(iulog,*) '  micro_mg_ngnst              = ', micro_mg_ngnst
     write(iulog,*) '  micro_mg_do_hail            = ', micro_mg_do_hail
     write(iulog,*) '  micro_mg_do_graupel         = ', micro_mg_do_graupel
     write(iulog,*) '  micro_do_massless_droplet_destroyer = ', micro_do_massless_droplet_destroyer
     write(iulog,*) '  micro_mg_nrcons             = ', micro_mg_nrcons
     write(iulog,*) '  micro_mg_nscons             = ', micro_mg_nscons
     write(iulog,*) '  micro_mg_nrnst              = ', micro_mg_nrnst
     write(iulog,*) '  micro_mg_nsnst              = ', micro_mg_nsnst
     write(iulog,*) '  micro_mg_evap_sed_off       = ', micro_mg_evap_sed_off
     write(iulog,*) '  micro_mg_icenuc_rh_off      = ', micro_mg_icenuc_rh_off
     write(iulog,*) '  micro_mg_icenuc_use_meyers  = ', micro_mg_icenuc_use_meyers
     write(iulog,*) '  micro_mg_evap_scl_ifs       = ', micro_mg_evap_scl_ifs
     write(iulog,*) '  micro_mg_evap_rhthrsh_ifs   = ', micro_mg_evap_rhthrsh_ifs
     write(iulog,*) '  micro_mg_rainfreeze_ifs     = ', micro_mg_rainfreeze_ifs
     write(iulog,*) '  micro_mg_ifs_sed            = ', micro_mg_ifs_sed
     write(iulog,*) '  micro_mg_precip_fall_corr     = ', micro_mg_precip_fall_corr
     write(iulog,*) '  micro_mg_implicit_fall     = ', micro_mg_implicit_fall
     write(iulog,*) '  micro_mg_accre_sees_auto     = ', micro_mg_accre_sees_auto
  end if

  ! Read in the emulated or tau namelist if needed
  if( trim(micro_mg_warm_rain) == 'emulated') then
     call stochastic_emulated_readnl(nlfile)
  else if (trim(micro_mg_warm_rain) == 'tau') then
     call stochastic_tau_readnl(nlfile)
  end if

contains

  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for MG microphysics: ", &
         micro_mg_version,".",micro_mg_sub_version
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_pumas_cam_readnl

!================================================================================================

subroutine micro_pumas_cam_register
   use cam_history_support, only: add_vert_coord, hist_dimension_values
   use cam_abortutils,      only: handle_allocate_error

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------

   logical :: prog_modal_aero
   logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics
   logical :: found

   integer :: i, ierr
   real(r8) :: all_levs(pver)

   allocate(trop_levs(pver-top_lev+1), stat=ierr)
   call handle_allocate_error(ierr, 'micro_pumas_cam_register', 'trop_levs')

   call phys_getopts(use_subcol_microp_out    = use_subcol_microp, &
                     prog_modal_aero_out      = prog_modal_aero)

   ! Register microphysics constituents and save indices.

   call cnst_add(cnst_names(1), mwh2o, cpair, 0._r8, ixcldliq, &
      longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
   call cnst_add(cnst_names(2), mwh2o, cpair, 0._r8, ixcldice, &
      longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

   call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
      longname='Grid box averaged cloud liquid number', is_convtran1=.true.)
   call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
      longname='Grid box averaged cloud ice number', is_convtran1=.true.)

   ! Add history coordinate for DDT nlev
   call hist_dimension_values('lev',all_levs, 1, pver, found)

   if (found) then
      trop_levs(1:pver-top_lev+1) = all_levs(top_lev:pver)
      call add_vert_coord('trop_cld_lev', pver-top_lev+1,                          &
            'troposphere hybrid level at midpoints (1000*(A+B))', 'hPa', trop_levs,  &
            positive='down' )
   else
      call endrun( "micro_pumas_cam_register: unable to find dimension field 'lev'")
   end if


! ---- Note is_convtran1 is set to .true.
   call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain, &
        longname='Grid box averaged rain amount', is_convtran1=.true.)
   call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixsnow, &
        longname='Grid box averaged snow amount', is_convtran1=.true.)
   call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain, &
        longname='Grid box averaged rain number', is_convtran1=.true.)
   call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixnumsnow, &
        longname='Grid box averaged snow number', is_convtran1=.true.)

   if (micro_mg_version > 2) then
         call cnst_add(cnst_names(9), mwh2o, cpair, 0._r8, ixgraupel, &
              longname='Grid box averaged graupel/hail amount', is_convtran1=.true.)
         call cnst_add(cnst_names(10), mwh2o, cpair, 0._r8, ixnumgraupel, &
              longname='Grid box averaged graupel/hail number', is_convtran1=.true.)
   end if

   ! Request physics buffer space for fields that persist across timesteps.

   call pbuf_add_field('CLDO','global',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)

   ! Physics buffer variables for convective cloud properties.

   call pbuf_add_field('QME',        'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
   call pbuf_add_field('PRAIN',      'physpkg',dtype_r8,(/pcols,pver/), prain_idx)
   call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)
   call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)
   call pbuf_add_field('BERGSO',     'physpkg',dtype_r8,(/pcols,pver/), bergso_idx)

   call pbuf_add_field('WSEDL',      'physpkg',dtype_r8,(/pcols,pver/), wsedl_idx)

   call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
   call pbuf_add_field('SADICE',     'physpkg',dtype_r8,(/pcols,pver/), sadice_idx)
   call pbuf_add_field('SADSNOW',    'physpkg',dtype_r8,(/pcols,pver/), sadsnow_idx)
   call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)

   ! Mitchell ice effective diameter for radiation
   call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
   ! Size distribution shape parameter for radiation
   call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
   ! Size distribution shape parameter for radiation
   call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

   ! Stratiform only in cloud ice water path for radiation
   call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)
   ! Stratiform in cloud liquid water path for radiation
   call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)

   ! Snow effective diameter for radiation
   call pbuf_add_field('DES',        'physpkg',dtype_r8,(/pcols,pver/), des_idx)
   ! In cloud snow water path for radiation
   call pbuf_add_field('ICSWP',      'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)
   ! Cloud fraction for liquid drops + snow
   call pbuf_add_field('CLDFSNOW ',  'physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx)

   if (micro_mg_version > 2) then
      ! Graupel effective diameter for radiation
      call pbuf_add_field('DEGRAU',        'physpkg',dtype_r8,(/pcols,pver/), degrau_idx)
      ! In cloud snow water path for radiation
      call pbuf_add_field('ICGRAUWP',      'physpkg',dtype_r8,(/pcols,pver/), icgrauwp_idx)
      ! Cloud fraction for liquid drops + graupel
      call pbuf_add_field('CLDFGRAU',      'physpkg',dtype_r8,(/pcols,pver/), cldfgrau_idx)
   end if

   if (prog_modal_aero) then
      call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)
   endif

   call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
   call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)


   ! Fields needed as inputs to COSP
   call pbuf_add_field('LS_MRPRC',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrprc_idx)
   call pbuf_add_field('LS_MRSNW',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrsnw_idx)
   call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
   call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
   call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
   call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

   ! CC_* Fields needed by Park macrophysics
   call pbuf_add_field('CC_T',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_t_idx)
   call pbuf_add_field('CC_qv',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qv_idx)
   call pbuf_add_field('CC_ql',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ql_idx)
   call pbuf_add_field('CC_qi',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qi_idx)
   call pbuf_add_field('CC_nl',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_nl_idx)
   call pbuf_add_field('CC_ni',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ni_idx)
   call pbuf_add_field('CC_qlst',  'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

   ! Fields for UNICON
   call pbuf_add_field('am_evp_st',  'global', dtype_r8, (/pcols,pver/), am_evp_st_idx)
   call pbuf_add_field('evprain_st', 'global', dtype_r8, (/pcols,pver/), evprain_st_idx)
   call pbuf_add_field('evpsnow_st', 'global', dtype_r8, (/pcols,pver/), evpsnow_st_idx)

   ! Register subcolumn pbuf fields
   if (use_subcol_microp) then
      ! Global pbuf fields
      call pbuf_register_subcol('CLDO',        'micro_pumas_cam_register', cldo_idx)

      ! CC_* Fields needed by Park macrophysics
      call pbuf_register_subcol('CC_T',        'micro_pumas_cam_register', cc_t_idx)
      call pbuf_register_subcol('CC_qv',       'micro_pumas_cam_register', cc_qv_idx)
      call pbuf_register_subcol('CC_ql',       'micro_pumas_cam_register', cc_ql_idx)
      call pbuf_register_subcol('CC_qi',       'micro_pumas_cam_register', cc_qi_idx)
      call pbuf_register_subcol('CC_nl',       'micro_pumas_cam_register', cc_nl_idx)
      call pbuf_register_subcol('CC_ni',       'micro_pumas_cam_register', cc_ni_idx)
      call pbuf_register_subcol('CC_qlst',     'micro_pumas_cam_register', cc_qlst_idx)

      ! Physpkg pbuf fields
      ! Physics buffer variables for convective cloud properties.

      call pbuf_register_subcol('QME',         'micro_pumas_cam_register', qme_idx)
      call pbuf_register_subcol('PRAIN',       'micro_pumas_cam_register', prain_idx)
      call pbuf_register_subcol('NEVAPR',      'micro_pumas_cam_register', nevapr_idx)
      call pbuf_register_subcol('PRER_EVAP',   'micro_pumas_cam_register', prer_evap_idx)
      call pbuf_register_subcol('BERGSO',      'micro_pumas_cam_register', bergso_idx)

      call pbuf_register_subcol('WSEDL',       'micro_pumas_cam_register', wsedl_idx)

      call pbuf_register_subcol('REI',         'micro_pumas_cam_register', rei_idx)
      call pbuf_register_subcol('SADICE',      'micro_pumas_cam_register', sadice_idx)
      call pbuf_register_subcol('SADSNOW',     'micro_pumas_cam_register', sadsnow_idx)
      call pbuf_register_subcol('REL',         'micro_pumas_cam_register', rel_idx)

      ! Mitchell ice effective diameter for radiation
      call pbuf_register_subcol('DEI',         'micro_pumas_cam_register', dei_idx)
      ! Size distribution shape parameter for radiation
      call pbuf_register_subcol('MU',          'micro_pumas_cam_register', mu_idx)
      ! Size distribution shape parameter for radiation
      call pbuf_register_subcol('LAMBDAC',     'micro_pumas_cam_register', lambdac_idx)

      ! Stratiform only in cloud ice water path for radiation
      call pbuf_register_subcol('ICIWPST',     'micro_pumas_cam_register', iciwpst_idx)
      ! Stratiform in cloud liquid water path for radiation
      call pbuf_register_subcol('ICLWPST',     'micro_pumas_cam_register', iclwpst_idx)

      ! Snow effective diameter for radiation
      call pbuf_register_subcol('DES',         'micro_pumas_cam_register', des_idx)
      ! In cloud snow water path for radiation
      call pbuf_register_subcol('ICSWP',       'micro_pumas_cam_register', icswp_idx)
      ! Cloud fraction for liquid drops + snow
      call pbuf_register_subcol('CLDFSNOW ',   'micro_pumas_cam_register', cldfsnow_idx)

      if (micro_mg_version > 2) then
         ! Graupel effective diameter for radiation
         call pbuf_register_subcol('DEGRAU',         'micro_pumas_cam_register', degrau_idx)
         ! In cloud snow water path for radiation
         call pbuf_register_subcol('ICGRAUWP',       'micro_pumas_cam_register', icgrauwp_idx)
         ! Cloud fraction for liquid drops + snow
         call pbuf_register_subcol('CLDFGRAU',   'micro_pumas_cam_register', cldfgrau_idx)
      end if

      if (prog_modal_aero) then
         call pbuf_register_subcol('RATE1_CW2PR_ST', 'micro_pumas_cam_register', rate1_cw2pr_st_idx)
      end if

      call pbuf_register_subcol('LS_FLXPRC',   'micro_pumas_cam_register', ls_flxprc_idx)
      call pbuf_register_subcol('LS_FLXSNW',   'micro_pumas_cam_register', ls_flxsnw_idx)

      ! Fields needed as inputs to COSP
      call pbuf_register_subcol('LS_MRPRC',    'micro_pumas_cam_register', ls_mrprc_idx)
      call pbuf_register_subcol('LS_MRSNW',    'micro_pumas_cam_register', ls_mrsnw_idx)
      call pbuf_register_subcol('LS_REFFRAIN', 'micro_pumas_cam_register', ls_reffrain_idx)
      call pbuf_register_subcol('LS_REFFSNOW', 'micro_pumas_cam_register', ls_reffsnow_idx)
      call pbuf_register_subcol('CV_REFFLIQ',  'micro_pumas_cam_register', cv_reffliq_idx)
      call pbuf_register_subcol('CV_REFFICE',  'micro_pumas_cam_register', cv_reffice_idx)
   end if

   ! Additional pbuf for CARMA interface
   if (.not. do_cldice) then
      call pbuf_add_field('TND_QSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_qsnow_idx)
      call pbuf_add_field('TND_NSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_nsnow_idx)
      call pbuf_add_field('RE_ICE',     'physpkg',dtype_r8,(/pcols,pver/), re_ice_idx)
   end if

   ! Precipitation efficiency fields across timesteps.
   call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
   call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
   call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

   ! SGS variability  -- These could be reset by CLUBB so they need to be grid only
   call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
   call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

   ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
   if (subcol_get_scheme() == 'SILHS') then
      call pbuf_add_field('QRAIN',   'global',dtype_r8,(/pcols,pver/), qrain_idx)
      call pbuf_add_field('QSNOW',   'global',dtype_r8,(/pcols,pver/), qsnow_idx)
      call pbuf_add_field('NRAIN',   'global',dtype_r8,(/pcols,pver/), nrain_idx)
      call pbuf_add_field('NSNOW',   'global',dtype_r8,(/pcols,pver/), nsnow_idx)

      ! Fields for subcol_SILHS hole filling
      ! Note -- hole filling is on the grid, so pbuf_register_setcols do not need to be called for these pbuf fields
      call pbuf_add_field('QCSEDTEN', 'physpkg', dtype_r8, (/pcols,pver/), qcsedten_idx)
      call pbuf_add_field('QRSEDTEN', 'physpkg', dtype_r8, (/pcols,pver/), qrsedten_idx)
      call pbuf_add_field('QISEDTEN', 'physpkg', dtype_r8, (/pcols,pver/), qisedten_idx)
      call pbuf_add_field('QSSEDTEN', 'physpkg', dtype_r8, (/pcols,pver/), qssedten_idx)
      call pbuf_add_field('VTRMC', 'physpkg', dtype_r8, (/pcols,pver/), vtrmc_idx)
      call pbuf_add_field('UMR', 'physpkg', dtype_r8, (/pcols,pver/), umr_idx)
      call pbuf_add_field('VTRMI', 'physpkg', dtype_r8, (/pcols,pver/), vtrmi_idx)
      call pbuf_add_field('UMS', 'physpkg', dtype_r8, (/pcols,pver/), ums_idx)
      call pbuf_add_field('QCSEVAP', 'physpkg', dtype_r8, (/pcols,pver/), qcsevap_idx)
      call pbuf_add_field('QISEVAP', 'physpkg', dtype_r8, (/pcols,pver/), qisevap_idx)
   end if

end subroutine micro_pumas_cam_register

!===============================================================================

function micro_pumas_cam_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: micro_pumas_cam_implements_cnst    ! return value

   !-----------------------------------------------------------------------

   micro_pumas_cam_implements_cnst = any(name == cnst_names)

end function micro_pumas_cam_implements_cnst

!===============================================================================

subroutine micro_pumas_cam_init_cnst(name, latvals, lonvals, mask, q)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name       ! constituent name
   real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev
   !-----------------------------------------------------------------------
   integer :: k

   if (micro_pumas_cam_implements_cnst(name)) then
     do k = 1, size(q, 2)
       where(mask)
         q(:, k) = 0.0_r8
       end where
     end do
   end if

end subroutine micro_pumas_cam_init_cnst

!===============================================================================

subroutine micro_pumas_cam_init(pbuf2d)
   use time_manager,   only: is_first_step
   use micro_pumas_utils, only: micro_pumas_utils_init
   use micro_pumas_v1, only: micro_mg_init3_0 => micro_pumas_init
   use stochastic_tau_cam, only:  stochastic_tau_init_cam
   use stochastic_emulated_cam, only:  stochastic_emulated_init_cam

   !-----------------------------------------------------------------------
   !
   ! Initialization for MG microphysics
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   integer :: m, mm
   logical :: history_amwg         ! output the variables used by the AMWG diag package
   logical :: history_budget       ! Output tendencies and state variables for CAM4
                                   ! temperature, water vapor, cloud ice and cloud
                                   ! liquid budgets.
   logical :: use_subcol_microp
   logical :: do_clubb_sgs
   integer :: budget_histfile      ! output history file number for budget fields
   integer :: ierr
   character(128) :: errstring     ! return status (non-blank for error return)

   character(len=cl) :: stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
                                       stochastic_emulated_filename_output_scale

   !-----------------------------------------------------------------------

   call phys_getopts(use_subcol_microp_out=use_subcol_microp, &
                     do_clubb_sgs_out     =do_clubb_sgs)

   if (do_clubb_sgs) then
     allow_sed_supersat = .false.
   else
     allow_sed_supersat = .true.
   endif

   if (masterproc) then
      write(iulog,"(A,I2,A,I2)") "Initializing MG version ",micro_mg_version,".",micro_mg_sub_version
      if (.not. do_cldliq) &
           write(iulog,*) "MG prognostic cloud liquid has been turned off via namelist."
      if (.not. do_cldice) &
           write(iulog,*) "MG prognostic cloud ice has been turned off via namelist."
      write(iulog,*) "Number of microphysics substeps is: ",num_steps
   end if

   ! Set constituent number for later loops.
   if(micro_mg_version == 2) then
         ncnst = 8
   else
         ncnst = 10
   end if

   ! If Machine learning is turned on, perform its initializations
   if (trim(micro_mg_warm_rain) == 'tau') then
      call stochastic_tau_init_cam()
   else if( trim(micro_mg_warm_rain) == 'emulated') then
      call stochastic_emulated_init_cam(stochastic_emulated_filename_quantile, &
                                       stochastic_emulated_filename_input_scale, &
                                       stochastic_emulated_filename_output_scale)
   end if

   call micro_mg_init3_0( &
           r8, gravit, rair, rh2o, cpair, &
           tmelt, latvap, latice, rhmini, &
           micro_mg_dcs,                  &
           micro_mg_do_hail,micro_mg_do_graupel, &
           microp_uniform, do_cldice, use_hetfrz_classnuc, &
           micro_mg_precip_frac_method, micro_mg_berg_eff_factor, &
           micro_mg_accre_enhan_fact , &
           micro_mg_autocon_fact , micro_mg_autocon_nd_exp, micro_mg_autocon_lwp_exp, micro_mg_homog_size, &
           micro_mg_vtrmi_factor, micro_mg_vtrms_factor, micro_mg_effi_factor, &
           micro_mg_iaccr_factor, micro_mg_max_nicons, &
           allow_sed_supersat, micro_mg_warm_rain, &
           micro_mg_evap_sed_off, micro_mg_icenuc_rh_off, micro_mg_icenuc_use_meyers, &
           micro_mg_evap_scl_ifs, micro_mg_evap_rhthrsh_ifs, &
           micro_mg_rainfreeze_ifs,  micro_mg_ifs_sed, micro_mg_precip_fall_corr,&
           micro_mg_accre_sees_auto, micro_mg_implicit_fall, &
           micro_mg_nccons, micro_mg_nicons, micro_mg_ncnst, &
           micro_mg_ninst, micro_mg_ngcons, micro_mg_ngnst, &
           micro_mg_nrcons,  micro_mg_nrnst, micro_mg_nscons, micro_mg_nsnst, &
           stochastic_emulated_filename_quantile, stochastic_emulated_filename_input_scale, &
           stochastic_emulated_filename_output_scale, iulog, errstring)

   call handle_errmsg(errstring, subname="micro_pumas_cam_init")

   ! Retrieve the index for water vapor
   call cnst_get_ind('Q', ixq)

   ! Register history variables
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow, ixgraupel /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm)                   )
         call addfld(sflxnam(mm),    horiz_only, 'A',   'kg/m2/s', trim(cnst_name(mm))//' surface flux')
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow, ixnumgraupel /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', cnst_longname(mm)                   )
         call addfld(sflxnam(mm),    horiz_only, 'A',   '1/m2/s', trim(cnst_name(mm))//' surface flux')
      else
         call endrun( "micro_pumas_cam_init: &
              &Could not call addfld for constituent with unknown units.")
      endif
   end do

   call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )

   call addfld(apcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
   call addfld(apcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics'  )
   call addfld(bpcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
   call addfld(bpcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics' )

   if (micro_mg_version > 2) then
      call addfld(apcnst(ixgraupel), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixgraupel))//' after physics'  )
      call addfld(bpcnst(ixgraupel), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixgraupel))//' before physics' )
   end if

   call addfld ('CME',        (/ 'lev' /), 'A', 'kg/kg/s',  'Rate of cond-evap within the cloud'                      )
   call addfld ('PRODPREC',   (/ 'lev' /), 'A', 'kg/kg/s',  'Rate of conversion of condensate to precip'              )
   call addfld ('EVAPPREC',   (/ 'lev' /), 'A', 'kg/kg/s',  'Rate of evaporation of falling precip'                   )
   call addfld ('EVAPSNOW',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Rate of evaporation of falling snow'                     )
   call addfld ('HPROGCLD',   (/ 'lev' /), 'A', 'W/kg'    , 'Heating from prognostic clouds'                          )
   call addfld ('FICE',       (/ 'lev' /), 'A', 'fraction', 'Fractional ice content within cloud'                     )
   call addfld ('CLDFSNOW',   (/ 'lev' /), 'A', '1',        'Cloud fraction adjusted for snow'                        )
   call addfld ('ICWMRST',    (/ 'lev' /), 'A', 'kg/kg',    'Prognostic in-stratus water mixing ratio'                )
   call addfld ('ICIMRST',    (/ 'lev' /), 'A', 'kg/kg',    'Prognostic in-stratus ice mixing ratio'                  )

   ! MG microphysics diagnostics
   call addfld ('QCSEVAP',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Rate of evaporation of falling cloud water'              )
   call addfld ('QISEVAP',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Rate of sublimation of falling cloud ice'                )
   call addfld ('QVRES',      (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Rate of residual condensation term'                      )
   call addfld ('CMEIOUT',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Rate of deposition/sublimation of cloud ice'             )
   call addfld ('VTRMC',      (/ 'trop_cld_lev' /), 'A', 'm/s',      'Mass-weighted cloud water fallspeed'                     )
   call addfld ('VTRMI',      (/ 'trop_cld_lev' /), 'A', 'm/s',      'Mass-weighted cloud ice fallspeed'                       )
   call addfld ('QCSEDTEN',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Cloud water mixing ratio tendency from sedimentation'    )
   call addfld ('QISEDTEN',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Cloud ice mixing ratio tendency from sedimentation'      )
   call addfld ('PRAO',       (/ 'lev' /), 'A', 'kg/kg/s',  'Accretion of cloud water by rain'                        )
   call addfld ('PRCO',       (/ 'lev' /), 'A', 'kg/kg/s',  'Autoconversion of cloud water'                           )
   call addfld ('MNUCCCO',    (/ 'lev' /), 'A', 'kg/kg/s',  'Immersion freezing of cloud water'                       )
   call addfld ('MNUCCTO',    (/ 'lev' /), 'A', 'kg/kg/s',  'Contact freezing of cloud water'                         )
   call addfld ('MNUCCDO',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Homogeneous and heterogeneous nucleation from vapor'     )
   call addfld ('MNUCCDOhet', (/ 'lev' /), 'A', 'kg/kg/s',  'Heterogeneous nucleation from vapor'                     )
   call addfld ('MSACWIO',    (/ 'lev' /), 'A', 'kg/kg/s',  'Conversion of cloud water from rime-splintering'         )
   call addfld ('PSACWSO',    (/ 'lev' /), 'A', 'kg/kg/s',  'Accretion of cloud water by snow'                        )
   call addfld ('BERGSO',     (/ 'lev' /), 'A', 'kg/kg/s',  'Conversion of cloud water to snow from bergeron'         )
   call addfld ('BERGO',      (/ 'lev' /), 'A', 'kg/kg/s',  'Conversion of cloud water to cloud ice from bergeron'    )
   call addfld ('MELTO',      (/ 'lev' /), 'A', 'kg/kg/s',  'Melting of cloud ice'                                    )
   call addfld ('MELTSTOT',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Melting of snow'                                    )
   call addfld ('MNUDEPO',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Deposition Nucleation'                                    )
   call addfld ('HOMOO',      (/ 'lev' /), 'A', 'kg/kg/s',  'Homogeneous freezing of cloud water'                     )
   call addfld ('QCRESO',     (/ 'lev' /), 'A', 'kg/kg/s',  'Residual condensation term for cloud water'              )
   call addfld ('PRCIO',      (/ 'lev' /), 'A', 'kg/kg/s',  'Autoconversion of cloud ice to snow'                     )
   call addfld ('PRAIO',      (/ 'lev' /), 'A', 'kg/kg/s',  'Accretion of cloud ice to snow'                          )
   call addfld ('QIRESO',     (/ 'lev' /), 'A', 'kg/kg/s',  'Residual deposition term for cloud ice'                  )
   call addfld ('MNUCCRO',    (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Heterogeneous freezing of rain to snow'                  )
   call addfld ('MNUCCRIO',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Heterogeneous freezing of rain to ice'                  )
   call addfld ('PRACSO',     (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Accretion of rain by snow'                               )
   call addfld ('VAPDEPSO',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Vapor deposition onto snow'                            )
   call addfld ('MELTSDT',    (/ 'trop_cld_lev' /), 'A', 'W/kg',     'Latent heating rate due to melting of snow'              )
   call addfld ('FRZRDT',     (/ 'trop_cld_lev' /), 'A', 'W/kg',     'Latent heating rate due to homogeneous freezing of rain' )
   call addfld ('QRSEDTEN',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s', 'Rain mixing ratio tendency from sedimentation'           )
   call addfld ('QSSEDTEN',   (/ 'trop_cld_lev' /), 'A', 'kg/kg/s', 'Snow mixing ratio tendency from sedimentation'           )
   call addfld ('NNUCCCO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Immersion freezing of cloud water')
   call addfld ('NNUCCTO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Contact freezing of cloud water')
   call addfld ('NNUCCDO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Ice nucleation')
   call addfld ('NNUDEPO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Deposition Nucleation')
   call addfld ('NHOMO',      (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Homogeneous freezing of cloud water')
   call addfld ('NNUCCRO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to heterogeneous freezing of rain to snow')
   call addfld ('NNUCCRIO',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Heterogeneous freezing of rain to ice')
   call addfld ('NSACWIO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Ice Multiplication- Rime-splintering')
   call addfld ('NPRAO',      (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Accretion of cloud water by rain')
   call addfld ('NPSACWSO',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Accretion of cloud water by snow')
   call addfld ('NPRAIO',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Accretion of cloud ice to snow')
   call addfld ('NPRACSO',    (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Accretion of rain by snow')
   call addfld ('NPRCO',      (/ 'trop_cld_lev' /), 'A', '#/kg/s', 'Number Tendency due to Autoconversion of cloud water [to rain]')
   call addfld ('NPRCIO',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Autoconversion of cloud ice to snow')
   call addfld ('NCSEDTEN',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to cloud liquid sedimentation')
   call addfld ('NISEDTEN',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to cloud ice sedimentation')
   call addfld ('NRSEDTEN',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to rain sedimentation')
   call addfld ('NSSEDTEN',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to snow sedimentation')
   call addfld ('NMELTO',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Melting of cloud ice ')
   call addfld ('NMELTS',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Melting of snow')

   if (trim(micro_mg_warm_rain) == 'kk2000') then
      call addfld ('qctend_KK2000',     (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'cloud liquid mass tendency due to autoconversion  & accretion from KK2000')
      call addfld ('nctend_KK2000',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'cloud number mass tendency due to autoconversion  & accretion from KK2000')
      call addfld ('qrtend_KK2000',     (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'rain mass tendency due to autoconversion  & accretion from KK2000')
      call addfld ('nrtend_KK2000',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'rain number tendency due to autoconversion  & accretion from KK2000')
   end if
   if (trim(micro_mg_warm_rain) == 'sb2001') then
      call addfld ('qctend_SB2001',     (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'cloud liquid mass tendency due to autoconversion  & accretion from SB2001')
      call addfld ('nctend_SB2001',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'cloud liquid number tendency due to autoconversion  & accretion from SB2001')
      call addfld ('qrtend_SB2001',     (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'rain mass tendency due to autoconversion  & accretion from SB2001')
      call addfld ('nrtend_SB2001',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'rain number tendency due to autoconversion  & accretion from SB2001')
   end if
   call addfld ('LAMC',    (/ 'trop_cld_lev' /), 'A', 'unitless',  'Size distribution parameter lambda for liquid'     )
   call addfld ('LAMR',    (/ 'trop_cld_lev' /), 'A', 'unitless',  'Size distribution parameter lambda for rain'   )
   call addfld ('PGAM',    (/ 'trop_cld_lev' /), 'A', 'unitless',  'Size distribution parameter mu (pgam) for liquid' )
   call addfld ('N0R',     (/ 'trop_cld_lev' /), 'A', 'unitless',  'Size distribution parameter n0 for rain' )

   if (micro_mg_version > 2) then
         call addfld ('NMELTG',     (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to Melting of graupel')
         call addfld ('NGSEDTEN',   (/ 'trop_cld_lev' /), 'A', '#/kg/s',  'Number Tendency due to graupel sedimentation')
         call addfld ('PSACRO',    (/ 'lev' /), 'A', 'kg/kg/s', 'Collisions between rain & snow (Graupel collecting snow)')
         call addfld ('PRACGO',    (/ 'lev' /), 'A', 'kg/kg/s',  'Change in q collection rain by graupel'             )
         call addfld ('PSACWGO',   (/ 'lev' /), 'A', 'kg/kg/s',  'Change in q collection droplets by graupel'         )
         call addfld ('PGSACWO',   (/ 'lev' /), 'A', 'kg/kg/s',  'Q conversion to graupel due to collection droplets by snow')
         call addfld ('PGRACSO',   (/ 'lev' /), 'A', 'kg/kg/s',  'Q conversion to graupel due to collection rain by snow')
         call addfld ('PRDGO',     (/ 'lev' /), 'A', 'kg/kg/s',  'Deposition of graupel')
         call addfld ('QMULTGO',   (/ 'lev' /), 'A', 'kg/kg/s',  'Q change due to ice mult droplets/graupel')
         call addfld ('QMULTRGO',  (/ 'lev' /), 'A', 'kg/kg/s',  'Q change due to ice mult rain/graupel')
         call addfld ('QGSEDTEN',  (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Graupel/Hail mixing ratio tendency from sedimentation')
         call addfld ('NPRACGO',   (/ 'lev' /), 'A', '#/kg/s',   'Change N collection rain by graupel')
         call addfld ('NSCNGO',    (/ 'lev' /), 'A', '#/kg/s',   'Change N conversion to graupel due to collection droplets by snow')
         call addfld ('NGRACSO',   (/ 'lev' /), 'A', '#/kg/s',   'Change N conversion to graupel due to collection rain by snow')
         call addfld ('NMULTGO',   (/ 'lev' /), 'A', '#/kg/s',  'Ice mult due to acc droplets by graupel ')
         call addfld ('NMULTRGO',  (/ 'lev' /), 'A', '#/kg/s',  'Ice mult due to acc rain by graupel')
         call addfld ('NPSACWGO',  (/ 'lev' /), 'A', '#/kg/s',   'Change N collection droplets by graupel')
         call addfld ('CLDFGRAU',  (/ 'lev' /), 'A', '1',        'Cloud fraction adjusted for graupel'                        )
         call addfld ('MELTGTOT',  (/ 'trop_cld_lev' /), 'A', 'kg/kg/s',  'Melting of graupel'                                    )

   end if

   ! History variables for CAM5 microphysics
   call addfld ('MPDT',       (/ 'lev' /), 'A', 'W/kg',     'Heating tendency - Morrison microphysics'                )
   call addfld ('MPDQ',       (/ 'lev' /), 'A', 'kg/kg/s',  'Q tendency - Morrison microphysics'                      )
   call addfld ('MPDLIQ',     (/ 'lev' /), 'A', 'kg/kg/s',  'CLDLIQ tendency - Morrison microphysics'                 )
   call addfld ('MPDICE',     (/ 'lev' /), 'A', 'kg/kg/s',  'CLDICE tendency - Morrison microphysics'                 )
   call addfld ('MPDNLIQ',    (/ 'lev' /), 'A', '1/kg/s',   'NUMLIQ tendency - Morrison microphysics'                 )
   call addfld ('MPDNICE',    (/ 'lev' /), 'A', '1/kg/s',   'NUMICE tendency - Morrison microphysics'                 )
   call addfld ('MPDW2V',     (/ 'lev' /), 'A', 'kg/kg/s',  'Water <--> Vapor tendency - Morrison microphysics'       )
   call addfld ('MPDW2I',     (/ 'lev' /), 'A', 'kg/kg/s',  'Water <--> Ice tendency - Morrison microphysics'         )
   call addfld ('MPDW2P',     (/ 'lev' /), 'A', 'kg/kg/s',  'Water <--> Precip tendency - Morrison microphysics'      )
   call addfld ('MPDI2V',     (/ 'lev' /), 'A', 'kg/kg/s',  'Ice <--> Vapor tendency - Morrison microphysics'         )
   call addfld ('MPDI2W',     (/ 'lev' /), 'A', 'kg/kg/s',  'Ice <--> Water tendency - Morrison microphysics'         )
   call addfld ('MPDI2P',     (/ 'lev' /), 'A', 'kg/kg/s',  'Ice <--> Precip tendency - Morrison microphysics'        )
   call addfld ('ICWNC',      (/ 'lev' /), 'A', 'm-3',      'Prognostic in-cloud water number conc'                   )
   call addfld ('ICINC',      (/ 'lev' /), 'A', 'm-3',      'Prognostic in-cloud ice number conc'                     )
   call addfld ('EFFLIQ_IND', (/ 'lev' /), 'A','Micron',    'Prognostic droplet effective radius (indirect effect)'   )
   call addfld ('CDNUMC',     horiz_only,  'A', '1/m2',     'Vertically-integrated droplet concentration'             )
   call addfld ('MPICLWPI',   horiz_only,  'A', 'kg/m2',    'Vertically-integrated &
        &in-cloud Initial Liquid WP (Before Micro)' )
   call addfld ('MPICIWPI',   horiz_only,  'A', 'kg/m2',    'Vertically-integrated &
        &in-cloud Initial Ice WP (Before Micro)'    )

   ! This is provided as an example on how to write out subcolumn output
   ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
   if (use_subcol_microp) then
      call addfld('FICE_SCOL', (/'psubcols','lev     '/), 'I', 'fraction', &
           'Sub-column fractional ice content within cloud', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('MPDICE_SCOL', (/'psubcols','lev     '/), 'I', 'kg/kg/s', &
           'Sub-column CLDICE tendency - Morrison microphysics', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('MPDLIQ_SCOL', (/'psubcols','lev     '/), 'I', 'kg/kg/s', &
           'Sub-column CLDLIQ tendency - Morrison microphysics', flag_xyfill=.true., fill_value=1.e30_r8)
   end if


   ! This is only if the coldpoint temperatures are being adjusted.
   ! NOTE: Some fields related to these and output later are added in tropopause.F90.
   if (micro_mg_adjust_cpt) then
     call addfld ('TROPF_TADJ', (/ 'lev' /), 'A', 'K',  'Temperatures after cold point adjustment'                    )
     call addfld ('TROPF_RHADJ', (/ 'lev' /), 'A', 'K', 'Relative Hunidity after cold point adjustment'               )
     call addfld ('TROPF_CDT',   horiz_only,  'A', 'K',  'Cold point temperature adjustment'                           )
     call addfld ('TROPF_CDZ',   horiz_only,  'A', 'm',  'Distance of coldpoint from coldest model level'              )
   end if


   ! Averaging for cloud particle number and size
   call addfld ('AWNC',        (/ 'lev' /),  'A', 'm-3',      'Average cloud water number conc'                                   )
   call addfld ('AWNI',        (/ 'lev' /),  'A', 'm-3',      'Average cloud ice number conc'                                     )
   call addfld ('AREL',        (/ 'lev' /),  'A', 'Micron',   'Average droplet effective radius'                                  )
   call addfld ('AREI',        (/ 'lev' /),  'A', 'Micron',   'Average ice effective radius'                                      )
   ! Frequency arrays for above
   call addfld ('FREQL',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of liquid'                                   )
   call addfld ('FREQI',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of ice'                                      )

   ! Average cloud top particle size and number (liq, ice) and frequency
   call addfld ('ACTREL',      horiz_only,   'A', 'Micron',   'Average Cloud Top droplet effective radius'                        )
   call addfld ('ACTREI',      horiz_only,   'A', 'Micron',   'Average Cloud Top ice effective radius'                            )
   call addfld ('ACTNL',       horiz_only,   'A', 'm-3',   'Average Cloud Top droplet number'                                  )
   call addfld ('ACTNI',       horiz_only,   'A', 'm-3',   'Average Cloud Top ice number'                                      )

   call addfld ('FCTL',        horiz_only,   'A', 'fraction', 'Fractional occurrence of cloud top liquid'                         )
   call addfld ('FCTI',        horiz_only,   'A', 'fraction', 'Fractional occurrence of cloud top ice'                            )

   ! New frequency arrays for mixed phase and supercooled liquid (only and mixed) for (a) Cloud Top and (b) everywhere..
   call addfld ('FREQM',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of mixed phase'                              )
   call addfld ('FREQSL',      (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of only supercooled liquid'                  )
   call addfld ('FREQSLM',     (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of super cooled liquid with ice'             )
   call addfld ('FCTM',        horiz_only,   'A', 'fraction', 'Fractional occurrence of cloud top mixed phase'                    )
   call addfld ('FCTSL',       horiz_only,   'A', 'fraction', 'Fractional occurrence of cloud top only supercooled liquid'        )
   call addfld ('FCTSLM',      horiz_only,   'A', 'fraction', 'Fractional occurrence of cloud top super cooled liquid with ice'   )

   call addfld ('LS_FLXPRC',   (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface rain+snow flux'                         )
   call addfld ('LS_FLXSNW',   (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface snow flux'                              )

   call addfld ('REL',         (/ 'lev' /),  'A', 'micron',   'MG REL stratiform cloud effective radius liquid'                   )
   call addfld ('REI',         (/ 'lev' /),  'A', 'micron',   'MG REI stratiform cloud effective radius ice'                      )
   call addfld ('LS_REFFRAIN', (/ 'lev' /),  'A', 'micron',   'ls stratiform rain effective radius'                               )
   call addfld ('LS_REFFSNOW', (/ 'lev' /),  'A', 'micron',   'ls stratiform snow effective radius'                               )
   call addfld ('CV_REFFLIQ',  (/ 'lev' /),  'A', 'micron',   'convective cloud liq effective radius'                             )
   call addfld ('CV_REFFICE',  (/ 'lev' /),  'A', 'micron',   'convective cloud ice effective radius'                             )
   call addfld ('MG_SADICE',   (/ 'lev' /),  'A', 'cm2/cm3',  'MG surface area density ice'                                       )
   call addfld ('MG_SADSNOW',  (/ 'lev' /),  'A', 'cm2/cm3',  'MG surface area density snow'                                       )

   ! diagnostic precip
   call addfld ('QRAIN',       (/ 'lev' /),  'A', 'kg/kg',    'Diagnostic grid-mean rain mixing ratio'                            )
   call addfld ('QSNOW',       (/ 'lev' /),  'A', 'kg/kg',    'Diagnostic grid-mean snow mixing ratio'                            )
   call addfld ('NRAIN',       (/ 'lev' /),  'A', 'm-3',      'Diagnostic grid-mean rain number conc'                             )
   call addfld ('NSNOW',       (/ 'lev' /),  'A', 'm-3',      'Diagnostic grid-mean snow number conc'                             )

   ! size of precip
   call addfld ('RERCLD',      (/ 'lev' /),  'A', 'm',         'Diagnostic effective radius of Liquid Cloud and Rain'             )
   call addfld ('DSNOW',       (/ 'lev' /),  'A', 'm',         'Diagnostic grid-mean snow diameter'                               )

   ! diagnostic radar reflectivity, cloud-averaged
   call addfld ('REFL',        (/ 'lev' /),  'A', 'DBz',      '94 GHz radar reflectivity'                                         )
   call addfld ('AREFL',       (/ 'lev' /),  'A', 'DBz',      'Average 94 GHz radar reflectivity'                                 )
   call addfld ('FREFL',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of radar reflectivity'                       )

   call addfld ('CSRFL',       (/ 'lev' /),  'A', 'DBz',      '94 GHz radar reflectivity (CloudSat thresholds)'                   )
   call addfld ('ACSRFL',      (/ 'lev' /),  'A', 'DBz',      'Average 94 GHz radar reflectivity (CloudSat thresholds)'           )
   call addfld ('FCSRFL',      (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of radar reflectivity (CloudSat thresholds)' )

   call addfld ('AREFLZ',      (/ 'lev' /),  'A', 'mm^6/m^3', 'Average 94 GHz radar reflectivity'                                 )

   ! 10cm (rain) radar reflectivity
   call addfld ('REFL10CM',    (/ 'lev' /),  'A', 'DBz',      '10cm (Rain) radar reflectivity (Dbz)'                              )
   call addfld ('REFLZ10CM',   (/ 'lev' /),  'A', 'mm^6/m^3', '10cm (Rain) radar reflectivity (Z units)'                          )

   ! Aerosol information
   call addfld ('NCAL',        (/ 'lev' /),  'A', '1/m3',     'Number Concentation Activated for Liquid'                          )
   call addfld ('NCAI',        (/ 'lev' /),  'A', '1/m3',     'Number Concentation Activated for Ice'                             )

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call addfld ('AQRAIN',      (/ 'lev' /),  'A', 'kg/kg',    'Average rain mixing ratio'                                         )
   call addfld ('AQSNOW',      (/ 'lev' /),  'A', 'kg/kg',    'Average snow mixing ratio'                                         )
   call addfld ('ANRAIN',      (/ 'lev' /),  'A', 'm-3',      'Average rain number conc'                                          )
   call addfld ('ANSNOW',      (/ 'lev' /),  'A', 'm-3',      'Average snow number conc'                                          )
   call addfld ('ADRAIN',      (/ 'lev' /),  'A', 'm',        'Average rain effective Diameter'                                   )
   call addfld ('ADSNOW',      (/ 'lev' /),  'A', 'm',        'Average snow effective Diameter'                                   )
   call addfld ('FREQR',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of rain'                                     )
   call addfld ('FREQS',       (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of snow'                                     )

   ! precipitation efficiency & other diagnostic fields
   call addfld('PE'    ,       horiz_only,   'A', '1',        'Stratiform Precipitation Efficiency  (precip/cmeliq)'              )
   call addfld('APRL'  ,       horiz_only,   'A', 'm/s',      'Average Stratiform Precip Rate over efficiency calculation'        )
   call addfld('PEFRAC',       horiz_only,   'A', '1',        'Fraction of timesteps precip efficiency reported'                  )
   call addfld('VPRCO' ,       horiz_only,   'A', 'kg/kg/s',  'Vertical average of autoconversion rate'                           )
   call addfld('VPRAO' ,       horiz_only,   'A', 'kg/kg/s',  'Vertical average of accretion rate'                                )
   call addfld('RACAU' ,       horiz_only,   'A', 'kg/kg/s',  'Accretion/autoconversion ratio from vertical average'              )

   call addfld('UMR', (/ 'trop_cld_lev' /), 'A',   'm/s', 'Mass-weighted rain  fallspeed'              )
   call addfld('UMS', (/ 'trop_cld_lev' /), 'A',   'm/s', 'Mass-weighted snow fallspeed'               )

   if (micro_mg_version > 2) then
      call addfld('UMG',    (/ 'trop_cld_lev' /), 'A',   'm/s', 'Mass-weighted graupel/hail  fallspeed'                    )
      call addfld ('FREQG', (/ 'lev' /),  'A', 'fraction', 'Fractional occurrence of Graupel'                     )
      call addfld ('LS_REFFGRAU', (/ 'lev' /),  'A', 'micron',   'ls stratiform graupel/hail effective radius'    )
      call addfld ('AQGRAU',      (/ 'lev' /),  'A', 'kg/kg',    'Average graupel/hail mixing ratio'              )
      call addfld ('ANGRAU',      (/ 'lev' /),  'A', 'm-3',      'Average graupel/hail number conc'               )
   end if


   ! qc limiter (only output in versions 1.5 and later)
   call addfld('QCRAT', (/ 'lev' /), 'A', 'fraction', 'Qc Limiter: Fraction of qc tendency applied')

   ! determine the add_default fields
   call phys_getopts(history_amwg_out           = history_amwg         , &
                     history_budget_out         = history_budget       , &
                     history_budget_histfile_num_out = budget_histfile)

   if (history_amwg) then
      call add_default ('FICE    ', 1, ' ')
      call add_default ('AQRAIN   ', 1, ' ')
      call add_default ('AQSNOW   ', 1, ' ')
      call add_default ('ANRAIN   ', 1, ' ')
      call add_default ('ANSNOW   ', 1, ' ')
      call add_default ('ADRAIN   ', 1, ' ')
      call add_default ('ADSNOW   ', 1, ' ')
      call add_default ('AREI     ', 1, ' ')
      call add_default ('AREL     ', 1, ' ')
      call add_default ('AWNC     ', 1, ' ')
      call add_default ('AWNI     ', 1, ' ')
      call add_default ('CDNUMC   ', 1, ' ')
      call add_default ('FREQR    ', 1, ' ')
      call add_default ('FREQS    ', 1, ' ')
      call add_default ('FREQL    ', 1, ' ')
      call add_default ('FREQI    ', 1, ' ')
      do m = 1, ncnst
         call cnst_get_ind(cnst_names(m), mm)
         call add_default(cnst_name(mm), 1, ' ')
      end do
   end if

   if ( history_budget ) then
      call add_default ('EVAPSNOW ', budget_histfile, ' ')
      call add_default ('EVAPPREC ', budget_histfile, ' ')
      call add_default ('QVRES    ', budget_histfile, ' ')
      call add_default ('QISEVAP  ', budget_histfile, ' ')
      call add_default ('QCSEVAP  ', budget_histfile, ' ')
      call add_default ('QISEDTEN ', budget_histfile, ' ')
      call add_default ('QCSEDTEN ', budget_histfile, ' ')
      call add_default ('QIRESO   ', budget_histfile, ' ')
      call add_default ('QCRESO   ', budget_histfile, ' ')
      call add_default ('QRSEDTEN ', budget_histfile, ' ')
      call add_default ('QSSEDTEN ', budget_histfile, ' ')
      call add_default ('PSACWSO  ', budget_histfile, ' ')
      call add_default ('PRCO     ', budget_histfile, ' ')
      call add_default ('PRCIO    ', budget_histfile, ' ')
      call add_default ('PRAO     ', budget_histfile, ' ')
      call add_default ('PRAIO    ', budget_histfile, ' ')
      call add_default ('PRACSO   ', budget_histfile, ' ')
      call add_default ('VAPDEPSO ', budget_histfile, ' ')
      call add_default ('MSACWIO  ', budget_histfile, ' ')
      call add_default ('MPDW2V   ', budget_histfile, ' ')
      call add_default ('MPDW2P   ', budget_histfile, ' ')
      call add_default ('MPDW2I   ', budget_histfile, ' ')
      call add_default ('MPDT     ', budget_histfile, ' ')
      call add_default ('MPDQ     ', budget_histfile, ' ')
      call add_default ('MPDLIQ   ', budget_histfile, ' ')
      call add_default ('MPDICE   ', budget_histfile, ' ')
      call add_default ('MPDI2W   ', budget_histfile, ' ')
      call add_default ('MPDI2V   ', budget_histfile, ' ')
      call add_default ('MPDI2P   ', budget_histfile, ' ')
      call add_default ('MNUCCTO  ', budget_histfile, ' ')
      call add_default ('MNUCCRO  ', budget_histfile, ' ')
      call add_default ('MNUCCRIO ', budget_histfile, ' ')
      call add_default ('MNUCCCO  ', budget_histfile, ' ')
      call add_default ('MELTSDT  ', budget_histfile, ' ')
      call add_default ('MELTO    ', budget_histfile, ' ')
      call add_default ('HOMOO    ', budget_histfile, ' ')
      call add_default ('FRZRDT   ', budget_histfile, ' ')
      call add_default ('CMEIOUT  ', budget_histfile, ' ')
      call add_default ('BERGSO   ', budget_histfile, ' ')
      call add_default ('BERGO    ', budget_histfile, ' ')
      call add_default ('MELTSTOT ', budget_histfile, ' ')
      call add_default ('MNUDEPO  ', budget_histfile, ' ')
      call add_default ('NNUCCCO  ', budget_histfile, ' ')
      call add_default ('NNUCCTO  ', budget_histfile, ' ')
      call add_default ('NNUCCDO  ', budget_histfile, ' ')
      call add_default ('NNUDEPO  ', budget_histfile, ' ')
      call add_default ('NHOMO    ', budget_histfile, ' ')
      call add_default ('NNUCCRO  ', budget_histfile, ' ')
      call add_default ('NNUCCRIO ', budget_histfile, ' ')
      call add_default ('NSACWIO  ', budget_histfile, ' ')
      call add_default ('NPRAO    ', budget_histfile, ' ')
      call add_default ('NPSACWSO ', budget_histfile, ' ')
      call add_default ('NPRAIO   ', budget_histfile, ' ')
      call add_default ('NPRACSO  ', budget_histfile, ' ')
      call add_default ('NPRCO    ', budget_histfile, ' ')
      call add_default ('NPRCIO   ', budget_histfile, ' ')
      call add_default ('NCSEDTEN ', budget_histfile, ' ')
      call add_default ('NISEDTEN ', budget_histfile, ' ')
      call add_default ('NRSEDTEN ', budget_histfile, ' ')
      call add_default ('NSSEDTEN ', budget_histfile, ' ')
      call add_default ('NMELTO   ', budget_histfile, ' ')
      call add_default ('NMELTS   ', budget_histfile, ' ')
      call add_default ('NCAL     ', budget_histfile, ' ')
      if (micro_mg_version > 2) then
         call add_default ('QGSEDTEN ', budget_histfile, ' ')
         call add_default ('PSACRO    ', budget_histfile, ' ')
         call add_default ('PRACGO    ', budget_histfile, ' ')
         call add_default ('PSACWGO   ', budget_histfile, ' ')
         call add_default ('PGSACWO   ', budget_histfile, ' ')
         call add_default ('PGRACSO   ', budget_histfile, ' ')
         call add_default ('PRDGO     ', budget_histfile, ' ')
         call add_default ('QMULTGO   ', budget_histfile, ' ')
         call add_default ('QMULTRGO  ', budget_histfile, ' ')
         call add_default ('MELTGTOT  ', budget_histfile, ' ')
         call add_default ('NPRACGO   ', budget_histfile, ' ')
         call add_default ('NSCNGO    ', budget_histfile, ' ')
         call add_default ('NGRACSO   ', budget_histfile, ' ')
         call add_default ('NMULTGO  ', budget_histfile, ' ')
         call add_default ('NMULTRGO  ', budget_histfile, ' ')
         call add_default ('NPSACWGO  ', budget_histfile, ' ')
         call add_default ('NGSEDTEN ', budget_histfile, ' ')
         call add_default ('NMELTG   ', budget_histfile, ' ')
      end if
      call add_default(cnst_name(ixcldliq), budget_histfile, ' ')
      call add_default(cnst_name(ixcldice), budget_histfile, ' ')
      call add_default(apcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(apcnst   (ixcldice), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldice), budget_histfile, ' ')
      call add_default(cnst_name(ixrain), budget_histfile, ' ')
      call add_default(cnst_name(ixsnow), budget_histfile, ' ')
      call add_default(apcnst   (ixrain), budget_histfile, ' ')
      call add_default(apcnst   (ixsnow), budget_histfile, ' ')
      call add_default(bpcnst   (ixrain), budget_histfile, ' ')
      call add_default(bpcnst   (ixsnow), budget_histfile, ' ')

      if (micro_mg_version > 2) then
         call add_default(cnst_name(ixgraupel), budget_histfile, ' ')
         call add_default(apcnst   (ixgraupel), budget_histfile, ' ')
         call add_default(bpcnst   (ixgraupel), budget_histfile, ' ')
      end if

   end if

   ! physics buffer indices
   ast_idx      = pbuf_get_index('AST')
   cld_idx      = pbuf_get_index('CLD')
   concld_idx   = pbuf_get_index('CONCLD')

   naai_idx     = pbuf_get_index('NAAI')
   naai_hom_idx = pbuf_get_index('NAAI_HOM')
   npccn_idx    = pbuf_get_index('NPCCN')
   rndst_idx    = pbuf_get_index('RNDST')
   nacon_idx    = pbuf_get_index('NACON')

   prec_str_idx = pbuf_get_index('PREC_STR')
   snow_str_idx = pbuf_get_index('SNOW_STR')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')

   cmeliq_idx = pbuf_get_index('CMELIQ')

   ! These fields may have been added, so don't abort if they have not been
   qsatfac_idx  = pbuf_get_index('QSATFAC', ierr)
   qrain_idx    = pbuf_get_index('QRAIN', ierr)
   qsnow_idx    = pbuf_get_index('QSNOW', ierr)
   nrain_idx    = pbuf_get_index('NRAIN', ierr)
   nsnow_idx    = pbuf_get_index('NSNOW', ierr)

  ! fields for heterogeneous freezing
  frzimm_idx = pbuf_get_index('FRZIMM', ierr)
  frzcnt_idx = pbuf_get_index('FRZCNT', ierr)
  frzdep_idx = pbuf_get_index('FRZDEP', ierr)

  ! Initialize physics buffer grid fields for accumulating precip and condensation
   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8)
      call pbuf_set_field(pbuf2d, acpr_idx,   0._r8)
      call pbuf_set_field(pbuf2d, acgcme_idx, 0._r8)
      call pbuf_set_field(pbuf2d, acnum_idx,  0)
      call pbuf_set_field(pbuf2d, relvar_idx, 2._r8)
      call pbuf_set_field(pbuf2d, accre_enhan_idx, 1._r8)
      call pbuf_set_field(pbuf2d, am_evp_st_idx,  0._r8)
      call pbuf_set_field(pbuf2d, evprain_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, evpsnow_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, prer_evap_idx,  0._r8)
      call pbuf_set_field(pbuf2d, bergso_idx, 0._r8)
      call pbuf_set_field(pbuf2d, icswp_idx, 0._r8)
      call pbuf_set_field(pbuf2d, cldfsnow_idx, 0._r8)
      call pbuf_set_field(pbuf2d, dei_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, des_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, mu_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, lambdac_idx, 0.0_r8)

      if (degrau_idx > 0) call pbuf_set_field(pbuf2d, degrau_idx, 0.0_r8)
      if (icgrauwp_idx > 0) call pbuf_set_field(pbuf2d, icgrauwp_idx, 0.0_r8)
      if (qrain_idx > 0)   call pbuf_set_field(pbuf2d, qrain_idx, 0._r8)
      if (qsnow_idx > 0)   call pbuf_set_field(pbuf2d, qsnow_idx, 0._r8)
      if (nrain_idx > 0)   call pbuf_set_field(pbuf2d, nrain_idx, 0._r8)
      if (nsnow_idx > 0)   call pbuf_set_field(pbuf2d, nsnow_idx, 0._r8)
      if (qcsedten_idx > 0)   call pbuf_set_field(pbuf2d, qcsedten_idx, 0._r8)
      if (qrsedten_idx > 0)   call pbuf_set_field(pbuf2d, qrsedten_idx, 0._r8)
      if (qisedten_idx > 0)   call pbuf_set_field(pbuf2d, qisedten_idx, 0._r8)
      if (qssedten_idx > 0)   call pbuf_set_field(pbuf2d, qssedten_idx, 0._r8)
      if (vtrmc_idx > 0)      call pbuf_set_field(pbuf2d, vtrmc_idx, 0._r8)
      if (umr_idx > 0)        call pbuf_set_field(pbuf2d, umr_idx, 0._r8)
      if (vtrmi_idx > 0)      call pbuf_set_field(pbuf2d, vtrmi_idx, 0._r8)
      if (ums_idx > 0)        call pbuf_set_field(pbuf2d, ums_idx, 0._r8)
      if (qcsevap_idx > 0)    call pbuf_set_field(pbuf2d, qcsevap_idx, 0._r8)
      if (qisevap_idx > 0)    call pbuf_set_field(pbuf2d, qisevap_idx, 0._r8)

      ! If sub-columns turned on, need to set the sub-column fields as well
      if (use_subcol_microp) then
         call pbuf_set_field(pbuf2d, cldo_idx,    0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_t_idx,    0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qv_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ql_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qi_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_nl_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ni_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qlst_idx, 0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, icswp_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cldfsnow_idx,0._r8, col_type=col_type_subcol)
      end if

   end if

end subroutine micro_pumas_cam_init

!===============================================================================

subroutine micro_pumas_cam_tend(state, ptend, dtime, pbuf)

   use micro_pumas_utils, only: size_dist_param_basic, size_dist_param_liq
   use micro_pumas_utils, only: mg_liq_props, mg_ice_props, avg_diameter
   use micro_pumas_utils, only: rhoi, rhosn, rhow, rhows, rhog, qsmall, mincld

   use micro_pumas_v1,    only: micro_pumas_tend

   use physics_buffer,  only: pbuf_col_type_index
   use subcol,          only: subcol_field_avg
   use tropopause,      only: tropopause_find, TROP_ALG_CPP, TROP_ALG_NONE, NOTFOUND
   use wv_saturation,   only: qsat
   use infnan,          only: nan, assignment(=)
   use cam_abortutils,  only: handle_allocate_error

   use stochastic_tau_cam, only: ncd

   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtime
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! Local variables

   type(proc_rates_type) :: proc_rates

   integer :: lchnk, ncol, psetcols, ngrdcol

   integer :: i, k, itim_old, it

   real(r8), parameter :: micron2meter = 1.e6_r8
   real(r8), parameter :: shapeparam = 1.e5_r8

   real(r8), pointer :: naai(:,:)      ! ice nucleation number
   real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
   real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
   real(r8), pointer :: rndst(:,:,:)
   real(r8), pointer :: nacon(:,:,:)
   real(r8), pointer :: am_evp_st_grid(:,:)    ! Evaporation area of stratiform precipitation. 0<= am_evp_st <=1.
   real(r8), pointer :: evprain_st_grid(:,:)   ! Evaporation rate of stratiform rain [kg/kg/s]
   real(r8), pointer :: evpsnow_st_grid(:,:)   ! Evaporation rate of stratiform snow [kg/kg/s]

   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

   real(r8), pointer :: ast(:,:)          ! Relative humidity cloud fraction
   real(r8), pointer :: qsatfac(:,:)      ! Subgrid cloud water saturation scaling factor.
   real(r8), pointer :: alst_mic(:,:)
   real(r8), pointer :: aist_mic(:,:)
   real(r8), pointer :: cldo(:,:)         ! Old cloud fraction
   real(r8), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
   real(r8), pointer :: prer_evap(:,:)    ! precipitation evaporation rate
   real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
   real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation
   real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
   real(r8), pointer :: dei(:,:)          ! Ice effective diameter (meters)
   real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
   real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
   real(r8), pointer :: des(:,:)          ! Snow effective diameter (m)
   real(r8), pointer :: degrau(:,:)       ! Graupel effective diameter (m)
   real(r8), pointer :: bergstot(:,:)     ! Conversion of cloud water to snow from bergeron

   real(r8) :: rho(state%psetcols,pver)
   real(r8) :: cldmax(state%psetcols,pver)

   real(r8)  :: rate1cld(state%psetcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

   real(r8)  :: tlat(state%psetcols,pver)
   real(r8)  :: qvlat(state%psetcols,pver)
   real(r8)  :: qcten(state%psetcols,pver)
   real(r8)  :: qiten(state%psetcols,pver)
   real(r8)  :: ncten(state%psetcols,pver)
   real(r8)  :: niten(state%psetcols,pver)

   real(r8)  :: qrten(state%psetcols,pver)
   real(r8)  :: qsten(state%psetcols,pver)
   real(r8)  :: nrten(state%psetcols,pver)
   real(r8)  :: nsten(state%psetcols,pver)
   real(r8)  :: qgten(state%psetcols,pver)
   real(r8)  :: ngten(state%psetcols,pver)

   real(r8)  :: prect(state%psetcols)
   real(r8)  :: preci(state%psetcols)
   real(r8)  :: am_evp_st(state%psetcols,pver)  ! Area over which precip evaporates
   real(r8)  :: cmeice(state%psetcols,pver)     ! Rate of cond-evap of ice within the cloud
   real(r8)  :: qsout(state%psetcols,pver)      ! Snow mixing ratio
   real(r8)  :: cflx(state%psetcols,pverp)      ! grid-box avg liq condensate flux (kg m^-2 s^-1)
   real(r8)  :: iflx(state%psetcols,pverp)      ! grid-box avg ice condensate flux (kg m^-2 s^-1)
   real(r8)  :: rflx(state%psetcols,pverp)      ! grid-box average rain flux (kg m^-2 s^-1)
   real(r8)  :: sflx(state%psetcols,pverp)      ! grid-box average snow flux (kg m^-2 s^-1)
   real(r8)  :: gflx(state%psetcols,pverp)      ! grid-box average snow flux (kg m^-2 s^-1)
   real(r8)  :: qrout(state%psetcols,pver)      ! Rain mixing ratio

   real(r8)  :: nrout(state%psetcols,pver)
   real(r8)  :: nsout(state%psetcols,pver)
   real(r8)  :: refl(state%psetcols,pver)    ! analytic radar reflectivity
   real(r8)  :: arefl(state%psetcols,pver)   ! average reflectivity will zero points outside valid range
   real(r8)  :: areflz(state%psetcols,pver)  ! average reflectivity in z.
   real(r8)  :: frefl(state%psetcols,pver)
   real(r8)  :: csrfl(state%psetcols,pver)   ! cloudsat reflectivity
   real(r8)  :: acsrfl(state%psetcols,pver)  ! cloudsat average
   real(r8)  :: fcsrfl(state%psetcols,pver)
   real(r8)  :: refl10cm(state%psetcols,pver)    ! analytic radar reflectivity
   real(r8)  :: reflz10cm(state%psetcols,pver)    ! analytic radar reflectivity Z
   real(r8)  :: rercld(state%psetcols,pver)  ! effective radius calculation for rain + cloud
   real(r8)  :: ncai(state%psetcols,pver)    ! output number conc of ice nuclei available (1/m3)
   real(r8)  :: ncal(state%psetcols,pver)    ! output number conc of CCN (1/m3)
   real(r8)  :: qrout2(state%psetcols,pver)
   real(r8)  :: qsout2(state%psetcols,pver)
   real(r8)  :: nrout2(state%psetcols,pver)
   real(r8)  :: nsout2(state%psetcols,pver)
   real(r8)  :: freqs(state%psetcols,pver)
   real(r8)  :: freqr(state%psetcols,pver)
   real(r8)  :: nfice(state%psetcols,pver)
   real(r8)  :: qcrat(state%psetcols,pver)   ! qc limiter ratio (1=no limit)

!Hail/Graupel Output
   real(r8)  :: freqg(state%psetcols,pver)
   real(r8)  :: qgout(state%psetcols,pver)
   real(r8)  :: ngout(state%psetcols,pver)
   real(r8)  :: dgout(state%psetcols,pver)
   real(r8)  :: qgout2(state%psetcols,pver)
   real(r8)  :: ngout2(state%psetcols,pver)
   real(r8)  :: dgout2(state%psetcols,pver)

   ! Dummy arrays for cases where we throw away the MG version and
   ! recalculate sizes on the CAM grid to avoid time/subcolumn averaging
   ! issues.
   real(r8) :: rel_fn_dum(state%ncol,pver)
   real(r8) :: dsout2_dum(state%ncol,pver)
   real(r8) :: drout_dum(state%ncol,pver)
   real(r8) :: reff_rain_dum(state%ncol,pver)
   real(r8) :: reff_snow_dum(state%ncol,pver)
   real(r8) :: reff_grau_dum(state%ncol,pver)   !not used for now or passed to COSP.
   real(r8), target :: nan_array(state%ncol,pver)   ! Array for NaN's

   ! Heterogeneous-only version of mnuccdtot.
   real(r8) :: mnuccdohet(state%psetcols,pver)

   ! physics buffer fields for COSP simulator
   real(r8), pointer :: mgflxprc(:,:)     ! MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgflxsnw(:,:)     ! MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgmrprc(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
   real(r8), pointer :: mgmrsnw(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
   real(r8), pointer :: mgreffrain_grid(:,:)   ! MG diagnostic rain effective radius (um)
   real(r8), pointer :: mgreffsnow_grid(:,:)   ! MG diagnostic snow effective radius (um)
   real(r8), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
   real(r8), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)

   ! physics buffer fields used with CARMA
   real(r8), pointer, dimension(:,:) :: tnd_qsnow    ! external tendency on snow mass (kg/kg/s)
   real(r8), pointer, dimension(:,:) :: tnd_nsnow    ! external tendency on snow number(#/kg/s)
   real(r8), pointer, dimension(:,:) :: re_ice       ! ice effective radius (m)

   real(r8), pointer :: rate1ord_cw2pr_st(:,:) ! 1st order rate for direct conversion of
                                               ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
   real(r8), pointer :: wsedl(:,:)        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]


   real(r8), pointer :: CC_T(:,:)         ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qv(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ql(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qi(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_nl(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ni(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qlst(:,:)      ! In-liquid stratus microphysical tendency

   ! variables for heterogeneous freezing
   real(r8), pointer :: frzimm(:,:)
   real(r8), pointer :: frzcnt(:,:)
   real(r8), pointer :: frzdep(:,:)

   real(r8), pointer :: qme(:,:)

   ! A local copy of state is used for diagnostic calculations
   type(physics_state) :: state_loc
   type(physics_ptend) :: ptend_loc

   real(r8) :: icecldf(state%psetcols,pver) ! Ice cloud fraction
   real(r8) :: liqcldf(state%psetcols,pver) ! Liquid cloud fraction (combined into cloud)

   real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
   real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)
   real(r8), pointer :: sadice(:,:)       ! Ice surface area density (cm2/cm3)
   real(r8), pointer :: sadsnow(:,:)      ! Snow surface area density (cm2/cm3)


   real(r8), pointer :: cmeliq(:,:)

   real(r8), pointer :: cld(:,:)          ! Total cloud fraction
   real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
   real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
   real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation
   real(r8), pointer :: cldfsnow(:,:)     ! Cloud fraction for liquid+snow
   real(r8), pointer :: icswp(:,:)        ! In-cloud snow water path

   real(r8), pointer :: cldfgrau(:,:)     ! Cloud fraction for liquid+snow
   real(r8), pointer :: icgrauwp(:,:)        ! In-cloud snow water path

   real(r8) :: icimrst(state%psetcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst(state%psetcols,pver) ! In stratus water mixing ratio
   real(r8) :: icinc(state%psetcols,pver)   ! In cloud ice number conc
   real(r8) :: icwnc(state%psetcols,pver)   ! In cloud water number conc

   real(r8) :: iclwpi(state%psetcols)       ! Vertically-integrated in-cloud Liquid WP before microphysics
   real(r8) :: iciwpi(state%psetcols)       ! Vertically-integrated in-cloud Ice WP before microphysics

   ! Averaging arrays for effective radius and number....
   real(r8) :: efiout_grid(pcols,pver)
   real(r8) :: efcout_grid(pcols,pver)
   real(r8) :: ncout_grid(pcols,pver)
   real(r8) :: niout_grid(pcols,pver)
   real(r8) :: freqi_grid(pcols,pver)
   real(r8) :: freql_grid(pcols,pver)

!  Averaging arrays for supercooled liquid
   real(r8) :: freqm_grid(pcols,pver)
   real(r8) :: freqsl_grid(pcols,pver)
   real(r8) :: freqslm_grid(pcols,pver)
   real(r8) :: fctm_grid(pcols)
   real(r8) :: fctsl_grid(pcols)
   real(r8) :: fctslm_grid(pcols)

   real(r8) :: cdnumc_grid(pcols)           ! Vertically-integrated droplet concentration
   real(r8) :: icimrst_grid_out(pcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst_grid_out(pcols,pver) ! In stratus water mixing ratio

   ! Cloud fraction used for precipitation.
   real(r8) :: cldmax_grid(pcols,pver)

   ! Average cloud top radius & number
   real(r8) :: ctrel_grid(pcols)
   real(r8) :: ctrei_grid(pcols)
   real(r8) :: ctnl_grid(pcols)
   real(r8) :: ctni_grid(pcols)
   real(r8) :: fcti_grid(pcols)
   real(r8) :: fctl_grid(pcols)

   real(r8) :: ftem_grid(pcols,pver)

   ! Variables for precip efficiency calculation
   real(r8) :: minlwp        ! LWP threshold

   real(r8), pointer, dimension(:) :: acprecl_grid ! accumulated precip across timesteps
   real(r8), pointer, dimension(:) :: acgcme_grid  ! accumulated condensation across timesteps
   integer,  pointer, dimension(:) :: acnum_grid   ! counter for # timesteps accumulated

   ! Variables for liquid water path and column condensation
   real(r8) :: tgliqwp_grid(pcols)   ! column liquid
   real(r8) :: tgcmeliq_grid(pcols)  ! column condensation rate (units)

   real(r8) :: pe_grid(pcols)        ! precip efficiency for output
   real(r8) :: pefrac_grid(pcols)    ! fraction of time precip efficiency is written out
   real(r8) :: tpr_grid(pcols)       ! average accumulated precipitation rate in pe calculation

   ! variables for autoconversion and accretion vertical averages
   real(r8) :: vprco_grid(pcols)     ! vertical average autoconversion
   real(r8) :: vprao_grid(pcols)     ! vertical average accretion
   real(r8) :: racau_grid(pcols)     ! ratio of vertical averages
   integer  :: cnt_grid(pcols)       ! counters

   logical  :: lq(pcnst)

   real(r8) :: icimrst_grid(pcols,pver) ! stratus ice mixing ratio - on grid
   real(r8) :: icwmrst_grid(pcols,pver) ! stratus water mixing ratio - on grid

   real(r8), pointer :: lambdac_grid(:,:)
   real(r8), pointer :: mu_grid(:,:)
   real(r8), pointer :: rel_grid(:,:)
   real(r8), pointer :: rei_grid(:,:)
   real(r8), pointer :: sadice_grid(:,:)
   real(r8), pointer :: sadsnow_grid(:,:)
   real(r8), pointer :: dei_grid(:,:)
   real(r8), pointer :: des_grid(:,:)
   real(r8), pointer :: iclwpst_grid(:,:)
   real(r8), pointer :: degrau_grid(:,:)

   real(r8) :: rho_grid(pcols,pver)
   real(r8) :: liqcldf_grid(pcols,pver)
   real(r8) :: qsout_grid(pcols,pver)
   real(r8) :: ncic_grid(pcols,pver)
   real(r8) :: niic_grid(pcols,pver)
   real(r8) :: rel_fn_grid(pcols,pver)    ! Ice effective drop size at fixed number (indirect effect) (microns) - on grid
   real(r8) :: qrout_grid(pcols,pver)
   real(r8) :: drout2_grid(pcols,pver)
   real(r8) :: dsout2_grid(pcols,pver)
   real(r8) :: nsout_grid(pcols,pver)
   real(r8) :: nrout_grid(pcols,pver)
   real(r8) :: reff_rain_grid(pcols,pver)
   real(r8) :: reff_snow_grid(pcols,pver)
   real(r8) :: reff_grau_grid(pcols,pver)
   real(r8) :: cld_grid(pcols,pver)
   real(r8) :: pdel_grid(pcols,pver)
   real(r8) :: prco_grid(pcols,pver)
   real(r8) :: prao_grid(pcols,pver)
   real(r8) :: icecldf_grid(pcols,pver)
   real(r8) :: icwnc_grid(pcols,pver)
   real(r8) :: icinc_grid(pcols,pver)
   real(r8) :: qcreso_grid(pcols,pver)
   real(r8) :: melto_grid(pcols,pver)
   real(r8) :: mnuccco_grid(pcols,pver)
   real(r8) :: mnuccto_grid(pcols,pver)
   real(r8) :: bergo_grid(pcols,pver)
   real(r8) :: homoo_grid(pcols,pver)
   real(r8) :: msacwio_grid(pcols,pver)
   real(r8) :: psacwso_grid(pcols,pver)
   real(r8) :: cmeiout_grid(pcols,pver)
   real(r8) :: qireso_grid(pcols,pver)
   real(r8) :: prcio_grid(pcols,pver)
   real(r8) :: praio_grid(pcols,pver)
   real(r8) :: psacro_grid(pcols,pver)
   real(r8) :: pracgo_grid(pcols,pver)
   real(r8) :: psacwgo_grid(pcols,pver)
   real(r8) :: pgsacwo_grid(pcols,pver)
   real(r8) :: pgracso_grid(pcols,pver)
   real(r8) :: prdgo_grid(pcols,pver)
   real(r8) :: qmultgo_grid(pcols,pver)
   real(r8) :: qmultrgo_grid(pcols,pver)
   real(r8) :: npracgo_grid(pcols,pver)
   real(r8) :: nscngo_grid(pcols,pver)
   real(r8) :: ngracso_grid(pcols,pver)
   real(r8) :: nmultgo_grid(pcols,pver)
   real(r8) :: nmultrgo_grid(pcols,pver)
   real(r8) :: npsacwgo_grid(pcols,pver)
   real(r8) :: qcsedtenout_grid(pcols,pver)
   real(r8) :: qrsedtenout_grid(pcols,pver)
   real(r8) :: qisedtenout_grid(pcols,pver)
   real(r8) :: qssedtenout_grid(pcols,pver)
   real(r8) :: vtrmcout_grid(pcols,pver)
   real(r8) :: umrout_grid(pcols,pver)
   real(r8) :: vtrmiout_grid(pcols,pver)
   real(r8) :: umsout_grid(pcols,pver)
   real(r8) :: qcsevapout_grid(pcols,pver)
   real(r8) :: qisevapout_grid(pcols,pver)

   real(r8) :: nc_grid(pcols,pver)
   real(r8) :: ni_grid(pcols,pver)
   real(r8) :: qr_grid(pcols,pver)
   real(r8) :: nr_grid(pcols,pver)
   real(r8) :: qs_grid(pcols,pver)
   real(r8) :: ns_grid(pcols,pver)
   real(r8) :: qg_grid(pcols,pver)
   real(r8) :: ng_grid(pcols,pver)

   real(r8) :: dgout2_grid(pcols,pver)

   real(r8) :: cp_rh(pcols,pver)
   real(r8) :: cp_t(pcols)
   real(r8) :: cp_z(pcols)
   real(r8) :: cp_dt(pcols)
   real(r8) :: cp_dz(pcols)
   integer  :: troplev(pcols)
   real(r8) :: es
   real(r8) :: qs

   real(r8) :: state_loc_graup(state%psetcols,pver)
   real(r8) :: state_loc_numgraup(state%psetcols,pver)

   real(r8), pointer :: cmeliq_grid(:,:)

   real(r8), pointer :: prec_str_grid(:)
   real(r8), pointer :: snow_str_grid(:)
   real(r8), pointer :: prec_pcw_grid(:)
   real(r8), pointer :: snow_pcw_grid(:)
   real(r8), pointer :: prec_sed_grid(:)
   real(r8), pointer :: snow_sed_grid(:)
   real(r8), pointer :: cldo_grid(:,:)
   real(r8), pointer :: nevapr_grid(:,:)
   real(r8), pointer :: prain_grid(:,:)
   real(r8), pointer :: mgflxprc_grid(:,:)
   real(r8), pointer :: mgflxsnw_grid(:,:)
   real(r8), pointer :: mgmrprc_grid(:,:)
   real(r8), pointer :: mgmrsnw_grid(:,:)
   real(r8), pointer :: cvreffliq_grid(:,:)
   real(r8), pointer :: cvreffice_grid(:,:)
   real(r8), pointer :: rate1ord_cw2pr_st_grid(:,:)
   real(r8), pointer :: wsedl_grid(:,:)
   real(r8), pointer :: CC_t_grid(:,:)
   real(r8), pointer :: CC_qv_grid(:,:)
   real(r8), pointer :: CC_ql_grid(:,:)
   real(r8), pointer :: CC_qi_grid(:,:)
   real(r8), pointer :: CC_nl_grid(:,:)
   real(r8), pointer :: CC_ni_grid(:,:)
   real(r8), pointer :: CC_qlst_grid(:,:)
   real(r8), pointer :: qme_grid(:,:)
   real(r8), pointer :: iciwpst_grid(:,:)
   real(r8), pointer :: icswp_grid(:,:)
   real(r8), pointer :: ast_grid(:,:)
   real(r8), pointer :: cldfsnow_grid(:,:)
   real(r8), pointer :: bergso_grid(:,:)

   real(r8), pointer :: icgrauwp_grid(:,:)
   real(r8), pointer :: cldfgrau_grid(:,:)

   real(r8), pointer :: qrout_grid_ptr(:,:)
   real(r8), pointer :: qsout_grid_ptr(:,:)
   real(r8), pointer :: nrout_grid_ptr(:,:)
   real(r8), pointer :: nsout_grid_ptr(:,:)
   real(r8), pointer :: qcsedtenout_grid_ptr(:,:)
   real(r8), pointer :: qrsedtenout_grid_ptr(:,:)
   real(r8), pointer :: qisedtenout_grid_ptr(:,:)
   real(r8), pointer :: qssedtenout_grid_ptr(:,:)
   real(r8), pointer :: vtrmcout_grid_ptr(:,:)
   real(r8), pointer :: umrout_grid_ptr(:,:)
   real(r8), pointer :: vtrmiout_grid_ptr(:,:)
   real(r8), pointer :: umsout_grid_ptr(:,:)
   real(r8), pointer :: qcsevapout_grid_ptr(:,:)
   real(r8), pointer :: qisevapout_grid_ptr(:,:)


   logical :: use_subcol_microp
   integer :: col_type ! Flag to store whether accessing grid or sub-columns in pbuf_get_field
   integer :: ierr
   integer :: nlev

   character(128) :: errstring   ! return status (non-blank for error return)

   ! For rrtmg optics. specified distribution.
   real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (meters)
   real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
   real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (meters)

   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   psetcols = state%psetcols
   ngrdcol  = state%ngrdcol
   itim_old = pbuf_old_tim_idx()
   nlev = pver - top_lev + 1

   nan_array = nan

   ! Allocate the proc_rates DDT
   ! IMPORTANT NOTE -- elements in proc_rates are dimensioned to the nlev dimension while
   !     all the other arrays in this routine are dimensioned pver.  This is required because
   !     PUMAS only gets the top_lev:pver array subsection, and the proc_rates arrays
   !     need to be the same levels.
   call proc_rates%allocate(ncol, nlev, ncd, micro_mg_warm_rain, errstring)

   call handle_errmsg(errstring, subname="micro_pumas_cam_tend")


   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

   ! Set the col_type flag to grid or subcolumn dependent on the value of use_subcol_microp
   call pbuf_col_type_index(use_subcol_microp, col_type=col_type)

   !-----------------------
   ! These physics buffer fields are read only and not set in this parameterization
   ! If these fields do not have subcolumn data, copy the grid to the subcolumn if subcolumns is turned on
   ! If subcolumns is not turned on, then these fields will be grid data

   call pbuf_get_field(pbuf, naai_idx,        naai,        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, naai_hom_idx,    naai_hom,    col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, npccn_idx,       npccn,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, rndst_idx,       rndst,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, nacon_idx,       nacon,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, relvar_idx,      relvar,      col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan, col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq,      col_type=col_type, copy_if_needed=use_subcol_microp)

   call pbuf_get_field(pbuf, cld_idx,         cld,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, concld_idx,      concld,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, ast_idx,         ast,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)

   if (.not. do_cldice) then
      ! If we are NOT prognosing ice and snow tendencies, then get them from the Pbuf
      call pbuf_get_field(pbuf, tnd_qsnow_idx,   tnd_qsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, tnd_nsnow_idx,   tnd_nsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, re_ice_idx,      re_ice,      col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      ! If we ARE prognosing tendencies, then just point to an array of NaN fields to have
      ! something for PUMAS to use in call
      tnd_qsnow => nan_array
      tnd_nsnow => nan_array
      re_ice => nan_array
   end if

   if (use_hetfrz_classnuc) then
      call pbuf_get_field(pbuf, frzimm_idx, frzimm, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzcnt_idx, frzcnt, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzdep_idx, frzdep, col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      ! Needed to satisfy gnu compiler with optional argument - set to an array of Nan fields
      frzimm => nan_array
      frzcnt => nan_array
      frzdep => nan_array
   end if

   if (qsatfac_idx > 0) then
      call pbuf_get_field(pbuf, qsatfac_idx, qsatfac, col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      allocate(qsatfac(ncol,pver),stat=ierr)
      call handle_allocate_error(ierr, 'micro_pumas_cam_tend', 'qsatfac')
      qsatfac = 1._r8
   end if

   ! initialize tendency variables
    preci  = 0._r8
    prect  = 0._r8


   !-----------------------
   ! These physics buffer fields are calculated and set in this parameterization
   ! If subcolumns is turned on, then these fields will be calculated on a subcolumn grid, otherwise they will be a normal grid

   call pbuf_get_field(pbuf, prec_str_idx,    prec_str,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_str_idx,    snow_str,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, nevapr_idx,      nevapr,      col_type=col_type)
   call pbuf_get_field(pbuf, prer_evap_idx,   prer_evap,   col_type=col_type)
   call pbuf_get_field(pbuf, prain_idx,       prain,       col_type=col_type)
   call pbuf_get_field(pbuf, dei_idx,         dei,         col_type=col_type)
   call pbuf_get_field(pbuf, mu_idx,          mu,          col_type=col_type)
   call pbuf_get_field(pbuf, lambdac_idx,     lambdac,     col_type=col_type)
   call pbuf_get_field(pbuf, des_idx,         des,         col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc,     col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw,     col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq,   col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice,   col_type=col_type)
   call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, icswp_idx,       icswp,       col_type=col_type)
   call pbuf_get_field(pbuf, rel_idx,         rel,         col_type=col_type)
   call pbuf_get_field(pbuf, rei_idx,         rei,         col_type=col_type)
   call pbuf_get_field(pbuf, sadice_idx,      sadice,      col_type=col_type)
   call pbuf_get_field(pbuf, sadsnow_idx,     sadsnow,     col_type=col_type)
   call pbuf_get_field(pbuf, wsedl_idx,       wsedl,       col_type=col_type)
   call pbuf_get_field(pbuf, qme_idx,         qme,         col_type=col_type)
   call pbuf_get_field(pbuf, bergso_idx,      bergstot,    col_type=col_type)

   ! Assign the pointer values to the non-pointer proc_rates element
   proc_rates%bergstot(:ncol,1:nlev) = bergstot(:ncol,top_lev:pver)

   if (degrau_idx > 0)   call pbuf_get_field(pbuf, degrau_idx,   degrau,   col_type=col_type)
   if (icgrauwp_idx > 0) call pbuf_get_field(pbuf, icgrauwp_idx, icgrauwp, col_type=col_type)
   if (cldfgrau_idx > 0) call pbuf_get_field(pbuf, cldfgrau_idx, cldfgrau, col_type=col_type)

   call pbuf_get_field(pbuf, cldo_idx,        cldo,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cldfsnow_idx,    cldfsnow, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_t_idx,        CC_t,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qv_idx,       CC_qv,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ql_idx,       CC_ql,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qi_idx,       CC_qi,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_nl_idx,       CC_nl,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ni_idx,       CC_ni,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qlst_idx,     CC_qlst,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)

   if (rate1_cw2pr_st_idx > 0) then
      call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st, col_type=col_type)
   end if

   if (qrain_idx > 0) call pbuf_get_field(pbuf, qrain_idx, qrout_grid_ptr)
   if (qsnow_idx > 0) call pbuf_get_field(pbuf, qsnow_idx, qsout_grid_ptr)
   if (nrain_idx > 0) call pbuf_get_field(pbuf, nrain_idx, nrout_grid_ptr)
   if (nsnow_idx > 0) call pbuf_get_field(pbuf, nsnow_idx, nsout_grid_ptr)
   if (qcsedten_idx > 0) call pbuf_get_field(pbuf, qcsedten_idx, qcsedtenout_grid_ptr)
   if (qrsedten_idx > 0) call pbuf_get_field(pbuf, qrsedten_idx, qrsedtenout_grid_ptr)
   if (qisedten_idx > 0) call pbuf_get_field(pbuf, qisedten_idx, qisedtenout_grid_ptr)
   if (qssedten_idx > 0) call pbuf_get_field(pbuf, qssedten_idx, qssedtenout_grid_ptr)
   if (vtrmc_idx > 0) call pbuf_get_field(pbuf, vtrmc_idx, vtrmcout_grid_ptr)
   if (umr_idx > 0) call pbuf_get_field(pbuf, umr_idx, umrout_grid_ptr)
   if (vtrmi_idx > 0) call pbuf_get_field(pbuf, vtrmi_idx, vtrmiout_grid_ptr)
   if (ums_idx > 0) call pbuf_get_field(pbuf, ums_idx, umsout_grid_ptr)
   if (qcsevap_idx > 0) call pbuf_get_field(pbuf, qcsevap_idx, qcsevapout_grid_ptr)
   if (qisevap_idx > 0) call pbuf_get_field(pbuf, qisevap_idx, qisevapout_grid_ptr)

   !-----------------------
   ! If subcolumns is turned on, all calculated fields which are on subcolumns
   ! need to be retrieved on the grid as well for storing averaged values

   if (use_subcol_microp) then
      call pbuf_get_field(pbuf, prec_str_idx,    prec_str_grid)
      call pbuf_get_field(pbuf, snow_str_idx,    snow_str_grid)
      call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw_grid)
      call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw_grid)
      call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed_grid)
      call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed_grid)
      call pbuf_get_field(pbuf, nevapr_idx,      nevapr_grid)
      call pbuf_get_field(pbuf, prain_idx,       prain_grid)
      call pbuf_get_field(pbuf, dei_idx,         dei_grid)
      call pbuf_get_field(pbuf, mu_idx,          mu_grid)
      call pbuf_get_field(pbuf, lambdac_idx,     lambdac_grid)
      call pbuf_get_field(pbuf, des_idx,         des_grid)
      call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc_grid)
      call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw_grid)
      call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc_grid)
      call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw_grid)
      call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq_grid)
      call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice_grid)
      call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst_grid)
      call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst_grid)
      call pbuf_get_field(pbuf, icswp_idx,       icswp_grid)
      call pbuf_get_field(pbuf, rel_idx,         rel_grid)
      call pbuf_get_field(pbuf, rei_idx,         rei_grid)
      call pbuf_get_field(pbuf, sadice_idx,      sadice_grid)
      call pbuf_get_field(pbuf, sadsnow_idx,     sadsnow_grid)
      call pbuf_get_field(pbuf, wsedl_idx,       wsedl_grid)
      call pbuf_get_field(pbuf, qme_idx,         qme_grid)
      call pbuf_get_field(pbuf, bergso_idx,      bergso_grid)
      if (degrau_idx > 0)   call pbuf_get_field(pbuf, degrau_idx,   degrau_grid)
      if (icgrauwp_idx > 0) call pbuf_get_field(pbuf, icgrauwp_idx, icgrauwp_grid)
      if (cldfgrau_idx > 0) call pbuf_get_field(pbuf, cldfgrau_idx, cldfgrau_grid)

      call pbuf_get_field(pbuf, cldo_idx,     cldo_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_t_idx,     CC_t_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qv_idx,    CC_qv_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ql_idx,    CC_ql_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qi_idx,    CC_qi_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_nl_idx,    CC_nl_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ni_idx,    CC_ni_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qlst_idx,  CC_qlst_grid,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

      if (rate1_cw2pr_st_idx > 0) then
         call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st_grid)
      end if

   else
      allocate(bergso_grid(pcols,pver), stat=ierr)
      call handle_allocate_error(ierr, 'micro_pumas_cam_tend', 'bergso_grid')
      bergso_grid(:,:) = 0._r8
   end if

   !-----------------------
   ! These are only on the grid regardless of whether subcolumns are turned on or not
   call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain_grid)
   call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow_grid)
   call pbuf_get_field(pbuf, acpr_idx,        acprecl_grid)
   call pbuf_get_field(pbuf, acgcme_idx,      acgcme_grid)
   call pbuf_get_field(pbuf, acnum_idx,       acnum_grid)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq_grid)
   call pbuf_get_field(pbuf, ast_idx,         ast_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, evprain_st_idx,  evprain_st_grid)
   call pbuf_get_field(pbuf, evpsnow_st_idx,  evpsnow_st_grid)
   call pbuf_get_field(pbuf, am_evp_st_idx,   am_evp_st_grid)

   !-------------------------------------------------------------------------------------
   ! Microphysics assumes 'liquid stratus frac = ice stratus frac
   !                      = max( liquid stratus frac, ice stratus frac )'.
   alst_mic => ast
   aist_mic => ast

   ! Output initial in-cloud LWP (before microphysics)

   iclwpi = 0._r8
   iciwpi = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         iclwpi(i) = iclwpi(i) + &
              min(state%q(i,k,ixcldliq) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
         iciwpi(i) = iciwpi(i) + &
              min(state%q(i,k,ixcldice) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
      end do
   end do

   cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

   ! Initialize local state from input.
   call physics_state_copy(state, state_loc)

   ! Because of the of limited vertical resolution, there can be a signifcant
   ! warm bias at the cold point tropopause, which can create a wet bias in the
   ! stratosphere. For the microphysics only, update the cold point temperature, with
   ! an estimate of the coldest point between the model layers.
   if (micro_mg_adjust_cpt) then
      cp_rh(:ncol, :pver)  = 0._r8
      cp_dt(:ncol)         = 0._r8
      cp_dz(:ncol)         = 0._r8

      call tropopause_find(state_loc, troplev, primary=TROP_ALG_CPP, backup=TROP_ALG_NONE, &
                           tropZ=cp_z, tropT=cp_t)

      do i = 1, ncol

         ! Update statistics and output results.
         if (troplev(i) .ne. NOTFOUND) then
            cp_dt(i) = cp_t(i) - state_loc%t(i,troplev(i))
            cp_dz(i) = cp_z(i) - state_loc%zm(i,troplev(i))

            ! NOTE: This change in temperature is just for the microphysics
            ! and should not be added to any tendencies or used to update
            ! any states
            state_loc%t(i,troplev(i)) = state_loc%t(i,troplev(i)) + cp_dt(i)
         end if
      end do

      ! Output all of the statistics related to the cold point
      ! tropopause adjustment. Th cold point information itself is
      ! output in tropopause.F90.
      call outfld("TROPF_TADJ", state_loc%t, pcols, lchnk)
      call outfld("TROPF_CDT",  cp_dt,       pcols, lchnk)
      call outfld("TROPF_CDZ",  cp_dz,       pcols, lchnk)
   end if

   ! Initialize ptend for output.
   lq = .false.
   lq(ixq) = .true.
   lq(ixcldliq) = .true.
   lq(ixcldice) = .true.
   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.
   lq(ixrain) = .true.
   lq(ixsnow) = .true.
   lq(ixnumrain) = .true.
   lq(ixnumsnow) = .true.
   if (micro_mg_version > 2) then
      lq(ixgraupel) = .true.
      lq(ixnumgraupel) = .true.
   end if

   ! the name 'cldwat' triggers special tests on cldliq
   ! and cldice in physics_update
   call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)

   if (micro_mg_version > 2) then
      state_loc_graup(:ncol,:) = state_loc%q(:ncol,:,ixgraupel)
      state_loc_numgraup(:ncol,:) = state_loc%q(:ncol,:,ixnumgraupel)
   else
      state_loc_graup(:ncol,:) = 0._r8
      state_loc_numgraup(:ncol,:) = 0._r8
   end if

   ! Zero out values above top_lev before passing into _tend for some pbuf variables that are inputs
   naai(:ncol,:top_lev-1) = 0._r8
   npccn(:ncol,:top_lev-1) = 0._r8

   ! The null value for qsatfac is 1, not zero
   qsatfac(:ncol,:top_lev-1) = 1._r8

   ! Zero out values above top_lev for all output variables
   ! Note that elements in proc_rates do not have the extra levels as they are dimensioned to be nlev instead of pver
   tlat(:ncol,:top_lev-1)=0._r8
   qvlat(:ncol,:top_lev-1)=0._r8
   qcten(:ncol,:top_lev-1)=0._r8
   qiten(:ncol,:top_lev-1)=0._r8
   ncten(:ncol,:top_lev-1)=0._r8
   niten(:ncol,:top_lev-1)=0._r8
   qrten(:ncol,:top_lev-1)=0._r8
   qsten(:ncol,:top_lev-1)=0._r8
   nrten(:ncol,:top_lev-1)=0._r8
   nsten(:ncol,:top_lev-1)=0._r8
   qgten(:ncol,:top_lev-1)=0._r8
   ngten(:ncol,:top_lev-1)=0._r8
   rel(:ncol,:top_lev-1)=0._r8
   rel_fn_dum(:ncol,:top_lev-1)=0._r8
   rei(:ncol,:top_lev-1)=0._r8
   sadice(:ncol,:top_lev-1)=0._r8
   sadsnow(:ncol,:top_lev-1)=0._r8
   prect(:ncol)=0._r8
   preci(:ncol)=0._r8
   nevapr(:ncol,:top_lev-1)=0._r8
   am_evp_st(:ncol,:top_lev-1)=0._r8
   prain(:ncol,:top_lev-1)=0._r8
   cmeice(:ncol,:top_lev-1)=0._r8
   dei(:ncol,:top_lev-1)=0._r8
   mu(:ncol,:top_lev-1)=0._r8
   lambdac(:ncol,:top_lev-1)=0._r8
   qsout(:ncol,:top_lev-1)=0._r8
   des(:ncol,:top_lev-1)=0._r8
   qgout(:ncol,:top_lev-1)=0._r8
   ngout(:ncol,:top_lev-1)=0._r8
   dgout(:ncol,:top_lev-1)=0._r8
   cflx(:ncol,:top_lev-1)=0._r8
   iflx(:ncol,:top_lev-1)=0._r8
   gflx(:ncol,:top_lev-1)=0._r8
   rflx(:ncol,:top_lev-1)=0._r8
   sflx(:ncol,:top_lev-1)=0._r8
   qrout(:ncol,:top_lev-1)=0._r8
   reff_rain_dum(:ncol,:top_lev-1)=0._r8
   reff_snow_dum(:ncol,:top_lev-1)=0._r8
   reff_grau_dum(:ncol,:top_lev-1)=0._r8
   nrout(:ncol,:top_lev-1)=0._r8
   nsout(:ncol,:top_lev-1)=0._r8
   refl(:ncol,:top_lev-1)=0._r8
   arefl(:ncol,:top_lev-1)=0._r8
   areflz(:ncol,:top_lev-1)=0._r8
   frefl(:ncol,:top_lev-1)=0._r8
   csrfl(:ncol,:top_lev-1)=0._r8
   acsrfl(:ncol,:top_lev-1)=0._r8
   fcsrfl(:ncol,:top_lev-1)=0._r8
   refl10cm(:ncol,:top_lev-1)=-9999._r8
   reflz10cm(:ncol,:top_lev-1)=0._r8
   rercld(:ncol,:top_lev-1)=0._r8
   ncai(:ncol,:top_lev-1)=0._r8
   ncal(:ncol,:top_lev-1)=0._r8
   qrout2(:ncol,:top_lev-1)=0._r8
   qsout2(:ncol,:top_lev-1)=0._r8
   nrout2(:ncol,:top_lev-1)=0._r8
   nsout2(:ncol,:top_lev-1)=0._r8
   qgout2(:ncol,:top_lev-1)=0._r8
   ngout2(:ncol,:top_lev-1)=0._r8
   dgout2(:ncol,:top_lev-1)=0._r8
   freqg(:ncol,:top_lev-1)=0._r8
   freqs(:ncol,:top_lev-1)=0._r8
   freqr(:ncol,:top_lev-1)=0._r8
   nfice(:ncol,:top_lev-1)=0._r8
   qcrat(:ncol,:top_lev-1)=0._r8
   tnd_qsnow(:ncol,:top_lev-1)=0._r8
   tnd_nsnow(:ncol,:top_lev-1)=0._r8
   re_ice(:ncol,:top_lev-1)=0._r8
   prer_evap(:ncol,:top_lev-1)=0._r8
   frzimm(:ncol,:top_lev-1)=0._r8
   frzcnt(:ncol,:top_lev-1)=0._r8
   frzdep(:ncol,:top_lev-1)=0._r8

   do it = 1, num_steps

     call micro_pumas_tend( &
              ncol,         nlev,           dtime/num_steps,&
              state_loc%t(:ncol,top_lev:),              state_loc%q(:ncol,top_lev:,ixq),            &
              state_loc%q(:ncol,top_lev:,ixcldliq),     state_loc%q(:ncol,top_lev:,ixcldice),          &
              state_loc%q(:ncol,top_lev:,ixnumliq),     state_loc%q(:ncol,top_lev:,ixnumice),       &
              state_loc%q(:ncol,top_lev:,ixrain),       state_loc%q(:ncol,top_lev:,ixsnow),         &
              state_loc%q(:ncol,top_lev:,ixnumrain),    state_loc%q(:ncol,top_lev:,ixnumsnow),      &
              state_loc_graup(:ncol,top_lev:),    state_loc_numgraup(:ncol,top_lev:),     &
              relvar(:ncol,top_lev:),         accre_enhan(:ncol,top_lev:),     &
              state_loc%pmid(:ncol,top_lev:),                state_loc%pdel(:ncol,top_lev:),  state_loc%pint(:ncol,top_lev:), &
              ast(:ncol,top_lev:), alst_mic(:ncol,top_lev:), aist_mic(:ncol,top_lev:), qsatfac(:ncol,top_lev:), &
              rate1cld(:ncol,top_lev:),                         &
              naai(:ncol,top_lev:),            npccn(:ncol,top_lev:),           &
              rndst(:ncol,top_lev:,:),    nacon(:ncol,top_lev:,:),           &
              tlat(:ncol,top_lev:),            qvlat(:ncol,top_lev:),           &
              qcten(:ncol,top_lev:),          qiten(:ncol,top_lev:),          &
              ncten(:ncol,top_lev:),          niten(:ncol,top_lev:),          &
              qrten(:ncol,top_lev:),          qsten(:ncol,top_lev:),          &
              nrten(:ncol,top_lev:),          nsten(:ncol,top_lev:),          &
              qgten(:ncol,top_lev:),          ngten(:ncol,top_lev:),          &
              rel(:ncol,top_lev:),     rel_fn_dum(:ncol,top_lev:),     rei(:ncol,top_lev:),     &
              sadice(:ncol,top_lev:),          sadsnow(:ncol,top_lev:),         &
              prect(:ncol),           preci(:ncol),           &
              nevapr(:ncol,top_lev:),          am_evp_st(:ncol,top_lev:),       &
              prain(:ncol,top_lev:),                   &
              cmeice(:ncol,top_lev:),          dei(:ncol,top_lev:),             &
              mu(:ncol,top_lev:),              lambdac(:ncol,top_lev:),         &
              qsout(:ncol,top_lev:),           des(:ncol,top_lev:),             &
              qgout(:ncol,top_lev:),   ngout(:ncol,top_lev:),   dgout(:ncol,top_lev:),   &
              cflx(:ncol,top_lev:),    iflx(:ncol,top_lev:),                    &
              gflx(:ncol,top_lev:),                                    &
              rflx(:ncol,top_lev:),    sflx(:ncol,top_lev:),    qrout(:ncol,top_lev:),   &
              reff_rain_dum(:ncol,top_lev:),          reff_snow_dum(:ncol,top_lev:),   reff_grau_dum(:ncol,top_lev:),       &
              nrout(:ncol,top_lev:),           nsout(:ncol,top_lev:),           &
              refl(:ncol,top_lev:),    arefl(:ncol,top_lev:),   areflz(:ncol,top_lev:),  &
              frefl(:ncol,top_lev:),   csrfl(:ncol,top_lev:),   acsrfl(:ncol,top_lev:),  &
              fcsrfl(:ncol,top_lev:),   &
              refl10cm(:ncol,top_lev:), reflz10cm(:ncol,top_lev:),    rercld(:ncol,top_lev:),          &
              ncai(:ncol,top_lev:),            ncal(:ncol,top_lev:),            &
              qrout2(:ncol,top_lev:),          qsout2(:ncol,top_lev:),          &
              nrout2(:ncol,top_lev:),          nsout2(:ncol,top_lev:),          &
              drout_dum(:ncol,top_lev:),              dsout2_dum(:ncol,top_lev:),             &
              qgout2(:ncol,top_lev:), ngout2(:ncol,top_lev:), dgout2(:ncol,top_lev:), freqg(:ncol,top_lev:),   &
              freqs(:ncol,top_lev:),           freqr(:ncol,top_lev:),           &
              nfice(:ncol,top_lev:),           qcrat(:ncol,top_lev:),           &
              proc_rates,                                                       &
              errstring, &
              tnd_qsnow(:ncol,top_lev:),tnd_nsnow(:ncol,top_lev:),re_ice(:ncol,top_lev:),&
              prer_evap(:ncol,top_lev:),                                     &
              frzimm(:ncol,top_lev:),  frzcnt(:ncol,top_lev:),  frzdep(:ncol,top_lev:)   )

      call handle_errmsg(errstring, subname="micro_pumas_cam_tend")

      call physics_ptend_init(ptend_loc, psetcols, "micro_pumas", &
                              ls=.true., lq=lq)

      ! Set local tendency.
      ptend_loc%s(:ncol,top_lev:) = tlat(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixq) = qvlat(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixcldliq) = qcten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixcldice) = qiten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumliq) = ncten(:ncol,top_lev:)

      if (do_cldice) then
         ptend_loc%q(:ncol,top_lev:,ixnumice) = niten(:ncol,top_lev:)
      else
         ! In this case, the tendency should be all 0.
         if (any(niten(:ncol,:) /= 0._r8)) then
              call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
              " but micro_pumas_tend has ice number tendencies.")
         end if
         ptend_loc%q(:ncol,:,ixnumice) = 0._r8
      end if

      ptend_loc%q(:ncol,top_lev:,ixrain) = qrten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixsnow) = qsten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumrain) = nrten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumsnow) = nsten(:ncol,top_lev:)

      if (micro_mg_version > 2) then
         ptend_loc%q(:ncol,top_lev:,ixgraupel) = qgten(:ncol,top_lev:)
         ptend_loc%q(:ncol,top_lev:,ixnumgraupel) = ngten(:ncol,top_lev:)
      end if

      ! Sum into overall ptend
      call physics_ptend_sum(ptend_loc, ptend, ncol)

      ! Update local state
      call physics_update(state_loc, ptend_loc, dtime/num_steps)

      if (trim(micro_mg_warm_rain) == 'tau') then
         proc_rates%amk_c(:ncol,:,:) = proc_rates%amk_c(:ncol,:,:)/num_steps
         proc_rates%ank_c(:ncol,:,:) = proc_rates%ank_c(:ncol,:,:)/num_steps
         proc_rates%amk_r(:ncol,:,:) = proc_rates%amk_r(:ncol,:,:)/num_steps
         proc_rates%ank_r(:ncol,:,:) = proc_rates%ank_r(:ncol,:,:)/num_steps
         proc_rates%amk(:ncol,:,:) = proc_rates%amk(:ncol,:,:)/num_steps
         proc_rates%ank(:ncol,:,:) = proc_rates%ank(:ncol,:,:)/num_steps
         proc_rates%amk_out(:ncol,:,:) = proc_rates%amk_out(:ncol,:,:)/num_steps
      end if

   end do

   ! Divide ptend by substeps.
   call physics_ptend_scale(ptend, 1._r8/num_steps, ncol)

   ! Check to make sure that the microphysics code is respecting the flags that control
   ! whether MG should be prognosing cloud ice and cloud liquid or not.
   if (.not. do_cldice) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldice) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_pumas_tend has ice mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumice) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_pumas_tend has ice number tendencies.")
   end if
   if (.not. do_cldliq) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldliq) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_pumas_tend has liquid mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumliq) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_pumas_tend has liquid number tendencies.")
   end if

   mnuccdohet = 0._r8
   do k=top_lev,pver
      do i=1,ncol
         if (naai(i,k) > 0._r8) then
            mnuccdohet(i,k) = proc_rates%mnuccdtot(i,k-top_lev+1) - (naai_hom(i,k)/naai(i,k))*proc_rates%mnuccdtot(i,k-top_lev+1)
         end if
      end do
   end do

   mgflxprc(:ncol,top_lev:pverp) = rflx(:ncol,top_lev:pverp) + sflx(:ncol,top_lev:pverp)
   mgflxsnw(:ncol,top_lev:pverp) = sflx(:ncol,top_lev:pverp)

   !add condensate fluxes for MG2 (ice and snow already added for MG1)
   if (micro_mg_version >= 2) then
      mgflxprc(:ncol,top_lev:pverp) = mgflxprc(:ncol,top_lev:pverp)+ iflx(:ncol,top_lev:pverp) + cflx(:ncol,top_lev:pverp)
      mgflxsnw(:ncol,top_lev:pverp) = mgflxsnw(:ncol,top_lev:pverp) + iflx(:ncol,top_lev:pverp)
   end if

   !add graupel fluxes for MG3 to snow flux
   if (micro_mg_version >= 3) then
      mgflxprc(:ncol,top_lev:pverp) = mgflxprc(:ncol,top_lev:pverp)+gflx(:ncol,top_lev:pverp)
      mgflxsnw(:ncol,top_lev:pverp) = mgflxsnw(:ncol,top_lev:pverp)+gflx(:ncol,top_lev:pverp)
   end if

   mgmrprc(:ncol,top_lev:pver) = qrout(:ncol,top_lev:pver) + qsout(:ncol,top_lev:pver)
   mgmrsnw(:ncol,top_lev:pver) = qsout(:ncol,top_lev:pver)

   !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
   !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
   cvreffliq(:ncol,top_lev:pver) = 9.0_r8
   cvreffice(:ncol,top_lev:pver) = 37.0_r8

   ! Reassign rate1 if modal aerosols
   if (rate1_cw2pr_st_idx > 0) then
      rate1ord_cw2pr_st(:ncol,top_lev:pver) = rate1cld(:ncol,top_lev:pver)
   end if

   ! Sedimentation velocity for liquid stratus cloud droplet
   wsedl(:ncol,top_lev:pver) = proc_rates%vtrmc(:ncol,1:nlev)

   ! Microphysical tendencies for use in the macrophysics at the next time step
   CC_T(:ncol,top_lev:pver)    = tlat(:ncol,top_lev:pver)/cpair
   CC_qv(:ncol,top_lev:pver)   = qvlat(:ncol,top_lev:pver)
   CC_ql(:ncol,top_lev:pver)   = qcten(:ncol,top_lev:pver)
   CC_qi(:ncol,top_lev:pver)   = qiten(:ncol,top_lev:pver)
   CC_nl(:ncol,top_lev:pver)   = ncten(:ncol,top_lev:pver)
   CC_ni(:ncol,top_lev:pver)   = niten(:ncol,top_lev:pver)
   CC_qlst(:ncol,top_lev:pver) = qcten(:ncol,top_lev:pver)/max(0.01_r8,alst_mic(:ncol,top_lev:pver))

   ! Net micro_pumas_cam condensation rate
   qme(:ncol,:top_lev-1) = 0._r8
   qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + proc_rates%cmeitot(:ncol,1:nlev)

   ! For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
   ! Other precip output variables are set to 0
   ! Do not subscript by ncol here, because in physpkg we divide the whole
   ! array and need to avoid an FPE due to uninitialized data.
   prec_pcw = prect
   snow_pcw = preci
   prec_sed = 0._r8
   snow_sed = 0._r8
   prec_str = prec_pcw + prec_sed
   snow_str = snow_pcw + snow_sed

   icecldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
   liqcldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)

   ! ------------------------------------------------------------ !
   ! Compute in cloud ice and liquid mixing ratios                !
   ! Note that 'iclwp, iciwp' are used for radiation computation. !
   ! ------------------------------------------------------------ !

   icinc = 0._r8
   icwnc = 0._r8
   iciwpst = 0._r8
   iclwpst = 0._r8
   icswp = 0._r8
   cldfsnow = 0._r8
   if (micro_mg_version > 2) then
      icgrauwp = 0._r8
      cldfgrau = 0._r8
   end if

   do k = top_lev, pver
      do i = 1, ncol
         ! Limits for in-cloud mixing ratios consistent with MG microphysics
         ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
         icimrst(i,k)   = min( state_loc%q(i,k,ixcldice) / max(mincld,icecldf(i,k)),0.005_r8 )
         icwmrst(i,k)   = min( state_loc%q(i,k,ixcldliq) / max(mincld,liqcldf(i,k)),0.005_r8 )
         icinc(i,k)     = state_loc%q(i,k,ixnumice) / max(mincld,icecldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         icwnc(i,k)     = state_loc%q(i,k,ixnumliq) / max(mincld,liqcldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         ! Calculate micro_pumas_cam cloud water paths in each layer
         ! Note: uses stratiform cloud fraction!
         iciwpst(i,k)   = min(state_loc%q(i,k,ixcldice)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
         iclwpst(i,k)   = min(state_loc%q(i,k,ixcldliq)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit

         ! ------------------------------ !
         ! Adjust cloud fraction for snow !
         ! ------------------------------ !
         cldfsnow(i,k) = cld(i,k)
         ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
         if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. &
            ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
            ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
            cldfsnow(i,k) = 0._r8
         end if
         ! If no cloud and snow, then set to 0.25
         if( ( cldfsnow(i,k) .le. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then
            cldfsnow(i,k) = 0.25_r8
         end if
         ! Calculate in-cloud snow water path
         icswp(i,k) = qsout(i,k) / max( mincld, cldfsnow(i,k) ) * state_loc%pdel(i,k) / gravit

         ! --------------------------------- !
         ! Adjust cloud fraction for graupel !
         ! --------------------------------- !
       if (micro_mg_version > 2) then
          cldfgrau(i,k) = cld(i,k)
         ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
          if( ( cldfgrau(i,k) .gt. 1.e-4_r8 ) .and. &
              ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
              ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
              cldfgrau(i,k) = 0._r8
           end if
         ! If no cloud and graupel, then set to 0.25
           if( ( cldfgrau(i,k) .le. 1.e-4_r8 ) .and. ( qgout(i,k) .gt. 1.e-9_r8 ) ) then
              cldfgrau(i,k) = 0.25_r8
           end if

         ! Calculate in-cloud snow water path
           icgrauwp(i,k) = qgout(i,k) / max( 1.e-2_r8, cldfgrau(i,k) ) * state_loc%pdel(i,k) / gravit
        end if

      end do
   end do

   ! Calculate cloud fraction for prognostic precip sizes.
   ! Cloud fraction for purposes of precipitation is maximum cloud
   ! fraction out of all the layers that the precipitation may be
   ! falling down from.
   cldmax(:ncol,top_lev:) = max(mincld, ast(:ncol,top_lev:))
   do k = top_lev+1, pver
      where (state_loc%q(:ncol,k-1,ixrain) >= qsmall .or. &
           state_loc%q(:ncol,k-1,ixsnow) >= qsmall)
         cldmax(:ncol,k) = max(cldmax(:ncol,k-1), cldmax(:ncol,k))
      end where
   end do

   !Copy pbuf field from proc_rates back to pbuf pointer
   bergstot(:ncol,top_lev:) = proc_rates%bergstot(:ncol,1:nlev)
   bergstot(:ncol,1:top_lev-1) = 0._r8

   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !
   ! All code from here to the end is on grid columns only  !
   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !

   ! Average the fields which are needed later in this paramterization to be on the grid
   if (use_subcol_microp) then
      call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
      call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
      call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
      call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
      call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
      call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
      call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
      call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
      call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)
      call subcol_field_avg(proc_rates%evapsnow,  ngrdcol, lchnk, evpsnow_st_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%bergstot,    ngrdcol, lchnk, bergso_grid(:,top_lev:))

      call subcol_field_avg(am_evp_st, ngrdcol, lchnk, am_evp_st_grid)

      ! Average fields which are not in pbuf
      call subcol_field_avg(qrout,     ngrdcol, lchnk, qrout_grid)
      call subcol_field_avg(qsout,     ngrdcol, lchnk, qsout_grid)
      call subcol_field_avg(nsout,     ngrdcol, lchnk, nsout_grid)
      call subcol_field_avg(nrout,     ngrdcol, lchnk, nrout_grid)
      call subcol_field_avg(cld,       ngrdcol, lchnk, cld_grid)
      call subcol_field_avg(proc_rates%qcrestot,    ngrdcol, lchnk, qcreso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%melttot,     ngrdcol, lchnk, melto_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%mnuccctot,   ngrdcol, lchnk, mnuccco_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%mnuccttot,   ngrdcol, lchnk, mnuccto_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%bergtot,     ngrdcol, lchnk, bergo_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%homotot,     ngrdcol, lchnk, homoo_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%msacwitot,   ngrdcol, lchnk, msacwio_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%psacwstot,   ngrdcol, lchnk, psacwso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%cmeitot,   ngrdcol, lchnk, cmeiout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qirestot,    ngrdcol, lchnk, qireso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%prcitot,     ngrdcol, lchnk, prcio_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%praitot,     ngrdcol, lchnk, praio_grid(:,top_lev:))
      call subcol_field_avg(icwmrst,   ngrdcol, lchnk, icwmrst_grid)
      call subcol_field_avg(icimrst,   ngrdcol, lchnk, icimrst_grid)
      call subcol_field_avg(liqcldf,   ngrdcol, lchnk, liqcldf_grid)
      call subcol_field_avg(icecldf,   ngrdcol, lchnk, icecldf_grid)
      call subcol_field_avg(icwnc,     ngrdcol, lchnk, icwnc_grid)
      call subcol_field_avg(icinc,     ngrdcol, lchnk, icinc_grid)
      call subcol_field_avg(state_loc%pdel,            ngrdcol, lchnk, pdel_grid)
      call subcol_field_avg(proc_rates%pratot,      ngrdcol, lchnk, prao_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%prctot,      ngrdcol, lchnk, prco_grid(:,top_lev:))

      call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, nc_grid(:,top_lev:))
      call subcol_field_avg(state_loc%q(:,:,ixnumice), ngrdcol, lchnk, ni_grid(:,top_lev:))

      call subcol_field_avg(proc_rates%qcsedten,  ngrdcol, lchnk, qcsedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qisedten,  ngrdcol, lchnk, qisedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%vtrmc,     ngrdcol, lchnk, vtrmcout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%vtrmi,     ngrdcol, lchnk, vtrmiout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qcsevap,  ngrdcol, lchnk, qcsevapout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qisevap,  ngrdcol, lchnk, qisevapout_grid(:,top_lev:))

      call subcol_field_avg(cldmax,    ngrdcol, lchnk, cldmax_grid)

      call subcol_field_avg(state_loc%q(:,:,ixrain),    ngrdcol, lchnk, qr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumrain), ngrdcol, lchnk, nr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixsnow),    ngrdcol, lchnk, qs_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumsnow), ngrdcol, lchnk, ns_grid)
      call subcol_field_avg(proc_rates%qrsedten,  ngrdcol, lchnk, qrsedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qssedten,  ngrdcol, lchnk, qssedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%umr,       ngrdcol, lchnk, umrout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%ums,       ngrdcol, lchnk, umsout_grid(:,top_lev:))

      if (micro_mg_version > 2) then
            call subcol_field_avg(state_loc%q(:,:,ixgraupel),    ngrdcol, lchnk, qg_grid)
            call subcol_field_avg(state_loc%q(:,:,ixnumgraupel), ngrdcol, lchnk, ng_grid)
            call subcol_field_avg(proc_rates%psacrtot,       ngrdcol, lchnk, psacro_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pracgtot,       ngrdcol, lchnk, pracgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%psacwgtot,      ngrdcol, lchnk, psacwgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pgsacwtot,      ngrdcol, lchnk, pgsacwo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pgracstot,      ngrdcol, lchnk, pgracso_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%prdgtot,        ngrdcol, lchnk, prdgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%qmultgtot,      ngrdcol, lchnk, qmultgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%qmultrgtot,     ngrdcol, lchnk, qmultrgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%npracgtot,      ngrdcol, lchnk, npracgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nscngtot,       ngrdcol, lchnk, nscngo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%ngracstot,      ngrdcol, lchnk, ngracso_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nmultgtot,      ngrdcol, lchnk, nmultgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nmultrgtot,     ngrdcol, lchnk, nmultrgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%npsacwgtot,     ngrdcol, lchnk, npsacwgo_grid(:,top_lev:))
      end if

   else
      qcreso_grid(:ncol,:top_lev-1)     = 0._r8
      melto_grid(:ncol,:top_lev-1)      = 0._r8
      mnuccco_grid(:ncol,:top_lev-1)    = 0._r8
      mnuccto_grid(:ncol,:top_lev-1)    = 0._r8
      bergo_grid(:ncol,:top_lev-1)      = 0._r8
      homoo_grid(:ncol,:top_lev-1)      = 0._r8
      msacwio_grid(:ncol,:top_lev-1)    = 0._r8
      psacwso_grid(:ncol,:top_lev-1)    = 0._r8
      cmeiout_grid(:ncol,:top_lev-1)    = 0._r8
      qireso_grid(:ncol,:top_lev-1)     = 0._r8
      prcio_grid(:ncol,:top_lev-1)      = 0._r8
      praio_grid(:ncol,:top_lev-1)      = 0._r8
      prao_grid(:ncol,:top_lev-1)       = 0._r8
      prco_grid(:ncol,:top_lev-1)       = 0._r8
      qcsedtenout_grid(:ncol,:top_lev-1) = 0._r8
      qisedtenout_grid(:ncol,:top_lev-1) = 0._r8
      vtrmcout_grid(:ncol,:top_lev-1)    = 0._r8
      vtrmiout_grid(:ncol,:top_lev-1)    = 0._r8
      qcsevapout_grid(:ncol,:top_lev-1) = 0._r8
      qisevapout_grid(:ncol,:top_lev-1) = 0._r8
      qrsedtenout_grid(:ncol,:top_lev-1) = 0._r8
      qssedtenout_grid(:ncol,:top_lev-1) = 0._r8
      umrout_grid(:ncol,:top_lev-1) = 0._r8
      umsout_grid(:ncol,:top_lev-1) = 0._r8
      psacro_grid(:ncol,:top_lev-1) = 0._r8
      pracgo_grid(:ncol,:top_lev-1) = 0._r8
      psacwgo_grid(:ncol,:top_lev-1) = 0._r8
      pgsacwo_grid(:ncol,:top_lev-1) = 0._r8
      pgracso_grid(:ncol,:top_lev-1) = 0._r8
      prdgo_grid(:ncol,:top_lev-1) = 0._r8
      qmultgo_grid(:ncol,:top_lev-1) = 0._r8
      qmultrgo_grid(:ncol,:top_lev-1) = 0._r8
      npracgo_grid(:ncol,:top_lev-1) = 0._r8
      nscngo_grid(:ncol,:top_lev-1) = 0._r8
      ngracso_grid(:ncol,:top_lev-1) = 0._r8
      nmultgo_grid(:ncol,:top_lev-1) = 0._r8
      nmultrgo_grid(:ncol,:top_lev-1) = 0._r8
      npsacwgo_grid(:ncol,:top_lev-1) = 0._r8
      bergso_grid(:ncol,:top_lev-1) = 0._r8

      ! These pbuf fields need to be assigned.  There is no corresponding subcol_field_avg
      ! as they are reset before being used, so it would be a needless calculation
      lambdac_grid    => lambdac
      mu_grid         => mu
      rel_grid        => rel
      rei_grid        => rei
      sadice_grid     => sadice
      sadsnow_grid    => sadsnow
      dei_grid        => dei
      des_grid        => des
      degrau_grid     => degrau

      ! fields already on grids, so just assign
      prec_str_grid   => prec_str
      iclwpst_grid    => iclwpst
      cvreffliq_grid  => cvreffliq
      cvreffice_grid  => cvreffice
      mgflxprc_grid   => mgflxprc
      mgflxsnw_grid   => mgflxsnw
      qme_grid        => qme
      nevapr_grid     => nevapr
      prain_grid      => prain

      bergso_grid(:ncol,top_lev:)    =  proc_rates%bergstot
      am_evp_st_grid  = am_evp_st

      evpsnow_st_grid(:ncol,top_lev:) = proc_rates%evapsnow
      qrout_grid      = qrout
      qsout_grid      = qsout
      nsout_grid      = nsout
      nrout_grid      = nrout
      cld_grid        = cld
      qcreso_grid(:ncol,top_lev:)     = proc_rates%qcrestot
      melto_grid(:ncol,top_lev:)      = proc_rates%melttot
      mnuccco_grid(:ncol,top_lev:)    = proc_rates%mnuccctot
      mnuccto_grid(:ncol,top_lev:)    = proc_rates%mnuccttot
      bergo_grid(:ncol,top_lev:)      = proc_rates%bergtot
      homoo_grid(:ncol,top_lev:)      = proc_rates%homotot
      msacwio_grid(:ncol,top_lev:)    = proc_rates%msacwitot
      psacwso_grid(:ncol,top_lev:)    = proc_rates%psacwstot
      cmeiout_grid(:ncol,top_lev:)    = proc_rates%cmeitot
      qireso_grid(:ncol,top_lev:)     = proc_rates%qirestot
      prcio_grid(:ncol,top_lev:)      = proc_rates%prcitot
      praio_grid(:ncol,top_lev:)      = proc_rates%praitot
      icwmrst_grid    = icwmrst
      icimrst_grid    = icimrst
      liqcldf_grid    = liqcldf
      icecldf_grid    = icecldf
      icwnc_grid      = icwnc
      icinc_grid      = icinc
      pdel_grid       = state_loc%pdel
      prao_grid(:ncol,top_lev:)       = proc_rates%pratot
      prco_grid(:ncol,top_lev:)       = proc_rates%prctot

      nc_grid = state_loc%q(:,:,ixnumliq)
      ni_grid = state_loc%q(:,:,ixnumice)

      qcsedtenout_grid(:ncol,top_lev:) = proc_rates%qcsedten
      qisedtenout_grid(:ncol,top_lev:) = proc_rates%qisedten
      vtrmcout_grid(:ncol,top_lev:)    = proc_rates%vtrmc
      vtrmiout_grid(:ncol,top_lev:)    = proc_rates%vtrmi
      qcsevapout_grid(:ncol,top_lev:) = proc_rates%qcsevap
      qisevapout_grid(:ncol,top_lev:) = proc_rates%qisevap

      cldmax_grid = cldmax

      qr_grid = state_loc%q(:,:,ixrain)
      nr_grid = state_loc%q(:,:,ixnumrain)
      qs_grid = state_loc%q(:,:,ixsnow)
      ns_grid = state_loc%q(:,:,ixnumsnow)
      qrsedtenout_grid(:ncol,top_lev:) = proc_rates%qrsedten
      qssedtenout_grid(:ncol,top_lev:) = proc_rates%qssedten
      umrout_grid(:ncol,top_lev:) = proc_rates%umr
      umsout_grid(:ncol,top_lev:) = proc_rates%ums

! Zero out terms for budgets if not mg3....
      psacwgo_grid = 0._r8
      pgsacwo_grid = 0._r8
      qmultgo_grid = 0._r8

      if (micro_mg_version > 2) then
            qg_grid = state_loc%q(:,:,ixgraupel)
            ng_grid = state_loc%q(:,:,ixnumgraupel)
            psacro_grid(:ncol,top_lev:) =     proc_rates%psacrtot
            pracgo_grid(:ncol,top_lev:) =     proc_rates%pracgtot
            psacwgo_grid(:ncol,top_lev:) =    proc_rates%psacwgtot
            pgsacwo_grid(:ncol,top_lev:) =    proc_rates%pgsacwtot
            pgracso_grid(:ncol,top_lev:) =    proc_rates%pgracstot
            prdgo_grid(:ncol,top_lev:) =      proc_rates%prdgtot
            qmultgo_grid(:ncol,top_lev:) =    proc_rates%qmultgtot
            qmultrgo_grid(:ncol,top_lev:) =   proc_rates%qmultrgtot
            npracgo_grid(:ncol,top_lev:) =   proc_rates%npracgtot
            nscngo_grid(:ncol,top_lev:) =   proc_rates%nscngtot
            ngracso_grid(:ncol,top_lev:) =   proc_rates%ngracstot
            nmultgo_grid(:ncol,top_lev:) =   proc_rates%nmultgtot
            nmultrgo_grid(:ncol,top_lev:) =   proc_rates%nmultrgtot
            npsacwgo_grid(:ncol,top_lev:) =   proc_rates%npsacwgtot
      end if


   end if

   ! If on subcolumns, average the rest of the pbuf fields which were modified on subcolumns but are not used further in
   ! this parameterization  (no need to assign in the non-subcolumn case -- the else step)
   if (use_subcol_microp) then
      call subcol_field_avg(snow_str,    ngrdcol, lchnk, snow_str_grid)
      call subcol_field_avg(prec_pcw,    ngrdcol, lchnk, prec_pcw_grid)
      call subcol_field_avg(snow_pcw,    ngrdcol, lchnk, snow_pcw_grid)
      call subcol_field_avg(prec_sed,    ngrdcol, lchnk, prec_sed_grid)
      call subcol_field_avg(snow_sed,    ngrdcol, lchnk, snow_sed_grid)
      call subcol_field_avg(cldo,        ngrdcol, lchnk, cldo_grid)
      call subcol_field_avg(mgmrprc,     ngrdcol, lchnk, mgmrprc_grid)
      call subcol_field_avg(mgmrsnw,     ngrdcol, lchnk, mgmrsnw_grid)
      call subcol_field_avg(wsedl,       ngrdcol, lchnk, wsedl_grid)
      call subcol_field_avg(cc_t,        ngrdcol, lchnk, cc_t_grid)
      call subcol_field_avg(cc_qv,       ngrdcol, lchnk, cc_qv_grid)
      call subcol_field_avg(cc_ql,       ngrdcol, lchnk, cc_ql_grid)
      call subcol_field_avg(cc_qi,       ngrdcol, lchnk, cc_qi_grid)
      call subcol_field_avg(cc_nl,       ngrdcol, lchnk, cc_nl_grid)
      call subcol_field_avg(cc_ni,       ngrdcol, lchnk, cc_ni_grid)
      call subcol_field_avg(cc_qlst,     ngrdcol, lchnk, cc_qlst_grid)
      call subcol_field_avg(iciwpst,     ngrdcol, lchnk, iciwpst_grid)
      call subcol_field_avg(icswp,       ngrdcol, lchnk, icswp_grid)
      call subcol_field_avg(cldfsnow,    ngrdcol, lchnk, cldfsnow_grid)

      if (micro_mg_version > 2) then
         call subcol_field_avg(icgrauwp,    ngrdcol, lchnk, icgrauwp_grid)
         call subcol_field_avg(cldfgrau,    ngrdcol, lchnk, cldfsnow_grid)
      end if

      if (rate1_cw2pr_st_idx > 0) then
         call subcol_field_avg(rate1ord_cw2pr_st,    ngrdcol, lchnk, rate1ord_cw2pr_st_grid)
      end if

   end if

   ! ------------------------------------- !
   ! Size distribution calculation         !
   ! ------------------------------------- !

   ! Calculate rho (on subcolumns if turned on) for size distribution
   ! parameter calculations and average it if needed
   !
   ! State instead of state_loc to preserve answers for MG1 (and in any
   ! case, it is unlikely to make much difference).
   rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
        (rair*state%t(:ncol,top_lev:))
   if (use_subcol_microp) then
      call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
   else
      rho_grid = rho
   end if

   ! Effective radius for cloud liquid, fixed number.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_fn_grid = 10._r8

   ncic_grid = 1.e8_r8

   do k = top_lev, pver
      !$acc data copyin  (mg_liq_props,icwmrst_grid(:ngrdcol,k),rho_grid(:ngrdcol,k)) &
      !$acc      copy    (ncic_grid(:ngrdcol,k)) &
      !$acc      copyout (mu_grid(:ngrdcol,k),lambdac_grid(:ngrdcol,k))
      call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,k), &
                               ncic_grid(:ngrdcol,k), rho_grid(:ngrdcol,k), &
                               mu_grid(:ngrdcol,k), lambdac_grid(:ngrdcol,k), ngrdcol)
      !$acc end data
   end do

   where (icwmrst_grid(:ngrdcol,top_lev:) > qsmall)
      rel_fn_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8)/ &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   end where

   ! Effective radius for cloud liquid, and size parameters
   ! mu_grid and lambdac_grid.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_grid = 10._r8

   ! Calculate ncic on the grid
   ncic_grid(:ngrdcol,top_lev:) = nc_grid(:ngrdcol,top_lev:) / &
        max(mincld,liqcldf_grid(:ngrdcol,top_lev:))

   do k = top_lev, pver
      !$acc data copyin  (mg_liq_props,icwmrst_grid(:ngrdcol,k), rho_grid(:ngrdcol,k)) &
      !$acc      copy    (ncic_grid(:ngrdcol,k)) &
      !$acc      copyout (mu_grid(:ngrdcol,k),lambdac_grid(:ngrdcol,k))
      call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,k), &
           ncic_grid(:ngrdcol,k), rho_grid(:ngrdcol,k), &
           mu_grid(:ngrdcol,k), lambdac_grid(:ngrdcol,k), ngrdcol)
      !$acc end data
   end do

   where (icwmrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rel_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8) / &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   elsewhere
      ! Deal with the fact that size_dist_param_liq sets mu_grid to -100
      ! wherever there is no cloud.
      mu_grid(:ngrdcol,top_lev:) = 0._r8
   end where

   ! Rain/Snow effective diameter.
   drout2_grid = 0._r8
   reff_rain_grid = 0._r8
   des_grid = 0._r8
   dsout2_grid = 0._r8
   reff_snow_grid = 0._r8
   reff_grau_grid = 0._r8

   ! Prognostic precipitation

   where (qr_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
      drout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
           qr_grid(:ngrdcol,top_lev:), &
           nr_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
           rho_grid(:ngrdcol,top_lev:), rhow)

      reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
           shapeparam * micron2meter
   end where

   where (qs_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
      dsout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
           qs_grid(:ngrdcol,top_lev:), &
           ns_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
           rho_grid(:ngrdcol,top_lev:), rhosn)

      des_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) *&
           3._r8 * rhosn/rhows

      reff_snow_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) * &
           shapeparam * micron2meter
   end where


! Graupel/Hail size distribution Placeholder
   if (micro_mg_version > 2) then
      degrau_grid = 0._r8
      where (qg_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         dgout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qg_grid(:ngrdcol,top_lev:), &
              ng_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhog)

         reff_grau_grid(:ngrdcol,top_lev:) = dgout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
         degrau_grid(:ngrdcol,top_lev:) = dgout2_grid(:ngrdcol,top_lev:) *&
              3._r8 * rhog/rhows
      end where
   end if

   ! Effective radius and diameter for cloud ice.
   rei_grid = 25._r8

   niic_grid(:ngrdcol,top_lev:) = ni_grid(:ngrdcol,top_lev:) / &
        max(mincld,icecldf_grid(:ngrdcol,top_lev:))

   do k = top_lev, pver
      !$acc data copyin  (mg_ice_props, icimrst_grid(:ngrdcol,k)) &
      !$acc      copy    (niic_grid(:ngrdcol,k)) &
      !$acc      copyout (rei_grid(:ngrdcol,k))
      call size_dist_param_basic(mg_ice_props,icimrst_grid(:ngrdcol,k), &
                                 niic_grid(:ngrdcol,k),rei_grid(:ngrdcol,k),ngrdcol)
      !$acc end data
   end do

   where (icimrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rei_grid(:ngrdcol,top_lev:) = 1.5_r8/rei_grid(:ngrdcol,top_lev:) &
           * 1.e6_r8
   elsewhere
      rei_grid(:ngrdcol,top_lev:) = 25._r8
   end where

   dei_grid = rei_grid * rhoi/rhows * 2._r8

   ! Limiters for low cloud fraction.
   do k = top_lev, pver
      do i = 1, ngrdcol
         ! Convert snow effective diameter to microns
         des_grid(i,k) = des_grid(i,k) * 1.e6_r8
         if ( ast_grid(i,k) < 1.e-4_r8 ) then
            mu_grid(i,k) = mucon
            lambdac_grid(i,k) = (mucon + 1._r8)/dcon
            dei_grid(i,k) = deicon
         end if
      end do
   end do

   mgreffrain_grid(:ngrdcol,top_lev:pver) = reff_rain_grid(:ngrdcol,top_lev:pver)
   mgreffsnow_grid(:ngrdcol,top_lev:pver) = reff_snow_grid(:ngrdcol,top_lev:pver)

   ! ------------------------------------- !
   ! Precipitation efficiency Calculation  !
   ! ------------------------------------- !

   !-----------------------------------------------------------------------
   ! Liquid water path

   ! Compute liquid water paths, and column condensation
   tgliqwp_grid(:ngrdcol) = 0._r8
   tgcmeliq_grid(:ngrdcol) = 0._r8
   do k = top_lev, pver
      do i = 1, ngrdcol
         tgliqwp_grid(i)  = tgliqwp_grid(i) + iclwpst_grid(i,k)*cld_grid(i,k)

         if (cmeliq_grid(i,k) > 1.e-12_r8) then
            !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
            tgcmeliq_grid(i) = tgcmeliq_grid(i) + cmeliq_grid(i,k) * &
                 (pdel_grid(i,k) / gravit) / rhoh2o
         end if
      end do
   end do

   ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
   ! this is 1ppmv of h2o in 10hpa
   ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

   !-----------------------------------------------------------------------
   ! precipitation efficiency calculation  (accumulate cme and precip)

   minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

   ! zero out precip efficiency and total averaged precip
   pe_grid(:ngrdcol)     = 0._r8
   tpr_grid(:ngrdcol)    = 0._r8
   pefrac_grid(:ngrdcol) = 0._r8

   ! accumulate precip and condensation
   do i = 1, ngrdcol

      acgcme_grid(i)  = acgcme_grid(i) + tgcmeliq_grid(i)
      acprecl_grid(i) = acprecl_grid(i) + prec_str_grid(i)
      acnum_grid(i)   = acnum_grid(i) + 1

      ! if LWP is zero, then 'end of cloud': calculate precip efficiency
      if (tgliqwp_grid(i) < minlwp) then
         if (acprecl_grid(i) > 5.e-8_r8) then
            tpr_grid(i) = max(acprecl_grid(i)/acnum_grid(i), 1.e-15_r8)
            if (acgcme_grid(i) > 1.e-10_r8) then
               pe_grid(i) = min(max(acprecl_grid(i)/acgcme_grid(i), 1.e-15_r8), 1.e5_r8)
               pefrac_grid(i) = 1._r8
            end if
         end if

         ! reset counters
!        if (pe_grid(i) /= 0._r8 .and. (pe_grid(i) < 1.e-8_r8 .or. pe_grid(i) > 1.e3_r8)) then
!           write (iulog,*) 'PE_grid:ANOMALY  pe_grid, acprecl_grid, acgcme_grid, tpr_grid, acnum_grid ', &
!                           pe_grid(i),acprecl_grid(i), acgcme_grid(i), tpr_grid(i), acnum_grid(i)
!        endif

         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
         acnum_grid(i)   = 0
      end if               ! end LWP zero conditional

      ! if never find any rain....(after 10^3 timesteps...)
      if (acnum_grid(i) > 1000) then
         acnum_grid(i)   = 0
         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
      end if

   end do

   !-----------------------------------------------------------------------
   ! vertical average of non-zero accretion, autoconversion and ratio.
   ! vars: vprco_grid(i),vprao_grid(i),racau_grid(i),cnt_grid

   vprao_grid = 0._r8
   cnt_grid = 0
   do k = top_lev, pver
      vprao_grid(:ngrdcol) = vprao_grid(:ngrdcol) + prao_grid(:ngrdcol,k)
      where (prao_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
   end do

   where (cnt_grid > 0) vprao_grid = vprao_grid/cnt_grid

   vprco_grid = 0._r8
   cnt_grid = 0
   do k = top_lev, pver
      vprco_grid(:ngrdcol) = vprco_grid(:ngrdcol) + prco_grid(:ngrdcol,k)
      where (prco_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
   end do

   where (cnt_grid > 0)
      vprco_grid = vprco_grid/cnt_grid
      racau_grid = vprao_grid/vprco_grid
   elsewhere
      racau_grid = 0._r8
   end where

   racau_grid = min(racau_grid, 1.e10_r8)

   ! --------------------- !
   ! History Output Fields !
   ! --------------------- !

   ! Column droplet concentration
   cdnumc_grid(:ngrdcol) = sum(nc_grid(:ngrdcol,top_lev:pver) * &
        pdel_grid(:ngrdcol,top_lev:pver)/gravit, dim=2)

   ! Averaging for new output fields
   efcout_grid      = 0._r8
   efiout_grid      = 0._r8
   ncout_grid       = 0._r8
   niout_grid       = 0._r8
   freql_grid       = 0._r8
   freqi_grid       = 0._r8
   icwmrst_grid_out = 0._r8
   icimrst_grid_out = 0._r8
   freqm_grid       = 0._r8
   freqsl_grid      = 0._r8
   freqslm_grid     = 0._r8

   do k = top_lev, pver
      do i = 1, ngrdcol
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 5.e-5_r8 ) then
            efcout_grid(i,k) = rel_grid(i,k) * liqcldf_grid(i,k)
            ncout_grid(i,k)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            freql_grid(i,k)  = liqcldf_grid(i,k)
            icwmrst_grid_out(i,k) = icwmrst_grid(i,k)
         end if
         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-6_r8 ) then
            efiout_grid(i,k) = rei_grid(i,k) * icecldf_grid(i,k)
            niout_grid(i,k)  = icinc_grid(i,k) * icecldf_grid(i,k)
            freqi_grid(i,k)  = icecldf_grid(i,k)
            icimrst_grid_out(i,k) = icimrst_grid(i,k)
         end if

         ! Supercooled liquid
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) > 0.01_r8 ) then
            freqm_grid(i,k)=min(liqcldf_grid(i,k),icecldf_grid(i,k))
         end if
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) < 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
            freqsl_grid(i,k)=liqcldf_grid(i,k)
         end if
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) > 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
            freqslm_grid(i,k)=liqcldf_grid(i,k)
         end if

      end do
   end do

   ! Cloud top effective radius and number.
   fcti_grid  = 0._r8
   fctl_grid  = 0._r8
   ctrel_grid = 0._r8
   ctrei_grid = 0._r8
   ctnl_grid  = 0._r8
   ctni_grid  = 0._r8
   fctm_grid  = 0._r8
   fctsl_grid = 0._r8
   fctslm_grid= 0._r8

   do i = 1, ngrdcol
      do k = top_lev, pver
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 1.e-7_r8 ) then
            ctrel_grid(i) = rel_grid(i,k) * liqcldf_grid(i,k)
            ctnl_grid(i)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            fctl_grid(i)  = liqcldf_grid(i,k)

            ! Cloud Top Mixed phase, supercooled liquid only and supercooled liquid mixed
            if (freqi_grid(i,k) > 0.01_r8) then
               fctm_grid(i)=min(liqcldf_grid(i,k),icecldf_grid(i,k))
            end if
            if (freqi_grid(i,k) < 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
               fctsl_grid(i)=liqcldf_grid(i,k)
            end if
            if (freqi_grid(i,k) > 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
               fctslm_grid(i)=liqcldf_grid(i,k)
            end if

            exit
         end if

         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-7_r8 ) then
            ctrei_grid(i) = rei_grid(i,k) * icecldf_grid(i,k)
            ctni_grid(i)  = icinc_grid(i,k) * icecldf_grid(i,k)
            fcti_grid(i)  = icecldf_grid(i,k)
            exit
         end if
      end do
   end do

   ! Evaporation of stratiform precipitation fields for UNICON
   evprain_st_grid(:ngrdcol,:pver) = nevapr_grid(:ngrdcol,:pver) - evpsnow_st_grid(:ngrdcol,:pver)
   do k = top_lev, pver
      do i = 1, ngrdcol
         evprain_st_grid(i,k) = max(evprain_st_grid(i,k), 0._r8)
         evpsnow_st_grid(i,k) = max(evpsnow_st_grid(i,k), 0._r8)
      end do
   end do

   ! Assign the values to the pbuf pointers if they exist in pbuf
   if (qrain_idx > 0)  qrout_grid_ptr = qrout_grid
   if (qsnow_idx > 0)  qsout_grid_ptr = qsout_grid
   if (nrain_idx > 0)  nrout_grid_ptr = nrout_grid
   if (nsnow_idx > 0)  nsout_grid_ptr = nsout_grid
   if (qcsedten_idx > 0) qcsedtenout_grid_ptr = qcsedtenout_grid
   if (qrsedten_idx > 0) qrsedtenout_grid_ptr = qrsedtenout_grid
   if (qisedten_idx > 0) qisedtenout_grid_ptr = qisedtenout_grid
   if (qssedten_idx > 0) qssedtenout_grid_ptr = qssedtenout_grid
   if (vtrmc_idx > 0)    vtrmcout_grid_ptr    = vtrmcout_grid
   if (umr_idx > 0)      umrout_grid_ptr      = umrout_grid
   if (vtrmi_idx > 0)    vtrmiout_grid_ptr    = vtrmiout_grid
   if (ums_idx > 0)      umsout_grid_ptr      = umsout_grid
   if (qcsevap_idx > 0 ) qcsevapout_grid_ptr  = qcsevapout_grid
   if (qisevap_idx > 0 ) qisevapout_grid_ptr  = qisevapout_grid

   ! --------------------------------------------- !
   ! General outfield calls for microphysics       !
   ! --------------------------------------------- !

   ! Output a handle of variables which are calculated on the fly

   ftem_grid = 0._r8

   ftem_grid(:ngrdcol,top_lev:pver) =  qcreso_grid(:ngrdcol,top_lev:pver)
   call outfld( 'MPDW2V', ftem_grid, pcols, lchnk)

   ftem_grid(:ngrdcol,top_lev:pver) =  melto_grid(:ngrdcol,top_lev:pver) - mnuccco_grid(:ngrdcol,top_lev:pver)&
        - mnuccto_grid(:ngrdcol,top_lev:pver) -  bergo_grid(:ngrdcol,top_lev:pver) - homoo_grid(:ngrdcol,top_lev:pver)&
        - msacwio_grid(:ngrdcol,top_lev:pver)
   call outfld( 'MPDW2I', ftem_grid, pcols, lchnk)

   if (micro_mg_version > 2) then
      ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
          - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)&
          - psacwgo_grid(:ngrdcol,top_lev:pver) - pgsacwo_grid(:ngrdcol,top_lev:pver)
   else
      ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
          - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)
   endif

   call outfld( 'MPDW2P', ftem_grid, pcols, lchnk)

   ftem_grid(:ngrdcol,top_lev:pver) =  cmeiout_grid(:ngrdcol,top_lev:pver) + qireso_grid(:ngrdcol,top_lev:pver)
   call outfld( 'MPDI2V', ftem_grid, pcols, lchnk)

   if (micro_mg_version > 2) then
      ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
          + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
          + msacwio_grid(:ngrdcol,top_lev:pver)&
          - qmultgo_grid(:ngrdcol,top_lev:pver)
   else
      ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
          + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
          + msacwio_grid(:ngrdcol,top_lev:pver)
   endif

   call outfld( 'MPDI2W', ftem_grid, pcols, lchnk)

   ftem_grid(:ngrdcol,top_lev:pver) = -prcio_grid(:ngrdcol,top_lev:pver) - praio_grid(:ngrdcol,top_lev:pver)
   call outfld( 'MPDI2P', ftem_grid, pcols, lchnk)

   ! Output fields which have not been averaged already, averaging if use_subcol_microp is true
   if (trim(micro_mg_warm_rain) == 'tau' .or. trim(micro_mg_warm_rain) == 'emulated') then
      call outfld('scale_qc',    proc_rates%scale_qc,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('scale_nc',    proc_rates%scale_nc,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('scale_qr',    proc_rates%scale_qr,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('scale_nr',    proc_rates%scale_nr,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('amk_c',       proc_rates%amk_c,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ank_c',       proc_rates%ank_c,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('amk_r',       proc_rates%amk_r,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ank_r',       proc_rates%ank_r,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('amk',         proc_rates%amk,         ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ank',         proc_rates%ank,         ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('amk_out',     proc_rates%amk_out,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ank_out',     proc_rates%ank_out,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QC_TAU_out',  proc_rates%qc_out_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NC_TAU_out',  proc_rates%nc_out_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QR_TAU_out',  proc_rates%qr_out_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NR_TAU_out',  proc_rates%nr_out_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qctend_TAU',  proc_rates%qctend_TAU,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nctend_TAU',  proc_rates%nctend_TAU,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qrtend_TAU',  proc_rates%qrtend_TAU,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nrtend_TAU',  proc_rates%nrtend_TAU,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('gmnnn_lmnnn_TAU',  proc_rates%gmnnn_lmnnn_TAU,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ML_fixer',     proc_rates%ML_fixer,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qc_fixer',     proc_rates%qc_fixer,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nc_fixer',     proc_rates%nc_fixer,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qr_fixer',     proc_rates%qr_fixer,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nr_fixer',     proc_rates%nr_fixer,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QC_TAU_in',   proc_rates%qc_in_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NC_TAU_in',   proc_rates%nc_in_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QR_TAU_in',   proc_rates%qr_in_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NR_TAU_in',   proc_rates%nr_in_TAU,      ncol, lchnk, avg_subcol_field=use_subcol_microp)
   end if

   if (trim(micro_mg_warm_rain) == 'sb2001') then
      call outfld('qctend_SB2001',  proc_rates%qctend_SB2001,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nctend_SB2001',  proc_rates%nctend_SB2001,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qrtend_SB2001',  proc_rates%qrtend_SB2001,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nrtend_SB2001',  proc_rates%nrtend_SB2001,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
   end if
   if (trim(micro_mg_warm_rain) == 'kk2000') then
      call outfld('qctend_KK2000',  proc_rates%qctend_KK2000,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nctend_KK2000',  proc_rates%nctend_KK2000,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('qrtend_KK2000',  proc_rates%qrtend_KK2000,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('nrtend_KK2000',  proc_rates%nrtend_KK2000,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
   end if
   call outfld('LAMC',  proc_rates%lamc_out,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('LAMR',  proc_rates%lamr_out,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('PGAM',  proc_rates%pgam_out,  ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('N0R',  proc_rates%n0r_out,  ncol, lchnk, avg_subcol_field=use_subcol_microp)

   call outfld('MPICLWPI',    iclwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPICIWPI',    iciwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('REFL',        refl,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AREFL',       arefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AREFLZ',      areflz,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREFL',       frefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('CSRFL',       csrfl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ACSRFL',      acsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FCSRFL',      fcsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('REFL10CM',    refl10cm,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('REFLZ10CM',   reflz10cm,   psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('RERCLD',      rercld,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('NCAL',        ncal,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('NCAI',        ncai,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AQRAIN',      qrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AQSNOW',      qsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ANRAIN',      nrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ANSNOW',      nsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREQR',       freqr,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREQS',       freqs,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDT',        tlat,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDQ',        qvlat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDLIQ',      qcten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDICE',      qiten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDNLIQ',     ncten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDNICE',     niten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('EVAPSNOW',    proc_rates%evapsnow,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEVAP',     proc_rates%qcsevap,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEVAP',     proc_rates%qisevap,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QVRES',       proc_rates%qvres,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMC',       proc_rates%vtrmc,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMI',       proc_rates%vtrmi,       ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEDTEN',    proc_rates%qcsedten,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEDTEN',    proc_rates%qisedten,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QRSEDTEN',    proc_rates%qrsedten,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QSSEDTEN',    proc_rates%qssedten,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCRIO',    proc_rates%mnuccritot,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUDEPO',     proc_rates%mnudeptot,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MELTSTOT',    proc_rates%meltstot,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCDO',     proc_rates%mnuccdtot,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCDOhet',  mnuccdohet,  psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCRO',     proc_rates%mnuccrtot,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('PRACSO',      proc_rates%pracstot ,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VAPDEPSO',    proc_rates%vapdepstot,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MELTSDT',     proc_rates%meltsdttot,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FRZRDT',      proc_rates%frzrdttot ,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FICE',        nfice,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('CLDFSNOW',    cldfsnow,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUCCCO',  proc_rates%nnuccctot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUCCTO',  proc_rates%nnuccttot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUCCDO',  proc_rates%nnuccdtot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUDEPO',  proc_rates%nnudeptot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NHOMO',    proc_rates%nhomotot   , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUCCRO',  proc_rates%nnuccrtot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NNUCCRIO', proc_rates%nnuccritot , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NSACWIO',  proc_rates%nsacwitot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPRAO',    proc_rates%npratot    , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPSACWSO', proc_rates%npsacwstot , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPRAIO',   proc_rates%npraitot   , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPRACSO',  proc_rates%npracstot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPRCO',    proc_rates%nprctot    , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NPRCIO',   proc_rates%nprcitot   , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NCSEDTEN', proc_rates%ncsedten , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NISEDTEN', proc_rates%nisedten , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NRSEDTEN', proc_rates%nrsedten , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NSSEDTEN', proc_rates%nssedten , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NMELTO',   proc_rates%nmelttot   , ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld ('NMELTS',   proc_rates%nmeltstot  , ncol, lchnk, avg_subcol_field=use_subcol_microp)

   call outfld('UMR',      proc_rates%umr,         ncol, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('UMS',      proc_rates%ums,         ncol, lchnk, avg_subcol_field=use_subcol_microp)

   call outfld('QCRAT',    qcrat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)

   if (micro_mg_version > 2) then
      call outfld('UMG',        proc_rates%umg,         ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QGSEDTEN',   proc_rates%qgsedten,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('FREQG',       freqg,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('AQGRAU',      qgout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('ANGRAU',      ngout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('CLDFGRAU',    cldfgrau,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('MELTGTOT',    proc_rates%meltgtot,    ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NMELTG',      proc_rates%nmeltgtot,     ncol, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('NGSEDTEN',    proc_rates%ngsedten ,   ncol, lchnk, avg_subcol_field=use_subcol_microp)

   end if

   ! Example subcolumn outfld call
   if (use_subcol_microp) then
      call outfld('FICE_SCOL',   nfice,       psubcols*pcols, lchnk)
      call outfld('MPDLIQ_SCOL', ptend%q(:,:,ixcldliq),       psubcols*pcols, lchnk)
      call outfld('MPDICE_SCOL', qiten,       psubcols*pcols, lchnk)
   end if

   ! Output fields which are already on the grid
   call outfld('QRAIN',       qrout_grid,       pcols, lchnk)
   call outfld('QSNOW',       qsout_grid,       pcols, lchnk)
   call outfld('NRAIN',       nrout_grid,       pcols, lchnk)
   call outfld('NSNOW',       nsout_grid,       pcols, lchnk)
   call outfld('CV_REFFLIQ',  cvreffliq_grid,   pcols, lchnk)
   call outfld('CV_REFFICE',  cvreffice_grid,   pcols, lchnk)
   call outfld('LS_FLXPRC',   mgflxprc_grid,    pcols, lchnk)
   call outfld('LS_FLXSNW',   mgflxsnw_grid,    pcols, lchnk)
   call outfld('CME',         qme_grid,         pcols, lchnk)
   call outfld('PRODPREC',    prain_grid,       pcols, lchnk)
   call outfld('EVAPPREC',    nevapr_grid,      pcols, lchnk)
   call outfld('QCRESO',      qcreso_grid,      pcols, lchnk)
   call outfld('LS_REFFRAIN', mgreffrain_grid,  pcols, lchnk)
   call outfld('LS_REFFSNOW', mgreffsnow_grid,  pcols, lchnk)
   call outfld('DSNOW',       des_grid,         pcols, lchnk)
   call outfld('ADRAIN',      drout2_grid,      pcols, lchnk)
   call outfld('ADSNOW',      dsout2_grid,      pcols, lchnk)
   call outfld('PE',          pe_grid,          pcols, lchnk)
   call outfld('PEFRAC',      pefrac_grid,      pcols, lchnk)
   call outfld('APRL',        tpr_grid,         pcols, lchnk)
   call outfld('VPRAO',       vprao_grid,       pcols, lchnk)
   call outfld('VPRCO',       vprco_grid,       pcols, lchnk)
   call outfld('RACAU',       racau_grid,       pcols, lchnk)
   call outfld('AREL',        efcout_grid,      pcols, lchnk)
   call outfld('AREI',        efiout_grid,      pcols, lchnk)
   call outfld('AWNC' ,       ncout_grid,       pcols, lchnk)
   call outfld('AWNI' ,       niout_grid,       pcols, lchnk)
   call outfld('FREQL',       freql_grid,       pcols, lchnk)
   call outfld('FREQI',       freqi_grid,       pcols, lchnk)
   call outfld('ACTREL',      ctrel_grid,       pcols, lchnk)
   call outfld('ACTREI',      ctrei_grid,       pcols, lchnk)
   call outfld('ACTNL',       ctnl_grid,        pcols, lchnk)
   call outfld('ACTNI',       ctni_grid,        pcols, lchnk)
   call outfld('FCTL',        fctl_grid,        pcols, lchnk)
   call outfld('FCTI',        fcti_grid,        pcols, lchnk)
   call outfld('ICINC',       icinc_grid,       pcols, lchnk)
   call outfld('ICWNC',       icwnc_grid,       pcols, lchnk)
   call outfld('EFFLIQ_IND',  rel_fn_grid,      pcols, lchnk)
   call outfld('CDNUMC',      cdnumc_grid,      pcols, lchnk)
   call outfld('REL',         rel_grid,         pcols, lchnk)
   call outfld('REI',         rei_grid,         pcols, lchnk)
   call outfld('MG_SADICE',   sadice_grid,      pcols, lchnk)
   call outfld('MG_SADSNOW',  sadsnow_grid,     pcols, lchnk)
   call outfld('ICIMRST',     icimrst_grid_out, pcols, lchnk)
   call outfld('ICWMRST',     icwmrst_grid_out, pcols, lchnk)
   call outfld('CMEIOUT',     cmeiout_grid,     pcols, lchnk)
   call outfld('PRAO',        prao_grid,        pcols, lchnk)
   call outfld('PRCO',        prco_grid,        pcols, lchnk)
   call outfld('MNUCCCO',     mnuccco_grid,     pcols, lchnk)
   call outfld('MNUCCTO',     mnuccto_grid,     pcols, lchnk)
   call outfld('MSACWIO',     msacwio_grid,     pcols, lchnk)
   call outfld('PSACWSO',     psacwso_grid,     pcols, lchnk)
   call outfld('BERGSO',      bergso_grid,      pcols, lchnk)
   call outfld('BERGO',       bergo_grid,       pcols, lchnk)
   call outfld('MELTO',       melto_grid,       pcols, lchnk)
   call outfld('HOMOO',       homoo_grid,       pcols, lchnk)
   call outfld('PRCIO',       prcio_grid,       pcols, lchnk)
   call outfld('PRAIO',       praio_grid,       pcols, lchnk)
   call outfld('QIRESO',      qireso_grid,      pcols, lchnk)
   call outfld('FREQM',       freqm_grid,       pcols, lchnk)
   call outfld('FREQSL',      freqsl_grid,      pcols, lchnk)
   call outfld('FREQSLM',     freqslm_grid,     pcols, lchnk)
   call outfld('FCTM',        fctm_grid,        pcols, lchnk)
   call outfld('FCTSL',       fctsl_grid,       pcols, lchnk)
   call outfld('FCTSLM',      fctslm_grid,      pcols, lchnk)

   if (micro_mg_version > 2) then
      call outfld('PRACGO',      pracgo_grid,      pcols, lchnk)
      call outfld('PSACRO',      psacro_grid,      pcols, lchnk)
      call outfld('PSACWGO',     psacwgo_grid,     pcols, lchnk)
      call outfld('PGSACWO',     pgsacwo_grid,     pcols, lchnk)
      call outfld('PGRACSO',     pgracso_grid,     pcols, lchnk)
      call outfld('PRDGO',       prdgo_grid,       pcols, lchnk)
      call outfld('QMULTGO',     qmultgo_grid,     pcols, lchnk)
      call outfld('QMULTRGO',    qmultrgo_grid,    pcols, lchnk)
      call outfld('LS_REFFGRAU', reff_grau_grid,  pcols, lchnk)
      call outfld ('NPRACGO',    npracgo_grid,  pcols, lchnk)
      call outfld ('NSCNGO',     nscngo_grid,  pcols, lchnk)
      call outfld ('NGRACSO',    ngracso_grid,  pcols, lchnk)
      call outfld ('NMULTGO',    nmultgo_grid,  pcols, lchnk)
      call outfld ('NMULTRGO',   nmultrgo_grid,  pcols, lchnk)
      call outfld ('NPSACWGO',   npsacwgo_grid,  pcols, lchnk)
   end if

   if (micro_mg_adjust_cpt) then
      cp_rh(:ncol, :pver)  = 0._r8

      do i = 1, ncol

         ! Calculate the RH including any T change that we make.
         do k = top_lev, pver
           call qsat(state_loc%t(i,k), state_loc%pmid(i,k), es, qs)
           cp_rh(i,k) = state_loc%q(i, k, ixq) / qs * 100._r8
         end do
      end do

      call outfld("TROPF_RHADJ", cp_rh,       pcols, lchnk)
   end if

   ! deallocate the temporary pbuf grid variable which was allocated if subcolumns are not used
   if (.not. use_subcol_microp) then
      deallocate(bergso_grid)
   end if

   ! deallocate the proc_rates DDT
   call proc_rates%deallocate(micro_mg_warm_rain)

   ! ptend_loc is deallocated in physics_update above
   call physics_state_dealloc(state_loc)

   if (qsatfac_idx <= 0) then
      deallocate(qsatfac)
   end if

end subroutine micro_pumas_cam_tend

subroutine massless_droplet_destroyer(ztodt, state,  ptend)

     ! This subroutine eradicates cloud droplets in grid boxes with no cloud
     ! mass.  This code is now expanded to remove massless rain drops, ice
     ! crystals, and snow flakes.
     !
     ! Note: qsmall, which is a small, positive number, is used as the
     !       threshold here instead of qmin, which is 0.  Some numbers that are
     !       supposed to have a value of 0, but don't because of numerical
     !       roundoff (especially after hole filling) will have small, positive
     !       values.  Using qsmall as the threshold here instead of qmin allows
     !       for unreasonable massless drop concentrations to be removed in
     !       those scenarios.

     use micro_pumas_utils,   only: qsmall
     use ref_pres,         only: top_lev => trop_cloud_top_lev

     implicit none

     ! Input Variables
     real(r8), intent(in)                  :: ztodt     ! model time increment
     type(physics_state), intent(in)       :: state     ! state for columns

     ! Input/Output Variables
     type(physics_ptend), intent(inout)    :: ptend     ! ptend for columns

     ! Local Variables
     integer :: icol, k

     !----- Begin Code -----

     ! Don't do anything if this option isn't enabled.
     if ( .not. micro_do_massless_droplet_destroyer ) return

     col_loop: do icol=1, state%ncol
       vert_loop: do k = top_lev, pver
         ! If updated qc (after microphysics) is zero, then ensure updated nc is also zero!!
         if ( state%q(icol,k,ixcldliq) + ztodt * ptend%q(icol,k,ixcldliq) < qsmall ) then
           ptend%lq(ixnumliq) = .true. ! This is probably already true, but it doesn't
                                       ! hurt to set it.
           ptend%q(icol,k,ixnumliq) = -(state%q(icol,k,ixnumliq) / ztodt)
         end if
         if ( ixnumrain > 0 ) then
            ! If updated qr (after microphysics) is zero, then ensure updated nr is also zero!!
            if ( state%q(icol,k,ixrain) + ztodt * ptend%q(icol,k,ixrain) < qsmall ) then
              ptend%lq(ixnumrain) = .true. ! This is probably already true, but it doesn't
                                           ! hurt to set it.
              ptend%q(icol,k,ixnumrain) = -(state%q(icol,k,ixnumrain) / ztodt)
            end if
         endif ! ixnumrain > 0
         ! If updated qi (after microphysics) is zero, then ensure updated ni is also zero!!
         if ( state%q(icol,k,ixcldice) + ztodt * ptend%q(icol,k,ixcldice) < qsmall ) then
           ptend%lq(ixnumice) = .true. ! This is probably already true, but it doesn't
                                       ! hurt to set it.
           ptend%q(icol,k,ixnumice) = -(state%q(icol,k,ixnumice) / ztodt)
         end if
         if ( ixnumsnow > 0 ) then
            ! If updated qs (after microphysics) is zero, then ensure updated ns is also zero!!
            if ( state%q(icol,k,ixsnow) + ztodt * ptend%q(icol,k,ixsnow) < qsmall ) then
              ptend%lq(ixnumsnow) = .true. ! This is probably already true, but it doesn't
                                           ! hurt to set it.
              ptend%q(icol,k,ixnumsnow) = -(state%q(icol,k,ixnumsnow) / ztodt)
            end if
         endif ! ixnumsnow > 0
       end do vert_loop
     end do col_loop

     return
end subroutine massless_droplet_destroyer

end module micro_pumas_cam
