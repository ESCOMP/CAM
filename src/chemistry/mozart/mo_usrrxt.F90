module mo_usrrxt

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_logfile,      only : iulog
  use ppgrid,           only : pver, pcols
  use cam_abortutils,   only : endrun

  implicit none

  private
  public :: usrrxt, usrrxt_inti, usrrxt_hrates
  public :: has_ice_trp_rxts

  save

  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_CO_OH_ndx
  integer :: usr_CO_OH_a_ndx
  integer :: usr_CO_OH_b_ndx
  integer :: usr_PAN_M_ndx
  integer :: usr_CH3COCH3_OH_ndx
  integer :: usr_MCO3_NO2_ndx
  integer :: usr_MPAN_M_ndx
  integer :: usr_XOOH_OH_ndx
  integer :: usr_SO2_OH_ndx
  integer :: usr_DMS_OH_ndx
  integer :: usr_HO2_aer_ndx
  integer :: usr_GLYOXAL_aer_ndx
  integer :: usr_N_O2_ndx
  integer :: usr_N2D_O2_ndx
  integer :: usr_N2D_e_ndx

  integer :: tag_NO2_NO3_ndx
  integer :: tag_NO2_OH_ndx
  integer :: tag_NO2_HO2_ndx
  integer :: tag_C2H4_OH_ndx
  integer :: tag_C3H6_OH_ndx
  integer :: tag_CH3CO3_NO2_ndx

!lke-TS1
  integer :: usr_PBZNIT_M_ndx
  integer :: tag_ACBZO2_NO2_ndx
  integer :: usr_ISOPNITA_aer_ndx
  integer :: usr_ISOPNITB_aer_ndx
  integer :: usr_ONITR_aer_ndx
  integer :: usr_HONITR_aer_ndx
  integer :: usr_TERPNIT_aer_ndx
  integer :: usr_NTERPOOH_aer_ndx
  integer :: usr_NC4CHO_aer_ndx
  integer :: usr_NC4CH2OH_aer_ndx
!TS2
 integer :: usr_ISOPZD1O2_ndx
 integer :: usr_ISOPZD4O2_ndx
 integer :: usr_ISOPFDN_aer_ndx
 integer :: usr_ISOPFNP_aer_ndx
 integer :: usr_ISOPN2B_aer_ndx
 integer :: usr_ISOPN1D_aer_ndx
 integer :: usr_ISOPN4D_aer_ndx
 integer :: usr_INOOHD_aer_ndx
 integer :: usr_INHEB_aer_ndx
 integer :: usr_INHED_aer_ndx
 integer :: usr_MACRN_aer_ndx
 integer :: usr_ISOPHFP_aer_ndx
 integer :: usr_IEPOX_aer_ndx
 integer :: usr_DHPMPAL_aer_ndx
 integer :: usr_ICHE_aer_ndx
 integer :: usr_ISOPFNC_aer_ndx
 integer :: usr_ISOPFDNC_aer_ndx
 integer :: usr_ISOPB1O2_NOn_ndx
 integer :: usr_ISOPB1O2_NOa_ndx
 integer :: usr_ISOPB4O2_NOn_ndx
 integer :: usr_ISOPB4O2_NOa_ndx
 integer :: usr_ISOPED1O2_NOn_ndx
 integer :: usr_ISOPED1O2_NOa_ndx
 integer :: usr_ISOPED4O2_NOn_ndx
 integer :: usr_ISOPED4O2_NOa_ndx
 integer :: usr_ISOPZD1O2_NOn_ndx
 integer :: usr_ISOPZD1O2_NOa_ndx
 integer :: usr_ISOPZD4O2_NOn_ndx
 integer :: usr_ISOPZD4O2_NOa_ndx
 integer :: usr_ISOPNO3_NOn_ndx
 integer :: usr_ISOPNO3_NOa_ndx
 integer :: usr_MVKO2_NOn_ndx
 integer :: usr_MVKO2_NOa_ndx
 integer :: usr_MACRO2_NOn_ndx
 integer :: usr_MACRO2_NOa_ndx
 integer :: usr_IEPOXOO_NOn_ndx
 integer :: usr_IEPOXOO_NOa_ndx
 integer :: usr_ISOPN1DO2_NOn_ndx
 integer :: usr_ISOPN1DO2_NOa_ndx
 integer :: usr_ISOPN2BO2_NOn_ndx
 integer :: usr_ISOPN2BO2_NOa_ndx
 integer :: usr_ISOPN3BO2_NOn_ndx
 integer :: usr_ISOPN3BO2_NOa_ndx
 integer :: usr_ISOPN4DO2_NOn_ndx
 integer :: usr_ISOPN4DO2_NOa_ndx
 integer :: usr_ISOPNBNO3O2_NOn_ndx
 integer :: usr_ISOPNBNO3O2_NOa_ndx
 integer :: usr_ISOPNOOHBO2_NOn_ndx
 integer :: usr_ISOPNOOHBO2_NOa_ndx
 integer :: usr_ISOPNOOHDO2_NOn_ndx
 integer :: usr_ISOPNOOHDO2_NOa_ndx
 integer :: usr_NC4CHOO2_NOn_ndx
 integer :: usr_NC4CHOO2_NOa_ndx
 integer :: tag_MCO3_NO2_ndx
 integer :: tag_TERPACO3_NO2_ndx
 integer :: usr_TERPAPAN_M_ndx
 integer :: tag_TERPA2CO3_NO2_ndx
 integer :: usr_TERPA2PAN_M_ndx
 integer :: tag_TERPA3CO3_NO2_ndx
 integer :: usr_TERPA3PAN_M_ndx
 integer :: usr_TERPNT_aer_ndx
 integer :: usr_TERPNT1_aer_ndx
 integer :: usr_TERPNPT_aer_ndx
 integer :: usr_TERPNPT1_aer_ndx
 integer :: usr_TERPFDN_aer_ndx
 integer :: usr_SQTN_aer_ndx
 integer :: usr_TERPHFN_aer_ndx
 integer :: usr_TERPDHDP_aer_ndx
 integer :: usr_TERPACID_aer_ndx
!
  integer :: usr_OA_O2_NDX
  integer :: usr_XNO2NO3_M_ndx
  integer :: usr_NO2XNO3_M_ndx
  integer :: usr_XHNO3_OH_ndx
  integer :: usr_XHO2NO2_M_ndx
  integer :: usr_XNO2NO3_aer_ndx
  integer :: usr_NO2XNO3_aer_ndx
  integer :: usr_XNO3_aer_ndx
  integer :: usr_XNO2_aer_ndx
  integer :: usr_XPAN_M_ndx
  integer :: usr_XMPAN_M_ndx
  integer :: usr_MCO3_XNO2_ndx

  integer :: usr_C2O3_NO2_ndx
  integer :: usr_C2H4_OH_ndx
  integer :: usr_XO2N_HO2_ndx
  integer :: usr_C2O3_XNO2_ndx

  integer :: tag_XO2N_NO_ndx
  integer :: tag_XO2_HO2_ndx
  integer :: tag_XO2_NO_ndx

  integer :: usr_O_O_ndx
  integer :: usr_CL2O2_M_ndx
  integer :: usr_SO3_H2O_ndx
  integer :: tag_CLO_CLO_M_ndx

  integer :: ion1_ndx, ion2_ndx, ion3_ndx, ion11_ndx
  integer :: elec1_ndx, elec2_ndx, elec3_ndx
  integer :: elec4_ndx, elec5_ndx, elec6_ndx
  integer :: het1_ndx

  integer, parameter :: nean = 3
  integer :: ean_ndx(nean)
  integer, parameter :: nrpe = 5
  integer :: rpe_ndx(nrpe)
  integer, parameter :: npir = 16
  integer :: pir_ndx(npir)
  integer, parameter :: nedn = 2
  integer :: edn_ndx(nedn)
  integer, parameter :: nnir = 13
  integer :: nir_ndx(nnir)
  integer, parameter :: niira = 112
  integer :: iira_ndx(niira)
  integer, parameter :: niirb = 14
  integer :: iirb_ndx(niirb)

  integer :: usr_clm_h2o_m_ndx, usr_clm_hcl_m_ndx
  integer :: usr_oh_co_ndx, het_no2_h2o_ndx, usr_oh_dms_ndx, aq_so2_h2o2_ndx, aq_so2_o3_ndx

  integer :: h2o_ndx
!
! jfl
!
  integer, parameter :: num_strat_tau = 22
  integer :: usr_strat_tau_ndx(num_strat_tau)
!
!lke++
  integer :: usr_COhc_OH_ndx
  integer :: usr_COme_OH_ndx
  integer :: usr_CO01_OH_ndx
  integer :: usr_CO02_OH_ndx
  integer :: usr_CO03_OH_ndx
  integer :: usr_CO04_OH_ndx
  integer :: usr_CO05_OH_ndx
  integer :: usr_CO06_OH_ndx
  integer :: usr_CO07_OH_ndx
  integer :: usr_CO08_OH_ndx
  integer :: usr_CO09_OH_ndx
  integer :: usr_CO10_OH_ndx
  integer :: usr_CO11_OH_ndx
  integer :: usr_CO12_OH_ndx
  integer :: usr_CO13_OH_ndx
  integer :: usr_CO14_OH_ndx
  integer :: usr_CO15_OH_ndx
  integer :: usr_CO16_OH_ndx
  integer :: usr_CO17_OH_ndx
  integer :: usr_CO18_OH_ndx
  integer :: usr_CO19_OH_ndx
  integer :: usr_CO20_OH_ndx
  integer :: usr_CO21_OH_ndx
  integer :: usr_CO22_OH_ndx
  integer :: usr_CO23_OH_ndx
  integer :: usr_CO24_OH_ndx
  integer :: usr_CO25_OH_ndx
  integer :: usr_CO26_OH_ndx
  integer :: usr_CO27_OH_ndx
  integer :: usr_CO28_OH_ndx
  integer :: usr_CO29_OH_ndx
  integer :: usr_CO30_OH_ndx
  integer :: usr_CO31_OH_ndx
  integer :: usr_CO32_OH_ndx
  integer :: usr_CO33_OH_ndx
  integer :: usr_CO34_OH_ndx
  integer :: usr_CO35_OH_ndx
  integer :: usr_CO36_OH_ndx
  integer :: usr_CO37_OH_ndx
  integer :: usr_CO38_OH_ndx
  integer :: usr_CO39_OH_ndx
  integer :: usr_CO40_OH_ndx
  integer :: usr_CO41_OH_ndx
  integer :: usr_CO42_OH_ndx
!lke--

! ================================================
! Halogen recycling on seasalts
! orginally from ordc
! ================================================
  integer :: het_ss_0_ndx
  integer :: het_ss_1_ndx
  integer :: het_ss_2_ndx
  integer :: het_ss_3_ndx
  integer :: het_ss_4_ndx
  integer :: het_ss_5_ndx
  integer :: het_ss_6_ndx
  integer :: het_ss_7_ndx
  integer :: het_ss_8_ndx
  integer :: het_ss_9_ndx
  integer :: het_ss_10_ndx
  integer :: het_ss_11_ndx
  integer :: het_ss_12_ndx

  integer :: ss_ixoy_2_ndx
  integer :: ss_ixoy_3_ndx
  integer :: ss_ixoy_4_ndx
  integer :: sslt1_ndx,sslt2_ndx,sslt3_ndx,sslt4_ndx

  integer :: id_hocl
  integer :: id_hcl
  integer :: id_hbr
  integer :: id_hi
  integer :: id_hobr
  integer :: id_hoi
  integer :: id_clono2
  integer :: id_brono2
  integer :: id_iono2
  integer :: id_n2o5

  integer :: usr_N2O5_aer1_ndx
  integer :: usr_N2O5_aer2_ndx
  integer :: usr_N2O5_HCL_ndx

  integer :: ice_trp_cl_1_ndx
  integer :: ice_trp_br_1_ndx
  integer :: ice_trp_i_1_ndx
  integer :: ice_trp_i_2_ndx
  integer :: ice_trp_i_3_ndx
  integer :: ice_trp_i_4_ndx
  integer :: ice_trp_hbr_5_ndx
  integer :: ice_trp_hbr_6_ndx
  integer :: ice_trp_hcl_5_ndx
  integer :: ice_trp_hcl_6_ndx
  integer :: ice_trp_hi_5_ndx
  integer :: ice_trp_hi_6_ndx

  integer :: ice_fr_hoi_ndx
  integer :: liq_fr_hoi_ndx
  integer :: ice_fr_hi_ndx
  integer :: liq_fr_hi_ndx
  integer :: ice_fr_iono2_ndx
  integer :: liq_fr_iono2_ndx
  integer :: ice_fr_brono2_ndx

  integer :: usr_IO_IO_a_ndx
  integer :: usr_IO_IO_b_ndx
  integer :: usr_IO_OIO_ndx
  integer :: usr_OIO_OIO_ndx
  integer :: usr_HOI_NO3_ndx
  integer :: usr_I2O2_a_ndx
  integer :: usr_I2O2_b_ndx
  integer :: usr_I2O4_ndx
  integer :: usr_IONO2_ndx

  logical, protected :: has_ice_trp_rxts
  logical :: has_het_ss_rxts
  logical :: has_ss_ixoy_rxts
  logical :: has_aerosols

  real(r8), parameter :: t0     = 300._r8                ! K
  real(r8), parameter :: trlim2 = 17._r8/3._r8           ! K
  real(r8), parameter :: trlim3 = 15._r8/3._r8           ! K

  logical :: has_ion_rxts, has_d_chem

contains

  subroutine usrrxt_inti
    !-----------------------------------------------------------------
    !        ... intialize the user reaction constants module
    !-----------------------------------------------------------------

    use mo_chem_utls,   only : get_rxt_ndx, get_spc_ndx
    use spmd_utils,     only : masterproc

    implicit none

    character(len=4) :: xchar
    character(len=32) :: rxtname
    integer :: i

!
! full tropospheric chemistry
!
    usr_O_O2_ndx         = get_rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = get_rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = get_rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = get_rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = get_rxt_ndx( 'usr_HO2NO2_M' )
    usr_N2O5_aer_ndx     = get_rxt_ndx( 'usr_N2O5_aer' )
    usr_N2O5_aer1_ndx    = get_rxt_ndx( 'usr_N2O5_aer1' )
    usr_N2O5_aer2_ndx    = get_rxt_ndx( 'usr_N2O5_aer2' )
    usr_N2O5_HCL_ndx     = get_rxt_ndx( 'usr_N2O5_HCL' )
    usr_NO3_aer_ndx      = get_rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = get_rxt_ndx( 'usr_NO2_aer' )
    usr_CO_OH_ndx        = get_rxt_ndx( 'usr_CO_OH' )
    usr_CO_OH_a_ndx      = get_rxt_ndx( 'usr_CO_OH_a' )
    usr_CO_OH_b_ndx      = get_rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = get_rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = get_rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = get_rxt_ndx( 'usr_MCO3_NO2' )
    tag_MCO3_NO2_ndx     = get_rxt_ndx( 'tag_MCO3_NO2' )
    if( tag_MCO3_NO2_ndx<0 .and. usr_MCO3_NO2_ndx>0 ) then
       tag_MCO3_NO2_ndx = usr_MCO3_NO2_ndx
    endif

    usr_MPAN_M_ndx       = get_rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = get_rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = get_rxt_ndx( 'usr_SO2_OH' )
    usr_N_O2_ndx         = get_rxt_ndx( 'usr_N_O2' )
    usr_N2D_O2_ndx       = get_rxt_ndx( 'usr_N2D_O2' )
    usr_N2D_e_ndx        = get_rxt_ndx( 'usr_N2D_e' )
    usr_DMS_OH_ndx       = get_rxt_ndx( 'usr_DMS_OH' )
    usr_HO2_aer_ndx      = get_rxt_ndx( 'usr_HO2_aer' )
    usr_GLYOXAL_aer_ndx  = get_rxt_ndx( 'usr_GLYOXAL_aer' )
 !
    tag_NO2_NO3_ndx      = get_rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_OH_ndx       = get_rxt_ndx( 'tag_NO2_OH' )
    tag_NO2_HO2_ndx      = get_rxt_ndx( 'tag_NO2_HO2' )
    tag_C2H4_OH_ndx      = get_rxt_ndx( 'tag_C2H4_OH' )
    tag_C3H6_OH_ndx      = get_rxt_ndx( 'tag_C3H6_OH' )
    tag_CH3CO3_NO2_ndx   = get_rxt_ndx( 'tag_CH3CO3_NO2' )
!lke-TS1
    usr_PBZNIT_M_ndx     = get_rxt_ndx( 'usr_PBZNIT_M' )
    tag_ACBZO2_NO2_ndx   = get_rxt_ndx( 'tag_ACBZO2_NO2' )
    usr_ISOPNITA_aer_ndx = get_rxt_ndx( 'usr_ISOPNITA_aer' )
    usr_ISOPNITB_aer_ndx = get_rxt_ndx( 'usr_ISOPNITB_aer' )
    usr_ONITR_aer_ndx    = get_rxt_ndx( 'usr_ONITR_aer' )
    usr_HONITR_aer_ndx   = get_rxt_ndx( 'usr_HONITR_aer' )
    usr_TERPNIT_aer_ndx  = get_rxt_ndx( 'usr_TERPNIT_aer' )
    usr_NTERPOOH_aer_ndx = get_rxt_ndx( 'usr_NTERPOOH_aer' )
    usr_NC4CHO_aer_ndx   = get_rxt_ndx( 'usr_NC4CHO_aer' )
    usr_NC4CH2OH_aer_ndx = get_rxt_ndx( 'usr_NC4CH2OH_aer' )
!TS2
    usr_ISOPZD1O2_ndx        = get_rxt_ndx( 'usr_ISOPZD1O2' )
    usr_ISOPZD4O2_ndx        = get_rxt_ndx( 'usr_ISOPZD4O2' )
    usr_ISOPFDN_aer_ndx      = get_rxt_ndx( 'usr_ISOPFDN_aer' )
    usr_ISOPFNP_aer_ndx      = get_rxt_ndx( 'usr_ISOPFNP_aer' )
    usr_ISOPN2B_aer_ndx      = get_rxt_ndx( 'usr_ISOPN2B_aer' )
    usr_ISOPN1D_aer_ndx      = get_rxt_ndx( 'usr_ISOPN1D_aer' )
    usr_ISOPN4D_aer_ndx      = get_rxt_ndx( 'usr_ISOPN4D_aer' )
    usr_INOOHD_aer_ndx       = get_rxt_ndx( 'usr_INOOHD_aer' )
    usr_INHEB_aer_ndx        = get_rxt_ndx( 'usr_INHEB_aer' )
    usr_INHED_aer_ndx        = get_rxt_ndx( 'usr_INHED_aer' )
    usr_MACRN_aer_ndx        = get_rxt_ndx( 'usr_MACRN_aer' )
    usr_ISOPHFP_aer_ndx      = get_rxt_ndx( 'usr_ISOPHFP_aer' )
    usr_IEPOX_aer_ndx        = get_rxt_ndx( 'usr_IEPOX_aer' )
    usr_DHPMPAL_aer_ndx      = get_rxt_ndx( 'usr_DHPMPAL_aer' )
    usr_ICHE_aer_ndx         = get_rxt_ndx( 'usr_ICHE_aer' )
    usr_ISOPFNC_aer_ndx      = get_rxt_ndx( 'usr_ISOPFNC_aer' )
    usr_ISOPFDNC_aer_ndx     = get_rxt_ndx( 'usr_ISOPFDNC_aer' )
    usr_ISOPB1O2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPB1O2_NOn' )
    usr_ISOPB1O2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPB1O2_NOa' )
    usr_ISOPB4O2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPB4O2_NOn' )
    usr_ISOPB4O2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPB4O2_NOa' )
    usr_ISOPED1O2_NOn_ndx    = get_rxt_ndx( 'usr_ISOPED1O2_NOn' )
    usr_ISOPED1O2_NOa_ndx    = get_rxt_ndx( 'usr_ISOPED1O2_NOa' )
    usr_ISOPED4O2_NOn_ndx    = get_rxt_ndx( 'usr_ISOPED4O2_NOn' )
    usr_ISOPED4O2_NOa_ndx    = get_rxt_ndx( 'usr_ISOPED4O2_NOa' )
    usr_ISOPZD1O2_NOn_ndx    = get_rxt_ndx( 'usr_ISOPZD1O2_NOn' )
    usr_ISOPZD1O2_NOa_ndx    = get_rxt_ndx( 'usr_ISOPZD1O2_NOa' )
    usr_ISOPZD4O2_NOn_ndx    = get_rxt_ndx( 'usr_ISOPZD4O2_NOn' )
    usr_ISOPZD4O2_NOa_ndx    = get_rxt_ndx( 'usr_ISOPZD4O2_NOa' )
    usr_ISOPNO3_NOn_ndx      = get_rxt_ndx( 'usr_ISOPNO3_NOn' )
    usr_ISOPNO3_NOa_ndx      = get_rxt_ndx( 'usr_ISOPNO3_NOa' )
    usr_MVKO2_NOn_ndx        = get_rxt_ndx( 'usr_MVKO2_NOn' )
    usr_MVKO2_NOa_ndx        = get_rxt_ndx( 'usr_MVKO2_NOa' )
    usr_MACRO2_NOn_ndx       = get_rxt_ndx( 'usr_MACRO2_NOn' )
    usr_MACRO2_NOa_ndx       = get_rxt_ndx( 'usr_MACRO2_NOa' )
    usr_IEPOXOO_NOn_ndx     = get_rxt_ndx( 'usr_IEPOXOO_NOn' )
    usr_IEPOXOO_NOa_ndx     = get_rxt_ndx( 'usr_IEPOXOO_NOa' )
    usr_ISOPN1DO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPN1DO2_NOn' )
    usr_ISOPN1DO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPN1DO2_NOa' )
    usr_ISOPN2BO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPN2BO2_NOn' )
    usr_ISOPN2BO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPN2BO2_NOa' )
    usr_ISOPN3BO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPN3BO2_NOn' )
    usr_ISOPN3BO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPN3BO2_NOa' )
    usr_ISOPN4DO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPN4DO2_NOn' )
    usr_ISOPN4DO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPN4DO2_NOa' )
    usr_ISOPNBNO3O2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPNBNO3O2_NOn' )
    usr_ISOPNBNO3O2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPNBNO3O2_NOa' )
    usr_ISOPNOOHBO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPNOOHBO2_NOn' )
    usr_ISOPNOOHBO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPNOOHBO2_NOa' )
    usr_ISOPNOOHDO2_NOn_ndx     = get_rxt_ndx( 'usr_ISOPNOOHDO2_NOn' )
    usr_ISOPNOOHDO2_NOa_ndx     = get_rxt_ndx( 'usr_ISOPNOOHDO2_NOa' )
    usr_NC4CHOO2_NOn_ndx     = get_rxt_ndx( 'usr_NC4CHOO2_NOn' )
    usr_NC4CHOO2_NOa_ndx     = get_rxt_ndx( 'usr_NC4CHOO2_NOa' )
    tag_TERPACO3_NO2_ndx  = get_rxt_ndx( 'tag_TERPACO3_NO2' )
    usr_TERPAPAN_M_ndx     = get_rxt_ndx( 'usr_TERPAPAN_M' )
    tag_TERPA2CO3_NO2_ndx  = get_rxt_ndx( 'tag_TERPA2CO3_NO2' )
    usr_TERPA2PAN_M_ndx     = get_rxt_ndx( 'usr_TERPA2PAN_M' )
    tag_TERPA3CO3_NO2_ndx  = get_rxt_ndx( 'tag_TERPA3CO3_NO2' )
    usr_TERPA3PAN_M_ndx     = get_rxt_ndx( 'usr_TERPA3PAN_M' )
    usr_TERPNT_aer_ndx       = get_rxt_ndx( 'usr_TERPNT_aer' )
    usr_TERPNT1_aer_ndx      = get_rxt_ndx( 'usr_TERPNT1_aer' )
    usr_TERPNPT_aer_ndx      = get_rxt_ndx( 'usr_TERPNPT_aer' )
    usr_TERPNPT1_aer_ndx     = get_rxt_ndx( 'usr_TERPNPT1_aer' )
    usr_TERPFDN_aer_ndx        = get_rxt_ndx( 'usr_TERPFDN_aer' )
    usr_SQTN_aer_ndx        = get_rxt_ndx( 'usr_SQTN_aer' )
    usr_TERPHFN_aer_ndx        = get_rxt_ndx( 'usr_TERPHFN_aer' )
    usr_TERPDHDP_aer_ndx        = get_rxt_ndx( 'usr_TERPDHDP_aer' )
    usr_TERPACID_aer_ndx        = get_rxt_ndx( 'usr_TERPACID_aer' )
 !
 ! additional reactions for O3A/XNO
 !
    usr_OA_O2_ndx        = get_rxt_ndx( 'usr_OA_O2' )
    usr_XNO2NO3_M_ndx    = get_rxt_ndx( 'usr_XNO2NO3_M' )
    usr_NO2XNO3_M_ndx    = get_rxt_ndx( 'usr_NO2XNO3_M' )
    usr_XNO2NO3_aer_ndx  = get_rxt_ndx( 'usr_XNO2NO3_aer' )
    usr_NO2XNO3_aer_ndx  = get_rxt_ndx( 'usr_NO2XNO3_aer' )
    usr_XHNO3_OH_ndx     = get_rxt_ndx( 'usr_XHNO3_OH' )
    usr_XNO3_aer_ndx     = get_rxt_ndx( 'usr_XNO3_aer' )
    usr_XNO2_aer_ndx     = get_rxt_ndx( 'usr_XNO2_aer' )
    usr_MCO3_XNO2_ndx    = get_rxt_ndx( 'usr_MCO3_XNO2' )
    usr_XPAN_M_ndx       = get_rxt_ndx( 'usr_XPAN_M' )
    usr_XMPAN_M_ndx      = get_rxt_ndx( 'usr_XMPAN_M' )
    usr_XHO2NO2_M_ndx    = get_rxt_ndx( 'usr_XHO2NO2_M' )
!
! reduced hydrocarbon chemistry
!
    usr_C2O3_NO2_ndx     = get_rxt_ndx( 'usr_C2O3_NO2' )
    usr_C2H4_OH_ndx      = get_rxt_ndx( 'usr_C2H4_OH' )
    usr_XO2N_HO2_ndx     = get_rxt_ndx( 'usr_XO2N_HO2' )
    usr_C2O3_XNO2_ndx    = get_rxt_ndx( 'usr_C2O3_XNO2' )
!
    tag_XO2N_NO_ndx      = get_rxt_ndx( 'tag_XO2N_NO' )
    tag_XO2_HO2_ndx      = get_rxt_ndx( 'tag_XO2_HO2' )
    tag_XO2_NO_ndx       = get_rxt_ndx( 'tag_XO2_NO' )
!
! stratospheric chemistry
!
    usr_O_O_ndx          = get_rxt_ndx( 'usr_O_O' )
    usr_CL2O2_M_ndx      = get_rxt_ndx( 'usr_CL2O2_M' )
    usr_SO3_H2O_ndx      = get_rxt_ndx( 'usr_SO3_H2O' )
!
    tag_CLO_CLO_M_ndx      = get_rxt_ndx( 'tag_CLO_CLO_M' )
    if (tag_CLO_CLO_M_ndx<0) then ! for backwards compatibility
       tag_CLO_CLO_M_ndx   = get_rxt_ndx( 'tag_CLO_CLO' )
    endif
!
! reactions to remove BAM aerosols in the stratosphere
!
    usr_strat_tau_ndx( 1) = get_rxt_ndx( 'usr_CB1_strat_tau' )
    usr_strat_tau_ndx( 2) = get_rxt_ndx( 'usr_CB2_strat_tau' )
    usr_strat_tau_ndx( 3) = get_rxt_ndx( 'usr_OC1_strat_tau' )
    usr_strat_tau_ndx( 4) = get_rxt_ndx( 'usr_OC2_strat_tau' )
    usr_strat_tau_ndx( 5) = get_rxt_ndx( 'usr_SO4_strat_tau' )
    usr_strat_tau_ndx( 6) = get_rxt_ndx( 'usr_SOA_strat_tau' )
    usr_strat_tau_ndx( 7) = get_rxt_ndx( 'usr_NH4_strat_tau' )
    usr_strat_tau_ndx( 8) = get_rxt_ndx( 'usr_NH4NO3_strat_tau' )
    usr_strat_tau_ndx( 9) = get_rxt_ndx( 'usr_SSLT01_strat_tau' )
    usr_strat_tau_ndx(10) = get_rxt_ndx( 'usr_SSLT02_strat_tau' )
    usr_strat_tau_ndx(11) = get_rxt_ndx( 'usr_SSLT03_strat_tau' )
    usr_strat_tau_ndx(12) = get_rxt_ndx( 'usr_SSLT04_strat_tau' )
    usr_strat_tau_ndx(13) = get_rxt_ndx( 'usr_DST01_strat_tau' )
    usr_strat_tau_ndx(14) = get_rxt_ndx( 'usr_DST02_strat_tau' )
    usr_strat_tau_ndx(15) = get_rxt_ndx( 'usr_DST03_strat_tau' )
    usr_strat_tau_ndx(16) = get_rxt_ndx( 'usr_DST04_strat_tau' )
    usr_strat_tau_ndx(17) = get_rxt_ndx( 'usr_SO2t_strat_tau' )
    usr_strat_tau_ndx(18) = get_rxt_ndx( 'usr_SOAI_strat_tau' )
    usr_strat_tau_ndx(19) = get_rxt_ndx( 'usr_SOAM_strat_tau' )
    usr_strat_tau_ndx(20) = get_rxt_ndx( 'usr_SOAB_strat_tau' )
    usr_strat_tau_ndx(21) = get_rxt_ndx( 'usr_SOAT_strat_tau' )
    usr_strat_tau_ndx(22) = get_rxt_ndx( 'usr_SOAX_strat_tau' )
!
! stratospheric aerosol chemistry
!
    het1_ndx             = get_rxt_ndx( 'het1' )
!
! ion chemistry
!
    ion1_ndx  = get_rxt_ndx( 'ion_Op_O2' )
    ion2_ndx  = get_rxt_ndx( 'ion_Op_N2' )
    ion3_ndx  = get_rxt_ndx( 'ion_N2p_Oa' )
    ion11_ndx = get_rxt_ndx( 'ion_N2p_Ob' )

    elec1_ndx  = get_rxt_ndx( 'elec1' )
    elec2_ndx  = get_rxt_ndx( 'elec2' )
    elec3_ndx  = get_rxt_ndx( 'elec3' )

! ================================================
! VSLS Halogen recycling on seasalts
! orginally from ordc
! ================================================
    het_ss_0_ndx         = get_rxt_ndx( 'het_ss_0' )
    het_ss_1_ndx         = get_rxt_ndx( 'het_ss_1' )
    het_ss_2_ndx         = get_rxt_ndx( 'het_ss_2' )
    het_ss_3_ndx         = get_rxt_ndx( 'het_ss_3' )
    het_ss_4_ndx         = get_rxt_ndx( 'het_ss_4' )
    het_ss_5_ndx         = get_rxt_ndx( 'het_ss_5' )
    het_ss_6_ndx         = get_rxt_ndx( 'het_ss_6' )
    het_ss_7_ndx         = get_rxt_ndx( 'het_ss_7' )
    het_ss_8_ndx         = get_rxt_ndx( 'het_ss_8' )
    het_ss_9_ndx         = get_rxt_ndx( 'het_ss_9'  )
    het_ss_10_ndx        = get_rxt_ndx( 'het_ss_10' )
    het_ss_11_ndx        = get_rxt_ndx( 'het_ss_11' )
    het_ss_12_ndx        = get_rxt_ndx( 'het_ss_12' )

    ss_ixoy_2_ndx        = get_rxt_ndx( 'ss_ixoy_2' )
    ss_ixoy_3_ndx        = get_rxt_ndx( 'ss_ixoy_3' )
    ss_ixoy_4_ndx        = get_rxt_ndx( 'ss_ixoy_4' )

!   sslt1_ndx = get_spc_ndx( 'SSLT01' )
!   sslt2_ndx = get_spc_ndx( 'SSLT02' )
!   sslt3_ndx = get_spc_ndx( 'SSLT03' )
!   sslt4_ndx = get_spc_ndx( 'SSLT04' )
    sslt1_ndx = get_spc_ndx( 'ncl_a1' )
    sslt2_ndx = get_spc_ndx( 'ncl_a2' )
    sslt3_ndx = get_spc_ndx( 'ncl_a3' )
    sslt4_ndx = get_spc_ndx( 'ncl_a4' )

    ice_trp_cl_1_ndx     = get_rxt_ndx( 'ice_trp_cl_1' )
    ice_trp_br_1_ndx     = get_rxt_ndx( 'ice_trp_br_1' )
    ice_trp_i_1_ndx      = get_rxt_ndx( 'ice_trp_i_1'  )
    ice_trp_i_2_ndx      = get_rxt_ndx( 'ice_trp_i_2' )
    ice_trp_i_3_ndx      = get_rxt_ndx( 'ice_trp_i_3' )
    ice_trp_i_4_ndx      = get_rxt_ndx( 'ice_trp_i_4' )
    ice_trp_hbr_5_ndx    = get_rxt_ndx( 'ice_trp_hbr_5' )
    ice_trp_hbr_6_ndx    = get_rxt_ndx( 'ice_trp_hbr_6' )
    ice_trp_hcl_5_ndx    = get_rxt_ndx( 'ice_trp_hcl_5' )
    ice_trp_hcl_6_ndx    = get_rxt_ndx( 'ice_trp_hcl_6' )
    ice_trp_hi_5_ndx     = get_rxt_ndx( 'ice_trp_hi_5'  )
    ice_trp_hi_6_ndx     = get_rxt_ndx( 'ice_trp_hi_6'  )

    ice_fr_hoi_ndx       = get_rxt_ndx( 'ice_fr_hoi' )
    liq_fr_hoi_ndx       = get_rxt_ndx( 'liq_fr_hoi' )
    ice_fr_hi_ndx        = get_rxt_ndx( 'ice_fr_hi' )
    liq_fr_hi_ndx        = get_rxt_ndx( 'liq_fr_hi' )
    ice_fr_iono2_ndx     = get_rxt_ndx( 'ice_fr_iono2' )
    liq_fr_iono2_ndx     = get_rxt_ndx( 'liq_fr_iono2' )
    ice_fr_brono2_ndx    = get_rxt_ndx( 'ice_fr_brono2' )

    usr_IO_IO_a_ndx      = get_rxt_ndx ( 'usr_IO_IO_a' )
    usr_IO_IO_b_ndx      = get_rxt_ndx ( 'usr_IO_IO_b' )
    usr_IO_OIO_ndx       = get_rxt_ndx ( 'usr_IO_OIO ' )
    usr_OIO_OIO_ndx      = get_rxt_ndx ( 'usr_OIO_OIO' )
    usr_HOI_NO3_ndx      = get_rxt_ndx ( 'usr_HOI_NO3' )
    usr_I2O2_a_ndx       = get_rxt_ndx ( 'usr_I2O2_a' )
    usr_I2O2_b_ndx       = get_rxt_ndx ( 'usr_I2O2_b' )
    usr_I2O4_ndx         = get_rxt_ndx ( 'usr_I2O4' )
    usr_IONO2_ndx        = get_rxt_ndx ( 'usr_IONO2' )

    id_hocl              = get_spc_ndx( 'HOCL' )
    id_hcl               = get_spc_ndx( 'HCL' )
    id_hbr               = get_spc_ndx( 'HBR' )
    id_hi                = get_spc_ndx( 'HI' )
    id_hobr              = get_spc_ndx( 'HOBR' )
    id_hoi               = get_spc_ndx( 'HOI' )
    id_clono2            = get_spc_ndx( 'CLONO2' )
    id_brono2            = get_spc_ndx( 'BRONO2' )
    id_iono2             = get_spc_ndx( 'IONO2' )

    id_n2o5              = get_spc_ndx( 'N2O5' )

    has_het_ss_rxts  = het_ss_0_ndx  > 0  .or.  het_ss_1_ndx  > 0  .or.  het_ss_2_ndx  > 0 .or.    &
                       het_ss_3_ndx  > 0  .or.  het_ss_4_ndx  > 0  .or.  het_ss_5_ndx  > 0 .or.    &
                       het_ss_6_ndx  > 0  .or.  het_ss_7_ndx  > 0  .or.  het_ss_8_ndx  > 0 .or.    &
                       het_ss_9_ndx  > 0  .or.  het_ss_10_ndx > 0  .or.  het_ss_11_ndx > 0 .or. het_ss_12_ndx > 0

    has_ss_ixoy_rxts = ss_ixoy_2_ndx  > 0  .or.  ss_ixoy_3_ndx  > 0  .or.  ss_ixoy_4_ndx  > 0

    has_ice_trp_rxts = ice_trp_cl_1_ndx > 0  .or.  ice_trp_br_1_ndx > 0  .or.  ice_trp_i_1_ndx  > 0  .or. &
                       ice_trp_i_2_ndx  > 0  .or.  ice_trp_i_3_ndx  > 0  .or.  ice_trp_i_4_ndx  > 0  .or. &
                       ice_trp_hbr_5_ndx > 0 .or.  ice_trp_hbr_6_ndx > 0 .or.  &
                       ice_trp_hcl_5_ndx > 0 .or.  ice_trp_hcl_6_ndx > 0 .or.  &
                       ice_trp_hi_5_ndx > 0  .or.  ice_trp_hi_6_ndx > 0  .or.  &
                       ice_fr_hoi_ndx > 0    .or.  liq_fr_hoi_ndx > 0    .or.  &
                       ice_fr_hi_ndx > 0     .or.  liq_fr_hi_ndx > 0     .or.  &
                       ice_fr_iono2_ndx > 0  .or.  liq_fr_iono2_ndx > 0  .or.  &
                       ice_fr_brono2_ndx > 0

    do i = 1,nean
      write (xchar,'(i4)') i
      rxtname = 'ean'//trim(adjustl(xchar))
      ean_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nrpe
      write (xchar,'(i4)') i
      rxtname = 'rpe'//trim(adjustl(xchar))
      rpe_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,npir
      write (xchar,'(i4)') i
      rxtname = 'pir'//trim(adjustl(xchar))
      pir_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nedn
      write (xchar,'(i4)') i
      rxtname = 'edn'//trim(adjustl(xchar))
      edn_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,nnir
      write (xchar,'(i4)') i
      rxtname = 'nir'//trim(adjustl(xchar))
      nir_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,niira
      write (xchar,'(i4)') i
      rxtname = 'iira'//trim(adjustl(xchar))
      iira_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    do i = 1,niirb
      write (xchar,'(i4)') i
      rxtname = 'iirb'//trim(adjustl(xchar))
      iirb_ndx(i) = get_rxt_ndx(trim(rxtname))
    enddo

    usr_clm_h2o_m_ndx = get_rxt_ndx( 'usr_CLm_H2O_M' )
    usr_clm_hcl_m_ndx = get_rxt_ndx( 'usr_CLm_HCL_M' )

    elec4_ndx  = get_rxt_ndx( 'Op2P_ea' )
    elec5_ndx  = get_rxt_ndx( 'Op2P_eb' )
    elec6_ndx  = get_rxt_ndx( 'Op2D_e' )

    has_ion_rxts = ion1_ndx>0 .and. ion2_ndx>0 .and. ion3_ndx>0 .and. elec1_ndx>0 &
                 .and. elec2_ndx>0 .and. elec3_ndx>0

    has_d_chem = &
         all(ean_ndx>0) .and. &
         all(rpe_ndx>0) .and. &
         all(pir_ndx>0) .and. &
         all(edn_ndx>0) .and. &
         all(nir_ndx>0) .and. &
         all(iira_ndx>0) .and. &
         all(iirb_ndx>0) .and. &
         usr_clm_h2o_m_ndx>0 .and. usr_clm_hcl_m_ndx>0

    h2o_ndx    = get_spc_ndx( 'H2O' )

    !
    ! llnl super fast
    !
    usr_oh_co_ndx  = get_rxt_ndx( 'usr_oh_co' )
    het_no2_h2o_ndx  = get_rxt_ndx( 'het_no2_h2o' )
    usr_oh_dms_ndx  = get_rxt_ndx( 'usr_oh_dms' )
    aq_so2_h2o2_ndx  = get_rxt_ndx( 'aq_so2_h2o2' )
    aq_so2_o3_ndx  = get_rxt_ndx( 'aq_so2_o3' )

!lke++
! CO tags
!
    usr_COhc_OH_ndx      = get_rxt_ndx( 'usr_COhc_OH' )
    usr_COme_OH_ndx      = get_rxt_ndx( 'usr_COme_OH' )
    usr_CO01_OH_ndx      = get_rxt_ndx( 'usr_CO01_OH' )
    usr_CO02_OH_ndx      = get_rxt_ndx( 'usr_CO02_OH' )
    usr_CO03_OH_ndx      = get_rxt_ndx( 'usr_CO03_OH' )
    usr_CO04_OH_ndx      = get_rxt_ndx( 'usr_CO04_OH' )
    usr_CO05_OH_ndx      = get_rxt_ndx( 'usr_CO05_OH' )
    usr_CO06_OH_ndx      = get_rxt_ndx( 'usr_CO06_OH' )
    usr_CO07_OH_ndx      = get_rxt_ndx( 'usr_CO07_OH' )
    usr_CO08_OH_ndx      = get_rxt_ndx( 'usr_CO08_OH' )
    usr_CO09_OH_ndx      = get_rxt_ndx( 'usr_CO09_OH' )
    usr_CO10_OH_ndx      = get_rxt_ndx( 'usr_CO10_OH' )
    usr_CO11_OH_ndx      = get_rxt_ndx( 'usr_CO11_OH' )
    usr_CO12_OH_ndx      = get_rxt_ndx( 'usr_CO12_OH' )
    usr_CO13_OH_ndx      = get_rxt_ndx( 'usr_CO13_OH' )
    usr_CO14_OH_ndx      = get_rxt_ndx( 'usr_CO14_OH' )
    usr_CO15_OH_ndx      = get_rxt_ndx( 'usr_CO15_OH' )
    usr_CO16_OH_ndx      = get_rxt_ndx( 'usr_CO16_OH' )
    usr_CO17_OH_ndx      = get_rxt_ndx( 'usr_CO17_OH' )
    usr_CO18_OH_ndx      = get_rxt_ndx( 'usr_CO18_OH' )
    usr_CO19_OH_ndx      = get_rxt_ndx( 'usr_CO19_OH' )
    usr_CO20_OH_ndx      = get_rxt_ndx( 'usr_CO20_OH' )
    usr_CO21_OH_ndx      = get_rxt_ndx( 'usr_CO21_OH' )
    usr_CO22_OH_ndx      = get_rxt_ndx( 'usr_CO22_OH' )
    usr_CO23_OH_ndx      = get_rxt_ndx( 'usr_CO23_OH' )
    usr_CO24_OH_ndx      = get_rxt_ndx( 'usr_CO24_OH' )
    usr_CO25_OH_ndx      = get_rxt_ndx( 'usr_CO25_OH' )
    usr_CO26_OH_ndx      = get_rxt_ndx( 'usr_CO26_OH' )
    usr_CO27_OH_ndx      = get_rxt_ndx( 'usr_CO27_OH' )
    usr_CO28_OH_ndx      = get_rxt_ndx( 'usr_CO28_OH' )
    usr_CO29_OH_ndx      = get_rxt_ndx( 'usr_CO29_OH' )
    usr_CO30_OH_ndx      = get_rxt_ndx( 'usr_CO30_OH' )
    usr_CO31_OH_ndx      = get_rxt_ndx( 'usr_CO31_OH' )
    usr_CO32_OH_ndx      = get_rxt_ndx( 'usr_CO32_OH' )
    usr_CO33_OH_ndx      = get_rxt_ndx( 'usr_CO33_OH' )
    usr_CO34_OH_ndx      = get_rxt_ndx( 'usr_CO34_OH' )
    usr_CO35_OH_ndx      = get_rxt_ndx( 'usr_CO35_OH' )
    usr_CO36_OH_ndx      = get_rxt_ndx( 'usr_CO36_OH' )
    usr_CO37_OH_ndx      = get_rxt_ndx( 'usr_CO37_OH' )
    usr_CO38_OH_ndx      = get_rxt_ndx( 'usr_CO38_OH' )
    usr_CO39_OH_ndx      = get_rxt_ndx( 'usr_CO39_OH' )
    usr_CO40_OH_ndx      = get_rxt_ndx( 'usr_CO40_OH' )
    usr_CO41_OH_ndx      = get_rxt_ndx( 'usr_CO41_OH' )
    usr_CO42_OH_ndx      = get_rxt_ndx( 'usr_CO42_OH' )
!lke--

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'usrrxt_inti: diagnostics '
       write(iulog,'(10i5)') usr_O_O2_ndx,usr_HO2_HO2_ndx,tag_NO2_NO3_ndx,usr_N2O5_M_ndx,tag_NO2_OH_ndx,usr_HNO3_OH_ndx &
                            ,tag_NO2_HO2_ndx,usr_HO2NO2_M_ndx,usr_N2O5_aer_ndx,usr_NO3_aer_ndx,usr_NO2_aer_ndx &
                            ,usr_CO_OH_b_ndx,tag_C2H4_OH_ndx,tag_C3H6_OH_ndx,tag_CH3CO3_NO2_ndx,usr_PAN_M_ndx,usr_CH3COCH3_OH_ndx &
                            ,usr_MCO3_NO2_ndx,usr_MPAN_M_ndx,usr_XOOH_OH_ndx,usr_SO2_OH_ndx,usr_N2D_O2_ndx,usr_N2D_e_ndx,usr_N_O2_ndx,usr_DMS_OH_ndx,usr_HO2_aer_ndx &
                            ,usr_GLYOXAL_aer_ndx,usr_ISOPNITA_aer_ndx,usr_ISOPNITB_aer_ndx,usr_ONITR_aer_ndx,usr_HONITR_aer_ndx &
                            ,usr_TERPNIT_aer_ndx,usr_NTERPOOH_aer_ndx,usr_NC4CHO_aer_ndx,usr_NC4CH2OH_aer_ndx,usr_ISOPZD1O2_ndx &
                            ,usr_ISOPZD4O2_ndx,usr_ISOPFDN_aer_ndx,usr_ISOPFNP_aer_ndx,usr_ISOPN2B_aer_ndx,usr_ISOPN1D_aer_ndx &
                            ,usr_ISOPN4D_aer_ndx,usr_INOOHD_aer_ndx,usr_INHEB_aer_ndx,usr_INHED_aer_ndx,usr_MACRN_aer_ndx &
                            ,usr_ISOPHFP_aer_ndx,usr_IEPOX_aer_ndx,usr_DHPMPAL_aer_ndx,usr_ISOPB1O2_NOn_ndx,usr_ISOPB1O2_NOa_ndx &
                            ,usr_ISOPB4O2_NOn_ndx,usr_ISOPB4O2_NOa_ndx,usr_ISOPED1O2_NOn_ndx,usr_ISOPED1O2_NOa_ndx &
                            ,usr_ISOPED4O2_NOn_ndx,usr_ISOPED4O2_NOa_ndx,usr_ISOPZD1O2_NOn_ndx,usr_ISOPZD1O2_NOa_ndx &
                            ,usr_ISOPZD4O2_NOn_ndx,usr_ISOPZD4O2_NOa_ndx,usr_ISOPNO3_NOn_ndx,usr_ISOPNO3_NOa_ndx &
                            ,usr_MVKO2_NOn_ndx,usr_MVKO2_NOa_ndx,usr_MACRO2_NOn_ndx,usr_MACRO2_NOa_ndx &
                            ,usr_IEPOXOO_NOn_ndx,usr_IEPOXOO_NOa_ndx,usr_ISOPN1DO2_NOn_ndx,usr_ISOPN1DO2_NOa_ndx &
                            ,usr_ISOPN2BO2_NOn_ndx,usr_ISOPN2BO2_NOa_ndx,usr_ISOPN3BO2_NOn_ndx,usr_ISOPN3BO2_NOa_ndx &
                            ,usr_ISOPN4DO2_NOn_ndx,usr_ISOPN4DO2_NOa_ndx,usr_ISOPNBNO3O2_NOn_ndx,usr_ISOPNBNO3O2_NOa_ndx &
                            ,usr_ISOPNOOHBO2_NOn_ndx,usr_ISOPNOOHBO2_NOa_ndx,usr_ISOPNOOHDO2_NOn_ndx,usr_ISOPNOOHDO2_NOa_ndx &
                            ,usr_NC4CHOO2_NOn_ndx,usr_NC4CHOO2_NOa_ndx,tag_MCO3_NO2_ndx &
                            ,usr_TERPNT_aer_ndx,tag_TERPA2CO3_NO2_ndx,usr_TERPA2PAN_M_ndx &
                            ,usr_TERPNT1_aer_ndx,usr_TERPNPT_aer_ndx,usr_TERPNPT1_aer_ndx,usr_TERPFDN_aer_ndx,usr_SQTN_aer_ndx &
                            ,usr_TERPHFN_aer_ndx,usr_TERPDHDP_aer_ndx,usr_TERPACID_aer_ndx,tag_TERPACO3_NO2_ndx &
                            ,usr_TERPAPAN_M_ndx,tag_TERPA3CO3_NO2_ndx, usr_TERPA3PAN_M_ndx,usr_ICHE_aer_ndx,usr_ISOPFNC_aer_ndx &
                            ,usr_ISOPFDNC_aer_ndx                                                     &
                            ,usr_IO_IO_a_ndx, usr_IO_IO_b_ndx, usr_IO_OIO_ndx, usr_OIO_OIO_ndx        &
                            ,usr_HOI_NO3_ndx                                                          &
                            ,usr_I2O2_a_ndx,  usr_I2O2_b_ndx,  usr_I2O4_ndx,   usr_IONO2_ndx          &
                            ,id_hocl,id_hcl,id_hbr,id_hi,id_hobr,id_hoi,id_clono2,id_brono2,id_iono2  &
                            ,usr_N2O5_HCL_ndx                                                         &
                            ,id_n2o5
       write(iulog,*) 'usrrxt_inti: has_het_ss_rxts : ',has_het_ss_rxts
       write(iulog,*) 'usrrxt_inti: has_ice_trp_rxts : ',has_ice_trp_rxts
       write(iulog,*) 'usrrxt_inti: has_ss_ixoy_rxts : ',has_ss_ixoy_rxts

    end if

  end subroutine usrrxt_inti

  subroutine usrrxt( state, rxt, temp, tempi, tempe, invariants, h2ovmr,                   &
                     pmid, m, sulfate, mmr, relhum, strato_sad,                            &
                     tropchemlev, dlat, ncol, sad_trop, reff_trop, cwat, mbar, pbuf,       &
                     sad_sslt, sad_sslt_eff,                                               &
                     clno2_yield,                                                          &
                     vmr, ocnfrac, icefrac,                                         &
                     sad_ice_trop, sad_liq_trop, sad_ice_trop_orig )

!-----------------------------------------------------------------
!        ... set the user specified reaction rates
!-----------------------------------------------------------------

    use mo_constants,    only : pi, avo => avogadro, boltz_cgs, rgas
    use chem_mods,       only : nfs, rxntot, gas_pcnst, inv_m_ndx=>indexm
    use mo_setinv,       only : inv_o2_ndx=>o2_ndx, inv_h2o_ndx=>h2o_ndx
    use physics_buffer,  only : physics_buffer_desc
    use physics_types,   only : physics_state
    use carma_flags_mod, only : carma_hetchem_feedback
    use aero_model,      only : aero_model_surfarea, n_supplemental_sad
    use radiative_aerosol,only : rad_aer_get_info
    use time_manager,     only : get_curr_calday
    use infnan,           only : nan
    use mo_slh_routines,  only : SSAdehal_ScalingFactor, SSAhno3_ScalingFactor, SSAn2o5_ScalingFactor, &
                                 ICEfraprx_ScalingFactor_I, ICEfraprx_ScalingFactor_Br, LIQfraprx_ScalingFactor_I

    implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    integer, intent(in)     :: tropchemlev(pcols)         ! trop/strat reaction separation vertical index
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    real(r8), intent(in)    :: temp(pcols,pver)           ! temperature (K); neutral temperature
    real(r8), intent(in)    :: tempi(pcols,pver)          ! ionic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: tempe(pcols,pver)          ! electronic temperature (K); only used if ion chemistry
    real(r8), intent(in)    :: m(ncol,pver)               ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(ncol,pver)         ! sulfate aerosol (mol/mol)
    real(r8), intent(in)    :: strato_sad(pcols,pver)     ! stratospheric aerosol sad (1/cm)
    real(r8), intent(in)    :: h2ovmr(ncol,pver)          ! water vapor (mol/mol)
    real(r8), intent(in)    :: relhum(ncol,pver)          ! relative humidity
    real(r8), intent(in)    :: pmid(pcols,pver)           ! midpoint pressure (Pa)
    real(r8), intent(in)    :: invariants(ncol,pver,nfs)  ! invariants density (/cm^3)
    real(r8), intent(in)    :: mmr(pcols,pver,gas_pcnst)  ! species concentrations (kg/kg)
    real(r8), intent(in)    :: cwat(ncol,pver) !PJC Condensed Water (liquid+ice) (kg/kg)
    real(r8), intent(in)    :: mbar(ncol,pver) !PJC Molar mass of air (g/mol)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)      ! gas phase rates
    real(r8), intent(out)   :: sad_trop(pcols,pver)       ! tropospheric surface area density (cm2/cm3)
    real(r8), intent(out)   :: reff_trop(pcols,pver)      ! tropospheric effective radius (cm)
    type(physics_state),    intent(in) :: state           ! Physics state variables

  ! --------------------------------
  ! VSLS Halogen chemistry incorporated
  ! --------------------------------
    real(r8), intent(in)    :: vmr(ncol ,pver,gas_pcnst)  ! species mixing ratio (mol/mol)
    real(r8), intent(in)    :: ocnfrac(pcols)             ! ocean fraction
    real(r8), intent(in)    :: icefrac(pcols)             ! sea-ice fraction
    real(r8), intent(in)    :: sad_ice_trop(ncol,pver)    ! surf area density of ice ( cm^2/cm^3 ) in Troposphere
    real(r8), intent(in)    :: sad_ice_trop_orig(ncol,pver)    ! surf area density of ice ( cm^2/cm^3 ) in Troposphere
    real(r8), intent(in)    :: sad_liq_trop(ncol,pver)    ! surf area density of liquid water ( cm^2/cm^3 ) in Troposphere
    real(r8), intent(inout) :: sad_sslt    (ncol,pver)    ! surf area density of ice ( cm^2/cm^3 ) in Troposphere
    real(r8), intent(inout) :: sad_sslt_eff(ncol,pver)    ! surf area density of ice ( cm^2/cm^3 ) in Troposphere
    real(r8), intent(inout) :: clno2_yield (ncol,pver)    ! originally implemented by jfl

    type(physics_buffer_desc), pointer :: pbuf(:)

!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------

    real(r8), parameter :: dg = 0.1_r8            ! mole diffusion =0.1 cm2/s (Dentener, 1993)

!-----------------------------------------------------------------
! 	... reaction probabilities for heterogeneous reactions
!-----------------------------------------------------------------
    real(r8), parameter :: gamma_n2o5 = 0.02_r8         ! JPL19
    real(r8), parameter :: gamma_ho2  = 0.10_r8         ! Gaubert et al., https://doi.org/10.5194/acp-20-14617-2020
    real(r8), parameter :: gamma_no2  = 8.0e-6_r8       ! Liu et al., Environ.Sci.&Tech, 53, 3517, 2019 doi:10.1021/acs.est.8b06367
    real(r8), parameter :: gamma_no3  = 0.002_r8        ! JPL19
    real(r8), parameter :: gamma_glyoxal  = 2.0e-4_r8   !  Washenfelder et al, JGR, 2011
!TS1 species
    real(r8), parameter :: gamma_isopnita  = 0.005_r8        ! from Fisher et al., ACP, 2016
    real(r8), parameter :: gamma_isopnitb  = 0.005_r8        !
    real(r8), parameter :: gamma_onitr     = 0.005_r8        !
    real(r8), parameter :: gamma_honitr    = 0.005_r8        !
    real(r8), parameter :: gamma_terpnit   = 0.01_r8         !
    real(r8), parameter :: gamma_nterpooh  = 0.01_r8         !
    real(r8), parameter :: gamma_nc4cho    = 0.02_r8        !
    real(r8), parameter :: gamma_nc4ch2oh  = 0.005_r8        !
!TS2 species
    real(r8), parameter :: gamma_isopfdn  = 0.1_r8          ! Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_isopfnp  = 0.1_r8          ! Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_isopn2b  = 0.02_r8        ! All isoprene nitrates Wolfe
    real(r8), parameter :: gamma_isopn1d  = 0.02_r8        !
    real(r8), parameter :: gamma_isopn4d  = 0.02_r8        !
    real(r8), parameter :: gamma_inoohd   = 0.02_r8        !
    real(r8), parameter :: gamma_inheb    = 0.02_r8       !Marais 2015 for IEPOX
    real(r8), parameter :: gamma_inhed    = 0.02_r8       !Marais 2015 for IEPOX
    real(r8), parameter :: gamma_macrn    = 0.02_r8        !
    real(r8), parameter :: gamma_isophfp  = 0.1_r8          !Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_iepox    = 0.0042_r8       !Marais 2015 for IEPOX
    real(r8), parameter :: gamma_dhpmpal  = 0.1_r8          !Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_iche     = 0.0042_r8       !Marais 2015 for IEPOX
    real(r8), parameter :: gamma_isopfnc  = 0.1_r8          ! Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_isopfdnc = 0.1_r8          ! Marais 2015 for C5-LVOC
    real(r8), parameter :: gamma_terpnt   = 0.02_r8        !
    real(r8), parameter :: gamma_terpnt1  = 0.02_r8        !
    real(r8), parameter :: gamma_terpnpt  = 0.02_r8        !
    real(r8), parameter :: gamma_terpnpt1 = 0.02_r8        !
    real(r8), parameter :: gamma_terpfdn  = 0.1_r8        !
    real(r8), parameter :: gamma_sqtn     = 0.1_r8        !
    real(r8), parameter :: gamma_terphfn  = 0.1_r8        !
    real(r8), parameter :: gamma_terpdhdp = 0.1_r8        !
    real(r8), parameter :: gamma_terpacid = 0.01_r8        !

    integer  ::  i, k
    integer  ::  l
    real(r8) ::  tp(ncol)                       ! 300/t
    real(r8) ::  tinv(ncol)                     ! 1/t
    real(r8) ::  ko(ncol)
    real(r8) ::  term1(ncol)
    real(r8) ::  term2(ncol)
    real(r8) ::  kinf(ncol)
    real(r8) ::  fc(ncol)
    real(r8) ::  xr(ncol)
    real(r8) ::  sur(ncol)
    real(r8) ::  sqrt_t(ncol)                   ! sqrt( temp )
    real(r8) ::  sqrt_t_58(ncol)                ! sqrt( temp / 58.)
    real(r8) ::  exp_fac(ncol)                  ! vector exponential
    real(r8) ::  lwc(ncol)
    real(r8) ::  ko_m(ncol)
    real(r8) ::  k0(ncol)
    real(r8) ::  kinf_m(ncol)
    real(r8) ::  o2(ncol)
    real(r8) ::  c_n2o5, c_ho2, c_no2, c_no3, c_glyoxal
!TS1 species
    real(r8) ::  c_isopnita, c_isopnitb, c_onitr, c_honitr, c_terpnit, c_nterpooh
    real(r8) ::  c_nc4cho, c_nc4ch2oh
!T2 species
    real(r8) ::  c_isopfdn, c_isopfnp, c_isopn2b, c_isopn1d, c_isopn4d, c_inoohd
    real(r8) ::  c_inheb, c_inhed, c_macrn, c_isophfp, c_iepox, c_dhpmpal
    real(r8) ::  c_iche, c_isopfnc, c_isopfdnc
    real(r8) ::  c_terpnt, c_terpnt1, c_terpnpt, c_terpnpt1, c_terpfdn, c_sqtn, c_terphfn, c_terpdhdp, c_terpacid

  ! --------------------------------------
  ! VSLS Halogen chemistry incorporated
  ! -------------------------------------
    real(r8), parameter :: gamma_brono2_ss    = 0.01_r8    ! originally 0.08_r8 in CESM1 (ordc) ! originaly in CESM2 0.020_r8, now x0.5 (jv)
    real(r8), parameter :: gamma_brno2_ss     = 0.005_r8   ! originally 0.04_r8 in CESM1 (ordc) ! originaly in CESM2 0.010_r8, now x0.5 (jv)
    real(r8), parameter :: gamma_hobr_ss      = 0.0125_r8  ! originally 0.1_r8  in CESM1 (ordc) ! originaly in CESM2 0.025_r8, now x0.5 (jv)
    real(r8), parameter :: gamma_clono2_ss    = 0.006_r8   ! originally 0.02_r8 in CESM1 (ordc)
    real(r8), parameter :: gamma_clno2_ss     = 0.006_r8   ! originally 0.02_r8 in CESM1 (ordc)
    real(r8), parameter :: gamma_hocl_ss      = 0.03_r8    ! originally 0.1_r8  in CESM1 (ordc)
    real(r8), parameter :: gamma_iono2_ss     = 0.003_r8   ! originally 0.01_r8 in CESM1 (ordc)
    real(r8), parameter :: gamma_ino2_ss      = 0.006_r8   ! originally 0.02_r8 in CESM1 (ordc)
    real(r8), parameter :: gamma_hoi_ss       = 0.0018_r8  ! originally 0.06_r8 in CESM1 (ordc)

    real(r8), parameter :: gamma_ss_ixoy_2    = 0.0025_r8  ! originally 0.01_r8 in CESM1 (rpf)
    real(r8), parameter :: gamma_ss_ixoy_3    = 0.0025_r8  ! originally 0.01_r8 in CESM1 (rpf)
    real(r8), parameter :: gamma_ss_ixoy_4    = 0.0025_r8  ! originally 0.01_r8 in CESM1 (rpf)

    real(r8), parameter :: gamma_hno3_ss      = 0.010_r8   ! originally 0.05_r8 in CESM1 (liqy)
    real(r8), parameter :: gamma_n2o5_ss      = 0.0143_r8  ! identical to CESM1 (jfl)
    real(r8), parameter :: gamma_hcln2o5_aerosol = 0.0143_r8

    real(r8), parameter :: gamma_clno2_sul  = 0.1_r8
    real(r8), parameter :: gamma_brno2_sul  = 0.1_r8
    real(r8), parameter :: gamma_iono2_sul  = 0.1_r8
    real(r8), parameter :: gamma_ino2_sul   = 0.1_r8
    real(r8), parameter :: gamma_brono2_sul = 0.2_r8
    real(r8), parameter :: gamma_hbr_wet    = 0.1_r8
    real(r8), parameter :: gamma_hcl_wet    = 0.1_r8
    real(r8), parameter :: gamma_hi_wet     = 0.1_r8

    !rpf: Now we implemented ScalingFactors in user_nl_cam (&slh_nl)
    real(r8), parameter :: gamma_fr_hoi_ice   = 3.0e-4_r8 ! originally 3.0e-4 for 2x2.5(26L) and 1.6e-4 for 1x1(56L) in CESM1 - Use ICEfraprx_ScalingFactor_I (rpf)
    real(r8), parameter :: gamma_fr_hoi_liq   = 3.0e-4_r8 ! originally 1.0e-4_r8 in CESM1 (rpf) - up to 3.0e-3_r8 in CESM2 - Use LIQfraprx_ScalingFactor_I (rpf)
    real(r8), parameter :: gamma_fr_iono2_ice = 0.02_r8   ! originally 0.005_r8 in CESM1 (rpf)
    real(r8), parameter :: gamma_fr_iono2_liq = 0.03_r8   ! identical to CESM1 (rpf)
    real(r8), parameter :: gamma_fr_hi_ice = 0.02_r8      ! identical to CESM1 (rpf)
    real(r8), parameter :: gamma_fr_hi_liq = 0.02_r8      ! identical to CESM1 (rpf)

    real(r8), parameter :: gamma_fr_brono2_ice = 1.0e-2_r8  ! in CESM1 BrONO2 washout was based on NEU routine (rpf) - from 1.0e-2_r8 up to 0.030_r8 in CESM2 - Use ICEfraprx_ScalingFactor_Br (rpf)

    real(r8), parameter :: gamma_hclhocl_ice = 0.2_r8     ! copied values from stratosphere
    real(r8), parameter :: gamma_hclhobr_ice = 0.3_r8     ! copied values from stratosphere
    real(r8), parameter :: gamma_hbrhocl_ice = 0.2_r8
    real(r8), parameter :: gamma_hbrhobr_ice = 0.12_r8
    real(r8), parameter :: gamma_hihocl_ice  = 0.12_r8
    real(r8), parameter :: gamma_hihobr_ice  = 0.12_r8
    real(r8), parameter :: gamma_hbrhoi_ice  = 0.12_r8
    real(r8), parameter :: gamma_hihoi_ice   = 0.12_r8
    real(r8), parameter :: gamma_hclhoi_ice  = 0.12_r8
    real(r8), parameter :: gamma_clono2_ice  = 0.3_r8
    real(r8), parameter :: gamma_brono2_ice  = 0.1_r8
    real(r8), parameter :: gamma_iono2_ice   = 0.1_r8
    real(r8), parameter :: gamma_hoi_ice     = 0.3_r8
    real(r8), parameter :: gamma_hbrhocl_add = 0.25_r8
    real(r8), parameter :: gamma_hbrhobr_add = 0.25_r8
    real(r8), parameter :: gamma_clono2_add  = 0.3_r8
    real(r8), parameter :: gamma_brono2_add  = 0.3_r8

    real(r8)            :: hclvmr, hcldeni
    real(r8)            :: hoclvmr, hocldeni
    real(r8)            :: hbrvmr, hbrdeni
    real(r8)            :: hobrvmr, hobrdeni
    real(r8)            :: clono2vmr, brono2vmr, iono2vmr
    real(r8)            :: hivmr, hoivmr, hoideni, hideni
    real(r8)            :: n2o5vmr, n2o5deni

    real(r8), parameter :: small = 1.e-16_r8

    real(r8)            :: latitude, DF     ! DF is the Depletion Factor (Yang et al., JGR 2005)
    real(r8), parameter :: rad2deg = 180.0_r8/pi  ! radians to degrees
    real(r8), parameter :: dfmax = 0.9_r8   ! originally 0.7_r8 in CESM1 (rpf)
    real(r8), parameter :: dfmin = 0.3_r8   ! originally 0.1_r8 in CESM1 (rpf)
    real(r8)            :: calday
    real(r8)            :: press_lev            ! Mid-point pressure (hPa)
    real(r8)            :: logical_sslt         ! Logical condition to avoid SSLT_recycling in the Stratosphere
    real(r8)            :: sad_sslt_mask        ! Logical condition to avoid SSLT_recycling above continents
    real(r8)            :: ssur_tot
    real(r8) ::  press_1d(ncol)                 ! Mid-point pressure (hPa)
    real(r8) ::  temp_1d (ncol)                 ! and temperature (K).
    real(r8) ::  work1 (ncol), work2 (ncol)     ! A & pre-exponential factors.
                                                ! Used for some user-defined iodine reacs.
                                                ! (some of them based on RRKM theory)
    real(r8) ::  amas
    !-----------------------------------------------------------------
    !	... density of sulfate aerosol
    !-----------------------------------------------------------------
    real(r8), parameter :: gam1 = 0.04_r8                 ! N2O5+SUL ->2HNO3
    real(r8), parameter :: wso4 = 98._r8
    real(r8), parameter :: den  = 1.15_r8                 ! each molecule of SO4(aer) density g/cm3
    !-------------------------------------------------
    ! 	... volume of sulfate particles
    !           assuming mean rm
    !           continient 0.05um  0.07um  0.09um
    !           ocean      0.09um  0.25um  0.37um
    !                      0.16um                  Blake JGR,7195, 1995
    !-------------------------------------------------
    real(r8), parameter :: rm1  = 0.16_r8*1.e-4_r8             ! mean radii in cm
    real(r8), parameter :: fare = 4._r8*pi*rm1*rm1             ! each mean particle(r=0.1u) area   cm2/cm3

    !-----------------------------------------------------------------------
    !        ... Aqueous phase sulfur quantities for SO2 + H2O2 and SO2 + O3
    !-----------------------------------------------------------------------
    real(r8), parameter  :: HENRY298_H2O2 =  7.45e+04_r8
    real(r8), parameter  :: H298_H2O2     = -1.45e+04_r8
    real(r8), parameter  :: HENRY298_SO2  =  1.23e+00_r8
    real(r8), parameter  :: H298_SO2      = -6.25e+03_r8
    real(r8), parameter  :: K298_SO2_HSO3 =  1.3e-02_r8
    real(r8), parameter  :: H298_SO2_HSO3 = -4.16e+03_r8
    real(r8), parameter  :: R_CONC        =  82.05e+00_r8 / avo
    real(r8), parameter  :: R_CAL         =  rgas * 0.239006e+00_r8
    real(r8), parameter  :: K_AQ          =  7.57e+07_r8
    real(r8), parameter  :: ER_AQ         =  4.43e+03_r8

    real(r8), parameter  :: HENRY298_O3   =  1.13e-02_r8
    real(r8), parameter  :: H298_O3       = -5.04e+03_r8
    real(r8), parameter  :: K298_HSO3_SO3 =  6.6e-08_r8
    real(r8), parameter  :: H298_HSO3_SO3 = -2.23e+03_r8
    real(r8), parameter  :: K0_AQ         =  2.4e+04_r8
    real(r8), parameter  :: ER0_AQ        =  0.0e+00_r8
    real(r8), parameter  :: K1_AQ         =  3.7e+05_r8
    real(r8), parameter  :: ER1_AQ        =  5.53e+03_r8
    real(r8), parameter  :: K2_AQ         =  1.5e+09_r8
    real(r8), parameter  :: ER2_AQ        =  5.28e+03_r8

    real(r8), parameter  :: pH            =  4.5e+00_r8

    real(r8), pointer :: sfc(:), dm_aer(:)
    integer :: ntot_amode, nbins, naero

    real(r8), pointer :: sfc_array(:,:,:), dm_array(:,:,:)
 !TS2
    real(r8) ::  aterm(ncol)
    real(r8) ::  natom
    real(r8) ::  nyield
    real(r8) ::  acorr
    real(r8) ::  exp_natom
    character(len=*), parameter :: subname = 'usrrxt'

    real(r8) :: total_sslt_mass
    real(r8) :: xmin,xmax,ymin,ymax,clno2_a,clno2_b
    real(r8) :: clno2_yield_Erin
    real(r8) :: dummy_n2o5_clno2_rate

    ! jfl: parameters for clno2 yield
    xmin = 5.e-10_r8
    xmax = 1.e-9_r8
    ymin = 0.0_r8
    ymax = 1.0_r8
    clno2_a = (ymax-ymin)/(xmax-xmin)
    clno2_b = ymax - clno2_a * xmax
    clno2_yield = 0._r8
    clno2_yield_Erin = 0.138_r8
    dummy_n2o5_clno2_rate = 0._r8

    ! get info about the modal aerosols
    ! get ntot_amode
    call rad_aer_get_info(0, nmodes=ntot_amode)
    call rad_aer_get_info(0, nbins=nbins)

    if (ntot_amode>0.and.nbins>0) then
      call endrun(subname // ':: ERROR running with MAM and CARMA simultaneously not supported.')
    end if

    calday = get_curr_calday()

    ! one slot per aerosol mode/bin, plus the slots for supplemental SADs from species not in rad_climate
    ! (bulk only: offline sulfate, ammonium nitrate, secondary organics)
    if (ntot_amode>0) then
       allocate(sfc_array(pcols,pver,ntot_amode+n_supplemental_sad), &
                dm_array (pcols,pver,ntot_amode+n_supplemental_sad) )
    else if (nbins>0) then
       allocate(sfc_array(pcols,pver,nbins+n_supplemental_sad), &
                dm_array (pcols,pver,nbins+n_supplemental_sad) )
    else
       call rad_aer_get_info(0, naero=naero)
       allocate(sfc_array(pcols,pver,naero+n_supplemental_sad), &
                dm_array (pcols,pver,naero+n_supplemental_sad) )
    endif

    sfc_array(:,:,:) = 0._r8
    dm_array (:,:,:) = 0._r8
    sad_trop (:,:)   = 0._r8
    reff_trop(:,:)   = 0._r8
    sad_sslt     (:,:) = 0._r8
    sad_sslt_eff (:,:) = 0._r8

    if( usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0 ) then

! sad_trop should be set outside of usrrxt ??
       if( carma_hetchem_feedback ) then
          sad_trop(:ncol,:pver)=strato_sad(:ncol,:pver)
       else
          call aero_model_surfarea( &
               state, relhum, pmid, temp, tropchemlev, &
               sfc_array, dm_array, sad_trop, reff_trop, sad_sslt )
       endif
    endif

    hoclvmr   = 0._r8
    hocldeni  = 0._r8
    hbrvmr    = 0._r8
    hbrdeni   = 0._r8
    hobrvmr   = 0._r8
    hobrdeni  = 0._r8
    clono2vmr = 0._r8
    brono2vmr = 0._r8
    iono2vmr  = 0._r8
    hoivmr    = 0._r8
    hoideni   = 0._r8
    !-----------------------------------------------------------------------
    !     	... intialize rate constants
    !-----------------------------------------------------------------------
    if ( ss_ixoy_2_ndx > 0 ) rxt(:,:,ss_ixoy_2_ndx) = 0._r8
    if ( ss_ixoy_3_ndx > 0 ) rxt(:,:,ss_ixoy_3_ndx) = 0._r8
    if ( ss_ixoy_4_ndx > 0 ) rxt(:,:,ss_ixoy_4_ndx) = 0._r8

    if ( ice_trp_i_2_ndx > 0 )   rxt(:,:,ice_trp_i_2_ndx)   = 0._r8
    if ( ice_trp_i_3_ndx > 0 )   rxt(:,:,ice_trp_i_3_ndx)   = 0._r8
    if ( ice_trp_i_4_ndx > 0 )   rxt(:,:,ice_trp_i_4_ndx)   = 0._r8

    if ( ice_trp_cl_1_ndx > 0 )  rxt(:,:,ice_trp_cl_1_ndx)  = 0._r8
    if ( ice_trp_br_1_ndx > 0 )  rxt(:,:,ice_trp_br_1_ndx)  = 0._r8
    if ( ice_trp_i_1_ndx > 0 )   rxt(:,:,ice_trp_i_1_ndx)   = 0._r8

    if ( ice_trp_hbr_5_ndx > 0 ) rxt(:,:,ice_trp_hbr_5_ndx) = 0._r8
    if ( ice_trp_hbr_6_ndx > 0 ) rxt(:,:,ice_trp_hbr_6_ndx) = 0._r8
    if ( ice_trp_hcl_5_ndx > 0 ) rxt(:,:,ice_trp_hcl_5_ndx) = 0._r8
    if ( ice_trp_hcl_6_ndx > 0 ) rxt(:,:,ice_trp_hcl_6_ndx) = 0._r8
    if ( ice_trp_hi_5_ndx > 0 )  rxt(:,:,ice_trp_hi_5_ndx)  = 0._r8
    if ( ice_trp_hi_6_ndx > 0 )  rxt(:,:,ice_trp_hi_6_ndx)  = 0._r8

    if ( ice_fr_hoi_ndx > 0 )    rxt(:,:,ice_fr_hoi_ndx)    = 0._r8
    if ( liq_fr_hoi_ndx > 0 )    rxt(:,:,liq_fr_hoi_ndx)    = 0._r8
    if ( ice_fr_hi_ndx > 0 )     rxt(:,:,ice_fr_hi_ndx)     = 0._r8
    if ( liq_fr_hi_ndx > 0 )     rxt(:,:,liq_fr_hi_ndx)     = 0._r8
    if ( ice_fr_iono2_ndx > 0 )  rxt(:,:,ice_fr_iono2_ndx)  = 0._r8
    if ( liq_fr_iono2_ndx > 0 )  rxt(:,:,liq_fr_iono2_ndx)  = 0._r8

    if ( ice_fr_brono2_ndx > 0 ) rxt(:,:,ice_fr_brono2_ndx) = 0._r8

    if (usr_N2O5_HCL_ndx > 0 )   rxt(:,:,usr_N2O5_HCL_ndx)  = 0._r8

    level_loop : do k = 1,pver
       tinv(:)           = 1._r8 / temp(:ncol,k)
       tp(:)             = 300._r8 * tinv(:)
       sqrt_t(:)         = sqrt( temp(:ncol,k) )
       sqrt_t_58(:)      = sqrt( temp(:ncol,k) / 58.0_r8 )
       press_1d(:)       = pmid (:ncol,k) / 100._r8  !hPa
       temp_1d (:)       = temp (:ncol,k)

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_O_O2_ndx > 0 ) then
          rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if
       if( usr_OA_O2_ndx > 0 ) then
          rxt(:,k,usr_OA_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
       end if

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
       if ( usr_O_O_ndx > 0 ) then
          rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )
       end if

!-----------------------------------------------------------------
! 	... cl2o2 + m -> 2*clo + m  (JPL15-10)
!-----------------------------------------------------------------
       if ( usr_CL2O2_M_ndx > 0 ) then
          if ( tag_CLO_CLO_M_ndx > 0 ) then
             ko(:)            = 2.16e-27_r8 * exp( 8537.0_r8* tinv(:) )
             rxt(:,k,usr_CL2O2_M_ndx) = rxt(:,k,tag_CLO_CLO_M_ndx)/ko(:)
          else
             rxt(:,k,usr_CL2O2_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!       ... so3 + 2*h2o --> h2so4 + h2o
!       Note: this reaction proceeds by the 2 intermediate steps below
!           so3 + h2o --> adduct
!           adduct + h2o --> h2so4 + h2o
!               (Lovejoy et al., JCP, pp. 19911-19916, 1996)
!	The first order rate constant used here is recommended by JPL 2011.
!	This rate involves the water vapor number density.
!-----------------------------------------------------------------

       if ( usr_SO3_H2O_ndx > 0 ) then
          call comp_exp( exp_fac, 6540.0_r8*tinv(:), ncol )
          if( h2o_ndx > 0 ) then
             fc(:) = 8.5e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 8.5e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_SO3_H2O_ndx) = 1.0e-20_r8 * fc(:)
       end if

!-----------------------------------------------------------------
!	... n2o5 + m --> no2 + no3 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_N2O5_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_N2O5_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 1.724138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_N2O5_M_ndx) = 0._r8
          end if
       end if
       if( usr_XNO2NO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_XNO2NO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) *1.724138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XNO2NO3_M_ndx) = 0._r8
          end if
       end if
       if( usr_NO2XNO3_M_ndx > 0 ) then
          if( tag_NO2_NO3_ndx > 0 ) then
             call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
             rxt(:,k,usr_NO2XNO3_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 1.724138e26_r8 * exp_fac(:)
          else
             rxt(:,k,usr_NO2XNO3_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	set rates for:
! 	... hno3 + oh --> no3 + h2o
!           ho2no2 + m --> ho2 + no2 + m
!-----------------------------------------------------------------
       if( usr_HNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_HNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_XHNO3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 1335._r8*tinv, ncol )
          ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
          call comp_exp( exp_fac, 2199._r8*tinv, ncol )
          ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          rxt(:,k,usr_XHNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
       end if
       if( usr_HO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_HO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_HO2NO2_M_ndx) = 0._r8
          end if
       end if
       if( usr_XHO2NO2_M_ndx > 0 ) then
          if( tag_NO2_HO2_ndx > 0 ) then
             call comp_exp( exp_fac, -10900._r8*tinv, ncol )
             rxt(:,k,usr_XHO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
          else
             rxt(:,k,usr_XHO2NO2_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
! 	... co + oh --> co2 + ho2 (new single reaction for combined branches [JPL19])
!-----------------------------------------------------------------
       if( usr_CO_OH_ndx > 0 ) then
         ko  (:)  = 6.9e-33_r8 * ( 298._r8 / temp(:ncol,k) )**(2.1_r8)
         kinf(:)  = 1.1e-12_r8 * ( 298._r8 / temp(:ncol,k) )**(-1.3_r8)

         term2(:) = (1 + (log10( ko(:)*m(:,k) / kinf(:) ))**2)**(-1)

         term1(:) = (kinf(:) * ko(:)*m(:,k)) / (kinf(:) + ko(:)*m(:,k)) * (0.6_r8)**term2(:)

         rxt(:ncol,k,usr_CO_OH_ndx) = term1(:) + 1.85e-13_r8 * exp(-65._r8/temp(:ncol,k)) * (1._r8 - term1(:)/kinf(:))

       end if
!-----------------------------------------------------------------
!       ... co + oh --> co2 + ho2     (combined branches - do not use with CO_OH_b)
!       note: for mechanisms prior to Dec 2022
!-----------------------------------------------------------------
       if( usr_CO_OH_a_ndx > 0 ) then
          rxt(:,k,usr_CO_OH_a_ndx) = 1.5e-13_r8 * &
               (1._r8 + 6.e-7_r8*boltz_cgs*m(:,k)*temp(:ncol,k))
       end if
!-----------------------------------------------------------------
! 	... co + oh --> co2 + h (second branch JPL15-10, with CO+OH+M)
!       note: for mechanisms prior to Dec 2022
!-----------------------------------------------------------------
       if( usr_CO_OH_b_ndx > 0 ) then
         kinf(:)  = 2.1e+09_r8 * (temp(:ncol,k)/ t0)**(6.1_r8)
         ko  (:)  = 1.5e-13_r8

         term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
         term2(:) = ko(:) / (1._r8 + term1(:))

         term1(:) = log10( term1(:) )
         term1(:) = 1.0_r8 / (1.0_r8 + term1(:)*term1(:))

         rxt(:ncol,k,usr_CO_OH_b_ndx) = term2(:) * (0.6_r8)**term1(:)
       end if

!-----------------------------------------------------------------
!       ... ho2 + ho2 --> h2o2
!       note: this rate involves the water vapor number density
!-----------------------------------------------------------------
       if( usr_HO2_HO2_ndx > 0 ) then

          call comp_exp( exp_fac, 460._r8*tinv, ncol )
          ko(:)   = 3.0e-13_r8 * exp_fac(:)
          call comp_exp( exp_fac, 920._r8*tinv, ncol )
          kinf(:) = 2.1e-33_r8 * m(:,k) * exp_fac(:)
          call comp_exp( exp_fac, 2200._r8*tinv, ncol )

          if( h2o_ndx > 0 ) then
             fc(:) = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
          else
             fc(:) = 1._r8 + 1.4e-21_r8 * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
          end if
          rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

       end if

!-----------------------------------------------------------------
!    	... mco3 + no2 -> mpan
!-----------------------------------------------------------------
       if( usr_MCO3_NO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_NO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if
       if( usr_MCO3_XNO2_ndx > 0 ) then
          rxt(:,k,usr_MCO3_XNO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
       end if

!-----------------------------------------------------------------
!	... pan + m --> ch3co3 + no2 + m (JPL15-10)
!-----------------------------------------------------------------
       call comp_exp( exp_fac, -14000._r8*tinv, ncol )
       if( usr_PAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_PAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_PAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XPAN_M_ndx > 0 ) then
          if( tag_CH3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XPAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XPAN_M_ndx) = 0._r8
          end if
       end if

!-----------------------------------------------------------------
!	... mpan + m --> mco3 + no2 + m (JPL15-10)
!-----------------------------------------------------------------
       if( usr_MPAN_M_ndx > 0 ) then
          if( tag_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_MPAN_M_ndx) = rxt(:,k,tag_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_MPAN_M_ndx) = 0._r8
          end if
       end if
       if( usr_XMPAN_M_ndx > 0 ) then
          if( tag_MCO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_XMPAN_M_ndx) = rxt(:,k,tag_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_XMPAN_M_ndx) = 0._r8
          end if
       end if

!lke-TS1
!-----------------------------------------------------------------
!       ... pbznit + m --> acbzo2 + no2 + m
!-----------------------------------------------------------------
       if( usr_PBZNIT_M_ndx > 0 ) then
          if( tag_ACBZO2_NO2_ndx > 0 ) then
             rxt(:,k,usr_PBZNIT_M_ndx) = rxt(:,k,tag_ACBZO2_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_PBZNIT_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!       ... TERPAPAN + m --> TERPACO3 + no2 + m
!-----------------------------------------------------------------
       if( usr_TERPAPAN_M_ndx > 0 ) then
          if( tag_TERPACO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_TERPAPAN_M_ndx) = rxt(:,k,tag_TERPACO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_TERPAPAN_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!       ... TERPA2PAN + m --> TERPA2CO3 + no2 + m
!-----------------------------------------------------------------
       if( usr_TERPA2PAN_M_ndx > 0 ) then
          if( tag_TERPA2CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_TERPA2PAN_M_ndx) = rxt(:,k,tag_TERPA2CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_TERPA2PAN_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!       ... TERPA3PAN + m --> TERPA3CO3 + no2 + m
!-----------------------------------------------------------------
       if( usr_TERPA3PAN_M_ndx > 0 ) then
          if( tag_TERPA3CO3_NO2_ndx > 0 ) then
             rxt(:,k,usr_TERPA3PAN_M_ndx) = rxt(:,k,tag_TERPA3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
          else
             rxt(:,k,usr_TERPA3PAN_M_ndx) = 0._r8
          end if
       end if
!-----------------------------------------------------------------
!       ... xooh + oh -> h2o + oh
!-----------------------------------------------------------------
       if( usr_XOOH_OH_ndx > 0 ) then
          call comp_exp( exp_fac, 253._r8*tinv, ncol )
          rxt(:,k,usr_XOOH_OH_ndx) = temp(:ncol,k)**2._r8 * 7.69e-17_r8 * exp_fac(:)
       end if

!-----------------------------------------------------------------
!       ... ch3coch3 + oh -> ro2 + h2o
!-----------------------------------------------------------------
       if( usr_CH3COCH3_OH_ndx > 0 ) then
          call comp_exp( exp_fac, -2000._r8*tinv, ncol )
          rxt(:,k,usr_CH3COCH3_OH_ndx) = 3.82e-11_r8 * exp_fac(:) + 1.33e-13_r8
       end if

!-----------------------------------------------------------------
!       ... N + O2 -> NO + O
! Rate coefficients calculated by J. Orlando, D. Kinnison, and J. Zhang (2024)
! from data published in:
! "Kinetics of the Reactions of N(4S) Atoms with O2 and CO2 over Wide Temperatures Ranges"
! Abel Fernandez, A. Goumri, and Arthur Fontijn, 1998
! https://doi.org/10.1021/jp972365k
!-----------------------------------------------------------------
       if( usr_N_O2_ndx > 0 ) then
          call comp_exp( exp_fac, -2557._r8*tinv, ncol )
          rxt(:,k,usr_N_O2_ndx) = 2.0e-18_r8 * temp(:ncol,k)**2.15_r8 * exp_fac(:)
       end if

!-----------------------------------------------------------------
!       ... N2D + O2 -> NO + O
! "On the rate coefficient of the N(2D)+O2->NO+O reaction in the terrestrial thermosphere"
! Duff, J.W., H. Dothe, and R. D. Sharma, 2003
! https://doi.org/10.1029/2002GL016720
!-----------------------------------------------------------------
       if( usr_N2D_O2_ndx > 0 ) then
          rxt(:,k,usr_N2D_O2_ndx) = 6.2e-12_r8 * temp(:ncol,k)/300.0_r8
       end if

!-----------------------------------------------------------------
!       ... N2D + e -> N + e  Roble, 1995
! "Energetics of the Mesosphere and Thermosphere"
! R. Roble, 1995
! https://doi.org/10.1029/GM087p0001
!-----------------------------------------------------------------
       if( usr_N2D_e_ndx > 0 ) then
          rxt(:,k,usr_N2D_e_ndx) = 3.6e-10_r8 * sqrt(tempe(:ncol,k)/300.0_r8)
       end if

!-----------------------------------------------------------------
!       ... DMS + OH  --> .5 * SO2
!       JPL15-10 (use [O2] = 0.21*[M])
!       k = 8.2E-39 * exp(5376/T) * [O2] / (1 + 1.05E-5 *([O2]/[M]) * exp(3644/T))
!-----------------------------------------------------------------
       if( usr_DMS_OH_ndx > 0 ) then
!          call comp_exp( exp_fac, 7460._r8*tinv, ncol )
!          ko(:) = 1._r8 + 5.5e-31_r8 * exp_fac * m(:,k) * 0.21_r8
!          call comp_exp( exp_fac, 7810._r8*tinv, ncol )
!          rxt(:,k,usr_DMS_OH_ndx) = 1.7e-42_r8 * exp_fac * m(:,k) * 0.21_r8 / ko(:)
          call comp_exp( exp_fac, 3644._r8*tinv, ncol )
          ko(:) = 1._r8 + 1.05e-5_r8 * exp_fac * 0.21_r8
          call comp_exp( exp_fac, 5376._r8*tinv, ncol )
          rxt(:,k,usr_DMS_OH_ndx) = 8.2e-39_r8 * exp_fac * m(:,k) * 0.21_r8 / ko(:)
       end if

!-----------------------------------------------------------------
!       ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
!-----------------------------------------------------------------
       if( usr_SO2_OH_ndx > 0 ) then
          fc(:) = 3.0e-31_r8 *(300._r8*tinv(:))**3.3_r8
          ko(:) = fc(:)*m(:,k)/(1._r8 + fc(:)*m(:,k)/1.5e-12_r8)
          rxt(:,k,usr_SO2_OH_ndx) = ko(:)*.6_r8**(1._r8 + (log10(fc(:)*m(:,k)/1.5e-12_r8))**2._r8)**(-1._r8)
       end if
!RHS TS2
!-----------------------------------------------------------------
!       ... ISOPZD1O2 --> HPALD etc. Wennberg 2018 for rate
!-----------------------------------------------------------------
       if( usr_ISOPZD1O2_ndx > 0 ) then
          call comp_exp( exp_fac, -12200._r8*tinv, ncol )
          ko(:) = 5.05e15_r8 * exp_fac(:)
          call comp_exp( exp_fac, 1.e8_r8*tinv**3._r8, ncol )
          rxt(:,k,usr_ISOPZD1O2_ndx) = ko(:)*exp_fac(:)
       end if
!-----------------------------------------------------------------
!       ... ISOPZD4O2 --> HPALD etc. Wennberg 2018 for rate
!-----------------------------------------------------------------
       if( usr_ISOPZD4O2_ndx > 0 ) then
          call comp_exp( exp_fac, -7160._r8*tinv, ncol )
          ko(:) = 2.22e9_r8 * exp_fac(:)
          call comp_exp( exp_fac, 1.e8_r8*tinv**3._r8, ncol )
          rxt(:,k,usr_ISOPZD4O2_ndx) = ko(:)*exp_fac(:)
       end if
!-----------------------------------------------------------------
!       ... ISOPB1O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPB1O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.14_r8)/0.14_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPB1O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPB1O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPB4O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPB4O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.13_r8)/0.13_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPB4O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPB4O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPED1O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPED1O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.12_r8)/0.12_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
        rxt(:,k,usr_ISOPED1O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
        rxt(:,k,usr_ISOPED1O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPED4O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPED4O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.12_r8)/0.12_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPED4O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPED4O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPZD1O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPZD1O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.12_r8)/0.12_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPZD1O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPZD1O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPZD4O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPZD4O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.12_r8)/0.12_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPZD4O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPZD4O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPNO3_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPNO3_NOn_ndx > 0 ) then
          nyield = (1._r8-0.135_r8)/0.135_r8
          natom = 9.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPNO3_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPNO3_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... MVKO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_MVKO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.04_r8)/0.04_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_MVKO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_MVKO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... MACRO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_MACRO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.06_r8)/0.06_r8
          natom = 6.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
        rxt(:,k,usr_MACRO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
        rxt(:,k,usr_MACRO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... IEPOXOO_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_IEPOXOO_NOn_ndx > 0 ) then
          nyield = (1._r8-0.025_r8)/0.025_r8
          natom = 8.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_IEPOXOO_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_IEPOXOO_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPN1DO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPN1DO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.084_r8)/0.084_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPN1DO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPN1DO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPN2BO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPN2BO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.065_r8)/0.065_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
        rxt(:,k,usr_ISOPN2BO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
        rxt(:,k,usr_ISOPN2BO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPN3BO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPN3BO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.053_r8)/0.053_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPN3BO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPN3BO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPN4DO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPN4DO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.165_r8)/0.165_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPN4DO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPN4DO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPNBNO3O2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPNBNO3O2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.203_r8)/0.203_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPNBNO3O2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPNBNO3O2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPNOOHBO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPNOOHBO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.141_r8)/0.141_r8
          natom = 12.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPNOOHBO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPNOOHBO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... ISOPNOOHDO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_ISOPNOOHDO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.045_r8)/0.045_r8
          natom = 12.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_ISOPNOOHDO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_ISOPNOOHDO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!-----------------------------------------------------------------
!       ... NC4CHOO2_NOn Temp/Pressure Dependent Nitrate Yield
!-----------------------------------------------------------------
       if( usr_NC4CHOO2_NOn_ndx > 0 ) then
          nyield = (1._r8-0.021_r8)/0.021_r8
          natom = 11.0_r8
          exp_natom = exp( natom )
          acorr = (2.0e-22_r8*exp_natom*2.45e19_r8)/(1._r8+((2.0e-22_r8* &
                      exp_natom*2.45e19_r8)/(0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*2.45e19_r8)/ &
                      (0.43_r8*(298._r8*(1._r8/293._r8))**8._r8)))**2._r8))
          aterm(:) = (2.0e-22_r8*exp_natom*m(:,k))/(1._r8+((2.0e-22_r8* &
                      exp_natom*m(:,k))/(0.43_r8*(298._r8*tinv(:))**8._r8)))* &
                      0.41_r8**(1._r8/(1._r8+(log10((2.0e-22_r8*exp_natom*m(:,k))/ &
                      (0.43_r8*(298._r8*tinv(:))**8._r8)))**2._r8))
          call comp_exp( exp_fac, 360._r8*tinv, ncol )
          rxt(:,k,usr_NC4CHOO2_NOn_ndx) = 2.7e-12_r8 * exp_fac(:)*aterm(:)/(aterm(:)+acorr*nyield)
          rxt(:,k,usr_NC4CHOO2_NOa_ndx) = 2.7e-12_r8 * exp_fac(:)*acorr*nyield/(aterm(:)+acorr*nyield)
       end if
!
! reduced hydrocarbon scheme
!
       if ( usr_C2O3_NO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_NO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2O3_XNO2_ndx > 0 ) then
          ko(:)   = 2.6e-28_r8 * m(:,k)
          kinf(:) = 1.2e-11_r8
          rxt(:,k,usr_C2O3_XNO2_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log10(ko/kinf))**2))
       end if
       if ( usr_C2H4_OH_ndx > 0 ) then
          ko(:)   = 1.0e-28_r8 * m(:,k)
          kinf(:) = 8.8e-12_r8
          rxt(:,k,usr_C2H4_OH_ndx) = (ko/(1._r8+ko/kinf)) * 0.6_r8**(1._r8/(1._r8+(log(ko/kinf))**2))
       end if
       if ( usr_XO2N_HO2_ndx > 0 ) then
          rxt(:,k,usr_XO2N_HO2_ndx) = rxt(:,k,tag_XO2N_NO_ndx)*rxt(:,k,tag_XO2_HO2_ndx)/(rxt(:,k,tag_XO2_NO_ndx)+1.e-36_r8)
       end if

!------------------------------------------------------------------------------------
! WSY: iodine chemistry incorporated. originally from ordc (May 09, 2012)
!      based on experimental and theoretical data from Juan Carlos Gomez Martin)
!------------------------------------------------------------------------------------
       if ( usr_IO_IO_a_ndx > 0 ) then
          rxt(:,k,usr_IO_IO_a_ndx) = 2.13e-11_r8 * exp (180.0_r8 * tinv(:)) * (1._r8 + exp (- press_1d(:) / 191.42_r8))
       endif

       if ( usr_IO_IO_b_ndx > 0 ) then
          rxt(:,k,usr_IO_IO_b_ndx) = 3.27e-11_r8 * exp (180.0_r8 * tinv(:)) * (1._r8 - 0.65_r8 * exp (-press_1d(:)/191.42_r8))
       endif

       if ( usr_IO_OIO_ndx > 0 ) then
          work1 (:) = 4.687e-10_r8   - 1.3855e-5_r8 * exp (- press_1d(:) * 0.75_r8 / 1.62265_r8)       &
                    + 5.51868e-10_r8 * exp (- press_1d(:) * 0.75_r8 / 199.328_r8)
          work2 (:) = 0.00331_r8 + 0.00514_r8   * exp(- press_1d(:) * 0.75_r8 / 325.68711_r8)          &
                    + 0.00444_r8 * exp (- press_1d(:) * 0.75_r8 / 40.81609_r8)
          rxt (:,k, usr_IO_OIO_ndx) = work1(:) * exp (- work2(:) * temp_1d(:))

          !There might be some negative values for low pressures. We set them
          !to a low value (done also for other iodine reactions below)
          where (rxt (:,k, usr_IO_OIO_ndx) < 0._r8) rxt (:,k, usr_IO_OIO_ndx) = 3.0e-11_r8
       endif

       if ( usr_OIO_OIO_ndx > 0 ) then
          work1 (:) = 1.1659e-9_r8  - 7.79644e-10_r8 * exp(-press_1d(:)*0.75_r8/22.09281_r8)           &
                    + 1.03779e-9_r8 * exp (- press_1d(:) * 0.75_r8 / 568.15381_r8)
          work2 (:) = 0.00813_r8    + 0.00382_r8 * exp(- press_1d(:) * 0.75_r8 / 45.57591_r8)          &
                    + 0.00643_r8    * exp (-press_1d(:) * 0.75_r8 / 417.95061_r8)
          rxt (:,k, usr_OIO_OIO_ndx) = work1(:) * exp (- work2(:) * temp_1d(:))
          where (rxt (:,k, usr_OIO_OIO_ndx) < 0._r8) rxt (:,k, usr_OIO_OIO_ndx) = 2.0e-11_r8
       endif

       ! 20150325-CAC update with HOI+NO3-->IO+HNO3
       if ( usr_HOI_NO3_ndx > 0 ) then
          rxt(:,k,usr_HOI_NO3_ndx) = 2.7e-12_r8 * (300.0_r8 * tinv(:))**2.66_r8
       endif

       if ( usr_I2O2_a_ndx > 0 ) then
          work1 (:) =  3.54288e10_r8 + 1.8523e11_r8 * press_1d(:) * 0.75_r8 -1.45435e8_r8 * (press_1d(:) * 0.75_r8)**2._r8  &
                    +  60799.4344_r8 * (press_1d(:)*0.75_r8)**3._r8
          work2 (:) = -9681.65989_r8 + 346.95538_r8 * exp (-press_1d(:) * 0.75_r8 / 343.25322_r8)                           &
                    +   251.78032_r8 * exp (-press_1d(:) * 0.75_r8 / 44.1466_r8)
          rxt (:,k, usr_I2O2_a_ndx) =  work1(:) * exp (work2(:) * tinv (:))
          where (rxt (:,k, usr_I2O2_a_ndx) < 0._r8) rxt (:,k, usr_I2O2_a_ndx) = 4.0e-6_r8
       endif

       if ( usr_I2O2_b_ndx > 0 ) then
          work1 (:) = 255335000000._r8 - 4418880000._r8 * press_1d(:) * 0.75_r8                                             &
                    + 85618600._r8 * (press_1d(:) * 0.75_r8)**2._r8 + 14218.81_r8 * (press_1d(:) * 0.75_r8)**3._r8
          work2 (:) = -11466.82304_r8 + 597.01334_r8  * exp (-press_1d(:) * 0.75_r8 / 1382.62325_r8)                        &
                    - 167.3391_r8 * exp (- press_1d(:) * 0.75_r8 / 43.75089_r8)
          rxt (:,k, usr_I2O2_b_ndx) =  work1(:) * exp (work2(:) * tinv (:))
          where (rxt (:,k, usr_I2O2_b_ndx) < 0._r8) rxt (:,k, usr_I2O2_b_ndx) = 4.5e-10_r8
       endif

      if ( usr_I2O4_ndx > 0 ) then
          work1 (:) = -1.92626e14_r8 + 4.67414e13_r8 * press_1d(:) * 0.75_r8                                                &
                    -  3.68651e8_r8  * (press_1d(:) * 0.75_r8)**2._r8 - 3.09109e6_r8 * (press_1d(:) * 0.75_r8)**3._r8
          work2 (:) = -12302.15294_r8 + 252.78367_r8 * exp (- press_1d(:) * 0.75_r8 / 46.12733_r8)                          &
                    +  437.62868_r8   * exp (-press_1d(:) * 0.75_r8 / 428.4413_r8)
          rxt (:,k, usr_I2O4_ndx) =  work1(:) * exp (work2(:) * tinv (:))
          where (rxt (:,k, usr_I2O4_ndx) < 0._r8) rxt (:,k, usr_I2O4_ndx) = 3.0e-9_r8
       endif

       if ( usr_IONO2_ndx > 0 ) then
          work1 (:) = -2.63544e13_r8  + 4.32845e12_r8 * (press_1d(:) * 0.75_r8)                                             &
                    +  3.73758e8_r8   * (press_1d(:) * 0.75_r8)**2._r8 - 628468.76313_r8 * (press_1d(:) * 0.75_r8)**3._r8
          work2 (:) = -13847.85015_r8 + 240.34465_r8 * exp (- press_1d(:) * 0.75_r8 / 49.27141_r8)                          &
                    + 451.35864_r8    * exp (- press_1d(:) * 0.75_r8 / 436.87605_r8)
          rxt (:,k, usr_IONO2_ndx) =  work1(:) * exp (work2(:) * tinv (:))
          where (rxt (:,k, usr_IONO2_ndx) < 0._r8) rxt (:,k, usr_IONO2_ndx) = 6.0e-13_r8
       endif

!
! hydrolysis reactions on wetted aerosols
!
       if( usr_NO2_aer_ndx > 0 .or. usr_N2O5_HCL_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 .or. usr_HO2_aer_ndx > 0  &
         .or. usr_GLYOXAL_aer_ndx > 0 &
         .or. has_het_ss_rxts .or. has_ice_trp_rxts .or. has_ss_ixoy_rxts ) then

          long_loop : do i = 1,ncol

             sfc    => sfc_array(i,k,:)
             dm_aer => dm_array(i,k,:)

             press_lev   = pmid(i,k) / 100._r8   !hPa

             if ( id_hocl>0 ) then
                hoclvmr   = vmr(i,k,id_hocl)
                hocldeni  = 1._r8/(hoclvmr*m(i,k))
             endif
             if ( id_hcl>0 ) then
                hclvmr    = vmr(i,k,id_hcl)
                hcldeni  = 1._r8/(hclvmr*m(i,k))
             endif
             if ( id_hbr>0 ) then
                hbrvmr    = vmr(i,k,id_hbr)
                hbrdeni  = 1._r8/(hbrvmr*m(i,k))
             endif
             if ( id_hi>0 ) then
                hivmr     = vmr(i,k,id_hi)
                hideni  = 1._r8/(hivmr*m(i,k))
             endif
             if ( id_hobr>0 ) then
                hobrvmr   = vmr(i,k,id_hobr)
                hobrdeni  = 1._r8/(hobrvmr*m(i,k))
             endif
             if ( id_hoi>0 ) then
                hoivmr  = vmr(i,k,id_hoi)
                hoideni = 1._r8/(hoivmr*m(i,k))
             endif
             if ( id_n2o5>0 ) then
                n2o5vmr    = vmr(i,k,id_n2o5)
                n2o5deni  = 1._r8/(n2o5vmr*m(i,k))
             endif
             if ( id_clono2>0 ) clono2vmr   = vmr(i,k,id_clono2)
             if ( id_brono2>0 ) brono2vmr   = vmr(i,k,id_brono2)
             if ( id_iono2>0 )  iono2vmr    = vmr(i,k,id_iono2)

             c_n2o5 = 1.40e3_r8 * sqrt_t(i)         ! mean molecular speed of n2o5
             c_no3  = 1.85e3_r8 * sqrt_t(i)         ! mean molecular speed of no3
             c_no2  = 2.15e3_r8 * sqrt_t(i)         ! mean molecular speed of no2
             c_ho2  = 2.53e3_r8 * sqrt_t(i)         ! mean molecular speed of ho2
             c_glyoxal = 1.455e4_r8 * sqrt_t_58(i)  ! mean molecular speed of ho2
             c_isopnita = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopnita
             c_isopnitb = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopnitb
             c_onitr    = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of onitr
             c_honitr   = 1.26e3_r8 * sqrt_t(i)         ! mean molecular speed of honitr
             c_terpnit  = 0.992e3_r8 * sqrt_t(i)        ! mean molecular speed of terpnit
             c_nterpooh = 0.957e3_r8 * sqrt_t(i)        ! mean molecular speed of nterpooh
             c_nc4cho   = 1.21e3_r8 * sqrt_t(i)         ! mean molecular speed of nc4cho
             c_nc4ch2oh = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of nc4ch2oh
             c_isopfdn  = 9.68e2_r8 * sqrt_t(i)         ! mean molecular speed of isopfdn
             c_isopfnp  = 1.04e3_r8 * sqrt_t(i)         ! mean molecular speed of isopfnp
             c_isopn2b  = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopn2
             c_isopn1d  = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopn1d
             c_isopn4d  = 1.20e3_r8 * sqrt_t(i)         ! mean molecular speed of isopn4d
             c_inoohd   = 1.14e3_r8 * sqrt_t(i)         ! mean molecular speed of inoohd
             c_inheb    = 1.14e3_r8 * sqrt_t(i)         ! mean molecular speed of inheb
             c_inhed    = 1.14e3_r8 * sqrt_t(i)         ! mean molecular speed of inhed
             c_macrn    = 1.19e3_r8 * sqrt_t(i)         ! mean molecular speed of macrn
             c_isophfp  = 1.19e3_r8 * sqrt_t(i)         ! mean molecular speed of isophfp
             c_iepox    = 1.34e3_r8 * sqrt_t(i)         ! mean molecular speed of iepox
             c_dhpmpal  = 1.24e3_r8 * sqrt_t(i)         ! mean molecular speed of dhpmpal
             c_iche     = 1.35e3_r8 * sqrt_t(i)         ! mean molecular speed of iche
             c_isopfnc  = 1.04e3_r8 * sqrt_t(i)         ! mean molecular speed of isopfnc
             c_isopfdnc = 9.72e2_r8 * sqrt_t(i)         ! mean molecular speed of isopfdnc
             c_terpnt   = 9.92e2_r8 * sqrt_t(i)         ! mean molecular speed of terpnt
             c_terpnt1  = 9.92e2_r8 * sqrt_t(i)         ! mean molecular speed of terpnt1
             c_terpnpt  = 9.57e2_r8 * sqrt_t(i)         ! mean molecular speed of terpnpt
             c_terpnpt1 = 9.57e2_r8 * sqrt_t(i)         ! mean molecular speed of terpnpt1
             c_terpfdn  = 8.48e2_r8 * sqrt_t(i)         ! mean molecular speed of terpfdn
             c_sqtn     = 8.64e2_r8 * sqrt_t(i)         ! mean molecular speed of sqtn
             c_terphfn  = 8.93e2_r8 * sqrt_t(i)         ! mean molecular speed of terphfn
             c_terpdhdp = 9.47e2_r8 * sqrt_t(i)         ! mean molecular speed of terpdhdp
             c_terpacid = 1.07e3_r8 * sqrt_t(i)         ! mean molecular speed of terpacid
             !-------------------------------------------------------------------------
             !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
             !    rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))
             !-------------------------------------------------------------------------
             !-------------------------------------------------------------------------
             ! 	... n2o5 -> 2 hno3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_N2O5_aer_ndx > 0 ) then
                if( .NOT. usr_N2O5_aer1_ndx > 0 .AND. (.NOT. usr_N2O5_aer2_ndx > 0) ) then
                   rxt(i,k,usr_N2O5_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
                else
       ! --------------------------------
       ! VSLS Halogen chemistry incorporated
       ! --------------------------------
!                  total_sslt_mass = mmr(i,k,sslt1_ndx) + mmr(i,k,sslt2_ndx) + mmr(i,k,sslt3_ndx) + mmr(i,k,sslt4_ndx)
                   total_sslt_mass = 0._r8
                   if (sslt1_ndx>0) total_sslt_mass = total_sslt_mass + mmr(i,k,sslt1_ndx)
                   if (sslt2_ndx>0) total_sslt_mass = total_sslt_mass + mmr(i,k,sslt2_ndx)
                   if (sslt3_ndx>0) total_sslt_mass = total_sslt_mass + mmr(i,k,sslt3_ndx)
                   if (sslt4_ndx>0) total_sslt_mass = total_sslt_mass + mmr(i,k,sslt4_ndx)

                   if ( total_sslt_mass <= xmin ) then
                      clno2_yield(i,k) = ymin
                   elseif ( total_sslt_mass >= xmax ) then
                      clno2_yield(i,k) = ymax
                   else
                      clno2_yield(i,k) = (ymax-ymin) * ((total_sslt_mass-xmin)/(xmax-xmin))**3 + ymin
                   end if

                   dummy_n2o5_clno2_rate      =                              hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 ) ! original
                   rxt(i,k,usr_N2O5_aer_ndx ) =                              dummy_n2o5_clno2_rate
                   rxt(i,k,usr_N2O5_aer1_ndx) = (1._r8 - clno2_yield_Erin) * dummy_n2o5_clno2_rate
                   if( n2o5vmr > small ) then
                      if( hclvmr > small ) then
                         if ( n2o5vmr > hclvmr ) then
                            rxt(i,k,usr_N2O5_aer2_ndx) =          clno2_yield_Erin  * dummy_n2o5_clno2_rate * n2o5deni
                         else
                            rxt(i,k,usr_N2O5_aer2_ndx) =          clno2_yield_Erin  * dummy_n2o5_clno2_rate * hcldeni
                         end if
                      else
                           rxt(i,k,usr_N2O5_aer1_ndx) = 1._r8 * dummy_n2o5_clno2_rate
                           rxt(i,k,usr_N2O5_aer2_ndx) = 0._r8
                      endif
                   else
                           rxt(i,k,usr_N2O5_aer1_ndx) = 1._r8 * dummy_n2o5_clno2_rate
                           rxt(i,k,usr_N2O5_aer2_ndx) = 0._r8
                   endif
                end if
             end if

             if( usr_XNO2NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             if( usr_NO2XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_n2o5, gamma_n2o5 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no3 -> hno3  (on sulfate, nh4no3, oc, soa)
             !-------------------------------------------------------------------------
             if( usr_NO3_aer_ndx > 0 ) then
                rxt(i,k,usr_NO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             if( usr_XNO3_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO3_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no3, gamma_no3 )
             end if
             !-------------------------------------------------------------------------
             ! 	... no2 -> 0.5 * (ho+no+hno3)  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NO2_aer_ndx > 0 ) then
                rxt(i,k,usr_NO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
             end if
             if( usr_XNO2_aer_ndx > 0 ) then
                rxt(i,k,usr_XNO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_no2, gamma_no2 )
             end if

             !-------------------------------------------------------------------------
             ! 	... ho2 -> 0.5 * h2o2  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_HO2_aer_ndx > 0 ) then
                rxt(i,k,usr_HO2_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_ho2, gamma_ho2 )
             end if
             !-------------------------------------------------------------------------
             !  ... glyoxal ->  soag1  (on sulfate, nh4no3, oc2, soa)
             ! first order uptake, Fuchs and Sutugin, 1971,  dCg = 1/4 * gamma * ! A * |v_mol| * Cg * dt
             !-------------------------------------------------------------------------
             if( usr_GLYOXAL_aer_ndx > 0 ) then
                rxt(i,k,usr_GLYOXAL_aer_ndx) = hetrxtrate_gly( sfc, c_glyoxal, gamma_glyoxal )
             end if
             !-------------------------------------------------------------------------
             ! 	... ISOPNITA -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPNITA_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPNITA_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopnita, gamma_isopnita )
             end if
             !-------------------------------------------------------------------------
             ! 	... ISOPNITB -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPNITB_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPNITB_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopnitb, gamma_isopnitb )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ONITR -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ONITR_aer_ndx > 0 ) then
                rxt(i,k,usr_ONITR_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_onitr, gamma_onitr )
             end if
             !-------------------------------------------------------------------------
             ! 	... HONITR -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_HONITR_aer_ndx > 0 ) then
                rxt(i,k,usr_HONITR_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_honitr, gamma_honitr )
             end if
             !-------------------------------------------------------------------------
             ! 	... TERPNIT -> HNO3  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNIT_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNIT_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnit, gamma_terpnit )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NTERPOOH -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NTERPOOH_aer_ndx > 0 ) then
                rxt(i,k,usr_NTERPOOH_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nterpooh, gamma_nterpooh )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NC4CHO -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NC4CHO_aer_ndx > 0 ) then
                rxt(i,k,usr_NC4CHO_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nc4cho, gamma_nc4cho )
             end if
             !-------------------------------------------------------------------------
             ! 	...  NC4CH2OH -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_NC4CH2OH_aer_ndx > 0 ) then
                rxt(i,k,usr_NC4CH2OH_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_nc4ch2oh, gamma_nc4ch2oh )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPFDN -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPFDN_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPFDN_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopfdn, gamma_isopfdn )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPFNP ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPFNP_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPFNP_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopfnp, gamma_isopfnp )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPN2B -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPN2B_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPN2B_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopn2b, gamma_isopn2b )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPN1D -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPN1D_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPN1D_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopn1d, gamma_isopn1d )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPN4D -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPN4D_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPN4D_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopn4d, gamma_isopn4d )
             end if
             !-------------------------------------------------------------------------
             ! 	...  INOOHD -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_INOOHD_aer_ndx > 0 ) then
                rxt(i,k,usr_INOOHD_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_inoohd, gamma_inoohd )
             end if
             !-------------------------------------------------------------------------
             ! 	...  INHEB -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_INHEB_aer_ndx > 0 ) then
                rxt(i,k,usr_INHEB_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_inheb, gamma_inheb )
             end if
             !-------------------------------------------------------------------------
             ! 	...  INHED -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_INHED_aer_ndx > 0 ) then
                rxt(i,k,usr_INHED_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_inhed, gamma_inhed )
             end if
             !-------------------------------------------------------------------------
             ! 	...  MACRN -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_MACRN_aer_ndx > 0 ) then
                rxt(i,k,usr_MACRN_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_macrn, gamma_macrn )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPHFP ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPHFP_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPHFP_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isophfp, gamma_isophfp )
             end if
             !-------------------------------------------------------------------------
             ! 	...  IEPOX ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_IEPOX_aer_ndx > 0 ) then
                rxt(i,k,usr_IEPOX_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_iepox, gamma_iepox )
             end if
             !-------------------------------------------------------------------------
             ! 	...  DHPMPAL ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_DHPMPAL_aer_ndx > 0 ) then
                rxt(i,k,usr_DHPMPAL_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_dhpmpal, gamma_dhpmpal )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ICHE ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ICHE_aer_ndx > 0 ) then
                rxt(i,k,usr_ICHE_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_iche, gamma_iche )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPFNC ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPFNC_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPFNC_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopfnc, gamma_isopfnc )
             end if
             !-------------------------------------------------------------------------
             ! 	...  ISOPFDNC ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_ISOPFDNC_aer_ndx > 0 ) then
                rxt(i,k,usr_ISOPFDNC_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_isopfdnc, gamma_isopfdnc )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPNT -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNT_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNT_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnt, gamma_terpnt )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPNT1 -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNT1_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNT1_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnt1, gamma_terpnt1 )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPNPT -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNPT_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNPT_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnpt, gamma_terpnpt )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPNPT1 -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPNPT1_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPNPT1_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpnpt1, gamma_terpnpt1 )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPFDN -> HNO3 (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPFDN_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPFDN_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpfdn, gamma_terpfdn )
             end if
             !-------------------------------------------------------------------------
             ! 	...  SQTN ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_SQTN_aer_ndx > 0 ) then
                rxt(i,k,usr_SQTN_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_sqtn, gamma_sqtn )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPHFN ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPHFN_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPHFN_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terphfn, gamma_terphfn )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPDHDP ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPDHDP_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPDHDP_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpdhdp, gamma_terpdhdp )
             end if
             !-------------------------------------------------------------------------
             ! 	...  TERPACID ->  (on sulfate, nh4no3, oc2, soa)
             !-------------------------------------------------------------------------
             if( usr_TERPACID_aer_ndx > 0 ) then
                rxt(i,k,usr_TERPACID_aer_ndx) = hetrxtrate( sfc, dm_aer, dg, c_terpacid, gamma_terpacid )
             end if

!-------------------------------------------------------------------------
!  ... WSY: VSLS chemistry incorporated
!-------------------------------------------------------------------------
             ssur_tot  = sad_sslt(i,k)
             ssur_tot  = min(3.0e-7_r8, ssur_tot)   ! ssur_tot is not longer used ... using sad_sslt_mask
             sad_sslt_mask  = sad_sslt(i,k)

             latitude = dlat(i)
             if ( ocnfrac(i) + icefrac(i) <= 0._r8 ) then
               if (latitude > -60.0_r8) then
                 ssur_tot      = 0._r8
                 sad_sslt_mask = 0._r8
               endif
             end if

             !rpf (Sept 18, 2012)_CAM-Chem_3.6.x_polar update: Apply Yang DF to cuped SSLT (ssur_tot) and not sad_sslt:
             ! The cup physical meaning is to avoid overestimation of SSLT production over the Southern Ocean
             ! The DF physical meaning is that the seasonal change in acidity modifies efficiency of recycling
             if (latitude > -30.0_r8) then
                DF = 0.5_r8
             else
                DF = dfmax + 0.5_r8  * (dfmin - dfmax) * (sin ((calday/(365.0_r8/2.0_r8) - 0.5_r8 ) * pi) + 1.0_r8 )
             endif

             !rpf (Feb 17, 2013): introduce logical condition to avoid SSLT_recycling in the UTLS and lower stratosphere
             if (press_lev < 300._r8) then
                logical_sslt = 0.0_r8
             else
                logical_sslt = 1.0_r8
             endif

             sad_sslt_eff(i,k) = sad_sslt_mask * DF * logical_sslt   ! for SSLT SAD Output

             !**********************************************************************************
             ! Here the uptake of halogen species is considered to be the rate limiting step
             ! of the recycling on seasalt aerosols (see Supplement of Ordonez et al., ACP, 2012).
             ! Therefore the reaction rates are given by the transfer coefs. k:
             !    k = (gamma / 4) * c * A
             !  where
             !    gamma = uptake coefficient
             !    A     = aerosol surface area
             !    c     = 100 * sqrt {(8*R*T)/pi*M} = root-mean square of molecular speed
             !            (the factor 100 converts to cm s-1)
             !
             ! Note that in the formulation below 100 * sqrt {(8*R)/pi*M} has already been calculated
             ! for the molecular masses M of the species uptaken
             !
             ! Unlike in CAMChem 3.6.x, the recycling of halogens on seasalt and any other tropospheric
             ! aerosols is now calculated even if the available aerosol surface area density is null.
             ! This is needed because the reaction rates are now initilized with NANs in gas_phase_chemdr.
             !
             ! In addition, for the heterog release of BR2 and BRCL we follow a similar approach
             ! to that of Yang et al., JGR, 2005:
             ! * Try to keep the total SAD of sea-salt (not the effective one)
             ! * Use bromine depletion factors (DF) of 0.5 North of 30 S and a sinusoidal formula
             !   to the south (average DF of 0.4 but with strong seasonality there)
             !**********************************************************************************

             if ( het_ss_0_ndx > 0 ) then
!              rxt(i,k,het_ss_0_ndx) = 0.25_r8 * gamma_brono2_ss * sad_sslt_mask * 1.22e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_0_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_brono2_ss * sad_sslt_mask * 1.22e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_1_ndx > 0 ) then
!              rxt(i,k,het_ss_1_ndx) = 0.25_r8 * gamma_brno2_ss  * sad_sslt_mask * 1.29e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_1_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_brno2_ss  * sad_sslt_mask * 1.29e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_2_ndx > 0 ) then
!              rxt(i,k,het_ss_2_ndx) = 0.25_r8 * gamma_hobr_ss   * sad_sslt_mask * 1.47e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_2_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_hobr_ss   * sad_sslt_mask * 1.47e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_3_ndx > 0 ) then
!              rxt(i,k,het_ss_3_ndx) = 0.25_r8 * gamma_clono2_ss * sad_sslt_mask * 1.47e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_3_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_clono2_ss * sad_sslt_mask * 1.47e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_4_ndx > 0 ) then
!              rxt(i,k,het_ss_4_ndx) = 0.25_r8 * gamma_clno2_ss  * sad_sslt_mask * 1.61e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_4_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_clno2_ss  * sad_sslt_mask * 1.61e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_5_ndx > 0 ) then
!              rxt(i,k,het_ss_5_ndx) = 0.25_r8 * gamma_hocl_ss   * sad_sslt_mask * 2.01e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_5_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_hocl_ss   * sad_sslt_mask * 2.01e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_6_ndx > 0 ) then
!              rxt(i,k,het_ss_6_ndx) = 0.25_r8 * gamma_iono2_ss  * sad_sslt_mask * 1.05e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_6_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_iono2_ss  * sad_sslt_mask * 1.05e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_7_ndx > 0 ) then
!              rxt(i,k,het_ss_7_ndx) = 0.25_r8 * gamma_ino2_ss   * sad_sslt_mask * 1.10e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_7_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_ino2_ss   * sad_sslt_mask * 1.10e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_8_ndx > 0 ) then
!              rxt(i,k,het_ss_8_ndx) = 0.25_r8 * gamma_hoi_ss    * sad_sslt_mask * 1.21e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_8_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_hoi_ss    * sad_sslt_mask * 1.21e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( het_ss_9_ndx > 0 ) then
!              rxt(i,k,het_ss_9_ndx) = 0.25_r8 * gamma_hno3_ss   * sad_sslt_mask * 1.83e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,het_ss_9_ndx) = SSAhno3_ScalingFactor  * 0.25_r8 * gamma_hno3_ss   * sad_sslt_mask * 1.83e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif


             if ( het_ss_10_ndx > 0 ) then
!              rxt(i,k,het_ss_10_ndx)=                               0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
               rxt(i,k,het_ss_10_ndx)= SSAn2o5_ScalingFactor  *                              0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
             endif
             if ( het_ss_11_ndx > 0 ) then
!              rxt(i,k,het_ss_11_ndx) = (1._r8 - clno2_yield(i,k)) * 0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
               rxt(i,k,het_ss_11_ndx) = SSAn2o5_ScalingFactor * (1._r8 - clno2_yield(i,k)) * 0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
             endif
             if ( het_ss_12_ndx > 0 ) then
!              rxt(i,k,het_ss_12_ndx) =           clno2_yield(i,k) * 0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
               rxt(i,k,het_ss_12_ndx) = SSAn2o5_ScalingFactor *           clno2_yield(i,k) * 0.25_r8 * gamma_n2o5_ss   * sad_sslt_mask * 1.4e3_r8  * sqrt(temp(i,k)) * logical_sslt
             endif


             if ( ss_ixoy_2_ndx > 0 ) then
!              rxt(i,k,ss_ixoy_2_ndx) = 0.25 * gamma_ss_ixoy_2 * sad_sslt_mask * 0.8607e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,ss_ixoy_2_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_ss_ixoy_2 * sad_sslt_mask * 0.8607e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( ss_ixoy_3_ndx > 0 ) then
!              rxt(i,k,ss_ixoy_3_ndx) = 0.25 * gamma_ss_ixoy_3 * sad_sslt_mask * 0.8376e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,ss_ixoy_3_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_ss_ixoy_3 * sad_sslt_mask * 0.8376e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif
             if ( ss_ixoy_4_ndx > 0 ) then
!              rxt(i,k,ss_ixoy_4_ndx) = 0.25 * gamma_ss_ixoy_4 * sad_sslt_mask * 0.8162e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
               rxt(i,k,ss_ixoy_4_ndx) = SSAdehal_ScalingFactor * 0.25_r8 * gamma_ss_ixoy_4 * sad_sslt_mask * 0.8162e3_r8 * sqrt(temp(i,k)) * DF * logical_sslt
             endif


             !*********************************************************************************************
             !  ICE_SAD reactions
             !*********************************************************************************************
             if ( ice_trp_hcl_5_ndx > 0 ) then     ! [ice_trp_hcl_5]    HOCL + HCL -> CL2 + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hoclvmr > small ) then
                      if( hclvmr > small ) then
                         ! I will aways consider velocity (c = 2.01e3) of HOCL (HOBR) to keep the mechanism similar to
                         ! how it was implemented in mo_strato_rates.F90
                         if ( hoclvmr > hclvmr ) then
                            rxt(i,k,ice_trp_hcl_5_ndx) = 0.25_r8 * gamma_hclhocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hocldeni
                         else
                            rxt(i,k,ice_trp_hcl_5_ndx) = 0.25_r8 * gamma_hclhocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hcldeni
                         end if
! This is another way of implementing, considering the limiting rate instead of the limiting vmr when computing rxt(i,k)
!                    half_rate_1 = 0.25 * gamma_hbrhobr_ice * sad_ice_trop(i,k) * 1.6e3_r8  * sqrt(temp(i,k)) * hobrdeni  ! limiting reactant is HBR
!                    half_rate_2 = 0.25 * gamma_hbrhobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hbrdeni   ! limiting reactant is HOBR
!                    rxt(i,k,het_ice_6_ndx) = min (half_rate_1,half_rate_2)
                      endif
                   endif
                endif
             endif
             if ( ice_trp_hcl_6_ndx > 0 ) then     ! [ice_trp_hcl_6]    HOBR + HCL -> BRCL + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hobrvmr > small ) then
                      if( hclvmr > small ) then
                         if ( hobrvmr > hclvmr ) then
                            rxt(i,k,ice_trp_hcl_6_ndx) = 0.25_r8 * gamma_hclhobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hobrdeni
                         else
                            rxt(i,k,ice_trp_hcl_6_ndx) = 0.25_r8 * gamma_hclhobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hcldeni
                         end if
                      endif
                   endif
                endif
             endif

             if ( usr_N2O5_HCL_ndx > 0 ) then     ! [usr_N2O5_HCL]    N2O5 + HCL -> CLNO2+HNO3
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( n2o5vmr > small ) then
                      if( hclvmr > small ) then
                         if ( n2o5vmr > hclvmr ) then
                            rxt(i,k,usr_N2O5_HCL_ndx) = clno2_yield_Erin * 0.25_r8 * gamma_hcln2o5_aerosol * sad_ice_trop(i,k) * 1.40e3_r8 * sqrt(temp(i,k)) * n2o5deni
                         else
                            rxt(i,k,usr_N2O5_HCL_ndx) = clno2_yield_Erin * 0.25_r8 * gamma_hcln2o5_aerosol * sad_ice_trop(i,k) * 1.40e3_r8 * sqrt(temp(i,k)) * hcldeni
                         end if
                      endif
                   endif
                endif
             endif

             if ( ice_trp_hbr_5_ndx > 0 ) then     ! [ice_trp_hbr_5]    HOCL + HBR -> BRCL + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hoclvmr > small ) then
                      if( hbrvmr > small ) then
                         if ( hoclvmr > hbrvmr ) then
                            rxt(i,k,ice_trp_hbr_5_ndx) = 0.25_r8 * gamma_hbrhocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hocldeni
                         else
                            rxt(i,k,ice_trp_hbr_5_ndx) = 0.25_r8 * gamma_hbrhocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hbrdeni
                         end if
                      endif
                   endif
                endif
             endif
             if ( ice_trp_hbr_6_ndx > 0 ) then     ! [ice_trp_hbr_6]    HOBR + HBR -> BR2 + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hobrvmr > small ) then
                      if( hbrvmr > small ) then
                         if ( hobrvmr > hbrvmr ) then
                            rxt(i,k,ice_trp_hbr_6_ndx) = 0.25_r8 * gamma_hbrhobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hobrdeni
                         else
                            rxt(i,k,ice_trp_hbr_6_ndx) = 0.25_r8 * gamma_hbrhobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hbrdeni
                         end if
                      endif
                   endif
                endif
             endif

             if ( ice_trp_hi_5_ndx > 0 ) then     ! [ice_trp_hi_5]    HOCL + HI  -> ICL + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hoclvmr > small ) then
                      if( hivmr > small ) then
                         if ( hoclvmr > hivmr ) then
                            rxt(i,k,ice_trp_hi_5_ndx) = 0.25_r8 * gamma_hihocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hocldeni
                         else
                            rxt(i,k,ice_trp_hi_5_ndx) = 0.25_r8 * gamma_hihocl_ice * sad_ice_trop(i,k) * 2.01e3_r8 * sqrt(temp(i,k)) * hideni
                         end if
                      endif
                   endif
                endif
             endif
             if ( ice_trp_hi_6_ndx > 0 ) then     ! [ice_trp_hi_6]    HOBR + HI  -> IBR + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hobrvmr > small ) then
                      if( hivmr > small ) then
                         if ( hobrvmr > hivmr ) then
                            rxt(i,k,ice_trp_hi_6_ndx) = 0.25_r8 * gamma_hihobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hobrdeni
                         else
                            rxt(i,k,ice_trp_hi_6_ndx) = 0.25_r8 * gamma_hihobr_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k)) * hideni
                         end if
                      endif
                   endif
                endif
             endif

             if ( ice_trp_i_2_ndx > 0 ) then     ! [ice_trp_i_2]    HOI + HCL -> ICL + H2O
             if ( sad_ice_trop(i,k) > 0 ) then
               if( hoivmr > small ) then
                 if( hclvmr > small ) then
                   if ( hoivmr > hclvmr ) then
                     rxt(i,k,ice_trp_i_2_ndx) = 0.25_r8 * gamma_hclhoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hoideni
                   else
                     rxt(i,k,ice_trp_i_2_ndx) = 0.25_r8 * gamma_hclhoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hcldeni
                   end if
                 endif
               endif
             endif
             endif
             if ( ice_trp_i_3_ndx > 0 ) then     ! [ice_trp_i_3]    HOI + HBR -> IBR + H2O
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hoivmr > small ) then
                      if( hbrvmr > small ) then
                         if ( hoivmr > hbrvmr ) then
                            rxt(i,k,ice_trp_i_3_ndx) = 0.25_r8 * gamma_hbrhoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hoideni
                         else
                            rxt(i,k,ice_trp_i_3_ndx) = 0.25_r8 * gamma_hbrhoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hbrdeni
                         end if
                      endif
                   endif
                endif
             endif
             if ( ice_trp_i_4_ndx > 0 ) then     ! [ice_trp_i_4]    HOI + HI -> I2 + H2O
             if ( sad_ice_trop(i,k) > 0 ) then
               if( hoivmr > small ) then
                 if( hivmr > small ) then
                   if ( hoivmr > hivmr ) then
                     rxt(i,k,ice_trp_i_4_ndx) = 0.25_r8 * gamma_hihoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hoideni
                   else
                     rxt(i,k,ice_trp_i_4_ndx) = 0.25_r8 * gamma_hihoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k)) * hideni
                   end if
                 endif
               endif
             endif
             endif

             if ( ice_trp_cl_1_ndx > 0 ) then     ! [ice_trp_cl_1]    CLONO2 -> HOCL + HNO3
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( clono2vmr > small ) then
                      rxt(i,k,ice_trp_cl_1_ndx) = 0.25_r8 * gamma_clono2_ice * sad_ice_trop(i,k) * 1.47e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( ice_trp_br_1_ndx > 0 ) then     ! [ice_trp_br_1]    BRONO2 -> HOBR + HNO3
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( brono2vmr > small ) then
                      rxt(i,k,ice_trp_br_1_ndx) = 0.25_r8 * gamma_brono2_ice * sad_ice_trop(i,k) * 1.22e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( ice_trp_i_1_ndx > 0 ) then     ! [ice_trp_i_1]    IONO2 -> HOI + HNO3   .or.   IONO2 -> deposition
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( iono2vmr > small ) then
                      rxt(i,k,ice_trp_i_1_ndx) = 0.25_r8 * gamma_iono2_ice * sad_ice_trop(i,k) * 1.06e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif

             !*********************************************************************************************
             !  HOI removal on liquid and ice in the Troposphere. Also applied to HI and IONO2
             !*********************************************************************************************
             if ( ice_fr_hoi_ndx > 0 ) then     ! [ice_fr_hoi]      HOI ->
                if ( sad_ice_trop(i,k) > 0._r8 ) then
                   if( hoivmr > small ) then
                    ! rxt(i,k,ice_fr_hoi_ndx) = 0.25_r8 * gamma_fr_hoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,ice_fr_hoi_ndx) = ICEfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_hoi_ice * sad_ice_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( liq_fr_hoi_ndx > 0 ) then     ! [ice_fr_hoi]      HOI ->
                if ( sad_liq_trop(i,k) > 0._r8 ) then
                   if( hoivmr > small ) then
                    ! rxt(i,k,liq_fr_hoi_ndx) = 0.25_r8 * gamma_fr_hoi_liq * sad_liq_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,liq_fr_hoi_ndx) = LIQfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_hoi_liq * sad_liq_trop(i,k) * 1.21e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( ice_fr_hi_ndx > 0 ) then     ! [ice_fr_hi]      HI ->
                if ( sad_ice_trop(i,k) > 0 ) then
                   if( hivmr > small ) then
                    ! rxt(i,k,ice_fr_hi_ndx) = 0.25 * gamma_fr_hi_ice * sad_ice_trop(i,k) * 1.287e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,ice_fr_hi_ndx) = ICEfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_hi_ice * sad_ice_trop(i,k) * 1.287e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( liq_fr_hi_ndx > 0 ) then     ! [liq_fr_hi]      HI ->
                if ( sad_liq_trop(i,k) > 0 ) then
                   if( hivmr > small ) then
                    ! rxt(i,k,liq_fr_hi_ndx) = 0.25 * gamma_fr_hi_liq * sad_liq_trop(i,k) * 1.287e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,liq_fr_hi_ndx) = LIQfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_hi_liq * sad_liq_trop(i,k) * 1.287e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( ice_fr_iono2_ndx > 0 ) then     ! [ice_fr_iono2]      IONO2 ->
                if ( sad_ice_trop(i,k) > 0 ) then
                   if( iono2vmr > small ) then
                    ! rxt(i,k,ice_fr_iono2_ndx) = 0.25 * gamma_fr_iono2_ice * sad_ice_trop(i,k) * 1.059e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,ice_fr_iono2_ndx) = ICEfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_iono2_ice * sad_ice_trop(i,k) * 1.059e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif
             if ( liq_fr_iono2_ndx > 0 ) then     ! [liq_fr_iono2]      IONO2 ->
                if ( sad_liq_trop(i,k) > 0 ) then
                   if( iono2vmr > small ) then
                    ! rxt(i,k,liq_fr_iono2_ndx) = 0.25 * gamma_fr_iono2_liq * sad_liq_trop(i,k) * 1.059e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,liq_fr_iono2_ndx) = LIQfraprx_ScalingFactor_I * 0.25_r8 * gamma_fr_iono2_liq * sad_liq_trop(i,k) * 1.059e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif

             !*********************************************************************************************
             !  BRONO2 removal on ice in the Troposphere
             !*********************************************************************************************
             if ( ice_fr_brono2_ndx > 0 ) then     ! [ice_fr_brono2]      BRONO2 ->
                if ( sad_ice_trop(i,k) > 0 ) then
                   if( brono2vmr > small ) then
                    ! rxt(i,k,ice_fr_brono2_ndx) = 0.25 * gamma_fr_brono2_ice * sad_ice_trop(i,k) * 1.22e3_r8 * sqrt(temp(i,k))
                      rxt(i,k,ice_fr_brono2_ndx) = ICEfraprx_ScalingFactor_Br * 0.25_r8 * gamma_fr_brono2_ice * sad_ice_trop(i,k) * 1.22e3_r8 * sqrt(temp(i,k))
                   endif
                endif
             endif


          end do long_loop
       end if

       ! LLNL super fast chem reaction rates

       !-----------------------------------------------------------------------
       !     ... CO + OH --> CO2 + HO2
       !-----------------------------------------------------------------------
       if ( usr_oh_co_ndx > 0 ) then
          ko(:)     = 5.9e-33_r8 * tp(:)**1.4_r8
          kinf(:)   = 1.1e-12_r8 * (temp(:ncol,k) / 300._r8)**1.3_r8
          ko_m(:)   = ko(:) * m(:,k)
          k0(:)     = 1.5e-13_r8 * (temp(:ncol,k) / 300._r8)**0.6_r8
          kinf_m(:) = (2.1e+09_r8 * (temp(:ncol,k) / 300._r8)**6.1_r8) / m(:,k)
          rxt(:,k,usr_oh_co_ndx) = (ko_m(:)/(1._r8+(ko_m(:)/kinf(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(ko_m(:)/kinf(:)))**2._r8)) + &
               (k0(:)/(1._r8+(k0(:)/kinf_m(:)))) * &
               0.6_r8**(1._r8/(1._r8+(log10(k0(:)/kinf_m(:)))**2._r8))
       endif
       !-----------------------------------------------------------------------
       !     ... NO2 + H2O --> 0.5 HONO + 0.5 HNO3
       !-----------------------------------------------------------------------
       if ( het_no2_h2o_ndx > 0 ) then
          rxt(:,k,het_no2_h2o_ndx) = 4.0e-24_r8
       endif
       !-----------------------------------------------------------------------
       !     ... DMS + OH --> 0.75 SO2 + 0.25 MSA
       !-----------------------------------------------------------------------
       if ( usr_oh_dms_ndx > 0 ) then
          o2(:ncol) = invariants(:ncol,k,inv_o2_ndx)
          rxt(:,k,usr_oh_dms_ndx) = 2.000e-10_r8 * exp(5820.0_r8 * tinv(:)) / &
               ((2.000e29_r8 / o2(:)) + exp(6280.0_r8 * tinv(:)))
       endif
       if ( aq_so2_h2o2_ndx > 0 .or. aq_so2_o3_ndx > 0 ) then
          lwc(:) = cwat(:ncol,k) * invariants(:ncol,k,inv_m_ndx) * mbar(:ncol,k) /avo !PJC convert kg/kg to g/cm3
          !-----------------------------------------------------------------------
          !     ... SO2 + H2O2 --> S(VI)
          !-----------------------------------------------------------------------
          if ( aq_so2_h2o2_ndx > 0 ) then
             rxt(:,k,aq_so2_h2o2_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  K_AQ * &

                  exp(ER_AQ * ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  HENRY298_SO2 * &
                  K298_SO2_HSO3 * &
                  HENRY298_H2O2 * &
                  exp(((H298_SO2 + H298_SO2_HSO3 + H298_H2O2) / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8 / &

                  (1.0e+00_r8 + 13.0e+00_r8 * 10.0e+00_r8**(-pH))
          endif
          !-----------------------------------------------------------------------
          !     ... SO2 + O3 --> S(VI)
          !-----------------------------------------------------------------------
          if (aq_so2_o3_ndx >0) then
             rxt(:,k,aq_so2_o3_ndx) = lwc(:) * 1.0e-03_r8 * avo * &
                  HENRY298_SO2 * exp((H298_SO2 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K0_AQ * exp(ER0_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) + &
                  K298_SO2_HSO3 * exp((H298_SO2_HSO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (K1_AQ * exp(ER1_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  10.0e+00_r8**(-pH) + K2_AQ * exp(ER2_AQ * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  K298_HSO3_SO3 * exp((H298_HSO3_SO3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) / &
                  (10.0e+00_r8**(-pH))**2.0e+00_r8) ) * &
                  HENRY298_O3 * exp((H298_O3 / R_CAL) * &
                  ((1.0e+00_r8 / 298.0e+00_r8) - tinv(:))) * &
                  (R_CONC * temp(:ncol,k))**2.0e+00_r8
          endif
       endif

    if ( has_d_chem ) then

        call comp_exp( exp_fac, -600._r8 * tinv, ncol )
        rxt(:,k,ean_ndx(1))  = 1.e-31_r8 * tp(:) * exp_fac(:)
        rxt(:,k,ean_ndx(2))  = 9.1e-12_r8 * tp(:)**(-1.46_r8)
        call comp_exp( exp_fac, -193._r8 * tinv, ncol )
        rxt(:,k,ean_ndx(3))  = (4.e-30_r8 * exp_fac(:)) * 0.21_r8

        rxt(:,k,rpe_ndx(1))  = 4.2e-6_r8 * tp(:)**0.5_r8
        rxt(:,k,rpe_ndx(2))  = 6.3e-7_r8 * tp(:)**0.5_r8
        rxt(:,k,rpe_ndx(3))  = 2.5e-6_r8 * tp(:)**0.1_r8
        rxt(:,k,rpe_ndx(4))  = 2.48e-6_r8 * tp(:)**0.76_r8
        rxt(:,k,rpe_ndx(5))  = 1.4e-6_r8 * tp(:)**0.4_r8

        rxt(:,k,pir_ndx(1)) = 4.e-30_r8 * tp(:)**2.93_r8
        rxt(:,k,pir_ndx(2))  = 4.6e-27_r8 * tp(:)**4._r8

        call comp_exp( exp_fac, -15900._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(3))  = (2.5e-2_r8 * tp(:)**5._r8) * exp_fac(:)
        rxt(:,k,pir_ndx(4))  = 2.3e-27_r8 * tp(:)**7.5_r8

        call comp_exp( exp_fac, -10272._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(5))  = (2.6e-3_r8 * tp(:)**8.5_r8) * exp_fac(:)
        rxt(:,k,pir_ndx(6))  = 3.6e-27_r8 * tp(:)**8.1_r8

        call comp_exp( exp_fac, -9000._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(7))  = (1.5e-1_r8 * tp(:)**9.1_r8) * exp_fac(:)
        rxt(:,k,pir_ndx(8))  = 4.6e-28_r8 * tp(:)**14._r8

        call comp_exp( exp_fac, -6400._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(9))  = (1.7e-3_r8 * tp(:)**15._r8) * exp_fac(:)
        rxt(:,k,pir_ndx(10)) = 1.35e-28_r8 * tp(:)**2.83_r8

        rxt(:,k,pir_ndx(11)) = 1.e-27_r8 * (308._r8 * tinv(:))**4.7_r8
        rxt(:,k,pir_ndx(12)) = rxt(:,k,pir_ndx(11))
        rxt(:,k,pir_ndx(13)) = 1.4e-29_r8 * tp(:)**4._r8

        call comp_exp( exp_fac, -3872._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(14)) = (3.4e-7_r8 * tp(:)**5._r8) * exp_fac(:)

        rxt(:,k,pir_ndx(15)) = 3.0e-31_r8 * tp(:)**4.3_r8
        call comp_exp( exp_fac, -2093._r8 * tinv, ncol )
        rxt(:,k,pir_ndx(16)) = (1.5e-8_r8 * tp(:)**4.3_r8) * exp_fac(:)

        rxt(:,k,edn_ndx(1)) = 3.1e-10_r8 * tp(:)**0.83_r8
        call comp_exp( exp_fac, -4990._r8 * tinv, ncol )
        rxt(:,k,edn_ndx(2)) = (1.9e-12_r8 * tp(:)**(-1.5_r8)) * exp_fac(:)

        rxt(:,k,nir_ndx(1)) = 1.05e-12_r8 * tp(:)**2.15_r8
        rxt(:,k,nir_ndx(2)) = 2.5e-11_r8 * tp(:)**0.79_r8
        rxt(:,k,nir_ndx(3)) = 7.5e-11_r8 * tp(:)**0.79_r8
        rxt(:,k,nir_ndx(4)) = rxt(:,k,nir_ndx(1))
        rxt(:,k,nir_ndx(5)) = 1.3e-11_r8 * tp(:)**1.64_r8
        rxt(:,k,nir_ndx(6)) = 3.3e-11_r8 * tp(:)**2.38_r8

        call comp_exp( exp_fac, -7300_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(7)) = (1.0e-3_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -7050_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(8)) = (7.2e-4_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -6800_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(9)) = (6.5e-3_r8 * tp(:)) * exp_fac(:)
        call comp_exp( exp_fac, -7600_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(10)) = (5.7e-4_r8 * tp(:)) * exp_fac(:)

        call comp_exp( exp_fac, -7150_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(11)) = (1.5e-2_r8 * tp(:)) * exp_fac(:)

        call comp_exp( exp_fac, -13130_r8 * tinv, ncol )
        rxt(:,k,nir_ndx(12)) = (6.0e-3_r8 * tp(:)) * exp_fac(:)
        rxt(:,k,nir_ndx(13)) = 5.22e-28_r8 * tp(:)**2.62_r8

        rxt(:,k,iira_ndx(1)) = 6.0e-8_r8 * tp(:)**.5_r8
        do i = 2,niira
          rxt(:,k,iira_ndx(i)) = rxt(:,k,iira_ndx(i-1))
        enddo

        rxt(:,k,iirb_ndx(1)) = 1.25e-25_r8 * tp(:)**4._r8
        do i = 2,niirb
          rxt(:,k,iirb_ndx(i)) = rxt(:,k,iirb_ndx(i-1))
        enddo

        call comp_exp( exp_fac, -6600._r8 * tinv, ncol )
        rxt(:,k,usr_clm_h2o_m_ndx) = 2.e-8_r8 * exp_fac(:)

        call comp_exp( exp_fac, -11926._r8 * tinv, ncol )
        rxt(:,k,usr_clm_hcl_m_ndx) =  tinv(:) * exp_fac(:)

     endif
    end do level_loop

!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
    if ( has_ion_rxts ) then
       level_loop2 : do k = 1,pver
 	   tp(:ncol)         = (2._r8*tempi(:ncol,k) + temp(:ncol,k)) / ( 3._r8 * t0 )
	   tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
	   rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
			 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
	   tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*temp(:ncol,k)) / t0
	   tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
	   rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
	   tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + temp(:ncol,k))
	   where( tp(:ncol) < trlim3 )
		  rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
		  rxt(:,k,ion11_ndx) = 1.e-11_r8 * tp(:)**.23_r8
           elsewhere
		  rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
	      rxt(:,k,ion11_ndx) = 3.6e-12_r8 / tp(:)**.41_r8
	   end where
	   tp(:ncol)          = t0 / tempe(:ncol,k)
	   rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
	   rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
	   where( tp(:ncol) < 4._r8 )
	      rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
	   elsewhere
	      rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
	   end where
	end do level_loop2
     endif

     ! quenching of O+(2P) and O+(2D) by e to produce O+
     ! See TABLE 1 of Roble (1995)
     ! drm 2015-07-27
     if (elec4_ndx > 0 .and. elec5_ndx > 0 .and. elec6_ndx > 0) then
         do k=1,pver
            tp(:ncol)          = sqrt(300._r8 / tempe(:ncol,k))
            rxt(:,k,elec4_ndx) = 1.5e-7_r8 * tp(:)
            rxt(:,k,elec5_ndx) = 4.0e-8_r8 * tp(:)
            rxt(:,k,elec6_ndx) = 6.6e-8_r8 * tp(:)
         end do
     endif

!-----------------------------------------------------------------
!	... tropospheric "aerosol" rate constants
!-----------------------------------------------------------------
     if ( het1_ndx > 0 .AND. (.NOT. usr_N2O5_aer_ndx > 0) ) then
         amas = 4._r8*pi*rm1**3*den/3._r8            ! each mean particle(r=0.1u) mass (g)
         do k = 1,pver
!-------------------------------------------------------------------------
! 	... estimate humidity effect on aerosols (from Shettle and Fenn, 1979)
!           xr is a factor of the increase aerosol radii with hum (hum=0., factor=1)
!-------------------------------------------------------------------------
            xr(:)     = .999151_r8 + relhum(:ncol,k)*(1.90445_r8 + relhum(:ncol,k)*(-6.35204_r8 + relhum(:ncol,k)*5.32061_r8))
!-------------------------------------------------------------------------
! 	... estimate sulfate particles surface area (cm2/cm3) in each grid
!-------------------------------------------------------------------------
            if ( carma_hetchem_feedback ) then
               sur(:ncol) = strato_sad(:ncol,k)
            else
               sur(:) = sulfate(:,k)*m(:,k)/avo*wso4 &              ! xform mixing ratio to g/cm3
                        / amas &                                    ! xform g/cm3 to num particles/cm3
                        * fare &                                    ! xform num particles/cm3 to cm2/cm3
                        * xr(:)*xr(:)                               ! humidity factor
            endif
!-----------------------------------------------------------------
!	... compute the "aerosol" reaction rates
!-----------------------------------------------------------------
!             k = gam * A * velo/4
!
!       where velo = sqrt[ 8*bk*T/pi/(w/av) ]
!             bk = 1.381e-16
!             av = 6.02e23
!             w  = 108 (n2o5)  HO2(33)  CH2O (30)  NH3(15)
!
!       so that velo = 1.40e3*sqrt(T)  (n2o5)   gama=0.1
!       so that velo = 2.53e3*sqrt(T)  (HO2)    gama>0.2
!       so that velo = 2.65e3*sqrt(T)  (CH2O)   gama>0.022
!       so that velo = 3.75e3*sqrt(T)  (NH3)    gama=0.4
!--------------------------------------------------------
!-----------------------------------------------------------------
!	... use this n2o5 -> 2*hno3 only in troposphere
!-----------------------------------------------------------------
	    rxt(:,k,het1_ndx) = rxt(:,k,het1_ndx) &
                                +.25_r8 * gam1 * sur(:) * 1.40e3_r8 * sqrt( temp(:ncol,k) )
         end do
      end if

!lke++
!-----------------------------------------------------------------
!      ... CO tags
!-----------------------------------------------------------------
      if( usr_CO_OH_b_ndx > 0 .and. usr_CO_OH_ndx < 0 ) then
         usr_CO_OH_ndx = usr_CO_OH_b_ndx
      end if
      if( usr_CO_OH_ndx > 0 ) then
         if( usr_COhc_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COhc_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_COme_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_COme_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO01_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO01_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO02_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO02_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO03_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO03_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO04_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO04_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO05_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO05_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO06_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO06_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO07_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO07_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO08_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO08_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO09_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO09_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO10_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO10_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO11_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO11_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO12_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO12_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO13_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO13_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO14_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO14_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO15_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO15_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO16_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO16_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO17_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO17_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO18_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO18_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO19_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO19_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO20_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO20_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO21_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO21_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO22_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO22_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO23_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO23_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO24_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO24_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO25_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO25_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO26_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO26_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO27_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO27_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO28_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO28_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO29_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO29_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO30_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO30_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO31_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO31_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO32_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO32_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO33_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO33_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO34_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO34_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO35_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO35_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO36_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO36_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO37_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO37_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO38_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO38_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO39_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO39_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO40_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO40_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO41_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO41_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
         if( usr_CO42_OH_ndx > 0 ) then
            rxt(:ncol,:,usr_CO42_OH_ndx) = rxt(:ncol,:,usr_CO_OH_ndx)
         end if
      end if
!lke--
!
! jfl : additional BAM removal reactions.  Zero out below the tropopause
!
      do l=1,num_strat_tau
!
         if ( usr_strat_tau_ndx(l) > 0 ) then
            do i=1,ncol
               rxt(i,tropchemlev(i)+1:pver,usr_strat_tau_ndx(l)) = 0._r8
            end do
         end if
!
      end do
!

      deallocate( sfc_array, dm_array )

  end subroutine usrrxt

      subroutine usrrxt_hrates( rxt, tempn, tempi, tempe, &
				h2ovmr, m, ncol, kbot )
!-----------------------------------------------------------------
!        ... set the user specified reaction rates for heating
!-----------------------------------------------------------------

      use shr_kind_mod,  only : r8 => shr_kind_r8
      use chem_mods,     only : rxntot
      use ppgrid,        only : pver, pcols

      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      integer, intent(in)     :: ncol                         ! number columns in chunk
      integer, intent(in)     :: kbot                         ! heating levels
      real(r8), intent(in)    :: tempn(pcols,pver)            ! neutral temperature (K)
      real(r8), intent(in)    :: tempi(pcols,pver)            ! ion temperature (K)
      real(r8), intent(in)    :: tempe(pcols,pver)            ! electron temperature (K)
      real(r8), intent(in)    :: m(ncol,pver)                 ! total atm density (1/cm^3)
      real(r8), intent(in)    :: h2ovmr(ncol,pver)            ! water vapor (vmr)
      real(r8), intent(inout) :: rxt(ncol,pver,rxntot)        ! gas phase rates

!-----------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------

      integer  ::  k
      real(r8), dimension(ncol) :: &
                   tp, &
                   tinv, &
                   ko, &
                   kinf, &
                   fc

!-----------------------------------------------------------------
!	... o + o2 + m --> o3 + m
!-----------------------------------------------------------------
      do k = 1,kbot
         tinv(:ncol)       = 1._r8 / tempn(:ncol,k)
         tp(:)             = 300._r8 * tinv(:)
         rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8

!-----------------------------------------------------------------
!	... o + o + m -> o2 + m
!-----------------------------------------------------------------
         rxt(:,k,usr_O_O_ndx) = 2.76e-34_r8 * exp( 720.0_r8*tinv(:) )

!-----------------------------------------------------------------
!	... ho2 + ho2 --> h2o2
!	Note: this rate involves the water vapor number density
!-----------------------------------------------------------------
         ko(:)   = 3.0e-13_r8  * exp( 460._r8*tinv(:) )
         kinf(:) = 2.1e-33_r8 * m(:,k) * exp( 920._r8*tinv(:) )
         fc(:)   = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp( 2200._r8*tinv(:) )
         rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)

      end do

!-----------------------------------------------------------------
! 	... the ionic rates
!-----------------------------------------------------------------
      if ( has_ion_rxts ) then
         level_loop2 :  do k = 1,kbot
            tp(:ncol)         = (2._r8*tempi(:ncol,k) + tempn(:ncol,k)) / ( 3._r8 * t0 )
            tp(:)             = max( min( tp(:),20._r8 ),1._r8 )
            rxt(:,k,ion1_ndx) = 2.82e-11_r8 + tp(:)*(-7.74e-12_r8 + tp(:)*(1.073e-12_r8  &
                 + tp(:)*(-5.17e-14_r8 + 9.65e-16_r8*tp(:))))
            tp(:ncol)         = (.6363_r8*tempi(:ncol,k) + .3637_r8*tempn(:ncol,k)) / t0
            tp(:)             = max( min( tp(:),trlim2 ),1._r8 )
            rxt(:,k,ion2_ndx) = 1.533e-12_r8 + tp(:)*(-5.92e-13_r8 + tp(:)*8.6e-14_r8)
            tp(:ncol)         = 2._r8 * t0 /(tempi(:ncol,k) + tempn(:ncol,k))
            where( tp(:ncol) < trlim3 )
               rxt(:,k,ion3_ndx)  = 1.4e-10_r8 * tp(:)**.44_r8
            elsewhere
               rxt(:,k,ion3_ndx)  = 5.2e-11_r8 / tp(:)**.2_r8
            endwhere
            tp(:ncol)          = t0 / tempe(:ncol,k)
            rxt(:,k,elec1_ndx) = 4.e-7_r8 * tp(:)**.85_r8
            rxt(:,k,elec3_ndx) = 1.8e-7_r8 * tp(:)**.39_r8
            where( tp(:ncol) < 4._r8 )
               rxt(:,k,elec2_ndx) = 2.7e-7_r8 * tp(:)**.7_r8
            elsewhere
               rxt(:,k,elec2_ndx) = 1.6e-7_r8 * tp(:)**.55_r8
            endwhere
         end do level_loop2
      endif
      end subroutine usrrxt_hrates

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  subroutine comp_exp( x, y, n )

    implicit none

    real(r8), intent(out) :: x(:)
    real(r8), intent(in)  :: y(:)
    integer,  intent(in)  :: n

#ifdef IBM
    call vexp( x, y, n )
#else
    x(:n) = exp( y(:n) )
#endif

  end subroutine comp_exp

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: dm_aer(:)
    real(r8), intent(in) :: dg_gas
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) = sfc(i) / (0.5_r8*dm_aer(i)/dg_gas + (4._r8/(c_gas*gamma_gas)))
    enddo

    rate = sum(rxt)

    deallocate(rxt)

  endfunction hetrxtrate

  !-------------------------------------------------------------------------
  !  Heterogeneous reaction rates for uptake of a glyoxal gas on an aerosol:
  !-------------------------------------------------------------------------
  function hetrxtrate_gly( sfc, c_gas, gamma_gas ) result(rate)

    real(r8), intent(in) :: sfc(:)
    real(r8), intent(in) :: c_gas
    real(r8), intent(in) :: gamma_gas
    real(r8) :: rate

    real(r8),allocatable :: rxt(:)
    integer :: n, i

    n = size(sfc)

    allocate(rxt(n))
    do i=1,n
       rxt(i) =  0.25_r8 * c_gas * sfc(i) * gamma_gas
    enddo

    rate = sum(rxt)

    deallocate(rxt)

  endfunction hetrxtrate_gly


end module mo_usrrxt
