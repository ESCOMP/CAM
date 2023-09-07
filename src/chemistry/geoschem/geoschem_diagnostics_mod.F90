MODULE GeosChem_Diagnostics_Mod

  ! CAM modules
  use cam_history,    only : fieldname_len
  use cam_logfile,    only : iulog
  use chem_mods,      only : gas_pcnst, map2chm, iFirstCnst
  use constituents,   only : pcnst
  use mo_tracname,    only : solsym
  use ppgrid,         only : begchunk, pver
  use shr_const_mod,  only : pi => shr_const_pi
  use shr_kind_mod,   only : r8 => shr_kind_r8, shr_kind_cl
  use spmd_utils,     only : MasterProc
  use string_utils,   only : to_upper

  ! GEOS-Chem modules
  use ErrCode_Mod,    only : GC_SUCCESS

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GC_Diagnostics_Init
  PUBLIC :: GC_Diagnostics_Calc
  PUBLIC :: wetdep_name, wtrate_name, dtchem_name

  CHARACTER(LEN=fieldname_len) :: srcnam(gas_pcnst)      ! Names of source/sink tendencies
  CHARACTER(LEN=fieldname_len) :: wetdep_name(gas_pcnst) ! Wet deposition tendencies
  CHARACTER(LEN=fieldname_len) :: wtrate_name(gas_pcnst) ! Column tendencies for wet dep
  CHARACTER(LEN=fieldname_len) :: dtchem_name(gas_pcnst) ! Chemical tendencies

  INTEGER :: aer_species(gas_pcnst)

  ! Chemical families
  INTEGER  :: NOx_species(3)
  INTEGER  :: NOy_species(62)
  INTEGER  :: HOx_species(4)
  INTEGER  :: ClOx_species(6)
  INTEGER  :: ClOy_species(11)
  INTEGER  :: tCly_species(30)
  INTEGER  :: BrOx_species(4)
  INTEGER  :: BrOy_species(9)
  INTEGER  :: tBry_species(18)
  INTEGER  :: SOx_species(2)
  INTEGER  :: NHx_species(2)
  INTEGER  :: TOTH_species(3)
  REAL(r8) :: NOx_MWs(3)
  REAL(r8) :: NOy_MWs(62)
  REAL(r8) :: HOx_MWs(4)
  REAL(r8) :: ClOx_MWs(6)
  REAL(r8) :: ClOy_MWs(11)
  REAL(r8) :: tCly_MWs(30)
  REAL(r8) :: BrOx_MWs(4)
  REAL(r8) :: BrOy_MWs(9)
  REAL(r8) :: tBry_MWs(18)
  REAL(r8) :: SOx_MWs(2)
  REAL(r8) :: NHx_MWs(2)
  REAL(r8) :: TOTH_MWs(3)

  REAL(r8), PARAMETER :: MW_NIT  = 62.01
  REAL(r8), PARAMETER :: MW_HNO3 = 63.01
  REAL(r8), PARAMETER :: MW_HCl  = 36.45
  REAL(r8), PARAMETER :: MW_H2O  = 18.02

  ! NOx species
  INTEGER :: i_NO, i_NO2, i_N
  ! NOy \ NOx species
  INTEGER :: i_BrNO2, i_BrNO3, i_ClNO2, i_ClNO3, i_ETHLN, i_ETNO3,      &
             i_HNO2, i_HNO3, i_HNO4, i_ICN, i_ICNOO, i_IDHNBOO,         &
             i_IDHNDOO1, i_IDN, i_IDNOO, i_IHN1, i_IHN2,                &
             i_IHN3, i_IHN4, i_IHPNBOO, i_IHPNDOO, i_INA, i_INO,        &
             i_INO2B, i_INO2D, i_INPB, i_INPD, i_IONO, i_IONO2,         &
             i_IPRNO3, i_ISOPNOO1, i_ISOPNOO2, i_ITCN, i_ITHN,          &
             i_MACRNO2, i_MCRHN, i_MCRHNB, i_MENO3, i_MONITS, i_MONITU, &
             i_MPAN, i_MPN, i_MVKN, i_N2O5, i_NO3, i_NPRNO3, i_OLND,    &
             i_OLNN, i_PAN, i_PPN, i_PRN1, i_PROPNN, i_PRPN, i_R4N1,    &
             i_R4N2, i_HONIT, i_IONITA, i_NIT, i_NITs
  ! HOx
  INTEGER :: i_H, i_OH, i_HO2, i_H2O2
  ! ClOx
  INTEGER :: i_Cl, i_ClO, i_HOCl, i_Cl2, i_Cl2O2, i_OClO
  ! tCly \ ClOx
  INTEGER :: i_ClOO, i_HCl, i_BrCl, i_ICl, i_H1211,                   &
             i_CFC115, i_CH3Cl, i_HCFC142b, i_HCFC22, i_CH2ICl,       &
             i_CFC114, i_CFC12, i_HCFC141b, i_HCFC123, i_CH2Cl2,      &
             i_CFC11, i_CH3CCl3, i_CHCl3, i_CCl4, i_CFC113, i_SALACL, &
             i_SALCCL !ClNO2, ClNO3 already defined in NOy_species
  ! BrOx
  INTEGER :: i_Br, i_BrO, i_HOBr !BrCl already defined in tCly_species
  ! Bry \ BrOx
  INTEGER :: i_HBr, i_IBr, i_Br2, i_CH3Br,                            &
             i_H1301, i_H2402, i_CH2Br2, i_CHBr3, i_BrSALA, i_BrSALC, &
             i_CH2IBr
             !BrNO2, BrNO3 already defined in NOy_speies
             !H1211 already defined in tCly_species
  ! SOx
  INTEGER :: i_SO2, i_SO4
  ! NHx
  INTEGER :: i_NH3, i_NH4
  ! TOTH
  INTEGER :: i_CH4, i_H2O, i_H2


  ! Index in solsym
  integer :: id_no,id_no3
  integer :: id_cfc11,id_cfc12
  integer :: id_ch4,id_h2o
  integer :: id_o,id_o2,id_h,id_n2o
  integer :: id_co2,id_o3,id_oh,id_ho2,id_so4_a1,id_so4_a2,id_so4_a3
  integer :: id_num_a2,id_num_a3,id_dst_a3,id_ncl_a3

! !REVISION HISTORY:
!  28 Oct 2020 - T. M. Fritz   - Initial version

CONTAINS

  SUBROUTINE GC_Diagnostics_Init( Input_Opt, State_Chm, State_Met )

    ! CAM modules 
    use cam_history,         only : addfld, add_default, horiz_only
    use constituents,        only : cnst_name, sflxnam, cnst_get_ind
    use mo_chem_utls,        only : get_spc_ndx
    use phys_control,        only : phys_getopts

    ! GEOS-Chem modules
    use Input_Opt_Mod,       only : OptInput
    use State_Chm_Mod,       only : ChmState
    use State_Met_Mod,       only : MetState
    use State_Diag_Mod,      only : get_TagInfo
    use Species_Mod,         only : Species
    use Registry_Mod,        only : MetaRegItem, RegItem
    use State_Chm_Mod,       only : Ind_
    use DryDep_Mod,          only : depName
    
    TYPE(OptInput),    INTENT(IN) :: Input_Opt   ! Input options
    TYPE(ChmState),    INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(MetState),    INTENT(IN) :: State_Met   ! Meteorology State object
    
    INTEGER                :: M, N, K, SM
    INTEGER                :: idx
    INTEGER                :: RC
    INTEGER                :: bulkaero_species(20)
    INTEGER                :: id_so4, id_nh4no3
    INTEGER                :: id_dst01, id_dst02, id_dst03, id_dst04
    INTEGER                :: id_sslt01, id_sslt02, id_sslt03, id_sslt04
    INTEGER                :: id_soa,  id_oc1, id_oc2, id_cb1, id_cb2
    INTEGER                :: id_soam,id_soai,id_soat,id_soab,id_soax
    INTEGER                :: id_bry, id_cly 
    INTEGER                :: history_budget_histfile_num    ! output history file number
                                                             ! for budget fields

    LOGICAL                :: Found
    LOGICAL                :: history_aerosol                ! Output the MAM aerosol
                                                             ! tendencies
    LOGICAL                :: history_chemistry
    LOGICAL                :: history_cesm_forcing
    LOGICAL                :: history_scwaccm_forcing
    LOGICAL                :: history_chemspecies_srf        ! Output the chemistry
                                                             ! constituents species
                                                             ! in the surface layer
    LOGICAL                :: history_dust
    LOGICAL                :: history_budget                 ! output tendencies and state
                                                             ! variables for CAM
                                                             ! temperature, water vapor, 
                                                             ! cloud ice and cloud
                                                             ! liquid budgets.

    CHARACTER(LEN=shr_kind_cl) :: SpcName
    CHARACTER(LEN=shr_kind_cl) :: tagName
    CHARACTER(LEN=shr_kind_cl) :: ThisLoc
    CHARACTER(LEN=shr_kind_cl) :: ErrMsg
    CHARACTER(LEN=2)           :: unit_basename  ! Units 'kg' or '1'

    ! Objects
    TYPE(Species),     POINTER :: SpcInfo
    TYPE(MetaRegItem), POINTER :: Current
    TYPE(RegItem    ), POINTER :: Item

    !=================================================================
    ! GC_Diagnostics_Init begins here!
    !=================================================================

    ! Initialize pointers
    SpcInfo => NULL()
    Current => NULL()
    Item    => NULL()

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    CALL phys_getopts( history_aerosol_out             = history_aerosol,             &
                       history_chemistry_out           = history_chemistry,           &
                       history_chemspecies_srf_out     = history_chemspecies_srf,     &
                       history_budget_out              = history_budget ,             &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_cesm_forcing_out        = history_cesm_forcing,        &
                       history_scwaccm_forcing_out     = history_scwaccm_forcing,     &
                       history_dust_out                = history_dust )

    id_no3     = get_spc_ndx( 'NO3',    ignore_case=.true. )
    id_o3      = get_spc_ndx( 'O3',     ignore_case=.true. )
    id_oh      = get_spc_ndx( 'OH',     ignore_case=.true. )
    id_ho2     = get_spc_ndx( 'HO2',    ignore_case=.true. )
    id_so4_a1  = get_spc_ndx( 'so4_a1', ignore_case=.true. )
    id_so4_a2  = get_spc_ndx( 'so4_a2', ignore_case=.true. )
    id_so4_a3  = get_spc_ndx( 'so4_a3', ignore_case=.true. )
    id_num_a2  = get_spc_ndx( 'num_a2', ignore_case=.true. )
    id_num_a3  = get_spc_ndx( 'num_a3', ignore_case=.true. )
    id_dst_a3  = get_spc_ndx( 'dst_a3', ignore_case=.true. )
    id_ncl_a3  = get_spc_ndx( 'ncl_a3', ignore_case=.true. )
    id_co2     = get_spc_ndx( 'CO2',    ignore_case=.true. )
    id_no      = get_spc_ndx( 'NO',     ignore_case=.true. )
    id_h       = get_spc_ndx( 'H',      ignore_case=.true. )
    id_o       = get_spc_ndx( 'O',      ignore_case=.true. )
    id_o2      = get_spc_ndx( 'O2',     ignore_case=.true. )
    id_ch4     = get_spc_ndx( 'CH4',    ignore_case=.true. )
    id_h2o     = get_spc_ndx( 'H2O',    ignore_case=.true. )
    id_n2o     = get_spc_ndx( 'N2O',    ignore_case=.true. )
    id_cfc11   = get_spc_ndx( 'CFC11',  ignore_case=.true. )
    id_cfc12   = get_spc_ndx( 'CFC12',  ignore_case=.true. )

    id_bry     = get_spc_ndx( 'BRY',    ignore_case=.true. )
    id_cly     = get_spc_ndx( 'CLY',    ignore_case=.true. )

    id_dst01   = get_spc_ndx( 'DST01',  ignore_case=.true. )
    id_dst02   = get_spc_ndx( 'DST02',  ignore_case=.true. )
    id_dst03   = get_spc_ndx( 'DST03',  ignore_case=.true. )
    id_dst04   = get_spc_ndx( 'DST04',  ignore_case=.true. )
    id_sslt01  = get_spc_ndx( 'SSLT01', ignore_case=.true. )
    id_sslt02  = get_spc_ndx( 'SSLT02', ignore_case=.true. )
    id_sslt03  = get_spc_ndx( 'SSLT03', ignore_case=.true. )
    id_sslt04  = get_spc_ndx( 'SSLT04', ignore_case=.true. )
    id_soa     = get_spc_ndx( 'SOA',    ignore_case=.true. )
    !id_so4     = get_spc_ndx( 'SO4',    ignore_case=.true. )i
    id_so4 = -1 ! Don't pick up GEOS-Chem's SO4!
    id_oc1     = get_spc_ndx( 'OC1',    ignore_case=.true. )
    id_oc2     = get_spc_ndx( 'OC2',    ignore_case=.true. )
    id_cb1     = get_spc_ndx( 'CB1',    ignore_case=.true. )
    id_cb2     = get_spc_ndx( 'CB2',    ignore_case=.true. )
    id_nh4no3  = get_spc_ndx( 'NH4NO3', ignore_case=.true. )
    id_soam    = get_spc_ndx( 'SOAM',   ignore_case=.true. )
    id_soai    = get_spc_ndx( 'SOAI',   ignore_case=.true. )
    id_soat    = get_spc_ndx( 'SOAT',   ignore_case=.true. )
    id_soab    = get_spc_ndx( 'SOAB',   ignore_case=.true. )
    id_soax    = get_spc_ndx( 'SOAX',   ignore_case=.true. )

    bulkaero_species(:) = -1
    bulkaero_species(1:20) = (/ id_dst01, id_dst02, id_dst03, id_dst04, &
                                id_sslt01, id_sslt02, id_sslt03, id_sslt04, &
                                id_soa, id_so4, id_oc1, id_oc2, id_cb1, id_cb2, id_nh4no3, &
                                id_soam,id_soai,id_soat,id_soab,id_soax /)
    aer_species(:) = -1
    n = 1
    do m = 1,gas_pcnst
       k=0
       if ( any(bulkaero_species(:)==m) ) k=1
       if ( k==0 ) k = index(trim(solsym(m)), '_a')
       if ( k==0 ) k = index(trim(solsym(m)), '_c')
       if ( k>0 ) then ! must be aerosol species
          aer_species(n) = m
          n = n+1
       endif
    enddo

    CALL Addfld( 'MASS', (/ 'lev' /), 'A', 'kg', 'Mass of grid box' )
    CALL Addfld( 'AREA', horiz_only,  'A', 'm2', 'Area of grid box' )
    CALL Addfld( 'HEIGHT', (/ 'ilev' /),'A','m', 'Geopotential height above surface at interfaces' )

    ! Note that constituents are already output by default
    ! Add all species as output fields if desired
    DO N = 1, gas_pcnst
       IF ( ANY( aer_species == N ) ) THEN
          SpcName = TRIM(solsym(N))
          unit_basename = 'kg'
          IF ( SpcName(1:3) == 'num' ) unit_basename = ' 1'
          CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', unit_basename//'/kg', &
             TRIM(SpcName)//' concentration' )
          CALL AddFld( TRIM(SpcName)//'_SRF', horiz_only, 'A', unit_basename//'/kg', &
             TRIM(SpcName)//' in bottom layer' )
       ELSE
          M = map2chm(N)
          SpcName = TRIM(solsym(N))
          CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
             TRIM(SpcName)//' volume mixing ratio')
          CALL AddFld( TRIM(SpcName)//'_SRF', horiz_only, 'A', 'mol/mol', &
             TRIM(SpcName)//' in bottom layer')
       ENDIF
       IF ( ( N /= id_cly ) .AND. ( N /= id_bry ) ) THEN
          IF ( history_aerosol .OR. history_chemistry ) THEN
             CALL Add_Default( TRIM(SpcName), 1, ' ' )
          ENDIF
          IF ( history_chemspecies_srf ) THEN
             CALL Add_Default( TRIM(SpcName)//'_SRF', 1, ' ' )
          ENDIF
       ENDIF

       IF ( history_cesm_forcing ) THEN
          IF ( N == id_o3     ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_oh     ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_no3    ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_ho2    ) CALL Add_Default( TRIM(SpcName), 1, ' ')

          IF ( N == id_o3     ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_so4_a1 ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_so4_a2 ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_so4_a3 ) CALL Add_Default( TRIM(SpcName), 8, ' ')

          IF ( N == id_num_a2 ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_num_a3 ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_dst_a3 ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_ncl_a3 ) CALL Add_Default( TRIM(SpcName), 8, ' ')

       ENDIF
       IF ( history_scwaccm_forcing ) THEN
          IF ( N == id_co2   ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_h     ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_no    ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_o     ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_o2    ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_o3    ) CALL Add_Default( TRIM(SpcName), 8, ' ')
          IF ( N == id_h2o   ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_ch4   ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_n2o   ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_cfc11 ) CALL Add_Default( TRIM(SpcName), 1, ' ')
          IF ( N == id_cfc12 ) CALL Add_Default( TRIM(SpcName), 1, ' ')
       ENDIF

       IF (history_dust .AND. (index(TRIM(SpcName),'dst_') > 0)) THEN
          CALL Add_Default( TRIM(SpcName), 1, ' ')
       ENDIF
    ENDDO

    IF ( Input_Opt%LDryD ) THEN
       DO N = 1, State_Chm%nDryDep
          SpcName = 'DV_'//to_upper(TRIM(depName(N)))
          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'm/s', &
            TRIM(SpcName)//' dry deposition velocity')
       ENDDO

       DO N = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          M = State_Chm%Map_Advect(N)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(M)%Info
          SpcName = 'DF_'//to_upper(TRIM(SpcInfo%Name))
          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
             TRIM(SpcName)//' dry deposition flux')
          IF ( history_chemistry ) THEN
             CALL Add_Default( TRIM(SpcName), 1, ' ' )
          ENDIF

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! Chemical tendencies and surface fluxes
    DO N = 1, gas_pcnst
       IF ( map2chm(N) > 0 ) THEN
          ! If this is a GEOS-Chem species then capitalize. This avoids
          ! issues where Br2 /= BR2
          srcnam(N) = 'CT_'//to_upper(TRIM(solsym(N))) ! chem tendency (source/sink)
       ELSE
          ! For MAM aerosols, keep as it is (i.e. bc_a1)
          srcnam(N) = 'CT_'//TRIM(solsym(N)) ! chem tendency (source/sink)
       ENDIF
       SpcName = srcnam(N)
       CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg/s', TRIM(SpcName)//' source/sink' )

       SpcName = TRIM(solsym(N))
       CALL cnst_get_ind( SpcName, M, abort=.false. )
       IF ( M > 0 ) THEN
          IF (sflxnam(M)(3:5) == 'num') THEN  ! name is in the form of "SF****"
             unit_basename = ' 1'
          ELSE
             unit_basename = 'kg'
          ENDIF
          SpcName = sflxnam(M)
          CALL Addfld ( TRIM(SpcName), horiz_only, 'A',  unit_basename//'/m2/s', &
             TRIM(solsym(N))//' surface flux')
          IF ( history_aerosol .OR. history_chemistry ) THEN 
             CALL Add_Default( TRIM(SpcName), 1, ' ' )
          ENDIF

          IF ( history_cesm_forcing ) THEN
             IF ( TRIM(SpcName(3:)) == 'NO' .OR. TRIM(SpcName(3:)) == 'NH3' ) THEN
                CALL Add_Default( TRIM(SpcName), 1, ' ' )
             ENDIF
          ENDIF
       ENDIF
    ENDDO

    ! Add chemical tendency of water vapor to water budget output
    IF ( history_budget ) THEN 
      CALL Add_Default ('CT_H2O' , history_budget_histfile_num, ' ')
    ENDIF

    ! Chemical tendencies
    DO N = 1, gas_pcnst
       M = map2chm(N)
       IF ( M > 0 ) THEN
          dtchem_name(N) = 'D'//to_upper(TRIM(solsym(N)))//'CHM'
       ELSE
          dtchem_name(N) = 'D'//TRIM(solsym(N))//'CHM'
       ENDIF
       SpcName = TRIM(dtchem_name(N))
       CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/s', &
          'net tendency from chemistry' )
    ENDDO

    i_NO          = Ind_('NO')
    i_NO2         = Ind_('NO2')
    i_N           = Ind_('N')
    i_BrNO2       = Ind_('BrNO2')
    i_BrNO3       = Ind_('BrNO3')
    i_ClNO2       = Ind_('ClNO2')
    i_ClNO3       = Ind_('ClNO3')
    i_ETHLN       = Ind_('ETHLN')
    i_ETNO3       = Ind_('ETNO3')
    i_HNO2        = Ind_('HNO2')
    i_HNO3        = Ind_('HNO3')
    i_HNO4        = Ind_('HNO4')
    i_ICN         = Ind_('ICN')
    i_ICNOO       = Ind_('ICNOO')
    i_IDHNBOO     = Ind_('IDHNBOO')
    i_IDHNDOO1    = Ind_('IDHNDOO1')
    i_IDN         = Ind_('IDN')
    i_IDNOO       = Ind_('IDNOO')
    i_IHN1        = Ind_('IHN1')
    i_IHN2        = Ind_('IHN2')
    i_IHN3        = Ind_('IHN3')
    i_IHN4        = Ind_('IHN4')
    i_IHPNBOO     = Ind_('IHPNBOO')
    i_IHPNDOO     = Ind_('IHPNDOO')
    i_INA         = Ind_('INA')
    i_INO         = Ind_('INO')
    i_INO2B       = Ind_('INO2B')
    i_INO2D       = Ind_('INO2D')
    i_INPB        = Ind_('INPB')
    i_INPD        = Ind_('INPD')
    i_IONO        = Ind_('IONO')
    i_IONO2       = Ind_('IONO2')
    i_IPRNO3      = Ind_('IPRNO3')
    i_ISOPNOO1    = Ind_('ISOPNOO1')
    i_ISOPNOO2    = Ind_('ISOPNOO2')
    i_ITCN        = Ind_('ITCN')
    i_ITHN        = Ind_('ITHN')
    i_MACRNO2     = Ind_('MACRNO2')
    i_MCRHN       = Ind_('MCRHN')
    i_MCRHNB      = Ind_('MCRHNB')
    i_MENO3       = Ind_('MENO3')
    i_MONITS      = Ind_('MONITS')
    i_MONITU      = Ind_('MONITU')
    i_MPAN        = Ind_('MPAN')
    i_MPN         = Ind_('MPN')
    i_MVKN        = Ind_('MVKN')
    i_N2O5        = Ind_('N2O5')
    i_NO3         = Ind_('NO3')
    i_NPRNO3      = Ind_('NPRNO3')
    i_OLND        = Ind_('OLND')
    i_OLNN        = Ind_('OLNN')
    i_PAN         = Ind_('PAN')
    i_PPN         = Ind_('PPN')
    i_PRN1        = Ind_('PRN1')
    i_PROPNN      = Ind_('PROPNN')
    i_PRPN        = Ind_('PRPN')
    i_R4N1        = Ind_('R4N1')
    i_R4N2        = Ind_('R4N2')
    i_HONIT       = Ind_('HONIT')
    i_IONITA      = Ind_('IONITA')
    i_NIT         = Ind_('NIT')
    i_NITs        = Ind_('NITs')
    i_H           = Ind_('H')
    i_OH          = Ind_('OH')
    i_HO2         = Ind_('HO2')
    i_H2O2        = Ind_('H2O2')
    i_Cl          = Ind_('Cl')
    i_ClO         = Ind_('ClO')
    i_HOCl        = Ind_('HOCl')
    i_Cl2         = Ind_('Cl2')
    i_Cl2O2       = Ind_('Cl2O2')
    i_OClO        = Ind_('OClO')
    i_ClOO        = Ind_('ClOO')
    i_HCl         = Ind_('HCl')
    i_ClNO2       = Ind_('ClNO2')
    i_ClNO3       = Ind_('ClNO3')
    i_BrCl        = Ind_('BrCl')
    i_ICl         = Ind_('ICl')
    i_H1211       = Ind_('H1211')
    i_CFC115      = Ind_('CFC115')
    i_CH3Cl       = Ind_('CH3Cl')
    i_HCFC142b    = Ind_('HCFC142b')
    i_HCFC22      = Ind_('HCFC22')
    i_CH2ICl      = Ind_('CH2ICl')
    i_CFC114      = Ind_('CFC114')
    i_CFC12       = Ind_('CFC12')
    i_HCFC141b    = Ind_('HCFC141b')
    i_HCFC123     = Ind_('HCFC123')
    i_CH2Cl2      = Ind_('CH2Cl2')
    i_CFC11       = Ind_('CFC11')
    i_CH3CCl3     = Ind_('CH3CCl3')
    i_CHCl3       = Ind_('CHCl3')
    i_CCl4        = Ind_('CCl4')
    i_CFC113      = Ind_('CFC113')
    i_SALACL      = Ind_('SALACL')
    i_SALCCL      = Ind_('SALCCL')
    i_Br          = Ind_('Br')
    i_BrO         = Ind_('BrO')
    i_BrCl        = Ind_('BrCl')
    i_HOBr        = Ind_('HOBr')
    i_HBr         = Ind_('HBr')
    i_BrNO2       = Ind_('BrNO2')
    i_BrNO3       = Ind_('BrNO3')
    i_IBr         = Ind_('IBr')
    i_Br2         = Ind_('Br2')
    i_CH3Br       = Ind_('CH3Br')
    i_H1211       = Ind_('H1211')
    i_H1301       = Ind_('H1301')
    i_H2402       = Ind_('H2402')
    i_CH2Br2      = Ind_('CH2Br2')
    i_CHBr3       = Ind_('CHBr3')
    i_BrSALA      = Ind_('BrSALA')
    i_BrSALC      = Ind_('BrSALC')
    i_CH2IBr      = Ind_('CH2IBr')
    i_SO2         = Ind_('SO2')
    i_SO4         = Ind_('SO4')
    i_NH3         = Ind_('NH3')
    i_NH4         = Ind_('NH4')
    i_CH4         = Ind_('CH4')
    i_H2O         = Ind_('H2O')
    i_H2          = Ind_('H2')

    NOx_species  = (/ i_N, i_NO, i_NO2 /)
    NOy_species  = (/ i_N, i_NO, i_NO2, i_BrNO2, i_BrNO3, i_ClNO2, i_ClNO3,&
                      i_ETHLN, i_ETNO3, i_HNO2, i_HNO3, i_HNO4, i_ICN,     &
                      i_ICNOO, i_IDHNBOO, i_IDHNDOO1, i_IDN,               &
                      i_IDNOO, i_IHN1, i_IHN2, i_IHN3, i_IHN4, i_IHPNBOO,  &
                      i_IHPNDOO, i_INA, i_INO, i_INO2B, i_INO2D, i_INPB,   &
                      i_INPD, i_IONO, i_IONO2, i_IPRNO3, i_ISOPNOO1,       &
                      i_ISOPNOO2, i_ITCN, i_ITHN, i_MACRNO2, i_MCRHN,      &
                      i_MCRHNB, i_MENO3, i_MONITS, i_MONITU, i_MPAN, i_MPN,&
                      i_MVKN, i_N2O5, i_NO3, i_NPRNO3, i_OLND, i_OLNN,     &
                      i_PAN, i_PPN, i_PRN1, i_PROPNN, i_PRPN, i_R4N1,      &
                      i_R4N2, i_HONIT, i_IONITA, i_NIT, i_NITs /)
    HOx_species  = (/ i_H, i_OH, i_HO2, i_H2O2 /)
    ClOx_species = (/ i_Cl, i_ClO, i_HOCl, i_Cl2, i_Cl2O2, i_OClO /)
    ClOy_species = (/ i_Cl, i_ClO, i_HOCl, i_Cl2, i_Cl2O2, i_OClO, &
                      i_HCl, i_ClNO3, i_BrCl, i_ICl, i_ClNO2 /)
    tCly_species = (/ i_Cl, i_ClO, i_HOCl, i_Cl2, i_Cl2O2, i_OClO, i_ClOO, &
                      i_HCl, i_ClNO2, i_ClNO3, i_BrCl, i_ICl, i_H1211,     &
                      i_CFC115, i_CH3Cl, i_HCFC142b, i_HCFC22, i_CH2ICl,   &
                      i_CFC114, i_CFC12, i_HCFC141b, i_HCFC123, i_CH2Cl2,  &
                      i_CFC11, i_CH3CCl3, i_CHCl3, i_CCl4, i_CFC113,       &
                      i_SALACL, i_SALCCL /)
    BrOx_species = (/ i_Br, i_BrO, i_BrCl, i_HOBr /)
    BrOy_species = (/ i_Br, i_BrO, i_BrCl, i_HOBr, i_HBr, i_BrNO2,         &
                      i_BrNO3, i_IBr, i_Br2 /)
    tBry_species = (/ i_Br, i_BrO, i_BrCl, i_HOBr, i_HBr, i_BrNO2,         &
                      i_BrNO3, i_IBr, i_Br2, i_CH3Br, i_H1211, i_H1301,    &
                      i_H2402, i_CH2Br2, i_CHBr3, i_BrSALA, i_BrSALC,      &
                      i_CH2IBr /)
    SOx_species  = (/ i_SO2, i_SO4 /)
    NHx_species  = (/ i_NH3, i_NH4 /)
    TOTH_species = (/ i_CH4, i_H2O, i_H2 /)

    DO N = 1, SIZE(NOx_species)
       idx = NOx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          NOx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          NOx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(NOy_species)
       idx = NOy_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          NOy_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          NOy_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(HOx_species)
       idx = HOx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          HOx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          HOx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(ClOx_species)
       idx = ClOx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          ClOx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          ClOx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(ClOy_species)
       idx = ClOy_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          ClOy_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          ClOy_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(tCly_species)
       idx = tCly_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          tCly_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          tCly_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(BrOx_species)
       idx = BrOx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          BrOx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          BrOx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(BrOy_species)
       idx = BrOy_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          BrOy_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          BrOy_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(tBry_species)
       idx = tBry_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          tBry_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          tBry_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(SOx_species)
       idx = SOx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          SOx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          SOx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(NHx_species)
       idx = NHx_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          NHx_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          NHx_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    DO N = 1, SIZE(TOTH_species)
       idx = TOTH_species(N)
       IF ( idx > 0 ) THEN
          SpcInfo => State_Chm%SpcData(idx)%Info
          TOTH_MWs(N) = REAL(SpcInfo%MW_g,r8)
          SpcInfo => NULL()
       ELSE
          TOTH_MWs(N) = -1.0e+00_r8
       ENDIF
    ENDDO

    IF ( ANY(NOx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "NOx indices: ", NOx_species
    ENDIF
    IF ( ANY(NOy_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "NOy indices: ", NOy_species
    ENDIF
    IF ( ANY(HOx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "HOx indices: ", HOx_species
    ENDIF
    IF ( ANY(ClOx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "ClOx indices: ", ClOx_species
    ENDIF
    IF ( ANY(ClOy_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "ClOy indices: ", ClOy_species
    ENDIF
    IF ( ANY(tCly_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "tCly indices: ", tCly_species
    ENDIF
    IF ( ANY(BrOx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "BrOx indices: ", BrOx_species
    ENDIF
    IF ( ANY(BrOy_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "BrOy indices: ", BrOy_species
    ENDIF
    IF ( ANY(tBry_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "tBry indices: ", tBry_species
    ENDIF
    IF ( ANY(SOx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "SOx indices: ", SOx_species
    ENDIF
    IF ( ANY(NHx_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "NHx indices: ", NHx_species
    ENDIF
    IF ( ANY(TOTH_species <= 0 ) ) THEN
       IF ( MasterProc ) Write(iulog,*) "TOTH indices: ", TOTH_species
    ENDIF

    CALL Addfld( 'NOX', (/ 'lev' /), 'A', 'mol/mol', &
       'NOx molar mixing ratio' )
    CALL Addfld( 'NOY', (/ 'lev' /), 'A', 'mol/mol', &
       'NOy molar mixing ratio' )
    CALL Addfld( 'NOY_mmr', (/ 'lev' /), 'A', 'kg/kg', &
       'NOy mass mixing ratio' )
    CALL Addfld( 'NOY_SRF', horiz_only, 'A', 'mol/mol', &
       'Surface NOy molar mixing ratio' )
    CALL Addfld( 'HOX', (/ 'lev' /), 'A', 'mol/mol', &
       'HOx molar mixing ratio' )
    CALL Addfld( 'CLOX', (/ 'lev' /), 'A', 'mol/mol', &
       'ClOx molar mixing ratio' )
    CALL Addfld( 'CLOY', (/ 'lev' /), 'A', 'mol/mol', &
       'Total inorganic chlorine (ClOy) molar mixing ratio' )
    CALL Addfld( 'TCLY', (/ 'lev' /), 'A', 'mol/mol', &
       'Total Cl molar mixing ratio' )
    CALL Addfld( 'BROX', (/ 'lev' /), 'A', 'mol/mol', &
       'BrOx molar mixing ratio' )
    CALL Addfld( 'BROY', (/ 'lev' /), 'A', 'mol/mol', &
       'Total inorganic bromine (BrOy) molar mixing ratio' )
    CALL Addfld( 'TBRY', (/ 'lev' /), 'A', 'mol/mol', &
       'Total Br molar mixing ratio' )
    CALL Addfld( 'SOX', (/ 'lev' /), 'A', 'mol/mol', &
       'SOx molar mixing ratio' )
    CALL Addfld( 'SOX_mmr', (/ 'lev' /), 'A', 'kg/kg', &
       'SOx mass mixing ratio' )
    CALL Addfld( 'NHX', (/ 'lev' /), 'A', 'mol/mol', &
       'NHx molar mixing ratio' )
    CALL Addfld( 'NHX_mmr', (/ 'lev' /), 'A', 'kg/kg', &
       'NHx mass mixing ratio' )
    CALL Addfld( 'TOTH', (/ 'lev' /), 'A', 'mol/mol', &
       'Total H2 molar mixing ratio' )

    CALL Addfld( 'SAD_STRAT',    (/ 'lev' /), 'I', 'cm2/cm3', 'Stratospheric aerosol SAD' )
    CALL Addfld( 'SAD_SULFC',    (/ 'lev' /), 'I', 'cm2/cm3', 'Chemical sulfate aerosol SAD' )
    CALL Addfld( 'SAD_PSC',      (/ 'lev' /), 'I', 'cm2/cm3', 'PSC aerosol SAD' )
    CALL Addfld( 'RAD_SULFC',    (/ 'lev' /), 'I', 'cm',      'Chemical sulfate radius' )
    CALL Addfld( 'RAD_PSC',      (/ 'lev' /), 'I', 'cm',      'PSC aerosol radius' )
    CALL Addfld( 'SAD_TROP',     (/ 'lev' /), 'I', 'cm2/cm3', 'Tropospheric aerosol SAD' )
    CALL Addfld( 'SAD_AERO',     (/ 'lev' /), 'I', 'cm2/cm3', 'Aerosol surface area density' )
    IF ( history_cesm_forcing ) THEN
       CALL Add_Default( 'SAD_AERO', 8, ' ' )
    ENDIF
    CALL Addfld( 'REFF_AERO',    (/ 'lev' /), 'I', 'cm',      'Aerosol effective radius')
    CALL Addfld( 'SULF_TROP',    (/ 'lev' /), 'I', 'cm2/cm3', 'Tropospheric sulfate area density')

    CALL Addfld( 'HNO3_TOTAL',   (/ 'lev' /), 'I', 'mol/mol', 'Total HNO3' )
    CALL Addfld( 'HNO3_STS',     (/ 'lev' /), 'I', 'mol/mol', 'STS condensed HNO3' )
    CALL Addfld( 'HNO3_NAT',     (/ 'lev' /), 'I', 'mol/mol', 'NAT condensed HNO3' )
    CALL Addfld( 'HNO3_GAS',     (/ 'lev' /), 'I', 'mol/mol', 'Gas phase HNO3' )
    CALL Addfld( 'H2O_GAS',      (/ 'lev' /), 'I', 'mol/mol', 'Gas phase H2O' )
    CALL Addfld( 'HCL_TOTAL',    (/ 'lev' /), 'I', 'mol/mol', 'Total HCl' )
    CALL Addfld( 'HCL_GAS',      (/ 'lev' /), 'I', 'mol/mol', 'Gas phase HCl' )
    CALL Addfld( 'HCL_STS',      (/ 'lev' /), 'I', 'mol/mol', 'STS condensend HCl' )

    CALL Addfld( 'SZA',          horiz_only,  'I', 'degrees', 'Solar Zenith Angle' )
    CALL Addfld( 'U_SRF',        horiz_only,  'I', 'm/s',     'Horizontal wind velocity' )
    CALL Addfld( 'V_SRF',        horiz_only,  'I', 'm/s',     'Vertical wind velocity' )
    CALL Addfld( 'Q_SRF',        horiz_only,  'I', 'kg/kg',   'Specific humidity' )

    CALL Addfld( 'CT_H2O_GHG',   (/ 'lev' /), 'A','kg/kg/s', 'ghg-chem h2o source/sink' )


    ! Cleanup
    Current => NULL()
    Item    => NULL()

  END SUBROUTINE GC_Diagnostics_Init
  
  SUBROUTINE GC_Diagnostics_Calc( Input_Opt,  State_Chm, State_Diag, &
                                  State_Grid, State_Met, cam_in, state, &
                                  mmr_tend,   LCHNK )

    ! CAM modules
    use cam_history,         only : outfld, hist_fld_active
    use camsrfexch,          only : cam_in_t
    use chem_mods,           only : adv_mass
    use constituents,        only : cnst_name, sflxnam
    use physconst,           only : MWDry
    use physics_types,       only : physics_state
    use spmd_utils,          only : MasterProc

    ! GEOS-Chem modules
    use CMN_Size_Mod,        only : NDUST
    use DryDep_Mod,          only : depName, Ndvzind
    use Input_Opt_Mod,       only : OptInput
    use Precision_Mod,       only : f8
    use Species_Mod,         only : Species
    use State_Chm_Mod,       only : ChmState
    use State_Diag_Mod,      only : DgnState, get_TagInfo
    use State_Grid_Mod,      only : GrdState
    use State_Met_Mod,       only : MetState
    use Registry_Mod,        only : MetaRegItem, RegItem, Registry_Lookup
    use UCX_Mod,             only : GET_STRAT_OPT

    TYPE(OptInput),      INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState),      INTENT(IN)    :: State_Diag  ! Diag State object
    TYPE(GrdState),      INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),      INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(cam_in_t),      INTENT(IN)    :: cam_in      ! import state
    TYPE(physics_state), INTENT(IN)    :: state       ! Physics state variables
    REAL(r8),            INTENT(IN)    :: mmr_tend(state%ncol,pver,gas_pcnst)
                                                      ! Net tendency from chemistry in kg/s
    INTEGER,             INTENT(IN)    :: LCHNK       ! Chunk number

    ! Integers
    INTEGER                :: I, J, L, M, N, ND, SM
    INTEGER                :: idx
    INTEGER                :: RC
    INTEGER                :: Rank              ! Size of data

    INTEGER                :: nY, nZ

    ! Logicals
    LOGICAL                :: Found
    LOGICAL                :: rootChunk
    LOGICAL                :: OnLevelEdges      ! Is the data defined
                                                !  on level edges (T/F)

    ! Strings
    CHARACTER(LEN=shr_kind_cl) :: ThisLoc
    CHARACTER(LEN=shr_kind_cl) :: ErrMsg
    CHARACTER(LEN=shr_kind_cl) :: SpcName
    CHARACTER(LEN=shr_kind_cl) :: tagName

    ! Real
    REAL(r8)               :: wgt
    REAL(r8)               :: MW
    REAL(r8)               :: RAER, REFF, SADSTRAT, XSASTRAT

    ! Arrays
    REAL(r8)               :: outTmp(State_Grid%nY,State_Grid%nZ)
    REAL(r8)               :: radTmp(State_Grid%nY,State_Grid%nZ)

    ! Floating-point data pointers (8-byte precision)
    REAL(f8),   POINTER    :: Ptr0d_8           ! 0D 8-byte data
    REAL(f8),   POINTER    :: Ptr1d_8(:    )    ! 1D 8-byte data
    REAL(f8),   POINTER    :: Ptr2d_8(:,:  )    ! 2D 8-byte data
    REAL(f8),   POINTER    :: Ptr3d_8(:,:,:)    ! 3D 8-byte data

    ! Objects
    TYPE(Species),     POINTER :: SpcInfo
    TYPE(MetaRegItem), POINTER :: Current
    TYPE(RegItem    ), POINTER :: Item

    !=================================================================
    ! GC_Diagnostics_Calc begins here!
    !=================================================================

    nY = State_Grid%nY
    nZ = State_Grid%nZ

    ! Initialize pointers
    SpcInfo => NULL()
    Current => NULL()
    Item    => NULL()

    ! For error trapping
    ErrMsg  = ''
    ThisLoc = ' -> at GC_Diagnostics_Calc (in chemistry/geoschem/geoschem_diagnostics_mod.F90)'

    ! Define rootChunk
    rootChunk = ( MasterProc.and.(LCHNK==BEGCHUNK) )

    CALL OutFld( 'AREA', State_Grid%Area_M2(1,:nY), nY, LCHNK)
    CALL OutFld( 'MASS', State_Met%AD(1,:nY,nZ:1:-1), nY, LCHNK)
    CALL Outfld( 'HEIGHT', state%zi(:nY,:), nY, LCHNK )

    ! ===============================================
    ! Diagnose chemical species (constituents and short-lived)
    ! ===============================================

    DO N = 1, gas_pcnst
       M = map2chm(N)
       SpcName = TRIM(solsym(N))
       outTmp = 0.0e+00_r8
       IF ( adv_mass(N) > 0.0e+00_r8 .AND. M /= 0 .AND. &
            (hist_fld_active(TRIM(SpcName)) .OR. hist_fld_active(TRIM(SpcName)//'_SRF')) ) THEN
          IF ( M > 0 ) THEN
             ! mol/mol
             outTmp(:nY,:) = REAL(State_Chm%Species(M)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / adv_mass(N)
          ELSEIF ( ANY( aer_species == N ) ) THEN
             ! kg/kg
             outTmp(:nY,:) = state%q(:nY,:nZ,-M)
          ELSE
             ! mol/mol
             outTmp(:nY,:) = state%q(:nY,:nZ,-M) * MWDry / adv_mass(N)
          ENDIF
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
          CALL OutFld( TRIM(SpcName)//'_SRF', outTmp(:nY,nZ), nY, LCHNK )
       ENDIF
    ENDDO

    ! ===============================================
    ! Diagnose chemical families (NOx, NOy, ...)
    ! ===============================================

    SpcName = 'NOX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(NOx_species)
          idx = NOx_species(N)
          MW  = NOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'NOY'
    IF ( hist_fld_active(TRIM(SpcName)) .OR. hist_fld_active(TRIM(SpcName)//'_SRF') ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(NOy_species)
          idx = NOy_species(N)
          MW  = NOy_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_N2O5 .OR. idx == i_IDN .OR. idx == i_IDNOO ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF
    SpcName = 'NOY_SRF'
    IF ( hist_fld_active(TRIM(SpcName)) ) CALL Outfld( TRIM(SpcName), outTmp(:nY,nZ), nY, LCHNK )

    SpcName = 'NOY_mmr'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(NOy_species)
          idx = NOy_species(N)
          MW  = NOy_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_N2O5 .OR. idx == i_IDN .OR. idx == i_IDNOO ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8)
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'HOX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(HOx_species)
          idx = HOx_species(N)
          MW  = HOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_H2O2 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'CLOX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(ClOx_species)
          idx = ClOx_species(N)
          MW  = ClOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_Cl2 .OR. idx == i_Cl2O2 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'CLOY'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(ClOy_species)
          idx = ClOy_species(N)
          MW  = ClOy_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_Cl2 .OR. idx == i_Cl2O2 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'TCLY'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(tCly_species)
          idx = tCly_species(N)
          MW  = tCly_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_Cl2 .OR. idx == i_Cl2O2 .OR. idx == i_CFC114 .OR. &
               idx == i_CFC12 .OR. idx == i_CH2Cl2 .OR. idx == i_HCFC123 .OR. &
               idx == i_HCFC141b ) THEN
             wgt = 2.0E+00_r8
          ELSEIF ( idx == i_CFC11 .OR. idx == i_CFC113 .OR. idx == i_CH3CCl3 .OR. &
                   idx == i_CHCl3 ) THEN
             wgt = 3.0E+00_r8
          ELSEIF ( idx == i_CCl4 ) THEN
             wgt = 4.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'BROX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(BrOx_species)
          idx = BrOx_species(N)
          MW  = BrOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'BROY'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(BrOy_species)
          idx = BrOy_species(N)
          MW  = BrOy_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_Br2 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'TBRY'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(tBry_species)
          idx = tBry_species(N)
          MW  = tBry_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_Br2 .OR. idx == i_H2402 .OR. idx == i_CH2Br2 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'SOX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(SOx_species)
          idx = SOx_species(N)
          MW  = SOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'SOX_mmr'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(SOx_species)
          idx = SOx_species(N)
          MW  = SOx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8)
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'NHX'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(NHx_species)
          idx = NHx_species(N)
          MW  = NHx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'NHX_mmr'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(NHx_species)
          idx = NHx_species(N)
          MW  = NHx_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8)
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    SpcName = 'TOTH'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp = 0.0e+00_r8
       DO N = 1, SIZE(TOTH_species)
          idx = TOTH_species(N)
          MW  = TOTH_MWs(N)
          IF ( idx <= 0 .OR. MW <= 0.0e+00_r8 ) CYCLE
          wgt = 1.0E+00_r8
          IF ( idx == i_CH4 ) THEN
             wgt = 2.0E+00_r8
          ENDIF
          outTmp(:nY,:) = outTmp(:nY,:) &
                        + wgt * REAL(State_Chm%Species(idx)%Conc(1,:nY,nZ:1:-1),r8) * MWDry / MW
       ENDDO
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    ! ===============================================
    ! Diagnose GEOS-Chem aerosol quantities
    ! ===============================================

    IF ( hist_fld_active('SAD_PSC') .OR. hist_fld_active('RAD_PSC') ) THEN
       outTmp = 0.0e+00_r8
       radTmp = 0.0e+00_r8
       DO J = 1, nY
       DO L = 1, nZ
          CALL GET_STRAT_OPT(State_Chm,1,J,L,1,RAER,REFF,SADSTRAT,XSASTRAT)
          outTmp(J,nZ+1-L) = SADSTRAT
          radTmp(J,nZ+1-L) = RAER
       ENDDO
       ENDDO
       CALL Outfld( 'SAD_PSC', outTmp(:nY,:), nY, LCHNK )
       CALL Outfld( 'RAD_PSC', radTmp(:nY,:), nY, LCHNK )
    ENDIF

    IF ( hist_fld_active('SAD_SULFC') .OR. hist_fld_active('RAD_SULFC') ) THEN
       outTmp = 0.0e+00_r8
       DO J = 1, nY
       DO L = 1, nZ
          CALL GET_STRAT_OPT(State_Chm,1,J,L,2,RAER,REFF,SADSTRAT,XSASTRAT)
          outTmp(J,nZ+1-L) = SADSTRAT
          radTmp(J,nZ+1-L) = RAER
       ENDDO
       ENDDO
       CALL Outfld( 'SAD_SULFC', outTmp(:nY,:), nY, LCHNK )
       CALL Outfld( 'RAD_SULFC', radTmp(:nY,:), nY, LCHNK )
    ENDIF

    IF ( hist_fld_active('SAD_AERO') .OR. hist_fld_active('SAD_TROP') ) THEN
       outTmp(:nY,:) = SUM(State_Chm%AeroArea(1,:nY,nZ:1:-1,:), DIM=3)
       CALL Outfld( 'SAD_AERO', outTmp(:nY,:), nY, LCHNK )
    ENDIF

    IF ( hist_fld_active('SAD_TROP') ) THEN
       DO J = 1, nY
       DO L = 1, nZ
          IF ( .NOT. State_Met%InTroposphere(1,J,nZ+1-L) ) THEN
             outTmp(J,L) = 0.0e+00_r8
          ENDIF
       ENDDO
       ENDDO
       CALL Outfld( 'SAD_TROP', outTmp(:nY,:), nY, LCHNK )
    ENDIF

    IF ( hist_fld_active('REFF_AERO') ) THEN
       !outTmp(:nY,:) = State_Chm%AeroRadi(1,:nY,nZ:1:-1,:)
       !CALL Outfld( 'REFF_AERO', outTmp(:nY,:), nY, LCHNK )
    ENDIF

    IF ( hist_fld_active('SULF_TROP') ) THEN
       outTmp(:nY,:) = State_Chm%AeroArea(1,:nY,nZ:1:-1,NDUST+1)
       CALL Outfld( 'SULF_TROP', outTmp(:nY,:), nY, LCHNK )
    ENDIF

    ! ===============================================
    ! Diagnose stratospheric quantities
    ! ===============================================

    outTmp(:nY,:) = State_Chm%Species(i_HNO3)%Conc(1,:nY,nZ:1:-1) * MWDry / MW_HNO3
    CALL Outfld( 'HNO3_GAS', outTmp(:nY,:), nY, LCHNK )

    ! TMMF, this requires to have access to the AERFRAC variable in ucx_mod.
    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,2)
    !CALL Outfld( 'HNO3_STS', outTmp(:nY,:), nY, LCHNK )

    outTmp = 0.0e+00_r8
    DO J = 1, nY
    DO L = 1, nZ
       IF ( State_Met%InTroposphere(1,J,nZ+1-L) ) CYCLE
       outTmp(J,L) = State_Chm%Species(i_NIT)%Conc(1,J,nZ+1-L) * MWDry / MW_NIT
    ENDDO
    ENDDO
    CALL Outfld( 'HNO3_NAT', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = outTmp(:nY,:)            + &
    !                AERFRAC(1,:nY,nZ:1:-1,2) + &
                    State_Chm%Species(i_HNO3)%Conc(1,:nY,nZ:1:-1) * MWDry / MW_HNO3
    CALL Outfld( 'HNO3_TOTAL', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = State_Chm%Species(i_H2O)%Conc(1,:nY,nZ:1:-1) * MWDry / MW_H2O
    CALL Outfld( 'H2O_GAS', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = State_Chm%Species(i_HCl)%Conc(1,:nY,nZ:1:-1) * MWDry / MW_HCl
    CALL Outfld( 'HCL_GAS', outTmp(:nY,:), nY, LCHNK )

    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,3)
    !CALL Outfld( 'HCL_STS', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = 0.0e+00_r8
    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,3)
    outTmp(:nY,:) = outTmp(:nY,:)            + &
                    State_Chm%Species(i_HCl)%Conc(1,:nY,nZ:1:-1) * MWDry / MW_HCl
    CALL Outfld( 'HCL_TOTAL', outTmp(:nY,:), nY, LCHNK )

    ! ===============================================
    ! Diagnose dry deposition velocities and fluxes
    ! ===============================================

    IF ( Input_Opt%LDryD ) THEN
       DO N = 1, State_Chm%nDryDep
          ND = NDVZIND(N)
          SpcName = 'DV_'//to_upper(TRIM(depName(N)))
          IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
          CALL OutFld( TRIM(SpcName), State_Chm%DryDepVel(1,:nY,ND), nY, LCHNK )
       ENDDO

       DO N = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          L = State_Chm%Map_Advect(N)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(L)%Info
          SpcName = 'DF_'//to_upper(TRIM(SpcInfo%Name))

          IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
          ! SurfaceFlux is Emissions - Drydep, but Emissions = 0, as it is applied
          ! externally
          CALL OutFld( TRIM(SpcName), -State_Chm%SurfaceFlux(1,:nY,N), nY, LCHNK )

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! ===============================================
    ! Diagnose surface fluxes (emissions - drydep)
    ! ===============================================

    DO N = iFirstCnst, pcnst
       SpcName = TRIM(sflxnam(N))
       IF ( TRIM(SpcName) == '' ) CYCLE
       IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
       CALL OutFld( TRIM(SpcName), cam_in%cflx(:nY,N), nY, LCHNK )
    ENDDO

    ! ===============================================
    ! Diagnose chemical tendencies
    ! ===============================================

    ! Chemical tendencies in kg/kg/s
    DO N = 1, gas_pcnst
       SpcName = TRIM(srcnam(N))
       IF ( TRIM(SpcName) == '' ) CYCLE
       IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
       CALL OutFld( TRIM(SpcName), mmr_tend(:nY,:nZ,N), nY, LCHNK )
    ENDDO

    ! Chemical tendencies in kg/s
    DO N = 1, gas_pcnst
       SpcName = TRIM(dtchem_name(N))
       IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
       outTmp  = 0.0e+0_r8
       outTmp(:nY,:nZ) = mmr_tend(:nY,:nZ,N) * REAL(State_Met%AD(1,:nY,nZ:1:-1),r8)
       CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDDO

    ! ===============================================
    ! Diagnose fields corresponding to State_Met
    ! ===============================================

    SpcName = 'SZA'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp(:nY,1) = ACOS(MIN(MAX(State_Met%SUNCOS(1,:nY),-1._r8),1._r8))/pi*180.e+0_r8
       CALL Outfld( TRIM(SpcName), outTmp(:nY,1) , nY, LCHNK )
    ENDIF

    SpcName = 'U_SRF'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp(:nY,:) = state%u(:nY,:)
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:) , nY, LCHNK )
    ENDIF

    SpcName = 'V_SRF'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp(:nY,:) = state%v(:nY,:)
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:) , nY, LCHNK )
    ENDIF

    SpcName = 'Q_SRF'
    IF ( hist_fld_active(TRIM(SpcName)) ) THEN
       outTmp(:nY,:) = State_Chm%Species(i_H2O)%Conc(1,:nY,nZ:1:-1)
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:) , nY, LCHNK )
    ENDIF

    ! Cleanup
    Current => NULL()
    Item    => NULL()
    Ptr0d_8 => NULL()
    Ptr1d_8 => NULL()
    Ptr2d_8 => NULL()
    Ptr3d_8 => NULL()

  END SUBROUTINE GC_Diagnostics_Calc

  END MODULE GeosChem_Diagnostics_Mod

