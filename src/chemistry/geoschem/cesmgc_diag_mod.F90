!------------------------------------------------------------------------------
!            "GEOS-Chem" chemistry diagnostics interface                      !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cesmgc_diag_mod.F90
!
! !DESCRIPTION: Module cesmgc\_diag\_mod contains routines which aim to
!  diagnose variables from GEOS-Chem
!\\
!\\
! !INTERFACE:
!
MODULE CESMGC_Diag_Mod
!
! !USES:
!
  USE SHR_KIND_MOD,        ONLY : r8 => shr_kind_r8
  USE SHR_CONST_MOD,       ONLY : pi => shr_const_pi
  USE CAM_HISTORY,         ONLY : fieldname_len
  USE CONSTITUENTS,        ONLY : pcnst
  USE CHEM_MODS,           ONLY : gas_pcnst, map2chm
  USE CHEM_MODS,           ONLY : iFirstCnst
  USE MO_TRACNAME,         ONLY : solsym
  USE SPMD_UTILS,          ONLY : MasterProc
  USE PPGRID,              ONLY : begchunk, pver
  USE CAM_LOGFILE,         ONLY : iulog
  USE STRING_UTILS,        ONLY : to_upper
  USE Error_Mod                                 ! For error checking
  USE ErrCode_Mod                               ! Error codes for success or failure

  IMPLICIT NONE

  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CESMGC_Diag_Init
  PUBLIC :: CESMGC_Diag_Calc
  PUBLIC :: wetdep_name, wtrate_name

  INTEGER                      :: nPhotol                ! Number of diagnosed photolytic reactions
  CHARACTER(LEN=fieldname_len) :: srcnam(gas_pcnst)      ! Names of source/sink tendencies
  CHARACTER(LEN=fieldname_len) :: wetdep_name(gas_pcnst) ! Wet deposition tendencies
  CHARACTER(LEN=fieldname_len) :: wtrate_name(gas_pcnst) ! Column tendencies for wet dep
  CHARACTER(LEN=fieldname_len) :: dtchem_name(gas_pcnst) ! Chemical tendencies
  CHARACTER(LEN=16)            :: sflxnam_loc(pcnst)     ! Names of surface fluxes

  ! Chemical families
  INTEGER  :: NOx_species(3)
  INTEGER  :: NOy_species(63)
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
  REAL(r8) :: NOy_MWs(64)
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
             i_R4N2, i_HONIT, i_IONITA, i_NIT, i_NITs, i_NH4
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
  INTEGER :: i_NH3 !NH4 already defined in NOy_species
  ! TOTH
  INTEGER :: i_CH4, i_H2O, i_H2
!
! !REVISION HISTORY:
!  28 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cesmgc_diag_init
!
! !DESCRIPTION: Subroutine CESMGC\_Diag\_Init declares the variables to
!  diagnosethe
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CESMGC_Diag_Init( Input_Opt, State_Chm, State_Met )
!
! !USES:
!
  USE Input_Opt_Mod,       ONLY : OptInput
  USE State_Chm_Mod,       ONLY : ChmState
  USE State_Met_Mod,       ONLY : MetState
  USE State_Diag_Mod,      ONLY : get_TagInfo
  USE Species_Mod,         ONLY : Species
  USE Registry_Mod,        ONLY : MetaRegItem, RegItem
  USE State_Chm_Mod,       ONLY : Ind_
  USE CONSTITUENTS,        ONLY : cnst_name, sflxnam
  USE CONSTITUENTS,        ONLY : cnst_get_ind
  USE CAM_HISTORY,         ONLY : addfld, add_default, horiz_only
  USE GAS_WETDEP_OPTS,     ONLY : gas_wetdep_method
  USE DRYDEP_MOD,          ONLY : depName
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN) :: Input_Opt   ! Input options
    TYPE(ChmState),    INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(MetState),    INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  20 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Integer
    INTEGER                :: M, N, SM
    INTEGER                :: idx
    INTEGER                :: RC

    ! Logical
    LOGICAL                :: Found
    LOGICAL                :: isWD

    ! Strings
    CHARACTER(LEN=255)     :: SpcName
    CHARACTER(LEN=255)     :: tagName
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=2)       :: unit_basename  ! Units 'kg' or '1'

    ! Objects
    TYPE(Species),     POINTER :: SpcInfo
    TYPE(MetaRegItem), POINTER :: Current
    TYPE(RegItem    ), POINTER :: Item

    !=================================================================
    ! CESMGC_Diag_Init begins here!
    !=================================================================

    ! Initialize pointers
    SpcInfo => NULL()
    Current => NULL()
    Item    => NULL()

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    CALL Addfld( 'MASS', (/ 'lev' /), 'A', 'kg', 'Mass of grid box' )
    CALL Addfld( 'AREA', horiz_only,  'A', 'm2', 'Area of grid box' )
    CALL Addfld( 'HEIGHT', (/ 'ilev' /),'A','m', 'Geopotential height above surface at interfaces' )

    ! Note that constituents are already output by default
    ! Add all species as output fields if desired
    DO N = 1, gas_pcnst
       M = map2chm(N)
       IF ( M > 0 ) THEN
          ! It's a GEOS-Chem species
          SpcName = to_upper(TRIM(solsym(N)))
          CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
             TRIM(SpcName)//' volume mixing ratio')
          CALL AddFld( TRIM(SpcName)//'_SRF', horiz_only, 'A', 'mol/mol', &
             TRIM(SpcName)//' in bottom layer')
          IF (TRIM(SpcName) == 'O3') CALL Add_Default( TRIM(SpcName), 2, ' ' )
       ELSE
          ! MAM aerosols
          SpcName = TRIM(solsym(N))
          unit_basename = 'kg'
          IF ( SpcName(1:3) == 'num' ) unit_basename = ' 1'
          CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', unit_basename//'/kg', &
             TRIM(SpcName)//' concentration' )
          CALL AddFld( TRIM(SpcName)//'_SRF', horiz_only, 'A', unit_basename//'/kg', &
             TRIM(SpcName)//' in bottom layer' )
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

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    sflxnam_loc(:) = ''
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
       IF ( map2chm(N) > 0 ) THEN
          !SpcName = to_upper(SpcName)
          IF ( State_Chm%SpcData(map2chm(N))%Info%Is_Gas .eqv. .False. ) THEN
             SpcName = 'GC_AER_'//to_upper(TRIM(SpcName))
          ENDIF
       ENDIF
       CALL cnst_get_ind( SpcName, M, abort=.false. )
       IF ( M > 0 ) THEN
          IF (sflxnam(M)(3:5) == 'num') then  ! name is in the form of "SF****"
             unit_basename = ' 1'
          ELSE
             unit_basename = 'kg'
          ENDIF
          IF ( map2chm(N) > 0 ) THEN
             IF ( State_Chm%SpcData(map2chm(N))%Info%Is_Gas .eqv. .True. ) THEN
                sflxnam_loc(M) = to_upper(sflxnam(M))
             ELSE
                ! Prevent from saving as GC_AER_...
                sflxnam_loc(M) = 'SF'//to_upper(cnst_name(M)(8:))
             ENDIF
          ELSE
             sflxnam_loc(M) = sflxnam(M)
          ENDIF
          SpcName = sflxnam_loc(M)
          CALL Addfld ( TRIM(SpcName), horiz_only, 'A',  unit_basename//'/m2/s', &
             TRIM(solsym(N))//' surface flux')
       ENDIF
    ENDDO

    IF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
       DO N = 1, gas_pcnst
          isWD = .False.
          wetdep_name(N) = 'DTWR_'//to_upper(TRIM(solsym(N)))
          wtrate_name(N) = 'WD_'//to_upper(TRIM(solsym(N)))

          M = map2chm(N)
          IF ( M > 0 ) THEN
             SpcInfo => State_Chm%SpcData(M)%Info
             isWD = SpcInfo%Is_WetDep

             ! Free pointer
             SpcInfo => NULL()
          ENDIF

          IF ( .NOT. isWD ) CYCLE

          SpcName = wetdep_name(N)
          CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg/s', &
             'wet removal tendency' )
          SpcName = wtrate_name(N)
          CALL Addfld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
             'vertical integrated wet deposition flux' )
       ENDDO
    ENDIF

    CALL get_TagInfo( Input_Opt = Input_Opt,  &
                      tagID     = 'PHO',      &
                      State_Chm = State_Chm,  &
                      Found     = Found,      &
                      RC        = RC,         &
                      nTags     = nPhotol    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Abnormal exit from routine "Get_TagInfo", could not '  // &
             ' get nTags!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DO M = 1, nPhotol
       CALL get_TagInfo( Input_Opt = Input_Opt,  &
                         tagID     = 'PHO',      &
                         State_Chm = State_Chm,  &
                         Found     = Found,      &
                         RC        = RC,         &
                         N         = M,          &
                         tagName   = tagName    )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Abnormal exit from routine "Get_TagInfo"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       SpcName = 'Jval_' // TRIM( tagName )
       CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', '1/s', &
          TRIM(tagName) // ' photolysis rate' )
    ENDDO
    ! Add Jval_O3O1D and Jval_O3O3P
    SpcName = 'Jval_O3O1D'
    CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', '1/s', &
       TRIM(tagName) // ' photolysis rate' )
    SpcName = 'Jval_O3O3P'
    CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', '1/s', &
       TRIM(tagName) // ' photolysis rate' )

    ! ==========================================
    ! Now add fields corresponding to State_Met
    ! ==========================================

    ! Copied from Headers/registry_mod.F90
    ! Point to the head node of the Registry
    Current => State_Met%Registry

    ! As long as the current node isn't NULL
    DO WHILE( ASSOCIATED( Current ) )

       ! Get the REGISTRY ITEM belonging to this node of the Registry
       Item => Current%Item

       ! Only print on the root CPU
       IF ( ASSOCIATED( Item ) ) THEN

          !IF (( TRIM(Item%FullName(1:8)) /= 'MET_XLAI' ) .AND. &
          !    ( TRIM(Item%FullName(1:8)) /= 'MET_IUSE' ) .AND. &
          !    ( TRIM(Item%FullName(1:9)) /= 'MET_ILAND' )) THEN
          !   IF ( TRIM(Item%DimNames) == 'xy' ) THEN
          !      CALL Addfld( TRIM( Item%FullName ), horiz_only, 'A', &
          !         TRIM( Item%Units ), TRIM( Item%Description ) )
          !   ELSE
          !      CALL Addfld( TRIM( Item%FullName ), (/ 'lev' /), 'A', &
          !         TRIM( Item%Units ), TRIM( Item%Description ) )
          !   ENDIF
          !ENDIF

       ENDIF

       ! Point to next node of the Registry
       Current => Current%Next

    ENDDO

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
                      i_R4N2, i_HONIT, i_IONITA, i_NIT, i_NITs, i_NH4 /)
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

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    Current => NULL()
    Item    => NULL()

  END SUBROUTINE CESMGC_Diag_Init
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cesmgc_diag_calc
!
! !DESCRIPTION: Subroutine CESMGC\_Diag\_Calc passes the diagnostics variable
!  to the CAM History routines
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CESMGC_Diag_Calc( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, cam_in, state, &
                               mmr_tend,   LCHNK )
!
! !USES:
!
  USE Input_Opt_Mod,       ONLY : OptInput
  USE State_Chm_Mod,       ONLY : ChmState
  USE State_Met_Mod,       ONLY : MetState
  USE State_Diag_Mod,      ONLY : DgnState
  USE State_Diag_Mod,      ONLY : get_TagInfo
  USE State_Grid_Mod,      ONLY : GrdState
  USE Species_Mod,         ONLY : Species
  USE Registry_Mod,        ONLY : MetaRegItem, RegItem
  USE Registry_Mod,        ONLY : Registry_Lookup
  USE Registry_Params_Mod
  USE PRECISION_MOD
  USE CHEM_MODS,           ONLY : adv_mass
  USE CAM_HISTORY,         ONLY : outfld, hist_fld_active
  USE CONSTITUENTS,        ONLY : cnst_name, sflxnam
  USE DRYDEP_MOD,          ONLY : depName, Ndvzind
  USE CAMSRFEXCH,          ONLY : cam_in_t
  USE PHYSICS_TYPES,       ONLY : physics_state
  USE SPMD_UTILS,          ONLY : MasterProc
  USE PHYSCONST,           ONLY : MWDry
  USE UCX_MOD,             ONLY : GET_STRAT_OPT!, AERFRAC
  USE CMN_SIZE_MOD,        ONLY : NDUST
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState),      INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(DgnState),      INTENT(IN)    :: State_Diag  ! Diag State object
    TYPE(GrdState),      INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),      INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(cam_in_t),      INTENT(IN)    :: cam_in      ! import state
    TYPE(physics_state), INTENT(IN)    :: state       ! Physics state variables
    REAL(r8),            INTENT(IN)    :: mmr_tend(state%ncol,pver,gas_pcnst)
                                                      ! Net tendency from chemistry in kg/s
    INTEGER,             INTENT(IN)    :: LCHNK       ! Chunk number
!
! !REVISION HISTORY:
!  20 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Integers
    INTEGER                :: I, J, L, M, N, ND, SM
    INTEGER                :: idx
    INTEGER                :: RC
    INTEGER                :: Source_KindVal    ! KIND value of data
    INTEGER                :: Output_KindVal    ! KIND value for output
    INTEGER                :: Rank              ! Size of data

    INTEGER                :: nY, nZ

    ! Logicals
    LOGICAL                :: Found
    LOGICAL                :: rootChunk
    LOGICAL                :: OnLevelEdges      ! Is the data defined
                                                !  on level edges (T/F)

    ! Strings
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: SpcName
    CHARACTER(LEN=255)     :: tagName

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
    ! CESMGC_Diag_Calc begins here!
    !=================================================================

    nY = State_Grid%nY
    nZ = State_Grid%nZ

    ! Initialize pointers
    SpcInfo => NULL()
    Current => NULL()
    Item    => NULL()

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at CESMGC_Diag_Calc (in chemistry/geoschem/cesmgc_diag_mod.F90)'

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
       IF ( M > 0 ) THEN
          ! It's a GEOS-Chem species
          SpcName = to_upper(TRIM(solsym(N)))
       ELSE
          ! MAM aerosols
          SpcName = TRIM(solsym(N))
       ENDIF
       outTmp = 0.0e+00_r8
       IF ( adv_mass(N) > 0.0e+00_r8 .AND. M /= 0 .AND. hist_fld_active(TRIM(SpcName)) ) THEN
          IF ( M > 0 ) THEN
             outTmp(:nY,:) = REAL(State_Chm%Species(1,:nY,nZ:1:-1,M),r8) * MWDry / adv_mass(N)
          ELSE
             outTmp(:nY,:) = state%q(:nY,:,-M)
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8)
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8)
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
          outTmp(:nY,:) = outTmp(:nY,:) + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8)
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
                        + wgt * REAL(State_Chm%Species(1,:nY,nZ:1:-1,idx),r8) * MWDry / MW
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
          CALL GET_STRAT_OPT(1,J,L,1,RAER,REFF,SADSTRAT,XSASTRAT)
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
          CALL GET_STRAT_OPT(1,J,L,2,RAER,REFF,SADSTRAT,XSASTRAT)
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

    outTmp(:nY,:) = State_Chm%Species(1,:nY,nZ:1:-1,i_HNO3) * MWDry / MW_HNO3
    CALL Outfld( 'HNO3_GAS', outTmp(:nY,:), nY, LCHNK )

    ! TMMF, this requires to have access to the AERFRAC variable in ucx_mod.
    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,2)
    !CALL Outfld( 'HNO3_STS', outTmp(:nY,:), nY, LCHNK )

    outTmp = 0.0e+00_r8
    DO J = 1, nY
    DO L = 1, nZ
       IF ( State_Met%InTroposphere(1,J,nZ+1-L) ) CYCLE
       outTmp(J,L) = State_Chm%Species(1,J,nZ+1-L,i_NIT) * MWDry / MW_NIT
    ENDDO
    ENDDO
    CALL Outfld( 'HNO3_NAT', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = outTmp(:nY,:)            + &
    !                AERFRAC(1,:nY,nZ:1:-1,2) + &
                    State_Chm%Species(1,:nY,nZ:1:-1,i_HNO3) * MWDry / MW_HNO3
    CALL Outfld( 'HNO3_TOTAL', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = State_Chm%Species(1,:nY,nZ:1:-1,i_H2O) * MWDry / MW_H2O
    CALL Outfld( 'H2O_GAS', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = State_Chm%Species(1,:nY,nZ:1:-1,i_HCl) * MWDry / MW_HCl
    CALL Outfld( 'HCL_GAS', outTmp(:nY,:), nY, LCHNK )

    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,3)
    !CALL Outfld( 'HCL_STS', outTmp(:nY,:), nY, LCHNK )

    outTmp(:nY,:) = 0.0e+00_r8
    !outTmp(:nY,:) = AERFRAC(1,:nY,nZ:1:-1,3)
    outTmp(:nY,:) = outTmp(:nY,:)            + &
                    State_Chm%Species(1,:nY,nZ:1:-1,i_HCl) * MWDry / MW_HCl
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
       SpcName = TRIM(sflxnam_loc(N))
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
       CALL OutFld( TRIM(SpcName), mmr_tend(:nY,:,N), nY, LCHNK )
    ENDDO

    ! Chemical tendencies in kg/s
    DO N = 1, gas_pcnst
       SpcName = TRIM(dtchem_name(N))
       IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
       outTmp  = 0.0e+0_r8
       DO J = 1, nY
       DO L = 1, nZ
          outTmp(J,L) = mmr_tend(J,L,N) * REAL(State_Met%AD(1,J,nZ+1-L),r8)
       ENDDO
       ENDDO
       CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDDO

    ! ===============================================
    ! Diagnose photolysis rates
    ! ===============================================

    IF ( ASSOCIATED(State_Diag%Jval) ) THEN
       DO M = 1, nPhotol
          CALL get_TagInfo( Input_Opt = Input_Opt,  &
                            tagID     = 'PHO',      &
                            State_Chm = State_Chm,  &
                            Found     = Found,      &
                            RC        = RC,         &
                            N         = M,          &
                            tagName   = tagName    )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Abnormal exit from routine "Get_TagInfo"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

          SpcName = 'Jval_' // TRIM( tagName )
          IF ( .NOT. hist_fld_active(TRIM(SpcName)) ) CYCLE
          DO J = 1, nY
          DO L = 1, nZ
             outTmp(J,nZ+1-L) = REAL(State_Diag%Jval(1,J,L,M),r8)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
       ENDDO
    ENDIF
    IF ( ASSOCIATED(State_Diag%JvalO3O1D) ) THEN
       SpcName = 'Jval_O3O1D'
       IF ( hist_fld_active(TRIM(SpcName)) ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             outTmp(J,nZ+1-L) = REAL(State_Diag%JvalO3O1D(1,J,L),r8)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
       ENDIF
    ENDIF
    IF ( ASSOCIATED(State_Diag%JvalO3O3P) ) THEN
       SpcName = 'Jval_O3O3P'
       IF ( hist_fld_active(TRIM(SpcName)) ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             outTmp(J,nZ+1-L) = REAL(State_Diag%JvalO3O3P(1,J,L),r8)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
       ENDIF
    ENDIF

    ! ===============================================
    ! Diagnose fields corresponding to State_Met
    ! ===============================================

    ! Copied from Headers/registry_mod.F90
    ! Point to the head node of the Registry
    Current => State_Met%Registry

    Source_KindVal = KINDVAL_F8
    Output_KindVal = KINDVAL_F8

    ! As long as the current node isn't NULL
    DO WHILE( ASSOCIATED( Current ) )

       ! Get the REGISTRY ITEM belonging to this node of the Registry
       Item => Current%Item

       ! Only print on the root CPU
       IF ( ASSOCIATED( Item ) ) THEN

          SpcName = TRIM(Item%FullName)
          IF (( TRIM(Item%FullName(1:8)) /= 'MET_XLAI' ) .AND. &
              ( TRIM(Item%FullName(1:8)) /= 'MET_IUSE' ) .AND. &
              ( TRIM(Item%FullName(1:9)) /= 'MET_ILAND' )) THEN
             CALL Registry_Lookup( am_I_Root      = Input_Opt%amIRoot,   &
                                   Registry       = State_Met%Registry,  &
                                   RegDict        = State_Met%RegDict,   &
                                   State          = State_Met%State,     &
                                   Variable       = Item%FullName,       &
                                   Source_KindVal = Source_KindVal,      &
                                   Output_KindVal = Output_KindVal,      &
                                   Rank           = Rank,                &
                                   OnLevelEdges   = OnLevelEdges,        &
                                   Ptr0d_8        = Ptr0d_8,             &
                                   Ptr1d_8        = Ptr1d_8,             &
                                   Ptr2d_8        = Ptr2d_8,             &
                                   Ptr3d_8        = Ptr3d_8,             &
                                   RC             = RC                  )

             !IF ( hist_fld_active(TRIM(SpcName)) ) THEN
             !   IF ( Source_KindVal /= KINDVAL_I4 ) THEN
             !      IF ( Rank == 2 ) THEN
             !         outTmp(:,nZ) = REAL(Ptr2d_8(1,:),r8)
             !         CALL Outfld( TRIM( Item%FullName ), outTmp(:,nZ), nY, LCHNK )
             !      ELSEIF ( Rank == 3 ) THEN
             !         ! For now, treat variables defined on level edges by ignoring top
             !         ! most layer
             !         DO J = 1, nY
             !         DO L = 1, nZ
             !            outTmp(J,nZ+1-L) = REAL(Ptr3d_8(1,J,L),r8)
             !         ENDDO
             !         ENDDO
             !         CALL Outfld( TRIM( Item%FullName ), outTmp, nY, LCHNK )
             !      ELSE
             !         IF ( rootChunk ) Write(iulog,*) " Item ", TRIM(Item%FullName), &
             !            " is of rank ", Rank, " and will not be diagnosed!"
             !      ENDIF
             !   ENDIF
             !ENDIF
          ENDIF

       ENDIF

       ! Point to next node of the Registry
       Current => Current%Next

    ENDDO

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
       outTmp(:nY,:) = State_Chm%Species(1,:nY,nZ:1:-1,i_H2O)
       CALL Outfld( TRIM(SpcName), outTmp(:nY,:) , nY, LCHNK )
    ENDIF

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    Current => NULL()
    Item    => NULL()
    Ptr0d_8 => NULL()
    Ptr1d_8 => NULL()
    Ptr2d_8 => NULL()
    Ptr3d_8 => NULL()

  END SUBROUTINE CESMGC_Diag_Calc
!EOC
!------------------------------------------------------------------------------
  END MODULE CESMGC_Diag_Mod

