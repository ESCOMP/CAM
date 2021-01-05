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
  USE CAM_HISTORY,         ONLY : fieldname_len
  USE CHEM_MODS,           ONLY : gas_pcnst
  USE SPMD_UTILS,          ONLY : MasterProc
  USE PPGRID,              ONLY : begchunk, pver
  USE CAM_LOGFILE,         ONLY : iulog
  USE Error_Mod                                 ! For error checking
  USE ErrCode_Mod                               ! Error codes for success or failure

  IMPLICIT NONE

  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CESMGC_Diag_Init
  PUBLIC :: CESMGC_Diag_Calc

  INTEGER                      :: nPhotol           ! Number of diagnosed photolytic reactions
  CHARACTER(LEN=FIELDNAME_LEN) :: srcnam(gas_pcnst) ! Names of source/sink tendencies
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
  USE CHEM_MODS,           ONLY : nTracers, nSls
  USE CHEM_MODS,           ONLY : tracerNames, tracerLongNames
  USE CHEM_MODS,           ONLY : slsNames, slsLongNames
  USE CHEM_MODS,           ONLY : iFirstCnst
  USE CHEM_MODS,           ONLY : map2chm
  USE MO_TRACNAME,         ONLY : solsym
  USE CONSTITUENTS,        ONLY : pcnst, cnst_name, sflxnam
  USE CONSTITUENTS,        ONLY : cnst_get_ind
  USE CAM_HISTORY,         ONLY : addfld, add_default, horiz_only
  USE GAS_WETDEP_OPTS,     ONLY : gas_wetdep_method
  USE DRYDEP_MOD,          ONLY : depName
#if defined( MODAL_AERO_4MODE )
  USE MODAL_AERO_DATA,     ONLY : ntot_amode, nspec_amode
  USE MODAL_AERO_DATA,     ONLY : xname_massptr
#endif
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
    INTEGER                :: I, L, M, N, SM
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
    DO I = 1, nTracers
       SpcName = TRIM(tracerNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
          TRIM(tracerLongNames(I))//' volume mixing ratio')
       IF (TRIM(SpcName) == 'O3') THEN
          CALL Add_Default ( TRIM(SpcName), 2, ' ')
       ENDIF
    ENDDO

    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          SpcName = TRIM(xname_massptr(SM,M))
          CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg', &
             TRIM(SpcName//' mass mixing ratio') )
       ENDDO
    ENDDO

    DO I = 1, nSls
       SpcName = TRIM(slsNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
          TRIM(slsLongNames(I))//' volume mixing ratio')
       !CALL Add_Default(TRIM(SpcName), 1, ' ')
    ENDDO

    IF ( Input_Opt%LDryD ) THEN
       DO I = 1, State_Chm%nDryDep
          SpcName = 'DV_'//TRIM(depName(I))
          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'm/s', &
            TRIM(SpcName)//' dry deposition velocity')
       ENDDO

       DO I = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          L = State_Chm%Map_Advect(I)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(L)%Info
          SpcName = 'DF_'//TRIM(SpcInfo%Name)

          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
             TRIM(SpcName)//' dry deposition flux')

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! Chemical tendencies and surface fluxes
    DO N = 1, gas_pcnst
       srcnam(N) = 'CT_'//TRIM(solsym(N)) ! chem tendency (source/sink)
       SpcName = srcnam(N)
       CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg/s', TRIM(SpcName)//' source/sink' )

       CALL cnst_get_ind(solsym(N), M, abort=.false. ) 
       IF ( M > 0 ) THEN
          IF (sflxnam(M)(3:5) == 'num') then  ! name is in the form of "SF****"
             unit_basename = ' 1'
          ELSE
             unit_basename = 'kg'  
          ENDIF
          CALL Addfld (sflxnam(M), horiz_only, 'A',  unit_basename//'/m2/s', &
             TRIM(solsym(N))//' surface flux')
       ENDIF
    ENDDO

    IF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
       DO N = 1, gas_pcnst
          M = map2chm(N)
          isWD = .False.
          IF ( M > 0 ) THEN
             SpcInfo => State_Chm%SpcData(M)%Info
             isWD = SpcInfo%Is_WetDep

             ! Free pointer
             SpcInfo => NULL()
          ENDIF

          IF ( .NOT. isWD ) CYCLE

          SpcName = 'DTWR_'//TRIM(solsym(N))
          CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg/s', &
             'wet removal tendency' )
          SpcName = 'WD_'//TRIM(solsym(N))
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
       SpcName = TRIM(solsym(N))
       CALL Addfld( 'D'//TRIM(SpcName)//'CHM', (/ 'lev' /), 'A', 'kg/s', &
          'net tendency from chemistry' )
    ENDDO

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
  USE MO_TRACNAME,         ONLY : solsym
  USE CAM_HISTORY,         ONLY : outfld
  USE CONSTITUENTS,        ONLY : pcnst, cnst_name, sflxnam
  USE CHEM_MODS,           ONLY : MWRatio
  USE CHEM_MODS,           ONLY : SlsMWRatio
  USE CHEM_MODS,           ONLY : tracerNames
  USE CHEM_MODS,           ONLY : slsNames
  USE CHEM_MODS,           ONLY : nSls
  USE CHEM_MODS,           ONLY : map2GC, map2Idx, map2GC_Sls
  USE CHEM_MODS,           ONLY : iFirstCnst
  USE DRYDEP_MOD,          ONLY : depName, Ndvzind
  USE CAMSRFEXCH,          ONLY : cam_in_t
  USE PHYSICS_TYPES,       ONLY : physics_state
  USE SPMD_UTILS,          ONLY : MasterProc
#if defined( MODAL_AERO_4MODE )
  USE MODAL_AERO_DATA,     ONLY : lmassptr_amode
  USE MODAL_AERO_DATA,     ONLY : ntot_amode, nspec_amode
  USE MODAL_AERO_DATA,     ONLY : xname_massptr
#endif
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

    ! Arrays
    REAL(r8)               :: outTmp(State_Grid%nY,State_Grid%nZ)
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

    CALL OutFld( 'AREA', State_Grid%Area_M2(1,1:nY), nY, LCHNK)
    CALL OutFld( 'MASS', State_Met%AD(1,1:nY,nZ:1:-1), nY, LCHNK)
    CALL Outfld( 'HEIGHT', state%zi(:nY,:), nY, LCHNK )

    ! Write diagnostic output
    DO N = iFirstCnst, pcnst
       M = map2GC(N)
       I = map2Idx(N)
       IF ( M > 0 ) THEN
          SpcName = tracerNames(I)
          outTmp  = 0.0e+0_r8
          DO J = 1, nY
          DO L = 1, nZ
             outTmp(J,nZ+1-L) = REAL(State_Chm%Species(1,J,L,M),r8) * MWRatio(I)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
       ENDIF
    ENDDO

    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          SpcName = TRIM(xname_massptr(SM,M))
          N = lmassptr_amode(SM,M)
          CALL OutFld( TRIM(SpcName), state%q(:nY,:,N), nY, LCHNK )
       ENDDO
    ENDDO

    DO N = 1, nSls
       SpcName = slsNames(N)
       outTmp  = 0.0e+0_r8
       M = map2GC_Sls(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             outTmp(J,nZ+1-L) = REAL(State_Chm%Species(1,J,L,M),r8) * SLSMWratio(N)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
       ENDIF
    ENDDO

    ! Dry deposition velocity and surface flux
    IF ( Input_Opt%LDryD ) THEN
       DO N = 1, State_Chm%nDryDep
          ND = NDVZIND(N)
          SpcName = 'DV_'//TRIM(depName(N))
          CALL OutFld( TRIM(SpcName), State_Chm%DryDepVel(1,:nY,ND), nY, LCHNK )
       ENDDO

       DO N = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          L = State_Chm%Map_Advect(N)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(L)%Info
          SpcName = 'DF_'//TRIM(SpcInfo%Name)

          ! SurfaceFlux is Emissions - Drydep, but Emissions = 0, as it is applied
          ! externally
          CALL OutFld( TRIM(SpcName), -State_Chm%SurfaceFlux(1,:nY,N), nY, LCHNK )

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! Surface fluxes (emissions - drydep)
    DO N = iFirstCnst, pcnst
       SpcName = sflxnam(N)
       CALL OutFld( TRIM(SpcName), cam_in%cflx(:nY,N), nY, LCHNK )
    ENDDO

    ! Chemical tendencies in kg/kg/s
    DO N = 1, gas_pcnst
       SpcName = TRIM(srcnam(N))
       CALL OutFld( TRIM(SpcName), mmr_tend(:nY,:,N), nY, LCHNK )
    ENDDO

    ! Chemical tendencies in kg/s
    DO N = 1, gas_pcnst
       SpcName = 'D'//TRIM(solsym(N))//'CHM'
       outTmp  = 0.0e+0_r8
       DO J = 1, nY
       DO L = 1, nZ
          outTmp(J,L) = mmr_tend(J,L,N) * REAL(State_Met%AD(1,J,nZ+1-L),r8)
       ENDDO
       ENDDO
       CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDDO

    ! Photolysis rates
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
       DO J = 1, nY
       DO L = 1, nZ
          outTmp(J,nZ+1-L) = REAL(State_Diag%JvalO3O1D(1,J,L),r8)
       ENDDO
       ENDDO
       CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF
    IF ( ASSOCIATED(State_Diag%JvalO3O3P) ) THEN
       SpcName = 'Jval_O3O3P'
       DO J = 1, nY
       DO L = 1, nZ
          outTmp(J,nZ+1-L) = REAL(State_Diag%JvalO3O3P(1,J,L),r8)
       ENDDO
       ENDDO
       CALL OutFld( TRIM(SpcName), outTmp(:nY,:), nY, LCHNK )
    ENDIF

    ! ===============================================
    ! Now diagnose fields corresponding to State_Met
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

             !IF ( Source_KindVal /= KINDVAL_I4 ) THEN
             !   IF ( Rank == 2 ) THEN
             !      outTmp(:,nZ) = REAL(Ptr2d_8(1,:),r8)
             !      CALL Outfld( TRIM( Item%FullName ), outTmp(:,nZ), nY, LCHNK )
             !   ELSEIF ( Rank == 3 ) THEN
             !      ! For now, treat variables defined on level edges by ignoring top
             !      ! most layer
             !      DO J = 1, nY
             !      DO L = 1, nZ
             !         outTmp(J,nZ+1-L) = REAL(Ptr3d_8(1,J,L),r8)
             !      ENDDO
             !      ENDDO
             !      CALL Outfld( TRIM( Item%FullName ), outTmp, nY, LCHNK )
             !   ELSE
             !      IF ( rootChunk ) Write(iulog,*) " Item ", TRIM(Item%FullName), &
             !         " is of rank ", Rank, " and will not be diagnosed!"
             !   ENDIF
             !ENDIF
          ENDIF

       ENDIF

       ! Point to next node of the Registry
       Current => Current%Next

    ENDDO

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

