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
  USE SPMD_UTILS,          ONLY : MasterProc
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

  ! Number of diagnosed photolytic reactions
  INTEGER :: nPhotol
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
  SUBROUTINE CESMGC_Diag_Init( Input_Opt, State_Chm )
!
! !USES:
!
  USE Input_Opt_Mod,       ONLY : OptInput
  USE State_Chm_Mod,       ONLY : ChmState
  USE State_Diag_Mod,      ONLY : get_TagInfo
  USE Species_Mod,         ONLY : Species
  USE CHEM_MODS,           ONLY : nTracers, nSls
  USE CHEM_MODS,           ONLY : tracerNames, tracerLongNames
  USE CHEM_MODS,           ONLY : slsNames, slsLongNames
  USE CHEM_MODS,           ONLY : gas_pcnst
  USE MO_TRACNAME,         ONLY : solsym
  USE CONSTITUENTS,        ONLY : pcnst, cnst_name
  USE CAM_HISTORY,         ONLY : addfld, add_default, horiz_only
  USE GAS_WETDEP_OPTS,     ONLY : gas_wetdep_method
  USE DRYDEP_MOD,          ONLY : depName
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  20 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Integer
    INTEGER                :: I, L, M, N
    INTEGER                :: RC

    ! Logical
    LOGICAL                :: Found

    ! Strings
    CHARACTER(LEN=255)     :: SpcName
    CHARACTER(LEN=255)     :: tagName
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=255)     :: ErrMsg

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! CESMGC_Diag_Init begins here!
    !=================================================================

    ! Initialize pointers
    SpcInfo   => NULL()

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default
    ! Add all species as output fields if desired
    DO I = 1, nTracers
       SpcName = TRIM(tracerNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
          TRIM(tracerLongNames(I))//' concentration')
       IF (TRIM(SpcName) == 'O3') THEN
          CALL Add_Default ( TRIM(SpcName), 1, ' ')
       ENDIF
    ENDDO

    DO I = 1, nSls
       SpcName = TRIM(slsNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', &
          TRIM(slsLongNames(I))//' concentration')
       !CALL Add_Default(TRIM(SpcName), 1, ' ')
    ENDDO

    IF ( Input_Opt%LDryD ) THEN
       DO I = 1, State_Chm%nDryDep
          SpcName = 'DepVel_'//TRIM(depName(I))
          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'm/s', &
            TRIM(SpcName)//' dry deposition velocity')
       ENDDO

       DO I = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          L = State_Chm%Map_Advect(I)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(L)%Info
          SpcName = 'DepFlux_'//TRIM(SpcInfo%Name)

          CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
             TRIM(SpcName)//' dry deposition flux')

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! Surface fluxes (emissions - drydep)
    DO I = 1, pcnst
       SpcName = 'SurfFlux_'//TRIM(cnst_name(I))
       CALL AddFld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
          TRIM(SpcName)//' surface flux')
    ENDDO

    IF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
       DO N = 1, gas_pcnst
          SpcName = 'DTWR_'//TRIM(solsym(N))
          CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/kg/s', &
             'wet removal tendency' )
          SpcName = 'WD_'//TRIM(solsym(N))
          CALL Addfld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
             'vertical integrated wet deposition flux' )
          SpcName = 'WDRATE_'//TRIM(solsym(N))
          CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'kg/s', &
             'wet removal rate' )
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
                               State_Grid, State_Met, cam_in,    LCHNK )
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
  USE MO_TRACNAME,         ONLY : solsym
  USE CAM_HISTORY,         ONLY : outfld
  USE CONSTITUENTS,        ONLY : pcnst, cnst_name
  USE CHEM_MODS,           ONLY : MWRatio
  USE CHEM_MODS,           ONLY : SlsMWRatio
  USE CHEM_MODS,           ONLY : tracerNames
  USE CHEM_MODS,           ONLY : slsNames
  USE CHEM_MODS,           ONLY : nSls
  USE CHEM_MODS,           ONLY : map2GC, map2Idx, map2GC_Sls
  USE DRYDEP_MOD,          ONLY : depName, Ndvzind
  USE CAMSRFEXCH,          ONLY : cam_in_t
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)    :: State_Diag  ! Diag State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(cam_in_t), INTENT(IN)    :: cam_in      ! import state
    INTEGER,        INTENT(IN)    :: LCHNK       ! Chunk number
!
! !REVISION HISTORY:
!  20 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Integers
    INTEGER                :: I, J, L, M, N, ND
    INTEGER                :: RC

    INTEGER                :: nY, nZ

    ! Logicals
    LOGICAL                :: Found

    ! Strings
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: SpcName
    CHARACTER(LEN=255)     :: tagName

    ! Arrays
    REAL(r8)               :: outTmp(State_Grid%nY,State_Grid%nZ)

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! CESMGC_Diag_Calc begins here!
    !=================================================================

    nY = State_Grid%nY
    nZ = State_Grid%nZ

    ! Initialize pointers
    SpcInfo   => NULL()

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at CESMGC_Diag_Calc (in chemistry/pp_geoschem/cesmgc_diag_mod.F90)'

    ! Write diagnostic output
    DO N = 1, pcnst
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
          SpcName = 'DepVel_'//TRIM(depName(N))
          CALL OutFld( TRIM(SpcName), State_Chm%DryDepVel(1,:nY,ND), nY, LCHNK )
       ENDDO

       DO N = 1, State_Chm%nAdvect
          ! Get the species ID from the advected species ID
          L = State_Chm%Map_Advect(N)

          ! Get info about this species from the species database
          SpcInfo => State_Chm%SpcData(L)%Info
          SpcName = 'DepFlux_'//TRIM(SpcInfo%Name)

          ! SurfaceFlux is Emissions - Drydep, but Emissions = 0, as it is applied
          ! externally
          CALL OutFld( TRIM(SpcName), -State_Chm%SurfaceFlux(1,:nY,N), nY, LCHNK )

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    ! Surface fluxes (emissions - drydep)
    DO N = 1, pcnst
       SpcName = 'SurfFlux_'//TRIM(cnst_name(N))
       CALL OutFld( TRIM(SpcName), cam_in%cflx(:nY,N), nY, LCHNK )
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

  END SUBROUTINE CESMGC_Diag_Calc
!EOC
!------------------------------------------------------------------------------
  END MODULE CESMGC_Diag_Mod

