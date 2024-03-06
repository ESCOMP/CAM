#define _ASSERT(cond,msg) if(.not.cond) then; print *, "assertion error: ", Iam, __LINE__; call endrun("assertion error - look above - in geoschem_history_mod.F90"); endif
#define _Iam_(name) character(len=255) :: Iam=name
#define __Iam__(name) integer :: STATUS; _Iam_(name)
! Above are compatibility shorthands to avoid excessive divergence from
! MAPL-based code. (hplin, 10/19/22)
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geoschem_history_mod.F90
!
! !DESCRIPTION: Module GeosChem\_History\_Mod interfaces between the CAM history
!  component, the HISTORY.rc configuration file, and the GEOS-Chem State registry.
!  This module is based off GCHP\_HistoryExports\_Mod originally developed by
!  Lizzie Lundgren for GCHP.
!\\
!\\
! !INTERFACE:
!
MODULE GeosChem_History_Mod
!
! !USES:
!
  ! CAM modules
  USE cam_abortutils,     ONLY : endrun

  ! GEOS-Chem modules
  USE DiagList_Mod,       ONLY : DgnItem, DgnList
  USE DiagList_Mod,       ONLY : Init_DiagList, Print_DiagList
  USE ErrCode_Mod,        ONLY : GC_SUCCESS, GC_FAILURE, GC_ERROR
  USE Precision_Mod,      ONLY : fp, f4, f8
  USE TaggedDiagList_Mod, ONLY : TaggedDgnList
  USE TaggedDiagList_Mod, ONLY : Init_TaggedDiagList, Print_TaggedDiagList
  
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HistoryExports_SetServices
  PUBLIC :: HistoryExports_SetDataPointers
  PUBLIC :: CopyGCStates2Exports
  PUBLIC :: Destroy_HistoryConfig
!
! !PRIVATE:
!
  PRIVATE :: Init_HistoryConfig
  PRIVATE :: Init_HistoryExport
  PRIVATE :: Init_HistoryExportsList
  PRIVATE :: Append_HistoryExportsList
  PRIVATE :: Check_HistoryExportsList
  PRIVATE :: Print_HistoryExportsList
 !
! !PUBLIC TYPES
!
  ! History Configuration Object
  TYPE, PUBLIC :: HistoryConfigObj

     CHARACTER(LEN=255)                   :: ROOT ! TODO: needed?
     CHARACTER(LEN=255)                   :: ConfigFileName
     LOGICAL                              :: ConfigFileRead
     TYPE(HistoryExportsListObj), POINTER :: HistoryExportsList
     TYPE(DgnList)                        :: DiagList
     TYPE(TaggedDgnList)                  :: TaggedDiagList

 END TYPE HistoryConfigObj
!
! !PRIVATE TYPES
!
  ! History Exports Linked List
  TYPE :: HistoryExportsListObj

     TYPE(HistoryExportObj), POINTER :: head
     INTEGER                         :: numExports

  END TYPE HistoryExportsListObj

  ! History Export Object
  TYPE :: HistoryExportObj

     CHARACTER(LEN=255)              :: name
     CHARACTER(LEN=255)              :: metadataID
     CHARACTER(LEN=255)              :: registryID
     CHARACTER(LEN=255)              :: long_name
     CHARACTER(LEN=255)              :: units
     INTEGER                         :: vloc
     INTEGER                         :: rank
     INTEGER                         :: type
     LOGICAL                         :: isMet
     LOGICAL                         :: isChem
     LOGICAL                         :: isDiag
     TYPE(HistoryExportObj), POINTER :: next

     ! Pointers to temporaries for CAM Export and GEOS-Chem State
     ! TODO: for now, include all possible data types in the registry.
     REAL(fp), POINTER :: GCStateData0d
     REAL(fp), POINTER :: GCStateData1d(:)
     REAL(fp), POINTER :: GCStateData2d(:,:)
     REAL(fp), POINTER :: GCStateData3d(:,:,:)
     REAL(f4), POINTER :: GCStateData0d_4
     REAL(f4), POINTER :: GCStateData1d_4(:)
     REAL(f4), POINTER :: GCStateData2d_4(:,:)
     REAL(f4), POINTER :: GCStateData3d_4(:,:,:)
     REAL(f8), POINTER :: GCStateData0d_8
     REAL(f8), POINTER :: GCStateData1d_8(:)
     REAL(f8), POINTER :: GCStateData2d_8(:,:)
     REAL(f8), POINTER :: GCStateData3d_8(:,:,:)
     INTEGER,  POINTER :: GCStateData0d_I
     INTEGER,  POINTER :: GCStateData1d_I(:)
     INTEGER,  POINTER :: GCStateData2d_I(:,:)
     INTEGER,  POINTER :: GCStateData3d_I(:,:,:)

  END TYPE HistoryExportObj
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version for GCHP/GEOS
!  19 Oct 2022 - H.P. Lin    - Adapted for CESM
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryConfig
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryConfig ( am_I_Root, HistoryConfig, configFile, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN) :: configFile
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
    INTEGER, INTENT(OUT)            :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_HistoryConfig (geoschem_history_mod.F90)')
    RC = GC_SUCCESS
    ALLOCATE(HistoryConfig)
    HistoryConfig%ROOT               =  ''
    HistoryConfig%ConfigFileName     =  TRIM(configFile)
    HistoryConfig%ConfigFileRead     =  .FALSE.

    CALL Init_DiagList( am_I_Root, configFile, HistoryConfig%DiagList, RC )
    IF ( RC == GC_FAILURE ) THEN
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    ! Optional debugging
    ! CALL Print_DiagList( am_I_Root, HistoryConfig%DiagList, RC )

    CALL Init_TaggedDiagList( am_I_Root, HistoryConfig%DiagList,  &
                              HistoryConfig%TaggedDiagList, RC   )
    IF ( RC == GC_FAILURE ) THEN
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    ! Optional debugging
    ! CALL Print_TaggedDiagList( am_I_Root, HistoryConfig%TaggedDiagList, RC )


    CALL Init_HistoryExportsList( am_I_Root, HistoryConfig, RC )
    IF ( RC == GC_FAILURE ) THEN
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

    ! Optional debugging
    ! CALL Print_HistoryExportsList( am_I_Root, HistoryConfig, RC )

  END SUBROUTINE Init_HistoryConfig
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryExportsList ( am_I_Root, HistoryConfig, RC )
!
! !USES:
!
    ! GEOS-Chem modules
    USE State_Chm_Mod,    ONLY : Get_Metadata_State_Chm
    USE State_Diag_Mod,   ONLY : Get_Metadata_State_Diag
    USE State_Met_Mod,    ONLY : Get_Metadata_State_Met
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N, rank, vloc, type
    CHARACTER(LEN=255)    :: ErrMsg, desc, units, tag
    LOGICAL               :: isMet, isChem, isDiag, found
    TYPE(HistoryExportObj),  POINTER :: NewHistExp
    TYPE(DgnItem),           POINTER :: current

    ! ================================================================
    ! Init_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Init_HistoryExportsList (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Init
    NewHistExp => NULL()

    ! Create HistoryExportsList object
    ALLOCATE(HistoryConfig%HistoryExportsList)
    HistoryConfig%HistoryExportsList%numExports = 0
    HistoryConfig%HistoryExportsList%head => NULL()

    ! Loop over entries in DiagList
    current => HistoryConfig%DiagList%head
    DO WHILE ( ASSOCIATED( current ) )

       ! Skip diagnostics handled by HEMCO, non-standard for GEOS,
       ! or species in the GCHP/GEOS internal state.
       ! See diaglist_mod.F90 for criteria for assigning diagnostic state.
       IF ( INDEX( current%state,  'HEMCO'    ) == 1 .OR. &
            INDEX( current%state,  'GEOS'     ) == 1 .OR. &
            INDEX( current%state,  'INTERNAL' ) == 1 ) THEN
          current => current%next
          CYCLE
       ENDIF

       ! Check history exports list to see if already added (unless wildcard)
       IF ( .NOT. current%isWildcard ) THEN
          CALL Check_HistoryExportsList( am_I_Root, current%name,           &
                                         HistoryConfig%HistoryExportsList,  &
                                         found, RC                         )
          IF ( found ) THEN
             current => current%next
             CYCLE
          ENDIF
       ENDIF

       ! Get metadata using metadataID and state
       ! If isTagged, then append to description
       ! If isWildcard, shouldn't get here
       ! The name of the export is simply name
       Found = .TRUE.
       isMet  = .FALSE.
       isChem = .FALSE.
       isDiag = .FALSE.
       IF ( TRIM(current%state) == 'MET' ) THEN
          isMet = .TRUE.
          CALL Get_Metadata_State_Met( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
          ! TODO: need to add found to outputs of get_metadata_state_met
       ELSEIF ( TRIM(current%state) == 'CHEM' ) THEN
          isCHEM = .TRUE.
          CALL Get_Metadata_State_Chm( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
       ELSEIF ( TRIM(current%state) == 'DIAG' ) THEN
          isDIAG = .TRUE.
          CALL Get_Metadata_State_Diag( am_I_Root, current%metadataID,     &
                                        Found, RC, desc=desc, units=units, &
                                        rank=rank, srcType=type, vloc=vloc )
       ELSE
          RC = GC_FAILURE
          ErrMsg = "Unknown state of item " // TRIM(current%name) // &
                   " in DiagList: " // TRIM(current%state)
          EXIT
       ENDIF

       IF ( .NOT. Found ) THEN
          RC = GC_FAILURE
          ErrMsg = "Metadata not found for " // TRIM(current%name) // &
                   " in state " // TRIM(current%state)
          EXIT
       ENDIF

       ! If wildcard is present
       IF ( current%isWildcard ) THEN
          ! Do nothing. This should never happen at this point since
          ! Init_DiagList will exit with an error if wildcard is
          ! encountered in HISTORY.rc while compiling with ESMF_.

          ! When it comes time to implement, create exports in a loop,
          ! either for all species or for advected species only. Include
          ! a check that the export was not already created. Loop over
          ! AdvNames if wildcard is ADV. Loop over SpecNames for all other
          ! cases, passing not found = OK so that not all are necessarily
          ! output. Later on, after species database is initialized, exports
          ! for only species in the specific wildcard will be associated
          ! with data and thus included in the output file.

          ! If the meantime, skip wildcards if it gets here.
          current => current%next
          CYCLE
       ENDIF

       ! If this item is for a specific tag, append description.
       ! This will need revisiting since there may be tag-dependent
       ! strings to append to long names
       IF ( current%isTagged ) THEN
          desc = TRIM(desc) // " for " // TRIM(current%tag)
       ENDIF

       ! Create a new HistoryExportObj object
       CALL Init_HistoryExport( am_I_Root, NewHistExp,         &
                                name=current%name,             &
                                metadataID=current%metadataID, &
                                registryID=current%registryID, &
                                long_name=desc,                &
                                units=units,                   &
                                vloc=vloc,                     &
                                rank=rank,                     &
                                type=type,                     &
                                isMet=isMet,                   &
                                isChem=isChem,                 &
                                isDiag=isDiag,                 &
                                RC=RC )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export init fail for " // TRIM(current%name)
          EXIT
       ENDIF

       ! Add new HistoryExportObj to linked list
       CALL Append_HistoryExportsList( am_I_Root,     NewHistExp, &
                                       HistoryConfig, RC       )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export append fail for " // TRIM(current%name)
          EXIT
       ENDIF

       ! Set up for next item in DiagList
       current => current%next

    ENDDO
    current => NULL()

    IF ( RC == GC_SUCCESS ) THEN
       HistoryConfig%ConfigFileRead = .TRUE.
    ELSE
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

  END SUBROUTINE Init_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryExport
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryExport ( am_I_Root,  NewHistExp, name,         &
                                  metadataID, registryID, long_name,    &
                                  units,      vloc,       rank,         &
                                  type,       isMet,      isChem,       &
                                  isDiag,     RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,                INTENT(IN) :: am_I_Root
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryExportObj), POINTER    :: NewHistExp
    CHARACTER(LEN=*),       OPTIONAL   :: name
    CHARACTER(LEN=*),       OPTIONAL   :: metadataID
    CHARACTER(LEN=*),       OPTIONAL   :: registryID
    CHARACTER(LEN=*),       OPTIONAL   :: long_name
    CHARACTER(LEN=*),       OPTIONAL   :: units
    INTEGER,                OPTIONAL   :: vloc
    INTEGER,                OPTIONAL   :: rank
    INTEGER,                OPTIONAL   :: type
    LOGICAL,                OPTIONAL   :: isMet
    LOGICAL,                OPTIONAL   :: isChem
    LOGICAL,                OPTIONAL   :: isDiag
    INTEGER,                OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_HistoryExport (geoschem_history_mod.F90)')
    RC = GC_SUCCESS
    ALLOCATE(NewHistExp)

    IF ( PRESENT( name ) ) THEN
       NewHistExp%name = TRIM(name)
    ELSE
       NewHistExp%name = ''
    ENDIF

    IF ( PRESENT( metaDataId ) ) THEN
       NewHistExp%metadataID  = TRIM(metadataID)
    ELSE
       NewHistExp%metadataID  = ''
    ENDIF

    IF ( PRESENT( registryId ) ) THEN
       NewHistExp%registryID = TRIM(registryID)
    ELSE
       NewHistExp%registryId = ''
    ENDIF

    IF ( PRESENT( long_name ) ) THEN
       NewHistExp%long_name = TRIM(long_name)
    ELSE
       NewHistExp%long_name = ''
    ENDIF

    IF ( PRESENT( units ) ) THEN
       NewHistExp%units = TRIM(units)
    ELSE
       NewHistExp%units = ''
    ENDIF

    IF ( PRESENT( vloc ) ) THEN
       NewHistExp%vloc = vloc
    ELSE
       NewHistExp%vloc = -1
    ENDIF

    IF ( PRESENT( rank ) ) THEN
       NewHistExp%rank = rank
    ELSE
       NewHistExp%rank = -1
    ENDIF

    IF ( PRESENT( type ) ) THEN
       NewHistExp%type = type
    ELSE
       NewHistExp%type = -1
    ENDIF

    IF ( PRESENT( isMet ) ) THEN
       NewHistExp%isMet = isMet
    ELSE
       NewHistExp%isMet = .FALSE.
    ENDIF

    IF ( PRESENT( isChem ) ) THEN
       NewHistExp%isChem = isChem
    ELSE
       NewHistExp%isChem = .FALSE.
    ENDIF

    IF ( PRESENT( isDiag ) ) THEN
       NewHistExp%isDiag = isDiag
    ELSE
       NewHistExp%isDiag = .FALSE.
    ENDIF

    NewHistExp%next            => NULL()
    NewHistExp%GCStateData0d   => NULL()
    NewHistExp%GCStateData1d   => NULL()
    NewHistExp%GCStateData2d   => NULL()
    NewHistExp%GCStateData3d   => NULL()
    NewHistExp%GCStateData0d_4 => NULL()
    NewHistExp%GCStateData1d_4 => NULL()
    NewHistExp%GCStateData2d_4 => NULL()
    NewHistExp%GCStateData3d_4 => NULL()
    NewHistExp%GCStateData0d_8 => NULL()
    NewHistExp%GCStateData1d_8 => NULL()
    NewHistExp%GCStateData2d_8 => NULL()
    NewHistExp%GCStateData3d_8 => NULL()
    NewHistExp%GCStateData0d_I => NULL()
    NewHistExp%GCStateData1d_I => NULL()
    NewHistExp%GCStateData2d_I => NULL()
    NewHistExp%GCStateData3d_I => NULL()

  END SUBROUTINE Init_HistoryExport
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Append_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Append_HistoryExportsList ( am_I_Root,     HistoryExport, &
                                         HistoryConfig, RC        )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(HistoryExportObj), POINTER :: HistoryExport
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj),  POINTER :: NewHistExp

    ! ================================================================
    ! Append_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Append_HistoryExportsList (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Add new object to the beginning of the linked list
    HistoryExport%next => HistoryConfig%HistoryExportsList%head
    HistoryConfig%HistoryExportsList%head => HistoryExport

    ! Update # of list items
    HistoryConfig%HistoryExportsList%numExports = &
         HistoryConfig%HistoryExportsList%numExports + 1

  END SUBROUTINE Append_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_HistoryExportsList ( am_I_Root, name,  &
                                        ExportsList, found, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)        :: am_I_Root
    CHARACTER(LEN=*),  INTENT(IN)        :: name
    TYPE(HistoryExportsListObj), POINTER :: ExportsList
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT)               :: found
    INTEGER, INTENT(OUT)               :: RC
!
! !REVISION HISTORY:
!  12 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    __Iam__('Check_HistoryExportsList (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Assume not found
    found = .False.

    current => ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%name == name ) THEN
          found = .TRUE.
          RETURN
       ENDIF
       current => current%next
    ENDDO
    current => NULL()

  END SUBROUTINE Check_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistoryExports_SetServices
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistoryExports_SetServices( am_I_Root, config_file, &
                                         HistoryConfig,  RC )
!
! !USES:
!
    ! CAM modules
    USE cam_history,         ONLY : addfld, add_default, horiz_only

    ! GEOS-Chem modules
    USE Registry_Params_Mod, ONLY : VLocationCenter, VLocationEdge    
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)    :: config_file
!
! !INPUT AND OUTPUT PARAMETERS:
!

!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config object
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version for GCHP/GEOS
!  19 Oct 2022 - H.P. Lin    - Adapted for CESM
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)                   :: ErrMsg
    TYPE(HistoryExportObj),      POINTER :: current

    ! ================================================================
    ! HistoryExports_SetServices begins here
    ! ================================================================

    ! For error handling (defines Iam and STATUS)
    __Iam__('HistoryExports_SetServices (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Create a config object if it does not already exist
    IF ( .NOT. ASSOCIATED(HistoryConfig) ) THEN
       CALL Init_HistoryConfig( am_I_Root, HistoryConfig, config_file, RC )
       IF ( RC == GC_FAILURE ) THEN
          _ASSERT(.FALSE., 'informative message here')
          RETURN
       ENDIF
    ENDIF

    ! Loop over the History Exports list to add one export per item
    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Adding history variables to CAM History State:"
    ENDIF
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) PRINT *, "adding export: ", TRIM(current%name)
       ! Create an export for this item
       IF ( current%rank == 3 ) THEN
          IF ( current%vloc == VLocationCenter ) THEN
            CALL addfld(trim(current%name),                                 &
                        (/'lev'/),                                          &
                        'I',                                                &
                        trim(current%units),                                &
                        trim(current%long_name)                            )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 3D export for " // TRIM(current%name)
             EXIT
          ENDIF
         ELSEIF ( current%vloc == VLocationEdge ) THEN
            CALL addfld(trim(current%name),                                 &
                        (/'ilev'/),                                         &
                        'I',                                                &
                        trim(current%units),                                &
                        trim(current%long_name)                            )
         ELSE
            IF ( am_I_Root ) THEN
               PRINT *, "Unknown vertical location for ", &
                        TRIM(current%name)
            ENDIF
         ENDIF
       ELSEIF ( current%rank == 2 ) THEN
          CALL addfld(trim(current%name),                                 &
                      horiz_only,                                         &
                      'I',                                                &
                      trim(current%units),                                &
                      trim(current%long_name)                            )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 2D export for " // TRIM(current%name)
             EXIT
          ENDIF
       ELSE
          RC = GC_FAILURE
          ErrMsg = "Problem adding export for " // TRIM(current%name) // &
                   ". Rank is only implemented for 2 or 3!"
          EXIT
       ENDIF

       current => current%next
    ENDDO
    current => NULL()

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

  END SUBROUTINE HistoryExports_SetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CopyGCStates2Exports
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CopyGCStates2Exports( am_I_Root, Input_Opt, State_Grid, HistoryConfig, LCHNK, RC )
!
! !USES:
!
    ! CAM modules
    USE cam_history,          ONLY : hist_fld_active, outfld
    USE shr_kind_mod,         ONLY : shr_kind_r8
    
    ! GEOS-Chem modules    
    USE HCO_Interface_GC_Mod, ONLY : HCOI_GC_WriteDiagn
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    INTEGER,        INTENT(IN)    :: LCHNK            ! Chunk number for CESM
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version for GCHP/GEOS
!  19 Oct 2022 - H.P. Lin    - Adapted for CESM
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: LMAX
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! Temporaries for CAM exports.
    ! Note that in CESM, State_Grid%NX is always length 1. (hplin, 11/16/22)
    REAL(shr_kind_r8)               :: outTmp_3D(State_Grid%NY, State_Grid%NZ)
    REAL(shr_kind_r8)               :: outTmp_2D(State_Grid%NY)

    ! ================================================================
    ! CopyGCStates2Exports begins here
    ! ================================================================
    __Iam__('CopyGCStates2Exports (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       ! Skip if not active
       if(.not. hist_fld_active(trim(current%name))) then
          current => current%next
          cycle
       endif

       ! if (am_I_Root) THEN
       !    print *, '  Copying ' // TRIM(current%name)
       ! endif
       IF ( current%rank == 2 ) THEN
          IF ( ASSOCIATED( current%GCStateData2d ) ) THEN
             outTmp_2D(1:State_Grid%NY) = current%GCStateData2d(1,1:State_Grid%NY)
          ELSE IF ( ASSOCIATED( current%GCStateData2d_4 ) ) THEN
             outTmp_2D(1:State_Grid%NY) = current%GCStateData2d_4(1,1:State_Grid%NY)
          ELSE IF ( ASSOCIATED( current%GCStateData2d_8 ) ) THEN
             outTmp_2D(1:State_Grid%NY) = current%GCStateData2d_8(1,1:State_Grid%NY)
          ELSE IF ( ASSOCIATED( current%GCStateData2d_I ) ) THEN
             ! Convert integer to float (integers not allowed in MAPL exports)
             outTmp_2D(1:State_Grid%NY) = FLOAT(current%GCStateData2d_I(1,1:State_Grid%NY))
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 2D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF

          ! Now call outfld to output for this chunk
          call outfld(trim(current%name),                              &
                      outTmp_2D,                                       &   ! Chunk width always 1
                      State_Grid%NY,                                   &
                      LCHNK                                           )
       ELSEIF ( current%rank == 3 ) THEN
          IF ( ASSOCIATED( current%GCStateData3d ) ) THEN
             outTmp_3D(1:State_Grid%NY, :) = current%GCStateData3d(1,1:State_Grid%NY,:)
          ELSE IF ( ASSOCIATED( current%GCStateData3d_4 ) ) THEN
             outTmp_3D(1:State_Grid%NY, :) = current%GCStateData3d_4(1,1:State_Grid%NY,:)
          ELSE IF ( ASSOCIATED( current%GCStateData3d_8 ) ) THEN
             outTmp_3D(1:State_Grid%NY, :) = current%GCStateData3d_8(1,1:State_Grid%NY,:)
          ELSE IF ( ASSOCIATED( current%GCStateData3d_I ) ) THEN
             outTmp_3D(1:State_Grid%NY, :) = FLOAT(current%GCStateData3d_I(1,1:State_Grid%NY,:))
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 3D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF
#if defined( MODEL_CESM )
          ! If using GEOS-5, flip the data vertically to match model
          ! convention
          ! Also do this in CESM. (hplin, 10/31/22)
          LMAX = SIZE(outTmp_3D, 2)
          outTmp_3D(:,1:LMAX) = outTmp_3D(:,LMAX:1:-1)
#endif

          ! Now call outfld to output for this chunk
          call outfld(trim(current%name),                              &
                      outTmp_3D,                                       &   ! Chunk width always 1. TOA is 1
                      State_Grid%NY,                                   &
                      LCHNK                                           )
       ENDIF

       current => current%next
    ENDDO
    current => NULL()

    ! Error handling
    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
  END SUBROUTINE CopyGCStates2Exports
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_HistoryExportsList( am_I_Root, HistoryConfig, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! Print_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Print_HistoryExportsList (geoschem_history_mod.F90)')
    RC = GC_SUCCESS

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    IF ( am_I_Root ) PRINT *, '==========================='
    IF ( am_I_Root ) PRINT *, 'History Exports List:'
    IF ( am_I_Root ) PRINT *, ' '
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) THEN
          PRINT *, "Name:        ",   TRIM(current%name)
          PRINT *, " MetadataID: ",   TRIM(current%metadataID)
          PRINT *, " RegistryID: ",   TRIM(current%registryID)
          PRINT *, " Long name:  ",   TRIM(current%long_name)
          PRINT *, " Units:      ",   TRIM(current%units)
          PRINT *, " Vert loc:   ",   current%vloc
          PRINT *, " Rank:       ",   current%rank
          PRINT *, " Type:       ",   current%type
          PRINT *, " isMet:      ",   current%isMet
          PRINT *, " isChem:     ",   current%isChem
          PRINT *, " isDiag:     ",   current%isDiag
          PRINT *, " "
       ENDIF
       current => current%next
    ENDDO
    IF ( am_I_Root ) PRINT *, '==========================='
    current => NULL()

  END SUBROUTINE Print_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistoryExports_SetDataPointers
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistoryExports_SetDataPointers( am_I_Root,    &
                                             HistoryConfig, State_Chm, &
                                             State_Grid,               &
                                             State_Diag,    State_Met, &
                                             RC                       )
!
! !USES:
!
    ! CAM modules
    USE cam_history,    ONLY : hist_fld_active

    ! GEOS-Chem modules
    USE Registry_Mod,   ONLY : Registry_Lookup
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)     :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config obj
    TYPE(GrdState),   INTENT(INOUT) :: State_Grid     ! Grid State obj
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry State obj
    TYPE(MetState),   INTENT(INOUT) :: State_Met      ! Meteorology State obj
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag     ! Diagnostics State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! HistoryExports_SetDataPointers begins here
    ! ================================================================
    __Iam__('HistoryExports_SetDataPointers')
    RC = GC_SUCCESS

    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Setting history variable pointers to GC and Export States"
    ENDIF

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       ! Skip if not active
       if(.not. hist_fld_active(trim(current%name))) then
          current => current%next
          cycle
       endif

       ! Get pointer to GC state data
       !IF ( am_I_Root ) WRITE(6,*) current%name
       IF ( current%isMET ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Met%Registry,      &
                                RegDict   = State_Met%RegDict,       &
                                State     = State_Met%State,         &
                                Variable  = current%registryID,      &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d_8, &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d_8, &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isChem ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Chm%Registry,      &
                                RegDict   = State_Chm%RegDict,       &
                                State     = State_Chm%State,         &
                                Variable  = current%registryID,      &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d_8, &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d_8, &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isDiag ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Diag%Registry,     &
                                RegDict   = State_Diag%RegDict,      &
                                State     = State_Diag%State,        &
                                Variable  = current%registryID,      &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d_8, &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d_8, &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ENDIF
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "Registry pointer not found for " // TRIM(current%name) // &
                   ". Check that the tag (e.g. species) is valid "         // &
                   "for this diagnostic."
          EXIT
       ENDIF

       !! debugging
       !IF ( Am_I_Root) THEN
       !   WRITE(6,*) TRIM(current%name)
       !ENDIF

       current => current%next
    ENDDO
    current => NULL()

    ! Optional debugging
    !WRITE(6,*) "hplin debug: after HistoryExports_SetDataPointers"
    !CALL Print_HistoryExportsList( am_I_Root, HistoryConfig, RC )

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

  END SUBROUTINE HistoryExports_SetDataPointers
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Destroy_HistoryConfig
!
! !DESCRIPTION: Subroutine Destroy_HistoryConfig deallocates a HistoryConfig
!  object and all of its member objects including the linked list of
!  HistoryExport objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Destroy_HistoryConfig ( am_I_Root, HistoryConfig, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root     ! root CPU?
    TYPE(HistoryConfigObj), POINTER   :: HistoryConfig
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current
    TYPE(HistoryExportObj), POINTER :: next

    ! ================================================================
    ! Destroy_HistoryConfig begins here
    ! ================================================================
    __Iam__('Destroy_HistoryConfig (geoschem_history_mod.F90)')

    current => NULL()
    next => NULL()

    ! Destroy each item in the linked list of HistoryExport objects
    current => HistoryConfig%HistoryExportsList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       _ASSERT( RC == GC_SUCCESS, 'informative message here' )
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO

    ! Deallocate the HistoryExportsList object
    DEALLOCATE( HistoryConfig%HistoryExportsList, STAT=RC )
    _ASSERT( RC == GC_SUCCESS, 'informative message here' )

    ! Deallocate the HistoryConfig object
    DEALLOCATE( HistoryConfig, STAT=RC )
    _ASSERT( RC == GC_SUCCESS, 'informative message here' )

    ! Final cleanup
    current => NULL()
    next    => NULL()

  END SUBROUTINE Destroy_HistoryConfig
!EOC
END MODULE GeosChem_History_Mod
