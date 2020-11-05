!------------------------------------------------------------------------------
!            "GEOS-Chem" chemistry emissions interface                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cesmgc_emissions_mod.F90
!
! !DESCRIPTION: Module cesmgc\_emissions\_mod contains routines which retrieve
!  emission fluxes from HEMCO and transfers it back to the CESM-GC interface
!\\
!\\
! !INTERFACE:
!
MODULE CESMGC_Emissions_Mod
!
! !USES:
!
  USE SHR_KIND_MOD,        ONLY : r8 => shr_kind_r8
  USE SPMD_UTILS,          ONLY : MasterProc
  USE CAM_LOGFILE,         ONLY : iulog

  IMPLICIT NONE

  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CESMGC_Emissions_Init
  PUBLIC  :: CESMGC_Emissions_Calc
  PUBLIC  :: CESMGC_Emissions_Final

  INTEGER :: iNO
!
! !REVISION HISTORY:
!  07 Oct 2020 - T. M. Fritz   - Initial version
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
! !IROUTINE: cesmgc_emissions_init
!
! !DESCRIPTION: Subroutine CESMGC\_Emissions\_Init initializes the emissions
!  routine
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CESMGC_Emissions_Init
!
! !USES:
!
    USE CONSTITUENTS,        ONLY : cnst_get_ind
!
! !INPUT PARAMETERS:
!
!
! !REVISION HISTORY:
!  07 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! CESMGC_Emissions_Init begins here!
    !=================================================================

    ! Get constituent index for NO
    CALL cnst_get_ind('NO', iNO, abort=.True.)

  END SUBROUTINE CESMGC_Emissions_Init
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cesmgc_emissions_calc
!
! !DESCRIPTION: Subroutine CESMGC\_Emissions\_Calc retrieves emission fluxes
!  from HEMCO and returns a 3-D array of emission flux to the CESM-GC
!  interface. On top of passing data, this routine handles a number of checks.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CESMGC_Emissions_Calc( state, hco_pbuf2d, State_Met, eflx )
!
! !USES:
!
    USE State_Met_Mod,       ONLY : MetState
    USE CONSTITUENTS,        ONLY : cnst_name, cnst_get_ind, cnst_mw, pcnst
    USE CHEM_MODS,           ONLY : tracerNames, nTracers, map2GCinv
    USE CAM_ABORTUTILS,      ONLY : endrun
    USE PHYSICS_TYPES,       ONLY : physics_state
    USE PHYSICS_BUFFER,      ONLY : pbuf_get_index, pbuf_get_chunk
    USE PHYSICS_BUFFER,      ONLY : physics_buffer_desc, pbuf_get_field
    USE PPGRID,              ONLY : pcols, pver, begchunk
    USE MO_LIGHTNING,        ONLY : prod_NO
    uSE PHYSCONSTANTS,       ONLY : AVO
!
! !INPUT PARAMETERS:
!
    TYPE(physics_state),                INTENT(IN)  :: state           ! Physics state variables
    TYPE(physics_buffer_desc), POINTER, INTENT(IN)  :: hco_pbuf2d(:,:) ! Pointer to 2-D pbuf
    TYPE(MetState),                     INTENT(IN)  :: State_Met       ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
     ! 3-D emissions in kg/m2/s
     REAL(r8),                          INTENT(OUT) :: eflx(pcols,pver,pcnst)
!
! !REVISION HISTORY:
!  07 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Integers
    INTEGER       :: LCHNK, nY, nZ
    INTEGER       :: M, N, J, L
    INTEGER       :: RC                    ! return code
    INTEGER       :: tmpIdx                ! pbuf field id

    ! Logical
    LOGICAL       :: rootChunk
    LOGICAL, SAVE :: FIRST = .True.

    ! Objects
    TYPE(physics_buffer_desc), POINTER :: pbuf_chnk(:) ! slice of pbuf in current chunk

    ! Real
    REAL(r8),                  POINTER :: pbuf_ik(:,:) ! pointer to pbuf data (/pcols,pver/)
    REAL(r8)                           :: SCALFAC      ! Multiplying factor

    ! Strings
    CHARACTER(LEN=255)                 :: fldname_ns   ! field name HCO_*

    !=================================================================
    ! CESMGC_Emissions_Calc begins here!
    !=================================================================

    ! Initialize pointers
    pbuf_chnk => NULL()
    pbuf_ik   => NULL()

    ! LCHNK: which chunk we have on this process
    LCHNK = state%LCHNK
    ! nY: number of atmospheric columns on this chunk
    nY    = state%NCOL
    nZ    = PVER
    rootChunk = ( MasterProc .AND. ( LCHNK.EQ.BEGCHUNK ) )

    ! Initialize emission flux
    eflx(:,:,:) = 0.0e+0_r8

    DO N = 1, nTracers

       fldname_ns = 'HCO_' // TRIM(tracerNames(N))
       tmpIdx = pbuf_get_index(fldname_ns, RC)

       IF ( tmpIdx < 0 ) THEN
          IF ( rootChunk ) Write(iulog,*) "CESMGC_Emissions_Calc: Field not found ", &
             TRIM(fldname_ns)
       ELSE
          ! This is already in chunk, retrieve it
          pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
          CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)

          IF ( .NOT. ASSOCIATED(pbuf_ik) ) THEN ! Sanity check
             CALL ENDRUN("CESMGC_Emissions_Calc: FATAL - tmpIdx > 0 but pbuf_ik not associated")
          ENDIF

          M = map2GCinv(N)

          IF ( M <= 0 ) CYCLE

          eflx(1:nY,:,M) = pbuf_ik(1:nY,:)

          ! Reset pointers
          pbuf_ik   => NULL()
          pbuf_chnk => NULL()

          IF ( MINVAL(eflx(:,:,M)) < 0.0e+00_r8 ) THEN
             Write(iulog,*) "CESMGC_Emissions_Calc: HEMCO emission flux is negative for ", &
                TRIM(cnst_name(M)), " with value ", MINVAL(eflx(:,:,M)), " at ", &
                MINLOC(eflx(:,:,M))
          ENDIF

          IF ( rootChunk .and. ( MAXVAL(eflx(1:nY,:,M)) > 0.0e+0_r8 ) ) THEN
             Write(iulog,'(a,a,a,a)') "CESMGC_Emissions_Calc: HEMCO flux ", &
                TRIM(fldname_ns), " added to ", TRIM(cnst_name(M))
             Write(iulog,'(a,a,E16.4)') "CESMGC_Emissions_Calc: Maximum flux ", &
                TRIM(fldname_ns), MAXVAL(eflx(1:nY,:,M))
          ENDIF
       ENDIF
    ENDDO

    ! Now add lightning emissions computed from lighning_no_prod
    M = iNO

    ! prod_NO is in atom N cm^-3 s^-1 <=> molec cm^-3 s^-1
    ! We need to convert this to kg NO/m2/s
    ! Multiply by AVO * MWNO * BXHEIGHT * 1.0E+06 
    !           = molec/mole * kg NO/mole * m * cm^3/m^3
    ! cnst_mw(M) is in g/mole
    SCALFAC = AVO * cnst_mw(M) * 1.0E-03 * 1.0E+06
    DO J = 1, nY
    DO L = 1, nZ
       eflx(J,L,M) = eflx(J,L,M)                      &
                   + prod_NO(J,L,LCHNK)               &
                     * State_Met%BXHEIGHT(1,J,nZ+1-L) &
                     * SCALFAC
    ENDDO
    ENDDO

    IF ( FIRST ) FIRST = .False.

  END SUBROUTINE CESMGC_Emissions_Calc
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cesmgc_emissions_final
!
! !DESCRIPTION: Subroutine CESMGC\_Emissions\_Final cleans up the module
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CESMGC_Emissions_Final
!
! !REVISION HISTORY:
!  07 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! CESMGC_Emissions_Final begins here!
    !=================================================================
!EOC
!------------------------------------------------------------------------------
  END SUBROUTINE CESMGC_Emissions_Final
!EOC
  END MODULE CESMGC_Emissions_Mod
