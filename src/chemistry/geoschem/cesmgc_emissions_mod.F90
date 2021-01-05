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
  USE CAM_ABORTUTILS,      ONLY : endrun
  USE SHR_MEGAN_MOD,       ONLY : shr_megan_mechcomps, shr_megan_mechcomps_n 
  USE CAM_LOGFILE,         ONLY : iulog

  IMPLICIT NONE

  PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CESMGC_Emissions_Init
  PUBLIC  :: CESMGC_Emissions_Calc
  PUBLIC  :: CESMGC_Emissions_Final

  ! Constituent number for NO
  INTEGER :: iNO

  ! Aerosol constituent number
  INTEGER :: iBC1
  INTEGER :: iBC4
  INTEGER :: iH2SO4
  INTEGER :: iSOA11
  INTEGER :: iSOA12
  INTEGER :: iSOA21
  INTEGER :: iSOA22
  INTEGER :: iSOA31
  INTEGER :: iSOA32
  INTEGER :: iSOA41
  INTEGER :: iSOA42
  INTEGER :: iSOA51
  INTEGER :: iSOA52
  INTEGER :: iPOM1
  INTEGER :: iPOM4

  INTEGER :: iBCPI
  INTEGER :: iBCPO
  INTEGER :: iOCPI
  INTEGER :: iOCPO
  INTEGER :: iSO4
  INTEGER :: iSOAS

  ! MEGAN Emissions
  INTEGER,  ALLOCATABLE :: megan_indices_map(:) 
  REAL(r8), ALLOCATABLE :: megan_wght_factors(:)

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
  SUBROUTINE CESMGC_Emissions_Init( lght_no_prd_factor )
!
! !USES:
!
    USE PHYSICS_TYPES,       ONLY : physics_state
    USE CONSTITUENTS,        ONLY : cnst_get_ind
    USE MO_CHEM_UTLS,        ONLY : get_spc_ndx
    USE CAM_HISTORY,         ONLY : addfld, add_default, horiz_only
    USE MO_LIGHTNING,        ONLY : lightning_inti
    USE FIRE_EMISSIONS,      ONLY : fire_emissions_init
    USE CHEM_MODS,           ONLY : adv_mass
    USE INFNAN,              ONLY : NaN, assignment(=)
!
! !INPUT PARAMETERS:
!
    REAL(r8),                INTENT(IN   ) :: lght_no_prd_factor ! Lightning scaling factor
!
! !REVISION HISTORY:
!  07 Oct 2020 - T. M. Fritz   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Integers
    INTEGER                :: IERR
    INTEGER                :: N, II

    ! Strings
    CHARACTER(LEN=255)     :: SpcName
    CHARACTER(LEN=255)     :: Description

    ! Real
    REAL(r8)               :: MW

    !=================================================================
    ! CESMGC_Emissions_Init begins here!
    !=================================================================

    ! Get constituent index for NO
    CALL cnst_get_ind('NO', iNO, abort=.True.)

#if defined( MODAL_AERO_4MODE )
    ! Get constituent index for aerosols
    CALL cnst_get_ind('bc_a1',       iBC1,   abort=.True.)
    CALL cnst_get_ind('bc_a4',       iBC4,   abort=.True.)
    CALL cnst_get_ind('soa1_a1',     iSOA11, abort=.True.)
    CALL cnst_get_ind('soa1_a2',     iSOA12, abort=.True.)
    CALL cnst_get_ind('soa2_a1',     iSOA21, abort=.True.)
    CALL cnst_get_ind('soa2_a2',     iSOA22, abort=.True.)
    CALL cnst_get_ind('soa3_a1',     iSOA31, abort=.True.)
    CALL cnst_get_ind('soa3_a2',     iSOA32, abort=.True.)
    CALL cnst_get_ind('soa4_a1',     iSOA41, abort=.True.)
    CALL cnst_get_ind('soa4_a2',     iSOA42, abort=.True.)
    CALL cnst_get_ind('soa5_a1',     iSOA51, abort=.True.)
    CALL cnst_get_ind('soa5_a2',     iSOA52, abort=.True.)
    CALL cnst_get_ind('H2SO4',       iH2SO4, abort=.True.)
    CALL cnst_get_ind('pom_a1',      iPOM1,  abort=.True.)
    CALL cnst_get_ind('pom_a4',      iPOM4,  abort=.True.)

    CALL cnst_get_ind('GC_AER_BCPI', iBCPI,  abort=.True.)
    CALL cnst_get_ind('GC_AER_BCPO', iBCPO,  abort=.True.)
    CALL cnst_get_ind('GC_AER_SOAS', iSOAS,  abort=.True.)
    CALL cnst_get_ind('GC_AER_SO4',  iSO4,   abort=.True.)
    CALL cnst_get_ind('GC_AER_OCPI', iOCPI,  abort=.True.)
    CALL cnst_get_ind('GC_AER_OCPO', iOCPO,  abort=.True.)
#endif

    !-----------------------------------------------------------------------
    !	... initialize the lightning module
    !-----------------------------------------------------------------------
    CALL lightning_inti(lght_no_prd_factor)

    !-----------------------------------------------------------------------
    ! ... MEGAN emissions
    !-----------------------------------------------------------------------
    IF ( shr_megan_mechcomps_n > 0 ) THEN

       ALLOCATE( megan_indices_map(shr_megan_mechcomps_n), STAT=IERR )
       IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating megan_indices_map')
       ALLOCATE( megan_wght_factors(shr_megan_mechcomps_n), STAT=IERR )
       IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating megan_wght_factors')
       megan_wght_factors(:) = NaN

       DO N = 1, shr_megan_mechcomps_n
          SpcName = TRIM(shr_megan_mechcomps(N)%name)
          ! Special handlings for GEOS-Chem species
          IF ( TRIM(SpcName) == 'MTERP' ) THEN
             SpcName = 'MTPA'
          ELSEIF ( TRIM(SpcName) == 'BCARY' ) THEN
             SpcName = 'None'
             MW      = 204.342600_r8 ! Taken from pp_trop_strat_mam4_vbs
          ELSEIF ( TRIM(SpcName) == 'CH3OH' ) THEN
             SpcName = 'MOH'
          ELSEIF ( TRIM(SpcName) == 'C2H5OH' ) THEN
             SpcName = 'EOH'
          ELSEIF ( TRIM(SpcName) == 'CH3CHO' ) THEN
             SpcName = 'ALD2'
          ELSEIF ( TRIM(SpcName) == 'CH3COOH' ) THEN
             SpcName = 'ACTA'
          ELSEIF ( TRIM(SpcName) == 'CH3COCH3' ) THEN
             SpcName = 'ACET'
          ELSEIF ( TRIM(SpcName) == 'HCN' ) THEN
             SpcName = 'None'
             MW      = 27.025140_r8 ! Taken from pp_trop_strat_mam4_vbs
          ELSEIF ( TRIM(SpcName) == 'C2H4' ) THEN
             SpcName = 'None'
             MW      = 28.051600_r8 ! Taken from pp_trop_strat_mam4_vbs
          ELSEIF ( TRIM(SpcName) == 'C3H6' ) THEN
             SpcName = 'PRPE'
          ELSEIF ( TRIM(SpcName) == 'BIGALK' ) THEN
             ! BIGALK = Pentane + Hexane + Heptane + Tricyclene
             SpcName = 'ALK4'
          ELSEIF ( TRIM(SpcName) == 'BIGENE' ) THEN
             ! BIGENE = butene (C4H8)
             SpcName = 'PRPE' ! Lumped >= C3 alkenes
          ELSEIF ( TRIM(SpcName) == 'TOLUENE' ) THEN
             SpcName = 'TOLU'
          ENDIF

          CALL cnst_get_ind (SpcName, megan_indices_map(N), abort=.False.)
          II = get_spc_ndx(SpcName)
          IF ( II > 0 ) THEN
             SpcName = TRIM(shr_megan_mechcomps(N)%name)
             megan_wght_factors(N) = adv_mass(II)*1.e-3_r8 ! kg/moles (to convert moles/m2/sec to kg/m2/sec)
             Description = TRIM(SpcName)//' MEGAN emissions flux (released as '//TRIM(SpcName)//' in GEOS-Chem)'
          ELSEIF ( TRIM(SpcName) == 'None' ) THEN
             SpcName = TRIM(shr_megan_mechcomps(N)%name)
             megan_wght_factors(N) = MW*1.e-3_r8 ! kg/moles
             IF ( MasterProc ) Write(iulog,*) " MEGAN ", TRIM(SpcName), &
                " emissions will be ignored as no species match in GEOS-Chem."
             Description = TRIM(SpcName)//' MEGAN emissions flux (not released in GEOS-Chem)'
          ELSE
             SpcName = TRIM(shr_megan_mechcomps(N)%name)
             CALL ENDRUN( 'chem_init: MEGAN compound not in chemistry mechanism : '//TRIM(SpcName))
          ENDIF

          ! MEGAN  history fields
          CALL Addfld( 'MEG_'//TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
               Description )

          !if (history_chemistry) then
          CALL Add_default('MEG_'//TRIM(SpcName), 1, ' ')
          !endif
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------
    ! ... Fire emissions
    !-----------------------------------------------------------------------
    CALL fire_emissions_init()

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
  SUBROUTINE CESMGC_Emissions_Calc( state, hco_pbuf2d, State_Met, cam_in, eflx, iStep )
!
! !USES:
!
    USE State_Met_Mod,       ONLY : MetState
    USE CAMSRFEXCH,          ONLY : cam_in_t
    USE CONSTITUENTS,        ONLY : cnst_name, cnst_get_ind, cnst_mw, pcnst
    USE CHEM_MODS,           ONLY : tracerNames, nTracers, map2GCinv
    USE PHYSICS_TYPES,       ONLY : physics_state
    USE PHYSICS_BUFFER,      ONLY : pbuf_get_index, pbuf_get_chunk
    USE PHYSICS_BUFFER,      ONLY : physics_buffer_desc, pbuf_get_field
    USE PPGRID,              ONLY : pcols, pver, begchunk
    USE CAM_HISTORY,         ONLY : outfld
    USE STRING_UTILS,        ONLY : to_upper

    ! Data from CLM
    USE CAM_CPL_INDICES,     ONLY : index_x2a_Fall_flxvoc

    ! Lightning emissions
    USE MO_LIGHTNING,        ONLY : prod_NO

    ! Fire emissions
    USE FIRE_EMISSIONS,      ONLY : fire_emissions_srf
    USE FIRE_EMISSIONS,      ONLY : fire_emissions_vrt

    ! Aerosol emissions
    USE AERO_MODEL,          ONLY : aero_model_emissions

    ! GEOS-Chem version of physical constants
    USE PHYSCONSTANTS,       ONLY : AVO
    ! CAM version of physical constants
    USE PHYSCONST,           ONLY : rga
!
! !INPUT PARAMETERS:
!
    TYPE(physics_state),                INTENT(IN   ) :: state           ! Physics state variables
    TYPE(physics_buffer_desc), POINTER, INTENT(IN   ) :: hco_pbuf2d(:,:) ! Pointer to 2-D pbuf
    TYPE(MetState),                     INTENT(IN   ) :: State_Met       ! Meteorology State object
    INTEGER,                            INTENT(IN   ) :: iStep
!
! !OUTPUT PARAMETERS:
!
     TYPE(cam_in_t),                    INTENT(INOUT) :: cam_in                 ! import state
     REAL(r8),                          INTENT(  OUT) :: eflx(pcols,pver,pcnst) ! 3-D emissions in kg/m2/s
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
    INTEGER                                :: LCHNK
    INTEGER                                :: nY, nZ
    INTEGER                                :: I, M, N, J, L
    INTEGER                                :: RC           ! return code
    INTEGER                                :: tmpIdx       ! pbuf field id

    ! Logical
    LOGICAL                                :: rootChunk

    ! Objects
    TYPE(physics_buffer_desc), POINTER     :: pbuf_chnk(:) ! slice of pbuf in current chunk

    ! Real
    REAL(r8),                      POINTER :: pbuf_ik(:,:)  ! pointer to pbuf data (/pcols,pver/)
    REAL(r8), DIMENSION(state%NCOL,PVER+1) :: zint          ! Interface geopotential in km
    REAL(r8), DIMENSION(state%NCOL)        :: zsurf         ! Surface height
    REAL(r8)                               :: SCALFAC       ! Multiplying factor
    REAL(r8)                               :: megflx(pcols) ! For MEGAN emissions
    REAL(r8), PARAMETER                    :: m2km  = 1.e-3_r8

    ! Strings
    CHARACTER(LEN=255)                     :: fldname_ns   ! field name HCO_*

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
          ! If previous field name was not found, try with capitalized version
          ! to_upper is required because pFe /= PFE
          fldname_ns = 'HCO_' // to_upper(TRIM(tracerNames(N)))
          tmpIdx = pbuf_get_index(fldname_ns, RC)
       ENDIF

       IF ( tmpIdx < 0 .OR. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,'(a,a)') " CESMGC_Emissions_Calc: Field not found ", &
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
             Write(iulog,*) " CESMGC_Emissions_Calc: HEMCO emission flux is negative for ", &
                TRIM(cnst_name(M)), " with value ", MINVAL(eflx(:,:,M)), " at ", &
                MINLOC(eflx(:,:,M))
          ENDIF

          IF ( rootChunk .and. ( MAXVAL(eflx(1:nY,:,M)) > 0.0e+0_r8 ) ) THEN
             Write(iulog,'(a,a,a,a)') " CESMGC_Emissions_Calc: HEMCO flux ", &
                TRIM(fldname_ns), " added to ", TRIM(cnst_name(M))
             Write(iulog,'(a,a,E16.4)') " CESMGC_Emissions_Calc: Maximum flux ", &
                TRIM(fldname_ns), MAXVAL(eflx(1:nY,:,M))
          ENDIF
       ENDIF
    ENDDO

    !-----------------------------------------------------------------------
    ! Lightning NO emissions
    !-----------------------------------------------------------------------
    M = iNO

    ! prod_NO is in atom N cm^-3 s^-1 <=> molec cm^-3 s^-1
    ! We need to convert this to kg NO/m2/s
    ! Multiply by MWNO * BXHEIGHT * 1.0E+06 / AVO
    !           = mole/molec * kg NO/mole * m * cm^3/m^3
    ! cnst_mw(M) is in g/mole
    SCALFAC = cnst_mw(M) * 1.0E-03 * 1.0E+06 / AVO
    DO J = 1, nY
    DO L = 1, nZ
       eflx(J,L,M) = eflx(J,L,M)                      &
                   + prod_NO(J,L,LCHNK)               &
                     * State_Met%BXHEIGHT(1,J,nZ+1-L) &
                     * SCALFAC
    ENDDO
    ENDDO

#if defined( MODAL_AERO_4MODE )
    !-----------------------------------------------------------------------
    ! Aerosol emissions (dust + seasalt) ...
    !-----------------------------------------------------------------------
    call aero_model_emissions( state, cam_in )

    ! Since GEOS-Chem DST* aerosols are inherited from MAM's DST, we do not
    ! need to feed MAM dust emissions into the GEOS-Chem DST* constituents
    ! Same thing applies for sea salt.

    ! However, emissions of other aerosols (black carbon, organic carbon,
    ! secondary organic aerosols and SO4), as read by HEMCO, need to be fed
    ! to MAM's aerosol emission flux
    eflx(:nY,:nZ,iBC1)   = eflx(:nY,:nZ,iBCPI)
    eflx(:nY,:nZ,iBCPI)  = 0.0e+00_r8
    eflx(:nY,:nZ,iBC4)   = eflx(:nY,:nZ,iBCPO)
    eflx(:nY,:nZ,iBCPO)  = 0.0e+00_r8

    eflx(:nY,:nZ,iH2SO4) = eflx(:nY,:nZ,iSO4)
    eflx(:nY,:nZ,iSO4)   = 0.0e+00_r8

    ! For SOA emission, split evently GEOS-Chem SOAS emission into each 
    ! VBS bin.
    eflx(:nY,:nZ,iSOA11) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA12) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA21) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA22) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA31) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA32) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA41) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA42) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA51) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOA52) = eflx(:nY,:nZ,iSOAS) / 10.0e+00_r8
    eflx(:nY,:nZ,iSOAS)  = 0.0e+00_r8

    eflx(:nY,:nZ,iPOM1)  = eflx(:nY,:nZ,iOCPI)
    eflx(:nY,:nZ,iOCPI)  = 0.0e+00_r8
    eflx(:nY,:nZ,iPOM4)  = eflx(:nY,:nZ,iOCPO)
    eflx(:nY,:nZ,iOCPO)  = 0.0e+00_r8
#endif

    !-----------------------------------------------------------------------
    ! MEGAN emissions ...
    !-----------------------------------------------------------------------

    IF ( index_x2a_Fall_flxvoc > 0 .AND. shr_megan_mechcomps_n > 0 ) THEN
       ! set MEGAN fluxes 
       DO N = 1, shr_megan_mechcomps_n
          DO I = 1, nY
             megflx(I) = -cam_in%meganflx(I,N) * megan_wght_factors(N)
          ENDDO
          IF ( ( megan_indices_map(N) > 0 ) .AND. ( megan_wght_factors(N) > 0.0e+00_r8 ) ) THEN
             DO I = 1, nY
                cam_in%cflx(I,megan_indices_map(N)) = cam_in%cflx(I,megan_indices_map(N)) &
                                                    + megflx(I)
             ENDDO
          ENDIF
          ! output MEGAN emis fluxes to history
          CALL Outfld('MEG_'//TRIM(shr_megan_mechcomps(N)%name), megflx(:nY), nY, LCHNK)
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------
    ! Fire surface emissions if not elevated forcing
    !-----------------------------------------------------------------------

    CALL fire_emissions_srf( LCHNK, nY, cam_in%fireflx, cam_in%cflx )

    !-----------------------------------------------------------------------
    ! Apply CLM emissions (for elevated forcing)
    !-----------------------------------------------------------------------

    ! Compute geopotential height in km (needed for vertical distribution of
    ! fire emissions
    zsurf(:nY) = rga * state%phis(:nY)
    DO L = 1, nZ
       zint(:nY,L) = m2km * ( state%zi(:nY,L) + zsurf(:nY) )
    ENDDO
    L = nZ+1
    zint(:nY,L) = m2km * ( state%zi(:nY,L) + zsurf(:nY) )

    ! Distributed fire emissions if elevated forcing
    ! extfrc is in molec/cm3/s
    ! TMMF - vertical distribution of fire emissions is not implemented yet
    !CALL fire_emissions_vrt( nY, LCHNK, zint, cam_in%fireflx, cam_in%fireztop, extfrc )

    !-----------------------------------------------------------------------
    ! Add near-surface emissions to surface flux boundary condition
    !-----------------------------------------------------------------------
    cam_in%cflx(1:nY,:) = cam_in%cflx(1:nY,:) + eflx(1:nY,nZ,:)
    eflx(1:nY,nZ,:)     = 0.0e+00_r8

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

    IF ( ALLOCATED( megan_indices_map  ) ) DEALLOCATE( megan_indices_map )
    IF ( ALLOCATED( megan_wght_factors ) ) DEALLOCATE( megan_wght_factors )

  END SUBROUTINE CESMGC_Emissions_Final
!EOC
!------------------------------------------------------------------------------
!EOC
  END MODULE CESMGC_Emissions_Mod
