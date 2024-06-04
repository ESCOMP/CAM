! Module geoschem_emissions_mod contains routines which retrieve
! emission fluxes from HEMCO and transfers it back to the CESM-GC interface
! 07 Oct 2020 - T. M. Fritz   - Initial version
MODULE GeosChem_Emissions_Mod

  ! CAM modules
  use cam_abortutils,      only : endrun
  use cam_logfile,         only : iulog
  use chem_mods,           only : iFirstCnst
  use constituents,        only : pcnst, cnst_name
  use shr_kind_mod,        only : r8 => shr_kind_r8, shr_kind_cl
  use shr_megan_mod,       only : shr_megan_mechcomps, shr_megan_mechcomps_n
  use spmd_utils,          only : MasterProc

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: GC_Emissions_Init
  PUBLIC  :: GC_Emissions_Calc
  PUBLIC  :: GC_Emissions_Final

  ! Constituent number for NO
  INTEGER :: iNO

  ! Aerosol constituent number
  INTEGER :: iBC1
  INTEGER :: iBC4
  INTEGER :: iH2SO4

  INTEGER :: iBCPI
  INTEGER :: iBCPO
  INTEGER :: iOCPI
  INTEGER :: iOCPO
  INTEGER :: iSO4

  ! MEGAN Emissions
  INTEGER,  ALLOCATABLE :: megan_indices_map(:) 
  REAL(r8), ALLOCATABLE :: megan_wght_factors(:)

  ! Cache for is_extfrc?
  LOGICAL,  ALLOCATABLE :: pcnst_is_extfrc(:) ! no idea why the indexing is not 1:gas_pcnst or why iFirstCnst can be < 0

CONTAINS

  SUBROUTINE GC_Emissions_Init( )

    ! CAM modules
    use cam_history,         only : addfld, add_default, horiz_only
    use chem_mods,           only : adv_mass
    use constituents,        only : cnst_get_ind
    use fire_emissions,      only : fire_emissions_init
    use infnan,              only : NaN, assignment(=)
    use mo_chem_utls,        only : get_spc_ndx, get_extfrc_ndx
    use phys_control,        only : phys_getopts
    use physics_types,       only : physics_state

    ! Integers
    INTEGER                :: IERR
    INTEGER                :: N, II

    ! Logicals
    LOGICAL                :: history_aerosol
    LOGICAL                :: history_chemistry
    LOGICAL                :: history_cesm_forcing

    ! Strings
    CHARACTER(LEN=shr_kind_cl) :: SpcName
    CHARACTER(LEN=shr_kind_cl) :: Description

    ! Real
    REAL(r8)               :: MW

    CALL phys_getopts( history_aerosol_out      = history_aerosol,   &
                       history_chemistry_out    = history_chemistry, &
                       history_cesm_forcing_out = history_cesm_forcing )

    ! Get constituent index for NO
    CALL cnst_get_ind('NO', iNO, abort=.True.)

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
          IF ( TRIM(SpcName) == 'HCN' ) THEN
             SpcName = 'None'
             MW      = 27.025140_r8 ! Taken from pp_trop_strat_mam4_vbs
          ELSEIF ( TRIM(SpcName) == 'C2H4' ) THEN
             SpcName = 'None'
             MW      = 28.051600_r8 ! Taken from pp_trop_strat_mam4_vbs
          ENDIF
          !IF ( TRIM(SpcName) == 'MTERP' ) THEN
          !   SpcName = 'MTPA'
          !ELSEIF ( TRIM(SpcName) == 'BCARY' ) THEN
          !   SpcName = 'None'
          !   MW      = 204.342600_r8 ! Taken from pp_trop_strat_mam4_vbs
          !ELSEIF ( TRIM(SpcName) == 'CH3OH' ) THEN
          !   SpcName = 'MOH'
          !ELSEIF ( TRIM(SpcName) == 'C2H5OH' ) THEN
          !   SpcName = 'EOH'
          !ELSEIF ( TRIM(SpcName) == 'CH3CHO' ) THEN
          !   SpcName = 'ALD2'
          !ELSEIF ( TRIM(SpcName) == 'CH3COOH' ) THEN
          !   SpcName = 'ACTA'
          !ELSEIF ( TRIM(SpcName) == 'CH3COCH3' ) THEN
          !   SpcName = 'ACET'
          !ELSEIF ( TRIM(SpcName) == 'HCN' ) THEN
          !   SpcName = 'None'
          !   MW      = 27.025140_r8 ! Taken from pp_trop_strat_mam4_vbs
          !ELSEIF ( TRIM(SpcName) == 'C2H4' ) THEN
          !   SpcName = 'None'
          !   MW      = 28.051600_r8 ! Taken from pp_trop_strat_mam4_vbs
          !ELSEIF ( TRIM(SpcName) == 'C3H6' ) THEN
          !   SpcName = 'PRPE'
          !ELSEIF ( TRIM(SpcName) == 'BIGALK' ) THEN
          !   ! BIGALK = Pentane + Hexane + Heptane + Tricyclene
          !   SpcName = 'ALK4'
          !ELSEIF ( TRIM(SpcName) == 'BIGENE' ) THEN
          !   ! BIGENE = butene (C4H8)
          !   SpcName = 'PRPE' ! Lumped >= C3 alkenes
          !ELSEIF ( TRIM(SpcName) == 'TOLUENE' ) THEN
          !   SpcName = 'TOLU'
          !ENDIF

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

          IF ( history_chemistry ) THEN
             CALL Add_default('MEG_'//TRIM(SpcName), 1, ' ')
          ENDIF
       ENDDO
    ENDIF

    DO N = iFirstCnst, pcnst
       SpcName = TRIM(cnst_name(N))//'_XFRC'
       CALL Addfld( TRIM(SpcName), (/ 'lev' /), 'A', 'molec/cm3/s', &
          'External forcing for '//TRIM(cnst_name(N)))
       SpcName = TRIM(cnst_name(N))//'_CLXF'
       CALL Addfld( TRIM(SpcName), horiz_only, 'A', 'molec/cm2/s', &
          'Vertically-integrated external forcing for '//TRIM(cnst_name(N)))
       IF ( history_aerosol .OR. history_chemistry ) THEN
          CALL Add_Default( TRIM(SpcName), 1, ' ' )
       ENDIF
       IF ( history_cesm_forcing .AND. TRIM(cnst_name(N)) == 'NO2' ) THEN
          CALL Add_Default( TRIM(SpcName), 1, ' ' )
       ENDIF
       SpcName = TRIM(cnst_name(N))//'_CMXF'
       CALL Addfld( TRIM(SpcName), horiz_only, 'A', 'kg/m2/s', &
          'Vertically-integrated external forcing for '//TRIM(cnst_name(N)))
       IF ( history_aerosol .OR. history_chemistry ) THEN
          CALL Add_Default( TRIM(SpcName), 1, ' ' )
       ENDIF
       IF ( history_cesm_forcing .AND. TRIM(cnst_name(N)) == 'NO2' ) THEN
          CALL Add_Default( TRIM(SpcName), 1, ' ' )
       ENDIF
    ENDDO

    CALL Addfld( 'NO_Lightning', (/ 'lev' /), 'A','molec/cm3/s', &
          'lightning NO source' )

    !-----------------------------------------------------------------------
    ! ... Fire emissions
    !-----------------------------------------------------------------------
    CALL fire_emissions_init()

    ! Initialize pcnst_is_extfrc cache to avoid lengthy lookups in future timesteps
    ! on the get_extfrc_ndx routine. (hplin 1/20/23)
    if(.not. allocated(pcnst_is_extfrc)) then
      allocate(pcnst_is_extfrc(pcnst - iFirstCnst + 1))
    endif
    do n = iFirstCnst, pcnst
       pcnst_is_extfrc(n - iFirstCnst + 1) = (get_extfrc_ndx(trim(cnst_name(n))) > 0)
    enddo

  END SUBROUTINE GC_Emissions_Init

  SUBROUTINE GC_Emissions_Calc( state, hco_pbuf2d, State_Met, cam_in, eflx, iStep )
    ! Subroutine GC_Emissions_Calc retrieves emission fluxes
    ! from HEMCO and returns a 3-D array of emission flux to the CESM-GC
    ! interface. On top of passing data, this routine handles a number of checks.

    ! CAM modules
    use aero_model,          only : aero_model_emissions ! Aerosol emissions
    use cam_history,         only : outfld
    use camsrfexch,          only : cam_in_t
    use constituents,        only : cnst_get_ind, cnst_mw
    use fire_emissions,      only : fire_emissions_srf, fire_emissions_vrt ! Fire emissions
    use mo_lightning,        only : prod_NO! Lightning emissions
    use physconst,           only : rga, avogad
    use physics_buffer,      only : pbuf_get_index, pbuf_get_chunk
    use physics_buffer,      only : physics_buffer_desc, pbuf_get_field
    use physics_types,       only : physics_state
    use ppgrid,              only : pcols, pver, begchunk
    use srf_field_check,     only : active_Fall_flxvoc ! MEGAN emissions
    use string_utils,        only : to_upper

    ! GEOS-Chem modules
    use PhysConstants,       only : AVO, PI
    use State_Met_Mod,       only : MetState    

    TYPE(physics_state),                INTENT(IN   ) :: state           ! Physics state variables
    TYPE(physics_buffer_desc), POINTER, INTENT(IN   ) :: hco_pbuf2d(:,:) ! Pointer to 2-D pbuf
    TYPE(MetState),                     INTENT(IN   ) :: State_Met       ! Meteorology State object
    INTEGER,                            INTENT(IN   ) :: iStep

    TYPE(cam_in_t),                    INTENT(INOUT) :: cam_in                 ! import state
     REAL(r8),                          INTENT(  OUT) :: eflx(pcols,pver,pcnst) ! 3-D emissions in kg/m2/s

    ! Integers
    INTEGER                                :: LCHNK
    INTEGER                                :: nY, nZ
    INTEGER                                :: J, L, N
    INTEGER                                :: RC             ! return code
    INTEGER                                :: tmpIdx         ! pbuf field id

    INTEGER                                :: id_O3, id_HNO3 ! Species IDs for reuse

    ! Logical
    LOGICAL                                :: rootChunk

    ! Objects
    TYPE(physics_buffer_desc), POINTER     :: pbuf_chnk(:)   ! slice of pbuf in current chunk

    ! Real
    REAL(r8),                      POINTER :: pbuf_ik(:,:)   ! pointer to pbuf data (/pcols,pver/)
    REAL(r8),                      POINTER :: pbuf_i(:)      ! pointer to 2-D (1-D in CAM) data (/pcols/)
    REAL(r8), DIMENSION(state%NCOL,PVER+1) :: zint           ! Interface geopotential in km
    REAL(r8), DIMENSION(state%NCOL)        :: zsurf          ! Surface height
    REAL(r8)                               :: SCALFAC        ! Multiplying factor
    REAL(r8)                               :: megflx(pcols)  ! For MEGAN emissions
    REAL(r8), PARAMETER                    :: m2km  = 1.e-3_r8

    ! Strings
    CHARACTER(LEN=shr_kind_cl)             :: SpcName
    CHARACTER(LEN=shr_kind_cl)             :: fldname_ns     ! field name HCO_*

    
    ! Initialize pointers
    pbuf_chnk => NULL()
    pbuf_ik   => NULL()
    pbuf_i    => NULL()

    ! LCHNK: which chunk we have on this process
    LCHNK = state%LCHNK
    ! nY: number of atmospheric columns on this chunk
    nY    = state%NCOL
    nZ    = PVER
    rootChunk = ( MasterProc .AND. ( LCHNK.EQ.BEGCHUNK ) )

    ! Initialize emission flux
    eflx(:,:,:) = 0.0e+0_r8

    DO N = iFirstCnst, pcnst
       fldname_ns = 'HCO_'//TRIM(cnst_name(N))
       tmpIdx = pbuf_get_index(fldname_ns, RC)

       IF ( tmpIdx < 0 .OR. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,'(a,a)') " GC_Emissions_Calc: Field not found ", &
             TRIM(fldname_ns)
       ELSE
          ! This is already in chunk, retrieve it
          pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)

          ! Check if we need to get 3-D, or 2-D data
          IF (pcnst_is_extfrc(N - iFirstCnst + 1)) THEN
             CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)

             IF ( .NOT. ASSOCIATED(pbuf_ik) ) THEN ! Sanity check
                CALL ENDRUN("GC_Emissions_Calc: FATAL - tmpIdx > 0 but pbuf_ik not associated (E-1)")
             ENDIF

             eflx(1:nY,:nZ,N) = pbuf_ik(1:nY,:nZ)

             ! Reset pointers
             pbuf_ik   => NULL()
          ELSE ! 2-D
             CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

             IF ( .NOT. ASSOCIATED(pbuf_i) ) THEN ! Sanity check
                CALL ENDRUN("GC_Emissions_Calc: FATAL - tmpIdx > 0 but pbuf_i not associated (E-2)")
             ENDIF

             ! note: write to nZ level here as this is surface
             eflx(1:nY,nZ,N) = pbuf_i(1:nY)

             ! Reset pointers
             pbuf_i    => NULL()
          ENDIF

          pbuf_chnk => NULL()

          !IF ( MINVAL(eflx(:nY,:nZ,N)) < 0.0e+00_r8 ) THEN
          !   Write(iulog,*) " GC_Emissions_Calc: HEMCO emission flux is negative for ", &
          !      TRIM(cnst_name(N)), " with value ", MINVAL(eflx(:nY,:nZ,N)), " at ", &
          !      MINLOC(eflx(:nY,:nZ,N))
          !ENDIF

          IF ( rootChunk .AND. (iStep == 2) .AND. ( MAXVAL(eflx(:nY,:nZ,N)) > 0.0e+0_r8 ) ) THEN
             ! Only print this once
             Write(iulog,'(a,a,a,a)') " GC_Emissions_Calc: HEMCO flux ", &
                TRIM(fldname_ns), " added to ", TRIM(cnst_name(N))
             Write(iulog,'(a,a,E16.4)') " GC_Emissions_Calc: Maximum flux ", &
                TRIM(fldname_ns), MAXVAL(eflx(:nY,:nZ,N))
          ENDIF
       ENDIF
    ENDDO

    !-----------------------------------------------------------------------
    ! Deposition fluxes from HEMCO
    !-----------------------------------------------------------------------

    ! Deposition velocities in HEMCO are now handled within HEMCO_CESM for a
    ! hardcoded list of species, primarily for the SeaFlux extension.
    ! This is not to be confused with dry deposition fluxes which are not
    ! handled by HEMCO.

    ! Part 2: Handle special deposition fluxes for the ParaNOx extension
    ! for PAR_O3_DEP and PAR_HNO3_DEP
    CALL cnst_get_ind('O3', id_O3)
    CALL cnst_get_ind('HNO3', id_HNO3)

    ! write(iulog,*) 'id_O3, cnst_name, id_HNO3, cnst_name', id_O3, cnst_name(id_O3), id_HNO3, cnst_name(id_HNO3)

    tmpIdx = pbuf_get_index('HCO_PAR_O3_DEP', RC)
    IF(tmpIdx < 0 .OR. ( iStep == 1 )) then
        ! No ParaNOx dep flux for O3
    ELSE
        pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
        CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

        IF ( .NOT. ASSOCIATED(pbuf_i) ) THEN ! Sanity check
           CALL ENDRUN("GC_Emissions_Calc: FATAL - tmpIdx > 0 but pbuf_i not associated (2)")
        ENDIF

        ! apply loss flux to surface (level nZ)
        eflx(1:NY,nZ,id_O3) = eflx(1:NY,nZ,id_O3) - pbuf_i(1:nY)

        !IF ( MINVAL(eflx(:nY,nZ,id_O3)) < 0.0e+00_r8 ) THEN
        !   Write(iulog,*) " GC_Emissions_Calc: HEMCO sfc flux after ParaNOx is negative for O3 with value ", MINVAL(eflx(:nY,:nZ,id_O3)), " at ", &
        !      MINLOC(eflx(:nY,nZ,id_O3))
        !ENDIF

        IF ( rootChunk .and. ( MINVAL(pbuf_i(1:nY)) < 0.0e+0_r8 ) ) THEN
           Write(iulog,'(a,a,a,a)') " GC_Emissions_Calc: HEMCO dflx(paranox) O3 added to ", TRIM(cnst_name(id_O3))
           Write(iulog,'(a,a,E16.4)') " GC_Emissions_Calc: Minval dflx(paranox), eflx(sfc) O3 ", MINVAL(pbuf_i(1:nY)), MINVAL(eflx(:nY,nZ,id_O3))
        ENDIF

        ! Reset pointers
        pbuf_i   => NULL()
        pbuf_chnk => NULL()
    ENDIF

    tmpIdx = pbuf_get_index('HCO_PAR_HNO3_DEP', RC)
    IF(tmpIdx < 0 .OR. ( iStep == 1 )) then
        ! No ParaNOx dep flux for HNO3
    ELSE
        pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
        CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

        IF ( .NOT. ASSOCIATED(pbuf_i) ) THEN ! Sanity check
           CALL ENDRUN("GC_Emissions_Calc: FATAL - tmpIdx > 0 but pbuf_i not associated (3)")
        ENDIF

        eflx(1:NY,nZ,id_HNO3) = eflx(1:NY,nZ,id_HNO3) - pbuf_i(1:nY)

        !IF ( MINVAL(eflx(:nY,nZ,id_HNO3)) < 0.0e+00_r8 ) THEN
        !   Write(iulog,*) " GC_Emissions_Calc: HEMCO sfc flux after ParaNOx is negative for HNO3 with value ", MINVAL(eflx(:nY,nZ,id_HNO3)), " at ", &
        !      MINLOC(eflx(:nY,nZ,id_HNO3))
        !ENDIF

        IF ( rootChunk .and. ( MINVAL(pbuf_i(1:nY)) < 0.0e+0_r8 ) ) THEN
           Write(iulog,'(a,a,a,a)') " GC_Emissions_Calc: HEMCO dflx(paranox) HNO3 added to ", TRIM(cnst_name(id_HNO3))
           Write(iulog,'(a,a,E16.4)') " GC_Emissions_Calc: Minval dflx(paranox), eflx(sfc) HNO3 ", MINVAL(pbuf_i(1:nY)), MINVAL(eflx(:nY,nZ,id_HNO3))
        ENDIF

        ! Reset pointers
        pbuf_i   => NULL()
        pbuf_chnk => NULL()
    ENDIF

#if defined( MODAL_AERO )

    !-----------------------------------------------------------------------
    ! Aerosol emissions (dust + seasalt) ...
    !-----------------------------------------------------------------------
    call aero_model_emissions( state, cam_in )

    ! Since GEOS-Chem DST* aerosols are inherited from MAM's DST, we do not
    ! need to feed MAM dust emissions into the GEOS-Chem DST* constituents
    ! Same thing applies for sea salt.

    ! HEMCO aerosol emissions are fed to MAM through the HEMCO_Config.rc
    ! where all GEOS-Chem aerosols (BCPI, BCPO, OCPI, OCPO, SO4) have been
    ! replaced with the corresponding MAM aerosols

#endif

    ! Output fields before lightning NO emissions are applied to eflx
    ! Make sure that we do not include surface emissions in the diagnostics!
    DO N = iFirstCnst, pcnst
       SpcName = TRIM(cnst_name(N))//'_XFRC'
       ! Convert from kg/m2/s to molec/cm3/s
       ! Note 1: cnst_mw is in kg/kmole
       ! Note 2: avogad is in molecules/kmole
       CALL Outfld( TRIM(SpcName), eflx(:nY,:nZ,N) / State_Met%BXHEIGHT(1,:nY,nZ:1:-1) * 1.0E-06 / cnst_mw(N) * avogad, nY, LCHNK )

       SpcName = TRIM(cnst_name(N))//'_CLXF'
       ! Convert from kg/m2/s to molec/cm2/s
       ! Note 1: cnst_mw is in kg/kmole
       ! Note 2: avogad is in molecules/kmole
       CALL Outfld( TRIM(SpcName), SUM(eflx(:nY,:nZ-1,N), DIM=2) * 1.0E-04 / cnst_mw(N) * avogad, nY, LCHNK )

       SpcName = TRIM(cnst_name(N))//'_CMXF'
       CALL Outfld( TRIM(SpcName), SUM(eflx(:nY,:nZ-1,N), DIM=2), nY, LCHNK )
    ENDDO

    !-----------------------------------------------------------------------
    ! Lightning NO emissions
    !-----------------------------------------------------------------------
    N = iNO

    ! prod_NO is in atom N cm^-3 s^-1 <=> molec cm^-3 s^-1
    ! We need to convert this to kg NO/m2/s
    ! Multiply by MWNO * BXHEIGHT * 1.0E+06 / AVO
    !           = mole/molec * kg NO/mole * m * cm^3/m^3
    ! cnst_mw(N) is in g/mole
    SCALFAC = cnst_mw(N) * 1.0E-03 * 1.0E+06 / AVO
    DO J = 1, nY
    DO L = 1, nZ
       eflx(J,L,N) = eflx(J,L,N)                      &
                   + prod_NO(J,L,LCHNK)               &
                     * State_Met%BXHEIGHT(1,J,nZ+1-L) &
                     * SCALFAC
    ENDDO
    ENDDO

    CALL Outfld( 'NO_Lightning', prod_NO(:nY,:nZ,LCHNK), nY, LCHNK )

    !-----------------------------------------------------------------------
    ! MEGAN emissions ...
    !-----------------------------------------------------------------------

    IF ( active_Fall_flxvoc .AND. shr_megan_mechcomps_n > 0 ) THEN
       ! set MEGAN fluxes 
       DO N = 1, shr_megan_mechcomps_n
          DO J = 1, nY
             megflx(J) = -cam_in%meganflx(J,N) * megan_wght_factors(N)
          ENDDO
          IF ( ( megan_indices_map(N) > 0 ) .AND. ( megan_wght_factors(N) > 0.0e+00_r8 ) ) THEN
             DO J = 1, nY
                cam_in%cflx(J,megan_indices_map(N)) = cam_in%cflx(J,megan_indices_map(N)) &
                                                    + megflx(J)
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

    ! Near-surface emissions are now emitted directly to GEOS-Chem Species array
    ! for consistency with CAM-chem implementation of HEMCO
    ! (but not with GEOS-Chem standalone, where fluxes are mixed by the turbulence routines)
    ! Refer to discussion here: https://github.com/ESCOMP/CAM/pull/560#discussion_r1084559191
    !
    ! To replicate old behavior, uncomment these two lines below:
    ! cam_in%cflx(1:nY,:) = cam_in%cflx(1:nY,:) + eflx(1:nY,nZ,:)
    ! eflx(1:nY,nZ,:)     = 0.0e+00_r8

  END SUBROUTINE GC_Emissions_Calc

  SUBROUTINE GC_Emissions_Final

    IF ( ALLOCATED( megan_indices_map  ) ) DEALLOCATE( megan_indices_map )
    IF ( ALLOCATED( megan_wght_factors ) ) DEALLOCATE( megan_wght_factors )

  END SUBROUTINE GC_Emissions_Final

END MODULE GeosChem_Emissions_Mod
