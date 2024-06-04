module chemistry

  ! CAM modules
  use cam_abortutils,      only : endrun
  use cam_logfile,         only : iulog
  use chem_mods,           only : nTracersMax, nTracers, tracerNames
  use chem_mods,           only : gas_pcnst, adv_mass, ref_MMR, iFirstCnst
  use chem_mods,           only : nSlsMax, nSls, slsNames, nSlvd, slvd_Lst
  use chem_mods,           only : nAerMax, nAer, aerNames, aerAdvMass
  use chem_mods,           only : map2GC, map2GCinv, map2GC_Sls
  use chem_mods,           only : mapCnst, map2chm, map2MAM4
  use constituents,        only : pcnst, cnst_add, cnst_get_ind, cnst_name
  use mo_tracname,         only : solsym
  use physics_buffer,      only : physics_buffer_desc
  use physics_types,       only : physics_state, physics_ptend, physics_ptend_init
  use ppgrid,              only : begchunk, endchunk, pcols, pver, pverp
  use shr_const_mod,       only : molw_dryair=>SHR_CONST_MWDAIR
  use shr_drydep_mod,      only : nddvels => n_drydep, drydep_list
  use shr_kind_mod,        only : r8 => shr_kind_r8, shr_kind_cl
  use spmd_utils,          only : MasterProc, myCPU=>Iam, nCPUs=>npes
  use string_utils,        only : to_upper
#if defined( MODAL_AERO )
  use modal_aero_data,     only : ntot_amode
#endif
  
  ! GEOS-Chem derived types
  USE DiagList_Mod,         ONLY : DgnList          ! Diagnostics list object
  use GeosChem_History_Mod, ONLY : HistoryConfigObj ! History diagnostic object
  USE Input_Opt_Mod,        ONLY : OptInput         ! Input Options
  USE Species_Mod,          ONLY : Species          ! Species object
  USE State_Chm_Mod,        ONLY : ChmState         ! Chemistry State object
  USE State_Diag_Mod,       ONLY : DgnState         ! Diagnostics State object
  USE State_Grid_Mod,       ONLY : GrdState         ! Grid State object
  USE State_Met_Mod,        ONLY : MetState         ! Meteorology State object
  USE TaggedDiagList_Mod,   ONLY : TaggedDgnList    ! Ragged diagnostics list

  ! GEOS-Chem utilities
  USE ErrCode_Mod,         ONLY : GC_SUCCESS, GC_FAILURE
  USE ErrCode_Mod,         ONLY : GC_Error, GC_CheckVar, GC_Warning
  USE Error_Mod,           ONLY : Error_Stop
  USE Precision_Mod,       ONLY : fp, f4                 ! Flexible precision

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! Public interfaces
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_readnl                    ! read chem namelist
  public :: chem_emissions
  public :: chem_timestep_init

  ! Location of valid geoschem_config.yml and species_database.yml
  ! Use local files in run folder
  CHARACTER(LEN=500) :: gcConfig = 'geoschem_config.yml'
  CHARACTER(LEN=500) :: speciesDB = 'species_database.yml'

  ! Location of chemistry input
  CHARACTER(LEN=shr_kind_cl) :: geoschem_cheminputs

  ! Debugging
  LOGICAL :: debug = .TRUE.

  ! Derived type objects
  TYPE(OptInput)                     :: Input_Opt       ! Input Options object
  TYPE(ChmState),ALLOCATABLE         :: State_Chm(:)    ! Chemistry State object
  TYPE(DgnState),ALLOCATABLE         :: State_Diag(:)   ! Diagnostics State object
  TYPE(GrdState),ALLOCATABLE         :: State_Grid(:)   ! Grid State object
  TYPE(MetState),ALLOCATABLE         :: State_Met(:)    ! Meteorology State object
  TYPE(DgnList )                     :: Diag_List       ! Diagnostics list object
  TYPE(TaggedDgnList )               :: TaggedDiag_List ! Tagged diagnostics list object
  TYPE(HistoryConfigObj), POINTER    :: HistoryConfig   ! HistoryConfig object for History diagn.
  type(physics_buffer_desc), POINTER :: hco_pbuf2d(:,:) ! Pointer to 2D pbuf

  ! Mimic code in sfcvmr_mod.F90
  TYPE :: SfcMrObj
     CHARACTER(LEN=63)         :: FldName        ! Field name
     INTEGER                   :: SpcID          ! ID in species database
     TYPE(SfcMrObj), POINTER   :: Next           ! Next element in list
  END TYPE SfcMrObj

  ! Heat of linked list with SfcMrObj objects
  TYPE(SfcMrObj),    POINTER   :: SfcMrHead => NULL()

  ! Field prefix
  CHARACTER(LEN=63), PARAMETER :: Prefix_SfcVMR = 'VMR_'

  ! Indices of critical species in GEOS-Chem
  INTEGER                    :: iH2O, iO3, iCO2, iSO4
  INTEGER                    :: iO, iH, iO2
  REAL(r8)                   :: MWO3
  ! Indices of critical species in the constituent list
  INTEGER                    :: cQ, cH2O, cH2SO4
  ! Indices of critical species in the solsym list
  INTEGER                    :: l_H2SO4, l_SO4
#if defined( MODAL_AERO )
  INTEGER, ALLOCATABLE       :: iSulf(:)
#endif

  ! Indices in the physics buffer
  INTEGER                    :: NDX_PBLH      ! PBL height [m]
  INTEGER                    :: NDX_FSDS      ! Downward shortwave flux at surface [W/m2]
  INTEGER                    :: NDX_CLDTOP    ! Cloud top height [index]
  INTEGER                    :: NDX_CLDFRC    ! Cloud fraction [-]
  INTEGER                    :: NDX_PRAIN     ! Rain production rate [kg/kg/s]
  INTEGER                    :: NDX_NEVAPR    ! Total rate of precipitation evaporation  [kg/kg/s]
  INTEGER                    :: NDX_LSFLXPRC  ! Large-scale precip. at interface (liq + snw) [kg/m2/s]
  INTEGER                    :: NDX_LSFLXSNW  ! Large-scale precip. at interface (snow only) [kg/m2/s]
  INTEGER                    :: NDX_CMFDQR    ! Convective total precip. production rate [kg/kg/s]

  ! Get constituent indices
  INTEGER                    :: ixCldLiq  ! Cloud liquid water
  INTEGER                    :: ixCldIce  ! Cloud ice
  INTEGER                    :: ixNDrop   ! Cloud droplet number index

  ! ghg
  LOGICAL                    :: ghg_chem = .false.  ! .true. => use ghg chem package
  CHARACTER(len=shr_kind_cl) :: bndtvg = ' '   ! pathname for greenhouse gas loss rate
  CHARACTER(len=shr_kind_cl) :: h2orates = ' ' ! pathname for greenhouse gas (lyman-alpha H2O loss)

  ! Strings
  CHARACTER(LEN=shr_kind_cl) :: ThisLoc
  CHARACTER(LEN=shr_kind_cl) :: ErrMsg

  ! For dry deposition
  character(len=shr_kind_cl) :: depvel_lnd_file = 'depvel_lnd_file'


contains

  !================================================================================================
  ! function chem_is
  !================================================================================================
  function chem_is (name) result (chem_name_is)

    ! CAM modules
    use string_utils, only : to_lower

    character(len=*), intent(in) :: name
    logical :: chem_name_is

    chem_name_is = (( to_lower(name) == 'geoschem'  ) .or. &
                    ( to_lower(name) == 'geos-chem' ))

  end function chem_is

  !================================================================================================
  ! subroutine chem_register
  !================================================================================================
  subroutine chem_register

    ! CAM modules
    use chem_mods,           only : drySpc_ndx
    use mo_chem_utls,        only : get_spc_ndx
    use physconst,           only : MWDry
    use physics_buffer,      only : pbuf_add_field, dtype_r8
    use short_lived_species, only : Register_Short_Lived_Species
#if defined( MODAL_AERO )
    use aero_model,          only : aero_model_register
    use modal_aero_data,     only : nspec_max
    use modal_aero_data,     only : ntot_amode, nspec_amode
    use rad_constituents,    only : rad_cnst_get_info
#endif

    ! GEOS-Chem interface modules in CAM
    use mo_sim_dat,          only : set_sim_dat

    ! GEOS-Chem modules
    use GC_Environment_Mod,  ONLY : GC_Init_Grid
    use Input_Opt_Mod,       only : Set_Input_Opt,  Cleanup_Input_Opt
    use State_Chm_Mod,       only : Init_State_Chm, Cleanup_State_Chm, Ind_
    use State_Grid_Mod,      only : Init_State_Grid, Cleanup_State_Grid

    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents for chemistry
    !
    !-----------------------------------------------------------------------
    ! Need to generate a temporary species database
    TYPE(ChmState)                 :: SC
    TYPE(GrdState)                 :: SG
    TYPE(OptInput)                 :: IO
    TYPE(Species), POINTER         :: ThisSpc

    INTEGER                        :: I, N, M, L
    INTEGER                        :: nIgnored
    INTEGER                        :: tmpIdx
    REAL(r8)                       :: cptmp
    REAL(r8)                       :: MWTmp
    REAL(r8)                       :: qmin
    REAL(r8)                       :: refmmr, refvmr
    REAL(r8), ALLOCATABLE          :: slvd_refmmr(:)
    CHARACTER(LEN=128)             :: mixtype
    CHARACTER(LEN=128)             :: molectype
    CHARACTER(LEN=128)             :: lngName
    CHARACTER(LEN=64)              :: cnstName
    CHARACTER(LEN=64)              :: trueName
    CHARACTER(LEN=64)              :: aerName
    LOGICAL                        :: camout
    LOGICAL                        :: ic_from_cam2
    LOGICAL                        :: has_fixed_ubc
    LOGICAL                        :: has_fixed_ubflx

    INTEGER                        :: RC, IERR

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at GEOS-Chem (in chemistry/geoschem/chemistry.F90)'

    ! Initialize pointer
    ThisSpc => NULL()

    if (debug .and. masterproc) write(iulog,'(a)') 'chem_register: registering advected constituents for GEOS-Chem chemistry'

    ! SDE 2018-05-02: This seems to get called before anything else
    ! that includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    ! hplin 2020-05-16: Call set_sim_dat to populate chemistry constituent information
    ! from mo_sim_dat.F90 in other places. This is needed for HEMCO_CESM.
    CALL Set_sim_dat()

    ! Prevent Reporting
    IO%amIRoot             = .False.
    IO%thisCpu             = MyCPU

    CALL Set_Input_Opt( am_I_Root = MasterProc, &
                        Input_Opt = IO,         &
                        RC        = RC       )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Could not generate reference input options object!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Options needed by Init_State_Chm
    IO%ITS_A_FULLCHEM_SIM  = .True.
    IO%LLinoz              = .True.
    IO%LPRT                = .False.
    IO%N_Advect            = nTracers
    DO I = 1, nTracers
       IO%AdvectSpc_Name(I) = TRIM(tracerNames(I))
    ENDDO
    IO%SALA_rEdge_um(1)    = 0.01e+0_fp
    IO%SALA_rEdge_um(2)    = 0.50e+0_fp
    IO%SALC_rEdge_um(1)    = 0.50e+0_fp
    IO%SALC_rEdge_um(2)    = 8.00e+0_fp

    IO%SpcDatabaseFile     = TRIM(speciesDB)

    CALL Init_State_Grid( Input_Opt  = IO,  &
                          State_Grid = SG,  &
                          RC         = RC  )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_State_Grid"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    SG%NX = 1
    SG%NY = 1
    SG%NZ = 1

    CALL GC_Init_Grid( Input_Opt  = IO,  &
                       State_Grid = SG,  &
                       RC         = RC  )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error in GC_Init_Grid"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Init_CMN_SIZE"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Init_State_Chm( Input_Opt  = IO,  &
                         State_Chm  = SC,  &
                         State_Grid = SG,  &
                         RC         = RC  )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Init_State_Chm"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    iFirstCnst = -1
    mapCnst    = -1
    map2GC     = -1
    map2GCinv  = -1
    map2chm    = -1
    ref_MMR(:) = 0.0e+0_r8

    ! nTracersMax must be # advected species in geoschem_config.yml (nTracers) plus
    ! # aerosols (nAer) plus 1 (for CO2). It is set in chem_mods.F90.
    DO I = 1, nTracersMax
       IF ( I .LE. nTracers ) THEN
          cnstName    = to_upper(TRIM(tracerNames(I)))
          trueName    = cnstName
          N           = Ind_(cnstName)
          ThisSpc     => SC%SpcData(N)%Info
          lngName     = TRIM(ThisSpc%FullName)
          MWTmp       = REAL(ThisSpc%MW_g,r8)
          refvmr      = REAL(ThisSpc%BackgroundVV,r8)
          refmmr      = refvmr / (MWDry / MWTmp)
          ! Make sure that solsym is following the list of tracers as listed
          ! geoschem_config.yml
          IF ( to_upper(TRIM(tracerNames(I))) /= to_upper(TRIM(solsym(I))) ) THEN
             Write(iulog,*) "tracerNames (", TRIM(tracerNames(I)), ") /= solsym (", &
                   TRIM(solsym(I)), ")"
             CALL ENDRUN('Solsym must be following GEOS-Chem tracer. Check geoschem/mo_sim.dat')
          ENDIF
          ! Nullify pointer
          ThisSpc => NULL()
       ELSEIF ( I .LE. (nTracers + nAer) ) THEN
          ! Add MAM4 aerosols
          cnstName    = TRIM(aerNames(I - nTracers))
          trueName    = cnstName
          lngName     = cnstName
          MWTmp       = aerAdvMass(I - nTracers)
          refmmr      = 1.0e-38_r8
       ELSEIF ( I .EQ. (nTracers + nAer + 1) ) THEN
          ! Add CO2 (which is not a GEOS-Chem tracer)
          cnstName    = 'CO2'
          trueName    = cnstName
          lngName     = cnstName
          MWTmp       = 44.009800_r8
          refmmr      = 1.0e-38_r8
       ELSE
          cnstName    = TRIM(tracerNames(I))
          trueName    = cnstName
          lngName     = cnstName
          MWTmp       = 1000.0e+0_r8 * (0.001e+0_r8)
          refmmr      = 1.0e-38_r8
       ENDIF

       ! dummy value for specific heat of constant pressure (Cp)
       cptmp = 666._r8
       ! minimum mixing ratio
       qmin = 1.e-38_r8
       ! mixing ratio type
       mixtype = 'dry'
       ! Used for ionospheric WACCM (WACCM-X)
       molectype = 'minor'
       ! Is an output field (?)
       camout = .false.
       ! Not true for O2(1-delta) or O2(1-sigma)
       ic_from_cam2  = .true.
       ! Use a fixed value at the upper boundary
       has_fixed_ubc = .false.
       ! Use a fixed flux condition at the upper boundary
       has_fixed_ubflx = .false.

       ! TMMF - 8/20/2020
       ! Note: I had to modify the IC file to rename variables such as
       ! CH3COCH3 into ACET. Using that new IC file, we can thus remove
       ! the unnecessary special handlings.
       ! Another option would have been to modify cnst_add and read_inidat
       ! to use a load_name the first time IC are read. Constituent names
       ! would be stored in cnst_name, while read_inidat would load from
       ! load_name. load_name would be an optional argument to cnst_add, such
       ! that, by default, load_name = cnst_name.
       ! However, this would be tricky to handle with restart files that
       ! would save cnst_name rather than load_name.

       ! Special handlings
       IF ( cnstName == 'HCHO' ) THEN
           cnstName = 'CH2O'
       !ELSEIF ( cnstName == 'HNO4' ) THEN
       !    cnstName = 'HO2NO2'
       !ELSEIF ( cnstName == 'HNO2' ) THEN
       !    cnstName = 'HONO'
       !ELSEIF ( cnstName == 'ACET' ) THEN
       !    cnstName = 'CH3COCH3'
       !ELSEIF ( cnstName == 'ALD2' ) THEN
       !    cnstName = 'CH3CHO'
       !ELSEIF ( cnstName == 'PRPE' ) THEN
       !    cnstName = 'C3H6'
       !ELSEIF ( cnstName == 'MP'   ) THEN
       !    cnstName = 'CH3OOH'
       !ELSEIF ( cnstName == 'HAC'  ) THEN
       !    cnstName = 'HYAC'
       !ELSEIF ( cnstName == 'GLYC' ) THEN
       !    cnstName = 'GLYALD'
       !ELSEIF ( cnstName == 'MAP' ) THEN
       !    cnstName = 'CH3COOOH'
       !ELSEIF ( cnstName == 'EOH' ) THEN
       !    cnstName = 'C2H5OH'
       !ELSEIF ( cnstName == 'MGLY' ) THEN
       !    cnstName = 'CH3COCHO'
       !ELSEIF ( cnstName == 'GLYX' ) THEN
       !    cnstName = 'GLYOXAL'
       !ELSEIF ( cnstName == 'ACTA' ) THEN
       !    cnstName = 'CH3COOH'
       !ELSEIF ( cnstName == 'TOLU' ) THEN
       !    cnstName = 'TOLUENE'
       ENDIF

       CALL cnst_add( cnstName, MWtmp, cptmp, qmin, N,        &
                      readiv=ic_from_cam2, mixtype=mixtype,   &
                      cam_outfld=camout, molectype=molectype, &
                      fixed_ubc=has_fixed_ubc,                &
                      fixed_ubflx=has_fixed_ubflx,            &
                      longname=TRIM(lngName)                 )

       IF ( iFirstCnst < 0 ) iFirstCnst = N

       ref_MMR(N) = refmmr

       ! Add to GC mapping. When starting a timestep, we will want to update the
       ! concentration of State_Chm(x)%Species(m)%Conc(1,iCol,iLev) with data from
       ! constituent n
       M = Ind_(TRIM(trueName))
       IF ( M > 0 ) THEN
          ! Map constituent onto GEOS-Chem tracer as indexed in State_Chm(LCHNK)%Species
          map2GC(N)    = M
          ! Map GEOS-Chem tracer onto constituent
          map2GCinv(M) = N
       ENDIF
       ! Map constituent onto chemically-active species (aka as indexed in solsym)
       M = get_spc_ndx(TRIM(trueName), ignore_case=.true.)
       IF ( M > 0 ) THEN
          mapCnst(N) = M
       ENDIF
    ENDDO

    ! Now unadvected species
    map2GC_Sls = 0
    ALLOCATE(slvd_refmmr(nslvd), STAT=IERR)
    slvd_refmmr(:) = 0.0e+0_r8
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate map2MAM4')
    DO I = 1, nSlvd
       N = Ind_(slsNames(I))
       IF ( N .GT. 0 ) THEN
          ThisSpc         => SC%SpcData(N)%Info
          MWTmp           = REAL(ThisSpc%MW_g,r8)
          refvmr          = REAL(ThisSpc%BackgroundVV,r8)
          lngName         = TRIM(ThisSpc%FullName)
          slvd_refmmr(I)  = refvmr / (MWDry / MWTmp)
          map2GC_Sls(I)   = N
          ThisSpc         => NULL()
       ENDIF
    ENDDO
    CALL Register_Short_Lived_Species(slvd_refmmr)
    DEALLOCATE(slvd_refmmr)
    ! More information:
    ! http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/node5.html

    if (debug .and. masterproc) write(iulog,'(a,i4,a)') 'chem_register: looping over gas_pcnst (length', gas_pcnst, ') to map solsym onto GEOS-Chem species'

    DO N = 1, gas_pcnst
       ! Map solsym onto GEOS-Chem species
       map2chm(N) = Ind_(TRIM(solsym(N)))
       IF ( map2chm(N) < 0 ) THEN
          ! This is not a GEOS-Chem species and we thus map to constituents list.
          ! Most likely, these will be MAM aerosols
          ! We store the index as the opposite to not confuse with GEOS-Chem
          ! indices.
          CALL cnst_get_ind(TRIM(solsym(N)), I, abort=.True.)
          map2chm(N) = -I
          if (debug .and. masterproc) write(iulog,'(a,a,a,I4,a,I4)') ' -> solsym species ', trim(solsym(N)), ' (index ', N, ') is not a GEOS-Chem species. Mapping to negative constituent index: ', map2chm(N)
       ELSE
          if (debug .and. masterproc) write(iulog,'(a,a,a,I4,a,I4)') ' -> solsym species ', trim(solsym(N)), ' (index ', N, ') mapped to GEOS-Chem species ', map2chm(N)
       ENDIF
    ENDDO
    ! Get constituent index of specific humidity
    CALL cnst_get_ind('Q',     cQ,     abort=.True.)
    CALL cnst_get_ind('H2O',   cH2O,   abort=.True.)
    CALL cnst_get_ind('H2SO4', cH2SO4, abort=.True.)
 
    !------------------------------------------------------------
    ! Get mapping between dry deposition species and species set
    !------------------------------------------------------------
    
    nIgnored = 0

    if (debug .and. masterproc) write(iulog,'(a,i4,a)') 'chem_register: looping over gas dry deposition list with ', nddvels, ' species'

    DO N = 1, nddvels

       ! The species names need to be convert to upper case as,
       ! for instance, BR2 != Br2
       drySpc_ndx(N) = get_spc_ndx( to_upper(drydep_list(N)), ignore_case=.true. )

       if (debug .and. masterproc) write(iulog,'(a,a,a,i4,a,i4)') ' -> species ', trim(drydep_list(N)), ' in dry deposition list at index ', N, ' maps to species in solsym at index ', drySpc_ndx(N)

       IF ( MasterProc .and. ( drySpc_ndx(N) < 0 ) ) THEN
          Write(iulog,'(a,a)') ' ## Ignoring dry deposition of ', &
                               TRIM(drydep_list(N))
          nIgnored = nIgnored + 1
       ENDIF
    ENDDO

    IF ( MasterProc .AND. ( nIgnored > 0 ) ) THEN
       Write(iulog,'(a,a)') ' The species listed above have dry', &
         ' deposition turned off for one of the following reasons:'
       Write(iulog,'(a)') '  - They are not present in the GEOS-Chem tracer list.'
       Write(iulog,'(a)') '  - They have a synonym (e.g. CH2O and HCHO).'
    ENDIF

#if defined( MODAL_AERO_4MODE )
    ! add fields to pbuf needed by aerosol models
    CALL aero_model_register()

    ! Mode                   | \sigma_g | Dry diameter (micrometers)
    ! -----------------------|----------|--------------------------
    ! a2 - Aitken mode       |   1.6    | 0.015 - 0.053
    ! a1 - Accumulation mode |   1.8    | 0.058 - 0.27
    ! a3 - Coarse mode       |   1.8    | 0.80  - 3.65
    ! a4 - Primary carbon    |   1.6    | 0.039 - 0.13
    ! -----------------------|----------|--------------------------
    ! Ref: Liu, Xiaohong, et al. "Toward a minimal representation of aerosols in
    ! climate models: Description and evaluation in the Community Atmosphere
    ! Model CAM5." Geoscientific Model Development 5.3 (2012): 709.

    ALLOCATE(map2MAM4(nspec_max,ntot_amode), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate map2MAM4')

    ALLOCATE(iSulf(ntot_amode), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate iSulf')

    ! Initialize indices
    map2MAM4(:,:) = -1
    iSulf(:)      = -1

    ! ewl notes: xname_massptr returns a name. The select case subsets characters? e.g. 1:3, 4:5, 5:6.
    ! so want to get a name give an L and M. Need anything else???

    DO M = 1, ntot_amode
       DO L = 1, nspec_amode(M)
          call rad_cnst_get_info(0,M,L,spec_name=aername)
          SELECT CASE ( to_upper(aername(:3)) )
             CASE ( 'BC_' )
                SELECT CASE ( to_upper(aername(4:5)) )
                   CASE ( 'A1' )
                       CALL cnst_get_ind( 'BCPI', map2MAM4(L,M) )
                   CASE ( 'A4' )
                       CALL cnst_get_ind( 'BCPO', map2MAM4(L,M) )
                END SELECT
             CASE ( 'DST' )
                SELECT CASE ( to_upper(aername(5:6)) )
                   ! DST1 - Dust aerosol, Reff = 0.7 micrometers
                   ! DST2 - Dust aerosol, Reff = 1.4 micrometers
                   ! DST3 - Dust aerosol, Reff = 2.4 micrometers
                   ! DST4 - Dust aerosol, Reff = 4.5 micrometers
                   CASE ( 'A1' )
                       CALL cnst_get_ind( 'DST1', map2MAM4(L,M) )
                   CASE ( 'A2' )
                       CALL cnst_get_ind( 'DST1', map2MAM4(L,M) )
                   CASE ( 'A3' )
                       CALL cnst_get_ind( 'DST4', map2MAM4(L,M) )
                END SELECT
             !CASE ( 'SOA' )
             !   CALL cnst_get_ind( 'SOAS', map2MAM4(L,M) )
             CASE ( 'SO4' )
                CALL cnst_get_ind( 'SO4', map2MAM4(L,M) )
                iSulf(M) = L
             CASE ( 'NCL' )
                SELECT CASE ( to_upper(aername(5:6)) )
                   ! SALA - Fine (0.01-0.05 micros) sea salt aerosol
                   ! SALC - Coarse (0.5-8 micros) sea salt aerosol
                   CASE ( 'A1' )
                      CALL cnst_get_ind( 'SALA', map2MAM4(L,M) )
                   CASE ( 'A2' )
                      CALL cnst_get_ind( 'SALA', map2MAM4(L,M) )
                   CASE ( 'A3' )
                      CALL cnst_get_ind( 'SALC', map2MAM4(L,M) )
                END SELECT
             CASE ( 'POM' )
                SELECT CASE ( to_upper(aername(5:6)) )
                   CASE ( 'A1' )
                      CALL cnst_get_ind( 'OCPI', map2MAM4(L,M) )
                   CASE ( 'A4' )
                      CALL cnst_get_ind( 'OCPO', map2MAM4(L,M) )
                END SELECT
          END SELECT
       ENDDO
    ENDDO

#endif

    ! Print summary
    IF ( MasterProc ) THEN
       Write(iulog,'(/, a)') '### Summary of GEOS-Chem species (end of chem_register): '
       Write(iulog,'( a)') REPEAT( '-', 50 )
       Write(iulog,'( a)') '+ List of advected species: '
       Write(iulog,100) 'ID', 'Tracer', ''!'Dry deposition (T/F)'
       DO N = 1, nTracers
          Write(iulog,120) N, TRIM(tracerNames(N))!, ANY(drySpc_ndx .eq. N)
       ENDDO
       IF ( nAer > 0 ) THEN
          Write(iulog,'(/, a)') '+ List of aerosols: '
          Write(iulog,110) 'ID', 'MAM4 Aerosol'
          DO N = 1, nAer
             Write(iulog,130) N, TRIM(aerNames(N))
          ENDDO
       ENDIF
       Write(iulog,'(/, a)') '+ List of short-lived species: '
       DO N = 1, nSls
          Write(iulog,130) N, TRIM(slsNames(N))
       ENDDO
    ENDIF

100 FORMAT( 1x, A3, 3x, A10, 1x, A25 )
110 FORMAT( 1x, A3, 3x, A15 )
!120 FORMAT( 1x, I3, 3x, A10, 1x, L15 )
120 FORMAT( 1x, I3, 3x, A10 )
130 FORMAT( 1x, I3, 3x, A10 )

    ! Clean up
    Call Cleanup_State_Chm ( SC, RC )
    Call Cleanup_State_Grid( SG, RC )
    Call Cleanup_Input_Opt ( IO, RC )

    ! Add data for HEMCO extensions to buffers
    call pbuf_add_field('HCO_IN_JNO2', 'global', dtype_r8, (/pcols/), tmpIdx)
    call pbuf_add_field('HCO_IN_JOH', 'global',  dtype_r8, (/pcols/), tmpIdx)

    if (debug .and. masterproc) write(iulog,'(a)') 'chem_register: advected constituent registration for GEOS-Chem chemistry complete '

  end subroutine chem_register

  !================================================================================================
  ! subroutine chem_readnl
  !================================================================================================
  subroutine chem_readnl(nlfile)

    ! CAM modules
    use cam_abortutils,  only : endrun
    use chem_mods,       only : drySpc_ndx
    use gas_wetdep_opts, only : gas_wetdep_readnl
    use gckpp_Model,     only : nSpec, Spc_Names
    use namelist_utils,  only : find_group_name
    use mo_lightning,    only : lightning_readnl
    use spmd_utils,      only : mpicom, masterprocid, mpi_success
    use spmd_utils,      only : mpi_character, mpi_integer, mpi_logical
    use units,           only : getunit, freeunit
#if defined( MODAL_AERO )
    use aero_model,      only : aero_model_readnl
    use dust_model,      only : dust_readnl
#endif
    ! For dry deposition on unstructured grids
    use mo_drydep,       only : drydep_srf_file

    ! args
    CHARACTER(LEN=*), INTENT(IN) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    INTEGER                      :: I, N
    INTEGER                      :: UNITN, IERR, RC
    CHARACTER(LEN=500)           :: line
    CHARACTER(LEN=63)            :: substrs(2)
    CHARACTER(LEN=*), PARAMETER  :: subname = 'chem_readnl'
    LOGICAL                      :: validSLS, v_bool

    namelist /chem_inparm/ depvel_lnd_file
    namelist /chem_inparm/ drydep_srf_file

    ! ghg chem
    namelist /chem_inparm/ bndtvg, h2orates, ghg_chem

    if (debug .and. masterproc) write(iulog,'(a)') 'chem_readnl: reading namelists for GEOS-Chem chemistry'

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ALLOCATE(drySpc_ndx(nddvels), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate drySpc_ndx')

#if defined( MODAL_AERO_4MODE )
    ! Get names and molar weights of aerosols in MAM4

    nAer = 33

    aerNames(:nAer) = (/ 'bc_a1          ','bc_a4          ','dst_a1         ', &
                         'dst_a2         ','dst_a3         ','ncl_a1         ', &
                         'ncl_a2         ','ncl_a3         ','num_a1         ', &
                         'num_a2         ','num_a3         ','num_a4         ', &
                         'pom_a1         ','pom_a4         ','so4_a1         ', &
                         'so4_a2         ','so4_a3         ','soa1_a1        ', &
                         'soa1_a2        ','soa2_a1        ','soa2_a2        ', &
                         'soa3_a1        ','soa3_a2        ','soa4_a1        ', &
                         'soa4_a2        ','soa5_a1        ','soa5_a2        ', &
                         'H2SO4          ','SOAG0          ','SOAG1          ', &
                         'SOAG2          ','SOAG3          ','SOAG4          ' /)

    aerAdvMass(:nAer) = (/ 12.011000_r8,    12.011000_r8,    135.064039_r8,  &
                           135.064039_r8,   135.064039_r8,   58.442468_r8,   &
                           58.442468_r8,    58.442468_r8,    1.007400_r8,    &
                           1.007400_r8,     1.007400_r8,     1.007400_r8,    &
                           12.011000_r8,    12.011000_r8,    115.107340_r8,  &
                           115.107340_r8,   115.107340_r8,   250.445000_r8,  &
                           250.445000_r8,   250.445000_r8,   250.445000_r8,  &
                           250.445000_r8,   250.445000_r8,   250.445000_r8,  &
                           250.445000_r8,   250.445000_r8,   250.445000_r8,  &
                           98.078400_r8,    250.445000_r8,   250.445000_r8,  &
                           250.445000_r8,   250.445000_r8,   250.445000_r8  /)

    CALL aero_model_readnl(nlfile)
    CALL dust_readnl(nlfile)
#endif

    DO I = (nAer+1), nAerMax
       aerNames(I)   = 'EMPTY_AER      '
       aerAdvMass(I) = -1.00_r8
    ENDDO

    CALL gas_wetdep_readnl(nlfile)

    CALL lightning_readnl(nlfile)

    CALL geoschem_readnl(nlfile)

    IF ( MasterProc ) THEN

       Write(iulog,'(/,a)') REPEAT( '=', 50 )
       Write(iulog,'(a)') REPEAT( '=', 50 )
       Write(iulog,'(a)') 'This is the GEOS-CHEM / CESM interface'
       Write(iulog,'(a)') REPEAT( '=', 50 )
       Write(iulog,'(a)') ' + Routines written by Thibaud M. Fritz'
       Write(iulog,'(a)') ' + Laboratory for Aviation and the Environment,'
       Write(iulog,'(a)') ' + Department of Aeronautics and Astronautics,'
       Write(iulog,'(a)') ' + Massachusetts Institute of Technology'
       Write(iulog,'(a)') REPEAT( '=', 50 )

       Write(iulog,'(/,a,/)') 'Now defining GEOS-Chem tracers and dry deposition mapping...'

       !----------------------------------------------------------
       ! Read GEOS-Chem advected species from geoschem_config.yml
       !----------------------------------------------------------

       unitn = getunit()

       OPEN( unitn, FILE=TRIM(gcConfig), STATUS='OLD', IOSTAT=IERR )
       IF (IERR .NE. 0) THEN
          CALL ENDRUN('chem_readnl: ERROR opening geoschem_config.yml')
       ENDIF

       ! Find the transported species section
       DO
          READ( unitn, '(a)', IOSTAT=IERR ) line
          IF ( IERR .NE. 0 ) CALL ENDRUN('chem_readnl: error finding adv spc list')
          LINE = ADJUSTL( ADJUSTR( LINE ) )
          IF ( INDEX( LINE, 'transported_species' ) > 0 ) EXIT
       ENDDO

       if (debug) write(iulog,'(a)') 'chem_readnl: reading advected species list from geoschem_config.yml'

       ! Read in all advected species names and add them to tracer names list
       nTracers = 0
       DO WHILE ( LEN_TRIM( line ) > 0 )
          READ(unitn,'(a)', IOSTAT=IERR) line
          IF ( IERR .NE. 0 ) CALL ENDRUN('chem_readnl: error setting adv spc list')
          line = ADJUSTL( ADJUSTR( line ) )
          IF ( INDEX( line, 'passive_species' ) > 0 ) EXIT
          IF ( INDEX( LINE, '-' ) > 0 ) THEN
             substrs(1) = LINE(3:)
             substrs(1) = ADJUSTL( ADJUSTR( substrs(1) ) )

             ! Remove quotes (i.e. 'NO' -> NO)
             I = INDEX( substrs(1), "'" )
             IF ( I > 0 ) THEN
                substrs(1) = substrs(1)(I+1:)
                I = INDEX( substrs(1), "'" )
                IF ( I > 0 ) substrs(1) = substrs(1)(1:I-1)
             ENDIF

             nTracers = nTracers + 1
             tracerNames(nTracers) = TRIM(substrs(1))

             write(iulog,'(a,i4,a,a)') ' ', nTracers, ' ', TRIM(substrs(1))
          ENDIF
       ENDDO
       CLOSE(unitn)
       CALL freeunit(unitn)

       ! Assign remaining tracers dummy names
       DO I = (nTracers+1), nTracersMax
          WRITE(tracerNames(I),'(a,I0.4)') 'GCTRC_', I
       ENDDO

       !----------------------------------------------------------
       ! Now go through the KPP mechanism and add any species not
       ! implemented by the tracer list in geoschem_config.yml
       !----------------------------------------------------------
       
       IF ( nSpec > nSlsMax ) THEN
          CALL ENDRUN('chem_readnl: too many species - increase nSlsmax')
       ENDIF

       if (debug .and. masterproc) write(iulog,'(a)') 'chem_readnl: getting non-advected (short-lived) species list from KPP'
       if (debug .and. masterproc) write(iulog,'(a)') 'NOTE: does not include CO2 even if CO2 is not advected'

       nSls = 0
       DO I = 1, nSpec
          ! Get the name of the species from KPP
          line = ADJUSTL(TRIM(Spc_Names(I)))
          ! Only add short-lived KPP species, except from CO2
          validSLS = (( .NOT. ANY(TRIM(line) .EQ. tracerNames) ) &
                        .AND. TRIM(line) /= 'CO2' )
          IF ( validSLS ) THEN
             ! Genuine new short-lived species
             nSls = nSls + 1
             slsNames(nSls) = TRIM(line)
             write(iulog,'(a,i4,a,a)') ' ', nSls, ' ', TRIM(slsNames(nSls))
          ENDIF
       ENDDO

       unitn = getunit()
       OPEN( unitn, FILE=TRIM(nlfile), STATUS='old', IOSTAT=IERR )
       IF (IERR .NE. 0) THEN
          CALL ENDRUN('chem_readnl: ERROR opening '//TRIM(nlfile))
       ENDIF

       CALL find_group_name(unitn, 'chem_inparm', STATUS=IERR)
       IF (IERR == 0) THEN
          READ(unitn, chem_inparm, IOSTAT=IERR)
          IF (IERR /= 0) THEN
             CALL endrun('chem_readnl: ERROR reading namelist chem_inparm')
          ENDIF
       ENDIF
       CLOSE(unitn)
       CALL freeunit(unitn)

    ENDIF

    !----------------------------------------------------------
    ! Broadcast to all processors
    !----------------------------------------------------------
    CALL mpi_bcast(nTracers, 1, mpi_integer, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: nTracers')
    ENDIF
    CALL mpi_bcast(tracerNames, LEN(tracerNames(1))*nTracersMax, mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: tracerNames')
    ENDIF
    CALL mpi_bcast(nSls, 1, mpi_integer, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: nSls')
    ENDIF
    CALL mpi_bcast(slsNames, LEN(slsNames(1))*nSlsMax, mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: slsNames')
    ENDIF

    ! Broadcast namelist variables
    CALL mpi_bcast(depvel_lnd_file, LEN(depvel_lnd_file), mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: depvel_lnd_file')
    ENDIF
    CALL mpi_bcast(drydep_srf_file, LEN(drydep_srf_file), mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: drydep_srf_file')
    ENDIF
    CALL mpi_bcast(ghg_chem, 1, mpi_logical, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: ghg_chem')
    ENDIF
    CALL mpi_bcast(bndtvg, LEN(bndtvg), mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: bndtvg')
    ENDIF
    CALL mpi_bcast(h2orates, LEN(h2orates), mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: h2orates')
    ENDIF

    IF ( nSls .NE. nSlvd ) THEN
       write(iulog,'(a,i4)') 'nSlvd in geoschem/chem_mods.F90 does not match # non-advected KPP species. Set nSlvd to ', nSls
       CALL ENDRUN('Failure while allocating slvd_Lst')
    ENDIF

    ALLOCATE(slvd_Lst(nSlvd), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating slvd_Lst')
    DO I = 1, nSls
       slvd_Lst(I) = TRIM(slsNames(I))
    ENDDO

    if (debug .and. masterproc) write(iulog,'(a)') 'chem_readnl: reading GEOS-Chem chemistry namelists complete'

  end subroutine chem_readnl

  !================================================================================================
  ! function chem_is_active
  !================================================================================================
  function chem_is_active()

    logical :: chem_is_active

    chem_is_active = .true.

  end function chem_is_active

  !================================================================================================
  ! function chem_implements_cnst
  !================================================================================================
  function chem_implements_cnst(name)
    ! Purpose: return true if specified constituent is implemented by this package
    ! Author: B. Eaton

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: name   ! constituent name
    LOGICAL :: chem_implements_cnst        ! return value
    INTEGER :: M

    chem_implements_cnst = .false.
    DO M = 1, gas_pcnst
       IF (TRIM(solsym(M)) .eq. TRIM(name)) THEN
          chem_implements_cnst = .true.
          EXIT
       ENDIF
    ENDDO

  end function chem_implements_cnst

  !================================================================================================
  ! subroutine chem_init
  !================================================================================================
  subroutine chem_init(phys_state, pbuf2d)
    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize GEOS-Chem parts (state objects, mainly)
    !          (and declare history variables)
    !
    !-----------------------------------------------------------------------

    ! CAM modules
    use cam_abortutils,        only : endrun
    use chem_mods,             only : map2GC_dryDep, drySpc_ndx
    use gas_wetdep_opts,       only : gas_wetdep_method
    use hycoef,                only : ps0, hyai, hybi, hyam
    use mo_chem_utls,          only : get_spc_ndx
    use mo_ghg_chem,           only : ghg_chem_init
    use mo_mean_mass,          only : init_mean_mass
    use mo_neu_wetdep,         only : neu_wetdep_init
    use mo_setinv,             only : setinv_inti
    use Phys_Grid,             only : get_Area_All_p
    use physics_buffer,        only : physics_buffer_desc, pbuf_get_index
    use spmd_utils,            only : mpicom, masterprocid, mpi_real8, mpi_success
    use tracer_cnst,           only : tracer_cnst_init
    use tracer_srcs,           only : tracer_srcs_init
#if defined( MODAL_AERO )
    use aero_model,            only : aero_model_init
    use mo_setsox,             only : sox_inti
    use mo_drydep,             only : drydep_inti
    use modal_aero_data,       only : ntot_amode, nspec_amode
#endif

    ! GEOS-Chem interface modules in CAM
    use geoschem_diagnostics_mod, only : GC_Diagnostics_Init
    use geoschem_emissions_mod,   only : GC_Emissions_Init
    use geoschem_history_mod,     only : HistoryExports_SetServices

    ! GEOS-Chem modules
    use Chemistry_Mod,         only : Init_Chemistry
    use DiagList_Mod,          only : Init_DiagList, Print_DiagList
    use Drydep_Mod,            only : depName, Ndvzind
    use Error_Mod,             only : Init_Error
    use GC_Environment_Mod,    only : GC_Init_Grid, GC_Init_StateObj
    use GC_Environment_Mod,    only : GC_Init_Extra, GC_Allocate_All
    use GC_Grid_Mod,           only : SetGridFromCtrEdges
    use Input_Mod,             only : Read_Input_File, Validate_Directories
    use Input_Opt_Mod,         only : Set_Input_Opt
    use isorropiaII_Mod,       only : Init_IsorropiaII
    use Linear_Chem_Mod,       only : Init_Linear_Chem
    use Linoz_Mod,             only : Linoz_Read
    use PhysConstants,         only : PI, PI_180, Re
    use Pressure_Mod,          only : Accept_External_ApBp
    use State_Chm_Mod,         only : Ind_
    use State_Grid_Mod,        only : Init_State_Grid, Cleanup_State_Grid
    use TaggedDiagList_Mod,    only : Init_TaggedDiagList, Print_TaggedDiagList
    use Time_Mod,              only : Accept_External_Date_Time
    use Ucx_Mod,               only : Init_Ucx
    use Vdiff_Mod,             only : Max_PblHt_For_Vdiff 

    TYPE(physics_state),                INTENT(IN   ) :: phys_state(BEGCHUNK:ENDCHUNK)
    TYPE(physics_buffer_desc), POINTER, INTENT(INOUT) :: pbuf2d(:,:)

    ! Local variables

    !----------------------------
    ! Scalars
    !----------------------------

    ! Integers
    INTEGER                :: LCHNK(BEGCHUNK:ENDCHUNK), NCOL(BEGCHUNK:ENDCHUNK)
    INTEGER                :: IWAIT, IERR
    INTEGER                :: nX, nY, nZ
    INTEGER                :: nStrat, nTrop
    INTEGER                :: I, J, L, N, M
    INTEGER                :: RC
    INTEGER                :: nLinoz

    ! Logicals
    LOGICAL                :: prtDebug
    LOGICAL                :: Found

    ! Strings
    CHARACTER(LEN=shr_kind_cl)  :: historyConfigFile
    CHARACTER(LEN=shr_kind_cl)  :: SpcName
    CHARACTER(LEN=*), PARAMETER :: subname = 'chem_init'

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! Grid setup
    REAL(fp)               :: lonVal,  latVal
    REAL(fp)               :: dLonFix, dLatFix
    REAL(f4), ALLOCATABLE  :: lonMidArr(:,:),  latMidArr(:,:)
    REAL(f4), ALLOCATABLE  :: lonEdgeArr(:,:), latEdgeArr(:,:)
    REAL(r8), ALLOCATABLE  :: linozData(:,:,:,:)

    ! Grid with largest number of columns
    TYPE(GrdState)         :: maxGrid   ! Grid State object

    REAL(r8), ALLOCATABLE  :: Col_Area(:)
    REAL(fp), ALLOCATABLE  :: Ap_CAM_Flip(:), Bp_CAM_Flip(:)

    !REAL(r8),      POINTER :: SlsPtr(:,:,:)

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at GEOS-Chem (in chemistry/geoschem/chemistry.F90)'

    ! Initialize pointers
    SpcInfo   => NULL()

    if (debug .and. masterproc) write(iulog,'(a)') 'chem_init: initializing GEOS-Chem chemistry'

    ! LCHNK: which chunks we have on this process
    LCHNK = phys_state%LCHNK
    ! NCOL: number of atmospheric columns for each chunk
    NCOL  = phys_state%NCOL

    ! The GEOS-Chem grids on every "chunk" will all be the same size, to avoid
    ! the possibility of having differently-sized chunks
    nX = 1
    !nY = MAXVAL(NCOL)
    nY = PCOLS
    nZ = PVER

    !! Add short lived speies to buffers
    !CALL Pbuf_add_field(Trim(SLSBuffer),'global',dtype_r8,(/PCOLS,PVER,nSls/),Sls_Pbf_Idx)
    !! Initialize
    !ALLOCATE(SlsPtr(PCOLS,PVER,BEGCHUNK:ENDCHUNK), STAT=IERR)
    !IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating SlsPtr')
    !SlsPtr(:,:,:) = 0.0e+0_r8
    !DO I=1,nSls
    !   SlsPtr(:,:,:) = sls_ref_MMR(I)
    !   CALL pbuf_set_field(pbuf2d,Sls_Pbf_Idx,SlsPtr,start=(/1,1,i/),kount=(/PCOLS,PVER,1/))
    !ENDDO
    !DEALLOCATE(SlsPtr)

    ! This ensures that each process allocates everything needed for its chunks
    ALLOCATE(State_Chm(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Chm')
    ALLOCATE(State_Diag(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Diag')
    ALLOCATE(State_Grid(BEGCHUNK:ENDCHUNK), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Grid')
    ALLOCATE(State_Met(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Met')

    ! Initialize fields of the Input Options object
    CALL Set_Input_Opt( am_I_Root = MasterProc, &
                        Input_Opt = Input_Opt,  &
                        RC        = RC         )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Find maximum tropopause level, set at 40 hPa (based on GEOS-Chem 72 and 47
    ! layer grids)
    nTrop = nZ
    DO WHILE ( hyam(nZ+1-nTrop) * ps0 < 4000.0 )
       nTrop = nTrop-1
    ENDDO
    ! Find stratopause level, defined at 1 hPa
    nStrat = nZ
    DO WHILE ( hyam(nZ+1-nStrat) * ps0 < 100.0 )
       nStrat = nStrat-1
    ENDDO

    ! Initialize grid with largest number of columns
    ! This is required as State_Grid(LCHNK) can have different
    ! number of columns, but GEOS-Chem arrays are defined based
    ! on State_Grid(BEGCHUNK).
    ! To go around this, we define all of GEOS-Chem arrays with
    ! size PCOLS x PVER, which is the largest possible number of
    ! grid cells. 
    CALL Init_State_Grid( Input_Opt  = Input_Opt, &
                          State_Grid = maxGrid,   &
                          RC         = RC        )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Init_State_Grid"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    maxGrid%NX = nX
    maxGrid%NY = nY
    maxGrid%NZ = nZ

    Input_Opt%thisCPU  = myCPU
    Input_Opt%amIRoot  = MasterProc

    CALL Read_Input_File( Input_Opt  = Input_Opt, &
                          State_Grid = maxGrid,   &
                          RC         = RC        )

    ! First setup directories
    Input_Opt%Chem_Inputs_Dir      = TRIM(geoschem_cheminputs)
    Input_Opt%SpcDatabaseFile      = TRIM(speciesDB)
    Input_Opt%FAST_JX_DIR          = TRIM(geoschem_cheminputs)//'FAST_JX/v2020-02/'

    !----------------------------------------------------------
    ! CESM-specific input flags
    !----------------------------------------------------------

    ! onlineAlbedo    -> True  (use CLM albedo)
    !                 -> False (read monthly-mean albedo from HEMCO)
    Input_Opt%onlineAlbedo           = .true.

    ! applyQtend: apply tendencies of water vapor to specific humidity
    Input_Opt%applyQtend             = .False.

    ! correctConvUTLS: Apply photolytic correction for convective scavenging of soluble tracers?
    Input_Opt%correctConvUTLS        = .true.

    IF ( .NOT. Input_Opt%LSOA ) THEN
       CALL ENDRUN('CESM2-GC requires the complex SOA option to be on!')
    ENDIF

    CALL Validate_Directories( Input_Opt, RC )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Validation_Directories"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt  = Input_Opt, &
                       State_Grid = maxGrid,   &
                       RC         = RC        )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Define more variables for maxGrid
    maxGrid%MaxTropLev  = nTrop
    maxGrid%MaxStratLev = nStrat
    maxGrid%MaxChemLev = maxGrid%MaxStratLev

    DO I = BEGCHUNK, ENDCHUNK

       ! Initialize fields of the Grid State object
       CALL Init_State_Grid( Input_Opt  = Input_Opt,      &
                             State_Grid = State_Grid(I),  &
                             RC         = RC             )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered within call to "Init_State_Grid"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       State_Grid(I)%NX = nX
       State_Grid(I)%NY = NCOL(I)
       State_Grid(I)%NZ = nZ

       ! Initialize GEOS-Chem horizontal grid structure
       CALL GC_Init_Grid( Input_Opt  = Input_Opt,      &
                          State_Grid = State_Grid(I),  &
                          RC         = RC             )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Define more variables for State_Grid
       State_Grid(I)%MaxTropLev  = nTrop
       State_Grid(I)%MaxStratLev = nStrat

       ! Set maximum number of levels in the chemistry grid
       State_Grid(I)%MaxChemLev = State_Grid(I)%MaxStratLev

    ENDDO

    ! Note - this is called AFTER chem_readnl, after X, and after
    ! every constituent has had its initial conditions read. Any
    ! constituent which is not found in the CAM restart file will
    ! then have already had a call to chem_implements_cnst, and will
    ! have then had a call to chem_init_cnst to set a default VMR
    ! Call the routine GC_Allocate_All (located in module file
    ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
    ! allocatable arrays used by GEOS-Chem.
    CALL GC_Allocate_All ( Input_Opt      = Input_Opt, &
                           State_Grid     = maxGrid,   &
                           RC             = RC        )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GC_Allocate_All"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Read in data for Linoz. All CPUs allocate one array to hold the data. Only
    ! the root CPU reads in the data; then we copy it out to a temporary array,
    ! broadcast to all other CPUs, and finally duplicate the data into every
    ! copy of Input_Opt
    IF ( Input_Opt%LLinoz ) THEN
       ! Allocate array for broadcast
       nLinoz = Input_Opt%Linoz_NLevels * &
                Input_Opt%Linoz_NLat    * &
                Input_Opt%Linoz_NMonths * &
                Input_Opt%Linoz_NFields
       ALLOCATE( linozData( Input_Opt%Linoz_NLevels,     &
                            Input_Opt%Linoz_NLat,        &
                            Input_Opt%Linoz_NMonths,     &
                            Input_Opt%Linoz_NFields  ), STAT=IERR)
       IF (IERR .NE. 0) CALL ENDRUN('Failure while allocating linozData')
       linozData = 0.0e+0_r8

       IF ( MasterProc ) THEN
          ! Read data in to Input_Opt%Linoz_TParm
          CALL Linoz_Read( Input_Opt, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Linoz_Read"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
          ! Copy the data to a temporary array
          linozData = REAL(Input_Opt%LINOZ_TPARM, r8)
       ENDIF
       CALL mpi_bcast(linozData, nLinoz, mpi_real8, masterprocid, mpicom, ierr)
       IF ( ierr /= mpi_success ) then
          CALL endrun(subname//': MPI_BCAST ERROR: linozData')
       ENDIF
       IF ( .NOT. MasterProc ) THEN
          Input_Opt%LINOZ_TPARM = REAL(linozData,fp)
       ENDIF
       IF ( ALLOCATED( linozData ) ) DEALLOCATE(linozData)
    ENDIF

    ! Note: The following calculations do not setup the gridcell areas.
    !       In any case, we will need to be constantly updating this grid
    !       to compensate for the "multiple chunks per processor" element
    ALLOCATE(lonMidArr(maxGrid%nX,maxGrid%nY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonMidArr')
    ALLOCATE(lonEdgeArr(maxGrid%nX+1,maxGrid%nY+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonEdgeArr')
    ALLOCATE(latMidArr(maxGrid%nX,maxGrid%nY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating latMidArr')
    ALLOCATE(latEdgeArr(maxGrid%nX+1,maxGrid%nY+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating latEdgeArr')

    ! We could try and get the data from CAM.. but the goal is to make this GC
    ! component completely grid independent. So for now, we set to arbitrary
    ! values
    ! TODO: This needs more refinement. For now, this generates identical
    ! State_Grid for all chunks
    DO L = BEGCHUNK, ENDCHUNK
       lonMidArr = 0.0e+0_f4
       latMidArr = 0.0e+0_f4
       dLonFix   = 360.0e+0_fp / REAL(nX,fp)
       dLatFix   = 180.0e+0_fp / REAL(NCOL(L),fp)
       DO I = 1, nX
          ! Center of box, assuming dateline edge
          lonVal = -180.0e+0_fp + (REAL(I-1,fp)*dLonFix)
          DO J = 1, NCOL(L)
             ! Center of box, assuming regular cells
             latVal = -90.0e+0_fp + (REAL(J-1,fp)*dLatFix)
             lonMidArr(I,J)  = REAL((lonVal + (0.5e+0_fp * dLonFix)) * PI_180, f4)
             latMidArr(I,J)  = REAL((latVal + (0.5e+0_fp * dLatFix)) * PI_180, f4)

             ! Edges of box, assuming regular cells
             lonEdgeArr(I,J) = REAL(lonVal * PI_180, f4)
             latEdgeArr(I,J) = REAL(latVal * PI_180, f4)
          ENDDO
          ! Edges of box, assuming regular cells
          lonEdgeArr(I,NCOL(L)+1)  = REAL((lonVal + dLonFix) * PI_180, f4)
          latEdgeArr(I,NCOL(L)+1)  = REAL((latVal + dLatFix) * PI_180, f4)
       ENDDO
       DO J = 1, NCOL(L)+1
          ! Edges of box, assuming regular cells
          latVal = -90.0e+0_fp + (REAL(J-1,fp)*dLatFix)
          lonEdgeArr(nX+1,J)  = REAL((lonVal + dLonFix) * PI_180, f4)
          latEdgeArr(nX+1,J)  = REAL((latVal) * PI_180, f4)
       ENDDO

       CALL SetGridFromCtrEdges( Input_Opt  = Input_Opt,     &
                                 State_Grid = State_Grid(L), &
                                 lonCtr     = lonMidArr,     &
                                 latCtr     = latMidArr,     &
                                 lonEdge    = lonEdgeArr,    &
                                 latEdge    = latEdgeArr,    &
                                 RC         = RC            )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "SetGridFromCtrEdges"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

    ENDDO
    IF ( ALLOCATED( lonMidArr  ) ) DEALLOCATE( lonMidArr  )
    IF ( ALLOCATED( latMidArr  ) ) DEALLOCATE( latMidArr  )
    IF ( ALLOCATED( lonEdgeArr ) ) DEALLOCATE( lonEdgeArr )
    IF ( ALLOCATED( latEdgeArr ) ) DEALLOCATE( latEdgeArr )

    ! Set the times held by "time_mod"
    CALL Accept_External_Date_Time( value_NYMDb = Input_Opt%NYMDb, &
                                    value_NHMSb = Input_Opt%NHMSb, &
                                    value_NYMDe = Input_Opt%NYMDe, &
                                    value_NHMSe = Input_Opt%NHMSe, &
                                    value_NYMD  = Input_Opt%NYMDb, &
                                    value_NHMS  = Input_Opt%NHMSb, &
                                    RC          = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Accept_External_Date_Time"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Start by setting some dummy timesteps
    CALL GC_Update_Timesteps(300.0E+0_r8)

    ! Initialize error module
    CALL Init_Error( Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Error"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set a flag to denote if we should print ND70 debug output
    prtDebug            = ( Input_Opt%LPRT .and. MasterProc )

    historyConfigFile = 'HISTORY.rc'
    ! This requires geoschem_config.yml and HISTORY.rc to be in the run directory
    ! This is the current way chosen to diagnose photolysis rates!
    CALL Init_DiagList( MasterProc, historyConfigFile, Diag_List, RC )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_DiagList"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Initialize the TaggedDiag_List (list of wildcards/tags per diagnostic)
    CALL Init_TaggedDiagList( Input_Opt%amIroot, Diag_List,  &
                              TaggedDiag_List,   RC         )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_TaggedDiagList"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( prtDebug ) THEN
       CALL Print_DiagList( Input_Opt%amIRoot, Diag_List, RC )
       CALL Print_TaggedDiagList( Input_Opt%amIRoot, TaggedDiag_List, RC )
    ENDIF

    ! There are actually two copies of the history configuration, one is contained
    ! within HistoryConfig to mimic the properties of GCHP.
    !
    ! The above original implementation is similar to GC-Classic and WRF-GC,
    ! and is used by geoschem_diagnostics_mod for lookups for certain diagnostic
    ! fields for compatibility with CAM-chem outputs.
    ! (hplin, 10/31/22)
    CALL HistoryExports_SetServices(am_I_Root     = masterproc,        &
                                    config_file   = historyConfigFile, &
                                    HistoryConfig = HistoryConfig,     &
                                    RC            = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HistoryExports_SetServices"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DO I = BEGCHUNK, ENDCHUNK
       Input_Opt%amIRoot = (MasterProc .AND. (I == BEGCHUNK))

       CALL GC_Init_StateObj( Diag_List       = Diag_List,       & ! Diagnostic list obj
                              TaggedDiag_List = TaggedDiag_List, & ! TaggedDiag list obj
                              Input_Opt       = Input_Opt,       & ! Input Options
                              State_Chm       = State_Chm(I),    & ! Chemistry State
                              State_Diag      = State_Diag(I),   & ! Diagnostics State
                              State_Grid      = maxGrid,         & ! Grid State
                              State_Met       = State_Met(I),    & ! Meteorology State
                              RC              = RC              )  ! Success or failure

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "GC_Init_StateObj"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Start with v/v dry (CAM standard)
       State_Chm(I)%Spc_Units = 'v/v dry'

    ENDDO
    Input_Opt%amIRoot = MasterProc

    CALL GC_Init_Extra( Diag_List  = Diag_List,            &  ! Diagnostic list obj
    &                   Input_Opt  = Input_Opt,            &  ! Input Options
    &                   State_Chm  = State_Chm(BEGCHUNK),  &  ! Chemistry State
    &                   State_Diag = State_Diag(BEGCHUNK), &  ! Diagnostics State
    &                   State_Grid = maxGrid,              &  ! Grid State
    &                   RC         = RC                   )   ! Success or failure

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "GC_Init_Extra"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( Input_Opt%LDryD ) THEN
       !----------------------------------------------------------
       ! Get mapping between CESM dry deposited species and the
       ! indices of State_Chm%DryDepVel. This needs to be done after
       ! Init_Drydep
       ! Thibaud M. Fritz - 04 Mar 2020
       !----------------------------------------------------------
       
       ALLOCATE(map2GC_dryDep(nddvels), STAT=IERR)
       IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate map2GC_dryDep')

       DO N = 1, nddvels
          ! Initialize index to -1
          map2GC_dryDep(N) = -1

          IF ( drySpc_ndx(N) > 0 ) THEN

             ! Convert to upper case
             SpcName = to_upper(drydep_list(N))

             DO I = 1, State_Chm(BEGCHUNK)%nDryDep
                IF ( TRIM( SpcName ) == TRIM( to_upper(depName(I)) ) ) THEN
                    map2GC_dryDep(N) = nDVZind(I)
                   EXIT
                ENDIF
             ENDDO

             ! Print out debug information
             IF ( masterProc ) THEN
                IF ( N == 1 ) Write(iulog,*) " ++ GEOS-Chem Dry deposition ++ "
                IF ( map2GC_dryDep(N) > 0 ) THEN
                    Write(iulog,*) " CESM species: ", TRIM(drydep_list(N)), &
                      ' is matched with ', depName(map2GC_dryDep(N))
                ELSE
                    Write(iulog,*) " CESM species: ", TRIM(drydep_list(N)), &
                      ' has no match'
                ENDIF
             ENDIF

          ENDIF
       ENDDO
    ENDIF

#if defined( MODAL_AERO )
    ! Initialize aqueous chem
    CALL SOx_inti()

    ! Initialize aerosols
    CALL aero_model_init( pbuf2d )

    ! Initialize drydep
    CALL drydep_inti( depvel_lnd_file )
#endif

    IF ( gas_wetdep_method == 'NEU' ) THEN
       ! Initialize MOZART's wet deposition
       CALL Neu_wetdep_init()
    ENDIF

    ! Set grid-cell area
    DO N = BEGCHUNK, ENDCHUNK
       ALLOCATE(Col_Area(State_Grid(N)%nY), STAT=IERR)
       IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Col_Area')

       CALL Get_Area_All_p(N, State_Grid(N)%nY, Col_Area)

       ! Set default value (in case of chunks with fewer columns)
       State_Grid(N)%Area_M2 = 1.0e+10_fp
       DO I = 1, State_Grid(N)%nX
       DO J = 1, State_Grid(N)%nY
          State_Grid(N)%Area_M2(I,J) = REAL(Col_Area(J) * Re**2,fp)
          State_Met(N)%Area_M2(I,J)  = State_Grid(N)%Area_M2(I,J)
       ENDDO
       ENDDO

       IF ( ALLOCATED( Col_Area ) ) DEALLOCATE(Col_Area)
    ENDDO

    ! Initialize (mostly unused) diagnostic arrays
    ! WARNING: This routine likely calls on modules which are currently
    ! excluded from the GC-CESM build (eg diag03)
    ! CALL Initialize( MasterProc, Input_Opt, 2, RC )
    ! CALL Initialize( Masterproc, Input_Opt, 3, RC )

    ! Get Ap and Bp from CAM at pressure edges
    ALLOCATE(Ap_CAM_Flip(nZ+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Ap_CAM_Flip')
    ALLOCATE(Bp_CAM_Flip(nZ+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Bp_CAM_Flip')

    Ap_CAM_Flip = 0.0e+0_fp
    Bp_CAM_Flip = 0.0e+0_fp
    DO I = 1, nZ+1
       Ap_CAM_Flip(I) = hyai(nZ+2-I) * ps0 * 0.01e+0_r8
       Bp_CAM_Flip(I) = hybi(nZ+2-I)
    ENDDO

    !-----------------------------------------------------------------
    ! Pass external Ap and Bp to GEOS-Chem's Pressure_Mod
    !-----------------------------------------------------------------
    CALL Accept_External_ApBp( State_Grid = maxGrid,     &  ! Grid State
                               ApIn       = Ap_CAM_Flip, &  ! "A" term for hybrid grid
                               BpIn       = Bp_CAM_Flip, &  ! "B" term for hybrid grid
                               RC         = RC          )   ! Success or failure

    ! Print vertical coordinates
    IF ( MasterProc ) THEN
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
       WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) Ap_CAM_Flip(1:maxGrid%nZ+1)
       WRITE( 6, '(a)'   )
       WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) Bp_CAM_Flip(1:maxGrid%nZ+1)
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
    ENDIF

    ! Trapping errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Accept_External_ApBp"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( ALLOCATED( Ap_CAM_Flip ) ) DEALLOCATE( Ap_CAM_Flip )
    IF ( ALLOCATED( Bp_CAM_Flip ) ) DEALLOCATE( Bp_CAM_Flip )

    ! Once the initial met fields have been read in, we need to find
    ! the maximum PBL level for the non-local mixing algorithm.
    CALL Max_PblHt_For_Vdiff( Input_Opt  = Input_Opt,            &
                              State_Grid = State_Grid(BEGCHUNK), &
                              State_Met  = State_Met(BEGCHUNK),  &
                              RC         = RC                   )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Max_PblHt_for_Vdiff"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    IF ( Input_Opt%Its_A_FullChem_Sim .OR. &
         Input_Opt%Its_An_Aerosol_Sim ) THEN
       ! This also initializes Fast-JX
       CALL Init_Chemistry( Input_Opt  = Input_Opt,            &
                            State_Chm  = State_Chm(BEGCHUNK),  &
                            State_Diag = State_Diag(BEGCHUNK), &
                            State_Grid = State_Grid(BEGCHUNK), &
                            RC         = RC                    )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_Chemistry"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF
    ENDIF

    ! hplin 3/3/23: note, since we moved UCX module variables to
    ! individual State_Chm variables, Init_UCX has to be called
    ! for all chunks (all State_Chm) to properly initialize all
    ! variables.
    IF ( Input_Opt%LChem ) THEN
        DO I = BEGCHUNK, ENDCHUNK
           CALL Init_UCX( Input_Opt  = Input_Opt,                &
                          State_Chm  = State_Chm(I),             &
                          State_Diag = State_Diag(I),            &
                          State_Grid = State_Grid(I)           )

           ! Because not all CPUs in the communicator have the same amount of chunks,
           ! it is only guaranteed that the first chunk in all CPUs can participate in
           ! MPI_bcast of the NOXCOEFF array. So only the root CPU & root chunk will
           ! read the NOXCOEFF array from disk, then broadcast to all other CPU's first
           ! chunks, then remaining chunks can be copied locally without MPI. (hplin, 10/17/23)
           IF( I == BEGCHUNK ) THEN
              CALL mpi_bcast( State_Chm(I)%NOXCOEFF, size(State_Chm(I)%NOXCOEFF), mpi_real8, masterprocid, mpicom, ierr )
              IF ( ierr /= mpi_success ) CALL endrun('Error in mpi_bcast of NOXCOEFF in first chunk')
           ELSE
              State_CHM(I)%NOXCOEFF = State_Chm(BEGCHUNK)%NOXCOEFF
           ENDIF
        ENDDO
    ENDIF

    IF ( Input_Opt%Linear_Chem ) THEN
        CALL Init_Linear_Chem( Input_Opt  = Input_Opt,           &
                               State_Chm  = State_Chm(BEGCHUNK), &
                               State_Met  = State_Met(BEGCHUNK), &
                               State_Grid = maxGrid,             &
                               RC         = RC                  )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Init_Linear_Chem"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF ( Input_Opt%LSSalt ) THEN
       CALL INIT_ISORROPIAII( State_Grid = maxGrid )
    ENDIF

    ! Get some indices
    iH2O  = Ind_('H2O')
    iO3   = Ind_('O3')
    iCO2  = Ind_('CO2')
    iSO4  = Ind_('SO4')
    ! The following indices are needed to compute invariants
    iO    = Ind_('O')
    iH    = Ind_('H')
    iO2   = Ind_('O2')

    ! This is used to compute overhead ozone column
    SpcInfo => State_Chm(BEGCHUNK)%SpcData(iO3)%Info
    MWO3    = REAL(SpcInfo%MW_g,r8)
    ! Free pointer
    SpcInfo => NULL()

    l_H2SO4 = get_spc_ndx('H2SO4', ignore_case=.true.)
    l_SO4   = get_spc_ndx('SO4', ignore_case=.true.)

    ! Get indices for physical fields in physics buffer
    NDX_PBLH     = pbuf_get_index('pblh'     )
    NDX_FSDS     = pbuf_get_index('FSDS'     )
    NDX_CLDTOP   = pbuf_get_index('CLDTOP'   )
    NDX_CLDFRC   = pbuf_get_index('CLD'      )
    NDX_PRAIN    = pbuf_get_index('PRAIN'    )
    NDX_NEVAPR   = pbuf_get_index('NEVAPR'   )
    NDX_LSFLXPRC = pbuf_get_index('LS_FLXPRC')
    NDX_LSFLXSNW = pbuf_get_index('LS_FLXSNW')
    NDX_CMFDQR   = pbuf_get_index('RPRDTOT'  )

    ! Get cloud water indices
    CALL cnst_get_ind( 'CLDLIQ', ixCldLiq)
    CALL cnst_get_ind( 'CLDICE', ixCldIce)
    CALL cnst_get_ind( 'NUMLIQ', ixNDrop, abort=.False.  )

    CALL init_mean_mass()
    CALL setinv_inti()

    !-----------------------------------------------------------------------
    !     ... initialize tracer modules
    !-----------------------------------------------------------------------
    CALL tracer_cnst_init()
    CALL tracer_srcs_init()

    IF ( ghg_chem ) THEN
       CALL ghg_chem_init(phys_state, bndtvg, h2orates)
    ENDIF

    ! Initialize diagnostics interface
    CALL GC_Diagnostics_Init( Input_Opt = Input_Opt,           &
                              State_Chm = State_Chm(BEGCHUNK), &
                              State_Met = State_Met(BEGCHUNK) )

    ! Initialize emissions interface
    CALL GC_Emissions_Init( )

    hco_pbuf2d => pbuf2d

    ! Cleanup
    Call Cleanup_State_Grid( maxGrid, RC )

    if (masterproc) write(iulog,'(a)') 'chem_init: GEOS-Chem chemistry initialization complete'

  end subroutine chem_init

  !================================================================================================
  ! chem_timestep_init
  !================================================================================================
  subroutine chem_timestep_init(phys_state, pbuf2d)

    ! CAM modules
    use mo_flbc,           only : flbc_chk
    use mo_ghg_chem,       only : ghg_chem_timestep_init
    use physics_buffer,    only : physics_buffer_desc
    
    TYPE(physics_state), INTENT(IN):: phys_state(begchunk:endchunk)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    ! Not sure what we would realistically do here rather than in tend

    !-----------------------------------------------------------------------
    ! Set fixed lower boundary timing factors
    !-----------------------------------------------------------------------
    CALL flbc_chk

    IF ( ghg_chem ) THEN
       CALL ghg_chem_timestep_init(phys_state)
    ENDIF

  end subroutine chem_timestep_init

  !================================================================================================
  ! subroutine gc_update_timesteps
  !================================================================================================
  subroutine gc_update_timesteps(DT)

    ! GEOS-Chem modules
    use Time_Mod,       only : Set_Timesteps

    REAL(r8), INTENT(IN) :: DT
    INTEGER              :: DT_MIN
    INTEGER, SAVE        :: DT_MIN_LAST = -1

    DT_MIN = NINT(DT)

    Input_Opt%TS_CHEM = DT_MIN
    Input_Opt%TS_EMIS = DT_MIN
    Input_Opt%TS_CONV = DT_MIN
    Input_Opt%TS_DYN  = DT_MIN
    Input_Opt%TS_RAD  = DT_MIN

    ! Only bother updating the module information if there's been a change
    IF (DT_MIN .NE. DT_MIN_LAST) THEN
        CALL Set_Timesteps( Input_Opt  =  Input_Opt, &
                            CHEMISTRY  =  DT_MIN,    &
                            EMISSION   =  DT_MIN,    &
                            DYNAMICS   =  DT_MIN,    &
                            UNIT_CONV  =  DT_MIN,    &
                            CONVECTION =  DT_MIN,    &
                            DIAGNOS    =  DT_MIN,    &
                            RADIATION  =  DT_MIN    )
        DT_MIN_LAST = DT_MIN
     ENDIF

  end subroutine gc_update_timesteps

  !================================================================================================
  ! subroutine geoschem_readnl
  !================================================================================================
  subroutine geoschem_readnl(nlfile)
    ! Purpose: reads the namelist from cam/src/control/runtime_opts

    ! CAM modules
    use spmd_utils,      only : mpicom, masterprocid, mpi_character, mpi_success
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'geoschem_readnl'

    namelist /geoschem_nl/ geoschem_cheminputs

    ! Read namelist
    IF ( MasterProc ) THEN
       unitn = getunit()
       OPEN( unitn, FILE=TRIM(nlfile), STATUS='old' )
       CALL find_group_name(unitn, 'geoschem_nl', STATUS=ierr)
       IF ( ierr == 0 ) THEN
          READ(unitn, geoschem_nl, IOSTAT=ierr)
          IF ( ierr /= 0 ) THEN
             CALL ENDRUN(subname // ':: ERROR reading namelist')
          ENDIF
       ENDIF
       CLOSE(unitn)
       CALL freeunit(unitn)
    ENDIF

    ! Broadcast namelist variables
    CALL mpi_bcast(geoschem_cheminputs, LEN(geoschem_cheminputs), mpi_character, masterprocid, mpicom, ierr)
    IF ( ierr /= mpi_success ) then
       CALL endrun(subname//': MPI_BCAST ERROR: geoschem_cheminputs')
    ENDIF

  end subroutine geoschem_readnl

  !================================================================================================
  ! subroutine chem_timestep_tend
  !================================================================================================
  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dT, pbuf, fh2o )

    ! CAM modules
    use cam_history,         only : outfld, hist_fld_active
    use camsrfexch,          only : cam_in_t, cam_out_t
    use chem_mods,           only : drySpc_ndx, map2GC_dryDep
    use chem_mods,           only : nfs, indexm, gas_pcnst
    use gas_wetdep_opts,     only : gas_wetdep_method
    use mo_chem_utls,        only : get_spc_ndx
    use mo_flbc,             only : flbc_set
    use mo_ghg_chem,         only : ghg_chem_set_flbc
    use mo_mean_mass,        only : set_mean_mass
    use mo_neu_wetdep,       only : neu_wetdep_tend
    use mo_setinv,           only : setinv
    use orbit,               only : zenith                         ! For computing SZA
    use physics_buffer,      only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use physics_buffer,      only : pbuf_get_chunk, pbuf_get_index
    use perf_mod,            only : t_startf, t_stopf
    use phys_grid,           only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use phys_grid,           only : get_area_all_p, get_lat_all_p, get_lon_all_p
    use physconst,           only : MWDry, Gravit
    use rad_constituents,    only : rad_cnst_get_info
    use short_lived_species, only : get_short_lived_species_gc, set_short_lived_species_gc
    use spmd_utils,          only : masterproc
    use time_manager,        only : Get_Curr_Calday, Get_Curr_Date ! For computing SZA
    use tropopause,          only : Tropopause_findChemTrop, Tropopause_Find
    use wv_saturation,       only : QSat
#if defined( MODAL_AERO )
    use aero_model,          only : aero_model_gasaerexch ! Aqueous chemistry and aerosol growth
    use modal_aero_data,     only : ntot_amode, nspec_amode
    use modal_aero_data,     only : nspec_max, nsoa
    use modal_aero_data,     only : lmassptr_amode, numptr_amode
    use modal_aero_data,     only : lptr_so4_a_amode
    use modal_aero_data,     only : lptr2_soa_a_amode, lptr2_soa_g_amode
#endif

    ! GEOS-Chem interface modules in CAM
    use GeosChem_Emissions_Mod,   only : GC_Emissions_Calc
    use GeosChem_Diagnostics_Mod, only : GC_Diagnostics_Calc, wetdep_name, wtrate_name
    use GeosChem_History_Mod,     only : HistoryExports_SetDataPointers, CopyGCStates2Exports

    ! GEOS-Chem modules
    use Aerosol_Mod,         only : Set_AerMass_Diagnostic
    use Calc_Met_Mod,        only : Set_Dry_Surface_Pressure, AirQnt
    use Chemistry_Mod,       only : Do_Chemistry
    use CMN_FJX_MOD,         only : ZPJ
    use CMN_Size_Mod,        only : NSURFTYPE, PTop
    use Diagnostics_Mod,     only : Zero_Diagnostics_StartOfTimestep, Set_Diagnostics_EndofTimestep
    use Drydep_Mod,          only : Do_Drydep, DEPNAME, NDVZIND, Update_DryDepFreq
    use FAST_JX_MOD,         only : RXN_NO2, RXN_O3_1
    use GC_Grid_Mod,         only : SetGridFromCtr
    use HCO_Interface_GC_Mod,only : Compute_Sflx_For_Vdiff
    use Linear_Chem_Mod,     only : TrID_GC, GC_Bry_TrID, NSCHEM
    use Linear_Chem_Mod,     only : BrPtrDay, BrPtrNight, PLVEC, GMI_OH
    use Olson_Landmap_Mod,   only : Compute_Olson_Landmap
    use Modis_LAI_Mod,       only : Compute_XLAI
    use PBL_Mix_Mod,         only : Compute_PBL_Height
    use PhysConstants,       only : PI, PI_180, g0, AVO, Re, g0_100
    use Pressure_Mod,        only : Set_Floating_Pressures, Accept_External_Pedge
    use State_Chm_Mod,       only : Ind_
    use State_Diag_Mod,      only : get_TagInfo
    use Time_Mod,            only : Accept_External_Date_Time
    use Toms_Mod,            only : Compute_Overhead_O3
    use UCX_Mod,             only : Set_H2O_Trac
    use Unitconv_Mod,        only : Convert_Spc_Units
    use Wetscav_Mod,         only : Setup_Wetscav

    REAL(r8),            INTENT(IN)    :: dT          ! Time step
    TYPE(physics_state), INTENT(IN)    :: state       ! Physics State variables
    TYPE(physics_ptend), INTENT(OUT)   :: ptend       ! indivdual parameterization tendencies
    TYPE(cam_in_t),      INTENT(INOUT) :: cam_in
    TYPE(cam_out_t),     INTENT(IN)    :: cam_out
    TYPE(physics_buffer_desc), POINTER :: pbuf(:)
    REAL(r8), OPTIONAL,  INTENT(OUT)   :: fh2o(PCOLS) ! h2o flux to balance source from chemistry

    ! Initial MMR for all species
    REAL(r8) :: MMR_Beg(PCOLS,PVER,MAXVAL(map2GC(:)))
    REAL(r8) :: MMR_End(PCOLS,PVER,MAXVAL(map2GC(:)))

    ! Logical to apply tendencies to mixing ratios
    LOGICAL :: lq(pcnst)

    ! Indexing
    INTEGER :: K, N, M, P, SM, ND
    INTEGER :: I, J, L, nX, nY, nZ

    INTEGER :: LCHNK, NCOL

    REAL(r8), DIMENSION(state%NCOL) :: &
        CSZA,                          &              ! Cosine of solar zenith angle
        CSZAmid,                       &              ! Cosine of solar zenith angle at the mid timestep
        Rlats, Rlons                                  ! Chunk latitudes and longitudes (radians)

    REAL(fp)          :: O3col(state%NCOL)            ! Overhead O3 column (DU)

    REAL(r8), POINTER :: PblH(:)                      ! PBL height on each chunk [m]
    REAL(r8), POINTER :: cldTop(:)                    ! Cloud top height [?]
    REAL(r8), POINTER :: cldFrc(:,:)                  ! Cloud fraction [-]
    REAL(r8), POINTER :: Fsds(:)                      ! Downward shortwave flux at surface [W/m2]
    REAL(r8), POINTER :: PRain(:,:)                   ! Total stratiform precip. prod. (rain + snow) [kg/kg/s]
    REAL(r8), POINTER :: NEvapr(:,:)                  ! Evaporation of total precipitation (rain + snow) [kg/kg/s]
    REAL(r8), POINTER :: LsFlxPrc(:,:)                ! Large-scale downward precip. flux at interface (rain + snow) [kg/m2/s]
    REAL(r8), POINTER :: LsFlxSnw(:,:)                ! Large-scale downward precip. flux at interface (snow only) [kg/m2/s]
    REAL(r8), POINTER :: cmfdqr(:,:)                  ! Total convective precip. prod. (rain + snow) [kg/kg/s]

    REAL(r8)          :: tmpMass
    REAL(r8)          :: cldW   (state%NCOL,PVER)     ! Cloud water (kg/kg)
    REAL(r8)          :: nCldWtr(state%NCOL,PVER)     ! Droplet number concentration (#/kg)

    REAL(r8)          :: relHum (state%NCOL,PVER)     ! Relative humidity [0-1]
    REAL(r8)          :: satV   (state%NCOL,PVER)     ! Work arrays
    REAL(r8)          :: satQ   (state%NCOL,PVER)     ! Work arrays
    REAL(r8)          :: qH2O   (state%NCOL,PVER)     ! Specific humidity [kg/kg]
    REAL(r8)          :: h2ovmr (state%NCOL,PVER)     ! H2O volume mixing ratio
    REAL(r8)          :: mBar   (state%NCOL,PVER)     ! Mean wet atmospheric mass [amu]
    REAL(r8)          :: invariants(state%NCOL,PVER,nfs)
    REAL(r8)          :: reaction_rates(1,1,1)        ! Reaction rates (unused)

    ! For aerosol formation
    REAL(r8)          :: del_h2so4_gasprod(state%NCOL,PVER)

    REAL(r8)          :: vmr0(state%NCOL,PVER,gas_pcnst)
    REAL(r8)          :: vmr1(state%NCOL,PVER,gas_pcnst)
    REAL(r8)          :: vmr2(state%NCOL,PVER,gas_pcnst)

    REAL(r8)          :: wetdepflx(pcols,pcnst)       ! Wet deposition fluxes (kg/m2/s)

#if defined( MODAL_AERO )
    REAL(r8)          :: binRatio(nspec_max,ntot_amode,state%NCOL,PVER)

    REAL(r8)          :: SO4_gasRatio(state%NCOL,PVER)

    ! For SOA mapping
    REAL(r8)          :: totMass(state%NCOL,PVER)
    REAL(r8)          :: bulkMass(state%NCOL,PVER)
    REAL(r8)          :: tmpMW_g
    CHARACTER(LEN=64) :: speciesName_1, speciesName_2, speciesName_3, speciesName_4
    INTEGER           :: speciesId_1, speciesId_2, speciesId_3, speciesId_4
    INTEGER           :: iMap, nMapping, iBin, binSOA_1, binSOA_2
    INTEGER           :: K1, K2, K3, K4
    LOGICAL           :: isSOA_aerosol

#endif

    ! For emissions
    REAL(r8)          :: eflx(pcols,pver,pcnst)       ! 3-D emissions in kg/m2/s

    ! For GEOS-Chem diagnostics
    REAL(r8)          :: mmr_tend(state%NCOL,PVER,gas_pcnst)
    REAL(r8)          :: wk_out(state%NCOL)
    LOGICAL           :: Found
    
    CHARACTER(LEN=shr_kind_cl) :: tagName

    REAL(r8), PARAMETER   :: zlnd  = 0.01_r8   ! Roughness length for soil [m]
    REAL(r8), PARAMETER   :: zslnd = 0.0024_r8 ! Roughness length for snow [m]
    REAL(r8), PARAMETER   :: zsice = 0.0400_r8 ! Roughness length for sea ice [m]
    REAL(r8), PARAMETER   :: zocn  = 0.0001_r8 ! Roughness length for oean [m]

    REAL(f4)      :: lonMidArr(1,PCOLS), latMidArr(1,PCOLS)
    INTEGER       :: iMaxLoc(1)

    REAL(r8)      :: Col_Area(state%NCOL)

    ! Intermediate arrays
    INTEGER      :: Trop_Lev (PCOLS)
    REAL(r8)     :: Trop_P   (PCOLS)
    REAL(r8)     :: Trop_T   (PCOLS)
    REAL(r8)     :: Trop_Ht  (PCOLS)
    REAL(r8)     :: SnowDepth(PCOLS)
    REAL(r8)     :: cld2D    (PCOLS)
    REAL(r8)     :: Z0       (PCOLS)
    REAL(r8)     :: Sd_Ice, Sd_Lnd, Sd_Avg, Frc_Ice

    ! Estimating cloud optical depth
    REAL(r8)     :: TauCli(PCOLS,PVER)
    REAL(r8)     :: TauClw(PCOLS,PVER)
    REAL(r8), PARAMETER :: re_m   = 1.0e-05_r8 ! Cloud drop radius in m
    REAL(r8), PARAMETER :: cldMin = 1.0e-02_r8 ! Minimum cloud cover
    REAL(r8), PARAMETER :: cnst   = 1.5e+00_r8 / (re_m * 1.0e+03_r8 * g0)

    ! Calculating SZA
    REAL(r8)      :: Calday

    CHARACTER(LEN=shr_kind_cl) :: SpcName
    CHARACTER(LEN=shr_kind_cl) :: Prefix, FieldName

    LOGICAL                 :: FND
    INTEGER                 :: SpcId
    TYPE(Species),  POINTER :: SpcInfo
    TYPE(SfcMrObj), POINTER :: iSfcMrObj

    CHARACTER(LEN=63)      :: OrigUnit

    REAL(r8)               :: SlsData(PCOLS, PVER, nSls)

    INTEGER                :: currYr, currMo, currDy, currTOD
    INTEGER                :: currYMD, currHMS, currHr, currMn, currSc
    REAL(f4)               :: currUTC

    TYPE(physics_buffer_desc), POINTER :: pbuf_chnk(:) ! slice of pbuf in chnk
    REAL(r8), POINTER      :: pbuf_ik(:,:)          ! ptr to pbuf data (/pcols,pver/)
    REAL(r8), POINTER      :: pbuf_i(:)             ! ptr to pbuf data (/pcols/) horizontal only (horiz_only)
    INTEGER                :: tmpIdx                ! pbuf field id

    INTEGER                :: TIM_NDX
    INTEGER                :: IERR

    INTEGER, SAVE          :: iStep = 0
    LOGICAL, SAVE          :: FIRST = .TRUE.
    LOGICAL                :: rootChunk
    LOGICAL                :: lastChunk
    INTEGER                :: RC


    ! Initialize pointers
    SpcInfo  => NULL()
    PblH     => NULL()
    cldTop   => NULL()
    cldFrc   => NULL()
    Fsds     => NULL()
    PRain    => NULL()
    NEvapr   => NULL()
    LsFlxPrc => NULL()
    LsFlxSnw => NULL()
    cmfdqr   => NULL()
    pbuf_chnk=> NULL()
    pbuf_ik  => NULL()
    pbuf_i   => NULL()

    ! LCHNK: which chunk we have on this process
    LCHNK = state%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    NCOL  = state%NCOL

    ! Root Chunk
    rootChunk = ( MasterProc .and. (LCHNK==BEGCHUNK) )
    ! Last Chunk
    lastChunk = ( MasterProc .and. (LCHNK==ENDCHUNK) )

    ! Count the number of steps which have passed
    IF ( LCHNK .EQ. BEGCHUNK ) iStep = iStep + 1

    ! Need to update the timesteps throughout the code
    CALL GC_Update_Timesteps(dT)

    ! For safety's sake
    PTop = state%pint(1,1)*0.01e+0_fp

    ! Need to be super careful that the module arrays are updated and correctly
    ! set. NOTE: First thing - you'll need to flip all the data vertically

    nX = 1
    nY = NCOL
    nZ = PVER

    ! Update the grid lat/lons since they are module variables
    ! Assume (!) that area hasn't changed for now, as GEOS-Chem will
    ! retrieve this from State_Met which is chunked
    !CALL get_rlat_all_p( LCHNK, NCOL, Rlats )
    !CALL get_rlon_all_p( LCHNK, NCOL, Rlons )
    Rlats(1:nY) = state%Lat(1:nY)
    Rlons(1:nY) = state%Lon(1:nY)

    lonMidArr = 0.0e+0_f4
    latMidArr = 0.0e+0_f4
    DO I = 1, nX
    DO J = 1, nY
       lonMidArr(I,J) = REAL(Rlons(J), f4)
       latMidArr(I,J) = REAL(Rlats(J), f4)
    ENDDO
    ENDDO

    ! Update the grid
    CALL SetGridFromCtr( Input_Opt  = Input_Opt,         &
                         State_Grid = State_Grid(LCHNK), &
                         lonCtr     = lonMidArr,         &
                         latCtr     = latMidArr,         &
                         RC         = RC )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "SetGridFromCtr"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set area
    CALL Get_Area_All_p( LCHNK, nY, Col_Area )

    ! Field      : AREA_M2
    ! Description: Grid box surface area
    ! Unit       : -
    ! Dimensions : nX, nY
    ! Note       : Set default value (in case of chunks with fewer columns)
    State_Grid(LCHNK)%Area_M2 = -1.0e+10_fp
    State_Met(LCHNK)%Area_M2  = -1.0e+10_fp
    State_Grid(LCHNK)%Area_M2(1,:nY) = REAL(Col_Area(:nY) * Re**2,fp)
    State_Met(LCHNK)%Area_M2(1,:nY)  = State_Grid(LCHNK)%Area_M2(1,:nY)

    ! 2. Copy tracers into State_Chm
    ! Data was received in kg/kg dry
    State_Chm(LCHNK)%Spc_Units = 'kg/kg dry'
    ! Initialize ALL State_Chm species data to zero, not just tracers
    DO N = 1, State_Chm(LCHNK)%nSpecies
       State_Chm(LCHNK)%Species(N)%Conc = 0.0e+0_fp
    ENDDO

    lq(:) = .False.

    ! Map and flip gaseous species
    MMR_Beg = 0.0e+0_r8
    MMR_End = 0.0e+0_r8
    DO N = 1, pcnst
       IF ( mapCnst(N) > 0 ) lq(N) = .True.
       M = map2GC(N)
       IF ( M <= 0 ) CYCLE
       MMR_Beg(:nY,:nZ,M) = state%q(:nY,nZ:1:-1,N)
       State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ) = REAL(MMR_Beg(:nY,:nZ,M),fp)
    ENDDO

    ! We need to let CAM know that 'H2O' and 'Q' are identical
    MMR_Beg(:nY,:nZ,iH2O) = state%q(:nY,nZ:1:-1,cQ)
    State_Chm(LCHNK)%Species(iH2O)%Conc(1,:nY,:nZ) = REAL(MMR_Beg(:nY,:nZ,iH2O),fp)

    ! Retrieve previous value of species data
    SlsData(:,:,:) = 0.0e+0_r8
    CALL get_short_lived_species_gc( SlsData, LCHNK, nY, pbuf )

    IF ( iStep == 1 ) THEN
       ! Retrieve list of species with surface boundary conditions (copied from
       ! sfcvmr_mod.F90)

       ! Head of linked list
       SfcMrHead => NULL()
       iSfcMrObj => NULL()
       SpcInfo   => NULL()

       ! Loop over all species
       DO N = 1, State_Chm(BEGCHUNK)%nSpecies
          ! Species information
          SpcInfo => State_Chm(BEGCHUNK)%SpcData(N)%Info

          ! Check if field exists (note: this needs to be less than 16
          ! characters long)
          FieldName = 'HCO_'//TRIM(Prefix_SfcVMR)//TRIM(to_upper(SpcInfo%Name))
          M = pbuf_get_index(FieldName, RC)
          IF ( M > 0 ) THEN

             ! Must have positive, non-zero MW
             IF ( SpcInfo%MW_g <= 0.0_fp ) THEN
                ErrMsg = 'Cannot use surface boundary condition for species '  &
                       // TRIM(SpcInfo%Name) // ' due to invalid MW!'
                CALL ENDRUN(TRIM(ErrMsg))
             ENDIF

             ! Create new object, add to list
             ALLOCATE( iSfcMrObj, STAT=RC )
             CALL GC_CheckVar( 'sfcvmr_mod.F90:iSfcMrObj', 0, RC )
             IF ( RC /= GC_SUCCESS ) CALL ENDRUN('Failure while allocating iSfcMrObj')

             iSfcMrObj%SpcID   =  N
             iSfcMrObj%FldName =  FieldName
             iSfcMrObj%Next    => SfcMrHead
             SfcMrHead         => iSfcMrObj
             IF ( rootChunk ) THEN
                WRITE( 6, 110 ) TRIM( SpcInfo%Name ), TRIM( iSfcMrObj%FldName )
 110            FORMAT( '--> ', a, ' will use prescribed surface boundary ',   &
                        'conditions from field ', a )
             ENDIF

             ! Free the pointer
             iSfcMrObj => NULL()
          ENDIF
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------
    !        ... Reset certain GEOS-Chem diagnostics at start of timestep
    !-----------------------------------------------------------------------
    CALL Zero_Diagnostics_StartOfTimestep( Input_Opt  = Input_Opt,         &
                                           State_Diag = State_Diag(LCHNK), &
                                           RC         = RC                )

    !-----------------------------------------------------------------------
    !        ... Set atmosphere mean mass
    !-----------------------------------------------------------------------
    ! This is not meant for simulations of the ionosphere. mBar will then just
    ! be set to mwdry and does not require to pass anything besides NCOL. We
    ! can then just past a dummy array as the second argument
    !CALL Set_mean_mass( NCOL, mmr, mBar )
    CALL Set_mean_mass( NCOL, vmr0, mBar )

    ! Map and flip gaseous short-lived species
    DO N = 1, nSls
       M = map2GC_Sls(N)
       IF ( M <= 0 ) CYCLE
       State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ) = REAL(SlsData(:nY,nZ:1:-1,N),fp)
    ENDDO

#if defined( MODAL_AERO )
    ! NOTE: GEOS-Chem bulk aerosol concentrations (BCPI, BCPO, SO4, ...) are ZEROED OUT
    ! here in order to be reconstructed from the modal concentrations.
    !
    ! This means that any changes to the BULK mass will be ignored between the end
    ! of the gas_phase_chemdr and the beginning of the next!!
    !
    ! First reset State_Chm%Species to zero out MAM-inherited GEOS-Chem aerosols
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M) ! Constituent index for GEOS-Chem
          IF ( P > 0 ) K = map2GC(P) ! Index in State_Chm
          IF ( K > 0 ) State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ) = 0.0e+00_fp
       ENDDO
    ENDDO

    ! Map and vertically flip aerosols
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M) ! Constituent index for GEOS-Chem
          IF ( P <= 0 ) CYCLE
          N = lmassptr_amode(SM,M)
          K = map2GC(P) ! Index in State_Chm
          ! /!\ MAM aerosols (with cnst index N) is mapped onto GEOS-Chem
          ! species (with cnst index P, which corresponds to index K in
          ! State_Chm)

          ! Multiple MAM4 bins are mapped to same GEOS-Chem species
          State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ) = State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ) &
                                                + REAL(state%q(:nY,nZ:1:-1,N),fp) *     &
                                                   adv_mass(mapCnst(P)) /               &
                                                   adv_mass(mapCnst(N))
       ENDDO
    ENDDO

    ! Compute ratios of bin to bulk mass
    !------------------------------------------------------------------------------------------
    ! Notes for the indices used here (hplin 3/3/23):
    !
    !   K = GEOS-Chem species index in State_Chm%Species(K).
    !   P = constituent index for BULK lumped tracer in GEOS-Chem (BCPI, BCPO, DST1, DST4, SO4, SALA, SALC, OCPI, OCPO)
    !   N = constituent index for MODAL tracer in MAM4 (bc_a1, bc_a4, ...)
    !      each combination of species and mode is described by (SM, M)
    !      SM = species (i.e., bc, dst, so4, ncl, pom) in mode M
    !       M = mode number
    !   constituent indices are used in state%q(column number,level number,constituent index)
    !   chemical tracer index (NOT constituent index) is used in mo_sim_dat, e.g., adv_mass(tracer index)
    !
    ! Mapping functions:                maps from...                    ...to
    !   mapCnst(constituent index)      constituent index               chemical tracer index
    !   lmassptr_amode(SM, M)           SM, M                           constituent index (modal)
    !   map2GC(bulk constituent index)  constituent index (bulk)        GEOS-Chem species index (bulk)
    !   map2MAM4(SM, M)                 SM, M (modal)                   constituent index (bulk)
    !                                                                   (map2MAM4 is a N to 1 operation)
    ! Query functions:
    !   xname_massptr(SM, M)            SM, M                           NAME of modal aer (bc_a1, bc_a4, ...)
    !------------------------------------------------------------------------------------------
    binRatio = 0.0e+00_r8
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M)
          IF ( P <= 0 ) CYCLE
          K = map2GC(P) ! Index in State_Chm
          N = lmassptr_amode(SM,M)
          IF ( N < 0 ) CYCLE
          DO J = 1, nY
          DO L = 1, nZ
             IF ( State_Chm(LCHNK)%Species(K)%Conc(1,J,nZ+1-L) > 0.0e+00_r8 ) THEN
                binRatio(SM,M,J,L) = state%q(J,L,N)              &
                   * adv_mass(mapCnst(P)) / adv_mass(mapCnst(N)) &
                   / REAL(State_Chm(LCHNK)%Species(K)%Conc(1,J,nZ+1-L), r8)
             ENDIF
          ENDDO
          ENDDO
          ! Overwrite MMR_Beg with value from MAM
          MMR_Beg(:nY,:nZ,K) = State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ)
       ENDDO
    ENDDO

    ! Deal with secondary organic aerosols (SOAs). This mapping is using the
    ! complex SOA option in GEOS-Chem. 
    ! MAM uses five volatility bins spanning saturation concentrations from 0.01
    ! to 100 ug/m3 (logarithmically). The complex SOA option has four volatility
    ! bins that 0.1 to 100 ug/m3. We lump the lowest two bins in CESM2 to the
    ! lowest bin in GEOS-Chem.
    !
    ! The mapping goes as follows:
    ! TSOA0 + ASOAN + SOAIE + SOAGX <- soa1_a* + soa2_a*
    ! TSOA1 + ASOA1 <- soa3_a*
    ! TSOA2 + ASOA2 <- soa4_a*
    ! TSOA3 + ASOA3 <- soa5_a*
    ! TSOG0         <- SOAG0 + SOAG1
    ! TSOG1 + ASOG1 <- SOAG2
    ! TSOG2 + ASOG2 <- SOAG3
    ! TSOG3 + ASOG3 <- SOAG4

    IF ( iStep > 1 ) THEN
       ! Do not perform this mapping on initialization as we first want to
       ! overwrite soa*_a* with the GEOS-Chem SOAs.
       nMapping = 8
       DO iMap = 1, nMapping
          speciesName_1 = ''
          speciesName_2 = ''
          speciesName_3 = ''
          speciesName_4 = ''
          IF ( iMap == 1 ) THEN
             binSOA_1 = 1
             binSOA_2 = 2
             speciesName_1 = 'TSOA0'
             speciesName_2 = 'ASOAN'
             speciesName_3 = 'SOAIE'
             speciesName_4 = 'SOAGX'
          ELSEIF ( iMap == 2 ) THEN
             binSOA_1 = 3
             binSOA_2 = 3
             speciesName_1 = 'TSOA1'
             speciesName_2 = 'ASOA1'
          ELSEIF ( iMap == 3 ) THEN
             binSOA_1 = 4
             binSOA_2 = 4
             speciesName_1 = 'TSOA2'
             speciesName_2 = 'ASOA2'
          ELSEIF ( iMap == 4 ) THEN
             binSOA_1 = 5
             binSOA_2 = 5
             speciesName_1 = 'TSOA3'
             speciesName_2 = 'ASOA3'
          ELSEIF ( iMap == 5 ) THEN
             binSOA_1 = 1
             binSOA_2 = 2
             speciesName_1 = 'TSOG0'
             speciesName_2 = 'TSOG0'
          ELSEIF ( iMap == 6 ) THEN
             binSOA_1 = 3
             binSOA_2 = 3
             speciesName_1 = 'TSOG1'
             speciesName_2 = 'ASOG1'
          ELSEIF ( iMap == 7 ) THEN
             binSOA_1 = 4
             binSOA_2 = 4
             speciesName_1 = 'TSOG2'
             speciesName_2 = 'ASOG2'
          ELSEIF ( iMap == 8 ) THEN
             binSOA_1 = 5
             binSOA_2 = 5
             speciesName_1 = 'TSOG3'
             speciesName_2 = 'ASOG3'
          ELSE
             CALL ENDRUN('Unknown SOA mapping!')
          ENDIF
          isSOA_aerosol = .False.
          IF ( iMap <= 4 ) isSOA_aerosol = .True.

          ! Compute total mass from GEOS-Chem species. This sets the ratio between
          ! speciesId_1 and speciesId_2
          totMass(:nY,:nZ) = 0.0e+00_r8

          CALL cnst_get_ind( speciesName_1, speciesId_1, abort=.True. )
          CALL cnst_get_ind( speciesName_2, speciesId_2, abort=.False. )
          CALL cnst_get_ind( speciesName_3, speciesId_3, abort=.False. )
          CALL cnst_get_ind( speciesName_4, speciesId_4, abort=.False. )
          IF ( speciesId_1 > 0 ) totMass(:nY,:nZ) = totMass(:nY,:nZ) + state%q(:nY,:nZ,speciesId_1)
          IF ( speciesId_2 > 0 ) totMass(:nY,:nZ) = totMass(:nY,:nZ) + state%q(:nY,:nZ,speciesId_2)
          IF ( speciesId_3 > 0 ) totMass(:nY,:nZ) = totMass(:nY,:nZ) + state%q(:nY,:nZ,speciesId_3)
          IF ( speciesId_4 > 0 ) totMass(:nY,:nZ) = totMass(:nY,:nZ) + state%q(:nY,:nZ,speciesId_4)

          ! Compute total bulk mass from MAM
          bulkMass(:nY,:nZ) = 0.0e+00_r8
          IF ( isSOA_aerosol ) THEN
             DO iBin = binSOA_1, binSOA_2
                DO M = 1, ntot_amode
                   N = lptr2_soa_a_amode(M,iBin)
                   IF ( N <= 0 ) CYCLE
                   tmpMW_g = adv_mass(mapCnst(N))
                   bulkMass(:nY,:nZ) = bulkMass(:nY,:nZ) + state%q(:nY,:nZ,N)
                ENDDO
             ENDDO
          ELSE
             DO iBin = binSOA_1, binSOA_2
                N = lptr2_soa_g_amode(iBin)
                IF ( N <= 0 ) CYCLE
                tmpMW_g = adv_mass(mapCnst(N))
                bulkMass(:nY,:nZ) = bulkMass(:nY,:nZ) + state%q(:nY,:nZ,N)
             ENDDO
          ENDIF

          K1 = Ind_(speciesName_1)
          K2 = Ind_(speciesName_2)
          K3 = Ind_(speciesName_3)
          K4 = Ind_(speciesName_4)
          DO J = 1, nY
          DO L = 1, nZ
             ! Total SOA aerosol masses from GC are available. Partition according to the ratio given in speciesId_N to totMass summed above.
             IF ( totMass(J,L) > 0.0e+00_r8 ) THEN
                IF ( K1 > 0 ) State_Chm(LCHNK)%Species(K1)%Conc(1,J,L) = state%q(J,nZ+1-L,speciesId_1) / totMass(J,nZ+1-L) * bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_1)) / tmpMW_g
                IF ( K2 > 0 ) State_Chm(LCHNK)%Species(K2)%Conc(1,J,L) = state%q(J,nZ+1-L,speciesId_2) / totMass(J,nZ+1-L) * bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_2)) / tmpMW_g
                IF ( K3 > 0 ) State_Chm(LCHNK)%Species(K3)%Conc(1,J,L) = state%q(J,nZ+1-L,speciesId_3) / totMass(J,nZ+1-L) * bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_3)) / tmpMW_g
                IF ( K4 > 0 ) State_Chm(LCHNK)%Species(K4)%Conc(1,J,L) = state%q(J,nZ+1-L,speciesId_4) / totMass(J,nZ+1-L) * bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_4)) / tmpMW_g
             ELSE
                ! Total SOA aerosol masses from GC are unknown. In this case partition the bulkMass by 1/2 to K1 and K2.
                IF ( K1 == K2 ) THEN
                   ! ... go in same bin. This actually does not exist in the partitioning above.
                   State_Chm(LCHNK)%Species(K1)%Conc(1,J,L) = bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_1)) / tmpMW_g
                ELSE
                   State_Chm(LCHNK)%Species(K1)%Conc(1,J,L) = bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_1)) / tmpMW_g / 2.0_r8
                   State_Chm(LCHNK)%Species(K2)%Conc(1,J,L) = bulkMass(J,nZ+1-L) * adv_mass(mapCnst(speciesId_2)) / tmpMW_g / 2.0_r8
                ENDIF
             ENDIF
          ENDDO
          ENDDO
          IF ( K1 > 0 ) MMR_Beg(:nY,:nZ,K1) = State_Chm(LCHNK)%Species(K1)%Conc(1,:nY,:nZ)
          IF ( K2 > 0 ) MMR_Beg(:nY,:nZ,K2) = State_Chm(LCHNK)%Species(K2)%Conc(1,:nY,:nZ)
          IF ( K3 > 0 ) MMR_Beg(:nY,:nZ,K3) = State_Chm(LCHNK)%Species(K3)%Conc(1,:nY,:nZ)
          IF ( K4 > 0 ) MMR_Beg(:nY,:nZ,K4) = State_Chm(LCHNK)%Species(K4)%Conc(1,:nY,:nZ)
       ENDDO
    ENDIF

    ! Add gas-phase H2SO4 to GEOS-Chem SO4 (which lumps SO4 aerosol and gaseous)
    K = iSO4
    N = cH2SO4
    IF ( K > 0 .AND. N > 0 .AND. l_SO4 > 0 ) THEN
       State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ) =                                     &
                                             State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ) &
                                             + REAL(state%q(:nY,nZ:1:-1,N),fp) *         &
                                                adv_mass(l_SO4) / adv_mass(mapCnst(N))
       ! SO4_gasRatio is in mol/mol
       SO4_gasRatio(:nY,:nZ) = state%q(:nY,:nZ,N)                      &
                             * adv_mass(l_SO4) / adv_mass(mapCnst(N))  &
                             / State_Chm(LCHNK)%Species(K)%Conc(1,:nY,nZ:1:-1)
       MMR_Beg(:nY,:nZ,K)    = State_Chm(LCHNK)%Species(K)%Conc(1,:nY,:nZ)
    ENDIF
#endif

    ! Convert mass fluxes to VMR as needed for MAM4 aerosols (these operate on vmr0 - initial and vmr1 - end of timestep)
    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          ! Is a GEOS-Chem species?
          vmr0(:nY,:nZ,N) = State_Chm(LCHNK)%Species(M)%Conc(1,:nY,nZ:1:-1) * &
                            MWDry / adv_mass(N)
          ! We'll substract concentrations after chemistry later
          mmr_tend(:nY,:nZ,N) = REAL(State_Chm(LCHNK)%Species(M)%Conc(1,:nY,nZ:1:-1),r8)
       ELSEIF ( M < 0 ) THEN
          ! Is a MAM4 species? Get VMR from state%q directly.
          vmr0(:nY,:nZ,N) = state%q(:nY,:nZ,-M) * &
                            MWDry / adv_mass(N)
          mmr_tend(:nY,:nZ,N) = state%q(:nY,:nZ,-M)
       ENDIF
    ENDDO

    ! If H2O tendencies are propagated to specific humidity, then make sure
    ! that Q actually applies tendencies
    IF ( Input_Opt%applyQtend ) lq(cQ) = .True.

    IF ( ghg_chem ) lq(1) = .True.

    ! Initialize tendency array
    CALL Physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)

    ! Reset chemical tendencies
    ptend%q(:,:,:) = 0.0e+0_r8

    ! Determine current date and time
    CALL Get_Curr_Date( yr  = currYr,  &
                        mon = currMo,  &
                        day = currDy,  &
                        tod = currTOD )

    currYMD = (currYr*1000) + (currMo*100) + (currDy)
    ! Deal with subdaily
    currUTC = REAL(currTOD,f4)/3600.0e+0_f4
    currSc  = 0
    currMn  = 0
    currHr  = 0
    DO WHILE (currTOD >= 3600)
       currTOD = currTOD - 3600
       currHr  = currHr + 1
    ENDDO
    DO WHILE (currTOD >= 60)
       currTOD = currTOD - 60
       currMn  = currMn + 1
    ENDDO
    currSc  = currTOD
    currHMS = (currHr*1000) + (currMn*100) + (currSc)

    ! Calculate COS(SZA)
    Calday = Get_Curr_Calday( INT(dT/2) )
    CALL Zenith( Calday, Rlats, Rlons, CSZAmid, nY )

    Calday = Get_Curr_Calday( )
    CALL Zenith( Calday, Rlats, Rlons, CSZA, nY )

    ! Get all required data from physics buffer
    TIM_NDX = pbuf_old_tim_idx()
    CALL pbuf_get_field( pbuf, NDX_PBLH,     PblH   )
    CALL pbuf_get_field( pbuf, NDX_FSDS,     Fsds   )
    CALL pbuf_get_field( pbuf, NDX_CLDTOP,   cldTop )
    CALL pbuf_get_field( pbuf, NDX_CLDFRC,   cldFrc,   START=(/1,1,TIM_NDX/), KOUNT=(/NCOL,PVER,1/) )
    CALL pbuf_get_field( pbuf, NDX_NEVAPR,   NEvapr,   START=(/1,1/),         KOUNT=(/NCOL,PVER/))
    CALL pbuf_get_field( pbuf, NDX_PRAIN,    PRain,    START=(/1,1/),         KOUNT=(/NCOL,PVER/))
    CALL pbuf_get_field( pbuf, NDX_LSFLXPRC, LsFlxPrc, START=(/1,1/),         KOUNT=(/NCOL,PVERP/))
    CALL pbuf_get_field( pbuf, NDX_LSFLXSNW, LsFlxSnw, START=(/1,1/),         KOUNT=(/NCOL,PVERP/))
    CALL pbuf_get_field( pbuf, NDX_CMFDQR,   cmfdqr,   START=(/1,1/),         KOUNT=(/NCOL,PVER/))

    ! Get VMR and MMR of H2O
    h2ovmr = 0.0e0_fp
    qH2O   = 0.0e0_fp
    ! Note MWDry = 28.966 g/mol
    DO J = 1, nY
    DO L = 1, nZ
       qH2O(J,L) = REAL(state%q(J,L,cQ),r8)
       ! Set GEOS-Chem's H2O mixing ratio to CAM's specific humidity 'q'
       State_Chm(LCHNK)%Species(iH2O)%Conc(1,J,nZ+1-L) = qH2O(J,L)
       h2ovmr(J,L) = qH2O(J,L) * MWDry / 18.016e+0_fp
    ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !        ... Set the "invariants"
    !-----------------------------------------------------------------------
    CALL Setinv( invariants, state%t(:,:), h2ovmr, vmr0, &
                 state%pmid(:,:), nY, LCHNK, pbuf )

    ! Calculate RH (range 0-1, note still level 1 = TOA)
    relHum(:,:) = 0.0e+0_r8
    CALL QSat(state%t(:nY,:), state%pmid(:nY,:), satV, satQ, state%NCOL,PVER)
    DO J = 1, nY
    DO L = 1, nZ
       relHum(J,L) = 0.622e+0_r8 * h2ovmr(J,L) / satQ(J,L)
       relHum(J,L) = MAX( 0.0e+0_r8, MIN( 1.0e+0_r8, relHum(J,L) ) )
    ENDDO
    ENDDO

    Z0 = 0.0e+0_r8
    DO J = 1, nY
       Z0(J) = cam_in%landFrac(J) * zlnd  &
             + cam_in%iceFrac(J)  * zsice &
             + cam_in%ocnFrac(J)  * zocn
       IF (( cam_in%snowhLand(J) > 0.01_r8 ) .OR. &
           ( cam_in%snowhIce(J)  > 0.01_r8 )) THEN
          ! Land is covered in snow
          Z0(J) = zslnd
       ENDIF
    ENDDO

    ! Estimate cloud liquid water content and OD
    TauCli = 0.0e+0_r8
    TauClw = 0.0e+0_r8

    cldW(:nY,:nZ) = state%q(:nY,:nZ,ixCldLiq) + state%q(:nY,:nZ,ixCldIce)
    IF ( ixNDrop > 0 ) nCldWtr(:nY,:nZ) = state%q(:nY,:nZ,ixNDrop)

    DO J = 1, nY
    DO L = nZ, 1, -1
       ! =================================================================
       ! ===========   Compute cloud optical depth based on   ============
       ! ===========     Liao et al. JGR, 104, 23697, 1999    ============
       ! =================================================================
       !
       ! Tau = 3/2 * LWC * dZ / ( \rho_w * r_e )
       ! dZ  = - dP / ( \rho_air * g )
       ! since Pint is ascending, we can neglect the minus sign
       !
       ! Tau = 3/2 * LWC * dP / ( \rho_air * r_e * \rho_w * g )
       ! LWC / \rho_air = Q
       !
       ! Tau    = 3/2 * Q * dP / ( r_e * rho_w * g )
       ! Tau(L) = 3/2 * Q(L) * (Pint(L+1) - Pint(L)) / (re * rho_w * g )
       ! Tau(L) = Q(L) * (Pint(L+1) - Pint(L)) * Cnst
       ! Then divide by cloud fraction to get the in-cloud optical depth

       ! Unit check:                    |
       ! Q    : [kg H2O/kg air]         |
       ! Pint : [Pa]=[kg air/m/s^2]     |
       ! re   : [m]                     |   = 1.0e-5
       ! rho_w: [kg H2O/m^3]            |   = 1.0e+3
       ! g    : [m/s^2]                 |   = 9.81
       IF ( cldFrc(J,L) > cldMin ) THEN
          TauClw(J,L) = state%q(J,L,ixCldLiq)               &
                      * (state%pint(J,L+1)-state%pint(J,L)) &
                      * cnst / cldFrc(J,L)
          TauClw(J,L) = MAX(TauClw(J,L), 0.0e+00_r8)
          TauCli(J,L) = state%q(J,L,ixCldIce)               &
                      * (state%pint(J,L+1)-state%pint(J,L)) &
                      * cnst / cldFrc(J,L)
          TauCli(J,L) = MAX(TauCli(J,L), 0.0e+00_r8)
       ENDIF
    ENDDO
    ENDDO

    ! Retrieve tropopause level
    Trop_Lev = 0.0e+0_r8
    CALL Tropopause_FindChemTrop(state, Trop_Lev)
    ! Back out the pressure
    Trop_P = 1000.0e+0_r8
    DO J = 1, nY
       Trop_P(J) = state%pmid(J,Trop_Lev(J)) * 0.01e+0_r8
    ENDDO

    ! Calculate snow depth
    snowDepth = 0.0e+0_r8
    DO J = 1, nY
       Sd_Ice  = MAX(0.0e+0_r8,cam_in%snowhIce(J))
       Sd_Lnd  = MAX(0.0e+0_r8,cam_in%snowhLand(J))
       Frc_Ice = MAX(0.0e+0_r8,cam_in%iceFrac(J))
       IF (Frc_Ice > 0.0e+0_r8) THEN
          Sd_Avg = (Sd_Lnd*(1.0e+0_r8 - Frc_Ice)) + (Sd_Ice * Frc_Ice)
       ELSE
          Sd_Avg = Sd_Lnd
       ENDIF
       snowDepth(J) = Sd_Avg
    ENDDO

    ! Field      : ALBD
    ! Description: Visible surface albedo
    ! Unit       : -
    ! Dimensions : nX, nY
    State_Met(LCHNK)%ALBD      (1,:nY) = cam_in%asdir(:nY)

    ! Field      : CLDFRC
    ! Description: Column cloud fraction
    ! Unit       : -
    ! Dimensions : nX, nY
    ! Note       : Estimate column cloud fraction as the maximum cloud
    !              fraction in the column (pessimistic assumption)
    DO J = 1, nY
       State_Met(LCHNK)%CLDFRC(1,J) = MAXVAL(cldFrc(J,:))
    ENDDO

    ! Field      : EFLUX, HFLUX
    ! Description: Latent heat flux, sensible heat flux
    ! Unit       : W/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%EFLUX     (1,:nY) = cam_in%Lhf(:nY)
    State_Met(LCHNK)%HFLUX     (1,:nY) = cam_in%Shf(:nY)

    ! Field      : LandTypeFrac
    ! Description: Olson fraction per type
    ! Unit       : - (between 0 and 1)
    ! Dimensions : nX, nY, NSURFTYPE
    ! Note       : Index 1 is water
    DO N = 1, NSURFTYPE
       Write(FieldName, '(a,i2.2)') 'HCO_LANDTYPE', N-1
       tmpIdx = pbuf_get_index(FieldName, rc)
       IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
       ELSE
          CALL pbuf_get_field(pbuf, tmpIdx, pbuf_i)
          DO J = 1, nY
             State_Met(LCHNK)%LandTypeFrac(1,J,N) = pbuf_i(J)
          ENDDO
          pbuf_i   => NULL()
       ENDIF

       Write(FieldName, '(a,i2.2)') 'HCO_XLAI', N-1
       tmpIdx = pbuf_get_index(FieldName, rc)
       IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
       ELSE
          CALL pbuf_get_field(pbuf, tmpIdx, pbuf_i)
          DO J = 1, nY
             State_Met(LCHNK)%XLAI_NATIVE(1,J,N) = pbuf_i(J)
          ENDDO
          pbuf_i   => NULL()
       ENDIF
    ENDDO

    ! Field      : FRCLND, FRLAND, FROCEAN, FRSEAICE, FRLAKE, FRLANDIC
    ! Description: Olson land fraction
    !              Fraction of land
    !              Fraction of ocean
    !              Fraction of sea ice
    !              Fraction of lake
    !              Fraction of land ice
    !              Fraction of snow
    ! Unit       : -
    ! Dimensions : nX, nY
    State_Met(LCHNK)%FRCLND    (1,:ny) = 1.e+0_fp - &
                    State_Met(LCHNK)%LandTypeFrac(1,:nY,1) ! Olson Land Fraction
    State_Met(LCHNK)%FRLAND    (1,:nY) = cam_in%landFrac(:nY)
    State_Met(LCHNK)%FROCEAN   (1,:nY) = cam_in%ocnFrac(:nY) + cam_in%iceFrac(:nY)
    State_Met(LCHNK)%FRSEAICE  (1,:nY) = cam_in%iceFrac(:nY)
    State_Met(LCHNK)%FRLAKE    (1,:nY) = 0.0e+0_fp
    State_Met(LCHNK)%FRLANDIC  (1,:nY) = 0.0e+0_fp
    State_Met(LCHNK)%FRSNO     (1,:nY) = 0.0e+0_fp

    ! Field      : GWETROOT, GWETTOP
    ! Description: Root and top soil moisture
    ! Unit       : -
    ! Dimensions : nX, nY
    State_Met(LCHNK)%GWETROOT  (1,:nY) = 0.0e+0_fp
    State_Met(LCHNK)%GWETTOP   (1,:nY) = 0.0e+0_fp

    ! Field      : LAI
    ! Description: Leaf area index
    ! Unit       : m^2/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%LAI       (1,:nY) = 0.0e+0_fp

    ! Field      : PARDR, PARDF
    ! Description: Direct and diffuse photosynthetically active radiation
    ! Unit       : W/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PARDR     (1,:nY) = 0.0e+0_fp
    State_Met(LCHNK)%PARDF     (1,:nY) = 0.0e+0_fp

    ! Field      : PBLH
    ! Description: PBL height
    ! Unit       : m
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PBLH      (1,:nY) = PblH(:nY)

    ! Field      : PHIS
    ! Description: Surface geopotential height
    ! Unit       : m
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PHIS      (1,:nY) = state%Phis(:nY)

    ! Field      : PRECANV, PRECCON, PRECLSC, PRECTOT
    ! Description: Anvil precipitation @ ground
    !              Convective precipitation @ ground
    !              Large-scale precipitation @ ground
    !              Total precipitation @ ground
    ! Unit       : kg/m^2/s
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PRECANV   (1,:nY) = 0.0e+0_fp
    State_Met(LCHNK)%PRECCON   (1,:nY) = cam_out%Precc(:nY)
    State_Met(LCHNK)%PRECLSC   (1,:nY) = cam_out%Precl(:nY)
    State_Met(LCHNK)%PRECTOT   (1,:nY) = cam_out%Precc(:nY) + cam_out%Precl(:nY)

    ! Field      : TROPP
    ! Description: Tropopause pressure
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TROPP     (1,:nY) = Trop_P(:nY)

    ! Field      : PS1_WET, PS2_WET
    ! Description: Wet surface pressure at start and end of timestep
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PS1_WET   (1,:nY) = state%ps(:nY)*0.01e+0_fp
    State_Met(LCHNK)%PS2_WET   (1,:nY) = state%ps(:nY)*0.01e+0_fp

    ! Field      : SLP
    ! Description: Sea level pressure
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%SLP       (1,:nY) = state%ps(:nY)*0.01e+0_fp

    ! Field      : TS
    ! Description: Surface temperature
    ! Unit       : K
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TS        (1,:nY) = cam_in%TS(:nY)

    ! Field      : TSKIN
    ! Description: Surface skin temperature
    ! Remarks    : NOT to be confused with TS (T at 2m) (hplin, 3/20/23)
    ! Unit       : K
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TSKIN     (1,:nY) = cam_in%SST(:nY)

    ! Field      : SWGDN
    ! Description: Incident radiation @ ground
    ! Unit       : W/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%SWGDN     (1,:nY) = fsds(:nY)

    ! Field      : SNODP, SNOMAS
    ! Description: Snow depth, snow mass
    ! Unit       : m, kg/m^2
    ! Dimensions : nX, nY
    ! Note       : Conversion from m to kg/m^2
    !              \rho_{ice} = 916.7 kg/m^3
    State_Met(LCHNK)%SNODP     (1,:nY) = snowDepth(:nY)
    State_Met(LCHNK)%SNOMAS    (1,:nY) = snowDepth(:nY) * 916.7e+0_r8

    ! Field      : SUNCOS, SUNCOSmid
    ! Description: COS(solar zenith angle) at current time and midpoint
    !              of chemistry timestep
    ! Unit       : -
    ! Dimensions : nX, nY
    State_Met(LCHNK)%SUNCOS    (1,:nY) = CSZA(:nY)
    State_Met(LCHNK)%SUNCOSmid (1,:nY) = CSZAmid(:nY)

    ! Field      : UVALBEDO
    ! Description: UV surface albedo
    ! Unit       : -
    ! Dimensions : nX, nY
    IF ( Input_Opt%onlineAlbedo ) THEN
       State_Met(LCHNK)%UVALBEDO(1,:nY) = cam_in%asdir(:nY)
    ELSE
       FieldName = 'HCO_UV_ALBEDO'
       tmpIdx = pbuf_get_index(FieldName, RC)
       IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
          State_Met(LCHNK)%UVALBEDO(1,:nY) = 0.0e+0_fp
       ELSE
          pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
          CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)
          State_Met(LCHNK)%UVALBEDO(1,:nY) = pbuf_i(:nY)
          pbuf_chnk => NULL()
          pbuf_i   => NULL()
       ENDIF
    ENDIF

    ! Field      : U10M, V10M
    ! Description: E/W and N/S wind speed @ 10m height
    ! Unit       : m/s
    ! Dimensions : nX, nY
    State_Met(LCHNK)%U10M      (1,:nY) = state%U(:nY,nZ)
    State_Met(LCHNK)%V10M      (1,:nY) = state%V(:nY,nZ)

    ! Field      : USTAR
    ! Description: Friction velocity
    ! Unit       : m/s
    ! Dimensions : nX, nY
    ! Note       : We here combine the land friction velocity (fv) with
    !              the ocean friction velocity (ustar)
    DO J = 1, nY
       State_Met(LCHNK)%USTAR     (1,J) =                       &
            cam_in%fv(J)    * ( cam_in%landFrac(J))             &
          + cam_in%uStar(J) * ( 1.0e+0_fp - cam_in%landFrac(J))
    ENDDO

    ! Field      : Z0
    ! Description: Surface roughness length
    ! Unit       : m
    ! Dimensions : nX, nY
    State_Met(LCHNK)%Z0        (1,:nY) = Z0(:nY)

    ! Field      : IODIDE
    ! Description: Surface iodide concentration
    ! Unit       : nM
    ! Dimensions : nX, nY
    FieldName = 'HCO_iodide'
    tmpIdx = pbuf_get_index(FieldName, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
       State_Chm(LCHNK)%IODIDE(1,:nY) = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)
       State_Chm(LCHNK)%IODIDE(1,:nY) = pbuf_i(:nY)
       pbuf_chnk => NULL()
       pbuf_i   => NULL()
    ENDIF

    ! Field      : SALINITY
    ! Description: Ocean salinity
    ! Unit       : PSU
    ! Dimensions : nX, nY
    ! Note       : Possibly get ocean salinity from POP?
    FieldName = 'HCO_salinity'
    tmpIdx = pbuf_get_index(FieldName, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
       State_Chm(LCHNK)%SALINITY(1,:nY) = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)
       State_Chm(LCHNK)%SALINITY(1,:nY) = pbuf_i(:nY)
       pbuf_chnk => NULL()
       pbuf_i   => NULL()
    ENDIF

    ! Field      : OMOC
    ! Description: OM/OC ratio
    ! Unit       : -
    ! Dimensions : nX, nY
    IF      ( currMo == 12 .or. currMo == 1  .or. currMo == 2  ) THEN
       FieldName = 'HCO_OMOC_DJF'
    ELSE IF ( currMo == 3  .or. currMo == 4  .or. currMo == 5  ) THEN
       FieldName = 'HCO_OMOC_MAM'
    ELSE IF ( currMo == 6  .or. currMo == 7  .or. currMo == 8  ) THEN
       FieldName = 'HCO_OMOC_JJA'
    ELSE IF ( currMo == 9  .or. currMo == 10 .or. currMo == 11 ) THEN
       FieldName = 'HCO_OMOC_SON'
    ENDIF
    tmpIdx = pbuf_get_index(FieldName, rc)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       ! there is an error here and the field was not found
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
    ELSE
       CALL pbuf_get_field(pbuf, tmpIdx, pbuf_i)
       DO J = 1, nY
          State_Chm(LCHNK)%OMOC(1,J) = pbuf_i(J)
       ENDDO
       pbuf_i   => NULL()
    ENDIF

    ! Three-dimensional fields on level edges
    DO J = 1, nY
    DO L = 1, nZ+1
       ! Field      : PEDGE
       ! Description: Wet air pressure at (vertical) level edges
       ! Unit       : hPa
       ! Dimensions : nX, nY, nZ+1
       State_Met(LCHNK)%PEDGE   (1,J,L) = state%pint(J,nZ+2-L)*0.01e+0_fp

       ! Field      : CMFMC
       ! Description: Upward moist convective mass flux
       ! Unit       : kg/m^2/s
       ! Dimensions : nX, nY, nZ+1
       State_Met(LCHNK)%CMFMC   (1,J,L) = 0.0e+0_fp

       ! Field      : PFICU, PFLCU
       ! Description: Downward flux of ice/liquid precipitation (convective)
       ! Unit       : kg/m^2/s
       ! Dimensions : nX, nY, nZ+1
       State_Met(LCHNK)%PFICU   (1,J,L) = 0.0e+0_fp
       State_Met(LCHNK)%PFLCU   (1,J,L) = 0.0e+0_fp

       ! Field      : PFILSAN, PFLLSAN
       ! Description: Downward flux of ice/liquid precipitation (Large-scale & anvil)
       ! Unit       : kg/m^2/s
       ! Dimensions : nX, nY, nZ+1
       State_Met(LCHNK)%PFILSAN (1,J,L) = LsFlxSnw(J,nZ+2-L) ! kg/m2/s
       State_Met(LCHNK)%PFLLSAN (1,J,L) = MAX(0.0e+0_fp,LsFlxPrc(J,nZ+2-L) - LsFlxSnw(J,nZ+2-L)) ! kg/m2/s
    ENDDO
    ENDDO

    DO J = 1, nY
       ! Field      : CLDTOPS
       ! Description: Max cloud top height
       ! Unit       : level
       ! Dimensions : nX, nY
       State_Met(LCHNK)%CLDTOPS(1,J) = nZ + 1 - NINT(cldTop(J))
    ENDDO

    ! Three-dimensional fields on level centers
    DO J = 1, nY
    DO L = 1, nZ
       ! Field      : U, V
       ! Description: E/W and N/S component of wind
       ! Unit       : m/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%U        (1,J,L) = state%U(J,nZ+1-L)
       State_Met(LCHNK)%V        (1,J,L) = state%V(J,nZ+1-L)

       ! Field      : OMEGA
       ! Description: Updraft velocity
       ! Unit       : Pa/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%OMEGA    (1,J,L) = state%Omega(J,nZ+1-L)

       ! Field      : CLDF
       ! Description: 3-D cloud fraction
       ! Unit       : -
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%CLDF     (1,J,L) = cldFrc(J,nZ+1-L)

       ! Field      : DTRAIN
       ! Description: Detrainment flux
       ! Unit       : kg/m^2/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%DTRAIN   (1,J,L) = 0.0e+0_fp ! Used in convection

       ! Field      : DQRCU
       ! Description: Convective precipitation production rate
       ! Unit       : kg/kg dry air/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%DQRCU    (1,J,L) = 0.0e+0_fp ! Used in convection

       ! Field      : DQRLSAN
       ! Description: Large-scale precipitation production rate
       ! Unit       : kg/kg dry air/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%DQRLSAN  (1,J,L) = PRain(J,nZ+1-L) ! kg/kg/s

       ! Field      : QI, QL
       ! Description: Cloud ice/water mixing ratio
       ! Unit       : kg/kg dry air
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%QI       (1,J,L) = MAX(1.0e-10_fp, state%q(J,nZ+1-L,ixCldIce)) ! kg ice / kg dry air
       State_Met(LCHNK)%QL       (1,J,L) = MAX(1.0e-10_fp, state%q(J,nZ+1-L,ixCldLiq)) ! kg water / kg dry air

       ! Field      : RH
       ! Description: Relative humidity
       ! Unit       : %
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%RH       (1,J,L) = relHum(J,nZ+1-L)  * 100.0e+0_fp

       ! Field      : TAUCLI, TAUCLW
       ! Description: Optical depth of ice/H2O clouds
       ! Unit       : -
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%TAUCLI   (1,J,L) = TauCli(J,nZ+1-L)
       State_Met(LCHNK)%TAUCLW   (1,J,L) = TauClw(J,nZ+1-L)

       ! Field      : REEVAPCN
       ! Description: Evaporation of convective precipitation
       !              (w/r/t dry air)
       ! Unit       : kg
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%REEVAPCN (1,J,L) = 0.0e+0_fp

       ! Field      : REEVAPLS
       ! Description: Evaporation of large-scale + anvil precipitation
       !              (w/r/t dry air)
       ! Unit       : kg/kg/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%REEVAPLS (1,J,L) = NEvapr(J,nZ+1-L) ! kg/kg/s

       ! Field      : SPHU1, SPHU2
       ! Description: Specific humidity at current and next timestep
       ! Unit       : g H2O/ kg air
       ! Dimensions : nX, nY, nZ
       ! Note       : Since we are using online meteorology, we do not have
       !              access to the data at the next time step
       !              Compute tendency in g H2O/kg air/s (tmmf, 1/13/20) ?
       State_Met(LCHNK)%SPHU1    (1,J,L) = qH2O(J,nZ+1-L)    * 1.0e+3_fp    ! g/kg
       State_Met(LCHNK)%SPHU2    (1,J,L) = qH2O(J,nZ+1-L)    * 1.0e+3_fp    ! g/kg

       ! Field      : TMPU1, TMPU2
       ! Description: Temperature at current and next timestep
       ! Unit       : K
       ! Dimensions : nX, nY, nZ
       ! Note       : Since we are using online meteorology, we do not have
       !              access to the data at the next time step
       !              Compute tendency in K/s (tmmf, 1/13/20) ?
       State_Met(LCHNK)%TMPU1    (1,J,L) = state%t(J,nZ+1-L)
       State_Met(LCHNK)%TMPU2    (1,J,L) = state%t(J,nZ+1-L)
    ENDDO
    ENDDO
    ! Note: Setting DQRLSAN to zero in the top layer prevents upcoming NaNs
    ! in the GEOS-Chem wet deposition routines. Given the altitude, it should
    ! be zero anyway, this is just to prevent any numerical artifacts from
    ! creeping in.
    State_Met(LCHNK)%DQRLSAN  (1,:nY,nZ) = 0.0e+00_fp

    ! Field      : T
    ! Description: Temperature at current time
    ! Unit       : K
    ! Dimensions : nX, nY, nZ
    ! Note       : Since we are using online meteorology, we do not have
    !              access to the data at the next time step
    !              Compute tendency in K/s (tmmf, 1/13/20) ?
    State_Met(LCHNK)%T    = (State_Met(LCHNK)%TMPU1 + State_Met(LCHNK)%TMPU2)*0.5e+0_fp

    ! Field      : SPHU
    ! Description: Specific humidity at current time
    ! Unit       : g H2O/ kg air
    ! Dimensions : nX, nY, nZ
    ! Note       : Since we are using online meteorology, we do not have
    !              access to the data at the next time step
    !              Compute tendency in g H2O/kg air/s (tmmf, 1/13/20) ?
    State_Met(LCHNK)%SPHU = (State_Met(LCHNK)%SPHU1 + State_Met(LCHNK)%SPHU2)*0.5e+0_fp

    ! Field      : OPTD
    ! Description: Total in-cloud optical depth (visible band)
    ! Unit       : -
    ! Dimensions : nX, nY, nZ
    State_Met(LCHNK)%OPTD =  State_Met(LCHNK)%TAUCLI + State_Met(LCHNK)%TAUCLW

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = currYMD,            &
                                    value_NHMS     = currHMS,            &
                                    value_YEAR     = currYr,             &
                                    value_MONTH    = currMo,             &
                                    value_DAY      = currDy,             &
                                    value_DAYOFYR  = INT(FLOOR(Calday)), &
                                    value_HOUR     = currHr,             &
                                    value_MINUTE   = currMn,             &
                                    value_HELAPSED = 0.0e+0_f4,          &
                                    value_UTC      = currUTC,            &
                                    RC             = RC    )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to update time in GEOS-Chem!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Accept_External_PEdge( State_Met  = State_Met(LCHNK),  &
                                State_Grid = State_Grid(LCHNK), &
                                RC         = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to update pressure edges!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Field      : PS1_DRY, PS2_DRY
    ! Description: Dry surface pressure at current and next timestep
    ! Unit       : hPa
    ! Dimensions : nX, nY, nZ+1
    ! Note       : 1. Use the CAM PSDry fields instead of using the
    !                 GEOS-Chem calculation
    !              2. As we are using online meteorology, we do not
    !                 have access to the fields at the next time step
    !                 Compute Pa/s tendency? (tmmf, 1/13/20)
    State_Met(LCHNK)%PS1_DRY (1,:nY) = state%PSDry(:nY) * 0.01e+0_fp
    State_Met(LCHNK)%PS2_DRY (1,:nY) = state%PSDry(:nY) * 0.01e+0_fp

    ! Field      : PSC2_WET, PSC2_DRY
    ! Description: Interpolated wet and dry surface pressure at the
    !              current time
    ! Unit       : hPa
    ! Dimensions : nX, nY, nZ+1
    ! Note       : As we are using online meteorology, we do not
    !              have access to the fields at the next time step
    !              Compute Pa/s tendency? (tmmf, 1/13/20)
    State_Met(LCHNK)%PSC2_WET = State_Met(LCHNK)%PS1_WET
    State_Met(LCHNK)%PSC2_DRY = State_Met(LCHNK)%PS1_DRY

    CALL Set_Floating_Pressures( State_Grid = State_Grid(LCHNK), &
                                 State_Met  = State_Met(LCHNK),  &
                                 RC         = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to set floating pressures!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set quantities of interest but do not change VMRs
    ! This function updates:
    !  ====================================================================
    !  (1)  PEDGE     : Moist air pressure at grid box bottom      [hPa]
    !  (2)  PEDGE_DRY : Dry air partial pressure at box bottom     [hPa]
    !  (3)  PMID      : Moist air pressure at grid box centroid    [hPa]
    !  (4)  PMID_DRY  : Dry air partial pressure at box centroid   [hPa]
    !  (5)  PMEAN     : Altitude-weighted mean moist air pressure  [hPa]
    !  (6)  PMEAN_DRY : Alt-weighted mean dry air partial pressure [hPa]
    !  (7)  DELP      : Delta-P extent of grid box                 [hPa]
    !                   (Same for both moist and dry air since we
    !                   assume constant water vapor pressure
    !                   across box)
    !  (8)  AIRDEN    : Mean grid box dry air density            [kg/m^3]
    !                   (defined as total dry air mass/box vol)
    !  (9)  AIRNUMDEN : Mean grid box dry air number density  [molec/m^3]
    !  (10) MAIRDEN   : Mean grid box moist air density          [kg/m^3]
    !                   (defined as total moist air mass/box vol)
    !  (11) AD        : Total dry air mass in grid box               [kg]
    !  (12) ADMOIST   : Total moist air mass in grid box             [kg]
    !  (13) BXHEIGHT  : Vertical height of grid box                   [m]
    !  (14) AIRVOL    : Volume of grid box                          [m^3]
    !  (15) MOISTMW   : Molecular weight of moist air in box      [g/mol]
    !  (16) IsLand    : Logical for grid cells over land              [-]
    !  (17) IsWater   : Logical for grid cells over water             [-]
    !  (18) IsIce     : Logical for grid cells over ice               [-]
    !  (19) IsSnow    : Logical for grid cells over snow              [-]
    !  (20) InTroposph: Logical for tropospheric grid cells           [-]
    !  (21) InStratMes: Logical for non-tropospheric grid cells       [-]
    !  (22) InStratosp: Logical for stratospheric grid cells          [-]
    !  (23) InChemGrid: Logical for chemistry grid cells              [-]
    !  (24) LocalSolar: Local solar time                              [-]
    !  (25) IsLocalNoo: Logical for local noon                        [-]
    !  (26) TropLev   : Maximum tropopause level                      [-]
    !  (27) TropHt    : Maximum tropopause height                    [km]
    !  ====================================================================
    CALL AirQnt( Input_Opt           = Input_Opt,         &
                 State_Chm           = State_Chm(LCHNK),  &
                 State_Grid          = State_Grid(LCHNK), &
                 State_Met           = State_Met(LCHNK),  &
                 RC                  = RC,                &
                 Update_Mixing_Ratio = .False. )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to calculate air properties!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! SDE 05/28/13: Set H2O to State_Chm tracer if relevant and,
    ! if LSETH2O=F and LACTIVEH2O=T, update specific humidity
    ! in the stratosphere
    !
    ! NOTE: Specific humidity may change in SET_H2O_TRAC and
    ! therefore this routine may call AIRQNT again to update
    ! air quantities and tracer concentrations (ewl, 10/28/15)
    IF ( Input_Opt%Its_A_Fullchem_Sim .and. iH2O > 0 ) THEN
       CALL Set_H2O_Trac( SETSTRAT   = Input_Opt%LSETH2O,          &
                          Input_Opt  = Input_Opt,                  &
                          State_Chm  = State_Chm(LCHNK),           &
                          State_Grid = State_Grid(LCHNK),          &
                          State_Met  = State_Met(LCHNK),           &
                          RC         = RC                         )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Set_H2O_Trac" #1!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Only force strat once if using UCX
       IF (Input_Opt%LSETH2O) Input_Opt%LSETH2O = .FALSE.
    ENDIF

    ! Do this after AirQnt, such that we overwrite GEOS-Chem isLand, isWater and
    ! isIce, which are based on albedo. Rather, we use CLM landFranc, ocnFrac
    ! and iceFrac. We also compute isSnow
    DO J = 1, nY
       iMaxLoc = MAXLOC( (/ State_Met(LCHNK)%FRLAND(1,J)   + &
                            State_Met(LCHNK)%FRLANDIC(1,J) + &
                            State_Met(LCHNK)%FRLAKE(1,J),    &
                            State_Met(LCHNK)%FRSEAICE(1,J),  &
                            State_Met(LCHNK)%FROCEAN(1,J)  - &
                            State_Met(LCHNK)%FRSEAICE(1,J) /) )
       IF ( iMaxLoc(1) == 3 ) iMaxLoc(1) = 0
       ! reset ocean to 0

       IF ( iMaxLoc(1) == 0 ) THEN
          State_Met(LCHNK)%isLand(1,J)  = .False.
          State_Met(LCHNK)%isWater(1,J) = .True.
          State_Met(LCHNK)%isIce(1,J)   = .False.
       ELSEIF ( iMaxLoc(1) == 1 ) THEN
          State_Met(LCHNK)%isLand(1,J)  = .True.
          State_Met(LCHNK)%isWater(1,J) = .False.
          State_Met(LCHNK)%isIce(1,J)   = .False.
       ELSEIF ( iMaxLoc(1) == 2 ) THEN
          State_Met(LCHNK)%isLand(1,J)  = .False.
          State_Met(LCHNK)%isWater(1,J) = .False.
          State_Met(LCHNK)%isIce(1,J)   = .True.
       ELSE
          Write(iulog,*) " iMaxLoc gets value: ", iMaxLoc
          ErrMsg = 'Failed to figure out land/water'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       State_Met(LCHNK)%isSnow(1,J) = &
                               ( State_Met(LCHNK)%FRSEAICE(1,J) > 0.0e+0_fp &
                                 .or. State_Met(LCHNK)%SNODP(1,J) > 0.01 )

    ENDDO

    ! Do this after AirQnt in order to use AIRDEN and BXHEIGHT
    DO J = 1, nY
       O3col(J) = 0.0e+0_fp
       DO L = 1, nZ
          O3col(J) = O3col(J) &
                   + State_Chm(LCHNK)%Species(iO3)%Conc(1,J,L) &
                      * State_Met(LCHNK)%AIRDEN(1,J,L)   &
                      * State_Met(LCHNK)%BXHEIGHT(1,J,L)
       ENDDO
       O3col(J) = O3col(J) * ( AVO / MWO3 ) / 1e+1_fp / 2.69e+16_fp
    ENDDO

    ! Field      : TO3
    ! Description: Total overhead ozone column
    ! Unit       : DU
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TO3       (1,:nY) = O3col(:nY)

    IF ( Input_Opt%Linear_Chem .AND. &
         State_Grid(LCHNK)%MaxChemLev /= State_Grid(LCHNK)%nZ ) THEN
       IF ( iStep == 1 ) THEN
          ALLOCATE( BrPtrDay  ( 6 ), STAT=IERR )
          IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating BrPtrDay')
          ALLOCATE( BrPtrNight( 6 ), STAT=IERR )
          IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating BrPtrNight')
          DO N = 1, 6
             ! Skip if species is not defined
             IF ( GC_Bry_TrID(N) <= 0 ) CYCLE

             ! Get Bry name
             SpcName = State_Chm(LCHNK)%SpcData(GC_Bry_TrID(N))%Info%Name

             ! Construct field name using Bry name
             PREFIX = 'GEOSCCM_'//TRIM(SpcName)

             ALLOCATE( BrPtrDay(N)%MR(1,PCOLS,nZ), STAT=IERR )
             IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating BrPtrDay%MR')
             ALLOCATE( BrPtrNight(N)%MR(1,PCOLS,nZ), STAT=IERR )
             IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating BrPtrNight%MR')

             ! Get pointer to this field. These are the mixing ratios (pptv).

             ! Day
             FieldName = TRIM(PREFIX) // '_DAY'
             tmpIdx = pbuf_get_index(FieldName, RC)
             IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
                IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
                BrPtrDay(N)%MR(1,:nY,nZ:1:-1) = 0.0e+0_f4
             ELSE
                pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
                CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
                BrPtrDay(N)%MR(1,:nY,nZ:1:-1) = REAL(pbuf_ik(:nY,:nZ), f4)
                pbuf_chnk => NULL()
                pbuf_ik   => NULL()
             ENDIF

             ! Night
             FieldName = TRIM(PREFIX) // '_NIGHT'
             tmpIdx = pbuf_get_index(FieldName, RC)
             IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
                IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
                BrPtrDay(N)%MR(1,:nY,nZ:1:-1) = 0.0e+0_f4
             ELSE
                pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
                CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
                BrPtrNight(N)%MR(1,:nY,nZ:1:-1) = REAL(pbuf_ik(:nY,:nZ), f4)
                pbuf_chnk => NULL()
                pbuf_ik   => NULL()
             ENDIF

          ENDDO

          DO N = 1,NSCHEM

             ! Get GEOS-Chem species index
             M = TrID_GC(N)

             ! Skip if species is not defined
             IF ( M <= 0 ) CYCLE

             ! Get species name
             SpcName = State_Chm(LCHNK)%SpcData(M)%Info%Name

             ! ---------------------------------------------------------------
             ! Get pointers to fields
             ! ---------------------------------------------------------------

             ! Production rates [v/v/s]
             FieldName = 'GMI_PROD_'//TRIM(SpcName)

             ALLOCATE( PLVEC(N)%PROD(1,PCOLS,nZ), STAT=IERR )
             IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating PLVEC%PROD')
             ALLOCATE( PLVEC(N)%LOSS(1,PCOLS,nZ), STAT=IERR )
             IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating PLVEC%PROD')

             ! Get pointer from HEMCO
             tmpIdx = pbuf_get_index(FieldName, RC)
             IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
                IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
                PLVEC(N)%PROD(1,:nY,nZ:1:-1) = 0.0e+0_f4
                FND = .False.
             ELSE
                pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
                CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
                PLVEC(N)%PROD(1,:nY,nZ:1:-1) = REAL(pbuf_ik(:nY,:nZ),f4)
                FND = .True.
                pbuf_chnk => NULL()
                pbuf_ik   => NULL()
             ENDIF

             ! Warning message
             IF ( .NOT. FND .AND. Input_Opt%amIRoot ) THEN
                ErrMsg = 'Cannot find archived production rates for '       // &
                          TRIM(SpcName) // ' - will use value of 0.0. '        // &
                         'To use archived rates, add the following field '      // &
                         'to the HEMCO configuration file: '// TRIM( FieldName )
                CALL GC_Warning( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Loss frequency [s-1]
             FieldName = 'GMI_LOSS_'//TRIM(SpcName)

             ! Get pointer from HEMCO
             tmpIdx = pbuf_get_index(FieldName, RC)
             IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
                IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
                PLVEC(N)%LOSS(1,:nY,nZ:1:-1) = 0.0e+0_f4
                FND = .False.
             ELSE
                pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
                CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
                PLVEC(N)%LOSS(1,:nY,nZ:1:-1) = REAL(pbuf_ik(:nY,:nZ), f4)
                FND = .True.
                pbuf_chnk => NULL()
                pbuf_ik   => NULL()
             ENDIF

             ! Warning message
             IF ( .NOT. FND .AND. Input_Opt%amIRoot ) THEN
                ErrMsg= 'Cannot find archived loss frequencies for '        // &
                        TRIM(SpcName) // ' - will use value of 0.0. '          // &
                        'To use archived rates, add the following field '       // &
                        'to the HEMCO configuration file: '//TRIM(FieldName)
                CALL GC_Warning( ErrMsg, RC, ThisLoc )
             ENDIF

          ENDDO !N

          ! Get pointer to GMI_OH

          ALLOCATE( GMI_OH(1,PCOLS,nZ), STAT=IERR )
          IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating GMI_OH')

          tmpIdx = pbuf_get_index(FieldName, RC)
          IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
             IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
             GMI_OH(1,:nY,nZ:1:-1) = 0.0e+0_f4
          ELSE
             pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
             CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
             GMI_OH(1,:nY,nZ:1:-1) = REAL(pbuf_ik(:nY,:nZ), f4)
             pbuf_chnk => NULL()
             pbuf_ik   => NULL()
          ENDIF
       ENDIF

    ENDIF

    ! This is not necessary as we prescribe CH4 surface mixing ratios
    ! through CAM.
    !! Prescribe methane surface concentrations throughout PBL
    !IF ( ITS_A_FULLCHEM_SIM .and. id_CH4 > 0 ) THEN
    !
    !   ! Set CH4 concentrations
    !   CALL SET_CH4( Input_Opt  = Input_Opt,         &
    !                 State_Chm  = State_Chm(LCHNK),  &
    !                 State_Diag = State_Diag(LCHNK), &
    !                 State_Grid = State_Grid(LCHNK), &
    !                 State_Met  = State_Met(LCHNK),  &
    !                 RC         = RC                )
    !
    !   ! Trap potential errors
    !   IF ( RC /= GC_SUCCESS ) THEN
    !      ErrMsg = 'Error encountered in call to "SET_CH4"!'
    !      CALL Error_Stop( ErrMsg, ThisLoc )
    !   ENDIF
    !ENDIF

    ! Eventually initialize/reset wetdep
    IF ( Input_Opt%LConv .OR. Input_Opt%LChem .OR. Input_Opt%LWetD ) THEN
        CALL Setup_WetScav( Input_Opt  = Input_Opt,         &
                            State_Chm  = State_Chm(LCHNK),  &
                            State_Grid = State_Grid(LCHNK), &
                            State_Met  = State_Met(LCHNK),  &
                            RC         = RC                )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Setup_WetScav"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    !==============================================================
    !     ***** C O M P U T E   P B L   H E I G H T  etc. *****
    !==============================================================
    ! Move this call from the PBL mixing routines because the PBL
    ! height is used by drydep and some of the emissions routines.
    ! (ckeller, 3/5/15)
    ! This function updates:
    !  ====================================================================
    !  (1)  InPbl      : Logical indicating if we are in the PBL    [-]
    !  (2)  PBL_TOP_L  : Number of layers in the PBL                [-]
    !  (3)  PBL_TOP_hPa: Pressure at the top of the PBL             [hPa]
    !  (4)  PBL_TOP_m  : PBL height                                 [m]
    !  (5)  PBL_THICK  : PBL thickness                              [hPa]
    !  (6)  F_OF_PBL   : Fraction of grid box within the PBL        [-]
    !  (7)  F_UNDER_PBLTOP: Fraction of grid box underneath the PBL top [-]
    !  (8)  PBL_MAX_L  : Model level where PBL top occurs           [-]
    !  ====================================================================
    CALL Compute_PBL_Height( Input_Opt  = Input_Opt,         &
                             State_Grid = State_Grid(LCHNK), &
                             State_Met  = State_Met(LCHNK),  &
                             RC         = RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Compute_PBL_Height"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    !--------------------------------------------------------------
    ! Test for emission timestep
    ! Now always do emissions here, even for full-mixing
    ! (ckeller, 3/5/15)
    !--------------------------------------------------------------
    !==================================================================
    !         ***** D R Y   D E P O S I T I O N *****
    !==================================================================
    !==================================================================
    ! Compute dry deposition velocities
    !
    ! CLM computes dry deposition velocities but only for gas-phase
    ! species and only over land. We therefore need to both pass the
    ! the CLM dry deposition velocities as well as compute them using
    ! the GEOS-Chem dry deposition module. If using the CLM velocities,
    ! then scale them with the ocean fraction; otherwise use GEOS-Chem
    ! computed velocities.
    !
    ! drydep_method must be set to DD_XLND.
    !
    !==================================================================
    !
    ! State_Chm expects dry deposition velocities in m/s, whereas
    ! CLM returns land deposition velocities in cm/s!
    !
    ! For now, dry deposition velocities are only computed for gases
    ! (which is what CLM deals with). Dry deposition for aerosols is
    ! work in progress.
    !
    ! Thibaud M. Fritz - 27 Feb 2020
    !==================================================================

    IF ( Input_Opt%LDryD ) THEN
       ! Compute the Olson landmap fields of State_Met
       ! (e.g. State_Met%IREG, State_Met%ILAND, etc.)
       CALL Compute_Olson_Landmap( Input_Opt  = Input_Opt,         &
                                   State_Grid = State_Grid(LCHNK), &
                                   State_Met  = State_Met(LCHNK),  &
                                   RC         = RC                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Olson_Landmap"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Compute State_Met%XLAI (for drydep) and State_Met%MODISLAI,
       ! which is the average LAI per grid box (for soil NOx emissions)
       CALL Compute_Xlai( Input_Opt  = Input_Opt,         &
                          State_Grid = State_Grid(LCHNK), &
                          State_Met  = State_Met(LCHNK),  &
                          RC         = RC                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Xlai"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Compute drydep velocities and update State_Chm%DryDepVel
       CALL Do_Drydep( Input_Opt  = Input_Opt,         &
                       State_Chm  = State_Chm(LCHNK),  &
                       State_Diag = State_Diag(LCHNK), &
                       State_Grid = State_Grid(LCHNK), &
                       State_Met  = State_Met(LCHNK),  &
                       RC         = RC                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Drydep"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       DO N = 1, nddvels

          !! Print debug
          !IF ( rootChunk ) THEN
          !    IF ( N == 1 ) THEN
          !    Write(iulog,*) "Number of GC dry deposition species = ", &
          !        SIZE(State_Chm(LCHNK)%DryDepVel(:,:,:),3)
          !    Write(iulog,*) "Number of CESM dry deposition species = ", &
          !        nddvels
          !    ENDIF
          !    Write(iulog,*) "N          = ", N
          !    Write(iulog,*) "drySpc_ndx = ", drySpc_ndx(N)
          !    Write(iulog,*) "GC index   = ", map2GC_dryDep(N)
          !    IF ( map2GC_dryDep(N) > 0 ) THEN
          !        Write(iulog,*) "GC name    = ", TRIM(DEPNAME(map2GC_dryDep(N)))
          !    ENDIF
          !    Write(iulog,*) "dry Species= ", TRIM(drydep_list(N))
          !    IF ( drySpc_ndx(N) > 0 ) THEN
          !        Write(iulog,*) "tracerName = ", TRIM(tracerNames(drySpc_ndx(N)))
          !    ENDIF
          !    Write(iulog,*) "CLM-depVel = ", &
          !  MAXVAL(cam_in%depvel(:nY,N)) * 1.0e-02_fp, " [m/s]"
          !    IF ( map2GC_dryDep(N) > 0 ) THEN
          !        Write(iulog,*) "GC-depVel  = ", &
          !  MAXVAL(State_Chm(LCHNK)%DryDepVel(1,:nY,map2GC_dryDep(N))), " [m/s]"
          !    ENDIF
          !ENDIF

          IF ( map2GC_dryDep(N) > 0 ) THEN
             ! State_Chm%DryDepVel is in m/s
             State_Chm(LCHNK)%DryDepVel(1,:nY,map2GC_dryDep(N)) = &
                ! This first bit corresponds to the dry deposition
                ! velocities over land as computed from CLM and
                ! converted to m/s. This is scaled by the fraction
                ! of land.
                  cam_in%depVel(:nY,N) * 1.0e-02_fp &
                   * MAX(0._fp, 1.0_fp - State_Met(LCHNK)%FROCEAN(1,:nY)) &
                ! This second bit corresponds to the dry deposition
                ! velocities over ocean and sea ice as computed from
                ! GEOS-Chem. This is scaled by the fraction of ocean
                ! and sea ice.
                + State_Chm(LCHNK)%DryDepVel(1,:nY,map2GC_dryDep(N)) &
                  * State_Met(LCHNK)%FROCEAN(1,:nY)
          ENDIF
       ENDDO

       CALL Update_DryDepFreq( Input_Opt  = Input_Opt,         &
                               State_Chm  = State_Chm(LCHNK),  &
                               State_Diag = State_Diag(LCHNK), &
                               State_Grid = State_Grid(LCHNK), &
                               State_Met  = State_Met(LCHNK),  &
                               RC         = RC                )

    ENDIF

    !===========================================================
    !      ***** M I X E D   L A Y E R   M I X I N G *****
    !===========================================================

    ! Updates from Bob Yantosca, 06/2020
    ! Compute the surface flux for the non-local mixing,
    ! (which means getting emissions & drydep from HEMCO)
    ! and store it in State_Chm%Surface_Flux
    ! 
    ! For CESM-GC, Surface_Flux will be equal to the opposite of the
    ! dry deposition flux since emissions are loaded externally
    ! ( SurfaceFlux = eflx - dflx = - dflx )
    IF ( Input_Opt%LTURB .and. Input_Opt%LNLPBL ) THEN
       CALL Compute_Sflx_For_Vdiff( Input_Opt  = Input_Opt,         &
                                    State_Chm  = State_Chm(LCHNK),  &
                                    State_Diag = State_Diag(LCHNK), &
                                    State_Grid = State_Grid(LCHNK), &
                                    State_Met  = State_Met(LCHNK),  &
                                    RC         = RC                )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Sflx_for_Vdiff"!'
          CALL Error_Stop( errMsg, thisLoc )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Get emissions from HEMCO + Lightning + Fire
    ! Add surface emissions to cam_in
    !-----------------------------------------------------------------------

    CALL GC_Emissions_Calc( state      = state,            &
                            hco_pbuf2d = hco_pbuf2d,       &
                            State_Met  = State_Met(LCHNK), &
                            cam_in     = cam_in,           &
                            eflx       = eflx,             &
                            iStep      = iStep            )

    !-----------------------------------------------------------------------
    ! Add dry deposition flux 
    ! (stored as SurfaceFlux = -dflx)
    !-----------------------------------------------------------------------

    IF ( Input_Opt%LDryD ) THEN
       DO ND = 1, State_Chm(BEGCHUNK)%nDryDep
          ! Get the species ID from the drydep ID
          N = State_Chm(BEGCHUNK)%Map_DryDep(ND)
          IF ( N <= 0 ) CYCLE
   
          M = map2GCinv(N)
          IF ( M <= 0 ) CYCLE
   
          cam_in%cflx(1:nY,M) = cam_in%cflx(1:nY,M) &
                              + State_Chm(LCHNK)%SurfaceFlux(1,1:nY,N)
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------
    ! Add non-surface emissions
    !-----------------------------------------------------------------------

    ! Use units of kg/m2 as State_Chm%Species to add emissions fluxes
    CALL Convert_Spc_Units( Input_Opt  = Input_Opt,         &
                            State_Chm  = State_Chm(LCHNK),  &
                            State_Grid = State_Grid(LCHNK), &
                            State_Met  = State_Met(LCHNK),  &
                            OutUnit    = 'kg/m2',           &
                            RC         = RC,                &
                            OrigUnit   = OrigUnit          )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DO N = 1, pcnst
       M = map2GC(N)
       IF ( M > 0 ) THEN
          ! Add to GEOS-Chem species
          State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ) = State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ) &
                                                + eflx(:nY,nZ:1:-1,N) * dT
       ELSEIF ( M < 0 ) THEN
          ! Add to constituent (mostly for MAM4 aerosols)
          ! Convert from kg/m2/s to kg/kg/s
          ptend%q(:nY,nZ:1:-1,N) = ptend%q(:nY,nZ:1:-1,N) &
                                 + eflx(:nY,nZ:1:-1,N)    &
                                   / ( g0_100 * State_Met(LCHNK)%DELP_DRY(1,:nY,:nZ) )
       ENDIF
    ENDDO

    ! Convert back to original unit
    CALL Convert_Spc_Units( Input_Opt  = Input_Opt,         &
                            State_Chm  = State_Chm(LCHNK),  &
                            State_Grid = State_Grid(LCHNK), &
                            State_Met  = State_Met(LCHNK),  &
                            OutUnit    = OrigUnit,          &
                            RC         = RC                )

    ! Convert State_Chm%Species back to original units
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    !==============================================================
    !               ***** C H E M I S T R Y *****
    !==============================================================

    call t_startf( 'chemdr' )

    ! Get the overhead column O3 for use with FAST-J
    IF ( Input_Opt%Its_A_FullChem_Sim .OR. &
         Input_Opt%Its_An_Aerosol_Sim ) THEN

       IF ( Input_Opt%LChem ) THEN
          CALL Compute_Overhead_O3( Input_Opt       = Input_Opt,                 &
                                    State_Grid      = State_Grid(LCHNK),         &
                                    State_Chm       = State_Chm(LCHNK),          &
                                    DAY             = currDy,                    &
                                    USE_O3_FROM_MET = Input_Opt%Use_O3_From_Met, &
                                    TO3             = State_Met(LCHNK)%TO3,      &
                                    RC              = RC                        )

          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Compute_Overhead_O3"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
       ENDIF
    ENDIF

    IF ( Input_Opt%Its_A_Fullchem_Sim .and. iH2O > 0 ) THEN
       CALL Set_H2O_Trac( SETSTRAT   = .False.               , &
                          Input_Opt  = Input_Opt,              &
                          State_Chm  = State_Chm(LCHNK),       &
                          State_Grid = State_Grid(LCHNK),      &
                          State_Met  = State_Met(LCHNK),       &
                          RC         = RC                     )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Set_H2O_Trac" #2!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF
    ENDIF

    ! Here, we apply surface mixing ratios for long-lived species
    ! (copied from sfcvmr_mod.F90)
    ! Loop over all objects
    iSfcMrObj => SfcMrHead
    DO WHILE( ASSOCIATED( iSfcMrObj ) )

       ! Get concentration for this species
       tmpIdx = pbuf_get_index(TRIM(iSfcMrObj%FldName), RC)
       IF ( tmpIdx < 0 .OR. (iStep == 1) ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(iSfcMrObj%FldName)
       ELSE
          CALL pbuf_get_field(pbuf, tmpIdx, pbuf_i)

          ! Set mixing ratio in PBL
          SpcInfo => State_Chm(LCHNK)%SpcData(iSfcMrObj%SpcID)%Info
          N = SpcInfo%ModelID
          IF ( N > 0 ) THEN
             DO L = 1, nZ
             DO J = 1, nY
                IF ( State_Met(LCHNK)%F_UNDER_PBLTOP(1,J,L) > 0.0_fp ) THEN
                   State_Chm(LCHNK)%Species(N)%Conc(1,J,L) =     &
                       ( pbuf_i(J) * 1.0e-9_fp       )     &
                     / ( MWDry      / SpcInfo%MW_g   )
                ENDIF  ! end selection of PBL boxes
             ENDDO
             ENDDO
          ENDIF
       ENDIF

       ! Point to next element in list
       iSfcMrObj => iSfcMrObj%Next
    ENDDO

    ! Reset photolysis rates
    ZPJ = 0.0e+0_r8

    ! Perform chemistry
    CALL Do_Chemistry( Input_Opt  = Input_Opt,         &
                       State_Chm  = State_Chm(LCHNK),  &
                       State_Diag = State_Diag(LCHNK), &
                       State_Grid = State_Grid(LCHNK), &
                       State_Met  = State_Met(LCHNK),  &
                       RC         = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Do_Chemistry"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! GEOS-Chem considers CO2 as a dead species and resets its concentration
    ! internally. Right after the call to `Do_Chemistry`, State_Chm%Species(iCO2)
    ! corresponds to the chemically-produced CO2. The real CO2 concentration
    ! is thus the concentration before chemistry + the chemically-produced CO2.
    State_Chm(LCHNK)%Species(iCO2)%Conc(1,:nY,:nZ) = State_Chm(LCHNK)%Species(iCO2)%Conc(1,:nY,:nZ) &
                                             + MMR_Beg(:nY,:nZ,iCO2)

    ! Make sure State_Chm(LCHNK) is back in kg/kg dry!
    IF ( TRIM(State_Chm(LCHNK)%Spc_Units) /= 'kg/kg dry' ) THEN
       Write(iulog,*) 'Current  unit = ', TRIM(State_Chm(LCHNK)%Spc_Units)
       Write(iulog,*) 'Expected unit = kg/ kg dry'
       CALL ENDRUN('Incorrect unit in GEOS-Chem State_Chm%Species')
    ENDIF

    call t_stopf( 'chemdr' )

    ! Save and write J-values to pbuf for HEMCO
    ! in HCO_IN_JNO2, HCO_IN_JOH
    FieldName = 'HCO_IN_JNO2'
    tmpIdx = pbuf_get_index(FieldName, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

       ! RXN_NO2: NO2 + hv --> NO  + O
       pbuf_i(:nY) = ZPJ(1,RXN_NO2,1,:nY)

       pbuf_chnk => NULL()
       pbuf_i    => NULL()
    ENDIF

    FieldName = 'HCO_IN_JOH'
    tmpIdx = pbuf_get_index(FieldName, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(FieldName)
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

       ! RXN_O3_1: O3  + hv --> O2  + O
       pbuf_i(:nY) = ZPJ(1,RXN_O3_1,1,:nY)
       pbuf_chnk => NULL()
       pbuf_i   => NULL()
    ENDIF

    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          vmr1(:nY,:nZ,N) = State_Chm(LCHNK)%Species(M)%Conc(1,:nY,nZ:1:-1) * &
                            MWDry / adv_mass(N)
       ELSEIF ( M < 0 ) THEN
          vmr1(:nY,:nZ,N) = state%q(:nY,:nZ,-M) * &
                            MWDry / adv_mass(N)
       ENDIF
    ENDDO

    !==============================================================
    ! ***** M A M   G A S - A E R O S O L   E X C H A N G E *****
    !==============================================================

#if defined( MODAL_AERO )
    ! Repartition SO4 into H2SO4 and so4_a*
    IF ( l_H2SO4 > 0 .AND. l_SO4 > 0 ) THEN
       P = l_H2SO4
       ! SO4_gasRatio is mol(SO4) (gaseous) / mol(SO4) (gaseous+aerosol)
       vmr1(:nY,:nZ,P) = SO4_gasRatio(:nY,:nZ) * vmr1(:nY,:nZ,l_SO4)
       ! binRatio is mol(SO4) (current bin) / mol(SO4) (all bins)
       DO M = 1, ntot_amode
          N = lptr_so4_a_amode(M)
          IF ( N <= 0 ) CYCLE
          P = mapCnst(N)
          vmr1(:nY,:nZ,P) = vmr1(:nY,:nZ,l_SO4)                &
                          * ( 1.0_r8 - SO4_gasRatio(:nY,:nZ) ) &
                          * binRatio(iSulf(M),M,:nY,:nZ)
       ENDDO
    ENDIF

    ! Amount of chemically-produced H2SO4 (mol/mol)
    ! This is archived from fullchem_mod.F90 using SO2 + OH rate from KPP (hplin, 1/25/23)
    del_h2so4_gasprod(:nY,:nZ) = State_Chm(LCHNK)%H2SO4_PRDR(1,:nY,nZ:1:-1)

    call aero_model_gasaerexch( loffset           = iFirstCnst - 1,         &
                                ncol              = NCOL,                   &
                                lchnk             = LCHNK,                  &
                                troplev           = Trop_Lev(:),            &
                                delt              = dT,                     &
                                reaction_rates    = reaction_rates,         &
                                tfld              = state%t(:,:),           &
                                pmid              = state%pmid(:,:),        &
                                pdel              = state%pdel(:,:),        &
                                mbar              = mBar,                   &
                                relhum            = relHum(:,:),            &
                                zm                = state%zm(:,:),          &
                                qh2o              = qH2O(:,:),              &
                                cwat              = cldW,                   &
                                cldfr             = cldFrc,                 &
                                cldnum            = nCldWtr,                &
                                airdens           = invariants(:,:,indexm), &
                                invariants        = invariants,             &
                                del_h2so4_gasprod = del_h2so4_gasprod,      &
                                vmr0              = vmr0,                   &
                                vmr               = vmr1,                   &
                                pbuf              = pbuf )

    ! Repartition MAM SOAs following mapping:
    ! TSOA0 + ASOAN + SOAIE + SOAGX -> soa1_a* + soa2_a*
    ! TSOA1 + ASOA1                 -> soa3_a*
    ! TSOA2 + ASOA2                 -> soa4_a*
    ! TSOA3 + ASOA3                 -> soa5_a*
    ! TSOG0                         -> SOAG0 + SOAG1
    ! TSOG1 + ASOG1                 -> SOAG2
    ! TSOG2 + ASOG2                 -> SOAG3
    ! TSOG3 + ASOG3                 -> SOAG4

    ! Deal with aerosol SOA species
    ! First deal with lowest two volatility bins
    ! Only map TOSA0 (K1) and ASOAN (K2) to soa1_ and soa2_, according to Fritz et al.
    ! SOAIE (K3) and SOAGX (K4) were mapped in the code but are inconsistent with the model description paper.
    speciesName_1 = 'TSOA0'
    speciesName_2 = 'ASOAN'
    speciesName_3 = 'SOAIE'
    speciesName_4 = 'SOAGX'
    K1 = get_spc_ndx(TRIM(speciesName_1), ignore_case=.true.)
    K2 = get_spc_ndx(TRIM(speciesName_2), ignore_case=.true.)
    K3 = get_spc_ndx(TRIM(speciesName_3), ignore_case=.true.)
    K4 = get_spc_ndx(TRIM(speciesName_4), ignore_case=.true.)
    bulkMass(:nY,:nZ) = 0.0e+00_r8
    DO iBin = 1, 2
       DO M = 1, ntot_amode
          N = lptr2_soa_a_amode(M,iBin)
          IF ( N <= 0 ) CYCLE
          bulkMass(:nY,:nZ) = bulkMass(:nY,:nZ) + state%q(:nY,:nZ,N)
       ENDDO
    ENDDO
    DO iBin = 1, 2
       DO M = 1, ntot_amode
          N = lptr2_soa_a_amode(M,iBin)
          IF ( N <= 0 ) CYCLE
          P = mapCnst(N)
          IF ( P > 0 .AND. K1 > 0 .AND. K2 > 0 ) THEN
             vmr1(:nY,:nZ,P) = state%q(:nY,:nZ,N) / bulkMass(:nY,:nZ) &
                             * (vmr1(:nY,:nZ,K1) + vmr1(:nY,:nZ,K2))
          ENDIF
       ENDDO
    ENDDO

    ! Now deal with other volatility bins
    DO iBin = 3, nsoa
       IF ( iBin == 3 ) THEN
          speciesName_1 = 'TSOA1'
          speciesName_2 = 'ASOA1'
       ELSEIF ( iBin == 4 ) THEN
          speciesName_1 = 'TSOA2'
          speciesName_2 = 'ASOA2'
       ELSEIF ( iBin == 5 ) THEN
          speciesName_1 = 'TSOA3'
          speciesName_2 = 'ASOA3'
       ENDIF
       K1 = get_spc_ndx(TRIM(speciesName_1), ignore_case=.true. )
       K2 = get_spc_ndx(TRIM(speciesName_2), ignore_case=.true. )
       bulkMass(:nY,:nZ) = 0.0e+00_r8
       DO M = 1, ntot_amode
          N = lptr2_soa_a_amode(M,iBin)
          IF ( N <= 0 ) CYCLE
          bulkMass(:nY,:nZ) = bulkMass(:nY,:nZ) + state%q(:nY,:nZ,N)
       ENDDO
       DO M = 1, ntot_amode
          N = lptr2_soa_a_amode(M,iBin)
          IF ( N <= 0 ) CYCLE
          P = mapCnst(N)
          IF ( P > 0 .AND. K1 > 0 .AND. K2 > 0 ) THEN
             vmr1(:nY,:nZ,P) = state%q(:nY,:nZ,N) / bulkMass(:nY,:nZ) &
                             * (vmr1(:nY,:nZ,K1) + vmr1(:nY,:nZ,K2))
          ENDIF
       ENDDO
    ENDDO

    ! Now deal with gaseous SOA species
    ! Deal with lowest two volatility bins - TSOG0 corresponds to SOAG0 and SOAG1
    speciesName_1 = 'TSOG0'
    K1 = get_spc_ndx(TRIM(speciesName_1), ignore_case=.true.)
    N = lptr2_soa_g_amode(1)
    P = mapCnst(N)
    !                                        current mode        other modes (this mapping was verified to be correct.)
    vmr1(:nY,:nZ,P) = vmr0(:nY,:nZ,P) / (vmr0(:nY,:nZ,P) + vmr0(:nY,:nZ,mapCnst(lptr2_soa_g_amode(2)))) &
                    * vmr1(:nY,:nZ,K1)
    N = lptr2_soa_g_amode(2)
    P = mapCnst(N)
    vmr1(:nY,:nZ,P) = vmr0(:nY,:nZ,P) / (vmr0(:nY,:nZ,P) + vmr0(:nY,:nZ,mapCnst(lptr2_soa_g_amode(1)))) &
                    * vmr1(:nY,:nZ,K1)

    ! Deal with other volatility bins
    DO iBin = 3, nsoa
       N = lptr2_soa_g_amode(iBin)
       P = mapCnst(N)
       IF ( iBin == 3 ) THEN
          speciesName_1 = 'TSOG1'
          speciesName_2 = 'ASOG1'
       ELSEIF ( iBin == 4 ) THEN
          speciesName_1 = 'TSOG2'
          speciesName_2 = 'ASOG2'
       ELSEIF ( iBin == 5 ) THEN
          speciesName_1 = 'TSOG3'
          speciesName_2 = 'ASOG3'
       ENDIF
       K1 = get_spc_ndx(TRIM(speciesName_1), ignore_case=.true.)
       K2 = get_spc_ndx(TRIM(speciesName_2), ignore_case=.true.)
       IF ( P > 0 .AND. K1 > 0 .AND. K2 > 0 ) vmr1(:nY,:nZ,P) = vmr1(:nY,:nZ,K1) + vmr1(:nY,:nZ,K2)
    ENDDO

#endif

    !==============================================================
    ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
    !==============================================================
    IF ( Input_Opt%LWetD ) THEN

       IF ( gas_wetdep_method == 'NEU' ) THEN
          CALL Neu_wetdep_tend( LCHNK       = LCHNK,      &
                                NCOL        = NCOL,       &
                                mmr         = state%q,    &
                                pmid        = state%pmid, &
                                pdel        = state%pdel, &
                                zint        = state%zi,   &
                                tfld        = state%t,    &
                                delt        = dT,         &
                                prain       = PRain,      &
                                nevapr      = NEvapr,     &
                                cld         = cldFrc,     &
                                cmfdqr      = cmfdqr,     &
                                wd_tend     = ptend%q,    &
                                wd_tend_int = wetdepflx  )
       ELSE
          ErrMsg = 'Unknown gas_wetdep_method '//TRIM(gas_wetdep_method)
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

    ENDIF

    !==============================================================
    ! ***** B O U N D A R Y   C O N D I T I O N S            *****
    !==============================================================
    ! Set boundary conditions of long-lived species (most likely
    ! CH4, OCS, N2O, CFC11, CFC12).
    ! Note: This will overwrite the UCX boundary conditions

    CALL flbc_set( vmr1(:nY,:nZ,:), nY, LCHNK, mapCnst )

    IF ( ghg_chem ) THEN
       CALL ghg_chem_set_flbc( vmr1, nY )
    ENDIF

    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M <= 0 ) CYCLE
       State_Chm(LCHNK)%Species(M)%Conc(1,:nY,nZ:1:-1) = vmr1(:nY,:nZ,N) * &
                        adv_mass(N) / MWDry
    ENDDO

    ! Make sure State_Chm(LCHNK) is back in kg/kg dry!
    IF ( TRIM(State_Chm(LCHNK)%Spc_Units) /= 'kg/kg dry' ) THEN
       Write(iulog,*) 'Current  unit = ', TRIM(State_Chm(LCHNK)%Spc_Units)
       Write(iulog,*) 'Expected unit = kg/ kg dry'
       CALL ENDRUN('Incorrect unit in GEOS-Chem State_Chm%Species')
    ENDIF

    ! Reset H2O MMR to the initial value (no chemistry tendency in H2O just yet)
    State_Chm(LCHNK)%Species(iH2O)%Conc(1,:,:) = MMR_Beg(:,:,iH2O)

    ! Store unadvected species data
    SlsData = 0.0e+0_r8
    DO N = 1, nSls
       M = map2GC_Sls(N)
       IF ( M <= 0 ) CYCLE
       SlsData(:nY,nZ:1:-1,N) = REAL(State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ),r8)
    ENDDO
    CALL set_short_lived_species_gc( SlsData, LCHNK, nY, pbuf )

    ! Apply tendencies to GEOS-Chem species
    DO N = 1, pcnst
       M = map2GC(N)
       IF ( M <= 0 ) CYCLE
       ! Add change in mass mixing ratio to tendencies.
       ! For NEU wet deposition, the wet removal rates are added to
       ! ptend.
       MMR_End(:nY,:nZ,M)     = REAL(State_Chm(LCHNK)%Species(M)%Conc(1,:nY,:nZ),r8)
       ptend%q(:nY,nZ:1:-1,N) = ptend%q(:nY,nZ:1:-1,N) &
                              + (MMR_End(:nY,:nZ,M)-MMR_Beg(:nY,:nZ,M))/dT
    ENDDO

#if defined( MODAL_AERO )
    ! Here apply tendencies to MAM aerosols
    ! Initial mass in bin SM is stored as state%q(N)
    ! Final mass in bin SM is stored as binRatio(SM,M) * State_Chm(P)
    !
    ! We decide to apply chemical tendencies to all MAM aerosols,
    ! except so4, for which the chemically-produced sulfate gets
    ! partitioned in aero_model_gasaerexch.
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          N = lmassptr_amode(SM,M)
          P = mapCnst(N)
          IF ( P <= 0 ) CYCLE
          ! Apply tendency from MAM gasaerexch
          ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N) &
                             + (vmr1(:nY,:nZ,P) - vmr0(:nY,:nZ,P))/dT &
                                  * adv_mass(P) / MWDry
          P = map2MAM4(SM,M)
          IF ( P <= 0 ) CYCLE
          K = map2GC(P)
          IF ( K <= 0 .or. K == iSO4 ) CYCLE
          ! Apply MAM4 chemical tendencies owing to GEOS-Chem aerosol processing
          ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N)                                  &
                             + (binRatio(SM,M,:nY,:nZ) *                           &
                                REAL(State_Chm(LCHNK)%Species(K)%Conc(1,:nY,nZ:1:-1),r8) &
                                  * adv_mass(mapCnst(N)) / adv_mass(mapCnst(P))    &
                                - state%q(:nY,:nZ,N))/dT
       ENDDO
       N = numptr_amode(M)
       P = mapCnst(N)
       IF ( P <= 0 ) CYCLE
       ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N) &
                          + (vmr1(:nY,:nZ,P) - vmr0(:nY,:nZ,P))/dT &
                               * adv_mass(P) / MWDry
    ENDDO
    N = cH2SO4
    P = l_H2SO4
    IF ( P > 0 ) THEN
       ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N) &
                          + (vmr1(:nY,:nZ,P) - vmr0(:nY,:nZ,P))/dT &
                               * adv_mass(P) / MWDry
    ENDIF
    DO iBin = 1, nsoa
       N = lptr2_soa_g_amode(iBin)
       P = mapCnst(N)
       IF ( P > 0 ) THEN
          ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N) &
                             + (vmr1(:nY,:nZ,P) - vmr0(:nY,:nZ,P))/dT &
                                  * adv_mass(P) / MWDry
       ENDIF
    ENDDO
#endif

    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          mmr_tend(:nY,:nZ,N) = ( REAL(State_Chm(LCHNK)%Species(M)%Conc(1,:nY,nZ:1:-1),r8) - mmr_tend(:nY,:nZ,N) ) / dT
       ELSEIF ( M < 0 ) THEN
          mmr_tend(:nY,:nZ,N) = ptend%q(:nY,:nZ,-M)
       ENDIF
    ENDDO

    IF ( Input_Opt%applyQtend ) THEN
       ! Apply GEOS-Chem's H2O mixing ratio tendency to CAM's specific humidity
       ! This requires to set lq(cQ) = lq(cH2O) ( = .True. )
       ptend%q(:,:,cQ) = ptend%q(:,:,cH2O)
    ENDIF

    CALL GC_Diagnostics_Calc( Input_Opt  = Input_Opt,         &
                              State_Chm  = State_Chm(LCHNK),  &
                              State_Diag = State_Diag(LCHNK), &
                              State_Grid = State_Grid(LCHNK), &
                              State_Met  = State_Met(LCHNK),  &
                              cam_in     = cam_in,            &
                              state      = state,             &
                              mmr_tend   = mmr_tend,          &
                              LCHNK      = LCHNK             )

    CALL Set_Diagnostics_EndofTimestep( Input_Opt  = Input_Opt,         &
                                        State_Chm  = State_Chm(LCHNK),  &
                                        State_Diag = State_Diag(LCHNK), &
                                        State_Grid = State_Grid(LCHNK), &
                                        State_Met  = State_Met(LCHNK),  &
                                        RC         = RC                )


    IF ( State_Diag(LCHNK)%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic(  Input_Opt  = Input_Opt,         &
                                     State_Chm  = State_Chm(LCHNK),  &
                                     State_Diag = State_Diag(LCHNK), &
                                     State_Grid = State_Grid(LCHNK), &
                                     State_Met  = State_Met(LCHNK),  &
                                     RC         = RC                )
    ENDIF

    ! Compute new GEOS-Chem diagnostics into CESM History (hplin, 10/31/22)
    ! Note that the containers (data pointers) actually need to be updated every time step,
    ! because the State_Chm(LCHNK) target changes. There is some registry lookup overhead
    ! but mitigated by a check to the history field activeness. (hplin, 11/1/22)
    CALL HistoryExports_SetDataPointers(rootChunk,            &
                                        HistoryConfig,        State_Chm(LCHNK), &
                                        State_Grid(LCHNK),                      &
                                        State_Diag(LCHNK),    State_Met(LCHNK), &
                                        RC)

    CALL CopyGCStates2Exports( am_I_Root     = rootChunk,         &
                               Input_Opt     = Input_Opt,         &
                               State_Grid    = State_Grid(LCHNK), &
                               HistoryConfig = HistoryConfig,     &
                               LCHNK         = LCHNK,             &
                               RC            = RC             )

    IF ( ghg_chem ) THEN
       ptend%lq(1) = .True.
       CALL outfld( 'CT_H2O_GHG', ptend%q(:,:,1), PCOLS, LCHNK )
    ENDIF

    !! Debug statements
    !! Ozone tendencies
    !IF ( rootChunk ) THEN
    !   Write(iulog,*) " MMR_Beg = ", MMR_Beg(1,:,iO3)
    !   Write(iulog,*) " MMR_End = ", MMR_End(1,:,iO3)
    !ENDIF

    IF (PRESENT(fh2o)) THEN
       fh2o(:nY) = 0.0e+0_r8
       !DO L = 1, nZ
       !   fh2o(:nY) = fh2o(:nY) + ptend%q(:nY,L,iH2O)*state%pdel(:nY,L)/Gravit
       !ENDDO
    ENDIF

    ! Nullify all pointers
    Nullify(PblH    )
    Nullify(Fsds    )
    Nullify(PRain   )
    Nullify(LsFlxSnw)
    Nullify(LsFlxPrc)
    Nullify(cldTop  )
    Nullify(cldFrc  )
    Nullify(NEvapr  )
    Nullify(cmfdqr  )

    IF ( rootChunk ) WRITE(iulog,*) 'GEOS-Chem Chemistry step ', iStep, ' completed'
    IF ( lastChunk ) WRITE(iulog,*) 'Chemistry completed on all chunks of root CPU'
    IF ( FIRST ) THEN
        FIRST = .false.
    ENDIF

  end subroutine chem_timestep_tend

  !================================================================================================
  ! subroutine chem_init_cnst
  !================================================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    CHARACTER(LEN=*), INTENT(IN)  :: name       !  constituent name
    REAL(r8),         INTENT(IN)  :: latvals(:) ! lat in degrees (NCOL)
    REAL(r8),         INTENT(IN)  :: lonvals(:) ! lon in degrees (NCOL)
    LOGICAL,          INTENT(IN)  :: mask(:)    ! Only initialize where .true.
    REAL(r8),         INTENT(OUT) :: q(:,:)     ! kg tracer/kg dry air (NCOL, PVER)
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    INTEGER  :: ilev, nlev, M
    REAL(r8) :: QTemp, Min_MMR

    nlev = SIZE(q, 2)

    ! Retrieve a "background value" for this from the database
    Min_MMR = 1.0e-38_r8
    CALL cnst_get_ind(TRIM(name), M, abort=.False.)
    IF ( M > 0 ) Min_MMR = ref_MMR(M)

    DO ilev = 1, nlev
       WHERE(mask)
          ! Set to the minimum mixing ratio
          q(:,ilev) = Min_MMR
       END WHERE
    ENDDO

  end subroutine chem_init_cnst

  !================================================================================================
  ! subroutine chem_final
  !================================================================================================
  subroutine chem_final

    ! CAM modules
    use short_lived_species,    only : short_lived_species_final

    ! GEOS-Chem interface modules in CAM
    use geoschem_emissions_mod, only : GC_Emissions_Final
    use geoschem_history_mod,   only : Destroy_HistoryConfig

    ! GEOS-Chem modules
    use Aerosol_Mod,     only : Cleanup_Aerosol
    use Carbon_Mod,      only : Cleanup_Carbon
    use CMN_FJX_Mod,     only : Cleanup_CMN_FJX
    use Drydep_Mod,      only : Cleanup_Drydep
    use Dust_Mod,        only : Cleanup_Dust
    use Error_Mod,       only : Cleanup_Error
    use Fullchem_Mod,    only : Cleanup_FullChem
    use Input_Opt_Mod,   only : Cleanup_Input_Opt
    use Linear_Chem_Mod, only : Cleanup_Linear_Chem
    use Pressure_Mod,    only : Cleanup_Pressure
    use Seasalt_Mod,     only : Cleanup_Seasalt
    use State_Chm_Mod,   only : Cleanup_State_Chm
    use State_Diag_Mod,  only : Cleanup_State_Diag
    use State_Grid_Mod,  only : Cleanup_State_Grid
    use State_Met_Mod,   only : Cleanup_State_Met
    use Sulfate_Mod,     only : Cleanup_Sulfate

    ! Local variables
    INTEGER  :: I, RC

    ! Destroy the history interface between GC States and CAM exports
    CALL Destroy_HistoryConfig(masterproc, HistoryConfig, RC)

    ! Finalize GEOS-Chem

    CALL Cleanup_Aerosol
    CALL Cleanup_Carbon
    CALL Cleanup_Drydep
    CALL Cleanup_Dust
    CALL Cleanup_FullChem( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_FullChem"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_Pressure
    CALL Cleanup_Seasalt
    CALL Cleanup_Sulfate
    CALL Cleanup_Linear_Chem

    CALL GC_Emissions_Final

    CALL short_lived_species_final()

    CALL Cleanup_CMN_FJX( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_FJX"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Cleanup Input_Opt
    CALL Cleanup_Input_Opt( Input_Opt, RC )

    ! Loop over each chunk and cleanup the variables
    DO I = BEGCHUNK, ENDCHUNK
       CALL Cleanup_State_Chm ( State_Chm(I),  RC )
       CALL Cleanup_State_Diag( State_Diag(I), RC )
       CALL Cleanup_State_Grid( State_Grid(I), RC )
       CALL Cleanup_State_Met ( State_Met(I),  RC )
    ENDDO
    CALL Cleanup_Error

    ! Finally deallocate state variables
    IF ( ALLOCATED( State_Chm  ) )    DEALLOCATE( State_Chm  )
    IF ( ALLOCATED( State_Diag ) )    DEALLOCATE( State_Diag )
    IF ( ALLOCATED( State_Grid ) )    DEALLOCATE( State_Grid )
    IF ( ALLOCATED( State_Met ) )     DEALLOCATE( State_Met )

    IF ( ALLOCATED( slvd_Lst     ) )  DEALLOCATE( slvd_Lst     )

    RETURN

  end subroutine chem_final

  !================================================================================================
  ! subroutine chem_init_restart
  !================================================================================================
  subroutine chem_init_restart(File)

    ! CAM modules
    use pio,              only : file_desc_t
    use tracer_cnst,      only : init_tracer_cnst_restart
    use tracer_srcs,      only : init_tracer_srcs_restart

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    WRITE(iulog,'(a)') 'chem_init_restart: init restarts for tracer sources and offline fields'

    !
    ! data for offline tracers
    !
    call init_tracer_cnst_restart(File)
    call init_tracer_srcs_restart(File)
    !call init_linoz_data_restart(File)

  end subroutine chem_init_restart

  !================================================================================================
  ! subroutine chem_write_restart
  !================================================================================================
  subroutine chem_write_restart( File )

    ! CAM modules
    use pio,         only : file_desc_t
    use tracer_cnst, only : write_tracer_cnst_restart
    use tracer_srcs, only : write_tracer_srcs_restart

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    WRITE(iulog,'(a)') 'chem_write_restart: writing restarts for tracer sources and offline fields'

    ! data for offline tracers
    call write_tracer_cnst_restart(File)
    call write_tracer_srcs_restart(File)

  end subroutine chem_write_restart

  !================================================================================================
  ! subroutine chem_read_restart
  !================================================================================================
  subroutine chem_read_restart( File )

    ! CAM modules
    use pio,         only : file_desc_t
    use tracer_cnst, only : read_tracer_cnst_restart
    use tracer_srcs, only : read_tracer_srcs_restart

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    WRITE(iulog,'(a)') 'GCCALL CHEM_READ_RESTART'

    ! data for offline tracers
    call read_tracer_cnst_restart(File)
    call read_tracer_srcs_restart(File)

  end subroutine chem_read_restart

  !================================================================================================
  ! subroutine chem_emissions
  !================================================================================================
  subroutine chem_emissions( state, cam_in, pbuf )

    ! CAM modules
    use camsrfexch,          only : cam_in_t
    use physics_buffer,      only : physics_buffer_desc

    TYPE(physics_state),    INTENT(IN)    :: state   ! Physics state variables
    TYPE(cam_in_t),         INTENT(INOUT) :: cam_in  ! import state
    TYPE(physics_buffer_desc), pointer    :: pbuf(:) ! Physics buffer in chunk, for HEMCO

    INTEGER :: M, N
    INTEGER :: nY
    LOGICAL :: rootChunk

    nY = state%NCOL    ! number of atmospheric columns on this chunk
    rootChunk = ( MasterProc .and. (state%LCHNK .eq. BEGCHUNK) )

    ! Reset surface fluxes
    DO M = iFirstCnst, pcnst
       !N = map2chm(M)
       !IF ( N > 0 ) cam_in%cflx(1:nY,N) = 0.0e+0_r8
       cam_in%cflx(1:nY,M) = 0.0e+0_r8
    ENDDO

  end subroutine chem_emissions

end module chemistry
