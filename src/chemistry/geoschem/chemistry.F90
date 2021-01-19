!================================================================================================
! This is the "GEOS-Chem" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only : r8 => shr_kind_r8, shr_kind_cl
  use physics_types,       only : physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,      only : physics_buffer_desc
  use ppgrid,              only : begchunk, endchunk, pcols
  use ppgrid,              only : pver, pverp
  use constituents,        only : pcnst, cnst_add, cnst_get_ind
  use constituents,        only : cnst_name
  use shr_const_mod,       only : molw_dryair=>SHR_CONST_MWDAIR
  use seq_drydep_mod,      only : nddvels => n_drydep, drydep_list
  use spmd_utils,          only : MasterProc, myCPU=>Iam, nCPUs=>npes
  use cam_logfile,         only : iulog
  use string_utils,        only : to_upper

  !--------------------------------------------------------------------
  ! Basic GEOS-Chem modules
  !--------------------------------------------------------------------
  USE DiagList_Mod,        ONLY : DgnList       ! Derived type for diagnostics list
  USE TaggedDiagList_Mod,  ONLY : TaggedDgnList ! Derived type for tagged diagnostics list
  USE Input_Opt_Mod,       ONLY : OptInput      ! Derived type for Input Options
  USE State_Chm_Mod,       ONLY : ChmState      ! Derived type for Chemistry State object
  USE State_Diag_Mod,      ONLY : DgnState      ! Derived type for Diagnostics State object
  USE State_Grid_Mod,      ONLY : GrdState      ! Derived type for Grid State object
  USE State_Met_Mod,       ONLY : MetState      ! Derived type for Meteorology State object
  USE Species_Mod,         ONLY : Species       ! Derived type for Species object
  USE GC_Environment_Mod                        ! Runtime GEOS-Chem environment
  USE ErrCode_Mod                               ! Error codes for success or failure
  USE Error_Mod                                 ! For error checking

  !-----------------------------------------------------------------
  ! Parameters to define floating-point variables
  !-----------------------------------------------------------------
  USE PRECISION_MOD,       ONLY : fp, f4     ! Flexible precision

  use chem_mods,           only : nSlvd, slvd_Lst, slvd_ref_MMR

  ! Exit routine in CAM
  use cam_abortutils,      only : endrun

  use chem_mods,           only : nTracersMax
  use chem_mods,           only : nTracers
  use chem_mods,           only : gas_pcnst
  use chem_mods,           only : tracerNames, tracerLongNames
  use chem_mods,           only : adv_mass
  use chem_mods,           only : mwRatio
  use chem_mods,           only : ref_MMR
  use chem_mods,           only : iFirstCnst
  use chem_mods,           only : nSlsMax
  use chem_mods,           only : nSls
  use chem_mods,           only : slsNames, slsLongNames
  use chem_mods,           only : sls_ref_MMR
  use chem_mods,           only : slsmwRatio
  use chem_mods,           only : nAerMax
  use chem_mods,           only : nAer
  use chem_mods,           only : aerNames
  use chem_mods,           only : aerAdvMass
  use chem_mods,           only : map2GC, map2GCinv
  use chem_mods,           only : map2GC_Sls
  use chem_mods,           only : map2chm
  use chem_mods,           only : map2Idx
  use chem_mods,           only : map2MAM4

  use mo_tracname,         only : solsym

  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  ! Public interfaces
  !
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

  ! Location of valid input.geos
  CHARACTER(LEN=500) :: inputGeosPath

  ! Location of valid species_database.yml
  CHARACTER(LEN=500) :: speciesDBPath

  ! Location of chemistry input (for now)
  CHARACTER(LEN=500) :: chemInputsDir

  !-----------------------------
  ! Derived type objects
  !-----------------------------
  TYPE(OptInput)             :: Input_Opt       ! Input Options object
  TYPE(ChmState),ALLOCATABLE :: State_Chm(:)    ! Chemistry State object
  TYPE(DgnState),ALLOCATABLE :: State_Diag(:)   ! Diagnostics State object
  TYPE(GrdState),ALLOCATABLE :: State_Grid(:)   ! Grid State object
  TYPE(MetState),ALLOCATABLE :: State_Met(:)    ! Meteorology State object
  TYPE(DgnList )             :: Diag_List       ! Diagnostics list object
  TYPE(TaggedDgnList )       :: TaggedDiag_List ! Tagged diagnostics list object

  type(physics_buffer_desc), pointer :: hco_pbuf2d(:,:)    ! Pointer to 2D pbuf

  ! Indices of critical species in GEOS-Chem
  INTEGER                    :: iH2O, iO3, iCH4, iCO, iNO, iOH
  INTEGER                    :: iO, iH, iO2, iPSO4
  REAL(r8)                   :: MWOH, MWPSO4, MWO3
  ! Indices of critical species in the constituent list
  INTEGER                    :: cQ, cH2O

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

  ! lightning
  REAL(r8)                   :: lght_no_prd_factor = 1._r8

  ! Strings
  CHARACTER(LEN=255)         :: ThisLoc
  CHARACTER(LEN=255)         :: ErrMsg

  ! Filenames to compute dry deposition velocities similarly to MOZART
  character(len=shr_kind_cl) :: clim_soilw_file = 'clim_soilw_file'
  character(len=shr_kind_cl) :: depvel_file     = ''
  character(len=shr_kind_cl) :: depvel_lnd_file = 'depvel_lnd_file'
  character(len=shr_kind_cl) :: season_wes_file = 'season_wes_file'

  character(len=shr_kind_cl) :: srf_emis_specifier(pcnst) = ''
  character(len=shr_kind_cl) :: ext_frc_specifier(pcnst) = ''

  character(len=24)          :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer                    :: srf_emis_cycle_yr  = 0
  integer                    :: srf_emis_fixed_ymd = 0
  integer                    :: srf_emis_fixed_tod = 0

  character(len=24)          :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer                    :: ext_frc_cycle_yr  = 0
  integer                    :: ext_frc_fixed_ymd = 0
  integer                    :: ext_frc_fixed_tod = 0


!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    use mo_chem_utls, only : utls_chem_is

    character(len=*), intent(in) :: name

    chem_is = utls_chem_is(name)

  end function chem_is

!================================================================================================

  subroutine chem_register

    use physics_buffer,      only : pbuf_add_field, dtype_r8
    use PhysConst,           only : MWDry

    use Short_Lived_Species, only : Register_Short_Lived_Species

    use State_Grid_Mod,      only : Init_State_Grid, Cleanup_State_Grid
    use State_Chm_Mod,       only : Init_State_Chm, Cleanup_State_Chm
    use State_Chm_Mod,       only : Ind_
    use Input_Opt_Mod,       only : Set_Input_Opt,  Cleanup_Input_Opt
    use CMN_SIZE_Mod,        only : Init_CMN_SIZE

    use mo_sim_dat,          only : set_sim_dat
    use mo_chem_utls,        only : get_spc_ndx
    use chem_mods,           only : drySpc_ndx
#if defined( MODAL_AERO_4MODE )
    use aero_model,          only : aero_model_register
    use modal_aero_data,     only : nspec_max
    use modal_aero_data,     only : ntot_amode, nspec_amode
    use modal_aero_data,     only : xname_massptr
#endif

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
    REAL(r8)                       :: cptmp
    REAL(r8)                       :: MWTmp
    REAL(r8)                       :: qmin
    REAL(r8)                       :: ref_VMR
    CHARACTER(LEN=128)             :: mixtype
    CHARACTER(LEN=128)             :: molectype
    CHARACTER(LEN=128)             :: lngName
    CHARACTER(LEN=64)              :: cnstName
    CHARACTER(LEN=64)              :: trueName
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

    ! SDE 2018-05-02: This seems to get called before anything else
    ! that includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    ! hplin 2020-05-16: Call set_sim_dat to populate chemistry constituent information
    ! from mo_sim_dat.F90 in other places. This is needed for HEMCO_CESM.
    CALL Set_sim_dat()

    ! Generate fake state_chm
    IO%Max_BPCH_Diag       = 1000
    IO%Max_AdvectSpc       = 500
    IO%Max_Families        = 250

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
    IO%LUCX                = .True.
    IO%LPRT                = .False.
    IO%N_Advect            = nTracers
    DO I = 1, nTracers
       IO%AdvectSpc_Name(I) = TRIM(tracerNames(I))
    ENDDO
    IO%SALA_rEdge_um(1)    = 0.01e+0_fp
    IO%SALA_rEdge_um(2)    = 0.50e+0_fp
    IO%SALC_rEdge_um(1)    = 0.50e+0_fp
    IO%SALC_rEdge_um(2)    = 8.00e+0_fp

    IO%SpcDatabaseFile     = TRIM(speciesDBPath)

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

    CALL Init_CMN_SIZE( Input_Opt = IO,  &
                        RC        = RC  )

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
    map2GC     = -1
    map2GCinv  = -1
    map2chm    = -1
    ref_MMR(:) = 0.0e+0_r8
    MWRatio(:) = 1.0e+0_r8
    tracerLongNames = ''

    DO I = 1, nTracersMax
       IF ( I .LE. nTracers ) THEN
          cnstName    = TRIM(tracerNames(I))
          trueName    = cnstName
          N           = Ind_(cnstName)
          ThisSpc     => SC%SpcData(N)%Info
          lngName     = TRIM(ThisSpc%FullName)
          MWTmp       = REAL(ThisSpc%MW_g,r8)
          ref_VMR     = REAL(ThisSpc%BackgroundVV,r8)
          ref_MMR(I)  = ref_VMR / (MWDry / MWTmp)
          ! This is required as we need to distinguish between MAM and GEOS-Chem aerosols
          ! (Both are included in aer_drydep_list)
          IF ( ThisSpc%Is_Gas .eqv. .False. ) THEN
             Write(cnstName, "(a,a)") 'GC_AER_', to_upper(TRIM(trueName))
          ENDIF
       ELSEIF ( I .LE. (nTracers + nAer) ) THEN
          ! Add MAM4 aerosols
          cnstName    = TRIM(aerNames(I - nTracers))
          lngName     = cnstName
          MWTmp       = aerAdvMass(I - nTracers)
          ref_MMR(I)  = 1.0e-38_r8
       ELSE
          cnstName    = TRIM(tracerNames(I))
          lngName     = cnstName
          MWTmp       = 1000.0e+0_r8 * (0.001e+0_r8)
          ref_MMR(I)  = 1.0e-38_r8
       ENDIF
       MWRatio(I) = MWDry/MWTmp
       tracerLongNames(I) = TRIM(lngName)

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
       If ( MasterProc ) Write(iulog,*) " Species = ", TRIM(cnstName)
       ! GEOS-Chem lumped species are not on restart file.
       ! Bromine, chlorine, iodine and halons species are missing
       ! from CESM restart file.
       ! These species will just be uniformily set to some low
       ! concentration.
       ! TMMF - 05/19/2020

       ! This is the call to add a "constituent"
       CALL cnst_add( cnstName, MWtmp, cptmp, qmin, N,        &
                      readiv=ic_from_cam2, mixtype=mixtype,   &
                      cam_outfld=camout, molectype=molectype, &
                      fixed_ubc=has_fixed_ubc,                &
                      fixed_ubflx=has_fixed_ubflx,            &
                      longname=TRIM(lngName)                 )

       IF ( iFirstCnst < 0 ) iFirstCnst = N

       ! Add to GC mapping. When starting a timestep, we will want to update the
       ! concentration of State_Chm(x)%Species(1,iCol,iLev,m) with data from
       ! constituent n
       M = Ind_(TRIM(tracerNames(I)))
       IF ( M > 0 ) THEN
          ! Map constituent onto GEOS-Chem tracer as indexed in State_Chm(LCHNK)%Species
          map2GC(N)    = M
          ! Map GEOS-Chem tracer onto constituent
          map2GCinv(M) = N
          ! Map constituent onto raw index
          map2Idx(N)   = I
       ENDIF
       ! Nullify pointer
       ThisSpc => NULL()
    ENDDO

    ! Now unadvected species
    map2GC_Sls = 0
    sls_ref_MMR(:) = 0.0e+0_r8
    SlsMWRatio(:)  = -1.0e+0_r8
    slsLongNames = ''
    DO I = 1, nSls
       N = Ind_(slsNames(I))
       IF ( N .GT. 0 ) THEN
          ThisSpc         => SC%SpcData(N)%Info
          MWTmp           = REAL(ThisSpc%MW_g,r8)
          ref_VMR         = REAL(ThisSpc%BackgroundVV,r8)
          lngName         = TRIM(ThisSpc%FullName)
          slsLongNames(I) = lngName
          sls_ref_MMR(I)  = ref_VMR / (MWDry / MWTmp)
          SlsMWRatio(I)   = MWDry / MWTmp
          map2GC_Sls(I)   = N
          ThisSpc         => NULL()
       ENDIF
    ENDDO

    ! Pass information to "short_lived_species" module
    slvd_ref_MMR(1:nSls) = sls_ref_MMR(1:nSls)
    CALL Register_Short_Lived_Species()
    ! More information:
    ! http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/node5.html

    DO N = 1, gas_pcnst
       ! Map solsym onto GEOS-Chem species
       map2chm(N) = Ind_(TRIM(solsym(N)))
       IF ( map2chm(N) < 0 ) THEN
          ! This is not a GEOS-Chem species and we thus map on constituents
          ! Most likely, these will be MAM aerosols
          ! We store the index as the opposite to not confuse with GEOS-Chem
          ! indices.
          CALL cnst_get_ind(TRIM(solsym(N)), I, abort=.True.)
          map2chm(N) = -I
       ENDIF
    ENDDO
    ! Get constituent index of specific humidity
    CALL cnst_get_ind('Q',   cQ,   abort=.True.)
    CALL cnst_get_ind('H2O', cH2O, abort=.True.)

    !==============================================================
    ! Get mapping between dry deposition species and species set
    !==============================================================

    nIgnored = 0

    DO N = 1, nddvels

       ! The species names need to be convert to upper case as,
       ! for instance, BR2 != Br2
       drySpc_ndx(N) = get_spc_ndx( to_upper(drydep_list(N)) )

       IF ( MasterProc .AND. ( drySpc_ndx(N) < 0 ) ) THEN
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

    ! Initialize indices
    map2MAM4(:,:) = -1

    DO M = 1, ntot_amode
       DO L = 1, nspec_amode(M)
          SELECT CASE ( to_upper(xname_massptr(L,M)(:3)) )
             CASE ( 'BC_' )
                SELECT CASE ( to_upper(xname_massptr(L,M)(4:5)) )
                   CASE ( 'A1' )
                       map2MAM4(L,M) = Ind_('BCPI')
                   CASE ( 'A4' )
                       map2MAM4(L,M) = Ind_('BCPO')
                END SELECT
             CASE ( 'DST' )
                SELECT CASE ( to_upper(xname_massptr(L,M)(5:6)) )
                   ! DST1 - Dust aerosol, Reff = 0.7 micrometers
                   ! DST2 - Dust aerosol, Reff = 1.4 micrometers
                   ! DST3 - Dust aerosol, Reff = 2.4 micrometers
                   ! DST4 - Dust aerosol, Reff = 4.5 micrometers
                   CASE ( 'A1' )
                       map2MAM4(L,M) = Ind_('DST1')
                   CASE ( 'A2' )
                       map2MAM4(L,M) = Ind_('DST1')
                   CASE ( 'A3' )
                       map2MAM4(L,M) = Ind_('DST4')
                END SELECT
             CASE ( 'SOA' )
                map2MAM4(L,M) = Ind_('SOAS')
             CASE ( 'SO4' )
                map2MAM4(L,M) = Ind_('SO4')
             CASE ( 'NCL' )
                SELECT CASE ( to_upper(xname_massptr(L,M)(5:6)) )
                   ! SALA - Fine (0.01-0.05 micros) sea salt aerosol
                   ! SALC - Coarse (0.5-8 micros) sea salt aerosol
                   CASE ( 'A1' )
                      map2MAM4(L,M) = Ind_('SALA')
                   CASE ( 'A2' )
                      map2MAM4(L,M) = Ind_('SALA')
                   CASE ( 'A3' )
                      map2MAM4(L,M) = Ind_('SALC')
                END SELECT
             CASE ( 'POM' )
                SELECT CASE ( to_upper(xname_massptr(L,M)(5:6)) )
                   CASE ( 'A1' )
                      map2MAM4(L,M) = Ind_('OCPI')
                   CASE ( 'A4' )
                      map2MAM4(L,M) = Ind_('OCPO')
                END SELECT
          END SELECT
       ENDDO
    ENDDO

#endif

    !==============================================================
    ! Print summary
    !==============================================================

    IF ( MasterProc ) THEN
       Write(iulog,'(/, a)') '### Summary of GEOS-Chem species: '
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

  end subroutine chem_register

  subroutine chem_readnl(nlfile)

    use cam_abortutils,  only : endrun
    use units,           only : getunit, freeunit
    use namelist_utils,  only : find_group_name
#if defined( MODAL_AERO_4MODE )
    use aero_model,      only : aero_model_readnl
    use dust_model,      only : dust_readnl
#endif
    use gas_wetdep_opts, only : gas_wetdep_readnl
#ifdef SPMD
    use mpishorthand
#endif
    use gckpp_Model,     only : nSpec, Spc_Names
    use chem_mods,       only : drySpc_ndx

    ! args
    CHARACTER(LEN=*), INTENT(IN) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    INTEGER                      :: I, N
    INTEGER                      :: UNITN, IERR
    CHARACTER(LEN=500)           :: line
    LOGICAL                      :: menuFound
    LOGICAL                      :: validSLS

    ! The following files are required to compute land maps, required to perform
    ! aerosol dry deposition
    namelist /chem_inparm/ clim_soilw_file,    &
                           depvel_file,        &
                           lght_no_prd_factor, &
                           depvel_lnd_file,    &
                           ext_frc_specifier,  &
                           ext_frc_type,       &
                           ext_frc_cycle_yr,   &
                           ext_frc_fixed_ymd,  &
                           ext_frc_fixed_tod,  &
                           season_wes_file,    &
                           srf_emis_specifier, &
                           srf_emis_cycle_yr,  &
                           srf_emis_fixed_ymd, &
                           srf_emis_fixed_tod, &
                           srf_emis_type

    inputGeosPath='/glade/u/home/fritzt/input.geos.template'
    speciesDBPath='/glade/u/home/fritzt/species_database.yml'
    chemInputsDir='/glade/p/univ/umit0034/ExtData/CHEM_INPUTS/'

    ALLOCATE(drySpc_ndx(nddvels), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failed to allocate drySpc_ndx')

#if defined( MODAL_AERO_4MODE )
    !==============================================================
    ! Get names and molar weights of aerosols in MAM4
    !==============================================================

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

       unitn = getunit()

       !==============================================================
       ! Opening input.geos and go to ADVECTED SPECIES MENU
       !==============================================================

       OPEN( unitn, FILE=TRIM(inputGeosPath), STATUS='OLD', IOSTAT=IERR )
       IF (IERR .NE. 0) THEN
          CALL ENDRUN('chem_readnl: ERROR opening input.geos')
       ENDIF

       ! Go to ADVECTED SPECIES MENU
       menuFound = .False.
       DO WHILE ( .NOT. menuFound )
          READ( unitn, '(a)', IOSTAT=IERR ) line
          IF ( IERR .NE. 0 ) THEN
              CALL ENDRUN('chem_readnl: ERROR finding advected species menu')
          ELSEIF ( INDEX(line, 'ADVECTED SPECIES MENU') > 0 ) THEN
              menuFound = .True.
          ENDIF
       ENDDO

       !==============================================================
       ! Read list of GEOS-Chem tracers
       !==============================================================

       DO
          ! Read line
          READ(unitn,'(26x,a)', IOSTAT=IERR) line

          IF ( INDEX( TRIM(line), '---' ) > 0 ) EXIT

          nTracers = nTracers + 1
          tracerNames(nTracers) = TRIM(line)

       ENDDO

       CLOSE(unitn)
       CALL freeunit(unitn)

       ! Assign remaining tracers dummy names
       DO I = (nTracers+1), nTracersMax
          WRITE(tracerNames(I),'(a,I0.4)') 'GCTRC_', I
       ENDDO

       !==============================================================
       ! Now go through the KPP mechanism and add any species not
       ! implemented by the tracer list in input.geos
       !==============================================================

       IF ( nSpec > nSlsMax ) THEN
          CALL ENDRUN('chem_readnl: too many species - increase nSlsmax')
       ENDIF

       nSls = 0
       DO I = 1, nSpec
          ! Get the name of the species from KPP
          line = ADJUSTL(TRIM(Spc_Names(I)))
          ! Only add this
          validSLS = ( .NOT. ANY(TRIM(line) .EQ. tracerNames) )
          IF (validSLS) THEN
             ! Genuine new short-lived species
             nSls = nSls + 1
             slsNames(nSls) = TRIM(line)
          ENDIF
       ENDDO

       !==============================================================

       unitn = getunit()
       OPEN( unitn, FILE=TRIM(nlfile), STATUS='old' )
       CALL find_group_name(unitn, 'chem_inparm', STATUS=IERR)
       IF (IERR == 0) THEN
          READ(unitn, chem_inparm, IOSTAT=IERR)
          IF (IERR /= 0) THEN
             CALL endrun('chem_readnl: ERROR reading namelist')
          ENDIF
       ENDIF
       CLOSE(unitn)
       CALL freeunit(unitn)

    ENDIF

    !==================================================================
    ! Broadcast to all processors
    !==================================================================

#if defined( SPMD )
    CALL MPIBCAST ( nTracers,    1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST ( tracerNames, LEN(tracerNames(1))*nTracersMax, MPICHAR, 0, MPICOM )
    CALL MPIBCAST ( nSls,        1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST ( slsNames,    LEN(slsNames(1))*nSlsMax,        MPICHAR, 0, MPICOM )

    ! The following files are required to compute land maps, required to perform
    ! aerosol dry deposition
    CALL MPIBCAST (depvel_lnd_file, LEN(depvel_lnd_file), MPICHAR, 0, MPICOM)
    CALL MPIBCAST (clim_soilw_file, LEN(clim_soilw_file), MPICHAR, 0, MPICOM)
    CALL MPIBCAST (season_wes_file, LEN(season_wes_file), MPICHAR, 0, MPICOM)

    CALL MPIBCAST (lght_no_prd_factor, 1,                                MPIR8,   0, MPICOM)
    CALL MPIBCAST (depvel_file,        LEN(depvel_file),                 MPICHAR, 0, MPICOM)
    CALL MPIBCAST (srf_emis_specifier, LEN(srf_emis_specifier(1))*pcnst, MPICHAR, 0, MPICOM)
    CALL MPIBCAST (srf_emis_type,      LEN(srf_emis_type),               MPICHAR, 0, MPICOM)
    CALL MPIBCAST (srf_emis_cycle_yr,  1,                                MPIINT,  0, MPICOM)
    CALL MPIBCAST (srf_emis_fixed_ymd, 1,                                MPIINT,  0, MPICOM)
    CALL MPIBCAST (srf_emis_fixed_tod, 1,                                MPIINT,  0, MPICOM)
    CALL MPIBCAST (ext_frc_specifier,  LEN(ext_frc_specifier(1))*pcnst,  MPICHAR, 0, MPICOM)
    CALL MPIBCAST (ext_frc_type,       LEN(ext_frc_type),                MPICHAR, 0, MPICOM)
    CALL MPIBCAST (ext_frc_cycle_yr,   1,                                MPIINT,  0, MPICOM)
    CALL MPIBCAST (ext_frc_fixed_ymd,  1,                                MPIINT,  0, MPICOM)
    CALL MPIBCAST (ext_frc_fixed_tod,  1,                                MPIINT,  0, MPICOM)
#endif

    ! Update "short_lived_species" arrays - will eventually unify these
    nSlvd = nSls
    ALLOCATE(slvd_Lst(nSlvd), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating slvd_Lst')
    ALLOCATE(slvd_ref_MMR(nSlvd), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating slvd_ref_MMR')
    DO I = 1, nSls
       slvd_Lst(I) = TRIM(slsNames(I))
    ENDDO

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------

    chem_is_active = .true.

  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    ! Author: B. Eaton
    !
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------Arguments---------------------------------

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

!===============================================================================

  subroutine chem_init(phys_state, pbuf2d)
    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize GEOS-Chem parts (state objects, mainly)
    !          (and declare history variables)
    !
    !-----------------------------------------------------------------------
    use physics_buffer,        only : physics_buffer_desc, pbuf_get_index
    use chem_mods,             only : map2GC_dryDep, drySpc_ndx

#ifdef SPMD
    use mpishorthand
#endif
    use cam_abortutils,        only : endrun
    use mo_chem_utls,          only : get_spc_ndx

    use Phys_Grid,             only : get_Area_All_p
    use hycoef,                only : ps0, hyai, hybi, hyam

    use seq_drydep_mod,        only : drydep_method, DD_XLND
    use gas_wetdep_opts,       only : gas_wetdep_method
    use mo_neu_wetdep,         only : neu_wetdep_init

#if defined( MODAL_AERO_4MODE )
    use aero_model,            only : aero_model_init
    use mo_setsox,             only : sox_inti
    use mo_drydep,             only : drydep_inti_landuse
    use modal_aero_data,       only : ntot_amode, nspec_amode
    use modal_aero_data,       only : xname_massptr
#endif

    use Input_Opt_Mod
    use State_Chm_Mod
    use State_Grid_Mod
    use State_Met_Mod
    use DiagList_Mod,          only : Init_DiagList, Print_DiagList
    use TaggedDiagList_Mod,    only : Init_TaggedDiagList, Print_TaggedDiagList
    use GC_Grid_Mod,           only : SetGridFromCtrEdges

    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,         only : PI, PI_180, Re

    use Time_Mod,              only : Accept_External_Date_Time
    use Linoz_Mod,             only : Linoz_Read

    use CMN_Size_Mod

    use Drydep_Mod,            only : depName, Ndvzind
    use Pressure_Mod,          only : Accept_External_ApBp
    use Chemistry_Mod,         only : Init_Chemistry
    use Ucx_Mod,               only : Init_Ucx
    use Strat_chem_Mod,        only : Init_Strat_Chem
    use isorropiaII_Mod,       only : Init_IsorropiaII
    use Input_mod,             only : Validate_Directories
    use Olson_Landmap_Mod
    use Vdiff_Mod

    use mo_setinv,             only : setinv_inti
    use mo_mean_mass,          only : init_mean_mass
    use tracer_cnst,           only : tracer_cnst_init
    use tracer_srcs,           only : tracer_srcs_init

    use CESMGC_Emissions_Mod,  only : CESMGC_Emissions_Init
    use CESMGC_Diag_Mod,       only : CESMGC_Diag_Init

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
    INTEGER                :: nStrat
    INTEGER                :: I, J, L, N, M
    INTEGER                :: RC
    INTEGER                :: nLinoz

    ! Logicals
    LOGICAL                :: prtDebug
    LOGICAL                :: Found

    ! Strings
    CHARACTER(LEN=255)     :: historyConfigFile
    CHARACTER(LEN=255)     :: SpcName

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

    ! LCHNK: which chunks we have on this process
    LCHNK = phys_state%LCHNK
    ! NCOL: number of atmospheric columns for each chunk
    NCOL  = phys_state%NCOL

    write(iulog,'(2(a,x,I6,x))') 'chem_init called on PE ', myCPU, ' of ', nCPUs

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

    ! Set some basic flags
    Input_Opt%LUCX      = .True.

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Validate_Directories( Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Validation_Directories"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

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

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt  = Input_Opt, &
                       State_Grid = maxGrid,   &
                       RC         = RC        )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Define more variables for maxGrid
    maxGrid%MaxTropLev  = nZ
    maxGrid%MaxStratLev = nStrat
    IF ( Input_Opt%LUCX ) THEN
       maxGrid%MaxChemLev = maxGrid%MaxStratLev
    ELSE
       maxGrid%MaxChemLev = maxGrid%MaxTropLev
    ENDIF


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
       State_Grid(I)%MaxTropLev  = nZ
       State_Grid(I)%MaxStratLev = nStrat

       ! Set maximum number of levels in the chemistry grid
       IF ( Input_Opt%LUCX ) THEN
          State_Grid(I)%MaxChemLev = State_Grid(I)%MaxStratLev
       ELSE
          State_Grid(I)%MaxChemLev = State_Grid(I)%MaxTropLev
       ENDIF

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

    Input_Opt%thisCPU  = myCPU
    Input_Opt%amIRoot  = MasterProc

    ! TODO: Mimic GEOS-Chem's reading of input options
    !IF (MasterProc) THEN
    !   CALL Read_Input_File( Input_Opt   = Input_Opt,     &
    !                         srcFile     = inputGeosPath, &
    !                         RC          = RC )
    !ENDIF
    !CALL <broadcast data to other CPUs>

    Input_Opt%DryRun               = .False.

    ! For now just hard-code it
    ! First setup directories
    Input_Opt%Chem_Inputs_Dir      = TRIM(chemInputsDir)

    Input_Opt%SpcDatabaseFile      = TRIM(speciesDBPath)

    ! Simulation menu
    Input_Opt%NYMDb                = 20000101
    Input_Opt%NHMSb                =   000000
    Input_Opt%NYMDe                = 20010101
    Input_Opt%NHMSe                =   000000

    ! Now READ_SIMULATION_MENU
    Input_Opt%ITS_A_CH4_SIM          = .False.
    Input_Opt%ITS_A_CO2_SIM          = .False.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .True.
    Input_Opt%ITS_A_MERCURY_SIM      = .False.
    Input_Opt%ITS_A_POPS_SIM         = .False.
    Input_Opt%ITS_A_RnPbBe_SIM       = .False.
    Input_Opt%ITS_A_TAGO3_SIM        = .False.
    Input_Opt%ITS_A_TAGCO_SIM        = .False.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .False.

    ! Now READ_ADVECTED_SPECIES_MENU
    Input_Opt%N_Advect               = nTracers
    IF (Input_Opt%N_Advect.GT.Input_Opt%Max_AdvectSpc) THEN
       CALL ENDRUN('Number of tracers exceeds max count')
    ENDIF
    ! Assign tracer names
    DO J = 1, Input_Opt%N_Advect
       Input_Opt%AdvectSpc_Name(J) = TRIM(tracerNames(J))
    ENDDO
    ! No tagged species
    Input_Opt%LSplit = .False.

    ! Now READ_TRANSPORT_MENU
    Input_Opt%LTran                  = .True.
    Input_Opt%LFill                  = .True.
    Input_Opt%TPCore_IOrd            = 3
    Input_Opt%TPCore_JOrd            = 3
    Input_Opt%TPCore_KOrd            = 3

    ! Now READ_PHOTOLYSIS_MENU
    Input_Opt%FAST_JX_DIR            ='/glade/p/univ/umit0034/ExtData/' // &
     'CHEM_INPUTS/FAST_JX/v2020-02/'

    ! Now READ_CONVECTION_MENU
    Input_Opt%LConv                  = .False.
    Input_Opt%LTurb                  = .True.
    Input_Opt%LNLPBL                 = .True.

    ! Now READ_EMISSIONS_MENU
    ! This menu is pointless in CESM-GC
    Input_Opt%LEmis                  = .False.
    Input_Opt%HCOConfigFile          = 'HEMCO_Config.rc'

    Input_Opt%LSoilNOx               = .True.

    ! Set surface VMRs - turn this off so that CAM can handle it
    Input_Opt%LCH4Emis               = .False.
    Input_Opt%LCH4SBC                = .False.

    ! Set initial conditions
    Input_Opt%LSetH2O                = .False. !TMMF

    ! Now READ_AEROSOL_MENU
    Input_Opt%LSulf                  = .True.
    Input_Opt%LMetalcatSO2           = .True.
    Input_Opt%LCarb                  = .True.
    Input_Opt%LBrC                   = .False.
    Input_Opt%LSOA                   = .False.
    Input_Opt%LSVPOA                 = .False.
    Input_Opt%LDust                  = .True.
    Input_Opt%LDstUp                 = .False.
    Input_Opt%LSSalt                 = .True.
    Input_Opt%SalA_rEdge_um(1)       = 0.01e+0_fp
    Input_Opt%SalA_rEdge_um(2)       = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(1)       = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(2)       = 8.00e+0_fp
    Input_Opt%LMPOA                  = .False.
    ! For now, disable solid PSCs and strat aerosol settling
    ! Our treatment of the stratosphere isn't really sophisticated
    ! enough to warrant it yet
    Input_Opt%LGravStrat             = .False.
    Input_Opt%LSolidPSC              = .False.
    Input_Opt%LHomNucNAT             = .False.
    Input_Opt%T_NAT_Supercool        = 3.0e+0_fp
    Input_Opt%P_Ice_Supersat         = 1.2e+0_fp
    Input_Opt%LPSCChem               = .True.
    Input_Opt%LStratOD               = .True.

    Input_Opt%LBCAE                  = .True.
    Input_Opt%BCAE_1                 = 1.5e+0_fp
    Input_Opt%BCAE_2                 = 1.0e+0_fp
    Input_Opt%hvAerNIT               = .False.
    Input_Opt%hvAerNIT_JNIT          = 0.0e+00_fp
    Input_Opt%hvAerNIT_JNITs         = 0.0e+00_fp
    Input_Opt%JNITChanA              = 66.667e+0_fp
    Input_Opt%JNITChanB              = 33.333e+0_fp

    ! Now READ_DEPOSITION_MENU
    Input_Opt%LDryD                  = .True.
    !==================================================================
    ! Add the following options:
    ! + GEOS-Chem computes ALL dry-deposition velocities
    ! + CLM computes land velocities. Velocities over ocean and ice are
    !   computed in a MOZART-like way
    ! + CLM computes land velocities. Velocities over ocean and ice are
    !   computed from GEOS-Chem
    !
    ! Note: What to do about aerosols? Who should compute the dry
    !       deposition velocities
    !
    ! Thibaud M. Fritz - 26 Feb 2020
    !==================================================================
    Input_Opt%LWetD                  = .True.
    Input_Opt%CO2_Effect             = .False.
    Input_Opt%CO2_Level              = 600.0_fp
    Input_Opt%CO2_Ref                = 390.0_fp

    ! Now READ_CHEMISTRY_MENU
    Input_Opt%LChem                  = .True.
    Input_Opt%LSChem                 = .False. ! .True. !TMMF
    Input_Opt%LLinoz                 = .True.
    Input_Opt%LSynoz                 = .True.
    Input_Opt%LUCX                   = .True.
    Input_Opt%LActiveH2O             = .True.
    Input_Opt%Use_Online_O3          = .True.
    ! Expect to get total overhead ozone, although it should not
    ! make too much of a difference since we want to use "full-UCX"
    Input_Opt%Use_O3_from_Met        = .True.
    Input_Opt%Use_TOMS_O3            = .False.
    Input_Opt%Gamma_HO2              = 0.2e+0_fp

    Input_Opt%LPRT                   = .False.

    !==================================================================
    ! CESM-specific input flags
    !==================================================================

    ! onlineAlbedo    -> True  (use CLM albedo)
    !                 -> False (read monthly-mean albedo from HEMCO)
    Input_Opt%onlineAlbedo           = .True.

    ! onlineLandTypes -> True  (use CLM landtypes)
    !                 -> False (read landtypes from HEMCO)
    Input_Opt%onlineLandTypes        = .True.

    ! ddVel_CLM       -> True  (use CLM dry deposition velocities)
    !                 -> False (let GEOS-Chem compute dry deposition velocities)
    Input_Opt%ddVel_CLM              = .False.

    ! applyQtend: apply tendencies of water vapor to specific humidity
    Input_Opt%applyQtend             = .False.

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
       IF (IERR.NE.0) CALL ENDRUN('Failure while allocating linozData')
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
#if defined( SPMD )
       CALL MPIBCAST( linozData, nLinoz, MPIR8, 0, MPICOM )
#endif
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
    ! This requires input.geos and HISTORY.rc to be in the run directory
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

#if defined ( MODAL_AERO )
    IF ( Input_Opt%LWetD ) THEN
       DO I = BEGCHUNK, ENDCHUNK
          DO N = 1, State_Chm(I)%nWetDep
             M = State_Chm(I)%Map_WetDep(N)
             SpcInfo => State_Chm(I)%SpcData(M)%Info
             SELECT CASE ( TRIM(SpcInfo%Name) )
             CASE ( 'BCPI', 'BCPO', 'DST1', 'DST2', 'DST3', 'DST4', &
                    'SOAS', 'SO4', 'SALA', 'SALC', 'OCPI' , 'OCPO' )
                SpcInfo%WD_ExternalDep = .True.
             CASE DEFAULT
                SpcInfo%WD_ExternalDep = .False.
             END SELECT
             SpcInfo => NULL()
          ENDDO
       ENDDO
    ENDIF
#endif

    IF ( Input_Opt%LDryD ) THEN
       !==============================================================
       ! Get mapping between CESM dry deposited species and the
       ! indices of State_Chm%DryDepVel. This needs to be done after
       ! Init_Drydep
       ! Thibaud M. Fritz - 04 Mar 2020
       !==============================================================

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

#if defined( MODAL_AERO_4MODE )
    ! Initialize aqueous chem
    CALL SOx_inti()

    ! Initialize aerosols
    CALL aero_model_init( pbuf2d )

    ! Initialize land maps for aerosol dry deposition
    IF ( drydep_method == DD_XLND ) THEN
       CALL drydep_inti_landuse( depvel_lnd_file, &
                                 clim_soilw_file )
    ELSE
       Write(iulog,'(a,a)') ' drydep_method is set to: ', TRIM(drydep_method)
       CALL ENDRUN('drydep_method must be DD_XLND to compute land maps for aerosol' // &
               ' dry deposition!')
    ENDIF
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

    IF ( Input_Opt%LChem .AND. &
         Input_Opt%LUCX ) THEN
       CALL Init_UCX( Input_Opt  = Input_Opt,            &
                      State_Chm  = State_Chm(BEGCHUNK),  &
                      State_Diag = State_Diag(BEGCHUNK), &
                      State_Grid = maxGrid              )
    ENDIF

    IF ( Input_Opt%LSCHEM ) THEN
        CALL Init_Strat_Chem( Input_Opt  = Input_Opt,           &
                              State_Chm  = State_Chm(BEGCHUNK), &
                              State_Met  = State_Met(BEGCHUNK), &
                              State_Grid = maxGrid,             &
                              RC         = RC                  )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Init_Strat_Chem"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF ( Input_Opt%LSSalt ) THEN
       CALL INIT_ISORROPIAII( State_Grid = maxGrid )
    ENDIF

    ! Get some indices
    iH2O  = Ind_('H2O')
    iOH   = Ind_('OH')
    iO3   = Ind_('O3')
    iCH4  = Ind_('CH4')
    iCO   = Ind_('CO')
    iNO   = Ind_('NO')
    ! The following indices are needed to compute invariants
    iO    = Ind_('O')
    iH    = Ind_('H')
    iO2   = Ind_('O2')

    ! This is used to compute gas-phase H2SO4 production
    SpcInfo => State_Chm(BEGCHUNK)%SpcData(iOH)%Info
    MWOH    = REAL(SpcInfo%MW_g,r8)
    ! Free pointer
    SpcInfo => NULL()

    ! This is used to compute gas-phase H2SO4 production
    iPSO4   = Ind_('PSO4')
    SpcInfo => State_Chm(BEGCHUNK)%SpcData(iPSO4)%Info
    MWPSO4  = REAL(SpcInfo%MW_g,r8)
    ! Free pointer
    SpcInfo => NULL()

    ! This is used to compute overhead ozone column
    SpcInfo => State_Chm(BEGCHUNK)%SpcData(iO3)%Info
    MWO3    = REAL(SpcInfo%MW_g,r8)
    ! Free pointer
    SpcInfo => NULL()

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

    ! Initialize diagnostics interface
    CALL CESMGC_Diag_Init( Input_Opt = Input_Opt,           &
                           State_Chm = State_Chm(BEGCHUNK), &
                           State_Met = State_Met(BEGCHUNK) )

    ! Initialize emissions interface
    CALL CESMGC_Emissions_Init( lght_no_prd_factor = lght_no_prd_factor )

    hco_pbuf2d => pbuf2d

    If ( MasterProc ) Write(iulog,*) "hco_pbuf2d now points to pbuf2d"

    ! Cleanup
    Call Cleanup_State_Grid( maxGrid, RC )

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer,   only: physics_buffer_desc

    TYPE(physics_state), INTENT(IN):: phys_state(begchunk:endchunk)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    ! Not sure what we would realistically do here rather than in tend

  end subroutine chem_timestep_init

!===============================================================================

  subroutine GC_Update_Timesteps(DT)

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
        IF (MasterProc) WRITE(iulog,'(a,F7.1,a)') ' --> GC: updating dt to ', DT, ' seconds'

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

  end subroutine

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dT, pbuf, fh2o )

    use physics_buffer,      only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use physics_buffer,      only : pbuf_get_chunk, pbuf_get_index
    use cam_history,         only : outfld, hist_fld_active
    use camsrfexch,          only : cam_in_t, cam_out_t

#ifdef SPMD
    use mpishorthand
#endif

    use phys_grid,           only : get_ncols_p, get_rlat_all_p, get_rlon_all_p

    use chem_mods,           only : drySpc_ndx, map2GC_dryDep
    use chem_mods,           only : nfs, indexm, gas_pcnst
    use mo_mean_mass,        only : set_mean_mass
    use mo_setinv,           only : setinv
    use mo_neu_wetdep,       only : neu_wetdep_tend
    use gas_wetdep_opts,     only : gas_wetdep_method
#if defined( MODAL_AERO_4MODE )
    use modal_aero_data,     only : ntot_amode, nspec_amode
    use modal_aero_data,     only : lmassptr_amode
    use modal_aero_data,     only : xname_massptr
#endif

    use Olson_Landmap_Mod,   only : Compute_Olson_Landmap
    use Modis_LAI_Mod,       only : Compute_XLAI
    use CMN_Size_Mod,        only : NSURFTYPE
    use Drydep_Mod,          only : Do_Drydep
    use Drydep_Mod,          only : DEPNAME, NDVZIND
    use Drydep_Mod,          only : Update_DryDepFreq

    use Calc_Met_Mod,        only : Set_Dry_Surface_Pressure
    use Calc_Met_Mod,        only : AirQnt
    use GC_Grid_Mod,         only : SetGridFromCtr
    use Pressure_Mod,        only : Set_Floating_Pressures
    use Pressure_Mod,        only : Accept_External_Pedge
    use Time_Mod,            only : Accept_External_Date_Time
    use Toms_Mod,            only : Compute_Overhead_O3
    use Chemistry_Mod,       only : Do_Chemistry
    use Wetscav_Mod,         only : Setup_Wetscav, Do_WetDep
    use CMN_Size_Mod,        only : PTop
    use PBL_Mix_Mod,         only : Compute_PBL_Height
    use UCX_Mod,             only : Set_H2O_Trac
    use CMN_FJX_MOD,         only : ZPJ
    use State_Diag_Mod,      only : get_TagInfo
    use Unitconv_Mod,        only : Convert_Spc_Units

    use CESMGC_Emissions_Mod,only : CESMGC_Emissions_Calc
    use CESMGC_Diag_Mod,     only : CESMGC_Diag_Calc
    use CESMGC_Diag_Mod,     only : wetdep_name, wtrate_name

    use Tropopause,          only : Tropopause_findChemTrop, Tropopause_Find
    use HCO_Utilities_GC_Mod  ! Utility routines for GC-HEMCO interface

    ! For calculating SZA
    use Orbit,               only : zenith
    use Time_Manager,        only : Get_Curr_Calday, Get_Curr_Date

    ! Calculating relative humidity
    use WV_Saturation,       only : QSat

    ! Grid area
    use Phys_Grid,           only : get_area_all_p, get_lat_all_p, get_lon_all_p

    use short_lived_species, only : get_short_lived_species
    use short_lived_species, only : set_short_lived_species

#if defined( MODAL_AERO )
    ! Aqueous chemistry and aerosol growth
    use aero_model,          only : aero_model_gasaerexch
#endif

    use rad_constituents,    only : rad_cnst_get_info

    ! GEOS-Chem version of physical constants
    use PhysConstants,       only : PI, PI_180, g0, AVO, Re, g0_100
    ! CAM version of physical constants
    use PhysConst,           only : MWDry, Gravit

    REAL(r8),            INTENT(IN)    :: dT          ! Time step
    TYPE(physics_state), INTENT(IN)    :: state       ! Physics State variables
    TYPE(physics_ptend), INTENT(OUT)   :: ptend       ! indivdual parameterization tendencies
    TYPE(cam_in_t),      INTENT(INOUT) :: cam_in
    TYPE(cam_out_t),     INTENT(IN)    :: cam_out
    TYPE(physics_buffer_desc), POINTER :: pbuf(:)
    REAL(r8), OPTIONAL,  INTENT(OUT)   :: fh2o(PCOLS) ! h2o flux to balance source from chemistry

    ! Initial MMR for all species
    REAL(r8) :: MMR_Beg(PCOLS,PVER,nSls+nTracers)
    REAL(r8) :: MMR_End(PCOLS,PVER,nSls+nTracers)

    ! Logical to apply tendencies to mixing ratios
    LOGICAL :: lq(pcnst)

    ! Indexing
    INTEGER :: N, M, P, SM, ND
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
    REAL(r8)          :: wetdepflx(pcols,pcnst)       ! Wet deposition fluxes (kg/m2/s)

#if defined( MODAL_AERO )
    REAL(r8)          :: binRatio(MAXVAL(nspec_amode(:)),ntot_amode,state%NCOL,PVER)
#endif

    ! For emissions
    REAL(r8)          :: eflx(pcols,pver,pcnst)       ! 3-D emissions in kg/m2/s

    ! For GEOS-Chem diagnostics
    REAL(r8)              :: mmr1(state%NCOL,PVER,gas_pcnst)
    REAL(r8)              :: mmr_tend(state%NCOL,PVER,gas_pcnst)
    REAL(r8)              :: wk_out(state%NCOL)
    LOGICAL               :: isWD
    LOGICAL               :: Found
    CHARACTER(LEN=255)    :: tagName

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
    REAL(r8)     :: cld(PCOLS,PVER)
    REAL(r8)     :: TauCli(PCOLS,PVER)
    REAL(r8)     :: TauClw(PCOLS,PVER)
    REAL(r8), PARAMETER :: re_m   = 1.0e-05_r8 ! Cloud drop radius in m
    REAL(r8), PARAMETER :: cldMin = 1.0e-02_r8 ! Minimum cloud cover
    REAL(r8), PARAMETER :: cnst   = 1.5e+00_r8 / (re_m * 1.0e+03_r8 * g0)

    ! Calculating SZA
    REAL(r8)      :: Calday

    CHARACTER(LEN=255)     :: SpcName
    INTEGER                :: SpcId
    TYPE(Species), POINTER :: SpcInfo

    CHARACTER(LEN=63)      :: OrigUnit

    REAL(r8)           :: SlsData(PCOLS, PVER, nSls)

    INTEGER            :: currYr, currMo, currDy, currTOD
    INTEGER            :: currYMD, currHMS, currHr, currMn, currSc
    REAL(f4)           :: currUTC
    LOGICAL            :: firstDay = .True.
    LOGICAL            :: newDay   = .False.
    LOGICAL            :: newMonth = .False.

    TYPE(physics_buffer_desc), POINTER :: pbuf_chnk(:) ! slice of pbuf in chnk
    REAL(r8), POINTER     :: pbuf_ik(:,:)          ! ptr to pbuf data (/pcols,pver/)
    INTEGER               :: tmpIdx                ! pbuf field id
    CHARACTER(LEN=255)    :: fldname_ns            ! field name

    INTEGER            :: TIM_NDX

    INTEGER, SAVE      :: iStep = 0
    LOGICAL            :: rootChunk
    LOGICAL            :: lastChunk
    INTEGER            :: RC


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

    ! LCHNK: which chunk we have on this process
    LCHNK = state%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    NCOL  = state%NCOL

    ! Root Chunk
    rootChunk = ( MasterProc .and. (LCHNK==BEGCHUNK) )
    ! Last Chunk
    lastChunk = ( MasterProc .and. (LCHNK==ENDCHUNK) )

    ! Count the number of steps which have passed
    IF (LCHNK.EQ.BEGCHUNK) iStep = iStep + 1

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
    DO J = 1, nY
       State_Grid(LCHNK)%Area_M2(1,J) = REAL(Col_Area(J) * Re**2,fp)
       State_Met(LCHNK)%Area_M2(1,J)  = State_Grid(LCHNK)%Area_M2(1,J)
    ENDDO

    ! 2. Copy tracers into State_Chm
    ! Data was received in kg/kg dry
    State_Chm(LCHNK)%Spc_Units = 'kg/kg dry'
    ! Initialize ALL State_Chm species data to zero, not just tracers
    State_Chm(LCHNK)%Species = 0.0e+0_fp

    lq(:) = .False.

    ! Map and flip gaseous species
    MMR_Beg = 0.0e+0_r8
    MMR_End = 0.0e+0_r8
    DO N = 1, pcnst
       M = map2GC(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             MMR_Beg(J,L,M) = state%q(J,nZ+1-L,N)
             State_Chm(LCHNK)%Species(1,J,L,M) = REAL(MMR_Beg(J,L,M),fp)
          ENDDO
          ENDDO
          lq(N) = .True.
       ENDIF
    ENDDO

    ! We need to let CAM know that 'H2O' and 'Q' are identical
    DO J = 1, nY
    DO L = 1, nZ
       MMR_Beg(J,L,iH2O) = state%q(J,nZ+1-L,cQ)
       State_Chm(LCHNK)%Species(1,J,L,iH2O) = REAL(MMR_Beg(J,L,iH2O),fp)
    ENDDO
    ENDDO

    ! Retrieve previous value of species data
    SlsData(:,:,:) = 0.0e+0_r8
    CALL get_short_lived_species( SlsData, LCHNK, nY, pbuf )

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
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             State_Chm(LCHNK)%Species(1,J,L,M) = REAL(SlsData(J,nZ+1-L,N),fp)
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          vmr0(:nY,:nZ,N) = State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,M) * &
                            MWDry / adv_mass(N)
          ! We'll substract concentrations after chemistry later
          mmr_tend(:nY,:nZ,N) = REAL(State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,M),r8)
       ELSEIF ( M < 0 ) THEN
          vmr0(:nY,:nZ,N) = state%q(:nY,:nZ,-M) * &
                            MWDry / adv_mass(N)
          mmr_tend(:nY,:nZ,N) = state%q(:nY,:nZ,-M)
       ENDIF
    ENDDO

#if defined( MODAL_AERO_4MODE )
    ! First reset State_Chm%Species to zero for aerosols
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M)
          IF ( P > 0 ) State_Chm(LCHNK)%Species(1,:nY,:nZ,P) = 0.0e+00_fp
       ENDDO
    ENDDO

    ! Map and flip aerosols
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          ! TMMF - Should there be a ratio of molar weights involved?
          P = map2MAM4(SM,M)
          IF ( P <= 0 ) CYCLE
          N = lmassptr_amode(SM,M)
          ! Multiple MAM4 bins are mapped to same GEOS-Chem species
          State_Chm(LCHNK)%Species(1,:nY,:nZ,P) = State_Chm(LCHNK)%Species(1,:nY,:nZ,P) &
                                                + REAL(state%q(:nY,nZ:1:-1,N),fp)
       ENDDO
    ENDDO
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M)
          IF ( P <= 0 ) CYCLE
          N = lmassptr_amode(SM,M)
          DO J = 1, nY
          DO L = 1, nZ
             IF ( State_Chm(LCHNK)%Species(1,J,nZ+1-L,P) > 0.0e+00_r8 ) THEN
                binRatio(SM,M,J,L) = REAL(state%q(J,L,N),r8) &
                   / State_Chm(LCHNK)%Species(1,J,nZ+1-L,P)
             ELSE
                binRatio(SM,M,J,L) = 0.0e+00_r8
             ENDIF
          ENDDO
          ENDDO
       ENDDO
    ENDDO
#endif

    ! If H2O tendencies are propagated to specific humidity, then make sure
    ! that Q actually applies tendencies
    IF ( Input_Opt%applyQtend ) lq(cQ) = .True.

    ! Initialize tendency array
    CALL Physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)

    ! Reset chemical tendencies
    ptend%q(:,:,:) = 0.0e+0_r8

    ! Determine current date and time
    CALL Get_Curr_Date( yr  = currYr,  &
                        mon = currMo,  &
                        day = currDy,  &
                        tod = currTOD )

    ! For now, force year to be 2000
    currYr  = 2000
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

    IF ( firstDay ) THEN
       newDay   = .True.
       newMonth = .True.
       firstDay = .False.
    ELSE IF ( currHMS < dT ) THEN
       newDay = .True.
       IF ( currDy == 1 ) THEN
          newMonth = .True.
       ELSE
          newMonth = .False.
       ENDIF
    ELSE
       newDay   = .False.
       newMonth = .False.
    ENDIF

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
       State_Chm(LCHNK)%Species(1,J,nZ+1-L,iH2O) = qH2O(J,L)
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
    CALL QSat(state%t(:nY,:), state%pmid(:nY,:), satV, satQ)
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

    ! Note: all using CAM vertical convention (1 = TOA)
    ! Calculation is based on that done for MOZART
    DO J = 1, nY
    DO L = nZ, 1, -1
       cldW(J,L) = state%q(J,L,ixCldLiq) + state%q(J,L,ixCldIce)
       ! Convert water mixing ratio [kg/kg] to water content [g/m^3]
       IF ( cldW(J,L) * state%pmid(J,L) / &
             (state%T(J,L) * 287.0e+00_r8) * 1.0e+03_r8 <= 0.01_r8 .AND. &
            cldFrc(J,L) /= 0.0e+00_r8 ) THEN
          cld(J,L) = 0.0e+00_r8
       ELSE
          cld(J,L) = cldFrc(J,L)
       ENDIF
       IF ( ixNDrop > 0 ) nCldWtr(J,L) = state%q(J,L,ixNDrop)
    ENDDO
    ENDDO

    DO J = 1, nY
       IF ( COUNT( cld(J,:nZ) > cldMin ) > 0 ) THEN
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
             !
             ! Unit check:                    |
             ! Q    : [kg H2O/kg air]         |
             ! Pint : [Pa]=[kg air/m/s^2]     |
             ! re   : [m]                     |   = 1.0e-5
             ! rho_w: [kg H2O/m^3]            |   = 1.0e+3
             ! g    : [m/s^2]                 |   = 9.81
             TauClw(J,L) = state%q(J,L,ixCldLiq)               &
                         * (state%pint(J,L+1)-state%pint(J,L)) &
                         * cnst
             TauClw(J,L) = MAX(TauClw(J,L), 0.0e+00_r8)
             TauCli(J,L) = state%q(J,L,ixCldIce)               &
                         * (state%pint(J,L+1)-state%pint(J,L)) &
                         * cnst
             TauCli(J,L) = MAX(TauCli(J,L), 0.0e+00_r8)
          ENDDO
       ENDIF
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
    IF ( Input_Opt%onlineLandTypes ) THEN
       ! Fill in water
       State_Met(LCHNK)%LandTypeFrac(1,:nY,1) = cam_in%ocnFrac(:nY)     &
                                              + cam_in%iceFrac(:nY)
       IF ( .NOT. Input_Opt%ddVel_CLM ) THEN
          CALL getLandTypes( cam_in,         &
                             nY,             &
                             State_Met(LCHNK) )
       ENDIF
    ELSE
       DO N = 1, NSURFTYPE
          Write(fldname_ns, '(a,i2.2)') 'HCO_LANDTYPE', N-1
          tmpIdx = pbuf_get_index(fldname_ns, rc)
          IF ( tmpIdx < 0 ) THEN
             ! there is an error here and the field was not found
             IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
          ELSE
             CALL pbuf_get_field(pbuf, tmpIdx, pbuf_ik)
             DO J = 1, nY
                State_Met(LCHNK)%LandTypeFrac(1,J,N) = pbuf_ik(J,nZ)
                ! 2-D data is stored in the 1st level of a
                ! 3-D array due to laziness
             ENDDO
             pbuf_ik   => NULL()
          ENDIF

          Write(fldname_ns, '(a,i2.2)') 'HCO_XLAI', N-1
          tmpIdx = pbuf_get_index(fldname_ns, rc)
          IF ( tmpIdx < 0 ) THEN
             ! there is an error here and the field was not found
             IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
          ELSE
             CALL pbuf_get_field(pbuf, tmpIdx, pbuf_ik)
             DO J = 1, nY
                State_Met(LCHNK)%XLAI_NATIVE(1,J,N) = pbuf_ik(J,nZ)
                ! 2-D data is stored in the 1st level of a
                ! 3-D array due to laziness
             ENDDO
             pbuf_ik   => NULL()
          ENDIF
       ENDDO
    ENDIF

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
    IF ( Input_Opt%onlineLandTypes ) THEN
       State_Met(LCHNK)%FRLAKE    (1,:nY) = cam_in%lwtgcell(:,3) + &
                                          cam_in%lwtgcell(:,4)
       State_Met(LCHNK)%FRLANDIC  (1,:nY) = cam_in%lwtgcell(:,2)
       State_Met(LCHNK)%FRSNO     (1,:nY) = 0.0e+0_fp
    ELSE
       State_Met(LCHNK)%FRLAKE    (1,:nY) = 0.0e+0_fp
       State_Met(LCHNK)%FRLANDIC  (1,:nY) = 0.0e+0_fp
       State_Met(LCHNK)%FRSNO     (1,:nY) = 0.0e+0_fp
    ENDIF

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

    ! Field      : TS, TSKIN
    ! Description: Surface temperature, surface skin temperature
    ! Unit       : K
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TS        (1,:nY) = cam_in%TS(:nY)
    State_Met(LCHNK)%TSKIN     (1,:nY) = cam_in%TS(:nY)

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
       fldname_ns = 'HCO_UV_ALBEDO'
       tmpIdx = pbuf_get_index(fldname_ns, RC)
       IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
          State_Met(LCHNK)%UVALBEDO(1,:nY) = 0.0e+0_fp
       ELSE
          pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
          CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
          State_Met(LCHNK)%UVALBEDO(1,:nY) = pbuf_ik(:nY,nZ)
          pbuf_chnk => NULL()
          pbuf_ik   => NULL()
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
    fldname_ns = 'HCO_iodide'
    tmpIdx = pbuf_get_index(fldname_ns, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
       State_Chm(LCHNK)%IODIDE(1,:nY)   = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
       State_Chm(LCHNK)%IODIDE(1,:nY) = pbuf_ik(:nY,nZ)
       pbuf_chnk => NULL()
       pbuf_ik   => NULL()
    ENDIF

    ! Field      : SALINITY
    ! Description: Ocean salinity
    ! Unit       : PSU
    ! Dimensions : nX, nY
    ! Note       : Possibly get ocean salinity from POP?
    fldname_ns = 'HCO_salinity'
    tmpIdx = pbuf_get_index(fldname_ns, RC)
    IF ( tmpIdx < 0 .or. ( iStep == 1 ) ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
       State_Chm(LCHNK)%SALINITY(1,:nY) = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
       State_Chm(LCHNK)%SALINITY(1,:nY) = pbuf_ik(:nY,nZ)
       pbuf_chnk => NULL()
       pbuf_ik   => NULL()
    ENDIF

    ! Field      : OMOC
    ! Description: OM/OC ratio
    ! Unit       : -
    ! Dimensions : nX, nY
    IF      ( currMo == 12 .or. currMo == 1  .or. currMo == 2  ) THEN
       fldname_ns = 'HCO_OMOC_DJF'
    ELSE IF ( currMo == 3  .or. currMo == 4  .or. currMo == 5  ) THEN
       fldname_ns = 'HCO_OMOC_MAM'
    ELSE IF ( currMo == 6  .or. currMo == 7  .or. currMo == 8  ) THEN
       fldname_ns = 'HCO_OMOC_JJA'
    ELSE IF ( currMo == 9  .or. currMo == 10 .or. currMo == 11 ) THEN
       fldname_ns = 'HCO_OMOC_SON'
    ENDIF
    tmpIdx = pbuf_get_index(fldname_ns, rc)
    IF ( tmpIdx < 0 ) THEN
       ! there is an error here and the field was not found
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
    ELSE
       CALL pbuf_get_field(pbuf, tmpIdx, pbuf_ik)
       DO J = 1, nY
          State_Chm(LCHNK)%OMOC(1,J) = pbuf_ik(J,nZ)
          ! 2-D data is stored in the 1st level of a
          ! 3-D array due to laziness
       ENDDO
       pbuf_ik   => NULL()
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
    ! if LUCX=T and LSETH2O=F and LACTIVEH2O=T, update specific humidity
    ! in the stratosphere
    !
    ! NOTE: Specific humidity may change in SET_H2O_TRAC and
    ! therefore this routine may call AIRQNT again to update
    ! air quantities and tracer concentrations (ewl, 10/28/15)
    IF ( Input_Opt%Its_A_Fullchem_Sim .and. iH2O > 0 ) THEN
       CALL Set_H2O_Trac( SETSTRAT   = ( ( .not. Input_Opt%LUCX )  &
                                         .or. Input_Opt%LSETH2O ), &
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

       ! Field      : LWI
       ! Description: Land/water indices
       ! Unit       : -
       ! Dimensions : nX, nY
       State_Met(LCHNK)%LWI(1,J) = FLOAT( iMaxLoc(1) )

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

       State_Met(LCHNK)%isSnow(1,J) = ( State_Met(LCHNK)%FRSEAICE(1,J) > 0.0e+0_fp &
                                   .or. State_Met(LCHNK)%SNODP(1,J) > 0.01 )

    ENDDO

    ! Do this after AirQnt in order to use AIRDEN and BXHEIGHT
    DO J = 1, nY
       O3col(J) = 0.0e+0_fp
       DO L = 1, nZ
          O3col(J) = O3col(J) &
                   + State_Chm(LCHNK)%Species(1,J,L,iO3) &
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

    !----------------------------------------------------------
    ! %%% GET SOME NON-EMISSIONS DATA FIELDS VIA HEMCO %%%
    !
    ! HEMCO can track non-emission data fields for chemistry
    ! simulations.  Put these subroutine calls after the
    ! call to EMISSIONS_RUN, so that the HEMCO data structure
    ! will be initialized. (bmy, 3/20/15)
    !
    ! HEMCO data list is now updated further above, so can
    ! take these calls out of the emissions sequence.
    ! (ckeller, 4/01/15)
    !----------------------------------------------------------
    !IF ( Input_Opt%LCHEM .and. newMonth ) THEN
    !
    !   ! The following only apply when photolysis is used,
    !   ! that is for fullchem or aerosol simulations.
    !   IF ( Input_Opt%Its_A_Fullchem_Sim  .or. Input_Opt%Its_An_Aerosol_Sim ) THEN
    !
    !      IF ( Input_Opt%USE_TOMS_O3 ) THEN
    !         ! Get TOMS overhead O3 columns for photolysis from
    !         ! the HEMCO data structure (bmy, 3/20/15)
    !         CALL Read_TOMS( Input_Opt = Input_Opt,  &
    !                         RC        = RC         )
    !
    !         ! Trap potential errors
    !         IF ( RC /= GC_SUCCESS ) THEN
    !            ErrMsg = 'Error encountered in "Read_TOMS"!'
    !            CALL Error_Stop( ErrMsg, ThisLoc )
    !         ENDIF
    !      ENDIF
    !
    !   ENDIF
    !
    !   ! Read data required for Hg2 gas-particle partitioning
    !   ! (H Amos, 25 Oct 2011)
    !   IF ( ITS_A_MERCURY_SIM ) THEN
    !      CALL Read_Hg2_Partitioning( Input_Opt  = Input_Opt,         &
    !                                  State_Grid = State_Grid(LCHNK), &
    !                                  State_Met  = State_Met(LCHNK),  &
    !                                  MONTH      = 1,                 & !TMMF
    !                                  RC         = RC                )
    !
    !      ! Trap potential errors
    !      IF ( RC /= GC_SUCCESS ) THEN
    !         ErrMsg =
    !            'Error encountered in "Read_Hg2_Partitioning"!'
    !         CALL Error_Stop( ErrMsg, ThisLoc )
    !      ENDIF
    !
    !   ENDIF
    !ENDIF

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
    ! CLM computes dry deposition velocities over land.
    ! We need to merge the land component passed through cam_in and
    ! the ocn/ice dry deposition velocities.
    !
    ! If using the CLM velocities, then use GEOS-Chem's dry deposition
    ! module to compute velocities and then scale them with the ocean
    ! fraction (Input_Opt%ddVel_CLM)
    !
    ! A second option would be to let GEOS-Chem compute dry deposition
    ! velocity, thus overwriting the input from CLM
    !
    ! drydep_method must be set to DD_XLND.
    !
    ! The GEOS-Chem option (.not. Input_Opt%ddVel_CLM) option coupled
    ! with Input_Opt%onlineLandTypes requires that CLM passes land
    ! type information (land type and leaf area index).
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

       IF ( Input_Opt%ddVel_CLM ) THEN
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
       ENDIF

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

    CALL CESMGC_Emissions_Calc( state      = state,            &
                                hco_pbuf2d = hco_pbuf2d,       &
                                State_Met  = State_Met(LCHNK), &
                                cam_in     = cam_in,           &
                                eflx       = eflx,             &
                                iStep      = iStep            )

    !-----------------------------------------------------------------------
    ! Add dry deposition flux 
    ! (stored as SurfaceFlux = -dflx)
    !-----------------------------------------------------------------------

    DO ND = 1, State_Chm(BEGCHUNK)%nDryDep
       ! Get the species ID from the drydep ID
       N = State_Chm(BEGCHUNK)%Map_DryDep(ND)
       IF ( N <= 0 ) CYCLE

       M = map2GCinv(N)
       IF ( M <= 0 ) CYCLE

       cam_in%cflx(1:nY,M) = cam_in%cflx(1:nY,M) &
                           + State_Chm(LCHNK)%SurfaceFlux(1,1:nY,N)
    ENDDO

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
          DO J = 1, nY
          DO L = 1, nZ
             State_Chm(LCHNK)%Species(1,J,L,M) = State_Chm(LCHNK)%Species(1,J,L,M) &
                                               + eflx(J,nZ+1-L,N) * dT
          ENDDO
          ENDDO
       ELSE
          ! Add to constituent (mostly for MAM4 aerosols)
          ! Convert from kg/m2/s to kg/kg/s
          DO J = 1, nY
          DO L = 1, nZ
             ptend%q(J,nZ+1-L,N) = ptend%q(J,nZ+1-L,N) &
                                 + eflx(J,nZ+1-L,N)    &
                                   / ( g0_100 * State_Met(LCHNK)%DELP_DRY(1,J,L) )
          ENDDO
          ENDDO
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
       CALL Set_H2O_Trac( SETSTRAT   = (.not. Input_Opt%LUCX), &
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

    ZPJ = 0.0e+0_r8
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

    !==============================================================
    ! ***** W E T   D E P O S I T I O N  (rainout + washout) *****
    !==============================================================
    IF ( Input_Opt%LWetD ) THEN

       ! Store mass mixing ratios before wet deposition is applied
       IF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
          DO N = 1, gas_pcnst
             isWD = .False.
             IF ( M > 0 ) THEN
                SpcInfo => State_Chm(BEGCHUNK)%SpcData(M)%Info
                isWD = SpcInfo%Is_WetDep

                ! Free pointer
                SpcInfo => NULL()
             ENDIF

             IF ( .NOT. isWD ) CYCLE

             IF ( hist_fld_active( TRIM(wetdep_name(N)) ) .OR. &
                  hist_fld_active( TRIM(wtrate_name(N)) ) ) THEN
                ! See definition of map2chm
                M = map2chm(N)
                IF ( M > 0 ) THEN
                   mmr1(:nY,:nZ,N) = State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,M)
                ENDIF
             ENDIF
          ENDDO
       ENDIF

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
       ELSEIF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
          ! TMMF, If we perform GEOS-Chem washout and rainout, we should turn
          ! it off for MAM4 aerosols, as this is done in physpkg.F90
          ! A way to do this would be to reapply vmr of MAM4 aerosols before
          ! wetdep is applied

          ! Do wet deposition
          ! Tracer concentration units are converted locally
          ! to [kg/m2] in wet deposition to enable calculations
          ! along the column (ewl, 9/18/15)
          CALL Do_WetDep( Input_Opt  = Input_Opt,         &
                          State_Chm  = State_Chm(LCHNK),  &
                          State_Diag = State_Diag(LCHNK), &
                          State_Grid = State_Grid(LCHNK), &
                          State_Met  = State_Met(LCHNK),  &
                          RC         = RC                )

          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_WetDep"!'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF

       ELSE
          ErrMsg = 'Unknown gas_wetdep_method ' //TRIM(gas_wetdep_method)
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       IF ( gas_wetdep_method == 'GEOS-CHEM' ) THEN
          DO N = 1, gas_pcnst
             isWD = .False.
             IF ( M > 0 ) THEN
                SpcInfo => State_Chm(BEGCHUNK)%SpcData(M)%Info
                isWD = SpcInfo%Is_WetDep

                ! Free pointer
                SpcInfo => NULL()
             ENDIF

             IF ( .NOT. isWD ) CYCLE

             IF ( hist_fld_active( TRIM(wetdep_name(N)) ) .OR. &
                  hist_fld_active( TRIM(wtrate_name(N)) ) ) THEN
                ! See definition of map2chm
                M = map2chm(N)
                IF ( M > 0 ) THEN
                   mmr1(:nY,:nZ,N) = State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,M) - mmr1(:nY,:nZ,N)
                ENDIF

                SpcName = wetdep_name(N)
                IF ( hist_fld_active(TRIM(SpcName)) ) THEN 
                   CALL Outfld( TRIM(SpcName), mmr1(1:nY,:nZ,N)/dT, nY, LCHNK )
                ENDIF

                SpcName = wtrate_name(N)
                IF ( hist_fld_active(TRIM(SpcName)) ) THEN 
                   wk_out = 0._r8
                   DO L = 1, nZ
                      wk_out(1:nY) = wk_out(1:nY) &
                                   + mmr1(1:nY,L,N)/dT                    * &
                                       State_Met(LCHNK)%AD(1,1:nY,nZ+1-L) / &
                                       State_Met(LCHNK)%Area_M2(1,1:nY)
                   ENDDO
                   CALL Outfld( TRIM(SpcName), wk_out(1:nY), nY, LCHNK )
                ENDIF
             ENDIF

             ! GEOS-Chem does not currently store HEFF, but calculates it
             ! internally. Some potential work around would be to add a
             ! SpcInfo%Heff variable.
          ENDDO
       ENDIF
    ENDIF

    !==============================================================
    ! ***** M A M   G A S - A E R O S O L   E X C H A N G E *****
    !==============================================================

#if defined( MODAL_AERO )
    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             vmr1(J,L,N) = State_Chm(LCHNK)%Species(1,J,nZ+1-L,M) * &
                MWDry / adv_mass(N)
          ENDDO
          ENDDO
       ELSEIF ( M < 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             vmr1(J,L,N) = state%q(J,L,-M) * &
                MWDry / adv_mass(N)
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    del_h2so4_gasprod = 0.0e+00_fp
    ! This needs to be in mol/mol over this timestep
    IF ( ( iPSO4 > 0 ) .and. ( MWPSO4 > 0.0e+00_fp ) ) THEN
       If ( rootChunk ) Write(iulog,*) " MAXVAL(PSO4) = ", &
          MAXVAL(State_Chm(LCHNK)%Species(1,:nY,:nZ,iPSO4))
       DO L = 1, nZ
          ! Convert from kg SO4/kg to mol/mol
          del_h2so4_gasprod(:nY,L) = &
          State_Chm(LCHNK)%Species(1,:nY,nZ+1-L,iPSO4) * MWDry / MWPSO4
       ENDDO
    ENDIF

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
#endif

    ! Make sure State_Chm(LCHNK) is back in kg/kg dry!
    IF ( TRIM(State_Chm(LCHNK)%Spc_Units) /= 'kg/kg dry' ) THEN
       Write(iulog,*) 'Current  unit = ', TRIM(State_Chm(LCHNK)%Spc_Units)
       Write(iulog,*) 'Expected unit = kg/ kg dry'
       CALL ENDRUN('Incorrect unit in GEOS-Chem State_Chm%Species')
    ENDIF

    ! Reset H2O MMR to the initial value (no chemistry tendency in H2O just yet)
    State_Chm(LCHNK)%Species(1,:,:,iH2O) = MMR_Beg(:,:,iH2O)

    ! Store unadvected species data
    SlsData = 0.0e+0_r8
    DO N = 1, nSls
       M = map2GC_Sls(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO L = 1, nZ
             SlsData(J,nZ+1-L,N) = REAL(State_Chm(LCHNK)%Species(1,J,L,M),r8)
          ENDDO
          ENDDO
       ENDIF
    ENDDO
    CALL set_short_lived_species( SlsData, LCHNK, nY, pbuf )

    DO N = 1, pcnst
       M = map2GC(N)
       IF ( M > 0 ) THEN
          ! Add change in mass mixing ratio to tendencies.
          ! For NEU wet deposition, the wet removal rates are added to
          ! ptend.
          DO J = 1, nY
          DO L = 1, nZ
             MMR_End (J,L,M) = REAL(State_Chm(LCHNK)%Species(1,J,L,M),r8)
             ptend%q(J,nZ+1-L,N) = ptend%q(J,nZ+1-L,N) &
                                 + (MMR_End(J,L,M)-MMR_Beg(J,L,M))/dT
          ENDDO
          ENDDO
       ENDIF
    ENDDO

#if defined( MODAL_AERO_4MODE )
    ! Here apply tendencies to MAM aerosols
    ! Initial mass in bin SM is stored as state%q(N)
    ! Final mass in bin SM is stored as binRatio(SM,M) * State_Chm(P)
    !
    ! We decide to apply chemical tendencies to all MAM aerosols,
    ! except so4, for which the chemically-produced sulfate gets
    ! partitioned in aero_model_gasaerexch
    DO M = 1, ntot_amode
       DO SM = 1, nspec_amode(M)
          P = map2MAM4(SM,M)
          IF ( P <= 0 .OR. to_upper(xname_massptr(SM,M)(:3)) == 'SO4' ) CYCLE
          N = lmassptr_amode(SM,M)
          ! Apply MAM4 chemical tendencies owing to GEOS-Chem aerosol processing
          ptend%q(:nY,:nZ,N) = ptend%q(:nY,:nZ,N) &
                             + (binRatio(SM,M,:nY,:nZ) * & 
                                REAL(State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,P),r8) &
                                - state%q(:nY,:nZ,N))/dT
       ENDDO
    ENDDO
#endif

    DO N = 1, gas_pcnst
       ! See definition of map2chm
       M = map2chm(N)
       IF ( M > 0 ) THEN
          mmr_tend(:nY,:nZ,N) = ( REAL(State_Chm(LCHNK)%Species(1,:nY,nZ:1:-1,M),r8) - mmr_tend(:nY,:nZ,N) ) / dT
       ELSEIF ( M < 0 ) THEN
          mmr_tend(:nY,:nZ,N) = ptend%q(:nY,:nZ,-M)
       ENDIF
    ENDDO

    IF ( Input_Opt%applyQtend ) THEN
       ! Apply GEOS-Chem's H2O mixing ratio tendency to CAM's specific humidity
       ! This requires to set lq(cQ) = lq(cH2O) ( = .True. )
       ptend%q(:,:,cQ) = ptend%q(:,:,cH2O)
    ENDIF

    CALL CESMGC_Diag_Calc( Input_Opt  = Input_Opt,         &
                           State_Chm  = State_Chm(LCHNK),  &
                           State_Diag = State_Diag(LCHNK), &
                           State_Grid = State_Grid(LCHNK), &
                           State_Met  = State_Met(LCHNK),  &
                           cam_in     = cam_in,            &
                           state      = state,             &
                           mmr_tend   = mmr_tend,          &
                           LCHNK      = LCHNK             )

    ! Debug statements
    ! Ozone tendencies
    IF ( rootChunk ) THEN
       Write(iulog,*) " MMR_Beg = ", MMR_Beg(1,:,iO3)
       Write(iulog,*) " MMR_End = ", MMR_End(1,:,iO3)
    ENDIF

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

    IF ( rootChunk ) WRITE(iulog,*) ' GEOS-Chem Chemistry step ', iStep, ' completed'
    IF ( lastChunk ) WRITE(iulog,*) ' Chemistry completed on all chunks completed of MasterProc'

  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    CHARACTER(LEN=*), INTENT(IN)  :: name       !  constituent name
    REAL(r8),         INTENT(IN)  :: latvals(:) ! lat in degrees (NCOL)
    REAL(r8),         INTENT(IN)  :: lonvals(:) ! lon in degrees (NCOL)
    LOGICAL,          INTENT(IN)  :: mask(:)    ! Only initialize where .true.
    REAL(r8),         INTENT(OUT) :: q(:,:)     ! kg tracer/kg dry air (NCOL, PVER)
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    INTEGER  :: iLev, NLEV, M
    REAL(r8) :: QTemp, Min_MMR

    NLEV = SIZE(q, 2)
    ! Retrieve a "background value" for this from the database
    Min_MMR = 1.0e-38_r8
    DO M = 1, gas_pcnst
       IF (TRIM(solsym(M)).eq.TRIM(name)) THEN
          Min_MMR = ref_MMR(M)
          EXIT
       ENDIF
    ENDDO

    DO iLev = 1, NLEV
       WHERE(mask)
          ! Set to the minimum mixing ratio
          q(:,iLev) = Min_MMR
       END WHERE
    ENDDO

  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final

    use Input_Opt_Mod,  only : Cleanup_Input_Opt
    use State_Chm_Mod,  only : Cleanup_State_Chm
    use State_Diag_Mod, only : Cleanup_State_Diag
    use State_Grid_Mod, only : Cleanup_State_Grid
    use State_Met_Mod,  only : Cleanup_State_Met
    use Error_Mod,      only : Cleanup_Error

    use FlexChem_Mod,   only : Cleanup_FlexChem
    use UCX_Mod,        only : Cleanup_UCX
    use Drydep_Mod,     only : Cleanup_Drydep
    use Carbon_Mod,     only : Cleanup_Carbon
    use Dust_Mod,       only : Cleanup_Dust
    use Seasalt_Mod,    only : Cleanup_Seasalt
    use Aerosol_Mod,    only : Cleanup_Aerosol
    use Sulfate_Mod,    only : Cleanup_Sulfate
    use Pressure_Mod,   only : Cleanup_Pressure
    use Strat_Chem_Mod, only : Cleanup_Strat_Chem

    use CMN_Size_Mod,   only : Cleanup_CMN_Size
    use CMN_FJX_Mod,    only : Cleanup_CMN_FJX

#ifdef BPCH_DIAG
    use CMN_O3_Mod,     only : Cleanup_CMN_O3
    ! Special: cleans up after NDXX_Setup
    use Diag_Mod,       only : Cleanup_Diag
#endif

    use CESMGC_Emissions_Mod, only: CESMGC_Emissions_Final

    ! Local variables
    INTEGER  :: I, RC

    ! Finalize GEOS-Chem

    CALL Cleanup_UCX
    CALL Cleanup_Aerosol
    CALL Cleanup_Carbon
    CALL Cleanup_Drydep
    CALL Cleanup_Dust
    CALL Cleanup_FlexChem( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_FlexChem"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_Pressure
    CALL Cleanup_Seasalt
    CALL Cleanup_Sulfate
    CALL Cleanup_Strat_Chem

    CALL CESMGC_Emissions_Final

    CALL Cleanup_CMN_SIZE( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_SIZE"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_CMN_FJX( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_FJX"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

#ifdef BPCH_DIAG
    CALL Cleanup_Diag

    ! Call extra cleanup routines, from modules in Headers/
    CALL Cleanup_CMN_O3( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_O3"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
#endif

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
    IF ( ALLOCATED( slvd_ref_MMR ) )  DEALLOCATE( slvd_ref_MMR )

    RETURN

  end subroutine chem_final
!===============================================================================
  subroutine chem_init_restart(File)
    use tracer_cnst,      only: init_tracer_cnst_restart
    use tracer_srcs,      only: init_tracer_srcs_restart
    use pio, only : file_desc_t

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT_RESTART'

    !
    ! data for offline tracers
    !
    call init_tracer_cnst_restart(File)
    call init_tracer_srcs_restart(File)
    !call init_linoz_data_restart(File)

  end subroutine chem_init_restart
!===============================================================================
  subroutine chem_write_restart( File )
    use tracer_cnst, only: write_tracer_cnst_restart
    use tracer_srcs, only: write_tracer_srcs_restart
    !use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    IF ( MasterProc ) WRITE(iulog,'(a)') 'GCCALL CHEM_WRITE_RESTART'
    !
    ! data for offline tracers
    !
    call write_tracer_cnst_restart(File)
    call write_tracer_srcs_restart(File)
    !call write_linoz_data_restart(File)
  end subroutine chem_write_restart
!===============================================================================
  subroutine chem_read_restart( File )
    use tracer_cnst, only: read_tracer_cnst_restart
    use tracer_srcs, only: read_tracer_srcs_restart
    !use linoz_data,  only: read_linoz_data_restart
    use pio, only : file_desc_t

    IMPLICIT NONE

    TYPE(file_desc_t) :: File

    IF ( MasterProc ) WRITE(iulog,'(a)') 'GCCALL CHEM_READ_RESTART'
    !
    ! data for offline tracers
    !
    call read_tracer_cnst_restart(File)
    call read_tracer_srcs_restart(File)
    !call read_linoz_data_restart(File)
  end subroutine chem_read_restart
!================================================================================
  subroutine chem_emissions( state, cam_in )

    use camsrfexch,          only : cam_in_t

    ! Arguments:

    TYPE(physics_state),    INTENT(IN)    :: state   ! Physics state variables
    TYPE(cam_in_t),         INTENT(INOUT) :: cam_in  ! import state

    INTEGER :: M, N
    INTEGER :: LCHNK, nY
    LOGICAL :: rootChunk


    ! LCHNK: which chunk we have on this process
    LCHNK = state%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    nY    = state%NCOL
    rootChunk = ( MasterProc.and.(LCHNK.EQ.BEGCHUNK) )

    !-----------------------------------------------------------------------
    ! Reset surface fluxes
    !-----------------------------------------------------------------------

    DO M = iFirstCnst, pcnst
       !N = map2chm(M)
       !IF ( N > 0 ) cam_in%cflx(1:nY,N) = 0.0e+0_r8
       cam_in%cflx(1:nY,M) = 0.0e+0_r8
    ENDDO

  end subroutine chem_emissions

end module chemistry
