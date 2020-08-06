!================================================================================================
! This is the "GEOS-Chem" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8, shr_kind_cl
  use physics_types,       only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,      only: physics_buffer_desc
  use ppgrid,              only: begchunk, endchunk, pcols
  use ppgrid,              only: pver, pverp
  use constituents,        only: pcnst, cnst_add, cnst_get_ind
  use shr_const_mod,       only: molw_dryair=>SHR_CONST_MWDAIR
  use seq_drydep_mod,      only : nddvels => n_drydep, drydep_list
  use spmd_utils,          only : MasterProc, myCPU=>Iam, nCPUs=>npes
  use cam_logfile,         only : iulog
  use string_utils,        only : to_upper

  !--------------------------------------------------------------------
  ! Basic GEOS-Chem modules
  !--------------------------------------------------------------------
  USE DiagList_Mod,        ONLY : DgnList    ! Derived type for diagnostics list
  USE Input_Opt_Mod,       ONLY : OptInput   ! Derived type for Input Options
  USE State_Chm_Mod,       ONLY : ChmState   ! Derived type for Chemistry State object
  USE State_Diag_Mod,      ONLY : DgnState   ! Derived type for Diagnostics State object
  USE State_Grid_Mod,      ONLY : GrdState   ! Derived type for Grid State object
  USE State_Met_Mod,       ONLY : MetState   ! Derived type for Meteorology State object
  USE ErrCode_Mod                            ! Error codes for success or failure
  USE Error_Mod                              ! For error checking

  !-----------------------------------------------------------------
  ! Parameters to define floating-point variables
  !-----------------------------------------------------------------
  USE PRECISION_MOD,       ONLY : fp, f4     ! Flexible precision

  use chem_mods,           only : nSlvd, slvd_Lst, slvd_ref_MMR

  ! Exit routine in CAM
  use cam_abortutils,      only : endrun

  use chem_mods,           only : nTracersMax
  use chem_mods,           only : nTracers
  use chem_mods,           only : tracerNames
  use chem_mods,           only : tracerLongNames
  use chem_mods,           only : adv_Mass
  use chem_mods,           only : mwRatio
  use chem_mods,           only : ref_MMR
  use chem_mods,           only : nSlsMax
  use chem_mods,           only : nSls
  use chem_mods,           only : slsNames
  use chem_mods,           only : slsLongNames
  use chem_mods,           only : sls_ref_MMR
  use chem_mods,           only : slsmwRatio
  use chem_mods,           only : nAerMax
  use chem_mods,           only : nAer
  use chem_mods,           only : aerNames
  use chem_mods,           only : aerAdvMass
  use chem_mods,           only : map2GC
  use chem_mods,           only : map2GC_Sls
  use chem_mods,           only : map2Idx
  use chem_mods,           only : map2MAM4

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
  TYPE(OptInput)             :: Input_Opt     ! Input Options object
  TYPE(ChmState),ALLOCATABLE :: State_Chm(:)  ! Chemistry State object
  TYPE(DgnState),ALLOCATABLE :: State_Diag(:) ! Diagnostics State object
  TYPE(GrdState),ALLOCATABLE :: State_Grid(:) ! Grid State object
  TYPE(MetState),ALLOCATABLE :: State_Met(:)  ! Meteorology State object
  TYPE(DgnList )             :: Diag_List     ! Diagnostics list object

  type(physics_buffer_desc), pointer :: hco_pbuf2d(:,:)  ! ptr to 2d pbuf

  ! Indices of critical species
  INTEGER                    :: iH2O, iO3, iCH4, iCO, iNO

  ! Indices in the physics buffer
  INTEGER                    :: NDX_PBLH      ! PBL height [m]
  INTEGER                    :: NDX_FSDS      ! Downward shortwave flux at surface [W/m2]
  INTEGER                    :: NDX_CLDTOP    ! Cloud top height [index]
  INTEGER                    :: NDX_CLDFRC    ! Cloud fraction [-]
  INTEGER                    :: NDX_PRAIN     ! Rain production rate [kg/kg/s]
  INTEGER                    :: NDX_NEVAPR    ! Total rate of precipitation evaporation  [kg/kg/s]
  INTEGER                    :: NDX_RPRDTOT   ! Convective total precip. production rate [kg/kg/s]
  INTEGER                    :: NDX_LSFLXPRC  ! Large-scale precip. at interface (liq + snw) [kg/m2/s]
  INTEGER                    :: NDX_LSFLXSNW  ! Large-scale precip. at interface (snow only) [kg/m2/s]

  ! Get constituent indices
  INTEGER :: ixCldLiq
  INTEGER :: ixCldIce

  ! Strings
  CHARACTER(LEN=255)              :: ThisLoc
  CHARACTER(LEN=255)              :: ErrMsg

#define ALLDDVEL_GEOSCHEM 1
#define OCNDDVEL_GEOSCHEM 0
#define OCNDDVEL_MOZART   0

! The following flags are only used if ALLDDVEL_GEOSCHEM is on
#define LANDTYPE_HEMCO    0
#define LANDTYPE_CLM      1

  ! Filenames to compute dry deposition velocities similarly to MOZART
  character(len=shr_kind_cl)  :: clim_soilw_file = 'clim_soilw_file'
  character(len=shr_kind_cl)  :: depvel_file     = ''
  character(len=shr_kind_cl)  :: depvel_lnd_file = 'depvel_lnd_file'
  character(len=shr_kind_cl)  :: season_wes_file = 'season_wes_file'

  character(len=shr_kind_cl) :: srf_emis_specifier(pcnst) = ''
  character(len=shr_kind_cl) :: ext_frc_specifier(pcnst) = ''

  character(len=24)  :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: srf_emis_cycle_yr  = 0
  integer            :: srf_emis_fixed_ymd = 0
  integer            :: srf_emis_fixed_tod = 0

  character(len=24)  :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' |  'INTERP_MISSING_MONTHS'
  integer            :: ext_frc_cycle_yr  = 0
  integer            :: ext_frc_fixed_ymd = 0
  integer            :: ext_frc_fixed_tod = 0

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
    use Species_Mod,         only : Species

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
    Type(ChmState)                 :: SC
    Type(GrdState)                 :: SG
    Type(OptInput)                 :: IO
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
    LOGICAL                        :: camout
    LOGICAL                        :: ic_from_cam2
    LOGICAL                        :: has_fixed_ubc
    LOGICAL                        :: has_fixed_ubflx

    INTEGER                        :: RC, IERR

    ! SDE 2018-05-02: This seems to get called before anything else
    ! that includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    IF ( MasterProc ) Write(iulog,'(a)') 'GCCALL CHEM_REGISTER'

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
    IO%SALA_rEdge_um(1)    = 0.01E+0_FP
    IO%SALA_rEdge_um(2)    = 0.50E+0_FP
    IO%SALC_rEdge_um(1)    = 0.50E+0_FP
    IO%SALC_rEdge_um(2)    = 8.00E+0_FP

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

    CALL Init_State_Chm( Input_Opt  = IO,  &
                         State_Chm  = SC,  &
                         State_Grid = SG,  &
                         RC         = RC  )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Init_State_Chm"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    map2GC = -1
    ref_MMR(:) = 0.0e+0_r8
    MWRatio(:) = 1.0e+0_r8
    tracerLongNames = ''

    DO I = 1, nTracersMax
       IF ( I .LE. nTracers ) THEN
           cnstName    = TRIM(tracerNames(I))
           N           = Ind_(cnstName)
           ThisSpc     => SC%SpcData(N)%Info
           lngName     = TRIM(ThisSpc%FullName)
           MWTmp       = REAL(ThisSpc%MW_g,r8)
           ref_VMR     = REAL(ThisSpc%BackgroundVV,r8)
           adv_Mass(I) = MWTmp
           ref_MMR(I)  = ref_VMR / (MWDry / MWTmp)
           IF ( ThisSpc%Is_Gas == .FALSE. ) THEN
               Write(cnstName, "(a,a)") 'GC_AER_', &
               to_upper(TRIM(tracerNames(I)))
               ! Aerosols that inherited from MAM do not need to be defined as a
               ! constituent
               ! For instance,
               ! - SOAGX is inherited from SOAG[0-4],
               !
               ! List of GEOS-Chem aerosols:
               ! DMS
               ! SO4
               ! SO4s
               ! MSA
               ! NH4
               ! NIT
               ! NITs
               ! BCPI
               ! OCPI
               ! BCPO
               ! OCPO
               ! DST[1-4]
               ! SALA
               ! SALC
               ! TSOA[0-3]
               ! ASOAN
               ! ASOA[1-3]
               ! SOAS
               ! SOAIE
               ! SOAME
               ! SOAGX
               ! SOAMG
               ! LVOCOA
               ! ISN1OA
               ! IONITA
               ! MONITA
               ! INDIOL
               ! BrSALA
               ! BrSALC
               ! ISALA
               ! ISALC
               ! AERI
               ! pFe
               ! TMMF, Update this
           ENDIF
       ELSEIF ( I .LE. (nTracers + nAer)) THEN
           ! Add MAM4 aerosols
           cnstName    = TRIM(aerNames(I - nTracers))
           lngName     = cnstName
           MWTmp       = aerAdvMass(I - nTracers)
           adv_Mass(I) = MWTmp
           ref_MMR(I)  = 1.0e-38_r8
       ELSE
           cnstName    = TRIM(tracerNames(I))
           lngName     = cnstName
           MWTmp       = 1000.0e+0_r8 * (0.001e+0_r8)
           adv_Mass(I) = MWTmp
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
       ! NOTE: In MOZART, this only gets called for tracers
       ! This is the call to add a "constituent"
       ! Special handlings
       IF ( cnstName == 'ACET' ) THEN
           cnstName = 'CH3COCH3'
       ELSEIF ( cnstName == 'ALD2' ) THEN
           cnstName = 'CH3CHO'
       ELSEIF ( cnstName == 'PRPE' ) THEN
           cnstName = 'C3H6'
       ELSEIF ( cnstName == 'HNO4' ) THEN
           cnstName = 'HO2NO2'
       ELSEIF ( cnstName == 'HNO2' ) THEN
           cnstName = 'HONO'
       ELSEIF ( cnstName == 'MP'   ) THEN
           cnstName = 'CH3OOH'
       ELSEIF ( cnstName == 'HAC'  ) THEN
           cnstName = 'HYAC'
       ELSEIF ( cnstName == 'GLYC' ) THEN
           cnstName = 'GLYALD'
       ELSEIF ( cnstName == 'MAP' ) THEN
           cnstName = 'CH3COOOH'
       ELSEIF ( cnstName == 'EOH' ) THEN
           cnstName = 'C2H5OH'
       ELSEIF ( cnstName == 'MGLY' ) THEN
           cnstName = 'CH3COCHO'
       ELSEIF ( cnstName == 'GLYX' ) THEN
           cnstName = 'GLYOXAL'
       ELSEIF ( cnstName == 'ACTA' ) THEN
           cnstName = 'CH3COOH'
       ELSEIF ( cnstName == 'TOLU' ) THEN
           cnstName = 'TOLUENE'
       ELSEIF ( cnstName == 'HCHO' ) THEN
           cnstName = 'CH2O'
       ENDIF
       If ( MasterProc ) Write(iulog,*) " Species = ", TRIM(cnstName)
       ! GEOS-Chem lumped species are not on restart file.
       ! Bromine, chlorine, iodine and halons species are missing
       ! from CESM restart file.
       ! These species will just be uniformily set to some low
       ! concentration.
       ! TMMF - 05/19/2020
       CALL cnst_add( cnstName, adv_Mass(I), cptmp, qmin, N,  &
                      readiv=ic_from_cam2, mixtype=mixtype,   &
                      cam_outfld=camout, molectype=molectype, &
                      fixed_ubc=has_fixed_ubc,                &
                      fixed_ubflx=has_fixed_ubflx,            &
                      longname=TRIM(lngName)                 )

       ! Add to GC mapping. When starting a timestep, we will want to update the
       ! concentration of State_Chm(x)%Species(1,iCol,iLev,m) with data from
       ! constituent n
       M = Ind_(TRIM(tracerNames(I)))
       IF ( M > 0 ) THEN
          map2GC(N)  = M
          map2Idx(N) = I
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
                map2MAM4(L,M) = Ind_('BCPI')
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
                SELECT CASE ( to_upper(xname_massptr(L,M)(4:4)) )
                   CASE ( '1' )
                      SELECT CASE ( to_upper(xname_massptr(L,M)(6:7)) )
                         CASE ( 'A1' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                         CASE ( 'A2' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                      END SELECT
                   CASE ( '2' )
                      SELECT CASE ( to_upper(xname_massptr(L,M)(6:7)) )
                         CASE ( 'A1' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                         CASE ( 'A2' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                      END SELECT
                   CASE ( '3' )
                      SELECT CASE ( to_upper(xname_massptr(L,M)(6:7)) )
                         CASE ( 'A1' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                         CASE ( 'A2' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                      END SELECT
                   CASE ( '4' )
                      SELECT CASE ( to_upper(xname_massptr(L,M)(6:7)) )
                         CASE ( 'A1' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                         CASE ( 'A2' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                      END SELECT
                   CASE ( '5' )
                      SELECT CASE ( to_upper(xname_massptr(L,M)(6:7)) )
                         CASE ( 'A1' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                         CASE ( 'A2' )
                            map2MAM4(L,M) = -1 !TMMF, Fill in
                      END SELECT
                END SELECT
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
                      map2MAM4(L,M) = -1 !TMMF, Fill in
                   CASE ( 'A4' )
                      map2MAM4(L,M) = -1 !TMMF, Fill in
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

    use cam_abortutils, only : endrun
    use units,          only : getunit, freeunit
    use namelist_utils, only : find_group_name
#if defined( MODAL_AERO_4MODE )
    use aero_model,     only : aero_model_readnl
    use dust_model,     only : dust_readnl
#endif
    use mpishorthand
    use gckpp_Model,    only : nSpec, Spc_Names
    use chem_mods,      only : drySpc_ndx

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
                           depvel_lnd_file,    &
                           ext_frc_cycle_yr,   &
                           ext_frc_specifier,  &
                           ext_frc_type,       &
                           season_wes_file,    &
                           srf_emis_cycle_yr,  &
                           srf_emis_specifier, &
                           srf_emis_type

    inputGeosPath='/glade/u/home/fritzt/input.geos.template'
    speciesDBPath='/glade/u/home/fritzt/species_database.yml'
    chemInputsDir='/glade/p/univ/umit0034/ExtData/CHEM_INPUTS/'


#if ( ALLDDVEL_GEOSCHEM + OCNDDVEL_GEOSCHEM + OCNDDVEL_MOZART != 1 )
    IF ( MasterProc ) THEN
        Write(iulog,'(/,a)') REPEAT( "=", 79 )
        Write(iulog,'(a)') " Preprocessor flags are not set correctly in chemistry.F90"
        Write(iulog,'(a)') " The user needs to decide how to compute dry deposition velocities"
        Write(iulog,'(a)') " Three options appear: "
        Write(iulog,'(a)') " + Let GEOS-Chem calculate all dry deposition velocities."
        Write(iulog,'(a)') "   Required setup:"
        Write(iulog,'(a)') "   ALLDDVEL_GEOSCHEM == 1"
        Write(iulog,'(a)') "   OCNDDVEL_GEOSCHEM == 0"
        Write(iulog,'(a)') "   OCNDDVEL_MOZART   == 0"
        Write(iulog,'(a)') " + Let CLM compute dry deposition velocities over land and let"
        Write(iulog,'(a)') "   GEOS-Chem compute velocities over ocean and ice"
        Write(iulog,'(a)') "   Required setup:"
        Write(iulog,'(a)') "   ALLDDVEL_GEOSCHEM == 0"
        Write(iulog,'(a)') "   OCNDDVEL_GEOSCHEM == 1"
        Write(iulog,'(a)') "   OCNDDVEL_MOZART   == 0"
        Write(iulog,'(a)') " + Let CLM compute dry deposition velocities over land and"
        Write(iulog,'(a)') "   compute velocities over ocean and ice in a similar way as"
        Write(iulog,'(a)') "   MOZART"
        Write(iulog,'(a)') "   Required setup:"
        Write(iulog,'(a)') "   ALLDDVEL_GEOSCHEM == 0"
        Write(iulog,'(a)') "   OCNDDVEL_GEOSCHEM == 0"
        Write(iulog,'(a)') "   OCNDDVEL_MOZART   == 1"
        Write(iulog,'(a)') REPEAT( "=", 79 )
        CALL ENDRUN('Incorrect definitions for dry deposition velocities')
    ENDIF
#endif
#if ( ALLDDVEL_GEOSCHEM && ( LANDTYPE_HEMCO + LANDTYPE_CLM != 1 ) )
    IF ( MasterProc ) THEN
        Write(iulog,'(/,a)') REPEAT( "=", 79 )
        Write(iulog,'(a)') REPEAT( "=", 79 )
        Write(iulog,'(a)') " Preprocessor flags are not set correctly in chemistry.F90"
        Write(iulog,'(a)') " Dry-deposition velocities are computed by GEOS-Chem"
        Write(iulog,'(a)') " The user needs to decide if land types should be from CLM or from HEMCO"
        CALL ENDRUN('Incorrect definitions for source of land type data')
    ENDIF
#endif

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

       Write(iulog,'(/,/, a)') 'Now defining GEOS-Chem tracers and dry deposition mapping...'

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
    CALL MPIBCAST( nTracers,    1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST( tracerNames, LEN(tracerNames(1))*nTracersMax, MPICHAR, 0, MPICOM )
    CALL MPIBCAST( nSls,        1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST( slsNames,    LEN(slsNames(1))*nSlsMax,        MPICHAR, 0, MPICOM )

    ! The following files are required to compute land maps, required to perform
    ! aerosol dry deposition
    CALL MPIBCAST(depvel_lnd_file, LEN(depvel_lnd_file), MPICHAR, 0, MPICOM)
    CALL MPIBCAST(clim_soilw_file, LEN(clim_soilw_file), MPICHAR, 0, MPICOM)
    CALL MPIBCAST(season_wes_file, LEN(season_wes_file), MPICHAR, 0, MPICOM)

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

    INTEGER :: I

    chem_implements_cnst = .false.

    DO I = 1, nTracers
       IF (TRIM(tracerNames(I)) .eq. TRIM(name)) THEN
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
    use physics_buffer, only : physics_buffer_desc, pbuf_get_index
    use cam_history,    only : addfld, add_default, horiz_only
    use chem_mods,      only : map2GC_dryDep, drySpc_ndx

    use mpishorthand
    use cam_abortutils, only : endrun

    use Phys_Grid,      only : get_Area_All_p
    use hycoef,         only : ps0, hyai, hybi

    use seq_drydep_mod, only : drydep_method, DD_XLND
#if ( OCNDDVEL_MOZART )
    use mo_drydep,      only : drydep_inti
#endif

#if defined( MODAL_AERO_4MODE )
    use aero_model,     only : aero_model_init
    use mo_drydep,      only : drydep_inti_landuse
    use modal_aero_data,only : ntot_amode, nspec_amode
    use modal_aero_data,only : xname_massptr
#endif

    use Input_Opt_Mod
    use State_Chm_Mod
    use State_Grid_Mod
    use State_Met_Mod
    use DiagList_Mod,   only : Init_DiagList, Print_DiagList
    use GC_Environment_Mod
    use GC_Grid_Mod,    only : SetGridFromCtrEdges

    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,  only : PI, PI_180, Re

    use Time_Mod,      only : Accept_External_Date_Time
    !use Time_Mod,      only : Set_Begin_Time,   Set_End_Time
    !use Time_Mod,      only : Set_Current_Time, Set_DiagB
    !use Transfer_Mod,  only : Init_Transfer
    use Linoz_Mod,     only : Linoz_Read

    use CMN_Size_Mod

    use Drydep_Mod,    only : Depname, Ndvzind
    use Pressure_Mod,  only : Accept_External_ApBp
    use Chemistry_Mod, only : Init_Chemistry
    use Ucx_Mod,       only : Init_Ucx
    use Input_mod,     only : Validate_Directories
#if   ( ALLDDVEL_GEOSCHEM && LANDTYPE_HEMCO )
    use Olson_Landmap_Mod
#endif
    use Mixing_Mod

    use GC_Emissions_Mod, only : GC_Emissions_Init

    TYPE(physics_state), INTENT(IN):: phys_state(BEGCHUNK:ENDCHUNK)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    ! Local variables

    !----------------------------
    ! Scalars
    !----------------------------

    ! Integers
    INTEGER               :: LCHNK(BEGCHUNK:ENDCHUNK), NCOL(BEGCHUNK:ENDCHUNK)
    INTEGER               :: IWAIT, IERR
    INTEGER               :: nX, nY, nZ
    INTEGER               :: iX, jY
    INTEGER               :: I, J, L, N
    INTEGER               :: RC
    INTEGER               :: nLinoz

    ! Logicals
    LOGICAL               :: prtDebug

    ! Strings
    CHARACTER(LEN=255)    :: historyConfigFile
    CHARACTER(LEN=255)    :: SpcName

    ! Grid setup
    REAL(fp)              :: lonVal,  latVal
    REAL(fp)              :: dLonFix, dLatFix
    REAL(f4), ALLOCATABLE :: lonMidArr(:,:),  latMidArr(:,:)
    REAL(f4), ALLOCATABLE :: lonEdgeArr(:,:), latEdgeArr(:,:)
    REAL(r8), ALLOCATABLE :: linozData(:,:,:,:)

    REAL(r8), ALLOCATABLE :: Col_Area(:)
    REAL(fp), ALLOCATABLE :: Ap_CAM_Flip(:), Bp_CAM_Flip(:)

    REAL(r8), POINTER     :: SlsPtr(:,:,:)


    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at GEOS-Chem (in chemistry/pp_geoschem/chemistry.F90)'

    ! LCHNK: which chunks we have on this process
    LCHNK = PHYS_STATE%LCHNK
    ! NCOL: number of atmospheric columns for each chunk
    NCOL  = PHYS_STATE%NCOL

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
       State_Grid(I)%NY = nY
       State_Grid(I)%NZ = nZ

       ! Initialize GEOS-Chem horizontal grid structure
       CALL GC_Init_Grid( Input_Opt  = Input_Opt,      &
                          State_Grid = State_Grid(I),  &
                          RC         = RC          )

       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
          CALL Error_Stop( ErrMsg, ThisLoc )
       ENDIF

       ! Define more variables for State_Grid
       ! TMMF, might need tweaking
       State_Grid(I)%MaxTropLev  = MIN(40, nZ)
       State_Grid(I)%MaxStratLev = MIN(59, nZ)

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
    CALL GC_Allocate_All ( Input_Opt      = Input_Opt,            &
                           State_Grid     = State_Grid(BEGCHUNK), &
                           value_I_Lo     = 1,                    &
                           value_J_Lo     = 1,                    &
                           value_I_Hi     = nX,                   &
                           value_J_Hi     = nY,                   &
                           value_IM       = nX,                   &
                           value_JM       = nY,                   &
                           value_LM       = nZ,                   &
                           value_IM_WORLD = nX,                   &
                           value_JM_WORLD = nY,                   &
                           value_LM_WORLD = nZ,                   &
                           value_LLSTRAT  = 59,                   & !TMMF
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
    ! For now, TMMF
    Input_Opt%LConv                  = .False.
    Input_Opt%LTurb                  = .True.
    Input_Opt%LNLPBL                 = .True.

    ! Now READ_EMISSIONS_MENU
    Input_Opt%LEmis                  = .False.
    Input_Opt%HCOConfigFile          = 'HEMCO_Config.rc'

    ! Set surface VMRs - turn this off so that CAM can handle it
    Input_Opt%LCH4Emis               = .False.
    Input_Opt%LCH4SBC                = .False.

    ! Set initial conditions
    Input_Opt%LSetH2O                = .True.

    ! Now READ_AEROSOL_MENU
    Input_Opt%LSulf               = .True.
    Input_Opt%LMetalcatSO2        = .True.
    Input_Opt%LCarb               = .True.
    Input_Opt%LBrC                = .False.
    Input_Opt%LSOA                = .False.
    Input_Opt%LSVPOA              = .False.
    Input_Opt%LDust               = .True.
    Input_Opt%LDstUp              = .False.
    Input_Opt%LSSalt              = .True.
    Input_Opt%SalA_rEdge_um(1)    = 0.01e+0_fp
    Input_Opt%SalA_rEdge_um(2)    = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(1)    = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(2)    = 8.00e+0_fp
    Input_Opt%LMPOA               = .False.
    ! For now, disable solid PSCs and strat aerosol settling
    ! Our treatment of the stratosphere isn't really sophisticated
    ! enough to warrant it yet
    Input_Opt%LGravStrat          = .False.
    Input_Opt%LSolidPSC           = .False.
    Input_Opt%LHomNucNAT          = .False.
    Input_Opt%T_NAT_Supercool     = 3.0e+0_fp
    Input_Opt%P_Ice_Supersat      = 1.2e+0_fp
    Input_Opt%LPSCChem            = .True.
    Input_Opt%LStratOD            = .True.
    Input_Opt%hvAerNIT            = .False.
    Input_Opt%hvAerNIT_JNIT       = .False.
    Input_Opt%hvAerNIT_JNITs      = .False.
    Input_Opt%JNITChanA           = 0e+0_fp
    Input_Opt%JNITChanB           = 0e+0_fp

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
    Input_Opt%CO2_Level              = 390.0_fp
    Input_Opt%CO2_Ref                = 390.0_fp

    ! Now READ_CHEMISTRY_MENU
    Input_Opt%LChem                  = .True.
    Input_Opt%LSChem                 = .False. ! .True. !TMMF
    Input_Opt%LLinoz                 = .True.
    Input_Opt%LSynoz                 = .True.
    Input_Opt%LUCX                   = .True.
    Input_Opt%LActiveH2O             = .True.
    Input_Opt%Use_Online_O3          = .True.
    ! Expect to get total overhead ozone, although it shouldn not
    ! make too much of a difference since we want to use "full-UCX"
    Input_Opt%Use_O3_from_Met        = .True.
    Input_Opt%Use_TOMS_O3            = .False.
    Input_Opt%Gamma_HO2              = 0.2e+0_fp

    Input_Opt%LPRT                   = .False.

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
       DEALLOCATE(linozData)
    ENDIF

    ! Note: The following calculations do not setup the gridcell areas.
    !       In any case, we will need to be constantly updating this grid
    !       to compensate for the "multiple chunks per processor" element
    ALLOCATE(lonMidArr(nX,nY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonMidArr')
    ALLOCATE(lonEdgeArr(nX+1,nY+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonEdgeArr')
    ALLOCATE(latMidArr(nX,nY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating latMidArr')
    ALLOCATE(latEdgeArr(nX+1,nY+1), STAT=IERR)
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
       dLatFix   = 180.0e+0_fp / REAL(nY,fp)
       DO I = 1, nX
          ! Center of box, assuming dateline edge
          lonVal = -180.0e+0_fp + (REAL(I-1,fp)*dLonFix)
          DO J = 1, nY
             ! Center of box, assuming regular cells
             latVal = -90.0e+0_fp + (REAL(J-1,fp)*dLatFix)
             lonMidArr(I,J)  = REAL((lonVal + (0.5e+0_fp * dLonFix)) * PI_180, f4)
             latMidArr(I,J)  = REAL((latVal + (0.5e+0_fp * dLatFix)) * PI_180, f4)

             ! Edges of box, assuming regular cells
             lonEdgeArr(I,J) = REAL(lonVal * PI_180, f4)
             latEdgeArr(I,J) = REAL(latVal * PI_180, f4)
          ENDDO
          ! Edges of box, assuming regular cells
          lonEdgeArr(I,nY+1)  = REAL((lonVal + dLonFix) * PI_180, f4)
          latEdgeArr(I,nY+1)  = REAL((latVal + dLatFix) * PI_180, f4)
       ENDDO
       DO J = 1, nY+1
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
    DEALLOCATE(lonMidArr)
    DEALLOCATE(latMidArr)
    DEALLOCATE(lonEdgeArr)
    DEALLOCATE(latEdgeArr)

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

    ! Debug output
    IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a READ_INPUT_FILE' )

    historyConfigFile = 'HISTORY.rc' ! InputOpt not yet initialized
    !TMMF need to pass input.geos path
    !CALL Init_DiagList( MasterProc, historyConfigFile, Diag_List, RC )
    !IF ( RC /= GC_SUCCESS ) THEN
    !   ErrMsg = 'Error encountered in "Init_DiagList"!'
    !   CALL Error_Stop( ErrMsg, ThisLoc )
    !ENDIF

    !!### Print diagnostic list if needed for debugging
    !IF ( prtDebug ) CALL Print_DiagList( Diag_List, RC )

    DO I = BEGCHUNK, ENDCHUNK
       Input_Opt%amIRoot = (MasterProc .AND. (I == BEGCHUNK))

       CALL GC_Init_StateObj( Diag_List  = Diag_List,     &  ! Diagnostic list obj
                              Input_Opt  = Input_Opt,     &  ! Input Options
                              State_Chm  = State_Chm(I),  &  ! Chemistry State
                              State_Diag = State_Diag(I), &  ! Diagnostics State
                              State_Grid = State_Grid(I), &  ! Grid State
                              State_Met  = State_Met(I),  &  ! Meteorology State
                              RC         = RC            )   ! Success or failure

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
    &                   State_Grid = State_Grid(BEGCHUNK), &  ! Grid State
    &                   RC         = RC                   )   ! Success or failure

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "GC_Init_Extra"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

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

#if ( OCNDDVEL_MOZART )
       !==============================================================
       ! The following line should only be called if we compute
       ! velocities over the ocean and ice in a MOZART-like way.
       ! Thibaud M. Fritz - 26 Feb 2020
       !==============================================================

       IF ( drydep_method == DD_XLND ) THEN
          CALL drydep_inti( depvel_lnd_file, &
                            clim_soilw_file, &
                            season_wes_file )
       ELSE
          Write(iulog,'(a,a)') ' drydep_method is set to: ', TRIM(drydep_method)
          CALL ENDRUN('drydep_method must be DD_XLND to compute dry deposition' // &
              ' velocities similarly to MOZART over ocean and ice!')
       ENDIF
#endif

    ENDIF

#if defined( MODAL_AERO_4MODE )
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

    ! Set grid-cell area
    DO I = BEGCHUNK, ENDCHUNK
       ALLOCATE(Col_Area(NCOL(I)), STAT=IERR)
       IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Col_Area')

       CALL Get_Area_All_p(I, NCOL(I), Col_Area)

       ! Set default value (in case of chunks with fewer columns)
       State_Grid(I)%Area_M2 = 1.0e+10_fp
       DO iX = 1, nX
       DO jY = 1, NCOL(I)
          State_Grid(I)%Area_M2(iX,jY) = REAL(Col_Area(jY) * Re**2,fp)
       ENDDO
       ENDDO

       DEALLOCATE(Col_Area)

       ! Copy to State_Met(I)%Area_M2
       State_Met(I)%Area_M2 = State_Grid(I)%Area_M2
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
    CALL Accept_External_ApBp( State_Grid = State_Grid(BEGCHUNK), &  ! Grid State
                               ApIn       = Ap_CAM_Flip,          &  ! "A" term for hybrid grid
                               BpIn       = Bp_CAM_Flip,          &  ! "B" term for hybrid grid
                               RC         = RC                   )   ! Success or failure

    ! Print vertical coordinates
    IF ( MasterProc ) THEN
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
       WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) Ap_CAM_Flip(1:State_Grid(BEGCHUNK)%NZ+1)
       WRITE( 6, '(a)'   )
       WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) Bp_CAM_Flip(1:State_Grid(BEGCHUNK)%NZ+1)
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
    ENDIF

    ! Trapping errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Accept_External_ApBp"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DEALLOCATE(Ap_CAM_Flip,Bp_CAM_Flip)

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
                      State_Grid = State_Grid(BEGCHUNK) )
    ENDIF

    ! Get the index of H2O
    iH2O = Ind_('H2O')
    iO3  = Ind_('O3')
    iCH4 = Ind_('CH4')
    iCO  = Ind_('CO')
    iNO  = Ind_('NO')

    ! Get indices for physical fields in physics buffer
    NDX_PBLH     = Pbuf_Get_Index('pblh'     )
    NDX_FSDS     = Pbuf_Get_Index('FSDS'     )
    NDX_CLDTOP   = Pbuf_Get_Index('CLDTOP'   )
    NDX_CLDFRC   = Pbuf_Get_Index('CLD'      )
    NDX_PRAIN    = Pbuf_Get_Index('PRAIN'    )
    NDX_NEVAPR   = Pbuf_Get_Index('NEVAPR'   )
    NDX_RPRDTOT  = Pbuf_Get_Index('RPRDTOT'  )
    NDX_LSFLXPRC = Pbuf_Get_Index('LS_FLXPRC')
    NDX_LSFLXSNW = Pbuf_Get_Index('LS_FLXSNW')

    ! Get cloud water indices
    CALL cnst_get_ind('CLDLIQ', ixCldLiq)
    CALL cnst_get_ind('CLDICE', ixCldIce)

    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default
    ! Add all species as output fields if desired
    DO I = 1, nTracers
       SpcName = TRIM(tracerNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', TRIM(tracerLongNames(I))//' concentration')
       IF (TRIM(SpcName) == 'O3') THEN
          CALL Add_Default ( TRIM(SpcName), 1, ' ')
       ENDIF
    ENDDO

    DO I =1, nSls
       SpcName = TRIM(slsNames(I))
       CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', TRIM(slsLongNames(I))//' concentration')
       !CALL Add_Default(TRIM(SpcName), 1, '')
    ENDDO

    ! Initialize emissions interface (this will eventually handle HEMCO)
    CALL GC_Emissions_Init

    !CALL AddFld ( 'BCPI', (/'lev'/), 'A', 'mole/mole', trim('BCPI')//' mixing ratio' )
    !CALL Add_Default ( 'BCPI',   1, ' ')

    hco_pbuf2d => pbuf2d

    If ( MasterProc ) Write(iulog,*) "hco_pbuf2d now points to pbuf2d"

    hco_pbuf2d => pbuf2d

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT'

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
                            RADIATION  =  DT_MIN      )
        DT_MIN_LAST = DT_MIN
     ENDIF

  end subroutine

!===============================================================================

  subroutine chem_timestep_tend( State, ptend, cam_in, cam_out, dT, pbuf,  fh2o )

    use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use physics_buffer,      only : pbuf_get_chunk, pbuf_get_index
    use cam_history,         only: outfld
    use camsrfexch,          only: cam_in_t, cam_out_t

    use phys_grid,           only: get_ncols_p, get_rlat_all_p, get_rlon_all_p

    use chem_mods,           only: drySpc_ndx, map2GC_dryDep
#if defined( MODAL_AERO_4MODE )
    use modal_aero_data,     only : ntot_amode, nspec_amode
    use modal_aero_data,     only : lmassptr_amode
#endif

    use Olson_Landmap_Mod,   only: Compute_Olson_Landmap
    use Modis_LAI_Mod,       only: Compute_XLAI
    use CMN_Size_Mod,        only: NSURFTYPE
#if   ( ALLDDVEL_GEOSCHEM || OCNDDVEL_GEOSCHEM )
    use Drydep_Mod,          only: Do_Drydep
#elif ( OCNDDVEL_MOZART )
    use mo_drydep,           only: drydep_update, drydep_fromlnd
#endif
    use Drydep_Mod,          only: DEPNAME
    use Drydep_Mod,          only: Update_DryDepSav
    use Mixing_Mod

    use Calc_Met_Mod,        only: Set_Dry_Surface_Pressure
    use Calc_Met_Mod,        only: AirQnt
    use GC_Grid_Mod,         only: SetGridFromCtr
    use Pressure_Mod,        only: Set_Floating_Pressures
    use Pressure_Mod,        only: Accept_External_Pedge
    use Time_Mod,            only: Accept_External_Date_Time
    use Strat_chem_Mod,      only: Init_Strat_Chem
    use Toms_Mod,            only: Compute_Overhead_O3
    use Chemistry_Mod,       only: Do_Chemistry
    use Wetscav_Mod,         only: Setup_Wetscav, Do_WetDep
    use CMN_Size_Mod,        only: PTop
    use PBL_Mix_Mod,         only: Compute_PBL_Height

    use Tropopause,          only: Tropopause_findChemTrop, Tropopause_Find

    ! For calculating SZA
    use Orbit,               only: zenith
    use Time_Manager,        only: Get_Curr_Calday, Get_Curr_Date

    ! Calculating relative humidity
    use WV_Saturation,       only: QSat
    use PhysConst,           only: MWDry

    ! Grid area
    use PhysConst,           only: Gravit
    use PhysConstants,       only: Re
    use Phys_Grid,           only: get_area_all_p, get_lat_all_p, get_lon_all_p

    use short_lived_species, only : get_short_lived_species
    use short_lived_species, only : set_short_lived_species

    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,       only: PI, PI_180, g0

    use rad_constituents,    only: rad_cnst_get_info

    REAL(r8),            INTENT(IN)    :: dT          ! Time step
    TYPE(physics_state), INTENT(IN)    :: State       ! Physics State variables
    TYPE(physics_ptend), INTENT(OUT)   :: ptend       ! indivdual parameterization tendencies
    TYPE(cam_in_t),      INTENT(INOUT) :: cam_in
    TYPE(cam_out_t),     INTENT(IN)    :: cam_out
    TYPE(physics_buffer_desc), POINTER :: pbuf(:)
    REAL(r8), OPTIONAL,  INTENT(OUT)   :: fh2o(PCOLS) ! h2o flux to balance source from chemistry

    ! Initial MMR for all species
    REAL(r8) :: MMR_Beg(PCOLS,PVER,nSls+nTracers)
    REAL(r8) :: MMR_End(PCOLS,PVER,nSls+nTracers)
    REAL(r8) :: MMR_TEnd(PCOLS,PVER,nSls+nTracers)

    ! Logical to apply tendencies to mixing ratios
    LOGICAL :: lq(pcnst)

    ! Indexing
    INTEGER :: I, J, K, L, N, M, P
    INTEGER :: nX, nY, nZ

    INTEGER :: LCHNK, NCOL

    REAL(r8), DIMENSION(State%NCOL) :: &
        CSZA,                          &              ! Cosine of solar zenith angle
        Zsurf,                         &              ! Surface height
        Rlats, Rlons                                  ! Chunk latitudes and longitudes (radians)

    REAL(r8), POINTER :: PblH(:)                      ! PBL height on each chunk [m]
    REAL(r8), POINTER :: cldTop(:)                    ! Cloud top height [?]
    REAL(r8), POINTER :: cldFrc(:,:)                  ! Cloud fraction [-]
    REAL(r8), POINTER :: Fsds(:)                      ! Downward shortwave flux at surface [W/m2]
    REAL(r8), POINTER :: PRain(:,:)                   ! Total stratiform precip. prod. (rain + snow) [kg/kg/s]
    REAL(r8), POINTER :: RprdTot(:,:)                 ! Total convective precip. prod. (rain + snow) [kg/kg/s]
    REAL(r8), POINTER :: NEvapr(:,:)                  ! Evaporation of total precipitation (rain + snow) [kg/kg/s]
    REAL(r8), POINTER :: LsFlxPrc(:,:)                ! Large-scale downward precip. flux at interface (rain + snow) [kg/m2/s]
    REAL(r8), POINTER :: LsFlxSnw(:,:)                ! Large-scale downward precip. flux at interface (snow only) [kg/m2/s]

    REAL(r8)          :: RelHum(State%NCOL, PVER)     ! Relative humidity [0-1]
    REAL(r8)          :: SatV  (State%NCOL, PVER)     ! Work arrays
    REAL(r8)          :: SatQ  (State%NCOL, PVER)     ! Work arrays
    REAL(r8)          :: qH2O  (State%NCOL, PVER)     ! Specific humidity [kg/kg]
    REAL(r8)          :: H2OVMR(State%NCOL, PVER)     ! H2O volume mixing ratio
#if ( OCNDDVEL_MOZART )
    REAL(r8)          :: windSpeed(State%NCOL)        ! Wind speed at ground level [m/s]
    REAL(r8)          :: potT(State%NCOL)             ! Potential temperature [K]

    INTEGER           :: latndx(PCOLS)
    INTEGER           :: lonndx(PCOLS)

    ! For MOZART's dry deposition over ocean and ice
    ! Deposition velocity (cm/s)
    REAL(r8)          :: MOZART_depVel(State%NCOL, nTracersMax)
    ! Deposition flux (/cm^2/s)
    REAL(r8)          :: MOZART_depFlx(State%NCOL, nTracersMax)
#endif
    REAL(r8), PARAMETER :: zlnd  = 0.01_r8   ! Roughness length for soil [m]
    REAL(r8), PARAMETER :: zslnd = 0.0024_r8 ! Roughness length for snow [m]
    REAL(r8), PARAMETER :: zsice = 0.0400_r8 ! Roughness length for sea ice [m]
    REAL(r8), PARAMETER :: zocn  = 0.0001_r8 ! Roughness length for oean [m]

    ! Because of strat chem
    LOGICAL, SAVE :: SCHEM_READY = .FALSE.

    REAL(f4)      :: lonMidArr(1,PCOLS), latMidArr(1,PCOLS)
    INTEGER       :: iMaxLoc(1)

    REAL(r8)      :: Col_Area(State%NCOL)

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

    ! For archiving
    CHARACTER(LEN=255) :: SpcName
    REAL(r8)           :: VMR(State%NCOL,PVER)

    REAL(r8)           :: SlsData(State%NCOL, PVER, nSls)

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
    INTEGER            :: RC

#if ( LANDTYPE_HEMCO )
    INTEGER            :: tmpIdx
    character(len=255) :: name
    real(r8), pointer  :: pbuf_ik(:,:)     ! Pointer to pbuf data  (/pcols,pver/)
#endif

    ! LCHNK: which chunk we have on this process
    LCHNK = State%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    NCOL  = State%NCOL

    ! Am I the first chunk on the first CPU?
    rootChunk = ( MasterProc.and.(LCHNK==BEGCHUNK) )

    ! Count the number of steps which have passed
    IF (LCHNK.EQ.BEGCHUNK) iStep = iStep + 1

    ! Need to update the timesteps throughout the code
    CALL GC_Update_Timesteps(dT)


    ! For safety's sake
    PTop = State%Pint(1,1)*0.01e+0_fp

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
    Rlats(1:nY) = State%Lat(1:nY)
    Rlons(1:nY) = State%Lon(1:nY)

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
    State_Grid(LCHNK)%Area_M2 = 1.0e+10_fp
    DO J = 1, nY
       State_Grid(LCHNK)%Area_M2(1,J) = REAL(Col_Area(J) * Re**2,fp)
    ENDDO
    State_Met(LCHNK)%Area_M2 = State_Grid(LCHNK)%Area_M2

    ! 2. Copy tracers into State_Chm
    ! Data was received in kg/kg dry
    State_Chm(LCHNK)%Spc_Units = 'kg/kg dry'
    ! Initialize ALL State_Chm species data to zero, not just tracers
    State_Chm(LCHNK)%Species = 0.0e+0_fp

    lq(:) = .FALSE.

    ! Map and flip gaseous species
    MMR_Beg = 0.0e+0_r8
    DO N = 1, pcnst
       M = map2GC(N)
       IF (M > 0) THEN
          DO J = 1, nY
          DO K = 1, nZ
             MMR_Beg(J,K,M) = State%q(J,nZ+1-K,N)
             State_Chm(LCHNK)%Species(1,J,K,M) = REAL(MMR_Beg(J,K,M),fp)
          ENDDO
          ENDDO
          lq(N) = .TRUE.
       ENDIF
    ENDDO

    ! Retrieve previous value of species data
    SlsData(:,:,:) = 0.0e+0_r8
    CALL get_short_lived_species( SlsData, LCHNK, nY, pbuf )

    ! Map and flip gaseous short-lived species
    DO N = 1, nSls
       M = map2GC_Sls(N)
       IF (M > 0) THEN
          DO J = 1, nY
          DO K = 1, nZ
             State_Chm(LCHNK)%Species(1,J,K,M) = REAL(SlsData(J,nZ+1-K,N),fp)
          ENDDO
          ENDDO
       ENDIF
    ENDDO

#if defined( MODAL_AERO_4MODE )
    ! TMMF - This needs more indices to MMR_Beg and MMR_End
    ! Map and flip aerosols
    DO M = 1, ntot_amode
       DO L = 1, nspec_amode(M)
          P = map2MAM4(L,M)
          N = lmassptr_amode(L,M)
          ! Multiple MAM4 bins are mapped to same GEOS-Chem species
          IF ( P > 0 ) THEN
             DO J = 1, nY
             DO K = 1, nZ
                !MMR_Beg(J,K,M) = State%q(J,nZ+1-K,N) 
                State_Chm(LCHNK)%Species(1,J,K,P) = State_Chm(LCHNK)%Species(1,J,K,P) &
                                                  + REAL(state%q(J,nZ+1-K,N),fp)
             ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
#endif

    ! Initialize tendency array
    CALL Physics_ptend_init(ptend, State%psetcols, 'chemistry', lq=lq)

    ! Calculate COS(SZA)
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
    CALL pbuf_get_field( pbuf, NDX_RPRDTOT,  RprdTot,  START=(/1,1/),         KOUNT=(/NCOL,PVER/))
    CALL pbuf_get_field( pbuf, NDX_LSFLXPRC, LsFlxPrc, START=(/1,1/),         KOUNT=(/NCOL,PVERP/))
    CALL pbuf_get_field( pbuf, NDX_LSFLXSNW, LsFlxSnw, START=(/1,1/),         KOUNT=(/NCOL,PVERP/))

    ! Get VMR and MMR of H2O
    H2OVMR = 0.0e0_fp
    qH2O   = 0.0e0_fp
    ! Note MWDRY = 28.966 g/mol
    DO J = 1, nY
    DO L = 1, nZ
       qH2O(J,L) = REAL(State_Chm(LCHNK)%Species(1,J,L,iH2O),r8)
       H2OVMR(J,L) = qH2O(J,L) * MWDry / 18.016e+0_fp
    ENDDO
    ENDDO

    ! Calculate RH (range 0-1, note still level 1 = TOA)
    relHum(:,:) = 0.0e+0_r8
    CALL QSat(State%T(:nY,:), State%Pmid(:nY,:), SatV, SatQ)
    DO J = 1, nY
    DO L = 1, nZ
       relHum(J,L) = 0.622e+0_r8 * H2OVMR(J,L) / SatQ(J,L)
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
       ! Convert water mixing ratio [kg/kg] to water content [g/m^3]
       IF ( ( State%Q(J,L,ixCldLiq) + State%Q(J,L,ixCldIce) ) * &
            State%PMid(J,L) / (State%T(J,L) * 287.0e+00_r8) * 1.0e+03_r8 <= 0.01_r8 .AND. &
            cldFrc(J,L) /= 0.0e+00_r8 ) THEN
          cld(J,L) = 0.0e+00_r8
       ELSE
          cld(J,L) = cldFrc(J,L)
       ENDIF
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
             ! Tau(K) = 3/2 * Q(K) * (Pint(K+1) - Pint(K)) / (re * rho_w * g )
             ! Tau(K) = Q(K) * (Pint(K+1) - Pint(K)) * Cnst
             !
             ! Unit check:                    |
             ! Q    : [kg H2O/kg air]         |
             ! Pint : [Pa]=[kg air/m/s^2]     |
             ! re   : [m]                     |   = 1.0e-5
             ! rho_w: [kg H2O/m^3]            |   = 1.0e+3
             ! g    : [m/s^2]                 |   = 9.81
             !
             TauClw(J,L) = State%Q(J,L,ixCldLiq)               &
                         * (State%Pint(J,L+1)-State%Pint(J,L)) &
                         * cnst
             TauClw(J,L) = MAX(TauClw(J,L), 0.0e+00_r8)
             TauCli(J,L) = State%Q(J,L,ixCldIce)               &
                         * (State%Pint(J,L+1)-State%Pint(J,L)) &
                         * cnst
             TauCli(J,L) = MAX(TauCli(J,L), 0.0e+00_r8)

          ENDDO
       ENDIF
    ENDDO

    ! Retrieve tropopause level
    Trop_Lev = 0.0e+0_r8
    CALL Tropopause_FindChemTrop(State, Trop_Lev)
    ! Back out the pressure
    Trop_P = 1000.0e+0_r8
    DO J = 1, nY
       Trop_P(J) = State%PMid(J,Trop_Lev(J)) * 0.01e+0_r8
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
    State_Met(LCHNK)%ALBD      (1,:) = cam_in%Asdir(:)

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
    State_Met(LCHNK)%EFLUX     (1,:) = cam_in%Lhf(:)
    State_Met(LCHNK)%HFLUX     (1,:) = cam_in%Shf(:)

    ! Field      : LandTypeFrac
    ! Description: Olson fraction per type
    ! Unit       : - (between 0 and 1)
    ! Dimensions : nX, nY, NSURFTYPE
    ! Note       : Index 1 is water
#if   ( LANDTYPE_CLM )
    ! Fill in water
    State_Met(LCHNK)%LandTypeFrac(1,:, 1) = cam_in%ocnFrac(:)     &
                                          + cam_in%iceFrac(:)
#if   ( ALLDDVEL_GEOSCHEM )
    CALL getLandTypes( cam_in,         &
                       nY,             &
                       State_Met(LCHNK) )
#endif
#elif ( LANDTYPE_HEMCO )
    DO N = 1, NSURFTYPE
       Write(name, '(a,i2.2)') 'HCO_LANDTYPE', N-1
       If ( MasterProc ) Write(iulog,*) " Getting ", TRIM(name)

       tmpIdx = pbuf_get_index(name, rc)
       IF ( tmpIdx < 0 ) THEN
          ! there is an error here and the field was not found
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(name)
       ELSE
          CALL pbuf_get_field(pbuf, tmpIdx, pbuf_ik)
          DO J = 1, nY
             State_Met(LCHNK)%LandTypeFrac(1,J,N) = pbuf_ik(J,nZ)
             ! 2-D data is stored in the 1st level of a
             ! 3-D array due to laziness
          ENDDO
       ENDIF

       Write(name, '(a,i2.2)') 'HCO_XLAI', N-1
       If ( MasterProc ) Write(iulog,*) " Getting ", TRIM(name)

       tmpIdx = pbuf_get_index(name, rc)
       IF ( tmpIdx < 0 ) THEN
          ! there is an error here and the field was not found
          IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(name)
       ELSE
          CALL pbuf_get_field(pbuf, tmpIdx, pbuf_ik)
          DO J = 1, nY
             State_Met(LCHNK)%XLAI_NATIVE(1,J,N) = pbuf_ik(J,nZ)
             ! 2-D data is stored in the 1st level of a
             ! 3-D array due to laziness
          ENDDO
       ENDIF
    ENDDO
#endif

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
    State_Met(LCHNK)%FRCLND    (1,:) = 1.e+0_fp - &
                    State_Met(LCHNK)%LandTypeFrac(1,:,1) ! Olson Land Fraction
    State_Met(LCHNK)%FRLAND    (1,:) = cam_in%landFrac(:)
    State_Met(LCHNK)%FROCEAN   (1,:) = cam_in%ocnFrac(:) + cam_in%iceFrac(:)
    State_Met(LCHNK)%FRSEAICE  (1,:) = cam_in%iceFrac(:)
#if   ( LANDTYPE_CLM )
    State_Met(LCHNK)%FRLAKE    (1,:) = cam_in%lwtgcell(:,3) + &
                                       cam_in%lwtgcell(:,4)
    State_Met(LCHNK)%FRLANDIC  (1,:) = cam_in%lwtgcell(:,2)
    State_Met(LCHNK)%FRSNO     (1,:) = 0.0e+0_fp
#else
    State_Met(LCHNK)%FRLAKE    (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%FRLANDIC  (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%FRSNO     (1,:) = 0.0e+0_fp
#endif

    ! Field      : GWETROOT, GWETTOP
    ! Description: Root and top soil moisture
    ! Unit       : -
    ! Dimensions : nX, nY
    State_Met(LCHNK)%GWETROOT  (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%GWETTOP   (1,:) = 0.0e+0_fp

    ! Field      : LAI
    ! Description: Leaf area index
    ! Unit       : m^2/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%LAI       (1,:) = 0.0e+0_fp

    ! Field      : PARDR, PARDF
    ! Description: Direct and diffuse photosynthetically active radiation
    ! Unit       : W/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PARDR     (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%PARDF     (1,:) = 0.0e+0_fp

    ! Field      : PBLH
    ! Description: PBL height
    ! Unit       : m
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PBLH      (1,:) = PblH(:nY)

    ! Field      : PHIS
    ! Description: Surface geopotential height
    ! Unit       : m
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PHIS      (1,:) = State%Phis(:)

    ! Field      : PRECANV, PRECCON, PRECLSC, PRECTOT
    ! Description: Anvil precipitation @ ground
    !              Convective precipitation @ ground
    !              Large-scale precipitation @ ground
    !              Total precipitation @ ground
    ! Unit       : kg/m^2/s
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PRECANV   (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%PRECCON   (1,:) = cam_out%Precc(:)
    State_Met(LCHNK)%PRECLSC   (1,:) = cam_out%Precl(:)
    State_Met(LCHNK)%PRECTOT   (1,:) = cam_out%Precc(:) + cam_out%Precl(:)

    ! Field      : TROPP
    ! Description: Tropopause pressure
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TROPP     (1,:) = Trop_P(:)

    ! Field      : PS1_WET, PS2_WET
    ! Description: Wet surface pressure at start and end of timestep
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%PS1_WET   (1,:) = State%ps(:)*0.01e+0_fp
    State_Met(LCHNK)%PS2_WET   (1,:) = State%ps(:)*0.01e+0_fp

    ! Field      : SLP
    ! Description: Sea level pressure
    ! Unit       : hPa
    ! Dimensions : nX, nY
    State_Met(LCHNK)%SLP       (1,:) = State%ps(:)*0.01e+0_fp

    ! Field      : TS, TSKIN
    ! Description: Surface temperature, surface skin temperature
    ! Unit       : K
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TS        (1,:) = cam_in%TS(:)
    State_Met(LCHNK)%TSKIN     (1,:) = cam_in%TS(:)

    ! Field      : SWGDN
    ! Description: Incident radiation @ ground
    ! Unit       : W/m^2
    ! Dimensions : nX, nY
    State_Met(LCHNK)%SWGDN     (1,:) = fsds(:)

    ! Field      : TO3
    ! Description: Total overhead ozone column
    ! Unit       : DU
    ! Dimensions : nX, nY
    State_Met(LCHNK)%TO3       (1,:) = 300.0e+0_fp ! TMMF

    ! Field      : SNODP, SNOMAS
    ! Description: Snow depth, snow mass
    ! Unit       : m, kg/m^2
    ! Dimensions : nX, nY
    ! Note       : Conversion from m to kg/m^2
    !              \rho_{ice} = 916.7 kg/m^3
    State_Met(LCHNK)%SNODP     (1,:) = snowDepth(:)
    State_Met(LCHNK)%SNOMAS    (1,:) = snowDepth(:) * 916.7e+0_r8

    ! Field      : SUNCOS, SUNCOSmid
    ! Description: COS(solar zenith angle) at current time and midpoint
    !              of chemistry timestep
    ! Unit       : -
    ! Dimensions : nX, nY
    ! Note       : Compute tendency in -/s (tmmf, 1/13/20) ?
    State_Met(LCHNK)%SUNCOS    (1,:) = CSZA(:)
    State_Met(LCHNK)%SUNCOSmid (1,:) = CSZA(:)

    ! Field      : U10M, V10M
    ! Description: E/W and N/S wind speed @ 10m height
    ! Unit       : m/s
    ! Dimensions : nX, nY
    State_Met(LCHNK)%U10M      (1,:) = State%U(:,nZ)
    State_Met(LCHNK)%V10M      (1,:) = State%V(:,nZ)

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
    State_Met(LCHNK)%Z0        (1,:) = Z0(:)

    ! Field      : IODIDE
    ! Description: Surface iodide concentration
    ! Unit       : nM
    ! Dimensions : nX, nY
    fldname_ns = 'HCO_iodide'
    tmpIdx = pbuf_get_index(fldname_ns, RC)
    IF ( tmpIdx < 0 ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
       State_Met(LCHNK)%IODIDE(1,:)   = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
       State_Met(LCHNK)%IODIDE(1,:) = pbuf_ik(:,nZ)
    ENDIF

    ! Field      : SALINITY
    ! Description: Ocean salinity
    ! Unit       : PSU
    ! Dimensions : nX, nY
    ! Note       : Possibly get ocean salinity from POP?
    fldname_ns = 'HCO_salinity'
    tmpIdx = pbuf_get_index(fldname_ns, RC)
    IF ( tmpIdx < 0 ) THEN
       IF ( rootChunk ) Write(iulog,*) "chem_timestep_tend: Field not found ", TRIM(fldname_ns)
       State_Met(LCHNK)%SALINITY(1,:) = 0.0e+0_fp
    ELSE
       pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
       CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)
       State_Met(LCHNK)%SALINITY(1,:) = pbuf_ik(:,nZ)
    ENDIF

    ! Three-dimensional fields on level edges
    DO J = 1, nY
    DO L = 1, nZ+1
       ! Field      : PEDGE
       ! Description: Wet air pressure at (vertical) level edges
       ! Unit       : hPa
       ! Dimensions : nX, nY, nZ+1
       State_Met(LCHNK)%PEDGE   (1,J,L) = State%Pint(J,nZ+2-L)*0.01e+0_fp

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
       State_Met(LCHNK)%PFILSAN (1,J,L) = LsFlxSnw(j,nZ+2-L) ! kg/m2/s
       State_Met(LCHNK)%PFLLSAN (1,J,L) = MAX(0.0e+0_fp,LsFlxPrc(J,nZ+2-L) - LsFlxSnw(J,nZ+2-L)) ! kg/m2/s
    ENDDO
    ENDDO

    DO J = 1, nY
       ! Field      : U, V
       ! Description: Max cloud top height
       ! Unit       : level
       ! Dimensions : nX, nY
       State_Met(LCHNK)%cldTops(1,J) = nZ + 1 - NINT(cldTop(J))
    ENDDO

    ! Three-dimensional fields on level centers
    DO J = 1, nY
    DO L = 1, nZ
       ! Field      : U, V
       ! Description: E/W and N/S component of wind
       ! Unit       : m/s
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%U        (1,J,L) = State%U(J,nZ+1-L)
       State_Met(LCHNK)%V        (1,J,L) = State%V(J,nZ+1-L)

       ! Field      : OMEGA
       ! Description: Updraft velocity
       ! Unit       : Pa/s
       ! Dimensions : nX, nY, nZ
       !State_Met(LCHNK)%OMEGA    (1,J,L) = State%Omega(J,nZ+1-L)

       ! Field      : CLDF
       ! Description: 3-D cloud fraction
       ! Unit       : -
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%CLDF     (1,J,L) = cldFrc(j,nZ+1-l)

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
       State_Met(LCHNK)%QI       (1,J,L) = MAX(1.0e-05_fp, State%Q(J,nZ+1-L,ixCldIce)) ! kg ice / kg dry air
       State_Met(LCHNK)%QL       (1,J,L) = MAX(1.0e-05_fp, State%Q(J,nZ+1-L,ixCldLiq)) ! kg water / kg dry air

       ! Field      : RH
       ! Description: Relative humidity
       ! Unit       : %
       ! Dimensions : nX, nY, nZ
       State_Met(LCHNK)%RH       (1,J,L) = RelHum(J,nZ+1-L)  * 100.0e+0_fp

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
       ! Unit       : kg
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
       State_Met(LCHNK)%TMPU1    (1,J,L) = State%T(J,nZ+1-L)
       State_Met(LCHNK)%TMPU2    (1,J,L) = State%T(J,nZ+1-L)
    ENDDO
    ENDDO

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

    ! Nullify all pointers
    Nullify(PblH    )
    Nullify(Fsds    )
    Nullify(PRain   )
    Nullify(LsFlxSnw)
    Nullify(LsFlxPrc)
    Nullify(cldTop  )
    Nullify(cldFrc  )
    Nullify(NEvapr  )
    Nullify(RprdTot )

    ! Field      : InChemGrid
    ! Description: Are we in the chemistry grid?
    ! Unit       : -
    ! Dimensions : nX, nY, nZ
    State_Met(LCHNK)%InChemGrid(:,:,:) = .True.

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

    CALL Accept_External_PEdge( State_Met = State_Met(LCHNK), &
                                RC        = RC               )

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
    State_Met(LCHNK)%PS1_DRY (1,:) = State%PSDry(:) * 0.01e+0_fp
    State_Met(LCHNK)%PS2_DRY (1,:) = State%PSDry(:) * 0.01e+0_fp

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
    !  (9)  MAIRDEN   : Mean grid box moist air density          [kg/m^3]
    !                   (defined as total moist air mass/box vol)
    !  (10) AD        : Total dry air mass in grid box             [kg]
    !  (11) ADMOIST   : Total moist air mass in grid box           [kg]
    !  (12) BXHEIGHT  : Vertical height of grid box                [m]
    !  (13) AIRVOL    : Volume of grid box                         [m^3]
    !  (14) MOISTMW   : Molecular weight of moist air in box     [g/mol]
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


    ! Initialize strat chem if not already done. This has to be done here because
    ! it needs to have non-zero values in State_Chm%AD, which only happens after
    ! the first call to AirQnt
    !IF ( (.not.SCHEM_READY) .and. Input_Opt%LSCHEM ) THEN
    IF ( (.not.SCHEM_READY) .and. .True. ) THEN !TMMF
        CALL Init_Strat_Chem( Input_Opt  = Input_Opt,         &
                              State_Chm  = State_Chm(LCHNK),  &
                              State_Met  = State_Met(LCHNK),  &
                              State_Grid = State_Grid(LCHNK), &
                              RC         = RC                )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Init_Strat_Chem"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
        SCHEM_READY = .True.
    ENDIF

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
    !IF ( LCHEM .and. newMonth ) THEN
    !
    !   ! The following only apply when photolysis is used,
    !   ! that is for fullchem or aerosol simulations.
    !   IF ( ITS_A_FULLCHEM_SIM  .or. ITS_AN_AEROSOL_SIM ) THEN
    !
    !      ! Copy UV Albedo data (for photolysis) into the
    !      ! State_Met%UVALBEDO field. (bmy, 3/20/15)
    !      CALL Get_UvAlbedo( Input_Opt = Input_Opt,        &
    !                         State_Met = State_Met(LCHNK), &
    !                         RC        = RC               )
    !
    !      ! Trap potential errors
    !      IF ( RC /= GC_SUCCESS ) THEN
    !         ErrMsg = 'Error encountered in "Get_UvAlbedo"!'
    !         CALL Error_Stop( ErrMsg, ThisLoc )
    !      ENDIF
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
    ! If using the CLM velocities, two options show up:
    ! 1. Compute dry deposition velocities over ocean and ice similarly
    !    to the way MOZART does it (OCNDDVEL_MOZART)
    ! 2. Use GEOS-Chem's dry deposition module to compute velocities
    !    and then scale them with the ocean fraction (OCNDDVEL_GEOSCHEM)
    !
    ! A third option would be to let GEOS-Chem compute dry deposition
    ! velocity (ALLDDVEL_GEOSCHEM), thus overwriting the input from CLM
    !
    ! drydep_method must be set to DD_XLND.
    !
    ! The following options are currently supported:
    ! - ALLDDVEL_GEOSCHEM
    ! - OCNDDVEL_GEOSCHEM
    ! - OCNDDVEL_MOZART
    !
    ! The ALLDDVEL_GEOSCHEM coupled with LANDTYPE_CLM requires that CLM
    ! passes land type information (land type and leaf area index).
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

#if   ( ALLDDVEL_GEOSCHEM || OCNDDVEL_GEOSCHEM )

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

#if   ( OCNDDVEL_GEOSCHEM )

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

#endif

#elif ( OCNDDVEL_MOZART )
       ! This routine updates the deposition velocities from CLM in the
       ! pointer lnd(LCHNK)%dvel as long as drydep_method == DD_XLND is
       ! True.
       CALL drydep_update( State, cam_in )

       windSpeed(:nY) = SQRT( State%U(:nY,nZ)*State%U(:nY,nZ) + &
                              State%V(:nY,nZ)*State%V(:nY,nZ)  )
       potT(:nY)      = State%T(:nY,nZ) * (1._fp + qH2O(:nY,nZ))

       CALL get_lat_all_p( LCHNK, nY, latndx )
       CALL get_lon_all_p( LCHNK, nY, lonndx )

       CALL drydep_fromlnd( ocnfrac      = cam_in%ocnfrac(:),             &
                            icefrac      = cam_in%icefrac(:),             &
                            ncdate       = currYMD,                       &
                            sfc_temp     = cam_in%TS(:),                  &
                            pressure_sfc = State%PS(:),                   &
                            wind_speed   = windSpeed(:),                  &
                            spec_hum     = qH2O(:,nZ),                    &
                            air_temp     = State%T(:,nZ),                 &
                            pressure_10m = State%PMid(:,nZ),              &
                            rain         = State_Met(LCHNK)%PRECTOT(1,:), &
                            snow         = cam_in%Snowhland(:),           &
                            solar_flux   = State_Met(LCHNK)%SWGDN(1,:),   &
                            dvelocity    = MOZART_depVel(:,:),            &
                            dflx         = MOZART_depFlx(:,:),            &
                            State_Chm    = State_Chm(LCHNK),              &
                            tv           = potT(:),                       &
                            soilw        = -99._fp,                       &
                            rh           = relHum(:,nZ),                  &
                            ncol         = nY,                            &
                            lonndx       = lonndx(:),                     &
                            latndx       = latndx(:),                     &
                            lchnk        = LCHNK                         )

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
          !    Write(iulog,*) "CLM-depVel    = ", &
          !  MAXVAL(cam_in%depvel(:nY,N)) * 1.0e-02_fp, " [m/s]", LCHNK
          !    IF ( drySpc_ndx(N) > 0 ) THEN
          !        Write(iulog,*) "Merged depVel = ", &
          !  MAXVAL(MOZART_depVel(:nY,drySpc_ndx(N))) * 1.0e-02_fp, " [m/s]", LCHNK
          !    ENDIF
          !ENDIF

          IF ( ( map2GC_dryDep(N) > 0 ) .AND. ( drySpc_ndx(N) > 0 ) ) THEN
              ! State_Chm%DryDepVel is in m/s
              State_Chm(LCHNK)%DryDepVel(1,:nY,map2GC_dryDep(N)) = &
                 MOZART_depVel(:nY,drySpc_ndx(N)) * 1.0e-02_fp
          ENDIF

       ENDDO

#else
       ! We should be in one of the cases above as any exceptions should be
       ! caught when running chem_readnl, but just for safety's safe:
       CALL ENDRUN('Incorrect definitions for dry deposition velocities')
#endif

       CALL Update_DryDepSav( Input_Opt  = Input_Opt,         &
                              State_Chm  = State_Chm(LCHNK),  &
                              State_Diag = State_Diag(LCHNK), &
                              State_Grid = State_Grid(LCHNK), &
                              State_Met  = State_Met(LCHNK),  &
                              RC         = RC                )

    ENDIF

    !===========================================================
    !      ***** M I X E D   L A Y E R   M I X I N G *****
    !===========================================================

    ! Note: mixing routine expects tracers in v/v
    ! DO_MIXING applies the tracer tendencies (dry deposition,
    ! emission rates) to the tracer arrays and performs PBL
    ! mixing.
    ! In the non-local PBL scheme, dry deposition and emission
    ! fluxes below the PBL are handled within the PBL mixing
    ! routine. Otherwise, tracer concentrations are first updated
    ! and the full-mixing is then applied.
    ! (ckeller, 3/5/15)
    ! NOTE: Tracer concentration units are converted locally
    ! to [v/v dry air] for mixing. Eventually mixing should
    ! be updated to use [kg/kg total air] (ewl, 9/18/15)
    !
    ! This requires HEMCO. For now comment out.
    ! Thibaud M. Fritz - 05/07/20
    !CALL Do_Mixing( Input_Opt  = Input_Opt,         &
    !                State_Chm  = State_Chm(LCHNK),  &
    !                State_Diag = State_Diag(LCHNK), &
    !                State_Grid = State_Grid(LCHNK), &
    !                State_Met  = State_Met(LCHNK),  &
    !                RC         = RC                )
    !
    !! Trap potential errors
    !IF ( RC /= GC_SUCCESS ) THEN
    !   ErrMsg = 'Error encountered in "Do_Mixing"!'
    !   CALL Error_Stop( ErrMsg, ThisLoc )
    !ENDIF

    !!===========================================================
    !!        ***** C L O U D   C O N V E C T I O N *****
    !!===========================================================
    !IF ( LCONV ) THEN
    !
    !   ! Call the appropriate convection routine
    !   ! NOTE: Tracer concentration units are converted locally
    !   ! to [kg/kg total air] for convection (ewl, 9/18/15)
    !   CALL Do_Convection( Input_Opt  = Input_Opt,         &
    !                       State_Chm  = State_Chm(LCHNK),  &
    !                       State_Diag = State_Diag(LCHNK), &
    !                       State_Grid = State_Grid(LCHNK), &
    !                       State_Met  = State_Met(LCHNK),  &
    !                       RC         = RC                )
    !
    !   ! Trap potential errors
    !   IF ( RC /= GC_SUCCESS ) THEN
    !      ErrMsg = 'Error encountered in "Do_Convection"!'
    !      CALL Error_Stop( ErrMsg, ThisLoc )
    !   ENDIF
    !ENDIF

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

        ! Do wet deposition
        ! NOTE: Tracer concentration units are converted locally
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

    ENDIF

    ! Make sure State_Chm(lchnk) is back in kg/kg dry!
    ! Reset H2O MMR to the initial value (no chemistry tendency in H2O just yet)
    State_Chm(LCHNK)%Species(1,:,:,iH2O) = MMR_Beg(:,:,iH2O)

    ! Store unadvected species data
    SlsData = 0.0e+0_r8
    DO N = 1, nSls
       M = map2GC_Sls(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO K = 1, nZ
             SlsData(J,nZ+1-K,N) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8)
          ENDDO
          ENDDO
       ENDIF
    ENDDO
    CALL set_short_lived_species( SlsData, LCHNK, nY, pbuf )

    ! Write diagnostic output
    DO N = 1, pcnst
       M = map2GC(N)
       I = map2Idx(N)
       IF ( M > 0 ) THEN
          SpcName = tracerNames(I)
          VMR     = 0.0e+0_r8
          DO J = 1, nY
          DO K = 1, nZ
             VMR(J,nZ+1-K) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8) * MWRatio(I)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), VMR(:nY,:), nY, LCHNK )
       ENDIF
    ENDDO

    DO N = 1, nSls
       SpcName = slsNames(n)
       VMR = 0.0e+0_r8
       M = map2GC_Sls(n)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO K = 1, nZ
             VMR(J,nZ+1-K) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8) * SLSMWratio(N)
          ENDDO
          ENDDO
          CALL OutFld( TRIM(SpcName), VMR(:nY,:), nY, LCHNK )
       ENDIF
    ENDDO

    ! NOTE: Re-flip all the arrays vertically or suffer the consequences
    ! ptend%q dimensions: [column, ?, species]
    ptend%q(:,:,:) = 0.0e+0_r8
    MMR_End = 0.0e+0_r8
    DO N = 1, pcnst
       M = map2GC(N)
       IF ( M > 0 ) THEN
          DO J = 1, nY
          DO K = 1, nZ
             MMR_End (J,K,M) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8)
             MMR_TEnd(J,K,M) = MMR_End(J,K,M) - MMR_Beg(J,K,M)
             ptend%q(J,nZ+1-K,N) = (MMR_End(J,K,M)-MMR_Beg(J,K,M))/dT
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    ! Debug statements
    ! Ozone tendencies
    IF ( rootChunk ) THEN
       Write(iulog,*) " MMR_Beg = ", MMR_Beg(1,:,2)
       Write(iulog,*) " MMR_End = ", MMR_End(1,:,2)
    ENDIF

    IF (PRESENT(fh2o)) THEN
       fh2o(:nY) = 0.0e+0_r8
       !DO K = 1, nZ
       !   fh2o(:nY) = fh2o(:nY) + Ptend%Q(:nY,K,iH2O)*State%Pdel(:nY,K)/Gravit
       !ENDDO
    ENDIF

    IF (rootChunk) WRITE(iulog,*) ' GEOS-Chem Chemistry step ', iStep, ' completed'

  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    CHARACTER(LEN=*), INTENT(IN)  :: name       !  constituent name
    REAL(r8),         INTENT(IN)  :: latvals(:) ! lat in degrees (NCOL)
    REAL(r8),         INTENT(IN)  :: lonvals(:) ! lon in degrees (NCOL)
    LOGICAL,          INTENT(IN)  :: mask(:)    ! Only initialize where .true.
    REAL(r8),         INTENT(OUT) :: q(:,:)     ! kg tracer/kg dry air (NCOL, PVER
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    INTEGER  :: iLev, NLEV, I
    REAL(r8) :: QTemp, Min_MMR

    NLEV = SIZE(q, 2)
    ! Retrieve a "background value" for this from the database
    Min_MMR = 1.0e-38_r8
    DO I = 1, nTracers
       IF (TRIM(tracerNames(I)).eq.TRIM(name)) THEN
          Min_MMR = ref_MMR(I)
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

    use GC_Emissions_Mod, only: GC_Emissions_Final

    INTEGER :: I, RC

    ! Finalize GEOS-Chem
    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_FINAL'

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

    CALL GC_Emissions_Final

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
    IF (ALLOCATED(State_Chm))     DEALLOCATE(State_Chm)
    IF (ALLOCATED(State_Diag))    DEALLOCATE(State_Diag)
    IF (ALLOCATED(State_Grid))    DEALLOCATE(State_Grid)
    IF (ALLOCATED(State_Met))     DEALLOCATE(State_Met)

    IF (ALLOCATED(slvd_Lst    ))  DEALLOCATE(slvd_Lst)
    IF (ALLOCATED(slvd_ref_MMR))  DEALLOCATE(slvd_ref_MMR)

    RETURN

  end subroutine chem_final
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    TYPE(file_desc_t) :: File

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT_RESTART'

    RETURN

  end subroutine chem_init_restart
!===============================================================================
  subroutine chem_write_restart( File )
    !use tracer_cnst, only: write_tracer_cnst_restart
    !use tracer_srcs, only: write_tracer_srcs_restart
    !use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t
    IMPLICIT NONE
    TYPE(file_desc_t) :: File

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_WRITE_RESTART'
    !
    ! data for offline tracers
    !
    !call write_tracer_cnst_restart(File)
    !call write_tracer_srcs_restart(File)
    !call write_linoz_data_restart(File)
  end subroutine chem_write_restart
!===============================================================================
  subroutine chem_read_restart( File )
    !use tracer_cnst, only: read_tracer_cnst_restart
    !use tracer_srcs, only: read_tracer_srcs_restart
    !use linoz_data,  only: read_linoz_data_restart

    use pio, only : file_desc_t
    IMPLICIT NONE
    TYPE(file_desc_t) :: File

    if (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_READ_RESTART'
    !
    ! data for offline tracers
    !
    !call read_tracer_cnst_restart(File)
    !call read_tracer_srcs_restart(File)
    !call read_linoz_data_restart(File)
  end subroutine chem_read_restart
!================================================================================
  subroutine chem_emissions( state, cam_in )
    use camsrfexch,       only: cam_in_t
    use physics_buffer, only : pbuf_get_chunk, pbuf_get_field, pbuf_get_index
    use ppgrid,         only : pver ! for vertical

    use constituents,      only: cnst_name
    use mo_chem_utls,        only : get_spc_ndx
    use chem_mods, only : tracerNames, nTracers
    use PhysConstants,    only: PI, PI_180

    ! Arguments:

    TYPE(physics_state),    INTENT(IN)    :: state   ! Physics state variables
    TYPE(cam_in_t),         INTENT(INOUT) :: cam_in  ! import state

    INTEGER :: M, N, I
    INTEGER :: LCHNK, NCOL
    LOGICAL :: rootChunk

    REAL(r8) :: sflx(State%NCOL,nTracers)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:) ! slice of pbuf in chnk
    real(r8), pointer     :: pbuf_ik(:,:)          ! ptr to pbuf data (/pcols,pver/)
    integer               :: tmpIdx                ! pbuf field id
    character(len=255)    :: fldname_ns            ! field name HCO_NH3
    integer               :: RC                    ! return code

    ! LCHNK: which chunk we have on this process
    LCHNK = State%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    NCOL  = State%NCOL
    rootChunk = ( MasterProc.and.(LCHNK.EQ.BEGCHUNK) )

    sflx(:,:) = 0.0e+0_r8

    DO N = 1, nTracers

       fldname_ns = 'HCO_' // TRIM(tracerNames(N))
       tmpIdx = pbuf_get_index(fldname_ns, RC)
       IF ( tmpIdx < 0 ) THEN
          IF ( rootChunk ) Write(iulog,*) "chem_emissions hemco: Field not found ", TRIM(fldname_ns)
       ELSE
          ! This is already in chunk, retrieve it
          pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, LCHNK)
          CALL pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)

          IF ( .NOT. ASSOCIATED(pbuf_ik) ) THEN ! Sanity check
             CALL ENDRUN("chem_emissions: FATAL - tmpIdx > 0 but pbuf_ik not associated")
          ENDIF

          ! For each column retrieve data from pbuf_ik(I,K)
          sflx(1:ncol,N) = pbuf_ik(1:ncol,pver) ! Only surface emissions for now,

          !TMMF, Replace this in chem_init
          M = -1
          DO I = 1, nTracers
             IF( TRIM( to_upper( TRIM(tracerNames(N)) ) ) == TRIM( to_upper( cnst_name(I)) ) ) THEN
                M = I
                EXIT
             ENDIF
          ENDDO

          IF ( M <= 0 ) CYCLE

          cam_in%cflx(1:ncol,M) = sflx(1:ncol,N)
          IF ( rootChunk .and. ( MAXVAL(sflx(1:ncol,N)) > 0.0e+0_fp ) ) THEN
             Write(iulog,'(a,a,E16.4,a,a)') "chem_emissions: debug added emiss for ", &
                TRIM(cnst_name(M)), MAXVAL(sflx(1:ncol,N)), " from ", TRIM(fldname_ns)
             Write(iulog,'(a,a,a,E16.4)') "chem_emissions: Total emission flux for ", &
                TRIM(cnst_name(M)), " is: ", MAXVAL(cam_in%cflx(1:ncol,M))
          ENDIF
       ENDIF
    ENDDO

  end subroutine chem_emissions

end module chemistry
