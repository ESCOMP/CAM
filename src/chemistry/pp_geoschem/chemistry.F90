!================================================================================================
! This is the "GEOS-Chem" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use physics_types,       only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,              only: begchunk, endchunk, pcols
  use ppgrid,              only: pver
  use constituents,        only: pcnst, cnst_add
  !use mo_gas_phase_chemdr, only: map2chm
  !use mo_constants,        only: pi
  use shr_const_mod,       only: molw_dryair=>SHR_CONST_MWDAIR
  !use mo_chem_utls,        only : get_spc_ndx
  !use chem_mods,           only : gas_pcnst, adv_mass
  !use mo_sim_dat, only: set_sim_dat
  use spmd_utils,          only : MasterProc, myCPU=>Iam, nCPUs=>npes
  use cam_logfile,         only : iulog

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

  ! Private data
  !===== SDE DEBUG =====
  integer, parameter :: ntracersmax = 200    ! Must be equal to nadv_chem
  integer            :: ntracers
  character(len=255) :: tracernames(ntracersmax)
  integer            :: indices(ntracersmax)
  real(r8)           :: adv_mass(ntracersmax)

  ! Short-lived species (i.e. not advected)
  integer, parameter :: nslsmax = 500        ! UNadvected species only
  integer            :: nsls
  character(len=255) :: slsnames(nslsmax)
  !===== SDE DEBUG =====

  ! Location of valid input.geos
  CHARACTER(LEN=500) :: inputGeosPath

  ! Location of chemistry input (for now)
  CHARACTER(LEN=500) :: chemInputsDir

  !-----------------------------
  ! Derived type objects
  !-----------------------------
  TYPE(OptInput)                  :: Input_Opt      ! Input Options object
  TYPE(ChmState),ALLOCATABLE      :: State_Chm(:)   ! Chemistry State object
  TYPE(DgnState),ALLOCATABLE      :: State_Diag(:)  ! Diagnostics State object
  TYPE(GrdState),ALLOCATABLE      :: State_Grid(:)  ! Grid State object
  TYPE(MetState),ALLOCATABLE      :: State_Met(:)   ! Meteorology State object
  TYPE(DgnList )                  :: Diag_List      ! Diagnostics list object

  ! Integer
  INTEGER                         :: RC

  ! Strings
  CHARACTER(LEN=255)              :: ThisLoc
  CHARACTER(LEN=255)              :: ErrMsg

!================================================================================================
contains
!================================================================================================

  LOGICAL function chem_is (NAME)

    CHARACTER(LEN=*), INTENT(IN) :: NAME

    chem_is = .false.
    IF (NAME == 'geoschem' ) THEN
       chem_is = .true.
    ENDIF
    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_IS'

  end function chem_is

!================================================================================================

  subroutine chem_register

    use physics_buffer, only : pbuf_add_field, dtype_r8
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents for chemistry
    !
    !-----------------------------------------------------------------------
    INTEGER            :: i, n
    REAL(r8)           :: cptmp
    REAL(r8)           :: qmin
    CHARACTER(LEN=128) :: mixtype
    CHARACTER(LEN=128) :: molectype
    CHARACTER(LEN=128) :: lng_name
    LOGICAL            :: camout
    LOGICAL            :: ic_from_cam2
    LOGICAL            :: has_fixed_ubc
    LOGICAL            :: has_fixed_ubflx
    ! SDE 2018-05-02: This seems to get called before anything else
    ! That includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    if (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_REGISTER'
    ! At the moment, we force nadv_chem=200 in the setup file
    DO I = 1, NTRACERSMAX
       ! TODO: Read input.geos in chem_readnl to get tracernames(1:ntracers)
       ! TODO: Get all other species properties here from species database
       ! Hard-code for now
       SELECT CASE (TRACERNAMES(I))
          CASE ('BCPI')
             lng_name = 'Hydrophilic black carbon'
             ! Molar mass (g/mol)
             adv_mass(i) = 1000.0e+0_r8 * (0.012e+0_r8)
          CASE ('OCS')
             lng_name = 'Carbonyl sulfide'
             ! Molar mass (g/mol)
             adv_mass(i) = 1000.0e+0_r8 * (0.060e+0_r8)
          CASE DEFAULT
             lng_name = tracernames(i)
             adv_mass(i) = 1000.0e+0_r8 * (0.001e+0_r8)
       END SELECT
       ! dummy value for specific heat of constant pressure (Cp)
       cptmp = 666._r8
       ! minimum mixing ratio
       qmin = 1.e-36_r8
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
       !write(tracernames(i),'(a,I0.4)') 'GCTRC_', i
       ! NOTE: In MOZART, this only gets called for tracers
       ! This is the call to add a "constituent"
       CALL cnst_add( trim(tracernames(i)), adv_mass(i), cptmp, qmin, n, &
                      readiv=ic_from_cam2, mixtype=mixtype, cam_outfld=camout, &
                      molectype=molectype, fixed_ubc=has_fixed_ubc, &
                      fixed_ubflx=has_fixed_ubflx, longname=trim(lng_name) )
    ENDDO

       ! MOZART uses this for short-lived species. Not certain exactly what it
       ! does, but note that the "ShortLivedSpecies" physics buffer already
       ! needs to have been initialized, which we haven't done. Physics buffers
       ! are fields which are available either across timesteps or for use to
       ! modules outside of chemistry
       ! More information:
       ! http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/node5.html
       !call pbuf_add_field('ShortLivedSpecies','global',dtype_r8,(/pcols,pver,nslvd/),pbf_idx)
       ! returned values
       !  n : mapping in CAM
       ! map2chm is a mozart variable
       !map2chm(n) = i
       !indices(i) = 0
       ! ===== SDE DEBUG =====

  end subroutine chem_register

  subroutine chem_readnl(nlfile)
    ! This is the FIRST routine to get called - so it should read in
    ! GEOS-Chem options from input.geos without actually doing any
    ! initialization

    use cam_abortutils, only : endrun
    use units,          only : getunit, freeunit
    use mpishorthand
    use gckpp_Model,    only : nspec, spc_names

    ! args
    CHARACTER(LEN=*), INTENT(IN) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    INTEGER :: I, UNITN, IERR
    CHARACTER(LEN=500) :: LINE
    logical :: menuFound, validSLS

    ! Set paths
    inputGeosPath='/n/scratchlfs/jacob_lab/elundgren/UT/runs/4x5_standard/input.geos.template'
    chemInputsDir='/n/holylfs/EXTERNAL_REPOS/GEOS-CHEM/gcgrid/gcdata/ExtData/CHEM_INPUTS/'

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_READNL'

    ! TODO: Read in input.geos and get species names
    IF (MasterProc) THEN
        UNITN = GETUNIT()
        OPEN( UNITN, FILE=TRIM(inputGeosPath), STATUS='old', IOSTAT=IERR )
        IF (IERR .NE. 0) THEN
            CALL ENDRUN('chem_readnl: ERROR opening input.geos')
        ENDIF

        ! Go to ADVECTED SPECIES MENU
        menuFound = .False.
        DO WHILE (.NOT.menuFound)
            READ( UNITN, '(a)', IOSTAT=IERR ) LINE
            IF (IERR.NE.0) THEN
                CALL ENDRUN('chem_readnl: ERROR finding advected species menu')
            ELSEIF (INDEX(LINE,'ADVECTED SPECIES MENU') > 0) then
                menuFound=.True.
            ENDIF
        ENDDO

        ! Skip first line
        READ(UNITN,'(a)',IOSTAT=IERR) LINE
        ! Read in tracer count
        READ(UNITN,'(26x,I)',IOSTAT=IERR) NTRACERS
        ! Skip divider line
        READ(UNITN,'(a)',IOSTAT=IERR) LINE
        ! Read in each tracer
        DO I=1,NTRACERS
            READ(UNITN,'(26x,a)',iostat=ierr) line
            tracernames(i) = trim(line)
        ENDDO
        CLOSE(UNITN)
        CALL FREEUNIT(UNITN)

        ! Assign remaining tracers dummy names
        DO I=(NTRACERS+1),NTRACERSMAX
            WRITE(TRACERNAMES(I),'(a,I0.4)') 'GCTRC_',I
        ENDDO

        ! Now go through the KPP mechanism and add any species not implemented by
        ! the tracer list in input.geos
        IF ( NSPEC > NSLSMAX ) THEN
            CALL ENDRUN('chem_readnl: too many species - increase nslsmax')
        ENDIF

        NSLS = 0
        DO I=1,NSPEC
            ! Get the name of the species from KPP
            LINE = ADJUSTL(TRIM(SPC_NAMES(I)))
            ! Only add this
            validSLS = ( (.NOT.ANY(TRIM(LINE).EQ.TRACERNAMES)).AND.&
                         (.NOT.(LINE(1:2) == 'RR')) )
            IF (validSLS) THEN
                ! Genuine new short-lived species
                NSLS = NSLS + 1
                SLSNAMES(NSLS) = TRIM(LINE)
                WRITE(iulog,'(a,I5,a,a)') ' --> GC species ',nsls, ': ', TRIM(LINE)
            ENDIF
        ENDDO
    ENDIF

    ! Broadcast to all processors
#if defined( SPMD )
    CALL MPIBCAST(NTRACERS,    1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST(TRACERNAMES, LEN(TRACERNAMES(1))*NTRACERSMAX, MPICHAR, 0, MPICOM )
    CALL MPIBCAST(NSLS,        1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST(SLSNAMES,    LEN(SLSNAMES(1))*NSLSMAX,        MPICHAR, 0, MPICOM )
#endif

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .true.

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_IS_ACTIVE'

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

    DO I = 1, NTRACERS
       IF (TRIM(TRACERNAMES(I)) .eq. TRIM(NAME)) THEN
          chem_implements_cnst = .true.
          EXIT
       ENDIF
    ENDDO

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_IMPLEMENTS_CNST'

  end function chem_implements_cnst

!===============================================================================

  subroutine chem_init(phys_state, pbuf2d)
    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize GEOS-Chem parts (state objects, mainly)
    !          (and declare history variables)
    !
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc
    use cam_history,    only: addfld, add_default, horiz_only

    use mpishorthand
    use cam_abortutils, only : endrun

    use Input_Opt_Mod
    use State_Chm_Mod
    use State_Grid_Mod
    use State_Met_Mod
    use DiagList_Mod,   only : Init_DiagList, Print_DiagList
    use GC_Environment_Mod
    use GC_Grid_Mod,    only : SetGridFromCtrEdges


    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,  only : PI, PI_180
    use PhysConstants,  only : Re

    use Phys_Grid,      only : get_Area_All_p
    use hycoef,         only : ps0, hyai, hybi

    use Time_Mod,      only : Accept_External_Date_Time
    !use Time_Mod,      only : Set_Begin_Time,   Set_End_Time
    !use Time_Mod,      only : Set_Current_Time, Set_DiagB
    !use Transfer_Mod,  only : Init_Transfer
    use Linoz_Mod,     only : Linoz_Read

    use CMN_Size_Mod

    use Drydep_Mod,    only : Init_Drydep
    use Carbon_Mod,    only : Init_Carbon
    use Dust_Mod,      only : Init_Dust
    use Seasalt_Mod,   only : Init_Seasalt
    use Sulfate_Mod,   only : Init_Sulfate
    use Aerosol_Mod,   only : Init_Aerosol
    use WetScav_Mod,   only : Init_WetScav
    use TOMS_Mod,      only : Init_TOMS
    use Pressure_Mod,  only : Init_Pressure, Accept_External_ApBp

    TYPE(physics_state), INTENT(IN):: phys_state(BEGCHUNK:ENDCHUNK)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    ! Local variables

    !----------------------------
    ! Scalars
    !----------------------------

    ! Integers
    INTEGER :: LCHNK(BEGCHUNK:ENDCHUNK), NCOL(BEGCHUNK:ENDCHUNK)
    INTEGER               :: IWAIT, IERR
    INTEGER               :: NX, NY, NZ
    INTEGER               :: IX, IY, IZ
    INTEGER               :: NLEV, I, J, L, RC
    INTEGER               :: NLINOZ

    ! Logicals
    LOGICAL               :: am_I_Root
    LOGICAL               :: prtDebug

    ! Strings
    CHARACTER(LEN=255)       :: historyConfigFile


    ! Grid setup
    REAL(fp)              :: lonVal,  latVal
    REAL(fp)              :: dLonFix, dLatFix
    REAL(f4), ALLOCATABLE :: lonMidArr(:,:),  latMidArr(:,:)
    REAL(f4), ALLOCATABLE :: lonEdgeArr(:,:), latEdgeArr(:,:)
    REAL(r8), ALLOCATABLE :: linozData(:,:,:,:)

    REAL(r8), ALLOCATABLE :: Col_Area(:)
    REAL(fp), ALLOCATABLE :: Ap_CAM_Flip(:), Bp_CAM_Flip(:)

    ! Assume a successful return until otherwise
    RC                      = GC_SUCCESS

    ! For error trapping
    ErrMsg                  = ''
    ThisLoc                 = ' -> at GEOS-Chem (in chemistry/pp_geoschem/chemistry.F90)'

    ! LCHNK: which chunks we have on this process
    LCHNK = PHYS_STATE%LCHNK
    ! NCOL: number of atmospheric columns for each chunk
    NCOL  = PHYS_STATE%NCOL
    ! NLEV: number of vertical levels
    NLEV  = PVER

    write(iulog,'(2(a,x,I6,x))') 'chem_init called on PE ', myCPU, ' of ', nCPUs

    ! The GEOS-Chem grids on every "chunk" will all be the same size, to avoid
    ! the possibility of having differently-sized chunks
    NX = 1
    NY = MAXVAL(NCOL)
    NZ = NLEV

    ! This ensures that each process allocates everything needed for its chunks
    ALLOCATE(State_Chm(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Chm')
    ALLOCATE(State_Diag(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Diag')
    ALLOCATE(State_Grid(BEGCHUNK:ENDCHUNK), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Grid')
    ALLOCATE(State_Met(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Met')

        ! Set some basic flags
    Input_Opt%Max_BPCH_Diag     = 1000
    Input_Opt%Max_AdvectSpc     = 500
    Input_Opt%Max_Families      = 250

    Input_Opt%Linoz_NLevels     = 25
    Input_Opt%Linoz_NLat        = 18
    Input_Opt%Linoz_NMonths     = 12
    Input_Opt%Linoz_NFields     = 7
    Input_Opt%RootCPU           = MasterProc

    DO I = BEGCHUNK, ENDCHUNK

        ! Only treat the first chunk as the "root"
        am_I_Root = ((I.EQ.BEGCHUNK) .and. MasterProc)

        ! Initialize fields of the Grid State object
        CALL Init_State_Grid( am_I_Root  = am_I_Root,      &
                              State_Grid = State_Grid(I),  &
                              RC         = RC         )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered within call to "Init_Grid_State"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        State_Grid(I)%NX = NX
        State_Grid(I)%NY = NY
        State_Grid(I)%NZ = NZ
        ! Initialize GEOS-Chem horizontal grid structure
        CALL GC_Init_Grid( am_I_Root  = am_I_Root,      &
                           Input_Opt  = Input_Opt,      &
                           State_Grid = State_Grid(I),  &
                           RC         = RC          )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
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
    CALL GC_Allocate_All ( am_I_Root      = MasterProc,           &
                           Input_Opt      = Input_Opt,            &
                           State_Grid     = State_Grid(BEGCHUNK), &
                           value_I_Lo     = 1,                    &
                           value_J_Lo     = 1,                    &
                           value_I_Hi     = NX,                   &
                           value_J_Hi     = NY,                   &
                           value_IM       = NX,                   &
                           value_JM       = NY,                   &
                           value_LM       = NZ,                   &
                           value_IM_WORLD = NX,                   &
                           value_JM_WORLD = NY,                   &
                           value_LM_WORLD = NZ,                   &
                           value_LLSTRAT  = 59,                   & !TMMF
                           RC             = RC        )
        IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GC_Allocate_All"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF


    ! TODO: Mimic GEOS-Chem's reading of input options
    !IF (MasterProc) THEN
    !   CALL Read_Input_File( am_I_Root   = .True., &
    !                         Input_Opt   = Input_Opt(BEGCHUNK), &
    !                         srcFile     = inputGeosPath,      &
    !                         RC          = RC )
    !ENDIF
    !CALL <broadcast data to other CPUs>

    ! For now just hard-code it
    ! First setup directories
    Input_Opt%Chem_Inputs_Dir      = TRIM(chemInputsDir)

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
    Input_Opt%N_Advect               = NTracers
    IF (Input_Opt%N_Advect.GT.Input_Opt%Max_AdvectSpc) THEN
        CALL ENDRUN('Number of tracers exceeds max count')
    ENDIF
    ! Assign tracer names
    DO J = 1, Input_Opt%N_Advect
        Input_Opt%AdvectSpc_Name(J) = TRIM(TRACERNAMES(J))
    ENDDO
    ! No tagged species
    Input_Opt%LSplit = .False.

    ! Now READ_TRANSPORT_MENU
    Input_Opt%LTran                  = .True.
    Input_Opt%LFill                  = .True.
    Input_Opt%TPCore_IOrd            = 3
    Input_Opt%TPCore_JOrd            = 3
    Input_Opt%TPCore_KOrd            = 3

    ! Now READ_CONVECTION_MENU
    ! For now, TMMF
    Input_Opt%LConv                  = .False.
    Input_Opt%LTurb                  = .False.
    Input_Opt%LNLPBL                 = .False.

    ! Now READ_EMISSIONS_MENU
    Input_Opt%LEmis                  = .False.
    Input_Opt%HCOConfigFile          = 'HEMCO_Config.rc'
    Input_Opt%LFix_PBL_Bro           = .False.

    ! Set surface VMRs - turn this off so that CAM can handle it
    Input_Opt%LCH4Emis               = .False.
    Input_Opt%LCH4SBC                = .False.
    Input_Opt%LOCSEmis               = .False.
    Input_Opt%LCFCEmis               = .False.
    Input_Opt%LClEmis                = .False.
    Input_Opt%LBrEmis                = .False.
    Input_Opt%LN2OEmis               = .False.
    Input_Opt%LBasicEmis             = .False.

    ! Set initial conditions
    Input_Opt%LSetH2O                = .True.

    ! CFC control
    Input_Opt%CFCYear                = 0

    ! Now READ_AEROSOL_MENU
    Input_Opt%LSulf               = .True.
    Input_Opt%LMetalcatSO2        = .True.
    Input_Opt%LCarb               = .True.
    Input_Opt%LBrC                = .False.
    Input_Opt%LSOA                = .True.
    Input_Opt%LSVPOA              = .False.
    Input_Opt%LOMOC               = .False.
    Input_Opt%LDust               = .True.
    Input_Opt%LDstUp              = .False.
    Input_Opt%LSSalt              = .True.
    Input_Opt%SalA_rEdge_um(1)    = 0.01e+0_fp
    Input_Opt%SalA_rEdge_um(2)    = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(1)    = 0.50e+0_fp
    Input_Opt%SalC_rEdge_um(2)    = 8.00e+0_fp
    Input_Opt%LMPOA               = .False.
    Input_Opt%LGravStrat          = .True.
    Input_Opt%LSolidPSC           = .True.
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
    ! Disable dry/wet dep for now
    Input_Opt%LDryD                  = .False.
    Input_Opt%LWetD                  = .False.
    Input_Opt%CO2_Effect             = .False.
    Input_Opt%CO2_Level              = 390.0_fp
    Input_Opt%CO2_Ref                = 390.0_fp

    ! Now READ_CHEMISTRY_MENU
    Input_Opt%LChem                  = .True.
    Input_Opt%LSChem                 = .True.
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

    ! Read in data for Linoz. All CPUs allocate one array to hold the data. Only
    ! the root CPU reads in the data; then we copy it out to a temporary array,
    ! broadcast to all other CPUs, and finally duplicate the data into every
    ! copy of Input_Opt
    IF (Input_Opt%LLinoz) THEN
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
            CALL Linoz_Read( MasterProc, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Linoz_Read"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
            ! Copy the data to a temporary array
            linozData = REAL(Input_Opt%LINOZ_TPARM,r8)
        ENDIF
#if defined( SPMD )
        CALL MPIBCAST(linozData, nLinoz, MPIR8, 0, MPICOM )
#endif
        IF ( .NOT. MasterProc ) THEN
            Input_Opt%LINOZ_TPARM = REAL(linozData,fp)
        ENDIF
        DEALLOCATE(linozData)
    ENDIF


    ! Note: The following calculations do not setup the gridcell areas.
    !       In any case, we will need to be constantly updating this grid
    !       to compensate for the "multiple chunks per processor" element
    ALLOCATE(lonMidArr(NX,NY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonMidArr')
    ALLOCATE(lonEdgeArr(NX+1,NY+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating lonEdgeArr')
    ALLOCATE(latMidArr(NX,NY), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating latMidArr')
    ALLOCATE(latEdgeArr(NX+1,NY+1), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating latEdgeArr')

    ! We could try and get the data from CAM.. but the goal is to make this GC
    ! component completely grid independent. So for now, we set to arbitrary
    ! values
    ! TODO: This needs more refinement. For now, this generates identical
    ! State_Grid for all chunks
    DO L = BEGCHUNK, ENDCHUNK
        lonMidArr = 0.0e+0_f4
        latMidArr = 0.0e+0_f4
        dLonFix   = 360.0e+0_fp / REAL(NX,fp)
        dLatFix   = 180.0e+0_fp / REAL(NY,fp)
        DO I = 1,NX
            ! Center of box, assuming dateline edge
            lonVal = -180.0e+0_fp + (REAL(I-1,fp)*dLonFix)
            DO J = 1,NY
                ! Center of box, assuming regular cells
                latVal = -90.0e+0_fp + (REAL(J-1,fp)*dLatFix)
                lonMidArr(I,J)  = REAL((lonVal + (0.5e+0_fp * dLonFix)) * PI_180, f4)
                latMidArr(I,J)  = REAL((latVal + (0.5e+0_fp * dLatFix)) * PI_180, f4)

                ! Edges of box, assuming regular cells
                lonEdgeArr(I,J) = REAL(lonVal * PI_180, f4)
                latEdgeArr(I,J) = REAL(latVal * PI_180, f4)
            ENDDO
            ! Edges of box, assuming regular cells
            lonEdgeArr(I,NY+1)  = REAL((lonVal + dLonFix) * PI_180, f4)
            latEdgeArr(I,NY+1)  = REAL((latVal + dLatFix) * PI_180, f4)
        ENDDO
        DO J = 1,NY+1
            ! Edges of box, assuming regular cells
            latVal = -90.0e+0_fp + (REAL(J-1,fp)*dLatFix)
            lonEdgeArr(NX+1,J)  = REAL((lonVal + dLonFix) * PI_180, f4)
            latEdgeArr(NX+1,J)  = REAL((latVal) * PI_180, f4)
    ENDDO

        CALL SetGridFromCtrEdges( am_I_Root  = MasterProc,    &
                                  State_Grid = State_Grid(L), &
                                  lonCtr     = lonMidArr,     &
                                  latCtr     = latMidArr,     &
                                  lonEdge    = lonEdgeArr,    &
                                  latEdge    = latEdgeArr,    &
                                  RC         = RC         )
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
    CALL Accept_External_Date_Time( am_I_Root   = MasterProc,                &
                                    value_NYMDb = Input_Opt%NYMDb, &
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
    CALL Init_Error( MasterProc, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Error"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set a flag to denote if we should print ND70 debug output
    prtDebug            = ( Input_Opt%LPRT .and. MasterProc )

    ! Debug output
    IF ( prtDebug ) CALL Debug_Msg( '### MAIN: a READ_INPUT_FILE' )

    historyConfigFile = 'HISTORY.rc' ! InputOpt not yet initialized
    CALL Init_DiagList( MasterProc, historyConfigFile, Diag_List, RC )
    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_State_Met"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    !### Print diagnostic list if needed for debugging
    IF ( prtDebug ) CALL Print_DiagList( am_I_Root, Diag_List, RC )

    DO I = BEGCHUNK, ENDCHUNK
        am_I_Root = (MasterProc .AND. (I == BEGCHUNK))

        CALL GC_Init_StateObj( am_I_Root  = am_I_Root,     &  ! Root CPU (Y/N)?
     &                         Diag_List  = Diag_List,     &  ! Diagnostic list obj
     &                         Input_Opt  = Input_Opt,     &  ! Input Options
     &                         State_Chm  = State_Chm(I),  &  ! Chemistry State
     &                         State_Diag = State_Diag(I), &  ! Diagnostics State
     &                         State_Grid = State_Grid(I), &  ! Grid State
     &                         State_Met  = State_Met(I),  &  ! Meteorology State
     &                         RC         = RC            )   ! Success or failure

        ! Trap potential errors
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "GC_Init_StateObj"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        ! Start with v/v dry (CAM standard)
        State_Chm(I)%Spc_Units = 'v/v dry'

    ENDDO

    ! Now replicate GC_Init_Extra
    IF ( Input_Opt%LDryD) THEN

        ! Setup for dry deposition
        CALL Init_Drydep( am_I_Root = MasterProc,            &
     &                    Input_Opt = Input_Opt,             &
     &                    State_Chm = State_Chm(BEGCHUNK),   &
     &                    State_Diag = State_Diag(BEGCHUNK), &
     &                    State_Grid = State_Grid(BEGCHUNK), &
     &                    RC         = RC                   )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Drydep"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    !=================================================================
    ! Call setup routines for wet deposition
    !
    ! We need to initialize the wetdep module if either wet
    ! deposition or convection is turned on, so that we can do the
    ! large-scale and convective scavenging.  Also initialize the
    ! wetdep module if both wetdep and convection are turned off,
    ! but chemistry is turned on.  The INIT_WETSCAV routine will also
    ! allocate the H2O2s and SO2s arrays that are referenced in the
    ! convection code. (bmy, 9/23/15)
    !=================================================================
    IF ( Input_Opt%LConv .OR. &
         Input_Opt%LWetD .OR. &
         Input_Opt%LChem ) THEN
        CALL Init_WetScav( am_I_Root  = MasterProc,           &
     &                     Input_Opt  = Input_Opt,            &
     &                     State_Chm  = State_Chm(BEGCHUNK),  &
     &                     State_Diag = State_Diag(BEGCHUNK), &
     &                     State_Grid = State_Grid(BEGCHUNK), &
     &                     RC         = RC                   )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_WetScav"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Call SET_VDIFF_VALUES so that we can pass several values from
    ! Input_Opt to the vdiff_mod.F90.  This replaces the functionality
    ! of logical_mod.F and tracer_mod.F..  This has to be called
    ! after the input.geos file has been read from disk.
    !-----------------------------------------------------------------
    !CALL Set_VDiff_Values( am_I_Root = MasterProc,           &
    !&                       Input_Opt = Input_Opt,           &
    !&                       State_Chm = State_Chm(BEGCHUNK), &
    !&                       RC        = RC )

    !&IF (RC /= GC_SUCCESS) THEN
    !    ErrMsg = 'Error encountered in "Set_VDiff_Values"!'
    !    CALL Error_Stop( ErrMsg, ThisLoc )
    !ENDIF

    !-----------------------------------------------------------------
    ! Initialize the GET_NDEP_MOD for soil NOx deposition (bmy, 6/17/16)
    !-----------------------------------------------------------------
    !CALL Init_Get_NDep( am_I_Root  = MasterProc,           &
    !&                   Input_Opt  = Input_Opt,            &
    !&                   State_Chm  = State_Chm(BEGCHUNK),  &
    !&                   State_Diag = State_Diag(BEGCHUNK), &
    !&                   RC         = RC                   )
    !
    !IF (RC /= GC_SUCCESS) THEN
    !    ErrMsg = 'Error encountered in "Init_Get_NDep"!'
    !    CALL Error_Stop( ErrMsg, ThisLoc )
    !ENDIF

    !-----------------------------------------------------------------
    ! Initialize "carbon_mod.F"
    !-----------------------------------------------------------------
    IF (Input_Opt%LCarb) THEN
        CALL Init_Carbon( am_I_Root = MasterProc,            &
     &                    Input_Opt = Input_Opt,             &
     &                    State_Chm = State_Chm(BEGCHUNK),   &
     &                    State_Diag = State_Diag(BEGCHUNK), &
     &                    State_Grid = State_Grid(BEGCHUNK), &
     &                    RC         = RC                   )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Carbon"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF (Input_Opt%LDust) THEN
        CALL Init_Dust( am_I_Root  = MasterProc,           &
     &                  Input_Opt  = Input_Opt,            &
     &                  State_Chm  = State_Chm(BEGCHUNK),  &
     &                  State_Diag = State_Diag(BEGCHUNK), &
     &                  State_Grid = State_Grid(BEGCHUNK), &
     &                  RC         = RC                    )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Dust"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF (Input_Opt%LSSalt) THEN
        CALL Init_Seasalt( am_I_Root  = MasterProc,           &
     &                     Input_Opt  = Input_Opt,            &
     &                     State_Chm  = State_Chm(BEGCHUNK),  &
     &                     State_Diag = State_Diag(BEGCHUNK), &
     &                     State_Grid = State_Grid(BEGCHUNK), &
     &                     RC         = RC                    )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Seasalt"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF (Input_Opt%LSulf) THEN
        CALL Init_Sulfate( am_I_Root  = MasterProc,           &
     &                     Input_Opt  = Input_Opt,            &
     &                     State_Chm  = State_Chm(BEGCHUNK),  &
     &                     State_Diag = State_Diag(BEGCHUNK), &
     &                     State_Grid = State_Grid(BEGCHUNK), &
     &                     RC         = RC                    )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Sulfate"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    IF (Input_Opt%LSulf.OR.Input_Opt%LCarb.OR.Input_Opt%LDust.OR.Input_Opt%LSSalt) THEN
        CALL Init_Aerosol( am_I_Root  = MasterProc,           &
     &                     Input_Opt  = Input_Opt,            &
     &                     State_Chm  = State_Chm(BEGCHUNK),  &
     &                     State_Diag = State_Diag(BEGCHUNK), &
     &                     State_Grid = State_Grid(BEGCHUNK), &
     &                     RC         = RC                    )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Aerosol"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    CALL Init_Toms( am_I_Root  = MasterProc,           &
     &              Input_Opt  = Input_Opt,            &
     &              State_Chm  = State_Chm(BEGCHUNK),  &
     &              State_Diag = State_Diag(BEGCHUNK), &
     &              State_Grid = State_Grid(BEGCHUNK), &
     &              RC         = RC                    )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_TOMS"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! This is a bare subroutine - no module
    CALL NDXX_Setup( MasterProc,           &
     &               Input_Opt,            &
     &               State_Chm(BEGCHUNK),  &
     &               State_Grid(BEGCHUNK), &
     &               RC                    )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_NDXX_Setup"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF


    ! Set grid-cell area
    DO I = BEGCHUNK, ENDCHUNK
        ALLOCATE(Col_Area(NCOL(I)), STAT=IERR)
        IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Col_Area')

        CALL Get_Area_All_p(I, NCOL(I), Col_Area)
        
        ! Set default value (in case of chunks with fewer columns)
        State_Grid(I)%Area_M2 = 1.0e+10_fp
        DO iX = 1, NX
            DO iY = 1, NCOL(I)
                State_Grid(I)%Area_M2(iX,iY) = REAL(Col_Area(iY) * Re**2,fp)
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
    DO I = 1, (nZ+1)
        Ap_CAM_Flip(I) = hyai(nZ+2-I) * ps0 * 0.01e+0_r8
        Bp_CAM_Flip(I) = hybi(nZ+2-I)
    ENDDO

    DO I = BEGCHUNK, ENDCHUNK
    
        !-----------------------------------------------------------------
        ! Initialize the hybrid pressure module.  Define Ap and Bp.
        !-----------------------------------------------------------------
        CALL Init_Pressure( am_I_Root  = MasterProc,    &  ! Root CPU (Y/N)?
        &                   State_Grid = State_Grid(I), &  ! Grid State
        &                   RC         = RC            )   ! Success or failure

        ! Trapping errors
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Init_Pressure"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        !-----------------------------------------------------------------
        ! Pass external Ap and Bp to GEOS-Chem's Pressure_Mod
        !-----------------------------------------------------------------
        CALL Accept_External_ApBp( am_I_Root  = MasterProc,    &  ! Root CPU (Y/N)?
        &                          State_Grid = State_Grid(I), &  ! Grid State
        &                          ApIn       = Ap_CAM_Flip,   &  ! "A" term for hybrid grid
        &                          BpIn       = Bp_CAM_Flip,   &  ! "B" term for hybrid grid
        &                          RC         = RC            )   ! Success or failure

        ! Trapping errors
        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Error encountered in "Accept_External_ApBp"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

    ENDDO
    DEALLOCATE(Ap_CAM_Flip,Bp_CAM_Flip)



    ! Init_FJX..
    ! Init_PBL_Mix...
    ! Init_Chemistry...
    ! Init_TOMS...
    ! Emissions_Init...
    ! Init_UCX...
    ! Convert_Spc_Units...



    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default
    CALL addfld ( 'BCPI', (/'lev'/), 'A', 'mole/mole', trim('BCPI')//' mixing ratio' )
    CALL add_default ( 'BCPI',   1, ' ')

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT'

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer,   only: physics_buffer_desc

    TYPE(physics_state), INTENT(IN):: phys_state(begchunk:endchunk)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_INIT'

    ! This is when we want to update State_Met and so on
    ! Note that here we have been passed MANY chunks

  end subroutine chem_timestep_init

!===============================================================================

  subroutine GC_Update_Timesteps(DT)

    use Time_Mod,       only : Set_Timesteps
    use cam_abortutils, only : endrun

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

        CALL Set_Timesteps( MasterProc,            &
                            CHEMISTRY  =  DT_MIN,  &
                            EMISSION   =  DT_MIN,  &
                            DYNAMICS   =  DT_MIN,  &
                            UNIT_CONV  =  DT_MIN,  &
                            CONVECTION =  DT_MIN,  &
                            DIAGNOS    =  DT_MIN,  &
                            RADIATION  =  DT_MIN    )
        DT_MIN_LAST = DT_MIN
     ENDIF

  end subroutine

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o )

    use physics_buffer,   only: physics_buffer_desc
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t

    REAL(r8),            INTENT(IN)    :: dt          ! time step
    TYPE(physics_state), INTENT(IN)    :: state       ! Physics state variables
    TYPE(physics_ptend), INTENT(OUT)   :: ptend       ! indivdual parameterization tendencies
    TYPE(cam_in_t),      INTENT(INOUT) :: cam_in
    TYPE(cam_out_t),     INTENT(IN)    :: cam_out
    TYPE(physics_buffer_desc), POINTER :: pbuf(:)
    REAL(r8), OPTIONAL,  INTENT(OUT)   :: fh2o(pcols) ! h2o flux to balance source from chemistry

    ! Mapping (?)
    logical :: lq(pcnst)
    integer :: n, m

    integer :: lchnk, ncol
    ! Here's where you'll call DO_CHEMISTRY
    ! NOTE: State_Met etc are in an ARRAY - so we will want to always pass
    ! State_Met%(lchnk) and so on
    ! lchnk: which chunk we have on this process
    lchnk = state%lchnk
    ! ncol: number of atmospheric columns on this chunk
    ncol  = state%ncol

   ! Need to update the timesteps throughout the code
    CALL GC_Update_Timesteps(DT)

    ! Need to be super careful that the module arrays are updated and correctly
    ! set. NOTE: First thing - you'll need to flip all the data vertically

    ! 1. Update State_Met etc for this timestep

    ! 2. Copy tracers into State_Chm

    !if (MasterProc) WRITE(iulog,*) ' --> TEND SIZE: ', size(state%ncol)
    !if (MasterProc) WRITE(iulog,'(a,2(x,I6))') ' --> TEND SIDE:  ', lbound(state%ncol),ubound(state%ncol)

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_TEND'

    ! NOTE: Re-flip all the arrays vertically or suffer the consequences
    lq(:) = .false.
    DO n=1,pcnst
        !m = map2chm(n)
        m=0
        IF (m > 0) lq(n) = .true.
    ENDDO

    CALL physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)
    ! ptend%q dimensions: [column, ?, species]
    !ptend%q(:ncol,:,:) = 0.0e+0_r8
    !ptend%q(:ncol,:,:) = 0.0e+0_r8
    ptend%q(:,:,:) = 0.0e+0_r8
    IF (present(fh2o)) fh2o(:) = 0.0e+0_r8

    RETURN

  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    CHARACTER(LEN=*), INTENT(IN)  :: name       !  constituent name
    REAL(r8),         INTENT(IN)  :: latvals(:) ! lat in degrees (ncol)
    REAL(r8),         INTENT(IN)  :: lonvals(:) ! lon in degrees (ncol)
    LOGICAL,          INTENT(IN)  :: mask(:)    ! Only initialize where .true.
    REAL(r8),         INTENT(OUT) :: q(:,:)     ! kg tracer/kg dry air (ncol, pver
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    INTEGER :: ILEV, NLEV

    if (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT_CNST'

    NLEV = SIZE(Q, 2)
    IF ( ANY( TRACERNAMES .EQ. NAME ) ) THEN
       DO ILEV=1,NLEV
          WHERE(MASK)
             ! Set to the minimum mixing ratio
             q(:,ILEV) = 1.0e-38_r8
          ENDWHERE
       ENDDO
    ENDIF

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
    use WetScav_Mod,    only : Cleanup_Wetscav
    use Carbon_Mod,     only : Cleanup_Carbon
    use Dust_Mod,       only : Cleanup_Dust
    use Seasalt_Mod,    only : Cleanup_Seasalt
    use Aerosol_Mod,    only : Cleanup_Aerosol
    use TOMS_Mod,       only : Cleanup_Toms
    use Sulfate_Mod,    only : Cleanup_Sulfate
    use Pressure_Mod,   only : Cleanup_Pressure

    ! Special: cleans up after NDXX_Setup
    use Diag_Mod,       only : Cleanup_Diag

    INTEGER :: I, RC
    LOGICAL :: am_I_Root

    ! Finalize GEOS-Chem
    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_FINAL'
    
    CALL Cleanup_UCX( MasterProc )
    CALL Cleanup_Aerosol
    CALL Cleanup_Carbon
    CALL Cleanup_Drydep
    CALL Cleanup_Dust
    CALL Cleanup_FlexChem( am_I_Root, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_FlexChem"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_Pressure
    CALL Cleanup_Seasalt
    CALL Cleanup_Sulfate
    CALL Cleanup_Toms( MasterProc, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_Toms"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_WetScav( MasterProc, RC)
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_WetScav"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_Diag

    ! Cleanup Input_Opt
    CALL Cleanup_Input_Opt( MasterProc, Input_Opt, RC )

    ! Loop over each chunk and cleanup the variables
    DO I = BEGCHUNK, ENDCHUNK
        am_I_Root = ((I.eq.BEGCHUNK) .and. MasterProc)

        CALL Cleanup_State_Chm ( am_I_Root, State_Chm(I),  RC )
        CALL Cleanup_State_Diag( am_I_Root, State_Diag(I), RC )
        CALL Cleanup_State_Grid( am_I_Root, State_Grid(I), RC )
        CALL Cleanup_State_Met ( am_I_Root, State_Met(I),  RC )
    ENDDO
    CALL Cleanup_Error

    ! Finally deallocate state variables
    IF (ALLOCATED(State_Chm))     DEALLOCATE(State_Chm)
    IF (ALLOCATED(State_Diag))    DEALLOCATE(State_Diag)
    IF (ALLOCATED(State_Grid))    DEALLOCATE(State_Grid)
    IF (ALLOCATED(State_Met))     DEALLOCATE(State_Met)

    IF (MasterProc) WRITE(iulog,'(a,4(x,L1))') ' --> DEALLOC CHECK : ', ALLOCATED(State_Chm), ALLOCATED(State_Diag), ALLOCATED(State_Grid), ALLOCATED(State_Met)

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
    use camsrfexch, only : cam_in_t

    ! Arguments:

    TYPE(physics_state),    INTENT(IN)    :: state   ! Physics state variables
    TYPE(cam_in_t),         INTENT(INOUT) :: cam_in  ! import state

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_EMISSIONS'

  end subroutine chem_emissions

end module chemistry
