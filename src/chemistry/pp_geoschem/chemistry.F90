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

  use Chem_Mods,           only : NSlvd, Slvd_Lst, Slvd_Ref_MMR

  ! Exit routine in CAM
  use cam_abortutils,      only : endrun

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
  integer, parameter :: NTracersMax = 200    ! Must be equal to nadv_chem
  integer            :: NTracers
  character(len=255) :: TracerNames(NTracersMax)
  character(len=255) :: TracerLongNames(NTracersMax)
  integer            :: Indices(NTracersMax)
  real(r8)           :: Adv_Mass(NTracersMax)
  real(r8)           :: MWRatio(NTracersMax)
  real(r8)           :: Ref_MMR(NTracersMax)

  ! Short-lived species (i.e. not advected)
  integer, parameter :: NSlsMax = 500        ! UNadvected species only
  integer            :: NSls
  character(len=255) :: SlsNames(NSlsMax)
  character(len=255) :: SlsLongNames(NSlsMax)
  real(r8)           :: Sls_Ref_MMR(NSlsMax)
  real(r8)           :: SLSMWRatio(NSlsMax)
  !===== SDE DEBUG =====

  ! Location of valid input.geos
  CHARACTER(LEN=500) :: inputGeosPath

  ! Location of chemistry input (for now)
  CHARACTER(LEN=500) :: chemInputsDir

  ! Mapping between constituents and GEOS-Chem tracers
  INTEGER :: Map2GC(pcnst)
  INTEGER :: Map2GC_Sls(NSlsMax)

  ! Mapping from constituents to raw index
  INTEGER :: Map2Idx(pcnst)

  !-----------------------------
  ! Derived type objects
  !-----------------------------
  TYPE(OptInput)                  :: Input_Opt      ! Input Options object
  TYPE(ChmState),ALLOCATABLE      :: State_Chm(:)   ! Chemistry State object
  TYPE(DgnState),ALLOCATABLE      :: State_Diag(:)  ! Diagnostics State object
  TYPE(GrdState),ALLOCATABLE      :: State_Grid(:)  ! Grid State object
  TYPE(MetState),ALLOCATABLE      :: State_Met(:)   ! Meteorology State object
  TYPE(DgnList )                  :: Diag_List      ! Diagnostics list object

  ! Indices of critical species
  INTEGER                         :: iH2O, iO3, iCH4, iCO

  ! Indices in the physics buffer
  INTEGER                         :: NDX_PBLH    ! PBL height [m]
  INTEGER                         :: NDX_FSDS    ! Downward shortwave flux at surface [W/m2]
  INTEGER                         :: NDX_CLDTOP  ! Cloud top height [?]
  INTEGER                         :: NDX_CLDFRC  ! Cloud fraction [?]
  INTEGER                         :: NDX_PRAIN   ! Rain production rate [?]

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
    use PhysConst,      only : MWDry

    use Short_Lived_Species, only : Register_Short_Lived_Species

    use State_Grid_Mod, only : Init_State_Grid, Cleanup_State_Grid
    use State_Chm_Mod,  only : Init_State_Chm, Cleanup_State_Chm
    use State_Chm_Mod,  only : Ind_
    use Input_Opt_Mod,  only : Set_Input_Opt,  Cleanup_Input_Opt
    use Species_Mod,    only : Species

    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents for chemistry
    !
    !-----------------------------------------------------------------------
    ! Need to generate a temporary species database - therefore temp State_Chm
    Type(ChmState)         :: SC
    Type(GrdState)         :: SG
    Type(OptInput)         :: IO
    TYPE(Species), POINTER :: ThisSpc

    INTEGER            :: I, N, M
    REAL(r8)           :: cptmp
    REAL(r8)           :: mwtmp
    REAL(r8)           :: qmin
    REAL(r8)           :: Ref_VMR
    CHARACTER(LEN=128) :: mixtype
    CHARACTER(LEN=128) :: molectype
    CHARACTER(LEN=128) :: Lng_Name
    LOGICAL            :: camout
    LOGICAL            :: ic_from_cam2
    LOGICAL            :: has_fixed_ubc
    LOGICAL            :: has_fixed_ubflx

    INTEGER            :: RC
    ! SDE 2018-05-02: This seems to get called before anything else
    ! That includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_REGISTER'

    ! Generate fake state_chm
    IO%Max_BPCH_Diag       = 1000
    IO%Max_AdvectSpc       = 500
    IO%Max_Families        = 250

    IO%RootCPU             = .False.

    CALL Set_Input_Opt( am_I_Root = .False., &
                        Input_Opt = IO,      &
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
        IO%AdvectSpc_Name(I) = TRIM(tracernames(I))
    ENDDO
    IO%SalA_rEdge_um(1)    = 0.01e+0_fp
    IO%SalA_rEdge_um(2)    = 0.50e+0_fp
    IO%SalC_rEdge_um(1)    = 0.50e+0_fp
    IO%SalC_rEdge_um(2)    = 8.00e+0_fp

    ! Prevent reporting
    IO%rootCPU             = .False.
    Input_Opt%myCPU        = myCPU

    CALL Init_State_Grid( am_I_Root  = .False.,       &
                          State_Grid = State_Grid(I), &
                          RC         = RC         )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Init_State_Grid"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    State_Grid(I)%NX = 1
    State_Grid(I)%NY = 1
    State_Grid(I)%NZ = 1

    CALL Init_State_Chm( am_I_Root  = .False., &
                         Input_Opt  = IO,      &
                         State_Chm  = SC,      &
                         State_Grid = SG,      &
                         RC         = RC      )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Init_State_Chm"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! At the moment, we force nadv_chem=200 in the setup file
    ! Default
    Map2GC = -1
    Ref_MMR(:) = 0.0e+0_r8
    MWRatio(:) = 1.0e+0_r8
    TracerLongNames = ''

    DO I = 1, NTRACERSMAX
        IF (I.LE.NTRACERS) THEN
            N           = Ind_(TracerNames(I))
            ThisSpc => SC%SpcData(N)%Info
            Lng_Name    = TRIM(ThisSpc%FullName)
            MWTmp       = REAL(ThisSpc%MW_g,r8)
            Ref_VMR     = REAL(ThisSpc%BackgroundVV,r8)
            Adv_Mass(I) = MWTmp
            Ref_MMR(I)  = Ref_VMR / (MWDry / MWTmp)
       ELSE
           Lng_Name    = TRIM(TracerNames(I))
           MWTmp       = 1000.0e+0_r8 * (0.001e+0_r8)
           Adv_Mass(I) = MWTmp
           Ref_MMR(I)  = 1.0e-38_r8
       ENDIF
       MWRatio(I) = MWDry/MWTmp
       TracerLongNames(I) = TRIM(Lng_Name)

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
       !write(tracernames(i),'(a,I0.4)') 'GCTRC_', i
       ! NOTE: In MOZART, this only gets called for tracers
       ! This is the call to add a "constituent"
       CALL cnst_add( trim(tracernames(i)), adv_mass(i), cptmp, qmin, n, &
                      readiv=ic_from_cam2, mixtype=mixtype, cam_outfld=camout, &
                      molectype=molectype, fixed_ubc=has_fixed_ubc, &
                      fixed_ubflx=has_fixed_ubflx, longname=trim(lng_name) )

       ! Add to GC mapping. When starting a timestep, we will want to update the
       ! concentration of State_Chm(x)%Species(1,iCol,iLev,m) with data from
       ! constituent n
       M = Ind_(TRIM(TracerNames(I)))
       IF ( M > 0 ) THEN
           Map2GC(N)  = M
           Map2Idx(N) = I
       ENDIF
       ! Nullify pointer
       ThisSpc => NULL()
    ENDDO

    ! Now unadvected species
    Map2GC_Sls = 0
    Sls_Ref_MMR(:) = 0.0e+0_r8
    SlsMWRatio(:)  = -1.0e+0_r8
    SlsLongNames = ''
    DO I = 1, NSls
        N = Ind_(SlsNames(I))
        IF (N.GT.0) THEN
            ThisSpc => SC%SpcData(N)%Info
            MWTmp           = REAL(ThisSpc%MW_g,r8)
            Ref_VMR         = REAL(ThisSpc%BackgroundVV,r8)
            Lng_Name        = TRIM(ThisSpc%FullName)
            SlsLongNames(I) = Lng_Name
            Sls_Ref_MMR(I)  = Ref_VMR / (MWDry / MWTmp)
            SlsMWRatio(I)   = MWDry / MWTmp
            Map2GC_Sls(I)   = N
            ThisSpc  => NULL()
        ENDIF
    ENDDO

    ! Pass information to "short_lived_species" module
    Slvd_Ref_MMR(1:NSls) = Sls_Ref_MMR(1:NSls)
    CALL Register_Short_Lived_Species()
       ! More information:
       ! http://www.cesm.ucar.edu/models/atm-cam/docs/phys-interface/node5.html

    ! Clean up
    Call Cleanup_State_Chm ( .False., SC, RC )
    Call Cleanup_State_Grid( .False., SG, RC )
    Call Cleanup_Input_Opt ( .False., IO, RC )

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
        IF ( NSPEC > NSlsMax ) THEN
            CALL ENDRUN('chem_readnl: too many species - increase NSlsmax')
        ENDIF

        NSls = 0
        DO I=1,NSPEC
            ! Get the name of the species from KPP
            LINE = ADJUSTL(TRIM(SPC_NAMES(I)))
            ! Only add this
            validSLS = ( (.NOT.ANY(TRIM(LINE).EQ.TRACERNAMES)).AND.&
                         (.NOT.(LINE(1:2) == 'RR')) )
            IF (validSLS) THEN
                ! Genuine new short-lived species
                NSls = NSls + 1
                SLSNAMES(NSls) = TRIM(LINE)
                WRITE(iulog,'(a,I5,a,a)') ' --> GC species ', NSls, ': ', TRIM(LINE)
            ENDIF
        ENDDO
    ENDIF

    ! Broadcast to all processors
#if defined( SPMD )
    CALL MPIBCAST(NTracers,    1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST(TracerNames, LEN(TracerNames(1))*NTracersMax, MPICHAR, 0, MPICOM )
    CALL MPIBCAST(NSls,        1,                               MPIINT,  0, MPICOM )
    CALL MPIBCAST(SlsNames,    LEN(SlsNames(1))*NSlsMax,        MPICHAR, 0, MPICOM )
#endif

    ! Update "short_lived_species" arrays - will eventually unify these
    NSlvd = NSls
    ALLOCATE(Slvd_Lst(NSlvd), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Slvd_Lst')
    ALLOCATE(Slvd_Ref_MMR(NSlvd), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Slvd_Ref_MMR')
    DO I=1,NSls
        Slvd_Lst(I) = TRIM(SlsNames(I))
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
    use physics_buffer, only: physics_buffer_desc, pbuf_get_index
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
    use Chemistry_Mod, only : Init_Chemistry
    use UCX_Mod,       only : Init_UCX

    use PBL_Mix_Mod,   only : Init_PBL_Mix

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
    LOGICAL               :: am_I_Root, rootChunk
    LOGICAL               :: prtDebug

    ! Strings
    CHARACTER(LEN=255)       :: historyConfigFile
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
    ! NLEV: number of vertical levels
    NLEV  = PVER

    write(iulog,'(2(a,x,I6,x))') 'chem_init called on PE ', myCPU, ' of ', nCPUs

    ! The GEOS-Chem grids on every "chunk" will all be the same size, to avoid
    ! the possibility of having differently-sized chunks
    NX = 1
    !NY = MAXVAL(NCOL)
    NY = PCOLS
    NZ = NLEV

    !! Add short lived speies to buffers
    !CALL Pbuf_add_field(Trim(SLSBuffer),'global',dtype_r8,(/PCOLS,PVER,NSls/),Sls_Pbf_Idx)
    !! Initialize
    !ALLOCATE(SlsPtr(PCOLS,PVER,BEGCHUNK:ENDCHUNK), STAT=IERR)
    !IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating SlsPtr')
    !SlsPtr(:,:,:) = 0.0e+0_r8
    !DO I=1,NSls
    !   SlsPtr(:,:,:) = Sls_Ref_MMR(I)
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

        ! Set some basic flags
    Input_Opt%Max_BPCH_Diag     = 1000
    Input_Opt%Max_AdvectSpc     = 500
    Input_Opt%Max_Families      = 250

    ! Initialize fields of the Input Options object
    CALL Set_Input_Opt( am_I_Root = MasterProc, &
                        Input_Opt = Input_Opt,  &
                        RC        = RC         )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DO I = BEGCHUNK, ENDCHUNK

        ! Only treat the first chunk as the "root"
        am_I_Root = ((I.EQ.BEGCHUNK) .and. MasterProc)

        ! Initialize fields of the Grid State object
        CALL Init_State_Grid( am_I_Root  = am_I_Root,      &
                              State_Grid = State_Grid(I),  &
                              RC         = RC         )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered within call to "Init_State_Grid"!'
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

    Input_Opt%myCPU    = myCPU
    Input_Opt%rootCPU  = MasterProc


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

    IF (Input_Opt%LChem) THEN
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

    CALL Init_PBL_Mix( am_I_Root  = MasterProc,           &
                       State_Grid = State_Grid(BEGCHUNK), &
                       RC         = RC                   )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered in "Init_PBL_Mix"!'
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

    IF (Input_Opt%Its_A_FullChem_Sim .OR. &
        Input_Opt%Its_An_Aerosol_Sim) THEN
        ! This also initializes Fast-JX
        CALL Init_Chemistry( am_I_Root  = MasterProc,           &
     &                       Input_Opt  = Input_Opt,            &
     &                       State_Chm  = State_Chm(BEGCHUNK),  &
     &                       State_Diag = State_Diag(BEGCHUNK), &
     &                       State_Grid = State_Grid(BEGCHUNK), &
     &                       RC         = RC                    )

        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Chemistry"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    ! Initialize HEMCO?
    !CALL EMISSIONS_INIT ( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
    !                      HcoConfig=HcoConfig )
    !ASSERT_(RC==GC_SUCCESS)

    IF (Input_Opt%LChem.and.Input_Opt%LUCX) THEN
        CALL Init_UCX( am_I_Root  = MasterProc,           &
     &                 Input_Opt  = Input_Opt,            &
     &                 State_Chm  = State_Chm(BEGCHUNK),  &
     &                 State_Diag = State_Diag(BEGCHUNK), &
     &                 State_Grid = State_Grid(BEGCHUNK) )
    ENDIF

    ! Get the index of H2O
    iH2O = Ind_('H2O')
    iO3  = Ind_('O3')
    iCH4 = Ind_('CH4')
    iCO  = Ind_('CO')

    ! Get indices for physical fields in physics buffer
    NDX_PBLH    = Pbuf_Get_Index('PblH'  )
    NDX_FSDS    = Pbuf_Get_Index('Fsds'  )
    NDX_CLDTOP  = Pbuf_Get_Index('CldTop')
    NDX_CLDFRC  = Pbuf_Get_Index('Cld'   )
    NDX_PRAIN   = Pbuf_Get_Index('PRain' )

    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default
    ! Add all species as output fields if desired
    DO I = 1, NTracers
        SpcName = TRIM(TracerNames(I))
        CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', TRIM(TracerLongNames(I))//' concentration')
        IF (TRIM(SpcName) == 'O3') THEN
            CALL Add_Default ( TRIM(SpcName), 1, ' ')
        ENDIF
    ENDDO
    DO I =1, NSls
        SpcName = TRIM(SlsNames(I))
        CALL AddFld( TRIM(SpcName), (/ 'lev' /), 'A', 'mol/mol', TRIM(SlsLongNames(I))//' concentration')
        !CALL Add_Default(TRIM(SpcName), 1, '')
    ENDDO
    !CALL AddFld ( 'BCPI', (/'lev'/), 'A', 'mole/mole', trim('BCPI')//' mixing ratio' )
    !CALL Add_Default ( 'BCPI',   1, ' ')

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

  subroutine chem_timestep_tend( State, ptend, cam_in, cam_out, dT, pbuf,  fh2o )

    use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t

    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p

    use Dao_Mod,          only: Set_Dry_Surface_Pressure
    use Dao_Mod,          only: AirQnt
    use GC_Grid_Mod,      only: SetGridFromCtr
    use Pressure_Mod,     only: Set_Floating_Pressures
    use Pressure_Mod,     only: Accept_External_Pedge
    use Time_Mod,         only: Accept_External_Date_Time
    use Strat_chem_Mod,   only: Init_Strat_Chem
    use Toms_Mod,         only: Compute_Overhead_O3
    use Chemistry_Mod,    only: Do_Chemistry
    use Wetscav_Mod,      only: Setup_Wetscav
    use CMN_Size_Mod,     only: PTop

    use Tropopause,          only: Tropopause_findChemTrop, Tropopause_Find

    ! For calculating SZA
    use Orbit,            only: zenith
    use Time_Manager,        only: Get_Curr_Calday, Get_Curr_Date

    ! Calculating relative humidity
    use WV_Saturation,    only: QSat
    use PhysConst,        only: MWDry

    ! Grid area
    use PhysConst,        only: Gravit
    use PhysConstants,    only: Re
    use Phys_Grid,        only: get_area_all_p

    use Short_Lived_Species, only : Get_Short_Lived_Species
    use Short_Lived_Species, only : Set_Short_Lived_Species

    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,    only: PI, PI_180

    REAL(r8),            INTENT(IN)    :: dT          ! Time step
    TYPE(physics_state), INTENT(IN)    :: State       ! Physics State variables
    TYPE(physics_ptend), INTENT(OUT)   :: ptend       ! indivdual parameterization tendencies
    TYPE(cam_in_t),      INTENT(INOUT) :: cam_in
    TYPE(cam_out_t),     INTENT(IN)    :: cam_out
    TYPE(physics_buffer_desc), POINTER :: pbuf(:)
    REAL(r8), OPTIONAL,  INTENT(OUT)   :: fh2o(PCOLS) ! h2o flux to balance source from chemistry

    ! Initial MMR for all species
    REAL(r8) :: MMR_Beg(PCOLS,PVER,NSls+NTracers)
    REAL(r8) :: MMR_End(PCOLS,PVER,NSls+NTracers)
    REAL(r8) :: MMR_TEnd(PCOLS,PVER,NSls+NTracers)


    ! Mapping (?)
    LOGICAL :: lq(pcnst)

    ! Indexing
    INTEGER :: I, J, K, L, N, M
    INTEGER :: NX, NY, NZ

    INTEGER :: LCHNK, NCOL

    REAL(r8), DIMENSION(State%NCOL) :: &
        CSZA,                          &           ! Cosine of solar zenith angle
        Zsurf,     &                               ! Surface height
        Rlats, Rlons                               ! Chunk latitudes and longitudes (radians)

    REAL(r8), POINTER :: PblH(:)                      ! PBL height on each chunk [m]
    REAL(r8), POINTER :: CldTop(:)                    ! Cloud top height 
    REAL(r8), POINTER :: CldFrc(:)                    ! Cloud fraction
    REAL(r8), POINTER :: Fsds(:)                      ! Downward shortwave flux at surface [W/m2]
    REAL(r8), POINTER :: PRain(:)                     ! Rain production rate
    REAL(r8)          :: RelHum(State%NCOL, PVER)     ! Relative humidity [0-1]
    REAL(r8)          :: SatV  (State%NCOL, PVER)     ! Work arrays
    REAL(r8)          :: SatQ  (State%NCOL, PVER)     ! Work arrays
    REAL(r8)          :: QH2O  (State%NCOL, PVER)     ! Specific humidity [kg/kg]
    REAL(r8)          :: H2OVMR(State%NCOL, PVER)     ! H2O volume mixing ratio

    ! Because of strat chem
    LOGICAL, SAVE :: SCHEM_READY = .FALSE.

    REAL(f4)      :: lonMidArr(1,PCOLS), latMidArr(1,PCOLS)
    INTEGER       :: iMaxLoc(1)

    REAL(r8)      :: Col_Area(State%NCOL)

    ! Intermediate arrays
    INTEGER      :: Trop_Lev(PCOLS)
    REAL(r8)     :: Trop_P(  PCOLS)
    REAL(r8)     :: Trop_T(  PCOLS)
    REAL(r8)     :: Trop_Ht( PCOLS)
    REAL(r8)     :: SnowDepth(PCOLS)
    REAL(r8)     :: Z0(PCOLS)
    REAL(r8)     :: Sd_Ice, Sd_Lnd, Sd_Avg, Frc_Ice

    ! Calculating SZA
    REAL(r8)      :: Calday

    ! For archiving
    CHARACTER(LEN=255) :: SpcName
    REAL(r8)           :: VMR(State%NCOL,PVER)
    REAL(r8)           :: MMR0, MMR1, Mass0, Mass1, AirMass, Mass10r, Mass10a
    REAL(r8)           :: MMR_Min, MMR_Max

    REAL(r8)           :: SlsData(State%NCOL, PVER, NSls)

    INTEGER      :: CurrYr, CurrMo, CurrDy, CurrTOD
    INTEGER      :: CurrYMD, CurrHMS, CurrHr, CurrMn, CurrSc
    REAL(f4)     :: CurrUTC

    INTEGER, SAVE :: iStep = 0
    LOGICAL       :: rootChunk
    INTEGER       :: RC

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

    NX = 1
    NY = NCOL

    ! Update the grid lat/lons since they are module variables
    ! Assume (!) that area hasn't changed for now, as GEOS-Chem will
    ! retrieve this from State_Met which is chunked
    !CALL get_rlat_all_p( LCHNK, NCOL, Rlats )
    !CALL get_rlon_all_p( LCHNK, NCOL, Rlons )
    Rlats(1:NCOL) = State%Lat(1:NCOL)
    Rlons(1:NCOL) = State%Lon(1:NCOL)

    lonMidArr = 0.0e+0_f4
    latMidArr = 0.0e+0_f4
    DO I = 1, NX
        DO J = 1, NY
            lonMidArr(I,J) = REAL(Rlons(J), f4)
            latMidArr(I,J) = REAL(Rlats(J), f4)
        ENDDO
    ENDDO

    ! Update the grid
    Call SetGridFromCtr( am_I_Root  = rootChunk,         &
                         State_Grid = State_Grid(LCHNK), &
                         lonCtr     = lonMidArr,         &
                         latCtr     = latMidArr,         &
                         RC         = RC )

    IF ( RC /= GC_SUCCESS ) THEN
        ErrMsg = 'Error encountered within call to "SetGridFromCtr"!'
        CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set area
    CALL Get_Area_All_p( LCHNK, NCOL, Col_Area )
    ! Set default value (in case of chunks with fewer columns)
    State_Grid(LCHNK)%Area_M2 = 1.0e+10_fp
    DO J = 1, NCOL
        State_Grid(LCHNK)%Area_M2(1,J) = REAL(Col_Area(J) * Re**2,fp)
    ENDDO
    State_Met(LCHNK)%Area_M2 = State_Grid(LCHNK)%Area_M2


    ! 2. Copy tracers into State_Chm
    ! Data was received in kg/kg dry
    State_Chm(LCHNK)%Spc_Units = 'kg/kg dry'
    ! Initialize ALL State_Chm species data to zero, not just tracers
    State_Chm(LCHNK)%Species = 0.0e+0_fp

    lq(:) = .FALSE.

    MMR_Beg = 0.0e+0_r8
    DO N = 1, pcnst
        M = Map2GC(N)
        IF (M > 0) THEN
            I = 1
            DO J = 1, NCOL
                DO K = 1, PVER
                    ! CURRENTLY KG/KG DRY
                    MMR_Beg(J,K,M) = State%q(J,PVER+1-K,N)
                    State_Chm(LCHNK)%Species(1,J,K,M) = REAL(MMR_Beg(J,K,M),fp)
                ENDDO
            ENDDO
            lq(N) = .TRUE.
        ENDIF
    ENDDO

    ! Retrieve previous value of species data
    SlsData(:,:,:) = 0.0e+0_r8
    CALL Get_Short_Lived_Species( SlsData, LCHNK, NCOL, Pbuf )

    ! Remap and flip them
    DO N = 1, NSls
        M = Map2GC_Sls(N)
        IF (M > 0) THEN
            DO J = 1, NCOL
                DO K = 1, PVER
                    State_Chm(LCHNK)%Species(1,J,K,M) = REAL(SlsData(J,PVER+1-K,N),fp)
                ENDDO
            ENDDO
        ENDIF
    ENDDO

    ! Initialize tendency array
    CALL Physics_ptend_init(ptend, State%psetcols, 'chemistry', lq=lq)

    ! Calculate COS(SZA)
    Calday = Get_Curr_Calday( )
    CALL Zenith( Calday, Rlats, Rlons, CSZA, NCOL )
    !CALL Outfld( 'SZA', SZA, NCOL, LCHNK )

    CALL Pbuf_Get_Field( Pbuf, NDX_PBLH,     PblH   )
    CALL Pbuf_Get_Field( Pbuf, NDX_PRAIN,    PRain  )
    CALL Pbuf_Get_Field( Pbuf, NDX_FSDS,     Fsds   )
    CALL Pbuf_Get_Field( Pbuf, NDX_CLDTOP,   CldTop )
    CALL Pbuf_Get_Field( Pbuf, NDX_CLDFRC,   CldFrc )

    ! Get VMR and MMR of H2O
    H2OVMR = 0.0e0_fp
    QH2O   = 0.0e0_fp
    ! Note MWDRY = 28.966 g/mol

    DO J = 1, NY
        DO L = 1, NZ
            QH2O(J,L) = REAL(State_Chm(LCHNK)%Species(1,J,K,iH2O),r8)
            H2OVMR(J,L) = QH2O(J,L) * MWDry / 18.016e+0_fp
        ENDDO
    ENDDO

    ! Calculate RH (range 0-1, note still level 1 = TOA)
    RELHUM(:,:) = 0.0e+0_r8
    CALL QSat(State%T(:NCOL,:), State%Pmid(:NCOL,:), SatV, SatQ)
    DO J = 1, NY
        DO L = 1, NZ
            RELHUM(J,L) = 0.622e+0_r8 * H2OVMR(J,L) / SatQ(J,L)
            RELHUM(J,L) = MAX( 0.0e+0_r8, MIN( 1.0e+0_r8, RELHUM(J,L) ) )
        ENDDO
    ENDDO

    ! Estimate roughness height VERY roughly
    Z0 = 0.0e+0_r8
    DO I = 1, NCOL 
        IF (Cam_in%LandFrac(I).GE.0.5e+0_r8) THEN
            Z0(I) = 0.035e+0_r8
        ELSE
            Z0(I) = 0.0001e+0_r8
        ENDIF
    ENDDO

    ! Retrieve tropopause level
    Trop_Lev = 0.0e+0_r8
    CALL Tropopause_FindChemTrop(State, Trop_Lev)
    ! Back out the pressure
    Trop_P = 1000.0e+0_r8
    DO I = 1, NCOL
        Trop_P(I) = State%PMid(I,Trop_Lev(I)) * 0.01e+0_r8
    ENDDO

   ! Calculate snow depth
   SnowDepth = 0.0e+0_r8
   DO I = 1, NCOL
       Sd_Ice  = MAX(0.0e+0_r8,Cam_in%Snowhice(I))
       Sd_Lnd  = MAX(0.0e+0_r8,Cam_in%Snowhland(I))
       Frc_Ice = MAX(0.0e+0_r8,Cam_in%Icefrac(I))
       IF (Frc_Ice > 0.0e+0_r8) THEN
           Sd_Avg = (Sd_Lnd*(1.0e+0_r8 - Frc_Ice)) + (Sd_Ice * Frc_Ice)
       ELSE
           Sd_Avg = Sd_Lnd
       ENDIF
       SnowDepth(I) = Sd_Avg
   ENDDO

    State_Met(LCHNK)%ALBD      (1,:) = Cam_in%Asdir(:)
    State_Met(LCHNK)%CLDFRC    (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%EFLUX     (1,:) = Cam_in%Lhf(:)
    State_Met(LCHNK)%HFLUX     (1,:) = Cam_in%Shf(:)
    State_Met(LCHNK)%FRCLND    (1,:) = 0.0e+0_fp ! Olson land fraction
    State_Met(LCHNK)%FRLAND    (1,:) = Cam_in%LandFrac(:)
    State_Met(LCHNK)%FROCEAN   (1,:) = Cam_in%OcnFrac(:) + Cam_in%IceFrac(:)
    State_Met(LCHNK)%FRSEAICE  (1,:) = Cam_in%IceFrac(:)
    State_Met(LCHNK)%FRLAKE    (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%FRLANDIC  (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%GWETROOT  (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%GWETTOP   (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%LAI       (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%PARDR     (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%PARDF     (1,:) = 0.0e+0_fp
    State_Met(LCHNK)%PBLH      (1,:) = PblH(:NCOL)
    State_Met(LCHNK)%PHIS      (1,:) = State%Phis(:)
    State_Met(LCHNK)%PRECANV   (1,:) = 0.0e+0_fp                           ! Not used
    State_Met(LCHNK)%PRECCON   (1,:) = Cam_Out%Precc(:)                    ! Convective precip
    State_Met(LCHNK)%PRECLSC   (1,:) = Cam_Out%Precl(:)                    ! "Stratiform" precip
    State_Met(LCHNK)%PRECTOT   (1,:) = Cam_Out%Precc(:) + Cam_Out%Precl(:) ! All precip
    State_Met(LCHNK)%TROPP     (1,:) = Trop_P(:)
    State_Met(LCHNK)%PS1_WET   (1,:) = State%ps(:)*0.01e+0_fp
    State_Met(LCHNK)%PS2_WET   (1,:) = State%ps(:)*0.01e+0_fp
    State_Met(LCHNK)%SLP       (1,:) = State%ps(:)*0.01e+0_fp
    State_Met(LCHNK)%TS        (1,:) = Cam_in%TS(:)
    State_Met(LCHNK)%TSKIN     (1,:) = Cam_in%TS(:)
    State_Met(LCHNK)%SWGDN     (1,:) = fsds(:)
    State_Met(LCHNK)%TO3       (1,:) = 300.0e+0_fp ! Dummy
    State_Met(LCHNK)%SNODP     (1,:) = snowdepth(:)
    State_Met(LCHNK)%SNOMAS    (1,:) = snowdepth(:) * 1000.0e+0_r8 ! m -> kg/m2 for ice w/rho ~ 1000 kg/m3
    State_Met(LCHNK)%SUNCOS    (1,:) = CSZA(:)
    State_Met(LCHNK)%SUNCOSmid (1,:) = CSZA(:)
    State_Met(LCHNK)%U10M      (1,:) = State%U(:,NZ)
    State_Met(LCHNK)%USTAR     (1,:) = Cam_In%UStar(:)
    State_Met(LCHNK)%V10M      (1,:) = State%V(:,NZ)
    State_Met(LCHNK)%Z0        (1,:) = 0.0e+0_fp

    DO J = 1, NY
        DO I = 1, NX
            iMaxLoc = MAXLOC( (/ State_Met(LCHNK)%FRLAND(I,J)   + &
                                 State_Met(LCHNK)%FRLANDIC(I,J) + &
                                 State_Met(LCHNK)%FRLAKE(I,J),    &
                                 State_Met(LCHNK)%FRSEAICE(I,J),  &
                                 State_Met(LCHNK)%FROCEAN(I,J)  - &
                                 State_Met(LCHNK)%FRSEAICE(I,J) /) )
            IF ( iMaxLoc(1) == 3 ) iMaxLoc(1) = 0
            ! reset ocean to 0
            State_Met(LCHNK)%LWI(I,J) = FLOAT( iMaxLoc(1) )
        ENDDO
    ENDDO

    ! Three-dimensional fields on level edges
    DO J = 1, NY
        DO L = 1, NZ+1 
            State_Met(LCHNK)%CMFMC   (1,J,L) = 0.0e+0_fp
            State_Met(LCHNK)%PFICU   (1,J,L) = 0.0e+0_fp
            State_Met(LCHNK)%PFILSAN (1,J,L) = 0.0e+0_fp
            State_Met(LCHNK)%PFLCU   (1,J,L) = 0.0e+0_fp
            State_Met(LCHNK)%PFLLSAN (1,J,L) = 0.0e+0_fp
            State_Met(LCHNK)%PEDGE   (1,J,L) = State%Pint(J,NZ+2-L)*0.01e+0_fp
        ENDDO
    ENDDO

    ! These are set later
    State_Met(LCHNK)%PS1_DRY (:,:) = 0.0e+0_fp
    State_Met(LCHNK)%PS2_DRY (:,:) = 0.0e+0_fp


    ! Calculate CLDTOPS (highest location of CMFMC in the column)
    I = 1
    DO J = 1, NY
        State_Met(lchnk)%CldTops(I,J) = 1
        DO L = NZ, 1, -1
            IF ( State_Met(LCHNK)%CMFMC(I,J,L) > 0.0e+0_fp ) THEN
                State_Met(LCHNK)%CLDTOPS(I,J) = L + 1
                EXIT
            ENDIF
        ENDDO
    ENDDO

    ! Three-dimensional fields on level centers
    DO J = 1, NY
        DO L = 1, NZ 
            State_Met(lchnk)%U        (1,J,L) = State%U(J,NZ+1-L)
            State_Met(lchnk)%V        (1,J,L) = State%V(J,NZ+1-L)
            !State_Met(lchnk)%OMEGA    (1,J,L) = State%Omega(J,NZ+1-L)
            State_Met(lchnk)%CLDF     (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%DTRAIN   (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%DQRCU    (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%DQRLSAN  (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%QI       (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%QL       (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%RH       (1,J,L) = RelHum(J,NZ+1-L)  * 100.0e+0_fp
            State_Met(lchnk)%TAUCLI   (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%TAUCLW   (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%REEVAPCN (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%REEVAPLS (1,J,L) = 0.0e+0_fp
            State_Met(lchnk)%SPHU1    (1,J,L) = QH2O(J,NZ+1-L)    * 1.0e+3_fp    ! g/kg
            State_Met(lchnk)%SPHU2    (1,J,L) = QH2O(J,NZ+1-L)    * 1.0e+3_fp    ! g/kg
            State_Met(lchnk)%TMPU1    (1,J,L) = State%T(J,NZ+1-L)
            State_Met(lchnk)%TMPU2    (1,J,L) = State%T(J,NZ+1-L)
        ENDDO
    ENDDO

    ! Derived fields
    State_Met(lchnk)%T    = (State_Met(lchnk)%TMPU1 + State_Met(lchnk)%TMPU2)*0.5e+0_fp
    State_Met(lchnk)%SPHU = (State_Met(lchnk)%SPHU1 + State_Met(lchnk)%SPHU2)*0.5e+0_fp

    ! Calculate total OD as liquid cloud OD + ice cloud OD
    State_Met(lchnk)%OPTD =  State_Met(lchnk)%TAUCLI + State_Met(lchnk)%TAUCLW

    ! Nullify all pointers
    NULLIFY(PblH   )
    NULLIFY(Fsds   )
    NULLIFY(PRain  )
    NULLIFY(CldTop )
    NULLIFY(CldFrc )

    ! Eventually initialize/reset wetdep
    IF ( Input_Opt%LConv .OR. Input_Opt%LChem .OR. Input_Opt%LWetD ) THEN
        CALL Setup_WetScav( am_I_Root  = rootChunk,         &
                            Input_Opt  = Input_Opt,         &
                            State_Chm  = State_Chm(LCHNK),  &
                            State_Grid = State_Grid(LCHNK), &
                            State_Met  = State_Met(LCHNK),  &
                            RC         = RC                )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Failed to set up wet scavenging!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF

    ! Determine current date and time
    CALL Get_Curr_Date(CurrYr,CurrMo,CurrDy,CurrTOD)
    ! For now, force year to be 2000
    CurrYr  = 2000
    CurrYMD = (CurrYr*1000) + (CurrMo*100) + (CurrDy)
    ! Deal with subdaily
    CurrUTC = REAL(CurrTOD,f4)/3600.0e+0_f4
    CurrSc  = 0
    CurrMn  = 0
    CurrHr  = 0
    DO WHILE (CurrTOD > 3600)
        CurrTOD = CurrTOD - 3600
        CurrHr  = CurrHr + 1
    ENDDO
    DO WHILE (CurrTOD > 60)
        CurrTOD = CurrTOD - 60
        CurrMn  = CurrMn + 1
    ENDDO
    CurrSc  = CurrTOD
    CurrHMS = (CurrHr*1000) + (CurrMn*100) + (CurrSc)


    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( am_I_Root      = rootChunk,  &
                                    value_NYMD     = CurrYMD,            &
                                    value_NHMS     = CurrHMS,            &
                                    value_YEAR     = CurrYr,             &
                                    value_MONTH    = CurrMo,             &
                                    value_DAY      = CurrDy,             &
                                    value_DAYOFYR  = INT(FLOOR(Calday)), &
                                    value_HOUR     = CurrHr,             &
                                    value_MINUTE   = CurrMn,             &
                                    value_HELAPSED = 0.0e+0_f4,    &
                                    value_UTC      = CurrUTC,            &
                                    RC             = RC    )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to update time in GEOS-Chem!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Accept_External_PEdge( am_I_Root = rootChunk,        &
                                State_Met = State_Met(LCHNK), &
                                RC        = RC               )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to update pressure edges!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Calculate State_Met etc for this timestep
    ! Use the CAM psdry fields instead of using the GC calculation
    !CALL Set_Dry_Surface_Pressure(State_Met(LCHNK), 1)
    State_Met(LCHNK)%PS1_DRY (1,:) = State%PSDry(:) * 0.01e+0_fp
    State_Met(LCHNK)%PS2_DRY (1,:) = State%PSDry(:) * 0.01e+0_fp

    ! Set surface pressures to match those in input
    State_Met(LCHNK)%PSC2_WET = State_Met(LCHNK)%PS1_WET
    State_Met(LCHNK)%PSC2_DRY = State_Met(LCHNK)%PS1_DRY
    CALL Set_Floating_Pressures( am_I_Root  = rootChunk,         &
                                 State_Grid = State_Grid(LCHNK), &
                                 State_Met  = State_Met(LCHNK),  &
                                 RC         = RC                )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to set floating pressures!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    ! Set quantities of interest but do not change VMRs
    CALL AirQnt( am_I_Root           = rootChunk,         &
                 Input_Opt           = Input_Opt,         &
                 State_Chm           = State_Chm(LCHNK),  &
                 State_Grid          = State_Grid(LCHNK), &
                 State_Met           = State_Met(LCHNK),  &
                 RC                  = RC,                &
                 Update_Mixing_Ratio = .False. )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Failed to calculate air properties!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF


    ! Initialize strat chem if not already done. This has to be done here because
    ! it needs to have non-zero values in State_Chm%AD, which only happens after
    ! the first call to AirQnt
    IF ( (.not.SCHEM_READY) .and. Input_Opt%LSCHEM ) THEN
        CALL Init_Strat_Chem( am_I_Root  = rootChunk,         &
                              Input_Opt  = Input_Opt,         &
                              State_Chm  = State_Chm(LCHNK),  &
                              State_Met  = State_Met(LCHNK),  &
                              State_Grid = State_Grid(LCHNK), &
                              RC         = RC                )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Could not initialize strat-chem!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
        SCHEM_READY = .True.
    ENDIF

    ! Run chemistry
    IF (Input_Opt%LChem) THEN
        CALL Compute_Overhead_O3( am_I_Root       = rootChunk,                 &
                                  State_Grid      = State_Grid(LCHNK),         &
                                  DAY             = 1,                         &
                                  USE_O3_FROM_MET = Input_Opt%Use_O3_From_Met, &
                                  TO3             = State_Met(LCHNK)%TO3 )

        CALL Do_Chemistry( am_I_Root  = rootChunk,         &
                           Input_Opt  = Input_Opt,         &
                           State_Chm  = State_Chm(LCHNK),  &
                           State_Diag = State_Diag(LCHNK), &
                           State_Grid = State_Grid(LCHNK), &
                           State_Met  = State_Met(LCHNK),  &
                           RC         = RC        )

        IF ( RC /= GC_SUCCESS ) THEN
           ErrMsg = 'Chemistry failed!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF
    ENDIF



    !IF (MasterProc) WRITE(iulog,*) ' --> TEND SIZE: ', size(State%NCOL)
    !IF (MasterProc) WRITE(iulog,'(a,2(x,I6))') ' --> TEND SIDE:  ', lbound(State%NCOL),ubound(State%NCOL)

    ! Make sure State_Chm(lchnk) is back in kg/kg dry!

    ! Reset H2O MMR to the initial value (no chemistry tendency in H2O just
    ! yet)
    State_Chm(LCHNK)%Species(1,:,:,iH2O) = MMR_Beg(:,:,iH2O)

    ! Store unadvected species data
    SlsData = 0.0e+0_r8
    DO N = 1, NSls
        M = MAP2GC_Sls(N)
        IF ( M > 0 ) THEN
            Mass1   = 0.0e+0_r8
            MMR_Min = 1.0e+9_r8
            MMR_Max = 0.0e+0_r8
            DO J = 1, NCOL
                DO K = 1, PVER
                    SlsData(J,PVER+1-K,N) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8)
                ENDDO
            ENDDO
        ENDIF
    ENDDO
    CALL Set_Short_Lived_Species( SlsData, LCHNK, NCOL, Pbuf )

    ! Write diagnostic output
    DO N=1, pcnst
        M = Map2GC(N)
        I = Map2IDX(N)
        IF ( M > 0 ) THEN
            SpcName = TracerNames(I)
            VMR     = 0.0e+0_r8
            Mass0   = 0.0e+0_r8
            Mass1   = 0.0e+0_r8
            DO J = 1, NCOL
                DO K = 1, PVER
                    AirMass         = REAL(State_Met(LCHNK)%AD(1,J,K),r8)
                    MMR0            = MMR_Beg(J,K,M)
                    MMR1            = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8)
                    VMR(J,PVER+1-K) = MMR1 * MWRatio(I)
                    Mass0           = Mass0 + (MMR0*AirMass)
                    Mass1           = Mass1 + (MMR1*AirMass)
                ENDDO
            ENDDO
            CALL OutFld( TRIM(SpcName), VMR(:NCOL,:), NCOL, LCHNK )
        ENDIF
    ENDDO
    DO N = 1, NSls
        SpcName = SlsNames(n)
        VMR = 0.0e+0_r8
        M = Map2GC_Sls(n)
        IF ( M > 0 ) THEN
            DO J = 1, NCOL
                DO K = 1, PVER
                    VMR(J,PVER+1-K) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8) * SLSMWratio(N)
                ENDDO
            ENDDO
            CALL OutFld( TRIM(SpcName), VMR(:NCOL,:), NCOL, LCHNK )
        ENDIF
    ENDDO

    ! NOTE: Re-flip all the arrays vertically or suffer the consequences
    ! ptend%q dimensions: [column, ?, species]
    Ptend%Q(:,:,:) = 0.0e+0_r8
    MMR_End = 0.0e+0_r8
    DO N = 1, pcnst
        M = Map2GC(N)
        IF (M > 0) THEN
            I = 1
            DO J = 1, NCOL
                DO K = 1, PVER
                    ! CURRENTLY KG/KG
                    MMR_End (J,K,M) = REAL(State_Chm(LCHNK)%Species(1,J,K,M),r8)
                    MMR_TEnd(J,K,M) = MMR_End(J,K,M) - MMR_Beg(J,K,M)
                    ptend%q(J,PVER+1-K,N) = (MMR_End(J,K,M)-MMR_Beg(J,K,M))/dT
                ENDDO
            ENDDO
        ENDIF
    ENDDO

    IF (PRESENT(fh2o)) THEN
        fh2o(:NCOL) = 0.0e+0_r8
        !DO K = 1, PVER
        !   fh2o(:NCOL) = fh2o(:NCOL) + Ptend%Q(:NCOL,K,iH2O)*State%Pdel(:NCOL,K)/Gravit
        !ENDDO
    ENDIF

    IF (rootChunk) WRITE(iulog,'(a)') ' GEOS-Chem chemistry step completed'
    RETURN

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

    INTEGER  :: ILEV, NLEV, I
    REAL(r8) :: QTemp, Min_MMR

    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_INIT_CNST'

    NLEV = SIZE(Q, 2)
    ! Retrieve a "background value" for this from the database
    Min_MMR = 1.0e-38_r8
    DO I = 1, NTracers
        IF (TRIM(TracerNames(I)).eq.TRIM(name)) THEN
            Min_MMR = Ref_MMR(i)
            EXIT
        ENDIF
    ENDDO

    DO ILEV=1,NLEV
       WHERE(MASK)
          ! Set to the minimum mixing ratio
          Q(:,ILEV) = Min_MMR
       ENDWHERE
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
    use WetScav_Mod,    only : Cleanup_Wetscav
    use Carbon_Mod,     only : Cleanup_Carbon
    use Dust_Mod,       only : Cleanup_Dust
    use Seasalt_Mod,    only : Cleanup_Seasalt
    use Aerosol_Mod,    only : Cleanup_Aerosol
    use TOMS_Mod,       only : Cleanup_Toms
    use Sulfate_Mod,    only : Cleanup_Sulfate
    use Pressure_Mod,   only : Cleanup_Pressure
    use Strat_Chem_Mod, only : Cleanup_Strat_Chem
    use PBL_Mix_Mod,    only : Cleanup_PBL_Mix

    use CMN_Size_Mod,   only : Cleanup_CMN_Size
    use CMN_O3_Mod,     only : Cleanup_CMN_O3
    use CMN_FJX_Mod,    only : Cleanup_CMN_FJX

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

    CALL Cleanup_PBL_Mix
    CALL Cleanup_Pressure
    CALL Cleanup_Seasalt
    CALL Cleanup_Sulfate
    CALL Cleanup_Strat_Chem
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

    ! Call extra cleanup routines, from modules in Headers/
    CALL Cleanup_CMN_O3( MasterProc, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_O3"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_CMN_SIZE( MasterProc, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_SIZE"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    CALL Cleanup_CMN_FJX( MasterProc, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Cleanup_CMN_FJX"!'
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

    IF (ALLOCATED(Slvd_Lst    ))  DEALLOCATE(Slvd_Lst)
    IF (ALLOCATED(Slvd_Ref_MMR))  DEALLOCATE(Slvd_Ref_MMR)

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

    INTEGER :: LCHNK, NCOL
    LOGICAL :: rootChunk

    ! LCHNK: which chunk we have on this process
    LCHNK = State%LCHNK
    ! NCOL: number of atmospheric columns on this chunk
    NCOL  = State%NCOL
    rootChunk = ( MasterProc.and.(LCHNK.EQ.BEGCHUNK) )

    IF (rootChunk) WRITE(iulog,'(a)') 'GCCALL CHEM_EMISSIONS'

  end subroutine chem_emissions

end module chemistry
