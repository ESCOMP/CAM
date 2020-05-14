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
  USE Input_Opt_Mod,       ONLY : OptInput   ! Derived type for Input Options
  USE State_Chm_Mod,       ONLY : ChmState   ! Derived type for Chemistry State object
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

  ! GEOS-Chem state variables
  TYPE(OptInput),ALLOCATABLE      :: Input_Opt(:)   ! Input Options object
  TYPE(ChmState),ALLOCATABLE      :: State_Chm(:)   ! Chemistry State object
  TYPE(GrdState),ALLOCATABLE      :: State_Grid(:)  ! Grid State object
  TYPE(MetState),ALLOCATABLE      :: State_Met(:)   ! Meteorology State object

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
    use GC_Environment_Mod
    use GC_Grid_Mod,    only : SetGridFromCtrEdges

    ! Use GEOS-Chem versions of physical constants
    use PhysConstants,  only : PI, PI_180

    use Time_Mod,      only : Accept_External_Date_Time
    !use Time_Mod,      only : Set_Begin_Time,   Set_End_Time
    !use Time_Mod,      only : Set_Current_Time, Set_DiagB
    use Transfer_Mod,  only : Init_Transfer
    use Linoz_Mod,     only : Linoz_Read

    TYPE(physics_state), INTENT(IN):: phys_state(BEGCHUNK:ENDCHUNK)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    ! Local variables
    INTEGER :: LCHNK(BEGCHUNK:ENDCHUNK), NCOL(BEGCHUNK:ENDCHUNK)
    INTEGER               :: IWAIT, IERR

    INTEGER               :: NX, NY, NZ
    INTEGER               :: NLEV, I, J, L, RC
    INTEGER               :: NLINOZ

    ! Grid setup
    REAL(fp)              :: lonVal,  latVal
    REAL(fp)              :: dLonFix, dLatFix
    REAL(f4), ALLOCATABLE :: lonMidArr(:,:),  latMidArr(:,:)
    REAL(f4), ALLOCATABLE :: lonEdgeArr(:,:), latEdgeArr(:,:)
    REAL(r8), ALLOCATABLE :: linozData(:,:,:,:)

    LOGICAL :: am_I_Root

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

    ! The GEOS-Chem grids on every "chunk" will all be the same size, to avoid
    ! the possibility of having differently-sized chunks
    NX = 1
    NY = MAXVAL(NCOL)
    NZ = NLEV

    ! This ensures that each process allocates everything needed for its chunks
    ALLOCATE(Input_Opt(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating Input_Opt')
    ALLOCATE(State_Chm(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Chm')
    ALLOCATE(State_Grid(BEGCHUNK:ENDCHUNK), STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Grid')
    ALLOCATE(State_Met(BEGCHUNK:ENDCHUNK) , STAT=IERR)
    IF ( IERR .NE. 0 ) CALL ENDRUN('Failure while allocating State_Met')

    DO I = BEGCHUNK, ENDCHUNK

        ! Only treat the first chunk as the "root" CPU
        am_I_Root = ((I.eq.BEGCHUNK) .and. MasterProc)

        ! Set some basic flags
        Input_Opt(I)%Max_BPCH_Diag     = 1000
        Input_Opt(I)%Max_AdvectSpc     = 500
        Input_Opt(I)%Max_Families      = 250

        Input_Opt(I)%Linoz_NLevels     = 25
        Input_Opt(I)%Linoz_NLat        = 18
        Input_Opt(I)%Linoz_NMonths     = 12
        Input_Opt(I)%Linoz_NFields     = 7
        Input_Opt(I)%RootCPU           = MasterProc

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
                           Input_Opt  = Input_Opt(I),   &
                           State_Grid = State_Grid(I),  &
                           RC         = RC          )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered within call to "GC_Init_Grid"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        ! Set default values
        CALL Set_Input_Opt( am_I_Root, Input_Opt(I), RC )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered within call to "Set_Input_Opt"!'
           CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        ! Each Input_Opt object should be indexed independently, but
        ! for now, we can just associate all of them with the same CPU
        Input_Opt(I)%myCPU = myCPU

    ENDDO

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
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%Chem_Inputs_Dir      = TRIM(chemInputsDir)
    ENDDO

    DO I = BEGCHUNK, ENDCHUNK
        ! Simulation menu
        Input_Opt(I)%NYMDb                = 20000101
        Input_Opt(I)%NHMSb                =   000000
        Input_Opt(I)%NYMDe                = 20010101
        Input_Opt(I)%NHMSe                =   000000
    ENDDO

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

        CALL Init_Transfer( State_Grid(L), 0, 0 )
    ENDDO
    DEALLOCATE(lonMidArr)
    DEALLOCATE(latMidArr)
    DEALLOCATE(lonEdgeArr)
    DEALLOCATE(latEdgeArr)


    ! Now READ_SIMULATION_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%ITS_A_CH4_SIM          = .False.
        Input_Opt(I)%ITS_A_CO2_SIM          = .False.
        Input_Opt(I)%ITS_A_FULLCHEM_SIM     = .True.
        Input_Opt(I)%ITS_A_MERCURY_SIM      = .False.
        Input_Opt(I)%ITS_A_POPS_SIM         = .False.
        Input_Opt(I)%ITS_A_RnPbBe_SIM       = .False.
        Input_Opt(I)%ITS_A_TAGO3_SIM        = .False.
        Input_Opt(I)%ITS_A_TAGCO_SIM        = .False.
        Input_Opt(I)%ITS_AN_AEROSOL_SIM     = .False.
    ENDDO

    ! Now READ_ADVECTED_SPECIES_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%N_Advect               = NTracers
        IF (Input_Opt(I)%N_Advect.GT.Input_Opt(I)%Max_AdvectSpc) THEN
            CALL ENDRUN('Number of tracers exceeds max count')
        ENDIF
        ! Assign tracer names
        DO J = 1, Input_Opt(I)%N_Advect
            Input_Opt(I)%AdvectSpc_Name(J) = TRIM(TRACERNAMES(J))
        ENDDO
        ! No tagged species
        Input_Opt(I)%LSplit = .False.
    ENDDO

    ! Now READ_TRANSPORT_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LTran                  = .True.
        Input_Opt(I)%LFill                  = .True.
        Input_Opt(I)%TPCore_IOrd            = 3
        Input_Opt(I)%TPCore_JOrd            = 3
        Input_Opt(I)%TPCore_KOrd            = 3
    ENDDO

    ! Now READ_CONVECTION_MENU
    ! For now, TMMF
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LConv                  = .False.
        Input_Opt(I)%LTurb                  = .False.
        Input_Opt(I)%LNLPBL                 = .False.
    ENDDO

    ! Now READ_EMISSIONS_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LEmis                  = .False.
        Input_Opt(I)%HCOConfigFile          = 'HEMCO_Config.rc'
        Input_Opt(I)%LFix_PBL_Bro           = .False.

        ! Set surface VMRs - turn this off so that CAM can handle it
        Input_Opt(I)%LCH4Emis               = .False.
        Input_Opt(I)%LCH4SBC                = .False.
        Input_Opt(I)%LOCSEmis               = .False.
        Input_Opt(I)%LCFCEmis               = .False.
        Input_Opt(I)%LClEmis                = .False.
        Input_Opt(I)%LBrEmis                = .False.
        Input_Opt(I)%LN2OEmis               = .False.
        Input_Opt(I)%LBasicEmis             = .False.

        ! Set initial conditions
        Input_Opt(I)%LSetH2O                = .True.

        ! CFC control
        Input_Opt(I)%CFCYear                = 0
    ENDDO

    ! Now READ_AEROSOL_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LSulf               = .True.
        Input_Opt(I)%LMetalcatSO2        = .True.
        Input_Opt(I)%LCarb               = .True.
        Input_Opt(I)%LBrC                = .False.
        Input_Opt(I)%LSOA                = .True.
        Input_Opt(I)%LSVPOA              = .False.
        Input_Opt(I)%LOMOC               = .False.
        Input_Opt(I)%LDust               = .True.
        Input_Opt(I)%LDstUp              = .False.
        Input_Opt(I)%LSSalt              = .True.
        Input_Opt(I)%SalA_rEdge_um(1)    = 0.01e+0_fp
        Input_Opt(I)%SalA_rEdge_um(2)    = 0.50e+0_fp
        Input_Opt(I)%SalC_rEdge_um(1)    = 0.50e+0_fp
        Input_Opt(I)%SalC_rEdge_um(2)    = 8.00e+0_fp
        Input_Opt(I)%LMPOA               = .False.
        Input_Opt(I)%LGravStrat          = .True.
        Input_Opt(I)%LSolidPSC           = .True.
        Input_Opt(I)%LHomNucNAT          = .False.
        Input_Opt(I)%T_NAT_Supercool     = 3.0e+0_fp
        Input_Opt(I)%P_Ice_Supersat      = 1.2e+0_fp
        Input_Opt(I)%LPSCChem            = .True.
        Input_Opt(I)%LStratOD            = .True.
        Input_Opt(I)%hvAerNIT            = .False.
        Input_Opt(I)%hvAerNIT_JNIT       = .False.
        Input_Opt(I)%hvAerNIT_JNITs      = .False.
        Input_Opt(I)%JNITChanA           = 0e+0_fp
        Input_Opt(I)%JNITChanB           = 0e+0_fp
    ENDDO

    ! Now READ_DEPOSITION_MENU
    ! Disable dry/wet dep for now
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LDryD                  = .False.
        Input_Opt(I)%LWetD                  = .False.
        Input_Opt(I)%CO2_Effect             = .False.
        Input_Opt(I)%CO2_Level              = 390.0_fp
        Input_Opt(I)%CO2_Ref                = 390.0_fp
    ENDDO

    ! Now READ_CHEMISTRY_MENU
    DO I = BEGCHUNK, ENDCHUNK
        Input_Opt(I)%LChem                  = .True.
        Input_Opt(I)%LSChem                 = .True.
        Input_Opt(I)%LLinoz                 = .True.
        Input_Opt(I)%LSynoz                 = .True.
        Input_Opt(I)%LUCX                   = .True.
        Input_Opt(I)%LActiveH2O             = .True.
        Input_Opt(I)%Use_Online_O3          = .True.
        ! Expect to get total overhead ozone, although it shouldn not
        ! make too much of a difference since we want to use "full-UCX"
        Input_Opt(I)%Use_O3_from_Met        = .True.
        Input_Opt(I)%Use_TOMS_O3            = .False.
        Input_Opt(I)%Gamma_HO2              = 0.2e+0_fp
    ENDDO

    ! Read in data for Linoz. All CPUs allocate one array to hold the data. Only
    ! the root CPU reads in the data; then we copy it out to a temporary array,
    ! broadcast to all other CPUs, and finally duplicate the data into every
    ! copy of Input_Opt
    IF (Input_Opt(BEGCHUNK)%LLinoz) THEN
        ! Allocate array for broadcast
        nLinoz = Input_Opt(BEGCHUNK)%Linoz_NLevels * &
                 Input_Opt(BEGCHUNK)%Linoz_NLat    * &
                 Input_Opt(BEGCHUNK)%Linoz_NMonths * &
                 Input_Opt(BEGCHUNK)%Linoz_NFields 
        ALLOCATE( linozData( Input_Opt(BEGCHUNK)%Linoz_NLevels,     &
                             Input_Opt(BEGCHUNK)%Linoz_NLat,        &
                             Input_Opt(BEGCHUNK)%Linoz_NMonths,     &
                             Input_Opt(BEGCHUNK)%Linoz_NFields  ), STAT=IERR)
        IF (IERR.NE.0) CALL ENDRUN('Failure while allocating linozData')
        linozData = 0.0e+0_r8

        IF ( MasterProc ) THEN
            ! Read data in to Input_Opt%Linoz_TParm
            CALL Linoz_Read( MasterProc, Input_Opt(BEGCHUNK), RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in "Linoz_Read"!'
               CALL Error_Stop( ErrMsg, ThisLoc )
            ENDIF
            ! Copy the data to a temporary array
            linozData = REAL(Input_Opt(BEGCHUNK)%LINOZ_TPARM,r8)
        ENDIF
#if defined( SPMD )
        CALL MPIBCAST(linozData, nLinoz, MPIR8, 0, MPICOM )
#endif
        ! Now copy the data to all other Input_Opt copies
        DO I=BEGCHUNK,ENDCHUNK
        Input_Opt(I)%LINOZ_TPARM = REAL(linozData,fp)
        ENDDO
        DEALLOCATE(linozData)
    ENDIF

    ! Note - this is called AFTER chem_readnl, after X, and after
    ! every constituent has had its initial conditions read. Any
    ! constituent which is not found in the CAM restart file will
    ! then have already had a call to chem_implements_cnst, and will
    ! have then had a call to chem_init_cnst to set a default VMR
    ! Call the routine GC_Allocate_All (located in module file
    ! GeosCore/gc_environment_mod.F90) to allocate all lat/lon
    ! allocatable arrays used by GEOS-Chem.
    CALL GC_Allocate_All ( am_I_Root      = MasterProc,           &
                           Input_Opt      = Input_Opt(BEGCHUNK),  &
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

    ! Set the times held by "time_mod"
    CALL Accept_External_Date_Time( am_I_Root   = MasterProc,                &
                                    value_NYMDb = Input_Opt(BEGCHUNK)%NYMDb, &
                                    value_NHMSb = Input_Opt(BEGCHUNK)%NHMSb, &
                                    value_NYMDe = Input_Opt(BEGCHUNK)%NYMDe, &
                                    value_NHMSe = Input_Opt(BEGCHUNK)%NHMSe, &
                                    value_NYMD  = Input_Opt(BEGCHUNK)%NYMDb, &
                                    value_NHMS  = Input_Opt(BEGCHUNK)%NHMSb, &
                                    RC          = RC                )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Accept_External_Date_Time"!'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    DO I = BEGCHUNK, ENDCHUNK
        ! Start by setting some dummy timesteps
        CALL GC_Update_Timesteps(I,300.0E+0_r8)
    ENDDO

    ! Initialize the state objects for each chunk
    DO I = BEGCHUNK, ENDCHUNK
        ! This reproduces GC_Init_Stateobj without the State_Diag object
        am_I_Root = (MasterProc .AND. (I == BEGCHUNK))
        CALL Init_State_Met( am_I_Root, State_Grid(I), State_Met(I), RC )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_State_Met"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        CALL Init_State_Chm( am_I_Root, Input_Opt(I),     &
                             State_Chm(I), State_Grid(I), &
                             RC )
        IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_State_Chm"!'
            CALL Error_Stop( ErrMsg, ThisLoc )
        ENDIF

        ! Now replicate GC_Init_Extra...

        ! Start with v/v dry (CAM standard)
        State_Chm(I)%Spc_Units = 'v/v dry'
    ENDDO
    ! Init_FJX..
    ! Init_Pressure...
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

  subroutine GC_Update_Timesteps(LCHNK,DT)

    use Time_Mod,       only : Set_Timesteps
    use cam_abortutils, only : endrun

    INTEGER,  INTENT(IN) :: LCHNK
    REAL(r8), INTENT(IN) :: DT
    INTEGER              :: DT_MIN
    INTEGER, SAVE        :: DT_MIN_LAST = -1

    DT_MIN = NINT(DT)

    Input_Opt(LCHNK)%TS_CHEM = DT_MIN
    Input_Opt(LCHNK)%TS_EMIS = DT_MIN
    Input_Opt(LCHNK)%TS_CONV = DT_MIN
    Input_Opt(LCHNK)%TS_DYN  = DT_MIN
    Input_Opt(LCHNK)%TS_RAD  = DT_MIN

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
    CALL GC_Update_Timesteps(LCHNK,DT)

    ! Need to be super careful that the module arrays are updated and correctly
    ! set

    ! 1. Update State_Met etc for this timestep

    !if (MasterProc) WRITE(iulog,*) ' --> TEND SIZE: ', size(state%ncol)
    !if (MasterProc) WRITE(iulog,'(a,2(x,I6))') ' --> TEND SIDE:  ', lbound(state%ncol),ubound(state%ncol)


    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_TEND'
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

    USE Input_Opt_Mod, ONLY : CLEANUP_INPUT_OPT
    USE State_Chm_Mod, ONLY : CLEANUP_STATE_CHM
    USE State_Met_Mod, ONLY : CLEANUP_STATE_MET
    Use UCX_Mod,       ONLY : CLEANUP_UCX

    INTEGER :: I, RC
    LOGICAL :: am_I_Root

    ! Finalize GEOS-Chem
    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_FINAL'
    
    ! Loop over each chunk and cleanup the independent variables
    DO I = BEGCHUNK, ENDCHUNK
        am_I_Root = ((I.eq.BEGCHUNK) .and. MasterProc)
        CALL CLEANUP_INPUT_OPT( am_I_Root, Input_Opt(I), RC )
        CALL CLEANUP_STATE_MET( am_I_Root, State_Met(I), RC )
        CALL CLEANUP_STATE_CHM( am_I_Root, State_Chm(I), RC )
    ENDDO

    ! Finally deallocate state variables
    IF (ALLOCATED(Input_Opt))     DEALLOCATE(Input_Opt)
    IF (ALLOCATED(State_Chm))     DEALLOCATE(State_Chm)
    IF (ALLOCATED(State_Grid))    DEALLOCATE(State_Grid)
    IF (ALLOCATED(State_Met))     DEALLOCATE(State_Met)

    IF (MasterProc) WRITE(iulog,'(a,4(x,L1))') ' --> DEALLOC CHECK : ', ALLOCATED(Input_Opt), ALLOCATED(State_Chm), ALLOCATED(State_Grid), ALLOCATED(State_Met)

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
