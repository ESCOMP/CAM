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
  use spmd_utils,          only : MasterProc
  use cam_logfile,         only : iulog

  use Input_Opt_Mod,       only : OptInput
  use State_Met_Mod,       only : MetState
  use State_Chm_Mod,       only : ChmState

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

  ! GEOS-Chem state variables
  Type(OptInput),Allocatable      :: Input_Opt(:)
  Type(MetState),Allocatable      :: State_Met(:)
  Type(ChmState),Allocatable      :: State_Chm(:)

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
    NTRACERS = 200
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

    inputGeosPath='/n/scratchlfs/jacob_lab/elundgren/UT/runs/4x5_standard/input.geos.template'

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

    use Input_Opt_Mod
    use State_Met_Mod
    use State_Chm_Mod
    use GC_Environment_Mod

    TYPE(physics_state), INTENT(IN):: phys_state(BEGCHUNK:ENDCHUNK)
    TYPE(physics_buffer_desc), POINTER :: pbuf2d(:,:)

    INTEGER :: LCHNK(BEGCHUNK:ENDCHUNK), NCOL(BEGCHUNK:ENDCHUNK)
    INTEGER :: IWAIT

    INTEGER :: IIPAR, JJPAR, LLPAR
    INTEGER :: NLEV, I, RC

    LOGICAL :: am_I_Root

    ! lchnk: which chunks we have on this process
    LCHNK = PHYS_STATE%LCHNK
    ! ncol: number of atmospheric columns for each chunk
    NCOL  = PHYS_STATE%NCOL
    ! nlev: number of vertical levels
    NLEV  = PVER

    ! This ensures that each process allocates everything needed for its chunks
    IF (.NOT.ALLOCATED(Input_Opt)) THEN
        ALLOCATE(Input_Opt(BEGCHUNK:ENDCHUNK))
        ALLOCATE(State_Met(BEGCHUNK:ENDCHUNK))
        ALLOCATE(State_Chm(BEGCHUNK:ENDCHUNK))
        IF (MasterProc) WRITE(iulog,'(a,3(x,L1))') ' --> ALLOC CHECK   : ', ALLOCATED(Input_Opt), ALLOCATED(State_Met), ALLOCATED(State_Chm)
    ENDIF
    WRITE(iulog,'(a,x,L1,2(x,I6))') ' --> SIZE  CHECK   : ', MasterProc, LBOUND(Input_Opt), UBOUND(Input_Opt)

    DO I = BEGCHUNK, ENDCHUNK

        ! Only treat the first chunk as the "root" CPU
        am_I_Root = ((I.eq.BEGCHUNK) .and. MasterProc)

        ! Set some basic flags
        Input_Opt(i)%Max_BPCH_Diag     = 1000
        Input_Opt(i)%Max_AdvectSpc     = 500
        Input_Opt(i)%Max_Families      = 250

        Input_Opt(i)%Linoz_NLevels     = 25
        Input_Opt(i)%Linoz_NLat        = 18
        Input_Opt(i)%Linoz_NMonths     = 12
        Input_Opt(i)%Linoz_NFields     = 7
        Input_Opt(i)%RootCPU           = am_I_Root

        IIPAR = 1
        JJPAR = NCOL(I)
        LLPAR = NLEV

    ENDDO

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
    REAL(r8),         INTENT(OUT) :: q(:,:)     ! kg tracer/kg dry air (ncol, plev
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

    INTEGER :: I, RC
    LOGICAL :: am_I_Root

    ! Finalize GEOS-Chem
    IF (MasterProc) WRITE(iulog,'(a)') 'GCCALL CHEM_FINAL'
    
    ! Loop over each chunk
    DO I = BEGCHUNK, ENDCHUNK
        am_I_Root = ((I.eq.BEGCHUNK) .and. MasterProc)
        CALL CLEANUP_INPUT_OPT( am_I_Root, Input_Opt(i), RC )
    ENDDO

    ! Finally deallocate state variables
    IF (ALLOCATED(Input_Opt))     DEALLOCATE(Input_Opt)
    IF (ALLOCATED(State_Met))     DEALLOCATE(State_Met)
    IF (ALLOCATED(State_Chm))     DEALLOCATE(State_Chm)

    IF (MasterProc) WRITE(iulog,'(a,3(x,L1))') ' --> DEALLOC CHECK : ', ALLOCATED(Input_Opt), ALLOCATED(State_Met), ALLOCATED(State_Chm)

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
