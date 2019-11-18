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
  use spmd_utils,          only : masterproc
  use cam_logfile,         only : iulog

  use Input_Opt_Mod,       only : OptInput
  use State_Met_Mod,       only : MetState
  use State_Chm_Mod,       only : ChmState

  implicit none
  private
  save
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
  character(len=500) :: inputGeosPath

  ! GEOS-Chem state variables
  Type(OptInput),Allocatable      :: Input_Opt(:)
  Type(MetState),Allocatable      :: State_Met(:)
  Type(ChmState),Allocatable      :: State_Chm(:)

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'geoschem' ) then
       chem_is = .true.
    end if
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_IS'

  end function chem_is

!================================================================================================

  subroutine chem_register

    use physics_buffer, only : pbuf_add_field, dtype_r8
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents for chemistry
    !
    !-----------------------------------------------------------------------
    integer            :: i, n
    real(r8)           :: cptmp
    real(r8)           :: qmin
    character(len=128) :: mixtype
    character(len=128) :: molectype
    character(len=128) :: lng_name
    logical            :: camout
    logical            :: ic_from_cam2
    logical            :: has_fixed_ubc
    logical            :: has_fixed_ubflx
    ! SDE 2018-05-02: This seems to get called before anything else
    ! That includes CHEM_INIT
    ! At this point, mozart calls SET_SIM_DAT, which is specified by each
    ! mechanism separately (ie mozart/chemistry.F90 calls the subroutine
    ! set_sim_dat which is in pp_[mechanism]/mo_sim_dat.F90. That sets a lot of
    ! data in other places, notably in "chem_mods"

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_REGISTER'
    ! At the moment, we force nadv_chem=200 in the setup file
    ntracers = 200
    do i = 1, ntracersmax
       ! TODO: Read input.geos in chem_readnl to get tracernames(1:ntracers)
       ! TODO: Get all other species properties here from species database
       ! Hard-code for now
       select case (tracernames(i))
          case ('BCPI')
             lng_name = 'Hydrophilic black carbon'
             ! Molar mass (g/mol)
             adv_mass(i) = 1000.0e+0_r8 * (0.012e+0_r8)
          case ('OCS')
             lng_name = 'Carbonyl sulfide'
             ! Molar mass (g/mol)
             adv_mass(i) = 1000.0e+0_r8 * (0.060e+0_r8)
          case default
             lng_name = tracernames(i)
             adv_mass(i) = 1000.0e+0_r8 * (0.001e+0_r8)
       end select
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
       call cnst_add( trim(tracernames(i)), adv_mass(i), cptmp, qmin, n, &
                      readiv=ic_from_cam2, mixtype=mixtype, cam_outfld=camout, &
                      molectype=molectype, fixed_ubc=has_fixed_ubc, &
                      fixed_ubflx=has_fixed_ubflx, longname=trim(lng_name) )
    end do

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
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: i, unitn, ierr
    character(len=500) :: line
    logical :: menuFound, validSLS

    inputGeosPath='/n/scratchlfs/jacob_lab/elundgren/UT/runs/4x5_standard/input.geos.template'

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_READNL'

    ! TODO: Read in input.geos and get species names
    if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(inputGeosPath), status='old', iostat=ierr )
        if (ierr .ne. 0) then
            call endrun('chem_readnl: ERROR opening input.geos')
        end if
        ! Go to ADVECTED SPECIES MENU
        menuFound = .False.
        Do While (.not.menuFound)
            read( unitn, '(a)', iostat=ierr ) line
            if (ierr.ne.0) then
                call endrun('chem_readnl: ERROR finding advected species menu')
            else  if (index(line,'ADVECTED SPECIES MENU') > 0) then
                menuFound=.True.
       end if
    end do
        ! Skip first line
        read(unitn,'(a)',iostat=ierr) line
        ! Read in tracer count
        read(unitn,'(26x,I)',iostat=ierr) ntracers
        ! Skip divider line
        read(unitn,'(a)',iostat=ierr) line
        ! Read in each tracer
        do i=1,ntracers
            read(unitn,'(26x,a)',iostat=ierr) line
            tracernames(i) = trim(line)
        end do
        close(unitn)
        call freeunit(unitn)
        ! Assign remaining tracers dummy names
    do i=(ntracers+1),ntracersmax
       write(tracernames(i),'(a,I0.4)') 'GCTRC_',i
    end do

        ! Now go through the KPP mechanism and add any species not implemented by
        ! the tracer list in input.geos
        if ( nspec > nslsmax ) then
            call endrun('chem_readnl: too many species - increase nslsmax')
        end If
        nsls = 0
        do i=1,nspec
            ! Get the name of the species from KPP
            line = adjustl(trim(spc_names(i)))
            ! Only add this
            validSLS = ( (.not.any(trim(line).eq.tracernames)).and.&
                         (.not.(line(1:2) == 'RR')) )
            if (validSLS) then
                ! Genuine new short-lived species
                nsls = nsls + 1
                slsnames(nsls) = trim(line)
                write(iulog,'(a,I5,a,a)') ' --> GC species ',nsls, ': ',trim(line)
            end if
        end do
    end if

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
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_IS_ACTIVE'
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
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    integer :: i

    chem_implements_cnst = .false.

    do i = 1, ntracers
       if (trim(tracernames(i)) .eq. trim(name)) then
          chem_implements_cnst = .true.
          exit
       end if
    end do

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_IMPLEMENTS_CNST'
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

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: lchnk(begchunk:endchunk), ncol(begchunk:endchunk)
    integer :: iwait

    integer :: iipar, jjpar, llpar
    integer :: nlev, i

    ! lchnk: which chunks we have on this process
    lchnk = phys_state%lchnk
    ! ncol: number of atmospheric columns for each chunk
    ncol  = phys_state%ncol
    ! nlev: number of vertical levels
    nlev  = pver

    ! This ensures that each process allocates everything needed for its chunks
    if (.not.allocated(Input_Opt)) then
        Allocate(Input_Opt(begchunk:endchunk))
        Allocate(State_Met(begchunk:endchunk))
        Allocate(State_Chm(begchunk:endchunk))
        if (masterproc) write(iulog,'(a,3(x,L1))') ' --> ALLOC CHECK   : ', Allocated(Input_Opt), Allocated(State_Met), Allocated(State_Chm)
    end if
    write(iulog,'(a,x,L1,2(x,I6))') ' --> SIZE  CHECK   : ', masterproc,lbound(Input_Opt),ubound(Input_Opt)

    Do i = begchunk, endchunk

        ! Set some basic flags
        Input_Opt(i)%Max_BPCH_Diag     = 1000
        Input_Opt(i)%Max_AdvectSpc     = 500
        Input_Opt(i)%Max_Families      = 250

        Input_Opt(i)%Linoz_NLevels     = 25
        Input_Opt(i)%Linoz_NLat        = 18
        Input_Opt(i)%Linoz_NMonths     = 12
        Input_Opt(i)%Linoz_NFields     = 7
        Input_Opt(i)%RootCPU           = ((i.eq.begchunk) .and. MasterProc)

        IIPAR = 1
        JJPAR = ncol(i)
        LLPAR = nlev

    end do

    ! Can add history output here too with the "addfld" & "add_default" routines
    ! Note that constituents are already output by default
    call addfld ( 'BCPI', (/'lev'/), 'A', 'mole/mole', trim('BCPI')//' mixing ratio' )
    call add_default ( 'BCPI',   1, ' ')
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT'

  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_INIT'

    ! This is when we want to update State_Met and so on
    ! Note that here we have been passed MANY chunks

  end subroutine chem_timestep_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o )

    use physics_buffer,   only: physics_buffer_desc
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t

    real(r8),            intent(in)    :: dt          ! time step
    type(physics_state), intent(in)    :: state       ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), optional,  intent(out)   :: fh2o(pcols) ! h2o flux to balance source from chemistry

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

    !if (masterproc) write(iulog,*) ' --> TEND SIZE: ', size(state%ncol)
    !if (masterproc) write(iulog,'(a,2(x,I6))') ' --> TEND SIDE:  ', lbound(state%ncol),ubound(state%ncol)


    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_TIMESTEP_TEND'
    lq(:) = .false.
    do n=1,pcnst
        !m = map2chm(n)
        m=0
        if (m > 0) lq(n) = .true.
    end do
    call physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)
    ! ptend%q dimensions: [column, ?, species]
    !ptend%q(:ncol,:,:) = 0.0e+0_r8
    !ptend%q(:ncol,:,:) = 0.0e+0_r8
    ptend%q(:,:,:) = 0.0e+0_r8
    if (present(fh2o)) fh2o(:) = 0.0e+0_r8

    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

    character(len=*), intent(in)  :: name       !  constituent name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (ncol, plev
    ! Used to initialize tracer fields if desired.
    ! Will need a simple mapping structure as well as the CAM tracer registration
    ! routines.

    integer :: ilev, nlev

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT_CNST'

    nlev = size(q, 2)
    if ( any( tracernames .eq. name ) ) then
       do ilev=1,nlev
          where(mask)
             ! Set to the minimum mixing ratio
             q(:,ilev) = 1.0e-38_r8
          end where
       end do
    end if

  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final

    ! Finalize GEOS-Chem
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_FINAL'
    If (allocated(Input_Opt))     Deallocate(Input_Opt)
    If (allocated(State_Met))     Deallocate(State_Met)
    If (allocated(State_Chm))     Deallocate(State_Chm)
    if (masterproc) write(iulog,'(a,3(x,L1))') ' --> DEALLOC CHECK : ', Allocated(Input_Opt), Allocated(State_Met), Allocated(State_Chm)

    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_INIT_RESTART'
    return
  end subroutine chem_init_restart
!===============================================================================
  subroutine chem_write_restart( File )
    !use tracer_cnst, only: write_tracer_cnst_restart
    !use tracer_srcs, only: write_tracer_srcs_restart
    !use linoz_data,  only: write_linoz_data_restart
    use pio, only : file_desc_t
    implicit none
    type(file_desc_t) :: File

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_WRITE_RESTART'
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
    implicit none
    type(file_desc_t) :: File

    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_READ_RESTART'
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

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state
    if (masterproc) write(iulog,'(a)') 'GCCALL CHEM_EMISSIONS'

  end subroutine chem_emissions

end module chemistry
