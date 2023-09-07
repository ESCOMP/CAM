module frierson_cam
!-----------------------------------------------------------------------
!
! Purpose: Implement idealized forcings described in
!          Frierson, et al. (2006), " A Gray-Radiation Aquaplanet
!          Moist GCM, Part I. Static Stability and Eddy Scale"
!          J. Atmos. Sci, Vol 63, 2548-2566.
!          doi: 10.1175/JAS3753.1
!
!============================================================================
  ! Useful modules
  !-------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: pi => shr_const_pi
  use physconst,      only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir
  use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
  use constituents,   only: pcnst
  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc, &
                            pbuf_set_field, pbuf_get_field
  use camsrfexch,     only: cam_in_t,cam_out_t
  use cam_history,    only: outfld
  use time_manager,   only: is_first_step
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use hycoef,         only: ps0, etamid
  use spmd_utils,     only: mpicom, mstrid=>masterprocid, mpi_real8

  use pio             ,only: file_desc_t, var_desc_t, io_desc_t, pio_double, pio_def_var
  use pio             ,only: pio_write_darray, pio_read_darray, pio_inq_varid
  use cam_grid_support,only: cam_grid_id, cam_grid_dimensions, cam_grid_get_decomp
  use shr_const_mod,   only: SHR_CONST_STEBOL, SHR_CONST_REARTH, SHR_CONST_KARMAN, SHR_CONST_TKTRIP

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !---------------------------------------------------------
  implicit none
  private
  save

  public :: frierson_register
  public :: frierson_readnl
  public :: frierson_init
  public :: frierson_condensate_tend
  public :: frierson_pbl_tend
  public :: frierson_radiative_tend
  public :: frierson_restart_init
  public :: frierson_restart_write
  public :: frierson_restart_read

  private :: frierson_surface_init

  ! PBL Configuatons
  !------------------
  integer,parameter :: PBL_FRIERSON = 0           ! Implementation of Frierson PBL
  integer,parameter :: PBL_USER     = 1           ! Optional call for user defined PBL

  ! Tags to identify optional model formulations
  !------------------------------------------------
  integer,parameter :: CONDENSATE_NONE       = 0  ! No Condensation, PRECL=0
  integer,parameter :: CONDENSATE_FRIERSON   = 1  ! Frierson condensation w/ re-evaporation
  integer,parameter :: CONDENSATE_TJ16       = 2  ! Condensation from TJ2016 model.
  integer,parameter :: CONDENSATE_USER       = 3  ! Optional user defined Condensation scheme

  integer,parameter :: RADIATION_FRIERSON    = 0  ! Frierson Gray radiation.
  integer,parameter :: RADIATION_USER        = 1  ! Optional user defined Radiation scheme

  ! Options selecting which PRECIP, PBL, RADIATION, etc.. formulations to use.
  !---------------------------------------------------------------------------------
  integer,parameter :: PBL_OPT          = PBL_FRIERSON
  integer,parameter :: CONDENSATE_OPT   = CONDENSATE_FRIERSON
  integer,parameter :: RADIATION_OPT    = RADIATION_FRIERSON

  ! Global Constants
  !---------------------
  real(r8),parameter :: frierson_T0     = SHR_CONST_TKTRIP  ! Reference Temperature for E0
  real(r8),parameter :: frierson_E0     = 610.78_r8     ! Saturation Vapor pressure @ T0
  real(r8),parameter :: frierson_Rs0    = 1360.0_r8     ! Solar Constant
  real(r8),parameter :: frierson_Erad   = SHR_CONST_REARTH  ! Earth Radius
  real(r8),parameter :: frierson_Karman = SHR_CONST_KARMAN  ! Von Karman constant
  real(r8),parameter :: frierson_Boltz  = SHR_CONST_STEBOL  ! Stefan-Boltzmann constant

  ! Some Physics buffer indices
  !-------------------------------
  integer :: prec_pcw_idx = 0
  integer :: prec_dp_idx  = 0
  integer :: relhum_idx   = 0

  ! Global values for Surface Temp, surface fluxes, and radiative heating
  !----------------------------------------------------------------------
  type(var_desc_t)    :: Tsurf_desc      ! Vardesc for restarts
  type(var_desc_t)    :: Qsurf_desc      ! Vardesc for restarts
  real(r8),allocatable :: Tsurf (:,:)     ! Surface Temp
  real(r8),allocatable :: Qsurf (:,:)     ! Surface Q
  real(r8),allocatable :: Fsolar(:,:)     ! Net Solar Heating
  real(r8),allocatable :: Fup   (:,:)     ! Upward Longwave flux
  real(r8),allocatable :: Fdown (:,:)     ! Downward Longwave flux
  real(r8),allocatable :: Fup_toa  (:,:)  ! Upward Longwave flux at TOA
  real(r8),allocatable :: Fdown_toa(:,:)  ! Downward Longwave flux at TOA
  real(r8),allocatable :: SHflux(:,:)     ! Sensible Heat flux
  real(r8),allocatable :: LHflux(:,:)     ! Latent Heat Flux
  real(r8),allocatable :: TUflux(:,:)     ! U stress momentum flux
  real(r8),allocatable :: TVflux(:,:)     ! V stress momentum flux
  real(r8),allocatable :: Cd    (:,:)     ! Surface Drag
  real(r8),allocatable :: clat  (:,:)     ! latitudes(radians) for columns
  real(r8),allocatable :: Fnet  (:,:)     ! Net Radiative Surface Heating
  real(r8),allocatable :: Fnet_toa(:,:)   ! Net Radiative Surface Heating at TOA

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Global Tuning values
  !------------------------
  real(r8) :: frierson_Wind_min   = unset_r8      ! Minimum wind threshold
  real(r8) :: frierson_Z0         = unset_r8      ! Roughness Length
  real(r8) :: frierson_Ri_c       = unset_r8      ! Crit. Richardson # for stable mixing
  real(r8) :: frierson_Fb         = unset_r8      ! Surface layer Fraction
  real(r8) :: frierson_Albedo     = unset_r8      ! Frierson Albedo
  real(r8) :: frierson_DeltaS     = unset_r8      ! Lat variation of shortwave radiation
  real(r8) :: frierson_Tau_eqtr   = unset_r8      ! Longwave optical depth at Equator
  real(r8) :: frierson_Tau_pole   = unset_r8      ! Longwave optical depth at poles.
  real(r8) :: frierson_LinFrac    = unset_r8      ! Stratosphere Linear optical depth param
  real(r8) :: frierson_C0         = unset_r8      ! Ocean mixed layer heat capacity
  real(r8) :: frierson_WetDryCoef = unset_r8      ! E0 Scale factor to control moisture
  real(r8) :: frierson_Tmin       = unset_r8      ! IC: Minimum sst (K)
  real(r8) :: frierson_Tdlt       = unset_r8      ! IC: eq-polar difference sst (K)
  real(r8) :: frierson_Twidth     = unset_r8      ! IC: Latitudinal width parameter for sst (degrees latitude)

contains
  !==============================================================================
  subroutine frierson_register()
    !
    ! frierson_register: Register physics buffer values
    !=====================================================================

    call pbuf_add_field('PREC_PCW','physpkg',dtype_r8, (/pcols/),     prec_pcw_idx)
    call pbuf_add_field('PREC_DP' ,'physpkg',dtype_r8, (/pcols/),     prec_dp_idx )
    call pbuf_add_field('RELHUM'  ,'physpkg',dtype_r8, (/pcols,pver/),relhum_idx  )

  end subroutine frierson_register
  !==============================================================================


  !==============================================================================
  subroutine frierson_readnl(nlfile)
    !
    ! frierson_readnl: Read in parameters controlling Frierson parameterizations.
    !=====================================================================
    use namelist_utils,only: find_group_name
    use units         ,only: getunit, freeunit
    !
    ! Passed Variables
    !------------------
    character(len=*),intent(in):: nlfile
    !
    ! Local Values
    !--------------
    integer:: ierr,unitn

    character(len=*), parameter :: sub = 'frierson_readnl'

    namelist /frierson_nl/ frierson_Wind_min, frierson_Z0        , frierson_Ri_c   , &
                           frierson_Fb      , frierson_Albedo    , frierson_DeltaS , &
                           frierson_Tau_eqtr, frierson_Tau_pole  , frierson_LinFrac, &
                           frierson_C0      , frierson_WetDryCoef, frierson_Tmin   , &
                           frierson_Tdlt    , frierson_Twidth

    ! Read in namelist values
    !-------------------------
    if(masterproc) then
      unitn = getunit()
      open(unitn,file=trim(nlfile),status='old')
      call find_group_name(unitn,'frierson_nl',status=ierr)
      if(ierr == 0) then
        read(unitn,frierson_nl,iostat=ierr)
        if(ierr /= 0) then
          call endrun(sub//': ERROR reading namelist')
        endif
      endif
      close(unitn)
      call freeunit(unitn)
    endif

    ! Broadcast namelist values
    !---------------------------
    call mpi_bcast(frierson_Wind_min  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Wind_min")
    call mpi_bcast(frierson_Z0        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Z0")
    call mpi_bcast(frierson_Ri_c      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Ri_c")
    call mpi_bcast(frierson_Fb        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Fb")
    call mpi_bcast(frierson_Albedo    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Albedo")
    call mpi_bcast(frierson_DeltaS    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_DeltaS")
    call mpi_bcast(frierson_Tau_eqtr  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Tau_eqtr")
    call mpi_bcast(frierson_Tau_pole  , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Tau_pole")
    call mpi_bcast(frierson_LinFrac   , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_LinFrac")
    call mpi_bcast(frierson_C0        , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_C0")
    call mpi_bcast(frierson_Tmin      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Tmin")
    call mpi_bcast(frierson_Tdlt      , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Tdlt")
    call mpi_bcast(frierson_Twidth    , 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_Twidth")
    call mpi_bcast(frierson_WetDryCoef, 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_WetDryCoef")

  end subroutine frierson_readnl
  !==============================================================================


  !==============================================================================
  subroutine frierson_init(phys_state,pbuf2d)
    !
    ! frierson_init: allocate space for global arrays and initialize values.
    !                Add variables to history outputs
    !=====================================================================
    use physics_types, only: physics_state
    use error_messages,only: alloc_err
    use cam_history,   only: addfld, add_default,horiz_only
    use phys_grid,     only: get_ncols_p, get_rlat_p
    use frierson,      only: frierson_set_const
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,pointer:: phys_state(:)
    type(physics_buffer_desc),pointer:: pbuf2d    (:,:)
    !
    ! Local Values
    !---------------
    integer :: istat,lchnk,icol,ncol
    real(r8):: adjusted_E0
    real(r8):: frierson_Rs

    ! Initialize constants in Frierson module
    !------------------------------------------
    adjusted_E0 = frierson_WetDryCoef*frierson_E0
    frierson_Rs = frierson_Rs0*(1._r8-frierson_Albedo)
    call frierson_set_const(gravit,cappa,rair,cpair,latvap,rh2o,epsilo,rhoh2o,zvir,ps0,etamid,     &
                            frierson_T0     ,adjusted_E0    ,frierson_Erad    ,frierson_Wind_min,  &
                            frierson_Z0     ,frierson_Ri_c  ,frierson_Karman  ,frierson_Fb      ,  &
                            frierson_Rs     ,frierson_DeltaS,frierson_Tau_eqtr,frierson_Tau_pole,  &
                            frierson_LinFrac,frierson_Boltz ,frierson_C0                           )

    ! Add values for history output
    !---------------------------------
    call addfld('gray_QRL'   ,(/'lev' /),'A','K/s'    ,'Longwave heating rate for gray atmosphere'       )
    call addfld('gray_QRS'   ,(/'lev' /),'A','K/s'    ,'Solar heating rate for gray atmosphere'          )
    call addfld('gray_DTCOND',(/'lev' /),'A','K/s'    ,'T tendency - gray atmosphere moist process'      )
    call addfld('gray_DQCOND',(/'lev' /),'A','kg/kg/s','Q tendency - gray atmosphere moist process'      )
    call addfld('gray_EVAPDT',(/'lev' /),'A','K/s'    ,'T tendency due to re-evaporation'                )
    call addfld('gray_EVAPDQ',(/'lev' /),'A','kg/kg/s','Q tendency due to re-evaporation'                )
    call addfld('gray_KVH'   ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (heat/moisture)')
    call addfld('gray_KVM'   ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (momentum)'     )
    call addfld('gray_VSE'   ,(/'lev' /),'A','K'      ,'VSE: (Tv + gZ/Cp)'                               )
    call addfld('gray_Zm'    ,(/'lev' /),'A','m'      ,'Geopotential height'                             )
    call addfld('gray_Rf'    ,(/'lev' /),'A','1'      ,'Bulk Richardson number (Frierson et al 2006, eq 16) / Ri_c' )
    call addfld('gray_DTV'   ,(/'lev' /),'A','K/s'    ,'T tendency due to vertical diffusion'            )
    call addfld('gray_DUV'   ,(/'lev' /),'A','m/s2'   ,'U tendency due to vertical diffusion'            )
    call addfld('gray_DVV'   ,(/'lev' /),'A','m/s2'   ,'V tendency due to vertical diffusion'            )
    call addfld('gray_VD01'  ,(/'lev' /),'A','kg/kg/s','Q tendency (vertical diffusion)'                 )
    call addfld('gray_PRECL' ,horiz_only,'A','m/s'    ,'Large-scale precipitation rate'                  )
    call addfld('gray_PRECC' ,horiz_only,'A','m/s'    ,'Convective precipitation rate'                   )
    call addfld('gray_Tsurf ',horiz_only,'I','K'      ,'Surface Temperature'                             )
    call addfld('gray_Qsurf ',horiz_only,'I','kg/kg'  ,'Surface Water Vapor'                             )
    call addfld('gray_Cdrag' ,horiz_only,'A','1'      ,'Surface Drag Coefficient'                        )
    call addfld('gray_Zpbl'  ,horiz_only,'I','m'      ,'PBL Height'                                      )
    call addfld('gray_SWflux',horiz_only,'I','W/m2'   ,'SW Solar Flux'                                   )
    call addfld('gray_LUflux',horiz_only,'I','W/m2'   ,'LW Upward Radiative Flux at Surface'             )
    call addfld('gray_LDflux',horiz_only,'I','W/m2'   ,'LW Downward Radiative Flux at Surface'           )
    call addfld('gray_LWflux',horiz_only,'I','W/m2'   ,'LW Net Radiative Flux at Surface'                )
    call addfld('gray_LUflux_TOA',horiz_only,'I','W/m2'   ,'LW Upward Radiative Flux at TOA'             )
    call addfld('gray_LDflux_TOA',horiz_only,'I','W/m2'   ,'LW Downward Radiative Flux at TOA'           )
    call addfld('gray_LWflux_TOA',horiz_only,'I','W/m2'   ,'LW Net Radiative Flux at TOA'                )
    call addfld('gray_SHflux',horiz_only,'I','W/m2'   , 'Sensible Heat Flux'        )
    call addfld('gray_LHflux',horiz_only,'I','W/m2'   , 'Latent Heat Flux'          )
    call addfld('gray_TauU'  ,horiz_only,'I','N/m2'   , 'U Surface Stress'          )
    call addfld('gray_TauV'  ,horiz_only,'I','N/m2'   , 'V Surface Stress'          )

    call add_default('gray_QRL'   ,1,' ')
    call add_default('gray_QRS'   ,1,' ')
    call add_default('gray_DTCOND',1,' ')
    call add_default('gray_DQCOND',1,' ')
    call add_default('gray_EVAPDT',1,' ')
    call add_default('gray_EVAPDQ',1,' ')
    call add_default('gray_KVH'   ,1,' ')
    call add_default('gray_KVM'   ,1,' ')
    call add_default('gray_VSE'   ,1,' ')
    call add_default('gray_Zm'    ,1,' ')
    call add_default('gray_Rf'    ,1,' ')
    call add_default('gray_DTV'   ,1,' ')
    call add_default('gray_DUV'   ,1,' ')
    call add_default('gray_DVV'   ,1,' ')
    call add_default('gray_VD01'  ,1,' ')
    call add_default('gray_PRECC' ,1,' ')
    call add_default('gray_PRECL' ,1,' ')
    call add_default('gray_Tsurf' ,1,' ')
    call add_default('gray_Qsurf' ,1,' ')
    call add_default('gray_Cdrag' ,1,' ')
    call add_default('gray_Zpbl'  ,1,' ')
    call add_default('gray_SWflux',1,' ')
    call add_default('gray_LUflux',1,' ')
    call add_default('gray_LDflux',1,' ')
    call add_default('gray_LWflux',1,' ')
    call add_default('gray_LUflux_TOA',1,' ')
    call add_default('gray_LDflux_TOA',1,' ')
    call add_default('gray_LWflux_TOA',1,' ')
    call add_default('gray_SHflux',1,' ')
    call add_default('gray_LHflux',1,' ')
    call add_default('gray_TauU'  ,1,' ')
    call add_default('gray_TauV'  ,1,' ')

    ! Allocate Global arrays
    !-------------------------
    allocate(Fsolar(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fsolar',pcols*(endchunk-begchunk+1))
    allocate(Fup   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fup'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fdown' ,pcols*(endchunk-begchunk+1))
    allocate(Fup_toa   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fup_toa'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown_toa (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fdown_toa' ,pcols*(endchunk-begchunk+1))
    allocate(Fnet(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fnet',pcols*(endchunk-begchunk+1))
    allocate(Fnet_toa(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fnet_toa'  ,pcols*(endchunk-begchunk+1))

    allocate(SHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','SHflux',pcols*(endchunk-begchunk+1))
    allocate(LHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','LHflux',pcols*(endchunk-begchunk+1))
    allocate(TUflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TUflux',pcols*(endchunk-begchunk+1))
    allocate(TVflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TVflux',pcols*(endchunk-begchunk+1))
    allocate(Cd    (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Cd'    ,pcols*(endchunk-begchunk+1))
    allocate(clat  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','clat'  ,pcols*(endchunk-begchunk+1))

    ! Initialize time indices and latitudes
    !----------------------------------------
    do lchnk = begchunk,endchunk
      ncol = get_ncols_p(lchnk)
      do icol = 1,ncol
        clat(icol,lchnk) = get_rlat_p(lchnk,icol)
      end do
    end do

    ! At first model step, initialize some values
    !-----------------------------------------------
    if(is_first_step()) then
      ! Initialize physics buffer values
      !----------------------------------
      call pbuf_set_field(pbuf2d, prec_pcw_idx, 0._r8)
      call pbuf_set_field(pbuf2d, prec_dp_idx , 0._r8)

      ! Allocate Surface fields
      !-------------------------
      allocate(Tsurf (pcols,begchunk:endchunk),stat=istat)
      call alloc_err(istat,'Frierson INIT','Tsurf' ,pcols*(endchunk-begchunk+1))
      allocate(Qsurf (pcols,begchunk:endchunk)  ,stat=istat)
      call alloc_err(istat,'Frierson INIT','Qsurf' ,pcols*(endchunk-begchunk+1))

      ! Initialize Surface temperatures and Q
      !-----------------------------------------------------------------------
      do lchnk = begchunk,endchunk
        ncol = get_ncols_p(lchnk)

        ! Set to reference values for initialization
        !------------------------------------------------------------
        phys_state(lchnk)%ps(:ncol) = ps0

        call frierson_surface_init(ncol,         clat(:ncol,lchnk), &
                                 phys_state(lchnk)%ps(:ncol),       &
                                                Tsurf(:ncol,lchnk), &
                                                Qsurf(:ncol,lchnk)  )
      end do
    endif

    ! Initialize radiation and flux values to 0.0
    !---------------------------------------------------------------------------
    do lchnk = begchunk,endchunk
      Fsolar(:,lchnk) = 0._r8
      Fup   (:,lchnk) = 0._r8
      Fdown (:,lchnk) = 0._r8
      Fup_toa   (:,lchnk) = 0._r8
      Fdown_toa (:,lchnk) = 0._r8
      SHflux(:,lchnk) = 0._r8
      LHflux(:,lchnk) = 0._r8
      TUflux(:,lchnk) = 0._r8
      TVflux(:,lchnk) = 0._r8
      Cd    (:,lchnk) = 0._r8
      Fnet  (:,lchnk) = 0._r8
    end do

    ! Informational Output
    !----------------------
    if(masterproc) then
      write(iulog,*) ' '
      write(iulog,*) '-----------------------------------------------------------'
      write(iulog,*) '  FRIERSON MODULE INITIALIZED WITH THE FOLLOWING SETTINGS: '
      write(iulog,*) '-----------------------------------------------------------'
      write(iulog,*) 'FRIERSON: gravit='    , gravit
      write(iulog,*) 'FRIERSON: cappa='     , cappa
      write(iulog,*) 'FRIERSON: rair ='     , rair
      write(iulog,*) 'FRIERSON: cpair='     , cpair
      write(iulog,*) 'FRIERSON: latvap='    , latvap
      write(iulog,*) 'FRIERSON: rh2o='      , rh2o
      write(iulog,*) 'FRIERSON: epsilo='    , epsilo
      write(iulog,*) 'FRIERSON: rhoh2o='    , rhoh2o
      write(iulog,*) 'FRIERSON: zvir='      , zvir
      write(iulog,*) 'FRIERSON: ps0='       , ps0
      write(iulog,*) 'FRIERSON: etamid='    , etamid
      write(iulog,*) 'FRIERSON: T0='        , frierson_T0
      write(iulog,*) 'FRIERSON: E0='        , frierson_E0
      write(iulog,*) 'FRIERSON: Erad='      , frierson_Erad
      write(iulog,*) 'FRIERSON: Wind_min='  , frierson_Wind_min
      write(iulog,*) 'FRIERSON: Z0='        , frierson_Z0
      write(iulog,*) 'FRIERSON: Ri_c='      , frierson_Ri_c
      write(iulog,*) 'FRIERSON: Karman='    , frierson_Karman
      write(iulog,*) 'FRIERSON: Fb='        , frierson_Fb
      write(iulog,*) 'FRIERSON: Rs0='       , frierson_Rs0
      write(iulog,*) 'FRIERSON: Albedo='    , frierson_Albedo
      write(iulog,*) 'FRIERSON: Rs='        , frierson_Rs
      write(iulog,*) 'FRIERSON: DeltaS='    , frierson_DeltaS
      write(iulog,*) 'FRIERSON: Tau_eqtr='  , frierson_Tau_eqtr
      write(iulog,*) 'FRIERSON: Tau_pole='  , frierson_Tau_pole
      write(iulog,*) 'FRIERSON: LinFrac='   , frierson_LinFrac
      write(iulog,*) 'FRIERSON: Boltz='     , frierson_Boltz
      write(iulog,*) 'FRIERSON: C0='        , frierson_C0
      write(iulog,*) 'FRIERSON: Tmin='      , frierson_Tmin
      write(iulog,*) 'FRIERSON: Tdlt='      , frierson_Tdlt
      write(iulog,*) 'FRIERSON: Twidth='    , frierson_Twidth
      write(iulog,*) 'FRIERSON: WetDryCoef=', frierson_WetDryCoef
      write(iulog,*) ' '
    endif

  end subroutine frierson_init
  !==============================================================================


  !==============================================================================


  !==============================================================================
  subroutine frierson_condensate_tend(state, ptend, ztodt, pbuf)
    !
    ! frierson_condensate_tend: Run the selected process to compute precipitation
    !                           due to large scale condensation.
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use frierson,     only: frierson_condensate_NONE,frierson_condensate
    use frierson,     only: frierson_condensate_USER,frierson_condensate_TJ16
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,intent(inout):: state
    real(r8)                 ,intent(in)   :: ztodt
    type(physics_ptend)      ,intent(out)  :: ptend
    type(physics_buffer_desc),pointer      :: pbuf(:)
    !
    ! Local Values
    !-----------------
    real(r8),pointer:: relhum  (:,:)
    real(r8),pointer:: prec_pcw(:)          ! large scale precip
    real(r8)        :: prec_cnv(state%ncol) ! Convective Precip
    real(r8)        :: evapdt(state%ncol, pver) ! T tendency due to re-evaporation of condensation
    real(r8)        :: evapdq(state%ncol, pver) ! Q tendency due to re-evaporation of condensation
    real(r8)        :: dtcond(state%ncol, pver) ! Temperature tendency due to condensation
    real(r8)        :: dqcond(state%ncol, pver) ! Q tendency due to condensation
    real(r8)        :: T     (state%ncol, pver) ! T temporary
    real(r8)        :: qv    (state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)                ! Calc tendencies?
    integer         :: lchnk                    ! chunk identifier
    integer         :: ncol                     ! number of atmospheric columns
    integer         :: k

    ! Set local copies of values
    !---------------------------------
    lchnk       = state%lchnk
    ncol        = state%ncol
    T (:ncol,:) = state%T(:ncol,:)
    qv(:ncol,:) = state%Q(:ncol,:,1)

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Frierson condensate', &
                                ls=.true., lu=.true., lv=.true., lq=lq)

    ! Get values from the physics buffer
    !------------------------------------
    call pbuf_get_field(pbuf,prec_pcw_idx,prec_pcw)
    call pbuf_get_field(pbuf,  relhum_idx,relhum  )

    ! Initialize values for condensate tendencies
    !---------------------------------------------
    do k = 1, pver
      dtcond(:ncol,k) = state%T(:ncol,k)
      dqcond(:ncol,k) = state%q(:ncol,k,1)
    end do

    ! Call the Selected condensation routine  ~~DEVO style~~
    !--------------------------------------------------------
    if(CONDENSATE_OPT == CONDENSATE_NONE) then
      prec_cnv(:ncol)   = 0._r8
      evapdt  (:ncol,:) = 0._r8
      evapdq  (:ncol,:) = 0._r8
      call frierson_condensate_NONE(ncol,pver,state%pmid(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT == CONDENSATE_FRIERSON) then
      prec_cnv(:ncol)   = 0._r8
      call frierson_condensate(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                               state%pdel(:ncol,:), &
                                                        T(:ncol,:), &
                                                       qv(:ncol,:), &
                                                   relhum(:ncol,:), &
                                                 prec_pcw(:ncol)  , &
                                                   evapdt(:ncol,:), &
                                                   evapdq(:ncol,:)  )
    elseif(CONDENSATE_OPT == CONDENSATE_TJ16) then
      prec_cnv(:ncol)   = 0._r8
      evapdt  (:ncol,:) = 0._r8
      evapdq  (:ncol,:) = 0._r8
      call frierson_condensate_TJ16(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT == CONDENSATE_USER) then
      prec_cnv(:ncol)   = 0._r8
      evapdt  (:ncol,:) = 0._r8
      evapdq  (:ncol,:) = 0._r8
      call frierson_condensate_USER(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    else
      ! ERROR: Unknown CONDENSATE_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown CONDENSATE_OPT=',CONDENSATE_OPT
      call endrun('frierson_condensate_tend() CONDENSATE_OPT ERROR')
    endif

    ! Back out temperature and specific humidity
    ! tendencies from updated fields
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k)   = (T (:,k)-state%T(:ncol,k)  )/ztodt*cpair
      ptend%q(:ncol,k,1) = (qv(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! Output condensate tendencies
    !------------------------------
    do k = 1, pver
      dtcond(:ncol,k) = (T (:ncol,k) - dtcond(:ncol,k))/ztodt
      dqcond(:ncol,k) = (qv(:ncol,k) - dqcond(:ncol,k))/ztodt
    end do
    call outfld('gray_EVAPDT',evapdt  ,ncol,lchnk)
    call outfld('gray_EVAPDQ',evapdq  ,ncol,lchnk)
    call outfld('gray_DTCOND',dtcond  ,ncol,lchnk)
    call outfld('gray_DQCOND',dqcond  ,ncol,lchnk)
    call outfld('gray_PRECL' ,prec_pcw,ncol,lchnk)
    call outfld('gray_PRECC' ,prec_cnv,ncol,lchnk)

  end subroutine frierson_condensate_tend
  !==============================================================================


  !============================================================================
  subroutine frierson_pbl_tend(state, ptend, ztodt, cam_in)
    !
    ! frierson_pbl_tend: Run the selected PBL process.
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use phys_grid,    only: get_rlat_all_p
    use frierson,     only: frierson_pbl,frierson_pbl_USER
    !
    ! Passed Variables
    !-------------------
    type(physics_state),intent(in)   :: state
    real(r8),           intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    !
    ! Local Values
    !----------------
    real(r8) :: T         (state%ncol,pver)   ! T temporary
    real(r8) :: qv        (state%ncol,pver)   ! Q temporary (specific humidity)
    real(r8) :: U         (state%ncol,pver)   ! U temporary
    real(r8) :: V         (state%ncol,pver)   ! V temporary
    real(r8) :: dqdt_vdiff(state%ncol,pver)   ! PBL Q vertical diffusion tend kg/kg/s
    real(r8) :: dtdt_vdiff(state%ncol,pver)   ! PBL T vertical diffusion tend  K/s
    real(r8) :: dudt_vdiff(state%ncol,pver)   ! PBL U vertical diffusion tend  m/s/s
    real(r8) :: dvdt_vdiff(state%ncol,pver)   ! PBL V vertical diffusion tend  m/s/s
    real(r8) :: Km        (state%ncol,pverp)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: Ke        (state%ncol,pverp)  ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: VSE       (state%ncol,pver)   ! Dry Static Energy divided by Cp (K)
    real(r8) :: Zm        (state%ncol,pver)   !
    real(r8) :: Zi        (state%ncol,pver)   !
    real(r8) :: Z_pbl     (state%ncol)        !
    real(r8) :: Rf        (state%ncol,pver)   !
    real(r8) :: Tsfc      (state%ncol)        ! Surface T
    real(r8) :: Qsfc      (state%ncol)        ! Surface Q (saturated)
    real(r8) :: Cdrag     (state%ncol)        ! Cdrag coef from surface calculation

    logical  :: lq        (pcnst)             ! Calc tendencies?
    real(r8) :: dTs       (state%ncol)
    real(r8) :: dUa       (state%ncol,pver)
    real(r8) :: dVa       (state%ncol,pver)
    real(r8) :: dTa       (state%ncol,pver)
    real(r8) :: dQa       (state%ncol,pver)
    integer  :: lchnk                        ! chunk identifier
    integer  :: ncol                         ! number of atmospheric columns
    integer  :: kk                           ! loop index

    ! Set local copies of values
    !---------------------------------
    lchnk              = state%lchnk
    ncol               = state%ncol
    Zm  (:ncol,:)      = state%zm  (:ncol,:)
    Zi  (:ncol,1:pver) = state%zi  (:ncol,1:pver)
    T   (:ncol,:)      = state%T   (:ncol,:)
    U   (:ncol,:)      = state%U   (:ncol,:)
    V   (:ncol,:)      = state%V   (:ncol,:)
    qv  (:ncol,:)      = state%Q   (:ncol,:,1)

    ! Initialize individual parameterization tendencies
    !-----------------------------------------------------
    lq    = .false.
    lq(1) = .true.
    call physics_ptend_init(ptend,state%psetcols,'Frierson pbl_tend',        &
                                       ls=.true., lu=.true., lv=.true., lq=lq)

    ! Call the Selected PBL routine
    !--------------------------------------------------------
    Tsfc(:ncol) = Tsurf(:ncol,lchnk)
    Qsfc(:ncol) = Qsurf(:ncol,lchnk)
    if(PBL_OPT == PBL_FRIERSON) then
      ! Call Frierson PBL scheme
      !--------------------------------------------------
      call frierson_pbl(ncol, pver, ztodt,state%pmid (:ncol,:),     &
                                          state%pint (:ncol,:),     &
                                                   Zm(:ncol,:),     &
                                                   Zi(:ncol,:),     &
                                             state%ps(:ncol)  ,     &
                                                 Tsfc(:ncol)  ,     &
                                                 Qsfc(:ncol)  ,     &
                                                    T(:ncol,:),     &
                                                    U(:ncol,:),     &
                                                    V(:ncol,:),     &
                                                   qv(:ncol,:),     &
                                               Fsolar(:ncol,lchnk), &
                                                Fdown(:ncol,lchnk), &
                                                Cdrag(:ncol)  ,     &
                                                   Km(:ncol,:),     &
                                                   Ke(:ncol,:),     &
                                                  VSE(:ncol,:),     &
                                                Z_pbl(:ncol)  ,     &
                                                   Rf(:ncol,:),     &
                                           dqdt_vdiff(:ncol,:),     &
                                           dtdt_vdiff(:ncol,:),     &
                                           dudt_vdiff(:ncol,:),     &
                                           dvdt_vdiff(:ncol,:),     &
                                               LHflux(:ncol,lchnk), &
                                               SHflux(:ncol,lchnk), &
                                               TUflux(:ncol,lchnk), &
                                               TVflux(:ncol,lchnk)  )
    elseif(PBL_OPT == PBL_USER) then
      ! Call USER implemented routine in frierson module
      !--------------------------------------------------
      call frierson_pbl_USER(ncol, pver, ztodt,state%pmid (:ncol,:),     &
                                               state%pint (:ncol,:),     &
                                                        Zm(:ncol,:),     &
                                                        Zi(:ncol,:),     &
                                                  state%ps(:ncol)  ,     &
                                                      Tsfc(:ncol)  ,     &
                                                      Qsfc(:ncol)  ,     &
                                                         T(:ncol,:),     &
                                                         U(:ncol,:),     &
                                                         V(:ncol,:),     &
                                                        qv(:ncol,:),     &
                                                    Fsolar(:ncol,lchnk), &
                                                     Fdown(:ncol,lchnk), &
                                                     Cdrag(:ncol)  ,     &
                                                        Km(:ncol,:),     &
                                                        Ke(:ncol,:),     &
                                                       VSE(:ncol,:),     &
                                                     Z_pbl(:ncol)  ,     &
                                                        Rf(:ncol,:),     &
                                                dqdt_vdiff(:ncol,:),     &
                                                dtdt_vdiff(:ncol,:),     &
                                                dudt_vdiff(:ncol,:),     &
                                                dvdt_vdiff(:ncol,:),     &
                                                    LHflux(:ncol,lchnk), &
                                                    SHflux(:ncol,lchnk), &
                                                    TUflux(:ncol,lchnk), &
                                                    TVflux(:ncol,lchnk)  )
    else
      ! ERROR: Unknown PBL_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown PBL_OPT=',PBL_OPT
      call endrun('frierson_pbl_tend() PBL_OPT ERROR')
    endif
    Tsurf(:ncol,lchnk) = Tsfc (:ncol)
    Qsurf(:ncol,lchnk) = Qsfc (:ncol)
    Cd   (:ncol,lchnk) = Cdrag(:ncol)

    ! Back out tendencies from updated fields
    !-----------------------------------------
    do kk = 1, pver
      ptend%s(:ncol,kk  ) = (T (:,kk)-state%T(:ncol,kk  ))/ztodt*cpair
      ptend%u(:ncol,kk  ) = (U (:,kk)-state%U(:ncol,kk  ))/ztodt
      ptend%v(:ncol,kk  ) = (V (:,kk)-state%V(:ncol,kk  ))/ztodt
      ptend%q(:ncol,kk,1) = (qv(:,kk)-state%q(:ncol,kk,1))/ztodt
    end do

    ! Archive diagnostic fields
    !----------------------------
    call outfld('gray_Tsurf' ,Tsurf(:ncol,lchnk) ,ncol,lchnk)
    call outfld('gray_Qsurf' ,Qsurf(:ncol,lchnk) ,ncol,lchnk)
    call outfld('gray_Cdrag' ,Cd   (:ncol,lchnk) ,ncol,lchnk)
    call outfld('gray_Zpbl'  ,Z_pbl              ,ncol,lchnk) !
    call outfld('gray_KVH'   ,Ke                 ,ncol,lchnk) ! Eddy diffusivity (heat and moisture,m2/s)
    call outfld('gray_KVM'   ,Km                 ,ncol,lchnk) ! Eddy diffusivity (momentum, m2/s)
    call outfld('gray_VSE'   ,VSE                ,ncol,lchnk) ! Virtual Dry Static Energy divided by Cp (K)
    call outfld('gray_Zm'    ,Zm                 ,ncol,lchnk) !
    call outfld('gray_Rf'    ,Rf                 ,ncol,lchnk) !
    call outfld('gray_DTV'   ,dtdt_vdiff         ,ncol,lchnk) ! PBL + surface flux T tendency (K/s)
    call outfld('gray_DUV'   ,dudt_vdiff         ,ncol,lchnk) ! PBL u tendency (m/s2)
    call outfld('gray_DVV'   ,dvdt_vdiff         ,ncol,lchnk) ! PBL v tendency (m/s2)
    call outfld('gray_VD01'  ,dqdt_vdiff         ,ncol,lchnk) ! PBL + surface flux Q tendency (kg/kg/s)
    call outfld('gray_SHflux',SHflux(:ncol,lchnk),ncol,lchnk) ! Sensible Heat Flux
    call outfld('gray_LHflux',LHflux(:ncol,lchnk),ncol,lchnk) ! Latent Heat Flux
    call outfld('gray_TauU'  ,TUflux(:ncol,lchnk),ncol,lchnk) ! U Surface Stress
    call outfld('gray_TauV'  ,TVflux(:ncol,lchnk),ncol,lchnk) ! V Surface Stress

  end subroutine frierson_pbl_tend
  !============================================================================


  !============================================================================
  subroutine frierson_radiative_tend(state, ptend, ztodt,cam_in,cam_out)
    !
    ! frierson_radiative_tend: Run the radiative process
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use phys_grid,    only: get_rlat_all_p
    use frierson,     only: frierson_radiation,frierson_radiation_USER
    !
    ! Passed Variables
    !------------------
    type(physics_state),intent(in)   :: state
    real(r8)           ,intent(in)   :: ztodt
    type(physics_ptend),intent(out)  :: ptend
    type(cam_in_t),     intent(inout):: cam_in
    type(cam_out_t),    intent(inout):: cam_out
    !
    ! Local Values
    !---------------
    real(r8):: T           (state%ncol,pver) ! T temporary
    real(r8):: qv          (state%ncol,pver) ! Q temporary
    real(r8):: dtdt_heating(state%ncol,pver) ! Longwave heating tendency K/s
    real(r8):: dtdt_solar  (state%ncol,pver) ! Shortwave heating tendency K/s
    real(r8):: Tsfc        (state%ncol)      ! Surface T
    real(r8):: Qsfc        (state%ncol)      ! Surface Q (saturated)
    logical :: lq(pcnst)                     ! Calc tendencies?
    integer :: lchnk                         ! chunk identifier
    integer :: ncol                          ! number of atmospheric columns
    integer :: k                             ! loop index

    ! Copy to local values
    !-------------------------------------------------
    lchnk         = state%lchnk
    ncol          = state%ncol
    T   (:ncol,:) = state%T(:ncol,:)
    qv  (:ncol,:) = state%Q(:ncol,:,1)

    !--------------------------------------
    Tsfc(:ncol)   = Tsurf(:ncol,lchnk)
    Qsfc(:ncol)   = Qsurf(:ncol,lchnk)

    ! initialize individual parameterization tendencies
    !---------------------------------------------------
    lq(:) = .false.
    call physics_ptend_init(ptend, state%psetcols, 'Frierson radiative_tend',   &
                                        ls=.true., lu=.false., lv=.false., lq=lq)

    ! Call the Selected radiative routine
    !--------------------------------------------------------
    if(RADIATION_OPT == RADIATION_FRIERSON) then
      call frierson_radiation(ncol,pver,ztodt,clat(:ncol,lchnk), &
                                        state%pint(:ncol,:),     &
                                        state%pmid(:ncol,:),     &
                                          state%ps(:ncol),       &
                                              Tsfc(:ncol),       &
                                              Qsfc(:ncol),       &
                                                 T(:ncol,:),     &
                                                qv(:ncol,:),     &
                                      dtdt_heating(:ncol,:),     &
                                            Fsolar(:ncol,lchnk), &
                                               Fup(:ncol,lchnk), &
                                             Fdown(:ncol,lchnk), &
                                           Fup_toa(:ncol,lchnk), &
                                         Fdown_toa(:ncol,lchnk)  )
      dtdt_solar(:ncol,:) = 0._r8
    elseif(RADIATION_OPT == RADIATION_USER) then
      call frierson_radiation_USER(ncol,pver,ztodt,clat(:ncol,lchnk), &
                                             state%pint(:ncol,:),     &
                                             state%pmid(:ncol,:),     &
                                               state%ps(:ncol),       &
                                                   Tsfc(:ncol),       &
                                                   Qsfc(:ncol),       &
                                                      T(:ncol,:),     &
                                                     qv(:ncol,:),     &
                                           dtdt_heating(:ncol,:),     &
                                                 Fsolar(:ncol,lchnk), &
                                                    Fup(:ncol,lchnk), &
                                                  Fdown(:ncol,lchnk), &
                                                Fup_toa(:ncol,lchnk), &
                                              Fdown_toa(:ncol,lchnk)  )
      dtdt_solar(:ncol,:) = 0._r8
    else
      ! ERROR: Unknown RADIATION_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown RADIATION_OPT=',RADIATION_OPT
      call endrun('frierson_pbl_tend() RADIATION_OPT ERROR')
    endif

    Fnet  (:ncol,lchnk)     = Fup(:ncol,lchnk)     - Fdown (:ncol,lchnk)
    Fnet_toa  (:ncol,lchnk) = Fup_toa(:ncol,lchnk) - Fdown_toa (:ncol,lchnk)

    ! Copy downward LW radiative heating values to cam_out%
    !---------------------------------------------------------
    cam_out%flwds(:ncol) = Fdown (:ncol,lchnk)
    cam_out%netsw(:ncol) = Fsolar(:ncol,lchnk)
    cam_out%sols (:ncol) = Fsolar(:ncol,lchnk)
    cam_out%solsd(:ncol) = Fsolar(:ncol,lchnk)
    cam_out%soll (:ncol) = Fsolar(:ncol,lchnk)
    cam_out%solld(:ncol) = Fsolar(:ncol,lchnk)

    ! Back out tendencies from updated T field
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k) = (T(:,k)-state%T(:ncol,k))/ztodt*cpair
    end do

    ! Archive T tendency from temperature relaxation (mimics radiation, K/s)
    !-----------------------------------------------------------------------
    call outfld('gray_QRL'   ,dtdt_heating, ncol,lchnk)
    call outfld('gray_QRS'   ,dtdt_solar  , ncol,lchnk)
    call outfld('gray_SWflux',Fsolar(:ncol,lchnk)      , ncol,lchnk)
    call outfld('gray_LUflux',Fup(:ncol,lchnk)         , ncol,lchnk)
    call outfld('gray_LDflux',Fdown(:ncol,lchnk)       , ncol,lchnk)
    call outfld('gray_LWflux',Fnet(:ncol,lchnk)        , ncol,lchnk)
    call outfld('gray_LUflux_TOA',Fup_toa(:ncol,lchnk)     , ncol,lchnk)
    call outfld('gray_LDflux_TOA',Fdown_toa(:ncol,lchnk)   , ncol,lchnk)
    call outfld('gray_LWflux_TOA',Fnet_toa(:ncol,lchnk)    , ncol,lchnk)

  end subroutine frierson_radiative_tend
  !============================================================================


  !=======================================================================
  subroutine frierson_surface_init(ncol, clat, PS, Tsfc, Qsfc)
    !
    !
    !==========================================================================
    !
    ! Passed variables
    !--------------------
    integer ,intent(in) :: ncol
    real(r8),intent(in) :: clat (ncol)
    real(r8),intent(in) :: PS   (ncol)
    real(r8),intent(out):: Tsfc(ncol)
    real(r8),intent(out):: Qsfc(ncol)
    !
    ! Local values
    !--------------
    integer :: ii
    real(r8):: T_width

    ! set SST profile
    !------------------
    T_width = frierson_Twidth*pi/180.0_r8
    do ii = 1, ncol
      Tsfc(ii) = frierson_Tmin + frierson_Tdlt*exp(-((clat(ii)/T_width)**2)/2.0_r8)
      Qsfc(ii) = epsilo*frierson_E0/PS(ii)                             &
                 *exp(-latvap/rh2o*((1._r8/Tsfc(ii))-1._r8/frierson_T0))
    end do

  end subroutine frierson_surface_init
  !=======================================================================


  !=======================================================================
  subroutine frierson_restart_init(File,hdimids,hdimcnt)
    !
    ! frierson_restart_init:
    !==========================================================================
    !
    ! Passed variables
    !--------------------
    type(file_desc_t),intent(inout):: File
    integer          ,intent(in)   :: hdimcnt
    integer          ,intent(in)   :: hdimids(1:hdimcnt)
    !
    ! Local values
    !--------------
    integer:: ierr

    ierr = pio_def_var(File,'Frierson_Tsfc',pio_double, hdimids, Tsurf_desc)
    if (ierr /= 0) then
       call endrun('frierson_restart_init: ERROR defining Frierson_Tsfc')
    end if

    ierr = pio_def_var(File,'Frierson_Qsfc',pio_double, hdimids, Qsurf_desc)
    if (ierr /= 0) then
       call endrun('frierson_restart_init: ERROR defining Frierson_Qsfc')
    end if

  end subroutine frierson_restart_init
  !=======================================================================


  !=======================================================================
  subroutine frierson_restart_write(File)
    !
    ! frierson_restart_write:
    !==========================================================================
    !
    ! Passed variables
    !--------------------
    type(file_desc_t),intent(inout):: File
    !
    ! Local values
    !--------------
    type(io_desc_t),pointer:: iodesc
    integer:: dims(3),gdims(3),nhdims
    integer:: physgrid
    integer:: ierr

    ! Get the iodesc for write calls
    !---------------------------------
    dims(1) = pcols
    dims(2) = endchunk - begchunk + 1
    physgrid = cam_grid_id('physgrid')
    call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
    call cam_grid_get_decomp(physgrid,  dims(1:2), gdims(1:nhdims), pio_double, iodesc)

    ! Write Surface values
    !---------------------
    call pio_write_darray(File, Tsurf_desc, iodesc, Tsurf, ierr)
    if (ierr /= 0) then
       call endrun('frierson_restart_write: ERROR writing Tsurf')
    end if

    call pio_write_darray(File, Qsurf_desc, iodesc, Qsurf, ierr)
    if (ierr /= 0) then
       call endrun('frierson_restart_write: ERROR writing Qsurf')
    end if

  end subroutine frierson_restart_write
  !=======================================================================


  !=======================================================================
  subroutine frierson_restart_read(File)
    !
    ! frierson_restart_read:
    !==========================================================================
    use error_messages,only: alloc_err
    !
    ! Passed variables
    !--------------------
    type(file_desc_t),intent(inout):: File
    !
    ! Local values
    !--------------
    type( io_desc_t),pointer:: iodesc
    type(var_desc_t)        :: vardesc
    integer:: dims(3),gdims(3),nhdims
    integer:: physgrid
    integer:: ierr

    ! Allocate space for the restart fields
    !-----------------------------------------
    allocate(Tsurf (pcols,begchunk:endchunk),stat=ierr)
    call alloc_err(ierr,'Frierson RESTART','Tsurf' ,pcols*(endchunk-begchunk+1))
    allocate(Qsurf (pcols,begchunk:endchunk)  ,stat=ierr)
    call alloc_err(ierr,'Frierson RESTART','Qsurf' ,pcols*(endchunk-begchunk+1))

    ! Get the iodesc for read calls
    !---------------------------------
    dims(1) = pcols
    dims(2) = endchunk - begchunk + 1
    physgrid = cam_grid_id('physgrid')
    call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)
    call cam_grid_get_decomp(physgrid,  dims(1:2), gdims(1:nhdims), pio_double, iodesc)

    ! Read Surface values
    !---------------------
    ierr = pio_inq_varid(File,'Frierson_Tsfc',vardesc)
    if (ierr /= 0) then
       call endrun('frierson_restart_read: ERROR PIO unable to find variable Frierson_Tsfc')
    end if

    call pio_read_darray(File, vardesc, iodesc, Tsurf, ierr)
    if (ierr /= 0) then
       call endrun('frierson_restart_read: ERROR PIO unable to read variable Tsurf')
    end if

    ierr = pio_inq_varid(File,'Frierson_Qsfc',vardesc)
    if (ierr /= 0) then
       call endrun('frierson_restart_read: ERROR PIO unable to find variable Frierson_Qsfc')
    end if

    call pio_read_darray(File, vardesc, iodesc, Qsurf, ierr)
    if (ierr /= 0) then
       call endrun('frierson_restart_read: ERROR PIO unable to read variable Qsurf')
    end if

  end subroutine frierson_restart_read
  !=======================================================================

end module frierson_cam

