module frierson_cam
!-----------------------------------------------------------------------
!
! Purpose: Implement idealized forcings described in 
!          Frierson, et al. (2006), " A Gray-Radiation Aquaplanet
!          Moist GCM, Part I. Static Stability and Eddy Scale"
!          J. Atmos. Sci, Vol 63, 2548-2566.
!
!============================================================================
  ! Useful modules
  !-------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: pi => shr_const_pi
  use physconst,      only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir
  use ppgrid,         only: pcols, pver, begchunk, endchunk
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
  use spmd_utils,     only: mpicom, mstrid=>masterprocid, mpi_logical, mpi_real8

  use pio             ,only: file_desc_t, var_desc_t, io_desc_t, pio_double, pio_def_var
  use pio             ,only: pio_write_darray, pio_read_darray, pio_inq_varid
  use cam_grid_support,only: cam_grid_id, cam_grid_dimensions, cam_grid_get_decomp
  use shr_const_mod,   only: SHR_CONST_STEBOL, SHR_CONST_REARTH, SHR_CONST_KARMAN

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure. 
  !---------------------------------------------------------
  implicit none
  private
  save

  public :: frierson_register
  public :: frierson_readnl
  public :: frierson_init
  public :: frierson_convection_tend
  public :: frierson_condensate_tend
  public :: frierson_pbl_tend
  public :: frierson_radiative_tend
  private:: frierson_surface_init
  public :: frierson_restart_init
  public :: frierson_restart_write
  public :: frierson_restart_read

  ! PBL Configuatons
  !------------------
  integer,parameter:: PBL_FRIERSON = 0           ! Implementation of Frierson PBL
  integer,parameter:: PBL_USER     = 1           ! Optional call for user defined PBL

  ! Tags to identify optional model formulations
  !------------------------------------------------
  integer,parameter:: CONVECTION_NONE       = 0  ! No convection, PRECC=0
  integer,parameter:: CONVECTION_FRIERSON   = 1  ! PFC todo: **Betts-Miller convection here..
  integer,parameter:: CONVECTION_USER       = 2  ! Optional user defined Convection scheme

  integer,parameter:: CONDENSATE_NONE       = 0  ! No Condensation, PRECL=0
  integer,parameter:: CONDENSATE_FRIERSON   = 1  ! Frierson condensation w/ re-evaporation
  integer,parameter:: CONDENSATE_TJ16       = 2  ! Consensation from TJ2016 model.
  integer,parameter:: CONDENSATE_USER       = 3  ! Optional user defined Consensation scheme

  integer,parameter:: RADIATION_FRIERSON    = 0  ! Frierson Gray radiation.
  integer,parameter:: RADIATION_USER        = 1  ! Optional user defined Radiation scheme

  ! Options selecting which PRECIP, PBL, RADIATION, etc.. formulations to use.
  !   PFC todo: Eventually change these to namelist options??
  !---------------------------------------------------------------------------------
  integer,parameter:: PBL_OPT          = PBL_FRIERSON
  integer,parameter:: CONVECTION_OPT   = CONVECTION_NONE
  integer,parameter:: CONDENSATE_OPT   = CONDENSATE_FRIERSON
  integer,parameter:: RADIATION_OPT    = RADIATION_FRIERSON

  ! Global Constants 
  !---------------------
  real(r8),parameter:: frierson_T0     = 273.16_r8     ! Reference Temperature for E0 
  real(r8),parameter:: frierson_E0     = 610.78_r8     ! Saturation Vapor pressure @ T0
  real(r8),parameter:: frierson_Rs0    = 1360.0_r8     ! Solar Constant  
  real(r8),parameter:: frierson_Erad   = SHR_CONST_REARTH  ! Earth Radius
  real(r8),parameter:: frierson_Karman = SHR_CONST_KARMAN  ! Von Karman constant
  real(r8),parameter:: frierson_Boltz  = SHR_CONST_STEBOL  ! Stefan-Boltzmann constant

  ! Some Physics buffer indicies
  !-------------------------------
  integer:: prec_pcw_idx = 0
  integer:: prec_dp_idx  = 0
  integer:: relhum_idx   = 0

  ! Global values for Surface Temp, surface fluxes, and radiative heating
  !----------------------------------------------------------------------
  type(var_desc_t)    :: Tsurf_desc      ! Vardesc for restarts
  type(var_desc_t)    :: Qsurf_desc      ! Vardesc for restarts
  real(r8),allocatable:: Tsurf (:,:)     ! Surface Temp 
  real(r8),allocatable:: Qsurf (:,:)     ! Surface Q
  real(r8),allocatable:: Fsolar(:,:)     ! Net Solar Heating
  real(r8),allocatable:: Fup   (:,:)     ! Upward Longwave heating
  real(r8),allocatable:: Fdown (:,:)     ! Downward Longwave heating
  real(r8),allocatable:: SHflux(:,:)     ! Sensible Heat flux
  real(r8),allocatable:: LHflux(:,:)     ! Latent Heat Flux
  real(r8),allocatable:: SWflux(:,:)     ! Surface Water flux
  real(r8),allocatable:: TUflux(:,:)     ! U momentum flux
  real(r8),allocatable:: TVflux(:,:)     ! V momentum flux
  real(r8),allocatable:: Evap  (:,:)     ! U momentum flux
  real(r8),allocatable:: Cd    (:,:)     ! V momentum flux
  real(r8),allocatable:: clat  (:,:)     ! latitudes(radians) for columns
  real(r8),allocatable:: Fnet  (:,:)     ! Net Radiative Surface Heating

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Global Tuning values
  !------------------------
  real(r8):: frierson_Wind_min   = unset_r8      ! Minimum wind threshold
  real(r8):: frierson_Z0         = unset_r8      ! Roughness Length
  real(r8):: frierson_Ri_c       = unset_r8      ! Crit. Richardson # for stable mixing
  real(r8):: frierson_Fb         = unset_r8      ! Surface layer Fraction
  real(r8):: frierson_Albedo     = unset_r8      ! Frierson Albeo
  real(r8):: frierson_DeltaS     = unset_r8      ! Lat variation of shortwave radiation
  real(r8):: frierson_Tau_eqtr   = unset_r8      ! Longwave optical depth at Equator
  real(r8):: frierson_Tau_pole   = unset_r8      ! Longwave optical depth at poles.
  real(r8):: frierson_LinFrac    = unset_r8      ! Stratosphere Linear optical depth param
  real(r8):: frierson_C0         = unset_r8      ! Ocean mixed layer heat capacity
  real(r8):: frierson_WetDryCoef = unset_r8      ! E0 Scale factor to control moisture
  real(r8):: frierson_Tmin       = unset_r8      ! IC: Minimum sst (K)
  real(r8):: frierson_Tdlt       = unset_r8      ! IC: eq-polar difference sst (K)
  real(r8):: frierson_Twidth     = unset_r8      ! IC: width parameter for sst (C)
  real(r8):: frierson_C_ocn      = unset_r8      ! CACQUESTION -- NEEDS description

contains
  !==============================================================================
  subroutine frierson_register()
    ! 
    ! frierson_register: Register physics buffer values
    !=====================================================================

    call pbuf_add_field('PREC_PCW','physpkg',dtype_r8, (/pcols/),     prec_pcw_idx)
    call pbuf_add_field('PREC_DP' ,'physpkg',dtype_r8, (/pcols/),     prec_dp_idx )
    call pbuf_add_field('RELHUM'  ,'physpkg',dtype_r8, (/pcols,pver/),relhum_idx  )

    ! End Routine
    !-------------
    return
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
                           frierson_Tdlt    , frierson_Twidth    , frierson_C_ocn

    ! Read in namelist values
    !-------------------------
    if(masterproc) then
      unitn = getunit()
      open(unitn,file=trim(nlfile),status='old')
      call find_group_name(unitn,'frierson_nl',status=ierr)
      if(ierr.eq.0) then
        read(unitn,frierson_nl,iostat=ierr)
        if(ierr.ne.0) then
          call endrun('frierson_readnl:: ERROR reading namelist')
        endif
      endif
      close(unitn)
      call freeunit(unitn)
    endif

    ! CACQUESTION - Either add checks or remove this comment
    ! Sanity Check namelist values
    !--------------------------------

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
    call mpi_bcast(frierson_C_ocn, 1, mpi_real8 , mstrid, mpicom, ierr)
    if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: frierson_C_ocn")

    ! End Routine
    !-------------
    return
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
    integer :: istat,lchnk,icol,ncol,ll
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
    call addfld('QRS',(/'lev'/),'A','K/s','Temperature tendency associated with the '//           &
                                          'relaxation toward the equilibrium temperature profile' )
    call addfld('KVH' ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (heat/moisture)')
    call addfld('KVM' ,(/'ilev'/),'A','m2/s'   ,'Vertical diffusion diffusivities (momentum)'     )
    call addfld('VSE' ,(/'lev' /),'A','K'      ,'VSE: (Tv + gZ/Cp)'                               )
    call addfld('Zm'  ,(/'lev' /),'A','m'      ,'ATM Layer Heights use in PBL'                    )
    call addfld('DTV' ,(/'lev' /),'A','K/s'    ,'T vertical diffusion'                            )
    call addfld('DUV' ,(/'lev' /),'A','m/s2'   ,'U vertical diffusion'                            )
    call addfld('DVV' ,(/'lev' /),'A','m/s2'   ,'V vertical diffusion'                            )
    call addfld('VD01',(/'lev' /),'A','kg/kg/s','Q tendency (vertical diffusion)'                 )
    call addfld('Cdrag',horiz_only,'A','n/a'    ,'Surface Drag'                                   )
    call addfld('Z_pbl',horiz_only,'I','m'      ,'PBL Height'                                     )
    call addfld('Rf'  ,(/'lev' /),'I','n/a'     ,'Another Richardson number / Ri_c'               )

    call addfld('R_Fsolar', horiz_only, 'I','W/m2', 'SW Solar Flux'             )
    call addfld('R_Fup'   , horiz_only, 'I','W/m2', 'LW Upward Radiative Flux'  )
    call addfld('R_Fdown' , horiz_only, 'I','W/m2', 'LW Downward Radiative Flux')
    call addfld('R_SHflux', horiz_only, 'I','W/m2', 'Sensible Heat Flux'        )
    call addfld('R_LHflux', horiz_only, 'I','W/m2', 'Latent Heat Flux'          )
    call addfld('R_Fnet  ', horiz_only, 'I','W/m2', 'Net Radiative Flux'        )
    call addfld('R_Tsurf ', horiz_only, 'I','K'   , 'Surface Temperature'       )
    call addfld('R_Qsurf ', horiz_only, 'I','kg/kg', 'Surface Water Vapor'      )
    call addfld('R_Cdrag' , horiz_only, 'I','n/a'  , 'Surface Drag'             )

    call add_default('QRS'  ,1,' ')
    call add_default('KVH'  ,1,' ')
    call add_default('KVM'  ,1,' ')
    call add_default('VSE'  ,1,' ')
    call add_default('Zm'   ,1,' ')
    call add_default('DTV'  ,1,' ')
    call add_default('DUV'  ,1,' ')
    call add_default('DVV'  ,1,' ')
    call add_default('VD01' ,1,' ')
    call add_default('Cdrag',1,' ')
    call add_default('Z_pbl',1,' ')
    call add_default('Rf'   ,1,' ')

    ! Allocate Global arrays
    !-------------------------
    allocate(Fsolar(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fsolar',pcols*(endchunk-begchunk+1))
    allocate(Fup   (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fup'   ,pcols*(endchunk-begchunk+1))
    allocate(Fdown (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fdown' ,pcols*(endchunk-begchunk+1))
    allocate(SHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','SHflux',pcols*(endchunk-begchunk+1))
    allocate(LHflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','LHflux',pcols*(endchunk-begchunk+1))
    allocate(SWflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','SWflux',pcols*(endchunk-begchunk+1))
    allocate(TUflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TUflux',pcols*(endchunk-begchunk+1))
    allocate(TVflux(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','TVflux',pcols*(endchunk-begchunk+1))
    allocate(Evap  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Evap'  ,pcols*(endchunk-begchunk+1))
    allocate(Cd    (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Cd'    ,pcols*(endchunk-begchunk+1))
    allocate(clat  (pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','clat'  ,pcols*(endchunk-begchunk+1))

    allocate(Fnet(pcols,begchunk:endchunk)  ,stat=istat)
    call alloc_err(istat,'Frierson INIT','Fnet',pcols*(endchunk-begchunk+1))

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

        ! phys_state PS values are Inf at this point.... HUH???
        ! Set to reference values for initialization
        !------------------------------------------------------------
        phys_state(lchnk)%ps(:ncol) = ps0
  
        call frierson_surface_init(ncol,         clat(:ncol,lchnk), &
                                 phys_state(lchnk)%ps(:ncol),       &
                                                Tsurf(:ncol,lchnk), &
                                                Qsurf(:ncol,lchnk)  )
      end do
    endif
      
    ! Initialize radition and flux values to 0.0  (Add Init from restart file???)
    !---------------------------------------------------------------------------
    do lchnk = begchunk,endchunk
      Fsolar(:,lchnk) = 0._r8
      Fup   (:,lchnk) = 0._r8
      Fdown (:,lchnk) = 0._r8
      SHflux(:,lchnk) = 0._r8
      LHflux(:,lchnk) = 0._r8
      SWflux(:,lchnk) = 0._r8
      TUflux(:,lchnk) = 0._r8
      TVflux(:,lchnk) = 0._r8
      Evap  (:,lchnk) = 0._r8
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

    ! End Routine
    !--------------
    return
  end subroutine frierson_init
  !==============================================================================


  !==============================================================================
  subroutine frierson_convection_tend(state, ptend, ztodt, pbuf)
    !
    ! frierson_convection_tend: Run the selected process to compute precipitation 
    !                           due to convection.
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
    use frierson,     only: frierson_convection_NONE,frierson_convection
    use frierson,     only: frierson_convection_USER
    !
    ! Passed Variables
    !------------------
    type(physics_state)      ,intent(in) :: state
    real(r8)                 ,intent(in) :: ztodt 
    type(physics_ptend)      ,intent(out):: ptend 
    type(physics_buffer_desc),pointer    :: pbuf(:)
    !
    ! Local Values
    !-----------------
    real(r8),pointer:: prec_dp (:)          ! convective precip
    real(r8),pointer:: relhum  (:,:)
    real(r8)        :: T (state%ncol, pver) ! T temporary
    real(r8)        :: qv(state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)            ! Calc tendencies?
    integer         :: lchnk                ! chunk identifier
    integer         :: ncol                 ! number of atmospheric columns
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
    call physics_ptend_init(ptend, state%psetcols, 'Frierson convection', &
                                ls=.true., lu=.true., lv=.true., lq=lq)

    ! Get values from the physics buffer
    !------------------------------------
    call pbuf_get_field(pbuf,prec_dp_idx ,prec_dp )
    call pbuf_get_field(pbuf,  relhum_idx,relhum  )

    ! Call the Selected convection routine
    !--------------------------------------------------------
    if(CONVECTION_OPT.eq.CONVECTION_NONE) then
      call frierson_convection_NONE(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                       prec_dp(:ncol)    )
    elseif(CONVECTION_OPT.eq.CONVECTION_FRIERSON) then
      call frierson_convection(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                               state%pdel(:ncol,:), &
                                                        T(:ncol,:), &
                                                       qv(:ncol,:), &
                                                   relhum(:ncol,:), &
                                                  prec_dp(:ncol)    )
    elseif(CONVECTION_OPT.eq.CONVECTION_USER) then
      call frierson_convection_USER(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                       prec_dp(:ncol)    )
    else
      ! ERROR: Unknown CONVECTION_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown CONVECTION_OPT=',CONVECTION_OPT
      call endrun('frierson_convection_tend() CONVECTION_OPT ERROR')
    endif

    ! Back out temperature and specific humidity 
    ! tendencies from updated fields
    !--------------------------------------------
    do k = 1, pver
      ptend%s(:ncol,k)   = (T (:,k)-state%T(:ncol,k)  )/ztodt*cpair
      ptend%q(:ncol,k,1) = (qv(:,k)-state%q(:ncol,k,1))/ztodt
    end do

    ! End Routine
    !--------------
    return
  end subroutine frierson_convection_tend
  !==============================================================================


  !==============================================================================
  subroutine frierson_condensate_tend(state, ptend, ztodt, pbuf)
    !
    ! frierson_condensate_tend: Run the selected process to compute precipitation 
    !                           due to large scale condensation.
    !=====================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
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
    real(r8),pointer:: prec_pcw(:)          ! large scale precip
    real(r8),pointer:: relhum  (:,:)
    real(r8)        :: T (state%ncol, pver) ! T temporary
    real(r8)        :: qv(state%ncol, pver) ! Q temporary
    logical         :: lq(pcnst)            ! Calc tendencies?
    integer         :: lchnk                ! chunk identifier
    integer         :: ncol                 ! number of atmospheric columns
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

    ! Call the Selected condensation routine  ~~DEVO style~~
    !--------------------------------------------------------
    if(CONDENSATE_OPT.eq.CONDENSATE_NONE) then
      call frierson_condensate_NONE(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT.eq.CONDENSATE_FRIERSON) then
      call frierson_condensate(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                               state%pdel(:ncol,:), &
                                                        T(:ncol,:), &
                                                       qv(:ncol,:), &
                                                   relhum(:ncol,:), &
                                                 prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT.eq.CONDENSATE_TJ16) then
      call frierson_condensate_TJ16(ncol,pver,ztodt,state%pmid(:ncol,:), &
                                                    state%pdel(:ncol,:), &
                                                             T(:ncol,:), &
                                                            qv(:ncol,:), &
                                                        relhum(:ncol,:), &
                                                      prec_pcw(:ncol)    )
    elseif(CONDENSATE_OPT.eq.CONDENSATE_USER) then
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

    ! End Routine
    !--------------
    return
  end subroutine frierson_condensate_tend
  !==============================================================================


  !============================================================================
  subroutine frierson_pbl_tend(state, ptend, ztodt, cam_in)
    !
    ! frierson_pbl_tend: Run the selected PBL process.
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
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
    real(r8) :: Km        (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
    real(r8) :: Ke        (state%ncol,pver+1) ! Eddy diffusivity at layer interfaces (m2/s)
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
    if(PBL_OPT.eq.PBL_FRIERSON) then
      ! Call Frierson PBL scheme
      !--------------------------------------------------
      call frierson_pbl(ncol, pver, ztodt,state%pmid (:ncol,:),     &
                                          state%pint (:ncol,:),     &
                                          state%rpdel(:ncol,:),     &
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
                                           dvdt_vdiff(:ncol,:)      )
    elseif(PBL_OPT.eq.PBL_USER) then
      ! Call USER implemented routine in frierson module
      !--------------------------------------------------
      call frierson_pbl_USER(ncol, pver, ztodt,state%pmid (:ncol,:),     &
                                               state%pint (:ncol,:),     &
                                               state%rpdel(:ncol,:),     &
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
                                                dvdt_vdiff(:ncol,:)      )
    else
      ! ERROR: Unknown PBL_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown PBL_OPT=',PBL_OPT
      call endrun('frierson_pbl_tend() PBL_OPT ERROR')
    endif
    Tsurf(:ncol,lchnk) = Tsfc(:ncol)
    Qsurf(:ncol,lchnk) = Qsfc(:ncol)

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
    call outfld('KVH'  ,Ke        ,ncol,lchnk) ! Eddy diffusivity (heat and moisture,m2/s)
    call outfld('KVM'  ,Km        ,ncol,lchnk) ! Eddy diffusivity (momentum, m2/s)
    call outfld('VSE'  ,VSE       ,ncol,lchnk) ! Virtual Dry Static Energy divided by Cp (K)
    call outfld('Zm'   ,Zm        ,ncol,lchnk) ! 
    call outfld('Z_pbl',Z_pbl     ,ncol,lchnk) ! 
    call outfld('Rf'   ,Rf        ,ncol,lchnk) ! 
    call outfld('DUV'  ,dudt_vdiff,ncol,lchnk) ! PBL u tendency (m/s2)
    call outfld('DVV'  ,dvdt_vdiff,ncol,lchnk) ! PBL v tendency (m/s2)
    call outfld('DTV'  ,dtdt_vdiff,ncol,lchnk) ! PBL + surface flux T tendency (K/s)
    call outfld('VD01' ,dqdt_vdiff,ncol,lchnk) ! PBL + surface flux Q tendency (kg/kg/s)
    call outfld('Cdrag',Cdrag     ,ncol,lchnk) ! 

    call outfld('R_Tsurf' , Tsurf (:ncol,lchnk),ncol,lchnk) 
    call outfld('R_Qsurf' , Qsurf (:ncol,lchnk),ncol,lchnk) 

    ! End Routine
    !--------------
    return
  end subroutine frierson_pbl_tend
  !============================================================================


  !============================================================================
  subroutine frierson_radiative_tend(state, ptend, ztodt,cam_in,cam_out)
    !
    ! frierson_radiative_tend: Run the radiatvie process
    !=========================================================================
    use physics_types,only: physics_state, physics_ptend
    use physics_types,only: physics_ptend_init
    use physconst,    only: cpair
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
    real(r8):: dtdt_heating(state%ncol,pver) ! temperature tendency from relaxation in K/s
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

    ! Set Tsfc() PFC?? Get rid of this!!!
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
    if(RADIATION_OPT.eq.RADIATION_FRIERSON) then
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
                                             Fdown(:ncol,lchnk)  )
    elseif(RADIATION_OPT.eq.RADIATION_USER) then
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
                                                  Fdown(:ncol,lchnk)  )
    else
      ! ERROR: Unknown RADIATION_OPT value
      !-------------------------------------
      write(iulog,*) 'ERROR: unknown RADIATION_OPT=',RADIATION_OPT
      call endrun('frierson_pbl_tend() RADIATION_OPT ERROR')
    endif

    Fnet  (:ncol,lchnk) = Fup(:ncol,lchnk) - Fdown (:ncol,lchnk)

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
    call outfld('QRS',dtdt_heating, ncol,lchnk) 

    ! End Routine
    !------------
    return
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

    ! End Routine
    !--------------
    return
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
    ierr = pio_def_var(File,'Frierson_Qsfc',pio_double, hdimids, Qsurf_desc)

    ! End Routine
    !--------------
    return
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
    call pio_write_darray(File, Qsurf_desc, iodesc, Qsurf, ierr)

    ! End Routine
    !--------------
    return
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
    call pio_read_darray(File, vardesc, iodesc, Tsurf, ierr)
    ierr = pio_inq_varid(File,'Frierson_Qsfc',vardesc)
    call pio_read_darray(File, vardesc, iodesc, Qsurf, ierr)

    ! End Routine
    !--------------
    return
  end subroutine frierson_restart_read
  !=======================================================================

end module frierson_cam

