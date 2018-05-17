
  module diffusion_solver

  !------------------------------------------------------------------------------------ !
  ! Module to solve vertical diffusion equations using a tri-diagonal solver.           !
  ! The module will also apply countergradient fluxes, and apply molecular              ! 
  ! diffusion for constituents.                                                         !
  !                                                                                     !
  ! Public interfaces :                                                                 ! 
  !    init_vdiff       initializes time independent coefficients                       !
  !    compute_vdiff    solves diffusion equations                                      !
  !    vdiff_selector   type for storing fields selected to be diffused                 !
  !    vdiff_select     selects fields to be diffused                                   !
  !    operator(.not.)  extends .not. to operate on type vdiff_selector                 !
  !    any              provides functionality of intrinsic any for type vdiff_selector !
  !                                                                                     !
  !------------------------------------ Code History ---------------------------------- !
  ! Initial subroutines :  B. Boville and others, 1991-2004                             !
  ! Modularization      :  J. McCaa, September 2004                                     !
  ! Most Recent Code    :  Sungsu Park, Aug. 2006, Dec. 2008, Jan. 2010.                !
  !------------------------------------------------------------------------------------ !

  implicit none
  private       
  save

  integer, parameter :: r8 = selected_real_kind(12)      ! 8 byte real

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public init_vdiff                                      ! Initialization
  public new_fieldlist_vdiff                             ! Returns an empty fieldlist
  public compute_vdiff                                   ! Full routine
  public vdiff_selector                                  ! Type for storing fields selected to be diffused
  public vdiff_select                                    ! Selects fields to be diffused
  public operator(.not.)                                 ! Extends .not. to operate on type vdiff_selector
  public any                                             ! Provides functionality of intrinsic any for type vdiff_selector
 
  ! Below stores logical array of fields to be diffused

  type vdiff_selector 
       private
       logical, allocatable, dimension(:) :: fields
  end type vdiff_selector

  ! Below extends .not. to operate on type vdiff_selector

  interface operator(.not.)
       module procedure not
  end interface

  ! Below provides functionality of intrinsic any for type vdiff_selector

  interface any                           
       module procedure my_any
  end interface

  ! ------------ !
  ! Private data !
  ! ------------ !

  ! Unit number for log output
  integer :: iulog = -1

  real(r8), private   :: cpair                           ! Specific heat of dry air
  real(r8), private   :: gravit                          ! Acceleration due to gravity
  real(r8), private   :: rair                            ! Gas constant for dry air

  logical,  private   :: do_iss                          ! Use implicit turbulent surface stress computation

  ! Parameters used for Turbulent Mountain Stress

  real(r8), parameter :: z0fac   = 0.025_r8              ! Factor determining z_0 from orographic standard deviation
  real(r8), parameter :: z0max   = 100._r8               ! Max value of z_0 for orography
  real(r8), parameter :: horomin = 10._r8                ! Min value of subgrid orographic height for mountain stress
  real(r8), parameter :: dv2min  = 0.01_r8               ! Minimum shear squared

  logical :: am_correction ! logical switch for AM correction 

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine init_vdiff( kind, iulog_in, rair_in, cpair_in, gravit_in, do_iss_in, &
                         am_correction_in, errstring )

    integer,              intent(in)  :: kind            ! Kind used for reals
    integer,              intent(in)  :: iulog_in        ! Unit number for log output.
    real(r8),             intent(in)  :: rair_in         ! Input gas constant for dry air
    real(r8),             intent(in)  :: cpair_in        ! Input heat capacity for dry air
    real(r8),             intent(in)  :: gravit_in       ! Input gravitational acceleration
    logical,              intent(in)  :: do_iss_in       ! Input ISS flag
    logical,              intent(in)  :: am_correction_in! for angular momentum conservation
    character(128),       intent(out) :: errstring       ! Output status
    
    errstring = ''
    iulog = iulog_in
    if( kind .ne. r8 ) then
        write(iulog,*) 'KIND of reals passed to init_vdiff -- exiting.'
        errstring = 'init_vdiff'
        return
    endif

    rair   = rair_in
    cpair  = cpair_in
    gravit = gravit_in
    do_iss = do_iss_in
    am_correction = am_correction_in

  end subroutine init_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  type(vdiff_selector) pure function new_fieldlist_vdiff(ncnst)

    integer,              intent(in)  :: ncnst           ! Number of constituents

    allocate( new_fieldlist_vdiff%fields( 3 + ncnst ) )
    new_fieldlist_vdiff%fields = .false.

  end function new_fieldlist_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine compute_vdiff( lchnk           ,                                                                   &
                            pcols           , pver               , ncnst         , ncol         , tint        , &
                            p               , t                  , rhoi          , ztodt        , taux        , &
                            tauy            , shflx              , cflx          , &
                            kvh             , kvm                , kvq           , cgs          , cgh         , &
                            zi              , ksrftms            , dragblj       , & 
                            qmincg          , fieldlist          , fieldlistm    , &
                            u               , v                  , q             , dse          ,               &
                            tautmsx         , tautmsy            , dtk           , topflx       , errstring   , &
                            tauresx         , tauresy            , itaures       , cpairv       , dse_top, &
                            do_molec_diff   , use_temperature_molec_diff, vd_lu_qdecomp, &
                            ubc_mmr, ubc_flux, kvt, pmid, &
                            cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx, nbot_molec, &
                            kq_scal, mw_fac)

    !-------------------------------------------------------------------------- !
    ! Driver routine to compute vertical diffusion of momentum, moisture, trace !
    ! constituents and dry static energy. The new temperature is computed from  !
    ! the diffused dry static energy.                                           ! 
    ! Turbulent diffusivities and boundary layer nonlocal transport terms are   !
    ! obtained from the turbulence module.                                      !
    !-------------------------------------------------------------------------- !

!    Used for CAM debugging.
!    use phys_debug_util,    only : phys_debug_col
!    use time_manager,       only : is_first_step, get_nstep

    use coords_1d, only: Coords1D
    use linear_1d_operators, only : BoundaryType, BoundaryFixedLayer, &
         BoundaryData, BoundaryFlux, TriDiagDecomp
    use vdiff_lu_solver,     only : fin_vol_lu_decomp
    use beljaars_drag_cam,   only : do_beljaars
    ! FIXME: This should not be needed
    use physconst, only: rairv
  
    use phys_control,        only : phys_getopts 
 
  ! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
  !                Also, vertical diffusion of cloud droplet number concentration and aerosol number
  !                concentration should be done very carefully in the future version.

    ! --------------- !
    ! Input Arguments !
    ! --------------- !

    integer,  intent(in)    :: lchnk                   
    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy'
                                                         ! is updated (1) or non-updated (0) in this subroutine.

    type(Coords1D), intent(in) :: p                      ! Pressure coordinates [ Pa ]
    real(r8), intent(in)    :: tint(pcols,pver+1)        ! Temperature [ K ]
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [ K ]
    real(r8), intent(in)    :: rhoi(pcols,pver+1)        ! Density of air at interfaces [ kg/m3 ]
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t [ s ]
    real(r8), intent(in)    :: taux(pcols)               ! Surface zonal      stress.
                                                         ! Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: tauy(pcols)               ! Surface meridional stress.
                                                         ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: shflx(pcols)              ! Surface sensible heat flux [ W/m2 ]
    real(r8), intent(in)    :: cflx(pcols,ncnst)         ! Surface constituent flux [ kg/m2/s ]
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights [ m ]
    real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
    real(r8), intent(in)    :: dragblj(pcols,pver)       ! Drag profile from Beljaars SGO form drag  > 0. [ 1/s ]
    real(r8), intent(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes
    real(r8), intent(in)    :: cpairv(pcols,pver)        ! Specific heat at constant pressure
    real(r8), intent(in)    :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]

    logical,  intent(in)    :: do_molec_diff             ! Flag indicating multiple constituent diffusivities
    logical,  intent(in)    :: use_temperature_molec_diff! Flag indicating that molecular diffusion should apply to temperature, not
                                                         ! dry static energy.

    type(vdiff_selector), intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse
    type(vdiff_selector), intent(in) :: fieldlistm       ! Array of flags selecting which fields for molecular diffusion

    ! Dry static energy top boundary condition.
    real(r8), intent(in) :: dse_top(pcols)

    real(r8), intent(in) :: kvm(pcols,pver+1)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    real(r8), intent(in) :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents
    real(r8), intent(in) :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
    real(r8), intent(in) :: cgh(pcols,pver+1)         ! Counter-gradient term for heat

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    real(r8), intent(inout) :: u(pcols,pver)             ! U wind. This input is the 'raw' input wind to
                                                         ! PBL scheme without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: v(pcols,pver)             ! V wind. This input is the 'raw' input wind to PBL scheme
                                                         ! without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: q(pcols,pver,ncnst)       ! Moisture and trace constituents [ kg/kg, #/kg ? ]
    real(r8), intent(inout) :: dse(pcols,pver)           ! Dry static energy [ J/kg ]

    real(r8), intent(inout) :: tauresx(pcols)            ! Input  : Reserved surface stress at previous time step
    real(r8), intent(inout) :: tauresy(pcols)            ! Output : Reserved surface stress at current  time step

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    real(r8), intent(out)   :: dtk(pcols,pver)           ! T tendency from KE dissipation
    real(r8), intent(out)   :: tautmsx(pcols)            ! Implicit zonal      turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: tautmsy(pcols)            ! Implicit meridional turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: topflx(pcols)             ! Molecular heat flux at the top interface
    character(128), intent(out) :: errstring             ! Output status

    ! ------------------ !
    ! Optional Arguments !
    ! ------------------ !

    ! The molecular diffusion module will likely change significantly in
    ! the future, and this module may directly depend on it after that.
    ! Until then, we have these highly specific interfaces hard-coded.

    optional :: vd_lu_qdecomp        ! Constituent-dependent molecular diffusivity routine

    interface
       function vd_lu_qdecomp( &
            pcols , pver   , ncol       , fixed_ubc  , mw     , &
            kv    , kq_scal, mw_facm    , dpidz_sq   , coords , &
            interface_boundary, molec_boundary, &
            tint  , ztodt  , nbot_molec , &
            lchnk , t          , m      , no_molec_decomp) result(decomp)
         import
         integer,  intent(in)    :: pcols
         integer,  intent(in)    :: pver
         integer,  intent(in)    :: ncol
         integer,  intent(in)    :: nbot_molec
         logical,  intent(in)    :: fixed_ubc
         real(r8), intent(in)    :: kv(pcols,pver+1)
         real(r8), intent(in)    :: kq_scal(pcols,pver+1)
         real(r8), intent(in)    :: mw
         real(r8), intent(in)    :: mw_facm(pcols,pver+1)
         real(r8), intent(in)    :: dpidz_sq(ncol,pver+1)
         type(Coords1D), intent(in) :: coords
         type(BoundaryType), intent(in) :: interface_boundary
         type(BoundaryType), intent(in) :: molec_boundary
         real(r8), intent(in)    :: tint(pcols,pver+1)
         real(r8), intent(in)    :: ztodt
         integer,  intent(in)    :: lchnk
         real(r8), intent(in)    :: t(pcols,pver)
         integer,  intent(in)    :: m
         type(TriDiagDecomp), intent(in) :: no_molec_decomp
         type(TriDiagDecomp) :: decomp
       end function vd_lu_qdecomp
    end interface

    real(r8), intent(in), optional :: ubc_mmr(pcols,ncnst) ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(in), optional :: ubc_flux(pcols,ncnst)      ! Upper boundary flux [ kg/s/m^2 ]

    real(r8), intent(in), optional :: kvt(pcols,pver+1) ! Kinematic molecular conductivity

    ! FIXME: This input should not be needed (and should not be passed in in vertical_diffusion).
    real(r8), intent(in), optional :: pmid(pcols,pver)

    real(r8), intent(in), optional :: cnst_mw(ncnst)          ! Molecular weight [ kg/kmole ]
    logical, intent(in), optional  :: cnst_fixed_ubc(ncnst)   ! Whether upper boundary condition is fixed
    logical, intent(in), optional  :: cnst_fixed_ubflx(ncnst) ! Whether upper boundary flux is a fixed non-zero value

    integer, intent(in), optional  :: nbot_molec              ! Bottom level where molecular diffusivity is applied

    ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(in), optional :: kq_scal(pcols,pver+1)
    ! Local sqrt(1/M_q + 1/M_d) for each constituent
    real(r8), intent(in), optional :: mw_fac(pcols,pver+1,ncnst)

    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    integer  :: i, k, m                                  ! Longitude, level, constituent indices
    logical  :: lqtst(pcols)                             ! Adjust vertical profiles

    ! LU decomposition information.
    type(TriDiagDecomp) :: decomp
    type(TriDiagDecomp) :: no_molec_decomp

    ! Square of derivative of pressure with height (on interfaces).
    real(r8) :: dpidz_sq(ncol,pver+1)

    ! Pressure coordinates over the molecular diffusion range only.
    type(Coords1D) :: p_molec

    ! Boundary layer objects
    type(BoundaryType) :: interface_boundary
    type(BoundaryType) :: molec_boundary

    real(r8) :: tmp1(pcols)                              ! Temporary storage
    real(r8) :: tmpi1(pcols,pver+1)                      ! Interface KE dissipation
    real(r8) :: tmpi2(pcols,pver+1)                      ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: keg_in(pcols,pver)                       ! KE on entry to subroutine
    real(r8) :: keg_out(pcols,pver)                      ! KE after U and V dissipation/diffusion
    real(r8) :: rrho(pcols)                              ! 1./bottom level density 

    real(r8) :: tautotx(pcols)                           ! Total surface stress ( zonal )
    real(r8) :: tautoty(pcols)                           ! Total surface stress ( meridional )

    real(r8) :: dinp_u(pcols,pver+1)                     ! Vertical difference at interfaces, input u
    real(r8) :: dinp_v(pcols,pver+1)                     ! Vertical difference at interfaces, input v
    real(r8) :: dout_u                                   ! Vertical difference at interfaces, output u
    real(r8) :: dout_v                                   ! Vertical difference at interfaces, output v

    real(r8) :: qtm(pcols,pver)                          ! Temporary copy of q

    real(r8) :: ws(pcols)                                ! Lowest-level wind speed [ m/s ]
    real(r8) :: tau(pcols)                               ! Turbulent surface stress ( not including mountain stress )
    real(r8) :: ksrfturb(pcols)                          ! Surface drag coefficient of 'normal' stress. > 0.
                                                         ! Virtual mass input per unit time per unit area [ kg/s/m2 ]
    real(r8) :: ksrf(pcols)                              ! Surface drag coefficient of 'normal' stress +
                                                         ! Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
    real(r8) :: usum_in(pcols)                           ! Vertical integral of input u-momentum. Total zonal
                                                         ! momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    real(r8) :: vsum_in(pcols)                           ! Vertical integral of input v-momentum. Total meridional
                                                         ! momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    real(r8) :: usum_mid(pcols)                          ! Vertical integral of u-momentum after adding explicit residual stress
    real(r8) :: vsum_mid(pcols)                          ! Vertical integral of v-momentum after adding explicit residual stress
    real(r8) :: usum_out(pcols)                          ! Vertical integral of u-momentum after doing implicit diffusion
    real(r8) :: vsum_out(pcols)                          ! Vertical integral of v-momentum after doing implicit diffusion
    real(r8) :: tauimpx(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: tauimpy(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: ramda                                    ! dt/timeres [ no unit ]

    real(r8) :: taubljx(pcols)                           ! recomputed explicit/residual beljaars stress
    real(r8) :: taubljy(pcols)                           ! recomputed explicit/residual beljaars stress

    ! Rate at which external (surface) stress damps wind speeds (1/s).
    real(r8) :: tau_damp_rate(ncol, pver)

    ! Combined molecular and eddy diffusion.
    real(r8) :: kv_total(pcols,pver+1)

    logical  :: use_spcam

    !--------------------------------
    ! Variables needed for WACCM-X
    !--------------------------------
    real(r8) :: ttemp(ncol,pver)			 ! temporary temperature array
    real(r8) :: ttemp0(ncol,pver)			 ! temporary temperature array

    ! ------------------------------------------------ !
    ! Parameters for implicit surface stress treatment !
    ! ------------------------------------------------ !

    real(r8), parameter :: wsmin = 1._r8                 ! Minimum sfc wind speed for estimating frictional
                                                         ! transfer velocity ksrf. [ m/s ]
    real(r8), parameter :: ksrfmin = 1.e-4_r8            ! Minimum surface drag coefficient [ kg/s/m^2 ]
    real(r8), parameter :: timeres = 7200._r8            ! Relaxation time scale of residual stress ( >= dt ) [ s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    call phys_getopts(use_spcam_out = use_spcam)

    errstring = ''
    if( ( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) .and. .not. diffuse(fieldlist,'s') ) then
          errstring = 'diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v'
          return
    end if

    !--------------------------------------- !
    ! Computation of Molecular Diffusivities !
    !--------------------------------------- !

  ! Modification : Why 'kvq' is not changed by molecular diffusion ? 

    if( do_molec_diff ) then

        if( (.not.present(vd_lu_qdecomp)) .or. (.not.present(kvt)) &
             .or. (.not. present(ubc_mmr)) .or. (.not. present(ubc_flux)) ) then
              errstring = 'compute_vdiff: do_molec_diff true but vd_lu_qdecomp or kvt missing'
              return
        endif

        p_molec = p%section([1, ncol], [1, nbot_molec])
        molec_boundary = BoundaryFixedLayer(p%del(:,nbot_molec+1))

    endif

    ! Boundary condition for a fixed concentration directly on a boundary
    ! interface (i.e. a boundary layer of size 0).
    interface_boundary = BoundaryFixedLayer(spread(0._r8, 1, ncol))

    ! Note that the *derivative* dp/dz is g*rho
    dpidz_sq = gravit*rhoi(:ncol,:)
    dpidz_sq = dpidz_sq * dpidz_sq

    rrho(:ncol) = rair  * t(:ncol,pver) / p%mid(:,pver)

    tmpi2(:ncol,1) = ztodt * dpidz_sq(:,1) / ( p%mid(:,1) - p%ifc(:,1) )
    tmpi2(:ncol,2:pver) = ztodt * dpidz_sq(:,2:pver) * p%rdst

    ! FIXME: The following four lines are kept in only to preserve answers;
    !        they really should be taken out completely.
    if (do_molec_diff) &
         tmpi2(:ncol,1) = ztodt * (gravit * rhoi(:ncol,1))**2 / ( pmid(:ncol,1) - p%ifc(:,1) )
    dpidz_sq(:,1) = gravit*(p%ifc(:,1) / (rairv(:ncol,1,lchnk)*t(:ncol,1)))
    dpidz_sq(:,1) = dpidz_sq(:,1)*dpidz_sq(:,1)

    tmp1(:ncol) = ztodt * gravit * p%rdel(:,pver)

    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !

    do k = 1, pver
       do i = 1, ncol
          keg_in(i,k) = 0.5_r8 * ( u(i,k)*u(i,k) + v(i,k)*v(i,k) )
       end do
    end do

    if( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) then

        ! Compute the vertical upward differences of the input u,v for KE dissipation
        ! at each interface.
        ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
        ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

        do i = 1, ncol
           dinp_u(i,1) = 0._r8
           dinp_v(i,1) = 0._r8
           dinp_u(i,pver+1) = -u(i,pver)
           dinp_v(i,pver+1) = -v(i,pver)
        end do
        do k = 2, pver
           do i = 1, ncol
              dinp_u(i,k) = u(i,k) - u(i,k-1)
              dinp_v(i,k) = v(i,k) - v(i,k-1)
           end do
        end do

       ! -------------------------------------------------------------- !
       ! Do 'Implicit Surface Stress' treatment for numerical stability !
       ! in the lowest model layer.                                     !
       ! -------------------------------------------------------------- !

       if( do_iss ) then

         ! Compute surface drag coefficient for implicit diffusion 
         ! including turbulent mountain stress. 

           do i = 1, ncol
              ws(i)       = max( sqrt( u(i,pver)**2._r8 + v(i,pver)**2._r8 ), wsmin )
              tau(i)      = sqrt( taux(i)**2._r8 + tauy(i)**2._r8 )
              ksrfturb(i) = max( tau(i) / ws(i), ksrfmin )
           end do
           ksrf(:ncol) = ksrfturb(:ncol) + ksrftms(:ncol)  ! Do all surface stress ( normal + tms ) implicitly

         ! Vertical integration of input momentum. 
         ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
         ! Note (u,v) are the raw input to the PBL scheme, not the
         ! provisionally-marched ones within the iteration loop of the PBL scheme.  

           do i = 1, ncol
              usum_in(i) = 0._r8
              vsum_in(i) = 0._r8
              do k = 1, pver
                 usum_in(i) = usum_in(i) + (1._r8/gravit)*u(i,k)*p%del(i,k)
                 vsum_in(i) = vsum_in(i) + (1._r8/gravit)*v(i,k)*p%del(i,k)
              end do
           end do              

         ! Add residual stress of previous time step explicitly into the lowest
         ! model layer with a relaxation time scale of 'timeres'.

           if (am_correction) then
              ! preserve time-mean torque 
              ramda         = 1._r8
           else
              ramda         = ztodt / timeres
           endif

           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*tauresx(:ncol)*ramda
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauresy(:ncol)*ramda

         ! Vertical integration of momentum after adding explicit residual stress
         ! into the lowest model layer.

           do i = 1, ncol
              usum_mid(i) = 0._r8
              vsum_mid(i) = 0._r8
              do k = 1, pver
                 usum_mid(i) = usum_mid(i) + (1._r8/gravit)*u(i,k)*p%del(i,k)
                 vsum_mid(i) = vsum_mid(i) + (1._r8/gravit)*v(i,k)*p%del(i,k)
              end do
           end do              

       else

         ! In this case, do 'turbulent mountain stress' implicitly, 
         ! but do 'normal turbulent stress' explicitly.
         ! In this case, there is no 'residual stress' as long as 'tms' is
         ! treated in a fully implicit way, which is true.

         ! 1. Do 'tms' implicitly

           ksrf(:ncol) = ksrftms(:ncol) 

         ! 2. Do 'normal stress' explicitly

           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*taux(:ncol)
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauy(:ncol)

       end if  ! End of 'do iss' ( implicit surface stress )

       ! --------------------------------------------------------------------------------------- !
       ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
       ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
       !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
       !                                   u(pver) : explicitly include 'residual normal' stress !
       !    For explicit 'normal' stress : ksrf = ksrftms                                        !
       !                                   u(pver) : explicitly include 'normal' stress          !
       ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
       ! --------------------------------------------------------------------------------------- !

       ! In most layers, no damping at all.
       tau_damp_rate = 0._r8

       ! Physical interpretation:
       ! ksrf is stress per unit wind speed.
       ! p%del / gravit is approximately the mass in the layer per unit of
       ! surface area.
       ! Therefore, gravit*ksrf/p%del is the acceleration of wind per unit
       ! wind speed, i.e. the rate at which wind is exponentially damped by
       ! surface stress.

       ! Beljaars et al SGO scheme incorporated here. It 
       ! appears as a "3D" tau_damp_rate specification.

       tau_damp_rate(:,pver) = -gravit*ksrf(:ncol)*p%rdel(:,pver)
       do k=1,pver
          tau_damp_rate(:,k) = tau_damp_rate(:,k) + dragblj(:ncol,k)
       end do

       decomp = fin_vol_lu_decomp(ztodt, p, &
            coef_q=tau_damp_rate, coef_q_diff=kvm(:ncol,:)*dpidz_sq)

       call decomp%left_div(u(:ncol,:))
       call decomp%left_div(v(:ncol,:))
       call decomp%finalize()

       ! ---------------------------------------------------------------------- !
       ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
       ! have been actually added into the atmosphere at the current time step. ! 
       ! Also, update residual stress, if required.                             !
       ! ---------------------------------------------------------------------- !

       do i = 1, ncol

          ! Compute the implicit 'tms' using the updated winds.
          ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
          ! that has been actually added into the atmosphere both for explicit
          ! and implicit approach. 

          tautmsx(i) = -ksrftms(i)*u(i,pver)
          tautmsy(i) = -ksrftms(i)*v(i,pver)

          ! We want to add vertically-integrated Beljaars drag to residual stress.
          ! So this has to be calculated locally. 
          ! We may want to rethink the residual drag calculation performed here on. (jtb)
          taubljx(i) = 0._r8
          taubljy(i) = 0._r8
          do k = 1, pver
             taubljx(i) = taubljx(i) + (1._r8/gravit)*dragblj(i,k)*u(i,k)*p%del(i,k)
             taubljy(i) = taubljy(i) + (1._r8/gravit)*dragblj(i,k)*v(i,k)*p%del(i,k)
          end do
        
          if( do_iss ) then

            ! Compute vertical integration of final horizontal momentum

              usum_out(i) = 0._r8
              vsum_out(i) = 0._r8
              do k = 1, pver
                 usum_out(i) = usum_out(i) + (1._r8/gravit)*u(i,k)*p%del(i,k)
                 vsum_out(i) = vsum_out(i) + (1._r8/gravit)*v(i,k)*p%del(i,k)
              end do

            ! Compute net stress added into the atmosphere at the current time step.
            ! Note that the difference between 'usum_in' and 'usum_out' are induced
            ! by 'explicit residual stress + implicit total stress' for implicit case, while
            ! by 'explicit normal   stress + implicit tms   stress' for explicit case. 
            ! Here, 'tautotx(i)' is net stress added into the air at the current time step.

              tauimpx(i) = ( usum_out(i) - usum_in(i) ) / ztodt
              tauimpy(i) = ( vsum_out(i) - vsum_in(i) ) / ztodt

              tautotx(i) = tauimpx(i) 
              tautoty(i) = tauimpy(i) 

            ! Compute residual stress and update if required.
            ! Note that the total stress we should have added at the current step is
            ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.

              if( itaures .eq. 1 ) then
                 tauresx(i) = taux(i) + tautmsx(i) + taubljx(i) + tauresx(i)- tauimpx(i)
                 tauresy(i) = tauy(i) + tautmsy(i) + taubljy(i) + tauresy(i)- tauimpy(i)
              endif

          else

             tautotx(i) = tautmsx(i) + taux(i)
             tautoty(i) = tautmsy(i) + tauy(i)
             tauresx(i) = 0._r8
             tauresy(i) = 0._r8

          end if  ! End of 'do_iss' if

       end do ! End of 'do i = 1, ncol' loop

       ! ------------------------------------ !
       ! Calculate kinetic energy dissipation !
       ! ------------------------------------ !       

     ! Modification : In future, this should be set exactly same as 
     !                the ones in the convection schemes 

       ! 1. Compute dissipation term at interfaces
       !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are 
       !    implicit stress that has been actually added. On the other hand,
       !    'dinp_u, dinp_v' were computed using non-diffused input wind.

     ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
     !                is correctly intended approach. I think so.

       k = pver + 1
       do i = 1, ncol
          tmpi1(i,1) = 0._r8
          tmpi1(i,k) = 0.5_r8 * ztodt * gravit * &
                       ( (-u(i,k-1) + dinp_u(i,k))*tautotx(i) + (-v(i,k-1) + dinp_v(i,k))*tautoty(i) )
       end do

       do k = 2, pver
          do i = 1, ncol
             dout_u = u(i,k) - u(i,k-1)
             dout_v = v(i,k) - v(i,k-1)
             tmpi1(i,k) = 0.25_r8 * tmpi2(i,k) * kvm(i,k) * &
                          ( dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k) )
          end do
       end do

       if (do_beljaars) then

          ! 2. Add Kinetic Energy change across dissipation to Static Energy
          do k = 1, pver
             do i = 1, ncol
                keg_out(i,k) = 0.5_r8 * ( u(i,k)*u(i,k) + v(i,k)*v(i,k) )
             end do
          end do
    
          do k = 1, pver
             do i = 1, ncol
                dtk(i,k) = keg_in(i,k) - keg_out(i,k)
                dse(i,k) = dse(i,k) + dtk(i,k) ! + dkeblj(i,k)
             end do
          end do

       else

          ! 2. Compute dissipation term at midpoints, add to dry static energy
          do k = 1, pver
             do i = 1, ncol
                dtk(i,k) = ( tmpi1(i,k+1) + tmpi1(i,k) ) * p%rdel(i,k)
                dse(i,k) = dse(i,k) + dtk(i,k)
             end do
          end do

       end if

    end if ! End of diffuse horizontal momentum, diffuse(fieldlist,'u') routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

  ! Modification : In future, we should diffuse the fully conservative 
  !                moist static energy,not the dry static energy.

    if( diffuse(fieldlist,'s') ) then
      if (.not. use_spcam) then

       ! Add counter-gradient to input static energy profiles

         do k = 1, pver
            dse(:ncol,k) = dse(:ncol,k) + ztodt * p%rdel(:,k) * gravit  *                &
                                        ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgh(:ncol,k+1) &
                                        - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgh(:ncol,k  ) )
         end do
       endif
       ! Add the explicit surface fluxes to the lowest layer
       dse(:ncol,pver) = dse(:ncol,pver) + tmp1(:ncol) * shflx(:ncol)

     ! Diffuse dry static energy

       !---------------------------------------------------
       ! Solve for temperature using thermal conductivity 
       !---------------------------------------------------
       if ( use_temperature_molec_diff ) then 
          !----------------------------------------------------------------------------------------------------
          ! In Extended WACCM, kvt is calculated rather kvh. This is because molecular diffusion operates on 
          ! temperature, while eddy diffusion operates on dse.  Also, pass in constituent dependent "constants"
          !----------------------------------------------------------------------------------------------------

          ! Boundary layer thickness of "0._r8" signifies that the boundary
          ! condition is defined directly on the top interface.
          decomp = fin_vol_lu_decomp(ztodt, p, &
               coef_q_diff=kvh(:ncol,:)*dpidz_sq, &
               upper_bndry=interface_boundary)

          if (.not. use_spcam) then
           call decomp%left_div(dse(:ncol,:), &
                l_cond=BoundaryData(dse_top(:ncol)))
          endif

          call decomp%finalize()

          ! Calculate flux at top interface

          ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

          topflx(:ncol) =  - kvh(:ncol,1) * tmpi2(:ncol,1) / (ztodt*gravit) * &
               ( dse(:ncol,1) - dse_top(:ncol) )

          decomp = fin_vol_lu_decomp(ztodt, p, &
               coef_q_diff=kvt(:ncol,:)*dpidz_sq, &
               coef_q_weight=cpairv(:ncol,:))

          ttemp0 = t(:ncol,:)
          ttemp = ttemp0

          ! upper boundary is zero flux for extended model
          if (.not. use_spcam) then
             call decomp%left_div(ttemp)
          end if

          call decomp%finalize()

          !-------------------------------------
          !  Update dry static energy
          !-------------------------------------
          do k = 1,pver
             dse(:ncol,k) = dse(:ncol,k) + &
                  cpairv(:ncol,k)*(ttemp(:,k) - ttemp0(:,k))
          enddo

       else

          if (do_molec_diff) then
             kv_total(:ncol,:) = kvh(:ncol,:) + kvt(:ncol,:)/cpair
          else
             kv_total(:ncol,:) = kvh(:ncol,:)
          end if

          ! Boundary layer thickness of "0._r8" signifies that the boundary
          ! condition is defined directly on the top interface.
          decomp = fin_vol_lu_decomp(ztodt, p, &
               coef_q_diff=kv_total(:ncol,:)*dpidz_sq, &
               upper_bndry=interface_boundary)

          if (.not. use_spcam) then
             call decomp%left_div(dse(:ncol,:), &
                  l_cond=BoundaryData(dse_top(:ncol)))
          end if

          call decomp%finalize()

          ! Calculate flux at top interface

          ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

          if( do_molec_diff ) then
             topflx(:ncol) =  - kv_total(:ncol,1) * tmpi2(:ncol,1) / (ztodt*gravit) * &
                  ( dse(:ncol,1) - dse_top(:ncol) )
          else
             topflx(:ncol) = 0._r8
          end if

       endif

    endif

    !---------------------------- !
    ! Diffuse Water Vapor Tracers !
    !---------------------------- !

  ! Modification : For aerosols, I need to use separate treatment 
  !                for aerosol mass and aerosol number. 

    ! Loop through constituents

    no_molec_decomp = fin_vol_lu_decomp(ztodt, p, &
         coef_q_diff=kvq(:ncol,:)*dpidz_sq)

    do m = 1, ncnst

       if( diffuse(fieldlist,'q',m) ) then
           if (.not. use_spcam) then

              ! Add the nonlocal transport terms to constituents in the PBL.
              ! Check for neg q's in each constituent and put the original vertical
              ! profile back if a neg value is found. A neg value implies that the
              ! quasi-equilibrium conditions assumed for the countergradient term are
              ! strongly violated.
   
              qtm(:ncol,:pver) = q(:ncol,:pver,m)
   
              do k = 1, pver
                 q(:ncol,k,m) = q(:ncol,k,m) + &
                                ztodt * p%rdel(:,k) * gravit  * ( cflx(:ncol,m) * rrho(:ncol) ) * &
                              ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgs(:ncol,k+1)                 &
                              - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgs(:ncol,k  ) )
              end do
              lqtst(:ncol) = all(q(:ncol,1:pver,m) >= qmincg(m), 2)
              do k = 1, pver
                 q(:ncol,k,m) = merge( q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol) )
              end do
           endif

           ! Add the explicit surface fluxes to the lowest layer

           q(:ncol,pver,m) = q(:ncol,pver,m) + tmp1(:ncol) * cflx(:ncol,m)

           ! Diffuse constituents.

           ! This is for solving molecular diffusion of minor species, thus, for WACCM-X, bypass O and O2 (major species)
           ! Major species diffusion is calculated separately.  -Hanli Liu

           if( do_molec_diff .and. diffuse(fieldlistm,'q',m)) then

              decomp = vd_lu_qdecomp( pcols , pver   , ncol              , cnst_fixed_ubc(m), cnst_mw(m), &
                   kvq   , kq_scal, mw_fac(:,:,m) ,dpidz_sq          , p_molec, &
                   interface_boundary, molec_boundary, &
                   tint  , ztodt  , nbot_molec       , &
                   lchnk , t                , m         , no_molec_decomp)

              ! This to calculate the upper boundary flux of H.    -Hanli Liu
              if ((cnst_fixed_ubflx(m))) then

                 ! ubc_flux is a flux of mass density through space, i.e.:
                 ! ubc_flux = rho_i * dz/dt = q_i * rho * dz/dt
                 ! For flux of mmr through pressure level, multiply by g:
                 ! q_i * rho * gravit * dz/dt = q_i * dp/dt

                 call decomp%left_div(q(:ncol,:,m), &
                      l_cond=BoundaryFlux( &
                      -gravit*ubc_flux(:ncol,m), ztodt, &
                      p%del(:,1)))

              else
                 call decomp%left_div(q(:ncol,:,m), &
                      l_cond=BoundaryData(ubc_mmr(:ncol,m)))
              end if

              call decomp%finalize()

           else

              if (.not. use_spcam) then
                 ! Currently, no ubc for constituents without molecular
                 ! diffusion (they cannot diffuse out the top of the model).
                 call no_molec_decomp%left_div(q(:ncol,:,m))
              end if

           end if

       end if
    end do

    call no_molec_decomp%finalize()

  end subroutine compute_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  character(128) function vdiff_select( fieldlist, name, qindex )
    ! --------------------------------------------------------------------- !
    ! This function sets the field with incoming name as one to be diffused !
    ! --------------------------------------------------------------------- !
    type(vdiff_selector), intent(inout)        :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    vdiff_select = ''
    select case (name)
    case ('u','U')
       fieldlist%fields(1) = .true.
    case ('v','V')
       fieldlist%fields(2) = .true.
    case ('s','S')
       fieldlist%fields(3) = .true.
    case ('q','Q')
       if( present(qindex) ) then
           fieldlist%fields(3 + qindex) = .true.
       else
           fieldlist%fields(4) = .true.
       endif
    case default
       write(vdiff_select,*) 'Bad argument to vdiff_index: ', name
    end select
    return
    
  end function vdiff_select

  type(vdiff_selector) function not(a)
    ! ------------------------------------------------------------- !
    ! This function extends .not. to operate on type vdiff_selector !
    ! ------------------------------------------------------------- !    
    type(vdiff_selector), intent(in)  :: a
    allocate(not%fields(size(a%fields)))
    not%fields = .not. a%fields
  end function not

  logical function my_any(a)
    ! -------------------------------------------------- !
    ! This function extends the intrinsic function 'any' ! 
    ! to operate on type vdiff_selector                  ! 
    ! -------------------------------------------------- !
    type(vdiff_selector), intent(in) :: a
    my_any = any(a%fields)
  end function my_any

  logical function diffuse(fieldlist,name,qindex)
    ! ---------------------------------------------------------------------------- !
    ! This function reports whether the field with incoming name is to be diffused !
    ! ---------------------------------------------------------------------------- !
    type(vdiff_selector), intent(in)           :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    select case (name)
    case ('u','U')
       diffuse = fieldlist%fields(1)
    case ('v','V')
       diffuse = fieldlist%fields(2)
    case ('s','S')
       diffuse = fieldlist%fields(3)
    case ('q','Q')
       if( present(qindex) ) then
           diffuse = fieldlist%fields(3 + qindex)
       else
           diffuse = fieldlist%fields(4)
       endif
    case default
       diffuse = .false.
    end select
    return
  end function diffuse

end module diffusion_solver
