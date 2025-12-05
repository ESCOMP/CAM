
  module diffusion_solver_cam

  !------------------------------------------------------------------------------------ !
  ! Module to solve vertical diffusion equations using a tri-diagonal solver.           !
  ! The module will also apply countergradient fluxes, and apply molecular              !
  ! diffusion for constituents.                                                         !
  !                                                                                     !
  ! Public interfaces :                                                                 !
  !    init_vdiff       initializes time independent coefficients                       !
  !    compute_vdiff    solves diffusion equations                                      !
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

  public :: init_vdiff                                      ! Initialization
  public :: compute_vdiff                                   ! Full routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  ! Unit number for log output
  integer :: iulog = -1

  real(r8), private   :: cpair                           ! Specific heat of dry air
  real(r8), private   :: gravit                          ! Acceleration due to gravity
  real(r8), private   :: rair                            ! Gas constant for dry air

  logical,  private   :: do_iss                          ! Use implicit turbulent surface stress computation
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

  subroutine compute_vdiff( ncol, pver, pverp, ncnst, tint        , &
                            p               , t                  , rhoi          , ztodt        , taux        , &
                            tauy            , shflx              , cflx          , &
                            kvh             , kvm                , kvq           , cgs          , cgh         , &
                            ksrftms         , dragblj         , &
                            qmincg          , &
                            do_diffusion_u_v, do_diffusion_s, &
                            do_diffusion_const, &
                            do_molecular_diffusion_const, &
                            u               , v                  , q             , dse          ,               &
                            tautmsx         , tautmsy            , dtk           , topflx       , errmsg   , &
                            tauresx         , tauresy            , itaures       , &
                            cpairv, rairv, mbarv, &
                            dse_top, &
                            do_beljaars, &
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

    use coords_1d, only: Coords1D
    use linear_1d_operators, only : BoundaryType, BoundaryFixedLayer, &
         BoundaryData, BoundaryFlux, TriDiagDecomp
    use vdiff_lu_solver,     only : fin_vol_lu_decomp
    use vertical_diffusion_solver, only : fin_vol_solve

    ! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
    !                Also, vertical diffusion of cloud droplet number concentration and aerosol number
    !                concentration should be done very carefully in the future version.

    ! Input Arguments
    integer,  intent(in) :: ncol             ! Number of atmospheric columns
    integer,  intent(in) :: pver
    integer,  intent(in) :: pverp
    integer,  intent(in) :: ncnst            ! # of constituents to diffuse. In eddy_diff, only wv. Others, pcnst.
    real(r8), intent(in) :: ztodt            ! delta-t [ s ]

    logical,  intent(in) :: do_diffusion_u_v                ! diffuse horizontal winds [flag]
    logical,  intent(in) :: do_diffusion_s                  ! diffuse dry static energy [flag]
    logical,  intent(in) :: do_diffusion_const(:)           ! diffuse constituents (size ncnst) [flag]
    logical,  intent(in) :: do_molecular_diffusion_const(:) ! molecular diffusion of constituents (size ncnst) [flag]

    logical,  intent(in) :: itaures          ! Whether 'tauresx,tauresy' is updated in this subroutine.

    type(Coords1D), intent(in) :: p          ! Pressure coordinates [ Pa ]
    real(r8), intent(in) :: t(:,:)           ! Temperature [ K ]
    real(r8), intent(in) :: tint(:,:)        ! Temperature at interfaces [ K ]
    real(r8), intent(in) :: rhoi(:,:)        ! Density of air at interfaces [ kg/m3 ]
    real(r8), intent(in) :: taux(:)          ! Surface zonal      stress.
                                             ! Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in) :: tauy(:)          ! Surface meridional stress.
                                             ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in) :: shflx(:)         ! Surface sensible heat flux [ W/m2 ]
    real(r8), intent(in) :: cflx(:,:)        ! Surface constituent flux [ kg/m2/s ] (ncol,ncnst)
    real(r8), intent(in) :: ksrftms(:)       ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
    real(r8), intent(in) :: dragblj(:,:)     ! Drag profile from Beljaars SGO form drag  > 0. [ 1/s ]
    real(r8), intent(in) :: qmincg(:)        ! Minimum constituent mixing ratios from cg fluxes, (ncnst)
    real(r8), intent(in) :: dse_top(:)       ! Dry static energy top boundary condition.
    real(r8), intent(in) :: kvh(:,:)         ! Eddy diffusivity for heat [ m2/s ], interfaces
    real(r8), intent(in) :: kvm(:,:)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ], interfaces
    real(r8), intent(in) :: kvq(:,:)         ! Eddy diffusivity for constituents, interfaces
    real(r8), intent(in) :: cgs(:,:)         ! Counter-gradient star [ cg/flux ], interfaces
    real(r8), intent(in) :: cgh(:,:)         ! Counter-gradient term for heat, interfaces
    real(r8), intent(in) :: rairv(:,:)       ! Composition dependent gas "constant"

    ! Input-Output Arguments
    real(r8), intent(inout) :: u(:,:)        ! U wind. This input is the 'raw' input wind to
                                             ! PBL scheme without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: v(:,:)        ! V wind. This input is the 'raw' input wind to PBL scheme
                                             ! without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: q(:,:,:)      ! Moisture and trace constituents [ kg/kg, #/kg ? ]
    real(r8), intent(inout) :: dse(:,:)      ! Dry static energy [ J/kg ]

    ! Input:  Reserved surface stress at previous time step
    ! Output: Reserved surface stress at current  time step
    real(r8), intent(inout) :: tauresx(:)
    real(r8), intent(inout) :: tauresy(:)

    ! Output Arguments
    real(r8), intent(out)   :: dtk(:,:)      ! T tendency from KE dissipation
    real(r8), intent(out)   :: tautmsx(:)    ! Implicit zonal      turbulent mountain surface stress
                                             ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: tautmsy(:)    ! Implicit meridional turbulent mountain surface stress
                                             ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: topflx(:)     ! Molecular heat flux at the top interface
    character(128), intent(out) :: errmsg    ! Output status

    ! Enable Beljaars drag?
    logical,  intent(in)    :: do_beljaars   ! Flag indicating Beljaars drag

    ! Molecular Diffusion Arguments - mostly optional
    logical,  intent(in)    :: do_molec_diff             ! Flag indicating multiple constituent diffusivities
    logical,  intent(in)    :: use_temperature_molec_diff! Flag indicating that molecular diffusion should apply to temperature, not
                                                         ! dry static energy.
    real(r8), intent(in)    :: cpairv(:,:)      ! Specific heat at constant pressure
    real(r8), intent(in)    :: mbarv(:,:)       ! Composition dependent atmosphere mean molar mass [kg mol-1]

    ! The molecular diffusion module will likely change significantly in
    ! the future, and this module may directly depend on it after that.
    ! Until then, we have these highly specific interfaces hard-coded.

    optional :: vd_lu_qdecomp        ! Constituent-dependent molecular diffusivity routine

    interface
       function vd_lu_qdecomp( &
            ncol , pver, fixed_ubc  , mw     , &
            kv    , kq_scal, mw_facm    , dpidz_sq   , p , &
            interface_boundary, molec_boundary, &
            tint  , ztodt  , nbot_molec , &
            t, m, no_molec_decomp, mbarv) result(decomp)
         import
         integer,  intent(in)    :: ncol
         integer,  intent(in)    :: pver
         integer,  intent(in)    :: nbot_molec
         logical,  intent(in)    :: fixed_ubc
         real(r8), intent(in)    :: kv(:,:) ! interfaces
         real(r8), intent(in)    :: kq_scal(:,:) ! interfaces
         real(r8), intent(in)    :: mw
         real(r8), intent(in)    :: mw_facm(:,:) ! interfaces
         real(r8), intent(in)    :: dpidz_sq(:,:) ! interfaces
         type(Coords1D), intent(in) :: p
         type(BoundaryType), intent(in) :: interface_boundary
         type(BoundaryType), intent(in) :: molec_boundary
         real(r8), intent(in)    :: tint(:,:) ! interfaces
         real(r8), intent(in)    :: ztodt
         real(r8), intent(in)    :: t(:,:)
         integer,  intent(in)    :: m
         type(TriDiagDecomp), intent(in) :: no_molec_decomp
         real(r8), intent(in)    :: mbarv(:,:)
         type(TriDiagDecomp) :: decomp
       end function vd_lu_qdecomp
    end interface

    real(r8), intent(in), optional :: ubc_mmr(:,:) ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(in), optional :: ubc_flux(:,:)      ! Upper boundary flux [ kg/s/m^2 ]

    real(r8), intent(in), optional :: kvt(:,:) ! Kinematic molecular conductivity

    ! FIXME: This input should not be needed (and should not be passed in in vertical_diffusion).
    real(r8), intent(in), optional :: pmid(:,:)

    real(r8), intent(in), optional  :: cnst_mw(:)          ! Molecular weight [ kg/kmole ]
    logical,  intent(in), optional  :: cnst_fixed_ubc(:)   ! Whether upper boundary condition is fixed
    logical,  intent(in), optional  :: cnst_fixed_ubflx(:) ! Whether upper boundary flux is a fixed non-zero value

    integer,  intent(in), optional  :: nbot_molec              ! Bottom level where molecular diffusivity is applied

    ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(in), optional :: kq_scal(:,:)
    ! Local sqrt(1/M_q + 1/M_d) for each constituent
    real(r8), intent(in), optional :: mw_fac(:,:,:)

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i, k, m                                  ! Longitude, level, constituent indices
    logical  :: lqtst(ncol)                             ! Adjust vertical profiles

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

    real(r8) :: tmp1(ncol)                              ! Temporary storage
    real(r8) :: tmpi1(ncol,pverp)                       ! Interface KE dissipation
    real(r8) :: tmpi2(ncol,pverp)                       ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: keg_in(ncol,pver)                       ! KE on entry to subroutine
    real(r8) :: keg_out(ncol,pver)                      ! KE after U and V dissipation/diffusion
    real(r8) :: rrho(ncol)                              ! 1./bottom level density

    real(r8) :: tautotx(ncol)                           ! Total surface stress ( zonal )
    real(r8) :: tautoty(ncol)                           ! Total surface stress ( meridional )

    real(r8) :: dinp_u(ncol,pverp)                      ! Vertical difference at interfaces, input u
    real(r8) :: dinp_v(ncol,pverp)                      ! Vertical difference at interfaces, input v
    real(r8) :: dout_u                                  ! Vertical difference at interfaces, output u
    real(r8) :: dout_v                                  ! Vertical difference at interfaces, output v

    real(r8) :: qtm(ncol,pver)                          ! Temporary copy of q

    real(r8) :: ws(ncol)                                ! Lowest-level wind speed [ m/s ]
    real(r8) :: tau(ncol)                               ! Turbulent surface stress ( not including mountain stress )
    real(r8) :: ksrfturb(ncol)                          ! Surface drag coefficient of 'normal' stress. > 0.
                                                        ! Virtual mass input per unit time per unit area [ kg/s/m2 ]
    real(r8) :: ksrf(ncol)                              ! Surface drag coefficient of 'normal' stress +
                                                        ! Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ]
    real(r8) :: usum_in(ncol)                           ! Vertical integral of input u-momentum. Total zonal
                                                        ! momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    real(r8) :: vsum_in(ncol)                           ! Vertical integral of input v-momentum. Total meridional
                                                        ! momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    real(r8) :: usum_mid(ncol)                          ! Vertical integral of u-momentum after adding explicit residual stress
    real(r8) :: vsum_mid(ncol)                          ! Vertical integral of v-momentum after adding explicit residual stress
    real(r8) :: usum_out(ncol)                          ! Vertical integral of u-momentum after doing implicit diffusion
    real(r8) :: vsum_out(ncol)                          ! Vertical integral of v-momentum after doing implicit diffusion
    real(r8) :: tauimpx(ncol)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: tauimpy(ncol)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: ramda                                   ! dt/timeres [ no unit ]

    real(r8) :: taubljx(ncol)                           ! recomputed explicit/residual beljaars stress
    real(r8) :: taubljy(ncol)                           ! recomputed explicit/residual beljaars stress

    ! Rate at which external (surface) stress damps wind speeds (1/s).
    real(r8) :: tau_damp_rate(ncol, pver)

    ! Combined molecular and eddy diffusion.
    real(r8) :: kv_total(ncol,pverp)

    !--------------------------------
    ! Variables needed for WACCM-X
    !--------------------------------
    real(r8) :: ttemp(ncol,pver)             ! temporary temperature array
    real(r8) :: ttemp0(ncol,pver)            ! temporary temperature array

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

    errmsg = ''
    if( do_diffusion_u_v .and. (.not. do_diffusion_s) ) then
        errmsg = 'compute_vdiff: must diffuse s if diffusing horizontal winds'
        return
    end if

    ! Check if logical array is of the expected size
    if ( size(do_diffusion_const) .ne. ncnst ) then
        write(errmsg,*) 'compute_vdiff: do_diffusion_const size ', size(do_diffusion_const), ' is not equal to ncnst ', ncnst
        return
    endif

    if ( size(do_molecular_diffusion_const) .ne. ncnst ) then
        write(errmsg,*) 'compute_vdiff: do_molecular_diffusion_const size ', size(do_molecular_diffusion_const), ' is not equal to ncnst ', ncnst
        return
    endif

    !--------------------------------------- !
    ! Computation of Molecular Diffusivities !
    !--------------------------------------- !

    ! Modification : Why 'kvq' is not changed by molecular diffusion ?

    if( do_molec_diff ) then

        if( (.not.present(vd_lu_qdecomp)) .or. (.not.present(kvt)) &
             .or. (.not. present(ubc_mmr)) .or. (.not. present(ubc_flux)) ) then
              errmsg = 'compute_vdiff: do_molec_diff true but vd_lu_qdecomp or kvt missing'
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

    ! FIXME: The following five lines are kept in only to preserve answers;
    !        they really should be taken out completely.
    if (do_molec_diff) then
      tmpi2(:ncol,1) = ztodt * (gravit * rhoi(:ncol,1))**2 / ( pmid(:ncol,1) - p%ifc(:,1) )
    endif
    dpidz_sq(:,1) = gravit*(p%ifc(:,1) / (rairv(:ncol,1)*t(:ncol,1)))
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

    if( do_diffusion_u_v ) then

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

       v(:ncol,:) = fin_vol_solve(ztodt, p, v(:ncol,:), ncol, pver, &
                         coef_q=tau_damp_rate,                      &
                         coef_q_diff=kvm(:ncol,:)*dpidz_sq(:ncol,:))

       u(:ncol,:) = fin_vol_solve(ztodt, p, u(:ncol,:), ncol, pver, &
                         coef_q=tau_damp_rate,                      &
                         coef_q_diff=kvm(:ncol,:)*dpidz_sq(:ncol,:))



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

              if( itaures ) then
                 tauresx(i) = taux(i) + tauresx(i) - tauimpx(i) + tautmsx(i) + taubljx(i)
                 tauresy(i) = tauy(i) + tauresy(i) - tauimpy(i) + tautmsy(i) + taubljy(i)
              endif

          else

             tautotx(i) = taux(i) + tautmsx(i)
             tautoty(i) = tauy(i) + tautmsy(i)
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

    end if ! End of diffuse horizontal momentum routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

    ! Modification : In future, we should diffuse the fully conservative
    !                moist static energy,not the dry static energy.

    if( do_diffusion_s ) then

       ! Add counter-gradient to input static energy profiles
       do k = 1, pver
          dse(:ncol,k) = dse(:ncol,k) + ztodt * p%rdel(:,k) * gravit  *      &
                         ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgh(:ncol,k+1) &
                         - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgh(:ncol,k  ) )
       end do

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

          dse(:ncol,:) = fin_vol_solve(ztodt, p, dse(:ncol,:), ncol, pver, &
                         coef_q_diff=kvh(:ncol,:)*dpidz_sq(:ncol,:),       &
                         upper_bndry=interface_boundary,                   &
                         l_cond=BoundaryData(dse_top(:ncol)))

          ! Calculate flux at top interface

          ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

          topflx(:ncol) =  - kvh(:ncol,1) * tmpi2(:ncol,1) / (ztodt*gravit) * &
               ( dse(:ncol,1) - dse_top(:ncol) )

          ttemp0 = t(:ncol,:)
          ttemp = ttemp0

          ! upper boundary is zero flux for extended model
          ttemp = fin_vol_solve(ztodt, p, ttemp, ncol, pver, &
                  coef_q_diff=kvt(:ncol,:)*dpidz_sq(:ncol,:),&
                  coef_q_weight=cpairv(:ncol,:))


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
          dse(:ncol,:) = fin_vol_solve(ztodt, p, dse(:ncol,:), ncol, pver, &
                         coef_q_diff=kv_total(:ncol,:)*dpidz_sq(:ncol,:),  &
                         upper_bndry=interface_boundary,                   &
                         l_cond=BoundaryData(dse_top(:ncol)))

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
         coef_q_diff=kvq(:ncol,:)*dpidz_sq(:ncol,:))

    const_diffuse_loop: do m = 1, ncnst
       if( .not. do_diffusion_const(m) ) then
         cycle const_diffuse_loop
       endif

       ! Add the nonlocal transport terms to constituents in the PBL.
       ! Check for neg q's in each constituent and put the original vertical
       ! profile back if a neg value is found. A neg value implies that the
       ! quasi-equilibrium conditions assumed for the countergradient term are
       ! strongly violated.
       qtm(:ncol,:pver) = q(:ncol,:pver,m)

       do k = 1, pver
          q(:ncol,k,m) = q(:ncol,k,m) + &
                         ztodt * p%rdel(:,k) * gravit  * ( cflx(:ncol,m) * rrho(:ncol) ) * &
                         ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgs(:ncol,k+1)               &
                         - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgs(:ncol,k  ) )
       end do
       lqtst(:ncol) = all(q(:ncol,1:pver,m) >= qmincg(m), 2)
       do k = 1, pver
          q(:ncol,k,m) = merge( q(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol) )
       end do

       ! Add the explicit surface fluxes to the lowest layer
       q(:ncol,pver,m) = q(:ncol,pver,m) + tmp1(:ncol) * cflx(:ncol,m)

       if( do_molec_diff .and. do_molecular_diffusion_const(m) ) then
         ! do molecular diffusion

         ! This is for solving molecular diffusion of minor species, thus, for WACCM-X, bypass O and O2 (major species)
         ! Major species diffusion is calculated separately.  -Hanli Liu
         decomp = vd_lu_qdecomp( &
           ncol = ncol, &
           pver = pver, &
           m = m, & ! constituent index
           fixed_ubc = cnst_fixed_ubc(m), &
           mw = cnst_mw(m), &
           kv = kvq(:ncol,:), & ! eddy diffusivity at interfaces
           kq_scal = kq_scal(:ncol,:), & ! molecular diffusivity at interfaces
           mw_facm = mw_fac(:ncol,:,m), & ! composition dependent sqrt(1/M_q + 1/M_d), at interfaces, this constituent
           dpidz_sq = dpidz_sq(:ncol,:), & ! square of vertical derivative of pint
           p = p_molec, &
           interface_boundary = interface_boundary, &
           molec_boundary = molec_boundary, &
           t = t(:ncol,:), & ! temperature, midpoints
           tint = tint(:ncol,:), & ! temperature, interfaces
           mbarv = mbarv(:ncol,:), & ! composition dependent atmosphere mean molar mass
           ztodt = ztodt, &
           nbot_molec = nbot_molec, &
           no_molec_decomp = no_molec_decomp)

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
         ! not doing molecular diffusion
         if (present(cnst_fixed_ubc)) then
            ! explicitly set mmr in top layer for cases where molecular diffusion is not active
            if (cnst_fixed_ubc(m)) then
               q(:ncol,1,m) = ubc_mmr(:ncol,m)
            endif
         end if
         call no_molec_decomp%left_div(q(:ncol,:,m))
       end if
    end do const_diffuse_loop

    call no_molec_decomp%finalize()

  end subroutine compute_vdiff

end module diffusion_solver_cam
