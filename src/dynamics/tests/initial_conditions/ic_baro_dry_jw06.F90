module ic_baro_dry_jw06
  !-----------------------------------------------------------------------
  !
  ! Purpose: Set idealized initial conditions for the Jablonowski and
  !          Williamson baroclinic instability test.
  !          References: 
  !          Jablonowski, C., and D. L. Williamson (2006), A Baroclinic Instability Test Case for 
  !              Atmospheric Model Dynamical Cores, Quart. J. Roy. Met. Soc., Vol. 132, 2943-2975
  !          Jablonowski, C., and D. L. Williamson (2006), A Baroclinic Wave Test Case for Dynamical 
  !              Cores of General Circulation Models: Model Intercomparisons, 
  !              NCAR Technical Note NCAR/TN-469+STR, Boulder, CO, 89 pp.
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  use physconst, only : rair, cpair, gravit, rearth, pi, omega
  use hycoef,    only : hyai, hybi, hyam, hybm, ps0

  implicit none
  private

  !=======================================================================
  !    JW06 Dry baroclinic wave test case parameters
  !=======================================================================
  real(r8), parameter, private ::             &
       eta_tropo              = 0.2_r8,       & ! tropopause level (hybrid vertical coordinate))
       u0                     = 35._r8,       & ! maximum jet speed 35 m/s
       T0                     = 288._r8,      & ! horizontal mean T at the surface
       p00                    = 1.e5_r8,      & ! surface pressure in Pa
       eta0                   = 0.252_r8,     & ! center of jets (hybrid vertical coordinate)
       radius                 = 10._r8,       & ! reciprocal radius of the perturbation without the Earth's radius 'a'
       perturbation_amplitude = 1._r8,        & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20._r8,       & ! longitudinal position, 20E
       perturbation_latitude  = 40._r8,       & ! latitudinal position, 40N
       eta_sfc                = 1._r8,        & ! hybrid value at the surface
       delta_T                = 480000._r8,   & ! in K, parameter for T mean calculation
       gamma                  = 0.005_r8        ! lapse rate in K/m
  real(r8) :: a_omega, exponent

  real(r8), parameter :: deg2rad = pi/180._r8   ! conversion to radians

  ! Public interface
  public :: bc_dry_jw06_set_ic

contains

  subroutine bc_dry_jw06_set_ic(vcoord, latvals, lonvals, U, V, T, PS, PHIS, &
                                Q, m_cnst, mask, verbose)
    use dyn_tests_utils, only: vc_moist_pressure, vc_dry_pressure, vc_height
    use constituents,    only: cnst_name
    use const_init,      only: cnst_init_default

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set baroclinic wave initial values for dynamics state variables
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer, intent(in)               :: vcoord
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
                                                    ! z_k for vccord 1)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(inout) :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use
    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    logical                           :: lu,lv,lt,lq,l3d_vars
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'BC_DRY_JW06_SET_IC'
    real(r8)                          :: tmp
    real(r8)                          :: r(size(latvals))
    real(r8)                          :: eta
    real(r8)                          :: factor
    real(r8)                          :: perturb_lon, perturb_lat
    real(r8)                          :: phi_vertical
    real(r8)                          :: u_wind(size(latvals))

       a_omega                = rearth*omega
       exponent               = rair*gamma/gravit

    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun(subname//': input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1

    !
    ! We do not yet handle height-based vertical coordinates
    if (vcoord == vc_height) then
      call endrun(subname//':  height-based vertical coordinate not currently supported')
    end if

    !
    !*******************************
    !
    ! initialize surface pressure
    !
    !*******************************
    !
    if (present(PS)) then
      where(mask_use)
        PS = p00
      end where

      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
    end if
    !
    !*******************************
    !
    ! Initialize PHIS
    !
    !*******************************
    !
    if (present(PHIS)) then
      tmp =  u0 * (cos((eta_sfc-eta0)*pi*0.5_r8))**1.5_r8
      where(mask_use)
        PHIS(:) = ((-2._r8*(sin(latvals(:)))**6 * ((cos(latvals(:)))**2 + 1._r8/3._r8) + 10._r8/63._r8)*tmp   &
                + (8._r8/5._r8*(cos(latvals(:)))**3 * ((sin(latvals(:)))**2 + 2._r8/3._r8) - pi/4._r8)*a_omega)*tmp
      end where
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
    end if
    !
    !*******************************
    !
    ! Initialize 3D vars
    !
    !
    !*******************************
    !
    lu = present(U)
    lv = present(V)
    lT = present(T)
    lq = present(Q)
    l3d_vars = lu .or. lv .or. lt .or.lq
    nlev = -1
    if (l3d_vars) then
      if (lu) nlev = size(U, 2)
      if (lv) nlev = size(V, 2)
      if (lt) nlev = size(T, 2)
      if (lq) nlev = size(Q, 2)

      if (lu) then
        do k = 1, nlev
          perturb_lon = perturbation_longitude * deg2rad
          perturb_lat = perturbation_latitude  * deg2rad
          phi_vertical = ((hyam(k)+hybm(k)) - eta0) *0.5_r8*pi
             where(mask_use)
               ! background wind
               u_wind(:) = (cos(phi_vertical))**1.5_r8 * 4._r8 * u0 * (sin(latvals(:)))**2 * (cos(latvals(:)))**2
               ! great circle distance without radius 'a'
               r(:)      = acos( sin(perturb_lat)*sin(latvals(:)) + cos(perturb_lat)*cos(latvals(:))*cos(lonvals(:)-perturb_lon))
               !  background + perturbation wind
               U(:,k)    = perturbation_amplitude*exp(- (r(:)*radius)**2 ) + u_wind(:)
             end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          U initialized by "',subname,'"'
        end if
      end if
      if (lv) then
        do k = 1, nlev
          where(mask_use)
            V(:,k) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          V initialized by "',subname,'"'
        end if
      end if
      if (lt) then
       do k = 1, nlev
         eta = hyam(k) + hybm(k)
            ! background temperature
            if (eta .ge. eta_tropo) then
              tmp = T0*eta**exponent
            else
              tmp = T0*eta**exponent + delta_T*(eta_tropo-eta)**5
            endif
            factor       = eta*pi*u0/rair
            phi_vertical = (eta - eta0) * 0.5_r8*pi
           where(mask_use)
            ! background temperature 'tmp' plus temperature deviation
            T(:,k)       = factor * 1.5_r8 * sin(phi_vertical) * (cos(phi_vertical))**0.5_r8 *                &
                         ((-2._r8*(sin(latvals(:)))**6 * ((cos(latvals(:)))**2 + 1._r8/3._r8) + 10._r8/63._r8)*              &
                         u0 * (cos(phi_vertical))**1.5_r8  + &
                         (8._r8/5._r8*(cos(latvals(:)))**3 * ((sin(latvals(:)))**2 + 2._r8/3._r8) - pi/4._r8)*a_omega*0.5_r8 ) + &
                         tmp
           end where
       enddo
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          T initialized by "',subname,'"'
        end if
      end if
      if (lq) then
        do k = 1, nlev
          where(mask_use)
            Q(:,k,1) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '         ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
        end if
      end if
    end if

    if (lq) then
      ncnst = size(m_cnst, 1)
      if ((vcoord == vc_moist_pressure) .or. (vcoord == vc_dry_pressure)) then
        do m = 2, ncnst
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end do
      end if
    end if

    deallocate(mask_use)

  end subroutine bc_dry_jw06_set_ic

end module ic_baro_dry_jw06
