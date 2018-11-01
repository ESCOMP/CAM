module TJ2016_cam
  !-----------------------------------------------------------------------
  !
  ! Purpose: Implement idealized moist Held-Suarez forcings
  !          Thatcher, D. R. and C. Jablonowski (2016)
  !          "A moist aquaplanet variant of the Held-Suarez test
  !           for atmospheric model dynamical cores",
  !          Geosci. Model Dev., Vol. 9, 1263-1292,
  !          doi:10.5194/gmd-9-1263-2016, 2016.
  !
  !-----------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use constituents,   only: pcnst

  use physics_buffer, only: dtype_r8, pbuf_add_field, physics_buffer_desc, &
                            pbuf_set_field, pbuf_get_field
  use camsrfexch,     only: cam_in_t

  use cam_history,    only: outfld
  use time_manager,   only: is_first_step

  implicit none
  private
  save

  public :: Thatcher_Jablonowski_register
  public :: Thatcher_Jablonowski_init
  public :: Thatcher_Jablonowski_precip_tend
  public :: Thatcher_Jablonowski_sfc_pbl_hs_tend

  integer :: prec_pcw_idx  = 0
  integer :: relhum_idx    = 0

!=======================================================================
CONTAINS
!=======================================================================

subroutine Thatcher_Jablonowski_register()

   call pbuf_add_field('PREC_PCW',  'physpkg', dtype_r8, (/pcols/),      prec_pcw_idx)
   call pbuf_add_field('RELHUM',    'physpkg', dtype_r8, (/pcols,pver/), relhum_idx)

end subroutine Thatcher_Jablonowski_register

!========================================================================================

  subroutine Thatcher_Jablonowski_init(pbuf2d)
    use cam_history,    only: addfld, add_default
    use physconst,      only: gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir
    use hycoef,         only: ps0, etamid
    use tj2016,         only: Thatcher_Jablonowski_set_const

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call Thatcher_Jablonowski_set_const(gravit, cappa, rair, cpair, latvap, rh2o, epsilo, rhoh2o, zvir, ps0, etamid)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
    call add_default('QRS', 1, ' ')
    ! These fields are added by vertical diffusion when full physics is used
    call addfld('KVH'        , (/ 'ilev' /), 'A', 'm2/s'      , 'Vertical diffusion diffusivities (heat/moisture)' )
    call add_default('KVH', 1, ' ')
    call addfld('KVM'        , (/ 'ilev' /), 'A', 'm2/s'      , 'Vertical diffusion diffusivities (momentum)'      )
    call add_default('KVM', 1, ' ')
    call addfld('DTV'         , (/ 'lev' /), 'A', 'K/s'       , 'T vertical diffusion'                             )
    call add_default('DTV', 1, ' ')
    call addfld('DUV'         , (/ 'lev' /), 'A', 'm/s2'      , 'U vertical diffusion'                             )
    call add_default('DUV', 1, ' ')
    call addfld('DVV'         , (/ 'lev' /), 'A', 'm/s2'      , 'V vertical diffusion'                             )
    call add_default('DVV', 1, ' ')
    call addfld ('VD01',        (/ 'lev' /), 'A', 'kg/kg/s'   , 'Q tendency (vertical diffusion)'                  )
    call add_default('VD01', 1, ' ')

    if (is_first_step()) then
       call pbuf_set_field(pbuf2d, prec_pcw_idx, 0._r8)
    end if

  end subroutine Thatcher_Jablonowski_init

!========================================================================================

  subroutine Thatcher_Jablonowski_precip_tend(state, ptend, ztodt, pbuf)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Run the precipitation process of the Thatcher-Jablonowski physics (see tj2016.F90)
    !
    !-----------------------------------------------------------------------
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use physconst,          only: cpair
    use TJ2016,             only: Thatcher_Jablonowski_precip

    ! arguments

    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt                         ! physics timestep

    type(physics_ptend), intent(out)   :: ptend                         ! Package tendencies

    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local variables

    real(r8), pointer :: prec_pcw(:)
    real(r8), pointer :: relhum(:,:)

    !---------------------------Local workspace-----------------------------
    !
    integer  :: lchnk                         ! chunk identifier
    integer  :: ncol                          ! number of atmospheric columns

    real(r8) :: T(state%ncol, pver)           ! T temporary
    real(r8) :: qv(state%ncol, pver)          ! Q temporary
    logical  :: lq(pcnst)                     ! Calc tendencies?
                                              ! output from parameterization
    real(r8) :: precc(state%ncol)             ! convective precip

    integer  :: k
    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! Gather temporary arrays
    T(:ncol, :)  = state%T(:ncol, :)
    qv(:ncol, :) = state%Q(:ncol, :, 1)

    ! initialize individual parameterization tendencies
    lq           = .false.
    lq(1)        = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Thatcher-Jablonowski precip',   &
         ls=.true., lu=.true., lv=.true., lq=lq)

    ! Thatcher_Jablonowski_precip physics arguments

    ! Input arguments
    ! ncol:   Number of columns
    ! pver:   Number of vertical levels
    ! ztodt:  Time step (s) (in)
    ! pmid:   Mid-point pressure (Pa)
    ! pdel:   Layer thickness (Pa)

    ! Input/Output arguments
    ! T:      Temperature (K)
    ! qv:     Specific humidity (kg/kg)

    ! Output arguments
    ! relhum: relative humidity (%)
    ! precl:  large-scale precipitation rate (m/s)
    ! precc:  convective precipitation rate (m/s) (optional process)

    call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
    call pbuf_get_field(pbuf, relhum_idx,   relhum)

    call Thatcher_Jablonowski_precip(ncol, pver, ztodt,                     &
         state%pmid(:ncol,:), state%pdel(:ncol,:),                          &
         T, qv, relhum(:ncol,:), prec_pcw(:ncol), precc)

    ! Back out temperature and specific humidity tendencies from updated fields
    do k = 1, pver
      ptend%s(:ncol,k)   = (T(:, k) - state%T(:ncol, k)) / ztodt * cpair
      ptend%q(:ncol,k,1) = (qv(:, k) - state%q(:ncol, k, 1)) / ztodt
    end do

 end subroutine Thatcher_Jablonowski_precip_tend

!========================================================================================

 subroutine Thatcher_Jablonowski_sfc_pbl_hs_tend(state, ptend, ztodt, cam_in)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Run the surface flux and PBL processes of the Thatcher-Jablonowski physics (moist Held-Suarez)
    !
    !-----------------------------------------------------------------------
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use physconst,          only: cpair
    use phys_grid,          only: get_rlat_all_p
    use TJ2016,             only: Thatcher_Jablonowski_sfc_pbl_hs

    ! Arguments
    type(physics_state), intent(in)    :: state
    real(r8),            intent(in)    :: ztodt
    type(physics_ptend), intent(out)   :: ptend
    type(cam_in_t),     intent(inout)  :: cam_in

    !---------------------------Local workspace-----------------------------

    integer  :: k                             ! loop index
    integer  :: lchnk                         ! chunk identifier
    integer  :: ncol                          ! number of atmospheric columns

    real(r8) :: clat(state%ncol)              ! latitudes(radians) for columns
    real(r8) :: lnpint(state%ncol, 2)         ! ln(int. press. (Pa))
    real(r8) :: T(state%ncol, pver)           ! T temporary
    real(r8) :: qv(state%ncol, pver)          ! Q temporary (specific humidity)
    real(r8) :: U(state%ncol, pver)           ! U temporary
    real(r8) :: V(state%ncol, pver)           ! V temporary
    logical  :: lq(pcnst)                     ! Calc tendencies?

    ! output from parameterization
    real(r8) :: dqdt_vdiff(state%ncol,pver)   ! Q tendency due to vertical PBL diffusion in kg/kg/s
    real(r8) :: dtdt_vdiff(state%ncol,pver)   ! T tendency due to vertical PBL diffusion in K/s
    real(r8) :: dtdt_heating(state%ncol,pver) ! temperature tendency from relaxation in K/s
    real(r8) :: Km(state%ncol,pver+1)         ! Eddy diffusivity at layer interfaces for boundary layer calculations (m2/s)
    real(r8) :: Ke(state%ncol,pver+1)         ! Eddy diffusivity at layer interfaces for boundary layer calculations (m2/s)
    !-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call get_rlat_all_p(lchnk, ncol, clat)

    ! Gather temporary arrays
    lnpint(:ncol, 1:2) = state%lnpint(:ncol,pver:pver+1)
    T(:ncol, :) = state%T(:ncol, :)
    U(:ncol, :) = state%U(:ncol, :)
    V(:ncol, :) = state%V(:ncol, :)
    qv(:ncol, :) = state%Q(:ncol, :, 1)

    ! initialize individual parameterization tendencies
    lq           = .false.
    lq(1)        = .true.
    call physics_ptend_init(ptend, state%psetcols, 'Thatcher-Jablonowski sfc_pbl_hs_tend',   &
         ls=.true., lu=.true., lv=.true., lq=lq)

    ! Thatcher_Jablonowski_sfc_pbl_hs physics arguments

    ! Input arguments
    ! ncol:   Number of columns
    ! pver:   Number of vertical levels
    ! ztodt:  Time step (s) (in)
    ! clat:   Column latitude
    ! PS:     Surface pressure
    ! pmid:   Mid-point pressure (Pa)
    ! pint:   Interface pressure (Pa)
    ! lnpint: ln(interface pressure (Pa)) at surface
    ! rpdel:  Reciprocal of layer thickness (Pa)

    ! Input/Output arguments
    ! T:      Temperature
    ! U:      Zonal wind
    ! V:      Meridional wind
    ! qv:     Specific humidity (moist mixing ratio)

    ! Output arguments
    ! cam_in%shf:   Surface sensible heat flux
    ! cam_in%lhf:   Surface latent heat flux
    ! cam_in%wsx:   Surface momentum flux in the zonal direction
    ! cam_in%wsy:   Surface momentum flux in the meridional direction
    ! cam_in%cflx:  surface water flux in kg/m2/s
    ! dqdt_vdiff:   Q tendency due to vertical diffusion (PBL)
    ! dtdt_vdiff:   T tendency due to vertical diffusion (PBL)
    ! dtdt_heating: Temperature tendency in K/s from relaxation
    ! Km:           Eddy diffusivity for boundary layer calculations
    ! Ke:           Eddy diffusivity for boundary layer calculations
    ! cam_in%sst:   Sea surface temperature K (varied by latitude)

    call Thatcher_Jablonowski_sfc_pbl_hs(ncol, pver, ztodt, clat,                 &
         state%ps(:ncol), state%pmid(:ncol,:), state%pint(:ncol,:), lnpint,       &
         state%rpdel(:ncol,:), T, U, V, qv, cam_in%shf(:ncol), cam_in%lhf(:ncol), &
         cam_in%wsx(:ncol), cam_in%wsy(:ncol), cam_in%cflx(:ncol,1), dqdt_vdiff,  &
         dtdt_vdiff, dtdt_heating, Km, Ke, cam_in%sst(:ncol))

    ! Back out tendencies from updated fields
    do k = 1, pver
      ptend%s(:ncol,k) = (T(:, k) - state%T(:ncol, k)) / ztodt * cpair
      ptend%u(:ncol,k) = (U(:, k) - state%U(:ncol, k)) / ztodt
      ptend%v(:ncol,k) = (V(:, k) - state%V(:ncol, k)) / ztodt
      ptend%q(:ncol,k,1) = (qv(:, k) - state%q(:ncol, k, 1)) / ztodt
    end do

    !===============================================================================
    ! Archive diagnostic fields
    !===============================================================================
    call outfld('QRS',  dtdt_heating, ncol,  lchnk) ! T tendency from temperature relaxation (mimics radiation, K/s)
    call outfld('KVH',  Ke,           ncol,  lchnk) ! Eddy diffusivity (heat and moisture, m2/s)
    call outfld('KVM',  Km,           ncol,  lchnk) ! Eddy diffusivity (momentum, m2/s)
    call outfld('DUV',  ptend%u,      pcols, lchnk) ! PBL u tendency (m/s2)
    call outfld('DVV',  ptend%v,      pcols, lchnk) ! PBL v tendency (m/s2)
    call outfld('DTV',  dtdt_vdiff,   ncol,  lchnk) ! PBL + surface flux T tendency (K/s)
    call outfld('VD01', dqdt_vdiff,   ncol,  lchnk) ! PBL + surface flux Q tendency (kg/kg/s)

 end subroutine Thatcher_Jablonowski_sfc_pbl_hs_tend

end module TJ2016_cam
