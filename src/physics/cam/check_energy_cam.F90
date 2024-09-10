! Shim for CCPP-ized check_energy routines in CAM
module check_energy_cam
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver
  use constituents,    only: pcnst
  use physics_types,   only: physics_state, physics_tend

  implicit none
  private

  public :: check_energy_cam_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_cam_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes

contains
  ! Compute initial values of energy and water integrals,
  ! zero cumulative tendencies
  subroutine check_energy_cam_timestep_init(state, tend, pbuf, col_type)
    use physics_buffer,  only: physics_buffer_desc, pbuf_set_field
    use cam_abortutils,  only: endrun
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height, vc_dry_pressure
    use physics_types,   only: phys_te_idx, dyn_te_idx
    use time_manager,    only: is_first_step
    use physconst,       only: cpair, rair
    use air_composition, only: cpairv, cp_or_cv_dycore

    ! To remove, pbuf indices in check_energy
    use check_energy,    only: teout_idx

    ! CCPP-ized subroutine
    use check_energy_chng, only: check_energy_chng_timestep_init

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag indicating whether using grid or subcolumns

    real(r8)  :: local_cp_phys(state%psetcols,pver)
    real(r8)  :: local_cp_or_cv_dycore(state%psetcols,pver)
    integer   :: lchnk     ! chunk identifier
    integer   :: ncol      ! number of atmospheric columns
    character(len=512) :: errmsg
    integer            :: errflg

    lchnk = state%lchnk
    ncol  = state%ncol

    ! The code below is split into not-subcolumns and subcolumns code, as there is different handling of the
    ! cp passed into the hydrostatic energy call. CAM-SIMA does not support subcolumns, so we keep this special
    ! handling inside this shim module. (hplin, 9/9/24)
    if(state%psetcols == pcols) then
        ! No subcolumns
        local_cp_phys(:,:) = cpairv(:,:,lchnk)
        local_cp_or_cv_dycore(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)
    else if (state%psetcols > pcols) then
        ! Subcolumns code
        ! Subcolumns specific error handling
        if(.not. all(cpairv(:,:,lchnk) == cpair)) then
            call endrun('check_energy_timestep_init: cpairv is not allowed to vary when subcolumns are turned on')
        endif

        local_cp_phys(1:ncol,:) = cpair

        if (vc_dycore == vc_height) then
            ! MPAS specific hydrostatic energy computation (internal energy)
            local_cp_or_cv_dycore(:ncol,:) = cpair-rair
        else if(vc_dycore == vc_dry_pressure) then
            ! SE specific hydrostatic energy (enthalpy)
            local_cp_or_cv_dycore(:ncol,:) = cpair
        else
            ! cp_or_cv is not used in the underlying subroutine, zero it out to be sure
            local_cp_or_cv_dycore(:ncol,:) = 0.0_r8
        endif
    end if

    ! Call CCPP-ized underlying subroutine.
    call check_energy_chng_timestep_init( &
        ncol            = ncol, &
        pver            = pver, &
        pcnst           = pcnst, &
        q               = state%q(1:ncol,1:pver,1:pcnst), &
        pdel            = state%pdel(1:ncol,1:pver), &
        u               = state%u(1:ncol,1:pver), &
        v               = state%v(1:ncol,1:pver), &
        T               = state%T(1:ncol,1:pver), &
        pintdry         = state%pintdry(1:ncol,1:pver), &
        phis            = state%phis(1:ncol), &
        zm              = state%zm(1:ncol,:), &
        temp_ini        = state%temp_ini(:ncol,:), &
        z_ini           = state%z_ini(:ncol,:), &
        cp_phys         = local_cp_phys(1:ncol,:), &
        cp_or_cv_dycore = local_cp_or_cv_dycore(1:ncol,:), &
        te_ini_phys     = state%te_ini(1:ncol,phys_te_idx), &
        te_ini_dyn      = state%te_ini(1:ncol,dyn_te_idx),  &
        tw_ini          = state%tw_ini(1:ncol),             &
        te_cur_phys     = state%te_cur(1:ncol,phys_te_idx), &
        te_cur_dyn      = state%te_cur(1:ncol,dyn_te_idx),  &
        tw_cur          = state%tw_cur(1:ncol),             &
        tend_te_tnd     = tend%te_tnd(1:ncol),              &
        tend_tw_tnd     = tend%tw_tnd(1:ncol),              &
        count           = state%count,                      &
        vc_physics      = vc_physics, &
        vc_dycore       = vc_dycore, &
        errmsg          = errmsg, &
        errflg          = errflg  &
    )

    ! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini(:,dyn_te_idx), col_type=col_type)
    end if

  end subroutine check_energy_cam_timestep_init

  ! Check that the energy and water change matches the boundary fluxes
  subroutine check_energy_cam_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)
    use cam_thermo,      only: get_hydrostatic_energy
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height, vc_dry_pressure
    use cam_abortutils,  only: endrun
    use physics_types,   only: phys_te_idx, dyn_te_idx
    use physconst,       only: cpair, rair
    use air_composition, only: cpairv, cp_or_cv_dycore

    ! CCPP-ized subroutine
    use check_energy_chng, only: check_energy_chng_run

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in) :: nstep                  ! current timestep number
    real(r8), intent(in) :: ztodt                  ! 2 delta t (model time increment)
    real(r8), intent(in) :: flx_vap(:)             ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in) :: flx_cnd(:)             ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in) :: flx_ice(:)             ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in) :: flx_sen(:)             ! (pcols) -boundary flux of sensible heat (w/m2)

    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i                                   ! column index
    real(r8)  :: local_cp_phys(state%psetcols,pver)
    real(r8)  :: local_cp_or_cv_dycore(state%psetcols,pver)
    real(r8)  :: scaling_dycore(state%ncol,pver)
    character(len=512) :: errmsg
    integer            :: errflg

    lchnk = state%lchnk
    ncol  = state%ncol

    if(state%psetcols == pcols) then
        ! No subcolumns
        local_cp_phys(:,:) = cpairv(:,:,lchnk)
        local_cp_or_cv_dycore(:ncol,:) = cp_or_cv_dycore(:ncol,:,lchnk)

        scaling_dycore(:ncol,:)  = cpairv(:,:,lchnk)/local_cp_or_cv_dycore(:ncol,:) ! cp/cv scaling
    elseif(state%psetcols > pcols) then
        ! Subcolumns
        if(.not. all(cpairv(:,:,:) == cpair)) then
            call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
        endif

        local_cp_phys(:,:) = cpair

        ! Note: cp_or_cv set above for pressure coordinate
        if (vc_dycore == vc_height) then
            ! compute cv if vertical coordinate is height: cv = cp - R
            local_cp_or_cv_dycore(:ncol,:) = cpair-rair
            scaling_dycore(:ncol,:)  = cpairv(:,:,lchnk)/local_cp_or_cv_dycore(:ncol,:) ! cp/cv scaling
        else if (vc_dycore == vc_dry_pressure) then
            ! SE specific hydrostatic energy
            local_cp_or_cv_dycore(:ncol,:) = cpair
            scaling_dycore(:ncol,:) = 1.0_r8
        else
            ! Moist pressure... use phys formula
            local_cp_or_cv_dycore(:ncol,:) = local_cp_phys(:ncol,:,lchnk)
            scaling_dycore(:ncol,:)  = cpairv(:,:,lchnk)/local_cp_or_cv_dycore(:ncol,:) ! cp/cv scaling
        end if
    endif

    ! Call CCPP-ized underlying subroutine.
    call check_energy_chng_run( &
        ncol            = ncol, &
        pver            = pver, &
        pcnst           = pcnst, &
        q               = state%q(1:ncol,1:pver,1:pcnst), &
        pdel            = state%pdel(1:ncol,1:pver), &
        u               = state%u(1:ncol,1:pver), &
        v               = state%v(1:ncol,1:pver), &
        T               = state%T(1:ncol,1:pver), &
        pintdry         = state%pintdry(1:ncol,1:pver), &
        phis            = state%phis(1:ncol), &
        zm              = state%zm(1:ncol,:), &
        cp_phys         = local_cp_phys(1:ncol,:), &
        cp_or_cv_dycore = local_cp_or_cv_dycore(1:ncol,:),  &
        scaling_dycore  = scaling_dycore(1:ncol,:),         &
        te_cur_phys     = state%te_cur(1:ncol,phys_te_idx), &
        te_cur_dyn      = state%te_cur(1:ncol,dyn_te_idx),  &
        tw_cur          = state%tw_cur(1:ncol),             &
        tend_te_tnd     = tend%te_tnd(1:ncol),              &
        tend_tw_tnd     = tend%tw_tnd(1:ncol),              &
        temp_ini        = state%temp_ini(:ncol,:),          &
        z_ini           = state%z_ini(:ncol,:),             &
        count           = state%count,                      &
        ztodt           = ztodt,                            &
        vc_physics      = vc_physics, &
        vc_dycore       = vc_dycore,  &
        errmsg          = errmsg, &
        errflg          = errflg  &
    )

  end subroutine check_energy_cam_chng
end module check_energy_cam