module kessler_cam

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc, pbuf_init_time, dtype_r8, &
                            pbuf_add_field, pbuf_get_field

  implicit none
  private
  save

  public :: kessler_register, kessler_cam_init, kessler_tend

  integer :: ixcldliq = -1 ! cloud liquid mixing ratio index
  integer :: ixrain   = -1 ! rain liquid mixing ratio index

  ! physics buffer indices
  integer :: prec_sed_idx  = 0
  integer :: relhum_idx    = 0

!========================================================================================
contains
!========================================================================================

  subroutine kessler_register()
    use physconst,      only: cpair, mwh2o
    use constituents,   only: cnst_add

    ! Add liquid constituents
    call cnst_add('CLDLIQ', mwh2o, cpair, 0._r8, ixcldliq,                    &
         longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
    call cnst_add('RAINQM', mwh2o, cpair, 0._r8, ixrain,                      &
         longname='Grid box averaged rain water amount', is_convtran1=.true.)

    call pbuf_add_field('PREC_SED', 'physpkg', dtype_r8, (/pcols/), prec_sed_idx)
    call pbuf_add_field('RELHUM',   'physpkg', dtype_r8, (/pcols,pver/), relhum_idx)

  end subroutine kessler_register

!========================================================================================

  subroutine kessler_cam_init()

    use physconst,      only: latvap, rhoh2o
    use ref_pres,       only: psurf_ref
    use constituents,   only: cnst_name, cnst_longname, bpcnst, apcnst
    use cam_history,    only: addfld, add_default
    use cam_abortutils, only: endrun

    use kessler,        only: kessler_init

    !
    !---------------------------Local workspace-----------------------------
    !

    integer                            :: errflg
    character(len=512)                 :: errmsg

    !
    !-----------------------------------------------------------------------
    !

    errflg = 0

    ! mass mixing ratios
    call addfld(cnst_name(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq)))
    call addfld(cnst_name(ixrain)  , (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixrain)))
    call addfld(bpcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq))//' (before physics)')
    call addfld(apcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq))//' (after physics)')

    call add_default(cnst_name(ixcldliq), 1, ' ')
    call add_default(cnst_name(ixrain),   1, ' ')

    ! Initialize Kessler with CAM physical constants
    call kessler_init(latvap, psurf_ref, rhoh2o, errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_cam_init error: Error returned from kessler_init: '//trim(errmsg))
    end if

  end subroutine kessler_cam_init

!========================================================================================

  subroutine kessler_tend(state, ptend, ztodt, pbuf)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Run Kessler physics (see kessler.F90)
    !
    !-----------------------------------------------------------------------
    use shr_kind_mod,     only: SHR_KIND_CM
    use ref_pres,         only: psurf_ref
    use physconst,        only: rh2o
    use air_composition,  only: cpairv, rairv
    use physics_types,    only: physics_state, physics_ptend
    use physics_types,    only: physics_ptend_init
    use constituents,     only: pcnst, cnst_name, cnst_type

    use kessler,          only: kessler_run
    use state_converters, only: wet_to_dry_water_vapor_run
    use state_converters, only: wet_to_dry_cloud_liquid_water_run
    use state_converters, only: wet_to_dry_rain_run
    use state_converters, only: dry_to_wet_water_vapor_run
    use state_converters, only: dry_to_wet_cloud_liquid_water_run
    use state_converters, only: dry_to_wet_rain_run
    use state_converters, only: temp_to_potential_temp_run
    use state_converters, only: calc_exner_run
    use state_converters, only: calc_dry_air_ideal_gas_density_run
    use state_converters, only: potential_temp_to_temp_run

    use cam_abortutils,   only: endrun
    use cam_history,      only: outfld


    ! arguments
    type(physics_state), intent(in)    :: state
    real(r8),            intent(in)    :: ztodt            ! physics timestep
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    !
    !---------------------------Local workspace-----------------------------
    !
    integer                            :: lchnk              ! chunk identifier
    integer                            :: ncol               ! number of atmospheric columns
    integer                            :: lyr_surf
    integer                            :: lyr_toa

    real(r8)                           :: zvirv(pcols,pver)  ! ratio of water vapor to dry air constants - 1
    real(r8)                           :: rho(pcols,pver)    ! Dry air density
    real(r8)                           :: pk(pcols,pver)     ! exner func.
    real(r8)                           :: th(pcols,pver)     ! Potential temp.
    real(r8)                           :: temp(pcols,pver)   ! temperature
    real(r8)                           :: qv(pcols,pver)     ! Water vapor mixing ratio wrt moist air
    real(r8)                           :: qc(pcols,pver)     ! Cloud water mixing ratio wrt moist air
    real(r8)                           :: qr(pcols,pver)     ! Rain mixing ratio wrt moist air
    real(r8)                           :: qv_dry(pcols,pver) ! Water vapor mixing ratio wrt dry air
    real(r8)                           :: qc_dry(pcols,pver) ! Cloud water mixing ratio wrt dry air
    real(r8)                           :: qr_dry(pcols,pver) ! Rain mixing ratio wrt dry air

    integer                            :: k,rk               ! vert. indices
    logical                            :: lq(pcnst)          ! Calc tendencies?
    character(len=SHR_KIND_CM)         :: errmsg             ! CCPP physics scheme error message

    integer                            :: errflg             ! CCPP physics scheme error flag

    character(len=64)                  :: scheme_name        ! CCPP physics scheme name (not used in CAM)

    real(r8), pointer                  :: prec_sed(:)        ! total precip from cloud sedimentation
    real(r8), pointer                  :: relhum(:,:)        ! relative humidity

    integer :: i

    !
    !-----------------------------------------------------------------------
    !
    errflg = 0
    lchnk = state%lchnk
    ncol  = state%ncol

    lyr_surf = pver
    lyr_toa = 1

    ! initialize individual parameterization tendencies
    lq           = .false.
    lq(1)        = .true.
    lq(ixcldliq) = .true.
    lq(ixrain)   = .true.
    call physics_ptend_init(ptend, state%psetcols, 'kessler',                 &
         ls=.true., lu=.true., lv=.true., lq=lq)

    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
    call pbuf_get_field(pbuf, relhum_idx,   relhum)

    do k = 1, pver
       ! Create temporaries for state variables changed by Kessler routine
       temp(:ncol,k) = state%t(:ncol,k)
       qv(:ncol,k)   = state%q(:ncol,k,1)
       qc(:ncol,k)   = state%q(:ncol,k,ixcldliq)
       qr(:ncol,k)   = state%q(:ncol,k,ixrain)

       !Also calculate gas constant ratio:
       zvirv(:ncol,k) = rh2o/rairv(:ncol,k, lchnk) - 1._r8
    end do

    ! Calculate Exner function:
    call calc_exner_run(ncol, pver, cpairv(1:ncol,:,lchnk), rairv(1:ncol,:,lchnk),   &
                        psurf_ref, state%pmid(1:ncol,:), pk(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from calc_exner_run: '//trim(errmsg))
    end if

    ! Calculate potential temperature:
    call temp_to_potential_temp_run(ncol, pver, temp(1:ncol,:), pk(1:ncol,:),        &
                                    th(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from temp_to_potential_temp_run: '//trim(errmsg))
    end if

    ! Calculate density using ideal gas law:
    call calc_dry_air_ideal_gas_density_run(ncol, pver, rairv(1:ncol,:,lchnk),       &
                                            state%pmiddry(1:ncol,:), temp(1:ncol,:), &
                                            rho(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from pres_to_density_dry_run: '//trim(errmsg))
    end if

    ! Convert moist air mixing ratios to dry air mixing ratios:
    !---------------------------------------------------------
    call wet_to_dry_water_vapor_run(ncol, pver, state%pdel(1:ncol,:),                &
                                    state%pdeldry(1:ncol,:), qv(1:ncol,:),           &
                                    qv_dry(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from wet_to_dry_water_vapor_run: '//trim(errmsg))
    end if

    call wet_to_dry_cloud_liquid_water_run(ncol, pver, state%pdel(1:ncol,:),         &
                                           state%pdeldry(1:ncol,:), qc(1:ncol,:),    &
                                           qc_dry(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from wet_to_dry_cloud_liquid_water_run: '//trim(errmsg))
    end if

    call wet_to_dry_rain_run(ncol, pver, state%pdel(1:ncol,:),                       &
                             state%pdeldry(1:ncol,:), qr(1:ncol,:),                  &
                             qr_dry(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from wet_to_dry_rain_run: '//trim(errmsg))
    end if
    !---------------------------------------------------------

    ! Run Kessler physics scheme:
    call kessler_run(ncol, pver, ztodt, lyr_surf, lyr_toa, cpairv(1:ncol,:,lchnk),   &
                     rairv(1:ncol,:,lchnk), rho(1:ncol,:), state%zm(1:ncol,:),       &
                     pk(1:ncol,:), th(1:ncol,:), qv_dry(1:ncol,:), qc_dry(1:ncol,:), &
                     qr_dry(1:ncol,:), prec_sed(1:ncol), relhum(1:ncol,:),         &
                     scheme_name, errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from kessler_run: '//trim(errmsg))
    end if

    ! Calculate air temperature from potential temperature:
    call potential_temp_to_temp_run(ncol, pver, th(1:ncol,:), pk(1:ncol,:),          &
                                    temp(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from potential_temp_to_temp_run: '//trim(errmsg))
    end if

    ! Convert dry air mixing ratios to moist air mixing ratios:
    !---------------------------------------------------------
    call dry_to_wet_water_vapor_run(ncol, pver, state%pdel(1:ncol,:),                &
                                    state%pdeldry(1:ncol,:), qv_dry(1:ncol,:),       &
                                    qv(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from dry_to_wet_water_vapor_run: '//trim(errmsg))
    end if

    call dry_to_wet_cloud_liquid_water_run(ncol, pver, state%pdel(1:ncol,:),         &
                                           state%pdeldry(1:ncol,:),                  &
                                           qc_dry(1:ncol,:), qc(1:ncol,:), errmsg,   &
                                           errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from dry_to_wet_cloud_liquid_water_run: '//trim(errmsg))
    end if

    call dry_to_wet_rain_run(ncol, pver, state%pdel(1:ncol,:),                       &
                             state%pdeldry(1:ncol,:), qr_dry(1:ncol,:),              &
                             qr(1:ncol,:), errmsg, errflg)
    if (errflg /=0) then
       call endrun('kessler_tend error: Error returned from dry_to_wet_rain_run: '//trim(errmsg))
    end if
    !---------------------------------------------------------

    ! Back out tendencies from updated fields
    do k = 1, pver
       ptend%s(:ncol,k)          = (th(:ncol,k)*pk(:ncol,k) - state%t(:ncol,k)) * cpairv(:ncol,k,lchnk) / ztodt
       ptend%q(:ncol,k,1)        = (qv(:ncol,k) - state%q(:ncol,k,1)) / ztodt
       ptend%q(:ncol,k,ixcldliq) = (qc(:ncol,k) - state%q(:ncol,k,ixcldliq)) / ztodt
       ptend%q(:ncol,k,ixrain)   = (qr(:ncol,k) - state%q(:ncol,k,ixrain)) / ztodt
    end do

    ! Output liquid tracers
    call outfld(cnst_name(ixcldliq), qc(:,pver:1:-1), pcols, lchnk)
    call outfld(cnst_name(ixrain  ), qr(:,pver:1:-1), pcols, lchnk)

  end subroutine kessler_tend
end module kessler_cam
