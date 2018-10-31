module kessler_cam

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver
  use physics_buffer, only: physics_buffer_desc, pbuf_init_time, dtype_r8, &
                            pbuf_add_field, pbuf_get_field

  implicit none
  private
  save

  public :: kessler_register, kessler_init, kessler_tend

  integer :: ixcldliq = -1 ! cloud liquid mixing ratio index
  integer :: ixrain   = -1 ! rain liquid mixing ratio index

  ! physics buffer indices
  integer :: prec_sed_idx  = 0

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

  end subroutine kessler_register

!========================================================================================

  subroutine kessler_init(pbuf2d)

    use physconst,      only: cpair, latvap,pstd, rair, rhoh2o
    use constituents,   only: cnst_name, cnst_longname, bpcnst, apcnst
    use cam_history,    only: addfld, add_default, horiz_only
    use kessler_mod,    only: kessler_set_const

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! mass mixing ratios
    call addfld(cnst_name(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq)))
    call addfld(cnst_name(ixrain)  , (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixrain)))
    call addfld(bpcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq))//' (before physics)')
    call addfld(apcnst(ixcldliq),    (/ 'lev' /), 'A', 'kg/kg', trim(cnst_longname(ixcldliq))//' (after physics)')

    call add_default(cnst_name(ixcldliq), 1, ' ')
    call add_default(cnst_name(ixrain),   1, ' ')

    ! Initialize Kessler with CAM physical constants
    call kessler_set_const(rair, cpair, latvap, pstd/100.0_r8, rhoh2o)

  end subroutine kessler_init

!========================================================================================

  subroutine kessler_tend(state, ptend, ztodt, pbuf)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Run Kessler physics (see kessler.F90)
    ! 
    !-----------------------------------------------------------------------
    use shr_kind_mod,       only: SHR_KIND_CM
    use physconst,          only: cpair, rair, zvir
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use constituents,       only: pcnst, cnst_name, cnst_type

    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use kessler_mod,        only: kessler


    ! arguments
    type(physics_state), intent(in)    :: state
    real(r8),            intent(in)    :: ztodt            ! physics timestep
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    !---------------------------Local workspace-----------------------------
    !
    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    real(r8)                           :: rho(pcols,pver)  ! Dry air density
    real(r8)                           :: exner(pcols,pver)! exner (CAM vertical order)
    real(r8)                           :: pk(pcols,pver)   ! exner func.
    real(r8)                           :: th(pcols,pver)   ! Potential temp.
    real(r8)                           :: qv(pcols,pver)   ! Water vapor
    real(r8)                           :: qc(pcols,pver)   ! Cloud water
    real(r8)                           :: qr(pcols,pver)   ! Rain water
    real(r8)                           :: z(pcols,pver)    ! height
    real(r8)                           :: wet_to_dry(pcols)! factor to convert from wet to dry mixing ratio 
    real(r8)                           :: dry_to_wet(pcols)! factor to convert from dry to wet mixing ratio 
    integer                            :: k,rk             ! vert. indices
    logical                            :: lq(pcnst)        ! Calc tendencies?
    character(len=SHR_KIND_CM)         :: errmsg

    real(r8), pointer                  :: prec_sed(:) ! total precip from cloud sedimentation

    integer :: i

    !
    !-----------------------------------------------------------------------
    !
    lchnk = state%lchnk
    ncol  = state%ncol

    ! initialize individual parameterization tendencies
    lq           = .false.
    lq(1)        = .true.
    lq(ixcldliq) = .true.
    lq(ixrain)   = .true.
    call physics_ptend_init(ptend, state%psetcols, 'kessler',                 &
         ls=.true., lu=.true., lv=.true., lq=lq)

    do k=1,pver
      exner(:ncol,k) = (state%pmid(:ncol,k)/1.e5_r8)**(rair/cpair)
    end do

    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)

    do k = 1, pver
      rk = pver - k + 1
      rho(:ncol,rk) = state%pmiddry(:ncol,k)/(rair*state%t(:ncol,k))
      pk(:ncol,rk) =  exner(:ncol,k) 
      ! Create temporaries for state variables changed by Kessler routine
      th(:ncol,rk) = state%t(:ncol,k) / exner(:ncol,k)
      z(:ncol,rk)  = state%zm(:ncol,k)
      qv(:ncol,rk) = state%q(:ncol,k,1)
      qc(:ncol,rk) = state%q(:ncol,k,ixcldliq)
      qr(:ncol,rk) = state%q(:ncol,k,ixrain)
      !
      ! mixing ratios are wet - convert to dry for Kessler physics
      !
      wet_to_dry(:ncol) = state%pdel(:ncol,k)/state%pdeldry(:ncol,k)

      if (cnst_type(1).eq.'wet')        qv(:ncol,rk) = wet_to_dry(:ncol)*qv(:ncol,rk)     
      if (cnst_type(ixcldliq).eq.'wet') qc(:ncol,rk) = wet_to_dry(:ncol)*qc(:ncol,rk)
      if (cnst_type(ixrain).eq.'wet')   qr(:ncol,rk) = wet_to_dry(:ncol)*qr(:ncol,rk)
    end do

    ! Kessler physics arguments
    ! ncol:  Number of columns
    ! nz:    Number of vertical levels
    ! dt:    Time step (s) (in)
    ! rho:   Dry air density (not mean state as in KW) (kg/m^3) (in)
    ! z:     Heights of thermo. levels in the grid column (m) (in)
    ! pk:    Exner function (p/p0)**(R/cp) (in)
    ! th:    Potential Temperature (K) (inout)
    ! qv:    Water vapor mixing ratio (gm/gm) (inout)
    ! qc:    Cloud water mixing ratio (gm/gm) (inout)
    ! qr:    Rain  water mixing ratio (gm/gm) (inout)
    ! precl: Precipitation rate (m_water / s) (out)
    ! errmsg: Error string if error found
    call kessler(ncol, pver, ztodt, rho(:ncol,:), z(:ncol,:), pk(:ncol,:),    &
         th(:ncol,:), qv(:ncol,:), qc(:ncol,:), qr(:ncol,:), prec_sed, errmsg)

    if (len_trim(errmsg) > 0) then
      call endrun(trim(errmsg))
    end if

    do k = 1, pver
      rk = pver - k + 1
      !
      ! mixing ratios are dry - convert to wet
      !
      dry_to_wet(:ncol) = state%pdeldry(:ncol,k)/state%pdel(:ncol,k)

      if (cnst_type(1).eq.'wet')        qv(:ncol,rk) = dry_to_wet(:ncol)*qv(:ncol,rk)     
      if (cnst_type(ixcldliq).eq.'wet') qc(:ncol,rk) = dry_to_wet(:ncol)*qc(:ncol,rk)
      if (cnst_type(ixrain).eq.'wet')   qr(:ncol,rk) = dry_to_wet(:ncol)*qr(:ncol,rk)
    end do


    ! Back out tendencies from updated fields
    do k = 1, pver
      rk = pver - k + 1
      ptend%s(:ncol,k) = (th(:ncol,rk)*exner(:ncol,k) - state%t(:ncol,k)) * cpair / ztodt
      ptend%q(:ncol,k,1) = (qv(:ncol,rk) - state%q(:ncol,k,1)) / ztodt
      ptend%q(:ncol,k,ixcldliq) = (qc(:ncol,rk) - state%q(:ncol,k,ixcldliq)) / ztodt
      ptend%q(:ncol,k,ixrain) = (qr(:ncol,rk) - state%q(:ncol,k,ixrain)) / ztodt
    end do

    ! Output liquid tracers
    call outfld(cnst_name(ixcldliq), qc(:,pver:1:-1), pcols, lchnk)
    call outfld(cnst_name(ixrain  ), qr(:,pver:1:-1), pcols, lchnk)

  end subroutine kessler_tend
end module kessler_cam
