module mcsp_cam
!---------------------------------------------------------------------------------
! CAM interface to the mesoscale coherent structure parameterization (MCSP).
! The portable scheme lives in atmospheric_physics (schemes/mcsp/mcsp.F90).
!
! MCSP redistributes the deep convective heating and moistening vertically, and can
! add momentum tendencies, to represent organized (mesoscale) convection. Its inputs
! are the core deep convective tendencies, taken before precipitation evaporation and
! momentum transport, and the deep convective cloud top. All three come through the
! physics buffer:
!   TTEND_DP_CORE, QTEND_DP_CORE  registered here when MCSP is active and filled by
!                                 the deep convection interface (zm_conv_intr)
!   CLDTOP                        set by convect_deep_tend
! MCSP is active only when at least one mcsp_nl coefficient is positive; otherwise
! nothing is registered or called.
!---------------------------------------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc

  implicit none
  private

  public :: mcsp_cam_readnl
  public :: mcsp_cam_register
  public :: mcsp_cam_init
  public :: mcsp_cam_tend

  ! True when any MCSP coefficient is positive; set by mcsp_cam_readnl
  logical, public, protected :: do_mcsp = .false.

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Namelist variables (group mcsp_nl).
  ! A coefficient of zero turns the corresponding tendency off, so zero is the
  ! documented default; the thresholds must be set when MCSP is active.
  real(r8) :: mcsp_heat_coeff       = 0._r8     ! heating coefficient [1]
  real(r8) :: mcsp_moisture_coeff   = 0._r8     ! moistening coefficient [1]
  real(r8) :: mcsp_uwind_coeff      = 0._r8     ! zonal wind coefficient [1]
  real(r8) :: mcsp_vwind_coeff      = 0._r8     ! meridional wind coefficient [1]
  real(r8) :: mcsp_storm_speed_pref = unset_r8  ! reference pressure of the storm-level zonal wind [Pa]
  real(r8) :: mcsp_conv_depth_min   = unset_r8  ! minimum pressure depth of deep convection [Pa]
  real(r8) :: mcsp_shear_min        = unset_r8  ! minimum low-level zonal wind shear magnitude [m s-1]

  ! Physics buffer indices
  integer :: cldtop_idx        = -1
  integer :: ttend_dp_core_idx = -1
  integer :: qtend_dp_core_idx = -1

!=========================================================================================
contains
!=========================================================================================

subroutine mcsp_cam_readnl(nlfile)

   use namelist_utils, only: find_group_name
   use spmd_utils,     only: mpicom, masterprocid, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'mcsp_cam_readnl'

   namelist /mcsp_nl/ mcsp_heat_coeff, mcsp_moisture_coeff, mcsp_uwind_coeff, mcsp_vwind_coeff, &
                      mcsp_storm_speed_pref, mcsp_conv_depth_min, mcsp_shear_min
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open(newunit=unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'mcsp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, mcsp_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist mcsp_nl')
         end if
      end if
      close(unitn)
   end if

   call mpi_bcast(mcsp_heat_coeff,       1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_heat_coeff')
   call mpi_bcast(mcsp_moisture_coeff,   1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_moisture_coeff')
   call mpi_bcast(mcsp_uwind_coeff,      1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_uwind_coeff')
   call mpi_bcast(mcsp_vwind_coeff,      1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_vwind_coeff')
   call mpi_bcast(mcsp_storm_speed_pref, 1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_storm_speed_pref')
   call mpi_bcast(mcsp_conv_depth_min,   1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_conv_depth_min')
   call mpi_bcast(mcsp_shear_min,        1, mpi_real8, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname // ': FATAL: mpi_bcast: mcsp_shear_min')

   ! Active when any coefficient is positive, matching the scheme's own test
   do_mcsp = mcsp_heat_coeff  > 0._r8 .or. mcsp_moisture_coeff > 0._r8 .or. &
             mcsp_uwind_coeff > 0._r8 .or. mcsp_vwind_coeff    > 0._r8

   if (do_mcsp) then
      if (mcsp_storm_speed_pref == unset_r8 .or. mcsp_conv_depth_min == unset_r8 .or. &
          mcsp_shear_min == unset_r8) then
         call endrun(subname // ': mcsp_storm_speed_pref, mcsp_conv_depth_min and ' // &
                     'mcsp_shear_min must be set when an MCSP coefficient is positive')
      end if
   end if

   if (masterproc) then
      write(iulog,*) subname // ': MCSP active = ', do_mcsp
   end if

end subroutine mcsp_cam_readnl

!=========================================================================================

subroutine mcsp_cam_register()

   use cam_history,    only: addfld, horiz_only
   use physics_buffer, only: pbuf_add_field, dtype_r8
   !-----------------------------------------------------------------------------

   if (.not. do_mcsp) return

   call addfld('MCSP_DT',         (/ 'lev' /), 'A', 'K/s',     'MCSP temperature tendency')
   call addfld('MCSP_DQ',         (/ 'lev' /), 'A', 'kg/kg/s', 'MCSP water vapor tendency')
   call addfld('MCSP_DU',         (/ 'lev' /), 'A', 'm/s2',    'MCSP zonal wind tendency')
   call addfld('MCSP_DV',         (/ 'lev' /), 'A', 'm/s2',    'MCSP meridional wind tendency')
   call addfld('MCSP_DT_max',     horiz_only,  'A', 'K/s',     'MCSP heating amplitude')
   call addfld('MCSP_freq',       horiz_only,  'A', '1',       'MCSP frequency of activation')
   call addfld('MCSP_shear',      horiz_only,  'A', 'm/s',     'MCSP low-level zonal wind shear')
   call addfld('MCSP_conv_depth', horiz_only,  'A', 'Pa',      'Deep convection pressure depth for MCSP')

   ! Core deep convective tendencies (before precipitation evaporation and momentum
   ! transport), filled by the deep convection interface when these fields exist
   call pbuf_add_field('TTEND_DP_CORE', 'physpkg', dtype_r8, (/pcols,pver/), ttend_dp_core_idx)
   call pbuf_add_field('QTEND_DP_CORE', 'physpkg', dtype_r8, (/pcols,pver/), qtend_dp_core_idx)

end subroutine mcsp_cam_register

!=========================================================================================

subroutine mcsp_cam_init()

   use mcsp,           only: mcsp_init
   use phys_control,   only: phys_getopts
   use physics_buffer, only: pbuf_get_index

   ! Local variables
   character(len=16)  :: deep_scheme
   character(len=512) :: errmsg
   integer            :: errflg
   !-----------------------------------------------------------------------------

   if (.not. do_mcsp) return

   ! Only ZM fills TTEND_DP_CORE and QTEND_DP_CORE
   call phys_getopts(deep_scheme_out=deep_scheme)
   if (deep_scheme /= 'ZM') then
      call endrun('mcsp_cam_init: MCSP requires deep_scheme = ZM, got ' // trim(deep_scheme))
   end if

   call mcsp_init(masterproc, iulog, &
                  mcsp_heat_coeff, mcsp_moisture_coeff, mcsp_uwind_coeff, mcsp_vwind_coeff, &
                  mcsp_storm_speed_pref, mcsp_conv_depth_min, mcsp_shear_min, &
                  errmsg, errflg)
   if (errflg /= 0) then
      call endrun('mcsp_cam_init: ' // trim(errmsg))
   end if

   ! Deep convective cloud top, set by convect_deep_tend before MCSP runs
   cldtop_idx = pbuf_get_index('CLDTOP')

end subroutine mcsp_cam_init

!=========================================================================================

subroutine mcsp_cam_tend(state, ptend, ztodt, pbuf)

   use mcsp,           only: mcsp_run
   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use physconst,      only: cpair, pi
   use constituents,   only: pcnst
   use cam_history,    only: outfld

   ! Arguments
   type(physics_state), intent(in)  :: state   ! physics state variables
   type(physics_ptend), intent(out) :: ptend   ! MCSP tendencies
   real(r8),            intent(in)  :: ztodt   ! physics time step [s]
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! Local variables
   integer  :: lchnk
   integer  :: ncol
   logical  :: lq(pcnst)

   real(r8), pointer :: jctop(:)            ! deep convective cloud top level as a real [index]
   real(r8), pointer :: ttend_dp_core(:,:)  ! core deep convective temperature tendency [K/s]
   real(r8), pointer :: qtend_dp_core(:,:)  ! core deep convective water vapor tendency [kg/kg/s]

   integer  :: jctop_int(pcols)             ! deep convective cloud top level [index]
   real(r8) :: mcsp_dt_out(pcols,pver)      ! diagnostic temperature tendency [K/s]
   real(r8) :: mcsp_dq_out(pcols,pver)      ! diagnostic water vapor tendency [kg/kg/s]
   real(r8) :: mcsp_du_out(pcols,pver)      ! diagnostic zonal wind tendency [m/s2]
   real(r8) :: mcsp_dv_out(pcols,pver)      ! diagnostic meridional wind tendency [m/s2]
   real(r8) :: mcsp_freq(pcols)             ! 1 where MCSP contributed a tendency [1]
   real(r8) :: mcsp_shear(pcols)            ! low-level zonal wind shear [m/s]
   real(r8) :: conv_depth(pcols)            ! pressure depth of deep convection [Pa]
   real(r8) :: mcsp_dt_max(pcols)           ! heating amplitude [K/s]

   character(len=512) :: errmsg
   integer            :: errflg
   !-----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   lq(:) = .false.
   lq(1) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'mcsp', ls=.true., lq=lq, lu=.true., lv=.true.)

   call pbuf_get_field(pbuf, cldtop_idx,        jctop)
   call pbuf_get_field(pbuf, ttend_dp_core_idx, ttend_dp_core)
   call pbuf_get_field(pbuf, qtend_dp_core_idx, qtend_dp_core)

   jctop_int(:)     = pver
   jctop_int(:ncol) = int(jctop(:ncol))

   mcsp_dt_out(:,:) = 0._r8
   mcsp_dq_out(:,:) = 0._r8
   mcsp_du_out(:,:) = 0._r8
   mcsp_dv_out(:,:) = 0._r8
   mcsp_freq(:)     = 0._r8
   mcsp_shear(:)    = 0._r8
   conv_depth(:)    = 0._r8
   mcsp_dt_max(:)   = 0._r8

   call mcsp_run(ncol, pver, pverp, cpair, pi, ztodt, jctop_int(:ncol),            &
                 state%pmid(:ncol,:), state%pint(:ncol,:), state%pdel(:ncol,:),     &
                 state%u(:ncol,:), state%v(:ncol,:),                                &
                 ttend_dp_core(:ncol,:), qtend_dp_core(:ncol,:),                    &
                 ptend%s(:ncol,:), ptend%q(:ncol,:,1), ptend%u(:ncol,:), ptend%v(:ncol,:), &
                 mcsp_dt_out(:ncol,:), mcsp_dq_out(:ncol,:),                        &
                 mcsp_du_out(:ncol,:), mcsp_dv_out(:ncol,:),                        &
                 mcsp_freq(:ncol), mcsp_shear(:ncol), conv_depth(:ncol), mcsp_dt_max(:ncol), &
                 errmsg, errflg)
   if (errflg /= 0) then
      call endrun('mcsp_cam_tend: ' // trim(errmsg))
   end if

   call outfld('MCSP_DT',         mcsp_dt_out, pcols, lchnk)
   call outfld('MCSP_DQ',         mcsp_dq_out, pcols, lchnk)
   call outfld('MCSP_DU',         mcsp_du_out, pcols, lchnk)
   call outfld('MCSP_DV',         mcsp_dv_out, pcols, lchnk)
   call outfld('MCSP_DT_max',     mcsp_dt_max, pcols, lchnk)
   call outfld('MCSP_freq',       mcsp_freq,   pcols, lchnk)
   call outfld('MCSP_shear',      mcsp_shear,  pcols, lchnk)
   call outfld('MCSP_conv_depth', conv_depth,  pcols, lchnk)

end subroutine mcsp_cam_tend

end module mcsp_cam
