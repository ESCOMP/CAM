module apply_iop_forcing_mod

use shr_kind_mod,   only:r8 => shr_kind_r8, i8 => shr_kind_i8
use pmgrid,         only:plev, plevp, plon
use constituents,   only:pcnst, cnst_get_ind, cnst_name
use physconst,      only:rair,cpair
use cam_logfile,    only:iulog
use hybvcoord_mod,  only: hvcoord_t
use scamMod,        only: use_3dfrc, single_column, have_u, have_v, divT3d, divq3d, divt, divq, &
                          wfld, uobs, vobs, tobs, qobs, plevs0, have_divt3d, have_divq3d, &
                          scm_relax_bot_p,scm_relax_linear,scm_relax_tau_bot_sec, &
                          scm_relax_tau_sec,scm_relax_tau_top_sec,scm_relax_top_p, &
                          scm_relaxation,scm_relax_fincl,qinitobs

use cam_abortutils, only: endrun
use string_utils,   only: to_upper

implicit none

public advance_iop_forcing
public advance_iop_nudging

!=========================================================================
contains
!=========================================================================

subroutine advance_iop_forcing(scm_dt, ps_in, &                   ! In
                    u_in, v_in, t_in, q_in, t_phys_frc, q_phys_frc, hvcoord, & ! In
                    u_update, v_update, t_update, q_update)       ! Out

!-----------------------------------------------------------------------
!
! Purpose:
! Apply large scale forcing for t, q, u, and v as provided by the
!   case IOP forcing file.
!
! Author:
! Original version: Adopted from CAM3.5/CAM5
! Updated version for E3SM: Peter Bogenschutz (bogenschutz1@llnl.gov)
!  and replaces the forecast.F90 routine in CAM3.5/CAM5/CAM6/E3SMv1/E3SMv2
!
!-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: ps_in             ! surface pressure [Pa]
  real(r8), intent(in) :: u_in(plev)        ! zonal wind [m/s]
  real(r8), intent(in) :: v_in(plev)        ! meridional wind [m/s]
  real(r8), intent(in) :: t_in(plev)        ! temperature [K]
  real(r8), intent(in) :: q_in(plev,pcnst)  ! q tracer array [units vary] already vertically advected
  real(r8), intent(in) :: t_phys_frc(plev)  ! temperature forcing from physics [K/s]
  real(r8), intent(in) :: q_phys_frc(plev,pcnst)  ! change in q due to physics.
  type (hvcoord_t), intent(in)   :: hvcoord
  real(r8), intent(in) :: scm_dt            ! model time step [s]

  ! Output arguments
  real(r8), intent(out) :: t_update(plev)      ! updated temperature [K]
  real(r8), intent(out) :: q_update(plev,pcnst)! updated q tracer array [units vary]
  real(r8), intent(out) :: u_update(plev)      ! updated zonal wind [m/s]
  real(r8), intent(out) :: v_update(plev)      ! updated meridional wind [m/s]

  ! Local variables
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
  real(r8) t_lsf(plev)       ! storage for temperature large scale forcing
  real(r8) q_lsf(plev,pcnst) ! storage for moisture large scale forcing
  real(r8) fac, t_expan

  integer i,k,m           ! longitude, level, constituent indices

  character(len=*), parameter :: subname = 'advance_iop_forcing'

  ! Get vertical level profiles
  call plevs0(plev, ps_in, hvcoord%ps0, hvcoord%hyam, hvcoord%hybm, hvcoord%hyai, hvcoord%hybi, pintm1 ,pmidm1 ,pdelm1)

  !  Advance T and Q due to large scale forcing
  if (use_3dfrc) then
    if(.not.(have_divt3d.and.have_divq3d)) call endrun(subname//": FATAL: divt3d and divq3d not available")
    t_lsf(:plev) = divt3d(:plev)
    q_lsf(:plev,:pcnst) = divq3d(:plev,:pcnst)
  else
    t_lsf(:plev) = divt(:plev)
    q_lsf(:plev,:pcnst) = divq(:plev,:pcnst)
  endif

  do k=1,plev
     ! Initialize thermal expansion term to zero.  This term is only
     !  considered if three dimensional forcing is not provided by IOP forcing file.
     t_expan = 0._r8

    if (.not. use_3dfrc) then
      t_expan = scm_dt*wfld(k)*t_in(k)*rair/(cpair*pmidm1(k))
    endif

     if (use_3dfrc) then
        do m=1,pcnst
           ! When using 3d dynamics tendencies, SCM skips the vertical advection step and thus
           ! q_in at this point has not had physics tendencies applied
           q_update(k,m) = q_in(k,m) + scm_dt*(q_phys_frc(k,m) + q_lsf(k,m))
        end do
        t_update(k)   = t_in(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k))
     else
        do m=1,pcnst
           ! When not using 3d dynamics tendencies, q_in at this point has had physics tend
           ! applied and has been vertically advected. Only horizontal dyn tend needed for forecast.
           q_update(k,m) = q_in(k,m) + scm_dt*q_lsf(k,m)
        end do
        t_update(k)   = t_in(k) + t_expan + scm_dt*t_lsf(k)
     end if
  end do

  !  Set U and V fields

  if ( have_v .and. have_u ) then
    do k=1,plev
      u_update(k) = uobs(k)
      v_update(k) = vobs(k)
    enddo
  endif

end subroutine advance_iop_forcing

!=========================================================================

subroutine advance_iop_nudging(ztodt, ps_in,                        &      ! In
                               tfcst, qfcst, ufcst, vfcst, hvcoord, &      ! Inout
                               relaxt, relaxq )                            ! Out

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Option to nudge t and q to observations as specified by the IOP file
  !-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: ztodt            ! model time step [s]
  real(r8), intent(in) :: ps_in            ! surface pressure [Pa]
  type (hvcoord_t), intent(in)   :: hvcoord

  ! Output arguments
  real(r8), intent(inout) :: tfcst(plev)  ! updated temperature [K]
  real(r8), intent(inout) :: qfcst(plon,plev,pcnst)  ! updated const field
  real(r8), intent(inout) :: ufcst(plev)  ! updated U wind
  real(r8), intent(inout) :: vfcst(plev)  ! updated V wind
  real(r8), intent(out) :: relaxt(plev)    ! relaxation of temperature [K/s]
  real(r8), intent(out) :: relaxq(plev)    ! relaxation of vapor [kg/kg/s]

  ! Local variables
  integer :: i, k, m
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)

  ! --------------------------- !
  ! For 'scm_relaxation' switch !
  ! --------------------------- !

  real(r8) rtau(plev)
  real(r8) relax_T(plev)
  real(r8) relax_u(plev)
  real(r8) relax_v(plev)
  real(r8) relax_q(plev,pcnst)
  ! +++BPM: allow linear relaxation profile
  real(r8) rslope ! [optional] slope for linear relaxation profile
  real(r8) rycept ! [optional] y-intercept for linear relaxtion profile
  logical scm_fincl_empty

  ! ------------------------------------------------------------------- !
  ! Relaxation to the observed or specified state                       !
  ! We should specify relaxation time scale ( rtau ) and                !
  ! target-relaxation state ( in the current case, either 'obs' or 0 )  !
  ! ------------------------------------------------------------------- !

  if ( .not. scm_relaxation) return

  call plevs0(plev, ps_in, hvcoord%ps0, hvcoord%hyam, hvcoord%hybm, hvcoord%hyai, hvcoord%hybi, pintm1 ,pmidm1 ,pdelm1)

  relax_T(:)             = 0._r8
  relax_u(:)             = 0._r8
  relax_v(:)             = 0._r8
  relax_q(:plev,:pcnst)  = 0._r8
  ! +++BPM: allow linear relaxation profile
  ! scm_relaxation is a logical from scamMod
  ! scm_relax_tau_top_sec and scm_relax_tau_bot_sec are the relaxation times at top and bottom of layer
  ! also defined in scamMod
  if ( scm_relax_linear ) then
     rslope = (scm_relax_top_p - scm_relax_bot_p)/(scm_relax_tau_top_sec - scm_relax_tau_bot_sec)
     rycept = scm_relax_tau_top_sec - (rslope*scm_relax_top_p)
  endif

  scm_fincl_empty=.true.
  do i=1,pcnst
     if (len_trim(scm_relax_fincl(i)) > 0) then
        scm_fincl_empty=.false.
        scm_relax_fincl(i)=trim(to_upper(scm_relax_fincl(i)))
     end if
  end do

  do k = 1, plev
     if ( pmidm1(k) <= scm_relax_bot_p.and.pmidm1(k) >= scm_relax_top_p ) then ! inside layer
        if (scm_relax_linear) then
           rtau(k) = rslope*pmidm1(k) + rycept ! linear regime
        else
           rtau(k)         = max( ztodt, scm_relax_tau_sec ) ! constant for whole layer / no relax outside
        endif
     else if  (scm_relax_linear .and. pmidm1(k) <= scm_relax_top_p ) then ! not linear => do nothing / linear => use upper value
        rtau(k) = scm_relax_tau_top_sec ! above layer keep rtau equal to the top
     endif
     ! +BPM: this can't be the best way...
     ! I put this in because if rtau doesn't get set above, then I don't want to do any relaxation in that layer.
     ! maybe the logic of this whole loop needs to be re-thinked.
     if (rtau(k) /= 0) then
        relax_T(k)      = -  ( tfcst(k)     - tobs(k) )    / rtau(k)
        relax_u(k)      = -  ( ufcst(k)     - uobs(k) )    / rtau(k)
        relax_v(k)      = -  ( vfcst(k)     - vobs(k) )    / rtau(k)
        relax_q(k,1)    = -  ( qfcst(1,k,1) - qobs(k) )    / rtau(k)
        do m = 2, pcnst
           relax_q(k,m) = -  ( qfcst(1,k,m) - qinitobs(k,m)   )    / rtau(k)
        enddo
        if (scm_fincl_empty .or. ANY(scm_relax_fincl(:) == 'T')) &
             tfcst(k)        =      tfcst(k)     + relax_T(k)   * ztodt
        if (scm_fincl_empty .or.ANY(scm_relax_fincl(:) == 'U')) &
             ufcst(k)        =      ufcst(k)     + relax_u(k)   * ztodt
        if (scm_fincl_empty .or. ANY(scm_relax_fincl(:) == 'V')) &
             vfcst(k)        =      vfcst(k)     + relax_v(k)   * ztodt
        do m = 1, pcnst
           if (scm_fincl_empty .or. ANY(scm_relax_fincl(:) == trim(to_upper(cnst_name(m)))) ) then
              qfcst(1,k,m) =      qfcst(1,k,m) + relax_q(k,m) * ztodt
           end if
        enddo
     end if
  enddo

end subroutine advance_iop_nudging

!-----------------------------------------------------------------------

end module apply_iop_forcing_mod
