module apply_iop_forcing_mod

use shr_kind_mod,   only:r8 => shr_kind_r8, i8 => shr_kind_i8
use pmgrid
use constituents,   only:pcnst, cnst_get_ind
use physconst,      only:rair,cpair
use cam_logfile,    only:iulog
use hybvcoord_mod,  only: hvcoord_t
use scamMod,        only: use_3dfrc, single_column, have_u, have_v, divT3d, divq3d, divt, divq, &
                          wfld, uobs, vobs, tobs, qobs, plevs0, have_divt3d, have_divq3d
use cam_abortutils, only: endrun
implicit none

public advance_iop_forcing
public advance_iop_nudging

!=========================================================================
contains
!=========================================================================

!$$subroutine advance_iop_forcing(scm_dt     ,hvcoord      ,psm1       ,  & ! In
!$$                               u3m1       ,u3m2         ,u_phys_frc ,  & ! In
!$$                               v3m1       ,v3m2         ,v_phys_frc ,  & ! In
!$$                               t3m1       ,t3m2         ,t_phys_frc ,  & ! In
!$$                               q3m1       ,q3m2         ,q_phys_frc ,  & ! In
!$$                               u3         ,v3           ,t3         ,  & ! Out
!$$                               q3 )                                      ! Out
!$$!-----------------------------------------------------------------------
!$$!
!$$! Purpose:
!$$! --------------------------------------------------------------------------- ! 
!$$!                                                                             !
!$$! Compute Time-Marched 'T, u, v, q' for SCAM by summing the 'physics',        !
!$$! 'horizontal advection', and 'vertical advection' tendencies.                ! 
!$$! This module is used only for SCAM.                                          ! 
!$$!                                                                             ! 
!$$! --------------------------------------------------------------------------- ! 
!$$!
!$$! Author:
!$$! Original version: Adopted from CAM3.5/CAM5
!$$! Substantially refactored by Sung Su Park
!$$! Updated for SE dycore: Peter Bogenschutz (bogenschutz1@llnl.gov) John Truesdale (jet.ucar.edu)
!$$!
!$$! --------------------------------------------------- 
!$$! x = t, u, v, q                                      
!$$! x3m1 : state variable used for computing 'forcing'  
!$$! x3m2 : initial state variable before time-marching  
!$$! x3   : final   state variable after  time-marching  
!$$!-----------------------------------------------------------------------
!$$
!$$! Input arguments
!$$  real(r8), intent(in) :: scm_dt            ! model time step [s]
!$$  type (hvcoord_t), intent(in)   :: hvcoord
!$$  real(r8), intent(in) :: psm1              ! surface pressure [Pa]
!$$  real(r8), intent(in) :: u3m1(plev)        ! Zonal wind [ m/s ]
!$$  real(r8), intent(in) :: u3m2(plev)        ! Zonal wind [ m/s ]
!$$  real(r8), intent(in) :: u_phys_frc(plev)        ! Zonal wind [ m/s ]
!$$  real(r8), intent(in) :: v3m1(plev)        ! Meridional wind [ m/s ]
!$$  real(r8), intent(in) :: v3m2(plev)        ! Meridional wind [ m/s ]
!$$  real(r8), intent(in) :: v_phys_frc(plev)        ! Meridional wind [ m/s ]
!$$  real(r8), intent(in) :: t3m1(plev)        ! temperature [K]
!$$  real(r8), intent(in) :: t3m2(plev)        ! temperature [K]
!$$  real(r8), intent(in) :: t_phys_frc(plev)  ! temperature forcing from physics [K/s]
!$$  real(r8), intent(in) :: q3m1(plev,pcnst)        ! Meridional wind [ m/s ]
!$$  real(r8), intent(in) :: q3m2(plev,pcnst)        ! Meridional wind [ m/s ]
!$$  real(r8), intent(in) :: q_phys_frc(plev,pcnst)  ! temperature forcing from physics [K/s]
!$$
!$$  ! Output arguments
!$$  real(r8), intent(out) :: t3(plev)      ! updated temperature [K]
!$$  real(r8), intent(out) :: q3(plev,pcnst)! updated q tracer array [units vary]
!$$  real(r8), intent(out) :: u3(plev)      ! updated zonal wind [m/s]
!$$  real(r8), intent(out) :: v3(plev)      ! updated meridional wind [m/s]
!$$
!$$  ! Local variables
!$$  real(r8) pmidm1(plev)  ! pressure at model levels
!$$  real(r8) pintm1(plevp) ! pressure at model interfaces
!$$  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)
!$$  real(r8) t_lsf(plev)       ! storage for temperature large scale forcing
!$$  real(r8) q_lsf(plev,pcnst) ! storage for moisture large scale forcing
!$$  real(r8) qten_phys(plev,pcnst) ! storage for moisture large scale forcing
!$$  real(r8) fac, t_expan
!$$
!$$  integer i,k,m           ! longitude, level, constituent indices
!$$  integer  ixcldliq, ixcldice, ixnumliq, ixnumice
!$$
!$$  ! Main Computation Begins Here                                 !
!$$  call cnst_get_ind( 'CLDLIQ', ixcldliq, abort=.false. )
!$$  call cnst_get_ind( 'CLDICE', ixcldice, abort=.false. )
!$$  call cnst_get_ind( 'NUMLIQ', ixnumliq, abort=.false. )
!$$  call cnst_get_ind( 'NUMICE', ixnumice, abort=.false. )
!$$
!$$  ! Calculate midpoint pressure levels                           !
!$$  call plevs0(plev    ,psm1   ,pintm1 ,pmidm1 ,pdelm1, hvcoord)
!$$
!$$
!$$   ! Extract physical tendencies of tracers q.
!$$   ! Note 'tten_phys, uten_phys, vten_phys' are already input.
!$$     qten_phys(:plev,:pcnst)     = ( q_phys_frc(:plev,:pcnst) - q3m2(:plev,:pcnst) ) / scm_dt
!$$
!$$
!$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$  !  Advance T and Q due to large scale forcing
!$$
!$$  if (use_3dfrc .or. have_divt3d) then
!$$    t_lsf(:plev) = divt3d(:plev)
!$$    q_lsf(:plev,:pcnst) = divq3d(:plev,:pcnst)
!$$  else
!$$    t_lsf(:plev) = divt(:plev)
!$$    q_lsf(:plev,:pcnst) = divq(:plev,:pcnst)
!$$  endif
!$$
!$$  do k=1,plev
!$$    ! Initialize thermal expansion term to zero.  This term is only
!$$    !  considered if using the preq-x dycore and if three dimensional
!$$    !  forcing is not provided by IOP forcing file.
!$$    t_expan = 0._r8
!$$
!$$!!$    if (.not. use_3dfrc) then
!$$!!$      t_expan = scm_dt*wfld(k)*t3m1(k)*rair/(cpair*pmidm1(k))
!$$!!$    endif
!$$!!$    t_update(k) = t3m1(k) + t_expan + scm_dt*(t_phys_frc(k) + t_lsf(k))
!$$    t3(k) = t3m1(k) + scm_dt*(t_phys_frc(k) + t_lsf(k))
!$$    do m=1,pcnst
!$$      q3(k,m) = scm_dt*(q_phys_frc(k,m) + q_lsf(k,m))
!$$    end do
!$$  enddo
!$$
!$$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$  !  Set U and V fields
!$$
!$$  if ( have_v .and. have_u ) then
!$$    do k=1,plev
!$$      u3(k) = uobs(k)
!$$      v3(k) = vobs(k)
!$$    enddo
!$$  endif
!$$
!$$end subroutine advance_iop_forcing
!$$
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
  real(r8), intent(in) :: q_phys_frc(plev,pcnst)  ! temperature forcing from physics [K/s]
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
  integer nlon

  character(len=*), parameter :: subname = 'advance_iop_forcing'

  !! Get vertical level profiles

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(plev    ,ps_in   ,pintm1 ,pmidm1 ,pdelm1, hvcoord)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set U and V fields

  if ( have_v .and. have_u ) then
    do k=1,plev
      u_update(k) = uobs(k)
      v_update(k) = vobs(k)
    enddo
  endif

end subroutine advance_iop_forcing

!=========================================================================

subroutine advance_iop_nudging(scm_dt, ps_in, t_in, q_in, hvcoord, &   ! In
                             t_update, q_update, relaxt, relaxq )      ! Out

!-----------------------------------------------------------------------
!
! Purpose:
! Option to nudge t and q to observations as specified by the IOP file
!-----------------------------------------------------------------------

  ! Input arguments
  real(r8), intent(in) :: scm_dt           ! model time step [s]
  real(r8), intent(in) :: ps_in            ! surface pressure [Pa]
  real(r8), intent(in) :: t_in(plev)       ! temperature [K]
  real(r8), intent(in) :: q_in(plev)       ! water vapor mixing ratio [kg/kg]
  type (hvcoord_t), intent(in)   :: hvcoord

  ! Output arguments
  real(r8), intent(out) :: t_update(plev)  ! updated temperature [K]
  real(r8), intent(out) :: q_update(plev)  ! updated water vapor [kg/kg]
  real(r8), intent(out) :: relaxt(plev)    ! relaxation of temperature [K/s]
  real(r8), intent(out) :: relaxq(plev)    ! relaxation of vapor [kg/kg/s]

  ! Local variables
  integer :: k, nlon
  real(r8) rtau(plev)
  real(r8) pmidm1(plev)  ! pressure at model levels
  real(r8) pintm1(plevp) ! pressure at model interfaces
  real(r8) pdelm1(plev)  ! pdel(k)   = pint  (k+1)-pint  (k)

  nlon = 1 ! number of columns for plevs0 routine
  call plevs0(plev    ,ps_in   ,pintm1 ,pmidm1 ,pdelm1, hvcoord)

  ! Set relaxation arrays to zero
  do k=1,plev
    relaxt(k) = 0.0_r8
    relaxq(k) = 0.0_r8
  end do

  do k=1,plev

!jt    if (pmidm1(k) .le. iop_nudge_tq_low*100._r8 .and. &
!jt      pmidm1(k) .ge. iop_nudge_tq_high*100._r8) then

      ! Set the relaxation time scale
!jt      rtau(k)   = iop_nudge_tscale
!jt      rtau(k)   = max(scm_dt,rtau(k))
      rtau(k)   = scm_dt
      relaxt(k) = -(t_update(k) - tobs(k))/rtau(k)
      relaxq(k) = -(q_update(k) - qobs(k))/rtau(k)

      t_update(k) = t_update(k) + relaxt(k)*scm_dt
      q_update(k) = q_update(k) + relaxq(k)*scm_dt

!jt    endif

  end do

end subroutine advance_iop_nudging

!-----------------------------------------------------------------------

end module apply_iop_forcing_mod


