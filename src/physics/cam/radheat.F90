
module radheat
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
!
! Original version: B.A. Boville
! Change weighting function for RRTMG: A J Conley
!-----------------------------------------------------------------------

! Use a cubic polynomial over the domain from minimum pressure to maximum pressure
! Cubic polynomial is chosen so that derivative is zero at minimum and maximum pressures
!   and is monotonically increasing from zero at minimum pressure to one at maximum pressure

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver
  use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
  use physconst,       only: cpair,mwco2
  use air_composition, only: cpairv
  use perf_mod
  use cam_logfile,     only: iulog

  implicit none
  private
  save

! Public interfaces
  public  &
       radheat_readnl,        &!
       radheat_register,      &!
       radheat_init,          &!
       radheat_timestep_init, &!
       radheat_tend            ! return net radiative heating

  public :: radheat_disable_waccm ! disable waccm heating in the upper atm

  real(r8), public :: p_top_for_equil_rad = 0._r8

! Private variables for merging heating rates
  real(r8):: qrs_wt(pver)             ! merge weight for cam solar heating
  real(r8):: qrl_wt(pver)             ! merge weight for cam long wave heating

  ! lw merge region
  ! highest altitude (lowest  pressure) of merge region (Pa)
  real(r8) :: min_pressure_lw= 5._r8
  ! lowest  altitude (highest pressure) of merge region (Pa)
  real(r8) :: max_pressure_lw=50._r8

  integer :: ntop_qrs_cam             ! top level for pure cam solar heating

!===============================================================================
contains
!===============================================================================
subroutine radheat_readnl(nlfile)
   ! Read radheat_nl namelist group
   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_real8
   use cam_abortutils,  only: endrun

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=32) :: errmsg
   character(len=*), parameter :: sub = 'radheat_readnl'

   namelist /radheat_nl/ p_top_for_equil_rad
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'radheat_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radheat_nl, iostat=ierr, iomsg=errmsg)
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist: '//trim(errmsg))
         end if
      end if
      close(unitn)
   end if

   call mpi_bcast(p_top_for_equil_rad, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: p_top_for_equil_rad")

end subroutine radheat_readnl

!================================================================================================

  subroutine radheat_register

  ! No options for this version of radheat; this is just a stub.

  end subroutine radheat_register

!================================================================================================
!================================================================================================
!================================================================================================

  subroutine radheat_init(pref_mid)

    
    use nlte_fomichev,    only: nlte_fomichev_init
    use phys_control,     only: phys_getopts

    ! args

    real(r8),          intent(in) :: pref_mid(pver) ! mid point reference pressure
    ! local vars
    real(r8) :: co2_mw, o1_mw, o2_mw, o3_mw, no_mw, n2_mw ! molecular weights

    real(r8) :: delta_merge_lw      ! range of merge region
    real(r8) :: midpoint_lw         ! midpoint of merge region
    real(r8) :: psh(pver)           ! pressure scale height
    integer  :: k
    logical :: camrt

    character(len=16) :: rad_pkg
    logical :: history_scwaccm_forcing
    logical :: history_waccm
    logical :: nlte_limit_co2

!-----------------------------------------------------------------------


    call phys_getopts(radiation_scheme_out=rad_pkg, &
                      history_waccm_out=history_waccm, &
                      history_scwaccm_forcing_out=history_scwaccm_forcing)
    camrt = rad_pkg == 'CAMRT' .or. rad_pkg == 'camrt'

    ! set max/min pressures for merging regions.

    if (camrt) then
       min_pressure_lw = 1e5_r8*exp(-10._r8)
       max_pressure_lw = 1e5_r8*exp(-8.57_r8)
    else
       min_pressure_lw = p_top_for_equil_rad
       max_pressure_lw = 50._r8
    endif

    delta_merge_lw = max_pressure_lw - min_pressure_lw

    midpoint_lw = (max_pressure_lw + min_pressure_lw)/2._r8

    do k=1,pver

       ! pressure scale heights for camrt merging (cam4)
       psh(k)=log(1e5_r8/pref_mid(k))

       if ( pref_mid(k) .le. min_pressure_lw  ) then
          qrl_wt(k)= 0._r8
       else if( pref_mid(k) .ge. max_pressure_lw) then
          qrl_wt(k)= 1._r8
       else
          if (camrt) then
             ! camrt
             qrl_wt(k) = 1._r8 - tanh( (psh(k) - 8.57_r8) / 0.71_r8 )
          else
             ! rrtmg
             qrl_wt(k) = 0.5_r8 + 1.5_r8*((pref_mid(k)-midpoint_lw)/delta_merge_lw) &
                       - 2._r8*((pref_mid(k)-midpoint_lw)/delta_merge_lw)**3._r8
          endif
       endif

    end do
    
    co2_mw = mwco2
    o1_mw  = 16._r8
    o2_mw  = 32._r8
    o3_mw  = 48._r8
    no_mw  = 30._r8
    n2_mw  = 28._r8
    nlte_limit_co2 = .true.
! Initialize Fomichev parameterization
    call nlte_fomichev_init (co2_mw, n2_mw, o1_mw, o2_mw, o3_mw, no_mw, nlte_limit_co2)
    

    ! determine upppermost level that is purely solar heating (no MLT chem heating)
    ntop_qrs_cam = 0
    do k=pver,1,-1
       if (qrs_wt(k)==1._r8) ntop_qrs_cam = k
    enddo

  end subroutine radheat_init

!================================================================================================

  subroutine radheat_timestep_init (state, pbuf2d)

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in):: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)



  end subroutine radheat_timestep_init

!================================================================================================

  subroutine radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
       fsnt, flns, flnt, asdir, net_flx)
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!
! This routine provides the waccm hook for computing nonLTE cooling and
! eUV heating.
!-----------------------------------------------------------------------

#if ( defined OFFLINE_DYN )
    use metdata, only: met_rlx, met_srf_feedback
#endif
    use cam_history,           only: outfld
    use nlte_fomichev,         only: nlte_fomichev_calc

    use physics_buffer,        only : physics_buffer_desc
    use constituents,          only: cnst_get_ind
    use calculate_net_heating, only: calculate_net_heating_run
    use cam_abortutils,        only: endrun

    use rad_constituents,      only: rad_cnst_get_gas

! Arguments
    type(physics_state), intent(in)  :: state             ! Physics state variables

    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_ptend), intent(out) :: ptend             ! individual parameterization tendencies
    real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating
    real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating
    real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
    real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
    real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
    real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
    real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
    real(r8),            intent(out) :: net_flx(pcols)

! Local variables
    integer  :: k 
    integer  :: ncol                                ! number of atmospheric columns
    integer  :: lchnk                               ! chunk identifier
    real(r8) :: qrl_mrg(pcols,pver)                 ! merged LW heating
    real(r8) :: qrl_mlt(pcols,pver)                 ! M/LT longwave heating rates
    real(r8) :: qrs_mrg(pcols,pver)                 ! merged SW heating

    real(r8) :: qrlfomichev(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: o3cool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: co2cool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: c2scool(pcols,pver) ! Fomichev cooling rate ! (K/s)
    real(r8) :: xco2mmr(pcols,pver) ! CO2
    real(r8) :: xo2mmr(pcols,pver)  ! O2
    real(r8) :: xo3mmr(pcols,pver)  ! O3
    real(r8) :: xommr(pcols,pver)   ! O
    real(r8) :: xn2mmr(pcols,pver)  ! N2

    real(r8), parameter :: N2_VMR       = 0.78084_r8 ! From US standard atmosphere
    real(r8), parameter :: N2_mass      = 28.0134_r8 ! g/mol
    real(r8), parameter :: mass_dry_air = 28.9647_r8 ! g/mol

    integer  :: icall
    real(r8), pointer  :: gas_mmr(:,:)
    character(len=512) :: errmsg
    integer            :: errflg

!-----------------------------------------------------------------------

    ncol  = state%ncol
    lchnk = state%lchnk
    call physics_ptend_init(ptend, state%psetcols, 'radheat', ls=.true.)
 
! Shortwave heating is RRTMG solar heating only
    qrs_mrg(:,:) = qrs(:,:)

   icall = 0

   xommr(:pcols,:pver) = 0._r8

   call rad_cnst_get_gas(icall,'O2   ', state, pbuf, gas_mmr)
   xo2mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

   call rad_cnst_get_gas(icall,'O3   ', state, pbuf, gas_mmr)
   xo3mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

! Using standard US N2_VMR(0.78084) converted to N2_MMR
   xn2mmr(:pcols,:pver) =  N2_VMR * (N2_mass / mass_dry_air)

   call rad_cnst_get_gas(icall,'CO2  ', state, pbuf, gas_mmr)
   xco2mmr(:pcols,:pver) = gas_mmr(:pcols,:pver)
   nullify(gas_mmr)

   call nlte_fomichev_calc (lchnk,ncol,state%pmid,state%pint,state%t, &
        xo2mmr,xommr,xo3mmr,xn2mmr,xco2mmr,qrlfomichev,co2cool,o3cool,c2scool)
   qrl_mlt = qrlfomichev

!  Merge cam long wave heating for lower atmosphere with M/LT (nlte) heating
   call merge_qrl (ncol, qrl, qrl_mlt, qrl_mrg)

   ! REMOVECAM no longer need once CAM is retired and pcols doesn't exist
   net_flx = 0._r8
   ptend%s = 0._r8
   ! END_REMOVECAM

#if ( defined OFFLINE_DYN )
   do k = 1,pver
     if (met_rlx(k) < 1._r8 .or. met_srf_feedback) then
       ptend%s(:ncol,k) = (qrs_mrg(:ncol,k) + qrl_mrg(:ncol,k))
     endif
   enddo
   call calculate_net_heating_run(ncol, ptend%s(:ncol,:), qrl_mrg(:ncol,:), qrs_mrg(:ncol,:), .true., &
           fsns(:ncol), fsnt(:ncol), flns(:ncol), flnt(:ncol), net_flx(:ncol), &
           errmsg, errflg)
#else
   call calculate_net_heating_run(ncol, ptend%s(:ncol,:), qrl_mrg(:ncol,:), qrs_mrg(:ncol,:), .false., &
           fsns(:ncol), fsnt(:ncol), flns(:ncol), flnt(:ncol), net_flx(:ncol), &
           errmsg, errflg)
#endif

  end subroutine radheat_tend

!================================================================================================
  subroutine radheat_disable_waccm()
  end subroutine radheat_disable_waccm

  subroutine merge_qrl (ncol, hcam, hmlt, hmrg)
!
!  Merges long wave heating rates
!
!-----------------Input arguments-----------------------------------
    integer ncol

    real(r8), intent(in)  :: hmlt(pcols,pver)               ! Upper atmosphere heating rates
    real(r8), intent(in)  :: hcam(pcols,pver)               ! CAM heating rate
    real(r8), intent(out) :: hmrg(pcols,pver)               ! merged heating rates

!-----------------Local workspace------------------------------------

    integer k

!--------------------------------------------------------------------

    do k = 1, pver
       hmrg(:ncol,k) = qrl_wt(k) * hcam(:ncol,k) + (1._r8-qrl_wt(k)) * hmlt(:ncol,k)
    end do

  end subroutine merge_qrl

end module radheat
