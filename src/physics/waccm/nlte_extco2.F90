!--------------------------------------------------------------------------------
! Ref:
!  Lopez-Puertas, M., Fabiano, F., Fomichev, V., Funke, B., and Marsh, D. R.:
!  An improved and extended parameterization of the CO2 15 um cooling in the middle
!  and upper atmosphere (CO2_cool_fort-1.0), Geosci. Model Dev., 17, 4401-4432,
!  https://doi.org/10.5194/gmd-17-4401-2024, 2024.
!--------------------------------------------------------------------------------
module nlte_extco2

  use ppgrid,             only: pcols, pver
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_logfile,        only: iulog
  use spmd_utils,         only: masterproc
  use cam_history,  only: add_default, addfld, outfld

  use co2cool, only: co2_nlte_cool

  implicit none

  private
  public :: nlte_extco2_init
  public :: nlte_extco2_hrate

  real(r8) :: o1_mw_inv  = -huge(1.0_r8)       ! O molecular weight (inverse)
  real(r8) :: o2_mw_inv  = -huge(1.0_r8)       ! O2 molecular weight (inverse)
  real(r8) :: co2_mw_inv = -huge(1.0_r8)       ! CO2 molecular weight (inverse)
  real(r8) :: n2_mw_inv  = -huge(1.0_r8)       ! N2 molecular weight (inverse)

  integer :: lev0 = -1

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine nlte_extco2_init(max_pressure_lw,co2_mw,n2_mw,o1_mw,o2_mw)
    use ref_pres, only: pref_mid

    ! Input variables
    real(r8), intent(in)  :: max_pressure_lw  ! Pa
    real(r8), intent(in) :: o1_mw             ! O molecular weight
    real(r8), intent(in) :: o2_mw             ! O2 molecular weight
    real(r8), intent(in) :: co2_mw            ! CO2 molecular weight
    real(r8), intent(in) :: n2_mw             ! N2 molecular weight

    co2_mw_inv = 1._r8/co2_mw
    o1_mw_inv  = 1._r8/o1_mw
    o2_mw_inv  = 1._r8/o2_mw
    n2_mw_inv  = 1._r8/n2_mw

    lev0 = 1
!    do while (pref_mid(lev0) < max_pressure_lw)
!       lev0 = lev0+1
!    end do
    lev0 = pver

    call addfld ('QCO2ext',   (/ 'lev' /), 'A','K/s','Extended CO2 cooling')
    call addfld ('TCO2ext',   (/ 'lev' /), 'A','K','Temp used in Extended CO2 cooling')
    call addfld ('PCO2ext',   (/ 'lev' /), 'A','hPa','Press used in Extended CO2 cooling')
    call addfld ('CO2_ext',   (/ 'lev' /), 'A','mol/mol','CO2 vmr used in Extended CO2 cooling')
    call addfld ('N2_ext',    (/ 'lev' /), 'A','mol/mol','N2 vmr used in Extended CO2 cooling')
    call addfld ('O_ext',     (/ 'lev' /), 'A','mol/mol','O vmr used in Extended CO2 cooling')
    call addfld ('O2_ext',    (/ 'lev' /), 'A','mol/mol','O2 vmr used in Extended CO2 cooling')

  end subroutine nlte_extco2_init


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine nlte_extco2_hrate(lchnk, ncol, temp, pres, co2mmr, n2mmr, ommr, o2mmr, &
                               co2cooling)
    use air_composition, only: mbarv

    ! Input variables
    integer, intent(in) :: ncol                          ! number of atmospheric columns
    integer, intent(in) :: lchnk                         ! chunk identifier

    real(r8), intent(in) :: temp(:,:)          ! K
    real(r8), intent(in) :: pres(:,:)          ! Pa
    real(r8), intent(in) :: co2mmr(:,:)        ! CO2 mass mixing ratio profile
    real(r8), intent(in) :: n2mmr(:,:)         ! N2 mass mixing ratio profile
    real(r8), intent(in) :: ommr(:,:)          ! O mass mixing ratio profile
    real(r8), intent(in) :: o2mmr(:,:)         ! O2 mass mixing ratio profile

    ! Output
    real(r8), intent(out) :: co2cooling(:,:)   ! K sec-1

    ! Local vars
    real(r8) :: co2vmr(pver)
    real(r8) :: ovmr(pver)
    real(r8) :: n2vmr(pver)
    real(r8) :: o2vmr(pver)
    real(r8) :: heatrate(pver) ! K/day
    real(r8) :: hPa(pver)

    real(r8) :: surf_temp ! K
    integer :: icol

    real(r8) :: tempout(pcols,pver)
    real(r8) :: presout(pcols,pver)
    real(r8) :: co2_out(pcols,pver)
    real(r8) :: n2_out(pcols,pver)
    real(r8) :: o_out(pcols,pver)
    real(r8) :: o2_out(pcols,pver)

    real(r8), parameter :: day_per_sec = 1._r8/86400._r8

    co2cooling = 0._r8

    do icol=1,ncol

       ! Convert to VMR from mmr
       co2vmr(:) = mbarv(icol,:,lchnk) * co2mmr(icol,:) * co2_mw_inv
       ovmr(:)   = mbarv(icol,:,lchnk) * ommr(icol,:) * o1_mw_inv
       n2vmr(:)  = mbarv(icol,:,lchnk) * n2mmr(icol,:) * n2_mw_inv
       o2vmr(:)  = mbarv(icol,:,lchnk) * o2mmr(icol,:) * o2_mw_inv
       hPa(:)    = pres(icol,:) * 1.e-2_r8 ! Pa --> hPa
       surf_temp = temp(icol,pver)

       tempout(icol,:) = temp(icol,:)
       presout(icol,:) = hPa(:)
       co2_out(icol,:) = co2vmr(:)
       o_out(icol,:)   = ovmr(:)
       o2_out(icol,:)  = o2vmr(:)
       n2_out(icol,:)  = n2vmr(:)

       heatrate(:) = 0._r8

       ! Units: Temperature in K, pressure in hPa, vmrs in mol/mol, heating rate in K/day
       call co2_nlte_cool(temp(icol,:), hPa, co2vmr, ovmr, o2vmr, n2vmr, lev0, &
                          surf_temp, heatrate) !   (K day-1)

       co2cooling(icol,:) = heatrate(:) * day_per_sec ! K day-1 --> K sec-1

    end do

    call outfld ('QCO2ext', co2cooling(:ncol,:), ncol, lchnk)
    call outfld ('TCO2ext', tempout(:ncol,:), ncol, lchnk)
    call outfld ('PCO2ext', presout(:ncol,:), ncol, lchnk)
    call outfld ('CO2_ext', co2_out(:ncol,:), ncol, lchnk)
    call outfld ('O_ext',   o_out(:ncol,:),   ncol, lchnk)
    call outfld ('O2_ext',  o2_out(:ncol,:),  ncol, lchnk)
    call outfld ('N2_ext',  n2_out(:ncol,:),  ncol, lchnk)

  end subroutine nlte_extco2_hrate



end module nlte_extco2
