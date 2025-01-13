module nlte_extco2

  use ppgrid,             only: pcols, pver
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_logfile,        only: iulog
  use spmd_utils,         only: masterproc

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
    do while (pref_mid(lev0) < max_pressure_lw)
       lev0 = lev0+1
    end do

  end subroutine nlte_extco2_init


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine nlte_extco2_hrate(lchnk, ncol, temp, pres, co2mmr, n2mmr, ommr, o2mmr, &
                                surf_temp, co2cooling)
    use air_composition, only: mbarv

    ! Input variables
    integer, intent(in) :: ncol                          ! number of atmospheric columns
    integer, intent(in) :: lchnk                         ! chunk identifier

    real(r8), intent(in) :: temp(pcols,pver)          ! K
    real(r8), intent(in) :: pres(pcols,pver)          ! hPa
    real(r8), intent(in) :: co2mmr(pcols,pver)        ! CO2 mass mixing ratio profile
    real(r8), intent(in) :: n2mmr(pcols,pver)         ! N2 mass mixing ratio profile
    real(r8), intent(in) :: ommr(pcols,pver)          ! O mass mixing ratio profile
    real(r8), intent(in) :: o2mmr(pcols,pver)         ! O2 mass mixing ratio profile
    real(r8), intent(in) :: surf_temp(pcols)          ! K

    ! Output
    real(r8), intent(out) :: co2cooling(pcols,pver)   ! K sec-1

    ! Local vars
    real(r8) :: co2vmr(pver)
    real(r8) :: ovmr(pver)
    real(r8) :: n2vmr(pver)
    real(r8) :: o2vmr(pver)
    real(r8) :: heatrate(pver)

    integer :: icol

    real(r8), parameter :: day_per_sec = 1._r8/86400._r8

    co2cooling = 0._r8

    do icol=1,ncol

       ! Convert to VMR from mmr
       co2vmr(:) = mbarv(icol,:,lchnk) * co2mmr(icol,:) * co2_mw_inv
       ovmr(:)   = mbarv(icol,:,lchnk) * ommr(icol,:) * o1_mw_inv
       n2vmr(:)  = mbarv(icol,:,lchnk) * n2mmr(icol,:) * n2_mw_inv
       o2vmr(:)  = mbarv(icol,:,lchnk) * o2mmr(icol,:) * o2_mw_inv

       heatrate(:) = 0._r8

       ! Units: Temperature in K, pressure in hPa, vmrs in mol/mol (not ppm), heating rate in K/day
       call co2_nlte_cool(temp(icol,:), pres(icol,:), co2vmr, ovmr, o2vmr, n2vmr, lev0, &
                          surf_temp(icol), heatrate) !   (K day-1)

       co2cooling(icol,:) = heatrate * day_per_sec ! K day-1 --> K sec-1

    end do

  end subroutine nlte_extco2_hrate



end module nlte_extco2
