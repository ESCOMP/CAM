module nlte_aliarms

!
! provides calculation of non-LTE heating rates by ALI-ARMS non-LTE code
!
  use ppgrid,             only: pcols, pver
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use cam_logfile,        only: iulog
  use spmd_utils,         only: masterproc

  implicit none
  private
  save

! Public interfaces
  public &
     nlte_aliarms_init, &
     nlte_aliarms_calc

  real(r8) :: max_pressure_aliarms = -huge(1.0_r8)  ! max_pressure_lw scaled bar

  real(r8) :: o1_mw_inv                        ! O molecular weight (inverse)
  real(r8) :: o2_mw_inv                        ! O2 molecular weight (inverse)
  real(r8) :: co2_mw_inv                       ! CO2 molecular weight (inverse)
  real(r8) :: n2_mw_inv                        ! N2 molecular weight (inverse)
contains

!-----------------------------------------------------------------
  subroutine nlte_aliarms_init(max_pressure_lw,co2_mw,n2_mw,o1_mw,o2_mw)
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use cam_history,  only: addfld

  real(r8), intent(in)  :: max_pressure_lw  ! Pa
  real(r8), intent(in) :: o1_mw             ! O molecular weight
  real(r8), intent(in) :: o2_mw             ! O2 molecular weight
  real(r8), intent(in) :: co2_mw            ! CO2 molecular weight
  real(r8), intent(in) :: n2_mw             ! N2 molecular weight

  if (masterproc) then
    write(iulog,*) 'init: ALI-ARMS non-LTE code'
  end if

  call addfld ('ALIARMS_Q',(/ 'lev' /), 'A','K/s','Non-LTE LW CO2 heating rate')

  ! Scale the max_pressure_aliarms to bar
  max_pressure_aliarms = max_pressure_lw * 1.e-05_r8

  co2_mw_inv = 1._r8/co2_mw
  o1_mw_inv  = 1._r8/o1_mw
  o2_mw_inv  = 1._r8/o2_mw
  n2_mw_inv  = 1._r8/n2_mw

  end subroutine nlte_aliarms_init

!-----------------------------------------------------------------
  subroutine nlte_aliarms_calc (lchnk,ncol,state_zm,pmid,t,xo2mmr,xommr,xn2mmr,xco2mmr,cool)
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use physconst,     only: mbarv
  use cam_history,   only: outfld
  use iso_c_binding, only: c_float, c_int
  use shr_infnan_mod, only: is_nan => shr_infnan_isnan
  use cam_abortutils,     only: endrun

! Input variables
  integer, intent(in) :: ncol                          ! number of atmospheric columns
  integer, intent(in) :: lchnk                         ! chunk identifier

  real(r8), intent(in) :: state_zm(pcols,pver)         ! model height (m)
  real(r8), intent(in) :: pmid(pcols,pver)             ! model pressure at mid-point (Pa)
  real(r8), intent(in) :: t(pcols,pver)                ! Neutral temperature (K)

  real(r8), intent(in) :: xco2mmr(pcols,pver)          ! CO2 mass mixing ratio profile
  real(r8), intent(in) :: xn2mmr(pcols,pver)           ! N2 mass mixing ratio profile
  real(r8), intent(in) :: xommr(pcols,pver)            ! O mass mixing ratio profile
  real(r8), intent(in) :: xo2mmr(pcols,pver)           ! O2 mass mixing ratio profile

! Output variables
  real(r8), intent(out) :: cool(pcols,pver)            ! CO2 NLTE cooling rate  (K/s)

! local variables

  real(c_float), dimension(pver) :: p, tn, zkm
  real(c_float), dimension(pver) :: co2_vmr, o_vmr, n2_vmr, o2_vmr
  real(c_float), dimension(pver) :: ali_cool

  integer(c_int) :: pver_c

  integer :: icol, iver, i, j

  character (len=160) :: errstring

  ! Interface to ali C routine
  interface
     subroutine ali_(zkm, p, tn, co2_vmr, o_vmr, n2_vmr, o2_vmr, ali_cool, pver_c) bind(c,name='ali_')
        use iso_c_binding, only: c_float, c_int
        real(c_float), dimension(*) :: p, tn, zkm                      ! (in) input pressure(Pa), temperature(K) and height (km)
        real(c_float), dimension(*) :: co2_vmr, o_vmr, n2_vmr, o2_vmr  ! (in) volume mixing ratios
        real(c_float), dimension(*) :: ali_cool                        ! (out) cooling rate (K/s)
        integer(c_int) :: pver_c
     end subroutine ali_
  end interface

  cool(:,:) = 0.0_r8

  do icol=1,ncol

      ali_cool(:) = 0.0_c_float
      zkm(:)      = 0.0_c_float
      tn(:)       = 0.0_c_float
      co2_vmr(:)  = 0.0_c_float
      o_vmr(:)    = 0.0_c_float
      n2_vmr(:)   = 0.0_c_float
      o2_vmr(:)   = 0.0_c_float

      p = pmid(icol,:)*1.0e-5_c_float ! convert pmid in Pa to bars
      pver_c=0
      do iver = 1,pver
        if (p(iver) < max_pressure_aliarms) then
          pver_c=pver_c+1
        else
          exit  ! Have gone past the maximum pressure
        end if
      end do
      zkm(:pver_c) = state_zm(icol,:pver_c)*1.e-3_c_float
      tn(:pver_c) = t(icol,:pver_c)

      ! Convert to VMR from mmr
      co2_vmr(:pver_c) = mbarv(icol,:pver_c ,lchnk) * xco2mmr(icol,:pver_c) * co2_mw_inv
      o_vmr(:pver_c) = mbarv(icol,:pver_c ,lchnk) * xommr(icol,:pver_c) * o1_mw_inv
      n2_vmr(:pver_c) = mbarv(icol,:pver_c ,lchnk) * xn2mmr(icol,:pver_c) * n2_mw_inv
      o2_vmr(:pver_c) = mbarv(icol,:pver_c ,lchnk) * xo2mmr(icol,:pver_c) * o2_mw_inv

      call ali(zkm, p, tn, co2_vmr, o_vmr, n2_vmr, o2_vmr, ali_cool, pver_c)

      cool(icol,:pver_c) = ali_cool(:pver_c)

  enddo

  ! Check the rates
  do j=1,pver
     do i=1,ncol
        if (is_nan(cool(i,j))) then
           write(errstring,*) 'nlte_aliarms_calc: Nan in qrlaliarms for chunk', lchnk
           call endrun (errstring)
        end if
     end do
  end do

  ! Check cool for any way out-of-bounds values
  if (any(cool(:ncol,:pver_c) > 0.2_r8) .or. any(cool(:ncol,:pver_c)<-1.0_r8)) then
     write(errstring,*) 'nlte_aliarms_calc: Cooling rate (cool) is greater than .2 or less than -1 K/s for chunk ', lchnk
     call endrun (errstring)
  end if

  call outfld ('ALIARMS_Q', cool, pcols, lchnk)

  end subroutine nlte_aliarms_calc

  end module nlte_aliarms

