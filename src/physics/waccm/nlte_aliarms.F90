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

  real(r8) :: max_pressure_aliarms   ! max_pressure_lw scaled bar

contains

!-----------------------------------------------------------------
  subroutine nlte_aliarms_init(max_pressure_lw)
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use cam_history,  only: addfld

  real(r8), intent(in)  :: max_pressure_lw  ! Pa

  if (masterproc) then
    write(iulog,*) 'init: ALI-ARMS non-LTE code'
  end if

  call addfld ('ALIARMS_Q',(/ 'lev' /), 'A','K/s','Non-LTE LW CO2 heating')

  ! Scale the max_pressure_aliarms to bar
  max_pressure_aliarms = max_pressure_lw * 1.e-05_r8

  end subroutine nlte_aliarms_init

!-----------------------------------------------------------------
  subroutine nlte_aliarms_calc (lchnk,ncol,state_zm,pmid,t,xo2,xo,xn2,xco2,cool)
!-----------------------------------------------------------------
!
!
!-----------------------------------------------------------------

  use cam_history,  only: outfld
  use iso_c_binding, only: c_float, c_int

! Input variables
  integer, intent(in) :: ncol                          ! number of atmospheric columns
  integer, intent(in) :: lchnk                         ! chunk identifier

  real(r8), intent(in) :: state_zm(pcols,pver)         ! model height
  real(r8), intent(in) :: pmid(pcols,pver)             ! model pressure at mid-point
  real(r8), intent(in) :: t(pcols,pver)                ! Neutral temperature (K)
  real(r8), intent(in) :: xco2(pcols,pver)             ! CO2 volume mixing ratio profile
  real(r8), intent(in) :: xn2(pcols,pver)              ! N2 volume mixing ratio profile
  real(r8), intent(in) :: xo(pcols,pver)               ! O volume mixing ratio profile
  real(r8), intent(in) :: xo2(pcols,pver)              ! O2 volume mixing ratio profile

! Output variables
  real(r8), intent(out) :: cool(pcols,pver)            ! CO2 NLTE cooling rate

! local variables

  real(c_float), dimension(pver) :: p, tn, zkm
  real(c_float), dimension(pver) :: co2_vmr, o_vmr, n2_vmr, o2_vmr
  real(c_float), dimension(pver) :: ali_cool

  integer(c_int) :: pver_c

  integer :: icol, iver

  ! Interface to ali C routine
  interface
     subroutine ali_(zkm, p, tn, co2_vmr, o_vmr, n2_vmr, o2_vmr, ali_cool, pver_c) bind(c,name='ali_')
        use iso_c_binding, only: c_float, c_int
        real(c_float), dimension(*) :: p, tn, zkm
        real(c_float), dimension(*) :: co2_vmr, o_vmr, n2_vmr, o2_vmr
        real(c_float), dimension(*) :: ali_cool
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

      co2_vmr(:pver_c) = xco2(icol,:pver_c)
      o_vmr(:pver_c) = xo(icol,:pver_c)
      n2_vmr(:pver_c) = xn2(icol,:pver_c)
      o2_vmr(:pver_c) = xo2(icol,:pver_c)

      call ali(zkm, p, tn, co2_vmr, o_vmr, n2_vmr, o2_vmr, ali_cool, pver_c)

      cool(icol,:pver_c) = ali_cool(:pver_c)

  enddo

  call outfld ('ALIARMS_Q', cool, pcols, lchnk)

  end subroutine nlte_aliarms_calc

  end module nlte_aliarms

