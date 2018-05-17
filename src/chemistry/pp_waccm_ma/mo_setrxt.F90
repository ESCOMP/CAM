
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,91) = 8.00e-14_r8
      rate(:,:,92) = 3.90e-17_r8
      rate(:,:,95) = 4.20e-13_r8
      rate(:,:,96) = 8.50e-2_r8
      rate(:,:,97) = 1.30e-16_r8
      rate(:,:,99) = 1.00e-20_r8
      rate(:,:,100) = 2.58e-04_r8
      rate(:,:,107) = 1.20e-10_r8
      rate(:,:,108) = 2.02e-10_r8
      rate(:,:,109) = 1.204e-10_r8
      rate(:,:,110) = 1.50e-10_r8
      rate(:,:,111) = 9.75e-11_r8
      rate(:,:,112) = 1.50e-11_r8
      rate(:,:,113) = 7.20e-11_r8
      rate(:,:,114) = 1.794e-10_r8
      rate(:,:,115) = 1.628e-10_r8
      rate(:,:,116) = 2.84e-10_r8
      rate(:,:,117) = 1.674e-10_r8
      rate(:,:,118) = 9.60e-11_r8
      rate(:,:,119) = 4.10e-11_r8
      rate(:,:,120) = 1.012e-10_r8
      rate(:,:,121) = 1.20e-10_r8
      rate(:,:,122) = 4.49e-10_r8
      rate(:,:,123) = 2.57e-10_r8
      rate(:,:,124) = 2.14e-11_r8
      rate(:,:,125) = 1.90e-10_r8
      rate(:,:,126) = 1.31e-10_r8
      rate(:,:,127) = 3.50e-11_r8
      rate(:,:,128) = 9.00e-12_r8
      rate(:,:,129) = 1.20e-10_r8
      rate(:,:,130) = 1.50e-10_r8
      rate(:,:,131) = 1.20e-10_r8
      rate(:,:,134) = 7.20e-11_r8
      rate(:,:,135) = 6.90e-12_r8
      rate(:,:,136) = 1.60e-12_r8
      rate(:,:,140) = 1.80e-12_r8
      rate(:,:,143) = 1.80e-12_r8
      rate(:,:,149) = 5.00e-12_r8
      rate(:,:,150) = 7.00e-13_r8
      rate(:,:,151) = 5.00e-11_r8
      rate(:,:,168) = 1.00e-11_r8
      rate(:,:,169) = 2.20e-11_r8
      rate(:,:,170) = 3.50e-12_r8
      rate(:,:,195) = 1.70e-13_r8
      rate(:,:,267) = 9.0e-10_r8
      rate(:,:,268) = 1.0e-10_r8
      rate(:,:,269) = 4.4e-10_r8
      rate(:,:,270) = 4.0e-10_r8
      rate(:,:,271) = 2.0e-10_r8
      rate(:,:,272) = 1.0e-12_r8
      rate(:,:,273) = 6.0e-11_r8
      rate(:,:,274) = 5.0e-16_r8
      rate(:,:,278) = 4.8e-10_r8
      rate(:,:,279) = 1.0e-10_r8
      rate(:,:,280) = 4.0e-10_r8
      rate(:,:,283) = 5.0e-12_r8
      rate(:,:,284) = 7.0e-10_r8
      rate(:,:,285) = 8.0e-10_r8
      rate(:,:,287) = 4.7e-2_r8
      rate(:,:,288) = 1.71e-1_r8
      rate(:,:,289) = 7.7e-5_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,89) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,93) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,94) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,98) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,101) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,102) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,103) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,104) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,105) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,106) = 4.63e-11_r8 * exp_fac(:,:)
      rate(:,:,133) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,137) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,138) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,139) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,205) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,142) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,144) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,145) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,213) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,241) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,146) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,148) = 1.40e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      rate(:,:,152) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,153) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,154) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,155) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,156) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,158) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,177) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,182) = 7.40e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,159) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,214) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,160) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,162) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,188) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,167) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,172) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,174) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,175) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,176) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,178) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,179) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,180) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,181) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,183) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,204) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,212) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,184) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,186) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,211) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,185) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,189) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,190) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,193) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,194) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,196) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,197) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,218) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,198) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,229) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,199) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,200) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,201) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,202) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,203) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,231) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,206) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,207) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,209) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,215) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,216) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,217) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,267) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,268) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,269) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,270) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,271) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,272) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,273) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,274) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,278) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,279) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,280) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,283) = 5.0e-12_r8 * exp_fac(:,:)
      rate(:,:,284) = 7.0e-10_r8 * exp_fac(:,:)
      rate(:,:,285) = 8.0e-10_r8 * exp_fac(:,:)
      rate(:,:,287) = 4.7e-2_r8 * exp_fac(:,:)
      rate(:,:,288) = 1.71e-1_r8 * exp_fac(:,:)
      rate(:,:,289) = 7.7e-5_r8 * exp_fac(:,:)
      rate(:,:,219) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      rate(:,:,220) = 6.00e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,221) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,222) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,223) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,224) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,227) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,238) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,225) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,226) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,228) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,230) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,232) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,233) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,236) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,237) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,239) = 2.80e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,240) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,132), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,141), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,157), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,161), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,163), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,165), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,171), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,187), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,191), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,208), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,235), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,91) = 8.00e-14_r8
      rate(:,:kbot,92) = 3.90e-17_r8
      rate(:,:kbot,97) = 1.30e-16_r8
      rate(:,:kbot,99) = 1.00e-20_r8
      rate(:,:kbot,135) = 6.90e-12_r8
      rate(:,:kbot,149) = 5.00e-12_r8
      rate(:,:kbot,150) = 7.00e-13_r8
      rate(:,:kbot,268) = 1.0e-10_r8
      rate(:,:kbot,269) = 4.4e-10_r8
      rate(:,:kbot,270) = 4.0e-10_r8
      rate(:,:kbot,271) = 2.0e-10_r8
      rate(:,:kbot,272) = 1.0e-12_r8
      rate(:,:kbot,273) = 6.0e-11_r8
      rate(:,:kbot,278) = 4.8e-10_r8
      rate(:,:kbot,279) = 1.0e-10_r8
      rate(:,:kbot,280) = 4.0e-10_r8
      rate(:,:kbot,283) = 5.0e-12_r8
      rate(:,:kbot,284) = 7.0e-10_r8
      rate(:,:kbot,285) = 8.0e-10_r8
      rate(:,:kbot,287) = 4.7e-2_r8
      rate(:,:kbot,288) = 1.71e-1_r8
      rate(:,:kbot,289) = 7.7e-5_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,89) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,93) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,94) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,98) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,101) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,102) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,103) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,133) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,137) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,138) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,139) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,145) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,146) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,152) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,153) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,158) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,159) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,160) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,132) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
