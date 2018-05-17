
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

      rate(:,:,121) = 8.00e-14_r8
      rate(:,:,122) = 3.90e-17_r8
      rate(:,:,125) = 4.20e-13_r8
      rate(:,:,126) = 8.50e-2_r8
      rate(:,:,127) = 1.30e-16_r8
      rate(:,:,129) = 1.00e-20_r8
      rate(:,:,130) = 2.58e-04_r8
      rate(:,:,137) = 1.20e-10_r8
      rate(:,:,138) = 2.02e-10_r8
      rate(:,:,139) = 1.204e-10_r8
      rate(:,:,140) = 1.50e-10_r8
      rate(:,:,141) = 9.75e-11_r8
      rate(:,:,142) = 1.50e-11_r8
      rate(:,:,143) = 7.20e-11_r8
      rate(:,:,144) = 1.794e-10_r8
      rate(:,:,145) = 1.628e-10_r8
      rate(:,:,146) = 2.84e-10_r8
      rate(:,:,147) = 1.674e-10_r8
      rate(:,:,148) = 9.60e-11_r8
      rate(:,:,149) = 4.10e-11_r8
      rate(:,:,150) = 1.012e-10_r8
      rate(:,:,151) = 1.20e-10_r8
      rate(:,:,152) = 4.49e-10_r8
      rate(:,:,153) = 2.57e-10_r8
      rate(:,:,154) = 2.14e-11_r8
      rate(:,:,155) = 1.90e-10_r8
      rate(:,:,156) = 1.31e-10_r8
      rate(:,:,157) = 3.50e-11_r8
      rate(:,:,158) = 9.00e-12_r8
      rate(:,:,159) = 1.20e-10_r8
      rate(:,:,160) = 1.50e-10_r8
      rate(:,:,161) = 1.20e-10_r8
      rate(:,:,165) = 7.20e-11_r8
      rate(:,:,166) = 6.90e-12_r8
      rate(:,:,167) = 1.60e-12_r8
      rate(:,:,171) = 1.80e-12_r8
      rate(:,:,174) = 1.80e-12_r8
      rate(:,:,182) = 5.00e-12_r8
      rate(:,:,183) = 7.00e-13_r8
      rate(:,:,184) = 5.00e-11_r8
      rate(:,:,201) = 1.00e-11_r8
      rate(:,:,202) = 2.20e-11_r8
      rate(:,:,203) = 3.50e-12_r8
      rate(:,:,228) = 1.70e-13_r8
      rate(:,:,279) = 4.50e-13_r8
      rate(:,:,291) = 1.00e-14_r8
      rate(:,:,294) = 7.00e-13_r8
      rate(:,:,297) = 2.00e-13_r8
      rate(:,:,298) = 6.80e-14_r8
      rate(:,:,307) = 1.00e-12_r8
      rate(:,:,308) = 1.00e-11_r8
      rate(:,:,309) = 1.15e-11_r8
      rate(:,:,312) = 4.00e-14_r8
      rate(:,:,329) = 3.00e-12_r8
      rate(:,:,332) = 6.80e-13_r8
      rate(:,:,333) = 5.40e-11_r8
      rate(:,:,345) = 2.40e-12_r8
      rate(:,:,348) = 1.40e-11_r8
      rate(:,:,351) = 5.00e-12_r8
      rate(:,:,363) = 2.40e-12_r8
      rate(:,:,367) = 1.40e-11_r8
      rate(:,:,369) = 2.40e-12_r8
      rate(:,:,371) = 3.50e-12_r8
      rate(:,:,372) = 4.50e-11_r8
      rate(:,:,379) = 2.40e-12_r8
      rate(:,:,389) = 3.00e-12_r8
      rate(:,:,390) = 1.00e-11_r8
      rate(:,:,394) = 2.3e-11_r8
      rate(:,:,406) = 7.10e-6_r8
      rate(:,:,407) = 7.10e-6_r8
      rate(:,:,409) = 6.34e-8_r8
      rate(:,:,410) = 6.34e-8_r8
      rate(:,:,411) = 6.34e-8_r8
      rate(:,:,412) = 6.34e-8_r8
      rate(:,:,413) = 6.34e-8_r8
      rate(:,:,414) = 6.34e-8_r8
      rate(:,:,415) = 6.34e-8_r8
      rate(:,:,416) = 6.34e-8_r8
      rate(:,:,417) = 6.34e-8_r8
      rate(:,:,418) = 6.34e-8_r8
      rate(:,:,419) = 6.34e-8_r8
      rate(:,:,420) = 6.34e-8_r8
      rate(:,:,421) = 6.34e-8_r8
      rate(:,:,422) = 6.34e-8_r8
      rate(:,:,423) = 6.34e-8_r8
      rate(:,:,424) = 6.34e-8_r8
      rate(:,:,425) = 6.34e-8_r8
      rate(:,:,426) = 6.34e-8_r8
      rate(:,:,427) = 6.34e-8_r8
      rate(:,:,428) = 6.34e-8_r8
      rate(:,:,431) = 6.60E-11_r8
      rate(:,:,432) = 2.30E-12_r8
      rate(:,:,433) = 1.20E-11_r8
      rate(:,:,437) = 1.40E-11_r8
      rate(:,:,438) = 2.80E-11_r8
      rate(:,:,439) = 5.70E-11_r8
      rate(:,:,440) = 1.90E-12_r8
      rate(:,:,467) = 9.0e-10_r8
      rate(:,:,468) = 1.0e-10_r8
      rate(:,:,469) = 4.4e-10_r8
      rate(:,:,470) = 4.0e-10_r8
      rate(:,:,471) = 2.0e-10_r8
      rate(:,:,472) = 1.0e-12_r8
      rate(:,:,473) = 6.0e-11_r8
      rate(:,:,474) = 5.0e-16_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,119) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,123) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,124) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,128) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,131) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,132) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,133) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,134) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,135) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,136) = 4.63e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,162) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,186) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,164) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,168) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,289) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,316) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,321) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,334) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,338) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,362) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,375) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,386) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,400) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,169) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,170) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,238) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,173) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,175) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,176) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,246) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,278) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,299) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,319) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,323) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,328) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,340) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,349) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,365) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,377) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,388) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,402) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,177) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,179) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,360) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,181) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      rate(:,:,185) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,187) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,188) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,189) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,191) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,210) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,215) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,302) = 8.10e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,192) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,247) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,193) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,195) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,221) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,200) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,205) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,207) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,208) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,209) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,211) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,212) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,213) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,214) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,216) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,237) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,245) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,217) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,219) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,244) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,218) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,222) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,223) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,226) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,227) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,229) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,230) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,251) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,231) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,262) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,232) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,233) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,234) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,235) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,236) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,264) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,239) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,240) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,243) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,242) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,248) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,249) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,250) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,279) = 4.50e-13_r8 * exp_fac(:,:)
      rate(:,:,291) = 1.00e-14_r8 * exp_fac(:,:)
      rate(:,:,294) = 7.00e-13_r8 * exp_fac(:,:)
      rate(:,:,297) = 2.00e-13_r8 * exp_fac(:,:)
      rate(:,:,298) = 6.80e-14_r8 * exp_fac(:,:)
      rate(:,:,307) = 1.00e-12_r8 * exp_fac(:,:)
      rate(:,:,308) = 1.00e-11_r8 * exp_fac(:,:)
      rate(:,:,309) = 1.15e-11_r8 * exp_fac(:,:)
      rate(:,:,312) = 4.00e-14_r8 * exp_fac(:,:)
      rate(:,:,329) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,332) = 6.80e-13_r8 * exp_fac(:,:)
      rate(:,:,333) = 5.40e-11_r8 * exp_fac(:,:)
      rate(:,:,345) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,348) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,351) = 5.00e-12_r8 * exp_fac(:,:)
      rate(:,:,363) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,367) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,369) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,371) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,372) = 4.50e-11_r8 * exp_fac(:,:)
      rate(:,:,379) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,389) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,390) = 1.00e-11_r8 * exp_fac(:,:)
      rate(:,:,394) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,406) = 7.10e-6_r8 * exp_fac(:,:)
      rate(:,:,407) = 7.10e-6_r8 * exp_fac(:,:)
      rate(:,:,409) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,410) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,411) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,412) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,413) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,414) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,415) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,416) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,417) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,418) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,419) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,420) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,421) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,422) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,423) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,424) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,425) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,426) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,427) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,428) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,431) = 6.60E-11_r8 * exp_fac(:,:)
      rate(:,:,432) = 2.30E-12_r8 * exp_fac(:,:)
      rate(:,:,433) = 1.20E-11_r8 * exp_fac(:,:)
      rate(:,:,437) = 1.40E-11_r8 * exp_fac(:,:)
      rate(:,:,438) = 2.80E-11_r8 * exp_fac(:,:)
      rate(:,:,439) = 5.70E-11_r8 * exp_fac(:,:)
      rate(:,:,440) = 1.90E-12_r8 * exp_fac(:,:)
      rate(:,:,467) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,468) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,469) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,470) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,471) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,472) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,473) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,474) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,252) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,253) = 6.00e-12_r8 * exp_fac(:,:)
      rate(:,:,347) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,366) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,381) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,254) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,255) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,256) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,257) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,260) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,271) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,258) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,259) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,261) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,263) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,265) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,266) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,269) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,270) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,272) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,273) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,325) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,274) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,275) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,276) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,277) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,280) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,281) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,282) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,290) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,296) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,317) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,322) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,326) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,339) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,346) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,364) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,370) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,376) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,380) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,387) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,392) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,395) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,401) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,285) = 7.20e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,287) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,292) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,293) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,295) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 350._r8 * itemp(:,:) )
      rate(:,:,300) = 4.63e-12_r8 * exp_fac(:,:)
      rate(:,:,393) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,396) = 2.6e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,301) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,314) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,304) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,352) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,305) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,306) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,327) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,353) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,310) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,315) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,318) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,320) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,330) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,331) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,373) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,335) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,336) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,337) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,341) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,374) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,342) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,343) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,344) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,350) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,368) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,378) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,354) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,355) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,359) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,361) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,382) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,383) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,385) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,391) = 2.3e-12_r8 * exp( -193._r8 * itemp(:,:) )
      rate(:,:,397) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,398) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,399) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,429) = 2.10E-11_r8 * exp( -2200.0_r8 * itemp(:,:) )
      rate(:,:,430) = 1.10E-13_r8 * exp( -1200.0_r8 * itemp(:,:) )
      rate(:,:,434) = 2.70E-11_r8 * exp( 335._r8 * itemp(:,:) )
      rate(:,:,435) = 1.25E-13_r8 * exp( -2190.0_r8 * itemp(:,:) )
      rate(:,:,436) = 3.40E-12_r8 * exp( -1100.0_r8 * itemp(:,:) )
      rate(:,:,443) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,445) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,163), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,172), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,180), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,190), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,194), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,196), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,198), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,204), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,220), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,224), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,241), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,268), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-30_r8 * itemp(:,:)**2.4_r8
      kinf(:,:) = 2.2e-10_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,283), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-29_r8 * itemp(:,:)**3.3_r8
      kinf(:,:) = 3.1e-10_r8 * itemp(:,:)
      call jpl( rate(1,1,284), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,286), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,288), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,303), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,313), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,358), m, 0.5_r8, ko, kinf, n )

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

      rate(:,:kbot,121) = 8.00e-14_r8
      rate(:,:kbot,122) = 3.90e-17_r8
      rate(:,:kbot,127) = 1.30e-16_r8
      rate(:,:kbot,129) = 1.00e-20_r8
      rate(:,:kbot,166) = 6.90e-12_r8
      rate(:,:kbot,182) = 5.00e-12_r8
      rate(:,:kbot,183) = 7.00e-13_r8
      rate(:,:kbot,468) = 1.0e-10_r8
      rate(:,:kbot,469) = 4.4e-10_r8
      rate(:,:kbot,470) = 4.0e-10_r8
      rate(:,:kbot,471) = 2.0e-10_r8
      rate(:,:kbot,472) = 1.0e-12_r8
      rate(:,:kbot,473) = 6.0e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,119) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,123) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,124) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,128) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,131) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,132) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,133) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,164) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,168) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,169) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,170) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,176) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,177) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,185) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,186) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,191) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,192) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,193) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,163) = wrk(:,:)



















      end subroutine setrxt_hrates

      end module mo_setrxt
