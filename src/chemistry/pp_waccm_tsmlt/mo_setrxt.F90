
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

      rate(:,:,118) = 8.00e-14_r8
      rate(:,:,119) = 3.90e-17_r8
      rate(:,:,122) = 4.20e-13_r8
      rate(:,:,123) = 8.50e-2_r8
      rate(:,:,124) = 1.30e-16_r8
      rate(:,:,126) = 1.00e-20_r8
      rate(:,:,127) = 2.58e-04_r8
      rate(:,:,134) = 1.20e-10_r8
      rate(:,:,135) = 2.02e-10_r8
      rate(:,:,136) = 1.204e-10_r8
      rate(:,:,137) = 1.50e-10_r8
      rate(:,:,138) = 9.75e-11_r8
      rate(:,:,139) = 1.50e-11_r8
      rate(:,:,140) = 7.20e-11_r8
      rate(:,:,141) = 1.794e-10_r8
      rate(:,:,142) = 1.628e-10_r8
      rate(:,:,143) = 2.84e-10_r8
      rate(:,:,144) = 1.674e-10_r8
      rate(:,:,145) = 9.60e-11_r8
      rate(:,:,146) = 4.10e-11_r8
      rate(:,:,147) = 1.012e-10_r8
      rate(:,:,148) = 1.20e-10_r8
      rate(:,:,149) = 4.49e-10_r8
      rate(:,:,150) = 2.57e-10_r8
      rate(:,:,151) = 2.14e-11_r8
      rate(:,:,152) = 1.90e-10_r8
      rate(:,:,153) = 1.31e-10_r8
      rate(:,:,154) = 3.50e-11_r8
      rate(:,:,155) = 9.00e-12_r8
      rate(:,:,156) = 1.20e-10_r8
      rate(:,:,157) = 1.50e-10_r8
      rate(:,:,158) = 1.20e-10_r8
      rate(:,:,162) = 7.20e-11_r8
      rate(:,:,163) = 6.90e-12_r8
      rate(:,:,164) = 1.60e-12_r8
      rate(:,:,168) = 1.80e-12_r8
      rate(:,:,171) = 1.80e-12_r8
      rate(:,:,179) = 5.00e-12_r8
      rate(:,:,180) = 7.00e-13_r8
      rate(:,:,181) = 5.00e-11_r8
      rate(:,:,198) = 1.00e-11_r8
      rate(:,:,199) = 2.20e-11_r8
      rate(:,:,200) = 3.50e-12_r8
      rate(:,:,225) = 1.70e-13_r8
      rate(:,:,276) = 4.50e-13_r8
      rate(:,:,288) = 1.00e-14_r8
      rate(:,:,291) = 7.00e-13_r8
      rate(:,:,294) = 2.00e-13_r8
      rate(:,:,295) = 6.80e-14_r8
      rate(:,:,304) = 1.00e-12_r8
      rate(:,:,305) = 1.00e-11_r8
      rate(:,:,306) = 1.15e-11_r8
      rate(:,:,309) = 4.00e-14_r8
      rate(:,:,326) = 3.00e-12_r8
      rate(:,:,329) = 6.80e-13_r8
      rate(:,:,330) = 5.40e-11_r8
      rate(:,:,342) = 2.40e-12_r8
      rate(:,:,345) = 1.40e-11_r8
      rate(:,:,348) = 5.00e-12_r8
      rate(:,:,360) = 2.40e-12_r8
      rate(:,:,364) = 1.40e-11_r8
      rate(:,:,366) = 2.40e-12_r8
      rate(:,:,368) = 3.50e-12_r8
      rate(:,:,369) = 4.50e-11_r8
      rate(:,:,376) = 2.40e-12_r8
      rate(:,:,386) = 3.00e-12_r8
      rate(:,:,387) = 1.00e-11_r8
      rate(:,:,391) = 2.3e-11_r8
      rate(:,:,403) = 7.10e-6_r8
      rate(:,:,409) = 7.10e-6_r8
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
      rate(:,:,429) = 6.34e-8_r8
      rate(:,:,430) = 6.34e-8_r8
      rate(:,:,431) = 6.34e-8_r8
      rate(:,:,453) = 9.0e-10_r8
      rate(:,:,454) = 1.0e-10_r8
      rate(:,:,455) = 4.4e-10_r8
      rate(:,:,456) = 4.0e-10_r8
      rate(:,:,457) = 2.0e-10_r8
      rate(:,:,458) = 1.0e-12_r8
      rate(:,:,459) = 6.0e-11_r8
      rate(:,:,460) = 5.0e-16_r8
      rate(:,:,464) = 2.31e-06_r8
      rate(:,:,465) = 2.31e-07_r8
      rate(:,:,466) = 2.31e-07_r8
      rate(:,:,467) = 4.63e-07_r8
      rate(:,:,468) = 4.63e-07_r8
      rate(:,:,469) = 2.31e-07_r8
      rate(:,:,470) = 1.29e-07_r8
      rate(:,:,471) = 1.29e-07_r8
      rate(:,:,472) = 1.29e-07_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,116) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,120) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,121) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,125) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,128) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,129) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:,130) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:,131) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,132) = 7.25e-11_r8 * exp_fac(:,:)
      rate(:,:,133) = 4.63e-11_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 100._r8 * itemp(:,:) )
      rate(:,:,159) = 7.70e-11_r8 * exp_fac(:,:)
      rate(:,:,183) = 2.10e-11_r8 * exp_fac(:,:)
      rate(:,:,161) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 180._r8 * itemp(:,:) )
      rate(:,:,165) = 1.80e-11_r8 * exp_fac(:,:)
      rate(:,:,286) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,313) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,318) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,331) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,335) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,359) = 4.40e-12_r8 * exp_fac(:,:)
      rate(:,:,372) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,383) = 4.20e-12_r8 * exp_fac(:,:)
      rate(:,:,397) = 4.2e-12_r8 * exp_fac(:,:)
      rate(:,:,166) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,167) = 4.80e-11_r8 * exp_fac(:,:)
      rate(:,:,235) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,170) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      rate(:,:,172) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,173) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,243) = 5.50e-12_r8 * exp_fac(:,:)
      rate(:,:,275) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,296) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,316) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,320) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,325) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,337) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,346) = 2.30e-11_r8 * exp_fac(:,:)
      rate(:,:,362) = 1.52e-11_r8 * exp_fac(:,:)
      rate(:,:,374) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,385) = 3.80e-12_r8 * exp_fac(:,:)
      rate(:,:,399) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,174) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2000._r8 * itemp(:,:) )
      rate(:,:,176) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,357) = 1.05e-14_r8 * exp_fac(:,:)
      rate(:,:,178) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:,:) )
      rate(:,:,182) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,184) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,185) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,186) = 1.45e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,188) = 3.30e-12_r8 * exp_fac(:,:)
      rate(:,:,207) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,212) = 7.40e-12_r8 * exp_fac(:,:)
      rate(:,:,299) = 8.10e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,189) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,244) = 5.80e-12_r8 * exp_fac(:,:)
      rate(:,:,190) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,192) = 1.20e-13_r8 * exp_fac(:,:)
      rate(:,:,218) = 3.00e-11_r8 * exp_fac(:,:)
      rate(:,:,197) = 1.50e-11_r8 * exp( 170._r8 * itemp(:,:) )
      rate(:,:,202) = 1.30e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,204) = 2.30e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,205) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,206) = 1.10e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,208) = 3.60e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,209) = 8.10e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,210) = 7.30e-12_r8 * exp( -1280._r8 * itemp(:,:) )
      rate(:,:,211) = 2.80e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,213) = 6.00e-13_r8 * exp_fac(:,:)
      rate(:,:,234) = 1.90e-11_r8 * exp_fac(:,:)
      rate(:,:,242) = 1.50e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,214) = 2.60e-12_r8 * exp_fac(:,:)
      rate(:,:,216) = 6.40e-12_r8 * exp_fac(:,:)
      rate(:,:,241) = 4.10e-13_r8 * exp_fac(:,:)
      rate(:,:,215) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,219) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,220) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      rate(:,:,223) = 1.80e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,224) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,226) = 3.40e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,227) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,248) = 1.40e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,228) = 3.60e-12_r8 * exp_fac(:,:)
      rate(:,:,259) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,229) = 1.20e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,230) = 6.50e-12_r8 * exp( 135._r8 * itemp(:,:) )
      rate(:,:,231) = 1.60e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,232) = 4.80e-12_r8 * exp( -310._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,233) = 1.70e-11_r8 * exp_fac(:,:)
      rate(:,:,261) = 6.30e-12_r8 * exp_fac(:,:)
      rate(:,:,236) = 4.50e-12_r8 * exp( 460._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,237) = 8.80e-12_r8 * exp_fac(:,:)
      rate(:,:,240) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,239) = 9.50e-13_r8 * exp( 550._r8 * itemp(:,:) )
      rate(:,:,245) = 1.20e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,246) = 1.90e-11_r8 * exp( 215._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,247) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,276) = 4.50e-13_r8 * exp_fac(:,:)
      rate(:,:,288) = 1.00e-14_r8 * exp_fac(:,:)
      rate(:,:,291) = 7.00e-13_r8 * exp_fac(:,:)
      rate(:,:,294) = 2.00e-13_r8 * exp_fac(:,:)
      rate(:,:,295) = 6.80e-14_r8 * exp_fac(:,:)
      rate(:,:,304) = 1.00e-12_r8 * exp_fac(:,:)
      rate(:,:,305) = 1.00e-11_r8 * exp_fac(:,:)
      rate(:,:,306) = 1.15e-11_r8 * exp_fac(:,:)
      rate(:,:,309) = 4.00e-14_r8 * exp_fac(:,:)
      rate(:,:,326) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,329) = 6.80e-13_r8 * exp_fac(:,:)
      rate(:,:,330) = 5.40e-11_r8 * exp_fac(:,:)
      rate(:,:,342) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,345) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,348) = 5.00e-12_r8 * exp_fac(:,:)
      rate(:,:,360) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,364) = 1.40e-11_r8 * exp_fac(:,:)
      rate(:,:,366) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,368) = 3.50e-12_r8 * exp_fac(:,:)
      rate(:,:,369) = 4.50e-11_r8 * exp_fac(:,:)
      rate(:,:,376) = 2.40e-12_r8 * exp_fac(:,:)
      rate(:,:,386) = 3.00e-12_r8 * exp_fac(:,:)
      rate(:,:,387) = 1.00e-11_r8 * exp_fac(:,:)
      rate(:,:,391) = 2.3e-11_r8 * exp_fac(:,:)
      rate(:,:,403) = 7.10e-6_r8 * exp_fac(:,:)
      rate(:,:,409) = 7.10e-6_r8 * exp_fac(:,:)
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
      rate(:,:,429) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,430) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,431) = 6.34e-8_r8 * exp_fac(:,:)
      rate(:,:,453) = 9.0e-10_r8 * exp_fac(:,:)
      rate(:,:,454) = 1.0e-10_r8 * exp_fac(:,:)
      rate(:,:,455) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,456) = 4.0e-10_r8 * exp_fac(:,:)
      rate(:,:,457) = 2.0e-10_r8 * exp_fac(:,:)
      rate(:,:,458) = 1.0e-12_r8 * exp_fac(:,:)
      rate(:,:,459) = 6.0e-11_r8 * exp_fac(:,:)
      rate(:,:,460) = 5.0e-16_r8 * exp_fac(:,:)
      rate(:,:,464) = 2.31e-06_r8 * exp_fac(:,:)
      rate(:,:,465) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,466) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,467) = 4.63e-07_r8 * exp_fac(:,:)
      rate(:,:,468) = 4.63e-07_r8 * exp_fac(:,:)
      rate(:,:,469) = 2.31e-07_r8 * exp_fac(:,:)
      rate(:,:,470) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,471) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,472) = 1.29e-07_r8 * exp_fac(:,:)
      rate(:,:,249) = 1.60e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 400._r8 * itemp(:,:) )
      rate(:,:,250) = 6.00e-12_r8 * exp_fac(:,:)
      rate(:,:,344) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,363) = 5.00e-13_r8 * exp_fac(:,:)
      rate(:,:,378) = 5.e-13_r8 * exp_fac(:,:)
      rate(:,:,251) = 2.17e-11_r8 * exp( -1130._r8 * itemp(:,:) )
      rate(:,:,252) = 2.40e-12_r8 * exp( -1250._r8 * itemp(:,:) )
      rate(:,:,253) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,254) = 1.05e-12_r8 * exp_fac(:,:)
      rate(:,:,257) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,268) = 3.40e-11_r8 * exp_fac(:,:)
      rate(:,:,255) = 2.35e-12_r8 * exp( -1300._r8 * itemp(:,:) )
      rate(:,:,256) = 1.40e-11_r8 * exp( -1030._r8 * itemp(:,:) )
      rate(:,:,258) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,260) = 1.35e-12_r8 * exp( -600._r8 * itemp(:,:) )
      rate(:,:,262) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,263) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )
      rate(:,:,266) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,267) = 5.50e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,269) = 9.7e-15_r8 * exp( 625._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 300._r8 * itemp(:,:) )
      rate(:,:,270) = 2.80e-12_r8 * exp_fac(:,:)
      rate(:,:,322) = 2.90e-12_r8 * exp_fac(:,:)
      rate(:,:,271) = 4.10e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,272) = 5.00e-13_r8 * exp( -424._r8 * itemp(:,:) )
      rate(:,:,273) = 1.90e-14_r8 * exp( 706._r8 * itemp(:,:) )
      rate(:,:,274) = 2.90e-12_r8 * exp( -345._r8 * itemp(:,:) )
      rate(:,:,277) = 2.40e12_r8 * exp( -7000._r8 * itemp(:,:) )
      rate(:,:,278) = 2.60e-12_r8 * exp( 265._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 700._r8 * itemp(:,:) )
      rate(:,:,279) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,287) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,293) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,314) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,319) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,323) = 8.60e-13_r8 * exp_fac(:,:)
      rate(:,:,336) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,343) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,361) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,367) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,373) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,377) = 8.00e-13_r8 * exp_fac(:,:)
      rate(:,:,384) = 7.50e-13_r8 * exp_fac(:,:)
      rate(:,:,389) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,392) = 1.4e-12_r8 * exp_fac(:,:)
      rate(:,:,398) = 7.5e-13_r8 * exp_fac(:,:)
      rate(:,:,282) = 7.20e-11_r8 * exp( -70._r8 * itemp(:,:) )
      rate(:,:,284) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:,:) )
      rate(:,:,289) = 1.60e11_r8 * exp( -4150._r8 * itemp(:,:) )
      rate(:,:,290) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:,:) )
      rate(:,:,292) = 2.60e-12_r8 * exp( 365._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 350._r8 * itemp(:,:) )
      rate(:,:,297) = 4.63e-12_r8 * exp_fac(:,:)
      rate(:,:,390) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,393) = 2.6e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1900._r8 * itemp(:,:) )
      rate(:,:,298) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,311) = 6.50e-15_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._r8 * itemp(:,:) )
      rate(:,:,301) = 4.30e-13_r8 * exp_fac(:,:)
      rate(:,:,349) = 4.30e-13_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 500._r8 * itemp(:,:) )
      rate(:,:,302) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,303) = 2.50e-12_r8 * exp_fac(:,:)
      rate(:,:,324) = 7.10e-13_r8 * exp_fac(:,:)
      rate(:,:,350) = 2.00e-12_r8 * exp_fac(:,:)
      rate(:,:,307) = 6.90e-12_r8 * exp( -230._r8 * itemp(:,:) )
      rate(:,:,312) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:,:) )
      rate(:,:,315) = 3.75e-13_r8 * exp( -40._r8 * itemp(:,:) )
      rate(:,:,317) = 8.70e-12_r8 * exp( -615._r8 * itemp(:,:) )
      rate(:,:,327) = 8.40e-13_r8 * exp( 830._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._r8 * itemp(:,:) )
      rate(:,:,328) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,370) = 1.40e-12_r8 * exp_fac(:,:)
      rate(:,:,332) = 4.13e-12_r8 * exp( 452._r8 * itemp(:,:) )
      rate(:,:,333) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:,:) )
      rate(:,:,334) = 2.30e-12_r8 * exp( -170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 175._r8 * itemp(:,:) )
      rate(:,:,338) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,371) = 1.86e-11_r8 * exp_fac(:,:)
      rate(:,:,339) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 360._r8 * itemp(:,:) )
      rate(:,:,340) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,341) = 1.30e-13_r8 * exp_fac(:,:)
      rate(:,:,347) = 5.30e-12_r8 * exp_fac(:,:)
      rate(:,:,365) = 2.70e-12_r8 * exp_fac(:,:)
      rate(:,:,375) = 2.7e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 530._r8 * itemp(:,:) )
      rate(:,:,351) = 4.60e-12_r8 * exp_fac(:,:)
      rate(:,:,352) = 2.30e-12_r8 * exp_fac(:,:)
      rate(:,:,356) = 2.54e-11_r8 * exp( 410._r8 * itemp(:,:) )
      rate(:,:,358) = 3.03e-12_r8 * exp( -446._r8 * itemp(:,:) )
      rate(:,:,379) = 1.3e-12_r8 * exp( 640._r8 * itemp(:,:) )
      rate(:,:,380) = 1.90e-12_r8 * exp( 190._r8 * itemp(:,:) )
      rate(:,:,382) = 1.70e-12_r8 * exp( 352._r8 * itemp(:,:) )
      rate(:,:,388) = 2.3e-12_r8 * exp( -193._r8 * itemp(:,:) )
      rate(:,:,394) = 1.2e-11_r8 * exp( 444._r8 * itemp(:,:) )
      rate(:,:,395) = 1.e-15_r8 * exp( -732._r8 * itemp(:,:) )
      rate(:,:,396) = 1.2e-12_r8 * exp( 490._r8 * itemp(:,:) )
      rate(:,:,405) = 9.60e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,407) = 1.90e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,408) = 1.70e-12_r8 * exp( -710._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,160), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.90e-31_r8 * itemp(:,:)**1.0_r8
      kinf(:,:) = 2.60e-11_r8
      call jpl( rate(1,1,169), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.28e-33_r8
      kinf(:,:) = 9.30e-15_r8 * itemp(:,:)**(-4.42_r8)
      call jpl( rate(1,1,177), m, 0.8_r8, ko, kinf, n )

      ko(:,:) = 9.00e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3.0e-11_r8
      call jpl( rate(1,1,187), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.50e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,191), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-30_r8 * itemp(:,:)**4.4_r8
      kinf(:,:) = 1.4e-12_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,193), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-30_r8 * itemp(:,:)**3.0_r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,195), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.00e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 2.9e-12_r8 * itemp(:,:)**1.1_r8
      call jpl( rate(1,1,201), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.80e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,217), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-32_r8 * itemp(:,:)**4.5_r8
      kinf(:,:) = 3.0e-12_r8 * itemp(:,:)**2.0_r8
      call jpl( rate(1,1,221), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,238), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.90e-33_r8 * itemp(:,:)**1.4_r8
      kinf(:,:) = 1.10e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,265), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.20e-30_r8 * itemp(:,:)**2.4_r8
      kinf(:,:) = 2.2e-10_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,280), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.60e-29_r8 * itemp(:,:)**3.3_r8
      kinf(:,:) = 3.1e-10_r8 * itemp(:,:)
      call jpl( rate(1,1,281), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.50e-30_r8
      kinf(:,:) = 8.3e-13_r8 * itemp(:,:)**(-2.0_r8)
      call jpl( rate(1,1,283), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.60e-29_r8 * itemp(:,:)**3.1_r8
      kinf(:,:) = 9.00e-12_r8 * itemp(:,:)**0.85_r8
      call jpl( rate(1,1,285), m, 0.48_r8, ko, kinf, n )

      ko(:,:) = 9.70e-29_r8 * itemp(:,:)**5.6_r8
      kinf(:,:) = 9.30e-12_r8 * itemp(:,:)**1.5_r8
      call jpl( rate(1,1,300), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,310), m, 0.5_r8, ko, kinf, n )

      ko(:,:) = 8.00e-27_r8 * itemp(:,:)**3.5_r8
      kinf(:,:) = 3.00e-11_r8
      call jpl( rate(1,1,355), m, 0.5_r8, ko, kinf, n )

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

      rate(:,:kbot,118) = 8.00e-14_r8
      rate(:,:kbot,119) = 3.90e-17_r8
      rate(:,:kbot,124) = 1.30e-16_r8
      rate(:,:kbot,126) = 1.00e-20_r8
      rate(:,:kbot,163) = 6.90e-12_r8
      rate(:,:kbot,179) = 5.00e-12_r8
      rate(:,:kbot,180) = 7.00e-13_r8
      rate(:,:kbot,454) = 1.0e-10_r8
      rate(:,:kbot,455) = 4.4e-10_r8
      rate(:,:kbot,456) = 4.0e-10_r8
      rate(:,:kbot,457) = 2.0e-10_r8
      rate(:,:kbot,458) = 1.0e-12_r8
      rate(:,:kbot,459) = 6.0e-11_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,116) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,120) = 1.80e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,121) = 3.50e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,125) = 3.60e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,128) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,129) = 3.135e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,130) = 1.65e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,161) = 1.40e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,165) = 1.80e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,166) = 1.70e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,167) = 4.80e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,173) = 3.00e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,174) = 1.00e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,182) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,183) = 2.10e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,188) = 3.30e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,189) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:,:) )
      rate(:,:kbot,190) = 5.10e-12_r8 * exp( 210._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 4.40e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,160) = wrk(:,:)



















      end subroutine setrxt_hrates

      end module mo_setrxt
