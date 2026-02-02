
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      real(r8)  :: itemp(ncol*pver)
      real(r8)  :: exp_fac(ncol*pver)
      real(r8)  :: ko(ncol*pver)
      real(r8)  :: kinf(ncol*pver)

      rate(:,124) = 1.29e-07_r8
      rate(:,125) = 1.2e-10_r8
      rate(:,129) = 1.2e-10_r8
      rate(:,130) = 1.2e-10_r8
      rate(:,136) = 6.9e-12_r8
      rate(:,137) = 7.2e-11_r8
      rate(:,138) = 1.6e-12_r8
      rate(:,144) = 1.8e-12_r8
      rate(:,148) = 1.8e-12_r8
      rate(:,160) = 3.5e-12_r8
      rate(:,162) = 1.3e-11_r8
      rate(:,163) = 2.2e-11_r8
      rate(:,164) = 5e-11_r8
      rate(:,199) = 1.7e-13_r8
      rate(:,201) = 2.607e-10_r8
      rate(:,202) = 9.75e-11_r8
      rate(:,203) = 2.07e-10_r8
      rate(:,204) = 2.088e-10_r8
      rate(:,205) = 1.17e-10_r8
      rate(:,206) = 4.644e-11_r8
      rate(:,207) = 1.204e-10_r8
      rate(:,208) = 9.9e-11_r8
      rate(:,209) = 3.3e-12_r8
      rate(:,228) = 4.5e-11_r8
      rate(:,229) = 4.62e-10_r8
      rate(:,230) = 1.2e-10_r8
      rate(:,231) = 9e-11_r8
      rate(:,232) = 3e-11_r8
      rate(:,237) = 2.14e-11_r8
      rate(:,238) = 1.9e-10_r8
      rate(:,251) = 2.57e-10_r8
      rate(:,252) = 1.8e-10_r8
      rate(:,253) = 1.794e-10_r8
      rate(:,254) = 1.3e-10_r8
      rate(:,255) = 7.65e-11_r8
      rate(:,268) = 4e-13_r8
      rate(:,272) = 1.31e-10_r8
      rate(:,273) = 3.5e-11_r8
      rate(:,274) = 9e-12_r8
      rate(:,281) = 6.8e-14_r8
      rate(:,282) = 2e-13_r8
      rate(:,297) = 1e-12_r8
      rate(:,301) = 1e-14_r8
      rate(:,302) = 1e-11_r8
      rate(:,303) = 1.15e-11_r8
      rate(:,304) = 4e-14_r8
      rate(:,317) = 3e-12_r8
      rate(:,318) = 6.7e-13_r8
      rate(:,328) = 3.5e-13_r8
      rate(:,329) = 5.4e-11_r8
      rate(:,332) = 2e-12_r8
      rate(:,333) = 1.4e-11_r8
      rate(:,336) = 2.4e-12_r8
      rate(:,347) = 5e-12_r8
      rate(:,357) = 2.2e-12_r8
      rate(:,359) = 6.7e-12_r8
      rate(:,362) = 3.5e-12_r8
      rate(:,365) = 1.3e-11_r8
      rate(:,366) = 1.4e-11_r8
      rate(:,370) = 2.4e-12_r8
      rate(:,371) = 1.4e-11_r8
      rate(:,376) = 2.4e-12_r8
      rate(:,377) = 4e-11_r8
      rate(:,378) = 4e-11_r8
      rate(:,380) = 1.4e-11_r8
      rate(:,384) = 2.4e-12_r8
      rate(:,385) = 4e-11_r8
      rate(:,389) = 7e-11_r8
      rate(:,390) = 1e-10_r8
      rate(:,395) = 2.4e-12_r8
      rate(:,410) = 4.7e-11_r8
      rate(:,423) = 2.1e-12_r8
      rate(:,424) = 2.8e-13_r8
      rate(:,432) = 1.7e-11_r8
      rate(:,438) = 8.4e-11_r8
      rate(:,440) = 1.9e-11_r8
      rate(:,441) = 1.2e-14_r8
      rate(:,442) = 2e-10_r8
      rate(:,449) = 2.4e-12_r8
      rate(:,450) = 2e-11_r8
      rate(:,454) = 2.3e-11_r8
      rate(:,455) = 2e-11_r8
      rate(:,459) = 3.3e-11_r8
      rate(:,460) = 1e-12_r8
      rate(:,461) = 5.7e-11_r8
      rate(:,462) = 3.4e-11_r8
      rate(:,467) = 2.3e-12_r8
      rate(:,469) = 1.2e-11_r8
      rate(:,470) = 5.7e-11_r8
      rate(:,471) = 2.8e-11_r8
      rate(:,472) = 6.6e-11_r8
      rate(:,473) = 1.4e-11_r8
      rate(:,476) = 1.9e-12_r8
      rate(:,488) = 6.34e-08_r8
      rate(:,494) = 1.9e-11_r8
      rate(:,497) = 1.2e-14_r8
      rate(:,498) = 2e-10_r8
      rate(:,509) = 1.34e-11_r8
      rate(:,515) = 1.34e-11_r8
      rate(:,520) = 1.7e-11_r8
      rate(:,540) = 2.31e-07_r8
      rate(:,541) = 2.31e-06_r8
      rate(:,542) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,126) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,128) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,131) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,134) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,135) = 1.4e-12_r8 * exp_fac(:)
      rate(:,386) = 1.05e-14_r8 * exp_fac(:)
      rate(:,505) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,140) = 3e-11_r8 * exp_fac(:)
      rate(:,226) = 5.5e-12_r8 * exp_fac(:)
      rate(:,265) = 3.8e-12_r8 * exp_fac(:)
      rate(:,286) = 3.8e-12_r8 * exp_fac(:)
      rate(:,313) = 3.8e-12_r8 * exp_fac(:)
      rate(:,321) = 3.8e-12_r8 * exp_fac(:)
      rate(:,325) = 3.8e-12_r8 * exp_fac(:)
      rate(:,341) = 2.3e-11_r8 * exp_fac(:)
      rate(:,351) = 3.8e-12_r8 * exp_fac(:)
      rate(:,361) = 3.8e-12_r8 * exp_fac(:)
      rate(:,388) = 1.52e-11_r8 * exp_fac(:)
      rate(:,396) = 1.52e-12_r8 * exp_fac(:)
      rate(:,402) = 3.8e-12_r8 * exp_fac(:)
      rate(:,405) = 3.8e-12_r8 * exp_fac(:)
      rate(:,409) = 3.8e-12_r8 * exp_fac(:)
      rate(:,425) = 3.8e-12_r8 * exp_fac(:)
      rate(:,429) = 3.8e-12_r8 * exp_fac(:)
      rate(:,435) = 3.8e-12_r8 * exp_fac(:)
      rate(:,439) = 3.8e-12_r8 * exp_fac(:)
      rate(:,141) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,142) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,143) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,145) = 4.8e-11_r8 * exp_fac(:)
      rate(:,224) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,146) = 1.8e-11_r8 * exp_fac(:)
      rate(:,299) = 4.2e-12_r8 * exp_fac(:)
      rate(:,320) = 4.2e-12_r8 * exp_fac(:)
      rate(:,349) = 4.2e-12_r8 * exp_fac(:)
      rate(:,369) = 4.4e-12_r8 * exp_fac(:)
      rate(:,375) = 4.4e-12_r8 * exp_fac(:)
      rate(:,448) = 4.2e-12_r8 * exp_fac(:)
      rate(:,453) = 4.2e-12_r8 * exp_fac(:)
      rate(:,458) = 4.2e-12_r8 * exp_fac(:)
      rate(:,147) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,151) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,152) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,153) = 2.9e-12_r8 * exp_fac(:)
      rate(:,154) = 1.45e-12_r8 * exp_fac(:)
      rate(:,155) = 1.45e-12_r8 * exp_fac(:)
      rate(:,156) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,157) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,158) = 1.2e-13_r8 * exp_fac(:)
      rate(:,184) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,161) = 1.7e-11_r8 * exp_fac(:)
      rate(:,259) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,165) = 3.44e-12_r8 * exp_fac(:)
      rate(:,217) = 2.3e-12_r8 * exp_fac(:)
      rate(:,220) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,166) = 3e-12_r8 * exp_fac(:)
      rate(:,225) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,168) = 7.26e-11_r8 * exp_fac(:)
      rate(:,169) = 4.64e-11_r8 * exp_fac(:)
      rate(:,176) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,177) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,178) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,179) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,180) = 1.4e-11_r8 * exp_fac(:)
      rate(:,194) = 7.4e-12_r8 * exp_fac(:)
      rate(:,295) = 8.1e-12_r8 * exp_fac(:)
      rate(:,181) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,182) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,183) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,185) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,186) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,187) = 2.6e-12_r8 * exp_fac(:)
      rate(:,188) = 6.4e-12_r8 * exp_fac(:)
      rate(:,218) = 4.1e-13_r8 * exp_fac(:)
      rate(:,398) = 7.5e-12_r8 * exp_fac(:)
      rate(:,412) = 7.5e-12_r8 * exp_fac(:)
      rate(:,415) = 7.5e-12_r8 * exp_fac(:)
      rate(:,418) = 7.5e-12_r8 * exp_fac(:)
      rate(:,189) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,191) = 3.6e-12_r8 * exp_fac(:)
      rate(:,240) = 2e-12_r8 * exp_fac(:)
      rate(:,192) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,193) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,195) = 6e-13_r8 * exp_fac(:)
      rate(:,215) = 1.5e-12_r8 * exp_fac(:)
      rate(:,223) = 1.9e-11_r8 * exp_fac(:)
      rate(:,196) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,197) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,198) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,200) = 3e-12_r8 * exp_fac(:)
      rate(:,234) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,212) = 1.7e-11_r8 * exp_fac(:)
      rate(:,239) = 6.3e-12_r8 * exp_fac(:)
      rate(:,213) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,214) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,216) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,219) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,222) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,227) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,233) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,235) = 1.4e-11_r8 * exp_fac(:)
      rate(:,237) = 2.14e-11_r8 * exp_fac(:)
      rate(:,238) = 1.9e-10_r8 * exp_fac(:)
      rate(:,251) = 2.57e-10_r8 * exp_fac(:)
      rate(:,252) = 1.8e-10_r8 * exp_fac(:)
      rate(:,253) = 1.794e-10_r8 * exp_fac(:)
      rate(:,254) = 1.3e-10_r8 * exp_fac(:)
      rate(:,255) = 7.65e-11_r8 * exp_fac(:)
      rate(:,268) = 4e-13_r8 * exp_fac(:)
      rate(:,272) = 1.31e-10_r8 * exp_fac(:)
      rate(:,273) = 3.5e-11_r8 * exp_fac(:)
      rate(:,274) = 9e-12_r8 * exp_fac(:)
      rate(:,281) = 6.8e-14_r8 * exp_fac(:)
      rate(:,282) = 2e-13_r8 * exp_fac(:)
      rate(:,297) = 1e-12_r8 * exp_fac(:)
      rate(:,301) = 1e-14_r8 * exp_fac(:)
      rate(:,302) = 1e-11_r8 * exp_fac(:)
      rate(:,303) = 1.15e-11_r8 * exp_fac(:)
      rate(:,304) = 4e-14_r8 * exp_fac(:)
      rate(:,317) = 3e-12_r8 * exp_fac(:)
      rate(:,318) = 6.7e-13_r8 * exp_fac(:)
      rate(:,328) = 3.5e-13_r8 * exp_fac(:)
      rate(:,329) = 5.4e-11_r8 * exp_fac(:)
      rate(:,332) = 2e-12_r8 * exp_fac(:)
      rate(:,333) = 1.4e-11_r8 * exp_fac(:)
      rate(:,336) = 2.4e-12_r8 * exp_fac(:)
      rate(:,347) = 5e-12_r8 * exp_fac(:)
      rate(:,357) = 2.2e-12_r8 * exp_fac(:)
      rate(:,359) = 6.7e-12_r8 * exp_fac(:)
      rate(:,362) = 3.5e-12_r8 * exp_fac(:)
      rate(:,365) = 1.3e-11_r8 * exp_fac(:)
      rate(:,366) = 1.4e-11_r8 * exp_fac(:)
      rate(:,370) = 2.4e-12_r8 * exp_fac(:)
      rate(:,371) = 1.4e-11_r8 * exp_fac(:)
      rate(:,376) = 2.4e-12_r8 * exp_fac(:)
      rate(:,377) = 4e-11_r8 * exp_fac(:)
      rate(:,378) = 4e-11_r8 * exp_fac(:)
      rate(:,380) = 1.4e-11_r8 * exp_fac(:)
      rate(:,384) = 2.4e-12_r8 * exp_fac(:)
      rate(:,385) = 4e-11_r8 * exp_fac(:)
      rate(:,389) = 7e-11_r8 * exp_fac(:)
      rate(:,390) = 1e-10_r8 * exp_fac(:)
      rate(:,395) = 2.4e-12_r8 * exp_fac(:)
      rate(:,410) = 4.7e-11_r8 * exp_fac(:)
      rate(:,423) = 2.1e-12_r8 * exp_fac(:)
      rate(:,424) = 2.8e-13_r8 * exp_fac(:)
      rate(:,432) = 1.7e-11_r8 * exp_fac(:)
      rate(:,438) = 8.4e-11_r8 * exp_fac(:)
      rate(:,440) = 1.9e-11_r8 * exp_fac(:)
      rate(:,441) = 1.2e-14_r8 * exp_fac(:)
      rate(:,442) = 2e-10_r8 * exp_fac(:)
      rate(:,449) = 2.4e-12_r8 * exp_fac(:)
      rate(:,450) = 2e-11_r8 * exp_fac(:)
      rate(:,454) = 2.3e-11_r8 * exp_fac(:)
      rate(:,455) = 2e-11_r8 * exp_fac(:)
      rate(:,459) = 3.3e-11_r8 * exp_fac(:)
      rate(:,460) = 1e-12_r8 * exp_fac(:)
      rate(:,461) = 5.7e-11_r8 * exp_fac(:)
      rate(:,462) = 3.4e-11_r8 * exp_fac(:)
      rate(:,467) = 2.3e-12_r8 * exp_fac(:)
      rate(:,469) = 1.2e-11_r8 * exp_fac(:)
      rate(:,470) = 5.7e-11_r8 * exp_fac(:)
      rate(:,471) = 2.8e-11_r8 * exp_fac(:)
      rate(:,472) = 6.6e-11_r8 * exp_fac(:)
      rate(:,473) = 1.4e-11_r8 * exp_fac(:)
      rate(:,476) = 1.9e-12_r8 * exp_fac(:)
      rate(:,488) = 6.34e-08_r8 * exp_fac(:)
      rate(:,494) = 1.9e-11_r8 * exp_fac(:)
      rate(:,497) = 1.2e-14_r8 * exp_fac(:)
      rate(:,498) = 2e-10_r8 * exp_fac(:)
      rate(:,509) = 1.34e-11_r8 * exp_fac(:)
      rate(:,515) = 1.34e-11_r8 * exp_fac(:)
      rate(:,520) = 1.7e-11_r8 * exp_fac(:)
      rate(:,540) = 2.31e-07_r8 * exp_fac(:)
      rate(:,541) = 2.31e-06_r8 * exp_fac(:)
      rate(:,542) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,236) = 6e-12_r8 * exp_fac(:)
      rate(:,334) = 5e-13_r8 * exp_fac(:)
      rate(:,367) = 5e-13_r8 * exp_fac(:)
      rate(:,372) = 5e-13_r8 * exp_fac(:)
      rate(:,381) = 5e-13_r8 * exp_fac(:)
      rate(:,392) = 5e-13_r8 * exp_fac(:)
      rate(:,241) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,242) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,243) = 1.64e-12_r8 * exp_fac(:)
      rate(:,353) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,244) = 2.03e-11_r8 * exp_fac(:)
      rate(:,475) = 3.4e-12_r8 * exp_fac(:)
      rate(:,245) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,246) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,247) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,248) = 1.25e-12_r8 * exp_fac(:)
      rate(:,258) = 3.4e-11_r8 * exp_fac(:)
      rate(:,249) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,250) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,256) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,257) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,260) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,261) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,262) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,263) = 2.8e-12_r8 * exp_fac(:)
      rate(:,324) = 2.9e-12_r8 * exp_fac(:)
      rate(:,264) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,266) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,269) = 7.5e-13_r8 * exp_fac(:)
      rate(:,283) = 7.5e-13_r8 * exp_fac(:)
      rate(:,298) = 7.5e-13_r8 * exp_fac(:)
      rate(:,319) = 7.5e-13_r8 * exp_fac(:)
      rate(:,323) = 8.6e-13_r8 * exp_fac(:)
      rate(:,335) = 8e-13_r8 * exp_fac(:)
      rate(:,348) = 7.5e-13_r8 * exp_fac(:)
      rate(:,358) = 7.5e-13_r8 * exp_fac(:)
      rate(:,368) = 8e-13_r8 * exp_fac(:)
      rate(:,373) = 8e-13_r8 * exp_fac(:)
      rate(:,382) = 8e-13_r8 * exp_fac(:)
      rate(:,393) = 8e-13_r8 * exp_fac(:)
      rate(:,400) = 7.5e-13_r8 * exp_fac(:)
      rate(:,404) = 7.5e-13_r8 * exp_fac(:)
      rate(:,407) = 7.5e-13_r8 * exp_fac(:)
      rate(:,420) = 7.5e-13_r8 * exp_fac(:)
      rate(:,427) = 7.5e-13_r8 * exp_fac(:)
      rate(:,433) = 7.5e-13_r8 * exp_fac(:)
      rate(:,436) = 7.5e-13_r8 * exp_fac(:)
      rate(:,447) = 7.5e-13_r8 * exp_fac(:)
      rate(:,452) = 7.5e-13_r8 * exp_fac(:)
      rate(:,457) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 7.5e-13_r8 * exp_fac(:)
      rate(:,507) = 7.5e-13_r8 * exp_fac(:)
      rate(:,517) = 7.5e-13_r8 * exp_fac(:)
      rate(:,521) = 7.5e-13_r8 * exp_fac(:)
      rate(:,270) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,271) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,275) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,280) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,284) = 2.6e-12_r8 * exp_fac(:)
      rate(:,401) = 2.6e-12_r8 * exp_fac(:)
      rate(:,406) = 2.6e-12_r8 * exp_fac(:)
      rate(:,408) = 2.6e-12_r8 * exp_fac(:)
      rate(:,421) = 2.6e-12_r8 * exp_fac(:)
      rate(:,428) = 2.6e-12_r8 * exp_fac(:)
      rate(:,434) = 2.6e-12_r8 * exp_fac(:)
      rate(:,437) = 2.6e-12_r8 * exp_fac(:)
      rate(:,501) = 2.6e-12_r8 * exp_fac(:)
      rate(:,508) = 2.6e-12_r8 * exp_fac(:)
      rate(:,518) = 2.6e-12_r8 * exp_fac(:)
      rate(:,522) = 2.6e-12_r8 * exp_fac(:)
      rate(:,285) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,287) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,288) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,289) = 1.4e-12_r8 * exp_fac(:)
      rate(:,309) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,290) = 4.63e-12_r8 * exp_fac(:)
      rate(:,504) = 2.7e-12_r8 * exp_fac(:)
      rate(:,291) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,292) = 2.9e-12_r8 * exp_fac(:)
      rate(:,293) = 2e-12_r8 * exp_fac(:)
      rate(:,322) = 7.1e-13_r8 * exp_fac(:)
      rate(:,343) = 2e-12_r8 * exp_fac(:)
      rate(:,446) = 2e-12_r8 * exp_fac(:)
      rate(:,451) = 2e-12_r8 * exp_fac(:)
      rate(:,456) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,294) = 4.3e-13_r8 * exp_fac(:)
      rate(:,344) = 4.3e-13_r8 * exp_fac(:)
      rate(:,397) = 4.3e-13_r8 * exp_fac(:)
      rate(:,411) = 4.3e-13_r8 * exp_fac(:)
      rate(:,414) = 4.3e-13_r8 * exp_fac(:)
      rate(:,417) = 4.3e-13_r8 * exp_fac(:)
      rate(:,296) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,300) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,308) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,310) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,311) = 1.41e-13_r8 * exp_fac(:)
      rate(:,495) = 2.75e-13_r8 * exp_fac(:)
      rate(:,503) = 2.12e-13_r8 * exp_fac(:)
      rate(:,511) = 2.6e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,312) = 2.7e-12_r8 * exp_fac(:)
      rate(:,337) = 2.7e-12_r8 * exp_fac(:)
      rate(:,338) = 1.3e-13_r8 * exp_fac(:)
      rate(:,340) = 9.6e-12_r8 * exp_fac(:)
      rate(:,346) = 5.3e-12_r8 * exp_fac(:)
      rate(:,383) = 2.7e-12_r8 * exp_fac(:)
      rate(:,394) = 2.7e-12_r8 * exp_fac(:)
      rate(:,496) = 2.7e-12_r8 * exp_fac(:)
      rate(:,512) = 2.7e-12_r8 * exp_fac(:)
      rate(:,314) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,315) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,316) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,330) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,331) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,339) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,342) = 4.6e-12_r8 * exp_fac(:)
      rate(:,345) = 2.3e-12_r8 * exp_fac(:)
      rate(:,350) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,354) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,360) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,363) = 1.86e-11_r8 * exp_fac(:)
      rate(:,364) = 1.86e-11_r8 * exp_fac(:)
      rate(:,374) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,379) = 3.03e-12_r8 * exp_fac(:)
      rate(:,502) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,387) = 2.54e-11_r8 * exp_fac(:)
      rate(:,506) = 2.54e-11_r8 * exp_fac(:)
      rate(:,391) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,399) = 2.3e-12_r8 * exp_fac(:)
      rate(:,499) = 2.3e-12_r8 * exp_fac(:)
      rate(:,403) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,422) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,430) = 1.7e-12_r8 * exp_fac(:)
      rate(:,516) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,443) = 1.2e-12_r8 * exp_fac(:)
      rate(:,510) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,444) = 6.3e-16_r8 * exp_fac(:)
      rate(:,513) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,445) = 1.2e-11_r8 * exp_fac(:)
      rate(:,514) = 1.2e-11_r8 * exp_fac(:)
      rate(:,463) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,464) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,465) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,466) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,474) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,477) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,480) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,139), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,149), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,159), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,167), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,171), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,172), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,190), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,210), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,221), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,267), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,277), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,278), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,279), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,305), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,306), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,326), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,352), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,355), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,413), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,416), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,419), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,426), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,468), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      integer   ::  k
      real(r8)  :: itemp(ncol*kbot)
      real(r8)  :: exp_fac(ncol*kbot)
      real(r8)  :: ko(ncol*kbot)
      real(r8)  :: kinf(ncol*kbot)
      real(r8)  :: wrk(ncol*kbot)
 
      n = ncol*kbot

      rate(:n,136) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,131) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,140) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,141) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,142) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,145) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,146) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,147) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,152) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,156) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,157) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,165) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,166) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,139) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
