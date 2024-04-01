
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

      rate(:,144) = 1.2e-10_r8
      rate(:,148) = 1.2e-10_r8
      rate(:,154) = 6.9e-12_r8
      rate(:,155) = 7.2e-11_r8
      rate(:,156) = 1.6e-12_r8
      rate(:,162) = 1.8e-12_r8
      rate(:,166) = 1.8e-12_r8
      rate(:,178) = 3.5e-12_r8
      rate(:,180) = 1.3e-11_r8
      rate(:,181) = 2.2e-11_r8
      rate(:,182) = 5e-11_r8
      rate(:,217) = 1.7e-13_r8
      rate(:,219) = 2.607e-10_r8
      rate(:,220) = 9.75e-11_r8
      rate(:,221) = 2.07e-10_r8
      rate(:,222) = 2.088e-10_r8
      rate(:,223) = 1.17e-10_r8
      rate(:,224) = 4.644e-11_r8
      rate(:,225) = 1.204e-10_r8
      rate(:,226) = 9.9e-11_r8
      rate(:,227) = 3.3e-12_r8
      rate(:,246) = 4.5e-11_r8
      rate(:,247) = 4.62e-10_r8
      rate(:,248) = 1.2e-10_r8
      rate(:,249) = 9e-11_r8
      rate(:,250) = 3e-11_r8
      rate(:,255) = 2.14e-11_r8
      rate(:,256) = 1.9e-10_r8
      rate(:,269) = 2.57e-10_r8
      rate(:,270) = 1.8e-10_r8
      rate(:,271) = 1.794e-10_r8
      rate(:,272) = 1.3e-10_r8
      rate(:,273) = 7.65e-11_r8
      rate(:,286) = 4e-13_r8
      rate(:,290) = 1.31e-10_r8
      rate(:,291) = 3.5e-11_r8
      rate(:,292) = 9e-12_r8
      rate(:,299) = 6.8e-14_r8
      rate(:,300) = 2e-13_r8
      rate(:,315) = 1e-12_r8
      rate(:,319) = 1e-14_r8
      rate(:,320) = 1e-11_r8
      rate(:,321) = 1.15e-11_r8
      rate(:,322) = 4e-14_r8
      rate(:,335) = 3e-12_r8
      rate(:,336) = 6.7e-13_r8
      rate(:,346) = 3.5e-13_r8
      rate(:,347) = 5.4e-11_r8
      rate(:,350) = 2e-12_r8
      rate(:,351) = 1.4e-11_r8
      rate(:,354) = 2.4e-12_r8
      rate(:,365) = 5e-12_r8
      rate(:,375) = 1.6e-12_r8
      rate(:,377) = 6.7e-12_r8
      rate(:,380) = 3.5e-12_r8
      rate(:,383) = 1.3e-11_r8
      rate(:,384) = 1.4e-11_r8
      rate(:,388) = 2.4e-12_r8
      rate(:,389) = 1.4e-11_r8
      rate(:,394) = 2.4e-12_r8
      rate(:,395) = 4e-11_r8
      rate(:,396) = 4e-11_r8
      rate(:,398) = 1.4e-11_r8
      rate(:,402) = 2.4e-12_r8
      rate(:,403) = 4e-11_r8
      rate(:,407) = 7e-11_r8
      rate(:,408) = 1e-10_r8
      rate(:,413) = 2.4e-12_r8
      rate(:,428) = 4.7e-11_r8
      rate(:,441) = 2.1e-12_r8
      rate(:,442) = 2.8e-13_r8
      rate(:,450) = 1.7e-11_r8
      rate(:,456) = 8.4e-11_r8
      rate(:,458) = 1.9e-11_r8
      rate(:,459) = 1.2e-14_r8
      rate(:,460) = 2e-10_r8
      rate(:,467) = 2.4e-12_r8
      rate(:,468) = 2e-11_r8
      rate(:,472) = 2.3e-11_r8
      rate(:,473) = 2e-11_r8
      rate(:,477) = 3.3e-11_r8
      rate(:,478) = 1e-12_r8
      rate(:,479) = 5.7e-11_r8
      rate(:,480) = 3.4e-11_r8
      rate(:,485) = 2.3e-12_r8
      rate(:,487) = 1.2e-11_r8
      rate(:,488) = 5.7e-11_r8
      rate(:,489) = 2.8e-11_r8
      rate(:,490) = 6.6e-11_r8
      rate(:,491) = 1.4e-11_r8
      rate(:,494) = 1.9e-12_r8
      rate(:,506) = 6.34e-08_r8
      rate(:,512) = 1.9e-11_r8
      rate(:,515) = 1.2e-14_r8
      rate(:,516) = 2e-10_r8
      rate(:,527) = 1.34e-11_r8
      rate(:,530) = 1.34e-11_r8
      rate(:,536) = 1.34e-11_r8
      rate(:,537) = 1.34e-11_r8
      rate(:,542) = 1.7e-11_r8
      rate(:,562) = 1.29e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,145) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,146) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,147) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,149) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,152) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,153) = 1.4e-12_r8 * exp_fac(:)
      rate(:,404) = 1.05e-14_r8 * exp_fac(:)
      rate(:,523) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,158) = 3e-11_r8 * exp_fac(:)
      rate(:,244) = 5.5e-12_r8 * exp_fac(:)
      rate(:,283) = 3.8e-12_r8 * exp_fac(:)
      rate(:,304) = 3.8e-12_r8 * exp_fac(:)
      rate(:,331) = 3.8e-12_r8 * exp_fac(:)
      rate(:,339) = 3.8e-12_r8 * exp_fac(:)
      rate(:,343) = 3.8e-12_r8 * exp_fac(:)
      rate(:,359) = 2.3e-11_r8 * exp_fac(:)
      rate(:,369) = 3.8e-12_r8 * exp_fac(:)
      rate(:,379) = 3.8e-12_r8 * exp_fac(:)
      rate(:,406) = 1.52e-11_r8 * exp_fac(:)
      rate(:,414) = 1.52e-12_r8 * exp_fac(:)
      rate(:,420) = 3.8e-12_r8 * exp_fac(:)
      rate(:,423) = 3.8e-12_r8 * exp_fac(:)
      rate(:,427) = 3.8e-12_r8 * exp_fac(:)
      rate(:,443) = 3.8e-12_r8 * exp_fac(:)
      rate(:,447) = 3.8e-12_r8 * exp_fac(:)
      rate(:,453) = 3.8e-12_r8 * exp_fac(:)
      rate(:,457) = 3.8e-12_r8 * exp_fac(:)
      rate(:,159) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,160) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,161) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,163) = 4.8e-11_r8 * exp_fac(:)
      rate(:,242) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,164) = 1.8e-11_r8 * exp_fac(:)
      rate(:,317) = 4.2e-12_r8 * exp_fac(:)
      rate(:,330) = 4.2e-12_r8 * exp_fac(:)
      rate(:,338) = 4.2e-12_r8 * exp_fac(:)
      rate(:,367) = 4.2e-12_r8 * exp_fac(:)
      rate(:,387) = 4.4e-12_r8 * exp_fac(:)
      rate(:,393) = 4.4e-12_r8 * exp_fac(:)
      rate(:,466) = 4.2e-12_r8 * exp_fac(:)
      rate(:,471) = 4.2e-12_r8 * exp_fac(:)
      rate(:,476) = 4.2e-12_r8 * exp_fac(:)
      rate(:,165) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,169) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,170) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,171) = 2.9e-12_r8 * exp_fac(:)
      rate(:,172) = 1.45e-12_r8 * exp_fac(:)
      rate(:,173) = 1.45e-12_r8 * exp_fac(:)
      rate(:,174) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,175) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,176) = 1.2e-13_r8 * exp_fac(:)
      rate(:,202) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,179) = 1.7e-11_r8 * exp_fac(:)
      rate(:,277) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,183) = 3.44e-12_r8 * exp_fac(:)
      rate(:,235) = 2.3e-12_r8 * exp_fac(:)
      rate(:,238) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,184) = 3e-12_r8 * exp_fac(:)
      rate(:,243) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,186) = 7.26e-11_r8 * exp_fac(:)
      rate(:,187) = 4.64e-11_r8 * exp_fac(:)
      rate(:,194) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,195) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,196) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,197) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,198) = 1.4e-11_r8 * exp_fac(:)
      rate(:,212) = 7.4e-12_r8 * exp_fac(:)
      rate(:,313) = 8.1e-12_r8 * exp_fac(:)
      rate(:,199) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,200) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,201) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,203) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,204) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,205) = 2.6e-12_r8 * exp_fac(:)
      rate(:,206) = 6.4e-12_r8 * exp_fac(:)
      rate(:,236) = 4.1e-13_r8 * exp_fac(:)
      rate(:,416) = 7.5e-12_r8 * exp_fac(:)
      rate(:,430) = 7.5e-12_r8 * exp_fac(:)
      rate(:,433) = 7.5e-12_r8 * exp_fac(:)
      rate(:,436) = 7.5e-12_r8 * exp_fac(:)
      rate(:,207) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,209) = 3.6e-12_r8 * exp_fac(:)
      rate(:,258) = 2e-12_r8 * exp_fac(:)
      rate(:,210) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,211) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,213) = 6e-13_r8 * exp_fac(:)
      rate(:,233) = 1.5e-12_r8 * exp_fac(:)
      rate(:,241) = 1.9e-11_r8 * exp_fac(:)
      rate(:,214) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,215) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,216) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,218) = 3e-12_r8 * exp_fac(:)
      rate(:,252) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,230) = 1.7e-11_r8 * exp_fac(:)
      rate(:,257) = 6.3e-12_r8 * exp_fac(:)
      rate(:,231) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,232) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,234) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,237) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,240) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,245) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,251) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,253) = 1.4e-11_r8 * exp_fac(:)
      rate(:,255) = 2.14e-11_r8 * exp_fac(:)
      rate(:,256) = 1.9e-10_r8 * exp_fac(:)
      rate(:,269) = 2.57e-10_r8 * exp_fac(:)
      rate(:,270) = 1.8e-10_r8 * exp_fac(:)
      rate(:,271) = 1.794e-10_r8 * exp_fac(:)
      rate(:,272) = 1.3e-10_r8 * exp_fac(:)
      rate(:,273) = 7.65e-11_r8 * exp_fac(:)
      rate(:,286) = 4e-13_r8 * exp_fac(:)
      rate(:,290) = 1.31e-10_r8 * exp_fac(:)
      rate(:,291) = 3.5e-11_r8 * exp_fac(:)
      rate(:,292) = 9e-12_r8 * exp_fac(:)
      rate(:,299) = 6.8e-14_r8 * exp_fac(:)
      rate(:,300) = 2e-13_r8 * exp_fac(:)
      rate(:,315) = 1e-12_r8 * exp_fac(:)
      rate(:,319) = 1e-14_r8 * exp_fac(:)
      rate(:,320) = 1e-11_r8 * exp_fac(:)
      rate(:,321) = 1.15e-11_r8 * exp_fac(:)
      rate(:,322) = 4e-14_r8 * exp_fac(:)
      rate(:,335) = 3e-12_r8 * exp_fac(:)
      rate(:,336) = 6.7e-13_r8 * exp_fac(:)
      rate(:,346) = 3.5e-13_r8 * exp_fac(:)
      rate(:,347) = 5.4e-11_r8 * exp_fac(:)
      rate(:,350) = 2e-12_r8 * exp_fac(:)
      rate(:,351) = 1.4e-11_r8 * exp_fac(:)
      rate(:,354) = 2.4e-12_r8 * exp_fac(:)
      rate(:,365) = 5e-12_r8 * exp_fac(:)
      rate(:,375) = 1.6e-12_r8 * exp_fac(:)
      rate(:,377) = 6.7e-12_r8 * exp_fac(:)
      rate(:,380) = 3.5e-12_r8 * exp_fac(:)
      rate(:,383) = 1.3e-11_r8 * exp_fac(:)
      rate(:,384) = 1.4e-11_r8 * exp_fac(:)
      rate(:,388) = 2.4e-12_r8 * exp_fac(:)
      rate(:,389) = 1.4e-11_r8 * exp_fac(:)
      rate(:,394) = 2.4e-12_r8 * exp_fac(:)
      rate(:,395) = 4e-11_r8 * exp_fac(:)
      rate(:,396) = 4e-11_r8 * exp_fac(:)
      rate(:,398) = 1.4e-11_r8 * exp_fac(:)
      rate(:,402) = 2.4e-12_r8 * exp_fac(:)
      rate(:,403) = 4e-11_r8 * exp_fac(:)
      rate(:,407) = 7e-11_r8 * exp_fac(:)
      rate(:,408) = 1e-10_r8 * exp_fac(:)
      rate(:,413) = 2.4e-12_r8 * exp_fac(:)
      rate(:,428) = 4.7e-11_r8 * exp_fac(:)
      rate(:,441) = 2.1e-12_r8 * exp_fac(:)
      rate(:,442) = 2.8e-13_r8 * exp_fac(:)
      rate(:,450) = 1.7e-11_r8 * exp_fac(:)
      rate(:,456) = 8.4e-11_r8 * exp_fac(:)
      rate(:,458) = 1.9e-11_r8 * exp_fac(:)
      rate(:,459) = 1.2e-14_r8 * exp_fac(:)
      rate(:,460) = 2e-10_r8 * exp_fac(:)
      rate(:,467) = 2.4e-12_r8 * exp_fac(:)
      rate(:,468) = 2e-11_r8 * exp_fac(:)
      rate(:,472) = 2.3e-11_r8 * exp_fac(:)
      rate(:,473) = 2e-11_r8 * exp_fac(:)
      rate(:,477) = 3.3e-11_r8 * exp_fac(:)
      rate(:,478) = 1e-12_r8 * exp_fac(:)
      rate(:,479) = 5.7e-11_r8 * exp_fac(:)
      rate(:,480) = 3.4e-11_r8 * exp_fac(:)
      rate(:,485) = 2.3e-12_r8 * exp_fac(:)
      rate(:,487) = 1.2e-11_r8 * exp_fac(:)
      rate(:,488) = 5.7e-11_r8 * exp_fac(:)
      rate(:,489) = 2.8e-11_r8 * exp_fac(:)
      rate(:,490) = 6.6e-11_r8 * exp_fac(:)
      rate(:,491) = 1.4e-11_r8 * exp_fac(:)
      rate(:,494) = 1.9e-12_r8 * exp_fac(:)
      rate(:,506) = 6.34e-08_r8 * exp_fac(:)
      rate(:,512) = 1.9e-11_r8 * exp_fac(:)
      rate(:,515) = 1.2e-14_r8 * exp_fac(:)
      rate(:,516) = 2e-10_r8 * exp_fac(:)
      rate(:,527) = 1.34e-11_r8 * exp_fac(:)
      rate(:,530) = 1.34e-11_r8 * exp_fac(:)
      rate(:,536) = 1.34e-11_r8 * exp_fac(:)
      rate(:,537) = 1.34e-11_r8 * exp_fac(:)
      rate(:,542) = 1.7e-11_r8 * exp_fac(:)
      rate(:,562) = 1.29e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,254) = 6e-12_r8 * exp_fac(:)
      rate(:,352) = 5e-13_r8 * exp_fac(:)
      rate(:,385) = 5e-13_r8 * exp_fac(:)
      rate(:,390) = 5e-13_r8 * exp_fac(:)
      rate(:,399) = 5e-13_r8 * exp_fac(:)
      rate(:,410) = 5e-13_r8 * exp_fac(:)
      rate(:,259) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,260) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,261) = 1.64e-12_r8 * exp_fac(:)
      rate(:,371) = 8.5e-16_r8 * exp_fac(:)
      rate(:,262) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,263) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,264) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,265) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,266) = 1.25e-12_r8 * exp_fac(:)
      rate(:,276) = 3.4e-11_r8 * exp_fac(:)
      rate(:,267) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,268) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,274) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,275) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,278) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,279) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,280) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,281) = 2.8e-12_r8 * exp_fac(:)
      rate(:,342) = 2.9e-12_r8 * exp_fac(:)
      rate(:,282) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,284) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,287) = 7.5e-13_r8 * exp_fac(:)
      rate(:,301) = 7.5e-13_r8 * exp_fac(:)
      rate(:,316) = 7.5e-13_r8 * exp_fac(:)
      rate(:,329) = 7.5e-13_r8 * exp_fac(:)
      rate(:,337) = 7.5e-13_r8 * exp_fac(:)
      rate(:,341) = 8.6e-13_r8 * exp_fac(:)
      rate(:,353) = 8e-13_r8 * exp_fac(:)
      rate(:,366) = 7.5e-13_r8 * exp_fac(:)
      rate(:,376) = 7.5e-13_r8 * exp_fac(:)
      rate(:,386) = 8e-13_r8 * exp_fac(:)
      rate(:,391) = 8e-13_r8 * exp_fac(:)
      rate(:,400) = 8e-13_r8 * exp_fac(:)
      rate(:,411) = 8e-13_r8 * exp_fac(:)
      rate(:,418) = 7.5e-13_r8 * exp_fac(:)
      rate(:,422) = 7.5e-13_r8 * exp_fac(:)
      rate(:,425) = 7.5e-13_r8 * exp_fac(:)
      rate(:,438) = 7.5e-13_r8 * exp_fac(:)
      rate(:,445) = 7.5e-13_r8 * exp_fac(:)
      rate(:,451) = 7.5e-13_r8 * exp_fac(:)
      rate(:,454) = 7.5e-13_r8 * exp_fac(:)
      rate(:,465) = 7.5e-13_r8 * exp_fac(:)
      rate(:,470) = 7.5e-13_r8 * exp_fac(:)
      rate(:,475) = 7.5e-13_r8 * exp_fac(:)
      rate(:,518) = 7.5e-13_r8 * exp_fac(:)
      rate(:,525) = 7.5e-13_r8 * exp_fac(:)
      rate(:,528) = 7.5e-13_r8 * exp_fac(:)
      rate(:,539) = 7.5e-13_r8 * exp_fac(:)
      rate(:,543) = 7.5e-13_r8 * exp_fac(:)
      rate(:,288) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,289) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,293) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,298) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,302) = 2.6e-12_r8 * exp_fac(:)
      rate(:,419) = 2.6e-12_r8 * exp_fac(:)
      rate(:,424) = 2.6e-12_r8 * exp_fac(:)
      rate(:,426) = 2.6e-12_r8 * exp_fac(:)
      rate(:,439) = 2.6e-12_r8 * exp_fac(:)
      rate(:,446) = 2.6e-12_r8 * exp_fac(:)
      rate(:,452) = 2.6e-12_r8 * exp_fac(:)
      rate(:,455) = 2.6e-12_r8 * exp_fac(:)
      rate(:,519) = 2.6e-12_r8 * exp_fac(:)
      rate(:,526) = 2.6e-12_r8 * exp_fac(:)
      rate(:,529) = 2.6e-12_r8 * exp_fac(:)
      rate(:,540) = 2.6e-12_r8 * exp_fac(:)
      rate(:,544) = 2.6e-12_r8 * exp_fac(:)
      rate(:,303) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,305) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,306) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,307) = 1.4e-12_r8 * exp_fac(:)
      rate(:,327) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,308) = 4.63e-12_r8 * exp_fac(:)
      rate(:,522) = 2.7e-12_r8 * exp_fac(:)
      rate(:,309) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,310) = 2.9e-12_r8 * exp_fac(:)
      rate(:,311) = 2e-12_r8 * exp_fac(:)
      rate(:,340) = 7.1e-13_r8 * exp_fac(:)
      rate(:,361) = 2e-12_r8 * exp_fac(:)
      rate(:,464) = 2e-12_r8 * exp_fac(:)
      rate(:,469) = 2e-12_r8 * exp_fac(:)
      rate(:,474) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,312) = 4.3e-13_r8 * exp_fac(:)
      rate(:,362) = 4.3e-13_r8 * exp_fac(:)
      rate(:,415) = 4.3e-13_r8 * exp_fac(:)
      rate(:,429) = 4.3e-13_r8 * exp_fac(:)
      rate(:,432) = 4.3e-13_r8 * exp_fac(:)
      rate(:,435) = 4.3e-13_r8 * exp_fac(:)
      rate(:,314) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,318) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,326) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,328) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,332) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,333) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,334) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,348) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,349) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,355) = 2.7e-12_r8 * exp_fac(:)
      rate(:,356) = 1.3e-13_r8 * exp_fac(:)
      rate(:,358) = 9.6e-12_r8 * exp_fac(:)
      rate(:,364) = 5.3e-12_r8 * exp_fac(:)
      rate(:,401) = 2.7e-12_r8 * exp_fac(:)
      rate(:,412) = 2.7e-12_r8 * exp_fac(:)
      rate(:,514) = 2.7e-12_r8 * exp_fac(:)
      rate(:,533) = 2.7e-12_r8 * exp_fac(:)
      rate(:,357) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,360) = 4.6e-12_r8 * exp_fac(:)
      rate(:,363) = 2.3e-12_r8 * exp_fac(:)
      rate(:,368) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,372) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,378) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,381) = 1.86e-11_r8 * exp_fac(:)
      rate(:,382) = 1.86e-11_r8 * exp_fac(:)
      rate(:,392) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,397) = 3.03e-12_r8 * exp_fac(:)
      rate(:,520) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,405) = 2.54e-11_r8 * exp_fac(:)
      rate(:,524) = 2.54e-11_r8 * exp_fac(:)
      rate(:,409) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,417) = 2.3e-12_r8 * exp_fac(:)
      rate(:,517) = 2.3e-12_r8 * exp_fac(:)
      rate(:,421) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,440) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,448) = 1.7e-12_r8 * exp_fac(:)
      rate(:,538) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,461) = 1.2e-12_r8 * exp_fac(:)
      rate(:,531) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,462) = 6.3e-16_r8 * exp_fac(:)
      rate(:,534) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,463) = 1.2e-11_r8 * exp_fac(:)
      rate(:,535) = 1.2e-11_r8 * exp_fac(:)
      rate(:,481) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,482) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,483) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,484) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,492) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,493) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,495) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,498) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,513) = 2.75e-13_r8 * exp_fac(:)
      rate(:,521) = 2.12e-13_r8 * exp_fac(:)
      rate(:,532) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,157), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,167), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,177), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,185), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,188), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,189), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,190), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,208), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,228), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,239), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,285), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,295), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,296), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,297), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,323), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,324), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,344), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,370), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,431), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,434), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,437), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,444), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,486), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,154) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,146) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,149) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,158) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,159) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,160) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,163) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,164) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,165) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,170) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,174) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,175) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,183) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,184) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,157) = wrk(:)























      end subroutine setrxt_hrates

      end module mo_setrxt
