
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

      rate(:,115) = 0.000258_r8
      rate(:,116) = 0.085_r8
      rate(:,117) = 1.2e-10_r8
      rate(:,122) = 1.2e-10_r8
      rate(:,123) = 1.2e-10_r8
      rate(:,124) = 1e-20_r8
      rate(:,125) = 1.3e-16_r8
      rate(:,127) = 4.2e-13_r8
      rate(:,129) = 8e-14_r8
      rate(:,130) = 3.9e-17_r8
      rate(:,137) = 6.9e-12_r8
      rate(:,138) = 7.2e-11_r8
      rate(:,139) = 1.6e-12_r8
      rate(:,145) = 1.8e-12_r8
      rate(:,149) = 1.8e-12_r8
      rate(:,152) = 1.06e-05_r8
      rate(:,154) = 7e-11_r8
      rate(:,155) = 7e-13_r8
      rate(:,163) = 3.5e-12_r8
      rate(:,165) = 1.3e-11_r8
      rate(:,166) = 2.2e-11_r8
      rate(:,167) = 5e-11_r8
      rate(:,205) = 1.7e-13_r8
      rate(:,207) = 2.607e-10_r8
      rate(:,208) = 9.75e-11_r8
      rate(:,209) = 2.07e-10_r8
      rate(:,210) = 2.088e-10_r8
      rate(:,211) = 1.17e-10_r8
      rate(:,212) = 4.644e-11_r8
      rate(:,213) = 1.204e-10_r8
      rate(:,214) = 9.9e-11_r8
      rate(:,215) = 3.3e-12_r8
      rate(:,234) = 4.5e-11_r8
      rate(:,235) = 4.62e-10_r8
      rate(:,236) = 1.2e-10_r8
      rate(:,237) = 9e-11_r8
      rate(:,238) = 3e-11_r8
      rate(:,243) = 2.14e-11_r8
      rate(:,244) = 1.9e-10_r8
      rate(:,257) = 2.57e-10_r8
      rate(:,258) = 1.8e-10_r8
      rate(:,259) = 1.794e-10_r8
      rate(:,260) = 1.3e-10_r8
      rate(:,261) = 7.65e-11_r8
      rate(:,272) = 1.31e-10_r8
      rate(:,273) = 3.5e-11_r8
      rate(:,274) = 9e-12_r8
      rate(:,278) = 6.8e-14_r8
      rate(:,279) = 2e-13_r8
      rate(:,293) = 1e-12_r8
      rate(:,297) = 1e-14_r8
      rate(:,298) = 1e-11_r8
      rate(:,299) = 1.15e-11_r8
      rate(:,300) = 4e-14_r8
      rate(:,313) = 3e-12_r8
      rate(:,314) = 6.7e-13_r8
      rate(:,324) = 1.4e-11_r8
      rate(:,327) = 2.4e-12_r8
      rate(:,338) = 5e-12_r8
      rate(:,344) = 3.5e-12_r8
      rate(:,349) = 2.4e-12_r8
      rate(:,350) = 1.4e-11_r8
      rate(:,354) = 2.4e-12_r8
      rate(:,359) = 4.5e-11_r8
      rate(:,364) = 2.4e-12_r8
      rate(:,373) = 2.3e-12_r8
      rate(:,375) = 1.2e-11_r8
      rate(:,376) = 5.7e-11_r8
      rate(:,377) = 2.8e-11_r8
      rate(:,378) = 6.6e-11_r8
      rate(:,379) = 1.4e-11_r8
      rate(:,382) = 1.9e-12_r8
      rate(:,389) = 6.34e-08_r8
      rate(:,393) = 1.157e-05_r8
      rate(:,411) = 0.047_r8
      rate(:,412) = 7.7e-05_r8
      rate(:,413) = 0.171_r8
      rate(:,417) = 6e-11_r8
      rate(:,420) = 1e-12_r8
      rate(:,421) = 4e-10_r8
      rate(:,422) = 2e-10_r8
      rate(:,423) = 1e-10_r8
      rate(:,424) = 5e-16_r8
      rate(:,425) = 4.4e-10_r8
      rate(:,426) = 9e-10_r8
      rate(:,428) = 1.3e-10_r8
      rate(:,431) = 8e-10_r8
      rate(:,432) = 5e-12_r8
      rate(:,433) = 7e-10_r8
      rate(:,436) = 4.8e-10_r8
      rate(:,437) = 1e-10_r8
      rate(:,438) = 4e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,118) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,119) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,120) = 2.64e-11_r8 * exp_fac(:)
      rate(:,121) = 6.6e-12_r8 * exp_fac(:)
      rate(:,126) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,128) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,131) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,132) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,135) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,136) = 1.4e-12_r8 * exp_fac(:)
      rate(:,355) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,141) = 3e-11_r8 * exp_fac(:)
      rate(:,232) = 5.5e-12_r8 * exp_fac(:)
      rate(:,270) = 3.8e-12_r8 * exp_fac(:)
      rate(:,283) = 3.8e-12_r8 * exp_fac(:)
      rate(:,309) = 3.8e-12_r8 * exp_fac(:)
      rate(:,317) = 3.8e-12_r8 * exp_fac(:)
      rate(:,321) = 3.8e-12_r8 * exp_fac(:)
      rate(:,332) = 2.3e-11_r8 * exp_fac(:)
      rate(:,357) = 1.52e-11_r8 * exp_fac(:)
      rate(:,365) = 1.52e-12_r8 * exp_fac(:)
      rate(:,142) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,143) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,144) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,146) = 4.8e-11_r8 * exp_fac(:)
      rate(:,230) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,147) = 1.8e-11_r8 * exp_fac(:)
      rate(:,295) = 4.2e-12_r8 * exp_fac(:)
      rate(:,316) = 4.2e-12_r8 * exp_fac(:)
      rate(:,353) = 4.4e-12_r8 * exp_fac(:)
      rate(:,148) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,153) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,156) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,157) = 2.9e-12_r8 * exp_fac(:)
      rate(:,158) = 1.45e-12_r8 * exp_fac(:)
      rate(:,159) = 1.45e-12_r8 * exp_fac(:)
      rate(:,160) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,161) = 1.2e-13_r8 * exp_fac(:)
      rate(:,190) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,164) = 1.7e-11_r8 * exp_fac(:)
      rate(:,264) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,168) = 3.44e-12_r8 * exp_fac(:)
      rate(:,223) = 2.3e-12_r8 * exp_fac(:)
      rate(:,226) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,169) = 3e-12_r8 * exp_fac(:)
      rate(:,231) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,171) = 7.26e-11_r8 * exp_fac(:)
      rate(:,172) = 4.64e-11_r8 * exp_fac(:)
      rate(:,182) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,183) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,184) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,185) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,186) = 1.4e-11_r8 * exp_fac(:)
      rate(:,200) = 7.4e-12_r8 * exp_fac(:)
      rate(:,291) = 8.1e-12_r8 * exp_fac(:)
      rate(:,187) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,188) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,189) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,191) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,192) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,193) = 2.6e-12_r8 * exp_fac(:)
      rate(:,194) = 6.4e-12_r8 * exp_fac(:)
      rate(:,224) = 4.1e-13_r8 * exp_fac(:)
      rate(:,195) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,197) = 3.6e-12_r8 * exp_fac(:)
      rate(:,246) = 2e-12_r8 * exp_fac(:)
      rate(:,198) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,199) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,201) = 6e-13_r8 * exp_fac(:)
      rate(:,221) = 1.5e-12_r8 * exp_fac(:)
      rate(:,229) = 1.9e-11_r8 * exp_fac(:)
      rate(:,202) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,203) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,204) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,206) = 3e-12_r8 * exp_fac(:)
      rate(:,240) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,218) = 1.7e-11_r8 * exp_fac(:)
      rate(:,245) = 6.3e-12_r8 * exp_fac(:)
      rate(:,219) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,220) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,222) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,225) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,228) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,233) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,239) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,241) = 1.4e-11_r8 * exp_fac(:)
      rate(:,243) = 2.14e-11_r8 * exp_fac(:)
      rate(:,244) = 1.9e-10_r8 * exp_fac(:)
      rate(:,257) = 2.57e-10_r8 * exp_fac(:)
      rate(:,258) = 1.8e-10_r8 * exp_fac(:)
      rate(:,259) = 1.794e-10_r8 * exp_fac(:)
      rate(:,260) = 1.3e-10_r8 * exp_fac(:)
      rate(:,261) = 7.65e-11_r8 * exp_fac(:)
      rate(:,272) = 1.31e-10_r8 * exp_fac(:)
      rate(:,273) = 3.5e-11_r8 * exp_fac(:)
      rate(:,274) = 9e-12_r8 * exp_fac(:)
      rate(:,278) = 6.8e-14_r8 * exp_fac(:)
      rate(:,279) = 2e-13_r8 * exp_fac(:)
      rate(:,293) = 1e-12_r8 * exp_fac(:)
      rate(:,297) = 1e-14_r8 * exp_fac(:)
      rate(:,298) = 1e-11_r8 * exp_fac(:)
      rate(:,299) = 1.15e-11_r8 * exp_fac(:)
      rate(:,300) = 4e-14_r8 * exp_fac(:)
      rate(:,313) = 3e-12_r8 * exp_fac(:)
      rate(:,314) = 6.7e-13_r8 * exp_fac(:)
      rate(:,324) = 1.4e-11_r8 * exp_fac(:)
      rate(:,327) = 2.4e-12_r8 * exp_fac(:)
      rate(:,338) = 5e-12_r8 * exp_fac(:)
      rate(:,344) = 3.5e-12_r8 * exp_fac(:)
      rate(:,349) = 2.4e-12_r8 * exp_fac(:)
      rate(:,350) = 1.4e-11_r8 * exp_fac(:)
      rate(:,354) = 2.4e-12_r8 * exp_fac(:)
      rate(:,359) = 4.5e-11_r8 * exp_fac(:)
      rate(:,364) = 2.4e-12_r8 * exp_fac(:)
      rate(:,373) = 2.3e-12_r8 * exp_fac(:)
      rate(:,375) = 1.2e-11_r8 * exp_fac(:)
      rate(:,376) = 5.7e-11_r8 * exp_fac(:)
      rate(:,377) = 2.8e-11_r8 * exp_fac(:)
      rate(:,378) = 6.6e-11_r8 * exp_fac(:)
      rate(:,379) = 1.4e-11_r8 * exp_fac(:)
      rate(:,382) = 1.9e-12_r8 * exp_fac(:)
      rate(:,389) = 6.34e-08_r8 * exp_fac(:)
      rate(:,393) = 1.157e-05_r8 * exp_fac(:)
      rate(:,411) = 0.047_r8 * exp_fac(:)
      rate(:,412) = 7.7e-05_r8 * exp_fac(:)
      rate(:,413) = 0.171_r8 * exp_fac(:)
      rate(:,417) = 6e-11_r8 * exp_fac(:)
      rate(:,420) = 1e-12_r8 * exp_fac(:)
      rate(:,421) = 4e-10_r8 * exp_fac(:)
      rate(:,422) = 2e-10_r8 * exp_fac(:)
      rate(:,423) = 1e-10_r8 * exp_fac(:)
      rate(:,424) = 5e-16_r8 * exp_fac(:)
      rate(:,425) = 4.4e-10_r8 * exp_fac(:)
      rate(:,426) = 9e-10_r8 * exp_fac(:)
      rate(:,428) = 1.3e-10_r8 * exp_fac(:)
      rate(:,431) = 8e-10_r8 * exp_fac(:)
      rate(:,432) = 5e-12_r8 * exp_fac(:)
      rate(:,433) = 7e-10_r8 * exp_fac(:)
      rate(:,436) = 4.8e-10_r8 * exp_fac(:)
      rate(:,437) = 1e-10_r8 * exp_fac(:)
      rate(:,438) = 4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,242) = 6e-12_r8 * exp_fac(:)
      rate(:,325) = 5e-13_r8 * exp_fac(:)
      rate(:,351) = 5e-13_r8 * exp_fac(:)
      rate(:,361) = 5e-13_r8 * exp_fac(:)
      rate(:,247) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,248) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,249) = 1.64e-12_r8 * exp_fac(:)
      rate(:,340) = 8.5e-16_r8 * exp_fac(:)
      rate(:,250) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,251) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,252) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,253) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,254) = 1.25e-12_r8 * exp_fac(:)
      rate(:,263) = 3.4e-11_r8 * exp_fac(:)
      rate(:,255) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,256) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,262) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,265) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,266) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,267) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,268) = 2.8e-12_r8 * exp_fac(:)
      rate(:,320) = 2.9e-12_r8 * exp_fac(:)
      rate(:,269) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,271) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,277) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,280) = 7.5e-13_r8 * exp_fac(:)
      rate(:,294) = 7.5e-13_r8 * exp_fac(:)
      rate(:,307) = 7.5e-13_r8 * exp_fac(:)
      rate(:,315) = 7.5e-13_r8 * exp_fac(:)
      rate(:,319) = 8.6e-13_r8 * exp_fac(:)
      rate(:,326) = 8e-13_r8 * exp_fac(:)
      rate(:,347) = 8e-13_r8 * exp_fac(:)
      rate(:,352) = 8e-13_r8 * exp_fac(:)
      rate(:,362) = 8e-13_r8 * exp_fac(:)
      rate(:,281) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,282) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,284) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,285) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,286) = 1.4e-12_r8 * exp_fac(:)
      rate(:,305) = 6.5e-15_r8 * exp_fac(:)
      rate(:,287) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,288) = 2.9e-12_r8 * exp_fac(:)
      rate(:,289) = 2e-12_r8 * exp_fac(:)
      rate(:,318) = 7.1e-13_r8 * exp_fac(:)
      rate(:,334) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,290) = 4.3e-13_r8 * exp_fac(:)
      rate(:,335) = 4.3e-13_r8 * exp_fac(:)
      rate(:,292) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,296) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,304) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,306) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,308) = 2.7e-12_r8 * exp_fac(:)
      rate(:,328) = 2.7e-12_r8 * exp_fac(:)
      rate(:,329) = 1.3e-13_r8 * exp_fac(:)
      rate(:,331) = 9.6e-12_r8 * exp_fac(:)
      rate(:,337) = 5.3e-12_r8 * exp_fac(:)
      rate(:,348) = 2.7e-12_r8 * exp_fac(:)
      rate(:,363) = 2.7e-12_r8 * exp_fac(:)
      rate(:,310) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,311) = 1.4e-12_r8 * exp_fac(:)
      rate(:,358) = 1.4e-12_r8 * exp_fac(:)
      rate(:,312) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,330) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,333) = 4.6e-12_r8 * exp_fac(:)
      rate(:,336) = 2.3e-12_r8 * exp_fac(:)
      rate(:,341) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,345) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,346) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,356) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,360) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,366) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,367) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,368) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,369) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,370) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,371) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,372) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,380) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,381) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,383) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,386) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,140), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,150), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,162), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,173), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,174), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,175), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,196), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,216), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,227), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,276), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,301), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,302), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,322), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,339), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,342), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,374), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,124) = 1e-20_r8
      rate(:n,125) = 1.3e-16_r8
      rate(:n,129) = 8e-14_r8
      rate(:n,130) = 3.9e-17_r8
      rate(:n,137) = 6.9e-12_r8
      rate(:n,154) = 7e-11_r8
      rate(:n,155) = 7e-13_r8
      rate(:n,411) = 0.047_r8
      rate(:n,412) = 7.7e-05_r8
      rate(:n,413) = 0.171_r8
      rate(:n,417) = 6e-11_r8
      rate(:n,420) = 1e-12_r8
      rate(:n,421) = 4e-10_r8
      rate(:n,422) = 2e-10_r8
      rate(:n,423) = 1e-10_r8
      rate(:n,425) = 4.4e-10_r8
      rate(:n,428) = 1.3e-10_r8
      rate(:n,431) = 8e-10_r8
      rate(:n,432) = 5e-12_r8
      rate(:n,433) = 7e-10_r8
      rate(:n,436) = 4.8e-10_r8
      rate(:n,437) = 1e-10_r8
      rate(:n,438) = 4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,119) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,120) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,121) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,126) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,128) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,131) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,132) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,141) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,142) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,143) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,146) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,147) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,148) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,156) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,160) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,168) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,169) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,140) = wrk(:)

















      end subroutine setrxt_hrates

      end module mo_setrxt
