
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
      rate(:,153) = 7e-13_r8
      rate(:,154) = 5e-12_r8
      rate(:,163) = 3.5e-12_r8
      rate(:,165) = 1.3e-11_r8
      rate(:,166) = 2.2e-11_r8
      rate(:,167) = 5e-11_r8
      rate(:,202) = 1.7e-13_r8
      rate(:,204) = 2.607e-10_r8
      rate(:,205) = 9.75e-11_r8
      rate(:,206) = 2.07e-10_r8
      rate(:,207) = 2.088e-10_r8
      rate(:,208) = 1.17e-10_r8
      rate(:,209) = 4.644e-11_r8
      rate(:,210) = 1.204e-10_r8
      rate(:,211) = 9.9e-11_r8
      rate(:,212) = 3.3e-12_r8
      rate(:,231) = 4.5e-11_r8
      rate(:,232) = 4.62e-10_r8
      rate(:,233) = 1.2e-10_r8
      rate(:,234) = 9e-11_r8
      rate(:,235) = 3e-11_r8
      rate(:,240) = 2.14e-11_r8
      rate(:,241) = 1.9e-10_r8
      rate(:,254) = 2.57e-10_r8
      rate(:,255) = 1.8e-10_r8
      rate(:,256) = 1.794e-10_r8
      rate(:,257) = 1.3e-10_r8
      rate(:,258) = 7.65e-11_r8
      rate(:,269) = 1.31e-10_r8
      rate(:,270) = 3.5e-11_r8
      rate(:,271) = 9e-12_r8
      rate(:,275) = 6.8e-14_r8
      rate(:,276) = 2e-13_r8
      rate(:,290) = 1e-12_r8
      rate(:,294) = 1e-14_r8
      rate(:,295) = 1e-11_r8
      rate(:,296) = 1.15e-11_r8
      rate(:,297) = 4e-14_r8
      rate(:,310) = 3e-12_r8
      rate(:,311) = 6.7e-13_r8
      rate(:,321) = 1.4e-11_r8
      rate(:,324) = 2.4e-12_r8
      rate(:,335) = 5e-12_r8
      rate(:,341) = 3.5e-12_r8
      rate(:,346) = 2.4e-12_r8
      rate(:,347) = 1.4e-11_r8
      rate(:,351) = 2.4e-12_r8
      rate(:,356) = 4.5e-11_r8
      rate(:,361) = 2.4e-12_r8
      rate(:,370) = 2.3e-12_r8
      rate(:,372) = 1.2e-11_r8
      rate(:,373) = 5.7e-11_r8
      rate(:,374) = 2.8e-11_r8
      rate(:,375) = 6.6e-11_r8
      rate(:,376) = 1.4e-11_r8
      rate(:,379) = 1.9e-12_r8
      rate(:,386) = 6.34e-08_r8
      rate(:,390) = 1.157e-05_r8
      rate(:,408) = 0.047_r8
      rate(:,409) = 7.7e-05_r8
      rate(:,410) = 0.171_r8
      rate(:,414) = 6e-11_r8
      rate(:,417) = 1e-12_r8
      rate(:,418) = 4e-10_r8
      rate(:,419) = 2e-10_r8
      rate(:,420) = 1e-10_r8
      rate(:,421) = 5e-16_r8
      rate(:,422) = 4.4e-10_r8
      rate(:,423) = 9e-10_r8
      rate(:,425) = 1.3e-10_r8
      rate(:,428) = 8e-10_r8
      rate(:,429) = 5e-12_r8
      rate(:,430) = 7e-10_r8
      rate(:,433) = 4.8e-10_r8
      rate(:,434) = 1e-10_r8
      rate(:,435) = 4e-10_r8
 
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
      rate(:,352) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,141) = 3e-11_r8 * exp_fac(:)
      rate(:,229) = 5.5e-12_r8 * exp_fac(:)
      rate(:,267) = 3.8e-12_r8 * exp_fac(:)
      rate(:,280) = 3.8e-12_r8 * exp_fac(:)
      rate(:,306) = 3.8e-12_r8 * exp_fac(:)
      rate(:,314) = 3.8e-12_r8 * exp_fac(:)
      rate(:,318) = 3.8e-12_r8 * exp_fac(:)
      rate(:,329) = 2.3e-11_r8 * exp_fac(:)
      rate(:,354) = 1.52e-11_r8 * exp_fac(:)
      rate(:,362) = 1.52e-12_r8 * exp_fac(:)
      rate(:,142) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,143) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,144) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,146) = 4.8e-11_r8 * exp_fac(:)
      rate(:,227) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,147) = 1.8e-11_r8 * exp_fac(:)
      rate(:,292) = 4.2e-12_r8 * exp_fac(:)
      rate(:,305) = 4.2e-12_r8 * exp_fac(:)
      rate(:,313) = 4.2e-12_r8 * exp_fac(:)
      rate(:,350) = 4.4e-12_r8 * exp_fac(:)
      rate(:,148) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,152) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,155) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,156) = 2.9e-12_r8 * exp_fac(:)
      rate(:,157) = 1.45e-12_r8 * exp_fac(:)
      rate(:,158) = 1.45e-12_r8 * exp_fac(:)
      rate(:,159) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,160) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,161) = 1.2e-13_r8 * exp_fac(:)
      rate(:,187) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,164) = 1.7e-11_r8 * exp_fac(:)
      rate(:,261) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,168) = 3.44e-12_r8 * exp_fac(:)
      rate(:,220) = 2.3e-12_r8 * exp_fac(:)
      rate(:,223) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,169) = 3e-12_r8 * exp_fac(:)
      rate(:,228) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,171) = 7.26e-11_r8 * exp_fac(:)
      rate(:,172) = 4.64e-11_r8 * exp_fac(:)
      rate(:,179) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,180) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,181) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,182) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,183) = 1.4e-11_r8 * exp_fac(:)
      rate(:,197) = 7.4e-12_r8 * exp_fac(:)
      rate(:,288) = 8.1e-12_r8 * exp_fac(:)
      rate(:,184) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,185) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,186) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,188) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,189) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,190) = 2.6e-12_r8 * exp_fac(:)
      rate(:,191) = 6.4e-12_r8 * exp_fac(:)
      rate(:,221) = 4.1e-13_r8 * exp_fac(:)
      rate(:,192) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,194) = 3.6e-12_r8 * exp_fac(:)
      rate(:,243) = 2e-12_r8 * exp_fac(:)
      rate(:,195) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,196) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,198) = 6e-13_r8 * exp_fac(:)
      rate(:,218) = 1.5e-12_r8 * exp_fac(:)
      rate(:,226) = 1.9e-11_r8 * exp_fac(:)
      rate(:,199) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,200) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,201) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,203) = 3e-12_r8 * exp_fac(:)
      rate(:,237) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,215) = 1.7e-11_r8 * exp_fac(:)
      rate(:,242) = 6.3e-12_r8 * exp_fac(:)
      rate(:,216) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,217) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,219) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,222) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,225) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,230) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,236) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,238) = 1.4e-11_r8 * exp_fac(:)
      rate(:,240) = 2.14e-11_r8 * exp_fac(:)
      rate(:,241) = 1.9e-10_r8 * exp_fac(:)
      rate(:,254) = 2.57e-10_r8 * exp_fac(:)
      rate(:,255) = 1.8e-10_r8 * exp_fac(:)
      rate(:,256) = 1.794e-10_r8 * exp_fac(:)
      rate(:,257) = 1.3e-10_r8 * exp_fac(:)
      rate(:,258) = 7.65e-11_r8 * exp_fac(:)
      rate(:,269) = 1.31e-10_r8 * exp_fac(:)
      rate(:,270) = 3.5e-11_r8 * exp_fac(:)
      rate(:,271) = 9e-12_r8 * exp_fac(:)
      rate(:,275) = 6.8e-14_r8 * exp_fac(:)
      rate(:,276) = 2e-13_r8 * exp_fac(:)
      rate(:,290) = 1e-12_r8 * exp_fac(:)
      rate(:,294) = 1e-14_r8 * exp_fac(:)
      rate(:,295) = 1e-11_r8 * exp_fac(:)
      rate(:,296) = 1.15e-11_r8 * exp_fac(:)
      rate(:,297) = 4e-14_r8 * exp_fac(:)
      rate(:,310) = 3e-12_r8 * exp_fac(:)
      rate(:,311) = 6.7e-13_r8 * exp_fac(:)
      rate(:,321) = 1.4e-11_r8 * exp_fac(:)
      rate(:,324) = 2.4e-12_r8 * exp_fac(:)
      rate(:,335) = 5e-12_r8 * exp_fac(:)
      rate(:,341) = 3.5e-12_r8 * exp_fac(:)
      rate(:,346) = 2.4e-12_r8 * exp_fac(:)
      rate(:,347) = 1.4e-11_r8 * exp_fac(:)
      rate(:,351) = 2.4e-12_r8 * exp_fac(:)
      rate(:,356) = 4.5e-11_r8 * exp_fac(:)
      rate(:,361) = 2.4e-12_r8 * exp_fac(:)
      rate(:,370) = 2.3e-12_r8 * exp_fac(:)
      rate(:,372) = 1.2e-11_r8 * exp_fac(:)
      rate(:,373) = 5.7e-11_r8 * exp_fac(:)
      rate(:,374) = 2.8e-11_r8 * exp_fac(:)
      rate(:,375) = 6.6e-11_r8 * exp_fac(:)
      rate(:,376) = 1.4e-11_r8 * exp_fac(:)
      rate(:,379) = 1.9e-12_r8 * exp_fac(:)
      rate(:,386) = 6.34e-08_r8 * exp_fac(:)
      rate(:,390) = 1.157e-05_r8 * exp_fac(:)
      rate(:,408) = 0.047_r8 * exp_fac(:)
      rate(:,409) = 7.7e-05_r8 * exp_fac(:)
      rate(:,410) = 0.171_r8 * exp_fac(:)
      rate(:,414) = 6e-11_r8 * exp_fac(:)
      rate(:,417) = 1e-12_r8 * exp_fac(:)
      rate(:,418) = 4e-10_r8 * exp_fac(:)
      rate(:,419) = 2e-10_r8 * exp_fac(:)
      rate(:,420) = 1e-10_r8 * exp_fac(:)
      rate(:,421) = 5e-16_r8 * exp_fac(:)
      rate(:,422) = 4.4e-10_r8 * exp_fac(:)
      rate(:,423) = 9e-10_r8 * exp_fac(:)
      rate(:,425) = 1.3e-10_r8 * exp_fac(:)
      rate(:,428) = 8e-10_r8 * exp_fac(:)
      rate(:,429) = 5e-12_r8 * exp_fac(:)
      rate(:,430) = 7e-10_r8 * exp_fac(:)
      rate(:,433) = 4.8e-10_r8 * exp_fac(:)
      rate(:,434) = 1e-10_r8 * exp_fac(:)
      rate(:,435) = 4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,239) = 6e-12_r8 * exp_fac(:)
      rate(:,322) = 5e-13_r8 * exp_fac(:)
      rate(:,348) = 5e-13_r8 * exp_fac(:)
      rate(:,358) = 5e-13_r8 * exp_fac(:)
      rate(:,244) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,245) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,246) = 1.64e-12_r8 * exp_fac(:)
      rate(:,337) = 8.5e-16_r8 * exp_fac(:)
      rate(:,247) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,248) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,249) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,250) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,251) = 1.25e-12_r8 * exp_fac(:)
      rate(:,260) = 3.4e-11_r8 * exp_fac(:)
      rate(:,252) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,253) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,259) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,262) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,263) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,264) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,265) = 2.8e-12_r8 * exp_fac(:)
      rate(:,317) = 2.9e-12_r8 * exp_fac(:)
      rate(:,266) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,268) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,274) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,277) = 7.5e-13_r8 * exp_fac(:)
      rate(:,291) = 7.5e-13_r8 * exp_fac(:)
      rate(:,304) = 7.5e-13_r8 * exp_fac(:)
      rate(:,312) = 7.5e-13_r8 * exp_fac(:)
      rate(:,316) = 8.6e-13_r8 * exp_fac(:)
      rate(:,323) = 8e-13_r8 * exp_fac(:)
      rate(:,344) = 8e-13_r8 * exp_fac(:)
      rate(:,349) = 8e-13_r8 * exp_fac(:)
      rate(:,359) = 8e-13_r8 * exp_fac(:)
      rate(:,278) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,279) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,281) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,282) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,283) = 1.4e-12_r8 * exp_fac(:)
      rate(:,302) = 6.5e-15_r8 * exp_fac(:)
      rate(:,284) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,285) = 2.9e-12_r8 * exp_fac(:)
      rate(:,286) = 2e-12_r8 * exp_fac(:)
      rate(:,315) = 7.1e-13_r8 * exp_fac(:)
      rate(:,331) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,287) = 4.3e-13_r8 * exp_fac(:)
      rate(:,332) = 4.3e-13_r8 * exp_fac(:)
      rate(:,289) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,293) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,301) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,303) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,307) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,308) = 1.4e-12_r8 * exp_fac(:)
      rate(:,355) = 1.4e-12_r8 * exp_fac(:)
      rate(:,309) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,325) = 2.7e-12_r8 * exp_fac(:)
      rate(:,326) = 1.3e-13_r8 * exp_fac(:)
      rate(:,328) = 9.6e-12_r8 * exp_fac(:)
      rate(:,334) = 5.3e-12_r8 * exp_fac(:)
      rate(:,345) = 2.7e-12_r8 * exp_fac(:)
      rate(:,360) = 2.7e-12_r8 * exp_fac(:)
      rate(:,327) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,330) = 4.6e-12_r8 * exp_fac(:)
      rate(:,333) = 2.3e-12_r8 * exp_fac(:)
      rate(:,338) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,342) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,343) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,353) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,357) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,363) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,364) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,365) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,366) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,367) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,368) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,369) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,377) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,378) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,380) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,383) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

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
      call jpl( rate(:,193), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,213), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,224), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,273), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,298), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,299), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,319), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,336), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,339), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,371), m, 0.6_r8, ko, kinf, n )

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
      rate(:n,153) = 7e-13_r8
      rate(:n,154) = 5e-12_r8
      rate(:n,408) = 0.047_r8
      rate(:n,409) = 7.7e-05_r8
      rate(:n,410) = 0.171_r8
      rate(:n,414) = 6e-11_r8
      rate(:n,417) = 1e-12_r8
      rate(:n,418) = 4e-10_r8
      rate(:n,419) = 2e-10_r8
      rate(:n,420) = 1e-10_r8
      rate(:n,422) = 4.4e-10_r8
      rate(:n,425) = 1.3e-10_r8
      rate(:n,428) = 8e-10_r8
      rate(:n,429) = 5e-12_r8
      rate(:n,430) = 7e-10_r8
      rate(:n,433) = 4.8e-10_r8
      rate(:n,434) = 1e-10_r8
      rate(:n,435) = 4e-10_r8
 
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
      rate(:n,155) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,159) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
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
