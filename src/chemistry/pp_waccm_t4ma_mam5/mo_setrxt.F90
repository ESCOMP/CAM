
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
      rate(:,156) = 5e-12_r8
      rate(:,164) = 3.5e-12_r8
      rate(:,166) = 1.3e-11_r8
      rate(:,167) = 2.2e-11_r8
      rate(:,168) = 5e-11_r8
      rate(:,206) = 1.7e-13_r8
      rate(:,208) = 2.607e-10_r8
      rate(:,209) = 9.75e-11_r8
      rate(:,210) = 2.07e-10_r8
      rate(:,211) = 2.088e-10_r8
      rate(:,212) = 1.17e-10_r8
      rate(:,213) = 4.644e-11_r8
      rate(:,214) = 1.204e-10_r8
      rate(:,215) = 9.9e-11_r8
      rate(:,216) = 3.3e-12_r8
      rate(:,235) = 4.5e-11_r8
      rate(:,236) = 4.62e-10_r8
      rate(:,237) = 1.2e-10_r8
      rate(:,238) = 9e-11_r8
      rate(:,239) = 3e-11_r8
      rate(:,244) = 2.14e-11_r8
      rate(:,245) = 1.9e-10_r8
      rate(:,258) = 2.57e-10_r8
      rate(:,259) = 1.8e-10_r8
      rate(:,260) = 1.794e-10_r8
      rate(:,261) = 1.3e-10_r8
      rate(:,262) = 7.65e-11_r8
      rate(:,273) = 1.31e-10_r8
      rate(:,274) = 3.5e-11_r8
      rate(:,275) = 9e-12_r8
      rate(:,279) = 6.8e-14_r8
      rate(:,280) = 2e-13_r8
      rate(:,294) = 1e-12_r8
      rate(:,298) = 1e-14_r8
      rate(:,299) = 1e-11_r8
      rate(:,300) = 1.15e-11_r8
      rate(:,301) = 4e-14_r8
      rate(:,314) = 3e-12_r8
      rate(:,315) = 6.7e-13_r8
      rate(:,325) = 1.4e-11_r8
      rate(:,328) = 2.4e-12_r8
      rate(:,339) = 5e-12_r8
      rate(:,345) = 3.5e-12_r8
      rate(:,350) = 2.4e-12_r8
      rate(:,351) = 1.4e-11_r8
      rate(:,355) = 2.4e-12_r8
      rate(:,360) = 4.5e-11_r8
      rate(:,365) = 2.4e-12_r8
      rate(:,374) = 2.3e-12_r8
      rate(:,376) = 1.2e-11_r8
      rate(:,377) = 5.7e-11_r8
      rate(:,378) = 2.8e-11_r8
      rate(:,379) = 6.6e-11_r8
      rate(:,380) = 1.4e-11_r8
      rate(:,383) = 1.9e-12_r8
      rate(:,390) = 6.34e-08_r8
      rate(:,394) = 1.157e-05_r8
      rate(:,412) = 0.047_r8
      rate(:,413) = 7.7e-05_r8
      rate(:,414) = 0.171_r8
      rate(:,418) = 6e-11_r8
      rate(:,421) = 1e-12_r8
      rate(:,422) = 4e-10_r8
      rate(:,423) = 2e-10_r8
      rate(:,424) = 1e-10_r8
      rate(:,425) = 5e-16_r8
      rate(:,426) = 4.4e-10_r8
      rate(:,427) = 9e-10_r8
      rate(:,429) = 1.3e-10_r8
      rate(:,432) = 8e-10_r8
      rate(:,433) = 5e-12_r8
      rate(:,434) = 7e-10_r8
      rate(:,437) = 4.8e-10_r8
      rate(:,438) = 1e-10_r8
      rate(:,439) = 4e-10_r8
 
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
      rate(:,356) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,141) = 3e-11_r8 * exp_fac(:)
      rate(:,233) = 5.5e-12_r8 * exp_fac(:)
      rate(:,271) = 3.8e-12_r8 * exp_fac(:)
      rate(:,284) = 3.8e-12_r8 * exp_fac(:)
      rate(:,310) = 3.8e-12_r8 * exp_fac(:)
      rate(:,318) = 3.8e-12_r8 * exp_fac(:)
      rate(:,322) = 3.8e-12_r8 * exp_fac(:)
      rate(:,333) = 2.3e-11_r8 * exp_fac(:)
      rate(:,358) = 1.52e-11_r8 * exp_fac(:)
      rate(:,366) = 1.52e-12_r8 * exp_fac(:)
      rate(:,142) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,143) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,144) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,146) = 4.8e-11_r8 * exp_fac(:)
      rate(:,231) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,147) = 1.8e-11_r8 * exp_fac(:)
      rate(:,296) = 4.2e-12_r8 * exp_fac(:)
      rate(:,317) = 4.2e-12_r8 * exp_fac(:)
      rate(:,354) = 4.4e-12_r8 * exp_fac(:)
      rate(:,148) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,153) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,157) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,158) = 2.9e-12_r8 * exp_fac(:)
      rate(:,159) = 1.45e-12_r8 * exp_fac(:)
      rate(:,160) = 1.45e-12_r8 * exp_fac(:)
      rate(:,161) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,162) = 1.2e-13_r8 * exp_fac(:)
      rate(:,191) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,165) = 1.7e-11_r8 * exp_fac(:)
      rate(:,265) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,169) = 3.44e-12_r8 * exp_fac(:)
      rate(:,224) = 2.3e-12_r8 * exp_fac(:)
      rate(:,227) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,170) = 3e-12_r8 * exp_fac(:)
      rate(:,232) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,172) = 7.26e-11_r8 * exp_fac(:)
      rate(:,173) = 4.64e-11_r8 * exp_fac(:)
      rate(:,183) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,184) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,185) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,186) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,187) = 1.4e-11_r8 * exp_fac(:)
      rate(:,201) = 7.4e-12_r8 * exp_fac(:)
      rate(:,292) = 8.1e-12_r8 * exp_fac(:)
      rate(:,188) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,189) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,190) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,192) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,193) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,194) = 2.6e-12_r8 * exp_fac(:)
      rate(:,195) = 6.4e-12_r8 * exp_fac(:)
      rate(:,225) = 4.1e-13_r8 * exp_fac(:)
      rate(:,196) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,198) = 3.6e-12_r8 * exp_fac(:)
      rate(:,247) = 2e-12_r8 * exp_fac(:)
      rate(:,199) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,200) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,202) = 6e-13_r8 * exp_fac(:)
      rate(:,222) = 1.5e-12_r8 * exp_fac(:)
      rate(:,230) = 1.9e-11_r8 * exp_fac(:)
      rate(:,203) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,204) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,205) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,207) = 3e-12_r8 * exp_fac(:)
      rate(:,241) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,219) = 1.7e-11_r8 * exp_fac(:)
      rate(:,246) = 6.3e-12_r8 * exp_fac(:)
      rate(:,220) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,221) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,223) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,226) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,229) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,234) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,240) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,242) = 1.4e-11_r8 * exp_fac(:)
      rate(:,244) = 2.14e-11_r8 * exp_fac(:)
      rate(:,245) = 1.9e-10_r8 * exp_fac(:)
      rate(:,258) = 2.57e-10_r8 * exp_fac(:)
      rate(:,259) = 1.8e-10_r8 * exp_fac(:)
      rate(:,260) = 1.794e-10_r8 * exp_fac(:)
      rate(:,261) = 1.3e-10_r8 * exp_fac(:)
      rate(:,262) = 7.65e-11_r8 * exp_fac(:)
      rate(:,273) = 1.31e-10_r8 * exp_fac(:)
      rate(:,274) = 3.5e-11_r8 * exp_fac(:)
      rate(:,275) = 9e-12_r8 * exp_fac(:)
      rate(:,279) = 6.8e-14_r8 * exp_fac(:)
      rate(:,280) = 2e-13_r8 * exp_fac(:)
      rate(:,294) = 1e-12_r8 * exp_fac(:)
      rate(:,298) = 1e-14_r8 * exp_fac(:)
      rate(:,299) = 1e-11_r8 * exp_fac(:)
      rate(:,300) = 1.15e-11_r8 * exp_fac(:)
      rate(:,301) = 4e-14_r8 * exp_fac(:)
      rate(:,314) = 3e-12_r8 * exp_fac(:)
      rate(:,315) = 6.7e-13_r8 * exp_fac(:)
      rate(:,325) = 1.4e-11_r8 * exp_fac(:)
      rate(:,328) = 2.4e-12_r8 * exp_fac(:)
      rate(:,339) = 5e-12_r8 * exp_fac(:)
      rate(:,345) = 3.5e-12_r8 * exp_fac(:)
      rate(:,350) = 2.4e-12_r8 * exp_fac(:)
      rate(:,351) = 1.4e-11_r8 * exp_fac(:)
      rate(:,355) = 2.4e-12_r8 * exp_fac(:)
      rate(:,360) = 4.5e-11_r8 * exp_fac(:)
      rate(:,365) = 2.4e-12_r8 * exp_fac(:)
      rate(:,374) = 2.3e-12_r8 * exp_fac(:)
      rate(:,376) = 1.2e-11_r8 * exp_fac(:)
      rate(:,377) = 5.7e-11_r8 * exp_fac(:)
      rate(:,378) = 2.8e-11_r8 * exp_fac(:)
      rate(:,379) = 6.6e-11_r8 * exp_fac(:)
      rate(:,380) = 1.4e-11_r8 * exp_fac(:)
      rate(:,383) = 1.9e-12_r8 * exp_fac(:)
      rate(:,390) = 6.34e-08_r8 * exp_fac(:)
      rate(:,394) = 1.157e-05_r8 * exp_fac(:)
      rate(:,412) = 0.047_r8 * exp_fac(:)
      rate(:,413) = 7.7e-05_r8 * exp_fac(:)
      rate(:,414) = 0.171_r8 * exp_fac(:)
      rate(:,418) = 6e-11_r8 * exp_fac(:)
      rate(:,421) = 1e-12_r8 * exp_fac(:)
      rate(:,422) = 4e-10_r8 * exp_fac(:)
      rate(:,423) = 2e-10_r8 * exp_fac(:)
      rate(:,424) = 1e-10_r8 * exp_fac(:)
      rate(:,425) = 5e-16_r8 * exp_fac(:)
      rate(:,426) = 4.4e-10_r8 * exp_fac(:)
      rate(:,427) = 9e-10_r8 * exp_fac(:)
      rate(:,429) = 1.3e-10_r8 * exp_fac(:)
      rate(:,432) = 8e-10_r8 * exp_fac(:)
      rate(:,433) = 5e-12_r8 * exp_fac(:)
      rate(:,434) = 7e-10_r8 * exp_fac(:)
      rate(:,437) = 4.8e-10_r8 * exp_fac(:)
      rate(:,438) = 1e-10_r8 * exp_fac(:)
      rate(:,439) = 4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,243) = 6e-12_r8 * exp_fac(:)
      rate(:,326) = 5e-13_r8 * exp_fac(:)
      rate(:,352) = 5e-13_r8 * exp_fac(:)
      rate(:,362) = 5e-13_r8 * exp_fac(:)
      rate(:,248) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,249) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,250) = 1.64e-12_r8 * exp_fac(:)
      rate(:,341) = 8.5e-16_r8 * exp_fac(:)
      rate(:,251) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,252) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,253) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,254) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,255) = 1.25e-12_r8 * exp_fac(:)
      rate(:,264) = 3.4e-11_r8 * exp_fac(:)
      rate(:,256) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,257) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,263) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,266) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,267) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,268) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,269) = 2.8e-12_r8 * exp_fac(:)
      rate(:,321) = 2.9e-12_r8 * exp_fac(:)
      rate(:,270) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,272) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,278) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,281) = 7.5e-13_r8 * exp_fac(:)
      rate(:,295) = 7.5e-13_r8 * exp_fac(:)
      rate(:,308) = 7.5e-13_r8 * exp_fac(:)
      rate(:,316) = 7.5e-13_r8 * exp_fac(:)
      rate(:,320) = 8.6e-13_r8 * exp_fac(:)
      rate(:,327) = 8e-13_r8 * exp_fac(:)
      rate(:,348) = 8e-13_r8 * exp_fac(:)
      rate(:,353) = 8e-13_r8 * exp_fac(:)
      rate(:,363) = 8e-13_r8 * exp_fac(:)
      rate(:,282) = 2.6e-12_r8 * exp( 365._r8 * itemp(:) )
      rate(:,283) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,285) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,286) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,287) = 1.4e-12_r8 * exp_fac(:)
      rate(:,306) = 6.5e-15_r8 * exp_fac(:)
      rate(:,288) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,289) = 2.9e-12_r8 * exp_fac(:)
      rate(:,290) = 2e-12_r8 * exp_fac(:)
      rate(:,319) = 7.1e-13_r8 * exp_fac(:)
      rate(:,335) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,291) = 4.3e-13_r8 * exp_fac(:)
      rate(:,336) = 4.3e-13_r8 * exp_fac(:)
      rate(:,293) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,297) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,305) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,307) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,309) = 2.7e-12_r8 * exp_fac(:)
      rate(:,329) = 2.7e-12_r8 * exp_fac(:)
      rate(:,330) = 1.3e-13_r8 * exp_fac(:)
      rate(:,332) = 9.6e-12_r8 * exp_fac(:)
      rate(:,338) = 5.3e-12_r8 * exp_fac(:)
      rate(:,349) = 2.7e-12_r8 * exp_fac(:)
      rate(:,364) = 2.7e-12_r8 * exp_fac(:)
      rate(:,311) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      exp_fac(:) = exp( -1860._r8 * itemp(:) )
      rate(:,312) = 1.4e-12_r8 * exp_fac(:)
      rate(:,359) = 1.4e-12_r8 * exp_fac(:)
      rate(:,313) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,331) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,334) = 4.6e-12_r8 * exp_fac(:)
      rate(:,337) = 2.3e-12_r8 * exp_fac(:)
      rate(:,342) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,346) = 1.86e-11_r8 * exp( 175._r8 * itemp(:) )
      rate(:,347) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,357) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,361) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,367) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,368) = 6.3e-16_r8 * exp( -580._r8 * itemp(:) )
      rate(:,369) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      rate(:,370) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,371) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,372) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,373) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,381) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,382) = 3.4e-12_r8 * exp( -1100._r8 * itemp(:) )
      rate(:,384) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,387) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

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
      call jpl( rate(:,163), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,171), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,174), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,175), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,176), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,197), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,217), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,228), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,277), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,302), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,303), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,323), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,340), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,343), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,375), m, 0.6_r8, ko, kinf, n )

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
      rate(:n,156) = 5e-12_r8
      rate(:n,412) = 0.047_r8
      rate(:n,413) = 7.7e-05_r8
      rate(:n,414) = 0.171_r8
      rate(:n,418) = 6e-11_r8
      rate(:n,421) = 1e-12_r8
      rate(:n,422) = 4e-10_r8
      rate(:n,423) = 2e-10_r8
      rate(:n,424) = 1e-10_r8
      rate(:n,426) = 4.4e-10_r8
      rate(:n,429) = 1.3e-10_r8
      rate(:n,432) = 8e-10_r8
      rate(:n,433) = 5e-12_r8
      rate(:n,434) = 7e-10_r8
      rate(:n,437) = 4.8e-10_r8
      rate(:n,438) = 1e-10_r8
      rate(:n,439) = 4e-10_r8
 
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
      rate(:n,157) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,161) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,169) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,170) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,140) = wrk(:)

















      end subroutine setrxt_hrates

      end module mo_setrxt
