
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
      rate(:,125) = 0.000258_r8
      rate(:,126) = 0.085_r8
      rate(:,127) = 1.2e-10_r8
      rate(:,132) = 1.2e-10_r8
      rate(:,133) = 1.2e-10_r8
      rate(:,134) = 1e-20_r8
      rate(:,135) = 1.3e-16_r8
      rate(:,137) = 4.2e-13_r8
      rate(:,139) = 8e-14_r8
      rate(:,140) = 3.9e-17_r8
      rate(:,147) = 6.9e-12_r8
      rate(:,148) = 7.2e-11_r8
      rate(:,149) = 1.6e-12_r8
      rate(:,155) = 1.8e-12_r8
      rate(:,159) = 1.8e-12_r8
      rate(:,171) = 3.5e-12_r8
      rate(:,173) = 1.3e-11_r8
      rate(:,174) = 2.2e-11_r8
      rate(:,175) = 5e-11_r8
      rate(:,210) = 1.7e-13_r8
      rate(:,212) = 2.607e-10_r8
      rate(:,213) = 9.75e-11_r8
      rate(:,214) = 2.07e-10_r8
      rate(:,215) = 2.088e-10_r8
      rate(:,216) = 1.17e-10_r8
      rate(:,217) = 4.644e-11_r8
      rate(:,218) = 1.204e-10_r8
      rate(:,219) = 9.9e-11_r8
      rate(:,220) = 3.3e-12_r8
      rate(:,239) = 4.5e-11_r8
      rate(:,240) = 4.62e-10_r8
      rate(:,241) = 1.2e-10_r8
      rate(:,242) = 9e-11_r8
      rate(:,243) = 3e-11_r8
      rate(:,248) = 2.14e-11_r8
      rate(:,249) = 1.9e-10_r8
      rate(:,262) = 2.57e-10_r8
      rate(:,263) = 1.8e-10_r8
      rate(:,264) = 1.794e-10_r8
      rate(:,265) = 1.3e-10_r8
      rate(:,266) = 7.65e-11_r8
      rate(:,279) = 4e-13_r8
      rate(:,283) = 1.31e-10_r8
      rate(:,284) = 3.5e-11_r8
      rate(:,285) = 9e-12_r8
      rate(:,292) = 6.8e-14_r8
      rate(:,293) = 2e-13_r8
      rate(:,308) = 1e-12_r8
      rate(:,312) = 1e-14_r8
      rate(:,313) = 1e-11_r8
      rate(:,314) = 1.15e-11_r8
      rate(:,315) = 4e-14_r8
      rate(:,328) = 3e-12_r8
      rate(:,329) = 6.7e-13_r8
      rate(:,339) = 3.5e-13_r8
      rate(:,340) = 5.4e-11_r8
      rate(:,343) = 2e-12_r8
      rate(:,344) = 1.4e-11_r8
      rate(:,347) = 2.4e-12_r8
      rate(:,358) = 5e-12_r8
      rate(:,368) = 2.2e-12_r8
      rate(:,370) = 6.7e-12_r8
      rate(:,373) = 3.5e-12_r8
      rate(:,376) = 1.3e-11_r8
      rate(:,377) = 1.4e-11_r8
      rate(:,381) = 2.4e-12_r8
      rate(:,382) = 1.4e-11_r8
      rate(:,387) = 2.4e-12_r8
      rate(:,388) = 4e-11_r8
      rate(:,389) = 4e-11_r8
      rate(:,391) = 1.4e-11_r8
      rate(:,395) = 2.4e-12_r8
      rate(:,396) = 4e-11_r8
      rate(:,400) = 7e-11_r8
      rate(:,401) = 1e-10_r8
      rate(:,406) = 2.4e-12_r8
      rate(:,421) = 4.7e-11_r8
      rate(:,434) = 2.1e-12_r8
      rate(:,435) = 2.8e-13_r8
      rate(:,443) = 1.7e-11_r8
      rate(:,449) = 8.4e-11_r8
      rate(:,451) = 1.9e-11_r8
      rate(:,452) = 1.2e-14_r8
      rate(:,453) = 2e-10_r8
      rate(:,460) = 2.4e-12_r8
      rate(:,461) = 2e-11_r8
      rate(:,465) = 2.3e-11_r8
      rate(:,466) = 2e-11_r8
      rate(:,470) = 3.3e-11_r8
      rate(:,471) = 1e-12_r8
      rate(:,472) = 5.7e-11_r8
      rate(:,473) = 3.4e-11_r8
      rate(:,478) = 2.3e-12_r8
      rate(:,480) = 1.2e-11_r8
      rate(:,481) = 5.7e-11_r8
      rate(:,482) = 2.8e-11_r8
      rate(:,483) = 6.6e-11_r8
      rate(:,484) = 1.4e-11_r8
      rate(:,487) = 1.9e-12_r8
      rate(:,499) = 6.34e-08_r8
      rate(:,505) = 1.9e-11_r8
      rate(:,508) = 1.2e-14_r8
      rate(:,509) = 2e-10_r8
      rate(:,520) = 1.34e-11_r8
      rate(:,526) = 1.34e-11_r8
      rate(:,531) = 1.7e-11_r8
      rate(:,551) = 2.31e-07_r8
      rate(:,552) = 2.31e-06_r8
      rate(:,553) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,128) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,129) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,130) = 2.64e-11_r8 * exp_fac(:)
      rate(:,131) = 6.6e-12_r8 * exp_fac(:)
      rate(:,136) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,138) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,141) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,142) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,145) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,146) = 1.4e-12_r8 * exp_fac(:)
      rate(:,397) = 1.05e-14_r8 * exp_fac(:)
      rate(:,516) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,151) = 3e-11_r8 * exp_fac(:)
      rate(:,237) = 5.5e-12_r8 * exp_fac(:)
      rate(:,276) = 3.8e-12_r8 * exp_fac(:)
      rate(:,297) = 3.8e-12_r8 * exp_fac(:)
      rate(:,324) = 3.8e-12_r8 * exp_fac(:)
      rate(:,332) = 3.8e-12_r8 * exp_fac(:)
      rate(:,336) = 3.8e-12_r8 * exp_fac(:)
      rate(:,352) = 2.3e-11_r8 * exp_fac(:)
      rate(:,362) = 3.8e-12_r8 * exp_fac(:)
      rate(:,372) = 3.8e-12_r8 * exp_fac(:)
      rate(:,399) = 1.52e-11_r8 * exp_fac(:)
      rate(:,407) = 1.52e-12_r8 * exp_fac(:)
      rate(:,413) = 3.8e-12_r8 * exp_fac(:)
      rate(:,416) = 3.8e-12_r8 * exp_fac(:)
      rate(:,420) = 3.8e-12_r8 * exp_fac(:)
      rate(:,436) = 3.8e-12_r8 * exp_fac(:)
      rate(:,440) = 3.8e-12_r8 * exp_fac(:)
      rate(:,446) = 3.8e-12_r8 * exp_fac(:)
      rate(:,450) = 3.8e-12_r8 * exp_fac(:)
      rate(:,152) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,153) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,154) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,156) = 4.8e-11_r8 * exp_fac(:)
      rate(:,235) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,157) = 1.8e-11_r8 * exp_fac(:)
      rate(:,310) = 4.2e-12_r8 * exp_fac(:)
      rate(:,331) = 4.2e-12_r8 * exp_fac(:)
      rate(:,360) = 4.2e-12_r8 * exp_fac(:)
      rate(:,380) = 4.4e-12_r8 * exp_fac(:)
      rate(:,386) = 4.4e-12_r8 * exp_fac(:)
      rate(:,459) = 4.2e-12_r8 * exp_fac(:)
      rate(:,464) = 4.2e-12_r8 * exp_fac(:)
      rate(:,469) = 4.2e-12_r8 * exp_fac(:)
      rate(:,158) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,162) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,163) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,164) = 2.9e-12_r8 * exp_fac(:)
      rate(:,165) = 1.45e-12_r8 * exp_fac(:)
      rate(:,166) = 1.45e-12_r8 * exp_fac(:)
      rate(:,167) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,168) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,169) = 1.2e-13_r8 * exp_fac(:)
      rate(:,195) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,172) = 1.7e-11_r8 * exp_fac(:)
      rate(:,270) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,176) = 3.44e-12_r8 * exp_fac(:)
      rate(:,228) = 2.3e-12_r8 * exp_fac(:)
      rate(:,231) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,177) = 3e-12_r8 * exp_fac(:)
      rate(:,236) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,179) = 7.26e-11_r8 * exp_fac(:)
      rate(:,180) = 4.64e-11_r8 * exp_fac(:)
      rate(:,187) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,188) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,189) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,190) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,191) = 1.4e-11_r8 * exp_fac(:)
      rate(:,205) = 7.4e-12_r8 * exp_fac(:)
      rate(:,306) = 8.1e-12_r8 * exp_fac(:)
      rate(:,192) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,193) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,194) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,196) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,197) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,198) = 2.6e-12_r8 * exp_fac(:)
      rate(:,199) = 6.4e-12_r8 * exp_fac(:)
      rate(:,229) = 4.1e-13_r8 * exp_fac(:)
      rate(:,409) = 7.5e-12_r8 * exp_fac(:)
      rate(:,423) = 7.5e-12_r8 * exp_fac(:)
      rate(:,426) = 7.5e-12_r8 * exp_fac(:)
      rate(:,429) = 7.5e-12_r8 * exp_fac(:)
      rate(:,200) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,202) = 3.6e-12_r8 * exp_fac(:)
      rate(:,251) = 2e-12_r8 * exp_fac(:)
      rate(:,203) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,204) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,206) = 6e-13_r8 * exp_fac(:)
      rate(:,226) = 1.5e-12_r8 * exp_fac(:)
      rate(:,234) = 1.9e-11_r8 * exp_fac(:)
      rate(:,207) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,208) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,209) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,211) = 3e-12_r8 * exp_fac(:)
      rate(:,245) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,223) = 1.7e-11_r8 * exp_fac(:)
      rate(:,250) = 6.3e-12_r8 * exp_fac(:)
      rate(:,224) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,225) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,227) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,230) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,233) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,238) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,244) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,246) = 1.4e-11_r8 * exp_fac(:)
      rate(:,248) = 2.14e-11_r8 * exp_fac(:)
      rate(:,249) = 1.9e-10_r8 * exp_fac(:)
      rate(:,262) = 2.57e-10_r8 * exp_fac(:)
      rate(:,263) = 1.8e-10_r8 * exp_fac(:)
      rate(:,264) = 1.794e-10_r8 * exp_fac(:)
      rate(:,265) = 1.3e-10_r8 * exp_fac(:)
      rate(:,266) = 7.65e-11_r8 * exp_fac(:)
      rate(:,279) = 4e-13_r8 * exp_fac(:)
      rate(:,283) = 1.31e-10_r8 * exp_fac(:)
      rate(:,284) = 3.5e-11_r8 * exp_fac(:)
      rate(:,285) = 9e-12_r8 * exp_fac(:)
      rate(:,292) = 6.8e-14_r8 * exp_fac(:)
      rate(:,293) = 2e-13_r8 * exp_fac(:)
      rate(:,308) = 1e-12_r8 * exp_fac(:)
      rate(:,312) = 1e-14_r8 * exp_fac(:)
      rate(:,313) = 1e-11_r8 * exp_fac(:)
      rate(:,314) = 1.15e-11_r8 * exp_fac(:)
      rate(:,315) = 4e-14_r8 * exp_fac(:)
      rate(:,328) = 3e-12_r8 * exp_fac(:)
      rate(:,329) = 6.7e-13_r8 * exp_fac(:)
      rate(:,339) = 3.5e-13_r8 * exp_fac(:)
      rate(:,340) = 5.4e-11_r8 * exp_fac(:)
      rate(:,343) = 2e-12_r8 * exp_fac(:)
      rate(:,344) = 1.4e-11_r8 * exp_fac(:)
      rate(:,347) = 2.4e-12_r8 * exp_fac(:)
      rate(:,358) = 5e-12_r8 * exp_fac(:)
      rate(:,368) = 2.2e-12_r8 * exp_fac(:)
      rate(:,370) = 6.7e-12_r8 * exp_fac(:)
      rate(:,373) = 3.5e-12_r8 * exp_fac(:)
      rate(:,376) = 1.3e-11_r8 * exp_fac(:)
      rate(:,377) = 1.4e-11_r8 * exp_fac(:)
      rate(:,381) = 2.4e-12_r8 * exp_fac(:)
      rate(:,382) = 1.4e-11_r8 * exp_fac(:)
      rate(:,387) = 2.4e-12_r8 * exp_fac(:)
      rate(:,388) = 4e-11_r8 * exp_fac(:)
      rate(:,389) = 4e-11_r8 * exp_fac(:)
      rate(:,391) = 1.4e-11_r8 * exp_fac(:)
      rate(:,395) = 2.4e-12_r8 * exp_fac(:)
      rate(:,396) = 4e-11_r8 * exp_fac(:)
      rate(:,400) = 7e-11_r8 * exp_fac(:)
      rate(:,401) = 1e-10_r8 * exp_fac(:)
      rate(:,406) = 2.4e-12_r8 * exp_fac(:)
      rate(:,421) = 4.7e-11_r8 * exp_fac(:)
      rate(:,434) = 2.1e-12_r8 * exp_fac(:)
      rate(:,435) = 2.8e-13_r8 * exp_fac(:)
      rate(:,443) = 1.7e-11_r8 * exp_fac(:)
      rate(:,449) = 8.4e-11_r8 * exp_fac(:)
      rate(:,451) = 1.9e-11_r8 * exp_fac(:)
      rate(:,452) = 1.2e-14_r8 * exp_fac(:)
      rate(:,453) = 2e-10_r8 * exp_fac(:)
      rate(:,460) = 2.4e-12_r8 * exp_fac(:)
      rate(:,461) = 2e-11_r8 * exp_fac(:)
      rate(:,465) = 2.3e-11_r8 * exp_fac(:)
      rate(:,466) = 2e-11_r8 * exp_fac(:)
      rate(:,470) = 3.3e-11_r8 * exp_fac(:)
      rate(:,471) = 1e-12_r8 * exp_fac(:)
      rate(:,472) = 5.7e-11_r8 * exp_fac(:)
      rate(:,473) = 3.4e-11_r8 * exp_fac(:)
      rate(:,478) = 2.3e-12_r8 * exp_fac(:)
      rate(:,480) = 1.2e-11_r8 * exp_fac(:)
      rate(:,481) = 5.7e-11_r8 * exp_fac(:)
      rate(:,482) = 2.8e-11_r8 * exp_fac(:)
      rate(:,483) = 6.6e-11_r8 * exp_fac(:)
      rate(:,484) = 1.4e-11_r8 * exp_fac(:)
      rate(:,487) = 1.9e-12_r8 * exp_fac(:)
      rate(:,499) = 6.34e-08_r8 * exp_fac(:)
      rate(:,505) = 1.9e-11_r8 * exp_fac(:)
      rate(:,508) = 1.2e-14_r8 * exp_fac(:)
      rate(:,509) = 2e-10_r8 * exp_fac(:)
      rate(:,520) = 1.34e-11_r8 * exp_fac(:)
      rate(:,526) = 1.34e-11_r8 * exp_fac(:)
      rate(:,531) = 1.7e-11_r8 * exp_fac(:)
      rate(:,551) = 2.31e-07_r8 * exp_fac(:)
      rate(:,552) = 2.31e-06_r8 * exp_fac(:)
      rate(:,553) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,247) = 6e-12_r8 * exp_fac(:)
      rate(:,345) = 5e-13_r8 * exp_fac(:)
      rate(:,378) = 5e-13_r8 * exp_fac(:)
      rate(:,383) = 5e-13_r8 * exp_fac(:)
      rate(:,392) = 5e-13_r8 * exp_fac(:)
      rate(:,403) = 5e-13_r8 * exp_fac(:)
      rate(:,252) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,253) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,254) = 1.64e-12_r8 * exp_fac(:)
      rate(:,364) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,255) = 2.03e-11_r8 * exp_fac(:)
      rate(:,486) = 3.4e-12_r8 * exp_fac(:)
      rate(:,256) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,257) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,258) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,259) = 1.25e-12_r8 * exp_fac(:)
      rate(:,269) = 3.4e-11_r8 * exp_fac(:)
      rate(:,260) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,261) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,267) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,268) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,271) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,272) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,273) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,274) = 2.8e-12_r8 * exp_fac(:)
      rate(:,335) = 2.9e-12_r8 * exp_fac(:)
      rate(:,275) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,277) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,280) = 7.5e-13_r8 * exp_fac(:)
      rate(:,294) = 7.5e-13_r8 * exp_fac(:)
      rate(:,309) = 7.5e-13_r8 * exp_fac(:)
      rate(:,330) = 7.5e-13_r8 * exp_fac(:)
      rate(:,334) = 8.6e-13_r8 * exp_fac(:)
      rate(:,346) = 8e-13_r8 * exp_fac(:)
      rate(:,359) = 7.5e-13_r8 * exp_fac(:)
      rate(:,369) = 7.5e-13_r8 * exp_fac(:)
      rate(:,379) = 8e-13_r8 * exp_fac(:)
      rate(:,384) = 8e-13_r8 * exp_fac(:)
      rate(:,393) = 8e-13_r8 * exp_fac(:)
      rate(:,404) = 8e-13_r8 * exp_fac(:)
      rate(:,411) = 7.5e-13_r8 * exp_fac(:)
      rate(:,415) = 7.5e-13_r8 * exp_fac(:)
      rate(:,418) = 7.5e-13_r8 * exp_fac(:)
      rate(:,431) = 7.5e-13_r8 * exp_fac(:)
      rate(:,438) = 7.5e-13_r8 * exp_fac(:)
      rate(:,444) = 7.5e-13_r8 * exp_fac(:)
      rate(:,447) = 7.5e-13_r8 * exp_fac(:)
      rate(:,458) = 7.5e-13_r8 * exp_fac(:)
      rate(:,463) = 7.5e-13_r8 * exp_fac(:)
      rate(:,468) = 7.5e-13_r8 * exp_fac(:)
      rate(:,511) = 7.5e-13_r8 * exp_fac(:)
      rate(:,518) = 7.5e-13_r8 * exp_fac(:)
      rate(:,528) = 7.5e-13_r8 * exp_fac(:)
      rate(:,532) = 7.5e-13_r8 * exp_fac(:)
      rate(:,281) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,282) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,286) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,291) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,295) = 2.6e-12_r8 * exp_fac(:)
      rate(:,412) = 2.6e-12_r8 * exp_fac(:)
      rate(:,417) = 2.6e-12_r8 * exp_fac(:)
      rate(:,419) = 2.6e-12_r8 * exp_fac(:)
      rate(:,432) = 2.6e-12_r8 * exp_fac(:)
      rate(:,439) = 2.6e-12_r8 * exp_fac(:)
      rate(:,445) = 2.6e-12_r8 * exp_fac(:)
      rate(:,448) = 2.6e-12_r8 * exp_fac(:)
      rate(:,512) = 2.6e-12_r8 * exp_fac(:)
      rate(:,519) = 2.6e-12_r8 * exp_fac(:)
      rate(:,529) = 2.6e-12_r8 * exp_fac(:)
      rate(:,533) = 2.6e-12_r8 * exp_fac(:)
      rate(:,296) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,298) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,299) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,300) = 1.4e-12_r8 * exp_fac(:)
      rate(:,320) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,301) = 4.63e-12_r8 * exp_fac(:)
      rate(:,515) = 2.7e-12_r8 * exp_fac(:)
      rate(:,302) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,303) = 2.9e-12_r8 * exp_fac(:)
      rate(:,304) = 2e-12_r8 * exp_fac(:)
      rate(:,333) = 7.1e-13_r8 * exp_fac(:)
      rate(:,354) = 2e-12_r8 * exp_fac(:)
      rate(:,457) = 2e-12_r8 * exp_fac(:)
      rate(:,462) = 2e-12_r8 * exp_fac(:)
      rate(:,467) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,305) = 4.3e-13_r8 * exp_fac(:)
      rate(:,355) = 4.3e-13_r8 * exp_fac(:)
      rate(:,408) = 4.3e-13_r8 * exp_fac(:)
      rate(:,422) = 4.3e-13_r8 * exp_fac(:)
      rate(:,425) = 4.3e-13_r8 * exp_fac(:)
      rate(:,428) = 4.3e-13_r8 * exp_fac(:)
      rate(:,307) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,311) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,319) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,321) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,322) = 1.41e-13_r8 * exp_fac(:)
      rate(:,506) = 2.75e-13_r8 * exp_fac(:)
      rate(:,514) = 2.12e-13_r8 * exp_fac(:)
      rate(:,522) = 2.6e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,323) = 2.7e-12_r8 * exp_fac(:)
      rate(:,348) = 2.7e-12_r8 * exp_fac(:)
      rate(:,349) = 1.3e-13_r8 * exp_fac(:)
      rate(:,351) = 9.6e-12_r8 * exp_fac(:)
      rate(:,357) = 5.3e-12_r8 * exp_fac(:)
      rate(:,394) = 2.7e-12_r8 * exp_fac(:)
      rate(:,405) = 2.7e-12_r8 * exp_fac(:)
      rate(:,507) = 2.7e-12_r8 * exp_fac(:)
      rate(:,523) = 2.7e-12_r8 * exp_fac(:)
      rate(:,325) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,326) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,327) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,341) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,342) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,350) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,353) = 4.6e-12_r8 * exp_fac(:)
      rate(:,356) = 2.3e-12_r8 * exp_fac(:)
      rate(:,361) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,365) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,371) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,374) = 1.86e-11_r8 * exp_fac(:)
      rate(:,375) = 1.86e-11_r8 * exp_fac(:)
      rate(:,385) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,390) = 3.03e-12_r8 * exp_fac(:)
      rate(:,513) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,398) = 2.54e-11_r8 * exp_fac(:)
      rate(:,517) = 2.54e-11_r8 * exp_fac(:)
      rate(:,402) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,410) = 2.3e-12_r8 * exp_fac(:)
      rate(:,510) = 2.3e-12_r8 * exp_fac(:)
      rate(:,414) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,433) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,441) = 1.7e-12_r8 * exp_fac(:)
      rate(:,527) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,454) = 1.2e-12_r8 * exp_fac(:)
      rate(:,521) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,455) = 6.3e-16_r8 * exp_fac(:)
      rate(:,524) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,456) = 1.2e-11_r8 * exp_fac(:)
      rate(:,525) = 1.2e-11_r8 * exp_fac(:)
      rate(:,474) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,475) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,476) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,477) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,485) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,488) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,491) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,150), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,160), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,178), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,181), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,182), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,183), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,201), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,221), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,232), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,278), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,288), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,289), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,290), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,316), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,317), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,337), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,363), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,366), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,424), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,427), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,430), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,437), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,479), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,134) = 1e-20_r8
      rate(:n,135) = 1.3e-16_r8
      rate(:n,139) = 8e-14_r8
      rate(:n,140) = 3.9e-17_r8
      rate(:n,147) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,129) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,130) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,131) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,136) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,138) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,141) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,142) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,151) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,152) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,153) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,156) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,157) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,158) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,163) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,167) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,168) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,176) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,177) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,150) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
