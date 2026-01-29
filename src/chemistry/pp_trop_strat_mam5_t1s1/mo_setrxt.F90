
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
      rate(:,135) = 6.9e-12_r8
      rate(:,136) = 7.2e-11_r8
      rate(:,137) = 1.6e-12_r8
      rate(:,143) = 1.8e-12_r8
      rate(:,147) = 1.8e-12_r8
      rate(:,159) = 3.5e-12_r8
      rate(:,161) = 1.3e-11_r8
      rate(:,162) = 2.2e-11_r8
      rate(:,163) = 5e-11_r8
      rate(:,198) = 1.7e-13_r8
      rate(:,200) = 2.607e-10_r8
      rate(:,201) = 9.75e-11_r8
      rate(:,202) = 2.07e-10_r8
      rate(:,203) = 2.088e-10_r8
      rate(:,204) = 1.17e-10_r8
      rate(:,205) = 4.644e-11_r8
      rate(:,206) = 1.204e-10_r8
      rate(:,207) = 9.9e-11_r8
      rate(:,208) = 3.3e-12_r8
      rate(:,227) = 4.5e-11_r8
      rate(:,228) = 4.62e-10_r8
      rate(:,229) = 1.2e-10_r8
      rate(:,230) = 9e-11_r8
      rate(:,231) = 3e-11_r8
      rate(:,236) = 2.14e-11_r8
      rate(:,237) = 1.9e-10_r8
      rate(:,250) = 2.57e-10_r8
      rate(:,251) = 1.8e-10_r8
      rate(:,252) = 1.794e-10_r8
      rate(:,253) = 1.3e-10_r8
      rate(:,254) = 7.65e-11_r8
      rate(:,267) = 4e-13_r8
      rate(:,271) = 1.31e-10_r8
      rate(:,272) = 3.5e-11_r8
      rate(:,273) = 9e-12_r8
      rate(:,280) = 6.8e-14_r8
      rate(:,281) = 2e-13_r8
      rate(:,296) = 1e-12_r8
      rate(:,300) = 1e-14_r8
      rate(:,301) = 1e-11_r8
      rate(:,302) = 1.15e-11_r8
      rate(:,303) = 4e-14_r8
      rate(:,316) = 3e-12_r8
      rate(:,317) = 6.7e-13_r8
      rate(:,327) = 3.5e-13_r8
      rate(:,328) = 5.4e-11_r8
      rate(:,331) = 2e-12_r8
      rate(:,332) = 1.4e-11_r8
      rate(:,335) = 2.4e-12_r8
      rate(:,346) = 5e-12_r8
      rate(:,356) = 1.6e-12_r8
      rate(:,358) = 6.7e-12_r8
      rate(:,361) = 3.5e-12_r8
      rate(:,364) = 1.3e-11_r8
      rate(:,365) = 1.4e-11_r8
      rate(:,369) = 2.4e-12_r8
      rate(:,370) = 1.4e-11_r8
      rate(:,375) = 2.4e-12_r8
      rate(:,376) = 4e-11_r8
      rate(:,377) = 4e-11_r8
      rate(:,379) = 1.4e-11_r8
      rate(:,383) = 2.4e-12_r8
      rate(:,384) = 4e-11_r8
      rate(:,388) = 7e-11_r8
      rate(:,389) = 1e-10_r8
      rate(:,394) = 2.4e-12_r8
      rate(:,409) = 4.7e-11_r8
      rate(:,422) = 2.1e-12_r8
      rate(:,423) = 2.8e-13_r8
      rate(:,431) = 1.7e-11_r8
      rate(:,437) = 8.4e-11_r8
      rate(:,439) = 1.9e-11_r8
      rate(:,440) = 1.2e-14_r8
      rate(:,441) = 2e-10_r8
      rate(:,448) = 2.4e-12_r8
      rate(:,449) = 2e-11_r8
      rate(:,453) = 2.3e-11_r8
      rate(:,454) = 2e-11_r8
      rate(:,458) = 3.3e-11_r8
      rate(:,459) = 1e-12_r8
      rate(:,460) = 5.7e-11_r8
      rate(:,461) = 3.4e-11_r8
      rate(:,466) = 2.3e-12_r8
      rate(:,468) = 1.2e-11_r8
      rate(:,469) = 5.7e-11_r8
      rate(:,470) = 2.8e-11_r8
      rate(:,471) = 6.6e-11_r8
      rate(:,472) = 1.4e-11_r8
      rate(:,475) = 1.9e-12_r8
      rate(:,488) = 6.34e-08_r8
      rate(:,494) = 1.9e-11_r8
      rate(:,497) = 1.2e-14_r8
      rate(:,498) = 2e-10_r8
      rate(:,509) = 1.34e-11_r8
      rate(:,515) = 1.34e-11_r8
      rate(:,519) = 1.7e-11_r8
      rate(:,539) = 2.31e-07_r8
      rate(:,540) = 2.31e-06_r8
      rate(:,541) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,126) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,128) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,130) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,133) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,134) = 1.4e-12_r8 * exp_fac(:)
      rate(:,385) = 1.05e-14_r8 * exp_fac(:)
      rate(:,505) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,139) = 3e-11_r8 * exp_fac(:)
      rate(:,225) = 5.5e-12_r8 * exp_fac(:)
      rate(:,264) = 3.8e-12_r8 * exp_fac(:)
      rate(:,285) = 3.8e-12_r8 * exp_fac(:)
      rate(:,312) = 3.8e-12_r8 * exp_fac(:)
      rate(:,320) = 3.8e-12_r8 * exp_fac(:)
      rate(:,324) = 3.8e-12_r8 * exp_fac(:)
      rate(:,340) = 2.3e-11_r8 * exp_fac(:)
      rate(:,350) = 3.8e-12_r8 * exp_fac(:)
      rate(:,360) = 3.8e-12_r8 * exp_fac(:)
      rate(:,387) = 1.52e-11_r8 * exp_fac(:)
      rate(:,395) = 1.52e-12_r8 * exp_fac(:)
      rate(:,401) = 3.8e-12_r8 * exp_fac(:)
      rate(:,404) = 3.8e-12_r8 * exp_fac(:)
      rate(:,408) = 3.8e-12_r8 * exp_fac(:)
      rate(:,424) = 3.8e-12_r8 * exp_fac(:)
      rate(:,428) = 3.8e-12_r8 * exp_fac(:)
      rate(:,434) = 3.8e-12_r8 * exp_fac(:)
      rate(:,438) = 3.8e-12_r8 * exp_fac(:)
      rate(:,140) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,141) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,142) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,144) = 4.8e-11_r8 * exp_fac(:)
      rate(:,223) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,145) = 1.8e-11_r8 * exp_fac(:)
      rate(:,298) = 4.2e-12_r8 * exp_fac(:)
      rate(:,311) = 4.2e-12_r8 * exp_fac(:)
      rate(:,319) = 4.2e-12_r8 * exp_fac(:)
      rate(:,348) = 4.2e-12_r8 * exp_fac(:)
      rate(:,368) = 4.4e-12_r8 * exp_fac(:)
      rate(:,374) = 4.4e-12_r8 * exp_fac(:)
      rate(:,447) = 4.2e-12_r8 * exp_fac(:)
      rate(:,452) = 4.2e-12_r8 * exp_fac(:)
      rate(:,457) = 4.2e-12_r8 * exp_fac(:)
      rate(:,146) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,150) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,151) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,152) = 2.9e-12_r8 * exp_fac(:)
      rate(:,153) = 1.45e-12_r8 * exp_fac(:)
      rate(:,154) = 1.45e-12_r8 * exp_fac(:)
      rate(:,155) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,156) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,157) = 1.2e-13_r8 * exp_fac(:)
      rate(:,183) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,160) = 1.7e-11_r8 * exp_fac(:)
      rate(:,258) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,164) = 3.44e-12_r8 * exp_fac(:)
      rate(:,216) = 2.3e-12_r8 * exp_fac(:)
      rate(:,219) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,165) = 3e-12_r8 * exp_fac(:)
      rate(:,224) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,167) = 7.26e-11_r8 * exp_fac(:)
      rate(:,168) = 4.64e-11_r8 * exp_fac(:)
      rate(:,175) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,176) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,177) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,178) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,179) = 1.4e-11_r8 * exp_fac(:)
      rate(:,193) = 7.4e-12_r8 * exp_fac(:)
      rate(:,294) = 8.1e-12_r8 * exp_fac(:)
      rate(:,180) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,181) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,182) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,184) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,185) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,186) = 2.6e-12_r8 * exp_fac(:)
      rate(:,187) = 6.4e-12_r8 * exp_fac(:)
      rate(:,217) = 4.1e-13_r8 * exp_fac(:)
      rate(:,397) = 7.5e-12_r8 * exp_fac(:)
      rate(:,411) = 7.5e-12_r8 * exp_fac(:)
      rate(:,414) = 7.5e-12_r8 * exp_fac(:)
      rate(:,417) = 7.5e-12_r8 * exp_fac(:)
      rate(:,188) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,190) = 3.6e-12_r8 * exp_fac(:)
      rate(:,239) = 2e-12_r8 * exp_fac(:)
      rate(:,191) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,192) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,194) = 6e-13_r8 * exp_fac(:)
      rate(:,214) = 1.5e-12_r8 * exp_fac(:)
      rate(:,222) = 1.9e-11_r8 * exp_fac(:)
      rate(:,195) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,196) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,197) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,199) = 3e-12_r8 * exp_fac(:)
      rate(:,233) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,211) = 1.7e-11_r8 * exp_fac(:)
      rate(:,238) = 6.3e-12_r8 * exp_fac(:)
      rate(:,212) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,213) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,215) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,218) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,221) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,226) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,232) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,234) = 1.4e-11_r8 * exp_fac(:)
      rate(:,236) = 2.14e-11_r8 * exp_fac(:)
      rate(:,237) = 1.9e-10_r8 * exp_fac(:)
      rate(:,250) = 2.57e-10_r8 * exp_fac(:)
      rate(:,251) = 1.8e-10_r8 * exp_fac(:)
      rate(:,252) = 1.794e-10_r8 * exp_fac(:)
      rate(:,253) = 1.3e-10_r8 * exp_fac(:)
      rate(:,254) = 7.65e-11_r8 * exp_fac(:)
      rate(:,267) = 4e-13_r8 * exp_fac(:)
      rate(:,271) = 1.31e-10_r8 * exp_fac(:)
      rate(:,272) = 3.5e-11_r8 * exp_fac(:)
      rate(:,273) = 9e-12_r8 * exp_fac(:)
      rate(:,280) = 6.8e-14_r8 * exp_fac(:)
      rate(:,281) = 2e-13_r8 * exp_fac(:)
      rate(:,296) = 1e-12_r8 * exp_fac(:)
      rate(:,300) = 1e-14_r8 * exp_fac(:)
      rate(:,301) = 1e-11_r8 * exp_fac(:)
      rate(:,302) = 1.15e-11_r8 * exp_fac(:)
      rate(:,303) = 4e-14_r8 * exp_fac(:)
      rate(:,316) = 3e-12_r8 * exp_fac(:)
      rate(:,317) = 6.7e-13_r8 * exp_fac(:)
      rate(:,327) = 3.5e-13_r8 * exp_fac(:)
      rate(:,328) = 5.4e-11_r8 * exp_fac(:)
      rate(:,331) = 2e-12_r8 * exp_fac(:)
      rate(:,332) = 1.4e-11_r8 * exp_fac(:)
      rate(:,335) = 2.4e-12_r8 * exp_fac(:)
      rate(:,346) = 5e-12_r8 * exp_fac(:)
      rate(:,356) = 1.6e-12_r8 * exp_fac(:)
      rate(:,358) = 6.7e-12_r8 * exp_fac(:)
      rate(:,361) = 3.5e-12_r8 * exp_fac(:)
      rate(:,364) = 1.3e-11_r8 * exp_fac(:)
      rate(:,365) = 1.4e-11_r8 * exp_fac(:)
      rate(:,369) = 2.4e-12_r8 * exp_fac(:)
      rate(:,370) = 1.4e-11_r8 * exp_fac(:)
      rate(:,375) = 2.4e-12_r8 * exp_fac(:)
      rate(:,376) = 4e-11_r8 * exp_fac(:)
      rate(:,377) = 4e-11_r8 * exp_fac(:)
      rate(:,379) = 1.4e-11_r8 * exp_fac(:)
      rate(:,383) = 2.4e-12_r8 * exp_fac(:)
      rate(:,384) = 4e-11_r8 * exp_fac(:)
      rate(:,388) = 7e-11_r8 * exp_fac(:)
      rate(:,389) = 1e-10_r8 * exp_fac(:)
      rate(:,394) = 2.4e-12_r8 * exp_fac(:)
      rate(:,409) = 4.7e-11_r8 * exp_fac(:)
      rate(:,422) = 2.1e-12_r8 * exp_fac(:)
      rate(:,423) = 2.8e-13_r8 * exp_fac(:)
      rate(:,431) = 1.7e-11_r8 * exp_fac(:)
      rate(:,437) = 8.4e-11_r8 * exp_fac(:)
      rate(:,439) = 1.9e-11_r8 * exp_fac(:)
      rate(:,440) = 1.2e-14_r8 * exp_fac(:)
      rate(:,441) = 2e-10_r8 * exp_fac(:)
      rate(:,448) = 2.4e-12_r8 * exp_fac(:)
      rate(:,449) = 2e-11_r8 * exp_fac(:)
      rate(:,453) = 2.3e-11_r8 * exp_fac(:)
      rate(:,454) = 2e-11_r8 * exp_fac(:)
      rate(:,458) = 3.3e-11_r8 * exp_fac(:)
      rate(:,459) = 1e-12_r8 * exp_fac(:)
      rate(:,460) = 5.7e-11_r8 * exp_fac(:)
      rate(:,461) = 3.4e-11_r8 * exp_fac(:)
      rate(:,466) = 2.3e-12_r8 * exp_fac(:)
      rate(:,468) = 1.2e-11_r8 * exp_fac(:)
      rate(:,469) = 5.7e-11_r8 * exp_fac(:)
      rate(:,470) = 2.8e-11_r8 * exp_fac(:)
      rate(:,471) = 6.6e-11_r8 * exp_fac(:)
      rate(:,472) = 1.4e-11_r8 * exp_fac(:)
      rate(:,475) = 1.9e-12_r8 * exp_fac(:)
      rate(:,488) = 6.34e-08_r8 * exp_fac(:)
      rate(:,494) = 1.9e-11_r8 * exp_fac(:)
      rate(:,497) = 1.2e-14_r8 * exp_fac(:)
      rate(:,498) = 2e-10_r8 * exp_fac(:)
      rate(:,509) = 1.34e-11_r8 * exp_fac(:)
      rate(:,515) = 1.34e-11_r8 * exp_fac(:)
      rate(:,519) = 1.7e-11_r8 * exp_fac(:)
      rate(:,539) = 2.31e-07_r8 * exp_fac(:)
      rate(:,540) = 2.31e-06_r8 * exp_fac(:)
      rate(:,541) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,235) = 6e-12_r8 * exp_fac(:)
      rate(:,333) = 5e-13_r8 * exp_fac(:)
      rate(:,366) = 5e-13_r8 * exp_fac(:)
      rate(:,371) = 5e-13_r8 * exp_fac(:)
      rate(:,380) = 5e-13_r8 * exp_fac(:)
      rate(:,391) = 5e-13_r8 * exp_fac(:)
      rate(:,240) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,241) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,242) = 1.64e-12_r8 * exp_fac(:)
      rate(:,352) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,243) = 2.03e-11_r8 * exp_fac(:)
      rate(:,474) = 3.4e-12_r8 * exp_fac(:)
      rate(:,244) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,245) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,246) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,247) = 1.25e-12_r8 * exp_fac(:)
      rate(:,257) = 3.4e-11_r8 * exp_fac(:)
      rate(:,248) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,249) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,255) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,256) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,259) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,260) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,261) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,262) = 2.8e-12_r8 * exp_fac(:)
      rate(:,323) = 2.9e-12_r8 * exp_fac(:)
      rate(:,263) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,265) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,268) = 7.5e-13_r8 * exp_fac(:)
      rate(:,282) = 7.5e-13_r8 * exp_fac(:)
      rate(:,297) = 7.5e-13_r8 * exp_fac(:)
      rate(:,310) = 7.5e-13_r8 * exp_fac(:)
      rate(:,318) = 7.5e-13_r8 * exp_fac(:)
      rate(:,322) = 8.6e-13_r8 * exp_fac(:)
      rate(:,334) = 8e-13_r8 * exp_fac(:)
      rate(:,347) = 7.5e-13_r8 * exp_fac(:)
      rate(:,357) = 7.5e-13_r8 * exp_fac(:)
      rate(:,367) = 8e-13_r8 * exp_fac(:)
      rate(:,372) = 8e-13_r8 * exp_fac(:)
      rate(:,381) = 8e-13_r8 * exp_fac(:)
      rate(:,392) = 8e-13_r8 * exp_fac(:)
      rate(:,399) = 7.5e-13_r8 * exp_fac(:)
      rate(:,403) = 7.5e-13_r8 * exp_fac(:)
      rate(:,406) = 7.5e-13_r8 * exp_fac(:)
      rate(:,419) = 7.5e-13_r8 * exp_fac(:)
      rate(:,426) = 7.5e-13_r8 * exp_fac(:)
      rate(:,432) = 7.5e-13_r8 * exp_fac(:)
      rate(:,435) = 7.5e-13_r8 * exp_fac(:)
      rate(:,446) = 7.5e-13_r8 * exp_fac(:)
      rate(:,451) = 7.5e-13_r8 * exp_fac(:)
      rate(:,456) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 7.5e-13_r8 * exp_fac(:)
      rate(:,507) = 7.5e-13_r8 * exp_fac(:)
      rate(:,517) = 7.5e-13_r8 * exp_fac(:)
      rate(:,520) = 7.5e-13_r8 * exp_fac(:)
      rate(:,269) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,270) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,274) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,279) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,283) = 2.6e-12_r8 * exp_fac(:)
      rate(:,400) = 2.6e-12_r8 * exp_fac(:)
      rate(:,405) = 2.6e-12_r8 * exp_fac(:)
      rate(:,407) = 2.6e-12_r8 * exp_fac(:)
      rate(:,420) = 2.6e-12_r8 * exp_fac(:)
      rate(:,427) = 2.6e-12_r8 * exp_fac(:)
      rate(:,433) = 2.6e-12_r8 * exp_fac(:)
      rate(:,436) = 2.6e-12_r8 * exp_fac(:)
      rate(:,501) = 2.6e-12_r8 * exp_fac(:)
      rate(:,508) = 2.6e-12_r8 * exp_fac(:)
      rate(:,518) = 2.6e-12_r8 * exp_fac(:)
      rate(:,521) = 2.6e-12_r8 * exp_fac(:)
      rate(:,284) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,286) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,287) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,288) = 1.4e-12_r8 * exp_fac(:)
      rate(:,308) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,289) = 4.63e-12_r8 * exp_fac(:)
      rate(:,504) = 2.7e-12_r8 * exp_fac(:)
      rate(:,290) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,291) = 2.9e-12_r8 * exp_fac(:)
      rate(:,292) = 2e-12_r8 * exp_fac(:)
      rate(:,321) = 7.1e-13_r8 * exp_fac(:)
      rate(:,342) = 2e-12_r8 * exp_fac(:)
      rate(:,445) = 2e-12_r8 * exp_fac(:)
      rate(:,450) = 2e-12_r8 * exp_fac(:)
      rate(:,455) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,293) = 4.3e-13_r8 * exp_fac(:)
      rate(:,343) = 4.3e-13_r8 * exp_fac(:)
      rate(:,396) = 4.3e-13_r8 * exp_fac(:)
      rate(:,410) = 4.3e-13_r8 * exp_fac(:)
      rate(:,413) = 4.3e-13_r8 * exp_fac(:)
      rate(:,416) = 4.3e-13_r8 * exp_fac(:)
      rate(:,295) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,299) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,307) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,309) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,313) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,314) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,315) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,329) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,330) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,336) = 2.7e-12_r8 * exp_fac(:)
      rate(:,337) = 1.3e-13_r8 * exp_fac(:)
      rate(:,339) = 9.6e-12_r8 * exp_fac(:)
      rate(:,345) = 5.3e-12_r8 * exp_fac(:)
      rate(:,382) = 2.7e-12_r8 * exp_fac(:)
      rate(:,393) = 2.7e-12_r8 * exp_fac(:)
      rate(:,496) = 2.7e-12_r8 * exp_fac(:)
      rate(:,512) = 2.7e-12_r8 * exp_fac(:)
      rate(:,338) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,341) = 4.6e-12_r8 * exp_fac(:)
      rate(:,344) = 2.3e-12_r8 * exp_fac(:)
      rate(:,349) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,353) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,359) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,362) = 1.86e-11_r8 * exp_fac(:)
      rate(:,363) = 1.86e-11_r8 * exp_fac(:)
      rate(:,373) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,378) = 3.03e-12_r8 * exp_fac(:)
      rate(:,502) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,386) = 2.54e-11_r8 * exp_fac(:)
      rate(:,506) = 2.54e-11_r8 * exp_fac(:)
      rate(:,390) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,398) = 2.3e-12_r8 * exp_fac(:)
      rate(:,499) = 2.3e-12_r8 * exp_fac(:)
      rate(:,402) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,421) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,429) = 1.7e-12_r8 * exp_fac(:)
      rate(:,516) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,442) = 1.2e-12_r8 * exp_fac(:)
      rate(:,510) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,443) = 6.3e-16_r8 * exp_fac(:)
      rate(:,513) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,444) = 1.2e-11_r8 * exp_fac(:)
      rate(:,514) = 1.2e-11_r8 * exp_fac(:)
      rate(:,462) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,463) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,464) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,465) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,473) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,476) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,479) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,495) = 2.75e-13_r8 * exp_fac(:)
      rate(:,503) = 2.12e-13_r8 * exp_fac(:)
      rate(:,511) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,138), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,148), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,158), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,166), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,169), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,170), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,171), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,189), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,220), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,266), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,276), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,277), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,278), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,304), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,305), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,325), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,351), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,354), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,412), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,415), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,418), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,425), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,467), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,135) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,127) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,130) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,139) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,140) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,141) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,144) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,145) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,146) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,151) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,155) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,156) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,164) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,165) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,138) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
