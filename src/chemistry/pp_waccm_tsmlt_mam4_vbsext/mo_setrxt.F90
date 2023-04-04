
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

      rate(:,171) = 0.000258_r8
      rate(:,172) = 0.085_r8
      rate(:,173) = 1.2e-10_r8
      rate(:,178) = 1.2e-10_r8
      rate(:,179) = 1e-20_r8
      rate(:,180) = 1.3e-16_r8
      rate(:,182) = 4.2e-13_r8
      rate(:,184) = 8e-14_r8
      rate(:,185) = 3.9e-17_r8
      rate(:,192) = 6.9e-12_r8
      rate(:,193) = 7.2e-11_r8
      rate(:,194) = 1.6e-12_r8
      rate(:,200) = 1.8e-12_r8
      rate(:,204) = 1.8e-12_r8
      rate(:,208) = 7e-13_r8
      rate(:,209) = 5e-12_r8
      rate(:,218) = 3.5e-12_r8
      rate(:,220) = 1.3e-11_r8
      rate(:,221) = 2.2e-11_r8
      rate(:,222) = 5e-11_r8
      rate(:,257) = 1.7e-13_r8
      rate(:,259) = 2.607e-10_r8
      rate(:,260) = 9.75e-11_r8
      rate(:,261) = 2.07e-10_r8
      rate(:,262) = 2.088e-10_r8
      rate(:,263) = 1.17e-10_r8
      rate(:,264) = 4.644e-11_r8
      rate(:,265) = 1.204e-10_r8
      rate(:,266) = 9.9e-11_r8
      rate(:,267) = 3.3e-12_r8
      rate(:,286) = 4.5e-11_r8
      rate(:,287) = 4.62e-10_r8
      rate(:,288) = 1.2e-10_r8
      rate(:,289) = 9e-11_r8
      rate(:,290) = 3e-11_r8
      rate(:,295) = 2.14e-11_r8
      rate(:,296) = 1.9e-10_r8
      rate(:,309) = 2.57e-10_r8
      rate(:,310) = 1.8e-10_r8
      rate(:,311) = 1.794e-10_r8
      rate(:,312) = 1.3e-10_r8
      rate(:,313) = 7.65e-11_r8
      rate(:,326) = 4e-13_r8
      rate(:,330) = 1.31e-10_r8
      rate(:,331) = 3.5e-11_r8
      rate(:,332) = 9e-12_r8
      rate(:,339) = 6.8e-14_r8
      rate(:,340) = 2e-13_r8
      rate(:,355) = 1e-12_r8
      rate(:,359) = 1e-14_r8
      rate(:,360) = 1e-11_r8
      rate(:,361) = 1.15e-11_r8
      rate(:,362) = 4e-14_r8
      rate(:,375) = 1.45e-10_r8
      rate(:,376) = 3e-12_r8
      rate(:,377) = 6.7e-13_r8
      rate(:,387) = 3.5e-13_r8
      rate(:,388) = 5.4e-11_r8
      rate(:,391) = 2e-12_r8
      rate(:,392) = 1.4e-11_r8
      rate(:,395) = 2.4e-12_r8
      rate(:,406) = 5e-12_r8
      rate(:,416) = 1.6e-12_r8
      rate(:,418) = 6.7e-12_r8
      rate(:,421) = 3.5e-12_r8
      rate(:,424) = 1.3e-11_r8
      rate(:,425) = 1.4e-11_r8
      rate(:,429) = 2.4e-12_r8
      rate(:,430) = 1.4e-11_r8
      rate(:,435) = 2.4e-12_r8
      rate(:,436) = 4e-11_r8
      rate(:,437) = 4e-11_r8
      rate(:,439) = 1.4e-11_r8
      rate(:,443) = 2.4e-12_r8
      rate(:,444) = 4e-11_r8
      rate(:,448) = 7e-11_r8
      rate(:,449) = 1e-10_r8
      rate(:,454) = 2.4e-12_r8
      rate(:,469) = 4.7e-11_r8
      rate(:,482) = 2.1e-12_r8
      rate(:,483) = 2.8e-13_r8
      rate(:,491) = 1.7e-11_r8
      rate(:,497) = 8.4e-11_r8
      rate(:,499) = 1.9e-11_r8
      rate(:,500) = 1.2e-14_r8
      rate(:,501) = 2e-10_r8
      rate(:,508) = 2.4e-12_r8
      rate(:,509) = 2e-11_r8
      rate(:,513) = 2.3e-11_r8
      rate(:,514) = 2e-11_r8
      rate(:,518) = 3.3e-11_r8
      rate(:,519) = 1e-12_r8
      rate(:,520) = 5.7e-11_r8
      rate(:,521) = 3.4e-11_r8
      rate(:,526) = 2.3e-12_r8
      rate(:,528) = 1.2e-11_r8
      rate(:,529) = 5.7e-11_r8
      rate(:,530) = 2.8e-11_r8
      rate(:,531) = 6.6e-11_r8
      rate(:,532) = 1.4e-11_r8
      rate(:,535) = 1.9e-12_r8
      rate(:,547) = 6.34e-08_r8
      rate(:,553) = 1.9e-11_r8
      rate(:,556) = 1.2e-14_r8
      rate(:,557) = 2e-10_r8
      rate(:,568) = 1.34e-11_r8
      rate(:,571) = 1.34e-11_r8
      rate(:,577) = 1.34e-11_r8
      rate(:,578) = 1.34e-11_r8
      rate(:,583) = 1.7e-11_r8
      rate(:,606) = 6e-11_r8
      rate(:,609) = 1e-12_r8
      rate(:,610) = 4e-10_r8
      rate(:,611) = 2e-10_r8
      rate(:,612) = 1e-10_r8
      rate(:,613) = 5e-16_r8
      rate(:,614) = 4.4e-10_r8
      rate(:,615) = 9e-10_r8
      rate(:,618) = 1.29e-07_r8
      rate(:,619) = 2.31e-07_r8
      rate(:,620) = 2.31e-06_r8
      rate(:,621) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,174) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,175) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,176) = 2.64e-11_r8 * exp_fac(:)
      rate(:,177) = 6.6e-12_r8 * exp_fac(:)
      rate(:,181) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,183) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,186) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,187) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,190) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,191) = 1.4e-12_r8 * exp_fac(:)
      rate(:,445) = 1.05e-14_r8 * exp_fac(:)
      rate(:,564) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,196) = 3e-11_r8 * exp_fac(:)
      rate(:,284) = 5.5e-12_r8 * exp_fac(:)
      rate(:,323) = 3.8e-12_r8 * exp_fac(:)
      rate(:,344) = 3.8e-12_r8 * exp_fac(:)
      rate(:,371) = 3.8e-12_r8 * exp_fac(:)
      rate(:,380) = 3.8e-12_r8 * exp_fac(:)
      rate(:,384) = 3.8e-12_r8 * exp_fac(:)
      rate(:,400) = 2.3e-11_r8 * exp_fac(:)
      rate(:,410) = 3.8e-12_r8 * exp_fac(:)
      rate(:,420) = 3.8e-12_r8 * exp_fac(:)
      rate(:,447) = 1.52e-11_r8 * exp_fac(:)
      rate(:,455) = 1.52e-12_r8 * exp_fac(:)
      rate(:,461) = 3.8e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,468) = 3.8e-12_r8 * exp_fac(:)
      rate(:,484) = 3.8e-12_r8 * exp_fac(:)
      rate(:,488) = 3.8e-12_r8 * exp_fac(:)
      rate(:,494) = 3.8e-12_r8 * exp_fac(:)
      rate(:,498) = 3.8e-12_r8 * exp_fac(:)
      rate(:,197) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,198) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,199) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,201) = 4.8e-11_r8 * exp_fac(:)
      rate(:,282) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,202) = 1.8e-11_r8 * exp_fac(:)
      rate(:,357) = 4.2e-12_r8 * exp_fac(:)
      rate(:,370) = 4.2e-12_r8 * exp_fac(:)
      rate(:,379) = 4.2e-12_r8 * exp_fac(:)
      rate(:,408) = 4.2e-12_r8 * exp_fac(:)
      rate(:,428) = 4.4e-12_r8 * exp_fac(:)
      rate(:,434) = 4.4e-12_r8 * exp_fac(:)
      rate(:,507) = 4.2e-12_r8 * exp_fac(:)
      rate(:,512) = 4.2e-12_r8 * exp_fac(:)
      rate(:,517) = 4.2e-12_r8 * exp_fac(:)
      rate(:,203) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,207) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,210) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,211) = 2.9e-12_r8 * exp_fac(:)
      rate(:,212) = 1.45e-12_r8 * exp_fac(:)
      rate(:,213) = 1.45e-12_r8 * exp_fac(:)
      rate(:,214) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,215) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,216) = 1.2e-13_r8 * exp_fac(:)
      rate(:,242) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,219) = 1.7e-11_r8 * exp_fac(:)
      rate(:,317) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,223) = 3.44e-12_r8 * exp_fac(:)
      rate(:,275) = 2.3e-12_r8 * exp_fac(:)
      rate(:,278) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,224) = 3e-12_r8 * exp_fac(:)
      rate(:,283) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,226) = 7.26e-11_r8 * exp_fac(:)
      rate(:,227) = 4.64e-11_r8 * exp_fac(:)
      rate(:,234) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,235) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,236) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,237) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,238) = 1.4e-11_r8 * exp_fac(:)
      rate(:,252) = 7.4e-12_r8 * exp_fac(:)
      rate(:,353) = 8.1e-12_r8 * exp_fac(:)
      rate(:,239) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,240) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,241) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,243) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,244) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,245) = 2.6e-12_r8 * exp_fac(:)
      rate(:,246) = 6.4e-12_r8 * exp_fac(:)
      rate(:,276) = 4.1e-13_r8 * exp_fac(:)
      rate(:,457) = 7.5e-12_r8 * exp_fac(:)
      rate(:,471) = 7.5e-12_r8 * exp_fac(:)
      rate(:,474) = 7.5e-12_r8 * exp_fac(:)
      rate(:,477) = 7.5e-12_r8 * exp_fac(:)
      rate(:,247) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,249) = 3.6e-12_r8 * exp_fac(:)
      rate(:,298) = 2e-12_r8 * exp_fac(:)
      rate(:,250) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,251) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,253) = 6e-13_r8 * exp_fac(:)
      rate(:,273) = 1.5e-12_r8 * exp_fac(:)
      rate(:,281) = 1.9e-11_r8 * exp_fac(:)
      rate(:,254) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,255) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,256) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,258) = 3e-12_r8 * exp_fac(:)
      rate(:,292) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,270) = 1.7e-11_r8 * exp_fac(:)
      rate(:,297) = 6.3e-12_r8 * exp_fac(:)
      rate(:,271) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,272) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,274) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,277) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,280) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,285) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,291) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,293) = 1.4e-11_r8 * exp_fac(:)
      rate(:,295) = 2.14e-11_r8 * exp_fac(:)
      rate(:,296) = 1.9e-10_r8 * exp_fac(:)
      rate(:,309) = 2.57e-10_r8 * exp_fac(:)
      rate(:,310) = 1.8e-10_r8 * exp_fac(:)
      rate(:,311) = 1.794e-10_r8 * exp_fac(:)
      rate(:,312) = 1.3e-10_r8 * exp_fac(:)
      rate(:,313) = 7.65e-11_r8 * exp_fac(:)
      rate(:,326) = 4e-13_r8 * exp_fac(:)
      rate(:,330) = 1.31e-10_r8 * exp_fac(:)
      rate(:,331) = 3.5e-11_r8 * exp_fac(:)
      rate(:,332) = 9e-12_r8 * exp_fac(:)
      rate(:,339) = 6.8e-14_r8 * exp_fac(:)
      rate(:,340) = 2e-13_r8 * exp_fac(:)
      rate(:,355) = 1e-12_r8 * exp_fac(:)
      rate(:,359) = 1e-14_r8 * exp_fac(:)
      rate(:,360) = 1e-11_r8 * exp_fac(:)
      rate(:,361) = 1.15e-11_r8 * exp_fac(:)
      rate(:,362) = 4e-14_r8 * exp_fac(:)
      rate(:,375) = 1.45e-10_r8 * exp_fac(:)
      rate(:,376) = 3e-12_r8 * exp_fac(:)
      rate(:,377) = 6.7e-13_r8 * exp_fac(:)
      rate(:,387) = 3.5e-13_r8 * exp_fac(:)
      rate(:,388) = 5.4e-11_r8 * exp_fac(:)
      rate(:,391) = 2e-12_r8 * exp_fac(:)
      rate(:,392) = 1.4e-11_r8 * exp_fac(:)
      rate(:,395) = 2.4e-12_r8 * exp_fac(:)
      rate(:,406) = 5e-12_r8 * exp_fac(:)
      rate(:,416) = 1.6e-12_r8 * exp_fac(:)
      rate(:,418) = 6.7e-12_r8 * exp_fac(:)
      rate(:,421) = 3.5e-12_r8 * exp_fac(:)
      rate(:,424) = 1.3e-11_r8 * exp_fac(:)
      rate(:,425) = 1.4e-11_r8 * exp_fac(:)
      rate(:,429) = 2.4e-12_r8 * exp_fac(:)
      rate(:,430) = 1.4e-11_r8 * exp_fac(:)
      rate(:,435) = 2.4e-12_r8 * exp_fac(:)
      rate(:,436) = 4e-11_r8 * exp_fac(:)
      rate(:,437) = 4e-11_r8 * exp_fac(:)
      rate(:,439) = 1.4e-11_r8 * exp_fac(:)
      rate(:,443) = 2.4e-12_r8 * exp_fac(:)
      rate(:,444) = 4e-11_r8 * exp_fac(:)
      rate(:,448) = 7e-11_r8 * exp_fac(:)
      rate(:,449) = 1e-10_r8 * exp_fac(:)
      rate(:,454) = 2.4e-12_r8 * exp_fac(:)
      rate(:,469) = 4.7e-11_r8 * exp_fac(:)
      rate(:,482) = 2.1e-12_r8 * exp_fac(:)
      rate(:,483) = 2.8e-13_r8 * exp_fac(:)
      rate(:,491) = 1.7e-11_r8 * exp_fac(:)
      rate(:,497) = 8.4e-11_r8 * exp_fac(:)
      rate(:,499) = 1.9e-11_r8 * exp_fac(:)
      rate(:,500) = 1.2e-14_r8 * exp_fac(:)
      rate(:,501) = 2e-10_r8 * exp_fac(:)
      rate(:,508) = 2.4e-12_r8 * exp_fac(:)
      rate(:,509) = 2e-11_r8 * exp_fac(:)
      rate(:,513) = 2.3e-11_r8 * exp_fac(:)
      rate(:,514) = 2e-11_r8 * exp_fac(:)
      rate(:,518) = 3.3e-11_r8 * exp_fac(:)
      rate(:,519) = 1e-12_r8 * exp_fac(:)
      rate(:,520) = 5.7e-11_r8 * exp_fac(:)
      rate(:,521) = 3.4e-11_r8 * exp_fac(:)
      rate(:,526) = 2.3e-12_r8 * exp_fac(:)
      rate(:,528) = 1.2e-11_r8 * exp_fac(:)
      rate(:,529) = 5.7e-11_r8 * exp_fac(:)
      rate(:,530) = 2.8e-11_r8 * exp_fac(:)
      rate(:,531) = 6.6e-11_r8 * exp_fac(:)
      rate(:,532) = 1.4e-11_r8 * exp_fac(:)
      rate(:,535) = 1.9e-12_r8 * exp_fac(:)
      rate(:,547) = 6.34e-08_r8 * exp_fac(:)
      rate(:,553) = 1.9e-11_r8 * exp_fac(:)
      rate(:,556) = 1.2e-14_r8 * exp_fac(:)
      rate(:,557) = 2e-10_r8 * exp_fac(:)
      rate(:,568) = 1.34e-11_r8 * exp_fac(:)
      rate(:,571) = 1.34e-11_r8 * exp_fac(:)
      rate(:,577) = 1.34e-11_r8 * exp_fac(:)
      rate(:,578) = 1.34e-11_r8 * exp_fac(:)
      rate(:,583) = 1.7e-11_r8 * exp_fac(:)
      rate(:,606) = 6e-11_r8 * exp_fac(:)
      rate(:,609) = 1e-12_r8 * exp_fac(:)
      rate(:,610) = 4e-10_r8 * exp_fac(:)
      rate(:,611) = 2e-10_r8 * exp_fac(:)
      rate(:,612) = 1e-10_r8 * exp_fac(:)
      rate(:,613) = 5e-16_r8 * exp_fac(:)
      rate(:,614) = 4.4e-10_r8 * exp_fac(:)
      rate(:,615) = 9e-10_r8 * exp_fac(:)
      rate(:,618) = 1.29e-07_r8 * exp_fac(:)
      rate(:,619) = 2.31e-07_r8 * exp_fac(:)
      rate(:,620) = 2.31e-06_r8 * exp_fac(:)
      rate(:,621) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,294) = 6e-12_r8 * exp_fac(:)
      rate(:,393) = 5e-13_r8 * exp_fac(:)
      rate(:,426) = 5e-13_r8 * exp_fac(:)
      rate(:,431) = 5e-13_r8 * exp_fac(:)
      rate(:,440) = 5e-13_r8 * exp_fac(:)
      rate(:,451) = 5e-13_r8 * exp_fac(:)
      rate(:,299) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,300) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,301) = 1.64e-12_r8 * exp_fac(:)
      rate(:,412) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,302) = 2.03e-11_r8 * exp_fac(:)
      rate(:,534) = 3.4e-12_r8 * exp_fac(:)
      rate(:,303) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,304) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,305) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,306) = 1.25e-12_r8 * exp_fac(:)
      rate(:,316) = 3.4e-11_r8 * exp_fac(:)
      rate(:,307) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,308) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,314) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,315) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,318) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,319) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,320) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,321) = 2.8e-12_r8 * exp_fac(:)
      rate(:,383) = 2.9e-12_r8 * exp_fac(:)
      rate(:,322) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,324) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,327) = 7.5e-13_r8 * exp_fac(:)
      rate(:,341) = 7.5e-13_r8 * exp_fac(:)
      rate(:,356) = 7.5e-13_r8 * exp_fac(:)
      rate(:,369) = 7.5e-13_r8 * exp_fac(:)
      rate(:,378) = 7.5e-13_r8 * exp_fac(:)
      rate(:,382) = 8.6e-13_r8 * exp_fac(:)
      rate(:,394) = 8e-13_r8 * exp_fac(:)
      rate(:,407) = 7.5e-13_r8 * exp_fac(:)
      rate(:,417) = 7.5e-13_r8 * exp_fac(:)
      rate(:,427) = 8e-13_r8 * exp_fac(:)
      rate(:,432) = 8e-13_r8 * exp_fac(:)
      rate(:,441) = 8e-13_r8 * exp_fac(:)
      rate(:,452) = 8e-13_r8 * exp_fac(:)
      rate(:,459) = 7.5e-13_r8 * exp_fac(:)
      rate(:,463) = 7.5e-13_r8 * exp_fac(:)
      rate(:,466) = 7.5e-13_r8 * exp_fac(:)
      rate(:,479) = 7.5e-13_r8 * exp_fac(:)
      rate(:,486) = 7.5e-13_r8 * exp_fac(:)
      rate(:,492) = 7.5e-13_r8 * exp_fac(:)
      rate(:,495) = 7.5e-13_r8 * exp_fac(:)
      rate(:,506) = 7.5e-13_r8 * exp_fac(:)
      rate(:,511) = 7.5e-13_r8 * exp_fac(:)
      rate(:,516) = 7.5e-13_r8 * exp_fac(:)
      rate(:,559) = 7.5e-13_r8 * exp_fac(:)
      rate(:,566) = 7.5e-13_r8 * exp_fac(:)
      rate(:,569) = 7.5e-13_r8 * exp_fac(:)
      rate(:,580) = 7.5e-13_r8 * exp_fac(:)
      rate(:,584) = 7.5e-13_r8 * exp_fac(:)
      rate(:,328) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,329) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,333) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,338) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,342) = 2.6e-12_r8 * exp_fac(:)
      rate(:,460) = 2.6e-12_r8 * exp_fac(:)
      rate(:,465) = 2.6e-12_r8 * exp_fac(:)
      rate(:,467) = 2.6e-12_r8 * exp_fac(:)
      rate(:,480) = 2.6e-12_r8 * exp_fac(:)
      rate(:,487) = 2.6e-12_r8 * exp_fac(:)
      rate(:,493) = 2.6e-12_r8 * exp_fac(:)
      rate(:,496) = 2.6e-12_r8 * exp_fac(:)
      rate(:,560) = 2.6e-12_r8 * exp_fac(:)
      rate(:,567) = 2.6e-12_r8 * exp_fac(:)
      rate(:,570) = 2.6e-12_r8 * exp_fac(:)
      rate(:,581) = 2.6e-12_r8 * exp_fac(:)
      rate(:,585) = 2.6e-12_r8 * exp_fac(:)
      rate(:,343) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,345) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,346) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,347) = 1.4e-12_r8 * exp_fac(:)
      rate(:,367) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,348) = 4.63e-12_r8 * exp_fac(:)
      rate(:,563) = 2.7e-12_r8 * exp_fac(:)
      rate(:,349) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,350) = 2.9e-12_r8 * exp_fac(:)
      rate(:,351) = 2e-12_r8 * exp_fac(:)
      rate(:,381) = 7.1e-13_r8 * exp_fac(:)
      rate(:,402) = 2e-12_r8 * exp_fac(:)
      rate(:,505) = 2e-12_r8 * exp_fac(:)
      rate(:,510) = 2e-12_r8 * exp_fac(:)
      rate(:,515) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,352) = 4.3e-13_r8 * exp_fac(:)
      rate(:,403) = 4.3e-13_r8 * exp_fac(:)
      rate(:,456) = 4.3e-13_r8 * exp_fac(:)
      rate(:,470) = 4.3e-13_r8 * exp_fac(:)
      rate(:,473) = 4.3e-13_r8 * exp_fac(:)
      rate(:,476) = 4.3e-13_r8 * exp_fac(:)
      rate(:,354) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,358) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,366) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,368) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,372) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,373) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,374) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,389) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,390) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,396) = 2.7e-12_r8 * exp_fac(:)
      rate(:,397) = 1.3e-13_r8 * exp_fac(:)
      rate(:,399) = 9.6e-12_r8 * exp_fac(:)
      rate(:,405) = 5.3e-12_r8 * exp_fac(:)
      rate(:,442) = 2.7e-12_r8 * exp_fac(:)
      rate(:,453) = 2.7e-12_r8 * exp_fac(:)
      rate(:,555) = 2.7e-12_r8 * exp_fac(:)
      rate(:,574) = 2.7e-12_r8 * exp_fac(:)
      rate(:,398) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,401) = 4.6e-12_r8 * exp_fac(:)
      rate(:,404) = 2.3e-12_r8 * exp_fac(:)
      rate(:,409) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,413) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,419) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,422) = 1.86e-11_r8 * exp_fac(:)
      rate(:,423) = 1.86e-11_r8 * exp_fac(:)
      rate(:,433) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,438) = 3.03e-12_r8 * exp_fac(:)
      rate(:,561) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,446) = 2.54e-11_r8 * exp_fac(:)
      rate(:,565) = 2.54e-11_r8 * exp_fac(:)
      rate(:,450) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,458) = 2.3e-12_r8 * exp_fac(:)
      rate(:,558) = 2.3e-12_r8 * exp_fac(:)
      rate(:,462) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,481) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,489) = 1.7e-12_r8 * exp_fac(:)
      rate(:,579) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,502) = 1.2e-12_r8 * exp_fac(:)
      rate(:,572) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,503) = 6.3e-16_r8 * exp_fac(:)
      rate(:,575) = 6.3e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,504) = 1.2e-11_r8 * exp_fac(:)
      rate(:,576) = 1.2e-11_r8 * exp_fac(:)
      rate(:,522) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,523) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,524) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,525) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,533) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,536) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,539) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,554) = 2.75e-13_r8 * exp_fac(:)
      rate(:,562) = 2.12e-13_r8 * exp_fac(:)
      rate(:,573) = 2.6e-13_r8 * exp_fac(:)

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,195), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,205), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,217), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,225), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,228), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,229), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,230), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,248), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,268), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,279), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,325), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,335), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,336), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,337), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,363), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,364), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,385), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,411), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,414), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,472), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,475), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,478), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,485), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,527), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,179) = 1e-20_r8
      rate(:n,180) = 1.3e-16_r8
      rate(:n,184) = 8e-14_r8
      rate(:n,185) = 3.9e-17_r8
      rate(:n,192) = 6.9e-12_r8
      rate(:n,208) = 7e-13_r8
      rate(:n,209) = 5e-12_r8
      rate(:n,606) = 6e-11_r8
      rate(:n,609) = 1e-12_r8
      rate(:n,610) = 4e-10_r8
      rate(:n,611) = 2e-10_r8
      rate(:n,612) = 1e-10_r8
      rate(:n,614) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,175) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,176) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,177) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,181) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,183) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,186) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,187) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,196) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,197) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,198) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,201) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,202) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,203) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,210) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,214) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,215) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,223) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,224) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,195) = wrk(:)
























      end subroutine setrxt_hrates

      end module mo_setrxt
