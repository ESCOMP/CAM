
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

      rate(:,180) = 0.000258_r8
      rate(:,181) = 0.085_r8
      rate(:,182) = 1.2e-10_r8
      rate(:,187) = 1.2e-10_r8
      rate(:,188) = 1e-20_r8
      rate(:,189) = 1.3e-16_r8
      rate(:,191) = 4.2e-13_r8
      rate(:,193) = 8e-14_r8
      rate(:,194) = 3.9e-17_r8
      rate(:,201) = 6.9e-12_r8
      rate(:,202) = 7.2e-11_r8
      rate(:,203) = 1.6e-12_r8
      rate(:,209) = 1.8e-12_r8
      rate(:,213) = 1.8e-12_r8
      rate(:,217) = 7e-13_r8
      rate(:,218) = 5e-12_r8
      rate(:,227) = 3.5e-12_r8
      rate(:,229) = 1.3e-11_r8
      rate(:,230) = 2.2e-11_r8
      rate(:,231) = 5e-11_r8
      rate(:,268) = 1.7e-13_r8
      rate(:,270) = 2.607e-10_r8
      rate(:,271) = 9.75e-11_r8
      rate(:,272) = 2.07e-10_r8
      rate(:,273) = 2.088e-10_r8
      rate(:,274) = 1.17e-10_r8
      rate(:,275) = 4.644e-11_r8
      rate(:,276) = 1.204e-10_r8
      rate(:,277) = 9.9e-11_r8
      rate(:,278) = 3.3e-12_r8
      rate(:,301) = 4.5e-11_r8
      rate(:,302) = 4.62e-10_r8
      rate(:,303) = 1.2e-10_r8
      rate(:,304) = 9e-11_r8
      rate(:,305) = 3e-11_r8
      rate(:,307) = 2e-13_r8
      rate(:,308) = 1.5e-12_r8
      rate(:,309) = 1.25e-10_r8
      rate(:,310) = 1.8e-10_r8
      rate(:,311) = 1.44e-11_r8
      rate(:,317) = 1e-10_r8
      rate(:,321) = 2.49e-11_r8
      rate(:,330) = 9e-12_r8
      rate(:,331) = 1.4e-10_r8
      rate(:,332) = 3.6e-16_r8
      rate(:,333) = 1e-10_r8
      rate(:,349) = 2.14e-11_r8
      rate(:,350) = 1.9e-10_r8
      rate(:,353) = 1.3e-12_r8
      rate(:,371) = 1.2e-12_r8
      rate(:,372) = 8e-13_r8
      rate(:,376) = 2.3e-12_r8
      rate(:,382) = 2.57e-10_r8
      rate(:,383) = 1.8e-10_r8
      rate(:,384) = 1.794e-10_r8
      rate(:,385) = 1.3e-10_r8
      rate(:,386) = 7.65e-11_r8
      rate(:,399) = 4e-13_r8
      rate(:,403) = 1.31e-10_r8
      rate(:,404) = 3.5e-11_r8
      rate(:,405) = 9e-12_r8
      rate(:,412) = 6.8e-14_r8
      rate(:,413) = 2e-13_r8
      rate(:,428) = 1e-12_r8
      rate(:,432) = 1e-14_r8
      rate(:,433) = 1e-11_r8
      rate(:,434) = 1.15e-11_r8
      rate(:,435) = 4e-14_r8
      rate(:,448) = 1.45e-10_r8
      rate(:,449) = 3e-12_r8
      rate(:,450) = 6.7e-13_r8
      rate(:,460) = 3.5e-13_r8
      rate(:,461) = 5.4e-11_r8
      rate(:,464) = 2e-12_r8
      rate(:,465) = 1.4e-11_r8
      rate(:,468) = 2.4e-12_r8
      rate(:,479) = 5e-12_r8
      rate(:,489) = 2.2e-12_r8
      rate(:,491) = 6.7e-12_r8
      rate(:,494) = 3.5e-12_r8
      rate(:,497) = 1.3e-11_r8
      rate(:,498) = 1.4e-11_r8
      rate(:,502) = 2.4e-12_r8
      rate(:,503) = 1.4e-11_r8
      rate(:,508) = 2.4e-12_r8
      rate(:,509) = 4e-11_r8
      rate(:,510) = 4e-11_r8
      rate(:,512) = 1.4e-11_r8
      rate(:,516) = 2.4e-12_r8
      rate(:,517) = 4e-11_r8
      rate(:,521) = 7e-11_r8
      rate(:,522) = 1e-10_r8
      rate(:,527) = 2.4e-12_r8
      rate(:,542) = 4.7e-11_r8
      rate(:,555) = 2.1e-12_r8
      rate(:,556) = 2.8e-13_r8
      rate(:,564) = 1.7e-11_r8
      rate(:,570) = 8.4e-11_r8
      rate(:,572) = 1.9e-11_r8
      rate(:,573) = 1.2e-14_r8
      rate(:,574) = 2e-10_r8
      rate(:,581) = 2.4e-12_r8
      rate(:,582) = 2e-11_r8
      rate(:,586) = 2.3e-11_r8
      rate(:,587) = 2e-11_r8
      rate(:,591) = 3.3e-11_r8
      rate(:,592) = 1e-12_r8
      rate(:,593) = 5.7e-11_r8
      rate(:,594) = 3.4e-11_r8
      rate(:,596) = 3.3e-10_r8
      rate(:,603) = 2.3e-12_r8
      rate(:,605) = 1.2e-11_r8
      rate(:,606) = 5.7e-11_r8
      rate(:,607) = 2.8e-11_r8
      rate(:,608) = 6.6e-11_r8
      rate(:,609) = 1.4e-11_r8
      rate(:,612) = 1.9e-12_r8
      rate(:,643) = 6.34e-08_r8
      rate(:,665) = 1.9e-11_r8
      rate(:,668) = 1.2e-14_r8
      rate(:,669) = 2e-10_r8
      rate(:,680) = 1.34e-11_r8
      rate(:,686) = 1.34e-11_r8
      rate(:,691) = 1.7e-11_r8
      rate(:,739) = 6e-11_r8
      rate(:,742) = 1e-12_r8
      rate(:,743) = 4e-10_r8
      rate(:,744) = 2e-10_r8
      rate(:,745) = 1e-10_r8
      rate(:,746) = 5e-16_r8
      rate(:,747) = 4.4e-10_r8
      rate(:,748) = 9e-10_r8
      rate(:,751) = 1.29e-07_r8
      rate(:,752) = 2.31e-07_r8
      rate(:,753) = 2.31e-06_r8
      rate(:,754) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,183) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,184) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,185) = 2.64e-11_r8 * exp_fac(:)
      rate(:,186) = 6.6e-12_r8 * exp_fac(:)
      rate(:,190) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,192) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,195) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,196) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,199) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,200) = 1.4e-12_r8 * exp_fac(:)
      rate(:,518) = 1.05e-14_r8 * exp_fac(:)
      rate(:,676) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,205) = 3e-11_r8 * exp_fac(:)
      rate(:,299) = 5.5e-12_r8 * exp_fac(:)
      rate(:,396) = 3.8e-12_r8 * exp_fac(:)
      rate(:,417) = 3.8e-12_r8 * exp_fac(:)
      rate(:,444) = 3.8e-12_r8 * exp_fac(:)
      rate(:,453) = 3.8e-12_r8 * exp_fac(:)
      rate(:,457) = 3.8e-12_r8 * exp_fac(:)
      rate(:,473) = 2.3e-11_r8 * exp_fac(:)
      rate(:,483) = 3.8e-12_r8 * exp_fac(:)
      rate(:,493) = 3.8e-12_r8 * exp_fac(:)
      rate(:,520) = 1.52e-11_r8 * exp_fac(:)
      rate(:,528) = 1.52e-12_r8 * exp_fac(:)
      rate(:,534) = 3.8e-12_r8 * exp_fac(:)
      rate(:,537) = 3.8e-12_r8 * exp_fac(:)
      rate(:,541) = 3.8e-12_r8 * exp_fac(:)
      rate(:,557) = 3.8e-12_r8 * exp_fac(:)
      rate(:,561) = 3.8e-12_r8 * exp_fac(:)
      rate(:,567) = 3.8e-12_r8 * exp_fac(:)
      rate(:,571) = 3.8e-12_r8 * exp_fac(:)
      rate(:,206) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,207) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,208) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,210) = 4.8e-11_r8 * exp_fac(:)
      rate(:,297) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,211) = 1.8e-11_r8 * exp_fac(:)
      rate(:,430) = 4.2e-12_r8 * exp_fac(:)
      rate(:,452) = 4.2e-12_r8 * exp_fac(:)
      rate(:,481) = 4.2e-12_r8 * exp_fac(:)
      rate(:,501) = 4.4e-12_r8 * exp_fac(:)
      rate(:,507) = 4.4e-12_r8 * exp_fac(:)
      rate(:,580) = 4.2e-12_r8 * exp_fac(:)
      rate(:,585) = 4.2e-12_r8 * exp_fac(:)
      rate(:,590) = 4.2e-12_r8 * exp_fac(:)
      rate(:,212) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,216) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,219) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,220) = 2.9e-12_r8 * exp_fac(:)
      rate(:,221) = 1.45e-12_r8 * exp_fac(:)
      rate(:,222) = 1.45e-12_r8 * exp_fac(:)
      rate(:,223) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,224) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,225) = 1.2e-13_r8 * exp_fac(:)
      rate(:,253) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,228) = 1.7e-11_r8 * exp_fac(:)
      rate(:,390) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,232) = 3.44e-12_r8 * exp_fac(:)
      rate(:,288) = 2.3e-12_r8 * exp_fac(:)
      rate(:,291) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,233) = 3e-12_r8 * exp_fac(:)
      rate(:,298) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,235) = 7.26e-11_r8 * exp_fac(:)
      rate(:,236) = 4.64e-11_r8 * exp_fac(:)
      rate(:,243) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,244) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,245) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,246) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,247) = 1.4e-11_r8 * exp_fac(:)
      rate(:,263) = 7.4e-12_r8 * exp_fac(:)
      rate(:,426) = 8.1e-12_r8 * exp_fac(:)
      rate(:,248) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,250) = 2.4e-12_r8 * exp( -1250._r8 * itemp(:) )
      rate(:,251) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,252) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,254) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,255) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,256) = 2.6e-12_r8 * exp_fac(:)
      rate(:,257) = 6.4e-12_r8 * exp_fac(:)
      rate(:,289) = 4.1e-13_r8 * exp_fac(:)
      rate(:,530) = 7.5e-12_r8 * exp_fac(:)
      rate(:,544) = 7.5e-12_r8 * exp_fac(:)
      rate(:,547) = 7.5e-12_r8 * exp_fac(:)
      rate(:,550) = 7.5e-12_r8 * exp_fac(:)
      rate(:,258) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,260) = 3.6e-12_r8 * exp_fac(:)
      rate(:,356) = 2e-12_r8 * exp_fac(:)
      rate(:,261) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,262) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,264) = 6e-13_r8 * exp_fac(:)
      rate(:,286) = 1.5e-12_r8 * exp_fac(:)
      rate(:,296) = 1.9e-11_r8 * exp_fac(:)
      rate(:,265) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,266) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,267) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,269) = 3e-12_r8 * exp_fac(:)
      rate(:,346) = 1.4e-10_r8 * exp_fac(:)
      rate(:,281) = 2.1e-11_r8 * exp( 240._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,282) = 1.7e-11_r8 * exp_fac(:)
      rate(:,355) = 6.3e-12_r8 * exp_fac(:)
      rate(:,283) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,285) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,287) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,290) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,292) = 1.78e-11_r8 * exp_fac(:)
      rate(:,415) = 2.6e-12_r8 * exp_fac(:)
      rate(:,533) = 2.6e-12_r8 * exp_fac(:)
      rate(:,538) = 2.6e-12_r8 * exp_fac(:)
      rate(:,540) = 2.6e-12_r8 * exp_fac(:)
      rate(:,553) = 2.6e-12_r8 * exp_fac(:)
      rate(:,560) = 2.6e-12_r8 * exp_fac(:)
      rate(:,566) = 2.6e-12_r8 * exp_fac(:)
      rate(:,569) = 2.6e-12_r8 * exp_fac(:)
      rate(:,672) = 2.6e-12_r8 * exp_fac(:)
      rate(:,679) = 2.6e-12_r8 * exp_fac(:)
      rate(:,689) = 2.6e-12_r8 * exp_fac(:)
      rate(:,693) = 2.6e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 215._r8 * itemp(:) )
      rate(:,293) = 6.28e-11_r8 * exp_fac(:)
      rate(:,295) = 1.9e-11_r8 * exp_fac(:)
      rate(:,300) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,306) = 1.3e-12_r8 * exp( -1830._r8 * itemp(:) )
      rate(:,312) = 1.5e-11_r8 * exp( -1090._r8 * itemp(:) )
      rate(:,313) = 9.1e-11_r8 * exp( -146._r8 * itemp(:) )
      rate(:,314) = 4.7e-13_r8 * exp( -1670._r8 * itemp(:) )
      rate(:,316) = 1.008e+18_r8 * exp( -13670._r8 * itemp(:) )
      rate(:,318) = 8.4e-11_r8 * exp( -2620._r8 * itemp(:) )
      rate(:,320) = 2.1e-11_r8 * exp( -830._r8 * itemp(:) )
      exp_fac(:) = exp( 510._r8 * itemp(:) )
      rate(:,322) = 3e-12_r8 * exp_fac(:)
      rate(:,323) = 1.2e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 280._r8 * itemp(:) )
      rate(:,324) = 2.585e-12_r8 * exp_fac(:)
      rate(:,325) = 1.175e-12_r8 * exp_fac(:)
      rate(:,326) = 9.4e-13_r8 * exp_fac(:)
      rate(:,327) = 1.3e-11_r8 * exp( 570._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,328) = 7.15e-12_r8 * exp_fac(:)
      rate(:,394) = 2.8e-12_r8 * exp_fac(:)
      rate(:,456) = 2.9e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,334) = 1.6e-11_r8 * exp_fac(:)
      rate(:,577) = 1.2e-11_r8 * exp_fac(:)
      rate(:,685) = 1.2e-11_r8 * exp_fac(:)
      rate(:,335) = 1.1e-12_r8 * exp( 542._r8 * itemp(:) )
      rate(:,345) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,347) = 1.4e-11_r8 * exp_fac(:)
      rate(:,349) = 2.14e-11_r8 * exp_fac(:)
      rate(:,350) = 1.9e-10_r8 * exp_fac(:)
      rate(:,353) = 1.3e-12_r8 * exp_fac(:)
      rate(:,371) = 1.2e-12_r8 * exp_fac(:)
      rate(:,372) = 8e-13_r8 * exp_fac(:)
      rate(:,376) = 2.3e-12_r8 * exp_fac(:)
      rate(:,382) = 2.57e-10_r8 * exp_fac(:)
      rate(:,383) = 1.8e-10_r8 * exp_fac(:)
      rate(:,384) = 1.794e-10_r8 * exp_fac(:)
      rate(:,385) = 1.3e-10_r8 * exp_fac(:)
      rate(:,386) = 7.65e-11_r8 * exp_fac(:)
      rate(:,399) = 4e-13_r8 * exp_fac(:)
      rate(:,403) = 1.31e-10_r8 * exp_fac(:)
      rate(:,404) = 3.5e-11_r8 * exp_fac(:)
      rate(:,405) = 9e-12_r8 * exp_fac(:)
      rate(:,412) = 6.8e-14_r8 * exp_fac(:)
      rate(:,413) = 2e-13_r8 * exp_fac(:)
      rate(:,428) = 1e-12_r8 * exp_fac(:)
      rate(:,432) = 1e-14_r8 * exp_fac(:)
      rate(:,433) = 1e-11_r8 * exp_fac(:)
      rate(:,434) = 1.15e-11_r8 * exp_fac(:)
      rate(:,435) = 4e-14_r8 * exp_fac(:)
      rate(:,448) = 1.45e-10_r8 * exp_fac(:)
      rate(:,449) = 3e-12_r8 * exp_fac(:)
      rate(:,450) = 6.7e-13_r8 * exp_fac(:)
      rate(:,460) = 3.5e-13_r8 * exp_fac(:)
      rate(:,461) = 5.4e-11_r8 * exp_fac(:)
      rate(:,464) = 2e-12_r8 * exp_fac(:)
      rate(:,465) = 1.4e-11_r8 * exp_fac(:)
      rate(:,468) = 2.4e-12_r8 * exp_fac(:)
      rate(:,479) = 5e-12_r8 * exp_fac(:)
      rate(:,489) = 2.2e-12_r8 * exp_fac(:)
      rate(:,491) = 6.7e-12_r8 * exp_fac(:)
      rate(:,494) = 3.5e-12_r8 * exp_fac(:)
      rate(:,497) = 1.3e-11_r8 * exp_fac(:)
      rate(:,498) = 1.4e-11_r8 * exp_fac(:)
      rate(:,502) = 2.4e-12_r8 * exp_fac(:)
      rate(:,503) = 1.4e-11_r8 * exp_fac(:)
      rate(:,508) = 2.4e-12_r8 * exp_fac(:)
      rate(:,509) = 4e-11_r8 * exp_fac(:)
      rate(:,510) = 4e-11_r8 * exp_fac(:)
      rate(:,512) = 1.4e-11_r8 * exp_fac(:)
      rate(:,516) = 2.4e-12_r8 * exp_fac(:)
      rate(:,517) = 4e-11_r8 * exp_fac(:)
      rate(:,521) = 7e-11_r8 * exp_fac(:)
      rate(:,522) = 1e-10_r8 * exp_fac(:)
      rate(:,527) = 2.4e-12_r8 * exp_fac(:)
      rate(:,542) = 4.7e-11_r8 * exp_fac(:)
      rate(:,555) = 2.1e-12_r8 * exp_fac(:)
      rate(:,556) = 2.8e-13_r8 * exp_fac(:)
      rate(:,564) = 1.7e-11_r8 * exp_fac(:)
      rate(:,570) = 8.4e-11_r8 * exp_fac(:)
      rate(:,572) = 1.9e-11_r8 * exp_fac(:)
      rate(:,573) = 1.2e-14_r8 * exp_fac(:)
      rate(:,574) = 2e-10_r8 * exp_fac(:)
      rate(:,581) = 2.4e-12_r8 * exp_fac(:)
      rate(:,582) = 2e-11_r8 * exp_fac(:)
      rate(:,586) = 2.3e-11_r8 * exp_fac(:)
      rate(:,587) = 2e-11_r8 * exp_fac(:)
      rate(:,591) = 3.3e-11_r8 * exp_fac(:)
      rate(:,592) = 1e-12_r8 * exp_fac(:)
      rate(:,593) = 5.7e-11_r8 * exp_fac(:)
      rate(:,594) = 3.4e-11_r8 * exp_fac(:)
      rate(:,596) = 3.3e-10_r8 * exp_fac(:)
      rate(:,603) = 2.3e-12_r8 * exp_fac(:)
      rate(:,605) = 1.2e-11_r8 * exp_fac(:)
      rate(:,606) = 5.7e-11_r8 * exp_fac(:)
      rate(:,607) = 2.8e-11_r8 * exp_fac(:)
      rate(:,608) = 6.6e-11_r8 * exp_fac(:)
      rate(:,609) = 1.4e-11_r8 * exp_fac(:)
      rate(:,612) = 1.9e-12_r8 * exp_fac(:)
      rate(:,643) = 6.34e-08_r8 * exp_fac(:)
      rate(:,665) = 1.9e-11_r8 * exp_fac(:)
      rate(:,668) = 1.2e-14_r8 * exp_fac(:)
      rate(:,669) = 2e-10_r8 * exp_fac(:)
      rate(:,680) = 1.34e-11_r8 * exp_fac(:)
      rate(:,686) = 1.34e-11_r8 * exp_fac(:)
      rate(:,691) = 1.7e-11_r8 * exp_fac(:)
      rate(:,739) = 6e-11_r8 * exp_fac(:)
      rate(:,742) = 1e-12_r8 * exp_fac(:)
      rate(:,743) = 4e-10_r8 * exp_fac(:)
      rate(:,744) = 2e-10_r8 * exp_fac(:)
      rate(:,745) = 1e-10_r8 * exp_fac(:)
      rate(:,746) = 5e-16_r8 * exp_fac(:)
      rate(:,747) = 4.4e-10_r8 * exp_fac(:)
      rate(:,748) = 9e-10_r8 * exp_fac(:)
      rate(:,751) = 1.29e-07_r8 * exp_fac(:)
      rate(:,752) = 2.31e-07_r8 * exp_fac(:)
      rate(:,753) = 2.31e-06_r8 * exp_fac(:)
      rate(:,754) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,348) = 6e-12_r8 * exp_fac(:)
      rate(:,466) = 5e-13_r8 * exp_fac(:)
      rate(:,499) = 5e-13_r8 * exp_fac(:)
      rate(:,504) = 5e-13_r8 * exp_fac(:)
      rate(:,513) = 5e-13_r8 * exp_fac(:)
      rate(:,524) = 5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( -990._r8 * itemp(:) )
      rate(:,352) = 4.7e-12_r8 * exp_fac(:)
      rate(:,377) = 3.3e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1150._r8 * itemp(:) )
      rate(:,354) = 1.14e-11_r8 * exp_fac(:)
      rate(:,361) = 1.42e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -880._r8 * itemp(:) )
      rate(:,357) = 2.1e-12_r8 * exp_fac(:)
      rate(:,359) = 1.92e-12_r8 * exp_fac(:)
      rate(:,358) = 7.4e-12_r8 * exp( -910._r8 * itemp(:) )
      rate(:,360) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,362) = 1.64e-12_r8 * exp_fac(:)
      rate(:,485) = 8.5e-16_r8 * exp_fac(:)
      rate(:,363) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,364) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,365) = 2.9e-11_r8 * exp( -1000._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,366) = 2.9e-12_r8 * exp_fac(:)
      rate(:,611) = 3.4e-12_r8 * exp_fac(:)
      rate(:,367) = 9e-13_r8 * exp( -420._r8 * itemp(:) )
      rate(:,368) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,369) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      rate(:,370) = 9.4e-13_r8 * exp( -510._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,373) = 3.92e-13_r8 * exp_fac(:)
      rate(:,374) = 1.68e-13_r8 * exp_fac(:)
      rate(:,400) = 7.5e-13_r8 * exp_fac(:)
      rate(:,414) = 7.5e-13_r8 * exp_fac(:)
      rate(:,429) = 7.5e-13_r8 * exp_fac(:)
      rate(:,451) = 7.5e-13_r8 * exp_fac(:)
      rate(:,455) = 8.6e-13_r8 * exp_fac(:)
      rate(:,467) = 8e-13_r8 * exp_fac(:)
      rate(:,480) = 7.5e-13_r8 * exp_fac(:)
      rate(:,490) = 7.5e-13_r8 * exp_fac(:)
      rate(:,500) = 8e-13_r8 * exp_fac(:)
      rate(:,505) = 8e-13_r8 * exp_fac(:)
      rate(:,514) = 8e-13_r8 * exp_fac(:)
      rate(:,525) = 8e-13_r8 * exp_fac(:)
      rate(:,532) = 7.5e-13_r8 * exp_fac(:)
      rate(:,536) = 7.5e-13_r8 * exp_fac(:)
      rate(:,539) = 7.5e-13_r8 * exp_fac(:)
      rate(:,552) = 7.5e-13_r8 * exp_fac(:)
      rate(:,559) = 7.5e-13_r8 * exp_fac(:)
      rate(:,565) = 7.5e-13_r8 * exp_fac(:)
      rate(:,568) = 7.5e-13_r8 * exp_fac(:)
      rate(:,579) = 7.5e-13_r8 * exp_fac(:)
      rate(:,584) = 7.5e-13_r8 * exp_fac(:)
      rate(:,589) = 7.5e-13_r8 * exp_fac(:)
      rate(:,671) = 7.5e-13_r8 * exp_fac(:)
      rate(:,678) = 7.5e-13_r8 * exp_fac(:)
      rate(:,688) = 7.5e-13_r8 * exp_fac(:)
      rate(:,692) = 7.5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,375) = 4.05e-12_r8 * exp_fac(:)
      rate(:,443) = 2.7e-12_r8 * exp_fac(:)
      rate(:,469) = 2.7e-12_r8 * exp_fac(:)
      rate(:,470) = 1.3e-13_r8 * exp_fac(:)
      rate(:,472) = 9.6e-12_r8 * exp_fac(:)
      rate(:,478) = 5.3e-12_r8 * exp_fac(:)
      rate(:,515) = 2.7e-12_r8 * exp_fac(:)
      rate(:,526) = 2.7e-12_r8 * exp_fac(:)
      rate(:,667) = 2.7e-12_r8 * exp_fac(:)
      rate(:,683) = 2.7e-12_r8 * exp_fac(:)
      rate(:,378) = 2.2e-12_r8 * exp( -920._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,379) = 1.25e-12_r8 * exp_fac(:)
      rate(:,389) = 3.4e-11_r8 * exp_fac(:)
      rate(:,380) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,381) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,387) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,388) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,391) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,392) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,393) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,395) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,397) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,401) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,402) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,406) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,411) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      rate(:,416) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,418) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,419) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,420) = 1.4e-12_r8 * exp_fac(:)
      rate(:,440) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,421) = 4.63e-12_r8 * exp_fac(:)
      rate(:,675) = 2.7e-12_r8 * exp_fac(:)
      rate(:,422) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,423) = 2.9e-12_r8 * exp_fac(:)
      rate(:,424) = 2e-12_r8 * exp_fac(:)
      rate(:,454) = 7.1e-13_r8 * exp_fac(:)
      rate(:,475) = 2e-12_r8 * exp_fac(:)
      rate(:,578) = 2e-12_r8 * exp_fac(:)
      rate(:,583) = 2e-12_r8 * exp_fac(:)
      rate(:,588) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,425) = 4.3e-13_r8 * exp_fac(:)
      rate(:,476) = 4.3e-13_r8 * exp_fac(:)
      rate(:,529) = 4.3e-13_r8 * exp_fac(:)
      rate(:,543) = 4.3e-13_r8 * exp_fac(:)
      rate(:,546) = 4.3e-13_r8 * exp_fac(:)
      rate(:,549) = 4.3e-13_r8 * exp_fac(:)
      rate(:,427) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,431) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,439) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,441) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,442) = 1.41e-13_r8 * exp_fac(:)
      rate(:,666) = 2.75e-13_r8 * exp_fac(:)
      rate(:,674) = 2.12e-13_r8 * exp_fac(:)
      rate(:,682) = 2.6e-13_r8 * exp_fac(:)
      rate(:,445) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,446) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,447) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,462) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,463) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,471) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,474) = 4.6e-12_r8 * exp_fac(:)
      rate(:,477) = 2.3e-12_r8 * exp_fac(:)
      rate(:,482) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,486) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,492) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,495) = 1.86e-11_r8 * exp_fac(:)
      rate(:,496) = 1.86e-11_r8 * exp_fac(:)
      rate(:,506) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,511) = 3.03e-12_r8 * exp_fac(:)
      rate(:,673) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,519) = 2.54e-11_r8 * exp_fac(:)
      rate(:,677) = 2.54e-11_r8 * exp_fac(:)
      rate(:,523) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,531) = 2.3e-12_r8 * exp_fac(:)
      rate(:,670) = 2.3e-12_r8 * exp_fac(:)
      rate(:,535) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,554) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,562) = 1.7e-12_r8 * exp_fac(:)
      rate(:,687) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,575) = 1.2e-12_r8 * exp_fac(:)
      rate(:,681) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,576) = 6.3e-16_r8 * exp_fac(:)
      rate(:,684) = 6.3e-16_r8 * exp_fac(:)
      rate(:,595) = 1e-14_r8 * exp( 950._r8 * itemp(:) )
      rate(:,597) = 9.4e-11_r8 * exp( 190._r8 * itemp(:) )
      rate(:,598) = 3.2e-13_r8 * exp( -925._r8 * itemp(:) )
      rate(:,599) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,600) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,601) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,602) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,610) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,613) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,633) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,204), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,214), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,226), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,234), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,237), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,238), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,239), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**2._r8
      kinf(:) = 1e-10_r8 * itemp(:)
      call jpl( rate(:,249), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,259), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,279), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.2e-31_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.7e-11_r8
      call jpl( rate(:,284), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,294), m, 0.6_r8, ko, kinf, n )

      ko(:) = 3e-31_r8 * itemp(:)**1._r8
      kinf(:) = 6.6e-11_r8
      call jpl( rate(:,315), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-32_r8 * itemp(:)**1._r8
      kinf(:) = 1.7e-11_r8
      call jpl( rate(:,319), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.5e-31_r8 * itemp(:)**3.5_r8
      kinf(:) = 7.6e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,329), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.4e-28_r8 * itemp(:)**8.5_r8
      kinf(:) = 4e-11_r8 * itemp(:)**1.2_r8
      call jpl( rate(:,351), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,398), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,408), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,409), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,410), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,436), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,437), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,458), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,484), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,487), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,545), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,548), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,551), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,558), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,604), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,188) = 1e-20_r8
      rate(:n,189) = 1.3e-16_r8
      rate(:n,193) = 8e-14_r8
      rate(:n,194) = 3.9e-17_r8
      rate(:n,201) = 6.9e-12_r8
      rate(:n,217) = 7e-13_r8
      rate(:n,218) = 5e-12_r8
      rate(:n,739) = 6e-11_r8
      rate(:n,742) = 1e-12_r8
      rate(:n,743) = 4e-10_r8
      rate(:n,744) = 2e-10_r8
      rate(:n,745) = 1e-10_r8
      rate(:n,747) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,184) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,185) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,186) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,190) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,192) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,195) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,196) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,205) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,206) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,207) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,210) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,211) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,212) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,219) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,223) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,224) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,232) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,233) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,204) = wrk(:)






























      end subroutine setrxt_hrates

      end module mo_setrxt
