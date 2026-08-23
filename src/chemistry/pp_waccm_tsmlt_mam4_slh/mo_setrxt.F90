
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

      rate(:,181) = 0.000258_r8
      rate(:,182) = 0.085_r8
      rate(:,183) = 1.2e-10_r8
      rate(:,188) = 1.2e-10_r8
      rate(:,189) = 1e-20_r8
      rate(:,190) = 1.3e-16_r8
      rate(:,192) = 4.2e-13_r8
      rate(:,194) = 8e-14_r8
      rate(:,195) = 3.9e-17_r8
      rate(:,202) = 6.9e-12_r8
      rate(:,203) = 7.2e-11_r8
      rate(:,204) = 1.6e-12_r8
      rate(:,210) = 1.8e-12_r8
      rate(:,214) = 1.8e-12_r8
      rate(:,218) = 7e-13_r8
      rate(:,219) = 5e-12_r8
      rate(:,228) = 3.5e-12_r8
      rate(:,230) = 1.3e-11_r8
      rate(:,231) = 2.2e-11_r8
      rate(:,232) = 5e-11_r8
      rate(:,269) = 1.7e-13_r8
      rate(:,271) = 2.607e-10_r8
      rate(:,272) = 9.75e-11_r8
      rate(:,273) = 2.07e-10_r8
      rate(:,274) = 2.088e-10_r8
      rate(:,275) = 1.17e-10_r8
      rate(:,276) = 4.644e-11_r8
      rate(:,277) = 1.204e-10_r8
      rate(:,278) = 9.9e-11_r8
      rate(:,279) = 3.3e-12_r8
      rate(:,302) = 4.5e-11_r8
      rate(:,303) = 4.62e-10_r8
      rate(:,304) = 1.2e-10_r8
      rate(:,305) = 9e-11_r8
      rate(:,306) = 3e-11_r8
      rate(:,308) = 2e-13_r8
      rate(:,309) = 1.5e-12_r8
      rate(:,310) = 1.25e-10_r8
      rate(:,311) = 1.8e-10_r8
      rate(:,312) = 1.44e-11_r8
      rate(:,318) = 1e-10_r8
      rate(:,322) = 2.49e-11_r8
      rate(:,331) = 9e-12_r8
      rate(:,332) = 1.4e-10_r8
      rate(:,333) = 3.6e-16_r8
      rate(:,334) = 1e-10_r8
      rate(:,350) = 2.14e-11_r8
      rate(:,351) = 1.9e-10_r8
      rate(:,354) = 1.3e-12_r8
      rate(:,372) = 1.2e-12_r8
      rate(:,373) = 8e-13_r8
      rate(:,377) = 2.3e-12_r8
      rate(:,383) = 2.57e-10_r8
      rate(:,384) = 1.8e-10_r8
      rate(:,385) = 1.794e-10_r8
      rate(:,386) = 1.3e-10_r8
      rate(:,387) = 7.65e-11_r8
      rate(:,400) = 4e-13_r8
      rate(:,404) = 1.31e-10_r8
      rate(:,405) = 3.5e-11_r8
      rate(:,406) = 9e-12_r8
      rate(:,413) = 6.8e-14_r8
      rate(:,414) = 2e-13_r8
      rate(:,429) = 1e-12_r8
      rate(:,433) = 1e-14_r8
      rate(:,434) = 1e-11_r8
      rate(:,435) = 1.15e-11_r8
      rate(:,436) = 4e-14_r8
      rate(:,449) = 1.45e-10_r8
      rate(:,450) = 3e-12_r8
      rate(:,451) = 6.7e-13_r8
      rate(:,461) = 3.5e-13_r8
      rate(:,462) = 5.4e-11_r8
      rate(:,465) = 2e-12_r8
      rate(:,466) = 1.4e-11_r8
      rate(:,469) = 2.4e-12_r8
      rate(:,480) = 5e-12_r8
      rate(:,490) = 2.2e-12_r8
      rate(:,492) = 6.7e-12_r8
      rate(:,495) = 3.5e-12_r8
      rate(:,498) = 1.3e-11_r8
      rate(:,499) = 1.4e-11_r8
      rate(:,503) = 2.4e-12_r8
      rate(:,504) = 1.4e-11_r8
      rate(:,509) = 2.4e-12_r8
      rate(:,510) = 4e-11_r8
      rate(:,511) = 4e-11_r8
      rate(:,513) = 1.4e-11_r8
      rate(:,517) = 2.4e-12_r8
      rate(:,518) = 4e-11_r8
      rate(:,522) = 7e-11_r8
      rate(:,523) = 1e-10_r8
      rate(:,528) = 2.4e-12_r8
      rate(:,543) = 4.7e-11_r8
      rate(:,556) = 2.1e-12_r8
      rate(:,557) = 2.8e-13_r8
      rate(:,565) = 1.7e-11_r8
      rate(:,571) = 8.4e-11_r8
      rate(:,573) = 1.9e-11_r8
      rate(:,574) = 1.2e-14_r8
      rate(:,575) = 2e-10_r8
      rate(:,582) = 2.4e-12_r8
      rate(:,583) = 2e-11_r8
      rate(:,587) = 2.3e-11_r8
      rate(:,588) = 2e-11_r8
      rate(:,592) = 3.3e-11_r8
      rate(:,593) = 1e-12_r8
      rate(:,594) = 5.7e-11_r8
      rate(:,595) = 3.4e-11_r8
      rate(:,597) = 3.3e-10_r8
      rate(:,604) = 2.3e-12_r8
      rate(:,606) = 1.2e-11_r8
      rate(:,607) = 5.7e-11_r8
      rate(:,608) = 2.8e-11_r8
      rate(:,609) = 6.6e-11_r8
      rate(:,610) = 1.4e-11_r8
      rate(:,613) = 1.9e-12_r8
      rate(:,644) = 6.34e-08_r8
      rate(:,666) = 1.9e-11_r8
      rate(:,669) = 1.2e-14_r8
      rate(:,670) = 2e-10_r8
      rate(:,681) = 1.34e-11_r8
      rate(:,687) = 1.34e-11_r8
      rate(:,692) = 1.7e-11_r8
      rate(:,740) = 6e-11_r8
      rate(:,743) = 1e-12_r8
      rate(:,744) = 4e-10_r8
      rate(:,745) = 2e-10_r8
      rate(:,746) = 1e-10_r8
      rate(:,747) = 5e-16_r8
      rate(:,748) = 4.4e-10_r8
      rate(:,749) = 9e-10_r8
      rate(:,752) = 1.29e-07_r8
      rate(:,753) = 2.31e-07_r8
      rate(:,754) = 2.31e-06_r8
      rate(:,755) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,184) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,185) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,186) = 2.64e-11_r8 * exp_fac(:)
      rate(:,187) = 6.6e-12_r8 * exp_fac(:)
      rate(:,191) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,193) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,196) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,197) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,200) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,201) = 1.4e-12_r8 * exp_fac(:)
      rate(:,519) = 1.05e-14_r8 * exp_fac(:)
      rate(:,677) = 1.05e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,206) = 3e-11_r8 * exp_fac(:)
      rate(:,300) = 5.5e-12_r8 * exp_fac(:)
      rate(:,397) = 3.8e-12_r8 * exp_fac(:)
      rate(:,418) = 3.8e-12_r8 * exp_fac(:)
      rate(:,445) = 3.8e-12_r8 * exp_fac(:)
      rate(:,454) = 3.8e-12_r8 * exp_fac(:)
      rate(:,458) = 3.8e-12_r8 * exp_fac(:)
      rate(:,474) = 2.3e-11_r8 * exp_fac(:)
      rate(:,484) = 3.8e-12_r8 * exp_fac(:)
      rate(:,494) = 3.8e-12_r8 * exp_fac(:)
      rate(:,521) = 1.52e-11_r8 * exp_fac(:)
      rate(:,529) = 1.52e-12_r8 * exp_fac(:)
      rate(:,535) = 3.8e-12_r8 * exp_fac(:)
      rate(:,538) = 3.8e-12_r8 * exp_fac(:)
      rate(:,542) = 3.8e-12_r8 * exp_fac(:)
      rate(:,558) = 3.8e-12_r8 * exp_fac(:)
      rate(:,562) = 3.8e-12_r8 * exp_fac(:)
      rate(:,568) = 3.8e-12_r8 * exp_fac(:)
      rate(:,572) = 3.8e-12_r8 * exp_fac(:)
      rate(:,207) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,208) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,209) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,211) = 4.8e-11_r8 * exp_fac(:)
      rate(:,298) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,212) = 1.8e-11_r8 * exp_fac(:)
      rate(:,431) = 4.2e-12_r8 * exp_fac(:)
      rate(:,453) = 4.2e-12_r8 * exp_fac(:)
      rate(:,482) = 4.2e-12_r8 * exp_fac(:)
      rate(:,502) = 4.4e-12_r8 * exp_fac(:)
      rate(:,508) = 4.4e-12_r8 * exp_fac(:)
      rate(:,581) = 4.2e-12_r8 * exp_fac(:)
      rate(:,586) = 4.2e-12_r8 * exp_fac(:)
      rate(:,591) = 4.2e-12_r8 * exp_fac(:)
      rate(:,213) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,217) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,220) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,221) = 2.9e-12_r8 * exp_fac(:)
      rate(:,222) = 1.45e-12_r8 * exp_fac(:)
      rate(:,223) = 1.45e-12_r8 * exp_fac(:)
      rate(:,224) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,225) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,226) = 1.2e-13_r8 * exp_fac(:)
      rate(:,254) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,229) = 1.7e-11_r8 * exp_fac(:)
      rate(:,391) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,233) = 3.44e-12_r8 * exp_fac(:)
      rate(:,289) = 2.3e-12_r8 * exp_fac(:)
      rate(:,292) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,234) = 3e-12_r8 * exp_fac(:)
      rate(:,299) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,236) = 7.26e-11_r8 * exp_fac(:)
      rate(:,237) = 4.64e-11_r8 * exp_fac(:)
      rate(:,244) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,245) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,246) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,247) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,248) = 1.4e-11_r8 * exp_fac(:)
      rate(:,264) = 7.4e-12_r8 * exp_fac(:)
      rate(:,427) = 8.1e-12_r8 * exp_fac(:)
      rate(:,249) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,251) = 2.4e-12_r8 * exp( -1250._r8 * itemp(:) )
      rate(:,252) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,253) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,255) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,256) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,257) = 2.6e-12_r8 * exp_fac(:)
      rate(:,258) = 6.4e-12_r8 * exp_fac(:)
      rate(:,290) = 4.1e-13_r8 * exp_fac(:)
      rate(:,531) = 7.5e-12_r8 * exp_fac(:)
      rate(:,545) = 7.5e-12_r8 * exp_fac(:)
      rate(:,548) = 7.5e-12_r8 * exp_fac(:)
      rate(:,551) = 7.5e-12_r8 * exp_fac(:)
      rate(:,259) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,261) = 3.6e-12_r8 * exp_fac(:)
      rate(:,357) = 2e-12_r8 * exp_fac(:)
      rate(:,262) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,263) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,265) = 6e-13_r8 * exp_fac(:)
      rate(:,287) = 1.5e-12_r8 * exp_fac(:)
      rate(:,297) = 1.9e-11_r8 * exp_fac(:)
      rate(:,266) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,267) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,268) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,270) = 3e-12_r8 * exp_fac(:)
      rate(:,347) = 1.4e-10_r8 * exp_fac(:)
      rate(:,282) = 2.1e-11_r8 * exp( 240._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,283) = 1.7e-11_r8 * exp_fac(:)
      rate(:,356) = 6.3e-12_r8 * exp_fac(:)
      rate(:,284) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,286) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,288) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,291) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,293) = 1.78e-11_r8 * exp_fac(:)
      rate(:,416) = 2.6e-12_r8 * exp_fac(:)
      rate(:,534) = 2.6e-12_r8 * exp_fac(:)
      rate(:,539) = 2.6e-12_r8 * exp_fac(:)
      rate(:,541) = 2.6e-12_r8 * exp_fac(:)
      rate(:,554) = 2.6e-12_r8 * exp_fac(:)
      rate(:,561) = 2.6e-12_r8 * exp_fac(:)
      rate(:,567) = 2.6e-12_r8 * exp_fac(:)
      rate(:,570) = 2.6e-12_r8 * exp_fac(:)
      rate(:,673) = 2.6e-12_r8 * exp_fac(:)
      rate(:,680) = 2.6e-12_r8 * exp_fac(:)
      rate(:,690) = 2.6e-12_r8 * exp_fac(:)
      rate(:,694) = 2.6e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 215._r8 * itemp(:) )
      rate(:,294) = 6.28e-11_r8 * exp_fac(:)
      rate(:,296) = 1.9e-11_r8 * exp_fac(:)
      rate(:,301) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,307) = 1.3e-12_r8 * exp( -1830._r8 * itemp(:) )
      rate(:,313) = 1.5e-11_r8 * exp( -1090._r8 * itemp(:) )
      rate(:,314) = 9.1e-11_r8 * exp( -146._r8 * itemp(:) )
      rate(:,315) = 4.7e-13_r8 * exp( -1670._r8 * itemp(:) )
      rate(:,317) = 1.008e+18_r8 * exp( -13670._r8 * itemp(:) )
      rate(:,319) = 8.4e-11_r8 * exp( -2620._r8 * itemp(:) )
      rate(:,321) = 2.1e-11_r8 * exp( -830._r8 * itemp(:) )
      exp_fac(:) = exp( 510._r8 * itemp(:) )
      rate(:,323) = 3e-12_r8 * exp_fac(:)
      rate(:,324) = 1.2e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 280._r8 * itemp(:) )
      rate(:,325) = 2.585e-12_r8 * exp_fac(:)
      rate(:,326) = 1.175e-12_r8 * exp_fac(:)
      rate(:,327) = 9.4e-13_r8 * exp_fac(:)
      rate(:,328) = 1.3e-11_r8 * exp( 570._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,329) = 7.15e-12_r8 * exp_fac(:)
      rate(:,395) = 2.8e-12_r8 * exp_fac(:)
      rate(:,457) = 2.9e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 440._r8 * itemp(:) )
      rate(:,335) = 1.6e-11_r8 * exp_fac(:)
      rate(:,578) = 1.2e-11_r8 * exp_fac(:)
      rate(:,686) = 1.2e-11_r8 * exp_fac(:)
      rate(:,336) = 1.1e-12_r8 * exp( 542._r8 * itemp(:) )
      rate(:,346) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,348) = 1.4e-11_r8 * exp_fac(:)
      rate(:,350) = 2.14e-11_r8 * exp_fac(:)
      rate(:,351) = 1.9e-10_r8 * exp_fac(:)
      rate(:,354) = 1.3e-12_r8 * exp_fac(:)
      rate(:,372) = 1.2e-12_r8 * exp_fac(:)
      rate(:,373) = 8e-13_r8 * exp_fac(:)
      rate(:,377) = 2.3e-12_r8 * exp_fac(:)
      rate(:,383) = 2.57e-10_r8 * exp_fac(:)
      rate(:,384) = 1.8e-10_r8 * exp_fac(:)
      rate(:,385) = 1.794e-10_r8 * exp_fac(:)
      rate(:,386) = 1.3e-10_r8 * exp_fac(:)
      rate(:,387) = 7.65e-11_r8 * exp_fac(:)
      rate(:,400) = 4e-13_r8 * exp_fac(:)
      rate(:,404) = 1.31e-10_r8 * exp_fac(:)
      rate(:,405) = 3.5e-11_r8 * exp_fac(:)
      rate(:,406) = 9e-12_r8 * exp_fac(:)
      rate(:,413) = 6.8e-14_r8 * exp_fac(:)
      rate(:,414) = 2e-13_r8 * exp_fac(:)
      rate(:,429) = 1e-12_r8 * exp_fac(:)
      rate(:,433) = 1e-14_r8 * exp_fac(:)
      rate(:,434) = 1e-11_r8 * exp_fac(:)
      rate(:,435) = 1.15e-11_r8 * exp_fac(:)
      rate(:,436) = 4e-14_r8 * exp_fac(:)
      rate(:,449) = 1.45e-10_r8 * exp_fac(:)
      rate(:,450) = 3e-12_r8 * exp_fac(:)
      rate(:,451) = 6.7e-13_r8 * exp_fac(:)
      rate(:,461) = 3.5e-13_r8 * exp_fac(:)
      rate(:,462) = 5.4e-11_r8 * exp_fac(:)
      rate(:,465) = 2e-12_r8 * exp_fac(:)
      rate(:,466) = 1.4e-11_r8 * exp_fac(:)
      rate(:,469) = 2.4e-12_r8 * exp_fac(:)
      rate(:,480) = 5e-12_r8 * exp_fac(:)
      rate(:,490) = 2.2e-12_r8 * exp_fac(:)
      rate(:,492) = 6.7e-12_r8 * exp_fac(:)
      rate(:,495) = 3.5e-12_r8 * exp_fac(:)
      rate(:,498) = 1.3e-11_r8 * exp_fac(:)
      rate(:,499) = 1.4e-11_r8 * exp_fac(:)
      rate(:,503) = 2.4e-12_r8 * exp_fac(:)
      rate(:,504) = 1.4e-11_r8 * exp_fac(:)
      rate(:,509) = 2.4e-12_r8 * exp_fac(:)
      rate(:,510) = 4e-11_r8 * exp_fac(:)
      rate(:,511) = 4e-11_r8 * exp_fac(:)
      rate(:,513) = 1.4e-11_r8 * exp_fac(:)
      rate(:,517) = 2.4e-12_r8 * exp_fac(:)
      rate(:,518) = 4e-11_r8 * exp_fac(:)
      rate(:,522) = 7e-11_r8 * exp_fac(:)
      rate(:,523) = 1e-10_r8 * exp_fac(:)
      rate(:,528) = 2.4e-12_r8 * exp_fac(:)
      rate(:,543) = 4.7e-11_r8 * exp_fac(:)
      rate(:,556) = 2.1e-12_r8 * exp_fac(:)
      rate(:,557) = 2.8e-13_r8 * exp_fac(:)
      rate(:,565) = 1.7e-11_r8 * exp_fac(:)
      rate(:,571) = 8.4e-11_r8 * exp_fac(:)
      rate(:,573) = 1.9e-11_r8 * exp_fac(:)
      rate(:,574) = 1.2e-14_r8 * exp_fac(:)
      rate(:,575) = 2e-10_r8 * exp_fac(:)
      rate(:,582) = 2.4e-12_r8 * exp_fac(:)
      rate(:,583) = 2e-11_r8 * exp_fac(:)
      rate(:,587) = 2.3e-11_r8 * exp_fac(:)
      rate(:,588) = 2e-11_r8 * exp_fac(:)
      rate(:,592) = 3.3e-11_r8 * exp_fac(:)
      rate(:,593) = 1e-12_r8 * exp_fac(:)
      rate(:,594) = 5.7e-11_r8 * exp_fac(:)
      rate(:,595) = 3.4e-11_r8 * exp_fac(:)
      rate(:,597) = 3.3e-10_r8 * exp_fac(:)
      rate(:,604) = 2.3e-12_r8 * exp_fac(:)
      rate(:,606) = 1.2e-11_r8 * exp_fac(:)
      rate(:,607) = 5.7e-11_r8 * exp_fac(:)
      rate(:,608) = 2.8e-11_r8 * exp_fac(:)
      rate(:,609) = 6.6e-11_r8 * exp_fac(:)
      rate(:,610) = 1.4e-11_r8 * exp_fac(:)
      rate(:,613) = 1.9e-12_r8 * exp_fac(:)
      rate(:,644) = 6.34e-08_r8 * exp_fac(:)
      rate(:,666) = 1.9e-11_r8 * exp_fac(:)
      rate(:,669) = 1.2e-14_r8 * exp_fac(:)
      rate(:,670) = 2e-10_r8 * exp_fac(:)
      rate(:,681) = 1.34e-11_r8 * exp_fac(:)
      rate(:,687) = 1.34e-11_r8 * exp_fac(:)
      rate(:,692) = 1.7e-11_r8 * exp_fac(:)
      rate(:,740) = 6e-11_r8 * exp_fac(:)
      rate(:,743) = 1e-12_r8 * exp_fac(:)
      rate(:,744) = 4e-10_r8 * exp_fac(:)
      rate(:,745) = 2e-10_r8 * exp_fac(:)
      rate(:,746) = 1e-10_r8 * exp_fac(:)
      rate(:,747) = 5e-16_r8 * exp_fac(:)
      rate(:,748) = 4.4e-10_r8 * exp_fac(:)
      rate(:,749) = 9e-10_r8 * exp_fac(:)
      rate(:,752) = 1.29e-07_r8 * exp_fac(:)
      rate(:,753) = 2.31e-07_r8 * exp_fac(:)
      rate(:,754) = 2.31e-06_r8 * exp_fac(:)
      rate(:,755) = 4.63e-07_r8 * exp_fac(:)
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,349) = 6e-12_r8 * exp_fac(:)
      rate(:,467) = 5e-13_r8 * exp_fac(:)
      rate(:,500) = 5e-13_r8 * exp_fac(:)
      rate(:,505) = 5e-13_r8 * exp_fac(:)
      rate(:,514) = 5e-13_r8 * exp_fac(:)
      rate(:,525) = 5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( -990._r8 * itemp(:) )
      rate(:,353) = 4.7e-12_r8 * exp_fac(:)
      rate(:,378) = 3.3e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1150._r8 * itemp(:) )
      rate(:,355) = 1.14e-11_r8 * exp_fac(:)
      rate(:,362) = 1.42e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -880._r8 * itemp(:) )
      rate(:,358) = 2.1e-12_r8 * exp_fac(:)
      rate(:,360) = 1.92e-12_r8 * exp_fac(:)
      rate(:,359) = 7.4e-12_r8 * exp( -910._r8 * itemp(:) )
      rate(:,361) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,363) = 1.64e-12_r8 * exp_fac(:)
      rate(:,486) = 8.5e-16_r8 * exp_fac(:)
      rate(:,364) = 2.03e-11_r8 * exp( -1110._r8 * itemp(:) )
      rate(:,365) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,366) = 2.9e-11_r8 * exp( -1000._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,367) = 2.9e-12_r8 * exp_fac(:)
      rate(:,612) = 3.4e-12_r8 * exp_fac(:)
      rate(:,368) = 9e-13_r8 * exp( -420._r8 * itemp(:) )
      rate(:,369) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,370) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      rate(:,371) = 9.4e-13_r8 * exp( -510._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,374) = 3.92e-13_r8 * exp_fac(:)
      rate(:,375) = 1.68e-13_r8 * exp_fac(:)
      rate(:,401) = 7.5e-13_r8 * exp_fac(:)
      rate(:,415) = 7.5e-13_r8 * exp_fac(:)
      rate(:,430) = 7.5e-13_r8 * exp_fac(:)
      rate(:,452) = 7.5e-13_r8 * exp_fac(:)
      rate(:,456) = 8.6e-13_r8 * exp_fac(:)
      rate(:,468) = 8e-13_r8 * exp_fac(:)
      rate(:,481) = 7.5e-13_r8 * exp_fac(:)
      rate(:,491) = 7.5e-13_r8 * exp_fac(:)
      rate(:,501) = 8e-13_r8 * exp_fac(:)
      rate(:,506) = 8e-13_r8 * exp_fac(:)
      rate(:,515) = 8e-13_r8 * exp_fac(:)
      rate(:,526) = 8e-13_r8 * exp_fac(:)
      rate(:,533) = 7.5e-13_r8 * exp_fac(:)
      rate(:,537) = 7.5e-13_r8 * exp_fac(:)
      rate(:,540) = 7.5e-13_r8 * exp_fac(:)
      rate(:,553) = 7.5e-13_r8 * exp_fac(:)
      rate(:,560) = 7.5e-13_r8 * exp_fac(:)
      rate(:,566) = 7.5e-13_r8 * exp_fac(:)
      rate(:,569) = 7.5e-13_r8 * exp_fac(:)
      rate(:,580) = 7.5e-13_r8 * exp_fac(:)
      rate(:,585) = 7.5e-13_r8 * exp_fac(:)
      rate(:,590) = 7.5e-13_r8 * exp_fac(:)
      rate(:,672) = 7.5e-13_r8 * exp_fac(:)
      rate(:,679) = 7.5e-13_r8 * exp_fac(:)
      rate(:,689) = 7.5e-13_r8 * exp_fac(:)
      rate(:,693) = 7.5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,376) = 4.05e-12_r8 * exp_fac(:)
      rate(:,444) = 2.7e-12_r8 * exp_fac(:)
      rate(:,470) = 2.7e-12_r8 * exp_fac(:)
      rate(:,471) = 1.3e-13_r8 * exp_fac(:)
      rate(:,473) = 9.6e-12_r8 * exp_fac(:)
      rate(:,479) = 5.3e-12_r8 * exp_fac(:)
      rate(:,516) = 2.7e-12_r8 * exp_fac(:)
      rate(:,527) = 2.7e-12_r8 * exp_fac(:)
      rate(:,668) = 2.7e-12_r8 * exp_fac(:)
      rate(:,684) = 2.7e-12_r8 * exp_fac(:)
      rate(:,379) = 2.2e-12_r8 * exp( -920._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,380) = 1.25e-12_r8 * exp_fac(:)
      rate(:,390) = 3.4e-11_r8 * exp_fac(:)
      rate(:,381) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,382) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,388) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,389) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,392) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,393) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,394) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,396) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,398) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,402) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,403) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,407) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,412) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      rate(:,417) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,419) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,420) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,421) = 1.4e-12_r8 * exp_fac(:)
      rate(:,441) = 6.5e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 350._r8 * itemp(:) )
      rate(:,422) = 4.63e-12_r8 * exp_fac(:)
      rate(:,676) = 2.7e-12_r8 * exp_fac(:)
      rate(:,423) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,424) = 2.9e-12_r8 * exp_fac(:)
      rate(:,425) = 2e-12_r8 * exp_fac(:)
      rate(:,455) = 7.1e-13_r8 * exp_fac(:)
      rate(:,476) = 2e-12_r8 * exp_fac(:)
      rate(:,579) = 2e-12_r8 * exp_fac(:)
      rate(:,584) = 2e-12_r8 * exp_fac(:)
      rate(:,589) = 2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,426) = 4.3e-13_r8 * exp_fac(:)
      rate(:,477) = 4.3e-13_r8 * exp_fac(:)
      rate(:,530) = 4.3e-13_r8 * exp_fac(:)
      rate(:,544) = 4.3e-13_r8 * exp_fac(:)
      rate(:,547) = 4.3e-13_r8 * exp_fac(:)
      rate(:,550) = 4.3e-13_r8 * exp_fac(:)
      rate(:,428) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,432) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,440) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,442) = 1e-13_r8 * exp( 557._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,443) = 1.41e-13_r8 * exp_fac(:)
      rate(:,667) = 2.75e-13_r8 * exp_fac(:)
      rate(:,675) = 2.12e-13_r8 * exp_fac(:)
      rate(:,683) = 2.6e-13_r8 * exp_fac(:)
      rate(:,446) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,447) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,448) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,463) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,464) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,472) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,475) = 4.6e-12_r8 * exp_fac(:)
      rate(:,478) = 2.3e-12_r8 * exp_fac(:)
      rate(:,483) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,487) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,493) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,496) = 1.86e-11_r8 * exp_fac(:)
      rate(:,497) = 1.86e-11_r8 * exp_fac(:)
      rate(:,507) = 1.6e+09_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( -446._r8 * itemp(:) )
      rate(:,512) = 3.03e-12_r8 * exp_fac(:)
      rate(:,674) = 3.03e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,520) = 2.54e-11_r8 * exp_fac(:)
      rate(:,678) = 2.54e-11_r8 * exp_fac(:)
      rate(:,524) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,532) = 2.3e-12_r8 * exp_fac(:)
      rate(:,671) = 2.3e-12_r8 * exp_fac(:)
      rate(:,536) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,555) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,563) = 1.7e-12_r8 * exp_fac(:)
      rate(:,688) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,576) = 1.2e-12_r8 * exp_fac(:)
      rate(:,682) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,577) = 6.3e-16_r8 * exp_fac(:)
      rate(:,685) = 6.3e-16_r8 * exp_fac(:)
      rate(:,596) = 1e-14_r8 * exp( 950._r8 * itemp(:) )
      rate(:,598) = 9.4e-11_r8 * exp( 190._r8 * itemp(:) )
      rate(:,599) = 3.2e-13_r8 * exp( -925._r8 * itemp(:) )
      rate(:,600) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,601) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,602) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,603) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,611) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,614) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,634) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,205), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,215), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,227), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,235), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,238), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,239), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,240), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**2._r8
      kinf(:) = 1e-10_r8 * itemp(:)
      call jpl( rate(:,250), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,260), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,280), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.2e-31_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.7e-11_r8
      call jpl( rate(:,285), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,295), m, 0.6_r8, ko, kinf, n )

      ko(:) = 3e-31_r8 * itemp(:)**1._r8
      kinf(:) = 6.6e-11_r8
      call jpl( rate(:,316), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-32_r8 * itemp(:)**1._r8
      kinf(:) = 1.7e-11_r8
      call jpl( rate(:,320), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.5e-31_r8 * itemp(:)**3.5_r8
      kinf(:) = 7.6e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,330), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.4e-28_r8 * itemp(:)**8.5_r8
      kinf(:) = 4e-11_r8 * itemp(:)**1.2_r8
      call jpl( rate(:,352), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,399), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,409), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,410), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,411), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,437), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,438), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,459), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,485), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,488), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,546), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,549), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,552), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,559), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,605), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,189) = 1e-20_r8
      rate(:n,190) = 1.3e-16_r8
      rate(:n,194) = 8e-14_r8
      rate(:n,195) = 3.9e-17_r8
      rate(:n,202) = 6.9e-12_r8
      rate(:n,218) = 7e-13_r8
      rate(:n,219) = 5e-12_r8
      rate(:n,740) = 6e-11_r8
      rate(:n,743) = 1e-12_r8
      rate(:n,744) = 4e-10_r8
      rate(:n,745) = 2e-10_r8
      rate(:n,746) = 1e-10_r8
      rate(:n,748) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,185) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,186) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,187) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,191) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,193) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,196) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,197) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,206) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,207) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,208) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,211) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,212) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,213) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,220) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,224) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,225) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,233) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,234) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,205) = wrk(:)






























      end subroutine setrxt_hrates

      end module mo_setrxt
