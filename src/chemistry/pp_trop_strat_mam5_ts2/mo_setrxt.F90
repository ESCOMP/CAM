
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

      rate(:,168) = 1.2e-10_r8
      rate(:,172) = 1.2e-10_r8
      rate(:,178) = 6.9e-12_r8
      rate(:,179) = 7.2e-11_r8
      rate(:,180) = 1.6e-12_r8
      rate(:,186) = 1.8e-12_r8
      rate(:,190) = 1.8e-12_r8
      rate(:,202) = 3.5e-12_r8
      rate(:,204) = 1.3e-11_r8
      rate(:,205) = 2.2e-11_r8
      rate(:,206) = 5e-11_r8
      rate(:,241) = 1.7e-13_r8
      rate(:,243) = 2.607e-10_r8
      rate(:,244) = 9.75e-11_r8
      rate(:,245) = 2.07e-10_r8
      rate(:,246) = 2.088e-10_r8
      rate(:,247) = 1.17e-10_r8
      rate(:,248) = 4.644e-11_r8
      rate(:,249) = 1.204e-10_r8
      rate(:,250) = 9.9e-11_r8
      rate(:,251) = 3.3e-12_r8
      rate(:,270) = 4.5e-11_r8
      rate(:,271) = 4.62e-10_r8
      rate(:,272) = 1.2e-10_r8
      rate(:,273) = 9e-11_r8
      rate(:,274) = 3e-11_r8
      rate(:,279) = 2.14e-11_r8
      rate(:,280) = 1.9e-10_r8
      rate(:,293) = 2.57e-10_r8
      rate(:,294) = 1.8e-10_r8
      rate(:,295) = 1.794e-10_r8
      rate(:,296) = 1.3e-10_r8
      rate(:,297) = 7.65e-11_r8
      rate(:,310) = 4e-13_r8
      rate(:,315) = 1.31e-10_r8
      rate(:,316) = 3.5e-11_r8
      rate(:,317) = 9e-12_r8
      rate(:,324) = 6.8e-14_r8
      rate(:,325) = 2e-13_r8
      rate(:,340) = 1e-12_r8
      rate(:,344) = 1e-14_r8
      rate(:,345) = 1e-11_r8
      rate(:,346) = 1.15e-11_r8
      rate(:,347) = 3.3e-11_r8
      rate(:,348) = 3.4e-12_r8
      rate(:,349) = 4e-14_r8
      rate(:,362) = 3e-12_r8
      rate(:,363) = 1.2e-11_r8
      rate(:,364) = 6.7e-13_r8
      rate(:,374) = 3.5e-13_r8
      rate(:,375) = 5.4e-11_r8
      rate(:,376) = 3.77e-11_r8
      rate(:,379) = 2e-12_r8
      rate(:,380) = 1.29e-11_r8
      rate(:,382) = 4.5e-14_r8
      rate(:,387) = 3.77e-11_r8
      rate(:,393) = 4e-12_r8
      rate(:,399) = 1.78e-12_r8
      rate(:,401) = 6.1e-13_r8
      rate(:,405) = 4.8e-11_r8
      rate(:,408) = 1.6e-12_r8
      rate(:,410) = 6.7e-12_r8
      rate(:,413) = 3.5e-12_r8
      rate(:,418) = 6.42e-11_r8
      rate(:,425) = 1.6e-13_r8
      rate(:,431) = 1.4e-12_r8
      rate(:,436) = 7.5e-13_r8
      rate(:,437) = 1.4e-13_r8
      rate(:,438) = 7.5e-13_r8
      rate(:,439) = 3.6e-13_r8
      rate(:,440) = 6.5e-13_r8
      rate(:,441) = 2.1e-13_r8
      rate(:,442) = 6.5e-13_r8
      rate(:,443) = 4.9e-13_r8
      rate(:,445) = 1.2e-12_r8
      rate(:,449) = 9.8e-13_r8
      rate(:,452) = 1.85e-11_r8
      rate(:,453) = 1.63e-12_r8
      rate(:,454) = 2.5e-11_r8
      rate(:,455) = 1.1e-11_r8
      rate(:,456) = 3.3e-11_r8
      rate(:,459) = 2.8e-17_r8
      rate(:,460) = 8e-11_r8
      rate(:,463) = 3e-11_r8
      rate(:,466) = 4.2e-11_r8
      rate(:,469) = 2.8e-17_r8
      rate(:,470) = 1.1e-10_r8
      rate(:,472) = 3.9e-11_r8
      rate(:,475) = 1.3e-12_r8
      rate(:,477) = 5e-12_r8
      rate(:,478) = 2.3e-12_r8
      rate(:,481) = 3.9e-11_r8
      rate(:,484) = 2.8e-17_r8
      rate(:,485) = 9.2e-11_r8
      rate(:,488) = 3.85e-11_r8
      rate(:,492) = 1.2e-12_r8
      rate(:,496) = 9.8e-13_r8
      rate(:,501) = 4.4e-18_r8
      rate(:,502) = 3.6e-11_r8
      rate(:,554) = 4.7e-11_r8
      rate(:,567) = 2.1e-12_r8
      rate(:,568) = 2.8e-13_r8
      rate(:,576) = 1.7e-11_r8
      rate(:,582) = 8.4e-11_r8
      rate(:,585) = 5.3e-13_r8
      rate(:,587) = 2e-12_r8
      rate(:,590) = 2.3e-12_r8
      rate(:,595) = 2e-12_r8
      rate(:,598) = 2.3e-12_r8
      rate(:,604) = 1.9e-11_r8
      rate(:,605) = 5.3e-13_r8
      rate(:,607) = 2e-12_r8
      rate(:,610) = 2.3e-12_r8
      rate(:,615) = 2e-12_r8
      rate(:,618) = 2.3e-12_r8
      rate(:,622) = 1.2e-14_r8
      rate(:,623) = 2e-10_r8
      rate(:,624) = 2.5e-12_r8
      rate(:,625) = 5.3e-13_r8
      rate(:,627) = 2e-12_r8
      rate(:,630) = 2.3e-12_r8
      rate(:,635) = 2e-12_r8
      rate(:,638) = 2.3e-12_r8
      rate(:,644) = 1.2e-11_r8
      rate(:,646) = 2e-12_r8
      rate(:,648) = 5.3e-13_r8
      rate(:,650) = 2.3e-12_r8
      rate(:,655) = 2e-12_r8
      rate(:,658) = 2.3e-12_r8
      rate(:,664) = 1.1e-11_r8
      rate(:,666) = 2e-12_r8
      rate(:,668) = 5.3e-13_r8
      rate(:,670) = 2.3e-12_r8
      rate(:,675) = 2e-12_r8
      rate(:,678) = 2.3e-12_r8
      rate(:,683) = 2.1e-10_r8
      rate(:,689) = 8.9e-11_r8
      rate(:,690) = 8.9e-11_r8
      rate(:,694) = 2e-12_r8
      rate(:,697) = 2.3e-12_r8
      rate(:,705) = 4e-12_r8
      rate(:,708) = 2e-14_r8
      rate(:,710) = 2e-12_r8
      rate(:,713) = 2.3e-12_r8
      rate(:,718) = 2.52e-11_r8
      rate(:,723) = 4e-12_r8
      rate(:,727) = 2e-14_r8
      rate(:,729) = 2e-12_r8
      rate(:,732) = 2.3e-12_r8
      rate(:,737) = 1.92e-11_r8
      rate(:,739) = 2e-12_r8
      rate(:,742) = 2.3e-12_r8
      rate(:,746) = 8.8e-12_r8
      rate(:,747) = 8.8e-12_r8
      rate(:,748) = 8.8e-12_r8
      rate(:,753) = 4e-12_r8
      rate(:,755) = 2e-14_r8
      rate(:,757) = 3.66e-12_r8
      rate(:,758) = 2.8e-11_r8
      rate(:,759) = 2.6e-13_r8
      rate(:,762) = 8.3e-18_r8
      rate(:,763) = 1.1e-10_r8
      rate(:,767) = 1.1e-16_r8
      rate(:,769) = 3.64e-12_r8
      rate(:,770) = 2.8e-11_r8
      rate(:,771) = 1.7e-11_r8
      rate(:,774) = 1.1e-10_r8
      rate(:,775) = 9.58e-12_r8
      rate(:,778) = 1.1e-10_r8
      rate(:,779) = 1.23e-11_r8
      rate(:,782) = 1.1e-10_r8
      rate(:,783) = 3.64e-12_r8
      rate(:,786) = 1.1e-10_r8
      rate(:,787) = 5.5e-12_r8
      rate(:,788) = 4.65e-11_r8
      rate(:,789) = 2.8e-11_r8
      rate(:,797) = 2.3e-12_r8
      rate(:,799) = 1.2e-11_r8
      rate(:,800) = 5.7e-11_r8
      rate(:,801) = 2.8e-11_r8
      rate(:,802) = 6.6e-11_r8
      rate(:,803) = 1.4e-11_r8
      rate(:,806) = 1.9e-12_r8
      rate(:,829) = 6.34e-08_r8
      rate(:,846) = 1.9e-11_r8
      rate(:,849) = 1.2e-14_r8
      rate(:,850) = 2e-10_r8
      rate(:,854) = 2.5e-12_r8
      rate(:,866) = 1.34e-11_r8
      rate(:,867) = 1.2e-11_r8
      rate(:,872) = 1.1e-11_r8
      rate(:,876) = 2.1e-10_r8
      rate(:,877) = 1.34e-11_r8
      rate(:,881) = 1.7e-11_r8
      rate(:,901) = 1.29e-07_r8
      rate(:,902) = 2.31e-07_r8
      rate(:,903) = 2.31e-06_r8
      rate(:,904) = 4.63e-07_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,169) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,170) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:,171) = 3.3e-11_r8 * exp( 55._r8 * itemp(:) )
      rate(:,173) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,176) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,177) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,182) = 3e-11_r8 * exp_fac(:)
      rate(:,268) = 5.5e-12_r8 * exp_fac(:)
      rate(:,307) = 3.8e-12_r8 * exp_fac(:)
      rate(:,329) = 3.8e-12_r8 * exp_fac(:)
      rate(:,358) = 3.8e-12_r8 * exp_fac(:)
      rate(:,367) = 3.8e-12_r8 * exp_fac(:)
      rate(:,371) = 3.8e-12_r8 * exp_fac(:)
      rate(:,397) = 3.8e-12_r8 * exp_fac(:)
      rate(:,412) = 3.8e-12_r8 * exp_fac(:)
      rate(:,489) = 5.53e-12_r8 * exp_fac(:)
      rate(:,546) = 3.8e-12_r8 * exp_fac(:)
      rate(:,549) = 3.8e-12_r8 * exp_fac(:)
      rate(:,553) = 3.8e-12_r8 * exp_fac(:)
      rate(:,569) = 3.8e-12_r8 * exp_fac(:)
      rate(:,573) = 3.8e-12_r8 * exp_fac(:)
      rate(:,579) = 3.8e-12_r8 * exp_fac(:)
      rate(:,583) = 3.8e-12_r8 * exp_fac(:)
      rate(:,183) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,184) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,185) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,187) = 4.8e-11_r8 * exp_fac(:)
      rate(:,266) = 1.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,188) = 1.8e-11_r8 * exp_fac(:)
      rate(:,342) = 4.2e-12_r8 * exp_fac(:)
      rate(:,357) = 4.2e-12_r8 * exp_fac(:)
      rate(:,366) = 4.2e-12_r8 * exp_fac(:)
      rate(:,395) = 4.2e-12_r8 * exp_fac(:)
      rate(:,189) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,193) = 4.5e-13_r8 * exp( 610._r8 * itemp(:) )
      rate(:,194) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,195) = 2.9e-12_r8 * exp_fac(:)
      rate(:,196) = 1.45e-12_r8 * exp_fac(:)
      rate(:,197) = 1.45e-12_r8 * exp_fac(:)
      rate(:,198) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:,199) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,200) = 1.2e-13_r8 * exp_fac(:)
      rate(:,226) = 3e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 125._r8 * itemp(:) )
      rate(:,203) = 1.7e-11_r8 * exp_fac(:)
      rate(:,301) = 5.5e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,207) = 3.44e-12_r8 * exp_fac(:)
      rate(:,259) = 2.3e-12_r8 * exp_fac(:)
      rate(:,262) = 8.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,208) = 3e-12_r8 * exp_fac(:)
      rate(:,267) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,210) = 7.26e-11_r8 * exp_fac(:)
      rate(:,211) = 4.64e-11_r8 * exp_fac(:)
      rate(:,218) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      exp_fac(:) = exp( -1270._r8 * itemp(:) )
      rate(:,219) = 7.1e-12_r8 * exp_fac(:)
      rate(:,642) = 1.35e-15_r8 * exp_fac(:)
      rate(:,857) = 1.35e-15_r8 * exp_fac(:)
      rate(:,220) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,221) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,222) = 1.4e-11_r8 * exp_fac(:)
      rate(:,236) = 7.4e-12_r8 * exp_fac(:)
      rate(:,338) = 8.1e-12_r8 * exp_fac(:)
      rate(:,392) = 8.1e-12_r8 * exp_fac(:)
      rate(:,704) = 8.1e-12_r8 * exp_fac(:)
      rate(:,722) = 8.1e-12_r8 * exp_fac(:)
      rate(:,752) = 8.1e-12_r8 * exp_fac(:)
      rate(:,223) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,224) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,225) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,227) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,228) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,229) = 2.6e-12_r8 * exp_fac(:)
      rate(:,230) = 6.4e-12_r8 * exp_fac(:)
      rate(:,260) = 4.1e-13_r8 * exp_fac(:)
      rate(:,542) = 7.5e-12_r8 * exp_fac(:)
      rate(:,556) = 7.5e-12_r8 * exp_fac(:)
      rate(:,559) = 7.5e-12_r8 * exp_fac(:)
      rate(:,562) = 7.5e-12_r8 * exp_fac(:)
      rate(:,231) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,233) = 3.6e-12_r8 * exp_fac(:)
      rate(:,282) = 2e-12_r8 * exp_fac(:)
      rate(:,234) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,235) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,237) = 6e-13_r8 * exp_fac(:)
      rate(:,257) = 1.5e-12_r8 * exp_fac(:)
      rate(:,265) = 1.9e-11_r8 * exp_fac(:)
      rate(:,238) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,239) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,240) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,242) = 3e-12_r8 * exp_fac(:)
      rate(:,276) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,254) = 1.7e-11_r8 * exp_fac(:)
      rate(:,281) = 6.3e-12_r8 * exp_fac(:)
      rate(:,255) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,256) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,258) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 460._r8 * itemp(:) )
      rate(:,261) = 4.5e-12_r8 * exp_fac(:)
      rate(:,643) = 1.62e-11_r8 * exp_fac(:)
      rate(:,858) = 1.62e-11_r8 * exp_fac(:)
      rate(:,264) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,269) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,275) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,277) = 1.4e-11_r8 * exp_fac(:)
      rate(:,279) = 2.14e-11_r8 * exp_fac(:)
      rate(:,280) = 1.9e-10_r8 * exp_fac(:)
      rate(:,293) = 2.57e-10_r8 * exp_fac(:)
      rate(:,294) = 1.8e-10_r8 * exp_fac(:)
      rate(:,295) = 1.794e-10_r8 * exp_fac(:)
      rate(:,296) = 1.3e-10_r8 * exp_fac(:)
      rate(:,297) = 7.65e-11_r8 * exp_fac(:)
      rate(:,310) = 4e-13_r8 * exp_fac(:)
      rate(:,315) = 1.31e-10_r8 * exp_fac(:)
      rate(:,316) = 3.5e-11_r8 * exp_fac(:)
      rate(:,317) = 9e-12_r8 * exp_fac(:)
      rate(:,324) = 6.8e-14_r8 * exp_fac(:)
      rate(:,325) = 2e-13_r8 * exp_fac(:)
      rate(:,340) = 1e-12_r8 * exp_fac(:)
      rate(:,344) = 1e-14_r8 * exp_fac(:)
      rate(:,345) = 1e-11_r8 * exp_fac(:)
      rate(:,346) = 1.15e-11_r8 * exp_fac(:)
      rate(:,347) = 3.3e-11_r8 * exp_fac(:)
      rate(:,348) = 3.4e-12_r8 * exp_fac(:)
      rate(:,349) = 4e-14_r8 * exp_fac(:)
      rate(:,362) = 3e-12_r8 * exp_fac(:)
      rate(:,363) = 1.2e-11_r8 * exp_fac(:)
      rate(:,364) = 6.7e-13_r8 * exp_fac(:)
      rate(:,374) = 3.5e-13_r8 * exp_fac(:)
      rate(:,375) = 5.4e-11_r8 * exp_fac(:)
      rate(:,376) = 3.77e-11_r8 * exp_fac(:)
      rate(:,379) = 2e-12_r8 * exp_fac(:)
      rate(:,380) = 1.29e-11_r8 * exp_fac(:)
      rate(:,382) = 4.5e-14_r8 * exp_fac(:)
      rate(:,387) = 3.77e-11_r8 * exp_fac(:)
      rate(:,393) = 4e-12_r8 * exp_fac(:)
      rate(:,399) = 1.78e-12_r8 * exp_fac(:)
      rate(:,401) = 6.1e-13_r8 * exp_fac(:)
      rate(:,405) = 4.8e-11_r8 * exp_fac(:)
      rate(:,408) = 1.6e-12_r8 * exp_fac(:)
      rate(:,410) = 6.7e-12_r8 * exp_fac(:)
      rate(:,413) = 3.5e-12_r8 * exp_fac(:)
      rate(:,418) = 6.42e-11_r8 * exp_fac(:)
      rate(:,425) = 1.6e-13_r8 * exp_fac(:)
      rate(:,431) = 1.4e-12_r8 * exp_fac(:)
      rate(:,436) = 7.5e-13_r8 * exp_fac(:)
      rate(:,437) = 1.4e-13_r8 * exp_fac(:)
      rate(:,438) = 7.5e-13_r8 * exp_fac(:)
      rate(:,439) = 3.6e-13_r8 * exp_fac(:)
      rate(:,440) = 6.5e-13_r8 * exp_fac(:)
      rate(:,441) = 2.1e-13_r8 * exp_fac(:)
      rate(:,442) = 6.5e-13_r8 * exp_fac(:)
      rate(:,443) = 4.9e-13_r8 * exp_fac(:)
      rate(:,445) = 1.2e-12_r8 * exp_fac(:)
      rate(:,449) = 9.8e-13_r8 * exp_fac(:)
      rate(:,452) = 1.85e-11_r8 * exp_fac(:)
      rate(:,453) = 1.63e-12_r8 * exp_fac(:)
      rate(:,454) = 2.5e-11_r8 * exp_fac(:)
      rate(:,455) = 1.1e-11_r8 * exp_fac(:)
      rate(:,456) = 3.3e-11_r8 * exp_fac(:)
      rate(:,459) = 2.8e-17_r8 * exp_fac(:)
      rate(:,460) = 8e-11_r8 * exp_fac(:)
      rate(:,463) = 3e-11_r8 * exp_fac(:)
      rate(:,466) = 4.2e-11_r8 * exp_fac(:)
      rate(:,469) = 2.8e-17_r8 * exp_fac(:)
      rate(:,470) = 1.1e-10_r8 * exp_fac(:)
      rate(:,472) = 3.9e-11_r8 * exp_fac(:)
      rate(:,475) = 1.3e-12_r8 * exp_fac(:)
      rate(:,477) = 5e-12_r8 * exp_fac(:)
      rate(:,478) = 2.3e-12_r8 * exp_fac(:)
      rate(:,481) = 3.9e-11_r8 * exp_fac(:)
      rate(:,484) = 2.8e-17_r8 * exp_fac(:)
      rate(:,485) = 9.2e-11_r8 * exp_fac(:)
      rate(:,488) = 3.85e-11_r8 * exp_fac(:)
      rate(:,492) = 1.2e-12_r8 * exp_fac(:)
      rate(:,496) = 9.8e-13_r8 * exp_fac(:)
      rate(:,501) = 4.4e-18_r8 * exp_fac(:)
      rate(:,502) = 3.6e-11_r8 * exp_fac(:)
      rate(:,554) = 4.7e-11_r8 * exp_fac(:)
      rate(:,567) = 2.1e-12_r8 * exp_fac(:)
      rate(:,568) = 2.8e-13_r8 * exp_fac(:)
      rate(:,576) = 1.7e-11_r8 * exp_fac(:)
      rate(:,582) = 8.4e-11_r8 * exp_fac(:)
      rate(:,585) = 5.3e-13_r8 * exp_fac(:)
      rate(:,587) = 2e-12_r8 * exp_fac(:)
      rate(:,590) = 2.3e-12_r8 * exp_fac(:)
      rate(:,595) = 2e-12_r8 * exp_fac(:)
      rate(:,598) = 2.3e-12_r8 * exp_fac(:)
      rate(:,604) = 1.9e-11_r8 * exp_fac(:)
      rate(:,605) = 5.3e-13_r8 * exp_fac(:)
      rate(:,607) = 2e-12_r8 * exp_fac(:)
      rate(:,610) = 2.3e-12_r8 * exp_fac(:)
      rate(:,615) = 2e-12_r8 * exp_fac(:)
      rate(:,618) = 2.3e-12_r8 * exp_fac(:)
      rate(:,622) = 1.2e-14_r8 * exp_fac(:)
      rate(:,623) = 2e-10_r8 * exp_fac(:)
      rate(:,624) = 2.5e-12_r8 * exp_fac(:)
      rate(:,625) = 5.3e-13_r8 * exp_fac(:)
      rate(:,627) = 2e-12_r8 * exp_fac(:)
      rate(:,630) = 2.3e-12_r8 * exp_fac(:)
      rate(:,635) = 2e-12_r8 * exp_fac(:)
      rate(:,638) = 2.3e-12_r8 * exp_fac(:)
      rate(:,644) = 1.2e-11_r8 * exp_fac(:)
      rate(:,646) = 2e-12_r8 * exp_fac(:)
      rate(:,648) = 5.3e-13_r8 * exp_fac(:)
      rate(:,650) = 2.3e-12_r8 * exp_fac(:)
      rate(:,655) = 2e-12_r8 * exp_fac(:)
      rate(:,658) = 2.3e-12_r8 * exp_fac(:)
      rate(:,664) = 1.1e-11_r8 * exp_fac(:)
      rate(:,666) = 2e-12_r8 * exp_fac(:)
      rate(:,668) = 5.3e-13_r8 * exp_fac(:)
      rate(:,670) = 2.3e-12_r8 * exp_fac(:)
      rate(:,675) = 2e-12_r8 * exp_fac(:)
      rate(:,678) = 2.3e-12_r8 * exp_fac(:)
      rate(:,683) = 2.1e-10_r8 * exp_fac(:)
      rate(:,689) = 8.9e-11_r8 * exp_fac(:)
      rate(:,690) = 8.9e-11_r8 * exp_fac(:)
      rate(:,694) = 2e-12_r8 * exp_fac(:)
      rate(:,697) = 2.3e-12_r8 * exp_fac(:)
      rate(:,705) = 4e-12_r8 * exp_fac(:)
      rate(:,708) = 2e-14_r8 * exp_fac(:)
      rate(:,710) = 2e-12_r8 * exp_fac(:)
      rate(:,713) = 2.3e-12_r8 * exp_fac(:)
      rate(:,718) = 2.52e-11_r8 * exp_fac(:)
      rate(:,723) = 4e-12_r8 * exp_fac(:)
      rate(:,727) = 2e-14_r8 * exp_fac(:)
      rate(:,729) = 2e-12_r8 * exp_fac(:)
      rate(:,732) = 2.3e-12_r8 * exp_fac(:)
      rate(:,737) = 1.92e-11_r8 * exp_fac(:)
      rate(:,739) = 2e-12_r8 * exp_fac(:)
      rate(:,742) = 2.3e-12_r8 * exp_fac(:)
      rate(:,746) = 8.8e-12_r8 * exp_fac(:)
      rate(:,747) = 8.8e-12_r8 * exp_fac(:)
      rate(:,748) = 8.8e-12_r8 * exp_fac(:)
      rate(:,753) = 4e-12_r8 * exp_fac(:)
      rate(:,755) = 2e-14_r8 * exp_fac(:)
      rate(:,757) = 3.66e-12_r8 * exp_fac(:)
      rate(:,758) = 2.8e-11_r8 * exp_fac(:)
      rate(:,759) = 2.6e-13_r8 * exp_fac(:)
      rate(:,762) = 8.3e-18_r8 * exp_fac(:)
      rate(:,763) = 1.1e-10_r8 * exp_fac(:)
      rate(:,767) = 1.1e-16_r8 * exp_fac(:)
      rate(:,769) = 3.64e-12_r8 * exp_fac(:)
      rate(:,770) = 2.8e-11_r8 * exp_fac(:)
      rate(:,771) = 1.7e-11_r8 * exp_fac(:)
      rate(:,774) = 1.1e-10_r8 * exp_fac(:)
      rate(:,775) = 9.58e-12_r8 * exp_fac(:)
      rate(:,778) = 1.1e-10_r8 * exp_fac(:)
      rate(:,779) = 1.23e-11_r8 * exp_fac(:)
      rate(:,782) = 1.1e-10_r8 * exp_fac(:)
      rate(:,783) = 3.64e-12_r8 * exp_fac(:)
      rate(:,786) = 1.1e-10_r8 * exp_fac(:)
      rate(:,787) = 5.5e-12_r8 * exp_fac(:)
      rate(:,788) = 4.65e-11_r8 * exp_fac(:)
      rate(:,789) = 2.8e-11_r8 * exp_fac(:)
      rate(:,797) = 2.3e-12_r8 * exp_fac(:)
      rate(:,799) = 1.2e-11_r8 * exp_fac(:)
      rate(:,800) = 5.7e-11_r8 * exp_fac(:)
      rate(:,801) = 2.8e-11_r8 * exp_fac(:)
      rate(:,802) = 6.6e-11_r8 * exp_fac(:)
      rate(:,803) = 1.4e-11_r8 * exp_fac(:)
      rate(:,806) = 1.9e-12_r8 * exp_fac(:)
      rate(:,829) = 6.34e-08_r8 * exp_fac(:)
      rate(:,846) = 1.9e-11_r8 * exp_fac(:)
      rate(:,849) = 1.2e-14_r8 * exp_fac(:)
      rate(:,850) = 2e-10_r8 * exp_fac(:)
      rate(:,854) = 2.5e-12_r8 * exp_fac(:)
      rate(:,866) = 1.34e-11_r8 * exp_fac(:)
      rate(:,867) = 1.2e-11_r8 * exp_fac(:)
      rate(:,872) = 1.1e-11_r8 * exp_fac(:)
      rate(:,876) = 2.1e-10_r8 * exp_fac(:)
      rate(:,877) = 1.34e-11_r8 * exp_fac(:)
      rate(:,881) = 1.7e-11_r8 * exp_fac(:)
      rate(:,901) = 1.29e-07_r8 * exp_fac(:)
      rate(:,902) = 2.31e-07_r8 * exp_fac(:)
      rate(:,903) = 2.31e-06_r8 * exp_fac(:)
      rate(:,904) = 4.63e-07_r8 * exp_fac(:)
      rate(:,278) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,283) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,284) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,285) = 1.64e-12_r8 * exp_fac(:)
      rate(:,403) = 8.5e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,286) = 2.03e-11_r8 * exp_fac(:)
      rate(:,805) = 3.4e-12_r8 * exp_fac(:)
      rate(:,287) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,288) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,289) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,290) = 1.25e-12_r8 * exp_fac(:)
      rate(:,300) = 3.4e-11_r8 * exp_fac(:)
      rate(:,291) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,292) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,298) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      rate(:,299) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,302) = 5e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,303) = 1.9e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,304) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,305) = 2.8e-12_r8 * exp_fac(:)
      rate(:,370) = 2.9e-12_r8 * exp_fac(:)
      rate(:,306) = 2.9e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,308) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,311) = 1.3e-12_r8 * exp_fac(:)
      rate(:,335) = 2.9e-12_r8 * exp_fac(:)
      rate(:,336) = 2e-12_r8 * exp_fac(:)
      rate(:,368) = 7.1e-13_r8 * exp_fac(:)
      rate(:,381) = 2e-12_r8 * exp_fac(:)
      rate(:,388) = 2.9e-12_r8 * exp_fac(:)
      rate(:,389) = 2e-12_r8 * exp_fac(:)
      rate(:,391) = 2.9e-12_r8 * exp_fac(:)
      rate(:,400) = 2e-12_r8 * exp_fac(:)
      rate(:,424) = 2e-12_r8 * exp_fac(:)
      rate(:,430) = 2e-12_r8 * exp_fac(:)
      rate(:,444) = 2e-12_r8 * exp_fac(:)
      rate(:,448) = 2e-12_r8 * exp_fac(:)
      rate(:,474) = 2e-12_r8 * exp_fac(:)
      rate(:,491) = 2e-12_r8 * exp_fac(:)
      rate(:,495) = 2e-12_r8 * exp_fac(:)
      rate(:,586) = 2e-12_r8 * exp_fac(:)
      rate(:,591) = 2e-12_r8 * exp_fac(:)
      rate(:,592) = 2e-12_r8 * exp_fac(:)
      rate(:,593) = 2e-12_r8 * exp_fac(:)
      rate(:,594) = 2e-12_r8 * exp_fac(:)
      rate(:,599) = 2e-12_r8 * exp_fac(:)
      rate(:,600) = 2e-12_r8 * exp_fac(:)
      rate(:,601) = 2e-12_r8 * exp_fac(:)
      rate(:,606) = 2e-12_r8 * exp_fac(:)
      rate(:,611) = 2e-12_r8 * exp_fac(:)
      rate(:,612) = 2e-12_r8 * exp_fac(:)
      rate(:,613) = 2e-12_r8 * exp_fac(:)
      rate(:,614) = 2e-12_r8 * exp_fac(:)
      rate(:,619) = 2e-12_r8 * exp_fac(:)
      rate(:,620) = 2e-12_r8 * exp_fac(:)
      rate(:,621) = 2e-12_r8 * exp_fac(:)
      rate(:,626) = 2e-12_r8 * exp_fac(:)
      rate(:,631) = 2e-12_r8 * exp_fac(:)
      rate(:,632) = 2e-12_r8 * exp_fac(:)
      rate(:,633) = 2e-12_r8 * exp_fac(:)
      rate(:,634) = 2e-12_r8 * exp_fac(:)
      rate(:,639) = 2e-12_r8 * exp_fac(:)
      rate(:,640) = 2e-12_r8 * exp_fac(:)
      rate(:,641) = 2e-12_r8 * exp_fac(:)
      rate(:,645) = 2e-12_r8 * exp_fac(:)
      rate(:,651) = 2e-12_r8 * exp_fac(:)
      rate(:,652) = 2e-12_r8 * exp_fac(:)
      rate(:,653) = 2e-12_r8 * exp_fac(:)
      rate(:,654) = 2e-12_r8 * exp_fac(:)
      rate(:,659) = 2e-12_r8 * exp_fac(:)
      rate(:,660) = 2e-12_r8 * exp_fac(:)
      rate(:,661) = 2e-12_r8 * exp_fac(:)
      rate(:,665) = 2e-12_r8 * exp_fac(:)
      rate(:,671) = 2e-12_r8 * exp_fac(:)
      rate(:,672) = 2e-12_r8 * exp_fac(:)
      rate(:,673) = 2e-12_r8 * exp_fac(:)
      rate(:,674) = 2e-12_r8 * exp_fac(:)
      rate(:,679) = 2e-12_r8 * exp_fac(:)
      rate(:,680) = 2e-12_r8 * exp_fac(:)
      rate(:,681) = 2e-12_r8 * exp_fac(:)
      rate(:,693) = 2e-12_r8 * exp_fac(:)
      rate(:,698) = 2e-12_r8 * exp_fac(:)
      rate(:,699) = 2e-12_r8 * exp_fac(:)
      rate(:,700) = 2e-12_r8 * exp_fac(:)
      rate(:,701) = 2.9e-12_r8 * exp_fac(:)
      rate(:,702) = 2e-12_r8 * exp_fac(:)
      rate(:,706) = 2.9e-12_r8 * exp_fac(:)
      rate(:,707) = 2.9e-12_r8 * exp_fac(:)
      rate(:,709) = 2e-12_r8 * exp_fac(:)
      rate(:,714) = 2e-12_r8 * exp_fac(:)
      rate(:,715) = 2e-12_r8 * exp_fac(:)
      rate(:,716) = 2e-12_r8 * exp_fac(:)
      rate(:,719) = 2.9e-12_r8 * exp_fac(:)
      rate(:,720) = 2e-12_r8 * exp_fac(:)
      rate(:,724) = 2.9e-12_r8 * exp_fac(:)
      rate(:,725) = 2.9e-12_r8 * exp_fac(:)
      rate(:,726) = 2.9e-12_r8 * exp_fac(:)
      rate(:,728) = 2e-12_r8 * exp_fac(:)
      rate(:,733) = 2e-12_r8 * exp_fac(:)
      rate(:,734) = 2e-12_r8 * exp_fac(:)
      rate(:,735) = 2e-12_r8 * exp_fac(:)
      rate(:,738) = 2e-12_r8 * exp_fac(:)
      rate(:,743) = 2e-12_r8 * exp_fac(:)
      rate(:,744) = 2e-12_r8 * exp_fac(:)
      rate(:,745) = 2e-12_r8 * exp_fac(:)
      rate(:,749) = 2.9e-12_r8 * exp_fac(:)
      rate(:,750) = 2e-12_r8 * exp_fac(:)
      rate(:,754) = 2.9e-12_r8 * exp_fac(:)
      rate(:,312) = 5.6e-15_r8 * exp( 2300._r8 * itemp(:) )
      rate(:,313) = 2.4e+12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,314) = 2.6e-12_r8 * exp( 265._r8 * itemp(:) )
      rate(:,318) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      rate(:,323) = 1.2e-14_r8 * exp( -2630._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,326) = 7.5e-13_r8 * exp_fac(:)
      rate(:,341) = 7.5e-13_r8 * exp_fac(:)
      rate(:,356) = 7.5e-13_r8 * exp_fac(:)
      rate(:,365) = 7.5e-13_r8 * exp_fac(:)
      rate(:,369) = 8.6e-13_r8 * exp_fac(:)
      rate(:,394) = 7.5e-13_r8 * exp_fac(:)
      rate(:,409) = 7.5e-13_r8 * exp_fac(:)
      rate(:,544) = 7.5e-13_r8 * exp_fac(:)
      rate(:,548) = 7.5e-13_r8 * exp_fac(:)
      rate(:,551) = 7.5e-13_r8 * exp_fac(:)
      rate(:,564) = 7.5e-13_r8 * exp_fac(:)
      rate(:,571) = 7.5e-13_r8 * exp_fac(:)
      rate(:,577) = 7.5e-13_r8 * exp_fac(:)
      rate(:,580) = 7.5e-13_r8 * exp_fac(:)
      rate(:,852) = 7.5e-13_r8 * exp_fac(:)
      rate(:,864) = 7.5e-13_r8 * exp_fac(:)
      rate(:,879) = 7.5e-13_r8 * exp_fac(:)
      rate(:,882) = 7.5e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,327) = 2.6e-12_r8 * exp_fac(:)
      rate(:,545) = 2.6e-12_r8 * exp_fac(:)
      rate(:,550) = 2.6e-12_r8 * exp_fac(:)
      rate(:,552) = 2.6e-12_r8 * exp_fac(:)
      rate(:,565) = 2.6e-12_r8 * exp_fac(:)
      rate(:,572) = 2.6e-12_r8 * exp_fac(:)
      rate(:,578) = 2.6e-12_r8 * exp_fac(:)
      rate(:,581) = 2.6e-12_r8 * exp_fac(:)
      rate(:,853) = 2.6e-12_r8 * exp_fac(:)
      rate(:,865) = 2.6e-12_r8 * exp_fac(:)
      rate(:,880) = 2.6e-12_r8 * exp_fac(:)
      rate(:,883) = 2.6e-12_r8 * exp_fac(:)
      rate(:,328) = 6.9e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,330) = 7.2e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,331) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,332) = 1.4e-12_r8 * exp_fac(:)
      rate(:,354) = 6.5e-15_r8 * exp_fac(:)
      rate(:,333) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      rate(:,334) = 7.8e-13_r8 * exp( -1050._r8 * itemp(:) )
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,337) = 4.3e-13_r8 * exp_fac(:)
      rate(:,390) = 4.3e-13_r8 * exp_fac(:)
      rate(:,541) = 4.3e-13_r8 * exp_fac(:)
      rate(:,555) = 4.3e-13_r8 * exp_fac(:)
      rate(:,558) = 4.3e-13_r8 * exp_fac(:)
      rate(:,561) = 4.3e-13_r8 * exp_fac(:)
      rate(:,703) = 4.3e-13_r8 * exp_fac(:)
      rate(:,721) = 4.3e-13_r8 * exp_fac(:)
      rate(:,751) = 4.3e-13_r8 * exp_fac(:)
      rate(:,339) = 3.15e-14_r8 * exp( 920._r8 * itemp(:) )
      rate(:,343) = 1.6e+11_r8 * exp( -4150._r8 * itemp(:) )
      rate(:,353) = 4.6e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,355) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,359) = 9.19e-12_r8 * exp( -630._r8 * itemp(:) )
      rate(:,360) = 1.4e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,361) = 8.4e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,377) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,378) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      exp_fac(:) = exp( 1300._r8 * itemp(:) )
      rate(:,383) = 2.11e-13_r8 * exp_fac(:)
      rate(:,402) = 2.11e-13_r8 * exp_fac(:)
      rate(:,421) = 2.38e-13_r8 * exp_fac(:)
      rate(:,426) = 2.12e-13_r8 * exp_fac(:)
      rate(:,432) = 2.12e-13_r8 * exp_fac(:)
      rate(:,446) = 2.12e-13_r8 * exp_fac(:)
      rate(:,450) = 2.12e-13_r8 * exp_fac(:)
      rate(:,457) = 2.6e-13_r8 * exp_fac(:)
      rate(:,461) = 2.6e-13_r8 * exp_fac(:)
      rate(:,464) = 2.6e-13_r8 * exp_fac(:)
      rate(:,467) = 2.6e-13_r8 * exp_fac(:)
      rate(:,471) = 2.6e-13_r8 * exp_fac(:)
      rate(:,476) = 2.47e-13_r8 * exp_fac(:)
      rate(:,479) = 2.64e-13_r8 * exp_fac(:)
      rate(:,482) = 2.64e-13_r8 * exp_fac(:)
      rate(:,493) = 2.12e-13_r8 * exp_fac(:)
      rate(:,497) = 2.12e-13_r8 * exp_fac(:)
      rate(:,499) = 2.6e-13_r8 * exp_fac(:)
      rate(:,588) = 2.71e-13_r8 * exp_fac(:)
      rate(:,596) = 2.6e-13_r8 * exp_fac(:)
      rate(:,608) = 2.78e-13_r8 * exp_fac(:)
      rate(:,616) = 2.75e-13_r8 * exp_fac(:)
      rate(:,628) = 2.71e-13_r8 * exp_fac(:)
      rate(:,636) = 2.6e-13_r8 * exp_fac(:)
      rate(:,647) = 2.71e-13_r8 * exp_fac(:)
      rate(:,656) = 2.6e-13_r8 * exp_fac(:)
      rate(:,667) = 2.71e-13_r8 * exp_fac(:)
      rate(:,676) = 2.6e-13_r8 * exp_fac(:)
      rate(:,687) = 2.71e-13_r8 * exp_fac(:)
      rate(:,691) = 2.71e-13_r8 * exp_fac(:)
      rate(:,695) = 2.54e-13_r8 * exp_fac(:)
      rate(:,711) = 2.62e-13_r8 * exp_fac(:)
      rate(:,730) = 2.66e-13_r8 * exp_fac(:)
      rate(:,740) = 2.51e-13_r8 * exp_fac(:)
      rate(:,760) = 2.68e-13_r8 * exp_fac(:)
      rate(:,765) = 2.47e-13_r8 * exp_fac(:)
      rate(:,772) = 2.76e-13_r8 * exp_fac(:)
      rate(:,776) = 2.76e-13_r8 * exp_fac(:)
      rate(:,780) = 2.75e-13_r8 * exp_fac(:)
      rate(:,784) = 2.75e-13_r8 * exp_fac(:)
      rate(:,842) = 2.6e-13_r8 * exp_fac(:)
      rate(:,847) = 2.75e-13_r8 * exp_fac(:)
      rate(:,855) = 2.6e-13_r8 * exp_fac(:)
      rate(:,860) = 2.12e-13_r8 * exp_fac(:)
      rate(:,868) = 2.6e-13_r8 * exp_fac(:)
      rate(:,873) = 2.6e-13_r8 * exp_fac(:)
      rate(:,384) = 2.9e+07_r8 * exp( -5297._r8 * itemp(:) )
      rate(:,385) = 1.5e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,386) = 9.6e-12_r8 * exp_fac(:)
      rate(:,589) = 2.7e-12_r8 * exp_fac(:)
      rate(:,597) = 2.7e-12_r8 * exp_fac(:)
      rate(:,609) = 2.7e-12_r8 * exp_fac(:)
      rate(:,617) = 2.7e-12_r8 * exp_fac(:)
      rate(:,629) = 2.7e-12_r8 * exp_fac(:)
      rate(:,637) = 2.7e-12_r8 * exp_fac(:)
      rate(:,649) = 2.7e-12_r8 * exp_fac(:)
      rate(:,657) = 2.7e-12_r8 * exp_fac(:)
      rate(:,669) = 2.7e-12_r8 * exp_fac(:)
      rate(:,677) = 2.7e-12_r8 * exp_fac(:)
      rate(:,688) = 2.7e-12_r8 * exp_fac(:)
      rate(:,692) = 2.7e-12_r8 * exp_fac(:)
      rate(:,696) = 2.7e-12_r8 * exp_fac(:)
      rate(:,712) = 2.7e-12_r8 * exp_fac(:)
      rate(:,731) = 2.7e-12_r8 * exp_fac(:)
      rate(:,741) = 2.7e-12_r8 * exp_fac(:)
      rate(:,761) = 2.7e-12_r8 * exp_fac(:)
      rate(:,766) = 2.7e-12_r8 * exp_fac(:)
      rate(:,773) = 2.7e-12_r8 * exp_fac(:)
      rate(:,777) = 2.7e-12_r8 * exp_fac(:)
      rate(:,781) = 2.7e-12_r8 * exp_fac(:)
      rate(:,785) = 2.7e-12_r8 * exp_fac(:)
      rate(:,843) = 2.7e-12_r8 * exp_fac(:)
      rate(:,848) = 2.7e-12_r8 * exp_fac(:)
      rate(:,856) = 2.7e-12_r8 * exp_fac(:)
      rate(:,861) = 2.7e-12_r8 * exp_fac(:)
      rate(:,869) = 2.7e-12_r8 * exp_fac(:)
      rate(:,874) = 2.7e-12_r8 * exp_fac(:)
      rate(:,396) = 2.3e-12_r8 * exp( -170._r8 * itemp(:) )
      rate(:,404) = 2.7e-12_r8 * exp( 580._r8 * itemp(:) )
      rate(:,411) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      exp_fac(:) = exp( 450._r8 * itemp(:) )
      rate(:,414) = 1.17e-11_r8 * exp_fac(:)
      rate(:,415) = 1.17e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 390._r8 * itemp(:) )
      rate(:,416) = 2.2e-11_r8 * exp_fac(:)
      rate(:,417) = 3.5e-11_r8 * exp_fac(:)
      rate(:,487) = 2.7e-11_r8 * exp_fac(:)
      rate(:,490) = 2.08e-11_r8 * exp_fac(:)
      rate(:,768) = 2.7e-11_r8 * exp_fac(:)
      rate(:,863) = 2.7e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( 410._r8 * itemp(:) )
      rate(:,419) = 9.85e-12_r8 * exp_fac(:)
      rate(:,603) = 1.34e-11_r8 * exp_fac(:)
      rate(:,845) = 1.34e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -400._r8 * itemp(:) )
      rate(:,420) = 4.43e-11_r8 * exp_fac(:)
      rate(:,422) = 4.43e-11_r8 * exp_fac(:)
      rate(:,423) = 3.22e-11_r8 * exp_fac(:)
      rate(:,427) = 1.04e+11_r8 * exp( -9746._r8 * itemp(:) )
      rate(:,428) = 2.24e+15_r8 * exp( -10865._r8 * itemp(:) )
      rate(:,429) = 2.22e+15_r8 * exp( -10355._r8 * itemp(:) )
      rate(:,433) = 1.88e+11_r8 * exp( -9752._r8 * itemp(:) )
      rate(:,434) = 2.49e+15_r8 * exp( -11112._r8 * itemp(:) )
      rate(:,435) = 2.49e+15_r8 * exp( -10890._r8 * itemp(:) )
      rate(:,447) = 1.83e+14_r8 * exp( -8930._r8 * itemp(:) )
      rate(:,451) = 2.08e+14_r8 * exp( -9400._r8 * itemp(:) )
      exp_fac(:) = exp( -10000._r8 * itemp(:) )
      rate(:,458) = 1.256e+13_r8 * exp_fac(:)
      rate(:,462) = 1.875e+13_r8 * exp_fac(:)
      rate(:,465) = 1.875e+13_r8 * exp_fac(:)
      rate(:,468) = 5.092e+12_r8 * exp_fac(:)
      rate(:,480) = 8.72e+12_r8 * exp_fac(:)
      rate(:,483) = 6.55e+12_r8 * exp_fac(:)
      exp_fac(:) = exp( -450._r8 * itemp(:) )
      rate(:,473) = 2.95e-12_r8 * exp_fac(:)
      rate(:,764) = 2.95e-12_r8 * exp_fac(:)
      rate(:,859) = 2.95e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1995._r8 * itemp(:) )
      rate(:,486) = 1.03e-14_r8 * exp_fac(:)
      rate(:,862) = 1.03e-14_r8 * exp_fac(:)
      rate(:,494) = 1.79e+14_r8 * exp( -8830._r8 * itemp(:) )
      rate(:,498) = 1.75e+14_r8 * exp( -9054._r8 * itemp(:) )
      rate(:,500) = 1e+07_r8 * exp( -5000._r8 * itemp(:) )
      exp_fac(:) = exp( -193._r8 * itemp(:) )
      rate(:,543) = 2.3e-12_r8 * exp_fac(:)
      rate(:,851) = 2.3e-12_r8 * exp_fac(:)
      rate(:,547) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,566) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      exp_fac(:) = exp( 352._r8 * itemp(:) )
      rate(:,574) = 1.7e-12_r8 * exp_fac(:)
      rate(:,878) = 1.7e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 490._r8 * itemp(:) )
      rate(:,584) = 1.2e-12_r8 * exp_fac(:)
      rate(:,841) = 1.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -640._r8 * itemp(:) )
      rate(:,602) = 8.05e-16_r8 * exp_fac(:)
      rate(:,844) = 8.05e-16_r8 * exp_fac(:)
      exp_fac(:) = exp( -770._r8 * itemp(:) )
      rate(:,662) = 2.8e-15_r8 * exp_fac(:)
      rate(:,870) = 2.8e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 470._r8 * itemp(:) )
      rate(:,663) = 3.41e-11_r8 * exp_fac(:)
      rate(:,871) = 3.41e-11_r8 * exp_fac(:)
      exp_fac(:) = exp( -520._r8 * itemp(:) )
      rate(:,682) = 2.65e-15_r8 * exp_fac(:)
      rate(:,875) = 2.65e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 600._r8 * itemp(:) )
      rate(:,717) = 5.2e-12_r8 * exp_fac(:)
      rate(:,736) = 5.2e-12_r8 * exp_fac(:)
      rate(:,756) = 5.2e-12_r8 * exp_fac(:)
      rate(:,793) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,794) = 1.1e-11_r8 * exp( -280._r8 * itemp(:) )
      rate(:,795) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,796) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,804) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,807) = 2.6e-11_r8 * exp( 330._r8 * itemp(:) )
      rate(:,810) = 1.7e-12_r8 * exp( -710._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( rate(:,181), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,191), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,201), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,209), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,212), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,213), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,214), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,232), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,252), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,263), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.1e-33_r8 * itemp(:)**1.5_r8
      kinf(:) = 9.8e-15_r8 * itemp(:)**(-4.6_r8)
      call jpl( rate(:,309), m, 0.8_r8, ko, kinf, n )

      ko(:) = 5.2e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,320), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.5e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2._r8)
      call jpl( rate(:,321), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.6e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,322), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.6e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,350), m, 0.48_r8, ko, kinf, n )

      ko(:) = 7.3e-29_r8 * itemp(:)**4.1_r8
      kinf(:) = 9.5e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,351), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,372), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,398), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,406), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,557), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,560), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,563), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,570), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,684), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,685), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.7e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.3e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,686), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.9e-31_r8 * itemp(:)**4.1_r8
      kinf(:) = 1.7e-12_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,798), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,178) = 6.9e-12_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,170) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      rate(:n,173) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,182) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,183) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,184) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,187) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,188) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,189) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,194) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,198) = 3.3e-12_r8 * exp( -3150._r8 * itemp(:) )
      rate(:n,199) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,207) = 3.44e-12_r8 * exp( 260._r8 * itemp(:) )
      rate(:n,208) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 5.3e-32_r8 * itemp(:)**1.8_r8
      kinf(:) = 9.5e-11_r8 * itemp(:)**(-0.4_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,181) = wrk(:)



























      end subroutine setrxt_hrates

      end module mo_setrxt
