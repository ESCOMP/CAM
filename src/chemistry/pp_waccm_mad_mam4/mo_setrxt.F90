
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

      rate(:,104) = 9.6e-10_r8
      rate(:,105) = 1.3e-09_r8
      rate(:,106) = 2e-29_r8
      rate(:,107) = 1e-27_r8
      rate(:,108) = 1.6e-09_r8
      rate(:,109) = 6e-12_r8
      rate(:,110) = 2.9e-12_r8
      rate(:,111) = 2.9e-11_r8
      rate(:,112) = 2e-10_r8
      rate(:,113) = 1e-10_r8
      rate(:,114) = 1e-10_r8
      rate(:,115) = 1e-11_r8
      rate(:,116) = 1.7e-10_r8
      rate(:,117) = 1e-28_r8
      rate(:,118) = 1e-28_r8
      rate(:,119) = 4e-11_r8
      rate(:,120) = 4e-11_r8
      rate(:,121) = 3.5e-12_r8
      rate(:,122) = 3.5e-12_r8
      rate(:,123) = 3.51e-10_r8
      rate(:,124) = 1.1e-10_r8
      rate(:,125) = 6e-15_r8
      rate(:,126) = 1e-10_r8
      rate(:,127) = 1e-10_r8
      rate(:,128) = 2.2e-10_r8
      rate(:,129) = 1.2e-09_r8
      rate(:,130) = 1.4e-10_r8
      rate(:,131) = 1.3e-10_r8
      rate(:,137) = 1.5e-06_r8
      rate(:,138) = 2e-09_r8
      rate(:,139) = 1e-09_r8
      rate(:,140) = 3.6e-06_r8
      rate(:,141) = 4e-12_r8
      rate(:,142) = 1e-09_r8
      rate(:,143) = 5e-06_r8
      rate(:,144) = 7e-12_r8
      rate(:,284) = 1e-10_r8
      rate(:,285) = 1e-10_r8
      rate(:,286) = 3e-10_r8
      rate(:,287) = 1.6e-28_r8
      rate(:,288) = 1.4e-09_r8
      rate(:,289) = 1.6e-09_r8
      rate(:,290) = 2e-13_r8
      rate(:,291) = 1.2e-10_r8
      rate(:,292) = 7e-10_r8
      rate(:,293) = 1.6e-28_r8
      rate(:,294) = 1.6e-09_r8
      rate(:,295) = 1.6e-28_r8
      rate(:,296) = 7e-10_r8
      rate(:,297) = 1e-12_r8
      rate(:,298) = 7.6e-10_r8
      rate(:,299) = 1.45e-26_r8
      rate(:,300) = 5e-12_r8
      rate(:,301) = 1e-13_r8
      rate(:,302) = 2e-06_r8
      rate(:,303) = 2e-06_r8
      rate(:,304) = 7e-11_r8
      rate(:,305) = 1.5e-06_r8
      rate(:,306) = 1e-09_r8
      rate(:,307) = 1.5e-06_r8
      rate(:,308) = 7e-12_r8
      rate(:,309) = 5e-10_r8
      rate(:,310) = 1e-10_r8
      rate(:,311) = 1e-09_r8
      rate(:,312) = 1e-09_r8
      rate(:,313) = 1e-10_r8
      rate(:,314) = 1e-10_r8
      rate(:,315) = 9.9e-30_r8
      rate(:,316) = 1.4e-09_r8
      rate(:,317) = 1.6e-09_r8
      rate(:,318) = 2.9e-09_r8
      rate(:,319) = 7e-10_r8
      rate(:,320) = 2e-10_r8
      rate(:,321) = 3.4e-31_r8
      rate(:,322) = 7.8e-10_r8
      rate(:,323) = 1.5e-10_r8
      rate(:,324) = 1.5e-10_r8
      rate(:,325) = 2e-06_r8
      rate(:,326) = 9e-10_r8
      rate(:,327) = 2.4e-10_r8
      rate(:,328) = 2.8e-28_r8
      rate(:,329) = 5.5e-10_r8
      rate(:,330) = 8.4e-10_r8
      rate(:,331) = 1e-10_r8
      rate(:,332) = 1e-10_r8
      rate(:,333) = 2.5e-10_r8
      rate(:,334) = 4.3e-10_r8
      rate(:,335) = 4e-10_r8
      rate(:,336) = 1.7e-09_r8
      rate(:,337) = 3e-10_r8
      rate(:,338) = 1.5e-10_r8
      rate(:,340) = 1e-10_r8
      rate(:,341) = 1e-10_r8
      rate(:,342) = 7.6e-28_r8
      rate(:,343) = 1.4e-09_r8
      rate(:,344) = 1e-09_r8
      rate(:,345) = 1.1e-09_r8
      rate(:,346) = 2e-10_r8
      rate(:,347) = 9e-10_r8
      rate(:,349) = 1e-10_r8
      rate(:,350) = 1e-10_r8
      rate(:,351) = 2e-28_r8
      rate(:,352) = 5.8e-10_r8
      rate(:,353) = 3.2e-11_r8
      rate(:,354) = 6e-13_r8
      rate(:,355) = 2e-09_r8
      rate(:,356) = 3.6e-09_r8
      rate(:,357) = 5e-13_r8
      rate(:,358) = 1e-09_r8
      rate(:,359) = 1.9e-10_r8
      rate(:,360) = 3e-10_r8
      rate(:,361) = 2.9e-31_r8
      rate(:,362) = 8e-10_r8
      rate(:,386) = 0.000258_r8
      rate(:,387) = 0.085_r8
      rate(:,388) = 1.2e-10_r8
      rate(:,393) = 1.2e-10_r8
      rate(:,394) = 1e-20_r8
      rate(:,395) = 1.3e-16_r8
      rate(:,397) = 4.2e-13_r8
      rate(:,399) = 8e-14_r8
      rate(:,400) = 3.9e-17_r8
      rate(:,407) = 6.9e-12_r8
      rate(:,408) = 7.2e-11_r8
      rate(:,409) = 1.6e-12_r8
      rate(:,415) = 1.8e-12_r8
      rate(:,419) = 1.8e-12_r8
      rate(:,423) = 7e-13_r8
      rate(:,424) = 5e-12_r8
      rate(:,433) = 3.5e-12_r8
      rate(:,435) = 1e-11_r8
      rate(:,436) = 2.2e-11_r8
      rate(:,437) = 5e-11_r8
      rate(:,472) = 1.7e-13_r8
      rate(:,474) = 2.607e-10_r8
      rate(:,475) = 9.75e-11_r8
      rate(:,476) = 2.07e-10_r8
      rate(:,477) = 2.088e-10_r8
      rate(:,478) = 1.17e-10_r8
      rate(:,479) = 4.644e-11_r8
      rate(:,480) = 1.204e-10_r8
      rate(:,481) = 9.9e-11_r8
      rate(:,482) = 3.3e-12_r8
      rate(:,501) = 4.5e-11_r8
      rate(:,502) = 4.62e-10_r8
      rate(:,503) = 1.2e-10_r8
      rate(:,504) = 9e-11_r8
      rate(:,505) = 3e-11_r8
      rate(:,510) = 2.14e-11_r8
      rate(:,511) = 1.9e-10_r8
      rate(:,524) = 2.57e-10_r8
      rate(:,525) = 1.8e-10_r8
      rate(:,526) = 1.794e-10_r8
      rate(:,527) = 1.3e-10_r8
      rate(:,528) = 7.65e-11_r8
      rate(:,537) = 1.31e-10_r8
      rate(:,538) = 3.5e-11_r8
      rate(:,539) = 9e-12_r8
      rate(:,543) = 2.3e-12_r8
      rate(:,544) = 1.2e-11_r8
      rate(:,545) = 5.7e-11_r8
      rate(:,546) = 2.8e-11_r8
      rate(:,547) = 6.6e-11_r8
      rate(:,548) = 1.4e-11_r8
      rate(:,551) = 1.9e-12_r8
      rate(:,582) = 6e-11_r8
      rate(:,585) = 1e-12_r8
      rate(:,586) = 4e-10_r8
      rate(:,587) = 2e-10_r8
      rate(:,588) = 1e-10_r8
      rate(:,589) = 5e-16_r8
      rate(:,590) = 4.4e-10_r8
      rate(:,591) = 9e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,339) = 1.8e-11_r8 * exp( 390._r8 * itemp(:) )
      rate(:,389) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,390) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,391) = 2.64e-11_r8 * exp_fac(:)
      rate(:,392) = 6.6e-12_r8 * exp_fac(:)
      rate(:,396) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,398) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,401) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,402) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,405) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,406) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,411) = 3e-11_r8 * exp_fac(:)
      rate(:,499) = 5.5e-12_r8 * exp_fac(:)
      rate(:,534) = 3.8e-12_r8 * exp_fac(:)
      rate(:,412) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,413) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,414) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,416) = 4.8e-11_r8 * exp_fac(:)
      rate(:,497) = 1.7e-11_r8 * exp_fac(:)
      rate(:,417) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,418) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,422) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,425) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,426) = 2.9e-12_r8 * exp_fac(:)
      rate(:,427) = 1.45e-12_r8 * exp_fac(:)
      rate(:,428) = 1.45e-12_r8 * exp_fac(:)
      rate(:,429) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,430) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,431) = 1.2e-13_r8 * exp_fac(:)
      rate(:,457) = 3e-11_r8 * exp_fac(:)
      rate(:,434) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,438) = 3.3e-12_r8 * exp_fac(:)
      rate(:,453) = 1.4e-11_r8 * exp_fac(:)
      rate(:,467) = 7.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,439) = 3e-12_r8 * exp_fac(:)
      rate(:,498) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,441) = 7.26e-11_r8 * exp_fac(:)
      rate(:,442) = 4.64e-11_r8 * exp_fac(:)
      rate(:,449) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,450) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,451) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,452) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,454) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,455) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,456) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,458) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,459) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,460) = 2.6e-12_r8 * exp_fac(:)
      rate(:,461) = 6.4e-12_r8 * exp_fac(:)
      rate(:,491) = 4.1e-13_r8 * exp_fac(:)
      rate(:,462) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,464) = 3.6e-12_r8 * exp_fac(:)
      rate(:,513) = 2e-12_r8 * exp_fac(:)
      rate(:,465) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,466) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,468) = 6e-13_r8 * exp_fac(:)
      rate(:,488) = 1.5e-12_r8 * exp_fac(:)
      rate(:,496) = 1.9e-11_r8 * exp_fac(:)
      rate(:,469) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,470) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,471) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,473) = 3e-12_r8 * exp_fac(:)
      rate(:,507) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,485) = 1.7e-11_r8 * exp_fac(:)
      rate(:,512) = 6.3e-12_r8 * exp_fac(:)
      rate(:,486) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,487) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,489) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,490) = 2.3e-12_r8 * exp_fac(:)
      rate(:,493) = 8.8e-12_r8 * exp_fac(:)
      rate(:,492) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,495) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,500) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,506) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,508) = 1.4e-11_r8 * exp_fac(:)
      rate(:,510) = 2.14e-11_r8 * exp_fac(:)
      rate(:,511) = 1.9e-10_r8 * exp_fac(:)
      rate(:,524) = 2.57e-10_r8 * exp_fac(:)
      rate(:,525) = 1.8e-10_r8 * exp_fac(:)
      rate(:,526) = 1.794e-10_r8 * exp_fac(:)
      rate(:,527) = 1.3e-10_r8 * exp_fac(:)
      rate(:,528) = 7.65e-11_r8 * exp_fac(:)
      rate(:,537) = 1.31e-10_r8 * exp_fac(:)
      rate(:,538) = 3.5e-11_r8 * exp_fac(:)
      rate(:,539) = 9e-12_r8 * exp_fac(:)
      rate(:,543) = 2.3e-12_r8 * exp_fac(:)
      rate(:,544) = 1.2e-11_r8 * exp_fac(:)
      rate(:,545) = 5.7e-11_r8 * exp_fac(:)
      rate(:,546) = 2.8e-11_r8 * exp_fac(:)
      rate(:,547) = 6.6e-11_r8 * exp_fac(:)
      rate(:,548) = 1.4e-11_r8 * exp_fac(:)
      rate(:,551) = 1.9e-12_r8 * exp_fac(:)
      rate(:,582) = 6e-11_r8 * exp_fac(:)
      rate(:,585) = 1e-12_r8 * exp_fac(:)
      rate(:,586) = 4e-10_r8 * exp_fac(:)
      rate(:,587) = 2e-10_r8 * exp_fac(:)
      rate(:,588) = 1e-10_r8 * exp_fac(:)
      rate(:,589) = 5e-16_r8 * exp_fac(:)
      rate(:,590) = 4.4e-10_r8 * exp_fac(:)
      rate(:,591) = 9e-10_r8 * exp_fac(:)
      rate(:,509) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,514) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,515) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,516) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,517) = 2.03e-11_r8 * exp_fac(:)
      rate(:,550) = 3.4e-12_r8 * exp_fac(:)
      rate(:,518) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,519) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,520) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,521) = 1.25e-12_r8 * exp_fac(:)
      rate(:,530) = 3.4e-11_r8 * exp_fac(:)
      rate(:,522) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,523) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,529) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,531) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,532) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,533) = 2.8e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,535) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,541) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,542) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,549) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,552) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,555) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,556) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 7e-31_r8 * itemp(:)**2.6_r8
      kinf(:) = 3.6e-11_r8 * itemp(:)**0.1_r8
      call jpl( rate(:,348), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,410), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,420), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,432), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,440), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,443), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,444), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,445), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,463), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,483), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,494), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,536), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,394) = 1e-20_r8
      rate(:n,395) = 1.3e-16_r8
      rate(:n,399) = 8e-14_r8
      rate(:n,400) = 3.9e-17_r8
      rate(:n,407) = 6.9e-12_r8
      rate(:n,423) = 7e-13_r8
      rate(:n,424) = 5e-12_r8
      rate(:n,582) = 6e-11_r8
      rate(:n,585) = 1e-12_r8
      rate(:n,586) = 4e-10_r8
      rate(:n,587) = 2e-10_r8
      rate(:n,588) = 1e-10_r8
      rate(:n,590) = 4.4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,390) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,391) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,392) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,396) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,398) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,401) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,402) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,411) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,412) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,413) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,416) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,417) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,418) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,425) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,429) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,430) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,438) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,439) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)


      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,410) = wrk(:)











      end subroutine setrxt_hrates

      end module mo_setrxt
