
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

      rate(:,105) = 9.6e-10_r8
      rate(:,106) = 1.3e-09_r8
      rate(:,107) = 2e-29_r8
      rate(:,108) = 1e-27_r8
      rate(:,109) = 1.6e-09_r8
      rate(:,110) = 6e-12_r8
      rate(:,111) = 2.9e-12_r8
      rate(:,112) = 2.9e-11_r8
      rate(:,113) = 2e-10_r8
      rate(:,114) = 1e-10_r8
      rate(:,115) = 1e-10_r8
      rate(:,116) = 1e-11_r8
      rate(:,117) = 1.7e-10_r8
      rate(:,118) = 1e-28_r8
      rate(:,119) = 1e-28_r8
      rate(:,120) = 4e-11_r8
      rate(:,121) = 4e-11_r8
      rate(:,122) = 3.5e-12_r8
      rate(:,123) = 3.5e-12_r8
      rate(:,124) = 3.51e-10_r8
      rate(:,125) = 1.1e-10_r8
      rate(:,126) = 6e-15_r8
      rate(:,127) = 1e-10_r8
      rate(:,128) = 1e-10_r8
      rate(:,129) = 2.2e-10_r8
      rate(:,130) = 1.2e-09_r8
      rate(:,131) = 1.4e-10_r8
      rate(:,132) = 1.3e-10_r8
      rate(:,138) = 1.5e-06_r8
      rate(:,139) = 2e-09_r8
      rate(:,140) = 1e-09_r8
      rate(:,141) = 3.6e-06_r8
      rate(:,142) = 4e-12_r8
      rate(:,143) = 1e-09_r8
      rate(:,144) = 5e-06_r8
      rate(:,145) = 7e-12_r8
      rate(:,285) = 1e-10_r8
      rate(:,286) = 1e-10_r8
      rate(:,287) = 3e-10_r8
      rate(:,288) = 1.6e-28_r8
      rate(:,289) = 1.4e-09_r8
      rate(:,290) = 1.6e-09_r8
      rate(:,291) = 2e-13_r8
      rate(:,292) = 1.2e-10_r8
      rate(:,293) = 7e-10_r8
      rate(:,294) = 1.6e-28_r8
      rate(:,295) = 1.6e-09_r8
      rate(:,296) = 1.6e-28_r8
      rate(:,297) = 7e-10_r8
      rate(:,298) = 1e-12_r8
      rate(:,299) = 7.6e-10_r8
      rate(:,300) = 1.45e-26_r8
      rate(:,301) = 5e-12_r8
      rate(:,302) = 1e-13_r8
      rate(:,303) = 2e-06_r8
      rate(:,304) = 2e-06_r8
      rate(:,305) = 7e-11_r8
      rate(:,306) = 1.5e-06_r8
      rate(:,307) = 1e-09_r8
      rate(:,308) = 1.5e-06_r8
      rate(:,309) = 7e-12_r8
      rate(:,310) = 5e-10_r8
      rate(:,311) = 1e-10_r8
      rate(:,312) = 1e-09_r8
      rate(:,313) = 1e-09_r8
      rate(:,314) = 1e-10_r8
      rate(:,315) = 1e-10_r8
      rate(:,316) = 9.9e-30_r8
      rate(:,317) = 1.4e-09_r8
      rate(:,318) = 1.6e-09_r8
      rate(:,319) = 2.9e-09_r8
      rate(:,320) = 7e-10_r8
      rate(:,321) = 2e-10_r8
      rate(:,322) = 3.4e-31_r8
      rate(:,323) = 7.8e-10_r8
      rate(:,324) = 1.5e-10_r8
      rate(:,325) = 1.5e-10_r8
      rate(:,326) = 2e-06_r8
      rate(:,327) = 9e-10_r8
      rate(:,328) = 2.4e-10_r8
      rate(:,329) = 2.8e-28_r8
      rate(:,330) = 5.5e-10_r8
      rate(:,331) = 8.4e-10_r8
      rate(:,332) = 1e-10_r8
      rate(:,333) = 1e-10_r8
      rate(:,334) = 2.5e-10_r8
      rate(:,335) = 4.3e-10_r8
      rate(:,336) = 4e-10_r8
      rate(:,337) = 1.7e-09_r8
      rate(:,338) = 3e-10_r8
      rate(:,339) = 1.5e-10_r8
      rate(:,341) = 1e-10_r8
      rate(:,342) = 1e-10_r8
      rate(:,343) = 7.6e-28_r8
      rate(:,344) = 1.4e-09_r8
      rate(:,345) = 1e-09_r8
      rate(:,346) = 1.1e-09_r8
      rate(:,347) = 2e-10_r8
      rate(:,348) = 9e-10_r8
      rate(:,350) = 1e-10_r8
      rate(:,351) = 1e-10_r8
      rate(:,352) = 2e-28_r8
      rate(:,353) = 5.8e-10_r8
      rate(:,354) = 3.2e-11_r8
      rate(:,355) = 6e-13_r8
      rate(:,356) = 2e-09_r8
      rate(:,357) = 3.6e-09_r8
      rate(:,358) = 5e-13_r8
      rate(:,359) = 1e-09_r8
      rate(:,360) = 1.9e-10_r8
      rate(:,361) = 3e-10_r8
      rate(:,362) = 2.9e-31_r8
      rate(:,363) = 8e-10_r8
      rate(:,387) = 0.000258_r8
      rate(:,388) = 0.085_r8
      rate(:,389) = 1.2e-10_r8
      rate(:,394) = 1.2e-10_r8
      rate(:,395) = 1e-20_r8
      rate(:,396) = 1.3e-16_r8
      rate(:,398) = 4.2e-13_r8
      rate(:,400) = 8e-14_r8
      rate(:,401) = 3.9e-17_r8
      rate(:,408) = 6.9e-12_r8
      rate(:,409) = 7.2e-11_r8
      rate(:,410) = 1.6e-12_r8
      rate(:,416) = 1.8e-12_r8
      rate(:,420) = 1.8e-12_r8
      rate(:,424) = 7e-13_r8
      rate(:,425) = 5e-12_r8
      rate(:,434) = 3.5e-12_r8
      rate(:,436) = 1e-11_r8
      rate(:,437) = 2.2e-11_r8
      rate(:,438) = 5e-11_r8
      rate(:,473) = 1.7e-13_r8
      rate(:,475) = 2.607e-10_r8
      rate(:,476) = 9.75e-11_r8
      rate(:,477) = 2.07e-10_r8
      rate(:,478) = 2.088e-10_r8
      rate(:,479) = 1.17e-10_r8
      rate(:,480) = 4.644e-11_r8
      rate(:,481) = 1.204e-10_r8
      rate(:,482) = 9.9e-11_r8
      rate(:,483) = 3.3e-12_r8
      rate(:,502) = 4.5e-11_r8
      rate(:,503) = 4.62e-10_r8
      rate(:,504) = 1.2e-10_r8
      rate(:,505) = 9e-11_r8
      rate(:,506) = 3e-11_r8
      rate(:,511) = 2.14e-11_r8
      rate(:,512) = 1.9e-10_r8
      rate(:,525) = 2.57e-10_r8
      rate(:,526) = 1.8e-10_r8
      rate(:,527) = 1.794e-10_r8
      rate(:,528) = 1.3e-10_r8
      rate(:,529) = 7.65e-11_r8
      rate(:,538) = 1.31e-10_r8
      rate(:,539) = 3.5e-11_r8
      rate(:,540) = 9e-12_r8
      rate(:,544) = 2.3e-12_r8
      rate(:,545) = 1.2e-11_r8
      rate(:,546) = 5.7e-11_r8
      rate(:,547) = 2.8e-11_r8
      rate(:,548) = 6.6e-11_r8
      rate(:,549) = 1.4e-11_r8
      rate(:,552) = 1.9e-12_r8
      rate(:,580) = 0.047_r8
      rate(:,581) = 7.7e-05_r8
      rate(:,582) = 0.171_r8
      rate(:,586) = 6e-11_r8
      rate(:,589) = 1e-12_r8
      rate(:,590) = 4e-10_r8
      rate(:,591) = 2e-10_r8
      rate(:,592) = 1e-10_r8
      rate(:,593) = 5e-16_r8
      rate(:,594) = 4.4e-10_r8
      rate(:,595) = 9e-10_r8
      rate(:,597) = 1.3e-10_r8
      rate(:,600) = 8e-10_r8
      rate(:,601) = 5e-12_r8
      rate(:,602) = 7e-10_r8
      rate(:,605) = 4.8e-10_r8
      rate(:,606) = 1e-10_r8
      rate(:,607) = 4e-10_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      rate(:,340) = 1.8e-11_r8 * exp( 390._r8 * itemp(:) )
      rate(:,390) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      rate(:,391) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,392) = 2.64e-11_r8 * exp_fac(:)
      rate(:,393) = 6.6e-12_r8 * exp_fac(:)
      rate(:,397) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,399) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,402) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,403) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:,406) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:) )
      rate(:,407) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,412) = 3e-11_r8 * exp_fac(:)
      rate(:,500) = 5.5e-12_r8 * exp_fac(:)
      rate(:,535) = 3.8e-12_r8 * exp_fac(:)
      rate(:,413) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:,414) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:,415) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:) )
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,417) = 4.8e-11_r8 * exp_fac(:)
      rate(:,498) = 1.7e-11_r8 * exp_fac(:)
      rate(:,418) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:,419) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:,423) = 1.3e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,426) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,427) = 2.9e-12_r8 * exp_fac(:)
      rate(:,428) = 1.45e-12_r8 * exp_fac(:)
      rate(:,429) = 1.45e-12_r8 * exp_fac(:)
      rate(:,430) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,431) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,432) = 1.2e-13_r8 * exp_fac(:)
      rate(:,458) = 3e-11_r8 * exp_fac(:)
      rate(:,435) = 1.5e-11_r8 * exp( 170._r8 * itemp(:) )
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,439) = 3.3e-12_r8 * exp_fac(:)
      rate(:,454) = 1.4e-11_r8 * exp_fac(:)
      rate(:,468) = 7.4e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,440) = 3e-12_r8 * exp_fac(:)
      rate(:,499) = 5.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,442) = 7.26e-11_r8 * exp_fac(:)
      rate(:,443) = 4.64e-11_r8 * exp_fac(:)
      rate(:,450) = 8.1e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,451) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,452) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,453) = 1.1e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,455) = 3.6e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,456) = 2.3e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,457) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,459) = 1e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,460) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:) )
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,461) = 2.6e-12_r8 * exp_fac(:)
      rate(:,462) = 6.4e-12_r8 * exp_fac(:)
      rate(:,492) = 4.1e-13_r8 * exp_fac(:)
      rate(:,463) = 6.5e-12_r8 * exp( 135._r8 * itemp(:) )
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,465) = 3.6e-12_r8 * exp_fac(:)
      rate(:,514) = 2e-12_r8 * exp_fac(:)
      rate(:,466) = 1.2e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,467) = 2.8e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,469) = 6e-13_r8 * exp_fac(:)
      rate(:,489) = 1.5e-12_r8 * exp_fac(:)
      rate(:,497) = 1.9e-11_r8 * exp_fac(:)
      rate(:,470) = 1e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,471) = 1.8e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,472) = 3.4e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,474) = 3e-12_r8 * exp_fac(:)
      rate(:,508) = 1.4e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,486) = 1.7e-11_r8 * exp_fac(:)
      rate(:,513) = 6.3e-12_r8 * exp_fac(:)
      rate(:,487) = 4.8e-12_r8 * exp( -310._r8 * itemp(:) )
      rate(:,488) = 1.6e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,490) = 9.5e-13_r8 * exp( 550._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,491) = 2.3e-12_r8 * exp_fac(:)
      rate(:,494) = 8.8e-12_r8 * exp_fac(:)
      rate(:,493) = 4.5e-12_r8 * exp( 460._r8 * itemp(:) )
      rate(:,496) = 1.9e-11_r8 * exp( 215._r8 * itemp(:) )
      rate(:,501) = 1.2e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,507) = 1.6e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,509) = 1.4e-11_r8 * exp_fac(:)
      rate(:,511) = 2.14e-11_r8 * exp_fac(:)
      rate(:,512) = 1.9e-10_r8 * exp_fac(:)
      rate(:,525) = 2.57e-10_r8 * exp_fac(:)
      rate(:,526) = 1.8e-10_r8 * exp_fac(:)
      rate(:,527) = 1.794e-10_r8 * exp_fac(:)
      rate(:,528) = 1.3e-10_r8 * exp_fac(:)
      rate(:,529) = 7.65e-11_r8 * exp_fac(:)
      rate(:,538) = 1.31e-10_r8 * exp_fac(:)
      rate(:,539) = 3.5e-11_r8 * exp_fac(:)
      rate(:,540) = 9e-12_r8 * exp_fac(:)
      rate(:,544) = 2.3e-12_r8 * exp_fac(:)
      rate(:,545) = 1.2e-11_r8 * exp_fac(:)
      rate(:,546) = 5.7e-11_r8 * exp_fac(:)
      rate(:,547) = 2.8e-11_r8 * exp_fac(:)
      rate(:,548) = 6.6e-11_r8 * exp_fac(:)
      rate(:,549) = 1.4e-11_r8 * exp_fac(:)
      rate(:,552) = 1.9e-12_r8 * exp_fac(:)
      rate(:,580) = 0.047_r8 * exp_fac(:)
      rate(:,581) = 7.7e-05_r8 * exp_fac(:)
      rate(:,582) = 0.171_r8 * exp_fac(:)
      rate(:,586) = 6e-11_r8 * exp_fac(:)
      rate(:,589) = 1e-12_r8 * exp_fac(:)
      rate(:,590) = 4e-10_r8 * exp_fac(:)
      rate(:,591) = 2e-10_r8 * exp_fac(:)
      rate(:,592) = 1e-10_r8 * exp_fac(:)
      rate(:,593) = 5e-16_r8 * exp_fac(:)
      rate(:,594) = 4.4e-10_r8 * exp_fac(:)
      rate(:,595) = 9e-10_r8 * exp_fac(:)
      rate(:,597) = 1.3e-10_r8 * exp_fac(:)
      rate(:,600) = 8e-10_r8 * exp_fac(:)
      rate(:,601) = 5e-12_r8 * exp_fac(:)
      rate(:,602) = 7e-10_r8 * exp_fac(:)
      rate(:,605) = 4.8e-10_r8 * exp_fac(:)
      rate(:,606) = 1e-10_r8 * exp_fac(:)
      rate(:,607) = 4e-10_r8 * exp_fac(:)
      rate(:,510) = 6e-12_r8 * exp( 400._r8 * itemp(:) )
      rate(:,515) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      rate(:,516) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,517) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:) )
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,518) = 2.03e-11_r8 * exp_fac(:)
      rate(:,551) = 3.4e-12_r8 * exp_fac(:)
      rate(:,519) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      rate(:,520) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,521) = 9e-13_r8 * exp( -360._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,522) = 1.25e-12_r8 * exp_fac(:)
      rate(:,531) = 3.4e-11_r8 * exp_fac(:)
      rate(:,523) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,524) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,530) = 6e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,532) = 5.5e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,533) = 4.1e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,534) = 2.8e-12_r8 * exp( 300._r8 * itemp(:) )
      rate(:,536) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,542) = 2.1e-11_r8 * exp( -2200._r8 * itemp(:) )
      rate(:,543) = 7.2e-14_r8 * exp( -1070._r8 * itemp(:) )
      rate(:,550) = 1.6e-13_r8 * exp( -2280._r8 * itemp(:) )
      rate(:,553) = 2.7e-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,556) = 1.9e-13_r8 * exp( 520._r8 * itemp(:) )
      rate(:,557) = 9.6e-12_r8 * exp( -234._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 7e-31_r8 * itemp(:)**2.6_r8
      kinf(:) = 3.6e-11_r8 * itemp(:)**0.1_r8
      call jpl( rate(:,349), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,411), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.9e-31_r8 * itemp(:)**1._r8
      kinf(:) = 2.6e-11_r8
      call jpl( rate(:,421), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.5e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,433), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3e-11_r8
      call jpl( rate(:,441), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,444), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.4e-30_r8 * itemp(:)**3._r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,445), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-30_r8 * itemp(:)**3._r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,446), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.8e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,464), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.9e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,484), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.2e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,495), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.9e-33_r8 * itemp(:)**1._r8
      kinf(:) = 1.1e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,537), m, 0.6_r8, ko, kinf, n )

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

      rate(:n,395) = 1e-20_r8
      rate(:n,396) = 1.3e-16_r8
      rate(:n,400) = 8e-14_r8
      rate(:n,401) = 3.9e-17_r8
      rate(:n,408) = 6.9e-12_r8
      rate(:n,424) = 7e-13_r8
      rate(:n,425) = 5e-12_r8
      rate(:n,580) = 0.047_r8
      rate(:n,581) = 7.7e-05_r8
      rate(:n,582) = 0.171_r8
      rate(:n,586) = 6e-11_r8
      rate(:n,589) = 1e-12_r8
      rate(:n,590) = 4e-10_r8
      rate(:n,591) = 2e-10_r8
      rate(:n,592) = 1e-10_r8
      rate(:n,594) = 4.4e-10_r8
      rate(:n,597) = 1.3e-10_r8
      rate(:n,600) = 8e-10_r8
      rate(:n,601) = 5e-12_r8
      rate(:n,602) = 7e-10_r8
      rate(:n,605) = 4.8e-10_r8
      rate(:n,606) = 1e-10_r8
      rate(:n,607) = 4e-10_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,391) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,392) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,393) = 6.6e-12_r8 * exp_fac(:)
      rate(:n,397) = 3.6e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,399) = 1.8e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,402) = 3.5e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,403) = 8e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,412) = 3e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,413) = 1e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,414) = 1.4e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,417) = 4.8e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,418) = 1.8e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,419) = 1.7e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,426) = 2.1e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,430) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,431) = 5.1e-12_r8 * exp( 210._r8 * itemp(:) )
      rate(:n,439) = 3.3e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,440) = 3e-12_r8 * exp( -1500._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)


      ko(:) = 4.4e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,411) = wrk(:)











      end subroutine setrxt_hrates

      end module mo_setrxt
