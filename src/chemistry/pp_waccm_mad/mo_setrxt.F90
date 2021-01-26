
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,100) = 9.6e-10_r8
      rate(:,:,101) = 1.3e-09_r8
      rate(:,:,102) = 2e-29_r8
      rate(:,:,103) = 1e-27_r8
      rate(:,:,104) = 1.6e-09_r8
      rate(:,:,105) = 6e-12_r8
      rate(:,:,106) = 2.9e-12_r8
      rate(:,:,107) = 2.9e-11_r8
      rate(:,:,108) = 2e-10_r8
      rate(:,:,109) = 1e-10_r8
      rate(:,:,110) = 1e-10_r8
      rate(:,:,111) = 1e-11_r8
      rate(:,:,112) = 1.7e-10_r8
      rate(:,:,113) = 1e-28_r8
      rate(:,:,114) = 1e-28_r8
      rate(:,:,115) = 4e-11_r8
      rate(:,:,116) = 4e-11_r8
      rate(:,:,117) = 3.5e-12_r8
      rate(:,:,118) = 3.5e-12_r8
      rate(:,:,119) = 3.51e-10_r8
      rate(:,:,120) = 1.1e-10_r8
      rate(:,:,121) = 6e-15_r8
      rate(:,:,122) = 1e-10_r8
      rate(:,:,123) = 1e-10_r8
      rate(:,:,124) = 2.2e-10_r8
      rate(:,:,125) = 1.2e-09_r8
      rate(:,:,126) = 1.4e-10_r8
      rate(:,:,127) = 1.3e-10_r8
      rate(:,:,133) = 1.5e-06_r8
      rate(:,:,134) = 2e-09_r8
      rate(:,:,135) = 1e-09_r8
      rate(:,:,136) = 3.6e-06_r8
      rate(:,:,137) = 4e-12_r8
      rate(:,:,138) = 1e-09_r8
      rate(:,:,139) = 5e-06_r8
      rate(:,:,140) = 7e-12_r8
      rate(:,:,280) = 1e-10_r8
      rate(:,:,281) = 1e-10_r8
      rate(:,:,282) = 3e-10_r8
      rate(:,:,283) = 1.6e-28_r8
      rate(:,:,284) = 1.4e-09_r8
      rate(:,:,285) = 1.6e-09_r8
      rate(:,:,286) = 2e-13_r8
      rate(:,:,287) = 1.2e-10_r8
      rate(:,:,288) = 7e-10_r8
      rate(:,:,289) = 1.6e-28_r8
      rate(:,:,290) = 1.6e-09_r8
      rate(:,:,291) = 1.6e-28_r8
      rate(:,:,292) = 7e-10_r8
      rate(:,:,293) = 1e-12_r8
      rate(:,:,294) = 7.6e-10_r8
      rate(:,:,295) = 1.45e-26_r8
      rate(:,:,296) = 5e-12_r8
      rate(:,:,297) = 1e-13_r8
      rate(:,:,298) = 2e-06_r8
      rate(:,:,299) = 2e-06_r8
      rate(:,:,300) = 7e-11_r8
      rate(:,:,301) = 1.5e-06_r8
      rate(:,:,302) = 1e-09_r8
      rate(:,:,303) = 1.5e-06_r8
      rate(:,:,304) = 7e-12_r8
      rate(:,:,305) = 5e-10_r8
      rate(:,:,306) = 1e-10_r8
      rate(:,:,307) = 1e-09_r8
      rate(:,:,308) = 1e-09_r8
      rate(:,:,309) = 1e-10_r8
      rate(:,:,310) = 1e-10_r8
      rate(:,:,311) = 9.9e-30_r8
      rate(:,:,312) = 1.4e-09_r8
      rate(:,:,313) = 1.6e-09_r8
      rate(:,:,314) = 2.9e-09_r8
      rate(:,:,315) = 7e-10_r8
      rate(:,:,316) = 2e-10_r8
      rate(:,:,317) = 3.4e-31_r8
      rate(:,:,318) = 7.8e-10_r8
      rate(:,:,319) = 1.5e-10_r8
      rate(:,:,320) = 1.5e-10_r8
      rate(:,:,321) = 2e-06_r8
      rate(:,:,322) = 9e-10_r8
      rate(:,:,323) = 2.4e-10_r8
      rate(:,:,324) = 2.8e-28_r8
      rate(:,:,325) = 5.5e-10_r8
      rate(:,:,326) = 8.4e-10_r8
      rate(:,:,327) = 1e-10_r8
      rate(:,:,328) = 1e-10_r8
      rate(:,:,329) = 2.5e-10_r8
      rate(:,:,330) = 4.3e-10_r8
      rate(:,:,331) = 4e-10_r8
      rate(:,:,332) = 1.7e-09_r8
      rate(:,:,333) = 3e-10_r8
      rate(:,:,334) = 1.5e-10_r8
      rate(:,:,336) = 1e-10_r8
      rate(:,:,337) = 1e-10_r8
      rate(:,:,338) = 7.6e-28_r8
      rate(:,:,339) = 1.4e-09_r8
      rate(:,:,340) = 1e-09_r8
      rate(:,:,341) = 1.1e-09_r8
      rate(:,:,342) = 2e-10_r8
      rate(:,:,343) = 9e-10_r8
      rate(:,:,345) = 1e-10_r8
      rate(:,:,346) = 1e-10_r8
      rate(:,:,347) = 2e-28_r8
      rate(:,:,348) = 5.8e-10_r8
      rate(:,:,349) = 3.2e-11_r8
      rate(:,:,350) = 6e-13_r8
      rate(:,:,351) = 2e-09_r8
      rate(:,:,352) = 3.6e-09_r8
      rate(:,:,353) = 5e-13_r8
      rate(:,:,354) = 1e-09_r8
      rate(:,:,355) = 1.9e-10_r8
      rate(:,:,356) = 3e-10_r8
      rate(:,:,357) = 2.9e-31_r8
      rate(:,:,358) = 8e-10_r8
      rate(:,:,382) = 0.000258_r8
      rate(:,:,383) = 0.085_r8
      rate(:,:,384) = 1.2e-10_r8
      rate(:,:,389) = 1.2e-10_r8
      rate(:,:,390) = 1e-20_r8
      rate(:,:,391) = 1.3e-16_r8
      rate(:,:,393) = 4.2e-13_r8
      rate(:,:,395) = 8e-14_r8
      rate(:,:,396) = 3.9e-17_r8
      rate(:,:,403) = 6.9e-12_r8
      rate(:,:,404) = 7.2e-11_r8
      rate(:,:,405) = 1.6e-12_r8
      rate(:,:,411) = 1.8e-12_r8
      rate(:,:,415) = 1.8e-12_r8
      rate(:,:,419) = 7e-13_r8
      rate(:,:,420) = 5e-12_r8
      rate(:,:,429) = 3.5e-12_r8
      rate(:,:,431) = 1e-11_r8
      rate(:,:,432) = 2.2e-11_r8
      rate(:,:,433) = 5e-11_r8
      rate(:,:,468) = 1.7e-13_r8
      rate(:,:,470) = 2.607e-10_r8
      rate(:,:,471) = 9.75e-11_r8
      rate(:,:,472) = 2.07e-10_r8
      rate(:,:,473) = 2.088e-10_r8
      rate(:,:,474) = 1.17e-10_r8
      rate(:,:,475) = 4.644e-11_r8
      rate(:,:,476) = 1.204e-10_r8
      rate(:,:,477) = 9.9e-11_r8
      rate(:,:,478) = 3.3e-12_r8
      rate(:,:,497) = 4.5e-11_r8
      rate(:,:,498) = 4.62e-10_r8
      rate(:,:,499) = 1.2e-10_r8
      rate(:,:,500) = 9e-11_r8
      rate(:,:,501) = 3e-11_r8
      rate(:,:,506) = 2.14e-11_r8
      rate(:,:,507) = 1.9e-10_r8
      rate(:,:,520) = 2.57e-10_r8
      rate(:,:,521) = 1.8e-10_r8
      rate(:,:,522) = 1.794e-10_r8
      rate(:,:,523) = 1.3e-10_r8
      rate(:,:,524) = 7.65e-11_r8
      rate(:,:,533) = 1.31e-10_r8
      rate(:,:,534) = 3.5e-11_r8
      rate(:,:,535) = 9e-12_r8
      rate(:,:,558) = 0.047_r8
      rate(:,:,559) = 7.7e-05_r8
      rate(:,:,560) = 0.171_r8
      rate(:,:,564) = 6e-11_r8
      rate(:,:,567) = 1e-12_r8
      rate(:,:,568) = 4e-10_r8
      rate(:,:,569) = 2e-10_r8
      rate(:,:,570) = 1e-10_r8
      rate(:,:,571) = 5e-16_r8
      rate(:,:,572) = 4.4e-10_r8
      rate(:,:,573) = 9e-10_r8
      rate(:,:,575) = 1.3e-10_r8
      rate(:,:,578) = 8e-10_r8
      rate(:,:,579) = 5e-12_r8
      rate(:,:,580) = 7e-10_r8
      rate(:,:,583) = 4.8e-10_r8
      rate(:,:,584) = 1e-10_r8
      rate(:,:,585) = 4e-10_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,335) = 1.8e-11_r8 * exp( 390._r8 * itemp(:,:) )
      rate(:,:,385) = 1.63e-10_r8 * exp( 60._r8 * itemp(:,:) )
      rate(:,:,386) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:,387) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:,388) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:,392) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:,394) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:,397) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:,398) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:,401) = 1.6e-11_r8 * exp( -4570._r8 * itemp(:,:) )
      rate(:,:,402) = 1.4e-12_r8 * exp( -2000._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 200._r8 * itemp(:,:) )
      rate(:,:,407) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,495) = 5.5e-12_r8 * exp_fac(:,:)
      rate(:,:,530) = 3.8e-12_r8 * exp_fac(:,:)
      rate(:,:,408) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:,409) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:,410) = 2.8e-12_r8 * exp( -1800._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 250._r8 * itemp(:,:) )
      rate(:,:,412) = 4.8e-11_r8 * exp_fac(:,:)
      rate(:,:,493) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,413) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:,414) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:,418) = 1.3e-12_r8 * exp( 380._r8 * itemp(:,:) )
      rate(:,:,421) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 220._r8 * itemp(:,:) )
      rate(:,:,422) = 2.9e-12_r8 * exp_fac(:,:)
      rate(:,:,423) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,424) = 1.45e-12_r8 * exp_fac(:,:)
      rate(:,:,425) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:,426) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._r8 * itemp(:,:) )
      rate(:,:,427) = 1.2e-13_r8 * exp_fac(:,:)
      rate(:,:,453) = 3e-11_r8 * exp_fac(:,:)
      rate(:,:,430) = 1.5e-11_r8 * exp( 170._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 270._r8 * itemp(:,:) )
      rate(:,:,434) = 3.3e-12_r8 * exp_fac(:,:)
      rate(:,:,449) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,463) = 7.4e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -1500._r8 * itemp(:,:) )
      rate(:,:,435) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,494) = 5.8e-12_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( 20._r8 * itemp(:,:) )
      rate(:,:,437) = 7.26e-11_r8 * exp_fac(:,:)
      rate(:,:,438) = 4.64e-11_r8 * exp_fac(:,:)
      rate(:,:,445) = 8.1e-11_r8 * exp( -30._r8 * itemp(:,:) )
      rate(:,:,446) = 7.1e-12_r8 * exp( -1270._r8 * itemp(:,:) )
      rate(:,:,447) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:,:) )
      rate(:,:,448) = 1.1e-11_r8 * exp( -980._r8 * itemp(:,:) )
      rate(:,:,450) = 3.6e-11_r8 * exp( -375._r8 * itemp(:,:) )
      rate(:,:,451) = 2.3e-11_r8 * exp( -200._r8 * itemp(:,:) )
      rate(:,:,452) = 3.3e-12_r8 * exp( -115._r8 * itemp(:,:) )
      rate(:,:,454) = 1e-12_r8 * exp( -1590._r8 * itemp(:,:) )
      rate(:,:,455) = 3.5e-13_r8 * exp( -1370._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 290._r8 * itemp(:,:) )
      rate(:,:,456) = 2.6e-12_r8 * exp_fac(:,:)
      rate(:,:,457) = 6.4e-12_r8 * exp_fac(:,:)
      rate(:,:,487) = 4.1e-13_r8 * exp_fac(:,:)
      rate(:,:,458) = 6.5e-12_r8 * exp( 135._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -840._r8 * itemp(:,:) )
      rate(:,:,460) = 3.6e-12_r8 * exp_fac(:,:)
      rate(:,:,509) = 2e-12_r8 * exp_fac(:,:)
      rate(:,:,461) = 1.2e-12_r8 * exp( -330._r8 * itemp(:,:) )
      rate(:,:,462) = 2.8e-11_r8 * exp( 85._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 230._r8 * itemp(:,:) )
      rate(:,:,464) = 6e-13_r8 * exp_fac(:,:)
      rate(:,:,484) = 1.5e-12_r8 * exp_fac(:,:)
      rate(:,:,492) = 1.9e-11_r8 * exp_fac(:,:)
      rate(:,:,465) = 1e-11_r8 * exp( -3300._r8 * itemp(:,:) )
      rate(:,:,466) = 1.8e-12_r8 * exp( -250._r8 * itemp(:,:) )
      rate(:,:,467) = 3.4e-12_r8 * exp( -130._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -500._r8 * itemp(:,:) )
      rate(:,:,469) = 3e-12_r8 * exp_fac(:,:)
      rate(:,:,503) = 1.4e-10_r8 * exp_fac(:,:)
      exp_fac(:,:) = exp( -800._r8 * itemp(:,:) )
      rate(:,:,481) = 1.7e-11_r8 * exp_fac(:,:)
      rate(:,:,508) = 6.3e-12_r8 * exp_fac(:,:)
      rate(:,:,482) = 4.8e-12_r8 * exp( -310._r8 * itemp(:,:) )
      rate(:,:,483) = 1.6e-11_r8 * exp( -780._r8 * itemp(:,:) )
      rate(:,:,485) = 9.5e-13_r8 * exp( 550._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 260._r8 * itemp(:,:) )
      rate(:,:,486) = 2.3e-12_r8 * exp_fac(:,:)
      rate(:,:,489) = 8.8e-12_r8 * exp_fac(:,:)
      rate(:,:,488) = 4.5e-12_r8 * exp( 460._r8 * itemp(:,:) )
      rate(:,:,491) = 1.9e-11_r8 * exp( 215._r8 * itemp(:,:) )
      rate(:,:,496) = 1.2e-10_r8 * exp( -430._r8 * itemp(:,:) )
      rate(:,:,502) = 1.6e-10_r8 * exp( -260._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 0._r8 * itemp(:,:) )
      rate(:,:,504) = 1.4e-11_r8 * exp_fac(:,:)
      rate(:,:,506) = 2.14e-11_r8 * exp_fac(:,:)
      rate(:,:,507) = 1.9e-10_r8 * exp_fac(:,:)
      rate(:,:,520) = 2.57e-10_r8 * exp_fac(:,:)
      rate(:,:,521) = 1.8e-10_r8 * exp_fac(:,:)
      rate(:,:,522) = 1.794e-10_r8 * exp_fac(:,:)
      rate(:,:,523) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,524) = 7.65e-11_r8 * exp_fac(:,:)
      rate(:,:,533) = 1.31e-10_r8 * exp_fac(:,:)
      rate(:,:,534) = 3.5e-11_r8 * exp_fac(:,:)
      rate(:,:,535) = 9e-12_r8 * exp_fac(:,:)
      rate(:,:,558) = 0.047_r8 * exp_fac(:,:)
      rate(:,:,559) = 7.7e-05_r8 * exp_fac(:,:)
      rate(:,:,560) = 0.171_r8 * exp_fac(:,:)
      rate(:,:,564) = 6e-11_r8 * exp_fac(:,:)
      rate(:,:,567) = 1e-12_r8 * exp_fac(:,:)
      rate(:,:,568) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,569) = 2e-10_r8 * exp_fac(:,:)
      rate(:,:,570) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,571) = 5e-16_r8 * exp_fac(:,:)
      rate(:,:,572) = 4.4e-10_r8 * exp_fac(:,:)
      rate(:,:,573) = 9e-10_r8 * exp_fac(:,:)
      rate(:,:,575) = 1.3e-10_r8 * exp_fac(:,:)
      rate(:,:,578) = 8e-10_r8 * exp_fac(:,:)
      rate(:,:,579) = 5e-12_r8 * exp_fac(:,:)
      rate(:,:,580) = 7e-10_r8 * exp_fac(:,:)
      rate(:,:,583) = 4.8e-10_r8 * exp_fac(:,:)
      rate(:,:,584) = 1e-10_r8 * exp_fac(:,:)
      rate(:,:,585) = 4e-10_r8 * exp_fac(:,:)
      rate(:,:,505) = 6e-12_r8 * exp( 400._r8 * itemp(:,:) )
      rate(:,:,510) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:,:) )
      rate(:,:,511) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:,:) )
      rate(:,:,512) = 1.64e-12_r8 * exp( -1520._r8 * itemp(:,:) )
      rate(:,:,513) = 2.03e-11_r8 * exp( -1100._r8 * itemp(:,:) )
      rate(:,:,514) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:,:) )
      rate(:,:,515) = 4.85e-12_r8 * exp( -850._r8 * itemp(:,:) )
      rate(:,:,516) = 9e-13_r8 * exp( -360._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._r8 * itemp(:,:) )
      rate(:,:,517) = 1.25e-12_r8 * exp_fac(:,:)
      rate(:,:,526) = 3.4e-11_r8 * exp_fac(:,:)
      rate(:,:,518) = 1.3e-12_r8 * exp( -1770._r8 * itemp(:,:) )
      rate(:,:,519) = 9.2e-13_r8 * exp( -1560._r8 * itemp(:,:) )
      rate(:,:,525) = 6e-13_r8 * exp( -2058._r8 * itemp(:,:) )
      rate(:,:,527) = 5.5e-12_r8 * exp( 125._r8 * itemp(:,:) )
      rate(:,:,528) = 4.1e-13_r8 * exp( 750._r8 * itemp(:,:) )
      rate(:,:,529) = 2.8e-12_r8 * exp( 300._r8 * itemp(:,:) )
      rate(:,:,531) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 7e-31_r8 * itemp(:,:)**2.6_r8
      kinf(:,:) = 3.6e-11_r8 * itemp(:,:)**0.1_r8
      call jpl( rate(1,1,344), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 4.4e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,406), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 6.9e-31_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 2.6e-11_r8
      call jpl( rate(1,1,416), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.5e-31_r8 * itemp(:,:)**1.8_r8
      kinf(:,:) = 2.2e-11_r8 * itemp(:,:)**0.7_r8
      call jpl( rate(1,1,428), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 9e-32_r8 * itemp(:,:)**1.5_r8
      kinf(:,:) = 3e-11_r8
      call jpl( rate(1,1,436), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 4e-12_r8 * itemp(:,:)**0.3_r8
      call jpl( rate(1,1,439), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 2.4e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 1.6e-12_r8 * itemp(:,:)**(-0.1_r8)
      call jpl( rate(1,1,440), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-30_r8 * itemp(:,:)**3._r8
      kinf(:,:) = 2.8e-11_r8
      call jpl( rate(1,1,441), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.8e-31_r8 * itemp(:,:)**3.4_r8
      kinf(:,:) = 1.5e-11_r8 * itemp(:,:)**1.9_r8
      call jpl( rate(1,1,459), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 1.9e-32_r8 * itemp(:,:)**3.6_r8
      kinf(:,:) = 3.7e-12_r8 * itemp(:,:)**1.6_r8
      call jpl( rate(1,1,479), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.2e-31_r8 * itemp(:,:)**3.2_r8
      kinf(:,:) = 6.9e-12_r8 * itemp(:,:)**2.9_r8
      call jpl( rate(1,1,490), m, 0.6_r8, ko, kinf, n )

      ko(:,:) = 5.9e-33_r8 * itemp(:,:)**1._r8
      kinf(:,:) = 1.1e-12_r8 * itemp(:,:)**(-1.3_r8)
      call jpl( rate(1,1,532), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)

      rate(:,:kbot,390) = 1e-20_r8
      rate(:,:kbot,391) = 1.3e-16_r8
      rate(:,:kbot,395) = 8e-14_r8
      rate(:,:kbot,396) = 3.9e-17_r8
      rate(:,:kbot,403) = 6.9e-12_r8
      rate(:,:kbot,419) = 7e-13_r8
      rate(:,:kbot,420) = 5e-12_r8
      rate(:,:kbot,558) = 0.047_r8
      rate(:,:kbot,559) = 7.7e-05_r8
      rate(:,:kbot,560) = 0.171_r8
      rate(:,:kbot,564) = 6e-11_r8
      rate(:,:kbot,567) = 1e-12_r8
      rate(:,:kbot,568) = 4e-10_r8
      rate(:,:kbot,569) = 2e-10_r8
      rate(:,:kbot,570) = 1e-10_r8
      rate(:,:kbot,572) = 4.4e-10_r8
      rate(:,:kbot,575) = 1.3e-10_r8
      rate(:,:kbot,578) = 8e-10_r8
      rate(:,:kbot,579) = 5e-12_r8
      rate(:,:kbot,580) = 7e-10_r8
      rate(:,:kbot,583) = 4.8e-10_r8
      rate(:,:kbot,584) = 1e-10_r8
      rate(:,:kbot,585) = 4e-10_r8
      itemp(:ncol,:kbot) = 1._r8 / temp(:ncol,:kbot)
      n = ncol*kbot
      rate(:,:kbot,386) = 2.15e-11_r8 * exp( 110._r8 * itemp(:,:) )
      exp_fac(:,:) = exp( 55._r8 * itemp(:,:) )
      rate(:,:kbot,387) = 2.64e-11_r8 * exp_fac(:,:)
      rate(:,:kbot,388) = 6.6e-12_r8 * exp_fac(:,:)
      rate(:,:kbot,392) = 3.6e-18_r8 * exp( -220._r8 * itemp(:,:) )
      rate(:,:kbot,394) = 1.8e-15_r8 * exp( 45._r8 * itemp(:,:) )
      rate(:,:kbot,397) = 3.5e-11_r8 * exp( -135._r8 * itemp(:,:) )
      rate(:,:kbot,398) = 8e-12_r8 * exp( -2060._r8 * itemp(:,:) )
      rate(:,:kbot,407) = 3e-11_r8 * exp( 200._r8 * itemp(:,:) )
      rate(:,:kbot,408) = 1e-14_r8 * exp( -490._r8 * itemp(:,:) )
      rate(:,:kbot,409) = 1.4e-10_r8 * exp( -470._r8 * itemp(:,:) )
      rate(:,:kbot,412) = 4.8e-11_r8 * exp( 250._r8 * itemp(:,:) )
      rate(:,:kbot,413) = 1.8e-11_r8 * exp( 180._r8 * itemp(:,:) )
      rate(:,:kbot,414) = 1.7e-12_r8 * exp( -940._r8 * itemp(:,:) )
      rate(:,:kbot,421) = 2.1e-11_r8 * exp( 100._r8 * itemp(:,:) )
      rate(:,:kbot,425) = 1.5e-11_r8 * exp( -3600._r8 * itemp(:,:) )
      rate(:,:kbot,426) = 5.1e-12_r8 * exp( 210._r8 * itemp(:,:) )
      rate(:,:kbot,434) = 3.3e-12_r8 * exp( 270._r8 * itemp(:,:) )
      rate(:,:kbot,435) = 3e-12_r8 * exp( -1500._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)


      ko(:,:) = 4.4e-32_r8 * itemp(:,:)**1.3_r8
      kinf(:,:) = 7.5e-11_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:,:kbot,406) = wrk(:,:)











      end subroutine setrxt_hrates

      end module mo_setrxt
