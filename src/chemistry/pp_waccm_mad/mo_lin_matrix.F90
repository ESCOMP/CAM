      module mo_lin_matrix
      use chem_mods, only: veclen
      private
      public :: linmat
      contains
      subroutine linmat01( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,16) = -( rxt(k,27) + het_rates(k,1) )
         mat(k,451) = -( rxt(k,28) + het_rates(k,2) )
         mat(k,79) = rxt(k,29)
         mat(k,76) = -( rxt(k,29) + rxt(k,30) + rxt(k,543) + rxt(k,546) + rxt(k,551) &
                 + het_rates(k,3) )
         mat(k,535) = -( rxt(k,21) + rxt(k,22) + het_rates(k,14) )
         mat(k,41) = rxt(k,23)
         mat(k,566) = rxt(k,534)*y(k,20) + rxt(k,535)*y(k,20)
         mat(k,252) = -( het_rates(k,18) )
         mat(k,1370) = rxt(k,446)*y(k,20)
         mat(k,175) = rxt(k,502)*y(k,20)
         mat(k,741) = rxt(k,531)*y(k,20)
         mat(k,561) = rxt(k,533)*y(k,20)
         mat(k,39) = -( rxt(k,23) + het_rates(k,19) )
         mat(k,10) = -( rxt(k,44) + het_rates(k,22) )
         mat(k,1) = -( rxt(k,45) + rxt(k,480) + het_rates(k,23) )
         mat(k,1520) = -( rxt(k,46) + het_rates(k,24) )
         mat(k,218) = rxt(k,48)
         mat(k,9) = rxt(k,60)
         mat(k,3) = 2.000_r8*rxt(k,480)
         mat(k,211) = -( rxt(k,47) + rxt(k,48) + rxt(k,545) + rxt(k,550) + rxt(k,556) &
                 + het_rates(k,25) )
         mat(k,142) = -( het_rates(k,27) )
         mat(k,530) = rxt(k,21) + rxt(k,22)
         mat(k,1067) = rxt(k,26) + rxt(k,62)
         mat(k,1366) = rxt(k,513)*y(k,17)
         mat(k,4) = -( rxt(k,49) + het_rates(k,29) )
         mat(k,555) = rxt(k,471)*y(k,6) + rxt(k,473)*y(k,9) + 2.000_r8*rxt(k,474)*y(k,10) &
                      + 2.000_r8*rxt(k,475)*y(k,11) + rxt(k,476)*y(k,12) &
                      + rxt(k,497)*y(k,7) + 2.000_r8*rxt(k,499)*y(k,34) &
                      + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39)
         mat(k,729) = rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39)
         mat(k,12) = -( rxt(k,50) + het_rates(k,30) )
         mat(k,556) = rxt(k,472)*y(k,8) + rxt(k,473)*y(k,9) + rxt(k,522)*y(k,37)
         mat(k,730) = rxt(k,517)*y(k,37)
         mat(k,1086) = -( rxt(k,26) + rxt(k,62) + het_rates(k,28) )
         mat(k,786) = rxt(k,63)
         mat(k,625) = rxt(k,65)
         mat(k,135) = rxt(k,364)
         mat(k,174) = -( rxt(k,502)*y(k,20) + het_rates(k,31) )
         mat(k,5) = 2.000_r8*rxt(k,49)
         mat(k,13) = rxt(k,50)
         mat(k,20) = rxt(k,57)
         mat(k,558) = rxt(k,475)*y(k,11) + rxt(k,497)*y(k,7)
         mat(k,1716) = -( het_rates(k,32) )
         mat(k,1814) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,552) = 2.000_r8*rxt(k,21)
         mat(k,43) = rxt(k,23)
         mat(k,188) = rxt(k,52)
         mat(k,1146) = rxt(k,56)
         mat(k,21) = rxt(k,57)
         mat(k,586) = rxt(k,534)*y(k,20)
         mat(k,413) = -( het_rates(k,33) )
         mat(k,1781) = rxt(k,1)
         mat(k,533) = rxt(k,22)
         mat(k,564) = rxt(k,535)*y(k,20)
         mat(k,84) = -( rxt(k,4) + het_rates(k,35) )
         mat(k,669) = .500_r8*rxt(k,537)
         mat(k,182) = -( rxt(k,52) + het_rates(k,36) )
         mat(k,1132) = -( rxt(k,56) + het_rates(k,40) )
         mat(k,341) = rxt(k,381)
         mat(k,1390) = rxt(k,446)*y(k,20) + rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + 2.000_r8*rxt(k,513)*y(k,17) + rxt(k,515)*y(k,21)
         mat(k,19) = -( rxt(k,57) + het_rates(k,41) )
         mat(k,173) = rxt(k,502)*y(k,20)
         mat(k,1173) = -( rxt(k,9) + het_rates(k,42) )
         mat(k,439) = rxt(k,270)
         mat(k,244) = 2.000_r8*rxt(k,538) + 2.000_r8*rxt(k,541) + 2.000_r8*rxt(k,544) &
                      + 2.000_r8*rxt(k,555)
         mat(k,1748) = .500_r8*rxt(k,539)
         mat(k,1426) = rxt(k,540)
         mat(k,80) = rxt(k,543) + rxt(k,546) + rxt(k,551)
         mat(k,214) = rxt(k,545) + rxt(k,550) + rxt(k,556)
         mat(k,680) = -( rxt(k,537) + het_rates(k,43) )
         mat(k,52) = rxt(k,11) + rxt(k,443)
         mat(k,1379) = rxt(k,510)*y(k,15) + rxt(k,513)*y(k,17)
         mat(k,748) = rxt(k,511)*y(k,15) + rxt(k,514)*y(k,17)
         mat(k,568) = rxt(k,534)*y(k,20)
         mat(k,51) = -( rxt(k,10) + rxt(k,11) + rxt(k,443) + het_rates(k,44) )
         mat(k,165) = -( rxt(k,58) + het_rates(k,45) )
         mat(k,77) = rxt(k,543) + rxt(k,546) + rxt(k,551)
         mat(k,191) = -( rxt(k,59) + het_rates(k,46) )
         mat(k,210) = rxt(k,545) + rxt(k,550) + rxt(k,556)
         mat(k,137) = -( rxt(k,12) + het_rates(k,47) )
         mat(k,265) = -( rxt(k,66) + het_rates(k,48) )
         mat(k,982) = rxt(k,17)
         mat(k,200) = rxt(k,574)
         mat(k,239) = -( rxt(k,14) + rxt(k,15) + rxt(k,444) + rxt(k,538) + rxt(k,541) &
                      + rxt(k,544) + rxt(k,555) + het_rates(k,50) )
         mat(k,1001) = -( rxt(k,16) + rxt(k,17) + het_rates(k,51) )
         mat(k,139) = rxt(k,12)
         mat(k,243) = rxt(k,15)
         mat(k,1744) = rxt(k,18) + .500_r8*rxt(k,539)
         mat(k,1422) = rxt(k,20)
         mat(k,1344) = rxt(k,571)
         mat(k,63) = rxt(k,583)
         mat(k,572) = 2.000_r8*rxt(k,437)*y(k,49)
         mat(k,1762) = -( rxt(k,18) + rxt(k,539) + het_rates(k,52) )
         mat(k,1187) = rxt(k,9)
         mat(k,56) = rxt(k,11) + rxt(k,443)
         mat(k,249) = rxt(k,14) + rxt(k,444)
         mat(k,1440) = rxt(k,19)
         mat(k,83) = rxt(k,29)
         mat(k,219) = rxt(k,48)
         mat(k,727) = rxt(k,75)
         mat(k,1432) = -( rxt(k,19) + rxt(k,20) + rxt(k,540) + het_rates(k,53) )
         mat(k,54) = rxt(k,10)
         mat(k,246) = rxt(k,14) + rxt(k,15) + rxt(k,444)
         mat(k,82) = rxt(k,30)
         mat(k,217) = rxt(k,47)
         mat(k,660) = rxt(k,76)
         mat(k,1230) = -( rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82) + het_rates(k,54) )
         mat(k,1802) = rxt(k,2)
         mat(k,1668) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,83) + 2.000_r8*rxt(k,84) &
                      + rxt(k,85) + rxt(k,87) + 2.000_r8*rxt(k,89) + rxt(k,90) + rxt(k,91) &
                      + rxt(k,92)
         mat(k,1268) = rxt(k,8)
         mat(k,245) = rxt(k,15)
         mat(k,1006) = rxt(k,17)
         mat(k,1749) = rxt(k,18)
         mat(k,1427) = rxt(k,19)
         mat(k,1089) = rxt(k,26) + rxt(k,62)
         mat(k,460) = rxt(k,28)
         mat(k,1513) = rxt(k,46)
         mat(k,8) = rxt(k,60)
         mat(k,1625) = rxt(k,99) + rxt(k,353)
         mat(k,576) = rxt(k,386)
         mat(k,68) = rxt(k,577)
         mat(k,64) = rxt(k,582)
         mat(k,1679) = -( rxt(k,5) + rxt(k,6) + rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86) + rxt(k,87) + rxt(k,88) + rxt(k,89) + rxt(k,90) &
                      + rxt(k,91) + rxt(k,92) + het_rates(k,55) )
         mat(k,1279) = rxt(k,8)
         mat(k,1438) = rxt(k,20)
         mat(k,1315) = rxt(k,93) + rxt(k,132)
         mat(k,407) = rxt(k,95)
         mat(k,104) = rxt(k,97)
         mat(k,236) = rxt(k,382) + rxt(k,390)
         mat(k,32) = rxt(k,383)
         mat(k,585) = rxt(k,438)*y(k,49)
         mat(k,1269) = -( rxt(k,7) + rxt(k,8) + het_rates(k,56) )
         mat(k,403) = rxt(k,96)
         mat(k,7) = -( rxt(k,60) + het_rates(k,57) )
         mat(k,322) = -( het_rates(k,59) )
         mat(k,17) = rxt(k,27)
         mat(k,450) = rxt(k,28)
         mat(k,78) = rxt(k,30)
         mat(k,183) = rxt(k,52)
         mat(k,166) = rxt(k,58)
         mat(k,563) = rxt(k,471)*y(k,6) + rxt(k,497)*y(k,7) + 3.000_r8*rxt(k,498)*y(k,21) &
                      + 2.000_r8*rxt(k,499)*y(k,34) + 2.000_r8*rxt(k,520)*y(k,13) &
                      + rxt(k,521)*y(k,15)
         mat(k,1371) = 2.000_r8*rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + 3.000_r8*rxt(k,515)*y(k,21)
         mat(k,743) = 2.000_r8*rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) &
                      + 3.000_r8*rxt(k,516)*y(k,21)
         mat(k,1396) = -( rxt(k,446)*y(k,20) + rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + rxt(k,513)*y(k,17) + rxt(k,515)*y(k,21) + het_rates(k,60) )
         mat(k,18) = rxt(k,27)
         mat(k,11) = 2.000_r8*rxt(k,44)
         mat(k,2) = 2.000_r8*rxt(k,45)
         mat(k,1517) = rxt(k,46)
         mat(k,216) = rxt(k,47)
         mat(k,15) = rxt(k,50)
         mat(k,1138) = rxt(k,56)
         mat(k,195) = rxt(k,59)
         mat(k,580) = 4.000_r8*rxt(k,470)*y(k,5) + rxt(k,471)*y(k,6) &
                      + 2.000_r8*rxt(k,472)*y(k,8) + 2.000_r8*rxt(k,473)*y(k,9) &
                      + 2.000_r8*rxt(k,474)*y(k,10) + rxt(k,475)*y(k,11) &
                      + 2.000_r8*rxt(k,476)*y(k,12) + rxt(k,522)*y(k,37) &
                      + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39)
         mat(k,762) = 3.000_r8*rxt(k,512)*y(k,16) + rxt(k,514)*y(k,17) &
                      + rxt(k,517)*y(k,37) + rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39)
         mat(k,591) = -( het_rates(k,61) )
         mat(k,351) = rxt(k,380)
         mat(k,335) = rxt(k,381)
         mat(k,350) = -( rxt(k,380) + het_rates(k,62) )
         mat(k,334) = -( rxt(k,381) + het_rates(k,63) )
         mat(k,221) = -( het_rates(k,64) )
         mat(k,779) = -( rxt(k,63) + het_rates(k,65) )
         mat(k,495) = rxt(k,64) + rxt(k,278)
         mat(k,305) = -( rxt(k,279) + het_rates(k,66) )
         mat(k,491) = -( rxt(k,64) + rxt(k,278) + het_rates(k,67) )
         mat(k,306) = rxt(k,279)
         mat(k,613) = -( rxt(k,65) + het_rates(k,68) )
         mat(k,958) = -( het_rates(k,69) )
         mat(k,1000) = rxt(k,16)
         mat(k,269) = rxt(k,66)
         mat(k,711) = rxt(k,75)
         mat(k,651) = rxt(k,76)
         mat(k,1224) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82)
         mat(k,1662) = rxt(k,83) + rxt(k,85) + rxt(k,86) + rxt(k,87) + rxt(k,88) &
                      + rxt(k,90) + rxt(k,91) + rxt(k,92)
         mat(k,1298) = rxt(k,93) + rxt(k,132)
         mat(k,399) = rxt(k,96)
         mat(k,1583) = rxt(k,98)
         mat(k,1619) = rxt(k,99) + rxt(k,353)
         mat(k,33) = -( het_rates(k,70) )
         mat(k,291) = -( het_rates(k,71) )
      end do
      end subroutine linmat01
      subroutine linmat02( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
      do k = 1,avec_len
         mat(k,105) = -( rxt(k,368) + het_rates(k,72) )
         mat(k,804) = rxt(k,370)
         mat(k,823) = -( rxt(k,370) + het_rates(k,73) )
         mat(k,1461) = rxt(k,372)
         mat(k,22) = -( het_rates(k,74) )
         mat(k,1476) = -( rxt(k,372) + het_rates(k,75) )
         mat(k,1562) = rxt(k,374)
         mat(k,26) = -( het_rates(k,76) )
         mat(k,1564) = -( rxt(k,374) + het_rates(k,77) )
         mat(k,148) = -( het_rates(k,78) )
         mat(k,106) = rxt(k,368)
         mat(k,71) = -( het_rates(k,79) )
         mat(k,112) = -( het_rates(k,80) )
         mat(k,66) = rxt(k,577)
         mat(k,61) = rxt(k,582)
         mat(k,706) = -( rxt(k,75) + het_rates(k,81) )
         mat(k,366) = rxt(k,268)
         mat(k,365) = -( rxt(k,268) + het_rates(k,82) )
         mat(k,646) = -( rxt(k,76) + het_rates(k,83) )
         mat(k,433) = rxt(k,270)
         mat(k,514) = rxt(k,277)
         mat(k,471) = -( rxt(k,269) + het_rates(k,84) )
         mat(k,513) = -( rxt(k,277) + het_rates(k,85) )
         mat(k,472) = rxt(k,269)
         mat(k,379) = -( het_rates(k,86) )
         mat(k,432) = -( rxt(k,270) + het_rates(k,87) )
         mat(k,870) = -( rxt(k,365) + het_rates(k,88) )
         mat(k,998) = rxt(k,16)
         mat(k,131) = rxt(k,364)
         mat(k,124) = rxt(k,366)
         mat(k,1341) = rxt(k,571)
         mat(k,203) = rxt(k,574)
         mat(k,914) = -( het_rates(k,89) )
         mat(k,45) = -( het_rates(k,90) )
         mat(k,130) = -( rxt(k,364) + het_rates(k,91) )
         mat(k,1044) = -( het_rates(k,92) )
         mat(k,122) = -( rxt(k,366) + het_rates(k,93) )
         mat(k,848) = rxt(k,365)
         mat(k,91) = -( het_rates(k,94) )
         mat(k,263) = rxt(k,66)
         mat(k,60) = rxt(k,583)
         mat(k,567) = -( rxt(k,386) + rxt(k,437)*y(k,49) + rxt(k,438)*y(k,49) &
                      + rxt(k,470)*y(k,5) + rxt(k,471)*y(k,6) + rxt(k,472)*y(k,8) &
                      + rxt(k,473)*y(k,9) + rxt(k,474)*y(k,10) + rxt(k,475)*y(k,11) &
                      + rxt(k,476)*y(k,12) + rxt(k,497)*y(k,7) + rxt(k,498)*y(k,21) &
                      + rxt(k,499)*y(k,34) + rxt(k,520)*y(k,13) + rxt(k,521)*y(k,15) &
                      + rxt(k,522)*y(k,37) + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39) &
                      + rxt(k,533)*y(k,20) + rxt(k,534)*y(k,20) + rxt(k,535)*y(k,20) &
                 + het_rates(k,95) )
         mat(k,1786) = rxt(k,1)
         mat(k,1652) = rxt(k,6)
         mat(k,1252) = rxt(k,7)
         mat(k,229) = -( rxt(k,382) + rxt(k,390) + het_rates(k,96) )
         mat(k,1246) = rxt(k,7)
         mat(k,31) = rxt(k,394)
         mat(k,30) = -( rxt(k,383) + rxt(k,394) + het_rates(k,97) )
         mat(k,1306) = -( rxt(k,93) + rxt(k,132) + het_rates(k,98) )
         mat(k,630) = rxt(k,65)
         mat(k,103) = rxt(k,97)
         mat(k,1352) = -( rxt(k,571) + het_rates(k,99) )
         mat(k,1671) = rxt(k,86) + rxt(k,88)
         mat(k,160) = rxt(k,94)
         mat(k,156) = -( rxt(k,94) + het_rates(k,100) )
         mat(k,395) = -( rxt(k,95) + rxt(k,96) + het_rates(k,101) )
         mat(k,98) = -( rxt(k,97) + het_rates(k,102) )
         mat(k,279) = -( het_rates(k,103) )
         mat(k,749) = -( rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) + rxt(k,512)*y(k,16) &
                      + rxt(k,514)*y(k,17) + rxt(k,516)*y(k,21) + rxt(k,517)*y(k,37) &
                      + rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39) + rxt(k,531)*y(k,20) &
                 + het_rates(k,104) )
         mat(k,1791) = rxt(k,3)
         mat(k,86) = 2.000_r8*rxt(k,4)
         mat(k,1163) = rxt(k,9)
         mat(k,53) = rxt(k,10)
         mat(k,138) = rxt(k,12)
         mat(k,42) = rxt(k,23)
         mat(k,168) = rxt(k,58)
         mat(k,192) = rxt(k,59)
         mat(k,1578) = rxt(k,98)
         mat(k,1738) = .500_r8*rxt(k,539)
         mat(k,569) = rxt(k,533)*y(k,20)
         mat(k,1598) = -( rxt(k,98) + het_rates(k,105) )
         mat(k,1635) = -( rxt(k,99) + rxt(k,353) + het_rates(k,106) )
         mat(k,799) = rxt(k,63)
         mat(k,406) = rxt(k,95)
         mat(k,199) = -( rxt(k,574) + het_rates(k,107) )
         mat(k,1202) = rxt(k,78) + rxt(k,79)
         mat(k,1647) = rxt(k,85) + rxt(k,87)
         mat(k,62) = rxt(k,558)
         mat(k,67) = rxt(k,559)
         mat(k,65) = -( rxt(k,559) + rxt(k,577) + het_rates(k,108) )
         mat(k,1191) = rxt(k,80) + rxt(k,81)
         mat(k,1642) = rxt(k,90) + rxt(k,91)
         mat(k,59) = rxt(k,560)
         mat(k,58) = -( rxt(k,558) + rxt(k,560) + rxt(k,582) + rxt(k,583) &
                 + het_rates(k,109) )
         mat(k,1190) = rxt(k,77) + rxt(k,82)
         mat(k,1641) = rxt(k,83) + rxt(k,92)
         mat(k,1816) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,110) )
         mat(k,509) = rxt(k,64) + rxt(k,278)
         mat(k,163) = rxt(k,94)
         mat(k,378) = rxt(k,268)
         mat(k,488) = rxt(k,269)
         mat(k,529) = rxt(k,277)
         mat(k,319) = rxt(k,279)
         mat(k,110) = rxt(k,368)
         mat(k,846) = rxt(k,370)
         mat(k,1484) = rxt(k,372)
         mat(k,1570) = rxt(k,374)
         mat(k,364) = rxt(k,380)
         mat(k,770) = rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) + rxt(k,512)*y(k,16) &
                      + rxt(k,514)*y(k,17) + rxt(k,519)*y(k,39) + rxt(k,531)*y(k,20)
      end do
      end subroutine linmat02
      subroutine linmat( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(in) :: het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
      call linmat01( avec_len, mat, y, rxt, het_rates )
      call linmat02( avec_len, mat, y, rxt, het_rates )
      end subroutine linmat
      end module mo_lin_matrix
