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
         mat(k,1061) = -( rxt(k,3) + rxt(k,4) + het_rates(k,1) )
         mat(k,464) = rxt(k,96)
         mat(k,1120) = -( rxt(k,67) + rxt(k,68) + rxt(k,69) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82) + het_rates(k,2) )
         mat(k,1206) = rxt(k,1) + 2.000_r8*rxt(k,2) + rxt(k,73) + rxt(k,74) + rxt(k,75) &
                      + 2.000_r8*rxt(k,78) + rxt(k,85) + rxt(k,86) + rxt(k,87) &
                      + 2.000_r8*rxt(k,90)
         mat(k,1062) = rxt(k,4)
         mat(k,1513) = rxt(k,6)
         mat(k,1559) = rxt(k,8)
         mat(k,288) = rxt(k,10)
         mat(k,1633) = rxt(k,12)
         mat(k,1343) = rxt(k,21)
         mat(k,1679) = rxt(k,24)
         mat(k,67) = rxt(k,25)
         mat(k,487) = rxt(k,32)
         mat(k,1427) = rxt(k,58) + rxt(k,92)
         mat(k,251) = rxt(k,62)
         mat(k,53) = rxt(k,63)
         mat(k,360) = rxt(k,65)
         mat(k,1284) = rxt(k,94) + rxt(k,367)
         mat(k,637) = rxt(k,118)
         mat(k,631) = -( rxt(k,118) + rxt(k,122)*y(k,4) + rxt(k,123)*y(k,4) &
                      + rxt(k,125)*y(k,37) + rxt(k,126)*y(k,38) + rxt(k,127)*y(k,39) &
                      + rxt(k,128)*y(k,47) + rxt(k,129)*y(k,48) + rxt(k,130)*y(k,40) &
                      + rxt(k,131)*y(k,45) + rxt(k,132)*y(k,46) + rxt(k,133)*y(k,41) &
                      + rxt(k,134)*y(k,36) + rxt(k,135)*y(k,44) + rxt(k,136)*y(k,43) &
                      + rxt(k,137)*y(k,49) + rxt(k,138)*y(k,50) + rxt(k,139)*y(k,51) &
                      + rxt(k,140)*y(k,52) + rxt(k,143)*y(k,12) + rxt(k,144)*y(k,12) &
                      + rxt(k,145)*y(k,12) + het_rates(k,134) )
         mat(k,1194) = rxt(k,1)
         mat(k,1050) = rxt(k,3)
         mat(k,1331) = rxt(k,20)
         mat(k,1208) = -( rxt(k,1) + rxt(k,2) + rxt(k,71) + rxt(k,73) + rxt(k,74) &
                      + rxt(k,75) + rxt(k,78) + rxt(k,83) + rxt(k,85) + rxt(k,86) &
                      + rxt(k,87) + rxt(k,90) + het_rates(k,3) )
         mat(k,1064) = rxt(k,4)
         mat(k,1635) = rxt(k,13)
         mat(k,1804) = rxt(k,95) + rxt(k,371)
         mat(k,467) = rxt(k,100)
         mat(k,126) = rxt(k,101)
         mat(k,63) = rxt(k,113)
         mat(k,277) = rxt(k,116) + rxt(k,117)
         mat(k,639) = rxt(k,123)*y(k,4)
         mat(k,61) = -( rxt(k,110) + rxt(k,113) + het_rates(k,132) )
         mat(k,273) = -( rxt(k,116) + rxt(k,117) + het_rates(k,133) )
         mat(k,1044) = rxt(k,3)
         mat(k,62) = rxt(k,110)
         mat(k,575) = -( het_rates(k,17) )
         mat(k,446) = rxt(k,18)
         mat(k,1329) = rxt(k,20)
         mat(k,629) = rxt(k,145)*y(k,12)
         mat(k,186) = -( het_rates(k,16) )
         mat(k,443) = rxt(k,17) + rxt(k,18)
         mat(k,1409) = rxt(k,58) + rxt(k,92)
         mat(k,71) = rxt(k,64)
         mat(k,1823) = rxt(k,240)*y(k,35)
         mat(k,1434) = -( rxt(k,58) + rxt(k,92) + het_rates(k,57) )
         mat(k,862) = rxt(k,102)
         mat(k,695) = rxt(k,103)
         mat(k,165) = rxt(k,361)
         mat(k,309) = -( rxt(k,70) + het_rates(k,5) )
         mat(k,1492) = rxt(k,6)
         mat(k,237) = rxt(k,300)
         mat(k,1522) = -( rxt(k,6) + rxt(k,7) + het_rates(k,6) )
         mat(k,1568) = rxt(k,8) + .500_r8*rxt(k,263)
         mat(k,291) = rxt(k,10)
         mat(k,1642) = rxt(k,13)
         mat(k,218) = rxt(k,66)
         mat(k,1775) = rxt(k,310)
         mat(k,644) = 2.000_r8*rxt(k,122)*y(k,4)
         mat(k,1569) = -( rxt(k,8) + rxt(k,263) + het_rates(k,7) )
         mat(k,292) = rxt(k,9) + rxt(k,183)
         mat(k,1730) = rxt(k,11)
         mat(k,1643) = rxt(k,12)
         mat(k,107) = rxt(k,15) + rxt(k,192)
         mat(k,260) = rxt(k,30)
         mat(k,119) = rxt(k,36)
         mat(k,785) = rxt(k,98)
         mat(k,819) = -( rxt(k,241)*y(k,35) + rxt(k,242)*y(k,42) + rxt(k,243)*y(k,40) &
                      + rxt(k,244)*y(k,36) + rxt(k,246)*y(k,45) + rxt(k,247)*y(k,46) &
                      + rxt(k,248)*y(k,52) + rxt(k,249)*y(k,51) + rxt(k,252)*y(k,12) &
                 + het_rates(k,87) )
         mat(k,1713) = rxt(k,11)
         mat(k,104) = rxt(k,14)
         mat(k,92) = rxt(k,16)
         mat(k,1336) = rxt(k,19)
         mat(k,131) = 2.000_r8*rxt(k,22)
         mat(k,229) = rxt(k,27)
         mat(k,195) = rxt(k,33)
         mat(k,216) = rxt(k,66)
         mat(k,1585) = rxt(k,97)
         mat(k,1552) = .500_r8*rxt(k,263)
         mat(k,633) = rxt(k,143)*y(k,12)
         mat(k,1645) = -( rxt(k,12) + rxt(k,13) + rxt(k,262) + het_rates(k,8) )
         mat(k,293) = rxt(k,9) + rxt(k,10) + rxt(k,183)
         mat(k,108) = rxt(k,14)
         mat(k,261) = rxt(k,29)
         mat(k,120) = rxt(k,35)
         mat(k,726) = rxt(k,99)
         mat(k,215) = -( rxt(k,66) + het_rates(k,20) )
         mat(k,1734) = -( rxt(k,11) + het_rates(k,9) )
         mat(k,294) = 2.000_r8*rxt(k,261) + 2.000_r8*rxt(k,282) + 2.000_r8*rxt(k,288) &
                      + 2.000_r8*rxt(k,293)
         mat(k,1647) = rxt(k,262)
         mat(k,1573) = .500_r8*rxt(k,263)
         mat(k,263) = rxt(k,283) + rxt(k,289) + rxt(k,294)
         mat(k,121) = rxt(k,284) + rxt(k,292) + rxt(k,295)
         mat(k,507) = rxt(k,462)
         mat(k,102) = -( rxt(k,14) + rxt(k,15) + rxt(k,192) + het_rates(k,10) )
         mat(k,283) = -( rxt(k,9) + rxt(k,10) + rxt(k,183) + rxt(k,261) + rxt(k,282) &
                      + rxt(k,288) + rxt(k,293) + het_rates(k,11) )
         mat(k,296) = -( het_rates(k,13) )
         mat(k,626) = rxt(k,143)*y(k,12)
         mat(k,1827) = rxt(k,199)*y(k,12)
         mat(k,202) = rxt(k,238)*y(k,12)
         mat(k,810) = rxt(k,252)*y(k,12)
         mat(k,89) = -( rxt(k,16) + het_rates(k,14) )
         mat(k,445) = -( rxt(k,17) + rxt(k,18) + het_rates(k,15) )
         mat(k,91) = rxt(k,16)
         mat(k,628) = rxt(k,144)*y(k,12) + rxt(k,145)*y(k,12)
         mat(k,931) = -( het_rates(k,18) )
         mat(k,93) = rxt(k,16)
         mat(k,451) = 2.000_r8*rxt(k,17)
         mat(k,1339) = rxt(k,19) + 2.000_r8*rxt(k,21)
         mat(k,1018) = rxt(k,28)
         mat(k,224) = rxt(k,34)
         mat(k,41) = rxt(k,57)
         mat(k,634) = rxt(k,144)*y(k,12)
         mat(k,744) = -( rxt(k,264) + het_rates(k,88) )
         mat(k,103) = rxt(k,15) + rxt(k,192)
         mat(k,632) = rxt(k,144)*y(k,12)
         mat(k,1835) = rxt(k,240)*y(k,35) + rxt(k,245)*y(k,36)
         mat(k,818) = rxt(k,241)*y(k,35) + rxt(k,244)*y(k,36)
         mat(k,129) = -( rxt(k,22) + het_rates(k,19) )
         mat(k,733) = .500_r8*rxt(k,264)
         mat(k,1348) = -( rxt(k,19) + rxt(k,20) + rxt(k,21) + het_rates(k,135) )
         mat(k,29) = rxt(k,61)
         mat(k,183) = rxt(k,93)
         mat(k,543) = rxt(k,104) + rxt(k,452)
         mat(k,140) = rxt(k,341)
         mat(k,1391) = rxt(k,343)
         mat(k,1475) = rxt(k,345)
         mat(k,983) = rxt(k,347)
         mat(k,406) = rxt(k,445)
         mat(k,562) = rxt(k,449)
         mat(k,375) = rxt(k,455)
         mat(k,420) = rxt(k,457)
         mat(k,521) = rxt(k,460)
         mat(k,828) = rxt(k,241)*y(k,35) + rxt(k,242)*y(k,42) + rxt(k,243)*y(k,40) &
                      + rxt(k,244)*y(k,36) + rxt(k,248)*y(k,52) + rxt(k,252)*y(k,12)
         mat(k,1861) = -( rxt(k,199)*y(k,12) + rxt(k,240)*y(k,35) + rxt(k,245)*y(k,36) &
                      + rxt(k,250)*y(k,52) + rxt(k,251)*y(k,51) + het_rates(k,85) )
         mat(k,31) = 2.000_r8*rxt(k,23)
         mat(k,1696) = rxt(k,24)
         mat(k,23) = 2.000_r8*rxt(k,26)
         mat(k,234) = rxt(k,27)
         mat(k,1039) = rxt(k,28)
         mat(k,264) = rxt(k,29)
         mat(k,38) = rxt(k,31)
         mat(k,35) = rxt(k,56)
         mat(k,652) = 2.000_r8*rxt(k,125)*y(k,37) + 2.000_r8*rxt(k,126)*y(k,38) &
                      + 2.000_r8*rxt(k,127)*y(k,39) + 2.000_r8*rxt(k,128)*y(k,47) &
                      + rxt(k,129)*y(k,48) + rxt(k,130)*y(k,40) + rxt(k,131)*y(k,45) &
                      + rxt(k,132)*y(k,46) + 4.000_r8*rxt(k,133)*y(k,41) &
                      + rxt(k,135)*y(k,44)
         mat(k,839) = rxt(k,241)*y(k,35) + 3.000_r8*rxt(k,242)*y(k,42) &
                      + rxt(k,243)*y(k,40) + rxt(k,246)*y(k,45) + rxt(k,247)*y(k,46)
         mat(k,30) = -( rxt(k,23) + het_rates(k,23) )
         mat(k,1692) = -( rxt(k,24) + het_rates(k,24) )
         mat(k,68) = rxt(k,25)
         mat(k,262) = rxt(k,30)
         mat(k,22) = 2.000_r8*rxt(k,211)
         mat(k,64) = -( rxt(k,25) + het_rates(k,25) )
         mat(k,21) = -( rxt(k,26) + rxt(k,211) + het_rates(k,26) )
         mat(k,1020) = -( rxt(k,28) + het_rates(k,27) )
         mat(k,387) = rxt(k,446)
         mat(k,1842) = rxt(k,199)*y(k,12) + 2.000_r8*rxt(k,240)*y(k,35) &
                      + rxt(k,245)*y(k,36) + rxt(k,250)*y(k,52) + rxt(k,251)*y(k,51)
         mat(k,228) = -( rxt(k,27) + het_rates(k,28) )
         mat(k,254) = rxt(k,283) + rxt(k,289) + rxt(k,294)
         mat(k,255) = -( rxt(k,29) + rxt(k,30) + rxt(k,283) + rxt(k,289) + rxt(k,294) &
                 + het_rates(k,29) )
         mat(k,36) = -( rxt(k,31) + het_rates(k,30) )
         mat(k,600) = -( het_rates(k,86) )
         mat(k,37) = rxt(k,31)
         mat(k,481) = rxt(k,32)
         mat(k,194) = rxt(k,33)
         mat(k,221) = rxt(k,34)
         mat(k,117) = rxt(k,35)
         mat(k,630) = rxt(k,134)*y(k,36) + rxt(k,135)*y(k,44) + rxt(k,136)*y(k,43) &
                      + 2.000_r8*rxt(k,137)*y(k,49) + 2.000_r8*rxt(k,138)*y(k,50) &
                      + 3.000_r8*rxt(k,139)*y(k,51) + 2.000_r8*rxt(k,140)*y(k,52)
         mat(k,816) = rxt(k,244)*y(k,36) + 2.000_r8*rxt(k,248)*y(k,52) &
                      + 3.000_r8*rxt(k,249)*y(k,51)
         mat(k,1830) = rxt(k,245)*y(k,36) + 2.000_r8*rxt(k,250)*y(k,52) &
                      + 3.000_r8*rxt(k,251)*y(k,51)
         mat(k,480) = -( rxt(k,32) + het_rates(k,31) )
         mat(k,116) = rxt(k,36)
         mat(k,220) = -( rxt(k,34) + het_rates(k,32) )
         mat(k,192) = -( rxt(k,33) + het_rates(k,33) )
         mat(k,115) = rxt(k,284) + rxt(k,292) + rxt(k,295)
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
         mat(k,114) = -( rxt(k,35) + rxt(k,36) + rxt(k,284) + rxt(k,292) + rxt(k,295) &
                 + het_rates(k,34) )
         mat(k,143) = -( het_rates(k,89) )
         mat(k,1781) = -( rxt(k,310) + het_rates(k,90) )
         mat(k,1221) = rxt(k,71) + rxt(k,83)
         mat(k,185) = rxt(k,93)
         mat(k,95) = -( het_rates(k,127) )
         mat(k,307) = rxt(k,70)
         mat(k,236) = -( rxt(k,300) + het_rates(k,128) )
         mat(k,1094) = rxt(k,67) + rxt(k,68) + rxt(k,69) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82)
         mat(k,1185) = rxt(k,73) + rxt(k,74) + rxt(k,75) + rxt(k,85) + rxt(k,86) &
                      + rxt(k,87)
         mat(k,1908) = -( rxt(k,362) + het_rates(k,129) )
         mat(k,1531) = rxt(k,7)
         mat(k,245) = rxt(k,300)
         mat(k,1784) = rxt(k,310)
         mat(k,167) = rxt(k,361)
         mat(k,160) = rxt(k,363)
         mat(k,109) = -( het_rates(k,131) )
         mat(k,324) = -( het_rates(k,91) )
         mat(k,177) = -( rxt(k,93) + het_rates(k,92) )
         mat(k,210) = -( het_rates(k,93) )
         mat(k,137) = rxt(k,341)
         mat(k,136) = -( rxt(k,341) + het_rates(k,94) )
         mat(k,1362) = rxt(k,343)
         mat(k,1392) = -( rxt(k,343) + het_rates(k,95) )
         mat(k,1476) = rxt(k,345)
         mat(k,1478) = -( rxt(k,345) + het_rates(k,96) )
         mat(k,986) = rxt(k,347)
         mat(k,975) = -( rxt(k,347) + het_rates(k,97) )
         mat(k,77) = -( het_rates(k,98) )
         mat(k,42) = -( het_rates(k,99) )
         mat(k,46) = -( het_rates(k,100) )
         mat(k,1251) = -( het_rates(k,101) )
         mat(k,893) = -( het_rates(k,102) )
         mat(k,83) = -( het_rates(k,103) )
         mat(k,161) = -( rxt(k,361) + het_rates(k,104) )
         mat(k,153) = -( rxt(k,363) + het_rates(k,105) )
         mat(k,1864) = rxt(k,362)
         mat(k,1288) = -( rxt(k,94) + rxt(k,367) + het_rates(k,106) )
         mat(k,468) = rxt(k,100)
         mat(k,859) = rxt(k,102)
         mat(k,1818) = -( rxt(k,95) + rxt(k,371) + het_rates(k,107) )
         mat(k,128) = rxt(k,101)
         mat(k,703) = rxt(k,103)
         mat(k,459) = -( rxt(k,96) + rxt(k,100) + het_rates(k,108) )
         mat(k,122) = -( rxt(k,101) + het_rates(k,109) )
         mat(k,849) = -( rxt(k,102) + het_rates(k,111) )
         mat(k,536) = rxt(k,104) + rxt(k,452)
         mat(k,677) = -( rxt(k,103) + het_rates(k,112) )
         mat(k,770) = -( rxt(k,98) + het_rates(k,113) )
         mat(k,414) = rxt(k,457)
         mat(k,710) = -( rxt(k,99) + het_rates(k,114) )
         mat(k,555) = rxt(k,449)
         mat(k,497) = rxt(k,462)
         mat(k,1603) = -( rxt(k,97) + het_rates(k,110) )
         mat(k,336) = -( het_rates(k,122) )
         mat(k,532) = -( rxt(k,104) + rxt(k,452) + het_rates(k,115) )
         mat(k,368) = rxt(k,455)
         mat(k,367) = -( rxt(k,455) + het_rates(k,116) )
         mat(k,413) = -( rxt(k,457) + het_rates(k,117) )
         mat(k,554) = -( rxt(k,449) + het_rates(k,118) )
         mat(k,513) = rxt(k,460)
         mat(k,512) = -( rxt(k,460) + het_rates(k,119) )
         mat(k,496) = -( rxt(k,462) + het_rates(k,120) )
         mat(k,427) = -( het_rates(k,121) )
         mat(k,655) = -( het_rates(k,123) )
         mat(k,399) = rxt(k,445)
         mat(k,383) = rxt(k,446)
         mat(k,265) = -( het_rates(k,124) )
         mat(k,398) = -( rxt(k,445) + het_rates(k,125) )
         mat(k,382) = -( rxt(k,446) + het_rates(k,126) )
         mat(k,1161) = -( het_rates(k,130) )
         mat(k,1514) = rxt(k,7)
         mat(k,1121) = rxt(k,67) + rxt(k,68) + rxt(k,69) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82)
         mat(k,314) = rxt(k,70)
         mat(k,1207) = rxt(k,71) + rxt(k,73) + rxt(k,74) + rxt(k,75) + rxt(k,83) &
                      + rxt(k,85) + rxt(k,86) + rxt(k,87)
         mat(k,1285) = rxt(k,94) + rxt(k,367)
         mat(k,1803) = rxt(k,95) + rxt(k,371)
         mat(k,466) = rxt(k,96)
         mat(k,1593) = rxt(k,97)
         mat(k,778) = rxt(k,98)
         mat(k,718) = rxt(k,99)
         mat(k,24) = -( rxt(k,55) + het_rates(k,53) )
         mat(k,620) = rxt(k,126)*y(k,38) + rxt(k,127)*y(k,39) &
                      + 2.000_r8*rxt(k,128)*y(k,47) + 2.000_r8*rxt(k,129)*y(k,48) &
                      + rxt(k,130)*y(k,40) + rxt(k,132)*y(k,46) + rxt(k,135)*y(k,44) &
                      + rxt(k,136)*y(k,43) + rxt(k,137)*y(k,49) &
                      + 2.000_r8*rxt(k,138)*y(k,50)
         mat(k,793) = rxt(k,243)*y(k,40) + rxt(k,247)*y(k,46)
         mat(k,32) = -( rxt(k,56) + het_rates(k,54) )
         mat(k,621) = rxt(k,125)*y(k,37) + rxt(k,127)*y(k,39) + rxt(k,131)*y(k,45)
         mat(k,794) = rxt(k,246)*y(k,45)
         mat(k,39) = -( rxt(k,57) + het_rates(k,55) )
         mat(k,200) = rxt(k,238)*y(k,12)
         mat(k,201) = -( rxt(k,238)*y(k,12) + het_rates(k,56) )
         mat(k,25) = 2.000_r8*rxt(k,55)
         mat(k,33) = rxt(k,56)
         mat(k,40) = rxt(k,57)
         mat(k,623) = rxt(k,129)*y(k,48) + rxt(k,136)*y(k,43)
         mat(k,69) = -( rxt(k,64) + het_rates(k,58) )
         mat(k,168) = -( het_rates(k,59) )
         mat(k,70) = rxt(k,64)
         mat(k,351) = rxt(k,65)
         mat(k,353) = -( rxt(k,65) + het_rates(k,60) )
         mat(k,248) = rxt(k,62)
         mat(k,247) = -( rxt(k,62) + het_rates(k,61) )
         mat(k,52) = rxt(k,63)
         mat(k,51) = -( rxt(k,63) + het_rates(k,62) )
         mat(k,28) = rxt(k,61)
         mat(k,27) = -( rxt(k,61) + het_rates(k,63) )
         mat(k,55) = -( het_rates(k,64) )
         mat(k,1) = -( het_rates(k,65) )
         mat(k,2) = -( het_rates(k,66) )
         mat(k,3) = -( het_rates(k,67) )
         mat(k,4) = -( het_rates(k,68) )
         mat(k,5) = -( het_rates(k,69) )
         mat(k,6) = -( het_rates(k,70) )
         mat(k,7) = -( het_rates(k,71) )
         mat(k,8) = -( het_rates(k,72) )
         mat(k,9) = -( het_rates(k,73) )
         mat(k,10) = -( het_rates(k,74) )
         mat(k,11) = -( het_rates(k,75) )
         mat(k,12) = -( het_rates(k,76) )
         mat(k,13) = -( het_rates(k,77) )
         mat(k,14) = -( het_rates(k,78) )
         mat(k,15) = -( het_rates(k,79) )
         mat(k,16) = -( het_rates(k,80) )
         mat(k,17) = -( het_rates(k,81) )
         mat(k,18) = -( het_rates(k,82) )
         mat(k,19) = -( het_rates(k,83) )
         mat(k,20) = -( het_rates(k,84) )
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
