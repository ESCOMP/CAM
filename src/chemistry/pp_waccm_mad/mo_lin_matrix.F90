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
         mat(k,1543) = -( rxt(k,3) + rxt(k,4) + het_rates(k,1) )
         mat(k,425) = rxt(k,92)
         mat(k,1817) = -( rxt(k,63) + rxt(k,64) + rxt(k,65) + rxt(k,76) + rxt(k,77) &
                      + rxt(k,78) + het_rates(k,2) )
         mat(k,1672) = rxt(k,1) + 2.000_r8*rxt(k,2) + rxt(k,69) + rxt(k,70) + rxt(k,71) &
                      + 2.000_r8*rxt(k,74) + rxt(k,81) + rxt(k,82) + rxt(k,83) &
                      + 2.000_r8*rxt(k,86)
         mat(k,1549) = rxt(k,4)
         mat(k,1173) = rxt(k,6)
         mat(k,1761) = rxt(k,8)
         mat(k,249) = rxt(k,10)
         mat(k,1629) = rxt(k,12)
         mat(k,1054) = rxt(k,21)
         mat(k,1716) = rxt(k,24)
         mat(k,6) = rxt(k,25)
         mat(k,308) = rxt(k,32)
         mat(k,1511) = rxt(k,58) + rxt(k,88)
         mat(k,958) = rxt(k,90) + rxt(k,358)
         mat(k,584) = rxt(k,114)
         mat(k,64) = rxt(k,293)
         mat(k,70) = rxt(k,300)
         mat(k,563) = -( rxt(k,114) + rxt(k,118)*y(k,7) + rxt(k,119)*y(k,7) &
                      + rxt(k,121)*y(k,45) + rxt(k,122)*y(k,46) + rxt(k,123)*y(k,47) &
                      + rxt(k,124)*y(k,55) + rxt(k,125)*y(k,56) + rxt(k,126)*y(k,48) &
                      + rxt(k,127)*y(k,53) + rxt(k,128)*y(k,54) + rxt(k,129)*y(k,49) &
                      + rxt(k,130)*y(k,44) + rxt(k,131)*y(k,52) + rxt(k,132)*y(k,51) &
                      + rxt(k,133)*y(k,57) + rxt(k,134)*y(k,58) + rxt(k,135)*y(k,59) &
                      + rxt(k,136)*y(k,60) + rxt(k,139)*y(k,15) + rxt(k,140)*y(k,15) &
                      + rxt(k,141)*y(k,15) + het_rates(k,3) )
         mat(k,1642) = rxt(k,1)
         mat(k,1519) = rxt(k,3)
         mat(k,1024) = rxt(k,20)
         mat(k,1669) = -( rxt(k,1) + rxt(k,2) + rxt(k,67) + rxt(k,69) + rxt(k,70) &
                      + rxt(k,71) + rxt(k,74) + rxt(k,79) + rxt(k,81) + rxt(k,82) &
                      + rxt(k,83) + rxt(k,86) + het_rates(k,4) )
         mat(k,1546) = rxt(k,4)
         mat(k,1626) = rxt(k,13)
         mat(k,832) = rxt(k,91) + rxt(k,362)
         mat(k,427) = rxt(k,96)
         mat(k,103) = rxt(k,97)
         mat(k,32) = rxt(k,109)
         mat(k,227) = rxt(k,112) + rxt(k,113)
         mat(k,581) = rxt(k,119)*y(k,7)
         mat(k,30) = -( rxt(k,106) + rxt(k,109) + het_rates(k,5) )
         mat(k,220) = -( rxt(k,112) + rxt(k,113) + het_rates(k,6) )
         mat(k,1513) = rxt(k,3)
         mat(k,31) = rxt(k,106)
         mat(k,532) = -( het_rates(k,20) )
         mat(k,403) = rxt(k,18)
         mat(k,1023) = rxt(k,20)
         mat(k,562) = rxt(k,141)*y(k,15)
         mat(k,148) = -( het_rates(k,19) )
         mat(k,399) = rxt(k,17) + rxt(k,18)
         mat(k,1475) = rxt(k,58) + rxt(k,88)
         mat(k,838) = rxt(k,236)*y(k,43)
         mat(k,1504) = -( rxt(k,58) + rxt(k,88) + het_rates(k,65) )
         mat(k,792) = rxt(k,98)
         mat(k,631) = rxt(k,99)
         mat(k,135) = rxt(k,352)
         mat(k,264) = -( rxt(k,66) + het_rates(k,8) )
         mat(k,1134) = rxt(k,6)
         mat(k,199) = rxt(k,279)
         mat(k,1158) = -( rxt(k,6) + rxt(k,7) + het_rates(k,9) )
         mat(k,1746) = rxt(k,8) + .500_r8*rxt(k,259)
         mat(k,243) = rxt(k,10)
         mat(k,1614) = rxt(k,13)
         mat(k,181) = rxt(k,62)
         mat(k,908) = rxt(k,289)
         mat(k,63) = rxt(k,294)
         mat(k,572) = 2.000_r8*rxt(k,118)*y(k,7)
         mat(k,1760) = -( rxt(k,8) + rxt(k,259) + het_rates(k,10) )
         mat(k,248) = rxt(k,9) + rxt(k,179)
         mat(k,1469) = rxt(k,11)
         mat(k,1628) = rxt(k,12)
         mat(k,57) = rxt(k,15) + rxt(k,188)
         mat(k,218) = rxt(k,30)
         mat(k,82) = rxt(k,36)
         mat(k,723) = rxt(k,94)
         mat(k,745) = -( rxt(k,237)*y(k,43) + rxt(k,238)*y(k,50) + rxt(k,239)*y(k,48) &
                      + rxt(k,240)*y(k,44) + rxt(k,242)*y(k,53) + rxt(k,243)*y(k,54) &
                      + rxt(k,244)*y(k,60) + rxt(k,245)*y(k,59) + rxt(k,248)*y(k,15) &
                 + het_rates(k,22) )
         mat(k,1445) = rxt(k,11)
         mat(k,53) = rxt(k,14)
         mat(k,36) = rxt(k,16)
         mat(k,1029) = rxt(k,19)
         mat(k,86) = 2.000_r8*rxt(k,22)
         mat(k,191) = rxt(k,27)
         mat(k,157) = rxt(k,33)
         mat(k,179) = rxt(k,62)
         mat(k,1105) = rxt(k,93)
         mat(k,1736) = .500_r8*rxt(k,259)
         mat(k,565) = rxt(k,139)*y(k,15)
         mat(k,1625) = -( rxt(k,12) + rxt(k,13) + rxt(k,258) + het_rates(k,11) )
         mat(k,247) = rxt(k,9) + rxt(k,10) + rxt(k,179)
         mat(k,55) = rxt(k,14)
         mat(k,216) = rxt(k,29)
         mat(k,81) = rxt(k,35)
         mat(k,659) = rxt(k,95)
         mat(k,178) = -( rxt(k,62) + het_rates(k,25) )
         mat(k,1462) = -( rxt(k,11) + het_rates(k,12) )
         mat(k,246) = 2.000_r8*rxt(k,257) + 2.000_r8*rxt(k,261) + 2.000_r8*rxt(k,267) &
                      + 2.000_r8*rxt(k,272)
         mat(k,1621) = rxt(k,258)
         mat(k,1753) = .500_r8*rxt(k,259)
         mat(k,215) = rxt(k,262) + rxt(k,268) + rxt(k,273)
         mat(k,80) = rxt(k,263) + rxt(k,271) + rxt(k,274)
         mat(k,441) = rxt(k,453)
         mat(k,51) = -( rxt(k,14) + rxt(k,15) + rxt(k,188) + het_rates(k,13) )
         mat(k,238) = -( rxt(k,9) + rxt(k,10) + rxt(k,179) + rxt(k,257) + rxt(k,261) &
                      + rxt(k,267) + rxt(k,272) + het_rates(k,14) )
         mat(k,251) = -( het_rates(k,16) )
         mat(k,558) = rxt(k,139)*y(k,15)
         mat(k,842) = rxt(k,195)*y(k,15)
         mat(k,164) = rxt(k,234)*y(k,15)
         mat(k,737) = rxt(k,248)*y(k,15)
         mat(k,33) = -( rxt(k,16) + het_rates(k,17) )
         mat(k,401) = -( rxt(k,17) + rxt(k,18) + het_rates(k,18) )
         mat(k,35) = rxt(k,16)
         mat(k,560) = rxt(k,140)*y(k,15) + rxt(k,141)*y(k,15)
         mat(k,1376) = -( het_rates(k,21) )
         mat(k,38) = rxt(k,16)
         mat(k,410) = 2.000_r8*rxt(k,17)
         mat(k,1044) = rxt(k,19) + 2.000_r8*rxt(k,21)
         mat(k,1297) = rxt(k,28)
         mat(k,176) = rxt(k,34)
         mat(k,21) = rxt(k,57)
         mat(k,575) = rxt(k,140)*y(k,15)
         mat(k,676) = -( rxt(k,260) + het_rates(k,23) )
         mat(k,52) = rxt(k,15) + rxt(k,188)
         mat(k,564) = rxt(k,140)*y(k,15)
         mat(k,850) = rxt(k,236)*y(k,43) + rxt(k,241)*y(k,44)
         mat(k,744) = rxt(k,237)*y(k,43) + rxt(k,240)*y(k,44)
         mat(k,84) = -( rxt(k,22) + het_rates(k,24) )
         mat(k,665) = .500_r8*rxt(k,260)
         mat(k,1036) = -( rxt(k,19) + rxt(k,20) + rxt(k,21) + het_rates(k,111) )
         mat(k,143) = rxt(k,89)
         mat(k,473) = rxt(k,100) + rxt(k,443)
         mat(k,107) = rxt(k,332)
         mat(k,1198) = rxt(k,334)
         mat(k,1332) = rxt(k,336)
         mat(k,1411) = rxt(k,338)
         mat(k,359) = rxt(k,436)
         mat(k,492) = rxt(k,440)
         mat(k,327) = rxt(k,446)
         mat(k,373) = rxt(k,448)
         mat(k,451) = rxt(k,451)
         mat(k,751) = rxt(k,237)*y(k,43) + rxt(k,238)*y(k,50) + rxt(k,239)*y(k,48) &
                      + rxt(k,240)*y(k,44) + rxt(k,244)*y(k,60) + rxt(k,248)*y(k,15)
         mat(k,855) = -( rxt(k,195)*y(k,15) + rxt(k,236)*y(k,43) + rxt(k,241)*y(k,44) &
                      + rxt(k,246)*y(k,60) + rxt(k,247)*y(k,59) + het_rates(k,29) )
         mat(k,11) = 2.000_r8*rxt(k,23)
         mat(k,1694) = rxt(k,24)
         mat(k,2) = 2.000_r8*rxt(k,26)
         mat(k,192) = rxt(k,27)
         mat(k,1285) = rxt(k,28)
         mat(k,212) = rxt(k,29)
         mat(k,18) = rxt(k,31)
         mat(k,15) = rxt(k,56)
         mat(k,567) = 2.000_r8*rxt(k,121)*y(k,45) + 2.000_r8*rxt(k,122)*y(k,46) &
                      + 2.000_r8*rxt(k,123)*y(k,47) + 2.000_r8*rxt(k,124)*y(k,55) &
                      + rxt(k,125)*y(k,56) + rxt(k,126)*y(k,48) + rxt(k,127)*y(k,53) &
                      + rxt(k,128)*y(k,54) + 4.000_r8*rxt(k,129)*y(k,49) &
                      + rxt(k,131)*y(k,52)
         mat(k,747) = rxt(k,237)*y(k,43) + 3.000_r8*rxt(k,238)*y(k,50) &
                      + rxt(k,239)*y(k,48) + rxt(k,242)*y(k,53) + rxt(k,243)*y(k,54)
         mat(k,10) = -( rxt(k,23) + het_rates(k,30) )
         mat(k,1714) = -( rxt(k,24) + het_rates(k,31) )
         mat(k,5) = rxt(k,25)
         mat(k,217) = rxt(k,30)
         mat(k,3) = 2.000_r8*rxt(k,207)
         mat(k,4) = -( rxt(k,25) + het_rates(k,32) )
         mat(k,1) = -( rxt(k,26) + rxt(k,207) + het_rates(k,33) )
         mat(k,1295) = -( rxt(k,28) + het_rates(k,34) )
         mat(k,348) = rxt(k,437)
         mat(k,865) = rxt(k,195)*y(k,15) + 2.000_r8*rxt(k,236)*y(k,43) &
                      + rxt(k,241)*y(k,44) + rxt(k,246)*y(k,60) + rxt(k,247)*y(k,59)
         mat(k,190) = -( rxt(k,27) + het_rates(k,35) )
         mat(k,209) = rxt(k,262) + rxt(k,268) + rxt(k,273)
         mat(k,210) = -( rxt(k,29) + rxt(k,30) + rxt(k,262) + rxt(k,268) + rxt(k,273) &
                 + het_rates(k,36) )
         mat(k,16) = -( rxt(k,31) + het_rates(k,37) )
         mat(k,509) = -( het_rates(k,38) )
         mat(k,17) = rxt(k,31)
         mat(k,296) = rxt(k,32)
         mat(k,156) = rxt(k,33)
         mat(k,172) = rxt(k,34)
         mat(k,79) = rxt(k,35)
         mat(k,561) = rxt(k,130)*y(k,44) + rxt(k,131)*y(k,52) + rxt(k,132)*y(k,51) &
                      + 2.000_r8*rxt(k,133)*y(k,57) + 2.000_r8*rxt(k,134)*y(k,58) &
                      + 3.000_r8*rxt(k,135)*y(k,59) + 2.000_r8*rxt(k,136)*y(k,60)
         mat(k,741) = rxt(k,240)*y(k,44) + 2.000_r8*rxt(k,244)*y(k,60) &
                      + 3.000_r8*rxt(k,245)*y(k,59)
         mat(k,844) = rxt(k,241)*y(k,44) + 2.000_r8*rxt(k,246)*y(k,60) &
                      + 3.000_r8*rxt(k,247)*y(k,59)
         mat(k,295) = -( rxt(k,32) + het_rates(k,39) )
         mat(k,78) = rxt(k,36)
         mat(k,171) = -( rxt(k,34) + het_rates(k,40) )
         mat(k,154) = -( rxt(k,33) + het_rates(k,41) )
         mat(k,77) = rxt(k,263) + rxt(k,271) + rxt(k,274)
         mat(k,76) = -( rxt(k,35) + rxt(k,36) + rxt(k,263) + rxt(k,271) + rxt(k,274) &
                 + het_rates(k,42) )
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
         mat(k,112) = -( het_rates(k,66) )
         mat(k,61) = rxt(k,293)
         mat(k,66) = rxt(k,300)
         mat(k,902) = -( rxt(k,289) + het_rates(k,67) )
         mat(k,1651) = rxt(k,67) + rxt(k,79)
         mat(k,142) = rxt(k,89)
         mat(k,91) = -( het_rates(k,68) )
         mat(k,262) = rxt(k,66)
         mat(k,60) = rxt(k,294)
         mat(k,198) = -( rxt(k,279) + het_rates(k,69) )
         mat(k,1775) = rxt(k,63) + rxt(k,76)
         mat(k,1637) = rxt(k,69) + rxt(k,81)
         mat(k,62) = rxt(k,302)
         mat(k,67) = rxt(k,304)
         mat(k,1249) = -( rxt(k,353) + het_rates(k,70) )
         mat(k,1160) = rxt(k,7)
         mat(k,203) = rxt(k,279)
         mat(k,910) = rxt(k,289)
         mat(k,134) = rxt(k,352)
         mat(k,127) = rxt(k,354)
         mat(k,71) = -( het_rates(k,72) )
         mat(k,1587) = -( het_rates(k,71) )
         mat(k,1168) = rxt(k,7)
         mat(k,1812) = rxt(k,63) + rxt(k,64) + rxt(k,65) + rxt(k,76) + rxt(k,77) &
                      + rxt(k,78)
         mat(k,272) = rxt(k,66)
         mat(k,1667) = rxt(k,67) + rxt(k,69) + rxt(k,70) + rxt(k,71) + rxt(k,79) &
                      + rxt(k,81) + rxt(k,82) + rxt(k,83)
         mat(k,953) = rxt(k,90) + rxt(k,358)
         mat(k,830) = rxt(k,91) + rxt(k,362)
         mat(k,426) = rxt(k,92)
         mat(k,1125) = rxt(k,93)
         mat(k,719) = rxt(k,94)
         mat(k,658) = rxt(k,95)
         mat(k,58) = -( rxt(k,293) + rxt(k,294) + rxt(k,302) + rxt(k,303) &
                 + het_rates(k,73) )
         mat(k,1763) = rxt(k,65) + rxt(k,78)
         mat(k,1631) = rxt(k,71) + rxt(k,83)
         mat(k,65) = -( rxt(k,300) + rxt(k,304) + het_rates(k,74) )
         mat(k,1764) = rxt(k,64) + rxt(k,77)
         mat(k,1632) = rxt(k,70) + rxt(k,82)
         mat(k,59) = rxt(k,303)
         mat(k,279) = -( het_rates(k,75) )
         mat(k,139) = -( rxt(k,89) + het_rates(k,76) )
         mat(k,184) = -( het_rates(k,77) )
         mat(k,106) = rxt(k,332)
         mat(k,105) = -( rxt(k,332) + het_rates(k,78) )
         mat(k,1174) = rxt(k,334)
         mat(k,1202) = -( rxt(k,334) + het_rates(k,79) )
         mat(k,1336) = rxt(k,336)
         mat(k,1339) = -( rxt(k,336) + het_rates(k,80) )
         mat(k,1418) = rxt(k,338)
         mat(k,1420) = -( rxt(k,338) + het_rates(k,81) )
         mat(k,39) = -( het_rates(k,82) )
         mat(k,22) = -( het_rates(k,83) )
         mat(k,26) = -( het_rates(k,84) )
         mat(k,981) = -( het_rates(k,85) )
         mat(k,1080) = -( het_rates(k,86) )
         mat(k,45) = -( het_rates(k,87) )
         mat(k,130) = -( rxt(k,352) + het_rates(k,88) )
         mat(k,122) = -( rxt(k,354) + het_rates(k,89) )
         mat(k,1218) = rxt(k,353)
         mat(k,938) = -( rxt(k,90) + rxt(k,358) + het_rates(k,90) )
         mat(k,420) = rxt(k,96)
         mat(k,779) = rxt(k,98)
         mat(k,812) = -( rxt(k,91) + rxt(k,362) + het_rates(k,91) )
         mat(k,101) = rxt(k,97)
         mat(k,615) = rxt(k,99)
         mat(k,415) = -( rxt(k,92) + rxt(k,96) + het_rates(k,92) )
         mat(k,98) = -( rxt(k,97) + het_rates(k,93) )
         mat(k,775) = -( rxt(k,98) + het_rates(k,95) )
         mat(k,470) = rxt(k,100) + rxt(k,443)
         mat(k,609) = -( rxt(k,99) + het_rates(k,96) )
         mat(k,702) = -( rxt(k,94) + het_rates(k,97) )
         mat(k,370) = rxt(k,448)
         mat(k,642) = -( rxt(k,95) + het_rates(k,98) )
         mat(k,489) = rxt(k,440)
         mat(k,431) = rxt(k,453)
         mat(k,1114) = -( rxt(k,93) + het_rates(k,94) )
         mat(k,309) = -( het_rates(k,106) )
         mat(k,466) = -( rxt(k,100) + rxt(k,443) + het_rates(k,99) )
         mat(k,324) = rxt(k,446)
         mat(k,323) = -( rxt(k,446) + het_rates(k,100) )
         mat(k,369) = -( rxt(k,448) + het_rates(k,101) )
         mat(k,488) = -( rxt(k,440) + het_rates(k,102) )
         mat(k,447) = rxt(k,451)
         mat(k,446) = -( rxt(k,451) + het_rates(k,103) )
         mat(k,430) = -( rxt(k,453) + het_rates(k,104) )
         mat(k,383) = -( het_rates(k,105) )
         mat(k,587) = -( het_rates(k,107) )
         mat(k,355) = rxt(k,436)
         mat(k,339) = rxt(k,437)
         mat(k,228) = -( het_rates(k,108) )
         mat(k,354) = -( rxt(k,436) + het_rates(k,109) )
         mat(k,338) = -( rxt(k,437) + het_rates(k,110) )
         mat(k,7) = -( rxt(k,55) + het_rates(k,61) )
         mat(k,552) = rxt(k,122)*y(k,46) + rxt(k,123)*y(k,47) &
                      + 2.000_r8*rxt(k,124)*y(k,55) + 2.000_r8*rxt(k,125)*y(k,56) &
                      + rxt(k,126)*y(k,48) + rxt(k,128)*y(k,54) + rxt(k,131)*y(k,52) &
                      + rxt(k,132)*y(k,51) + rxt(k,133)*y(k,57) &
                      + 2.000_r8*rxt(k,134)*y(k,58)
         mat(k,725) = rxt(k,239)*y(k,48) + rxt(k,243)*y(k,54)
         mat(k,12) = -( rxt(k,56) + het_rates(k,62) )
         mat(k,553) = rxt(k,121)*y(k,45) + rxt(k,123)*y(k,47) + rxt(k,127)*y(k,53)
         mat(k,726) = rxt(k,242)*y(k,53)
         mat(k,19) = -( rxt(k,57) + het_rates(k,63) )
         mat(k,162) = rxt(k,234)*y(k,15)
         mat(k,163) = -( rxt(k,234)*y(k,15) + het_rates(k,64) )
         mat(k,8) = 2.000_r8*rxt(k,55)
         mat(k,13) = rxt(k,56)
         mat(k,20) = rxt(k,57)
         mat(k,555) = rxt(k,125)*y(k,56) + rxt(k,132)*y(k,51)
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
