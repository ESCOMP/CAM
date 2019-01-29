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
         mat(k,1) = -( het_rates(k,1) )
         mat(k,2) = -( het_rates(k,2) )
         mat(k,39) = -( rxt(k,27) + het_rates(k,3) )
         mat(k,584) = -( rxt(k,28) + het_rates(k,4) )
         mat(k,141) = rxt(k,29)
         mat(k,138) = -( rxt(k,29) + rxt(k,30) + rxt(k,565) + rxt(k,568) + rxt(k,573) &
                 + het_rates(k,5) )
         mat(k,606) = -( rxt(k,21) + rxt(k,22) + het_rates(k,16) )
         mat(k,79) = rxt(k,23)
         mat(k,636) = rxt(k,539)*y(k,22) + rxt(k,540)*y(k,22)
         mat(k,324) = -( het_rates(k,20) )
         mat(k,1450) = rxt(k,451)*y(k,22)
         mat(k,199) = rxt(k,507)*y(k,22)
         mat(k,786) = rxt(k,536)*y(k,22)
         mat(k,631) = rxt(k,538)*y(k,22)
         mat(k,77) = -( rxt(k,23) + het_rates(k,21) )
         mat(k,30) = -( rxt(k,44) + het_rates(k,24) )
         mat(k,21) = -( rxt(k,45) + rxt(k,485) + het_rates(k,25) )
         mat(k,1761) = -( rxt(k,46) + het_rates(k,26) )
         mat(k,285) = rxt(k,48)
         mat(k,68) = rxt(k,60)
         mat(k,23) = 2.000_r8*rxt(k,485)
         mat(k,277) = -( rxt(k,47) + rxt(k,48) + rxt(k,567) + rxt(k,572) + rxt(k,578) &
                 + het_rates(k,27) )
         mat(k,134) = -( het_rates(k,29) )
         mat(k,601) = rxt(k,21) + rxt(k,22)
         mat(k,70) = rxt(k,101)
         mat(k,1445) = rxt(k,518)*y(k,19)
         mat(k,261) = rxt(k,595)*y(k,30)
         mat(k,27) = -( rxt(k,49) + het_rates(k,31) )
         mat(k,625) = rxt(k,476)*y(k,8) + rxt(k,478)*y(k,11) + 2.000_r8*rxt(k,479)*y(k,12) &
                      + 2.000_r8*rxt(k,480)*y(k,13) + rxt(k,481)*y(k,14) &
                      + rxt(k,502)*y(k,9) + 2.000_r8*rxt(k,504)*y(k,40) &
                      + rxt(k,528)*y(k,45) + rxt(k,529)*y(k,46)
         mat(k,768) = rxt(k,523)*y(k,45) + rxt(k,524)*y(k,46)
         mat(k,32) = -( rxt(k,50) + het_rates(k,32) )
         mat(k,626) = rxt(k,477)*y(k,10) + rxt(k,478)*y(k,11) + rxt(k,527)*y(k,44)
         mat(k,769) = rxt(k,522)*y(k,44)
         mat(k,58) = -( het_rates(k,33) )
         mat(k,3) = -( het_rates(k,34) )
         mat(k,4) = -( het_rates(k,35) )
         mat(k,5) = -( het_rates(k,36) )
         mat(k,198) = -( rxt(k,507)*y(k,22) + het_rates(k,37) )
         mat(k,28) = 2.000_r8*rxt(k,49)
         mat(k,33) = rxt(k,50)
         mat(k,37) = rxt(k,57)
         mat(k,628) = rxt(k,480)*y(k,13) + rxt(k,502)*y(k,9)
         mat(k,1797) = -( het_rates(k,38) )
         mat(k,1855) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,623) = 2.000_r8*rxt(k,21)
         mat(k,81) = rxt(k,23)
         mat(k,244) = rxt(k,52)
         mat(k,1175) = rxt(k,56)
         mat(k,38) = rxt(k,57)
         mat(k,656) = rxt(k,539)*y(k,22)
         mat(k,484) = -( het_rates(k,39) )
         mat(k,1821) = rxt(k,1)
         mat(k,604) = rxt(k,22)
         mat(k,634) = rxt(k,540)*y(k,22)
         mat(k,146) = -( rxt(k,4) + het_rates(k,41) )
         mat(k,736) = .500_r8*rxt(k,559)
         mat(k,24) = -( rxt(k,100) + het_rates(k,42) )
         mat(k,238) = -( rxt(k,52) + het_rates(k,43) )
         mat(k,1160) = -( rxt(k,56) + het_rates(k,47) )
         mat(k,399) = rxt(k,386)
         mat(k,1469) = rxt(k,451)*y(k,22) + rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) &
                      + 2.000_r8*rxt(k,518)*y(k,19) + rxt(k,520)*y(k,23)
         mat(k,36) = -( rxt(k,57) + het_rates(k,48) )
         mat(k,197) = rxt(k,507)*y(k,22)
         mat(k,1292) = -( rxt(k,9) + het_rates(k,49) )
         mat(k,511) = rxt(k,275)
         mat(k,303) = 2.000_r8*rxt(k,560) + 2.000_r8*rxt(k,563) + 2.000_r8*rxt(k,566) &
                      + 2.000_r8*rxt(k,577)
         mat(k,997) = .500_r8*rxt(k,561)
         mat(k,1119) = rxt(k,562)
         mat(k,144) = rxt(k,565) + rxt(k,568) + rxt(k,573)
         mat(k,282) = rxt(k,567) + rxt(k,572) + rxt(k,578)
         mat(k,101) = -( rxt(k,10) + rxt(k,11) + rxt(k,448) + het_rates(k,50) )
         mat(k,207) = -( rxt(k,58) + het_rates(k,51) )
         mat(k,139) = rxt(k,565) + rxt(k,568) + rxt(k,573)
         mat(k,231) = -( rxt(k,59) + het_rates(k,52) )
         mat(k,276) = rxt(k,567) + rxt(k,572) + rxt(k,578)
         mat(k,175) = -( rxt(k,12) + het_rates(k,53) )
         mat(k,311) = -( rxt(k,66) + het_rates(k,54) )
         mat(k,1488) = rxt(k,17)
         mat(k,265) = rxt(k,596)
         mat(k,297) = -( rxt(k,14) + rxt(k,15) + rxt(k,449) + rxt(k,560) + rxt(k,563) &
                      + rxt(k,566) + rxt(k,577) + het_rates(k,56) )
         mat(k,6) = -( het_rates(k,57) )
         mat(k,7) = -( het_rates(k,58) )
         mat(k,8) = -( het_rates(k,59) )
         mat(k,1518) = -( rxt(k,16) + rxt(k,17) + het_rates(k,60) )
         mat(k,178) = rxt(k,12)
         mat(k,305) = rxt(k,15)
         mat(k,1002) = rxt(k,18) + .500_r8*rxt(k,561)
         mat(k,1124) = rxt(k,20)
         mat(k,1376) = rxt(k,593)
         mat(k,127) = rxt(k,606)
         mat(k,651) = 2.000_r8*rxt(k,442)*y(k,55)
         mat(k,990) = -( rxt(k,18) + rxt(k,561) + het_rates(k,61) )
         mat(k,1285) = rxt(k,9)
         mat(k,104) = rxt(k,11) + rxt(k,448)
         mat(k,301) = rxt(k,14) + rxt(k,449)
         mat(k,1112) = rxt(k,19)
         mat(k,142) = rxt(k,29)
         mat(k,279) = rxt(k,48)
         mat(k,855) = rxt(k,75)
         mat(k,1115) = -( rxt(k,19) + rxt(k,20) + rxt(k,562) + het_rates(k,62) )
         mat(k,105) = rxt(k,10)
         mat(k,302) = rxt(k,14) + rxt(k,15) + rxt(k,449)
         mat(k,143) = rxt(k,30)
         mat(k,280) = rxt(k,47)
         mat(k,720) = rxt(k,76)
         mat(k,9) = -( het_rates(k,63) )
         mat(k,10) = -( het_rates(k,64) )
         mat(k,11) = -( het_rates(k,65) )
         mat(k,12) = -( het_rates(k,66) )
         mat(k,1433) = -( rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82) + het_rates(k,67) )
         mat(k,1846) = rxt(k,2)
         mat(k,1256) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,83) + rxt(k,85) + rxt(k,87) &
                      + 2.000_r8*rxt(k,88) + 2.000_r8*rxt(k,89) + rxt(k,90) + rxt(k,91) &
                      + rxt(k,92)
         mat(k,1207) = rxt(k,8)
         mat(k,304) = rxt(k,15)
         mat(k,1516) = rxt(k,17)
         mat(k,1000) = rxt(k,18)
         mat(k,1122) = rxt(k,19)
         mat(k,595) = rxt(k,28)
         mat(k,1753) = rxt(k,46)
         mat(k,67) = rxt(k,60)
         mat(k,1667) = rxt(k,99) + rxt(k,358)
         mat(k,387) = rxt(k,102)
         mat(k,259) = rxt(k,103)
         mat(k,48) = rxt(k,104)
         mat(k,649) = rxt(k,391)
         mat(k,133) = rxt(k,600)
         mat(k,126) = rxt(k,605)
         mat(k,1252) = -( rxt(k,5) + rxt(k,6) + rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86) + rxt(k,87) + rxt(k,88) + rxt(k,89) + rxt(k,90) &
                      + rxt(k,91) + rxt(k,92) + het_rates(k,68) )
         mat(k,1203) = rxt(k,8)
         mat(k,1118) = rxt(k,20)
         mat(k,1326) = rxt(k,93) + rxt(k,137)
         mat(k,444) = rxt(k,95) + rxt(k,330)*y(k,30)
         mat(k,98) = rxt(k,97) + rxt(k,335)*y(k,30)
         mat(k,249) = rxt(k,387) + rxt(k,395)
         mat(k,44) = rxt(k,388)
         mat(k,645) = rxt(k,443)*y(k,55)
         mat(k,1202) = -( rxt(k,7) + rxt(k,8) + het_rates(k,69) )
         mat(k,443) = rxt(k,96)
         mat(k,64) = -( rxt(k,60) + het_rates(k,70) )
         mat(k,69) = -( rxt(k,101) + het_rates(k,71) )
         mat(k,13) = -( het_rates(k,72) )
         mat(k,14) = -( het_rates(k,73) )
         mat(k,168) = -( het_rates(k,74) )
         mat(k,71) = rxt(k,101)
         mat(k,377) = rxt(k,102)
         mat(k,379) = -( rxt(k,102) + het_rates(k,76) )
         mat(k,256) = rxt(k,103)
         mat(k,255) = -( rxt(k,103) + het_rates(k,77) )
         mat(k,47) = rxt(k,104)
         mat(k,46) = -( rxt(k,104) + het_rates(k,78) )
         mat(k,25) = rxt(k,100)
         mat(k,15) = -( het_rates(k,79) )
         mat(k,16) = -( het_rates(k,80) )
         mat(k,17) = -( het_rates(k,81) )
         mat(k,18) = -( het_rates(k,82) )
         mat(k,19) = -( het_rates(k,83) )
         mat(k,20) = -( het_rates(k,84) )
         mat(k,454) = -( het_rates(k,85) )
         mat(k,40) = rxt(k,27)
         mat(k,583) = rxt(k,28)
         mat(k,140) = rxt(k,30)
         mat(k,239) = rxt(k,52)
         mat(k,208) = rxt(k,58)
         mat(k,633) = rxt(k,476)*y(k,8) + rxt(k,502)*y(k,9) + 3.000_r8*rxt(k,503)*y(k,23) &
                      + 2.000_r8*rxt(k,504)*y(k,40) + 2.000_r8*rxt(k,525)*y(k,15) &
                      + rxt(k,526)*y(k,17)
         mat(k,1451) = 2.000_r8*rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) &
                      + 3.000_r8*rxt(k,520)*y(k,23)
         mat(k,788) = 2.000_r8*rxt(k,514)*y(k,15) + rxt(k,516)*y(k,17) &
                      + 3.000_r8*rxt(k,521)*y(k,23)
         mat(k,1476) = -( rxt(k,451)*y(k,22) + rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) &
                      + rxt(k,518)*y(k,19) + rxt(k,520)*y(k,23) + het_rates(k,86) )
         mat(k,41) = rxt(k,27)
         mat(k,31) = 2.000_r8*rxt(k,44)
         mat(k,22) = 2.000_r8*rxt(k,45)
         mat(k,1754) = rxt(k,46)
         mat(k,284) = rxt(k,47)
         mat(k,35) = rxt(k,50)
         mat(k,1167) = rxt(k,56)
         mat(k,235) = rxt(k,59)
         mat(k,650) = 4.000_r8*rxt(k,475)*y(k,7) + rxt(k,476)*y(k,8) &
                      + 2.000_r8*rxt(k,477)*y(k,10) + 2.000_r8*rxt(k,478)*y(k,11) &
                      + 2.000_r8*rxt(k,479)*y(k,12) + rxt(k,480)*y(k,13) &
                      + 2.000_r8*rxt(k,481)*y(k,14) + rxt(k,527)*y(k,44) &
                      + rxt(k,528)*y(k,45) + rxt(k,529)*y(k,46)
         mat(k,807) = 3.000_r8*rxt(k,517)*y(k,18) + rxt(k,519)*y(k,19) &
                      + rxt(k,522)*y(k,44) + rxt(k,523)*y(k,45) + rxt(k,524)*y(k,46)
         mat(k,690) = -( het_rates(k,87) )
         mat(k,410) = rxt(k,385)
         mat(k,394) = rxt(k,386)
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
         mat(k,409) = -( rxt(k,385) + het_rates(k,88) )
         mat(k,393) = -( rxt(k,386) + het_rates(k,89) )
         mat(k,287) = -( het_rates(k,90) )
         mat(k,822) = -( rxt(k,63) + het_rates(k,91) )
         mat(k,542) = rxt(k,64) + rxt(k,283)
         mat(k,440) = rxt(k,330)*y(k,30)
         mat(k,1653) = rxt(k,352)*y(k,30)
         mat(k,362) = -( rxt(k,284) + het_rates(k,92) )
         mat(k,539) = -( rxt(k,64) + rxt(k,283) + het_rates(k,93) )
         mat(k,363) = rxt(k,284)
         mat(k,660) = -( rxt(k,65) + het_rates(k,94) )
         mat(k,1311) = rxt(k,316)*y(k,30)
         mat(k,97) = rxt(k,335)*y(k,30)
         mat(k,1561) = -( het_rates(k,95) )
         mat(k,1519) = rxt(k,16)
         mat(k,321) = rxt(k,66)
         mat(k,868) = rxt(k,75)
         mat(k,729) = rxt(k,76)
         mat(k,1436) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82)
         mat(k,1259) = rxt(k,83) + rxt(k,84) + rxt(k,85) + rxt(k,86) + rxt(k,87) &
                      + rxt(k,90) + rxt(k,91) + rxt(k,92)
         mat(k,1333) = rxt(k,93) + rxt(k,137)
         mat(k,448) = rxt(k,96)
         mat(k,1635) = rxt(k,98)
         mat(k,1670) = rxt(k,99) + rxt(k,358)
         mat(k,83) = -( het_rates(k,96) )
         mat(k,335) = -( het_rates(k,97) )
         mat(k,1612) = rxt(k,343)*y(k,30)
         mat(k,747) = -( rxt(k,559) + het_rates(k,98) )
         mat(k,102) = rxt(k,11) + rxt(k,448)
         mat(k,1459) = rxt(k,515)*y(k,17) + rxt(k,518)*y(k,19)
         mat(k,793) = rxt(k,516)*y(k,17) + rxt(k,519)*y(k,19)
         mat(k,638) = rxt(k,539)*y(k,22)
         mat(k,162) = -( rxt(k,373) + het_rates(k,99) )
         mat(k,1011) = rxt(k,375)
         mat(k,1033) = -( rxt(k,375) + het_rates(k,100) )
         mat(k,1700) = rxt(k,377)
         mat(k,50) = -( het_rates(k,101) )
         mat(k,1716) = -( rxt(k,377) + het_rates(k,102) )
         mat(k,1607) = rxt(k,379)
         mat(k,54) = -( het_rates(k,103) )
         mat(k,1604) = -( rxt(k,379) + het_rates(k,104) )
         mat(k,181) = -( het_rates(k,105) )
         mat(k,163) = rxt(k,373)
         mat(k,215) = -( het_rates(k,106) )
         mat(k,153) = -( het_rates(k,107) )
         mat(k,129) = rxt(k,600)
         mat(k,123) = rxt(k,605)
         mat(k,852) = -( rxt(k,75) + het_rates(k,108) )
         mat(k,425) = rxt(k,273)
         mat(k,424) = -( rxt(k,273) + het_rates(k,109) )
         mat(k,713) = -( rxt(k,76) + het_rates(k,110) )
         mat(k,504) = rxt(k,275)
         mat(k,561) = rxt(k,282)
         mat(k,519) = -( rxt(k,274) + het_rates(k,111) )
         mat(k,560) = -( rxt(k,282) + het_rates(k,112) )
         mat(k,520) = rxt(k,274)
         mat(k,465) = -( het_rates(k,113) )
         mat(k,503) = -( rxt(k,275) + het_rates(k,114) )
         mat(k,900) = -( rxt(k,370) + rxt(k,368)*y(k,30) + het_rates(k,115) )
         mat(k,1504) = rxt(k,16)
         mat(k,116) = rxt(k,369)
         mat(k,110) = rxt(k,371)
         mat(k,1362) = rxt(k,593)
         mat(k,269) = rxt(k,596)
         mat(k,943) = -( het_rates(k,116) )
         mat(k,89) = -( het_rates(k,117) )
         mat(k,115) = -( rxt(k,369) + het_rates(k,118) )
         mat(k,109) = rxt(k,312)*y(k,30)
         mat(k,877) = rxt(k,368)*y(k,30)
         mat(k,1075) = -( het_rates(k,119) )
         mat(k,108) = -( rxt(k,371) + rxt(k,312)*y(k,30) + het_rates(k,120) )
         mat(k,876) = rxt(k,370)
         mat(k,223) = -( het_rates(k,121) )
         mat(k,309) = rxt(k,66)
         mat(k,124) = rxt(k,606)
         mat(k,637) = -( rxt(k,391) + rxt(k,442)*y(k,55) + rxt(k,443)*y(k,55) &
                      + rxt(k,475)*y(k,7) + rxt(k,476)*y(k,8) + rxt(k,477)*y(k,10) &
                      + rxt(k,478)*y(k,11) + rxt(k,479)*y(k,12) + rxt(k,480)*y(k,13) &
                      + rxt(k,481)*y(k,14) + rxt(k,502)*y(k,9) + rxt(k,503)*y(k,23) &
                      + rxt(k,504)*y(k,40) + rxt(k,525)*y(k,15) + rxt(k,526)*y(k,17) &
                      + rxt(k,527)*y(k,44) + rxt(k,528)*y(k,45) + rxt(k,529)*y(k,46) &
                      + rxt(k,538)*y(k,22) + rxt(k,539)*y(k,22) + rxt(k,540)*y(k,22) &
                 + het_rates(k,122) )
         mat(k,1827) = rxt(k,1)
         mat(k,1236) = rxt(k,6)
         mat(k,1187) = rxt(k,7)
         mat(k,246) = -( rxt(k,387) + rxt(k,395) + het_rates(k,123) )
         mat(k,1179) = rxt(k,7)
         mat(k,43) = rxt(k,399) + rxt(k,398)*y(k,30)
         mat(k,42) = -( rxt(k,388) + rxt(k,399) + rxt(k,398)*y(k,30) + het_rates(k,124) )
         mat(k,1328) = -( rxt(k,93) + rxt(k,137) + rxt(k,316)*y(k,30) + het_rates(k,125) &
       )
         mat(k,675) = rxt(k,65)
         mat(k,99) = rxt(k,97)
         mat(k,1373) = -( rxt(k,593) + het_rates(k,126) )
         mat(k,1255) = rxt(k,84) + rxt(k,86)
         mat(k,193) = rxt(k,94)
         mat(k,271) = rxt(k,595)*y(k,30)
         mat(k,189) = -( rxt(k,94) + het_rates(k,127) )
         mat(k,438) = -( rxt(k,95) + rxt(k,96) + rxt(k,330)*y(k,30) + het_rates(k,128) )
         mat(k,95) = -( rxt(k,97) + rxt(k,335)*y(k,30) + het_rates(k,129) )
         mat(k,350) = -( het_rates(k,130) )
         mat(k,794) = -( rxt(k,514)*y(k,15) + rxt(k,516)*y(k,17) + rxt(k,517)*y(k,18) &
                      + rxt(k,519)*y(k,19) + rxt(k,521)*y(k,23) + rxt(k,522)*y(k,44) &
                      + rxt(k,523)*y(k,45) + rxt(k,524)*y(k,46) + rxt(k,536)*y(k,22) &
                 + het_rates(k,131) )
         mat(k,1831) = rxt(k,3)
         mat(k,148) = 2.000_r8*rxt(k,4)
         mat(k,1280) = rxt(k,9)
         mat(k,103) = rxt(k,10)
         mat(k,176) = rxt(k,12)
         mat(k,80) = rxt(k,23)
         mat(k,210) = rxt(k,58)
         mat(k,232) = rxt(k,59)
         mat(k,1617) = rxt(k,98)
         mat(k,985) = .500_r8*rxt(k,561)
         mat(k,639) = rxt(k,538)*y(k,22)
         mat(k,1637) = -( rxt(k,98) + rxt(k,343)*y(k,30) + het_rates(k,132) )
         mat(k,1673) = -( rxt(k,99) + rxt(k,358) + rxt(k,352)*y(k,30) + het_rates(k,133) &
       )
         mat(k,842) = rxt(k,63)
         mat(k,450) = rxt(k,95)
         mat(k,264) = -( rxt(k,596) + rxt(k,595)*y(k,30) + het_rates(k,134) )
         mat(k,1402) = rxt(k,77) + rxt(k,81)
         mat(k,1228) = rxt(k,85) + rxt(k,87)
         mat(k,125) = rxt(k,580)
         mat(k,130) = rxt(k,581)
         mat(k,128) = -( rxt(k,581) + rxt(k,600) + het_rates(k,135) )
         mat(k,1389) = rxt(k,78) + rxt(k,82)
         mat(k,1221) = rxt(k,83) + rxt(k,92)
         mat(k,122) = rxt(k,582)
         mat(k,121) = -( rxt(k,580) + rxt(k,582) + rxt(k,605) + rxt(k,606) &
                 + het_rates(k,136) )
         mat(k,1388) = rxt(k,79) + rxt(k,80)
         mat(k,1220) = rxt(k,90) + rxt(k,91)
         mat(k,1856) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,137) )
         mat(k,556) = rxt(k,64) + rxt(k,283)
         mat(k,196) = rxt(k,94)
         mat(k,26) = rxt(k,100)
         mat(k,437) = rxt(k,273)
         mat(k,536) = rxt(k,274)
         mat(k,576) = rxt(k,282)
         mat(k,375) = rxt(k,284)
         mat(k,167) = rxt(k,373)
         mat(k,1052) = rxt(k,375)
         mat(k,1719) = rxt(k,377)
         mat(k,1610) = rxt(k,379)
         mat(k,423) = rxt(k,385)
         mat(k,814) = rxt(k,514)*y(k,15) + rxt(k,516)*y(k,17) + rxt(k,517)*y(k,18) &
                      + rxt(k,519)*y(k,19) + rxt(k,524)*y(k,46) + rxt(k,536)*y(k,22)
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
