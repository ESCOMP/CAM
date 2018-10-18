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
         mat(k,565) = -( rxt(k,28) + het_rates(k,4) )
         mat(k,140) = rxt(k,29)
         mat(k,137) = -( rxt(k,29) + rxt(k,30) + rxt(k,564) + rxt(k,567) + rxt(k,572) &
                 + het_rates(k,5) )
         mat(k,587) = -( rxt(k,21) + rxt(k,22) + het_rates(k,16) )
         mat(k,79) = rxt(k,23)
         mat(k,617) = rxt(k,538)*y(k,22) + rxt(k,539)*y(k,22)
         mat(k,305) = -( het_rates(k,20) )
         mat(k,873) = rxt(k,450)*y(k,22)
         mat(k,210) = rxt(k,506)*y(k,22)
         mat(k,767) = rxt(k,535)*y(k,22)
         mat(k,612) = rxt(k,537)*y(k,22)
         mat(k,77) = -( rxt(k,23) + het_rates(k,21) )
         mat(k,30) = -( rxt(k,44) + het_rates(k,24) )
         mat(k,21) = -( rxt(k,45) + rxt(k,484) + het_rates(k,25) )
         mat(k,1105) = -( rxt(k,46) + het_rates(k,26) )
         mat(k,262) = rxt(k,48)
         mat(k,67) = rxt(k,60)
         mat(k,23) = 2.000_r8*rxt(k,484)
         mat(k,259) = -( rxt(k,47) + rxt(k,48) + rxt(k,566) + rxt(k,571) + rxt(k,577) &
                 + het_rates(k,27) )
         mat(k,133) = -( het_rates(k,29) )
         mat(k,582) = rxt(k,21) + rxt(k,22)
         mat(k,70) = rxt(k,100)
         mat(k,868) = rxt(k,517)*y(k,19)
         mat(k,198) = rxt(k,591)*y(k,30)
         mat(k,27) = -( rxt(k,49) + het_rates(k,31) )
         mat(k,606) = rxt(k,475)*y(k,8) + rxt(k,477)*y(k,11) + 2.000_r8*rxt(k,478)*y(k,12) &
                      + 2.000_r8*rxt(k,479)*y(k,13) + rxt(k,480)*y(k,14) &
                      + rxt(k,501)*y(k,9) + 2.000_r8*rxt(k,503)*y(k,40) &
                      + rxt(k,527)*y(k,45) + rxt(k,528)*y(k,46)
         mat(k,749) = rxt(k,522)*y(k,45) + rxt(k,523)*y(k,46)
         mat(k,32) = -( rxt(k,50) + het_rates(k,32) )
         mat(k,607) = rxt(k,476)*y(k,10) + rxt(k,477)*y(k,11) + rxt(k,526)*y(k,44)
         mat(k,750) = rxt(k,521)*y(k,44)
         mat(k,58) = -( het_rates(k,33) )
         mat(k,3) = -( het_rates(k,34) )
         mat(k,4) = -( het_rates(k,35) )
         mat(k,5) = -( het_rates(k,36) )
         mat(k,209) = -( rxt(k,506)*y(k,22) + het_rates(k,37) )
         mat(k,28) = 2.000_r8*rxt(k,49)
         mat(k,33) = rxt(k,50)
         mat(k,37) = rxt(k,57)
         mat(k,609) = rxt(k,479)*y(k,13) + rxt(k,501)*y(k,9)
         mat(k,1223) = -( het_rates(k,38) )
         mat(k,1814) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,596) = 2.000_r8*rxt(k,21)
         mat(k,81) = rxt(k,23)
         mat(k,231) = rxt(k,52)
         mat(k,1314) = rxt(k,56)
         mat(k,38) = rxt(k,57)
         mat(k,626) = rxt(k,538)*y(k,22)
         mat(k,465) = -( het_rates(k,39) )
         mat(k,1793) = rxt(k,1)
         mat(k,585) = rxt(k,22)
         mat(k,615) = rxt(k,539)*y(k,22)
         mat(k,145) = -( rxt(k,4) + het_rates(k,41) )
         mat(k,717) = .500_r8*rxt(k,558)
         mat(k,24) = -( rxt(k,99) + het_rates(k,42) )
         mat(k,226) = -( rxt(k,52) + het_rates(k,43) )
         mat(k,1316) = -( rxt(k,56) + het_rates(k,47) )
         mat(k,383) = rxt(k,385)
         mat(k,896) = rxt(k,450)*y(k,22) + rxt(k,512)*y(k,15) + rxt(k,514)*y(k,17) &
                      + 2.000_r8*rxt(k,517)*y(k,19) + rxt(k,519)*y(k,23)
         mat(k,36) = -( rxt(k,57) + het_rates(k,48) )
         mat(k,208) = rxt(k,506)*y(k,22)
         mat(k,1356) = -( rxt(k,9) + het_rates(k,49) )
         mat(k,492) = rxt(k,274)
         mat(k,285) = 2.000_r8*rxt(k,559) + 2.000_r8*rxt(k,562) + 2.000_r8*rxt(k,565) &
                      + 2.000_r8*rxt(k,576)
         mat(k,1481) = .500_r8*rxt(k,560)
         mat(k,1191) = rxt(k,561)
         mat(k,142) = rxt(k,564) + rxt(k,567) + rxt(k,572)
         mat(k,265) = rxt(k,566) + rxt(k,571) + rxt(k,577)
         mat(k,101) = -( rxt(k,10) + rxt(k,11) + rxt(k,447) + het_rates(k,50) )
         mat(k,218) = -( rxt(k,58) + het_rates(k,51) )
         mat(k,138) = rxt(k,564) + rxt(k,567) + rxt(k,572)
         mat(k,235) = -( rxt(k,59) + het_rates(k,52) )
         mat(k,258) = rxt(k,566) + rxt(k,571) + rxt(k,577)
         mat(k,176) = -( rxt(k,12) + het_rates(k,53) )
         mat(k,293) = -( rxt(k,65) + het_rates(k,54) )
         mat(k,1537) = rxt(k,17)
         mat(k,200) = rxt(k,592)
         mat(k,279) = -( rxt(k,14) + rxt(k,15) + rxt(k,448) + rxt(k,559) + rxt(k,562) &
                      + rxt(k,565) + rxt(k,576) + het_rates(k,56) )
         mat(k,6) = -( het_rates(k,57) )
         mat(k,7) = -( het_rates(k,58) )
         mat(k,8) = -( het_rates(k,59) )
         mat(k,1569) = -( rxt(k,16) + rxt(k,17) + het_rates(k,60) )
         mat(k,179) = rxt(k,12)
         mat(k,288) = rxt(k,15)
         mat(k,1486) = rxt(k,18) + .500_r8*rxt(k,560)
         mat(k,1196) = rxt(k,20)
         mat(k,1440) = rxt(k,589)
         mat(k,633) = 2.000_r8*rxt(k,441)*y(k,55)
         mat(k,1484) = -( rxt(k,18) + rxt(k,560) + het_rates(k,61) )
         mat(k,1359) = rxt(k,9)
         mat(k,106) = rxt(k,11) + rxt(k,447)
         mat(k,286) = rxt(k,14) + rxt(k,448)
         mat(k,1194) = rxt(k,19)
         mat(k,143) = rxt(k,29)
         mat(k,266) = rxt(k,48)
         mat(k,1653) = rxt(k,74)
         mat(k,1187) = -( rxt(k,19) + rxt(k,20) + rxt(k,561) + het_rates(k,62) )
         mat(k,104) = rxt(k,10)
         mat(k,284) = rxt(k,14) + rxt(k,15) + rxt(k,448)
         mat(k,141) = rxt(k,30)
         mat(k,263) = rxt(k,47)
         mat(k,702) = rxt(k,75)
         mat(k,9) = -( het_rates(k,63) )
         mat(k,10) = -( het_rates(k,64) )
         mat(k,11) = -( het_rates(k,65) )
         mat(k,12) = -( het_rates(k,66) )
         mat(k,1627) = -( rxt(k,76) + rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) &
                      + rxt(k,81) + het_rates(k,67) )
         mat(k,1823) = rxt(k,2)
         mat(k,1279) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,82) + rxt(k,84) &
                      + 2.000_r8*rxt(k,85) + rxt(k,86) + rxt(k,88) + 2.000_r8*rxt(k,89) &
                      + rxt(k,90) + rxt(k,91)
         mat(k,944) = rxt(k,8)
         mat(k,289) = rxt(k,15)
         mat(k,1570) = rxt(k,17)
         mat(k,1487) = rxt(k,18)
         mat(k,1197) = rxt(k,19)
         mat(k,580) = rxt(k,28)
         mat(k,1117) = rxt(k,46)
         mat(k,68) = rxt(k,60)
         mat(k,1723) = rxt(k,98) + rxt(k,357)
         mat(k,372) = rxt(k,101)
         mat(k,255) = rxt(k,102)
         mat(k,48) = rxt(k,103)
         mat(k,634) = rxt(k,390)
         mat(k,1271) = -( rxt(k,5) + rxt(k,6) + rxt(k,82) + rxt(k,83) + rxt(k,84) &
                      + rxt(k,85) + rxt(k,86) + rxt(k,87) + rxt(k,88) + rxt(k,89) &
                      + rxt(k,90) + rxt(k,91) + het_rates(k,68) )
         mat(k,936) = rxt(k,8)
         mat(k,1189) = rxt(k,20)
         mat(k,1389) = rxt(k,92) + rxt(k,136)
         mat(k,425) = rxt(k,95) + rxt(k,329)*y(k,30)
         mat(k,98) = rxt(k,96) + rxt(k,334)*y(k,30)
         mat(k,246) = rxt(k,386) + rxt(k,394)
         mat(k,44) = rxt(k,387)
         mat(k,627) = rxt(k,442)*y(k,55)
         mat(k,928) = -( rxt(k,7) + rxt(k,8) + het_rates(k,69) )
         mat(k,423) = rxt(k,94)
         mat(k,64) = -( rxt(k,60) + het_rates(k,70) )
         mat(k,69) = -( rxt(k,100) + het_rates(k,71) )
         mat(k,13) = -( het_rates(k,72) )
         mat(k,14) = -( het_rates(k,73) )
         mat(k,169) = -( het_rates(k,74) )
         mat(k,71) = rxt(k,100)
         mat(k,358) = rxt(k,101)
         mat(k,360) = -( rxt(k,101) + het_rates(k,75) )
         mat(k,252) = rxt(k,102)
         mat(k,251) = -( rxt(k,102) + het_rates(k,76) )
         mat(k,47) = rxt(k,103)
         mat(k,46) = -( rxt(k,103) + het_rates(k,77) )
         mat(k,25) = rxt(k,99)
         mat(k,15) = -( het_rates(k,78) )
         mat(k,16) = -( het_rates(k,79) )
         mat(k,17) = -( het_rates(k,80) )
         mat(k,18) = -( het_rates(k,81) )
         mat(k,19) = -( het_rates(k,82) )
         mat(k,20) = -( het_rates(k,83) )
         mat(k,451) = -( het_rates(k,84) )
         mat(k,40) = rxt(k,27)
         mat(k,564) = rxt(k,28)
         mat(k,139) = rxt(k,30)
         mat(k,227) = rxt(k,52)
         mat(k,219) = rxt(k,58)
         mat(k,614) = rxt(k,475)*y(k,8) + rxt(k,501)*y(k,9) + 3.000_r8*rxt(k,502)*y(k,23) &
                      + 2.000_r8*rxt(k,503)*y(k,40) + 2.000_r8*rxt(k,524)*y(k,15) &
                      + rxt(k,525)*y(k,17)
         mat(k,874) = 2.000_r8*rxt(k,512)*y(k,15) + rxt(k,514)*y(k,17) &
                      + 3.000_r8*rxt(k,519)*y(k,23)
         mat(k,769) = 2.000_r8*rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) &
                      + 3.000_r8*rxt(k,520)*y(k,23)
         mat(k,886) = -( rxt(k,450)*y(k,22) + rxt(k,512)*y(k,15) + rxt(k,514)*y(k,17) &
                      + rxt(k,517)*y(k,19) + rxt(k,519)*y(k,23) + het_rates(k,85) )
         mat(k,41) = rxt(k,27)
         mat(k,31) = 2.000_r8*rxt(k,44)
         mat(k,22) = 2.000_r8*rxt(k,45)
         mat(k,1100) = rxt(k,46)
         mat(k,261) = rxt(k,47)
         mat(k,35) = rxt(k,50)
         mat(k,1306) = rxt(k,56)
         mat(k,237) = rxt(k,59)
         mat(k,622) = 4.000_r8*rxt(k,474)*y(k,7) + rxt(k,475)*y(k,8) &
                      + 2.000_r8*rxt(k,476)*y(k,10) + 2.000_r8*rxt(k,477)*y(k,11) &
                      + 2.000_r8*rxt(k,478)*y(k,12) + rxt(k,479)*y(k,13) &
                      + 2.000_r8*rxt(k,480)*y(k,14) + rxt(k,526)*y(k,44) &
                      + rxt(k,527)*y(k,45) + rxt(k,528)*y(k,46)
         mat(k,777) = 3.000_r8*rxt(k,516)*y(k,18) + rxt(k,518)*y(k,19) &
                      + rxt(k,521)*y(k,44) + rxt(k,522)*y(k,45) + rxt(k,523)*y(k,46)
         mat(k,671) = -( het_rates(k,86) )
         mat(k,391) = rxt(k,384)
         mat(k,375) = rxt(k,385)
         mat(k,390) = -( rxt(k,384) + het_rates(k,87) )
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
         mat(k,374) = -( rxt(k,385) + het_rates(k,88) )
         mat(k,269) = -( het_rates(k,89) )
         mat(k,803) = -( rxt(k,62) + het_rates(k,90) )
         mat(k,523) = rxt(k,63) + rxt(k,282)
         mat(k,421) = rxt(k,329)*y(k,30)
         mat(k,1704) = rxt(k,351)*y(k,30)
         mat(k,343) = -( rxt(k,283) + het_rates(k,91) )
         mat(k,520) = -( rxt(k,63) + rxt(k,282) + het_rates(k,92) )
         mat(k,344) = rxt(k,283)
         mat(k,641) = -( rxt(k,64) + het_rates(k,93) )
         mat(k,1373) = rxt(k,315)*y(k,30)
         mat(k,97) = rxt(k,334)*y(k,30)
         mat(k,844) = -( het_rates(k,94) )
         mat(k,1552) = rxt(k,16)
         mat(k,296) = rxt(k,65)
         mat(k,1638) = rxt(k,74)
         mat(k,695) = rxt(k,75)
         mat(k,1609) = rxt(k,76) + rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) &
                      + rxt(k,81)
         mat(k,1261) = rxt(k,82) + rxt(k,83) + rxt(k,84) + rxt(k,86) + rxt(k,87) &
                      + rxt(k,88) + rxt(k,90) + rxt(k,91)
         mat(k,1379) = rxt(k,92) + rxt(k,136)
         mat(k,422) = rxt(k,94)
         mat(k,1670) = rxt(k,97)
         mat(k,1705) = rxt(k,98) + rxt(k,357)
         mat(k,83) = -( het_rates(k,95) )
         mat(k,316) = -( het_rates(k,96) )
         mat(k,1663) = rxt(k,342)*y(k,30)
         mat(k,728) = -( rxt(k,558) + het_rates(k,97) )
         mat(k,102) = rxt(k,11) + rxt(k,447)
         mat(k,882) = rxt(k,514)*y(k,17) + rxt(k,517)*y(k,19)
         mat(k,774) = rxt(k,515)*y(k,17) + rxt(k,518)*y(k,19)
         mat(k,619) = rxt(k,538)*y(k,22)
         mat(k,152) = -( rxt(k,372) + het_rates(k,98) )
         mat(k,1729) = rxt(k,374)
         mat(k,1769) = -( rxt(k,374) + het_rates(k,99) )
         mat(k,1533) = rxt(k,376)
         mat(k,50) = -( het_rates(k,100) )
         mat(k,1527) = -( rxt(k,376) + het_rates(k,101) )
         mat(k,1071) = rxt(k,378)
         mat(k,54) = -( het_rates(k,102) )
         mat(k,1060) = -( rxt(k,378) + het_rates(k,103) )
         mat(k,182) = -( het_rates(k,104) )
         mat(k,153) = rxt(k,372)
         mat(k,128) = -( het_rates(k,105) )
         mat(k,159) = -( het_rates(k,106) )
         mat(k,1657) = -( rxt(k,74) + het_rates(k,107) )
         mat(k,416) = rxt(k,272)
         mat(k,405) = -( rxt(k,272) + het_rates(k,108) )
         mat(k,694) = -( rxt(k,75) + het_rates(k,109) )
         mat(k,485) = rxt(k,274)
         mat(k,542) = rxt(k,281)
         mat(k,500) = -( rxt(k,273) + het_rates(k,110) )
         mat(k,541) = -( rxt(k,281) + het_rates(k,111) )
         mat(k,501) = rxt(k,273)
         mat(k,433) = -( het_rates(k,112) )
         mat(k,484) = -( rxt(k,274) + het_rates(k,113) )
         mat(k,974) = -( rxt(k,369) + rxt(k,367)*y(k,30) + het_rates(k,114) )
         mat(k,1555) = rxt(k,16)
         mat(k,117) = rxt(k,368)
         mat(k,111) = rxt(k,370)
         mat(k,1426) = rxt(k,589)
         mat(k,203) = rxt(k,592)
         mat(k,1017) = -( het_rates(k,115) )
         mat(k,89) = -( het_rates(k,116) )
         mat(k,115) = -( rxt(k,368) + het_rates(k,117) )
         mat(k,109) = rxt(k,311)*y(k,30)
         mat(k,951) = rxt(k,367)*y(k,30)
         mat(k,1147) = -( het_rates(k,118) )
         mat(k,108) = -( rxt(k,370) + rxt(k,311)*y(k,30) + het_rates(k,119) )
         mat(k,950) = rxt(k,369)
         mat(k,121) = -( het_rates(k,120) )
         mat(k,291) = rxt(k,65)
         mat(k,618) = -( rxt(k,390) + rxt(k,441)*y(k,55) + rxt(k,442)*y(k,55) &
                      + rxt(k,474)*y(k,7) + rxt(k,475)*y(k,8) + rxt(k,476)*y(k,10) &
                      + rxt(k,477)*y(k,11) + rxt(k,478)*y(k,12) + rxt(k,479)*y(k,13) &
                      + rxt(k,480)*y(k,14) + rxt(k,501)*y(k,9) + rxt(k,502)*y(k,23) &
                      + rxt(k,503)*y(k,40) + rxt(k,524)*y(k,15) + rxt(k,525)*y(k,17) &
                      + rxt(k,526)*y(k,44) + rxt(k,527)*y(k,45) + rxt(k,528)*y(k,46) &
                      + rxt(k,537)*y(k,22) + rxt(k,538)*y(k,22) + rxt(k,539)*y(k,22) &
                 + het_rates(k,121) )
         mat(k,1799) = rxt(k,1)
         mat(k,1254) = rxt(k,6)
         mat(k,919) = rxt(k,7)
         mat(k,242) = -( rxt(k,386) + rxt(k,394) + het_rates(k,122) )
         mat(k,911) = rxt(k,7)
         mat(k,43) = rxt(k,398) + rxt(k,397)*y(k,30)
         mat(k,42) = -( rxt(k,387) + rxt(k,398) + rxt(k,397)*y(k,30) + het_rates(k,123) )
         mat(k,1392) = -( rxt(k,92) + rxt(k,136) + rxt(k,315)*y(k,30) + het_rates(k,124) &
       )
         mat(k,658) = rxt(k,64)
         mat(k,99) = rxt(k,96)
         mat(k,1437) = -( rxt(k,589) + het_rates(k,125) )
         mat(k,1275) = rxt(k,83) + rxt(k,87)
         mat(k,195) = rxt(k,93)
         mat(k,206) = rxt(k,591)*y(k,30)
         mat(k,190) = -( rxt(k,93) + het_rates(k,126) )
         mat(k,419) = -( rxt(k,94) + rxt(k,95) + rxt(k,329)*y(k,30) + het_rates(k,127) )
         mat(k,95) = -( rxt(k,96) + rxt(k,334)*y(k,30) + het_rates(k,128) )
         mat(k,331) = -( het_rates(k,129) )
         mat(k,775) = -( rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) + rxt(k,516)*y(k,18) &
                      + rxt(k,518)*y(k,19) + rxt(k,520)*y(k,23) + rxt(k,521)*y(k,44) &
                      + rxt(k,522)*y(k,45) + rxt(k,523)*y(k,46) + rxt(k,535)*y(k,22) &
                 + het_rates(k,130) )
         mat(k,1803) = rxt(k,3)
         mat(k,147) = 2.000_r8*rxt(k,4)
         mat(k,1342) = rxt(k,9)
         mat(k,103) = rxt(k,10)
         mat(k,177) = rxt(k,12)
         mat(k,80) = rxt(k,23)
         mat(k,221) = rxt(k,58)
         mat(k,236) = rxt(k,59)
         mat(k,1668) = rxt(k,97)
         mat(k,1467) = .500_r8*rxt(k,560)
         mat(k,620) = rxt(k,537)*y(k,22)
         mat(k,1690) = -( rxt(k,97) + rxt(k,342)*y(k,30) + het_rates(k,131) )
         mat(k,1726) = -( rxt(k,98) + rxt(k,357) + rxt(k,351)*y(k,30) + het_rates(k,132) &
       )
         mat(k,825) = rxt(k,62)
         mat(k,432) = rxt(k,95)
         mat(k,199) = -( rxt(k,592) + rxt(k,591)*y(k,30) + het_rates(k,133) )
         mat(k,1586) = rxt(k,76) + rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) &
                      + rxt(k,81)
         mat(k,1244) = rxt(k,82) + rxt(k,84) + rxt(k,86) + rxt(k,88) + rxt(k,90) &
                      + rxt(k,91)
         mat(k,1828) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,134) )
         mat(k,537) = rxt(k,63) + rxt(k,282)
         mat(k,197) = rxt(k,93)
         mat(k,26) = rxt(k,99)
         mat(k,418) = rxt(k,272)
         mat(k,517) = rxt(k,273)
         mat(k,557) = rxt(k,281)
         mat(k,356) = rxt(k,283)
         mat(k,157) = rxt(k,372)
         mat(k,1770) = rxt(k,374)
         mat(k,1534) = rxt(k,376)
         mat(k,1078) = rxt(k,378)
         mat(k,404) = rxt(k,384)
         mat(k,795) = rxt(k,513)*y(k,15) + rxt(k,515)*y(k,17) + rxt(k,516)*y(k,18) &
                      + rxt(k,518)*y(k,19) + rxt(k,523)*y(k,46) + rxt(k,535)*y(k,22)
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
