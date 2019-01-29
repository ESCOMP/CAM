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
         mat(k,457) = -( rxt(k,28) + het_rates(k,2) )
         mat(k,74) = rxt(k,29)
         mat(k,71) = -( rxt(k,29) + rxt(k,30) + rxt(k,543) + rxt(k,546) + rxt(k,551) &
                 + het_rates(k,3) )
         mat(k,541) = -( rxt(k,21) + rxt(k,22) + het_rates(k,14) )
         mat(k,41) = rxt(k,23)
         mat(k,572) = rxt(k,534)*y(k,20) + rxt(k,535)*y(k,20)
         mat(k,257) = -( het_rates(k,18) )
         mat(k,1658) = rxt(k,446)*y(k,20)
         mat(k,169) = rxt(k,502)*y(k,20)
         mat(k,747) = rxt(k,531)*y(k,20)
         mat(k,567) = rxt(k,533)*y(k,20)
         mat(k,39) = -( rxt(k,23) + het_rates(k,19) )
         mat(k,10) = -( rxt(k,44) + het_rates(k,22) )
         mat(k,1) = -( rxt(k,45) + rxt(k,480) + het_rates(k,23) )
         mat(k,1532) = -( rxt(k,46) + het_rates(k,24) )
         mat(k,208) = rxt(k,48)
         mat(k,9) = rxt(k,60)
         mat(k,2) = 2.000_r8*rxt(k,480)
         mat(k,201) = -( rxt(k,47) + rxt(k,48) + rxt(k,545) + rxt(k,550) + rxt(k,556) &
                 + het_rates(k,25) )
         mat(k,128) = -( het_rates(k,27) )
         mat(k,536) = rxt(k,21) + rxt(k,22)
         mat(k,1376) = rxt(k,26) + rxt(k,62)
         mat(k,1654) = rxt(k,513)*y(k,17)
         mat(k,4) = -( rxt(k,49) + het_rates(k,29) )
         mat(k,561) = rxt(k,471)*y(k,6) + rxt(k,473)*y(k,9) + 2.000_r8*rxt(k,474)*y(k,10) &
                      + 2.000_r8*rxt(k,475)*y(k,11) + rxt(k,476)*y(k,12) &
                      + rxt(k,497)*y(k,7) + 2.000_r8*rxt(k,499)*y(k,34) &
                      + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39)
         mat(k,735) = rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39)
         mat(k,12) = -( rxt(k,50) + het_rates(k,30) )
         mat(k,562) = rxt(k,472)*y(k,8) + rxt(k,473)*y(k,9) + rxt(k,522)*y(k,37)
         mat(k,736) = rxt(k,517)*y(k,37)
         mat(k,1403) = -( rxt(k,26) + rxt(k,62) + het_rates(k,28) )
         mat(k,799) = rxt(k,63)
         mat(k,637) = rxt(k,65)
         mat(k,120) = rxt(k,364)
         mat(k,168) = -( rxt(k,502)*y(k,20) + het_rates(k,31) )
         mat(k,5) = 2.000_r8*rxt(k,49)
         mat(k,13) = rxt(k,50)
         mat(k,20) = rxt(k,57)
         mat(k,564) = rxt(k,475)*y(k,11) + rxt(k,497)*y(k,7)
         mat(k,1728) = -( het_rates(k,32) )
         mat(k,1823) = 2.000_r8*rxt(k,2) + rxt(k,3)
         mat(k,558) = 2.000_r8*rxt(k,21)
         mat(k,43) = rxt(k,23)
         mat(k,197) = rxt(k,52)
         mat(k,1494) = rxt(k,56)
         mat(k,21) = rxt(k,57)
         mat(k,592) = rxt(k,534)*y(k,20)
         mat(k,419) = -( het_rates(k,33) )
         mat(k,1790) = rxt(k,1)
         mat(k,539) = rxt(k,22)
         mat(k,570) = rxt(k,535)*y(k,20)
         mat(k,79) = -( rxt(k,4) + het_rates(k,35) )
         mat(k,675) = .500_r8*rxt(k,537)
         mat(k,191) = -( rxt(k,52) + het_rates(k,36) )
         mat(k,1488) = -( rxt(k,56) + het_rates(k,40) )
         mat(k,349) = rxt(k,381)
         mat(k,1686) = rxt(k,446)*y(k,20) + rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + 2.000_r8*rxt(k,513)*y(k,17) + rxt(k,515)*y(k,21)
         mat(k,19) = -( rxt(k,57) + het_rates(k,41) )
         mat(k,167) = rxt(k,502)*y(k,20)
         mat(k,1107) = -( rxt(k,9) + het_rates(k,42) )
         mat(k,444) = rxt(k,270)
         mat(k,234) = 2.000_r8*rxt(k,538) + 2.000_r8*rxt(k,541) + 2.000_r8*rxt(k,544) &
                      + 2.000_r8*rxt(k,555)
         mat(k,1152) = .500_r8*rxt(k,539)
         mat(k,1230) = rxt(k,540)
         mat(k,76) = rxt(k,543) + rxt(k,546) + rxt(k,551)
         mat(k,204) = rxt(k,545) + rxt(k,550) + rxt(k,556)
         mat(k,686) = -( rxt(k,537) + het_rates(k,43) )
         mat(k,52) = rxt(k,11) + rxt(k,443)
         mat(k,1667) = rxt(k,510)*y(k,15) + rxt(k,513)*y(k,17)
         mat(k,754) = rxt(k,511)*y(k,15) + rxt(k,514)*y(k,17)
         mat(k,574) = rxt(k,534)*y(k,20)
         mat(k,51) = -( rxt(k,10) + rxt(k,11) + rxt(k,443) + het_rates(k,44) )
         mat(k,151) = -( rxt(k,58) + het_rates(k,45) )
         mat(k,72) = rxt(k,543) + rxt(k,546) + rxt(k,551)
         mat(k,184) = -( rxt(k,59) + het_rates(k,46) )
         mat(k,200) = rxt(k,545) + rxt(k,550) + rxt(k,556)
         mat(k,123) = -( rxt(k,12) + het_rates(k,47) )
         mat(k,270) = -( rxt(k,66) + het_rates(k,48) )
         mat(k,1734) = rxt(k,17)
         mat(k,245) = rxt(k,574)
         mat(k,229) = -( rxt(k,14) + rxt(k,15) + rxt(k,444) + rxt(k,538) + rxt(k,541) &
                      + rxt(k,544) + rxt(k,555) + het_rates(k,50) )
         mat(k,1771) = -( rxt(k,16) + rxt(k,17) + het_rates(k,51) )
         mat(k,126) = rxt(k,12)
         mat(k,239) = rxt(k,15)
         mat(k,1168) = rxt(k,18) + .500_r8*rxt(k,539)
         mat(k,1246) = rxt(k,20)
         mat(k,1370) = rxt(k,571)
         mat(k,64) = rxt(k,584)
         mat(k,593) = 2.000_r8*rxt(k,437)*y(k,49)
         mat(k,1153) = -( rxt(k,18) + rxt(k,539) + het_rates(k,52) )
         mat(k,1108) = rxt(k,9)
         mat(k,54) = rxt(k,11) + rxt(k,443)
         mat(k,235) = rxt(k,14) + rxt(k,444)
         mat(k,1231) = rxt(k,19)
         mat(k,77) = rxt(k,29)
         mat(k,205) = rxt(k,48)
         mat(k,721) = rxt(k,75)
         mat(k,1233) = -( rxt(k,19) + rxt(k,20) + rxt(k,540) + het_rates(k,53) )
         mat(k,55) = rxt(k,10)
         mat(k,237) = rxt(k,14) + rxt(k,15) + rxt(k,444)
         mat(k,78) = rxt(k,30)
         mat(k,206) = rxt(k,47)
         mat(k,663) = rxt(k,76)
         mat(k,979) = -( rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82) + het_rates(k,54) )
         mat(k,1805) = rxt(k,2)
         mat(k,1270) = 2.000_r8*rxt(k,5) + rxt(k,6) + rxt(k,83) + 2.000_r8*rxt(k,84) &
                      + rxt(k,85) + rxt(k,87) + 2.000_r8*rxt(k,89) + rxt(k,90) + rxt(k,91) &
                      + rxt(k,92)
         mat(k,1431) = rxt(k,8)
         mat(k,233) = rxt(k,15)
         mat(k,1752) = rxt(k,17)
         mat(k,1149) = rxt(k,18)
         mat(k,1227) = rxt(k,19)
         mat(k,1393) = rxt(k,26) + rxt(k,62)
         mat(k,462) = rxt(k,28)
         mat(k,1519) = rxt(k,46)
         mat(k,8) = rxt(k,60)
         mat(k,1631) = rxt(k,99) + rxt(k,353)
         mat(k,577) = rxt(k,386)
         mat(k,68) = rxt(k,578)
         mat(k,63) = rxt(k,583)
         mat(k,1277) = -( rxt(k,5) + rxt(k,6) + rxt(k,83) + rxt(k,84) + rxt(k,85) &
                      + rxt(k,86) + rxt(k,87) + rxt(k,88) + rxt(k,89) + rxt(k,90) &
                      + rxt(k,91) + rxt(k,92) + het_rates(k,55) )
         mat(k,1438) = rxt(k,8)
         mat(k,1234) = rxt(k,20)
         mat(k,1313) = rxt(k,93) + rxt(k,132)
         mat(k,408) = rxt(k,95)
         mat(k,99) = rxt(k,97)
         mat(k,222) = rxt(k,382) + rxt(k,390)
         mat(k,32) = rxt(k,383)
         mat(k,582) = rxt(k,438)*y(k,49)
         mat(k,1442) = -( rxt(k,7) + rxt(k,8) + het_rates(k,56) )
         mat(k,411) = rxt(k,96)
         mat(k,7) = -( rxt(k,60) + het_rates(k,57) )
         mat(k,328) = -( het_rates(k,59) )
         mat(k,17) = rxt(k,27)
         mat(k,456) = rxt(k,28)
         mat(k,73) = rxt(k,30)
         mat(k,192) = rxt(k,52)
         mat(k,152) = rxt(k,58)
         mat(k,569) = rxt(k,471)*y(k,6) + rxt(k,497)*y(k,7) + 3.000_r8*rxt(k,498)*y(k,21) &
                      + 2.000_r8*rxt(k,499)*y(k,34) + 2.000_r8*rxt(k,520)*y(k,13) &
                      + rxt(k,521)*y(k,15)
         mat(k,1659) = 2.000_r8*rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + 3.000_r8*rxt(k,515)*y(k,21)
         mat(k,749) = 2.000_r8*rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) &
                      + 3.000_r8*rxt(k,516)*y(k,21)
         mat(k,1691) = -( rxt(k,446)*y(k,20) + rxt(k,508)*y(k,13) + rxt(k,510)*y(k,15) &
                      + rxt(k,513)*y(k,17) + rxt(k,515)*y(k,21) + het_rates(k,60) )
         mat(k,18) = rxt(k,27)
         mat(k,11) = 2.000_r8*rxt(k,44)
         mat(k,3) = 2.000_r8*rxt(k,45)
         mat(k,1536) = rxt(k,46)
         mat(k,209) = rxt(k,47)
         mat(k,15) = rxt(k,50)
         mat(k,1493) = rxt(k,56)
         mat(k,189) = rxt(k,59)
         mat(k,591) = 4.000_r8*rxt(k,470)*y(k,5) + rxt(k,471)*y(k,6) &
                      + 2.000_r8*rxt(k,472)*y(k,8) + 2.000_r8*rxt(k,473)*y(k,9) &
                      + 2.000_r8*rxt(k,474)*y(k,10) + rxt(k,475)*y(k,11) &
                      + 2.000_r8*rxt(k,476)*y(k,12) + rxt(k,522)*y(k,37) &
                      + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39)
         mat(k,773) = 3.000_r8*rxt(k,512)*y(k,16) + rxt(k,514)*y(k,17) &
                      + rxt(k,517)*y(k,37) + rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39)
         mat(k,597) = -( het_rates(k,61) )
         mat(k,357) = rxt(k,380)
         mat(k,341) = rxt(k,381)
         mat(k,356) = -( rxt(k,380) + het_rates(k,62) )
         mat(k,340) = -( rxt(k,381) + het_rates(k,63) )
         mat(k,211) = -( het_rates(k,64) )
         mat(k,785) = -( rxt(k,63) + het_rates(k,65) )
         mat(k,501) = rxt(k,64) + rxt(k,278)
         mat(k,311) = -( rxt(k,279) + het_rates(k,66) )
         mat(k,497) = -( rxt(k,64) + rxt(k,278) + het_rates(k,67) )
         mat(k,312) = rxt(k,279)
         mat(k,619) = -( rxt(k,65) + het_rates(k,68) )
         mat(k,1023) = -( het_rates(k,69) )
         mat(k,1753) = rxt(k,16)
         mat(k,276) = rxt(k,66)
         mat(k,718) = rxt(k,75)
         mat(k,658) = rxt(k,76)
         mat(k,980) = rxt(k,77) + rxt(k,78) + rxt(k,79) + rxt(k,80) + rxt(k,81) &
                      + rxt(k,82)
         mat(k,1271) = rxt(k,83) + rxt(k,85) + rxt(k,86) + rxt(k,87) + rxt(k,88) &
                      + rxt(k,90) + rxt(k,91) + rxt(k,92)
         mat(k,1307) = rxt(k,93) + rxt(k,132)
         mat(k,406) = rxt(k,96)
         mat(k,1596) = rxt(k,98)
         mat(k,1632) = rxt(k,99) + rxt(k,353)
         mat(k,33) = -( het_rates(k,70) )
         mat(k,297) = -( het_rates(k,71) )
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
         mat(k,102) = -( rxt(k,368) + het_rates(k,72) )
         mat(k,810) = rxt(k,370)
         mat(k,829) = -( rxt(k,370) + het_rates(k,73) )
         mat(k,1189) = rxt(k,372)
         mat(k,22) = -( het_rates(k,74) )
         mat(k,1197) = -( rxt(k,372) + het_rates(k,75) )
         mat(k,1567) = rxt(k,374)
         mat(k,26) = -( het_rates(k,76) )
         mat(k,1576) = -( rxt(k,374) + het_rates(k,77) )
         mat(k,134) = -( het_rates(k,78) )
         mat(k,103) = rxt(k,368)
         mat(k,159) = -( het_rates(k,79) )
         mat(k,86) = -( het_rates(k,80) )
         mat(k,66) = rxt(k,578)
         mat(k,60) = rxt(k,583)
         mat(k,712) = -( rxt(k,75) + het_rates(k,81) )
         mat(k,372) = rxt(k,268)
         mat(k,371) = -( rxt(k,268) + het_rates(k,82) )
         mat(k,652) = -( rxt(k,76) + het_rates(k,83) )
         mat(k,439) = rxt(k,270)
         mat(k,520) = rxt(k,277)
         mat(k,477) = -( rxt(k,269) + het_rates(k,84) )
         mat(k,519) = -( rxt(k,277) + het_rates(k,85) )
         mat(k,478) = rxt(k,269)
         mat(k,385) = -( het_rates(k,86) )
         mat(k,438) = -( rxt(k,270) + het_rates(k,87) )
         mat(k,878) = -( rxt(k,365) + het_rates(k,88) )
         mat(k,1750) = rxt(k,16)
         mat(k,117) = rxt(k,364)
         mat(k,110) = rxt(k,366)
         mat(k,1349) = rxt(k,571)
         mat(k,249) = rxt(k,574)
         mat(k,922) = -( het_rates(k,89) )
         mat(k,45) = -( het_rates(k,90) )
         mat(k,116) = -( rxt(k,364) + het_rates(k,91) )
         mat(k,1066) = -( het_rates(k,92) )
         mat(k,108) = -( rxt(k,366) + het_rates(k,93) )
         mat(k,853) = rxt(k,365)
         mat(k,176) = -( het_rates(k,94) )
         mat(k,268) = rxt(k,66)
         mat(k,61) = rxt(k,584)
         mat(k,573) = -( rxt(k,386) + rxt(k,437)*y(k,49) + rxt(k,438)*y(k,49) &
                      + rxt(k,470)*y(k,5) + rxt(k,471)*y(k,6) + rxt(k,472)*y(k,8) &
                      + rxt(k,473)*y(k,9) + rxt(k,474)*y(k,10) + rxt(k,475)*y(k,11) &
                      + rxt(k,476)*y(k,12) + rxt(k,497)*y(k,7) + rxt(k,498)*y(k,21) &
                      + rxt(k,499)*y(k,34) + rxt(k,520)*y(k,13) + rxt(k,521)*y(k,15) &
                      + rxt(k,522)*y(k,37) + rxt(k,523)*y(k,38) + rxt(k,524)*y(k,39) &
                      + rxt(k,533)*y(k,20) + rxt(k,534)*y(k,20) + rxt(k,535)*y(k,20) &
                 + het_rates(k,95) )
         mat(k,1795) = rxt(k,1)
         mat(k,1260) = rxt(k,6)
         mat(k,1421) = rxt(k,7)
         mat(k,219) = -( rxt(k,382) + rxt(k,390) + het_rates(k,96) )
         mat(k,1415) = rxt(k,7)
         mat(k,31) = rxt(k,394)
         mat(k,30) = -( rxt(k,383) + rxt(k,394) + het_rates(k,97) )
         mat(k,1314) = -( rxt(k,93) + rxt(k,132) + het_rates(k,98) )
         mat(k,635) = rxt(k,65)
         mat(k,100) = rxt(k,97)
         mat(k,1360) = -( rxt(k,571) + het_rates(k,99) )
         mat(k,1279) = rxt(k,86) + rxt(k,88)
         mat(k,147) = rxt(k,94)
         mat(k,142) = -( rxt(k,94) + het_rates(k,100) )
         mat(k,401) = -( rxt(k,95) + rxt(k,96) + het_rates(k,101) )
         mat(k,95) = -( rxt(k,97) + het_rates(k,102) )
         mat(k,285) = -( het_rates(k,103) )
         mat(k,755) = -( rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) + rxt(k,512)*y(k,16) &
                      + rxt(k,514)*y(k,17) + rxt(k,516)*y(k,21) + rxt(k,517)*y(k,37) &
                      + rxt(k,518)*y(k,38) + rxt(k,519)*y(k,39) + rxt(k,531)*y(k,20) &
                 + het_rates(k,104) )
         mat(k,1800) = rxt(k,3)
         mat(k,81) = 2.000_r8*rxt(k,4)
         mat(k,1099) = rxt(k,9)
         mat(k,53) = rxt(k,10)
         mat(k,124) = rxt(k,12)
         mat(k,42) = rxt(k,23)
         mat(k,154) = rxt(k,58)
         mat(k,185) = rxt(k,59)
         mat(k,1590) = rxt(k,98)
         mat(k,1144) = .500_r8*rxt(k,539)
         mat(k,575) = rxt(k,533)*y(k,20)
         mat(k,1610) = -( rxt(k,98) + het_rates(k,105) )
         mat(k,1647) = -( rxt(k,99) + rxt(k,353) + het_rates(k,106) )
         mat(k,805) = rxt(k,63)
         mat(k,413) = rxt(k,95)
         mat(k,244) = -( rxt(k,574) + het_rates(k,107) )
         mat(k,960) = rxt(k,78) + rxt(k,79)
         mat(k,1256) = rxt(k,85) + rxt(k,87)
         mat(k,62) = rxt(k,558)
         mat(k,67) = rxt(k,559)
         mat(k,65) = -( rxt(k,559) + rxt(k,578) + het_rates(k,108) )
         mat(k,946) = rxt(k,80) + rxt(k,81)
         mat(k,1250) = rxt(k,90) + rxt(k,91)
         mat(k,59) = rxt(k,560)
         mat(k,58) = -( rxt(k,558) + rxt(k,560) + rxt(k,583) + rxt(k,584) &
                 + het_rates(k,109) )
         mat(k,945) = rxt(k,77) + rxt(k,82)
         mat(k,1249) = rxt(k,83) + rxt(k,92)
         mat(k,1825) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,110) )
         mat(k,515) = rxt(k,64) + rxt(k,278)
         mat(k,149) = rxt(k,94)
         mat(k,384) = rxt(k,268)
         mat(k,494) = rxt(k,269)
         mat(k,535) = rxt(k,277)
         mat(k,325) = rxt(k,279)
         mat(k,107) = rxt(k,368)
         mat(k,852) = rxt(k,370)
         mat(k,1212) = rxt(k,372)
         mat(k,1582) = rxt(k,374)
         mat(k,370) = rxt(k,380)
         mat(k,776) = rxt(k,509)*y(k,13) + rxt(k,511)*y(k,15) + rxt(k,512)*y(k,16) &
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
