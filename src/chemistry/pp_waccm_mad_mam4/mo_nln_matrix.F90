      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only: veclen
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,577) = rxt(k,492)*y(k,26)
         mat(k,1722) = rxt(k,492)*y(k,4)
         mat(k,1134) = (rxt(k,571)+rxt(k,576))*y(k,51)
         mat(k,206) = (rxt(k,571)+rxt(k,576))*y(k,47)
         mat(k,584) = -(4._r8*rxt(k,489)*y(k,4) + (rxt(k,490) + rxt(k,491) + rxt(k,492) &
                      ) * y(k,26) + rxt(k,493)*y(k,98) + rxt(k,494)*y(k,60) + rxt(k,495) &
                      *y(k,61) + rxt(k,497)*y(k,67) + rxt(k,498)*y(k,131) + rxt(k,546) &
                      *y(k,76))
         mat(k,1731) = -(rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,4)
         mat(k,744) = -rxt(k,493)*y(k,4)
         mat(k,1495) = -rxt(k,494)*y(k,4)
         mat(k,979) = -rxt(k,495)*y(k,4)
         mat(k,1411) = -rxt(k,497)*y(k,4)
         mat(k,790) = -rxt(k,498)*y(k,4)
         mat(k,381) = -rxt(k,546)*y(k,4)
         mat(k,141) = rxt(k,496)*y(k,67)
         mat(k,240) = rxt(k,506)*y(k,122)
         mat(k,209) = rxt(k,501)*y(k,67)
         mat(k,1411) = mat(k,1411) + rxt(k,496)*y(k,5) + rxt(k,501)*y(k,51)
         mat(k,1185) = rxt(k,488)*y(k,85)
         mat(k,455) = rxt(k,488)*y(k,69)
         mat(k,635) = rxt(k,506)*y(k,43)
         mat(k,138) = -(rxt(k,496)*y(k,67))
         mat(k,1391) = -rxt(k,496)*y(k,5)
         mat(k,579) = rxt(k,495)*y(k,61)
         mat(k,966) = rxt(k,495)*y(k,4)
         mat(k,606) = -(rxt(k,450)*y(k,86) + rxt(k,486)*y(k,85) + rxt(k,530)*y(k,62) &
                      + rxt(k,531)*y(k,67) + rxt(k,532)*y(k,131))
         mat(k,1454) = -rxt(k,450)*y(k,16)
         mat(k,456) = -rxt(k,486)*y(k,16)
         mat(k,1104) = -rxt(k,530)*y(k,16)
         mat(k,1412) = -rxt(k,531)*y(k,16)
         mat(k,791) = -rxt(k,532)*y(k,16)
         mat(k,325) = rxt(k,457)*y(k,26) + rxt(k,534)*y(k,60)
         mat(k,79) = .300_r8*rxt(k,535)*y(k,131)
         mat(k,1732) = rxt(k,457)*y(k,20)
         mat(k,1496) = rxt(k,534)*y(k,20)
         mat(k,791) = mat(k,791) + .300_r8*rxt(k,535)*y(k,21)
         mat(k,324) = -(rxt(k,457)*y(k,26) + rxt(k,533)*y(k,98) + rxt(k,534)*y(k,60))
         mat(k,1728) = -rxt(k,457)*y(k,20)
         mat(k,741) = -rxt(k,533)*y(k,20)
         mat(k,1489) = -rxt(k,534)*y(k,20)
         mat(k,78) = .700_r8*rxt(k,535)*y(k,131)
         mat(k,786) = .700_r8*rxt(k,535)*y(k,21)
         mat(k,77) = -(rxt(k,535)*y(k,131))
         mat(k,773) = -rxt(k,535)*y(k,21)
         mat(k,323) = rxt(k,533)*y(k,98)
         mat(k,734) = rxt(k,533)*y(k,20)
         mat(k,1721) = 2.000_r8*rxt(k,459)*y(k,26)
         mat(k,275) = (rxt(k,569)+rxt(k,574)+rxt(k,579))*y(k,47) + rxt(k,463)*y(k,86)
         mat(k,1133) = (rxt(k,569)+rxt(k,574)+rxt(k,579))*y(k,27) + (rxt(k,564) &
                       +rxt(k,570)+rxt(k,575))*y(k,52)
         mat(k,230) = (rxt(k,564)+rxt(k,570)+rxt(k,575))*y(k,47)
         mat(k,1444) = rxt(k,463)*y(k,27)
         mat(k,1720) = 2.000_r8*rxt(k,484)*y(k,26)
         mat(k,1761) = -(rxt(k,116)*y(k,91) + rxt(k,128)*y(k,94) + rxt(k,286)*y(k,108) &
                      + rxt(k,315)*y(k,125) + rxt(k,342)*y(k,132) + rxt(k,351) &
                      *y(k,133) + rxt(k,457)*y(k,20) + (4._r8*rxt(k,458) &
                      + 4._r8*rxt(k,459) + 4._r8*rxt(k,460) + 4._r8*rxt(k,484) &
                      ) * y(k,26) + rxt(k,461)*y(k,98) + rxt(k,462)*y(k,60) + rxt(k,464) &
                      *y(k,61) + rxt(k,467)*y(k,67) + (rxt(k,468) + rxt(k,469) &
                      ) * y(k,131) + (rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,4) &
                      + rxt(k,547)*y(k,76))
         mat(k,844) = -rxt(k,116)*y(k,26)
         mat(k,685) = -rxt(k,128)*y(k,26)
         mat(k,873) = -rxt(k,286)*y(k,26)
         mat(k,1338) = -rxt(k,315)*y(k,26)
         mat(k,1640) = -rxt(k,342)*y(k,26)
         mat(k,1675) = -rxt(k,351)*y(k,26)
         mat(k,332) = -rxt(k,457)*y(k,26)
         mat(k,765) = -rxt(k,461)*y(k,26)
         mat(k,1524) = -rxt(k,462)*y(k,26)
         mat(k,1008) = -rxt(k,464)*y(k,26)
         mat(k,1441) = -rxt(k,467)*y(k,26)
         mat(k,812) = -(rxt(k,468) + rxt(k,469)) * y(k,26)
         mat(k,598) = -(rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,26)
         mat(k,390) = -rxt(k,547)*y(k,26)
         mat(k,285) = rxt(k,465)*y(k,67)
         mat(k,1174) = rxt(k,483)*y(k,122)
         mat(k,236) = rxt(k,473)*y(k,67) + rxt(k,472)*y(k,86) + rxt(k,474)*y(k,131)
         mat(k,1441) = mat(k,1441) + rxt(k,465)*y(k,27) + rxt(k,473)*y(k,52)
         mat(k,1215) = rxt(k,456)*y(k,86)
         mat(k,68) = rxt(k,552)*y(k,76)
         mat(k,390) = mat(k,390) + rxt(k,552)*y(k,70)
         mat(k,1483) = rxt(k,472)*y(k,52) + rxt(k,456)*y(k,69) + rxt(k,455)*y(k,98)
         mat(k,765) = mat(k,765) + rxt(k,455)*y(k,86)
         mat(k,655) = rxt(k,483)*y(k,47)
         mat(k,812) = mat(k,812) + rxt(k,474)*y(k,52)
         mat(k,277) = -(rxt(k,463)*y(k,86) + rxt(k,465)*y(k,67) + rxt(k,466)*y(k,131) &
                      + (rxt(k,569) + rxt(k,574) + rxt(k,579)) * y(k,47))
         mat(k,1448) = -rxt(k,463)*y(k,27)
         mat(k,1403) = -rxt(k,465)*y(k,27)
         mat(k,784) = -rxt(k,466)*y(k,27)
         mat(k,1138) = -(rxt(k,569) + rxt(k,574) + rxt(k,579)) * y(k,27)
         mat(k,1726) = rxt(k,464)*y(k,61)
         mat(k,969) = rxt(k,464)*y(k,26)
         mat(k,134) = -((rxt(k,537) + rxt(k,541)) * y(k,131))
         mat(k,775) = -(rxt(k,537) + rxt(k,541)) * y(k,29)
         mat(k,601) = rxt(k,530)*y(k,62) + rxt(k,531)*y(k,67) + rxt(k,486)*y(k,85) &
                      + rxt(k,450)*y(k,86) + rxt(k,532)*y(k,131)
         mat(k,1095) = rxt(k,530)*y(k,16)
         mat(k,1390) = rxt(k,531)*y(k,16) + rxt(k,542)*y(k,71)
         mat(k,70) = rxt(k,542)*y(k,67) + rxt(k,543)*y(k,131)
         mat(k,452) = rxt(k,486)*y(k,16)
         mat(k,1445) = rxt(k,450)*y(k,16)
         mat(k,775) = mat(k,775) + rxt(k,532)*y(k,16) + rxt(k,543)*y(k,71)
         mat(k,27) = -(rxt(k,511)*y(k,122))
         mat(k,625) = -rxt(k,511)*y(k,31)
         mat(k,32) = -(rxt(k,512)*y(k,122))
         mat(k,626) = -rxt(k,512)*y(k,32)
         mat(k,58) = -(rxt(k,556)*y(k,62) + (rxt(k,557) + rxt(k,558)) * y(k,131))
         mat(k,1094) = -rxt(k,556)*y(k,33)
         mat(k,771) = -(rxt(k,557) + rxt(k,558)) * y(k,33)
         mat(k,198) = -(rxt(k,508)*y(k,39) + rxt(k,509)*y(k,137) + rxt(k,510)*y(k,49))
         mat(k,482) = -rxt(k,508)*y(k,37)
         mat(k,1812) = -rxt(k,509)*y(k,37)
         mat(k,1269) = -rxt(k,510)*y(k,37)
         mat(k,28) = 2.000_r8*rxt(k,511)*y(k,122)
         mat(k,33) = rxt(k,512)*y(k,122)
         mat(k,628) = 2.000_r8*rxt(k,511)*y(k,31) + rxt(k,512)*y(k,32)
         mat(k,1797) = -(rxt(k,105)*y(k,87) + rxt(k,117)*y(k,91) + rxt(k,129)*y(k,94) &
                      + rxt(k,287)*y(k,108) + rxt(k,309)*y(k,119) + rxt(k,317) &
                      *y(k,125) + rxt(k,331)*y(k,128) + rxt(k,344)*y(k,132) + (rxt(k,408) &
                      + rxt(k,409) + rxt(k,410)) * y(k,98) + rxt(k,411)*y(k,68) &
                      + rxt(k,414)*y(k,69))
         mat(k,707) = -rxt(k,105)*y(k,38)
         mat(k,845) = -rxt(k,117)*y(k,38)
         mat(k,686) = -rxt(k,129)*y(k,38)
         mat(k,874) = -rxt(k,287)*y(k,38)
         mat(k,1092) = -rxt(k,309)*y(k,38)
         mat(k,1339) = -rxt(k,317)*y(k,38)
         mat(k,451) = -rxt(k,331)*y(k,38)
         mat(k,1641) = -rxt(k,344)*y(k,38)
         mat(k,766) = -(rxt(k,408) + rxt(k,409) + rxt(k,410)) * y(k,38)
         mat(k,1265) = -rxt(k,411)*y(k,38)
         mat(k,1216) = -rxt(k,414)*y(k,38)
         mat(k,623) = rxt(k,532)*y(k,131)
         mat(k,137) = rxt(k,541)*y(k,131)
         mat(k,204) = rxt(k,508)*y(k,39)
         mat(k,501) = rxt(k,508)*y(k,37) + rxt(k,406)*y(k,67) + rxt(k,452)*y(k,86) &
                      + rxt(k,389)*y(k,122) + rxt(k,415)*y(k,131) + rxt(k,354) &
                      *y(k,133)
         mat(k,244) = rxt(k,506)*y(k,122)
         mat(k,1175) = rxt(k,483)*y(k,122)
         mat(k,322) = rxt(k,438)*y(k,131)
         mat(k,1442) = rxt(k,406)*y(k,39) + rxt(k,418)*y(k,131)
         mat(k,76) = rxt(k,543)*y(k,131)
         mat(k,174) = rxt(k,548)*y(k,131)
         mat(k,391) = rxt(k,553)*y(k,131)
         mat(k,1484) = rxt(k,452)*y(k,39)
         mat(k,707) = mat(k,707) + rxt(k,195)*y(k,100) + rxt(k,147)*y(k,102) &
                      + rxt(k,177)*y(k,104)
         mat(k,422) = rxt(k,199)*y(k,100) + rxt(k,164)*y(k,102) + rxt(k,182)*y(k,104)
         mat(k,406) = rxt(k,187)*y(k,100) + rxt(k,181)*y(k,102) + rxt(k,169)*y(k,104)
         mat(k,845) = mat(k,845) + rxt(k,186)*y(k,100) + (rxt(k,170)+rxt(k,258)) &
                      *y(k,102) + (rxt(k,168)+rxt(k,265))*y(k,104)
         mat(k,374) = rxt(k,194)*y(k,100) + (rxt(k,247)+rxt(k,271))*y(k,102) + ( &
                      + rxt(k,176)+rxt(k,259))*y(k,104)
         mat(k,555) = rxt(k,196)*y(k,100) + (rxt(k,158)+rxt(k,260))*y(k,102) + ( &
                      + rxt(k,178)+rxt(k,261))*y(k,104)
         mat(k,686) = mat(k,686) + rxt(k,191)*y(k,100) + rxt(k,225)*y(k,102) &
                      + rxt(k,174)*y(k,104)
         mat(k,1567) = rxt(k,138)*y(k,96) + rxt(k,382)*y(k,99) + rxt(k,383)*y(k,100) &
                      + rxt(k,141)*y(k,102) + rxt(k,144)*y(k,104) + rxt(k,381) &
                      *y(k,105)
         mat(k,87) = rxt(k,138)*y(k,95)
         mat(k,346) = rxt(k,189)*y(k,100) + rxt(k,203)*y(k,102) + rxt(k,172)*y(k,104)
         mat(k,166) = rxt(k,382)*y(k,95)
         mat(k,1051) = rxt(k,195)*y(k,87) + rxt(k,199)*y(k,88) + rxt(k,187)*y(k,89) &
                      + rxt(k,186)*y(k,91) + rxt(k,194)*y(k,92) + rxt(k,196)*y(k,93) &
                      + rxt(k,191)*y(k,94) + rxt(k,383)*y(k,95) + rxt(k,189)*y(k,97) &
                      + rxt(k,201)*y(k,108) + rxt(k,197)*y(k,109) + rxt(k,200) &
                      *y(k,111) + rxt(k,193)*y(k,112) + rxt(k,198)*y(k,113) &
                      + rxt(k,190)*y(k,125)
         mat(k,1718) = rxt(k,147)*y(k,87) + rxt(k,164)*y(k,88) + rxt(k,181)*y(k,89) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,91) + (rxt(k,247)+rxt(k,271)) &
                      *y(k,92) + (rxt(k,158)+rxt(k,260))*y(k,93) + rxt(k,225)*y(k,94) &
                      + rxt(k,141)*y(k,95) + rxt(k,203)*y(k,97) + rxt(k,166)*y(k,108) &
                      + rxt(k,162)*y(k,109) + rxt(k,165)*y(k,111) + (rxt(k,236) &
                       +rxt(k,262))*y(k,112) + rxt(k,163)*y(k,113) + rxt(k,214) &
                      *y(k,125)
         mat(k,1609) = rxt(k,177)*y(k,87) + rxt(k,182)*y(k,88) + rxt(k,169)*y(k,89) + ( &
                      + rxt(k,168)+rxt(k,265))*y(k,91) + (rxt(k,176)+rxt(k,259)) &
                      *y(k,92) + (rxt(k,178)+rxt(k,261))*y(k,93) + rxt(k,174)*y(k,94) &
                      + rxt(k,144)*y(k,95) + rxt(k,172)*y(k,97) + rxt(k,184)*y(k,108) &
                      + rxt(k,179)*y(k,109) + rxt(k,183)*y(k,111) + (rxt(k,175) &
                       +rxt(k,263))*y(k,112) + rxt(k,180)*y(k,113) + rxt(k,173) &
                      *y(k,125)
         mat(k,184) = rxt(k,381)*y(k,95)
         mat(k,874) = mat(k,874) + rxt(k,201)*y(k,100) + rxt(k,166)*y(k,102) &
                      + rxt(k,184)*y(k,104)
         mat(k,436) = rxt(k,197)*y(k,100) + rxt(k,162)*y(k,102) + rxt(k,179)*y(k,104)
         mat(k,535) = rxt(k,200)*y(k,100) + rxt(k,165)*y(k,102) + rxt(k,183)*y(k,104)
         mat(k,575) = rxt(k,193)*y(k,100) + (rxt(k,236)+rxt(k,262))*y(k,102) + ( &
                      + rxt(k,175)+rxt(k,263))*y(k,104)
         mat(k,479) = rxt(k,198)*y(k,100) + rxt(k,163)*y(k,102) + rxt(k,180)*y(k,104)
         mat(k,656) = rxt(k,389)*y(k,39) + rxt(k,506)*y(k,43) + rxt(k,483)*y(k,47)
         mat(k,1339) = mat(k,1339) + rxt(k,190)*y(k,100) + rxt(k,214)*y(k,102) &
                      + rxt(k,173)*y(k,104)
         mat(k,813) = rxt(k,532)*y(k,16) + rxt(k,541)*y(k,29) + rxt(k,415)*y(k,39) &
                      + rxt(k,438)*y(k,54) + rxt(k,418)*y(k,67) + rxt(k,543)*y(k,71) &
                      + rxt(k,548)*y(k,74) + rxt(k,553)*y(k,76)
         mat(k,1676) = rxt(k,354)*y(k,39)
      end do
      end subroutine nlnmat01
      subroutine nlnmat02( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,484) = -((rxt(k,353) + rxt(k,354)) * y(k,133) + rxt(k,389)*y(k,122) &
                      + rxt(k,406)*y(k,67) + rxt(k,415)*y(k,131) + rxt(k,452)*y(k,86) &
                      + rxt(k,508)*y(k,37))
         mat(k,1646) = -(rxt(k,353) + rxt(k,354)) * y(k,39)
         mat(k,634) = -rxt(k,389)*y(k,39)
         mat(k,1410) = -rxt(k,406)*y(k,39)
         mat(k,789) = -rxt(k,415)*y(k,39)
         mat(k,1452) = -rxt(k,452)*y(k,39)
         mat(k,200) = -rxt(k,508)*y(k,39)
         mat(k,1766) = rxt(k,408)*y(k,98)
         mat(k,743) = rxt(k,408)*y(k,38)
         mat(k,146) = -(rxt(k,407)*y(k,67) + rxt(k,416)*y(k,131) + rxt(k,453)*y(k,86))
         mat(k,1392) = -rxt(k,407)*y(k,41)
         mat(k,776) = -rxt(k,416)*y(k,41)
         mat(k,1446) = -rxt(k,453)*y(k,41)
         mat(k,736) = 2.000_r8*rxt(k,422)*y(k,98)
         mat(k,776) = mat(k,776) + 2.000_r8*rxt(k,421)*y(k,131)
         mat(k,45) = rxt(k,555)*y(k,137)
         mat(k,1799) = rxt(k,555)*y(k,78)
         mat(k,238) = -(rxt(k,499)*y(k,67) + rxt(k,500)*y(k,131) + (rxt(k,505) &
                      + rxt(k,506)) * y(k,122))
         mat(k,1399) = -rxt(k,499)*y(k,43)
         mat(k,782) = -rxt(k,500)*y(k,43)
         mat(k,629) = -(rxt(k,505) + rxt(k,506)) * y(k,43)
         mat(k,602) = rxt(k,486)*y(k,85)
         mat(k,453) = rxt(k,486)*y(k,16) + rxt(k,487)*y(k,98)
         mat(k,740) = rxt(k,487)*y(k,85)
         mat(k,1160) = -(rxt(k,106)*y(k,88) + rxt(k,108)*y(k,87) + rxt(k,130)*y(k,94) &
                      + (rxt(k,276) + rxt(k,298)) * y(k,110) + rxt(k,289)*y(k,108) &
                      + rxt(k,318)*y(k,125) + rxt(k,345)*y(k,132) + rxt(k,356) &
                      *y(k,133) + rxt(k,470)*y(k,67) + rxt(k,471)*y(k,131) + (rxt(k,482) &
                      + rxt(k,483)) * y(k,122) + (rxt(k,564) + rxt(k,570) + rxt(k,575) &
                      ) * y(k,52) + (rxt(k,569) + rxt(k,574) + rxt(k,579)) * y(k,27) &
                      + (rxt(k,571) + rxt(k,576)) * y(k,51))
         mat(k,415) = -rxt(k,106)*y(k,47)
         mat(k,698) = -rxt(k,108)*y(k,47)
         mat(k,672) = -rxt(k,130)*y(k,47)
         mat(k,721) = -(rxt(k,276) + rxt(k,298)) * y(k,47)
         mat(k,859) = -rxt(k,289)*y(k,47)
         mat(k,1324) = -rxt(k,318)*y(k,47)
         mat(k,1626) = -rxt(k,345)*y(k,47)
         mat(k,1661) = -rxt(k,356)*y(k,47)
         mat(k,1427) = -rxt(k,470)*y(k,47)
         mat(k,800) = -rxt(k,471)*y(k,47)
         mat(k,643) = -(rxt(k,482) + rxt(k,483)) * y(k,47)
         mat(k,233) = -(rxt(k,564) + rxt(k,570) + rxt(k,575)) * y(k,47)
         mat(k,281) = -(rxt(k,569) + rxt(k,574) + rxt(k,579)) * y(k,47)
         mat(k,211) = -(rxt(k,571) + rxt(k,576)) * y(k,47)
         mat(k,612) = rxt(k,450)*y(k,86)
         mat(k,1747) = rxt(k,469)*y(k,131)
         mat(k,1782) = rxt(k,105)*y(k,87)
         mat(k,491) = rxt(k,452)*y(k,86)
         mat(k,149) = rxt(k,453)*y(k,86)
         mat(k,1289) = rxt(k,109)*y(k,87) + rxt(k,299)*y(k,113)
         mat(k,233) = mat(k,233) + rxt(k,472)*y(k,86)
         mat(k,1469) = rxt(k,450)*y(k,16) + rxt(k,452)*y(k,39) + rxt(k,453)*y(k,41) &
                      + rxt(k,472)*y(k,52) + rxt(k,454)*y(k,98)
         mat(k,698) = mat(k,698) + rxt(k,105)*y(k,38) + rxt(k,109)*y(k,49)
         mat(k,399) = rxt(k,187)*y(k,100) + (rxt(k,181)+2.000_r8*rxt(k,267))*y(k,102) + ( &
                      + rxt(k,169)+2.000_r8*rxt(k,268))*y(k,104) + rxt(k,240)*y(k,115) &
                      + rxt(k,222)*y(k,116) + rxt(k,205)*y(k,119) + rxt(k,257) &
                      *y(k,126)
         mat(k,753) = rxt(k,454)*y(k,86)
         mat(k,1036) = rxt(k,187)*y(k,89) + rxt(k,198)*y(k,113)
         mat(k,1703) = (rxt(k,181)+2.000_r8*rxt(k,267))*y(k,89) + rxt(k,163)*y(k,113)
         mat(k,1594) = (rxt(k,169)+2.000_r8*rxt(k,268))*y(k,89) + rxt(k,180)*y(k,113)
         mat(k,472) = rxt(k,299)*y(k,49) + rxt(k,198)*y(k,100) + rxt(k,163)*y(k,102) &
                      + rxt(k,180)*y(k,104) + rxt(k,251)*y(k,115) + rxt(k,233) &
                      *y(k,116) + rxt(k,216)*y(k,119) + rxt(k,157)*y(k,126)
         mat(k,906) = rxt(k,240)*y(k,89) + rxt(k,251)*y(k,113)
         mat(k,948) = rxt(k,222)*y(k,89) + rxt(k,233)*y(k,113)
         mat(k,1077) = rxt(k,205)*y(k,89) + rxt(k,216)*y(k,113)
         mat(k,1368) = rxt(k,257)*y(k,89) + rxt(k,157)*y(k,113)
         mat(k,800) = mat(k,800) + rxt(k,469)*y(k,26)
         mat(k,197) = rxt(k,508)*y(k,39) + rxt(k,510)*y(k,49) + rxt(k,509)*y(k,137)
         mat(k,481) = rxt(k,508)*y(k,37)
         mat(k,1267) = rxt(k,510)*y(k,37)
         mat(k,1800) = rxt(k,509)*y(k,37)
         mat(k,1292) = -(rxt(k,109)*y(k,87) + rxt(k,124)*y(k,91) + rxt(k,290)*y(k,108) &
                      + rxt(k,295)*y(k,112) + rxt(k,299)*y(k,113) + rxt(k,300) &
                      *y(k,110) + rxt(k,319)*y(k,125) + rxt(k,357)*y(k,133) + rxt(k,447) &
                      *y(k,131) + rxt(k,510)*y(k,37))
         mat(k,700) = -rxt(k,109)*y(k,49)
         mat(k,833) = -rxt(k,124)*y(k,49)
         mat(k,862) = -rxt(k,290)*y(k,49)
         mat(k,569) = -rxt(k,295)*y(k,49)
         mat(k,474) = -rxt(k,299)*y(k,49)
         mat(k,724) = -rxt(k,300)*y(k,49)
         mat(k,1327) = -rxt(k,319)*y(k,49)
         mat(k,1664) = -rxt(k,357)*y(k,49)
         mat(k,803) = -rxt(k,447)*y(k,49)
         mat(k,203) = -rxt(k,510)*y(k,49)
         mat(k,615) = rxt(k,530)*y(k,62)
         mat(k,282) = (rxt(k,569)+rxt(k,574)+rxt(k,579))*y(k,47)
         mat(k,63) = rxt(k,556)*y(k,62)
         mat(k,1163) = (rxt(k,569)+rxt(k,574)+rxt(k,579))*y(k,27) + rxt(k,298) &
                      *y(k,110)
         mat(k,303) = rxt(k,142)*y(k,102) + rxt(k,145)*y(k,104) + rxt(k,293)*y(k,111) &
                      + rxt(k,297)*y(k,112)
         mat(k,997) = rxt(k,446)*y(k,131)
         mat(k,1119) = rxt(k,530)*y(k,16) + rxt(k,556)*y(k,33)
         mat(k,1039) = rxt(k,188)*y(k,110) + 2.000_r8*rxt(k,185)*y(k,114)
         mat(k,51) = rxt(k,140)*y(k,137)
         mat(k,1706) = rxt(k,142)*y(k,56) + (rxt(k,192)+rxt(k,264))*y(k,110) + ( &
                      + 2.000_r8*rxt(k,146)+2.000_r8*rxt(k,269))*y(k,114)
         mat(k,55) = rxt(k,143)*y(k,137)
         mat(k,1597) = rxt(k,145)*y(k,56) + (rxt(k,171)+rxt(k,266))*y(k,110) + ( &
                      + 2.000_r8*rxt(k,167)+2.000_r8*rxt(k,270))*y(k,114)
         mat(k,724) = mat(k,724) + rxt(k,298)*y(k,47) + rxt(k,188)*y(k,100) + ( &
                      + rxt(k,192)+rxt(k,264))*y(k,102) + (rxt(k,171)+rxt(k,266)) &
                      *y(k,104)
         mat(k,529) = rxt(k,293)*y(k,56)
         mat(k,569) = mat(k,569) + rxt(k,297)*y(k,56)
         mat(k,511) = 2.000_r8*rxt(k,185)*y(k,100) + (2.000_r8*rxt(k,146) &
                       +2.000_r8*rxt(k,269))*y(k,102) + (2.000_r8*rxt(k,167) &
                       +2.000_r8*rxt(k,270))*y(k,104) + rxt(k,238)*y(k,115) &
                      + rxt(k,220)*y(k,116) + rxt(k,202)*y(k,119) + rxt(k,255) &
                      *y(k,126)
         mat(k,909) = rxt(k,238)*y(k,114)
         mat(k,951) = rxt(k,220)*y(k,114)
         mat(k,1080) = rxt(k,202)*y(k,114)
         mat(k,1371) = rxt(k,255)*y(k,114)
         mat(k,803) = mat(k,803) + rxt(k,446)*y(k,61)
         mat(k,1843) = rxt(k,140)*y(k,101) + rxt(k,143)*y(k,103)
         mat(k,101) = -(rxt(k,423)*y(k,131))
         mat(k,774) = -rxt(k,423)*y(k,50)
         mat(k,965) = rxt(k,444)*y(k,98)
         mat(k,735) = rxt(k,444)*y(k,61)
         mat(k,207) = -(rxt(k,501)*y(k,67) + (rxt(k,571) + rxt(k,576)) * y(k,47))
         mat(k,1395) = -rxt(k,501)*y(k,51)
         mat(k,1136) = -(rxt(k,571) + rxt(k,576)) * y(k,51)
         mat(k,580) = rxt(k,493)*y(k,98)
         mat(k,738) = rxt(k,493)*y(k,4)
         mat(k,231) = -(rxt(k,472)*y(k,86) + rxt(k,473)*y(k,67) + rxt(k,474)*y(k,131) &
                      + (rxt(k,564) + rxt(k,570) + rxt(k,575)) * y(k,47))
         mat(k,1447) = -rxt(k,472)*y(k,52)
         mat(k,1398) = -rxt(k,473)*y(k,52)
         mat(k,781) = -rxt(k,474)*y(k,52)
         mat(k,1137) = -(rxt(k,564) + rxt(k,570) + rxt(k,575)) * y(k,52)
         mat(k,1724) = rxt(k,461)*y(k,98)
         mat(k,276) = rxt(k,466)*y(k,131)
         mat(k,739) = rxt(k,461)*y(k,26)
         mat(k,781) = mat(k,781) + rxt(k,466)*y(k,27)
         mat(k,175) = -(rxt(k,340)*y(k,131))
         mat(k,778) = -rxt(k,340)*y(k,53)
         mat(k,1135) = rxt(k,289)*y(k,108)
         mat(k,1268) = rxt(k,290)*y(k,108)
         mat(k,1486) = rxt(k,349)*y(k,131)
         mat(k,847) = rxt(k,289)*y(k,47) + rxt(k,290)*y(k,49)
         mat(k,90) = rxt(k,305)*y(k,137)
         mat(k,778) = mat(k,778) + rxt(k,349)*y(k,60)
         mat(k,1809) = rxt(k,305)*y(k,117)
         mat(k,311) = -(rxt(k,426)*y(k,60) + (rxt(k,427) + rxt(k,428) + rxt(k,429) &
                      ) * y(k,61) + rxt(k,430)*y(k,68) + rxt(k,438)*y(k,131) + rxt(k,592) &
                      *y(k,126))
         mat(k,1488) = -rxt(k,426)*y(k,54)
         mat(k,971) = -(rxt(k,427) + rxt(k,428) + rxt(k,429)) * y(k,54)
         mat(k,1229) = -rxt(k,430)*y(k,54)
         mat(k,785) = -rxt(k,438)*y(k,54)
         mat(k,1342) = -rxt(k,592)*y(k,54)
         mat(k,1405) = rxt(k,424)*y(k,106) + rxt(k,589)*y(k,121)
         mat(k,1229) = mat(k,1229) + rxt(k,590)*y(k,121)
         mat(k,1541) = 1.100_r8*rxt(k,585)*y(k,107) + .200_r8*rxt(k,583)*y(k,115)
         mat(k,218) = rxt(k,424)*y(k,67)
         mat(k,156) = 1.100_r8*rxt(k,585)*y(k,95)
         mat(k,881) = .200_r8*rxt(k,583)*y(k,95)
         mat(k,225) = rxt(k,589)*y(k,67) + rxt(k,590)*y(k,68)
         mat(k,297) = -(rxt(k,142)*y(k,102) + rxt(k,145)*y(k,104) + rxt(k,293) &
                      *y(k,111) + rxt(k,297)*y(k,112))
         mat(k,1679) = -rxt(k,142)*y(k,56)
         mat(k,1570) = -rxt(k,145)*y(k,56)
         mat(k,517) = -rxt(k,293)*y(k,56)
         mat(k,557) = -rxt(k,297)*y(k,56)
         mat(k,970) = rxt(k,445)*y(k,62)
         mat(k,1097) = rxt(k,445)*y(k,61)
         mat(k,1518) = -((rxt(k,111) + rxt(k,112)) * y(k,90) + (rxt(k,122) + rxt(k,123) &
                      ) * y(k,93) + rxt(k,136)*y(k,133) + (rxt(k,272) + rxt(k,279) &
                      ) * y(k,128) + rxt(k,280)*y(k,91) + rxt(k,349)*y(k,131) &
                      + rxt(k,426)*y(k,54) + rxt(k,435)*y(k,62) + rxt(k,439)*y(k,98) &
                      + rxt(k,440)*y(k,69) + rxt(k,441)*y(k,67) + rxt(k,462)*y(k,26) &
                      + rxt(k,494)*y(k,4) + rxt(k,534)*y(k,20) + rxt(k,594)*y(k,126))
         mat(k,294) = -(rxt(k,111) + rxt(k,112)) * y(k,60)
         mat(k,552) = -(rxt(k,122) + rxt(k,123)) * y(k,60)
         mat(k,1669) = -rxt(k,136)*y(k,60)
         mat(k,447) = -(rxt(k,272) + rxt(k,279)) * y(k,60)
         mat(k,838) = -rxt(k,280)*y(k,60)
         mat(k,808) = -rxt(k,349)*y(k,60)
         mat(k,320) = -rxt(k,426)*y(k,60)
         mat(k,1124) = -rxt(k,435)*y(k,60)
         mat(k,761) = -rxt(k,439)*y(k,60)
         mat(k,1209) = -rxt(k,440)*y(k,60)
         mat(k,1435) = -rxt(k,441)*y(k,60)
         mat(k,1755) = -rxt(k,462)*y(k,60)
         mat(k,597) = -rxt(k,494)*y(k,60)
         mat(k,331) = -rxt(k,534)*y(k,60)
         mat(k,1376) = -rxt(k,594)*y(k,60)
         mat(k,1790) = rxt(k,287)*y(k,108) + rxt(k,309)*y(k,119)
         mat(k,320) = mat(k,320) + 2.000_r8*rxt(k,428)*y(k,61) + rxt(k,430)*y(k,68) &
                      + rxt(k,438)*y(k,131)
         mat(k,1002) = 2.000_r8*rxt(k,428)*y(k,54) + rxt(k,431)*y(k,67) + rxt(k,549) &
                      *y(k,76) + rxt(k,291)*y(k,108)
         mat(k,1435) = mat(k,1435) + rxt(k,431)*y(k,61)
         mat(k,1258) = rxt(k,430)*y(k,54) + rxt(k,425)*y(k,106)
         mat(k,389) = rxt(k,549)*y(k,61)
         mat(k,703) = rxt(k,248)*y(k,115) + rxt(k,230)*y(k,116) + rxt(k,212)*y(k,119)
         mat(k,419) = rxt(k,252)*y(k,115) + rxt(k,234)*y(k,116) + rxt(k,217)*y(k,119)
         mat(k,403) = rxt(k,240)*y(k,115) + rxt(k,222)*y(k,116) + rxt(k,205)*y(k,119)
         mat(k,838) = mat(k,838) + rxt(k,239)*y(k,115) + rxt(k,221)*y(k,116) &
                      + rxt(k,204)*y(k,119)
         mat(k,371) = rxt(k,246)*y(k,115) + rxt(k,229)*y(k,116) + rxt(k,211)*y(k,119)
         mat(k,552) = mat(k,552) + rxt(k,249)*y(k,115) + rxt(k,231)*y(k,116) &
                      + rxt(k,213)*y(k,119)
         mat(k,679) = rxt(k,244)*y(k,115) + rxt(k,227)*y(k,116) + rxt(k,209)*y(k,119)
         mat(k,1560) = rxt(k,303)*y(k,116) + rxt(k,304)*y(k,117) + rxt(k,306)*y(k,118) &
                      + rxt(k,308)*y(k,119) + rxt(k,384)*y(k,120)
         mat(k,343) = rxt(k,242)*y(k,115) + rxt(k,224)*y(k,116) + rxt(k,207)*y(k,119)
         mat(k,222) = rxt(k,425)*y(k,68)
         mat(k,867) = rxt(k,287)*y(k,38) + rxt(k,291)*y(k,61) + rxt(k,254)*y(k,115) &
                      + rxt(k,237)*y(k,116) + rxt(k,219)*y(k,119)
         mat(k,433) = rxt(k,250)*y(k,115) + rxt(k,232)*y(k,116) + rxt(k,215)*y(k,119)
         mat(k,728) = rxt(k,241)*y(k,115) + rxt(k,223)*y(k,116) + rxt(k,206)*y(k,119)
         mat(k,532) = rxt(k,253)*y(k,115) + rxt(k,235)*y(k,116) + rxt(k,218)*y(k,119)
         mat(k,572) = rxt(k,245)*y(k,115) + rxt(k,228)*y(k,116) + rxt(k,210)*y(k,119)
         mat(k,476) = rxt(k,251)*y(k,115) + rxt(k,233)*y(k,116) + rxt(k,216)*y(k,119)
         mat(k,513) = rxt(k,238)*y(k,115) + rxt(k,220)*y(k,116) + rxt(k,202)*y(k,119)
         mat(k,914) = rxt(k,248)*y(k,87) + rxt(k,252)*y(k,88) + rxt(k,240)*y(k,89) &
                      + rxt(k,239)*y(k,91) + rxt(k,246)*y(k,92) + rxt(k,249)*y(k,93) &
                      + rxt(k,244)*y(k,94) + rxt(k,242)*y(k,97) + rxt(k,254)*y(k,108) &
                      + rxt(k,250)*y(k,109) + rxt(k,241)*y(k,110) + rxt(k,253) &
                      *y(k,111) + rxt(k,245)*y(k,112) + rxt(k,251)*y(k,113) &
                      + rxt(k,238)*y(k,114) + rxt(k,243)*y(k,125)
         mat(k,956) = rxt(k,230)*y(k,87) + rxt(k,234)*y(k,88) + rxt(k,222)*y(k,89) &
                      + rxt(k,221)*y(k,91) + rxt(k,229)*y(k,92) + rxt(k,231)*y(k,93) &
                      + rxt(k,227)*y(k,94) + rxt(k,303)*y(k,95) + rxt(k,224)*y(k,97) &
                      + rxt(k,237)*y(k,108) + rxt(k,232)*y(k,109) + rxt(k,223) &
                      *y(k,110) + rxt(k,235)*y(k,111) + rxt(k,228)*y(k,112) &
                      + rxt(k,233)*y(k,113) + rxt(k,220)*y(k,114) + rxt(k,226) &
                      *y(k,125)
         mat(k,92) = rxt(k,304)*y(k,95)
         mat(k,118) = rxt(k,306)*y(k,95)
         mat(k,1085) = rxt(k,309)*y(k,38) + rxt(k,212)*y(k,87) + rxt(k,217)*y(k,88) &
                      + rxt(k,205)*y(k,89) + rxt(k,204)*y(k,91) + rxt(k,211)*y(k,92) &
                      + rxt(k,213)*y(k,93) + rxt(k,209)*y(k,94) + rxt(k,308)*y(k,95) &
                      + rxt(k,207)*y(k,97) + rxt(k,219)*y(k,108) + rxt(k,215)*y(k,109) &
                      + rxt(k,206)*y(k,110) + rxt(k,218)*y(k,111) + rxt(k,210) &
                      *y(k,112) + rxt(k,216)*y(k,113) + rxt(k,202)*y(k,114) &
                      + rxt(k,208)*y(k,125)
         mat(k,112) = rxt(k,384)*y(k,95)
         mat(k,1332) = rxt(k,243)*y(k,115) + rxt(k,226)*y(k,116) + rxt(k,208)*y(k,119)
         mat(k,808) = mat(k,808) + rxt(k,438)*y(k,54)
      end do
      end subroutine nlnmat02
      subroutine nlnmat03( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,990) = -(rxt(k,110)*y(k,87) + (rxt(k,120) + rxt(k,121)) * y(k,93) &
                      + (rxt(k,277) + rxt(k,278)) * y(k,128) + rxt(k,281)*y(k,91) &
                      + rxt(k,291)*y(k,108) + rxt(k,320)*y(k,125) + rxt(k,346) &
                      *y(k,132) + rxt(k,359)*y(k,133) + (rxt(k,427) + rxt(k,428) &
                      + rxt(k,429)) * y(k,54) + (rxt(k,431) + rxt(k,433)) * y(k,67) &
                      + rxt(k,432)*y(k,69) + rxt(k,444)*y(k,98) + rxt(k,445)*y(k,62) &
                      + rxt(k,446)*y(k,131) + rxt(k,464)*y(k,26) + rxt(k,495)*y(k,4) &
                      + rxt(k,549)*y(k,76))
         mat(k,695) = -rxt(k,110)*y(k,61)
         mat(k,546) = -(rxt(k,120) + rxt(k,121)) * y(k,61)
         mat(k,442) = -(rxt(k,277) + rxt(k,278)) * y(k,61)
         mat(k,826) = -rxt(k,281)*y(k,61)
         mat(k,855) = -rxt(k,291)*y(k,61)
         mat(k,1320) = -rxt(k,320)*y(k,61)
         mat(k,1622) = -rxt(k,346)*y(k,61)
         mat(k,1657) = -rxt(k,359)*y(k,61)
         mat(k,316) = -(rxt(k,427) + rxt(k,428) + rxt(k,429)) * y(k,61)
         mat(k,1423) = -(rxt(k,431) + rxt(k,433)) * y(k,61)
         mat(k,1197) = -rxt(k,432)*y(k,61)
         mat(k,749) = -rxt(k,444)*y(k,61)
         mat(k,1112) = -rxt(k,445)*y(k,61)
         mat(k,796) = -rxt(k,446)*y(k,61)
         mat(k,1743) = -rxt(k,464)*y(k,61)
         mat(k,589) = -rxt(k,495)*y(k,61)
         mat(k,384) = -rxt(k,549)*y(k,61)
         mat(k,589) = mat(k,589) + rxt(k,494)*y(k,60)
         mat(k,328) = rxt(k,534)*y(k,60)
         mat(k,1743) = mat(k,1743) + rxt(k,462)*y(k,60)
         mat(k,104) = rxt(k,423)*y(k,131)
         mat(k,177) = rxt(k,340)*y(k,131)
         mat(k,1506) = rxt(k,494)*y(k,4) + rxt(k,534)*y(k,20) + rxt(k,462)*y(k,26) &
                      + 2.000_r8*rxt(k,435)*y(k,62) + rxt(k,441)*y(k,67) + rxt(k,440) &
                      *y(k,69) + rxt(k,112)*y(k,90) + rxt(k,439)*y(k,98) + rxt(k,136) &
                      *y(k,133)
         mat(k,1112) = mat(k,1112) + 2.000_r8*rxt(k,435)*y(k,60) + rxt(k,436)*y(k,67) &
                      + rxt(k,434)*y(k,98) + rxt(k,437)*y(k,131)
         mat(k,1423) = mat(k,1423) + rxt(k,441)*y(k,60) + rxt(k,436)*y(k,62)
         mat(k,1197) = mat(k,1197) + rxt(k,440)*y(k,60)
         mat(k,1465) = rxt(k,285)*y(k,108)
         mat(k,290) = rxt(k,112)*y(k,60)
         mat(k,749) = mat(k,749) + rxt(k,439)*y(k,60) + rxt(k,434)*y(k,62)
         mat(k,1032) = rxt(k,201)*y(k,108) + rxt(k,197)*y(k,109)
         mat(k,1699) = rxt(k,166)*y(k,108) + rxt(k,162)*y(k,109)
         mat(k,1590) = rxt(k,184)*y(k,108) + rxt(k,179)*y(k,109)
         mat(k,855) = mat(k,855) + rxt(k,285)*y(k,86) + rxt(k,201)*y(k,100) &
                      + rxt(k,166)*y(k,102) + rxt(k,184)*y(k,104) + rxt(k,254) &
                      *y(k,115) + rxt(k,237)*y(k,116) + rxt(k,219)*y(k,119) &
                      + rxt(k,161)*y(k,126)
         mat(k,428) = rxt(k,197)*y(k,100) + rxt(k,162)*y(k,102) + rxt(k,179)*y(k,104) &
                      + rxt(k,250)*y(k,115) + rxt(k,232)*y(k,116) + rxt(k,215) &
                      *y(k,119) + rxt(k,156)*y(k,126)
         mat(k,902) = rxt(k,254)*y(k,108) + rxt(k,250)*y(k,109)
         mat(k,944) = rxt(k,237)*y(k,108) + rxt(k,232)*y(k,109)
         mat(k,1073) = rxt(k,219)*y(k,108) + rxt(k,215)*y(k,109) + rxt(k,311)*y(k,131)
         mat(k,1364) = rxt(k,161)*y(k,108) + rxt(k,156)*y(k,109)
         mat(k,796) = mat(k,796) + rxt(k,423)*y(k,50) + rxt(k,340)*y(k,53) &
                      + rxt(k,437)*y(k,62) + rxt(k,311)*y(k,119)
         mat(k,1657) = mat(k,1657) + rxt(k,136)*y(k,60)
         mat(k,1115) = -(rxt(k,434)*y(k,98) + rxt(k,435)*y(k,60) + rxt(k,436)*y(k,67) &
                      + rxt(k,437)*y(k,131) + rxt(k,445)*y(k,61) + rxt(k,530)*y(k,16) &
                      + rxt(k,556)*y(k,33))
         mat(k,752) = -rxt(k,434)*y(k,62)
         mat(k,1509) = -rxt(k,435)*y(k,62)
         mat(k,1426) = -rxt(k,436)*y(k,62)
         mat(k,799) = -rxt(k,437)*y(k,62)
         mat(k,993) = -rxt(k,445)*y(k,62)
         mat(k,611) = -rxt(k,530)*y(k,62)
         mat(k,62) = -rxt(k,556)*y(k,62)
         mat(k,143) = rxt(k,496)*y(k,67)
         mat(k,1746) = rxt(k,286)*y(k,108)
         mat(k,280) = rxt(k,465)*y(k,67) + rxt(k,463)*y(k,86) + rxt(k,466)*y(k,131)
         mat(k,202) = rxt(k,510)*y(k,49)
         mat(k,1288) = rxt(k,510)*y(k,37) + rxt(k,447)*y(k,131)
         mat(k,993) = mat(k,993) + rxt(k,433)*y(k,67) + rxt(k,432)*y(k,69)
         mat(k,1426) = mat(k,1426) + rxt(k,496)*y(k,5) + rxt(k,465)*y(k,27) &
                      + rxt(k,433)*y(k,61)
         mat(k,1200) = rxt(k,432)*y(k,61)
         mat(k,1468) = rxt(k,463)*y(k,27)
         mat(k,752) = mat(k,752) + rxt(k,310)*y(k,119)
         mat(k,1035) = rxt(k,200)*y(k,111) + rxt(k,193)*y(k,112) + rxt(k,198)*y(k,113)
         mat(k,1702) = rxt(k,165)*y(k,111) + (rxt(k,236)+rxt(k,262))*y(k,112) &
                      + rxt(k,163)*y(k,113)
         mat(k,1593) = rxt(k,183)*y(k,111) + (rxt(k,175)+rxt(k,263))*y(k,112) &
                      + rxt(k,180)*y(k,113)
         mat(k,858) = rxt(k,286)*y(k,26)
         mat(k,720) = rxt(k,241)*y(k,115) + rxt(k,223)*y(k,116) + rxt(k,206)*y(k,119) &
                      + rxt(k,148)*y(k,126)
         mat(k,527) = rxt(k,200)*y(k,100) + rxt(k,165)*y(k,102) + rxt(k,183)*y(k,104) &
                      + rxt(k,253)*y(k,115) + rxt(k,235)*y(k,116) + rxt(k,218) &
                      *y(k,119) + rxt(k,160)*y(k,126)
         mat(k,567) = rxt(k,193)*y(k,100) + (rxt(k,236)+rxt(k,262))*y(k,102) + ( &
                      + rxt(k,175)+rxt(k,263))*y(k,104) + rxt(k,245)*y(k,115) &
                      + rxt(k,228)*y(k,116) + rxt(k,210)*y(k,119) + rxt(k,152) &
                      *y(k,126)
         mat(k,471) = rxt(k,198)*y(k,100) + rxt(k,163)*y(k,102) + rxt(k,180)*y(k,104) &
                      + rxt(k,251)*y(k,115) + rxt(k,233)*y(k,116) + rxt(k,216) &
                      *y(k,119) + rxt(k,157)*y(k,126)
         mat(k,509) = rxt(k,238)*y(k,115) + rxt(k,220)*y(k,116) + rxt(k,202)*y(k,119) &
                      + rxt(k,255)*y(k,126)
         mat(k,905) = rxt(k,241)*y(k,110) + rxt(k,253)*y(k,111) + rxt(k,245)*y(k,112) &
                      + rxt(k,251)*y(k,113) + rxt(k,238)*y(k,114)
         mat(k,947) = rxt(k,223)*y(k,110) + rxt(k,235)*y(k,111) + rxt(k,228)*y(k,112) &
                      + rxt(k,233)*y(k,113) + rxt(k,220)*y(k,114)
         mat(k,1076) = rxt(k,310)*y(k,98) + rxt(k,206)*y(k,110) + rxt(k,218)*y(k,111) &
                      + rxt(k,210)*y(k,112) + rxt(k,216)*y(k,113) + rxt(k,202) &
                      *y(k,114)
         mat(k,1367) = rxt(k,148)*y(k,110) + rxt(k,160)*y(k,111) + rxt(k,152)*y(k,112) &
                      + rxt(k,157)*y(k,113) + rxt(k,255)*y(k,114)
         mat(k,799) = mat(k,799) + rxt(k,466)*y(k,27) + rxt(k,447)*y(k,49)
         mat(k,1433) = -(rxt(k,113)*y(k,90) + rxt(k,125)*y(k,91) + rxt(k,131)*y(k,94) &
                      + rxt(k,301)*y(k,110) + (rxt(k,324) + rxt(k,325)) * y(k,125) &
                      + (rxt(k,333) + rxt(k,334)) * y(k,128) + rxt(k,336)*y(k,129) &
                      + rxt(k,338)*y(k,130) + rxt(k,347)*y(k,132) + rxt(k,360) &
                      *y(k,133) + rxt(k,403)*y(k,69) + 4._r8*rxt(k,404)*y(k,67) &
                      + rxt(k,405)*y(k,68) + rxt(k,406)*y(k,39) + rxt(k,407)*y(k,41) &
                      + rxt(k,412)*y(k,98) + rxt(k,418)*y(k,131) + (rxt(k,431) &
                      + rxt(k,433)) * y(k,61) + rxt(k,436)*y(k,62) + rxt(k,441) &
                      *y(k,60) + rxt(k,465)*y(k,27) + rxt(k,467)*y(k,26) + rxt(k,470) &
                      *y(k,47) + rxt(k,473)*y(k,52) + rxt(k,496)*y(k,5) + rxt(k,497) &
                      *y(k,4) + rxt(k,499)*y(k,43) + rxt(k,501)*y(k,51) + rxt(k,531) &
                      *y(k,16) + rxt(k,542)*y(k,71) + (rxt(k,587) + rxt(k,588) &
                      ) * y(k,107) + rxt(k,589)*y(k,121))
         mat(k,292) = -rxt(k,113)*y(k,67)
         mat(k,836) = -rxt(k,125)*y(k,67)
         mat(k,677) = -rxt(k,131)*y(k,67)
         mat(k,726) = -rxt(k,301)*y(k,67)
         mat(k,1330) = -(rxt(k,324) + rxt(k,325)) * y(k,67)
         mat(k,446) = -(rxt(k,333) + rxt(k,334)) * y(k,67)
         mat(k,100) = -rxt(k,336)*y(k,67)
         mat(k,357) = -rxt(k,338)*y(k,67)
         mat(k,1632) = -rxt(k,347)*y(k,67)
         mat(k,1667) = -rxt(k,360)*y(k,67)
         mat(k,1207) = -rxt(k,403)*y(k,67)
         mat(k,1256) = -rxt(k,405)*y(k,67)
         mat(k,494) = -rxt(k,406)*y(k,67)
         mat(k,150) = -rxt(k,407)*y(k,67)
         mat(k,759) = -rxt(k,412)*y(k,67)
         mat(k,806) = -rxt(k,418)*y(k,67)
         mat(k,1000) = -(rxt(k,431) + rxt(k,433)) * y(k,67)
         mat(k,1122) = -rxt(k,436)*y(k,67)
         mat(k,1516) = -rxt(k,441)*y(k,67)
         mat(k,283) = -rxt(k,465)*y(k,67)
         mat(k,1753) = -rxt(k,467)*y(k,67)
         mat(k,1166) = -rxt(k,470)*y(k,67)
         mat(k,234) = -rxt(k,473)*y(k,67)
         mat(k,145) = -rxt(k,496)*y(k,67)
         mat(k,595) = -rxt(k,497)*y(k,67)
         mat(k,243) = -rxt(k,499)*y(k,67)
         mat(k,212) = -rxt(k,501)*y(k,67)
         mat(k,616) = -rxt(k,531)*y(k,67)
         mat(k,75) = -rxt(k,542)*y(k,67)
         mat(k,160) = -(rxt(k,587) + rxt(k,588)) * y(k,67)
         mat(k,229) = -rxt(k,589)*y(k,67)
         mat(k,1788) = rxt(k,410)*y(k,98)
         mat(k,319) = rxt(k,426)*y(k,60) + rxt(k,427)*y(k,61) + rxt(k,430)*y(k,68) &
                      + rxt(k,592)*y(k,126)
         mat(k,1516) = mat(k,1516) + rxt(k,426)*y(k,54) + rxt(k,272)*y(k,128)
         mat(k,1000) = mat(k,1000) + rxt(k,427)*y(k,54) + rxt(k,359)*y(k,133)
         mat(k,1256) = mat(k,1256) + rxt(k,430)*y(k,54) + rxt(k,544)*y(k,74) &
                      + rxt(k,550)*y(k,76) + rxt(k,591)*y(k,121) + (rxt(k,392) &
                       +rxt(k,393))*y(k,122) + rxt(k,598)*y(k,134) + rxt(k,602) &
                      *y(k,135)
         mat(k,1207) = mat(k,1207) + rxt(k,363)*y(k,133)
         mat(k,173) = rxt(k,544)*y(k,68)
         mat(k,387) = rxt(k,550)*y(k,68)
         mat(k,1475) = rxt(k,114)*y(k,91) + rxt(k,350)*y(k,133)
         mat(k,836) = mat(k,836) + rxt(k,114)*y(k,86) + rxt(k,186)*y(k,100) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,102) + (rxt(k,168)+rxt(k,265)) &
                      *y(k,104) + rxt(k,239)*y(k,115) + rxt(k,221)*y(k,116) &
                      + rxt(k,204)*y(k,119) + rxt(k,256)*y(k,126)
         mat(k,370) = rxt(k,194)*y(k,100) + (rxt(k,247)+rxt(k,271))*y(k,102) + ( &
                      + rxt(k,176)+rxt(k,259))*y(k,104) + rxt(k,246)*y(k,115) &
                      + rxt(k,229)*y(k,116) + rxt(k,211)*y(k,119) + rxt(k,153) &
                      *y(k,126)
         mat(k,551) = rxt(k,196)*y(k,100) + (rxt(k,158)+rxt(k,260))*y(k,102) + ( &
                      + rxt(k,178)+rxt(k,261))*y(k,104) + rxt(k,249)*y(k,115) &
                      + rxt(k,231)*y(k,116) + rxt(k,213)*y(k,119) + rxt(k,155) &
                      *y(k,126)
         mat(k,1558) = rxt(k,583)*y(k,115) + 1.150_r8*rxt(k,584)*y(k,126)
         mat(k,759) = mat(k,759) + rxt(k,410)*y(k,38)
         mat(k,1042) = rxt(k,186)*y(k,91) + rxt(k,194)*y(k,92) + rxt(k,196)*y(k,93)
         mat(k,1709) = (rxt(k,170)+rxt(k,258))*y(k,91) + (rxt(k,247)+rxt(k,271)) &
                      *y(k,92) + (rxt(k,158)+rxt(k,260))*y(k,93)
         mat(k,1600) = (rxt(k,168)+rxt(k,265))*y(k,91) + (rxt(k,176)+rxt(k,259)) &
                      *y(k,92) + (rxt(k,178)+rxt(k,261))*y(k,93)
         mat(k,221) = rxt(k,597)*y(k,134)
         mat(k,912) = rxt(k,239)*y(k,91) + rxt(k,246)*y(k,92) + rxt(k,249)*y(k,93) &
                      + rxt(k,583)*y(k,95)
         mat(k,954) = rxt(k,221)*y(k,91) + rxt(k,229)*y(k,92) + rxt(k,231)*y(k,93)
         mat(k,1083) = rxt(k,204)*y(k,91) + rxt(k,211)*y(k,92) + rxt(k,213)*y(k,93)
         mat(k,229) = mat(k,229) + rxt(k,591)*y(k,68)
         mat(k,649) = (rxt(k,392)+rxt(k,393))*y(k,68)
         mat(k,1374) = rxt(k,592)*y(k,54) + rxt(k,256)*y(k,91) + rxt(k,153)*y(k,92) &
                      + rxt(k,155)*y(k,93) + 1.150_r8*rxt(k,584)*y(k,95)
         mat(k,446) = mat(k,446) + rxt(k,272)*y(k,60)
         mat(k,806) = mat(k,806) + 2.000_r8*rxt(k,420)*y(k,131)
         mat(k,1667) = mat(k,1667) + rxt(k,359)*y(k,61) + rxt(k,363)*y(k,69) &
                      + rxt(k,350)*y(k,86)
         mat(k,272) = rxt(k,598)*y(k,68) + rxt(k,597)*y(k,106)
         mat(k,133) = rxt(k,602)*y(k,68)
         mat(k,1252) = -(rxt(k,126)*y(k,91) + (rxt(k,133) + rxt(k,135)) * y(k,95) &
                      + rxt(k,322)*y(k,125) + rxt(k,362)*y(k,133) + rxt(k,364) &
                      *y(k,126) + rxt(k,392)*y(k,122) + rxt(k,397)*y(k,123) + rxt(k,405) &
                      *y(k,67) + rxt(k,411)*y(k,38) + rxt(k,425)*y(k,106) + rxt(k,430) &
                      *y(k,54) + rxt(k,544)*y(k,74) + rxt(k,550)*y(k,76) + rxt(k,586) &
                      *y(k,107) + (rxt(k,590) + rxt(k,591)) * y(k,121) + rxt(k,598) &
                      *y(k,134) + rxt(k,602)*y(k,135))
         mat(k,832) = -rxt(k,126)*y(k,68)
         mat(k,1554) = -(rxt(k,133) + rxt(k,135)) * y(k,68)
         mat(k,1326) = -rxt(k,322)*y(k,68)
         mat(k,1663) = -rxt(k,362)*y(k,68)
         mat(k,1370) = -rxt(k,364)*y(k,68)
         mat(k,645) = -rxt(k,392)*y(k,68)
         mat(k,249) = -rxt(k,397)*y(k,68)
         mat(k,1429) = -rxt(k,405)*y(k,68)
         mat(k,1784) = -rxt(k,411)*y(k,68)
         mat(k,220) = -rxt(k,425)*y(k,68)
         mat(k,317) = -rxt(k,430)*y(k,68)
         mat(k,172) = -rxt(k,544)*y(k,68)
         mat(k,386) = -rxt(k,550)*y(k,68)
         mat(k,158) = -rxt(k,586)*y(k,68)
         mat(k,227) = -(rxt(k,590) + rxt(k,591)) * y(k,68)
         mat(k,270) = -rxt(k,598)*y(k,68)
         mat(k,131) = -rxt(k,602)*y(k,68)
         mat(k,593) = 2.000_r8*rxt(k,489)*y(k,4) + (rxt(k,491)+rxt(k,492))*y(k,26) &
                      + rxt(k,497)*y(k,67) + rxt(k,493)*y(k,98)
         mat(k,329) = rxt(k,533)*y(k,98)
         mat(k,1749) = (rxt(k,491)+rxt(k,492))*y(k,4) + (2.000_r8*rxt(k,458) &
                       +2.000_r8*rxt(k,459))*y(k,26) + rxt(k,467)*y(k,67) + rxt(k,116) &
                      *y(k,91) + rxt(k,128)*y(k,94) + rxt(k,461)*y(k,98) + rxt(k,315) &
                      *y(k,125) + rxt(k,469)*y(k,131) + rxt(k,351)*y(k,133)
         mat(k,1784) = mat(k,1784) + rxt(k,414)*y(k,69) + rxt(k,408)*y(k,98) &
                      + rxt(k,331)*y(k,128)
         mat(k,106) = rxt(k,423)*y(k,131)
         mat(k,317) = mat(k,317) + rxt(k,429)*y(k,61)
         mat(k,1512) = rxt(k,440)*y(k,69) + rxt(k,594)*y(k,126) + rxt(k,279)*y(k,128)
         mat(k,996) = rxt(k,429)*y(k,54) + rxt(k,431)*y(k,67) + rxt(k,432)*y(k,69) &
                      + rxt(k,320)*y(k,125) + rxt(k,277)*y(k,128)
         mat(k,1118) = rxt(k,436)*y(k,67) + rxt(k,434)*y(k,98)
         mat(k,1429) = mat(k,1429) + rxt(k,497)*y(k,4) + rxt(k,467)*y(k,26) &
                      + rxt(k,431)*y(k,61) + rxt(k,436)*y(k,62) + 2.000_r8*rxt(k,404) &
                      *y(k,67) + 2.000_r8*rxt(k,403)*y(k,69) + rxt(k,113)*y(k,90) &
                      + rxt(k,131)*y(k,94) + rxt(k,412)*y(k,98) + rxt(k,301)*y(k,110) &
                      + rxt(k,396)*y(k,123) + rxt(k,325)*y(k,125) + ( &
                      + 2.000_r8*rxt(k,333)+rxt(k,334))*y(k,128) + rxt(k,336)*y(k,129) &
                      + rxt(k,418)*y(k,131) + rxt(k,360)*y(k,133)
         mat(k,1252) = mat(k,1252) + 2.000_r8*rxt(k,397)*y(k,123)
         mat(k,1203) = rxt(k,414)*y(k,38) + rxt(k,440)*y(k,60) + rxt(k,432)*y(k,61) &
                      + 2.000_r8*rxt(k,403)*y(k,67) + rxt(k,545)*y(k,74) + rxt(k,551) &
                      *y(k,76) + rxt(k,488)*y(k,85) + rxt(k,456)*y(k,86) + rxt(k,132) &
                      *y(k,94) + rxt(k,134)*y(k,95) + 2.000_r8*rxt(k,413)*y(k,98) &
                      + rxt(k,292)*y(k,108) + 2.000_r8*rxt(k,302)*y(k,110) &
                      + 2.000_r8*rxt(k,394)*y(k,122) + rxt(k,323)*y(k,125) &
                      + 3.000_r8*rxt(k,332)*y(k,128) + rxt(k,419)*y(k,131)
         mat(k,172) = mat(k,172) + rxt(k,545)*y(k,69)
         mat(k,386) = mat(k,386) + rxt(k,551)*y(k,69)
         mat(k,461) = rxt(k,488)*y(k,69) + rxt(k,487)*y(k,98)
         mat(k,1471) = rxt(k,456)*y(k,69) + rxt(k,127)*y(k,94) + rxt(k,454)*y(k,98) &
                      + rxt(k,314)*y(k,125)
         mat(k,699) = rxt(k,154)*y(k,126)
         mat(k,416) = rxt(k,159)*y(k,126)
         mat(k,400) = rxt(k,257)*y(k,126)
         mat(k,291) = rxt(k,113)*y(k,67)
         mat(k,832) = mat(k,832) + rxt(k,116)*y(k,26) + rxt(k,256)*y(k,126)
         mat(k,368) = rxt(k,153)*y(k,126)
         mat(k,549) = rxt(k,155)*y(k,126)
         mat(k,674) = rxt(k,128)*y(k,26) + rxt(k,131)*y(k,67) + rxt(k,132)*y(k,69) &
                      + rxt(k,127)*y(k,86) + rxt(k,191)*y(k,100) + rxt(k,225)*y(k,102) &
                      + rxt(k,174)*y(k,104) + rxt(k,244)*y(k,115) + rxt(k,227) &
                      *y(k,116) + rxt(k,209)*y(k,119) + 2.000_r8*rxt(k,151)*y(k,126)
         mat(k,1554) = mat(k,1554) + rxt(k,134)*y(k,69) + rxt(k,326)*y(k,127) &
                      + 2.000_r8*rxt(k,380)*y(k,130)
         mat(k,341) = rxt(k,149)*y(k,126)
         mat(k,755) = rxt(k,493)*y(k,4) + rxt(k,533)*y(k,20) + rxt(k,461)*y(k,26) &
                      + rxt(k,408)*y(k,38) + rxt(k,434)*y(k,62) + rxt(k,412)*y(k,67) &
                      + 2.000_r8*rxt(k,413)*y(k,69) + rxt(k,487)*y(k,85) + rxt(k,454) &
                      *y(k,86) + 2.000_r8*rxt(k,422)*y(k,98) + rxt(k,417)*y(k,131)
         mat(k,1038) = rxt(k,191)*y(k,94) + rxt(k,190)*y(k,125)
         mat(k,1705) = rxt(k,225)*y(k,94) + rxt(k,214)*y(k,125)
         mat(k,1596) = rxt(k,174)*y(k,94) + rxt(k,173)*y(k,125)
         mat(k,861) = rxt(k,292)*y(k,69) + rxt(k,161)*y(k,126)
         mat(k,431) = rxt(k,156)*y(k,126)
         mat(k,723) = rxt(k,301)*y(k,67) + 2.000_r8*rxt(k,302)*y(k,69) + rxt(k,148) &
                      *y(k,126)
         mat(k,528) = rxt(k,160)*y(k,126)
         mat(k,568) = rxt(k,152)*y(k,126)
         mat(k,473) = rxt(k,157)*y(k,126)
         mat(k,510) = rxt(k,255)*y(k,126)
         mat(k,908) = rxt(k,244)*y(k,94) + rxt(k,243)*y(k,125)
         mat(k,950) = rxt(k,227)*y(k,94) + rxt(k,226)*y(k,125)
         mat(k,1079) = rxt(k,209)*y(k,94) + rxt(k,208)*y(k,125)
         mat(k,645) = mat(k,645) + 2.000_r8*rxt(k,394)*y(k,69)
         mat(k,249) = mat(k,249) + rxt(k,396)*y(k,67) + 2.000_r8*rxt(k,397)*y(k,68) &
                      + 2.000_r8*rxt(k,321)*y(k,125) + 2.000_r8*rxt(k,339)*y(k,130)
         mat(k,1326) = mat(k,1326) + rxt(k,315)*y(k,26) + rxt(k,320)*y(k,61) &
                      + rxt(k,325)*y(k,67) + rxt(k,323)*y(k,69) + rxt(k,314)*y(k,86) &
                      + rxt(k,190)*y(k,100) + rxt(k,214)*y(k,102) + rxt(k,173) &
                      *y(k,104) + rxt(k,243)*y(k,115) + rxt(k,226)*y(k,116) &
                      + rxt(k,208)*y(k,119) + 2.000_r8*rxt(k,321)*y(k,123)
         mat(k,1370) = mat(k,1370) + rxt(k,594)*y(k,60) + rxt(k,154)*y(k,87) &
                      + rxt(k,159)*y(k,88) + rxt(k,257)*y(k,89) + rxt(k,256)*y(k,91) &
                      + rxt(k,153)*y(k,92) + rxt(k,155)*y(k,93) + 2.000_r8*rxt(k,151) &
                      *y(k,94) + rxt(k,149)*y(k,97) + rxt(k,161)*y(k,108) + rxt(k,156) &
                      *y(k,109) + rxt(k,148)*y(k,110) + rxt(k,160)*y(k,111) &
                      + rxt(k,152)*y(k,112) + rxt(k,157)*y(k,113) + rxt(k,255) &
                      *y(k,114)
         mat(k,192) = rxt(k,326)*y(k,95) + (rxt(k,327)+rxt(k,328))*y(k,137)
         mat(k,444) = rxt(k,331)*y(k,38) + rxt(k,279)*y(k,60) + rxt(k,277)*y(k,61) + ( &
                      + 2.000_r8*rxt(k,333)+rxt(k,334))*y(k,67) + 3.000_r8*rxt(k,332) &
                      *y(k,69)
         mat(k,98) = rxt(k,336)*y(k,67)
         mat(k,354) = 2.000_r8*rxt(k,380)*y(k,95) + 2.000_r8*rxt(k,339)*y(k,123) &
                      + rxt(k,337)*y(k,137)
         mat(k,802) = rxt(k,469)*y(k,26) + rxt(k,423)*y(k,50) + rxt(k,418)*y(k,67) &
                      + rxt(k,419)*y(k,69) + rxt(k,417)*y(k,98)
         mat(k,1663) = mat(k,1663) + rxt(k,351)*y(k,26) + rxt(k,360)*y(k,67)
         mat(k,1842) = (rxt(k,327)+rxt(k,328))*y(k,127) + rxt(k,337)*y(k,130)
      end do
      end subroutine nlnmat03
      subroutine nlnmat04( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,1202) = -(rxt(k,132)*y(k,94) + rxt(k,134)*y(k,95) + rxt(k,292)*y(k,108) &
                      + rxt(k,302)*y(k,110) + rxt(k,323)*y(k,125) + rxt(k,332) &
                      *y(k,128) + rxt(k,348)*y(k,132) + rxt(k,363)*y(k,133) + rxt(k,394) &
                      *y(k,122) + rxt(k,403)*y(k,67) + rxt(k,413)*y(k,98) + rxt(k,414) &
                      *y(k,38) + rxt(k,419)*y(k,131) + rxt(k,432)*y(k,61) + rxt(k,440) &
                      *y(k,60) + rxt(k,456)*y(k,86) + rxt(k,488)*y(k,85) + rxt(k,545) &
                      *y(k,74) + rxt(k,551)*y(k,76))
         mat(k,673) = -rxt(k,132)*y(k,69)
         mat(k,1553) = -rxt(k,134)*y(k,69)
         mat(k,860) = -rxt(k,292)*y(k,69)
         mat(k,722) = -rxt(k,302)*y(k,69)
         mat(k,1325) = -rxt(k,323)*y(k,69)
         mat(k,443) = -rxt(k,332)*y(k,69)
         mat(k,1627) = -rxt(k,348)*y(k,69)
         mat(k,1662) = -rxt(k,363)*y(k,69)
         mat(k,644) = -rxt(k,394)*y(k,69)
         mat(k,1428) = -rxt(k,403)*y(k,69)
         mat(k,754) = -rxt(k,413)*y(k,69)
         mat(k,1783) = -rxt(k,414)*y(k,69)
         mat(k,801) = -rxt(k,419)*y(k,69)
         mat(k,995) = -rxt(k,432)*y(k,69)
         mat(k,1511) = -rxt(k,440)*y(k,69)
         mat(k,1470) = -rxt(k,456)*y(k,69)
         mat(k,460) = -rxt(k,488)*y(k,69)
         mat(k,171) = -rxt(k,545)*y(k,69)
         mat(k,385) = -rxt(k,551)*y(k,69)
         mat(k,995) = mat(k,995) + rxt(k,278)*y(k,128)
         mat(k,1428) = mat(k,1428) + rxt(k,405)*y(k,68) + rxt(k,324)*y(k,125) &
                      + rxt(k,338)*y(k,130)
         mat(k,1251) = rxt(k,405)*y(k,67)
         mat(k,248) = rxt(k,361)*y(k,133)
         mat(k,1325) = mat(k,1325) + rxt(k,324)*y(k,67)
         mat(k,443) = mat(k,443) + rxt(k,278)*y(k,61)
         mat(k,353) = rxt(k,338)*y(k,67)
         mat(k,1662) = mat(k,1662) + rxt(k,361)*y(k,123)
         mat(k,64) = -(rxt(k,552)*y(k,76))
         mat(k,376) = -rxt(k,552)*y(k,70)
         mat(k,578) = rxt(k,490)*y(k,26)
         mat(k,1723) = rxt(k,490)*y(k,4) + 2.000_r8*rxt(k,460)*y(k,26)
         mat(k,69) = -(rxt(k,542)*y(k,67) + rxt(k,543)*y(k,131))
         mat(k,1386) = -rxt(k,542)*y(k,71)
         mat(k,772) = -rxt(k,543)*y(k,71)
         mat(k,168) = -(rxt(k,544)*y(k,68) + rxt(k,545)*y(k,69) + rxt(k,548)*y(k,131))
         mat(k,1223) = -rxt(k,544)*y(k,74)
         mat(k,1178) = -rxt(k,545)*y(k,74)
         mat(k,777) = -rxt(k,548)*y(k,74)
         mat(k,379) = -(rxt(k,546)*y(k,4) + rxt(k,547)*y(k,26) + rxt(k,549)*y(k,61) &
                      + rxt(k,550)*y(k,68) + rxt(k,551)*y(k,69) + rxt(k,552)*y(k,70) &
                      + rxt(k,553)*y(k,131))
         mat(k,582) = -rxt(k,546)*y(k,76)
         mat(k,1729) = -rxt(k,547)*y(k,76)
         mat(k,972) = -rxt(k,549)*y(k,76)
         mat(k,1231) = -rxt(k,550)*y(k,76)
         mat(k,1182) = -rxt(k,551)*y(k,76)
         mat(k,66) = -rxt(k,552)*y(k,76)
         mat(k,787) = -rxt(k,553)*y(k,76)
         mat(k,1407) = rxt(k,542)*y(k,71)
         mat(k,1231) = mat(k,1231) + rxt(k,544)*y(k,74)
         mat(k,1182) = mat(k,1182) + rxt(k,545)*y(k,74)
         mat(k,73) = rxt(k,542)*y(k,67)
         mat(k,169) = rxt(k,544)*y(k,68) + rxt(k,545)*y(k,69) + rxt(k,548)*y(k,131)
         mat(k,787) = mat(k,787) + rxt(k,548)*y(k,74)
         mat(k,255) = -(rxt(k,554)*y(k,131))
         mat(k,783) = -rxt(k,554)*y(k,77)
         mat(k,581) = rxt(k,546)*y(k,76)
         mat(k,1725) = rxt(k,547)*y(k,76)
         mat(k,59) = rxt(k,556)*y(k,62) + (rxt(k,557)+.500_r8*rxt(k,558))*y(k,131)
         mat(k,968) = rxt(k,549)*y(k,76)
         mat(k,1096) = rxt(k,556)*y(k,33)
         mat(k,1227) = rxt(k,550)*y(k,76)
         mat(k,1180) = rxt(k,551)*y(k,76)
         mat(k,65) = rxt(k,552)*y(k,76)
         mat(k,72) = rxt(k,543)*y(k,131)
         mat(k,378) = rxt(k,546)*y(k,4) + rxt(k,547)*y(k,26) + rxt(k,549)*y(k,61) &
                      + rxt(k,550)*y(k,68) + rxt(k,551)*y(k,69) + rxt(k,552)*y(k,70) &
                      + rxt(k,553)*y(k,131)
         mat(k,783) = mat(k,783) + (rxt(k,557)+.500_r8*rxt(k,558))*y(k,33) &
                      + rxt(k,543)*y(k,71) + rxt(k,553)*y(k,76)
         mat(k,46) = -(rxt(k,555)*y(k,137))
         mat(k,1801) = -rxt(k,555)*y(k,78)
         mat(k,254) = rxt(k,554)*y(k,131)
         mat(k,770) = rxt(k,554)*y(k,77)
         mat(k,454) = -(rxt(k,486)*y(k,16) + rxt(k,487)*y(k,98) + rxt(k,488)*y(k,69))
         mat(k,603) = -rxt(k,486)*y(k,85)
         mat(k,742) = -rxt(k,487)*y(k,85)
         mat(k,1184) = -rxt(k,488)*y(k,85)
         mat(k,583) = 4.000_r8*rxt(k,489)*y(k,4) + (rxt(k,490)+rxt(k,491))*y(k,26) &
                      + rxt(k,494)*y(k,60) + rxt(k,497)*y(k,67) + rxt(k,546)*y(k,76) &
                      + rxt(k,498)*y(k,131)
         mat(k,1730) = (rxt(k,490)+rxt(k,491))*y(k,4)
         mat(k,239) = rxt(k,499)*y(k,67) + rxt(k,505)*y(k,122) + rxt(k,500)*y(k,131)
         mat(k,1492) = rxt(k,494)*y(k,4)
         mat(k,1409) = rxt(k,497)*y(k,4) + rxt(k,499)*y(k,43)
         mat(k,380) = rxt(k,546)*y(k,4)
         mat(k,633) = rxt(k,505)*y(k,43)
         mat(k,788) = rxt(k,498)*y(k,4) + rxt(k,500)*y(k,43)
         mat(k,1476) = -((rxt(k,114) + rxt(k,115)) * y(k,91) + rxt(k,127)*y(k,94) &
                      + rxt(k,285)*y(k,108) + rxt(k,314)*y(k,125) + rxt(k,341) &
                      *y(k,132) + rxt(k,350)*y(k,133) + rxt(k,450)*y(k,16) + rxt(k,452) &
                      *y(k,39) + rxt(k,453)*y(k,41) + (rxt(k,454) + rxt(k,455) &
                      ) * y(k,98) + rxt(k,456)*y(k,69) + rxt(k,463)*y(k,27) + rxt(k,472) &
                      *y(k,52))
         mat(k,837) = -(rxt(k,114) + rxt(k,115)) * y(k,86)
         mat(k,678) = -rxt(k,127)*y(k,86)
         mat(k,866) = -rxt(k,285)*y(k,86)
         mat(k,1331) = -rxt(k,314)*y(k,86)
         mat(k,1633) = -rxt(k,341)*y(k,86)
         mat(k,1668) = -rxt(k,350)*y(k,86)
         mat(k,617) = -rxt(k,450)*y(k,86)
         mat(k,495) = -rxt(k,452)*y(k,86)
         mat(k,151) = -rxt(k,453)*y(k,86)
         mat(k,760) = -(rxt(k,454) + rxt(k,455)) * y(k,86)
         mat(k,1208) = -rxt(k,456)*y(k,86)
         mat(k,284) = -rxt(k,463)*y(k,86)
         mat(k,235) = -rxt(k,472)*y(k,86)
         mat(k,596) = rxt(k,491)*y(k,26)
         mat(k,330) = rxt(k,457)*y(k,26)
         mat(k,1754) = rxt(k,491)*y(k,4) + rxt(k,457)*y(k,20) + (4.000_r8*rxt(k,458) &
                       +2.000_r8*rxt(k,460))*y(k,26) + rxt(k,462)*y(k,60) + rxt(k,467) &
                      *y(k,67) + rxt(k,547)*y(k,76) + rxt(k,468)*y(k,131)
         mat(k,35) = rxt(k,512)*y(k,122)
         mat(k,1167) = rxt(k,470)*y(k,67) + rxt(k,482)*y(k,122) + rxt(k,471)*y(k,131)
         mat(k,1517) = rxt(k,462)*y(k,26) + rxt(k,111)*y(k,90)
         mat(k,1001) = rxt(k,110)*y(k,87)
         mat(k,1434) = rxt(k,467)*y(k,26) + rxt(k,470)*y(k,47)
         mat(k,388) = rxt(k,547)*y(k,26)
         mat(k,702) = rxt(k,110)*y(k,61) + rxt(k,195)*y(k,100) + rxt(k,147)*y(k,102) &
                      + rxt(k,177)*y(k,104) + rxt(k,248)*y(k,115) + rxt(k,230) &
                      *y(k,116) + rxt(k,212)*y(k,119) + rxt(k,154)*y(k,126)
         mat(k,418) = rxt(k,199)*y(k,100) + rxt(k,164)*y(k,102) + rxt(k,182)*y(k,104) &
                      + rxt(k,252)*y(k,115) + rxt(k,234)*y(k,116) + rxt(k,217) &
                      *y(k,119) + rxt(k,159)*y(k,126)
         mat(k,402) = rxt(k,187)*y(k,100) + rxt(k,181)*y(k,102) + rxt(k,169)*y(k,104) &
                      + rxt(k,240)*y(k,115) + rxt(k,222)*y(k,116) + rxt(k,205) &
                      *y(k,119) + rxt(k,257)*y(k,126)
         mat(k,293) = rxt(k,111)*y(k,60)
         mat(k,1043) = rxt(k,195)*y(k,87) + rxt(k,199)*y(k,88) + rxt(k,187)*y(k,89)
         mat(k,1710) = rxt(k,147)*y(k,87) + rxt(k,164)*y(k,88) + rxt(k,181)*y(k,89)
         mat(k,1601) = rxt(k,177)*y(k,87) + rxt(k,182)*y(k,88) + rxt(k,169)*y(k,89)
         mat(k,913) = rxt(k,248)*y(k,87) + rxt(k,252)*y(k,88) + rxt(k,240)*y(k,89)
         mat(k,955) = rxt(k,230)*y(k,87) + rxt(k,234)*y(k,88) + rxt(k,222)*y(k,89)
         mat(k,1084) = rxt(k,212)*y(k,87) + rxt(k,217)*y(k,88) + rxt(k,205)*y(k,89)
         mat(k,650) = rxt(k,512)*y(k,32) + rxt(k,482)*y(k,47)
         mat(k,1375) = rxt(k,154)*y(k,87) + rxt(k,159)*y(k,88) + rxt(k,257)*y(k,89)
         mat(k,807) = rxt(k,468)*y(k,26) + rxt(k,471)*y(k,47)
         mat(k,690) = -(rxt(k,105)*y(k,38) + rxt(k,107)*y(k,137) + rxt(k,108)*y(k,47) &
                      + rxt(k,109)*y(k,49) + rxt(k,110)*y(k,61) + rxt(k,147)*y(k,102) &
                      + rxt(k,154)*y(k,126) + rxt(k,177)*y(k,104) + rxt(k,195) &
                      *y(k,100) + rxt(k,212)*y(k,119) + rxt(k,230)*y(k,116) + rxt(k,248) &
                      *y(k,115))
         mat(k,1770) = -rxt(k,105)*y(k,87)
         mat(k,1828) = -rxt(k,107)*y(k,87)
         mat(k,1148) = -rxt(k,108)*y(k,87)
         mat(k,1277) = -rxt(k,109)*y(k,87)
         mat(k,982) = -rxt(k,110)*y(k,87)
         mat(k,1691) = -rxt(k,147)*y(k,87)
         mat(k,1356) = -rxt(k,154)*y(k,87)
         mat(k,1582) = -rxt(k,177)*y(k,87)
         mat(k,1024) = -rxt(k,195)*y(k,87)
         mat(k,1065) = -rxt(k,212)*y(k,87)
         mat(k,936) = -rxt(k,230)*y(k,87)
         mat(k,894) = -rxt(k,248)*y(k,87)
         mat(k,1735) = rxt(k,116)*y(k,91) + rxt(k,286)*y(k,108) + rxt(k,351)*y(k,133)
         mat(k,1148) = mat(k,1148) + rxt(k,130)*y(k,94) + rxt(k,289)*y(k,108) &
                      + rxt(k,298)*y(k,110) + rxt(k,318)*y(k,125) + rxt(k,345) &
                      *y(k,132) + rxt(k,356)*y(k,133)
         mat(k,1498) = rxt(k,112)*y(k,90)
         mat(k,1415) = rxt(k,113)*y(k,90)
         mat(k,1457) = rxt(k,114)*y(k,91) + rxt(k,127)*y(k,94) + rxt(k,285)*y(k,108) &
                      + rxt(k,314)*y(k,125) + rxt(k,341)*y(k,132) + rxt(k,350) &
                      *y(k,133)
         mat(k,288) = rxt(k,112)*y(k,60) + rxt(k,113)*y(k,67)
         mat(k,819) = rxt(k,116)*y(k,26) + rxt(k,114)*y(k,86)
         mat(k,661) = rxt(k,130)*y(k,47) + rxt(k,127)*y(k,86)
         mat(k,849) = rxt(k,286)*y(k,26) + rxt(k,289)*y(k,47) + rxt(k,285)*y(k,86)
         mat(k,712) = rxt(k,298)*y(k,47)
         mat(k,1312) = rxt(k,318)*y(k,47) + rxt(k,314)*y(k,86)
         mat(k,1614) = rxt(k,345)*y(k,47) + rxt(k,341)*y(k,86)
         mat(k,1649) = rxt(k,351)*y(k,26) + rxt(k,356)*y(k,47) + rxt(k,350)*y(k,86)
         mat(k,409) = -(rxt(k,106)*y(k,47) + rxt(k,159)*y(k,126) + rxt(k,164)*y(k,102) &
                      + rxt(k,182)*y(k,104) + rxt(k,199)*y(k,100) + rxt(k,217) &
                      *y(k,119) + rxt(k,234)*y(k,116) + rxt(k,252)*y(k,115))
         mat(k,1140) = -rxt(k,106)*y(k,88)
         mat(k,1347) = -rxt(k,159)*y(k,88)
         mat(k,1683) = -rxt(k,164)*y(k,88)
         mat(k,1574) = -rxt(k,182)*y(k,88)
         mat(k,1016) = -rxt(k,199)*y(k,88)
         mat(k,1057) = -rxt(k,217)*y(k,88)
         mat(k,928) = -rxt(k,234)*y(k,88)
         mat(k,885) = -rxt(k,252)*y(k,88)
         mat(k,689) = rxt(k,107)*y(k,137)
         mat(k,1818) = rxt(k,107)*y(k,87)
         mat(k,393) = -((rxt(k,169) + rxt(k,268)) * y(k,104) + (rxt(k,181) + rxt(k,267) &
                      ) * y(k,102) + rxt(k,187)*y(k,100) + rxt(k,205)*y(k,119) &
                      + rxt(k,222)*y(k,116) + rxt(k,240)*y(k,115) + rxt(k,257) &
                      *y(k,126))
         mat(k,1573) = -(rxt(k,169) + rxt(k,268)) * y(k,89)
         mat(k,1682) = -(rxt(k,181) + rxt(k,267)) * y(k,89)
         mat(k,1015) = -rxt(k,187)*y(k,89)
         mat(k,1056) = -rxt(k,205)*y(k,89)
         mat(k,927) = -rxt(k,222)*y(k,89)
         mat(k,884) = -rxt(k,240)*y(k,89)
         mat(k,1346) = -rxt(k,257)*y(k,89)
         mat(k,1139) = rxt(k,108)*y(k,87) + rxt(k,106)*y(k,88)
         mat(k,688) = rxt(k,108)*y(k,47)
         mat(k,408) = rxt(k,106)*y(k,47)
      end do
      end subroutine nlnmat04
      subroutine nlnmat05( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,287) = -((rxt(k,111) + rxt(k,112)) * y(k,60) + rxt(k,113)*y(k,67))
         mat(k,1487) = -(rxt(k,111) + rxt(k,112)) * y(k,90)
         mat(k,1404) = -rxt(k,113)*y(k,90)
         mat(k,1727) = rxt(k,128)*y(k,94) + rxt(k,315)*y(k,125) + rxt(k,342)*y(k,132)
         mat(k,1449) = rxt(k,115)*y(k,91)
         mat(k,815) = rxt(k,115)*y(k,86)
         mat(k,658) = rxt(k,128)*y(k,26)
         mat(k,1308) = rxt(k,315)*y(k,26)
         mat(k,1611) = rxt(k,342)*y(k,26)
         mat(k,822) = -((rxt(k,114) + rxt(k,115)) * y(k,86) + rxt(k,116)*y(k,26) &
                      + rxt(k,117)*y(k,38) + rxt(k,119)*y(k,137) + rxt(k,124)*y(k,49) &
                      + rxt(k,125)*y(k,67) + rxt(k,126)*y(k,68) + (rxt(k,168) &
                      + rxt(k,265)) * y(k,104) + (rxt(k,170) + rxt(k,258)) * y(k,102) &
                      + rxt(k,186)*y(k,100) + rxt(k,204)*y(k,119) + rxt(k,221) &
                      *y(k,116) + rxt(k,239)*y(k,115) + rxt(k,256)*y(k,126) + rxt(k,280) &
                      *y(k,60) + rxt(k,281)*y(k,61))
         mat(k,1461) = -(rxt(k,114) + rxt(k,115)) * y(k,91)
         mat(k,1739) = -rxt(k,116)*y(k,91)
         mat(k,1774) = -rxt(k,117)*y(k,91)
         mat(k,1832) = -rxt(k,119)*y(k,91)
         mat(k,1281) = -rxt(k,124)*y(k,91)
         mat(k,1419) = -rxt(k,125)*y(k,91)
         mat(k,1242) = -rxt(k,126)*y(k,91)
         mat(k,1586) = -(rxt(k,168) + rxt(k,265)) * y(k,91)
         mat(k,1695) = -(rxt(k,170) + rxt(k,258)) * y(k,91)
         mat(k,1028) = -rxt(k,186)*y(k,91)
         mat(k,1069) = -rxt(k,204)*y(k,91)
         mat(k,940) = -rxt(k,221)*y(k,91)
         mat(k,898) = -rxt(k,239)*y(k,91)
         mat(k,1360) = -rxt(k,256)*y(k,91)
         mat(k,1502) = -rxt(k,280)*y(k,91)
         mat(k,986) = -rxt(k,281)*y(k,91)
         mat(k,1774) = mat(k,1774) + rxt(k,129)*y(k,94)
         mat(k,1419) = mat(k,1419) + rxt(k,131)*y(k,94)
         mat(k,665) = rxt(k,129)*y(k,38) + rxt(k,131)*y(k,67)
         mat(k,362) = -(rxt(k,153)*y(k,126) + (rxt(k,176) + rxt(k,259)) * y(k,104) &
                      + rxt(k,194)*y(k,100) + rxt(k,211)*y(k,119) + rxt(k,229) &
                      *y(k,116) + rxt(k,246)*y(k,115) + (rxt(k,247) + rxt(k,271) &
                      ) * y(k,102))
         mat(k,1345) = -rxt(k,153)*y(k,92)
         mat(k,1572) = -(rxt(k,176) + rxt(k,259)) * y(k,92)
         mat(k,1014) = -rxt(k,194)*y(k,92)
         mat(k,1055) = -rxt(k,211)*y(k,92)
         mat(k,926) = -rxt(k,229)*y(k,92)
         mat(k,883) = -rxt(k,246)*y(k,92)
         mat(k,1681) = -(rxt(k,247) + rxt(k,271)) * y(k,92)
         mat(k,537) = rxt(k,118)*y(k,137)
         mat(k,1816) = rxt(k,118)*y(k,93)
         mat(k,539) = -(rxt(k,118)*y(k,137) + (rxt(k,120) + rxt(k,121)) * y(k,61) &
                      + (rxt(k,122) + rxt(k,123)) * y(k,60) + rxt(k,155)*y(k,126) &
                      + (rxt(k,158) + rxt(k,260)) * y(k,102) + (rxt(k,178) + rxt(k,261) &
                      ) * y(k,104) + rxt(k,196)*y(k,100) + rxt(k,213)*y(k,119) &
                      + rxt(k,231)*y(k,116) + rxt(k,249)*y(k,115))
         mat(k,1823) = -rxt(k,118)*y(k,93)
         mat(k,977) = -(rxt(k,120) + rxt(k,121)) * y(k,93)
         mat(k,1493) = -(rxt(k,122) + rxt(k,123)) * y(k,93)
         mat(k,1352) = -rxt(k,155)*y(k,93)
         mat(k,1688) = -(rxt(k,158) + rxt(k,260)) * y(k,93)
         mat(k,1579) = -(rxt(k,178) + rxt(k,261)) * y(k,93)
         mat(k,1021) = -rxt(k,196)*y(k,93)
         mat(k,1062) = -rxt(k,213)*y(k,93)
         mat(k,933) = -rxt(k,231)*y(k,93)
         mat(k,890) = -rxt(k,249)*y(k,93)
         mat(k,817) = rxt(k,119)*y(k,137)
         mat(k,1823) = mat(k,1823) + rxt(k,119)*y(k,91)
         mat(k,660) = -(rxt(k,127)*y(k,86) + rxt(k,128)*y(k,26) + rxt(k,129)*y(k,38) &
                      + rxt(k,130)*y(k,47) + rxt(k,131)*y(k,67) + rxt(k,132)*y(k,69) &
                      + rxt(k,151)*y(k,126) + rxt(k,174)*y(k,104) + rxt(k,191) &
                      *y(k,100) + rxt(k,209)*y(k,119) + rxt(k,225)*y(k,102) + rxt(k,227) &
                      *y(k,116) + rxt(k,244)*y(k,115))
         mat(k,1456) = -rxt(k,127)*y(k,94)
         mat(k,1734) = -rxt(k,128)*y(k,94)
         mat(k,1769) = -rxt(k,129)*y(k,94)
         mat(k,1147) = -rxt(k,130)*y(k,94)
         mat(k,1414) = -rxt(k,131)*y(k,94)
         mat(k,1188) = -rxt(k,132)*y(k,94)
         mat(k,1355) = -rxt(k,151)*y(k,94)
         mat(k,1581) = -rxt(k,174)*y(k,94)
         mat(k,1023) = -rxt(k,191)*y(k,94)
         mat(k,1064) = -rxt(k,209)*y(k,94)
         mat(k,1690) = -rxt(k,225)*y(k,94)
         mat(k,935) = -rxt(k,227)*y(k,94)
         mat(k,893) = -rxt(k,244)*y(k,94)
         mat(k,1561) = -((rxt(k,133) + rxt(k,135)) * y(k,68) + rxt(k,134)*y(k,69) &
                      + rxt(k,138)*y(k,96) + rxt(k,141)*y(k,102) + rxt(k,144)*y(k,104) &
                      + rxt(k,303)*y(k,116) + rxt(k,304)*y(k,117) + rxt(k,306) &
                      *y(k,118) + rxt(k,308)*y(k,119) + rxt(k,326)*y(k,127) + rxt(k,380) &
                      *y(k,130) + rxt(k,381)*y(k,105) + rxt(k,382)*y(k,99) + rxt(k,383) &
                      *y(k,100) + rxt(k,384)*y(k,120) + rxt(k,583)*y(k,115) + rxt(k,584) &
                      *y(k,126) + rxt(k,585)*y(k,107))
         mat(k,1259) = -(rxt(k,133) + rxt(k,135)) * y(k,95)
         mat(k,1210) = -rxt(k,134)*y(k,95)
         mat(k,86) = -rxt(k,138)*y(k,95)
         mat(k,1712) = -rxt(k,141)*y(k,95)
         mat(k,1603) = -rxt(k,144)*y(k,95)
         mat(k,957) = -rxt(k,303)*y(k,95)
         mat(k,93) = -rxt(k,304)*y(k,95)
         mat(k,119) = -rxt(k,306)*y(k,95)
         mat(k,1086) = -rxt(k,308)*y(k,95)
         mat(k,194) = -rxt(k,326)*y(k,95)
         mat(k,358) = -rxt(k,380)*y(k,95)
         mat(k,183) = -rxt(k,381)*y(k,95)
         mat(k,165) = -rxt(k,382)*y(k,95)
         mat(k,1045) = -rxt(k,383)*y(k,95)
         mat(k,113) = -rxt(k,384)*y(k,95)
         mat(k,915) = -rxt(k,583)*y(k,95)
         mat(k,1377) = -rxt(k,584)*y(k,95)
         mat(k,161) = -rxt(k,585)*y(k,95)
         mat(k,1791) = rxt(k,105)*y(k,87) + rxt(k,317)*y(k,125) + rxt(k,344)*y(k,132)
         mat(k,497) = rxt(k,353)*y(k,133)
         mat(k,1519) = rxt(k,136)*y(k,133)
         mat(k,1436) = rxt(k,324)*y(k,125) + rxt(k,333)*y(k,128) + rxt(k,347)*y(k,132) &
                      + rxt(k,360)*y(k,133)
         mat(k,1210) = mat(k,1210) + rxt(k,332)*y(k,128)
         mat(k,704) = rxt(k,105)*y(k,38)
         mat(k,252) = rxt(k,321)*y(k,125) + rxt(k,361)*y(k,133)
         mat(k,1333) = rxt(k,317)*y(k,38) + rxt(k,324)*y(k,67) + rxt(k,321)*y(k,123)
         mat(k,448) = rxt(k,333)*y(k,67) + rxt(k,332)*y(k,69)
         mat(k,1635) = rxt(k,344)*y(k,38) + rxt(k,347)*y(k,67)
         mat(k,1670) = rxt(k,353)*y(k,39) + rxt(k,136)*y(k,60) + rxt(k,360)*y(k,67) &
                      + rxt(k,361)*y(k,123)
         mat(k,83) = -(rxt(k,138)*y(k,95) + rxt(k,139)*y(k,137))
         mat(k,1527) = -rxt(k,138)*y(k,96)
         mat(k,1804) = -rxt(k,139)*y(k,96)
         mat(k,186) = rxt(k,327)*y(k,137)
         mat(k,1804) = mat(k,1804) + rxt(k,327)*y(k,127)
         mat(k,335) = -(rxt(k,149)*y(k,126) + rxt(k,172)*y(k,104) + rxt(k,189) &
                      *y(k,100) + rxt(k,203)*y(k,102) + rxt(k,207)*y(k,119) + rxt(k,224) &
                      *y(k,116) + rxt(k,242)*y(k,115))
         mat(k,1343) = -rxt(k,149)*y(k,97)
         mat(k,1571) = -rxt(k,172)*y(k,97)
         mat(k,1013) = -rxt(k,189)*y(k,97)
         mat(k,1680) = -rxt(k,203)*y(k,97)
         mat(k,1054) = -rxt(k,207)*y(k,97)
         mat(k,925) = -rxt(k,224)*y(k,97)
         mat(k,882) = -rxt(k,242)*y(k,97)
         mat(k,747) = -(rxt(k,310)*y(k,119) + (rxt(k,408) + rxt(k,409) + rxt(k,410) &
                      ) * y(k,38) + rxt(k,412)*y(k,67) + rxt(k,413)*y(k,69) + rxt(k,417) &
                      *y(k,131) + 4._r8*rxt(k,422)*y(k,98) + rxt(k,434)*y(k,62) &
                      + rxt(k,439)*y(k,60) + rxt(k,444)*y(k,61) + (rxt(k,454) &
                      + rxt(k,455)) * y(k,86) + rxt(k,461)*y(k,26) + rxt(k,487) &
                      *y(k,85) + rxt(k,493)*y(k,4) + rxt(k,533)*y(k,20))
         mat(k,1067) = -rxt(k,310)*y(k,98)
         mat(k,1772) = -(rxt(k,408) + rxt(k,409) + rxt(k,410)) * y(k,98)
         mat(k,1417) = -rxt(k,412)*y(k,98)
         mat(k,1191) = -rxt(k,413)*y(k,98)
         mat(k,793) = -rxt(k,417)*y(k,98)
         mat(k,1107) = -rxt(k,434)*y(k,98)
         mat(k,1500) = -rxt(k,439)*y(k,98)
         mat(k,984) = -rxt(k,444)*y(k,98)
         mat(k,1459) = -(rxt(k,454) + rxt(k,455)) * y(k,98)
         mat(k,1737) = -rxt(k,461)*y(k,98)
         mat(k,458) = -rxt(k,487)*y(k,98)
         mat(k,587) = -rxt(k,493)*y(k,98)
         mat(k,326) = -rxt(k,533)*y(k,98)
         mat(k,587) = mat(k,587) + rxt(k,498)*y(k,131)
         mat(k,608) = rxt(k,530)*y(k,62) + rxt(k,531)*y(k,67) + rxt(k,486)*y(k,85) &
                      + rxt(k,450)*y(k,86)
         mat(k,326) = mat(k,326) + rxt(k,457)*y(k,26) + rxt(k,534)*y(k,60)
         mat(k,1737) = mat(k,1737) + rxt(k,457)*y(k,20) + rxt(k,468)*y(k,131)
         mat(k,135) = rxt(k,537)*y(k,131)
         mat(k,60) = .500_r8*rxt(k,558)*y(k,131)
         mat(k,1772) = mat(k,1772) + rxt(k,411)*y(k,68) + rxt(k,317)*y(k,125)
         mat(k,147) = rxt(k,407)*y(k,67) + rxt(k,453)*y(k,86) + rxt(k,416)*y(k,131)
         mat(k,1150) = rxt(k,130)*y(k,94) + rxt(k,318)*y(k,125)
         mat(k,1279) = rxt(k,319)*y(k,125)
         mat(k,1500) = mat(k,1500) + rxt(k,534)*y(k,20)
         mat(k,1107) = mat(k,1107) + rxt(k,530)*y(k,16) + rxt(k,437)*y(k,131)
         mat(k,1417) = mat(k,1417) + rxt(k,531)*y(k,16) + rxt(k,407)*y(k,41) &
                      + rxt(k,347)*y(k,132)
         mat(k,1240) = rxt(k,411)*y(k,38)
         mat(k,1191) = mat(k,1191) + rxt(k,419)*y(k,131)
         mat(k,257) = rxt(k,554)*y(k,131)
         mat(k,458) = mat(k,458) + rxt(k,486)*y(k,16)
         mat(k,1459) = mat(k,1459) + rxt(k,450)*y(k,16) + rxt(k,453)*y(k,41)
         mat(k,663) = rxt(k,130)*y(k,47)
         mat(k,1314) = rxt(k,317)*y(k,38) + rxt(k,318)*y(k,47) + rxt(k,319)*y(k,49)
         mat(k,793) = mat(k,793) + rxt(k,498)*y(k,4) + rxt(k,468)*y(k,26) + rxt(k,537) &
                      *y(k,29) + .500_r8*rxt(k,558)*y(k,33) + rxt(k,416)*y(k,41) &
                      + rxt(k,437)*y(k,62) + rxt(k,419)*y(k,69) + rxt(k,554)*y(k,77)
         mat(k,1616) = rxt(k,347)*y(k,67)
         mat(k,162) = -(rxt(k,374)*y(k,137) + rxt(k,382)*y(k,95))
         mat(k,1808) = -rxt(k,374)*y(k,99)
         mat(k,1534) = -rxt(k,382)*y(k,99)
         mat(k,84) = rxt(k,139)*y(k,137)
         mat(k,180) = rxt(k,372)*y(k,137)
         mat(k,1808) = mat(k,1808) + rxt(k,139)*y(k,96) + rxt(k,372)*y(k,105)
         mat(k,1033) = -(rxt(k,185)*y(k,114) + rxt(k,186)*y(k,91) + rxt(k,187)*y(k,89) &
                      + rxt(k,188)*y(k,110) + rxt(k,189)*y(k,97) + rxt(k,190)*y(k,125) &
                      + rxt(k,191)*y(k,94) + rxt(k,193)*y(k,112) + rxt(k,194)*y(k,92) &
                      + rxt(k,195)*y(k,87) + rxt(k,196)*y(k,93) + rxt(k,197)*y(k,109) &
                      + rxt(k,198)*y(k,113) + rxt(k,199)*y(k,88) + rxt(k,200)*y(k,111) &
                      + rxt(k,201)*y(k,108) + rxt(k,376)*y(k,137) + rxt(k,383)*y(k,95))
         mat(k,507) = -rxt(k,185)*y(k,100)
         mat(k,827) = -rxt(k,186)*y(k,100)
         mat(k,397) = -rxt(k,187)*y(k,100)
         mat(k,718) = -rxt(k,188)*y(k,100)
         mat(k,339) = -rxt(k,189)*y(k,100)
         mat(k,1321) = -rxt(k,190)*y(k,100)
         mat(k,670) = -rxt(k,191)*y(k,100)
         mat(k,565) = -rxt(k,193)*y(k,100)
         mat(k,366) = -rxt(k,194)*y(k,100)
         mat(k,696) = -rxt(k,195)*y(k,100)
         mat(k,547) = -rxt(k,196)*y(k,100)
         mat(k,429) = -rxt(k,197)*y(k,100)
         mat(k,469) = -rxt(k,198)*y(k,100)
         mat(k,413) = -rxt(k,199)*y(k,100)
         mat(k,525) = -rxt(k,200)*y(k,100)
         mat(k,856) = -rxt(k,201)*y(k,100)
         mat(k,1837) = -rxt(k,376)*y(k,100)
         mat(k,1549) = -rxt(k,383)*y(k,100)
         mat(k,164) = rxt(k,374)*y(k,137)
         mat(k,91) = rxt(k,305)*y(k,137)
         mat(k,1837) = mat(k,1837) + rxt(k,374)*y(k,99) + rxt(k,305)*y(k,117)
      end do
      end subroutine nlnmat05
      subroutine nlnmat06( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,50) = -(rxt(k,140)*y(k,137))
         mat(k,1802) = -rxt(k,140)*y(k,101)
         mat(k,295) = rxt(k,142)*y(k,102)
         mat(k,1678) = rxt(k,142)*y(k,56)
         mat(k,1716) = -(rxt(k,141)*y(k,95) + rxt(k,142)*y(k,56) + (rxt(k,146) &
                      + rxt(k,269)) * y(k,114) + rxt(k,147)*y(k,87) + (rxt(k,158) &
                      + rxt(k,260)) * y(k,93) + rxt(k,162)*y(k,109) + rxt(k,163) &
                      *y(k,113) + rxt(k,164)*y(k,88) + rxt(k,165)*y(k,111) + rxt(k,166) &
                      *y(k,108) + (rxt(k,170) + rxt(k,258)) * y(k,91) + (rxt(k,181) &
                      + rxt(k,267)) * y(k,89) + (rxt(k,192) + rxt(k,264)) * y(k,110) &
                      + rxt(k,203)*y(k,97) + rxt(k,214)*y(k,125) + rxt(k,225)*y(k,94) &
                      + (rxt(k,236) + rxt(k,262)) * y(k,112) + (rxt(k,247) + rxt(k,271) &
                      ) * y(k,92) + rxt(k,378)*y(k,137))
         mat(k,1565) = -rxt(k,141)*y(k,102)
         mat(k,307) = -rxt(k,142)*y(k,102)
         mat(k,515) = -(rxt(k,146) + rxt(k,269)) * y(k,102)
         mat(k,706) = -rxt(k,147)*y(k,102)
         mat(k,554) = -(rxt(k,158) + rxt(k,260)) * y(k,102)
         mat(k,435) = -rxt(k,162)*y(k,102)
         mat(k,478) = -rxt(k,163)*y(k,102)
         mat(k,421) = -rxt(k,164)*y(k,102)
         mat(k,534) = -rxt(k,165)*y(k,102)
         mat(k,872) = -rxt(k,166)*y(k,102)
         mat(k,843) = -(rxt(k,170) + rxt(k,258)) * y(k,102)
         mat(k,405) = -(rxt(k,181) + rxt(k,267)) * y(k,102)
         mat(k,731) = -(rxt(k,192) + rxt(k,264)) * y(k,102)
         mat(k,345) = -rxt(k,203)*y(k,102)
         mat(k,1337) = -rxt(k,214)*y(k,102)
         mat(k,684) = -rxt(k,225)*y(k,102)
         mat(k,574) = -(rxt(k,236) + rxt(k,262)) * y(k,102)
         mat(k,373) = -(rxt(k,247) + rxt(k,271)) * y(k,102)
         mat(k,1853) = -rxt(k,378)*y(k,102)
         mat(k,1049) = rxt(k,376)*y(k,137)
         mat(k,52) = rxt(k,140)*y(k,137)
         mat(k,1853) = mat(k,1853) + rxt(k,376)*y(k,100) + rxt(k,140)*y(k,101)
         mat(k,54) = -(rxt(k,143)*y(k,137))
         mat(k,1803) = -rxt(k,143)*y(k,103)
         mat(k,296) = rxt(k,145)*y(k,104)
         mat(k,1569) = rxt(k,145)*y(k,56)
         mat(k,1604) = -(rxt(k,144)*y(k,95) + rxt(k,145)*y(k,56) + (rxt(k,167) &
                      + rxt(k,270)) * y(k,114) + (rxt(k,168) + rxt(k,265)) * y(k,91) &
                      + (rxt(k,169) + rxt(k,268)) * y(k,89) + (rxt(k,171) + rxt(k,266) &
                      ) * y(k,110) + rxt(k,172)*y(k,97) + rxt(k,173)*y(k,125) &
                      + rxt(k,174)*y(k,94) + (rxt(k,175) + rxt(k,263)) * y(k,112) &
                      + (rxt(k,176) + rxt(k,259)) * y(k,92) + rxt(k,177)*y(k,87) &
                      + (rxt(k,178) + rxt(k,261)) * y(k,93) + rxt(k,179)*y(k,109) &
                      + rxt(k,180)*y(k,113) + rxt(k,182)*y(k,88) + rxt(k,183)*y(k,111) &
                      + rxt(k,184)*y(k,108))
         mat(k,1562) = -rxt(k,144)*y(k,104)
         mat(k,306) = -rxt(k,145)*y(k,104)
         mat(k,514) = -(rxt(k,167) + rxt(k,270)) * y(k,104)
         mat(k,840) = -(rxt(k,168) + rxt(k,265)) * y(k,104)
         mat(k,404) = -(rxt(k,169) + rxt(k,268)) * y(k,104)
         mat(k,730) = -(rxt(k,171) + rxt(k,266)) * y(k,104)
         mat(k,344) = -rxt(k,172)*y(k,104)
         mat(k,1334) = -rxt(k,173)*y(k,104)
         mat(k,681) = -rxt(k,174)*y(k,104)
         mat(k,573) = -(rxt(k,175) + rxt(k,263)) * y(k,104)
         mat(k,372) = -(rxt(k,176) + rxt(k,259)) * y(k,104)
         mat(k,705) = -rxt(k,177)*y(k,104)
         mat(k,553) = -(rxt(k,178) + rxt(k,261)) * y(k,104)
         mat(k,434) = -rxt(k,179)*y(k,104)
         mat(k,477) = -rxt(k,180)*y(k,104)
         mat(k,420) = -rxt(k,182)*y(k,104)
         mat(k,533) = -rxt(k,183)*y(k,104)
         mat(k,869) = -rxt(k,184)*y(k,104)
         mat(k,1713) = rxt(k,378)*y(k,137)
         mat(k,56) = rxt(k,143)*y(k,137)
         mat(k,1850) = rxt(k,378)*y(k,102) + rxt(k,143)*y(k,103)
         mat(k,181) = -(rxt(k,372)*y(k,137) + rxt(k,381)*y(k,95))
         mat(k,1810) = -rxt(k,372)*y(k,105)
         mat(k,1536) = -rxt(k,381)*y(k,105)
         mat(k,1764) = rxt(k,309)*y(k,119)
         mat(k,737) = rxt(k,310)*y(k,119)
         mat(k,1053) = rxt(k,309)*y(k,38) + rxt(k,310)*y(k,98) + rxt(k,311)*y(k,131)
         mat(k,188) = rxt(k,328)*y(k,137)
         mat(k,779) = rxt(k,311)*y(k,119)
         mat(k,1810) = mat(k,1810) + rxt(k,328)*y(k,127)
         mat(k,215) = -(rxt(k,424)*y(k,67) + rxt(k,425)*y(k,68) + rxt(k,597)*y(k,134))
         mat(k,1396) = -rxt(k,424)*y(k,106)
         mat(k,1224) = -rxt(k,425)*y(k,106)
         mat(k,262) = -rxt(k,597)*y(k,106)
         mat(k,1396) = mat(k,1396) + rxt(k,587)*y(k,107)
         mat(k,1538) = .900_r8*rxt(k,585)*y(k,107) + .800_r8*rxt(k,583)*y(k,115)
         mat(k,154) = rxt(k,587)*y(k,67) + .900_r8*rxt(k,585)*y(k,95)
         mat(k,878) = .800_r8*rxt(k,583)*y(k,95)
         mat(k,153) = -(rxt(k,585)*y(k,95) + rxt(k,586)*y(k,68) + (rxt(k,587) &
                      + rxt(k,588)) * y(k,67))
         mat(k,1533) = -rxt(k,585)*y(k,107)
         mat(k,1222) = -rxt(k,586)*y(k,107)
         mat(k,1393) = -(rxt(k,587) + rxt(k,588)) * y(k,107)
         mat(k,852) = -(rxt(k,161)*y(k,126) + rxt(k,166)*y(k,102) + rxt(k,184) &
                      *y(k,104) + rxt(k,201)*y(k,100) + rxt(k,219)*y(k,119) + rxt(k,237) &
                      *y(k,116) + rxt(k,254)*y(k,115) + rxt(k,285)*y(k,86) + rxt(k,286) &
                      *y(k,26) + rxt(k,287)*y(k,38) + rxt(k,288)*y(k,137) + rxt(k,289) &
                      *y(k,47) + rxt(k,290)*y(k,49) + rxt(k,291)*y(k,61) + rxt(k,292) &
                      *y(k,69))
         mat(k,1361) = -rxt(k,161)*y(k,108)
         mat(k,1696) = -rxt(k,166)*y(k,108)
         mat(k,1587) = -rxt(k,184)*y(k,108)
         mat(k,1029) = -rxt(k,201)*y(k,108)
         mat(k,1070) = -rxt(k,219)*y(k,108)
         mat(k,941) = -rxt(k,237)*y(k,108)
         mat(k,899) = -rxt(k,254)*y(k,108)
         mat(k,1462) = -rxt(k,285)*y(k,108)
         mat(k,1740) = -rxt(k,286)*y(k,108)
         mat(k,1775) = -rxt(k,287)*y(k,108)
         mat(k,1833) = -rxt(k,288)*y(k,108)
         mat(k,1153) = -rxt(k,289)*y(k,108)
         mat(k,1282) = -rxt(k,290)*y(k,108)
         mat(k,987) = -rxt(k,291)*y(k,108)
         mat(k,1194) = -rxt(k,292)*y(k,108)
         mat(k,1503) = rxt(k,111)*y(k,90) + rxt(k,280)*y(k,91) + rxt(k,123)*y(k,93) &
                      + rxt(k,279)*y(k,128)
         mat(k,987) = mat(k,987) + rxt(k,110)*y(k,87) + rxt(k,320)*y(k,125) &
                      + rxt(k,278)*y(k,128) + rxt(k,346)*y(k,132) + rxt(k,359) &
                      *y(k,133)
         mat(k,1420) = rxt(k,301)*y(k,110)
         mat(k,1194) = mat(k,1194) + rxt(k,302)*y(k,110)
         mat(k,692) = rxt(k,110)*y(k,61)
         mat(k,289) = rxt(k,111)*y(k,60)
         mat(k,823) = rxt(k,280)*y(k,60)
         mat(k,543) = rxt(k,123)*y(k,60)
         mat(k,714) = rxt(k,301)*y(k,67) + rxt(k,302)*y(k,69)
         mat(k,1317) = rxt(k,320)*y(k,61)
         mat(k,441) = rxt(k,279)*y(k,60) + rxt(k,278)*y(k,61)
         mat(k,1619) = rxt(k,346)*y(k,61)
         mat(k,1654) = rxt(k,359)*y(k,61)
         mat(k,424) = -(rxt(k,156)*y(k,126) + rxt(k,162)*y(k,102) + rxt(k,179) &
                      *y(k,104) + rxt(k,197)*y(k,100) + rxt(k,215)*y(k,119) + rxt(k,232) &
                      *y(k,116) + rxt(k,250)*y(k,115))
         mat(k,1348) = -rxt(k,156)*y(k,109)
         mat(k,1684) = -rxt(k,162)*y(k,109)
         mat(k,1575) = -rxt(k,179)*y(k,109)
         mat(k,1017) = -rxt(k,197)*y(k,109)
         mat(k,1058) = -rxt(k,215)*y(k,109)
         mat(k,929) = -rxt(k,232)*y(k,109)
         mat(k,886) = -rxt(k,250)*y(k,109)
         mat(k,1490) = rxt(k,122)*y(k,93)
         mat(k,538) = rxt(k,122)*y(k,60)
         mat(k,848) = rxt(k,288)*y(k,137)
         mat(k,1819) = rxt(k,288)*y(k,108)
         mat(k,713) = -(rxt(k,148)*y(k,126) + (rxt(k,171) + rxt(k,266)) * y(k,104) &
                      + rxt(k,188)*y(k,100) + (rxt(k,192) + rxt(k,264)) * y(k,102) &
                      + rxt(k,206)*y(k,119) + rxt(k,223)*y(k,116) + rxt(k,241) &
                      *y(k,115) + (rxt(k,276) + rxt(k,298)) * y(k,47) + rxt(k,296) &
                      *y(k,137) + rxt(k,300)*y(k,49) + rxt(k,301)*y(k,67) + rxt(k,302) &
                      *y(k,69))
         mat(k,1357) = -rxt(k,148)*y(k,110)
         mat(k,1583) = -(rxt(k,171) + rxt(k,266)) * y(k,110)
         mat(k,1025) = -rxt(k,188)*y(k,110)
         mat(k,1692) = -(rxt(k,192) + rxt(k,264)) * y(k,110)
         mat(k,1066) = -rxt(k,206)*y(k,110)
         mat(k,937) = -rxt(k,223)*y(k,110)
         mat(k,895) = -rxt(k,241)*y(k,110)
         mat(k,1149) = -(rxt(k,276) + rxt(k,298)) * y(k,110)
         mat(k,1829) = -rxt(k,296)*y(k,110)
         mat(k,1278) = -rxt(k,300)*y(k,110)
         mat(k,1416) = -rxt(k,301)*y(k,110)
         mat(k,1190) = -rxt(k,302)*y(k,110)
         mat(k,1278) = mat(k,1278) + rxt(k,109)*y(k,87) + rxt(k,124)*y(k,91) &
                      + rxt(k,290)*y(k,108) + rxt(k,319)*y(k,125) + rxt(k,357) &
                      *y(k,133)
         mat(k,1499) = rxt(k,272)*y(k,128)
         mat(k,983) = rxt(k,281)*y(k,91) + rxt(k,120)*y(k,93) + rxt(k,291)*y(k,108) &
                      + rxt(k,277)*y(k,128)
         mat(k,1190) = mat(k,1190) + rxt(k,292)*y(k,108)
         mat(k,691) = rxt(k,109)*y(k,49)
         mat(k,820) = rxt(k,124)*y(k,49) + rxt(k,281)*y(k,61)
         mat(k,541) = rxt(k,120)*y(k,61)
         mat(k,850) = rxt(k,290)*y(k,49) + rxt(k,291)*y(k,61) + rxt(k,292)*y(k,69)
         mat(k,1313) = rxt(k,319)*y(k,49)
         mat(k,439) = rxt(k,272)*y(k,60) + rxt(k,277)*y(k,61)
         mat(k,1650) = rxt(k,357)*y(k,49)
         mat(k,519) = -(rxt(k,160)*y(k,126) + rxt(k,165)*y(k,102) + rxt(k,183) &
                      *y(k,104) + rxt(k,200)*y(k,100) + rxt(k,218)*y(k,119) + rxt(k,235) &
                      *y(k,116) + rxt(k,253)*y(k,115) + rxt(k,293)*y(k,56))
         mat(k,1351) = -rxt(k,160)*y(k,111)
         mat(k,1687) = -rxt(k,165)*y(k,111)
         mat(k,1578) = -rxt(k,183)*y(k,111)
         mat(k,1020) = -rxt(k,200)*y(k,111)
         mat(k,1061) = -rxt(k,218)*y(k,111)
         mat(k,932) = -rxt(k,235)*y(k,111)
         mat(k,889) = -rxt(k,253)*y(k,111)
         mat(k,299) = -rxt(k,293)*y(k,111)
         mat(k,559) = rxt(k,294)*y(k,137)
         mat(k,1822) = rxt(k,294)*y(k,112)
         mat(k,560) = -(rxt(k,152)*y(k,126) + (rxt(k,175) + rxt(k,263)) * y(k,104) &
                      + rxt(k,193)*y(k,100) + rxt(k,210)*y(k,119) + rxt(k,228) &
                      *y(k,116) + (rxt(k,236) + rxt(k,262)) * y(k,102) + rxt(k,245) &
                      *y(k,115) + rxt(k,294)*y(k,137) + rxt(k,295)*y(k,49) + rxt(k,297) &
                      *y(k,56))
         mat(k,1353) = -rxt(k,152)*y(k,112)
         mat(k,1580) = -(rxt(k,175) + rxt(k,263)) * y(k,112)
         mat(k,1022) = -rxt(k,193)*y(k,112)
         mat(k,1063) = -rxt(k,210)*y(k,112)
         mat(k,934) = -rxt(k,228)*y(k,112)
         mat(k,1689) = -(rxt(k,236) + rxt(k,262)) * y(k,112)
         mat(k,891) = -rxt(k,245)*y(k,112)
         mat(k,1824) = -rxt(k,294)*y(k,112)
         mat(k,1274) = -rxt(k,295)*y(k,112)
         mat(k,300) = -rxt(k,297)*y(k,112)
         mat(k,978) = rxt(k,121)*y(k,93)
         mat(k,540) = rxt(k,121)*y(k,61)
         mat(k,711) = rxt(k,296)*y(k,137)
         mat(k,1824) = mat(k,1824) + rxt(k,296)*y(k,110)
      end do
      end subroutine nlnmat06
      subroutine nlnmat07( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,465) = -(rxt(k,157)*y(k,126) + rxt(k,163)*y(k,102) + rxt(k,180) &
                      *y(k,104) + rxt(k,198)*y(k,100) + rxt(k,216)*y(k,119) + rxt(k,233) &
                      *y(k,116) + rxt(k,251)*y(k,115) + rxt(k,299)*y(k,49))
         mat(k,1349) = -rxt(k,157)*y(k,113)
         mat(k,1685) = -rxt(k,163)*y(k,113)
         mat(k,1576) = -rxt(k,180)*y(k,113)
         mat(k,1018) = -rxt(k,198)*y(k,113)
         mat(k,1059) = -rxt(k,216)*y(k,113)
         mat(k,930) = -rxt(k,233)*y(k,113)
         mat(k,887) = -rxt(k,251)*y(k,113)
         mat(k,1271) = -rxt(k,299)*y(k,113)
         mat(k,1142) = rxt(k,276)*y(k,110)
         mat(k,709) = rxt(k,276)*y(k,47)
         mat(k,503) = -((rxt(k,146) + rxt(k,269)) * y(k,102) + (rxt(k,167) + rxt(k,270) &
                      ) * y(k,104) + rxt(k,185)*y(k,100) + rxt(k,202)*y(k,119) &
                      + rxt(k,220)*y(k,116) + rxt(k,238)*y(k,115) + rxt(k,255) &
                      *y(k,126))
         mat(k,1686) = -(rxt(k,146) + rxt(k,269)) * y(k,114)
         mat(k,1577) = -(rxt(k,167) + rxt(k,270)) * y(k,114)
         mat(k,1019) = -rxt(k,185)*y(k,114)
         mat(k,1060) = -rxt(k,202)*y(k,114)
         mat(k,931) = -rxt(k,220)*y(k,114)
         mat(k,888) = -rxt(k,238)*y(k,114)
         mat(k,1350) = -rxt(k,255)*y(k,114)
         mat(k,1273) = rxt(k,300)*y(k,110) + rxt(k,295)*y(k,112) + rxt(k,299)*y(k,113)
         mat(k,298) = rxt(k,293)*y(k,111) + rxt(k,297)*y(k,112)
         mat(k,710) = rxt(k,300)*y(k,49)
         mat(k,518) = rxt(k,293)*y(k,56)
         mat(k,558) = rxt(k,295)*y(k,49) + rxt(k,297)*y(k,56)
         mat(k,466) = rxt(k,299)*y(k,49)
         mat(k,900) = -(rxt(k,238)*y(k,114) + rxt(k,239)*y(k,91) + rxt(k,240)*y(k,89) &
                      + rxt(k,241)*y(k,110) + rxt(k,242)*y(k,97) + rxt(k,243)*y(k,125) &
                      + rxt(k,244)*y(k,94) + rxt(k,245)*y(k,112) + rxt(k,246)*y(k,92) &
                      + rxt(k,248)*y(k,87) + rxt(k,249)*y(k,93) + rxt(k,250)*y(k,109) &
                      + rxt(k,251)*y(k,113) + rxt(k,252)*y(k,88) + rxt(k,253)*y(k,111) &
                      + rxt(k,254)*y(k,108) + rxt(k,365)*y(k,137) + rxt(k,583)*y(k,95))
         mat(k,505) = -rxt(k,238)*y(k,115)
         mat(k,824) = -rxt(k,239)*y(k,115)
         mat(k,395) = -rxt(k,240)*y(k,115)
         mat(k,715) = -rxt(k,241)*y(k,115)
         mat(k,337) = -rxt(k,242)*y(k,115)
         mat(k,1318) = -rxt(k,243)*y(k,115)
         mat(k,667) = -rxt(k,244)*y(k,115)
         mat(k,562) = -rxt(k,245)*y(k,115)
         mat(k,364) = -rxt(k,246)*y(k,115)
         mat(k,693) = -rxt(k,248)*y(k,115)
         mat(k,544) = -rxt(k,249)*y(k,115)
         mat(k,426) = -rxt(k,250)*y(k,115)
         mat(k,467) = -rxt(k,251)*y(k,115)
         mat(k,411) = -rxt(k,252)*y(k,115)
         mat(k,522) = -rxt(k,253)*y(k,115)
         mat(k,853) = -rxt(k,254)*y(k,115)
         mat(k,1834) = -rxt(k,365)*y(k,115)
         mat(k,1546) = -rxt(k,583)*y(k,115)
         mat(k,315) = rxt(k,592)*y(k,126)
         mat(k,1504) = rxt(k,594)*y(k,126)
         mat(k,1421) = rxt(k,587)*y(k,107)
         mat(k,1244) = rxt(k,591)*y(k,121)
         mat(k,157) = rxt(k,587)*y(k,67)
         mat(k,226) = rxt(k,591)*y(k,68)
         mat(k,1362) = rxt(k,592)*y(k,54) + rxt(k,594)*y(k,60)
         mat(k,943) = -(rxt(k,220)*y(k,114) + rxt(k,221)*y(k,91) + rxt(k,222)*y(k,89) &
                      + rxt(k,223)*y(k,110) + rxt(k,224)*y(k,97) + rxt(k,226)*y(k,125) &
                      + rxt(k,227)*y(k,94) + rxt(k,228)*y(k,112) + rxt(k,229)*y(k,92) &
                      + rxt(k,230)*y(k,87) + rxt(k,231)*y(k,93) + rxt(k,232)*y(k,109) &
                      + rxt(k,233)*y(k,113) + rxt(k,234)*y(k,88) + rxt(k,235)*y(k,111) &
                      + rxt(k,237)*y(k,108) + rxt(k,303)*y(k,95) + rxt(k,367)*y(k,137))
         mat(k,506) = -rxt(k,220)*y(k,116)
         mat(k,825) = -rxt(k,221)*y(k,116)
         mat(k,396) = -rxt(k,222)*y(k,116)
         mat(k,716) = -rxt(k,223)*y(k,116)
         mat(k,338) = -rxt(k,224)*y(k,116)
         mat(k,1319) = -rxt(k,226)*y(k,116)
         mat(k,668) = -rxt(k,227)*y(k,116)
         mat(k,563) = -rxt(k,228)*y(k,116)
         mat(k,365) = -rxt(k,229)*y(k,116)
         mat(k,694) = -rxt(k,230)*y(k,116)
         mat(k,545) = -rxt(k,231)*y(k,116)
         mat(k,427) = -rxt(k,232)*y(k,116)
         mat(k,468) = -rxt(k,233)*y(k,116)
         mat(k,412) = -rxt(k,234)*y(k,116)
         mat(k,523) = -rxt(k,235)*y(k,116)
         mat(k,854) = -rxt(k,237)*y(k,116)
         mat(k,1547) = -rxt(k,303)*y(k,116)
         mat(k,1835) = -rxt(k,367)*y(k,116)
         mat(k,1072) = rxt(k,366)*y(k,137)
         mat(k,1835) = mat(k,1835) + rxt(k,366)*y(k,119)
         mat(k,89) = -(rxt(k,304)*y(k,95) + rxt(k,305)*y(k,137))
         mat(k,1528) = -rxt(k,304)*y(k,117)
         mat(k,1805) = -rxt(k,305)*y(k,117)
         mat(k,923) = rxt(k,367)*y(k,137)
         mat(k,1805) = mat(k,1805) + rxt(k,367)*y(k,116)
         mat(k,115) = -(rxt(k,306)*y(k,95) + rxt(k,307)*y(k,137))
         mat(k,1530) = -rxt(k,306)*y(k,118)
         mat(k,1807) = -rxt(k,307)*y(k,118)
         mat(k,1075) = -(rxt(k,202)*y(k,114) + rxt(k,204)*y(k,91) + rxt(k,205)*y(k,89) &
                      + rxt(k,206)*y(k,110) + rxt(k,207)*y(k,97) + rxt(k,208)*y(k,125) &
                      + rxt(k,209)*y(k,94) + rxt(k,210)*y(k,112) + rxt(k,211)*y(k,92) &
                      + rxt(k,212)*y(k,87) + rxt(k,213)*y(k,93) + rxt(k,215)*y(k,109) &
                      + rxt(k,216)*y(k,113) + rxt(k,217)*y(k,88) + rxt(k,218)*y(k,111) &
                      + rxt(k,219)*y(k,108) + rxt(k,308)*y(k,95) + rxt(k,309)*y(k,38) &
                      + rxt(k,310)*y(k,98) + rxt(k,311)*y(k,131) + rxt(k,366)*y(k,137))
         mat(k,508) = -rxt(k,202)*y(k,119)
         mat(k,828) = -rxt(k,204)*y(k,119)
         mat(k,398) = -rxt(k,205)*y(k,119)
         mat(k,719) = -rxt(k,206)*y(k,119)
         mat(k,340) = -rxt(k,207)*y(k,119)
         mat(k,1322) = -rxt(k,208)*y(k,119)
         mat(k,671) = -rxt(k,209)*y(k,119)
         mat(k,566) = -rxt(k,210)*y(k,119)
         mat(k,367) = -rxt(k,211)*y(k,119)
         mat(k,697) = -rxt(k,212)*y(k,119)
         mat(k,548) = -rxt(k,213)*y(k,119)
         mat(k,430) = -rxt(k,215)*y(k,119)
         mat(k,470) = -rxt(k,216)*y(k,119)
         mat(k,414) = -rxt(k,217)*y(k,119)
         mat(k,526) = -rxt(k,218)*y(k,119)
         mat(k,857) = -rxt(k,219)*y(k,119)
         mat(k,1550) = -rxt(k,308)*y(k,119)
         mat(k,1780) = -rxt(k,309)*y(k,119)
         mat(k,751) = -rxt(k,310)*y(k,119)
         mat(k,798) = -rxt(k,311)*y(k,119)
         mat(k,1838) = -rxt(k,366)*y(k,119)
         mat(k,904) = rxt(k,365)*y(k,137)
         mat(k,117) = rxt(k,307)*y(k,137)
         mat(k,111) = rxt(k,313)*y(k,137)
         mat(k,1838) = mat(k,1838) + rxt(k,365)*y(k,115) + rxt(k,307)*y(k,118) &
                      + rxt(k,313)*y(k,120)
         mat(k,108) = -(rxt(k,313)*y(k,137) + rxt(k,384)*y(k,95))
         mat(k,1806) = -rxt(k,313)*y(k,120)
         mat(k,1529) = -rxt(k,384)*y(k,120)
         mat(k,223) = -(rxt(k,589)*y(k,67) + (rxt(k,590) + rxt(k,591)) * y(k,68))
         mat(k,1397) = -rxt(k,589)*y(k,121)
         mat(k,1225) = -(rxt(k,590) + rxt(k,591)) * y(k,121)
         mat(k,216) = rxt(k,597)*y(k,134)
         mat(k,263) = rxt(k,597)*y(k,106)
         mat(k,637) = -(rxt(k,389)*y(k,39) + rxt(k,390)*y(k,137) + (rxt(k,392) &
                      + rxt(k,393)) * y(k,68) + rxt(k,394)*y(k,69) + (rxt(k,482) &
                      + rxt(k,483)) * y(k,47) + (rxt(k,505) + rxt(k,506)) * y(k,43) &
                      + rxt(k,511)*y(k,31) + rxt(k,512)*y(k,32))
         mat(k,486) = -rxt(k,389)*y(k,122)
         mat(k,1827) = -rxt(k,390)*y(k,122)
         mat(k,1236) = -(rxt(k,392) + rxt(k,393)) * y(k,122)
         mat(k,1187) = -rxt(k,394)*y(k,122)
         mat(k,1146) = -(rxt(k,482) + rxt(k,483)) * y(k,122)
         mat(k,241) = -(rxt(k,505) + rxt(k,506)) * y(k,122)
         mat(k,29) = -rxt(k,511)*y(k,122)
         mat(k,34) = -rxt(k,512)*y(k,122)
         mat(k,1236) = mat(k,1236) + rxt(k,425)*y(k,106)
         mat(k,1543) = .850_r8*rxt(k,584)*y(k,126)
         mat(k,219) = rxt(k,425)*y(k,68)
         mat(k,1354) = .850_r8*rxt(k,584)*y(k,95)
         mat(k,246) = -(rxt(k,321)*y(k,125) + rxt(k,339)*y(k,130) + rxt(k,361) &
                      *y(k,133) + rxt(k,396)*y(k,67) + rxt(k,397)*y(k,68))
         mat(k,1307) = -rxt(k,321)*y(k,123)
         mat(k,349) = -rxt(k,339)*y(k,123)
         mat(k,1643) = -rxt(k,361)*y(k,123)
         mat(k,1400) = -rxt(k,396)*y(k,123)
         mat(k,1226) = -rxt(k,397)*y(k,123)
         mat(k,1400) = mat(k,1400) + rxt(k,400)*y(k,124)
         mat(k,1226) = mat(k,1226) + rxt(k,401)*y(k,124)
         mat(k,1179) = rxt(k,402)*y(k,124)
         mat(k,43) = rxt(k,400)*y(k,67) + rxt(k,401)*y(k,68) + rxt(k,402)*y(k,69)
         mat(k,42) = -(rxt(k,400)*y(k,67) + rxt(k,401)*y(k,68) + rxt(k,402)*y(k,69))
         mat(k,1385) = -rxt(k,400)*y(k,124)
         mat(k,1218) = -rxt(k,401)*y(k,124)
         mat(k,1177) = -rxt(k,402)*y(k,124)
         mat(k,1218) = mat(k,1218) + rxt(k,392)*y(k,122)
         mat(k,627) = rxt(k,392)*y(k,68)
         mat(k,1328) = -(rxt(k,150)*y(k,126) + rxt(k,173)*y(k,104) + rxt(k,190) &
                      *y(k,100) + rxt(k,208)*y(k,119) + rxt(k,214)*y(k,102) + rxt(k,226) &
                      *y(k,116) + rxt(k,243)*y(k,115) + rxt(k,314)*y(k,86) + rxt(k,315) &
                      *y(k,26) + rxt(k,317)*y(k,38) + rxt(k,318)*y(k,47) + rxt(k,319) &
                      *y(k,49) + rxt(k,320)*y(k,61) + rxt(k,321)*y(k,123) + rxt(k,322) &
                      *y(k,68) + rxt(k,323)*y(k,69) + (rxt(k,324) + rxt(k,325) &
                      ) * y(k,67))
         mat(k,1372) = -rxt(k,150)*y(k,125)
         mat(k,1598) = -rxt(k,173)*y(k,125)
         mat(k,1040) = -rxt(k,190)*y(k,125)
         mat(k,1081) = -rxt(k,208)*y(k,125)
         mat(k,1707) = -rxt(k,214)*y(k,125)
         mat(k,952) = -rxt(k,226)*y(k,125)
         mat(k,910) = -rxt(k,243)*y(k,125)
         mat(k,1473) = -rxt(k,314)*y(k,125)
         mat(k,1751) = -rxt(k,315)*y(k,125)
         mat(k,1786) = -rxt(k,317)*y(k,125)
         mat(k,1164) = -rxt(k,318)*y(k,125)
         mat(k,1293) = -rxt(k,319)*y(k,125)
         mat(k,998) = -rxt(k,320)*y(k,125)
         mat(k,250) = -rxt(k,321)*y(k,125)
         mat(k,1254) = -rxt(k,322)*y(k,125)
         mat(k,1205) = -rxt(k,323)*y(k,125)
         mat(k,1431) = -(rxt(k,324) + rxt(k,325)) * y(k,125)
         mat(k,1431) = mat(k,1431) + rxt(k,125)*y(k,91) + rxt(k,334)*y(k,128)
         mat(k,1254) = mat(k,1254) + (rxt(k,133)+rxt(k,135))*y(k,95)
         mat(k,834) = rxt(k,125)*y(k,67)
         mat(k,1556) = (rxt(k,133)+rxt(k,135))*y(k,68)
         mat(k,445) = rxt(k,334)*y(k,67)
         mat(k,1373) = -(rxt(k,148)*y(k,110) + rxt(k,149)*y(k,97) + rxt(k,150) &
                      *y(k,125) + rxt(k,151)*y(k,94) + rxt(k,152)*y(k,112) + rxt(k,153) &
                      *y(k,92) + rxt(k,154)*y(k,87) + rxt(k,155)*y(k,93) + rxt(k,156) &
                      *y(k,109) + rxt(k,157)*y(k,113) + rxt(k,159)*y(k,88) + rxt(k,160) &
                      *y(k,111) + rxt(k,161)*y(k,108) + rxt(k,255)*y(k,114) + rxt(k,256) &
                      *y(k,91) + rxt(k,257)*y(k,89) + rxt(k,329)*y(k,137) + rxt(k,364) &
                      *y(k,68) + rxt(k,584)*y(k,95) + rxt(k,592)*y(k,54) + rxt(k,594) &
                      *y(k,60))
         mat(k,725) = -rxt(k,148)*y(k,126)
         mat(k,342) = -rxt(k,149)*y(k,126)
         mat(k,1329) = -rxt(k,150)*y(k,126)
         mat(k,676) = -rxt(k,151)*y(k,126)
         mat(k,570) = -rxt(k,152)*y(k,126)
         mat(k,369) = -rxt(k,153)*y(k,126)
         mat(k,701) = -rxt(k,154)*y(k,126)
         mat(k,550) = -rxt(k,155)*y(k,126)
         mat(k,432) = -rxt(k,156)*y(k,126)
         mat(k,475) = -rxt(k,157)*y(k,126)
         mat(k,417) = -rxt(k,159)*y(k,126)
         mat(k,530) = -rxt(k,160)*y(k,126)
         mat(k,864) = -rxt(k,161)*y(k,126)
         mat(k,512) = -rxt(k,255)*y(k,126)
         mat(k,835) = -rxt(k,256)*y(k,126)
         mat(k,401) = -rxt(k,257)*y(k,126)
         mat(k,1845) = -rxt(k,329)*y(k,126)
         mat(k,1255) = -rxt(k,364)*y(k,126)
         mat(k,1557) = -rxt(k,584)*y(k,126)
         mat(k,318) = -rxt(k,592)*y(k,126)
         mat(k,1515) = -rxt(k,594)*y(k,126)
         mat(k,1432) = rxt(k,338)*y(k,130)
         mat(k,1255) = mat(k,1255) + rxt(k,586)*y(k,107) + rxt(k,590)*y(k,121) &
                      + rxt(k,598)*y(k,134) + rxt(k,602)*y(k,135)
         mat(k,159) = rxt(k,586)*y(k,68)
         mat(k,228) = rxt(k,590)*y(k,68)
         mat(k,251) = rxt(k,339)*y(k,130)
         mat(k,1329) = mat(k,1329) + 2.000_r8*rxt(k,150)*y(k,126)
         mat(k,1373) = mat(k,1373) + 2.000_r8*rxt(k,150)*y(k,125)
         mat(k,356) = rxt(k,338)*y(k,67) + rxt(k,339)*y(k,123)
         mat(k,271) = rxt(k,598)*y(k,68)
         mat(k,132) = rxt(k,602)*y(k,68)
      end do
      end subroutine nlnmat07
      subroutine nlnmat08( avec_len, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k,189) = -(rxt(k,326)*y(k,95) + (rxt(k,327) + rxt(k,328)) * y(k,137))
         mat(k,1537) = -rxt(k,326)*y(k,127)
         mat(k,1811) = -(rxt(k,327) + rxt(k,328)) * y(k,127)
         mat(k,1341) = rxt(k,329)*y(k,137)
         mat(k,348) = rxt(k,337)*y(k,137)
         mat(k,1811) = mat(k,1811) + rxt(k,329)*y(k,126) + rxt(k,337)*y(k,130)
         mat(k,438) = -((rxt(k,272) + rxt(k,279)) * y(k,60) + (rxt(k,277) + rxt(k,278) &
                      ) * y(k,61) + rxt(k,331)*y(k,38) + rxt(k,332)*y(k,69) + (rxt(k,333) &
                      + rxt(k,334)) * y(k,67))
         mat(k,1491) = -(rxt(k,272) + rxt(k,279)) * y(k,128)
         mat(k,973) = -(rxt(k,277) + rxt(k,278)) * y(k,128)
         mat(k,1765) = -rxt(k,331)*y(k,128)
         mat(k,1183) = -rxt(k,332)*y(k,128)
         mat(k,1408) = -(rxt(k,333) + rxt(k,334)) * y(k,128)
         mat(k,1408) = mat(k,1408) + rxt(k,336)*y(k,129)
         mat(k,1232) = rxt(k,126)*y(k,91) + rxt(k,362)*y(k,133)
         mat(k,1183) = mat(k,1183) + rxt(k,132)*y(k,94) + rxt(k,323)*y(k,125) &
                      + rxt(k,348)*y(k,132) + rxt(k,363)*y(k,133)
         mat(k,816) = rxt(k,126)*y(k,68)
         mat(k,659) = rxt(k,132)*y(k,69)
         mat(k,1310) = rxt(k,323)*y(k,69)
         mat(k,96) = rxt(k,336)*y(k,67)
         mat(k,1613) = rxt(k,348)*y(k,69)
         mat(k,1645) = rxt(k,362)*y(k,68) + rxt(k,363)*y(k,69)
         mat(k,95) = -(rxt(k,336)*y(k,67))
         mat(k,1387) = -rxt(k,336)*y(k,129)
         mat(k,1219) = rxt(k,322)*y(k,125)
         mat(k,1306) = rxt(k,322)*y(k,68)
         mat(k,350) = -(rxt(k,337)*y(k,137) + rxt(k,338)*y(k,67) + rxt(k,339)*y(k,123) &
                      + rxt(k,380)*y(k,95))
         mat(k,1815) = -rxt(k,337)*y(k,130)
         mat(k,1406) = -rxt(k,338)*y(k,130)
         mat(k,247) = -rxt(k,339)*y(k,130)
         mat(k,1542) = -rxt(k,380)*y(k,130)
         mat(k,1230) = rxt(k,364)*y(k,126)
         mat(k,1344) = rxt(k,364)*y(k,68)
         mat(k,794) = -(rxt(k,311)*y(k,119) + rxt(k,340)*y(k,53) + rxt(k,349)*y(k,60) &
                      + rxt(k,415)*y(k,39) + rxt(k,416)*y(k,41) + rxt(k,417)*y(k,98) &
                      + rxt(k,418)*y(k,67) + rxt(k,419)*y(k,69) + (4._r8*rxt(k,420) &
                      + 4._r8*rxt(k,421)) * y(k,131) + rxt(k,423)*y(k,50) + rxt(k,437) &
                      *y(k,62) + rxt(k,438)*y(k,54) + rxt(k,446)*y(k,61) + rxt(k,447) &
                      *y(k,49) + rxt(k,466)*y(k,27) + (rxt(k,468) + rxt(k,469) &
                      ) * y(k,26) + rxt(k,471)*y(k,47) + rxt(k,474)*y(k,52) + rxt(k,498) &
                      *y(k,4) + rxt(k,500)*y(k,43) + rxt(k,532)*y(k,16) + rxt(k,535) &
                      *y(k,21) + (rxt(k,537) + rxt(k,541)) * y(k,29) + rxt(k,543) &
                      *y(k,71) + rxt(k,548)*y(k,74) + rxt(k,553)*y(k,76) + rxt(k,554) &
                      *y(k,77) + (rxt(k,557) + rxt(k,558)) * y(k,33))
         mat(k,1068) = -rxt(k,311)*y(k,131)
         mat(k,176) = -rxt(k,340)*y(k,131)
         mat(k,1501) = -rxt(k,349)*y(k,131)
         mat(k,488) = -rxt(k,415)*y(k,131)
         mat(k,148) = -rxt(k,416)*y(k,131)
         mat(k,748) = -rxt(k,417)*y(k,131)
         mat(k,1418) = -rxt(k,418)*y(k,131)
         mat(k,1192) = -rxt(k,419)*y(k,131)
         mat(k,103) = -rxt(k,423)*y(k,131)
         mat(k,1108) = -rxt(k,437)*y(k,131)
         mat(k,314) = -rxt(k,438)*y(k,131)
         mat(k,985) = -rxt(k,446)*y(k,131)
         mat(k,1280) = -rxt(k,447)*y(k,131)
         mat(k,278) = -rxt(k,466)*y(k,131)
         mat(k,1738) = -(rxt(k,468) + rxt(k,469)) * y(k,131)
         mat(k,1151) = -rxt(k,471)*y(k,131)
         mat(k,232) = -rxt(k,474)*y(k,131)
         mat(k,588) = -rxt(k,498)*y(k,131)
         mat(k,242) = -rxt(k,500)*y(k,131)
         mat(k,609) = -rxt(k,532)*y(k,131)
         mat(k,80) = -rxt(k,535)*y(k,131)
         mat(k,136) = -(rxt(k,537) + rxt(k,541)) * y(k,131)
         mat(k,74) = -rxt(k,543)*y(k,131)
         mat(k,170) = -rxt(k,548)*y(k,131)
         mat(k,383) = -rxt(k,553)*y(k,131)
         mat(k,258) = -rxt(k,554)*y(k,131)
         mat(k,61) = -(rxt(k,557) + rxt(k,558)) * y(k,131)
         mat(k,609) = mat(k,609) + rxt(k,531)*y(k,67)
         mat(k,80) = mat(k,80) + .300_r8*rxt(k,535)*y(k,131)
         mat(k,1738) = mat(k,1738) + rxt(k,342)*y(k,132)
         mat(k,201) = rxt(k,509)*y(k,137)
         mat(k,1773) = rxt(k,414)*y(k,69) + rxt(k,129)*y(k,94) + 2.000_r8*rxt(k,409) &
                      *y(k,98)
         mat(k,488) = mat(k,488) + rxt(k,406)*y(k,67) + rxt(k,389)*y(k,122)
         mat(k,148) = mat(k,148) + rxt(k,407)*y(k,67)
         mat(k,242) = mat(k,242) + rxt(k,499)*y(k,67) + rxt(k,505)*y(k,122)
         mat(k,1151) = mat(k,1151) + rxt(k,470)*y(k,67) + rxt(k,482)*y(k,122) &
                      + rxt(k,356)*y(k,133)
         mat(k,1280) = mat(k,1280) + rxt(k,124)*y(k,91) + rxt(k,357)*y(k,133)
         mat(k,210) = rxt(k,501)*y(k,67)
         mat(k,232) = mat(k,232) + rxt(k,473)*y(k,67)
         mat(k,1501) = mat(k,1501) + rxt(k,439)*y(k,98)
         mat(k,985) = mat(k,985) + rxt(k,346)*y(k,132)
         mat(k,1108) = mat(k,1108) + rxt(k,434)*y(k,98)
         mat(k,1418) = mat(k,1418) + rxt(k,531)*y(k,16) + rxt(k,406)*y(k,39) &
                      + rxt(k,407)*y(k,41) + rxt(k,499)*y(k,43) + rxt(k,470)*y(k,47) &
                      + rxt(k,501)*y(k,51) + rxt(k,473)*y(k,52) + rxt(k,412)*y(k,98)
         mat(k,1192) = mat(k,1192) + rxt(k,414)*y(k,38) + rxt(k,413)*y(k,98) &
                      + rxt(k,348)*y(k,132)
         mat(k,1460) = rxt(k,455)*y(k,98) + rxt(k,341)*y(k,132)
         mat(k,821) = rxt(k,124)*y(k,49)
         mat(k,664) = rxt(k,129)*y(k,38)
         mat(k,1545) = rxt(k,138)*y(k,96)
         mat(k,85) = rxt(k,138)*y(k,95) + rxt(k,139)*y(k,137)
         mat(k,336) = rxt(k,189)*y(k,100) + rxt(k,203)*y(k,102) + rxt(k,172)*y(k,104) &
                      + rxt(k,242)*y(k,115) + rxt(k,224)*y(k,116) + rxt(k,207) &
                      *y(k,119) + rxt(k,149)*y(k,126)
         mat(k,748) = mat(k,748) + 2.000_r8*rxt(k,409)*y(k,38) + rxt(k,439)*y(k,60) &
                      + rxt(k,434)*y(k,62) + rxt(k,412)*y(k,67) + rxt(k,413)*y(k,69) &
                      + rxt(k,455)*y(k,86)
         mat(k,1027) = rxt(k,189)*y(k,97)
         mat(k,1694) = rxt(k,203)*y(k,97)
         mat(k,1585) = rxt(k,172)*y(k,97)
         mat(k,897) = rxt(k,242)*y(k,97)
         mat(k,939) = rxt(k,224)*y(k,97)
         mat(k,1068) = mat(k,1068) + rxt(k,207)*y(k,97)
         mat(k,639) = rxt(k,389)*y(k,39) + rxt(k,505)*y(k,43) + rxt(k,482)*y(k,47) &
                      + 2.000_r8*rxt(k,390)*y(k,137)
         mat(k,1359) = rxt(k,149)*y(k,97)
         mat(k,190) = rxt(k,328)*y(k,137)
         mat(k,794) = mat(k,794) + .300_r8*rxt(k,535)*y(k,21)
         mat(k,1617) = rxt(k,342)*y(k,26) + rxt(k,346)*y(k,61) + rxt(k,348)*y(k,69) &
                      + rxt(k,341)*y(k,86)
         mat(k,1652) = rxt(k,356)*y(k,47) + rxt(k,357)*y(k,49) + rxt(k,355)*y(k,137)
         mat(k,1831) = rxt(k,509)*y(k,37) + rxt(k,139)*y(k,96) + 2.000_r8*rxt(k,390) &
                      *y(k,122) + rxt(k,328)*y(k,127) + rxt(k,355)*y(k,133)
         mat(k,1637) = -(rxt(k,341)*y(k,86) + rxt(k,342)*y(k,26) + rxt(k,344)*y(k,38) &
                      + rxt(k,345)*y(k,47) + rxt(k,346)*y(k,61) + rxt(k,347)*y(k,67) &
                      + rxt(k,348)*y(k,69))
         mat(k,1480) = -rxt(k,341)*y(k,132)
         mat(k,1758) = -rxt(k,342)*y(k,132)
         mat(k,1793) = -rxt(k,344)*y(k,132)
         mat(k,1171) = -rxt(k,345)*y(k,132)
         mat(k,1005) = -rxt(k,346)*y(k,132)
         mat(k,1438) = -rxt(k,347)*y(k,132)
         mat(k,1212) = -rxt(k,348)*y(k,132)
         mat(k,1793) = mat(k,1793) + rxt(k,117)*y(k,91) + rxt(k,287)*y(k,108) &
                      + rxt(k,331)*y(k,128)
         mat(k,498) = rxt(k,354)*y(k,133)
         mat(k,841) = rxt(k,117)*y(k,38)
         mat(k,870) = rxt(k,287)*y(k,38)
         mat(k,449) = rxt(k,331)*y(k,38)
         mat(k,1672) = rxt(k,354)*y(k,39) + rxt(k,355)*y(k,137)
         mat(k,1851) = rxt(k,355)*y(k,133)
         mat(k,1673) = -(rxt(k,136)*y(k,60) + rxt(k,350)*y(k,86) + rxt(k,351)*y(k,26) &
                      + (rxt(k,353) + rxt(k,354)) * y(k,39) + rxt(k,355)*y(k,137) &
                      + rxt(k,356)*y(k,47) + rxt(k,357)*y(k,49) + rxt(k,359)*y(k,61) &
                      + rxt(k,360)*y(k,67) + rxt(k,361)*y(k,123) + rxt(k,362)*y(k,68) &
                      + rxt(k,363)*y(k,69))
         mat(k,1522) = -rxt(k,136)*y(k,133)
         mat(k,1481) = -rxt(k,350)*y(k,133)
         mat(k,1759) = -rxt(k,351)*y(k,133)
         mat(k,499) = -(rxt(k,353) + rxt(k,354)) * y(k,133)
         mat(k,1852) = -rxt(k,355)*y(k,133)
         mat(k,1172) = -rxt(k,356)*y(k,133)
         mat(k,1301) = -rxt(k,357)*y(k,133)
         mat(k,1006) = -rxt(k,359)*y(k,133)
         mat(k,1439) = -rxt(k,360)*y(k,133)
         mat(k,253) = -rxt(k,361)*y(k,133)
         mat(k,1262) = -rxt(k,362)*y(k,133)
         mat(k,1213) = -rxt(k,363)*y(k,133)
         mat(k,1439) = mat(k,1439) + rxt(k,325)*y(k,125)
         mat(k,1213) = mat(k,1213) + rxt(k,134)*y(k,95)
         mat(k,1564) = rxt(k,134)*y(k,69)
         mat(k,1336) = rxt(k,325)*y(k,67)
         mat(k,264) = -(rxt(k,597)*y(k,106) + rxt(k,598)*y(k,68))
         mat(k,217) = -rxt(k,597)*y(k,134)
         mat(k,1228) = -rxt(k,598)*y(k,134)
         mat(k,1402) = rxt(k,588)*y(k,107) + rxt(k,589)*y(k,121) + rxt(k,601)*y(k,135) &
                      + rxt(k,607)*y(k,136)
         mat(k,1540) = rxt(k,599)*y(k,135) + rxt(k,604)*y(k,136)
         mat(k,155) = rxt(k,588)*y(k,67)
         mat(k,224) = rxt(k,589)*y(k,67)
         mat(k,130) = rxt(k,601)*y(k,67) + rxt(k,599)*y(k,95)
         mat(k,125) = rxt(k,607)*y(k,67) + rxt(k,604)*y(k,95)
         mat(k,128) = -(rxt(k,599)*y(k,95) + rxt(k,601)*y(k,67) + rxt(k,602)*y(k,68))
         mat(k,1532) = -rxt(k,599)*y(k,135)
         mat(k,1389) = -rxt(k,601)*y(k,135)
         mat(k,1221) = -rxt(k,602)*y(k,135)
         mat(k,1532) = mat(k,1532) + rxt(k,603)*y(k,136)
         mat(k,122) = rxt(k,603)*y(k,95)
         mat(k,121) = -((rxt(k,603) + rxt(k,604)) * y(k,95) + rxt(k,607)*y(k,67))
         mat(k,1531) = -(rxt(k,603) + rxt(k,604)) * y(k,136)
         mat(k,1388) = -rxt(k,607)*y(k,136)
         mat(k,1856) = -(rxt(k,107)*y(k,87) + rxt(k,118)*y(k,93) + rxt(k,119)*y(k,91) &
                      + rxt(k,139)*y(k,96) + rxt(k,140)*y(k,101) + rxt(k,143)*y(k,103) &
                      + rxt(k,288)*y(k,108) + rxt(k,294)*y(k,112) + rxt(k,296) &
                      *y(k,110) + rxt(k,305)*y(k,117) + rxt(k,307)*y(k,118) + rxt(k,313) &
                      *y(k,120) + (rxt(k,327) + rxt(k,328)) * y(k,127) + rxt(k,329) &
                      *y(k,126) + rxt(k,337)*y(k,130) + rxt(k,355)*y(k,133) + rxt(k,365) &
                      *y(k,115) + rxt(k,366)*y(k,119) + rxt(k,367)*y(k,116) + rxt(k,372) &
                      *y(k,105) + rxt(k,374)*y(k,99) + rxt(k,376)*y(k,100) + rxt(k,378) &
                      *y(k,102) + rxt(k,390)*y(k,122) + rxt(k,509)*y(k,37) + rxt(k,555) &
                      *y(k,78))
         mat(k,708) = -rxt(k,107)*y(k,137)
         mat(k,556) = -rxt(k,118)*y(k,137)
         mat(k,846) = -rxt(k,119)*y(k,137)
         mat(k,88) = -rxt(k,139)*y(k,137)
         mat(k,53) = -rxt(k,140)*y(k,137)
         mat(k,57) = -rxt(k,143)*y(k,137)
         mat(k,875) = -rxt(k,288)*y(k,137)
         mat(k,576) = -rxt(k,294)*y(k,137)
         mat(k,733) = -rxt(k,296)*y(k,137)
         mat(k,94) = -rxt(k,305)*y(k,137)
         mat(k,120) = -rxt(k,307)*y(k,137)
         mat(k,114) = -rxt(k,313)*y(k,137)
         mat(k,196) = -(rxt(k,327) + rxt(k,328)) * y(k,137)
         mat(k,1384) = -rxt(k,329)*y(k,137)
         mat(k,361) = -rxt(k,337)*y(k,137)
         mat(k,1677) = -rxt(k,355)*y(k,137)
         mat(k,922) = -rxt(k,365)*y(k,137)
         mat(k,1093) = -rxt(k,366)*y(k,137)
         mat(k,964) = -rxt(k,367)*y(k,137)
         mat(k,185) = -rxt(k,372)*y(k,137)
         mat(k,167) = -rxt(k,374)*y(k,137)
         mat(k,1052) = -rxt(k,376)*y(k,137)
         mat(k,1719) = -rxt(k,378)*y(k,137)
         mat(k,657) = -rxt(k,390)*y(k,137)
         mat(k,205) = -rxt(k,509)*y(k,137)
         mat(k,49) = -rxt(k,555)*y(k,137)
         mat(k,624) = rxt(k,532)*y(k,131)
         mat(k,82) = rxt(k,535)*y(k,131)
         mat(k,1798) = rxt(k,410)*y(k,98) + rxt(k,344)*y(k,132)
         mat(k,502) = rxt(k,415)*y(k,131) + rxt(k,353)*y(k,133)
         mat(k,152) = rxt(k,416)*y(k,131)
         mat(k,245) = rxt(k,500)*y(k,131)
         mat(k,1176) = (rxt(k,571)+rxt(k,576))*y(k,51) + (rxt(k,564)+rxt(k,570) &
                       +rxt(k,575))*y(k,52) + rxt(k,106)*y(k,88) + rxt(k,471)*y(k,131) &
                      + rxt(k,345)*y(k,132)
         mat(k,1305) = rxt(k,295)*y(k,112) + rxt(k,447)*y(k,131)
         mat(k,107) = rxt(k,423)*y(k,131)
         mat(k,214) = (rxt(k,571)+rxt(k,576))*y(k,47)
         mat(k,237) = (rxt(k,564)+rxt(k,570)+rxt(k,575))*y(k,47) + rxt(k,474)*y(k,131)
         mat(k,179) = rxt(k,340)*y(k,131)
         mat(k,308) = rxt(k,293)*y(k,111)
         mat(k,1526) = rxt(k,123)*y(k,93)
         mat(k,1010) = rxt(k,120)*y(k,93)
         mat(k,708) = mat(k,708) + 3.000_r8*rxt(k,195)*y(k,100) + 4.000_r8*rxt(k,147) &
                      *y(k,102) + 5.000_r8*rxt(k,177)*y(k,104) + 2.000_r8*rxt(k,230) &
                      *y(k,116) + rxt(k,212)*y(k,119)
         mat(k,423) = rxt(k,106)*y(k,47) + 4.000_r8*rxt(k,199)*y(k,100) &
                      + 5.000_r8*rxt(k,164)*y(k,102) + 6.000_r8*rxt(k,182)*y(k,104) &
                      + rxt(k,252)*y(k,115) + 3.000_r8*rxt(k,234)*y(k,116) &
                      + 2.000_r8*rxt(k,217)*y(k,119) + rxt(k,159)*y(k,126)
         mat(k,407) = 3.000_r8*rxt(k,187)*y(k,100) + (4.000_r8*rxt(k,181) &
                       +4.000_r8*rxt(k,267))*y(k,102) + (5.000_r8*rxt(k,169) &
                       +5.000_r8*rxt(k,268))*y(k,104) + 2.000_r8*rxt(k,222)*y(k,116) &
                      + rxt(k,205)*y(k,119)
         mat(k,846) = mat(k,846) + 3.000_r8*rxt(k,186)*y(k,100) + (4.000_r8*rxt(k,170) &
                       +4.000_r8*rxt(k,258))*y(k,102) + (5.000_r8*rxt(k,168) &
                       +5.000_r8*rxt(k,265))*y(k,104) + 2.000_r8*rxt(k,221)*y(k,116) &
                      + rxt(k,204)*y(k,119)
         mat(k,375) = 5.000_r8*rxt(k,194)*y(k,100) + (6.000_r8*rxt(k,247) &
                       +6.000_r8*rxt(k,271))*y(k,102) + (7.000_r8*rxt(k,176) &
                       +7.000_r8*rxt(k,259))*y(k,104) + 2.000_r8*rxt(k,246)*y(k,115) &
                      + 4.000_r8*rxt(k,229)*y(k,116) + 3.000_r8*rxt(k,211)*y(k,119) &
                      + 2.000_r8*rxt(k,153)*y(k,126)
         mat(k,556) = mat(k,556) + rxt(k,123)*y(k,60) + rxt(k,120)*y(k,61) &
                      + 4.000_r8*rxt(k,196)*y(k,100) + (5.000_r8*rxt(k,158) &
                       +5.000_r8*rxt(k,260))*y(k,102) + (6.000_r8*rxt(k,178) &
                       +6.000_r8*rxt(k,261))*y(k,104) + rxt(k,249)*y(k,115) &
                      + 3.000_r8*rxt(k,231)*y(k,116) + 2.000_r8*rxt(k,213)*y(k,119) &
                      + rxt(k,155)*y(k,126)
         mat(k,687) = 3.000_r8*rxt(k,191)*y(k,100) + 4.000_r8*rxt(k,225)*y(k,102) &
                      + 5.000_r8*rxt(k,174)*y(k,104) + 2.000_r8*rxt(k,227)*y(k,116) &
                      + rxt(k,209)*y(k,119)
         mat(k,1568) = rxt(k,138)*y(k,96) + 2.000_r8*rxt(k,382)*y(k,99) &
                      + 3.000_r8*rxt(k,383)*y(k,100) + 4.000_r8*rxt(k,141)*y(k,102) &
                      + 5.000_r8*rxt(k,144)*y(k,104) + rxt(k,381)*y(k,105) &
                      + 2.000_r8*rxt(k,303)*y(k,116) + 3.000_r8*rxt(k,304)*y(k,117) &
                      + rxt(k,308)*y(k,119) + rxt(k,326)*y(k,127)
         mat(k,88) = mat(k,88) + rxt(k,138)*y(k,95)
         mat(k,347) = 3.000_r8*rxt(k,189)*y(k,100) + 4.000_r8*rxt(k,203)*y(k,102) &
                      + 5.000_r8*rxt(k,172)*y(k,104) + 2.000_r8*rxt(k,224)*y(k,116) &
                      + rxt(k,207)*y(k,119)
         mat(k,767) = rxt(k,410)*y(k,38) + rxt(k,417)*y(k,131)
         mat(k,167) = mat(k,167) + 2.000_r8*rxt(k,382)*y(k,95)
         mat(k,1052) = mat(k,1052) + 3.000_r8*rxt(k,195)*y(k,87) + 4.000_r8*rxt(k,199) &
                      *y(k,88) + 3.000_r8*rxt(k,187)*y(k,89) + 3.000_r8*rxt(k,186) &
                      *y(k,91) + 5.000_r8*rxt(k,194)*y(k,92) + 4.000_r8*rxt(k,196) &
                      *y(k,93) + 3.000_r8*rxt(k,191)*y(k,94) + 3.000_r8*rxt(k,383) &
                      *y(k,95) + 3.000_r8*rxt(k,189)*y(k,97) + 3.000_r8*rxt(k,201) &
                      *y(k,108) + 4.000_r8*rxt(k,197)*y(k,109) + 3.000_r8*rxt(k,188) &
                      *y(k,110) + 5.000_r8*rxt(k,200)*y(k,111) + 4.000_r8*rxt(k,193) &
                      *y(k,112) + 3.000_r8*rxt(k,198)*y(k,113) + 3.000_r8*rxt(k,185) &
                      *y(k,114) + 3.000_r8*rxt(k,190)*y(k,125)
         mat(k,1719) = mat(k,1719) + 4.000_r8*rxt(k,147)*y(k,87) + 5.000_r8*rxt(k,164) &
                      *y(k,88) + (4.000_r8*rxt(k,181)+4.000_r8*rxt(k,267))*y(k,89) + ( &
                      + 4.000_r8*rxt(k,170)+4.000_r8*rxt(k,258))*y(k,91) + ( &
                      + 6.000_r8*rxt(k,247)+6.000_r8*rxt(k,271))*y(k,92) + ( &
                      + 5.000_r8*rxt(k,158)+5.000_r8*rxt(k,260))*y(k,93) &
                      + 4.000_r8*rxt(k,225)*y(k,94) + 4.000_r8*rxt(k,141)*y(k,95) &
                      + 4.000_r8*rxt(k,203)*y(k,97) + 4.000_r8*rxt(k,166)*y(k,108) &
                      + 5.000_r8*rxt(k,162)*y(k,109) + (4.000_r8*rxt(k,192) &
                       +4.000_r8*rxt(k,264))*y(k,110) + 6.000_r8*rxt(k,165)*y(k,111) + ( &
                      + 5.000_r8*rxt(k,236)+5.000_r8*rxt(k,262))*y(k,112) &
                      + 4.000_r8*rxt(k,163)*y(k,113) + (4.000_r8*rxt(k,146) &
                       +4.000_r8*rxt(k,269))*y(k,114) + 4.000_r8*rxt(k,214)*y(k,125)
         mat(k,1610) = 5.000_r8*rxt(k,177)*y(k,87) + 6.000_r8*rxt(k,182)*y(k,88) + ( &
                      + 5.000_r8*rxt(k,169)+5.000_r8*rxt(k,268))*y(k,89) + ( &
                      + 5.000_r8*rxt(k,168)+5.000_r8*rxt(k,265))*y(k,91) + ( &
                      + 7.000_r8*rxt(k,176)+7.000_r8*rxt(k,259))*y(k,92) + ( &
                      + 6.000_r8*rxt(k,178)+6.000_r8*rxt(k,261))*y(k,93) &
                      + 5.000_r8*rxt(k,174)*y(k,94) + 5.000_r8*rxt(k,144)*y(k,95) &
                      + 5.000_r8*rxt(k,172)*y(k,97) + 5.000_r8*rxt(k,184)*y(k,108) &
                      + 6.000_r8*rxt(k,179)*y(k,109) + (5.000_r8*rxt(k,171) &
                       +5.000_r8*rxt(k,266))*y(k,110) + 7.000_r8*rxt(k,183)*y(k,111) + ( &
                      + 6.000_r8*rxt(k,175)+6.000_r8*rxt(k,263))*y(k,112) &
                      + 5.000_r8*rxt(k,180)*y(k,113) + (5.000_r8*rxt(k,167) &
                       +5.000_r8*rxt(k,270))*y(k,114) + 5.000_r8*rxt(k,173)*y(k,125)
         mat(k,185) = mat(k,185) + rxt(k,381)*y(k,95)
         mat(k,875) = mat(k,875) + 3.000_r8*rxt(k,201)*y(k,100) + 4.000_r8*rxt(k,166) &
                      *y(k,102) + 5.000_r8*rxt(k,184)*y(k,104) + 2.000_r8*rxt(k,237) &
                      *y(k,116) + rxt(k,219)*y(k,119)
         mat(k,437) = 4.000_r8*rxt(k,197)*y(k,100) + 5.000_r8*rxt(k,162)*y(k,102) &
                      + 6.000_r8*rxt(k,179)*y(k,104) + rxt(k,250)*y(k,115) &
                      + 3.000_r8*rxt(k,232)*y(k,116) + 2.000_r8*rxt(k,215)*y(k,119) &
                      + rxt(k,156)*y(k,126)
         mat(k,733) = mat(k,733) + 3.000_r8*rxt(k,188)*y(k,100) + (4.000_r8*rxt(k,192) &
                       +4.000_r8*rxt(k,264))*y(k,102) + (5.000_r8*rxt(k,171) &
                       +5.000_r8*rxt(k,266))*y(k,104) + 2.000_r8*rxt(k,223)*y(k,116) &
                      + rxt(k,206)*y(k,119)
         mat(k,536) = rxt(k,293)*y(k,56) + 5.000_r8*rxt(k,200)*y(k,100) &
                      + 6.000_r8*rxt(k,165)*y(k,102) + 7.000_r8*rxt(k,183)*y(k,104) &
                      + 2.000_r8*rxt(k,253)*y(k,115) + 4.000_r8*rxt(k,235)*y(k,116) &
                      + 3.000_r8*rxt(k,218)*y(k,119) + 2.000_r8*rxt(k,160)*y(k,126)
         mat(k,576) = mat(k,576) + rxt(k,295)*y(k,49) + 4.000_r8*rxt(k,193)*y(k,100) + ( &
                      + 5.000_r8*rxt(k,236)+5.000_r8*rxt(k,262))*y(k,102) + ( &
                      + 6.000_r8*rxt(k,175)+6.000_r8*rxt(k,263))*y(k,104) + rxt(k,245) &
                      *y(k,115) + 3.000_r8*rxt(k,228)*y(k,116) + 2.000_r8*rxt(k,210) &
                      *y(k,119) + rxt(k,152)*y(k,126)
         mat(k,480) = 3.000_r8*rxt(k,198)*y(k,100) + 4.000_r8*rxt(k,163)*y(k,102) &
                      + 5.000_r8*rxt(k,180)*y(k,104) + 2.000_r8*rxt(k,233)*y(k,116) &
                      + rxt(k,216)*y(k,119)
         mat(k,516) = 3.000_r8*rxt(k,185)*y(k,100) + (4.000_r8*rxt(k,146) &
                       +4.000_r8*rxt(k,269))*y(k,102) + (5.000_r8*rxt(k,167) &
                       +5.000_r8*rxt(k,270))*y(k,104) + 2.000_r8*rxt(k,220)*y(k,116) &
                      + rxt(k,202)*y(k,119)
         mat(k,922) = mat(k,922) + rxt(k,252)*y(k,88) + 2.000_r8*rxt(k,246)*y(k,92) &
                      + rxt(k,249)*y(k,93) + rxt(k,250)*y(k,109) + 2.000_r8*rxt(k,253) &
                      *y(k,111) + rxt(k,245)*y(k,112)
         mat(k,964) = mat(k,964) + 2.000_r8*rxt(k,230)*y(k,87) + 3.000_r8*rxt(k,234) &
                      *y(k,88) + 2.000_r8*rxt(k,222)*y(k,89) + 2.000_r8*rxt(k,221) &
                      *y(k,91) + 4.000_r8*rxt(k,229)*y(k,92) + 3.000_r8*rxt(k,231) &
                      *y(k,93) + 2.000_r8*rxt(k,227)*y(k,94) + 2.000_r8*rxt(k,303) &
                      *y(k,95) + 2.000_r8*rxt(k,224)*y(k,97) + 2.000_r8*rxt(k,237) &
                      *y(k,108) + 3.000_r8*rxt(k,232)*y(k,109) + 2.000_r8*rxt(k,223) &
                      *y(k,110) + 4.000_r8*rxt(k,235)*y(k,111) + 3.000_r8*rxt(k,228) &
                      *y(k,112) + 2.000_r8*rxt(k,233)*y(k,113) + 2.000_r8*rxt(k,220) &
                      *y(k,114) + 2.000_r8*rxt(k,226)*y(k,125)
         mat(k,94) = mat(k,94) + 3.000_r8*rxt(k,304)*y(k,95)
         mat(k,1093) = mat(k,1093) + rxt(k,212)*y(k,87) + 2.000_r8*rxt(k,217)*y(k,88) &
                      + rxt(k,205)*y(k,89) + rxt(k,204)*y(k,91) + 3.000_r8*rxt(k,211) &
                      *y(k,92) + 2.000_r8*rxt(k,213)*y(k,93) + rxt(k,209)*y(k,94) &
                      + rxt(k,308)*y(k,95) + rxt(k,207)*y(k,97) + rxt(k,219)*y(k,108) &
                      + 2.000_r8*rxt(k,215)*y(k,109) + rxt(k,206)*y(k,110) &
                      + 3.000_r8*rxt(k,218)*y(k,111) + 2.000_r8*rxt(k,210)*y(k,112) &
                      + rxt(k,216)*y(k,113) + rxt(k,202)*y(k,114) + rxt(k,208) &
                      *y(k,125)
         mat(k,1340) = 3.000_r8*rxt(k,190)*y(k,100) + 4.000_r8*rxt(k,214)*y(k,102) &
                      + 5.000_r8*rxt(k,173)*y(k,104) + 2.000_r8*rxt(k,226)*y(k,116) &
                      + rxt(k,208)*y(k,119)
         mat(k,1384) = mat(k,1384) + rxt(k,159)*y(k,88) + 2.000_r8*rxt(k,153)*y(k,92) &
                      + rxt(k,155)*y(k,93) + rxt(k,156)*y(k,109) + 2.000_r8*rxt(k,160) &
                      *y(k,111) + rxt(k,152)*y(k,112)
         mat(k,196) = mat(k,196) + rxt(k,326)*y(k,95)
         mat(k,814) = rxt(k,532)*y(k,16) + rxt(k,535)*y(k,21) + rxt(k,415)*y(k,39) &
                      + rxt(k,416)*y(k,41) + rxt(k,500)*y(k,43) + rxt(k,471)*y(k,47) &
                      + rxt(k,447)*y(k,49) + rxt(k,423)*y(k,50) + rxt(k,474)*y(k,52) &
                      + rxt(k,340)*y(k,53) + rxt(k,417)*y(k,98) + 2.000_r8*rxt(k,420) &
                      *y(k,131)
         mat(k,1642) = rxt(k,344)*y(k,38) + rxt(k,345)*y(k,47)
         mat(k,1677) = mat(k,1677) + rxt(k,353)*y(k,39)
      end do
      end subroutine nlnmat08
      subroutine nlnmat_finit( avec_len, mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: dti(veclen)
      real(r8), intent(in) :: lmat(veclen,nzcnt)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = 1,avec_len
         mat(k, 1) = lmat(k, 1)
         mat(k, 2) = lmat(k, 2)
         mat(k, 3) = lmat(k, 3)
         mat(k, 4) = lmat(k, 4)
         mat(k, 5) = lmat(k, 5)
         mat(k, 6) = lmat(k, 6)
         mat(k, 7) = lmat(k, 7)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = lmat(k, 10)
         mat(k, 11) = lmat(k, 11)
         mat(k, 12) = lmat(k, 12)
         mat(k, 13) = lmat(k, 13)
         mat(k, 14) = lmat(k, 14)
         mat(k, 15) = lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = lmat(k, 18)
         mat(k, 19) = lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = lmat(k, 22)
         mat(k, 23) = lmat(k, 23)
         mat(k, 24) = lmat(k, 24)
         mat(k, 25) = lmat(k, 25)
         mat(k, 26) = lmat(k, 26)
         mat(k, 27) = mat(k, 27) + lmat(k, 27)
         mat(k, 28) = mat(k, 28) + lmat(k, 28)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 35) = mat(k, 35) + lmat(k, 35)
         mat(k, 36) = lmat(k, 36)
         mat(k, 37) = lmat(k, 37)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 44) = lmat(k, 44)
         mat(k, 46) = mat(k, 46) + lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 64) = mat(k, 64) + lmat(k, 64)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 83) = mat(k, 83) + lmat(k, 83)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 99) = lmat(k, 99)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = lmat(k, 116)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = lmat(k, 144)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 163) = lmat(k, 163)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 178) = lmat(k, 178)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 193) = lmat(k, 193)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 198) = mat(k, 198) + lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 259) = lmat(k, 259)
         mat(k, 261) = lmat(k, 261)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 265) = lmat(k, 265)
         mat(k, 269) = lmat(k, 269)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 282) = mat(k, 282) + lmat(k, 282)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 287) = mat(k, 287) + lmat(k, 287)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 301) = lmat(k, 301)
         mat(k, 302) = lmat(k, 302)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 304) = lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 309) = lmat(k, 309)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 321) = lmat(k, 321)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 363) = lmat(k, 363)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 377) = lmat(k, 377)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = lmat(k, 394)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 440) = lmat(k, 440)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 450) = lmat(k, 450)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 484) = mat(k, 484) + lmat(k, 484)
         mat(k, 503) = mat(k, 503) + lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 542) = lmat(k, 542)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 560) = mat(k, 560) + lmat(k, 560)
         mat(k, 561) = lmat(k, 561)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 604) = lmat(k, 604)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 631) = lmat(k, 631)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 634) = mat(k, 634) + lmat(k, 634)
         mat(k, 636) = lmat(k, 636)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 638) = lmat(k, 638)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 675) = lmat(k, 675)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 713) = mat(k, 713) + lmat(k, 713)
         mat(k, 720) = mat(k, 720) + lmat(k, 720)
         mat(k, 729) = lmat(k, 729)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 747) = mat(k, 747) + lmat(k, 747)
         mat(k, 768) = lmat(k, 768)
         mat(k, 769) = lmat(k, 769)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 788) = mat(k, 788) + lmat(k, 788)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 794) = mat(k, 794) + lmat(k, 794)
         mat(k, 807) = mat(k, 807) + lmat(k, 807)
         mat(k, 814) = mat(k, 814) + lmat(k, 814)
         mat(k, 822) = mat(k, 822) + lmat(k, 822)
         mat(k, 842) = lmat(k, 842)
         mat(k, 852) = mat(k, 852) + lmat(k, 852)
         mat(k, 855) = mat(k, 855) + lmat(k, 855)
         mat(k, 868) = lmat(k, 868)
         mat(k, 876) = lmat(k, 876)
         mat(k, 877) = lmat(k, 877)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k, 990) = mat(k, 990) + lmat(k, 990)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k,1000) = mat(k,1000) + lmat(k,1000)
         mat(k,1002) = mat(k,1002) + lmat(k,1002)
         mat(k,1011) = lmat(k,1011)
         mat(k,1033) = mat(k,1033) + lmat(k,1033)
         mat(k,1052) = mat(k,1052) + lmat(k,1052)
         mat(k,1075) = mat(k,1075) + lmat(k,1075)
         mat(k,1112) = mat(k,1112) + lmat(k,1112)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1118) = mat(k,1118) + lmat(k,1118)
         mat(k,1119) = mat(k,1119) + lmat(k,1119)
         mat(k,1122) = mat(k,1122) + lmat(k,1122)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1160) = mat(k,1160) + lmat(k,1160)
         mat(k,1167) = mat(k,1167) + lmat(k,1167)
         mat(k,1175) = mat(k,1175) + lmat(k,1175)
         mat(k,1179) = mat(k,1179) + lmat(k,1179)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1203) = mat(k,1203) + lmat(k,1203)
         mat(k,1207) = mat(k,1207) + lmat(k,1207)
         mat(k,1220) = lmat(k,1220)
         mat(k,1221) = mat(k,1221) + lmat(k,1221)
         mat(k,1228) = mat(k,1228) + lmat(k,1228)
         mat(k,1236) = mat(k,1236) + lmat(k,1236)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1255) = mat(k,1255) + lmat(k,1255)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1259) = mat(k,1259) + lmat(k,1259)
         mat(k,1280) = mat(k,1280) + lmat(k,1280)
         mat(k,1285) = lmat(k,1285)
         mat(k,1292) = mat(k,1292) + lmat(k,1292)
         mat(k,1311) = lmat(k,1311)
         mat(k,1326) = mat(k,1326) + lmat(k,1326)
         mat(k,1328) = mat(k,1328) + lmat(k,1328)
         mat(k,1333) = mat(k,1333) + lmat(k,1333)
         mat(k,1362) = mat(k,1362) + lmat(k,1362)
         mat(k,1373) = mat(k,1373) + lmat(k,1373)
         mat(k,1376) = mat(k,1376) + lmat(k,1376)
         mat(k,1388) = mat(k,1388) + lmat(k,1388)
         mat(k,1389) = mat(k,1389) + lmat(k,1389)
         mat(k,1402) = mat(k,1402) + lmat(k,1402)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1436) = mat(k,1436) + lmat(k,1436)
         mat(k,1445) = mat(k,1445) + lmat(k,1445)
         mat(k,1450) = lmat(k,1450)
         mat(k,1451) = lmat(k,1451)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1469) = mat(k,1469) + lmat(k,1469)
         mat(k,1476) = mat(k,1476) + lmat(k,1476)
         mat(k,1488) = mat(k,1488) + lmat(k,1488)
         mat(k,1504) = mat(k,1504) + lmat(k,1504)
         mat(k,1516) = mat(k,1516) + lmat(k,1516)
         mat(k,1518) = mat(k,1518) + lmat(k,1518)
         mat(k,1519) = mat(k,1519) + lmat(k,1519)
         mat(k,1561) = mat(k,1561) + lmat(k,1561)
         mat(k,1604) = mat(k,1604) + lmat(k,1604)
         mat(k,1607) = lmat(k,1607)
         mat(k,1610) = mat(k,1610) + lmat(k,1610)
         mat(k,1612) = lmat(k,1612)
         mat(k,1617) = mat(k,1617) + lmat(k,1617)
         mat(k,1635) = mat(k,1635) + lmat(k,1635)
         mat(k,1637) = mat(k,1637) + lmat(k,1637)
         mat(k,1653) = lmat(k,1653)
         mat(k,1667) = mat(k,1667) + lmat(k,1667)
         mat(k,1670) = mat(k,1670) + lmat(k,1670)
         mat(k,1673) = mat(k,1673) + lmat(k,1673)
         mat(k,1700) = lmat(k,1700)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1753) = mat(k,1753) + lmat(k,1753)
         mat(k,1754) = mat(k,1754) + lmat(k,1754)
         mat(k,1761) = mat(k,1761) + lmat(k,1761)
         mat(k,1797) = mat(k,1797) + lmat(k,1797)
         mat(k,1821) = lmat(k,1821)
         mat(k,1827) = mat(k,1827) + lmat(k,1827)
         mat(k,1831) = mat(k,1831) + lmat(k,1831)
         mat(k,1846) = lmat(k,1846)
         mat(k,1855) = lmat(k,1855)
         mat(k,1856) = mat(k,1856) + lmat(k,1856)
         mat(k, 182) = 0._r8
         mat(k, 187) = 0._r8
         mat(k, 191) = 0._r8
         mat(k, 195) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 260) = 0._r8
         mat(k, 266) = 0._r8
         mat(k, 267) = 0._r8
         mat(k, 268) = 0._r8
         mat(k, 273) = 0._r8
         mat(k, 274) = 0._r8
         mat(k, 286) = 0._r8
         mat(k, 310) = 0._r8
         mat(k, 312) = 0._r8
         mat(k, 313) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 334) = 0._r8
         mat(k, 351) = 0._r8
         mat(k, 352) = 0._r8
         mat(k, 355) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 360) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 392) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 464) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 490) = 0._r8
         mat(k, 492) = 0._r8
         mat(k, 493) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 524) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 564) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 594) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 600) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 610) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 618) = 0._r8
         mat(k, 619) = 0._r8
         mat(k, 620) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 630) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 646) = 0._r8
         mat(k, 647) = 0._r8
         mat(k, 648) = 0._r8
         mat(k, 652) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 680) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 763) = 0._r8
         mat(k, 764) = 0._r8
         mat(k, 780) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 805) = 0._r8
         mat(k, 809) = 0._r8
         mat(k, 810) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 830) = 0._r8
         mat(k, 831) = 0._r8
         mat(k, 839) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 863) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 903) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 949) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 975) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 988) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 991) = 0._r8
         mat(k, 992) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1098) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1100) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1120) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1132) = 0._r8
         mat(k,1141) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1266) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1363) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1477) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1485) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1507) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1548) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1584) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1589) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1602) = 0._r8
         mat(k,1605) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1608) = 0._r8
         mat(k,1615) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1623) = 0._r8
         mat(k,1624) = 0._r8
         mat(k,1625) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1629) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1631) = 0._r8
         mat(k,1634) = 0._r8
         mat(k,1636) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1656) = 0._r8
         mat(k,1658) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1693) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1701) = 0._r8
         mat(k,1704) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1744) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1748) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1777) = 0._r8
         mat(k,1778) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1789) = 0._r8
         mat(k,1792) = 0._r8
         mat(k,1794) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1813) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1820) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1840) = 0._r8
         mat(k,1841) = 0._r8
         mat(k,1844) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1848) = 0._r8
         mat(k,1849) = 0._r8
         mat(k,1854) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 3) = mat(k, 3) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 6) = mat(k, 6) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 8) = mat(k, 8) - dti(k)
         mat(k, 9) = mat(k, 9) - dti(k)
         mat(k, 10) = mat(k, 10) - dti(k)
         mat(k, 11) = mat(k, 11) - dti(k)
         mat(k, 12) = mat(k, 12) - dti(k)
         mat(k, 13) = mat(k, 13) - dti(k)
         mat(k, 14) = mat(k, 14) - dti(k)
         mat(k, 15) = mat(k, 15) - dti(k)
         mat(k, 16) = mat(k, 16) - dti(k)
         mat(k, 17) = mat(k, 17) - dti(k)
         mat(k, 18) = mat(k, 18) - dti(k)
         mat(k, 19) = mat(k, 19) - dti(k)
         mat(k, 20) = mat(k, 20) - dti(k)
         mat(k, 21) = mat(k, 21) - dti(k)
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 42) = mat(k, 42) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 54) = mat(k, 54) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 69) = mat(k, 69) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 83) = mat(k, 83) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 121) = mat(k, 121) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 198) = mat(k, 198) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 215) = mat(k, 215) - dti(k)
         mat(k, 223) = mat(k, 223) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 255) = mat(k, 255) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 277) = mat(k, 277) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 335) = mat(k, 335) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 484) = mat(k, 484) - dti(k)
         mat(k, 503) = mat(k, 503) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 539) = mat(k, 539) - dti(k)
         mat(k, 560) = mat(k, 560) - dti(k)
         mat(k, 584) = mat(k, 584) - dti(k)
         mat(k, 606) = mat(k, 606) - dti(k)
         mat(k, 637) = mat(k, 637) - dti(k)
         mat(k, 660) = mat(k, 660) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 713) = mat(k, 713) - dti(k)
         mat(k, 747) = mat(k, 747) - dti(k)
         mat(k, 794) = mat(k, 794) - dti(k)
         mat(k, 822) = mat(k, 822) - dti(k)
         mat(k, 852) = mat(k, 852) - dti(k)
         mat(k, 900) = mat(k, 900) - dti(k)
         mat(k, 943) = mat(k, 943) - dti(k)
         mat(k, 990) = mat(k, 990) - dti(k)
         mat(k,1033) = mat(k,1033) - dti(k)
         mat(k,1075) = mat(k,1075) - dti(k)
         mat(k,1115) = mat(k,1115) - dti(k)
         mat(k,1160) = mat(k,1160) - dti(k)
         mat(k,1202) = mat(k,1202) - dti(k)
         mat(k,1252) = mat(k,1252) - dti(k)
         mat(k,1292) = mat(k,1292) - dti(k)
         mat(k,1328) = mat(k,1328) - dti(k)
         mat(k,1373) = mat(k,1373) - dti(k)
         mat(k,1433) = mat(k,1433) - dti(k)
         mat(k,1476) = mat(k,1476) - dti(k)
         mat(k,1518) = mat(k,1518) - dti(k)
         mat(k,1561) = mat(k,1561) - dti(k)
         mat(k,1604) = mat(k,1604) - dti(k)
         mat(k,1637) = mat(k,1637) - dti(k)
         mat(k,1673) = mat(k,1673) - dti(k)
         mat(k,1716) = mat(k,1716) - dti(k)
         mat(k,1761) = mat(k,1761) - dti(k)
         mat(k,1797) = mat(k,1797) - dti(k)
         mat(k,1856) = mat(k,1856) - dti(k)
      end do
      end subroutine nlnmat_finit
      subroutine nlnmat( avec_len, mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: avec_len
      real(r8), intent(in) :: dti(veclen)
      real(r8), intent(in) :: lmat(veclen,nzcnt)
      real(r8), intent(in) :: y(veclen,gas_pcnst)
      real(r8), intent(in) :: rxt(veclen,rxntot)
      real(r8), intent(inout) :: mat(veclen,nzcnt)
      call nlnmat01( avec_len, mat, y, rxt )
      call nlnmat02( avec_len, mat, y, rxt )
      call nlnmat03( avec_len, mat, y, rxt )
      call nlnmat04( avec_len, mat, y, rxt )
      call nlnmat05( avec_len, mat, y, rxt )
      call nlnmat06( avec_len, mat, y, rxt )
      call nlnmat07( avec_len, mat, y, rxt )
      call nlnmat08( avec_len, mat, y, rxt )
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
