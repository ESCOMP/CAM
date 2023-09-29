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
         mat(k,582) = -(rxt(k,486)*y(k,17) + rxt(k,487)*y(k,101) + rxt(k,488)*y(k,73))
         mat(k,735) = -rxt(k,486)*y(k,3)
         mat(k,895) = -rxt(k,487)*y(k,3)
         mat(k,1909) = -rxt(k,488)*y(k,3)
         mat(k,714) = 4.000_r8*rxt(k,489)*y(k,5) + (rxt(k,490)+rxt(k,491))*y(k,28) &
                      + rxt(k,494)*y(k,62) + rxt(k,497)*y(k,70) + rxt(k,548)*y(k,80) &
                      + rxt(k,498)*y(k,133)
         mat(k,57) = rxt(k,476)*y(k,71)
         mat(k,63) = rxt(k,502)*y(k,71)
         mat(k,187) = 2.000_r8*rxt(k,513)*y(k,25) + 2.000_r8*rxt(k,525)*y(k,71) &
                      + 2.000_r8*rxt(k,514)*y(k,133)
         mat(k,223) = rxt(k,515)*y(k,25) + rxt(k,526)*y(k,71) + rxt(k,516)*y(k,133)
         mat(k,162) = 3.000_r8*rxt(k,520)*y(k,25) + 3.000_r8*rxt(k,503)*y(k,71) &
                      + 3.000_r8*rxt(k,521)*y(k,133)
         mat(k,1466) = 2.000_r8*rxt(k,513)*y(k,16) + rxt(k,515)*y(k,18) &
                      + 3.000_r8*rxt(k,520)*y(k,24)
         mat(k,1092) = (rxt(k,490)+rxt(k,491))*y(k,5)
         mat(k,34) = 2.000_r8*rxt(k,504)*y(k,71)
         mat(k,325) = rxt(k,499)*y(k,70) + rxt(k,505)*y(k,71) + rxt(k,500)*y(k,133)
         mat(k,1949) = rxt(k,494)*y(k,5)
         mat(k,1368) = rxt(k,497)*y(k,5) + rxt(k,499)*y(k,45)
         mat(k,783) = rxt(k,476)*y(k,9) + rxt(k,502)*y(k,10) + 2.000_r8*rxt(k,525) &
                      *y(k,16) + rxt(k,526)*y(k,18) + 3.000_r8*rxt(k,503)*y(k,24) &
                      + 2.000_r8*rxt(k,504)*y(k,42) + rxt(k,505)*y(k,45)
         mat(k,477) = rxt(k,548)*y(k,5)
         mat(k,978) = rxt(k,498)*y(k,5) + 2.000_r8*rxt(k,514)*y(k,16) + rxt(k,516) &
                      *y(k,18) + 3.000_r8*rxt(k,521)*y(k,24) + rxt(k,500)*y(k,45)
         mat(k,708) = rxt(k,492)*y(k,28)
         mat(k,1084) = rxt(k,492)*y(k,5)
         mat(k,1987) = (rxt(k,570)+rxt(k,575))*y(k,53)
         mat(k,308) = (rxt(k,570)+rxt(k,575))*y(k,49)
         mat(k,715) = -(4._r8*rxt(k,489)*y(k,5) + (rxt(k,490) + rxt(k,491) + rxt(k,492) &
                      ) * y(k,28) + rxt(k,493)*y(k,101) + rxt(k,494)*y(k,62) + rxt(k,495) &
                      *y(k,63) + rxt(k,497)*y(k,70) + rxt(k,498)*y(k,133) + rxt(k,548) &
                      *y(k,80))
         mat(k,1093) = -(rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,5)
         mat(k,897) = -rxt(k,493)*y(k,5)
         mat(k,1953) = -rxt(k,494)*y(k,5)
         mat(k,1801) = -rxt(k,495)*y(k,5)
         mat(k,1371) = -rxt(k,497)*y(k,5)
         mat(k,980) = -rxt(k,498)*y(k,5)
         mat(k,478) = -rxt(k,548)*y(k,5)
         mat(k,583) = rxt(k,488)*y(k,73)
         mat(k,209) = rxt(k,496)*y(k,70)
         mat(k,326) = rxt(k,506)*y(k,71)
         mat(k,311) = rxt(k,501)*y(k,70)
         mat(k,1371) = mat(k,1371) + rxt(k,496)*y(k,6) + rxt(k,501)*y(k,53)
         mat(k,785) = rxt(k,506)*y(k,45)
         mat(k,1911) = rxt(k,488)*y(k,3)
         mat(k,206) = -(rxt(k,496)*y(k,70))
         mat(k,1349) = -rxt(k,496)*y(k,6)
         mat(k,710) = rxt(k,495)*y(k,63)
         mat(k,1788) = rxt(k,495)*y(k,5)
         mat(k,29) = -(rxt(k,475)*y(k,71))
         mat(k,758) = -rxt(k,475)*y(k,8)
         mat(k,55) = -(rxt(k,476)*y(k,71))
         mat(k,763) = -rxt(k,476)*y(k,9)
         mat(k,60) = -(rxt(k,502)*y(k,71))
         mat(k,764) = -rxt(k,502)*y(k,10)
         mat(k,36) = -(rxt(k,477)*y(k,71))
         mat(k,760) = -rxt(k,477)*y(k,11)
         mat(k,65) = -(rxt(k,478)*y(k,71))
         mat(k,765) = -rxt(k,478)*y(k,12)
         mat(k,40) = -(rxt(k,479)*y(k,71))
         mat(k,761) = -rxt(k,479)*y(k,13)
         mat(k,70) = -(rxt(k,480)*y(k,71))
         mat(k,766) = -rxt(k,480)*y(k,14)
         mat(k,44) = -(rxt(k,481)*y(k,71))
         mat(k,762) = -rxt(k,481)*y(k,15)
         mat(k,186) = -(rxt(k,513)*y(k,25) + rxt(k,514)*y(k,133) + rxt(k,525)*y(k,71))
         mat(k,1457) = -rxt(k,513)*y(k,16)
         mat(k,962) = -rxt(k,514)*y(k,16)
         mat(k,775) = -rxt(k,525)*y(k,16)
         mat(k,738) = -(rxt(k,450)*y(k,25) + rxt(k,486)*y(k,3) + rxt(k,530)*y(k,64) &
                      + rxt(k,531)*y(k,70) + rxt(k,532)*y(k,133))
         mat(k,1469) = -rxt(k,450)*y(k,17)
         mat(k,584) = -rxt(k,486)*y(k,17)
         mat(k,1720) = -rxt(k,530)*y(k,17)
         mat(k,1372) = -rxt(k,531)*y(k,17)
         mat(k,981) = -rxt(k,532)*y(k,17)
         mat(k,433) = rxt(k,457)*y(k,28) + rxt(k,534)*y(k,62)
         mat(k,157) = .300_r8*rxt(k,535)*y(k,133)
         mat(k,419) = (rxt(k,538)+rxt(k,539))*y(k,71)
         mat(k,1094) = rxt(k,457)*y(k,21)
         mat(k,1954) = rxt(k,534)*y(k,21)
         mat(k,786) = (rxt(k,538)+rxt(k,539))*y(k,23)
         mat(k,981) = mat(k,981) + .300_r8*rxt(k,535)*y(k,22)
         mat(k,221) = -(rxt(k,515)*y(k,25) + rxt(k,516)*y(k,133) + rxt(k,526)*y(k,71))
         mat(k,1459) = -rxt(k,515)*y(k,18)
         mat(k,964) = -rxt(k,516)*y(k,18)
         mat(k,776) = -rxt(k,526)*y(k,18)
         mat(k,48) = -(rxt(k,517)*y(k,133))
         mat(k,949) = -rxt(k,517)*y(k,19)
         mat(k,147) = -(rxt(k,518)*y(k,25) + rxt(k,519)*y(k,133))
         mat(k,1455) = -rxt(k,518)*y(k,20)
         mat(k,958) = -rxt(k,519)*y(k,20)
         mat(k,432) = -(rxt(k,457)*y(k,28) + rxt(k,533)*y(k,101) + rxt(k,534)*y(k,62))
         mat(k,1090) = -rxt(k,457)*y(k,21)
         mat(k,894) = -rxt(k,533)*y(k,21)
         mat(k,1946) = -rxt(k,534)*y(k,21)
         mat(k,156) = .700_r8*rxt(k,535)*y(k,133)
         mat(k,417) = rxt(k,451)*y(k,25) + rxt(k,507)*y(k,39) + rxt(k,537)*y(k,71) &
                      + rxt(k,536)*y(k,133)
         mat(k,1465) = rxt(k,451)*y(k,23)
         mat(k,343) = rxt(k,507)*y(k,23)
         mat(k,781) = rxt(k,537)*y(k,23)
         mat(k,975) = .700_r8*rxt(k,535)*y(k,22) + rxt(k,536)*y(k,23)
         mat(k,155) = -(rxt(k,535)*y(k,133))
         mat(k,959) = -rxt(k,535)*y(k,22)
         mat(k,431) = rxt(k,533)*y(k,101)
         mat(k,887) = rxt(k,533)*y(k,21)
         mat(k,416) = -(rxt(k,451)*y(k,25) + rxt(k,507)*y(k,39) + rxt(k,536)*y(k,133) &
                      + (rxt(k,537) + rxt(k,538) + rxt(k,539)) * y(k,71))
         mat(k,1464) = -rxt(k,451)*y(k,23)
         mat(k,342) = -rxt(k,507)*y(k,23)
         mat(k,974) = -rxt(k,536)*y(k,23)
         mat(k,780) = -(rxt(k,537) + rxt(k,538) + rxt(k,539)) * y(k,23)
         mat(k,161) = -(rxt(k,503)*y(k,71) + rxt(k,520)*y(k,25) + rxt(k,521)*y(k,133))
         mat(k,774) = -rxt(k,503)*y(k,24)
         mat(k,1456) = -rxt(k,520)*y(k,24)
         mat(k,960) = -rxt(k,521)*y(k,24)
         mat(k,1487) = -((rxt(k,114) + rxt(k,115)) * y(k,94) + rxt(k,127)*y(k,97) &
                      + rxt(k,285)*y(k,111) + rxt(k,314)*y(k,127) + rxt(k,341) &
                      *y(k,134) + rxt(k,350)*y(k,135) + rxt(k,450)*y(k,17) + rxt(k,451) &
                      *y(k,23) + rxt(k,452)*y(k,41) + rxt(k,453)*y(k,43) + (rxt(k,454) &
                      + rxt(k,455)) * y(k,101) + rxt(k,456)*y(k,73) + rxt(k,463) &
                      *y(k,29) + rxt(k,472)*y(k,54) + rxt(k,513)*y(k,16) + rxt(k,515) &
                      *y(k,18) + rxt(k,518)*y(k,20) + rxt(k,520)*y(k,24))
         mat(k,1024) = -(rxt(k,114) + rxt(k,115)) * y(k,25)
         mat(k,848) = -rxt(k,127)*y(k,25)
         mat(k,936) = -rxt(k,285)*y(k,25)
         mat(k,1564) = -rxt(k,314)*y(k,25)
         mat(k,1852) = -rxt(k,341)*y(k,25)
         mat(k,1888) = -rxt(k,350)*y(k,25)
         mat(k,747) = -rxt(k,450)*y(k,25)
         mat(k,425) = -rxt(k,451)*y(k,25)
         mat(k,657) = -rxt(k,452)*y(k,25)
         mat(k,218) = -rxt(k,453)*y(k,25)
         mat(k,909) = -(rxt(k,454) + rxt(k,455)) * y(k,25)
         mat(k,1929) = -rxt(k,456)*y(k,25)
         mat(k,363) = -rxt(k,463)*y(k,25)
         mat(k,337) = -rxt(k,472)*y(k,25)
         mat(k,190) = -rxt(k,513)*y(k,25)
         mat(k,227) = -rxt(k,515)*y(k,25)
         mat(k,152) = -rxt(k,518)*y(k,25)
         mat(k,165) = -rxt(k,520)*y(k,25)
         mat(k,724) = rxt(k,491)*y(k,28)
         mat(k,31) = 4.000_r8*rxt(k,475)*y(k,71)
         mat(k,59) = rxt(k,476)*y(k,71)
         mat(k,39) = 2.000_r8*rxt(k,477)*y(k,71)
         mat(k,69) = 2.000_r8*rxt(k,478)*y(k,71)
         mat(k,43) = 2.000_r8*rxt(k,479)*y(k,71)
         mat(k,74) = rxt(k,480)*y(k,71)
         mat(k,47) = 2.000_r8*rxt(k,481)*y(k,71)
         mat(k,50) = 3.000_r8*rxt(k,517)*y(k,133)
         mat(k,152) = mat(k,152) + rxt(k,519)*y(k,133)
         mat(k,438) = rxt(k,457)*y(k,28)
         mat(k,1112) = rxt(k,491)*y(k,5) + rxt(k,457)*y(k,21) + (4.000_r8*rxt(k,458) &
                       +2.000_r8*rxt(k,460))*y(k,28) + rxt(k,462)*y(k,62) + rxt(k,467) &
                      *y(k,70) + rxt(k,549)*y(k,80) + rxt(k,468)*y(k,133)
         mat(k,88) = rxt(k,512)*y(k,71)
         mat(k,84) = rxt(k,527)*y(k,71) + rxt(k,522)*y(k,133)
         mat(k,94) = rxt(k,528)*y(k,71) + rxt(k,523)*y(k,133)
         mat(k,120) = rxt(k,529)*y(k,71) + rxt(k,524)*y(k,133)
         mat(k,2016) = rxt(k,470)*y(k,70) + rxt(k,482)*y(k,71) + rxt(k,471)*y(k,133)
         mat(k,1971) = rxt(k,462)*y(k,28) + rxt(k,111)*y(k,93)
         mat(k,1819) = rxt(k,110)*y(k,90)
         mat(k,1390) = rxt(k,467)*y(k,28) + rxt(k,470)*y(k,49)
         mat(k,796) = 4.000_r8*rxt(k,475)*y(k,8) + rxt(k,476)*y(k,9) &
                      + 2.000_r8*rxt(k,477)*y(k,11) + 2.000_r8*rxt(k,478)*y(k,12) &
                      + 2.000_r8*rxt(k,479)*y(k,13) + rxt(k,480)*y(k,14) &
                      + 2.000_r8*rxt(k,481)*y(k,15) + rxt(k,512)*y(k,34) + rxt(k,527) &
                      *y(k,46) + rxt(k,528)*y(k,47) + rxt(k,529)*y(k,48) + rxt(k,482) &
                      *y(k,49)
         mat(k,484) = rxt(k,549)*y(k,28)
         mat(k,820) = rxt(k,110)*y(k,63) + rxt(k,195)*y(k,103) + rxt(k,147)*y(k,105) &
                      + rxt(k,177)*y(k,107) + rxt(k,248)*y(k,118) + rxt(k,230) &
                      *y(k,119) + rxt(k,212)*y(k,122) + rxt(k,154)*y(k,128)
         mat(k,542) = rxt(k,199)*y(k,103) + rxt(k,164)*y(k,105) + rxt(k,182)*y(k,107) &
                      + rxt(k,252)*y(k,118) + rxt(k,234)*y(k,119) + rxt(k,217) &
                      *y(k,122) + rxt(k,159)*y(k,128)
         mat(k,526) = rxt(k,187)*y(k,103) + rxt(k,181)*y(k,105) + rxt(k,169)*y(k,107) &
                      + rxt(k,240)*y(k,118) + rxt(k,222)*y(k,119) + rxt(k,205) &
                      *y(k,122) + rxt(k,257)*y(k,128)
         mat(k,374) = rxt(k,111)*y(k,62)
         mat(k,1246) = rxt(k,195)*y(k,90) + rxt(k,199)*y(k,91) + rxt(k,187)*y(k,92)
         mat(k,1695) = rxt(k,147)*y(k,90) + rxt(k,164)*y(k,91) + rxt(k,181)*y(k,92)
         mat(k,1652) = rxt(k,177)*y(k,90) + rxt(k,182)*y(k,91) + rxt(k,169)*y(k,92)
         mat(k,1160) = rxt(k,248)*y(k,90) + rxt(k,252)*y(k,91) + rxt(k,240)*y(k,92)
         mat(k,1203) = rxt(k,230)*y(k,90) + rxt(k,234)*y(k,91) + rxt(k,222)*y(k,92)
         mat(k,1330) = rxt(k,212)*y(k,90) + rxt(k,217)*y(k,91) + rxt(k,205)*y(k,92)
         mat(k,1609) = rxt(k,154)*y(k,90) + rxt(k,159)*y(k,91) + rxt(k,257)*y(k,92)
         mat(k,993) = 3.000_r8*rxt(k,517)*y(k,19) + rxt(k,519)*y(k,20) + rxt(k,468) &
                      *y(k,28) + rxt(k,522)*y(k,46) + rxt(k,523)*y(k,47) + rxt(k,524) &
                      *y(k,48) + rxt(k,471)*y(k,49)
         mat(k,1454) = rxt(k,463)*y(k,29)
         mat(k,1083) = 2.000_r8*rxt(k,459)*y(k,28)
         mat(k,357) = rxt(k,463)*y(k,25) + (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,49)
         mat(k,1986) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,29) + (rxt(k,563) &
                       +rxt(k,569)+rxt(k,574))*y(k,54)
         mat(k,332) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,49)
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
         mat(k,1082) = 2.000_r8*rxt(k,484)*y(k,28)
         mat(k,1104) = -(rxt(k,116)*y(k,94) + rxt(k,128)*y(k,97) + rxt(k,286)*y(k,111) &
                      + rxt(k,315)*y(k,127) + rxt(k,342)*y(k,134) + rxt(k,351) &
                      *y(k,135) + rxt(k,457)*y(k,21) + (4._r8*rxt(k,458) &
                      + 4._r8*rxt(k,459) + 4._r8*rxt(k,460) + 4._r8*rxt(k,484) &
                      ) * y(k,28) + rxt(k,461)*y(k,101) + rxt(k,462)*y(k,62) + rxt(k,464) &
                      *y(k,63) + rxt(k,467)*y(k,70) + (rxt(k,468) + rxt(k,469) &
                      ) * y(k,133) + (rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,5) &
                      + rxt(k,549)*y(k,80))
         mat(k,1016) = -rxt(k,116)*y(k,28)
         mat(k,840) = -rxt(k,128)*y(k,28)
         mat(k,929) = -rxt(k,286)*y(k,28)
         mat(k,1556) = -rxt(k,315)*y(k,28)
         mat(k,1844) = -rxt(k,342)*y(k,28)
         mat(k,1880) = -rxt(k,351)*y(k,28)
         mat(k,436) = -rxt(k,457)*y(k,28)
         mat(k,903) = -rxt(k,461)*y(k,28)
         mat(k,1963) = -rxt(k,462)*y(k,28)
         mat(k,1811) = -rxt(k,464)*y(k,28)
         mat(k,1382) = -rxt(k,467)*y(k,28)
         mat(k,986) = -(rxt(k,468) + rxt(k,469)) * y(k,28)
         mat(k,720) = -(rxt(k,490) + rxt(k,491) + rxt(k,492)) * y(k,28)
         mat(k,481) = -rxt(k,549)*y(k,28)
         mat(k,1479) = rxt(k,472)*y(k,54) + rxt(k,456)*y(k,73) + rxt(k,455)*y(k,101)
         mat(k,361) = rxt(k,465)*y(k,70)
         mat(k,2008) = rxt(k,483)*y(k,71)
         mat(k,335) = rxt(k,472)*y(k,25) + rxt(k,473)*y(k,70) + rxt(k,474)*y(k,133)
         mat(k,1382) = mat(k,1382) + rxt(k,465)*y(k,29) + rxt(k,473)*y(k,54)
         mat(k,791) = rxt(k,483)*y(k,49)
         mat(k,1921) = rxt(k,456)*y(k,25)
         mat(k,137) = rxt(k,554)*y(k,80)
         mat(k,481) = mat(k,481) + rxt(k,554)*y(k,74)
         mat(k,903) = mat(k,903) + rxt(k,455)*y(k,25)
         mat(k,986) = mat(k,986) + rxt(k,474)*y(k,54)
         mat(k,359) = -(rxt(k,463)*y(k,25) + rxt(k,465)*y(k,70) + rxt(k,466)*y(k,133) &
                      + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,49))
         mat(k,1462) = -rxt(k,463)*y(k,29)
         mat(k,1361) = -rxt(k,465)*y(k,29)
         mat(k,973) = -rxt(k,466)*y(k,29)
         mat(k,1991) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,29)
         mat(k,1088) = rxt(k,464)*y(k,63)
         mat(k,1791) = rxt(k,464)*y(k,28)
         mat(k,252) = -(rxt(k,540)*y(k,133))
         mat(k,965) = -rxt(k,540)*y(k,31)
         mat(k,580) = rxt(k,486)*y(k,17)
         mat(k,733) = rxt(k,486)*y(k,3) + rxt(k,450)*y(k,25) + rxt(k,530)*y(k,64) &
                      + rxt(k,531)*y(k,70) + rxt(k,532)*y(k,133)
         mat(k,148) = rxt(k,518)*y(k,25)
         mat(k,1460) = rxt(k,450)*y(k,17) + rxt(k,518)*y(k,20)
         mat(k,1263) = rxt(k,594)*y(k,136)
         mat(k,1711) = rxt(k,530)*y(k,17)
         mat(k,1353) = rxt(k,531)*y(k,17) + rxt(k,543)*y(k,75)
         mat(k,140) = rxt(k,543)*y(k,70) + rxt(k,544)*y(k,133)
         mat(k,965) = mat(k,965) + rxt(k,532)*y(k,17) + rxt(k,544)*y(k,75)
         mat(k,385) = rxt(k,594)*y(k,32)
         mat(k,1284) = -(rxt(k,312)*y(k,123) + rxt(k,316)*y(k,127) + rxt(k,330) &
                      *y(k,130) + rxt(k,335)*y(k,131) + rxt(k,343)*y(k,134) + rxt(k,352) &
                      *y(k,135) + rxt(k,368)*y(k,118) + rxt(k,594)*y(k,136))
         mat(k,260) = -rxt(k,312)*y(k,32)
         mat(k,1560) = -rxt(k,316)*y(k,32)
         mat(k,599) = -rxt(k,330)*y(k,32)
         mat(k,242) = -rxt(k,335)*y(k,32)
         mat(k,1848) = -rxt(k,343)*y(k,32)
         mat(k,1884) = -rxt(k,352)*y(k,32)
         mat(k,1156) = -rxt(k,368)*y(k,32)
         mat(k,394) = -rxt(k,594)*y(k,32)
         mat(k,1483) = (rxt(k,114)+rxt(k,115))*y(k,94) + rxt(k,127)*y(k,97)
         mat(k,1108) = rxt(k,116)*y(k,94) + rxt(k,128)*y(k,97)
         mat(k,255) = rxt(k,540)*y(k,133)
         mat(k,1767) = rxt(k,117)*y(k,94)
         mat(k,2012) = rxt(k,130)*y(k,97)
         mat(k,1524) = rxt(k,124)*y(k,94)
         mat(k,1967) = rxt(k,280)*y(k,94) + (rxt(k,122)+rxt(k,123))*y(k,96)
         mat(k,1815) = rxt(k,281)*y(k,94) + (rxt(k,120)+rxt(k,121))*y(k,96)
         mat(k,1386) = rxt(k,125)*y(k,94)
         mat(k,1435) = rxt(k,126)*y(k,94)
         mat(k,1925) = rxt(k,132)*y(k,97)
         mat(k,1020) = (rxt(k,114)+rxt(k,115))*y(k,25) + rxt(k,116)*y(k,28) &
                      + rxt(k,117)*y(k,40) + rxt(k,124)*y(k,51) + rxt(k,280)*y(k,62) &
                      + rxt(k,281)*y(k,63) + rxt(k,125)*y(k,70) + rxt(k,126)*y(k,72) &
                      + rxt(k,186)*y(k,103) + (rxt(k,170)+rxt(k,258))*y(k,105) + ( &
                      + rxt(k,168)+rxt(k,265))*y(k,107) + rxt(k,239)*y(k,118) &
                      + rxt(k,221)*y(k,119) + rxt(k,204)*y(k,122) + rxt(k,256) &
                      *y(k,128)
         mat(k,509) = rxt(k,194)*y(k,103) + (rxt(k,247)+rxt(k,271))*y(k,105) + ( &
                      + rxt(k,176)+rxt(k,259))*y(k,107) + rxt(k,246)*y(k,118) &
                      + rxt(k,229)*y(k,119) + rxt(k,211)*y(k,122) + rxt(k,153) &
                      *y(k,128)
         mat(k,677) = (rxt(k,122)+rxt(k,123))*y(k,62) + (rxt(k,120)+rxt(k,121)) &
                      *y(k,63) + rxt(k,196)*y(k,103) + (rxt(k,158)+rxt(k,260)) &
                      *y(k,105) + (rxt(k,178)+rxt(k,261))*y(k,107) + rxt(k,249) &
                      *y(k,118) + rxt(k,231)*y(k,119) + rxt(k,213)*y(k,122) &
                      + rxt(k,155)*y(k,128)
         mat(k,844) = rxt(k,127)*y(k,25) + rxt(k,128)*y(k,28) + rxt(k,130)*y(k,49) &
                      + rxt(k,132)*y(k,73) + rxt(k,191)*y(k,103) + rxt(k,225)*y(k,105) &
                      + rxt(k,174)*y(k,107) + rxt(k,244)*y(k,118) + rxt(k,227) &
                      *y(k,119) + rxt(k,209)*y(k,122) + rxt(k,151)*y(k,128)
         mat(k,1063) = rxt(k,306)*y(k,121)
         mat(k,495) = rxt(k,189)*y(k,103) + rxt(k,203)*y(k,105) + rxt(k,172)*y(k,107) &
                      + rxt(k,242)*y(k,118) + rxt(k,224)*y(k,119) + rxt(k,207) &
                      *y(k,122) + rxt(k,149)*y(k,128)
         mat(k,1242) = rxt(k,186)*y(k,94) + rxt(k,194)*y(k,95) + rxt(k,196)*y(k,96) &
                      + rxt(k,191)*y(k,97) + rxt(k,189)*y(k,100)
         mat(k,1691) = (rxt(k,170)+rxt(k,258))*y(k,94) + (rxt(k,247)+rxt(k,271)) &
                      *y(k,95) + (rxt(k,158)+rxt(k,260))*y(k,96) + rxt(k,225)*y(k,97) &
                      + rxt(k,203)*y(k,100)
         mat(k,1648) = (rxt(k,168)+rxt(k,265))*y(k,94) + (rxt(k,176)+rxt(k,259)) &
                      *y(k,95) + (rxt(k,178)+rxt(k,261))*y(k,96) + rxt(k,174)*y(k,97) &
                      + rxt(k,172)*y(k,100)
         mat(k,1156) = mat(k,1156) + rxt(k,239)*y(k,94) + rxt(k,246)*y(k,95) &
                      + rxt(k,249)*y(k,96) + rxt(k,244)*y(k,97) + rxt(k,242)*y(k,100)
         mat(k,1199) = rxt(k,221)*y(k,94) + rxt(k,229)*y(k,95) + rxt(k,231)*y(k,96) &
                      + rxt(k,227)*y(k,97) + rxt(k,224)*y(k,100)
         mat(k,267) = rxt(k,306)*y(k,98) + rxt(k,307)*y(k,139)
         mat(k,1326) = rxt(k,204)*y(k,94) + rxt(k,211)*y(k,95) + rxt(k,213)*y(k,96) &
                      + rxt(k,209)*y(k,97) + rxt(k,207)*y(k,100)
         mat(k,1605) = rxt(k,256)*y(k,94) + rxt(k,153)*y(k,95) + rxt(k,155)*y(k,96) &
                      + rxt(k,151)*y(k,97) + rxt(k,149)*y(k,100)
         mat(k,989) = rxt(k,540)*y(k,31)
         mat(k,2072) = rxt(k,307)*y(k,121)
         mat(k,77) = -(rxt(k,511)*y(k,71))
         mat(k,767) = -rxt(k,511)*y(k,33)
         mat(k,56) = rxt(k,476)*y(k,71)
         mat(k,61) = rxt(k,502)*y(k,71)
         mat(k,66) = rxt(k,478)*y(k,71)
         mat(k,41) = 2.000_r8*rxt(k,479)*y(k,71)
         mat(k,71) = 2.000_r8*rxt(k,480)*y(k,71)
         mat(k,45) = rxt(k,481)*y(k,71)
         mat(k,33) = 2.000_r8*rxt(k,504)*y(k,71)
         mat(k,89) = rxt(k,528)*y(k,71) + rxt(k,523)*y(k,133)
         mat(k,115) = rxt(k,529)*y(k,71) + rxt(k,524)*y(k,133)
         mat(k,767) = mat(k,767) + rxt(k,476)*y(k,9) + rxt(k,502)*y(k,10) + rxt(k,478) &
                      *y(k,12) + 2.000_r8*rxt(k,479)*y(k,13) + 2.000_r8*rxt(k,480) &
                      *y(k,14) + rxt(k,481)*y(k,15) + 2.000_r8*rxt(k,504)*y(k,42) &
                      + rxt(k,528)*y(k,47) + rxt(k,529)*y(k,48)
         mat(k,950) = rxt(k,523)*y(k,47) + rxt(k,524)*y(k,48)
         mat(k,85) = -(rxt(k,512)*y(k,71))
         mat(k,769) = -rxt(k,512)*y(k,34)
         mat(k,37) = rxt(k,477)*y(k,71)
         mat(k,67) = rxt(k,478)*y(k,71)
         mat(k,81) = rxt(k,527)*y(k,71) + rxt(k,522)*y(k,133)
         mat(k,769) = mat(k,769) + rxt(k,477)*y(k,11) + rxt(k,478)*y(k,12) &
                      + rxt(k,527)*y(k,46)
         mat(k,952) = rxt(k,522)*y(k,46)
         mat(k,125) = -(rxt(k,541)*y(k,64) + (rxt(k,542) + rxt(k,556)) * y(k,133))
         mat(k,1710) = -rxt(k,541)*y(k,35)
         mat(k,956) = -(rxt(k,542) + rxt(k,556)) * y(k,35)
         mat(k,341) = -(rxt(k,507)*y(k,23) + rxt(k,508)*y(k,41) + rxt(k,509)*y(k,139) &
                      + rxt(k,510)*y(k,51))
         mat(k,415) = -rxt(k,507)*y(k,39)
         mat(k,644) = -rxt(k,508)*y(k,39)
         mat(k,2044) = -rxt(k,509)*y(k,39)
         mat(k,1504) = -rxt(k,510)*y(k,39)
         mat(k,62) = rxt(k,502)*y(k,71)
         mat(k,72) = rxt(k,480)*y(k,71)
         mat(k,78) = 2.000_r8*rxt(k,511)*y(k,71)
         mat(k,86) = rxt(k,512)*y(k,71)
         mat(k,778) = rxt(k,502)*y(k,10) + rxt(k,480)*y(k,14) + 2.000_r8*rxt(k,511) &
                      *y(k,33) + rxt(k,512)*y(k,34)
         mat(k,1778) = -(rxt(k,105)*y(k,90) + rxt(k,117)*y(k,94) + rxt(k,129)*y(k,97) &
                      + rxt(k,287)*y(k,111) + rxt(k,309)*y(k,122) + rxt(k,317) &
                      *y(k,127) + rxt(k,331)*y(k,130) + rxt(k,344)*y(k,134) + (rxt(k,408) &
                      + rxt(k,409) + rxt(k,410)) * y(k,101) + rxt(k,411)*y(k,72) &
                      + rxt(k,414)*y(k,73))
         mat(k,825) = -rxt(k,105)*y(k,40)
         mat(k,1031) = -rxt(k,117)*y(k,40)
         mat(k,854) = -rxt(k,129)*y(k,40)
         mat(k,942) = -rxt(k,287)*y(k,40)
         mat(k,1337) = -rxt(k,309)*y(k,40)
         mat(k,1571) = -rxt(k,317)*y(k,40)
         mat(k,603) = -rxt(k,331)*y(k,40)
         mat(k,1859) = -rxt(k,344)*y(k,40)
         mat(k,914) = -(rxt(k,408) + rxt(k,409) + rxt(k,410)) * y(k,40)
         mat(k,1446) = -rxt(k,411)*y(k,40)
         mat(k,1936) = -rxt(k,414)*y(k,40)
         mat(k,750) = rxt(k,532)*y(k,133)
         mat(k,428) = rxt(k,538)*y(k,71)
         mat(k,1494) = rxt(k,452)*y(k,41)
         mat(k,348) = rxt(k,508)*y(k,41)
         mat(k,660) = rxt(k,452)*y(k,25) + rxt(k,508)*y(k,39) + rxt(k,406)*y(k,70) &
                      + rxt(k,389)*y(k,71) + rxt(k,415)*y(k,133) + rxt(k,354)*y(k,135)
         mat(k,330) = rxt(k,506)*y(k,71)
         mat(k,2023) = rxt(k,483)*y(k,71)
         mat(k,456) = rxt(k,438)*y(k,133)
         mat(k,1397) = rxt(k,406)*y(k,41) + rxt(k,418)*y(k,133)
         mat(k,801) = rxt(k,538)*y(k,23) + rxt(k,389)*y(k,41) + rxt(k,506)*y(k,45) &
                      + rxt(k,483)*y(k,49)
         mat(k,146) = rxt(k,544)*y(k,133)
         mat(k,276) = rxt(k,550)*y(k,133)
         mat(k,485) = rxt(k,555)*y(k,133)
         mat(k,825) = mat(k,825) + rxt(k,195)*y(k,103) + rxt(k,147)*y(k,105) &
                      + rxt(k,177)*y(k,107)
         mat(k,546) = rxt(k,199)*y(k,103) + rxt(k,164)*y(k,105) + rxt(k,182)*y(k,107)
         mat(k,530) = rxt(k,187)*y(k,103) + rxt(k,181)*y(k,105) + rxt(k,169)*y(k,107)
         mat(k,1031) = mat(k,1031) + rxt(k,186)*y(k,103) + (rxt(k,170)+rxt(k,258)) &
                      *y(k,105) + (rxt(k,168)+rxt(k,265))*y(k,107)
         mat(k,516) = rxt(k,194)*y(k,103) + (rxt(k,247)+rxt(k,271))*y(k,105) + ( &
                      + rxt(k,176)+rxt(k,259))*y(k,107)
         mat(k,684) = rxt(k,196)*y(k,103) + (rxt(k,158)+rxt(k,260))*y(k,105) + ( &
                      + rxt(k,178)+rxt(k,261))*y(k,107)
         mat(k,854) = mat(k,854) + rxt(k,191)*y(k,103) + rxt(k,225)*y(k,105) &
                      + rxt(k,174)*y(k,107)
         mat(k,1074) = rxt(k,138)*y(k,99) + rxt(k,382)*y(k,102) + rxt(k,383)*y(k,103) &
                      + rxt(k,141)*y(k,105) + rxt(k,144)*y(k,107) + rxt(k,381) &
                      *y(k,108)
         mat(k,171) = rxt(k,138)*y(k,98)
         mat(k,501) = rxt(k,189)*y(k,103) + rxt(k,203)*y(k,105) + rxt(k,172)*y(k,107)
         mat(k,250) = rxt(k,382)*y(k,98)
         mat(k,1253) = rxt(k,195)*y(k,90) + rxt(k,199)*y(k,91) + rxt(k,187)*y(k,92) &
                      + rxt(k,186)*y(k,94) + rxt(k,194)*y(k,95) + rxt(k,196)*y(k,96) &
                      + rxt(k,191)*y(k,97) + rxt(k,383)*y(k,98) + rxt(k,189)*y(k,100) &
                      + rxt(k,201)*y(k,111) + rxt(k,197)*y(k,112) + rxt(k,200) &
                      *y(k,114) + rxt(k,193)*y(k,115) + rxt(k,198)*y(k,116) &
                      + rxt(k,190)*y(k,127)
         mat(k,1702) = rxt(k,147)*y(k,90) + rxt(k,164)*y(k,91) + rxt(k,181)*y(k,92) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,94) + (rxt(k,247)+rxt(k,271)) &
                      *y(k,95) + (rxt(k,158)+rxt(k,260))*y(k,96) + rxt(k,225)*y(k,97) &
                      + rxt(k,141)*y(k,98) + rxt(k,203)*y(k,100) + rxt(k,166)*y(k,111) &
                      + rxt(k,162)*y(k,112) + rxt(k,165)*y(k,114) + (rxt(k,236) &
                       +rxt(k,262))*y(k,115) + rxt(k,163)*y(k,116) + rxt(k,214) &
                      *y(k,127)
         mat(k,1659) = rxt(k,177)*y(k,90) + rxt(k,182)*y(k,91) + rxt(k,169)*y(k,92) + ( &
                      + rxt(k,168)+rxt(k,265))*y(k,94) + (rxt(k,176)+rxt(k,259)) &
                      *y(k,95) + (rxt(k,178)+rxt(k,261))*y(k,96) + rxt(k,174)*y(k,97) &
                      + rxt(k,144)*y(k,98) + rxt(k,172)*y(k,100) + rxt(k,184)*y(k,111) &
                      + rxt(k,179)*y(k,112) + rxt(k,183)*y(k,114) + (rxt(k,175) &
                       +rxt(k,263))*y(k,115) + rxt(k,180)*y(k,116) + rxt(k,173) &
                      *y(k,127)
         mat(k,287) = rxt(k,381)*y(k,98)
         mat(k,942) = mat(k,942) + rxt(k,201)*y(k,103) + rxt(k,166)*y(k,105) &
                      + rxt(k,184)*y(k,107)
         mat(k,560) = rxt(k,197)*y(k,103) + rxt(k,162)*y(k,105) + rxt(k,179)*y(k,107)
         mat(k,639) = rxt(k,200)*y(k,103) + rxt(k,165)*y(k,105) + rxt(k,183)*y(k,107)
         mat(k,704) = rxt(k,193)*y(k,103) + (rxt(k,236)+rxt(k,262))*y(k,105) + ( &
                      + rxt(k,175)+rxt(k,263))*y(k,107)
         mat(k,576) = rxt(k,198)*y(k,103) + rxt(k,163)*y(k,105) + rxt(k,180)*y(k,107)
         mat(k,1571) = mat(k,1571) + rxt(k,190)*y(k,103) + rxt(k,214)*y(k,105) &
                      + rxt(k,173)*y(k,107)
         mat(k,998) = rxt(k,532)*y(k,17) + rxt(k,415)*y(k,41) + rxt(k,438)*y(k,56) &
                      + rxt(k,418)*y(k,70) + rxt(k,544)*y(k,75) + rxt(k,550)*y(k,78) &
                      + rxt(k,555)*y(k,80)
         mat(k,1895) = rxt(k,354)*y(k,41)
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
         mat(k,647) = -((rxt(k,353) + rxt(k,354)) * y(k,135) + rxt(k,389)*y(k,71) &
                      + rxt(k,406)*y(k,70) + rxt(k,415)*y(k,133) + rxt(k,452)*y(k,25) &
                      + rxt(k,508)*y(k,39))
         mat(k,1870) = -(rxt(k,353) + rxt(k,354)) * y(k,41)
         mat(k,784) = -rxt(k,389)*y(k,41)
         mat(k,1370) = -rxt(k,406)*y(k,41)
         mat(k,979) = -rxt(k,415)*y(k,41)
         mat(k,1467) = -rxt(k,452)*y(k,41)
         mat(k,344) = -rxt(k,508)*y(k,41)
         mat(k,418) = rxt(k,539)*y(k,71)
         mat(k,1752) = rxt(k,408)*y(k,101)
         mat(k,784) = mat(k,784) + rxt(k,539)*y(k,23)
         mat(k,896) = rxt(k,408)*y(k,40)
         mat(k,32) = -(rxt(k,504)*y(k,71))
         mat(k,759) = -rxt(k,504)*y(k,42)
         mat(k,214) = -(rxt(k,407)*y(k,70) + rxt(k,416)*y(k,133) + rxt(k,453)*y(k,25))
         mat(k,1350) = -rxt(k,407)*y(k,43)
         mat(k,963) = -rxt(k,416)*y(k,43)
         mat(k,1458) = -rxt(k,453)*y(k,43)
         mat(k,889) = 2.000_r8*rxt(k,422)*y(k,101)
         mat(k,963) = mat(k,963) + 2.000_r8*rxt(k,421)*y(k,133)
         mat(k,102) = rxt(k,557)*y(k,139)
         mat(k,2031) = rxt(k,557)*y(k,82)
         mat(k,324) = -(rxt(k,499)*y(k,70) + rxt(k,500)*y(k,133) + (rxt(k,505) &
                      + rxt(k,506)) * y(k,71))
         mat(k,1358) = -rxt(k,499)*y(k,45)
         mat(k,969) = -rxt(k,500)*y(k,45)
         mat(k,777) = -(rxt(k,505) + rxt(k,506)) * y(k,45)
         mat(k,581) = rxt(k,486)*y(k,17) + rxt(k,487)*y(k,101)
         mat(k,734) = rxt(k,486)*y(k,3)
         mat(k,892) = rxt(k,487)*y(k,3)
         mat(k,80) = -(rxt(k,522)*y(k,133) + rxt(k,527)*y(k,71))
         mat(k,951) = -rxt(k,522)*y(k,46)
         mat(k,768) = -rxt(k,527)*y(k,46)
         mat(k,90) = -(rxt(k,523)*y(k,133) + rxt(k,528)*y(k,71))
         mat(k,953) = -rxt(k,523)*y(k,47)
         mat(k,770) = -rxt(k,528)*y(k,47)
         mat(k,116) = -(rxt(k,524)*y(k,133) + rxt(k,529)*y(k,71))
         mat(k,955) = -rxt(k,524)*y(k,48)
         mat(k,772) = -rxt(k,529)*y(k,48)
         mat(k,2029) = -(rxt(k,106)*y(k,91) + rxt(k,108)*y(k,90) + rxt(k,130)*y(k,97) &
                      + (rxt(k,276) + rxt(k,298)) * y(k,113) + rxt(k,289)*y(k,111) &
                      + rxt(k,318)*y(k,127) + rxt(k,345)*y(k,134) + rxt(k,356) &
                      *y(k,135) + rxt(k,470)*y(k,70) + rxt(k,471)*y(k,133) + (rxt(k,482) &
                      + rxt(k,483)) * y(k,71) + (rxt(k,563) + rxt(k,569) + rxt(k,574) &
                      ) * y(k,54) + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,29) &
                      + (rxt(k,570) + rxt(k,575)) * y(k,53))
         mat(k,548) = -rxt(k,106)*y(k,49)
         mat(k,828) = -rxt(k,108)*y(k,49)
         mat(k,860) = -rxt(k,130)*y(k,49)
         mat(k,885) = -(rxt(k,276) + rxt(k,298)) * y(k,49)
         mat(k,947) = -rxt(k,289)*y(k,49)
         mat(k,1577) = -rxt(k,318)*y(k,49)
         mat(k,1865) = -rxt(k,345)*y(k,49)
         mat(k,1901) = -rxt(k,356)*y(k,49)
         mat(k,1403) = -rxt(k,470)*y(k,49)
         mat(k,1004) = -rxt(k,471)*y(k,49)
         mat(k,807) = -(rxt(k,482) + rxt(k,483)) * y(k,49)
         mat(k,338) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,49)
         mat(k,367) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,49)
         mat(k,315) = -(rxt(k,570) + rxt(k,575)) * y(k,49)
         mat(k,191) = rxt(k,513)*y(k,25)
         mat(k,756) = rxt(k,450)*y(k,25)
         mat(k,228) = rxt(k,515)*y(k,25)
         mat(k,153) = 2.000_r8*rxt(k,518)*y(k,25)
         mat(k,429) = rxt(k,451)*y(k,25)
         mat(k,166) = rxt(k,520)*y(k,25)
         mat(k,1500) = rxt(k,513)*y(k,16) + rxt(k,450)*y(k,17) + rxt(k,515)*y(k,18) &
                      + 2.000_r8*rxt(k,518)*y(k,20) + rxt(k,451)*y(k,23) + rxt(k,520) &
                      *y(k,24) + rxt(k,452)*y(k,41) + rxt(k,453)*y(k,43) + rxt(k,472) &
                      *y(k,54) + rxt(k,454)*y(k,101)
         mat(k,1125) = rxt(k,469)*y(k,133)
         mat(k,1784) = rxt(k,105)*y(k,90)
         mat(k,665) = rxt(k,452)*y(k,25)
         mat(k,219) = rxt(k,453)*y(k,25)
         mat(k,1541) = rxt(k,109)*y(k,90) + rxt(k,299)*y(k,116)
         mat(k,338) = mat(k,338) + rxt(k,472)*y(k,25)
         mat(k,828) = mat(k,828) + rxt(k,105)*y(k,40) + rxt(k,109)*y(k,51)
         mat(k,532) = rxt(k,187)*y(k,103) + (rxt(k,181)+2.000_r8*rxt(k,267))*y(k,105) + ( &
                      + rxt(k,169)+2.000_r8*rxt(k,268))*y(k,107) + rxt(k,240)*y(k,118) &
                      + rxt(k,222)*y(k,119) + rxt(k,205)*y(k,122) + rxt(k,257) &
                      *y(k,128)
         mat(k,920) = rxt(k,454)*y(k,25)
         mat(k,1259) = rxt(k,187)*y(k,92) + rxt(k,198)*y(k,116)
         mat(k,1708) = (rxt(k,181)+2.000_r8*rxt(k,267))*y(k,92) + rxt(k,163)*y(k,116)
         mat(k,1665) = (rxt(k,169)+2.000_r8*rxt(k,268))*y(k,92) + rxt(k,180)*y(k,116)
         mat(k,578) = rxt(k,299)*y(k,51) + rxt(k,198)*y(k,103) + rxt(k,163)*y(k,105) &
                      + rxt(k,180)*y(k,107) + rxt(k,251)*y(k,118) + rxt(k,233) &
                      *y(k,119) + rxt(k,216)*y(k,122) + rxt(k,157)*y(k,128)
         mat(k,1173) = rxt(k,240)*y(k,92) + rxt(k,251)*y(k,116)
         mat(k,1216) = rxt(k,222)*y(k,92) + rxt(k,233)*y(k,116)
         mat(k,1343) = rxt(k,205)*y(k,92) + rxt(k,216)*y(k,116)
         mat(k,1622) = rxt(k,257)*y(k,92) + rxt(k,157)*y(k,116)
         mat(k,1004) = mat(k,1004) + rxt(k,469)*y(k,28)
         mat(k,413) = rxt(k,507)*y(k,39)
         mat(k,340) = rxt(k,507)*y(k,23) + rxt(k,508)*y(k,41) + rxt(k,510)*y(k,51) &
                      + rxt(k,509)*y(k,139)
         mat(k,643) = rxt(k,508)*y(k,39)
         mat(k,1502) = rxt(k,510)*y(k,39)
         mat(k,2035) = rxt(k,509)*y(k,39)
         mat(k,1529) = -(rxt(k,109)*y(k,90) + rxt(k,124)*y(k,94) + rxt(k,290)*y(k,111) &
                      + rxt(k,295)*y(k,115) + rxt(k,299)*y(k,116) + rxt(k,300) &
                      *y(k,113) + rxt(k,319)*y(k,127) + rxt(k,357)*y(k,135) + rxt(k,447) &
                      *y(k,133) + rxt(k,510)*y(k,39))
         mat(k,821) = -rxt(k,109)*y(k,51)
         mat(k,1025) = -rxt(k,124)*y(k,51)
         mat(k,937) = -rxt(k,290)*y(k,51)
         mat(k,699) = -rxt(k,295)*y(k,51)
         mat(k,571) = -rxt(k,299)*y(k,51)
         mat(k,876) = -rxt(k,300)*y(k,51)
         mat(k,1565) = -rxt(k,319)*y(k,51)
         mat(k,1889) = -rxt(k,357)*y(k,51)
         mat(k,994) = -rxt(k,447)*y(k,51)
         mat(k,346) = -rxt(k,510)*y(k,51)
         mat(k,748) = rxt(k,530)*y(k,64)
         mat(k,364) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,49)
         mat(k,129) = rxt(k,541)*y(k,64)
         mat(k,2017) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,29) + rxt(k,298) &
                      *y(k,113)
         mat(k,406) = rxt(k,142)*y(k,105) + rxt(k,145)*y(k,107) + rxt(k,293)*y(k,114) &
                      + rxt(k,297)*y(k,115)
         mat(k,1820) = rxt(k,446)*y(k,133)
         mat(k,1736) = rxt(k,530)*y(k,17) + rxt(k,541)*y(k,35)
         mat(k,1247) = rxt(k,188)*y(k,113) + 2.000_r8*rxt(k,185)*y(k,117)
         mat(k,108) = rxt(k,140)*y(k,139)
         mat(k,1696) = rxt(k,142)*y(k,58) + (rxt(k,192)+rxt(k,264))*y(k,113) + ( &
                      + 2.000_r8*rxt(k,146)+2.000_r8*rxt(k,269))*y(k,117)
         mat(k,112) = rxt(k,143)*y(k,139)
         mat(k,1653) = rxt(k,145)*y(k,58) + (rxt(k,171)+rxt(k,266))*y(k,113) + ( &
                      + 2.000_r8*rxt(k,167)+2.000_r8*rxt(k,270))*y(k,117)
         mat(k,876) = mat(k,876) + rxt(k,298)*y(k,49) + rxt(k,188)*y(k,103) + ( &
                      + rxt(k,192)+rxt(k,264))*y(k,105) + (rxt(k,171)+rxt(k,266)) &
                      *y(k,107)
         mat(k,634) = rxt(k,293)*y(k,58)
         mat(k,699) = mat(k,699) + rxt(k,297)*y(k,58)
         mat(k,616) = 2.000_r8*rxt(k,185)*y(k,103) + (2.000_r8*rxt(k,146) &
                       +2.000_r8*rxt(k,269))*y(k,105) + (2.000_r8*rxt(k,167) &
                       +2.000_r8*rxt(k,270))*y(k,107) + rxt(k,238)*y(k,118) &
                      + rxt(k,220)*y(k,119) + rxt(k,202)*y(k,122) + rxt(k,255) &
                      *y(k,128)
         mat(k,1161) = rxt(k,238)*y(k,117)
         mat(k,1204) = rxt(k,220)*y(k,117)
         mat(k,1331) = rxt(k,202)*y(k,117)
         mat(k,1610) = rxt(k,255)*y(k,117)
         mat(k,994) = mat(k,994) + rxt(k,446)*y(k,63)
         mat(k,2077) = rxt(k,140)*y(k,104) + rxt(k,143)*y(k,106)
         mat(k,179) = -(rxt(k,423)*y(k,133))
         mat(k,961) = -rxt(k,423)*y(k,52)
         mat(k,1787) = rxt(k,444)*y(k,101)
         mat(k,888) = rxt(k,444)*y(k,63)
         mat(k,309) = -(rxt(k,501)*y(k,70) + (rxt(k,570) + rxt(k,575)) * y(k,49))
         mat(k,1356) = -rxt(k,501)*y(k,53)
         mat(k,1989) = -(rxt(k,570) + rxt(k,575)) * y(k,53)
         mat(k,711) = rxt(k,493)*y(k,101)
         mat(k,891) = rxt(k,493)*y(k,5)
         mat(k,333) = -(rxt(k,472)*y(k,25) + rxt(k,473)*y(k,70) + rxt(k,474)*y(k,133) &
                      + (rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,49))
         mat(k,1461) = -rxt(k,472)*y(k,54)
         mat(k,1359) = -rxt(k,473)*y(k,54)
         mat(k,970) = -rxt(k,474)*y(k,54)
         mat(k,1990) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,54)
         mat(k,1086) = rxt(k,461)*y(k,101)
         mat(k,358) = rxt(k,466)*y(k,133)
         mat(k,893) = rxt(k,461)*y(k,28)
         mat(k,970) = mat(k,970) + rxt(k,466)*y(k,29)
         mat(k,278) = -(rxt(k,340)*y(k,133))
         mat(k,967) = -rxt(k,340)*y(k,55)
         mat(k,1988) = rxt(k,289)*y(k,111)
         mat(k,1503) = rxt(k,290)*y(k,111)
         mat(k,1944) = rxt(k,349)*y(k,133)
         mat(k,922) = rxt(k,289)*y(k,49) + rxt(k,290)*y(k,51)
         mat(k,174) = rxt(k,305)*y(k,139)
         mat(k,967) = mat(k,967) + rxt(k,349)*y(k,62)
         mat(k,2041) = rxt(k,305)*y(k,120)
         mat(k,446) = -(rxt(k,426)*y(k,62) + (rxt(k,427) + rxt(k,428) + rxt(k,429) &
                      ) * y(k,63) + rxt(k,430)*y(k,72) + rxt(k,438)*y(k,133) + rxt(k,591) &
                      *y(k,128))
         mat(k,1947) = -rxt(k,426)*y(k,56)
         mat(k,1793) = -(rxt(k,427) + rxt(k,428) + rxt(k,429)) * y(k,56)
         mat(k,1416) = -rxt(k,430)*y(k,56)
         mat(k,976) = -rxt(k,438)*y(k,56)
         mat(k,1580) = -rxt(k,591)*y(k,56)
         mat(k,1365) = rxt(k,424)*y(k,109) + rxt(k,588)*y(k,124)
         mat(k,1416) = mat(k,1416) + rxt(k,589)*y(k,124)
         mat(k,1053) = 1.100_r8*rxt(k,584)*y(k,110) + .200_r8*rxt(k,582)*y(k,118)
         mat(k,303) = rxt(k,424)*y(k,70)
         mat(k,233) = 1.100_r8*rxt(k,584)*y(k,98)
         mat(k,1132) = .200_r8*rxt(k,582)*y(k,98)
         mat(k,319) = rxt(k,588)*y(k,70) + rxt(k,589)*y(k,72)
         mat(k,98) = -((rxt(k,442) + rxt(k,443)) * y(k,71))
         mat(k,771) = -(rxt(k,442) + rxt(k,443)) * y(k,57)
         mat(k,443) = rxt(k,427)*y(k,63)
         mat(k,1786) = rxt(k,427)*y(k,56)
         mat(k,401) = -(rxt(k,142)*y(k,105) + rxt(k,145)*y(k,107) + rxt(k,293) &
                      *y(k,114) + rxt(k,297)*y(k,115))
         mat(k,1668) = -rxt(k,142)*y(k,58)
         mat(k,1625) = -rxt(k,145)*y(k,58)
         mat(k,623) = -rxt(k,293)*y(k,58)
         mat(k,688) = -rxt(k,297)*y(k,58)
         mat(k,1792) = rxt(k,445)*y(k,64)
         mat(k,1713) = rxt(k,445)*y(k,63)
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
         mat(k,1983) = -((rxt(k,111) + rxt(k,112)) * y(k,93) + (rxt(k,122) + rxt(k,123) &
                      ) * y(k,96) + rxt(k,136)*y(k,135) + (rxt(k,272) + rxt(k,279) &
                      ) * y(k,130) + rxt(k,280)*y(k,94) + rxt(k,349)*y(k,133) &
                      + rxt(k,426)*y(k,56) + rxt(k,435)*y(k,64) + rxt(k,439)*y(k,101) &
                      + rxt(k,440)*y(k,73) + rxt(k,441)*y(k,70) + rxt(k,462)*y(k,28) &
                      + rxt(k,494)*y(k,5) + rxt(k,534)*y(k,21) + rxt(k,593)*y(k,128))
         mat(k,376) = -(rxt(k,111) + rxt(k,112)) * y(k,62)
         mat(k,686) = -(rxt(k,122) + rxt(k,123)) * y(k,62)
         mat(k,1900) = -rxt(k,136)*y(k,62)
         mat(k,608) = -(rxt(k,272) + rxt(k,279)) * y(k,62)
         mat(k,1036) = -rxt(k,280)*y(k,62)
         mat(k,1003) = -rxt(k,349)*y(k,62)
         mat(k,458) = -rxt(k,426)*y(k,62)
         mat(k,1747) = -rxt(k,435)*y(k,62)
         mat(k,919) = -rxt(k,439)*y(k,62)
         mat(k,1941) = -rxt(k,440)*y(k,62)
         mat(k,1402) = -rxt(k,441)*y(k,62)
         mat(k,1124) = -rxt(k,462)*y(k,62)
         mat(k,730) = -rxt(k,494)*y(k,62)
         mat(k,441) = -rxt(k,534)*y(k,62)
         mat(k,1621) = -rxt(k,593)*y(k,62)
         mat(k,1783) = rxt(k,287)*y(k,111) + rxt(k,309)*y(k,122)
         mat(k,458) = mat(k,458) + 2.000_r8*rxt(k,428)*y(k,63) + rxt(k,430)*y(k,72) &
                      + rxt(k,438)*y(k,133)
         mat(k,101) = 2.000_r8*rxt(k,442)*y(k,71)
         mat(k,1831) = 2.000_r8*rxt(k,428)*y(k,56) + rxt(k,431)*y(k,70) + rxt(k,551) &
                      *y(k,80) + rxt(k,291)*y(k,111)
         mat(k,1402) = mat(k,1402) + rxt(k,431)*y(k,63)
         mat(k,806) = 2.000_r8*rxt(k,442)*y(k,57)
         mat(k,1451) = rxt(k,430)*y(k,56) + rxt(k,425)*y(k,109)
         mat(k,488) = rxt(k,551)*y(k,63)
         mat(k,827) = rxt(k,248)*y(k,118) + rxt(k,230)*y(k,119) + rxt(k,212)*y(k,122)
         mat(k,547) = rxt(k,252)*y(k,118) + rxt(k,234)*y(k,119) + rxt(k,217)*y(k,122)
         mat(k,531) = rxt(k,240)*y(k,118) + rxt(k,222)*y(k,119) + rxt(k,205)*y(k,122)
         mat(k,1036) = mat(k,1036) + rxt(k,239)*y(k,118) + rxt(k,221)*y(k,119) &
                      + rxt(k,204)*y(k,122)
         mat(k,517) = rxt(k,246)*y(k,118) + rxt(k,229)*y(k,119) + rxt(k,211)*y(k,122)
         mat(k,686) = mat(k,686) + rxt(k,249)*y(k,118) + rxt(k,231)*y(k,119) &
                      + rxt(k,213)*y(k,122)
         mat(k,859) = rxt(k,244)*y(k,118) + rxt(k,227)*y(k,119) + rxt(k,209)*y(k,122)
         mat(k,1079) = rxt(k,303)*y(k,119) + rxt(k,304)*y(k,120) + rxt(k,306)*y(k,121) &
                      + rxt(k,308)*y(k,122) + rxt(k,384)*y(k,123)
         mat(k,502) = rxt(k,242)*y(k,118) + rxt(k,224)*y(k,119) + rxt(k,207)*y(k,122)
         mat(k,307) = rxt(k,425)*y(k,72)
         mat(k,946) = rxt(k,287)*y(k,40) + rxt(k,291)*y(k,63) + rxt(k,254)*y(k,118) &
                      + rxt(k,237)*y(k,119) + rxt(k,219)*y(k,122)
         mat(k,562) = rxt(k,250)*y(k,118) + rxt(k,232)*y(k,119) + rxt(k,215)*y(k,122)
         mat(k,884) = rxt(k,241)*y(k,118) + rxt(k,223)*y(k,119) + rxt(k,206)*y(k,122)
         mat(k,641) = rxt(k,253)*y(k,118) + rxt(k,235)*y(k,119) + rxt(k,218)*y(k,122)
         mat(k,706) = rxt(k,245)*y(k,118) + rxt(k,228)*y(k,119) + rxt(k,210)*y(k,122)
         mat(k,577) = rxt(k,251)*y(k,118) + rxt(k,233)*y(k,119) + rxt(k,216)*y(k,122)
         mat(k,621) = rxt(k,238)*y(k,118) + rxt(k,220)*y(k,119) + rxt(k,202)*y(k,122)
         mat(k,1172) = rxt(k,248)*y(k,90) + rxt(k,252)*y(k,91) + rxt(k,240)*y(k,92) &
                      + rxt(k,239)*y(k,94) + rxt(k,246)*y(k,95) + rxt(k,249)*y(k,96) &
                      + rxt(k,244)*y(k,97) + rxt(k,242)*y(k,100) + rxt(k,254)*y(k,111) &
                      + rxt(k,250)*y(k,112) + rxt(k,241)*y(k,113) + rxt(k,253) &
                      *y(k,114) + rxt(k,245)*y(k,115) + rxt(k,251)*y(k,116) &
                      + rxt(k,238)*y(k,117) + rxt(k,243)*y(k,127)
         mat(k,1215) = rxt(k,230)*y(k,90) + rxt(k,234)*y(k,91) + rxt(k,222)*y(k,92) &
                      + rxt(k,221)*y(k,94) + rxt(k,229)*y(k,95) + rxt(k,231)*y(k,96) &
                      + rxt(k,227)*y(k,97) + rxt(k,303)*y(k,98) + rxt(k,224)*y(k,100) &
                      + rxt(k,237)*y(k,111) + rxt(k,232)*y(k,112) + rxt(k,223) &
                      *y(k,113) + rxt(k,235)*y(k,114) + rxt(k,228)*y(k,115) &
                      + rxt(k,233)*y(k,116) + rxt(k,220)*y(k,117) + rxt(k,226) &
                      *y(k,127)
         mat(k,177) = rxt(k,304)*y(k,98)
         mat(k,269) = rxt(k,306)*y(k,98)
         mat(k,1342) = rxt(k,309)*y(k,40) + rxt(k,212)*y(k,90) + rxt(k,217)*y(k,91) &
                      + rxt(k,205)*y(k,92) + rxt(k,204)*y(k,94) + rxt(k,211)*y(k,95) &
                      + rxt(k,213)*y(k,96) + rxt(k,209)*y(k,97) + rxt(k,308)*y(k,98) &
                      + rxt(k,207)*y(k,100) + rxt(k,219)*y(k,111) + rxt(k,215) &
                      *y(k,112) + rxt(k,206)*y(k,113) + rxt(k,218)*y(k,114) &
                      + rxt(k,210)*y(k,115) + rxt(k,216)*y(k,116) + rxt(k,202) &
                      *y(k,117) + rxt(k,208)*y(k,127)
         mat(k,262) = rxt(k,384)*y(k,98)
         mat(k,1576) = rxt(k,243)*y(k,118) + rxt(k,226)*y(k,119) + rxt(k,208)*y(k,122)
         mat(k,1003) = mat(k,1003) + rxt(k,438)*y(k,56)
         mat(k,1827) = -(rxt(k,110)*y(k,90) + (rxt(k,120) + rxt(k,121)) * y(k,96) &
                      + (rxt(k,277) + rxt(k,278)) * y(k,130) + rxt(k,281)*y(k,94) &
                      + rxt(k,291)*y(k,111) + rxt(k,320)*y(k,127) + rxt(k,346) &
                      *y(k,134) + rxt(k,359)*y(k,135) + (rxt(k,427) + rxt(k,428) &
                      + rxt(k,429)) * y(k,56) + (rxt(k,431) + rxt(k,433)) * y(k,70) &
                      + rxt(k,432)*y(k,73) + rxt(k,444)*y(k,101) + rxt(k,445)*y(k,64) &
                      + rxt(k,446)*y(k,133) + rxt(k,464)*y(k,28) + rxt(k,495)*y(k,5) &
                      + rxt(k,551)*y(k,80))
         mat(k,826) = -rxt(k,110)*y(k,63)
         mat(k,685) = -(rxt(k,120) + rxt(k,121)) * y(k,63)
         mat(k,604) = -(rxt(k,277) + rxt(k,278)) * y(k,63)
         mat(k,1032) = -rxt(k,281)*y(k,63)
         mat(k,943) = -rxt(k,291)*y(k,63)
         mat(k,1572) = -rxt(k,320)*y(k,63)
         mat(k,1860) = -rxt(k,346)*y(k,63)
         mat(k,1896) = -rxt(k,359)*y(k,63)
         mat(k,457) = -(rxt(k,427) + rxt(k,428) + rxt(k,429)) * y(k,63)
         mat(k,1398) = -(rxt(k,431) + rxt(k,433)) * y(k,63)
         mat(k,1937) = -rxt(k,432)*y(k,63)
         mat(k,915) = -rxt(k,444)*y(k,63)
         mat(k,1743) = -rxt(k,445)*y(k,63)
         mat(k,999) = -rxt(k,446)*y(k,63)
         mat(k,1120) = -rxt(k,464)*y(k,63)
         mat(k,728) = -rxt(k,495)*y(k,63)
         mat(k,486) = -rxt(k,551)*y(k,63)
         mat(k,728) = mat(k,728) + rxt(k,494)*y(k,62)
         mat(k,440) = rxt(k,534)*y(k,62)
         mat(k,1495) = rxt(k,285)*y(k,111)
         mat(k,1120) = mat(k,1120) + rxt(k,462)*y(k,62)
         mat(k,184) = rxt(k,423)*y(k,133)
         mat(k,280) = rxt(k,340)*y(k,133)
         mat(k,1979) = rxt(k,494)*y(k,5) + rxt(k,534)*y(k,21) + rxt(k,462)*y(k,28) &
                      + 2.000_r8*rxt(k,435)*y(k,64) + rxt(k,441)*y(k,70) + rxt(k,440) &
                      *y(k,73) + rxt(k,112)*y(k,93) + rxt(k,439)*y(k,101) + rxt(k,136) &
                      *y(k,135)
         mat(k,1743) = mat(k,1743) + 2.000_r8*rxt(k,435)*y(k,62) + rxt(k,436)*y(k,70) &
                      + rxt(k,434)*y(k,101) + rxt(k,437)*y(k,133)
         mat(k,1398) = mat(k,1398) + rxt(k,441)*y(k,62) + rxt(k,436)*y(k,64)
         mat(k,1937) = mat(k,1937) + rxt(k,440)*y(k,62)
         mat(k,375) = rxt(k,112)*y(k,62)
         mat(k,915) = mat(k,915) + rxt(k,439)*y(k,62) + rxt(k,434)*y(k,64)
         mat(k,1254) = rxt(k,201)*y(k,111) + rxt(k,197)*y(k,112)
         mat(k,1703) = rxt(k,166)*y(k,111) + rxt(k,162)*y(k,112)
         mat(k,1660) = rxt(k,184)*y(k,111) + rxt(k,179)*y(k,112)
         mat(k,943) = mat(k,943) + rxt(k,285)*y(k,25) + rxt(k,201)*y(k,103) &
                      + rxt(k,166)*y(k,105) + rxt(k,184)*y(k,107) + rxt(k,254) &
                      *y(k,118) + rxt(k,237)*y(k,119) + rxt(k,219)*y(k,122) &
                      + rxt(k,161)*y(k,128)
         mat(k,561) = rxt(k,197)*y(k,103) + rxt(k,162)*y(k,105) + rxt(k,179)*y(k,107) &
                      + rxt(k,250)*y(k,118) + rxt(k,232)*y(k,119) + rxt(k,215) &
                      *y(k,122) + rxt(k,156)*y(k,128)
         mat(k,1168) = rxt(k,254)*y(k,111) + rxt(k,250)*y(k,112)
         mat(k,1211) = rxt(k,237)*y(k,111) + rxt(k,232)*y(k,112)
         mat(k,1338) = rxt(k,219)*y(k,111) + rxt(k,215)*y(k,112) + rxt(k,311)*y(k,133)
         mat(k,1617) = rxt(k,161)*y(k,111) + rxt(k,156)*y(k,112)
         mat(k,999) = mat(k,999) + rxt(k,423)*y(k,52) + rxt(k,340)*y(k,55) &
                      + rxt(k,437)*y(k,64) + rxt(k,311)*y(k,122)
         mat(k,1896) = mat(k,1896) + rxt(k,136)*y(k,62)
         mat(k,1741) = -(rxt(k,434)*y(k,101) + rxt(k,435)*y(k,62) + rxt(k,436)*y(k,70) &
                      + rxt(k,437)*y(k,133) + rxt(k,445)*y(k,63) + rxt(k,530)*y(k,17) &
                      + rxt(k,541)*y(k,35))
         mat(k,913) = -rxt(k,434)*y(k,64)
         mat(k,1977) = -rxt(k,435)*y(k,64)
         mat(k,1396) = -rxt(k,436)*y(k,64)
         mat(k,997) = -rxt(k,437)*y(k,64)
         mat(k,1825) = -rxt(k,445)*y(k,64)
         mat(k,749) = -rxt(k,530)*y(k,64)
         mat(k,130) = -rxt(k,541)*y(k,64)
         mat(k,212) = rxt(k,496)*y(k,70)
         mat(k,1493) = rxt(k,463)*y(k,29)
         mat(k,1118) = rxt(k,286)*y(k,111)
         mat(k,365) = rxt(k,463)*y(k,25) + rxt(k,465)*y(k,70) + rxt(k,466)*y(k,133)
         mat(k,347) = rxt(k,510)*y(k,51)
         mat(k,1534) = rxt(k,510)*y(k,39) + rxt(k,447)*y(k,133)
         mat(k,1825) = mat(k,1825) + rxt(k,433)*y(k,70) + rxt(k,432)*y(k,73)
         mat(k,1396) = mat(k,1396) + rxt(k,496)*y(k,6) + rxt(k,465)*y(k,29) &
                      + rxt(k,433)*y(k,63)
         mat(k,1935) = rxt(k,432)*y(k,63)
         mat(k,913) = mat(k,913) + rxt(k,310)*y(k,122)
         mat(k,1252) = rxt(k,200)*y(k,114) + rxt(k,193)*y(k,115) + rxt(k,198)*y(k,116)
         mat(k,1701) = rxt(k,165)*y(k,114) + (rxt(k,236)+rxt(k,262))*y(k,115) &
                      + rxt(k,163)*y(k,116)
         mat(k,1658) = rxt(k,183)*y(k,114) + (rxt(k,175)+rxt(k,263))*y(k,115) &
                      + rxt(k,180)*y(k,116)
         mat(k,941) = rxt(k,286)*y(k,28)
         mat(k,880) = rxt(k,241)*y(k,118) + rxt(k,223)*y(k,119) + rxt(k,206)*y(k,122) &
                      + rxt(k,148)*y(k,128)
         mat(k,638) = rxt(k,200)*y(k,103) + rxt(k,165)*y(k,105) + rxt(k,183)*y(k,107) &
                      + rxt(k,253)*y(k,118) + rxt(k,235)*y(k,119) + rxt(k,218) &
                      *y(k,122) + rxt(k,160)*y(k,128)
         mat(k,703) = rxt(k,193)*y(k,103) + (rxt(k,236)+rxt(k,262))*y(k,105) + ( &
                      + rxt(k,175)+rxt(k,263))*y(k,107) + rxt(k,245)*y(k,118) &
                      + rxt(k,228)*y(k,119) + rxt(k,210)*y(k,122) + rxt(k,152) &
                      *y(k,128)
         mat(k,575) = rxt(k,198)*y(k,103) + rxt(k,163)*y(k,105) + rxt(k,180)*y(k,107) &
                      + rxt(k,251)*y(k,118) + rxt(k,233)*y(k,119) + rxt(k,216) &
                      *y(k,122) + rxt(k,157)*y(k,128)
         mat(k,620) = rxt(k,238)*y(k,118) + rxt(k,220)*y(k,119) + rxt(k,202)*y(k,122) &
                      + rxt(k,255)*y(k,128)
         mat(k,1166) = rxt(k,241)*y(k,113) + rxt(k,253)*y(k,114) + rxt(k,245)*y(k,115) &
                      + rxt(k,251)*y(k,116) + rxt(k,238)*y(k,117)
         mat(k,1209) = rxt(k,223)*y(k,113) + rxt(k,235)*y(k,114) + rxt(k,228)*y(k,115) &
                      + rxt(k,233)*y(k,116) + rxt(k,220)*y(k,117)
         mat(k,1336) = rxt(k,310)*y(k,101) + rxt(k,206)*y(k,113) + rxt(k,218)*y(k,114) &
                      + rxt(k,210)*y(k,115) + rxt(k,216)*y(k,116) + rxt(k,202) &
                      *y(k,117)
         mat(k,1615) = rxt(k,148)*y(k,113) + rxt(k,160)*y(k,114) + rxt(k,152)*y(k,115) &
                      + rxt(k,157)*y(k,116) + rxt(k,255)*y(k,117)
         mat(k,997) = mat(k,997) + rxt(k,466)*y(k,29) + rxt(k,447)*y(k,51)
         mat(k,1388) = -(rxt(k,113)*y(k,93) + rxt(k,125)*y(k,94) + rxt(k,131)*y(k,97) &
                      + rxt(k,301)*y(k,113) + (rxt(k,324) + rxt(k,325)) * y(k,127) &
                      + (rxt(k,333) + rxt(k,334)) * y(k,130) + rxt(k,336)*y(k,131) &
                      + rxt(k,338)*y(k,132) + rxt(k,347)*y(k,134) + rxt(k,360) &
                      *y(k,135) + rxt(k,403)*y(k,73) + 4._r8*rxt(k,404)*y(k,70) &
                      + rxt(k,405)*y(k,72) + rxt(k,406)*y(k,41) + rxt(k,407)*y(k,43) &
                      + rxt(k,412)*y(k,101) + rxt(k,418)*y(k,133) + (rxt(k,431) &
                      + rxt(k,433)) * y(k,63) + rxt(k,436)*y(k,64) + rxt(k,441) &
                      *y(k,62) + rxt(k,465)*y(k,29) + rxt(k,467)*y(k,28) + rxt(k,470) &
                      *y(k,49) + rxt(k,473)*y(k,54) + rxt(k,496)*y(k,6) + rxt(k,497) &
                      *y(k,5) + rxt(k,499)*y(k,45) + rxt(k,501)*y(k,53) + rxt(k,531) &
                      *y(k,17) + rxt(k,543)*y(k,75) + (rxt(k,586) + rxt(k,587) &
                      ) * y(k,110) + rxt(k,588)*y(k,124))
         mat(k,372) = -rxt(k,113)*y(k,70)
         mat(k,1022) = -rxt(k,125)*y(k,70)
         mat(k,846) = -rxt(k,131)*y(k,70)
         mat(k,873) = -rxt(k,301)*y(k,70)
         mat(k,1562) = -(rxt(k,324) + rxt(k,325)) * y(k,70)
         mat(k,600) = -(rxt(k,333) + rxt(k,334)) * y(k,70)
         mat(k,243) = -rxt(k,336)*y(k,70)
         mat(k,465) = -rxt(k,338)*y(k,70)
         mat(k,1850) = -rxt(k,347)*y(k,70)
         mat(k,1886) = -rxt(k,360)*y(k,70)
         mat(k,1927) = -rxt(k,403)*y(k,70)
         mat(k,1437) = -rxt(k,405)*y(k,70)
         mat(k,655) = -rxt(k,406)*y(k,70)
         mat(k,217) = -rxt(k,407)*y(k,70)
         mat(k,907) = -rxt(k,412)*y(k,70)
         mat(k,991) = -rxt(k,418)*y(k,70)
         mat(k,1817) = -(rxt(k,431) + rxt(k,433)) * y(k,70)
         mat(k,1733) = -rxt(k,436)*y(k,70)
         mat(k,1969) = -rxt(k,441)*y(k,70)
         mat(k,362) = -rxt(k,465)*y(k,70)
         mat(k,1110) = -rxt(k,467)*y(k,70)
         mat(k,2014) = -rxt(k,470)*y(k,70)
         mat(k,336) = -rxt(k,473)*y(k,70)
         mat(k,210) = -rxt(k,496)*y(k,70)
         mat(k,722) = -rxt(k,497)*y(k,70)
         mat(k,329) = -rxt(k,499)*y(k,70)
         mat(k,313) = -rxt(k,501)*y(k,70)
         mat(k,745) = -rxt(k,531)*y(k,70)
         mat(k,145) = -rxt(k,543)*y(k,70)
         mat(k,236) = -(rxt(k,586) + rxt(k,587)) * y(k,70)
         mat(k,321) = -rxt(k,588)*y(k,70)
         mat(k,1485) = rxt(k,114)*y(k,94) + rxt(k,350)*y(k,135)
         mat(k,1769) = rxt(k,410)*y(k,101)
         mat(k,453) = rxt(k,426)*y(k,62) + rxt(k,427)*y(k,63) + rxt(k,430)*y(k,72) &
                      + rxt(k,591)*y(k,128)
         mat(k,1969) = mat(k,1969) + rxt(k,426)*y(k,56) + rxt(k,272)*y(k,130)
         mat(k,1817) = mat(k,1817) + rxt(k,427)*y(k,56) + rxt(k,359)*y(k,135)
         mat(k,794) = (rxt(k,392)+rxt(k,393))*y(k,72)
         mat(k,1437) = mat(k,1437) + rxt(k,430)*y(k,56) + (rxt(k,392)+rxt(k,393)) &
                      *y(k,71) + rxt(k,545)*y(k,78) + rxt(k,552)*y(k,80) + rxt(k,590) &
                      *y(k,124) + rxt(k,597)*y(k,136) + rxt(k,601)*y(k,137)
         mat(k,1927) = mat(k,1927) + rxt(k,363)*y(k,135)
         mat(k,274) = rxt(k,545)*y(k,72)
         mat(k,482) = rxt(k,552)*y(k,72)
         mat(k,1022) = mat(k,1022) + rxt(k,114)*y(k,25) + rxt(k,186)*y(k,103) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,105) + (rxt(k,168)+rxt(k,265)) &
                      *y(k,107) + rxt(k,239)*y(k,118) + rxt(k,221)*y(k,119) &
                      + rxt(k,204)*y(k,122) + rxt(k,256)*y(k,128)
         mat(k,511) = rxt(k,194)*y(k,103) + (rxt(k,247)+rxt(k,271))*y(k,105) + ( &
                      + rxt(k,176)+rxt(k,259))*y(k,107) + rxt(k,246)*y(k,118) &
                      + rxt(k,229)*y(k,119) + rxt(k,211)*y(k,122) + rxt(k,153) &
                      *y(k,128)
         mat(k,679) = rxt(k,196)*y(k,103) + (rxt(k,158)+rxt(k,260))*y(k,105) + ( &
                      + rxt(k,178)+rxt(k,261))*y(k,107) + rxt(k,249)*y(k,118) &
                      + rxt(k,231)*y(k,119) + rxt(k,213)*y(k,122) + rxt(k,155) &
                      *y(k,128)
         mat(k,1065) = rxt(k,582)*y(k,118) + 1.150_r8*rxt(k,583)*y(k,128)
         mat(k,907) = mat(k,907) + rxt(k,410)*y(k,40)
         mat(k,1244) = rxt(k,186)*y(k,94) + rxt(k,194)*y(k,95) + rxt(k,196)*y(k,96)
         mat(k,1693) = (rxt(k,170)+rxt(k,258))*y(k,94) + (rxt(k,247)+rxt(k,271)) &
                      *y(k,95) + (rxt(k,158)+rxt(k,260))*y(k,96)
         mat(k,1650) = (rxt(k,168)+rxt(k,265))*y(k,94) + (rxt(k,176)+rxt(k,259)) &
                      *y(k,95) + (rxt(k,178)+rxt(k,261))*y(k,96)
         mat(k,305) = rxt(k,596)*y(k,136)
         mat(k,1158) = rxt(k,239)*y(k,94) + rxt(k,246)*y(k,95) + rxt(k,249)*y(k,96) &
                      + rxt(k,582)*y(k,98)
         mat(k,1201) = rxt(k,221)*y(k,94) + rxt(k,229)*y(k,95) + rxt(k,231)*y(k,96)
         mat(k,1328) = rxt(k,204)*y(k,94) + rxt(k,211)*y(k,95) + rxt(k,213)*y(k,96)
         mat(k,321) = mat(k,321) + rxt(k,590)*y(k,72)
         mat(k,1607) = rxt(k,591)*y(k,56) + rxt(k,256)*y(k,94) + rxt(k,153)*y(k,95) &
                      + rxt(k,155)*y(k,96) + 1.150_r8*rxt(k,583)*y(k,98)
         mat(k,600) = mat(k,600) + rxt(k,272)*y(k,62)
         mat(k,991) = mat(k,991) + 2.000_r8*rxt(k,420)*y(k,133)
         mat(k,1886) = mat(k,1886) + rxt(k,350)*y(k,25) + rxt(k,359)*y(k,63) &
                      + rxt(k,363)*y(k,73)
         mat(k,395) = rxt(k,597)*y(k,72) + rxt(k,596)*y(k,109)
         mat(k,203) = rxt(k,601)*y(k,72)
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
         mat(k,787) = -(rxt(k,389)*y(k,41) + rxt(k,390)*y(k,139) + (rxt(k,392) &
                      + rxt(k,393)) * y(k,72) + rxt(k,394)*y(k,73) + (rxt(k,442) &
                      + rxt(k,443)) * y(k,57) + rxt(k,475)*y(k,8) + rxt(k,476)*y(k,9) &
                      + rxt(k,477)*y(k,11) + rxt(k,478)*y(k,12) + rxt(k,479)*y(k,13) &
                      + rxt(k,480)*y(k,14) + rxt(k,481)*y(k,15) + (rxt(k,482) &
                      + rxt(k,483)) * y(k,49) + rxt(k,502)*y(k,10) + rxt(k,503) &
                      *y(k,24) + rxt(k,504)*y(k,42) + (rxt(k,505) + rxt(k,506) &
                      ) * y(k,45) + rxt(k,511)*y(k,33) + rxt(k,512)*y(k,34) + rxt(k,525) &
                      *y(k,16) + rxt(k,526)*y(k,18) + rxt(k,527)*y(k,46) + rxt(k,528) &
                      *y(k,47) + rxt(k,529)*y(k,48) + (rxt(k,537) + rxt(k,538) &
                      + rxt(k,539)) * y(k,23))
         mat(k,649) = -rxt(k,389)*y(k,71)
         mat(k,2060) = -rxt(k,390)*y(k,71)
         mat(k,1423) = -(rxt(k,392) + rxt(k,393)) * y(k,71)
         mat(k,1913) = -rxt(k,394)*y(k,71)
         mat(k,99) = -(rxt(k,442) + rxt(k,443)) * y(k,71)
         mat(k,30) = -rxt(k,475)*y(k,71)
         mat(k,58) = -rxt(k,476)*y(k,71)
         mat(k,38) = -rxt(k,477)*y(k,71)
         mat(k,68) = -rxt(k,478)*y(k,71)
         mat(k,42) = -rxt(k,479)*y(k,71)
         mat(k,73) = -rxt(k,480)*y(k,71)
         mat(k,46) = -rxt(k,481)*y(k,71)
         mat(k,1999) = -(rxt(k,482) + rxt(k,483)) * y(k,71)
         mat(k,64) = -rxt(k,502)*y(k,71)
         mat(k,163) = -rxt(k,503)*y(k,71)
         mat(k,35) = -rxt(k,504)*y(k,71)
         mat(k,327) = -(rxt(k,505) + rxt(k,506)) * y(k,71)
         mat(k,79) = -rxt(k,511)*y(k,71)
         mat(k,87) = -rxt(k,512)*y(k,71)
         mat(k,188) = -rxt(k,525)*y(k,71)
         mat(k,224) = -rxt(k,526)*y(k,71)
         mat(k,82) = -rxt(k,527)*y(k,71)
         mat(k,92) = -rxt(k,528)*y(k,71)
         mat(k,118) = -rxt(k,529)*y(k,71)
         mat(k,420) = -(rxt(k,537) + rxt(k,538) + rxt(k,539)) * y(k,71)
         mat(k,1423) = mat(k,1423) + rxt(k,425)*y(k,109)
         mat(k,1055) = .850_r8*rxt(k,583)*y(k,128)
         mat(k,304) = rxt(k,425)*y(k,72)
         mat(k,1592) = .850_r8*rxt(k,583)*y(k,98)
         mat(k,1438) = -(rxt(k,126)*y(k,94) + (rxt(k,133) + rxt(k,135)) * y(k,98) &
                      + rxt(k,322)*y(k,127) + rxt(k,362)*y(k,135) + rxt(k,364) &
                      *y(k,128) + rxt(k,392)*y(k,71) + rxt(k,397)*y(k,125) + rxt(k,405) &
                      *y(k,70) + rxt(k,411)*y(k,40) + rxt(k,425)*y(k,109) + rxt(k,430) &
                      *y(k,56) + rxt(k,545)*y(k,78) + rxt(k,552)*y(k,80) + rxt(k,585) &
                      *y(k,110) + (rxt(k,589) + rxt(k,590)) * y(k,124) + rxt(k,597) &
                      *y(k,136) + rxt(k,601)*y(k,137))
         mat(k,1023) = -rxt(k,126)*y(k,72)
         mat(k,1066) = -(rxt(k,133) + rxt(k,135)) * y(k,72)
         mat(k,1563) = -rxt(k,322)*y(k,72)
         mat(k,1887) = -rxt(k,362)*y(k,72)
         mat(k,1608) = -rxt(k,364)*y(k,72)
         mat(k,795) = -rxt(k,392)*y(k,72)
         mat(k,380) = -rxt(k,397)*y(k,72)
         mat(k,1389) = -rxt(k,405)*y(k,72)
         mat(k,1770) = -rxt(k,411)*y(k,72)
         mat(k,306) = -rxt(k,425)*y(k,72)
         mat(k,454) = -rxt(k,430)*y(k,72)
         mat(k,275) = -rxt(k,545)*y(k,72)
         mat(k,483) = -rxt(k,552)*y(k,72)
         mat(k,237) = -rxt(k,585)*y(k,72)
         mat(k,322) = -(rxt(k,589) + rxt(k,590)) * y(k,72)
         mat(k,396) = -rxt(k,597)*y(k,72)
         mat(k,204) = -rxt(k,601)*y(k,72)
         mat(k,590) = rxt(k,488)*y(k,73) + rxt(k,487)*y(k,101)
         mat(k,723) = 2.000_r8*rxt(k,489)*y(k,5) + (rxt(k,491)+rxt(k,492))*y(k,28) &
                      + rxt(k,497)*y(k,70) + rxt(k,493)*y(k,101)
         mat(k,437) = rxt(k,533)*y(k,101)
         mat(k,1486) = rxt(k,456)*y(k,73) + rxt(k,127)*y(k,97) + rxt(k,454)*y(k,101) &
                      + rxt(k,314)*y(k,127)
         mat(k,1111) = (rxt(k,491)+rxt(k,492))*y(k,5) + (2.000_r8*rxt(k,458) &
                       +2.000_r8*rxt(k,459))*y(k,28) + rxt(k,467)*y(k,70) + rxt(k,116) &
                      *y(k,94) + rxt(k,128)*y(k,97) + rxt(k,461)*y(k,101) + rxt(k,315) &
                      *y(k,127) + rxt(k,469)*y(k,133) + rxt(k,351)*y(k,135)
         mat(k,1287) = rxt(k,330)*y(k,130) + rxt(k,335)*y(k,131)
         mat(k,1770) = mat(k,1770) + rxt(k,414)*y(k,73) + rxt(k,408)*y(k,101) &
                      + rxt(k,331)*y(k,130)
         mat(k,182) = rxt(k,423)*y(k,133)
         mat(k,454) = mat(k,454) + rxt(k,429)*y(k,63)
         mat(k,100) = rxt(k,443)*y(k,71)
         mat(k,1970) = rxt(k,440)*y(k,73) + rxt(k,593)*y(k,128) + rxt(k,279)*y(k,130)
         mat(k,1818) = rxt(k,429)*y(k,56) + rxt(k,431)*y(k,70) + rxt(k,432)*y(k,73) &
                      + rxt(k,320)*y(k,127) + rxt(k,277)*y(k,130)
         mat(k,1734) = rxt(k,436)*y(k,70) + rxt(k,434)*y(k,101)
         mat(k,1389) = mat(k,1389) + rxt(k,497)*y(k,5) + rxt(k,467)*y(k,28) &
                      + rxt(k,431)*y(k,63) + rxt(k,436)*y(k,64) + 2.000_r8*rxt(k,404) &
                      *y(k,70) + 2.000_r8*rxt(k,403)*y(k,73) + rxt(k,113)*y(k,93) &
                      + rxt(k,131)*y(k,97) + rxt(k,412)*y(k,101) + rxt(k,301)*y(k,113) &
                      + rxt(k,396)*y(k,125) + rxt(k,325)*y(k,127) + ( &
                      + 2.000_r8*rxt(k,333)+rxt(k,334))*y(k,130) + rxt(k,336)*y(k,131) &
                      + rxt(k,418)*y(k,133) + rxt(k,360)*y(k,135)
         mat(k,795) = mat(k,795) + rxt(k,443)*y(k,57) + 2.000_r8*rxt(k,394)*y(k,73)
         mat(k,1438) = mat(k,1438) + 2.000_r8*rxt(k,397)*y(k,125)
         mat(k,1928) = rxt(k,488)*y(k,3) + rxt(k,456)*y(k,25) + rxt(k,414)*y(k,40) &
                      + rxt(k,440)*y(k,62) + rxt(k,432)*y(k,63) + 2.000_r8*rxt(k,403) &
                      *y(k,70) + 2.000_r8*rxt(k,394)*y(k,71) + rxt(k,547)*y(k,78) &
                      + rxt(k,553)*y(k,80) + rxt(k,132)*y(k,97) + rxt(k,134)*y(k,98) &
                      + 2.000_r8*rxt(k,413)*y(k,101) + rxt(k,292)*y(k,111) &
                      + 2.000_r8*rxt(k,302)*y(k,113) + rxt(k,323)*y(k,127) &
                      + 3.000_r8*rxt(k,332)*y(k,130) + rxt(k,419)*y(k,133)
         mat(k,275) = mat(k,275) + rxt(k,547)*y(k,73)
         mat(k,483) = mat(k,483) + rxt(k,553)*y(k,73)
         mat(k,819) = rxt(k,154)*y(k,128)
         mat(k,541) = rxt(k,159)*y(k,128)
         mat(k,525) = rxt(k,257)*y(k,128)
         mat(k,373) = rxt(k,113)*y(k,70)
         mat(k,1023) = mat(k,1023) + rxt(k,116)*y(k,28) + rxt(k,256)*y(k,128)
         mat(k,512) = rxt(k,153)*y(k,128)
         mat(k,680) = rxt(k,155)*y(k,128)
         mat(k,847) = rxt(k,127)*y(k,25) + rxt(k,128)*y(k,28) + rxt(k,131)*y(k,70) &
                      + rxt(k,132)*y(k,73) + rxt(k,191)*y(k,103) + rxt(k,225)*y(k,105) &
                      + rxt(k,174)*y(k,107) + rxt(k,244)*y(k,118) + rxt(k,227) &
                      *y(k,119) + rxt(k,209)*y(k,122) + 2.000_r8*rxt(k,151)*y(k,128)
         mat(k,1066) = mat(k,1066) + rxt(k,134)*y(k,73) + rxt(k,326)*y(k,129) &
                      + 2.000_r8*rxt(k,380)*y(k,132)
         mat(k,497) = rxt(k,149)*y(k,128)
         mat(k,908) = rxt(k,487)*y(k,3) + rxt(k,493)*y(k,5) + rxt(k,533)*y(k,21) &
                      + rxt(k,454)*y(k,25) + rxt(k,461)*y(k,28) + rxt(k,408)*y(k,40) &
                      + rxt(k,434)*y(k,64) + rxt(k,412)*y(k,70) + 2.000_r8*rxt(k,413) &
                      *y(k,73) + 2.000_r8*rxt(k,422)*y(k,101) + rxt(k,417)*y(k,133)
         mat(k,1245) = rxt(k,191)*y(k,97) + rxt(k,190)*y(k,127)
         mat(k,1694) = rxt(k,225)*y(k,97) + rxt(k,214)*y(k,127)
         mat(k,1651) = rxt(k,174)*y(k,97) + rxt(k,173)*y(k,127)
         mat(k,935) = rxt(k,292)*y(k,73) + rxt(k,161)*y(k,128)
         mat(k,556) = rxt(k,156)*y(k,128)
         mat(k,874) = rxt(k,301)*y(k,70) + 2.000_r8*rxt(k,302)*y(k,73) + rxt(k,148) &
                      *y(k,128)
         mat(k,633) = rxt(k,160)*y(k,128)
         mat(k,698) = rxt(k,152)*y(k,128)
         mat(k,570) = rxt(k,157)*y(k,128)
         mat(k,615) = rxt(k,255)*y(k,128)
         mat(k,1159) = rxt(k,244)*y(k,97) + rxt(k,243)*y(k,127)
         mat(k,1202) = rxt(k,227)*y(k,97) + rxt(k,226)*y(k,127)
         mat(k,1329) = rxt(k,209)*y(k,97) + rxt(k,208)*y(k,127)
         mat(k,380) = mat(k,380) + rxt(k,396)*y(k,70) + 2.000_r8*rxt(k,397)*y(k,72) &
                      + 2.000_r8*rxt(k,321)*y(k,127) + 2.000_r8*rxt(k,339)*y(k,132)
         mat(k,1563) = mat(k,1563) + rxt(k,314)*y(k,25) + rxt(k,315)*y(k,28) &
                      + rxt(k,320)*y(k,63) + rxt(k,325)*y(k,70) + rxt(k,323)*y(k,73) &
                      + rxt(k,190)*y(k,103) + rxt(k,214)*y(k,105) + rxt(k,173) &
                      *y(k,107) + rxt(k,243)*y(k,118) + rxt(k,226)*y(k,119) &
                      + rxt(k,208)*y(k,122) + 2.000_r8*rxt(k,321)*y(k,125)
         mat(k,1608) = mat(k,1608) + rxt(k,593)*y(k,62) + rxt(k,154)*y(k,90) &
                      + rxt(k,159)*y(k,91) + rxt(k,257)*y(k,92) + rxt(k,256)*y(k,94) &
                      + rxt(k,153)*y(k,95) + rxt(k,155)*y(k,96) + 2.000_r8*rxt(k,151) &
                      *y(k,97) + rxt(k,149)*y(k,100) + rxt(k,161)*y(k,111) &
                      + rxt(k,156)*y(k,112) + rxt(k,148)*y(k,113) + rxt(k,160) &
                      *y(k,114) + rxt(k,152)*y(k,115) + rxt(k,157)*y(k,116) &
                      + rxt(k,255)*y(k,117)
         mat(k,296) = rxt(k,326)*y(k,98) + (rxt(k,327)+rxt(k,328))*y(k,139)
         mat(k,601) = rxt(k,330)*y(k,32) + rxt(k,331)*y(k,40) + rxt(k,279)*y(k,62) &
                      + rxt(k,277)*y(k,63) + (2.000_r8*rxt(k,333)+rxt(k,334))*y(k,70) &
                      + 3.000_r8*rxt(k,332)*y(k,73)
         mat(k,244) = rxt(k,335)*y(k,32) + rxt(k,336)*y(k,70)
         mat(k,466) = 2.000_r8*rxt(k,380)*y(k,98) + 2.000_r8*rxt(k,339)*y(k,125) &
                      + rxt(k,337)*y(k,139)
         mat(k,992) = rxt(k,469)*y(k,28) + rxt(k,423)*y(k,52) + rxt(k,418)*y(k,70) &
                      + rxt(k,419)*y(k,73) + rxt(k,417)*y(k,101)
         mat(k,1887) = mat(k,1887) + rxt(k,351)*y(k,28) + rxt(k,360)*y(k,70)
         mat(k,2075) = (rxt(k,327)+rxt(k,328))*y(k,129) + rxt(k,337)*y(k,132)
         mat(k,1940) = -(rxt(k,132)*y(k,97) + rxt(k,134)*y(k,98) + rxt(k,292)*y(k,111) &
                      + rxt(k,302)*y(k,113) + rxt(k,323)*y(k,127) + rxt(k,332) &
                      *y(k,130) + rxt(k,348)*y(k,134) + rxt(k,363)*y(k,135) + rxt(k,394) &
                      *y(k,71) + rxt(k,403)*y(k,70) + rxt(k,413)*y(k,101) + rxt(k,414) &
                      *y(k,40) + rxt(k,419)*y(k,133) + rxt(k,432)*y(k,63) + rxt(k,440) &
                      *y(k,62) + rxt(k,456)*y(k,25) + rxt(k,488)*y(k,3) + rxt(k,547) &
                      *y(k,78) + rxt(k,553)*y(k,80))
         mat(k,858) = -rxt(k,132)*y(k,73)
         mat(k,1078) = -rxt(k,134)*y(k,73)
         mat(k,945) = -rxt(k,292)*y(k,73)
         mat(k,883) = -rxt(k,302)*y(k,73)
         mat(k,1575) = -rxt(k,323)*y(k,73)
         mat(k,607) = -rxt(k,332)*y(k,73)
         mat(k,1863) = -rxt(k,348)*y(k,73)
         mat(k,1899) = -rxt(k,363)*y(k,73)
         mat(k,805) = -rxt(k,394)*y(k,73)
         mat(k,1401) = -rxt(k,403)*y(k,73)
         mat(k,918) = -rxt(k,413)*y(k,73)
         mat(k,1782) = -rxt(k,414)*y(k,73)
         mat(k,1002) = -rxt(k,419)*y(k,73)
         mat(k,1830) = -rxt(k,432)*y(k,73)
         mat(k,1982) = -rxt(k,440)*y(k,73)
         mat(k,1498) = -rxt(k,456)*y(k,73)
         mat(k,592) = -rxt(k,488)*y(k,73)
         mat(k,277) = -rxt(k,547)*y(k,73)
         mat(k,487) = -rxt(k,553)*y(k,73)
         mat(k,1830) = mat(k,1830) + rxt(k,278)*y(k,130)
         mat(k,1401) = mat(k,1401) + rxt(k,405)*y(k,72) + rxt(k,324)*y(k,127) &
                      + rxt(k,338)*y(k,132)
         mat(k,1450) = rxt(k,405)*y(k,70)
         mat(k,384) = rxt(k,361)*y(k,135)
         mat(k,1575) = mat(k,1575) + rxt(k,324)*y(k,70)
         mat(k,607) = mat(k,607) + rxt(k,278)*y(k,63)
         mat(k,471) = rxt(k,338)*y(k,70)
         mat(k,1899) = mat(k,1899) + rxt(k,361)*y(k,125)
         mat(k,134) = -(rxt(k,554)*y(k,80))
         mat(k,473) = -rxt(k,554)*y(k,74)
         mat(k,709) = rxt(k,490)*y(k,28)
         mat(k,1085) = rxt(k,490)*y(k,5) + 2.000_r8*rxt(k,460)*y(k,28)
         mat(k,139) = -(rxt(k,543)*y(k,70) + rxt(k,544)*y(k,133))
         mat(k,1346) = -rxt(k,543)*y(k,75)
         mat(k,957) = -rxt(k,544)*y(k,75)
         mat(k,271) = -(rxt(k,545)*y(k,72) + rxt(k,547)*y(k,73) + rxt(k,550)*y(k,133))
         mat(k,1410) = -rxt(k,545)*y(k,78)
         mat(k,1904) = -rxt(k,547)*y(k,78)
         mat(k,966) = -rxt(k,550)*y(k,78)
         mat(k,476) = -(rxt(k,548)*y(k,5) + rxt(k,549)*y(k,28) + rxt(k,551)*y(k,63) &
                      + rxt(k,552)*y(k,72) + rxt(k,553)*y(k,73) + rxt(k,554)*y(k,74) &
                      + rxt(k,555)*y(k,133))
         mat(k,713) = -rxt(k,548)*y(k,80)
         mat(k,1091) = -rxt(k,549)*y(k,80)
         mat(k,1794) = -rxt(k,551)*y(k,80)
         mat(k,1418) = -rxt(k,552)*y(k,80)
         mat(k,1908) = -rxt(k,553)*y(k,80)
         mat(k,136) = -rxt(k,554)*y(k,80)
         mat(k,977) = -rxt(k,555)*y(k,80)
         mat(k,1367) = rxt(k,543)*y(k,75)
         mat(k,1418) = mat(k,1418) + rxt(k,545)*y(k,78)
         mat(k,1908) = mat(k,1908) + rxt(k,547)*y(k,78)
         mat(k,143) = rxt(k,543)*y(k,70)
         mat(k,272) = rxt(k,545)*y(k,72) + rxt(k,547)*y(k,73) + rxt(k,550)*y(k,133)
         mat(k,977) = mat(k,977) + rxt(k,550)*y(k,78)
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
         mat(k,351) = -(rxt(k,546)*y(k,133))
         mat(k,972) = -rxt(k,546)*y(k,81)
         mat(k,712) = rxt(k,548)*y(k,80)
         mat(k,1087) = rxt(k,549)*y(k,80)
         mat(k,126) = rxt(k,541)*y(k,64) + (rxt(k,542)+.500_r8*rxt(k,556))*y(k,133)
         mat(k,1790) = rxt(k,551)*y(k,80)
         mat(k,1712) = rxt(k,541)*y(k,35)
         mat(k,1413) = rxt(k,552)*y(k,80)
         mat(k,1905) = rxt(k,553)*y(k,80)
         mat(k,135) = rxt(k,554)*y(k,80)
         mat(k,142) = rxt(k,544)*y(k,133)
         mat(k,475) = rxt(k,548)*y(k,5) + rxt(k,549)*y(k,28) + rxt(k,551)*y(k,63) &
                      + rxt(k,552)*y(k,72) + rxt(k,553)*y(k,73) + rxt(k,554)*y(k,74) &
                      + rxt(k,555)*y(k,133)
         mat(k,972) = mat(k,972) + (rxt(k,542)+.500_r8*rxt(k,556))*y(k,35) &
                      + rxt(k,544)*y(k,75) + rxt(k,555)*y(k,80)
         mat(k,103) = -(rxt(k,557)*y(k,139))
         mat(k,2032) = -rxt(k,557)*y(k,82)
         mat(k,350) = rxt(k,546)*y(k,133)
         mat(k,954) = rxt(k,546)*y(k,81)
         mat(k,811) = -(rxt(k,105)*y(k,40) + rxt(k,107)*y(k,139) + rxt(k,108)*y(k,49) &
                      + rxt(k,109)*y(k,51) + rxt(k,110)*y(k,63) + rxt(k,147)*y(k,105) &
                      + rxt(k,154)*y(k,128) + rxt(k,177)*y(k,107) + rxt(k,195) &
                      *y(k,103) + rxt(k,212)*y(k,122) + rxt(k,230)*y(k,119) + rxt(k,248) &
                      *y(k,118))
         mat(k,1755) = -rxt(k,105)*y(k,90)
         mat(k,2061) = -rxt(k,107)*y(k,90)
         mat(k,2000) = -rxt(k,108)*y(k,90)
         mat(k,1513) = -rxt(k,109)*y(k,90)
         mat(k,1804) = -rxt(k,110)*y(k,90)
         mat(k,1679) = -rxt(k,147)*y(k,90)
         mat(k,1593) = -rxt(k,154)*y(k,90)
         mat(k,1636) = -rxt(k,177)*y(k,90)
         mat(k,1230) = -rxt(k,195)*y(k,90)
         mat(k,1314) = -rxt(k,212)*y(k,90)
         mat(k,1187) = -rxt(k,230)*y(k,90)
         mat(k,1144) = -rxt(k,248)*y(k,90)
         mat(k,1471) = rxt(k,114)*y(k,94) + rxt(k,127)*y(k,97) + rxt(k,285)*y(k,111) &
                      + rxt(k,314)*y(k,127) + rxt(k,341)*y(k,134) + rxt(k,350) &
                      *y(k,135)
         mat(k,1096) = rxt(k,116)*y(k,94) + rxt(k,286)*y(k,111) + rxt(k,351)*y(k,135)
         mat(k,2000) = mat(k,2000) + rxt(k,130)*y(k,97) + rxt(k,289)*y(k,111) &
                      + rxt(k,298)*y(k,113) + rxt(k,318)*y(k,127) + rxt(k,345) &
                      *y(k,134) + rxt(k,356)*y(k,135)
         mat(k,1956) = rxt(k,112)*y(k,93)
         mat(k,1374) = rxt(k,113)*y(k,93)
         mat(k,370) = rxt(k,112)*y(k,62) + rxt(k,113)*y(k,70)
         mat(k,1010) = rxt(k,114)*y(k,25) + rxt(k,116)*y(k,28)
         mat(k,832) = rxt(k,127)*y(k,25) + rxt(k,130)*y(k,49)
         mat(k,924) = rxt(k,285)*y(k,25) + rxt(k,286)*y(k,28) + rxt(k,289)*y(k,49)
         mat(k,865) = rxt(k,298)*y(k,49)
         mat(k,1548) = rxt(k,314)*y(k,25) + rxt(k,318)*y(k,49)
         mat(k,1837) = rxt(k,341)*y(k,25) + rxt(k,345)*y(k,49)
         mat(k,1873) = rxt(k,350)*y(k,25) + rxt(k,351)*y(k,28) + rxt(k,356)*y(k,49)
         mat(k,535) = -(rxt(k,106)*y(k,49) + rxt(k,159)*y(k,128) + rxt(k,164)*y(k,105) &
                      + rxt(k,182)*y(k,107) + rxt(k,199)*y(k,103) + rxt(k,217) &
                      *y(k,122) + rxt(k,234)*y(k,119) + rxt(k,252)*y(k,118))
         mat(k,1993) = -rxt(k,106)*y(k,91)
         mat(k,1585) = -rxt(k,159)*y(k,91)
         mat(k,1672) = -rxt(k,164)*y(k,91)
         mat(k,1629) = -rxt(k,182)*y(k,91)
         mat(k,1223) = -rxt(k,199)*y(k,91)
         mat(k,1307) = -rxt(k,217)*y(k,91)
         mat(k,1180) = -rxt(k,234)*y(k,91)
         mat(k,1136) = -rxt(k,252)*y(k,91)
         mat(k,810) = rxt(k,107)*y(k,139)
         mat(k,2051) = rxt(k,107)*y(k,90)
         mat(k,519) = -((rxt(k,169) + rxt(k,268)) * y(k,107) + (rxt(k,181) + rxt(k,267) &
                      ) * y(k,105) + rxt(k,187)*y(k,103) + rxt(k,205)*y(k,122) &
                      + rxt(k,222)*y(k,119) + rxt(k,240)*y(k,118) + rxt(k,257) &
                      *y(k,128))
         mat(k,1628) = -(rxt(k,169) + rxt(k,268)) * y(k,92)
         mat(k,1671) = -(rxt(k,181) + rxt(k,267)) * y(k,92)
         mat(k,1222) = -rxt(k,187)*y(k,92)
         mat(k,1306) = -rxt(k,205)*y(k,92)
         mat(k,1179) = -rxt(k,222)*y(k,92)
         mat(k,1135) = -rxt(k,240)*y(k,92)
         mat(k,1584) = -rxt(k,257)*y(k,92)
         mat(k,1992) = rxt(k,108)*y(k,90) + rxt(k,106)*y(k,91)
         mat(k,809) = rxt(k,108)*y(k,49)
         mat(k,534) = rxt(k,106)*y(k,49)
         mat(k,369) = -((rxt(k,111) + rxt(k,112)) * y(k,62) + rxt(k,113)*y(k,70))
         mat(k,1945) = -(rxt(k,111) + rxt(k,112)) * y(k,93)
         mat(k,1362) = -rxt(k,113)*y(k,93)
         mat(k,1463) = rxt(k,115)*y(k,94)
         mat(k,1089) = rxt(k,128)*y(k,97) + rxt(k,315)*y(k,127) + rxt(k,342)*y(k,134)
         mat(k,1006) = rxt(k,115)*y(k,25)
         mat(k,830) = rxt(k,128)*y(k,28)
         mat(k,1544) = rxt(k,315)*y(k,28)
         mat(k,1834) = rxt(k,342)*y(k,28)
         mat(k,1014) = -((rxt(k,114) + rxt(k,115)) * y(k,25) + rxt(k,116)*y(k,28) &
                      + rxt(k,117)*y(k,40) + rxt(k,119)*y(k,139) + rxt(k,124)*y(k,51) &
                      + rxt(k,125)*y(k,70) + rxt(k,126)*y(k,72) + (rxt(k,168) &
                      + rxt(k,265)) * y(k,107) + (rxt(k,170) + rxt(k,258)) * y(k,105) &
                      + rxt(k,186)*y(k,103) + rxt(k,204)*y(k,122) + rxt(k,221) &
                      *y(k,119) + rxt(k,239)*y(k,118) + rxt(k,256)*y(k,128) + rxt(k,280) &
                      *y(k,62) + rxt(k,281)*y(k,63))
         mat(k,1477) = -(rxt(k,114) + rxt(k,115)) * y(k,94)
         mat(k,1102) = -rxt(k,116)*y(k,94)
         mat(k,1761) = -rxt(k,117)*y(k,94)
         mat(k,2066) = -rxt(k,119)*y(k,94)
         mat(k,1518) = -rxt(k,124)*y(k,94)
         mat(k,1380) = -rxt(k,125)*y(k,94)
         mat(k,1429) = -rxt(k,126)*y(k,94)
         mat(k,1642) = -(rxt(k,168) + rxt(k,265)) * y(k,94)
         mat(k,1685) = -(rxt(k,170) + rxt(k,258)) * y(k,94)
         mat(k,1236) = -rxt(k,186)*y(k,94)
         mat(k,1320) = -rxt(k,204)*y(k,94)
         mat(k,1193) = -rxt(k,221)*y(k,94)
         mat(k,1150) = -rxt(k,239)*y(k,94)
         mat(k,1599) = -rxt(k,256)*y(k,94)
         mat(k,1961) = -rxt(k,280)*y(k,94)
         mat(k,1809) = -rxt(k,281)*y(k,94)
         mat(k,1278) = rxt(k,330)*y(k,130) + rxt(k,352)*y(k,135)
         mat(k,1761) = mat(k,1761) + rxt(k,129)*y(k,97)
         mat(k,1380) = mat(k,1380) + rxt(k,131)*y(k,97)
         mat(k,838) = rxt(k,129)*y(k,40) + rxt(k,131)*y(k,70)
         mat(k,597) = rxt(k,330)*y(k,32)
         mat(k,1878) = rxt(k,352)*y(k,32)
         mat(k,504) = -(rxt(k,153)*y(k,128) + (rxt(k,176) + rxt(k,259)) * y(k,107) &
                      + rxt(k,194)*y(k,103) + rxt(k,211)*y(k,122) + rxt(k,229) &
                      *y(k,119) + rxt(k,246)*y(k,118) + (rxt(k,247) + rxt(k,271) &
                      ) * y(k,105))
         mat(k,1583) = -rxt(k,153)*y(k,95)
         mat(k,1627) = -(rxt(k,176) + rxt(k,259)) * y(k,95)
         mat(k,1221) = -rxt(k,194)*y(k,95)
         mat(k,1305) = -rxt(k,211)*y(k,95)
         mat(k,1178) = -rxt(k,229)*y(k,95)
         mat(k,1134) = -rxt(k,246)*y(k,95)
         mat(k,1670) = -(rxt(k,247) + rxt(k,271)) * y(k,95)
         mat(k,667) = rxt(k,118)*y(k,139)
         mat(k,2050) = rxt(k,118)*y(k,96)
         mat(k,669) = -(rxt(k,118)*y(k,139) + (rxt(k,120) + rxt(k,121)) * y(k,63) &
                      + (rxt(k,122) + rxt(k,123)) * y(k,62) + rxt(k,155)*y(k,128) &
                      + (rxt(k,158) + rxt(k,260)) * y(k,105) + (rxt(k,178) + rxt(k,261) &
                      ) * y(k,107) + rxt(k,196)*y(k,103) + rxt(k,213)*y(k,122) &
                      + rxt(k,231)*y(k,119) + rxt(k,249)*y(k,118))
         mat(k,2056) = -rxt(k,118)*y(k,96)
         mat(k,1799) = -(rxt(k,120) + rxt(k,121)) * y(k,96)
         mat(k,1951) = -(rxt(k,122) + rxt(k,123)) * y(k,96)
         mat(k,1590) = -rxt(k,155)*y(k,96)
         mat(k,1677) = -(rxt(k,158) + rxt(k,260)) * y(k,96)
         mat(k,1634) = -(rxt(k,178) + rxt(k,261)) * y(k,96)
         mat(k,1228) = -rxt(k,196)*y(k,96)
         mat(k,1312) = -rxt(k,213)*y(k,96)
         mat(k,1185) = -rxt(k,231)*y(k,96)
         mat(k,1141) = -rxt(k,249)*y(k,96)
         mat(k,1008) = rxt(k,119)*y(k,139)
         mat(k,2056) = mat(k,2056) + rxt(k,119)*y(k,94)
         mat(k,833) = -(rxt(k,127)*y(k,25) + rxt(k,128)*y(k,28) + rxt(k,129)*y(k,40) &
                      + rxt(k,130)*y(k,49) + rxt(k,131)*y(k,70) + rxt(k,132)*y(k,73) &
                      + rxt(k,151)*y(k,128) + rxt(k,174)*y(k,107) + rxt(k,191) &
                      *y(k,103) + rxt(k,209)*y(k,122) + rxt(k,225)*y(k,105) + rxt(k,227) &
                      *y(k,119) + rxt(k,244)*y(k,118))
         mat(k,1472) = -rxt(k,127)*y(k,97)
         mat(k,1097) = -rxt(k,128)*y(k,97)
         mat(k,1756) = -rxt(k,129)*y(k,97)
         mat(k,2001) = -rxt(k,130)*y(k,97)
         mat(k,1375) = -rxt(k,131)*y(k,97)
         mat(k,1914) = -rxt(k,132)*y(k,97)
         mat(k,1594) = -rxt(k,151)*y(k,97)
         mat(k,1637) = -rxt(k,174)*y(k,97)
         mat(k,1231) = -rxt(k,191)*y(k,97)
         mat(k,1315) = -rxt(k,209)*y(k,97)
         mat(k,1680) = -rxt(k,225)*y(k,97)
         mat(k,1188) = -rxt(k,227)*y(k,97)
         mat(k,1145) = -rxt(k,244)*y(k,97)
         mat(k,1273) = rxt(k,316)*y(k,127) + rxt(k,335)*y(k,131)
         mat(k,1549) = rxt(k,316)*y(k,32)
         mat(k,241) = rxt(k,335)*y(k,32)
         mat(k,1058) = -((rxt(k,133) + rxt(k,135)) * y(k,72) + rxt(k,134)*y(k,73) &
                      + rxt(k,138)*y(k,99) + rxt(k,141)*y(k,105) + rxt(k,144)*y(k,107) &
                      + rxt(k,303)*y(k,119) + rxt(k,304)*y(k,120) + rxt(k,306) &
                      *y(k,121) + rxt(k,308)*y(k,122) + rxt(k,326)*y(k,129) + rxt(k,380) &
                      *y(k,132) + rxt(k,381)*y(k,108) + rxt(k,382)*y(k,102) + rxt(k,383) &
                      *y(k,103) + rxt(k,384)*y(k,123) + rxt(k,582)*y(k,118) + rxt(k,583) &
                      *y(k,128) + rxt(k,584)*y(k,110))
         mat(k,1430) = -(rxt(k,133) + rxt(k,135)) * y(k,98)
         mat(k,1920) = -rxt(k,134)*y(k,98)
         mat(k,170) = -rxt(k,138)*y(k,98)
         mat(k,1686) = -rxt(k,141)*y(k,98)
         mat(k,1643) = -rxt(k,144)*y(k,98)
         mat(k,1194) = -rxt(k,303)*y(k,98)
         mat(k,175) = -rxt(k,304)*y(k,98)
         mat(k,265) = -rxt(k,306)*y(k,98)
         mat(k,1321) = -rxt(k,308)*y(k,98)
         mat(k,294) = -rxt(k,326)*y(k,98)
         mat(k,463) = -rxt(k,380)*y(k,98)
         mat(k,285) = -rxt(k,381)*y(k,98)
         mat(k,248) = -rxt(k,382)*y(k,98)
         mat(k,1237) = -rxt(k,383)*y(k,98)
         mat(k,258) = -rxt(k,384)*y(k,98)
         mat(k,1151) = -rxt(k,582)*y(k,98)
         mat(k,1600) = -rxt(k,583)*y(k,98)
         mat(k,234) = -rxt(k,584)*y(k,98)
         mat(k,1762) = rxt(k,105)*y(k,90) + rxt(k,317)*y(k,127) + rxt(k,344)*y(k,134)
         mat(k,652) = rxt(k,353)*y(k,135)
         mat(k,1962) = rxt(k,136)*y(k,135)
         mat(k,1381) = rxt(k,324)*y(k,127) + rxt(k,333)*y(k,130) + rxt(k,347)*y(k,134) &
                      + rxt(k,360)*y(k,135)
         mat(k,1920) = mat(k,1920) + rxt(k,332)*y(k,130)
         mat(k,814) = rxt(k,105)*y(k,40)
         mat(k,379) = rxt(k,321)*y(k,127) + rxt(k,361)*y(k,135)
         mat(k,1555) = rxt(k,317)*y(k,40) + rxt(k,324)*y(k,70) + rxt(k,321)*y(k,125)
         mat(k,598) = rxt(k,333)*y(k,70) + rxt(k,332)*y(k,73)
         mat(k,1843) = rxt(k,344)*y(k,40) + rxt(k,347)*y(k,70)
         mat(k,1879) = rxt(k,353)*y(k,41) + rxt(k,136)*y(k,62) + rxt(k,360)*y(k,70) &
                      + rxt(k,361)*y(k,125)
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
         mat(k,167) = -(rxt(k,138)*y(k,98) + rxt(k,139)*y(k,139))
         mat(k,1039) = -rxt(k,138)*y(k,99)
         mat(k,2036) = -rxt(k,139)*y(k,99)
         mat(k,289) = rxt(k,327)*y(k,139)
         mat(k,2036) = mat(k,2036) + rxt(k,327)*y(k,129)
         mat(k,490) = -(rxt(k,149)*y(k,128) + rxt(k,172)*y(k,107) + rxt(k,189) &
                      *y(k,103) + rxt(k,203)*y(k,105) + rxt(k,207)*y(k,122) + rxt(k,224) &
                      *y(k,119) + rxt(k,242)*y(k,118))
         mat(k,1582) = -rxt(k,149)*y(k,100)
         mat(k,1626) = -rxt(k,172)*y(k,100)
         mat(k,1220) = -rxt(k,189)*y(k,100)
         mat(k,1669) = -rxt(k,203)*y(k,100)
         mat(k,1304) = -rxt(k,207)*y(k,100)
         mat(k,1177) = -rxt(k,224)*y(k,100)
         mat(k,1133) = -rxt(k,242)*y(k,100)
         mat(k,1270) = rxt(k,343)*y(k,134)
         mat(k,1835) = rxt(k,343)*y(k,32)
         mat(k,900) = -(rxt(k,310)*y(k,122) + (rxt(k,408) + rxt(k,409) + rxt(k,410) &
                      ) * y(k,40) + rxt(k,412)*y(k,70) + rxt(k,413)*y(k,73) + rxt(k,417) &
                      *y(k,133) + 4._r8*rxt(k,422)*y(k,101) + rxt(k,434)*y(k,64) &
                      + rxt(k,439)*y(k,62) + rxt(k,444)*y(k,63) + (rxt(k,454) &
                      + rxt(k,455)) * y(k,25) + rxt(k,461)*y(k,28) + rxt(k,487)*y(k,3) &
                      + rxt(k,493)*y(k,5) + rxt(k,533)*y(k,21))
         mat(k,1317) = -rxt(k,310)*y(k,101)
         mat(k,1758) = -(rxt(k,408) + rxt(k,409) + rxt(k,410)) * y(k,101)
         mat(k,1377) = -rxt(k,412)*y(k,101)
         mat(k,1916) = -rxt(k,413)*y(k,101)
         mat(k,983) = -rxt(k,417)*y(k,101)
         mat(k,1723) = -rxt(k,434)*y(k,101)
         mat(k,1958) = -rxt(k,439)*y(k,101)
         mat(k,1806) = -rxt(k,444)*y(k,101)
         mat(k,1474) = -(rxt(k,454) + rxt(k,455)) * y(k,101)
         mat(k,1099) = -rxt(k,461)*y(k,101)
         mat(k,586) = -rxt(k,487)*y(k,101)
         mat(k,718) = -rxt(k,493)*y(k,101)
         mat(k,434) = -rxt(k,533)*y(k,101)
         mat(k,586) = mat(k,586) + rxt(k,486)*y(k,17)
         mat(k,718) = mat(k,718) + rxt(k,498)*y(k,133)
         mat(k,740) = rxt(k,486)*y(k,3) + rxt(k,450)*y(k,25) + rxt(k,530)*y(k,64) &
                      + rxt(k,531)*y(k,70)
         mat(k,225) = rxt(k,515)*y(k,25) + rxt(k,516)*y(k,133)
         mat(k,150) = rxt(k,518)*y(k,25) + rxt(k,519)*y(k,133)
         mat(k,434) = mat(k,434) + rxt(k,457)*y(k,28) + rxt(k,534)*y(k,62)
         mat(k,421) = rxt(k,538)*y(k,71)
         mat(k,1474) = mat(k,1474) + rxt(k,450)*y(k,17) + rxt(k,515)*y(k,18) &
                      + rxt(k,518)*y(k,20) + rxt(k,453)*y(k,43)
         mat(k,1099) = mat(k,1099) + rxt(k,457)*y(k,21) + rxt(k,468)*y(k,133)
         mat(k,253) = rxt(k,540)*y(k,133)
         mat(k,127) = .500_r8*rxt(k,556)*y(k,133)
         mat(k,1758) = mat(k,1758) + rxt(k,411)*y(k,72) + rxt(k,317)*y(k,127)
         mat(k,215) = rxt(k,453)*y(k,25) + rxt(k,407)*y(k,70) + rxt(k,416)*y(k,133)
         mat(k,2003) = rxt(k,130)*y(k,97) + rxt(k,318)*y(k,127)
         mat(k,1515) = rxt(k,319)*y(k,127)
         mat(k,1958) = mat(k,1958) + rxt(k,534)*y(k,21)
         mat(k,1723) = mat(k,1723) + rxt(k,530)*y(k,17) + rxt(k,437)*y(k,133)
         mat(k,1377) = mat(k,1377) + rxt(k,531)*y(k,17) + rxt(k,407)*y(k,43) &
                      + rxt(k,347)*y(k,134)
         mat(k,788) = rxt(k,538)*y(k,23)
         mat(k,1426) = rxt(k,411)*y(k,40)
         mat(k,1916) = mat(k,1916) + rxt(k,419)*y(k,133)
         mat(k,353) = rxt(k,546)*y(k,133)
         mat(k,835) = rxt(k,130)*y(k,49)
         mat(k,1551) = rxt(k,317)*y(k,40) + rxt(k,318)*y(k,49) + rxt(k,319)*y(k,51)
         mat(k,983) = mat(k,983) + rxt(k,498)*y(k,5) + rxt(k,516)*y(k,18) + rxt(k,519) &
                      *y(k,20) + rxt(k,468)*y(k,28) + rxt(k,540)*y(k,31) &
                      + .500_r8*rxt(k,556)*y(k,35) + rxt(k,416)*y(k,43) + rxt(k,437) &
                      *y(k,64) + rxt(k,419)*y(k,73) + rxt(k,546)*y(k,81)
         mat(k,1839) = rxt(k,347)*y(k,70)
         mat(k,246) = -(rxt(k,374)*y(k,139) + rxt(k,382)*y(k,98))
         mat(k,2038) = -rxt(k,374)*y(k,102)
         mat(k,1044) = -rxt(k,382)*y(k,102)
         mat(k,168) = rxt(k,139)*y(k,139)
         mat(k,283) = rxt(k,372)*y(k,139)
         mat(k,2038) = mat(k,2038) + rxt(k,139)*y(k,99) + rxt(k,372)*y(k,108)
         mat(k,1241) = -(rxt(k,185)*y(k,117) + rxt(k,186)*y(k,94) + rxt(k,187)*y(k,92) &
                      + rxt(k,188)*y(k,113) + rxt(k,189)*y(k,100) + rxt(k,190) &
                      *y(k,127) + rxt(k,191)*y(k,97) + rxt(k,193)*y(k,115) + rxt(k,194) &
                      *y(k,95) + rxt(k,195)*y(k,90) + rxt(k,196)*y(k,96) + rxt(k,197) &
                      *y(k,112) + rxt(k,198)*y(k,116) + rxt(k,199)*y(k,91) + rxt(k,200) &
                      *y(k,114) + rxt(k,201)*y(k,111) + rxt(k,376)*y(k,139) + rxt(k,383) &
                      *y(k,98))
         mat(k,613) = -rxt(k,185)*y(k,103)
         mat(k,1019) = -rxt(k,186)*y(k,103)
         mat(k,523) = -rxt(k,187)*y(k,103)
         mat(k,871) = -rxt(k,188)*y(k,103)
         mat(k,494) = -rxt(k,189)*y(k,103)
         mat(k,1559) = -rxt(k,190)*y(k,103)
         mat(k,843) = -rxt(k,191)*y(k,103)
         mat(k,695) = -rxt(k,193)*y(k,103)
         mat(k,508) = -rxt(k,194)*y(k,103)
         mat(k,817) = -rxt(k,195)*y(k,103)
         mat(k,676) = -rxt(k,196)*y(k,103)
         mat(k,554) = -rxt(k,197)*y(k,103)
         mat(k,568) = -rxt(k,198)*y(k,103)
         mat(k,539) = -rxt(k,199)*y(k,103)
         mat(k,630) = -rxt(k,200)*y(k,103)
         mat(k,932) = -rxt(k,201)*y(k,103)
         mat(k,2071) = -rxt(k,376)*y(k,103)
         mat(k,1062) = -rxt(k,383)*y(k,103)
         mat(k,249) = rxt(k,374)*y(k,139)
         mat(k,176) = rxt(k,305)*y(k,139)
         mat(k,2071) = mat(k,2071) + rxt(k,374)*y(k,102) + rxt(k,305)*y(k,120)
         mat(k,107) = -(rxt(k,140)*y(k,139))
         mat(k,2033) = -rxt(k,140)*y(k,104)
         mat(k,399) = rxt(k,142)*y(k,105)
         mat(k,1667) = rxt(k,142)*y(k,58)
         mat(k,1700) = -(rxt(k,141)*y(k,98) + rxt(k,142)*y(k,58) + (rxt(k,146) &
                      + rxt(k,269)) * y(k,117) + rxt(k,147)*y(k,90) + (rxt(k,158) &
                      + rxt(k,260)) * y(k,96) + rxt(k,162)*y(k,112) + rxt(k,163) &
                      *y(k,116) + rxt(k,164)*y(k,91) + rxt(k,165)*y(k,114) + rxt(k,166) &
                      *y(k,111) + (rxt(k,170) + rxt(k,258)) * y(k,94) + (rxt(k,181) &
                      + rxt(k,267)) * y(k,92) + (rxt(k,192) + rxt(k,264)) * y(k,113) &
                      + rxt(k,203)*y(k,100) + rxt(k,214)*y(k,127) + rxt(k,225)*y(k,97) &
                      + (rxt(k,236) + rxt(k,262)) * y(k,115) + (rxt(k,247) + rxt(k,271) &
                      ) * y(k,95) + rxt(k,378)*y(k,139))
         mat(k,1072) = -rxt(k,141)*y(k,105)
         mat(k,408) = -rxt(k,142)*y(k,105)
         mat(k,619) = -(rxt(k,146) + rxt(k,269)) * y(k,105)
         mat(k,824) = -rxt(k,147)*y(k,105)
         mat(k,683) = -(rxt(k,158) + rxt(k,260)) * y(k,105)
         mat(k,559) = -rxt(k,162)*y(k,105)
         mat(k,574) = -rxt(k,163)*y(k,105)
         mat(k,545) = -rxt(k,164)*y(k,105)
         mat(k,637) = -rxt(k,165)*y(k,105)
         mat(k,940) = -rxt(k,166)*y(k,105)
         mat(k,1029) = -(rxt(k,170) + rxt(k,258)) * y(k,105)
         mat(k,529) = -(rxt(k,181) + rxt(k,267)) * y(k,105)
         mat(k,879) = -(rxt(k,192) + rxt(k,264)) * y(k,105)
         mat(k,500) = -rxt(k,203)*y(k,105)
         mat(k,1569) = -rxt(k,214)*y(k,105)
         mat(k,853) = -rxt(k,225)*y(k,105)
         mat(k,702) = -(rxt(k,236) + rxt(k,262)) * y(k,105)
         mat(k,515) = -(rxt(k,247) + rxt(k,271)) * y(k,105)
         mat(k,2081) = -rxt(k,378)*y(k,105)
         mat(k,1251) = rxt(k,376)*y(k,139)
         mat(k,109) = rxt(k,140)*y(k,139)
         mat(k,2081) = mat(k,2081) + rxt(k,376)*y(k,103) + rxt(k,140)*y(k,104)
         mat(k,111) = -(rxt(k,143)*y(k,139))
         mat(k,2034) = -rxt(k,143)*y(k,106)
         mat(k,400) = rxt(k,145)*y(k,107)
         mat(k,1624) = rxt(k,145)*y(k,58)
         mat(k,1656) = -(rxt(k,144)*y(k,98) + rxt(k,145)*y(k,58) + (rxt(k,167) &
                      + rxt(k,270)) * y(k,117) + (rxt(k,168) + rxt(k,265)) * y(k,94) &
                      + (rxt(k,169) + rxt(k,268)) * y(k,92) + (rxt(k,171) + rxt(k,266) &
                      ) * y(k,113) + rxt(k,172)*y(k,100) + rxt(k,173)*y(k,127) &
                      + rxt(k,174)*y(k,97) + (rxt(k,175) + rxt(k,263)) * y(k,115) &
                      + (rxt(k,176) + rxt(k,259)) * y(k,95) + rxt(k,177)*y(k,90) &
                      + (rxt(k,178) + rxt(k,261)) * y(k,96) + rxt(k,179)*y(k,112) &
                      + rxt(k,180)*y(k,116) + rxt(k,182)*y(k,91) + rxt(k,183)*y(k,114) &
                      + rxt(k,184)*y(k,111))
         mat(k,1071) = -rxt(k,144)*y(k,107)
         mat(k,407) = -rxt(k,145)*y(k,107)
         mat(k,618) = -(rxt(k,167) + rxt(k,270)) * y(k,107)
         mat(k,1028) = -(rxt(k,168) + rxt(k,265)) * y(k,107)
         mat(k,528) = -(rxt(k,169) + rxt(k,268)) * y(k,107)
         mat(k,878) = -(rxt(k,171) + rxt(k,266)) * y(k,107)
         mat(k,499) = -rxt(k,172)*y(k,107)
         mat(k,1568) = -rxt(k,173)*y(k,107)
         mat(k,852) = -rxt(k,174)*y(k,107)
         mat(k,701) = -(rxt(k,175) + rxt(k,263)) * y(k,107)
         mat(k,514) = -(rxt(k,176) + rxt(k,259)) * y(k,107)
         mat(k,823) = -rxt(k,177)*y(k,107)
         mat(k,682) = -(rxt(k,178) + rxt(k,261)) * y(k,107)
         mat(k,558) = -rxt(k,179)*y(k,107)
         mat(k,573) = -rxt(k,180)*y(k,107)
         mat(k,544) = -rxt(k,182)*y(k,107)
         mat(k,636) = -rxt(k,183)*y(k,107)
         mat(k,939) = -rxt(k,184)*y(k,107)
         mat(k,1699) = rxt(k,378)*y(k,139)
         mat(k,113) = rxt(k,143)*y(k,139)
         mat(k,2080) = rxt(k,378)*y(k,105) + rxt(k,143)*y(k,106)
         mat(k,284) = -(rxt(k,372)*y(k,139) + rxt(k,381)*y(k,98))
         mat(k,2042) = -rxt(k,372)*y(k,108)
         mat(k,1048) = -rxt(k,381)*y(k,108)
         mat(k,1750) = rxt(k,309)*y(k,122)
         mat(k,890) = rxt(k,310)*y(k,122)
         mat(k,1303) = rxt(k,309)*y(k,40) + rxt(k,310)*y(k,101) + rxt(k,311)*y(k,133)
         mat(k,291) = rxt(k,328)*y(k,139)
         mat(k,968) = rxt(k,311)*y(k,122)
         mat(k,2042) = mat(k,2042) + rxt(k,328)*y(k,129)
         mat(k,300) = -(rxt(k,424)*y(k,70) + rxt(k,425)*y(k,72) + rxt(k,596)*y(k,136))
         mat(k,1355) = -rxt(k,424)*y(k,109)
         mat(k,1411) = -rxt(k,425)*y(k,109)
         mat(k,386) = -rxt(k,596)*y(k,109)
         mat(k,1355) = mat(k,1355) + rxt(k,586)*y(k,110)
         mat(k,1050) = .900_r8*rxt(k,584)*y(k,110) + .800_r8*rxt(k,582)*y(k,118)
         mat(k,231) = rxt(k,586)*y(k,70) + .900_r8*rxt(k,584)*y(k,98)
         mat(k,1129) = .800_r8*rxt(k,582)*y(k,98)
         mat(k,230) = -(rxt(k,584)*y(k,98) + rxt(k,585)*y(k,72) + (rxt(k,586) &
                      + rxt(k,587)) * y(k,70))
         mat(k,1043) = -rxt(k,584)*y(k,110)
         mat(k,1408) = -rxt(k,585)*y(k,110)
         mat(k,1351) = -(rxt(k,586) + rxt(k,587)) * y(k,110)
         mat(k,926) = -(rxt(k,161)*y(k,128) + rxt(k,166)*y(k,105) + rxt(k,184) &
                      *y(k,107) + rxt(k,201)*y(k,103) + rxt(k,219)*y(k,122) + rxt(k,237) &
                      *y(k,119) + rxt(k,254)*y(k,118) + rxt(k,285)*y(k,25) + rxt(k,286) &
                      *y(k,28) + rxt(k,287)*y(k,40) + rxt(k,288)*y(k,139) + rxt(k,289) &
                      *y(k,49) + rxt(k,290)*y(k,51) + rxt(k,291)*y(k,63) + rxt(k,292) &
                      *y(k,73))
         mat(k,1597) = -rxt(k,161)*y(k,111)
         mat(k,1683) = -rxt(k,166)*y(k,111)
         mat(k,1640) = -rxt(k,184)*y(k,111)
         mat(k,1234) = -rxt(k,201)*y(k,111)
         mat(k,1318) = -rxt(k,219)*y(k,111)
         mat(k,1191) = -rxt(k,237)*y(k,111)
         mat(k,1148) = -rxt(k,254)*y(k,111)
         mat(k,1475) = -rxt(k,285)*y(k,111)
         mat(k,1100) = -rxt(k,286)*y(k,111)
         mat(k,1759) = -rxt(k,287)*y(k,111)
         mat(k,2064) = -rxt(k,288)*y(k,111)
         mat(k,2004) = -rxt(k,289)*y(k,111)
         mat(k,1516) = -rxt(k,290)*y(k,111)
         mat(k,1807) = -rxt(k,291)*y(k,111)
         mat(k,1917) = -rxt(k,292)*y(k,111)
         mat(k,1959) = rxt(k,111)*y(k,93) + rxt(k,280)*y(k,94) + rxt(k,123)*y(k,96) &
                      + rxt(k,279)*y(k,130)
         mat(k,1807) = mat(k,1807) + rxt(k,110)*y(k,90) + rxt(k,320)*y(k,127) &
                      + rxt(k,278)*y(k,130) + rxt(k,346)*y(k,134) + rxt(k,359) &
                      *y(k,135)
         mat(k,1378) = rxt(k,301)*y(k,113)
         mat(k,1917) = mat(k,1917) + rxt(k,302)*y(k,113)
         mat(k,813) = rxt(k,110)*y(k,63)
         mat(k,371) = rxt(k,111)*y(k,62)
         mat(k,1012) = rxt(k,280)*y(k,62)
         mat(k,672) = rxt(k,123)*y(k,62)
         mat(k,867) = rxt(k,301)*y(k,70) + rxt(k,302)*y(k,73)
         mat(k,1552) = rxt(k,320)*y(k,63)
         mat(k,596) = rxt(k,279)*y(k,62) + rxt(k,278)*y(k,63)
         mat(k,1840) = rxt(k,346)*y(k,63)
         mat(k,1876) = rxt(k,359)*y(k,63)
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
         mat(k,550) = -(rxt(k,156)*y(k,128) + rxt(k,162)*y(k,105) + rxt(k,179) &
                      *y(k,107) + rxt(k,197)*y(k,103) + rxt(k,215)*y(k,122) + rxt(k,232) &
                      *y(k,119) + rxt(k,250)*y(k,118))
         mat(k,1586) = -rxt(k,156)*y(k,112)
         mat(k,1673) = -rxt(k,162)*y(k,112)
         mat(k,1630) = -rxt(k,179)*y(k,112)
         mat(k,1224) = -rxt(k,197)*y(k,112)
         mat(k,1308) = -rxt(k,215)*y(k,112)
         mat(k,1181) = -rxt(k,232)*y(k,112)
         mat(k,1137) = -rxt(k,250)*y(k,112)
         mat(k,1948) = rxt(k,122)*y(k,96)
         mat(k,668) = rxt(k,122)*y(k,62)
         mat(k,923) = rxt(k,288)*y(k,139)
         mat(k,2052) = rxt(k,288)*y(k,111)
         mat(k,866) = -(rxt(k,148)*y(k,128) + (rxt(k,171) + rxt(k,266)) * y(k,107) &
                      + rxt(k,188)*y(k,103) + (rxt(k,192) + rxt(k,264)) * y(k,105) &
                      + rxt(k,206)*y(k,122) + rxt(k,223)*y(k,119) + rxt(k,241) &
                      *y(k,118) + (rxt(k,276) + rxt(k,298)) * y(k,49) + rxt(k,296) &
                      *y(k,139) + rxt(k,300)*y(k,51) + rxt(k,301)*y(k,70) + rxt(k,302) &
                      *y(k,73))
         mat(k,1595) = -rxt(k,148)*y(k,113)
         mat(k,1638) = -(rxt(k,171) + rxt(k,266)) * y(k,113)
         mat(k,1232) = -rxt(k,188)*y(k,113)
         mat(k,1681) = -(rxt(k,192) + rxt(k,264)) * y(k,113)
         mat(k,1316) = -rxt(k,206)*y(k,113)
         mat(k,1189) = -rxt(k,223)*y(k,113)
         mat(k,1146) = -rxt(k,241)*y(k,113)
         mat(k,2002) = -(rxt(k,276) + rxt(k,298)) * y(k,113)
         mat(k,2062) = -rxt(k,296)*y(k,113)
         mat(k,1514) = -rxt(k,300)*y(k,113)
         mat(k,1376) = -rxt(k,301)*y(k,113)
         mat(k,1915) = -rxt(k,302)*y(k,113)
         mat(k,1514) = mat(k,1514) + rxt(k,109)*y(k,90) + rxt(k,124)*y(k,94) &
                      + rxt(k,290)*y(k,111) + rxt(k,319)*y(k,127) + rxt(k,357) &
                      *y(k,135)
         mat(k,1957) = rxt(k,272)*y(k,130)
         mat(k,1805) = rxt(k,281)*y(k,94) + rxt(k,120)*y(k,96) + rxt(k,291)*y(k,111) &
                      + rxt(k,277)*y(k,130)
         mat(k,1915) = mat(k,1915) + rxt(k,292)*y(k,111)
         mat(k,812) = rxt(k,109)*y(k,51)
         mat(k,1011) = rxt(k,124)*y(k,51) + rxt(k,281)*y(k,63)
         mat(k,671) = rxt(k,120)*y(k,63)
         mat(k,925) = rxt(k,290)*y(k,51) + rxt(k,291)*y(k,63) + rxt(k,292)*y(k,73)
         mat(k,1550) = rxt(k,319)*y(k,51)
         mat(k,595) = rxt(k,272)*y(k,62) + rxt(k,277)*y(k,63)
         mat(k,1874) = rxt(k,357)*y(k,51)
         mat(k,625) = -(rxt(k,160)*y(k,128) + rxt(k,165)*y(k,105) + rxt(k,183) &
                      *y(k,107) + rxt(k,200)*y(k,103) + rxt(k,218)*y(k,122) + rxt(k,235) &
                      *y(k,119) + rxt(k,253)*y(k,118) + rxt(k,293)*y(k,58))
         mat(k,1589) = -rxt(k,160)*y(k,114)
         mat(k,1676) = -rxt(k,165)*y(k,114)
         mat(k,1633) = -rxt(k,183)*y(k,114)
         mat(k,1227) = -rxt(k,200)*y(k,114)
         mat(k,1311) = -rxt(k,218)*y(k,114)
         mat(k,1184) = -rxt(k,235)*y(k,114)
         mat(k,1140) = -rxt(k,253)*y(k,114)
         mat(k,403) = -rxt(k,293)*y(k,114)
         mat(k,690) = rxt(k,294)*y(k,139)
         mat(k,2054) = rxt(k,294)*y(k,115)
         mat(k,691) = -(rxt(k,152)*y(k,128) + (rxt(k,175) + rxt(k,263)) * y(k,107) &
                      + rxt(k,193)*y(k,103) + rxt(k,210)*y(k,122) + rxt(k,228) &
                      *y(k,119) + (rxt(k,236) + rxt(k,262)) * y(k,105) + rxt(k,245) &
                      *y(k,118) + rxt(k,294)*y(k,139) + rxt(k,295)*y(k,51) + rxt(k,297) &
                      *y(k,58))
         mat(k,1591) = -rxt(k,152)*y(k,115)
         mat(k,1635) = -(rxt(k,175) + rxt(k,263)) * y(k,115)
         mat(k,1229) = -rxt(k,193)*y(k,115)
         mat(k,1313) = -rxt(k,210)*y(k,115)
         mat(k,1186) = -rxt(k,228)*y(k,115)
         mat(k,1678) = -(rxt(k,236) + rxt(k,262)) * y(k,115)
         mat(k,1142) = -rxt(k,245)*y(k,115)
         mat(k,2057) = -rxt(k,294)*y(k,115)
         mat(k,1510) = -rxt(k,295)*y(k,115)
         mat(k,404) = -rxt(k,297)*y(k,115)
         mat(k,1800) = rxt(k,121)*y(k,96)
         mat(k,670) = rxt(k,121)*y(k,63)
         mat(k,864) = rxt(k,296)*y(k,139)
         mat(k,2057) = mat(k,2057) + rxt(k,296)*y(k,113)
         mat(k,564) = -(rxt(k,157)*y(k,128) + rxt(k,163)*y(k,105) + rxt(k,180) &
                      *y(k,107) + rxt(k,198)*y(k,103) + rxt(k,216)*y(k,122) + rxt(k,233) &
                      *y(k,119) + rxt(k,251)*y(k,118) + rxt(k,299)*y(k,51))
         mat(k,1587) = -rxt(k,157)*y(k,116)
         mat(k,1674) = -rxt(k,163)*y(k,116)
         mat(k,1631) = -rxt(k,180)*y(k,116)
         mat(k,1225) = -rxt(k,198)*y(k,116)
         mat(k,1309) = -rxt(k,216)*y(k,116)
         mat(k,1182) = -rxt(k,233)*y(k,116)
         mat(k,1138) = -rxt(k,251)*y(k,116)
         mat(k,1507) = -rxt(k,299)*y(k,116)
         mat(k,1994) = rxt(k,276)*y(k,113)
         mat(k,862) = rxt(k,276)*y(k,49)
         mat(k,609) = -((rxt(k,146) + rxt(k,269)) * y(k,105) + (rxt(k,167) + rxt(k,270) &
                      ) * y(k,107) + rxt(k,185)*y(k,103) + rxt(k,202)*y(k,122) &
                      + rxt(k,220)*y(k,119) + rxt(k,238)*y(k,118) + rxt(k,255) &
                      *y(k,128))
         mat(k,1675) = -(rxt(k,146) + rxt(k,269)) * y(k,117)
         mat(k,1632) = -(rxt(k,167) + rxt(k,270)) * y(k,117)
         mat(k,1226) = -rxt(k,185)*y(k,117)
         mat(k,1310) = -rxt(k,202)*y(k,117)
         mat(k,1183) = -rxt(k,220)*y(k,117)
         mat(k,1139) = -rxt(k,238)*y(k,117)
         mat(k,1588) = -rxt(k,255)*y(k,117)
         mat(k,1508) = rxt(k,300)*y(k,113) + rxt(k,295)*y(k,115) + rxt(k,299)*y(k,116)
         mat(k,402) = rxt(k,293)*y(k,114) + rxt(k,297)*y(k,115)
         mat(k,863) = rxt(k,300)*y(k,51)
         mat(k,624) = rxt(k,293)*y(k,58)
         mat(k,689) = rxt(k,295)*y(k,51) + rxt(k,297)*y(k,58)
         mat(k,565) = rxt(k,299)*y(k,51)
         mat(k,1153) = -(rxt(k,238)*y(k,117) + rxt(k,239)*y(k,94) + rxt(k,240)*y(k,92) &
                      + rxt(k,241)*y(k,113) + rxt(k,242)*y(k,100) + rxt(k,243) &
                      *y(k,127) + rxt(k,244)*y(k,97) + rxt(k,245)*y(k,115) + rxt(k,246) &
                      *y(k,95) + rxt(k,248)*y(k,90) + rxt(k,249)*y(k,96) + rxt(k,250) &
                      *y(k,112) + rxt(k,251)*y(k,116) + rxt(k,252)*y(k,91) + rxt(k,253) &
                      *y(k,114) + rxt(k,254)*y(k,111) + rxt(k,365)*y(k,139) + rxt(k,368) &
                      *y(k,32) + rxt(k,582)*y(k,98))
         mat(k,611) = -rxt(k,238)*y(k,118)
         mat(k,1017) = -rxt(k,239)*y(k,118)
         mat(k,521) = -rxt(k,240)*y(k,118)
         mat(k,869) = -rxt(k,241)*y(k,118)
         mat(k,492) = -rxt(k,242)*y(k,118)
         mat(k,1557) = -rxt(k,243)*y(k,118)
         mat(k,841) = -rxt(k,244)*y(k,118)
         mat(k,693) = -rxt(k,245)*y(k,118)
         mat(k,506) = -rxt(k,246)*y(k,118)
         mat(k,815) = -rxt(k,248)*y(k,118)
         mat(k,674) = -rxt(k,249)*y(k,118)
         mat(k,552) = -rxt(k,250)*y(k,118)
         mat(k,566) = -rxt(k,251)*y(k,118)
         mat(k,537) = -rxt(k,252)*y(k,118)
         mat(k,628) = -rxt(k,253)*y(k,118)
         mat(k,930) = -rxt(k,254)*y(k,118)
         mat(k,2069) = -rxt(k,365)*y(k,118)
         mat(k,1281) = -rxt(k,368)*y(k,118)
         mat(k,1060) = -rxt(k,582)*y(k,118)
         mat(k,451) = rxt(k,591)*y(k,128)
         mat(k,1964) = rxt(k,593)*y(k,128)
         mat(k,1383) = rxt(k,586)*y(k,110)
         mat(k,1432) = rxt(k,590)*y(k,124)
         mat(k,235) = rxt(k,586)*y(k,70)
         mat(k,320) = rxt(k,590)*y(k,72)
         mat(k,1602) = rxt(k,591)*y(k,56) + rxt(k,593)*y(k,62)
         mat(k,1197) = -(rxt(k,220)*y(k,117) + rxt(k,221)*y(k,94) + rxt(k,222)*y(k,92) &
                      + rxt(k,223)*y(k,113) + rxt(k,224)*y(k,100) + rxt(k,226) &
                      *y(k,127) + rxt(k,227)*y(k,97) + rxt(k,228)*y(k,115) + rxt(k,229) &
                      *y(k,95) + rxt(k,230)*y(k,90) + rxt(k,231)*y(k,96) + rxt(k,232) &
                      *y(k,112) + rxt(k,233)*y(k,116) + rxt(k,234)*y(k,91) + rxt(k,235) &
                      *y(k,114) + rxt(k,237)*y(k,111) + rxt(k,303)*y(k,98) + rxt(k,367) &
                      *y(k,139))
         mat(k,612) = -rxt(k,220)*y(k,119)
         mat(k,1018) = -rxt(k,221)*y(k,119)
         mat(k,522) = -rxt(k,222)*y(k,119)
         mat(k,870) = -rxt(k,223)*y(k,119)
         mat(k,493) = -rxt(k,224)*y(k,119)
         mat(k,1558) = -rxt(k,226)*y(k,119)
         mat(k,842) = -rxt(k,227)*y(k,119)
         mat(k,694) = -rxt(k,228)*y(k,119)
         mat(k,507) = -rxt(k,229)*y(k,119)
         mat(k,816) = -rxt(k,230)*y(k,119)
         mat(k,675) = -rxt(k,231)*y(k,119)
         mat(k,553) = -rxt(k,232)*y(k,119)
         mat(k,567) = -rxt(k,233)*y(k,119)
         mat(k,538) = -rxt(k,234)*y(k,119)
         mat(k,629) = -rxt(k,235)*y(k,119)
         mat(k,931) = -rxt(k,237)*y(k,119)
         mat(k,1061) = -rxt(k,303)*y(k,119)
         mat(k,2070) = -rxt(k,367)*y(k,119)
         mat(k,1324) = rxt(k,366)*y(k,139)
         mat(k,2070) = mat(k,2070) + rxt(k,366)*y(k,122)
         mat(k,173) = -(rxt(k,304)*y(k,98) + rxt(k,305)*y(k,139))
         mat(k,1040) = -rxt(k,304)*y(k,120)
         mat(k,2037) = -rxt(k,305)*y(k,120)
         mat(k,1175) = rxt(k,367)*y(k,139)
         mat(k,2037) = mat(k,2037) + rxt(k,367)*y(k,119)
         mat(k,264) = -(rxt(k,306)*y(k,98) + rxt(k,307)*y(k,139))
         mat(k,1046) = -rxt(k,306)*y(k,121)
         mat(k,2040) = -rxt(k,307)*y(k,121)
         mat(k,1265) = rxt(k,368)*y(k,118) + rxt(k,312)*y(k,123)
         mat(k,1128) = rxt(k,368)*y(k,32)
         mat(k,257) = rxt(k,312)*y(k,32)
         mat(k,1327) = -(rxt(k,202)*y(k,117) + rxt(k,204)*y(k,94) + rxt(k,205)*y(k,92) &
                      + rxt(k,206)*y(k,113) + rxt(k,207)*y(k,100) + rxt(k,208) &
                      *y(k,127) + rxt(k,209)*y(k,97) + rxt(k,210)*y(k,115) + rxt(k,211) &
                      *y(k,95) + rxt(k,212)*y(k,90) + rxt(k,213)*y(k,96) + rxt(k,215) &
                      *y(k,112) + rxt(k,216)*y(k,116) + rxt(k,217)*y(k,91) + rxt(k,218) &
                      *y(k,114) + rxt(k,219)*y(k,111) + rxt(k,308)*y(k,98) + rxt(k,309) &
                      *y(k,40) + rxt(k,310)*y(k,101) + rxt(k,311)*y(k,133) + rxt(k,366) &
                      *y(k,139))
         mat(k,614) = -rxt(k,202)*y(k,122)
         mat(k,1021) = -rxt(k,204)*y(k,122)
         mat(k,524) = -rxt(k,205)*y(k,122)
         mat(k,872) = -rxt(k,206)*y(k,122)
         mat(k,496) = -rxt(k,207)*y(k,122)
         mat(k,1561) = -rxt(k,208)*y(k,122)
         mat(k,845) = -rxt(k,209)*y(k,122)
         mat(k,696) = -rxt(k,210)*y(k,122)
         mat(k,510) = -rxt(k,211)*y(k,122)
         mat(k,818) = -rxt(k,212)*y(k,122)
         mat(k,678) = -rxt(k,213)*y(k,122)
         mat(k,555) = -rxt(k,215)*y(k,122)
         mat(k,569) = -rxt(k,216)*y(k,122)
         mat(k,540) = -rxt(k,217)*y(k,122)
         mat(k,631) = -rxt(k,218)*y(k,122)
         mat(k,933) = -rxt(k,219)*y(k,122)
         mat(k,1064) = -rxt(k,308)*y(k,122)
         mat(k,1768) = -rxt(k,309)*y(k,122)
         mat(k,906) = -rxt(k,310)*y(k,122)
         mat(k,990) = -rxt(k,311)*y(k,122)
         mat(k,2073) = -rxt(k,366)*y(k,122)
         mat(k,1157) = rxt(k,365)*y(k,139)
         mat(k,268) = rxt(k,307)*y(k,139)
         mat(k,261) = rxt(k,313)*y(k,139)
         mat(k,2073) = mat(k,2073) + rxt(k,365)*y(k,118) + rxt(k,307)*y(k,121) &
                      + rxt(k,313)*y(k,123)
      end do
      end subroutine nlnmat08
      subroutine nlnmat09( avec_len, mat, y, rxt )
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
         mat(k,256) = -(rxt(k,312)*y(k,32) + rxt(k,313)*y(k,139) + rxt(k,384)*y(k,98))
         mat(k,1264) = -rxt(k,312)*y(k,123)
         mat(k,2039) = -rxt(k,313)*y(k,123)
         mat(k,1045) = -rxt(k,384)*y(k,123)
         mat(k,317) = -(rxt(k,588)*y(k,70) + (rxt(k,589) + rxt(k,590)) * y(k,72))
         mat(k,1357) = -rxt(k,588)*y(k,124)
         mat(k,1412) = -(rxt(k,589) + rxt(k,590)) * y(k,124)
         mat(k,301) = rxt(k,596)*y(k,136)
         mat(k,387) = rxt(k,596)*y(k,109)
         mat(k,377) = -(rxt(k,321)*y(k,127) + rxt(k,339)*y(k,132) + rxt(k,361) &
                      *y(k,135) + rxt(k,396)*y(k,70) + rxt(k,397)*y(k,72))
         mat(k,1545) = -rxt(k,321)*y(k,125)
         mat(k,460) = -rxt(k,339)*y(k,125)
         mat(k,1867) = -rxt(k,361)*y(k,125)
         mat(k,1363) = -rxt(k,396)*y(k,125)
         mat(k,1414) = -rxt(k,397)*y(k,125)
         mat(k,1266) = rxt(k,398)*y(k,126)
         mat(k,1363) = mat(k,1363) + rxt(k,400)*y(k,126)
         mat(k,1414) = mat(k,1414) + rxt(k,401)*y(k,126)
         mat(k,1906) = rxt(k,402)*y(k,126)
         mat(k,132) = rxt(k,398)*y(k,32) + rxt(k,400)*y(k,70) + rxt(k,401)*y(k,72) &
                      + rxt(k,402)*y(k,73)
         mat(k,131) = -(rxt(k,398)*y(k,32) + rxt(k,400)*y(k,70) + rxt(k,401)*y(k,72) &
                      + rxt(k,402)*y(k,73))
         mat(k,1261) = -rxt(k,398)*y(k,126)
         mat(k,1345) = -rxt(k,400)*y(k,126)
         mat(k,1405) = -rxt(k,401)*y(k,126)
         mat(k,1903) = -rxt(k,402)*y(k,126)
         mat(k,773) = rxt(k,392)*y(k,72)
         mat(k,1405) = mat(k,1405) + rxt(k,392)*y(k,71)
         mat(k,1566) = -(rxt(k,150)*y(k,128) + rxt(k,173)*y(k,107) + rxt(k,190) &
                      *y(k,103) + rxt(k,208)*y(k,122) + rxt(k,214)*y(k,105) + rxt(k,226) &
                      *y(k,119) + rxt(k,243)*y(k,118) + rxt(k,314)*y(k,25) + rxt(k,315) &
                      *y(k,28) + rxt(k,316)*y(k,32) + rxt(k,317)*y(k,40) + rxt(k,318) &
                      *y(k,49) + rxt(k,319)*y(k,51) + rxt(k,320)*y(k,63) + rxt(k,321) &
                      *y(k,125) + rxt(k,322)*y(k,72) + rxt(k,323)*y(k,73) + (rxt(k,324) &
                      + rxt(k,325)) * y(k,70))
         mat(k,1611) = -rxt(k,150)*y(k,127)
         mat(k,1654) = -rxt(k,173)*y(k,127)
         mat(k,1248) = -rxt(k,190)*y(k,127)
         mat(k,1332) = -rxt(k,208)*y(k,127)
         mat(k,1697) = -rxt(k,214)*y(k,127)
         mat(k,1205) = -rxt(k,226)*y(k,127)
         mat(k,1162) = -rxt(k,243)*y(k,127)
         mat(k,1489) = -rxt(k,314)*y(k,127)
         mat(k,1114) = -rxt(k,315)*y(k,127)
         mat(k,1290) = -rxt(k,316)*y(k,127)
         mat(k,1773) = -rxt(k,317)*y(k,127)
         mat(k,2018) = -rxt(k,318)*y(k,127)
         mat(k,1530) = -rxt(k,319)*y(k,127)
         mat(k,1821) = -rxt(k,320)*y(k,127)
         mat(k,381) = -rxt(k,321)*y(k,127)
         mat(k,1441) = -rxt(k,322)*y(k,127)
         mat(k,1931) = -rxt(k,323)*y(k,127)
         mat(k,1392) = -(rxt(k,324) + rxt(k,325)) * y(k,127)
         mat(k,1392) = mat(k,1392) + rxt(k,125)*y(k,94) + rxt(k,334)*y(k,130)
         mat(k,1441) = mat(k,1441) + (rxt(k,133)+rxt(k,135))*y(k,98)
         mat(k,1026) = rxt(k,125)*y(k,70)
         mat(k,1069) = (rxt(k,133)+rxt(k,135))*y(k,72)
         mat(k,602) = rxt(k,334)*y(k,70)
         mat(k,1612) = -(rxt(k,148)*y(k,113) + rxt(k,149)*y(k,100) + rxt(k,150) &
                      *y(k,127) + rxt(k,151)*y(k,97) + rxt(k,152)*y(k,115) + rxt(k,153) &
                      *y(k,95) + rxt(k,154)*y(k,90) + rxt(k,155)*y(k,96) + rxt(k,156) &
                      *y(k,112) + rxt(k,157)*y(k,116) + rxt(k,159)*y(k,91) + rxt(k,160) &
                      *y(k,114) + rxt(k,161)*y(k,111) + rxt(k,255)*y(k,117) + rxt(k,256) &
                      *y(k,94) + rxt(k,257)*y(k,92) + rxt(k,329)*y(k,139) + rxt(k,364) &
                      *y(k,72) + rxt(k,583)*y(k,98) + rxt(k,591)*y(k,56) + rxt(k,593) &
                      *y(k,62))
         mat(k,877) = -rxt(k,148)*y(k,128)
         mat(k,498) = -rxt(k,149)*y(k,128)
         mat(k,1567) = -rxt(k,150)*y(k,128)
         mat(k,851) = -rxt(k,151)*y(k,128)
         mat(k,700) = -rxt(k,152)*y(k,128)
         mat(k,513) = -rxt(k,153)*y(k,128)
         mat(k,822) = -rxt(k,154)*y(k,128)
         mat(k,681) = -rxt(k,155)*y(k,128)
         mat(k,557) = -rxt(k,156)*y(k,128)
         mat(k,572) = -rxt(k,157)*y(k,128)
         mat(k,543) = -rxt(k,159)*y(k,128)
         mat(k,635) = -rxt(k,160)*y(k,128)
         mat(k,938) = -rxt(k,161)*y(k,128)
         mat(k,617) = -rxt(k,255)*y(k,128)
         mat(k,1027) = -rxt(k,256)*y(k,128)
         mat(k,527) = -rxt(k,257)*y(k,128)
         mat(k,2079) = -rxt(k,329)*y(k,128)
         mat(k,1442) = -rxt(k,364)*y(k,128)
         mat(k,1070) = -rxt(k,583)*y(k,128)
         mat(k,455) = -rxt(k,591)*y(k,128)
         mat(k,1974) = -rxt(k,593)*y(k,128)
         mat(k,1291) = rxt(k,594)*y(k,136)
         mat(k,1393) = rxt(k,338)*y(k,132)
         mat(k,1442) = mat(k,1442) + rxt(k,585)*y(k,110) + rxt(k,589)*y(k,124) &
                      + rxt(k,597)*y(k,136) + rxt(k,601)*y(k,137)
         mat(k,238) = rxt(k,585)*y(k,72)
         mat(k,323) = rxt(k,589)*y(k,72)
         mat(k,382) = rxt(k,339)*y(k,132)
         mat(k,1567) = mat(k,1567) + 2.000_r8*rxt(k,150)*y(k,128)
         mat(k,1612) = mat(k,1612) + 2.000_r8*rxt(k,150)*y(k,127)
         mat(k,468) = rxt(k,338)*y(k,70) + rxt(k,339)*y(k,125)
         mat(k,397) = rxt(k,594)*y(k,32) + rxt(k,597)*y(k,72)
         mat(k,205) = rxt(k,601)*y(k,72)
         mat(k,292) = -(rxt(k,326)*y(k,98) + (rxt(k,327) + rxt(k,328)) * y(k,139))
         mat(k,1049) = -rxt(k,326)*y(k,129)
         mat(k,2043) = -(rxt(k,327) + rxt(k,328)) * y(k,129)
         mat(k,1579) = rxt(k,329)*y(k,139)
         mat(k,459) = rxt(k,337)*y(k,139)
         mat(k,2043) = mat(k,2043) + rxt(k,329)*y(k,128) + rxt(k,337)*y(k,132)
         mat(k,594) = -((rxt(k,272) + rxt(k,279)) * y(k,62) + (rxt(k,277) + rxt(k,278) &
                      ) * y(k,63) + rxt(k,330)*y(k,32) + rxt(k,331)*y(k,40) + rxt(k,332) &
                      *y(k,73) + (rxt(k,333) + rxt(k,334)) * y(k,70))
         mat(k,1950) = -(rxt(k,272) + rxt(k,279)) * y(k,130)
         mat(k,1796) = -(rxt(k,277) + rxt(k,278)) * y(k,130)
         mat(k,1271) = -rxt(k,330)*y(k,130)
         mat(k,1751) = -rxt(k,331)*y(k,130)
         mat(k,1910) = -rxt(k,332)*y(k,130)
         mat(k,1369) = -(rxt(k,333) + rxt(k,334)) * y(k,130)
         mat(k,1369) = mat(k,1369) + rxt(k,336)*y(k,131)
         mat(k,1420) = rxt(k,126)*y(k,94) + rxt(k,362)*y(k,135)
         mat(k,1910) = mat(k,1910) + rxt(k,132)*y(k,97) + rxt(k,323)*y(k,127) &
                      + rxt(k,348)*y(k,134) + rxt(k,363)*y(k,135)
         mat(k,1007) = rxt(k,126)*y(k,72)
         mat(k,831) = rxt(k,132)*y(k,73)
         mat(k,1547) = rxt(k,323)*y(k,73)
         mat(k,240) = rxt(k,336)*y(k,70)
         mat(k,1836) = rxt(k,348)*y(k,73)
         mat(k,1869) = rxt(k,362)*y(k,72) + rxt(k,363)*y(k,73)
         mat(k,239) = -(rxt(k,335)*y(k,32) + rxt(k,336)*y(k,70))
         mat(k,1262) = -rxt(k,335)*y(k,131)
         mat(k,1352) = -rxt(k,336)*y(k,131)
         mat(k,1409) = rxt(k,322)*y(k,127)
         mat(k,1543) = rxt(k,322)*y(k,72)
         mat(k,461) = -(rxt(k,337)*y(k,139) + rxt(k,338)*y(k,70) + rxt(k,339)*y(k,125) &
                      + rxt(k,380)*y(k,98))
         mat(k,2048) = -rxt(k,337)*y(k,132)
         mat(k,1366) = -rxt(k,338)*y(k,132)
         mat(k,378) = -rxt(k,339)*y(k,132)
         mat(k,1054) = -rxt(k,380)*y(k,132)
         mat(k,1417) = rxt(k,364)*y(k,128)
         mat(k,1581) = rxt(k,364)*y(k,72)
         mat(k,984) = -(rxt(k,311)*y(k,122) + rxt(k,340)*y(k,55) + rxt(k,349)*y(k,62) &
                      + rxt(k,415)*y(k,41) + rxt(k,416)*y(k,43) + rxt(k,417)*y(k,101) &
                      + rxt(k,418)*y(k,70) + rxt(k,419)*y(k,73) + (4._r8*rxt(k,420) &
                      + 4._r8*rxt(k,421)) * y(k,133) + rxt(k,423)*y(k,52) + rxt(k,437) &
                      *y(k,64) + rxt(k,438)*y(k,56) + rxt(k,446)*y(k,63) + rxt(k,447) &
                      *y(k,51) + rxt(k,466)*y(k,29) + (rxt(k,468) + rxt(k,469) &
                      ) * y(k,28) + rxt(k,471)*y(k,49) + rxt(k,474)*y(k,54) + rxt(k,498) &
                      *y(k,5) + rxt(k,500)*y(k,45) + rxt(k,514)*y(k,16) + rxt(k,516) &
                      *y(k,18) + rxt(k,517)*y(k,19) + rxt(k,519)*y(k,20) + rxt(k,521) &
                      *y(k,24) + rxt(k,522)*y(k,46) + rxt(k,523)*y(k,47) + rxt(k,524) &
                      *y(k,48) + rxt(k,532)*y(k,17) + rxt(k,535)*y(k,22) + rxt(k,536) &
                      *y(k,23) + rxt(k,540)*y(k,31) + (rxt(k,542) + rxt(k,556) &
                      ) * y(k,35) + rxt(k,544)*y(k,75) + rxt(k,546)*y(k,81) + rxt(k,550) &
                      *y(k,78) + rxt(k,555)*y(k,80))
         mat(k,1319) = -rxt(k,311)*y(k,133)
         mat(k,279) = -rxt(k,340)*y(k,133)
         mat(k,1960) = -rxt(k,349)*y(k,133)
         mat(k,651) = -rxt(k,415)*y(k,133)
         mat(k,216) = -rxt(k,416)*y(k,133)
         mat(k,901) = -rxt(k,417)*y(k,133)
         mat(k,1379) = -rxt(k,418)*y(k,133)
         mat(k,1918) = -rxt(k,419)*y(k,133)
         mat(k,181) = -rxt(k,423)*y(k,133)
         mat(k,1725) = -rxt(k,437)*y(k,133)
         mat(k,449) = -rxt(k,438)*y(k,133)
         mat(k,1808) = -rxt(k,446)*y(k,133)
         mat(k,1517) = -rxt(k,447)*y(k,133)
         mat(k,360) = -rxt(k,466)*y(k,133)
         mat(k,1101) = -(rxt(k,468) + rxt(k,469)) * y(k,133)
         mat(k,2005) = -rxt(k,471)*y(k,133)
         mat(k,334) = -rxt(k,474)*y(k,133)
         mat(k,719) = -rxt(k,498)*y(k,133)
         mat(k,328) = -rxt(k,500)*y(k,133)
         mat(k,189) = -rxt(k,514)*y(k,133)
         mat(k,226) = -rxt(k,516)*y(k,133)
         mat(k,49) = -rxt(k,517)*y(k,133)
         mat(k,151) = -rxt(k,519)*y(k,133)
         mat(k,164) = -rxt(k,521)*y(k,133)
         mat(k,83) = -rxt(k,522)*y(k,133)
         mat(k,93) = -rxt(k,523)*y(k,133)
         mat(k,119) = -rxt(k,524)*y(k,133)
         mat(k,741) = -rxt(k,532)*y(k,133)
         mat(k,158) = -rxt(k,535)*y(k,133)
         mat(k,422) = -rxt(k,536)*y(k,133)
         mat(k,254) = -rxt(k,540)*y(k,133)
         mat(k,128) = -(rxt(k,542) + rxt(k,556)) * y(k,133)
         mat(k,144) = -rxt(k,544)*y(k,133)
         mat(k,354) = -rxt(k,546)*y(k,133)
         mat(k,273) = -rxt(k,550)*y(k,133)
         mat(k,480) = -rxt(k,555)*y(k,133)
         mat(k,741) = mat(k,741) + rxt(k,531)*y(k,70)
         mat(k,158) = mat(k,158) + .300_r8*rxt(k,535)*y(k,133)
         mat(k,422) = mat(k,422) + rxt(k,537)*y(k,71)
         mat(k,1476) = rxt(k,455)*y(k,101) + rxt(k,341)*y(k,134)
         mat(k,1101) = mat(k,1101) + rxt(k,342)*y(k,134)
         mat(k,345) = rxt(k,509)*y(k,139)
         mat(k,1760) = rxt(k,414)*y(k,73) + rxt(k,129)*y(k,97) + 2.000_r8*rxt(k,409) &
                      *y(k,101)
         mat(k,651) = mat(k,651) + rxt(k,406)*y(k,70) + rxt(k,389)*y(k,71)
         mat(k,216) = mat(k,216) + rxt(k,407)*y(k,70)
         mat(k,328) = mat(k,328) + rxt(k,499)*y(k,70) + rxt(k,505)*y(k,71)
         mat(k,2005) = mat(k,2005) + rxt(k,470)*y(k,70) + rxt(k,482)*y(k,71) &
                      + rxt(k,356)*y(k,135)
         mat(k,1517) = mat(k,1517) + rxt(k,124)*y(k,94) + rxt(k,357)*y(k,135)
         mat(k,312) = rxt(k,501)*y(k,70)
         mat(k,334) = mat(k,334) + rxt(k,473)*y(k,70)
         mat(k,1960) = mat(k,1960) + rxt(k,439)*y(k,101)
         mat(k,1808) = mat(k,1808) + rxt(k,346)*y(k,134)
         mat(k,1725) = mat(k,1725) + rxt(k,434)*y(k,101)
         mat(k,1379) = mat(k,1379) + rxt(k,531)*y(k,17) + rxt(k,406)*y(k,41) &
                      + rxt(k,407)*y(k,43) + rxt(k,499)*y(k,45) + rxt(k,470)*y(k,49) &
                      + rxt(k,501)*y(k,53) + rxt(k,473)*y(k,54) + rxt(k,412)*y(k,101)
         mat(k,789) = rxt(k,537)*y(k,23) + rxt(k,389)*y(k,41) + rxt(k,505)*y(k,45) &
                      + rxt(k,482)*y(k,49) + 2.000_r8*rxt(k,390)*y(k,139)
         mat(k,1918) = mat(k,1918) + rxt(k,414)*y(k,40) + rxt(k,413)*y(k,101) &
                      + rxt(k,348)*y(k,134)
         mat(k,1013) = rxt(k,124)*y(k,51)
         mat(k,837) = rxt(k,129)*y(k,40)
         mat(k,1057) = rxt(k,138)*y(k,99)
         mat(k,169) = rxt(k,138)*y(k,98) + rxt(k,139)*y(k,139)
         mat(k,491) = rxt(k,189)*y(k,103) + rxt(k,203)*y(k,105) + rxt(k,172)*y(k,107) &
                      + rxt(k,242)*y(k,118) + rxt(k,224)*y(k,119) + rxt(k,207) &
                      *y(k,122) + rxt(k,149)*y(k,128)
         mat(k,901) = mat(k,901) + rxt(k,455)*y(k,25) + 2.000_r8*rxt(k,409)*y(k,40) &
                      + rxt(k,439)*y(k,62) + rxt(k,434)*y(k,64) + rxt(k,412)*y(k,70) &
                      + rxt(k,413)*y(k,73)
         mat(k,1235) = rxt(k,189)*y(k,100)
         mat(k,1684) = rxt(k,203)*y(k,100)
         mat(k,1641) = rxt(k,172)*y(k,100)
         mat(k,1149) = rxt(k,242)*y(k,100)
         mat(k,1192) = rxt(k,224)*y(k,100)
         mat(k,1319) = mat(k,1319) + rxt(k,207)*y(k,100)
         mat(k,1598) = rxt(k,149)*y(k,100)
         mat(k,293) = rxt(k,328)*y(k,139)
         mat(k,984) = mat(k,984) + .300_r8*rxt(k,535)*y(k,22)
         mat(k,1841) = rxt(k,341)*y(k,25) + rxt(k,342)*y(k,28) + rxt(k,346)*y(k,63) &
                      + rxt(k,348)*y(k,73)
         mat(k,1877) = rxt(k,356)*y(k,49) + rxt(k,357)*y(k,51) + rxt(k,355)*y(k,139)
         mat(k,2065) = rxt(k,509)*y(k,39) + 2.000_r8*rxt(k,390)*y(k,71) + rxt(k,139) &
                      *y(k,99) + rxt(k,328)*y(k,129) + rxt(k,355)*y(k,135)
      end do
      end subroutine nlnmat09
      subroutine nlnmat10( avec_len, mat, y, rxt )
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
         mat(k,1861) = -(rxt(k,341)*y(k,25) + rxt(k,342)*y(k,28) + rxt(k,343)*y(k,32) &
                      + rxt(k,344)*y(k,40) + rxt(k,345)*y(k,49) + rxt(k,346)*y(k,63) &
                      + rxt(k,347)*y(k,70) + rxt(k,348)*y(k,73))
         mat(k,1496) = -rxt(k,341)*y(k,134)
         mat(k,1121) = -rxt(k,342)*y(k,134)
         mat(k,1297) = -rxt(k,343)*y(k,134)
         mat(k,1780) = -rxt(k,344)*y(k,134)
         mat(k,2025) = -rxt(k,345)*y(k,134)
         mat(k,1828) = -rxt(k,346)*y(k,134)
         mat(k,1399) = -rxt(k,347)*y(k,134)
         mat(k,1938) = -rxt(k,348)*y(k,134)
         mat(k,1780) = mat(k,1780) + rxt(k,117)*y(k,94) + rxt(k,287)*y(k,111) &
                      + rxt(k,331)*y(k,130)
         mat(k,662) = rxt(k,354)*y(k,135)
         mat(k,1033) = rxt(k,117)*y(k,40)
         mat(k,944) = rxt(k,287)*y(k,40)
         mat(k,605) = rxt(k,331)*y(k,40)
         mat(k,1897) = rxt(k,354)*y(k,41) + rxt(k,355)*y(k,139)
         mat(k,2085) = rxt(k,355)*y(k,135)
         mat(k,1898) = -(rxt(k,136)*y(k,62) + rxt(k,350)*y(k,25) + rxt(k,351)*y(k,28) &
                      + rxt(k,352)*y(k,32) + (rxt(k,353) + rxt(k,354)) * y(k,41) &
                      + rxt(k,355)*y(k,139) + rxt(k,356)*y(k,49) + rxt(k,357)*y(k,51) &
                      + rxt(k,359)*y(k,63) + rxt(k,360)*y(k,70) + rxt(k,361)*y(k,125) &
                      + rxt(k,362)*y(k,72) + rxt(k,363)*y(k,73))
         mat(k,1981) = -rxt(k,136)*y(k,135)
         mat(k,1497) = -rxt(k,350)*y(k,135)
         mat(k,1122) = -rxt(k,351)*y(k,135)
         mat(k,1298) = -rxt(k,352)*y(k,135)
         mat(k,663) = -(rxt(k,353) + rxt(k,354)) * y(k,135)
         mat(k,2086) = -rxt(k,355)*y(k,135)
         mat(k,2026) = -rxt(k,356)*y(k,135)
         mat(k,1538) = -rxt(k,357)*y(k,135)
         mat(k,1829) = -rxt(k,359)*y(k,135)
         mat(k,1400) = -rxt(k,360)*y(k,135)
         mat(k,383) = -rxt(k,361)*y(k,135)
         mat(k,1449) = -rxt(k,362)*y(k,135)
         mat(k,1939) = -rxt(k,363)*y(k,135)
         mat(k,1400) = mat(k,1400) + rxt(k,325)*y(k,127)
         mat(k,1939) = mat(k,1939) + rxt(k,134)*y(k,98)
         mat(k,1077) = rxt(k,134)*y(k,73)
         mat(k,1574) = rxt(k,325)*y(k,70)
         mat(k,388) = -(rxt(k,594)*y(k,32) + rxt(k,596)*y(k,109) + rxt(k,597)*y(k,72))
         mat(k,1267) = -rxt(k,594)*y(k,136)
         mat(k,302) = -rxt(k,596)*y(k,136)
         mat(k,1415) = -rxt(k,597)*y(k,136)
         mat(k,1364) = rxt(k,587)*y(k,110) + rxt(k,588)*y(k,124) + rxt(k,600)*y(k,137) &
                      + rxt(k,606)*y(k,138)
         mat(k,1052) = rxt(k,598)*y(k,137) + rxt(k,603)*y(k,138)
         mat(k,232) = rxt(k,587)*y(k,70)
         mat(k,318) = rxt(k,588)*y(k,70)
         mat(k,202) = rxt(k,600)*y(k,70) + rxt(k,598)*y(k,98)
         mat(k,197) = rxt(k,606)*y(k,70) + rxt(k,603)*y(k,98)
         mat(k,200) = -(rxt(k,598)*y(k,98) + rxt(k,600)*y(k,70) + rxt(k,601)*y(k,72))
         mat(k,1042) = -rxt(k,598)*y(k,137)
         mat(k,1348) = -rxt(k,600)*y(k,137)
         mat(k,1407) = -rxt(k,601)*y(k,137)
         mat(k,1042) = mat(k,1042) + rxt(k,602)*y(k,138)
         mat(k,194) = rxt(k,602)*y(k,98)
         mat(k,193) = -((rxt(k,602) + rxt(k,603)) * y(k,98) + rxt(k,606)*y(k,70))
         mat(k,1041) = -(rxt(k,602) + rxt(k,603)) * y(k,138)
         mat(k,1347) = -rxt(k,606)*y(k,138)
         mat(k,2090) = -(rxt(k,107)*y(k,90) + rxt(k,118)*y(k,96) + rxt(k,119)*y(k,94) &
                      + rxt(k,139)*y(k,99) + rxt(k,140)*y(k,104) + rxt(k,143)*y(k,106) &
                      + rxt(k,288)*y(k,111) + rxt(k,294)*y(k,115) + rxt(k,296) &
                      *y(k,113) + rxt(k,305)*y(k,120) + rxt(k,307)*y(k,121) + rxt(k,313) &
                      *y(k,123) + (rxt(k,327) + rxt(k,328)) * y(k,129) + rxt(k,329) &
                      *y(k,128) + rxt(k,337)*y(k,132) + rxt(k,355)*y(k,135) + rxt(k,365) &
                      *y(k,118) + rxt(k,366)*y(k,122) + rxt(k,367)*y(k,119) + rxt(k,372) &
                      *y(k,108) + rxt(k,374)*y(k,102) + rxt(k,376)*y(k,103) + rxt(k,378) &
                      *y(k,105) + rxt(k,390)*y(k,71) + rxt(k,509)*y(k,39) + rxt(k,557) &
                      *y(k,82))
         mat(k,829) = -rxt(k,107)*y(k,139)
         mat(k,687) = -rxt(k,118)*y(k,139)
         mat(k,1038) = -rxt(k,119)*y(k,139)
         mat(k,172) = -rxt(k,139)*y(k,139)
         mat(k,110) = -rxt(k,140)*y(k,139)
         mat(k,114) = -rxt(k,143)*y(k,139)
         mat(k,948) = -rxt(k,288)*y(k,139)
         mat(k,707) = -rxt(k,294)*y(k,139)
         mat(k,886) = -rxt(k,296)*y(k,139)
         mat(k,178) = -rxt(k,305)*y(k,139)
         mat(k,270) = -rxt(k,307)*y(k,139)
         mat(k,263) = -rxt(k,313)*y(k,139)
         mat(k,299) = -(rxt(k,327) + rxt(k,328)) * y(k,139)
         mat(k,1623) = -rxt(k,329)*y(k,139)
         mat(k,472) = -rxt(k,337)*y(k,139)
         mat(k,1902) = -rxt(k,355)*y(k,139)
         mat(k,1174) = -rxt(k,365)*y(k,139)
         mat(k,1344) = -rxt(k,366)*y(k,139)
         mat(k,1217) = -rxt(k,367)*y(k,139)
         mat(k,288) = -rxt(k,372)*y(k,139)
         mat(k,251) = -rxt(k,374)*y(k,139)
         mat(k,1260) = -rxt(k,376)*y(k,139)
         mat(k,1709) = -rxt(k,378)*y(k,139)
         mat(k,808) = -rxt(k,390)*y(k,139)
         mat(k,349) = -rxt(k,509)*y(k,139)
         mat(k,106) = -rxt(k,557)*y(k,139)
         mat(k,192) = rxt(k,514)*y(k,133)
         mat(k,757) = rxt(k,532)*y(k,133)
         mat(k,229) = rxt(k,516)*y(k,133)
         mat(k,51) = rxt(k,517)*y(k,133)
         mat(k,154) = rxt(k,519)*y(k,133)
         mat(k,160) = rxt(k,535)*y(k,133)
         mat(k,430) = rxt(k,536)*y(k,133)
         mat(k,1785) = rxt(k,410)*y(k,101) + rxt(k,344)*y(k,134)
         mat(k,666) = rxt(k,415)*y(k,133) + rxt(k,353)*y(k,135)
         mat(k,220) = rxt(k,416)*y(k,133)
         mat(k,331) = rxt(k,500)*y(k,133)
         mat(k,121) = rxt(k,524)*y(k,133)
         mat(k,2030) = (rxt(k,570)+rxt(k,575))*y(k,53) + (rxt(k,563)+rxt(k,569) &
                       +rxt(k,574))*y(k,54) + rxt(k,106)*y(k,91) + rxt(k,471)*y(k,133) &
                      + rxt(k,345)*y(k,134)
         mat(k,1542) = rxt(k,295)*y(k,115) + rxt(k,447)*y(k,133)
         mat(k,185) = rxt(k,423)*y(k,133)
         mat(k,316) = (rxt(k,570)+rxt(k,575))*y(k,49)
         mat(k,339) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,49) + rxt(k,474)*y(k,133)
         mat(k,282) = rxt(k,340)*y(k,133)
         mat(k,412) = rxt(k,293)*y(k,114)
         mat(k,1985) = rxt(k,123)*y(k,96)
         mat(k,1833) = rxt(k,120)*y(k,96)
         mat(k,829) = mat(k,829) + 3.000_r8*rxt(k,195)*y(k,103) + 4.000_r8*rxt(k,147) &
                      *y(k,105) + 5.000_r8*rxt(k,177)*y(k,107) + 2.000_r8*rxt(k,230) &
                      *y(k,119) + rxt(k,212)*y(k,122)
         mat(k,549) = rxt(k,106)*y(k,49) + 4.000_r8*rxt(k,199)*y(k,103) &
                      + 5.000_r8*rxt(k,164)*y(k,105) + 6.000_r8*rxt(k,182)*y(k,107) &
                      + rxt(k,252)*y(k,118) + 3.000_r8*rxt(k,234)*y(k,119) &
                      + 2.000_r8*rxt(k,217)*y(k,122) + rxt(k,159)*y(k,128)
         mat(k,533) = 3.000_r8*rxt(k,187)*y(k,103) + (4.000_r8*rxt(k,181) &
                       +4.000_r8*rxt(k,267))*y(k,105) + (5.000_r8*rxt(k,169) &
                       +5.000_r8*rxt(k,268))*y(k,107) + 2.000_r8*rxt(k,222)*y(k,119) &
                      + rxt(k,205)*y(k,122)
         mat(k,1038) = mat(k,1038) + 3.000_r8*rxt(k,186)*y(k,103) + ( &
                      + 4.000_r8*rxt(k,170)+4.000_r8*rxt(k,258))*y(k,105) + ( &
                      + 5.000_r8*rxt(k,168)+5.000_r8*rxt(k,265))*y(k,107) &
                      + 2.000_r8*rxt(k,221)*y(k,119) + rxt(k,204)*y(k,122)
         mat(k,518) = 5.000_r8*rxt(k,194)*y(k,103) + (6.000_r8*rxt(k,247) &
                       +6.000_r8*rxt(k,271))*y(k,105) + (7.000_r8*rxt(k,176) &
                       +7.000_r8*rxt(k,259))*y(k,107) + 2.000_r8*rxt(k,246)*y(k,118) &
                      + 4.000_r8*rxt(k,229)*y(k,119) + 3.000_r8*rxt(k,211)*y(k,122) &
                      + 2.000_r8*rxt(k,153)*y(k,128)
         mat(k,687) = mat(k,687) + rxt(k,123)*y(k,62) + rxt(k,120)*y(k,63) &
                      + 4.000_r8*rxt(k,196)*y(k,103) + (5.000_r8*rxt(k,158) &
                       +5.000_r8*rxt(k,260))*y(k,105) + (6.000_r8*rxt(k,178) &
                       +6.000_r8*rxt(k,261))*y(k,107) + rxt(k,249)*y(k,118) &
                      + 3.000_r8*rxt(k,231)*y(k,119) + 2.000_r8*rxt(k,213)*y(k,122) &
                      + rxt(k,155)*y(k,128)
         mat(k,861) = 3.000_r8*rxt(k,191)*y(k,103) + 4.000_r8*rxt(k,225)*y(k,105) &
                      + 5.000_r8*rxt(k,174)*y(k,107) + 2.000_r8*rxt(k,227)*y(k,119) &
                      + rxt(k,209)*y(k,122)
         mat(k,1081) = rxt(k,138)*y(k,99) + 2.000_r8*rxt(k,382)*y(k,102) &
                      + 3.000_r8*rxt(k,383)*y(k,103) + 4.000_r8*rxt(k,141)*y(k,105) &
                      + 5.000_r8*rxt(k,144)*y(k,107) + rxt(k,381)*y(k,108) &
                      + 2.000_r8*rxt(k,303)*y(k,119) + 3.000_r8*rxt(k,304)*y(k,120) &
                      + rxt(k,308)*y(k,122) + rxt(k,326)*y(k,129)
         mat(k,172) = mat(k,172) + rxt(k,138)*y(k,98)
         mat(k,503) = 3.000_r8*rxt(k,189)*y(k,103) + 4.000_r8*rxt(k,203)*y(k,105) &
                      + 5.000_r8*rxt(k,172)*y(k,107) + 2.000_r8*rxt(k,224)*y(k,119) &
                      + rxt(k,207)*y(k,122)
         mat(k,921) = rxt(k,410)*y(k,40) + rxt(k,417)*y(k,133)
         mat(k,251) = mat(k,251) + 2.000_r8*rxt(k,382)*y(k,98)
         mat(k,1260) = mat(k,1260) + 3.000_r8*rxt(k,195)*y(k,90) + 4.000_r8*rxt(k,199) &
                      *y(k,91) + 3.000_r8*rxt(k,187)*y(k,92) + 3.000_r8*rxt(k,186) &
                      *y(k,94) + 5.000_r8*rxt(k,194)*y(k,95) + 4.000_r8*rxt(k,196) &
                      *y(k,96) + 3.000_r8*rxt(k,191)*y(k,97) + 3.000_r8*rxt(k,383) &
                      *y(k,98) + 3.000_r8*rxt(k,189)*y(k,100) + 3.000_r8*rxt(k,201) &
                      *y(k,111) + 4.000_r8*rxt(k,197)*y(k,112) + 3.000_r8*rxt(k,188) &
                      *y(k,113) + 5.000_r8*rxt(k,200)*y(k,114) + 4.000_r8*rxt(k,193) &
                      *y(k,115) + 3.000_r8*rxt(k,198)*y(k,116) + 3.000_r8*rxt(k,185) &
                      *y(k,117) + 3.000_r8*rxt(k,190)*y(k,127)
         mat(k,1709) = mat(k,1709) + 4.000_r8*rxt(k,147)*y(k,90) + 5.000_r8*rxt(k,164) &
                      *y(k,91) + (4.000_r8*rxt(k,181)+4.000_r8*rxt(k,267))*y(k,92) + ( &
                      + 4.000_r8*rxt(k,170)+4.000_r8*rxt(k,258))*y(k,94) + ( &
                      + 6.000_r8*rxt(k,247)+6.000_r8*rxt(k,271))*y(k,95) + ( &
                      + 5.000_r8*rxt(k,158)+5.000_r8*rxt(k,260))*y(k,96) &
                      + 4.000_r8*rxt(k,225)*y(k,97) + 4.000_r8*rxt(k,141)*y(k,98) &
                      + 4.000_r8*rxt(k,203)*y(k,100) + 4.000_r8*rxt(k,166)*y(k,111) &
                      + 5.000_r8*rxt(k,162)*y(k,112) + (4.000_r8*rxt(k,192) &
                       +4.000_r8*rxt(k,264))*y(k,113) + 6.000_r8*rxt(k,165)*y(k,114) + ( &
                      + 5.000_r8*rxt(k,236)+5.000_r8*rxt(k,262))*y(k,115) &
                      + 4.000_r8*rxt(k,163)*y(k,116) + (4.000_r8*rxt(k,146) &
                       +4.000_r8*rxt(k,269))*y(k,117) + 4.000_r8*rxt(k,214)*y(k,127)
         mat(k,1666) = 5.000_r8*rxt(k,177)*y(k,90) + 6.000_r8*rxt(k,182)*y(k,91) + ( &
                      + 5.000_r8*rxt(k,169)+5.000_r8*rxt(k,268))*y(k,92) + ( &
                      + 5.000_r8*rxt(k,168)+5.000_r8*rxt(k,265))*y(k,94) + ( &
                      + 7.000_r8*rxt(k,176)+7.000_r8*rxt(k,259))*y(k,95) + ( &
                      + 6.000_r8*rxt(k,178)+6.000_r8*rxt(k,261))*y(k,96) &
                      + 5.000_r8*rxt(k,174)*y(k,97) + 5.000_r8*rxt(k,144)*y(k,98) &
                      + 5.000_r8*rxt(k,172)*y(k,100) + 5.000_r8*rxt(k,184)*y(k,111) &
                      + 6.000_r8*rxt(k,179)*y(k,112) + (5.000_r8*rxt(k,171) &
                       +5.000_r8*rxt(k,266))*y(k,113) + 7.000_r8*rxt(k,183)*y(k,114) + ( &
                      + 6.000_r8*rxt(k,175)+6.000_r8*rxt(k,263))*y(k,115) &
                      + 5.000_r8*rxt(k,180)*y(k,116) + (5.000_r8*rxt(k,167) &
                       +5.000_r8*rxt(k,270))*y(k,117) + 5.000_r8*rxt(k,173)*y(k,127)
         mat(k,288) = mat(k,288) + rxt(k,381)*y(k,98)
         mat(k,948) = mat(k,948) + 3.000_r8*rxt(k,201)*y(k,103) + 4.000_r8*rxt(k,166) &
                      *y(k,105) + 5.000_r8*rxt(k,184)*y(k,107) + 2.000_r8*rxt(k,237) &
                      *y(k,119) + rxt(k,219)*y(k,122)
         mat(k,563) = 4.000_r8*rxt(k,197)*y(k,103) + 5.000_r8*rxt(k,162)*y(k,105) &
                      + 6.000_r8*rxt(k,179)*y(k,107) + rxt(k,250)*y(k,118) &
                      + 3.000_r8*rxt(k,232)*y(k,119) + 2.000_r8*rxt(k,215)*y(k,122) &
                      + rxt(k,156)*y(k,128)
         mat(k,886) = mat(k,886) + 3.000_r8*rxt(k,188)*y(k,103) + (4.000_r8*rxt(k,192) &
                       +4.000_r8*rxt(k,264))*y(k,105) + (5.000_r8*rxt(k,171) &
                       +5.000_r8*rxt(k,266))*y(k,107) + 2.000_r8*rxt(k,223)*y(k,119) &
                      + rxt(k,206)*y(k,122)
         mat(k,642) = rxt(k,293)*y(k,58) + 5.000_r8*rxt(k,200)*y(k,103) &
                      + 6.000_r8*rxt(k,165)*y(k,105) + 7.000_r8*rxt(k,183)*y(k,107) &
                      + 2.000_r8*rxt(k,253)*y(k,118) + 4.000_r8*rxt(k,235)*y(k,119) &
                      + 3.000_r8*rxt(k,218)*y(k,122) + 2.000_r8*rxt(k,160)*y(k,128)
         mat(k,707) = mat(k,707) + rxt(k,295)*y(k,51) + 4.000_r8*rxt(k,193)*y(k,103) + ( &
                      + 5.000_r8*rxt(k,236)+5.000_r8*rxt(k,262))*y(k,105) + ( &
                      + 6.000_r8*rxt(k,175)+6.000_r8*rxt(k,263))*y(k,107) + rxt(k,245) &
                      *y(k,118) + 3.000_r8*rxt(k,228)*y(k,119) + 2.000_r8*rxt(k,210) &
                      *y(k,122) + rxt(k,152)*y(k,128)
         mat(k,579) = 3.000_r8*rxt(k,198)*y(k,103) + 4.000_r8*rxt(k,163)*y(k,105) &
                      + 5.000_r8*rxt(k,180)*y(k,107) + 2.000_r8*rxt(k,233)*y(k,119) &
                      + rxt(k,216)*y(k,122)
         mat(k,622) = 3.000_r8*rxt(k,185)*y(k,103) + (4.000_r8*rxt(k,146) &
                       +4.000_r8*rxt(k,269))*y(k,105) + (5.000_r8*rxt(k,167) &
                       +5.000_r8*rxt(k,270))*y(k,107) + 2.000_r8*rxt(k,220)*y(k,119) &
                      + rxt(k,202)*y(k,122)
         mat(k,1174) = mat(k,1174) + rxt(k,252)*y(k,91) + 2.000_r8*rxt(k,246)*y(k,95) &
                      + rxt(k,249)*y(k,96) + rxt(k,250)*y(k,112) + 2.000_r8*rxt(k,253) &
                      *y(k,114) + rxt(k,245)*y(k,115)
         mat(k,1217) = mat(k,1217) + 2.000_r8*rxt(k,230)*y(k,90) + 3.000_r8*rxt(k,234) &
                      *y(k,91) + 2.000_r8*rxt(k,222)*y(k,92) + 2.000_r8*rxt(k,221) &
                      *y(k,94) + 4.000_r8*rxt(k,229)*y(k,95) + 3.000_r8*rxt(k,231) &
                      *y(k,96) + 2.000_r8*rxt(k,227)*y(k,97) + 2.000_r8*rxt(k,303) &
                      *y(k,98) + 2.000_r8*rxt(k,224)*y(k,100) + 2.000_r8*rxt(k,237) &
                      *y(k,111) + 3.000_r8*rxt(k,232)*y(k,112) + 2.000_r8*rxt(k,223) &
                      *y(k,113) + 4.000_r8*rxt(k,235)*y(k,114) + 3.000_r8*rxt(k,228) &
                      *y(k,115) + 2.000_r8*rxt(k,233)*y(k,116) + 2.000_r8*rxt(k,220) &
                      *y(k,117) + 2.000_r8*rxt(k,226)*y(k,127)
         mat(k,178) = mat(k,178) + 3.000_r8*rxt(k,304)*y(k,98)
         mat(k,1344) = mat(k,1344) + rxt(k,212)*y(k,90) + 2.000_r8*rxt(k,217)*y(k,91) &
                      + rxt(k,205)*y(k,92) + rxt(k,204)*y(k,94) + 3.000_r8*rxt(k,211) &
                      *y(k,95) + 2.000_r8*rxt(k,213)*y(k,96) + rxt(k,209)*y(k,97) &
                      + rxt(k,308)*y(k,98) + rxt(k,207)*y(k,100) + rxt(k,219)*y(k,111) &
                      + 2.000_r8*rxt(k,215)*y(k,112) + rxt(k,206)*y(k,113) &
                      + 3.000_r8*rxt(k,218)*y(k,114) + 2.000_r8*rxt(k,210)*y(k,115) &
                      + rxt(k,216)*y(k,116) + rxt(k,202)*y(k,117) + rxt(k,208) &
                      *y(k,127)
         mat(k,1578) = 3.000_r8*rxt(k,190)*y(k,103) + 4.000_r8*rxt(k,214)*y(k,105) &
                      + 5.000_r8*rxt(k,173)*y(k,107) + 2.000_r8*rxt(k,226)*y(k,119) &
                      + rxt(k,208)*y(k,122)
         mat(k,1623) = mat(k,1623) + rxt(k,159)*y(k,91) + 2.000_r8*rxt(k,153)*y(k,95) &
                      + rxt(k,155)*y(k,96) + rxt(k,156)*y(k,112) + 2.000_r8*rxt(k,160) &
                      *y(k,114) + rxt(k,152)*y(k,115)
         mat(k,299) = mat(k,299) + rxt(k,326)*y(k,98)
         mat(k,1005) = rxt(k,514)*y(k,16) + rxt(k,532)*y(k,17) + rxt(k,516)*y(k,18) &
                      + rxt(k,517)*y(k,19) + rxt(k,519)*y(k,20) + rxt(k,535)*y(k,22) &
                      + rxt(k,536)*y(k,23) + rxt(k,415)*y(k,41) + rxt(k,416)*y(k,43) &
                      + rxt(k,500)*y(k,45) + rxt(k,524)*y(k,48) + rxt(k,471)*y(k,49) &
                      + rxt(k,447)*y(k,51) + rxt(k,423)*y(k,52) + rxt(k,474)*y(k,54) &
                      + rxt(k,340)*y(k,55) + rxt(k,417)*y(k,101) + 2.000_r8*rxt(k,420) &
                      *y(k,133)
         mat(k,1866) = rxt(k,344)*y(k,40) + rxt(k,345)*y(k,49)
         mat(k,1902) = mat(k,1902) + rxt(k,353)*y(k,41)
      end do
      end subroutine nlnmat10
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
         mat(k, 27) = lmat(k, 27)
         mat(k, 28) = lmat(k, 28)
         mat(k, 29) = mat(k, 29) + lmat(k, 29)
         mat(k, 31) = mat(k, 31) + lmat(k, 31)
         mat(k, 32) = mat(k, 32) + lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 34) = mat(k, 34) + lmat(k, 34)
         mat(k, 36) = mat(k, 36) + lmat(k, 36)
         mat(k, 37) = mat(k, 37) + lmat(k, 37)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 40) = mat(k, 40) + lmat(k, 40)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 43) = mat(k, 43) + lmat(k, 43)
         mat(k, 44) = mat(k, 44) + lmat(k, 44)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 69) = mat(k, 69) + lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 90) = mat(k, 90) + lmat(k, 90)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = lmat(k, 105)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = lmat(k, 122)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 141) = lmat(k, 141)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = lmat(k, 149)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 195) = lmat(k, 195)
         mat(k, 196) = lmat(k, 196)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = mat(k, 202) + lmat(k, 202)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = lmat(k, 211)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 259) = lmat(k, 259)
         mat(k, 264) = mat(k, 264) + lmat(k, 264)
         mat(k, 266) = lmat(k, 266)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 281) = lmat(k, 281)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 297) = lmat(k, 297)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 309) = mat(k, 309) + lmat(k, 309)
         mat(k, 310) = lmat(k, 310)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 352) = lmat(k, 352)
         mat(k, 355) = lmat(k, 355)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 380) = mat(k, 380) + lmat(k, 380)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 389) = lmat(k, 389)
         mat(k, 393) = lmat(k, 393)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 409) = lmat(k, 409)
         mat(k, 410) = lmat(k, 410)
         mat(k, 411) = lmat(k, 411)
         mat(k, 414) = lmat(k, 414)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 444) = lmat(k, 444)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 450) = lmat(k, 450)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 474) = lmat(k, 474)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 536) = lmat(k, 536)
         mat(k, 549) = mat(k, 549) + lmat(k, 549)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 551) = lmat(k, 551)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 594) = mat(k, 594) + lmat(k, 594)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 606) = lmat(k, 606)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 610) = lmat(k, 610)
         mat(k, 616) = mat(k, 616) + lmat(k, 616)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 626) = lmat(k, 626)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 647) = mat(k, 647) + lmat(k, 647)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 673) = lmat(k, 673)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 691) = mat(k, 691) + lmat(k, 691)
         mat(k, 692) = lmat(k, 692)
         mat(k, 707) = mat(k, 707) + lmat(k, 707)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 722) = mat(k, 722) + lmat(k, 722)
         mat(k, 733) = mat(k, 733) + lmat(k, 733)
         mat(k, 736) = lmat(k, 736)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 787) = mat(k, 787) + lmat(k, 787)
         mat(k, 794) = mat(k, 794) + lmat(k, 794)
         mat(k, 811) = mat(k, 811) + lmat(k, 811)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 844) = mat(k, 844) + lmat(k, 844)
         mat(k, 850) = lmat(k, 850)
         mat(k, 866) = mat(k, 866) + lmat(k, 866)
         mat(k, 868) = lmat(k, 868)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 926) = mat(k, 926) + lmat(k, 926)
         mat(k, 928) = lmat(k, 928)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 984) = mat(k, 984) + lmat(k, 984)
         mat(k,1014) = mat(k,1014) + lmat(k,1014)
         mat(k,1020) = mat(k,1020) + lmat(k,1020)
         mat(k,1034) = lmat(k,1034)
         mat(k,1058) = mat(k,1058) + lmat(k,1058)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1112) = mat(k,1112) + lmat(k,1112)
         mat(k,1127) = lmat(k,1127)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1197) = mat(k,1197) + lmat(k,1197)
         mat(k,1218) = lmat(k,1218)
         mat(k,1241) = mat(k,1241) + lmat(k,1241)
         mat(k,1260) = mat(k,1260) + lmat(k,1260)
         mat(k,1263) = mat(k,1263) + lmat(k,1263)
         mat(k,1284) = mat(k,1284) + lmat(k,1284)
         mat(k,1286) = lmat(k,1286)
         mat(k,1327) = mat(k,1327) + lmat(k,1327)
         mat(k,1347) = mat(k,1347) + lmat(k,1347)
         mat(k,1348) = mat(k,1348) + lmat(k,1348)
         mat(k,1364) = mat(k,1364) + lmat(k,1364)
         mat(k,1381) = mat(k,1381) + lmat(k,1381)
         mat(k,1388) = mat(k,1388) + lmat(k,1388)
         mat(k,1406) = lmat(k,1406)
         mat(k,1407) = mat(k,1407) + lmat(k,1407)
         mat(k,1415) = mat(k,1415) + lmat(k,1415)
         mat(k,1423) = mat(k,1423) + lmat(k,1423)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1437) = mat(k,1437) + lmat(k,1437)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1487) = mat(k,1487) + lmat(k,1487)
         mat(k,1517) = mat(k,1517) + lmat(k,1517)
         mat(k,1529) = mat(k,1529) + lmat(k,1529)
         mat(k,1536) = lmat(k,1536)
         mat(k,1555) = mat(k,1555) + lmat(k,1555)
         mat(k,1563) = mat(k,1563) + lmat(k,1563)
         mat(k,1566) = mat(k,1566) + lmat(k,1566)
         mat(k,1602) = mat(k,1602) + lmat(k,1602)
         mat(k,1612) = mat(k,1612) + lmat(k,1612)
         mat(k,1621) = mat(k,1621) + lmat(k,1621)
         mat(k,1656) = mat(k,1656) + lmat(k,1656)
         mat(k,1657) = lmat(k,1657)
         mat(k,1666) = mat(k,1666) + lmat(k,1666)
         mat(k,1690) = lmat(k,1690)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1709) = mat(k,1709) + lmat(k,1709)
         mat(k,1733) = mat(k,1733) + lmat(k,1733)
         mat(k,1734) = mat(k,1734) + lmat(k,1734)
         mat(k,1736) = mat(k,1736) + lmat(k,1736)
         mat(k,1741) = mat(k,1741) + lmat(k,1741)
         mat(k,1743) = mat(k,1743) + lmat(k,1743)
         mat(k,1747) = mat(k,1747) + lmat(k,1747)
         mat(k,1778) = mat(k,1778) + lmat(k,1778)
         mat(k,1808) = mat(k,1808) + lmat(k,1808)
         mat(k,1817) = mat(k,1817) + lmat(k,1817)
         mat(k,1820) = mat(k,1820) + lmat(k,1820)
         mat(k,1827) = mat(k,1827) + lmat(k,1827)
         mat(k,1831) = mat(k,1831) + lmat(k,1831)
         mat(k,1841) = mat(k,1841) + lmat(k,1841)
         mat(k,1843) = mat(k,1843) + lmat(k,1843)
         mat(k,1861) = mat(k,1861) + lmat(k,1861)
         mat(k,1879) = mat(k,1879) + lmat(k,1879)
         mat(k,1886) = mat(k,1886) + lmat(k,1886)
         mat(k,1898) = mat(k,1898) + lmat(k,1898)
         mat(k,1906) = mat(k,1906) + lmat(k,1906)
         mat(k,1913) = mat(k,1913) + lmat(k,1913)
         mat(k,1927) = mat(k,1927) + lmat(k,1927)
         mat(k,1928) = mat(k,1928) + lmat(k,1928)
         mat(k,1940) = mat(k,1940) + lmat(k,1940)
         mat(k,1947) = mat(k,1947) + lmat(k,1947)
         mat(k,1962) = mat(k,1962) + lmat(k,1962)
         mat(k,1964) = mat(k,1964) + lmat(k,1964)
         mat(k,1969) = mat(k,1969) + lmat(k,1969)
         mat(k,1983) = mat(k,1983) + lmat(k,1983)
         mat(k,2016) = mat(k,2016) + lmat(k,2016)
         mat(k,2023) = mat(k,2023) + lmat(k,2023)
         mat(k,2029) = mat(k,2029) + lmat(k,2029)
         mat(k,2055) = lmat(k,2055)
         mat(k,2060) = mat(k,2060) + lmat(k,2060)
         mat(k,2065) = mat(k,2065) + lmat(k,2065)
         mat(k,2074) = lmat(k,2074)
         mat(k,2083) = lmat(k,2083)
         mat(k,2090) = mat(k,2090) + lmat(k,2090)
         mat(k, 91) = 0._r8
         mat(k, 117) = 0._r8
         mat(k, 286) = 0._r8
         mat(k, 290) = 0._r8
         mat(k, 295) = 0._r8
         mat(k, 298) = 0._r8
         mat(k, 314) = 0._r8
         mat(k, 356) = 0._r8
         mat(k, 368) = 0._r8
         mat(k, 390) = 0._r8
         mat(k, 391) = 0._r8
         mat(k, 392) = 0._r8
         mat(k, 398) = 0._r8
         mat(k, 426) = 0._r8
         mat(k, 427) = 0._r8
         mat(k, 435) = 0._r8
         mat(k, 439) = 0._r8
         mat(k, 442) = 0._r8
         mat(k, 445) = 0._r8
         mat(k, 447) = 0._r8
         mat(k, 448) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 464) = 0._r8
         mat(k, 467) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 470) = 0._r8
         mat(k, 479) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 588) = 0._r8
         mat(k, 589) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 645) = 0._r8
         mat(k, 646) = 0._r8
         mat(k, 648) = 0._r8
         mat(k, 650) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 726) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 779) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 793) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 803) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 834) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 839) = 0._r8
         mat(k, 849) = 0._r8
         mat(k, 855) = 0._r8
         mat(k, 856) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 881) = 0._r8
         mat(k, 882) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 911) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 987) = 0._r8
         mat(k, 988) = 0._r8
         mat(k, 995) = 0._r8
         mat(k, 996) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1035) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1056) = 0._r8
         mat(k,1059) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1098) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1210) = 0._r8
         mat(k,1212) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1233) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1334) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1444) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1447) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1481) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1526) = 0._r8
         mat(k,1527) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1531) = 0._r8
         mat(k,1532) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1540) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1573) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1601) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1619) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1649) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1664) = 0._r8
         mat(k,1682) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1704) = 0._r8
         mat(k,1705) = 0._r8
         mat(k,1706) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1716) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1722) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1729) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1731) = 0._r8
         mat(k,1732) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1739) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1744) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1748) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1764) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1775) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1777) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1789) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1797) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1802) = 0._r8
         mat(k,1803) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1813) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1832) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1842) = 0._r8
         mat(k,1845) = 0._r8
         mat(k,1846) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1849) = 0._r8
         mat(k,1851) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1855) = 0._r8
         mat(k,1856) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1858) = 0._r8
         mat(k,1862) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1871) = 0._r8
         mat(k,1872) = 0._r8
         mat(k,1875) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1882) = 0._r8
         mat(k,1883) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1890) = 0._r8
         mat(k,1891) = 0._r8
         mat(k,1892) = 0._r8
         mat(k,1893) = 0._r8
         mat(k,1894) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1912) = 0._r8
         mat(k,1919) = 0._r8
         mat(k,1922) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1934) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1966) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1973) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,1998) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2007) = 0._r8
         mat(k,2009) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2013) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2024) = 0._r8
         mat(k,2027) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2049) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2068) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2082) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2089) = 0._r8
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
         mat(k, 22) = mat(k, 22) - dti(k)
         mat(k, 23) = mat(k, 23) - dti(k)
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 25) = mat(k, 25) - dti(k)
         mat(k, 26) = mat(k, 26) - dti(k)
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 40) = mat(k, 40) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 90) = mat(k, 90) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 161) = mat(k, 161) - dti(k)
         mat(k, 167) = mat(k, 167) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 186) = mat(k, 186) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 256) = mat(k, 256) - dti(k)
         mat(k, 264) = mat(k, 264) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 309) = mat(k, 309) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 341) = mat(k, 341) - dti(k)
         mat(k, 351) = mat(k, 351) - dti(k)
         mat(k, 359) = mat(k, 359) - dti(k)
         mat(k, 369) = mat(k, 369) - dti(k)
         mat(k, 377) = mat(k, 377) - dti(k)
         mat(k, 388) = mat(k, 388) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 416) = mat(k, 416) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 446) = mat(k, 446) - dti(k)
         mat(k, 461) = mat(k, 461) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 490) = mat(k, 490) - dti(k)
         mat(k, 504) = mat(k, 504) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 564) = mat(k, 564) - dti(k)
         mat(k, 582) = mat(k, 582) - dti(k)
         mat(k, 594) = mat(k, 594) - dti(k)
         mat(k, 609) = mat(k, 609) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 647) = mat(k, 647) - dti(k)
         mat(k, 669) = mat(k, 669) - dti(k)
         mat(k, 691) = mat(k, 691) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 787) = mat(k, 787) - dti(k)
         mat(k, 811) = mat(k, 811) - dti(k)
         mat(k, 833) = mat(k, 833) - dti(k)
         mat(k, 866) = mat(k, 866) - dti(k)
         mat(k, 900) = mat(k, 900) - dti(k)
         mat(k, 926) = mat(k, 926) - dti(k)
         mat(k, 984) = mat(k, 984) - dti(k)
         mat(k,1014) = mat(k,1014) - dti(k)
         mat(k,1058) = mat(k,1058) - dti(k)
         mat(k,1104) = mat(k,1104) - dti(k)
         mat(k,1153) = mat(k,1153) - dti(k)
         mat(k,1197) = mat(k,1197) - dti(k)
         mat(k,1241) = mat(k,1241) - dti(k)
         mat(k,1284) = mat(k,1284) - dti(k)
         mat(k,1327) = mat(k,1327) - dti(k)
         mat(k,1388) = mat(k,1388) - dti(k)
         mat(k,1438) = mat(k,1438) - dti(k)
         mat(k,1487) = mat(k,1487) - dti(k)
         mat(k,1529) = mat(k,1529) - dti(k)
         mat(k,1566) = mat(k,1566) - dti(k)
         mat(k,1612) = mat(k,1612) - dti(k)
         mat(k,1656) = mat(k,1656) - dti(k)
         mat(k,1700) = mat(k,1700) - dti(k)
         mat(k,1741) = mat(k,1741) - dti(k)
         mat(k,1778) = mat(k,1778) - dti(k)
         mat(k,1827) = mat(k,1827) - dti(k)
         mat(k,1861) = mat(k,1861) - dti(k)
         mat(k,1898) = mat(k,1898) - dti(k)
         mat(k,1940) = mat(k,1940) - dti(k)
         mat(k,1983) = mat(k,1983) - dti(k)
         mat(k,2029) = mat(k,2029) - dti(k)
         mat(k,2090) = mat(k,2090) - dti(k)
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
      call nlnmat09( avec_len, mat, y, rxt )
      call nlnmat10( avec_len, mat, y, rxt )
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
