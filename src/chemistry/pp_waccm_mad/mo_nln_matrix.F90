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
         mat(k,453) = rxt(k,487)*y(k,24)
         mat(k,1500) = rxt(k,487)*y(k,2)
         mat(k,1453) = (rxt(k,549)+rxt(k,554))*y(k,45)
         mat(k,150) = (rxt(k,549)+rxt(k,554))*y(k,40)
         mat(k,457) = -(4._r8*rxt(k,484)*y(k,2) + (rxt(k,485) + rxt(k,486) + rxt(k,487) &
                      ) * y(k,24) + rxt(k,488)*y(k,43) + rxt(k,489)*y(k,51) + rxt(k,490) &
                      *y(k,52) + rxt(k,492)*y(k,54) + rxt(k,493)*y(k,104))
         mat(k,1506) = -(rxt(k,485) + rxt(k,486) + rxt(k,487)) * y(k,2)
         mat(k,683) = -rxt(k,488)*y(k,2)
         mat(k,1738) = -rxt(k,489)*y(k,2)
         mat(k,1134) = -rxt(k,490)*y(k,2)
         mat(k,966) = -rxt(k,492)*y(k,2)
         mat(k,751) = -rxt(k,493)*y(k,2)
         mat(k,74) = rxt(k,491)*y(k,54)
         mat(k,193) = rxt(k,501)*y(k,95)
         mat(k,153) = rxt(k,496)*y(k,54)
         mat(k,966) = mat(k,966) + rxt(k,491)*y(k,3) + rxt(k,496)*y(k,45)
         mat(k,1419) = rxt(k,483)*y(k,59)
         mat(k,329) = rxt(k,483)*y(k,56)
         mat(k,571) = rxt(k,501)*y(k,36)
         mat(k,71) = -(rxt(k,491)*y(k,54))
         mat(k,947) = -rxt(k,491)*y(k,3)
         mat(k,454) = rxt(k,490)*y(k,52)
         mat(k,1126) = rxt(k,490)*y(k,2)
         mat(k,541) = -(rxt(k,445)*y(k,60) + rxt(k,481)*y(k,59) + rxt(k,525)*y(k,53) &
                      + rxt(k,526)*y(k,54) + rxt(k,527)*y(k,104))
         mat(k,1662) = -rxt(k,445)*y(k,14)
         mat(k,330) = -rxt(k,481)*y(k,14)
         mat(k,1218) = -rxt(k,525)*y(k,14)
         mat(k,967) = -rxt(k,526)*y(k,14)
         mat(k,752) = -rxt(k,527)*y(k,14)
         mat(k,258) = rxt(k,452)*y(k,24) + rxt(k,529)*y(k,51)
         mat(k,41) = .300_r8*rxt(k,530)*y(k,104)
         mat(k,1507) = rxt(k,452)*y(k,18)
         mat(k,1741) = rxt(k,529)*y(k,18)
         mat(k,752) = mat(k,752) + .300_r8*rxt(k,530)*y(k,19)
         mat(k,257) = -(rxt(k,452)*y(k,24) + rxt(k,528)*y(k,43) + rxt(k,529)*y(k,51))
         mat(k,1504) = -rxt(k,452)*y(k,18)
         mat(k,680) = -rxt(k,528)*y(k,18)
         mat(k,1733) = -rxt(k,529)*y(k,18)
         mat(k,40) = .700_r8*rxt(k,530)*y(k,104)
         mat(k,747) = .700_r8*rxt(k,530)*y(k,19)
         mat(k,39) = -(rxt(k,530)*y(k,104))
         mat(k,737) = -rxt(k,530)*y(k,19)
         mat(k,256) = rxt(k,528)*y(k,43)
         mat(k,673) = rxt(k,528)*y(k,18)
         mat(k,1499) = 2.000_r8*rxt(k,454)*y(k,24)
         mat(k,199) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,40) + rxt(k,458)*y(k,60)
         mat(k,1452) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,25) + (rxt(k,542) &
                       +rxt(k,548)+rxt(k,553))*y(k,46)
         mat(k,183) = (rxt(k,542)+rxt(k,548)+rxt(k,553))*y(k,40)
         mat(k,1652) = rxt(k,458)*y(k,25)
         mat(k,1497) = 2.000_r8*rxt(k,479)*y(k,24)
         mat(k,1532) = -(rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68) + rxt(k,281)*y(k,81) &
                      + rxt(k,310)*y(k,98) + rxt(k,337)*y(k,105) + rxt(k,346)*y(k,106) &
                      + rxt(k,452)*y(k,18) + (4._r8*rxt(k,453) + 4._r8*rxt(k,454) &
                      + 4._r8*rxt(k,455) + 4._r8*rxt(k,479)) * y(k,24) + rxt(k,456) &
                      *y(k,43) + rxt(k,457)*y(k,51) + rxt(k,459)*y(k,52) + rxt(k,462) &
                      *y(k,54) + (rxt(k,463) + rxt(k,464)) * y(k,104) + (rxt(k,485) &
                      + rxt(k,486) + rxt(k,487)) * y(k,2))
         mat(k,802) = -rxt(k,111)*y(k,24)
         mat(k,640) = -rxt(k,123)*y(k,24)
         mat(k,728) = -rxt(k,281)*y(k,24)
         mat(k,1319) = -rxt(k,310)*y(k,24)
         mat(k,1608) = -rxt(k,337)*y(k,24)
         mat(k,1644) = -rxt(k,346)*y(k,24)
         mat(k,263) = -rxt(k,452)*y(k,24)
         mat(k,701) = -rxt(k,456)*y(k,24)
         mat(k,1765) = -rxt(k,457)*y(k,24)
         mat(k,1162) = -rxt(k,459)*y(k,24)
         mat(k,992) = -rxt(k,462)*y(k,24)
         mat(k,770) = -(rxt(k,463) + rxt(k,464)) * y(k,24)
         mat(k,470) = -(rxt(k,485) + rxt(k,486) + rxt(k,487)) * y(k,24)
         mat(k,208) = rxt(k,460)*y(k,54)
         mat(k,1489) = rxt(k,478)*y(k,95)
         mat(k,701) = mat(k,701) + rxt(k,450)*y(k,60)
         mat(k,188) = rxt(k,468)*y(k,54) + rxt(k,467)*y(k,60) + rxt(k,469)*y(k,104)
         mat(k,992) = mat(k,992) + rxt(k,460)*y(k,25) + rxt(k,468)*y(k,46)
         mat(k,1444) = rxt(k,451)*y(k,60)
         mat(k,1687) = rxt(k,450)*y(k,43) + rxt(k,467)*y(k,46) + rxt(k,451)*y(k,56)
         mat(k,588) = rxt(k,478)*y(k,40)
         mat(k,770) = mat(k,770) + rxt(k,469)*y(k,46)
         mat(k,201) = -(rxt(k,458)*y(k,60) + rxt(k,460)*y(k,54) + rxt(k,461)*y(k,104) &
                      + (rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,40))
         mat(k,1656) = -rxt(k,458)*y(k,25)
         mat(k,957) = -rxt(k,460)*y(k,25)
         mat(k,746) = -rxt(k,461)*y(k,25)
         mat(k,1457) = -(rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,25)
         mat(k,1502) = rxt(k,459)*y(k,52)
         mat(k,1128) = rxt(k,459)*y(k,24)
         mat(k,128) = -((rxt(k,532) + rxt(k,536)) * y(k,104))
         mat(k,741) = -(rxt(k,532) + rxt(k,536)) * y(k,27)
         mat(k,536) = rxt(k,525)*y(k,53) + rxt(k,526)*y(k,54) + rxt(k,481)*y(k,59) &
                      + rxt(k,445)*y(k,60) + rxt(k,527)*y(k,104)
         mat(k,1376) = rxt(k,573)*y(k,107)
         mat(k,1213) = rxt(k,525)*y(k,14)
         mat(k,951) = rxt(k,526)*y(k,14)
         mat(k,326) = rxt(k,481)*y(k,14)
         mat(k,1654) = rxt(k,445)*y(k,14)
         mat(k,741) = mat(k,741) + rxt(k,527)*y(k,14)
         mat(k,241) = rxt(k,573)*y(k,28)
         mat(k,4) = -(rxt(k,506)*y(k,95))
         mat(k,561) = -rxt(k,506)*y(k,29)
         mat(k,12) = -(rxt(k,507)*y(k,95))
         mat(k,562) = -rxt(k,507)*y(k,30)
         mat(k,1403) = -(rxt(k,307)*y(k,93) + rxt(k,311)*y(k,98) + rxt(k,325)*y(k,101) &
                      + rxt(k,330)*y(k,102) + rxt(k,338)*y(k,105) + rxt(k,347) &
                      *y(k,106) + rxt(k,363)*y(k,88) + rxt(k,573)*y(k,107))
         mat(k,113) = -rxt(k,307)*y(k,28)
         mat(k,1316) = -rxt(k,311)*y(k,28)
         mat(k,410) = -rxt(k,325)*y(k,28)
         mat(k,101) = -rxt(k,330)*y(k,28)
         mat(k,1605) = -rxt(k,338)*y(k,28)
         mat(k,1641) = -rxt(k,347)*y(k,28)
         mat(k,890) = -rxt(k,363)*y(k,28)
         mat(k,253) = -rxt(k,573)*y(k,28)
         mat(k,1529) = rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68)
         mat(k,131) = (rxt(k,532)+rxt(k,536))*y(k,104)
         mat(k,1720) = rxt(k,112)*y(k,65)
         mat(k,1486) = rxt(k,125)*y(k,68)
         mat(k,1114) = rxt(k,119)*y(k,65)
         mat(k,1762) = rxt(k,275)*y(k,65) + (rxt(k,117)+rxt(k,118))*y(k,67)
         mat(k,1159) = rxt(k,276)*y(k,65) + (rxt(k,115)+rxt(k,116))*y(k,67)
         mat(k,989) = rxt(k,120)*y(k,65)
         mat(k,1280) = rxt(k,121)*y(k,65)
         mat(k,1441) = rxt(k,127)*y(k,68)
         mat(k,1684) = (rxt(k,109)+rxt(k,110))*y(k,65) + rxt(k,122)*y(k,68)
         mat(k,799) = rxt(k,111)*y(k,24) + rxt(k,112)*y(k,32) + rxt(k,119)*y(k,42) &
                      + rxt(k,275)*y(k,51) + rxt(k,276)*y(k,52) + rxt(k,120)*y(k,54) &
                      + rxt(k,121)*y(k,55) + (rxt(k,109)+rxt(k,110))*y(k,60) &
                      + rxt(k,181)*y(k,73) + (rxt(k,165)+rxt(k,253))*y(k,75) + ( &
                      + rxt(k,163)+rxt(k,260))*y(k,77) + rxt(k,234)*y(k,88) &
                      + rxt(k,216)*y(k,89) + rxt(k,199)*y(k,92) + rxt(k,251)*y(k,99)
         mat(k,321) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77) + rxt(k,241)*y(k,88) &
                      + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92) + rxt(k,148)*y(k,99)
         mat(k,511) = (rxt(k,117)+rxt(k,118))*y(k,51) + (rxt(k,115)+rxt(k,116)) &
                      *y(k,52) + rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + rxt(k,226)*y(k,89) + rxt(k,208)*y(k,92) + rxt(k,150)*y(k,99)
         mat(k,637) = rxt(k,123)*y(k,24) + rxt(k,125)*y(k,40) + rxt(k,127)*y(k,56) &
                      + rxt(k,122)*y(k,60) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77) + rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92) + rxt(k,146)*y(k,99)
         mat(k,1032) = rxt(k,301)*y(k,91)
         mat(k,306) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77) &
                      + rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92) &
                      + rxt(k,144)*y(k,99)
         mat(k,842) = rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67) &
                      + rxt(k,186)*y(k,68) + rxt(k,184)*y(k,71)
         mat(k,1202) = (rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67) + rxt(k,220)*y(k,68) &
                      + rxt(k,198)*y(k,71)
         mat(k,1572) = (rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67) + rxt(k,169)*y(k,68) &
                      + rxt(k,167)*y(k,71)
         mat(k,890) = mat(k,890) + rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) &
                      + rxt(k,244)*y(k,67) + rxt(k,239)*y(k,68) + rxt(k,237)*y(k,71)
         mat(k,933) = rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67) &
                      + rxt(k,222)*y(k,68) + rxt(k,219)*y(k,71)
         mat(k,120) = rxt(k,301)*y(k,69) + rxt(k,302)*y(k,110)
         mat(k,1074) = rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) + rxt(k,208)*y(k,67) &
                      + rxt(k,204)*y(k,68) + rxt(k,202)*y(k,71)
         mat(k,1361) = rxt(k,251)*y(k,65) + rxt(k,148)*y(k,66) + rxt(k,150)*y(k,67) &
                      + rxt(k,146)*y(k,68) + rxt(k,144)*y(k,71)
         mat(k,767) = (rxt(k,532)+rxt(k,536))*y(k,27)
         mat(k,1815) = rxt(k,302)*y(k,91)
         mat(k,168) = -(rxt(k,503)*y(k,33) + rxt(k,504)*y(k,110) + rxt(k,505)*y(k,42))
         mat(k,417) = -rxt(k,503)*y(k,31)
         mat(k,1784) = -rxt(k,504)*y(k,31)
         mat(k,1087) = -rxt(k,505)*y(k,31)
         mat(k,5) = 2.000_r8*rxt(k,506)*y(k,95)
         mat(k,13) = rxt(k,507)*y(k,95)
         mat(k,564) = 2.000_r8*rxt(k,506)*y(k,29) + rxt(k,507)*y(k,30)
         mat(k,1728) = -(rxt(k,100)*y(k,61) + rxt(k,112)*y(k,65) + rxt(k,124)*y(k,68) &
                      + rxt(k,282)*y(k,81) + rxt(k,304)*y(k,92) + rxt(k,312)*y(k,98) &
                      + rxt(k,326)*y(k,101) + rxt(k,339)*y(k,105) + (rxt(k,403) &
                      + rxt(k,404) + rxt(k,405)) * y(k,43) + rxt(k,406)*y(k,55) &
                      + rxt(k,409)*y(k,56))
         mat(k,613) = -rxt(k,100)*y(k,32)
         mat(k,807) = -rxt(k,112)*y(k,32)
         mat(k,645) = -rxt(k,124)*y(k,32)
         mat(k,732) = -rxt(k,282)*y(k,32)
         mat(k,1082) = -rxt(k,304)*y(k,32)
         mat(k,1324) = -rxt(k,312)*y(k,32)
         mat(k,414) = -rxt(k,326)*y(k,32)
         mat(k,1613) = -rxt(k,339)*y(k,32)
         mat(k,705) = -(rxt(k,403) + rxt(k,404) + rxt(k,405)) * y(k,32)
         mat(k,1288) = -rxt(k,406)*y(k,32)
         mat(k,1449) = -rxt(k,409)*y(k,32)
         mat(k,558) = rxt(k,527)*y(k,104)
         mat(k,132) = rxt(k,536)*y(k,104)
         mat(k,174) = rxt(k,503)*y(k,33)
         mat(k,435) = rxt(k,503)*y(k,31) + rxt(k,401)*y(k,54) + rxt(k,447)*y(k,60) &
                      + rxt(k,384)*y(k,95) + rxt(k,410)*y(k,104) + rxt(k,349)*y(k,106)
         mat(k,197) = rxt(k,501)*y(k,95)
         mat(k,1494) = rxt(k,478)*y(k,95)
         mat(k,281) = rxt(k,433)*y(k,104)
         mat(k,997) = rxt(k,401)*y(k,33) + rxt(k,413)*y(k,104)
         mat(k,1692) = rxt(k,447)*y(k,33)
         mat(k,613) = mat(k,613) + rxt(k,190)*y(k,73) + rxt(k,142)*y(k,75) &
                      + rxt(k,172)*y(k,77)
         mat(k,368) = rxt(k,194)*y(k,73) + rxt(k,159)*y(k,75) + rxt(k,177)*y(k,77)
         mat(k,352) = rxt(k,182)*y(k,73) + rxt(k,176)*y(k,75) + rxt(k,164)*y(k,77)
         mat(k,807) = mat(k,807) + rxt(k,181)*y(k,73) + (rxt(k,165)+rxt(k,253)) &
                      *y(k,75) + (rxt(k,163)+rxt(k,260))*y(k,77)
         mat(k,323) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77)
         mat(k,513) = rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77)
         mat(k,645) = mat(k,645) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77)
         mat(k,1040) = rxt(k,133)*y(k,70) + rxt(k,377)*y(k,72) + rxt(k,378)*y(k,73) &
                      + rxt(k,136)*y(k,75) + rxt(k,139)*y(k,77) + rxt(k,376)*y(k,78)
         mat(k,37) = rxt(k,133)*y(k,69)
         mat(k,308) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77)
         mat(k,106) = rxt(k,377)*y(k,69)
         mat(k,850) = rxt(k,190)*y(k,61) + rxt(k,194)*y(k,62) + rxt(k,182)*y(k,63) &
                      + rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67) &
                      + rxt(k,186)*y(k,68) + rxt(k,378)*y(k,69) + rxt(k,184)*y(k,71) &
                      + rxt(k,196)*y(k,81) + rxt(k,192)*y(k,82) + rxt(k,195)*y(k,84) &
                      + rxt(k,188)*y(k,85) + rxt(k,193)*y(k,86) + rxt(k,185)*y(k,98)
         mat(k,1210) = rxt(k,142)*y(k,61) + rxt(k,159)*y(k,62) + rxt(k,176)*y(k,63) + ( &
                      + rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67) + rxt(k,220)*y(k,68) &
                      + rxt(k,136)*y(k,69) + rxt(k,198)*y(k,71) + rxt(k,161)*y(k,81) &
                      + rxt(k,157)*y(k,82) + rxt(k,160)*y(k,84) + (rxt(k,231) &
                       +rxt(k,257))*y(k,85) + rxt(k,158)*y(k,86) + rxt(k,209)*y(k,98)
         mat(k,1580) = rxt(k,172)*y(k,61) + rxt(k,177)*y(k,62) + rxt(k,164)*y(k,63) + ( &
                      + rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67) + rxt(k,169)*y(k,68) &
                      + rxt(k,139)*y(k,69) + rxt(k,167)*y(k,71) + rxt(k,179)*y(k,81) &
                      + rxt(k,174)*y(k,82) + rxt(k,178)*y(k,84) + (rxt(k,170) &
                       +rxt(k,258))*y(k,85) + rxt(k,175)*y(k,86) + rxt(k,168)*y(k,98)
         mat(k,137) = rxt(k,376)*y(k,69)
         mat(k,732) = mat(k,732) + rxt(k,196)*y(k,73) + rxt(k,161)*y(k,75) &
                      + rxt(k,179)*y(k,77)
         mat(k,382) = rxt(k,192)*y(k,73) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77)
         mat(k,492) = rxt(k,195)*y(k,73) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77)
         mat(k,533) = rxt(k,188)*y(k,73) + (rxt(k,231)+rxt(k,257))*y(k,75) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,77)
         mat(k,398) = rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77)
         mat(k,592) = rxt(k,384)*y(k,33) + rxt(k,501)*y(k,36) + rxt(k,478)*y(k,40)
         mat(k,1324) = mat(k,1324) + rxt(k,185)*y(k,73) + rxt(k,209)*y(k,75) &
                      + rxt(k,168)*y(k,77)
         mat(k,774) = rxt(k,527)*y(k,14) + rxt(k,536)*y(k,27) + rxt(k,410)*y(k,33) &
                      + rxt(k,433)*y(k,48) + rxt(k,413)*y(k,54)
         mat(k,1649) = rxt(k,349)*y(k,33)
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
         mat(k,419) = -((rxt(k,348) + rxt(k,349)) * y(k,106) + rxt(k,384)*y(k,95) &
                      + rxt(k,401)*y(k,54) + rxt(k,410)*y(k,104) + rxt(k,447)*y(k,60) &
                      + rxt(k,503)*y(k,31))
         mat(k,1619) = -(rxt(k,348) + rxt(k,349)) * y(k,33)
         mat(k,570) = -rxt(k,384)*y(k,33)
         mat(k,965) = -rxt(k,401)*y(k,33)
         mat(k,750) = -rxt(k,410)*y(k,33)
         mat(k,1660) = -rxt(k,447)*y(k,33)
         mat(k,170) = -rxt(k,503)*y(k,33)
         mat(k,1697) = rxt(k,403)*y(k,43)
         mat(k,682) = rxt(k,403)*y(k,32)
         mat(k,79) = -(rxt(k,402)*y(k,54) + rxt(k,411)*y(k,104) + rxt(k,448)*y(k,60))
         mat(k,948) = -rxt(k,402)*y(k,35)
         mat(k,739) = -rxt(k,411)*y(k,35)
         mat(k,1653) = -rxt(k,448)*y(k,35)
         mat(k,675) = 2.000_r8*rxt(k,417)*y(k,43)
         mat(k,739) = mat(k,739) + 2.000_r8*rxt(k,416)*y(k,104)
         mat(k,191) = -(rxt(k,494)*y(k,54) + rxt(k,495)*y(k,104) + (rxt(k,500) &
                      + rxt(k,501)) * y(k,95))
         mat(k,956) = -rxt(k,494)*y(k,36)
         mat(k,745) = -rxt(k,495)*y(k,36)
         mat(k,565) = -(rxt(k,500) + rxt(k,501)) * y(k,36)
         mat(k,537) = rxt(k,481)*y(k,59)
         mat(k,679) = rxt(k,482)*y(k,59)
         mat(k,327) = rxt(k,481)*y(k,14) + rxt(k,482)*y(k,43)
         mat(k,1488) = -(rxt(k,101)*y(k,62) + rxt(k,103)*y(k,61) + rxt(k,125)*y(k,68) &
                      + (rxt(k,271) + rxt(k,293)) * y(k,83) + rxt(k,284)*y(k,81) &
                      + rxt(k,313)*y(k,98) + rxt(k,340)*y(k,105) + rxt(k,351)*y(k,106) &
                      + rxt(k,465)*y(k,54) + rxt(k,466)*y(k,104) + (rxt(k,477) &
                      + rxt(k,478)) * y(k,95) + (rxt(k,542) + rxt(k,548) + rxt(k,553) &
                      ) * y(k,46) + (rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,25) &
                      + (rxt(k,549) + rxt(k,554)) * y(k,45))
         mat(k,365) = -rxt(k,101)*y(k,40)
         mat(k,610) = -rxt(k,103)*y(k,40)
         mat(k,639) = -rxt(k,125)*y(k,40)
         mat(k,667) = -(rxt(k,271) + rxt(k,293)) * y(k,40)
         mat(k,727) = -rxt(k,284)*y(k,40)
         mat(k,1318) = -rxt(k,313)*y(k,40)
         mat(k,1607) = -rxt(k,340)*y(k,40)
         mat(k,1643) = -rxt(k,351)*y(k,40)
         mat(k,991) = -rxt(k,465)*y(k,40)
         mat(k,769) = -rxt(k,466)*y(k,40)
         mat(k,587) = -(rxt(k,477) + rxt(k,478)) * y(k,40)
         mat(k,187) = -(rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,40)
         mat(k,207) = -(rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,40)
         mat(k,156) = -(rxt(k,549) + rxt(k,554)) * y(k,40)
         mat(k,553) = rxt(k,445)*y(k,60)
         mat(k,1531) = rxt(k,464)*y(k,104)
         mat(k,1722) = rxt(k,100)*y(k,61)
         mat(k,430) = rxt(k,447)*y(k,60)
         mat(k,83) = rxt(k,448)*y(k,60)
         mat(k,1116) = rxt(k,104)*y(k,61) + rxt(k,294)*y(k,86)
         mat(k,700) = rxt(k,449)*y(k,60)
         mat(k,187) = mat(k,187) + rxt(k,467)*y(k,60)
         mat(k,1686) = rxt(k,445)*y(k,14) + rxt(k,447)*y(k,33) + rxt(k,448)*y(k,35) &
                      + rxt(k,449)*y(k,43) + rxt(k,467)*y(k,46)
         mat(k,610) = mat(k,610) + rxt(k,100)*y(k,32) + rxt(k,104)*y(k,42)
         mat(k,349) = rxt(k,182)*y(k,73) + (rxt(k,176)+2.000_r8*rxt(k,262))*y(k,75) + ( &
                      + rxt(k,164)+2.000_r8*rxt(k,263))*y(k,77) + rxt(k,235)*y(k,88) &
                      + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92) + rxt(k,252)*y(k,99)
         mat(k,844) = rxt(k,182)*y(k,63) + rxt(k,193)*y(k,86)
         mat(k,1204) = (rxt(k,176)+2.000_r8*rxt(k,262))*y(k,63) + rxt(k,158)*y(k,86)
         mat(k,1574) = (rxt(k,164)+2.000_r8*rxt(k,263))*y(k,63) + rxt(k,175)*y(k,86)
         mat(k,396) = rxt(k,294)*y(k,42) + rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) &
                      + rxt(k,175)*y(k,77) + rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) &
                      + rxt(k,211)*y(k,92) + rxt(k,152)*y(k,99)
         mat(k,892) = rxt(k,235)*y(k,63) + rxt(k,246)*y(k,86)
         mat(k,935) = rxt(k,217)*y(k,63) + rxt(k,228)*y(k,86)
         mat(k,1076) = rxt(k,200)*y(k,63) + rxt(k,211)*y(k,86)
         mat(k,1363) = rxt(k,252)*y(k,63) + rxt(k,152)*y(k,86)
         mat(k,769) = mat(k,769) + rxt(k,464)*y(k,24)
         mat(k,167) = rxt(k,503)*y(k,33) + rxt(k,505)*y(k,42) + rxt(k,504)*y(k,110)
         mat(k,416) = rxt(k,503)*y(k,31)
         mat(k,1085) = rxt(k,505)*y(k,31)
         mat(k,1773) = rxt(k,504)*y(k,31)
         mat(k,1107) = -(rxt(k,104)*y(k,61) + rxt(k,119)*y(k,65) + rxt(k,285)*y(k,81) &
                      + rxt(k,290)*y(k,85) + rxt(k,294)*y(k,86) + rxt(k,295)*y(k,83) &
                      + rxt(k,314)*y(k,98) + rxt(k,352)*y(k,106) + rxt(k,442)*y(k,104) &
                      + rxt(k,505)*y(k,31))
         mat(k,605) = -rxt(k,104)*y(k,42)
         mat(k,792) = -rxt(k,119)*y(k,42)
         mat(k,720) = -rxt(k,285)*y(k,42)
         mat(k,526) = -rxt(k,290)*y(k,42)
         mat(k,391) = -rxt(k,294)*y(k,42)
         mat(k,660) = -rxt(k,295)*y(k,42)
         mat(k,1309) = -rxt(k,314)*y(k,42)
         mat(k,1634) = -rxt(k,352)*y(k,42)
         mat(k,761) = -rxt(k,442)*y(k,42)
         mat(k,172) = -rxt(k,505)*y(k,42)
         mat(k,547) = rxt(k,525)*y(k,53)
         mat(k,204) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,40)
         mat(k,1479) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,25) + rxt(k,293)*y(k,83)
         mat(k,234) = rxt(k,137)*y(k,75) + rxt(k,140)*y(k,77) + rxt(k,288)*y(k,84) &
                      + rxt(k,292)*y(k,85)
         mat(k,1152) = rxt(k,441)*y(k,104)
         mat(k,1230) = rxt(k,525)*y(k,14)
         mat(k,835) = rxt(k,183)*y(k,83) + 2.000_r8*rxt(k,180)*y(k,87)
         mat(k,23) = rxt(k,135)*y(k,110)
         mat(k,1195) = rxt(k,137)*y(k,50) + (rxt(k,187)+rxt(k,259))*y(k,83) + ( &
                      + 2.000_r8*rxt(k,141)+2.000_r8*rxt(k,264))*y(k,87)
         mat(k,27) = rxt(k,138)*y(k,110)
         mat(k,1565) = rxt(k,140)*y(k,50) + (rxt(k,166)+rxt(k,261))*y(k,83) + ( &
                      + 2.000_r8*rxt(k,162)+2.000_r8*rxt(k,265))*y(k,87)
         mat(k,660) = mat(k,660) + rxt(k,293)*y(k,40) + rxt(k,183)*y(k,73) + ( &
                      + rxt(k,187)+rxt(k,259))*y(k,75) + (rxt(k,166)+rxt(k,261)) &
                      *y(k,77)
         mat(k,485) = rxt(k,288)*y(k,50)
         mat(k,526) = mat(k,526) + rxt(k,292)*y(k,50)
         mat(k,444) = 2.000_r8*rxt(k,180)*y(k,73) + (2.000_r8*rxt(k,141) &
                       +2.000_r8*rxt(k,264))*y(k,75) + (2.000_r8*rxt(k,162) &
                       +2.000_r8*rxt(k,265))*y(k,77) + rxt(k,233)*y(k,88) + rxt(k,215) &
                      *y(k,89) + rxt(k,197)*y(k,92) + rxt(k,250)*y(k,99)
         mat(k,883) = rxt(k,233)*y(k,87)
         mat(k,926) = rxt(k,215)*y(k,87)
         mat(k,1067) = rxt(k,197)*y(k,87)
         mat(k,1354) = rxt(k,250)*y(k,87)
         mat(k,761) = mat(k,761) + rxt(k,441)*y(k,52)
         mat(k,1808) = rxt(k,135)*y(k,74) + rxt(k,138)*y(k,76)
         mat(k,686) = -(rxt(k,305)*y(k,92) + (rxt(k,403) + rxt(k,404) + rxt(k,405) &
                      ) * y(k,32) + rxt(k,407)*y(k,54) + rxt(k,408)*y(k,56) + rxt(k,412) &
                      *y(k,104) + 4._r8*rxt(k,417)*y(k,43) + rxt(k,429)*y(k,53) &
                      + rxt(k,434)*y(k,51) + rxt(k,439)*y(k,52) + (rxt(k,449) &
                      + rxt(k,450)) * y(k,60) + rxt(k,456)*y(k,24) + rxt(k,482) &
                      *y(k,59) + rxt(k,488)*y(k,2) + rxt(k,528)*y(k,18))
         mat(k,1057) = -rxt(k,305)*y(k,43)
         mat(k,1703) = -(rxt(k,403) + rxt(k,404) + rxt(k,405)) * y(k,43)
         mat(k,972) = -rxt(k,407)*y(k,43)
         mat(k,1424) = -rxt(k,408)*y(k,43)
         mat(k,754) = -rxt(k,412)*y(k,43)
         mat(k,1221) = -rxt(k,429)*y(k,43)
         mat(k,1745) = -rxt(k,434)*y(k,43)
         mat(k,1142) = -rxt(k,439)*y(k,43)
         mat(k,1667) = -(rxt(k,449) + rxt(k,450)) * y(k,43)
         mat(k,1512) = -rxt(k,456)*y(k,43)
         mat(k,332) = -rxt(k,482)*y(k,43)
         mat(k,460) = -rxt(k,488)*y(k,43)
         mat(k,259) = -rxt(k,528)*y(k,43)
         mat(k,460) = mat(k,460) + rxt(k,493)*y(k,104)
         mat(k,543) = rxt(k,525)*y(k,53) + rxt(k,526)*y(k,54) + rxt(k,481)*y(k,59) &
                      + rxt(k,445)*y(k,60)
         mat(k,259) = mat(k,259) + rxt(k,452)*y(k,24) + rxt(k,529)*y(k,51)
         mat(k,1512) = mat(k,1512) + rxt(k,452)*y(k,18) + rxt(k,463)*y(k,104)
         mat(k,129) = rxt(k,532)*y(k,104)
         mat(k,1703) = mat(k,1703) + rxt(k,406)*y(k,55) + rxt(k,312)*y(k,98)
         mat(k,80) = rxt(k,402)*y(k,54) + rxt(k,448)*y(k,60) + rxt(k,411)*y(k,104)
         mat(k,1469) = rxt(k,125)*y(k,68) + rxt(k,313)*y(k,98)
         mat(k,1097) = rxt(k,314)*y(k,98)
         mat(k,1745) = mat(k,1745) + rxt(k,529)*y(k,18)
         mat(k,1221) = mat(k,1221) + rxt(k,525)*y(k,14) + rxt(k,432)*y(k,104)
         mat(k,972) = mat(k,972) + rxt(k,526)*y(k,14) + rxt(k,402)*y(k,35) &
                      + rxt(k,342)*y(k,105)
         mat(k,1263) = rxt(k,406)*y(k,32)
         mat(k,1424) = mat(k,1424) + rxt(k,414)*y(k,104)
         mat(k,332) = mat(k,332) + rxt(k,481)*y(k,14)
         mat(k,1667) = mat(k,1667) + rxt(k,445)*y(k,14) + rxt(k,448)*y(k,35)
         mat(k,621) = rxt(k,125)*y(k,40)
         mat(k,1299) = rxt(k,312)*y(k,32) + rxt(k,313)*y(k,40) + rxt(k,314)*y(k,42)
         mat(k,754) = mat(k,754) + rxt(k,493)*y(k,2) + rxt(k,463)*y(k,24) + rxt(k,532) &
                      *y(k,27) + rxt(k,411)*y(k,35) + rxt(k,432)*y(k,53) + rxt(k,414) &
                      *y(k,56)
         mat(k,1588) = rxt(k,342)*y(k,54)
         mat(k,51) = -(rxt(k,418)*y(k,104))
         mat(k,738) = -rxt(k,418)*y(k,44)
         mat(k,674) = rxt(k,439)*y(k,52)
         mat(k,1125) = rxt(k,439)*y(k,43)
         mat(k,151) = -(rxt(k,496)*y(k,54) + (rxt(k,549) + rxt(k,554)) * y(k,40))
         mat(k,952) = -rxt(k,496)*y(k,45)
         mat(k,1455) = -(rxt(k,549) + rxt(k,554)) * y(k,45)
         mat(k,455) = rxt(k,488)*y(k,43)
         mat(k,677) = rxt(k,488)*y(k,2)
         mat(k,184) = -(rxt(k,467)*y(k,60) + rxt(k,468)*y(k,54) + rxt(k,469)*y(k,104) &
                      + (rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,40))
         mat(k,1655) = -rxt(k,467)*y(k,46)
         mat(k,955) = -rxt(k,468)*y(k,46)
         mat(k,744) = -rxt(k,469)*y(k,46)
         mat(k,1456) = -(rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,46)
         mat(k,1501) = rxt(k,456)*y(k,43)
         mat(k,200) = rxt(k,461)*y(k,104)
         mat(k,678) = rxt(k,456)*y(k,24)
         mat(k,744) = mat(k,744) + rxt(k,461)*y(k,25)
         mat(k,123) = -(rxt(k,335)*y(k,104))
         mat(k,740) = -rxt(k,335)*y(k,47)
         mat(k,1454) = rxt(k,284)*y(k,81)
         mat(k,1086) = rxt(k,285)*y(k,81)
         mat(k,1731) = rxt(k,344)*y(k,104)
         mat(k,708) = rxt(k,284)*y(k,40) + rxt(k,285)*y(k,42)
         mat(k,46) = rxt(k,300)*y(k,110)
         mat(k,740) = mat(k,740) + rxt(k,344)*y(k,51)
         mat(k,1781) = rxt(k,300)*y(k,90)
         mat(k,270) = -(rxt(k,421)*y(k,51) + (rxt(k,422) + rxt(k,423) + rxt(k,424) &
                      ) * y(k,52) + rxt(k,425)*y(k,55) + rxt(k,433)*y(k,104) + rxt(k,570) &
                      *y(k,99))
         mat(k,1734) = -rxt(k,421)*y(k,48)
         mat(k,1130) = -(rxt(k,422) + rxt(k,423) + rxt(k,424)) * y(k,48)
         mat(k,1257) = -rxt(k,425)*y(k,48)
         mat(k,748) = -rxt(k,433)*y(k,48)
         mat(k,1328) = -rxt(k,570)*y(k,48)
         mat(k,961) = rxt(k,419)*y(k,79) + rxt(k,567)*y(k,94)
         mat(k,1257) = mat(k,1257) + rxt(k,568)*y(k,94)
         mat(k,1014) = 1.100_r8*rxt(k,563)*y(k,80) + .200_r8*rxt(k,561)*y(k,88)
         mat(k,162) = rxt(k,419)*y(k,54)
         mat(k,89) = 1.100_r8*rxt(k,563)*y(k,69)
         mat(k,858) = .200_r8*rxt(k,561)*y(k,69)
         mat(k,178) = rxt(k,567)*y(k,54) + rxt(k,568)*y(k,55)
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
         mat(k,229) = -(rxt(k,137)*y(k,75) + rxt(k,140)*y(k,77) + rxt(k,288)*y(k,84) &
                      + rxt(k,292)*y(k,85))
         mat(k,1171) = -rxt(k,137)*y(k,50)
         mat(k,1541) = -rxt(k,140)*y(k,50)
         mat(k,475) = -rxt(k,288)*y(k,50)
         mat(k,516) = -rxt(k,292)*y(k,50)
         mat(k,1129) = rxt(k,440)*y(k,53)
         mat(k,1214) = rxt(k,440)*y(k,52)
         mat(k,1771) = -((rxt(k,106) + rxt(k,107)) * y(k,64) + (rxt(k,117) + rxt(k,118) &
                      ) * y(k,67) + rxt(k,131)*y(k,106) + (rxt(k,267) + rxt(k,274) &
                      ) * y(k,101) + rxt(k,275)*y(k,65) + rxt(k,344)*y(k,104) &
                      + rxt(k,421)*y(k,48) + rxt(k,430)*y(k,53) + rxt(k,434)*y(k,43) &
                      + rxt(k,435)*y(k,56) + rxt(k,436)*y(k,54) + rxt(k,457)*y(k,24) &
                      + rxt(k,489)*y(k,2) + rxt(k,529)*y(k,18) + rxt(k,572)*y(k,99))
         mat(k,218) = -(rxt(k,106) + rxt(k,107)) * y(k,51)
         mat(k,514) = -(rxt(k,117) + rxt(k,118)) * y(k,51)
         mat(k,1650) = -rxt(k,131)*y(k,51)
         mat(k,415) = -(rxt(k,267) + rxt(k,274)) * y(k,51)
         mat(k,808) = -rxt(k,275)*y(k,51)
         mat(k,775) = -rxt(k,344)*y(k,51)
         mat(k,282) = -rxt(k,421)*y(k,51)
         mat(k,1246) = -rxt(k,430)*y(k,51)
         mat(k,706) = -rxt(k,434)*y(k,51)
         mat(k,1450) = -rxt(k,435)*y(k,51)
         mat(k,998) = -rxt(k,436)*y(k,51)
         mat(k,1538) = -rxt(k,457)*y(k,51)
         mat(k,473) = -rxt(k,489)*y(k,51)
         mat(k,266) = -rxt(k,529)*y(k,51)
         mat(k,1370) = -rxt(k,572)*y(k,51)
         mat(k,1729) = rxt(k,282)*y(k,81) + rxt(k,304)*y(k,92)
         mat(k,282) = mat(k,282) + 2.000_r8*rxt(k,423)*y(k,52) + rxt(k,425)*y(k,55) &
                      + rxt(k,433)*y(k,104)
         mat(k,1168) = 2.000_r8*rxt(k,423)*y(k,48) + rxt(k,426)*y(k,54) + rxt(k,286) &
                      *y(k,81)
         mat(k,998) = mat(k,998) + rxt(k,426)*y(k,52)
         mat(k,1289) = rxt(k,425)*y(k,48) + rxt(k,420)*y(k,79)
         mat(k,614) = rxt(k,243)*y(k,88) + rxt(k,225)*y(k,89) + rxt(k,207)*y(k,92)
         mat(k,369) = rxt(k,247)*y(k,88) + rxt(k,229)*y(k,89) + rxt(k,212)*y(k,92)
         mat(k,353) = rxt(k,235)*y(k,88) + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92)
         mat(k,808) = mat(k,808) + rxt(k,234)*y(k,88) + rxt(k,216)*y(k,89) &
                      + rxt(k,199)*y(k,92)
         mat(k,324) = rxt(k,241)*y(k,88) + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92)
         mat(k,514) = mat(k,514) + rxt(k,244)*y(k,88) + rxt(k,226)*y(k,89) &
                      + rxt(k,208)*y(k,92)
         mat(k,646) = rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) + rxt(k,204)*y(k,92)
         mat(k,1041) = rxt(k,298)*y(k,89) + rxt(k,299)*y(k,90) + rxt(k,301)*y(k,91) &
                      + rxt(k,303)*y(k,92) + rxt(k,379)*y(k,93)
         mat(k,309) = rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92)
         mat(k,166) = rxt(k,420)*y(k,55)
         mat(k,733) = rxt(k,282)*y(k,32) + rxt(k,286)*y(k,52) + rxt(k,249)*y(k,88) &
                      + rxt(k,232)*y(k,89) + rxt(k,214)*y(k,92)
         mat(k,383) = rxt(k,245)*y(k,88) + rxt(k,227)*y(k,89) + rxt(k,210)*y(k,92)
         mat(k,671) = rxt(k,236)*y(k,88) + rxt(k,218)*y(k,89) + rxt(k,201)*y(k,92)
         mat(k,493) = rxt(k,248)*y(k,88) + rxt(k,230)*y(k,89) + rxt(k,213)*y(k,92)
         mat(k,534) = rxt(k,240)*y(k,88) + rxt(k,223)*y(k,89) + rxt(k,205)*y(k,92)
         mat(k,399) = rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) + rxt(k,211)*y(k,92)
         mat(k,450) = rxt(k,233)*y(k,88) + rxt(k,215)*y(k,89) + rxt(k,197)*y(k,92)
         mat(k,899) = rxt(k,243)*y(k,61) + rxt(k,247)*y(k,62) + rxt(k,235)*y(k,63) &
                      + rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) + rxt(k,244)*y(k,67) &
                      + rxt(k,239)*y(k,68) + rxt(k,237)*y(k,71) + rxt(k,249)*y(k,81) &
                      + rxt(k,245)*y(k,82) + rxt(k,236)*y(k,83) + rxt(k,248)*y(k,84) &
                      + rxt(k,240)*y(k,85) + rxt(k,246)*y(k,86) + rxt(k,233)*y(k,87) &
                      + rxt(k,238)*y(k,98)
         mat(k,942) = rxt(k,225)*y(k,61) + rxt(k,229)*y(k,62) + rxt(k,217)*y(k,63) &
                      + rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67) &
                      + rxt(k,222)*y(k,68) + rxt(k,298)*y(k,69) + rxt(k,219)*y(k,71) &
                      + rxt(k,232)*y(k,81) + rxt(k,227)*y(k,82) + rxt(k,218)*y(k,83) &
                      + rxt(k,230)*y(k,84) + rxt(k,223)*y(k,85) + rxt(k,228)*y(k,86) &
                      + rxt(k,215)*y(k,87) + rxt(k,221)*y(k,98)
         mat(k,49) = rxt(k,299)*y(k,69)
         mat(k,121) = rxt(k,301)*y(k,69)
         mat(k,1083) = rxt(k,304)*y(k,32) + rxt(k,207)*y(k,61) + rxt(k,212)*y(k,62) &
                      + rxt(k,200)*y(k,63) + rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) &
                      + rxt(k,208)*y(k,67) + rxt(k,204)*y(k,68) + rxt(k,303)*y(k,69) &
                      + rxt(k,202)*y(k,71) + rxt(k,214)*y(k,81) + rxt(k,210)*y(k,82) &
                      + rxt(k,201)*y(k,83) + rxt(k,213)*y(k,84) + rxt(k,205)*y(k,85) &
                      + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87) + rxt(k,203)*y(k,98)
         mat(k,114) = rxt(k,379)*y(k,69)
         mat(k,1325) = rxt(k,238)*y(k,88) + rxt(k,221)*y(k,89) + rxt(k,203)*y(k,92)
         mat(k,775) = mat(k,775) + rxt(k,433)*y(k,48)
         mat(k,1153) = -(rxt(k,105)*y(k,61) + (rxt(k,115) + rxt(k,116)) * y(k,67) &
                      + (rxt(k,272) + rxt(k,273)) * y(k,101) + rxt(k,276)*y(k,65) &
                      + rxt(k,286)*y(k,81) + rxt(k,315)*y(k,98) + rxt(k,341)*y(k,105) &
                      + rxt(k,354)*y(k,106) + (rxt(k,422) + rxt(k,423) + rxt(k,424) &
                      ) * y(k,48) + (rxt(k,426) + rxt(k,428)) * y(k,54) + rxt(k,427) &
                      *y(k,56) + rxt(k,439)*y(k,43) + rxt(k,440)*y(k,53) + rxt(k,441) &
                      *y(k,104) + rxt(k,459)*y(k,24) + rxt(k,490)*y(k,2))
         mat(k,606) = -rxt(k,105)*y(k,52)
         mat(k,507) = -(rxt(k,115) + rxt(k,116)) * y(k,52)
         mat(k,407) = -(rxt(k,272) + rxt(k,273)) * y(k,52)
         mat(k,793) = -rxt(k,276)*y(k,52)
         mat(k,721) = -rxt(k,286)*y(k,52)
         mat(k,1310) = -rxt(k,315)*y(k,52)
         mat(k,1599) = -rxt(k,341)*y(k,52)
         mat(k,1635) = -rxt(k,354)*y(k,52)
         mat(k,277) = -(rxt(k,422) + rxt(k,423) + rxt(k,424)) * y(k,52)
         mat(k,983) = -(rxt(k,426) + rxt(k,428)) * y(k,52)
         mat(k,1435) = -rxt(k,427)*y(k,52)
         mat(k,693) = -rxt(k,439)*y(k,52)
         mat(k,1231) = -rxt(k,440)*y(k,52)
         mat(k,762) = -rxt(k,441)*y(k,52)
         mat(k,1523) = -rxt(k,459)*y(k,52)
         mat(k,464) = -rxt(k,490)*y(k,52)
         mat(k,464) = mat(k,464) + rxt(k,489)*y(k,51)
         mat(k,261) = rxt(k,529)*y(k,51)
         mat(k,1523) = mat(k,1523) + rxt(k,457)*y(k,51)
         mat(k,693) = mat(k,693) + rxt(k,434)*y(k,51) + rxt(k,429)*y(k,53)
         mat(k,54) = rxt(k,418)*y(k,104)
         mat(k,125) = rxt(k,335)*y(k,104)
         mat(k,1756) = rxt(k,489)*y(k,2) + rxt(k,529)*y(k,18) + rxt(k,457)*y(k,24) &
                      + rxt(k,434)*y(k,43) + 2.000_r8*rxt(k,430)*y(k,53) + rxt(k,436) &
                      *y(k,54) + rxt(k,435)*y(k,56) + rxt(k,107)*y(k,64) + rxt(k,131) &
                      *y(k,106)
         mat(k,1231) = mat(k,1231) + rxt(k,429)*y(k,43) + 2.000_r8*rxt(k,430)*y(k,51) &
                      + rxt(k,431)*y(k,54) + rxt(k,432)*y(k,104)
         mat(k,983) = mat(k,983) + rxt(k,436)*y(k,51) + rxt(k,431)*y(k,53)
         mat(k,1435) = mat(k,1435) + rxt(k,435)*y(k,51)
         mat(k,1678) = rxt(k,280)*y(k,81)
         mat(k,215) = rxt(k,107)*y(k,51)
         mat(k,836) = rxt(k,196)*y(k,81) + rxt(k,192)*y(k,82)
         mat(k,1196) = rxt(k,161)*y(k,81) + rxt(k,157)*y(k,82)
         mat(k,1566) = rxt(k,179)*y(k,81) + rxt(k,174)*y(k,82)
         mat(k,721) = mat(k,721) + rxt(k,280)*y(k,60) + rxt(k,196)*y(k,73) &
                      + rxt(k,161)*y(k,75) + rxt(k,179)*y(k,77) + rxt(k,249)*y(k,88) &
                      + rxt(k,232)*y(k,89) + rxt(k,214)*y(k,92) + rxt(k,156)*y(k,99)
         mat(k,377) = rxt(k,192)*y(k,73) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77) &
                      + rxt(k,245)*y(k,88) + rxt(k,227)*y(k,89) + rxt(k,210)*y(k,92) &
                      + rxt(k,151)*y(k,99)
         mat(k,884) = rxt(k,249)*y(k,81) + rxt(k,245)*y(k,82)
         mat(k,927) = rxt(k,232)*y(k,81) + rxt(k,227)*y(k,82)
         mat(k,1068) = rxt(k,214)*y(k,81) + rxt(k,210)*y(k,82) + rxt(k,306)*y(k,104)
         mat(k,1355) = rxt(k,156)*y(k,81) + rxt(k,151)*y(k,82)
         mat(k,762) = mat(k,762) + rxt(k,418)*y(k,44) + rxt(k,335)*y(k,47) &
                      + rxt(k,432)*y(k,53) + rxt(k,306)*y(k,92)
         mat(k,1635) = mat(k,1635) + rxt(k,131)*y(k,51)
         mat(k,1233) = -(rxt(k,429)*y(k,43) + rxt(k,430)*y(k,51) + rxt(k,431)*y(k,54) &
                      + rxt(k,432)*y(k,104) + rxt(k,440)*y(k,52) + rxt(k,525)*y(k,14))
         mat(k,694) = -rxt(k,429)*y(k,53)
         mat(k,1758) = -rxt(k,430)*y(k,53)
         mat(k,985) = -rxt(k,431)*y(k,53)
         mat(k,763) = -rxt(k,432)*y(k,53)
         mat(k,1155) = -rxt(k,440)*y(k,53)
         mat(k,549) = -rxt(k,525)*y(k,53)
         mat(k,78) = rxt(k,491)*y(k,54)
         mat(k,1525) = rxt(k,281)*y(k,81)
         mat(k,206) = rxt(k,460)*y(k,54) + rxt(k,458)*y(k,60) + rxt(k,461)*y(k,104)
         mat(k,173) = rxt(k,505)*y(k,42)
         mat(k,1110) = rxt(k,505)*y(k,31) + rxt(k,442)*y(k,104)
         mat(k,694) = mat(k,694) + rxt(k,305)*y(k,92)
         mat(k,1155) = mat(k,1155) + rxt(k,428)*y(k,54) + rxt(k,427)*y(k,56)
         mat(k,985) = mat(k,985) + rxt(k,491)*y(k,3) + rxt(k,460)*y(k,25) + rxt(k,428) &
                      *y(k,52)
         mat(k,1437) = rxt(k,427)*y(k,52)
         mat(k,1680) = rxt(k,458)*y(k,25)
         mat(k,838) = rxt(k,195)*y(k,84) + rxt(k,188)*y(k,85) + rxt(k,193)*y(k,86)
         mat(k,1198) = rxt(k,160)*y(k,84) + (rxt(k,231)+rxt(k,257))*y(k,85) &
                      + rxt(k,158)*y(k,86)
         mat(k,1568) = rxt(k,178)*y(k,84) + (rxt(k,170)+rxt(k,258))*y(k,85) &
                      + rxt(k,175)*y(k,86)
         mat(k,723) = rxt(k,281)*y(k,24)
         mat(k,663) = rxt(k,236)*y(k,88) + rxt(k,218)*y(k,89) + rxt(k,201)*y(k,92) &
                      + rxt(k,143)*y(k,99)
         mat(k,488) = rxt(k,195)*y(k,73) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77) &
                      + rxt(k,248)*y(k,88) + rxt(k,230)*y(k,89) + rxt(k,213)*y(k,92) &
                      + rxt(k,155)*y(k,99)
         mat(k,529) = rxt(k,188)*y(k,73) + (rxt(k,231)+rxt(k,257))*y(k,75) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,77) + rxt(k,240)*y(k,88) &
                      + rxt(k,223)*y(k,89) + rxt(k,205)*y(k,92) + rxt(k,147)*y(k,99)
         mat(k,393) = rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77) &
                      + rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) + rxt(k,211)*y(k,92) &
                      + rxt(k,152)*y(k,99)
         mat(k,446) = rxt(k,233)*y(k,88) + rxt(k,215)*y(k,89) + rxt(k,197)*y(k,92) &
                      + rxt(k,250)*y(k,99)
         mat(k,886) = rxt(k,236)*y(k,83) + rxt(k,248)*y(k,84) + rxt(k,240)*y(k,85) &
                      + rxt(k,246)*y(k,86) + rxt(k,233)*y(k,87)
         mat(k,929) = rxt(k,218)*y(k,83) + rxt(k,230)*y(k,84) + rxt(k,223)*y(k,85) &
                      + rxt(k,228)*y(k,86) + rxt(k,215)*y(k,87)
         mat(k,1070) = rxt(k,305)*y(k,43) + rxt(k,201)*y(k,83) + rxt(k,213)*y(k,84) &
                      + rxt(k,205)*y(k,85) + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87)
         mat(k,1357) = rxt(k,143)*y(k,83) + rxt(k,155)*y(k,84) + rxt(k,147)*y(k,85) &
                      + rxt(k,152)*y(k,86) + rxt(k,250)*y(k,87)
         mat(k,763) = mat(k,763) + rxt(k,461)*y(k,25) + rxt(k,442)*y(k,42)
         mat(k,979) = -(rxt(k,108)*y(k,64) + rxt(k,120)*y(k,65) + rxt(k,126)*y(k,68) &
                      + rxt(k,296)*y(k,83) + (rxt(k,319) + rxt(k,320)) * y(k,98) &
                      + (rxt(k,328) + rxt(k,329)) * y(k,101) + rxt(k,331)*y(k,102) &
                      + rxt(k,333)*y(k,103) + rxt(k,342)*y(k,105) + rxt(k,355) &
                      *y(k,106) + rxt(k,398)*y(k,56) + 4._r8*rxt(k,399)*y(k,54) &
                      + rxt(k,400)*y(k,55) + rxt(k,401)*y(k,33) + rxt(k,402)*y(k,35) &
                      + rxt(k,407)*y(k,43) + rxt(k,413)*y(k,104) + (rxt(k,426) &
                      + rxt(k,428)) * y(k,52) + rxt(k,431)*y(k,53) + rxt(k,436) &
                      *y(k,51) + rxt(k,460)*y(k,25) + rxt(k,462)*y(k,24) + rxt(k,465) &
                      *y(k,40) + rxt(k,468)*y(k,46) + rxt(k,491)*y(k,3) + rxt(k,492) &
                      *y(k,2) + rxt(k,494)*y(k,36) + rxt(k,496)*y(k,45) + rxt(k,526) &
                      *y(k,14) + (rxt(k,565) + rxt(k,566)) * y(k,80) + rxt(k,567) &
                      *y(k,94))
         mat(k,214) = -rxt(k,108)*y(k,54)
         mat(k,789) = -rxt(k,120)*y(k,54)
         mat(k,628) = -rxt(k,126)*y(k,54)
         mat(k,657) = -rxt(k,296)*y(k,54)
         mat(k,1306) = -(rxt(k,319) + rxt(k,320)) * y(k,54)
         mat(k,405) = -(rxt(k,328) + rxt(k,329)) * y(k,54)
         mat(k,98) = -rxt(k,331)*y(k,54)
         mat(k,288) = -rxt(k,333)*y(k,54)
         mat(k,1595) = -rxt(k,342)*y(k,54)
         mat(k,1631) = -rxt(k,355)*y(k,54)
         mat(k,1431) = -rxt(k,398)*y(k,54)
         mat(k,1270) = -rxt(k,400)*y(k,54)
         mat(k,424) = -rxt(k,401)*y(k,54)
         mat(k,82) = -rxt(k,402)*y(k,54)
         mat(k,689) = -rxt(k,407)*y(k,54)
         mat(k,758) = -rxt(k,413)*y(k,54)
         mat(k,1149) = -(rxt(k,426) + rxt(k,428)) * y(k,54)
         mat(k,1227) = -rxt(k,431)*y(k,54)
         mat(k,1752) = -rxt(k,436)*y(k,54)
         mat(k,203) = -rxt(k,460)*y(k,54)
         mat(k,1519) = -rxt(k,462)*y(k,54)
         mat(k,1476) = -rxt(k,465)*y(k,54)
         mat(k,186) = -rxt(k,468)*y(k,54)
         mat(k,75) = -rxt(k,491)*y(k,54)
         mat(k,462) = -rxt(k,492)*y(k,54)
         mat(k,196) = -rxt(k,494)*y(k,54)
         mat(k,155) = -rxt(k,496)*y(k,54)
         mat(k,545) = -rxt(k,526)*y(k,54)
         mat(k,91) = -(rxt(k,565) + rxt(k,566)) * y(k,54)
         mat(k,180) = -rxt(k,567)*y(k,54)
         mat(k,1710) = rxt(k,405)*y(k,43)
         mat(k,689) = mat(k,689) + rxt(k,405)*y(k,32)
         mat(k,275) = rxt(k,421)*y(k,51) + rxt(k,422)*y(k,52) + rxt(k,425)*y(k,55) &
                      + rxt(k,570)*y(k,99)
         mat(k,1752) = mat(k,1752) + rxt(k,421)*y(k,48) + rxt(k,267)*y(k,101)
         mat(k,1149) = mat(k,1149) + rxt(k,422)*y(k,48) + rxt(k,354)*y(k,106)
         mat(k,1270) = mat(k,1270) + rxt(k,425)*y(k,48) + rxt(k,569)*y(k,94) + ( &
                      + rxt(k,387)+rxt(k,388))*y(k,95) + rxt(k,576)*y(k,107) &
                      + rxt(k,580)*y(k,108)
         mat(k,1431) = mat(k,1431) + rxt(k,358)*y(k,106)
         mat(k,1674) = rxt(k,109)*y(k,65) + rxt(k,345)*y(k,106)
         mat(k,789) = mat(k,789) + rxt(k,109)*y(k,60) + rxt(k,181)*y(k,73) + ( &
                      + rxt(k,165)+rxt(k,253))*y(k,75) + (rxt(k,163)+rxt(k,260)) &
                      *y(k,77) + rxt(k,234)*y(k,88) + rxt(k,216)*y(k,89) + rxt(k,199) &
                      *y(k,92) + rxt(k,251)*y(k,99)
         mat(k,316) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77) + rxt(k,241)*y(k,88) &
                      + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92) + rxt(k,148)*y(k,99)
         mat(k,505) = rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + rxt(k,226)*y(k,89) + rxt(k,208)*y(k,92) + rxt(k,150)*y(k,99)
         mat(k,1022) = rxt(k,561)*y(k,88) + 1.150_r8*rxt(k,562)*y(k,99)
         mat(k,832) = rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67)
         mat(k,1192) = (rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67)
         mat(k,1562) = (rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67)
         mat(k,164) = rxt(k,575)*y(k,107)
         mat(k,880) = rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) + rxt(k,244)*y(k,67) &
                      + rxt(k,561)*y(k,69)
         mat(k,923) = rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67)
         mat(k,1064) = rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) + rxt(k,208)*y(k,67)
         mat(k,180) = mat(k,180) + rxt(k,569)*y(k,55)
         mat(k,577) = (rxt(k,387)+rxt(k,388))*y(k,55)
         mat(k,1351) = rxt(k,570)*y(k,48) + rxt(k,251)*y(k,65) + rxt(k,148)*y(k,66) &
                      + rxt(k,150)*y(k,67) + 1.150_r8*rxt(k,562)*y(k,69)
         mat(k,405) = mat(k,405) + rxt(k,267)*y(k,51)
         mat(k,758) = mat(k,758) + 2.000_r8*rxt(k,415)*y(k,104)
         mat(k,1631) = mat(k,1631) + rxt(k,354)*y(k,52) + rxt(k,358)*y(k,56) &
                      + rxt(k,345)*y(k,60)
         mat(k,250) = rxt(k,576)*y(k,55) + rxt(k,575)*y(k,79)
         mat(k,68) = rxt(k,580)*y(k,55)
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
         mat(k,1277) = -(rxt(k,121)*y(k,65) + (rxt(k,128) + rxt(k,130)) * y(k,69) &
                      + rxt(k,317)*y(k,98) + rxt(k,357)*y(k,106) + rxt(k,359)*y(k,99) &
                      + rxt(k,387)*y(k,95) + rxt(k,392)*y(k,96) + rxt(k,400)*y(k,54) &
                      + rxt(k,406)*y(k,32) + rxt(k,420)*y(k,79) + rxt(k,425)*y(k,48) &
                      + rxt(k,564)*y(k,80) + (rxt(k,568) + rxt(k,569)) * y(k,94) &
                      + rxt(k,576)*y(k,107) + rxt(k,580)*y(k,108))
         mat(k,796) = -rxt(k,121)*y(k,55)
         mat(k,1029) = -(rxt(k,128) + rxt(k,130)) * y(k,55)
         mat(k,1313) = -rxt(k,317)*y(k,55)
         mat(k,1638) = -rxt(k,357)*y(k,55)
         mat(k,1358) = -rxt(k,359)*y(k,55)
         mat(k,582) = -rxt(k,387)*y(k,55)
         mat(k,222) = -rxt(k,392)*y(k,55)
         mat(k,986) = -rxt(k,400)*y(k,55)
         mat(k,1717) = -rxt(k,406)*y(k,55)
         mat(k,165) = -rxt(k,420)*y(k,55)
         mat(k,278) = -rxt(k,425)*y(k,55)
         mat(k,93) = -rxt(k,564)*y(k,55)
         mat(k,181) = -(rxt(k,568) + rxt(k,569)) * y(k,55)
         mat(k,251) = -rxt(k,576)*y(k,55)
         mat(k,69) = -rxt(k,580)*y(k,55)
         mat(k,466) = 2.000_r8*rxt(k,484)*y(k,2) + (rxt(k,486)+rxt(k,487))*y(k,24) &
                      + rxt(k,488)*y(k,43) + rxt(k,492)*y(k,54)
         mat(k,262) = rxt(k,528)*y(k,43)
         mat(k,1526) = (rxt(k,486)+rxt(k,487))*y(k,2) + (2.000_r8*rxt(k,453) &
                       +2.000_r8*rxt(k,454))*y(k,24) + rxt(k,456)*y(k,43) + rxt(k,462) &
                      *y(k,54) + rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68) + rxt(k,310) &
                      *y(k,98) + rxt(k,464)*y(k,104) + rxt(k,346)*y(k,106)
         mat(k,1400) = rxt(k,325)*y(k,101) + rxt(k,330)*y(k,102)
         mat(k,1717) = mat(k,1717) + rxt(k,403)*y(k,43) + rxt(k,409)*y(k,56) &
                      + rxt(k,326)*y(k,101)
         mat(k,695) = rxt(k,488)*y(k,2) + rxt(k,528)*y(k,18) + rxt(k,456)*y(k,24) &
                      + rxt(k,403)*y(k,32) + 2.000_r8*rxt(k,417)*y(k,43) + rxt(k,429) &
                      *y(k,53) + rxt(k,407)*y(k,54) + 2.000_r8*rxt(k,408)*y(k,56) &
                      + rxt(k,482)*y(k,59) + rxt(k,449)*y(k,60) + rxt(k,412)*y(k,104)
         mat(k,56) = rxt(k,418)*y(k,104)
         mat(k,278) = mat(k,278) + rxt(k,424)*y(k,52)
         mat(k,1759) = rxt(k,435)*y(k,56) + rxt(k,572)*y(k,99) + rxt(k,274)*y(k,101)
         mat(k,1156) = rxt(k,424)*y(k,48) + rxt(k,426)*y(k,54) + rxt(k,427)*y(k,56) &
                      + rxt(k,315)*y(k,98) + rxt(k,272)*y(k,101)
         mat(k,1234) = rxt(k,429)*y(k,43) + rxt(k,431)*y(k,54)
         mat(k,986) = mat(k,986) + rxt(k,492)*y(k,2) + rxt(k,462)*y(k,24) + rxt(k,407) &
                      *y(k,43) + rxt(k,426)*y(k,52) + rxt(k,431)*y(k,53) &
                      + 2.000_r8*rxt(k,399)*y(k,54) + 2.000_r8*rxt(k,398)*y(k,56) &
                      + rxt(k,108)*y(k,64) + rxt(k,126)*y(k,68) + rxt(k,296)*y(k,83) &
                      + rxt(k,391)*y(k,96) + rxt(k,320)*y(k,98) + (2.000_r8*rxt(k,328) &
                       +rxt(k,329))*y(k,101) + rxt(k,331)*y(k,102) + rxt(k,413) &
                      *y(k,104) + rxt(k,355)*y(k,106)
         mat(k,1277) = mat(k,1277) + 2.000_r8*rxt(k,392)*y(k,96)
         mat(k,1438) = rxt(k,409)*y(k,32) + 2.000_r8*rxt(k,408)*y(k,43) + rxt(k,435) &
                      *y(k,51) + rxt(k,427)*y(k,52) + 2.000_r8*rxt(k,398)*y(k,54) &
                      + rxt(k,483)*y(k,59) + rxt(k,451)*y(k,60) + rxt(k,127)*y(k,68) &
                      + rxt(k,129)*y(k,69) + rxt(k,287)*y(k,81) + 2.000_r8*rxt(k,297) &
                      *y(k,83) + 2.000_r8*rxt(k,389)*y(k,95) + rxt(k,318)*y(k,98) &
                      + 3.000_r8*rxt(k,327)*y(k,101) + rxt(k,414)*y(k,104)
         mat(k,335) = rxt(k,482)*y(k,43) + rxt(k,483)*y(k,56)
         mat(k,1681) = rxt(k,449)*y(k,43) + rxt(k,451)*y(k,56) + rxt(k,122)*y(k,68) &
                      + rxt(k,309)*y(k,98)
         mat(k,608) = rxt(k,149)*y(k,99)
         mat(k,363) = rxt(k,154)*y(k,99)
         mat(k,347) = rxt(k,252)*y(k,99)
         mat(k,216) = rxt(k,108)*y(k,54)
         mat(k,796) = mat(k,796) + rxt(k,111)*y(k,24) + rxt(k,251)*y(k,99)
         mat(k,319) = rxt(k,148)*y(k,99)
         mat(k,509) = rxt(k,150)*y(k,99)
         mat(k,634) = rxt(k,123)*y(k,24) + rxt(k,126)*y(k,54) + rxt(k,127)*y(k,56) &
                      + rxt(k,122)*y(k,60) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77) + rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92) + 2.000_r8*rxt(k,146)*y(k,99)
         mat(k,1029) = mat(k,1029) + rxt(k,129)*y(k,56) + rxt(k,321)*y(k,100) &
                      + 2.000_r8*rxt(k,375)*y(k,103)
         mat(k,304) = rxt(k,144)*y(k,99)
         mat(k,839) = rxt(k,186)*y(k,68) + rxt(k,185)*y(k,98)
         mat(k,1199) = rxt(k,220)*y(k,68) + rxt(k,209)*y(k,98)
         mat(k,1569) = rxt(k,169)*y(k,68) + rxt(k,168)*y(k,98)
         mat(k,724) = rxt(k,287)*y(k,56) + rxt(k,156)*y(k,99)
         mat(k,379) = rxt(k,151)*y(k,99)
         mat(k,664) = rxt(k,296)*y(k,54) + 2.000_r8*rxt(k,297)*y(k,56) + rxt(k,143) &
                      *y(k,99)
         mat(k,489) = rxt(k,155)*y(k,99)
         mat(k,530) = rxt(k,147)*y(k,99)
         mat(k,394) = rxt(k,152)*y(k,99)
         mat(k,447) = rxt(k,250)*y(k,99)
         mat(k,887) = rxt(k,239)*y(k,68) + rxt(k,238)*y(k,98)
         mat(k,930) = rxt(k,222)*y(k,68) + rxt(k,221)*y(k,98)
         mat(k,1071) = rxt(k,204)*y(k,68) + rxt(k,203)*y(k,98)
         mat(k,582) = mat(k,582) + 2.000_r8*rxt(k,389)*y(k,56)
         mat(k,222) = mat(k,222) + rxt(k,391)*y(k,54) + 2.000_r8*rxt(k,392)*y(k,55) &
                      + 2.000_r8*rxt(k,316)*y(k,98) + 2.000_r8*rxt(k,334)*y(k,103)
         mat(k,1313) = mat(k,1313) + rxt(k,310)*y(k,24) + rxt(k,315)*y(k,52) &
                      + rxt(k,320)*y(k,54) + rxt(k,318)*y(k,56) + rxt(k,309)*y(k,60) &
                      + rxt(k,185)*y(k,73) + rxt(k,209)*y(k,75) + rxt(k,168)*y(k,77) &
                      + rxt(k,238)*y(k,88) + rxt(k,221)*y(k,89) + rxt(k,203)*y(k,92) &
                      + 2.000_r8*rxt(k,316)*y(k,96)
         mat(k,1358) = mat(k,1358) + rxt(k,572)*y(k,51) + rxt(k,149)*y(k,61) &
                      + rxt(k,154)*y(k,62) + rxt(k,252)*y(k,63) + rxt(k,251)*y(k,65) &
                      + rxt(k,148)*y(k,66) + rxt(k,150)*y(k,67) + 2.000_r8*rxt(k,146) &
                      *y(k,68) + rxt(k,144)*y(k,71) + rxt(k,156)*y(k,81) + rxt(k,151) &
                      *y(k,82) + rxt(k,143)*y(k,83) + rxt(k,155)*y(k,84) + rxt(k,147) &
                      *y(k,85) + rxt(k,152)*y(k,86) + rxt(k,250)*y(k,87)
         mat(k,146) = rxt(k,321)*y(k,69) + (rxt(k,322)+rxt(k,323))*y(k,110)
         mat(k,408) = rxt(k,325)*y(k,28) + rxt(k,326)*y(k,32) + rxt(k,274)*y(k,51) &
                      + rxt(k,272)*y(k,52) + (2.000_r8*rxt(k,328)+rxt(k,329))*y(k,54) &
                      + 3.000_r8*rxt(k,327)*y(k,56)
         mat(k,99) = rxt(k,330)*y(k,28) + rxt(k,331)*y(k,54)
         mat(k,290) = 2.000_r8*rxt(k,375)*y(k,69) + 2.000_r8*rxt(k,334)*y(k,96) &
                      + rxt(k,332)*y(k,110)
         mat(k,764) = rxt(k,464)*y(k,24) + rxt(k,412)*y(k,43) + rxt(k,418)*y(k,44) &
                      + rxt(k,413)*y(k,54) + rxt(k,414)*y(k,56)
         mat(k,1638) = mat(k,1638) + rxt(k,346)*y(k,24) + rxt(k,355)*y(k,54)
         mat(k,1812) = (rxt(k,322)+rxt(k,323))*y(k,100) + rxt(k,332)*y(k,103)
         mat(k,1442) = -(rxt(k,127)*y(k,68) + rxt(k,129)*y(k,69) + rxt(k,287)*y(k,81) &
                      + rxt(k,297)*y(k,83) + rxt(k,318)*y(k,98) + rxt(k,327)*y(k,101) &
                      + rxt(k,343)*y(k,105) + rxt(k,358)*y(k,106) + rxt(k,389)*y(k,95) &
                      + rxt(k,398)*y(k,54) + rxt(k,408)*y(k,43) + rxt(k,409)*y(k,32) &
                      + rxt(k,414)*y(k,104) + rxt(k,427)*y(k,52) + rxt(k,435)*y(k,51) &
                      + rxt(k,451)*y(k,60) + rxt(k,483)*y(k,59))
         mat(k,638) = -rxt(k,127)*y(k,56)
         mat(k,1033) = -rxt(k,129)*y(k,56)
         mat(k,726) = -rxt(k,287)*y(k,56)
         mat(k,666) = -rxt(k,297)*y(k,56)
         mat(k,1317) = -rxt(k,318)*y(k,56)
         mat(k,411) = -rxt(k,327)*y(k,56)
         mat(k,1606) = -rxt(k,343)*y(k,56)
         mat(k,1642) = -rxt(k,358)*y(k,56)
         mat(k,586) = -rxt(k,389)*y(k,56)
         mat(k,990) = -rxt(k,398)*y(k,56)
         mat(k,699) = -rxt(k,408)*y(k,56)
         mat(k,1721) = -rxt(k,409)*y(k,56)
         mat(k,768) = -rxt(k,414)*y(k,56)
         mat(k,1160) = -rxt(k,427)*y(k,56)
         mat(k,1763) = -rxt(k,435)*y(k,56)
         mat(k,1685) = -rxt(k,451)*y(k,56)
         mat(k,337) = -rxt(k,483)*y(k,56)
         mat(k,1160) = mat(k,1160) + rxt(k,273)*y(k,101)
         mat(k,990) = mat(k,990) + rxt(k,400)*y(k,55) + rxt(k,319)*y(k,98) &
                      + rxt(k,333)*y(k,103)
         mat(k,1281) = rxt(k,400)*y(k,54)
         mat(k,225) = rxt(k,356)*y(k,106)
         mat(k,1317) = mat(k,1317) + rxt(k,319)*y(k,54)
         mat(k,411) = mat(k,411) + rxt(k,273)*y(k,52)
         mat(k,293) = rxt(k,333)*y(k,54)
         mat(k,1642) = mat(k,1642) + rxt(k,356)*y(k,96)
         mat(k,452) = rxt(k,485)*y(k,24)
         mat(k,1498) = rxt(k,485)*y(k,2) + 2.000_r8*rxt(k,455)*y(k,24)
         mat(k,328) = -(rxt(k,481)*y(k,14) + rxt(k,482)*y(k,43) + rxt(k,483)*y(k,56))
         mat(k,538) = -rxt(k,481)*y(k,59)
         mat(k,681) = -rxt(k,482)*y(k,59)
         mat(k,1417) = -rxt(k,483)*y(k,59)
         mat(k,456) = 4.000_r8*rxt(k,484)*y(k,2) + (rxt(k,485)+rxt(k,486))*y(k,24) &
                      + rxt(k,489)*y(k,51) + rxt(k,492)*y(k,54) + rxt(k,493)*y(k,104)
         mat(k,1505) = (rxt(k,485)+rxt(k,486))*y(k,2)
         mat(k,192) = rxt(k,494)*y(k,54) + rxt(k,500)*y(k,95) + rxt(k,495)*y(k,104)
         mat(k,1735) = rxt(k,489)*y(k,2)
         mat(k,963) = rxt(k,492)*y(k,2) + rxt(k,494)*y(k,36)
         mat(k,569) = rxt(k,500)*y(k,36)
         mat(k,749) = rxt(k,493)*y(k,2) + rxt(k,495)*y(k,36)
         mat(k,1691) = -((rxt(k,109) + rxt(k,110)) * y(k,65) + rxt(k,122)*y(k,68) &
                      + rxt(k,280)*y(k,81) + rxt(k,309)*y(k,98) + rxt(k,336)*y(k,105) &
                      + rxt(k,345)*y(k,106) + rxt(k,445)*y(k,14) + rxt(k,447)*y(k,33) &
                      + rxt(k,448)*y(k,35) + (rxt(k,449) + rxt(k,450)) * y(k,43) &
                      + rxt(k,451)*y(k,56) + rxt(k,458)*y(k,25) + rxt(k,467)*y(k,46))
         mat(k,806) = -(rxt(k,109) + rxt(k,110)) * y(k,60)
         mat(k,644) = -rxt(k,122)*y(k,60)
         mat(k,731) = -rxt(k,280)*y(k,60)
         mat(k,1323) = -rxt(k,309)*y(k,60)
         mat(k,1612) = -rxt(k,336)*y(k,60)
         mat(k,1648) = -rxt(k,345)*y(k,60)
         mat(k,557) = -rxt(k,445)*y(k,60)
         mat(k,434) = -rxt(k,447)*y(k,60)
         mat(k,84) = -rxt(k,448)*y(k,60)
         mat(k,704) = -(rxt(k,449) + rxt(k,450)) * y(k,60)
         mat(k,1448) = -rxt(k,451)*y(k,60)
         mat(k,209) = -rxt(k,458)*y(k,60)
         mat(k,189) = -rxt(k,467)*y(k,60)
         mat(k,471) = rxt(k,486)*y(k,24)
         mat(k,264) = rxt(k,452)*y(k,24)
         mat(k,1536) = rxt(k,486)*y(k,2) + rxt(k,452)*y(k,18) + (4.000_r8*rxt(k,453) &
                       +2.000_r8*rxt(k,455))*y(k,24) + rxt(k,457)*y(k,51) + rxt(k,462) &
                      *y(k,54) + rxt(k,463)*y(k,104)
         mat(k,15) = rxt(k,507)*y(k,95)
         mat(k,1493) = rxt(k,465)*y(k,54) + rxt(k,477)*y(k,95) + rxt(k,466)*y(k,104)
         mat(k,1769) = rxt(k,457)*y(k,24) + rxt(k,106)*y(k,64)
         mat(k,1166) = rxt(k,105)*y(k,61)
         mat(k,996) = rxt(k,462)*y(k,24) + rxt(k,465)*y(k,40)
         mat(k,612) = rxt(k,105)*y(k,52) + rxt(k,190)*y(k,73) + rxt(k,142)*y(k,75) &
                      + rxt(k,172)*y(k,77) + rxt(k,243)*y(k,88) + rxt(k,225)*y(k,89) &
                      + rxt(k,207)*y(k,92) + rxt(k,149)*y(k,99)
         mat(k,367) = rxt(k,194)*y(k,73) + rxt(k,159)*y(k,75) + rxt(k,177)*y(k,77) &
                      + rxt(k,247)*y(k,88) + rxt(k,229)*y(k,89) + rxt(k,212)*y(k,92) &
                      + rxt(k,154)*y(k,99)
         mat(k,351) = rxt(k,182)*y(k,73) + rxt(k,176)*y(k,75) + rxt(k,164)*y(k,77) &
                      + rxt(k,235)*y(k,88) + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92) &
                      + rxt(k,252)*y(k,99)
         mat(k,217) = rxt(k,106)*y(k,51)
         mat(k,849) = rxt(k,190)*y(k,61) + rxt(k,194)*y(k,62) + rxt(k,182)*y(k,63)
         mat(k,1209) = rxt(k,142)*y(k,61) + rxt(k,159)*y(k,62) + rxt(k,176)*y(k,63)
         mat(k,1579) = rxt(k,172)*y(k,61) + rxt(k,177)*y(k,62) + rxt(k,164)*y(k,63)
         mat(k,897) = rxt(k,243)*y(k,61) + rxt(k,247)*y(k,62) + rxt(k,235)*y(k,63)
         mat(k,940) = rxt(k,225)*y(k,61) + rxt(k,229)*y(k,62) + rxt(k,217)*y(k,63)
         mat(k,1081) = rxt(k,207)*y(k,61) + rxt(k,212)*y(k,62) + rxt(k,200)*y(k,63)
         mat(k,591) = rxt(k,507)*y(k,30) + rxt(k,477)*y(k,40)
         mat(k,1368) = rxt(k,149)*y(k,61) + rxt(k,154)*y(k,62) + rxt(k,252)*y(k,63)
         mat(k,773) = rxt(k,463)*y(k,24) + rxt(k,466)*y(k,40)
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
         mat(k,597) = -(rxt(k,100)*y(k,32) + rxt(k,102)*y(k,110) + rxt(k,103)*y(k,40) &
                      + rxt(k,104)*y(k,42) + rxt(k,105)*y(k,52) + rxt(k,142)*y(k,75) &
                      + rxt(k,149)*y(k,99) + rxt(k,172)*y(k,77) + rxt(k,190)*y(k,73) &
                      + rxt(k,207)*y(k,92) + rxt(k,225)*y(k,89) + rxt(k,243)*y(k,88))
         mat(k,1700) = -rxt(k,100)*y(k,61)
         mat(k,1796) = -rxt(k,102)*y(k,61)
         mat(k,1466) = -rxt(k,103)*y(k,61)
         mat(k,1095) = -rxt(k,104)*y(k,61)
         mat(k,1140) = -rxt(k,105)*y(k,61)
         mat(k,1182) = -rxt(k,142)*y(k,61)
         mat(k,1341) = -rxt(k,149)*y(k,61)
         mat(k,1552) = -rxt(k,172)*y(k,61)
         mat(k,822) = -rxt(k,190)*y(k,61)
         mat(k,1054) = -rxt(k,207)*y(k,61)
         mat(k,913) = -rxt(k,225)*y(k,61)
         mat(k,870) = -rxt(k,243)*y(k,61)
         mat(k,1509) = rxt(k,111)*y(k,65) + rxt(k,281)*y(k,81) + rxt(k,346)*y(k,106)
         mat(k,1466) = mat(k,1466) + rxt(k,125)*y(k,68) + rxt(k,284)*y(k,81) &
                      + rxt(k,293)*y(k,83) + rxt(k,313)*y(k,98) + rxt(k,340)*y(k,105) &
                      + rxt(k,351)*y(k,106)
         mat(k,1743) = rxt(k,107)*y(k,64)
         mat(k,969) = rxt(k,108)*y(k,64)
         mat(k,1664) = rxt(k,109)*y(k,65) + rxt(k,122)*y(k,68) + rxt(k,280)*y(k,81) &
                      + rxt(k,309)*y(k,98) + rxt(k,336)*y(k,105) + rxt(k,345)*y(k,106)
         mat(k,212) = rxt(k,107)*y(k,51) + rxt(k,108)*y(k,54)
         mat(k,781) = rxt(k,111)*y(k,24) + rxt(k,109)*y(k,60)
         mat(k,618) = rxt(k,125)*y(k,40) + rxt(k,122)*y(k,60)
         mat(k,710) = rxt(k,281)*y(k,24) + rxt(k,284)*y(k,40) + rxt(k,280)*y(k,60)
         mat(k,651) = rxt(k,293)*y(k,40)
         mat(k,1296) = rxt(k,313)*y(k,40) + rxt(k,309)*y(k,60)
         mat(k,1586) = rxt(k,340)*y(k,40) + rxt(k,336)*y(k,60)
         mat(k,1622) = rxt(k,346)*y(k,24) + rxt(k,351)*y(k,40) + rxt(k,345)*y(k,60)
         mat(k,356) = -(rxt(k,101)*y(k,40) + rxt(k,154)*y(k,99) + rxt(k,159)*y(k,75) &
                      + rxt(k,177)*y(k,77) + rxt(k,194)*y(k,73) + rxt(k,212)*y(k,92) &
                      + rxt(k,229)*y(k,89) + rxt(k,247)*y(k,88))
         mat(k,1460) = -rxt(k,101)*y(k,62)
         mat(k,1333) = -rxt(k,154)*y(k,62)
         mat(k,1175) = -rxt(k,159)*y(k,62)
         mat(k,1545) = -rxt(k,177)*y(k,62)
         mat(k,815) = -rxt(k,194)*y(k,62)
         mat(k,1047) = -rxt(k,212)*y(k,62)
         mat(k,906) = -rxt(k,229)*y(k,62)
         mat(k,862) = -rxt(k,247)*y(k,62)
         mat(k,596) = rxt(k,102)*y(k,110)
         mat(k,1788) = rxt(k,102)*y(k,61)
         mat(k,340) = -((rxt(k,164) + rxt(k,263)) * y(k,77) + (rxt(k,176) + rxt(k,262) &
                      ) * y(k,75) + rxt(k,182)*y(k,73) + rxt(k,200)*y(k,92) + rxt(k,217) &
                      *y(k,89) + rxt(k,235)*y(k,88) + rxt(k,252)*y(k,99))
         mat(k,1544) = -(rxt(k,164) + rxt(k,263)) * y(k,63)
         mat(k,1174) = -(rxt(k,176) + rxt(k,262)) * y(k,63)
         mat(k,814) = -rxt(k,182)*y(k,63)
         mat(k,1046) = -rxt(k,200)*y(k,63)
         mat(k,905) = -rxt(k,217)*y(k,63)
         mat(k,861) = -rxt(k,235)*y(k,63)
         mat(k,1332) = -rxt(k,252)*y(k,63)
         mat(k,1459) = rxt(k,103)*y(k,61) + rxt(k,101)*y(k,62)
         mat(k,595) = rxt(k,103)*y(k,40)
         mat(k,355) = rxt(k,101)*y(k,40)
         mat(k,211) = -((rxt(k,106) + rxt(k,107)) * y(k,51) + rxt(k,108)*y(k,54))
         mat(k,1732) = -(rxt(k,106) + rxt(k,107)) * y(k,64)
         mat(k,958) = -rxt(k,108)*y(k,64)
         mat(k,1503) = rxt(k,123)*y(k,68) + rxt(k,310)*y(k,98) + rxt(k,337)*y(k,105)
         mat(k,1657) = rxt(k,110)*y(k,65)
         mat(k,777) = rxt(k,110)*y(k,60)
         mat(k,616) = rxt(k,123)*y(k,24)
         mat(k,1292) = rxt(k,310)*y(k,24)
         mat(k,1583) = rxt(k,337)*y(k,24)
         mat(k,785) = -((rxt(k,109) + rxt(k,110)) * y(k,60) + rxt(k,111)*y(k,24) &
                      + rxt(k,112)*y(k,32) + rxt(k,114)*y(k,110) + rxt(k,119)*y(k,42) &
                      + rxt(k,120)*y(k,54) + rxt(k,121)*y(k,55) + (rxt(k,163) &
                      + rxt(k,260)) * y(k,77) + (rxt(k,165) + rxt(k,253)) * y(k,75) &
                      + rxt(k,181)*y(k,73) + rxt(k,199)*y(k,92) + rxt(k,216)*y(k,89) &
                      + rxt(k,234)*y(k,88) + rxt(k,251)*y(k,99) + rxt(k,275)*y(k,51) &
                      + rxt(k,276)*y(k,52))
         mat(k,1670) = -(rxt(k,109) + rxt(k,110)) * y(k,65)
         mat(k,1515) = -rxt(k,111)*y(k,65)
         mat(k,1706) = -rxt(k,112)*y(k,65)
         mat(k,1801) = -rxt(k,114)*y(k,65)
         mat(k,1100) = -rxt(k,119)*y(k,65)
         mat(k,975) = -rxt(k,120)*y(k,65)
         mat(k,1266) = -rxt(k,121)*y(k,65)
         mat(k,1558) = -(rxt(k,163) + rxt(k,260)) * y(k,65)
         mat(k,1188) = -(rxt(k,165) + rxt(k,253)) * y(k,65)
         mat(k,828) = -rxt(k,181)*y(k,65)
         mat(k,1060) = -rxt(k,199)*y(k,65)
         mat(k,919) = -rxt(k,216)*y(k,65)
         mat(k,876) = -rxt(k,234)*y(k,65)
         mat(k,1347) = -rxt(k,251)*y(k,65)
         mat(k,1748) = -rxt(k,275)*y(k,65)
         mat(k,1145) = -rxt(k,276)*y(k,65)
         mat(k,1389) = rxt(k,325)*y(k,101) + rxt(k,347)*y(k,106)
         mat(k,1706) = mat(k,1706) + rxt(k,124)*y(k,68)
         mat(k,975) = mat(k,975) + rxt(k,126)*y(k,68)
         mat(k,624) = rxt(k,124)*y(k,32) + rxt(k,126)*y(k,54)
         mat(k,404) = rxt(k,325)*y(k,28)
         mat(k,1627) = rxt(k,347)*y(k,28)
         mat(k,311) = -(rxt(k,148)*y(k,99) + (rxt(k,171) + rxt(k,254)) * y(k,77) &
                      + rxt(k,189)*y(k,73) + rxt(k,206)*y(k,92) + rxt(k,224)*y(k,89) &
                      + rxt(k,241)*y(k,88) + (rxt(k,242) + rxt(k,266)) * y(k,75))
         mat(k,1331) = -rxt(k,148)*y(k,66)
         mat(k,1543) = -(rxt(k,171) + rxt(k,254)) * y(k,66)
         mat(k,813) = -rxt(k,189)*y(k,66)
         mat(k,1045) = -rxt(k,206)*y(k,66)
         mat(k,904) = -rxt(k,224)*y(k,66)
         mat(k,860) = -rxt(k,241)*y(k,66)
         mat(k,1173) = -(rxt(k,242) + rxt(k,266)) * y(k,66)
         mat(k,495) = rxt(k,113)*y(k,110)
         mat(k,1787) = rxt(k,113)*y(k,67)
         mat(k,497) = -(rxt(k,113)*y(k,110) + (rxt(k,115) + rxt(k,116)) * y(k,52) &
                      + (rxt(k,117) + rxt(k,118)) * y(k,51) + rxt(k,150)*y(k,99) &
                      + (rxt(k,153) + rxt(k,255)) * y(k,75) + (rxt(k,173) + rxt(k,256) &
                      ) * y(k,77) + rxt(k,191)*y(k,73) + rxt(k,208)*y(k,92) + rxt(k,226) &
                      *y(k,89) + rxt(k,244)*y(k,88))
         mat(k,1792) = -rxt(k,113)*y(k,67)
         mat(k,1136) = -(rxt(k,115) + rxt(k,116)) * y(k,67)
         mat(k,1739) = -(rxt(k,117) + rxt(k,118)) * y(k,67)
         mat(k,1338) = -rxt(k,150)*y(k,67)
         mat(k,1180) = -(rxt(k,153) + rxt(k,255)) * y(k,67)
         mat(k,1550) = -(rxt(k,173) + rxt(k,256)) * y(k,67)
         mat(k,820) = -rxt(k,191)*y(k,67)
         mat(k,1052) = -rxt(k,208)*y(k,67)
         mat(k,911) = -rxt(k,226)*y(k,67)
         mat(k,867) = -rxt(k,244)*y(k,67)
         mat(k,779) = rxt(k,114)*y(k,110)
         mat(k,1792) = mat(k,1792) + rxt(k,114)*y(k,65)
         mat(k,619) = -(rxt(k,122)*y(k,60) + rxt(k,123)*y(k,24) + rxt(k,124)*y(k,32) &
                      + rxt(k,125)*y(k,40) + rxt(k,126)*y(k,54) + rxt(k,127)*y(k,56) &
                      + rxt(k,146)*y(k,99) + rxt(k,169)*y(k,77) + rxt(k,186)*y(k,73) &
                      + rxt(k,204)*y(k,92) + rxt(k,220)*y(k,75) + rxt(k,222)*y(k,89) &
                      + rxt(k,239)*y(k,88))
         mat(k,1665) = -rxt(k,122)*y(k,68)
         mat(k,1510) = -rxt(k,123)*y(k,68)
         mat(k,1701) = -rxt(k,124)*y(k,68)
         mat(k,1467) = -rxt(k,125)*y(k,68)
         mat(k,970) = -rxt(k,126)*y(k,68)
         mat(k,1422) = -rxt(k,127)*y(k,68)
         mat(k,1342) = -rxt(k,146)*y(k,68)
         mat(k,1553) = -rxt(k,169)*y(k,68)
         mat(k,823) = -rxt(k,186)*y(k,68)
         mat(k,1055) = -rxt(k,204)*y(k,68)
         mat(k,1183) = -rxt(k,220)*y(k,68)
         mat(k,914) = -rxt(k,222)*y(k,68)
         mat(k,871) = -rxt(k,239)*y(k,68)
         mat(k,1384) = rxt(k,311)*y(k,98) + rxt(k,330)*y(k,102)
         mat(k,1297) = rxt(k,311)*y(k,28)
         mat(k,97) = rxt(k,330)*y(k,28)
         mat(k,1023) = -((rxt(k,128) + rxt(k,130)) * y(k,55) + rxt(k,129)*y(k,56) &
                      + rxt(k,133)*y(k,70) + rxt(k,136)*y(k,75) + rxt(k,139)*y(k,77) &
                      + rxt(k,298)*y(k,89) + rxt(k,299)*y(k,90) + rxt(k,301)*y(k,91) &
                      + rxt(k,303)*y(k,92) + rxt(k,321)*y(k,100) + rxt(k,375)*y(k,103) &
                      + rxt(k,376)*y(k,78) + rxt(k,377)*y(k,72) + rxt(k,378)*y(k,73) &
                      + rxt(k,379)*y(k,93) + rxt(k,561)*y(k,88) + rxt(k,562)*y(k,99) &
                      + rxt(k,563)*y(k,80))
         mat(k,1271) = -(rxt(k,128) + rxt(k,130)) * y(k,69)
         mat(k,1432) = -rxt(k,129)*y(k,69)
         mat(k,36) = -rxt(k,133)*y(k,69)
         mat(k,1193) = -rxt(k,136)*y(k,69)
         mat(k,1563) = -rxt(k,139)*y(k,69)
         mat(k,924) = -rxt(k,298)*y(k,69)
         mat(k,48) = -rxt(k,299)*y(k,69)
         mat(k,118) = -rxt(k,301)*y(k,69)
         mat(k,1065) = -rxt(k,303)*y(k,69)
         mat(k,145) = -rxt(k,321)*y(k,69)
         mat(k,289) = -rxt(k,375)*y(k,69)
         mat(k,136) = -rxt(k,376)*y(k,69)
         mat(k,105) = -rxt(k,377)*y(k,69)
         mat(k,833) = -rxt(k,378)*y(k,69)
         mat(k,111) = -rxt(k,379)*y(k,69)
         mat(k,881) = -rxt(k,561)*y(k,69)
         mat(k,1352) = -rxt(k,562)*y(k,69)
         mat(k,92) = -rxt(k,563)*y(k,69)
         mat(k,1711) = rxt(k,100)*y(k,61) + rxt(k,312)*y(k,98) + rxt(k,339)*y(k,105)
         mat(k,425) = rxt(k,348)*y(k,106)
         mat(k,1753) = rxt(k,131)*y(k,106)
         mat(k,980) = rxt(k,319)*y(k,98) + rxt(k,328)*y(k,101) + rxt(k,342)*y(k,105) &
                      + rxt(k,355)*y(k,106)
         mat(k,1432) = mat(k,1432) + rxt(k,327)*y(k,101)
         mat(k,603) = rxt(k,100)*y(k,32)
         mat(k,221) = rxt(k,316)*y(k,98) + rxt(k,356)*y(k,106)
         mat(k,1307) = rxt(k,312)*y(k,32) + rxt(k,319)*y(k,54) + rxt(k,316)*y(k,96)
         mat(k,406) = rxt(k,328)*y(k,54) + rxt(k,327)*y(k,56)
         mat(k,1596) = rxt(k,339)*y(k,32) + rxt(k,342)*y(k,54)
         mat(k,1632) = rxt(k,348)*y(k,33) + rxt(k,131)*y(k,51) + rxt(k,355)*y(k,54) &
                      + rxt(k,356)*y(k,96)
         mat(k,33) = -(rxt(k,133)*y(k,69) + rxt(k,134)*y(k,110))
         mat(k,1000) = -rxt(k,133)*y(k,70)
         mat(k,1776) = -rxt(k,134)*y(k,70)
         mat(k,139) = rxt(k,322)*y(k,110)
         mat(k,1776) = mat(k,1776) + rxt(k,322)*y(k,100)
         mat(k,297) = -(rxt(k,144)*y(k,99) + rxt(k,167)*y(k,77) + rxt(k,184)*y(k,73) &
                      + rxt(k,198)*y(k,75) + rxt(k,202)*y(k,92) + rxt(k,219)*y(k,89) &
                      + rxt(k,237)*y(k,88))
         mat(k,1330) = -rxt(k,144)*y(k,71)
         mat(k,1542) = -rxt(k,167)*y(k,71)
         mat(k,812) = -rxt(k,184)*y(k,71)
         mat(k,1172) = -rxt(k,198)*y(k,71)
         mat(k,1044) = -rxt(k,202)*y(k,71)
         mat(k,903) = -rxt(k,219)*y(k,71)
         mat(k,859) = -rxt(k,237)*y(k,71)
         mat(k,1381) = rxt(k,338)*y(k,105)
         mat(k,1584) = rxt(k,338)*y(k,28)
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
         mat(k,102) = -(rxt(k,369)*y(k,110) + rxt(k,377)*y(k,69))
         mat(k,1778) = -rxt(k,369)*y(k,72)
         mat(k,1005) = -rxt(k,377)*y(k,72)
         mat(k,34) = rxt(k,134)*y(k,110)
         mat(k,133) = rxt(k,367)*y(k,110)
         mat(k,1778) = mat(k,1778) + rxt(k,134)*y(k,70) + rxt(k,367)*y(k,78)
         mat(k,829) = -(rxt(k,180)*y(k,87) + rxt(k,181)*y(k,65) + rxt(k,182)*y(k,63) &
                      + rxt(k,183)*y(k,83) + rxt(k,184)*y(k,71) + rxt(k,185)*y(k,98) &
                      + rxt(k,186)*y(k,68) + rxt(k,188)*y(k,85) + rxt(k,189)*y(k,66) &
                      + rxt(k,190)*y(k,61) + rxt(k,191)*y(k,67) + rxt(k,192)*y(k,82) &
                      + rxt(k,193)*y(k,86) + rxt(k,194)*y(k,62) + rxt(k,195)*y(k,84) &
                      + rxt(k,196)*y(k,81) + rxt(k,371)*y(k,110) + rxt(k,378)*y(k,69))
         mat(k,440) = -rxt(k,180)*y(k,73)
         mat(k,786) = -rxt(k,181)*y(k,73)
         mat(k,342) = -rxt(k,182)*y(k,73)
         mat(k,654) = -rxt(k,183)*y(k,73)
         mat(k,299) = -rxt(k,184)*y(k,73)
         mat(k,1303) = -rxt(k,185)*y(k,73)
         mat(k,625) = -rxt(k,186)*y(k,73)
         mat(k,521) = -rxt(k,188)*y(k,73)
         mat(k,313) = -rxt(k,189)*y(k,73)
         mat(k,600) = -rxt(k,190)*y(k,73)
         mat(k,502) = -rxt(k,191)*y(k,73)
         mat(k,373) = -rxt(k,192)*y(k,73)
         mat(k,387) = -rxt(k,193)*y(k,73)
         mat(k,358) = -rxt(k,194)*y(k,73)
         mat(k,480) = -rxt(k,195)*y(k,73)
         mat(k,714) = -rxt(k,196)*y(k,73)
         mat(k,1802) = -rxt(k,371)*y(k,73)
         mat(k,1019) = -rxt(k,378)*y(k,73)
         mat(k,104) = rxt(k,369)*y(k,110)
         mat(k,47) = rxt(k,300)*y(k,110)
         mat(k,1802) = mat(k,1802) + rxt(k,369)*y(k,72) + rxt(k,300)*y(k,90)
         mat(k,22) = -(rxt(k,135)*y(k,110))
         mat(k,1774) = -rxt(k,135)*y(k,74)
         mat(k,227) = rxt(k,137)*y(k,75)
         mat(k,1170) = rxt(k,137)*y(k,50)
         mat(k,1197) = -(rxt(k,136)*y(k,69) + rxt(k,137)*y(k,50) + (rxt(k,141) &
                      + rxt(k,264)) * y(k,87) + rxt(k,142)*y(k,61) + (rxt(k,153) &
                      + rxt(k,255)) * y(k,67) + rxt(k,157)*y(k,82) + rxt(k,158) &
                      *y(k,86) + rxt(k,159)*y(k,62) + rxt(k,160)*y(k,84) + rxt(k,161) &
                      *y(k,81) + (rxt(k,165) + rxt(k,253)) * y(k,65) + (rxt(k,176) &
                      + rxt(k,262)) * y(k,63) + (rxt(k,187) + rxt(k,259)) * y(k,83) &
                      + rxt(k,198)*y(k,71) + rxt(k,209)*y(k,98) + rxt(k,220)*y(k,68) &
                      + (rxt(k,231) + rxt(k,257)) * y(k,85) + (rxt(k,242) + rxt(k,266) &
                      ) * y(k,66) + rxt(k,373)*y(k,110))
         mat(k,1027) = -rxt(k,136)*y(k,75)
         mat(k,236) = -rxt(k,137)*y(k,75)
         mat(k,445) = -(rxt(k,141) + rxt(k,264)) * y(k,75)
         mat(k,607) = -rxt(k,142)*y(k,75)
         mat(k,508) = -(rxt(k,153) + rxt(k,255)) * y(k,75)
         mat(k,378) = -rxt(k,157)*y(k,75)
         mat(k,392) = -rxt(k,158)*y(k,75)
         mat(k,362) = -rxt(k,159)*y(k,75)
         mat(k,487) = -rxt(k,160)*y(k,75)
         mat(k,722) = -rxt(k,161)*y(k,75)
         mat(k,794) = -(rxt(k,165) + rxt(k,253)) * y(k,75)
         mat(k,346) = -(rxt(k,176) + rxt(k,262)) * y(k,75)
         mat(k,662) = -(rxt(k,187) + rxt(k,259)) * y(k,75)
         mat(k,303) = -rxt(k,198)*y(k,75)
         mat(k,1311) = -rxt(k,209)*y(k,75)
         mat(k,633) = -rxt(k,220)*y(k,75)
         mat(k,528) = -(rxt(k,231) + rxt(k,257)) * y(k,75)
         mat(k,318) = -(rxt(k,242) + rxt(k,266)) * y(k,75)
         mat(k,1810) = -rxt(k,373)*y(k,75)
         mat(k,837) = rxt(k,371)*y(k,110)
         mat(k,24) = rxt(k,135)*y(k,110)
         mat(k,1810) = mat(k,1810) + rxt(k,371)*y(k,73) + rxt(k,135)*y(k,74)
         mat(k,26) = -(rxt(k,138)*y(k,110))
         mat(k,1775) = -rxt(k,138)*y(k,76)
         mat(k,228) = rxt(k,140)*y(k,77)
         mat(k,1540) = rxt(k,140)*y(k,50)
         mat(k,1576) = -(rxt(k,139)*y(k,69) + rxt(k,140)*y(k,50) + (rxt(k,162) &
                      + rxt(k,265)) * y(k,87) + (rxt(k,163) + rxt(k,260)) * y(k,65) &
                      + (rxt(k,164) + rxt(k,263)) * y(k,63) + (rxt(k,166) + rxt(k,261) &
                      ) * y(k,83) + rxt(k,167)*y(k,71) + rxt(k,168)*y(k,98) + rxt(k,169) &
                      *y(k,68) + (rxt(k,170) + rxt(k,258)) * y(k,85) + (rxt(k,171) &
                      + rxt(k,254)) * y(k,66) + rxt(k,172)*y(k,61) + (rxt(k,173) &
                      + rxt(k,256)) * y(k,67) + rxt(k,174)*y(k,82) + rxt(k,175) &
                      *y(k,86) + rxt(k,177)*y(k,62) + rxt(k,178)*y(k,84) + rxt(k,179) &
                      *y(k,81))
         mat(k,1036) = -rxt(k,139)*y(k,77)
         mat(k,238) = -rxt(k,140)*y(k,77)
         mat(k,449) = -(rxt(k,162) + rxt(k,265)) * y(k,77)
         mat(k,803) = -(rxt(k,163) + rxt(k,260)) * y(k,77)
         mat(k,350) = -(rxt(k,164) + rxt(k,263)) * y(k,77)
         mat(k,668) = -(rxt(k,166) + rxt(k,261)) * y(k,77)
         mat(k,307) = -rxt(k,167)*y(k,77)
         mat(k,1320) = -rxt(k,168)*y(k,77)
         mat(k,641) = -rxt(k,169)*y(k,77)
         mat(k,532) = -(rxt(k,170) + rxt(k,258)) * y(k,77)
         mat(k,322) = -(rxt(k,171) + rxt(k,254)) * y(k,77)
         mat(k,611) = -rxt(k,172)*y(k,77)
         mat(k,512) = -(rxt(k,173) + rxt(k,256)) * y(k,77)
         mat(k,381) = -rxt(k,174)*y(k,77)
         mat(k,397) = -rxt(k,175)*y(k,77)
         mat(k,366) = -rxt(k,177)*y(k,77)
         mat(k,491) = -rxt(k,178)*y(k,77)
         mat(k,729) = -rxt(k,179)*y(k,77)
         mat(k,1206) = rxt(k,373)*y(k,110)
         mat(k,28) = rxt(k,138)*y(k,110)
         mat(k,1819) = rxt(k,373)*y(k,75) + rxt(k,138)*y(k,76)
         mat(k,134) = -(rxt(k,367)*y(k,110) + rxt(k,376)*y(k,69))
         mat(k,1782) = -rxt(k,367)*y(k,78)
         mat(k,1009) = -rxt(k,376)*y(k,78)
         mat(k,1695) = rxt(k,304)*y(k,92)
         mat(k,676) = rxt(k,305)*y(k,92)
         mat(k,1043) = rxt(k,304)*y(k,32) + rxt(k,305)*y(k,43) + rxt(k,306)*y(k,104)
         mat(k,141) = rxt(k,323)*y(k,110)
         mat(k,742) = rxt(k,306)*y(k,92)
         mat(k,1782) = mat(k,1782) + rxt(k,323)*y(k,100)
         mat(k,159) = -(rxt(k,419)*y(k,54) + rxt(k,420)*y(k,55) + rxt(k,575)*y(k,107))
         mat(k,953) = -rxt(k,419)*y(k,79)
         mat(k,1253) = -rxt(k,420)*y(k,79)
         mat(k,242) = -rxt(k,575)*y(k,79)
         mat(k,953) = mat(k,953) + rxt(k,565)*y(k,80)
         mat(k,1011) = .900_r8*rxt(k,563)*y(k,80) + .800_r8*rxt(k,561)*y(k,88)
         mat(k,87) = rxt(k,565)*y(k,54) + .900_r8*rxt(k,563)*y(k,69)
         mat(k,855) = .800_r8*rxt(k,561)*y(k,69)
         mat(k,86) = -(rxt(k,563)*y(k,69) + rxt(k,564)*y(k,55) + (rxt(k,565) + rxt(k,566) &
                      ) * y(k,54))
         mat(k,1004) = -rxt(k,563)*y(k,80)
         mat(k,1251) = -rxt(k,564)*y(k,80)
         mat(k,949) = -(rxt(k,565) + rxt(k,566)) * y(k,80)
         mat(k,712) = -(rxt(k,156)*y(k,99) + rxt(k,161)*y(k,75) + rxt(k,179)*y(k,77) &
                      + rxt(k,196)*y(k,73) + rxt(k,214)*y(k,92) + rxt(k,232)*y(k,89) &
                      + rxt(k,249)*y(k,88) + rxt(k,280)*y(k,60) + rxt(k,281)*y(k,24) &
                      + rxt(k,282)*y(k,32) + rxt(k,283)*y(k,110) + rxt(k,284)*y(k,40) &
                      + rxt(k,285)*y(k,42) + rxt(k,286)*y(k,52) + rxt(k,287)*y(k,56))
         mat(k,1345) = -rxt(k,156)*y(k,81)
         mat(k,1186) = -rxt(k,161)*y(k,81)
         mat(k,1556) = -rxt(k,179)*y(k,81)
         mat(k,826) = -rxt(k,196)*y(k,81)
         mat(k,1058) = -rxt(k,214)*y(k,81)
         mat(k,917) = -rxt(k,232)*y(k,81)
         mat(k,874) = -rxt(k,249)*y(k,81)
         mat(k,1668) = -rxt(k,280)*y(k,81)
         mat(k,1513) = -rxt(k,281)*y(k,81)
         mat(k,1704) = -rxt(k,282)*y(k,81)
         mat(k,1799) = -rxt(k,283)*y(k,81)
         mat(k,1470) = -rxt(k,284)*y(k,81)
         mat(k,1098) = -rxt(k,285)*y(k,81)
         mat(k,1143) = -rxt(k,286)*y(k,81)
         mat(k,1425) = -rxt(k,287)*y(k,81)
         mat(k,1746) = rxt(k,106)*y(k,64) + rxt(k,275)*y(k,65) + rxt(k,118)*y(k,67) &
                      + rxt(k,274)*y(k,101)
         mat(k,1143) = mat(k,1143) + rxt(k,105)*y(k,61) + rxt(k,315)*y(k,98) &
                      + rxt(k,273)*y(k,101) + rxt(k,341)*y(k,105) + rxt(k,354) &
                      *y(k,106)
         mat(k,973) = rxt(k,296)*y(k,83)
         mat(k,1425) = mat(k,1425) + rxt(k,297)*y(k,83)
         mat(k,599) = rxt(k,105)*y(k,52)
         mat(k,213) = rxt(k,106)*y(k,51)
         mat(k,783) = rxt(k,275)*y(k,51)
         mat(k,500) = rxt(k,118)*y(k,51)
         mat(k,653) = rxt(k,296)*y(k,54) + rxt(k,297)*y(k,56)
         mat(k,1300) = rxt(k,315)*y(k,52)
         mat(k,403) = rxt(k,274)*y(k,51) + rxt(k,273)*y(k,52)
         mat(k,1589) = rxt(k,341)*y(k,52)
         mat(k,1625) = rxt(k,354)*y(k,52)
         mat(k,371) = -(rxt(k,151)*y(k,99) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77) &
                      + rxt(k,192)*y(k,73) + rxt(k,210)*y(k,92) + rxt(k,227)*y(k,89) &
                      + rxt(k,245)*y(k,88))
         mat(k,1334) = -rxt(k,151)*y(k,82)
         mat(k,1176) = -rxt(k,157)*y(k,82)
         mat(k,1546) = -rxt(k,174)*y(k,82)
         mat(k,816) = -rxt(k,192)*y(k,82)
         mat(k,1048) = -rxt(k,210)*y(k,82)
         mat(k,907) = -rxt(k,227)*y(k,82)
         mat(k,863) = -rxt(k,245)*y(k,82)
         mat(k,1736) = rxt(k,117)*y(k,67)
         mat(k,496) = rxt(k,117)*y(k,51)
         mat(k,709) = rxt(k,283)*y(k,110)
         mat(k,1789) = rxt(k,283)*y(k,81)
         mat(k,652) = -(rxt(k,143)*y(k,99) + (rxt(k,166) + rxt(k,261)) * y(k,77) &
                      + rxt(k,183)*y(k,73) + (rxt(k,187) + rxt(k,259)) * y(k,75) &
                      + rxt(k,201)*y(k,92) + rxt(k,218)*y(k,89) + rxt(k,236)*y(k,88) &
                      + (rxt(k,271) + rxt(k,293)) * y(k,40) + rxt(k,291)*y(k,110) &
                      + rxt(k,295)*y(k,42) + rxt(k,296)*y(k,54) + rxt(k,297)*y(k,56))
         mat(k,1343) = -rxt(k,143)*y(k,83)
         mat(k,1554) = -(rxt(k,166) + rxt(k,261)) * y(k,83)
         mat(k,824) = -rxt(k,183)*y(k,83)
         mat(k,1184) = -(rxt(k,187) + rxt(k,259)) * y(k,83)
         mat(k,1056) = -rxt(k,201)*y(k,83)
         mat(k,915) = -rxt(k,218)*y(k,83)
         mat(k,872) = -rxt(k,236)*y(k,83)
         mat(k,1468) = -(rxt(k,271) + rxt(k,293)) * y(k,83)
         mat(k,1797) = -rxt(k,291)*y(k,83)
         mat(k,1096) = -rxt(k,295)*y(k,83)
         mat(k,971) = -rxt(k,296)*y(k,83)
         mat(k,1423) = -rxt(k,297)*y(k,83)
         mat(k,1096) = mat(k,1096) + rxt(k,104)*y(k,61) + rxt(k,119)*y(k,65) &
                      + rxt(k,285)*y(k,81) + rxt(k,314)*y(k,98) + rxt(k,352)*y(k,106)
         mat(k,1744) = rxt(k,267)*y(k,101)
         mat(k,1141) = rxt(k,276)*y(k,65) + rxt(k,115)*y(k,67) + rxt(k,286)*y(k,81) &
                      + rxt(k,272)*y(k,101)
         mat(k,1423) = mat(k,1423) + rxt(k,287)*y(k,81)
         mat(k,598) = rxt(k,104)*y(k,42)
         mat(k,782) = rxt(k,119)*y(k,42) + rxt(k,276)*y(k,52)
         mat(k,499) = rxt(k,115)*y(k,52)
         mat(k,711) = rxt(k,285)*y(k,42) + rxt(k,286)*y(k,52) + rxt(k,287)*y(k,56)
         mat(k,1298) = rxt(k,314)*y(k,42)
         mat(k,402) = rxt(k,267)*y(k,51) + rxt(k,272)*y(k,52)
         mat(k,1623) = rxt(k,352)*y(k,42)
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
         mat(k,477) = -(rxt(k,155)*y(k,99) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77) &
                      + rxt(k,195)*y(k,73) + rxt(k,213)*y(k,92) + rxt(k,230)*y(k,89) &
                      + rxt(k,248)*y(k,88) + rxt(k,288)*y(k,50))
         mat(k,1337) = -rxt(k,155)*y(k,84)
         mat(k,1179) = -rxt(k,160)*y(k,84)
         mat(k,1549) = -rxt(k,178)*y(k,84)
         mat(k,819) = -rxt(k,195)*y(k,84)
         mat(k,1051) = -rxt(k,213)*y(k,84)
         mat(k,910) = -rxt(k,230)*y(k,84)
         mat(k,866) = -rxt(k,248)*y(k,84)
         mat(k,231) = -rxt(k,288)*y(k,84)
         mat(k,518) = rxt(k,289)*y(k,110)
         mat(k,1791) = rxt(k,289)*y(k,85)
         mat(k,519) = -(rxt(k,147)*y(k,99) + (rxt(k,170) + rxt(k,258)) * y(k,77) &
                      + rxt(k,188)*y(k,73) + rxt(k,205)*y(k,92) + rxt(k,223)*y(k,89) &
                      + (rxt(k,231) + rxt(k,257)) * y(k,75) + rxt(k,240)*y(k,88) &
                      + rxt(k,289)*y(k,110) + rxt(k,290)*y(k,42) + rxt(k,292)*y(k,50))
         mat(k,1339) = -rxt(k,147)*y(k,85)
         mat(k,1551) = -(rxt(k,170) + rxt(k,258)) * y(k,85)
         mat(k,821) = -rxt(k,188)*y(k,85)
         mat(k,1053) = -rxt(k,205)*y(k,85)
         mat(k,912) = -rxt(k,223)*y(k,85)
         mat(k,1181) = -(rxt(k,231) + rxt(k,257)) * y(k,85)
         mat(k,868) = -rxt(k,240)*y(k,85)
         mat(k,1793) = -rxt(k,289)*y(k,85)
         mat(k,1092) = -rxt(k,290)*y(k,85)
         mat(k,232) = -rxt(k,292)*y(k,85)
         mat(k,1137) = rxt(k,116)*y(k,67)
         mat(k,498) = rxt(k,116)*y(k,52)
         mat(k,650) = rxt(k,291)*y(k,110)
         mat(k,1793) = mat(k,1793) + rxt(k,291)*y(k,83)
         mat(k,385) = -(rxt(k,152)*y(k,99) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77) &
                      + rxt(k,193)*y(k,73) + rxt(k,211)*y(k,92) + rxt(k,228)*y(k,89) &
                      + rxt(k,246)*y(k,88) + rxt(k,294)*y(k,42))
         mat(k,1335) = -rxt(k,152)*y(k,86)
         mat(k,1177) = -rxt(k,158)*y(k,86)
         mat(k,1547) = -rxt(k,175)*y(k,86)
         mat(k,817) = -rxt(k,193)*y(k,86)
         mat(k,1049) = -rxt(k,211)*y(k,86)
         mat(k,908) = -rxt(k,228)*y(k,86)
         mat(k,864) = -rxt(k,246)*y(k,86)
         mat(k,1089) = -rxt(k,294)*y(k,86)
         mat(k,1461) = rxt(k,271)*y(k,83)
         mat(k,648) = rxt(k,271)*y(k,40)
         mat(k,438) = -((rxt(k,141) + rxt(k,264)) * y(k,75) + (rxt(k,162) + rxt(k,265) &
                      ) * y(k,77) + rxt(k,180)*y(k,73) + rxt(k,197)*y(k,92) + rxt(k,215) &
                      *y(k,89) + rxt(k,233)*y(k,88) + rxt(k,250)*y(k,99))
         mat(k,1178) = -(rxt(k,141) + rxt(k,264)) * y(k,87)
         mat(k,1548) = -(rxt(k,162) + rxt(k,265)) * y(k,87)
         mat(k,818) = -rxt(k,180)*y(k,87)
         mat(k,1050) = -rxt(k,197)*y(k,87)
         mat(k,909) = -rxt(k,215)*y(k,87)
         mat(k,865) = -rxt(k,233)*y(k,87)
         mat(k,1336) = -rxt(k,250)*y(k,87)
         mat(k,1091) = rxt(k,295)*y(k,83) + rxt(k,290)*y(k,85) + rxt(k,294)*y(k,86)
         mat(k,230) = rxt(k,288)*y(k,84) + rxt(k,292)*y(k,85)
         mat(k,649) = rxt(k,295)*y(k,42)
         mat(k,476) = rxt(k,288)*y(k,50)
         mat(k,517) = rxt(k,290)*y(k,42) + rxt(k,292)*y(k,50)
         mat(k,386) = rxt(k,294)*y(k,42)
         mat(k,878) = -(rxt(k,233)*y(k,87) + rxt(k,234)*y(k,65) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,83) + rxt(k,237)*y(k,71) + rxt(k,238)*y(k,98) &
                      + rxt(k,239)*y(k,68) + rxt(k,240)*y(k,85) + rxt(k,241)*y(k,66) &
                      + rxt(k,243)*y(k,61) + rxt(k,244)*y(k,67) + rxt(k,245)*y(k,82) &
                      + rxt(k,246)*y(k,86) + rxt(k,247)*y(k,62) + rxt(k,248)*y(k,84) &
                      + rxt(k,249)*y(k,81) + rxt(k,360)*y(k,110) + rxt(k,363)*y(k,28) &
                      + rxt(k,561)*y(k,69))
         mat(k,441) = -rxt(k,233)*y(k,88)
         mat(k,787) = -rxt(k,234)*y(k,88)
         mat(k,343) = -rxt(k,235)*y(k,88)
         mat(k,655) = -rxt(k,236)*y(k,88)
         mat(k,300) = -rxt(k,237)*y(k,88)
         mat(k,1304) = -rxt(k,238)*y(k,88)
         mat(k,626) = -rxt(k,239)*y(k,88)
         mat(k,522) = -rxt(k,240)*y(k,88)
         mat(k,314) = -rxt(k,241)*y(k,88)
         mat(k,601) = -rxt(k,243)*y(k,88)
         mat(k,503) = -rxt(k,244)*y(k,88)
         mat(k,374) = -rxt(k,245)*y(k,88)
         mat(k,388) = -rxt(k,246)*y(k,88)
         mat(k,359) = -rxt(k,247)*y(k,88)
         mat(k,481) = -rxt(k,248)*y(k,88)
         mat(k,715) = -rxt(k,249)*y(k,88)
         mat(k,1803) = -rxt(k,360)*y(k,88)
         mat(k,1391) = -rxt(k,363)*y(k,88)
         mat(k,1020) = -rxt(k,561)*y(k,88)
         mat(k,274) = rxt(k,570)*y(k,99)
         mat(k,1750) = rxt(k,572)*y(k,99)
         mat(k,977) = rxt(k,565)*y(k,80)
         mat(k,1268) = rxt(k,569)*y(k,94)
         mat(k,90) = rxt(k,565)*y(k,54)
         mat(k,179) = rxt(k,569)*y(k,55)
         mat(k,1349) = rxt(k,570)*y(k,48) + rxt(k,572)*y(k,51)
         mat(k,922) = -(rxt(k,215)*y(k,87) + rxt(k,216)*y(k,65) + rxt(k,217)*y(k,63) &
                      + rxt(k,218)*y(k,83) + rxt(k,219)*y(k,71) + rxt(k,221)*y(k,98) &
                      + rxt(k,222)*y(k,68) + rxt(k,223)*y(k,85) + rxt(k,224)*y(k,66) &
                      + rxt(k,225)*y(k,61) + rxt(k,226)*y(k,67) + rxt(k,227)*y(k,82) &
                      + rxt(k,228)*y(k,86) + rxt(k,229)*y(k,62) + rxt(k,230)*y(k,84) &
                      + rxt(k,232)*y(k,81) + rxt(k,298)*y(k,69) + rxt(k,362)*y(k,110))
         mat(k,442) = -rxt(k,215)*y(k,89)
         mat(k,788) = -rxt(k,216)*y(k,89)
         mat(k,344) = -rxt(k,217)*y(k,89)
         mat(k,656) = -rxt(k,218)*y(k,89)
         mat(k,301) = -rxt(k,219)*y(k,89)
         mat(k,1305) = -rxt(k,221)*y(k,89)
         mat(k,627) = -rxt(k,222)*y(k,89)
         mat(k,523) = -rxt(k,223)*y(k,89)
         mat(k,315) = -rxt(k,224)*y(k,89)
         mat(k,602) = -rxt(k,225)*y(k,89)
         mat(k,504) = -rxt(k,226)*y(k,89)
         mat(k,375) = -rxt(k,227)*y(k,89)
         mat(k,389) = -rxt(k,228)*y(k,89)
         mat(k,360) = -rxt(k,229)*y(k,89)
         mat(k,482) = -rxt(k,230)*y(k,89)
         mat(k,716) = -rxt(k,232)*y(k,89)
         mat(k,1021) = -rxt(k,298)*y(k,89)
         mat(k,1804) = -rxt(k,362)*y(k,89)
         mat(k,1063) = rxt(k,361)*y(k,110)
         mat(k,1804) = mat(k,1804) + rxt(k,361)*y(k,92)
         mat(k,45) = -(rxt(k,299)*y(k,69) + rxt(k,300)*y(k,110))
         mat(k,1001) = -rxt(k,299)*y(k,90)
         mat(k,1777) = -rxt(k,300)*y(k,90)
         mat(k,901) = rxt(k,362)*y(k,110)
         mat(k,1777) = mat(k,1777) + rxt(k,362)*y(k,89)
         mat(k,116) = -(rxt(k,301)*y(k,69) + rxt(k,302)*y(k,110))
         mat(k,1007) = -rxt(k,301)*y(k,91)
         mat(k,1780) = -rxt(k,302)*y(k,91)
         mat(k,1375) = rxt(k,363)*y(k,88) + rxt(k,307)*y(k,93)
         mat(k,854) = rxt(k,363)*y(k,28)
         mat(k,109) = rxt(k,307)*y(k,28)
         mat(k,1066) = -(rxt(k,197)*y(k,87) + rxt(k,199)*y(k,65) + rxt(k,200)*y(k,63) &
                      + rxt(k,201)*y(k,83) + rxt(k,202)*y(k,71) + rxt(k,203)*y(k,98) &
                      + rxt(k,204)*y(k,68) + rxt(k,205)*y(k,85) + rxt(k,206)*y(k,66) &
                      + rxt(k,207)*y(k,61) + rxt(k,208)*y(k,67) + rxt(k,210)*y(k,82) &
                      + rxt(k,211)*y(k,86) + rxt(k,212)*y(k,62) + rxt(k,213)*y(k,84) &
                      + rxt(k,214)*y(k,81) + rxt(k,303)*y(k,69) + rxt(k,304)*y(k,32) &
                      + rxt(k,305)*y(k,43) + rxt(k,306)*y(k,104) + rxt(k,361)*y(k,110))
         mat(k,443) = -rxt(k,197)*y(k,92)
         mat(k,791) = -rxt(k,199)*y(k,92)
         mat(k,345) = -rxt(k,200)*y(k,92)
         mat(k,659) = -rxt(k,201)*y(k,92)
         mat(k,302) = -rxt(k,202)*y(k,92)
         mat(k,1308) = -rxt(k,203)*y(k,92)
         mat(k,630) = -rxt(k,204)*y(k,92)
         mat(k,525) = -rxt(k,205)*y(k,92)
         mat(k,317) = -rxt(k,206)*y(k,92)
         mat(k,604) = -rxt(k,207)*y(k,92)
         mat(k,506) = -rxt(k,208)*y(k,92)
         mat(k,376) = -rxt(k,210)*y(k,92)
         mat(k,390) = -rxt(k,211)*y(k,92)
         mat(k,361) = -rxt(k,212)*y(k,92)
         mat(k,484) = -rxt(k,213)*y(k,92)
         mat(k,719) = -rxt(k,214)*y(k,92)
         mat(k,1024) = -rxt(k,303)*y(k,92)
         mat(k,1712) = -rxt(k,304)*y(k,92)
         mat(k,691) = -rxt(k,305)*y(k,92)
         mat(k,760) = -rxt(k,306)*y(k,92)
         mat(k,1807) = -rxt(k,361)*y(k,92)
         mat(k,882) = rxt(k,360)*y(k,110)
         mat(k,119) = rxt(k,302)*y(k,110)
         mat(k,112) = rxt(k,308)*y(k,110)
         mat(k,1807) = mat(k,1807) + rxt(k,360)*y(k,88) + rxt(k,302)*y(k,91) &
                      + rxt(k,308)*y(k,93)
         mat(k,108) = -(rxt(k,307)*y(k,28) + rxt(k,308)*y(k,110) + rxt(k,379)*y(k,69))
         mat(k,1374) = -rxt(k,307)*y(k,93)
         mat(k,1779) = -rxt(k,308)*y(k,93)
         mat(k,1006) = -rxt(k,379)*y(k,93)
         mat(k,176) = -(rxt(k,567)*y(k,54) + (rxt(k,568) + rxt(k,569)) * y(k,55))
         mat(k,954) = -rxt(k,567)*y(k,94)
         mat(k,1254) = -(rxt(k,568) + rxt(k,569)) * y(k,94)
         mat(k,160) = rxt(k,575)*y(k,107)
         mat(k,243) = rxt(k,575)*y(k,79)
         mat(k,573) = -(rxt(k,384)*y(k,33) + rxt(k,385)*y(k,110) + (rxt(k,387) &
                      + rxt(k,388)) * y(k,55) + rxt(k,389)*y(k,56) + (rxt(k,477) &
                      + rxt(k,478)) * y(k,40) + (rxt(k,500) + rxt(k,501)) * y(k,36) &
                      + rxt(k,506)*y(k,29) + rxt(k,507)*y(k,30))
         mat(k,421) = -rxt(k,384)*y(k,95)
         mat(k,1795) = -rxt(k,385)*y(k,95)
         mat(k,1260) = -(rxt(k,387) + rxt(k,388)) * y(k,95)
         mat(k,1421) = -rxt(k,389)*y(k,95)
         mat(k,1465) = -(rxt(k,477) + rxt(k,478)) * y(k,95)
         mat(k,194) = -(rxt(k,500) + rxt(k,501)) * y(k,95)
         mat(k,6) = -rxt(k,506)*y(k,95)
         mat(k,14) = -rxt(k,507)*y(k,95)
         mat(k,1260) = mat(k,1260) + rxt(k,420)*y(k,79)
         mat(k,1016) = .850_r8*rxt(k,562)*y(k,99)
         mat(k,163) = rxt(k,420)*y(k,55)
         mat(k,1340) = .850_r8*rxt(k,562)*y(k,69)
         mat(k,219) = -(rxt(k,316)*y(k,98) + rxt(k,334)*y(k,103) + rxt(k,356)*y(k,106) &
                      + rxt(k,391)*y(k,54) + rxt(k,392)*y(k,55))
         mat(k,1293) = -rxt(k,316)*y(k,96)
         mat(k,284) = -rxt(k,334)*y(k,96)
         mat(k,1616) = -rxt(k,356)*y(k,96)
         mat(k,959) = -rxt(k,391)*y(k,96)
         mat(k,1255) = -rxt(k,392)*y(k,96)
         mat(k,1377) = rxt(k,393)*y(k,97)
         mat(k,959) = mat(k,959) + rxt(k,395)*y(k,97)
         mat(k,1255) = mat(k,1255) + rxt(k,396)*y(k,97)
         mat(k,1415) = rxt(k,397)*y(k,97)
         mat(k,31) = rxt(k,393)*y(k,28) + rxt(k,395)*y(k,54) + rxt(k,396)*y(k,55) &
                      + rxt(k,397)*y(k,56)
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
         mat(k,30) = -(rxt(k,393)*y(k,28) + rxt(k,395)*y(k,54) + rxt(k,396)*y(k,55) &
                      + rxt(k,397)*y(k,56))
         mat(k,1372) = -rxt(k,393)*y(k,97)
         mat(k,944) = -rxt(k,395)*y(k,97)
         mat(k,1248) = -rxt(k,396)*y(k,97)
         mat(k,1414) = -rxt(k,397)*y(k,97)
         mat(k,1248) = mat(k,1248) + rxt(k,387)*y(k,95)
         mat(k,563) = rxt(k,387)*y(k,55)
         mat(k,1314) = -(rxt(k,145)*y(k,99) + rxt(k,168)*y(k,77) + rxt(k,185)*y(k,73) &
                      + rxt(k,203)*y(k,92) + rxt(k,209)*y(k,75) + rxt(k,221)*y(k,89) &
                      + rxt(k,238)*y(k,88) + rxt(k,309)*y(k,60) + rxt(k,310)*y(k,24) &
                      + rxt(k,311)*y(k,28) + rxt(k,312)*y(k,32) + rxt(k,313)*y(k,40) &
                      + rxt(k,314)*y(k,42) + rxt(k,315)*y(k,52) + rxt(k,316)*y(k,96) &
                      + rxt(k,317)*y(k,55) + rxt(k,318)*y(k,56) + (rxt(k,319) &
                      + rxt(k,320)) * y(k,54))
         mat(k,1359) = -rxt(k,145)*y(k,98)
         mat(k,1570) = -rxt(k,168)*y(k,98)
         mat(k,840) = -rxt(k,185)*y(k,98)
         mat(k,1072) = -rxt(k,203)*y(k,98)
         mat(k,1200) = -rxt(k,209)*y(k,98)
         mat(k,931) = -rxt(k,221)*y(k,98)
         mat(k,888) = -rxt(k,238)*y(k,98)
         mat(k,1682) = -rxt(k,309)*y(k,98)
         mat(k,1527) = -rxt(k,310)*y(k,98)
         mat(k,1401) = -rxt(k,311)*y(k,98)
         mat(k,1718) = -rxt(k,312)*y(k,98)
         mat(k,1484) = -rxt(k,313)*y(k,98)
         mat(k,1112) = -rxt(k,314)*y(k,98)
         mat(k,1157) = -rxt(k,315)*y(k,98)
         mat(k,223) = -rxt(k,316)*y(k,98)
         mat(k,1278) = -rxt(k,317)*y(k,98)
         mat(k,1439) = -rxt(k,318)*y(k,98)
         mat(k,987) = -(rxt(k,319) + rxt(k,320)) * y(k,98)
         mat(k,987) = mat(k,987) + rxt(k,120)*y(k,65) + rxt(k,329)*y(k,101)
         mat(k,1278) = mat(k,1278) + (rxt(k,128)+rxt(k,130))*y(k,69)
         mat(k,797) = rxt(k,120)*y(k,54)
         mat(k,1030) = (rxt(k,128)+rxt(k,130))*y(k,55)
         mat(k,409) = rxt(k,329)*y(k,54)
         mat(k,1360) = -(rxt(k,143)*y(k,83) + rxt(k,144)*y(k,71) + rxt(k,145)*y(k,98) &
                      + rxt(k,146)*y(k,68) + rxt(k,147)*y(k,85) + rxt(k,148)*y(k,66) &
                      + rxt(k,149)*y(k,61) + rxt(k,150)*y(k,67) + rxt(k,151)*y(k,82) &
                      + rxt(k,152)*y(k,86) + rxt(k,154)*y(k,62) + rxt(k,155)*y(k,84) &
                      + rxt(k,156)*y(k,81) + rxt(k,250)*y(k,87) + rxt(k,251)*y(k,65) &
                      + rxt(k,252)*y(k,63) + rxt(k,324)*y(k,110) + rxt(k,359)*y(k,55) &
                      + rxt(k,562)*y(k,69) + rxt(k,570)*y(k,48) + rxt(k,572)*y(k,51))
         mat(k,665) = -rxt(k,143)*y(k,99)
         mat(k,305) = -rxt(k,144)*y(k,99)
         mat(k,1315) = -rxt(k,145)*y(k,99)
         mat(k,636) = -rxt(k,146)*y(k,99)
         mat(k,531) = -rxt(k,147)*y(k,99)
         mat(k,320) = -rxt(k,148)*y(k,99)
         mat(k,609) = -rxt(k,149)*y(k,99)
         mat(k,510) = -rxt(k,150)*y(k,99)
         mat(k,380) = -rxt(k,151)*y(k,99)
         mat(k,395) = -rxt(k,152)*y(k,99)
         mat(k,364) = -rxt(k,154)*y(k,99)
         mat(k,490) = -rxt(k,155)*y(k,99)
         mat(k,725) = -rxt(k,156)*y(k,99)
         mat(k,448) = -rxt(k,250)*y(k,99)
         mat(k,798) = -rxt(k,251)*y(k,99)
         mat(k,348) = -rxt(k,252)*y(k,99)
         mat(k,1814) = -rxt(k,324)*y(k,99)
         mat(k,1279) = -rxt(k,359)*y(k,99)
         mat(k,1031) = -rxt(k,562)*y(k,99)
         mat(k,279) = -rxt(k,570)*y(k,99)
         mat(k,1761) = -rxt(k,572)*y(k,99)
         mat(k,1402) = rxt(k,573)*y(k,107)
         mat(k,988) = rxt(k,333)*y(k,103)
         mat(k,1279) = mat(k,1279) + rxt(k,564)*y(k,80) + rxt(k,568)*y(k,94) &
                      + rxt(k,576)*y(k,107) + rxt(k,580)*y(k,108)
         mat(k,94) = rxt(k,564)*y(k,55)
         mat(k,182) = rxt(k,568)*y(k,55)
         mat(k,224) = rxt(k,334)*y(k,103)
         mat(k,1315) = mat(k,1315) + 2.000_r8*rxt(k,145)*y(k,99)
         mat(k,1360) = mat(k,1360) + 2.000_r8*rxt(k,145)*y(k,98)
         mat(k,292) = rxt(k,333)*y(k,54) + rxt(k,334)*y(k,96)
         mat(k,252) = rxt(k,573)*y(k,28) + rxt(k,576)*y(k,55)
         mat(k,70) = rxt(k,580)*y(k,55)
         mat(k,142) = -(rxt(k,321)*y(k,69) + (rxt(k,322) + rxt(k,323)) * y(k,110))
         mat(k,1010) = -rxt(k,321)*y(k,100)
         mat(k,1783) = -(rxt(k,322) + rxt(k,323)) * y(k,100)
         mat(k,1327) = rxt(k,324)*y(k,110)
         mat(k,283) = rxt(k,332)*y(k,110)
         mat(k,1783) = mat(k,1783) + rxt(k,324)*y(k,99) + rxt(k,332)*y(k,103)
         mat(k,401) = -((rxt(k,267) + rxt(k,274)) * y(k,51) + (rxt(k,272) + rxt(k,273) &
                      ) * y(k,52) + rxt(k,325)*y(k,28) + rxt(k,326)*y(k,32) + rxt(k,327) &
                      *y(k,56) + (rxt(k,328) + rxt(k,329)) * y(k,54))
         mat(k,1737) = -(rxt(k,267) + rxt(k,274)) * y(k,101)
         mat(k,1132) = -(rxt(k,272) + rxt(k,273)) * y(k,101)
         mat(k,1382) = -rxt(k,325)*y(k,101)
         mat(k,1696) = -rxt(k,326)*y(k,101)
         mat(k,1418) = -rxt(k,327)*y(k,101)
         mat(k,964) = -(rxt(k,328) + rxt(k,329)) * y(k,101)
         mat(k,964) = mat(k,964) + rxt(k,331)*y(k,102)
         mat(k,1259) = rxt(k,121)*y(k,65) + rxt(k,357)*y(k,106)
         mat(k,1418) = mat(k,1418) + rxt(k,127)*y(k,68) + rxt(k,318)*y(k,98) &
                      + rxt(k,343)*y(k,105) + rxt(k,358)*y(k,106)
         mat(k,778) = rxt(k,121)*y(k,55)
         mat(k,617) = rxt(k,127)*y(k,56)
         mat(k,1295) = rxt(k,318)*y(k,56)
         mat(k,96) = rxt(k,331)*y(k,54)
         mat(k,1585) = rxt(k,343)*y(k,56)
         mat(k,1618) = rxt(k,357)*y(k,55) + rxt(k,358)*y(k,56)
         mat(k,95) = -(rxt(k,330)*y(k,28) + rxt(k,331)*y(k,54))
         mat(k,1373) = -rxt(k,330)*y(k,102)
         mat(k,950) = -rxt(k,331)*y(k,102)
         mat(k,1252) = rxt(k,317)*y(k,98)
         mat(k,1291) = rxt(k,317)*y(k,55)
         mat(k,285) = -(rxt(k,332)*y(k,110) + rxt(k,333)*y(k,54) + rxt(k,334)*y(k,96) &
                      + rxt(k,375)*y(k,69))
         mat(k,1786) = -rxt(k,332)*y(k,103)
         mat(k,962) = -rxt(k,333)*y(k,103)
         mat(k,220) = -rxt(k,334)*y(k,103)
         mat(k,1015) = -rxt(k,375)*y(k,103)
         mat(k,1258) = rxt(k,359)*y(k,99)
         mat(k,1329) = rxt(k,359)*y(k,55)
         mat(k,755) = -(rxt(k,306)*y(k,92) + rxt(k,335)*y(k,47) + rxt(k,344)*y(k,51) &
                      + rxt(k,410)*y(k,33) + rxt(k,411)*y(k,35) + rxt(k,412)*y(k,43) &
                      + rxt(k,413)*y(k,54) + rxt(k,414)*y(k,56) + (4._r8*rxt(k,415) &
                      + 4._r8*rxt(k,416)) * y(k,104) + rxt(k,418)*y(k,44) + rxt(k,432) &
                      *y(k,53) + rxt(k,433)*y(k,48) + rxt(k,441)*y(k,52) + rxt(k,442) &
                      *y(k,42) + rxt(k,461)*y(k,25) + (rxt(k,463) + rxt(k,464) &
                      ) * y(k,24) + rxt(k,466)*y(k,40) + rxt(k,469)*y(k,46) + rxt(k,493) &
                      *y(k,2) + rxt(k,495)*y(k,36) + rxt(k,527)*y(k,14) + rxt(k,530) &
                      *y(k,19) + (rxt(k,532) + rxt(k,536)) * y(k,27))
         mat(k,1059) = -rxt(k,306)*y(k,104)
         mat(k,124) = -rxt(k,335)*y(k,104)
         mat(k,1747) = -rxt(k,344)*y(k,104)
         mat(k,423) = -rxt(k,410)*y(k,104)
         mat(k,81) = -rxt(k,411)*y(k,104)
         mat(k,687) = -rxt(k,412)*y(k,104)
         mat(k,974) = -rxt(k,413)*y(k,104)
         mat(k,1426) = -rxt(k,414)*y(k,104)
         mat(k,53) = -rxt(k,418)*y(k,104)
         mat(k,1223) = -rxt(k,432)*y(k,104)
         mat(k,273) = -rxt(k,433)*y(k,104)
         mat(k,1144) = -rxt(k,441)*y(k,104)
         mat(k,1099) = -rxt(k,442)*y(k,104)
         mat(k,202) = -rxt(k,461)*y(k,104)
         mat(k,1514) = -(rxt(k,463) + rxt(k,464)) * y(k,104)
         mat(k,1471) = -rxt(k,466)*y(k,104)
         mat(k,185) = -rxt(k,469)*y(k,104)
         mat(k,461) = -rxt(k,493)*y(k,104)
         mat(k,195) = -rxt(k,495)*y(k,104)
         mat(k,544) = -rxt(k,527)*y(k,104)
         mat(k,42) = -rxt(k,530)*y(k,104)
         mat(k,130) = -(rxt(k,532) + rxt(k,536)) * y(k,104)
         mat(k,544) = mat(k,544) + rxt(k,526)*y(k,54)
         mat(k,42) = mat(k,42) + .300_r8*rxt(k,530)*y(k,104)
         mat(k,1514) = mat(k,1514) + rxt(k,337)*y(k,105)
         mat(k,171) = rxt(k,504)*y(k,110)
         mat(k,1705) = 2.000_r8*rxt(k,404)*y(k,43) + rxt(k,409)*y(k,56) + rxt(k,124) &
                      *y(k,68)
         mat(k,423) = mat(k,423) + rxt(k,401)*y(k,54) + rxt(k,384)*y(k,95)
         mat(k,81) = mat(k,81) + rxt(k,402)*y(k,54)
         mat(k,195) = mat(k,195) + rxt(k,494)*y(k,54) + rxt(k,500)*y(k,95)
         mat(k,1471) = mat(k,1471) + rxt(k,465)*y(k,54) + rxt(k,477)*y(k,95) &
                      + rxt(k,351)*y(k,106)
         mat(k,1099) = mat(k,1099) + rxt(k,119)*y(k,65) + rxt(k,352)*y(k,106)
         mat(k,687) = mat(k,687) + 2.000_r8*rxt(k,404)*y(k,32) + rxt(k,434)*y(k,51) &
                      + rxt(k,429)*y(k,53) + rxt(k,407)*y(k,54) + rxt(k,408)*y(k,56) &
                      + rxt(k,450)*y(k,60)
         mat(k,154) = rxt(k,496)*y(k,54)
         mat(k,185) = mat(k,185) + rxt(k,468)*y(k,54)
         mat(k,1747) = mat(k,1747) + rxt(k,434)*y(k,43)
         mat(k,1144) = mat(k,1144) + rxt(k,341)*y(k,105)
         mat(k,1223) = mat(k,1223) + rxt(k,429)*y(k,43)
         mat(k,974) = mat(k,974) + rxt(k,526)*y(k,14) + rxt(k,401)*y(k,33) &
                      + rxt(k,402)*y(k,35) + rxt(k,494)*y(k,36) + rxt(k,465)*y(k,40) &
                      + rxt(k,407)*y(k,43) + rxt(k,496)*y(k,45) + rxt(k,468)*y(k,46)
         mat(k,1426) = mat(k,1426) + rxt(k,409)*y(k,32) + rxt(k,408)*y(k,43) &
                      + rxt(k,343)*y(k,105)
         mat(k,1669) = rxt(k,450)*y(k,43) + rxt(k,336)*y(k,105)
         mat(k,784) = rxt(k,119)*y(k,42)
         mat(k,623) = rxt(k,124)*y(k,32)
         mat(k,1018) = rxt(k,133)*y(k,70)
         mat(k,35) = rxt(k,133)*y(k,69) + rxt(k,134)*y(k,110)
         mat(k,298) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77) &
                      + rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92) &
                      + rxt(k,144)*y(k,99)
         mat(k,827) = rxt(k,184)*y(k,71)
         mat(k,1187) = rxt(k,198)*y(k,71)
         mat(k,1557) = rxt(k,167)*y(k,71)
         mat(k,875) = rxt(k,237)*y(k,71)
         mat(k,918) = rxt(k,219)*y(k,71)
         mat(k,1059) = mat(k,1059) + rxt(k,202)*y(k,71)
         mat(k,575) = rxt(k,384)*y(k,33) + rxt(k,500)*y(k,36) + rxt(k,477)*y(k,40) &
                      + 2.000_r8*rxt(k,385)*y(k,110)
         mat(k,1346) = rxt(k,144)*y(k,71)
         mat(k,143) = rxt(k,323)*y(k,110)
         mat(k,755) = mat(k,755) + .300_r8*rxt(k,530)*y(k,19)
         mat(k,1590) = rxt(k,337)*y(k,24) + rxt(k,341)*y(k,52) + rxt(k,343)*y(k,56) &
                      + rxt(k,336)*y(k,60)
         mat(k,1626) = rxt(k,351)*y(k,40) + rxt(k,352)*y(k,42) + rxt(k,350)*y(k,110)
         mat(k,1800) = rxt(k,504)*y(k,31) + rxt(k,134)*y(k,70) + 2.000_r8*rxt(k,385) &
                      *y(k,95) + rxt(k,323)*y(k,100) + rxt(k,350)*y(k,106)
         mat(k,1610) = -(rxt(k,336)*y(k,60) + rxt(k,337)*y(k,24) + rxt(k,338)*y(k,28) &
                      + rxt(k,339)*y(k,32) + rxt(k,340)*y(k,40) + rxt(k,341)*y(k,52) &
                      + rxt(k,342)*y(k,54) + rxt(k,343)*y(k,56))
         mat(k,1689) = -rxt(k,336)*y(k,105)
         mat(k,1534) = -rxt(k,337)*y(k,105)
         mat(k,1408) = -rxt(k,338)*y(k,105)
         mat(k,1725) = -rxt(k,339)*y(k,105)
         mat(k,1491) = -rxt(k,340)*y(k,105)
         mat(k,1164) = -rxt(k,341)*y(k,105)
         mat(k,994) = -rxt(k,342)*y(k,105)
         mat(k,1446) = -rxt(k,343)*y(k,105)
         mat(k,1725) = mat(k,1725) + rxt(k,112)*y(k,65) + rxt(k,282)*y(k,81) &
                      + rxt(k,326)*y(k,101)
         mat(k,432) = rxt(k,349)*y(k,106)
         mat(k,804) = rxt(k,112)*y(k,32)
         mat(k,730) = rxt(k,282)*y(k,32)
         mat(k,412) = rxt(k,326)*y(k,32)
         mat(k,1646) = rxt(k,349)*y(k,33) + rxt(k,350)*y(k,110)
         mat(k,1820) = rxt(k,350)*y(k,106)
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
         mat(k,1647) = -(rxt(k,131)*y(k,51) + rxt(k,345)*y(k,60) + rxt(k,346)*y(k,24) &
                      + rxt(k,347)*y(k,28) + (rxt(k,348) + rxt(k,349)) * y(k,33) &
                      + rxt(k,350)*y(k,110) + rxt(k,351)*y(k,40) + rxt(k,352)*y(k,42) &
                      + rxt(k,354)*y(k,52) + rxt(k,355)*y(k,54) + rxt(k,356)*y(k,96) &
                      + rxt(k,357)*y(k,55) + rxt(k,358)*y(k,56))
         mat(k,1768) = -rxt(k,131)*y(k,106)
         mat(k,1690) = -rxt(k,345)*y(k,106)
         mat(k,1535) = -rxt(k,346)*y(k,106)
         mat(k,1409) = -rxt(k,347)*y(k,106)
         mat(k,433) = -(rxt(k,348) + rxt(k,349)) * y(k,106)
         mat(k,1821) = -rxt(k,350)*y(k,106)
         mat(k,1492) = -rxt(k,351)*y(k,106)
         mat(k,1120) = -rxt(k,352)*y(k,106)
         mat(k,1165) = -rxt(k,354)*y(k,106)
         mat(k,995) = -rxt(k,355)*y(k,106)
         mat(k,226) = -rxt(k,356)*y(k,106)
         mat(k,1286) = -rxt(k,357)*y(k,106)
         mat(k,1447) = -rxt(k,358)*y(k,106)
         mat(k,995) = mat(k,995) + rxt(k,320)*y(k,98)
         mat(k,1447) = mat(k,1447) + rxt(k,129)*y(k,69)
         mat(k,1038) = rxt(k,129)*y(k,56)
         mat(k,1322) = rxt(k,320)*y(k,54)
         mat(k,244) = -(rxt(k,573)*y(k,28) + rxt(k,575)*y(k,79) + rxt(k,576)*y(k,55))
         mat(k,1378) = -rxt(k,573)*y(k,107)
         mat(k,161) = -rxt(k,575)*y(k,107)
         mat(k,1256) = -rxt(k,576)*y(k,107)
         mat(k,960) = rxt(k,566)*y(k,80) + rxt(k,567)*y(k,94) + rxt(k,579)*y(k,108) &
                      + rxt(k,585)*y(k,109)
         mat(k,1013) = rxt(k,577)*y(k,108) + rxt(k,582)*y(k,109)
         mat(k,88) = rxt(k,566)*y(k,54)
         mat(k,177) = rxt(k,567)*y(k,54)
         mat(k,67) = rxt(k,579)*y(k,54) + rxt(k,577)*y(k,69)
         mat(k,62) = rxt(k,585)*y(k,54) + rxt(k,582)*y(k,69)
         mat(k,65) = -(rxt(k,577)*y(k,69) + rxt(k,579)*y(k,54) + rxt(k,580)*y(k,55))
         mat(k,1003) = -rxt(k,577)*y(k,108)
         mat(k,946) = -rxt(k,579)*y(k,108)
         mat(k,1250) = -rxt(k,580)*y(k,108)
         mat(k,1003) = mat(k,1003) + rxt(k,581)*y(k,109)
         mat(k,59) = rxt(k,581)*y(k,69)
         mat(k,58) = -((rxt(k,581) + rxt(k,582)) * y(k,69) + rxt(k,585)*y(k,54))
         mat(k,1002) = -(rxt(k,581) + rxt(k,582)) * y(k,109)
         mat(k,945) = -rxt(k,585)*y(k,109)
         mat(k,1825) = -(rxt(k,102)*y(k,61) + rxt(k,113)*y(k,67) + rxt(k,114)*y(k,65) &
                      + rxt(k,134)*y(k,70) + rxt(k,135)*y(k,74) + rxt(k,138)*y(k,76) &
                      + rxt(k,283)*y(k,81) + rxt(k,289)*y(k,85) + rxt(k,291)*y(k,83) &
                      + rxt(k,300)*y(k,90) + rxt(k,302)*y(k,91) + rxt(k,308)*y(k,93) &
                      + (rxt(k,322) + rxt(k,323)) * y(k,100) + rxt(k,324)*y(k,99) &
                      + rxt(k,332)*y(k,103) + rxt(k,350)*y(k,106) + rxt(k,360)*y(k,88) &
                      + rxt(k,361)*y(k,92) + rxt(k,362)*y(k,89) + rxt(k,367)*y(k,78) &
                      + rxt(k,369)*y(k,72) + rxt(k,371)*y(k,73) + rxt(k,373)*y(k,75) &
                      + rxt(k,385)*y(k,95) + rxt(k,504)*y(k,31))
         mat(k,615) = -rxt(k,102)*y(k,110)
         mat(k,515) = -rxt(k,113)*y(k,110)
         mat(k,809) = -rxt(k,114)*y(k,110)
         mat(k,38) = -rxt(k,134)*y(k,110)
         mat(k,25) = -rxt(k,135)*y(k,110)
         mat(k,29) = -rxt(k,138)*y(k,110)
         mat(k,734) = -rxt(k,283)*y(k,110)
         mat(k,535) = -rxt(k,289)*y(k,110)
         mat(k,672) = -rxt(k,291)*y(k,110)
         mat(k,50) = -rxt(k,300)*y(k,110)
         mat(k,122) = -rxt(k,302)*y(k,110)
         mat(k,115) = -rxt(k,308)*y(k,110)
         mat(k,149) = -(rxt(k,322) + rxt(k,323)) * y(k,110)
         mat(k,1371) = -rxt(k,324)*y(k,110)
         mat(k,296) = -rxt(k,332)*y(k,110)
         mat(k,1651) = -rxt(k,350)*y(k,110)
         mat(k,900) = -rxt(k,360)*y(k,110)
         mat(k,1084) = -rxt(k,361)*y(k,110)
         mat(k,943) = -rxt(k,362)*y(k,110)
         mat(k,138) = -rxt(k,367)*y(k,110)
         mat(k,107) = -rxt(k,369)*y(k,110)
         mat(k,852) = -rxt(k,371)*y(k,110)
         mat(k,1212) = -rxt(k,373)*y(k,110)
         mat(k,594) = -rxt(k,385)*y(k,110)
         mat(k,175) = -rxt(k,504)*y(k,110)
         mat(k,560) = rxt(k,527)*y(k,104)
         mat(k,44) = rxt(k,530)*y(k,104)
         mat(k,1730) = rxt(k,405)*y(k,43) + rxt(k,339)*y(k,105)
         mat(k,437) = rxt(k,410)*y(k,104) + rxt(k,348)*y(k,106)
         mat(k,85) = rxt(k,411)*y(k,104)
         mat(k,198) = rxt(k,495)*y(k,104)
         mat(k,1496) = (rxt(k,549)+rxt(k,554))*y(k,45) + (rxt(k,542)+rxt(k,548) &
                       +rxt(k,553))*y(k,46) + rxt(k,101)*y(k,62) + rxt(k,466)*y(k,104) &
                      + rxt(k,340)*y(k,105)
         mat(k,1124) = rxt(k,290)*y(k,85) + rxt(k,442)*y(k,104)
         mat(k,707) = rxt(k,405)*y(k,32) + rxt(k,412)*y(k,104)
         mat(k,57) = rxt(k,418)*y(k,104)
         mat(k,158) = (rxt(k,549)+rxt(k,554))*y(k,40)
         mat(k,190) = (rxt(k,542)+rxt(k,548)+rxt(k,553))*y(k,40) + rxt(k,469)*y(k,104)
         mat(k,127) = rxt(k,335)*y(k,104)
         mat(k,240) = rxt(k,288)*y(k,84)
         mat(k,1772) = rxt(k,118)*y(k,67)
         mat(k,1169) = rxt(k,115)*y(k,67)
         mat(k,615) = mat(k,615) + 3.000_r8*rxt(k,190)*y(k,73) + 4.000_r8*rxt(k,142) &
                      *y(k,75) + 5.000_r8*rxt(k,172)*y(k,77) + 2.000_r8*rxt(k,225) &
                      *y(k,89) + rxt(k,207)*y(k,92)
         mat(k,370) = rxt(k,101)*y(k,40) + 4.000_r8*rxt(k,194)*y(k,73) &
                      + 5.000_r8*rxt(k,159)*y(k,75) + 6.000_r8*rxt(k,177)*y(k,77) &
                      + rxt(k,247)*y(k,88) + 3.000_r8*rxt(k,229)*y(k,89) &
                      + 2.000_r8*rxt(k,212)*y(k,92) + rxt(k,154)*y(k,99)
         mat(k,354) = 3.000_r8*rxt(k,182)*y(k,73) + (4.000_r8*rxt(k,176) &
                       +4.000_r8*rxt(k,262))*y(k,75) + (5.000_r8*rxt(k,164) &
                       +5.000_r8*rxt(k,263))*y(k,77) + 2.000_r8*rxt(k,217)*y(k,89) &
                      + rxt(k,200)*y(k,92)
         mat(k,809) = mat(k,809) + 3.000_r8*rxt(k,181)*y(k,73) + (4.000_r8*rxt(k,165) &
                       +4.000_r8*rxt(k,253))*y(k,75) + (5.000_r8*rxt(k,163) &
                       +5.000_r8*rxt(k,260))*y(k,77) + 2.000_r8*rxt(k,216)*y(k,89) &
                      + rxt(k,199)*y(k,92)
         mat(k,325) = 5.000_r8*rxt(k,189)*y(k,73) + (6.000_r8*rxt(k,242) &
                       +6.000_r8*rxt(k,266))*y(k,75) + (7.000_r8*rxt(k,171) &
                       +7.000_r8*rxt(k,254))*y(k,77) + 2.000_r8*rxt(k,241)*y(k,88) &
                      + 4.000_r8*rxt(k,224)*y(k,89) + 3.000_r8*rxt(k,206)*y(k,92) &
                      + 2.000_r8*rxt(k,148)*y(k,99)
         mat(k,515) = mat(k,515) + rxt(k,118)*y(k,51) + rxt(k,115)*y(k,52) &
                      + 4.000_r8*rxt(k,191)*y(k,73) + (5.000_r8*rxt(k,153) &
                       +5.000_r8*rxt(k,255))*y(k,75) + (6.000_r8*rxt(k,173) &
                       +6.000_r8*rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + 3.000_r8*rxt(k,226)*y(k,89) + 2.000_r8*rxt(k,208)*y(k,92) &
                      + rxt(k,150)*y(k,99)
         mat(k,647) = 3.000_r8*rxt(k,186)*y(k,73) + 4.000_r8*rxt(k,220)*y(k,75) &
                      + 5.000_r8*rxt(k,169)*y(k,77) + 2.000_r8*rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92)
         mat(k,1042) = rxt(k,133)*y(k,70) + 2.000_r8*rxt(k,377)*y(k,72) &
                      + 3.000_r8*rxt(k,378)*y(k,73) + 4.000_r8*rxt(k,136)*y(k,75) &
                      + 5.000_r8*rxt(k,139)*y(k,77) + rxt(k,376)*y(k,78) &
                      + 2.000_r8*rxt(k,298)*y(k,89) + 3.000_r8*rxt(k,299)*y(k,90) &
                      + rxt(k,303)*y(k,92) + rxt(k,321)*y(k,100)
         mat(k,38) = mat(k,38) + rxt(k,133)*y(k,69)
         mat(k,310) = 3.000_r8*rxt(k,184)*y(k,73) + 4.000_r8*rxt(k,198)*y(k,75) &
                      + 5.000_r8*rxt(k,167)*y(k,77) + 2.000_r8*rxt(k,219)*y(k,89) &
                      + rxt(k,202)*y(k,92)
         mat(k,107) = mat(k,107) + 2.000_r8*rxt(k,377)*y(k,69)
         mat(k,852) = mat(k,852) + 3.000_r8*rxt(k,190)*y(k,61) + 4.000_r8*rxt(k,194) &
                      *y(k,62) + 3.000_r8*rxt(k,182)*y(k,63) + 3.000_r8*rxt(k,181) &
                      *y(k,65) + 5.000_r8*rxt(k,189)*y(k,66) + 4.000_r8*rxt(k,191) &
                      *y(k,67) + 3.000_r8*rxt(k,186)*y(k,68) + 3.000_r8*rxt(k,378) &
                      *y(k,69) + 3.000_r8*rxt(k,184)*y(k,71) + 3.000_r8*rxt(k,196) &
                      *y(k,81) + 4.000_r8*rxt(k,192)*y(k,82) + 3.000_r8*rxt(k,183) &
                      *y(k,83) + 5.000_r8*rxt(k,195)*y(k,84) + 4.000_r8*rxt(k,188) &
                      *y(k,85) + 3.000_r8*rxt(k,193)*y(k,86) + 3.000_r8*rxt(k,180) &
                      *y(k,87) + 3.000_r8*rxt(k,185)*y(k,98)
         mat(k,1212) = mat(k,1212) + 4.000_r8*rxt(k,142)*y(k,61) + 5.000_r8*rxt(k,159) &
                      *y(k,62) + (4.000_r8*rxt(k,176)+4.000_r8*rxt(k,262))*y(k,63) + ( &
                      + 4.000_r8*rxt(k,165)+4.000_r8*rxt(k,253))*y(k,65) + ( &
                      + 6.000_r8*rxt(k,242)+6.000_r8*rxt(k,266))*y(k,66) + ( &
                      + 5.000_r8*rxt(k,153)+5.000_r8*rxt(k,255))*y(k,67) &
                      + 4.000_r8*rxt(k,220)*y(k,68) + 4.000_r8*rxt(k,136)*y(k,69) &
                      + 4.000_r8*rxt(k,198)*y(k,71) + 4.000_r8*rxt(k,161)*y(k,81) &
                      + 5.000_r8*rxt(k,157)*y(k,82) + (4.000_r8*rxt(k,187) &
                       +4.000_r8*rxt(k,259))*y(k,83) + 6.000_r8*rxt(k,160)*y(k,84) + ( &
                      + 5.000_r8*rxt(k,231)+5.000_r8*rxt(k,257))*y(k,85) &
                      + 4.000_r8*rxt(k,158)*y(k,86) + (4.000_r8*rxt(k,141) &
                       +4.000_r8*rxt(k,264))*y(k,87) + 4.000_r8*rxt(k,209)*y(k,98)
         mat(k,1582) = 5.000_r8*rxt(k,172)*y(k,61) + 6.000_r8*rxt(k,177)*y(k,62) + ( &
                      + 5.000_r8*rxt(k,164)+5.000_r8*rxt(k,263))*y(k,63) + ( &
                      + 5.000_r8*rxt(k,163)+5.000_r8*rxt(k,260))*y(k,65) + ( &
                      + 7.000_r8*rxt(k,171)+7.000_r8*rxt(k,254))*y(k,66) + ( &
                      + 6.000_r8*rxt(k,173)+6.000_r8*rxt(k,256))*y(k,67) &
                      + 5.000_r8*rxt(k,169)*y(k,68) + 5.000_r8*rxt(k,139)*y(k,69) &
                      + 5.000_r8*rxt(k,167)*y(k,71) + 5.000_r8*rxt(k,179)*y(k,81) &
                      + 6.000_r8*rxt(k,174)*y(k,82) + (5.000_r8*rxt(k,166) &
                       +5.000_r8*rxt(k,261))*y(k,83) + 7.000_r8*rxt(k,178)*y(k,84) + ( &
                      + 6.000_r8*rxt(k,170)+6.000_r8*rxt(k,258))*y(k,85) &
                      + 5.000_r8*rxt(k,175)*y(k,86) + (5.000_r8*rxt(k,162) &
                       +5.000_r8*rxt(k,265))*y(k,87) + 5.000_r8*rxt(k,168)*y(k,98)
         mat(k,138) = mat(k,138) + rxt(k,376)*y(k,69)
         mat(k,734) = mat(k,734) + 3.000_r8*rxt(k,196)*y(k,73) + 4.000_r8*rxt(k,161) &
                      *y(k,75) + 5.000_r8*rxt(k,179)*y(k,77) + 2.000_r8*rxt(k,232) &
                      *y(k,89) + rxt(k,214)*y(k,92)
         mat(k,384) = 4.000_r8*rxt(k,192)*y(k,73) + 5.000_r8*rxt(k,157)*y(k,75) &
                      + 6.000_r8*rxt(k,174)*y(k,77) + rxt(k,245)*y(k,88) &
                      + 3.000_r8*rxt(k,227)*y(k,89) + 2.000_r8*rxt(k,210)*y(k,92) &
                      + rxt(k,151)*y(k,99)
         mat(k,672) = mat(k,672) + 3.000_r8*rxt(k,183)*y(k,73) + (4.000_r8*rxt(k,187) &
                       +4.000_r8*rxt(k,259))*y(k,75) + (5.000_r8*rxt(k,166) &
                       +5.000_r8*rxt(k,261))*y(k,77) + 2.000_r8*rxt(k,218)*y(k,89) &
                      + rxt(k,201)*y(k,92)
         mat(k,494) = rxt(k,288)*y(k,50) + 5.000_r8*rxt(k,195)*y(k,73) &
                      + 6.000_r8*rxt(k,160)*y(k,75) + 7.000_r8*rxt(k,178)*y(k,77) &
                      + 2.000_r8*rxt(k,248)*y(k,88) + 4.000_r8*rxt(k,230)*y(k,89) &
                      + 3.000_r8*rxt(k,213)*y(k,92) + 2.000_r8*rxt(k,155)*y(k,99)
         mat(k,535) = mat(k,535) + rxt(k,290)*y(k,42) + 4.000_r8*rxt(k,188)*y(k,73) + ( &
                      + 5.000_r8*rxt(k,231)+5.000_r8*rxt(k,257))*y(k,75) + ( &
                      + 6.000_r8*rxt(k,170)+6.000_r8*rxt(k,258))*y(k,77) + rxt(k,240) &
                      *y(k,88) + 3.000_r8*rxt(k,223)*y(k,89) + 2.000_r8*rxt(k,205) &
                      *y(k,92) + rxt(k,147)*y(k,99)
         mat(k,400) = 3.000_r8*rxt(k,193)*y(k,73) + 4.000_r8*rxt(k,158)*y(k,75) &
                      + 5.000_r8*rxt(k,175)*y(k,77) + 2.000_r8*rxt(k,228)*y(k,89) &
                      + rxt(k,211)*y(k,92)
         mat(k,451) = 3.000_r8*rxt(k,180)*y(k,73) + (4.000_r8*rxt(k,141) &
                       +4.000_r8*rxt(k,264))*y(k,75) + (5.000_r8*rxt(k,162) &
                       +5.000_r8*rxt(k,265))*y(k,77) + 2.000_r8*rxt(k,215)*y(k,89) &
                      + rxt(k,197)*y(k,92)
         mat(k,900) = mat(k,900) + rxt(k,247)*y(k,62) + 2.000_r8*rxt(k,241)*y(k,66) &
                      + rxt(k,244)*y(k,67) + rxt(k,245)*y(k,82) + 2.000_r8*rxt(k,248) &
                      *y(k,84) + rxt(k,240)*y(k,85)
         mat(k,943) = mat(k,943) + 2.000_r8*rxt(k,225)*y(k,61) + 3.000_r8*rxt(k,229) &
                      *y(k,62) + 2.000_r8*rxt(k,217)*y(k,63) + 2.000_r8*rxt(k,216) &
                      *y(k,65) + 4.000_r8*rxt(k,224)*y(k,66) + 3.000_r8*rxt(k,226) &
                      *y(k,67) + 2.000_r8*rxt(k,222)*y(k,68) + 2.000_r8*rxt(k,298) &
                      *y(k,69) + 2.000_r8*rxt(k,219)*y(k,71) + 2.000_r8*rxt(k,232) &
                      *y(k,81) + 3.000_r8*rxt(k,227)*y(k,82) + 2.000_r8*rxt(k,218) &
                      *y(k,83) + 4.000_r8*rxt(k,230)*y(k,84) + 3.000_r8*rxt(k,223) &
                      *y(k,85) + 2.000_r8*rxt(k,228)*y(k,86) + 2.000_r8*rxt(k,215) &
                      *y(k,87) + 2.000_r8*rxt(k,221)*y(k,98)
         mat(k,50) = mat(k,50) + 3.000_r8*rxt(k,299)*y(k,69)
         mat(k,1084) = mat(k,1084) + rxt(k,207)*y(k,61) + 2.000_r8*rxt(k,212)*y(k,62) &
                      + rxt(k,200)*y(k,63) + rxt(k,199)*y(k,65) + 3.000_r8*rxt(k,206) &
                      *y(k,66) + 2.000_r8*rxt(k,208)*y(k,67) + rxt(k,204)*y(k,68) &
                      + rxt(k,303)*y(k,69) + rxt(k,202)*y(k,71) + rxt(k,214)*y(k,81) &
                      + 2.000_r8*rxt(k,210)*y(k,82) + rxt(k,201)*y(k,83) &
                      + 3.000_r8*rxt(k,213)*y(k,84) + 2.000_r8*rxt(k,205)*y(k,85) &
                      + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87) + rxt(k,203)*y(k,98)
         mat(k,1326) = 3.000_r8*rxt(k,185)*y(k,73) + 4.000_r8*rxt(k,209)*y(k,75) &
                      + 5.000_r8*rxt(k,168)*y(k,77) + 2.000_r8*rxt(k,221)*y(k,89) &
                      + rxt(k,203)*y(k,92)
         mat(k,1371) = mat(k,1371) + rxt(k,154)*y(k,62) + 2.000_r8*rxt(k,148)*y(k,66) &
                      + rxt(k,150)*y(k,67) + rxt(k,151)*y(k,82) + 2.000_r8*rxt(k,155) &
                      *y(k,84) + rxt(k,147)*y(k,85)
         mat(k,149) = mat(k,149) + rxt(k,321)*y(k,69)
         mat(k,776) = rxt(k,527)*y(k,14) + rxt(k,530)*y(k,19) + rxt(k,410)*y(k,33) &
                      + rxt(k,411)*y(k,35) + rxt(k,495)*y(k,36) + rxt(k,466)*y(k,40) &
                      + rxt(k,442)*y(k,42) + rxt(k,412)*y(k,43) + rxt(k,418)*y(k,44) &
                      + rxt(k,469)*y(k,46) + rxt(k,335)*y(k,47) + 2.000_r8*rxt(k,415) &
                      *y(k,104)
         mat(k,1615) = rxt(k,339)*y(k,32) + rxt(k,340)*y(k,40)
         mat(k,1651) = mat(k,1651) + rxt(k,348)*y(k,33)
      end do
      end subroutine nlnmat09
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
         mat(k, 4) = mat(k, 4) + lmat(k, 4)
         mat(k, 5) = mat(k, 5) + lmat(k, 5)
         mat(k, 7) = lmat(k, 7)
         mat(k, 8) = lmat(k, 8)
         mat(k, 9) = lmat(k, 9)
         mat(k, 10) = lmat(k, 10)
         mat(k, 11) = lmat(k, 11)
         mat(k, 12) = mat(k, 12) + lmat(k, 12)
         mat(k, 13) = mat(k, 13) + lmat(k, 13)
         mat(k, 15) = mat(k, 15) + lmat(k, 15)
         mat(k, 16) = lmat(k, 16)
         mat(k, 17) = lmat(k, 17)
         mat(k, 18) = lmat(k, 18)
         mat(k, 19) = lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = mat(k, 22) + lmat(k, 22)
         mat(k, 26) = mat(k, 26) + lmat(k, 26)
         mat(k, 30) = mat(k, 30) + lmat(k, 30)
         mat(k, 31) = mat(k, 31) + lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = mat(k, 33) + lmat(k, 33)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 42) = mat(k, 42) + lmat(k, 42)
         mat(k, 43) = lmat(k, 43)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 52) = lmat(k, 52)
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 55) = lmat(k, 55)
         mat(k, 58) = mat(k, 58) + lmat(k, 58)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = lmat(k, 103)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 110) = lmat(k, 110)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = lmat(k, 117)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 126) = lmat(k, 126)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 147) = lmat(k, 147)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 169) = lmat(k, 169)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 197) = mat(k, 197) + lmat(k, 197)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 201) = mat(k, 201) + lmat(k, 201)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 237) = lmat(k, 237)
         mat(k, 239) = lmat(k, 239)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 249) = lmat(k, 249)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 268) = lmat(k, 268)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 276) = lmat(k, 276)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 312) = lmat(k, 312)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 413) = lmat(k, 413)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 439) = lmat(k, 439)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 462) = mat(k, 462) + lmat(k, 462)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 494) = mat(k, 494) + lmat(k, 494)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 501) = lmat(k, 501)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 539) = lmat(k, 539)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 562) = mat(k, 562) + lmat(k, 562)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 567) = lmat(k, 567)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = mat(k, 570) + lmat(k, 570)
         mat(k, 572) = lmat(k, 572)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 574) = lmat(k, 574)
         mat(k, 575) = mat(k, 575) + lmat(k, 575)
         mat(k, 577) = mat(k, 577) + lmat(k, 577)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 597) = mat(k, 597) + lmat(k, 597)
         mat(k, 619) = mat(k, 619) + lmat(k, 619)
         mat(k, 635) = lmat(k, 635)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 658) = lmat(k, 658)
         mat(k, 663) = mat(k, 663) + lmat(k, 663)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 712) = mat(k, 712) + lmat(k, 712)
         mat(k, 718) = lmat(k, 718)
         mat(k, 721) = mat(k, 721) + lmat(k, 721)
         mat(k, 735) = lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 747) = mat(k, 747) + lmat(k, 747)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 755) = mat(k, 755) + lmat(k, 755)
         mat(k, 773) = mat(k, 773) + lmat(k, 773)
         mat(k, 776) = mat(k, 776) + lmat(k, 776)
         mat(k, 785) = mat(k, 785) + lmat(k, 785)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 805) = lmat(k, 805)
         mat(k, 810) = lmat(k, 810)
         mat(k, 829) = mat(k, 829) + lmat(k, 829)
         mat(k, 852) = mat(k, 852) + lmat(k, 852)
         mat(k, 853) = lmat(k, 853)
         mat(k, 878) = mat(k, 878) + lmat(k, 878)
         mat(k, 922) = mat(k, 922) + lmat(k, 922)
         mat(k, 945) = mat(k, 945) + lmat(k, 945)
         mat(k, 946) = mat(k, 946) + lmat(k, 946)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k,1023) = mat(k,1023) + lmat(k,1023)
         mat(k,1066) = mat(k,1066) + lmat(k,1066)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1107) = mat(k,1107) + lmat(k,1107)
         mat(k,1108) = lmat(k,1108)
         mat(k,1144) = mat(k,1144) + lmat(k,1144)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1168) = mat(k,1168) + lmat(k,1168)
         mat(k,1189) = lmat(k,1189)
         mat(k,1197) = mat(k,1197) + lmat(k,1197)
         mat(k,1212) = mat(k,1212) + lmat(k,1212)
         mat(k,1227) = mat(k,1227) + lmat(k,1227)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1249) = lmat(k,1249)
         mat(k,1250) = mat(k,1250) + lmat(k,1250)
         mat(k,1256) = mat(k,1256) + lmat(k,1256)
         mat(k,1260) = mat(k,1260) + lmat(k,1260)
         mat(k,1270) = mat(k,1270) + lmat(k,1270)
         mat(k,1271) = mat(k,1271) + lmat(k,1271)
         mat(k,1277) = mat(k,1277) + lmat(k,1277)
         mat(k,1279) = mat(k,1279) + lmat(k,1279)
         mat(k,1307) = mat(k,1307) + lmat(k,1307)
         mat(k,1313) = mat(k,1313) + lmat(k,1313)
         mat(k,1314) = mat(k,1314) + lmat(k,1314)
         mat(k,1349) = mat(k,1349) + lmat(k,1349)
         mat(k,1360) = mat(k,1360) + lmat(k,1360)
         mat(k,1370) = mat(k,1370) + lmat(k,1370)
         mat(k,1376) = mat(k,1376) + lmat(k,1376)
         mat(k,1393) = lmat(k,1393)
         mat(k,1403) = mat(k,1403) + lmat(k,1403)
         mat(k,1415) = mat(k,1415) + lmat(k,1415)
         mat(k,1421) = mat(k,1421) + lmat(k,1421)
         mat(k,1431) = mat(k,1431) + lmat(k,1431)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1488) = mat(k,1488) + lmat(k,1488)
         mat(k,1493) = mat(k,1493) + lmat(k,1493)
         mat(k,1494) = mat(k,1494) + lmat(k,1494)
         mat(k,1519) = mat(k,1519) + lmat(k,1519)
         mat(k,1532) = mat(k,1532) + lmat(k,1532)
         mat(k,1536) = mat(k,1536) + lmat(k,1536)
         mat(k,1567) = lmat(k,1567)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1582) = mat(k,1582) + lmat(k,1582)
         mat(k,1590) = mat(k,1590) + lmat(k,1590)
         mat(k,1596) = mat(k,1596) + lmat(k,1596)
         mat(k,1610) = mat(k,1610) + lmat(k,1610)
         mat(k,1631) = mat(k,1631) + lmat(k,1631)
         mat(k,1632) = mat(k,1632) + lmat(k,1632)
         mat(k,1647) = mat(k,1647) + lmat(k,1647)
         mat(k,1654) = mat(k,1654) + lmat(k,1654)
         mat(k,1658) = lmat(k,1658)
         mat(k,1659) = lmat(k,1659)
         mat(k,1667) = mat(k,1667) + lmat(k,1667)
         mat(k,1686) = mat(k,1686) + lmat(k,1686)
         mat(k,1691) = mat(k,1691) + lmat(k,1691)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1734) = mat(k,1734) + lmat(k,1734)
         mat(k,1750) = mat(k,1750) + lmat(k,1750)
         mat(k,1752) = mat(k,1752) + lmat(k,1752)
         mat(k,1753) = mat(k,1753) + lmat(k,1753)
         mat(k,1771) = mat(k,1771) + lmat(k,1771)
         mat(k,1790) = lmat(k,1790)
         mat(k,1795) = mat(k,1795) + lmat(k,1795)
         mat(k,1800) = mat(k,1800) + lmat(k,1800)
         mat(k,1805) = lmat(k,1805)
         mat(k,1823) = lmat(k,1823)
         mat(k,1825) = mat(k,1825) + lmat(k,1825)
         mat(k, 135) = 0._r8
         mat(k, 140) = 0._r8
         mat(k, 144) = 0._r8
         mat(k, 148) = 0._r8
         mat(k, 157) = 0._r8
         mat(k, 210) = 0._r8
         mat(k, 246) = 0._r8
         mat(k, 247) = 0._r8
         mat(k, 248) = 0._r8
         mat(k, 254) = 0._r8
         mat(k, 255) = 0._r8
         mat(k, 260) = 0._r8
         mat(k, 265) = 0._r8
         mat(k, 267) = 0._r8
         mat(k, 269) = 0._r8
         mat(k, 271) = 0._r8
         mat(k, 272) = 0._r8
         mat(k, 280) = 0._r8
         mat(k, 286) = 0._r8
         mat(k, 287) = 0._r8
         mat(k, 291) = 0._r8
         mat(k, 294) = 0._r8
         mat(k, 295) = 0._r8
         mat(k, 331) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 334) = 0._r8
         mat(k, 336) = 0._r8
         mat(k, 338) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 418) = 0._r8
         mat(k, 420) = 0._r8
         mat(k, 422) = 0._r8
         mat(k, 426) = 0._r8
         mat(k, 427) = 0._r8
         mat(k, 428) = 0._r8
         mat(k, 429) = 0._r8
         mat(k, 431) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 458) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 465) = 0._r8
         mat(k, 467) = 0._r8
         mat(k, 468) = 0._r8
         mat(k, 469) = 0._r8
         mat(k, 472) = 0._r8
         mat(k, 474) = 0._r8
         mat(k, 479) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 524) = 0._r8
         mat(k, 527) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 542) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 550) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 552) = 0._r8
         mat(k, 554) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 556) = 0._r8
         mat(k, 559) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 576) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 580) = 0._r8
         mat(k, 581) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 585) = 0._r8
         mat(k, 589) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 620) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 685) = 0._r8
         mat(k, 688) = 0._r8
         mat(k, 690) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 696) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 717) = 0._r8
         mat(k, 743) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 759) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 766) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 780) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 801) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 830) = 0._r8
         mat(k, 831) = 0._r8
         mat(k, 834) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 848) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 856) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 873) = 0._r8
         mat(k, 877) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 889) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 896) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 921) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 928) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 978) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1035) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1062) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1079) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1105) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1118) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1222) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1225) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1273) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1350) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1369) = 0._r8
         mat(k,1379) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1411) = 0._r8
         mat(k,1412) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1420) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1477) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1481) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1483) = 0._r8
         mat(k,1485) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1517) = 0._r8
         mat(k,1518) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1564) = 0._r8
         mat(k,1571) = 0._r8
         mat(k,1573) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1577) = 0._r8
         mat(k,1578) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1587) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1593) = 0._r8
         mat(k,1594) = 0._r8
         mat(k,1597) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1601) = 0._r8
         mat(k,1602) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1611) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1617) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1624) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1629) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1636) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1639) = 0._r8
         mat(k,1640) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1671) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1673) = 0._r8
         mat(k,1675) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1679) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1693) = 0._r8
         mat(k,1694) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1709) = 0._r8
         mat(k,1713) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1716) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1723) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1726) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1764) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1794) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1809) = 0._r8
         mat(k,1811) = 0._r8
         mat(k,1813) = 0._r8
         mat(k,1816) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1818) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1824) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 7) = mat(k, 7) - dti(k)
         mat(k, 10) = mat(k, 10) - dti(k)
         mat(k, 12) = mat(k, 12) - dti(k)
         mat(k, 16) = mat(k, 16) - dti(k)
         mat(k, 19) = mat(k, 19) - dti(k)
         mat(k, 22) = mat(k, 22) - dti(k)
         mat(k, 26) = mat(k, 26) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 51) = mat(k, 51) - dti(k)
         mat(k, 58) = mat(k, 58) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 116) = mat(k, 116) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 151) = mat(k, 151) - dti(k)
         mat(k, 159) = mat(k, 159) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 176) = mat(k, 176) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 191) = mat(k, 191) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 211) = mat(k, 211) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 328) = mat(k, 328) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 356) = mat(k, 356) - dti(k)
         mat(k, 371) = mat(k, 371) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 457) = mat(k, 457) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 497) = mat(k, 497) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 573) = mat(k, 573) - dti(k)
         mat(k, 597) = mat(k, 597) - dti(k)
         mat(k, 619) = mat(k, 619) - dti(k)
         mat(k, 652) = mat(k, 652) - dti(k)
         mat(k, 686) = mat(k, 686) - dti(k)
         mat(k, 712) = mat(k, 712) - dti(k)
         mat(k, 755) = mat(k, 755) - dti(k)
         mat(k, 785) = mat(k, 785) - dti(k)
         mat(k, 829) = mat(k, 829) - dti(k)
         mat(k, 878) = mat(k, 878) - dti(k)
         mat(k, 922) = mat(k, 922) - dti(k)
         mat(k, 979) = mat(k, 979) - dti(k)
         mat(k,1023) = mat(k,1023) - dti(k)
         mat(k,1066) = mat(k,1066) - dti(k)
         mat(k,1107) = mat(k,1107) - dti(k)
         mat(k,1153) = mat(k,1153) - dti(k)
         mat(k,1197) = mat(k,1197) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1277) = mat(k,1277) - dti(k)
         mat(k,1314) = mat(k,1314) - dti(k)
         mat(k,1360) = mat(k,1360) - dti(k)
         mat(k,1403) = mat(k,1403) - dti(k)
         mat(k,1442) = mat(k,1442) - dti(k)
         mat(k,1488) = mat(k,1488) - dti(k)
         mat(k,1532) = mat(k,1532) - dti(k)
         mat(k,1576) = mat(k,1576) - dti(k)
         mat(k,1610) = mat(k,1610) - dti(k)
         mat(k,1647) = mat(k,1647) - dti(k)
         mat(k,1691) = mat(k,1691) - dti(k)
         mat(k,1728) = mat(k,1728) - dti(k)
         mat(k,1771) = mat(k,1771) - dti(k)
         mat(k,1825) = mat(k,1825) - dti(k)
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
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
