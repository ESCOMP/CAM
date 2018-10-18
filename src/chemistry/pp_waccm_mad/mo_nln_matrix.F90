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
         mat(k,447) = rxt(k,487)*y(k,24)
         mat(k,1488) = rxt(k,487)*y(k,2)
         mat(k,1105) = (rxt(k,549)+rxt(k,554))*y(k,45)
         mat(k,164) = (rxt(k,549)+rxt(k,554))*y(k,40)
         mat(k,451) = -(4._r8*rxt(k,484)*y(k,2) + (rxt(k,485) + rxt(k,486) + rxt(k,487) &
                      ) * y(k,24) + rxt(k,488)*y(k,43) + rxt(k,489)*y(k,51) + rxt(k,490) &
                      *y(k,52) + rxt(k,492)*y(k,54) + rxt(k,493)*y(k,104))
         mat(k,1494) = -(rxt(k,485) + rxt(k,486) + rxt(k,487)) * y(k,2)
         mat(k,677) = -rxt(k,488)*y(k,2)
         mat(k,986) = -rxt(k,489)*y(k,2)
         mat(k,1728) = -rxt(k,490)*y(k,2)
         mat(k,1211) = -rxt(k,492)*y(k,2)
         mat(k,745) = -rxt(k,493)*y(k,2)
         mat(k,79) = rxt(k,491)*y(k,54)
         mat(k,184) = rxt(k,501)*y(k,95)
         mat(k,167) = rxt(k,496)*y(k,54)
         mat(k,1211) = mat(k,1211) + rxt(k,491)*y(k,3) + rxt(k,496)*y(k,45)
         mat(k,1250) = rxt(k,483)*y(k,59)
         mat(k,323) = rxt(k,483)*y(k,56)
         mat(k,565) = rxt(k,501)*y(k,36)
         mat(k,76) = -(rxt(k,491)*y(k,54))
         mat(k,1193) = -rxt(k,491)*y(k,3)
         mat(k,448) = rxt(k,490)*y(k,52)
         mat(k,1720) = rxt(k,490)*y(k,2)
         mat(k,535) = -(rxt(k,445)*y(k,60) + rxt(k,481)*y(k,59) + rxt(k,525)*y(k,53) &
                      + rxt(k,526)*y(k,54) + rxt(k,527)*y(k,104))
         mat(k,1374) = -rxt(k,445)*y(k,14)
         mat(k,324) = -rxt(k,481)*y(k,14)
         mat(k,1412) = -rxt(k,525)*y(k,14)
         mat(k,1212) = -rxt(k,526)*y(k,14)
         mat(k,746) = -rxt(k,527)*y(k,14)
         mat(k,253) = rxt(k,452)*y(k,24) + rxt(k,529)*y(k,51)
         mat(k,41) = .300_r8*rxt(k,530)*y(k,104)
         mat(k,1495) = rxt(k,452)*y(k,18)
         mat(k,989) = rxt(k,529)*y(k,18)
         mat(k,746) = mat(k,746) + .300_r8*rxt(k,530)*y(k,19)
         mat(k,252) = -(rxt(k,452)*y(k,24) + rxt(k,528)*y(k,43) + rxt(k,529)*y(k,51))
         mat(k,1492) = -rxt(k,452)*y(k,18)
         mat(k,674) = -rxt(k,528)*y(k,18)
         mat(k,981) = -rxt(k,529)*y(k,18)
         mat(k,40) = .700_r8*rxt(k,530)*y(k,104)
         mat(k,741) = .700_r8*rxt(k,530)*y(k,19)
         mat(k,39) = -(rxt(k,530)*y(k,104))
         mat(k,731) = -rxt(k,530)*y(k,19)
         mat(k,251) = rxt(k,528)*y(k,43)
         mat(k,667) = rxt(k,528)*y(k,18)
         mat(k,1487) = 2.000_r8*rxt(k,454)*y(k,24)
         mat(k,209) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,40) + rxt(k,458)*y(k,60)
         mat(k,1104) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,25) + (rxt(k,542) &
                       +rxt(k,548)+rxt(k,553))*y(k,46)
         mat(k,190) = (rxt(k,542)+rxt(k,548)+rxt(k,553))*y(k,40)
         mat(k,1364) = rxt(k,458)*y(k,25)
         mat(k,1485) = 2.000_r8*rxt(k,479)*y(k,24)
         mat(k,1520) = -(rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68) + rxt(k,281)*y(k,81) &
                      + rxt(k,310)*y(k,98) + rxt(k,337)*y(k,105) + rxt(k,346)*y(k,106) &
                      + rxt(k,452)*y(k,18) + (4._r8*rxt(k,453) + 4._r8*rxt(k,454) &
                      + 4._r8*rxt(k,455) + 4._r8*rxt(k,479)) * y(k,24) + rxt(k,456) &
                      *y(k,43) + rxt(k,457)*y(k,51) + rxt(k,459)*y(k,52) + rxt(k,462) &
                      *y(k,54) + (rxt(k,463) + rxt(k,464)) * y(k,104) + (rxt(k,485) &
                      + rxt(k,486) + rxt(k,487)) * y(k,2))
         mat(k,796) = -rxt(k,111)*y(k,24)
         mat(k,634) = -rxt(k,123)*y(k,24)
         mat(k,722) = -rxt(k,281)*y(k,24)
         mat(k,1311) = -rxt(k,310)*y(k,24)
         mat(k,1596) = -rxt(k,337)*y(k,24)
         mat(k,1632) = -rxt(k,346)*y(k,24)
         mat(k,258) = -rxt(k,452)*y(k,24)
         mat(k,695) = -rxt(k,456)*y(k,24)
         mat(k,1013) = -rxt(k,457)*y(k,24)
         mat(k,1756) = -rxt(k,459)*y(k,24)
         mat(k,1237) = -rxt(k,462)*y(k,24)
         mat(k,764) = -(rxt(k,463) + rxt(k,464)) * y(k,24)
         mat(k,464) = -(rxt(k,485) + rxt(k,486) + rxt(k,487)) * y(k,24)
         mat(k,218) = rxt(k,460)*y(k,54)
         mat(k,1141) = rxt(k,478)*y(k,95)
         mat(k,695) = mat(k,695) + rxt(k,450)*y(k,60)
         mat(k,196) = rxt(k,468)*y(k,54) + rxt(k,467)*y(k,60) + rxt(k,469)*y(k,104)
         mat(k,1237) = mat(k,1237) + rxt(k,460)*y(k,25) + rxt(k,468)*y(k,46)
         mat(k,1275) = rxt(k,451)*y(k,60)
         mat(k,1399) = rxt(k,450)*y(k,43) + rxt(k,467)*y(k,46) + rxt(k,451)*y(k,56)
         mat(k,582) = rxt(k,478)*y(k,40)
         mat(k,764) = mat(k,764) + rxt(k,469)*y(k,46)
         mat(k,211) = -(rxt(k,458)*y(k,60) + rxt(k,460)*y(k,54) + rxt(k,461)*y(k,104) &
                      + (rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,40))
         mat(k,1368) = -rxt(k,458)*y(k,25)
         mat(k,1203) = -rxt(k,460)*y(k,25)
         mat(k,740) = -rxt(k,461)*y(k,25)
         mat(k,1109) = -(rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,25)
         mat(k,1490) = rxt(k,459)*y(k,52)
         mat(k,1722) = rxt(k,459)*y(k,24)
         mat(k,142) = -((rxt(k,532) + rxt(k,536)) * y(k,104))
         mat(k,735) = -(rxt(k,532) + rxt(k,536)) * y(k,27)
         mat(k,530) = rxt(k,525)*y(k,53) + rxt(k,526)*y(k,54) + rxt(k,481)*y(k,59) &
                      + rxt(k,445)*y(k,60) + rxt(k,527)*y(k,104)
         mat(k,1067) = rxt(k,573)*y(k,107)
         mat(k,1407) = rxt(k,525)*y(k,14)
         mat(k,1198) = rxt(k,526)*y(k,14)
         mat(k,320) = rxt(k,481)*y(k,14)
         mat(k,1366) = rxt(k,445)*y(k,14)
         mat(k,735) = mat(k,735) + rxt(k,527)*y(k,14)
         mat(k,198) = rxt(k,573)*y(k,28)
         mat(k,4) = -(rxt(k,506)*y(k,95))
         mat(k,555) = -rxt(k,506)*y(k,29)
         mat(k,12) = -(rxt(k,507)*y(k,95))
         mat(k,556) = -rxt(k,507)*y(k,30)
         mat(k,1086) = -(rxt(k,307)*y(k,93) + rxt(k,311)*y(k,98) + rxt(k,325)*y(k,101) &
                      + rxt(k,330)*y(k,102) + rxt(k,338)*y(k,105) + rxt(k,347) &
                      *y(k,106) + rxt(k,363)*y(k,88) + rxt(k,573)*y(k,107))
         mat(k,128) = -rxt(k,307)*y(k,28)
         mat(k,1301) = -rxt(k,311)*y(k,28)
         mat(k,401) = -rxt(k,325)*y(k,28)
         mat(k,101) = -rxt(k,330)*y(k,28)
         mat(k,1586) = -rxt(k,338)*y(k,28)
         mat(k,1622) = -rxt(k,347)*y(k,28)
         mat(k,875) = -rxt(k,363)*y(k,28)
         mat(k,204) = -rxt(k,573)*y(k,28)
         mat(k,1510) = rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68)
         mat(k,145) = (rxt(k,532)+rxt(k,536))*y(k,104)
         mat(k,1701) = rxt(k,112)*y(k,65)
         mat(k,1131) = rxt(k,125)*y(k,68)
         mat(k,1171) = rxt(k,119)*y(k,65)
         mat(k,1003) = rxt(k,275)*y(k,65) + (rxt(k,117)+rxt(k,118))*y(k,67)
         mat(k,1746) = rxt(k,276)*y(k,65) + (rxt(k,115)+rxt(k,116))*y(k,67)
         mat(k,1227) = rxt(k,120)*y(k,65)
         mat(k,1665) = rxt(k,121)*y(k,65)
         mat(k,1265) = rxt(k,127)*y(k,68)
         mat(k,1389) = (rxt(k,109)+rxt(k,110))*y(k,65) + rxt(k,122)*y(k,68)
         mat(k,786) = rxt(k,111)*y(k,24) + rxt(k,112)*y(k,32) + rxt(k,119)*y(k,42) &
                      + rxt(k,275)*y(k,51) + rxt(k,276)*y(k,52) + rxt(k,120)*y(k,54) &
                      + rxt(k,121)*y(k,55) + (rxt(k,109)+rxt(k,110))*y(k,60) &
                      + rxt(k,181)*y(k,73) + (rxt(k,165)+rxt(k,253))*y(k,75) + ( &
                      + rxt(k,163)+rxt(k,260))*y(k,77) + rxt(k,234)*y(k,88) &
                      + rxt(k,216)*y(k,89) + rxt(k,199)*y(k,92) + rxt(k,251)*y(k,99)
         mat(k,312) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77) + rxt(k,241)*y(k,88) &
                      + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92) + rxt(k,148)*y(k,99)
         mat(k,501) = (rxt(k,117)+rxt(k,118))*y(k,51) + (rxt(k,115)+rxt(k,116)) &
                      *y(k,52) + rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + rxt(k,226)*y(k,89) + rxt(k,208)*y(k,92) + rxt(k,150)*y(k,99)
         mat(k,625) = rxt(k,123)*y(k,24) + rxt(k,125)*y(k,40) + rxt(k,127)*y(k,56) &
                      + rxt(k,122)*y(k,60) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77) + rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92) + rxt(k,146)*y(k,99)
         mat(k,961) = rxt(k,301)*y(k,91)
         mat(k,298) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77) &
                      + rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92) &
                      + rxt(k,144)*y(k,99)
         mat(k,829) = rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67) &
                      + rxt(k,186)*y(k,68) + rxt(k,184)*y(k,71)
         mat(k,1467) = (rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67) + rxt(k,220)*y(k,68) &
                      + rxt(k,198)*y(k,71)
         mat(k,1553) = (rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67) + rxt(k,169)*y(k,68) &
                      + rxt(k,167)*y(k,71)
         mat(k,875) = mat(k,875) + rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) &
                      + rxt(k,244)*y(k,67) + rxt(k,239)*y(k,68) + rxt(k,237)*y(k,71)
         mat(k,918) = rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67) &
                      + rxt(k,222)*y(k,68) + rxt(k,219)*y(k,71)
         mat(k,135) = rxt(k,301)*y(k,69) + rxt(k,302)*y(k,110)
         mat(k,1045) = rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) + rxt(k,208)*y(k,67) &
                      + rxt(k,204)*y(k,68) + rxt(k,202)*y(k,71)
         mat(k,1346) = rxt(k,251)*y(k,65) + rxt(k,148)*y(k,66) + rxt(k,150)*y(k,67) &
                      + rxt(k,146)*y(k,68) + rxt(k,144)*y(k,71)
         mat(k,755) = (rxt(k,532)+rxt(k,536))*y(k,27)
         mat(k,1799) = rxt(k,302)*y(k,91)
         mat(k,174) = -(rxt(k,503)*y(k,33) + rxt(k,504)*y(k,110) + rxt(k,505)*y(k,42))
         mat(k,411) = -rxt(k,503)*y(k,31)
         mat(k,1775) = -rxt(k,504)*y(k,31)
         mat(k,1151) = -rxt(k,505)*y(k,31)
         mat(k,5) = 2.000_r8*rxt(k,506)*y(k,95)
         mat(k,13) = rxt(k,507)*y(k,95)
         mat(k,558) = 2.000_r8*rxt(k,506)*y(k,29) + rxt(k,507)*y(k,30)
         mat(k,1716) = -(rxt(k,100)*y(k,61) + rxt(k,112)*y(k,65) + rxt(k,124)*y(k,68) &
                      + rxt(k,282)*y(k,81) + rxt(k,304)*y(k,92) + rxt(k,312)*y(k,98) &
                      + rxt(k,326)*y(k,101) + rxt(k,339)*y(k,105) + (rxt(k,403) &
                      + rxt(k,404) + rxt(k,405)) * y(k,43) + rxt(k,406)*y(k,55) &
                      + rxt(k,409)*y(k,56))
         mat(k,607) = -rxt(k,100)*y(k,32)
         mat(k,801) = -rxt(k,112)*y(k,32)
         mat(k,639) = -rxt(k,124)*y(k,32)
         mat(k,726) = -rxt(k,282)*y(k,32)
         mat(k,1060) = -rxt(k,304)*y(k,32)
         mat(k,1316) = -rxt(k,312)*y(k,32)
         mat(k,408) = -rxt(k,326)*y(k,32)
         mat(k,1601) = -rxt(k,339)*y(k,32)
         mat(k,699) = -(rxt(k,403) + rxt(k,404) + rxt(k,405)) * y(k,32)
         mat(k,1680) = -rxt(k,406)*y(k,32)
         mat(k,1280) = -rxt(k,409)*y(k,32)
         mat(k,552) = rxt(k,527)*y(k,104)
         mat(k,146) = rxt(k,536)*y(k,104)
         mat(k,180) = rxt(k,503)*y(k,33)
         mat(k,429) = rxt(k,503)*y(k,31) + rxt(k,401)*y(k,54) + rxt(k,447)*y(k,60) &
                      + rxt(k,384)*y(k,95) + rxt(k,410)*y(k,104) + rxt(k,349)*y(k,106)
         mat(k,188) = rxt(k,501)*y(k,95)
         mat(k,1146) = rxt(k,478)*y(k,95)
         mat(k,275) = rxt(k,433)*y(k,104)
         mat(k,1242) = rxt(k,401)*y(k,33) + rxt(k,413)*y(k,104)
         mat(k,1404) = rxt(k,447)*y(k,33)
         mat(k,607) = mat(k,607) + rxt(k,190)*y(k,73) + rxt(k,142)*y(k,75) &
                      + rxt(k,172)*y(k,77)
         mat(k,363) = rxt(k,194)*y(k,73) + rxt(k,159)*y(k,75) + rxt(k,177)*y(k,77)
         mat(k,347) = rxt(k,182)*y(k,73) + rxt(k,176)*y(k,75) + rxt(k,164)*y(k,77)
         mat(k,801) = mat(k,801) + rxt(k,181)*y(k,73) + (rxt(k,165)+rxt(k,253)) &
                      *y(k,75) + (rxt(k,163)+rxt(k,260))*y(k,77)
         mat(k,318) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77)
         mat(k,507) = rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77)
         mat(k,639) = mat(k,639) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77)
         mat(k,976) = rxt(k,133)*y(k,70) + rxt(k,377)*y(k,72) + rxt(k,378)*y(k,73) &
                      + rxt(k,136)*y(k,75) + rxt(k,139)*y(k,77) + rxt(k,376)*y(k,78)
         mat(k,37) = rxt(k,133)*y(k,69)
         mat(k,303) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77)
         mat(k,109) = rxt(k,377)*y(k,69)
         mat(k,844) = rxt(k,190)*y(k,61) + rxt(k,194)*y(k,62) + rxt(k,182)*y(k,63) &
                      + rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67) &
                      + rxt(k,186)*y(k,68) + rxt(k,378)*y(k,69) + rxt(k,184)*y(k,71) &
                      + rxt(k,196)*y(k,81) + rxt(k,192)*y(k,82) + rxt(k,195)*y(k,84) &
                      + rxt(k,188)*y(k,85) + rxt(k,193)*y(k,86) + rxt(k,185)*y(k,98)
         mat(k,1482) = rxt(k,142)*y(k,61) + rxt(k,159)*y(k,62) + rxt(k,176)*y(k,63) + ( &
                      + rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67) + rxt(k,220)*y(k,68) &
                      + rxt(k,136)*y(k,69) + rxt(k,198)*y(k,71) + rxt(k,161)*y(k,81) &
                      + rxt(k,157)*y(k,82) + rxt(k,160)*y(k,84) + (rxt(k,231) &
                       +rxt(k,257))*y(k,85) + rxt(k,158)*y(k,86) + rxt(k,209)*y(k,98)
         mat(k,1568) = rxt(k,172)*y(k,61) + rxt(k,177)*y(k,62) + rxt(k,164)*y(k,63) + ( &
                      + rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67) + rxt(k,169)*y(k,68) &
                      + rxt(k,139)*y(k,69) + rxt(k,167)*y(k,71) + rxt(k,179)*y(k,81) &
                      + rxt(k,174)*y(k,82) + rxt(k,178)*y(k,84) + (rxt(k,170) &
                       +rxt(k,258))*y(k,85) + rxt(k,175)*y(k,86) + rxt(k,168)*y(k,98)
         mat(k,151) = rxt(k,376)*y(k,69)
         mat(k,726) = mat(k,726) + rxt(k,196)*y(k,73) + rxt(k,161)*y(k,75) &
                      + rxt(k,179)*y(k,77)
         mat(k,376) = rxt(k,192)*y(k,73) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77)
         mat(k,486) = rxt(k,195)*y(k,73) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77)
         mat(k,527) = rxt(k,188)*y(k,73) + (rxt(k,231)+rxt(k,257))*y(k,75) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,77)
         mat(k,393) = rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77)
         mat(k,586) = rxt(k,384)*y(k,33) + rxt(k,501)*y(k,36) + rxt(k,478)*y(k,40)
         mat(k,1316) = mat(k,1316) + rxt(k,185)*y(k,73) + rxt(k,209)*y(k,75) &
                      + rxt(k,168)*y(k,77)
         mat(k,768) = rxt(k,527)*y(k,14) + rxt(k,536)*y(k,27) + rxt(k,410)*y(k,33) &
                      + rxt(k,433)*y(k,48) + rxt(k,413)*y(k,54)
         mat(k,1637) = rxt(k,349)*y(k,33)
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
         mat(k,413) = -((rxt(k,348) + rxt(k,349)) * y(k,106) + rxt(k,384)*y(k,95) &
                      + rxt(k,401)*y(k,54) + rxt(k,410)*y(k,104) + rxt(k,447)*y(k,60) &
                      + rxt(k,503)*y(k,31))
         mat(k,1607) = -(rxt(k,348) + rxt(k,349)) * y(k,33)
         mat(k,564) = -rxt(k,384)*y(k,33)
         mat(k,1210) = -rxt(k,401)*y(k,33)
         mat(k,744) = -rxt(k,410)*y(k,33)
         mat(k,1372) = -rxt(k,447)*y(k,33)
         mat(k,176) = -rxt(k,503)*y(k,33)
         mat(k,1685) = rxt(k,403)*y(k,43)
         mat(k,676) = rxt(k,403)*y(k,32)
         mat(k,84) = -(rxt(k,402)*y(k,54) + rxt(k,411)*y(k,104) + rxt(k,448)*y(k,60))
         mat(k,1194) = -rxt(k,402)*y(k,35)
         mat(k,733) = -rxt(k,411)*y(k,35)
         mat(k,1365) = -rxt(k,448)*y(k,35)
         mat(k,669) = 2.000_r8*rxt(k,417)*y(k,43)
         mat(k,733) = mat(k,733) + 2.000_r8*rxt(k,416)*y(k,104)
         mat(k,182) = -(rxt(k,494)*y(k,54) + rxt(k,495)*y(k,104) + (rxt(k,500) &
                      + rxt(k,501)) * y(k,95))
         mat(k,1200) = -rxt(k,494)*y(k,36)
         mat(k,738) = -rxt(k,495)*y(k,36)
         mat(k,559) = -(rxt(k,500) + rxt(k,501)) * y(k,36)
         mat(k,531) = rxt(k,481)*y(k,59)
         mat(k,672) = rxt(k,482)*y(k,59)
         mat(k,321) = rxt(k,481)*y(k,14) + rxt(k,482)*y(k,43)
         mat(k,1132) = -(rxt(k,101)*y(k,62) + rxt(k,103)*y(k,61) + rxt(k,125)*y(k,68) &
                      + (rxt(k,271) + rxt(k,293)) * y(k,83) + rxt(k,284)*y(k,81) &
                      + rxt(k,313)*y(k,98) + rxt(k,340)*y(k,105) + rxt(k,351)*y(k,106) &
                      + rxt(k,465)*y(k,54) + rxt(k,466)*y(k,104) + (rxt(k,477) &
                      + rxt(k,478)) * y(k,95) + (rxt(k,542) + rxt(k,548) + rxt(k,553) &
                      ) * y(k,46) + (rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,25) &
                      + (rxt(k,549) + rxt(k,554)) * y(k,45))
         mat(k,357) = -rxt(k,101)*y(k,40)
         mat(k,600) = -rxt(k,103)*y(k,40)
         mat(k,626) = -rxt(k,125)*y(k,40)
         mat(k,654) = -(rxt(k,271) + rxt(k,293)) * y(k,40)
         mat(k,714) = -rxt(k,284)*y(k,40)
         mat(k,1302) = -rxt(k,313)*y(k,40)
         mat(k,1587) = -rxt(k,340)*y(k,40)
         mat(k,1623) = -rxt(k,351)*y(k,40)
         mat(k,1228) = -rxt(k,465)*y(k,40)
         mat(k,756) = -rxt(k,466)*y(k,40)
         mat(k,574) = -(rxt(k,477) + rxt(k,478)) * y(k,40)
         mat(k,193) = -(rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,40)
         mat(k,213) = -(rxt(k,547) + rxt(k,552) + rxt(k,557)) * y(k,40)
         mat(k,169) = -(rxt(k,549) + rxt(k,554)) * y(k,40)
         mat(k,542) = rxt(k,445)*y(k,60)
         mat(k,1511) = rxt(k,464)*y(k,104)
         mat(k,1702) = rxt(k,100)*y(k,61)
         mat(k,420) = rxt(k,447)*y(k,60)
         mat(k,87) = rxt(k,448)*y(k,60)
         mat(k,1172) = rxt(k,104)*y(k,61) + rxt(k,294)*y(k,86)
         mat(k,687) = rxt(k,449)*y(k,60)
         mat(k,193) = mat(k,193) + rxt(k,467)*y(k,60)
         mat(k,1390) = rxt(k,445)*y(k,14) + rxt(k,447)*y(k,33) + rxt(k,448)*y(k,35) &
                      + rxt(k,449)*y(k,43) + rxt(k,467)*y(k,46)
         mat(k,600) = mat(k,600) + rxt(k,100)*y(k,32) + rxt(k,104)*y(k,42)
         mat(k,341) = rxt(k,182)*y(k,73) + (rxt(k,176)+2.000_r8*rxt(k,262))*y(k,75) + ( &
                      + rxt(k,164)+2.000_r8*rxt(k,263))*y(k,77) + rxt(k,235)*y(k,88) &
                      + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92) + rxt(k,252)*y(k,99)
         mat(k,830) = rxt(k,182)*y(k,63) + rxt(k,193)*y(k,86)
         mat(k,1468) = (rxt(k,176)+2.000_r8*rxt(k,262))*y(k,63) + rxt(k,158)*y(k,86)
         mat(k,1554) = (rxt(k,164)+2.000_r8*rxt(k,263))*y(k,63) + rxt(k,175)*y(k,86)
         mat(k,386) = rxt(k,294)*y(k,42) + rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) &
                      + rxt(k,175)*y(k,77) + rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) &
                      + rxt(k,211)*y(k,92) + rxt(k,152)*y(k,99)
         mat(k,876) = rxt(k,235)*y(k,63) + rxt(k,246)*y(k,86)
         mat(k,919) = rxt(k,217)*y(k,63) + rxt(k,228)*y(k,86)
         mat(k,1046) = rxt(k,200)*y(k,63) + rxt(k,211)*y(k,86)
         mat(k,1347) = rxt(k,252)*y(k,63) + rxt(k,152)*y(k,86)
         mat(k,756) = mat(k,756) + rxt(k,464)*y(k,24)
         mat(k,173) = rxt(k,503)*y(k,33) + rxt(k,505)*y(k,42) + rxt(k,504)*y(k,110)
         mat(k,410) = rxt(k,503)*y(k,31)
         mat(k,1149) = rxt(k,505)*y(k,31)
         mat(k,1764) = rxt(k,504)*y(k,31)
         mat(k,1173) = -(rxt(k,104)*y(k,61) + rxt(k,119)*y(k,65) + rxt(k,285)*y(k,81) &
                      + rxt(k,290)*y(k,85) + rxt(k,294)*y(k,86) + rxt(k,295)*y(k,83) &
                      + rxt(k,314)*y(k,98) + rxt(k,352)*y(k,106) + rxt(k,442)*y(k,104) &
                      + rxt(k,505)*y(k,31))
         mat(k,601) = -rxt(k,104)*y(k,42)
         mat(k,788) = -rxt(k,119)*y(k,42)
         mat(k,715) = -rxt(k,285)*y(k,42)
         mat(k,520) = -rxt(k,290)*y(k,42)
         mat(k,387) = -rxt(k,294)*y(k,42)
         mat(k,655) = -rxt(k,295)*y(k,42)
         mat(k,1303) = -rxt(k,314)*y(k,42)
         mat(k,1624) = -rxt(k,352)*y(k,42)
         mat(k,757) = -rxt(k,442)*y(k,42)
         mat(k,178) = -rxt(k,505)*y(k,42)
         mat(k,543) = rxt(k,525)*y(k,53)
         mat(k,214) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,40)
         mat(k,1133) = (rxt(k,547)+rxt(k,552)+rxt(k,557))*y(k,25) + rxt(k,293)*y(k,83)
         mat(k,244) = rxt(k,137)*y(k,75) + rxt(k,140)*y(k,77) + rxt(k,288)*y(k,84) &
                      + rxt(k,292)*y(k,85)
         mat(k,1748) = rxt(k,441)*y(k,104)
         mat(k,1426) = rxt(k,525)*y(k,14)
         mat(k,831) = rxt(k,183)*y(k,83) + 2.000_r8*rxt(k,180)*y(k,87)
         mat(k,23) = rxt(k,135)*y(k,110)
         mat(k,1469) = rxt(k,137)*y(k,50) + (rxt(k,187)+rxt(k,259))*y(k,83) + ( &
                      + 2.000_r8*rxt(k,141)+2.000_r8*rxt(k,264))*y(k,87)
         mat(k,27) = rxt(k,138)*y(k,110)
         mat(k,1555) = rxt(k,140)*y(k,50) + (rxt(k,166)+rxt(k,261))*y(k,83) + ( &
                      + 2.000_r8*rxt(k,162)+2.000_r8*rxt(k,265))*y(k,87)
         mat(k,655) = mat(k,655) + rxt(k,293)*y(k,40) + rxt(k,183)*y(k,73) + ( &
                      + rxt(k,187)+rxt(k,259))*y(k,75) + (rxt(k,166)+rxt(k,261)) &
                      *y(k,77)
         mat(k,479) = rxt(k,288)*y(k,50)
         mat(k,520) = mat(k,520) + rxt(k,292)*y(k,50)
         mat(k,439) = 2.000_r8*rxt(k,180)*y(k,73) + (2.000_r8*rxt(k,141) &
                       +2.000_r8*rxt(k,264))*y(k,75) + (2.000_r8*rxt(k,162) &
                       +2.000_r8*rxt(k,265))*y(k,77) + rxt(k,233)*y(k,88) + rxt(k,215) &
                      *y(k,89) + rxt(k,197)*y(k,92) + rxt(k,250)*y(k,99)
         mat(k,877) = rxt(k,233)*y(k,87)
         mat(k,920) = rxt(k,215)*y(k,87)
         mat(k,1047) = rxt(k,197)*y(k,87)
         mat(k,1348) = rxt(k,250)*y(k,87)
         mat(k,757) = mat(k,757) + rxt(k,441)*y(k,52)
         mat(k,1801) = rxt(k,135)*y(k,74) + rxt(k,138)*y(k,76)
         mat(k,680) = -(rxt(k,305)*y(k,92) + (rxt(k,403) + rxt(k,404) + rxt(k,405) &
                      ) * y(k,32) + rxt(k,407)*y(k,54) + rxt(k,408)*y(k,56) + rxt(k,412) &
                      *y(k,104) + 4._r8*rxt(k,417)*y(k,43) + rxt(k,429)*y(k,53) &
                      + rxt(k,434)*y(k,51) + rxt(k,439)*y(k,52) + (rxt(k,449) &
                      + rxt(k,450)) * y(k,60) + rxt(k,456)*y(k,24) + rxt(k,482) &
                      *y(k,59) + rxt(k,488)*y(k,2) + rxt(k,528)*y(k,18))
         mat(k,1035) = -rxt(k,305)*y(k,43)
         mat(k,1691) = -(rxt(k,403) + rxt(k,404) + rxt(k,405)) * y(k,43)
         mat(k,1217) = -rxt(k,407)*y(k,43)
         mat(k,1255) = -rxt(k,408)*y(k,43)
         mat(k,748) = -rxt(k,412)*y(k,43)
         mat(k,1415) = -rxt(k,429)*y(k,43)
         mat(k,993) = -rxt(k,434)*y(k,43)
         mat(k,1736) = -rxt(k,439)*y(k,43)
         mat(k,1379) = -(rxt(k,449) + rxt(k,450)) * y(k,43)
         mat(k,1500) = -rxt(k,456)*y(k,43)
         mat(k,326) = -rxt(k,482)*y(k,43)
         mat(k,454) = -rxt(k,488)*y(k,43)
         mat(k,254) = -rxt(k,528)*y(k,43)
         mat(k,454) = mat(k,454) + rxt(k,493)*y(k,104)
         mat(k,537) = rxt(k,525)*y(k,53) + rxt(k,526)*y(k,54) + rxt(k,481)*y(k,59) &
                      + rxt(k,445)*y(k,60)
         mat(k,254) = mat(k,254) + rxt(k,452)*y(k,24) + rxt(k,529)*y(k,51)
         mat(k,1500) = mat(k,1500) + rxt(k,452)*y(k,18) + rxt(k,463)*y(k,104)
         mat(k,143) = rxt(k,532)*y(k,104)
         mat(k,1691) = mat(k,1691) + rxt(k,406)*y(k,55) + rxt(k,312)*y(k,98)
         mat(k,85) = rxt(k,402)*y(k,54) + rxt(k,448)*y(k,60) + rxt(k,411)*y(k,104)
         mat(k,1121) = rxt(k,125)*y(k,68) + rxt(k,313)*y(k,98)
         mat(k,1161) = rxt(k,314)*y(k,98)
         mat(k,993) = mat(k,993) + rxt(k,529)*y(k,18)
         mat(k,1415) = mat(k,1415) + rxt(k,525)*y(k,14) + rxt(k,432)*y(k,104)
         mat(k,1217) = mat(k,1217) + rxt(k,526)*y(k,14) + rxt(k,402)*y(k,35) &
                      + rxt(k,342)*y(k,105)
         mat(k,1655) = rxt(k,406)*y(k,32)
         mat(k,1255) = mat(k,1255) + rxt(k,414)*y(k,104)
         mat(k,326) = mat(k,326) + rxt(k,481)*y(k,14)
         mat(k,1379) = mat(k,1379) + rxt(k,445)*y(k,14) + rxt(k,448)*y(k,35)
         mat(k,615) = rxt(k,125)*y(k,40)
         mat(k,1291) = rxt(k,312)*y(k,32) + rxt(k,313)*y(k,40) + rxt(k,314)*y(k,42)
         mat(k,748) = mat(k,748) + rxt(k,493)*y(k,2) + rxt(k,463)*y(k,24) + rxt(k,532) &
                      *y(k,27) + rxt(k,411)*y(k,35) + rxt(k,432)*y(k,53) + rxt(k,414) &
                      *y(k,56)
         mat(k,1576) = rxt(k,342)*y(k,54)
         mat(k,51) = -(rxt(k,418)*y(k,104))
         mat(k,732) = -rxt(k,418)*y(k,44)
         mat(k,668) = rxt(k,439)*y(k,52)
         mat(k,1719) = rxt(k,439)*y(k,43)
         mat(k,165) = -(rxt(k,496)*y(k,54) + (rxt(k,549) + rxt(k,554)) * y(k,40))
         mat(k,1199) = -rxt(k,496)*y(k,45)
         mat(k,1107) = -(rxt(k,549) + rxt(k,554)) * y(k,45)
         mat(k,449) = rxt(k,488)*y(k,43)
         mat(k,671) = rxt(k,488)*y(k,2)
         mat(k,191) = -(rxt(k,467)*y(k,60) + rxt(k,468)*y(k,54) + rxt(k,469)*y(k,104) &
                      + (rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,40))
         mat(k,1367) = -rxt(k,467)*y(k,46)
         mat(k,1201) = -rxt(k,468)*y(k,46)
         mat(k,739) = -rxt(k,469)*y(k,46)
         mat(k,1108) = -(rxt(k,542) + rxt(k,548) + rxt(k,553)) * y(k,46)
         mat(k,1489) = rxt(k,456)*y(k,43)
         mat(k,210) = rxt(k,461)*y(k,104)
         mat(k,673) = rxt(k,456)*y(k,24)
         mat(k,739) = mat(k,739) + rxt(k,461)*y(k,25)
         mat(k,137) = -(rxt(k,335)*y(k,104))
         mat(k,734) = -rxt(k,335)*y(k,47)
         mat(k,1106) = rxt(k,284)*y(k,81)
         mat(k,1150) = rxt(k,285)*y(k,81)
         mat(k,979) = rxt(k,344)*y(k,104)
         mat(k,702) = rxt(k,284)*y(k,40) + rxt(k,285)*y(k,42)
         mat(k,46) = rxt(k,300)*y(k,110)
         mat(k,734) = mat(k,734) + rxt(k,344)*y(k,51)
         mat(k,1772) = rxt(k,300)*y(k,90)
         mat(k,265) = -(rxt(k,421)*y(k,51) + (rxt(k,422) + rxt(k,423) + rxt(k,424) &
                      ) * y(k,52) + rxt(k,425)*y(k,55) + rxt(k,433)*y(k,104) + rxt(k,570) &
                      *y(k,99))
         mat(k,982) = -rxt(k,421)*y(k,48)
         mat(k,1724) = -(rxt(k,422) + rxt(k,423) + rxt(k,424)) * y(k,48)
         mat(k,1649) = -rxt(k,425)*y(k,48)
         mat(k,742) = -rxt(k,433)*y(k,48)
         mat(k,1320) = -rxt(k,570)*y(k,48)
         mat(k,1206) = rxt(k,419)*y(k,79) + rxt(k,567)*y(k,94)
         mat(k,1649) = mat(k,1649) + rxt(k,568)*y(k,94)
         mat(k,950) = 1.100_r8*rxt(k,563)*y(k,80) + .200_r8*rxt(k,561)*y(k,88)
         mat(k,72) = rxt(k,419)*y(k,54)
         mat(k,114) = 1.100_r8*rxt(k,563)*y(k,69)
         mat(k,850) = .200_r8*rxt(k,561)*y(k,69)
         mat(k,93) = rxt(k,567)*y(k,54) + rxt(k,568)*y(k,55)
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
         mat(k,239) = -(rxt(k,137)*y(k,75) + rxt(k,140)*y(k,77) + rxt(k,288)*y(k,84) &
                      + rxt(k,292)*y(k,85))
         mat(k,1443) = -rxt(k,137)*y(k,50)
         mat(k,1529) = -rxt(k,140)*y(k,50)
         mat(k,469) = -rxt(k,288)*y(k,50)
         mat(k,510) = -rxt(k,292)*y(k,50)
         mat(k,1723) = rxt(k,440)*y(k,53)
         mat(k,1408) = rxt(k,440)*y(k,52)
         mat(k,1001) = -((rxt(k,106) + rxt(k,107)) * y(k,64) + (rxt(k,117) + rxt(k,118) &
                      ) * y(k,67) + rxt(k,131)*y(k,106) + (rxt(k,267) + rxt(k,274) &
                      ) * y(k,101) + rxt(k,275)*y(k,65) + rxt(k,344)*y(k,104) &
                      + rxt(k,421)*y(k,48) + rxt(k,430)*y(k,53) + rxt(k,434)*y(k,43) &
                      + rxt(k,435)*y(k,56) + rxt(k,436)*y(k,54) + rxt(k,457)*y(k,24) &
                      + rxt(k,489)*y(k,2) + rxt(k,529)*y(k,18) + rxt(k,572)*y(k,99))
         mat(k,224) = -(rxt(k,106) + rxt(k,107)) * y(k,51)
         mat(k,499) = -(rxt(k,117) + rxt(k,118)) * y(k,51)
         mat(k,1620) = -rxt(k,131)*y(k,51)
         mat(k,400) = -(rxt(k,267) + rxt(k,274)) * y(k,51)
         mat(k,784) = -rxt(k,275)*y(k,51)
         mat(k,753) = -rxt(k,344)*y(k,51)
         mat(k,270) = -rxt(k,421)*y(k,51)
         mat(k,1422) = -rxt(k,430)*y(k,51)
         mat(k,684) = -rxt(k,434)*y(k,51)
         mat(k,1263) = -rxt(k,435)*y(k,51)
         mat(k,1225) = -rxt(k,436)*y(k,51)
         mat(k,1508) = -rxt(k,457)*y(k,51)
         mat(k,456) = -rxt(k,489)*y(k,51)
         mat(k,256) = -rxt(k,529)*y(k,51)
         mat(k,1344) = -rxt(k,572)*y(k,51)
         mat(k,1699) = rxt(k,282)*y(k,81) + rxt(k,304)*y(k,92)
         mat(k,270) = mat(k,270) + 2.000_r8*rxt(k,423)*y(k,52) + rxt(k,425)*y(k,55) &
                      + rxt(k,433)*y(k,104)
         mat(k,1744) = 2.000_r8*rxt(k,423)*y(k,48) + rxt(k,426)*y(k,54) + rxt(k,286) &
                      *y(k,81)
         mat(k,1225) = mat(k,1225) + rxt(k,426)*y(k,52)
         mat(k,1663) = rxt(k,425)*y(k,48) + rxt(k,420)*y(k,79)
         mat(k,598) = rxt(k,243)*y(k,88) + rxt(k,225)*y(k,89) + rxt(k,207)*y(k,92)
         mat(k,355) = rxt(k,247)*y(k,88) + rxt(k,229)*y(k,89) + rxt(k,212)*y(k,92)
         mat(k,339) = rxt(k,235)*y(k,88) + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92)
         mat(k,784) = mat(k,784) + rxt(k,234)*y(k,88) + rxt(k,216)*y(k,89) &
                      + rxt(k,199)*y(k,92)
         mat(k,310) = rxt(k,241)*y(k,88) + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92)
         mat(k,499) = mat(k,499) + rxt(k,244)*y(k,88) + rxt(k,226)*y(k,89) &
                      + rxt(k,208)*y(k,92)
         mat(k,623) = rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) + rxt(k,204)*y(k,92)
         mat(k,959) = rxt(k,298)*y(k,89) + rxt(k,299)*y(k,90) + rxt(k,301)*y(k,91) &
                      + rxt(k,303)*y(k,92) + rxt(k,379)*y(k,93)
         mat(k,296) = rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92)
         mat(k,74) = rxt(k,420)*y(k,55)
         mat(k,712) = rxt(k,282)*y(k,32) + rxt(k,286)*y(k,52) + rxt(k,249)*y(k,88) &
                      + rxt(k,232)*y(k,89) + rxt(k,214)*y(k,92)
         mat(k,370) = rxt(k,245)*y(k,88) + rxt(k,227)*y(k,89) + rxt(k,210)*y(k,92)
         mat(k,652) = rxt(k,236)*y(k,88) + rxt(k,218)*y(k,89) + rxt(k,201)*y(k,92)
         mat(k,477) = rxt(k,248)*y(k,88) + rxt(k,230)*y(k,89) + rxt(k,213)*y(k,92)
         mat(k,518) = rxt(k,240)*y(k,88) + rxt(k,223)*y(k,89) + rxt(k,205)*y(k,92)
         mat(k,384) = rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) + rxt(k,211)*y(k,92)
         mat(k,437) = rxt(k,233)*y(k,88) + rxt(k,215)*y(k,89) + rxt(k,197)*y(k,92)
         mat(k,873) = rxt(k,243)*y(k,61) + rxt(k,247)*y(k,62) + rxt(k,235)*y(k,63) &
                      + rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) + rxt(k,244)*y(k,67) &
                      + rxt(k,239)*y(k,68) + rxt(k,237)*y(k,71) + rxt(k,249)*y(k,81) &
                      + rxt(k,245)*y(k,82) + rxt(k,236)*y(k,83) + rxt(k,248)*y(k,84) &
                      + rxt(k,240)*y(k,85) + rxt(k,246)*y(k,86) + rxt(k,233)*y(k,87) &
                      + rxt(k,238)*y(k,98)
         mat(k,916) = rxt(k,225)*y(k,61) + rxt(k,229)*y(k,62) + rxt(k,217)*y(k,63) &
                      + rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67) &
                      + rxt(k,222)*y(k,68) + rxt(k,298)*y(k,69) + rxt(k,219)*y(k,71) &
                      + rxt(k,232)*y(k,81) + rxt(k,227)*y(k,82) + rxt(k,218)*y(k,83) &
                      + rxt(k,230)*y(k,84) + rxt(k,223)*y(k,85) + rxt(k,228)*y(k,86) &
                      + rxt(k,215)*y(k,87) + rxt(k,221)*y(k,98)
         mat(k,49) = rxt(k,299)*y(k,69)
         mat(k,133) = rxt(k,301)*y(k,69)
         mat(k,1043) = rxt(k,304)*y(k,32) + rxt(k,207)*y(k,61) + rxt(k,212)*y(k,62) &
                      + rxt(k,200)*y(k,63) + rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) &
                      + rxt(k,208)*y(k,67) + rxt(k,204)*y(k,68) + rxt(k,303)*y(k,69) &
                      + rxt(k,202)*y(k,71) + rxt(k,214)*y(k,81) + rxt(k,210)*y(k,82) &
                      + rxt(k,201)*y(k,83) + rxt(k,213)*y(k,84) + rxt(k,205)*y(k,85) &
                      + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87) + rxt(k,203)*y(k,98)
         mat(k,126) = rxt(k,379)*y(k,69)
         mat(k,1299) = rxt(k,238)*y(k,88) + rxt(k,221)*y(k,89) + rxt(k,203)*y(k,92)
         mat(k,753) = mat(k,753) + rxt(k,433)*y(k,48)
         mat(k,1762) = -(rxt(k,105)*y(k,61) + (rxt(k,115) + rxt(k,116)) * y(k,67) &
                      + (rxt(k,272) + rxt(k,273)) * y(k,101) + rxt(k,276)*y(k,65) &
                      + rxt(k,286)*y(k,81) + rxt(k,315)*y(k,98) + rxt(k,341)*y(k,105) &
                      + rxt(k,354)*y(k,106) + (rxt(k,422) + rxt(k,423) + rxt(k,424) &
                      ) * y(k,48) + (rxt(k,426) + rxt(k,428)) * y(k,54) + rxt(k,427) &
                      *y(k,56) + rxt(k,439)*y(k,43) + rxt(k,440)*y(k,53) + rxt(k,441) &
                      *y(k,104) + rxt(k,459)*y(k,24) + rxt(k,490)*y(k,2))
         mat(k,608) = -rxt(k,105)*y(k,52)
         mat(k,508) = -(rxt(k,115) + rxt(k,116)) * y(k,52)
         mat(k,409) = -(rxt(k,272) + rxt(k,273)) * y(k,52)
         mat(k,802) = -rxt(k,276)*y(k,52)
         mat(k,727) = -rxt(k,286)*y(k,52)
         mat(k,1317) = -rxt(k,315)*y(k,52)
         mat(k,1602) = -rxt(k,341)*y(k,52)
         mat(k,1638) = -rxt(k,354)*y(k,52)
         mat(k,276) = -(rxt(k,422) + rxt(k,423) + rxt(k,424)) * y(k,52)
         mat(k,1243) = -(rxt(k,426) + rxt(k,428)) * y(k,52)
         mat(k,1281) = -rxt(k,427)*y(k,52)
         mat(k,700) = -rxt(k,439)*y(k,52)
         mat(k,1440) = -rxt(k,440)*y(k,52)
         mat(k,769) = -rxt(k,441)*y(k,52)
         mat(k,1526) = -rxt(k,459)*y(k,52)
         mat(k,467) = -rxt(k,490)*y(k,52)
         mat(k,467) = mat(k,467) + rxt(k,489)*y(k,51)
         mat(k,261) = rxt(k,529)*y(k,51)
         mat(k,1526) = mat(k,1526) + rxt(k,457)*y(k,51)
         mat(k,700) = mat(k,700) + rxt(k,434)*y(k,51) + rxt(k,429)*y(k,53)
         mat(k,56) = rxt(k,418)*y(k,104)
         mat(k,140) = rxt(k,335)*y(k,104)
         mat(k,1019) = rxt(k,489)*y(k,2) + rxt(k,529)*y(k,18) + rxt(k,457)*y(k,24) &
                      + rxt(k,434)*y(k,43) + 2.000_r8*rxt(k,430)*y(k,53) + rxt(k,436) &
                      *y(k,54) + rxt(k,435)*y(k,56) + rxt(k,107)*y(k,64) + rxt(k,131) &
                      *y(k,106)
         mat(k,1440) = mat(k,1440) + rxt(k,429)*y(k,43) + 2.000_r8*rxt(k,430)*y(k,51) &
                      + rxt(k,431)*y(k,54) + rxt(k,432)*y(k,104)
         mat(k,1243) = mat(k,1243) + rxt(k,436)*y(k,51) + rxt(k,431)*y(k,53)
         mat(k,1281) = mat(k,1281) + rxt(k,435)*y(k,51)
         mat(k,1405) = rxt(k,280)*y(k,81)
         mat(k,228) = rxt(k,107)*y(k,51)
         mat(k,845) = rxt(k,196)*y(k,81) + rxt(k,192)*y(k,82)
         mat(k,1483) = rxt(k,161)*y(k,81) + rxt(k,157)*y(k,82)
         mat(k,1569) = rxt(k,179)*y(k,81) + rxt(k,174)*y(k,82)
         mat(k,727) = mat(k,727) + rxt(k,280)*y(k,60) + rxt(k,196)*y(k,73) &
                      + rxt(k,161)*y(k,75) + rxt(k,179)*y(k,77) + rxt(k,249)*y(k,88) &
                      + rxt(k,232)*y(k,89) + rxt(k,214)*y(k,92) + rxt(k,156)*y(k,99)
         mat(k,377) = rxt(k,192)*y(k,73) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77) &
                      + rxt(k,245)*y(k,88) + rxt(k,227)*y(k,89) + rxt(k,210)*y(k,92) &
                      + rxt(k,151)*y(k,99)
         mat(k,891) = rxt(k,249)*y(k,81) + rxt(k,245)*y(k,82)
         mat(k,934) = rxt(k,232)*y(k,81) + rxt(k,227)*y(k,82)
         mat(k,1061) = rxt(k,214)*y(k,81) + rxt(k,210)*y(k,82) + rxt(k,306)*y(k,104)
         mat(k,1362) = rxt(k,156)*y(k,81) + rxt(k,151)*y(k,82)
         mat(k,769) = mat(k,769) + rxt(k,418)*y(k,44) + rxt(k,335)*y(k,47) &
                      + rxt(k,432)*y(k,53) + rxt(k,306)*y(k,92)
         mat(k,1638) = mat(k,1638) + rxt(k,131)*y(k,51)
         mat(k,1432) = -(rxt(k,429)*y(k,43) + rxt(k,430)*y(k,51) + rxt(k,431)*y(k,54) &
                      + rxt(k,432)*y(k,104) + rxt(k,440)*y(k,52) + rxt(k,525)*y(k,14))
         mat(k,694) = -rxt(k,429)*y(k,53)
         mat(k,1011) = -rxt(k,430)*y(k,53)
         mat(k,1235) = -rxt(k,431)*y(k,53)
         mat(k,763) = -rxt(k,432)*y(k,53)
         mat(k,1754) = -rxt(k,440)*y(k,53)
         mat(k,547) = -rxt(k,525)*y(k,53)
         mat(k,82) = rxt(k,491)*y(k,54)
         mat(k,1518) = rxt(k,281)*y(k,81)
         mat(k,217) = rxt(k,460)*y(k,54) + rxt(k,458)*y(k,60) + rxt(k,461)*y(k,104)
         mat(k,179) = rxt(k,505)*y(k,42)
         mat(k,1179) = rxt(k,505)*y(k,31) + rxt(k,442)*y(k,104)
         mat(k,694) = mat(k,694) + rxt(k,305)*y(k,92)
         mat(k,1754) = mat(k,1754) + rxt(k,428)*y(k,54) + rxt(k,427)*y(k,56)
         mat(k,1235) = mat(k,1235) + rxt(k,491)*y(k,3) + rxt(k,460)*y(k,25) &
                      + rxt(k,428)*y(k,52)
         mat(k,1273) = rxt(k,427)*y(k,52)
         mat(k,1397) = rxt(k,458)*y(k,25)
         mat(k,837) = rxt(k,195)*y(k,84) + rxt(k,188)*y(k,85) + rxt(k,193)*y(k,86)
         mat(k,1475) = rxt(k,160)*y(k,84) + (rxt(k,231)+rxt(k,257))*y(k,85) &
                      + rxt(k,158)*y(k,86)
         mat(k,1561) = rxt(k,178)*y(k,84) + (rxt(k,170)+rxt(k,258))*y(k,85) &
                      + rxt(k,175)*y(k,86)
         mat(k,720) = rxt(k,281)*y(k,24)
         mat(k,660) = rxt(k,236)*y(k,88) + rxt(k,218)*y(k,89) + rxt(k,201)*y(k,92) &
                      + rxt(k,143)*y(k,99)
         mat(k,482) = rxt(k,195)*y(k,73) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77) &
                      + rxt(k,248)*y(k,88) + rxt(k,230)*y(k,89) + rxt(k,213)*y(k,92) &
                      + rxt(k,155)*y(k,99)
         mat(k,523) = rxt(k,188)*y(k,73) + (rxt(k,231)+rxt(k,257))*y(k,75) + ( &
                      + rxt(k,170)+rxt(k,258))*y(k,77) + rxt(k,240)*y(k,88) &
                      + rxt(k,223)*y(k,89) + rxt(k,205)*y(k,92) + rxt(k,147)*y(k,99)
         mat(k,389) = rxt(k,193)*y(k,73) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77) &
                      + rxt(k,246)*y(k,88) + rxt(k,228)*y(k,89) + rxt(k,211)*y(k,92) &
                      + rxt(k,152)*y(k,99)
         mat(k,441) = rxt(k,233)*y(k,88) + rxt(k,215)*y(k,89) + rxt(k,197)*y(k,92) &
                      + rxt(k,250)*y(k,99)
         mat(k,883) = rxt(k,236)*y(k,83) + rxt(k,248)*y(k,84) + rxt(k,240)*y(k,85) &
                      + rxt(k,246)*y(k,86) + rxt(k,233)*y(k,87)
         mat(k,926) = rxt(k,218)*y(k,83) + rxt(k,230)*y(k,84) + rxt(k,223)*y(k,85) &
                      + rxt(k,228)*y(k,86) + rxt(k,215)*y(k,87)
         mat(k,1053) = rxt(k,305)*y(k,43) + rxt(k,201)*y(k,83) + rxt(k,213)*y(k,84) &
                      + rxt(k,205)*y(k,85) + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87)
         mat(k,1354) = rxt(k,143)*y(k,83) + rxt(k,155)*y(k,84) + rxt(k,147)*y(k,85) &
                      + rxt(k,152)*y(k,86) + rxt(k,250)*y(k,87)
         mat(k,763) = mat(k,763) + rxt(k,461)*y(k,25) + rxt(k,442)*y(k,42)
         mat(k,1230) = -(rxt(k,108)*y(k,64) + rxt(k,120)*y(k,65) + rxt(k,126)*y(k,68) &
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
         mat(k,225) = -rxt(k,108)*y(k,54)
         mat(k,789) = -rxt(k,120)*y(k,54)
         mat(k,628) = -rxt(k,126)*y(k,54)
         mat(k,656) = -rxt(k,296)*y(k,54)
         mat(k,1304) = -(rxt(k,319) + rxt(k,320)) * y(k,54)
         mat(k,402) = -(rxt(k,328) + rxt(k,329)) * y(k,54)
         mat(k,102) = -rxt(k,331)*y(k,54)
         mat(k,283) = -rxt(k,333)*y(k,54)
         mat(k,1589) = -rxt(k,342)*y(k,54)
         mat(k,1625) = -rxt(k,355)*y(k,54)
         mat(k,1268) = -rxt(k,398)*y(k,54)
         mat(k,1668) = -rxt(k,400)*y(k,54)
         mat(k,422) = -rxt(k,401)*y(k,54)
         mat(k,88) = -rxt(k,402)*y(k,54)
         mat(k,689) = -rxt(k,407)*y(k,54)
         mat(k,758) = -rxt(k,413)*y(k,54)
         mat(k,1749) = -(rxt(k,426) + rxt(k,428)) * y(k,54)
         mat(k,1427) = -rxt(k,431)*y(k,54)
         mat(k,1006) = -rxt(k,436)*y(k,54)
         mat(k,215) = -rxt(k,460)*y(k,54)
         mat(k,1513) = -rxt(k,462)*y(k,54)
         mat(k,1134) = -rxt(k,465)*y(k,54)
         mat(k,194) = -rxt(k,468)*y(k,54)
         mat(k,81) = -rxt(k,491)*y(k,54)
         mat(k,460) = -rxt(k,492)*y(k,54)
         mat(k,187) = -rxt(k,494)*y(k,54)
         mat(k,170) = -rxt(k,496)*y(k,54)
         mat(k,544) = -rxt(k,526)*y(k,54)
         mat(k,119) = -(rxt(k,565) + rxt(k,566)) * y(k,54)
         mat(k,95) = -rxt(k,567)*y(k,54)
         mat(k,1704) = rxt(k,405)*y(k,43)
         mat(k,689) = mat(k,689) + rxt(k,405)*y(k,32)
         mat(k,272) = rxt(k,421)*y(k,51) + rxt(k,422)*y(k,52) + rxt(k,425)*y(k,55) &
                      + rxt(k,570)*y(k,99)
         mat(k,1006) = mat(k,1006) + rxt(k,421)*y(k,48) + rxt(k,267)*y(k,101)
         mat(k,1749) = mat(k,1749) + rxt(k,422)*y(k,48) + rxt(k,354)*y(k,106)
         mat(k,1668) = mat(k,1668) + rxt(k,425)*y(k,48) + rxt(k,569)*y(k,94) + ( &
                      + rxt(k,387)+rxt(k,388))*y(k,95) + rxt(k,575)*y(k,107) &
                      + rxt(k,579)*y(k,108)
         mat(k,1268) = mat(k,1268) + rxt(k,358)*y(k,106)
         mat(k,1392) = rxt(k,109)*y(k,65) + rxt(k,345)*y(k,106)
         mat(k,789) = mat(k,789) + rxt(k,109)*y(k,60) + rxt(k,181)*y(k,73) + ( &
                      + rxt(k,165)+rxt(k,253))*y(k,75) + (rxt(k,163)+rxt(k,260)) &
                      *y(k,77) + rxt(k,234)*y(k,88) + rxt(k,216)*y(k,89) + rxt(k,199) &
                      *y(k,92) + rxt(k,251)*y(k,99)
         mat(k,313) = rxt(k,189)*y(k,73) + (rxt(k,242)+rxt(k,266))*y(k,75) + ( &
                      + rxt(k,171)+rxt(k,254))*y(k,77) + rxt(k,241)*y(k,88) &
                      + rxt(k,224)*y(k,89) + rxt(k,206)*y(k,92) + rxt(k,148)*y(k,99)
         mat(k,502) = rxt(k,191)*y(k,73) + (rxt(k,153)+rxt(k,255))*y(k,75) + ( &
                      + rxt(k,173)+rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + rxt(k,226)*y(k,89) + rxt(k,208)*y(k,92) + rxt(k,150)*y(k,99)
         mat(k,964) = rxt(k,561)*y(k,88) + 1.150_r8*rxt(k,562)*y(k,99)
         mat(k,832) = rxt(k,181)*y(k,65) + rxt(k,189)*y(k,66) + rxt(k,191)*y(k,67)
         mat(k,1470) = (rxt(k,165)+rxt(k,253))*y(k,65) + (rxt(k,242)+rxt(k,266)) &
                      *y(k,66) + (rxt(k,153)+rxt(k,255))*y(k,67)
         mat(k,1556) = (rxt(k,163)+rxt(k,260))*y(k,65) + (rxt(k,171)+rxt(k,254)) &
                      *y(k,66) + (rxt(k,173)+rxt(k,256))*y(k,67)
         mat(k,878) = rxt(k,234)*y(k,65) + rxt(k,241)*y(k,66) + rxt(k,244)*y(k,67) &
                      + rxt(k,561)*y(k,69)
         mat(k,921) = rxt(k,216)*y(k,65) + rxt(k,224)*y(k,66) + rxt(k,226)*y(k,67)
         mat(k,1048) = rxt(k,199)*y(k,65) + rxt(k,206)*y(k,66) + rxt(k,208)*y(k,67)
         mat(k,95) = mat(k,95) + rxt(k,569)*y(k,55)
         mat(k,576) = (rxt(k,387)+rxt(k,388))*y(k,55)
         mat(k,1349) = rxt(k,570)*y(k,48) + rxt(k,251)*y(k,65) + rxt(k,148)*y(k,66) &
                      + rxt(k,150)*y(k,67) + 1.150_r8*rxt(k,562)*y(k,69)
         mat(k,402) = mat(k,402) + rxt(k,267)*y(k,51)
         mat(k,758) = mat(k,758) + 2.000_r8*rxt(k,415)*y(k,104)
         mat(k,1625) = mat(k,1625) + rxt(k,354)*y(k,52) + rxt(k,358)*y(k,56) &
                      + rxt(k,345)*y(k,60)
         mat(k,205) = rxt(k,575)*y(k,55)
         mat(k,68) = rxt(k,579)*y(k,55)
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
         mat(k,1679) = -(rxt(k,121)*y(k,65) + (rxt(k,128) + rxt(k,130)) * y(k,69) &
                      + rxt(k,317)*y(k,98) + rxt(k,357)*y(k,106) + rxt(k,359)*y(k,99) &
                      + rxt(k,387)*y(k,95) + rxt(k,392)*y(k,96) + rxt(k,400)*y(k,54) &
                      + rxt(k,406)*y(k,32) + rxt(k,420)*y(k,79) + rxt(k,425)*y(k,48) &
                      + rxt(k,564)*y(k,80) + (rxt(k,568) + rxt(k,569)) * y(k,94) &
                      + rxt(k,575)*y(k,107) + rxt(k,579)*y(k,108))
         mat(k,800) = -rxt(k,121)*y(k,55)
         mat(k,975) = -(rxt(k,128) + rxt(k,130)) * y(k,55)
         mat(k,1315) = -rxt(k,317)*y(k,55)
         mat(k,1636) = -rxt(k,357)*y(k,55)
         mat(k,1360) = -rxt(k,359)*y(k,55)
         mat(k,585) = -rxt(k,387)*y(k,55)
         mat(k,236) = -rxt(k,392)*y(k,55)
         mat(k,1241) = -rxt(k,400)*y(k,55)
         mat(k,1715) = -rxt(k,406)*y(k,55)
         mat(k,75) = -rxt(k,420)*y(k,55)
         mat(k,274) = -rxt(k,425)*y(k,55)
         mat(k,121) = -rxt(k,564)*y(k,55)
         mat(k,97) = -(rxt(k,568) + rxt(k,569)) * y(k,55)
         mat(k,207) = -rxt(k,575)*y(k,55)
         mat(k,70) = -rxt(k,579)*y(k,55)
         mat(k,465) = 2.000_r8*rxt(k,484)*y(k,2) + (rxt(k,486)+rxt(k,487))*y(k,24) &
                      + rxt(k,488)*y(k,43) + rxt(k,492)*y(k,54)
         mat(k,259) = rxt(k,528)*y(k,43)
         mat(k,1524) = (rxt(k,486)+rxt(k,487))*y(k,2) + (2.000_r8*rxt(k,453) &
                       +2.000_r8*rxt(k,454))*y(k,24) + rxt(k,456)*y(k,43) + rxt(k,462) &
                      *y(k,54) + rxt(k,111)*y(k,65) + rxt(k,123)*y(k,68) + rxt(k,310) &
                      *y(k,98) + rxt(k,464)*y(k,104) + rxt(k,346)*y(k,106)
         mat(k,1100) = rxt(k,325)*y(k,101) + rxt(k,330)*y(k,102)
         mat(k,1715) = mat(k,1715) + rxt(k,403)*y(k,43) + rxt(k,409)*y(k,56) &
                      + rxt(k,326)*y(k,101)
         mat(k,698) = rxt(k,488)*y(k,2) + rxt(k,528)*y(k,18) + rxt(k,456)*y(k,24) &
                      + rxt(k,403)*y(k,32) + 2.000_r8*rxt(k,417)*y(k,43) + rxt(k,429) &
                      *y(k,53) + rxt(k,407)*y(k,54) + 2.000_r8*rxt(k,408)*y(k,56) &
                      + rxt(k,482)*y(k,59) + rxt(k,449)*y(k,60) + rxt(k,412)*y(k,104)
         mat(k,55) = rxt(k,418)*y(k,104)
         mat(k,274) = mat(k,274) + rxt(k,424)*y(k,52)
         mat(k,1017) = rxt(k,435)*y(k,56) + rxt(k,572)*y(k,99) + rxt(k,274)*y(k,101)
         mat(k,1760) = rxt(k,424)*y(k,48) + rxt(k,426)*y(k,54) + rxt(k,427)*y(k,56) &
                      + rxt(k,315)*y(k,98) + rxt(k,272)*y(k,101)
         mat(k,1438) = rxt(k,429)*y(k,43) + rxt(k,431)*y(k,54)
         mat(k,1241) = mat(k,1241) + rxt(k,492)*y(k,2) + rxt(k,462)*y(k,24) &
                      + rxt(k,407)*y(k,43) + rxt(k,426)*y(k,52) + rxt(k,431)*y(k,53) &
                      + 2.000_r8*rxt(k,399)*y(k,54) + 2.000_r8*rxt(k,398)*y(k,56) &
                      + rxt(k,108)*y(k,64) + rxt(k,126)*y(k,68) + rxt(k,296)*y(k,83) &
                      + rxt(k,391)*y(k,96) + rxt(k,320)*y(k,98) + (2.000_r8*rxt(k,328) &
                       +rxt(k,329))*y(k,101) + rxt(k,331)*y(k,102) + rxt(k,413) &
                      *y(k,104) + rxt(k,355)*y(k,106)
         mat(k,1679) = mat(k,1679) + 2.000_r8*rxt(k,392)*y(k,96)
         mat(k,1279) = rxt(k,409)*y(k,32) + 2.000_r8*rxt(k,408)*y(k,43) + rxt(k,435) &
                      *y(k,51) + rxt(k,427)*y(k,52) + 2.000_r8*rxt(k,398)*y(k,54) &
                      + rxt(k,483)*y(k,59) + rxt(k,451)*y(k,60) + rxt(k,127)*y(k,68) &
                      + rxt(k,129)*y(k,69) + rxt(k,287)*y(k,81) + 2.000_r8*rxt(k,297) &
                      *y(k,83) + 2.000_r8*rxt(k,389)*y(k,95) + rxt(k,318)*y(k,98) &
                      + 3.000_r8*rxt(k,327)*y(k,101) + rxt(k,414)*y(k,104)
         mat(k,331) = rxt(k,482)*y(k,43) + rxt(k,483)*y(k,56)
         mat(k,1403) = rxt(k,449)*y(k,43) + rxt(k,451)*y(k,56) + rxt(k,122)*y(k,68) &
                      + rxt(k,309)*y(k,98)
         mat(k,606) = rxt(k,149)*y(k,99)
         mat(k,362) = rxt(k,154)*y(k,99)
         mat(k,346) = rxt(k,252)*y(k,99)
         mat(k,227) = rxt(k,108)*y(k,54)
         mat(k,800) = mat(k,800) + rxt(k,111)*y(k,24) + rxt(k,251)*y(k,99)
         mat(k,317) = rxt(k,148)*y(k,99)
         mat(k,506) = rxt(k,150)*y(k,99)
         mat(k,638) = rxt(k,123)*y(k,24) + rxt(k,126)*y(k,54) + rxt(k,127)*y(k,56) &
                      + rxt(k,122)*y(k,60) + rxt(k,186)*y(k,73) + rxt(k,220)*y(k,75) &
                      + rxt(k,169)*y(k,77) + rxt(k,239)*y(k,88) + rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92) + 2.000_r8*rxt(k,146)*y(k,99)
         mat(k,975) = mat(k,975) + rxt(k,129)*y(k,56) + rxt(k,321)*y(k,100) &
                      + 2.000_r8*rxt(k,375)*y(k,103)
         mat(k,302) = rxt(k,144)*y(k,99)
         mat(k,843) = rxt(k,186)*y(k,68) + rxt(k,185)*y(k,98)
         mat(k,1481) = rxt(k,220)*y(k,68) + rxt(k,209)*y(k,98)
         mat(k,1567) = rxt(k,169)*y(k,68) + rxt(k,168)*y(k,98)
         mat(k,725) = rxt(k,287)*y(k,56) + rxt(k,156)*y(k,99)
         mat(k,375) = rxt(k,151)*y(k,99)
         mat(k,663) = rxt(k,296)*y(k,54) + 2.000_r8*rxt(k,297)*y(k,56) + rxt(k,143) &
                      *y(k,99)
         mat(k,485) = rxt(k,155)*y(k,99)
         mat(k,526) = rxt(k,147)*y(k,99)
         mat(k,392) = rxt(k,152)*y(k,99)
         mat(k,444) = rxt(k,250)*y(k,99)
         mat(k,889) = rxt(k,239)*y(k,68) + rxt(k,238)*y(k,98)
         mat(k,932) = rxt(k,222)*y(k,68) + rxt(k,221)*y(k,98)
         mat(k,1059) = rxt(k,204)*y(k,68) + rxt(k,203)*y(k,98)
         mat(k,585) = mat(k,585) + 2.000_r8*rxt(k,389)*y(k,56)
         mat(k,236) = mat(k,236) + rxt(k,391)*y(k,54) + 2.000_r8*rxt(k,392)*y(k,55) &
                      + 2.000_r8*rxt(k,316)*y(k,98) + 2.000_r8*rxt(k,334)*y(k,103)
         mat(k,1315) = mat(k,1315) + rxt(k,310)*y(k,24) + rxt(k,315)*y(k,52) &
                      + rxt(k,320)*y(k,54) + rxt(k,318)*y(k,56) + rxt(k,309)*y(k,60) &
                      + rxt(k,185)*y(k,73) + rxt(k,209)*y(k,75) + rxt(k,168)*y(k,77) &
                      + rxt(k,238)*y(k,88) + rxt(k,221)*y(k,89) + rxt(k,203)*y(k,92) &
                      + 2.000_r8*rxt(k,316)*y(k,96)
         mat(k,1360) = mat(k,1360) + rxt(k,572)*y(k,51) + rxt(k,149)*y(k,61) &
                      + rxt(k,154)*y(k,62) + rxt(k,252)*y(k,63) + rxt(k,251)*y(k,65) &
                      + rxt(k,148)*y(k,66) + rxt(k,150)*y(k,67) + 2.000_r8*rxt(k,146) &
                      *y(k,68) + rxt(k,144)*y(k,71) + rxt(k,156)*y(k,81) + rxt(k,151) &
                      *y(k,82) + rxt(k,143)*y(k,83) + rxt(k,155)*y(k,84) + rxt(k,147) &
                      *y(k,85) + rxt(k,152)*y(k,86) + rxt(k,250)*y(k,87)
         mat(k,161) = rxt(k,321)*y(k,69) + (rxt(k,322)+rxt(k,323))*y(k,110)
         mat(k,407) = rxt(k,325)*y(k,28) + rxt(k,326)*y(k,32) + rxt(k,274)*y(k,51) &
                      + rxt(k,272)*y(k,52) + (2.000_r8*rxt(k,328)+rxt(k,329))*y(k,54) &
                      + 3.000_r8*rxt(k,327)*y(k,56)
         mat(k,104) = rxt(k,330)*y(k,28) + rxt(k,331)*y(k,54)
         mat(k,288) = 2.000_r8*rxt(k,375)*y(k,69) + 2.000_r8*rxt(k,334)*y(k,96) &
                      + rxt(k,332)*y(k,110)
         mat(k,767) = rxt(k,464)*y(k,24) + rxt(k,412)*y(k,43) + rxt(k,418)*y(k,44) &
                      + rxt(k,413)*y(k,54) + rxt(k,414)*y(k,56)
         mat(k,1636) = mat(k,1636) + rxt(k,346)*y(k,24) + rxt(k,355)*y(k,54)
         mat(k,1813) = (rxt(k,322)+rxt(k,323))*y(k,100) + rxt(k,332)*y(k,103)
         mat(k,1269) = -(rxt(k,127)*y(k,68) + rxt(k,129)*y(k,69) + rxt(k,287)*y(k,81) &
                      + rxt(k,297)*y(k,83) + rxt(k,318)*y(k,98) + rxt(k,327)*y(k,101) &
                      + rxt(k,343)*y(k,105) + rxt(k,358)*y(k,106) + rxt(k,389)*y(k,95) &
                      + rxt(k,398)*y(k,54) + rxt(k,408)*y(k,43) + rxt(k,409)*y(k,32) &
                      + rxt(k,414)*y(k,104) + rxt(k,427)*y(k,52) + rxt(k,435)*y(k,51) &
                      + rxt(k,451)*y(k,60) + rxt(k,483)*y(k,59))
         mat(k,629) = -rxt(k,127)*y(k,56)
         mat(k,965) = -rxt(k,129)*y(k,56)
         mat(k,717) = -rxt(k,287)*y(k,56)
         mat(k,657) = -rxt(k,297)*y(k,56)
         mat(k,1305) = -rxt(k,318)*y(k,56)
         mat(k,403) = -rxt(k,327)*y(k,56)
         mat(k,1590) = -rxt(k,343)*y(k,56)
         mat(k,1626) = -rxt(k,358)*y(k,56)
         mat(k,577) = -rxt(k,389)*y(k,56)
         mat(k,1231) = -rxt(k,398)*y(k,56)
         mat(k,690) = -rxt(k,408)*y(k,56)
         mat(k,1705) = -rxt(k,409)*y(k,56)
         mat(k,759) = -rxt(k,414)*y(k,56)
         mat(k,1750) = -rxt(k,427)*y(k,56)
         mat(k,1007) = -rxt(k,435)*y(k,56)
         mat(k,1393) = -rxt(k,451)*y(k,56)
         mat(k,330) = -rxt(k,483)*y(k,56)
         mat(k,1750) = mat(k,1750) + rxt(k,273)*y(k,101)
         mat(k,1231) = mat(k,1231) + rxt(k,400)*y(k,55) + rxt(k,319)*y(k,98) &
                      + rxt(k,333)*y(k,103)
         mat(k,1669) = rxt(k,400)*y(k,54)
         mat(k,232) = rxt(k,356)*y(k,106)
         mat(k,1305) = mat(k,1305) + rxt(k,319)*y(k,54)
         mat(k,403) = mat(k,403) + rxt(k,273)*y(k,52)
         mat(k,284) = rxt(k,333)*y(k,54)
         mat(k,1626) = mat(k,1626) + rxt(k,356)*y(k,96)
         mat(k,446) = rxt(k,485)*y(k,24)
         mat(k,1486) = rxt(k,485)*y(k,2) + 2.000_r8*rxt(k,455)*y(k,24)
         mat(k,322) = -(rxt(k,481)*y(k,14) + rxt(k,482)*y(k,43) + rxt(k,483)*y(k,56))
         mat(k,532) = -rxt(k,481)*y(k,59)
         mat(k,675) = -rxt(k,482)*y(k,59)
         mat(k,1248) = -rxt(k,483)*y(k,59)
         mat(k,450) = 4.000_r8*rxt(k,484)*y(k,2) + (rxt(k,485)+rxt(k,486))*y(k,24) &
                      + rxt(k,489)*y(k,51) + rxt(k,492)*y(k,54) + rxt(k,493)*y(k,104)
         mat(k,1493) = (rxt(k,485)+rxt(k,486))*y(k,2)
         mat(k,183) = rxt(k,494)*y(k,54) + rxt(k,500)*y(k,95) + rxt(k,495)*y(k,104)
         mat(k,983) = rxt(k,489)*y(k,2)
         mat(k,1208) = rxt(k,492)*y(k,2) + rxt(k,494)*y(k,36)
         mat(k,563) = rxt(k,500)*y(k,36)
         mat(k,743) = rxt(k,493)*y(k,2) + rxt(k,495)*y(k,36)
         mat(k,1396) = -((rxt(k,109) + rxt(k,110)) * y(k,65) + rxt(k,122)*y(k,68) &
                      + rxt(k,280)*y(k,81) + rxt(k,309)*y(k,98) + rxt(k,336)*y(k,105) &
                      + rxt(k,345)*y(k,106) + rxt(k,445)*y(k,14) + rxt(k,447)*y(k,33) &
                      + rxt(k,448)*y(k,35) + (rxt(k,449) + rxt(k,450)) * y(k,43) &
                      + rxt(k,451)*y(k,56) + rxt(k,458)*y(k,25) + rxt(k,467)*y(k,46))
         mat(k,793) = -(rxt(k,109) + rxt(k,110)) * y(k,60)
         mat(k,632) = -rxt(k,122)*y(k,60)
         mat(k,719) = -rxt(k,280)*y(k,60)
         mat(k,1308) = -rxt(k,309)*y(k,60)
         mat(k,1593) = -rxt(k,336)*y(k,60)
         mat(k,1629) = -rxt(k,345)*y(k,60)
         mat(k,546) = -rxt(k,445)*y(k,60)
         mat(k,423) = -rxt(k,447)*y(k,60)
         mat(k,89) = -rxt(k,448)*y(k,60)
         mat(k,693) = -(rxt(k,449) + rxt(k,450)) * y(k,60)
         mat(k,1272) = -rxt(k,451)*y(k,60)
         mat(k,216) = -rxt(k,458)*y(k,60)
         mat(k,195) = -rxt(k,467)*y(k,60)
         mat(k,462) = rxt(k,486)*y(k,24)
         mat(k,257) = rxt(k,452)*y(k,24)
         mat(k,1517) = rxt(k,486)*y(k,2) + rxt(k,452)*y(k,18) + (4.000_r8*rxt(k,453) &
                       +2.000_r8*rxt(k,455))*y(k,24) + rxt(k,457)*y(k,51) + rxt(k,462) &
                      *y(k,54) + rxt(k,463)*y(k,104)
         mat(k,15) = rxt(k,507)*y(k,95)
         mat(k,1138) = rxt(k,465)*y(k,54) + rxt(k,477)*y(k,95) + rxt(k,466)*y(k,104)
         mat(k,1010) = rxt(k,457)*y(k,24) + rxt(k,106)*y(k,64)
         mat(k,1753) = rxt(k,105)*y(k,61)
         mat(k,1234) = rxt(k,462)*y(k,24) + rxt(k,465)*y(k,40)
         mat(k,603) = rxt(k,105)*y(k,52) + rxt(k,190)*y(k,73) + rxt(k,142)*y(k,75) &
                      + rxt(k,172)*y(k,77) + rxt(k,243)*y(k,88) + rxt(k,225)*y(k,89) &
                      + rxt(k,207)*y(k,92) + rxt(k,149)*y(k,99)
         mat(k,359) = rxt(k,194)*y(k,73) + rxt(k,159)*y(k,75) + rxt(k,177)*y(k,77) &
                      + rxt(k,247)*y(k,88) + rxt(k,229)*y(k,89) + rxt(k,212)*y(k,92) &
                      + rxt(k,154)*y(k,99)
         mat(k,343) = rxt(k,182)*y(k,73) + rxt(k,176)*y(k,75) + rxt(k,164)*y(k,77) &
                      + rxt(k,235)*y(k,88) + rxt(k,217)*y(k,89) + rxt(k,200)*y(k,92) &
                      + rxt(k,252)*y(k,99)
         mat(k,226) = rxt(k,106)*y(k,51)
         mat(k,836) = rxt(k,190)*y(k,61) + rxt(k,194)*y(k,62) + rxt(k,182)*y(k,63)
         mat(k,1474) = rxt(k,142)*y(k,61) + rxt(k,159)*y(k,62) + rxt(k,176)*y(k,63)
         mat(k,1560) = rxt(k,172)*y(k,61) + rxt(k,177)*y(k,62) + rxt(k,164)*y(k,63)
         mat(k,882) = rxt(k,243)*y(k,61) + rxt(k,247)*y(k,62) + rxt(k,235)*y(k,63)
         mat(k,925) = rxt(k,225)*y(k,61) + rxt(k,229)*y(k,62) + rxt(k,217)*y(k,63)
         mat(k,1052) = rxt(k,207)*y(k,61) + rxt(k,212)*y(k,62) + rxt(k,200)*y(k,63)
         mat(k,580) = rxt(k,507)*y(k,30) + rxt(k,477)*y(k,40)
         mat(k,1353) = rxt(k,149)*y(k,61) + rxt(k,154)*y(k,62) + rxt(k,252)*y(k,63)
         mat(k,762) = rxt(k,463)*y(k,24) + rxt(k,466)*y(k,40)
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
         mat(k,591) = -(rxt(k,100)*y(k,32) + rxt(k,102)*y(k,110) + rxt(k,103)*y(k,40) &
                      + rxt(k,104)*y(k,42) + rxt(k,105)*y(k,52) + rxt(k,142)*y(k,75) &
                      + rxt(k,149)*y(k,99) + rxt(k,172)*y(k,77) + rxt(k,190)*y(k,73) &
                      + rxt(k,207)*y(k,92) + rxt(k,225)*y(k,89) + rxt(k,243)*y(k,88))
         mat(k,1688) = -rxt(k,100)*y(k,61)
         mat(k,1787) = -rxt(k,102)*y(k,61)
         mat(k,1118) = -rxt(k,103)*y(k,61)
         mat(k,1159) = -rxt(k,104)*y(k,61)
         mat(k,1734) = -rxt(k,105)*y(k,61)
         mat(k,1454) = -rxt(k,142)*y(k,61)
         mat(k,1333) = -rxt(k,149)*y(k,61)
         mat(k,1540) = -rxt(k,172)*y(k,61)
         mat(k,816) = -rxt(k,190)*y(k,61)
         mat(k,1032) = -rxt(k,207)*y(k,61)
         mat(k,905) = -rxt(k,225)*y(k,61)
         mat(k,862) = -rxt(k,243)*y(k,61)
         mat(k,1497) = rxt(k,111)*y(k,65) + rxt(k,281)*y(k,81) + rxt(k,346)*y(k,106)
         mat(k,1118) = mat(k,1118) + rxt(k,125)*y(k,68) + rxt(k,284)*y(k,81) &
                      + rxt(k,293)*y(k,83) + rxt(k,313)*y(k,98) + rxt(k,340)*y(k,105) &
                      + rxt(k,351)*y(k,106)
         mat(k,991) = rxt(k,107)*y(k,64)
         mat(k,1214) = rxt(k,108)*y(k,64)
         mat(k,1376) = rxt(k,109)*y(k,65) + rxt(k,122)*y(k,68) + rxt(k,280)*y(k,81) &
                      + rxt(k,309)*y(k,98) + rxt(k,336)*y(k,105) + rxt(k,345)*y(k,106)
         mat(k,222) = rxt(k,107)*y(k,51) + rxt(k,108)*y(k,54)
         mat(k,775) = rxt(k,111)*y(k,24) + rxt(k,109)*y(k,60)
         mat(k,612) = rxt(k,125)*y(k,40) + rxt(k,122)*y(k,60)
         mat(k,704) = rxt(k,281)*y(k,24) + rxt(k,284)*y(k,40) + rxt(k,280)*y(k,60)
         mat(k,645) = rxt(k,293)*y(k,40)
         mat(k,1288) = rxt(k,313)*y(k,40) + rxt(k,309)*y(k,60)
         mat(k,1574) = rxt(k,340)*y(k,40) + rxt(k,336)*y(k,60)
         mat(k,1610) = rxt(k,346)*y(k,24) + rxt(k,351)*y(k,40) + rxt(k,345)*y(k,60)
         mat(k,350) = -(rxt(k,101)*y(k,40) + rxt(k,154)*y(k,99) + rxt(k,159)*y(k,75) &
                      + rxt(k,177)*y(k,77) + rxt(k,194)*y(k,73) + rxt(k,212)*y(k,92) &
                      + rxt(k,229)*y(k,89) + rxt(k,247)*y(k,88))
         mat(k,1112) = -rxt(k,101)*y(k,62)
         mat(k,1325) = -rxt(k,154)*y(k,62)
         mat(k,1447) = -rxt(k,159)*y(k,62)
         mat(k,1533) = -rxt(k,177)*y(k,62)
         mat(k,809) = -rxt(k,194)*y(k,62)
         mat(k,1025) = -rxt(k,212)*y(k,62)
         mat(k,898) = -rxt(k,229)*y(k,62)
         mat(k,854) = -rxt(k,247)*y(k,62)
         mat(k,590) = rxt(k,102)*y(k,110)
         mat(k,1779) = rxt(k,102)*y(k,61)
         mat(k,334) = -((rxt(k,164) + rxt(k,263)) * y(k,77) + (rxt(k,176) + rxt(k,262) &
                      ) * y(k,75) + rxt(k,182)*y(k,73) + rxt(k,200)*y(k,92) + rxt(k,217) &
                      *y(k,89) + rxt(k,235)*y(k,88) + rxt(k,252)*y(k,99))
         mat(k,1532) = -(rxt(k,164) + rxt(k,263)) * y(k,63)
         mat(k,1446) = -(rxt(k,176) + rxt(k,262)) * y(k,63)
         mat(k,808) = -rxt(k,182)*y(k,63)
         mat(k,1024) = -rxt(k,200)*y(k,63)
         mat(k,897) = -rxt(k,217)*y(k,63)
         mat(k,853) = -rxt(k,235)*y(k,63)
         mat(k,1324) = -rxt(k,252)*y(k,63)
         mat(k,1111) = rxt(k,103)*y(k,61) + rxt(k,101)*y(k,62)
         mat(k,589) = rxt(k,103)*y(k,40)
         mat(k,349) = rxt(k,101)*y(k,40)
         mat(k,221) = -((rxt(k,106) + rxt(k,107)) * y(k,51) + rxt(k,108)*y(k,54))
         mat(k,980) = -(rxt(k,106) + rxt(k,107)) * y(k,64)
         mat(k,1204) = -rxt(k,108)*y(k,64)
         mat(k,1491) = rxt(k,123)*y(k,68) + rxt(k,310)*y(k,98) + rxt(k,337)*y(k,105)
         mat(k,1369) = rxt(k,110)*y(k,65)
         mat(k,771) = rxt(k,110)*y(k,60)
         mat(k,610) = rxt(k,123)*y(k,24)
         mat(k,1284) = rxt(k,310)*y(k,24)
         mat(k,1571) = rxt(k,337)*y(k,24)
         mat(k,779) = -((rxt(k,109) + rxt(k,110)) * y(k,60) + rxt(k,111)*y(k,24) &
                      + rxt(k,112)*y(k,32) + rxt(k,114)*y(k,110) + rxt(k,119)*y(k,42) &
                      + rxt(k,120)*y(k,54) + rxt(k,121)*y(k,55) + (rxt(k,163) &
                      + rxt(k,260)) * y(k,77) + (rxt(k,165) + rxt(k,253)) * y(k,75) &
                      + rxt(k,181)*y(k,73) + rxt(k,199)*y(k,92) + rxt(k,216)*y(k,89) &
                      + rxt(k,234)*y(k,88) + rxt(k,251)*y(k,99) + rxt(k,275)*y(k,51) &
                      + rxt(k,276)*y(k,52))
         mat(k,1382) = -(rxt(k,109) + rxt(k,110)) * y(k,65)
         mat(k,1503) = -rxt(k,111)*y(k,65)
         mat(k,1694) = -rxt(k,112)*y(k,65)
         mat(k,1792) = -rxt(k,114)*y(k,65)
         mat(k,1164) = -rxt(k,119)*y(k,65)
         mat(k,1220) = -rxt(k,120)*y(k,65)
         mat(k,1658) = -rxt(k,121)*y(k,65)
         mat(k,1546) = -(rxt(k,163) + rxt(k,260)) * y(k,65)
         mat(k,1460) = -(rxt(k,165) + rxt(k,253)) * y(k,65)
         mat(k,822) = -rxt(k,181)*y(k,65)
         mat(k,1038) = -rxt(k,199)*y(k,65)
         mat(k,911) = -rxt(k,216)*y(k,65)
         mat(k,868) = -rxt(k,234)*y(k,65)
         mat(k,1339) = -rxt(k,251)*y(k,65)
         mat(k,996) = -rxt(k,275)*y(k,65)
         mat(k,1739) = -rxt(k,276)*y(k,65)
         mat(k,1079) = rxt(k,325)*y(k,101) + rxt(k,347)*y(k,106)
         mat(k,1694) = mat(k,1694) + rxt(k,124)*y(k,68)
         mat(k,1220) = mat(k,1220) + rxt(k,126)*y(k,68)
         mat(k,618) = rxt(k,124)*y(k,32) + rxt(k,126)*y(k,54)
         mat(k,398) = rxt(k,325)*y(k,28)
         mat(k,1615) = rxt(k,347)*y(k,28)
         mat(k,305) = -(rxt(k,148)*y(k,99) + (rxt(k,171) + rxt(k,254)) * y(k,77) &
                      + rxt(k,189)*y(k,73) + rxt(k,206)*y(k,92) + rxt(k,224)*y(k,89) &
                      + rxt(k,241)*y(k,88) + (rxt(k,242) + rxt(k,266)) * y(k,75))
         mat(k,1323) = -rxt(k,148)*y(k,66)
         mat(k,1531) = -(rxt(k,171) + rxt(k,254)) * y(k,66)
         mat(k,807) = -rxt(k,189)*y(k,66)
         mat(k,1023) = -rxt(k,206)*y(k,66)
         mat(k,896) = -rxt(k,224)*y(k,66)
         mat(k,852) = -rxt(k,241)*y(k,66)
         mat(k,1445) = -(rxt(k,242) + rxt(k,266)) * y(k,66)
         mat(k,489) = rxt(k,113)*y(k,110)
         mat(k,1778) = rxt(k,113)*y(k,67)
         mat(k,491) = -(rxt(k,113)*y(k,110) + (rxt(k,115) + rxt(k,116)) * y(k,52) &
                      + (rxt(k,117) + rxt(k,118)) * y(k,51) + rxt(k,150)*y(k,99) &
                      + (rxt(k,153) + rxt(k,255)) * y(k,75) + (rxt(k,173) + rxt(k,256) &
                      ) * y(k,77) + rxt(k,191)*y(k,73) + rxt(k,208)*y(k,92) + rxt(k,226) &
                      *y(k,89) + rxt(k,244)*y(k,88))
         mat(k,1783) = -rxt(k,113)*y(k,67)
         mat(k,1730) = -(rxt(k,115) + rxt(k,116)) * y(k,67)
         mat(k,987) = -(rxt(k,117) + rxt(k,118)) * y(k,67)
         mat(k,1330) = -rxt(k,150)*y(k,67)
         mat(k,1452) = -(rxt(k,153) + rxt(k,255)) * y(k,67)
         mat(k,1538) = -(rxt(k,173) + rxt(k,256)) * y(k,67)
         mat(k,814) = -rxt(k,191)*y(k,67)
         mat(k,1030) = -rxt(k,208)*y(k,67)
         mat(k,903) = -rxt(k,226)*y(k,67)
         mat(k,859) = -rxt(k,244)*y(k,67)
         mat(k,773) = rxt(k,114)*y(k,110)
         mat(k,1783) = mat(k,1783) + rxt(k,114)*y(k,65)
         mat(k,613) = -(rxt(k,122)*y(k,60) + rxt(k,123)*y(k,24) + rxt(k,124)*y(k,32) &
                      + rxt(k,125)*y(k,40) + rxt(k,126)*y(k,54) + rxt(k,127)*y(k,56) &
                      + rxt(k,146)*y(k,99) + rxt(k,169)*y(k,77) + rxt(k,186)*y(k,73) &
                      + rxt(k,204)*y(k,92) + rxt(k,220)*y(k,75) + rxt(k,222)*y(k,89) &
                      + rxt(k,239)*y(k,88))
         mat(k,1377) = -rxt(k,122)*y(k,68)
         mat(k,1498) = -rxt(k,123)*y(k,68)
         mat(k,1689) = -rxt(k,124)*y(k,68)
         mat(k,1119) = -rxt(k,125)*y(k,68)
         mat(k,1215) = -rxt(k,126)*y(k,68)
         mat(k,1253) = -rxt(k,127)*y(k,68)
         mat(k,1334) = -rxt(k,146)*y(k,68)
         mat(k,1541) = -rxt(k,169)*y(k,68)
         mat(k,817) = -rxt(k,186)*y(k,68)
         mat(k,1033) = -rxt(k,204)*y(k,68)
         mat(k,1455) = -rxt(k,220)*y(k,68)
         mat(k,906) = -rxt(k,222)*y(k,68)
         mat(k,863) = -rxt(k,239)*y(k,68)
         mat(k,1074) = rxt(k,311)*y(k,98) + rxt(k,330)*y(k,102)
         mat(k,1289) = rxt(k,311)*y(k,28)
         mat(k,100) = rxt(k,330)*y(k,28)
         mat(k,958) = -((rxt(k,128) + rxt(k,130)) * y(k,55) + rxt(k,129)*y(k,56) &
                      + rxt(k,133)*y(k,70) + rxt(k,136)*y(k,75) + rxt(k,139)*y(k,77) &
                      + rxt(k,298)*y(k,89) + rxt(k,299)*y(k,90) + rxt(k,301)*y(k,91) &
                      + rxt(k,303)*y(k,92) + rxt(k,321)*y(k,100) + rxt(k,375)*y(k,103) &
                      + rxt(k,376)*y(k,78) + rxt(k,377)*y(k,72) + rxt(k,378)*y(k,73) &
                      + rxt(k,379)*y(k,93) + rxt(k,561)*y(k,88) + rxt(k,562)*y(k,99) &
                      + rxt(k,563)*y(k,80))
         mat(k,1662) = -(rxt(k,128) + rxt(k,130)) * y(k,69)
         mat(k,1262) = -rxt(k,129)*y(k,69)
         mat(k,36) = -rxt(k,133)*y(k,69)
         mat(k,1464) = -rxt(k,136)*y(k,69)
         mat(k,1550) = -rxt(k,139)*y(k,69)
         mat(k,915) = -rxt(k,298)*y(k,69)
         mat(k,48) = -rxt(k,299)*y(k,69)
         mat(k,132) = -rxt(k,301)*y(k,69)
         mat(k,1042) = -rxt(k,303)*y(k,69)
         mat(k,159) = -rxt(k,321)*y(k,69)
         mat(k,282) = -rxt(k,375)*y(k,69)
         mat(k,150) = -rxt(k,376)*y(k,69)
         mat(k,108) = -rxt(k,377)*y(k,69)
         mat(k,826) = -rxt(k,378)*y(k,69)
         mat(k,125) = -rxt(k,379)*y(k,69)
         mat(k,872) = -rxt(k,561)*y(k,69)
         mat(k,1343) = -rxt(k,562)*y(k,69)
         mat(k,117) = -rxt(k,563)*y(k,69)
         mat(k,1698) = rxt(k,100)*y(k,61) + rxt(k,312)*y(k,98) + rxt(k,339)*y(k,105)
         mat(k,418) = rxt(k,348)*y(k,106)
         mat(k,1000) = rxt(k,131)*y(k,106)
         mat(k,1224) = rxt(k,319)*y(k,98) + rxt(k,328)*y(k,101) + rxt(k,342)*y(k,105) &
                      + rxt(k,355)*y(k,106)
         mat(k,1262) = mat(k,1262) + rxt(k,327)*y(k,101)
         mat(k,597) = rxt(k,100)*y(k,32)
         mat(k,231) = rxt(k,316)*y(k,98) + rxt(k,356)*y(k,106)
         mat(k,1298) = rxt(k,312)*y(k,32) + rxt(k,319)*y(k,54) + rxt(k,316)*y(k,96)
         mat(k,399) = rxt(k,328)*y(k,54) + rxt(k,327)*y(k,56)
         mat(k,1583) = rxt(k,339)*y(k,32) + rxt(k,342)*y(k,54)
         mat(k,1619) = rxt(k,348)*y(k,33) + rxt(k,131)*y(k,51) + rxt(k,355)*y(k,54) &
                      + rxt(k,356)*y(k,96)
         mat(k,33) = -(rxt(k,133)*y(k,69) + rxt(k,134)*y(k,110))
         mat(k,936) = -rxt(k,133)*y(k,70)
         mat(k,1767) = -rxt(k,134)*y(k,70)
         mat(k,153) = rxt(k,322)*y(k,110)
         mat(k,1767) = mat(k,1767) + rxt(k,322)*y(k,100)
         mat(k,291) = -(rxt(k,144)*y(k,99) + rxt(k,167)*y(k,77) + rxt(k,184)*y(k,73) &
                      + rxt(k,198)*y(k,75) + rxt(k,202)*y(k,92) + rxt(k,219)*y(k,89) &
                      + rxt(k,237)*y(k,88))
         mat(k,1322) = -rxt(k,144)*y(k,71)
         mat(k,1530) = -rxt(k,167)*y(k,71)
         mat(k,806) = -rxt(k,184)*y(k,71)
         mat(k,1444) = -rxt(k,198)*y(k,71)
         mat(k,1022) = -rxt(k,202)*y(k,71)
         mat(k,895) = -rxt(k,219)*y(k,71)
         mat(k,851) = -rxt(k,237)*y(k,71)
         mat(k,1072) = rxt(k,338)*y(k,105)
         mat(k,1572) = rxt(k,338)*y(k,28)
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
         mat(k,105) = -(rxt(k,369)*y(k,110) + rxt(k,377)*y(k,69))
         mat(k,1769) = -rxt(k,369)*y(k,72)
         mat(k,942) = -rxt(k,377)*y(k,72)
         mat(k,34) = rxt(k,134)*y(k,110)
         mat(k,147) = rxt(k,367)*y(k,110)
         mat(k,1769) = mat(k,1769) + rxt(k,134)*y(k,70) + rxt(k,367)*y(k,78)
         mat(k,823) = -(rxt(k,180)*y(k,87) + rxt(k,181)*y(k,65) + rxt(k,182)*y(k,63) &
                      + rxt(k,183)*y(k,83) + rxt(k,184)*y(k,71) + rxt(k,185)*y(k,98) &
                      + rxt(k,186)*y(k,68) + rxt(k,188)*y(k,85) + rxt(k,189)*y(k,66) &
                      + rxt(k,190)*y(k,61) + rxt(k,191)*y(k,67) + rxt(k,192)*y(k,82) &
                      + rxt(k,193)*y(k,86) + rxt(k,194)*y(k,62) + rxt(k,195)*y(k,84) &
                      + rxt(k,196)*y(k,81) + rxt(k,371)*y(k,110) + rxt(k,378)*y(k,69))
         mat(k,434) = -rxt(k,180)*y(k,73)
         mat(k,780) = -rxt(k,181)*y(k,73)
         mat(k,336) = -rxt(k,182)*y(k,73)
         mat(k,648) = -rxt(k,183)*y(k,73)
         mat(k,293) = -rxt(k,184)*y(k,73)
         mat(k,1295) = -rxt(k,185)*y(k,73)
         mat(k,619) = -rxt(k,186)*y(k,73)
         mat(k,515) = -rxt(k,188)*y(k,73)
         mat(k,307) = -rxt(k,189)*y(k,73)
         mat(k,594) = -rxt(k,190)*y(k,73)
         mat(k,496) = -rxt(k,191)*y(k,73)
         mat(k,367) = -rxt(k,192)*y(k,73)
         mat(k,381) = -rxt(k,193)*y(k,73)
         mat(k,352) = -rxt(k,194)*y(k,73)
         mat(k,474) = -rxt(k,195)*y(k,73)
         mat(k,708) = -rxt(k,196)*y(k,73)
         mat(k,1793) = -rxt(k,371)*y(k,73)
         mat(k,955) = -rxt(k,378)*y(k,73)
         mat(k,107) = rxt(k,369)*y(k,110)
         mat(k,47) = rxt(k,300)*y(k,110)
         mat(k,1793) = mat(k,1793) + rxt(k,369)*y(k,72) + rxt(k,300)*y(k,90)
         mat(k,22) = -(rxt(k,135)*y(k,110))
         mat(k,1765) = -rxt(k,135)*y(k,74)
         mat(k,237) = rxt(k,137)*y(k,75)
         mat(k,1442) = rxt(k,137)*y(k,50)
         mat(k,1476) = -(rxt(k,136)*y(k,69) + rxt(k,137)*y(k,50) + (rxt(k,141) &
                      + rxt(k,264)) * y(k,87) + rxt(k,142)*y(k,61) + (rxt(k,153) &
                      + rxt(k,255)) * y(k,67) + rxt(k,157)*y(k,82) + rxt(k,158) &
                      *y(k,86) + rxt(k,159)*y(k,62) + rxt(k,160)*y(k,84) + rxt(k,161) &
                      *y(k,81) + (rxt(k,165) + rxt(k,253)) * y(k,65) + (rxt(k,176) &
                      + rxt(k,262)) * y(k,63) + (rxt(k,187) + rxt(k,259)) * y(k,83) &
                      + rxt(k,198)*y(k,71) + rxt(k,209)*y(k,98) + rxt(k,220)*y(k,68) &
                      + (rxt(k,231) + rxt(k,257)) * y(k,85) + (rxt(k,242) + rxt(k,266) &
                      ) * y(k,66) + rxt(k,373)*y(k,110))
         mat(k,970) = -rxt(k,136)*y(k,75)
         mat(k,247) = -rxt(k,137)*y(k,75)
         mat(k,442) = -(rxt(k,141) + rxt(k,264)) * y(k,75)
         mat(k,604) = -rxt(k,142)*y(k,75)
         mat(k,504) = -(rxt(k,153) + rxt(k,255)) * y(k,75)
         mat(k,373) = -rxt(k,157)*y(k,75)
         mat(k,390) = -rxt(k,158)*y(k,75)
         mat(k,360) = -rxt(k,159)*y(k,75)
         mat(k,483) = -rxt(k,160)*y(k,75)
         mat(k,721) = -rxt(k,161)*y(k,75)
         mat(k,795) = -(rxt(k,165) + rxt(k,253)) * y(k,75)
         mat(k,344) = -(rxt(k,176) + rxt(k,262)) * y(k,75)
         mat(k,661) = -(rxt(k,187) + rxt(k,259)) * y(k,75)
         mat(k,300) = -rxt(k,198)*y(k,75)
         mat(k,1310) = -rxt(k,209)*y(k,75)
         mat(k,633) = -rxt(k,220)*y(k,75)
         mat(k,524) = -(rxt(k,231) + rxt(k,257)) * y(k,75)
         mat(k,315) = -(rxt(k,242) + rxt(k,266)) * y(k,75)
         mat(k,1808) = -rxt(k,373)*y(k,75)
         mat(k,838) = rxt(k,371)*y(k,110)
         mat(k,24) = rxt(k,135)*y(k,110)
         mat(k,1808) = mat(k,1808) + rxt(k,371)*y(k,73) + rxt(k,135)*y(k,74)
         mat(k,26) = -(rxt(k,138)*y(k,110))
         mat(k,1766) = -rxt(k,138)*y(k,76)
         mat(k,238) = rxt(k,140)*y(k,77)
         mat(k,1528) = rxt(k,140)*y(k,50)
         mat(k,1564) = -(rxt(k,139)*y(k,69) + rxt(k,140)*y(k,50) + (rxt(k,162) &
                      + rxt(k,265)) * y(k,87) + (rxt(k,163) + rxt(k,260)) * y(k,65) &
                      + (rxt(k,164) + rxt(k,263)) * y(k,63) + (rxt(k,166) + rxt(k,261) &
                      ) * y(k,83) + rxt(k,167)*y(k,71) + rxt(k,168)*y(k,98) + rxt(k,169) &
                      *y(k,68) + (rxt(k,170) + rxt(k,258)) * y(k,85) + (rxt(k,171) &
                      + rxt(k,254)) * y(k,66) + rxt(k,172)*y(k,61) + (rxt(k,173) &
                      + rxt(k,256)) * y(k,67) + rxt(k,174)*y(k,82) + rxt(k,175) &
                      *y(k,86) + rxt(k,177)*y(k,62) + rxt(k,178)*y(k,84) + rxt(k,179) &
                      *y(k,81))
         mat(k,972) = -rxt(k,139)*y(k,77)
         mat(k,248) = -rxt(k,140)*y(k,77)
         mat(k,443) = -(rxt(k,162) + rxt(k,265)) * y(k,77)
         mat(k,797) = -(rxt(k,163) + rxt(k,260)) * y(k,77)
         mat(k,345) = -(rxt(k,164) + rxt(k,263)) * y(k,77)
         mat(k,662) = -(rxt(k,166) + rxt(k,261)) * y(k,77)
         mat(k,301) = -rxt(k,167)*y(k,77)
         mat(k,1312) = -rxt(k,168)*y(k,77)
         mat(k,635) = -rxt(k,169)*y(k,77)
         mat(k,525) = -(rxt(k,170) + rxt(k,258)) * y(k,77)
         mat(k,316) = -(rxt(k,171) + rxt(k,254)) * y(k,77)
         mat(k,605) = -rxt(k,172)*y(k,77)
         mat(k,505) = -(rxt(k,173) + rxt(k,256)) * y(k,77)
         mat(k,374) = -rxt(k,174)*y(k,77)
         mat(k,391) = -rxt(k,175)*y(k,77)
         mat(k,361) = -rxt(k,177)*y(k,77)
         mat(k,484) = -rxt(k,178)*y(k,77)
         mat(k,723) = -rxt(k,179)*y(k,77)
         mat(k,1478) = rxt(k,373)*y(k,110)
         mat(k,28) = rxt(k,138)*y(k,110)
         mat(k,1810) = rxt(k,373)*y(k,75) + rxt(k,138)*y(k,76)
         mat(k,148) = -(rxt(k,367)*y(k,110) + rxt(k,376)*y(k,69))
         mat(k,1773) = -rxt(k,367)*y(k,78)
         mat(k,947) = -rxt(k,376)*y(k,78)
         mat(k,1683) = rxt(k,304)*y(k,92)
         mat(k,670) = rxt(k,305)*y(k,92)
         mat(k,1021) = rxt(k,304)*y(k,32) + rxt(k,305)*y(k,43) + rxt(k,306)*y(k,104)
         mat(k,155) = rxt(k,323)*y(k,110)
         mat(k,736) = rxt(k,306)*y(k,92)
         mat(k,1773) = mat(k,1773) + rxt(k,323)*y(k,100)
         mat(k,71) = -(rxt(k,419)*y(k,54) + rxt(k,420)*y(k,55))
         mat(k,1192) = -rxt(k,419)*y(k,79)
         mat(k,1643) = -rxt(k,420)*y(k,79)
         mat(k,1192) = mat(k,1192) + rxt(k,565)*y(k,80)
         mat(k,940) = .900_r8*rxt(k,563)*y(k,80) + .800_r8*rxt(k,561)*y(k,88)
         mat(k,111) = rxt(k,565)*y(k,54) + .900_r8*rxt(k,563)*y(k,69)
         mat(k,847) = .800_r8*rxt(k,561)*y(k,69)
         mat(k,112) = -(rxt(k,563)*y(k,69) + rxt(k,564)*y(k,55) + (rxt(k,565) &
                      + rxt(k,566)) * y(k,54))
         mat(k,943) = -rxt(k,563)*y(k,80)
         mat(k,1646) = -rxt(k,564)*y(k,80)
         mat(k,1197) = -(rxt(k,565) + rxt(k,566)) * y(k,80)
         mat(k,706) = -(rxt(k,156)*y(k,99) + rxt(k,161)*y(k,75) + rxt(k,179)*y(k,77) &
                      + rxt(k,196)*y(k,73) + rxt(k,214)*y(k,92) + rxt(k,232)*y(k,89) &
                      + rxt(k,249)*y(k,88) + rxt(k,280)*y(k,60) + rxt(k,281)*y(k,24) &
                      + rxt(k,282)*y(k,32) + rxt(k,283)*y(k,110) + rxt(k,284)*y(k,40) &
                      + rxt(k,285)*y(k,42) + rxt(k,286)*y(k,52) + rxt(k,287)*y(k,56))
         mat(k,1337) = -rxt(k,156)*y(k,81)
         mat(k,1458) = -rxt(k,161)*y(k,81)
         mat(k,1544) = -rxt(k,179)*y(k,81)
         mat(k,820) = -rxt(k,196)*y(k,81)
         mat(k,1036) = -rxt(k,214)*y(k,81)
         mat(k,909) = -rxt(k,232)*y(k,81)
         mat(k,866) = -rxt(k,249)*y(k,81)
         mat(k,1380) = -rxt(k,280)*y(k,81)
         mat(k,1501) = -rxt(k,281)*y(k,81)
         mat(k,1692) = -rxt(k,282)*y(k,81)
         mat(k,1790) = -rxt(k,283)*y(k,81)
         mat(k,1122) = -rxt(k,284)*y(k,81)
         mat(k,1162) = -rxt(k,285)*y(k,81)
         mat(k,1737) = -rxt(k,286)*y(k,81)
         mat(k,1256) = -rxt(k,287)*y(k,81)
         mat(k,994) = rxt(k,106)*y(k,64) + rxt(k,275)*y(k,65) + rxt(k,118)*y(k,67) &
                      + rxt(k,274)*y(k,101)
         mat(k,1737) = mat(k,1737) + rxt(k,105)*y(k,61) + rxt(k,315)*y(k,98) &
                      + rxt(k,273)*y(k,101) + rxt(k,341)*y(k,105) + rxt(k,354) &
                      *y(k,106)
         mat(k,1218) = rxt(k,296)*y(k,83)
         mat(k,1256) = mat(k,1256) + rxt(k,297)*y(k,83)
         mat(k,593) = rxt(k,105)*y(k,52)
         mat(k,223) = rxt(k,106)*y(k,51)
         mat(k,777) = rxt(k,275)*y(k,51)
         mat(k,494) = rxt(k,118)*y(k,51)
         mat(k,647) = rxt(k,296)*y(k,54) + rxt(k,297)*y(k,56)
         mat(k,1292) = rxt(k,315)*y(k,52)
         mat(k,397) = rxt(k,274)*y(k,51) + rxt(k,273)*y(k,52)
         mat(k,1577) = rxt(k,341)*y(k,52)
         mat(k,1613) = rxt(k,354)*y(k,52)
         mat(k,365) = -(rxt(k,151)*y(k,99) + rxt(k,157)*y(k,75) + rxt(k,174)*y(k,77) &
                      + rxt(k,192)*y(k,73) + rxt(k,210)*y(k,92) + rxt(k,227)*y(k,89) &
                      + rxt(k,245)*y(k,88))
         mat(k,1326) = -rxt(k,151)*y(k,82)
         mat(k,1448) = -rxt(k,157)*y(k,82)
         mat(k,1534) = -rxt(k,174)*y(k,82)
         mat(k,810) = -rxt(k,192)*y(k,82)
         mat(k,1026) = -rxt(k,210)*y(k,82)
         mat(k,899) = -rxt(k,227)*y(k,82)
         mat(k,855) = -rxt(k,245)*y(k,82)
         mat(k,984) = rxt(k,117)*y(k,67)
         mat(k,490) = rxt(k,117)*y(k,51)
         mat(k,703) = rxt(k,283)*y(k,110)
         mat(k,1780) = rxt(k,283)*y(k,81)
         mat(k,646) = -(rxt(k,143)*y(k,99) + (rxt(k,166) + rxt(k,261)) * y(k,77) &
                      + rxt(k,183)*y(k,73) + (rxt(k,187) + rxt(k,259)) * y(k,75) &
                      + rxt(k,201)*y(k,92) + rxt(k,218)*y(k,89) + rxt(k,236)*y(k,88) &
                      + (rxt(k,271) + rxt(k,293)) * y(k,40) + rxt(k,291)*y(k,110) &
                      + rxt(k,295)*y(k,42) + rxt(k,296)*y(k,54) + rxt(k,297)*y(k,56))
         mat(k,1335) = -rxt(k,143)*y(k,83)
         mat(k,1542) = -(rxt(k,166) + rxt(k,261)) * y(k,83)
         mat(k,818) = -rxt(k,183)*y(k,83)
         mat(k,1456) = -(rxt(k,187) + rxt(k,259)) * y(k,83)
         mat(k,1034) = -rxt(k,201)*y(k,83)
         mat(k,907) = -rxt(k,218)*y(k,83)
         mat(k,864) = -rxt(k,236)*y(k,83)
         mat(k,1120) = -(rxt(k,271) + rxt(k,293)) * y(k,83)
         mat(k,1788) = -rxt(k,291)*y(k,83)
         mat(k,1160) = -rxt(k,295)*y(k,83)
         mat(k,1216) = -rxt(k,296)*y(k,83)
         mat(k,1254) = -rxt(k,297)*y(k,83)
         mat(k,1160) = mat(k,1160) + rxt(k,104)*y(k,61) + rxt(k,119)*y(k,65) &
                      + rxt(k,285)*y(k,81) + rxt(k,314)*y(k,98) + rxt(k,352)*y(k,106)
         mat(k,992) = rxt(k,267)*y(k,101)
         mat(k,1735) = rxt(k,276)*y(k,65) + rxt(k,115)*y(k,67) + rxt(k,286)*y(k,81) &
                      + rxt(k,272)*y(k,101)
         mat(k,1254) = mat(k,1254) + rxt(k,287)*y(k,81)
         mat(k,592) = rxt(k,104)*y(k,42)
         mat(k,776) = rxt(k,119)*y(k,42) + rxt(k,276)*y(k,52)
         mat(k,493) = rxt(k,115)*y(k,52)
         mat(k,705) = rxt(k,285)*y(k,42) + rxt(k,286)*y(k,52) + rxt(k,287)*y(k,56)
         mat(k,1290) = rxt(k,314)*y(k,42)
         mat(k,396) = rxt(k,267)*y(k,51) + rxt(k,272)*y(k,52)
         mat(k,1611) = rxt(k,352)*y(k,42)
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
         mat(k,471) = -(rxt(k,155)*y(k,99) + rxt(k,160)*y(k,75) + rxt(k,178)*y(k,77) &
                      + rxt(k,195)*y(k,73) + rxt(k,213)*y(k,92) + rxt(k,230)*y(k,89) &
                      + rxt(k,248)*y(k,88) + rxt(k,288)*y(k,50))
         mat(k,1329) = -rxt(k,155)*y(k,84)
         mat(k,1451) = -rxt(k,160)*y(k,84)
         mat(k,1537) = -rxt(k,178)*y(k,84)
         mat(k,813) = -rxt(k,195)*y(k,84)
         mat(k,1029) = -rxt(k,213)*y(k,84)
         mat(k,902) = -rxt(k,230)*y(k,84)
         mat(k,858) = -rxt(k,248)*y(k,84)
         mat(k,241) = -rxt(k,288)*y(k,84)
         mat(k,512) = rxt(k,289)*y(k,110)
         mat(k,1782) = rxt(k,289)*y(k,85)
         mat(k,513) = -(rxt(k,147)*y(k,99) + (rxt(k,170) + rxt(k,258)) * y(k,77) &
                      + rxt(k,188)*y(k,73) + rxt(k,205)*y(k,92) + rxt(k,223)*y(k,89) &
                      + (rxt(k,231) + rxt(k,257)) * y(k,75) + rxt(k,240)*y(k,88) &
                      + rxt(k,289)*y(k,110) + rxt(k,290)*y(k,42) + rxt(k,292)*y(k,50))
         mat(k,1331) = -rxt(k,147)*y(k,85)
         mat(k,1539) = -(rxt(k,170) + rxt(k,258)) * y(k,85)
         mat(k,815) = -rxt(k,188)*y(k,85)
         mat(k,1031) = -rxt(k,205)*y(k,85)
         mat(k,904) = -rxt(k,223)*y(k,85)
         mat(k,1453) = -(rxt(k,231) + rxt(k,257)) * y(k,85)
         mat(k,860) = -rxt(k,240)*y(k,85)
         mat(k,1784) = -rxt(k,289)*y(k,85)
         mat(k,1156) = -rxt(k,290)*y(k,85)
         mat(k,242) = -rxt(k,292)*y(k,85)
         mat(k,1731) = rxt(k,116)*y(k,67)
         mat(k,492) = rxt(k,116)*y(k,52)
         mat(k,644) = rxt(k,291)*y(k,110)
         mat(k,1784) = mat(k,1784) + rxt(k,291)*y(k,83)
         mat(k,379) = -(rxt(k,152)*y(k,99) + rxt(k,158)*y(k,75) + rxt(k,175)*y(k,77) &
                      + rxt(k,193)*y(k,73) + rxt(k,211)*y(k,92) + rxt(k,228)*y(k,89) &
                      + rxt(k,246)*y(k,88) + rxt(k,294)*y(k,42))
         mat(k,1327) = -rxt(k,152)*y(k,86)
         mat(k,1449) = -rxt(k,158)*y(k,86)
         mat(k,1535) = -rxt(k,175)*y(k,86)
         mat(k,811) = -rxt(k,193)*y(k,86)
         mat(k,1027) = -rxt(k,211)*y(k,86)
         mat(k,900) = -rxt(k,228)*y(k,86)
         mat(k,856) = -rxt(k,246)*y(k,86)
         mat(k,1153) = -rxt(k,294)*y(k,86)
         mat(k,1113) = rxt(k,271)*y(k,83)
         mat(k,642) = rxt(k,271)*y(k,40)
         mat(k,432) = -((rxt(k,141) + rxt(k,264)) * y(k,75) + (rxt(k,162) + rxt(k,265) &
                      ) * y(k,77) + rxt(k,180)*y(k,73) + rxt(k,197)*y(k,92) + rxt(k,215) &
                      *y(k,89) + rxt(k,233)*y(k,88) + rxt(k,250)*y(k,99))
         mat(k,1450) = -(rxt(k,141) + rxt(k,264)) * y(k,87)
         mat(k,1536) = -(rxt(k,162) + rxt(k,265)) * y(k,87)
         mat(k,812) = -rxt(k,180)*y(k,87)
         mat(k,1028) = -rxt(k,197)*y(k,87)
         mat(k,901) = -rxt(k,215)*y(k,87)
         mat(k,857) = -rxt(k,233)*y(k,87)
         mat(k,1328) = -rxt(k,250)*y(k,87)
         mat(k,1155) = rxt(k,295)*y(k,83) + rxt(k,290)*y(k,85) + rxt(k,294)*y(k,86)
         mat(k,240) = rxt(k,288)*y(k,84) + rxt(k,292)*y(k,85)
         mat(k,643) = rxt(k,295)*y(k,42)
         mat(k,470) = rxt(k,288)*y(k,50)
         mat(k,511) = rxt(k,290)*y(k,42) + rxt(k,292)*y(k,50)
         mat(k,380) = rxt(k,294)*y(k,42)
         mat(k,870) = -(rxt(k,233)*y(k,87) + rxt(k,234)*y(k,65) + rxt(k,235)*y(k,63) &
                      + rxt(k,236)*y(k,83) + rxt(k,237)*y(k,71) + rxt(k,238)*y(k,98) &
                      + rxt(k,239)*y(k,68) + rxt(k,240)*y(k,85) + rxt(k,241)*y(k,66) &
                      + rxt(k,243)*y(k,61) + rxt(k,244)*y(k,67) + rxt(k,245)*y(k,82) &
                      + rxt(k,246)*y(k,86) + rxt(k,247)*y(k,62) + rxt(k,248)*y(k,84) &
                      + rxt(k,249)*y(k,81) + rxt(k,360)*y(k,110) + rxt(k,363)*y(k,28) &
                      + rxt(k,561)*y(k,69))
         mat(k,435) = -rxt(k,233)*y(k,88)
         mat(k,781) = -rxt(k,234)*y(k,88)
         mat(k,337) = -rxt(k,235)*y(k,88)
         mat(k,649) = -rxt(k,236)*y(k,88)
         mat(k,294) = -rxt(k,237)*y(k,88)
         mat(k,1296) = -rxt(k,238)*y(k,88)
         mat(k,620) = -rxt(k,239)*y(k,88)
         mat(k,516) = -rxt(k,240)*y(k,88)
         mat(k,308) = -rxt(k,241)*y(k,88)
         mat(k,595) = -rxt(k,243)*y(k,88)
         mat(k,497) = -rxt(k,244)*y(k,88)
         mat(k,368) = -rxt(k,245)*y(k,88)
         mat(k,382) = -rxt(k,246)*y(k,88)
         mat(k,353) = -rxt(k,247)*y(k,88)
         mat(k,475) = -rxt(k,248)*y(k,88)
         mat(k,709) = -rxt(k,249)*y(k,88)
         mat(k,1794) = -rxt(k,360)*y(k,88)
         mat(k,1081) = -rxt(k,363)*y(k,88)
         mat(k,956) = -rxt(k,561)*y(k,88)
         mat(k,268) = rxt(k,570)*y(k,99)
         mat(k,998) = rxt(k,572)*y(k,99)
         mat(k,1222) = rxt(k,565)*y(k,80)
         mat(k,1660) = rxt(k,569)*y(k,94)
         mat(k,116) = rxt(k,565)*y(k,54)
         mat(k,94) = rxt(k,569)*y(k,55)
         mat(k,1341) = rxt(k,570)*y(k,48) + rxt(k,572)*y(k,51)
         mat(k,914) = -(rxt(k,215)*y(k,87) + rxt(k,216)*y(k,65) + rxt(k,217)*y(k,63) &
                      + rxt(k,218)*y(k,83) + rxt(k,219)*y(k,71) + rxt(k,221)*y(k,98) &
                      + rxt(k,222)*y(k,68) + rxt(k,223)*y(k,85) + rxt(k,224)*y(k,66) &
                      + rxt(k,225)*y(k,61) + rxt(k,226)*y(k,67) + rxt(k,227)*y(k,82) &
                      + rxt(k,228)*y(k,86) + rxt(k,229)*y(k,62) + rxt(k,230)*y(k,84) &
                      + rxt(k,232)*y(k,81) + rxt(k,298)*y(k,69) + rxt(k,362)*y(k,110))
         mat(k,436) = -rxt(k,215)*y(k,89)
         mat(k,782) = -rxt(k,216)*y(k,89)
         mat(k,338) = -rxt(k,217)*y(k,89)
         mat(k,650) = -rxt(k,218)*y(k,89)
         mat(k,295) = -rxt(k,219)*y(k,89)
         mat(k,1297) = -rxt(k,221)*y(k,89)
         mat(k,621) = -rxt(k,222)*y(k,89)
         mat(k,517) = -rxt(k,223)*y(k,89)
         mat(k,309) = -rxt(k,224)*y(k,89)
         mat(k,596) = -rxt(k,225)*y(k,89)
         mat(k,498) = -rxt(k,226)*y(k,89)
         mat(k,369) = -rxt(k,227)*y(k,89)
         mat(k,383) = -rxt(k,228)*y(k,89)
         mat(k,354) = -rxt(k,229)*y(k,89)
         mat(k,476) = -rxt(k,230)*y(k,89)
         mat(k,710) = -rxt(k,232)*y(k,89)
         mat(k,957) = -rxt(k,298)*y(k,89)
         mat(k,1795) = -rxt(k,362)*y(k,89)
         mat(k,1041) = rxt(k,361)*y(k,110)
         mat(k,1795) = mat(k,1795) + rxt(k,361)*y(k,92)
         mat(k,45) = -(rxt(k,299)*y(k,69) + rxt(k,300)*y(k,110))
         mat(k,937) = -rxt(k,299)*y(k,90)
         mat(k,1768) = -rxt(k,300)*y(k,90)
         mat(k,893) = rxt(k,362)*y(k,110)
         mat(k,1768) = mat(k,1768) + rxt(k,362)*y(k,89)
         mat(k,130) = -(rxt(k,301)*y(k,69) + rxt(k,302)*y(k,110))
         mat(k,945) = -rxt(k,301)*y(k,91)
         mat(k,1771) = -rxt(k,302)*y(k,91)
         mat(k,1066) = rxt(k,363)*y(k,88) + rxt(k,307)*y(k,93)
         mat(k,849) = rxt(k,363)*y(k,28)
         mat(k,123) = rxt(k,307)*y(k,28)
         mat(k,1044) = -(rxt(k,197)*y(k,87) + rxt(k,199)*y(k,65) + rxt(k,200)*y(k,63) &
                      + rxt(k,201)*y(k,83) + rxt(k,202)*y(k,71) + rxt(k,203)*y(k,98) &
                      + rxt(k,204)*y(k,68) + rxt(k,205)*y(k,85) + rxt(k,206)*y(k,66) &
                      + rxt(k,207)*y(k,61) + rxt(k,208)*y(k,67) + rxt(k,210)*y(k,82) &
                      + rxt(k,211)*y(k,86) + rxt(k,212)*y(k,62) + rxt(k,213)*y(k,84) &
                      + rxt(k,214)*y(k,81) + rxt(k,303)*y(k,69) + rxt(k,304)*y(k,32) &
                      + rxt(k,305)*y(k,43) + rxt(k,306)*y(k,104) + rxt(k,361)*y(k,110))
         mat(k,438) = -rxt(k,197)*y(k,92)
         mat(k,785) = -rxt(k,199)*y(k,92)
         mat(k,340) = -rxt(k,200)*y(k,92)
         mat(k,653) = -rxt(k,201)*y(k,92)
         mat(k,297) = -rxt(k,202)*y(k,92)
         mat(k,1300) = -rxt(k,203)*y(k,92)
         mat(k,624) = -rxt(k,204)*y(k,92)
         mat(k,519) = -rxt(k,205)*y(k,92)
         mat(k,311) = -rxt(k,206)*y(k,92)
         mat(k,599) = -rxt(k,207)*y(k,92)
         mat(k,500) = -rxt(k,208)*y(k,92)
         mat(k,371) = -rxt(k,210)*y(k,92)
         mat(k,385) = -rxt(k,211)*y(k,92)
         mat(k,356) = -rxt(k,212)*y(k,92)
         mat(k,478) = -rxt(k,213)*y(k,92)
         mat(k,713) = -rxt(k,214)*y(k,92)
         mat(k,960) = -rxt(k,303)*y(k,92)
         mat(k,1700) = -rxt(k,304)*y(k,92)
         mat(k,685) = -rxt(k,305)*y(k,92)
         mat(k,754) = -rxt(k,306)*y(k,92)
         mat(k,1798) = -rxt(k,361)*y(k,92)
         mat(k,874) = rxt(k,360)*y(k,110)
         mat(k,134) = rxt(k,302)*y(k,110)
         mat(k,127) = rxt(k,308)*y(k,110)
         mat(k,1798) = mat(k,1798) + rxt(k,360)*y(k,88) + rxt(k,302)*y(k,91) &
                      + rxt(k,308)*y(k,93)
         mat(k,122) = -(rxt(k,307)*y(k,28) + rxt(k,308)*y(k,110) + rxt(k,379)*y(k,69))
         mat(k,1065) = -rxt(k,307)*y(k,93)
         mat(k,1770) = -rxt(k,308)*y(k,93)
         mat(k,944) = -rxt(k,379)*y(k,93)
         mat(k,91) = -(rxt(k,567)*y(k,54) + (rxt(k,568) + rxt(k,569)) * y(k,55))
         mat(k,1195) = -rxt(k,567)*y(k,94)
         mat(k,1644) = -(rxt(k,568) + rxt(k,569)) * y(k,94)
         mat(k,567) = -(rxt(k,384)*y(k,33) + rxt(k,385)*y(k,110) + (rxt(k,387) &
                      + rxt(k,388)) * y(k,55) + rxt(k,389)*y(k,56) + (rxt(k,477) &
                      + rxt(k,478)) * y(k,40) + (rxt(k,500) + rxt(k,501)) * y(k,36) &
                      + rxt(k,506)*y(k,29) + rxt(k,507)*y(k,30))
         mat(k,415) = -rxt(k,384)*y(k,95)
         mat(k,1786) = -rxt(k,385)*y(k,95)
         mat(k,1652) = -(rxt(k,387) + rxt(k,388)) * y(k,95)
         mat(k,1252) = -rxt(k,389)*y(k,95)
         mat(k,1117) = -(rxt(k,477) + rxt(k,478)) * y(k,95)
         mat(k,185) = -(rxt(k,500) + rxt(k,501)) * y(k,95)
         mat(k,6) = -rxt(k,506)*y(k,95)
         mat(k,14) = -rxt(k,507)*y(k,95)
         mat(k,1652) = mat(k,1652) + rxt(k,420)*y(k,79)
         mat(k,952) = .850_r8*rxt(k,562)*y(k,99)
         mat(k,73) = rxt(k,420)*y(k,55)
         mat(k,1332) = .850_r8*rxt(k,562)*y(k,69)
         mat(k,229) = -(rxt(k,316)*y(k,98) + rxt(k,334)*y(k,103) + rxt(k,356)*y(k,106) &
                      + rxt(k,391)*y(k,54) + rxt(k,392)*y(k,55))
         mat(k,1285) = -rxt(k,316)*y(k,96)
         mat(k,278) = -rxt(k,334)*y(k,96)
         mat(k,1604) = -rxt(k,356)*y(k,96)
         mat(k,1205) = -rxt(k,391)*y(k,96)
         mat(k,1648) = -rxt(k,392)*y(k,96)
         mat(k,1069) = rxt(k,393)*y(k,97)
         mat(k,1205) = mat(k,1205) + rxt(k,395)*y(k,97)
         mat(k,1648) = mat(k,1648) + rxt(k,396)*y(k,97)
         mat(k,1246) = rxt(k,397)*y(k,97)
         mat(k,31) = rxt(k,393)*y(k,28) + rxt(k,395)*y(k,54) + rxt(k,396)*y(k,55) &
                      + rxt(k,397)*y(k,56)
         mat(k,30) = -(rxt(k,393)*y(k,28) + rxt(k,395)*y(k,54) + rxt(k,396)*y(k,55) &
                      + rxt(k,397)*y(k,56))
         mat(k,1063) = -rxt(k,393)*y(k,97)
         mat(k,1189) = -rxt(k,395)*y(k,97)
         mat(k,1640) = -rxt(k,396)*y(k,97)
         mat(k,1245) = -rxt(k,397)*y(k,97)
         mat(k,1640) = mat(k,1640) + rxt(k,387)*y(k,95)
         mat(k,557) = rxt(k,387)*y(k,55)
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
         mat(k,1306) = -(rxt(k,145)*y(k,99) + rxt(k,168)*y(k,77) + rxt(k,185)*y(k,73) &
                      + rxt(k,203)*y(k,92) + rxt(k,209)*y(k,75) + rxt(k,221)*y(k,89) &
                      + rxt(k,238)*y(k,88) + rxt(k,309)*y(k,60) + rxt(k,310)*y(k,24) &
                      + rxt(k,311)*y(k,28) + rxt(k,312)*y(k,32) + rxt(k,313)*y(k,40) &
                      + rxt(k,314)*y(k,42) + rxt(k,315)*y(k,52) + rxt(k,316)*y(k,96) &
                      + rxt(k,317)*y(k,55) + rxt(k,318)*y(k,56) + (rxt(k,319) &
                      + rxt(k,320)) * y(k,54))
         mat(k,1351) = -rxt(k,145)*y(k,98)
         mat(k,1558) = -rxt(k,168)*y(k,98)
         mat(k,834) = -rxt(k,185)*y(k,98)
         mat(k,1050) = -rxt(k,203)*y(k,98)
         mat(k,1472) = -rxt(k,209)*y(k,98)
         mat(k,923) = -rxt(k,221)*y(k,98)
         mat(k,880) = -rxt(k,238)*y(k,98)
         mat(k,1394) = -rxt(k,309)*y(k,98)
         mat(k,1515) = -rxt(k,310)*y(k,98)
         mat(k,1091) = -rxt(k,311)*y(k,98)
         mat(k,1706) = -rxt(k,312)*y(k,98)
         mat(k,1136) = -rxt(k,313)*y(k,98)
         mat(k,1176) = -rxt(k,314)*y(k,98)
         mat(k,1751) = -rxt(k,315)*y(k,98)
         mat(k,233) = -rxt(k,316)*y(k,98)
         mat(k,1670) = -rxt(k,317)*y(k,98)
         mat(k,1270) = -rxt(k,318)*y(k,98)
         mat(k,1232) = -(rxt(k,319) + rxt(k,320)) * y(k,98)
         mat(k,1232) = mat(k,1232) + rxt(k,120)*y(k,65) + rxt(k,329)*y(k,101)
         mat(k,1670) = mat(k,1670) + (rxt(k,128)+rxt(k,130))*y(k,69)
         mat(k,791) = rxt(k,120)*y(k,54)
         mat(k,966) = (rxt(k,128)+rxt(k,130))*y(k,55)
         mat(k,404) = rxt(k,329)*y(k,54)
         mat(k,1352) = -(rxt(k,143)*y(k,83) + rxt(k,144)*y(k,71) + rxt(k,145)*y(k,98) &
                      + rxt(k,146)*y(k,68) + rxt(k,147)*y(k,85) + rxt(k,148)*y(k,66) &
                      + rxt(k,149)*y(k,61) + rxt(k,150)*y(k,67) + rxt(k,151)*y(k,82) &
                      + rxt(k,152)*y(k,86) + rxt(k,154)*y(k,62) + rxt(k,155)*y(k,84) &
                      + rxt(k,156)*y(k,81) + rxt(k,250)*y(k,87) + rxt(k,251)*y(k,65) &
                      + rxt(k,252)*y(k,63) + rxt(k,324)*y(k,110) + rxt(k,359)*y(k,55) &
                      + rxt(k,562)*y(k,69) + rxt(k,570)*y(k,48) + rxt(k,572)*y(k,51))
         mat(k,658) = -rxt(k,143)*y(k,99)
         mat(k,299) = -rxt(k,144)*y(k,99)
         mat(k,1307) = -rxt(k,145)*y(k,99)
         mat(k,631) = -rxt(k,146)*y(k,99)
         mat(k,522) = -rxt(k,147)*y(k,99)
         mat(k,314) = -rxt(k,148)*y(k,99)
         mat(k,602) = -rxt(k,149)*y(k,99)
         mat(k,503) = -rxt(k,150)*y(k,99)
         mat(k,372) = -rxt(k,151)*y(k,99)
         mat(k,388) = -rxt(k,152)*y(k,99)
         mat(k,358) = -rxt(k,154)*y(k,99)
         mat(k,481) = -rxt(k,155)*y(k,99)
         mat(k,718) = -rxt(k,156)*y(k,99)
         mat(k,440) = -rxt(k,250)*y(k,99)
         mat(k,792) = -rxt(k,251)*y(k,99)
         mat(k,342) = -rxt(k,252)*y(k,99)
         mat(k,1805) = -rxt(k,324)*y(k,99)
         mat(k,1671) = -rxt(k,359)*y(k,99)
         mat(k,967) = -rxt(k,562)*y(k,99)
         mat(k,273) = -rxt(k,570)*y(k,99)
         mat(k,1009) = -rxt(k,572)*y(k,99)
         mat(k,1092) = rxt(k,573)*y(k,107)
         mat(k,1233) = rxt(k,333)*y(k,103)
         mat(k,1671) = mat(k,1671) + rxt(k,564)*y(k,80) + rxt(k,568)*y(k,94) &
                      + rxt(k,575)*y(k,107) + rxt(k,579)*y(k,108)
         mat(k,120) = rxt(k,564)*y(k,55)
         mat(k,96) = rxt(k,568)*y(k,55)
         mat(k,234) = rxt(k,334)*y(k,103)
         mat(k,1307) = mat(k,1307) + 2.000_r8*rxt(k,145)*y(k,99)
         mat(k,1352) = mat(k,1352) + 2.000_r8*rxt(k,145)*y(k,98)
         mat(k,286) = rxt(k,333)*y(k,54) + rxt(k,334)*y(k,96)
         mat(k,206) = rxt(k,573)*y(k,28) + rxt(k,575)*y(k,55)
         mat(k,69) = rxt(k,579)*y(k,55)
         mat(k,156) = -(rxt(k,321)*y(k,69) + (rxt(k,322) + rxt(k,323)) * y(k,110))
         mat(k,948) = -rxt(k,321)*y(k,100)
         mat(k,1774) = -(rxt(k,322) + rxt(k,323)) * y(k,100)
         mat(k,1319) = rxt(k,324)*y(k,110)
         mat(k,277) = rxt(k,332)*y(k,110)
         mat(k,1774) = mat(k,1774) + rxt(k,324)*y(k,99) + rxt(k,332)*y(k,103)
         mat(k,395) = -((rxt(k,267) + rxt(k,274)) * y(k,51) + (rxt(k,272) + rxt(k,273) &
                      ) * y(k,52) + rxt(k,325)*y(k,28) + rxt(k,326)*y(k,32) + rxt(k,327) &
                      *y(k,56) + (rxt(k,328) + rxt(k,329)) * y(k,54))
         mat(k,985) = -(rxt(k,267) + rxt(k,274)) * y(k,101)
         mat(k,1726) = -(rxt(k,272) + rxt(k,273)) * y(k,101)
         mat(k,1073) = -rxt(k,325)*y(k,101)
         mat(k,1684) = -rxt(k,326)*y(k,101)
         mat(k,1249) = -rxt(k,327)*y(k,101)
         mat(k,1209) = -(rxt(k,328) + rxt(k,329)) * y(k,101)
         mat(k,1209) = mat(k,1209) + rxt(k,331)*y(k,102)
         mat(k,1651) = rxt(k,121)*y(k,65) + rxt(k,357)*y(k,106)
         mat(k,1249) = mat(k,1249) + rxt(k,127)*y(k,68) + rxt(k,318)*y(k,98) &
                      + rxt(k,343)*y(k,105) + rxt(k,358)*y(k,106)
         mat(k,772) = rxt(k,121)*y(k,55)
         mat(k,611) = rxt(k,127)*y(k,56)
         mat(k,1287) = rxt(k,318)*y(k,56)
         mat(k,99) = rxt(k,331)*y(k,54)
         mat(k,1573) = rxt(k,343)*y(k,56)
         mat(k,1606) = rxt(k,357)*y(k,55) + rxt(k,358)*y(k,56)
         mat(k,98) = -(rxt(k,330)*y(k,28) + rxt(k,331)*y(k,54))
         mat(k,1064) = -rxt(k,330)*y(k,102)
         mat(k,1196) = -rxt(k,331)*y(k,102)
         mat(k,1645) = rxt(k,317)*y(k,98)
         mat(k,1283) = rxt(k,317)*y(k,55)
         mat(k,279) = -(rxt(k,332)*y(k,110) + rxt(k,333)*y(k,54) + rxt(k,334)*y(k,96) &
                      + rxt(k,375)*y(k,69))
         mat(k,1777) = -rxt(k,332)*y(k,103)
         mat(k,1207) = -rxt(k,333)*y(k,103)
         mat(k,230) = -rxt(k,334)*y(k,103)
         mat(k,951) = -rxt(k,375)*y(k,103)
         mat(k,1650) = rxt(k,359)*y(k,99)
         mat(k,1321) = rxt(k,359)*y(k,55)
         mat(k,749) = -(rxt(k,306)*y(k,92) + rxt(k,335)*y(k,47) + rxt(k,344)*y(k,51) &
                      + rxt(k,410)*y(k,33) + rxt(k,411)*y(k,35) + rxt(k,412)*y(k,43) &
                      + rxt(k,413)*y(k,54) + rxt(k,414)*y(k,56) + (4._r8*rxt(k,415) &
                      + 4._r8*rxt(k,416)) * y(k,104) + rxt(k,418)*y(k,44) + rxt(k,432) &
                      *y(k,53) + rxt(k,433)*y(k,48) + rxt(k,441)*y(k,52) + rxt(k,442) &
                      *y(k,42) + rxt(k,461)*y(k,25) + (rxt(k,463) + rxt(k,464) &
                      ) * y(k,24) + rxt(k,466)*y(k,40) + rxt(k,469)*y(k,46) + rxt(k,493) &
                      *y(k,2) + rxt(k,495)*y(k,36) + rxt(k,527)*y(k,14) + rxt(k,530) &
                      *y(k,19) + (rxt(k,532) + rxt(k,536)) * y(k,27))
         mat(k,1037) = -rxt(k,306)*y(k,104)
         mat(k,138) = -rxt(k,335)*y(k,104)
         mat(k,995) = -rxt(k,344)*y(k,104)
         mat(k,417) = -rxt(k,410)*y(k,104)
         mat(k,86) = -rxt(k,411)*y(k,104)
         mat(k,681) = -rxt(k,412)*y(k,104)
         mat(k,1219) = -rxt(k,413)*y(k,104)
         mat(k,1257) = -rxt(k,414)*y(k,104)
         mat(k,53) = -rxt(k,418)*y(k,104)
         mat(k,1417) = -rxt(k,432)*y(k,104)
         mat(k,267) = -rxt(k,433)*y(k,104)
         mat(k,1738) = -rxt(k,441)*y(k,104)
         mat(k,1163) = -rxt(k,442)*y(k,104)
         mat(k,212) = -rxt(k,461)*y(k,104)
         mat(k,1502) = -(rxt(k,463) + rxt(k,464)) * y(k,104)
         mat(k,1123) = -rxt(k,466)*y(k,104)
         mat(k,192) = -rxt(k,469)*y(k,104)
         mat(k,455) = -rxt(k,493)*y(k,104)
         mat(k,186) = -rxt(k,495)*y(k,104)
         mat(k,538) = -rxt(k,527)*y(k,104)
         mat(k,42) = -rxt(k,530)*y(k,104)
         mat(k,144) = -(rxt(k,532) + rxt(k,536)) * y(k,104)
         mat(k,538) = mat(k,538) + rxt(k,526)*y(k,54)
         mat(k,42) = mat(k,42) + .300_r8*rxt(k,530)*y(k,104)
         mat(k,1502) = mat(k,1502) + rxt(k,337)*y(k,105)
         mat(k,177) = rxt(k,504)*y(k,110)
         mat(k,1693) = 2.000_r8*rxt(k,404)*y(k,43) + rxt(k,409)*y(k,56) + rxt(k,124) &
                      *y(k,68)
         mat(k,417) = mat(k,417) + rxt(k,401)*y(k,54) + rxt(k,384)*y(k,95)
         mat(k,86) = mat(k,86) + rxt(k,402)*y(k,54)
         mat(k,186) = mat(k,186) + rxt(k,494)*y(k,54) + rxt(k,500)*y(k,95)
         mat(k,1123) = mat(k,1123) + rxt(k,465)*y(k,54) + rxt(k,477)*y(k,95) &
                      + rxt(k,351)*y(k,106)
         mat(k,1163) = mat(k,1163) + rxt(k,119)*y(k,65) + rxt(k,352)*y(k,106)
         mat(k,681) = mat(k,681) + 2.000_r8*rxt(k,404)*y(k,32) + rxt(k,434)*y(k,51) &
                      + rxt(k,429)*y(k,53) + rxt(k,407)*y(k,54) + rxt(k,408)*y(k,56) &
                      + rxt(k,450)*y(k,60)
         mat(k,168) = rxt(k,496)*y(k,54)
         mat(k,192) = mat(k,192) + rxt(k,468)*y(k,54)
         mat(k,995) = mat(k,995) + rxt(k,434)*y(k,43)
         mat(k,1738) = mat(k,1738) + rxt(k,341)*y(k,105)
         mat(k,1417) = mat(k,1417) + rxt(k,429)*y(k,43)
         mat(k,1219) = mat(k,1219) + rxt(k,526)*y(k,14) + rxt(k,401)*y(k,33) &
                      + rxt(k,402)*y(k,35) + rxt(k,494)*y(k,36) + rxt(k,465)*y(k,40) &
                      + rxt(k,407)*y(k,43) + rxt(k,496)*y(k,45) + rxt(k,468)*y(k,46)
         mat(k,1257) = mat(k,1257) + rxt(k,409)*y(k,32) + rxt(k,408)*y(k,43) &
                      + rxt(k,343)*y(k,105)
         mat(k,1381) = rxt(k,450)*y(k,43) + rxt(k,336)*y(k,105)
         mat(k,778) = rxt(k,119)*y(k,42)
         mat(k,617) = rxt(k,124)*y(k,32)
         mat(k,954) = rxt(k,133)*y(k,70)
         mat(k,35) = rxt(k,133)*y(k,69) + rxt(k,134)*y(k,110)
         mat(k,292) = rxt(k,184)*y(k,73) + rxt(k,198)*y(k,75) + rxt(k,167)*y(k,77) &
                      + rxt(k,237)*y(k,88) + rxt(k,219)*y(k,89) + rxt(k,202)*y(k,92) &
                      + rxt(k,144)*y(k,99)
         mat(k,821) = rxt(k,184)*y(k,71)
         mat(k,1459) = rxt(k,198)*y(k,71)
         mat(k,1545) = rxt(k,167)*y(k,71)
         mat(k,867) = rxt(k,237)*y(k,71)
         mat(k,910) = rxt(k,219)*y(k,71)
         mat(k,1037) = mat(k,1037) + rxt(k,202)*y(k,71)
         mat(k,569) = rxt(k,384)*y(k,33) + rxt(k,500)*y(k,36) + rxt(k,477)*y(k,40) &
                      + 2.000_r8*rxt(k,385)*y(k,110)
         mat(k,1338) = rxt(k,144)*y(k,71)
         mat(k,157) = rxt(k,323)*y(k,110)
         mat(k,749) = mat(k,749) + .300_r8*rxt(k,530)*y(k,19)
         mat(k,1578) = rxt(k,337)*y(k,24) + rxt(k,341)*y(k,52) + rxt(k,343)*y(k,56) &
                      + rxt(k,336)*y(k,60)
         mat(k,1614) = rxt(k,351)*y(k,40) + rxt(k,352)*y(k,42) + rxt(k,350)*y(k,110)
         mat(k,1791) = rxt(k,504)*y(k,31) + rxt(k,134)*y(k,70) + 2.000_r8*rxt(k,385) &
                      *y(k,95) + rxt(k,323)*y(k,100) + rxt(k,350)*y(k,106)
         mat(k,1598) = -(rxt(k,336)*y(k,60) + rxt(k,337)*y(k,24) + rxt(k,338)*y(k,28) &
                      + rxt(k,339)*y(k,32) + rxt(k,340)*y(k,40) + rxt(k,341)*y(k,52) &
                      + rxt(k,342)*y(k,54) + rxt(k,343)*y(k,56))
         mat(k,1401) = -rxt(k,336)*y(k,105)
         mat(k,1522) = -rxt(k,337)*y(k,105)
         mat(k,1098) = -rxt(k,338)*y(k,105)
         mat(k,1713) = -rxt(k,339)*y(k,105)
         mat(k,1143) = -rxt(k,340)*y(k,105)
         mat(k,1758) = -rxt(k,341)*y(k,105)
         mat(k,1239) = -rxt(k,342)*y(k,105)
         mat(k,1277) = -rxt(k,343)*y(k,105)
         mat(k,1713) = mat(k,1713) + rxt(k,112)*y(k,65) + rxt(k,282)*y(k,81) &
                      + rxt(k,326)*y(k,101)
         mat(k,426) = rxt(k,349)*y(k,106)
         mat(k,798) = rxt(k,112)*y(k,32)
         mat(k,724) = rxt(k,282)*y(k,32)
         mat(k,405) = rxt(k,326)*y(k,32)
         mat(k,1634) = rxt(k,349)*y(k,33) + rxt(k,350)*y(k,110)
         mat(k,1811) = rxt(k,350)*y(k,106)
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
         mat(k,1635) = -(rxt(k,131)*y(k,51) + rxt(k,345)*y(k,60) + rxt(k,346)*y(k,24) &
                      + rxt(k,347)*y(k,28) + (rxt(k,348) + rxt(k,349)) * y(k,33) &
                      + rxt(k,350)*y(k,110) + rxt(k,351)*y(k,40) + rxt(k,352)*y(k,42) &
                      + rxt(k,354)*y(k,52) + rxt(k,355)*y(k,54) + rxt(k,356)*y(k,96) &
                      + rxt(k,357)*y(k,55) + rxt(k,358)*y(k,56))
         mat(k,1016) = -rxt(k,131)*y(k,106)
         mat(k,1402) = -rxt(k,345)*y(k,106)
         mat(k,1523) = -rxt(k,346)*y(k,106)
         mat(k,1099) = -rxt(k,347)*y(k,106)
         mat(k,427) = -(rxt(k,348) + rxt(k,349)) * y(k,106)
         mat(k,1812) = -rxt(k,350)*y(k,106)
         mat(k,1144) = -rxt(k,351)*y(k,106)
         mat(k,1184) = -rxt(k,352)*y(k,106)
         mat(k,1759) = -rxt(k,354)*y(k,106)
         mat(k,1240) = -rxt(k,355)*y(k,106)
         mat(k,235) = -rxt(k,356)*y(k,106)
         mat(k,1678) = -rxt(k,357)*y(k,106)
         mat(k,1278) = -rxt(k,358)*y(k,106)
         mat(k,1240) = mat(k,1240) + rxt(k,320)*y(k,98)
         mat(k,1278) = mat(k,1278) + rxt(k,129)*y(k,69)
         mat(k,974) = rxt(k,129)*y(k,56)
         mat(k,1314) = rxt(k,320)*y(k,54)
         mat(k,199) = -(rxt(k,573)*y(k,28) + rxt(k,575)*y(k,55))
         mat(k,1068) = -rxt(k,573)*y(k,107)
         mat(k,1647) = -rxt(k,575)*y(k,107)
         mat(k,1202) = rxt(k,566)*y(k,80) + rxt(k,567)*y(k,94) + rxt(k,578)*y(k,108) &
                      + rxt(k,584)*y(k,109)
         mat(k,949) = rxt(k,576)*y(k,108) + rxt(k,581)*y(k,109)
         mat(k,113) = rxt(k,566)*y(k,54)
         mat(k,92) = rxt(k,567)*y(k,54)
         mat(k,67) = rxt(k,578)*y(k,54) + rxt(k,576)*y(k,69)
         mat(k,62) = rxt(k,584)*y(k,54) + rxt(k,581)*y(k,69)
         mat(k,65) = -(rxt(k,576)*y(k,69) + rxt(k,578)*y(k,54) + rxt(k,579)*y(k,55))
         mat(k,939) = -rxt(k,576)*y(k,108)
         mat(k,1191) = -rxt(k,578)*y(k,108)
         mat(k,1642) = -rxt(k,579)*y(k,108)
         mat(k,939) = mat(k,939) + rxt(k,580)*y(k,109)
         mat(k,59) = rxt(k,580)*y(k,69)
         mat(k,58) = -((rxt(k,580) + rxt(k,581)) * y(k,69) + rxt(k,584)*y(k,54))
         mat(k,938) = -(rxt(k,580) + rxt(k,581)) * y(k,109)
         mat(k,1190) = -rxt(k,584)*y(k,109)
         mat(k,1816) = -(rxt(k,102)*y(k,61) + rxt(k,113)*y(k,67) + rxt(k,114)*y(k,65) &
                      + rxt(k,134)*y(k,70) + rxt(k,135)*y(k,74) + rxt(k,138)*y(k,76) &
                      + rxt(k,283)*y(k,81) + rxt(k,289)*y(k,85) + rxt(k,291)*y(k,83) &
                      + rxt(k,300)*y(k,90) + rxt(k,302)*y(k,91) + rxt(k,308)*y(k,93) &
                      + (rxt(k,322) + rxt(k,323)) * y(k,100) + rxt(k,324)*y(k,99) &
                      + rxt(k,332)*y(k,103) + rxt(k,350)*y(k,106) + rxt(k,360)*y(k,88) &
                      + rxt(k,361)*y(k,92) + rxt(k,362)*y(k,89) + rxt(k,367)*y(k,78) &
                      + rxt(k,369)*y(k,72) + rxt(k,371)*y(k,73) + rxt(k,373)*y(k,75) &
                      + rxt(k,385)*y(k,95) + rxt(k,504)*y(k,31))
         mat(k,609) = -rxt(k,102)*y(k,110)
         mat(k,509) = -rxt(k,113)*y(k,110)
         mat(k,803) = -rxt(k,114)*y(k,110)
         mat(k,38) = -rxt(k,134)*y(k,110)
         mat(k,25) = -rxt(k,135)*y(k,110)
         mat(k,29) = -rxt(k,138)*y(k,110)
         mat(k,728) = -rxt(k,283)*y(k,110)
         mat(k,529) = -rxt(k,289)*y(k,110)
         mat(k,666) = -rxt(k,291)*y(k,110)
         mat(k,50) = -rxt(k,300)*y(k,110)
         mat(k,136) = -rxt(k,302)*y(k,110)
         mat(k,129) = -rxt(k,308)*y(k,110)
         mat(k,163) = -(rxt(k,322) + rxt(k,323)) * y(k,110)
         mat(k,1363) = -rxt(k,324)*y(k,110)
         mat(k,290) = -rxt(k,332)*y(k,110)
         mat(k,1639) = -rxt(k,350)*y(k,110)
         mat(k,892) = -rxt(k,360)*y(k,110)
         mat(k,1062) = -rxt(k,361)*y(k,110)
         mat(k,935) = -rxt(k,362)*y(k,110)
         mat(k,152) = -rxt(k,367)*y(k,110)
         mat(k,110) = -rxt(k,369)*y(k,110)
         mat(k,846) = -rxt(k,371)*y(k,110)
         mat(k,1484) = -rxt(k,373)*y(k,110)
         mat(k,588) = -rxt(k,385)*y(k,110)
         mat(k,181) = -rxt(k,504)*y(k,110)
         mat(k,554) = rxt(k,527)*y(k,104)
         mat(k,44) = rxt(k,530)*y(k,104)
         mat(k,1718) = rxt(k,405)*y(k,43) + rxt(k,339)*y(k,105)
         mat(k,431) = rxt(k,410)*y(k,104) + rxt(k,348)*y(k,106)
         mat(k,90) = rxt(k,411)*y(k,104)
         mat(k,189) = rxt(k,495)*y(k,104)
         mat(k,1148) = (rxt(k,549)+rxt(k,554))*y(k,45) + (rxt(k,542)+rxt(k,548) &
                       +rxt(k,553))*y(k,46) + rxt(k,101)*y(k,62) + rxt(k,466)*y(k,104) &
                      + rxt(k,340)*y(k,105)
         mat(k,1188) = rxt(k,290)*y(k,85) + rxt(k,442)*y(k,104)
         mat(k,701) = rxt(k,405)*y(k,32) + rxt(k,412)*y(k,104)
         mat(k,57) = rxt(k,418)*y(k,104)
         mat(k,172) = (rxt(k,549)+rxt(k,554))*y(k,40)
         mat(k,197) = (rxt(k,542)+rxt(k,548)+rxt(k,553))*y(k,40) + rxt(k,469)*y(k,104)
         mat(k,141) = rxt(k,335)*y(k,104)
         mat(k,250) = rxt(k,288)*y(k,84)
         mat(k,1020) = rxt(k,118)*y(k,67)
         mat(k,1763) = rxt(k,115)*y(k,67)
         mat(k,609) = mat(k,609) + 3.000_r8*rxt(k,190)*y(k,73) + 4.000_r8*rxt(k,142) &
                      *y(k,75) + 5.000_r8*rxt(k,172)*y(k,77) + 2.000_r8*rxt(k,225) &
                      *y(k,89) + rxt(k,207)*y(k,92)
         mat(k,364) = rxt(k,101)*y(k,40) + 4.000_r8*rxt(k,194)*y(k,73) &
                      + 5.000_r8*rxt(k,159)*y(k,75) + 6.000_r8*rxt(k,177)*y(k,77) &
                      + rxt(k,247)*y(k,88) + 3.000_r8*rxt(k,229)*y(k,89) &
                      + 2.000_r8*rxt(k,212)*y(k,92) + rxt(k,154)*y(k,99)
         mat(k,348) = 3.000_r8*rxt(k,182)*y(k,73) + (4.000_r8*rxt(k,176) &
                       +4.000_r8*rxt(k,262))*y(k,75) + (5.000_r8*rxt(k,164) &
                       +5.000_r8*rxt(k,263))*y(k,77) + 2.000_r8*rxt(k,217)*y(k,89) &
                      + rxt(k,200)*y(k,92)
         mat(k,803) = mat(k,803) + 3.000_r8*rxt(k,181)*y(k,73) + (4.000_r8*rxt(k,165) &
                       +4.000_r8*rxt(k,253))*y(k,75) + (5.000_r8*rxt(k,163) &
                       +5.000_r8*rxt(k,260))*y(k,77) + 2.000_r8*rxt(k,216)*y(k,89) &
                      + rxt(k,199)*y(k,92)
         mat(k,319) = 5.000_r8*rxt(k,189)*y(k,73) + (6.000_r8*rxt(k,242) &
                       +6.000_r8*rxt(k,266))*y(k,75) + (7.000_r8*rxt(k,171) &
                       +7.000_r8*rxt(k,254))*y(k,77) + 2.000_r8*rxt(k,241)*y(k,88) &
                      + 4.000_r8*rxt(k,224)*y(k,89) + 3.000_r8*rxt(k,206)*y(k,92) &
                      + 2.000_r8*rxt(k,148)*y(k,99)
         mat(k,509) = mat(k,509) + rxt(k,118)*y(k,51) + rxt(k,115)*y(k,52) &
                      + 4.000_r8*rxt(k,191)*y(k,73) + (5.000_r8*rxt(k,153) &
                       +5.000_r8*rxt(k,255))*y(k,75) + (6.000_r8*rxt(k,173) &
                       +6.000_r8*rxt(k,256))*y(k,77) + rxt(k,244)*y(k,88) &
                      + 3.000_r8*rxt(k,226)*y(k,89) + 2.000_r8*rxt(k,208)*y(k,92) &
                      + rxt(k,150)*y(k,99)
         mat(k,641) = 3.000_r8*rxt(k,186)*y(k,73) + 4.000_r8*rxt(k,220)*y(k,75) &
                      + 5.000_r8*rxt(k,169)*y(k,77) + 2.000_r8*rxt(k,222)*y(k,89) &
                      + rxt(k,204)*y(k,92)
         mat(k,978) = rxt(k,133)*y(k,70) + 2.000_r8*rxt(k,377)*y(k,72) &
                      + 3.000_r8*rxt(k,378)*y(k,73) + 4.000_r8*rxt(k,136)*y(k,75) &
                      + 5.000_r8*rxt(k,139)*y(k,77) + rxt(k,376)*y(k,78) &
                      + 2.000_r8*rxt(k,298)*y(k,89) + 3.000_r8*rxt(k,299)*y(k,90) &
                      + rxt(k,303)*y(k,92) + rxt(k,321)*y(k,100)
         mat(k,38) = mat(k,38) + rxt(k,133)*y(k,69)
         mat(k,304) = 3.000_r8*rxt(k,184)*y(k,73) + 4.000_r8*rxt(k,198)*y(k,75) &
                      + 5.000_r8*rxt(k,167)*y(k,77) + 2.000_r8*rxt(k,219)*y(k,89) &
                      + rxt(k,202)*y(k,92)
         mat(k,110) = mat(k,110) + 2.000_r8*rxt(k,377)*y(k,69)
         mat(k,846) = mat(k,846) + 3.000_r8*rxt(k,190)*y(k,61) + 4.000_r8*rxt(k,194) &
                      *y(k,62) + 3.000_r8*rxt(k,182)*y(k,63) + 3.000_r8*rxt(k,181) &
                      *y(k,65) + 5.000_r8*rxt(k,189)*y(k,66) + 4.000_r8*rxt(k,191) &
                      *y(k,67) + 3.000_r8*rxt(k,186)*y(k,68) + 3.000_r8*rxt(k,378) &
                      *y(k,69) + 3.000_r8*rxt(k,184)*y(k,71) + 3.000_r8*rxt(k,196) &
                      *y(k,81) + 4.000_r8*rxt(k,192)*y(k,82) + 3.000_r8*rxt(k,183) &
                      *y(k,83) + 5.000_r8*rxt(k,195)*y(k,84) + 4.000_r8*rxt(k,188) &
                      *y(k,85) + 3.000_r8*rxt(k,193)*y(k,86) + 3.000_r8*rxt(k,180) &
                      *y(k,87) + 3.000_r8*rxt(k,185)*y(k,98)
         mat(k,1484) = mat(k,1484) + 4.000_r8*rxt(k,142)*y(k,61) + 5.000_r8*rxt(k,159) &
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
         mat(k,1570) = 5.000_r8*rxt(k,172)*y(k,61) + 6.000_r8*rxt(k,177)*y(k,62) + ( &
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
         mat(k,152) = mat(k,152) + rxt(k,376)*y(k,69)
         mat(k,728) = mat(k,728) + 3.000_r8*rxt(k,196)*y(k,73) + 4.000_r8*rxt(k,161) &
                      *y(k,75) + 5.000_r8*rxt(k,179)*y(k,77) + 2.000_r8*rxt(k,232) &
                      *y(k,89) + rxt(k,214)*y(k,92)
         mat(k,378) = 4.000_r8*rxt(k,192)*y(k,73) + 5.000_r8*rxt(k,157)*y(k,75) &
                      + 6.000_r8*rxt(k,174)*y(k,77) + rxt(k,245)*y(k,88) &
                      + 3.000_r8*rxt(k,227)*y(k,89) + 2.000_r8*rxt(k,210)*y(k,92) &
                      + rxt(k,151)*y(k,99)
         mat(k,666) = mat(k,666) + 3.000_r8*rxt(k,183)*y(k,73) + (4.000_r8*rxt(k,187) &
                       +4.000_r8*rxt(k,259))*y(k,75) + (5.000_r8*rxt(k,166) &
                       +5.000_r8*rxt(k,261))*y(k,77) + 2.000_r8*rxt(k,218)*y(k,89) &
                      + rxt(k,201)*y(k,92)
         mat(k,488) = rxt(k,288)*y(k,50) + 5.000_r8*rxt(k,195)*y(k,73) &
                      + 6.000_r8*rxt(k,160)*y(k,75) + 7.000_r8*rxt(k,178)*y(k,77) &
                      + 2.000_r8*rxt(k,248)*y(k,88) + 4.000_r8*rxt(k,230)*y(k,89) &
                      + 3.000_r8*rxt(k,213)*y(k,92) + 2.000_r8*rxt(k,155)*y(k,99)
         mat(k,529) = mat(k,529) + rxt(k,290)*y(k,42) + 4.000_r8*rxt(k,188)*y(k,73) + ( &
                      + 5.000_r8*rxt(k,231)+5.000_r8*rxt(k,257))*y(k,75) + ( &
                      + 6.000_r8*rxt(k,170)+6.000_r8*rxt(k,258))*y(k,77) + rxt(k,240) &
                      *y(k,88) + 3.000_r8*rxt(k,223)*y(k,89) + 2.000_r8*rxt(k,205) &
                      *y(k,92) + rxt(k,147)*y(k,99)
         mat(k,394) = 3.000_r8*rxt(k,193)*y(k,73) + 4.000_r8*rxt(k,158)*y(k,75) &
                      + 5.000_r8*rxt(k,175)*y(k,77) + 2.000_r8*rxt(k,228)*y(k,89) &
                      + rxt(k,211)*y(k,92)
         mat(k,445) = 3.000_r8*rxt(k,180)*y(k,73) + (4.000_r8*rxt(k,141) &
                       +4.000_r8*rxt(k,264))*y(k,75) + (5.000_r8*rxt(k,162) &
                       +5.000_r8*rxt(k,265))*y(k,77) + 2.000_r8*rxt(k,215)*y(k,89) &
                      + rxt(k,197)*y(k,92)
         mat(k,892) = mat(k,892) + rxt(k,247)*y(k,62) + 2.000_r8*rxt(k,241)*y(k,66) &
                      + rxt(k,244)*y(k,67) + rxt(k,245)*y(k,82) + 2.000_r8*rxt(k,248) &
                      *y(k,84) + rxt(k,240)*y(k,85)
         mat(k,935) = mat(k,935) + 2.000_r8*rxt(k,225)*y(k,61) + 3.000_r8*rxt(k,229) &
                      *y(k,62) + 2.000_r8*rxt(k,217)*y(k,63) + 2.000_r8*rxt(k,216) &
                      *y(k,65) + 4.000_r8*rxt(k,224)*y(k,66) + 3.000_r8*rxt(k,226) &
                      *y(k,67) + 2.000_r8*rxt(k,222)*y(k,68) + 2.000_r8*rxt(k,298) &
                      *y(k,69) + 2.000_r8*rxt(k,219)*y(k,71) + 2.000_r8*rxt(k,232) &
                      *y(k,81) + 3.000_r8*rxt(k,227)*y(k,82) + 2.000_r8*rxt(k,218) &
                      *y(k,83) + 4.000_r8*rxt(k,230)*y(k,84) + 3.000_r8*rxt(k,223) &
                      *y(k,85) + 2.000_r8*rxt(k,228)*y(k,86) + 2.000_r8*rxt(k,215) &
                      *y(k,87) + 2.000_r8*rxt(k,221)*y(k,98)
         mat(k,50) = mat(k,50) + 3.000_r8*rxt(k,299)*y(k,69)
         mat(k,1062) = mat(k,1062) + rxt(k,207)*y(k,61) + 2.000_r8*rxt(k,212)*y(k,62) &
                      + rxt(k,200)*y(k,63) + rxt(k,199)*y(k,65) + 3.000_r8*rxt(k,206) &
                      *y(k,66) + 2.000_r8*rxt(k,208)*y(k,67) + rxt(k,204)*y(k,68) &
                      + rxt(k,303)*y(k,69) + rxt(k,202)*y(k,71) + rxt(k,214)*y(k,81) &
                      + 2.000_r8*rxt(k,210)*y(k,82) + rxt(k,201)*y(k,83) &
                      + 3.000_r8*rxt(k,213)*y(k,84) + 2.000_r8*rxt(k,205)*y(k,85) &
                      + rxt(k,211)*y(k,86) + rxt(k,197)*y(k,87) + rxt(k,203)*y(k,98)
         mat(k,1318) = 3.000_r8*rxt(k,185)*y(k,73) + 4.000_r8*rxt(k,209)*y(k,75) &
                      + 5.000_r8*rxt(k,168)*y(k,77) + 2.000_r8*rxt(k,221)*y(k,89) &
                      + rxt(k,203)*y(k,92)
         mat(k,1363) = mat(k,1363) + rxt(k,154)*y(k,62) + 2.000_r8*rxt(k,148)*y(k,66) &
                      + rxt(k,150)*y(k,67) + rxt(k,151)*y(k,82) + 2.000_r8*rxt(k,155) &
                      *y(k,84) + rxt(k,147)*y(k,85)
         mat(k,163) = mat(k,163) + rxt(k,321)*y(k,69)
         mat(k,770) = rxt(k,527)*y(k,14) + rxt(k,530)*y(k,19) + rxt(k,410)*y(k,33) &
                      + rxt(k,411)*y(k,35) + rxt(k,495)*y(k,36) + rxt(k,466)*y(k,40) &
                      + rxt(k,442)*y(k,42) + rxt(k,412)*y(k,43) + rxt(k,418)*y(k,44) &
                      + rxt(k,469)*y(k,46) + rxt(k,335)*y(k,47) + 2.000_r8*rxt(k,415) &
                      *y(k,104)
         mat(k,1603) = rxt(k,339)*y(k,32) + rxt(k,340)*y(k,40)
         mat(k,1639) = mat(k,1639) + rxt(k,348)*y(k,33)
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
         mat(k, 54) = lmat(k, 54)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
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
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = lmat(k, 80)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = mat(k, 84) + lmat(k, 84)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 91) = mat(k, 91) + lmat(k, 91)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 103) = lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 106) = lmat(k, 106)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 124) = lmat(k, 124)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 135) = mat(k, 135) + lmat(k, 135)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 160) = lmat(k, 160)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 165) = mat(k, 165) + lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 203) = lmat(k, 203)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 219) = lmat(k, 219)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 249) = lmat(k, 249)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 263) = lmat(k, 263)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 269) = lmat(k, 269)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 335) = lmat(k, 335)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 378) = mat(k, 378) + lmat(k, 378)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 471) = mat(k, 471) + lmat(k, 471)
         mat(k, 472) = lmat(k, 472)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 495) = lmat(k, 495)
         mat(k, 509) = mat(k, 509) + lmat(k, 509)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 514) = lmat(k, 514)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 533) = lmat(k, 533)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 552) = mat(k, 552) + lmat(k, 552)
         mat(k, 555) = mat(k, 555) + lmat(k, 555)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 561) = lmat(k, 561)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 566) = lmat(k, 566)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 572) = lmat(k, 572)
         mat(k, 576) = mat(k, 576) + lmat(k, 576)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 586) = mat(k, 586) + lmat(k, 586)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 630) = lmat(k, 630)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 651) = lmat(k, 651)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 680) = mat(k, 680) + lmat(k, 680)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 711) = lmat(k, 711)
         mat(k, 727) = mat(k, 727) + lmat(k, 727)
         mat(k, 729) = lmat(k, 729)
         mat(k, 730) = lmat(k, 730)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 743) = mat(k, 743) + lmat(k, 743)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 762) = mat(k, 762) + lmat(k, 762)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 779) = mat(k, 779) + lmat(k, 779)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 799) = lmat(k, 799)
         mat(k, 804) = lmat(k, 804)
         mat(k, 823) = mat(k, 823) + lmat(k, 823)
         mat(k, 846) = mat(k, 846) + lmat(k, 846)
         mat(k, 848) = lmat(k, 848)
         mat(k, 870) = mat(k, 870) + lmat(k, 870)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 998) = mat(k, 998) + lmat(k, 998)
         mat(k,1000) = mat(k,1000) + lmat(k,1000)
         mat(k,1001) = mat(k,1001) + lmat(k,1001)
         mat(k,1006) = mat(k,1006) + lmat(k,1006)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1067) = mat(k,1067) + lmat(k,1067)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1089) = lmat(k,1089)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1138) = mat(k,1138) + lmat(k,1138)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1163) = mat(k,1163) + lmat(k,1163)
         mat(k,1173) = mat(k,1173) + lmat(k,1173)
         mat(k,1187) = lmat(k,1187)
         mat(k,1190) = mat(k,1190) + lmat(k,1190)
         mat(k,1191) = mat(k,1191) + lmat(k,1191)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1224) = mat(k,1224) + lmat(k,1224)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1268) = mat(k,1268) + lmat(k,1268)
         mat(k,1269) = mat(k,1269) + lmat(k,1269)
         mat(k,1279) = mat(k,1279) + lmat(k,1279)
         mat(k,1298) = mat(k,1298) + lmat(k,1298)
         mat(k,1306) = mat(k,1306) + lmat(k,1306)
         mat(k,1315) = mat(k,1315) + lmat(k,1315)
         mat(k,1341) = mat(k,1341) + lmat(k,1341)
         mat(k,1344) = mat(k,1344) + lmat(k,1344)
         mat(k,1352) = mat(k,1352) + lmat(k,1352)
         mat(k,1366) = mat(k,1366) + lmat(k,1366)
         mat(k,1370) = lmat(k,1370)
         mat(k,1371) = lmat(k,1371)
         mat(k,1379) = mat(k,1379) + lmat(k,1379)
         mat(k,1390) = mat(k,1390) + lmat(k,1390)
         mat(k,1396) = mat(k,1396) + lmat(k,1396)
         mat(k,1422) = mat(k,1422) + lmat(k,1422)
         mat(k,1426) = mat(k,1426) + lmat(k,1426)
         mat(k,1427) = mat(k,1427) + lmat(k,1427)
         mat(k,1432) = mat(k,1432) + lmat(k,1432)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1440) = mat(k,1440) + lmat(k,1440)
         mat(k,1461) = lmat(k,1461)
         mat(k,1476) = mat(k,1476) + lmat(k,1476)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1513) = mat(k,1513) + lmat(k,1513)
         mat(k,1517) = mat(k,1517) + lmat(k,1517)
         mat(k,1520) = mat(k,1520) + lmat(k,1520)
         mat(k,1562) = lmat(k,1562)
         mat(k,1564) = mat(k,1564) + lmat(k,1564)
         mat(k,1570) = mat(k,1570) + lmat(k,1570)
         mat(k,1578) = mat(k,1578) + lmat(k,1578)
         mat(k,1583) = mat(k,1583) + lmat(k,1583)
         mat(k,1598) = mat(k,1598) + lmat(k,1598)
         mat(k,1619) = mat(k,1619) + lmat(k,1619)
         mat(k,1625) = mat(k,1625) + lmat(k,1625)
         mat(k,1635) = mat(k,1635) + lmat(k,1635)
         mat(k,1641) = lmat(k,1641)
         mat(k,1642) = mat(k,1642) + lmat(k,1642)
         mat(k,1647) = mat(k,1647) + lmat(k,1647)
         mat(k,1652) = mat(k,1652) + lmat(k,1652)
         mat(k,1662) = mat(k,1662) + lmat(k,1662)
         mat(k,1668) = mat(k,1668) + lmat(k,1668)
         mat(k,1671) = mat(k,1671) + lmat(k,1671)
         mat(k,1679) = mat(k,1679) + lmat(k,1679)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1738) = mat(k,1738) + lmat(k,1738)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1748) = mat(k,1748) + lmat(k,1748)
         mat(k,1749) = mat(k,1749) + lmat(k,1749)
         mat(k,1762) = mat(k,1762) + lmat(k,1762)
         mat(k,1781) = lmat(k,1781)
         mat(k,1786) = mat(k,1786) + lmat(k,1786)
         mat(k,1791) = mat(k,1791) + lmat(k,1791)
         mat(k,1802) = lmat(k,1802)
         mat(k,1814) = lmat(k,1814)
         mat(k,1816) = mat(k,1816) + lmat(k,1816)
         mat(k, 115) = 0._r8
         mat(k, 118) = 0._r8
         mat(k, 149) = 0._r8
         mat(k, 154) = 0._r8
         mat(k, 158) = 0._r8
         mat(k, 162) = 0._r8
         mat(k, 171) = 0._r8
         mat(k, 201) = 0._r8
         mat(k, 202) = 0._r8
         mat(k, 208) = 0._r8
         mat(k, 220) = 0._r8
         mat(k, 255) = 0._r8
         mat(k, 260) = 0._r8
         mat(k, 262) = 0._r8
         mat(k, 264) = 0._r8
         mat(k, 266) = 0._r8
         mat(k, 271) = 0._r8
         mat(k, 280) = 0._r8
         mat(k, 281) = 0._r8
         mat(k, 285) = 0._r8
         mat(k, 287) = 0._r8
         mat(k, 289) = 0._r8
         mat(k, 325) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 328) = 0._r8
         mat(k, 329) = 0._r8
         mat(k, 332) = 0._r8
         mat(k, 333) = 0._r8
         mat(k, 412) = 0._r8
         mat(k, 414) = 0._r8
         mat(k, 416) = 0._r8
         mat(k, 419) = 0._r8
         mat(k, 421) = 0._r8
         mat(k, 424) = 0._r8
         mat(k, 425) = 0._r8
         mat(k, 428) = 0._r8
         mat(k, 430) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 453) = 0._r8
         mat(k, 457) = 0._r8
         mat(k, 458) = 0._r8
         mat(k, 459) = 0._r8
         mat(k, 461) = 0._r8
         mat(k, 463) = 0._r8
         mat(k, 466) = 0._r8
         mat(k, 468) = 0._r8
         mat(k, 473) = 0._r8
         mat(k, 480) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 528) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 536) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 541) = 0._r8
         mat(k, 545) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 550) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 570) = 0._r8
         mat(k, 571) = 0._r8
         mat(k, 573) = 0._r8
         mat(k, 575) = 0._r8
         mat(k, 578) = 0._r8
         mat(k, 579) = 0._r8
         mat(k, 581) = 0._r8
         mat(k, 583) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 587) = 0._r8
         mat(k, 614) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 636) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 678) = 0._r8
         mat(k, 679) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 686) = 0._r8
         mat(k, 688) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 696) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 716) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 747) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 765) = 0._r8
         mat(k, 766) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 783) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 794) = 0._r8
         mat(k, 805) = 0._r8
         mat(k, 819) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 827) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 833) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 839) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 841) = 0._r8
         mat(k, 842) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 871) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 881) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 890) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 908) = 0._r8
         mat(k, 912) = 0._r8
         mat(k, 913) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 922) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 928) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 931) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 977) = 0._r8
         mat(k, 988) = 0._r8
         mat(k, 990) = 0._r8
         mat(k, 997) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1005) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1054) = 0._r8
         mat(k,1055) = 0._r8
         mat(k,1056) = 0._r8
         mat(k,1057) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1097) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1124) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1266) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1350) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1387) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1411) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1420) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1477) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1507) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1527) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1548) = 0._r8
         mat(k,1549) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1557) = 0._r8
         mat(k,1559) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1579) = 0._r8
         mat(k,1580) = 0._r8
         mat(k,1581) = 0._r8
         mat(k,1582) = 0._r8
         mat(k,1584) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1592) = 0._r8
         mat(k,1594) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1597) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1605) = 0._r8
         mat(k,1608) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1612) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1617) = 0._r8
         mat(k,1618) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1627) = 0._r8
         mat(k,1628) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1631) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1656) = 0._r8
         mat(k,1657) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1664) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1667) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1673) = 0._r8
         mat(k,1674) = 0._r8
         mat(k,1675) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1682) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1690) = 0._r8
         mat(k,1695) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1703) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1708) = 0._r8
         mat(k,1709) = 0._r8
         mat(k,1710) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1729) = 0._r8
         mat(k,1732) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1789) = 0._r8
         mat(k,1796) = 0._r8
         mat(k,1797) = 0._r8
         mat(k,1800) = 0._r8
         mat(k,1803) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1807) = 0._r8
         mat(k,1809) = 0._r8
         mat(k,1815) = 0._r8
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
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 84) = mat(k, 84) - dti(k)
         mat(k, 91) = mat(k, 91) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 137) = mat(k, 137) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 156) = mat(k, 156) - dti(k)
         mat(k, 165) = mat(k, 165) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 191) = mat(k, 191) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 211) = mat(k, 211) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 322) = mat(k, 322) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 365) = mat(k, 365) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 413) = mat(k, 413) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 471) = mat(k, 471) - dti(k)
         mat(k, 491) = mat(k, 491) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 613) = mat(k, 613) - dti(k)
         mat(k, 646) = mat(k, 646) - dti(k)
         mat(k, 680) = mat(k, 680) - dti(k)
         mat(k, 706) = mat(k, 706) - dti(k)
         mat(k, 749) = mat(k, 749) - dti(k)
         mat(k, 779) = mat(k, 779) - dti(k)
         mat(k, 823) = mat(k, 823) - dti(k)
         mat(k, 870) = mat(k, 870) - dti(k)
         mat(k, 914) = mat(k, 914) - dti(k)
         mat(k, 958) = mat(k, 958) - dti(k)
         mat(k,1001) = mat(k,1001) - dti(k)
         mat(k,1044) = mat(k,1044) - dti(k)
         mat(k,1086) = mat(k,1086) - dti(k)
         mat(k,1132) = mat(k,1132) - dti(k)
         mat(k,1173) = mat(k,1173) - dti(k)
         mat(k,1230) = mat(k,1230) - dti(k)
         mat(k,1269) = mat(k,1269) - dti(k)
         mat(k,1306) = mat(k,1306) - dti(k)
         mat(k,1352) = mat(k,1352) - dti(k)
         mat(k,1396) = mat(k,1396) - dti(k)
         mat(k,1432) = mat(k,1432) - dti(k)
         mat(k,1476) = mat(k,1476) - dti(k)
         mat(k,1520) = mat(k,1520) - dti(k)
         mat(k,1564) = mat(k,1564) - dti(k)
         mat(k,1598) = mat(k,1598) - dti(k)
         mat(k,1635) = mat(k,1635) - dti(k)
         mat(k,1679) = mat(k,1679) - dti(k)
         mat(k,1716) = mat(k,1716) - dti(k)
         mat(k,1762) = mat(k,1762) - dti(k)
         mat(k,1816) = mat(k,1816) - dti(k)
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
