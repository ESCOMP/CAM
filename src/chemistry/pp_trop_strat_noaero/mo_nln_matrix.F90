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
         mat(k,532) = -(rxt(k,346)*y(k,190))
         mat(k,1442) = -rxt(k,346)*y(k,1)
         mat(k,1755) = rxt(k,349)*y(k,162)
         mat(k,805) = rxt(k,349)*y(k,116)
         mat(k,521) = -(rxt(k,350)*y(k,190))
         mat(k,1441) = -rxt(k,350)*y(k,2)
         mat(k,804) = rxt(k,347)*y(k,176)
         mat(k,1635) = rxt(k,347)*y(k,162)
         mat(k,759) = -(rxt(k,429)*y(k,118) + rxt(k,430)*y(k,122) + rxt(k,431) &
                      *y(k,190))
         mat(k,1980) = -rxt(k,429)*y(k,4)
         mat(k,1828) = -rxt(k,430)*y(k,4)
         mat(k,1464) = -rxt(k,431)*y(k,4)
         mat(k,88) = -(rxt(k,388)*y(k,190))
         mat(k,1377) = -rxt(k,388)*y(k,5)
         mat(k,313) = -(rxt(k,391)*y(k,190))
         mat(k,1414) = -rxt(k,391)*y(k,6)
         mat(k,373) = rxt(k,389)*y(k,176)
         mat(k,1616) = rxt(k,389)*y(k,164)
         mat(k,89) = .120_r8*rxt(k,388)*y(k,190)
         mat(k,1378) = .120_r8*rxt(k,388)*y(k,5)
         mat(k,756) = .100_r8*rxt(k,430)*y(k,122)
         mat(k,783) = .100_r8*rxt(k,433)*y(k,122)
         mat(k,1817) = .100_r8*rxt(k,430)*y(k,4) + .100_r8*rxt(k,433)*y(k,105)
         mat(k,1743) = .500_r8*rxt(k,390)*y(k,164) + .200_r8*rxt(k,417)*y(k,196) &
                      + .060_r8*rxt(k,423)*y(k,199)
         mat(k,374) = .500_r8*rxt(k,390)*y(k,116)
         mat(k,581) = .200_r8*rxt(k,417)*y(k,116)
         mat(k,605) = .060_r8*rxt(k,423)*y(k,116)
         mat(k,1736) = .200_r8*rxt(k,417)*y(k,196) + .200_r8*rxt(k,423)*y(k,199)
         mat(k,580) = .200_r8*rxt(k,417)*y(k,116)
         mat(k,603) = .200_r8*rxt(k,423)*y(k,116)
         mat(k,1752) = .200_r8*rxt(k,417)*y(k,196) + .150_r8*rxt(k,423)*y(k,199)
         mat(k,583) = .200_r8*rxt(k,417)*y(k,116)
         mat(k,606) = .150_r8*rxt(k,423)*y(k,116)
         mat(k,1738) = .210_r8*rxt(k,423)*y(k,199)
         mat(k,604) = .210_r8*rxt(k,423)*y(k,116)
         mat(k,154) = -(rxt(k,351)*y(k,190))
         mat(k,1389) = -rxt(k,351)*y(k,13)
         mat(k,755) = .050_r8*rxt(k,430)*y(k,122)
         mat(k,782) = .050_r8*rxt(k,433)*y(k,122)
         mat(k,1816) = .050_r8*rxt(k,430)*y(k,4) + .050_r8*rxt(k,433)*y(k,105)
         mat(k,257) = -(rxt(k,317)*y(k,118) + rxt(k,318)*y(k,190))
         mat(k,1974) = -rxt(k,317)*y(k,14)
         mat(k,1405) = -rxt(k,318)*y(k,14)
         mat(k,1247) = -(rxt(k,200)*y(k,40) + rxt(k,201)*y(k,176) + rxt(k,202) &
                      *y(k,122))
         mat(k,1541) = -rxt(k,200)*y(k,15)
         mat(k,1677) = -rxt(k,201)*y(k,15)
         mat(k,1853) = -rxt(k,202)*y(k,15)
         mat(k,1518) = 4.000_r8*rxt(k,203)*y(k,17) + (rxt(k,204)+rxt(k,205))*y(k,57) &
                      + rxt(k,208)*y(k,116) + rxt(k,211)*y(k,121) + rxt(k,458) &
                      *y(k,136) + rxt(k,212)*y(k,190)
         mat(k,1703) = (rxt(k,204)+rxt(k,205))*y(k,17)
         mat(k,695) = rxt(k,213)*y(k,121) + rxt(k,219)*y(k,189) + rxt(k,214)*y(k,190)
         mat(k,1793) = rxt(k,208)*y(k,17)
         mat(k,1571) = rxt(k,211)*y(k,17) + rxt(k,213)*y(k,76)
         mat(k,1084) = rxt(k,458)*y(k,17)
         mat(k,1340) = rxt(k,219)*y(k,76)
         mat(k,1494) = rxt(k,212)*y(k,17) + rxt(k,214)*y(k,76)
         mat(k,1512) = rxt(k,206)*y(k,57)
         mat(k,1697) = rxt(k,206)*y(k,17)
         mat(k,1913) = (rxt(k,520)+rxt(k,525))*y(k,86)
         mat(k,647) = (rxt(k,520)+rxt(k,525))*y(k,80)
         mat(k,1524) = -(4._r8*rxt(k,203)*y(k,17) + (rxt(k,204) + rxt(k,205) + rxt(k,206) &
                      ) * y(k,57) + rxt(k,207)*y(k,176) + rxt(k,208)*y(k,116) &
                      + rxt(k,209)*y(k,117) + rxt(k,211)*y(k,121) + rxt(k,212) &
                      *y(k,190) + rxt(k,458)*y(k,136))
         mat(k,1709) = -(rxt(k,204) + rxt(k,205) + rxt(k,206)) * y(k,17)
         mat(k,1683) = -rxt(k,207)*y(k,17)
         mat(k,1799) = -rxt(k,208)*y(k,17)
         mat(k,1900) = -rxt(k,209)*y(k,17)
         mat(k,1577) = -rxt(k,211)*y(k,17)
         mat(k,1500) = -rxt(k,212)*y(k,17)
         mat(k,1087) = -rxt(k,458)*y(k,17)
         mat(k,1251) = rxt(k,202)*y(k,122)
         mat(k,418) = rxt(k,210)*y(k,121)
         mat(k,699) = rxt(k,220)*y(k,189)
         mat(k,651) = rxt(k,215)*y(k,121)
         mat(k,1577) = mat(k,1577) + rxt(k,210)*y(k,18) + rxt(k,215)*y(k,86)
         mat(k,1859) = rxt(k,202)*y(k,15)
         mat(k,1346) = rxt(k,220)*y(k,76)
         mat(k,414) = -(rxt(k,210)*y(k,121))
         mat(k,1561) = -rxt(k,210)*y(k,18)
         mat(k,1514) = rxt(k,209)*y(k,117)
         mat(k,1878) = rxt(k,209)*y(k,17)
         mat(k,157) = -(rxt(k,392)*y(k,190))
         mat(k,1390) = -rxt(k,392)*y(k,20)
         mat(k,1733) = rxt(k,395)*y(k,166)
         mat(k,331) = rxt(k,395)*y(k,116)
         mat(k,239) = -(rxt(k,394)*y(k,190))
         mat(k,1403) = -rxt(k,394)*y(k,21)
         mat(k,332) = rxt(k,393)*y(k,176)
         mat(k,1609) = rxt(k,393)*y(k,166)
         mat(k,201) = -(rxt(k,266)*y(k,54) + rxt(k,267)*y(k,190))
         mat(k,1937) = -rxt(k,266)*y(k,22)
         mat(k,1397) = -rxt(k,267)*y(k,22)
         mat(k,442) = -(rxt(k,268)*y(k,54) + rxt(k,269)*y(k,122) + rxt(k,294)*y(k,190))
         mat(k,1938) = -rxt(k,268)*y(k,23)
         mat(k,1821) = -rxt(k,269)*y(k,23)
         mat(k,1431) = -rxt(k,294)*y(k,23)
         mat(k,166) = -(rxt(k,274)*y(k,190))
         mat(k,1392) = -rxt(k,274)*y(k,24)
         mat(k,683) = .800_r8*rxt(k,270)*y(k,167) + .200_r8*rxt(k,271)*y(k,171)
         mat(k,1285) = .200_r8*rxt(k,271)*y(k,167)
         mat(k,213) = -(rxt(k,275)*y(k,190))
         mat(k,1399) = -rxt(k,275)*y(k,25)
         mat(k,684) = rxt(k,272)*y(k,176)
         mat(k,1605) = rxt(k,272)*y(k,167)
         mat(k,195) = -(rxt(k,276)*y(k,54) + rxt(k,277)*y(k,190))
         mat(k,1936) = -rxt(k,276)*y(k,26)
         mat(k,1396) = -rxt(k,277)*y(k,26)
         mat(k,840) = -(rxt(k,297)*y(k,118) + rxt(k,298)*y(k,122) + rxt(k,315) &
                      *y(k,190))
         mat(k,1984) = -rxt(k,297)*y(k,27)
         mat(k,1832) = -rxt(k,298)*y(k,27)
         mat(k,1469) = -rxt(k,315)*y(k,27)
         mat(k,707) = .130_r8*rxt(k,375)*y(k,122)
         mat(k,1832) = mat(k,1832) + .130_r8*rxt(k,375)*y(k,93)
         mat(k,307) = -(rxt(k,302)*y(k,190))
         mat(k,1413) = -rxt(k,302)*y(k,28)
         mat(k,660) = rxt(k,300)*y(k,176)
         mat(k,1615) = rxt(k,300)*y(k,168)
         mat(k,66) = -(rxt(k,303)*y(k,190))
         mat(k,1374) = -rxt(k,303)*y(k,29)
         mat(k,170) = -(rxt(k,398)*y(k,190))
         mat(k,1393) = -rxt(k,398)*y(k,30)
         mat(k,512) = rxt(k,396)*y(k,176)
         mat(k,1603) = rxt(k,396)*y(k,169)
         mat(k,1548) = -(rxt(k,164)*y(k,54) + rxt(k,200)*y(k,15) + rxt(k,244)*y(k,176) &
                      + rxt(k,245)*y(k,118) + rxt(k,246)*y(k,121) + rxt(k,247) &
                      *y(k,190))
         mat(k,1958) = -rxt(k,164)*y(k,40)
         mat(k,1252) = -rxt(k,200)*y(k,40)
         mat(k,1684) = -rxt(k,244)*y(k,40)
         mat(k,2015) = -rxt(k,245)*y(k,40)
         mat(k,1578) = -rxt(k,246)*y(k,40)
         mat(k,1501) = -rxt(k,247)*y(k,40)
         mat(k,539) = .400_r8*rxt(k,346)*y(k,190)
         mat(k,771) = .340_r8*rxt(k,430)*y(k,122)
         mat(k,262) = .500_r8*rxt(k,317)*y(k,118)
         mat(k,447) = rxt(k,269)*y(k,122)
         mat(k,848) = .500_r8*rxt(k,298)*y(k,122)
         mat(k,405) = .500_r8*rxt(k,286)*y(k,190)
         mat(k,658) = rxt(k,252)*y(k,190)
         mat(k,299) = .300_r8*rxt(k,253)*y(k,190)
         mat(k,1710) = rxt(k,171)*y(k,171)
         mat(k,867) = .800_r8*rxt(k,291)*y(k,190)
         mat(k,717) = .910_r8*rxt(k,375)*y(k,122)
         mat(k,479) = .300_r8*rxt(k,366)*y(k,190)
         mat(k,1053) = .800_r8*rxt(k,370)*y(k,171)
         mat(k,1067) = .120_r8*rxt(k,328)*y(k,122)
         mat(k,462) = .500_r8*rxt(k,341)*y(k,190)
         mat(k,798) = .340_r8*rxt(k,433)*y(k,122)
         mat(k,1142) = .600_r8*rxt(k,342)*y(k,122)
         mat(k,1800) = .100_r8*rxt(k,348)*y(k,162) + rxt(k,251)*y(k,171) &
                      + .500_r8*rxt(k,319)*y(k,173) + .500_r8*rxt(k,288)*y(k,175) &
                      + .920_r8*rxt(k,358)*y(k,178) + .250_r8*rxt(k,326)*y(k,182) &
                      + rxt(k,335)*y(k,184) + rxt(k,309)*y(k,192) + rxt(k,313) &
                      *y(k,193) + .340_r8*rxt(k,442)*y(k,194) + .320_r8*rxt(k,447) &
                      *y(k,195) + .250_r8*rxt(k,383)*y(k,198)
         mat(k,2015) = mat(k,2015) + .500_r8*rxt(k,317)*y(k,14) + rxt(k,359)*y(k,178) &
                      + .250_r8*rxt(k,325)*y(k,182) + rxt(k,336)*y(k,184)
         mat(k,1860) = .340_r8*rxt(k,430)*y(k,4) + rxt(k,269)*y(k,23) &
                      + .500_r8*rxt(k,298)*y(k,27) + .910_r8*rxt(k,375)*y(k,93) &
                      + .120_r8*rxt(k,328)*y(k,100) + .340_r8*rxt(k,433)*y(k,105) &
                      + .600_r8*rxt(k,342)*y(k,106)
         mat(k,359) = rxt(k,293)*y(k,190)
         mat(k,891) = .680_r8*rxt(k,451)*y(k,190)
         mat(k,814) = .100_r8*rxt(k,348)*y(k,116)
         mat(k,690) = .700_r8*rxt(k,271)*y(k,171)
         mat(k,666) = rxt(k,299)*y(k,171)
         mat(k,1238) = rxt(k,282)*y(k,171) + rxt(k,355)*y(k,178) + .250_r8*rxt(k,322) &
                      *y(k,182) + rxt(k,331)*y(k,184) + .250_r8*rxt(k,380)*y(k,198)
         mat(k,1323) = rxt(k,171)*y(k,57) + .800_r8*rxt(k,370)*y(k,96) + rxt(k,251) &
                      *y(k,116) + .700_r8*rxt(k,271)*y(k,167) + rxt(k,299)*y(k,168) &
                      + rxt(k,282)*y(k,170) + (4.000_r8*rxt(k,248)+2.000_r8*rxt(k,249)) &
                      *y(k,171) + 1.500_r8*rxt(k,356)*y(k,178) + .750_r8*rxt(k,361) &
                      *y(k,179) + .880_r8*rxt(k,323)*y(k,182) + 2.000_r8*rxt(k,332) &
                      *y(k,184) + .750_r8*rxt(k,435)*y(k,188) + .800_r8*rxt(k,311) &
                      *y(k,193) + .930_r8*rxt(k,440)*y(k,194) + .950_r8*rxt(k,445) &
                      *y(k,195) + .800_r8*rxt(k,381)*y(k,198)
         mat(k,454) = .500_r8*rxt(k,319)*y(k,116)
         mat(k,570) = .500_r8*rxt(k,288)*y(k,116)
         mat(k,1684) = mat(k,1684) + .450_r8*rxt(k,333)*y(k,184) + .150_r8*rxt(k,312) &
                      *y(k,193)
         mat(k,1190) = .920_r8*rxt(k,358)*y(k,116) + rxt(k,359)*y(k,118) + rxt(k,355) &
                      *y(k,170) + 1.500_r8*rxt(k,356)*y(k,171)
         mat(k,1122) = .750_r8*rxt(k,361)*y(k,171)
         mat(k,1164) = .250_r8*rxt(k,326)*y(k,116) + .250_r8*rxt(k,325)*y(k,118) &
                      + .250_r8*rxt(k,322)*y(k,170) + .880_r8*rxt(k,323)*y(k,171)
         mat(k,1208) = rxt(k,335)*y(k,116) + rxt(k,336)*y(k,118) + rxt(k,331)*y(k,170) &
                      + 2.000_r8*rxt(k,332)*y(k,171) + .450_r8*rxt(k,333)*y(k,176) &
                      + 4.000_r8*rxt(k,334)*y(k,184)
         mat(k,986) = .750_r8*rxt(k,435)*y(k,171)
         mat(k,1501) = mat(k,1501) + .400_r8*rxt(k,346)*y(k,1) + .500_r8*rxt(k,286) &
                      *y(k,49) + rxt(k,252)*y(k,50) + .300_r8*rxt(k,253)*y(k,51) &
                      + .800_r8*rxt(k,291)*y(k,69) + .300_r8*rxt(k,366)*y(k,94) &
                      + .500_r8*rxt(k,341)*y(k,104) + rxt(k,293)*y(k,127) &
                      + .680_r8*rxt(k,451)*y(k,151)
         mat(k,633) = rxt(k,309)*y(k,116)
         mat(k,1000) = rxt(k,313)*y(k,116) + .800_r8*rxt(k,311)*y(k,171) &
                      + .150_r8*rxt(k,312)*y(k,176)
         mat(k,967) = .340_r8*rxt(k,442)*y(k,116) + .930_r8*rxt(k,440)*y(k,171)
         mat(k,947) = .320_r8*rxt(k,447)*y(k,116) + .950_r8*rxt(k,445)*y(k,171)
         mat(k,1017) = .250_r8*rxt(k,383)*y(k,116) + .250_r8*rxt(k,380)*y(k,170) &
                      + .800_r8*rxt(k,381)*y(k,171)
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
         mat(k,894) = -(rxt(k,278)*y(k,118) + rxt(k,279)*y(k,190))
         mat(k,1989) = -rxt(k,278)*y(k,43)
         mat(k,1474) = -rxt(k,279)*y(k,43)
         mat(k,536) = .800_r8*rxt(k,346)*y(k,190)
         mat(k,260) = rxt(k,317)*y(k,118)
         mat(k,167) = rxt(k,274)*y(k,190)
         mat(k,215) = .500_r8*rxt(k,275)*y(k,190)
         mat(k,841) = .500_r8*rxt(k,298)*y(k,122)
         mat(k,1131) = .100_r8*rxt(k,342)*y(k,122)
         mat(k,1775) = .400_r8*rxt(k,348)*y(k,162) + rxt(k,273)*y(k,167) &
                      + .270_r8*rxt(k,301)*y(k,168) + rxt(k,319)*y(k,173) + rxt(k,338) &
                      *y(k,186) + rxt(k,309)*y(k,192)
         mat(k,1989) = mat(k,1989) + rxt(k,317)*y(k,14)
         mat(k,1836) = .500_r8*rxt(k,298)*y(k,27) + .100_r8*rxt(k,342)*y(k,106)
         mat(k,810) = .400_r8*rxt(k,348)*y(k,116)
         mat(k,687) = rxt(k,273)*y(k,116) + 3.200_r8*rxt(k,270)*y(k,167) &
                      + .800_r8*rxt(k,271)*y(k,171)
         mat(k,663) = .270_r8*rxt(k,301)*y(k,116)
         mat(k,1301) = .800_r8*rxt(k,271)*y(k,167)
         mat(k,452) = rxt(k,319)*y(k,116)
         mat(k,1659) = .200_r8*rxt(k,337)*y(k,186)
         mat(k,544) = rxt(k,338)*y(k,116) + .200_r8*rxt(k,337)*y(k,176)
         mat(k,1474) = mat(k,1474) + .800_r8*rxt(k,346)*y(k,1) + rxt(k,274)*y(k,24) &
                      + .500_r8*rxt(k,275)*y(k,25)
         mat(k,630) = rxt(k,309)*y(k,116)
         mat(k,57) = -(rxt(k,280)*y(k,190))
         mat(k,1372) = -rxt(k,280)*y(k,45)
         mat(k,818) = -(rxt(k,316)*y(k,190))
         mat(k,1467) = -rxt(k,316)*y(k,46)
         mat(k,535) = .800_r8*rxt(k,346)*y(k,190)
         mat(k,761) = .520_r8*rxt(k,430)*y(k,122)
         mat(k,259) = .500_r8*rxt(k,317)*y(k,118)
         mat(k,788) = .520_r8*rxt(k,433)*y(k,122)
         mat(k,1770) = .250_r8*rxt(k,348)*y(k,162) + .820_r8*rxt(k,301)*y(k,168) &
                      + .500_r8*rxt(k,319)*y(k,173) + .270_r8*rxt(k,442)*y(k,194) &
                      + .040_r8*rxt(k,447)*y(k,195)
         mat(k,1983) = .500_r8*rxt(k,317)*y(k,14)
         mat(k,1831) = .520_r8*rxt(k,430)*y(k,4) + .520_r8*rxt(k,433)*y(k,105)
         mat(k,884) = .500_r8*rxt(k,451)*y(k,190)
         mat(k,809) = .250_r8*rxt(k,348)*y(k,116)
         mat(k,662) = .820_r8*rxt(k,301)*y(k,116) + .820_r8*rxt(k,299)*y(k,171)
         mat(k,1296) = .820_r8*rxt(k,299)*y(k,168) + .150_r8*rxt(k,440)*y(k,194) &
                      + .025_r8*rxt(k,445)*y(k,195)
         mat(k,451) = .500_r8*rxt(k,319)*y(k,116)
         mat(k,1467) = mat(k,1467) + .800_r8*rxt(k,346)*y(k,1) + .500_r8*rxt(k,451) &
                      *y(k,151)
         mat(k,956) = .270_r8*rxt(k,442)*y(k,116) + .150_r8*rxt(k,440)*y(k,171)
         mat(k,934) = .040_r8*rxt(k,447)*y(k,116) + .025_r8*rxt(k,445)*y(k,171)
         mat(k,1072) = -(rxt(k,304)*y(k,118) + rxt(k,305)*y(k,190))
         mat(k,2000) = -rxt(k,304)*y(k,47)
         mat(k,1486) = -rxt(k,305)*y(k,47)
         mat(k,926) = rxt(k,306)*y(k,190)
         mat(k,1061) = .880_r8*rxt(k,328)*y(k,122)
         mat(k,1134) = .500_r8*rxt(k,342)*y(k,122)
         mat(k,1786) = .170_r8*rxt(k,401)*y(k,172) + .050_r8*rxt(k,364)*y(k,179) &
                      + .250_r8*rxt(k,326)*y(k,182) + .170_r8*rxt(k,407)*y(k,185) &
                      + .400_r8*rxt(k,417)*y(k,196) + .250_r8*rxt(k,383)*y(k,198) &
                      + .540_r8*rxt(k,423)*y(k,199) + .510_r8*rxt(k,426)*y(k,201)
         mat(k,2000) = mat(k,2000) + .050_r8*rxt(k,365)*y(k,179) + .250_r8*rxt(k,325) &
                      *y(k,182) + .250_r8*rxt(k,384)*y(k,198)
         mat(k,730) = rxt(k,307)*y(k,190)
         mat(k,1845) = .880_r8*rxt(k,328)*y(k,100) + .500_r8*rxt(k,342)*y(k,106)
         mat(k,1227) = .250_r8*rxt(k,322)*y(k,182) + .250_r8*rxt(k,380)*y(k,198)
         mat(k,1311) = .240_r8*rxt(k,323)*y(k,182) + .500_r8*rxt(k,311)*y(k,193) &
                      + .100_r8*rxt(k,381)*y(k,198)
         mat(k,622) = .170_r8*rxt(k,401)*y(k,116) + .070_r8*rxt(k,400)*y(k,176)
         mat(k,1670) = .070_r8*rxt(k,400)*y(k,172) + .070_r8*rxt(k,406)*y(k,185)
         mat(k,1112) = .050_r8*rxt(k,364)*y(k,116) + .050_r8*rxt(k,365)*y(k,118)
         mat(k,1156) = .250_r8*rxt(k,326)*y(k,116) + .250_r8*rxt(k,325)*y(k,118) &
                      + .250_r8*rxt(k,322)*y(k,170) + .240_r8*rxt(k,323)*y(k,171)
         mat(k,743) = .170_r8*rxt(k,407)*y(k,116) + .070_r8*rxt(k,406)*y(k,176)
         mat(k,1486) = mat(k,1486) + rxt(k,306)*y(k,90) + rxt(k,307)*y(k,119)
         mat(k,996) = .500_r8*rxt(k,311)*y(k,171)
         mat(k,590) = .400_r8*rxt(k,417)*y(k,116)
         mat(k,1012) = .250_r8*rxt(k,383)*y(k,116) + .250_r8*rxt(k,384)*y(k,118) &
                      + .250_r8*rxt(k,380)*y(k,170) + .100_r8*rxt(k,381)*y(k,171)
         mat(k,614) = .540_r8*rxt(k,423)*y(k,116)
         mat(k,385) = .510_r8*rxt(k,426)*y(k,116)
         mat(k,430) = -(rxt(k,285)*y(k,190))
         mat(k,1429) = -rxt(k,285)*y(k,48)
         mat(k,836) = .120_r8*rxt(k,298)*y(k,122)
         mat(k,1820) = .120_r8*rxt(k,298)*y(k,27)
         mat(k,1218) = .100_r8*rxt(k,282)*y(k,171) + .150_r8*rxt(k,283)*y(k,176)
         mat(k,1289) = .100_r8*rxt(k,282)*y(k,170)
         mat(k,1629) = .150_r8*rxt(k,283)*y(k,170) + .150_r8*rxt(k,333)*y(k,184)
         mat(k,1198) = .150_r8*rxt(k,333)*y(k,176)
         mat(k,401) = -(rxt(k,286)*y(k,190))
         mat(k,1426) = -rxt(k,286)*y(k,49)
         mat(k,1217) = .400_r8*rxt(k,283)*y(k,176)
         mat(k,1626) = .400_r8*rxt(k,283)*y(k,170) + .400_r8*rxt(k,333)*y(k,184)
         mat(k,1197) = .400_r8*rxt(k,333)*y(k,176)
         mat(k,656) = -(rxt(k,252)*y(k,190))
         mat(k,1453) = -rxt(k,252)*y(k,50)
         mat(k,1038) = .200_r8*rxt(k,370)*y(k,171)
         mat(k,685) = .300_r8*rxt(k,271)*y(k,171)
         mat(k,1291) = .200_r8*rxt(k,370)*y(k,96) + .300_r8*rxt(k,271)*y(k,167) &
                      + 2.000_r8*rxt(k,249)*y(k,171) + .250_r8*rxt(k,356)*y(k,178) &
                      + .250_r8*rxt(k,361)*y(k,179) + .250_r8*rxt(k,323)*y(k,182) &
                      + .250_r8*rxt(k,435)*y(k,188) + .500_r8*rxt(k,311)*y(k,193) &
                      + .250_r8*rxt(k,440)*y(k,194) + .250_r8*rxt(k,445)*y(k,195) &
                      + .300_r8*rxt(k,381)*y(k,198)
         mat(k,1172) = .250_r8*rxt(k,356)*y(k,171)
         mat(k,1101) = .250_r8*rxt(k,361)*y(k,171)
         mat(k,1150) = .250_r8*rxt(k,323)*y(k,171)
         mat(k,974) = .250_r8*rxt(k,435)*y(k,171)
         mat(k,993) = .500_r8*rxt(k,311)*y(k,171)
         mat(k,955) = .250_r8*rxt(k,440)*y(k,171)
         mat(k,933) = .250_r8*rxt(k,445)*y(k,171)
         mat(k,1006) = .300_r8*rxt(k,381)*y(k,171)
         mat(k,295) = -(rxt(k,253)*y(k,190))
         mat(k,1411) = -rxt(k,253)*y(k,51)
         mat(k,1288) = rxt(k,250)*y(k,176)
         mat(k,1613) = rxt(k,250)*y(k,171)
         mat(k,1966) = -(rxt(k,164)*y(k,40) + rxt(k,166)*y(k,72) + rxt(k,167)*y(k,74) &
                      + (rxt(k,168) + rxt(k,169)) * y(k,176) + rxt(k,170)*y(k,122) &
                      + rxt(k,177)*y(k,58) + rxt(k,186)*y(k,87) + rxt(k,276)*y(k,26))
         mat(k,1556) = -rxt(k,164)*y(k,54)
         mat(k,1032) = -rxt(k,166)*y(k,54)
         mat(k,471) = -rxt(k,167)*y(k,54)
         mat(k,1692) = -(rxt(k,168) + rxt(k,169)) * y(k,54)
         mat(k,1868) = -rxt(k,170)*y(k,54)
         mat(k,833) = -rxt(k,177)*y(k,54)
         mat(k,681) = -rxt(k,186)*y(k,54)
         mat(k,199) = -rxt(k,276)*y(k,54)
         mat(k,1533) = rxt(k,205)*y(k,57)
         mat(k,1718) = rxt(k,205)*y(k,17) + (4.000_r8*rxt(k,172)+2.000_r8*rxt(k,174)) &
                      *y(k,57) + rxt(k,176)*y(k,116) + rxt(k,181)*y(k,121) &
                      + rxt(k,459)*y(k,136) + rxt(k,171)*y(k,171) + rxt(k,182) &
                      *y(k,190)
         mat(k,106) = rxt(k,226)*y(k,189)
         mat(k,1932) = rxt(k,184)*y(k,121) + rxt(k,196)*y(k,189) + rxt(k,185)*y(k,190)
         mat(k,1808) = rxt(k,176)*y(k,57)
         mat(k,1586) = rxt(k,181)*y(k,57) + rxt(k,184)*y(k,80)
         mat(k,1094) = rxt(k,459)*y(k,57)
         mat(k,1331) = rxt(k,171)*y(k,57)
         mat(k,1355) = rxt(k,226)*y(k,63) + rxt(k,196)*y(k,80)
         mat(k,1509) = rxt(k,182)*y(k,57) + rxt(k,185)*y(k,80)
         mat(k,1935) = rxt(k,177)*y(k,58)
         mat(k,1696) = 2.000_r8*rxt(k,173)*y(k,57)
         mat(k,824) = rxt(k,177)*y(k,54) + (rxt(k,518)+rxt(k,523)+rxt(k,528))*y(k,80)
         mat(k,1912) = (rxt(k,518)+rxt(k,523)+rxt(k,528))*y(k,58) + (rxt(k,513) &
                       +rxt(k,519)+rxt(k,524))*y(k,87)
         mat(k,675) = (rxt(k,513)+rxt(k,519)+rxt(k,524))*y(k,80)
         mat(k,1695) = 2.000_r8*rxt(k,198)*y(k,57)
         mat(k,1713) = -(rxt(k,171)*y(k,171) + (4._r8*rxt(k,172) + 4._r8*rxt(k,173) &
                      + 4._r8*rxt(k,174) + 4._r8*rxt(k,198)) * y(k,57) + rxt(k,175) &
                      *y(k,176) + rxt(k,176)*y(k,116) + rxt(k,178)*y(k,117) + rxt(k,181) &
                      *y(k,121) + (rxt(k,182) + rxt(k,183)) * y(k,190) + (rxt(k,204) &
                      + rxt(k,205) + rxt(k,206)) * y(k,17) + rxt(k,459)*y(k,136))
         mat(k,1326) = -rxt(k,171)*y(k,57)
         mat(k,1687) = -rxt(k,175)*y(k,57)
         mat(k,1803) = -rxt(k,176)*y(k,57)
         mat(k,1904) = -rxt(k,178)*y(k,57)
         mat(k,1581) = -rxt(k,181)*y(k,57)
         mat(k,1504) = -(rxt(k,182) + rxt(k,183)) * y(k,57)
         mat(k,1528) = -(rxt(k,204) + rxt(k,205) + rxt(k,206)) * y(k,57)
         mat(k,1090) = -rxt(k,459)*y(k,57)
         mat(k,1961) = rxt(k,186)*y(k,87) + rxt(k,170)*y(k,122) + rxt(k,169)*y(k,176)
         mat(k,830) = rxt(k,179)*y(k,121)
         mat(k,1927) = rxt(k,197)*y(k,189)
         mat(k,679) = rxt(k,186)*y(k,54) + rxt(k,187)*y(k,121) + rxt(k,188)*y(k,190)
         mat(k,1581) = mat(k,1581) + rxt(k,179)*y(k,58) + rxt(k,187)*y(k,87)
         mat(k,1863) = rxt(k,170)*y(k,54)
         mat(k,232) = rxt(k,464)*y(k,136)
         mat(k,1090) = mat(k,1090) + rxt(k,464)*y(k,124)
         mat(k,1687) = mat(k,1687) + rxt(k,169)*y(k,54)
         mat(k,1350) = rxt(k,197)*y(k,80)
         mat(k,1504) = mat(k,1504) + rxt(k,188)*y(k,87)
         mat(k,826) = -(rxt(k,177)*y(k,54) + rxt(k,179)*y(k,121) + rxt(k,180)*y(k,190) &
                      + (rxt(k,518) + rxt(k,523) + rxt(k,528)) * y(k,80))
         mat(k,1945) = -rxt(k,177)*y(k,58)
         mat(k,1567) = -rxt(k,179)*y(k,58)
         mat(k,1468) = -rxt(k,180)*y(k,58)
         mat(k,1916) = -(rxt(k,518) + rxt(k,523) + rxt(k,528)) * y(k,58)
         mat(k,1701) = rxt(k,178)*y(k,117)
         mat(k,1887) = rxt(k,178)*y(k,57)
         mat(k,903) = -((rxt(k,255) + rxt(k,265)) * y(k,190))
         mat(k,1475) = -(rxt(k,255) + rxt(k,265)) * y(k,60)
         mat(k,764) = .230_r8*rxt(k,430)*y(k,122)
         mat(k,1246) = rxt(k,200)*y(k,40)
         mat(k,204) = .350_r8*rxt(k,267)*y(k,190)
         mat(k,445) = .630_r8*rxt(k,269)*y(k,122)
         mat(k,842) = .560_r8*rxt(k,298)*y(k,122)
         mat(k,1539) = rxt(k,200)*y(k,15) + rxt(k,164)*y(k,54) + rxt(k,245)*y(k,118) &
                      + rxt(k,246)*y(k,121) + rxt(k,247)*y(k,190)
         mat(k,1071) = rxt(k,304)*y(k,118) + rxt(k,305)*y(k,190)
         mat(k,1948) = rxt(k,164)*y(k,40)
         mat(k,737) = rxt(k,292)*y(k,190)
         mat(k,708) = .620_r8*rxt(k,375)*y(k,122)
         mat(k,1059) = .650_r8*rxt(k,328)*y(k,122)
         mat(k,791) = .230_r8*rxt(k,433)*y(k,122)
         mat(k,1132) = .560_r8*rxt(k,342)*y(k,122)
         mat(k,1776) = .170_r8*rxt(k,401)*y(k,172) + .220_r8*rxt(k,326)*y(k,182) &
                      + .400_r8*rxt(k,404)*y(k,183) + .350_r8*rxt(k,407)*y(k,185) &
                      + .225_r8*rxt(k,442)*y(k,194) + .250_r8*rxt(k,383)*y(k,198)
         mat(k,1990) = rxt(k,245)*y(k,40) + rxt(k,304)*y(k,47) + .220_r8*rxt(k,325) &
                      *y(k,182) + .500_r8*rxt(k,384)*y(k,198)
         mat(k,1568) = rxt(k,246)*y(k,40) + rxt(k,454)*y(k,125)
         mat(k,1837) = .230_r8*rxt(k,430)*y(k,4) + .630_r8*rxt(k,269)*y(k,23) &
                      + .560_r8*rxt(k,298)*y(k,27) + .620_r8*rxt(k,375)*y(k,93) &
                      + .650_r8*rxt(k,328)*y(k,100) + .230_r8*rxt(k,433)*y(k,105) &
                      + .560_r8*rxt(k,342)*y(k,106)
         mat(k,252) = rxt(k,454)*y(k,121) + rxt(k,455)*y(k,190)
         mat(k,886) = .700_r8*rxt(k,451)*y(k,190)
         mat(k,1222) = .220_r8*rxt(k,322)*y(k,182) + .250_r8*rxt(k,380)*y(k,198)
         mat(k,1302) = .110_r8*rxt(k,323)*y(k,182) + .125_r8*rxt(k,440)*y(k,194) &
                      + .200_r8*rxt(k,381)*y(k,198)
         mat(k,621) = .170_r8*rxt(k,401)*y(k,116) + .070_r8*rxt(k,400)*y(k,176)
         mat(k,1660) = .070_r8*rxt(k,400)*y(k,172) + .160_r8*rxt(k,403)*y(k,183) &
                      + .140_r8*rxt(k,406)*y(k,185)
         mat(k,1152) = .220_r8*rxt(k,326)*y(k,116) + .220_r8*rxt(k,325)*y(k,118) &
                      + .220_r8*rxt(k,322)*y(k,170) + .110_r8*rxt(k,323)*y(k,171)
         mat(k,576) = .400_r8*rxt(k,404)*y(k,116) + .160_r8*rxt(k,403)*y(k,176)
         mat(k,742) = .350_r8*rxt(k,407)*y(k,116) + .140_r8*rxt(k,406)*y(k,176)
         mat(k,1475) = mat(k,1475) + .350_r8*rxt(k,267)*y(k,22) + rxt(k,247)*y(k,40) &
                      + rxt(k,305)*y(k,47) + rxt(k,292)*y(k,70) + rxt(k,455)*y(k,125) &
                      + .700_r8*rxt(k,451)*y(k,151)
         mat(k,959) = .225_r8*rxt(k,442)*y(k,116) + .125_r8*rxt(k,440)*y(k,171)
         mat(k,1009) = .250_r8*rxt(k,383)*y(k,116) + .500_r8*rxt(k,384)*y(k,118) &
                      + .250_r8*rxt(k,380)*y(k,170) + .200_r8*rxt(k,381)*y(k,171)
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
         mat(k,70) = -(rxt(k,225)*y(k,189))
         mat(k,1334) = -rxt(k,225)*y(k,62)
         mat(k,103) = -(rxt(k,226)*y(k,189))
         mat(k,1335) = -rxt(k,226)*y(k,63)
         mat(k,119) = -(rxt(k,399)*y(k,190))
         mat(k,1383) = -rxt(k,399)*y(k,64)
         mat(k,113) = .180_r8*rxt(k,419)*y(k,190)
         mat(k,1383) = mat(k,1383) + .180_r8*rxt(k,419)*y(k,153)
         mat(k,183) = -(rxt(k,452)*y(k,118) + (rxt(k,453) + rxt(k,466)) * y(k,190))
         mat(k,1971) = -rxt(k,452)*y(k,65)
         mat(k,1395) = -(rxt(k,453) + rxt(k,466)) * y(k,65)
         mat(k,565) = rxt(k,287)*y(k,176)
         mat(k,1601) = rxt(k,287)*y(k,175)
         mat(k,639) = -(rxt(k,222)*y(k,72) + rxt(k,223)*y(k,202) + rxt(k,224)*y(k,84))
         mat(k,1023) = -rxt(k,222)*y(k,68)
         mat(k,2029) = -rxt(k,223)*y(k,68)
         mat(k,1258) = -rxt(k,224)*y(k,68)
         mat(k,71) = 2.000_r8*rxt(k,225)*y(k,189)
         mat(k,104) = rxt(k,226)*y(k,189)
         mat(k,1337) = 2.000_r8*rxt(k,225)*y(k,62) + rxt(k,226)*y(k,63)
         mat(k,864) = -(rxt(k,291)*y(k,190))
         mat(k,1471) = -rxt(k,291)*y(k,69)
         mat(k,474) = .700_r8*rxt(k,366)*y(k,190)
         mat(k,436) = .500_r8*rxt(k,367)*y(k,190)
         mat(k,285) = rxt(k,378)*y(k,190)
         mat(k,1772) = .050_r8*rxt(k,364)*y(k,179) + .530_r8*rxt(k,326)*y(k,182) &
                      + .225_r8*rxt(k,442)*y(k,194) + .250_r8*rxt(k,383)*y(k,198)
         mat(k,1986) = .050_r8*rxt(k,365)*y(k,179) + .530_r8*rxt(k,325)*y(k,182) &
                      + .250_r8*rxt(k,384)*y(k,198)
         mat(k,1220) = .530_r8*rxt(k,322)*y(k,182) + .250_r8*rxt(k,380)*y(k,198)
         mat(k,1298) = .260_r8*rxt(k,323)*y(k,182) + .125_r8*rxt(k,440)*y(k,194) &
                      + .100_r8*rxt(k,381)*y(k,198)
         mat(k,1105) = .050_r8*rxt(k,364)*y(k,116) + .050_r8*rxt(k,365)*y(k,118)
         mat(k,1151) = .530_r8*rxt(k,326)*y(k,116) + .530_r8*rxt(k,325)*y(k,118) &
                      + .530_r8*rxt(k,322)*y(k,170) + .260_r8*rxt(k,323)*y(k,171)
         mat(k,1471) = mat(k,1471) + .700_r8*rxt(k,366)*y(k,94) + .500_r8*rxt(k,367) &
                      *y(k,95) + rxt(k,378)*y(k,110)
         mat(k,957) = .225_r8*rxt(k,442)*y(k,116) + .125_r8*rxt(k,440)*y(k,171)
         mat(k,1008) = .250_r8*rxt(k,383)*y(k,116) + .250_r8*rxt(k,384)*y(k,118) &
                      + .250_r8*rxt(k,380)*y(k,170) + .100_r8*rxt(k,381)*y(k,171)
         mat(k,736) = -(rxt(k,292)*y(k,190))
         mat(k,1462) = -rxt(k,292)*y(k,70)
         mat(k,203) = .650_r8*rxt(k,267)*y(k,190)
         mat(k,863) = .200_r8*rxt(k,291)*y(k,190)
         mat(k,871) = rxt(k,379)*y(k,190)
         mat(k,1767) = rxt(k,390)*y(k,164) + .050_r8*rxt(k,364)*y(k,179) &
                      + .400_r8*rxt(k,404)*y(k,183) + .170_r8*rxt(k,407)*y(k,185) &
                      + .700_r8*rxt(k,410)*y(k,191) + .600_r8*rxt(k,417)*y(k,196) &
                      + .250_r8*rxt(k,383)*y(k,198) + .340_r8*rxt(k,423)*y(k,199) &
                      + .170_r8*rxt(k,426)*y(k,201)
         mat(k,1979) = .050_r8*rxt(k,365)*y(k,179) + .250_r8*rxt(k,384)*y(k,198)
         mat(k,377) = rxt(k,390)*y(k,116)
         mat(k,1219) = .250_r8*rxt(k,380)*y(k,198)
         mat(k,1295) = .100_r8*rxt(k,381)*y(k,198)
         mat(k,1652) = .160_r8*rxt(k,403)*y(k,183) + .070_r8*rxt(k,406)*y(k,185)
         mat(k,1103) = .050_r8*rxt(k,364)*y(k,116) + .050_r8*rxt(k,365)*y(k,118)
         mat(k,575) = .400_r8*rxt(k,404)*y(k,116) + .160_r8*rxt(k,403)*y(k,176)
         mat(k,740) = .170_r8*rxt(k,407)*y(k,116) + .070_r8*rxt(k,406)*y(k,176)
         mat(k,1462) = mat(k,1462) + .650_r8*rxt(k,267)*y(k,22) + .200_r8*rxt(k,291) &
                      *y(k,69) + rxt(k,379)*y(k,111)
         mat(k,347) = .700_r8*rxt(k,410)*y(k,116)
         mat(k,587) = .600_r8*rxt(k,417)*y(k,116)
         mat(k,1007) = .250_r8*rxt(k,383)*y(k,116) + .250_r8*rxt(k,384)*y(k,118) &
                      + .250_r8*rxt(k,380)*y(k,170) + .100_r8*rxt(k,381)*y(k,171)
         mat(k,611) = .340_r8*rxt(k,423)*y(k,116)
         mat(k,384) = .170_r8*rxt(k,426)*y(k,116)
         mat(k,1273) = -((rxt(k,124) + rxt(k,125) + rxt(k,126)) * y(k,176) + rxt(k,130) &
                      *y(k,122))
         mat(k,1679) = -(rxt(k,124) + rxt(k,125) + rxt(k,126)) * y(k,71)
         mat(k,1855) = -rxt(k,130)*y(k,71)
         mat(k,1543) = rxt(k,247)*y(k,190)
         mat(k,1953) = rxt(k,166)*y(k,72)
         mat(k,904) = rxt(k,265)*y(k,190)
         mat(k,642) = rxt(k,222)*y(k,72)
         mat(k,1026) = rxt(k,166)*y(k,54) + rxt(k,222)*y(k,68) + rxt(k,122)*y(k,121) &
                      + rxt(k,114)*y(k,189) + rxt(k,131)*y(k,190)
         mat(k,696) = rxt(k,220)*y(k,189)
         mat(k,1919) = rxt(k,197)*y(k,189)
         mat(k,278) = rxt(k,152)*y(k,190)
         mat(k,1573) = rxt(k,122)*y(k,72) + rxt(k,134)*y(k,190)
         mat(k,254) = rxt(k,455)*y(k,190)
         mat(k,392) = rxt(k,460)*y(k,190)
         mat(k,1085) = rxt(k,465)*y(k,190)
         mat(k,1342) = rxt(k,114)*y(k,72) + rxt(k,220)*y(k,76) + rxt(k,197)*y(k,80)
         mat(k,1496) = rxt(k,247)*y(k,40) + rxt(k,265)*y(k,60) + rxt(k,131)*y(k,72) &
                      + rxt(k,152)*y(k,107) + rxt(k,134)*y(k,121) + rxt(k,455) &
                      *y(k,125) + rxt(k,460)*y(k,134) + rxt(k,465)*y(k,136)
         mat(k,1024) = -(rxt(k,114)*y(k,189) + rxt(k,122)*y(k,121) + rxt(k,131) &
                      *y(k,190) + rxt(k,166)*y(k,54) + rxt(k,222)*y(k,68))
         mat(k,1339) = -rxt(k,114)*y(k,72)
         mat(k,1569) = -rxt(k,122)*y(k,72)
         mat(k,1483) = -rxt(k,131)*y(k,72)
         mat(k,1949) = -rxt(k,166)*y(k,72)
         mat(k,640) = -rxt(k,222)*y(k,72)
         mat(k,1271) = rxt(k,124)*y(k,176)
         mat(k,1667) = rxt(k,124)*y(k,71)
         mat(k,466) = -(rxt(k,123)*y(k,121) + rxt(k,132)*y(k,190) + rxt(k,167)*y(k,54))
         mat(k,1562) = -rxt(k,123)*y(k,74)
         mat(k,1434) = -rxt(k,132)*y(k,74)
         mat(k,1939) = -rxt(k,167)*y(k,74)
         mat(k,1630) = 2.000_r8*rxt(k,138)*y(k,176)
         mat(k,1434) = mat(k,1434) + 2.000_r8*rxt(k,137)*y(k,190)
         mat(k,174) = rxt(k,468)*y(k,202)
         mat(k,2026) = rxt(k,468)*y(k,138)
         mat(k,694) = -(rxt(k,213)*y(k,121) + rxt(k,214)*y(k,190) + (rxt(k,219) &
                      + rxt(k,220)) * y(k,189))
         mat(k,1565) = -rxt(k,213)*y(k,76)
         mat(k,1458) = -rxt(k,214)*y(k,76)
         mat(k,1338) = -(rxt(k,219) + rxt(k,220)) * y(k,76)
         mat(k,1245) = rxt(k,200)*y(k,40) + rxt(k,201)*y(k,176)
         mat(k,1538) = rxt(k,200)*y(k,15)
         mat(k,1650) = rxt(k,201)*y(k,15)
         mat(k,1931) = -(rxt(k,184)*y(k,121) + rxt(k,185)*y(k,190) + (rxt(k,196) &
                      + rxt(k,197)) * y(k,189) + (rxt(k,513) + rxt(k,519) + rxt(k,524) &
                      ) * y(k,87) + (rxt(k,518) + rxt(k,523) + rxt(k,528)) * y(k,58) &
                      + (rxt(k,520) + rxt(k,525)) * y(k,86))
         mat(k,1585) = -rxt(k,184)*y(k,80)
         mat(k,1508) = -rxt(k,185)*y(k,80)
         mat(k,1354) = -(rxt(k,196) + rxt(k,197)) * y(k,80)
         mat(k,680) = -(rxt(k,513) + rxt(k,519) + rxt(k,524)) * y(k,80)
         mat(k,832) = -(rxt(k,518) + rxt(k,523) + rxt(k,528)) * y(k,80)
         mat(k,653) = -(rxt(k,520) + rxt(k,525)) * y(k,80)
         mat(k,198) = rxt(k,276)*y(k,54)
         mat(k,1555) = rxt(k,164)*y(k,54)
         mat(k,1965) = rxt(k,276)*y(k,26) + rxt(k,164)*y(k,40) + rxt(k,166)*y(k,72) &
                      + rxt(k,167)*y(k,74) + rxt(k,186)*y(k,87) + rxt(k,168)*y(k,176)
         mat(k,1717) = rxt(k,183)*y(k,190)
         mat(k,1031) = rxt(k,166)*y(k,54)
         mat(k,470) = rxt(k,167)*y(k,54)
         mat(k,680) = mat(k,680) + rxt(k,186)*y(k,54)
         mat(k,1691) = rxt(k,168)*y(k,54)
         mat(k,1508) = mat(k,1508) + rxt(k,183)*y(k,57)
         mat(k,107) = -(rxt(k,256)*y(k,190) + rxt(k,264)*y(k,189))
         mat(k,1381) = -rxt(k,256)*y(k,81)
         mat(k,1336) = -rxt(k,264)*y(k,81)
         mat(k,671) = -(rxt(k,257)*y(k,190))
         mat(k,1455) = -rxt(k,257)*y(k,82)
         mat(k,757) = .050_r8*rxt(k,430)*y(k,122)
         mat(k,202) = .350_r8*rxt(k,267)*y(k,190)
         mat(k,444) = .370_r8*rxt(k,269)*y(k,122)
         mat(k,838) = .120_r8*rxt(k,298)*y(k,122)
         mat(k,705) = .110_r8*rxt(k,375)*y(k,122)
         mat(k,1058) = .330_r8*rxt(k,328)*y(k,122)
         mat(k,784) = .050_r8*rxt(k,433)*y(k,122)
         mat(k,1129) = .120_r8*rxt(k,342)*y(k,122)
         mat(k,1764) = rxt(k,260)*y(k,177)
         mat(k,1824) = .050_r8*rxt(k,430)*y(k,4) + .370_r8*rxt(k,269)*y(k,23) &
                      + .120_r8*rxt(k,298)*y(k,27) + .110_r8*rxt(k,375)*y(k,93) &
                      + .330_r8*rxt(k,328)*y(k,100) + .050_r8*rxt(k,433)*y(k,105) &
                      + .120_r8*rxt(k,342)*y(k,106)
         mat(k,1647) = rxt(k,258)*y(k,177)
         mat(k,340) = rxt(k,260)*y(k,116) + rxt(k,258)*y(k,176)
         mat(k,1455) = mat(k,1455) + .350_r8*rxt(k,267)*y(k,22)
         mat(k,638) = rxt(k,222)*y(k,72) + rxt(k,224)*y(k,84) + rxt(k,223)*y(k,202)
         mat(k,1022) = rxt(k,222)*y(k,68)
         mat(k,1257) = rxt(k,224)*y(k,68)
         mat(k,2027) = rxt(k,223)*y(k,68)
         mat(k,1260) = -(rxt(k,161)*y(k,190) + rxt(k,224)*y(k,68))
         mat(k,1495) = -rxt(k,161)*y(k,84)
         mat(k,641) = -rxt(k,224)*y(k,84)
         mat(k,1542) = rxt(k,245)*y(k,118)
         mat(k,897) = rxt(k,278)*y(k,118)
         mat(k,1074) = rxt(k,304)*y(k,118)
         mat(k,827) = (rxt(k,518)+rxt(k,523)+rxt(k,528))*y(k,80)
         mat(k,185) = rxt(k,452)*y(k,118)
         mat(k,1918) = (rxt(k,518)+rxt(k,523)+rxt(k,528))*y(k,58)
         mat(k,1895) = rxt(k,160)*y(k,190)
         mat(k,2009) = rxt(k,245)*y(k,40) + rxt(k,278)*y(k,43) + rxt(k,304)*y(k,47) &
                      + rxt(k,452)*y(k,65)
         mat(k,1495) = mat(k,1495) + rxt(k,160)*y(k,117)
         mat(k,265) = -(rxt(k,139)*y(k,190))
         mat(k,1406) = -rxt(k,139)*y(k,85)
         mat(k,1873) = rxt(k,158)*y(k,176)
         mat(k,1610) = rxt(k,158)*y(k,117)
         mat(k,648) = -(rxt(k,215)*y(k,121) + (rxt(k,520) + rxt(k,525)) * y(k,80))
         mat(k,1563) = -rxt(k,215)*y(k,86)
         mat(k,1914) = -(rxt(k,520) + rxt(k,525)) * y(k,86)
         mat(k,1515) = rxt(k,207)*y(k,176)
         mat(k,1645) = rxt(k,207)*y(k,17)
         mat(k,676) = -(rxt(k,186)*y(k,54) + rxt(k,187)*y(k,121) + rxt(k,188)*y(k,190) &
                      + (rxt(k,513) + rxt(k,519) + rxt(k,524)) * y(k,80))
         mat(k,1942) = -rxt(k,186)*y(k,87)
         mat(k,1564) = -rxt(k,187)*y(k,87)
         mat(k,1456) = -rxt(k,188)*y(k,87)
         mat(k,1915) = -(rxt(k,513) + rxt(k,519) + rxt(k,524)) * y(k,87)
         mat(k,1699) = rxt(k,175)*y(k,176)
         mat(k,825) = rxt(k,180)*y(k,190)
         mat(k,1648) = rxt(k,175)*y(k,57)
         mat(k,1456) = mat(k,1456) + rxt(k,180)*y(k,58)
         mat(k,912) = -(rxt(k,321)*y(k,190))
         mat(k,1476) = -rxt(k,321)*y(k,88)
         mat(k,475) = .300_r8*rxt(k,366)*y(k,190)
         mat(k,437) = .500_r8*rxt(k,367)*y(k,190)
         mat(k,1777) = rxt(k,320)*y(k,173) + rxt(k,327)*y(k,182)
         mat(k,453) = rxt(k,320)*y(k,116)
         mat(k,1153) = rxt(k,327)*y(k,116)
         mat(k,1476) = mat(k,1476) + .300_r8*rxt(k,366)*y(k,94) + .500_r8*rxt(k,367) &
                      *y(k,95)
         mat(k,149) = -(rxt(k,352)*y(k,190))
         mat(k,1388) = -rxt(k,352)*y(k,89)
         mat(k,925) = -(rxt(k,306)*y(k,190))
         mat(k,1477) = -rxt(k,306)*y(k,90)
         mat(k,476) = .700_r8*rxt(k,366)*y(k,190)
         mat(k,438) = .500_r8*rxt(k,367)*y(k,190)
         mat(k,459) = .500_r8*rxt(k,341)*y(k,190)
         mat(k,1778) = .050_r8*rxt(k,364)*y(k,179) + .220_r8*rxt(k,326)*y(k,182) &
                      + .250_r8*rxt(k,383)*y(k,198)
         mat(k,1992) = .050_r8*rxt(k,365)*y(k,179) + .220_r8*rxt(k,325)*y(k,182) &
                      + .250_r8*rxt(k,384)*y(k,198)
         mat(k,425) = .500_r8*rxt(k,310)*y(k,190)
         mat(k,1223) = .220_r8*rxt(k,322)*y(k,182) + .250_r8*rxt(k,380)*y(k,198)
         mat(k,1303) = .230_r8*rxt(k,323)*y(k,182) + .200_r8*rxt(k,311)*y(k,193) &
                      + .100_r8*rxt(k,381)*y(k,198)
         mat(k,1108) = .050_r8*rxt(k,364)*y(k,116) + .050_r8*rxt(k,365)*y(k,118)
         mat(k,1154) = .220_r8*rxt(k,326)*y(k,116) + .220_r8*rxt(k,325)*y(k,118) &
                      + .220_r8*rxt(k,322)*y(k,170) + .230_r8*rxt(k,323)*y(k,171)
         mat(k,1477) = mat(k,1477) + .700_r8*rxt(k,366)*y(k,94) + .500_r8*rxt(k,367) &
                      *y(k,95) + .500_r8*rxt(k,341)*y(k,104) + .500_r8*rxt(k,310) &
                      *y(k,132)
         mat(k,994) = .200_r8*rxt(k,311)*y(k,171)
         mat(k,1010) = .250_r8*rxt(k,383)*y(k,116) + .250_r8*rxt(k,384)*y(k,118) &
                      + .250_r8*rxt(k,380)*y(k,170) + .100_r8*rxt(k,381)*y(k,171)
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
         mat(k,210) = -(rxt(k,353)*y(k,190))
         mat(k,1398) = -rxt(k,353)*y(k,91)
         mat(k,1737) = .870_r8*rxt(k,364)*y(k,179)
         mat(k,1973) = .950_r8*rxt(k,365)*y(k,179)
         mat(k,1215) = rxt(k,360)*y(k,179)
         mat(k,1286) = .750_r8*rxt(k,361)*y(k,179)
         mat(k,1097) = .870_r8*rxt(k,364)*y(k,116) + .950_r8*rxt(k,365)*y(k,118) &
                      + rxt(k,360)*y(k,170) + .750_r8*rxt(k,361)*y(k,171)
         mat(k,79) = -(rxt(k,354)*y(k,190))
         mat(k,1376) = -rxt(k,354)*y(k,92)
         mat(k,595) = .600_r8*rxt(k,377)*y(k,190)
         mat(k,1376) = mat(k,1376) + .600_r8*rxt(k,377)*y(k,98)
         mat(k,706) = -(rxt(k,368)*y(k,118) + rxt(k,375)*y(k,122) + rxt(k,376) &
                      *y(k,190))
         mat(k,1976) = -rxt(k,368)*y(k,93)
         mat(k,1825) = -rxt(k,375)*y(k,93)
         mat(k,1459) = -rxt(k,376)*y(k,93)
         mat(k,473) = -(rxt(k,366)*y(k,190))
         mat(k,1435) = -rxt(k,366)*y(k,94)
         mat(k,1751) = .080_r8*rxt(k,358)*y(k,178)
         mat(k,1170) = .080_r8*rxt(k,358)*y(k,116)
         mat(k,434) = -(rxt(k,367)*y(k,190))
         mat(k,1430) = -rxt(k,367)*y(k,95)
         mat(k,1749) = .080_r8*rxt(k,364)*y(k,179)
         mat(k,1098) = .080_r8*rxt(k,364)*y(k,116)
         mat(k,1044) = -(rxt(k,369)*y(k,170) + rxt(k,370)*y(k,171) + rxt(k,371) &
                      *y(k,176) + rxt(k,372)*y(k,116) + rxt(k,373)*y(k,118))
         mat(k,1225) = -rxt(k,369)*y(k,96)
         mat(k,1309) = -rxt(k,370)*y(k,96)
         mat(k,1668) = -rxt(k,371)*y(k,96)
         mat(k,1784) = -rxt(k,372)*y(k,96)
         mat(k,1998) = -rxt(k,373)*y(k,96)
         mat(k,709) = rxt(k,368)*y(k,118)
         mat(k,1998) = mat(k,1998) + rxt(k,368)*y(k,93)
         mat(k,271) = -(rxt(k,374)*y(k,190))
         mat(k,1407) = -rxt(k,374)*y(k,97)
         mat(k,1035) = rxt(k,371)*y(k,176)
         mat(k,1611) = rxt(k,371)*y(k,96)
         mat(k,596) = -(rxt(k,377)*y(k,190))
         mat(k,1448) = -rxt(k,377)*y(k,98)
         mat(k,1641) = rxt(k,357)*y(k,178) + rxt(k,362)*y(k,179)
         mat(k,1171) = rxt(k,357)*y(k,176)
         mat(k,1100) = rxt(k,362)*y(k,176)
         mat(k,38) = -(rxt(k,499)*y(k,190))
         mat(k,1368) = -rxt(k,499)*y(k,99)
         mat(k,1060) = -(rxt(k,328)*y(k,122) + rxt(k,329)*y(k,190))
         mat(k,1844) = -rxt(k,328)*y(k,100)
         mat(k,1485) = -rxt(k,329)*y(k,100)
         mat(k,710) = .300_r8*rxt(k,375)*y(k,122)
         mat(k,1785) = .360_r8*rxt(k,358)*y(k,178)
         mat(k,1999) = .400_r8*rxt(k,359)*y(k,178)
         mat(k,1844) = mat(k,1844) + .300_r8*rxt(k,375)*y(k,93)
         mat(k,1226) = .390_r8*rxt(k,355)*y(k,178)
         mat(k,1310) = .310_r8*rxt(k,356)*y(k,178)
         mat(k,1179) = .360_r8*rxt(k,358)*y(k,116) + .400_r8*rxt(k,359)*y(k,118) &
                      + .390_r8*rxt(k,355)*y(k,170) + .310_r8*rxt(k,356)*y(k,171)
         mat(k,218) = -(rxt(k,330)*y(k,190))
         mat(k,1400) = -rxt(k,330)*y(k,101)
         mat(k,1606) = rxt(k,324)*y(k,182)
         mat(k,1149) = rxt(k,324)*y(k,176)
         mat(k,396) = -(rxt(k,339)*y(k,190))
         mat(k,1425) = -rxt(k,339)*y(k,102)
         mat(k,1747) = .800_r8*rxt(k,348)*y(k,162)
         mat(k,803) = .800_r8*rxt(k,348)*y(k,116)
         mat(k,223) = -(rxt(k,340)*y(k,190))
         mat(k,1401) = -rxt(k,340)*y(k,103)
         mat(k,1607) = .800_r8*rxt(k,337)*y(k,186)
         mat(k,542) = .800_r8*rxt(k,337)*y(k,176)
         mat(k,458) = -(rxt(k,341)*y(k,190))
         mat(k,1433) = -rxt(k,341)*y(k,104)
         mat(k,1879) = rxt(k,344)*y(k,184)
         mat(k,1199) = rxt(k,344)*y(k,117)
         mat(k,786) = -(rxt(k,432)*y(k,118) + rxt(k,433)*y(k,122) + rxt(k,434) &
                      *y(k,190))
         mat(k,1981) = -rxt(k,432)*y(k,105)
         mat(k,1829) = -rxt(k,433)*y(k,105)
         mat(k,1465) = -rxt(k,434)*y(k,105)
         mat(k,1135) = -(rxt(k,342)*y(k,122) + rxt(k,343)*y(k,190))
         mat(k,1848) = -rxt(k,342)*y(k,106)
         mat(k,1489) = -rxt(k,343)*y(k,106)
         mat(k,712) = .200_r8*rxt(k,375)*y(k,122)
         mat(k,1788) = .560_r8*rxt(k,358)*y(k,178)
         mat(k,2003) = .600_r8*rxt(k,359)*y(k,178)
         mat(k,1848) = mat(k,1848) + .200_r8*rxt(k,375)*y(k,93)
         mat(k,1229) = .610_r8*rxt(k,355)*y(k,178)
         mat(k,1313) = .440_r8*rxt(k,356)*y(k,178)
         mat(k,1181) = .560_r8*rxt(k,358)*y(k,116) + .600_r8*rxt(k,359)*y(k,118) &
                      + .610_r8*rxt(k,355)*y(k,170) + .440_r8*rxt(k,356)*y(k,171)
         mat(k,277) = -(rxt(k,140)*y(k,116) + (rxt(k,141) + rxt(k,142) + rxt(k,143) &
                      ) * y(k,117) + rxt(k,152)*y(k,190))
         mat(k,1739) = -rxt(k,140)*y(k,107)
         mat(k,1874) = -(rxt(k,141) + rxt(k,142) + rxt(k,143)) * y(k,107)
         mat(k,1408) = -rxt(k,152)*y(k,107)
         mat(k,1872) = rxt(k,159)*y(k,118)
         mat(k,1972) = rxt(k,159)*y(k,117)
         mat(k,283) = -(rxt(k,378)*y(k,190))
         mat(k,1409) = -rxt(k,378)*y(k,110)
         mat(k,1036) = .200_r8*rxt(k,370)*y(k,171)
         mat(k,1287) = .200_r8*rxt(k,370)*y(k,96)
         mat(k,873) = -(rxt(k,379)*y(k,190))
         mat(k,1472) = -rxt(k,379)*y(k,111)
         mat(k,1041) = rxt(k,372)*y(k,116) + rxt(k,373)*y(k,118) + rxt(k,369)*y(k,170) &
                      + .800_r8*rxt(k,370)*y(k,171)
         mat(k,1773) = rxt(k,372)*y(k,96)
         mat(k,1987) = rxt(k,373)*y(k,96)
         mat(k,1221) = rxt(k,369)*y(k,96)
         mat(k,1299) = .800_r8*rxt(k,370)*y(k,96)
         mat(k,63) = -(rxt(k,469)*y(k,190))
         mat(k,1373) = -rxt(k,469)*y(k,112)
         mat(k,1804) = -(rxt(k,140)*y(k,107) + rxt(k,149)*y(k,118) + rxt(k,153) &
                      *y(k,176) + rxt(k,154)*y(k,122) + rxt(k,155)*y(k,121) + rxt(k,176) &
                      *y(k,57) + rxt(k,208)*y(k,17) + rxt(k,251)*y(k,171) + rxt(k,260) &
                      *y(k,177) + rxt(k,273)*y(k,167) + rxt(k,284)*y(k,170) + rxt(k,288) &
                      *y(k,175) + rxt(k,301)*y(k,168) + rxt(k,309)*y(k,192) + rxt(k,313) &
                      *y(k,193) + (rxt(k,319) + rxt(k,320)) * y(k,173) + (rxt(k,326) &
                      + rxt(k,327)) * y(k,182) + rxt(k,335)*y(k,184) + rxt(k,338) &
                      *y(k,186) + (rxt(k,348) + rxt(k,349)) * y(k,162) + rxt(k,358) &
                      *y(k,178) + rxt(k,364)*y(k,179) + rxt(k,372)*y(k,96) + rxt(k,383) &
                      *y(k,198) + rxt(k,387)*y(k,161) + rxt(k,390)*y(k,164) + rxt(k,395) &
                      *y(k,166) + rxt(k,397)*y(k,169) + rxt(k,401)*y(k,172) + rxt(k,404) &
                      *y(k,183) + rxt(k,407)*y(k,185) + rxt(k,410)*y(k,191) + rxt(k,417) &
                      *y(k,196) + rxt(k,423)*y(k,199) + rxt(k,426)*y(k,201) + rxt(k,437) &
                      *y(k,188) + rxt(k,442)*y(k,194) + rxt(k,447)*y(k,195))
         mat(k,281) = -rxt(k,140)*y(k,116)
         mat(k,2019) = -rxt(k,149)*y(k,116)
         mat(k,1688) = -rxt(k,153)*y(k,116)
         mat(k,1864) = -rxt(k,154)*y(k,116)
         mat(k,1582) = -rxt(k,155)*y(k,116)
         mat(k,1714) = -rxt(k,176)*y(k,116)
         mat(k,1529) = -rxt(k,208)*y(k,116)
         mat(k,1327) = -rxt(k,251)*y(k,116)
         mat(k,343) = -rxt(k,260)*y(k,116)
         mat(k,692) = -rxt(k,273)*y(k,116)
         mat(k,1240) = -rxt(k,284)*y(k,116)
         mat(k,572) = -rxt(k,288)*y(k,116)
         mat(k,668) = -rxt(k,301)*y(k,116)
         mat(k,635) = -rxt(k,309)*y(k,116)
         mat(k,1002) = -rxt(k,313)*y(k,116)
         mat(k,456) = -(rxt(k,319) + rxt(k,320)) * y(k,116)
         mat(k,1166) = -(rxt(k,326) + rxt(k,327)) * y(k,116)
         mat(k,1210) = -rxt(k,335)*y(k,116)
         mat(k,548) = -rxt(k,338)*y(k,116)
         mat(k,816) = -(rxt(k,348) + rxt(k,349)) * y(k,116)
         mat(k,1192) = -rxt(k,358)*y(k,116)
         mat(k,1124) = -rxt(k,364)*y(k,116)
         mat(k,1055) = -rxt(k,372)*y(k,116)
         mat(k,1019) = -rxt(k,383)*y(k,116)
         mat(k,412) = -rxt(k,387)*y(k,116)
         mat(k,380) = -rxt(k,390)*y(k,116)
         mat(k,337) = -rxt(k,395)*y(k,116)
         mat(k,517) = -rxt(k,397)*y(k,116)
         mat(k,626) = -rxt(k,401)*y(k,116)
         mat(k,578) = -rxt(k,404)*y(k,116)
         mat(k,747) = -rxt(k,407)*y(k,116)
         mat(k,350) = -rxt(k,410)*y(k,116)
         mat(k,593) = -rxt(k,417)*y(k,116)
         mat(k,618) = -rxt(k,423)*y(k,116)
         mat(k,388) = -rxt(k,426)*y(k,116)
         mat(k,988) = -rxt(k,437)*y(k,116)
         mat(k,969) = -rxt(k,442)*y(k,116)
         mat(k,949) = -rxt(k,447)*y(k,116)
         mat(k,281) = mat(k,281) + 2.000_r8*rxt(k,142)*y(k,117) + rxt(k,152)*y(k,190)
         mat(k,1905) = 2.000_r8*rxt(k,142)*y(k,107) + rxt(k,145)*y(k,121) + rxt(k,461) &
                      *y(k,136)
         mat(k,1582) = mat(k,1582) + rxt(k,145)*y(k,117)
         mat(k,1091) = rxt(k,461)*y(k,117)
         mat(k,1505) = rxt(k,152)*y(k,107)
         mat(k,1907) = -((rxt(k,141) + rxt(k,142) + rxt(k,143)) * y(k,107) + (rxt(k,145) &
                      + rxt(k,147)) * y(k,121) + rxt(k,146)*y(k,122) + rxt(k,158) &
                      *y(k,176) + rxt(k,159)*y(k,118) + rxt(k,160)*y(k,190) + rxt(k,178) &
                      *y(k,57) + rxt(k,209)*y(k,17) + rxt(k,295)*y(k,170) + rxt(k,344) &
                      *y(k,184) + rxt(k,402)*y(k,172) + rxt(k,405)*y(k,183) + rxt(k,408) &
                      *y(k,185) + rxt(k,412)*y(k,129) + rxt(k,415)*y(k,161) + rxt(k,461) &
                      *y(k,136))
         mat(k,282) = -(rxt(k,141) + rxt(k,142) + rxt(k,143)) * y(k,117)
         mat(k,1584) = -(rxt(k,145) + rxt(k,147)) * y(k,117)
         mat(k,1866) = -rxt(k,146)*y(k,117)
         mat(k,1690) = -rxt(k,158)*y(k,117)
         mat(k,2021) = -rxt(k,159)*y(k,117)
         mat(k,1507) = -rxt(k,160)*y(k,117)
         mat(k,1716) = -rxt(k,178)*y(k,117)
         mat(k,1531) = -rxt(k,209)*y(k,117)
         mat(k,1242) = -rxt(k,295)*y(k,117)
         mat(k,1212) = -rxt(k,344)*y(k,117)
         mat(k,627) = -rxt(k,402)*y(k,117)
         mat(k,579) = -rxt(k,405)*y(k,117)
         mat(k,748) = -rxt(k,408)*y(k,117)
         mat(k,365) = -rxt(k,412)*y(k,117)
         mat(k,413) = -rxt(k,415)*y(k,117)
         mat(k,1093) = -rxt(k,461)*y(k,117)
         mat(k,541) = rxt(k,346)*y(k,190)
         mat(k,263) = rxt(k,317)*y(k,118)
         mat(k,1531) = mat(k,1531) + rxt(k,208)*y(k,116)
         mat(k,1716) = mat(k,1716) + rxt(k,176)*y(k,116)
         mat(k,268) = rxt(k,139)*y(k,190)
         mat(k,481) = .700_r8*rxt(k,366)*y(k,190)
         mat(k,1056) = rxt(k,372)*y(k,116) + rxt(k,373)*y(k,118)
         mat(k,1806) = rxt(k,208)*y(k,17) + rxt(k,176)*y(k,57) + rxt(k,372)*y(k,96) &
                      + 2.000_r8*rxt(k,149)*y(k,118) + rxt(k,155)*y(k,121) &
                      + rxt(k,154)*y(k,122) + rxt(k,387)*y(k,161) + rxt(k,348) &
                      *y(k,162) + rxt(k,390)*y(k,164) + rxt(k,395)*y(k,166) &
                      + rxt(k,273)*y(k,167) + rxt(k,301)*y(k,168) + rxt(k,397) &
                      *y(k,169) + rxt(k,284)*y(k,170) + rxt(k,251)*y(k,171) &
                      + rxt(k,401)*y(k,172) + rxt(k,319)*y(k,173) + rxt(k,288) &
                      *y(k,175) + rxt(k,153)*y(k,176) + rxt(k,260)*y(k,177) &
                      + .920_r8*rxt(k,358)*y(k,178) + .920_r8*rxt(k,364)*y(k,179) &
                      + rxt(k,326)*y(k,182) + rxt(k,404)*y(k,183) + rxt(k,335) &
                      *y(k,184) + rxt(k,407)*y(k,185) + rxt(k,338)*y(k,186) &
                      + 1.600_r8*rxt(k,437)*y(k,188) + rxt(k,410)*y(k,191) &
                      + rxt(k,309)*y(k,192) + rxt(k,313)*y(k,193) + .900_r8*rxt(k,442) &
                      *y(k,194) + .800_r8*rxt(k,447)*y(k,195) + rxt(k,417)*y(k,196) &
                      + rxt(k,383)*y(k,198) + rxt(k,423)*y(k,199) + rxt(k,426) &
                      *y(k,201)
         mat(k,2021) = mat(k,2021) + rxt(k,317)*y(k,14) + rxt(k,373)*y(k,96) &
                      + 2.000_r8*rxt(k,149)*y(k,116) + rxt(k,150)*y(k,121) &
                      + rxt(k,148)*y(k,176) + rxt(k,359)*y(k,178) + rxt(k,365) &
                      *y(k,179) + rxt(k,325)*y(k,182) + rxt(k,336)*y(k,184) &
                      + 2.000_r8*rxt(k,438)*y(k,188) + rxt(k,151)*y(k,190) &
                      + rxt(k,384)*y(k,198)
         mat(k,734) = rxt(k,307)*y(k,190)
         mat(k,1584) = mat(k,1584) + rxt(k,155)*y(k,116) + rxt(k,150)*y(k,118)
         mat(k,1866) = mat(k,1866) + rxt(k,154)*y(k,116)
         mat(k,511) = rxt(k,444)*y(k,190)
         mat(k,413) = mat(k,413) + rxt(k,387)*y(k,116)
         mat(k,817) = rxt(k,348)*y(k,116)
         mat(k,381) = rxt(k,390)*y(k,116)
         mat(k,338) = rxt(k,395)*y(k,116)
         mat(k,693) = rxt(k,273)*y(k,116)
         mat(k,669) = rxt(k,301)*y(k,116)
         mat(k,519) = rxt(k,397)*y(k,116)
         mat(k,1242) = mat(k,1242) + rxt(k,284)*y(k,116)
         mat(k,1329) = rxt(k,251)*y(k,116) + .500_r8*rxt(k,435)*y(k,188)
         mat(k,627) = mat(k,627) + rxt(k,401)*y(k,116)
         mat(k,457) = rxt(k,319)*y(k,116)
         mat(k,573) = rxt(k,288)*y(k,116)
         mat(k,1690) = mat(k,1690) + rxt(k,153)*y(k,116) + rxt(k,148)*y(k,118)
         mat(k,344) = rxt(k,260)*y(k,116)
         mat(k,1194) = .920_r8*rxt(k,358)*y(k,116) + rxt(k,359)*y(k,118)
         mat(k,1126) = .920_r8*rxt(k,364)*y(k,116) + rxt(k,365)*y(k,118)
         mat(k,1167) = rxt(k,326)*y(k,116) + rxt(k,325)*y(k,118)
         mat(k,579) = mat(k,579) + rxt(k,404)*y(k,116)
         mat(k,1212) = mat(k,1212) + rxt(k,335)*y(k,116) + rxt(k,336)*y(k,118)
         mat(k,748) = mat(k,748) + rxt(k,407)*y(k,116)
         mat(k,549) = rxt(k,338)*y(k,116)
         mat(k,989) = 1.600_r8*rxt(k,437)*y(k,116) + 2.000_r8*rxt(k,438)*y(k,118) &
                      + .500_r8*rxt(k,435)*y(k,171)
         mat(k,1507) = mat(k,1507) + rxt(k,346)*y(k,1) + rxt(k,139)*y(k,85) &
                      + .700_r8*rxt(k,366)*y(k,94) + rxt(k,151)*y(k,118) + rxt(k,307) &
                      *y(k,119) + rxt(k,444)*y(k,148)
         mat(k,351) = rxt(k,410)*y(k,116)
         mat(k,636) = rxt(k,309)*y(k,116)
         mat(k,1003) = rxt(k,313)*y(k,116)
         mat(k,970) = .900_r8*rxt(k,442)*y(k,116)
         mat(k,950) = .800_r8*rxt(k,447)*y(k,116)
         mat(k,594) = rxt(k,417)*y(k,116)
         mat(k,1020) = rxt(k,383)*y(k,116) + rxt(k,384)*y(k,118)
         mat(k,619) = rxt(k,423)*y(k,116)
         mat(k,389) = rxt(k,426)*y(k,116)
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
         mat(k,2024) = -(rxt(k,148)*y(k,176) + rxt(k,149)*y(k,116) + rxt(k,150) &
                      *y(k,121) + rxt(k,151)*y(k,190) + rxt(k,159)*y(k,117) + rxt(k,245) &
                      *y(k,40) + rxt(k,278)*y(k,43) + rxt(k,297)*y(k,27) + rxt(k,304) &
                      *y(k,47) + rxt(k,317)*y(k,14) + rxt(k,325)*y(k,182) + rxt(k,336) &
                      *y(k,184) + rxt(k,359)*y(k,178) + rxt(k,365)*y(k,179) + rxt(k,368) &
                      *y(k,93) + rxt(k,373)*y(k,96) + rxt(k,384)*y(k,198) + rxt(k,429) &
                      *y(k,4) + rxt(k,432)*y(k,105) + rxt(k,438)*y(k,188) + rxt(k,449) &
                      *y(k,150) + rxt(k,452)*y(k,65))
         mat(k,1693) = -rxt(k,148)*y(k,118)
         mat(k,1809) = -rxt(k,149)*y(k,118)
         mat(k,1587) = -rxt(k,150)*y(k,118)
         mat(k,1510) = -rxt(k,151)*y(k,118)
         mat(k,1910) = -rxt(k,159)*y(k,118)
         mat(k,1557) = -rxt(k,245)*y(k,118)
         mat(k,901) = -rxt(k,278)*y(k,118)
         mat(k,853) = -rxt(k,297)*y(k,118)
         mat(k,1078) = -rxt(k,304)*y(k,118)
         mat(k,264) = -rxt(k,317)*y(k,118)
         mat(k,1168) = -rxt(k,325)*y(k,118)
         mat(k,1213) = -rxt(k,336)*y(k,118)
         mat(k,1195) = -rxt(k,359)*y(k,118)
         mat(k,1127) = -rxt(k,365)*y(k,118)
         mat(k,720) = -rxt(k,368)*y(k,118)
         mat(k,1057) = -rxt(k,373)*y(k,118)
         mat(k,1021) = -rxt(k,384)*y(k,118)
         mat(k,774) = -rxt(k,429)*y(k,118)
         mat(k,801) = -rxt(k,432)*y(k,118)
         mat(k,990) = -rxt(k,438)*y(k,118)
         mat(k,862) = -rxt(k,449)*y(k,118)
         mat(k,188) = -rxt(k,452)*y(k,118)
         mat(k,421) = rxt(k,210)*y(k,121)
         mat(k,1967) = rxt(k,177)*y(k,58)
         mat(k,834) = rxt(k,177)*y(k,54) + rxt(k,179)*y(k,121) + rxt(k,180)*y(k,190)
         mat(k,645) = rxt(k,224)*y(k,84)
         mat(k,1269) = rxt(k,224)*y(k,68) + rxt(k,161)*y(k,190)
         mat(k,465) = .500_r8*rxt(k,341)*y(k,190)
         mat(k,1910) = mat(k,1910) + rxt(k,147)*y(k,121) + rxt(k,146)*y(k,122)
         mat(k,1587) = mat(k,1587) + rxt(k,210)*y(k,18) + rxt(k,179)*y(k,58) &
                      + rxt(k,147)*y(k,117)
         mat(k,1869) = rxt(k,146)*y(k,117)
         mat(k,361) = rxt(k,293)*y(k,190)
         mat(k,1510) = mat(k,1510) + rxt(k,180)*y(k,58) + rxt(k,161)*y(k,84) &
                      + .500_r8*rxt(k,341)*y(k,104) + rxt(k,293)*y(k,127)
         mat(k,729) = -(rxt(k,307)*y(k,190))
         mat(k,1461) = -rxt(k,307)*y(k,119)
         mat(k,839) = rxt(k,297)*y(k,118)
         mat(k,435) = .500_r8*rxt(k,367)*y(k,190)
         mat(k,273) = rxt(k,374)*y(k,190)
         mat(k,284) = rxt(k,378)*y(k,190)
         mat(k,870) = rxt(k,379)*y(k,190)
         mat(k,1978) = rxt(k,297)*y(k,27)
         mat(k,1461) = mat(k,1461) + .500_r8*rxt(k,367)*y(k,95) + rxt(k,374)*y(k,97) &
                      + rxt(k,378)*y(k,110) + rxt(k,379)*y(k,111)
         mat(k,289) = -(rxt(k,439)*y(k,190))
         mat(k,1410) = -rxt(k,439)*y(k,120)
         mat(k,1612) = rxt(k,436)*y(k,188)
         mat(k,972) = rxt(k,436)*y(k,176)
         mat(k,1579) = -(rxt(k,119)*y(k,122) + 4._r8*rxt(k,120)*y(k,121) + rxt(k,122) &
                      *y(k,72) + rxt(k,123)*y(k,74) + rxt(k,128)*y(k,176) + rxt(k,134) &
                      *y(k,190) + (rxt(k,145) + rxt(k,147)) * y(k,117) + rxt(k,150) &
                      *y(k,118) + rxt(k,155)*y(k,116) + rxt(k,179)*y(k,58) + rxt(k,181) &
                      *y(k,57) + rxt(k,184)*y(k,80) + rxt(k,187)*y(k,87) + rxt(k,210) &
                      *y(k,18) + rxt(k,211)*y(k,17) + rxt(k,213)*y(k,76) + rxt(k,215) &
                      *y(k,86) + rxt(k,246)*y(k,40) + rxt(k,454)*y(k,125))
         mat(k,1861) = -rxt(k,119)*y(k,121)
         mat(k,1030) = -rxt(k,122)*y(k,121)
         mat(k,468) = -rxt(k,123)*y(k,121)
         mat(k,1685) = -rxt(k,128)*y(k,121)
         mat(k,1502) = -rxt(k,134)*y(k,121)
         mat(k,1902) = -(rxt(k,145) + rxt(k,147)) * y(k,121)
         mat(k,2016) = -rxt(k,150)*y(k,121)
         mat(k,1801) = -rxt(k,155)*y(k,121)
         mat(k,829) = -rxt(k,179)*y(k,121)
         mat(k,1711) = -rxt(k,181)*y(k,121)
         mat(k,1925) = -rxt(k,184)*y(k,121)
         mat(k,678) = -rxt(k,187)*y(k,121)
         mat(k,419) = -rxt(k,210)*y(k,121)
         mat(k,1526) = -rxt(k,211)*y(k,121)
         mat(k,700) = -rxt(k,213)*y(k,121)
         mat(k,652) = -rxt(k,215)*y(k,121)
         mat(k,1549) = -rxt(k,246)*y(k,121)
         mat(k,256) = -rxt(k,454)*y(k,121)
         mat(k,1277) = rxt(k,126)*y(k,176)
         mat(k,280) = rxt(k,140)*y(k,116) + rxt(k,141)*y(k,117)
         mat(k,1801) = mat(k,1801) + rxt(k,140)*y(k,107)
         mat(k,1902) = mat(k,1902) + rxt(k,141)*y(k,107)
         mat(k,1685) = mat(k,1685) + rxt(k,126)*y(k,71)
         mat(k,1502) = mat(k,1502) + 2.000_r8*rxt(k,136)*y(k,190)
         mat(k,1865) = -(rxt(k,118)*y(k,189) + rxt(k,119)*y(k,121) + rxt(k,129) &
                      *y(k,176) + rxt(k,130)*y(k,71) + rxt(k,135)*y(k,190) + rxt(k,146) &
                      *y(k,117) + rxt(k,154)*y(k,116) + rxt(k,170)*y(k,54) + rxt(k,202) &
                      *y(k,15) + rxt(k,269)*y(k,23) + rxt(k,298)*y(k,27) + rxt(k,328) &
                      *y(k,100) + rxt(k,342)*y(k,106) + rxt(k,375)*y(k,93) + rxt(k,413) &
                      *y(k,129) + rxt(k,430)*y(k,4) + rxt(k,433)*y(k,105) + rxt(k,457) &
                      *y(k,134) + rxt(k,463)*y(k,136))
         mat(k,1352) = -rxt(k,118)*y(k,122)
         mat(k,1583) = -rxt(k,119)*y(k,122)
         mat(k,1689) = -rxt(k,129)*y(k,122)
         mat(k,1279) = -rxt(k,130)*y(k,122)
         mat(k,1506) = -rxt(k,135)*y(k,122)
         mat(k,1906) = -rxt(k,146)*y(k,122)
         mat(k,1805) = -rxt(k,154)*y(k,122)
         mat(k,1963) = -rxt(k,170)*y(k,122)
         mat(k,1255) = -rxt(k,202)*y(k,122)
         mat(k,449) = -rxt(k,269)*y(k,122)
         mat(k,851) = -rxt(k,298)*y(k,122)
         mat(k,1069) = -rxt(k,328)*y(k,122)
         mat(k,1145) = -rxt(k,342)*y(k,122)
         mat(k,719) = -rxt(k,375)*y(k,122)
         mat(k,364) = -rxt(k,413)*y(k,122)
         mat(k,773) = -rxt(k,430)*y(k,122)
         mat(k,800) = -rxt(k,433)*y(k,122)
         mat(k,395) = -rxt(k,457)*y(k,122)
         mat(k,1092) = -rxt(k,463)*y(k,122)
         mat(k,1241) = .150_r8*rxt(k,283)*y(k,176)
         mat(k,1689) = mat(k,1689) + .150_r8*rxt(k,283)*y(k,170) + .150_r8*rxt(k,333) &
                      *y(k,184)
         mat(k,1211) = .150_r8*rxt(k,333)*y(k,176)
         mat(k,228) = -(rxt(k,464)*y(k,136))
         mat(k,1080) = -rxt(k,464)*y(k,124)
         mat(k,1513) = rxt(k,204)*y(k,57)
         mat(k,1698) = rxt(k,204)*y(k,17) + 2.000_r8*rxt(k,174)*y(k,57)
         mat(k,249) = -(rxt(k,454)*y(k,121) + rxt(k,455)*y(k,190))
         mat(k,1559) = -rxt(k,454)*y(k,125)
         mat(k,1404) = -rxt(k,455)*y(k,125)
         mat(k,907) = rxt(k,321)*y(k,190)
         mat(k,1734) = .100_r8*rxt(k,442)*y(k,194)
         mat(k,1391) = rxt(k,321)*y(k,88)
         mat(k,953) = .100_r8*rxt(k,442)*y(k,116)
         mat(k,355) = -(rxt(k,293)*y(k,190))
         mat(k,1419) = -rxt(k,293)*y(k,127)
         mat(k,1875) = rxt(k,295)*y(k,170)
         mat(k,1216) = rxt(k,295)*y(k,117)
         mat(k,1871) = rxt(k,415)*y(k,161)
         mat(k,407) = rxt(k,415)*y(k,117)
         mat(k,362) = -(rxt(k,412)*y(k,117) + rxt(k,413)*y(k,122))
         mat(k,1876) = -rxt(k,412)*y(k,129)
         mat(k,1818) = -rxt(k,413)*y(k,129)
         mat(k,121) = .070_r8*rxt(k,399)*y(k,190)
         mat(k,1744) = rxt(k,397)*y(k,169)
         mat(k,100) = .060_r8*rxt(k,411)*y(k,190)
         mat(k,142) = .070_r8*rxt(k,427)*y(k,190)
         mat(k,513) = rxt(k,397)*y(k,116)
         mat(k,1420) = .070_r8*rxt(k,399)*y(k,64) + .060_r8*rxt(k,411)*y(k,130) &
                      + .070_r8*rxt(k,427)*y(k,157)
         mat(k,98) = -(rxt(k,411)*y(k,190))
         mat(k,1379) = -rxt(k,411)*y(k,130)
         mat(k,90) = .530_r8*rxt(k,388)*y(k,190)
         mat(k,1379) = mat(k,1379) + .530_r8*rxt(k,388)*y(k,5)
         mat(k,233) = -(rxt(k,414)*y(k,190))
         mat(k,1402) = -rxt(k,414)*y(k,131)
         mat(k,1608) = rxt(k,409)*y(k,191)
         mat(k,345) = rxt(k,409)*y(k,176)
         mat(k,422) = -(rxt(k,310)*y(k,190))
         mat(k,1428) = -rxt(k,310)*y(k,132)
         mat(k,1628) = rxt(k,308)*y(k,192)
         mat(k,628) = rxt(k,308)*y(k,176)
         mat(k,301) = -(rxt(k,314)*y(k,190))
         mat(k,1412) = -rxt(k,314)*y(k,133)
         mat(k,1614) = .850_r8*rxt(k,312)*y(k,193)
         mat(k,992) = .850_r8*rxt(k,312)*y(k,176)
         mat(k,390) = -(rxt(k,457)*y(k,122) + rxt(k,460)*y(k,190))
         mat(k,1819) = -rxt(k,457)*y(k,134)
         mat(k,1424) = -rxt(k,460)*y(k,134)
         mat(k,1083) = -(rxt(k,458)*y(k,17) + rxt(k,459)*y(k,57) + rxt(k,461)*y(k,117) &
                      + rxt(k,463)*y(k,122) + rxt(k,464)*y(k,124) + rxt(k,465) &
                      *y(k,190))
         mat(k,1517) = -rxt(k,458)*y(k,136)
         mat(k,1702) = -rxt(k,459)*y(k,136)
         mat(k,1891) = -rxt(k,461)*y(k,136)
         mat(k,1846) = -rxt(k,463)*y(k,136)
         mat(k,230) = -rxt(k,464)*y(k,136)
         mat(k,1487) = -rxt(k,465)*y(k,136)
         mat(k,1570) = rxt(k,454)*y(k,125)
         mat(k,1846) = mat(k,1846) + rxt(k,457)*y(k,134)
         mat(k,253) = rxt(k,454)*y(k,121)
         mat(k,391) = rxt(k,457)*y(k,122) + rxt(k,460)*y(k,190)
         mat(k,1487) = mat(k,1487) + rxt(k,460)*y(k,134)
         mat(k,723) = -(rxt(k,467)*y(k,190))
         mat(k,1460) = -rxt(k,467)*y(k,137)
         mat(k,1516) = rxt(k,458)*y(k,136)
         mat(k,1700) = rxt(k,459)*y(k,136)
         mat(k,184) = rxt(k,452)*y(k,118) + (rxt(k,453)+.500_r8*rxt(k,466))*y(k,190)
         mat(k,1884) = rxt(k,461)*y(k,136)
         mat(k,1977) = rxt(k,452)*y(k,65)
         mat(k,1826) = rxt(k,463)*y(k,136)
         mat(k,229) = rxt(k,464)*y(k,136)
         mat(k,251) = rxt(k,455)*y(k,190)
         mat(k,1082) = rxt(k,458)*y(k,17) + rxt(k,459)*y(k,57) + rxt(k,461)*y(k,117) &
                      + rxt(k,463)*y(k,122) + rxt(k,464)*y(k,124) + rxt(k,465) &
                      *y(k,190)
         mat(k,1460) = mat(k,1460) + (rxt(k,453)+.500_r8*rxt(k,466))*y(k,65) &
                      + rxt(k,455)*y(k,125) + rxt(k,465)*y(k,136)
         mat(k,175) = -(rxt(k,468)*y(k,202))
         mat(k,2028) = -rxt(k,468)*y(k,138)
         mat(k,722) = rxt(k,467)*y(k,190)
         mat(k,1394) = rxt(k,467)*y(k,137)
         mat(k,749) = .2202005_r8*rxt(k,487)*y(k,122)
         mat(k,776) = .0508005_r8*rxt(k,503)*y(k,122)
         mat(k,1721) = .1279005_r8*rxt(k,486)*y(k,163) + .0097005_r8*rxt(k,491) &
                      *y(k,165) + .0003005_r8*rxt(k,494)*y(k,180) &
                      + .1056005_r8*rxt(k,498)*y(k,181) + .0245005_r8*rxt(k,502) &
                      *y(k,187) + .0154005_r8*rxt(k,508)*y(k,197) &
                      + .0063005_r8*rxt(k,511)*y(k,200)
         mat(k,1811) = .2202005_r8*rxt(k,487)*y(k,4) + .0508005_r8*rxt(k,503)*y(k,105)
         mat(k,7) = .5931005_r8*rxt(k,505)*y(k,190)
         mat(k,13) = .1279005_r8*rxt(k,486)*y(k,116) + .2202005_r8*rxt(k,485)*y(k,176)
         mat(k,19) = .0097005_r8*rxt(k,491)*y(k,116) + .0023005_r8*rxt(k,490)*y(k,176)
         mat(k,1589) = .2202005_r8*rxt(k,485)*y(k,163) + .0023005_r8*rxt(k,490) &
                      *y(k,165) + .0031005_r8*rxt(k,493)*y(k,180) &
                      + .2381005_r8*rxt(k,497)*y(k,181) + .0508005_r8*rxt(k,501) &
                      *y(k,187) + .1364005_r8*rxt(k,507)*y(k,197) &
                      + .1677005_r8*rxt(k,510)*y(k,200)
         mat(k,25) = .0003005_r8*rxt(k,494)*y(k,116) + .0031005_r8*rxt(k,493)*y(k,176)
         mat(k,31) = .1056005_r8*rxt(k,498)*y(k,116) + .2381005_r8*rxt(k,497)*y(k,176)
         mat(k,39) = .0245005_r8*rxt(k,502)*y(k,116) + .0508005_r8*rxt(k,501)*y(k,176)
         mat(k,1358) = .5931005_r8*rxt(k,505)*y(k,145)
         mat(k,45) = .0154005_r8*rxt(k,508)*y(k,116) + .1364005_r8*rxt(k,507)*y(k,176)
         mat(k,51) = .0063005_r8*rxt(k,511)*y(k,116) + .1677005_r8*rxt(k,510)*y(k,176)
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
         mat(k,750) = .2067005_r8*rxt(k,487)*y(k,122)
         mat(k,777) = .1149005_r8*rxt(k,503)*y(k,122)
         mat(k,1722) = .1792005_r8*rxt(k,486)*y(k,163) + .0034005_r8*rxt(k,491) &
                      *y(k,165) + .0003005_r8*rxt(k,494)*y(k,180) &
                      + .1026005_r8*rxt(k,498)*y(k,181) + .0082005_r8*rxt(k,502) &
                      *y(k,187) + .0452005_r8*rxt(k,508)*y(k,197) &
                      + .0237005_r8*rxt(k,511)*y(k,200)
         mat(k,1812) = .2067005_r8*rxt(k,487)*y(k,4) + .1149005_r8*rxt(k,503)*y(k,105)
         mat(k,8) = .1534005_r8*rxt(k,505)*y(k,190)
         mat(k,14) = .1792005_r8*rxt(k,486)*y(k,116) + .2067005_r8*rxt(k,485)*y(k,176)
         mat(k,20) = .0034005_r8*rxt(k,491)*y(k,116) + .0008005_r8*rxt(k,490)*y(k,176)
         mat(k,1590) = .2067005_r8*rxt(k,485)*y(k,163) + .0008005_r8*rxt(k,490) &
                      *y(k,165) + .0035005_r8*rxt(k,493)*y(k,180) &
                      + .1308005_r8*rxt(k,497)*y(k,181) + .1149005_r8*rxt(k,501) &
                      *y(k,187) + .0101005_r8*rxt(k,507)*y(k,197) &
                      + .0174005_r8*rxt(k,510)*y(k,200)
         mat(k,26) = .0003005_r8*rxt(k,494)*y(k,116) + .0035005_r8*rxt(k,493)*y(k,176)
         mat(k,32) = .1026005_r8*rxt(k,498)*y(k,116) + .1308005_r8*rxt(k,497)*y(k,176)
         mat(k,40) = .0082005_r8*rxt(k,502)*y(k,116) + .1149005_r8*rxt(k,501)*y(k,176)
         mat(k,1359) = .1534005_r8*rxt(k,505)*y(k,145)
         mat(k,46) = .0452005_r8*rxt(k,508)*y(k,116) + .0101005_r8*rxt(k,507)*y(k,176)
         mat(k,52) = .0237005_r8*rxt(k,511)*y(k,116) + .0174005_r8*rxt(k,510)*y(k,176)
         mat(k,751) = .0653005_r8*rxt(k,487)*y(k,122)
         mat(k,778) = .0348005_r8*rxt(k,503)*y(k,122)
         mat(k,1723) = .0676005_r8*rxt(k,486)*y(k,163) + .1579005_r8*rxt(k,491) &
                      *y(k,165) + .0073005_r8*rxt(k,494)*y(k,180) &
                      + .0521005_r8*rxt(k,498)*y(k,181) + .0772005_r8*rxt(k,502) &
                      *y(k,187) + .0966005_r8*rxt(k,508)*y(k,197) &
                      + .0025005_r8*rxt(k,511)*y(k,200)
         mat(k,1813) = .0653005_r8*rxt(k,487)*y(k,4) + .0348005_r8*rxt(k,503)*y(k,105)
         mat(k,9) = .0459005_r8*rxt(k,505)*y(k,190)
         mat(k,15) = .0676005_r8*rxt(k,486)*y(k,116) + .0653005_r8*rxt(k,485)*y(k,176)
         mat(k,21) = .1579005_r8*rxt(k,491)*y(k,116) + .0843005_r8*rxt(k,490)*y(k,176)
         mat(k,1591) = .0653005_r8*rxt(k,485)*y(k,163) + .0843005_r8*rxt(k,490) &
                      *y(k,165) + .0003005_r8*rxt(k,493)*y(k,180) &
                      + .0348005_r8*rxt(k,497)*y(k,181) + .0348005_r8*rxt(k,501) &
                      *y(k,187) + .0763005_r8*rxt(k,507)*y(k,197) + .086_r8*rxt(k,510) &
                      *y(k,200)
         mat(k,27) = .0073005_r8*rxt(k,494)*y(k,116) + .0003005_r8*rxt(k,493)*y(k,176)
         mat(k,33) = .0521005_r8*rxt(k,498)*y(k,116) + .0348005_r8*rxt(k,497)*y(k,176)
         mat(k,41) = .0772005_r8*rxt(k,502)*y(k,116) + .0348005_r8*rxt(k,501)*y(k,176)
         mat(k,1360) = .0459005_r8*rxt(k,505)*y(k,145)
         mat(k,47) = .0966005_r8*rxt(k,508)*y(k,116) + .0763005_r8*rxt(k,507)*y(k,176)
         mat(k,53) = .0025005_r8*rxt(k,511)*y(k,116) + .086_r8*rxt(k,510)*y(k,176)
         mat(k,752) = .1749305_r8*rxt(k,484)*y(k,118) + .1284005_r8*rxt(k,487) &
                      *y(k,122)
         mat(k,702) = .0590245_r8*rxt(k,492)*y(k,118) + .0033005_r8*rxt(k,495) &
                      *y(k,122)
         mat(k,779) = .1749305_r8*rxt(k,500)*y(k,118) + .0554005_r8*rxt(k,503) &
                      *y(k,122)
         mat(k,1724) = .079_r8*rxt(k,486)*y(k,163) + .0059005_r8*rxt(k,491)*y(k,165) &
                      + .0057005_r8*rxt(k,494)*y(k,180) + .0143005_r8*rxt(k,498) &
                      *y(k,181) + .0332005_r8*rxt(k,502)*y(k,187) &
                      + .0073005_r8*rxt(k,508)*y(k,197) + .011_r8*rxt(k,511)*y(k,200)
         mat(k,1969) = .1749305_r8*rxt(k,484)*y(k,4) + .0590245_r8*rxt(k,492)*y(k,93) &
                      + .1749305_r8*rxt(k,500)*y(k,105)
         mat(k,1814) = .1284005_r8*rxt(k,487)*y(k,4) + .0033005_r8*rxt(k,495)*y(k,93) &
                      + .0554005_r8*rxt(k,503)*y(k,105)
         mat(k,10) = .0085005_r8*rxt(k,505)*y(k,190)
         mat(k,16) = .079_r8*rxt(k,486)*y(k,116) + .1284005_r8*rxt(k,485)*y(k,176)
         mat(k,22) = .0059005_r8*rxt(k,491)*y(k,116) + .0443005_r8*rxt(k,490)*y(k,176)
         mat(k,1592) = .1284005_r8*rxt(k,485)*y(k,163) + .0443005_r8*rxt(k,490) &
                      *y(k,165) + .0271005_r8*rxt(k,493)*y(k,180) &
                      + .0076005_r8*rxt(k,497)*y(k,181) + .0554005_r8*rxt(k,501) &
                      *y(k,187) + .2157005_r8*rxt(k,507)*y(k,197) &
                      + .0512005_r8*rxt(k,510)*y(k,200)
         mat(k,28) = .0057005_r8*rxt(k,494)*y(k,116) + .0271005_r8*rxt(k,493)*y(k,176)
         mat(k,34) = .0143005_r8*rxt(k,498)*y(k,116) + .0076005_r8*rxt(k,497)*y(k,176)
         mat(k,42) = .0332005_r8*rxt(k,502)*y(k,116) + .0554005_r8*rxt(k,501)*y(k,176)
         mat(k,1361) = .0085005_r8*rxt(k,505)*y(k,145)
         mat(k,48) = .0073005_r8*rxt(k,508)*y(k,116) + .2157005_r8*rxt(k,507)*y(k,176)
         mat(k,54) = .011_r8*rxt(k,511)*y(k,116) + .0512005_r8*rxt(k,510)*y(k,176)
         mat(k,753) = .5901905_r8*rxt(k,484)*y(k,118) + .114_r8*rxt(k,487)*y(k,122)
         mat(k,703) = .0250245_r8*rxt(k,492)*y(k,118)
         mat(k,780) = .5901905_r8*rxt(k,500)*y(k,118) + .1278005_r8*rxt(k,503) &
                      *y(k,122)
         mat(k,1725) = .1254005_r8*rxt(k,486)*y(k,163) + .0536005_r8*rxt(k,491) &
                      *y(k,165) + .0623005_r8*rxt(k,494)*y(k,180) &
                      + .0166005_r8*rxt(k,498)*y(k,181) + .130_r8*rxt(k,502)*y(k,187) &
                      + .238_r8*rxt(k,508)*y(k,197) + .1185005_r8*rxt(k,511)*y(k,200)
         mat(k,1970) = .5901905_r8*rxt(k,484)*y(k,4) + .0250245_r8*rxt(k,492)*y(k,93) &
                      + .5901905_r8*rxt(k,500)*y(k,105)
         mat(k,1815) = .114_r8*rxt(k,487)*y(k,4) + .1278005_r8*rxt(k,503)*y(k,105)
         mat(k,11) = .0128005_r8*rxt(k,505)*y(k,190)
         mat(k,17) = .1254005_r8*rxt(k,486)*y(k,116) + .114_r8*rxt(k,485)*y(k,176)
         mat(k,23) = .0536005_r8*rxt(k,491)*y(k,116) + .1621005_r8*rxt(k,490)*y(k,176)
         mat(k,1593) = .114_r8*rxt(k,485)*y(k,163) + .1621005_r8*rxt(k,490)*y(k,165) &
                      + .0474005_r8*rxt(k,493)*y(k,180) + .0113005_r8*rxt(k,497) &
                      *y(k,181) + .1278005_r8*rxt(k,501)*y(k,187) &
                      + .0738005_r8*rxt(k,507)*y(k,197) + .1598005_r8*rxt(k,510) &
                      *y(k,200)
         mat(k,29) = .0623005_r8*rxt(k,494)*y(k,116) + .0474005_r8*rxt(k,493)*y(k,176)
         mat(k,35) = .0166005_r8*rxt(k,498)*y(k,116) + .0113005_r8*rxt(k,497)*y(k,176)
         mat(k,43) = .130_r8*rxt(k,502)*y(k,116) + .1278005_r8*rxt(k,501)*y(k,176)
         mat(k,1362) = .0128005_r8*rxt(k,505)*y(k,145)
         mat(k,49) = .238_r8*rxt(k,508)*y(k,116) + .0738005_r8*rxt(k,507)*y(k,176)
         mat(k,55) = .1185005_r8*rxt(k,511)*y(k,116) + .1598005_r8*rxt(k,510)*y(k,176)
         mat(k,12) = -(rxt(k,505)*y(k,190))
         mat(k,1363) = -rxt(k,505)*y(k,145)
         mat(k,114) = .100_r8*rxt(k,419)*y(k,190)
         mat(k,132) = .230_r8*rxt(k,421)*y(k,190)
         mat(k,1384) = .100_r8*rxt(k,419)*y(k,153) + .230_r8*rxt(k,421)*y(k,155)
         mat(k,482) = -(rxt(k,443)*y(k,190))
         mat(k,1436) = -rxt(k,443)*y(k,147)
         mat(k,1631) = rxt(k,441)*y(k,194)
         mat(k,954) = rxt(k,441)*y(k,176)
         mat(k,506) = -(rxt(k,444)*y(k,190))
         mat(k,1439) = -rxt(k,444)*y(k,148)
         mat(k,1753) = .200_r8*rxt(k,437)*y(k,188) + .200_r8*rxt(k,447)*y(k,195)
         mat(k,1290) = .500_r8*rxt(k,435)*y(k,188)
         mat(k,973) = .200_r8*rxt(k,437)*y(k,116) + .500_r8*rxt(k,435)*y(k,171)
         mat(k,932) = .200_r8*rxt(k,447)*y(k,116)
         mat(k,366) = -(rxt(k,448)*y(k,190))
         mat(k,1421) = -rxt(k,448)*y(k,149)
         mat(k,1623) = rxt(k,446)*y(k,195)
         mat(k,931) = rxt(k,446)*y(k,176)
         mat(k,855) = -(rxt(k,449)*y(k,118) + rxt(k,450)*y(k,190))
         mat(k,1985) = -rxt(k,449)*y(k,150)
         mat(k,1470) = -rxt(k,450)*y(k,150)
         mat(k,762) = .330_r8*rxt(k,430)*y(k,122)
         mat(k,789) = .330_r8*rxt(k,433)*y(k,122)
         mat(k,1771) = .800_r8*rxt(k,437)*y(k,188) + .800_r8*rxt(k,447)*y(k,195)
         mat(k,1985) = mat(k,1985) + rxt(k,438)*y(k,188)
         mat(k,1833) = .330_r8*rxt(k,430)*y(k,4) + .330_r8*rxt(k,433)*y(k,105)
         mat(k,507) = rxt(k,444)*y(k,190)
         mat(k,1297) = .500_r8*rxt(k,435)*y(k,188) + rxt(k,445)*y(k,195)
         mat(k,975) = .800_r8*rxt(k,437)*y(k,116) + rxt(k,438)*y(k,118) &
                      + .500_r8*rxt(k,435)*y(k,171)
         mat(k,1470) = mat(k,1470) + rxt(k,444)*y(k,148)
         mat(k,935) = .800_r8*rxt(k,447)*y(k,116) + rxt(k,445)*y(k,171)
         mat(k,885) = -(rxt(k,451)*y(k,190))
         mat(k,1473) = -rxt(k,451)*y(k,151)
         mat(k,763) = .300_r8*rxt(k,430)*y(k,122)
         mat(k,790) = .300_r8*rxt(k,433)*y(k,122)
         mat(k,1774) = .900_r8*rxt(k,442)*y(k,194)
         mat(k,1835) = .300_r8*rxt(k,430)*y(k,4) + .300_r8*rxt(k,433)*y(k,105)
         mat(k,1300) = rxt(k,440)*y(k,194)
         mat(k,958) = .900_r8*rxt(k,442)*y(k,116) + rxt(k,440)*y(k,171)
         mat(k,493) = -(rxt(k,418)*y(k,190))
         mat(k,1437) = -rxt(k,418)*y(k,152)
         mat(k,1632) = rxt(k,416)*y(k,196)
         mat(k,582) = rxt(k,416)*y(k,176)
         mat(k,112) = -(rxt(k,419)*y(k,190))
         mat(k,1382) = -rxt(k,419)*y(k,153)
         mat(k,128) = -(rxt(k,385)*y(k,190))
         mat(k,1385) = -rxt(k,385)*y(k,154)
         mat(k,1602) = rxt(k,382)*y(k,198)
         mat(k,1005) = rxt(k,382)*y(k,176)
         mat(k,133) = -(rxt(k,421)*y(k,190))
         mat(k,1386) = -rxt(k,421)*y(k,155)
         mat(k,554) = -(rxt(k,424)*y(k,190))
         mat(k,1444) = -rxt(k,424)*y(k,156)
         mat(k,1637) = rxt(k,422)*y(k,199)
         mat(k,607) = rxt(k,422)*y(k,176)
         mat(k,141) = -(rxt(k,427)*y(k,190))
         mat(k,1387) = -rxt(k,427)*y(k,157)
         mat(k,134) = .150_r8*rxt(k,421)*y(k,190)
         mat(k,1387) = mat(k,1387) + .150_r8*rxt(k,421)*y(k,155)
         mat(k,325) = -(rxt(k,428)*y(k,190))
         mat(k,1415) = -rxt(k,428)*y(k,158)
         mat(k,1617) = rxt(k,425)*y(k,201)
         mat(k,382) = rxt(k,425)*y(k,176)
         mat(k,408) = -(rxt(k,386)*y(k,176) + rxt(k,387)*y(k,116) + rxt(k,415) &
                      *y(k,117))
         mat(k,1627) = -rxt(k,386)*y(k,161)
         mat(k,1748) = -rxt(k,387)*y(k,161)
         mat(k,1877) = -rxt(k,415)*y(k,161)
         mat(k,158) = rxt(k,392)*y(k,190)
         mat(k,1427) = rxt(k,392)*y(k,20)
         mat(k,808) = -(rxt(k,347)*y(k,176) + (rxt(k,348) + rxt(k,349)) * y(k,116))
         mat(k,1654) = -rxt(k,347)*y(k,162)
         mat(k,1769) = -(rxt(k,348) + rxt(k,349)) * y(k,162)
         mat(k,524) = rxt(k,350)*y(k,190)
         mat(k,155) = rxt(k,351)*y(k,190)
         mat(k,1466) = rxt(k,350)*y(k,2) + rxt(k,351)*y(k,13)
         mat(k,18) = -(rxt(k,485)*y(k,176) + rxt(k,486)*y(k,116))
         mat(k,1594) = -rxt(k,485)*y(k,163)
         mat(k,1726) = -rxt(k,486)*y(k,163)
         mat(k,754) = rxt(k,488)*y(k,190)
         mat(k,1364) = rxt(k,488)*y(k,4)
         mat(k,375) = -(rxt(k,389)*y(k,176) + rxt(k,390)*y(k,116))
         mat(k,1624) = -rxt(k,389)*y(k,164)
         mat(k,1745) = -rxt(k,390)*y(k,164)
         mat(k,91) = .350_r8*rxt(k,388)*y(k,190)
         mat(k,315) = rxt(k,391)*y(k,190)
         mat(k,1422) = .350_r8*rxt(k,388)*y(k,5) + rxt(k,391)*y(k,6)
         mat(k,24) = -(rxt(k,490)*y(k,176) + rxt(k,491)*y(k,116))
         mat(k,1595) = -rxt(k,490)*y(k,165)
         mat(k,1727) = -rxt(k,491)*y(k,165)
         mat(k,87) = rxt(k,489)*y(k,190)
         mat(k,1365) = rxt(k,489)*y(k,5)
         mat(k,333) = -(rxt(k,393)*y(k,176) + rxt(k,395)*y(k,116))
         mat(k,1618) = -rxt(k,393)*y(k,166)
         mat(k,1740) = -rxt(k,395)*y(k,166)
         mat(k,240) = rxt(k,394)*y(k,190)
         mat(k,115) = .070_r8*rxt(k,419)*y(k,190)
         mat(k,135) = .060_r8*rxt(k,421)*y(k,190)
         mat(k,1416) = rxt(k,394)*y(k,21) + .070_r8*rxt(k,419)*y(k,153) &
                      + .060_r8*rxt(k,421)*y(k,155)
         mat(k,686) = -(4._r8*rxt(k,270)*y(k,167) + rxt(k,271)*y(k,171) + rxt(k,272) &
                      *y(k,176) + rxt(k,273)*y(k,116))
         mat(k,1293) = -rxt(k,271)*y(k,167)
         mat(k,1649) = -rxt(k,272)*y(k,167)
         mat(k,1765) = -rxt(k,273)*y(k,167)
         mat(k,214) = .500_r8*rxt(k,275)*y(k,190)
         mat(k,196) = rxt(k,276)*y(k,54) + rxt(k,277)*y(k,190)
         mat(k,1943) = rxt(k,276)*y(k,26)
         mat(k,1457) = .500_r8*rxt(k,275)*y(k,25) + rxt(k,277)*y(k,26)
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
         mat(k,661) = -(rxt(k,299)*y(k,171) + rxt(k,300)*y(k,176) + rxt(k,301) &
                      *y(k,116))
         mat(k,1292) = -rxt(k,299)*y(k,168)
         mat(k,1646) = -rxt(k,300)*y(k,168)
         mat(k,1763) = -rxt(k,301)*y(k,168)
         mat(k,308) = rxt(k,302)*y(k,190)
         mat(k,67) = rxt(k,303)*y(k,190)
         mat(k,1454) = rxt(k,302)*y(k,28) + rxt(k,303)*y(k,29)
         mat(k,514) = -(rxt(k,396)*y(k,176) + rxt(k,397)*y(k,116))
         mat(k,1634) = -rxt(k,396)*y(k,169)
         mat(k,1754) = -rxt(k,397)*y(k,169)
         mat(k,172) = rxt(k,398)*y(k,190)
         mat(k,1754) = mat(k,1754) + rxt(k,387)*y(k,161)
         mat(k,1822) = rxt(k,413)*y(k,129)
         mat(k,363) = rxt(k,413)*y(k,122)
         mat(k,409) = rxt(k,387)*y(k,116) + .400_r8*rxt(k,386)*y(k,176)
         mat(k,1634) = mat(k,1634) + .400_r8*rxt(k,386)*y(k,161)
         mat(k,1440) = rxt(k,398)*y(k,30)
         mat(k,1233) = -(4._r8*rxt(k,281)*y(k,170) + rxt(k,282)*y(k,171) + rxt(k,283) &
                      *y(k,176) + rxt(k,284)*y(k,116) + rxt(k,295)*y(k,117) + rxt(k,322) &
                      *y(k,182) + rxt(k,355)*y(k,178) + rxt(k,360)*y(k,179) + rxt(k,369) &
                      *y(k,96) + rxt(k,380)*y(k,198))
         mat(k,1317) = -rxt(k,282)*y(k,170)
         mat(k,1676) = -rxt(k,283)*y(k,170)
         mat(k,1792) = -rxt(k,284)*y(k,170)
         mat(k,1893) = -rxt(k,295)*y(k,170)
         mat(k,1159) = -rxt(k,322)*y(k,170)
         mat(k,1185) = -rxt(k,355)*y(k,170)
         mat(k,1117) = -rxt(k,360)*y(k,170)
         mat(k,1048) = -rxt(k,369)*y(k,170)
         mat(k,1013) = -rxt(k,380)*y(k,170)
         mat(k,769) = .060_r8*rxt(k,430)*y(k,122)
         mat(k,896) = rxt(k,278)*y(k,118) + rxt(k,279)*y(k,190)
         mat(k,1073) = rxt(k,304)*y(k,118) + rxt(k,305)*y(k,190)
         mat(k,402) = .500_r8*rxt(k,286)*y(k,190)
         mat(k,714) = .080_r8*rxt(k,375)*y(k,122)
         mat(k,1064) = .100_r8*rxt(k,328)*y(k,122)
         mat(k,796) = .060_r8*rxt(k,433)*y(k,122)
         mat(k,1137) = .280_r8*rxt(k,342)*y(k,122)
         mat(k,1792) = mat(k,1792) + .530_r8*rxt(k,326)*y(k,182) + rxt(k,335)*y(k,184) &
                      + rxt(k,338)*y(k,186) + rxt(k,313)*y(k,193)
         mat(k,2007) = rxt(k,278)*y(k,43) + rxt(k,304)*y(k,47) + .530_r8*rxt(k,325) &
                      *y(k,182) + rxt(k,336)*y(k,184)
         mat(k,1852) = .060_r8*rxt(k,430)*y(k,4) + .080_r8*rxt(k,375)*y(k,93) &
                      + .100_r8*rxt(k,328)*y(k,100) + .060_r8*rxt(k,433)*y(k,105) &
                      + .280_r8*rxt(k,342)*y(k,106)
         mat(k,888) = .650_r8*rxt(k,451)*y(k,190)
         mat(k,1233) = mat(k,1233) + .530_r8*rxt(k,322)*y(k,182)
         mat(k,1317) = mat(k,1317) + .260_r8*rxt(k,323)*y(k,182) + rxt(k,332)*y(k,184) &
                      + .300_r8*rxt(k,311)*y(k,193)
         mat(k,1676) = mat(k,1676) + .450_r8*rxt(k,333)*y(k,184) + .200_r8*rxt(k,337) &
                      *y(k,186) + .150_r8*rxt(k,312)*y(k,193)
         mat(k,1159) = mat(k,1159) + .530_r8*rxt(k,326)*y(k,116) + .530_r8*rxt(k,325) &
                      *y(k,118) + .530_r8*rxt(k,322)*y(k,170) + .260_r8*rxt(k,323) &
                      *y(k,171)
         mat(k,1203) = rxt(k,335)*y(k,116) + rxt(k,336)*y(k,118) + rxt(k,332)*y(k,171) &
                      + .450_r8*rxt(k,333)*y(k,176) + 4.000_r8*rxt(k,334)*y(k,184)
         mat(k,545) = rxt(k,338)*y(k,116) + .200_r8*rxt(k,337)*y(k,176)
         mat(k,1493) = rxt(k,279)*y(k,43) + rxt(k,305)*y(k,47) + .500_r8*rxt(k,286) &
                      *y(k,49) + .650_r8*rxt(k,451)*y(k,151)
         mat(k,997) = rxt(k,313)*y(k,116) + .300_r8*rxt(k,311)*y(k,171) &
                      + .150_r8*rxt(k,312)*y(k,176)
         mat(k,1320) = -(rxt(k,171)*y(k,57) + (4._r8*rxt(k,248) + 4._r8*rxt(k,249) &
                      ) * y(k,171) + rxt(k,250)*y(k,176) + rxt(k,251)*y(k,116) &
                      + rxt(k,271)*y(k,167) + rxt(k,282)*y(k,170) + rxt(k,299) &
                      *y(k,168) + rxt(k,311)*y(k,193) + rxt(k,323)*y(k,182) + rxt(k,332) &
                      *y(k,184) + rxt(k,356)*y(k,178) + rxt(k,361)*y(k,179) + rxt(k,370) &
                      *y(k,96) + rxt(k,381)*y(k,198) + rxt(k,435)*y(k,188) + rxt(k,440) &
                      *y(k,194) + rxt(k,445)*y(k,195))
         mat(k,1706) = -rxt(k,171)*y(k,171)
         mat(k,1680) = -rxt(k,250)*y(k,171)
         mat(k,1796) = -rxt(k,251)*y(k,171)
         mat(k,688) = -rxt(k,271)*y(k,171)
         mat(k,1236) = -rxt(k,282)*y(k,171)
         mat(k,664) = -rxt(k,299)*y(k,171)
         mat(k,998) = -rxt(k,311)*y(k,171)
         mat(k,1162) = -rxt(k,323)*y(k,171)
         mat(k,1206) = -rxt(k,332)*y(k,171)
         mat(k,1188) = -rxt(k,356)*y(k,171)
         mat(k,1120) = -rxt(k,361)*y(k,171)
         mat(k,1051) = -rxt(k,370)*y(k,171)
         mat(k,1015) = -rxt(k,381)*y(k,171)
         mat(k,984) = -rxt(k,435)*y(k,171)
         mat(k,965) = -rxt(k,440)*y(k,171)
         mat(k,945) = -rxt(k,445)*y(k,171)
         mat(k,846) = .280_r8*rxt(k,298)*y(k,122)
         mat(k,431) = rxt(k,285)*y(k,190)
         mat(k,297) = .700_r8*rxt(k,253)*y(k,190)
         mat(k,715) = .050_r8*rxt(k,375)*y(k,122)
         mat(k,1051) = mat(k,1051) + rxt(k,369)*y(k,170)
         mat(k,1796) = mat(k,1796) + rxt(k,284)*y(k,170) + .830_r8*rxt(k,401)*y(k,172) &
                      + .170_r8*rxt(k,407)*y(k,185)
         mat(k,1856) = .280_r8*rxt(k,298)*y(k,27) + .050_r8*rxt(k,375)*y(k,93)
         mat(k,1236) = mat(k,1236) + rxt(k,369)*y(k,96) + rxt(k,284)*y(k,116) &
                      + 4.000_r8*rxt(k,281)*y(k,170) + .900_r8*rxt(k,282)*y(k,171) &
                      + .450_r8*rxt(k,283)*y(k,176) + rxt(k,355)*y(k,178) + rxt(k,360) &
                      *y(k,179) + rxt(k,322)*y(k,182) + rxt(k,331)*y(k,184) &
                      + rxt(k,380)*y(k,198)
         mat(k,1320) = mat(k,1320) + .900_r8*rxt(k,282)*y(k,170)
         mat(k,623) = .830_r8*rxt(k,401)*y(k,116) + .330_r8*rxt(k,400)*y(k,176)
         mat(k,1680) = mat(k,1680) + .450_r8*rxt(k,283)*y(k,170) + .330_r8*rxt(k,400) &
                      *y(k,172) + .070_r8*rxt(k,406)*y(k,185)
         mat(k,1188) = mat(k,1188) + rxt(k,355)*y(k,170)
         mat(k,1120) = mat(k,1120) + rxt(k,360)*y(k,170)
         mat(k,1162) = mat(k,1162) + rxt(k,322)*y(k,170)
         mat(k,1206) = mat(k,1206) + rxt(k,331)*y(k,170)
         mat(k,744) = .170_r8*rxt(k,407)*y(k,116) + .070_r8*rxt(k,406)*y(k,176)
         mat(k,1497) = rxt(k,285)*y(k,48) + .700_r8*rxt(k,253)*y(k,51)
         mat(k,1015) = mat(k,1015) + rxt(k,380)*y(k,170)
         mat(k,620) = -(rxt(k,400)*y(k,176) + rxt(k,401)*y(k,116) + rxt(k,402) &
                      *y(k,117))
         mat(k,1643) = -rxt(k,400)*y(k,172)
         mat(k,1761) = -rxt(k,401)*y(k,172)
         mat(k,1882) = -rxt(k,402)*y(k,172)
         mat(k,450) = -((rxt(k,319) + rxt(k,320)) * y(k,116))
         mat(k,1750) = -(rxt(k,319) + rxt(k,320)) * y(k,173)
         mat(k,258) = rxt(k,318)*y(k,190)
         mat(k,1432) = rxt(k,318)*y(k,14)
         mat(k,1735) = .750_r8*rxt(k,288)*y(k,175)
         mat(k,566) = .750_r8*rxt(k,288)*y(k,116)
         mat(k,567) = -(rxt(k,287)*y(k,176) + rxt(k,288)*y(k,116))
         mat(k,1638) = -rxt(k,287)*y(k,175)
         mat(k,1757) = -rxt(k,288)*y(k,175)
         mat(k,443) = rxt(k,294)*y(k,190)
         mat(k,1445) = rxt(k,294)*y(k,23)
         mat(k,1686) = -((rxt(k,124) + rxt(k,125) + rxt(k,126)) * y(k,71) + rxt(k,128) &
                      *y(k,121) + rxt(k,129)*y(k,122) + rxt(k,133)*y(k,190) &
                      + 4._r8*rxt(k,138)*y(k,176) + rxt(k,148)*y(k,118) + rxt(k,153) &
                      *y(k,116) + rxt(k,158)*y(k,117) + (rxt(k,168) + rxt(k,169) &
                      ) * y(k,54) + rxt(k,175)*y(k,57) + rxt(k,201)*y(k,15) + rxt(k,207) &
                      *y(k,17) + rxt(k,244)*y(k,40) + rxt(k,250)*y(k,171) + rxt(k,258) &
                      *y(k,177) + rxt(k,272)*y(k,167) + rxt(k,283)*y(k,170) + rxt(k,287) &
                      *y(k,175) + rxt(k,300)*y(k,168) + rxt(k,308)*y(k,192) + rxt(k,312) &
                      *y(k,193) + rxt(k,324)*y(k,182) + rxt(k,333)*y(k,184) + rxt(k,337) &
                      *y(k,186) + rxt(k,347)*y(k,162) + rxt(k,357)*y(k,178) + rxt(k,362) &
                      *y(k,179) + rxt(k,371)*y(k,96) + rxt(k,382)*y(k,198) + rxt(k,386) &
                      *y(k,161) + rxt(k,389)*y(k,164) + rxt(k,393)*y(k,166) + rxt(k,396) &
                      *y(k,169) + rxt(k,400)*y(k,172) + rxt(k,403)*y(k,183) + rxt(k,406) &
                      *y(k,185) + rxt(k,409)*y(k,191) + rxt(k,416)*y(k,196) + rxt(k,422) &
                      *y(k,199) + rxt(k,425)*y(k,201) + rxt(k,436)*y(k,188) + rxt(k,441) &
                      *y(k,194) + rxt(k,446)*y(k,195))
         mat(k,1278) = -(rxt(k,124) + rxt(k,125) + rxt(k,126)) * y(k,176)
         mat(k,1580) = -rxt(k,128)*y(k,176)
         mat(k,1862) = -rxt(k,129)*y(k,176)
         mat(k,1503) = -rxt(k,133)*y(k,176)
         mat(k,2017) = -rxt(k,148)*y(k,176)
         mat(k,1802) = -rxt(k,153)*y(k,176)
         mat(k,1903) = -rxt(k,158)*y(k,176)
         mat(k,1960) = -(rxt(k,168) + rxt(k,169)) * y(k,176)
         mat(k,1712) = -rxt(k,175)*y(k,176)
         mat(k,1254) = -rxt(k,201)*y(k,176)
         mat(k,1527) = -rxt(k,207)*y(k,176)
         mat(k,1550) = -rxt(k,244)*y(k,176)
         mat(k,1325) = -rxt(k,250)*y(k,176)
         mat(k,342) = -rxt(k,258)*y(k,176)
         mat(k,691) = -rxt(k,272)*y(k,176)
         mat(k,1239) = -rxt(k,283)*y(k,176)
         mat(k,571) = -rxt(k,287)*y(k,176)
         mat(k,667) = -rxt(k,300)*y(k,176)
         mat(k,634) = -rxt(k,308)*y(k,176)
         mat(k,1001) = -rxt(k,312)*y(k,176)
         mat(k,1165) = -rxt(k,324)*y(k,176)
         mat(k,1209) = -rxt(k,333)*y(k,176)
         mat(k,547) = -rxt(k,337)*y(k,176)
         mat(k,815) = -rxt(k,347)*y(k,176)
         mat(k,1191) = -rxt(k,357)*y(k,176)
         mat(k,1123) = -rxt(k,362)*y(k,176)
         mat(k,1054) = -rxt(k,371)*y(k,176)
         mat(k,1018) = -rxt(k,382)*y(k,176)
         mat(k,411) = -rxt(k,386)*y(k,176)
         mat(k,379) = -rxt(k,389)*y(k,176)
         mat(k,336) = -rxt(k,393)*y(k,176)
         mat(k,516) = -rxt(k,396)*y(k,176)
         mat(k,625) = -rxt(k,400)*y(k,176)
         mat(k,577) = -rxt(k,403)*y(k,176)
         mat(k,746) = -rxt(k,406)*y(k,176)
         mat(k,349) = -rxt(k,409)*y(k,176)
         mat(k,592) = -rxt(k,416)*y(k,176)
         mat(k,617) = -rxt(k,422)*y(k,176)
         mat(k,387) = -rxt(k,425)*y(k,176)
         mat(k,987) = -rxt(k,436)*y(k,176)
         mat(k,968) = -rxt(k,441)*y(k,176)
         mat(k,948) = -rxt(k,446)*y(k,176)
         mat(k,772) = .570_r8*rxt(k,430)*y(k,122)
         mat(k,93) = .650_r8*rxt(k,388)*y(k,190)
         mat(k,1254) = mat(k,1254) + rxt(k,200)*y(k,40)
         mat(k,1527) = mat(k,1527) + rxt(k,212)*y(k,190)
         mat(k,206) = .350_r8*rxt(k,267)*y(k,190)
         mat(k,448) = .130_r8*rxt(k,269)*y(k,122)
         mat(k,169) = rxt(k,274)*y(k,190)
         mat(k,849) = .280_r8*rxt(k,298)*y(k,122)
         mat(k,1550) = mat(k,1550) + rxt(k,200)*y(k,15) + rxt(k,164)*y(k,54) &
                      + rxt(k,245)*y(k,118) + rxt(k,246)*y(k,121)
         mat(k,59) = rxt(k,280)*y(k,190)
         mat(k,659) = rxt(k,252)*y(k,190)
         mat(k,1960) = mat(k,1960) + rxt(k,164)*y(k,40) + rxt(k,167)*y(k,74)
         mat(k,1712) = mat(k,1712) + rxt(k,171)*y(k,171) + rxt(k,182)*y(k,190)
         mat(k,906) = rxt(k,255)*y(k,190)
         mat(k,123) = .730_r8*rxt(k,399)*y(k,190)
         mat(k,187) = .500_r8*rxt(k,466)*y(k,190)
         mat(k,868) = rxt(k,291)*y(k,190)
         mat(k,739) = rxt(k,292)*y(k,190)
         mat(k,469) = rxt(k,167)*y(k,54) + rxt(k,123)*y(k,121) + rxt(k,132)*y(k,190)
         mat(k,110) = rxt(k,256)*y(k,190)
         mat(k,673) = rxt(k,257)*y(k,190)
         mat(k,921) = rxt(k,321)*y(k,190)
         mat(k,930) = rxt(k,306)*y(k,190)
         mat(k,718) = .370_r8*rxt(k,375)*y(k,122)
         mat(k,480) = .300_r8*rxt(k,366)*y(k,190)
         mat(k,441) = rxt(k,367)*y(k,190)
         mat(k,1054) = mat(k,1054) + rxt(k,372)*y(k,116) + rxt(k,373)*y(k,118) &
                      + rxt(k,369)*y(k,170) + 1.200_r8*rxt(k,370)*y(k,171)
         mat(k,275) = rxt(k,374)*y(k,190)
         mat(k,1068) = .140_r8*rxt(k,328)*y(k,122)
         mat(k,222) = .200_r8*rxt(k,330)*y(k,190)
         mat(k,463) = .500_r8*rxt(k,341)*y(k,190)
         mat(k,799) = .570_r8*rxt(k,433)*y(k,122)
         mat(k,1143) = .280_r8*rxt(k,342)*y(k,122)
         mat(k,288) = rxt(k,378)*y(k,190)
         mat(k,881) = rxt(k,379)*y(k,190)
         mat(k,1802) = mat(k,1802) + rxt(k,372)*y(k,96) + rxt(k,348)*y(k,162) &
                      + rxt(k,390)*y(k,164) + rxt(k,395)*y(k,166) + rxt(k,273) &
                      *y(k,167) + rxt(k,301)*y(k,168) + rxt(k,251)*y(k,171) &
                      + .170_r8*rxt(k,401)*y(k,172) + rxt(k,319)*y(k,173) &
                      + .250_r8*rxt(k,288)*y(k,175) + rxt(k,260)*y(k,177) &
                      + .920_r8*rxt(k,358)*y(k,178) + .920_r8*rxt(k,364)*y(k,179) &
                      + .470_r8*rxt(k,326)*y(k,182) + .400_r8*rxt(k,404)*y(k,183) &
                      + .830_r8*rxt(k,407)*y(k,185) + rxt(k,410)*y(k,191) + rxt(k,309) &
                      *y(k,192) + .900_r8*rxt(k,442)*y(k,194) + .800_r8*rxt(k,447) &
                      *y(k,195) + rxt(k,417)*y(k,196) + rxt(k,383)*y(k,198) &
                      + rxt(k,423)*y(k,199) + rxt(k,426)*y(k,201)
         mat(k,2017) = mat(k,2017) + rxt(k,245)*y(k,40) + rxt(k,373)*y(k,96) &
                      + rxt(k,359)*y(k,178) + rxt(k,365)*y(k,179) + .470_r8*rxt(k,325) &
                      *y(k,182) + rxt(k,151)*y(k,190) + rxt(k,384)*y(k,198)
         mat(k,1580) = mat(k,1580) + rxt(k,246)*y(k,40) + rxt(k,123)*y(k,74)
         mat(k,1862) = mat(k,1862) + .570_r8*rxt(k,430)*y(k,4) + .130_r8*rxt(k,269) &
                      *y(k,23) + .280_r8*rxt(k,298)*y(k,27) + .370_r8*rxt(k,375) &
                      *y(k,93) + .140_r8*rxt(k,328)*y(k,100) + .570_r8*rxt(k,433) &
                      *y(k,105) + .280_r8*rxt(k,342)*y(k,106) + rxt(k,135)*y(k,190)
         mat(k,102) = .800_r8*rxt(k,411)*y(k,190)
         mat(k,727) = rxt(k,467)*y(k,190)
         mat(k,892) = .200_r8*rxt(k,451)*y(k,190)
         mat(k,118) = .280_r8*rxt(k,419)*y(k,190)
         mat(k,140) = .380_r8*rxt(k,421)*y(k,190)
         mat(k,145) = .630_r8*rxt(k,427)*y(k,190)
         mat(k,815) = mat(k,815) + rxt(k,348)*y(k,116)
         mat(k,379) = mat(k,379) + rxt(k,390)*y(k,116)
         mat(k,336) = mat(k,336) + rxt(k,395)*y(k,116)
         mat(k,691) = mat(k,691) + rxt(k,273)*y(k,116) + 2.400_r8*rxt(k,270)*y(k,167) &
                      + rxt(k,271)*y(k,171)
         mat(k,667) = mat(k,667) + rxt(k,301)*y(k,116) + rxt(k,299)*y(k,171)
         mat(k,1239) = mat(k,1239) + rxt(k,369)*y(k,96) + .900_r8*rxt(k,282)*y(k,171) &
                      + rxt(k,355)*y(k,178) + rxt(k,360)*y(k,179) + .470_r8*rxt(k,322) &
                      *y(k,182) + rxt(k,380)*y(k,198)
         mat(k,1325) = mat(k,1325) + rxt(k,171)*y(k,57) + 1.200_r8*rxt(k,370)*y(k,96) &
                      + rxt(k,251)*y(k,116) + rxt(k,271)*y(k,167) + rxt(k,299) &
                      *y(k,168) + .900_r8*rxt(k,282)*y(k,170) + 4.000_r8*rxt(k,248) &
                      *y(k,171) + rxt(k,356)*y(k,178) + rxt(k,361)*y(k,179) &
                      + .730_r8*rxt(k,323)*y(k,182) + rxt(k,332)*y(k,184) &
                      + .500_r8*rxt(k,435)*y(k,188) + .300_r8*rxt(k,311)*y(k,193) &
                      + rxt(k,440)*y(k,194) + rxt(k,445)*y(k,195) + .800_r8*rxt(k,381) &
                      *y(k,198)
         mat(k,625) = mat(k,625) + .170_r8*rxt(k,401)*y(k,116) + .070_r8*rxt(k,400) &
                      *y(k,176)
         mat(k,455) = rxt(k,319)*y(k,116)
         mat(k,571) = mat(k,571) + .250_r8*rxt(k,288)*y(k,116)
         mat(k,1686) = mat(k,1686) + .070_r8*rxt(k,400)*y(k,172) + .160_r8*rxt(k,403) &
                      *y(k,183) + .330_r8*rxt(k,406)*y(k,185)
         mat(k,342) = mat(k,342) + rxt(k,260)*y(k,116)
         mat(k,1191) = mat(k,1191) + .920_r8*rxt(k,358)*y(k,116) + rxt(k,359)*y(k,118) &
                      + rxt(k,355)*y(k,170) + rxt(k,356)*y(k,171)
         mat(k,1123) = mat(k,1123) + .920_r8*rxt(k,364)*y(k,116) + rxt(k,365)*y(k,118) &
                      + rxt(k,360)*y(k,170) + rxt(k,361)*y(k,171)
         mat(k,1165) = mat(k,1165) + .470_r8*rxt(k,326)*y(k,116) + .470_r8*rxt(k,325) &
                      *y(k,118) + .470_r8*rxt(k,322)*y(k,170) + .730_r8*rxt(k,323) &
                      *y(k,171)
         mat(k,577) = mat(k,577) + .400_r8*rxt(k,404)*y(k,116) + .160_r8*rxt(k,403) &
                      *y(k,176)
         mat(k,1209) = mat(k,1209) + rxt(k,332)*y(k,171)
         mat(k,746) = mat(k,746) + .830_r8*rxt(k,407)*y(k,116) + .330_r8*rxt(k,406) &
                      *y(k,176)
         mat(k,987) = mat(k,987) + .500_r8*rxt(k,435)*y(k,171)
         mat(k,1503) = mat(k,1503) + .650_r8*rxt(k,388)*y(k,5) + rxt(k,212)*y(k,17) &
                      + .350_r8*rxt(k,267)*y(k,22) + rxt(k,274)*y(k,24) + rxt(k,280) &
                      *y(k,45) + rxt(k,252)*y(k,50) + rxt(k,182)*y(k,57) + rxt(k,255) &
                      *y(k,60) + .730_r8*rxt(k,399)*y(k,64) + .500_r8*rxt(k,466) &
                      *y(k,65) + rxt(k,291)*y(k,69) + rxt(k,292)*y(k,70) + rxt(k,132) &
                      *y(k,74) + rxt(k,256)*y(k,81) + rxt(k,257)*y(k,82) + rxt(k,321) &
                      *y(k,88) + rxt(k,306)*y(k,90) + .300_r8*rxt(k,366)*y(k,94) &
                      + rxt(k,367)*y(k,95) + rxt(k,374)*y(k,97) + .200_r8*rxt(k,330) &
                      *y(k,101) + .500_r8*rxt(k,341)*y(k,104) + rxt(k,378)*y(k,110) &
                      + rxt(k,379)*y(k,111) + rxt(k,151)*y(k,118) + rxt(k,135) &
                      *y(k,122) + .800_r8*rxt(k,411)*y(k,130) + rxt(k,467)*y(k,137) &
                      + .200_r8*rxt(k,451)*y(k,151) + .280_r8*rxt(k,419)*y(k,153) &
                      + .380_r8*rxt(k,421)*y(k,155) + .630_r8*rxt(k,427)*y(k,157)
         mat(k,349) = mat(k,349) + rxt(k,410)*y(k,116)
         mat(k,634) = mat(k,634) + rxt(k,309)*y(k,116)
         mat(k,1001) = mat(k,1001) + .300_r8*rxt(k,311)*y(k,171)
         mat(k,968) = mat(k,968) + .900_r8*rxt(k,442)*y(k,116) + rxt(k,440)*y(k,171)
         mat(k,948) = mat(k,948) + .800_r8*rxt(k,447)*y(k,116) + rxt(k,445)*y(k,171)
         mat(k,592) = mat(k,592) + rxt(k,417)*y(k,116)
         mat(k,1018) = mat(k,1018) + rxt(k,383)*y(k,116) + rxt(k,384)*y(k,118) &
                      + rxt(k,380)*y(k,170) + .800_r8*rxt(k,381)*y(k,171)
         mat(k,617) = mat(k,617) + rxt(k,423)*y(k,116)
         mat(k,387) = mat(k,387) + rxt(k,426)*y(k,116)
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
         mat(k,339) = -(rxt(k,258)*y(k,176) + rxt(k,260)*y(k,116))
         mat(k,1619) = -rxt(k,258)*y(k,177)
         mat(k,1741) = -rxt(k,260)*y(k,177)
         mat(k,1536) = rxt(k,244)*y(k,176)
         mat(k,1619) = mat(k,1619) + rxt(k,244)*y(k,40)
         mat(k,1183) = -(rxt(k,355)*y(k,170) + rxt(k,356)*y(k,171) + rxt(k,357) &
                      *y(k,176) + rxt(k,358)*y(k,116) + rxt(k,359)*y(k,118))
         mat(k,1231) = -rxt(k,355)*y(k,178)
         mat(k,1315) = -rxt(k,356)*y(k,178)
         mat(k,1674) = -rxt(k,357)*y(k,178)
         mat(k,1790) = -rxt(k,358)*y(k,178)
         mat(k,2005) = -rxt(k,359)*y(k,178)
         mat(k,713) = .600_r8*rxt(k,376)*y(k,190)
         mat(k,1491) = .600_r8*rxt(k,376)*y(k,93)
         mat(k,1113) = -(rxt(k,360)*y(k,170) + rxt(k,361)*y(k,171) + rxt(k,362) &
                      *y(k,176) + rxt(k,364)*y(k,116) + rxt(k,365)*y(k,118))
         mat(k,1228) = -rxt(k,360)*y(k,179)
         mat(k,1312) = -rxt(k,361)*y(k,179)
         mat(k,1671) = -rxt(k,362)*y(k,179)
         mat(k,1787) = -rxt(k,364)*y(k,179)
         mat(k,2002) = -rxt(k,365)*y(k,179)
         mat(k,711) = .400_r8*rxt(k,376)*y(k,190)
         mat(k,1488) = .400_r8*rxt(k,376)*y(k,93)
         mat(k,30) = -(rxt(k,493)*y(k,176) + rxt(k,494)*y(k,116))
         mat(k,1596) = -rxt(k,493)*y(k,180)
         mat(k,1728) = -rxt(k,494)*y(k,180)
         mat(k,704) = rxt(k,496)*y(k,190)
         mat(k,1366) = rxt(k,496)*y(k,93)
         mat(k,36) = -(rxt(k,497)*y(k,176) + rxt(k,498)*y(k,116))
         mat(k,1597) = -rxt(k,497)*y(k,181)
         mat(k,1729) = -rxt(k,498)*y(k,181)
         mat(k,37) = rxt(k,499)*y(k,190)
         mat(k,1367) = rxt(k,499)*y(k,99)
         mat(k,1157) = -(rxt(k,322)*y(k,170) + rxt(k,323)*y(k,171) + rxt(k,324) &
                      *y(k,176) + rxt(k,325)*y(k,118) + (rxt(k,326) + rxt(k,327) &
                      ) * y(k,116))
         mat(k,1230) = -rxt(k,322)*y(k,182)
         mat(k,1314) = -rxt(k,323)*y(k,182)
         mat(k,1673) = -rxt(k,324)*y(k,182)
         mat(k,2004) = -rxt(k,325)*y(k,182)
         mat(k,1789) = -(rxt(k,326) + rxt(k,327)) * y(k,182)
         mat(k,1062) = .500_r8*rxt(k,329)*y(k,190)
         mat(k,219) = .200_r8*rxt(k,330)*y(k,190)
         mat(k,1136) = rxt(k,343)*y(k,190)
         mat(k,1490) = .500_r8*rxt(k,329)*y(k,100) + .200_r8*rxt(k,330)*y(k,101) &
                      + rxt(k,343)*y(k,106)
         mat(k,574) = -(rxt(k,403)*y(k,176) + rxt(k,404)*y(k,116) + rxt(k,405) &
                      *y(k,117))
         mat(k,1639) = -rxt(k,403)*y(k,183)
         mat(k,1758) = -rxt(k,404)*y(k,183)
         mat(k,1881) = -rxt(k,405)*y(k,183)
         mat(k,1202) = -(rxt(k,331)*y(k,170) + rxt(k,332)*y(k,171) + rxt(k,333) &
                      *y(k,176) + 4._r8*rxt(k,334)*y(k,184) + rxt(k,335)*y(k,116) &
                      + rxt(k,336)*y(k,118) + rxt(k,344)*y(k,117))
         mat(k,1232) = -rxt(k,331)*y(k,184)
         mat(k,1316) = -rxt(k,332)*y(k,184)
         mat(k,1675) = -rxt(k,333)*y(k,184)
         mat(k,1791) = -rxt(k,335)*y(k,184)
         mat(k,2006) = -rxt(k,336)*y(k,184)
         mat(k,1892) = -rxt(k,344)*y(k,184)
         mat(k,1063) = .500_r8*rxt(k,329)*y(k,190)
         mat(k,220) = .500_r8*rxt(k,330)*y(k,190)
         mat(k,1492) = .500_r8*rxt(k,329)*y(k,100) + .500_r8*rxt(k,330)*y(k,101)
         mat(k,741) = -(rxt(k,406)*y(k,176) + rxt(k,407)*y(k,116) + rxt(k,408) &
                      *y(k,117))
         mat(k,1653) = -rxt(k,406)*y(k,185)
         mat(k,1768) = -rxt(k,407)*y(k,185)
         mat(k,1886) = -rxt(k,408)*y(k,185)
         mat(k,543) = -(rxt(k,337)*y(k,176) + rxt(k,338)*y(k,116))
         mat(k,1636) = -rxt(k,337)*y(k,186)
         mat(k,1756) = -rxt(k,338)*y(k,186)
         mat(k,397) = rxt(k,339)*y(k,190)
         mat(k,224) = rxt(k,340)*y(k,190)
         mat(k,1443) = rxt(k,339)*y(k,102) + rxt(k,340)*y(k,103)
         mat(k,44) = -(rxt(k,501)*y(k,176) + rxt(k,502)*y(k,116))
         mat(k,1598) = -rxt(k,501)*y(k,187)
         mat(k,1730) = -rxt(k,502)*y(k,187)
         mat(k,781) = rxt(k,504)*y(k,190)
         mat(k,1369) = rxt(k,504)*y(k,105)
         mat(k,979) = -(rxt(k,435)*y(k,171) + rxt(k,436)*y(k,176) + rxt(k,437) &
                      *y(k,116) + rxt(k,438)*y(k,118))
         mat(k,1306) = -rxt(k,435)*y(k,188)
         mat(k,1664) = -rxt(k,436)*y(k,188)
         mat(k,1781) = -rxt(k,437)*y(k,188)
         mat(k,1995) = -rxt(k,438)*y(k,188)
         mat(k,766) = rxt(k,429)*y(k,118)
         mat(k,793) = rxt(k,432)*y(k,118)
         mat(k,1995) = mat(k,1995) + rxt(k,429)*y(k,4) + rxt(k,432)*y(k,105) &
                      + .500_r8*rxt(k,449)*y(k,150)
         mat(k,291) = rxt(k,439)*y(k,190)
         mat(k,859) = .500_r8*rxt(k,449)*y(k,118)
         mat(k,1480) = rxt(k,439)*y(k,120)
         mat(k,1344) = -(rxt(k,114)*y(k,72) + rxt(k,115)*y(k,202) + rxt(k,118) &
                      *y(k,122) + (rxt(k,196) + rxt(k,197)) * y(k,80) + (rxt(k,219) &
                      + rxt(k,220)) * y(k,76) + rxt(k,225)*y(k,62) + rxt(k,226) &
                      *y(k,63) + rxt(k,264)*y(k,81))
         mat(k,1028) = -rxt(k,114)*y(k,189)
         mat(k,2037) = -rxt(k,115)*y(k,189)
         mat(k,1857) = -rxt(k,118)*y(k,189)
         mat(k,1921) = -(rxt(k,196) + rxt(k,197)) * y(k,189)
         mat(k,697) = -(rxt(k,219) + rxt(k,220)) * y(k,189)
         mat(k,72) = -rxt(k,225)*y(k,189)
         mat(k,105) = -rxt(k,226)*y(k,189)
         mat(k,108) = -rxt(k,264)*y(k,189)
         mat(k,1499) = -(rxt(k,131)*y(k,72) + rxt(k,132)*y(k,74) + rxt(k,133)*y(k,176) &
                      + rxt(k,134)*y(k,121) + rxt(k,135)*y(k,122) + (4._r8*rxt(k,136) &
                      + 4._r8*rxt(k,137)) * y(k,190) + rxt(k,139)*y(k,85) + rxt(k,151) &
                      *y(k,118) + rxt(k,152)*y(k,107) + rxt(k,160)*y(k,117) + rxt(k,161) &
                      *y(k,84) + rxt(k,180)*y(k,58) + (rxt(k,182) + rxt(k,183) &
                      ) * y(k,57) + rxt(k,185)*y(k,80) + rxt(k,188)*y(k,87) + rxt(k,212) &
                      *y(k,17) + rxt(k,214)*y(k,76) + rxt(k,247)*y(k,40) + rxt(k,252) &
                      *y(k,50) + rxt(k,253)*y(k,51) + (rxt(k,255) + rxt(k,265) &
                      ) * y(k,60) + rxt(k,256)*y(k,81) + rxt(k,257)*y(k,82) + rxt(k,267) &
                      *y(k,22) + rxt(k,274)*y(k,24) + rxt(k,275)*y(k,25) + rxt(k,277) &
                      *y(k,26) + rxt(k,279)*y(k,43) + rxt(k,280)*y(k,45) + rxt(k,285) &
                      *y(k,48) + rxt(k,286)*y(k,49) + rxt(k,291)*y(k,69) + rxt(k,292) &
                      *y(k,70) + rxt(k,293)*y(k,127) + rxt(k,294)*y(k,23) + rxt(k,302) &
                      *y(k,28) + rxt(k,303)*y(k,29) + rxt(k,305)*y(k,47) + rxt(k,306) &
                      *y(k,90) + rxt(k,307)*y(k,119) + rxt(k,310)*y(k,132) + rxt(k,314) &
                      *y(k,133) + rxt(k,315)*y(k,27) + rxt(k,316)*y(k,46) + rxt(k,318) &
                      *y(k,14) + rxt(k,321)*y(k,88) + rxt(k,329)*y(k,100) + rxt(k,330) &
                      *y(k,101) + rxt(k,339)*y(k,102) + rxt(k,340)*y(k,103) + rxt(k,341) &
                      *y(k,104) + rxt(k,343)*y(k,106) + rxt(k,346)*y(k,1) + rxt(k,350) &
                      *y(k,2) + rxt(k,351)*y(k,13) + rxt(k,352)*y(k,89) + rxt(k,353) &
                      *y(k,91) + rxt(k,354)*y(k,92) + rxt(k,366)*y(k,94) + rxt(k,367) &
                      *y(k,95) + rxt(k,374)*y(k,97) + rxt(k,376)*y(k,93) + rxt(k,377) &
                      *y(k,98) + rxt(k,378)*y(k,110) + rxt(k,379)*y(k,111) + rxt(k,385) &
                      *y(k,154) + rxt(k,388)*y(k,5) + rxt(k,391)*y(k,6) + rxt(k,392) &
                      *y(k,20) + rxt(k,394)*y(k,21) + rxt(k,398)*y(k,30) + rxt(k,399) &
                      *y(k,64) + rxt(k,411)*y(k,130) + rxt(k,414)*y(k,131) + rxt(k,418) &
                      *y(k,152) + rxt(k,419)*y(k,153) + rxt(k,421)*y(k,155) + rxt(k,424) &
                      *y(k,156) + rxt(k,427)*y(k,157) + rxt(k,428)*y(k,158) + rxt(k,431) &
                      *y(k,4) + rxt(k,434)*y(k,105) + rxt(k,439)*y(k,120) + rxt(k,443) &
                      *y(k,147) + rxt(k,444)*y(k,148) + rxt(k,448)*y(k,149) + rxt(k,450) &
                      *y(k,150) + rxt(k,451)*y(k,151) + (rxt(k,453) + rxt(k,466) &
                      ) * y(k,65) + rxt(k,455)*y(k,125) + rxt(k,460)*y(k,134) &
                      + rxt(k,465)*y(k,136) + rxt(k,467)*y(k,137) + rxt(k,469) &
                      *y(k,112))
         mat(k,1029) = -rxt(k,131)*y(k,190)
         mat(k,467) = -rxt(k,132)*y(k,190)
         mat(k,1682) = -rxt(k,133)*y(k,190)
         mat(k,1576) = -rxt(k,134)*y(k,190)
         mat(k,1858) = -rxt(k,135)*y(k,190)
         mat(k,266) = -rxt(k,139)*y(k,190)
         mat(k,2013) = -rxt(k,151)*y(k,190)
         mat(k,279) = -rxt(k,152)*y(k,190)
         mat(k,1899) = -rxt(k,160)*y(k,190)
         mat(k,1264) = -rxt(k,161)*y(k,190)
         mat(k,828) = -rxt(k,180)*y(k,190)
         mat(k,1708) = -(rxt(k,182) + rxt(k,183)) * y(k,190)
         mat(k,1922) = -rxt(k,185)*y(k,190)
         mat(k,677) = -rxt(k,188)*y(k,190)
         mat(k,1523) = -rxt(k,212)*y(k,190)
         mat(k,698) = -rxt(k,214)*y(k,190)
         mat(k,1546) = -rxt(k,247)*y(k,190)
         mat(k,657) = -rxt(k,252)*y(k,190)
         mat(k,298) = -rxt(k,253)*y(k,190)
         mat(k,905) = -(rxt(k,255) + rxt(k,265)) * y(k,190)
         mat(k,109) = -rxt(k,256)*y(k,190)
         mat(k,672) = -rxt(k,257)*y(k,190)
         mat(k,205) = -rxt(k,267)*y(k,190)
         mat(k,168) = -rxt(k,274)*y(k,190)
         mat(k,216) = -rxt(k,275)*y(k,190)
         mat(k,197) = -rxt(k,277)*y(k,190)
         mat(k,899) = -rxt(k,279)*y(k,190)
         mat(k,58) = -rxt(k,280)*y(k,190)
         mat(k,432) = -rxt(k,285)*y(k,190)
         mat(k,404) = -rxt(k,286)*y(k,190)
         mat(k,866) = -rxt(k,291)*y(k,190)
         mat(k,738) = -rxt(k,292)*y(k,190)
         mat(k,358) = -rxt(k,293)*y(k,190)
         mat(k,446) = -rxt(k,294)*y(k,190)
         mat(k,310) = -rxt(k,302)*y(k,190)
         mat(k,68) = -rxt(k,303)*y(k,190)
         mat(k,1076) = -rxt(k,305)*y(k,190)
         mat(k,928) = -rxt(k,306)*y(k,190)
         mat(k,732) = -rxt(k,307)*y(k,190)
         mat(k,426) = -rxt(k,310)*y(k,190)
         mat(k,304) = -rxt(k,314)*y(k,190)
         mat(k,847) = -rxt(k,315)*y(k,190)
         mat(k,822) = -rxt(k,316)*y(k,190)
         mat(k,261) = -rxt(k,318)*y(k,190)
         mat(k,919) = -rxt(k,321)*y(k,190)
         mat(k,1066) = -rxt(k,329)*y(k,190)
         mat(k,221) = -rxt(k,330)*y(k,190)
         mat(k,400) = -rxt(k,339)*y(k,190)
         mat(k,227) = -rxt(k,340)*y(k,190)
         mat(k,461) = -rxt(k,341)*y(k,190)
         mat(k,1141) = -rxt(k,343)*y(k,190)
         mat(k,538) = -rxt(k,346)*y(k,190)
         mat(k,528) = -rxt(k,350)*y(k,190)
         mat(k,156) = -rxt(k,351)*y(k,190)
         mat(k,152) = -rxt(k,352)*y(k,190)
         mat(k,212) = -rxt(k,353)*y(k,190)
         mat(k,81) = -rxt(k,354)*y(k,190)
         mat(k,478) = -rxt(k,366)*y(k,190)
         mat(k,440) = -rxt(k,367)*y(k,190)
         mat(k,274) = -rxt(k,374)*y(k,190)
         mat(k,716) = -rxt(k,376)*y(k,190)
         mat(k,600) = -rxt(k,377)*y(k,190)
         mat(k,287) = -rxt(k,378)*y(k,190)
         mat(k,879) = -rxt(k,379)*y(k,190)
         mat(k,130) = -rxt(k,385)*y(k,190)
         mat(k,92) = -rxt(k,388)*y(k,190)
         mat(k,317) = -rxt(k,391)*y(k,190)
         mat(k,159) = -rxt(k,392)*y(k,190)
         mat(k,242) = -rxt(k,394)*y(k,190)
         mat(k,173) = -rxt(k,398)*y(k,190)
         mat(k,122) = -rxt(k,399)*y(k,190)
         mat(k,101) = -rxt(k,411)*y(k,190)
         mat(k,236) = -rxt(k,414)*y(k,190)
         mat(k,500) = -rxt(k,418)*y(k,190)
         mat(k,117) = -rxt(k,419)*y(k,190)
         mat(k,139) = -rxt(k,421)*y(k,190)
         mat(k,563) = -rxt(k,424)*y(k,190)
         mat(k,144) = -rxt(k,427)*y(k,190)
         mat(k,329) = -rxt(k,428)*y(k,190)
         mat(k,770) = -rxt(k,431)*y(k,190)
         mat(k,797) = -rxt(k,434)*y(k,190)
         mat(k,293) = -rxt(k,439)*y(k,190)
         mat(k,488) = -rxt(k,443)*y(k,190)
         mat(k,509) = -rxt(k,444)*y(k,190)
         mat(k,370) = -rxt(k,448)*y(k,190)
         mat(k,860) = -rxt(k,450)*y(k,190)
         mat(k,890) = -rxt(k,451)*y(k,190)
         mat(k,186) = -(rxt(k,453) + rxt(k,466)) * y(k,190)
         mat(k,255) = -rxt(k,455)*y(k,190)
         mat(k,393) = -rxt(k,460)*y(k,190)
         mat(k,1086) = -rxt(k,465)*y(k,190)
         mat(k,725) = -rxt(k,467)*y(k,190)
         mat(k,64) = -rxt(k,469)*y(k,190)
         mat(k,770) = mat(k,770) + .630_r8*rxt(k,430)*y(k,122)
         mat(k,205) = mat(k,205) + .650_r8*rxt(k,267)*y(k,190)
         mat(k,446) = mat(k,446) + .130_r8*rxt(k,269)*y(k,122)
         mat(k,216) = mat(k,216) + .500_r8*rxt(k,275)*y(k,190)
         mat(k,847) = mat(k,847) + .360_r8*rxt(k,298)*y(k,122)
         mat(k,1546) = mat(k,1546) + rxt(k,246)*y(k,121)
         mat(k,298) = mat(k,298) + .300_r8*rxt(k,253)*y(k,190)
         mat(k,1956) = rxt(k,169)*y(k,176)
         mat(k,644) = rxt(k,223)*y(k,202)
         mat(k,1276) = rxt(k,130)*y(k,122) + 2.000_r8*rxt(k,125)*y(k,176)
         mat(k,1029) = mat(k,1029) + rxt(k,122)*y(k,121) + rxt(k,114)*y(k,189)
         mat(k,467) = mat(k,467) + rxt(k,123)*y(k,121)
         mat(k,698) = mat(k,698) + rxt(k,213)*y(k,121) + rxt(k,219)*y(k,189)
         mat(k,1922) = mat(k,1922) + rxt(k,184)*y(k,121) + rxt(k,196)*y(k,189)
         mat(k,109) = mat(k,109) + rxt(k,264)*y(k,189)
         mat(k,650) = rxt(k,215)*y(k,121)
         mat(k,677) = mat(k,677) + rxt(k,187)*y(k,121)
         mat(k,716) = mat(k,716) + .320_r8*rxt(k,375)*y(k,122)
         mat(k,600) = mat(k,600) + .600_r8*rxt(k,377)*y(k,190)
         mat(k,1066) = mat(k,1066) + .240_r8*rxt(k,328)*y(k,122)
         mat(k,221) = mat(k,221) + .100_r8*rxt(k,330)*y(k,190)
         mat(k,797) = mat(k,797) + .630_r8*rxt(k,433)*y(k,122)
         mat(k,1141) = mat(k,1141) + .360_r8*rxt(k,342)*y(k,122)
         mat(k,1798) = rxt(k,153)*y(k,176)
         mat(k,2013) = mat(k,2013) + rxt(k,148)*y(k,176)
         mat(k,1576) = mat(k,1576) + rxt(k,246)*y(k,40) + rxt(k,122)*y(k,72) &
                      + rxt(k,123)*y(k,74) + rxt(k,213)*y(k,76) + rxt(k,184)*y(k,80) &
                      + rxt(k,215)*y(k,86) + rxt(k,187)*y(k,87) + rxt(k,128)*y(k,176)
         mat(k,1858) = mat(k,1858) + .630_r8*rxt(k,430)*y(k,4) + .130_r8*rxt(k,269) &
                      *y(k,23) + .360_r8*rxt(k,298)*y(k,27) + rxt(k,130)*y(k,71) &
                      + .320_r8*rxt(k,375)*y(k,93) + .240_r8*rxt(k,328)*y(k,100) &
                      + .630_r8*rxt(k,433)*y(k,105) + .360_r8*rxt(k,342)*y(k,106) &
                      + rxt(k,129)*y(k,176)
         mat(k,426) = mat(k,426) + .500_r8*rxt(k,310)*y(k,190)
         mat(k,130) = mat(k,130) + .500_r8*rxt(k,385)*y(k,190)
         mat(k,410) = .400_r8*rxt(k,386)*y(k,176)
         mat(k,1237) = .450_r8*rxt(k,283)*y(k,176)
         mat(k,624) = .400_r8*rxt(k,400)*y(k,176)
         mat(k,1682) = mat(k,1682) + rxt(k,169)*y(k,54) + 2.000_r8*rxt(k,125)*y(k,71) &
                      + rxt(k,153)*y(k,116) + rxt(k,148)*y(k,118) + rxt(k,128) &
                      *y(k,121) + rxt(k,129)*y(k,122) + .400_r8*rxt(k,386)*y(k,161) &
                      + .450_r8*rxt(k,283)*y(k,170) + .400_r8*rxt(k,400)*y(k,172) &
                      + .450_r8*rxt(k,333)*y(k,184) + .400_r8*rxt(k,406)*y(k,185) &
                      + .200_r8*rxt(k,337)*y(k,186) + .150_r8*rxt(k,312)*y(k,193)
         mat(k,1207) = .450_r8*rxt(k,333)*y(k,176)
         mat(k,745) = .400_r8*rxt(k,406)*y(k,176)
         mat(k,546) = .200_r8*rxt(k,337)*y(k,176)
         mat(k,1345) = rxt(k,114)*y(k,72) + rxt(k,219)*y(k,76) + rxt(k,196)*y(k,80) &
                      + rxt(k,264)*y(k,81) + 2.000_r8*rxt(k,115)*y(k,202)
         mat(k,1499) = mat(k,1499) + .650_r8*rxt(k,267)*y(k,22) + .500_r8*rxt(k,275) &
                      *y(k,25) + .300_r8*rxt(k,253)*y(k,51) + .600_r8*rxt(k,377) &
                      *y(k,98) + .100_r8*rxt(k,330)*y(k,101) + .500_r8*rxt(k,310) &
                      *y(k,132) + .500_r8*rxt(k,385)*y(k,154)
         mat(k,999) = .150_r8*rxt(k,312)*y(k,176)
         mat(k,2038) = rxt(k,223)*y(k,68) + 2.000_r8*rxt(k,115)*y(k,189)
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
         mat(k,346) = -(rxt(k,409)*y(k,176) + rxt(k,410)*y(k,116))
         mat(k,1620) = -rxt(k,409)*y(k,191)
         mat(k,1742) = -rxt(k,410)*y(k,191)
         mat(k,120) = .200_r8*rxt(k,399)*y(k,190)
         mat(k,99) = .140_r8*rxt(k,411)*y(k,190)
         mat(k,234) = rxt(k,414)*y(k,190)
         mat(k,1417) = .200_r8*rxt(k,399)*y(k,64) + .140_r8*rxt(k,411)*y(k,130) &
                      + rxt(k,414)*y(k,131)
         mat(k,629) = -(rxt(k,308)*y(k,176) + rxt(k,309)*y(k,116))
         mat(k,1644) = -rxt(k,308)*y(k,192)
         mat(k,1762) = -rxt(k,309)*y(k,192)
         mat(k,837) = rxt(k,315)*y(k,190)
         mat(k,423) = .500_r8*rxt(k,310)*y(k,190)
         mat(k,1451) = rxt(k,315)*y(k,27) + .500_r8*rxt(k,310)*y(k,132)
         mat(k,995) = -(rxt(k,311)*y(k,171) + rxt(k,312)*y(k,176) + rxt(k,313) &
                      *y(k,116))
         mat(k,1307) = -rxt(k,311)*y(k,193)
         mat(k,1665) = -rxt(k,312)*y(k,193)
         mat(k,1782) = -rxt(k,313)*y(k,193)
         mat(k,767) = .060_r8*rxt(k,430)*y(k,122)
         mat(k,819) = rxt(k,316)*y(k,190)
         mat(k,794) = .060_r8*rxt(k,433)*y(k,122)
         mat(k,1842) = .060_r8*rxt(k,430)*y(k,4) + .060_r8*rxt(k,433)*y(k,105)
         mat(k,302) = rxt(k,314)*y(k,190)
         mat(k,887) = .150_r8*rxt(k,451)*y(k,190)
         mat(k,1481) = rxt(k,316)*y(k,46) + rxt(k,314)*y(k,133) + .150_r8*rxt(k,451) &
                      *y(k,151)
         mat(k,960) = -(rxt(k,440)*y(k,171) + rxt(k,441)*y(k,176) + rxt(k,442) &
                      *y(k,116))
         mat(k,1305) = -rxt(k,440)*y(k,194)
         mat(k,1663) = -rxt(k,441)*y(k,194)
         mat(k,1780) = -rxt(k,442)*y(k,194)
         mat(k,1994) = .500_r8*rxt(k,449)*y(k,150)
         mat(k,487) = rxt(k,443)*y(k,190)
         mat(k,858) = .500_r8*rxt(k,449)*y(k,118) + rxt(k,450)*y(k,190)
         mat(k,1479) = rxt(k,443)*y(k,147) + rxt(k,450)*y(k,150)
         mat(k,938) = -(rxt(k,445)*y(k,171) + rxt(k,446)*y(k,176) + rxt(k,447) &
                      *y(k,116))
         mat(k,1304) = -rxt(k,445)*y(k,195)
         mat(k,1662) = -rxt(k,446)*y(k,195)
         mat(k,1779) = -rxt(k,447)*y(k,195)
         mat(k,765) = rxt(k,431)*y(k,190)
         mat(k,792) = rxt(k,434)*y(k,190)
         mat(k,369) = rxt(k,448)*y(k,190)
         mat(k,1478) = rxt(k,431)*y(k,4) + rxt(k,434)*y(k,105) + rxt(k,448)*y(k,149)
         mat(k,585) = -(rxt(k,416)*y(k,176) + rxt(k,417)*y(k,116))
         mat(k,1640) = -rxt(k,416)*y(k,196)
         mat(k,1759) = -rxt(k,417)*y(k,196)
         mat(k,496) = rxt(k,418)*y(k,190)
         mat(k,116) = .650_r8*rxt(k,419)*y(k,190)
         mat(k,1447) = rxt(k,418)*y(k,152) + .650_r8*rxt(k,419)*y(k,153)
         mat(k,50) = -(rxt(k,507)*y(k,176) + rxt(k,508)*y(k,116))
         mat(k,1599) = -rxt(k,507)*y(k,197)
         mat(k,1731) = -rxt(k,508)*y(k,197)
         mat(k,111) = rxt(k,506)*y(k,190)
         mat(k,1370) = rxt(k,506)*y(k,153)
         mat(k,1011) = -(rxt(k,380)*y(k,170) + rxt(k,381)*y(k,171) + rxt(k,382) &
                      *y(k,176) + rxt(k,383)*y(k,116) + rxt(k,384)*y(k,118))
         mat(k,1224) = -rxt(k,380)*y(k,198)
         mat(k,1308) = -rxt(k,381)*y(k,198)
         mat(k,1666) = -rxt(k,382)*y(k,198)
         mat(k,1783) = -rxt(k,383)*y(k,198)
         mat(k,1997) = -rxt(k,384)*y(k,198)
         mat(k,151) = rxt(k,352)*y(k,190)
         mat(k,211) = rxt(k,353)*y(k,190)
         mat(k,80) = rxt(k,354)*y(k,190)
         mat(k,597) = .400_r8*rxt(k,377)*y(k,190)
         mat(k,129) = .500_r8*rxt(k,385)*y(k,190)
         mat(k,1482) = rxt(k,352)*y(k,89) + rxt(k,353)*y(k,91) + rxt(k,354)*y(k,92) &
                      + .400_r8*rxt(k,377)*y(k,98) + .500_r8*rxt(k,385)*y(k,154)
         mat(k,609) = -(rxt(k,422)*y(k,176) + rxt(k,423)*y(k,116))
         mat(k,1642) = -rxt(k,422)*y(k,199)
         mat(k,1760) = -rxt(k,423)*y(k,199)
         mat(k,136) = .560_r8*rxt(k,421)*y(k,190)
         mat(k,556) = rxt(k,424)*y(k,190)
         mat(k,1449) = .560_r8*rxt(k,421)*y(k,155) + rxt(k,424)*y(k,156)
         mat(k,56) = -(rxt(k,510)*y(k,176) + rxt(k,511)*y(k,116))
         mat(k,1600) = -rxt(k,510)*y(k,200)
         mat(k,1732) = -rxt(k,511)*y(k,200)
         mat(k,131) = rxt(k,509)*y(k,190)
         mat(k,1371) = rxt(k,509)*y(k,155)
         mat(k,383) = -(rxt(k,425)*y(k,176) + rxt(k,426)*y(k,116))
         mat(k,1625) = -rxt(k,425)*y(k,201)
         mat(k,1746) = -rxt(k,426)*y(k,201)
         mat(k,143) = .300_r8*rxt(k,427)*y(k,190)
         mat(k,326) = rxt(k,428)*y(k,190)
         mat(k,1423) = .300_r8*rxt(k,427)*y(k,157) + rxt(k,428)*y(k,158)
         mat(k,2050) = -(rxt(k,115)*y(k,189) + rxt(k,223)*y(k,68) + rxt(k,468) &
                      *y(k,138))
         mat(k,1357) = -rxt(k,115)*y(k,202)
         mat(k,646) = -rxt(k,223)*y(k,202)
         mat(k,178) = -rxt(k,468)*y(k,202)
         mat(k,200) = rxt(k,277)*y(k,190)
         mat(k,312) = rxt(k,302)*y(k,190)
         mat(k,69) = rxt(k,303)*y(k,190)
         mat(k,1558) = rxt(k,247)*y(k,190)
         mat(k,902) = rxt(k,279)*y(k,190)
         mat(k,823) = rxt(k,316)*y(k,190)
         mat(k,1079) = rxt(k,305)*y(k,190)
         mat(k,433) = rxt(k,285)*y(k,190)
         mat(k,406) = rxt(k,286)*y(k,190)
         mat(k,300) = rxt(k,253)*y(k,190)
         mat(k,1284) = rxt(k,126)*y(k,176)
         mat(k,1034) = rxt(k,131)*y(k,190)
         mat(k,472) = rxt(k,132)*y(k,190)
         mat(k,701) = rxt(k,214)*y(k,190)
         mat(k,1934) = (rxt(k,520)+rxt(k,525))*y(k,86) + (rxt(k,513)+rxt(k,519) &
                       +rxt(k,524))*y(k,87) + rxt(k,185)*y(k,190)
         mat(k,674) = rxt(k,257)*y(k,190)
         mat(k,1270) = rxt(k,161)*y(k,190)
         mat(k,270) = rxt(k,139)*y(k,190)
         mat(k,655) = (rxt(k,520)+rxt(k,525))*y(k,80)
         mat(k,682) = (rxt(k,513)+rxt(k,519)+rxt(k,524))*y(k,80) + rxt(k,188)*y(k,190)
         mat(k,1070) = .500_r8*rxt(k,329)*y(k,190)
         mat(k,65) = rxt(k,469)*y(k,190)
         mat(k,429) = rxt(k,310)*y(k,190)
         mat(k,306) = rxt(k,314)*y(k,190)
         mat(k,1694) = rxt(k,126)*y(k,71) + rxt(k,133)*y(k,190)
         mat(k,1511) = rxt(k,277)*y(k,26) + rxt(k,302)*y(k,28) + rxt(k,303)*y(k,29) &
                      + rxt(k,247)*y(k,40) + rxt(k,279)*y(k,43) + rxt(k,316)*y(k,46) &
                      + rxt(k,305)*y(k,47) + rxt(k,285)*y(k,48) + rxt(k,286)*y(k,49) &
                      + rxt(k,253)*y(k,51) + rxt(k,131)*y(k,72) + rxt(k,132)*y(k,74) &
                      + rxt(k,214)*y(k,76) + rxt(k,185)*y(k,80) + rxt(k,257)*y(k,82) &
                      + rxt(k,161)*y(k,84) + rxt(k,139)*y(k,85) + rxt(k,188)*y(k,87) &
                      + .500_r8*rxt(k,329)*y(k,100) + rxt(k,469)*y(k,112) + rxt(k,310) &
                      *y(k,132) + rxt(k,314)*y(k,133) + rxt(k,133)*y(k,176) &
                      + 2.000_r8*rxt(k,136)*y(k,190)
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
         mat(k, 4) = lmat(k, 4)
         mat(k, 5) = lmat(k, 5)
         mat(k, 6) = lmat(k, 6)
         mat(k, 12) = mat(k, 12) + lmat(k, 12)
         mat(k, 18) = mat(k, 18) + lmat(k, 18)
         mat(k, 24) = mat(k, 24) + lmat(k, 24)
         mat(k, 30) = mat(k, 30) + lmat(k, 30)
         mat(k, 36) = mat(k, 36) + lmat(k, 36)
         mat(k, 38) = mat(k, 38) + lmat(k, 38)
         mat(k, 44) = mat(k, 44) + lmat(k, 44)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 60) = lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 71) = mat(k, 71) + lmat(k, 71)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = lmat(k, 76)
         mat(k, 77) = lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 82) = lmat(k, 82)
         mat(k, 83) = lmat(k, 83)
         mat(k, 84) = lmat(k, 84)
         mat(k, 85) = lmat(k, 85)
         mat(k, 86) = lmat(k, 86)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 106) = mat(k, 106) + lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 146) = lmat(k, 146)
         mat(k, 147) = lmat(k, 147)
         mat(k, 148) = lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 160) = lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 163) = lmat(k, 163)
         mat(k, 164) = lmat(k, 164)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = mat(k, 166) + lmat(k, 166)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = lmat(k, 194)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 201) = mat(k, 201) + lmat(k, 201)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 213) = mat(k, 213) + lmat(k, 213)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 237) = lmat(k, 237)
         mat(k, 238) = lmat(k, 238)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 280) = mat(k, 280) + lmat(k, 280)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 286) = lmat(k, 286)
         mat(k, 289) = mat(k, 289) + lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 295) = mat(k, 295) + lmat(k, 295)
         mat(k, 296) = lmat(k, 296)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 303) = lmat(k, 303)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 307) = mat(k, 307) + lmat(k, 307)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 313) = mat(k, 313) + lmat(k, 313)
         mat(k, 314) = lmat(k, 314)
         mat(k, 316) = lmat(k, 316)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = lmat(k, 328)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 330) = lmat(k, 330)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 352) = lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 356) = lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 360) = lmat(k, 360)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 367) = lmat(k, 367)
         mat(k, 368) = lmat(k, 368)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 383) = mat(k, 383) + lmat(k, 383)
         mat(k, 390) = mat(k, 390) + lmat(k, 390)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 394) = lmat(k, 394)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 403) = lmat(k, 403)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 414) = mat(k, 414) + lmat(k, 414)
         mat(k, 415) = lmat(k, 415)
         mat(k, 416) = lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = mat(k, 418) + lmat(k, 418)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 424) = lmat(k, 424)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = lmat(k, 428)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 434) = mat(k, 434) + lmat(k, 434)
         mat(k, 439) = lmat(k, 439)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 460) = lmat(k, 460)
         mat(k, 464) = lmat(k, 464)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 477) = lmat(k, 477)
         mat(k, 482) = mat(k, 482) + lmat(k, 482)
         mat(k, 483) = lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 498) = lmat(k, 498)
         mat(k, 499) = lmat(k, 499)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 501) = lmat(k, 501)
         mat(k, 502) = lmat(k, 502)
         mat(k, 503) = lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = mat(k, 506) + lmat(k, 506)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 510) = lmat(k, 510)
         mat(k, 511) = mat(k, 511) + lmat(k, 511)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 525) = lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 540) = lmat(k, 540)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 550) = lmat(k, 550)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = lmat(k, 552)
         mat(k, 553) = lmat(k, 553)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 558) = lmat(k, 558)
         mat(k, 561) = lmat(k, 561)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 564) = lmat(k, 564)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 585) = mat(k, 585) + lmat(k, 585)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 598) = lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 601) = lmat(k, 601)
         mat(k, 602) = lmat(k, 602)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 629) = mat(k, 629) + lmat(k, 629)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 639) = mat(k, 639) + lmat(k, 639)
         mat(k, 643) = lmat(k, 643)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 696) = mat(k, 696) + lmat(k, 696)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 723) = mat(k, 723) + lmat(k, 723)
         mat(k, 724) = lmat(k, 724)
         mat(k, 726) = lmat(k, 726)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 731) = lmat(k, 731)
         mat(k, 733) = lmat(k, 733)
         mat(k, 734) = mat(k, 734) + lmat(k, 734)
         mat(k, 735) = lmat(k, 735)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 737) = mat(k, 737) + lmat(k, 737)
         mat(k, 739) = mat(k, 739) + lmat(k, 739)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 786) = mat(k, 786) + lmat(k, 786)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 818) = mat(k, 818) + lmat(k, 818)
         mat(k, 820) = lmat(k, 820)
         mat(k, 821) = lmat(k, 821)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 830) = mat(k, 830) + lmat(k, 830)
         mat(k, 831) = lmat(k, 831)
         mat(k, 833) = mat(k, 833) + lmat(k, 833)
         mat(k, 834) = mat(k, 834) + lmat(k, 834)
         mat(k, 840) = mat(k, 840) + lmat(k, 840)
         mat(k, 855) = mat(k, 855) + lmat(k, 855)
         mat(k, 856) = lmat(k, 856)
         mat(k, 857) = lmat(k, 857)
         mat(k, 861) = lmat(k, 861)
         mat(k, 864) = mat(k, 864) + lmat(k, 864)
         mat(k, 865) = lmat(k, 865)
         mat(k, 867) = mat(k, 867) + lmat(k, 867)
         mat(k, 868) = mat(k, 868) + lmat(k, 868)
         mat(k, 869) = lmat(k, 869)
         mat(k, 873) = mat(k, 873) + lmat(k, 873)
         mat(k, 877) = lmat(k, 877)
         mat(k, 881) = mat(k, 881) + lmat(k, 881)
         mat(k, 883) = lmat(k, 883)
         mat(k, 884) = mat(k, 884) + lmat(k, 884)
         mat(k, 885) = mat(k, 885) + lmat(k, 885)
         mat(k, 886) = mat(k, 886) + lmat(k, 886)
         mat(k, 887) = mat(k, 887) + lmat(k, 887)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 891) = mat(k, 891) + lmat(k, 891)
         mat(k, 892) = mat(k, 892) + lmat(k, 892)
         mat(k, 894) = mat(k, 894) + lmat(k, 894)
         mat(k, 895) = lmat(k, 895)
         mat(k, 898) = lmat(k, 898)
         mat(k, 900) = lmat(k, 900)
         mat(k, 903) = mat(k, 903) + lmat(k, 903)
         mat(k, 908) = lmat(k, 908)
         mat(k, 909) = lmat(k, 909)
         mat(k, 910) = lmat(k, 910)
         mat(k, 911) = lmat(k, 911)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 913) = lmat(k, 913)
         mat(k, 915) = lmat(k, 915)
         mat(k, 916) = lmat(k, 916)
         mat(k, 920) = lmat(k, 920)
         mat(k, 921) = mat(k, 921) + lmat(k, 921)
         mat(k, 922) = lmat(k, 922)
         mat(k, 925) = mat(k, 925) + lmat(k, 925)
         mat(k, 927) = lmat(k, 927)
         mat(k, 929) = lmat(k, 929)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 938) = mat(k, 938) + lmat(k, 938)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 995) = mat(k, 995) + lmat(k, 995)
         mat(k,1011) = mat(k,1011) + lmat(k,1011)
         mat(k,1024) = mat(k,1024) + lmat(k,1024)
         mat(k,1044) = mat(k,1044) + lmat(k,1044)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1067) = mat(k,1067) + lmat(k,1067)
         mat(k,1068) = mat(k,1068) + lmat(k,1068)
         mat(k,1071) = mat(k,1071) + lmat(k,1071)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1073) = mat(k,1073) + lmat(k,1073)
         mat(k,1077) = lmat(k,1077)
         mat(k,1081) = lmat(k,1081)
         mat(k,1082) = mat(k,1082) + lmat(k,1082)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1088) = lmat(k,1088)
         mat(k,1096) = lmat(k,1096)
         mat(k,1113) = mat(k,1113) + lmat(k,1113)
         mat(k,1123) = mat(k,1123) + lmat(k,1123)
         mat(k,1130) = lmat(k,1130)
         mat(k,1132) = mat(k,1132) + lmat(k,1132)
         mat(k,1135) = mat(k,1135) + lmat(k,1135)
         mat(k,1137) = mat(k,1137) + lmat(k,1137)
         mat(k,1140) = lmat(k,1140)
         mat(k,1157) = mat(k,1157) + lmat(k,1157)
         mat(k,1183) = mat(k,1183) + lmat(k,1183)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1247) = mat(k,1247) + lmat(k,1247)
         mat(k,1260) = mat(k,1260) + lmat(k,1260)
         mat(k,1264) = mat(k,1264) + lmat(k,1264)
         mat(k,1266) = lmat(k,1266)
         mat(k,1273) = mat(k,1273) + lmat(k,1273)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1320) = mat(k,1320) + lmat(k,1320)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1335) = mat(k,1335) + lmat(k,1335)
         mat(k,1337) = mat(k,1337) + lmat(k,1337)
         mat(k,1339) = mat(k,1339) + lmat(k,1339)
         mat(k,1340) = mat(k,1340) + lmat(k,1340)
         mat(k,1342) = mat(k,1342) + lmat(k,1342)
         mat(k,1343) = lmat(k,1343)
         mat(k,1344) = mat(k,1344) + lmat(k,1344)
         mat(k,1345) = mat(k,1345) + lmat(k,1345)
         mat(k,1347) = lmat(k,1347)
         mat(k,1348) = lmat(k,1348)
         mat(k,1349) = lmat(k,1349)
         mat(k,1351) = lmat(k,1351)
         mat(k,1355) = mat(k,1355) + lmat(k,1355)
         mat(k,1375) = lmat(k,1375)
         mat(k,1380) = lmat(k,1380)
         mat(k,1494) = mat(k,1494) + lmat(k,1494)
         mat(k,1497) = mat(k,1497) + lmat(k,1497)
         mat(k,1499) = mat(k,1499) + lmat(k,1499)
         mat(k,1503) = mat(k,1503) + lmat(k,1503)
         mat(k,1509) = mat(k,1509) + lmat(k,1509)
         mat(k,1511) = mat(k,1511) + lmat(k,1511)
         mat(k,1518) = mat(k,1518) + lmat(k,1518)
         mat(k,1524) = mat(k,1524) + lmat(k,1524)
         mat(k,1526) = mat(k,1526) + lmat(k,1526)
         mat(k,1539) = mat(k,1539) + lmat(k,1539)
         mat(k,1540) = lmat(k,1540)
         mat(k,1543) = mat(k,1543) + lmat(k,1543)
         mat(k,1548) = mat(k,1548) + lmat(k,1548)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1583) = mat(k,1583) + lmat(k,1583)
         mat(k,1686) = mat(k,1686) + lmat(k,1686)
         mat(k,1694) = mat(k,1694) + lmat(k,1694)
         mat(k,1711) = mat(k,1711) + lmat(k,1711)
         mat(k,1713) = mat(k,1713) + lmat(k,1713)
         mat(k,1718) = mat(k,1718) + lmat(k,1718)
         mat(k,1739) = mat(k,1739) + lmat(k,1739)
         mat(k,1801) = mat(k,1801) + lmat(k,1801)
         mat(k,1804) = mat(k,1804) + lmat(k,1804)
         mat(k,1857) = mat(k,1857) + lmat(k,1857)
         mat(k,1861) = mat(k,1861) + lmat(k,1861)
         mat(k,1865) = mat(k,1865) + lmat(k,1865)
         mat(k,1895) = mat(k,1895) + lmat(k,1895)
         mat(k,1899) = mat(k,1899) + lmat(k,1899)
         mat(k,1902) = mat(k,1902) + lmat(k,1902)
         mat(k,1905) = mat(k,1905) + lmat(k,1905)
         mat(k,1907) = mat(k,1907) + lmat(k,1907)
         mat(k,1919) = mat(k,1919) + lmat(k,1919)
         mat(k,1931) = mat(k,1931) + lmat(k,1931)
         mat(k,1932) = mat(k,1932) + lmat(k,1932)
         mat(k,1948) = mat(k,1948) + lmat(k,1948)
         mat(k,1951) = lmat(k,1951)
         mat(k,1954) = lmat(k,1954)
         mat(k,1960) = mat(k,1960) + lmat(k,1960)
         mat(k,1965) = mat(k,1965) + lmat(k,1965)
         mat(k,1966) = mat(k,1966) + lmat(k,1966)
         mat(k,2009) = mat(k,2009) + lmat(k,2009)
         mat(k,2016) = mat(k,2016) + lmat(k,2016)
         mat(k,2019) = mat(k,2019) + lmat(k,2019)
         mat(k,2021) = mat(k,2021) + lmat(k,2021)
         mat(k,2024) = mat(k,2024) + lmat(k,2024)
         mat(k,2031) = lmat(k,2031)
         mat(k,2035) = lmat(k,2035)
         mat(k,2037) = mat(k,2037) + lmat(k,2037)
         mat(k,2038) = mat(k,2038) + lmat(k,2038)
         mat(k,2041) = lmat(k,2041)
         mat(k,2050) = mat(k,2050) + lmat(k,2050)
         mat(k, 137) = 0._r8
         mat(k, 138) = 0._r8
         mat(k, 241) = 0._r8
         mat(k, 334) = 0._r8
         mat(k, 335) = 0._r8
         mat(k, 348) = 0._r8
         mat(k, 376) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 386) = 0._r8
         mat(k, 495) = 0._r8
         mat(k, 497) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 518) = 0._r8
         mat(k, 522) = 0._r8
         mat(k, 523) = 0._r8
         mat(k, 527) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 534) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 559) = 0._r8
         mat(k, 560) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 584) = 0._r8
         mat(k, 586) = 0._r8
         mat(k, 588) = 0._r8
         mat(k, 589) = 0._r8
         mat(k, 591) = 0._r8
         mat(k, 608) = 0._r8
         mat(k, 610) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 613) = 0._r8
         mat(k, 615) = 0._r8
         mat(k, 616) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 654) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 728) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 768) = 0._r8
         mat(k, 775) = 0._r8
         mat(k, 785) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 806) = 0._r8
         mat(k, 807) = 0._r8
         mat(k, 811) = 0._r8
         mat(k, 812) = 0._r8
         mat(k, 813) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 852) = 0._r8
         mat(k, 854) = 0._r8
         mat(k, 872) = 0._r8
         mat(k, 874) = 0._r8
         mat(k, 875) = 0._r8
         mat(k, 876) = 0._r8
         mat(k, 878) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 882) = 0._r8
         mat(k, 889) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 914) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 923) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 940) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 977) = 0._r8
         mat(k, 978) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 983) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 991) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1045) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1052) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1118) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1234) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1532) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1545) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1566) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1574) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1588) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1621) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1656) = 0._r8
         mat(k,1657) = 0._r8
         mat(k,1658) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1672) = 0._r8
         mat(k,1678) = 0._r8
         mat(k,1681) = 0._r8
         mat(k,1704) = 0._r8
         mat(k,1705) = 0._r8
         mat(k,1707) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1794) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1797) = 0._r8
         mat(k,1807) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1827) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1834) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1840) = 0._r8
         mat(k,1841) = 0._r8
         mat(k,1843) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1849) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1851) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1867) = 0._r8
         mat(k,1870) = 0._r8
         mat(k,1880) = 0._r8
         mat(k,1883) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1888) = 0._r8
         mat(k,1889) = 0._r8
         mat(k,1890) = 0._r8
         mat(k,1894) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1897) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1901) = 0._r8
         mat(k,1908) = 0._r8
         mat(k,1909) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1924) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1940) = 0._r8
         mat(k,1941) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1964) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,2001) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2023) = 0._r8
         mat(k,2025) = 0._r8
         mat(k,2030) = 0._r8
         mat(k,2032) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2036) = 0._r8
         mat(k,2039) = 0._r8
         mat(k,2040) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2044) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2049) = 0._r8
         mat(k, 1) = mat(k, 1) - dti(k)
         mat(k, 2) = mat(k, 2) - dti(k)
         mat(k, 3) = mat(k, 3) - dti(k)
         mat(k, 4) = mat(k, 4) - dti(k)
         mat(k, 5) = mat(k, 5) - dti(k)
         mat(k, 6) = mat(k, 6) - dti(k)
         mat(k, 12) = mat(k, 12) - dti(k)
         mat(k, 18) = mat(k, 18) - dti(k)
         mat(k, 24) = mat(k, 24) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 56) = mat(k, 56) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 63) = mat(k, 63) - dti(k)
         mat(k, 66) = mat(k, 66) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 82) = mat(k, 82) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 88) = mat(k, 88) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 146) = mat(k, 146) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 195) = mat(k, 195) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 207) = mat(k, 207) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 213) = mat(k, 213) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 223) = mat(k, 223) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 249) = mat(k, 249) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 277) = mat(k, 277) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 289) = mat(k, 289) - dti(k)
         mat(k, 295) = mat(k, 295) - dti(k)
         mat(k, 301) = mat(k, 301) - dti(k)
         mat(k, 307) = mat(k, 307) - dti(k)
         mat(k, 313) = mat(k, 313) - dti(k)
         mat(k, 319) = mat(k, 319) - dti(k)
         mat(k, 325) = mat(k, 325) - dti(k)
         mat(k, 333) = mat(k, 333) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 346) = mat(k, 346) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 355) = mat(k, 355) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 366) = mat(k, 366) - dti(k)
         mat(k, 375) = mat(k, 375) - dti(k)
         mat(k, 383) = mat(k, 383) - dti(k)
         mat(k, 390) = mat(k, 390) - dti(k)
         mat(k, 396) = mat(k, 396) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 408) = mat(k, 408) - dti(k)
         mat(k, 414) = mat(k, 414) - dti(k)
         mat(k, 422) = mat(k, 422) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 434) = mat(k, 434) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 450) = mat(k, 450) - dti(k)
         mat(k, 458) = mat(k, 458) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 482) = mat(k, 482) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 506) = mat(k, 506) - dti(k)
         mat(k, 514) = mat(k, 514) - dti(k)
         mat(k, 521) = mat(k, 521) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 543) = mat(k, 543) - dti(k)
         mat(k, 554) = mat(k, 554) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 585) = mat(k, 585) - dti(k)
         mat(k, 596) = mat(k, 596) - dti(k)
         mat(k, 609) = mat(k, 609) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 629) = mat(k, 629) - dti(k)
         mat(k, 639) = mat(k, 639) - dti(k)
         mat(k, 648) = mat(k, 648) - dti(k)
         mat(k, 656) = mat(k, 656) - dti(k)
         mat(k, 661) = mat(k, 661) - dti(k)
         mat(k, 671) = mat(k, 671) - dti(k)
         mat(k, 676) = mat(k, 676) - dti(k)
         mat(k, 686) = mat(k, 686) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 706) = mat(k, 706) - dti(k)
         mat(k, 723) = mat(k, 723) - dti(k)
         mat(k, 729) = mat(k, 729) - dti(k)
         mat(k, 736) = mat(k, 736) - dti(k)
         mat(k, 741) = mat(k, 741) - dti(k)
         mat(k, 759) = mat(k, 759) - dti(k)
         mat(k, 786) = mat(k, 786) - dti(k)
         mat(k, 808) = mat(k, 808) - dti(k)
         mat(k, 818) = mat(k, 818) - dti(k)
         mat(k, 826) = mat(k, 826) - dti(k)
         mat(k, 840) = mat(k, 840) - dti(k)
         mat(k, 855) = mat(k, 855) - dti(k)
         mat(k, 864) = mat(k, 864) - dti(k)
         mat(k, 873) = mat(k, 873) - dti(k)
         mat(k, 885) = mat(k, 885) - dti(k)
         mat(k, 894) = mat(k, 894) - dti(k)
         mat(k, 903) = mat(k, 903) - dti(k)
         mat(k, 912) = mat(k, 912) - dti(k)
         mat(k, 925) = mat(k, 925) - dti(k)
         mat(k, 938) = mat(k, 938) - dti(k)
         mat(k, 960) = mat(k, 960) - dti(k)
         mat(k, 979) = mat(k, 979) - dti(k)
         mat(k, 995) = mat(k, 995) - dti(k)
         mat(k,1011) = mat(k,1011) - dti(k)
         mat(k,1024) = mat(k,1024) - dti(k)
         mat(k,1044) = mat(k,1044) - dti(k)
         mat(k,1060) = mat(k,1060) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1083) = mat(k,1083) - dti(k)
         mat(k,1113) = mat(k,1113) - dti(k)
         mat(k,1135) = mat(k,1135) - dti(k)
         mat(k,1157) = mat(k,1157) - dti(k)
         mat(k,1183) = mat(k,1183) - dti(k)
         mat(k,1202) = mat(k,1202) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1247) = mat(k,1247) - dti(k)
         mat(k,1260) = mat(k,1260) - dti(k)
         mat(k,1273) = mat(k,1273) - dti(k)
         mat(k,1320) = mat(k,1320) - dti(k)
         mat(k,1344) = mat(k,1344) - dti(k)
         mat(k,1499) = mat(k,1499) - dti(k)
         mat(k,1524) = mat(k,1524) - dti(k)
         mat(k,1548) = mat(k,1548) - dti(k)
         mat(k,1579) = mat(k,1579) - dti(k)
         mat(k,1686) = mat(k,1686) - dti(k)
         mat(k,1713) = mat(k,1713) - dti(k)
         mat(k,1804) = mat(k,1804) - dti(k)
         mat(k,1865) = mat(k,1865) - dti(k)
         mat(k,1907) = mat(k,1907) - dti(k)
         mat(k,1931) = mat(k,1931) - dti(k)
         mat(k,1966) = mat(k,1966) - dti(k)
         mat(k,2024) = mat(k,2024) - dti(k)
         mat(k,2050) = mat(k,2050) - dti(k)
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
