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
         mat(k,536) = -(rxt(k,396)*y(k,221))
         mat(k,1657) = -rxt(k,396)*y(k,1)
         mat(k,1462) = rxt(k,399)*y(k,190)
         mat(k,876) = rxt(k,399)*y(k,124)
         mat(k,547) = -(rxt(k,400)*y(k,221))
         mat(k,1658) = -rxt(k,400)*y(k,2)
         mat(k,877) = rxt(k,397)*y(k,203)
         mat(k,2057) = rxt(k,397)*y(k,190)
         mat(k,858) = -(rxt(k,479)*y(k,126) + rxt(k,480)*y(k,135) + rxt(k,481) &
                      *y(k,221))
         mat(k,1869) = -rxt(k,479)*y(k,6)
         mat(k,1750) = -rxt(k,480)*y(k,6)
         mat(k,1682) = -rxt(k,481)*y(k,6)
         mat(k,86) = -(rxt(k,438)*y(k,221))
         mat(k,1595) = -rxt(k,438)*y(k,7)
         mat(k,273) = -(rxt(k,441)*y(k,221))
         mat(k,1625) = -rxt(k,441)*y(k,8)
         mat(k,375) = rxt(k,439)*y(k,203)
         mat(k,2031) = rxt(k,439)*y(k,191)
         mat(k,87) = .120_r8*rxt(k,438)*y(k,221)
         mat(k,1596) = .120_r8*rxt(k,438)*y(k,7)
         mat(k,855) = .100_r8*rxt(k,480)*y(k,135)
         mat(k,817) = .100_r8*rxt(k,483)*y(k,135)
         mat(k,1739) = .100_r8*rxt(k,480)*y(k,6) + .100_r8*rxt(k,483)*y(k,110)
         mat(k,1449) = .500_r8*rxt(k,440)*y(k,191) + .200_r8*rxt(k,467)*y(k,228) &
                      + .060_r8*rxt(k,473)*y(k,230)
         mat(k,376) = .500_r8*rxt(k,440)*y(k,124)
         mat(k,614) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,630) = .060_r8*rxt(k,473)*y(k,124)
         mat(k,1443) = .200_r8*rxt(k,467)*y(k,228) + .200_r8*rxt(k,473)*y(k,230)
         mat(k,613) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,628) = .200_r8*rxt(k,473)*y(k,124)
         mat(k,1458) = .200_r8*rxt(k,467)*y(k,228) + .150_r8*rxt(k,473)*y(k,230)
         mat(k,616) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,631) = .150_r8*rxt(k,473)*y(k,124)
         mat(k,1445) = .210_r8*rxt(k,473)*y(k,230)
         mat(k,629) = .210_r8*rxt(k,473)*y(k,124)
         mat(k,163) = -(rxt(k,401)*y(k,221))
         mat(k,1608) = -rxt(k,401)*y(k,15)
         mat(k,854) = .050_r8*rxt(k,480)*y(k,135)
         mat(k,816) = .050_r8*rxt(k,483)*y(k,135)
         mat(k,1738) = .050_r8*rxt(k,480)*y(k,6) + .050_r8*rxt(k,483)*y(k,110)
         mat(k,259) = -(rxt(k,367)*y(k,126) + rxt(k,368)*y(k,221))
         mat(k,1863) = -rxt(k,367)*y(k,16)
         mat(k,1623) = -rxt(k,368)*y(k,16)
         mat(k,1352) = -(rxt(k,250)*y(k,42) + rxt(k,251)*y(k,203) + rxt(k,252) &
                      *y(k,135))
         mat(k,1840) = -rxt(k,250)*y(k,17)
         mat(k,2100) = -rxt(k,251)*y(k,17)
         mat(k,1775) = -rxt(k,252)*y(k,17)
         mat(k,2006) = 4.000_r8*rxt(k,253)*y(k,19) + (rxt(k,254)+rxt(k,255))*y(k,59) &
                      + rxt(k,258)*y(k,124) + rxt(k,261)*y(k,133) + rxt(k,508) &
                      *y(k,151) + rxt(k,262)*y(k,221)
         mat(k,2127) = (rxt(k,254)+rxt(k,255))*y(k,19)
         mat(k,754) = rxt(k,263)*y(k,133) + rxt(k,269)*y(k,217) + rxt(k,264)*y(k,221)
         mat(k,1504) = rxt(k,258)*y(k,19)
         mat(k,1816) = rxt(k,261)*y(k,19) + rxt(k,263)*y(k,81)
         mat(k,1319) = rxt(k,508)*y(k,19)
         mat(k,1565) = rxt(k,269)*y(k,81)
         mat(k,1713) = rxt(k,262)*y(k,19) + rxt(k,264)*y(k,81)
         mat(k,1999) = rxt(k,256)*y(k,59)
         mat(k,2120) = rxt(k,256)*y(k,19)
         mat(k,1333) = (rxt(k,556)+rxt(k,561))*y(k,91)
         mat(k,653) = (rxt(k,556)+rxt(k,561))*y(k,85)
         mat(k,2019) = -(4._r8*rxt(k,253)*y(k,19) + (rxt(k,254) + rxt(k,255) + rxt(k,256) &
                      ) * y(k,59) + rxt(k,257)*y(k,203) + rxt(k,258)*y(k,124) &
                      + rxt(k,259)*y(k,125) + rxt(k,261)*y(k,133) + rxt(k,262) &
                      *y(k,221) + rxt(k,508)*y(k,151))
         mat(k,2141) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,19)
         mat(k,2114) = -rxt(k,257)*y(k,19)
         mat(k,1518) = -rxt(k,258)*y(k,19)
         mat(k,1995) = -rxt(k,259)*y(k,19)
         mat(k,1830) = -rxt(k,261)*y(k,19)
         mat(k,1727) = -rxt(k,262)*y(k,19)
         mat(k,1328) = -rxt(k,508)*y(k,19)
         mat(k,1360) = rxt(k,252)*y(k,135)
         mat(k,445) = rxt(k,260)*y(k,133)
         mat(k,759) = rxt(k,270)*y(k,217)
         mat(k,660) = rxt(k,265)*y(k,133)
         mat(k,1830) = mat(k,1830) + rxt(k,260)*y(k,20) + rxt(k,265)*y(k,91)
         mat(k,1789) = rxt(k,252)*y(k,17)
         mat(k,1579) = rxt(k,270)*y(k,81)
         mat(k,438) = -(rxt(k,260)*y(k,133))
         mat(k,1798) = -rxt(k,260)*y(k,20)
         mat(k,2001) = rxt(k,259)*y(k,125)
         mat(k,1963) = rxt(k,259)*y(k,19)
         mat(k,169) = -(rxt(k,442)*y(k,221))
         mat(k,1609) = -rxt(k,442)*y(k,22)
         mat(k,1442) = rxt(k,445)*y(k,192)
         mat(k,321) = rxt(k,445)*y(k,124)
         mat(k,241) = -(rxt(k,444)*y(k,221))
         mat(k,1620) = -rxt(k,444)*y(k,23)
         mat(k,322) = rxt(k,443)*y(k,203)
         mat(k,2029) = rxt(k,443)*y(k,192)
         mat(k,200) = -(rxt(k,316)*y(k,56) + rxt(k,317)*y(k,221))
         mat(k,1523) = -rxt(k,316)*y(k,24)
         mat(k,1614) = -rxt(k,317)*y(k,24)
         mat(k,450) = -(rxt(k,318)*y(k,56) + rxt(k,319)*y(k,135) + rxt(k,344)*y(k,221))
         mat(k,1525) = -rxt(k,318)*y(k,25)
         mat(k,1742) = -rxt(k,319)*y(k,25)
         mat(k,1647) = -rxt(k,344)*y(k,25)
         mat(k,177) = -(rxt(k,324)*y(k,221))
         mat(k,1611) = -rxt(k,324)*y(k,26)
         mat(k,837) = .800_r8*rxt(k,320)*y(k,193) + .200_r8*rxt(k,321)*y(k,197)
         mat(k,1363) = .200_r8*rxt(k,321)*y(k,193)
         mat(k,246) = -(rxt(k,325)*y(k,221))
         mat(k,1621) = -rxt(k,325)*y(k,27)
         mat(k,838) = rxt(k,322)*y(k,203)
         mat(k,2030) = rxt(k,322)*y(k,193)
         mat(k,206) = -(rxt(k,326)*y(k,56) + rxt(k,327)*y(k,221))
         mat(k,1524) = -rxt(k,326)*y(k,28)
         mat(k,1615) = -rxt(k,327)*y(k,28)
         mat(k,943) = -(rxt(k,347)*y(k,126) + rxt(k,348)*y(k,135) + rxt(k,365) &
                      *y(k,221))
         mat(k,1875) = -rxt(k,347)*y(k,29)
         mat(k,1755) = -rxt(k,348)*y(k,29)
         mat(k,1689) = -rxt(k,365)*y(k,29)
         mat(k,782) = .130_r8*rxt(k,425)*y(k,135)
         mat(k,1755) = mat(k,1755) + .130_r8*rxt(k,425)*y(k,98)
         mat(k,303) = -(rxt(k,352)*y(k,221))
         mat(k,1629) = -rxt(k,352)*y(k,30)
         mat(k,727) = rxt(k,350)*y(k,203)
         mat(k,2035) = rxt(k,350)*y(k,194)
         mat(k,55) = -(rxt(k,353)*y(k,221))
         mat(k,1592) = -rxt(k,353)*y(k,31)
         mat(k,181) = -(rxt(k,448)*y(k,221))
         mat(k,1612) = -rxt(k,448)*y(k,32)
         mat(k,527) = rxt(k,446)*y(k,203)
         mat(k,2025) = rxt(k,446)*y(k,195)
         mat(k,1849) = -(rxt(k,214)*y(k,56) + rxt(k,250)*y(k,17) + rxt(k,294)*y(k,203) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,133) + rxt(k,297) &
                      *y(k,221))
         mat(k,1548) = -rxt(k,214)*y(k,42)
         mat(k,1358) = -rxt(k,250)*y(k,42)
         mat(k,2109) = -rxt(k,294)*y(k,42)
         mat(k,1906) = -rxt(k,295)*y(k,42)
         mat(k,1825) = -rxt(k,296)*y(k,42)
         mat(k,1722) = -rxt(k,297)*y(k,42)
         mat(k,543) = .400_r8*rxt(k,396)*y(k,221)
         mat(k,871) = .340_r8*rxt(k,480)*y(k,135)
         mat(k,264) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,456) = rxt(k,319)*y(k,135)
         mat(k,954) = .500_r8*rxt(k,348)*y(k,135)
         mat(k,410) = .500_r8*rxt(k,336)*y(k,221)
         mat(k,710) = rxt(k,302)*y(k,221)
         mat(k,312) = .300_r8*rxt(k,303)*y(k,221)
         mat(k,2136) = rxt(k,221)*y(k,197)
         mat(k,963) = .800_r8*rxt(k,341)*y(k,221)
         mat(k,793) = .910_r8*rxt(k,425)*y(k,135)
         mat(k,510) = .300_r8*rxt(k,416)*y(k,221)
         mat(k,1132) = .800_r8*rxt(k,420)*y(k,197)
         mat(k,1147) = .120_r8*rxt(k,378)*y(k,135)
         mat(k,470) = .500_r8*rxt(k,391)*y(k,221)
         mat(k,833) = .340_r8*rxt(k,483)*y(k,135)
         mat(k,1259) = .600_r8*rxt(k,392)*y(k,135)
         mat(k,1513) = .100_r8*rxt(k,398)*y(k,190) + rxt(k,301)*y(k,197) &
                      + .500_r8*rxt(k,369)*y(k,200) + .500_r8*rxt(k,338)*y(k,202) &
                      + .920_r8*rxt(k,408)*y(k,205) + .250_r8*rxt(k,376)*y(k,207) &
                      + rxt(k,385)*y(k,209) + rxt(k,359)*y(k,224) + rxt(k,363) &
                      *y(k,225) + .340_r8*rxt(k,492)*y(k,226) + .320_r8*rxt(k,497) &
                      *y(k,227) + .250_r8*rxt(k,433)*y(k,229)
         mat(k,1906) = mat(k,1906) + .500_r8*rxt(k,367)*y(k,16) + rxt(k,409)*y(k,205) &
                      + .250_r8*rxt(k,375)*y(k,207) + rxt(k,386)*y(k,209)
         mat(k,1784) = .340_r8*rxt(k,480)*y(k,6) + rxt(k,319)*y(k,25) &
                      + .500_r8*rxt(k,348)*y(k,29) + .910_r8*rxt(k,425)*y(k,98) &
                      + .120_r8*rxt(k,378)*y(k,105) + .340_r8*rxt(k,483)*y(k,110) &
                      + .600_r8*rxt(k,392)*y(k,111)
         mat(k,354) = rxt(k,343)*y(k,221)
         mat(k,988) = .680_r8*rxt(k,501)*y(k,221)
         mat(k,888) = .100_r8*rxt(k,398)*y(k,124)
         mat(k,846) = .700_r8*rxt(k,321)*y(k,197)
         mat(k,735) = rxt(k,349)*y(k,197)
         mat(k,1308) = rxt(k,332)*y(k,197) + rxt(k,405)*y(k,205) + .250_r8*rxt(k,372) &
                      *y(k,207) + rxt(k,381)*y(k,209) + .250_r8*rxt(k,430)*y(k,229)
         mat(k,1402) = rxt(k,221)*y(k,59) + .800_r8*rxt(k,420)*y(k,101) + rxt(k,301) &
                      *y(k,124) + .700_r8*rxt(k,321)*y(k,193) + rxt(k,349)*y(k,194) &
                      + rxt(k,332)*y(k,196) + (4.000_r8*rxt(k,298)+2.000_r8*rxt(k,299)) &
                      *y(k,197) + 1.500_r8*rxt(k,406)*y(k,205) + .750_r8*rxt(k,411) &
                      *y(k,206) + .880_r8*rxt(k,373)*y(k,207) + 2.000_r8*rxt(k,382) &
                      *y(k,209) + .750_r8*rxt(k,485)*y(k,216) + .800_r8*rxt(k,361) &
                      *y(k,225) + .930_r8*rxt(k,490)*y(k,226) + .950_r8*rxt(k,495) &
                      *y(k,227) + .800_r8*rxt(k,431)*y(k,229)
         mat(k,463) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,669) = .500_r8*rxt(k,338)*y(k,124)
         mat(k,2109) = mat(k,2109) + .450_r8*rxt(k,383)*y(k,209) + .150_r8*rxt(k,362) &
                      *y(k,225)
         mat(k,1181) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) + rxt(k,405) &
                      *y(k,196) + 1.500_r8*rxt(k,406)*y(k,197)
         mat(k,1215) = .750_r8*rxt(k,411)*y(k,197)
         mat(k,1237) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,196) + .880_r8*rxt(k,373)*y(k,197)
         mat(k,1277) = rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126) + rxt(k,381)*y(k,196) &
                      + 2.000_r8*rxt(k,382)*y(k,197) + .450_r8*rxt(k,383)*y(k,203) &
                      + 4.000_r8*rxt(k,384)*y(k,209)
         mat(k,1055) = .750_r8*rxt(k,485)*y(k,197)
         mat(k,1722) = mat(k,1722) + .400_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,336) &
                      *y(k,51) + rxt(k,302)*y(k,52) + .300_r8*rxt(k,303)*y(k,53) &
                      + .800_r8*rxt(k,341)*y(k,74) + .300_r8*rxt(k,416)*y(k,99) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,343)*y(k,140) &
                      + .680_r8*rxt(k,501)*y(k,179)
         mat(k,723) = rxt(k,359)*y(k,124)
         mat(k,1071) = rxt(k,363)*y(k,124) + .800_r8*rxt(k,361)*y(k,197) &
                      + .150_r8*rxt(k,362)*y(k,203)
         mat(k,1036) = .340_r8*rxt(k,492)*y(k,124) + .930_r8*rxt(k,490)*y(k,197)
         mat(k,1016) = .320_r8*rxt(k,497)*y(k,124) + .950_r8*rxt(k,495)*y(k,197)
         mat(k,1097) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,430)*y(k,196) &
                      + .800_r8*rxt(k,431)*y(k,197)
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
         mat(k,1076) = -(rxt(k,328)*y(k,126) + rxt(k,329)*y(k,221))
         mat(k,1885) = -rxt(k,328)*y(k,45)
         mat(k,1699) = -rxt(k,329)*y(k,45)
         mat(k,540) = .800_r8*rxt(k,396)*y(k,221)
         mat(k,262) = rxt(k,367)*y(k,126)
         mat(k,178) = rxt(k,324)*y(k,221)
         mat(k,248) = .500_r8*rxt(k,325)*y(k,221)
         mat(k,946) = .500_r8*rxt(k,348)*y(k,135)
         mat(k,1248) = .100_r8*rxt(k,392)*y(k,135)
         mat(k,1493) = .400_r8*rxt(k,398)*y(k,190) + rxt(k,323)*y(k,193) &
                      + .270_r8*rxt(k,351)*y(k,194) + rxt(k,369)*y(k,200) + rxt(k,388) &
                      *y(k,211) + rxt(k,359)*y(k,224)
         mat(k,1885) = mat(k,1885) + rxt(k,367)*y(k,16)
         mat(k,1764) = .500_r8*rxt(k,348)*y(k,29) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,882) = .400_r8*rxt(k,398)*y(k,124)
         mat(k,841) = rxt(k,323)*y(k,124) + 3.200_r8*rxt(k,320)*y(k,193) &
                      + .800_r8*rxt(k,321)*y(k,197)
         mat(k,730) = .270_r8*rxt(k,351)*y(k,124)
         mat(k,1385) = .800_r8*rxt(k,321)*y(k,193)
         mat(k,461) = rxt(k,369)*y(k,124)
         mat(k,2087) = .200_r8*rxt(k,387)*y(k,211)
         mat(k,559) = rxt(k,388)*y(k,124) + .200_r8*rxt(k,387)*y(k,203)
         mat(k,1699) = mat(k,1699) + .800_r8*rxt(k,396)*y(k,1) + rxt(k,324)*y(k,26) &
                      + .500_r8*rxt(k,325)*y(k,27)
         mat(k,719) = rxt(k,359)*y(k,124)
         mat(k,52) = -(rxt(k,330)*y(k,221))
         mat(k,1591) = -rxt(k,330)*y(k,47)
         mat(k,898) = -(rxt(k,366)*y(k,221))
         mat(k,1685) = -rxt(k,366)*y(k,48)
         mat(k,539) = .800_r8*rxt(k,396)*y(k,221)
         mat(k,860) = .520_r8*rxt(k,480)*y(k,135)
         mat(k,261) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,822) = .520_r8*rxt(k,483)*y(k,135)
         mat(k,1481) = .250_r8*rxt(k,398)*y(k,190) + .820_r8*rxt(k,351)*y(k,194) &
                      + .500_r8*rxt(k,369)*y(k,200) + .270_r8*rxt(k,492)*y(k,226) &
                      + .040_r8*rxt(k,497)*y(k,227)
         mat(k,1872) = .500_r8*rxt(k,367)*y(k,16)
         mat(k,1753) = .520_r8*rxt(k,480)*y(k,6) + .520_r8*rxt(k,483)*y(k,110)
         mat(k,981) = .500_r8*rxt(k,501)*y(k,221)
         mat(k,881) = .250_r8*rxt(k,398)*y(k,124)
         mat(k,729) = .820_r8*rxt(k,351)*y(k,124) + .820_r8*rxt(k,349)*y(k,197)
         mat(k,1374) = .820_r8*rxt(k,349)*y(k,194) + .150_r8*rxt(k,490)*y(k,226) &
                      + .025_r8*rxt(k,495)*y(k,227)
         mat(k,459) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,1685) = mat(k,1685) + .800_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,501) &
                      *y(k,179)
         mat(k,1026) = .270_r8*rxt(k,492)*y(k,124) + .150_r8*rxt(k,490)*y(k,197)
         mat(k,1004) = .040_r8*rxt(k,497)*y(k,124) + .025_r8*rxt(k,495)*y(k,197)
         mat(k,1152) = -(rxt(k,354)*y(k,126) + rxt(k,355)*y(k,221))
         mat(k,1889) = -rxt(k,354)*y(k,49)
         mat(k,1704) = -rxt(k,355)*y(k,49)
         mat(k,992) = rxt(k,356)*y(k,221)
         mat(k,1141) = .880_r8*rxt(k,378)*y(k,135)
         mat(k,1249) = .500_r8*rxt(k,392)*y(k,135)
         mat(k,1497) = .170_r8*rxt(k,451)*y(k,198) + .050_r8*rxt(k,414)*y(k,206) &
                      + .250_r8*rxt(k,376)*y(k,207) + .170_r8*rxt(k,457)*y(k,210) &
                      + .400_r8*rxt(k,467)*y(k,228) + .250_r8*rxt(k,433)*y(k,229) &
                      + .540_r8*rxt(k,473)*y(k,230) + .510_r8*rxt(k,476)*y(k,231)
         mat(k,1889) = mat(k,1889) + .050_r8*rxt(k,415)*y(k,206) + .250_r8*rxt(k,375) &
                      *y(k,207) + .250_r8*rxt(k,434)*y(k,229)
         mat(k,770) = rxt(k,357)*y(k,221)
         mat(k,1767) = .880_r8*rxt(k,378)*y(k,105) + .500_r8*rxt(k,392)*y(k,111)
         mat(k,1296) = .250_r8*rxt(k,372)*y(k,207) + .250_r8*rxt(k,430)*y(k,229)
         mat(k,1389) = .240_r8*rxt(k,373)*y(k,207) + .500_r8*rxt(k,361)*y(k,225) &
                      + .100_r8*rxt(k,431)*y(k,229)
         mat(k,647) = .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450)*y(k,203)
         mat(k,2092) = .070_r8*rxt(k,450)*y(k,198) + .070_r8*rxt(k,456)*y(k,210)
         mat(k,1205) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1230) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,196) + .240_r8*rxt(k,373)*y(k,197)
         mat(k,805) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,203)
         mat(k,1704) = mat(k,1704) + rxt(k,356)*y(k,95) + rxt(k,357)*y(k,127)
         mat(k,1066) = .500_r8*rxt(k,361)*y(k,197)
         mat(k,623) = .400_r8*rxt(k,467)*y(k,124)
         mat(k,1092) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,196) + .100_r8*rxt(k,431)*y(k,197)
         mat(k,639) = .540_r8*rxt(k,473)*y(k,124)
         mat(k,394) = .510_r8*rxt(k,476)*y(k,124)
         mat(k,446) = -(rxt(k,335)*y(k,221))
         mat(k,1646) = -rxt(k,335)*y(k,50)
         mat(k,939) = .120_r8*rxt(k,348)*y(k,135)
         mat(k,1741) = .120_r8*rxt(k,348)*y(k,29)
         mat(k,1287) = .100_r8*rxt(k,332)*y(k,197) + .150_r8*rxt(k,333)*y(k,203)
         mat(k,1367) = .100_r8*rxt(k,332)*y(k,196)
         mat(k,2051) = .150_r8*rxt(k,333)*y(k,196) + .150_r8*rxt(k,383)*y(k,209)
         mat(k,1267) = .150_r8*rxt(k,383)*y(k,203)
         mat(k,406) = -(rxt(k,336)*y(k,221))
         mat(k,1642) = -rxt(k,336)*y(k,51)
         mat(k,1286) = .400_r8*rxt(k,333)*y(k,203)
         mat(k,2049) = .400_r8*rxt(k,333)*y(k,196) + .400_r8*rxt(k,383)*y(k,209)
         mat(k,1266) = .400_r8*rxt(k,383)*y(k,203)
         mat(k,708) = -(rxt(k,302)*y(k,221))
         mat(k,1669) = -rxt(k,302)*y(k,52)
         mat(k,1118) = .200_r8*rxt(k,420)*y(k,197)
         mat(k,839) = .300_r8*rxt(k,321)*y(k,197)
         mat(k,1369) = .200_r8*rxt(k,420)*y(k,101) + .300_r8*rxt(k,321)*y(k,193) &
                      + 2.000_r8*rxt(k,299)*y(k,197) + .250_r8*rxt(k,406)*y(k,205) &
                      + .250_r8*rxt(k,411)*y(k,206) + .250_r8*rxt(k,373)*y(k,207) &
                      + .250_r8*rxt(k,485)*y(k,216) + .500_r8*rxt(k,361)*y(k,225) &
                      + .250_r8*rxt(k,490)*y(k,226) + .250_r8*rxt(k,495)*y(k,227) &
                      + .300_r8*rxt(k,431)*y(k,229)
         mat(k,1162) = .250_r8*rxt(k,406)*y(k,197)
         mat(k,1193) = .250_r8*rxt(k,411)*y(k,197)
         mat(k,1223) = .250_r8*rxt(k,373)*y(k,197)
         mat(k,1044) = .250_r8*rxt(k,485)*y(k,197)
         mat(k,1063) = .500_r8*rxt(k,361)*y(k,197)
         mat(k,1025) = .250_r8*rxt(k,490)*y(k,197)
         mat(k,1003) = .250_r8*rxt(k,495)*y(k,197)
         mat(k,1086) = .300_r8*rxt(k,431)*y(k,197)
         mat(k,309) = -(rxt(k,303)*y(k,221))
         mat(k,1630) = -rxt(k,303)*y(k,53)
         mat(k,1366) = rxt(k,300)*y(k,203)
         mat(k,2036) = rxt(k,300)*y(k,197)
         mat(k,1543) = -(rxt(k,214)*y(k,42) + rxt(k,216)*y(k,77) + rxt(k,217)*y(k,79) &
                      + (rxt(k,218) + rxt(k,219)) * y(k,203) + rxt(k,220)*y(k,135) &
                      + rxt(k,227)*y(k,60) + rxt(k,236)*y(k,92) + rxt(k,326)*y(k,28))
         mat(k,1844) = -rxt(k,214)*y(k,56)
         mat(k,1107) = -rxt(k,216)*y(k,56)
         mat(k,476) = -rxt(k,217)*y(k,56)
         mat(k,2104) = -(rxt(k,218) + rxt(k,219)) * y(k,56)
         mat(k,1779) = -rxt(k,220)*y(k,56)
         mat(k,908) = -rxt(k,227)*y(k,56)
         mat(k,764) = -rxt(k,236)*y(k,56)
         mat(k,209) = -rxt(k,326)*y(k,56)
         mat(k,2009) = rxt(k,255)*y(k,59)
         mat(k,2131) = rxt(k,255)*y(k,19) + (4.000_r8*rxt(k,222)+2.000_r8*rxt(k,224)) &
                      *y(k,59) + rxt(k,226)*y(k,124) + rxt(k,231)*y(k,133) &
                      + rxt(k,509)*y(k,151) + rxt(k,221)*y(k,197) + rxt(k,232) &
                      *y(k,221)
         mat(k,103) = rxt(k,276)*y(k,217)
         mat(k,1339) = rxt(k,234)*y(k,133) + rxt(k,246)*y(k,217) + rxt(k,235)*y(k,221)
         mat(k,1508) = rxt(k,226)*y(k,59)
         mat(k,1820) = rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85)
         mat(k,1322) = rxt(k,509)*y(k,59)
         mat(k,1399) = rxt(k,221)*y(k,59)
         mat(k,1569) = rxt(k,276)*y(k,65) + rxt(k,246)*y(k,85)
         mat(k,1717) = rxt(k,232)*y(k,59) + rxt(k,235)*y(k,85)
         mat(k,1522) = rxt(k,227)*y(k,60)
         mat(k,2119) = 2.000_r8*rxt(k,223)*y(k,59)
         mat(k,904) = rxt(k,227)*y(k,56) + (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,85)
         mat(k,1332) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,60) + (rxt(k,549) &
                       +rxt(k,555)+rxt(k,560))*y(k,92)
         mat(k,761) = (rxt(k,549)+rxt(k,555)+rxt(k,560))*y(k,85)
         mat(k,2118) = 2.000_r8*rxt(k,248)*y(k,59)
         mat(k,2143) = -(rxt(k,221)*y(k,197) + (4._r8*rxt(k,222) + 4._r8*rxt(k,223) &
                      + 4._r8*rxt(k,224) + 4._r8*rxt(k,248)) * y(k,59) + rxt(k,225) &
                      *y(k,203) + rxt(k,226)*y(k,124) + rxt(k,228)*y(k,125) + rxt(k,231) &
                      *y(k,133) + (rxt(k,232) + rxt(k,233)) * y(k,221) + (rxt(k,254) &
                      + rxt(k,255) + rxt(k,256)) * y(k,19) + rxt(k,509)*y(k,151))
         mat(k,1408) = -rxt(k,221)*y(k,59)
         mat(k,2116) = -rxt(k,225)*y(k,59)
         mat(k,1520) = -rxt(k,226)*y(k,59)
         mat(k,1997) = -rxt(k,228)*y(k,59)
         mat(k,1832) = -rxt(k,231)*y(k,59)
         mat(k,1729) = -(rxt(k,232) + rxt(k,233)) * y(k,59)
         mat(k,2021) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,59)
         mat(k,1330) = -rxt(k,509)*y(k,59)
         mat(k,1555) = rxt(k,236)*y(k,92) + rxt(k,220)*y(k,135) + rxt(k,219)*y(k,203)
         mat(k,914) = rxt(k,229)*y(k,133)
         mat(k,1348) = rxt(k,247)*y(k,217)
         mat(k,767) = rxt(k,236)*y(k,56) + rxt(k,237)*y(k,133) + rxt(k,238)*y(k,221)
         mat(k,1832) = mat(k,1832) + rxt(k,229)*y(k,60) + rxt(k,237)*y(k,92)
         mat(k,1791) = rxt(k,220)*y(k,56)
         mat(k,234) = rxt(k,514)*y(k,151)
         mat(k,1330) = mat(k,1330) + rxt(k,514)*y(k,137)
         mat(k,2116) = mat(k,2116) + rxt(k,219)*y(k,56)
         mat(k,1581) = rxt(k,247)*y(k,85)
         mat(k,1729) = mat(k,1729) + rxt(k,238)*y(k,92)
         mat(k,906) = -(rxt(k,227)*y(k,56) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,221) &
                      + (rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,85))
         mat(k,1532) = -rxt(k,227)*y(k,60)
         mat(k,1811) = -rxt(k,229)*y(k,60)
         mat(k,1686) = -rxt(k,230)*y(k,60)
         mat(k,1336) = -(rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,60)
         mat(k,2124) = rxt(k,228)*y(k,125)
         mat(k,1973) = rxt(k,228)*y(k,59)
         mat(k,997) = -((rxt(k,305) + rxt(k,315)) * y(k,221))
         mat(k,1694) = -(rxt(k,305) + rxt(k,315)) * y(k,62)
         mat(k,863) = .230_r8*rxt(k,480)*y(k,135)
         mat(k,1351) = rxt(k,250)*y(k,42)
         mat(k,203) = .350_r8*rxt(k,317)*y(k,221)
         mat(k,453) = .630_r8*rxt(k,319)*y(k,135)
         mat(k,945) = .560_r8*rxt(k,348)*y(k,135)
         mat(k,1837) = rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) + rxt(k,295)*y(k,126) &
                      + rxt(k,296)*y(k,133) + rxt(k,297)*y(k,221)
         mat(k,1151) = rxt(k,354)*y(k,126) + rxt(k,355)*y(k,221)
         mat(k,1534) = rxt(k,214)*y(k,42)
         mat(k,799) = rxt(k,342)*y(k,221)
         mat(k,783) = .620_r8*rxt(k,425)*y(k,135)
         mat(k,1139) = .650_r8*rxt(k,378)*y(k,135)
         mat(k,825) = .230_r8*rxt(k,483)*y(k,135)
         mat(k,1247) = .560_r8*rxt(k,392)*y(k,135)
         mat(k,1488) = .170_r8*rxt(k,451)*y(k,198) + .220_r8*rxt(k,376)*y(k,207) &
                      + .400_r8*rxt(k,454)*y(k,208) + .350_r8*rxt(k,457)*y(k,210) &
                      + .225_r8*rxt(k,492)*y(k,226) + .250_r8*rxt(k,433)*y(k,229)
         mat(k,1880) = rxt(k,295)*y(k,42) + rxt(k,354)*y(k,49) + .220_r8*rxt(k,375) &
                      *y(k,207) + .500_r8*rxt(k,434)*y(k,229)
         mat(k,1812) = rxt(k,296)*y(k,42) + rxt(k,504)*y(k,138)
         mat(k,1759) = .230_r8*rxt(k,480)*y(k,6) + .630_r8*rxt(k,319)*y(k,25) &
                      + .560_r8*rxt(k,348)*y(k,29) + .620_r8*rxt(k,425)*y(k,98) &
                      + .650_r8*rxt(k,378)*y(k,105) + .230_r8*rxt(k,483)*y(k,110) &
                      + .560_r8*rxt(k,392)*y(k,111)
         mat(k,254) = rxt(k,504)*y(k,133) + rxt(k,505)*y(k,221)
         mat(k,983) = .700_r8*rxt(k,501)*y(k,221)
         mat(k,1292) = .220_r8*rxt(k,372)*y(k,207) + .250_r8*rxt(k,430)*y(k,229)
         mat(k,1380) = .110_r8*rxt(k,373)*y(k,207) + .125_r8*rxt(k,490)*y(k,226) &
                      + .200_r8*rxt(k,431)*y(k,229)
         mat(k,646) = .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450)*y(k,203)
         mat(k,2082) = .070_r8*rxt(k,450)*y(k,198) + .160_r8*rxt(k,453)*y(k,208) &
                      + .140_r8*rxt(k,456)*y(k,210)
         mat(k,1227) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,196) + .110_r8*rxt(k,373)*y(k,197)
         mat(k,609) = .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453)*y(k,203)
         mat(k,804) = .350_r8*rxt(k,457)*y(k,124) + .140_r8*rxt(k,456)*y(k,203)
         mat(k,1694) = mat(k,1694) + .350_r8*rxt(k,317)*y(k,24) + rxt(k,297)*y(k,42) &
                      + rxt(k,355)*y(k,49) + rxt(k,342)*y(k,75) + rxt(k,505)*y(k,138) &
                      + .700_r8*rxt(k,501)*y(k,179)
         mat(k,1029) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,197)
         mat(k,1090) = .250_r8*rxt(k,433)*y(k,124) + .500_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,196) + .200_r8*rxt(k,431)*y(k,197)
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
         mat(k,59) = -(rxt(k,275)*y(k,217))
         mat(k,1557) = -rxt(k,275)*y(k,64)
         mat(k,101) = -(rxt(k,276)*y(k,217))
         mat(k,1559) = -rxt(k,276)*y(k,65)
         mat(k,121) = -(rxt(k,449)*y(k,221))
         mat(k,1601) = -rxt(k,449)*y(k,66)
         mat(k,115) = .180_r8*rxt(k,469)*y(k,221)
         mat(k,1601) = mat(k,1601) + .180_r8*rxt(k,469)*y(k,181)
         mat(k,191) = -(rxt(k,502)*y(k,126) + (rxt(k,503) + rxt(k,516)) * y(k,221))
         mat(k,1861) = -rxt(k,502)*y(k,67)
         mat(k,1613) = -(rxt(k,503) + rxt(k,516)) * y(k,67)
         mat(k,662) = rxt(k,337)*y(k,203)
         mat(k,2023) = rxt(k,337)*y(k,202)
         mat(k,673) = -(rxt(k,272)*y(k,77) + rxt(k,273)*y(k,232) + rxt(k,274)*y(k,89))
         mat(k,1103) = -rxt(k,272)*y(k,73)
         mat(k,2148) = -rxt(k,273)*y(k,73)
         mat(k,1936) = -rxt(k,274)*y(k,73)
         mat(k,60) = 2.000_r8*rxt(k,275)*y(k,217)
         mat(k,102) = rxt(k,276)*y(k,217)
         mat(k,1561) = 2.000_r8*rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65)
         mat(k,960) = -(rxt(k,341)*y(k,221))
         mat(k,1690) = -rxt(k,341)*y(k,74)
         mat(k,506) = .700_r8*rxt(k,416)*y(k,221)
         mat(k,424) = .500_r8*rxt(k,417)*y(k,221)
         mat(k,269) = rxt(k,428)*y(k,221)
         mat(k,1484) = .050_r8*rxt(k,414)*y(k,206) + .530_r8*rxt(k,376)*y(k,207) &
                      + .225_r8*rxt(k,492)*y(k,226) + .250_r8*rxt(k,433)*y(k,229)
         mat(k,1876) = .050_r8*rxt(k,415)*y(k,206) + .530_r8*rxt(k,375)*y(k,207) &
                      + .250_r8*rxt(k,434)*y(k,229)
         mat(k,1423) = rxt(k,340)*y(k,201)
         mat(k,1290) = .530_r8*rxt(k,372)*y(k,207) + .250_r8*rxt(k,430)*y(k,229)
         mat(k,1377) = .260_r8*rxt(k,373)*y(k,207) + .125_r8*rxt(k,490)*y(k,226) &
                      + .100_r8*rxt(k,431)*y(k,229)
         mat(k,346) = rxt(k,340)*y(k,134)
         mat(k,1197) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1224) = .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375)*y(k,126) &
                      + .530_r8*rxt(k,372)*y(k,196) + .260_r8*rxt(k,373)*y(k,197)
         mat(k,1690) = mat(k,1690) + .700_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
                      *y(k,100) + rxt(k,428)*y(k,115)
         mat(k,1027) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,197)
         mat(k,1088) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,196) + .100_r8*rxt(k,431)*y(k,197)
         mat(k,798) = -(rxt(k,342)*y(k,221))
         mat(k,1678) = -rxt(k,342)*y(k,75)
         mat(k,202) = .650_r8*rxt(k,317)*y(k,221)
         mat(k,959) = .200_r8*rxt(k,341)*y(k,221)
         mat(k,926) = rxt(k,429)*y(k,221)
         mat(k,1477) = rxt(k,440)*y(k,191) + .050_r8*rxt(k,414)*y(k,206) &
                      + .400_r8*rxt(k,454)*y(k,208) + .170_r8*rxt(k,457)*y(k,210) &
                      + .700_r8*rxt(k,460)*y(k,223) + .600_r8*rxt(k,467)*y(k,228) &
                      + .250_r8*rxt(k,433)*y(k,229) + .340_r8*rxt(k,473)*y(k,230) &
                      + .170_r8*rxt(k,476)*y(k,231)
         mat(k,1867) = .050_r8*rxt(k,415)*y(k,206) + .250_r8*rxt(k,434)*y(k,229)
         mat(k,379) = rxt(k,440)*y(k,124)
         mat(k,1288) = .250_r8*rxt(k,430)*y(k,229)
         mat(k,1372) = .100_r8*rxt(k,431)*y(k,229)
         mat(k,2073) = .160_r8*rxt(k,453)*y(k,208) + .070_r8*rxt(k,456)*y(k,210)
         mat(k,1195) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,608) = .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453)*y(k,203)
         mat(k,802) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,203)
         mat(k,1678) = mat(k,1678) + .650_r8*rxt(k,317)*y(k,24) + .200_r8*rxt(k,341) &
                      *y(k,74) + rxt(k,429)*y(k,116)
         mat(k,337) = .700_r8*rxt(k,460)*y(k,124)
         mat(k,620) = .600_r8*rxt(k,467)*y(k,124)
         mat(k,1087) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,196) + .100_r8*rxt(k,431)*y(k,197)
         mat(k,636) = .340_r8*rxt(k,473)*y(k,124)
         mat(k,393) = .170_r8*rxt(k,476)*y(k,124)
         mat(k,1928) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,203) + rxt(k,175) &
                      *y(k,134) + rxt(k,178)*y(k,135))
         mat(k,2111) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76)
         mat(k,1435) = -rxt(k,175)*y(k,76)
         mat(k,1786) = -rxt(k,178)*y(k,76)
         mat(k,1851) = rxt(k,297)*y(k,221)
         mat(k,1550) = rxt(k,216)*y(k,77)
         mat(k,999) = rxt(k,315)*y(k,221)
         mat(k,678) = rxt(k,272)*y(k,77)
         mat(k,1112) = rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73) + rxt(k,170)*y(k,133) &
                      + rxt(k,153)*y(k,217) + rxt(k,179)*y(k,221)
         mat(k,758) = rxt(k,270)*y(k,217)
         mat(k,1344) = rxt(k,247)*y(k,217)
         mat(k,751) = rxt(k,202)*y(k,221)
         mat(k,1827) = rxt(k,170)*y(k,77) + rxt(k,182)*y(k,221)
         mat(k,258) = rxt(k,505)*y(k,221)
         mat(k,606) = rxt(k,510)*y(k,221)
         mat(k,1326) = rxt(k,515)*y(k,221)
         mat(k,1576) = rxt(k,153)*y(k,77) + rxt(k,270)*y(k,81) + rxt(k,247)*y(k,85)
         mat(k,1724) = rxt(k,297)*y(k,42) + rxt(k,315)*y(k,62) + rxt(k,179)*y(k,77) &
                      + rxt(k,202)*y(k,112) + rxt(k,182)*y(k,133) + rxt(k,505) &
                      *y(k,138) + rxt(k,510)*y(k,149) + rxt(k,515)*y(k,151)
         mat(k,1104) = -(rxt(k,153)*y(k,217) + rxt(k,170)*y(k,133) + rxt(k,179) &
                      *y(k,221) + rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73))
         mat(k,1563) = -rxt(k,153)*y(k,77)
         mat(k,1813) = -rxt(k,170)*y(k,77)
         mat(k,1701) = -rxt(k,179)*y(k,77)
         mat(k,1536) = -rxt(k,216)*y(k,77)
         mat(k,674) = -rxt(k,272)*y(k,77)
         mat(k,1915) = rxt(k,172)*y(k,203)
         mat(k,2089) = rxt(k,172)*y(k,76)
         mat(k,474) = -(rxt(k,171)*y(k,133) + rxt(k,180)*y(k,221) + rxt(k,217)*y(k,56))
         mat(k,1799) = -rxt(k,171)*y(k,79)
         mat(k,1650) = -rxt(k,180)*y(k,79)
         mat(k,1526) = -rxt(k,217)*y(k,79)
         mat(k,2052) = 2.000_r8*rxt(k,186)*y(k,203)
         mat(k,1650) = mat(k,1650) + 2.000_r8*rxt(k,185)*y(k,221)
         mat(k,172) = rxt(k,518)*y(k,232)
         mat(k,2145) = rxt(k,518)*y(k,153)
         mat(k,753) = -(rxt(k,263)*y(k,133) + rxt(k,264)*y(k,221) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,217))
         mat(k,1808) = -rxt(k,263)*y(k,81)
         mat(k,1674) = -rxt(k,264)*y(k,81)
         mat(k,1562) = -(rxt(k,269) + rxt(k,270)) * y(k,81)
         mat(k,1350) = rxt(k,250)*y(k,42) + rxt(k,251)*y(k,203)
         mat(k,1836) = rxt(k,250)*y(k,17)
         mat(k,2070) = rxt(k,251)*y(k,17)
         mat(k,1337) = -(rxt(k,234)*y(k,133) + rxt(k,235)*y(k,221) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,217) + (rxt(k,549) + rxt(k,555) + rxt(k,560) &
                      ) * y(k,92) + (rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,60) &
                      + (rxt(k,556) + rxt(k,561)) * y(k,91))
         mat(k,1815) = -rxt(k,234)*y(k,85)
         mat(k,1712) = -rxt(k,235)*y(k,85)
         mat(k,1564) = -(rxt(k,246) + rxt(k,247)) * y(k,85)
         mat(k,763) = -(rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,85)
         mat(k,907) = -(rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,85)
         mat(k,655) = -(rxt(k,556) + rxt(k,561)) * y(k,85)
         mat(k,208) = rxt(k,326)*y(k,56)
         mat(k,1839) = rxt(k,214)*y(k,56)
         mat(k,1538) = rxt(k,326)*y(k,28) + rxt(k,214)*y(k,42) + rxt(k,216)*y(k,77) &
                      + rxt(k,217)*y(k,79) + rxt(k,236)*y(k,92) + rxt(k,218)*y(k,203)
         mat(k,2126) = rxt(k,233)*y(k,221)
         mat(k,1105) = rxt(k,216)*y(k,56)
         mat(k,475) = rxt(k,217)*y(k,56)
         mat(k,763) = mat(k,763) + rxt(k,236)*y(k,56)
         mat(k,2099) = rxt(k,218)*y(k,56)
         mat(k,1712) = mat(k,1712) + rxt(k,233)*y(k,59)
         mat(k,105) = -(rxt(k,306)*y(k,221) + rxt(k,314)*y(k,217))
         mat(k,1599) = -rxt(k,306)*y(k,86)
         mat(k,1560) = -rxt(k,314)*y(k,86)
         mat(k,712) = -(rxt(k,307)*y(k,221))
         mat(k,1670) = -rxt(k,307)*y(k,87)
         mat(k,856) = .050_r8*rxt(k,480)*y(k,135)
         mat(k,201) = .350_r8*rxt(k,317)*y(k,221)
         mat(k,452) = .370_r8*rxt(k,319)*y(k,135)
         mat(k,940) = .120_r8*rxt(k,348)*y(k,135)
         mat(k,780) = .110_r8*rxt(k,425)*y(k,135)
         mat(k,1138) = .330_r8*rxt(k,378)*y(k,135)
         mat(k,818) = .050_r8*rxt(k,483)*y(k,135)
         mat(k,1244) = .120_r8*rxt(k,392)*y(k,135)
         mat(k,1472) = rxt(k,310)*y(k,204)
         mat(k,1746) = .050_r8*rxt(k,480)*y(k,6) + .370_r8*rxt(k,319)*y(k,25) &
                      + .120_r8*rxt(k,348)*y(k,29) + .110_r8*rxt(k,425)*y(k,98) &
                      + .330_r8*rxt(k,378)*y(k,105) + .050_r8*rxt(k,483)*y(k,110) &
                      + .120_r8*rxt(k,392)*y(k,111)
         mat(k,2067) = rxt(k,308)*y(k,204)
         mat(k,330) = rxt(k,310)*y(k,124) + rxt(k,308)*y(k,203)
         mat(k,1670) = mat(k,1670) + .350_r8*rxt(k,317)*y(k,24)
         mat(k,672) = rxt(k,272)*y(k,77) + rxt(k,274)*y(k,89) + rxt(k,273)*y(k,232)
         mat(k,1102) = rxt(k,272)*y(k,73)
         mat(k,1935) = rxt(k,274)*y(k,73)
         mat(k,2146) = rxt(k,273)*y(k,73)
         mat(k,1951) = -(rxt(k,211)*y(k,221) + rxt(k,274)*y(k,73))
         mat(k,1725) = -rxt(k,211)*y(k,89)
         mat(k,679) = -rxt(k,274)*y(k,89)
         mat(k,1852) = rxt(k,295)*y(k,126)
         mat(k,1082) = rxt(k,328)*y(k,126)
         mat(k,1157) = rxt(k,354)*y(k,126)
         mat(k,912) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,85)
         mat(k,195) = rxt(k,502)*y(k,126)
         mat(k,1345) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,60)
         mat(k,1993) = rxt(k,210)*y(k,221)
         mat(k,1909) = rxt(k,295)*y(k,42) + rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) &
                      + rxt(k,502)*y(k,67)
         mat(k,1725) = mat(k,1725) + rxt(k,210)*y(k,125)
         mat(k,361) = -(rxt(k,187)*y(k,221))
         mat(k,1637) = -rxt(k,187)*y(k,90)
         mat(k,1961) = rxt(k,208)*y(k,203)
         mat(k,2044) = rxt(k,208)*y(k,125)
         mat(k,654) = -(rxt(k,265)*y(k,133) + (rxt(k,556) + rxt(k,561)) * y(k,85))
         mat(k,1803) = -rxt(k,265)*y(k,91)
         mat(k,1334) = -(rxt(k,556) + rxt(k,561)) * y(k,91)
         mat(k,2002) = rxt(k,257)*y(k,203)
         mat(k,2065) = rxt(k,257)*y(k,19)
         mat(k,762) = -(rxt(k,236)*y(k,56) + rxt(k,237)*y(k,133) + rxt(k,238)*y(k,221) &
                      + (rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,85))
         mat(k,1529) = -rxt(k,236)*y(k,92)
         mat(k,1809) = -rxt(k,237)*y(k,92)
         mat(k,1675) = -rxt(k,238)*y(k,92)
         mat(k,1335) = -(rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,92)
         mat(k,2122) = rxt(k,225)*y(k,203)
         mat(k,905) = rxt(k,230)*y(k,221)
         mat(k,2071) = rxt(k,225)*y(k,59)
         mat(k,1675) = mat(k,1675) + rxt(k,230)*y(k,60)
         mat(k,968) = -(rxt(k,371)*y(k,221))
         mat(k,1691) = -rxt(k,371)*y(k,93)
         mat(k,507) = .300_r8*rxt(k,416)*y(k,221)
         mat(k,425) = .500_r8*rxt(k,417)*y(k,221)
         mat(k,1485) = rxt(k,370)*y(k,200) + rxt(k,377)*y(k,207)
         mat(k,460) = rxt(k,370)*y(k,124)
         mat(k,1225) = rxt(k,377)*y(k,124)
         mat(k,1691) = mat(k,1691) + .300_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
                      *y(k,100)
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
         mat(k,155) = -(rxt(k,402)*y(k,221))
         mat(k,1606) = -rxt(k,402)*y(k,94)
         mat(k,991) = -(rxt(k,356)*y(k,221))
         mat(k,1693) = -rxt(k,356)*y(k,95)
         mat(k,508) = .700_r8*rxt(k,416)*y(k,221)
         mat(k,426) = .500_r8*rxt(k,417)*y(k,221)
         mat(k,467) = .500_r8*rxt(k,391)*y(k,221)
         mat(k,1487) = .050_r8*rxt(k,414)*y(k,206) + .220_r8*rxt(k,376)*y(k,207) &
                      + .250_r8*rxt(k,433)*y(k,229)
         mat(k,1879) = .050_r8*rxt(k,415)*y(k,206) + .220_r8*rxt(k,375)*y(k,207) &
                      + .250_r8*rxt(k,434)*y(k,229)
         mat(k,432) = .500_r8*rxt(k,360)*y(k,221)
         mat(k,1291) = .220_r8*rxt(k,372)*y(k,207) + .250_r8*rxt(k,430)*y(k,229)
         mat(k,1379) = .230_r8*rxt(k,373)*y(k,207) + .200_r8*rxt(k,361)*y(k,225) &
                      + .100_r8*rxt(k,431)*y(k,229)
         mat(k,1199) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1226) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,196) + .230_r8*rxt(k,373)*y(k,197)
         mat(k,1693) = mat(k,1693) + .700_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
                      *y(k,100) + .500_r8*rxt(k,391)*y(k,109) + .500_r8*rxt(k,360) &
                      *y(k,147)
         mat(k,1064) = .200_r8*rxt(k,361)*y(k,197)
         mat(k,1089) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,196) + .100_r8*rxt(k,431)*y(k,197)
         mat(k,212) = -(rxt(k,403)*y(k,221))
         mat(k,1616) = -rxt(k,403)*y(k,96)
         mat(k,1444) = .870_r8*rxt(k,414)*y(k,206)
         mat(k,1862) = .950_r8*rxt(k,415)*y(k,206)
         mat(k,1284) = rxt(k,410)*y(k,206)
         mat(k,1364) = .750_r8*rxt(k,411)*y(k,206)
         mat(k,1189) = .870_r8*rxt(k,414)*y(k,124) + .950_r8*rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,196) + .750_r8*rxt(k,411)*y(k,197)
         mat(k,68) = -(rxt(k,404)*y(k,221))
         mat(k,1594) = -rxt(k,404)*y(k,97)
         mat(k,577) = .600_r8*rxt(k,427)*y(k,221)
         mat(k,1594) = mat(k,1594) + .600_r8*rxt(k,427)*y(k,103)
         mat(k,781) = -(rxt(k,418)*y(k,126) + rxt(k,425)*y(k,135) + rxt(k,426) &
                      *y(k,221))
         mat(k,1866) = -rxt(k,418)*y(k,98)
         mat(k,1747) = -rxt(k,425)*y(k,98)
         mat(k,1677) = -rxt(k,426)*y(k,98)
         mat(k,505) = -(rxt(k,416)*y(k,221))
         mat(k,1654) = -rxt(k,416)*y(k,99)
         mat(k,1459) = .080_r8*rxt(k,408)*y(k,205)
         mat(k,1160) = .080_r8*rxt(k,408)*y(k,124)
         mat(k,422) = -(rxt(k,417)*y(k,221))
         mat(k,1644) = -rxt(k,417)*y(k,100)
         mat(k,1456) = .080_r8*rxt(k,414)*y(k,206)
         mat(k,1190) = .080_r8*rxt(k,414)*y(k,124)
         mat(k,1124) = -(rxt(k,419)*y(k,196) + rxt(k,420)*y(k,197) + rxt(k,421) &
                      *y(k,203) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126))
         mat(k,1294) = -rxt(k,419)*y(k,101)
         mat(k,1387) = -rxt(k,420)*y(k,101)
         mat(k,2090) = -rxt(k,421)*y(k,101)
         mat(k,1495) = -rxt(k,422)*y(k,101)
         mat(k,1887) = -rxt(k,423)*y(k,101)
         mat(k,784) = rxt(k,418)*y(k,126)
         mat(k,1887) = mat(k,1887) + rxt(k,418)*y(k,98)
         mat(k,297) = -(rxt(k,424)*y(k,221))
         mat(k,1628) = -rxt(k,424)*y(k,102)
         mat(k,1116) = rxt(k,421)*y(k,203)
         mat(k,2034) = rxt(k,421)*y(k,101)
         mat(k,578) = -(rxt(k,427)*y(k,221))
         mat(k,1660) = -rxt(k,427)*y(k,103)
         mat(k,2059) = rxt(k,407)*y(k,205) + rxt(k,412)*y(k,206)
         mat(k,1161) = rxt(k,407)*y(k,203)
         mat(k,1192) = rxt(k,412)*y(k,203)
         mat(k,39) = -(rxt(k,541)*y(k,221))
         mat(k,1588) = -rxt(k,541)*y(k,104)
         mat(k,1140) = -(rxt(k,378)*y(k,135) + rxt(k,379)*y(k,221))
         mat(k,1766) = -rxt(k,378)*y(k,105)
         mat(k,1703) = -rxt(k,379)*y(k,105)
         mat(k,785) = .300_r8*rxt(k,425)*y(k,135)
         mat(k,1496) = .360_r8*rxt(k,408)*y(k,205)
         mat(k,1888) = .400_r8*rxt(k,409)*y(k,205)
         mat(k,1766) = mat(k,1766) + .300_r8*rxt(k,425)*y(k,98)
         mat(k,1295) = .390_r8*rxt(k,405)*y(k,205)
         mat(k,1388) = .310_r8*rxt(k,406)*y(k,205)
         mat(k,1170) = .360_r8*rxt(k,408)*y(k,124) + .400_r8*rxt(k,409)*y(k,126) &
                      + .390_r8*rxt(k,405)*y(k,196) + .310_r8*rxt(k,406)*y(k,197)
         mat(k,215) = -(rxt(k,380)*y(k,221))
         mat(k,1617) = -rxt(k,380)*y(k,106)
         mat(k,2026) = rxt(k,374)*y(k,207)
         mat(k,1222) = rxt(k,374)*y(k,203)
         mat(k,417) = -(rxt(k,389)*y(k,221))
         mat(k,1643) = -rxt(k,389)*y(k,107)
         mat(k,1455) = .800_r8*rxt(k,398)*y(k,190)
         mat(k,875) = .800_r8*rxt(k,398)*y(k,124)
         mat(k,220) = -(rxt(k,390)*y(k,221))
         mat(k,1618) = -rxt(k,390)*y(k,108)
         mat(k,2027) = .800_r8*rxt(k,387)*y(k,211)
         mat(k,557) = .800_r8*rxt(k,387)*y(k,203)
         mat(k,466) = -(rxt(k,391)*y(k,221))
         mat(k,1649) = -rxt(k,391)*y(k,109)
         mat(k,1964) = rxt(k,394)*y(k,209)
         mat(k,1268) = rxt(k,394)*y(k,125)
         mat(k,820) = -(rxt(k,482)*y(k,126) + rxt(k,483)*y(k,135) + rxt(k,484) &
                      *y(k,221))
         mat(k,1868) = -rxt(k,482)*y(k,110)
         mat(k,1749) = -rxt(k,483)*y(k,110)
         mat(k,1680) = -rxt(k,484)*y(k,110)
         mat(k,1251) = -(rxt(k,392)*y(k,135) + rxt(k,393)*y(k,221))
         mat(k,1771) = -rxt(k,392)*y(k,111)
         mat(k,1708) = -rxt(k,393)*y(k,111)
         mat(k,788) = .200_r8*rxt(k,425)*y(k,135)
         mat(k,1501) = .560_r8*rxt(k,408)*y(k,205)
         mat(k,1893) = .600_r8*rxt(k,409)*y(k,205)
         mat(k,1771) = mat(k,1771) + .200_r8*rxt(k,425)*y(k,98)
         mat(k,1300) = .610_r8*rxt(k,405)*y(k,205)
         mat(k,1393) = .440_r8*rxt(k,406)*y(k,205)
         mat(k,1174) = .560_r8*rxt(k,408)*y(k,124) + .600_r8*rxt(k,409)*y(k,126) &
                      + .610_r8*rxt(k,405)*y(k,196) + .440_r8*rxt(k,406)*y(k,197)
         mat(k,744) = -(rxt(k,190)*y(k,124) + (rxt(k,191) + rxt(k,192) + rxt(k,193) &
                      ) * y(k,125) + rxt(k,194)*y(k,134) + rxt(k,202)*y(k,221) &
                      + rxt(k,574)*y(k,220))
         mat(k,1475) = -rxt(k,190)*y(k,112)
         mat(k,1969) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112)
         mat(k,1421) = -rxt(k,194)*y(k,112)
         mat(k,1673) = -rxt(k,202)*y(k,112)
         mat(k,684) = -rxt(k,574)*y(k,112)
         mat(k,1807) = rxt(k,188)*y(k,212) + rxt(k,571)*y(k,215)
         mat(k,1421) = mat(k,1421) + rxt(k,572)*y(k,215)
         mat(k,702) = 1.100_r8*rxt(k,567)*y(k,213) + .200_r8*rxt(k,565)*y(k,214)
         mat(k,413) = rxt(k,188)*y(k,133)
         mat(k,571) = 1.100_r8*rxt(k,567)*y(k,199)
         mat(k,692) = .200_r8*rxt(k,565)*y(k,199)
         mat(k,388) = rxt(k,571)*y(k,133) + rxt(k,572)*y(k,134)
         mat(k,1958) = rxt(k,209)*y(k,126)
         mat(k,1860) = rxt(k,209)*y(k,125)
         mat(k,267) = -(rxt(k,428)*y(k,221))
         mat(k,1624) = -rxt(k,428)*y(k,115)
         mat(k,1115) = .200_r8*rxt(k,420)*y(k,197)
         mat(k,1365) = .200_r8*rxt(k,420)*y(k,101)
         mat(k,928) = -(rxt(k,429)*y(k,221))
         mat(k,1688) = -rxt(k,429)*y(k,116)
         mat(k,1120) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) + rxt(k,419)*y(k,196) &
                      + .800_r8*rxt(k,420)*y(k,197)
         mat(k,1483) = rxt(k,422)*y(k,101)
         mat(k,1874) = rxt(k,423)*y(k,101)
         mat(k,1289) = rxt(k,419)*y(k,101)
         mat(k,1376) = .800_r8*rxt(k,420)*y(k,101)
         mat(k,49) = -(rxt(k,519)*y(k,221))
         mat(k,1590) = -rxt(k,519)*y(k,120)
         mat(k,1507) = -(rxt(k,190)*y(k,112) + rxt(k,199)*y(k,126) + rxt(k,203) &
                      *y(k,203) + rxt(k,204)*y(k,135) + rxt(k,205)*y(k,133) + rxt(k,226) &
                      *y(k,59) + rxt(k,258)*y(k,19) + rxt(k,301)*y(k,197) + rxt(k,310) &
                      *y(k,204) + rxt(k,323)*y(k,193) + rxt(k,334)*y(k,196) + rxt(k,338) &
                      *y(k,202) + rxt(k,351)*y(k,194) + rxt(k,359)*y(k,224) + rxt(k,363) &
                      *y(k,225) + (rxt(k,369) + rxt(k,370)) * y(k,200) + (rxt(k,376) &
                      + rxt(k,377)) * y(k,207) + rxt(k,385)*y(k,209) + rxt(k,388) &
                      *y(k,211) + (rxt(k,398) + rxt(k,399)) * y(k,190) + rxt(k,408) &
                      *y(k,205) + rxt(k,414)*y(k,206) + rxt(k,422)*y(k,101) + rxt(k,433) &
                      *y(k,229) + rxt(k,437)*y(k,189) + rxt(k,440)*y(k,191) + rxt(k,445) &
                      *y(k,192) + rxt(k,447)*y(k,195) + rxt(k,451)*y(k,198) + rxt(k,454) &
                      *y(k,208) + rxt(k,457)*y(k,210) + rxt(k,460)*y(k,223) + rxt(k,467) &
                      *y(k,228) + rxt(k,473)*y(k,230) + rxt(k,476)*y(k,231) + rxt(k,487) &
                      *y(k,216) + rxt(k,492)*y(k,226) + rxt(k,497)*y(k,227) + rxt(k,576) &
                      *y(k,220))
         mat(k,747) = -rxt(k,190)*y(k,124)
         mat(k,1900) = -rxt(k,199)*y(k,124)
         mat(k,2103) = -rxt(k,203)*y(k,124)
         mat(k,1778) = -rxt(k,204)*y(k,124)
         mat(k,1819) = -rxt(k,205)*y(k,124)
         mat(k,2130) = -rxt(k,226)*y(k,124)
         mat(k,2008) = -rxt(k,258)*y(k,124)
         mat(k,1398) = -rxt(k,301)*y(k,124)
         mat(k,331) = -rxt(k,310)*y(k,124)
         mat(k,844) = -rxt(k,323)*y(k,124)
         mat(k,1305) = -rxt(k,334)*y(k,124)
         mat(k,667) = -rxt(k,338)*y(k,124)
         mat(k,733) = -rxt(k,351)*y(k,124)
         mat(k,721) = -rxt(k,359)*y(k,124)
         mat(k,1069) = -rxt(k,363)*y(k,124)
         mat(k,462) = -(rxt(k,369) + rxt(k,370)) * y(k,124)
         mat(k,1235) = -(rxt(k,376) + rxt(k,377)) * y(k,124)
         mat(k,1274) = -rxt(k,385)*y(k,124)
         mat(k,561) = -rxt(k,388)*y(k,124)
         mat(k,886) = -(rxt(k,398) + rxt(k,399)) * y(k,124)
         mat(k,1178) = -rxt(k,408)*y(k,124)
         mat(k,1212) = -rxt(k,414)*y(k,124)
         mat(k,1130) = -rxt(k,422)*y(k,124)
         mat(k,1095) = -rxt(k,433)*y(k,124)
         mat(k,402) = -rxt(k,437)*y(k,124)
         mat(k,380) = -rxt(k,440)*y(k,124)
         mat(k,325) = -rxt(k,445)*y(k,124)
         mat(k,530) = -rxt(k,447)*y(k,124)
         mat(k,649) = -rxt(k,451)*y(k,124)
         mat(k,610) = -rxt(k,454)*y(k,124)
         mat(k,807) = -rxt(k,457)*y(k,124)
         mat(k,338) = -rxt(k,460)*y(k,124)
         mat(k,624) = -rxt(k,467)*y(k,124)
         mat(k,641) = -rxt(k,473)*y(k,124)
         mat(k,395) = -rxt(k,476)*y(k,124)
         mat(k,1053) = -rxt(k,487)*y(k,124)
         mat(k,1034) = -rxt(k,492)*y(k,124)
         mat(k,1014) = -rxt(k,497)*y(k,124)
         mat(k,686) = -rxt(k,576)*y(k,124)
         mat(k,747) = mat(k,747) + 2.000_r8*rxt(k,192)*y(k,125) + rxt(k,194)*y(k,134) &
                      + rxt(k,202)*y(k,221)
         mat(k,1984) = 2.000_r8*rxt(k,192)*y(k,112) + rxt(k,195)*y(k,133) + rxt(k,511) &
                      *y(k,151)
         mat(k,1819) = mat(k,1819) + rxt(k,195)*y(k,125)
         mat(k,1428) = rxt(k,194)*y(k,112) + rxt(k,189)*y(k,212)
         mat(k,1321) = rxt(k,511)*y(k,125)
         mat(k,415) = rxt(k,189)*y(k,134)
         mat(k,1716) = rxt(k,202)*y(k,112)
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
         mat(k,1994) = -((rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112) + (rxt(k,195) &
                      + rxt(k,197)) * y(k,133) + rxt(k,196)*y(k,135) + rxt(k,208) &
                      *y(k,203) + rxt(k,209)*y(k,126) + rxt(k,210)*y(k,221) + rxt(k,228) &
                      *y(k,59) + rxt(k,259)*y(k,19) + rxt(k,345)*y(k,196) + rxt(k,394) &
                      *y(k,209) + rxt(k,452)*y(k,198) + rxt(k,455)*y(k,208) + rxt(k,458) &
                      *y(k,210) + rxt(k,462)*y(k,142) + rxt(k,465)*y(k,189) + rxt(k,511) &
                      *y(k,151))
         mat(k,752) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,125)
         mat(k,1829) = -(rxt(k,195) + rxt(k,197)) * y(k,125)
         mat(k,1788) = -rxt(k,196)*y(k,125)
         mat(k,2113) = -rxt(k,208)*y(k,125)
         mat(k,1910) = -rxt(k,209)*y(k,125)
         mat(k,1726) = -rxt(k,210)*y(k,125)
         mat(k,2140) = -rxt(k,228)*y(k,125)
         mat(k,2018) = -rxt(k,259)*y(k,125)
         mat(k,1312) = -rxt(k,345)*y(k,125)
         mat(k,1281) = -rxt(k,394)*y(k,125)
         mat(k,651) = -rxt(k,452)*y(k,125)
         mat(k,611) = -rxt(k,455)*y(k,125)
         mat(k,809) = -rxt(k,458)*y(k,125)
         mat(k,360) = -rxt(k,462)*y(k,125)
         mat(k,404) = -rxt(k,465)*y(k,125)
         mat(k,1327) = -rxt(k,511)*y(k,125)
         mat(k,544) = rxt(k,396)*y(k,221)
         mat(k,266) = rxt(k,367)*y(k,126)
         mat(k,2018) = mat(k,2018) + rxt(k,258)*y(k,124)
         mat(k,2140) = mat(k,2140) + rxt(k,226)*y(k,124)
         mat(k,365) = rxt(k,187)*y(k,221)
         mat(k,512) = .700_r8*rxt(k,416)*y(k,221)
         mat(k,1136) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126)
         mat(k,1517) = rxt(k,258)*y(k,19) + rxt(k,226)*y(k,59) + rxt(k,422)*y(k,101) &
                      + 2.000_r8*rxt(k,199)*y(k,126) + rxt(k,205)*y(k,133) &
                      + rxt(k,204)*y(k,135) + rxt(k,437)*y(k,189) + rxt(k,398) &
                      *y(k,190) + rxt(k,440)*y(k,191) + rxt(k,445)*y(k,192) &
                      + rxt(k,323)*y(k,193) + rxt(k,351)*y(k,194) + rxt(k,447) &
                      *y(k,195) + rxt(k,334)*y(k,196) + rxt(k,301)*y(k,197) &
                      + rxt(k,451)*y(k,198) + rxt(k,369)*y(k,200) + rxt(k,338) &
                      *y(k,202) + rxt(k,203)*y(k,203) + rxt(k,310)*y(k,204) &
                      + .920_r8*rxt(k,408)*y(k,205) + .920_r8*rxt(k,414)*y(k,206) &
                      + rxt(k,376)*y(k,207) + rxt(k,454)*y(k,208) + rxt(k,385) &
                      *y(k,209) + rxt(k,457)*y(k,210) + rxt(k,388)*y(k,211) &
                      + 1.600_r8*rxt(k,487)*y(k,216) + rxt(k,460)*y(k,223) &
                      + rxt(k,359)*y(k,224) + rxt(k,363)*y(k,225) + .900_r8*rxt(k,492) &
                      *y(k,226) + .800_r8*rxt(k,497)*y(k,227) + rxt(k,467)*y(k,228) &
                      + rxt(k,433)*y(k,229) + rxt(k,473)*y(k,230) + rxt(k,476) &
                      *y(k,231)
         mat(k,1910) = mat(k,1910) + rxt(k,367)*y(k,16) + rxt(k,423)*y(k,101) &
                      + 2.000_r8*rxt(k,199)*y(k,124) + rxt(k,200)*y(k,133) &
                      + rxt(k,198)*y(k,203) + rxt(k,409)*y(k,205) + rxt(k,415) &
                      *y(k,206) + rxt(k,375)*y(k,207) + rxt(k,386)*y(k,209) &
                      + 2.000_r8*rxt(k,488)*y(k,216) + rxt(k,201)*y(k,221) &
                      + rxt(k,434)*y(k,229)
         mat(k,774) = rxt(k,357)*y(k,221)
         mat(k,1829) = mat(k,1829) + rxt(k,205)*y(k,124) + rxt(k,200)*y(k,126)
         mat(k,1788) = mat(k,1788) + rxt(k,204)*y(k,124)
         mat(k,525) = rxt(k,494)*y(k,221)
         mat(k,404) = mat(k,404) + rxt(k,437)*y(k,124)
         mat(k,889) = rxt(k,398)*y(k,124)
         mat(k,382) = rxt(k,440)*y(k,124)
         mat(k,327) = rxt(k,445)*y(k,124)
         mat(k,847) = rxt(k,323)*y(k,124)
         mat(k,736) = rxt(k,351)*y(k,124)
         mat(k,533) = rxt(k,447)*y(k,124)
         mat(k,1312) = mat(k,1312) + rxt(k,334)*y(k,124)
         mat(k,1406) = rxt(k,301)*y(k,124) + .500_r8*rxt(k,485)*y(k,216)
         mat(k,651) = mat(k,651) + rxt(k,451)*y(k,124)
         mat(k,464) = rxt(k,369)*y(k,124)
         mat(k,670) = rxt(k,338)*y(k,124)
         mat(k,2113) = mat(k,2113) + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126)
         mat(k,333) = rxt(k,310)*y(k,124)
         mat(k,1185) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126)
         mat(k,1219) = .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126)
         mat(k,1241) = rxt(k,376)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,611) = mat(k,611) + rxt(k,454)*y(k,124)
         mat(k,1281) = mat(k,1281) + rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126)
         mat(k,809) = mat(k,809) + rxt(k,457)*y(k,124)
         mat(k,563) = rxt(k,388)*y(k,124)
         mat(k,1059) = 1.600_r8*rxt(k,487)*y(k,124) + 2.000_r8*rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,197)
         mat(k,1726) = mat(k,1726) + rxt(k,396)*y(k,1) + rxt(k,187)*y(k,90) &
                      + .700_r8*rxt(k,416)*y(k,99) + rxt(k,201)*y(k,126) + rxt(k,357) &
                      *y(k,127) + rxt(k,494)*y(k,176)
         mat(k,340) = rxt(k,460)*y(k,124)
         mat(k,724) = rxt(k,359)*y(k,124)
         mat(k,1072) = rxt(k,363)*y(k,124)
         mat(k,1039) = .900_r8*rxt(k,492)*y(k,124)
         mat(k,1020) = .800_r8*rxt(k,497)*y(k,124)
         mat(k,626) = rxt(k,467)*y(k,124)
         mat(k,1100) = rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126)
         mat(k,643) = rxt(k,473)*y(k,124)
         mat(k,397) = rxt(k,476)*y(k,124)
         mat(k,1907) = -(rxt(k,198)*y(k,203) + rxt(k,199)*y(k,124) + rxt(k,200) &
                      *y(k,133) + rxt(k,201)*y(k,221) + rxt(k,209)*y(k,125) + rxt(k,295) &
                      *y(k,42) + rxt(k,328)*y(k,45) + rxt(k,347)*y(k,29) + rxt(k,354) &
                      *y(k,49) + rxt(k,367)*y(k,16) + rxt(k,375)*y(k,207) + rxt(k,386) &
                      *y(k,209) + rxt(k,409)*y(k,205) + rxt(k,415)*y(k,206) + rxt(k,418) &
                      *y(k,98) + rxt(k,423)*y(k,101) + rxt(k,434)*y(k,229) + rxt(k,479) &
                      *y(k,6) + rxt(k,482)*y(k,110) + rxt(k,488)*y(k,216) + rxt(k,499) &
                      *y(k,178) + rxt(k,502)*y(k,67))
         mat(k,2110) = -rxt(k,198)*y(k,126)
         mat(k,1514) = -rxt(k,199)*y(k,126)
         mat(k,1826) = -rxt(k,200)*y(k,126)
         mat(k,1723) = -rxt(k,201)*y(k,126)
         mat(k,1991) = -rxt(k,209)*y(k,126)
         mat(k,1850) = -rxt(k,295)*y(k,126)
         mat(k,1080) = -rxt(k,328)*y(k,126)
         mat(k,955) = -rxt(k,347)*y(k,126)
         mat(k,1155) = -rxt(k,354)*y(k,126)
         mat(k,265) = -rxt(k,367)*y(k,126)
         mat(k,1238) = -rxt(k,375)*y(k,126)
         mat(k,1278) = -rxt(k,386)*y(k,126)
         mat(k,1182) = -rxt(k,409)*y(k,126)
         mat(k,1216) = -rxt(k,415)*y(k,126)
         mat(k,794) = -rxt(k,418)*y(k,126)
         mat(k,1133) = -rxt(k,423)*y(k,126)
         mat(k,1098) = -rxt(k,434)*y(k,126)
         mat(k,872) = -rxt(k,479)*y(k,126)
         mat(k,834) = -rxt(k,482)*y(k,126)
         mat(k,1056) = -rxt(k,488)*y(k,126)
         mat(k,922) = -rxt(k,499)*y(k,126)
         mat(k,194) = -rxt(k,502)*y(k,126)
         mat(k,442) = rxt(k,260)*y(k,133)
         mat(k,1549) = rxt(k,227)*y(k,60)
         mat(k,911) = rxt(k,227)*y(k,56) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,221)
         mat(k,677) = rxt(k,274)*y(k,89)
         mat(k,1949) = rxt(k,274)*y(k,73) + rxt(k,211)*y(k,221)
         mat(k,471) = .500_r8*rxt(k,391)*y(k,221)
         mat(k,1991) = mat(k,1991) + rxt(k,197)*y(k,133) + rxt(k,196)*y(k,135)
         mat(k,1826) = mat(k,1826) + rxt(k,260)*y(k,20) + rxt(k,229)*y(k,60) &
                      + rxt(k,197)*y(k,125)
         mat(k,1785) = rxt(k,196)*y(k,125)
         mat(k,355) = rxt(k,343)*y(k,221)
         mat(k,1723) = mat(k,1723) + rxt(k,230)*y(k,60) + rxt(k,211)*y(k,89) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,343)*y(k,140)
         mat(k,769) = -(rxt(k,357)*y(k,221))
         mat(k,1676) = -rxt(k,357)*y(k,127)
         mat(k,942) = rxt(k,347)*y(k,126)
         mat(k,423) = .500_r8*rxt(k,417)*y(k,221)
         mat(k,299) = rxt(k,424)*y(k,221)
         mat(k,268) = rxt(k,428)*y(k,221)
         mat(k,925) = rxt(k,429)*y(k,221)
         mat(k,1865) = rxt(k,347)*y(k,29)
         mat(k,1676) = mat(k,1676) + .500_r8*rxt(k,417)*y(k,100) + rxt(k,424)*y(k,102) &
                      + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116)
         mat(k,285) = -(rxt(k,489)*y(k,221))
         mat(k,1626) = -rxt(k,489)*y(k,128)
         mat(k,2032) = rxt(k,486)*y(k,216)
         mat(k,1042) = rxt(k,486)*y(k,203)
         mat(k,1824) = -(rxt(k,167)*y(k,135) + 4._r8*rxt(k,168)*y(k,133) + rxt(k,169) &
                      *y(k,134) + rxt(k,170)*y(k,77) + rxt(k,171)*y(k,79) + rxt(k,176) &
                      *y(k,203) + rxt(k,182)*y(k,221) + (rxt(k,195) + rxt(k,197) &
                      ) * y(k,125) + rxt(k,200)*y(k,126) + rxt(k,205)*y(k,124) &
                      + rxt(k,229)*y(k,60) + rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85) &
                      + rxt(k,237)*y(k,92) + rxt(k,260)*y(k,20) + rxt(k,261)*y(k,19) &
                      + rxt(k,263)*y(k,81) + rxt(k,265)*y(k,91) + rxt(k,296)*y(k,42) &
                      + rxt(k,504)*y(k,138) + (rxt(k,569) + rxt(k,570)) * y(k,213) &
                      + rxt(k,571)*y(k,215))
         mat(k,1783) = -rxt(k,167)*y(k,133)
         mat(k,1433) = -rxt(k,169)*y(k,133)
         mat(k,1110) = -rxt(k,170)*y(k,133)
         mat(k,478) = -rxt(k,171)*y(k,133)
         mat(k,2108) = -rxt(k,176)*y(k,133)
         mat(k,1721) = -rxt(k,182)*y(k,133)
         mat(k,1989) = -(rxt(k,195) + rxt(k,197)) * y(k,133)
         mat(k,1905) = -rxt(k,200)*y(k,133)
         mat(k,1512) = -rxt(k,205)*y(k,133)
         mat(k,910) = -rxt(k,229)*y(k,133)
         mat(k,2135) = -rxt(k,231)*y(k,133)
         mat(k,1342) = -rxt(k,234)*y(k,133)
         mat(k,766) = -rxt(k,237)*y(k,133)
         mat(k,441) = -rxt(k,260)*y(k,133)
         mat(k,2013) = -rxt(k,261)*y(k,133)
         mat(k,757) = -rxt(k,263)*y(k,133)
         mat(k,659) = -rxt(k,265)*y(k,133)
         mat(k,1848) = -rxt(k,296)*y(k,133)
         mat(k,257) = -rxt(k,504)*y(k,133)
         mat(k,576) = -(rxt(k,569) + rxt(k,570)) * y(k,133)
         mat(k,390) = -rxt(k,571)*y(k,133)
         mat(k,1925) = rxt(k,174)*y(k,203)
         mat(k,750) = rxt(k,190)*y(k,124) + rxt(k,191)*y(k,125) + rxt(k,194)*y(k,134) &
                      + rxt(k,574)*y(k,220)
         mat(k,1512) = mat(k,1512) + rxt(k,190)*y(k,112)
         mat(k,1989) = mat(k,1989) + rxt(k,191)*y(k,112)
         mat(k,1433) = mat(k,1433) + rxt(k,194)*y(k,112) + rxt(k,506)*y(k,149) &
                      + rxt(k,512)*y(k,151) + rxt(k,573)*y(k,215) + (rxt(k,156) &
                       +rxt(k,157))*y(k,217) + rxt(k,579)*y(k,222)
         mat(k,605) = rxt(k,506)*y(k,134)
         mat(k,1325) = rxt(k,512)*y(k,134)
         mat(k,707) = rxt(k,565)*y(k,214) + 1.150_r8*rxt(k,566)*y(k,220)
         mat(k,2108) = mat(k,2108) + rxt(k,174)*y(k,76)
         mat(k,696) = rxt(k,565)*y(k,199)
         mat(k,390) = mat(k,390) + rxt(k,573)*y(k,134)
         mat(k,1573) = (rxt(k,156)+rxt(k,157))*y(k,134)
         mat(k,688) = rxt(k,574)*y(k,112) + 1.150_r8*rxt(k,566)*y(k,199)
         mat(k,1721) = mat(k,1721) + 2.000_r8*rxt(k,184)*y(k,221)
         mat(k,520) = rxt(k,579)*y(k,134)
         mat(k,1427) = -(rxt(k,156)*y(k,217) + rxt(k,161)*y(k,218) + rxt(k,169) &
                      *y(k,133) + rxt(k,175)*y(k,76) + rxt(k,189)*y(k,212) + rxt(k,194) &
                      *y(k,112) + rxt(k,340)*y(k,201) + rxt(k,506)*y(k,149) + rxt(k,512) &
                      *y(k,151) + rxt(k,568)*y(k,213) + (rxt(k,572) + rxt(k,573) &
                      ) * y(k,215) + rxt(k,579)*y(k,222))
         mat(k,1567) = -rxt(k,156)*y(k,134)
         mat(k,75) = -rxt(k,161)*y(k,134)
         mat(k,1818) = -rxt(k,169)*y(k,134)
         mat(k,1919) = -rxt(k,175)*y(k,134)
         mat(k,414) = -rxt(k,189)*y(k,134)
         mat(k,746) = -rxt(k,194)*y(k,134)
         mat(k,347) = -rxt(k,340)*y(k,134)
         mat(k,602) = -rxt(k,506)*y(k,134)
         mat(k,1320) = -rxt(k,512)*y(k,134)
         mat(k,573) = -rxt(k,568)*y(k,134)
         mat(k,389) = -(rxt(k,572) + rxt(k,573)) * y(k,134)
         mat(k,519) = -rxt(k,579)*y(k,134)
         mat(k,1353) = rxt(k,252)*y(k,135) + rxt(k,251)*y(k,203)
         mat(k,2007) = 2.000_r8*rxt(k,253)*y(k,19) + (rxt(k,255)+rxt(k,256))*y(k,59) &
                      + rxt(k,261)*y(k,133) + rxt(k,257)*y(k,203)
         mat(k,1541) = rxt(k,220)*y(k,135) + rxt(k,218)*y(k,203)
         mat(k,2129) = (rxt(k,255)+rxt(k,256))*y(k,19) + (2.000_r8*rxt(k,222) &
                       +2.000_r8*rxt(k,223))*y(k,59) + rxt(k,231)*y(k,133) &
                      + rxt(k,225)*y(k,203) + rxt(k,233)*y(k,221)
         mat(k,1919) = mat(k,1919) + rxt(k,178)*y(k,135) + rxt(k,172)*y(k,203)
         mat(k,362) = rxt(k,187)*y(k,221)
         mat(k,746) = mat(k,746) + rxt(k,193)*y(k,125)
         mat(k,1506) = rxt(k,204)*y(k,135) + rxt(k,576)*y(k,220)
         mat(k,1983) = rxt(k,193)*y(k,112) + rxt(k,195)*y(k,133) + rxt(k,196)*y(k,135)
         mat(k,1899) = rxt(k,200)*y(k,133) + rxt(k,198)*y(k,203)
         mat(k,1818) = mat(k,1818) + rxt(k,261)*y(k,19) + rxt(k,231)*y(k,59) &
                      + rxt(k,195)*y(k,125) + rxt(k,200)*y(k,126) &
                      + 2.000_r8*rxt(k,168)*y(k,133) + 2.000_r8*rxt(k,167)*y(k,135) &
                      + rxt(k,176)*y(k,203) + rxt(k,160)*y(k,218) + rxt(k,182) &
                      *y(k,221)
         mat(k,1427) = mat(k,1427) + 2.000_r8*rxt(k,161)*y(k,218)
         mat(k,1777) = rxt(k,252)*y(k,17) + rxt(k,220)*y(k,56) + rxt(k,178)*y(k,76) &
                      + rxt(k,204)*y(k,124) + rxt(k,196)*y(k,125) &
                      + 2.000_r8*rxt(k,167)*y(k,133) + rxt(k,507)*y(k,149) &
                      + rxt(k,513)*y(k,151) + 2.000_r8*rxt(k,177)*y(k,203) &
                      + 2.000_r8*rxt(k,158)*y(k,217) + rxt(k,183)*y(k,221)
         mat(k,602) = mat(k,602) + rxt(k,507)*y(k,135)
         mat(k,1320) = mat(k,1320) + rxt(k,513)*y(k,135)
         mat(k,843) = rxt(k,322)*y(k,203)
         mat(k,732) = rxt(k,350)*y(k,203)
         mat(k,1397) = rxt(k,300)*y(k,203)
         mat(k,2102) = rxt(k,251)*y(k,17) + rxt(k,257)*y(k,19) + rxt(k,218)*y(k,56) &
                      + rxt(k,225)*y(k,59) + rxt(k,172)*y(k,76) + rxt(k,198)*y(k,126) &
                      + rxt(k,176)*y(k,133) + 2.000_r8*rxt(k,177)*y(k,135) &
                      + rxt(k,322)*y(k,193) + rxt(k,350)*y(k,194) + rxt(k,300) &
                      *y(k,197) + 2.000_r8*rxt(k,186)*y(k,203) + rxt(k,181)*y(k,221) &
                      + rxt(k,358)*y(k,224)
         mat(k,1567) = mat(k,1567) + 2.000_r8*rxt(k,158)*y(k,135)
         mat(k,75) = mat(k,75) + rxt(k,160)*y(k,133) + 2.000_r8*rxt(k,161)*y(k,134)
         mat(k,685) = rxt(k,576)*y(k,124)
         mat(k,1715) = rxt(k,233)*y(k,59) + rxt(k,187)*y(k,90) + rxt(k,182)*y(k,133) &
                      + rxt(k,183)*y(k,135) + rxt(k,181)*y(k,203)
         mat(k,720) = rxt(k,358)*y(k,203)
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
         mat(k,1782) = -(rxt(k,158)*y(k,217) + rxt(k,167)*y(k,133) + rxt(k,177) &
                      *y(k,203) + rxt(k,178)*y(k,76) + rxt(k,183)*y(k,221) + rxt(k,196) &
                      *y(k,125) + rxt(k,204)*y(k,124) + rxt(k,220)*y(k,56) + rxt(k,252) &
                      *y(k,17) + rxt(k,319)*y(k,25) + rxt(k,348)*y(k,29) + rxt(k,378) &
                      *y(k,105) + rxt(k,392)*y(k,111) + rxt(k,425)*y(k,98) + rxt(k,463) &
                      *y(k,142) + rxt(k,480)*y(k,6) + rxt(k,483)*y(k,110) + rxt(k,507) &
                      *y(k,149) + rxt(k,513)*y(k,151))
         mat(k,1572) = -rxt(k,158)*y(k,135)
         mat(k,1823) = -rxt(k,167)*y(k,135)
         mat(k,2107) = -rxt(k,177)*y(k,135)
         mat(k,1924) = -rxt(k,178)*y(k,135)
         mat(k,1720) = -rxt(k,183)*y(k,135)
         mat(k,1988) = -rxt(k,196)*y(k,135)
         mat(k,1511) = -rxt(k,204)*y(k,135)
         mat(k,1546) = -rxt(k,220)*y(k,135)
         mat(k,1356) = -rxt(k,252)*y(k,135)
         mat(k,455) = -rxt(k,319)*y(k,135)
         mat(k,953) = -rxt(k,348)*y(k,135)
         mat(k,1146) = -rxt(k,378)*y(k,135)
         mat(k,1258) = -rxt(k,392)*y(k,135)
         mat(k,792) = -rxt(k,425)*y(k,135)
         mat(k,359) = -rxt(k,463)*y(k,135)
         mat(k,870) = -rxt(k,480)*y(k,135)
         mat(k,832) = -rxt(k,483)*y(k,135)
         mat(k,604) = -rxt(k,507)*y(k,135)
         mat(k,1324) = -rxt(k,513)*y(k,135)
         mat(k,1823) = mat(k,1823) + rxt(k,169)*y(k,134)
         mat(k,1432) = rxt(k,169)*y(k,133)
         mat(k,1307) = .150_r8*rxt(k,333)*y(k,203)
         mat(k,2107) = mat(k,2107) + .150_r8*rxt(k,333)*y(k,196) + .150_r8*rxt(k,383) &
                      *y(k,209)
         mat(k,1276) = .150_r8*rxt(k,383)*y(k,203)
         mat(k,230) = -(rxt(k,514)*y(k,151))
         mat(k,1315) = -rxt(k,514)*y(k,137)
         mat(k,2000) = rxt(k,254)*y(k,59)
         mat(k,2121) = rxt(k,254)*y(k,19) + 2.000_r8*rxt(k,224)*y(k,59)
         mat(k,251) = -(rxt(k,504)*y(k,133) + rxt(k,505)*y(k,221))
         mat(k,1795) = -rxt(k,504)*y(k,138)
         mat(k,1622) = -rxt(k,505)*y(k,138)
         mat(k,965) = rxt(k,371)*y(k,221)
         mat(k,1441) = .100_r8*rxt(k,492)*y(k,226)
         mat(k,1607) = rxt(k,371)*y(k,93)
         mat(k,1023) = .100_r8*rxt(k,492)*y(k,124)
         mat(k,350) = -(rxt(k,343)*y(k,221))
         mat(k,1635) = -rxt(k,343)*y(k,140)
         mat(k,1959) = rxt(k,345)*y(k,196)
         mat(k,1285) = rxt(k,345)*y(k,125)
         mat(k,1957) = rxt(k,465)*y(k,189)
         mat(k,399) = rxt(k,465)*y(k,125)
         mat(k,357) = -(rxt(k,462)*y(k,125) + rxt(k,463)*y(k,135))
         mat(k,1960) = -rxt(k,462)*y(k,142)
         mat(k,1740) = -rxt(k,463)*y(k,142)
         mat(k,123) = .070_r8*rxt(k,449)*y(k,221)
         mat(k,1451) = rxt(k,447)*y(k,195)
         mat(k,98) = .060_r8*rxt(k,461)*y(k,221)
         mat(k,148) = .070_r8*rxt(k,477)*y(k,221)
         mat(k,528) = rxt(k,447)*y(k,124)
         mat(k,1636) = .070_r8*rxt(k,449)*y(k,66) + .060_r8*rxt(k,461)*y(k,143) &
                      + .070_r8*rxt(k,477)*y(k,185)
         mat(k,96) = -(rxt(k,461)*y(k,221))
         mat(k,1597) = -rxt(k,461)*y(k,143)
         mat(k,88) = .530_r8*rxt(k,438)*y(k,221)
         mat(k,1597) = mat(k,1597) + .530_r8*rxt(k,438)*y(k,7)
         mat(k,235) = -(rxt(k,464)*y(k,221))
         mat(k,1619) = -rxt(k,464)*y(k,144)
         mat(k,2028) = rxt(k,459)*y(k,223)
         mat(k,335) = rxt(k,459)*y(k,203)
         mat(k,430) = -(rxt(k,360)*y(k,221))
         mat(k,1645) = -rxt(k,360)*y(k,147)
         mat(k,2050) = rxt(k,358)*y(k,224)
         mat(k,716) = rxt(k,358)*y(k,203)
         mat(k,291) = -(rxt(k,364)*y(k,221))
         mat(k,1627) = -rxt(k,364)*y(k,148)
         mat(k,2033) = .850_r8*rxt(k,362)*y(k,225)
         mat(k,1062) = .850_r8*rxt(k,362)*y(k,203)
         mat(k,600) = -(rxt(k,506)*y(k,134) + rxt(k,507)*y(k,135) + rxt(k,510) &
                      *y(k,221))
         mat(k,1417) = -rxt(k,506)*y(k,149)
         mat(k,1744) = -rxt(k,507)*y(k,149)
         mat(k,1662) = -rxt(k,510)*y(k,149)
         mat(k,1318) = -(rxt(k,508)*y(k,19) + rxt(k,509)*y(k,59) + rxt(k,511)*y(k,125) &
                      + rxt(k,512)*y(k,134) + rxt(k,513)*y(k,135) + rxt(k,514) &
                      *y(k,137) + rxt(k,515)*y(k,221))
         mat(k,2004) = -rxt(k,508)*y(k,151)
         mat(k,2125) = -rxt(k,509)*y(k,151)
         mat(k,1979) = -rxt(k,511)*y(k,151)
         mat(k,1425) = -rxt(k,512)*y(k,151)
         mat(k,1774) = -rxt(k,513)*y(k,151)
         mat(k,232) = -rxt(k,514)*y(k,151)
         mat(k,1711) = -rxt(k,515)*y(k,151)
         mat(k,1814) = rxt(k,504)*y(k,138)
         mat(k,1425) = mat(k,1425) + rxt(k,506)*y(k,149)
         mat(k,1774) = mat(k,1774) + rxt(k,507)*y(k,149)
         mat(k,255) = rxt(k,504)*y(k,133)
         mat(k,601) = rxt(k,506)*y(k,134) + rxt(k,507)*y(k,135) + rxt(k,510)*y(k,221)
         mat(k,1711) = mat(k,1711) + rxt(k,510)*y(k,149)
         mat(k,892) = -(rxt(k,517)*y(k,221))
         mat(k,1684) = -rxt(k,517)*y(k,152)
         mat(k,2003) = rxt(k,508)*y(k,151)
         mat(k,2123) = rxt(k,509)*y(k,151)
         mat(k,192) = rxt(k,502)*y(k,126) + (rxt(k,503)+.500_r8*rxt(k,516))*y(k,221)
         mat(k,1972) = rxt(k,511)*y(k,151)
         mat(k,1871) = rxt(k,502)*y(k,67)
         mat(k,1422) = rxt(k,512)*y(k,151)
         mat(k,1752) = rxt(k,513)*y(k,151)
         mat(k,231) = rxt(k,514)*y(k,151)
         mat(k,253) = rxt(k,505)*y(k,221)
         mat(k,1317) = rxt(k,508)*y(k,19) + rxt(k,509)*y(k,59) + rxt(k,511)*y(k,125) &
                      + rxt(k,512)*y(k,134) + rxt(k,513)*y(k,135) + rxt(k,514) &
                      *y(k,137) + rxt(k,515)*y(k,221)
         mat(k,1684) = mat(k,1684) + (rxt(k,503)+.500_r8*rxt(k,516))*y(k,67) &
                      + rxt(k,505)*y(k,138) + rxt(k,515)*y(k,151)
         mat(k,173) = -(rxt(k,518)*y(k,232))
         mat(k,2147) = -rxt(k,518)*y(k,153)
         mat(k,891) = rxt(k,517)*y(k,221)
         mat(k,1610) = rxt(k,517)*y(k,152)
         mat(k,849) = .2202005_r8*rxt(k,535)*y(k,135) + .2202005_r8*rxt(k,536) &
                      *y(k,221)
         mat(k,81) = .0023005_r8*rxt(k,537)*y(k,221)
         mat(k,775) = .0031005_r8*rxt(k,540)*y(k,221)
         mat(k,34) = .2381005_r8*rxt(k,541)*y(k,221)
         mat(k,811) = .0508005_r8*rxt(k,543)*y(k,135) + .0508005_r8*rxt(k,544) &
                      *y(k,221)
         mat(k,1731) = .2202005_r8*rxt(k,535)*y(k,6) + .0508005_r8*rxt(k,543)*y(k,110)
         mat(k,40) = .5931005_r8*rxt(k,545)*y(k,221)
         mat(k,109) = .1364005_r8*rxt(k,546)*y(k,221)
         mat(k,133) = .1677005_r8*rxt(k,547)*y(k,221)
         mat(k,1583) = .2202005_r8*rxt(k,536)*y(k,6) + .0023005_r8*rxt(k,537)*y(k,7) &
                      + .0031005_r8*rxt(k,540)*y(k,98) + .2381005_r8*rxt(k,541) &
                      *y(k,104) + .0508005_r8*rxt(k,544)*y(k,110) &
                      + .5931005_r8*rxt(k,545)*y(k,173) + .1364005_r8*rxt(k,546) &
                      *y(k,181) + .1677005_r8*rxt(k,547)*y(k,183)
         mat(k,850) = .2067005_r8*rxt(k,535)*y(k,135) + .2067005_r8*rxt(k,536) &
                      *y(k,221)
         mat(k,82) = .0008005_r8*rxt(k,537)*y(k,221)
         mat(k,776) = .0035005_r8*rxt(k,540)*y(k,221)
         mat(k,35) = .1308005_r8*rxt(k,541)*y(k,221)
         mat(k,812) = .1149005_r8*rxt(k,543)*y(k,135) + .1149005_r8*rxt(k,544) &
                      *y(k,221)
         mat(k,1732) = .2067005_r8*rxt(k,535)*y(k,6) + .1149005_r8*rxt(k,543)*y(k,110)
         mat(k,41) = .1534005_r8*rxt(k,545)*y(k,221)
         mat(k,110) = .0101005_r8*rxt(k,546)*y(k,221)
         mat(k,134) = .0174005_r8*rxt(k,547)*y(k,221)
         mat(k,1584) = .2067005_r8*rxt(k,536)*y(k,6) + .0008005_r8*rxt(k,537)*y(k,7) &
                      + .0035005_r8*rxt(k,540)*y(k,98) + .1308005_r8*rxt(k,541) &
                      *y(k,104) + .1149005_r8*rxt(k,544)*y(k,110) &
                      + .1534005_r8*rxt(k,545)*y(k,173) + .0101005_r8*rxt(k,546) &
                      *y(k,181) + .0174005_r8*rxt(k,547)*y(k,183)
         mat(k,851) = .0653005_r8*rxt(k,535)*y(k,135) + .0653005_r8*rxt(k,536) &
                      *y(k,221)
         mat(k,83) = .0843005_r8*rxt(k,537)*y(k,221)
         mat(k,777) = .0003005_r8*rxt(k,540)*y(k,221)
         mat(k,36) = .0348005_r8*rxt(k,541)*y(k,221)
         mat(k,813) = .0348005_r8*rxt(k,543)*y(k,135) + .0348005_r8*rxt(k,544) &
                      *y(k,221)
         mat(k,1733) = .0653005_r8*rxt(k,535)*y(k,6) + .0348005_r8*rxt(k,543)*y(k,110)
         mat(k,42) = .0459005_r8*rxt(k,545)*y(k,221)
         mat(k,111) = .0763005_r8*rxt(k,546)*y(k,221)
         mat(k,135) = .086_r8*rxt(k,547)*y(k,221)
         mat(k,1585) = .0653005_r8*rxt(k,536)*y(k,6) + .0843005_r8*rxt(k,537)*y(k,7) &
                      + .0003005_r8*rxt(k,540)*y(k,98) + .0348005_r8*rxt(k,541) &
                      *y(k,104) + .0348005_r8*rxt(k,544)*y(k,110) &
                      + .0459005_r8*rxt(k,545)*y(k,173) + .0763005_r8*rxt(k,546) &
                      *y(k,181) + .086_r8*rxt(k,547)*y(k,183)
         mat(k,852) = .1749305_r8*rxt(k,534)*y(k,126) + .1284005_r8*rxt(k,535) &
                      *y(k,135) + .1284005_r8*rxt(k,536)*y(k,221)
         mat(k,84) = .0443005_r8*rxt(k,537)*y(k,221)
         mat(k,778) = .0590245_r8*rxt(k,538)*y(k,126) + .0033005_r8*rxt(k,539) &
                      *y(k,135) + .0271005_r8*rxt(k,540)*y(k,221)
         mat(k,37) = .0076005_r8*rxt(k,541)*y(k,221)
         mat(k,814) = .1749305_r8*rxt(k,542)*y(k,126) + .0554005_r8*rxt(k,543) &
                      *y(k,135) + .0554005_r8*rxt(k,544)*y(k,221)
         mat(k,1858) = .1749305_r8*rxt(k,534)*y(k,6) + .0590245_r8*rxt(k,538)*y(k,98) &
                      + .1749305_r8*rxt(k,542)*y(k,110)
         mat(k,1734) = .1284005_r8*rxt(k,535)*y(k,6) + .0033005_r8*rxt(k,539)*y(k,98) &
                      + .0554005_r8*rxt(k,543)*y(k,110)
         mat(k,43) = .0085005_r8*rxt(k,545)*y(k,221)
         mat(k,112) = .2157005_r8*rxt(k,546)*y(k,221)
         mat(k,136) = .0512005_r8*rxt(k,547)*y(k,221)
         mat(k,1586) = .1284005_r8*rxt(k,536)*y(k,6) + .0443005_r8*rxt(k,537)*y(k,7) &
                      + .0271005_r8*rxt(k,540)*y(k,98) + .0076005_r8*rxt(k,541) &
                      *y(k,104) + .0554005_r8*rxt(k,544)*y(k,110) &
                      + .0085005_r8*rxt(k,545)*y(k,173) + .2157005_r8*rxt(k,546) &
                      *y(k,181) + .0512005_r8*rxt(k,547)*y(k,183)
         mat(k,853) = .5901905_r8*rxt(k,534)*y(k,126) + .114_r8*rxt(k,535)*y(k,135) &
                      + .114_r8*rxt(k,536)*y(k,221)
         mat(k,85) = .1621005_r8*rxt(k,537)*y(k,221)
         mat(k,779) = .0250245_r8*rxt(k,538)*y(k,126) + .0474005_r8*rxt(k,540) &
                      *y(k,221)
         mat(k,38) = .0113005_r8*rxt(k,541)*y(k,221)
         mat(k,815) = .5901905_r8*rxt(k,542)*y(k,126) + .1278005_r8*rxt(k,543) &
                      *y(k,135) + .1278005_r8*rxt(k,544)*y(k,221)
         mat(k,1859) = .5901905_r8*rxt(k,534)*y(k,6) + .0250245_r8*rxt(k,538)*y(k,98) &
                      + .5901905_r8*rxt(k,542)*y(k,110)
         mat(k,1735) = .114_r8*rxt(k,535)*y(k,6) + .1278005_r8*rxt(k,543)*y(k,110)
         mat(k,44) = .0128005_r8*rxt(k,545)*y(k,221)
         mat(k,113) = .0738005_r8*rxt(k,546)*y(k,221)
         mat(k,137) = .1598005_r8*rxt(k,547)*y(k,221)
         mat(k,1587) = .114_r8*rxt(k,536)*y(k,6) + .1621005_r8*rxt(k,537)*y(k,7) &
                      + .0474005_r8*rxt(k,540)*y(k,98) + .0113005_r8*rxt(k,541) &
                      *y(k,104) + .1278005_r8*rxt(k,544)*y(k,110) &
                      + .0128005_r8*rxt(k,545)*y(k,173) + .0738005_r8*rxt(k,546) &
                      *y(k,181) + .1598005_r8*rxt(k,547)*y(k,183)
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
         mat(k,45) = -(rxt(k,545)*y(k,221))
         mat(k,1589) = -rxt(k,545)*y(k,173)
         mat(k,116) = .100_r8*rxt(k,469)*y(k,221)
         mat(k,138) = .230_r8*rxt(k,471)*y(k,221)
         mat(k,1602) = .100_r8*rxt(k,469)*y(k,181) + .230_r8*rxt(k,471)*y(k,183)
         mat(k,481) = -(rxt(k,493)*y(k,221))
         mat(k,1651) = -rxt(k,493)*y(k,175)
         mat(k,2053) = rxt(k,491)*y(k,226)
         mat(k,1024) = rxt(k,491)*y(k,203)
         mat(k,521) = -(rxt(k,494)*y(k,221))
         mat(k,1655) = -rxt(k,494)*y(k,176)
         mat(k,1460) = .200_r8*rxt(k,487)*y(k,216) + .200_r8*rxt(k,497)*y(k,227)
         mat(k,1368) = .500_r8*rxt(k,485)*y(k,216)
         mat(k,1043) = .200_r8*rxt(k,487)*y(k,124) + .500_r8*rxt(k,485)*y(k,197)
         mat(k,1002) = .200_r8*rxt(k,497)*y(k,124)
         mat(k,368) = -(rxt(k,498)*y(k,221))
         mat(k,1638) = -rxt(k,498)*y(k,177)
         mat(k,2045) = rxt(k,496)*y(k,227)
         mat(k,1001) = rxt(k,496)*y(k,203)
         mat(k,916) = -(rxt(k,499)*y(k,126) + rxt(k,500)*y(k,221))
         mat(k,1873) = -rxt(k,499)*y(k,178)
         mat(k,1687) = -rxt(k,500)*y(k,178)
         mat(k,861) = .330_r8*rxt(k,480)*y(k,135)
         mat(k,823) = .330_r8*rxt(k,483)*y(k,135)
         mat(k,1482) = .800_r8*rxt(k,487)*y(k,216) + .800_r8*rxt(k,497)*y(k,227)
         mat(k,1873) = mat(k,1873) + rxt(k,488)*y(k,216)
         mat(k,1754) = .330_r8*rxt(k,480)*y(k,6) + .330_r8*rxt(k,483)*y(k,110)
         mat(k,522) = rxt(k,494)*y(k,221)
         mat(k,1375) = .500_r8*rxt(k,485)*y(k,216) + rxt(k,495)*y(k,227)
         mat(k,1045) = .800_r8*rxt(k,487)*y(k,124) + rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,197)
         mat(k,1687) = mat(k,1687) + rxt(k,494)*y(k,176)
         mat(k,1005) = .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,197)
         mat(k,982) = -(rxt(k,501)*y(k,221))
         mat(k,1692) = -rxt(k,501)*y(k,179)
         mat(k,862) = .300_r8*rxt(k,480)*y(k,135)
         mat(k,824) = .300_r8*rxt(k,483)*y(k,135)
         mat(k,1486) = .900_r8*rxt(k,492)*y(k,226)
         mat(k,1757) = .300_r8*rxt(k,480)*y(k,6) + .300_r8*rxt(k,483)*y(k,110)
         mat(k,1378) = rxt(k,490)*y(k,226)
         mat(k,1028) = .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,197)
         mat(k,492) = -(rxt(k,468)*y(k,221))
         mat(k,1652) = -rxt(k,468)*y(k,180)
         mat(k,2054) = rxt(k,466)*y(k,228)
         mat(k,615) = rxt(k,466)*y(k,203)
         mat(k,114) = -(rxt(k,469)*y(k,221))
         mat(k,1600) = -rxt(k,469)*y(k,181)
         mat(k,130) = -(rxt(k,435)*y(k,221))
         mat(k,1603) = -rxt(k,435)*y(k,182)
         mat(k,2024) = rxt(k,432)*y(k,229)
         mat(k,1085) = rxt(k,432)*y(k,203)
         mat(k,139) = -(rxt(k,471)*y(k,221))
         mat(k,1604) = -rxt(k,471)*y(k,183)
         mat(k,589) = -(rxt(k,474)*y(k,221))
         mat(k,1661) = -rxt(k,474)*y(k,184)
         mat(k,2060) = rxt(k,472)*y(k,230)
         mat(k,632) = rxt(k,472)*y(k,203)
         mat(k,147) = -(rxt(k,477)*y(k,221))
         mat(k,1605) = -rxt(k,477)*y(k,185)
         mat(k,140) = .150_r8*rxt(k,471)*y(k,221)
         mat(k,1605) = mat(k,1605) + .150_r8*rxt(k,471)*y(k,183)
         mat(k,315) = -(rxt(k,478)*y(k,221))
         mat(k,1631) = -rxt(k,478)*y(k,186)
         mat(k,2037) = rxt(k,475)*y(k,231)
         mat(k,391) = rxt(k,475)*y(k,203)
         mat(k,400) = -(rxt(k,436)*y(k,203) + rxt(k,437)*y(k,124) + rxt(k,465) &
                      *y(k,125))
         mat(k,2048) = -rxt(k,436)*y(k,189)
         mat(k,1454) = -rxt(k,437)*y(k,189)
         mat(k,1962) = -rxt(k,465)*y(k,189)
         mat(k,170) = rxt(k,442)*y(k,221)
         mat(k,1641) = rxt(k,442)*y(k,22)
         mat(k,880) = -(rxt(k,397)*y(k,203) + (rxt(k,398) + rxt(k,399)) * y(k,124))
         mat(k,2076) = -rxt(k,397)*y(k,190)
         mat(k,1480) = -(rxt(k,398) + rxt(k,399)) * y(k,190)
         mat(k,550) = rxt(k,400)*y(k,221)
         mat(k,164) = rxt(k,401)*y(k,221)
         mat(k,1683) = rxt(k,400)*y(k,2) + rxt(k,401)*y(k,15)
         mat(k,377) = -(rxt(k,439)*y(k,203) + rxt(k,440)*y(k,124))
         mat(k,2046) = -rxt(k,439)*y(k,191)
         mat(k,1452) = -rxt(k,440)*y(k,191)
         mat(k,89) = .350_r8*rxt(k,438)*y(k,221)
         mat(k,275) = rxt(k,441)*y(k,221)
         mat(k,1639) = .350_r8*rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8)
         mat(k,323) = -(rxt(k,443)*y(k,203) + rxt(k,445)*y(k,124))
         mat(k,2038) = -rxt(k,443)*y(k,192)
         mat(k,1446) = -rxt(k,445)*y(k,192)
         mat(k,242) = rxt(k,444)*y(k,221)
         mat(k,117) = .070_r8*rxt(k,469)*y(k,221)
         mat(k,141) = .060_r8*rxt(k,471)*y(k,221)
         mat(k,1632) = rxt(k,444)*y(k,23) + .070_r8*rxt(k,469)*y(k,181) &
                      + .060_r8*rxt(k,471)*y(k,183)
         mat(k,840) = -(4._r8*rxt(k,320)*y(k,193) + rxt(k,321)*y(k,197) + rxt(k,322) &
                      *y(k,203) + rxt(k,323)*y(k,124))
         mat(k,1373) = -rxt(k,321)*y(k,193)
         mat(k,2075) = -rxt(k,322)*y(k,193)
         mat(k,1479) = -rxt(k,323)*y(k,193)
         mat(k,247) = .500_r8*rxt(k,325)*y(k,221)
         mat(k,207) = rxt(k,326)*y(k,56) + rxt(k,327)*y(k,221)
         mat(k,1531) = rxt(k,326)*y(k,28)
         mat(k,1681) = .500_r8*rxt(k,325)*y(k,27) + rxt(k,327)*y(k,28)
         mat(k,728) = -(rxt(k,349)*y(k,197) + rxt(k,350)*y(k,203) + rxt(k,351) &
                      *y(k,124))
         mat(k,1370) = -rxt(k,349)*y(k,194)
         mat(k,2069) = -rxt(k,350)*y(k,194)
         mat(k,1474) = -rxt(k,351)*y(k,194)
         mat(k,304) = rxt(k,352)*y(k,221)
         mat(k,56) = rxt(k,353)*y(k,221)
         mat(k,1672) = rxt(k,352)*y(k,30) + rxt(k,353)*y(k,31)
         mat(k,529) = -(rxt(k,446)*y(k,203) + rxt(k,447)*y(k,124))
         mat(k,2056) = -rxt(k,446)*y(k,195)
         mat(k,1461) = -rxt(k,447)*y(k,195)
         mat(k,183) = rxt(k,448)*y(k,221)
         mat(k,1461) = mat(k,1461) + rxt(k,437)*y(k,189)
         mat(k,1743) = rxt(k,463)*y(k,142)
         mat(k,358) = rxt(k,463)*y(k,135)
         mat(k,401) = rxt(k,437)*y(k,124) + .400_r8*rxt(k,436)*y(k,203)
         mat(k,2056) = mat(k,2056) + .400_r8*rxt(k,436)*y(k,189)
         mat(k,1656) = rxt(k,448)*y(k,32)
         mat(k,1302) = -(4._r8*rxt(k,331)*y(k,196) + rxt(k,332)*y(k,197) + rxt(k,333) &
                      *y(k,203) + rxt(k,334)*y(k,124) + rxt(k,345)*y(k,125) + rxt(k,372) &
                      *y(k,207) + rxt(k,405)*y(k,205) + rxt(k,410)*y(k,206) + rxt(k,419) &
                      *y(k,101) + rxt(k,430)*y(k,229))
         mat(k,1395) = -rxt(k,332)*y(k,196)
         mat(k,2098) = -rxt(k,333)*y(k,196)
         mat(k,1503) = -rxt(k,334)*y(k,196)
         mat(k,1978) = -rxt(k,345)*y(k,196)
         mat(k,1233) = -rxt(k,372)*y(k,196)
         mat(k,1176) = -rxt(k,405)*y(k,196)
         mat(k,1210) = -rxt(k,410)*y(k,196)
         mat(k,1128) = -rxt(k,419)*y(k,196)
         mat(k,1093) = -rxt(k,430)*y(k,196)
         mat(k,868) = .060_r8*rxt(k,480)*y(k,135)
         mat(k,1077) = rxt(k,328)*y(k,126) + rxt(k,329)*y(k,221)
         mat(k,1153) = rxt(k,354)*y(k,126) + rxt(k,355)*y(k,221)
         mat(k,407) = .500_r8*rxt(k,336)*y(k,221)
         mat(k,789) = .080_r8*rxt(k,425)*y(k,135)
         mat(k,1144) = .100_r8*rxt(k,378)*y(k,135)
         mat(k,830) = .060_r8*rxt(k,483)*y(k,135)
         mat(k,1253) = .280_r8*rxt(k,392)*y(k,135)
         mat(k,1503) = mat(k,1503) + .530_r8*rxt(k,376)*y(k,207) + rxt(k,385)*y(k,209) &
                      + rxt(k,388)*y(k,211) + rxt(k,363)*y(k,225)
         mat(k,1895) = rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) + .530_r8*rxt(k,375) &
                      *y(k,207) + rxt(k,386)*y(k,209)
         mat(k,1773) = .060_r8*rxt(k,480)*y(k,6) + .080_r8*rxt(k,425)*y(k,98) &
                      + .100_r8*rxt(k,378)*y(k,105) + .060_r8*rxt(k,483)*y(k,110) &
                      + .280_r8*rxt(k,392)*y(k,111)
         mat(k,985) = .650_r8*rxt(k,501)*y(k,221)
         mat(k,1302) = mat(k,1302) + .530_r8*rxt(k,372)*y(k,207)
         mat(k,1395) = mat(k,1395) + .260_r8*rxt(k,373)*y(k,207) + rxt(k,382)*y(k,209) &
                      + .300_r8*rxt(k,361)*y(k,225)
         mat(k,2098) = mat(k,2098) + .450_r8*rxt(k,383)*y(k,209) + .200_r8*rxt(k,387) &
                      *y(k,211) + .150_r8*rxt(k,362)*y(k,225)
         mat(k,1233) = mat(k,1233) + .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375) &
                      *y(k,126) + .530_r8*rxt(k,372)*y(k,196) + .260_r8*rxt(k,373) &
                      *y(k,197)
         mat(k,1272) = rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126) + rxt(k,382)*y(k,197) &
                      + .450_r8*rxt(k,383)*y(k,203) + 4.000_r8*rxt(k,384)*y(k,209)
         mat(k,560) = rxt(k,388)*y(k,124) + .200_r8*rxt(k,387)*y(k,203)
         mat(k,1710) = rxt(k,329)*y(k,45) + rxt(k,355)*y(k,49) + .500_r8*rxt(k,336) &
                      *y(k,51) + .650_r8*rxt(k,501)*y(k,179)
         mat(k,1067) = rxt(k,363)*y(k,124) + .300_r8*rxt(k,361)*y(k,197) &
                      + .150_r8*rxt(k,362)*y(k,203)
         mat(k,1396) = -(rxt(k,221)*y(k,59) + (4._r8*rxt(k,298) + 4._r8*rxt(k,299) &
                      ) * y(k,197) + rxt(k,300)*y(k,203) + rxt(k,301)*y(k,124) &
                      + rxt(k,321)*y(k,193) + rxt(k,332)*y(k,196) + rxt(k,349) &
                      *y(k,194) + rxt(k,361)*y(k,225) + rxt(k,373)*y(k,207) + rxt(k,382) &
                      *y(k,209) + rxt(k,406)*y(k,205) + rxt(k,411)*y(k,206) + rxt(k,420) &
                      *y(k,101) + rxt(k,431)*y(k,229) + rxt(k,485)*y(k,216) + rxt(k,490) &
                      *y(k,226) + rxt(k,495)*y(k,227))
         mat(k,2128) = -rxt(k,221)*y(k,197)
         mat(k,2101) = -rxt(k,300)*y(k,197)
         mat(k,1505) = -rxt(k,301)*y(k,197)
         mat(k,842) = -rxt(k,321)*y(k,197)
         mat(k,1303) = -rxt(k,332)*y(k,197)
         mat(k,731) = -rxt(k,349)*y(k,197)
         mat(k,1068) = -rxt(k,361)*y(k,197)
         mat(k,1234) = -rxt(k,373)*y(k,197)
         mat(k,1273) = -rxt(k,382)*y(k,197)
         mat(k,1177) = -rxt(k,406)*y(k,197)
         mat(k,1211) = -rxt(k,411)*y(k,197)
         mat(k,1129) = -rxt(k,420)*y(k,197)
         mat(k,1094) = -rxt(k,431)*y(k,197)
         mat(k,1052) = -rxt(k,485)*y(k,197)
         mat(k,1033) = -rxt(k,490)*y(k,197)
         mat(k,1013) = -rxt(k,495)*y(k,197)
         mat(k,949) = .280_r8*rxt(k,348)*y(k,135)
         mat(k,447) = rxt(k,335)*y(k,221)
         mat(k,310) = .700_r8*rxt(k,303)*y(k,221)
         mat(k,790) = .050_r8*rxt(k,425)*y(k,135)
         mat(k,1129) = mat(k,1129) + rxt(k,419)*y(k,196)
         mat(k,1505) = mat(k,1505) + rxt(k,334)*y(k,196) + .830_r8*rxt(k,451)*y(k,198) &
                      + .170_r8*rxt(k,457)*y(k,210)
         mat(k,1776) = .280_r8*rxt(k,348)*y(k,29) + .050_r8*rxt(k,425)*y(k,98)
         mat(k,1303) = mat(k,1303) + rxt(k,419)*y(k,101) + rxt(k,334)*y(k,124) &
                      + 4.000_r8*rxt(k,331)*y(k,196) + .900_r8*rxt(k,332)*y(k,197) &
                      + .450_r8*rxt(k,333)*y(k,203) + rxt(k,405)*y(k,205) + rxt(k,410) &
                      *y(k,206) + rxt(k,372)*y(k,207) + rxt(k,381)*y(k,209) &
                      + rxt(k,430)*y(k,229)
         mat(k,1396) = mat(k,1396) + .900_r8*rxt(k,332)*y(k,196)
         mat(k,648) = .830_r8*rxt(k,451)*y(k,124) + .330_r8*rxt(k,450)*y(k,203)
         mat(k,2101) = mat(k,2101) + .450_r8*rxt(k,333)*y(k,196) + .330_r8*rxt(k,450) &
                      *y(k,198) + .070_r8*rxt(k,456)*y(k,210)
         mat(k,1177) = mat(k,1177) + rxt(k,405)*y(k,196)
         mat(k,1211) = mat(k,1211) + rxt(k,410)*y(k,196)
         mat(k,1234) = mat(k,1234) + rxt(k,372)*y(k,196)
         mat(k,1273) = mat(k,1273) + rxt(k,381)*y(k,196)
         mat(k,806) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,203)
         mat(k,1714) = rxt(k,335)*y(k,50) + .700_r8*rxt(k,303)*y(k,53)
         mat(k,1094) = mat(k,1094) + rxt(k,430)*y(k,196)
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
         mat(k,645) = -(rxt(k,450)*y(k,203) + rxt(k,451)*y(k,124) + rxt(k,452) &
                      *y(k,125))
         mat(k,2064) = -rxt(k,450)*y(k,198)
         mat(k,1467) = -rxt(k,451)*y(k,198)
         mat(k,1967) = -rxt(k,452)*y(k,198)
         mat(k,701) = -(rxt(k,565)*y(k,214) + rxt(k,566)*y(k,220) + rxt(k,567) &
                      *y(k,213))
         mat(k,691) = -rxt(k,565)*y(k,199)
         mat(k,683) = -rxt(k,566)*y(k,199)
         mat(k,570) = -rxt(k,567)*y(k,199)
         mat(k,458) = -((rxt(k,369) + rxt(k,370)) * y(k,124))
         mat(k,1457) = -(rxt(k,369) + rxt(k,370)) * y(k,200)
         mat(k,260) = rxt(k,368)*y(k,221)
         mat(k,1648) = rxt(k,368)*y(k,16)
         mat(k,345) = -(rxt(k,340)*y(k,134))
         mat(k,1412) = -rxt(k,340)*y(k,201)
         mat(k,1450) = .750_r8*rxt(k,338)*y(k,202)
         mat(k,663) = .750_r8*rxt(k,338)*y(k,124)
         mat(k,664) = -(rxt(k,337)*y(k,203) + rxt(k,338)*y(k,124))
         mat(k,2066) = -rxt(k,337)*y(k,202)
         mat(k,1468) = -rxt(k,338)*y(k,202)
         mat(k,451) = rxt(k,344)*y(k,221)
         mat(k,1667) = rxt(k,344)*y(k,25)
         mat(k,2115) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76) + rxt(k,176) &
                      *y(k,133) + rxt(k,177)*y(k,135) + rxt(k,181)*y(k,221) &
                      + 4._r8*rxt(k,186)*y(k,203) + rxt(k,198)*y(k,126) + rxt(k,203) &
                      *y(k,124) + rxt(k,208)*y(k,125) + (rxt(k,218) + rxt(k,219) &
                      ) * y(k,56) + rxt(k,225)*y(k,59) + rxt(k,251)*y(k,17) + rxt(k,257) &
                      *y(k,19) + rxt(k,294)*y(k,42) + rxt(k,300)*y(k,197) + rxt(k,308) &
                      *y(k,204) + rxt(k,322)*y(k,193) + rxt(k,333)*y(k,196) + rxt(k,337) &
                      *y(k,202) + rxt(k,350)*y(k,194) + rxt(k,358)*y(k,224) + rxt(k,362) &
                      *y(k,225) + rxt(k,374)*y(k,207) + rxt(k,383)*y(k,209) + rxt(k,387) &
                      *y(k,211) + rxt(k,397)*y(k,190) + rxt(k,407)*y(k,205) + rxt(k,412) &
                      *y(k,206) + rxt(k,421)*y(k,101) + rxt(k,432)*y(k,229) + rxt(k,436) &
                      *y(k,189) + rxt(k,439)*y(k,191) + rxt(k,443)*y(k,192) + rxt(k,446) &
                      *y(k,195) + rxt(k,450)*y(k,198) + rxt(k,453)*y(k,208) + rxt(k,456) &
                      *y(k,210) + rxt(k,459)*y(k,223) + rxt(k,466)*y(k,228) + rxt(k,472) &
                      *y(k,230) + rxt(k,475)*y(k,231) + rxt(k,486)*y(k,216) + rxt(k,491) &
                      *y(k,226) + rxt(k,496)*y(k,227))
         mat(k,1932) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,203)
         mat(k,1831) = -rxt(k,176)*y(k,203)
         mat(k,1790) = -rxt(k,177)*y(k,203)
         mat(k,1728) = -rxt(k,181)*y(k,203)
         mat(k,1912) = -rxt(k,198)*y(k,203)
         mat(k,1519) = -rxt(k,203)*y(k,203)
         mat(k,1996) = -rxt(k,208)*y(k,203)
         mat(k,1554) = -(rxt(k,218) + rxt(k,219)) * y(k,203)
         mat(k,2142) = -rxt(k,225)*y(k,203)
         mat(k,1361) = -rxt(k,251)*y(k,203)
         mat(k,2020) = -rxt(k,257)*y(k,203)
         mat(k,1855) = -rxt(k,294)*y(k,203)
         mat(k,1407) = -rxt(k,300)*y(k,203)
         mat(k,334) = -rxt(k,308)*y(k,203)
         mat(k,848) = -rxt(k,322)*y(k,203)
         mat(k,1313) = -rxt(k,333)*y(k,203)
         mat(k,671) = -rxt(k,337)*y(k,203)
         mat(k,737) = -rxt(k,350)*y(k,203)
         mat(k,725) = -rxt(k,358)*y(k,203)
         mat(k,1073) = -rxt(k,362)*y(k,203)
         mat(k,1242) = -rxt(k,374)*y(k,203)
         mat(k,1282) = -rxt(k,383)*y(k,203)
         mat(k,564) = -rxt(k,387)*y(k,203)
         mat(k,890) = -rxt(k,397)*y(k,203)
         mat(k,1186) = -rxt(k,407)*y(k,203)
         mat(k,1220) = -rxt(k,412)*y(k,203)
         mat(k,1137) = -rxt(k,421)*y(k,203)
         mat(k,1101) = -rxt(k,432)*y(k,203)
         mat(k,405) = -rxt(k,436)*y(k,203)
         mat(k,383) = -rxt(k,439)*y(k,203)
         mat(k,328) = -rxt(k,443)*y(k,203)
         mat(k,534) = -rxt(k,446)*y(k,203)
         mat(k,652) = -rxt(k,450)*y(k,203)
         mat(k,612) = -rxt(k,453)*y(k,203)
         mat(k,810) = -rxt(k,456)*y(k,203)
         mat(k,341) = -rxt(k,459)*y(k,203)
         mat(k,627) = -rxt(k,466)*y(k,203)
         mat(k,644) = -rxt(k,472)*y(k,203)
         mat(k,398) = -rxt(k,475)*y(k,203)
         mat(k,1060) = -rxt(k,486)*y(k,203)
         mat(k,1040) = -rxt(k,491)*y(k,203)
         mat(k,1021) = -rxt(k,496)*y(k,203)
         mat(k,873) = .570_r8*rxt(k,480)*y(k,135)
         mat(k,91) = .650_r8*rxt(k,438)*y(k,221)
         mat(k,1361) = mat(k,1361) + rxt(k,250)*y(k,42)
         mat(k,2020) = mat(k,2020) + rxt(k,262)*y(k,221)
         mat(k,205) = .350_r8*rxt(k,317)*y(k,221)
         mat(k,457) = .130_r8*rxt(k,319)*y(k,135)
         mat(k,180) = rxt(k,324)*y(k,221)
         mat(k,957) = .280_r8*rxt(k,348)*y(k,135)
         mat(k,1855) = mat(k,1855) + rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,133)
         mat(k,54) = rxt(k,330)*y(k,221)
         mat(k,711) = rxt(k,302)*y(k,221)
         mat(k,1554) = mat(k,1554) + rxt(k,214)*y(k,42) + rxt(k,217)*y(k,79)
         mat(k,2142) = mat(k,2142) + rxt(k,221)*y(k,197) + rxt(k,232)*y(k,221)
         mat(k,1000) = rxt(k,305)*y(k,221)
         mat(k,125) = .730_r8*rxt(k,449)*y(k,221)
         mat(k,196) = .500_r8*rxt(k,516)*y(k,221)
         mat(k,964) = rxt(k,341)*y(k,221)
         mat(k,801) = rxt(k,342)*y(k,221)
         mat(k,1932) = mat(k,1932) + rxt(k,175)*y(k,134)
         mat(k,479) = rxt(k,217)*y(k,56) + rxt(k,171)*y(k,133) + rxt(k,180)*y(k,221)
         mat(k,108) = rxt(k,306)*y(k,221)
         mat(k,714) = rxt(k,307)*y(k,221)
         mat(k,979) = rxt(k,371)*y(k,221)
         mat(k,996) = rxt(k,356)*y(k,221)
         mat(k,795) = .370_r8*rxt(k,425)*y(k,135)
         mat(k,513) = .300_r8*rxt(k,416)*y(k,221)
         mat(k,429) = rxt(k,417)*y(k,221)
         mat(k,1137) = mat(k,1137) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) &
                      + rxt(k,419)*y(k,196) + 1.200_r8*rxt(k,420)*y(k,197)
         mat(k,302) = rxt(k,424)*y(k,221)
         mat(k,1149) = .140_r8*rxt(k,378)*y(k,135)
         mat(k,219) = .200_r8*rxt(k,380)*y(k,221)
         mat(k,473) = .500_r8*rxt(k,391)*y(k,221)
         mat(k,835) = .570_r8*rxt(k,483)*y(k,135)
         mat(k,1264) = .280_r8*rxt(k,392)*y(k,135)
         mat(k,272) = rxt(k,428)*y(k,221)
         mat(k,938) = rxt(k,429)*y(k,221)
         mat(k,1519) = mat(k,1519) + rxt(k,422)*y(k,101) + rxt(k,398)*y(k,190) &
                      + rxt(k,440)*y(k,191) + rxt(k,445)*y(k,192) + rxt(k,323) &
                      *y(k,193) + rxt(k,351)*y(k,194) + rxt(k,301)*y(k,197) &
                      + .170_r8*rxt(k,451)*y(k,198) + rxt(k,369)*y(k,200) &
                      + .250_r8*rxt(k,338)*y(k,202) + rxt(k,310)*y(k,204) &
                      + .920_r8*rxt(k,408)*y(k,205) + .920_r8*rxt(k,414)*y(k,206) &
                      + .470_r8*rxt(k,376)*y(k,207) + .400_r8*rxt(k,454)*y(k,208) &
                      + .830_r8*rxt(k,457)*y(k,210) + rxt(k,460)*y(k,223) + rxt(k,359) &
                      *y(k,224) + .900_r8*rxt(k,492)*y(k,226) + .800_r8*rxt(k,497) &
                      *y(k,227) + rxt(k,467)*y(k,228) + rxt(k,433)*y(k,229) &
                      + rxt(k,473)*y(k,230) + rxt(k,476)*y(k,231)
         mat(k,1912) = mat(k,1912) + rxt(k,295)*y(k,42) + rxt(k,423)*y(k,101) &
                      + rxt(k,409)*y(k,205) + rxt(k,415)*y(k,206) + .470_r8*rxt(k,375) &
                      *y(k,207) + rxt(k,201)*y(k,221) + rxt(k,434)*y(k,229)
         mat(k,1831) = mat(k,1831) + rxt(k,296)*y(k,42) + rxt(k,171)*y(k,79)
         mat(k,1438) = rxt(k,175)*y(k,76) + rxt(k,340)*y(k,201)
         mat(k,1790) = mat(k,1790) + .570_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,319) &
                      *y(k,25) + .280_r8*rxt(k,348)*y(k,29) + .370_r8*rxt(k,425) &
                      *y(k,98) + .140_r8*rxt(k,378)*y(k,105) + .570_r8*rxt(k,483) &
                      *y(k,110) + .280_r8*rxt(k,392)*y(k,111) + rxt(k,183)*y(k,221)
         mat(k,100) = .800_r8*rxt(k,461)*y(k,221)
         mat(k,896) = rxt(k,517)*y(k,221)
         mat(k,989) = .200_r8*rxt(k,501)*y(k,221)
         mat(k,120) = .280_r8*rxt(k,469)*y(k,221)
         mat(k,146) = .380_r8*rxt(k,471)*y(k,221)
         mat(k,151) = .630_r8*rxt(k,477)*y(k,221)
         mat(k,890) = mat(k,890) + rxt(k,398)*y(k,124)
         mat(k,383) = mat(k,383) + rxt(k,440)*y(k,124)
         mat(k,328) = mat(k,328) + rxt(k,445)*y(k,124)
         mat(k,848) = mat(k,848) + rxt(k,323)*y(k,124) + 2.400_r8*rxt(k,320)*y(k,193) &
                      + rxt(k,321)*y(k,197)
         mat(k,737) = mat(k,737) + rxt(k,351)*y(k,124) + rxt(k,349)*y(k,197)
         mat(k,1313) = mat(k,1313) + rxt(k,419)*y(k,101) + .900_r8*rxt(k,332)*y(k,197) &
                      + rxt(k,405)*y(k,205) + rxt(k,410)*y(k,206) + .470_r8*rxt(k,372) &
                      *y(k,207) + rxt(k,430)*y(k,229)
         mat(k,1407) = mat(k,1407) + rxt(k,221)*y(k,59) + 1.200_r8*rxt(k,420)*y(k,101) &
                      + rxt(k,301)*y(k,124) + rxt(k,321)*y(k,193) + rxt(k,349) &
                      *y(k,194) + .900_r8*rxt(k,332)*y(k,196) + 4.000_r8*rxt(k,298) &
                      *y(k,197) + rxt(k,406)*y(k,205) + rxt(k,411)*y(k,206) &
                      + .730_r8*rxt(k,373)*y(k,207) + rxt(k,382)*y(k,209) &
                      + .500_r8*rxt(k,485)*y(k,216) + .300_r8*rxt(k,361)*y(k,225) &
                      + rxt(k,490)*y(k,226) + rxt(k,495)*y(k,227) + .800_r8*rxt(k,431) &
                      *y(k,229)
         mat(k,652) = mat(k,652) + .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450) &
                      *y(k,203)
         mat(k,465) = rxt(k,369)*y(k,124)
         mat(k,349) = rxt(k,340)*y(k,134)
         mat(k,671) = mat(k,671) + .250_r8*rxt(k,338)*y(k,124)
         mat(k,2115) = mat(k,2115) + .070_r8*rxt(k,450)*y(k,198) + .160_r8*rxt(k,453) &
                      *y(k,208) + .330_r8*rxt(k,456)*y(k,210)
         mat(k,334) = mat(k,334) + rxt(k,310)*y(k,124)
         mat(k,1186) = mat(k,1186) + .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) &
                      + rxt(k,405)*y(k,196) + rxt(k,406)*y(k,197)
         mat(k,1220) = mat(k,1220) + .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,196) + rxt(k,411)*y(k,197)
         mat(k,1242) = mat(k,1242) + .470_r8*rxt(k,376)*y(k,124) + .470_r8*rxt(k,375) &
                      *y(k,126) + .470_r8*rxt(k,372)*y(k,196) + .730_r8*rxt(k,373) &
                      *y(k,197)
         mat(k,612) = mat(k,612) + .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453) &
                      *y(k,203)
         mat(k,1282) = mat(k,1282) + rxt(k,382)*y(k,197)
         mat(k,810) = mat(k,810) + .830_r8*rxt(k,457)*y(k,124) + .330_r8*rxt(k,456) &
                      *y(k,203)
         mat(k,1060) = mat(k,1060) + .500_r8*rxt(k,485)*y(k,197)
         mat(k,1728) = mat(k,1728) + .650_r8*rxt(k,438)*y(k,7) + rxt(k,262)*y(k,19) &
                      + .350_r8*rxt(k,317)*y(k,24) + rxt(k,324)*y(k,26) + rxt(k,330) &
                      *y(k,47) + rxt(k,302)*y(k,52) + rxt(k,232)*y(k,59) + rxt(k,305) &
                      *y(k,62) + .730_r8*rxt(k,449)*y(k,66) + .500_r8*rxt(k,516) &
                      *y(k,67) + rxt(k,341)*y(k,74) + rxt(k,342)*y(k,75) + rxt(k,180) &
                      *y(k,79) + rxt(k,306)*y(k,86) + rxt(k,307)*y(k,87) + rxt(k,371) &
                      *y(k,93) + rxt(k,356)*y(k,95) + .300_r8*rxt(k,416)*y(k,99) &
                      + rxt(k,417)*y(k,100) + rxt(k,424)*y(k,102) + .200_r8*rxt(k,380) &
                      *y(k,106) + .500_r8*rxt(k,391)*y(k,109) + rxt(k,428)*y(k,115) &
                      + rxt(k,429)*y(k,116) + rxt(k,201)*y(k,126) + rxt(k,183) &
                      *y(k,135) + .800_r8*rxt(k,461)*y(k,143) + rxt(k,517)*y(k,152) &
                      + .200_r8*rxt(k,501)*y(k,179) + .280_r8*rxt(k,469)*y(k,181) &
                      + .380_r8*rxt(k,471)*y(k,183) + .630_r8*rxt(k,477)*y(k,185)
         mat(k,341) = mat(k,341) + rxt(k,460)*y(k,124)
         mat(k,725) = mat(k,725) + rxt(k,359)*y(k,124)
         mat(k,1073) = mat(k,1073) + .300_r8*rxt(k,361)*y(k,197)
         mat(k,1040) = mat(k,1040) + .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,197)
         mat(k,1021) = mat(k,1021) + .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,197)
         mat(k,627) = mat(k,627) + rxt(k,467)*y(k,124)
         mat(k,1101) = mat(k,1101) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126) &
                      + rxt(k,430)*y(k,196) + .800_r8*rxt(k,431)*y(k,197)
         mat(k,644) = mat(k,644) + rxt(k,473)*y(k,124)
         mat(k,398) = mat(k,398) + rxt(k,476)*y(k,124)
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
         mat(k,329) = -(rxt(k,308)*y(k,203) + rxt(k,310)*y(k,124))
         mat(k,2039) = -rxt(k,308)*y(k,204)
         mat(k,1447) = -rxt(k,310)*y(k,204)
         mat(k,1834) = rxt(k,294)*y(k,203)
         mat(k,2039) = mat(k,2039) + rxt(k,294)*y(k,42)
         mat(k,1172) = -(rxt(k,405)*y(k,196) + rxt(k,406)*y(k,197) + rxt(k,407) &
                      *y(k,203) + rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126))
         mat(k,1297) = -rxt(k,405)*y(k,205)
         mat(k,1390) = -rxt(k,406)*y(k,205)
         mat(k,2093) = -rxt(k,407)*y(k,205)
         mat(k,1498) = -rxt(k,408)*y(k,205)
         mat(k,1890) = -rxt(k,409)*y(k,205)
         mat(k,786) = .600_r8*rxt(k,426)*y(k,221)
         mat(k,1705) = .600_r8*rxt(k,426)*y(k,98)
         mat(k,1206) = -(rxt(k,410)*y(k,196) + rxt(k,411)*y(k,197) + rxt(k,412) &
                      *y(k,203) + rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126))
         mat(k,1298) = -rxt(k,410)*y(k,206)
         mat(k,1391) = -rxt(k,411)*y(k,206)
         mat(k,2094) = -rxt(k,412)*y(k,206)
         mat(k,1499) = -rxt(k,414)*y(k,206)
         mat(k,1891) = -rxt(k,415)*y(k,206)
         mat(k,787) = .400_r8*rxt(k,426)*y(k,221)
         mat(k,1706) = .400_r8*rxt(k,426)*y(k,98)
         mat(k,1231) = -(rxt(k,372)*y(k,196) + rxt(k,373)*y(k,197) + rxt(k,374) &
                      *y(k,203) + rxt(k,375)*y(k,126) + (rxt(k,376) + rxt(k,377) &
                      ) * y(k,124))
         mat(k,1299) = -rxt(k,372)*y(k,207)
         mat(k,1392) = -rxt(k,373)*y(k,207)
         mat(k,2095) = -rxt(k,374)*y(k,207)
         mat(k,1892) = -rxt(k,375)*y(k,207)
         mat(k,1500) = -(rxt(k,376) + rxt(k,377)) * y(k,207)
         mat(k,1142) = .500_r8*rxt(k,379)*y(k,221)
         mat(k,216) = .200_r8*rxt(k,380)*y(k,221)
         mat(k,1250) = rxt(k,393)*y(k,221)
         mat(k,1707) = .500_r8*rxt(k,379)*y(k,105) + .200_r8*rxt(k,380)*y(k,106) &
                      + rxt(k,393)*y(k,111)
         mat(k,607) = -(rxt(k,453)*y(k,203) + rxt(k,454)*y(k,124) + rxt(k,455) &
                      *y(k,125))
         mat(k,2061) = -rxt(k,453)*y(k,208)
         mat(k,1464) = -rxt(k,454)*y(k,208)
         mat(k,1966) = -rxt(k,455)*y(k,208)
         mat(k,1271) = -(rxt(k,381)*y(k,196) + rxt(k,382)*y(k,197) + rxt(k,383) &
                      *y(k,203) + 4._r8*rxt(k,384)*y(k,209) + rxt(k,385)*y(k,124) &
                      + rxt(k,386)*y(k,126) + rxt(k,394)*y(k,125))
         mat(k,1301) = -rxt(k,381)*y(k,209)
         mat(k,1394) = -rxt(k,382)*y(k,209)
         mat(k,2097) = -rxt(k,383)*y(k,209)
         mat(k,1502) = -rxt(k,385)*y(k,209)
         mat(k,1894) = -rxt(k,386)*y(k,209)
         mat(k,1977) = -rxt(k,394)*y(k,209)
         mat(k,1143) = .500_r8*rxt(k,379)*y(k,221)
         mat(k,217) = .500_r8*rxt(k,380)*y(k,221)
         mat(k,1709) = .500_r8*rxt(k,379)*y(k,105) + .500_r8*rxt(k,380)*y(k,106)
         mat(k,803) = -(rxt(k,456)*y(k,203) + rxt(k,457)*y(k,124) + rxt(k,458) &
                      *y(k,125))
         mat(k,2074) = -rxt(k,456)*y(k,210)
         mat(k,1478) = -rxt(k,457)*y(k,210)
         mat(k,1971) = -rxt(k,458)*y(k,210)
         mat(k,558) = -(rxt(k,387)*y(k,203) + rxt(k,388)*y(k,124))
         mat(k,2058) = -rxt(k,387)*y(k,211)
         mat(k,1463) = -rxt(k,388)*y(k,211)
         mat(k,418) = rxt(k,389)*y(k,221)
         mat(k,221) = rxt(k,390)*y(k,221)
         mat(k,1659) = rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108)
         mat(k,412) = -(rxt(k,188)*y(k,133) + rxt(k,189)*y(k,134))
         mat(k,1797) = -rxt(k,188)*y(k,212)
         mat(k,1414) = -rxt(k,189)*y(k,212)
         mat(k,1797) = mat(k,1797) + rxt(k,569)*y(k,213)
         mat(k,697) = .900_r8*rxt(k,567)*y(k,213) + .800_r8*rxt(k,565)*y(k,214)
         mat(k,565) = rxt(k,569)*y(k,133) + .900_r8*rxt(k,567)*y(k,199)
         mat(k,689) = .800_r8*rxt(k,565)*y(k,199)
         mat(k,567) = -(rxt(k,567)*y(k,199) + rxt(k,568)*y(k,134) + (rxt(k,569) &
                      + rxt(k,570)) * y(k,133))
         mat(k,698) = -rxt(k,567)*y(k,213)
         mat(k,1416) = -rxt(k,568)*y(k,213)
         mat(k,1801) = -(rxt(k,569) + rxt(k,570)) * y(k,213)
         mat(k,690) = -(rxt(k,565)*y(k,199))
         mat(k,700) = -rxt(k,565)*y(k,214)
         mat(k,742) = rxt(k,574)*y(k,220)
         mat(k,1470) = rxt(k,576)*y(k,220)
         mat(k,1805) = rxt(k,569)*y(k,213)
         mat(k,1419) = rxt(k,573)*y(k,215)
         mat(k,569) = rxt(k,569)*y(k,133)
         mat(k,387) = rxt(k,573)*y(k,134)
         mat(k,682) = rxt(k,574)*y(k,112) + rxt(k,576)*y(k,124)
         mat(k,384) = -(rxt(k,571)*y(k,133) + (rxt(k,572) + rxt(k,573)) * y(k,134))
         mat(k,1796) = -rxt(k,571)*y(k,215)
         mat(k,1413) = -(rxt(k,572) + rxt(k,573)) * y(k,215)
         mat(k,1049) = -(rxt(k,485)*y(k,197) + rxt(k,486)*y(k,203) + rxt(k,487) &
                      *y(k,124) + rxt(k,488)*y(k,126))
         mat(k,1383) = -rxt(k,485)*y(k,216)
         mat(k,2085) = -rxt(k,486)*y(k,216)
         mat(k,1491) = -rxt(k,487)*y(k,216)
         mat(k,1883) = -rxt(k,488)*y(k,216)
         mat(k,865) = rxt(k,479)*y(k,126)
         mat(k,827) = rxt(k,482)*y(k,126)
         mat(k,1883) = mat(k,1883) + rxt(k,479)*y(k,6) + rxt(k,482)*y(k,110) &
                      + .500_r8*rxt(k,499)*y(k,178)
         mat(k,287) = rxt(k,489)*y(k,221)
         mat(k,920) = .500_r8*rxt(k,499)*y(k,126)
         mat(k,1697) = rxt(k,489)*y(k,128)
         mat(k,1570) = -(rxt(k,153)*y(k,77) + rxt(k,154)*y(k,232) + (rxt(k,156) &
                      + rxt(k,157)) * y(k,134) + rxt(k,158)*y(k,135) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,85) + (rxt(k,269) + rxt(k,270)) * y(k,81) &
                      + rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65) + rxt(k,314)*y(k,86))
         mat(k,1108) = -rxt(k,153)*y(k,217)
         mat(k,2158) = -rxt(k,154)*y(k,217)
         mat(k,1430) = -(rxt(k,156) + rxt(k,157)) * y(k,217)
         mat(k,1780) = -rxt(k,158)*y(k,217)
         mat(k,1340) = -(rxt(k,246) + rxt(k,247)) * y(k,217)
         mat(k,755) = -(rxt(k,269) + rxt(k,270)) * y(k,217)
         mat(k,61) = -rxt(k,275)*y(k,217)
         mat(k,104) = -rxt(k,276)*y(k,217)
         mat(k,106) = -rxt(k,314)*y(k,217)
         mat(k,1430) = mat(k,1430) + rxt(k,189)*y(k,212)
         mat(k,706) = .850_r8*rxt(k,566)*y(k,220)
         mat(k,416) = rxt(k,189)*y(k,134)
         mat(k,687) = .850_r8*rxt(k,566)*y(k,199)
         mat(k,74) = -(rxt(k,160)*y(k,133) + rxt(k,161)*y(k,134))
         mat(k,1793) = -rxt(k,160)*y(k,218)
         mat(k,1410) = -rxt(k,161)*y(k,218)
         mat(k,1793) = mat(k,1793) + rxt(k,164)*y(k,219)
         mat(k,1410) = mat(k,1410) + rxt(k,165)*y(k,219)
         mat(k,1736) = rxt(k,166)*y(k,219)
         mat(k,76) = rxt(k,164)*y(k,133) + rxt(k,165)*y(k,134) + rxt(k,166)*y(k,135)
         mat(k,77) = -(rxt(k,164)*y(k,133) + rxt(k,165)*y(k,134) + rxt(k,166)*y(k,135))
         mat(k,1794) = -rxt(k,164)*y(k,219)
         mat(k,1411) = -rxt(k,165)*y(k,219)
         mat(k,1737) = -rxt(k,166)*y(k,219)
         mat(k,1411) = mat(k,1411) + rxt(k,156)*y(k,217)
         mat(k,1558) = rxt(k,156)*y(k,134)
         mat(k,681) = -(rxt(k,566)*y(k,199) + rxt(k,574)*y(k,112) + rxt(k,576) &
                      *y(k,124))
         mat(k,699) = -rxt(k,566)*y(k,220)
         mat(k,741) = -rxt(k,574)*y(k,220)
         mat(k,1469) = -rxt(k,576)*y(k,220)
         mat(k,1418) = rxt(k,568)*y(k,213) + rxt(k,572)*y(k,215) + rxt(k,579)*y(k,222)
         mat(k,568) = rxt(k,568)*y(k,134)
         mat(k,386) = rxt(k,572)*y(k,134)
         mat(k,515) = rxt(k,579)*y(k,134)
         mat(k,1719) = -(rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,181)*y(k,203) &
                      + rxt(k,182)*y(k,133) + rxt(k,183)*y(k,135) + (4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185)) * y(k,221) + rxt(k,187)*y(k,90) + rxt(k,201) &
                      *y(k,126) + rxt(k,202)*y(k,112) + rxt(k,210)*y(k,125) + rxt(k,211) &
                      *y(k,89) + rxt(k,230)*y(k,60) + (rxt(k,232) + rxt(k,233) &
                      ) * y(k,59) + rxt(k,235)*y(k,85) + rxt(k,238)*y(k,92) + rxt(k,262) &
                      *y(k,19) + rxt(k,264)*y(k,81) + rxt(k,297)*y(k,42) + rxt(k,302) &
                      *y(k,52) + rxt(k,303)*y(k,53) + (rxt(k,305) + rxt(k,315) &
                      ) * y(k,62) + rxt(k,306)*y(k,86) + rxt(k,307)*y(k,87) + rxt(k,317) &
                      *y(k,24) + rxt(k,324)*y(k,26) + rxt(k,325)*y(k,27) + rxt(k,327) &
                      *y(k,28) + rxt(k,329)*y(k,45) + rxt(k,330)*y(k,47) + rxt(k,335) &
                      *y(k,50) + rxt(k,336)*y(k,51) + rxt(k,341)*y(k,74) + rxt(k,342) &
                      *y(k,75) + rxt(k,343)*y(k,140) + rxt(k,344)*y(k,25) + rxt(k,352) &
                      *y(k,30) + rxt(k,353)*y(k,31) + rxt(k,355)*y(k,49) + rxt(k,356) &
                      *y(k,95) + rxt(k,357)*y(k,127) + rxt(k,360)*y(k,147) + rxt(k,364) &
                      *y(k,148) + rxt(k,365)*y(k,29) + rxt(k,366)*y(k,48) + rxt(k,368) &
                      *y(k,16) + rxt(k,371)*y(k,93) + rxt(k,379)*y(k,105) + rxt(k,380) &
                      *y(k,106) + rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108) + rxt(k,391) &
                      *y(k,109) + rxt(k,393)*y(k,111) + rxt(k,396)*y(k,1) + rxt(k,400) &
                      *y(k,2) + rxt(k,401)*y(k,15) + rxt(k,402)*y(k,94) + rxt(k,403) &
                      *y(k,96) + rxt(k,404)*y(k,97) + rxt(k,416)*y(k,99) + rxt(k,417) &
                      *y(k,100) + rxt(k,424)*y(k,102) + rxt(k,426)*y(k,98) + rxt(k,427) &
                      *y(k,103) + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116) + rxt(k,435) &
                      *y(k,182) + rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8) + rxt(k,442) &
                      *y(k,22) + rxt(k,444)*y(k,23) + rxt(k,448)*y(k,32) + rxt(k,449) &
                      *y(k,66) + rxt(k,461)*y(k,143) + rxt(k,464)*y(k,144) + rxt(k,468) &
                      *y(k,180) + rxt(k,469)*y(k,181) + rxt(k,471)*y(k,183) + rxt(k,474) &
                      *y(k,184) + rxt(k,477)*y(k,185) + rxt(k,478)*y(k,186) + rxt(k,481) &
                      *y(k,6) + rxt(k,484)*y(k,110) + rxt(k,489)*y(k,128) + rxt(k,493) &
                      *y(k,175) + rxt(k,494)*y(k,176) + rxt(k,498)*y(k,177) + rxt(k,500) &
                      *y(k,178) + rxt(k,501)*y(k,179) + (rxt(k,503) + rxt(k,516) &
                      ) * y(k,67) + rxt(k,505)*y(k,138) + rxt(k,510)*y(k,149) &
                      + rxt(k,515)*y(k,151) + rxt(k,517)*y(k,152) + rxt(k,519) &
                      *y(k,120))
         mat(k,1109) = -rxt(k,179)*y(k,221)
         mat(k,477) = -rxt(k,180)*y(k,221)
         mat(k,2106) = -rxt(k,181)*y(k,221)
         mat(k,1822) = -rxt(k,182)*y(k,221)
         mat(k,1781) = -rxt(k,183)*y(k,221)
         mat(k,363) = -rxt(k,187)*y(k,221)
         mat(k,1903) = -rxt(k,201)*y(k,221)
         mat(k,749) = -rxt(k,202)*y(k,221)
         mat(k,1987) = -rxt(k,210)*y(k,221)
         mat(k,1945) = -rxt(k,211)*y(k,221)
         mat(k,909) = -rxt(k,230)*y(k,221)
         mat(k,2133) = -(rxt(k,232) + rxt(k,233)) * y(k,221)
         mat(k,1341) = -rxt(k,235)*y(k,221)
         mat(k,765) = -rxt(k,238)*y(k,221)
         mat(k,2011) = -rxt(k,262)*y(k,221)
         mat(k,756) = -rxt(k,264)*y(k,221)
         mat(k,1846) = -rxt(k,297)*y(k,221)
         mat(k,709) = -rxt(k,302)*y(k,221)
         mat(k,311) = -rxt(k,303)*y(k,221)
         mat(k,998) = -(rxt(k,305) + rxt(k,315)) * y(k,221)
         mat(k,107) = -rxt(k,306)*y(k,221)
         mat(k,713) = -rxt(k,307)*y(k,221)
         mat(k,204) = -rxt(k,317)*y(k,221)
         mat(k,179) = -rxt(k,324)*y(k,221)
         mat(k,249) = -rxt(k,325)*y(k,221)
         mat(k,210) = -rxt(k,327)*y(k,221)
         mat(k,1079) = -rxt(k,329)*y(k,221)
         mat(k,53) = -rxt(k,330)*y(k,221)
         mat(k,448) = -rxt(k,335)*y(k,221)
         mat(k,409) = -rxt(k,336)*y(k,221)
         mat(k,962) = -rxt(k,341)*y(k,221)
         mat(k,800) = -rxt(k,342)*y(k,221)
         mat(k,353) = -rxt(k,343)*y(k,221)
         mat(k,454) = -rxt(k,344)*y(k,221)
         mat(k,306) = -rxt(k,352)*y(k,221)
         mat(k,57) = -rxt(k,353)*y(k,221)
         mat(k,1154) = -rxt(k,355)*y(k,221)
         mat(k,994) = -rxt(k,356)*y(k,221)
         mat(k,772) = -rxt(k,357)*y(k,221)
         mat(k,434) = -rxt(k,360)*y(k,221)
         mat(k,294) = -rxt(k,364)*y(k,221)
         mat(k,952) = -rxt(k,365)*y(k,221)
         mat(k,902) = -rxt(k,366)*y(k,221)
         mat(k,263) = -rxt(k,368)*y(k,221)
         mat(k,975) = -rxt(k,371)*y(k,221)
         mat(k,1145) = -rxt(k,379)*y(k,221)
         mat(k,218) = -rxt(k,380)*y(k,221)
         mat(k,421) = -rxt(k,389)*y(k,221)
         mat(k,224) = -rxt(k,390)*y(k,221)
         mat(k,469) = -rxt(k,391)*y(k,221)
         mat(k,1257) = -rxt(k,393)*y(k,221)
         mat(k,542) = -rxt(k,396)*y(k,221)
         mat(k,554) = -rxt(k,400)*y(k,221)
         mat(k,165) = -rxt(k,401)*y(k,221)
         mat(k,158) = -rxt(k,402)*y(k,221)
         mat(k,214) = -rxt(k,403)*y(k,221)
         mat(k,70) = -rxt(k,404)*y(k,221)
         mat(k,509) = -rxt(k,416)*y(k,221)
         mat(k,427) = -rxt(k,417)*y(k,221)
         mat(k,300) = -rxt(k,424)*y(k,221)
         mat(k,791) = -rxt(k,426)*y(k,221)
         mat(k,582) = -rxt(k,427)*y(k,221)
         mat(k,270) = -rxt(k,428)*y(k,221)
         mat(k,934) = -rxt(k,429)*y(k,221)
         mat(k,132) = -rxt(k,435)*y(k,221)
         mat(k,90) = -rxt(k,438)*y(k,221)
         mat(k,277) = -rxt(k,441)*y(k,221)
         mat(k,171) = -rxt(k,442)*y(k,221)
         mat(k,244) = -rxt(k,444)*y(k,221)
         mat(k,184) = -rxt(k,448)*y(k,221)
         mat(k,124) = -rxt(k,449)*y(k,221)
         mat(k,99) = -rxt(k,461)*y(k,221)
         mat(k,238) = -rxt(k,464)*y(k,221)
         mat(k,499) = -rxt(k,468)*y(k,221)
         mat(k,119) = -rxt(k,469)*y(k,221)
         mat(k,145) = -rxt(k,471)*y(k,221)
         mat(k,598) = -rxt(k,474)*y(k,221)
         mat(k,150) = -rxt(k,477)*y(k,221)
         mat(k,319) = -rxt(k,478)*y(k,221)
         mat(k,869) = -rxt(k,481)*y(k,221)
         mat(k,831) = -rxt(k,484)*y(k,221)
         mat(k,288) = -rxt(k,489)*y(k,221)
         mat(k,487) = -rxt(k,493)*y(k,221)
         mat(k,523) = -rxt(k,494)*y(k,221)
         mat(k,372) = -rxt(k,498)*y(k,221)
         mat(k,921) = -rxt(k,500)*y(k,221)
         mat(k,987) = -rxt(k,501)*y(k,221)
         mat(k,193) = -(rxt(k,503) + rxt(k,516)) * y(k,221)
         mat(k,256) = -rxt(k,505)*y(k,221)
         mat(k,603) = -rxt(k,510)*y(k,221)
         mat(k,1323) = -rxt(k,515)*y(k,221)
         mat(k,894) = -rxt(k,517)*y(k,221)
         mat(k,50) = -rxt(k,519)*y(k,221)
         mat(k,869) = mat(k,869) + .630_r8*rxt(k,480)*y(k,135)
         mat(k,204) = mat(k,204) + .650_r8*rxt(k,317)*y(k,221)
         mat(k,454) = mat(k,454) + .130_r8*rxt(k,319)*y(k,135)
         mat(k,249) = mat(k,249) + .500_r8*rxt(k,325)*y(k,221)
         mat(k,952) = mat(k,952) + .360_r8*rxt(k,348)*y(k,135)
         mat(k,1846) = mat(k,1846) + rxt(k,296)*y(k,133)
         mat(k,311) = mat(k,311) + .300_r8*rxt(k,303)*y(k,221)
         mat(k,1545) = rxt(k,219)*y(k,203)
         mat(k,676) = rxt(k,273)*y(k,232)
         mat(k,1923) = rxt(k,178)*y(k,135) + 2.000_r8*rxt(k,173)*y(k,203)
         mat(k,1109) = mat(k,1109) + rxt(k,170)*y(k,133) + rxt(k,153)*y(k,217)
         mat(k,477) = mat(k,477) + rxt(k,171)*y(k,133)
         mat(k,756) = mat(k,756) + rxt(k,263)*y(k,133) + rxt(k,269)*y(k,217)
         mat(k,1341) = mat(k,1341) + rxt(k,234)*y(k,133) + rxt(k,246)*y(k,217)
         mat(k,107) = mat(k,107) + rxt(k,314)*y(k,217)
         mat(k,658) = rxt(k,265)*y(k,133)
         mat(k,765) = mat(k,765) + rxt(k,237)*y(k,133)
         mat(k,791) = mat(k,791) + .320_r8*rxt(k,425)*y(k,135)
         mat(k,582) = mat(k,582) + .600_r8*rxt(k,427)*y(k,221)
         mat(k,1145) = mat(k,1145) + .240_r8*rxt(k,378)*y(k,135)
         mat(k,218) = mat(k,218) + .100_r8*rxt(k,380)*y(k,221)
         mat(k,831) = mat(k,831) + .630_r8*rxt(k,483)*y(k,135)
         mat(k,1257) = mat(k,1257) + .360_r8*rxt(k,392)*y(k,135)
         mat(k,1510) = rxt(k,203)*y(k,203)
         mat(k,1903) = mat(k,1903) + rxt(k,198)*y(k,203)
         mat(k,1822) = mat(k,1822) + rxt(k,296)*y(k,42) + rxt(k,170)*y(k,77) &
                      + rxt(k,171)*y(k,79) + rxt(k,263)*y(k,81) + rxt(k,234)*y(k,85) &
                      + rxt(k,265)*y(k,91) + rxt(k,237)*y(k,92) + rxt(k,176)*y(k,203)
         mat(k,1781) = mat(k,1781) + .630_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,319) &
                      *y(k,25) + .360_r8*rxt(k,348)*y(k,29) + rxt(k,178)*y(k,76) &
                      + .320_r8*rxt(k,425)*y(k,98) + .240_r8*rxt(k,378)*y(k,105) &
                      + .630_r8*rxt(k,483)*y(k,110) + .360_r8*rxt(k,392)*y(k,111) &
                      + rxt(k,177)*y(k,203)
         mat(k,434) = mat(k,434) + .500_r8*rxt(k,360)*y(k,221)
         mat(k,132) = mat(k,132) + .500_r8*rxt(k,435)*y(k,221)
         mat(k,403) = .400_r8*rxt(k,436)*y(k,203)
         mat(k,1306) = .450_r8*rxt(k,333)*y(k,203)
         mat(k,650) = .400_r8*rxt(k,450)*y(k,203)
         mat(k,2106) = mat(k,2106) + rxt(k,219)*y(k,56) + 2.000_r8*rxt(k,173)*y(k,76) &
                      + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126) + rxt(k,176) &
                      *y(k,133) + rxt(k,177)*y(k,135) + .400_r8*rxt(k,436)*y(k,189) &
                      + .450_r8*rxt(k,333)*y(k,196) + .400_r8*rxt(k,450)*y(k,198) &
                      + .450_r8*rxt(k,383)*y(k,209) + .400_r8*rxt(k,456)*y(k,210) &
                      + .200_r8*rxt(k,387)*y(k,211) + .150_r8*rxt(k,362)*y(k,225)
         mat(k,1275) = .450_r8*rxt(k,383)*y(k,203)
         mat(k,808) = .400_r8*rxt(k,456)*y(k,203)
         mat(k,562) = .200_r8*rxt(k,387)*y(k,203)
         mat(k,1571) = rxt(k,153)*y(k,77) + rxt(k,269)*y(k,81) + rxt(k,246)*y(k,85) &
                      + rxt(k,314)*y(k,86) + 2.000_r8*rxt(k,154)*y(k,232)
         mat(k,1719) = mat(k,1719) + .650_r8*rxt(k,317)*y(k,24) + .500_r8*rxt(k,325) &
                      *y(k,27) + .300_r8*rxt(k,303)*y(k,53) + .600_r8*rxt(k,427) &
                      *y(k,103) + .100_r8*rxt(k,380)*y(k,106) + .500_r8*rxt(k,360) &
                      *y(k,147) + .500_r8*rxt(k,435)*y(k,182)
         mat(k,1070) = .150_r8*rxt(k,362)*y(k,203)
         mat(k,2159) = rxt(k,273)*y(k,73) + 2.000_r8*rxt(k,154)*y(k,217)
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
         mat(k,514) = -(rxt(k,579)*y(k,134))
         mat(k,1415) = -rxt(k,579)*y(k,222)
         mat(k,1800) = rxt(k,570)*y(k,213) + rxt(k,571)*y(k,215)
         mat(k,566) = rxt(k,570)*y(k,133)
         mat(k,385) = rxt(k,571)*y(k,133)
         mat(k,336) = -(rxt(k,459)*y(k,203) + rxt(k,460)*y(k,124))
         mat(k,2040) = -rxt(k,459)*y(k,223)
         mat(k,1448) = -rxt(k,460)*y(k,223)
         mat(k,122) = .200_r8*rxt(k,449)*y(k,221)
         mat(k,97) = .140_r8*rxt(k,461)*y(k,221)
         mat(k,236) = rxt(k,464)*y(k,221)
         mat(k,1633) = .200_r8*rxt(k,449)*y(k,66) + .140_r8*rxt(k,461)*y(k,143) &
                      + rxt(k,464)*y(k,144)
         mat(k,717) = -(rxt(k,358)*y(k,203) + rxt(k,359)*y(k,124))
         mat(k,2068) = -rxt(k,358)*y(k,224)
         mat(k,1473) = -rxt(k,359)*y(k,224)
         mat(k,941) = rxt(k,365)*y(k,221)
         mat(k,431) = .500_r8*rxt(k,360)*y(k,221)
         mat(k,1671) = rxt(k,365)*y(k,29) + .500_r8*rxt(k,360)*y(k,147)
         mat(k,1065) = -(rxt(k,361)*y(k,197) + rxt(k,362)*y(k,203) + rxt(k,363) &
                      *y(k,124))
         mat(k,1384) = -rxt(k,361)*y(k,225)
         mat(k,2086) = -rxt(k,362)*y(k,225)
         mat(k,1492) = -rxt(k,363)*y(k,225)
         mat(k,866) = .060_r8*rxt(k,480)*y(k,135)
         mat(k,899) = rxt(k,366)*y(k,221)
         mat(k,828) = .060_r8*rxt(k,483)*y(k,135)
         mat(k,1763) = .060_r8*rxt(k,480)*y(k,6) + .060_r8*rxt(k,483)*y(k,110)
         mat(k,292) = rxt(k,364)*y(k,221)
         mat(k,984) = .150_r8*rxt(k,501)*y(k,221)
         mat(k,1698) = rxt(k,366)*y(k,48) + rxt(k,364)*y(k,148) + .150_r8*rxt(k,501) &
                      *y(k,179)
         mat(k,1030) = -(rxt(k,490)*y(k,197) + rxt(k,491)*y(k,203) + rxt(k,492) &
                      *y(k,124))
         mat(k,1382) = -rxt(k,490)*y(k,226)
         mat(k,2084) = -rxt(k,491)*y(k,226)
         mat(k,1490) = -rxt(k,492)*y(k,226)
         mat(k,1882) = .500_r8*rxt(k,499)*y(k,178)
         mat(k,486) = rxt(k,493)*y(k,221)
         mat(k,919) = .500_r8*rxt(k,499)*y(k,126) + rxt(k,500)*y(k,221)
         mat(k,1696) = rxt(k,493)*y(k,175) + rxt(k,500)*y(k,178)
         mat(k,1008) = -(rxt(k,495)*y(k,197) + rxt(k,496)*y(k,203) + rxt(k,497) &
                      *y(k,124))
         mat(k,1381) = -rxt(k,495)*y(k,227)
         mat(k,2083) = -rxt(k,496)*y(k,227)
         mat(k,1489) = -rxt(k,497)*y(k,227)
         mat(k,864) = rxt(k,481)*y(k,221)
         mat(k,826) = rxt(k,484)*y(k,221)
         mat(k,371) = rxt(k,498)*y(k,221)
         mat(k,1695) = rxt(k,481)*y(k,6) + rxt(k,484)*y(k,110) + rxt(k,498)*y(k,177)
         mat(k,618) = -(rxt(k,466)*y(k,203) + rxt(k,467)*y(k,124))
         mat(k,2062) = -rxt(k,466)*y(k,228)
         mat(k,1465) = -rxt(k,467)*y(k,228)
         mat(k,495) = rxt(k,468)*y(k,221)
         mat(k,118) = .650_r8*rxt(k,469)*y(k,221)
         mat(k,1664) = rxt(k,468)*y(k,180) + .650_r8*rxt(k,469)*y(k,181)
         mat(k,1091) = -(rxt(k,430)*y(k,196) + rxt(k,431)*y(k,197) + rxt(k,432) &
                      *y(k,203) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126))
         mat(k,1293) = -rxt(k,430)*y(k,229)
         mat(k,1386) = -rxt(k,431)*y(k,229)
         mat(k,2088) = -rxt(k,432)*y(k,229)
         mat(k,1494) = -rxt(k,433)*y(k,229)
         mat(k,1886) = -rxt(k,434)*y(k,229)
         mat(k,157) = rxt(k,402)*y(k,221)
         mat(k,213) = rxt(k,403)*y(k,221)
         mat(k,69) = rxt(k,404)*y(k,221)
         mat(k,579) = .400_r8*rxt(k,427)*y(k,221)
         mat(k,131) = .500_r8*rxt(k,435)*y(k,221)
         mat(k,1700) = rxt(k,402)*y(k,94) + rxt(k,403)*y(k,96) + rxt(k,404)*y(k,97) &
                      + .400_r8*rxt(k,427)*y(k,103) + .500_r8*rxt(k,435)*y(k,182)
         mat(k,634) = -(rxt(k,472)*y(k,203) + rxt(k,473)*y(k,124))
         mat(k,2063) = -rxt(k,472)*y(k,230)
         mat(k,1466) = -rxt(k,473)*y(k,230)
         mat(k,142) = .560_r8*rxt(k,471)*y(k,221)
         mat(k,591) = rxt(k,474)*y(k,221)
         mat(k,1665) = .560_r8*rxt(k,471)*y(k,183) + rxt(k,474)*y(k,184)
         mat(k,392) = -(rxt(k,475)*y(k,203) + rxt(k,476)*y(k,124))
         mat(k,2047) = -rxt(k,475)*y(k,231)
         mat(k,1453) = -rxt(k,476)*y(k,231)
         mat(k,149) = .300_r8*rxt(k,477)*y(k,221)
         mat(k,316) = rxt(k,478)*y(k,221)
         mat(k,1640) = .300_r8*rxt(k,477)*y(k,185) + rxt(k,478)*y(k,186)
         mat(k,2170) = -(rxt(k,154)*y(k,217) + rxt(k,273)*y(k,73) + rxt(k,518) &
                      *y(k,153))
         mat(k,1582) = -rxt(k,154)*y(k,232)
         mat(k,680) = -rxt(k,273)*y(k,232)
         mat(k,176) = -rxt(k,518)*y(k,232)
         mat(k,211) = rxt(k,327)*y(k,221)
         mat(k,308) = rxt(k,352)*y(k,221)
         mat(k,58) = rxt(k,353)*y(k,221)
         mat(k,1857) = rxt(k,297)*y(k,221)
         mat(k,1084) = rxt(k,329)*y(k,221)
         mat(k,903) = rxt(k,366)*y(k,221)
         mat(k,1159) = rxt(k,355)*y(k,221)
         mat(k,449) = rxt(k,335)*y(k,221)
         mat(k,411) = rxt(k,336)*y(k,221)
         mat(k,314) = rxt(k,303)*y(k,221)
         mat(k,1934) = rxt(k,174)*y(k,203)
         mat(k,1114) = rxt(k,179)*y(k,221)
         mat(k,480) = rxt(k,180)*y(k,221)
         mat(k,760) = rxt(k,264)*y(k,221)
         mat(k,1349) = (rxt(k,556)+rxt(k,561))*y(k,91) + (rxt(k,549)+rxt(k,555) &
                       +rxt(k,560))*y(k,92) + rxt(k,235)*y(k,221)
         mat(k,715) = rxt(k,307)*y(k,221)
         mat(k,1956) = rxt(k,211)*y(k,221)
         mat(k,367) = rxt(k,187)*y(k,221)
         mat(k,661) = (rxt(k,556)+rxt(k,561))*y(k,85)
         mat(k,768) = (rxt(k,549)+rxt(k,555)+rxt(k,560))*y(k,85) + rxt(k,238)*y(k,221)
         mat(k,1150) = .500_r8*rxt(k,379)*y(k,221)
         mat(k,51) = rxt(k,519)*y(k,221)
         mat(k,437) = rxt(k,360)*y(k,221)
         mat(k,296) = rxt(k,364)*y(k,221)
         mat(k,2117) = rxt(k,174)*y(k,76) + rxt(k,181)*y(k,221)
         mat(k,1730) = rxt(k,327)*y(k,28) + rxt(k,352)*y(k,30) + rxt(k,353)*y(k,31) &
                      + rxt(k,297)*y(k,42) + rxt(k,329)*y(k,45) + rxt(k,366)*y(k,48) &
                      + rxt(k,355)*y(k,49) + rxt(k,335)*y(k,50) + rxt(k,336)*y(k,51) &
                      + rxt(k,303)*y(k,53) + rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) &
                      + rxt(k,264)*y(k,81) + rxt(k,235)*y(k,85) + rxt(k,307)*y(k,87) &
                      + rxt(k,211)*y(k,89) + rxt(k,187)*y(k,90) + rxt(k,238)*y(k,92) &
                      + .500_r8*rxt(k,379)*y(k,105) + rxt(k,519)*y(k,120) + rxt(k,360) &
                      *y(k,147) + rxt(k,364)*y(k,148) + rxt(k,181)*y(k,203) &
                      + 2.000_r8*rxt(k,184)*y(k,221)
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
         mat(k, 29) = lmat(k, 29)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = lmat(k, 33)
         mat(k, 39) = mat(k, 39) + lmat(k, 39)
         mat(k, 45) = mat(k, 45) + lmat(k, 45)
         mat(k, 46) = lmat(k, 46)
         mat(k, 47) = lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 52) = mat(k, 52) + lmat(k, 52)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = lmat(k, 63)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 71) = lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = lmat(k, 78)
         mat(k, 79) = lmat(k, 79)
         mat(k, 80) = lmat(k, 80)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = lmat(k, 159)
         mat(k, 160) = lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = mat(k, 169) + lmat(k, 169)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = mat(k, 191) + lmat(k, 191)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 212) = mat(k, 212) + lmat(k, 212)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = lmat(k, 223)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 229) = lmat(k, 229)
         mat(k, 230) = mat(k, 230) + lmat(k, 230)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 237) = lmat(k, 237)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = lmat(k, 252)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 271) = lmat(k, 271)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = mat(k, 277) + lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 288) = mat(k, 288) + lmat(k, 288)
         mat(k, 289) = lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 297) = mat(k, 297) + lmat(k, 297)
         mat(k, 298) = lmat(k, 298)
         mat(k, 301) = lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 305) = lmat(k, 305)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 307) = lmat(k, 307)
         mat(k, 309) = mat(k, 309) + lmat(k, 309)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 313) = lmat(k, 313)
         mat(k, 315) = mat(k, 315) + lmat(k, 315)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 323) = mat(k, 323) + lmat(k, 323)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 332) = lmat(k, 332)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 348) = lmat(k, 348)
         mat(k, 349) = mat(k, 349) + lmat(k, 349)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = lmat(k, 352)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 356) = lmat(k, 356)
         mat(k, 357) = mat(k, 357) + lmat(k, 357)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 364) = lmat(k, 364)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = lmat(k, 370)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 373) = lmat(k, 373)
         mat(k, 374) = lmat(k, 374)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 400) = mat(k, 400) + lmat(k, 400)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 419) = lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 428) = lmat(k, 428)
         mat(k, 430) = mat(k, 430) + lmat(k, 430)
         mat(k, 433) = lmat(k, 433)
         mat(k, 434) = mat(k, 434) + lmat(k, 434)
         mat(k, 435) = lmat(k, 435)
         mat(k, 436) = lmat(k, 436)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 439) = lmat(k, 439)
         mat(k, 440) = lmat(k, 440)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 443) = lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 468) = lmat(k, 468)
         mat(k, 472) = lmat(k, 472)
         mat(k, 474) = mat(k, 474) + lmat(k, 474)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 483) = lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 493) = lmat(k, 493)
         mat(k, 497) = lmat(k, 497)
         mat(k, 498) = lmat(k, 498)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 500) = lmat(k, 500)
         mat(k, 501) = lmat(k, 501)
         mat(k, 502) = lmat(k, 502)
         mat(k, 503) = lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 511) = lmat(k, 511)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 517) = lmat(k, 517)
         mat(k, 518) = lmat(k, 518)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 524) = lmat(k, 524)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 535) = lmat(k, 535)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = lmat(k, 546)
         mat(k, 547) = mat(k, 547) + lmat(k, 547)
         mat(k, 551) = lmat(k, 551)
         mat(k, 552) = lmat(k, 552)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 555) = lmat(k, 555)
         mat(k, 556) = lmat(k, 556)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 580) = lmat(k, 580)
         mat(k, 581) = lmat(k, 581)
         mat(k, 583) = lmat(k, 583)
         mat(k, 584) = lmat(k, 584)
         mat(k, 585) = lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 593) = lmat(k, 593)
         mat(k, 596) = lmat(k, 596)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = mat(k, 600) + lmat(k, 600)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 634) = mat(k, 634) + lmat(k, 634)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 654) = mat(k, 654) + lmat(k, 654)
         mat(k, 656) = lmat(k, 656)
         mat(k, 658) = mat(k, 658) + lmat(k, 658)
         mat(k, 664) = mat(k, 664) + lmat(k, 664)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 675) = lmat(k, 675)
         mat(k, 681) = mat(k, 681) + lmat(k, 681)
         mat(k, 682) = mat(k, 682) + lmat(k, 682)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 701) = mat(k, 701) + lmat(k, 701)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 712) = mat(k, 712) + lmat(k, 712)
         mat(k, 717) = mat(k, 717) + lmat(k, 717)
         mat(k, 728) = mat(k, 728) + lmat(k, 728)
         mat(k, 739) = lmat(k, 739)
         mat(k, 743) = lmat(k, 743)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 753) = mat(k, 753) + lmat(k, 753)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 758) = mat(k, 758) + lmat(k, 758)
         mat(k, 762) = mat(k, 762) + lmat(k, 762)
         mat(k, 764) = mat(k, 764) + lmat(k, 764)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 769) = mat(k, 769) + lmat(k, 769)
         mat(k, 771) = lmat(k, 771)
         mat(k, 773) = lmat(k, 773)
         mat(k, 774) = mat(k, 774) + lmat(k, 774)
         mat(k, 781) = mat(k, 781) + lmat(k, 781)
         mat(k, 797) = lmat(k, 797)
         mat(k, 798) = mat(k, 798) + lmat(k, 798)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 801) = mat(k, 801) + lmat(k, 801)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 820) = mat(k, 820) + lmat(k, 820)
         mat(k, 840) = mat(k, 840) + lmat(k, 840)
         mat(k, 858) = mat(k, 858) + lmat(k, 858)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 892) = mat(k, 892) + lmat(k, 892)
         mat(k, 893) = lmat(k, 893)
         mat(k, 895) = lmat(k, 895)
         mat(k, 898) = mat(k, 898) + lmat(k, 898)
         mat(k, 900) = lmat(k, 900)
         mat(k, 901) = lmat(k, 901)
         mat(k, 905) = mat(k, 905) + lmat(k, 905)
         mat(k, 906) = mat(k, 906) + lmat(k, 906)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 913) = lmat(k, 913)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 917) = lmat(k, 917)
         mat(k, 918) = lmat(k, 918)
         mat(k, 923) = lmat(k, 923)
         mat(k, 924) = lmat(k, 924)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 936) = lmat(k, 936)
         mat(k, 937) = lmat(k, 937)
         mat(k, 938) = mat(k, 938) + lmat(k, 938)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 961) = lmat(k, 961)
         mat(k, 963) = mat(k, 963) + lmat(k, 963)
         mat(k, 964) = mat(k, 964) + lmat(k, 964)
         mat(k, 966) = lmat(k, 966)
         mat(k, 967) = lmat(k, 967)
         mat(k, 968) = mat(k, 968) + lmat(k, 968)
         mat(k, 969) = lmat(k, 969)
         mat(k, 970) = lmat(k, 970)
         mat(k, 972) = lmat(k, 972)
         mat(k, 973) = lmat(k, 973)
         mat(k, 976) = lmat(k, 976)
         mat(k, 977) = lmat(k, 977)
         mat(k, 978) = lmat(k, 978)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 981) = mat(k, 981) + lmat(k, 981)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 983) = mat(k, 983) + lmat(k, 983)
         mat(k, 984) = mat(k, 984) + lmat(k, 984)
         mat(k, 985) = mat(k, 985) + lmat(k, 985)
         mat(k, 988) = mat(k, 988) + lmat(k, 988)
         mat(k, 989) = mat(k, 989) + lmat(k, 989)
         mat(k, 991) = mat(k, 991) + lmat(k, 991)
         mat(k, 993) = lmat(k, 993)
         mat(k, 995) = lmat(k, 995)
         mat(k, 996) = mat(k, 996) + lmat(k, 996)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k,1008) = mat(k,1008) + lmat(k,1008)
         mat(k,1030) = mat(k,1030) + lmat(k,1030)
         mat(k,1049) = mat(k,1049) + lmat(k,1049)
         mat(k,1065) = mat(k,1065) + lmat(k,1065)
         mat(k,1075) = lmat(k,1075)
         mat(k,1076) = mat(k,1076) + lmat(k,1076)
         mat(k,1078) = lmat(k,1078)
         mat(k,1083) = lmat(k,1083)
         mat(k,1091) = mat(k,1091) + lmat(k,1091)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1139) = mat(k,1139) + lmat(k,1139)
         mat(k,1140) = mat(k,1140) + lmat(k,1140)
         mat(k,1143) = mat(k,1143) + lmat(k,1143)
         mat(k,1144) = mat(k,1144) + lmat(k,1144)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1151) = mat(k,1151) + lmat(k,1151)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1158) = lmat(k,1158)
         mat(k,1172) = mat(k,1172) + lmat(k,1172)
         mat(k,1188) = lmat(k,1188)
         mat(k,1206) = mat(k,1206) + lmat(k,1206)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1245) = lmat(k,1245)
         mat(k,1247) = mat(k,1247) + lmat(k,1247)
         mat(k,1251) = mat(k,1251) + lmat(k,1251)
         mat(k,1253) = mat(k,1253) + lmat(k,1253)
         mat(k,1254) = lmat(k,1254)
         mat(k,1271) = mat(k,1271) + lmat(k,1271)
         mat(k,1302) = mat(k,1302) + lmat(k,1302)
         mat(k,1316) = lmat(k,1316)
         mat(k,1318) = mat(k,1318) + lmat(k,1318)
         mat(k,1325) = mat(k,1325) + lmat(k,1325)
         mat(k,1337) = mat(k,1337) + lmat(k,1337)
         mat(k,1339) = mat(k,1339) + lmat(k,1339)
         mat(k,1344) = mat(k,1344) + lmat(k,1344)
         mat(k,1352) = mat(k,1352) + lmat(k,1352)
         mat(k,1396) = mat(k,1396) + lmat(k,1396)
         mat(k,1415) = mat(k,1415) + lmat(k,1415)
         mat(k,1418) = mat(k,1418) + lmat(k,1418)
         mat(k,1420) = lmat(k,1420)
         mat(k,1427) = mat(k,1427) + lmat(k,1427)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1470) = mat(k,1470) + lmat(k,1470)
         mat(k,1471) = lmat(k,1471)
         mat(k,1475) = mat(k,1475) + lmat(k,1475)
         mat(k,1507) = mat(k,1507) + lmat(k,1507)
         mat(k,1512) = mat(k,1512) + lmat(k,1512)
         mat(k,1534) = mat(k,1534) + lmat(k,1534)
         mat(k,1538) = mat(k,1538) + lmat(k,1538)
         mat(k,1539) = lmat(k,1539)
         mat(k,1540) = lmat(k,1540)
         mat(k,1543) = mat(k,1543) + lmat(k,1543)
         mat(k,1554) = mat(k,1554) + lmat(k,1554)
         mat(k,1557) = mat(k,1557) + lmat(k,1557)
         mat(k,1559) = mat(k,1559) + lmat(k,1559)
         mat(k,1561) = mat(k,1561) + lmat(k,1561)
         mat(k,1563) = mat(k,1563) + lmat(k,1563)
         mat(k,1565) = mat(k,1565) + lmat(k,1565)
         mat(k,1566) = lmat(k,1566)
         mat(k,1567) = mat(k,1567) + lmat(k,1567)
         mat(k,1568) = lmat(k,1568)
         mat(k,1569) = mat(k,1569) + lmat(k,1569)
         mat(k,1570) = mat(k,1570) + lmat(k,1570)
         mat(k,1571) = mat(k,1571) + lmat(k,1571)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1574) = lmat(k,1574)
         mat(k,1576) = mat(k,1576) + lmat(k,1576)
         mat(k,1580) = lmat(k,1580)
         mat(k,1593) = lmat(k,1593)
         mat(k,1598) = lmat(k,1598)
         mat(k,1713) = mat(k,1713) + lmat(k,1713)
         mat(k,1714) = mat(k,1714) + lmat(k,1714)
         mat(k,1717) = mat(k,1717) + lmat(k,1717)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1730) = mat(k,1730) + lmat(k,1730)
         mat(k,1736) = mat(k,1736) + lmat(k,1736)
         mat(k,1777) = mat(k,1777) + lmat(k,1777)
         mat(k,1780) = mat(k,1780) + lmat(k,1780)
         mat(k,1782) = mat(k,1782) + lmat(k,1782)
         mat(k,1783) = mat(k,1783) + lmat(k,1783)
         mat(k,1800) = mat(k,1800) + lmat(k,1800)
         mat(k,1806) = lmat(k,1806)
         mat(k,1824) = mat(k,1824) + lmat(k,1824)
         mat(k,1837) = mat(k,1837) + lmat(k,1837)
         mat(k,1838) = lmat(k,1838)
         mat(k,1849) = mat(k,1849) + lmat(k,1849)
         mat(k,1851) = mat(k,1851) + lmat(k,1851)
         mat(k,1899) = mat(k,1899) + lmat(k,1899)
         mat(k,1900) = mat(k,1900) + lmat(k,1900)
         mat(k,1905) = mat(k,1905) + lmat(k,1905)
         mat(k,1907) = mat(k,1907) + lmat(k,1907)
         mat(k,1909) = mat(k,1909) + lmat(k,1909)
         mat(k,1910) = mat(k,1910) + lmat(k,1910)
         mat(k,1928) = mat(k,1928) + lmat(k,1928)
         mat(k,1945) = mat(k,1945) + lmat(k,1945)
         mat(k,1951) = mat(k,1951) + lmat(k,1951)
         mat(k,1952) = lmat(k,1952)
         mat(k,1984) = mat(k,1984) + lmat(k,1984)
         mat(k,1987) = mat(k,1987) + lmat(k,1987)
         mat(k,1989) = mat(k,1989) + lmat(k,1989)
         mat(k,1993) = mat(k,1993) + lmat(k,1993)
         mat(k,1994) = mat(k,1994) + lmat(k,1994)
         mat(k,2006) = mat(k,2006) + lmat(k,2006)
         mat(k,2013) = mat(k,2013) + lmat(k,2013)
         mat(k,2019) = mat(k,2019) + lmat(k,2019)
         mat(k,2052) = mat(k,2052) + lmat(k,2052)
         mat(k,2115) = mat(k,2115) + lmat(k,2115)
         mat(k,2131) = mat(k,2131) + lmat(k,2131)
         mat(k,2135) = mat(k,2135) + lmat(k,2135)
         mat(k,2143) = mat(k,2143) + lmat(k,2143)
         mat(k,2150) = lmat(k,2150)
         mat(k,2158) = mat(k,2158) + lmat(k,2158)
         mat(k,2159) = mat(k,2159) + lmat(k,2159)
         mat(k,2161) = lmat(k,2161)
         mat(k,2164) = lmat(k,2164)
         mat(k,2170) = mat(k,2170) + lmat(k,2170)
         mat(k, 143) = 0._r8
         mat(k, 144) = 0._r8
         mat(k, 243) = 0._r8
         mat(k, 324) = 0._r8
         mat(k, 326) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 381) = 0._r8
         mat(k, 396) = 0._r8
         mat(k, 494) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 532) = 0._r8
         mat(k, 537) = 0._r8
         mat(k, 538) = 0._r8
         mat(k, 541) = 0._r8
         mat(k, 548) = 0._r8
         mat(k, 549) = 0._r8
         mat(k, 553) = 0._r8
         mat(k, 572) = 0._r8
         mat(k, 574) = 0._r8
         mat(k, 575) = 0._r8
         mat(k, 590) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 594) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 617) = 0._r8
         mat(k, 619) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 622) = 0._r8
         mat(k, 625) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 694) = 0._r8
         mat(k, 695) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 704) = 0._r8
         mat(k, 705) = 0._r8
         mat(k, 718) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 726) = 0._r8
         mat(k, 734) = 0._r8
         mat(k, 738) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 745) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 819) = 0._r8
         mat(k, 821) = 0._r8
         mat(k, 829) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 845) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 874) = 0._r8
         mat(k, 878) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 884) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 930) = 0._r8
         mat(k, 931) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 944) = 0._r8
         mat(k, 947) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 956) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 974) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 986) = 0._r8
         mat(k, 990) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1011) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1019) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1032) = 0._r8
         mat(k,1035) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1054) = 0._r8
         mat(k,1057) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1061) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1081) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1106) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1131) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1171) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1183) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1202) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1240) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1362) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1527) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1577) = 0._r8
         mat(k,1578) = 0._r8
         mat(k,1634) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1668) = 0._r8
         mat(k,1679) = 0._r8
         mat(k,1702) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1745) = 0._r8
         mat(k,1748) = 0._r8
         mat(k,1751) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1760) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1769) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1792) = 0._r8
         mat(k,1802) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1821) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1833) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1841) = 0._r8
         mat(k,1842) = 0._r8
         mat(k,1843) = 0._r8
         mat(k,1845) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1854) = 0._r8
         mat(k,1856) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1870) = 0._r8
         mat(k,1877) = 0._r8
         mat(k,1878) = 0._r8
         mat(k,1881) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1896) = 0._r8
         mat(k,1897) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1901) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1904) = 0._r8
         mat(k,1908) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1921) = 0._r8
         mat(k,1922) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1927) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1931) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1937) = 0._r8
         mat(k,1938) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1940) = 0._r8
         mat(k,1941) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1998) = 0._r8
         mat(k,2005) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2016) = 0._r8
         mat(k,2017) = 0._r8
         mat(k,2022) = 0._r8
         mat(k,2041) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2055) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2080) = 0._r8
         mat(k,2081) = 0._r8
         mat(k,2091) = 0._r8
         mat(k,2096) = 0._r8
         mat(k,2105) = 0._r8
         mat(k,2112) = 0._r8
         mat(k,2132) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2139) = 0._r8
         mat(k,2144) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2153) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2160) = 0._r8
         mat(k,2162) = 0._r8
         mat(k,2163) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2166) = 0._r8
         mat(k,2167) = 0._r8
         mat(k,2168) = 0._r8
         mat(k,2169) = 0._r8
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
         mat(k, 27) = mat(k, 27) - dti(k)
         mat(k, 28) = mat(k, 28) - dti(k)
         mat(k, 29) = mat(k, 29) - dti(k)
         mat(k, 30) = mat(k, 30) - dti(k)
         mat(k, 31) = mat(k, 31) - dti(k)
         mat(k, 32) = mat(k, 32) - dti(k)
         mat(k, 33) = mat(k, 33) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 45) = mat(k, 45) - dti(k)
         mat(k, 46) = mat(k, 46) - dti(k)
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 52) = mat(k, 52) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 59) = mat(k, 59) - dti(k)
         mat(k, 62) = mat(k, 62) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 68) = mat(k, 68) - dti(k)
         mat(k, 71) = mat(k, 71) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 77) = mat(k, 77) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 121) = mat(k, 121) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 155) = mat(k, 155) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 191) = mat(k, 191) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 212) = mat(k, 212) - dti(k)
         mat(k, 215) = mat(k, 215) - dti(k)
         mat(k, 220) = mat(k, 220) - dti(k)
         mat(k, 225) = mat(k, 225) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 251) = mat(k, 251) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 267) = mat(k, 267) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 297) = mat(k, 297) - dti(k)
         mat(k, 303) = mat(k, 303) - dti(k)
         mat(k, 309) = mat(k, 309) - dti(k)
         mat(k, 315) = mat(k, 315) - dti(k)
         mat(k, 323) = mat(k, 323) - dti(k)
         mat(k, 329) = mat(k, 329) - dti(k)
         mat(k, 336) = mat(k, 336) - dti(k)
         mat(k, 342) = mat(k, 342) - dti(k)
         mat(k, 345) = mat(k, 345) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 357) = mat(k, 357) - dti(k)
         mat(k, 361) = mat(k, 361) - dti(k)
         mat(k, 368) = mat(k, 368) - dti(k)
         mat(k, 377) = mat(k, 377) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 392) = mat(k, 392) - dti(k)
         mat(k, 400) = mat(k, 400) - dti(k)
         mat(k, 406) = mat(k, 406) - dti(k)
         mat(k, 412) = mat(k, 412) - dti(k)
         mat(k, 417) = mat(k, 417) - dti(k)
         mat(k, 422) = mat(k, 422) - dti(k)
         mat(k, 430) = mat(k, 430) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 446) = mat(k, 446) - dti(k)
         mat(k, 450) = mat(k, 450) - dti(k)
         mat(k, 458) = mat(k, 458) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 474) = mat(k, 474) - dti(k)
         mat(k, 481) = mat(k, 481) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 501) = mat(k, 501) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 514) = mat(k, 514) - dti(k)
         mat(k, 521) = mat(k, 521) - dti(k)
         mat(k, 529) = mat(k, 529) - dti(k)
         mat(k, 536) = mat(k, 536) - dti(k)
         mat(k, 547) = mat(k, 547) - dti(k)
         mat(k, 558) = mat(k, 558) - dti(k)
         mat(k, 567) = mat(k, 567) - dti(k)
         mat(k, 578) = mat(k, 578) - dti(k)
         mat(k, 589) = mat(k, 589) - dti(k)
         mat(k, 600) = mat(k, 600) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 618) = mat(k, 618) - dti(k)
         mat(k, 634) = mat(k, 634) - dti(k)
         mat(k, 645) = mat(k, 645) - dti(k)
         mat(k, 654) = mat(k, 654) - dti(k)
         mat(k, 664) = mat(k, 664) - dti(k)
         mat(k, 673) = mat(k, 673) - dti(k)
         mat(k, 681) = mat(k, 681) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 701) = mat(k, 701) - dti(k)
         mat(k, 708) = mat(k, 708) - dti(k)
         mat(k, 712) = mat(k, 712) - dti(k)
         mat(k, 717) = mat(k, 717) - dti(k)
         mat(k, 728) = mat(k, 728) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 753) = mat(k, 753) - dti(k)
         mat(k, 762) = mat(k, 762) - dti(k)
         mat(k, 769) = mat(k, 769) - dti(k)
         mat(k, 781) = mat(k, 781) - dti(k)
         mat(k, 798) = mat(k, 798) - dti(k)
         mat(k, 803) = mat(k, 803) - dti(k)
         mat(k, 820) = mat(k, 820) - dti(k)
         mat(k, 840) = mat(k, 840) - dti(k)
         mat(k, 858) = mat(k, 858) - dti(k)
         mat(k, 880) = mat(k, 880) - dti(k)
         mat(k, 892) = mat(k, 892) - dti(k)
         mat(k, 898) = mat(k, 898) - dti(k)
         mat(k, 906) = mat(k, 906) - dti(k)
         mat(k, 916) = mat(k, 916) - dti(k)
         mat(k, 928) = mat(k, 928) - dti(k)
         mat(k, 943) = mat(k, 943) - dti(k)
         mat(k, 960) = mat(k, 960) - dti(k)
         mat(k, 968) = mat(k, 968) - dti(k)
         mat(k, 982) = mat(k, 982) - dti(k)
         mat(k, 991) = mat(k, 991) - dti(k)
         mat(k, 997) = mat(k, 997) - dti(k)
         mat(k,1008) = mat(k,1008) - dti(k)
         mat(k,1030) = mat(k,1030) - dti(k)
         mat(k,1049) = mat(k,1049) - dti(k)
         mat(k,1065) = mat(k,1065) - dti(k)
         mat(k,1076) = mat(k,1076) - dti(k)
         mat(k,1091) = mat(k,1091) - dti(k)
         mat(k,1104) = mat(k,1104) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1140) = mat(k,1140) - dti(k)
         mat(k,1152) = mat(k,1152) - dti(k)
         mat(k,1172) = mat(k,1172) - dti(k)
         mat(k,1206) = mat(k,1206) - dti(k)
         mat(k,1231) = mat(k,1231) - dti(k)
         mat(k,1251) = mat(k,1251) - dti(k)
         mat(k,1271) = mat(k,1271) - dti(k)
         mat(k,1302) = mat(k,1302) - dti(k)
         mat(k,1318) = mat(k,1318) - dti(k)
         mat(k,1337) = mat(k,1337) - dti(k)
         mat(k,1352) = mat(k,1352) - dti(k)
         mat(k,1396) = mat(k,1396) - dti(k)
         mat(k,1427) = mat(k,1427) - dti(k)
         mat(k,1507) = mat(k,1507) - dti(k)
         mat(k,1543) = mat(k,1543) - dti(k)
         mat(k,1570) = mat(k,1570) - dti(k)
         mat(k,1719) = mat(k,1719) - dti(k)
         mat(k,1782) = mat(k,1782) - dti(k)
         mat(k,1824) = mat(k,1824) - dti(k)
         mat(k,1849) = mat(k,1849) - dti(k)
         mat(k,1907) = mat(k,1907) - dti(k)
         mat(k,1928) = mat(k,1928) - dti(k)
         mat(k,1951) = mat(k,1951) - dti(k)
         mat(k,1994) = mat(k,1994) - dti(k)
         mat(k,2019) = mat(k,2019) - dti(k)
         mat(k,2115) = mat(k,2115) - dti(k)
         mat(k,2143) = mat(k,2143) - dti(k)
         mat(k,2170) = mat(k,2170) - dti(k)
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
