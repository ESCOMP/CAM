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
         mat(k,640) = -(rxt(k,396)*y(k,226))
         mat(k,1734) = -rxt(k,396)*y(k,1)
         mat(k,2020) = rxt(k,399)*y(k,190)
         mat(k,1036) = rxt(k,399)*y(k,124)
         mat(k,685) = -(rxt(k,400)*y(k,226))
         mat(k,1737) = -rxt(k,400)*y(k,2)
         mat(k,2147) = rxt(k,397)*y(k,190)
         mat(k,1037) = rxt(k,397)*y(k,90)
         mat(k,988) = -(rxt(k,479)*y(k,126) + rxt(k,480)*y(k,135) + rxt(k,481) &
                      *y(k,226))
         mat(k,1869) = -rxt(k,479)*y(k,6)
         mat(k,2375) = -rxt(k,480)*y(k,6)
         mat(k,1762) = -rxt(k,481)*y(k,6)
         mat(k,158) = -(rxt(k,438)*y(k,226))
         mat(k,1666) = -rxt(k,438)*y(k,7)
         mat(k,401) = -(rxt(k,441)*y(k,226))
         mat(k,1704) = -rxt(k,441)*y(k,8)
         mat(k,2125) = rxt(k,439)*y(k,192)
         mat(k,484) = rxt(k,439)*y(k,90)
         mat(k,159) = .120_r8*rxt(k,438)*y(k,226)
         mat(k,1667) = .120_r8*rxt(k,438)*y(k,7)
         mat(k,986) = .100_r8*rxt(k,480)*y(k,135)
         mat(k,1014) = .100_r8*rxt(k,483)*y(k,135)
         mat(k,2365) = .100_r8*rxt(k,480)*y(k,6) + .100_r8*rxt(k,483)*y(k,110)
         mat(k,2007) = .500_r8*rxt(k,440)*y(k,192) + .200_r8*rxt(k,467)*y(k,233) &
                      + .060_r8*rxt(k,473)*y(k,236)
         mat(k,485) = .500_r8*rxt(k,440)*y(k,124)
         mat(k,745) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,761) = .060_r8*rxt(k,473)*y(k,124)
         mat(k,2001) = .200_r8*rxt(k,467)*y(k,233) + .200_r8*rxt(k,473)*y(k,236)
         mat(k,744) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,759) = .200_r8*rxt(k,473)*y(k,124)
         mat(k,2017) = .200_r8*rxt(k,467)*y(k,233) + .150_r8*rxt(k,473)*y(k,236)
         mat(k,746) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,762) = .150_r8*rxt(k,473)*y(k,124)
         mat(k,2003) = .210_r8*rxt(k,473)*y(k,236)
         mat(k,760) = .210_r8*rxt(k,473)*y(k,124)
         mat(k,233) = -(rxt(k,401)*y(k,226))
         mat(k,1679) = -rxt(k,401)*y(k,15)
         mat(k,985) = .050_r8*rxt(k,480)*y(k,135)
         mat(k,1013) = .050_r8*rxt(k,483)*y(k,135)
         mat(k,2364) = .050_r8*rxt(k,480)*y(k,6) + .050_r8*rxt(k,483)*y(k,110)
         mat(k,353) = -(rxt(k,367)*y(k,126) + rxt(k,368)*y(k,226))
         mat(k,1863) = -rxt(k,367)*y(k,16)
         mat(k,1698) = -rxt(k,368)*y(k,16)
         mat(k,1509) = -(rxt(k,250)*y(k,42) + rxt(k,251)*y(k,90) + rxt(k,252)*y(k,135))
         mat(k,2338) = -rxt(k,250)*y(k,17)
         mat(k,2192) = -rxt(k,251)*y(k,17)
         mat(k,2402) = -rxt(k,252)*y(k,17)
         mat(k,1561) = 4.000_r8*rxt(k,253)*y(k,19) + (rxt(k,254)+rxt(k,255))*y(k,59) &
                      + rxt(k,258)*y(k,124) + rxt(k,261)*y(k,133) + rxt(k,509) &
                      *y(k,151) + rxt(k,262)*y(k,226)
         mat(k,139) = rxt(k,240)*y(k,222)
         mat(k,145) = rxt(k,266)*y(k,222)
         mat(k,472) = 2.000_r8*rxt(k,277)*y(k,56) + 2.000_r8*rxt(k,289)*y(k,222) &
                      + 2.000_r8*rxt(k,278)*y(k,226)
         mat(k,604) = rxt(k,279)*y(k,56) + rxt(k,290)*y(k,222) + rxt(k,280)*y(k,226)
         mat(k,447) = 3.000_r8*rxt(k,284)*y(k,56) + 3.000_r8*rxt(k,267)*y(k,222) &
                      + 3.000_r8*rxt(k,285)*y(k,226)
         mat(k,1945) = 2.000_r8*rxt(k,277)*y(k,41) + rxt(k,279)*y(k,43) &
                      + 3.000_r8*rxt(k,284)*y(k,55)
         mat(k,1587) = (rxt(k,254)+rxt(k,255))*y(k,19)
         mat(k,107) = 2.000_r8*rxt(k,268)*y(k,222)
         mat(k,827) = rxt(k,263)*y(k,133) + rxt(k,269)*y(k,222) + rxt(k,264)*y(k,226)
         mat(k,2063) = rxt(k,258)*y(k,19)
         mat(k,2312) = rxt(k,261)*y(k,19) + rxt(k,263)*y(k,81)
         mat(k,1475) = rxt(k,509)*y(k,19)
         mat(k,1629) = rxt(k,240)*y(k,34) + rxt(k,266)*y(k,35) + 2.000_r8*rxt(k,289) &
                      *y(k,41) + rxt(k,290)*y(k,43) + 3.000_r8*rxt(k,267)*y(k,55) &
                      + 2.000_r8*rxt(k,268)*y(k,78) + rxt(k,269)*y(k,81)
         mat(k,1794) = rxt(k,262)*y(k,19) + 2.000_r8*rxt(k,278)*y(k,41) + rxt(k,280) &
                      *y(k,43) + 3.000_r8*rxt(k,285)*y(k,55) + rxt(k,264)*y(k,81)
         mat(k,1554) = rxt(k,256)*y(k,59)
         mat(k,1580) = rxt(k,256)*y(k,19)
         mat(k,1489) = (rxt(k,570)+rxt(k,575))*y(k,92)
         mat(k,784) = (rxt(k,570)+rxt(k,575))*y(k,85)
         mat(k,1563) = -(4._r8*rxt(k,253)*y(k,19) + (rxt(k,254) + rxt(k,255) + rxt(k,256) &
                      ) * y(k,59) + rxt(k,257)*y(k,90) + rxt(k,258)*y(k,124) + rxt(k,259) &
                      *y(k,125) + rxt(k,261)*y(k,133) + rxt(k,262)*y(k,226) + rxt(k,509) &
                      *y(k,151))
         mat(k,1589) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,19)
         mat(k,2194) = -rxt(k,257)*y(k,19)
         mat(k,2065) = -rxt(k,258)*y(k,19)
         mat(k,1841) = -rxt(k,259)*y(k,19)
         mat(k,2314) = -rxt(k,261)*y(k,19)
         mat(k,1796) = -rxt(k,262)*y(k,19)
         mat(k,1477) = -rxt(k,509)*y(k,19)
         mat(k,1511) = rxt(k,252)*y(k,135)
         mat(k,567) = rxt(k,260)*y(k,133)
         mat(k,828) = rxt(k,270)*y(k,222)
         mat(k,788) = rxt(k,265)*y(k,133)
         mat(k,2314) = mat(k,2314) + rxt(k,260)*y(k,20) + rxt(k,265)*y(k,92)
         mat(k,2404) = rxt(k,252)*y(k,17)
         mat(k,1631) = rxt(k,270)*y(k,81)
         mat(k,564) = -(rxt(k,260)*y(k,133))
         mat(k,2293) = -rxt(k,260)*y(k,20)
         mat(k,1556) = rxt(k,259)*y(k,125)
         mat(k,1820) = rxt(k,259)*y(k,19)
         mat(k,245) = -(rxt(k,442)*y(k,226))
         mat(k,1682) = -rxt(k,442)*y(k,22)
         mat(k,2000) = rxt(k,445)*y(k,194)
         mat(k,431) = rxt(k,445)*y(k,124)
         mat(k,343) = -(rxt(k,444)*y(k,226))
         mat(k,1696) = -rxt(k,444)*y(k,23)
         mat(k,2121) = rxt(k,443)*y(k,194)
         mat(k,432) = rxt(k,443)*y(k,90)
         mat(k,290) = -(rxt(k,315)*y(k,56) + rxt(k,316)*y(k,226))
         mat(k,1919) = -rxt(k,315)*y(k,24)
         mat(k,1688) = -rxt(k,316)*y(k,24)
         mat(k,548) = -(rxt(k,317)*y(k,56) + rxt(k,318)*y(k,135) + rxt(k,343)*y(k,226))
         mat(k,1925) = -rxt(k,317)*y(k,25)
         mat(k,2367) = -rxt(k,318)*y(k,25)
         mat(k,1723) = -rxt(k,343)*y(k,25)
         mat(k,263) = -(rxt(k,323)*y(k,226))
         mat(k,1685) = -rxt(k,323)*y(k,26)
         mat(k,896) = .800_r8*rxt(k,319)*y(k,195) + .200_r8*rxt(k,320)*y(k,199)
         mat(k,2211) = .200_r8*rxt(k,320)*y(k,195)
         mat(k,348) = -(rxt(k,324)*y(k,226))
         mat(k,1697) = -rxt(k,324)*y(k,27)
         mat(k,2122) = rxt(k,321)*y(k,195)
         mat(k,897) = rxt(k,321)*y(k,90)
         mat(k,296) = -(rxt(k,325)*y(k,56) + rxt(k,326)*y(k,226))
         mat(k,1920) = -rxt(k,325)*y(k,28)
         mat(k,1689) = -rxt(k,326)*y(k,28)
         mat(k,1131) = -(rxt(k,346)*y(k,126) + rxt(k,347)*y(k,135) + rxt(k,365) &
                      *y(k,226))
         mat(k,1879) = -rxt(k,346)*y(k,29)
         mat(k,2384) = -rxt(k,347)*y(k,29)
         mat(k,1772) = -rxt(k,365)*y(k,29)
         mat(k,882) = .130_r8*rxt(k,425)*y(k,135)
         mat(k,2384) = mat(k,2384) + .130_r8*rxt(k,425)*y(k,99)
         mat(k,413) = -(rxt(k,351)*y(k,226))
         mat(k,1706) = -rxt(k,351)*y(k,30)
         mat(k,2127) = rxt(k,349)*y(k,196)
         mat(k,932) = rxt(k,349)*y(k,90)
         mat(k,302) = -(rxt(k,352)*y(k,226) + rxt(k,355)*y(k,56))
         mat(k,1690) = -rxt(k,352)*y(k,31)
         mat(k,1921) = -rxt(k,355)*y(k,31)
         mat(k,267) = -(rxt(k,448)*y(k,226))
         mat(k,1686) = -rxt(k,448)*y(k,32)
         mat(k,2117) = rxt(k,446)*y(k,197)
         mat(k,631) = rxt(k,446)*y(k,90)
         mat(k,99) = -(rxt(k,239)*y(k,222))
         mat(k,1605) = -rxt(k,239)*y(k,33)
         mat(k,137) = -(rxt(k,240)*y(k,222))
         mat(k,1610) = -rxt(k,240)*y(k,34)
         mat(k,142) = -(rxt(k,266)*y(k,222))
         mat(k,1611) = -rxt(k,266)*y(k,35)
         mat(k,109) = -(rxt(k,241)*y(k,222))
         mat(k,1607) = -rxt(k,241)*y(k,36)
         mat(k,147) = -(rxt(k,242)*y(k,222))
         mat(k,1612) = -rxt(k,242)*y(k,37)
         mat(k,113) = -(rxt(k,243)*y(k,222))
         mat(k,1608) = -rxt(k,243)*y(k,38)
         mat(k,152) = -(rxt(k,244)*y(k,222))
         mat(k,1613) = -rxt(k,244)*y(k,39)
         mat(k,117) = -(rxt(k,245)*y(k,222))
         mat(k,1609) = -rxt(k,245)*y(k,40)
         mat(k,470) = -(rxt(k,277)*y(k,56) + rxt(k,278)*y(k,226) + rxt(k,289)*y(k,222))
         mat(k,1924) = -rxt(k,277)*y(k,41)
         mat(k,1714) = -rxt(k,278)*y(k,41)
         mat(k,1623) = -rxt(k,289)*y(k,41)
         mat(k,2354) = -(rxt(k,214)*y(k,56) + rxt(k,250)*y(k,17) + rxt(k,294)*y(k,90) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,133) + rxt(k,297) &
                      *y(k,226))
         mat(k,1961) = -rxt(k,214)*y(k,42)
         mat(k,1519) = -rxt(k,250)*y(k,42)
         mat(k,2208) = -rxt(k,294)*y(k,42)
         mat(k,1915) = -rxt(k,295)*y(k,42)
         mat(k,2328) = -rxt(k,296)*y(k,42)
         mat(k,1810) = -rxt(k,297)*y(k,42)
         mat(k,649) = .400_r8*rxt(k,396)*y(k,226)
         mat(k,1004) = .340_r8*rxt(k,480)*y(k,135)
         mat(k,360) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,554) = rxt(k,318)*y(k,135)
         mat(k,1146) = .500_r8*rxt(k,347)*y(k,135)
         mat(k,623) = .500_r8*rxt(k,335)*y(k,226)
         mat(k,814) = rxt(k,302)*y(k,226)
         mat(k,393) = .300_r8*rxt(k,303)*y(k,226)
         mat(k,2285) = (rxt(k,311)+rxt(k,312))*y(k,222)
         mat(k,1602) = rxt(k,221)*y(k,199)
         mat(k,1168) = .800_r8*rxt(k,340)*y(k,226)
         mat(k,2208) = mat(k,2208) + .450_r8*rxt(k,383)*y(k,213) + .150_r8*rxt(k,362) &
                      *y(k,230)
         mat(k,894) = .910_r8*rxt(k,425)*y(k,135)
         mat(k,597) = .300_r8*rxt(k,416)*y(k,226)
         mat(k,1275) = .120_r8*rxt(k,378)*y(k,135)
         mat(k,588) = .500_r8*rxt(k,391)*y(k,226)
         mat(k,1032) = .340_r8*rxt(k,483)*y(k,135)
         mat(k,1387) = .600_r8*rxt(k,392)*y(k,135)
         mat(k,2079) = .100_r8*rxt(k,398)*y(k,190) + rxt(k,301)*y(k,199) &
                      + .500_r8*rxt(k,369)*y(k,202) + .500_r8*rxt(k,337)*y(k,204) &
                      + .920_r8*rxt(k,408)*y(k,206) + .250_r8*rxt(k,376)*y(k,211) &
                      + rxt(k,385)*y(k,213) + rxt(k,359)*y(k,229) + rxt(k,363) &
                      *y(k,230) + .340_r8*rxt(k,492)*y(k,231) + .320_r8*rxt(k,497) &
                      *y(k,232) + .250_r8*rxt(k,433)*y(k,235)
         mat(k,1915) = mat(k,1915) + .500_r8*rxt(k,367)*y(k,16) + rxt(k,409)*y(k,206) &
                      + .250_r8*rxt(k,375)*y(k,211) + rxt(k,386)*y(k,213)
         mat(k,2418) = .340_r8*rxt(k,480)*y(k,6) + rxt(k,318)*y(k,25) &
                      + .500_r8*rxt(k,347)*y(k,29) + .910_r8*rxt(k,425)*y(k,99) &
                      + .120_r8*rxt(k,378)*y(k,105) + .340_r8*rxt(k,483)*y(k,110) &
                      + .600_r8*rxt(k,392)*y(k,111)
         mat(k,539) = rxt(k,342)*y(k,226)
         mat(k,1109) = .680_r8*rxt(k,501)*y(k,226)
         mat(k,1050) = .100_r8*rxt(k,398)*y(k,124)
         mat(k,907) = .700_r8*rxt(k,320)*y(k,199)
         mat(k,942) = rxt(k,348)*y(k,199)
         mat(k,1437) = rxt(k,331)*y(k,199) + rxt(k,405)*y(k,206) + .250_r8*rxt(k,372) &
                      *y(k,211) + rxt(k,381)*y(k,213) + .250_r8*rxt(k,430)*y(k,235)
         mat(k,2260) = rxt(k,221)*y(k,59) + rxt(k,301)*y(k,124) + .700_r8*rxt(k,320) &
                      *y(k,195) + rxt(k,348)*y(k,196) + rxt(k,331)*y(k,198) + ( &
                      + 4.000_r8*rxt(k,298)+2.000_r8*rxt(k,299))*y(k,199) &
                      + 1.500_r8*rxt(k,406)*y(k,206) + .750_r8*rxt(k,411)*y(k,207) &
                      + .800_r8*rxt(k,420)*y(k,208) + .880_r8*rxt(k,373)*y(k,211) &
                      + 2.000_r8*rxt(k,382)*y(k,213) + .750_r8*rxt(k,485)*y(k,221) &
                      + .800_r8*rxt(k,361)*y(k,230) + .930_r8*rxt(k,490)*y(k,231) &
                      + .950_r8*rxt(k,495)*y(k,232) + .800_r8*rxt(k,431)*y(k,235)
         mat(k,579) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,802) = .500_r8*rxt(k,337)*y(k,124)
         mat(k,1310) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) + rxt(k,405) &
                      *y(k,198) + 1.500_r8*rxt(k,406)*y(k,199)
         mat(k,1343) = .750_r8*rxt(k,411)*y(k,199)
         mat(k,1264) = .800_r8*rxt(k,420)*y(k,199)
         mat(k,1365) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,198) + .880_r8*rxt(k,373)*y(k,199)
         mat(k,1405) = .450_r8*rxt(k,383)*y(k,90) + rxt(k,385)*y(k,124) + rxt(k,386) &
                      *y(k,126) + rxt(k,381)*y(k,198) + 2.000_r8*rxt(k,382)*y(k,199) &
                      + 4.000_r8*rxt(k,384)*y(k,213)
         mat(k,1099) = .750_r8*rxt(k,485)*y(k,199)
         mat(k,1645) = (rxt(k,311)+rxt(k,312))*y(k,54)
         mat(k,1810) = mat(k,1810) + .400_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,335) &
                      *y(k,51) + rxt(k,302)*y(k,52) + .300_r8*rxt(k,303)*y(k,53) &
                      + .800_r8*rxt(k,340)*y(k,74) + .300_r8*rxt(k,416)*y(k,100) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,342)*y(k,140) &
                      + .680_r8*rxt(k,501)*y(k,179)
         mat(k,824) = rxt(k,359)*y(k,124)
         mat(k,1223) = .150_r8*rxt(k,362)*y(k,90) + rxt(k,363)*y(k,124) &
                      + .800_r8*rxt(k,361)*y(k,199)
         mat(k,1186) = .340_r8*rxt(k,492)*y(k,124) + .930_r8*rxt(k,490)*y(k,199)
         mat(k,1069) = .320_r8*rxt(k,497)*y(k,124) + .950_r8*rxt(k,495)*y(k,199)
         mat(k,1241) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,430)*y(k,198) &
                      + .800_r8*rxt(k,431)*y(k,199)
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
         mat(k,602) = -(rxt(k,279)*y(k,56) + rxt(k,280)*y(k,226) + rxt(k,290)*y(k,222))
         mat(k,1926) = -rxt(k,279)*y(k,43)
         mat(k,1729) = -rxt(k,280)*y(k,43)
         mat(k,1624) = -rxt(k,290)*y(k,43)
         mat(k,121) = -(rxt(k,281)*y(k,226))
         mat(k,1664) = -rxt(k,281)*y(k,44)
         mat(k,1149) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,226))
         mat(k,1880) = -rxt(k,327)*y(k,45)
         mat(k,1773) = -rxt(k,328)*y(k,45)
         mat(k,644) = .800_r8*rxt(k,396)*y(k,226)
         mat(k,356) = rxt(k,367)*y(k,126)
         mat(k,264) = rxt(k,323)*y(k,226)
         mat(k,350) = .500_r8*rxt(k,324)*y(k,226)
         mat(k,1132) = .500_r8*rxt(k,347)*y(k,135)
         mat(k,2173) = .200_r8*rxt(k,387)*y(k,215)
         mat(k,1369) = .100_r8*rxt(k,392)*y(k,135)
         mat(k,2045) = .400_r8*rxt(k,398)*y(k,190) + rxt(k,322)*y(k,195) &
                      + .270_r8*rxt(k,350)*y(k,196) + rxt(k,369)*y(k,202) + rxt(k,388) &
                      *y(k,215) + rxt(k,359)*y(k,229)
         mat(k,1880) = mat(k,1880) + rxt(k,367)*y(k,16)
         mat(k,2385) = .500_r8*rxt(k,347)*y(k,29) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,1042) = .400_r8*rxt(k,398)*y(k,124)
         mat(k,900) = rxt(k,322)*y(k,124) + 3.200_r8*rxt(k,319)*y(k,195) &
                      + .800_r8*rxt(k,320)*y(k,199)
         mat(k,935) = .270_r8*rxt(k,350)*y(k,124)
         mat(k,2228) = .800_r8*rxt(k,320)*y(k,195)
         mat(k,574) = rxt(k,369)*y(k,124)
         mat(k,697) = .200_r8*rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124)
         mat(k,1773) = mat(k,1773) + .800_r8*rxt(k,396)*y(k,1) + rxt(k,323)*y(k,26) &
                      + .500_r8*rxt(k,324)*y(k,27)
         mat(k,817) = rxt(k,359)*y(k,124)
         mat(k,369) = -(rxt(k,282)*y(k,56) + rxt(k,283)*y(k,226))
         mat(k,1922) = -rxt(k,282)*y(k,46)
         mat(k,1700) = -rxt(k,283)*y(k,46)
         mat(k,102) = -(rxt(k,329)*y(k,226))
         mat(k,1663) = -rxt(k,329)*y(k,47)
         mat(k,1078) = -(rxt(k,366)*y(k,226))
         mat(k,1768) = -rxt(k,366)*y(k,48)
         mat(k,643) = .800_r8*rxt(k,396)*y(k,226)
         mat(k,993) = .520_r8*rxt(k,480)*y(k,135)
         mat(k,355) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,1021) = .520_r8*rxt(k,483)*y(k,135)
         mat(k,2041) = .250_r8*rxt(k,398)*y(k,190) + .820_r8*rxt(k,350)*y(k,196) &
                      + .500_r8*rxt(k,369)*y(k,202) + .270_r8*rxt(k,492)*y(k,231) &
                      + .040_r8*rxt(k,497)*y(k,232)
         mat(k,1875) = .500_r8*rxt(k,367)*y(k,16)
         mat(k,2381) = .520_r8*rxt(k,480)*y(k,6) + .520_r8*rxt(k,483)*y(k,110)
         mat(k,1100) = .500_r8*rxt(k,501)*y(k,226)
         mat(k,1041) = .250_r8*rxt(k,398)*y(k,124)
         mat(k,934) = .820_r8*rxt(k,350)*y(k,124) + .820_r8*rxt(k,348)*y(k,199)
         mat(k,2224) = .820_r8*rxt(k,348)*y(k,196) + .150_r8*rxt(k,490)*y(k,231) &
                      + .025_r8*rxt(k,495)*y(k,232)
         mat(k,573) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,1768) = mat(k,1768) + .800_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,501) &
                      *y(k,179)
         mat(k,1172) = .270_r8*rxt(k,492)*y(k,124) + .150_r8*rxt(k,490)*y(k,199)
         mat(k,1062) = .040_r8*rxt(k,497)*y(k,124) + .025_r8*rxt(k,495)*y(k,199)
         mat(k,1279) = -(rxt(k,353)*y(k,126) + rxt(k,354)*y(k,226))
         mat(k,1890) = -rxt(k,353)*y(k,49)
         mat(k,1783) = -rxt(k,354)*y(k,49)
         mat(k,2182) = .070_r8*rxt(k,450)*y(k,200) + .070_r8*rxt(k,456)*y(k,214)
         mat(k,1207) = rxt(k,356)*y(k,226)
         mat(k,1268) = .880_r8*rxt(k,378)*y(k,135)
         mat(k,1372) = .500_r8*rxt(k,392)*y(k,135)
         mat(k,2055) = .170_r8*rxt(k,451)*y(k,200) + .050_r8*rxt(k,414)*y(k,207) &
                      + .250_r8*rxt(k,376)*y(k,211) + .170_r8*rxt(k,457)*y(k,214) &
                      + .400_r8*rxt(k,467)*y(k,233) + .250_r8*rxt(k,433)*y(k,235) &
                      + .540_r8*rxt(k,473)*y(k,236) + .510_r8*rxt(k,476)*y(k,238)
         mat(k,1890) = mat(k,1890) + .050_r8*rxt(k,415)*y(k,207) + .250_r8*rxt(k,375) &
                      *y(k,211) + .250_r8*rxt(k,434)*y(k,235)
         mat(k,872) = rxt(k,357)*y(k,226)
         mat(k,2393) = .880_r8*rxt(k,378)*y(k,105) + .500_r8*rxt(k,392)*y(k,111)
         mat(k,1420) = .250_r8*rxt(k,372)*y(k,211) + .250_r8*rxt(k,430)*y(k,235)
         mat(k,2237) = .240_r8*rxt(k,373)*y(k,211) + .500_r8*rxt(k,361)*y(k,230) &
                      + .100_r8*rxt(k,431)*y(k,235)
         mat(k,778) = .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451)*y(k,124)
         mat(k,1329) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1353) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,198) + .240_r8*rxt(k,373)*y(k,199)
         mat(k,916) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,1783) = mat(k,1783) + rxt(k,356)*y(k,96) + rxt(k,357)*y(k,127)
         mat(k,1216) = .500_r8*rxt(k,361)*y(k,199)
         mat(k,754) = .400_r8*rxt(k,467)*y(k,124)
         mat(k,1232) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,770) = .540_r8*rxt(k,473)*y(k,124)
         mat(k,510) = .510_r8*rxt(k,476)*y(k,124)
         mat(k,703) = -(rxt(k,334)*y(k,226))
         mat(k,1739) = -rxt(k,334)*y(k,50)
         mat(k,1127) = .120_r8*rxt(k,347)*y(k,135)
         mat(k,2149) = .150_r8*rxt(k,332)*y(k,198) + .150_r8*rxt(k,383)*y(k,213)
         mat(k,2369) = .120_r8*rxt(k,347)*y(k,29)
         mat(k,1411) = .150_r8*rxt(k,332)*y(k,90) + .100_r8*rxt(k,331)*y(k,199)
         mat(k,2216) = .100_r8*rxt(k,331)*y(k,198)
         mat(k,1392) = .150_r8*rxt(k,383)*y(k,90)
         mat(k,618) = -(rxt(k,335)*y(k,226))
         mat(k,1731) = -rxt(k,335)*y(k,51)
         mat(k,2143) = .400_r8*rxt(k,332)*y(k,198) + .400_r8*rxt(k,383)*y(k,213)
         mat(k,1410) = .400_r8*rxt(k,332)*y(k,90)
         mat(k,1391) = .400_r8*rxt(k,383)*y(k,90)
         mat(k,811) = -(rxt(k,302)*y(k,226))
         mat(k,1748) = -rxt(k,302)*y(k,52)
         mat(k,898) = .300_r8*rxt(k,320)*y(k,199)
         mat(k,2217) = .300_r8*rxt(k,320)*y(k,195) + 2.000_r8*rxt(k,299)*y(k,199) &
                      + .250_r8*rxt(k,406)*y(k,206) + .250_r8*rxt(k,411)*y(k,207) &
                      + .200_r8*rxt(k,420)*y(k,208) + .250_r8*rxt(k,373)*y(k,211) &
                      + .250_r8*rxt(k,485)*y(k,221) + .500_r8*rxt(k,361)*y(k,230) &
                      + .250_r8*rxt(k,490)*y(k,231) + .250_r8*rxt(k,495)*y(k,232) &
                      + .300_r8*rxt(k,431)*y(k,235)
         mat(k,1289) = .250_r8*rxt(k,406)*y(k,199)
         mat(k,1318) = .250_r8*rxt(k,411)*y(k,199)
         mat(k,1245) = .200_r8*rxt(k,420)*y(k,199)
         mat(k,1347) = .250_r8*rxt(k,373)*y(k,199)
         mat(k,1086) = .250_r8*rxt(k,485)*y(k,199)
         mat(k,1213) = .500_r8*rxt(k,361)*y(k,199)
         mat(k,1171) = .250_r8*rxt(k,490)*y(k,199)
         mat(k,1059) = .250_r8*rxt(k,495)*y(k,199)
         mat(k,1226) = .300_r8*rxt(k,431)*y(k,199)
         mat(k,389) = -(rxt(k,303)*y(k,226))
         mat(k,1702) = -rxt(k,303)*y(k,53)
         mat(k,2123) = rxt(k,300)*y(k,199)
         mat(k,2214) = rxt(k,300)*y(k,90)
         mat(k,2283) = -(rxt(k,215)*y(k,56) + rxt(k,271)*y(k,73) + rxt(k,304)*y(k,226) &
                      + (rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,222))
         mat(k,1959) = -rxt(k,215)*y(k,54)
         mat(k,930) = -rxt(k,271)*y(k,54)
         mat(k,1808) = -rxt(k,304)*y(k,54)
         mat(k,1643) = -(rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,54)
         mat(k,1145) = .100_r8*rxt(k,347)*y(k,135)
         mat(k,2416) = .100_r8*rxt(k,347)*y(k,29)
         mat(k,445) = -(rxt(k,267)*y(k,222) + rxt(k,284)*y(k,56) + rxt(k,285)*y(k,226))
         mat(k,1622) = -rxt(k,267)*y(k,55)
         mat(k,1923) = -rxt(k,284)*y(k,55)
         mat(k,1710) = -rxt(k,285)*y(k,55)
         mat(k,1953) = -(rxt(k,214)*y(k,42) + rxt(k,215)*y(k,54) + rxt(k,216)*y(k,77) &
                      + rxt(k,217)*y(k,79) + (rxt(k,218) + rxt(k,219)) * y(k,90) &
                      + rxt(k,220)*y(k,135) + rxt(k,227)*y(k,60) + rxt(k,236)*y(k,93) &
                      + rxt(k,277)*y(k,41) + rxt(k,279)*y(k,43) + rxt(k,282)*y(k,46) &
                      + rxt(k,284)*y(k,55) + rxt(k,325)*y(k,28) + rxt(k,355)*y(k,31))
         mat(k,2346) = -rxt(k,214)*y(k,56)
         mat(k,2277) = -rxt(k,215)*y(k,56)
         mat(k,1464) = -rxt(k,216)*y(k,56)
         mat(k,614) = -rxt(k,217)*y(k,56)
         mat(k,2200) = -(rxt(k,218) + rxt(k,219)) * y(k,56)
         mat(k,2410) = -rxt(k,220)*y(k,56)
         mat(k,959) = -rxt(k,227)*y(k,56)
         mat(k,839) = -rxt(k,236)*y(k,56)
         mat(k,475) = -rxt(k,277)*y(k,56)
         mat(k,607) = -rxt(k,279)*y(k,56)
         mat(k,373) = -rxt(k,282)*y(k,56)
         mat(k,450) = -rxt(k,284)*y(k,56)
         mat(k,300) = -rxt(k,325)*y(k,56)
         mat(k,306) = -rxt(k,355)*y(k,56)
         mat(k,1569) = rxt(k,255)*y(k,59)
         mat(k,101) = 4.000_r8*rxt(k,239)*y(k,222)
         mat(k,141) = rxt(k,240)*y(k,222)
         mat(k,112) = 2.000_r8*rxt(k,241)*y(k,222)
         mat(k,151) = 2.000_r8*rxt(k,242)*y(k,222)
         mat(k,116) = 2.000_r8*rxt(k,243)*y(k,222)
         mat(k,156) = rxt(k,244)*y(k,222)
         mat(k,120) = 2.000_r8*rxt(k,245)*y(k,222)
         mat(k,123) = 3.000_r8*rxt(k,281)*y(k,226)
         mat(k,373) = mat(k,373) + rxt(k,283)*y(k,226)
         mat(k,1595) = rxt(k,255)*y(k,19) + (4.000_r8*rxt(k,222)+2.000_r8*rxt(k,224)) &
                      *y(k,59) + rxt(k,226)*y(k,124) + rxt(k,231)*y(k,133) &
                      + rxt(k,510)*y(k,151) + rxt(k,221)*y(k,199) + rxt(k,232) &
                      *y(k,226)
         mat(k,227) = rxt(k,276)*y(k,222)
         mat(k,223) = rxt(k,291)*y(k,222) + rxt(k,286)*y(k,226)
         mat(k,253) = rxt(k,292)*y(k,222) + rxt(k,287)*y(k,226)
         mat(k,276) = rxt(k,293)*y(k,222) + rxt(k,288)*y(k,226)
         mat(k,1501) = rxt(k,234)*y(k,133) + rxt(k,246)*y(k,222) + rxt(k,235)*y(k,226)
         mat(k,2071) = rxt(k,226)*y(k,59)
         mat(k,2320) = rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85)
         mat(k,1481) = rxt(k,510)*y(k,59)
         mat(k,2252) = rxt(k,221)*y(k,59)
         mat(k,1637) = 4.000_r8*rxt(k,239)*y(k,33) + rxt(k,240)*y(k,34) &
                      + 2.000_r8*rxt(k,241)*y(k,36) + 2.000_r8*rxt(k,242)*y(k,37) &
                      + 2.000_r8*rxt(k,243)*y(k,38) + rxt(k,244)*y(k,39) &
                      + 2.000_r8*rxt(k,245)*y(k,40) + rxt(k,276)*y(k,65) + rxt(k,291) &
                      *y(k,82) + rxt(k,292)*y(k,83) + rxt(k,293)*y(k,84) + rxt(k,246) &
                      *y(k,85)
         mat(k,1802) = 3.000_r8*rxt(k,281)*y(k,44) + rxt(k,283)*y(k,46) + rxt(k,232) &
                      *y(k,59) + rxt(k,286)*y(k,82) + rxt(k,287)*y(k,83) + rxt(k,288) &
                      *y(k,84) + rxt(k,235)*y(k,85)
         mat(k,1918) = rxt(k,227)*y(k,60)
         mat(k,1579) = 2.000_r8*rxt(k,223)*y(k,59)
         mat(k,951) = rxt(k,227)*y(k,56) + (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,85)
         mat(k,1488) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,60) + (rxt(k,563) &
                       +rxt(k,569)+rxt(k,574))*y(k,93)
         mat(k,834) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,85)
         mat(k,1578) = 2.000_r8*rxt(k,248)*y(k,59)
         mat(k,1590) = -(rxt(k,221)*y(k,199) + (4._r8*rxt(k,222) + 4._r8*rxt(k,223) &
                      + 4._r8*rxt(k,224) + 4._r8*rxt(k,248)) * y(k,59) + rxt(k,225) &
                      *y(k,90) + rxt(k,226)*y(k,124) + rxt(k,228)*y(k,125) + rxt(k,231) &
                      *y(k,133) + (rxt(k,232) + rxt(k,233)) * y(k,226) + (rxt(k,254) &
                      + rxt(k,255) + rxt(k,256)) * y(k,19) + rxt(k,510)*y(k,151))
         mat(k,2247) = -rxt(k,221)*y(k,59)
         mat(k,2195) = -rxt(k,225)*y(k,59)
         mat(k,2066) = -rxt(k,226)*y(k,59)
         mat(k,1842) = -rxt(k,228)*y(k,59)
         mat(k,2315) = -rxt(k,231)*y(k,59)
         mat(k,1797) = -(rxt(k,232) + rxt(k,233)) * y(k,59)
         mat(k,1564) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,59)
         mat(k,1478) = -rxt(k,510)*y(k,59)
         mat(k,1948) = rxt(k,219)*y(k,90) + rxt(k,236)*y(k,93) + rxt(k,220)*y(k,135)
         mat(k,955) = rxt(k,229)*y(k,133)
         mat(k,1496) = rxt(k,247)*y(k,222)
         mat(k,2195) = mat(k,2195) + rxt(k,219)*y(k,56)
         mat(k,837) = rxt(k,236)*y(k,56) + rxt(k,237)*y(k,133) + rxt(k,238)*y(k,226)
         mat(k,2315) = mat(k,2315) + rxt(k,229)*y(k,60) + rxt(k,237)*y(k,93)
         mat(k,2405) = rxt(k,220)*y(k,56)
         mat(k,335) = rxt(k,515)*y(k,151)
         mat(k,1478) = mat(k,1478) + rxt(k,515)*y(k,137)
         mat(k,1632) = rxt(k,247)*y(k,85)
         mat(k,1797) = mat(k,1797) + rxt(k,238)*y(k,93)
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
         mat(k,953) = -(rxt(k,227)*y(k,56) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,226) &
                      + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,85))
         mat(k,1933) = -rxt(k,227)*y(k,60)
         mat(k,2305) = -rxt(k,229)*y(k,60)
         mat(k,1760) = -rxt(k,230)*y(k,60)
         mat(k,1492) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,60)
         mat(k,1584) = rxt(k,228)*y(k,125)
         mat(k,1828) = rxt(k,228)*y(k,59)
         mat(k,1158) = -(rxt(k,314)*y(k,226))
         mat(k,1774) = -rxt(k,314)*y(k,62)
         mat(k,996) = .230_r8*rxt(k,480)*y(k,135)
         mat(k,1507) = rxt(k,250)*y(k,42)
         mat(k,293) = .350_r8*rxt(k,316)*y(k,226)
         mat(k,551) = .630_r8*rxt(k,318)*y(k,135)
         mat(k,1133) = .560_r8*rxt(k,347)*y(k,135)
         mat(k,2334) = rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) + rxt(k,295)*y(k,126) &
                      + rxt(k,296)*y(k,133) + rxt(k,297)*y(k,226)
         mat(k,370) = rxt(k,282)*y(k,56)
         mat(k,1278) = rxt(k,353)*y(k,126) + rxt(k,354)*y(k,226)
         mat(k,1937) = rxt(k,214)*y(k,42) + rxt(k,282)*y(k,46)
         mat(k,1447) = rxt(k,591)*y(k,227)
         mat(k,1053) = rxt(k,341)*y(k,226)
         mat(k,2174) = .070_r8*rxt(k,450)*y(k,200) + .160_r8*rxt(k,453)*y(k,212) &
                      + .140_r8*rxt(k,456)*y(k,214)
         mat(k,883) = .620_r8*rxt(k,425)*y(k,135)
         mat(k,1266) = .650_r8*rxt(k,378)*y(k,135)
         mat(k,1024) = .230_r8*rxt(k,483)*y(k,135)
         mat(k,1370) = .560_r8*rxt(k,392)*y(k,135)
         mat(k,2046) = .170_r8*rxt(k,451)*y(k,200) + .220_r8*rxt(k,376)*y(k,211) &
                      + .400_r8*rxt(k,454)*y(k,212) + .350_r8*rxt(k,457)*y(k,214) &
                      + .225_r8*rxt(k,492)*y(k,231) + .250_r8*rxt(k,433)*y(k,235)
         mat(k,1881) = rxt(k,295)*y(k,42) + rxt(k,353)*y(k,49) + .220_r8*rxt(k,375) &
                      *y(k,211) + .500_r8*rxt(k,434)*y(k,235)
         mat(k,2307) = rxt(k,296)*y(k,42) + rxt(k,504)*y(k,138)
         mat(k,2386) = .230_r8*rxt(k,480)*y(k,6) + .630_r8*rxt(k,318)*y(k,25) &
                      + .560_r8*rxt(k,347)*y(k,29) + .620_r8*rxt(k,425)*y(k,99) &
                      + .650_r8*rxt(k,378)*y(k,105) + .230_r8*rxt(k,483)*y(k,110) &
                      + .560_r8*rxt(k,392)*y(k,111)
         mat(k,364) = rxt(k,504)*y(k,133) + rxt(k,505)*y(k,226)
         mat(k,1102) = .700_r8*rxt(k,501)*y(k,226)
         mat(k,1414) = .220_r8*rxt(k,372)*y(k,211) + .250_r8*rxt(k,430)*y(k,235)
         mat(k,2229) = .110_r8*rxt(k,373)*y(k,211) + .125_r8*rxt(k,490)*y(k,231) &
                      + .200_r8*rxt(k,431)*y(k,235)
         mat(k,777) = .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451)*y(k,124)
         mat(k,1348) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,198) + .110_r8*rxt(k,373)*y(k,199)
         mat(k,740) = .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454)*y(k,124)
         mat(k,915) = .140_r8*rxt(k,456)*y(k,90) + .350_r8*rxt(k,457)*y(k,124)
         mat(k,1774) = mat(k,1774) + .350_r8*rxt(k,316)*y(k,24) + rxt(k,297)*y(k,42) &
                      + rxt(k,354)*y(k,49) + rxt(k,341)*y(k,75) + rxt(k,505)*y(k,138) &
                      + .700_r8*rxt(k,501)*y(k,179)
         mat(k,807) = rxt(k,591)*y(k,63)
         mat(k,1174) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,199)
         mat(k,1228) = .250_r8*rxt(k,433)*y(k,124) + .500_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .200_r8*rxt(k,431)*y(k,199)
         mat(k,1448) = -(rxt(k,591)*y(k,227))
         mat(k,808) = -rxt(k,591)*y(k,63)
         mat(k,1000) = .270_r8*rxt(k,480)*y(k,135)
         mat(k,1137) = .200_r8*rxt(k,347)*y(k,135)
         mat(k,704) = rxt(k,334)*y(k,226)
         mat(k,620) = .500_r8*rxt(k,335)*y(k,226)
         mat(k,1159) = rxt(k,314)*y(k,226)
         mat(k,1165) = .800_r8*rxt(k,340)*y(k,226)
         mat(k,1054) = rxt(k,341)*y(k,226)
         mat(k,909) = rxt(k,306)*y(k,226)
         mat(k,2189) = .450_r8*rxt(k,383)*y(k,213)
         mat(k,583) = .500_r8*rxt(k,391)*y(k,226)
         mat(k,1028) = .270_r8*rxt(k,483)*y(k,135)
         mat(k,1377) = .100_r8*rxt(k,392)*y(k,135)
         mat(k,2062) = rxt(k,333)*y(k,198) + .900_r8*rxt(k,492)*y(k,231)
         mat(k,2400) = .270_r8*rxt(k,480)*y(k,6) + .200_r8*rxt(k,347)*y(k,29) &
                      + .270_r8*rxt(k,483)*y(k,110) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,1105) = 1.800_r8*rxt(k,501)*y(k,226)
         mat(k,1427) = rxt(k,333)*y(k,124) + 4.000_r8*rxt(k,330)*y(k,198) &
                      + .900_r8*rxt(k,331)*y(k,199) + rxt(k,405)*y(k,206) &
                      + 2.000_r8*rxt(k,381)*y(k,213) + rxt(k,430)*y(k,235)
         mat(k,2244) = .900_r8*rxt(k,331)*y(k,198) + rxt(k,382)*y(k,213) &
                      + .500_r8*rxt(k,490)*y(k,231)
         mat(k,1302) = rxt(k,405)*y(k,198)
         mat(k,1397) = .450_r8*rxt(k,383)*y(k,90) + 2.000_r8*rxt(k,381)*y(k,198) &
                      + rxt(k,382)*y(k,199) + 4.000_r8*rxt(k,384)*y(k,213)
         mat(k,1790) = rxt(k,334)*y(k,50) + .500_r8*rxt(k,335)*y(k,51) + rxt(k,314) &
                      *y(k,62) + .800_r8*rxt(k,340)*y(k,74) + rxt(k,341)*y(k,75) &
                      + rxt(k,306)*y(k,87) + .500_r8*rxt(k,391)*y(k,109) &
                      + 1.800_r8*rxt(k,501)*y(k,179)
         mat(k,1179) = .900_r8*rxt(k,492)*y(k,124) + .500_r8*rxt(k,490)*y(k,199)
         mat(k,1234) = rxt(k,430)*y(k,198)
         mat(k,242) = -(rxt(k,275)*y(k,222))
         mat(k,1618) = -rxt(k,275)*y(k,64)
         mat(k,138) = rxt(k,240)*y(k,222)
         mat(k,143) = rxt(k,266)*y(k,222)
         mat(k,149) = rxt(k,242)*y(k,222)
         mat(k,114) = 2.000_r8*rxt(k,243)*y(k,222)
         mat(k,153) = 2.000_r8*rxt(k,244)*y(k,222)
         mat(k,118) = rxt(k,245)*y(k,222)
         mat(k,106) = 2.000_r8*rxt(k,268)*y(k,222)
         mat(k,248) = rxt(k,292)*y(k,222) + rxt(k,287)*y(k,226)
         mat(k,271) = rxt(k,293)*y(k,222) + rxt(k,288)*y(k,226)
         mat(k,1618) = mat(k,1618) + rxt(k,240)*y(k,34) + rxt(k,266)*y(k,35) &
                      + rxt(k,242)*y(k,37) + 2.000_r8*rxt(k,243)*y(k,38) &
                      + 2.000_r8*rxt(k,244)*y(k,39) + rxt(k,245)*y(k,40) &
                      + 2.000_r8*rxt(k,268)*y(k,78) + rxt(k,292)*y(k,83) + rxt(k,293) &
                      *y(k,84)
         mat(k,1681) = rxt(k,287)*y(k,83) + rxt(k,288)*y(k,84)
         mat(k,224) = -(rxt(k,276)*y(k,222))
         mat(k,1617) = -rxt(k,276)*y(k,65)
         mat(k,110) = rxt(k,241)*y(k,222)
         mat(k,148) = rxt(k,242)*y(k,222)
         mat(k,220) = rxt(k,291)*y(k,222) + rxt(k,286)*y(k,226)
         mat(k,1617) = mat(k,1617) + rxt(k,241)*y(k,36) + rxt(k,242)*y(k,37) &
                      + rxt(k,291)*y(k,82)
         mat(k,1677) = rxt(k,286)*y(k,82)
         mat(k,192) = -(rxt(k,449)*y(k,226))
         mat(k,1671) = -rxt(k,449)*y(k,66)
         mat(k,186) = .180_r8*rxt(k,469)*y(k,226)
         mat(k,1671) = mat(k,1671) + .180_r8*rxt(k,469)*y(k,181)
         mat(k,308) = -(rxt(k,502)*y(k,126) + (rxt(k,503) + rxt(k,517)) * y(k,226))
         mat(k,1861) = -rxt(k,502)*y(k,67)
         mat(k,1691) = -(rxt(k,503) + rxt(k,517)) * y(k,67)
         mat(k,2115) = rxt(k,336)*y(k,204)
         mat(k,793) = rxt(k,336)*y(k,90)
         mat(k,923) = -(rxt(k,271)*y(k,54) + rxt(k,272)*y(k,77) + rxt(k,273)*y(k,239) &
                      + rxt(k,274)*y(k,89))
         mat(k,2264) = -rxt(k,271)*y(k,73)
         mat(k,1458) = -rxt(k,272)*y(k,73)
         mat(k,2424) = -rxt(k,273)*y(k,73)
         mat(k,1965) = -rxt(k,274)*y(k,73)
         mat(k,144) = rxt(k,266)*y(k,222)
         mat(k,154) = rxt(k,244)*y(k,222)
         mat(k,243) = 2.000_r8*rxt(k,275)*y(k,222)
         mat(k,225) = rxt(k,276)*y(k,222)
         mat(k,1626) = rxt(k,266)*y(k,35) + rxt(k,244)*y(k,39) + 2.000_r8*rxt(k,275) &
                      *y(k,64) + rxt(k,276)*y(k,65)
         mat(k,1164) = -(rxt(k,340)*y(k,226))
         mat(k,1775) = -rxt(k,340)*y(k,74)
         mat(k,590) = .700_r8*rxt(k,416)*y(k,226)
         mat(k,558) = .500_r8*rxt(k,417)*y(k,226)
         mat(k,379) = rxt(k,428)*y(k,226)
         mat(k,2047) = .050_r8*rxt(k,414)*y(k,207) + .530_r8*rxt(k,376)*y(k,211) &
                      + .225_r8*rxt(k,492)*y(k,231) + .250_r8*rxt(k,433)*y(k,235)
         mat(k,1882) = .050_r8*rxt(k,415)*y(k,207) + .530_r8*rxt(k,375)*y(k,211) &
                      + .250_r8*rxt(k,434)*y(k,235)
         mat(k,1536) = rxt(k,339)*y(k,203)
         mat(k,1415) = .530_r8*rxt(k,372)*y(k,211) + .250_r8*rxt(k,430)*y(k,235)
         mat(k,2230) = .260_r8*rxt(k,373)*y(k,211) + .125_r8*rxt(k,490)*y(k,231) &
                      + .100_r8*rxt(k,431)*y(k,235)
         mat(k,462) = rxt(k,339)*y(k,134)
         mat(k,1323) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1349) = .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375)*y(k,126) &
                      + .530_r8*rxt(k,372)*y(k,198) + .260_r8*rxt(k,373)*y(k,199)
         mat(k,1775) = mat(k,1775) + .700_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101) + rxt(k,428)*y(k,115)
         mat(k,1175) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,199)
         mat(k,1229) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,1052) = -(rxt(k,341)*y(k,226))
         mat(k,1765) = -rxt(k,341)*y(k,75)
         mat(k,292) = .650_r8*rxt(k,316)*y(k,226)
         mat(k,1162) = .200_r8*rxt(k,340)*y(k,226)
         mat(k,2167) = .160_r8*rxt(k,453)*y(k,212) + .070_r8*rxt(k,456)*y(k,214)
         mat(k,1114) = rxt(k,429)*y(k,226)
         mat(k,2038) = rxt(k,440)*y(k,192) + .050_r8*rxt(k,414)*y(k,207) &
                      + .400_r8*rxt(k,454)*y(k,212) + .170_r8*rxt(k,457)*y(k,214) &
                      + .700_r8*rxt(k,460)*y(k,228) + .600_r8*rxt(k,467)*y(k,233) &
                      + .250_r8*rxt(k,433)*y(k,235) + .340_r8*rxt(k,473)*y(k,236) &
                      + .170_r8*rxt(k,476)*y(k,238)
         mat(k,1872) = .050_r8*rxt(k,415)*y(k,207) + .250_r8*rxt(k,434)*y(k,235)
         mat(k,488) = rxt(k,440)*y(k,124)
         mat(k,1412) = .250_r8*rxt(k,430)*y(k,235)
         mat(k,2221) = .100_r8*rxt(k,431)*y(k,235)
         mat(k,1321) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,739) = .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454)*y(k,124)
         mat(k,914) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,1765) = mat(k,1765) + .650_r8*rxt(k,316)*y(k,24) + .200_r8*rxt(k,340) &
                      *y(k,74) + rxt(k,429)*y(k,116)
         mat(k,453) = .700_r8*rxt(k,460)*y(k,124)
         mat(k,752) = .600_r8*rxt(k,467)*y(k,124)
         mat(k,1227) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,768) = .340_r8*rxt(k,473)*y(k,124)
         mat(k,509) = .170_r8*rxt(k,476)*y(k,124)
         mat(k,2095) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,90) + rxt(k,175) &
                      *y(k,134) + rxt(k,178)*y(k,135))
         mat(k,2203) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76)
         mat(k,1548) = -rxt(k,175)*y(k,76)
         mat(k,2413) = -rxt(k,178)*y(k,76)
         mat(k,2349) = rxt(k,297)*y(k,226)
         mat(k,2280) = rxt(k,311)*y(k,222)
         mat(k,1956) = rxt(k,216)*y(k,77)
         mat(k,928) = rxt(k,272)*y(k,77)
         mat(k,1466) = rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73) + rxt(k,170)*y(k,133) &
                      + rxt(k,153)*y(k,222) + rxt(k,179)*y(k,226)
         mat(k,831) = rxt(k,270)*y(k,222)
         mat(k,1503) = rxt(k,247)*y(k,222)
         mat(k,977) = rxt(k,202)*y(k,226)
         mat(k,2323) = rxt(k,170)*y(k,77) + rxt(k,182)*y(k,226)
         mat(k,367) = rxt(k,505)*y(k,226)
         mat(k,712) = rxt(k,511)*y(k,226)
         mat(k,1483) = rxt(k,516)*y(k,226)
         mat(k,1640) = rxt(k,311)*y(k,54) + rxt(k,153)*y(k,77) + rxt(k,270)*y(k,81) &
                      + rxt(k,247)*y(k,85)
         mat(k,1805) = rxt(k,297)*y(k,42) + rxt(k,179)*y(k,77) + rxt(k,202)*y(k,112) &
                      + rxt(k,182)*y(k,133) + rxt(k,505)*y(k,138) + rxt(k,511) &
                      *y(k,149) + rxt(k,516)*y(k,151)
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
         mat(k,1459) = -(rxt(k,153)*y(k,222) + rxt(k,170)*y(k,133) + rxt(k,179) &
                      *y(k,226) + rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73))
         mat(k,1627) = -rxt(k,153)*y(k,77)
         mat(k,2309) = -rxt(k,170)*y(k,77)
         mat(k,1791) = -rxt(k,179)*y(k,77)
         mat(k,1943) = -rxt(k,216)*y(k,77)
         mat(k,924) = -rxt(k,272)*y(k,77)
         mat(k,2267) = rxt(k,312)*y(k,222)
         mat(k,2082) = rxt(k,172)*y(k,90)
         mat(k,2190) = rxt(k,172)*y(k,76)
         mat(k,1627) = mat(k,1627) + rxt(k,312)*y(k,54)
         mat(k,105) = -(rxt(k,268)*y(k,222))
         mat(k,1606) = -rxt(k,268)*y(k,78)
         mat(k,611) = -(rxt(k,171)*y(k,133) + rxt(k,180)*y(k,226) + rxt(k,217)*y(k,56))
         mat(k,2294) = -rxt(k,171)*y(k,79)
         mat(k,1730) = -rxt(k,180)*y(k,79)
         mat(k,1927) = -rxt(k,217)*y(k,79)
         mat(k,2142) = 2.000_r8*rxt(k,186)*y(k,90)
         mat(k,1730) = mat(k,1730) + 2.000_r8*rxt(k,185)*y(k,226)
         mat(k,258) = rxt(k,518)*y(k,239)
         mat(k,2421) = rxt(k,518)*y(k,153)
         mat(k,826) = -(rxt(k,263)*y(k,133) + rxt(k,264)*y(k,226) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,222))
         mat(k,2299) = -rxt(k,263)*y(k,81)
         mat(k,1750) = -rxt(k,264)*y(k,81)
         mat(k,1625) = -(rxt(k,269) + rxt(k,270)) * y(k,81)
         mat(k,1506) = rxt(k,250)*y(k,42) + rxt(k,251)*y(k,90)
         mat(k,2332) = rxt(k,250)*y(k,17)
         mat(k,2159) = rxt(k,251)*y(k,17)
         mat(k,219) = -(rxt(k,286)*y(k,226) + rxt(k,291)*y(k,222))
         mat(k,1676) = -rxt(k,286)*y(k,82)
         mat(k,1616) = -rxt(k,291)*y(k,82)
         mat(k,249) = -(rxt(k,287)*y(k,226) + rxt(k,292)*y(k,222))
         mat(k,1683) = -rxt(k,287)*y(k,83)
         mat(k,1619) = -rxt(k,292)*y(k,83)
         mat(k,272) = -(rxt(k,288)*y(k,226) + rxt(k,293)*y(k,222))
         mat(k,1687) = -rxt(k,288)*y(k,84)
         mat(k,1621) = -rxt(k,293)*y(k,84)
         mat(k,1493) = -(rxt(k,234)*y(k,133) + rxt(k,235)*y(k,226) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,222) + (rxt(k,563) + rxt(k,569) + rxt(k,574) &
                      ) * y(k,93) + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,60) &
                      + (rxt(k,570) + rxt(k,575)) * y(k,92))
         mat(k,2311) = -rxt(k,234)*y(k,85)
         mat(k,1793) = -rxt(k,235)*y(k,85)
         mat(k,1628) = -(rxt(k,246) + rxt(k,247)) * y(k,85)
         mat(k,836) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,85)
         mat(k,954) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,85)
         mat(k,786) = -(rxt(k,570) + rxt(k,575)) * y(k,85)
         mat(k,298) = rxt(k,325)*y(k,56)
         mat(k,304) = rxt(k,355)*y(k,56)
         mat(k,471) = rxt(k,277)*y(k,56)
         mat(k,2337) = rxt(k,214)*y(k,56)
         mat(k,603) = rxt(k,279)*y(k,56)
         mat(k,371) = 2.000_r8*rxt(k,282)*y(k,56)
         mat(k,2268) = rxt(k,215)*y(k,56)
         mat(k,446) = rxt(k,284)*y(k,56)
         mat(k,1944) = rxt(k,325)*y(k,28) + rxt(k,355)*y(k,31) + rxt(k,277)*y(k,41) &
                      + rxt(k,214)*y(k,42) + rxt(k,279)*y(k,43) + 2.000_r8*rxt(k,282) &
                      *y(k,46) + rxt(k,215)*y(k,54) + rxt(k,284)*y(k,55) + rxt(k,216) &
                      *y(k,77) + rxt(k,217)*y(k,79) + rxt(k,218)*y(k,90) + rxt(k,236) &
                      *y(k,93)
         mat(k,1586) = rxt(k,233)*y(k,226)
         mat(k,1460) = rxt(k,216)*y(k,56)
         mat(k,612) = rxt(k,217)*y(k,56)
         mat(k,2191) = rxt(k,218)*y(k,56)
         mat(k,836) = mat(k,836) + rxt(k,236)*y(k,56)
         mat(k,1793) = mat(k,1793) + rxt(k,233)*y(k,59)
         mat(k,180) = -(rxt(k,305)*y(k,226) + rxt(k,313)*y(k,222))
         mat(k,1669) = -rxt(k,305)*y(k,86)
         mat(k,1615) = -rxt(k,313)*y(k,86)
         mat(k,908) = -(rxt(k,306)*y(k,226))
         mat(k,1755) = -rxt(k,306)*y(k,87)
         mat(k,987) = .050_r8*rxt(k,480)*y(k,135)
         mat(k,291) = .350_r8*rxt(k,316)*y(k,226)
         mat(k,550) = .370_r8*rxt(k,318)*y(k,135)
         mat(k,1130) = .120_r8*rxt(k,347)*y(k,135)
         mat(k,2163) = rxt(k,307)*y(k,205)
         mat(k,881) = .110_r8*rxt(k,425)*y(k,135)
         mat(k,1265) = .330_r8*rxt(k,378)*y(k,135)
         mat(k,1015) = .050_r8*rxt(k,483)*y(k,135)
         mat(k,1367) = .120_r8*rxt(k,392)*y(k,135)
         mat(k,2033) = rxt(k,309)*y(k,205)
         mat(k,2373) = .050_r8*rxt(k,480)*y(k,6) + .370_r8*rxt(k,318)*y(k,25) &
                      + .120_r8*rxt(k,347)*y(k,29) + .110_r8*rxt(k,425)*y(k,99) &
                      + .330_r8*rxt(k,378)*y(k,105) + .050_r8*rxt(k,483)*y(k,110) &
                      + .120_r8*rxt(k,392)*y(k,111)
         mat(k,440) = rxt(k,307)*y(k,90) + rxt(k,309)*y(k,124)
         mat(k,1755) = mat(k,1755) + .350_r8*rxt(k,316)*y(k,24)
         mat(k,2263) = rxt(k,271)*y(k,73)
         mat(k,922) = rxt(k,271)*y(k,54) + rxt(k,272)*y(k,77) + rxt(k,274)*y(k,89) &
                      + rxt(k,273)*y(k,239)
         mat(k,1457) = rxt(k,272)*y(k,73)
         mat(k,1964) = rxt(k,274)*y(k,73)
         mat(k,2423) = rxt(k,273)*y(k,73)
         mat(k,1977) = -(rxt(k,211)*y(k,226) + rxt(k,274)*y(k,73))
         mat(k,1803) = -rxt(k,211)*y(k,89)
         mat(k,927) = -rxt(k,274)*y(k,89)
         mat(k,2347) = rxt(k,295)*y(k,126)
         mat(k,1154) = rxt(k,327)*y(k,126)
         mat(k,1284) = rxt(k,353)*y(k,126)
         mat(k,960) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,85)
         mat(k,312) = rxt(k,502)*y(k,126)
         mat(k,1502) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,60)
         mat(k,1848) = rxt(k,210)*y(k,226)
         mat(k,1908) = rxt(k,295)*y(k,42) + rxt(k,327)*y(k,45) + rxt(k,353)*y(k,49) &
                      + rxt(k,502)*y(k,67)
         mat(k,1803) = mat(k,1803) + rxt(k,210)*y(k,125)
         mat(k,2204) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76) + rxt(k,176) &
                      *y(k,133) + rxt(k,177)*y(k,135) + rxt(k,181)*y(k,226) &
                      + 4._r8*rxt(k,186)*y(k,90) + rxt(k,198)*y(k,126) + rxt(k,203) &
                      *y(k,124) + rxt(k,208)*y(k,125) + (rxt(k,218) + rxt(k,219) &
                      ) * y(k,56) + rxt(k,225)*y(k,59) + rxt(k,251)*y(k,17) + rxt(k,257) &
                      *y(k,19) + rxt(k,294)*y(k,42) + rxt(k,300)*y(k,199) + rxt(k,307) &
                      *y(k,205) + rxt(k,321)*y(k,195) + rxt(k,332)*y(k,198) + rxt(k,336) &
                      *y(k,204) + rxt(k,349)*y(k,196) + rxt(k,358)*y(k,229) + rxt(k,362) &
                      *y(k,230) + rxt(k,374)*y(k,211) + rxt(k,383)*y(k,213) + rxt(k,387) &
                      *y(k,215) + rxt(k,397)*y(k,190) + rxt(k,407)*y(k,206) + rxt(k,412) &
                      *y(k,207) + rxt(k,421)*y(k,208) + rxt(k,432)*y(k,235) + rxt(k,436) &
                      *y(k,189) + rxt(k,439)*y(k,192) + rxt(k,443)*y(k,194) + rxt(k,446) &
                      *y(k,197) + rxt(k,450)*y(k,200) + rxt(k,453)*y(k,212) + rxt(k,456) &
                      *y(k,214) + rxt(k,459)*y(k,228) + rxt(k,466)*y(k,233) + rxt(k,472) &
                      *y(k,236) + rxt(k,475)*y(k,238) + rxt(k,486)*y(k,221) + rxt(k,491) &
                      *y(k,231) + rxt(k,496)*y(k,232))
         mat(k,2096) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,90)
         mat(k,2324) = -rxt(k,176)*y(k,90)
         mat(k,2414) = -rxt(k,177)*y(k,90)
         mat(k,1806) = -rxt(k,181)*y(k,90)
         mat(k,1911) = -rxt(k,198)*y(k,90)
         mat(k,2075) = -rxt(k,203)*y(k,90)
         mat(k,1851) = -rxt(k,208)*y(k,90)
         mat(k,1957) = -(rxt(k,218) + rxt(k,219)) * y(k,90)
         mat(k,1599) = -rxt(k,225)*y(k,90)
         mat(k,1517) = -rxt(k,251)*y(k,90)
         mat(k,1573) = -rxt(k,257)*y(k,90)
         mat(k,2350) = -rxt(k,294)*y(k,90)
         mat(k,2256) = -rxt(k,300)*y(k,90)
         mat(k,443) = -rxt(k,307)*y(k,90)
         mat(k,905) = -rxt(k,321)*y(k,90)
         mat(k,1434) = -rxt(k,332)*y(k,90)
         mat(k,801) = -rxt(k,336)*y(k,90)
         mat(k,940) = -rxt(k,349)*y(k,90)
         mat(k,823) = -rxt(k,358)*y(k,90)
         mat(k,1221) = -rxt(k,362)*y(k,90)
         mat(k,1363) = -rxt(k,374)*y(k,90)
         mat(k,1403) = -rxt(k,383)*y(k,90)
         mat(k,702) = -rxt(k,387)*y(k,90)
         mat(k,1048) = -rxt(k,397)*y(k,90)
         mat(k,1308) = -rxt(k,407)*y(k,90)
         mat(k,1341) = -rxt(k,412)*y(k,90)
         mat(k,1262) = -rxt(k,421)*y(k,90)
         mat(k,1239) = -rxt(k,432)*y(k,90)
         mat(k,526) = -rxt(k,436)*y(k,90)
         mat(k,492) = -rxt(k,439)*y(k,90)
         mat(k,438) = -rxt(k,443)*y(k,90)
         mat(k,637) = -rxt(k,446)*y(k,90)
         mat(k,782) = -rxt(k,450)*y(k,90)
         mat(k,743) = -rxt(k,453)*y(k,90)
         mat(k,920) = -rxt(k,456)*y(k,90)
         mat(k,457) = -rxt(k,459)*y(k,90)
         mat(k,758) = -rxt(k,466)*y(k,90)
         mat(k,775) = -rxt(k,472)*y(k,90)
         mat(k,514) = -rxt(k,475)*y(k,90)
         mat(k,1097) = -rxt(k,486)*y(k,90)
         mat(k,1184) = -rxt(k,491)*y(k,90)
         mat(k,1067) = -rxt(k,496)*y(k,90)
         mat(k,1003) = .570_r8*rxt(k,480)*y(k,135)
         mat(k,163) = .650_r8*rxt(k,438)*y(k,226)
         mat(k,1517) = mat(k,1517) + rxt(k,250)*y(k,42)
         mat(k,1573) = mat(k,1573) + rxt(k,262)*y(k,226)
         mat(k,295) = .350_r8*rxt(k,316)*y(k,226)
         mat(k,553) = .130_r8*rxt(k,318)*y(k,135)
         mat(k,266) = rxt(k,323)*y(k,226)
         mat(k,1143) = .280_r8*rxt(k,347)*y(k,135)
         mat(k,2350) = mat(k,2350) + rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,133)
         mat(k,608) = rxt(k,279)*y(k,56) + rxt(k,280)*y(k,226)
         mat(k,374) = rxt(k,282)*y(k,56) + rxt(k,283)*y(k,226)
         mat(k,104) = rxt(k,329)*y(k,226)
         mat(k,813) = rxt(k,302)*y(k,226)
         mat(k,2281) = rxt(k,311)*y(k,222)
         mat(k,1957) = mat(k,1957) + rxt(k,214)*y(k,42) + rxt(k,279)*y(k,43) &
                      + rxt(k,282)*y(k,46) + rxt(k,217)*y(k,79)
         mat(k,1599) = mat(k,1599) + rxt(k,221)*y(k,199) + rxt(k,232)*y(k,226)
         mat(k,1161) = rxt(k,314)*y(k,226)
         mat(k,196) = .730_r8*rxt(k,449)*y(k,226)
         mat(k,313) = .500_r8*rxt(k,517)*y(k,226)
         mat(k,1167) = rxt(k,340)*y(k,226)
         mat(k,1056) = rxt(k,341)*y(k,226)
         mat(k,2096) = mat(k,2096) + rxt(k,175)*y(k,134)
         mat(k,615) = rxt(k,217)*y(k,56) + rxt(k,171)*y(k,133) + rxt(k,180)*y(k,226)
         mat(k,183) = rxt(k,305)*y(k,226)
         mat(k,911) = rxt(k,306)*y(k,226)
         mat(k,2204) = mat(k,2204) + .070_r8*rxt(k,450)*y(k,200) + .160_r8*rxt(k,453) &
                      *y(k,212) + .330_r8*rxt(k,456)*y(k,214)
         mat(k,1202) = rxt(k,371)*y(k,226)
         mat(k,1210) = rxt(k,356)*y(k,226)
         mat(k,892) = .370_r8*rxt(k,425)*y(k,135)
         mat(k,596) = .300_r8*rxt(k,416)*y(k,226)
         mat(k,563) = rxt(k,417)*y(k,226)
         mat(k,424) = rxt(k,424)*y(k,226)
         mat(k,1274) = .140_r8*rxt(k,378)*y(k,135)
         mat(k,318) = .200_r8*rxt(k,380)*y(k,226)
         mat(k,587) = .500_r8*rxt(k,391)*y(k,226)
         mat(k,1031) = .570_r8*rxt(k,483)*y(k,135)
         mat(k,1384) = .280_r8*rxt(k,392)*y(k,135)
         mat(k,382) = rxt(k,428)*y(k,226)
         mat(k,1124) = rxt(k,429)*y(k,226)
         mat(k,2075) = mat(k,2075) + rxt(k,398)*y(k,190) + rxt(k,440)*y(k,192) &
                      + rxt(k,445)*y(k,194) + rxt(k,322)*y(k,195) + rxt(k,350) &
                      *y(k,196) + rxt(k,301)*y(k,199) + .170_r8*rxt(k,451)*y(k,200) &
                      + rxt(k,369)*y(k,202) + .250_r8*rxt(k,337)*y(k,204) + rxt(k,309) &
                      *y(k,205) + .920_r8*rxt(k,408)*y(k,206) + .920_r8*rxt(k,414) &
                      *y(k,207) + rxt(k,422)*y(k,208) + .470_r8*rxt(k,376)*y(k,211) &
                      + .400_r8*rxt(k,454)*y(k,212) + .830_r8*rxt(k,457)*y(k,214) &
                      + rxt(k,460)*y(k,228) + rxt(k,359)*y(k,229) + .900_r8*rxt(k,492) &
                      *y(k,231) + .800_r8*rxt(k,497)*y(k,232) + rxt(k,467)*y(k,233) &
                      + rxt(k,433)*y(k,235) + rxt(k,473)*y(k,236) + rxt(k,476) &
                      *y(k,238)
         mat(k,1911) = mat(k,1911) + rxt(k,295)*y(k,42) + rxt(k,409)*y(k,206) &
                      + rxt(k,415)*y(k,207) + rxt(k,423)*y(k,208) + .470_r8*rxt(k,375) &
                      *y(k,211) + rxt(k,201)*y(k,226) + rxt(k,434)*y(k,235)
         mat(k,2324) = mat(k,2324) + rxt(k,296)*y(k,42) + rxt(k,171)*y(k,79)
         mat(k,1549) = rxt(k,175)*y(k,76) + rxt(k,339)*y(k,203)
         mat(k,2414) = mat(k,2414) + .570_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,318) &
                      *y(k,25) + .280_r8*rxt(k,347)*y(k,29) + .370_r8*rxt(k,425) &
                      *y(k,99) + .140_r8*rxt(k,378)*y(k,105) + .570_r8*rxt(k,483) &
                      *y(k,110) + .280_r8*rxt(k,392)*y(k,111) + rxt(k,183)*y(k,226)
         mat(k,172) = .800_r8*rxt(k,461)*y(k,226)
         mat(k,948) = rxt(k,507)*y(k,226)
         mat(k,1107) = .200_r8*rxt(k,501)*y(k,226)
         mat(k,191) = .280_r8*rxt(k,469)*y(k,226)
         mat(k,213) = .380_r8*rxt(k,471)*y(k,226)
         mat(k,218) = .630_r8*rxt(k,477)*y(k,226)
         mat(k,1048) = mat(k,1048) + rxt(k,398)*y(k,124)
         mat(k,492) = mat(k,492) + rxt(k,440)*y(k,124)
         mat(k,438) = mat(k,438) + rxt(k,445)*y(k,124)
         mat(k,905) = mat(k,905) + rxt(k,322)*y(k,124) + 2.400_r8*rxt(k,319)*y(k,195) &
                      + rxt(k,320)*y(k,199)
         mat(k,940) = mat(k,940) + rxt(k,350)*y(k,124) + rxt(k,348)*y(k,199)
         mat(k,1434) = mat(k,1434) + .900_r8*rxt(k,331)*y(k,199) + rxt(k,405)*y(k,206) &
                      + rxt(k,410)*y(k,207) + rxt(k,419)*y(k,208) + .470_r8*rxt(k,372) &
                      *y(k,211) + rxt(k,430)*y(k,235)
         mat(k,2256) = mat(k,2256) + rxt(k,221)*y(k,59) + rxt(k,301)*y(k,124) &
                      + rxt(k,320)*y(k,195) + rxt(k,348)*y(k,196) + .900_r8*rxt(k,331) &
                      *y(k,198) + 4.000_r8*rxt(k,298)*y(k,199) + rxt(k,406)*y(k,206) &
                      + rxt(k,411)*y(k,207) + 1.200_r8*rxt(k,420)*y(k,208) &
                      + .730_r8*rxt(k,373)*y(k,211) + rxt(k,382)*y(k,213) &
                      + .500_r8*rxt(k,485)*y(k,221) + .300_r8*rxt(k,361)*y(k,230) &
                      + rxt(k,490)*y(k,231) + rxt(k,495)*y(k,232) + .800_r8*rxt(k,431) &
                      *y(k,235)
         mat(k,782) = mat(k,782) + .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451) &
                      *y(k,124)
         mat(k,578) = rxt(k,369)*y(k,124)
         mat(k,464) = rxt(k,339)*y(k,134)
         mat(k,801) = mat(k,801) + .250_r8*rxt(k,337)*y(k,124)
         mat(k,443) = mat(k,443) + rxt(k,309)*y(k,124)
         mat(k,1308) = mat(k,1308) + .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) &
                      + rxt(k,405)*y(k,198) + rxt(k,406)*y(k,199)
         mat(k,1341) = mat(k,1341) + .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,198) + rxt(k,411)*y(k,199)
         mat(k,1262) = mat(k,1262) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) &
                      + rxt(k,419)*y(k,198) + 1.200_r8*rxt(k,420)*y(k,199)
         mat(k,1363) = mat(k,1363) + .470_r8*rxt(k,376)*y(k,124) + .470_r8*rxt(k,375) &
                      *y(k,126) + .470_r8*rxt(k,372)*y(k,198) + .730_r8*rxt(k,373) &
                      *y(k,199)
         mat(k,743) = mat(k,743) + .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454) &
                      *y(k,124)
         mat(k,1403) = mat(k,1403) + rxt(k,382)*y(k,199)
         mat(k,920) = mat(k,920) + .330_r8*rxt(k,456)*y(k,90) + .830_r8*rxt(k,457) &
                      *y(k,124)
         mat(k,1097) = mat(k,1097) + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1641) = rxt(k,311)*y(k,54)
         mat(k,1806) = mat(k,1806) + .650_r8*rxt(k,438)*y(k,7) + rxt(k,262)*y(k,19) &
                      + .350_r8*rxt(k,316)*y(k,24) + rxt(k,323)*y(k,26) + rxt(k,280) &
                      *y(k,43) + rxt(k,283)*y(k,46) + rxt(k,329)*y(k,47) + rxt(k,302) &
                      *y(k,52) + rxt(k,232)*y(k,59) + rxt(k,314)*y(k,62) &
                      + .730_r8*rxt(k,449)*y(k,66) + .500_r8*rxt(k,517)*y(k,67) &
                      + rxt(k,340)*y(k,74) + rxt(k,341)*y(k,75) + rxt(k,180)*y(k,79) &
                      + rxt(k,305)*y(k,86) + rxt(k,306)*y(k,87) + rxt(k,371)*y(k,94) &
                      + rxt(k,356)*y(k,96) + .300_r8*rxt(k,416)*y(k,100) + rxt(k,417) &
                      *y(k,101) + rxt(k,424)*y(k,102) + .200_r8*rxt(k,380)*y(k,106) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,428)*y(k,115) + rxt(k,429) &
                      *y(k,116) + rxt(k,201)*y(k,126) + rxt(k,183)*y(k,135) &
                      + .800_r8*rxt(k,461)*y(k,143) + rxt(k,507)*y(k,152) &
                      + .200_r8*rxt(k,501)*y(k,179) + .280_r8*rxt(k,469)*y(k,181) &
                      + .380_r8*rxt(k,471)*y(k,183) + .630_r8*rxt(k,477)*y(k,185)
         mat(k,457) = mat(k,457) + rxt(k,460)*y(k,124)
         mat(k,823) = mat(k,823) + rxt(k,359)*y(k,124)
         mat(k,1221) = mat(k,1221) + .300_r8*rxt(k,361)*y(k,199)
         mat(k,1184) = mat(k,1184) + .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,199)
         mat(k,1067) = mat(k,1067) + .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,199)
         mat(k,758) = mat(k,758) + rxt(k,467)*y(k,124)
         mat(k,1239) = mat(k,1239) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126) &
                      + rxt(k,430)*y(k,198) + .800_r8*rxt(k,431)*y(k,199)
         mat(k,775) = mat(k,775) + rxt(k,473)*y(k,124)
         mat(k,514) = mat(k,514) + rxt(k,476)*y(k,124)
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
         mat(k,493) = -(rxt(k,187)*y(k,226))
         mat(k,1717) = -rxt(k,187)*y(k,91)
         mat(k,2138) = rxt(k,208)*y(k,125)
         mat(k,1817) = rxt(k,208)*y(k,90)
         mat(k,785) = -(rxt(k,265)*y(k,133) + (rxt(k,570) + rxt(k,575)) * y(k,85))
         mat(k,2297) = -rxt(k,265)*y(k,92)
         mat(k,1490) = -(rxt(k,570) + rxt(k,575)) * y(k,92)
         mat(k,1557) = rxt(k,257)*y(k,90)
         mat(k,2156) = rxt(k,257)*y(k,19)
         mat(k,835) = -(rxt(k,236)*y(k,56) + rxt(k,237)*y(k,133) + rxt(k,238)*y(k,226) &
                      + (rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,85))
         mat(k,1929) = -rxt(k,236)*y(k,93)
         mat(k,2300) = -rxt(k,237)*y(k,93)
         mat(k,1751) = -rxt(k,238)*y(k,93)
         mat(k,1491) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,93)
         mat(k,1582) = rxt(k,225)*y(k,90)
         mat(k,952) = rxt(k,230)*y(k,226)
         mat(k,2160) = rxt(k,225)*y(k,59)
         mat(k,1751) = mat(k,1751) + rxt(k,230)*y(k,60)
         mat(k,1193) = -(rxt(k,371)*y(k,226))
         mat(k,1777) = -rxt(k,371)*y(k,94)
         mat(k,591) = .300_r8*rxt(k,416)*y(k,226)
         mat(k,559) = .500_r8*rxt(k,417)*y(k,226)
         mat(k,2049) = rxt(k,370)*y(k,202) + rxt(k,377)*y(k,211)
         mat(k,575) = rxt(k,370)*y(k,124)
         mat(k,1350) = rxt(k,377)*y(k,124)
         mat(k,1777) = mat(k,1777) + .300_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101)
         mat(k,228) = -(rxt(k,402)*y(k,226))
         mat(k,1678) = -rxt(k,402)*y(k,95)
         mat(k,1206) = -(rxt(k,356)*y(k,226))
         mat(k,1778) = -rxt(k,356)*y(k,96)
         mat(k,592) = .700_r8*rxt(k,416)*y(k,226)
         mat(k,560) = .500_r8*rxt(k,417)*y(k,226)
         mat(k,581) = .500_r8*rxt(k,391)*y(k,226)
         mat(k,2050) = .050_r8*rxt(k,414)*y(k,207) + .220_r8*rxt(k,376)*y(k,211) &
                      + .250_r8*rxt(k,433)*y(k,235)
         mat(k,1885) = .050_r8*rxt(k,415)*y(k,207) + .220_r8*rxt(k,375)*y(k,211) &
                      + .250_r8*rxt(k,434)*y(k,235)
         mat(k,543) = .500_r8*rxt(k,360)*y(k,226)
         mat(k,1416) = .220_r8*rxt(k,372)*y(k,211) + .250_r8*rxt(k,430)*y(k,235)
         mat(k,2232) = .230_r8*rxt(k,373)*y(k,211) + .200_r8*rxt(k,361)*y(k,230) &
                      + .100_r8*rxt(k,431)*y(k,235)
         mat(k,1325) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1351) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,198) + .230_r8*rxt(k,373)*y(k,199)
         mat(k,1778) = mat(k,1778) + .700_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101) + .500_r8*rxt(k,391)*y(k,109) + .500_r8*rxt(k,360) &
                      *y(k,147)
         mat(k,1214) = .200_r8*rxt(k,361)*y(k,199)
         mat(k,1230) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,324) = -(rxt(k,403)*y(k,226))
         mat(k,1694) = -rxt(k,403)*y(k,97)
         mat(k,2002) = .870_r8*rxt(k,414)*y(k,207)
         mat(k,1862) = .950_r8*rxt(k,415)*y(k,207)
         mat(k,1408) = rxt(k,410)*y(k,207)
         mat(k,2212) = .750_r8*rxt(k,411)*y(k,207)
         mat(k,1314) = .870_r8*rxt(k,414)*y(k,124) + .950_r8*rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,198) + .750_r8*rxt(k,411)*y(k,199)
         mat(k,131) = -(rxt(k,404)*y(k,226))
         mat(k,1665) = -rxt(k,404)*y(k,98)
         mat(k,730) = .600_r8*rxt(k,427)*y(k,226)
         mat(k,1665) = mat(k,1665) + .600_r8*rxt(k,427)*y(k,103)
         mat(k,880) = -(rxt(k,418)*y(k,126) + rxt(k,425)*y(k,135) + rxt(k,426) &
                      *y(k,226))
         mat(k,1866) = -rxt(k,418)*y(k,99)
         mat(k,2372) = -rxt(k,425)*y(k,99)
         mat(k,1753) = -rxt(k,426)*y(k,99)
         mat(k,589) = -(rxt(k,416)*y(k,226))
         mat(k,1727) = -rxt(k,416)*y(k,100)
         mat(k,2016) = .080_r8*rxt(k,408)*y(k,206)
         mat(k,1287) = .080_r8*rxt(k,408)*y(k,124)
         mat(k,556) = -(rxt(k,417)*y(k,226))
         mat(k,1724) = -rxt(k,417)*y(k,101)
         mat(k,2014) = .080_r8*rxt(k,414)*y(k,207)
         mat(k,1315) = .080_r8*rxt(k,414)*y(k,124)
         mat(k,419) = -(rxt(k,424)*y(k,226))
         mat(k,1707) = -rxt(k,424)*y(k,102)
         mat(k,2128) = rxt(k,421)*y(k,208)
         mat(k,1243) = rxt(k,421)*y(k,90)
         mat(k,731) = -(rxt(k,427)*y(k,226))
         mat(k,1742) = -rxt(k,427)*y(k,103)
         mat(k,2151) = rxt(k,407)*y(k,206) + rxt(k,412)*y(k,207)
         mat(k,1288) = rxt(k,407)*y(k,90)
         mat(k,1317) = rxt(k,412)*y(k,90)
         mat(k,74) = -(rxt(k,549)*y(k,226))
         mat(k,1658) = -rxt(k,549)*y(k,104)
         mat(k,1267) = -(rxt(k,378)*y(k,135) + rxt(k,379)*y(k,226))
         mat(k,2392) = -rxt(k,378)*y(k,105)
         mat(k,1782) = -rxt(k,379)*y(k,105)
         mat(k,885) = .300_r8*rxt(k,425)*y(k,135)
         mat(k,2054) = .360_r8*rxt(k,408)*y(k,206)
         mat(k,1889) = .400_r8*rxt(k,409)*y(k,206)
         mat(k,2392) = mat(k,2392) + .300_r8*rxt(k,425)*y(k,99)
         mat(k,1419) = .390_r8*rxt(k,405)*y(k,206)
         mat(k,2236) = .310_r8*rxt(k,406)*y(k,206)
         mat(k,1295) = .360_r8*rxt(k,408)*y(k,124) + .400_r8*rxt(k,409)*y(k,126) &
                      + .390_r8*rxt(k,405)*y(k,198) + .310_r8*rxt(k,406)*y(k,199)
         mat(k,314) = -(rxt(k,380)*y(k,226))
         mat(k,1692) = -rxt(k,380)*y(k,106)
         mat(k,2118) = rxt(k,374)*y(k,211)
         mat(k,1346) = rxt(k,374)*y(k,90)
         mat(k,515) = -(rxt(k,389)*y(k,226))
         mat(k,1719) = -rxt(k,389)*y(k,107)
         mat(k,2012) = .800_r8*rxt(k,398)*y(k,190)
         mat(k,1035) = .800_r8*rxt(k,398)*y(k,124)
         mat(k,319) = -(rxt(k,390)*y(k,226))
         mat(k,1693) = -rxt(k,390)*y(k,108)
         mat(k,2119) = .800_r8*rxt(k,387)*y(k,215)
         mat(k,695) = .800_r8*rxt(k,387)*y(k,90)
         mat(k,580) = -(rxt(k,391)*y(k,226))
         mat(k,1726) = -rxt(k,391)*y(k,109)
         mat(k,1821) = rxt(k,394)*y(k,213)
         mat(k,1390) = rxt(k,394)*y(k,125)
         mat(k,1016) = -(rxt(k,482)*y(k,126) + rxt(k,483)*y(k,135) + rxt(k,484) &
                      *y(k,226))
         mat(k,1870) = -rxt(k,482)*y(k,110)
         mat(k,2376) = -rxt(k,483)*y(k,110)
         mat(k,1763) = -rxt(k,484)*y(k,110)
         mat(k,1374) = -(rxt(k,392)*y(k,135) + rxt(k,393)*y(k,226))
         mat(k,2397) = -rxt(k,392)*y(k,111)
         mat(k,1787) = -rxt(k,393)*y(k,111)
         mat(k,888) = .200_r8*rxt(k,425)*y(k,135)
         mat(k,2059) = .560_r8*rxt(k,408)*y(k,206)
         mat(k,1894) = .600_r8*rxt(k,409)*y(k,206)
         mat(k,2397) = mat(k,2397) + .200_r8*rxt(k,425)*y(k,99)
         mat(k,1424) = .610_r8*rxt(k,405)*y(k,206)
         mat(k,2241) = .440_r8*rxt(k,406)*y(k,206)
         mat(k,1299) = .560_r8*rxt(k,408)*y(k,124) + .600_r8*rxt(k,409)*y(k,126) &
                      + .610_r8*rxt(k,405)*y(k,198) + .440_r8*rxt(k,406)*y(k,199)
         mat(k,969) = -(rxt(k,190)*y(k,124) + (rxt(k,191) + rxt(k,192) + rxt(k,193) &
                      ) * y(k,125) + rxt(k,194)*y(k,134) + rxt(k,202)*y(k,226) &
                      + rxt(k,588)*y(k,225))
         mat(k,2036) = -rxt(k,190)*y(k,112)
         mat(k,1829) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112)
         mat(k,1534) = -rxt(k,194)*y(k,112)
         mat(k,1761) = -rxt(k,202)*y(k,112)
         mat(k,853) = -rxt(k,588)*y(k,112)
         mat(k,2306) = rxt(k,188)*y(k,217) + rxt(k,585)*y(k,220)
         mat(k,1534) = mat(k,1534) + rxt(k,586)*y(k,220)
         mat(k,864) = 1.100_r8*rxt(k,581)*y(k,218) + .200_r8*rxt(k,579)*y(k,219)
         mat(k,528) = rxt(k,188)*y(k,133)
         mat(k,679) = 1.100_r8*rxt(k,581)*y(k,201)
         mat(k,845) = .200_r8*rxt(k,579)*y(k,201)
         mat(k,504) = rxt(k,585)*y(k,133) + rxt(k,586)*y(k,134)
         mat(k,254) = -((rxt(k,206) + rxt(k,207)) * y(k,222))
         mat(k,1620) = -(rxt(k,206) + rxt(k,207)) * y(k,113)
         mat(k,963) = rxt(k,191)*y(k,125)
         mat(k,1814) = rxt(k,191)*y(k,112)
         mat(k,1815) = rxt(k,209)*y(k,126)
         mat(k,1860) = rxt(k,209)*y(k,125)
         mat(k,377) = -(rxt(k,428)*y(k,226))
         mat(k,1701) = -rxt(k,428)*y(k,115)
         mat(k,2213) = .200_r8*rxt(k,420)*y(k,208)
         mat(k,1242) = .200_r8*rxt(k,420)*y(k,199)
         mat(k,1115) = -(rxt(k,429)*y(k,226))
         mat(k,1771) = -rxt(k,429)*y(k,116)
         mat(k,2044) = rxt(k,422)*y(k,208)
         mat(k,1878) = rxt(k,423)*y(k,208)
         mat(k,1413) = rxt(k,419)*y(k,208)
         mat(k,2227) = .800_r8*rxt(k,420)*y(k,208)
         mat(k,1247) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) + rxt(k,419)*y(k,198) &
                      + .800_r8*rxt(k,420)*y(k,199)
         mat(k,96) = -(rxt(k,519)*y(k,226))
         mat(k,1662) = -rxt(k,519)*y(k,120)
         mat(k,2073) = -(rxt(k,190)*y(k,112) + rxt(k,199)*y(k,126) + rxt(k,203) &
                      *y(k,90) + rxt(k,204)*y(k,135) + rxt(k,205)*y(k,133) + rxt(k,226) &
                      *y(k,59) + rxt(k,258)*y(k,19) + rxt(k,301)*y(k,199) + rxt(k,309) &
                      *y(k,205) + rxt(k,322)*y(k,195) + rxt(k,333)*y(k,198) + rxt(k,337) &
                      *y(k,204) + rxt(k,350)*y(k,196) + rxt(k,359)*y(k,229) + rxt(k,363) &
                      *y(k,230) + (rxt(k,369) + rxt(k,370)) * y(k,202) + (rxt(k,376) &
                      + rxt(k,377)) * y(k,211) + rxt(k,385)*y(k,213) + rxt(k,388) &
                      *y(k,215) + (rxt(k,398) + rxt(k,399)) * y(k,190) + rxt(k,408) &
                      *y(k,206) + rxt(k,414)*y(k,207) + rxt(k,422)*y(k,208) + rxt(k,433) &
                      *y(k,235) + rxt(k,437)*y(k,189) + rxt(k,440)*y(k,192) + rxt(k,445) &
                      *y(k,194) + rxt(k,447)*y(k,197) + rxt(k,451)*y(k,200) + rxt(k,454) &
                      *y(k,212) + rxt(k,457)*y(k,214) + rxt(k,460)*y(k,228) + rxt(k,467) &
                      *y(k,233) + rxt(k,473)*y(k,236) + rxt(k,476)*y(k,238) + rxt(k,487) &
                      *y(k,221) + rxt(k,492)*y(k,231) + rxt(k,497)*y(k,232) + rxt(k,590) &
                      *y(k,225))
         mat(k,976) = -rxt(k,190)*y(k,124)
         mat(k,1909) = -rxt(k,199)*y(k,124)
         mat(k,2202) = -rxt(k,203)*y(k,124)
         mat(k,2412) = -rxt(k,204)*y(k,124)
         mat(k,2322) = -rxt(k,205)*y(k,124)
         mat(k,1597) = -rxt(k,226)*y(k,124)
         mat(k,1571) = -rxt(k,258)*y(k,124)
         mat(k,2254) = -rxt(k,301)*y(k,124)
         mat(k,442) = -rxt(k,309)*y(k,124)
         mat(k,904) = -rxt(k,322)*y(k,124)
         mat(k,1433) = -rxt(k,333)*y(k,124)
         mat(k,800) = -rxt(k,337)*y(k,124)
         mat(k,939) = -rxt(k,350)*y(k,124)
         mat(k,822) = -rxt(k,359)*y(k,124)
         mat(k,1220) = -rxt(k,363)*y(k,124)
         mat(k,577) = -(rxt(k,369) + rxt(k,370)) * y(k,124)
         mat(k,1362) = -(rxt(k,376) + rxt(k,377)) * y(k,124)
         mat(k,1402) = -rxt(k,385)*y(k,124)
         mat(k,701) = -rxt(k,388)*y(k,124)
         mat(k,1047) = -(rxt(k,398) + rxt(k,399)) * y(k,124)
         mat(k,1307) = -rxt(k,408)*y(k,124)
         mat(k,1340) = -rxt(k,414)*y(k,124)
         mat(k,1261) = -rxt(k,422)*y(k,124)
         mat(k,1238) = -rxt(k,433)*y(k,124)
         mat(k,525) = -rxt(k,437)*y(k,124)
         mat(k,491) = -rxt(k,440)*y(k,124)
         mat(k,437) = -rxt(k,445)*y(k,124)
         mat(k,636) = -rxt(k,447)*y(k,124)
         mat(k,781) = -rxt(k,451)*y(k,124)
         mat(k,742) = -rxt(k,454)*y(k,124)
         mat(k,919) = -rxt(k,457)*y(k,124)
         mat(k,456) = -rxt(k,460)*y(k,124)
         mat(k,757) = -rxt(k,467)*y(k,124)
         mat(k,774) = -rxt(k,473)*y(k,124)
         mat(k,513) = -rxt(k,476)*y(k,124)
         mat(k,1096) = -rxt(k,487)*y(k,124)
         mat(k,1183) = -rxt(k,492)*y(k,124)
         mat(k,1066) = -rxt(k,497)*y(k,124)
         mat(k,856) = -rxt(k,590)*y(k,124)
         mat(k,976) = mat(k,976) + 2.000_r8*rxt(k,192)*y(k,125) + rxt(k,194)*y(k,134) &
                      + rxt(k,202)*y(k,226)
         mat(k,257) = 2.000_r8*rxt(k,206)*y(k,222)
         mat(k,1849) = 2.000_r8*rxt(k,192)*y(k,112) + rxt(k,195)*y(k,133) + rxt(k,512) &
                      *y(k,151)
         mat(k,2322) = mat(k,2322) + rxt(k,195)*y(k,125)
         mat(k,1547) = rxt(k,194)*y(k,112) + rxt(k,189)*y(k,217)
         mat(k,1482) = rxt(k,512)*y(k,125)
         mat(k,531) = rxt(k,189)*y(k,134)
         mat(k,1639) = 2.000_r8*rxt(k,206)*y(k,113)
         mat(k,1804) = rxt(k,202)*y(k,112)
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
         mat(k,1845) = -((rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112) + (rxt(k,195) &
                      + rxt(k,197)) * y(k,133) + rxt(k,196)*y(k,135) + rxt(k,208) &
                      *y(k,90) + rxt(k,209)*y(k,126) + rxt(k,210)*y(k,226) + rxt(k,228) &
                      *y(k,59) + rxt(k,259)*y(k,19) + rxt(k,344)*y(k,198) + rxt(k,394) &
                      *y(k,213) + rxt(k,452)*y(k,200) + rxt(k,455)*y(k,212) + rxt(k,458) &
                      *y(k,214) + rxt(k,462)*y(k,142) + rxt(k,465)*y(k,189) + rxt(k,512) &
                      *y(k,151))
         mat(k,975) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,125)
         mat(k,2318) = -(rxt(k,195) + rxt(k,197)) * y(k,125)
         mat(k,2408) = -rxt(k,196)*y(k,125)
         mat(k,2198) = -rxt(k,208)*y(k,125)
         mat(k,1905) = -rxt(k,209)*y(k,125)
         mat(k,1800) = -rxt(k,210)*y(k,125)
         mat(k,1593) = -rxt(k,228)*y(k,125)
         mat(k,1567) = -rxt(k,259)*y(k,125)
         mat(k,1430) = -rxt(k,344)*y(k,125)
         mat(k,1399) = -rxt(k,394)*y(k,125)
         mat(k,780) = -rxt(k,452)*y(k,125)
         mat(k,741) = -rxt(k,455)*y(k,125)
         mat(k,918) = -rxt(k,458)*y(k,125)
         mat(k,468) = -rxt(k,462)*y(k,125)
         mat(k,524) = -rxt(k,465)*y(k,125)
         mat(k,1480) = -rxt(k,512)*y(k,125)
         mat(k,647) = rxt(k,396)*y(k,226)
         mat(k,358) = rxt(k,367)*y(k,126)
         mat(k,1567) = mat(k,1567) + rxt(k,258)*y(k,124)
         mat(k,1593) = mat(k,1593) + rxt(k,226)*y(k,124)
         mat(k,2198) = mat(k,2198) + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126)
         mat(k,496) = rxt(k,187)*y(k,226)
         mat(k,594) = .700_r8*rxt(k,416)*y(k,226)
         mat(k,2069) = rxt(k,258)*y(k,19) + rxt(k,226)*y(k,59) + rxt(k,203)*y(k,90) &
                      + 2.000_r8*rxt(k,199)*y(k,126) + rxt(k,205)*y(k,133) &
                      + rxt(k,204)*y(k,135) + rxt(k,437)*y(k,189) + rxt(k,398) &
                      *y(k,190) + rxt(k,440)*y(k,192) + rxt(k,445)*y(k,194) &
                      + rxt(k,322)*y(k,195) + rxt(k,350)*y(k,196) + rxt(k,447) &
                      *y(k,197) + rxt(k,333)*y(k,198) + rxt(k,301)*y(k,199) &
                      + rxt(k,451)*y(k,200) + rxt(k,369)*y(k,202) + rxt(k,337) &
                      *y(k,204) + rxt(k,309)*y(k,205) + .920_r8*rxt(k,408)*y(k,206) &
                      + .920_r8*rxt(k,414)*y(k,207) + rxt(k,422)*y(k,208) + rxt(k,376) &
                      *y(k,211) + rxt(k,454)*y(k,212) + rxt(k,385)*y(k,213) &
                      + rxt(k,457)*y(k,214) + rxt(k,388)*y(k,215) &
                      + 1.600_r8*rxt(k,487)*y(k,221) + rxt(k,460)*y(k,228) &
                      + rxt(k,359)*y(k,229) + rxt(k,363)*y(k,230) + .900_r8*rxt(k,492) &
                      *y(k,231) + .800_r8*rxt(k,497)*y(k,232) + rxt(k,467)*y(k,233) &
                      + rxt(k,433)*y(k,235) + rxt(k,473)*y(k,236) + rxt(k,476) &
                      *y(k,238)
         mat(k,1905) = mat(k,1905) + rxt(k,367)*y(k,16) + rxt(k,198)*y(k,90) &
                      + 2.000_r8*rxt(k,199)*y(k,124) + rxt(k,200)*y(k,133) &
                      + rxt(k,409)*y(k,206) + rxt(k,415)*y(k,207) + rxt(k,423) &
                      *y(k,208) + rxt(k,375)*y(k,211) + rxt(k,386)*y(k,213) &
                      + 2.000_r8*rxt(k,488)*y(k,221) + rxt(k,201)*y(k,226) &
                      + rxt(k,434)*y(k,235)
         mat(k,875) = rxt(k,357)*y(k,226)
         mat(k,2318) = mat(k,2318) + rxt(k,205)*y(k,124) + rxt(k,200)*y(k,126)
         mat(k,2408) = mat(k,2408) + rxt(k,204)*y(k,124)
         mat(k,628) = rxt(k,494)*y(k,226)
         mat(k,524) = mat(k,524) + rxt(k,437)*y(k,124)
         mat(k,1046) = rxt(k,398)*y(k,124)
         mat(k,490) = rxt(k,440)*y(k,124)
         mat(k,436) = rxt(k,445)*y(k,124)
         mat(k,903) = rxt(k,322)*y(k,124)
         mat(k,938) = rxt(k,350)*y(k,124)
         mat(k,635) = rxt(k,447)*y(k,124)
         mat(k,1430) = mat(k,1430) + rxt(k,333)*y(k,124)
         mat(k,2250) = rxt(k,301)*y(k,124) + .500_r8*rxt(k,485)*y(k,221)
         mat(k,780) = mat(k,780) + rxt(k,451)*y(k,124)
         mat(k,576) = rxt(k,369)*y(k,124)
         mat(k,799) = rxt(k,337)*y(k,124)
         mat(k,441) = rxt(k,309)*y(k,124)
         mat(k,1304) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126)
         mat(k,1337) = .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126)
         mat(k,1258) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126)
         mat(k,1359) = rxt(k,376)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,741) = mat(k,741) + rxt(k,454)*y(k,124)
         mat(k,1399) = mat(k,1399) + rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126)
         mat(k,918) = mat(k,918) + rxt(k,457)*y(k,124)
         mat(k,700) = rxt(k,388)*y(k,124)
         mat(k,1093) = 1.600_r8*rxt(k,487)*y(k,124) + 2.000_r8*rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1800) = mat(k,1800) + rxt(k,396)*y(k,1) + rxt(k,187)*y(k,91) &
                      + .700_r8*rxt(k,416)*y(k,100) + rxt(k,201)*y(k,126) + rxt(k,357) &
                      *y(k,127) + rxt(k,494)*y(k,176)
         mat(k,455) = rxt(k,460)*y(k,124)
         mat(k,821) = rxt(k,359)*y(k,124)
         mat(k,1219) = rxt(k,363)*y(k,124)
         mat(k,1181) = .900_r8*rxt(k,492)*y(k,124)
         mat(k,1064) = .800_r8*rxt(k,497)*y(k,124)
         mat(k,756) = rxt(k,467)*y(k,124)
         mat(k,1236) = rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126)
         mat(k,773) = rxt(k,473)*y(k,124)
         mat(k,512) = rxt(k,476)*y(k,124)
         mat(k,1906) = -(rxt(k,198)*y(k,90) + rxt(k,199)*y(k,124) + rxt(k,200) &
                      *y(k,133) + rxt(k,201)*y(k,226) + rxt(k,209)*y(k,125) + rxt(k,295) &
                      *y(k,42) + rxt(k,327)*y(k,45) + rxt(k,346)*y(k,29) + rxt(k,353) &
                      *y(k,49) + rxt(k,367)*y(k,16) + rxt(k,375)*y(k,211) + rxt(k,386) &
                      *y(k,213) + rxt(k,409)*y(k,206) + rxt(k,415)*y(k,207) + rxt(k,418) &
                      *y(k,99) + rxt(k,423)*y(k,208) + rxt(k,434)*y(k,235) + rxt(k,479) &
                      *y(k,6) + rxt(k,482)*y(k,110) + rxt(k,488)*y(k,221) + rxt(k,499) &
                      *y(k,178) + rxt(k,502)*y(k,67))
         mat(k,2199) = -rxt(k,198)*y(k,126)
         mat(k,2070) = -rxt(k,199)*y(k,126)
         mat(k,2319) = -rxt(k,200)*y(k,126)
         mat(k,1801) = -rxt(k,201)*y(k,126)
         mat(k,1846) = -rxt(k,209)*y(k,126)
         mat(k,2345) = -rxt(k,295)*y(k,126)
         mat(k,1153) = -rxt(k,327)*y(k,126)
         mat(k,1141) = -rxt(k,346)*y(k,126)
         mat(k,1283) = -rxt(k,353)*y(k,126)
         mat(k,359) = -rxt(k,367)*y(k,126)
         mat(k,1360) = -rxt(k,375)*y(k,126)
         mat(k,1400) = -rxt(k,386)*y(k,126)
         mat(k,1305) = -rxt(k,409)*y(k,126)
         mat(k,1338) = -rxt(k,415)*y(k,126)
         mat(k,891) = -rxt(k,418)*y(k,126)
         mat(k,1259) = -rxt(k,423)*y(k,126)
         mat(k,1237) = -rxt(k,434)*y(k,126)
         mat(k,1002) = -rxt(k,479)*y(k,126)
         mat(k,1030) = -rxt(k,482)*y(k,126)
         mat(k,1094) = -rxt(k,488)*y(k,126)
         mat(k,1076) = -rxt(k,499)*y(k,126)
         mat(k,311) = -rxt(k,502)*y(k,126)
         mat(k,569) = rxt(k,260)*y(k,133)
         mat(k,1952) = rxt(k,227)*y(k,60)
         mat(k,958) = rxt(k,227)*y(k,56) + rxt(k,229)*y(k,133) + rxt(k,230)*y(k,226)
         mat(k,926) = rxt(k,274)*y(k,89)
         mat(k,1975) = rxt(k,274)*y(k,73) + rxt(k,211)*y(k,226)
         mat(k,586) = .500_r8*rxt(k,391)*y(k,226)
         mat(k,1846) = mat(k,1846) + rxt(k,197)*y(k,133) + rxt(k,196)*y(k,135)
         mat(k,2319) = mat(k,2319) + rxt(k,260)*y(k,20) + rxt(k,229)*y(k,60) &
                      + rxt(k,197)*y(k,125)
         mat(k,2409) = rxt(k,196)*y(k,125)
         mat(k,537) = rxt(k,342)*y(k,226)
         mat(k,1801) = mat(k,1801) + rxt(k,230)*y(k,60) + rxt(k,211)*y(k,89) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,342)*y(k,140)
         mat(k,871) = -(rxt(k,357)*y(k,226))
         mat(k,1752) = -rxt(k,357)*y(k,127)
         mat(k,1129) = rxt(k,346)*y(k,126)
         mat(k,557) = .500_r8*rxt(k,417)*y(k,226)
         mat(k,421) = rxt(k,424)*y(k,226)
         mat(k,378) = rxt(k,428)*y(k,226)
         mat(k,1112) = rxt(k,429)*y(k,226)
         mat(k,1865) = rxt(k,346)*y(k,29)
         mat(k,1752) = mat(k,1752) + .500_r8*rxt(k,417)*y(k,101) + rxt(k,424)*y(k,102) &
                      + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116)
         mat(k,395) = -(rxt(k,489)*y(k,226))
         mat(k,1703) = -rxt(k,489)*y(k,128)
         mat(k,2124) = rxt(k,486)*y(k,221)
         mat(k,1084) = rxt(k,486)*y(k,90)
         mat(k,2327) = -(rxt(k,167)*y(k,135) + 4._r8*rxt(k,168)*y(k,133) + rxt(k,169) &
                      *y(k,134) + rxt(k,170)*y(k,77) + rxt(k,171)*y(k,79) + rxt(k,176) &
                      *y(k,90) + rxt(k,182)*y(k,226) + (rxt(k,195) + rxt(k,197) &
                      ) * y(k,125) + rxt(k,200)*y(k,126) + rxt(k,205)*y(k,124) &
                      + rxt(k,229)*y(k,60) + rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85) &
                      + rxt(k,237)*y(k,93) + rxt(k,260)*y(k,20) + rxt(k,261)*y(k,19) &
                      + rxt(k,263)*y(k,81) + rxt(k,265)*y(k,92) + rxt(k,296)*y(k,42) &
                      + rxt(k,504)*y(k,138) + (rxt(k,583) + rxt(k,584)) * y(k,218) &
                      + rxt(k,585)*y(k,220))
         mat(k,2417) = -rxt(k,167)*y(k,133)
         mat(k,1550) = -rxt(k,169)*y(k,133)
         mat(k,1469) = -rxt(k,170)*y(k,133)
         mat(k,616) = -rxt(k,171)*y(k,133)
         mat(k,2207) = -rxt(k,176)*y(k,133)
         mat(k,1809) = -rxt(k,182)*y(k,133)
         mat(k,1854) = -(rxt(k,195) + rxt(k,197)) * y(k,133)
         mat(k,1914) = -rxt(k,200)*y(k,133)
         mat(k,2078) = -rxt(k,205)*y(k,133)
         mat(k,961) = -rxt(k,229)*y(k,133)
         mat(k,1601) = -rxt(k,231)*y(k,133)
         mat(k,1504) = -rxt(k,234)*y(k,133)
         mat(k,840) = -rxt(k,237)*y(k,133)
         mat(k,571) = -rxt(k,260)*y(k,133)
         mat(k,1574) = -rxt(k,261)*y(k,133)
         mat(k,832) = -rxt(k,263)*y(k,133)
         mat(k,791) = -rxt(k,265)*y(k,133)
         mat(k,2353) = -rxt(k,296)*y(k,133)
         mat(k,368) = -rxt(k,504)*y(k,133)
         mat(k,683) = -(rxt(k,583) + rxt(k,584)) * y(k,133)
         mat(k,506) = -rxt(k,585)*y(k,133)
         mat(k,2099) = rxt(k,174)*y(k,90)
         mat(k,2207) = mat(k,2207) + rxt(k,174)*y(k,76)
         mat(k,978) = rxt(k,190)*y(k,124) + rxt(k,191)*y(k,125) + rxt(k,194)*y(k,134) &
                      + rxt(k,588)*y(k,225)
         mat(k,2078) = mat(k,2078) + rxt(k,190)*y(k,112)
         mat(k,1854) = mat(k,1854) + rxt(k,191)*y(k,112)
         mat(k,1550) = mat(k,1550) + rxt(k,194)*y(k,112) + rxt(k,506)*y(k,149) &
                      + rxt(k,513)*y(k,151) + rxt(k,587)*y(k,220) + (rxt(k,156) &
                       +rxt(k,157))*y(k,222) + rxt(k,593)*y(k,227)
         mat(k,713) = rxt(k,506)*y(k,134)
         mat(k,1485) = rxt(k,513)*y(k,134)
         mat(k,870) = rxt(k,579)*y(k,219) + 1.150_r8*rxt(k,580)*y(k,225)
         mat(k,849) = rxt(k,579)*y(k,201)
         mat(k,506) = mat(k,506) + rxt(k,587)*y(k,134)
         mat(k,1644) = (rxt(k,156)+rxt(k,157))*y(k,134)
         mat(k,857) = rxt(k,588)*y(k,112) + 1.150_r8*rxt(k,580)*y(k,201)
         mat(k,1809) = mat(k,1809) + 2.000_r8*rxt(k,184)*y(k,226)
         mat(k,810) = rxt(k,593)*y(k,134)
         mat(k,1540) = -(rxt(k,156)*y(k,222) + rxt(k,161)*y(k,223) + rxt(k,169) &
                      *y(k,133) + rxt(k,175)*y(k,76) + rxt(k,189)*y(k,217) + rxt(k,194) &
                      *y(k,112) + rxt(k,339)*y(k,203) + rxt(k,506)*y(k,149) + rxt(k,513) &
                      *y(k,151) + rxt(k,582)*y(k,218) + (rxt(k,586) + rxt(k,587) &
                      ) * y(k,220) + rxt(k,593)*y(k,227))
         mat(k,1630) = -rxt(k,156)*y(k,134)
         mat(k,176) = -rxt(k,161)*y(k,134)
         mat(k,2313) = -rxt(k,169)*y(k,134)
         mat(k,2085) = -rxt(k,175)*y(k,134)
         mat(k,529) = -rxt(k,189)*y(k,134)
         mat(k,972) = -rxt(k,194)*y(k,134)
         mat(k,463) = -rxt(k,339)*y(k,134)
         mat(k,710) = -rxt(k,506)*y(k,134)
         mat(k,1476) = -rxt(k,513)*y(k,134)
         mat(k,680) = -rxt(k,582)*y(k,134)
         mat(k,505) = -(rxt(k,586) + rxt(k,587)) * y(k,134)
         mat(k,809) = -rxt(k,593)*y(k,134)
         mat(k,1510) = rxt(k,251)*y(k,90) + rxt(k,252)*y(k,135)
         mat(k,1562) = 2.000_r8*rxt(k,253)*y(k,19) + (rxt(k,255)+rxt(k,256))*y(k,59) &
                      + rxt(k,257)*y(k,90) + rxt(k,261)*y(k,133)
         mat(k,1946) = rxt(k,218)*y(k,90) + rxt(k,220)*y(k,135)
         mat(k,1588) = (rxt(k,255)+rxt(k,256))*y(k,19) + (2.000_r8*rxt(k,222) &
                       +2.000_r8*rxt(k,223))*y(k,59) + rxt(k,225)*y(k,90) + rxt(k,231) &
                      *y(k,133) + rxt(k,233)*y(k,226)
         mat(k,2085) = mat(k,2085) + rxt(k,172)*y(k,90) + rxt(k,178)*y(k,135)
         mat(k,2193) = rxt(k,251)*y(k,17) + rxt(k,257)*y(k,19) + rxt(k,218)*y(k,56) &
                      + rxt(k,225)*y(k,59) + rxt(k,172)*y(k,76) + 2.000_r8*rxt(k,186) &
                      *y(k,90) + rxt(k,198)*y(k,126) + rxt(k,176)*y(k,133) &
                      + 2.000_r8*rxt(k,177)*y(k,135) + rxt(k,321)*y(k,195) &
                      + rxt(k,349)*y(k,196) + rxt(k,300)*y(k,199) + rxt(k,181) &
                      *y(k,226) + rxt(k,358)*y(k,229)
         mat(k,494) = rxt(k,187)*y(k,226)
         mat(k,972) = mat(k,972) + rxt(k,193)*y(k,125)
         mat(k,255) = rxt(k,207)*y(k,222)
         mat(k,2064) = rxt(k,204)*y(k,135) + rxt(k,590)*y(k,225)
         mat(k,1840) = rxt(k,193)*y(k,112) + rxt(k,195)*y(k,133) + rxt(k,196)*y(k,135)
         mat(k,1900) = rxt(k,198)*y(k,90) + rxt(k,200)*y(k,133)
         mat(k,2313) = mat(k,2313) + rxt(k,261)*y(k,19) + rxt(k,231)*y(k,59) &
                      + rxt(k,176)*y(k,90) + rxt(k,195)*y(k,125) + rxt(k,200)*y(k,126) &
                      + 2.000_r8*rxt(k,168)*y(k,133) + 2.000_r8*rxt(k,167)*y(k,135) &
                      + rxt(k,160)*y(k,223) + rxt(k,182)*y(k,226)
         mat(k,1540) = mat(k,1540) + 2.000_r8*rxt(k,161)*y(k,223)
         mat(k,2403) = rxt(k,252)*y(k,17) + rxt(k,220)*y(k,56) + rxt(k,178)*y(k,76) &
                      + 2.000_r8*rxt(k,177)*y(k,90) + rxt(k,204)*y(k,124) + rxt(k,196) &
                      *y(k,125) + 2.000_r8*rxt(k,167)*y(k,133) + rxt(k,508)*y(k,149) &
                      + rxt(k,514)*y(k,151) + 2.000_r8*rxt(k,158)*y(k,222) &
                      + rxt(k,183)*y(k,226)
         mat(k,710) = mat(k,710) + rxt(k,508)*y(k,135)
         mat(k,1476) = mat(k,1476) + rxt(k,514)*y(k,135)
         mat(k,901) = rxt(k,321)*y(k,90)
         mat(k,936) = rxt(k,349)*y(k,90)
         mat(k,2245) = rxt(k,300)*y(k,90)
         mat(k,1630) = mat(k,1630) + rxt(k,207)*y(k,113) + 2.000_r8*rxt(k,158) &
                      *y(k,135)
         mat(k,176) = mat(k,176) + rxt(k,160)*y(k,133) + 2.000_r8*rxt(k,161)*y(k,134)
         mat(k,854) = rxt(k,590)*y(k,124)
         mat(k,1795) = rxt(k,233)*y(k,59) + rxt(k,181)*y(k,90) + rxt(k,187)*y(k,91) &
                      + rxt(k,182)*y(k,133) + rxt(k,183)*y(k,135)
         mat(k,819) = rxt(k,358)*y(k,90)
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
         mat(k,2419) = -(rxt(k,158)*y(k,222) + rxt(k,167)*y(k,133) + rxt(k,177) &
                      *y(k,90) + rxt(k,178)*y(k,76) + rxt(k,183)*y(k,226) + rxt(k,196) &
                      *y(k,125) + rxt(k,204)*y(k,124) + rxt(k,220)*y(k,56) + rxt(k,252) &
                      *y(k,17) + rxt(k,318)*y(k,25) + rxt(k,347)*y(k,29) + rxt(k,378) &
                      *y(k,105) + rxt(k,392)*y(k,111) + rxt(k,425)*y(k,99) + rxt(k,463) &
                      *y(k,142) + rxt(k,480)*y(k,6) + rxt(k,483)*y(k,110) + rxt(k,508) &
                      *y(k,149) + rxt(k,514)*y(k,151))
         mat(k,1646) = -rxt(k,158)*y(k,135)
         mat(k,2329) = -rxt(k,167)*y(k,135)
         mat(k,2209) = -rxt(k,177)*y(k,135)
         mat(k,2101) = -rxt(k,178)*y(k,135)
         mat(k,1811) = -rxt(k,183)*y(k,135)
         mat(k,1856) = -rxt(k,196)*y(k,135)
         mat(k,2080) = -rxt(k,204)*y(k,135)
         mat(k,1962) = -rxt(k,220)*y(k,135)
         mat(k,1520) = -rxt(k,252)*y(k,135)
         mat(k,555) = -rxt(k,318)*y(k,135)
         mat(k,1147) = -rxt(k,347)*y(k,135)
         mat(k,1276) = -rxt(k,378)*y(k,135)
         mat(k,1388) = -rxt(k,392)*y(k,135)
         mat(k,895) = -rxt(k,425)*y(k,135)
         mat(k,469) = -rxt(k,463)*y(k,135)
         mat(k,1005) = -rxt(k,480)*y(k,135)
         mat(k,1033) = -rxt(k,483)*y(k,135)
         mat(k,714) = -rxt(k,508)*y(k,135)
         mat(k,1486) = -rxt(k,514)*y(k,135)
         mat(k,2209) = mat(k,2209) + .150_r8*rxt(k,332)*y(k,198) + .150_r8*rxt(k,383) &
                      *y(k,213)
         mat(k,2329) = mat(k,2329) + rxt(k,169)*y(k,134)
         mat(k,1552) = rxt(k,169)*y(k,133)
         mat(k,1438) = .150_r8*rxt(k,332)*y(k,90)
         mat(k,1406) = .150_r8*rxt(k,383)*y(k,90)
         mat(k,332) = -(rxt(k,515)*y(k,151))
         mat(k,1471) = -rxt(k,515)*y(k,137)
         mat(k,1555) = rxt(k,254)*y(k,59)
         mat(k,1581) = rxt(k,254)*y(k,19) + 2.000_r8*rxt(k,224)*y(k,59)
         mat(k,361) = -(rxt(k,504)*y(k,133) + rxt(k,505)*y(k,226))
         mat(k,2290) = -rxt(k,504)*y(k,138)
         mat(k,1699) = -rxt(k,505)*y(k,138)
         mat(k,1188) = rxt(k,371)*y(k,226)
         mat(k,1999) = .100_r8*rxt(k,492)*y(k,231)
         mat(k,1680) = rxt(k,371)*y(k,94)
         mat(k,1169) = .100_r8*rxt(k,492)*y(k,124)
         mat(k,532) = -(rxt(k,342)*y(k,226))
         mat(k,1721) = -rxt(k,342)*y(k,140)
         mat(k,1819) = rxt(k,344)*y(k,198)
         mat(k,1409) = rxt(k,344)*y(k,125)
         mat(k,1813) = rxt(k,465)*y(k,189)
         mat(k,520) = rxt(k,465)*y(k,125)
         mat(k,466) = -(rxt(k,462)*y(k,125) + rxt(k,463)*y(k,135))
         mat(k,1816) = -rxt(k,462)*y(k,142)
         mat(k,2366) = -rxt(k,463)*y(k,142)
         mat(k,194) = .070_r8*rxt(k,449)*y(k,226)
         mat(k,2009) = rxt(k,447)*y(k,197)
         mat(k,170) = .060_r8*rxt(k,461)*y(k,226)
         mat(k,215) = .070_r8*rxt(k,477)*y(k,226)
         mat(k,632) = rxt(k,447)*y(k,124)
         mat(k,1713) = .070_r8*rxt(k,449)*y(k,66) + .060_r8*rxt(k,461)*y(k,143) &
                      + .070_r8*rxt(k,477)*y(k,185)
         mat(k,168) = -(rxt(k,461)*y(k,226))
         mat(k,1668) = -rxt(k,461)*y(k,143)
         mat(k,160) = .530_r8*rxt(k,438)*y(k,226)
         mat(k,1668) = mat(k,1668) + .530_r8*rxt(k,438)*y(k,7)
         mat(k,337) = -(rxt(k,464)*y(k,226))
         mat(k,1695) = -rxt(k,464)*y(k,144)
         mat(k,2120) = rxt(k,459)*y(k,228)
         mat(k,451) = rxt(k,459)*y(k,90)
         mat(k,540) = -(rxt(k,360)*y(k,226))
         mat(k,1722) = -rxt(k,360)*y(k,147)
         mat(k,2141) = rxt(k,358)*y(k,229)
         mat(k,815) = rxt(k,358)*y(k,90)
         mat(k,407) = -(rxt(k,364)*y(k,226))
         mat(k,1705) = -rxt(k,364)*y(k,148)
         mat(k,2126) = .850_r8*rxt(k,362)*y(k,230)
         mat(k,1212) = .850_r8*rxt(k,362)*y(k,90)
         mat(k,708) = -(rxt(k,506)*y(k,134) + rxt(k,508)*y(k,135) + rxt(k,511) &
                      *y(k,226))
         mat(k,1528) = -rxt(k,506)*y(k,149)
         mat(k,2370) = -rxt(k,508)*y(k,149)
         mat(k,1740) = -rxt(k,511)*y(k,149)
         mat(k,1474) = -(rxt(k,509)*y(k,19) + rxt(k,510)*y(k,59) + rxt(k,512)*y(k,125) &
                      + rxt(k,513)*y(k,134) + rxt(k,514)*y(k,135) + rxt(k,515) &
                      *y(k,137) + rxt(k,516)*y(k,226))
         mat(k,1559) = -rxt(k,509)*y(k,151)
         mat(k,1585) = -rxt(k,510)*y(k,151)
         mat(k,1837) = -rxt(k,512)*y(k,151)
         mat(k,1538) = -rxt(k,513)*y(k,151)
         mat(k,2401) = -rxt(k,514)*y(k,151)
         mat(k,334) = -rxt(k,515)*y(k,151)
         mat(k,1792) = -rxt(k,516)*y(k,151)
         mat(k,2310) = rxt(k,504)*y(k,138)
         mat(k,1538) = mat(k,1538) + rxt(k,506)*y(k,149)
         mat(k,2401) = mat(k,2401) + rxt(k,508)*y(k,149)
         mat(k,365) = rxt(k,504)*y(k,133)
         mat(k,709) = rxt(k,506)*y(k,134) + rxt(k,508)*y(k,135) + rxt(k,511)*y(k,226)
         mat(k,1792) = mat(k,1792) + rxt(k,511)*y(k,149)
         mat(k,945) = -(rxt(k,507)*y(k,226))
         mat(k,1759) = -rxt(k,507)*y(k,152)
         mat(k,1558) = rxt(k,509)*y(k,151)
         mat(k,1583) = rxt(k,510)*y(k,151)
         mat(k,309) = rxt(k,502)*y(k,126) + (rxt(k,503)+.500_r8*rxt(k,517))*y(k,226)
         mat(k,1827) = rxt(k,512)*y(k,151)
         mat(k,1868) = rxt(k,502)*y(k,67)
         mat(k,1533) = rxt(k,513)*y(k,151)
         mat(k,2374) = rxt(k,514)*y(k,151)
         mat(k,333) = rxt(k,515)*y(k,151)
         mat(k,363) = rxt(k,505)*y(k,226)
         mat(k,1473) = rxt(k,509)*y(k,19) + rxt(k,510)*y(k,59) + rxt(k,512)*y(k,125) &
                      + rxt(k,513)*y(k,134) + rxt(k,514)*y(k,135) + rxt(k,515) &
                      *y(k,137) + rxt(k,516)*y(k,226)
         mat(k,1759) = mat(k,1759) + (rxt(k,503)+.500_r8*rxt(k,517))*y(k,67) &
                      + rxt(k,505)*y(k,138) + rxt(k,516)*y(k,151)
         mat(k,259) = -(rxt(k,518)*y(k,239))
         mat(k,2422) = -rxt(k,518)*y(k,153)
         mat(k,944) = rxt(k,507)*y(k,226)
         mat(k,1684) = rxt(k,507)*y(k,152)
         mat(k,979) = .2202005_r8*rxt(k,537)*y(k,135)
         mat(k,2103) = .2202005_r8*rxt(k,535)*y(k,191) + .0023005_r8*rxt(k,540) &
                      *y(k,193) + .0031005_r8*rxt(k,543)*y(k,209) &
                      + .2381005_r8*rxt(k,547)*y(k,210) + .0508005_r8*rxt(k,551) &
                      *y(k,216) + .1364005_r8*rxt(k,557)*y(k,234) &
                      + .1677005_r8*rxt(k,560)*y(k,237)
         mat(k,1007) = .0508005_r8*rxt(k,553)*y(k,135)
         mat(k,1987) = .1279005_r8*rxt(k,536)*y(k,191) + .0097005_r8*rxt(k,541) &
                      *y(k,193) + .0003005_r8*rxt(k,544)*y(k,209) &
                      + .1056005_r8*rxt(k,548)*y(k,210) + .0245005_r8*rxt(k,552) &
                      *y(k,216) + .0154005_r8*rxt(k,558)*y(k,234) &
                      + .0063005_r8*rxt(k,561)*y(k,237)
         mat(k,2357) = .2202005_r8*rxt(k,537)*y(k,6) + .0508005_r8*rxt(k,553)*y(k,110)
         mat(k,43) = .5931005_r8*rxt(k,555)*y(k,226)
         mat(k,49) = .2202005_r8*rxt(k,535)*y(k,90) + .1279005_r8*rxt(k,536)*y(k,124)
         mat(k,55) = .0023005_r8*rxt(k,540)*y(k,90) + .0097005_r8*rxt(k,541)*y(k,124)
         mat(k,61) = .0031005_r8*rxt(k,543)*y(k,90) + .0003005_r8*rxt(k,544)*y(k,124)
         mat(k,67) = .2381005_r8*rxt(k,547)*y(k,90) + .1056005_r8*rxt(k,548)*y(k,124)
         mat(k,75) = .0508005_r8*rxt(k,551)*y(k,90) + .0245005_r8*rxt(k,552)*y(k,124)
         mat(k,1648) = .5931005_r8*rxt(k,555)*y(k,173)
         mat(k,81) = .1364005_r8*rxt(k,557)*y(k,90) + .0154005_r8*rxt(k,558)*y(k,124)
         mat(k,87) = .1677005_r8*rxt(k,560)*y(k,90) + .0063005_r8*rxt(k,561)*y(k,124)
         mat(k,980) = .2067005_r8*rxt(k,537)*y(k,135)
         mat(k,2104) = .2067005_r8*rxt(k,535)*y(k,191) + .0008005_r8*rxt(k,540) &
                      *y(k,193) + .0035005_r8*rxt(k,543)*y(k,209) &
                      + .1308005_r8*rxt(k,547)*y(k,210) + .1149005_r8*rxt(k,551) &
                      *y(k,216) + .0101005_r8*rxt(k,557)*y(k,234) &
                      + .0174005_r8*rxt(k,560)*y(k,237)
         mat(k,1008) = .1149005_r8*rxt(k,553)*y(k,135)
         mat(k,1988) = .1792005_r8*rxt(k,536)*y(k,191) + .0034005_r8*rxt(k,541) &
                      *y(k,193) + .0003005_r8*rxt(k,544)*y(k,209) &
                      + .1026005_r8*rxt(k,548)*y(k,210) + .0082005_r8*rxt(k,552) &
                      *y(k,216) + .0452005_r8*rxt(k,558)*y(k,234) &
                      + .0237005_r8*rxt(k,561)*y(k,237)
         mat(k,2358) = .2067005_r8*rxt(k,537)*y(k,6) + .1149005_r8*rxt(k,553)*y(k,110)
         mat(k,44) = .1534005_r8*rxt(k,555)*y(k,226)
         mat(k,50) = .2067005_r8*rxt(k,535)*y(k,90) + .1792005_r8*rxt(k,536)*y(k,124)
         mat(k,56) = .0008005_r8*rxt(k,540)*y(k,90) + .0034005_r8*rxt(k,541)*y(k,124)
         mat(k,62) = .0035005_r8*rxt(k,543)*y(k,90) + .0003005_r8*rxt(k,544)*y(k,124)
         mat(k,68) = .1308005_r8*rxt(k,547)*y(k,90) + .1026005_r8*rxt(k,548)*y(k,124)
         mat(k,76) = .1149005_r8*rxt(k,551)*y(k,90) + .0082005_r8*rxt(k,552)*y(k,124)
         mat(k,1649) = .1534005_r8*rxt(k,555)*y(k,173)
         mat(k,82) = .0101005_r8*rxt(k,557)*y(k,90) + .0452005_r8*rxt(k,558)*y(k,124)
         mat(k,88) = .0174005_r8*rxt(k,560)*y(k,90) + .0237005_r8*rxt(k,561)*y(k,124)
         mat(k,981) = .0653005_r8*rxt(k,537)*y(k,135)
         mat(k,2105) = .0653005_r8*rxt(k,535)*y(k,191) + .0843005_r8*rxt(k,540) &
                      *y(k,193) + .0003005_r8*rxt(k,543)*y(k,209) &
                      + .0348005_r8*rxt(k,547)*y(k,210) + .0348005_r8*rxt(k,551) &
                      *y(k,216) + .0763005_r8*rxt(k,557)*y(k,234) + .086_r8*rxt(k,560) &
                      *y(k,237)
         mat(k,1009) = .0348005_r8*rxt(k,553)*y(k,135)
         mat(k,1989) = .0676005_r8*rxt(k,536)*y(k,191) + .1579005_r8*rxt(k,541) &
                      *y(k,193) + .0073005_r8*rxt(k,544)*y(k,209) &
                      + .0521005_r8*rxt(k,548)*y(k,210) + .0772005_r8*rxt(k,552) &
                      *y(k,216) + .0966005_r8*rxt(k,558)*y(k,234) &
                      + .0025005_r8*rxt(k,561)*y(k,237)
         mat(k,2359) = .0653005_r8*rxt(k,537)*y(k,6) + .0348005_r8*rxt(k,553)*y(k,110)
         mat(k,45) = .0459005_r8*rxt(k,555)*y(k,226)
         mat(k,51) = .0653005_r8*rxt(k,535)*y(k,90) + .0676005_r8*rxt(k,536)*y(k,124)
         mat(k,57) = .0843005_r8*rxt(k,540)*y(k,90) + .1579005_r8*rxt(k,541)*y(k,124)
         mat(k,63) = .0003005_r8*rxt(k,543)*y(k,90) + .0073005_r8*rxt(k,544)*y(k,124)
         mat(k,69) = .0348005_r8*rxt(k,547)*y(k,90) + .0521005_r8*rxt(k,548)*y(k,124)
         mat(k,77) = .0348005_r8*rxt(k,551)*y(k,90) + .0772005_r8*rxt(k,552)*y(k,124)
         mat(k,1650) = .0459005_r8*rxt(k,555)*y(k,173)
         mat(k,83) = .0763005_r8*rxt(k,557)*y(k,90) + .0966005_r8*rxt(k,558)*y(k,124)
         mat(k,89) = .086_r8*rxt(k,560)*y(k,90) + .0025005_r8*rxt(k,561)*y(k,124)
         mat(k,982) = .1749305_r8*rxt(k,534)*y(k,126) + .1284005_r8*rxt(k,537) &
                      *y(k,135)
         mat(k,2106) = .1284005_r8*rxt(k,535)*y(k,191) + .0443005_r8*rxt(k,540) &
                      *y(k,193) + .0271005_r8*rxt(k,543)*y(k,209) &
                      + .0076005_r8*rxt(k,547)*y(k,210) + .0554005_r8*rxt(k,551) &
                      *y(k,216) + .2157005_r8*rxt(k,557)*y(k,234) &
                      + .0512005_r8*rxt(k,560)*y(k,237)
         mat(k,877) = .0590245_r8*rxt(k,542)*y(k,126) + .0033005_r8*rxt(k,545) &
                      *y(k,135)
         mat(k,1010) = .1749305_r8*rxt(k,550)*y(k,126) + .0554005_r8*rxt(k,553) &
                      *y(k,135)
         mat(k,1990) = .079_r8*rxt(k,536)*y(k,191) + .0059005_r8*rxt(k,541)*y(k,193) &
                      + .0057005_r8*rxt(k,544)*y(k,209) + .0143005_r8*rxt(k,548) &
                      *y(k,210) + .0332005_r8*rxt(k,552)*y(k,216) &
                      + .0073005_r8*rxt(k,558)*y(k,234) + .011_r8*rxt(k,561)*y(k,237)
         mat(k,1858) = .1749305_r8*rxt(k,534)*y(k,6) + .0590245_r8*rxt(k,542)*y(k,99) &
                      + .1749305_r8*rxt(k,550)*y(k,110)
         mat(k,2360) = .1284005_r8*rxt(k,537)*y(k,6) + .0033005_r8*rxt(k,545)*y(k,99) &
                      + .0554005_r8*rxt(k,553)*y(k,110)
         mat(k,46) = .0085005_r8*rxt(k,555)*y(k,226)
         mat(k,52) = .1284005_r8*rxt(k,535)*y(k,90) + .079_r8*rxt(k,536)*y(k,124)
         mat(k,58) = .0443005_r8*rxt(k,540)*y(k,90) + .0059005_r8*rxt(k,541)*y(k,124)
         mat(k,64) = .0271005_r8*rxt(k,543)*y(k,90) + .0057005_r8*rxt(k,544)*y(k,124)
         mat(k,70) = .0076005_r8*rxt(k,547)*y(k,90) + .0143005_r8*rxt(k,548)*y(k,124)
         mat(k,78) = .0554005_r8*rxt(k,551)*y(k,90) + .0332005_r8*rxt(k,552)*y(k,124)
         mat(k,1651) = .0085005_r8*rxt(k,555)*y(k,173)
         mat(k,84) = .2157005_r8*rxt(k,557)*y(k,90) + .0073005_r8*rxt(k,558)*y(k,124)
         mat(k,90) = .0512005_r8*rxt(k,560)*y(k,90) + .011_r8*rxt(k,561)*y(k,124)
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
         mat(k,983) = .5901905_r8*rxt(k,534)*y(k,126) + .114_r8*rxt(k,537)*y(k,135)
         mat(k,2107) = .114_r8*rxt(k,535)*y(k,191) + .1621005_r8*rxt(k,540)*y(k,193) &
                      + .0474005_r8*rxt(k,543)*y(k,209) + .0113005_r8*rxt(k,547) &
                      *y(k,210) + .1278005_r8*rxt(k,551)*y(k,216) &
                      + .0738005_r8*rxt(k,557)*y(k,234) + .1598005_r8*rxt(k,560) &
                      *y(k,237)
         mat(k,878) = .0250245_r8*rxt(k,542)*y(k,126)
         mat(k,1011) = .5901905_r8*rxt(k,550)*y(k,126) + .1278005_r8*rxt(k,553) &
                      *y(k,135)
         mat(k,1991) = .1254005_r8*rxt(k,536)*y(k,191) + .0536005_r8*rxt(k,541) &
                      *y(k,193) + .0623005_r8*rxt(k,544)*y(k,209) &
                      + .0166005_r8*rxt(k,548)*y(k,210) + .130_r8*rxt(k,552)*y(k,216) &
                      + .238_r8*rxt(k,558)*y(k,234) + .1185005_r8*rxt(k,561)*y(k,237)
         mat(k,1859) = .5901905_r8*rxt(k,534)*y(k,6) + .0250245_r8*rxt(k,542)*y(k,99) &
                      + .5901905_r8*rxt(k,550)*y(k,110)
         mat(k,2361) = .114_r8*rxt(k,537)*y(k,6) + .1278005_r8*rxt(k,553)*y(k,110)
         mat(k,47) = .0128005_r8*rxt(k,555)*y(k,226)
         mat(k,53) = .114_r8*rxt(k,535)*y(k,90) + .1254005_r8*rxt(k,536)*y(k,124)
         mat(k,59) = .1621005_r8*rxt(k,540)*y(k,90) + .0536005_r8*rxt(k,541)*y(k,124)
         mat(k,65) = .0474005_r8*rxt(k,543)*y(k,90) + .0623005_r8*rxt(k,544)*y(k,124)
         mat(k,71) = .0113005_r8*rxt(k,547)*y(k,90) + .0166005_r8*rxt(k,548)*y(k,124)
         mat(k,79) = .1278005_r8*rxt(k,551)*y(k,90) + .130_r8*rxt(k,552)*y(k,124)
         mat(k,1652) = .0128005_r8*rxt(k,555)*y(k,173)
         mat(k,85) = .0738005_r8*rxt(k,557)*y(k,90) + .238_r8*rxt(k,558)*y(k,124)
         mat(k,91) = .1598005_r8*rxt(k,560)*y(k,90) + .1185005_r8*rxt(k,561)*y(k,124)
         mat(k,48) = -(rxt(k,555)*y(k,226))
         mat(k,1653) = -rxt(k,555)*y(k,173)
         mat(k,187) = .100_r8*rxt(k,469)*y(k,226)
         mat(k,205) = .230_r8*rxt(k,471)*y(k,226)
         mat(k,1672) = .100_r8*rxt(k,469)*y(k,181) + .230_r8*rxt(k,471)*y(k,183)
         mat(k,650) = -(rxt(k,493)*y(k,226))
         mat(k,1735) = -rxt(k,493)*y(k,175)
         mat(k,2145) = rxt(k,491)*y(k,231)
         mat(k,1170) = rxt(k,491)*y(k,90)
         mat(k,625) = -(rxt(k,494)*y(k,226))
         mat(k,1732) = -rxt(k,494)*y(k,176)
         mat(k,2018) = .200_r8*rxt(k,487)*y(k,221) + .200_r8*rxt(k,497)*y(k,232)
         mat(k,2215) = .500_r8*rxt(k,485)*y(k,221)
         mat(k,1085) = .200_r8*rxt(k,487)*y(k,124) + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1058) = .200_r8*rxt(k,497)*y(k,124)
         mat(k,477) = -(rxt(k,498)*y(k,226))
         mat(k,1715) = -rxt(k,498)*y(k,177)
         mat(k,2136) = rxt(k,496)*y(k,232)
         mat(k,1057) = rxt(k,496)*y(k,90)
         mat(k,1070) = -(rxt(k,499)*y(k,126) + rxt(k,500)*y(k,226))
         mat(k,1874) = -rxt(k,499)*y(k,178)
         mat(k,1767) = -rxt(k,500)*y(k,178)
         mat(k,992) = .330_r8*rxt(k,480)*y(k,135)
         mat(k,1020) = .330_r8*rxt(k,483)*y(k,135)
         mat(k,2040) = .800_r8*rxt(k,487)*y(k,221) + .800_r8*rxt(k,497)*y(k,232)
         mat(k,1874) = mat(k,1874) + rxt(k,488)*y(k,221)
         mat(k,2380) = .330_r8*rxt(k,480)*y(k,6) + .330_r8*rxt(k,483)*y(k,110)
         mat(k,626) = rxt(k,494)*y(k,226)
         mat(k,2223) = .500_r8*rxt(k,485)*y(k,221) + rxt(k,495)*y(k,232)
         mat(k,1087) = .800_r8*rxt(k,487)*y(k,124) + rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1767) = mat(k,1767) + rxt(k,494)*y(k,176)
         mat(k,1061) = .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,199)
         mat(k,1101) = -(rxt(k,501)*y(k,226))
         mat(k,1770) = -rxt(k,501)*y(k,179)
         mat(k,995) = .300_r8*rxt(k,480)*y(k,135)
         mat(k,1023) = .300_r8*rxt(k,483)*y(k,135)
         mat(k,2043) = .900_r8*rxt(k,492)*y(k,231)
         mat(k,2383) = .300_r8*rxt(k,480)*y(k,6) + .300_r8*rxt(k,483)*y(k,110)
         mat(k,2226) = rxt(k,490)*y(k,231)
         mat(k,1173) = .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,199)
         mat(k,663) = -(rxt(k,468)*y(k,226))
         mat(k,1736) = -rxt(k,468)*y(k,180)
         mat(k,2146) = rxt(k,466)*y(k,233)
         mat(k,747) = rxt(k,466)*y(k,90)
         mat(k,185) = -(rxt(k,469)*y(k,226))
         mat(k,1670) = -rxt(k,469)*y(k,181)
         mat(k,201) = -(rxt(k,435)*y(k,226))
         mat(k,1673) = -rxt(k,435)*y(k,182)
         mat(k,2116) = rxt(k,432)*y(k,235)
         mat(k,1225) = rxt(k,432)*y(k,90)
         mat(k,206) = -(rxt(k,471)*y(k,226))
         mat(k,1674) = -rxt(k,471)*y(k,183)
         mat(k,719) = -(rxt(k,474)*y(k,226))
         mat(k,1741) = -rxt(k,474)*y(k,184)
         mat(k,2150) = rxt(k,472)*y(k,236)
         mat(k,763) = rxt(k,472)*y(k,90)
         mat(k,214) = -(rxt(k,477)*y(k,226))
         mat(k,1675) = -rxt(k,477)*y(k,185)
         mat(k,207) = .150_r8*rxt(k,471)*y(k,226)
         mat(k,1675) = mat(k,1675) + .150_r8*rxt(k,471)*y(k,183)
         mat(k,425) = -(rxt(k,478)*y(k,226))
         mat(k,1708) = -rxt(k,478)*y(k,186)
         mat(k,2129) = rxt(k,475)*y(k,238)
         mat(k,507) = rxt(k,475)*y(k,90)
         mat(k,521) = -(rxt(k,436)*y(k,90) + rxt(k,437)*y(k,124) + rxt(k,465)*y(k,125))
         mat(k,2140) = -rxt(k,436)*y(k,189)
         mat(k,2013) = -rxt(k,437)*y(k,189)
         mat(k,1818) = -rxt(k,465)*y(k,189)
         mat(k,246) = rxt(k,442)*y(k,226)
         mat(k,1720) = rxt(k,442)*y(k,22)
         mat(k,1040) = -(rxt(k,397)*y(k,90) + (rxt(k,398) + rxt(k,399)) * y(k,124))
         mat(k,2166) = -rxt(k,397)*y(k,190)
         mat(k,2037) = -(rxt(k,398) + rxt(k,399)) * y(k,190)
         mat(k,688) = rxt(k,400)*y(k,226)
         mat(k,234) = rxt(k,401)*y(k,226)
         mat(k,1764) = rxt(k,400)*y(k,2) + rxt(k,401)*y(k,15)
         mat(k,54) = -(rxt(k,535)*y(k,90) + rxt(k,536)*y(k,124))
         mat(k,2108) = -rxt(k,535)*y(k,191)
         mat(k,1992) = -rxt(k,536)*y(k,191)
         mat(k,984) = rxt(k,538)*y(k,226)
         mat(k,1654) = rxt(k,538)*y(k,6)
         mat(k,486) = -(rxt(k,439)*y(k,90) + rxt(k,440)*y(k,124))
         mat(k,2137) = -rxt(k,439)*y(k,192)
         mat(k,2010) = -rxt(k,440)*y(k,192)
         mat(k,161) = .350_r8*rxt(k,438)*y(k,226)
         mat(k,403) = rxt(k,441)*y(k,226)
         mat(k,1716) = .350_r8*rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8)
         mat(k,60) = -(rxt(k,540)*y(k,90) + rxt(k,541)*y(k,124))
         mat(k,2109) = -rxt(k,540)*y(k,193)
         mat(k,1993) = -rxt(k,541)*y(k,193)
         mat(k,157) = rxt(k,539)*y(k,226)
         mat(k,1655) = rxt(k,539)*y(k,7)
         mat(k,433) = -(rxt(k,443)*y(k,90) + rxt(k,445)*y(k,124))
         mat(k,2130) = -rxt(k,443)*y(k,194)
         mat(k,2004) = -rxt(k,445)*y(k,194)
         mat(k,344) = rxt(k,444)*y(k,226)
         mat(k,188) = .070_r8*rxt(k,469)*y(k,226)
         mat(k,208) = .060_r8*rxt(k,471)*y(k,226)
         mat(k,1709) = rxt(k,444)*y(k,23) + .070_r8*rxt(k,469)*y(k,181) &
                      + .060_r8*rxt(k,471)*y(k,183)
         mat(k,899) = -(4._r8*rxt(k,319)*y(k,195) + rxt(k,320)*y(k,199) + rxt(k,321) &
                      *y(k,90) + rxt(k,322)*y(k,124))
         mat(k,2219) = -rxt(k,320)*y(k,195)
         mat(k,2162) = -rxt(k,321)*y(k,195)
         mat(k,2032) = -rxt(k,322)*y(k,195)
         mat(k,349) = .500_r8*rxt(k,324)*y(k,226)
         mat(k,297) = rxt(k,325)*y(k,56) + rxt(k,326)*y(k,226)
         mat(k,1930) = rxt(k,325)*y(k,28)
         mat(k,1754) = .500_r8*rxt(k,324)*y(k,27) + rxt(k,326)*y(k,28)
         mat(k,933) = -(rxt(k,348)*y(k,199) + rxt(k,349)*y(k,90) + rxt(k,350)*y(k,124))
         mat(k,2220) = -rxt(k,348)*y(k,196)
         mat(k,2165) = -rxt(k,349)*y(k,196)
         mat(k,2035) = -rxt(k,350)*y(k,196)
         mat(k,414) = rxt(k,351)*y(k,226)
         mat(k,303) = rxt(k,355)*y(k,56) + rxt(k,352)*y(k,226)
         mat(k,1932) = rxt(k,355)*y(k,31)
         mat(k,1758) = rxt(k,351)*y(k,30) + rxt(k,352)*y(k,31)
         mat(k,633) = -(rxt(k,446)*y(k,90) + rxt(k,447)*y(k,124))
         mat(k,2144) = -rxt(k,446)*y(k,197)
         mat(k,2019) = -rxt(k,447)*y(k,197)
         mat(k,269) = rxt(k,448)*y(k,226)
         mat(k,2144) = mat(k,2144) + .400_r8*rxt(k,436)*y(k,189)
         mat(k,2019) = mat(k,2019) + rxt(k,437)*y(k,189)
         mat(k,2368) = rxt(k,463)*y(k,142)
         mat(k,467) = rxt(k,463)*y(k,135)
         mat(k,522) = .400_r8*rxt(k,436)*y(k,90) + rxt(k,437)*y(k,124)
         mat(k,1733) = rxt(k,448)*y(k,32)
         mat(k,1426) = -(4._r8*rxt(k,330)*y(k,198) + rxt(k,331)*y(k,199) + rxt(k,332) &
                      *y(k,90) + rxt(k,333)*y(k,124) + rxt(k,344)*y(k,125) + rxt(k,372) &
                      *y(k,211) + rxt(k,405)*y(k,206) + rxt(k,410)*y(k,207) + rxt(k,419) &
                      *y(k,208) + rxt(k,430)*y(k,235))
         mat(k,2243) = -rxt(k,331)*y(k,198)
         mat(k,2188) = -rxt(k,332)*y(k,198)
         mat(k,2061) = -rxt(k,333)*y(k,198)
         mat(k,1835) = -rxt(k,344)*y(k,198)
         mat(k,1356) = -rxt(k,372)*y(k,198)
         mat(k,1301) = -rxt(k,405)*y(k,198)
         mat(k,1334) = -rxt(k,410)*y(k,198)
         mat(k,1255) = -rxt(k,419)*y(k,198)
         mat(k,1233) = -rxt(k,430)*y(k,198)
         mat(k,999) = .060_r8*rxt(k,480)*y(k,135)
         mat(k,1151) = rxt(k,327)*y(k,126) + rxt(k,328)*y(k,226)
         mat(k,1280) = rxt(k,353)*y(k,126) + rxt(k,354)*y(k,226)
         mat(k,619) = .500_r8*rxt(k,335)*y(k,226)
         mat(k,2188) = mat(k,2188) + .450_r8*rxt(k,383)*y(k,213) + .200_r8*rxt(k,387) &
                      *y(k,215) + .150_r8*rxt(k,362)*y(k,230)
         mat(k,889) = .080_r8*rxt(k,425)*y(k,135)
         mat(k,1271) = .100_r8*rxt(k,378)*y(k,135)
         mat(k,1027) = .060_r8*rxt(k,483)*y(k,135)
         mat(k,1376) = .280_r8*rxt(k,392)*y(k,135)
         mat(k,2061) = mat(k,2061) + .530_r8*rxt(k,376)*y(k,211) + rxt(k,385)*y(k,213) &
                      + rxt(k,388)*y(k,215) + rxt(k,363)*y(k,230)
         mat(k,1896) = rxt(k,327)*y(k,45) + rxt(k,353)*y(k,49) + .530_r8*rxt(k,375) &
                      *y(k,211) + rxt(k,386)*y(k,213)
         mat(k,2399) = .060_r8*rxt(k,480)*y(k,6) + .080_r8*rxt(k,425)*y(k,99) &
                      + .100_r8*rxt(k,378)*y(k,105) + .060_r8*rxt(k,483)*y(k,110) &
                      + .280_r8*rxt(k,392)*y(k,111)
         mat(k,1104) = .650_r8*rxt(k,501)*y(k,226)
         mat(k,1426) = mat(k,1426) + .530_r8*rxt(k,372)*y(k,211)
         mat(k,2243) = mat(k,2243) + .260_r8*rxt(k,373)*y(k,211) + rxt(k,382)*y(k,213) &
                      + .300_r8*rxt(k,361)*y(k,230)
         mat(k,1356) = mat(k,1356) + .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375) &
                      *y(k,126) + .530_r8*rxt(k,372)*y(k,198) + .260_r8*rxt(k,373) &
                      *y(k,199)
         mat(k,1396) = .450_r8*rxt(k,383)*y(k,90) + rxt(k,385)*y(k,124) + rxt(k,386) &
                      *y(k,126) + rxt(k,382)*y(k,199) + 4.000_r8*rxt(k,384)*y(k,213)
         mat(k,698) = .200_r8*rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124)
         mat(k,1789) = rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) + .500_r8*rxt(k,335) &
                      *y(k,51) + .650_r8*rxt(k,501)*y(k,179)
         mat(k,1217) = .150_r8*rxt(k,362)*y(k,90) + rxt(k,363)*y(k,124) &
                      + .300_r8*rxt(k,361)*y(k,199)
         mat(k,2257) = -(rxt(k,221)*y(k,59) + (4._r8*rxt(k,298) + 4._r8*rxt(k,299) &
                      ) * y(k,199) + rxt(k,300)*y(k,90) + rxt(k,301)*y(k,124) &
                      + rxt(k,320)*y(k,195) + rxt(k,331)*y(k,198) + rxt(k,348) &
                      *y(k,196) + rxt(k,361)*y(k,230) + rxt(k,373)*y(k,211) + rxt(k,382) &
                      *y(k,213) + rxt(k,406)*y(k,206) + rxt(k,411)*y(k,207) + rxt(k,420) &
                      *y(k,208) + rxt(k,431)*y(k,235) + rxt(k,485)*y(k,221) + rxt(k,490) &
                      *y(k,231) + rxt(k,495)*y(k,232))
         mat(k,1600) = -rxt(k,221)*y(k,199)
         mat(k,2205) = -rxt(k,300)*y(k,199)
         mat(k,2076) = -rxt(k,301)*y(k,199)
         mat(k,906) = -rxt(k,320)*y(k,199)
         mat(k,1435) = -rxt(k,331)*y(k,199)
         mat(k,941) = -rxt(k,348)*y(k,199)
         mat(k,1222) = -rxt(k,361)*y(k,199)
         mat(k,1364) = -rxt(k,373)*y(k,199)
         mat(k,1404) = -rxt(k,382)*y(k,199)
         mat(k,1309) = -rxt(k,406)*y(k,199)
         mat(k,1342) = -rxt(k,411)*y(k,199)
         mat(k,1263) = -rxt(k,420)*y(k,199)
         mat(k,1240) = -rxt(k,431)*y(k,199)
         mat(k,1098) = -rxt(k,485)*y(k,199)
         mat(k,1185) = -rxt(k,490)*y(k,199)
         mat(k,1068) = -rxt(k,495)*y(k,199)
         mat(k,1144) = .280_r8*rxt(k,347)*y(k,135)
         mat(k,706) = rxt(k,334)*y(k,226)
         mat(k,392) = .700_r8*rxt(k,303)*y(k,226)
         mat(k,2282) = rxt(k,215)*y(k,56) + rxt(k,271)*y(k,73) + rxt(k,310)*y(k,222) &
                      + rxt(k,304)*y(k,226)
         mat(k,1958) = rxt(k,215)*y(k,54)
         mat(k,929) = rxt(k,271)*y(k,54)
         mat(k,2205) = mat(k,2205) + .450_r8*rxt(k,332)*y(k,198) + .330_r8*rxt(k,450) &
                      *y(k,200) + .070_r8*rxt(k,456)*y(k,214)
         mat(k,893) = .050_r8*rxt(k,425)*y(k,135)
         mat(k,2076) = mat(k,2076) + rxt(k,333)*y(k,198) + .830_r8*rxt(k,451)*y(k,200) &
                      + .170_r8*rxt(k,457)*y(k,214)
         mat(k,2415) = .280_r8*rxt(k,347)*y(k,29) + .050_r8*rxt(k,425)*y(k,99)
         mat(k,1435) = mat(k,1435) + .450_r8*rxt(k,332)*y(k,90) + rxt(k,333)*y(k,124) &
                      + 4.000_r8*rxt(k,330)*y(k,198) + .900_r8*rxt(k,331)*y(k,199) &
                      + rxt(k,405)*y(k,206) + rxt(k,410)*y(k,207) + rxt(k,419) &
                      *y(k,208) + rxt(k,372)*y(k,211) + rxt(k,381)*y(k,213) &
                      + rxt(k,430)*y(k,235)
         mat(k,2257) = mat(k,2257) + .900_r8*rxt(k,331)*y(k,198)
         mat(k,783) = .330_r8*rxt(k,450)*y(k,90) + .830_r8*rxt(k,451)*y(k,124)
         mat(k,1309) = mat(k,1309) + rxt(k,405)*y(k,198)
         mat(k,1342) = mat(k,1342) + rxt(k,410)*y(k,198)
         mat(k,1263) = mat(k,1263) + rxt(k,419)*y(k,198)
         mat(k,1364) = mat(k,1364) + rxt(k,372)*y(k,198)
         mat(k,1404) = mat(k,1404) + rxt(k,381)*y(k,198)
         mat(k,921) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,1642) = rxt(k,310)*y(k,54)
         mat(k,1807) = rxt(k,334)*y(k,50) + .700_r8*rxt(k,303)*y(k,53) + rxt(k,304) &
                      *y(k,54)
         mat(k,1240) = mat(k,1240) + rxt(k,430)*y(k,198)
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
         mat(k,776) = -(rxt(k,450)*y(k,90) + rxt(k,451)*y(k,124) + rxt(k,452)*y(k,125))
         mat(k,2155) = -rxt(k,450)*y(k,200)
         mat(k,2025) = -rxt(k,451)*y(k,200)
         mat(k,1824) = -rxt(k,452)*y(k,200)
         mat(k,863) = -(rxt(k,579)*y(k,219) + rxt(k,580)*y(k,225) + rxt(k,581) &
                      *y(k,218))
         mat(k,844) = -rxt(k,579)*y(k,201)
         mat(k,852) = -rxt(k,580)*y(k,201)
         mat(k,678) = -rxt(k,581)*y(k,201)
         mat(k,572) = -((rxt(k,369) + rxt(k,370)) * y(k,124))
         mat(k,2015) = -(rxt(k,369) + rxt(k,370)) * y(k,202)
         mat(k,354) = rxt(k,368)*y(k,226)
         mat(k,1725) = rxt(k,368)*y(k,16)
         mat(k,461) = -(rxt(k,339)*y(k,134))
         mat(k,1524) = -rxt(k,339)*y(k,203)
         mat(k,2008) = .750_r8*rxt(k,337)*y(k,204)
         mat(k,794) = .750_r8*rxt(k,337)*y(k,124)
         mat(k,795) = -(rxt(k,336)*y(k,90) + rxt(k,337)*y(k,124))
         mat(k,2157) = -rxt(k,336)*y(k,204)
         mat(k,2026) = -rxt(k,337)*y(k,204)
         mat(k,549) = rxt(k,343)*y(k,226)
         mat(k,1747) = rxt(k,343)*y(k,25)
         mat(k,439) = -(rxt(k,307)*y(k,90) + rxt(k,309)*y(k,124))
         mat(k,2131) = -rxt(k,307)*y(k,205)
         mat(k,2005) = -rxt(k,309)*y(k,205)
         mat(k,2331) = rxt(k,294)*y(k,90)
         mat(k,2131) = mat(k,2131) + rxt(k,294)*y(k,42)
         mat(k,1297) = -(rxt(k,405)*y(k,198) + rxt(k,406)*y(k,199) + rxt(k,407) &
                      *y(k,90) + rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126))
         mat(k,1421) = -rxt(k,405)*y(k,206)
         mat(k,2238) = -rxt(k,406)*y(k,206)
         mat(k,2183) = -rxt(k,407)*y(k,206)
         mat(k,2056) = -rxt(k,408)*y(k,206)
         mat(k,1891) = -rxt(k,409)*y(k,206)
         mat(k,886) = .600_r8*rxt(k,426)*y(k,226)
         mat(k,1784) = .600_r8*rxt(k,426)*y(k,99)
         mat(k,1330) = -(rxt(k,410)*y(k,198) + rxt(k,411)*y(k,199) + rxt(k,412) &
                      *y(k,90) + rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126))
         mat(k,1422) = -rxt(k,410)*y(k,207)
         mat(k,2239) = -rxt(k,411)*y(k,207)
         mat(k,2184) = -rxt(k,412)*y(k,207)
         mat(k,2057) = -rxt(k,414)*y(k,207)
         mat(k,1892) = -rxt(k,415)*y(k,207)
         mat(k,887) = .400_r8*rxt(k,426)*y(k,226)
         mat(k,1785) = .400_r8*rxt(k,426)*y(k,99)
         mat(k,1251) = -(rxt(k,419)*y(k,198) + rxt(k,420)*y(k,199) + rxt(k,421) &
                      *y(k,90) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126))
         mat(k,1418) = -rxt(k,419)*y(k,208)
         mat(k,2235) = -rxt(k,420)*y(k,208)
         mat(k,2180) = -rxt(k,421)*y(k,208)
         mat(k,2053) = -rxt(k,422)*y(k,208)
         mat(k,1888) = -rxt(k,423)*y(k,208)
         mat(k,884) = rxt(k,418)*y(k,126)
         mat(k,1888) = mat(k,1888) + rxt(k,418)*y(k,99)
         mat(k,66) = -(rxt(k,543)*y(k,90) + rxt(k,544)*y(k,124))
         mat(k,2110) = -rxt(k,543)*y(k,209)
         mat(k,1994) = -rxt(k,544)*y(k,209)
         mat(k,879) = rxt(k,546)*y(k,226)
         mat(k,1656) = rxt(k,546)*y(k,99)
         mat(k,72) = -(rxt(k,547)*y(k,90) + rxt(k,548)*y(k,124))
         mat(k,2111) = -rxt(k,547)*y(k,210)
         mat(k,1995) = -rxt(k,548)*y(k,210)
         mat(k,73) = rxt(k,549)*y(k,226)
         mat(k,1657) = rxt(k,549)*y(k,104)
         mat(k,1354) = -(rxt(k,372)*y(k,198) + rxt(k,373)*y(k,199) + rxt(k,374) &
                      *y(k,90) + rxt(k,375)*y(k,126) + (rxt(k,376) + rxt(k,377) &
                      ) * y(k,124))
         mat(k,1423) = -rxt(k,372)*y(k,211)
         mat(k,2240) = -rxt(k,373)*y(k,211)
         mat(k,2185) = -rxt(k,374)*y(k,211)
         mat(k,1893) = -rxt(k,375)*y(k,211)
         mat(k,2058) = -(rxt(k,376) + rxt(k,377)) * y(k,211)
         mat(k,1269) = .500_r8*rxt(k,379)*y(k,226)
         mat(k,315) = .200_r8*rxt(k,380)*y(k,226)
         mat(k,1373) = rxt(k,393)*y(k,226)
         mat(k,1786) = .500_r8*rxt(k,379)*y(k,105) + .200_r8*rxt(k,380)*y(k,106) &
                      + rxt(k,393)*y(k,111)
         mat(k,738) = -(rxt(k,453)*y(k,90) + rxt(k,454)*y(k,124) + rxt(k,455)*y(k,125))
         mat(k,2152) = -rxt(k,453)*y(k,212)
         mat(k,2022) = -rxt(k,454)*y(k,212)
         mat(k,1823) = -rxt(k,455)*y(k,212)
         mat(k,1395) = -(rxt(k,381)*y(k,198) + rxt(k,382)*y(k,199) + rxt(k,383) &
                      *y(k,90) + 4._r8*rxt(k,384)*y(k,213) + rxt(k,385)*y(k,124) &
                      + rxt(k,386)*y(k,126) + rxt(k,394)*y(k,125))
         mat(k,1425) = -rxt(k,381)*y(k,213)
         mat(k,2242) = -rxt(k,382)*y(k,213)
         mat(k,2187) = -rxt(k,383)*y(k,213)
         mat(k,2060) = -rxt(k,385)*y(k,213)
         mat(k,1895) = -rxt(k,386)*y(k,213)
         mat(k,1834) = -rxt(k,394)*y(k,213)
         mat(k,1270) = .500_r8*rxt(k,379)*y(k,226)
         mat(k,316) = .500_r8*rxt(k,380)*y(k,226)
         mat(k,1788) = .500_r8*rxt(k,379)*y(k,105) + .500_r8*rxt(k,380)*y(k,106)
         mat(k,913) = -(rxt(k,456)*y(k,90) + rxt(k,457)*y(k,124) + rxt(k,458)*y(k,125))
         mat(k,2164) = -rxt(k,456)*y(k,214)
         mat(k,2034) = -rxt(k,457)*y(k,214)
         mat(k,1826) = -rxt(k,458)*y(k,214)
         mat(k,696) = -(rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124))
         mat(k,2148) = -rxt(k,387)*y(k,215)
         mat(k,2021) = -rxt(k,388)*y(k,215)
         mat(k,516) = rxt(k,389)*y(k,226)
         mat(k,320) = rxt(k,390)*y(k,226)
         mat(k,1738) = rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108)
         mat(k,80) = -(rxt(k,551)*y(k,90) + rxt(k,552)*y(k,124))
         mat(k,2112) = -rxt(k,551)*y(k,216)
         mat(k,1996) = -rxt(k,552)*y(k,216)
         mat(k,1012) = rxt(k,554)*y(k,226)
         mat(k,1659) = rxt(k,554)*y(k,110)
         mat(k,527) = -(rxt(k,188)*y(k,133) + rxt(k,189)*y(k,134))
         mat(k,2292) = -rxt(k,188)*y(k,217)
         mat(k,1526) = -rxt(k,189)*y(k,217)
         mat(k,2292) = mat(k,2292) + rxt(k,583)*y(k,218)
         mat(k,858) = .900_r8*rxt(k,581)*y(k,218) + .800_r8*rxt(k,579)*y(k,219)
         mat(k,673) = rxt(k,583)*y(k,133) + .900_r8*rxt(k,581)*y(k,201)
         mat(k,842) = .800_r8*rxt(k,579)*y(k,201)
         mat(k,674) = -(rxt(k,581)*y(k,201) + rxt(k,582)*y(k,134) + (rxt(k,583) &
                      + rxt(k,584)) * y(k,133))
         mat(k,859) = -rxt(k,581)*y(k,218)
         mat(k,1527) = -rxt(k,582)*y(k,218)
         mat(k,2295) = -(rxt(k,583) + rxt(k,584)) * y(k,218)
         mat(k,843) = -(rxt(k,579)*y(k,201))
         mat(k,861) = -rxt(k,579)*y(k,219)
         mat(k,966) = rxt(k,588)*y(k,225)
         mat(k,2028) = rxt(k,590)*y(k,225)
         mat(k,2301) = rxt(k,583)*y(k,218)
         mat(k,1530) = rxt(k,587)*y(k,220)
         mat(k,676) = rxt(k,583)*y(k,133)
         mat(k,502) = rxt(k,587)*y(k,134)
         mat(k,850) = rxt(k,588)*y(k,112) + rxt(k,590)*y(k,124)
         mat(k,500) = -(rxt(k,585)*y(k,133) + (rxt(k,586) + rxt(k,587)) * y(k,134))
         mat(k,2291) = -rxt(k,585)*y(k,220)
         mat(k,1525) = -(rxt(k,586) + rxt(k,587)) * y(k,220)
         mat(k,1088) = -(rxt(k,485)*y(k,199) + rxt(k,486)*y(k,90) + rxt(k,487) &
                      *y(k,124) + rxt(k,488)*y(k,126))
         mat(k,2225) = -rxt(k,485)*y(k,221)
         mat(k,2171) = -rxt(k,486)*y(k,221)
         mat(k,2042) = -rxt(k,487)*y(k,221)
         mat(k,1876) = -rxt(k,488)*y(k,221)
         mat(k,994) = rxt(k,479)*y(k,126)
         mat(k,1022) = rxt(k,482)*y(k,126)
         mat(k,1876) = mat(k,1876) + rxt(k,479)*y(k,6) + rxt(k,482)*y(k,110) &
                      + .500_r8*rxt(k,499)*y(k,178)
         mat(k,397) = rxt(k,489)*y(k,226)
         mat(k,1071) = .500_r8*rxt(k,499)*y(k,126)
         mat(k,1769) = rxt(k,489)*y(k,128)
         mat(k,1633) = -(rxt(k,153)*y(k,77) + rxt(k,154)*y(k,239) + (rxt(k,156) &
                      + rxt(k,157)) * y(k,134) + rxt(k,158)*y(k,135) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,113) + rxt(k,239)*y(k,33) + rxt(k,240) &
                      *y(k,34) + rxt(k,241)*y(k,36) + rxt(k,242)*y(k,37) + rxt(k,243) &
                      *y(k,38) + rxt(k,244)*y(k,39) + rxt(k,245)*y(k,40) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,85) + rxt(k,266)*y(k,35) + rxt(k,267) &
                      *y(k,55) + rxt(k,268)*y(k,78) + (rxt(k,269) + rxt(k,270) &
                      ) * y(k,81) + rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65) + rxt(k,289) &
                      *y(k,41) + rxt(k,290)*y(k,43) + rxt(k,291)*y(k,82) + rxt(k,292) &
                      *y(k,83) + rxt(k,293)*y(k,84) + (rxt(k,310) + rxt(k,311) &
                      + rxt(k,312)) * y(k,54) + rxt(k,313)*y(k,86))
         mat(k,1461) = -rxt(k,153)*y(k,222)
         mat(k,2433) = -rxt(k,154)*y(k,222)
         mat(k,1543) = -(rxt(k,156) + rxt(k,157)) * y(k,222)
         mat(k,2406) = -rxt(k,158)*y(k,222)
         mat(k,256) = -(rxt(k,206) + rxt(k,207)) * y(k,222)
         mat(k,100) = -rxt(k,239)*y(k,222)
         mat(k,140) = -rxt(k,240)*y(k,222)
         mat(k,111) = -rxt(k,241)*y(k,222)
         mat(k,150) = -rxt(k,242)*y(k,222)
         mat(k,115) = -rxt(k,243)*y(k,222)
         mat(k,155) = -rxt(k,244)*y(k,222)
         mat(k,119) = -rxt(k,245)*y(k,222)
         mat(k,1497) = -(rxt(k,246) + rxt(k,247)) * y(k,222)
         mat(k,146) = -rxt(k,266)*y(k,222)
         mat(k,448) = -rxt(k,267)*y(k,222)
         mat(k,108) = -rxt(k,268)*y(k,222)
         mat(k,829) = -(rxt(k,269) + rxt(k,270)) * y(k,222)
         mat(k,244) = -rxt(k,275)*y(k,222)
         mat(k,226) = -rxt(k,276)*y(k,222)
         mat(k,473) = -rxt(k,289)*y(k,222)
         mat(k,605) = -rxt(k,290)*y(k,222)
         mat(k,221) = -rxt(k,291)*y(k,222)
         mat(k,251) = -rxt(k,292)*y(k,222)
         mat(k,274) = -rxt(k,293)*y(k,222)
         mat(k,2273) = -(rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,222)
         mat(k,181) = -rxt(k,313)*y(k,222)
         mat(k,1543) = mat(k,1543) + rxt(k,189)*y(k,217)
         mat(k,868) = .850_r8*rxt(k,580)*y(k,225)
         mat(k,530) = rxt(k,189)*y(k,134)
         mat(k,855) = .850_r8*rxt(k,580)*y(k,201)
         mat(k,175) = -(rxt(k,160)*y(k,133) + rxt(k,161)*y(k,134))
         mat(k,2288) = -rxt(k,160)*y(k,223)
         mat(k,1522) = -rxt(k,161)*y(k,223)
         mat(k,1440) = rxt(k,162)*y(k,224)
         mat(k,2288) = mat(k,2288) + rxt(k,164)*y(k,224)
         mat(k,1522) = mat(k,1522) + rxt(k,165)*y(k,224)
         mat(k,2362) = rxt(k,166)*y(k,224)
         mat(k,177) = rxt(k,162)*y(k,63) + rxt(k,164)*y(k,133) + rxt(k,165)*y(k,134) &
                      + rxt(k,166)*y(k,135)
         mat(k,178) = -(rxt(k,162)*y(k,63) + rxt(k,164)*y(k,133) + rxt(k,165)*y(k,134) &
                      + rxt(k,166)*y(k,135))
         mat(k,1441) = -rxt(k,162)*y(k,224)
         mat(k,2289) = -rxt(k,164)*y(k,224)
         mat(k,1523) = -rxt(k,165)*y(k,224)
         mat(k,2363) = -rxt(k,166)*y(k,224)
         mat(k,1523) = mat(k,1523) + rxt(k,156)*y(k,222)
         mat(k,1614) = rxt(k,156)*y(k,134)
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
         mat(k,851) = -(rxt(k,580)*y(k,201) + rxt(k,588)*y(k,112) + rxt(k,590) &
                      *y(k,124))
         mat(k,862) = -rxt(k,580)*y(k,225)
         mat(k,967) = -rxt(k,588)*y(k,225)
         mat(k,2029) = -rxt(k,590)*y(k,225)
         mat(k,1444) = rxt(k,591)*y(k,227)
         mat(k,1531) = rxt(k,582)*y(k,218) + rxt(k,586)*y(k,220) + rxt(k,593)*y(k,227)
         mat(k,677) = rxt(k,582)*y(k,134)
         mat(k,503) = rxt(k,586)*y(k,134)
         mat(k,805) = rxt(k,591)*y(k,63) + rxt(k,593)*y(k,134)
         mat(k,1799) = -(rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,181)*y(k,90) &
                      + rxt(k,182)*y(k,133) + rxt(k,183)*y(k,135) + (4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185)) * y(k,226) + rxt(k,187)*y(k,91) + rxt(k,201) &
                      *y(k,126) + rxt(k,202)*y(k,112) + rxt(k,210)*y(k,125) + rxt(k,211) &
                      *y(k,89) + rxt(k,230)*y(k,60) + (rxt(k,232) + rxt(k,233) &
                      ) * y(k,59) + rxt(k,235)*y(k,85) + rxt(k,238)*y(k,93) + rxt(k,262) &
                      *y(k,19) + rxt(k,264)*y(k,81) + rxt(k,278)*y(k,41) + rxt(k,280) &
                      *y(k,43) + rxt(k,281)*y(k,44) + rxt(k,283)*y(k,46) + rxt(k,285) &
                      *y(k,55) + rxt(k,286)*y(k,82) + rxt(k,287)*y(k,83) + rxt(k,288) &
                      *y(k,84) + rxt(k,297)*y(k,42) + rxt(k,302)*y(k,52) + rxt(k,303) &
                      *y(k,53) + rxt(k,304)*y(k,54) + rxt(k,305)*y(k,86) + rxt(k,306) &
                      *y(k,87) + rxt(k,314)*y(k,62) + rxt(k,316)*y(k,24) + rxt(k,323) &
                      *y(k,26) + rxt(k,324)*y(k,27) + rxt(k,326)*y(k,28) + rxt(k,328) &
                      *y(k,45) + rxt(k,329)*y(k,47) + rxt(k,334)*y(k,50) + rxt(k,335) &
                      *y(k,51) + rxt(k,340)*y(k,74) + rxt(k,341)*y(k,75) + rxt(k,342) &
                      *y(k,140) + rxt(k,343)*y(k,25) + rxt(k,351)*y(k,30) + rxt(k,352) &
                      *y(k,31) + rxt(k,354)*y(k,49) + rxt(k,356)*y(k,96) + rxt(k,357) &
                      *y(k,127) + rxt(k,360)*y(k,147) + rxt(k,364)*y(k,148) + rxt(k,365) &
                      *y(k,29) + rxt(k,366)*y(k,48) + rxt(k,368)*y(k,16) + rxt(k,371) &
                      *y(k,94) + rxt(k,379)*y(k,105) + rxt(k,380)*y(k,106) + rxt(k,389) &
                      *y(k,107) + rxt(k,390)*y(k,108) + rxt(k,391)*y(k,109) + rxt(k,393) &
                      *y(k,111) + rxt(k,396)*y(k,1) + rxt(k,400)*y(k,2) + rxt(k,401) &
                      *y(k,15) + rxt(k,402)*y(k,95) + rxt(k,403)*y(k,97) + rxt(k,404) &
                      *y(k,98) + rxt(k,416)*y(k,100) + rxt(k,417)*y(k,101) + rxt(k,424) &
                      *y(k,102) + rxt(k,426)*y(k,99) + rxt(k,427)*y(k,103) + rxt(k,428) &
                      *y(k,115) + rxt(k,429)*y(k,116) + rxt(k,435)*y(k,182) + rxt(k,438) &
                      *y(k,7) + rxt(k,441)*y(k,8) + rxt(k,442)*y(k,22) + rxt(k,444) &
                      *y(k,23) + rxt(k,448)*y(k,32) + rxt(k,449)*y(k,66) + rxt(k,461) &
                      *y(k,143) + rxt(k,464)*y(k,144) + rxt(k,468)*y(k,180) + rxt(k,469) &
                      *y(k,181) + rxt(k,471)*y(k,183) + rxt(k,474)*y(k,184) + rxt(k,477) &
                      *y(k,185) + rxt(k,478)*y(k,186) + rxt(k,481)*y(k,6) + rxt(k,484) &
                      *y(k,110) + rxt(k,489)*y(k,128) + rxt(k,493)*y(k,175) + rxt(k,494) &
                      *y(k,176) + rxt(k,498)*y(k,177) + rxt(k,500)*y(k,178) + rxt(k,501) &
                      *y(k,179) + (rxt(k,503) + rxt(k,517)) * y(k,67) + rxt(k,505) &
                      *y(k,138) + rxt(k,507)*y(k,152) + rxt(k,511)*y(k,149) + rxt(k,516) &
                      *y(k,151) + rxt(k,519)*y(k,120))
         mat(k,1462) = -rxt(k,179)*y(k,226)
         mat(k,613) = -rxt(k,180)*y(k,226)
         mat(k,2197) = -rxt(k,181)*y(k,226)
         mat(k,2317) = -rxt(k,182)*y(k,226)
         mat(k,2407) = -rxt(k,183)*y(k,226)
         mat(k,495) = -rxt(k,187)*y(k,226)
         mat(k,1904) = -rxt(k,201)*y(k,226)
         mat(k,974) = -rxt(k,202)*y(k,226)
         mat(k,1844) = -rxt(k,210)*y(k,226)
         mat(k,1973) = -rxt(k,211)*y(k,226)
         mat(k,956) = -rxt(k,230)*y(k,226)
         mat(k,1592) = -(rxt(k,232) + rxt(k,233)) * y(k,226)
         mat(k,1498) = -rxt(k,235)*y(k,226)
         mat(k,838) = -rxt(k,238)*y(k,226)
         mat(k,1566) = -rxt(k,262)*y(k,226)
         mat(k,830) = -rxt(k,264)*y(k,226)
         mat(k,474) = -rxt(k,278)*y(k,226)
         mat(k,606) = -rxt(k,280)*y(k,226)
         mat(k,122) = -rxt(k,281)*y(k,226)
         mat(k,372) = -rxt(k,283)*y(k,226)
         mat(k,449) = -rxt(k,285)*y(k,226)
         mat(k,222) = -rxt(k,286)*y(k,226)
         mat(k,252) = -rxt(k,287)*y(k,226)
         mat(k,275) = -rxt(k,288)*y(k,226)
         mat(k,2343) = -rxt(k,297)*y(k,226)
         mat(k,812) = -rxt(k,302)*y(k,226)
         mat(k,390) = -rxt(k,303)*y(k,226)
         mat(k,2274) = -rxt(k,304)*y(k,226)
         mat(k,182) = -rxt(k,305)*y(k,226)
         mat(k,910) = -rxt(k,306)*y(k,226)
         mat(k,1160) = -rxt(k,314)*y(k,226)
         mat(k,294) = -rxt(k,316)*y(k,226)
         mat(k,265) = -rxt(k,323)*y(k,226)
         mat(k,351) = -rxt(k,324)*y(k,226)
         mat(k,299) = -rxt(k,326)*y(k,226)
         mat(k,1152) = -rxt(k,328)*y(k,226)
         mat(k,103) = -rxt(k,329)*y(k,226)
         mat(k,705) = -rxt(k,334)*y(k,226)
         mat(k,621) = -rxt(k,335)*y(k,226)
         mat(k,1166) = -rxt(k,340)*y(k,226)
         mat(k,1055) = -rxt(k,341)*y(k,226)
         mat(k,535) = -rxt(k,342)*y(k,226)
         mat(k,552) = -rxt(k,343)*y(k,226)
         mat(k,416) = -rxt(k,351)*y(k,226)
         mat(k,305) = -rxt(k,352)*y(k,226)
         mat(k,1282) = -rxt(k,354)*y(k,226)
         mat(k,1209) = -rxt(k,356)*y(k,226)
         mat(k,874) = -rxt(k,357)*y(k,226)
         mat(k,544) = -rxt(k,360)*y(k,226)
         mat(k,410) = -rxt(k,364)*y(k,226)
         mat(k,1139) = -rxt(k,365)*y(k,226)
         mat(k,1081) = -rxt(k,366)*y(k,226)
         mat(k,357) = -rxt(k,368)*y(k,226)
         mat(k,1198) = -rxt(k,371)*y(k,226)
         mat(k,1273) = -rxt(k,379)*y(k,226)
         mat(k,317) = -rxt(k,380)*y(k,226)
         mat(k,519) = -rxt(k,389)*y(k,226)
         mat(k,323) = -rxt(k,390)*y(k,226)
         mat(k,584) = -rxt(k,391)*y(k,226)
         mat(k,1379) = -rxt(k,393)*y(k,226)
         mat(k,646) = -rxt(k,396)*y(k,226)
         mat(k,692) = -rxt(k,400)*y(k,226)
         mat(k,235) = -rxt(k,401)*y(k,226)
         mat(k,231) = -rxt(k,402)*y(k,226)
         mat(k,326) = -rxt(k,403)*y(k,226)
         mat(k,133) = -rxt(k,404)*y(k,226)
         mat(k,593) = -rxt(k,416)*y(k,226)
         mat(k,561) = -rxt(k,417)*y(k,226)
         mat(k,422) = -rxt(k,424)*y(k,226)
         mat(k,890) = -rxt(k,426)*y(k,226)
         mat(k,735) = -rxt(k,427)*y(k,226)
         mat(k,380) = -rxt(k,428)*y(k,226)
         mat(k,1120) = -rxt(k,429)*y(k,226)
         mat(k,203) = -rxt(k,435)*y(k,226)
         mat(k,162) = -rxt(k,438)*y(k,226)
         mat(k,405) = -rxt(k,441)*y(k,226)
         mat(k,247) = -rxt(k,442)*y(k,226)
         mat(k,346) = -rxt(k,444)*y(k,226)
         mat(k,270) = -rxt(k,448)*y(k,226)
         mat(k,195) = -rxt(k,449)*y(k,226)
         mat(k,171) = -rxt(k,461)*y(k,226)
         mat(k,340) = -rxt(k,464)*y(k,226)
         mat(k,671) = -rxt(k,468)*y(k,226)
         mat(k,190) = -rxt(k,469)*y(k,226)
         mat(k,212) = -rxt(k,471)*y(k,226)
         mat(k,728) = -rxt(k,474)*y(k,226)
         mat(k,217) = -rxt(k,477)*y(k,226)
         mat(k,429) = -rxt(k,478)*y(k,226)
         mat(k,1001) = -rxt(k,481)*y(k,226)
         mat(k,1029) = -rxt(k,484)*y(k,226)
         mat(k,398) = -rxt(k,489)*y(k,226)
         mat(k,657) = -rxt(k,493)*y(k,226)
         mat(k,627) = -rxt(k,494)*y(k,226)
         mat(k,481) = -rxt(k,498)*y(k,226)
         mat(k,1075) = -rxt(k,500)*y(k,226)
         mat(k,1106) = -rxt(k,501)*y(k,226)
         mat(k,310) = -(rxt(k,503) + rxt(k,517)) * y(k,226)
         mat(k,366) = -rxt(k,505)*y(k,226)
         mat(k,947) = -rxt(k,507)*y(k,226)
         mat(k,711) = -rxt(k,511)*y(k,226)
         mat(k,1479) = -rxt(k,516)*y(k,226)
         mat(k,97) = -rxt(k,519)*y(k,226)
         mat(k,1001) = mat(k,1001) + .630_r8*rxt(k,480)*y(k,135)
         mat(k,294) = mat(k,294) + .650_r8*rxt(k,316)*y(k,226)
         mat(k,552) = mat(k,552) + .130_r8*rxt(k,318)*y(k,135)
         mat(k,351) = mat(k,351) + .500_r8*rxt(k,324)*y(k,226)
         mat(k,1139) = mat(k,1139) + .360_r8*rxt(k,347)*y(k,135)
         mat(k,2343) = mat(k,2343) + rxt(k,296)*y(k,133)
         mat(k,390) = mat(k,390) + .300_r8*rxt(k,303)*y(k,226)
         mat(k,2274) = mat(k,2274) + rxt(k,310)*y(k,222)
         mat(k,1950) = rxt(k,219)*y(k,90)
         mat(k,925) = rxt(k,273)*y(k,239)
         mat(k,2089) = 2.000_r8*rxt(k,173)*y(k,90) + rxt(k,178)*y(k,135)
         mat(k,1462) = mat(k,1462) + rxt(k,170)*y(k,133) + rxt(k,153)*y(k,222)
         mat(k,613) = mat(k,613) + rxt(k,171)*y(k,133)
         mat(k,830) = mat(k,830) + rxt(k,263)*y(k,133) + rxt(k,269)*y(k,222)
         mat(k,1498) = mat(k,1498) + rxt(k,234)*y(k,133) + rxt(k,246)*y(k,222)
         mat(k,182) = mat(k,182) + rxt(k,313)*y(k,222)
         mat(k,2197) = mat(k,2197) + rxt(k,219)*y(k,56) + 2.000_r8*rxt(k,173)*y(k,76) &
                      + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126) + rxt(k,176) &
                      *y(k,133) + rxt(k,177)*y(k,135) + .400_r8*rxt(k,436)*y(k,189) &
                      + .450_r8*rxt(k,332)*y(k,198) + .400_r8*rxt(k,450)*y(k,200) &
                      + .450_r8*rxt(k,383)*y(k,213) + .400_r8*rxt(k,456)*y(k,214) &
                      + .200_r8*rxt(k,387)*y(k,215) + .150_r8*rxt(k,362)*y(k,230)
         mat(k,789) = rxt(k,265)*y(k,133)
         mat(k,838) = mat(k,838) + rxt(k,237)*y(k,133)
         mat(k,890) = mat(k,890) + .320_r8*rxt(k,425)*y(k,135)
         mat(k,735) = mat(k,735) + .600_r8*rxt(k,427)*y(k,226)
         mat(k,1273) = mat(k,1273) + .240_r8*rxt(k,378)*y(k,135)
         mat(k,317) = mat(k,317) + .100_r8*rxt(k,380)*y(k,226)
         mat(k,1029) = mat(k,1029) + .630_r8*rxt(k,483)*y(k,135)
         mat(k,1379) = mat(k,1379) + .360_r8*rxt(k,392)*y(k,135)
         mat(k,2068) = rxt(k,203)*y(k,90)
         mat(k,1904) = mat(k,1904) + rxt(k,198)*y(k,90)
         mat(k,2317) = mat(k,2317) + rxt(k,296)*y(k,42) + rxt(k,170)*y(k,77) &
                      + rxt(k,171)*y(k,79) + rxt(k,263)*y(k,81) + rxt(k,234)*y(k,85) &
                      + rxt(k,176)*y(k,90) + rxt(k,265)*y(k,92) + rxt(k,237)*y(k,93)
         mat(k,2407) = mat(k,2407) + .630_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,318) &
                      *y(k,25) + .360_r8*rxt(k,347)*y(k,29) + rxt(k,178)*y(k,76) &
                      + rxt(k,177)*y(k,90) + .320_r8*rxt(k,425)*y(k,99) &
                      + .240_r8*rxt(k,378)*y(k,105) + .630_r8*rxt(k,483)*y(k,110) &
                      + .360_r8*rxt(k,392)*y(k,111)
         mat(k,544) = mat(k,544) + .500_r8*rxt(k,360)*y(k,226)
         mat(k,203) = mat(k,203) + .500_r8*rxt(k,435)*y(k,226)
         mat(k,523) = .400_r8*rxt(k,436)*y(k,90)
         mat(k,1429) = .450_r8*rxt(k,332)*y(k,90)
         mat(k,779) = .400_r8*rxt(k,450)*y(k,90)
         mat(k,1398) = .450_r8*rxt(k,383)*y(k,90)
         mat(k,917) = .400_r8*rxt(k,456)*y(k,90)
         mat(k,699) = .200_r8*rxt(k,387)*y(k,90)
         mat(k,1634) = rxt(k,310)*y(k,54) + rxt(k,153)*y(k,77) + rxt(k,269)*y(k,81) &
                      + rxt(k,246)*y(k,85) + rxt(k,313)*y(k,86) + 2.000_r8*rxt(k,154) &
                      *y(k,239)
         mat(k,1799) = mat(k,1799) + .650_r8*rxt(k,316)*y(k,24) + .500_r8*rxt(k,324) &
                      *y(k,27) + .300_r8*rxt(k,303)*y(k,53) + .600_r8*rxt(k,427) &
                      *y(k,103) + .100_r8*rxt(k,380)*y(k,106) + .500_r8*rxt(k,360) &
                      *y(k,147) + .500_r8*rxt(k,435)*y(k,182)
         mat(k,1218) = .150_r8*rxt(k,362)*y(k,90)
         mat(k,2434) = rxt(k,273)*y(k,73) + 2.000_r8*rxt(k,154)*y(k,222)
      end do
      end subroutine nlnmat10
      subroutine nlnmat11( avec_len, mat, y, rxt )
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
         mat(k,803) = -(rxt(k,591)*y(k,63) + rxt(k,593)*y(k,134))
         mat(k,1442) = -rxt(k,591)*y(k,227)
         mat(k,1529) = -rxt(k,593)*y(k,227)
         mat(k,2298) = rxt(k,584)*y(k,218) + rxt(k,585)*y(k,220)
         mat(k,675) = rxt(k,584)*y(k,133)
         mat(k,501) = rxt(k,585)*y(k,133)
         mat(k,452) = -(rxt(k,459)*y(k,90) + rxt(k,460)*y(k,124))
         mat(k,2132) = -rxt(k,459)*y(k,228)
         mat(k,2006) = -rxt(k,460)*y(k,228)
         mat(k,193) = .200_r8*rxt(k,449)*y(k,226)
         mat(k,169) = .140_r8*rxt(k,461)*y(k,226)
         mat(k,338) = rxt(k,464)*y(k,226)
         mat(k,1711) = .200_r8*rxt(k,449)*y(k,66) + .140_r8*rxt(k,461)*y(k,143) &
                      + rxt(k,464)*y(k,144)
         mat(k,816) = -(rxt(k,358)*y(k,90) + rxt(k,359)*y(k,124))
         mat(k,2158) = -rxt(k,358)*y(k,229)
         mat(k,2027) = -rxt(k,359)*y(k,229)
         mat(k,1128) = rxt(k,365)*y(k,226)
         mat(k,541) = .500_r8*rxt(k,360)*y(k,226)
         mat(k,1749) = rxt(k,365)*y(k,29) + .500_r8*rxt(k,360)*y(k,147)
         mat(k,1215) = -(rxt(k,361)*y(k,199) + rxt(k,362)*y(k,90) + rxt(k,363) &
                      *y(k,124))
         mat(k,2233) = -rxt(k,361)*y(k,230)
         mat(k,2178) = -rxt(k,362)*y(k,230)
         mat(k,2051) = -rxt(k,363)*y(k,230)
         mat(k,997) = .060_r8*rxt(k,480)*y(k,135)
         mat(k,1079) = rxt(k,366)*y(k,226)
         mat(k,1025) = .060_r8*rxt(k,483)*y(k,135)
         mat(k,2390) = .060_r8*rxt(k,480)*y(k,6) + .060_r8*rxt(k,483)*y(k,110)
         mat(k,408) = rxt(k,364)*y(k,226)
         mat(k,1103) = .150_r8*rxt(k,501)*y(k,226)
         mat(k,1779) = rxt(k,366)*y(k,48) + rxt(k,364)*y(k,148) + .150_r8*rxt(k,501) &
                      *y(k,179)
         mat(k,1176) = -(rxt(k,490)*y(k,199) + rxt(k,491)*y(k,90) + rxt(k,492) &
                      *y(k,124))
         mat(k,2231) = -rxt(k,490)*y(k,231)
         mat(k,2176) = -rxt(k,491)*y(k,231)
         mat(k,2048) = -rxt(k,492)*y(k,231)
         mat(k,1883) = .500_r8*rxt(k,499)*y(k,178)
         mat(k,655) = rxt(k,493)*y(k,226)
         mat(k,1074) = .500_r8*rxt(k,499)*y(k,126) + rxt(k,500)*y(k,226)
         mat(k,1776) = rxt(k,493)*y(k,175) + rxt(k,500)*y(k,178)
         mat(k,1060) = -(rxt(k,495)*y(k,199) + rxt(k,496)*y(k,90) + rxt(k,497) &
                      *y(k,124))
         mat(k,2222) = -rxt(k,495)*y(k,232)
         mat(k,2168) = -rxt(k,496)*y(k,232)
         mat(k,2039) = -rxt(k,497)*y(k,232)
         mat(k,991) = rxt(k,481)*y(k,226)
         mat(k,1019) = rxt(k,484)*y(k,226)
         mat(k,478) = rxt(k,498)*y(k,226)
         mat(k,1766) = rxt(k,481)*y(k,6) + rxt(k,484)*y(k,110) + rxt(k,498)*y(k,177)
         mat(k,749) = -(rxt(k,466)*y(k,90) + rxt(k,467)*y(k,124))
         mat(k,2153) = -rxt(k,466)*y(k,233)
         mat(k,2023) = -rxt(k,467)*y(k,233)
         mat(k,665) = rxt(k,468)*y(k,226)
         mat(k,189) = .650_r8*rxt(k,469)*y(k,226)
         mat(k,1744) = rxt(k,468)*y(k,180) + .650_r8*rxt(k,469)*y(k,181)
         mat(k,86) = -(rxt(k,557)*y(k,90) + rxt(k,558)*y(k,124))
         mat(k,2113) = -rxt(k,557)*y(k,234)
         mat(k,1997) = -rxt(k,558)*y(k,234)
         mat(k,184) = rxt(k,556)*y(k,226)
         mat(k,1660) = rxt(k,556)*y(k,181)
         mat(k,1231) = -(rxt(k,430)*y(k,198) + rxt(k,431)*y(k,199) + rxt(k,432) &
                      *y(k,90) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126))
         mat(k,1417) = -rxt(k,430)*y(k,235)
         mat(k,2234) = -rxt(k,431)*y(k,235)
         mat(k,2179) = -rxt(k,432)*y(k,235)
         mat(k,2052) = -rxt(k,433)*y(k,235)
         mat(k,1887) = -rxt(k,434)*y(k,235)
         mat(k,230) = rxt(k,402)*y(k,226)
         mat(k,325) = rxt(k,403)*y(k,226)
         mat(k,132) = rxt(k,404)*y(k,226)
         mat(k,732) = .400_r8*rxt(k,427)*y(k,226)
         mat(k,202) = .500_r8*rxt(k,435)*y(k,226)
         mat(k,1780) = rxt(k,402)*y(k,95) + rxt(k,403)*y(k,97) + rxt(k,404)*y(k,98) &
                      + .400_r8*rxt(k,427)*y(k,103) + .500_r8*rxt(k,435)*y(k,182)
         mat(k,765) = -(rxt(k,472)*y(k,90) + rxt(k,473)*y(k,124))
         mat(k,2154) = -rxt(k,472)*y(k,236)
         mat(k,2024) = -rxt(k,473)*y(k,236)
         mat(k,209) = .560_r8*rxt(k,471)*y(k,226)
         mat(k,721) = rxt(k,474)*y(k,226)
         mat(k,1745) = .560_r8*rxt(k,471)*y(k,183) + rxt(k,474)*y(k,184)
         mat(k,92) = -(rxt(k,560)*y(k,90) + rxt(k,561)*y(k,124))
         mat(k,2114) = -rxt(k,560)*y(k,237)
         mat(k,1998) = -rxt(k,561)*y(k,237)
         mat(k,204) = rxt(k,559)*y(k,226)
         mat(k,1661) = rxt(k,559)*y(k,183)
         mat(k,508) = -(rxt(k,475)*y(k,90) + rxt(k,476)*y(k,124))
         mat(k,2139) = -rxt(k,475)*y(k,238)
         mat(k,2011) = -rxt(k,476)*y(k,238)
         mat(k,216) = .300_r8*rxt(k,477)*y(k,226)
         mat(k,426) = rxt(k,478)*y(k,226)
         mat(k,1718) = .300_r8*rxt(k,477)*y(k,185) + rxt(k,478)*y(k,186)
         mat(k,2447) = -(rxt(k,154)*y(k,222) + rxt(k,273)*y(k,73) + rxt(k,518) &
                      *y(k,153))
         mat(k,1647) = -rxt(k,154)*y(k,239)
         mat(k,931) = -rxt(k,273)*y(k,239)
         mat(k,262) = -rxt(k,518)*y(k,239)
         mat(k,301) = rxt(k,326)*y(k,226)
         mat(k,418) = rxt(k,351)*y(k,226)
         mat(k,307) = rxt(k,352)*y(k,226)
         mat(k,476) = rxt(k,278)*y(k,226)
         mat(k,2356) = rxt(k,297)*y(k,226)
         mat(k,610) = rxt(k,280)*y(k,226)
         mat(k,124) = rxt(k,281)*y(k,226)
         mat(k,1157) = rxt(k,328)*y(k,226)
         mat(k,376) = rxt(k,283)*y(k,226)
         mat(k,1083) = rxt(k,366)*y(k,226)
         mat(k,1286) = rxt(k,354)*y(k,226)
         mat(k,707) = rxt(k,334)*y(k,226)
         mat(k,624) = rxt(k,335)*y(k,226)
         mat(k,394) = rxt(k,303)*y(k,226)
         mat(k,2287) = rxt(k,304)*y(k,226)
         mat(k,2102) = rxt(k,174)*y(k,90)
         mat(k,1470) = rxt(k,179)*y(k,226)
         mat(k,617) = rxt(k,180)*y(k,226)
         mat(k,833) = rxt(k,264)*y(k,226)
         mat(k,277) = rxt(k,288)*y(k,226)
         mat(k,1505) = (rxt(k,570)+rxt(k,575))*y(k,92) + (rxt(k,563)+rxt(k,569) &
                       +rxt(k,574))*y(k,93) + rxt(k,235)*y(k,226)
         mat(k,912) = rxt(k,306)*y(k,226)
         mat(k,1986) = rxt(k,211)*y(k,226)
         mat(k,2210) = rxt(k,174)*y(k,76) + rxt(k,181)*y(k,226)
         mat(k,499) = rxt(k,187)*y(k,226)
         mat(k,792) = (rxt(k,570)+rxt(k,575))*y(k,85)
         mat(k,841) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,85) + rxt(k,238)*y(k,226)
         mat(k,1277) = .500_r8*rxt(k,379)*y(k,226)
         mat(k,98) = rxt(k,519)*y(k,226)
         mat(k,547) = rxt(k,360)*y(k,226)
         mat(k,412) = rxt(k,364)*y(k,226)
         mat(k,1812) = rxt(k,326)*y(k,28) + rxt(k,351)*y(k,30) + rxt(k,352)*y(k,31) &
                      + rxt(k,278)*y(k,41) + rxt(k,297)*y(k,42) + rxt(k,280)*y(k,43) &
                      + rxt(k,281)*y(k,44) + rxt(k,328)*y(k,45) + rxt(k,283)*y(k,46) &
                      + rxt(k,366)*y(k,48) + rxt(k,354)*y(k,49) + rxt(k,334)*y(k,50) &
                      + rxt(k,335)*y(k,51) + rxt(k,303)*y(k,53) + rxt(k,304)*y(k,54) &
                      + rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,264)*y(k,81) &
                      + rxt(k,288)*y(k,84) + rxt(k,235)*y(k,85) + rxt(k,306)*y(k,87) &
                      + rxt(k,211)*y(k,89) + rxt(k,181)*y(k,90) + rxt(k,187)*y(k,91) &
                      + rxt(k,238)*y(k,93) + .500_r8*rxt(k,379)*y(k,105) + rxt(k,519) &
                      *y(k,120) + rxt(k,360)*y(k,147) + rxt(k,364)*y(k,148) &
                      + 2.000_r8*rxt(k,184)*y(k,226)
      end do
      end subroutine nlnmat11
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
         mat(k, 34) = lmat(k, 34)
         mat(k, 35) = lmat(k, 35)
         mat(k, 36) = lmat(k, 36)
         mat(k, 37) = lmat(k, 37)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = lmat(k, 42)
         mat(k, 48) = mat(k, 48) + lmat(k, 48)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 60) = mat(k, 60) + lmat(k, 60)
         mat(k, 66) = mat(k, 66) + lmat(k, 66)
         mat(k, 72) = mat(k, 72) + lmat(k, 72)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 92) = mat(k, 92) + lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 99) = mat(k, 99) + lmat(k, 99)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 106) = mat(k, 106) + lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = lmat(k, 136)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 138) = mat(k, 138) + lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 164) = lmat(k, 164)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = mat(k, 176) + lmat(k, 176)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 192) = mat(k, 192) + lmat(k, 192)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = mat(k, 201) + lmat(k, 201)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 229) = lmat(k, 229)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = lmat(k, 237)
         mat(k, 238) = lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 260) = lmat(k, 260)
         mat(k, 261) = lmat(k, 261)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = lmat(k, 285)
         mat(k, 286) = lmat(k, 286)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = lmat(k, 288)
         mat(k, 289) = lmat(k, 289)
         mat(k, 290) = mat(k, 290) + lmat(k, 290)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = mat(k, 323) + lmat(k, 323)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 327) = lmat(k, 327)
         mat(k, 328) = lmat(k, 328)
         mat(k, 329) = lmat(k, 329)
         mat(k, 330) = lmat(k, 330)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 352) = lmat(k, 352)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 362) = lmat(k, 362)
         mat(k, 364) = mat(k, 364) + lmat(k, 364)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 375) = lmat(k, 375)
         mat(k, 377) = mat(k, 377) + lmat(k, 377)
         mat(k, 381) = lmat(k, 381)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = lmat(k, 387)
         mat(k, 388) = lmat(k, 388)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = mat(k, 390) + lmat(k, 390)
         mat(k, 391) = lmat(k, 391)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 399) = lmat(k, 399)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 404) = lmat(k, 404)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 406) = lmat(k, 406)
         mat(k, 407) = mat(k, 407) + lmat(k, 407)
         mat(k, 409) = lmat(k, 409)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 411) = lmat(k, 411)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 415) = lmat(k, 415)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 427) = lmat(k, 427)
         mat(k, 428) = lmat(k, 428)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 430) = lmat(k, 430)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 445) = mat(k, 445) + lmat(k, 445)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 452) = mat(k, 452) + lmat(k, 452)
         mat(k, 458) = lmat(k, 458)
         mat(k, 459) = lmat(k, 459)
         mat(k, 460) = lmat(k, 460)
         mat(k, 461) = mat(k, 461) + lmat(k, 461)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 465) = lmat(k, 465)
         mat(k, 466) = mat(k, 466) + lmat(k, 466)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 477) = mat(k, 477) + lmat(k, 477)
         mat(k, 479) = lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 483) = lmat(k, 483)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 497) = lmat(k, 497)
         mat(k, 498) = lmat(k, 498)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 517) = lmat(k, 517)
         mat(k, 518) = lmat(k, 518)
         mat(k, 521) = mat(k, 521) + lmat(k, 521)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 532) = mat(k, 532) + lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = lmat(k, 534)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 538) = lmat(k, 538)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 542) = lmat(k, 542)
         mat(k, 544) = mat(k, 544) + lmat(k, 544)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = lmat(k, 546)
         mat(k, 548) = mat(k, 548) + lmat(k, 548)
         mat(k, 556) = mat(k, 556) + lmat(k, 556)
         mat(k, 562) = lmat(k, 562)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = lmat(k, 566)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = lmat(k, 570)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 582) = lmat(k, 582)
         mat(k, 585) = lmat(k, 585)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 595) = lmat(k, 595)
         mat(k, 598) = lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = lmat(k, 600)
         mat(k, 601) = lmat(k, 601)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 609) = lmat(k, 609)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 621) = mat(k, 621) + lmat(k, 621)
         mat(k, 622) = lmat(k, 622)
         mat(k, 625) = mat(k, 625) + lmat(k, 625)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 629) = lmat(k, 629)
         mat(k, 630) = lmat(k, 630)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 639) = lmat(k, 639)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 644) = mat(k, 644) + lmat(k, 644)
         mat(k, 647) = mat(k, 647) + lmat(k, 647)
         mat(k, 648) = lmat(k, 648)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 656) = lmat(k, 656)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 658) = lmat(k, 658)
         mat(k, 659) = lmat(k, 659)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = lmat(k, 662)
         mat(k, 663) = mat(k, 663) + lmat(k, 663)
         mat(k, 668) = lmat(k, 668)
         mat(k, 670) = lmat(k, 670)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 672) = lmat(k, 672)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 684) = lmat(k, 684)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 689) = lmat(k, 689)
         mat(k, 690) = lmat(k, 690)
         mat(k, 692) = mat(k, 692) + lmat(k, 692)
         mat(k, 693) = lmat(k, 693)
         mat(k, 694) = lmat(k, 694)
         mat(k, 696) = mat(k, 696) + lmat(k, 696)
         mat(k, 703) = mat(k, 703) + lmat(k, 703)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 715) = lmat(k, 715)
         mat(k, 716) = lmat(k, 716)
         mat(k, 717) = lmat(k, 717)
         mat(k, 718) = lmat(k, 718)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 724) = lmat(k, 724)
         mat(k, 726) = lmat(k, 726)
         mat(k, 728) = mat(k, 728) + lmat(k, 728)
         mat(k, 729) = lmat(k, 729)
         mat(k, 731) = mat(k, 731) + lmat(k, 731)
         mat(k, 733) = lmat(k, 733)
         mat(k, 734) = lmat(k, 734)
         mat(k, 735) = mat(k, 735) + lmat(k, 735)
         mat(k, 736) = lmat(k, 736)
         mat(k, 737) = lmat(k, 737)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 776) = mat(k, 776) + lmat(k, 776)
         mat(k, 785) = mat(k, 785) + lmat(k, 785)
         mat(k, 787) = lmat(k, 787)
         mat(k, 789) = mat(k, 789) + lmat(k, 789)
         mat(k, 795) = mat(k, 795) + lmat(k, 795)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 804) = lmat(k, 804)
         mat(k, 806) = lmat(k, 806)
         mat(k, 811) = mat(k, 811) + lmat(k, 811)
         mat(k, 816) = mat(k, 816) + lmat(k, 816)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 827) = mat(k, 827) + lmat(k, 827)
         mat(k, 831) = mat(k, 831) + lmat(k, 831)
         mat(k, 835) = mat(k, 835) + lmat(k, 835)
         mat(k, 838) = mat(k, 838) + lmat(k, 838)
         mat(k, 839) = mat(k, 839) + lmat(k, 839)
         mat(k, 843) = mat(k, 843) + lmat(k, 843)
         mat(k, 850) = mat(k, 850) + lmat(k, 850)
         mat(k, 851) = mat(k, 851) + lmat(k, 851)
         mat(k, 856) = mat(k, 856) + lmat(k, 856)
         mat(k, 863) = mat(k, 863) + lmat(k, 863)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 873) = lmat(k, 873)
         mat(k, 875) = mat(k, 875) + lmat(k, 875)
         mat(k, 876) = lmat(k, 876)
         mat(k, 880) = mat(k, 880) + lmat(k, 880)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 923) = mat(k, 923) + lmat(k, 923)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 945) = mat(k, 945) + lmat(k, 945)
         mat(k, 946) = lmat(k, 946)
         mat(k, 949) = lmat(k, 949)
         mat(k, 952) = mat(k, 952) + lmat(k, 952)
         mat(k, 953) = mat(k, 953) + lmat(k, 953)
         mat(k, 955) = mat(k, 955) + lmat(k, 955)
         mat(k, 957) = lmat(k, 957)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 960) = mat(k, 960) + lmat(k, 960)
         mat(k, 964) = lmat(k, 964)
         mat(k, 968) = lmat(k, 968)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 988) = mat(k, 988) + lmat(k, 988)
         mat(k,1016) = mat(k,1016) + lmat(k,1016)
         mat(k,1040) = mat(k,1040) + lmat(k,1040)
         mat(k,1051) = lmat(k,1051)
         mat(k,1052) = mat(k,1052) + lmat(k,1052)
         mat(k,1053) = mat(k,1053) + lmat(k,1053)
         mat(k,1056) = mat(k,1056) + lmat(k,1056)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1070) = mat(k,1070) + lmat(k,1070)
         mat(k,1072) = lmat(k,1072)
         mat(k,1073) = lmat(k,1073)
         mat(k,1077) = lmat(k,1077)
         mat(k,1078) = mat(k,1078) + lmat(k,1078)
         mat(k,1080) = lmat(k,1080)
         mat(k,1082) = lmat(k,1082)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1105) = mat(k,1105) + lmat(k,1105)
         mat(k,1107) = mat(k,1107) + lmat(k,1107)
         mat(k,1109) = mat(k,1109) + lmat(k,1109)
         mat(k,1111) = lmat(k,1111)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1121) = lmat(k,1121)
         mat(k,1122) = lmat(k,1122)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1131) = mat(k,1131) + lmat(k,1131)
         mat(k,1149) = mat(k,1149) + lmat(k,1149)
         mat(k,1150) = lmat(k,1150)
         mat(k,1155) = lmat(k,1155)
         mat(k,1156) = lmat(k,1156)
         mat(k,1158) = mat(k,1158) + lmat(k,1158)
         mat(k,1163) = lmat(k,1163)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1167) = mat(k,1167) + lmat(k,1167)
         mat(k,1168) = mat(k,1168) + lmat(k,1168)
         mat(k,1176) = mat(k,1176) + lmat(k,1176)
         mat(k,1189) = lmat(k,1189)
         mat(k,1190) = lmat(k,1190)
         mat(k,1191) = lmat(k,1191)
         mat(k,1192) = lmat(k,1192)
         mat(k,1193) = mat(k,1193) + lmat(k,1193)
         mat(k,1194) = lmat(k,1194)
         mat(k,1196) = lmat(k,1196)
         mat(k,1199) = lmat(k,1199)
         mat(k,1201) = lmat(k,1201)
         mat(k,1202) = mat(k,1202) + lmat(k,1202)
         mat(k,1204) = lmat(k,1204)
         mat(k,1206) = mat(k,1206) + lmat(k,1206)
         mat(k,1208) = lmat(k,1208)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1211) = lmat(k,1211)
         mat(k,1215) = mat(k,1215) + lmat(k,1215)
         mat(k,1231) = mat(k,1231) + lmat(k,1231)
         mat(k,1251) = mat(k,1251) + lmat(k,1251)
         mat(k,1266) = mat(k,1266) + lmat(k,1266)
         mat(k,1267) = mat(k,1267) + lmat(k,1267)
         mat(k,1270) = mat(k,1270) + lmat(k,1270)
         mat(k,1271) = mat(k,1271) + lmat(k,1271)
         mat(k,1274) = mat(k,1274) + lmat(k,1274)
         mat(k,1275) = mat(k,1275) + lmat(k,1275)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1279) = mat(k,1279) + lmat(k,1279)
         mat(k,1280) = mat(k,1280) + lmat(k,1280)
         mat(k,1285) = lmat(k,1285)
         mat(k,1297) = mat(k,1297) + lmat(k,1297)
         mat(k,1313) = lmat(k,1313)
         mat(k,1330) = mat(k,1330) + lmat(k,1330)
         mat(k,1341) = mat(k,1341) + lmat(k,1341)
         mat(k,1354) = mat(k,1354) + lmat(k,1354)
         mat(k,1368) = lmat(k,1368)
         mat(k,1370) = mat(k,1370) + lmat(k,1370)
         mat(k,1374) = mat(k,1374) + lmat(k,1374)
         mat(k,1376) = mat(k,1376) + lmat(k,1376)
         mat(k,1385) = lmat(k,1385)
         mat(k,1395) = mat(k,1395) + lmat(k,1395)
         mat(k,1426) = mat(k,1426) + lmat(k,1426)
         mat(k,1447) = mat(k,1447) + lmat(k,1447)
         mat(k,1448) = mat(k,1448) + lmat(k,1448)
         mat(k,1456) = lmat(k,1456)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1472) = lmat(k,1472)
         mat(k,1474) = mat(k,1474) + lmat(k,1474)
         mat(k,1485) = mat(k,1485) + lmat(k,1485)
         mat(k,1493) = mat(k,1493) + lmat(k,1493)
         mat(k,1501) = mat(k,1501) + lmat(k,1501)
         mat(k,1503) = mat(k,1503) + lmat(k,1503)
         mat(k,1509) = mat(k,1509) + lmat(k,1509)
         mat(k,1529) = mat(k,1529) + lmat(k,1529)
         mat(k,1531) = mat(k,1531) + lmat(k,1531)
         mat(k,1532) = lmat(k,1532)
         mat(k,1540) = mat(k,1540) + lmat(k,1540)
         mat(k,1543) = mat(k,1543) + lmat(k,1543)
         mat(k,1550) = mat(k,1550) + lmat(k,1550)
         mat(k,1561) = mat(k,1561) + lmat(k,1561)
         mat(k,1563) = mat(k,1563) + lmat(k,1563)
         mat(k,1574) = mat(k,1574) + lmat(k,1574)
         mat(k,1590) = mat(k,1590) + lmat(k,1590)
         mat(k,1595) = mat(k,1595) + lmat(k,1595)
         mat(k,1601) = mat(k,1601) + lmat(k,1601)
         mat(k,1633) = mat(k,1633) + lmat(k,1633)
         mat(k,1644) = mat(k,1644) + lmat(k,1644)
         mat(k,1799) = mat(k,1799) + lmat(k,1799)
         mat(k,1844) = mat(k,1844) + lmat(k,1844)
         mat(k,1845) = mat(k,1845) + lmat(k,1845)
         mat(k,1848) = mat(k,1848) + lmat(k,1848)
         mat(k,1849) = mat(k,1849) + lmat(k,1849)
         mat(k,1854) = mat(k,1854) + lmat(k,1854)
         mat(k,1900) = mat(k,1900) + lmat(k,1900)
         mat(k,1905) = mat(k,1905) + lmat(k,1905)
         mat(k,1906) = mat(k,1906) + lmat(k,1906)
         mat(k,1908) = mat(k,1908) + lmat(k,1908)
         mat(k,1909) = mat(k,1909) + lmat(k,1909)
         mat(k,1914) = mat(k,1914) + lmat(k,1914)
         mat(k,1953) = mat(k,1953) + lmat(k,1953)
         mat(k,1973) = mat(k,1973) + lmat(k,1973)
         mat(k,1974) = lmat(k,1974)
         mat(k,1977) = mat(k,1977) + lmat(k,1977)
         mat(k,2028) = mat(k,2028) + lmat(k,2028)
         mat(k,2030) = lmat(k,2030)
         mat(k,2036) = mat(k,2036) + lmat(k,2036)
         mat(k,2073) = mat(k,2073) + lmat(k,2073)
         mat(k,2078) = mat(k,2078) + lmat(k,2078)
         mat(k,2095) = mat(k,2095) + lmat(k,2095)
         mat(k,2204) = mat(k,2204) + lmat(k,2204)
         mat(k,2210) = mat(k,2210) + lmat(k,2210)
         mat(k,2257) = mat(k,2257) + lmat(k,2257)
         mat(k,2265) = lmat(k,2265)
         mat(k,2266) = lmat(k,2266)
         mat(k,2267) = mat(k,2267) + lmat(k,2267)
         mat(k,2274) = mat(k,2274) + lmat(k,2274)
         mat(k,2280) = mat(k,2280) + lmat(k,2280)
         mat(k,2282) = mat(k,2282) + lmat(k,2282)
         mat(k,2283) = mat(k,2283) + lmat(k,2283)
         mat(k,2284) = lmat(k,2284)
         mat(k,2285) = mat(k,2285) + lmat(k,2285)
         mat(k,2287) = mat(k,2287) + lmat(k,2287)
         mat(k,2298) = mat(k,2298) + lmat(k,2298)
         mat(k,2303) = lmat(k,2303)
         mat(k,2327) = mat(k,2327) + lmat(k,2327)
         mat(k,2334) = mat(k,2334) + lmat(k,2334)
         mat(k,2336) = lmat(k,2336)
         mat(k,2349) = mat(k,2349) + lmat(k,2349)
         mat(k,2354) = mat(k,2354) + lmat(k,2354)
         mat(k,2362) = mat(k,2362) + lmat(k,2362)
         mat(k,2403) = mat(k,2403) + lmat(k,2403)
         mat(k,2406) = mat(k,2406) + lmat(k,2406)
         mat(k,2417) = mat(k,2417) + lmat(k,2417)
         mat(k,2419) = mat(k,2419) + lmat(k,2419)
         mat(k,2426) = lmat(k,2426)
         mat(k,2433) = mat(k,2433) + lmat(k,2433)
         mat(k,2434) = mat(k,2434) + lmat(k,2434)
         mat(k,2440) = lmat(k,2440)
         mat(k,2444) = lmat(k,2444)
         mat(k,2447) = mat(k,2447) + lmat(k,2447)
         mat(k, 210) = 0._r8
         mat(k, 211) = 0._r8
         mat(k, 250) = 0._r8
         mat(k, 273) = 0._r8
         mat(k, 345) = 0._r8
         mat(k, 434) = 0._r8
         mat(k, 435) = 0._r8
         mat(k, 454) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 489) = 0._r8
         mat(k, 511) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 645) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 681) = 0._r8
         mat(k, 682) = 0._r8
         mat(k, 686) = 0._r8
         mat(k, 687) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 720) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 723) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 764) = 0._r8
         mat(k, 766) = 0._r8
         mat(k, 767) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 790) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 820) = 0._r8
         mat(k, 825) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 847) = 0._r8
         mat(k, 848) = 0._r8
         mat(k, 860) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 902) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 943) = 0._r8
         mat(k, 950) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 971) = 0._r8
         mat(k, 973) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 990) = 0._r8
         mat(k, 998) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1045) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1092) = 0._r8
         mat(k,1095) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1110) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1118) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1319) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1327) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1332) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1371) = 0._r8
         mat(k,1375) = 0._r8
         mat(k,1378) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1382) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1394) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1407) = 0._r8
         mat(k,1428) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1432) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1508) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1518) = 0._r8
         mat(k,1521) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1545) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1551) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1560) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1568) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1576) = 0._r8
         mat(k,1577) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1594) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1636) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1728) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1831) = 0._r8
         mat(k,1832) = 0._r8
         mat(k,1833) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1843) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1852) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1855) = 0._r8
         mat(k,1857) = 0._r8
         mat(k,1864) = 0._r8
         mat(k,1867) = 0._r8
         mat(k,1871) = 0._r8
         mat(k,1873) = 0._r8
         mat(k,1877) = 0._r8
         mat(k,1884) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1897) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1899) = 0._r8
         mat(k,1901) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1903) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1910) = 0._r8
         mat(k,1912) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1931) = 0._r8
         mat(k,1934) = 0._r8
         mat(k,1935) = 0._r8
         mat(k,1936) = 0._r8
         mat(k,1938) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1940) = 0._r8
         mat(k,1941) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1949) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1963) = 0._r8
         mat(k,1966) = 0._r8
         mat(k,1967) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1979) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2081) = 0._r8
         mat(k,2083) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2090) = 0._r8
         mat(k,2091) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2093) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2100) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2161) = 0._r8
         mat(k,2169) = 0._r8
         mat(k,2170) = 0._r8
         mat(k,2172) = 0._r8
         mat(k,2175) = 0._r8
         mat(k,2177) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2196) = 0._r8
         mat(k,2201) = 0._r8
         mat(k,2206) = 0._r8
         mat(k,2218) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2248) = 0._r8
         mat(k,2249) = 0._r8
         mat(k,2251) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2255) = 0._r8
         mat(k,2258) = 0._r8
         mat(k,2259) = 0._r8
         mat(k,2261) = 0._r8
         mat(k,2262) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2270) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2272) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2296) = 0._r8
         mat(k,2302) = 0._r8
         mat(k,2304) = 0._r8
         mat(k,2308) = 0._r8
         mat(k,2316) = 0._r8
         mat(k,2321) = 0._r8
         mat(k,2325) = 0._r8
         mat(k,2326) = 0._r8
         mat(k,2330) = 0._r8
         mat(k,2333) = 0._r8
         mat(k,2335) = 0._r8
         mat(k,2339) = 0._r8
         mat(k,2340) = 0._r8
         mat(k,2341) = 0._r8
         mat(k,2342) = 0._r8
         mat(k,2344) = 0._r8
         mat(k,2348) = 0._r8
         mat(k,2351) = 0._r8
         mat(k,2352) = 0._r8
         mat(k,2355) = 0._r8
         mat(k,2371) = 0._r8
         mat(k,2377) = 0._r8
         mat(k,2378) = 0._r8
         mat(k,2379) = 0._r8
         mat(k,2382) = 0._r8
         mat(k,2387) = 0._r8
         mat(k,2388) = 0._r8
         mat(k,2389) = 0._r8
         mat(k,2391) = 0._r8
         mat(k,2394) = 0._r8
         mat(k,2395) = 0._r8
         mat(k,2396) = 0._r8
         mat(k,2398) = 0._r8
         mat(k,2411) = 0._r8
         mat(k,2420) = 0._r8
         mat(k,2425) = 0._r8
         mat(k,2427) = 0._r8
         mat(k,2428) = 0._r8
         mat(k,2429) = 0._r8
         mat(k,2430) = 0._r8
         mat(k,2431) = 0._r8
         mat(k,2432) = 0._r8
         mat(k,2435) = 0._r8
         mat(k,2436) = 0._r8
         mat(k,2437) = 0._r8
         mat(k,2438) = 0._r8
         mat(k,2439) = 0._r8
         mat(k,2441) = 0._r8
         mat(k,2442) = 0._r8
         mat(k,2443) = 0._r8
         mat(k,2445) = 0._r8
         mat(k,2446) = 0._r8
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
         mat(k, 34) = mat(k, 34) - dti(k)
         mat(k, 35) = mat(k, 35) - dti(k)
         mat(k, 36) = mat(k, 36) - dti(k)
         mat(k, 37) = mat(k, 37) - dti(k)
         mat(k, 38) = mat(k, 38) - dti(k)
         mat(k, 39) = mat(k, 39) - dti(k)
         mat(k, 40) = mat(k, 40) - dti(k)
         mat(k, 41) = mat(k, 41) - dti(k)
         mat(k, 42) = mat(k, 42) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 54) = mat(k, 54) - dti(k)
         mat(k, 60) = mat(k, 60) - dti(k)
         mat(k, 66) = mat(k, 66) - dti(k)
         mat(k, 72) = mat(k, 72) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 80) = mat(k, 80) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 109) = mat(k, 109) - dti(k)
         mat(k, 113) = mat(k, 113) - dti(k)
         mat(k, 117) = mat(k, 117) - dti(k)
         mat(k, 121) = mat(k, 121) - dti(k)
         mat(k, 125) = mat(k, 125) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 137) = mat(k, 137) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 192) = mat(k, 192) - dti(k)
         mat(k, 197) = mat(k, 197) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 206) = mat(k, 206) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 236) = mat(k, 236) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 245) = mat(k, 245) - dti(k)
         mat(k, 249) = mat(k, 249) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 263) = mat(k, 263) - dti(k)
         mat(k, 267) = mat(k, 267) - dti(k)
         mat(k, 272) = mat(k, 272) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 296) = mat(k, 296) - dti(k)
         mat(k, 302) = mat(k, 302) - dti(k)
         mat(k, 308) = mat(k, 308) - dti(k)
         mat(k, 314) = mat(k, 314) - dti(k)
         mat(k, 319) = mat(k, 319) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 332) = mat(k, 332) - dti(k)
         mat(k, 337) = mat(k, 337) - dti(k)
         mat(k, 343) = mat(k, 343) - dti(k)
         mat(k, 348) = mat(k, 348) - dti(k)
         mat(k, 353) = mat(k, 353) - dti(k)
         mat(k, 361) = mat(k, 361) - dti(k)
         mat(k, 369) = mat(k, 369) - dti(k)
         mat(k, 377) = mat(k, 377) - dti(k)
         mat(k, 383) = mat(k, 383) - dti(k)
         mat(k, 389) = mat(k, 389) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 407) = mat(k, 407) - dti(k)
         mat(k, 413) = mat(k, 413) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 425) = mat(k, 425) - dti(k)
         mat(k, 433) = mat(k, 433) - dti(k)
         mat(k, 439) = mat(k, 439) - dti(k)
         mat(k, 445) = mat(k, 445) - dti(k)
         mat(k, 452) = mat(k, 452) - dti(k)
         mat(k, 458) = mat(k, 458) - dti(k)
         mat(k, 461) = mat(k, 461) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 470) = mat(k, 470) - dti(k)
         mat(k, 477) = mat(k, 477) - dti(k)
         mat(k, 486) = mat(k, 486) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 500) = mat(k, 500) - dti(k)
         mat(k, 508) = mat(k, 508) - dti(k)
         mat(k, 515) = mat(k, 515) - dti(k)
         mat(k, 521) = mat(k, 521) - dti(k)
         mat(k, 527) = mat(k, 527) - dti(k)
         mat(k, 532) = mat(k, 532) - dti(k)
         mat(k, 540) = mat(k, 540) - dti(k)
         mat(k, 548) = mat(k, 548) - dti(k)
         mat(k, 556) = mat(k, 556) - dti(k)
         mat(k, 564) = mat(k, 564) - dti(k)
         mat(k, 572) = mat(k, 572) - dti(k)
         mat(k, 580) = mat(k, 580) - dti(k)
         mat(k, 589) = mat(k, 589) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 602) = mat(k, 602) - dti(k)
         mat(k, 611) = mat(k, 611) - dti(k)
         mat(k, 618) = mat(k, 618) - dti(k)
         mat(k, 625) = mat(k, 625) - dti(k)
         mat(k, 633) = mat(k, 633) - dti(k)
         mat(k, 640) = mat(k, 640) - dti(k)
         mat(k, 650) = mat(k, 650) - dti(k)
         mat(k, 663) = mat(k, 663) - dti(k)
         mat(k, 674) = mat(k, 674) - dti(k)
         mat(k, 685) = mat(k, 685) - dti(k)
         mat(k, 696) = mat(k, 696) - dti(k)
         mat(k, 703) = mat(k, 703) - dti(k)
         mat(k, 708) = mat(k, 708) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 731) = mat(k, 731) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 749) = mat(k, 749) - dti(k)
         mat(k, 765) = mat(k, 765) - dti(k)
         mat(k, 776) = mat(k, 776) - dti(k)
         mat(k, 785) = mat(k, 785) - dti(k)
         mat(k, 795) = mat(k, 795) - dti(k)
         mat(k, 803) = mat(k, 803) - dti(k)
         mat(k, 811) = mat(k, 811) - dti(k)
         mat(k, 816) = mat(k, 816) - dti(k)
         mat(k, 826) = mat(k, 826) - dti(k)
         mat(k, 835) = mat(k, 835) - dti(k)
         mat(k, 843) = mat(k, 843) - dti(k)
         mat(k, 851) = mat(k, 851) - dti(k)
         mat(k, 863) = mat(k, 863) - dti(k)
         mat(k, 871) = mat(k, 871) - dti(k)
         mat(k, 880) = mat(k, 880) - dti(k)
         mat(k, 899) = mat(k, 899) - dti(k)
         mat(k, 908) = mat(k, 908) - dti(k)
         mat(k, 913) = mat(k, 913) - dti(k)
         mat(k, 923) = mat(k, 923) - dti(k)
         mat(k, 933) = mat(k, 933) - dti(k)
         mat(k, 945) = mat(k, 945) - dti(k)
         mat(k, 953) = mat(k, 953) - dti(k)
         mat(k, 969) = mat(k, 969) - dti(k)
         mat(k, 988) = mat(k, 988) - dti(k)
         mat(k,1016) = mat(k,1016) - dti(k)
         mat(k,1040) = mat(k,1040) - dti(k)
         mat(k,1052) = mat(k,1052) - dti(k)
         mat(k,1060) = mat(k,1060) - dti(k)
         mat(k,1070) = mat(k,1070) - dti(k)
         mat(k,1078) = mat(k,1078) - dti(k)
         mat(k,1088) = mat(k,1088) - dti(k)
         mat(k,1101) = mat(k,1101) - dti(k)
         mat(k,1115) = mat(k,1115) - dti(k)
         mat(k,1131) = mat(k,1131) - dti(k)
         mat(k,1149) = mat(k,1149) - dti(k)
         mat(k,1158) = mat(k,1158) - dti(k)
         mat(k,1164) = mat(k,1164) - dti(k)
         mat(k,1176) = mat(k,1176) - dti(k)
         mat(k,1193) = mat(k,1193) - dti(k)
         mat(k,1206) = mat(k,1206) - dti(k)
         mat(k,1215) = mat(k,1215) - dti(k)
         mat(k,1231) = mat(k,1231) - dti(k)
         mat(k,1251) = mat(k,1251) - dti(k)
         mat(k,1267) = mat(k,1267) - dti(k)
         mat(k,1279) = mat(k,1279) - dti(k)
         mat(k,1297) = mat(k,1297) - dti(k)
         mat(k,1330) = mat(k,1330) - dti(k)
         mat(k,1354) = mat(k,1354) - dti(k)
         mat(k,1374) = mat(k,1374) - dti(k)
         mat(k,1395) = mat(k,1395) - dti(k)
         mat(k,1426) = mat(k,1426) - dti(k)
         mat(k,1448) = mat(k,1448) - dti(k)
         mat(k,1459) = mat(k,1459) - dti(k)
         mat(k,1474) = mat(k,1474) - dti(k)
         mat(k,1493) = mat(k,1493) - dti(k)
         mat(k,1509) = mat(k,1509) - dti(k)
         mat(k,1540) = mat(k,1540) - dti(k)
         mat(k,1563) = mat(k,1563) - dti(k)
         mat(k,1590) = mat(k,1590) - dti(k)
         mat(k,1633) = mat(k,1633) - dti(k)
         mat(k,1799) = mat(k,1799) - dti(k)
         mat(k,1845) = mat(k,1845) - dti(k)
         mat(k,1906) = mat(k,1906) - dti(k)
         mat(k,1953) = mat(k,1953) - dti(k)
         mat(k,1977) = mat(k,1977) - dti(k)
         mat(k,2073) = mat(k,2073) - dti(k)
         mat(k,2095) = mat(k,2095) - dti(k)
         mat(k,2204) = mat(k,2204) - dti(k)
         mat(k,2257) = mat(k,2257) - dti(k)
         mat(k,2283) = mat(k,2283) - dti(k)
         mat(k,2327) = mat(k,2327) - dti(k)
         mat(k,2354) = mat(k,2354) - dti(k)
         mat(k,2419) = mat(k,2419) - dti(k)
         mat(k,2447) = mat(k,2447) - dti(k)
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
      call nlnmat11( avec_len, mat, y, rxt )
      call nlnmat_finit( avec_len, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
