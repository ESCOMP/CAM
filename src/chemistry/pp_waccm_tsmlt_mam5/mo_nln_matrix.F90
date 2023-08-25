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
         mat(k,642) = -(rxt(k,396)*y(k,228))
         mat(k,1753) = -rxt(k,396)*y(k,1)
         mat(k,1865) = rxt(k,399)*y(k,192)
         mat(k,1038) = rxt(k,399)*y(k,124)
         mat(k,687) = -(rxt(k,400)*y(k,228))
         mat(k,1756) = -rxt(k,400)*y(k,2)
         mat(k,2307) = rxt(k,397)*y(k,192)
         mat(k,1039) = rxt(k,397)*y(k,90)
         mat(k,974) = -(rxt(k,479)*y(k,126) + rxt(k,480)*y(k,136) + rxt(k,481) &
                      *y(k,228))
         mat(k,1618) = -rxt(k,479)*y(k,6)
         mat(k,2192) = -rxt(k,480)*y(k,6)
         mat(k,1780) = -rxt(k,481)*y(k,6)
         mat(k,160) = -(rxt(k,438)*y(k,228))
         mat(k,1685) = -rxt(k,438)*y(k,7)
         mat(k,415) = -(rxt(k,441)*y(k,228))
         mat(k,1725) = -rxt(k,441)*y(k,8)
         mat(k,2287) = rxt(k,439)*y(k,194)
         mat(k,493) = rxt(k,439)*y(k,90)
         mat(k,161) = .120_r8*rxt(k,438)*y(k,228)
         mat(k,1686) = .120_r8*rxt(k,438)*y(k,7)
         mat(k,972) = .100_r8*rxt(k,480)*y(k,136)
         mat(k,1016) = .100_r8*rxt(k,483)*y(k,136)
         mat(k,2182) = .100_r8*rxt(k,480)*y(k,6) + .100_r8*rxt(k,483)*y(k,110)
         mat(k,1853) = .500_r8*rxt(k,440)*y(k,194) + .200_r8*rxt(k,467)*y(k,235) &
                      + .060_r8*rxt(k,473)*y(k,238)
         mat(k,494) = .500_r8*rxt(k,440)*y(k,124)
         mat(k,747) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,763) = .060_r8*rxt(k,473)*y(k,124)
         mat(k,1846) = .200_r8*rxt(k,467)*y(k,235) + .200_r8*rxt(k,473)*y(k,238)
         mat(k,746) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,761) = .200_r8*rxt(k,473)*y(k,124)
         mat(k,1862) = .200_r8*rxt(k,467)*y(k,235) + .150_r8*rxt(k,473)*y(k,238)
         mat(k,748) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,764) = .150_r8*rxt(k,473)*y(k,124)
         mat(k,1848) = .210_r8*rxt(k,473)*y(k,238)
         mat(k,762) = .210_r8*rxt(k,473)*y(k,124)
         mat(k,226) = -(rxt(k,401)*y(k,228))
         mat(k,1696) = -rxt(k,401)*y(k,15)
         mat(k,971) = .050_r8*rxt(k,480)*y(k,136)
         mat(k,1015) = .050_r8*rxt(k,483)*y(k,136)
         mat(k,2181) = .050_r8*rxt(k,480)*y(k,6) + .050_r8*rxt(k,483)*y(k,110)
         mat(k,355) = -(rxt(k,367)*y(k,126) + rxt(k,368)*y(k,228))
         mat(k,1612) = -rxt(k,367)*y(k,16)
         mat(k,1717) = -rxt(k,368)*y(k,16)
         mat(k,1511) = -(rxt(k,250)*y(k,42) + rxt(k,251)*y(k,90) + rxt(k,252)*y(k,136))
         mat(k,1979) = -rxt(k,250)*y(k,17)
         mat(k,2352) = -rxt(k,251)*y(k,17)
         mat(k,2219) = -rxt(k,252)*y(k,17)
         mat(k,1563) = 4.000_r8*rxt(k,253)*y(k,19) + (rxt(k,254)+rxt(k,255))*y(k,59) &
                      + rxt(k,258)*y(k,124) + rxt(k,261)*y(k,134) + rxt(k,509) &
                      *y(k,152) + rxt(k,262)*y(k,228)
         mat(k,141) = rxt(k,240)*y(k,224)
         mat(k,147) = rxt(k,266)*y(k,224)
         mat(k,481) = 2.000_r8*rxt(k,277)*y(k,56) + 2.000_r8*rxt(k,289)*y(k,224) &
                      + 2.000_r8*rxt(k,278)*y(k,228)
         mat(k,604) = rxt(k,279)*y(k,56) + rxt(k,290)*y(k,224) + rxt(k,280)*y(k,228)
         mat(k,387) = 3.000_r8*rxt(k,284)*y(k,56) + 3.000_r8*rxt(k,267)*y(k,224) &
                      + 3.000_r8*rxt(k,285)*y(k,228)
         mat(k,2155) = 2.000_r8*rxt(k,277)*y(k,41) + rxt(k,279)*y(k,43) &
                      + 3.000_r8*rxt(k,284)*y(k,55)
         mat(k,1589) = (rxt(k,254)+rxt(k,255))*y(k,19)
         mat(k,109) = 2.000_r8*rxt(k,268)*y(k,224)
         mat(k,829) = rxt(k,263)*y(k,134) + rxt(k,269)*y(k,224) + rxt(k,264)*y(k,228)
         mat(k,1908) = rxt(k,258)*y(k,19)
         mat(k,2088) = rxt(k,261)*y(k,19) + rxt(k,263)*y(k,81)
         mat(k,1477) = rxt(k,509)*y(k,19)
         mat(k,2022) = rxt(k,240)*y(k,34) + rxt(k,266)*y(k,35) + 2.000_r8*rxt(k,289) &
                      *y(k,41) + rxt(k,290)*y(k,43) + 3.000_r8*rxt(k,267)*y(k,55) &
                      + 2.000_r8*rxt(k,268)*y(k,78) + rxt(k,269)*y(k,81)
         mat(k,1813) = rxt(k,262)*y(k,19) + 2.000_r8*rxt(k,278)*y(k,41) + rxt(k,280) &
                      *y(k,43) + 3.000_r8*rxt(k,285)*y(k,55) + rxt(k,264)*y(k,81)
         mat(k,1556) = rxt(k,256)*y(k,59)
         mat(k,1582) = rxt(k,256)*y(k,19)
         mat(k,1491) = (rxt(k,570)+rxt(k,575))*y(k,92)
         mat(k,786) = (rxt(k,570)+rxt(k,575))*y(k,85)
         mat(k,1565) = -(4._r8*rxt(k,253)*y(k,19) + (rxt(k,254) + rxt(k,255) + rxt(k,256) &
                      ) * y(k,59) + rxt(k,257)*y(k,90) + rxt(k,258)*y(k,124) + rxt(k,259) &
                      *y(k,125) + rxt(k,261)*y(k,134) + rxt(k,262)*y(k,228) + rxt(k,509) &
                      *y(k,152))
         mat(k,1591) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,19)
         mat(k,2354) = -rxt(k,257)*y(k,19)
         mat(k,1910) = -rxt(k,258)*y(k,19)
         mat(k,1955) = -rxt(k,259)*y(k,19)
         mat(k,2090) = -rxt(k,261)*y(k,19)
         mat(k,1815) = -rxt(k,262)*y(k,19)
         mat(k,1479) = -rxt(k,509)*y(k,19)
         mat(k,1513) = rxt(k,252)*y(k,136)
         mat(k,569) = rxt(k,260)*y(k,134)
         mat(k,830) = rxt(k,270)*y(k,224)
         mat(k,790) = rxt(k,265)*y(k,134)
         mat(k,2090) = mat(k,2090) + rxt(k,260)*y(k,20) + rxt(k,265)*y(k,92)
         mat(k,2221) = rxt(k,252)*y(k,17)
         mat(k,2024) = rxt(k,270)*y(k,81)
         mat(k,566) = -(rxt(k,260)*y(k,134))
         mat(k,2069) = -rxt(k,260)*y(k,20)
         mat(k,1558) = rxt(k,259)*y(k,125)
         mat(k,1934) = rxt(k,259)*y(k,19)
         mat(k,235) = -(rxt(k,442)*y(k,228))
         mat(k,1698) = -rxt(k,442)*y(k,22)
         mat(k,1845) = rxt(k,445)*y(k,196)
         mat(k,433) = rxt(k,445)*y(k,124)
         mat(k,345) = -(rxt(k,444)*y(k,228))
         mat(k,1715) = -rxt(k,444)*y(k,23)
         mat(k,2281) = rxt(k,443)*y(k,196)
         mat(k,434) = rxt(k,443)*y(k,90)
         mat(k,285) = -(rxt(k,315)*y(k,56) + rxt(k,316)*y(k,228))
         mat(k,2129) = -rxt(k,315)*y(k,24)
         mat(k,1706) = -rxt(k,316)*y(k,24)
         mat(k,550) = -(rxt(k,317)*y(k,56) + rxt(k,318)*y(k,136) + rxt(k,343)*y(k,228))
         mat(k,2135) = -rxt(k,317)*y(k,25)
         mat(k,2184) = -rxt(k,318)*y(k,25)
         mat(k,1742) = -rxt(k,343)*y(k,25)
         mat(k,265) = -(rxt(k,323)*y(k,228))
         mat(k,1704) = -rxt(k,323)*y(k,26)
         mat(k,898) = .800_r8*rxt(k,319)*y(k,197) + .200_r8*rxt(k,320)*y(k,201)
         mat(k,2371) = .200_r8*rxt(k,320)*y(k,197)
         mat(k,350) = -(rxt(k,324)*y(k,228))
         mat(k,1716) = -rxt(k,324)*y(k,27)
         mat(k,2282) = rxt(k,321)*y(k,197)
         mat(k,899) = rxt(k,321)*y(k,90)
         mat(k,298) = -(rxt(k,325)*y(k,56) + rxt(k,326)*y(k,228))
         mat(k,2130) = -rxt(k,325)*y(k,28)
         mat(k,1708) = -rxt(k,326)*y(k,28)
         mat(k,1133) = -(rxt(k,346)*y(k,126) + rxt(k,347)*y(k,136) + rxt(k,365) &
                      *y(k,228))
         mat(k,1628) = -rxt(k,346)*y(k,29)
         mat(k,2201) = -rxt(k,347)*y(k,29)
         mat(k,1791) = -rxt(k,365)*y(k,29)
         mat(k,878) = .130_r8*rxt(k,425)*y(k,136)
         mat(k,2201) = mat(k,2201) + .130_r8*rxt(k,425)*y(k,99)
         mat(k,409) = -(rxt(k,351)*y(k,228))
         mat(k,1724) = -rxt(k,351)*y(k,30)
         mat(k,2286) = rxt(k,349)*y(k,198)
         mat(k,934) = rxt(k,349)*y(k,90)
         mat(k,304) = -(rxt(k,352)*y(k,228) + rxt(k,355)*y(k,56))
         mat(k,1709) = -rxt(k,352)*y(k,31)
         mat(k,2131) = -rxt(k,355)*y(k,31)
         mat(k,269) = -(rxt(k,448)*y(k,228))
         mat(k,1705) = -rxt(k,448)*y(k,32)
         mat(k,2277) = rxt(k,446)*y(k,199)
         mat(k,633) = rxt(k,446)*y(k,90)
         mat(k,101) = -(rxt(k,239)*y(k,224))
         mat(k,1998) = -rxt(k,239)*y(k,33)
         mat(k,139) = -(rxt(k,240)*y(k,224))
         mat(k,2003) = -rxt(k,240)*y(k,34)
         mat(k,144) = -(rxt(k,266)*y(k,224))
         mat(k,2004) = -rxt(k,266)*y(k,35)
         mat(k,111) = -(rxt(k,241)*y(k,224))
         mat(k,2000) = -rxt(k,241)*y(k,36)
         mat(k,149) = -(rxt(k,242)*y(k,224))
         mat(k,2005) = -rxt(k,242)*y(k,37)
         mat(k,115) = -(rxt(k,243)*y(k,224))
         mat(k,2001) = -rxt(k,243)*y(k,38)
         mat(k,154) = -(rxt(k,244)*y(k,224))
         mat(k,2006) = -rxt(k,244)*y(k,39)
         mat(k,119) = -(rxt(k,245)*y(k,224))
         mat(k,2002) = -rxt(k,245)*y(k,40)
         mat(k,479) = -(rxt(k,277)*y(k,56) + rxt(k,278)*y(k,228) + rxt(k,289)*y(k,224))
         mat(k,2134) = -rxt(k,277)*y(k,41)
         mat(k,1734) = -rxt(k,278)*y(k,41)
         mat(k,2016) = -rxt(k,289)*y(k,41)
         mat(k,1987) = -(rxt(k,214)*y(k,56) + rxt(k,250)*y(k,17) + rxt(k,294)*y(k,90) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,134) + rxt(k,297) &
                      *y(k,228))
         mat(k,2163) = -rxt(k,214)*y(k,42)
         mat(k,1517) = -rxt(k,250)*y(k,42)
         mat(k,2360) = -rxt(k,294)*y(k,42)
         mat(k,1656) = -rxt(k,295)*y(k,42)
         mat(k,2096) = -rxt(k,296)*y(k,42)
         mat(k,1821) = -rxt(k,297)*y(k,42)
         mat(k,650) = .400_r8*rxt(k,396)*y(k,228)
         mat(k,989) = .340_r8*rxt(k,480)*y(k,136)
         mat(k,362) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,555) = rxt(k,318)*y(k,136)
         mat(k,1145) = .500_r8*rxt(k,347)*y(k,136)
         mat(k,624) = .500_r8*rxt(k,335)*y(k,228)
         mat(k,826) = rxt(k,302)*y(k,228)
         mat(k,456) = .300_r8*rxt(k,303)*y(k,228)
         mat(k,2252) = (rxt(k,311)+rxt(k,312))*y(k,224)
         mat(k,1597) = rxt(k,221)*y(k,201)
         mat(k,1169) = .800_r8*rxt(k,340)*y(k,228)
         mat(k,2360) = mat(k,2360) + .450_r8*rxt(k,383)*y(k,215) + .150_r8*rxt(k,362) &
                      *y(k,232)
         mat(k,888) = .910_r8*rxt(k,425)*y(k,136)
         mat(k,595) = .300_r8*rxt(k,416)*y(k,228)
         mat(k,1276) = .120_r8*rxt(k,378)*y(k,136)
         mat(k,618) = .500_r8*rxt(k,391)*y(k,228)
         mat(k,1033) = .340_r8*rxt(k,483)*y(k,136)
         mat(k,1385) = .600_r8*rxt(k,392)*y(k,136)
         mat(k,1916) = .100_r8*rxt(k,398)*y(k,192) + rxt(k,301)*y(k,201) &
                      + .500_r8*rxt(k,369)*y(k,204) + .500_r8*rxt(k,337)*y(k,206) &
                      + .920_r8*rxt(k,408)*y(k,208) + .250_r8*rxt(k,376)*y(k,213) &
                      + rxt(k,385)*y(k,215) + rxt(k,359)*y(k,231) + rxt(k,363) &
                      *y(k,232) + .340_r8*rxt(k,492)*y(k,233) + .320_r8*rxt(k,497) &
                      *y(k,234) + .250_r8*rxt(k,433)*y(k,237)
         mat(k,1656) = mat(k,1656) + .500_r8*rxt(k,367)*y(k,16) + rxt(k,409)*y(k,208) &
                      + .250_r8*rxt(k,375)*y(k,213) + rxt(k,386)*y(k,215)
         mat(k,2227) = .340_r8*rxt(k,480)*y(k,6) + rxt(k,318)*y(k,25) &
                      + .500_r8*rxt(k,347)*y(k,29) + .910_r8*rxt(k,425)*y(k,99) &
                      + .120_r8*rxt(k,378)*y(k,105) + .340_r8*rxt(k,483)*y(k,110) &
                      + .600_r8*rxt(k,392)*y(k,111)
         mat(k,540) = rxt(k,342)*y(k,228)
         mat(k,1125) = .680_r8*rxt(k,501)*y(k,228)
         mat(k,1050) = .100_r8*rxt(k,398)*y(k,124)
         mat(k,907) = .700_r8*rxt(k,320)*y(k,201)
         mat(k,942) = rxt(k,348)*y(k,201)
         mat(k,1435) = rxt(k,331)*y(k,201) + rxt(k,405)*y(k,208) + .250_r8*rxt(k,372) &
                      *y(k,213) + rxt(k,381)*y(k,215) + .250_r8*rxt(k,430)*y(k,237)
         mat(k,2412) = rxt(k,221)*y(k,59) + rxt(k,301)*y(k,124) + .700_r8*rxt(k,320) &
                      *y(k,197) + rxt(k,348)*y(k,198) + rxt(k,331)*y(k,200) + ( &
                      + 4.000_r8*rxt(k,298)+2.000_r8*rxt(k,299))*y(k,201) &
                      + 1.500_r8*rxt(k,406)*y(k,208) + .750_r8*rxt(k,411)*y(k,209) &
                      + .800_r8*rxt(k,420)*y(k,210) + .880_r8*rxt(k,373)*y(k,213) &
                      + 2.000_r8*rxt(k,382)*y(k,215) + .750_r8*rxt(k,485)*y(k,223) &
                      + .800_r8*rxt(k,361)*y(k,232) + .930_r8*rxt(k,490)*y(k,233) &
                      + .950_r8*rxt(k,495)*y(k,234) + .800_r8*rxt(k,431)*y(k,237)
         mat(k,580) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,803) = .500_r8*rxt(k,337)*y(k,124)
         mat(k,1309) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) + rxt(k,405) &
                      *y(k,200) + 1.500_r8*rxt(k,406)*y(k,201)
         mat(k,1342) = .750_r8*rxt(k,411)*y(k,201)
         mat(k,1263) = .800_r8*rxt(k,420)*y(k,201)
         mat(k,1364) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,200) + .880_r8*rxt(k,373)*y(k,201)
         mat(k,1404) = .450_r8*rxt(k,383)*y(k,90) + rxt(k,385)*y(k,124) + rxt(k,386) &
                      *y(k,126) + rxt(k,381)*y(k,200) + 2.000_r8*rxt(k,382)*y(k,201) &
                      + 4.000_r8*rxt(k,384)*y(k,215)
         mat(k,1114) = .750_r8*rxt(k,485)*y(k,201)
         mat(k,2030) = (rxt(k,311)+rxt(k,312))*y(k,54)
         mat(k,1821) = mat(k,1821) + .400_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,335) &
                      *y(k,51) + rxt(k,302)*y(k,52) + .300_r8*rxt(k,303)*y(k,53) &
                      + .800_r8*rxt(k,340)*y(k,74) + .300_r8*rxt(k,416)*y(k,100) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,342)*y(k,141) &
                      + .680_r8*rxt(k,501)*y(k,181)
         mat(k,821) = rxt(k,359)*y(k,124)
         mat(k,1223) = .150_r8*rxt(k,362)*y(k,90) + rxt(k,363)*y(k,124) &
                      + .800_r8*rxt(k,361)*y(k,201)
         mat(k,1185) = .340_r8*rxt(k,492)*y(k,124) + .930_r8*rxt(k,490)*y(k,201)
         mat(k,1068) = .320_r8*rxt(k,497)*y(k,124) + .950_r8*rxt(k,495)*y(k,201)
         mat(k,1241) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,430)*y(k,200) &
                      + .800_r8*rxt(k,431)*y(k,201)
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
         mat(k,602) = -(rxt(k,279)*y(k,56) + rxt(k,280)*y(k,228) + rxt(k,290)*y(k,224))
         mat(k,2137) = -rxt(k,279)*y(k,43)
         mat(k,1748) = -rxt(k,280)*y(k,43)
         mat(k,2017) = -rxt(k,290)*y(k,43)
         mat(k,123) = -(rxt(k,281)*y(k,228))
         mat(k,1683) = -rxt(k,281)*y(k,44)
         mat(k,1151) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,228))
         mat(k,1629) = -rxt(k,327)*y(k,45)
         mat(k,1792) = -rxt(k,328)*y(k,45)
         mat(k,646) = .800_r8*rxt(k,396)*y(k,228)
         mat(k,358) = rxt(k,367)*y(k,126)
         mat(k,266) = rxt(k,323)*y(k,228)
         mat(k,352) = .500_r8*rxt(k,324)*y(k,228)
         mat(k,1134) = .500_r8*rxt(k,347)*y(k,136)
         mat(k,2333) = .200_r8*rxt(k,387)*y(k,217)
         mat(k,1371) = .100_r8*rxt(k,392)*y(k,136)
         mat(k,1890) = .400_r8*rxt(k,398)*y(k,192) + rxt(k,322)*y(k,197) &
                      + .270_r8*rxt(k,350)*y(k,198) + rxt(k,369)*y(k,204) + rxt(k,388) &
                      *y(k,217) + rxt(k,359)*y(k,231)
         mat(k,1629) = mat(k,1629) + rxt(k,367)*y(k,16)
         mat(k,2202) = .500_r8*rxt(k,347)*y(k,29) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,1044) = .400_r8*rxt(k,398)*y(k,124)
         mat(k,902) = rxt(k,322)*y(k,124) + 3.200_r8*rxt(k,319)*y(k,197) &
                      + .800_r8*rxt(k,320)*y(k,201)
         mat(k,937) = .270_r8*rxt(k,350)*y(k,124)
         mat(k,2388) = .800_r8*rxt(k,320)*y(k,197)
         mat(k,576) = rxt(k,369)*y(k,124)
         mat(k,699) = .200_r8*rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124)
         mat(k,1792) = mat(k,1792) + .800_r8*rxt(k,396)*y(k,1) + rxt(k,323)*y(k,26) &
                      + .500_r8*rxt(k,324)*y(k,27)
         mat(k,815) = rxt(k,359)*y(k,124)
         mat(k,371) = -(rxt(k,282)*y(k,56) + rxt(k,283)*y(k,228))
         mat(k,2132) = -rxt(k,282)*y(k,46)
         mat(k,1719) = -rxt(k,283)*y(k,46)
         mat(k,104) = -(rxt(k,329)*y(k,228))
         mat(k,1682) = -rxt(k,329)*y(k,47)
         mat(k,1080) = -(rxt(k,366)*y(k,228))
         mat(k,1787) = -rxt(k,366)*y(k,48)
         mat(k,645) = .800_r8*rxt(k,396)*y(k,228)
         mat(k,979) = .520_r8*rxt(k,480)*y(k,136)
         mat(k,357) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,1023) = .520_r8*rxt(k,483)*y(k,136)
         mat(k,1886) = .250_r8*rxt(k,398)*y(k,192) + .820_r8*rxt(k,350)*y(k,198) &
                      + .500_r8*rxt(k,369)*y(k,204) + .270_r8*rxt(k,492)*y(k,233) &
                      + .040_r8*rxt(k,497)*y(k,234)
         mat(k,1624) = .500_r8*rxt(k,367)*y(k,16)
         mat(k,2198) = .520_r8*rxt(k,480)*y(k,6) + .520_r8*rxt(k,483)*y(k,110)
         mat(k,1118) = .500_r8*rxt(k,501)*y(k,228)
         mat(k,1043) = .250_r8*rxt(k,398)*y(k,124)
         mat(k,936) = .820_r8*rxt(k,350)*y(k,124) + .820_r8*rxt(k,348)*y(k,201)
         mat(k,2384) = .820_r8*rxt(k,348)*y(k,198) + .150_r8*rxt(k,490)*y(k,233) &
                      + .025_r8*rxt(k,495)*y(k,234)
         mat(k,575) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,1787) = mat(k,1787) + .800_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,501) &
                      *y(k,181)
         mat(k,1174) = .270_r8*rxt(k,492)*y(k,124) + .150_r8*rxt(k,490)*y(k,201)
         mat(k,1064) = .040_r8*rxt(k,497)*y(k,124) + .025_r8*rxt(k,495)*y(k,201)
         mat(k,1281) = -(rxt(k,353)*y(k,126) + rxt(k,354)*y(k,228))
         mat(k,1639) = -rxt(k,353)*y(k,49)
         mat(k,1802) = -rxt(k,354)*y(k,49)
         mat(k,2342) = .070_r8*rxt(k,450)*y(k,202) + .070_r8*rxt(k,456)*y(k,216)
         mat(k,1209) = rxt(k,356)*y(k,228)
         mat(k,1270) = .880_r8*rxt(k,378)*y(k,136)
         mat(k,1374) = .500_r8*rxt(k,392)*y(k,136)
         mat(k,1900) = .170_r8*rxt(k,451)*y(k,202) + .050_r8*rxt(k,414)*y(k,209) &
                      + .250_r8*rxt(k,376)*y(k,213) + .170_r8*rxt(k,457)*y(k,216) &
                      + .400_r8*rxt(k,467)*y(k,235) + .250_r8*rxt(k,433)*y(k,237) &
                      + .540_r8*rxt(k,473)*y(k,238) + .510_r8*rxt(k,476)*y(k,240)
         mat(k,1639) = mat(k,1639) + .050_r8*rxt(k,415)*y(k,209) + .250_r8*rxt(k,375) &
                      *y(k,213) + .250_r8*rxt(k,434)*y(k,237)
         mat(k,893) = rxt(k,357)*y(k,228)
         mat(k,2210) = .880_r8*rxt(k,378)*y(k,105) + .500_r8*rxt(k,392)*y(k,111)
         mat(k,1422) = .250_r8*rxt(k,372)*y(k,213) + .250_r8*rxt(k,430)*y(k,237)
         mat(k,2397) = .240_r8*rxt(k,373)*y(k,213) + .500_r8*rxt(k,361)*y(k,232) &
                      + .100_r8*rxt(k,431)*y(k,237)
         mat(k,780) = .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451)*y(k,124)
         mat(k,1331) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1355) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,200) + .240_r8*rxt(k,373)*y(k,201)
         mat(k,913) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,1802) = mat(k,1802) + rxt(k,356)*y(k,96) + rxt(k,357)*y(k,127)
         mat(k,1218) = .500_r8*rxt(k,361)*y(k,201)
         mat(k,756) = .400_r8*rxt(k,467)*y(k,124)
         mat(k,1234) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,200) + .100_r8*rxt(k,431)*y(k,201)
         mat(k,772) = .540_r8*rxt(k,473)*y(k,124)
         mat(k,512) = .510_r8*rxt(k,476)*y(k,124)
         mat(k,705) = -(rxt(k,334)*y(k,228))
         mat(k,1758) = -rxt(k,334)*y(k,50)
         mat(k,1129) = .120_r8*rxt(k,347)*y(k,136)
         mat(k,2309) = .150_r8*rxt(k,332)*y(k,200) + .150_r8*rxt(k,383)*y(k,215)
         mat(k,2186) = .120_r8*rxt(k,347)*y(k,29)
         mat(k,1413) = .150_r8*rxt(k,332)*y(k,90) + .100_r8*rxt(k,331)*y(k,201)
         mat(k,2376) = .100_r8*rxt(k,331)*y(k,200)
         mat(k,1394) = .150_r8*rxt(k,383)*y(k,90)
         mat(k,620) = -(rxt(k,335)*y(k,228))
         mat(k,1750) = -rxt(k,335)*y(k,51)
         mat(k,2303) = .400_r8*rxt(k,332)*y(k,200) + .400_r8*rxt(k,383)*y(k,215)
         mat(k,1412) = .400_r8*rxt(k,332)*y(k,90)
         mat(k,1393) = .400_r8*rxt(k,383)*y(k,90)
         mat(k,824) = -(rxt(k,302)*y(k,228))
         mat(k,1768) = -rxt(k,302)*y(k,52)
         mat(k,900) = .300_r8*rxt(k,320)*y(k,201)
         mat(k,2377) = .300_r8*rxt(k,320)*y(k,197) + 2.000_r8*rxt(k,299)*y(k,201) &
                      + .250_r8*rxt(k,406)*y(k,208) + .250_r8*rxt(k,411)*y(k,209) &
                      + .200_r8*rxt(k,420)*y(k,210) + .250_r8*rxt(k,373)*y(k,213) &
                      + .250_r8*rxt(k,485)*y(k,223) + .500_r8*rxt(k,361)*y(k,232) &
                      + .250_r8*rxt(k,490)*y(k,233) + .250_r8*rxt(k,495)*y(k,234) &
                      + .300_r8*rxt(k,431)*y(k,237)
         mat(k,1291) = .250_r8*rxt(k,406)*y(k,201)
         mat(k,1320) = .250_r8*rxt(k,411)*y(k,201)
         mat(k,1247) = .200_r8*rxt(k,420)*y(k,201)
         mat(k,1349) = .250_r8*rxt(k,373)*y(k,201)
         mat(k,1104) = .250_r8*rxt(k,485)*y(k,201)
         mat(k,1215) = .500_r8*rxt(k,361)*y(k,201)
         mat(k,1173) = .250_r8*rxt(k,490)*y(k,201)
         mat(k,1061) = .250_r8*rxt(k,495)*y(k,201)
         mat(k,1228) = .300_r8*rxt(k,431)*y(k,201)
         mat(k,454) = -(rxt(k,303)*y(k,228))
         mat(k,1730) = -rxt(k,303)*y(k,53)
         mat(k,2292) = rxt(k,300)*y(k,201)
         mat(k,2374) = rxt(k,300)*y(k,90)
         mat(k,2259) = -(rxt(k,215)*y(k,56) + rxt(k,271)*y(k,73) + rxt(k,304)*y(k,228) &
                      + (rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,224))
         mat(k,2170) = -rxt(k,215)*y(k,54)
         mat(k,931) = -rxt(k,271)*y(k,54)
         mat(k,1828) = -rxt(k,304)*y(k,54)
         mat(k,2037) = -(rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,54)
         mat(k,1147) = .100_r8*rxt(k,347)*y(k,136)
         mat(k,2234) = .100_r8*rxt(k,347)*y(k,29)
         mat(k,385) = -(rxt(k,267)*y(k,224) + rxt(k,284)*y(k,56) + rxt(k,285)*y(k,228))
         mat(k,2015) = -rxt(k,267)*y(k,55)
         mat(k,2133) = -rxt(k,284)*y(k,55)
         mat(k,1720) = -rxt(k,285)*y(k,55)
         mat(k,2168) = -(rxt(k,214)*y(k,42) + rxt(k,215)*y(k,54) + rxt(k,216)*y(k,77) &
                      + rxt(k,217)*y(k,79) + (rxt(k,218) + rxt(k,219)) * y(k,90) &
                      + rxt(k,220)*y(k,136) + rxt(k,227)*y(k,60) + rxt(k,236)*y(k,93) &
                      + rxt(k,277)*y(k,41) + rxt(k,279)*y(k,43) + rxt(k,282)*y(k,46) &
                      + rxt(k,284)*y(k,55) + rxt(k,325)*y(k,28) + rxt(k,355)*y(k,31))
         mat(k,1992) = -rxt(k,214)*y(k,56)
         mat(k,2257) = -rxt(k,215)*y(k,56)
         mat(k,1469) = -rxt(k,216)*y(k,56)
         mat(k,586) = -rxt(k,217)*y(k,56)
         mat(k,2365) = -(rxt(k,218) + rxt(k,219)) * y(k,56)
         mat(k,2232) = -rxt(k,220)*y(k,56)
         mat(k,963) = -rxt(k,227)*y(k,56)
         mat(k,842) = -rxt(k,236)*y(k,56)
         mat(k,484) = -rxt(k,277)*y(k,56)
         mat(k,607) = -rxt(k,279)*y(k,56)
         mat(k,375) = -rxt(k,282)*y(k,56)
         mat(k,390) = -rxt(k,284)*y(k,56)
         mat(k,302) = -rxt(k,325)*y(k,56)
         mat(k,308) = -rxt(k,355)*y(k,56)
         mat(k,1576) = rxt(k,255)*y(k,59)
         mat(k,103) = 4.000_r8*rxt(k,239)*y(k,224)
         mat(k,143) = rxt(k,240)*y(k,224)
         mat(k,114) = 2.000_r8*rxt(k,241)*y(k,224)
         mat(k,153) = 2.000_r8*rxt(k,242)*y(k,224)
         mat(k,118) = 2.000_r8*rxt(k,243)*y(k,224)
         mat(k,158) = rxt(k,244)*y(k,224)
         mat(k,122) = 2.000_r8*rxt(k,245)*y(k,224)
         mat(k,125) = 3.000_r8*rxt(k,281)*y(k,228)
         mat(k,375) = mat(k,375) + rxt(k,283)*y(k,228)
         mat(k,1602) = rxt(k,255)*y(k,19) + (4.000_r8*rxt(k,222)+2.000_r8*rxt(k,224)) &
                      *y(k,59) + rxt(k,226)*y(k,124) + rxt(k,231)*y(k,134) &
                      + rxt(k,510)*y(k,152) + rxt(k,221)*y(k,201) + rxt(k,232) &
                      *y(k,228)
         mat(k,249) = rxt(k,276)*y(k,224)
         mat(k,245) = rxt(k,291)*y(k,224) + rxt(k,286)*y(k,228)
         mat(k,255) = rxt(k,292)*y(k,224) + rxt(k,287)*y(k,228)
         mat(k,296) = rxt(k,293)*y(k,224) + rxt(k,288)*y(k,228)
         mat(k,1506) = rxt(k,234)*y(k,134) + rxt(k,246)*y(k,224) + rxt(k,235)*y(k,228)
         mat(k,1921) = rxt(k,226)*y(k,59)
         mat(k,2101) = rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85)
         mat(k,1486) = rxt(k,510)*y(k,59)
         mat(k,2417) = rxt(k,221)*y(k,59)
         mat(k,2035) = 4.000_r8*rxt(k,239)*y(k,33) + rxt(k,240)*y(k,34) &
                      + 2.000_r8*rxt(k,241)*y(k,36) + 2.000_r8*rxt(k,242)*y(k,37) &
                      + 2.000_r8*rxt(k,243)*y(k,38) + rxt(k,244)*y(k,39) &
                      + 2.000_r8*rxt(k,245)*y(k,40) + rxt(k,276)*y(k,65) + rxt(k,291) &
                      *y(k,82) + rxt(k,292)*y(k,83) + rxt(k,293)*y(k,84) + rxt(k,246) &
                      *y(k,85)
         mat(k,1826) = 3.000_r8*rxt(k,281)*y(k,44) + rxt(k,283)*y(k,46) + rxt(k,232) &
                      *y(k,59) + rxt(k,286)*y(k,82) + rxt(k,287)*y(k,83) + rxt(k,288) &
                      *y(k,84) + rxt(k,235)*y(k,85)
         mat(k,2128) = rxt(k,227)*y(k,60)
         mat(k,1581) = 2.000_r8*rxt(k,223)*y(k,59)
         mat(k,953) = rxt(k,227)*y(k,56) + (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,85)
         mat(k,1490) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,60) + (rxt(k,563) &
                       +rxt(k,569)+rxt(k,574))*y(k,93)
         mat(k,836) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,85)
         mat(k,1580) = 2.000_r8*rxt(k,248)*y(k,59)
         mat(k,1592) = -(rxt(k,221)*y(k,201) + (4._r8*rxt(k,222) + 4._r8*rxt(k,223) &
                      + 4._r8*rxt(k,224) + 4._r8*rxt(k,248)) * y(k,59) + rxt(k,225) &
                      *y(k,90) + rxt(k,226)*y(k,124) + rxt(k,228)*y(k,125) + rxt(k,231) &
                      *y(k,134) + (rxt(k,232) + rxt(k,233)) * y(k,228) + (rxt(k,254) &
                      + rxt(k,255) + rxt(k,256)) * y(k,19) + rxt(k,510)*y(k,152))
         mat(k,2407) = -rxt(k,221)*y(k,59)
         mat(k,2355) = -rxt(k,225)*y(k,59)
         mat(k,1911) = -rxt(k,226)*y(k,59)
         mat(k,1956) = -rxt(k,228)*y(k,59)
         mat(k,2091) = -rxt(k,231)*y(k,59)
         mat(k,1816) = -(rxt(k,232) + rxt(k,233)) * y(k,59)
         mat(k,1566) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,59)
         mat(k,1480) = -rxt(k,510)*y(k,59)
         mat(k,2158) = rxt(k,219)*y(k,90) + rxt(k,236)*y(k,93) + rxt(k,220)*y(k,136)
         mat(k,957) = rxt(k,229)*y(k,134)
         mat(k,1498) = rxt(k,247)*y(k,224)
         mat(k,2355) = mat(k,2355) + rxt(k,219)*y(k,56)
         mat(k,839) = rxt(k,236)*y(k,56) + rxt(k,237)*y(k,134) + rxt(k,238)*y(k,228)
         mat(k,2091) = mat(k,2091) + rxt(k,229)*y(k,60) + rxt(k,237)*y(k,93)
         mat(k,2222) = rxt(k,220)*y(k,56)
         mat(k,337) = rxt(k,515)*y(k,152)
         mat(k,1480) = mat(k,1480) + rxt(k,515)*y(k,138)
         mat(k,2025) = rxt(k,247)*y(k,85)
         mat(k,1816) = mat(k,1816) + rxt(k,238)*y(k,93)
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
         mat(k,955) = -(rxt(k,227)*y(k,56) + rxt(k,229)*y(k,134) + rxt(k,230)*y(k,228) &
                      + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,85))
         mat(k,2143) = -rxt(k,227)*y(k,60)
         mat(k,2081) = -rxt(k,229)*y(k,60)
         mat(k,1779) = -rxt(k,230)*y(k,60)
         mat(k,1494) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,60)
         mat(k,1586) = rxt(k,228)*y(k,125)
         mat(k,1942) = rxt(k,228)*y(k,59)
         mat(k,1160) = -(rxt(k,314)*y(k,228))
         mat(k,1793) = -rxt(k,314)*y(k,62)
         mat(k,982) = .230_r8*rxt(k,480)*y(k,136)
         mat(k,1509) = rxt(k,250)*y(k,42)
         mat(k,288) = .350_r8*rxt(k,316)*y(k,228)
         mat(k,553) = .630_r8*rxt(k,318)*y(k,136)
         mat(k,1135) = .560_r8*rxt(k,347)*y(k,136)
         mat(k,1975) = rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) + rxt(k,295)*y(k,126) &
                      + rxt(k,296)*y(k,134) + rxt(k,297)*y(k,228)
         mat(k,372) = rxt(k,282)*y(k,56)
         mat(k,1280) = rxt(k,353)*y(k,126) + rxt(k,354)*y(k,228)
         mat(k,2147) = rxt(k,214)*y(k,42) + rxt(k,282)*y(k,46)
         mat(k,1449) = rxt(k,591)*y(k,229)
         mat(k,1055) = rxt(k,341)*y(k,228)
         mat(k,2334) = .070_r8*rxt(k,450)*y(k,202) + .160_r8*rxt(k,453)*y(k,214) &
                      + .140_r8*rxt(k,456)*y(k,216)
         mat(k,879) = .620_r8*rxt(k,425)*y(k,136)
         mat(k,1268) = .650_r8*rxt(k,378)*y(k,136)
         mat(k,1026) = .230_r8*rxt(k,483)*y(k,136)
         mat(k,1372) = .560_r8*rxt(k,392)*y(k,136)
         mat(k,1891) = .170_r8*rxt(k,451)*y(k,202) + .220_r8*rxt(k,376)*y(k,213) &
                      + .400_r8*rxt(k,454)*y(k,214) + .350_r8*rxt(k,457)*y(k,216) &
                      + .225_r8*rxt(k,492)*y(k,233) + .250_r8*rxt(k,433)*y(k,237)
         mat(k,1630) = rxt(k,295)*y(k,42) + rxt(k,353)*y(k,49) + .220_r8*rxt(k,375) &
                      *y(k,213) + .500_r8*rxt(k,434)*y(k,237)
         mat(k,2083) = rxt(k,296)*y(k,42) + rxt(k,504)*y(k,139)
         mat(k,2203) = .230_r8*rxt(k,480)*y(k,6) + .630_r8*rxt(k,318)*y(k,25) &
                      + .560_r8*rxt(k,347)*y(k,29) + .620_r8*rxt(k,425)*y(k,99) &
                      + .650_r8*rxt(k,378)*y(k,105) + .230_r8*rxt(k,483)*y(k,110) &
                      + .560_r8*rxt(k,392)*y(k,111)
         mat(k,366) = rxt(k,504)*y(k,134) + rxt(k,505)*y(k,228)
         mat(k,1120) = .700_r8*rxt(k,501)*y(k,228)
         mat(k,1416) = .220_r8*rxt(k,372)*y(k,213) + .250_r8*rxt(k,430)*y(k,237)
         mat(k,2389) = .110_r8*rxt(k,373)*y(k,213) + .125_r8*rxt(k,490)*y(k,233) &
                      + .200_r8*rxt(k,431)*y(k,237)
         mat(k,779) = .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451)*y(k,124)
         mat(k,1350) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,200) + .110_r8*rxt(k,373)*y(k,201)
         mat(k,742) = .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454)*y(k,124)
         mat(k,912) = .140_r8*rxt(k,456)*y(k,90) + .350_r8*rxt(k,457)*y(k,124)
         mat(k,1793) = mat(k,1793) + .350_r8*rxt(k,316)*y(k,24) + rxt(k,297)*y(k,42) &
                      + rxt(k,354)*y(k,49) + rxt(k,341)*y(k,75) + rxt(k,505)*y(k,139) &
                      + .700_r8*rxt(k,501)*y(k,181)
         mat(k,809) = rxt(k,591)*y(k,63)
         mat(k,1176) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,201)
         mat(k,1230) = .250_r8*rxt(k,433)*y(k,124) + .500_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,200) + .200_r8*rxt(k,431)*y(k,201)
         mat(k,1450) = -(rxt(k,591)*y(k,229))
         mat(k,810) = -rxt(k,591)*y(k,63)
         mat(k,986) = .270_r8*rxt(k,480)*y(k,136)
         mat(k,1139) = .200_r8*rxt(k,347)*y(k,136)
         mat(k,706) = rxt(k,334)*y(k,228)
         mat(k,622) = .500_r8*rxt(k,335)*y(k,228)
         mat(k,1161) = rxt(k,314)*y(k,228)
         mat(k,1167) = .800_r8*rxt(k,340)*y(k,228)
         mat(k,1056) = rxt(k,341)*y(k,228)
         mat(k,920) = rxt(k,306)*y(k,228)
         mat(k,2349) = .450_r8*rxt(k,383)*y(k,215)
         mat(k,614) = .500_r8*rxt(k,391)*y(k,228)
         mat(k,1030) = .270_r8*rxt(k,483)*y(k,136)
         mat(k,1379) = .100_r8*rxt(k,392)*y(k,136)
         mat(k,1907) = rxt(k,333)*y(k,200) + .900_r8*rxt(k,492)*y(k,233)
         mat(k,2217) = .270_r8*rxt(k,480)*y(k,6) + .200_r8*rxt(k,347)*y(k,29) &
                      + .270_r8*rxt(k,483)*y(k,110) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,1123) = 1.800_r8*rxt(k,501)*y(k,228)
         mat(k,1429) = rxt(k,333)*y(k,124) + 4.000_r8*rxt(k,330)*y(k,200) &
                      + .900_r8*rxt(k,331)*y(k,201) + rxt(k,405)*y(k,208) &
                      + 2.000_r8*rxt(k,381)*y(k,215) + rxt(k,430)*y(k,237)
         mat(k,2404) = .900_r8*rxt(k,331)*y(k,200) + rxt(k,382)*y(k,215) &
                      + .500_r8*rxt(k,490)*y(k,233)
         mat(k,1304) = rxt(k,405)*y(k,200)
         mat(k,1399) = .450_r8*rxt(k,383)*y(k,90) + 2.000_r8*rxt(k,381)*y(k,200) &
                      + rxt(k,382)*y(k,201) + 4.000_r8*rxt(k,384)*y(k,215)
         mat(k,1809) = rxt(k,334)*y(k,50) + .500_r8*rxt(k,335)*y(k,51) + rxt(k,314) &
                      *y(k,62) + .800_r8*rxt(k,340)*y(k,74) + rxt(k,341)*y(k,75) &
                      + rxt(k,306)*y(k,87) + .500_r8*rxt(k,391)*y(k,109) &
                      + 1.800_r8*rxt(k,501)*y(k,181)
         mat(k,1181) = .900_r8*rxt(k,492)*y(k,124) + .500_r8*rxt(k,490)*y(k,201)
         mat(k,1236) = rxt(k,430)*y(k,200)
         mat(k,238) = -(rxt(k,275)*y(k,224))
         mat(k,2009) = -rxt(k,275)*y(k,64)
         mat(k,140) = rxt(k,240)*y(k,224)
         mat(k,145) = rxt(k,266)*y(k,224)
         mat(k,150) = rxt(k,242)*y(k,224)
         mat(k,116) = 2.000_r8*rxt(k,243)*y(k,224)
         mat(k,155) = 2.000_r8*rxt(k,244)*y(k,224)
         mat(k,120) = rxt(k,245)*y(k,224)
         mat(k,108) = 2.000_r8*rxt(k,268)*y(k,224)
         mat(k,250) = rxt(k,292)*y(k,224) + rxt(k,287)*y(k,228)
         mat(k,291) = rxt(k,293)*y(k,224) + rxt(k,288)*y(k,228)
         mat(k,2009) = mat(k,2009) + rxt(k,240)*y(k,34) + rxt(k,266)*y(k,35) &
                      + rxt(k,242)*y(k,37) + 2.000_r8*rxt(k,243)*y(k,38) &
                      + 2.000_r8*rxt(k,244)*y(k,39) + rxt(k,245)*y(k,40) &
                      + 2.000_r8*rxt(k,268)*y(k,78) + rxt(k,292)*y(k,83) + rxt(k,293) &
                      *y(k,84)
         mat(k,1699) = rxt(k,287)*y(k,83) + rxt(k,288)*y(k,84)
         mat(k,246) = -(rxt(k,276)*y(k,224))
         mat(k,2011) = -rxt(k,276)*y(k,65)
         mat(k,112) = rxt(k,241)*y(k,224)
         mat(k,151) = rxt(k,242)*y(k,224)
         mat(k,242) = rxt(k,291)*y(k,224) + rxt(k,286)*y(k,228)
         mat(k,2011) = mat(k,2011) + rxt(k,241)*y(k,36) + rxt(k,242)*y(k,37) &
                      + rxt(k,291)*y(k,82)
         mat(k,1701) = rxt(k,286)*y(k,82)
         mat(k,194) = -(rxt(k,449)*y(k,228))
         mat(k,1690) = -rxt(k,449)*y(k,66)
         mat(k,188) = .180_r8*rxt(k,469)*y(k,228)
         mat(k,1690) = mat(k,1690) + .180_r8*rxt(k,469)*y(k,183)
         mat(k,310) = -(rxt(k,502)*y(k,126) + (rxt(k,503) + rxt(k,517)) * y(k,228))
         mat(k,1610) = -rxt(k,502)*y(k,67)
         mat(k,1710) = -(rxt(k,503) + rxt(k,517)) * y(k,67)
         mat(k,2275) = rxt(k,336)*y(k,206)
         mat(k,795) = rxt(k,336)*y(k,90)
         mat(k,925) = -(rxt(k,271)*y(k,54) + rxt(k,272)*y(k,77) + rxt(k,273)*y(k,241) &
                      + rxt(k,274)*y(k,89))
         mat(k,2239) = -rxt(k,271)*y(k,73)
         mat(k,1460) = -rxt(k,272)*y(k,73)
         mat(k,2426) = -rxt(k,273)*y(k,73)
         mat(k,2042) = -rxt(k,274)*y(k,73)
         mat(k,146) = rxt(k,266)*y(k,224)
         mat(k,156) = rxt(k,244)*y(k,224)
         mat(k,239) = 2.000_r8*rxt(k,275)*y(k,224)
         mat(k,247) = rxt(k,276)*y(k,224)
         mat(k,2019) = rxt(k,266)*y(k,35) + rxt(k,244)*y(k,39) + 2.000_r8*rxt(k,275) &
                      *y(k,64) + rxt(k,276)*y(k,65)
         mat(k,1166) = -(rxt(k,340)*y(k,228))
         mat(k,1794) = -rxt(k,340)*y(k,74)
         mat(k,590) = .700_r8*rxt(k,416)*y(k,228)
         mat(k,560) = .500_r8*rxt(k,417)*y(k,228)
         mat(k,429) = rxt(k,428)*y(k,228)
         mat(k,1892) = .050_r8*rxt(k,414)*y(k,209) + .530_r8*rxt(k,376)*y(k,213) &
                      + .225_r8*rxt(k,492)*y(k,233) + .250_r8*rxt(k,433)*y(k,237)
         mat(k,1631) = .050_r8*rxt(k,415)*y(k,209) + .530_r8*rxt(k,375)*y(k,213) &
                      + .250_r8*rxt(k,434)*y(k,237)
         mat(k,1538) = rxt(k,339)*y(k,205)
         mat(k,1417) = .530_r8*rxt(k,372)*y(k,213) + .250_r8*rxt(k,430)*y(k,237)
         mat(k,2390) = .260_r8*rxt(k,373)*y(k,213) + .125_r8*rxt(k,490)*y(k,233) &
                      + .100_r8*rxt(k,431)*y(k,237)
         mat(k,461) = rxt(k,339)*y(k,135)
         mat(k,1325) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1351) = .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375)*y(k,126) &
                      + .530_r8*rxt(k,372)*y(k,200) + .260_r8*rxt(k,373)*y(k,201)
         mat(k,1794) = mat(k,1794) + .700_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101) + rxt(k,428)*y(k,115)
         mat(k,1177) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,201)
         mat(k,1231) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,200) + .100_r8*rxt(k,431)*y(k,201)
         mat(k,1054) = -(rxt(k,341)*y(k,228))
         mat(k,1784) = -rxt(k,341)*y(k,75)
         mat(k,287) = .650_r8*rxt(k,316)*y(k,228)
         mat(k,1164) = .200_r8*rxt(k,340)*y(k,228)
         mat(k,2327) = .160_r8*rxt(k,453)*y(k,214) + .070_r8*rxt(k,456)*y(k,216)
         mat(k,1089) = rxt(k,429)*y(k,228)
         mat(k,1883) = rxt(k,440)*y(k,194) + .050_r8*rxt(k,414)*y(k,209) &
                      + .400_r8*rxt(k,454)*y(k,214) + .170_r8*rxt(k,457)*y(k,216) &
                      + .700_r8*rxt(k,460)*y(k,230) + .600_r8*rxt(k,467)*y(k,235) &
                      + .250_r8*rxt(k,433)*y(k,237) + .340_r8*rxt(k,473)*y(k,238) &
                      + .170_r8*rxt(k,476)*y(k,240)
         mat(k,1621) = .050_r8*rxt(k,415)*y(k,209) + .250_r8*rxt(k,434)*y(k,237)
         mat(k,497) = rxt(k,440)*y(k,124)
         mat(k,1414) = .250_r8*rxt(k,430)*y(k,237)
         mat(k,2381) = .100_r8*rxt(k,431)*y(k,237)
         mat(k,1323) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,741) = .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454)*y(k,124)
         mat(k,911) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,1784) = mat(k,1784) + .650_r8*rxt(k,316)*y(k,24) + .200_r8*rxt(k,340) &
                      *y(k,74) + rxt(k,429)*y(k,116)
         mat(k,449) = .700_r8*rxt(k,460)*y(k,124)
         mat(k,754) = .600_r8*rxt(k,467)*y(k,124)
         mat(k,1229) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,200) + .100_r8*rxt(k,431)*y(k,201)
         mat(k,770) = .340_r8*rxt(k,473)*y(k,124)
         mat(k,511) = .170_r8*rxt(k,476)*y(k,124)
         mat(k,2121) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,90) + rxt(k,175) &
                      *y(k,135) + rxt(k,178)*y(k,136))
         mat(k,2364) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76)
         mat(k,1551) = -rxt(k,175)*y(k,76)
         mat(k,2231) = -rxt(k,178)*y(k,76)
         mat(k,1991) = rxt(k,297)*y(k,228)
         mat(k,2256) = rxt(k,311)*y(k,224)
         mat(k,2167) = rxt(k,216)*y(k,77)
         mat(k,930) = rxt(k,272)*y(k,77)
         mat(k,1468) = rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73) + rxt(k,170)*y(k,134) &
                      + rxt(k,153)*y(k,224) + rxt(k,179)*y(k,228)
         mat(k,834) = rxt(k,270)*y(k,224)
         mat(k,1505) = rxt(k,247)*y(k,224)
         mat(k,1008) = rxt(k,202)*y(k,228)
         mat(k,2100) = rxt(k,170)*y(k,77) + rxt(k,182)*y(k,228)
         mat(k,370) = rxt(k,505)*y(k,228)
         mat(k,723) = rxt(k,511)*y(k,228)
         mat(k,1485) = rxt(k,516)*y(k,228)
         mat(k,2034) = rxt(k,311)*y(k,54) + rxt(k,153)*y(k,77) + rxt(k,270)*y(k,81) &
                      + rxt(k,247)*y(k,85)
         mat(k,1825) = rxt(k,297)*y(k,42) + rxt(k,179)*y(k,77) + rxt(k,202)*y(k,112) &
                      + rxt(k,182)*y(k,134) + rxt(k,505)*y(k,139) + rxt(k,511) &
                      *y(k,150) + rxt(k,516)*y(k,152)
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
         mat(k,1461) = -(rxt(k,153)*y(k,224) + rxt(k,170)*y(k,134) + rxt(k,179) &
                      *y(k,228) + rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73))
         mat(k,2020) = -rxt(k,153)*y(k,77)
         mat(k,2085) = -rxt(k,170)*y(k,77)
         mat(k,1810) = -rxt(k,179)*y(k,77)
         mat(k,2153) = -rxt(k,216)*y(k,77)
         mat(k,926) = -rxt(k,272)*y(k,77)
         mat(k,2242) = rxt(k,312)*y(k,224)
         mat(k,2107) = rxt(k,172)*y(k,90)
         mat(k,2350) = rxt(k,172)*y(k,76)
         mat(k,2020) = mat(k,2020) + rxt(k,312)*y(k,54)
         mat(k,107) = -(rxt(k,268)*y(k,224))
         mat(k,1999) = -rxt(k,268)*y(k,78)
         mat(k,582) = -(rxt(k,171)*y(k,134) + rxt(k,180)*y(k,228) + rxt(k,217)*y(k,56))
         mat(k,2070) = -rxt(k,171)*y(k,79)
         mat(k,1745) = -rxt(k,180)*y(k,79)
         mat(k,2136) = -rxt(k,217)*y(k,79)
         mat(k,2302) = 2.000_r8*rxt(k,186)*y(k,90)
         mat(k,1745) = mat(k,1745) + 2.000_r8*rxt(k,185)*y(k,228)
         mat(k,260) = rxt(k,518)*y(k,241)
         mat(k,2423) = rxt(k,518)*y(k,154)
         mat(k,828) = -(rxt(k,263)*y(k,134) + rxt(k,264)*y(k,228) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,224))
         mat(k,2075) = -rxt(k,263)*y(k,81)
         mat(k,1769) = -rxt(k,264)*y(k,81)
         mat(k,2018) = -(rxt(k,269) + rxt(k,270)) * y(k,81)
         mat(k,1508) = rxt(k,250)*y(k,42) + rxt(k,251)*y(k,90)
         mat(k,1973) = rxt(k,250)*y(k,17)
         mat(k,2319) = rxt(k,251)*y(k,17)
         mat(k,241) = -(rxt(k,286)*y(k,228) + rxt(k,291)*y(k,224))
         mat(k,1700) = -rxt(k,286)*y(k,82)
         mat(k,2010) = -rxt(k,291)*y(k,82)
         mat(k,251) = -(rxt(k,287)*y(k,228) + rxt(k,292)*y(k,224))
         mat(k,1702) = -rxt(k,287)*y(k,83)
         mat(k,2012) = -rxt(k,292)*y(k,83)
         mat(k,292) = -(rxt(k,288)*y(k,228) + rxt(k,293)*y(k,224))
         mat(k,1707) = -rxt(k,288)*y(k,84)
         mat(k,2014) = -rxt(k,293)*y(k,84)
         mat(k,1495) = -(rxt(k,234)*y(k,134) + rxt(k,235)*y(k,228) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,224) + (rxt(k,563) + rxt(k,569) + rxt(k,574) &
                      ) * y(k,93) + (rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,60) &
                      + (rxt(k,570) + rxt(k,575)) * y(k,92))
         mat(k,2087) = -rxt(k,234)*y(k,85)
         mat(k,1812) = -rxt(k,235)*y(k,85)
         mat(k,2021) = -(rxt(k,246) + rxt(k,247)) * y(k,85)
         mat(k,838) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,85)
         mat(k,956) = -(rxt(k,568) + rxt(k,573) + rxt(k,578)) * y(k,85)
         mat(k,788) = -(rxt(k,570) + rxt(k,575)) * y(k,85)
         mat(k,300) = rxt(k,325)*y(k,56)
         mat(k,306) = rxt(k,355)*y(k,56)
         mat(k,480) = rxt(k,277)*y(k,56)
         mat(k,1978) = rxt(k,214)*y(k,56)
         mat(k,603) = rxt(k,279)*y(k,56)
         mat(k,373) = 2.000_r8*rxt(k,282)*y(k,56)
         mat(k,2243) = rxt(k,215)*y(k,56)
         mat(k,386) = rxt(k,284)*y(k,56)
         mat(k,2154) = rxt(k,325)*y(k,28) + rxt(k,355)*y(k,31) + rxt(k,277)*y(k,41) &
                      + rxt(k,214)*y(k,42) + rxt(k,279)*y(k,43) + 2.000_r8*rxt(k,282) &
                      *y(k,46) + rxt(k,215)*y(k,54) + rxt(k,284)*y(k,55) + rxt(k,216) &
                      *y(k,77) + rxt(k,217)*y(k,79) + rxt(k,218)*y(k,90) + rxt(k,236) &
                      *y(k,93)
         mat(k,1588) = rxt(k,233)*y(k,228)
         mat(k,1462) = rxt(k,216)*y(k,56)
         mat(k,583) = rxt(k,217)*y(k,56)
         mat(k,2351) = rxt(k,218)*y(k,56)
         mat(k,838) = mat(k,838) + rxt(k,236)*y(k,56)
         mat(k,1812) = mat(k,1812) + rxt(k,233)*y(k,59)
         mat(k,182) = -(rxt(k,305)*y(k,228) + rxt(k,313)*y(k,224))
         mat(k,1688) = -rxt(k,305)*y(k,86)
         mat(k,2008) = -rxt(k,313)*y(k,86)
         mat(k,919) = -(rxt(k,306)*y(k,228))
         mat(k,1775) = -rxt(k,306)*y(k,87)
         mat(k,973) = .050_r8*rxt(k,480)*y(k,136)
         mat(k,286) = .350_r8*rxt(k,316)*y(k,228)
         mat(k,552) = .370_r8*rxt(k,318)*y(k,136)
         mat(k,1132) = .120_r8*rxt(k,347)*y(k,136)
         mat(k,2324) = rxt(k,307)*y(k,207)
         mat(k,877) = .110_r8*rxt(k,425)*y(k,136)
         mat(k,1267) = .330_r8*rxt(k,378)*y(k,136)
         mat(k,1017) = .050_r8*rxt(k,483)*y(k,136)
         mat(k,1369) = .120_r8*rxt(k,392)*y(k,136)
         mat(k,1879) = rxt(k,309)*y(k,207)
         mat(k,2190) = .050_r8*rxt(k,480)*y(k,6) + .370_r8*rxt(k,318)*y(k,25) &
                      + .120_r8*rxt(k,347)*y(k,29) + .110_r8*rxt(k,425)*y(k,99) &
                      + .330_r8*rxt(k,378)*y(k,105) + .050_r8*rxt(k,483)*y(k,110) &
                      + .120_r8*rxt(k,392)*y(k,111)
         mat(k,442) = rxt(k,307)*y(k,90) + rxt(k,309)*y(k,124)
         mat(k,1775) = mat(k,1775) + .350_r8*rxt(k,316)*y(k,24)
         mat(k,2238) = rxt(k,271)*y(k,73)
         mat(k,924) = rxt(k,271)*y(k,54) + rxt(k,272)*y(k,77) + rxt(k,274)*y(k,89) &
                      + rxt(k,273)*y(k,241)
         mat(k,1459) = rxt(k,272)*y(k,73)
         mat(k,2041) = rxt(k,274)*y(k,73)
         mat(k,2425) = rxt(k,273)*y(k,73)
         mat(k,2055) = -(rxt(k,211)*y(k,228) + rxt(k,274)*y(k,73))
         mat(k,1823) = -rxt(k,211)*y(k,89)
         mat(k,929) = -rxt(k,274)*y(k,89)
         mat(k,1989) = rxt(k,295)*y(k,126)
         mat(k,1156) = rxt(k,327)*y(k,126)
         mat(k,1286) = rxt(k,353)*y(k,126)
         mat(k,961) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,85)
         mat(k,314) = rxt(k,502)*y(k,126)
         mat(k,1503) = (rxt(k,568)+rxt(k,573)+rxt(k,578))*y(k,60)
         mat(k,1963) = rxt(k,210)*y(k,228)
         mat(k,1658) = rxt(k,295)*y(k,42) + rxt(k,327)*y(k,45) + rxt(k,353)*y(k,49) &
                      + rxt(k,502)*y(k,67)
         mat(k,1823) = mat(k,1823) + rxt(k,210)*y(k,125)
         mat(k,2368) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76) + rxt(k,176) &
                      *y(k,134) + rxt(k,177)*y(k,136) + rxt(k,181)*y(k,228) &
                      + 4._r8*rxt(k,186)*y(k,90) + rxt(k,198)*y(k,126) + rxt(k,203) &
                      *y(k,124) + rxt(k,208)*y(k,125) + (rxt(k,218) + rxt(k,219) &
                      ) * y(k,56) + rxt(k,225)*y(k,59) + rxt(k,251)*y(k,17) + rxt(k,257) &
                      *y(k,19) + rxt(k,294)*y(k,42) + rxt(k,300)*y(k,201) + rxt(k,307) &
                      *y(k,207) + rxt(k,321)*y(k,197) + rxt(k,332)*y(k,200) + rxt(k,336) &
                      *y(k,206) + rxt(k,349)*y(k,198) + rxt(k,358)*y(k,231) + rxt(k,362) &
                      *y(k,232) + rxt(k,374)*y(k,213) + rxt(k,383)*y(k,215) + rxt(k,387) &
                      *y(k,217) + rxt(k,397)*y(k,192) + rxt(k,407)*y(k,208) + rxt(k,412) &
                      *y(k,209) + rxt(k,421)*y(k,210) + rxt(k,432)*y(k,237) + rxt(k,436) &
                      *y(k,191) + rxt(k,439)*y(k,194) + rxt(k,443)*y(k,196) + rxt(k,446) &
                      *y(k,199) + rxt(k,450)*y(k,202) + rxt(k,453)*y(k,214) + rxt(k,456) &
                      *y(k,216) + rxt(k,459)*y(k,230) + rxt(k,466)*y(k,235) + rxt(k,472) &
                      *y(k,238) + rxt(k,475)*y(k,240) + rxt(k,486)*y(k,223) + rxt(k,491) &
                      *y(k,233) + rxt(k,496)*y(k,234))
         mat(k,2125) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,90)
         mat(k,2104) = -rxt(k,176)*y(k,90)
         mat(k,2235) = -rxt(k,177)*y(k,90)
         mat(k,1829) = -rxt(k,181)*y(k,90)
         mat(k,1664) = -rxt(k,198)*y(k,90)
         mat(k,1924) = -rxt(k,203)*y(k,90)
         mat(k,1969) = -rxt(k,208)*y(k,90)
         mat(k,2171) = -(rxt(k,218) + rxt(k,219)) * y(k,90)
         mat(k,1604) = -rxt(k,225)*y(k,90)
         mat(k,1522) = -rxt(k,251)*y(k,90)
         mat(k,1578) = -rxt(k,257)*y(k,90)
         mat(k,1995) = -rxt(k,294)*y(k,90)
         mat(k,2420) = -rxt(k,300)*y(k,90)
         mat(k,446) = -rxt(k,307)*y(k,90)
         mat(k,908) = -rxt(k,321)*y(k,90)
         mat(k,1439) = -rxt(k,332)*y(k,90)
         mat(k,804) = -rxt(k,336)*y(k,90)
         mat(k,943) = -rxt(k,349)*y(k,90)
         mat(k,822) = -rxt(k,358)*y(k,90)
         mat(k,1224) = -rxt(k,362)*y(k,90)
         mat(k,1366) = -rxt(k,374)*y(k,90)
         mat(k,1407) = -rxt(k,383)*y(k,90)
         mat(k,704) = -rxt(k,387)*y(k,90)
         mat(k,1051) = -rxt(k,397)*y(k,90)
         mat(k,1312) = -rxt(k,407)*y(k,90)
         mat(k,1345) = -rxt(k,412)*y(k,90)
         mat(k,1265) = -rxt(k,421)*y(k,90)
         mat(k,1242) = -rxt(k,432)*y(k,90)
         mat(k,528) = -rxt(k,436)*y(k,90)
         mat(k,501) = -rxt(k,439)*y(k,90)
         mat(k,440) = -rxt(k,443)*y(k,90)
         mat(k,640) = -rxt(k,446)*y(k,90)
         mat(k,784) = -rxt(k,450)*y(k,90)
         mat(k,745) = -rxt(k,453)*y(k,90)
         mat(k,917) = -rxt(k,456)*y(k,90)
         mat(k,453) = -rxt(k,459)*y(k,90)
         mat(k,760) = -rxt(k,466)*y(k,90)
         mat(k,777) = -rxt(k,472)*y(k,90)
         mat(k,516) = -rxt(k,475)*y(k,90)
         mat(k,1116) = -rxt(k,486)*y(k,90)
         mat(k,1187) = -rxt(k,491)*y(k,90)
         mat(k,1070) = -rxt(k,496)*y(k,90)
         mat(k,991) = .570_r8*rxt(k,480)*y(k,136)
         mat(k,165) = .650_r8*rxt(k,438)*y(k,228)
         mat(k,1522) = mat(k,1522) + rxt(k,250)*y(k,42)
         mat(k,1578) = mat(k,1578) + rxt(k,262)*y(k,228)
         mat(k,290) = .350_r8*rxt(k,316)*y(k,228)
         mat(k,557) = .130_r8*rxt(k,318)*y(k,136)
         mat(k,268) = rxt(k,323)*y(k,228)
         mat(k,1148) = .280_r8*rxt(k,347)*y(k,136)
         mat(k,1995) = mat(k,1995) + rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,134)
         mat(k,608) = rxt(k,279)*y(k,56) + rxt(k,280)*y(k,228)
         mat(k,376) = rxt(k,282)*y(k,56) + rxt(k,283)*y(k,228)
         mat(k,106) = rxt(k,329)*y(k,228)
         mat(k,827) = rxt(k,302)*y(k,228)
         mat(k,2260) = rxt(k,311)*y(k,224)
         mat(k,2171) = mat(k,2171) + rxt(k,214)*y(k,42) + rxt(k,279)*y(k,43) &
                      + rxt(k,282)*y(k,46) + rxt(k,217)*y(k,79)
         mat(k,1604) = mat(k,1604) + rxt(k,221)*y(k,201) + rxt(k,232)*y(k,228)
         mat(k,1163) = rxt(k,314)*y(k,228)
         mat(k,198) = .730_r8*rxt(k,449)*y(k,228)
         mat(k,315) = .500_r8*rxt(k,517)*y(k,228)
         mat(k,1170) = rxt(k,340)*y(k,228)
         mat(k,1058) = rxt(k,341)*y(k,228)
         mat(k,2125) = mat(k,2125) + rxt(k,175)*y(k,135)
         mat(k,587) = rxt(k,217)*y(k,56) + rxt(k,171)*y(k,134) + rxt(k,180)*y(k,228)
         mat(k,185) = rxt(k,305)*y(k,228)
         mat(k,922) = rxt(k,306)*y(k,228)
         mat(k,2368) = mat(k,2368) + .070_r8*rxt(k,450)*y(k,202) + .160_r8*rxt(k,453) &
                      *y(k,214) + .330_r8*rxt(k,456)*y(k,216)
         mat(k,1205) = rxt(k,371)*y(k,228)
         mat(k,1213) = rxt(k,356)*y(k,228)
         mat(k,890) = .370_r8*rxt(k,425)*y(k,136)
         mat(k,597) = .300_r8*rxt(k,416)*y(k,228)
         mat(k,565) = rxt(k,417)*y(k,228)
         mat(k,408) = rxt(k,424)*y(k,228)
         mat(k,1278) = .140_r8*rxt(k,378)*y(k,136)
         mat(k,320) = .200_r8*rxt(k,380)*y(k,228)
         mat(k,619) = .500_r8*rxt(k,391)*y(k,228)
         mat(k,1035) = .570_r8*rxt(k,483)*y(k,136)
         mat(k,1389) = .280_r8*rxt(k,392)*y(k,136)
         mat(k,432) = rxt(k,428)*y(k,228)
         mat(k,1100) = rxt(k,429)*y(k,228)
         mat(k,1924) = mat(k,1924) + rxt(k,398)*y(k,192) + rxt(k,440)*y(k,194) &
                      + rxt(k,445)*y(k,196) + rxt(k,322)*y(k,197) + rxt(k,350) &
                      *y(k,198) + rxt(k,301)*y(k,201) + .170_r8*rxt(k,451)*y(k,202) &
                      + rxt(k,369)*y(k,204) + .250_r8*rxt(k,337)*y(k,206) + rxt(k,309) &
                      *y(k,207) + .920_r8*rxt(k,408)*y(k,208) + .920_r8*rxt(k,414) &
                      *y(k,209) + rxt(k,422)*y(k,210) + .470_r8*rxt(k,376)*y(k,213) &
                      + .400_r8*rxt(k,454)*y(k,214) + .830_r8*rxt(k,457)*y(k,216) &
                      + rxt(k,460)*y(k,230) + rxt(k,359)*y(k,231) + .900_r8*rxt(k,492) &
                      *y(k,233) + .800_r8*rxt(k,497)*y(k,234) + rxt(k,467)*y(k,235) &
                      + rxt(k,433)*y(k,237) + rxt(k,473)*y(k,238) + rxt(k,476) &
                      *y(k,240)
         mat(k,1664) = mat(k,1664) + rxt(k,295)*y(k,42) + rxt(k,409)*y(k,208) &
                      + rxt(k,415)*y(k,209) + rxt(k,423)*y(k,210) + .470_r8*rxt(k,375) &
                      *y(k,213) + rxt(k,201)*y(k,228) + rxt(k,434)*y(k,237)
         mat(k,2104) = mat(k,2104) + rxt(k,296)*y(k,42) + rxt(k,171)*y(k,79)
         mat(k,1554) = rxt(k,175)*y(k,76) + rxt(k,339)*y(k,205)
         mat(k,2235) = mat(k,2235) + .570_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,318) &
                      *y(k,25) + .280_r8*rxt(k,347)*y(k,29) + .370_r8*rxt(k,425) &
                      *y(k,99) + .140_r8*rxt(k,378)*y(k,105) + .570_r8*rxt(k,483) &
                      *y(k,110) + .280_r8*rxt(k,392)*y(k,111) + rxt(k,183)*y(k,228)
         mat(k,174) = .800_r8*rxt(k,461)*y(k,228)
         mat(k,951) = rxt(k,507)*y(k,228)
         mat(k,1126) = .200_r8*rxt(k,501)*y(k,228)
         mat(k,193) = .280_r8*rxt(k,469)*y(k,228)
         mat(k,215) = .380_r8*rxt(k,471)*y(k,228)
         mat(k,220) = .630_r8*rxt(k,477)*y(k,228)
         mat(k,1051) = mat(k,1051) + rxt(k,398)*y(k,124)
         mat(k,501) = mat(k,501) + rxt(k,440)*y(k,124)
         mat(k,440) = mat(k,440) + rxt(k,445)*y(k,124)
         mat(k,908) = mat(k,908) + rxt(k,322)*y(k,124) + 2.400_r8*rxt(k,319)*y(k,197) &
                      + rxt(k,320)*y(k,201)
         mat(k,943) = mat(k,943) + rxt(k,350)*y(k,124) + rxt(k,348)*y(k,201)
         mat(k,1439) = mat(k,1439) + .900_r8*rxt(k,331)*y(k,201) + rxt(k,405)*y(k,208) &
                      + rxt(k,410)*y(k,209) + rxt(k,419)*y(k,210) + .470_r8*rxt(k,372) &
                      *y(k,213) + rxt(k,430)*y(k,237)
         mat(k,2420) = mat(k,2420) + rxt(k,221)*y(k,59) + rxt(k,301)*y(k,124) &
                      + rxt(k,320)*y(k,197) + rxt(k,348)*y(k,198) + .900_r8*rxt(k,331) &
                      *y(k,200) + 4.000_r8*rxt(k,298)*y(k,201) + rxt(k,406)*y(k,208) &
                      + rxt(k,411)*y(k,209) + 1.200_r8*rxt(k,420)*y(k,210) &
                      + .730_r8*rxt(k,373)*y(k,213) + rxt(k,382)*y(k,215) &
                      + .500_r8*rxt(k,485)*y(k,223) + .300_r8*rxt(k,361)*y(k,232) &
                      + rxt(k,490)*y(k,233) + rxt(k,495)*y(k,234) + .800_r8*rxt(k,431) &
                      *y(k,237)
         mat(k,784) = mat(k,784) + .070_r8*rxt(k,450)*y(k,90) + .170_r8*rxt(k,451) &
                      *y(k,124)
         mat(k,581) = rxt(k,369)*y(k,124)
         mat(k,464) = rxt(k,339)*y(k,135)
         mat(k,804) = mat(k,804) + .250_r8*rxt(k,337)*y(k,124)
         mat(k,446) = mat(k,446) + rxt(k,309)*y(k,124)
         mat(k,1312) = mat(k,1312) + .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) &
                      + rxt(k,405)*y(k,200) + rxt(k,406)*y(k,201)
         mat(k,1345) = mat(k,1345) + .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,200) + rxt(k,411)*y(k,201)
         mat(k,1265) = mat(k,1265) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) &
                      + rxt(k,419)*y(k,200) + 1.200_r8*rxt(k,420)*y(k,201)
         mat(k,1366) = mat(k,1366) + .470_r8*rxt(k,376)*y(k,124) + .470_r8*rxt(k,375) &
                      *y(k,126) + .470_r8*rxt(k,372)*y(k,200) + .730_r8*rxt(k,373) &
                      *y(k,201)
         mat(k,745) = mat(k,745) + .160_r8*rxt(k,453)*y(k,90) + .400_r8*rxt(k,454) &
                      *y(k,124)
         mat(k,1407) = mat(k,1407) + rxt(k,382)*y(k,201)
         mat(k,917) = mat(k,917) + .330_r8*rxt(k,456)*y(k,90) + .830_r8*rxt(k,457) &
                      *y(k,124)
         mat(k,1116) = mat(k,1116) + .500_r8*rxt(k,485)*y(k,201)
         mat(k,2038) = rxt(k,311)*y(k,54)
         mat(k,1829) = mat(k,1829) + .650_r8*rxt(k,438)*y(k,7) + rxt(k,262)*y(k,19) &
                      + .350_r8*rxt(k,316)*y(k,24) + rxt(k,323)*y(k,26) + rxt(k,280) &
                      *y(k,43) + rxt(k,283)*y(k,46) + rxt(k,329)*y(k,47) + rxt(k,302) &
                      *y(k,52) + rxt(k,232)*y(k,59) + rxt(k,314)*y(k,62) &
                      + .730_r8*rxt(k,449)*y(k,66) + .500_r8*rxt(k,517)*y(k,67) &
                      + rxt(k,340)*y(k,74) + rxt(k,341)*y(k,75) + rxt(k,180)*y(k,79) &
                      + rxt(k,305)*y(k,86) + rxt(k,306)*y(k,87) + rxt(k,371)*y(k,94) &
                      + rxt(k,356)*y(k,96) + .300_r8*rxt(k,416)*y(k,100) + rxt(k,417) &
                      *y(k,101) + rxt(k,424)*y(k,102) + .200_r8*rxt(k,380)*y(k,106) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,428)*y(k,115) + rxt(k,429) &
                      *y(k,116) + rxt(k,201)*y(k,126) + rxt(k,183)*y(k,136) &
                      + .800_r8*rxt(k,461)*y(k,144) + rxt(k,507)*y(k,153) &
                      + .200_r8*rxt(k,501)*y(k,181) + .280_r8*rxt(k,469)*y(k,183) &
                      + .380_r8*rxt(k,471)*y(k,185) + .630_r8*rxt(k,477)*y(k,187)
         mat(k,453) = mat(k,453) + rxt(k,460)*y(k,124)
         mat(k,822) = mat(k,822) + rxt(k,359)*y(k,124)
         mat(k,1224) = mat(k,1224) + .300_r8*rxt(k,361)*y(k,201)
         mat(k,1187) = mat(k,1187) + .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,201)
         mat(k,1070) = mat(k,1070) + .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,201)
         mat(k,760) = mat(k,760) + rxt(k,467)*y(k,124)
         mat(k,1242) = mat(k,1242) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126) &
                      + rxt(k,430)*y(k,200) + .800_r8*rxt(k,431)*y(k,201)
         mat(k,777) = mat(k,777) + rxt(k,473)*y(k,124)
         mat(k,516) = mat(k,516) + rxt(k,476)*y(k,124)
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
         mat(k,472) = -(rxt(k,187)*y(k,228))
         mat(k,1733) = -rxt(k,187)*y(k,91)
         mat(k,2296) = rxt(k,208)*y(k,125)
         mat(k,1931) = rxt(k,208)*y(k,90)
         mat(k,787) = -(rxt(k,265)*y(k,134) + (rxt(k,570) + rxt(k,575)) * y(k,85))
         mat(k,2073) = -rxt(k,265)*y(k,92)
         mat(k,1492) = -(rxt(k,570) + rxt(k,575)) * y(k,92)
         mat(k,1559) = rxt(k,257)*y(k,90)
         mat(k,2316) = rxt(k,257)*y(k,19)
         mat(k,837) = -(rxt(k,236)*y(k,56) + rxt(k,237)*y(k,134) + rxt(k,238)*y(k,228) &
                      + (rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,85))
         mat(k,2139) = -rxt(k,236)*y(k,93)
         mat(k,2076) = -rxt(k,237)*y(k,93)
         mat(k,1770) = -rxt(k,238)*y(k,93)
         mat(k,1493) = -(rxt(k,563) + rxt(k,569) + rxt(k,574)) * y(k,93)
         mat(k,1584) = rxt(k,225)*y(k,90)
         mat(k,954) = rxt(k,230)*y(k,228)
         mat(k,2320) = rxt(k,225)*y(k,59)
         mat(k,1770) = mat(k,1770) + rxt(k,230)*y(k,60)
         mat(k,1195) = -(rxt(k,371)*y(k,228))
         mat(k,1796) = -rxt(k,371)*y(k,94)
         mat(k,591) = .300_r8*rxt(k,416)*y(k,228)
         mat(k,561) = .500_r8*rxt(k,417)*y(k,228)
         mat(k,1894) = rxt(k,370)*y(k,204) + rxt(k,377)*y(k,213)
         mat(k,577) = rxt(k,370)*y(k,124)
         mat(k,1352) = rxt(k,377)*y(k,124)
         mat(k,1796) = mat(k,1796) + .300_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101)
         mat(k,221) = -(rxt(k,402)*y(k,228))
         mat(k,1695) = -rxt(k,402)*y(k,95)
         mat(k,1208) = -(rxt(k,356)*y(k,228))
         mat(k,1797) = -rxt(k,356)*y(k,96)
         mat(k,592) = .700_r8*rxt(k,416)*y(k,228)
         mat(k,562) = .500_r8*rxt(k,417)*y(k,228)
         mat(k,612) = .500_r8*rxt(k,391)*y(k,228)
         mat(k,1895) = .050_r8*rxt(k,414)*y(k,209) + .220_r8*rxt(k,376)*y(k,213) &
                      + .250_r8*rxt(k,433)*y(k,237)
         mat(k,1634) = .050_r8*rxt(k,415)*y(k,209) + .220_r8*rxt(k,375)*y(k,213) &
                      + .250_r8*rxt(k,434)*y(k,237)
         mat(k,545) = .500_r8*rxt(k,360)*y(k,228)
         mat(k,1418) = .220_r8*rxt(k,372)*y(k,213) + .250_r8*rxt(k,430)*y(k,237)
         mat(k,2392) = .230_r8*rxt(k,373)*y(k,213) + .200_r8*rxt(k,361)*y(k,232) &
                      + .100_r8*rxt(k,431)*y(k,237)
         mat(k,1327) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1353) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,200) + .230_r8*rxt(k,373)*y(k,201)
         mat(k,1797) = mat(k,1797) + .700_r8*rxt(k,416)*y(k,100) + .500_r8*rxt(k,417) &
                      *y(k,101) + .500_r8*rxt(k,391)*y(k,109) + .500_r8*rxt(k,360) &
                      *y(k,148)
         mat(k,1216) = .200_r8*rxt(k,361)*y(k,201)
         mat(k,1232) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,200) + .100_r8*rxt(k,431)*y(k,201)
         mat(k,326) = -(rxt(k,403)*y(k,228))
         mat(k,1713) = -rxt(k,403)*y(k,97)
         mat(k,1847) = .870_r8*rxt(k,414)*y(k,209)
         mat(k,1611) = .950_r8*rxt(k,415)*y(k,209)
         mat(k,1410) = rxt(k,410)*y(k,209)
         mat(k,2372) = .750_r8*rxt(k,411)*y(k,209)
         mat(k,1316) = .870_r8*rxt(k,414)*y(k,124) + .950_r8*rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,200) + .750_r8*rxt(k,411)*y(k,201)
         mat(k,133) = -(rxt(k,404)*y(k,228))
         mat(k,1684) = -rxt(k,404)*y(k,98)
         mat(k,710) = .600_r8*rxt(k,427)*y(k,228)
         mat(k,1684) = mat(k,1684) + .600_r8*rxt(k,427)*y(k,103)
         mat(k,876) = -(rxt(k,418)*y(k,126) + rxt(k,425)*y(k,136) + rxt(k,426) &
                      *y(k,228))
         mat(k,1614) = -rxt(k,418)*y(k,99)
         mat(k,2189) = -rxt(k,425)*y(k,99)
         mat(k,1771) = -rxt(k,426)*y(k,99)
         mat(k,589) = -(rxt(k,416)*y(k,228))
         mat(k,1746) = -rxt(k,416)*y(k,100)
         mat(k,1861) = .080_r8*rxt(k,408)*y(k,208)
         mat(k,1289) = .080_r8*rxt(k,408)*y(k,124)
         mat(k,558) = -(rxt(k,417)*y(k,228))
         mat(k,1743) = -rxt(k,417)*y(k,101)
         mat(k,1859) = .080_r8*rxt(k,414)*y(k,209)
         mat(k,1317) = .080_r8*rxt(k,414)*y(k,124)
         mat(k,403) = -(rxt(k,424)*y(k,228))
         mat(k,1723) = -rxt(k,424)*y(k,102)
         mat(k,2285) = rxt(k,421)*y(k,210)
         mat(k,1244) = rxt(k,421)*y(k,90)
         mat(k,711) = -(rxt(k,427)*y(k,228))
         mat(k,1759) = -rxt(k,427)*y(k,103)
         mat(k,2310) = rxt(k,407)*y(k,208) + rxt(k,412)*y(k,209)
         mat(k,1290) = rxt(k,407)*y(k,90)
         mat(k,1319) = rxt(k,412)*y(k,90)
         mat(k,76) = -(rxt(k,549)*y(k,228))
         mat(k,1677) = -rxt(k,549)*y(k,104)
         mat(k,1269) = -(rxt(k,378)*y(k,136) + rxt(k,379)*y(k,228))
         mat(k,2209) = -rxt(k,378)*y(k,105)
         mat(k,1801) = -rxt(k,379)*y(k,105)
         mat(k,881) = .300_r8*rxt(k,425)*y(k,136)
         mat(k,1899) = .360_r8*rxt(k,408)*y(k,208)
         mat(k,1638) = .400_r8*rxt(k,409)*y(k,208)
         mat(k,2209) = mat(k,2209) + .300_r8*rxt(k,425)*y(k,99)
         mat(k,1421) = .390_r8*rxt(k,405)*y(k,208)
         mat(k,2396) = .310_r8*rxt(k,406)*y(k,208)
         mat(k,1297) = .360_r8*rxt(k,408)*y(k,124) + .400_r8*rxt(k,409)*y(k,126) &
                      + .390_r8*rxt(k,405)*y(k,200) + .310_r8*rxt(k,406)*y(k,201)
         mat(k,316) = -(rxt(k,380)*y(k,228))
         mat(k,1711) = -rxt(k,380)*y(k,106)
         mat(k,2278) = rxt(k,374)*y(k,213)
         mat(k,1348) = rxt(k,374)*y(k,90)
         mat(k,517) = -(rxt(k,389)*y(k,228))
         mat(k,1738) = -rxt(k,389)*y(k,107)
         mat(k,1857) = .800_r8*rxt(k,398)*y(k,192)
         mat(k,1037) = .800_r8*rxt(k,398)*y(k,124)
         mat(k,321) = -(rxt(k,390)*y(k,228))
         mat(k,1712) = -rxt(k,390)*y(k,108)
         mat(k,2279) = .800_r8*rxt(k,387)*y(k,217)
         mat(k,697) = .800_r8*rxt(k,387)*y(k,90)
         mat(k,611) = -(rxt(k,391)*y(k,228))
         mat(k,1749) = -rxt(k,391)*y(k,109)
         mat(k,1935) = rxt(k,394)*y(k,215)
         mat(k,1392) = rxt(k,394)*y(k,125)
         mat(k,1018) = -(rxt(k,482)*y(k,126) + rxt(k,483)*y(k,136) + rxt(k,484) &
                      *y(k,228))
         mat(k,1619) = -rxt(k,482)*y(k,110)
         mat(k,2193) = -rxt(k,483)*y(k,110)
         mat(k,1782) = -rxt(k,484)*y(k,110)
         mat(k,1376) = -(rxt(k,392)*y(k,136) + rxt(k,393)*y(k,228))
         mat(k,2214) = -rxt(k,392)*y(k,111)
         mat(k,1806) = -rxt(k,393)*y(k,111)
         mat(k,884) = .200_r8*rxt(k,425)*y(k,136)
         mat(k,1904) = .560_r8*rxt(k,408)*y(k,208)
         mat(k,1643) = .600_r8*rxt(k,409)*y(k,208)
         mat(k,2214) = mat(k,2214) + .200_r8*rxt(k,425)*y(k,99)
         mat(k,1426) = .610_r8*rxt(k,405)*y(k,208)
         mat(k,2401) = .440_r8*rxt(k,406)*y(k,208)
         mat(k,1301) = .560_r8*rxt(k,408)*y(k,124) + .600_r8*rxt(k,409)*y(k,126) &
                      + .610_r8*rxt(k,405)*y(k,200) + .440_r8*rxt(k,406)*y(k,201)
         mat(k,999) = -(rxt(k,190)*y(k,124) + (rxt(k,191) + rxt(k,192) + rxt(k,193) &
                      ) * y(k,125) + rxt(k,194)*y(k,135) + rxt(k,202)*y(k,228) &
                      + rxt(k,588)*y(k,227))
         mat(k,1881) = -rxt(k,190)*y(k,112)
         mat(k,1943) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112)
         mat(k,1536) = -rxt(k,194)*y(k,112)
         mat(k,1781) = -rxt(k,202)*y(k,112)
         mat(k,855) = -rxt(k,588)*y(k,112)
         mat(k,2082) = rxt(k,188)*y(k,219) + rxt(k,585)*y(k,222)
         mat(k,1536) = mat(k,1536) + rxt(k,586)*y(k,222)
         mat(k,866) = 1.100_r8*rxt(k,581)*y(k,220) + .200_r8*rxt(k,579)*y(k,221)
         mat(k,530) = rxt(k,188)*y(k,134)
         mat(k,681) = 1.100_r8*rxt(k,581)*y(k,203)
         mat(k,847) = .200_r8*rxt(k,579)*y(k,203)
         mat(k,506) = rxt(k,585)*y(k,134) + rxt(k,586)*y(k,135)
         mat(k,256) = -((rxt(k,206) + rxt(k,207)) * y(k,224))
         mat(k,2013) = -(rxt(k,206) + rxt(k,207)) * y(k,113)
         mat(k,993) = rxt(k,191)*y(k,125)
         mat(k,1928) = rxt(k,191)*y(k,112)
         mat(k,1929) = rxt(k,209)*y(k,126)
         mat(k,1609) = rxt(k,209)*y(k,125)
         mat(k,427) = -(rxt(k,428)*y(k,228))
         mat(k,1727) = -rxt(k,428)*y(k,115)
         mat(k,2373) = .200_r8*rxt(k,420)*y(k,210)
         mat(k,1245) = .200_r8*rxt(k,420)*y(k,201)
         mat(k,1090) = -(rxt(k,429)*y(k,228))
         mat(k,1788) = -rxt(k,429)*y(k,116)
         mat(k,1887) = rxt(k,422)*y(k,210)
         mat(k,1625) = rxt(k,423)*y(k,210)
         mat(k,1415) = rxt(k,419)*y(k,210)
         mat(k,2385) = .800_r8*rxt(k,420)*y(k,210)
         mat(k,1249) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) + rxt(k,419)*y(k,200) &
                      + .800_r8*rxt(k,420)*y(k,201)
         mat(k,98) = -(rxt(k,519)*y(k,228))
         mat(k,1681) = -rxt(k,519)*y(k,120)
         mat(k,1914) = -(rxt(k,190)*y(k,112) + rxt(k,199)*y(k,126) + rxt(k,203) &
                      *y(k,90) + rxt(k,204)*y(k,136) + rxt(k,205)*y(k,134) + rxt(k,226) &
                      *y(k,59) + rxt(k,258)*y(k,19) + rxt(k,301)*y(k,201) + rxt(k,309) &
                      *y(k,207) + rxt(k,322)*y(k,197) + rxt(k,333)*y(k,200) + rxt(k,337) &
                      *y(k,206) + rxt(k,350)*y(k,198) + rxt(k,359)*y(k,231) + rxt(k,363) &
                      *y(k,232) + (rxt(k,369) + rxt(k,370)) * y(k,204) + (rxt(k,376) &
                      + rxt(k,377)) * y(k,213) + rxt(k,385)*y(k,215) + rxt(k,388) &
                      *y(k,217) + (rxt(k,398) + rxt(k,399)) * y(k,192) + rxt(k,408) &
                      *y(k,208) + rxt(k,414)*y(k,209) + rxt(k,422)*y(k,210) + rxt(k,433) &
                      *y(k,237) + rxt(k,437)*y(k,191) + rxt(k,440)*y(k,194) + rxt(k,445) &
                      *y(k,196) + rxt(k,447)*y(k,199) + rxt(k,451)*y(k,202) + rxt(k,454) &
                      *y(k,214) + rxt(k,457)*y(k,216) + rxt(k,460)*y(k,230) + rxt(k,467) &
                      *y(k,235) + rxt(k,473)*y(k,238) + rxt(k,476)*y(k,240) + rxt(k,487) &
                      *y(k,223) + rxt(k,492)*y(k,233) + rxt(k,497)*y(k,234) + rxt(k,590) &
                      *y(k,227))
         mat(k,1004) = -rxt(k,190)*y(k,124)
         mat(k,1654) = -rxt(k,199)*y(k,124)
         mat(k,2358) = -rxt(k,203)*y(k,124)
         mat(k,2225) = -rxt(k,204)*y(k,124)
         mat(k,2094) = -rxt(k,205)*y(k,124)
         mat(k,1595) = -rxt(k,226)*y(k,124)
         mat(k,1569) = -rxt(k,258)*y(k,124)
         mat(k,2410) = -rxt(k,301)*y(k,124)
         mat(k,443) = -rxt(k,309)*y(k,124)
         mat(k,905) = -rxt(k,322)*y(k,124)
         mat(k,1433) = -rxt(k,333)*y(k,124)
         mat(k,801) = -rxt(k,337)*y(k,124)
         mat(k,940) = -rxt(k,350)*y(k,124)
         mat(k,819) = -rxt(k,359)*y(k,124)
         mat(k,1221) = -rxt(k,363)*y(k,124)
         mat(k,578) = -(rxt(k,369) + rxt(k,370)) * y(k,124)
         mat(k,1362) = -(rxt(k,376) + rxt(k,377)) * y(k,124)
         mat(k,1402) = -rxt(k,385)*y(k,124)
         mat(k,702) = -rxt(k,388)*y(k,124)
         mat(k,1048) = -(rxt(k,398) + rxt(k,399)) * y(k,124)
         mat(k,1307) = -rxt(k,408)*y(k,124)
         mat(k,1340) = -rxt(k,414)*y(k,124)
         mat(k,1261) = -rxt(k,422)*y(k,124)
         mat(k,1239) = -rxt(k,433)*y(k,124)
         mat(k,526) = -rxt(k,437)*y(k,124)
         mat(k,499) = -rxt(k,440)*y(k,124)
         mat(k,438) = -rxt(k,445)*y(k,124)
         mat(k,637) = -rxt(k,447)*y(k,124)
         mat(k,782) = -rxt(k,451)*y(k,124)
         mat(k,743) = -rxt(k,454)*y(k,124)
         mat(k,915) = -rxt(k,457)*y(k,124)
         mat(k,451) = -rxt(k,460)*y(k,124)
         mat(k,758) = -rxt(k,467)*y(k,124)
         mat(k,775) = -rxt(k,473)*y(k,124)
         mat(k,514) = -rxt(k,476)*y(k,124)
         mat(k,1112) = -rxt(k,487)*y(k,124)
         mat(k,1183) = -rxt(k,492)*y(k,124)
         mat(k,1066) = -rxt(k,497)*y(k,124)
         mat(k,857) = -rxt(k,590)*y(k,124)
         mat(k,1004) = mat(k,1004) + 2.000_r8*rxt(k,192)*y(k,125) + rxt(k,194) &
                      *y(k,135) + rxt(k,202)*y(k,228)
         mat(k,258) = 2.000_r8*rxt(k,206)*y(k,224)
         mat(k,1959) = 2.000_r8*rxt(k,192)*y(k,112) + rxt(k,195)*y(k,134) + rxt(k,512) &
                      *y(k,152)
         mat(k,2094) = mat(k,2094) + rxt(k,195)*y(k,125)
         mat(k,1546) = rxt(k,194)*y(k,112) + rxt(k,189)*y(k,219)
         mat(k,1482) = rxt(k,512)*y(k,125)
         mat(k,532) = rxt(k,189)*y(k,135)
         mat(k,2028) = 2.000_r8*rxt(k,206)*y(k,113)
         mat(k,1819) = rxt(k,202)*y(k,112)
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
         mat(k,1960) = -((rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112) + (rxt(k,195) &
                      + rxt(k,197)) * y(k,134) + rxt(k,196)*y(k,136) + rxt(k,208) &
                      *y(k,90) + rxt(k,209)*y(k,126) + rxt(k,210)*y(k,228) + rxt(k,228) &
                      *y(k,59) + rxt(k,259)*y(k,19) + rxt(k,344)*y(k,200) + rxt(k,394) &
                      *y(k,215) + rxt(k,452)*y(k,202) + rxt(k,455)*y(k,214) + rxt(k,458) &
                      *y(k,216) + rxt(k,462)*y(k,143) + rxt(k,465)*y(k,191) + rxt(k,512) &
                      *y(k,152))
         mat(k,1005) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,125)
         mat(k,2095) = -(rxt(k,195) + rxt(k,197)) * y(k,125)
         mat(k,2226) = -rxt(k,196)*y(k,125)
         mat(k,2359) = -rxt(k,208)*y(k,125)
         mat(k,1655) = -rxt(k,209)*y(k,125)
         mat(k,1820) = -rxt(k,210)*y(k,125)
         mat(k,1596) = -rxt(k,228)*y(k,125)
         mat(k,1570) = -rxt(k,259)*y(k,125)
         mat(k,1434) = -rxt(k,344)*y(k,125)
         mat(k,1403) = -rxt(k,394)*y(k,125)
         mat(k,783) = -rxt(k,452)*y(k,125)
         mat(k,744) = -rxt(k,455)*y(k,125)
         mat(k,916) = -rxt(k,458)*y(k,125)
         mat(k,470) = -rxt(k,462)*y(k,125)
         mat(k,527) = -rxt(k,465)*y(k,125)
         mat(k,1483) = -rxt(k,512)*y(k,125)
         mat(k,649) = rxt(k,396)*y(k,228)
         mat(k,361) = rxt(k,367)*y(k,126)
         mat(k,1570) = mat(k,1570) + rxt(k,258)*y(k,124)
         mat(k,1596) = mat(k,1596) + rxt(k,226)*y(k,124)
         mat(k,2359) = mat(k,2359) + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126)
         mat(k,476) = rxt(k,187)*y(k,228)
         mat(k,594) = .700_r8*rxt(k,416)*y(k,228)
         mat(k,1915) = rxt(k,258)*y(k,19) + rxt(k,226)*y(k,59) + rxt(k,203)*y(k,90) &
                      + 2.000_r8*rxt(k,199)*y(k,126) + rxt(k,205)*y(k,134) &
                      + rxt(k,204)*y(k,136) + rxt(k,437)*y(k,191) + rxt(k,398) &
                      *y(k,192) + rxt(k,440)*y(k,194) + rxt(k,445)*y(k,196) &
                      + rxt(k,322)*y(k,197) + rxt(k,350)*y(k,198) + rxt(k,447) &
                      *y(k,199) + rxt(k,333)*y(k,200) + rxt(k,301)*y(k,201) &
                      + rxt(k,451)*y(k,202) + rxt(k,369)*y(k,204) + rxt(k,337) &
                      *y(k,206) + rxt(k,309)*y(k,207) + .920_r8*rxt(k,408)*y(k,208) &
                      + .920_r8*rxt(k,414)*y(k,209) + rxt(k,422)*y(k,210) + rxt(k,376) &
                      *y(k,213) + rxt(k,454)*y(k,214) + rxt(k,385)*y(k,215) &
                      + rxt(k,457)*y(k,216) + rxt(k,388)*y(k,217) &
                      + 1.600_r8*rxt(k,487)*y(k,223) + rxt(k,460)*y(k,230) &
                      + rxt(k,359)*y(k,231) + rxt(k,363)*y(k,232) + .900_r8*rxt(k,492) &
                      *y(k,233) + .800_r8*rxt(k,497)*y(k,234) + rxt(k,467)*y(k,235) &
                      + rxt(k,433)*y(k,237) + rxt(k,473)*y(k,238) + rxt(k,476) &
                      *y(k,240)
         mat(k,1655) = mat(k,1655) + rxt(k,367)*y(k,16) + rxt(k,198)*y(k,90) &
                      + 2.000_r8*rxt(k,199)*y(k,124) + rxt(k,200)*y(k,134) &
                      + rxt(k,409)*y(k,208) + rxt(k,415)*y(k,209) + rxt(k,423) &
                      *y(k,210) + rxt(k,375)*y(k,213) + rxt(k,386)*y(k,215) &
                      + 2.000_r8*rxt(k,488)*y(k,223) + rxt(k,201)*y(k,228) &
                      + rxt(k,434)*y(k,237)
         mat(k,896) = rxt(k,357)*y(k,228)
         mat(k,2095) = mat(k,2095) + rxt(k,205)*y(k,124) + rxt(k,200)*y(k,126)
         mat(k,2226) = mat(k,2226) + rxt(k,204)*y(k,124)
         mat(k,630) = rxt(k,494)*y(k,228)
         mat(k,527) = mat(k,527) + rxt(k,437)*y(k,124)
         mat(k,1049) = rxt(k,398)*y(k,124)
         mat(k,500) = rxt(k,440)*y(k,124)
         mat(k,439) = rxt(k,445)*y(k,124)
         mat(k,906) = rxt(k,322)*y(k,124)
         mat(k,941) = rxt(k,350)*y(k,124)
         mat(k,638) = rxt(k,447)*y(k,124)
         mat(k,1434) = mat(k,1434) + rxt(k,333)*y(k,124)
         mat(k,2411) = rxt(k,301)*y(k,124) + .500_r8*rxt(k,485)*y(k,223)
         mat(k,783) = mat(k,783) + rxt(k,451)*y(k,124)
         mat(k,579) = rxt(k,369)*y(k,124)
         mat(k,802) = rxt(k,337)*y(k,124)
         mat(k,444) = rxt(k,309)*y(k,124)
         mat(k,1308) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126)
         mat(k,1341) = .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126)
         mat(k,1262) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126)
         mat(k,1363) = rxt(k,376)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,744) = mat(k,744) + rxt(k,454)*y(k,124)
         mat(k,1403) = mat(k,1403) + rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126)
         mat(k,916) = mat(k,916) + rxt(k,457)*y(k,124)
         mat(k,703) = rxt(k,388)*y(k,124)
         mat(k,1113) = 1.600_r8*rxt(k,487)*y(k,124) + 2.000_r8*rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,201)
         mat(k,1820) = mat(k,1820) + rxt(k,396)*y(k,1) + rxt(k,187)*y(k,91) &
                      + .700_r8*rxt(k,416)*y(k,100) + rxt(k,201)*y(k,126) + rxt(k,357) &
                      *y(k,127) + rxt(k,494)*y(k,178)
         mat(k,452) = rxt(k,460)*y(k,124)
         mat(k,820) = rxt(k,359)*y(k,124)
         mat(k,1222) = rxt(k,363)*y(k,124)
         mat(k,1184) = .900_r8*rxt(k,492)*y(k,124)
         mat(k,1067) = .800_r8*rxt(k,497)*y(k,124)
         mat(k,759) = rxt(k,467)*y(k,124)
         mat(k,1240) = rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126)
         mat(k,776) = rxt(k,473)*y(k,124)
         mat(k,515) = rxt(k,476)*y(k,124)
         mat(k,1652) = -(rxt(k,198)*y(k,90) + rxt(k,199)*y(k,124) + rxt(k,200) &
                      *y(k,134) + rxt(k,201)*y(k,228) + rxt(k,209)*y(k,125) + rxt(k,295) &
                      *y(k,42) + rxt(k,327)*y(k,45) + rxt(k,346)*y(k,29) + rxt(k,353) &
                      *y(k,49) + rxt(k,367)*y(k,16) + rxt(k,375)*y(k,213) + rxt(k,386) &
                      *y(k,215) + rxt(k,409)*y(k,208) + rxt(k,415)*y(k,209) + rxt(k,418) &
                      *y(k,99) + rxt(k,423)*y(k,210) + rxt(k,434)*y(k,237) + rxt(k,479) &
                      *y(k,6) + rxt(k,482)*y(k,110) + rxt(k,488)*y(k,223) + rxt(k,499) &
                      *y(k,180) + rxt(k,502)*y(k,67))
         mat(k,2356) = -rxt(k,198)*y(k,126)
         mat(k,1912) = -rxt(k,199)*y(k,126)
         mat(k,2092) = -rxt(k,200)*y(k,126)
         mat(k,1817) = -rxt(k,201)*y(k,126)
         mat(k,1957) = -rxt(k,209)*y(k,126)
         mat(k,1983) = -rxt(k,295)*y(k,126)
         mat(k,1154) = -rxt(k,327)*y(k,126)
         mat(k,1141) = -rxt(k,346)*y(k,126)
         mat(k,1284) = -rxt(k,353)*y(k,126)
         mat(k,359) = -rxt(k,367)*y(k,126)
         mat(k,1360) = -rxt(k,375)*y(k,126)
         mat(k,1400) = -rxt(k,386)*y(k,126)
         mat(k,1305) = -rxt(k,409)*y(k,126)
         mat(k,1338) = -rxt(k,415)*y(k,126)
         mat(k,886) = -rxt(k,418)*y(k,126)
         mat(k,1259) = -rxt(k,423)*y(k,126)
         mat(k,1237) = -rxt(k,434)*y(k,126)
         mat(k,987) = -rxt(k,479)*y(k,126)
         mat(k,1031) = -rxt(k,482)*y(k,126)
         mat(k,1110) = -rxt(k,488)*y(k,126)
         mat(k,1077) = -rxt(k,499)*y(k,126)
         mat(k,312) = -rxt(k,502)*y(k,126)
         mat(k,570) = rxt(k,260)*y(k,134)
         mat(k,2159) = rxt(k,227)*y(k,60)
         mat(k,958) = rxt(k,227)*y(k,56) + rxt(k,229)*y(k,134) + rxt(k,230)*y(k,228)
         mat(k,927) = rxt(k,274)*y(k,89)
         mat(k,2049) = rxt(k,274)*y(k,73) + rxt(k,211)*y(k,228)
         mat(k,615) = .500_r8*rxt(k,391)*y(k,228)
         mat(k,1957) = mat(k,1957) + rxt(k,197)*y(k,134) + rxt(k,196)*y(k,136)
         mat(k,2092) = mat(k,2092) + rxt(k,260)*y(k,20) + rxt(k,229)*y(k,60) &
                      + rxt(k,197)*y(k,125)
         mat(k,2223) = rxt(k,196)*y(k,125)
         mat(k,537) = rxt(k,342)*y(k,228)
         mat(k,1817) = mat(k,1817) + rxt(k,230)*y(k,60) + rxt(k,211)*y(k,89) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,342)*y(k,141)
         mat(k,892) = -(rxt(k,357)*y(k,228))
         mat(k,1772) = -rxt(k,357)*y(k,127)
         mat(k,1131) = rxt(k,346)*y(k,126)
         mat(k,559) = .500_r8*rxt(k,417)*y(k,228)
         mat(k,405) = rxt(k,424)*y(k,228)
         mat(k,428) = rxt(k,428)*y(k,228)
         mat(k,1087) = rxt(k,429)*y(k,228)
         mat(k,1615) = rxt(k,346)*y(k,29)
         mat(k,1772) = mat(k,1772) + .500_r8*rxt(k,417)*y(k,101) + rxt(k,424)*y(k,102) &
                      + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116)
         mat(k,391) = -(rxt(k,489)*y(k,228))
         mat(k,1721) = -rxt(k,489)*y(k,128)
         mat(k,2283) = rxt(k,486)*y(k,223)
         mat(k,1102) = rxt(k,486)*y(k,90)
         mat(k,2099) = -(rxt(k,167)*y(k,136) + 4._r8*rxt(k,168)*y(k,134) + rxt(k,169) &
                      *y(k,135) + rxt(k,170)*y(k,77) + rxt(k,171)*y(k,79) + rxt(k,176) &
                      *y(k,90) + rxt(k,182)*y(k,228) + (rxt(k,195) + rxt(k,197) &
                      ) * y(k,125) + rxt(k,200)*y(k,126) + rxt(k,205)*y(k,124) &
                      + rxt(k,229)*y(k,60) + rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85) &
                      + rxt(k,237)*y(k,93) + rxt(k,260)*y(k,20) + rxt(k,261)*y(k,19) &
                      + rxt(k,263)*y(k,81) + rxt(k,265)*y(k,92) + rxt(k,296)*y(k,42) &
                      + rxt(k,504)*y(k,139) + (rxt(k,583) + rxt(k,584)) * y(k,220) &
                      + rxt(k,585)*y(k,222))
         mat(k,2230) = -rxt(k,167)*y(k,134)
         mat(k,1550) = -rxt(k,169)*y(k,134)
         mat(k,1467) = -rxt(k,170)*y(k,134)
         mat(k,585) = -rxt(k,171)*y(k,134)
         mat(k,2363) = -rxt(k,176)*y(k,134)
         mat(k,1824) = -rxt(k,182)*y(k,134)
         mat(k,1964) = -(rxt(k,195) + rxt(k,197)) * y(k,134)
         mat(k,1659) = -rxt(k,200)*y(k,134)
         mat(k,1919) = -rxt(k,205)*y(k,134)
         mat(k,962) = -rxt(k,229)*y(k,134)
         mat(k,1600) = -rxt(k,231)*y(k,134)
         mat(k,1504) = -rxt(k,234)*y(k,134)
         mat(k,841) = -rxt(k,237)*y(k,134)
         mat(k,573) = -rxt(k,260)*y(k,134)
         mat(k,1574) = -rxt(k,261)*y(k,134)
         mat(k,833) = -rxt(k,263)*y(k,134)
         mat(k,792) = -rxt(k,265)*y(k,134)
         mat(k,1990) = -rxt(k,296)*y(k,134)
         mat(k,369) = -rxt(k,504)*y(k,134)
         mat(k,685) = -(rxt(k,583) + rxt(k,584)) * y(k,134)
         mat(k,508) = -rxt(k,585)*y(k,134)
         mat(k,2120) = rxt(k,174)*y(k,90)
         mat(k,2363) = mat(k,2363) + rxt(k,174)*y(k,76)
         mat(k,1007) = rxt(k,190)*y(k,124) + rxt(k,191)*y(k,125) + rxt(k,194)*y(k,135) &
                      + rxt(k,588)*y(k,227)
         mat(k,1919) = mat(k,1919) + rxt(k,190)*y(k,112)
         mat(k,1964) = mat(k,1964) + rxt(k,191)*y(k,112)
         mat(k,1550) = mat(k,1550) + rxt(k,194)*y(k,112) + rxt(k,506)*y(k,150) &
                      + rxt(k,513)*y(k,152) + rxt(k,587)*y(k,222) + (rxt(k,156) &
                       +rxt(k,157))*y(k,224) + rxt(k,593)*y(k,229)
         mat(k,722) = rxt(k,506)*y(k,135)
         mat(k,1484) = rxt(k,513)*y(k,135)
         mat(k,872) = rxt(k,579)*y(k,221) + 1.150_r8*rxt(k,580)*y(k,227)
         mat(k,851) = rxt(k,579)*y(k,203)
         mat(k,508) = mat(k,508) + rxt(k,587)*y(k,135)
         mat(k,2033) = (rxt(k,156)+rxt(k,157))*y(k,135)
         mat(k,859) = rxt(k,588)*y(k,112) + 1.150_r8*rxt(k,580)*y(k,203)
         mat(k,1824) = mat(k,1824) + 2.000_r8*rxt(k,184)*y(k,228)
         mat(k,812) = rxt(k,593)*y(k,135)
         mat(k,1542) = -(rxt(k,156)*y(k,224) + rxt(k,161)*y(k,225) + rxt(k,169) &
                      *y(k,134) + rxt(k,175)*y(k,76) + rxt(k,189)*y(k,219) + rxt(k,194) &
                      *y(k,112) + rxt(k,339)*y(k,205) + rxt(k,506)*y(k,150) + rxt(k,513) &
                      *y(k,152) + rxt(k,582)*y(k,220) + (rxt(k,586) + rxt(k,587) &
                      ) * y(k,222) + rxt(k,593)*y(k,229))
         mat(k,2023) = -rxt(k,156)*y(k,135)
         mat(k,178) = -rxt(k,161)*y(k,135)
         mat(k,2089) = -rxt(k,169)*y(k,135)
         mat(k,2110) = -rxt(k,175)*y(k,135)
         mat(k,531) = -rxt(k,189)*y(k,135)
         mat(k,1002) = -rxt(k,194)*y(k,135)
         mat(k,462) = -rxt(k,339)*y(k,135)
         mat(k,720) = -rxt(k,506)*y(k,135)
         mat(k,1478) = -rxt(k,513)*y(k,135)
         mat(k,682) = -rxt(k,582)*y(k,135)
         mat(k,507) = -(rxt(k,586) + rxt(k,587)) * y(k,135)
         mat(k,811) = -rxt(k,593)*y(k,135)
         mat(k,1512) = rxt(k,251)*y(k,90) + rxt(k,252)*y(k,136)
         mat(k,1564) = 2.000_r8*rxt(k,253)*y(k,19) + (rxt(k,255)+rxt(k,256))*y(k,59) &
                      + rxt(k,257)*y(k,90) + rxt(k,261)*y(k,134)
         mat(k,2156) = rxt(k,218)*y(k,90) + rxt(k,220)*y(k,136)
         mat(k,1590) = (rxt(k,255)+rxt(k,256))*y(k,19) + (2.000_r8*rxt(k,222) &
                       +2.000_r8*rxt(k,223))*y(k,59) + rxt(k,225)*y(k,90) + rxt(k,231) &
                      *y(k,134) + rxt(k,233)*y(k,228)
         mat(k,2110) = mat(k,2110) + rxt(k,172)*y(k,90) + rxt(k,178)*y(k,136)
         mat(k,2353) = rxt(k,251)*y(k,17) + rxt(k,257)*y(k,19) + rxt(k,218)*y(k,56) &
                      + rxt(k,225)*y(k,59) + rxt(k,172)*y(k,76) + 2.000_r8*rxt(k,186) &
                      *y(k,90) + rxt(k,198)*y(k,126) + rxt(k,176)*y(k,134) &
                      + 2.000_r8*rxt(k,177)*y(k,136) + rxt(k,321)*y(k,197) &
                      + rxt(k,349)*y(k,198) + rxt(k,300)*y(k,201) + rxt(k,181) &
                      *y(k,228) + rxt(k,358)*y(k,231)
         mat(k,473) = rxt(k,187)*y(k,228)
         mat(k,1002) = mat(k,1002) + rxt(k,193)*y(k,125)
         mat(k,257) = rxt(k,207)*y(k,224)
         mat(k,1909) = rxt(k,204)*y(k,136) + rxt(k,590)*y(k,227)
         mat(k,1954) = rxt(k,193)*y(k,112) + rxt(k,195)*y(k,134) + rxt(k,196)*y(k,136)
         mat(k,1649) = rxt(k,198)*y(k,90) + rxt(k,200)*y(k,134)
         mat(k,2089) = mat(k,2089) + rxt(k,261)*y(k,19) + rxt(k,231)*y(k,59) &
                      + rxt(k,176)*y(k,90) + rxt(k,195)*y(k,125) + rxt(k,200)*y(k,126) &
                      + 2.000_r8*rxt(k,168)*y(k,134) + 2.000_r8*rxt(k,167)*y(k,136) &
                      + rxt(k,160)*y(k,225) + rxt(k,182)*y(k,228)
         mat(k,1542) = mat(k,1542) + 2.000_r8*rxt(k,161)*y(k,225)
         mat(k,2220) = rxt(k,252)*y(k,17) + rxt(k,220)*y(k,56) + rxt(k,178)*y(k,76) &
                      + 2.000_r8*rxt(k,177)*y(k,90) + rxt(k,204)*y(k,124) + rxt(k,196) &
                      *y(k,125) + 2.000_r8*rxt(k,167)*y(k,134) + rxt(k,508)*y(k,150) &
                      + rxt(k,514)*y(k,152) + 2.000_r8*rxt(k,158)*y(k,224) &
                      + rxt(k,183)*y(k,228)
         mat(k,720) = mat(k,720) + rxt(k,508)*y(k,136)
         mat(k,1478) = mat(k,1478) + rxt(k,514)*y(k,136)
         mat(k,903) = rxt(k,321)*y(k,90)
         mat(k,938) = rxt(k,349)*y(k,90)
         mat(k,2405) = rxt(k,300)*y(k,90)
         mat(k,2023) = mat(k,2023) + rxt(k,207)*y(k,113) + 2.000_r8*rxt(k,158) &
                      *y(k,136)
         mat(k,178) = mat(k,178) + rxt(k,160)*y(k,134) + 2.000_r8*rxt(k,161)*y(k,135)
         mat(k,856) = rxt(k,590)*y(k,124)
         mat(k,1814) = rxt(k,233)*y(k,59) + rxt(k,181)*y(k,90) + rxt(k,187)*y(k,91) &
                      + rxt(k,182)*y(k,134) + rxt(k,183)*y(k,136)
         mat(k,817) = rxt(k,358)*y(k,90)
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
         mat(k,2233) = -(rxt(k,158)*y(k,224) + rxt(k,167)*y(k,134) + rxt(k,177) &
                      *y(k,90) + rxt(k,178)*y(k,76) + rxt(k,183)*y(k,228) + rxt(k,196) &
                      *y(k,125) + rxt(k,204)*y(k,124) + rxt(k,220)*y(k,56) + rxt(k,252) &
                      *y(k,17) + rxt(k,318)*y(k,25) + rxt(k,347)*y(k,29) + rxt(k,378) &
                      *y(k,105) + rxt(k,392)*y(k,111) + rxt(k,425)*y(k,99) + rxt(k,463) &
                      *y(k,143) + rxt(k,480)*y(k,6) + rxt(k,483)*y(k,110) + rxt(k,508) &
                      *y(k,150) + rxt(k,514)*y(k,152))
         mat(k,2036) = -rxt(k,158)*y(k,136)
         mat(k,2102) = -rxt(k,167)*y(k,136)
         mat(k,2366) = -rxt(k,177)*y(k,136)
         mat(k,2123) = -rxt(k,178)*y(k,136)
         mat(k,1827) = -rxt(k,183)*y(k,136)
         mat(k,1967) = -rxt(k,196)*y(k,136)
         mat(k,1922) = -rxt(k,204)*y(k,136)
         mat(k,2169) = -rxt(k,220)*y(k,136)
         mat(k,1521) = -rxt(k,252)*y(k,136)
         mat(k,556) = -rxt(k,318)*y(k,136)
         mat(k,1146) = -rxt(k,347)*y(k,136)
         mat(k,1277) = -rxt(k,378)*y(k,136)
         mat(k,1387) = -rxt(k,392)*y(k,136)
         mat(k,889) = -rxt(k,425)*y(k,136)
         mat(k,471) = -rxt(k,463)*y(k,136)
         mat(k,990) = -rxt(k,480)*y(k,136)
         mat(k,1034) = -rxt(k,483)*y(k,136)
         mat(k,724) = -rxt(k,508)*y(k,136)
         mat(k,1487) = -rxt(k,514)*y(k,136)
         mat(k,2366) = mat(k,2366) + .150_r8*rxt(k,332)*y(k,200) + .150_r8*rxt(k,383) &
                      *y(k,215)
         mat(k,2102) = mat(k,2102) + rxt(k,169)*y(k,135)
         mat(k,1553) = rxt(k,169)*y(k,134)
         mat(k,1437) = .150_r8*rxt(k,332)*y(k,90)
         mat(k,1406) = .150_r8*rxt(k,383)*y(k,90)
         mat(k,334) = -(rxt(k,515)*y(k,152))
         mat(k,1473) = -rxt(k,515)*y(k,138)
         mat(k,1557) = rxt(k,254)*y(k,59)
         mat(k,1583) = rxt(k,254)*y(k,19) + 2.000_r8*rxt(k,224)*y(k,59)
         mat(k,363) = -(rxt(k,504)*y(k,134) + rxt(k,505)*y(k,228))
         mat(k,2066) = -rxt(k,504)*y(k,139)
         mat(k,1718) = -rxt(k,505)*y(k,139)
         mat(k,1190) = rxt(k,371)*y(k,228)
         mat(k,1844) = .100_r8*rxt(k,492)*y(k,233)
         mat(k,1697) = rxt(k,371)*y(k,94)
         mat(k,1171) = .100_r8*rxt(k,492)*y(k,124)
         mat(k,534) = -(rxt(k,342)*y(k,228))
         mat(k,1740) = -rxt(k,342)*y(k,141)
         mat(k,1933) = rxt(k,344)*y(k,200)
         mat(k,1411) = rxt(k,344)*y(k,125)
         mat(k,1927) = rxt(k,465)*y(k,191)
         mat(k,522) = rxt(k,465)*y(k,125)
         mat(k,468) = -(rxt(k,462)*y(k,125) + rxt(k,463)*y(k,136))
         mat(k,1930) = -rxt(k,462)*y(k,143)
         mat(k,2183) = -rxt(k,463)*y(k,143)
         mat(k,196) = .070_r8*rxt(k,449)*y(k,228)
         mat(k,1854) = rxt(k,447)*y(k,199)
         mat(k,172) = .060_r8*rxt(k,461)*y(k,228)
         mat(k,217) = .070_r8*rxt(k,477)*y(k,228)
         mat(k,634) = rxt(k,447)*y(k,124)
         mat(k,1732) = .070_r8*rxt(k,449)*y(k,66) + .060_r8*rxt(k,461)*y(k,144) &
                      + .070_r8*rxt(k,477)*y(k,187)
         mat(k,170) = -(rxt(k,461)*y(k,228))
         mat(k,1687) = -rxt(k,461)*y(k,144)
         mat(k,162) = .530_r8*rxt(k,438)*y(k,228)
         mat(k,1687) = mat(k,1687) + .530_r8*rxt(k,438)*y(k,7)
         mat(k,339) = -(rxt(k,464)*y(k,228))
         mat(k,1714) = -rxt(k,464)*y(k,145)
         mat(k,2280) = rxt(k,459)*y(k,230)
         mat(k,447) = rxt(k,459)*y(k,90)
         mat(k,542) = -(rxt(k,360)*y(k,228))
         mat(k,1741) = -rxt(k,360)*y(k,148)
         mat(k,2301) = rxt(k,358)*y(k,231)
         mat(k,813) = rxt(k,358)*y(k,90)
         mat(k,397) = -(rxt(k,364)*y(k,228))
         mat(k,1722) = -rxt(k,364)*y(k,149)
         mat(k,2284) = .850_r8*rxt(k,362)*y(k,232)
         mat(k,1214) = .850_r8*rxt(k,362)*y(k,90)
         mat(k,718) = -(rxt(k,506)*y(k,135) + rxt(k,508)*y(k,136) + rxt(k,511) &
                      *y(k,228))
         mat(k,1530) = -rxt(k,506)*y(k,150)
         mat(k,2187) = -rxt(k,508)*y(k,150)
         mat(k,1760) = -rxt(k,511)*y(k,150)
         mat(k,1476) = -(rxt(k,509)*y(k,19) + rxt(k,510)*y(k,59) + rxt(k,512)*y(k,125) &
                      + rxt(k,513)*y(k,135) + rxt(k,514)*y(k,136) + rxt(k,515) &
                      *y(k,138) + rxt(k,516)*y(k,228))
         mat(k,1561) = -rxt(k,509)*y(k,152)
         mat(k,1587) = -rxt(k,510)*y(k,152)
         mat(k,1951) = -rxt(k,512)*y(k,152)
         mat(k,1540) = -rxt(k,513)*y(k,152)
         mat(k,2218) = -rxt(k,514)*y(k,152)
         mat(k,336) = -rxt(k,515)*y(k,152)
         mat(k,1811) = -rxt(k,516)*y(k,152)
         mat(k,2086) = rxt(k,504)*y(k,139)
         mat(k,1540) = mat(k,1540) + rxt(k,506)*y(k,150)
         mat(k,2218) = mat(k,2218) + rxt(k,508)*y(k,150)
         mat(k,367) = rxt(k,504)*y(k,134)
         mat(k,719) = rxt(k,506)*y(k,135) + rxt(k,508)*y(k,136) + rxt(k,511)*y(k,228)
         mat(k,1811) = mat(k,1811) + rxt(k,511)*y(k,150)
         mat(k,947) = -(rxt(k,507)*y(k,228))
         mat(k,1778) = -rxt(k,507)*y(k,153)
         mat(k,1560) = rxt(k,509)*y(k,152)
         mat(k,1585) = rxt(k,510)*y(k,152)
         mat(k,311) = rxt(k,502)*y(k,126) + (rxt(k,503)+.500_r8*rxt(k,517))*y(k,228)
         mat(k,1941) = rxt(k,512)*y(k,152)
         mat(k,1617) = rxt(k,502)*y(k,67)
         mat(k,1535) = rxt(k,513)*y(k,152)
         mat(k,2191) = rxt(k,514)*y(k,152)
         mat(k,335) = rxt(k,515)*y(k,152)
         mat(k,365) = rxt(k,505)*y(k,228)
         mat(k,1475) = rxt(k,509)*y(k,19) + rxt(k,510)*y(k,59) + rxt(k,512)*y(k,125) &
                      + rxt(k,513)*y(k,135) + rxt(k,514)*y(k,136) + rxt(k,515) &
                      *y(k,138) + rxt(k,516)*y(k,228)
         mat(k,1778) = mat(k,1778) + (rxt(k,503)+.500_r8*rxt(k,517))*y(k,67) &
                      + rxt(k,505)*y(k,139) + rxt(k,516)*y(k,152)
         mat(k,261) = -(rxt(k,518)*y(k,241))
         mat(k,2424) = -rxt(k,518)*y(k,154)
         mat(k,946) = rxt(k,507)*y(k,228)
         mat(k,1703) = rxt(k,507)*y(k,153)
         mat(k,965) = .2202005_r8*rxt(k,537)*y(k,136)
         mat(k,2263) = .2202005_r8*rxt(k,535)*y(k,193) + .0023005_r8*rxt(k,540) &
                      *y(k,195) + .0031005_r8*rxt(k,543)*y(k,211) &
                      + .2381005_r8*rxt(k,547)*y(k,212) + .0508005_r8*rxt(k,551) &
                      *y(k,218) + .1364005_r8*rxt(k,557)*y(k,236) &
                      + .1677005_r8*rxt(k,560)*y(k,239)
         mat(k,1009) = .0508005_r8*rxt(k,553)*y(k,136)
         mat(k,1832) = .1279005_r8*rxt(k,536)*y(k,193) + .0097005_r8*rxt(k,541) &
                      *y(k,195) + .0003005_r8*rxt(k,544)*y(k,211) &
                      + .1056005_r8*rxt(k,548)*y(k,212) + .0245005_r8*rxt(k,552) &
                      *y(k,218) + .0154005_r8*rxt(k,558)*y(k,236) &
                      + .0063005_r8*rxt(k,561)*y(k,239)
         mat(k,2174) = .2202005_r8*rxt(k,537)*y(k,6) + .0508005_r8*rxt(k,553)*y(k,110)
         mat(k,45) = .5931005_r8*rxt(k,555)*y(k,228)
         mat(k,51) = .2202005_r8*rxt(k,535)*y(k,90) + .1279005_r8*rxt(k,536)*y(k,124)
         mat(k,57) = .0023005_r8*rxt(k,540)*y(k,90) + .0097005_r8*rxt(k,541)*y(k,124)
         mat(k,63) = .0031005_r8*rxt(k,543)*y(k,90) + .0003005_r8*rxt(k,544)*y(k,124)
         mat(k,69) = .2381005_r8*rxt(k,547)*y(k,90) + .1056005_r8*rxt(k,548)*y(k,124)
         mat(k,77) = .0508005_r8*rxt(k,551)*y(k,90) + .0245005_r8*rxt(k,552)*y(k,124)
         mat(k,1667) = .5931005_r8*rxt(k,555)*y(k,175)
         mat(k,83) = .1364005_r8*rxt(k,557)*y(k,90) + .0154005_r8*rxt(k,558)*y(k,124)
         mat(k,89) = .1677005_r8*rxt(k,560)*y(k,90) + .0063005_r8*rxt(k,561)*y(k,124)
         mat(k,966) = .2067005_r8*rxt(k,537)*y(k,136)
         mat(k,2264) = .2067005_r8*rxt(k,535)*y(k,193) + .0008005_r8*rxt(k,540) &
                      *y(k,195) + .0035005_r8*rxt(k,543)*y(k,211) &
                      + .1308005_r8*rxt(k,547)*y(k,212) + .1149005_r8*rxt(k,551) &
                      *y(k,218) + .0101005_r8*rxt(k,557)*y(k,236) &
                      + .0174005_r8*rxt(k,560)*y(k,239)
         mat(k,1010) = .1149005_r8*rxt(k,553)*y(k,136)
         mat(k,1833) = .1792005_r8*rxt(k,536)*y(k,193) + .0034005_r8*rxt(k,541) &
                      *y(k,195) + .0003005_r8*rxt(k,544)*y(k,211) &
                      + .1026005_r8*rxt(k,548)*y(k,212) + .0082005_r8*rxt(k,552) &
                      *y(k,218) + .0452005_r8*rxt(k,558)*y(k,236) &
                      + .0237005_r8*rxt(k,561)*y(k,239)
         mat(k,2175) = .2067005_r8*rxt(k,537)*y(k,6) + .1149005_r8*rxt(k,553)*y(k,110)
         mat(k,46) = .1534005_r8*rxt(k,555)*y(k,228)
         mat(k,52) = .2067005_r8*rxt(k,535)*y(k,90) + .1792005_r8*rxt(k,536)*y(k,124)
         mat(k,58) = .0008005_r8*rxt(k,540)*y(k,90) + .0034005_r8*rxt(k,541)*y(k,124)
         mat(k,64) = .0035005_r8*rxt(k,543)*y(k,90) + .0003005_r8*rxt(k,544)*y(k,124)
         mat(k,70) = .1308005_r8*rxt(k,547)*y(k,90) + .1026005_r8*rxt(k,548)*y(k,124)
         mat(k,78) = .1149005_r8*rxt(k,551)*y(k,90) + .0082005_r8*rxt(k,552)*y(k,124)
         mat(k,1668) = .1534005_r8*rxt(k,555)*y(k,175)
         mat(k,84) = .0101005_r8*rxt(k,557)*y(k,90) + .0452005_r8*rxt(k,558)*y(k,124)
         mat(k,90) = .0174005_r8*rxt(k,560)*y(k,90) + .0237005_r8*rxt(k,561)*y(k,124)
         mat(k,967) = .0653005_r8*rxt(k,537)*y(k,136)
         mat(k,2265) = .0653005_r8*rxt(k,535)*y(k,193) + .0843005_r8*rxt(k,540) &
                      *y(k,195) + .0003005_r8*rxt(k,543)*y(k,211) &
                      + .0348005_r8*rxt(k,547)*y(k,212) + .0348005_r8*rxt(k,551) &
                      *y(k,218) + .0763005_r8*rxt(k,557)*y(k,236) + .086_r8*rxt(k,560) &
                      *y(k,239)
         mat(k,1011) = .0348005_r8*rxt(k,553)*y(k,136)
         mat(k,1834) = .0676005_r8*rxt(k,536)*y(k,193) + .1579005_r8*rxt(k,541) &
                      *y(k,195) + .0073005_r8*rxt(k,544)*y(k,211) &
                      + .0521005_r8*rxt(k,548)*y(k,212) + .0772005_r8*rxt(k,552) &
                      *y(k,218) + .0966005_r8*rxt(k,558)*y(k,236) &
                      + .0025005_r8*rxt(k,561)*y(k,239)
         mat(k,2176) = .0653005_r8*rxt(k,537)*y(k,6) + .0348005_r8*rxt(k,553)*y(k,110)
         mat(k,47) = .0459005_r8*rxt(k,555)*y(k,228)
         mat(k,53) = .0653005_r8*rxt(k,535)*y(k,90) + .0676005_r8*rxt(k,536)*y(k,124)
         mat(k,59) = .0843005_r8*rxt(k,540)*y(k,90) + .1579005_r8*rxt(k,541)*y(k,124)
         mat(k,65) = .0003005_r8*rxt(k,543)*y(k,90) + .0073005_r8*rxt(k,544)*y(k,124)
         mat(k,71) = .0348005_r8*rxt(k,547)*y(k,90) + .0521005_r8*rxt(k,548)*y(k,124)
         mat(k,79) = .0348005_r8*rxt(k,551)*y(k,90) + .0772005_r8*rxt(k,552)*y(k,124)
         mat(k,1669) = .0459005_r8*rxt(k,555)*y(k,175)
         mat(k,85) = .0763005_r8*rxt(k,557)*y(k,90) + .0966005_r8*rxt(k,558)*y(k,124)
         mat(k,91) = .086_r8*rxt(k,560)*y(k,90) + .0025005_r8*rxt(k,561)*y(k,124)
         mat(k,968) = .1749305_r8*rxt(k,534)*y(k,126) + .1284005_r8*rxt(k,537) &
                      *y(k,136)
         mat(k,2266) = .1284005_r8*rxt(k,535)*y(k,193) + .0443005_r8*rxt(k,540) &
                      *y(k,195) + .0271005_r8*rxt(k,543)*y(k,211) &
                      + .0076005_r8*rxt(k,547)*y(k,212) + .0554005_r8*rxt(k,551) &
                      *y(k,218) + .2157005_r8*rxt(k,557)*y(k,236) &
                      + .0512005_r8*rxt(k,560)*y(k,239)
         mat(k,873) = .0590245_r8*rxt(k,542)*y(k,126) + .0033005_r8*rxt(k,545) &
                      *y(k,136)
         mat(k,1012) = .1749305_r8*rxt(k,550)*y(k,126) + .0554005_r8*rxt(k,553) &
                      *y(k,136)
         mat(k,1835) = .079_r8*rxt(k,536)*y(k,193) + .0059005_r8*rxt(k,541)*y(k,195) &
                      + .0057005_r8*rxt(k,544)*y(k,211) + .0143005_r8*rxt(k,548) &
                      *y(k,212) + .0332005_r8*rxt(k,552)*y(k,218) &
                      + .0073005_r8*rxt(k,558)*y(k,236) + .011_r8*rxt(k,561)*y(k,239)
         mat(k,1607) = .1749305_r8*rxt(k,534)*y(k,6) + .0590245_r8*rxt(k,542)*y(k,99) &
                      + .1749305_r8*rxt(k,550)*y(k,110)
         mat(k,2177) = .1284005_r8*rxt(k,537)*y(k,6) + .0033005_r8*rxt(k,545)*y(k,99) &
                      + .0554005_r8*rxt(k,553)*y(k,110)
         mat(k,48) = .0085005_r8*rxt(k,555)*y(k,228)
         mat(k,54) = .1284005_r8*rxt(k,535)*y(k,90) + .079_r8*rxt(k,536)*y(k,124)
         mat(k,60) = .0443005_r8*rxt(k,540)*y(k,90) + .0059005_r8*rxt(k,541)*y(k,124)
         mat(k,66) = .0271005_r8*rxt(k,543)*y(k,90) + .0057005_r8*rxt(k,544)*y(k,124)
         mat(k,72) = .0076005_r8*rxt(k,547)*y(k,90) + .0143005_r8*rxt(k,548)*y(k,124)
         mat(k,80) = .0554005_r8*rxt(k,551)*y(k,90) + .0332005_r8*rxt(k,552)*y(k,124)
         mat(k,1670) = .0085005_r8*rxt(k,555)*y(k,175)
         mat(k,86) = .2157005_r8*rxt(k,557)*y(k,90) + .0073005_r8*rxt(k,558)*y(k,124)
         mat(k,92) = .0512005_r8*rxt(k,560)*y(k,90) + .011_r8*rxt(k,561)*y(k,124)
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
         mat(k,969) = .5901905_r8*rxt(k,534)*y(k,126) + .114_r8*rxt(k,537)*y(k,136)
         mat(k,2267) = .114_r8*rxt(k,535)*y(k,193) + .1621005_r8*rxt(k,540)*y(k,195) &
                      + .0474005_r8*rxt(k,543)*y(k,211) + .0113005_r8*rxt(k,547) &
                      *y(k,212) + .1278005_r8*rxt(k,551)*y(k,218) &
                      + .0738005_r8*rxt(k,557)*y(k,236) + .1598005_r8*rxt(k,560) &
                      *y(k,239)
         mat(k,874) = .0250245_r8*rxt(k,542)*y(k,126)
         mat(k,1013) = .5901905_r8*rxt(k,550)*y(k,126) + .1278005_r8*rxt(k,553) &
                      *y(k,136)
         mat(k,1836) = .1254005_r8*rxt(k,536)*y(k,193) + .0536005_r8*rxt(k,541) &
                      *y(k,195) + .0623005_r8*rxt(k,544)*y(k,211) &
                      + .0166005_r8*rxt(k,548)*y(k,212) + .130_r8*rxt(k,552)*y(k,218) &
                      + .238_r8*rxt(k,558)*y(k,236) + .1185005_r8*rxt(k,561)*y(k,239)
         mat(k,1608) = .5901905_r8*rxt(k,534)*y(k,6) + .0250245_r8*rxt(k,542)*y(k,99) &
                      + .5901905_r8*rxt(k,550)*y(k,110)
         mat(k,2178) = .114_r8*rxt(k,537)*y(k,6) + .1278005_r8*rxt(k,553)*y(k,110)
         mat(k,49) = .0128005_r8*rxt(k,555)*y(k,228)
         mat(k,55) = .114_r8*rxt(k,535)*y(k,90) + .1254005_r8*rxt(k,536)*y(k,124)
         mat(k,61) = .1621005_r8*rxt(k,540)*y(k,90) + .0536005_r8*rxt(k,541)*y(k,124)
         mat(k,67) = .0474005_r8*rxt(k,543)*y(k,90) + .0623005_r8*rxt(k,544)*y(k,124)
         mat(k,73) = .0113005_r8*rxt(k,547)*y(k,90) + .0166005_r8*rxt(k,548)*y(k,124)
         mat(k,81) = .1278005_r8*rxt(k,551)*y(k,90) + .130_r8*rxt(k,552)*y(k,124)
         mat(k,1671) = .0128005_r8*rxt(k,555)*y(k,175)
         mat(k,87) = .0738005_r8*rxt(k,557)*y(k,90) + .238_r8*rxt(k,558)*y(k,124)
         mat(k,93) = .1598005_r8*rxt(k,560)*y(k,90) + .1185005_r8*rxt(k,561)*y(k,124)
         mat(k,50) = -(rxt(k,555)*y(k,228))
         mat(k,1672) = -rxt(k,555)*y(k,175)
         mat(k,189) = .100_r8*rxt(k,469)*y(k,228)
         mat(k,207) = .230_r8*rxt(k,471)*y(k,228)
         mat(k,1691) = .100_r8*rxt(k,469)*y(k,183) + .230_r8*rxt(k,471)*y(k,185)
         mat(k,652) = -(rxt(k,493)*y(k,228))
         mat(k,1754) = -rxt(k,493)*y(k,177)
         mat(k,2305) = rxt(k,491)*y(k,233)
         mat(k,1172) = rxt(k,491)*y(k,90)
         mat(k,627) = -(rxt(k,494)*y(k,228))
         mat(k,1751) = -rxt(k,494)*y(k,178)
         mat(k,1863) = .200_r8*rxt(k,487)*y(k,223) + .200_r8*rxt(k,497)*y(k,234)
         mat(k,2375) = .500_r8*rxt(k,485)*y(k,223)
         mat(k,1103) = .200_r8*rxt(k,487)*y(k,124) + .500_r8*rxt(k,485)*y(k,201)
         mat(k,1060) = .200_r8*rxt(k,497)*y(k,124)
         mat(k,486) = -(rxt(k,498)*y(k,228))
         mat(k,1735) = -rxt(k,498)*y(k,179)
         mat(k,2297) = rxt(k,496)*y(k,234)
         mat(k,1059) = rxt(k,496)*y(k,90)
         mat(k,1072) = -(rxt(k,499)*y(k,126) + rxt(k,500)*y(k,228))
         mat(k,1623) = -rxt(k,499)*y(k,180)
         mat(k,1786) = -rxt(k,500)*y(k,180)
         mat(k,978) = .330_r8*rxt(k,480)*y(k,136)
         mat(k,1022) = .330_r8*rxt(k,483)*y(k,136)
         mat(k,1885) = .800_r8*rxt(k,487)*y(k,223) + .800_r8*rxt(k,497)*y(k,234)
         mat(k,1623) = mat(k,1623) + rxt(k,488)*y(k,223)
         mat(k,2197) = .330_r8*rxt(k,480)*y(k,6) + .330_r8*rxt(k,483)*y(k,110)
         mat(k,628) = rxt(k,494)*y(k,228)
         mat(k,2383) = .500_r8*rxt(k,485)*y(k,223) + rxt(k,495)*y(k,234)
         mat(k,1105) = .800_r8*rxt(k,487)*y(k,124) + rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,201)
         mat(k,1786) = mat(k,1786) + rxt(k,494)*y(k,178)
         mat(k,1063) = .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,201)
         mat(k,1119) = -(rxt(k,501)*y(k,228))
         mat(k,1790) = -rxt(k,501)*y(k,181)
         mat(k,981) = .300_r8*rxt(k,480)*y(k,136)
         mat(k,1025) = .300_r8*rxt(k,483)*y(k,136)
         mat(k,1889) = .900_r8*rxt(k,492)*y(k,233)
         mat(k,2200) = .300_r8*rxt(k,480)*y(k,6) + .300_r8*rxt(k,483)*y(k,110)
         mat(k,2387) = rxt(k,490)*y(k,233)
         mat(k,1175) = .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,201)
         mat(k,665) = -(rxt(k,468)*y(k,228))
         mat(k,1755) = -rxt(k,468)*y(k,182)
         mat(k,2306) = rxt(k,466)*y(k,235)
         mat(k,749) = rxt(k,466)*y(k,90)
         mat(k,187) = -(rxt(k,469)*y(k,228))
         mat(k,1689) = -rxt(k,469)*y(k,183)
         mat(k,203) = -(rxt(k,435)*y(k,228))
         mat(k,1692) = -rxt(k,435)*y(k,184)
         mat(k,2276) = rxt(k,432)*y(k,237)
         mat(k,1227) = rxt(k,432)*y(k,90)
         mat(k,208) = -(rxt(k,471)*y(k,228))
         mat(k,1693) = -rxt(k,471)*y(k,185)
         mat(k,729) = -(rxt(k,474)*y(k,228))
         mat(k,1761) = -rxt(k,474)*y(k,186)
         mat(k,2311) = rxt(k,472)*y(k,238)
         mat(k,765) = rxt(k,472)*y(k,90)
         mat(k,216) = -(rxt(k,477)*y(k,228))
         mat(k,1694) = -rxt(k,477)*y(k,187)
         mat(k,209) = .150_r8*rxt(k,471)*y(k,228)
         mat(k,1694) = mat(k,1694) + .150_r8*rxt(k,471)*y(k,185)
         mat(k,421) = -(rxt(k,478)*y(k,228))
         mat(k,1726) = -rxt(k,478)*y(k,188)
         mat(k,2288) = rxt(k,475)*y(k,240)
         mat(k,509) = rxt(k,475)*y(k,90)
         mat(k,523) = -(rxt(k,436)*y(k,90) + rxt(k,437)*y(k,124) + rxt(k,465)*y(k,125))
         mat(k,2300) = -rxt(k,436)*y(k,191)
         mat(k,1858) = -rxt(k,437)*y(k,191)
         mat(k,1932) = -rxt(k,465)*y(k,191)
         mat(k,236) = rxt(k,442)*y(k,228)
         mat(k,1739) = rxt(k,442)*y(k,22)
         mat(k,1042) = -(rxt(k,397)*y(k,90) + (rxt(k,398) + rxt(k,399)) * y(k,124))
         mat(k,2326) = -rxt(k,397)*y(k,192)
         mat(k,1882) = -(rxt(k,398) + rxt(k,399)) * y(k,192)
         mat(k,690) = rxt(k,400)*y(k,228)
         mat(k,227) = rxt(k,401)*y(k,228)
         mat(k,1783) = rxt(k,400)*y(k,2) + rxt(k,401)*y(k,15)
         mat(k,56) = -(rxt(k,535)*y(k,90) + rxt(k,536)*y(k,124))
         mat(k,2268) = -rxt(k,535)*y(k,193)
         mat(k,1837) = -rxt(k,536)*y(k,193)
         mat(k,970) = rxt(k,538)*y(k,228)
         mat(k,1673) = rxt(k,538)*y(k,6)
         mat(k,495) = -(rxt(k,439)*y(k,90) + rxt(k,440)*y(k,124))
         mat(k,2298) = -rxt(k,439)*y(k,194)
         mat(k,1855) = -rxt(k,440)*y(k,194)
         mat(k,163) = .350_r8*rxt(k,438)*y(k,228)
         mat(k,417) = rxt(k,441)*y(k,228)
         mat(k,1736) = .350_r8*rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8)
         mat(k,62) = -(rxt(k,540)*y(k,90) + rxt(k,541)*y(k,124))
         mat(k,2269) = -rxt(k,540)*y(k,195)
         mat(k,1838) = -rxt(k,541)*y(k,195)
         mat(k,159) = rxt(k,539)*y(k,228)
         mat(k,1674) = rxt(k,539)*y(k,7)
         mat(k,435) = -(rxt(k,443)*y(k,90) + rxt(k,445)*y(k,124))
         mat(k,2289) = -rxt(k,443)*y(k,196)
         mat(k,1849) = -rxt(k,445)*y(k,196)
         mat(k,346) = rxt(k,444)*y(k,228)
         mat(k,190) = .070_r8*rxt(k,469)*y(k,228)
         mat(k,210) = .060_r8*rxt(k,471)*y(k,228)
         mat(k,1728) = rxt(k,444)*y(k,23) + .070_r8*rxt(k,469)*y(k,183) &
                      + .060_r8*rxt(k,471)*y(k,185)
         mat(k,901) = -(4._r8*rxt(k,319)*y(k,197) + rxt(k,320)*y(k,201) + rxt(k,321) &
                      *y(k,90) + rxt(k,322)*y(k,124))
         mat(k,2379) = -rxt(k,320)*y(k,197)
         mat(k,2322) = -rxt(k,321)*y(k,197)
         mat(k,1877) = -rxt(k,322)*y(k,197)
         mat(k,351) = .500_r8*rxt(k,324)*y(k,228)
         mat(k,299) = rxt(k,325)*y(k,56) + rxt(k,326)*y(k,228)
         mat(k,2140) = rxt(k,325)*y(k,28)
         mat(k,1773) = .500_r8*rxt(k,324)*y(k,27) + rxt(k,326)*y(k,28)
         mat(k,935) = -(rxt(k,348)*y(k,201) + rxt(k,349)*y(k,90) + rxt(k,350)*y(k,124))
         mat(k,2380) = -rxt(k,348)*y(k,198)
         mat(k,2325) = -rxt(k,349)*y(k,198)
         mat(k,1880) = -rxt(k,350)*y(k,198)
         mat(k,410) = rxt(k,351)*y(k,228)
         mat(k,305) = rxt(k,355)*y(k,56) + rxt(k,352)*y(k,228)
         mat(k,2142) = rxt(k,355)*y(k,31)
         mat(k,1777) = rxt(k,351)*y(k,30) + rxt(k,352)*y(k,31)
         mat(k,635) = -(rxt(k,446)*y(k,90) + rxt(k,447)*y(k,124))
         mat(k,2304) = -rxt(k,446)*y(k,199)
         mat(k,1864) = -rxt(k,447)*y(k,199)
         mat(k,271) = rxt(k,448)*y(k,228)
         mat(k,2304) = mat(k,2304) + .400_r8*rxt(k,436)*y(k,191)
         mat(k,1864) = mat(k,1864) + rxt(k,437)*y(k,191)
         mat(k,2185) = rxt(k,463)*y(k,143)
         mat(k,469) = rxt(k,463)*y(k,136)
         mat(k,524) = .400_r8*rxt(k,436)*y(k,90) + rxt(k,437)*y(k,124)
         mat(k,1752) = rxt(k,448)*y(k,32)
         mat(k,1428) = -(4._r8*rxt(k,330)*y(k,200) + rxt(k,331)*y(k,201) + rxt(k,332) &
                      *y(k,90) + rxt(k,333)*y(k,124) + rxt(k,344)*y(k,125) + rxt(k,372) &
                      *y(k,213) + rxt(k,405)*y(k,208) + rxt(k,410)*y(k,209) + rxt(k,419) &
                      *y(k,210) + rxt(k,430)*y(k,237))
         mat(k,2403) = -rxt(k,331)*y(k,200)
         mat(k,2348) = -rxt(k,332)*y(k,200)
         mat(k,1906) = -rxt(k,333)*y(k,200)
         mat(k,1949) = -rxt(k,344)*y(k,200)
         mat(k,1358) = -rxt(k,372)*y(k,200)
         mat(k,1303) = -rxt(k,405)*y(k,200)
         mat(k,1336) = -rxt(k,410)*y(k,200)
         mat(k,1257) = -rxt(k,419)*y(k,200)
         mat(k,1235) = -rxt(k,430)*y(k,200)
         mat(k,985) = .060_r8*rxt(k,480)*y(k,136)
         mat(k,1153) = rxt(k,327)*y(k,126) + rxt(k,328)*y(k,228)
         mat(k,1282) = rxt(k,353)*y(k,126) + rxt(k,354)*y(k,228)
         mat(k,621) = .500_r8*rxt(k,335)*y(k,228)
         mat(k,2348) = mat(k,2348) + .450_r8*rxt(k,383)*y(k,215) + .200_r8*rxt(k,387) &
                      *y(k,217) + .150_r8*rxt(k,362)*y(k,232)
         mat(k,885) = .080_r8*rxt(k,425)*y(k,136)
         mat(k,1273) = .100_r8*rxt(k,378)*y(k,136)
         mat(k,1029) = .060_r8*rxt(k,483)*y(k,136)
         mat(k,1378) = .280_r8*rxt(k,392)*y(k,136)
         mat(k,1906) = mat(k,1906) + .530_r8*rxt(k,376)*y(k,213) + rxt(k,385)*y(k,215) &
                      + rxt(k,388)*y(k,217) + rxt(k,363)*y(k,232)
         mat(k,1645) = rxt(k,327)*y(k,45) + rxt(k,353)*y(k,49) + .530_r8*rxt(k,375) &
                      *y(k,213) + rxt(k,386)*y(k,215)
         mat(k,2216) = .060_r8*rxt(k,480)*y(k,6) + .080_r8*rxt(k,425)*y(k,99) &
                      + .100_r8*rxt(k,378)*y(k,105) + .060_r8*rxt(k,483)*y(k,110) &
                      + .280_r8*rxt(k,392)*y(k,111)
         mat(k,1122) = .650_r8*rxt(k,501)*y(k,228)
         mat(k,1428) = mat(k,1428) + .530_r8*rxt(k,372)*y(k,213)
         mat(k,2403) = mat(k,2403) + .260_r8*rxt(k,373)*y(k,213) + rxt(k,382)*y(k,215) &
                      + .300_r8*rxt(k,361)*y(k,232)
         mat(k,1358) = mat(k,1358) + .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375) &
                      *y(k,126) + .530_r8*rxt(k,372)*y(k,200) + .260_r8*rxt(k,373) &
                      *y(k,201)
         mat(k,1398) = .450_r8*rxt(k,383)*y(k,90) + rxt(k,385)*y(k,124) + rxt(k,386) &
                      *y(k,126) + rxt(k,382)*y(k,201) + 4.000_r8*rxt(k,384)*y(k,215)
         mat(k,700) = .200_r8*rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124)
         mat(k,1808) = rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) + .500_r8*rxt(k,335) &
                      *y(k,51) + .650_r8*rxt(k,501)*y(k,181)
         mat(k,1219) = .150_r8*rxt(k,362)*y(k,90) + rxt(k,363)*y(k,124) &
                      + .300_r8*rxt(k,361)*y(k,201)
         mat(k,2421) = -(rxt(k,221)*y(k,59) + (4._r8*rxt(k,298) + 4._r8*rxt(k,299) &
                      ) * y(k,201) + rxt(k,300)*y(k,90) + rxt(k,301)*y(k,124) &
                      + rxt(k,320)*y(k,197) + rxt(k,331)*y(k,200) + rxt(k,348) &
                      *y(k,198) + rxt(k,361)*y(k,232) + rxt(k,373)*y(k,213) + rxt(k,382) &
                      *y(k,215) + rxt(k,406)*y(k,208) + rxt(k,411)*y(k,209) + rxt(k,420) &
                      *y(k,210) + rxt(k,431)*y(k,237) + rxt(k,485)*y(k,223) + rxt(k,490) &
                      *y(k,233) + rxt(k,495)*y(k,234))
         mat(k,1605) = -rxt(k,221)*y(k,201)
         mat(k,2369) = -rxt(k,300)*y(k,201)
         mat(k,1925) = -rxt(k,301)*y(k,201)
         mat(k,909) = -rxt(k,320)*y(k,201)
         mat(k,1440) = -rxt(k,331)*y(k,201)
         mat(k,944) = -rxt(k,348)*y(k,201)
         mat(k,1225) = -rxt(k,361)*y(k,201)
         mat(k,1367) = -rxt(k,373)*y(k,201)
         mat(k,1408) = -rxt(k,382)*y(k,201)
         mat(k,1313) = -rxt(k,406)*y(k,201)
         mat(k,1346) = -rxt(k,411)*y(k,201)
         mat(k,1266) = -rxt(k,420)*y(k,201)
         mat(k,1243) = -rxt(k,431)*y(k,201)
         mat(k,1117) = -rxt(k,485)*y(k,201)
         mat(k,1188) = -rxt(k,490)*y(k,201)
         mat(k,1071) = -rxt(k,495)*y(k,201)
         mat(k,1149) = .280_r8*rxt(k,347)*y(k,136)
         mat(k,708) = rxt(k,334)*y(k,228)
         mat(k,458) = .700_r8*rxt(k,303)*y(k,228)
         mat(k,2261) = rxt(k,215)*y(k,56) + rxt(k,271)*y(k,73) + rxt(k,310)*y(k,224) &
                      + rxt(k,304)*y(k,228)
         mat(k,2172) = rxt(k,215)*y(k,54)
         mat(k,932) = rxt(k,271)*y(k,54)
         mat(k,2369) = mat(k,2369) + .450_r8*rxt(k,332)*y(k,200) + .330_r8*rxt(k,450) &
                      *y(k,202) + .070_r8*rxt(k,456)*y(k,216)
         mat(k,891) = .050_r8*rxt(k,425)*y(k,136)
         mat(k,1925) = mat(k,1925) + rxt(k,333)*y(k,200) + .830_r8*rxt(k,451)*y(k,202) &
                      + .170_r8*rxt(k,457)*y(k,216)
         mat(k,2236) = .280_r8*rxt(k,347)*y(k,29) + .050_r8*rxt(k,425)*y(k,99)
         mat(k,1440) = mat(k,1440) + .450_r8*rxt(k,332)*y(k,90) + rxt(k,333)*y(k,124) &
                      + 4.000_r8*rxt(k,330)*y(k,200) + .900_r8*rxt(k,331)*y(k,201) &
                      + rxt(k,405)*y(k,208) + rxt(k,410)*y(k,209) + rxt(k,419) &
                      *y(k,210) + rxt(k,372)*y(k,213) + rxt(k,381)*y(k,215) &
                      + rxt(k,430)*y(k,237)
         mat(k,2421) = mat(k,2421) + .900_r8*rxt(k,331)*y(k,200)
         mat(k,785) = .330_r8*rxt(k,450)*y(k,90) + .830_r8*rxt(k,451)*y(k,124)
         mat(k,1313) = mat(k,1313) + rxt(k,405)*y(k,200)
         mat(k,1346) = mat(k,1346) + rxt(k,410)*y(k,200)
         mat(k,1266) = mat(k,1266) + rxt(k,419)*y(k,200)
         mat(k,1367) = mat(k,1367) + rxt(k,372)*y(k,200)
         mat(k,1408) = mat(k,1408) + rxt(k,381)*y(k,200)
         mat(k,918) = .070_r8*rxt(k,456)*y(k,90) + .170_r8*rxt(k,457)*y(k,124)
         mat(k,2039) = rxt(k,310)*y(k,54)
         mat(k,1830) = rxt(k,334)*y(k,50) + .700_r8*rxt(k,303)*y(k,53) + rxt(k,304) &
                      *y(k,54)
         mat(k,1243) = mat(k,1243) + rxt(k,430)*y(k,200)
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
         mat(k,778) = -(rxt(k,450)*y(k,90) + rxt(k,451)*y(k,124) + rxt(k,452)*y(k,125))
         mat(k,2315) = -rxt(k,450)*y(k,202)
         mat(k,1870) = -rxt(k,451)*y(k,202)
         mat(k,1938) = -rxt(k,452)*y(k,202)
         mat(k,865) = -(rxt(k,579)*y(k,221) + rxt(k,580)*y(k,227) + rxt(k,581) &
                      *y(k,220))
         mat(k,846) = -rxt(k,579)*y(k,203)
         mat(k,854) = -rxt(k,580)*y(k,203)
         mat(k,680) = -rxt(k,581)*y(k,203)
         mat(k,574) = -((rxt(k,369) + rxt(k,370)) * y(k,124))
         mat(k,1860) = -(rxt(k,369) + rxt(k,370)) * y(k,204)
         mat(k,356) = rxt(k,368)*y(k,228)
         mat(k,1744) = rxt(k,368)*y(k,16)
         mat(k,460) = -(rxt(k,339)*y(k,135))
         mat(k,1526) = -rxt(k,339)*y(k,205)
         mat(k,1852) = .750_r8*rxt(k,337)*y(k,206)
         mat(k,796) = .750_r8*rxt(k,337)*y(k,124)
         mat(k,797) = -(rxt(k,336)*y(k,90) + rxt(k,337)*y(k,124))
         mat(k,2317) = -rxt(k,336)*y(k,206)
         mat(k,1871) = -rxt(k,337)*y(k,206)
         mat(k,551) = rxt(k,343)*y(k,228)
         mat(k,1766) = rxt(k,343)*y(k,25)
         mat(k,441) = -(rxt(k,307)*y(k,90) + rxt(k,309)*y(k,124))
         mat(k,2290) = -rxt(k,307)*y(k,207)
         mat(k,1850) = -rxt(k,309)*y(k,207)
         mat(k,1972) = rxt(k,294)*y(k,90)
         mat(k,2290) = mat(k,2290) + rxt(k,294)*y(k,42)
         mat(k,1299) = -(rxt(k,405)*y(k,200) + rxt(k,406)*y(k,201) + rxt(k,407) &
                      *y(k,90) + rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126))
         mat(k,1423) = -rxt(k,405)*y(k,208)
         mat(k,2398) = -rxt(k,406)*y(k,208)
         mat(k,2343) = -rxt(k,407)*y(k,208)
         mat(k,1901) = -rxt(k,408)*y(k,208)
         mat(k,1640) = -rxt(k,409)*y(k,208)
         mat(k,882) = .600_r8*rxt(k,426)*y(k,228)
         mat(k,1803) = .600_r8*rxt(k,426)*y(k,99)
         mat(k,1332) = -(rxt(k,410)*y(k,200) + rxt(k,411)*y(k,201) + rxt(k,412) &
                      *y(k,90) + rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126))
         mat(k,1424) = -rxt(k,410)*y(k,209)
         mat(k,2399) = -rxt(k,411)*y(k,209)
         mat(k,2344) = -rxt(k,412)*y(k,209)
         mat(k,1902) = -rxt(k,414)*y(k,209)
         mat(k,1641) = -rxt(k,415)*y(k,209)
         mat(k,883) = .400_r8*rxt(k,426)*y(k,228)
         mat(k,1804) = .400_r8*rxt(k,426)*y(k,99)
         mat(k,1253) = -(rxt(k,419)*y(k,200) + rxt(k,420)*y(k,201) + rxt(k,421) &
                      *y(k,90) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126))
         mat(k,1420) = -rxt(k,419)*y(k,210)
         mat(k,2395) = -rxt(k,420)*y(k,210)
         mat(k,2340) = -rxt(k,421)*y(k,210)
         mat(k,1898) = -rxt(k,422)*y(k,210)
         mat(k,1637) = -rxt(k,423)*y(k,210)
         mat(k,880) = rxt(k,418)*y(k,126)
         mat(k,1637) = mat(k,1637) + rxt(k,418)*y(k,99)
         mat(k,68) = -(rxt(k,543)*y(k,90) + rxt(k,544)*y(k,124))
         mat(k,2270) = -rxt(k,543)*y(k,211)
         mat(k,1839) = -rxt(k,544)*y(k,211)
         mat(k,875) = rxt(k,546)*y(k,228)
         mat(k,1675) = rxt(k,546)*y(k,99)
         mat(k,74) = -(rxt(k,547)*y(k,90) + rxt(k,548)*y(k,124))
         mat(k,2271) = -rxt(k,547)*y(k,212)
         mat(k,1840) = -rxt(k,548)*y(k,212)
         mat(k,75) = rxt(k,549)*y(k,228)
         mat(k,1676) = rxt(k,549)*y(k,104)
         mat(k,1356) = -(rxt(k,372)*y(k,200) + rxt(k,373)*y(k,201) + rxt(k,374) &
                      *y(k,90) + rxt(k,375)*y(k,126) + (rxt(k,376) + rxt(k,377) &
                      ) * y(k,124))
         mat(k,1425) = -rxt(k,372)*y(k,213)
         mat(k,2400) = -rxt(k,373)*y(k,213)
         mat(k,2345) = -rxt(k,374)*y(k,213)
         mat(k,1642) = -rxt(k,375)*y(k,213)
         mat(k,1903) = -(rxt(k,376) + rxt(k,377)) * y(k,213)
         mat(k,1271) = .500_r8*rxt(k,379)*y(k,228)
         mat(k,317) = .200_r8*rxt(k,380)*y(k,228)
         mat(k,1375) = rxt(k,393)*y(k,228)
         mat(k,1805) = .500_r8*rxt(k,379)*y(k,105) + .200_r8*rxt(k,380)*y(k,106) &
                      + rxt(k,393)*y(k,111)
         mat(k,740) = -(rxt(k,453)*y(k,90) + rxt(k,454)*y(k,124) + rxt(k,455)*y(k,125))
         mat(k,2312) = -rxt(k,453)*y(k,214)
         mat(k,1867) = -rxt(k,454)*y(k,214)
         mat(k,1937) = -rxt(k,455)*y(k,214)
         mat(k,1397) = -(rxt(k,381)*y(k,200) + rxt(k,382)*y(k,201) + rxt(k,383) &
                      *y(k,90) + 4._r8*rxt(k,384)*y(k,215) + rxt(k,385)*y(k,124) &
                      + rxt(k,386)*y(k,126) + rxt(k,394)*y(k,125))
         mat(k,1427) = -rxt(k,381)*y(k,215)
         mat(k,2402) = -rxt(k,382)*y(k,215)
         mat(k,2347) = -rxt(k,383)*y(k,215)
         mat(k,1905) = -rxt(k,385)*y(k,215)
         mat(k,1644) = -rxt(k,386)*y(k,215)
         mat(k,1948) = -rxt(k,394)*y(k,215)
         mat(k,1272) = .500_r8*rxt(k,379)*y(k,228)
         mat(k,318) = .500_r8*rxt(k,380)*y(k,228)
         mat(k,1807) = .500_r8*rxt(k,379)*y(k,105) + .500_r8*rxt(k,380)*y(k,106)
         mat(k,910) = -(rxt(k,456)*y(k,90) + rxt(k,457)*y(k,124) + rxt(k,458)*y(k,125))
         mat(k,2323) = -rxt(k,456)*y(k,216)
         mat(k,1878) = -rxt(k,457)*y(k,216)
         mat(k,1940) = -rxt(k,458)*y(k,216)
         mat(k,698) = -(rxt(k,387)*y(k,90) + rxt(k,388)*y(k,124))
         mat(k,2308) = -rxt(k,387)*y(k,217)
         mat(k,1866) = -rxt(k,388)*y(k,217)
         mat(k,518) = rxt(k,389)*y(k,228)
         mat(k,322) = rxt(k,390)*y(k,228)
         mat(k,1757) = rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108)
         mat(k,82) = -(rxt(k,551)*y(k,90) + rxt(k,552)*y(k,124))
         mat(k,2272) = -rxt(k,551)*y(k,218)
         mat(k,1841) = -rxt(k,552)*y(k,218)
         mat(k,1014) = rxt(k,554)*y(k,228)
         mat(k,1678) = rxt(k,554)*y(k,110)
         mat(k,529) = -(rxt(k,188)*y(k,134) + rxt(k,189)*y(k,135))
         mat(k,2068) = -rxt(k,188)*y(k,219)
         mat(k,1528) = -rxt(k,189)*y(k,219)
         mat(k,2068) = mat(k,2068) + rxt(k,583)*y(k,220)
         mat(k,860) = .900_r8*rxt(k,581)*y(k,220) + .800_r8*rxt(k,579)*y(k,221)
         mat(k,675) = rxt(k,583)*y(k,134) + .900_r8*rxt(k,581)*y(k,203)
         mat(k,844) = .800_r8*rxt(k,579)*y(k,203)
         mat(k,676) = -(rxt(k,581)*y(k,203) + rxt(k,582)*y(k,135) + (rxt(k,583) &
                      + rxt(k,584)) * y(k,134))
         mat(k,861) = -rxt(k,581)*y(k,220)
         mat(k,1529) = -rxt(k,582)*y(k,220)
         mat(k,2071) = -(rxt(k,583) + rxt(k,584)) * y(k,220)
         mat(k,845) = -(rxt(k,579)*y(k,203))
         mat(k,863) = -rxt(k,579)*y(k,221)
         mat(k,996) = rxt(k,588)*y(k,227)
         mat(k,1873) = rxt(k,590)*y(k,227)
         mat(k,2077) = rxt(k,583)*y(k,220)
         mat(k,1532) = rxt(k,587)*y(k,222)
         mat(k,678) = rxt(k,583)*y(k,134)
         mat(k,504) = rxt(k,587)*y(k,135)
         mat(k,852) = rxt(k,588)*y(k,112) + rxt(k,590)*y(k,124)
         mat(k,502) = -(rxt(k,585)*y(k,134) + (rxt(k,586) + rxt(k,587)) * y(k,135))
         mat(k,2067) = -rxt(k,585)*y(k,222)
         mat(k,1527) = -(rxt(k,586) + rxt(k,587)) * y(k,222)
         mat(k,1106) = -(rxt(k,485)*y(k,201) + rxt(k,486)*y(k,90) + rxt(k,487) &
                      *y(k,124) + rxt(k,488)*y(k,126))
         mat(k,2386) = -rxt(k,485)*y(k,223)
         mat(k,2331) = -rxt(k,486)*y(k,223)
         mat(k,1888) = -rxt(k,487)*y(k,223)
         mat(k,1626) = -rxt(k,488)*y(k,223)
         mat(k,980) = rxt(k,479)*y(k,126)
         mat(k,1024) = rxt(k,482)*y(k,126)
         mat(k,1626) = mat(k,1626) + rxt(k,479)*y(k,6) + rxt(k,482)*y(k,110) &
                      + .500_r8*rxt(k,499)*y(k,180)
         mat(k,393) = rxt(k,489)*y(k,228)
         mat(k,1073) = .500_r8*rxt(k,499)*y(k,126)
         mat(k,1789) = rxt(k,489)*y(k,128)
         mat(k,2031) = -(rxt(k,153)*y(k,77) + rxt(k,154)*y(k,241) + (rxt(k,156) &
                      + rxt(k,157)) * y(k,135) + rxt(k,158)*y(k,136) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,113) + rxt(k,239)*y(k,33) + rxt(k,240) &
                      *y(k,34) + rxt(k,241)*y(k,36) + rxt(k,242)*y(k,37) + rxt(k,243) &
                      *y(k,38) + rxt(k,244)*y(k,39) + rxt(k,245)*y(k,40) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,85) + rxt(k,266)*y(k,35) + rxt(k,267) &
                      *y(k,55) + rxt(k,268)*y(k,78) + (rxt(k,269) + rxt(k,270) &
                      ) * y(k,81) + rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65) + rxt(k,289) &
                      *y(k,41) + rxt(k,290)*y(k,43) + rxt(k,291)*y(k,82) + rxt(k,292) &
                      *y(k,83) + rxt(k,293)*y(k,84) + (rxt(k,310) + rxt(k,311) &
                      + rxt(k,312)) * y(k,54) + rxt(k,313)*y(k,86))
         mat(k,1465) = -rxt(k,153)*y(k,224)
         mat(k,2440) = -rxt(k,154)*y(k,224)
         mat(k,1549) = -(rxt(k,156) + rxt(k,157)) * y(k,224)
         mat(k,2228) = -rxt(k,158)*y(k,224)
         mat(k,259) = -(rxt(k,206) + rxt(k,207)) * y(k,224)
         mat(k,102) = -rxt(k,239)*y(k,224)
         mat(k,142) = -rxt(k,240)*y(k,224)
         mat(k,113) = -rxt(k,241)*y(k,224)
         mat(k,152) = -rxt(k,242)*y(k,224)
         mat(k,117) = -rxt(k,243)*y(k,224)
         mat(k,157) = -rxt(k,244)*y(k,224)
         mat(k,121) = -rxt(k,245)*y(k,224)
         mat(k,1502) = -(rxt(k,246) + rxt(k,247)) * y(k,224)
         mat(k,148) = -rxt(k,266)*y(k,224)
         mat(k,389) = -rxt(k,267)*y(k,224)
         mat(k,110) = -rxt(k,268)*y(k,224)
         mat(k,832) = -(rxt(k,269) + rxt(k,270)) * y(k,224)
         mat(k,240) = -rxt(k,275)*y(k,224)
         mat(k,248) = -rxt(k,276)*y(k,224)
         mat(k,483) = -rxt(k,289)*y(k,224)
         mat(k,606) = -rxt(k,290)*y(k,224)
         mat(k,244) = -rxt(k,291)*y(k,224)
         mat(k,254) = -rxt(k,292)*y(k,224)
         mat(k,295) = -rxt(k,293)*y(k,224)
         mat(k,2253) = -(rxt(k,310) + rxt(k,311) + rxt(k,312)) * y(k,224)
         mat(k,184) = -rxt(k,313)*y(k,224)
         mat(k,1549) = mat(k,1549) + rxt(k,189)*y(k,219)
         mat(k,871) = .850_r8*rxt(k,580)*y(k,227)
         mat(k,533) = rxt(k,189)*y(k,135)
         mat(k,858) = .850_r8*rxt(k,580)*y(k,203)
         mat(k,177) = -(rxt(k,160)*y(k,134) + rxt(k,161)*y(k,135))
         mat(k,2064) = -rxt(k,160)*y(k,225)
         mat(k,1524) = -rxt(k,161)*y(k,225)
         mat(k,1442) = rxt(k,162)*y(k,226)
         mat(k,2064) = mat(k,2064) + rxt(k,164)*y(k,226)
         mat(k,1524) = mat(k,1524) + rxt(k,165)*y(k,226)
         mat(k,2179) = rxt(k,166)*y(k,226)
         mat(k,179) = rxt(k,162)*y(k,63) + rxt(k,164)*y(k,134) + rxt(k,165)*y(k,135) &
                      + rxt(k,166)*y(k,136)
         mat(k,180) = -(rxt(k,162)*y(k,63) + rxt(k,164)*y(k,134) + rxt(k,165)*y(k,135) &
                      + rxt(k,166)*y(k,136))
         mat(k,1443) = -rxt(k,162)*y(k,226)
         mat(k,2065) = -rxt(k,164)*y(k,226)
         mat(k,1525) = -rxt(k,165)*y(k,226)
         mat(k,2180) = -rxt(k,166)*y(k,226)
         mat(k,1525) = mat(k,1525) + rxt(k,156)*y(k,224)
         mat(k,2007) = rxt(k,156)*y(k,135)
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
         mat(k,853) = -(rxt(k,580)*y(k,203) + rxt(k,588)*y(k,112) + rxt(k,590) &
                      *y(k,124))
         mat(k,864) = -rxt(k,580)*y(k,227)
         mat(k,997) = -rxt(k,588)*y(k,227)
         mat(k,1874) = -rxt(k,590)*y(k,227)
         mat(k,1446) = rxt(k,591)*y(k,229)
         mat(k,1533) = rxt(k,582)*y(k,220) + rxt(k,586)*y(k,222) + rxt(k,593)*y(k,229)
         mat(k,679) = rxt(k,582)*y(k,135)
         mat(k,505) = rxt(k,586)*y(k,135)
         mat(k,807) = rxt(k,591)*y(k,63) + rxt(k,593)*y(k,135)
         mat(k,1818) = -(rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,181)*y(k,90) &
                      + rxt(k,182)*y(k,134) + rxt(k,183)*y(k,136) + (4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185)) * y(k,228) + rxt(k,187)*y(k,91) + rxt(k,201) &
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
                      *y(k,141) + rxt(k,343)*y(k,25) + rxt(k,351)*y(k,30) + rxt(k,352) &
                      *y(k,31) + rxt(k,354)*y(k,49) + rxt(k,356)*y(k,96) + rxt(k,357) &
                      *y(k,127) + rxt(k,360)*y(k,148) + rxt(k,364)*y(k,149) + rxt(k,365) &
                      *y(k,29) + rxt(k,366)*y(k,48) + rxt(k,368)*y(k,16) + rxt(k,371) &
                      *y(k,94) + rxt(k,379)*y(k,105) + rxt(k,380)*y(k,106) + rxt(k,389) &
                      *y(k,107) + rxt(k,390)*y(k,108) + rxt(k,391)*y(k,109) + rxt(k,393) &
                      *y(k,111) + rxt(k,396)*y(k,1) + rxt(k,400)*y(k,2) + rxt(k,401) &
                      *y(k,15) + rxt(k,402)*y(k,95) + rxt(k,403)*y(k,97) + rxt(k,404) &
                      *y(k,98) + rxt(k,416)*y(k,100) + rxt(k,417)*y(k,101) + rxt(k,424) &
                      *y(k,102) + rxt(k,426)*y(k,99) + rxt(k,427)*y(k,103) + rxt(k,428) &
                      *y(k,115) + rxt(k,429)*y(k,116) + rxt(k,435)*y(k,184) + rxt(k,438) &
                      *y(k,7) + rxt(k,441)*y(k,8) + rxt(k,442)*y(k,22) + rxt(k,444) &
                      *y(k,23) + rxt(k,448)*y(k,32) + rxt(k,449)*y(k,66) + rxt(k,461) &
                      *y(k,144) + rxt(k,464)*y(k,145) + rxt(k,468)*y(k,182) + rxt(k,469) &
                      *y(k,183) + rxt(k,471)*y(k,185) + rxt(k,474)*y(k,186) + rxt(k,477) &
                      *y(k,187) + rxt(k,478)*y(k,188) + rxt(k,481)*y(k,6) + rxt(k,484) &
                      *y(k,110) + rxt(k,489)*y(k,128) + rxt(k,493)*y(k,177) + rxt(k,494) &
                      *y(k,178) + rxt(k,498)*y(k,179) + rxt(k,500)*y(k,180) + rxt(k,501) &
                      *y(k,181) + (rxt(k,503) + rxt(k,517)) * y(k,67) + rxt(k,505) &
                      *y(k,139) + rxt(k,507)*y(k,153) + rxt(k,511)*y(k,150) + rxt(k,516) &
                      *y(k,152) + rxt(k,519)*y(k,120))
         mat(k,1464) = -rxt(k,179)*y(k,228)
         mat(k,584) = -rxt(k,180)*y(k,228)
         mat(k,2357) = -rxt(k,181)*y(k,228)
         mat(k,2093) = -rxt(k,182)*y(k,228)
         mat(k,2224) = -rxt(k,183)*y(k,228)
         mat(k,475) = -rxt(k,187)*y(k,228)
         mat(k,1653) = -rxt(k,201)*y(k,228)
         mat(k,1003) = -rxt(k,202)*y(k,228)
         mat(k,1958) = -rxt(k,210)*y(k,228)
         mat(k,2050) = -rxt(k,211)*y(k,228)
         mat(k,959) = -rxt(k,230)*y(k,228)
         mat(k,1594) = -(rxt(k,232) + rxt(k,233)) * y(k,228)
         mat(k,1500) = -rxt(k,235)*y(k,228)
         mat(k,840) = -rxt(k,238)*y(k,228)
         mat(k,1568) = -rxt(k,262)*y(k,228)
         mat(k,831) = -rxt(k,264)*y(k,228)
         mat(k,482) = -rxt(k,278)*y(k,228)
         mat(k,605) = -rxt(k,280)*y(k,228)
         mat(k,124) = -rxt(k,281)*y(k,228)
         mat(k,374) = -rxt(k,283)*y(k,228)
         mat(k,388) = -rxt(k,285)*y(k,228)
         mat(k,243) = -rxt(k,286)*y(k,228)
         mat(k,253) = -rxt(k,287)*y(k,228)
         mat(k,294) = -rxt(k,288)*y(k,228)
         mat(k,1984) = -rxt(k,297)*y(k,228)
         mat(k,825) = -rxt(k,302)*y(k,228)
         mat(k,455) = -rxt(k,303)*y(k,228)
         mat(k,2249) = -rxt(k,304)*y(k,228)
         mat(k,183) = -rxt(k,305)*y(k,228)
         mat(k,921) = -rxt(k,306)*y(k,228)
         mat(k,1162) = -rxt(k,314)*y(k,228)
         mat(k,289) = -rxt(k,316)*y(k,228)
         mat(k,267) = -rxt(k,323)*y(k,228)
         mat(k,353) = -rxt(k,324)*y(k,228)
         mat(k,301) = -rxt(k,326)*y(k,228)
         mat(k,1155) = -rxt(k,328)*y(k,228)
         mat(k,105) = -rxt(k,329)*y(k,228)
         mat(k,707) = -rxt(k,334)*y(k,228)
         mat(k,623) = -rxt(k,335)*y(k,228)
         mat(k,1168) = -rxt(k,340)*y(k,228)
         mat(k,1057) = -rxt(k,341)*y(k,228)
         mat(k,538) = -rxt(k,342)*y(k,228)
         mat(k,554) = -rxt(k,343)*y(k,228)
         mat(k,412) = -rxt(k,351)*y(k,228)
         mat(k,307) = -rxt(k,352)*y(k,228)
         mat(k,1285) = -rxt(k,354)*y(k,228)
         mat(k,1211) = -rxt(k,356)*y(k,228)
         mat(k,895) = -rxt(k,357)*y(k,228)
         mat(k,546) = -rxt(k,360)*y(k,228)
         mat(k,400) = -rxt(k,364)*y(k,228)
         mat(k,1142) = -rxt(k,365)*y(k,228)
         mat(k,1083) = -rxt(k,366)*y(k,228)
         mat(k,360) = -rxt(k,368)*y(k,228)
         mat(k,1201) = -rxt(k,371)*y(k,228)
         mat(k,1275) = -rxt(k,379)*y(k,228)
         mat(k,319) = -rxt(k,380)*y(k,228)
         mat(k,521) = -rxt(k,389)*y(k,228)
         mat(k,325) = -rxt(k,390)*y(k,228)
         mat(k,616) = -rxt(k,391)*y(k,228)
         mat(k,1382) = -rxt(k,393)*y(k,228)
         mat(k,648) = -rxt(k,396)*y(k,228)
         mat(k,694) = -rxt(k,400)*y(k,228)
         mat(k,228) = -rxt(k,401)*y(k,228)
         mat(k,224) = -rxt(k,402)*y(k,228)
         mat(k,328) = -rxt(k,403)*y(k,228)
         mat(k,135) = -rxt(k,404)*y(k,228)
         mat(k,593) = -rxt(k,416)*y(k,228)
         mat(k,563) = -rxt(k,417)*y(k,228)
         mat(k,406) = -rxt(k,424)*y(k,228)
         mat(k,887) = -rxt(k,426)*y(k,228)
         mat(k,715) = -rxt(k,427)*y(k,228)
         mat(k,430) = -rxt(k,428)*y(k,228)
         mat(k,1095) = -rxt(k,429)*y(k,228)
         mat(k,205) = -rxt(k,435)*y(k,228)
         mat(k,164) = -rxt(k,438)*y(k,228)
         mat(k,419) = -rxt(k,441)*y(k,228)
         mat(k,237) = -rxt(k,442)*y(k,228)
         mat(k,348) = -rxt(k,444)*y(k,228)
         mat(k,272) = -rxt(k,448)*y(k,228)
         mat(k,197) = -rxt(k,449)*y(k,228)
         mat(k,173) = -rxt(k,461)*y(k,228)
         mat(k,342) = -rxt(k,464)*y(k,228)
         mat(k,673) = -rxt(k,468)*y(k,228)
         mat(k,192) = -rxt(k,469)*y(k,228)
         mat(k,214) = -rxt(k,471)*y(k,228)
         mat(k,738) = -rxt(k,474)*y(k,228)
         mat(k,219) = -rxt(k,477)*y(k,228)
         mat(k,425) = -rxt(k,478)*y(k,228)
         mat(k,988) = -rxt(k,481)*y(k,228)
         mat(k,1032) = -rxt(k,484)*y(k,228)
         mat(k,394) = -rxt(k,489)*y(k,228)
         mat(k,659) = -rxt(k,493)*y(k,228)
         mat(k,629) = -rxt(k,494)*y(k,228)
         mat(k,490) = -rxt(k,498)*y(k,228)
         mat(k,1078) = -rxt(k,500)*y(k,228)
         mat(k,1124) = -rxt(k,501)*y(k,228)
         mat(k,313) = -(rxt(k,503) + rxt(k,517)) * y(k,228)
         mat(k,368) = -rxt(k,505)*y(k,228)
         mat(k,949) = -rxt(k,507)*y(k,228)
         mat(k,721) = -rxt(k,511)*y(k,228)
         mat(k,1481) = -rxt(k,516)*y(k,228)
         mat(k,99) = -rxt(k,519)*y(k,228)
         mat(k,988) = mat(k,988) + .630_r8*rxt(k,480)*y(k,136)
         mat(k,289) = mat(k,289) + .650_r8*rxt(k,316)*y(k,228)
         mat(k,554) = mat(k,554) + .130_r8*rxt(k,318)*y(k,136)
         mat(k,353) = mat(k,353) + .500_r8*rxt(k,324)*y(k,228)
         mat(k,1142) = mat(k,1142) + .360_r8*rxt(k,347)*y(k,136)
         mat(k,1984) = mat(k,1984) + rxt(k,296)*y(k,134)
         mat(k,455) = mat(k,455) + .300_r8*rxt(k,303)*y(k,228)
         mat(k,2249) = mat(k,2249) + rxt(k,310)*y(k,224)
         mat(k,2160) = rxt(k,219)*y(k,90)
         mat(k,928) = rxt(k,273)*y(k,241)
         mat(k,2114) = 2.000_r8*rxt(k,173)*y(k,90) + rxt(k,178)*y(k,136)
         mat(k,1464) = mat(k,1464) + rxt(k,170)*y(k,134) + rxt(k,153)*y(k,224)
         mat(k,584) = mat(k,584) + rxt(k,171)*y(k,134)
         mat(k,831) = mat(k,831) + rxt(k,263)*y(k,134) + rxt(k,269)*y(k,224)
         mat(k,1500) = mat(k,1500) + rxt(k,234)*y(k,134) + rxt(k,246)*y(k,224)
         mat(k,183) = mat(k,183) + rxt(k,313)*y(k,224)
         mat(k,2357) = mat(k,2357) + rxt(k,219)*y(k,56) + 2.000_r8*rxt(k,173)*y(k,76) &
                      + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126) + rxt(k,176) &
                      *y(k,134) + rxt(k,177)*y(k,136) + .400_r8*rxt(k,436)*y(k,191) &
                      + .450_r8*rxt(k,332)*y(k,200) + .400_r8*rxt(k,450)*y(k,202) &
                      + .450_r8*rxt(k,383)*y(k,215) + .400_r8*rxt(k,456)*y(k,216) &
                      + .200_r8*rxt(k,387)*y(k,217) + .150_r8*rxt(k,362)*y(k,232)
         mat(k,791) = rxt(k,265)*y(k,134)
         mat(k,840) = mat(k,840) + rxt(k,237)*y(k,134)
         mat(k,887) = mat(k,887) + .320_r8*rxt(k,425)*y(k,136)
         mat(k,715) = mat(k,715) + .600_r8*rxt(k,427)*y(k,228)
         mat(k,1275) = mat(k,1275) + .240_r8*rxt(k,378)*y(k,136)
         mat(k,319) = mat(k,319) + .100_r8*rxt(k,380)*y(k,228)
         mat(k,1032) = mat(k,1032) + .630_r8*rxt(k,483)*y(k,136)
         mat(k,1382) = mat(k,1382) + .360_r8*rxt(k,392)*y(k,136)
         mat(k,1913) = rxt(k,203)*y(k,90)
         mat(k,1653) = mat(k,1653) + rxt(k,198)*y(k,90)
         mat(k,2093) = mat(k,2093) + rxt(k,296)*y(k,42) + rxt(k,170)*y(k,77) &
                      + rxt(k,171)*y(k,79) + rxt(k,263)*y(k,81) + rxt(k,234)*y(k,85) &
                      + rxt(k,176)*y(k,90) + rxt(k,265)*y(k,92) + rxt(k,237)*y(k,93)
         mat(k,2224) = mat(k,2224) + .630_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,318) &
                      *y(k,25) + .360_r8*rxt(k,347)*y(k,29) + rxt(k,178)*y(k,76) &
                      + rxt(k,177)*y(k,90) + .320_r8*rxt(k,425)*y(k,99) &
                      + .240_r8*rxt(k,378)*y(k,105) + .630_r8*rxt(k,483)*y(k,110) &
                      + .360_r8*rxt(k,392)*y(k,111)
         mat(k,546) = mat(k,546) + .500_r8*rxt(k,360)*y(k,228)
         mat(k,205) = mat(k,205) + .500_r8*rxt(k,435)*y(k,228)
         mat(k,525) = .400_r8*rxt(k,436)*y(k,90)
         mat(k,1432) = .450_r8*rxt(k,332)*y(k,90)
         mat(k,781) = .400_r8*rxt(k,450)*y(k,90)
         mat(k,1401) = .450_r8*rxt(k,383)*y(k,90)
         mat(k,914) = .400_r8*rxt(k,456)*y(k,90)
         mat(k,701) = .200_r8*rxt(k,387)*y(k,90)
         mat(k,2027) = rxt(k,310)*y(k,54) + rxt(k,153)*y(k,77) + rxt(k,269)*y(k,81) &
                      + rxt(k,246)*y(k,85) + rxt(k,313)*y(k,86) + 2.000_r8*rxt(k,154) &
                      *y(k,241)
         mat(k,1818) = mat(k,1818) + .650_r8*rxt(k,316)*y(k,24) + .500_r8*rxt(k,324) &
                      *y(k,27) + .300_r8*rxt(k,303)*y(k,53) + .600_r8*rxt(k,427) &
                      *y(k,103) + .100_r8*rxt(k,380)*y(k,106) + .500_r8*rxt(k,360) &
                      *y(k,148) + .500_r8*rxt(k,435)*y(k,184)
         mat(k,1220) = .150_r8*rxt(k,362)*y(k,90)
         mat(k,2436) = rxt(k,273)*y(k,73) + 2.000_r8*rxt(k,154)*y(k,224)
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
         mat(k,805) = -(rxt(k,591)*y(k,63) + rxt(k,593)*y(k,135))
         mat(k,1444) = -rxt(k,591)*y(k,229)
         mat(k,1531) = -rxt(k,593)*y(k,229)
         mat(k,2074) = rxt(k,584)*y(k,220) + rxt(k,585)*y(k,222)
         mat(k,677) = rxt(k,584)*y(k,134)
         mat(k,503) = rxt(k,585)*y(k,134)
         mat(k,448) = -(rxt(k,459)*y(k,90) + rxt(k,460)*y(k,124))
         mat(k,2291) = -rxt(k,459)*y(k,230)
         mat(k,1851) = -rxt(k,460)*y(k,230)
         mat(k,195) = .200_r8*rxt(k,449)*y(k,228)
         mat(k,171) = .140_r8*rxt(k,461)*y(k,228)
         mat(k,340) = rxt(k,464)*y(k,228)
         mat(k,1729) = .200_r8*rxt(k,449)*y(k,66) + .140_r8*rxt(k,461)*y(k,144) &
                      + rxt(k,464)*y(k,145)
         mat(k,814) = -(rxt(k,358)*y(k,90) + rxt(k,359)*y(k,124))
         mat(k,2318) = -rxt(k,358)*y(k,231)
         mat(k,1872) = -rxt(k,359)*y(k,231)
         mat(k,1130) = rxt(k,365)*y(k,228)
         mat(k,543) = .500_r8*rxt(k,360)*y(k,228)
         mat(k,1767) = rxt(k,365)*y(k,29) + .500_r8*rxt(k,360)*y(k,148)
         mat(k,1217) = -(rxt(k,361)*y(k,201) + rxt(k,362)*y(k,90) + rxt(k,363) &
                      *y(k,124))
         mat(k,2393) = -rxt(k,361)*y(k,232)
         mat(k,2338) = -rxt(k,362)*y(k,232)
         mat(k,1896) = -rxt(k,363)*y(k,232)
         mat(k,983) = .060_r8*rxt(k,480)*y(k,136)
         mat(k,1081) = rxt(k,366)*y(k,228)
         mat(k,1027) = .060_r8*rxt(k,483)*y(k,136)
         mat(k,2207) = .060_r8*rxt(k,480)*y(k,6) + .060_r8*rxt(k,483)*y(k,110)
         mat(k,398) = rxt(k,364)*y(k,228)
         mat(k,1121) = .150_r8*rxt(k,501)*y(k,228)
         mat(k,1798) = rxt(k,366)*y(k,48) + rxt(k,364)*y(k,149) + .150_r8*rxt(k,501) &
                      *y(k,181)
         mat(k,1178) = -(rxt(k,490)*y(k,201) + rxt(k,491)*y(k,90) + rxt(k,492) &
                      *y(k,124))
         mat(k,2391) = -rxt(k,490)*y(k,233)
         mat(k,2336) = -rxt(k,491)*y(k,233)
         mat(k,1893) = -rxt(k,492)*y(k,233)
         mat(k,1632) = .500_r8*rxt(k,499)*y(k,180)
         mat(k,657) = rxt(k,493)*y(k,228)
         mat(k,1076) = .500_r8*rxt(k,499)*y(k,126) + rxt(k,500)*y(k,228)
         mat(k,1795) = rxt(k,493)*y(k,177) + rxt(k,500)*y(k,180)
         mat(k,1062) = -(rxt(k,495)*y(k,201) + rxt(k,496)*y(k,90) + rxt(k,497) &
                      *y(k,124))
         mat(k,2382) = -rxt(k,495)*y(k,234)
         mat(k,2328) = -rxt(k,496)*y(k,234)
         mat(k,1884) = -rxt(k,497)*y(k,234)
         mat(k,977) = rxt(k,481)*y(k,228)
         mat(k,1021) = rxt(k,484)*y(k,228)
         mat(k,487) = rxt(k,498)*y(k,228)
         mat(k,1785) = rxt(k,481)*y(k,6) + rxt(k,484)*y(k,110) + rxt(k,498)*y(k,179)
         mat(k,751) = -(rxt(k,466)*y(k,90) + rxt(k,467)*y(k,124))
         mat(k,2313) = -rxt(k,466)*y(k,235)
         mat(k,1868) = -rxt(k,467)*y(k,235)
         mat(k,667) = rxt(k,468)*y(k,228)
         mat(k,191) = .650_r8*rxt(k,469)*y(k,228)
         mat(k,1763) = rxt(k,468)*y(k,182) + .650_r8*rxt(k,469)*y(k,183)
         mat(k,88) = -(rxt(k,557)*y(k,90) + rxt(k,558)*y(k,124))
         mat(k,2273) = -rxt(k,557)*y(k,236)
         mat(k,1842) = -rxt(k,558)*y(k,236)
         mat(k,186) = rxt(k,556)*y(k,228)
         mat(k,1679) = rxt(k,556)*y(k,183)
         mat(k,1233) = -(rxt(k,430)*y(k,200) + rxt(k,431)*y(k,201) + rxt(k,432) &
                      *y(k,90) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126))
         mat(k,1419) = -rxt(k,430)*y(k,237)
         mat(k,2394) = -rxt(k,431)*y(k,237)
         mat(k,2339) = -rxt(k,432)*y(k,237)
         mat(k,1897) = -rxt(k,433)*y(k,237)
         mat(k,1636) = -rxt(k,434)*y(k,237)
         mat(k,223) = rxt(k,402)*y(k,228)
         mat(k,327) = rxt(k,403)*y(k,228)
         mat(k,134) = rxt(k,404)*y(k,228)
         mat(k,712) = .400_r8*rxt(k,427)*y(k,228)
         mat(k,204) = .500_r8*rxt(k,435)*y(k,228)
         mat(k,1799) = rxt(k,402)*y(k,95) + rxt(k,403)*y(k,97) + rxt(k,404)*y(k,98) &
                      + .400_r8*rxt(k,427)*y(k,103) + .500_r8*rxt(k,435)*y(k,184)
         mat(k,767) = -(rxt(k,472)*y(k,90) + rxt(k,473)*y(k,124))
         mat(k,2314) = -rxt(k,472)*y(k,238)
         mat(k,1869) = -rxt(k,473)*y(k,238)
         mat(k,211) = .560_r8*rxt(k,471)*y(k,228)
         mat(k,731) = rxt(k,474)*y(k,228)
         mat(k,1764) = .560_r8*rxt(k,471)*y(k,185) + rxt(k,474)*y(k,186)
         mat(k,94) = -(rxt(k,560)*y(k,90) + rxt(k,561)*y(k,124))
         mat(k,2274) = -rxt(k,560)*y(k,239)
         mat(k,1843) = -rxt(k,561)*y(k,239)
         mat(k,206) = rxt(k,559)*y(k,228)
         mat(k,1680) = rxt(k,559)*y(k,185)
         mat(k,510) = -(rxt(k,475)*y(k,90) + rxt(k,476)*y(k,124))
         mat(k,2299) = -rxt(k,475)*y(k,240)
         mat(k,1856) = -rxt(k,476)*y(k,240)
         mat(k,218) = .300_r8*rxt(k,477)*y(k,228)
         mat(k,422) = rxt(k,478)*y(k,228)
         mat(k,1737) = .300_r8*rxt(k,477)*y(k,187) + rxt(k,478)*y(k,188)
         mat(k,2449) = -(rxt(k,154)*y(k,224) + rxt(k,273)*y(k,73) + rxt(k,518) &
                      *y(k,154))
         mat(k,2040) = -rxt(k,154)*y(k,241)
         mat(k,933) = -rxt(k,273)*y(k,241)
         mat(k,264) = -rxt(k,518)*y(k,241)
         mat(k,303) = rxt(k,326)*y(k,228)
         mat(k,414) = rxt(k,351)*y(k,228)
         mat(k,309) = rxt(k,352)*y(k,228)
         mat(k,485) = rxt(k,278)*y(k,228)
         mat(k,1997) = rxt(k,297)*y(k,228)
         mat(k,610) = rxt(k,280)*y(k,228)
         mat(k,126) = rxt(k,281)*y(k,228)
         mat(k,1159) = rxt(k,328)*y(k,228)
         mat(k,378) = rxt(k,283)*y(k,228)
         mat(k,1085) = rxt(k,366)*y(k,228)
         mat(k,1288) = rxt(k,354)*y(k,228)
         mat(k,709) = rxt(k,334)*y(k,228)
         mat(k,626) = rxt(k,335)*y(k,228)
         mat(k,459) = rxt(k,303)*y(k,228)
         mat(k,2262) = rxt(k,304)*y(k,228)
         mat(k,2127) = rxt(k,174)*y(k,90)
         mat(k,1472) = rxt(k,179)*y(k,228)
         mat(k,588) = rxt(k,180)*y(k,228)
         mat(k,835) = rxt(k,264)*y(k,228)
         mat(k,297) = rxt(k,288)*y(k,228)
         mat(k,1507) = (rxt(k,570)+rxt(k,575))*y(k,92) + (rxt(k,563)+rxt(k,569) &
                       +rxt(k,574))*y(k,93) + rxt(k,235)*y(k,228)
         mat(k,923) = rxt(k,306)*y(k,228)
         mat(k,2063) = rxt(k,211)*y(k,228)
         mat(k,2370) = rxt(k,174)*y(k,76) + rxt(k,181)*y(k,228)
         mat(k,478) = rxt(k,187)*y(k,228)
         mat(k,794) = (rxt(k,570)+rxt(k,575))*y(k,85)
         mat(k,843) = (rxt(k,563)+rxt(k,569)+rxt(k,574))*y(k,85) + rxt(k,238)*y(k,228)
         mat(k,1279) = .500_r8*rxt(k,379)*y(k,228)
         mat(k,100) = rxt(k,519)*y(k,228)
         mat(k,549) = rxt(k,360)*y(k,228)
         mat(k,402) = rxt(k,364)*y(k,228)
         mat(k,1831) = rxt(k,326)*y(k,28) + rxt(k,351)*y(k,30) + rxt(k,352)*y(k,31) &
                      + rxt(k,278)*y(k,41) + rxt(k,297)*y(k,42) + rxt(k,280)*y(k,43) &
                      + rxt(k,281)*y(k,44) + rxt(k,328)*y(k,45) + rxt(k,283)*y(k,46) &
                      + rxt(k,366)*y(k,48) + rxt(k,354)*y(k,49) + rxt(k,334)*y(k,50) &
                      + rxt(k,335)*y(k,51) + rxt(k,303)*y(k,53) + rxt(k,304)*y(k,54) &
                      + rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,264)*y(k,81) &
                      + rxt(k,288)*y(k,84) + rxt(k,235)*y(k,85) + rxt(k,306)*y(k,87) &
                      + rxt(k,211)*y(k,89) + rxt(k,181)*y(k,90) + rxt(k,187)*y(k,91) &
                      + rxt(k,238)*y(k,93) + .500_r8*rxt(k,379)*y(k,105) + rxt(k,519) &
                      *y(k,120) + rxt(k,360)*y(k,148) + rxt(k,364)*y(k,149) &
                      + 2.000_r8*rxt(k,184)*y(k,228)
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
         mat(k, 43) = lmat(k, 43)
         mat(k, 44) = lmat(k, 44)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 68) = mat(k, 68) + lmat(k, 68)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 82) = mat(k, 82) + lmat(k, 82)
         mat(k, 88) = mat(k, 88) + lmat(k, 88)
         mat(k, 94) = mat(k, 94) + lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = lmat(k, 97)
         mat(k, 98) = mat(k, 98) + lmat(k, 98)
         mat(k, 101) = mat(k, 101) + lmat(k, 101)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 104) = mat(k, 104) + lmat(k, 104)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = lmat(k, 128)
         mat(k, 129) = lmat(k, 129)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 136) = lmat(k, 136)
         mat(k, 137) = lmat(k, 137)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = mat(k, 139) + lmat(k, 139)
         mat(k, 140) = mat(k, 140) + lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = lmat(k, 169)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = mat(k, 182) + lmat(k, 182)
         mat(k, 187) = mat(k, 187) + lmat(k, 187)
         mat(k, 194) = mat(k, 194) + lmat(k, 194)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = mat(k, 203) + lmat(k, 203)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 224) = mat(k, 224) + lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 229) = lmat(k, 229)
         mat(k, 230) = lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 238) = mat(k, 238) + lmat(k, 238)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = lmat(k, 262)
         mat(k, 263) = lmat(k, 263)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 272) = mat(k, 272) + lmat(k, 272)
         mat(k, 273) = lmat(k, 273)
         mat(k, 274) = lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 291) = mat(k, 291) + lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 329) = lmat(k, 329)
         mat(k, 330) = lmat(k, 330)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = mat(k, 339) + lmat(k, 339)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 343) = lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = mat(k, 353) + lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 364) = lmat(k, 364)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 377) = lmat(k, 377)
         mat(k, 379) = lmat(k, 379)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 395) = lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 399) = lmat(k, 399)
         mat(k, 400) = mat(k, 400) + lmat(k, 400)
         mat(k, 401) = lmat(k, 401)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 404) = lmat(k, 404)
         mat(k, 407) = lmat(k, 407)
         mat(k, 408) = mat(k, 408) + lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 411) = lmat(k, 411)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 413) = lmat(k, 413)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 416) = lmat(k, 416)
         mat(k, 418) = lmat(k, 418)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 431) = lmat(k, 431)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 445) = lmat(k, 445)
         mat(k, 446) = mat(k, 446) + lmat(k, 446)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = lmat(k, 457)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 465) = lmat(k, 465)
         mat(k, 466) = lmat(k, 466)
         mat(k, 467) = lmat(k, 467)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 472) = mat(k, 472) + lmat(k, 472)
         mat(k, 474) = lmat(k, 474)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 477) = lmat(k, 477)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 486) = mat(k, 486) + lmat(k, 486)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = mat(k, 490) + lmat(k, 490)
         mat(k, 491) = lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 502) = mat(k, 502) + lmat(k, 502)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 519) = lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 529) = mat(k, 529) + lmat(k, 529)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 535) = lmat(k, 535)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 539) = lmat(k, 539)
         mat(k, 541) = lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 544) = lmat(k, 544)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 547) = lmat(k, 547)
         mat(k, 548) = lmat(k, 548)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 564) = lmat(k, 564)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 567) = lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 569) = mat(k, 569) + lmat(k, 569)
         mat(k, 570) = mat(k, 570) + lmat(k, 570)
         mat(k, 571) = lmat(k, 571)
         mat(k, 572) = lmat(k, 572)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 584) = mat(k, 584) + lmat(k, 584)
         mat(k, 589) = mat(k, 589) + lmat(k, 589)
         mat(k, 596) = lmat(k, 596)
         mat(k, 598) = lmat(k, 598)
         mat(k, 599) = lmat(k, 599)
         mat(k, 600) = lmat(k, 600)
         mat(k, 601) = lmat(k, 601)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 609) = lmat(k, 609)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 613) = lmat(k, 613)
         mat(k, 617) = lmat(k, 617)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 625) = lmat(k, 625)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 631) = lmat(k, 631)
         mat(k, 632) = lmat(k, 632)
         mat(k, 635) = mat(k, 635) + lmat(k, 635)
         mat(k, 641) = lmat(k, 641)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 645) = mat(k, 645) + lmat(k, 645)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 650) = mat(k, 650) + lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = lmat(k, 655)
         mat(k, 656) = lmat(k, 656)
         mat(k, 658) = lmat(k, 658)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = lmat(k, 661)
         mat(k, 662) = lmat(k, 662)
         mat(k, 663) = lmat(k, 663)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 670) = lmat(k, 670)
         mat(k, 672) = lmat(k, 672)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 674) = lmat(k, 674)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 686) = lmat(k, 686)
         mat(k, 687) = mat(k, 687) + lmat(k, 687)
         mat(k, 691) = lmat(k, 691)
         mat(k, 692) = lmat(k, 692)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 695) = lmat(k, 695)
         mat(k, 696) = lmat(k, 696)
         mat(k, 698) = mat(k, 698) + lmat(k, 698)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 711) = mat(k, 711) + lmat(k, 711)
         mat(k, 713) = lmat(k, 713)
         mat(k, 714) = lmat(k, 714)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 716) = lmat(k, 716)
         mat(k, 717) = lmat(k, 717)
         mat(k, 718) = mat(k, 718) + lmat(k, 718)
         mat(k, 725) = lmat(k, 725)
         mat(k, 726) = lmat(k, 726)
         mat(k, 727) = lmat(k, 727)
         mat(k, 728) = lmat(k, 728)
         mat(k, 729) = mat(k, 729) + lmat(k, 729)
         mat(k, 734) = lmat(k, 734)
         mat(k, 736) = lmat(k, 736)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 739) = lmat(k, 739)
         mat(k, 740) = mat(k, 740) + lmat(k, 740)
         mat(k, 751) = mat(k, 751) + lmat(k, 751)
         mat(k, 767) = mat(k, 767) + lmat(k, 767)
         mat(k, 778) = mat(k, 778) + lmat(k, 778)
         mat(k, 787) = mat(k, 787) + lmat(k, 787)
         mat(k, 789) = lmat(k, 789)
         mat(k, 791) = mat(k, 791) + lmat(k, 791)
         mat(k, 797) = mat(k, 797) + lmat(k, 797)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 806) = lmat(k, 806)
         mat(k, 808) = lmat(k, 808)
         mat(k, 814) = mat(k, 814) + lmat(k, 814)
         mat(k, 824) = mat(k, 824) + lmat(k, 824)
         mat(k, 828) = mat(k, 828) + lmat(k, 828)
         mat(k, 829) = mat(k, 829) + lmat(k, 829)
         mat(k, 834) = mat(k, 834) + lmat(k, 834)
         mat(k, 837) = mat(k, 837) + lmat(k, 837)
         mat(k, 840) = mat(k, 840) + lmat(k, 840)
         mat(k, 842) = mat(k, 842) + lmat(k, 842)
         mat(k, 845) = mat(k, 845) + lmat(k, 845)
         mat(k, 852) = mat(k, 852) + lmat(k, 852)
         mat(k, 853) = mat(k, 853) + lmat(k, 853)
         mat(k, 857) = mat(k, 857) + lmat(k, 857)
         mat(k, 865) = mat(k, 865) + lmat(k, 865)
         mat(k, 876) = mat(k, 876) + lmat(k, 876)
         mat(k, 892) = mat(k, 892) + lmat(k, 892)
         mat(k, 894) = lmat(k, 894)
         mat(k, 896) = mat(k, 896) + lmat(k, 896)
         mat(k, 897) = lmat(k, 897)
         mat(k, 901) = mat(k, 901) + lmat(k, 901)
         mat(k, 910) = mat(k, 910) + lmat(k, 910)
         mat(k, 919) = mat(k, 919) + lmat(k, 919)
         mat(k, 925) = mat(k, 925) + lmat(k, 925)
         mat(k, 935) = mat(k, 935) + lmat(k, 935)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 948) = lmat(k, 948)
         mat(k, 950) = lmat(k, 950)
         mat(k, 954) = mat(k, 954) + lmat(k, 954)
         mat(k, 955) = mat(k, 955) + lmat(k, 955)
         mat(k, 957) = mat(k, 957) + lmat(k, 957)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 960) = lmat(k, 960)
         mat(k, 961) = mat(k, 961) + lmat(k, 961)
         mat(k, 963) = mat(k, 963) + lmat(k, 963)
         mat(k, 974) = mat(k, 974) + lmat(k, 974)
         mat(k, 994) = lmat(k, 994)
         mat(k, 998) = lmat(k, 998)
         mat(k, 999) = mat(k, 999) + lmat(k, 999)
         mat(k,1018) = mat(k,1018) + lmat(k,1018)
         mat(k,1042) = mat(k,1042) + lmat(k,1042)
         mat(k,1053) = lmat(k,1053)
         mat(k,1054) = mat(k,1054) + lmat(k,1054)
         mat(k,1055) = mat(k,1055) + lmat(k,1055)
         mat(k,1058) = mat(k,1058) + lmat(k,1058)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1074) = lmat(k,1074)
         mat(k,1075) = lmat(k,1075)
         mat(k,1079) = lmat(k,1079)
         mat(k,1080) = mat(k,1080) + lmat(k,1080)
         mat(k,1082) = lmat(k,1082)
         mat(k,1084) = lmat(k,1084)
         mat(k,1086) = lmat(k,1086)
         mat(k,1090) = mat(k,1090) + lmat(k,1090)
         mat(k,1097) = lmat(k,1097)
         mat(k,1099) = lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1106) = mat(k,1106) + lmat(k,1106)
         mat(k,1118) = mat(k,1118) + lmat(k,1118)
         mat(k,1119) = mat(k,1119) + lmat(k,1119)
         mat(k,1120) = mat(k,1120) + lmat(k,1120)
         mat(k,1121) = mat(k,1121) + lmat(k,1121)
         mat(k,1122) = mat(k,1122) + lmat(k,1122)
         mat(k,1123) = mat(k,1123) + lmat(k,1123)
         mat(k,1125) = mat(k,1125) + lmat(k,1125)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1133) = mat(k,1133) + lmat(k,1133)
         mat(k,1151) = mat(k,1151) + lmat(k,1151)
         mat(k,1152) = lmat(k,1152)
         mat(k,1157) = lmat(k,1157)
         mat(k,1158) = lmat(k,1158)
         mat(k,1160) = mat(k,1160) + lmat(k,1160)
         mat(k,1165) = lmat(k,1165)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1169) = mat(k,1169) + lmat(k,1169)
         mat(k,1170) = mat(k,1170) + lmat(k,1170)
         mat(k,1178) = mat(k,1178) + lmat(k,1178)
         mat(k,1191) = lmat(k,1191)
         mat(k,1192) = lmat(k,1192)
         mat(k,1193) = lmat(k,1193)
         mat(k,1194) = lmat(k,1194)
         mat(k,1195) = mat(k,1195) + lmat(k,1195)
         mat(k,1196) = lmat(k,1196)
         mat(k,1198) = lmat(k,1198)
         mat(k,1202) = lmat(k,1202)
         mat(k,1203) = lmat(k,1203)
         mat(k,1204) = lmat(k,1204)
         mat(k,1205) = mat(k,1205) + lmat(k,1205)
         mat(k,1208) = mat(k,1208) + lmat(k,1208)
         mat(k,1210) = lmat(k,1210)
         mat(k,1212) = lmat(k,1212)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1253) = mat(k,1253) + lmat(k,1253)
         mat(k,1268) = mat(k,1268) + lmat(k,1268)
         mat(k,1269) = mat(k,1269) + lmat(k,1269)
         mat(k,1272) = mat(k,1272) + lmat(k,1272)
         mat(k,1273) = mat(k,1273) + lmat(k,1273)
         mat(k,1276) = mat(k,1276) + lmat(k,1276)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1280) = mat(k,1280) + lmat(k,1280)
         mat(k,1281) = mat(k,1281) + lmat(k,1281)
         mat(k,1282) = mat(k,1282) + lmat(k,1282)
         mat(k,1287) = lmat(k,1287)
         mat(k,1299) = mat(k,1299) + lmat(k,1299)
         mat(k,1315) = lmat(k,1315)
         mat(k,1332) = mat(k,1332) + lmat(k,1332)
         mat(k,1345) = mat(k,1345) + lmat(k,1345)
         mat(k,1356) = mat(k,1356) + lmat(k,1356)
         mat(k,1370) = lmat(k,1370)
         mat(k,1372) = mat(k,1372) + lmat(k,1372)
         mat(k,1376) = mat(k,1376) + lmat(k,1376)
         mat(k,1378) = mat(k,1378) + lmat(k,1378)
         mat(k,1390) = lmat(k,1390)
         mat(k,1397) = mat(k,1397) + lmat(k,1397)
         mat(k,1428) = mat(k,1428) + lmat(k,1428)
         mat(k,1449) = mat(k,1449) + lmat(k,1449)
         mat(k,1450) = mat(k,1450) + lmat(k,1450)
         mat(k,1456) = lmat(k,1456)
         mat(k,1461) = mat(k,1461) + lmat(k,1461)
         mat(k,1474) = lmat(k,1474)
         mat(k,1476) = mat(k,1476) + lmat(k,1476)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1495) = mat(k,1495) + lmat(k,1495)
         mat(k,1505) = mat(k,1505) + lmat(k,1505)
         mat(k,1506) = mat(k,1506) + lmat(k,1506)
         mat(k,1511) = mat(k,1511) + lmat(k,1511)
         mat(k,1531) = mat(k,1531) + lmat(k,1531)
         mat(k,1533) = mat(k,1533) + lmat(k,1533)
         mat(k,1534) = lmat(k,1534)
         mat(k,1542) = mat(k,1542) + lmat(k,1542)
         mat(k,1549) = mat(k,1549) + lmat(k,1549)
         mat(k,1550) = mat(k,1550) + lmat(k,1550)
         mat(k,1563) = mat(k,1563) + lmat(k,1563)
         mat(k,1565) = mat(k,1565) + lmat(k,1565)
         mat(k,1574) = mat(k,1574) + lmat(k,1574)
         mat(k,1592) = mat(k,1592) + lmat(k,1592)
         mat(k,1600) = mat(k,1600) + lmat(k,1600)
         mat(k,1602) = mat(k,1602) + lmat(k,1602)
         mat(k,1649) = mat(k,1649) + lmat(k,1649)
         mat(k,1652) = mat(k,1652) + lmat(k,1652)
         mat(k,1654) = mat(k,1654) + lmat(k,1654)
         mat(k,1655) = mat(k,1655) + lmat(k,1655)
         mat(k,1658) = mat(k,1658) + lmat(k,1658)
         mat(k,1659) = mat(k,1659) + lmat(k,1659)
         mat(k,1818) = mat(k,1818) + lmat(k,1818)
         mat(k,1873) = mat(k,1873) + lmat(k,1873)
         mat(k,1875) = lmat(k,1875)
         mat(k,1881) = mat(k,1881) + lmat(k,1881)
         mat(k,1914) = mat(k,1914) + lmat(k,1914)
         mat(k,1919) = mat(k,1919) + lmat(k,1919)
         mat(k,1958) = mat(k,1958) + lmat(k,1958)
         mat(k,1959) = mat(k,1959) + lmat(k,1959)
         mat(k,1960) = mat(k,1960) + lmat(k,1960)
         mat(k,1963) = mat(k,1963) + lmat(k,1963)
         mat(k,1964) = mat(k,1964) + lmat(k,1964)
         mat(k,1975) = mat(k,1975) + lmat(k,1975)
         mat(k,1977) = lmat(k,1977)
         mat(k,1987) = mat(k,1987) + lmat(k,1987)
         mat(k,1991) = mat(k,1991) + lmat(k,1991)
         mat(k,2031) = mat(k,2031) + lmat(k,2031)
         mat(k,2033) = mat(k,2033) + lmat(k,2033)
         mat(k,2050) = mat(k,2050) + lmat(k,2050)
         mat(k,2052) = lmat(k,2052)
         mat(k,2055) = mat(k,2055) + lmat(k,2055)
         mat(k,2074) = mat(k,2074) + lmat(k,2074)
         mat(k,2079) = lmat(k,2079)
         mat(k,2099) = mat(k,2099) + lmat(k,2099)
         mat(k,2121) = mat(k,2121) + lmat(k,2121)
         mat(k,2168) = mat(k,2168) + lmat(k,2168)
         mat(k,2179) = mat(k,2179) + lmat(k,2179)
         mat(k,2220) = mat(k,2220) + lmat(k,2220)
         mat(k,2228) = mat(k,2228) + lmat(k,2228)
         mat(k,2230) = mat(k,2230) + lmat(k,2230)
         mat(k,2233) = mat(k,2233) + lmat(k,2233)
         mat(k,2240) = lmat(k,2240)
         mat(k,2241) = lmat(k,2241)
         mat(k,2242) = mat(k,2242) + lmat(k,2242)
         mat(k,2249) = mat(k,2249) + lmat(k,2249)
         mat(k,2252) = mat(k,2252) + lmat(k,2252)
         mat(k,2255) = lmat(k,2255)
         mat(k,2256) = mat(k,2256) + lmat(k,2256)
         mat(k,2259) = mat(k,2259) + lmat(k,2259)
         mat(k,2261) = mat(k,2261) + lmat(k,2261)
         mat(k,2262) = mat(k,2262) + lmat(k,2262)
         mat(k,2368) = mat(k,2368) + lmat(k,2368)
         mat(k,2370) = mat(k,2370) + lmat(k,2370)
         mat(k,2421) = mat(k,2421) + lmat(k,2421)
         mat(k,2428) = lmat(k,2428)
         mat(k,2436) = mat(k,2436) + lmat(k,2436)
         mat(k,2440) = mat(k,2440) + lmat(k,2440)
         mat(k,2442) = lmat(k,2442)
         mat(k,2443) = lmat(k,2443)
         mat(k,2449) = mat(k,2449) + lmat(k,2449)
         mat(k, 212) = 0._r8
         mat(k, 213) = 0._r8
         mat(k, 252) = 0._r8
         mat(k, 293) = 0._r8
         mat(k, 347) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 450) = 0._r8
         mat(k, 496) = 0._r8
         mat(k, 498) = 0._r8
         mat(k, 513) = 0._r8
         mat(k, 636) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 644) = 0._r8
         mat(k, 647) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 671) = 0._r8
         mat(k, 683) = 0._r8
         mat(k, 684) = 0._r8
         mat(k, 688) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 693) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 766) = 0._r8
         mat(k, 768) = 0._r8
         mat(k, 769) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 793) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 799) = 0._r8
         mat(k, 800) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 818) = 0._r8
         mat(k, 823) = 0._r8
         mat(k, 848) = 0._r8
         mat(k, 849) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 862) = 0._r8
         mat(k, 867) = 0._r8
         mat(k, 868) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 870) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 939) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 975) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 992) = 0._r8
         mat(k, 995) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1006) = 0._r8
         mat(k,1019) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1045) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1052) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1088) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1092) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1096) = 0._r8
         mat(k,1098) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1107) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1136) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1306) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1321) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1329) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1333) = 0._r8
         mat(k,1334) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1339) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1377) = 0._r8
         mat(k,1380) = 0._r8
         mat(k,1381) = 0._r8
         mat(k,1383) = 0._r8
         mat(k,1384) = 0._r8
         mat(k,1386) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1391) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1430) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1436) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1441) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1447) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1489) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1501) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1516) = 0._r8
         mat(k,1518) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1520) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1543) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1545) = 0._r8
         mat(k,1547) = 0._r8
         mat(k,1548) = 0._r8
         mat(k,1552) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1571) = 0._r8
         mat(k,1572) = 0._r8
         mat(k,1573) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1577) = 0._r8
         mat(k,1579) = 0._r8
         mat(k,1593) = 0._r8
         mat(k,1598) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1601) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1627) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1635) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1648) = 0._r8
         mat(k,1650) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1657) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1662) = 0._r8
         mat(k,1663) = 0._r8
         mat(k,1665) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1731) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1800) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1876) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1936) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1961) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1965) = 0._r8
         mat(k,1966) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1981) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1996) = 0._r8
         mat(k,2026) = 0._r8
         mat(k,2029) = 0._r8
         mat(k,2032) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2044) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2051) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2054) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2059) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2062) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2080) = 0._r8
         mat(k,2084) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2098) = 0._r8
         mat(k,2103) = 0._r8
         mat(k,2105) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2112) = 0._r8
         mat(k,2113) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2116) = 0._r8
         mat(k,2117) = 0._r8
         mat(k,2118) = 0._r8
         mat(k,2119) = 0._r8
         mat(k,2122) = 0._r8
         mat(k,2124) = 0._r8
         mat(k,2126) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2144) = 0._r8
         mat(k,2145) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2161) = 0._r8
         mat(k,2162) = 0._r8
         mat(k,2164) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2166) = 0._r8
         mat(k,2173) = 0._r8
         mat(k,2188) = 0._r8
         mat(k,2194) = 0._r8
         mat(k,2195) = 0._r8
         mat(k,2196) = 0._r8
         mat(k,2199) = 0._r8
         mat(k,2204) = 0._r8
         mat(k,2205) = 0._r8
         mat(k,2206) = 0._r8
         mat(k,2208) = 0._r8
         mat(k,2211) = 0._r8
         mat(k,2212) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2229) = 0._r8
         mat(k,2237) = 0._r8
         mat(k,2244) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2247) = 0._r8
         mat(k,2248) = 0._r8
         mat(k,2250) = 0._r8
         mat(k,2251) = 0._r8
         mat(k,2254) = 0._r8
         mat(k,2258) = 0._r8
         mat(k,2293) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,2295) = 0._r8
         mat(k,2321) = 0._r8
         mat(k,2329) = 0._r8
         mat(k,2330) = 0._r8
         mat(k,2332) = 0._r8
         mat(k,2335) = 0._r8
         mat(k,2337) = 0._r8
         mat(k,2341) = 0._r8
         mat(k,2346) = 0._r8
         mat(k,2361) = 0._r8
         mat(k,2362) = 0._r8
         mat(k,2367) = 0._r8
         mat(k,2378) = 0._r8
         mat(k,2406) = 0._r8
         mat(k,2408) = 0._r8
         mat(k,2409) = 0._r8
         mat(k,2413) = 0._r8
         mat(k,2414) = 0._r8
         mat(k,2415) = 0._r8
         mat(k,2416) = 0._r8
         mat(k,2418) = 0._r8
         mat(k,2419) = 0._r8
         mat(k,2422) = 0._r8
         mat(k,2427) = 0._r8
         mat(k,2429) = 0._r8
         mat(k,2430) = 0._r8
         mat(k,2431) = 0._r8
         mat(k,2432) = 0._r8
         mat(k,2433) = 0._r8
         mat(k,2434) = 0._r8
         mat(k,2435) = 0._r8
         mat(k,2437) = 0._r8
         mat(k,2438) = 0._r8
         mat(k,2439) = 0._r8
         mat(k,2441) = 0._r8
         mat(k,2444) = 0._r8
         mat(k,2445) = 0._r8
         mat(k,2446) = 0._r8
         mat(k,2447) = 0._r8
         mat(k,2448) = 0._r8
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
         mat(k, 43) = mat(k, 43) - dti(k)
         mat(k, 44) = mat(k, 44) - dti(k)
         mat(k, 50) = mat(k, 50) - dti(k)
         mat(k, 56) = mat(k, 56) - dti(k)
         mat(k, 62) = mat(k, 62) - dti(k)
         mat(k, 68) = mat(k, 68) - dti(k)
         mat(k, 74) = mat(k, 74) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 82) = mat(k, 82) - dti(k)
         mat(k, 88) = mat(k, 88) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 98) = mat(k, 98) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 104) = mat(k, 104) - dti(k)
         mat(k, 107) = mat(k, 107) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 127) = mat(k, 127) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 160) = mat(k, 160) - dti(k)
         mat(k, 166) = mat(k, 166) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 177) = mat(k, 177) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 182) = mat(k, 182) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 194) = mat(k, 194) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 203) = mat(k, 203) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 216) = mat(k, 216) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 226) = mat(k, 226) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 235) = mat(k, 235) - dti(k)
         mat(k, 238) = mat(k, 238) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 246) = mat(k, 246) - dti(k)
         mat(k, 251) = mat(k, 251) - dti(k)
         mat(k, 256) = mat(k, 256) - dti(k)
         mat(k, 261) = mat(k, 261) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 282) = mat(k, 282) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 298) = mat(k, 298) - dti(k)
         mat(k, 304) = mat(k, 304) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 326) = mat(k, 326) - dti(k)
         mat(k, 329) = mat(k, 329) - dti(k)
         mat(k, 334) = mat(k, 334) - dti(k)
         mat(k, 339) = mat(k, 339) - dti(k)
         mat(k, 345) = mat(k, 345) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 355) = mat(k, 355) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 371) = mat(k, 371) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 397) = mat(k, 397) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
         mat(k, 427) = mat(k, 427) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 448) = mat(k, 448) - dti(k)
         mat(k, 454) = mat(k, 454) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 468) = mat(k, 468) - dti(k)
         mat(k, 472) = mat(k, 472) - dti(k)
         mat(k, 479) = mat(k, 479) - dti(k)
         mat(k, 486) = mat(k, 486) - dti(k)
         mat(k, 495) = mat(k, 495) - dti(k)
         mat(k, 502) = mat(k, 502) - dti(k)
         mat(k, 510) = mat(k, 510) - dti(k)
         mat(k, 517) = mat(k, 517) - dti(k)
         mat(k, 523) = mat(k, 523) - dti(k)
         mat(k, 529) = mat(k, 529) - dti(k)
         mat(k, 534) = mat(k, 534) - dti(k)
         mat(k, 542) = mat(k, 542) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 558) = mat(k, 558) - dti(k)
         mat(k, 566) = mat(k, 566) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 582) = mat(k, 582) - dti(k)
         mat(k, 589) = mat(k, 589) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 602) = mat(k, 602) - dti(k)
         mat(k, 611) = mat(k, 611) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 635) = mat(k, 635) - dti(k)
         mat(k, 642) = mat(k, 642) - dti(k)
         mat(k, 652) = mat(k, 652) - dti(k)
         mat(k, 665) = mat(k, 665) - dti(k)
         mat(k, 676) = mat(k, 676) - dti(k)
         mat(k, 687) = mat(k, 687) - dti(k)
         mat(k, 698) = mat(k, 698) - dti(k)
         mat(k, 705) = mat(k, 705) - dti(k)
         mat(k, 711) = mat(k, 711) - dti(k)
         mat(k, 718) = mat(k, 718) - dti(k)
         mat(k, 729) = mat(k, 729) - dti(k)
         mat(k, 740) = mat(k, 740) - dti(k)
         mat(k, 751) = mat(k, 751) - dti(k)
         mat(k, 767) = mat(k, 767) - dti(k)
         mat(k, 778) = mat(k, 778) - dti(k)
         mat(k, 787) = mat(k, 787) - dti(k)
         mat(k, 797) = mat(k, 797) - dti(k)
         mat(k, 805) = mat(k, 805) - dti(k)
         mat(k, 814) = mat(k, 814) - dti(k)
         mat(k, 824) = mat(k, 824) - dti(k)
         mat(k, 828) = mat(k, 828) - dti(k)
         mat(k, 837) = mat(k, 837) - dti(k)
         mat(k, 845) = mat(k, 845) - dti(k)
         mat(k, 853) = mat(k, 853) - dti(k)
         mat(k, 865) = mat(k, 865) - dti(k)
         mat(k, 876) = mat(k, 876) - dti(k)
         mat(k, 892) = mat(k, 892) - dti(k)
         mat(k, 901) = mat(k, 901) - dti(k)
         mat(k, 910) = mat(k, 910) - dti(k)
         mat(k, 919) = mat(k, 919) - dti(k)
         mat(k, 925) = mat(k, 925) - dti(k)
         mat(k, 935) = mat(k, 935) - dti(k)
         mat(k, 947) = mat(k, 947) - dti(k)
         mat(k, 955) = mat(k, 955) - dti(k)
         mat(k, 974) = mat(k, 974) - dti(k)
         mat(k, 999) = mat(k, 999) - dti(k)
         mat(k,1018) = mat(k,1018) - dti(k)
         mat(k,1042) = mat(k,1042) - dti(k)
         mat(k,1054) = mat(k,1054) - dti(k)
         mat(k,1062) = mat(k,1062) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1080) = mat(k,1080) - dti(k)
         mat(k,1090) = mat(k,1090) - dti(k)
         mat(k,1106) = mat(k,1106) - dti(k)
         mat(k,1119) = mat(k,1119) - dti(k)
         mat(k,1133) = mat(k,1133) - dti(k)
         mat(k,1151) = mat(k,1151) - dti(k)
         mat(k,1160) = mat(k,1160) - dti(k)
         mat(k,1166) = mat(k,1166) - dti(k)
         mat(k,1178) = mat(k,1178) - dti(k)
         mat(k,1195) = mat(k,1195) - dti(k)
         mat(k,1208) = mat(k,1208) - dti(k)
         mat(k,1217) = mat(k,1217) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1253) = mat(k,1253) - dti(k)
         mat(k,1269) = mat(k,1269) - dti(k)
         mat(k,1281) = mat(k,1281) - dti(k)
         mat(k,1299) = mat(k,1299) - dti(k)
         mat(k,1332) = mat(k,1332) - dti(k)
         mat(k,1356) = mat(k,1356) - dti(k)
         mat(k,1376) = mat(k,1376) - dti(k)
         mat(k,1397) = mat(k,1397) - dti(k)
         mat(k,1428) = mat(k,1428) - dti(k)
         mat(k,1450) = mat(k,1450) - dti(k)
         mat(k,1461) = mat(k,1461) - dti(k)
         mat(k,1476) = mat(k,1476) - dti(k)
         mat(k,1495) = mat(k,1495) - dti(k)
         mat(k,1511) = mat(k,1511) - dti(k)
         mat(k,1542) = mat(k,1542) - dti(k)
         mat(k,1565) = mat(k,1565) - dti(k)
         mat(k,1592) = mat(k,1592) - dti(k)
         mat(k,1652) = mat(k,1652) - dti(k)
         mat(k,1818) = mat(k,1818) - dti(k)
         mat(k,1914) = mat(k,1914) - dti(k)
         mat(k,1960) = mat(k,1960) - dti(k)
         mat(k,1987) = mat(k,1987) - dti(k)
         mat(k,2031) = mat(k,2031) - dti(k)
         mat(k,2055) = mat(k,2055) - dti(k)
         mat(k,2099) = mat(k,2099) - dti(k)
         mat(k,2121) = mat(k,2121) - dti(k)
         mat(k,2168) = mat(k,2168) - dti(k)
         mat(k,2233) = mat(k,2233) - dti(k)
         mat(k,2259) = mat(k,2259) - dti(k)
         mat(k,2368) = mat(k,2368) - dti(k)
         mat(k,2421) = mat(k,2421) - dti(k)
         mat(k,2449) = mat(k,2449) - dti(k)
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
