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
         mat(k,538) = -(rxt(k,396)*y(k,223))
         mat(k,1660) = -rxt(k,396)*y(k,1)
         mat(k,1823) = rxt(k,399)*y(k,192)
         mat(k,878) = rxt(k,399)*y(k,124)
         mat(k,561) = -(rxt(k,400)*y(k,223))
         mat(k,1661) = -rxt(k,400)*y(k,2)
         mat(k,879) = rxt(k,397)*y(k,205)
         mat(k,2035) = rxt(k,397)*y(k,192)
         mat(k,834) = -(rxt(k,479)*y(k,126) + rxt(k,480)*y(k,136) + rxt(k,481) &
                      *y(k,223))
         mat(k,1450) = -rxt(k,479)*y(k,6)
         mat(k,1542) = -rxt(k,480)*y(k,6)
         mat(k,1684) = -rxt(k,481)*y(k,6)
         mat(k,86) = -(rxt(k,438)*y(k,223))
         mat(k,1598) = -rxt(k,438)*y(k,7)
         mat(k,275) = -(rxt(k,441)*y(k,223))
         mat(k,1628) = -rxt(k,441)*y(k,8)
         mat(k,377) = rxt(k,439)*y(k,205)
         mat(k,2009) = rxt(k,439)*y(k,193)
         mat(k,87) = .120_r8*rxt(k,438)*y(k,223)
         mat(k,1599) = .120_r8*rxt(k,438)*y(k,7)
         mat(k,831) = .100_r8*rxt(k,480)*y(k,136)
         mat(k,857) = .100_r8*rxt(k,483)*y(k,136)
         mat(k,1532) = .100_r8*rxt(k,480)*y(k,6) + .100_r8*rxt(k,483)*y(k,110)
         mat(k,1810) = .500_r8*rxt(k,440)*y(k,193) + .200_r8*rxt(k,467)*y(k,230) &
                      + .060_r8*rxt(k,473)*y(k,232)
         mat(k,378) = .500_r8*rxt(k,440)*y(k,124)
         mat(k,616) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,632) = .060_r8*rxt(k,473)*y(k,124)
         mat(k,1804) = .200_r8*rxt(k,467)*y(k,230) + .200_r8*rxt(k,473)*y(k,232)
         mat(k,615) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,630) = .200_r8*rxt(k,473)*y(k,124)
         mat(k,1820) = .200_r8*rxt(k,467)*y(k,230) + .150_r8*rxt(k,473)*y(k,232)
         mat(k,618) = .200_r8*rxt(k,467)*y(k,124)
         mat(k,633) = .150_r8*rxt(k,473)*y(k,124)
         mat(k,1805) = .210_r8*rxt(k,473)*y(k,232)
         mat(k,631) = .210_r8*rxt(k,473)*y(k,124)
         mat(k,162) = -(rxt(k,401)*y(k,223))
         mat(k,1610) = -rxt(k,401)*y(k,15)
         mat(k,830) = .050_r8*rxt(k,480)*y(k,136)
         mat(k,856) = .050_r8*rxt(k,483)*y(k,136)
         mat(k,1531) = .050_r8*rxt(k,480)*y(k,6) + .050_r8*rxt(k,483)*y(k,110)
         mat(k,261) = -(rxt(k,367)*y(k,126) + rxt(k,368)*y(k,223))
         mat(k,1445) = -rxt(k,367)*y(k,16)
         mat(k,1626) = -rxt(k,368)*y(k,16)
         mat(k,1351) = -(rxt(k,250)*y(k,42) + rxt(k,251)*y(k,205) + rxt(k,252) &
                      *y(k,136))
         mat(k,2102) = -rxt(k,250)*y(k,17)
         mat(k,2078) = -rxt(k,251)*y(k,17)
         mat(k,1568) = -rxt(k,252)*y(k,17)
         mat(k,2127) = 4.000_r8*rxt(k,253)*y(k,19) + (rxt(k,254)+rxt(k,255))*y(k,59) &
                      + rxt(k,258)*y(k,124) + rxt(k,261)*y(k,134) + rxt(k,508) &
                      *y(k,152) + rxt(k,262)*y(k,223)
         mat(k,1506) = (rxt(k,254)+rxt(k,255))*y(k,19)
         mat(k,764) = rxt(k,263)*y(k,134) + rxt(k,269)*y(k,219) + rxt(k,264)*y(k,223)
         mat(k,1865) = rxt(k,258)*y(k,19)
         mat(k,1906) = rxt(k,261)*y(k,19) + rxt(k,263)*y(k,81)
         mat(k,1318) = rxt(k,508)*y(k,19)
         mat(k,1784) = rxt(k,269)*y(k,81)
         mat(k,1716) = rxt(k,262)*y(k,19) + rxt(k,264)*y(k,81)
         mat(k,2120) = rxt(k,256)*y(k,59)
         mat(k,1499) = rxt(k,256)*y(k,19)
         mat(k,1332) = (rxt(k,556)+rxt(k,561))*y(k,91)
         mat(k,701) = (rxt(k,556)+rxt(k,561))*y(k,85)
         mat(k,2142) = -(4._r8*rxt(k,253)*y(k,19) + (rxt(k,254) + rxt(k,255) + rxt(k,256) &
                      ) * y(k,59) + rxt(k,257)*y(k,205) + rxt(k,258)*y(k,124) &
                      + rxt(k,259)*y(k,125) + rxt(k,261)*y(k,134) + rxt(k,262) &
                      *y(k,223) + rxt(k,508)*y(k,152))
         mat(k,1522) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,19)
         mat(k,2094) = -rxt(k,257)*y(k,19)
         mat(k,1881) = -rxt(k,258)*y(k,19)
         mat(k,1999) = -rxt(k,259)*y(k,19)
         mat(k,1922) = -rxt(k,261)*y(k,19)
         mat(k,1732) = -rxt(k,262)*y(k,19)
         mat(k,1329) = -rxt(k,508)*y(k,19)
         mat(k,1360) = rxt(k,252)*y(k,136)
         mat(k,447) = rxt(k,260)*y(k,134)
         mat(k,769) = rxt(k,270)*y(k,219)
         mat(k,708) = rxt(k,265)*y(k,134)
         mat(k,1922) = mat(k,1922) + rxt(k,260)*y(k,20) + rxt(k,265)*y(k,91)
         mat(k,1584) = rxt(k,252)*y(k,17)
         mat(k,1800) = rxt(k,270)*y(k,81)
         mat(k,440) = -(rxt(k,260)*y(k,134))
         mat(k,1888) = -rxt(k,260)*y(k,20)
         mat(k,2122) = rxt(k,259)*y(k,125)
         mat(k,1965) = rxt(k,259)*y(k,19)
         mat(k,171) = -(rxt(k,442)*y(k,223))
         mat(k,1612) = -rxt(k,442)*y(k,22)
         mat(k,1803) = rxt(k,445)*y(k,194)
         mat(k,323) = rxt(k,445)*y(k,124)
         mat(k,240) = -(rxt(k,444)*y(k,223))
         mat(k,1622) = -rxt(k,444)*y(k,23)
         mat(k,324) = rxt(k,443)*y(k,205)
         mat(k,2007) = rxt(k,443)*y(k,194)
         mat(k,202) = -(rxt(k,316)*y(k,56) + rxt(k,317)*y(k,223))
         mat(k,1925) = -rxt(k,316)*y(k,24)
         mat(k,1617) = -rxt(k,317)*y(k,24)
         mat(k,452) = -(rxt(k,318)*y(k,56) + rxt(k,319)*y(k,136) + rxt(k,344)*y(k,223))
         mat(k,1927) = -rxt(k,318)*y(k,25)
         mat(k,1535) = -rxt(k,319)*y(k,25)
         mat(k,1650) = -rxt(k,344)*y(k,25)
         mat(k,179) = -(rxt(k,324)*y(k,223))
         mat(k,1614) = -rxt(k,324)*y(k,26)
         mat(k,799) = .800_r8*rxt(k,320)*y(k,195) + .200_r8*rxt(k,321)*y(k,199)
         mat(k,1362) = .200_r8*rxt(k,321)*y(k,195)
         mat(k,248) = -(rxt(k,325)*y(k,223))
         mat(k,1624) = -rxt(k,325)*y(k,27)
         mat(k,800) = rxt(k,322)*y(k,205)
         mat(k,2008) = rxt(k,322)*y(k,195)
         mat(k,208) = -(rxt(k,326)*y(k,56) + rxt(k,327)*y(k,223))
         mat(k,1926) = -rxt(k,326)*y(k,28)
         mat(k,1618) = -rxt(k,327)*y(k,28)
         mat(k,930) = -(rxt(k,347)*y(k,126) + rxt(k,348)*y(k,136) + rxt(k,365) &
                      *y(k,223))
         mat(k,1456) = -rxt(k,347)*y(k,29)
         mat(k,1548) = -rxt(k,348)*y(k,29)
         mat(k,1691) = -rxt(k,365)*y(k,29)
         mat(k,778) = .130_r8*rxt(k,425)*y(k,136)
         mat(k,1548) = mat(k,1548) + .130_r8*rxt(k,425)*y(k,98)
         mat(k,311) = -(rxt(k,352)*y(k,223))
         mat(k,1633) = -rxt(k,352)*y(k,30)
         mat(k,743) = rxt(k,350)*y(k,205)
         mat(k,2014) = rxt(k,350)*y(k,196)
         mat(k,57) = -(rxt(k,353)*y(k,223))
         mat(k,1595) = -rxt(k,353)*y(k,31)
         mat(k,183) = -(rxt(k,448)*y(k,223))
         mat(k,1615) = -rxt(k,448)*y(k,32)
         mat(k,529) = rxt(k,446)*y(k,205)
         mat(k,2003) = rxt(k,446)*y(k,197)
         mat(k,2117) = -(rxt(k,214)*y(k,56) + rxt(k,250)*y(k,17) + rxt(k,294)*y(k,205) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,134) + rxt(k,297) &
                      *y(k,223))
         mat(k,1956) = -rxt(k,214)*y(k,42)
         mat(k,1359) = -rxt(k,250)*y(k,42)
         mat(k,2093) = -rxt(k,294)*y(k,42)
         mat(k,1494) = -rxt(k,295)*y(k,42)
         mat(k,1921) = -rxt(k,296)*y(k,42)
         mat(k,1731) = -rxt(k,297)*y(k,42)
         mat(k,547) = .400_r8*rxt(k,396)*y(k,223)
         mat(k,849) = .340_r8*rxt(k,480)*y(k,136)
         mat(k,268) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,459) = rxt(k,319)*y(k,136)
         mat(k,944) = .500_r8*rxt(k,348)*y(k,136)
         mat(k,410) = .500_r8*rxt(k,336)*y(k,223)
         mat(k,713) = rxt(k,302)*y(k,223)
         mat(k,285) = .300_r8*rxt(k,303)*y(k,223)
         mat(k,1521) = rxt(k,221)*y(k,199)
         mat(k,966) = .800_r8*rxt(k,341)*y(k,223)
         mat(k,791) = .910_r8*rxt(k,425)*y(k,136)
         mat(k,484) = .300_r8*rxt(k,416)*y(k,223)
         mat(k,1139) = .800_r8*rxt(k,420)*y(k,199)
         mat(k,1151) = .120_r8*rxt(k,378)*y(k,136)
         mat(k,467) = .500_r8*rxt(k,391)*y(k,223)
         mat(k,875) = .340_r8*rxt(k,483)*y(k,136)
         mat(k,1263) = .600_r8*rxt(k,392)*y(k,136)
         mat(k,1880) = .100_r8*rxt(k,398)*y(k,192) + rxt(k,301)*y(k,199) &
                      + .500_r8*rxt(k,369)*y(k,202) + .500_r8*rxt(k,338)*y(k,204) &
                      + .920_r8*rxt(k,408)*y(k,207) + .250_r8*rxt(k,376)*y(k,209) &
                      + rxt(k,385)*y(k,211) + rxt(k,359)*y(k,226) + rxt(k,363) &
                      *y(k,227) + .340_r8*rxt(k,492)*y(k,228) + .320_r8*rxt(k,497) &
                      *y(k,229) + .250_r8*rxt(k,433)*y(k,231)
         mat(k,1494) = mat(k,1494) + .500_r8*rxt(k,367)*y(k,16) + rxt(k,409)*y(k,207) &
                      + .250_r8*rxt(k,375)*y(k,209) + rxt(k,386)*y(k,211)
         mat(k,1583) = .340_r8*rxt(k,480)*y(k,6) + rxt(k,319)*y(k,25) &
                      + .500_r8*rxt(k,348)*y(k,29) + .910_r8*rxt(k,425)*y(k,98) &
                      + .120_r8*rxt(k,378)*y(k,105) + .340_r8*rxt(k,483)*y(k,110) &
                      + .600_r8*rxt(k,392)*y(k,111)
         mat(k,358) = rxt(k,343)*y(k,223)
         mat(k,975) = .680_r8*rxt(k,501)*y(k,223)
         mat(k,892) = .100_r8*rxt(k,398)*y(k,124)
         mat(k,810) = .700_r8*rxt(k,321)*y(k,199)
         mat(k,753) = rxt(k,349)*y(k,199)
         mat(k,1312) = rxt(k,332)*y(k,199) + rxt(k,405)*y(k,207) + .250_r8*rxt(k,372) &
                      *y(k,209) + rxt(k,381)*y(k,211) + .250_r8*rxt(k,430)*y(k,231)
         mat(k,1407) = rxt(k,221)*y(k,59) + .800_r8*rxt(k,420)*y(k,101) + rxt(k,301) &
                      *y(k,124) + .700_r8*rxt(k,321)*y(k,195) + rxt(k,349)*y(k,196) &
                      + rxt(k,332)*y(k,198) + (4.000_r8*rxt(k,298)+2.000_r8*rxt(k,299)) &
                      *y(k,199) + 1.500_r8*rxt(k,406)*y(k,207) + .750_r8*rxt(k,411) &
                      *y(k,208) + .880_r8*rxt(k,373)*y(k,209) + 2.000_r8*rxt(k,382) &
                      *y(k,211) + .750_r8*rxt(k,485)*y(k,218) + .800_r8*rxt(k,361) &
                      *y(k,227) + .930_r8*rxt(k,490)*y(k,228) + .950_r8*rxt(k,495) &
                      *y(k,229) + .800_r8*rxt(k,431)*y(k,231)
         mat(k,475) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,673) = .500_r8*rxt(k,338)*y(k,124)
         mat(k,2093) = mat(k,2093) + .450_r8*rxt(k,383)*y(k,211) + .150_r8*rxt(k,362) &
                      *y(k,227)
         mat(k,1220) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) + rxt(k,405) &
                      *y(k,198) + 1.500_r8*rxt(k,406)*y(k,199)
         mat(k,1193) = .750_r8*rxt(k,411)*y(k,199)
         mat(k,1241) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,198) + .880_r8*rxt(k,373)*y(k,199)
         mat(k,1281) = rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126) + rxt(k,381)*y(k,198) &
                      + 2.000_r8*rxt(k,382)*y(k,199) + .450_r8*rxt(k,383)*y(k,205) &
                      + 4.000_r8*rxt(k,384)*y(k,211)
         mat(k,1072) = .750_r8*rxt(k,485)*y(k,199)
         mat(k,1731) = mat(k,1731) + .400_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,336) &
                      *y(k,51) + rxt(k,302)*y(k,52) + .300_r8*rxt(k,303)*y(k,53) &
                      + .800_r8*rxt(k,341)*y(k,74) + .300_r8*rxt(k,416)*y(k,99) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,343)*y(k,141) &
                      + .680_r8*rxt(k,501)*y(k,181)
         mat(k,727) = rxt(k,359)*y(k,124)
         mat(k,1085) = rxt(k,363)*y(k,124) + .800_r8*rxt(k,361)*y(k,199) &
                      + .150_r8*rxt(k,362)*y(k,205)
         mat(k,1052) = .340_r8*rxt(k,492)*y(k,124) + .930_r8*rxt(k,490)*y(k,199)
         mat(k,1033) = .320_r8*rxt(k,497)*y(k,124) + .950_r8*rxt(k,495)*y(k,199)
         mat(k,1116) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,430)*y(k,198) &
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
         mat(k,977) = -(rxt(k,328)*y(k,126) + rxt(k,329)*y(k,223))
         mat(k,1460) = -rxt(k,328)*y(k,45)
         mat(k,1695) = -rxt(k,329)*y(k,45)
         mat(k,542) = .800_r8*rxt(k,396)*y(k,223)
         mat(k,264) = rxt(k,367)*y(k,126)
         mat(k,180) = rxt(k,324)*y(k,223)
         mat(k,250) = .500_r8*rxt(k,325)*y(k,223)
         mat(k,931) = .500_r8*rxt(k,348)*y(k,136)
         mat(k,1245) = .100_r8*rxt(k,392)*y(k,136)
         mat(k,1847) = .400_r8*rxt(k,398)*y(k,192) + rxt(k,323)*y(k,195) &
                      + .270_r8*rxt(k,351)*y(k,196) + rxt(k,369)*y(k,202) + rxt(k,388) &
                      *y(k,213) + rxt(k,359)*y(k,226)
         mat(k,1460) = mat(k,1460) + rxt(k,367)*y(k,16)
         mat(k,1551) = .500_r8*rxt(k,348)*y(k,29) + .100_r8*rxt(k,392)*y(k,111)
         mat(k,884) = .400_r8*rxt(k,398)*y(k,124)
         mat(k,803) = rxt(k,323)*y(k,124) + 3.200_r8*rxt(k,320)*y(k,195) &
                      + .800_r8*rxt(k,321)*y(k,199)
         mat(k,746) = .270_r8*rxt(k,351)*y(k,124)
         mat(k,1378) = .800_r8*rxt(k,321)*y(k,195)
         mat(k,470) = rxt(k,369)*y(k,124)
         mat(k,2059) = .200_r8*rxt(k,387)*y(k,213)
         mat(k,573) = rxt(k,388)*y(k,124) + .200_r8*rxt(k,387)*y(k,205)
         mat(k,1695) = mat(k,1695) + .800_r8*rxt(k,396)*y(k,1) + rxt(k,324)*y(k,26) &
                      + .500_r8*rxt(k,325)*y(k,27)
         mat(k,720) = rxt(k,359)*y(k,124)
         mat(k,54) = -(rxt(k,330)*y(k,223))
         mat(k,1594) = -rxt(k,330)*y(k,47)
         mat(k,900) = -(rxt(k,366)*y(k,223))
         mat(k,1688) = -rxt(k,366)*y(k,48)
         mat(k,541) = .800_r8*rxt(k,396)*y(k,223)
         mat(k,836) = .520_r8*rxt(k,480)*y(k,136)
         mat(k,263) = .500_r8*rxt(k,367)*y(k,126)
         mat(k,862) = .520_r8*rxt(k,483)*y(k,136)
         mat(k,1842) = .250_r8*rxt(k,398)*y(k,192) + .820_r8*rxt(k,351)*y(k,196) &
                      + .500_r8*rxt(k,369)*y(k,202) + .270_r8*rxt(k,492)*y(k,228) &
                      + .040_r8*rxt(k,497)*y(k,229)
         mat(k,1454) = .500_r8*rxt(k,367)*y(k,16)
         mat(k,1546) = .520_r8*rxt(k,480)*y(k,6) + .520_r8*rxt(k,483)*y(k,110)
         mat(k,967) = .500_r8*rxt(k,501)*y(k,223)
         mat(k,883) = .250_r8*rxt(k,398)*y(k,124)
         mat(k,745) = .820_r8*rxt(k,351)*y(k,124) + .820_r8*rxt(k,349)*y(k,199)
         mat(k,1373) = .820_r8*rxt(k,349)*y(k,196) + .150_r8*rxt(k,490)*y(k,228) &
                      + .025_r8*rxt(k,495)*y(k,229)
         mat(k,469) = .500_r8*rxt(k,369)*y(k,124)
         mat(k,1688) = mat(k,1688) + .800_r8*rxt(k,396)*y(k,1) + .500_r8*rxt(k,501) &
                      *y(k,181)
         mat(k,1038) = .270_r8*rxt(k,492)*y(k,124) + .150_r8*rxt(k,490)*y(k,199)
         mat(k,1016) = .040_r8*rxt(k,497)*y(k,124) + .025_r8*rxt(k,495)*y(k,199)
         mat(k,1154) = -(rxt(k,354)*y(k,126) + rxt(k,355)*y(k,223))
         mat(k,1471) = -rxt(k,354)*y(k,49)
         mat(k,1707) = -rxt(k,355)*y(k,49)
         mat(k,1008) = rxt(k,356)*y(k,223)
         mat(k,1143) = .880_r8*rxt(k,378)*y(k,136)
         mat(k,1248) = .500_r8*rxt(k,392)*y(k,136)
         mat(k,1858) = .170_r8*rxt(k,451)*y(k,200) + .050_r8*rxt(k,414)*y(k,208) &
                      + .250_r8*rxt(k,376)*y(k,209) + .170_r8*rxt(k,457)*y(k,212) &
                      + .400_r8*rxt(k,467)*y(k,230) + .250_r8*rxt(k,433)*y(k,231) &
                      + .540_r8*rxt(k,473)*y(k,232) + .510_r8*rxt(k,476)*y(k,233)
         mat(k,1471) = mat(k,1471) + .050_r8*rxt(k,415)*y(k,208) + .250_r8*rxt(k,375) &
                      *y(k,209) + .250_r8*rxt(k,434)*y(k,231)
         mat(k,794) = rxt(k,357)*y(k,223)
         mat(k,1560) = .880_r8*rxt(k,378)*y(k,105) + .500_r8*rxt(k,392)*y(k,111)
         mat(k,1295) = .250_r8*rxt(k,372)*y(k,209) + .250_r8*rxt(k,430)*y(k,231)
         mat(k,1388) = .240_r8*rxt(k,373)*y(k,209) + .500_r8*rxt(k,361)*y(k,227) &
                      + .100_r8*rxt(k,431)*y(k,231)
         mat(k,649) = .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450)*y(k,205)
         mat(k,2070) = .070_r8*rxt(k,450)*y(k,200) + .070_r8*rxt(k,456)*y(k,212)
         mat(k,1178) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1229) = .250_r8*rxt(k,376)*y(k,124) + .250_r8*rxt(k,375)*y(k,126) &
                      + .250_r8*rxt(k,372)*y(k,198) + .240_r8*rxt(k,373)*y(k,199)
         mat(k,819) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,205)
         mat(k,1707) = mat(k,1707) + rxt(k,356)*y(k,95) + rxt(k,357)*y(k,127)
         mat(k,1078) = .500_r8*rxt(k,361)*y(k,199)
         mat(k,625) = .400_r8*rxt(k,467)*y(k,124)
         mat(k,1107) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,641) = .540_r8*rxt(k,473)*y(k,124)
         mat(k,396) = .510_r8*rxt(k,476)*y(k,124)
         mat(k,448) = -(rxt(k,335)*y(k,223))
         mat(k,1649) = -rxt(k,335)*y(k,50)
         mat(k,926) = .120_r8*rxt(k,348)*y(k,136)
         mat(k,1534) = .120_r8*rxt(k,348)*y(k,29)
         mat(k,1286) = .100_r8*rxt(k,332)*y(k,199) + .150_r8*rxt(k,333)*y(k,205)
         mat(k,1366) = .100_r8*rxt(k,332)*y(k,198)
         mat(k,2029) = .150_r8*rxt(k,333)*y(k,198) + .150_r8*rxt(k,383)*y(k,211)
         mat(k,1266) = .150_r8*rxt(k,383)*y(k,205)
         mat(k,406) = -(rxt(k,336)*y(k,223))
         mat(k,1645) = -rxt(k,336)*y(k,51)
         mat(k,1285) = .400_r8*rxt(k,333)*y(k,205)
         mat(k,2026) = .400_r8*rxt(k,333)*y(k,198) + .400_r8*rxt(k,383)*y(k,211)
         mat(k,1265) = .400_r8*rxt(k,383)*y(k,205)
         mat(k,710) = -(rxt(k,302)*y(k,223))
         mat(k,1672) = -rxt(k,302)*y(k,52)
         mat(k,1120) = .200_r8*rxt(k,420)*y(k,199)
         mat(k,801) = .300_r8*rxt(k,321)*y(k,199)
         mat(k,1368) = .200_r8*rxt(k,420)*y(k,101) + .300_r8*rxt(k,321)*y(k,195) &
                      + 2.000_r8*rxt(k,299)*y(k,199) + .250_r8*rxt(k,406)*y(k,207) &
                      + .250_r8*rxt(k,411)*y(k,208) + .250_r8*rxt(k,373)*y(k,209) &
                      + .250_r8*rxt(k,485)*y(k,218) + .500_r8*rxt(k,361)*y(k,227) &
                      + .250_r8*rxt(k,490)*y(k,228) + .250_r8*rxt(k,495)*y(k,229) &
                      + .300_r8*rxt(k,431)*y(k,231)
         mat(k,1197) = .250_r8*rxt(k,406)*y(k,199)
         mat(k,1167) = .250_r8*rxt(k,411)*y(k,199)
         mat(k,1223) = .250_r8*rxt(k,373)*y(k,199)
         mat(k,1056) = .250_r8*rxt(k,485)*y(k,199)
         mat(k,1075) = .500_r8*rxt(k,361)*y(k,199)
         mat(k,1037) = .250_r8*rxt(k,490)*y(k,199)
         mat(k,1015) = .250_r8*rxt(k,495)*y(k,199)
         mat(k,1101) = .300_r8*rxt(k,431)*y(k,199)
         mat(k,281) = -(rxt(k,303)*y(k,223))
         mat(k,1629) = -rxt(k,303)*y(k,53)
         mat(k,1365) = rxt(k,300)*y(k,205)
         mat(k,2010) = rxt(k,300)*y(k,199)
         mat(k,1953) = -(rxt(k,214)*y(k,42) + rxt(k,216)*y(k,77) + rxt(k,217)*y(k,79) &
                      + (rxt(k,218) + rxt(k,219)) * y(k,205) + rxt(k,220)*y(k,136) &
                      + rxt(k,227)*y(k,60) + rxt(k,236)*y(k,92) + rxt(k,326)*y(k,28))
         mat(k,2114) = -rxt(k,214)*y(k,56)
         mat(k,1098) = -rxt(k,216)*y(k,56)
         mat(k,509) = -rxt(k,217)*y(k,56)
         mat(k,2090) = -(rxt(k,218) + rxt(k,219)) * y(k,56)
         mat(k,1580) = -rxt(k,220)*y(k,56)
         mat(k,915) = -rxt(k,227)*y(k,56)
         mat(k,761) = -rxt(k,236)*y(k,56)
         mat(k,212) = -rxt(k,326)*y(k,56)
         mat(k,2138) = rxt(k,255)*y(k,59)
         mat(k,1518) = rxt(k,255)*y(k,19) + (4.000_r8*rxt(k,222)+2.000_r8*rxt(k,224)) &
                      *y(k,59) + rxt(k,226)*y(k,124) + rxt(k,231)*y(k,134) &
                      + rxt(k,509)*y(k,152) + rxt(k,221)*y(k,199) + rxt(k,232) &
                      *y(k,223)
         mat(k,134) = rxt(k,276)*y(k,219)
         mat(k,1345) = rxt(k,234)*y(k,134) + rxt(k,246)*y(k,219) + rxt(k,235)*y(k,223)
         mat(k,1877) = rxt(k,226)*y(k,59)
         mat(k,1918) = rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85)
         mat(k,1326) = rxt(k,509)*y(k,59)
         mat(k,1404) = rxt(k,221)*y(k,59)
         mat(k,1796) = rxt(k,276)*y(k,65) + rxt(k,246)*y(k,85)
         mat(k,1728) = rxt(k,232)*y(k,59) + rxt(k,235)*y(k,85)
         mat(k,1924) = rxt(k,227)*y(k,60)
         mat(k,1498) = 2.000_r8*rxt(k,223)*y(k,59)
         mat(k,906) = rxt(k,227)*y(k,56) + (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,85)
         mat(k,1331) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,60) + (rxt(k,549) &
                       +rxt(k,555)+rxt(k,560))*y(k,92)
         mat(k,755) = (rxt(k,549)+rxt(k,555)+rxt(k,560))*y(k,85)
         mat(k,1497) = 2.000_r8*rxt(k,248)*y(k,59)
         mat(k,1510) = -(rxt(k,221)*y(k,199) + (4._r8*rxt(k,222) + 4._r8*rxt(k,223) &
                      + 4._r8*rxt(k,224) + 4._r8*rxt(k,248)) * y(k,59) + rxt(k,225) &
                      *y(k,205) + rxt(k,226)*y(k,124) + rxt(k,228)*y(k,125) + rxt(k,231) &
                      *y(k,134) + (rxt(k,232) + rxt(k,233)) * y(k,223) + (rxt(k,254) &
                      + rxt(k,255) + rxt(k,256)) * y(k,19) + rxt(k,509)*y(k,152))
         mat(k,1398) = -rxt(k,221)*y(k,59)
         mat(k,2082) = -rxt(k,225)*y(k,59)
         mat(k,1869) = -rxt(k,226)*y(k,59)
         mat(k,1987) = -rxt(k,228)*y(k,59)
         mat(k,1910) = -rxt(k,231)*y(k,59)
         mat(k,1720) = -(rxt(k,232) + rxt(k,233)) * y(k,59)
         mat(k,2130) = -(rxt(k,254) + rxt(k,255) + rxt(k,256)) * y(k,59)
         mat(k,1320) = -rxt(k,509)*y(k,59)
         mat(k,1945) = rxt(k,236)*y(k,92) + rxt(k,220)*y(k,136) + rxt(k,219)*y(k,205)
         mat(k,911) = rxt(k,229)*y(k,134)
         mat(k,1339) = rxt(k,247)*y(k,219)
         mat(k,758) = rxt(k,236)*y(k,56) + rxt(k,237)*y(k,134) + rxt(k,238)*y(k,223)
         mat(k,1910) = mat(k,1910) + rxt(k,229)*y(k,60) + rxt(k,237)*y(k,92)
         mat(k,1572) = rxt(k,220)*y(k,56)
         mat(k,232) = rxt(k,514)*y(k,152)
         mat(k,1320) = mat(k,1320) + rxt(k,514)*y(k,138)
         mat(k,2082) = mat(k,2082) + rxt(k,219)*y(k,56)
         mat(k,1788) = rxt(k,247)*y(k,85)
         mat(k,1720) = mat(k,1720) + rxt(k,238)*y(k,92)
         mat(k,908) = -(rxt(k,227)*y(k,56) + rxt(k,229)*y(k,134) + rxt(k,230)*y(k,223) &
                      + (rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,85))
         mat(k,1934) = -rxt(k,227)*y(k,60)
         mat(k,1901) = -rxt(k,229)*y(k,60)
         mat(k,1689) = -rxt(k,230)*y(k,60)
         mat(k,1335) = -(rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,60)
         mat(k,1503) = rxt(k,228)*y(k,125)
         mat(k,1975) = rxt(k,228)*y(k,59)
         mat(k,1003) = -((rxt(k,305) + rxt(k,315)) * y(k,223))
         mat(k,1697) = -(rxt(k,305) + rxt(k,315)) * y(k,62)
         mat(k,839) = .230_r8*rxt(k,480)*y(k,136)
         mat(k,1350) = rxt(k,250)*y(k,42)
         mat(k,205) = .350_r8*rxt(k,317)*y(k,223)
         mat(k,455) = .630_r8*rxt(k,319)*y(k,136)
         mat(k,932) = .560_r8*rxt(k,348)*y(k,136)
         mat(k,2099) = rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) + rxt(k,295)*y(k,126) &
                      + rxt(k,296)*y(k,134) + rxt(k,297)*y(k,223)
         mat(k,1153) = rxt(k,354)*y(k,126) + rxt(k,355)*y(k,223)
         mat(k,1937) = rxt(k,214)*y(k,42)
         mat(k,813) = rxt(k,342)*y(k,223)
         mat(k,779) = .620_r8*rxt(k,425)*y(k,136)
         mat(k,1141) = .650_r8*rxt(k,378)*y(k,136)
         mat(k,865) = .230_r8*rxt(k,483)*y(k,136)
         mat(k,1246) = .560_r8*rxt(k,392)*y(k,136)
         mat(k,1849) = .170_r8*rxt(k,451)*y(k,200) + .220_r8*rxt(k,376)*y(k,209) &
                      + .400_r8*rxt(k,454)*y(k,210) + .350_r8*rxt(k,457)*y(k,212) &
                      + .225_r8*rxt(k,492)*y(k,228) + .250_r8*rxt(k,433)*y(k,231)
         mat(k,1462) = rxt(k,295)*y(k,42) + rxt(k,354)*y(k,49) + .220_r8*rxt(k,375) &
                      *y(k,209) + .500_r8*rxt(k,434)*y(k,231)
         mat(k,1902) = rxt(k,296)*y(k,42) + rxt(k,504)*y(k,139)
         mat(k,1552) = .230_r8*rxt(k,480)*y(k,6) + .630_r8*rxt(k,319)*y(k,25) &
                      + .560_r8*rxt(k,348)*y(k,29) + .620_r8*rxt(k,425)*y(k,98) &
                      + .650_r8*rxt(k,378)*y(k,105) + .230_r8*rxt(k,483)*y(k,110) &
                      + .560_r8*rxt(k,392)*y(k,111)
         mat(k,256) = rxt(k,504)*y(k,134) + rxt(k,505)*y(k,223)
         mat(k,969) = .700_r8*rxt(k,501)*y(k,223)
         mat(k,1290) = .220_r8*rxt(k,372)*y(k,209) + .250_r8*rxt(k,430)*y(k,231)
         mat(k,1379) = .110_r8*rxt(k,373)*y(k,209) + .125_r8*rxt(k,490)*y(k,228) &
                      + .200_r8*rxt(k,431)*y(k,231)
         mat(k,648) = .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450)*y(k,205)
         mat(k,2060) = .070_r8*rxt(k,450)*y(k,200) + .160_r8*rxt(k,453)*y(k,210) &
                      + .140_r8*rxt(k,456)*y(k,212)
         mat(k,1226) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,198) + .110_r8*rxt(k,373)*y(k,199)
         mat(k,611) = .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453)*y(k,205)
         mat(k,818) = .350_r8*rxt(k,457)*y(k,124) + .140_r8*rxt(k,456)*y(k,205)
         mat(k,1697) = mat(k,1697) + .350_r8*rxt(k,317)*y(k,24) + rxt(k,297)*y(k,42) &
                      + rxt(k,355)*y(k,49) + rxt(k,342)*y(k,75) + rxt(k,505)*y(k,139) &
                      + .700_r8*rxt(k,501)*y(k,181)
         mat(k,1041) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,199)
         mat(k,1104) = .250_r8*rxt(k,433)*y(k,124) + .500_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .200_r8*rxt(k,431)*y(k,199)
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
         mat(k,61) = -(rxt(k,275)*y(k,219))
         mat(k,1776) = -rxt(k,275)*y(k,64)
         mat(k,131) = -(rxt(k,276)*y(k,219))
         mat(k,1779) = -rxt(k,276)*y(k,65)
         mat(k,119) = -(rxt(k,449)*y(k,223))
         mat(k,1603) = -rxt(k,449)*y(k,66)
         mat(k,113) = .180_r8*rxt(k,469)*y(k,223)
         mat(k,1603) = mat(k,1603) + .180_r8*rxt(k,469)*y(k,183)
         mat(k,193) = -(rxt(k,502)*y(k,126) + (rxt(k,503) + rxt(k,516)) * y(k,223))
         mat(k,1443) = -rxt(k,502)*y(k,67)
         mat(k,1616) = -(rxt(k,503) + rxt(k,516)) * y(k,67)
         mat(k,664) = rxt(k,337)*y(k,205)
         mat(k,2001) = rxt(k,337)*y(k,204)
         mat(k,656) = -(rxt(k,272)*y(k,77) + rxt(k,273)*y(k,234) + rxt(k,274)*y(k,89))
         mat(k,1088) = -rxt(k,272)*y(k,73)
         mat(k,2147) = -rxt(k,273)*y(k,73)
         mat(k,1735) = -rxt(k,274)*y(k,73)
         mat(k,62) = 2.000_r8*rxt(k,275)*y(k,219)
         mat(k,132) = rxt(k,276)*y(k,219)
         mat(k,1780) = 2.000_r8*rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65)
         mat(k,962) = -(rxt(k,341)*y(k,223))
         mat(k,1693) = -rxt(k,341)*y(k,74)
         mat(k,477) = .700_r8*rxt(k,416)*y(k,223)
         mat(k,426) = .500_r8*rxt(k,417)*y(k,223)
         mat(k,271) = rxt(k,428)*y(k,223)
         mat(k,1845) = .050_r8*rxt(k,414)*y(k,208) + .530_r8*rxt(k,376)*y(k,209) &
                      + .225_r8*rxt(k,492)*y(k,228) + .250_r8*rxt(k,433)*y(k,231)
         mat(k,1458) = .050_r8*rxt(k,415)*y(k,208) + .530_r8*rxt(k,375)*y(k,209) &
                      + .250_r8*rxt(k,434)*y(k,231)
         mat(k,1422) = rxt(k,340)*y(k,203)
         mat(k,1289) = .530_r8*rxt(k,372)*y(k,209) + .250_r8*rxt(k,430)*y(k,231)
         mat(k,1376) = .260_r8*rxt(k,373)*y(k,209) + .125_r8*rxt(k,490)*y(k,228) &
                      + .100_r8*rxt(k,431)*y(k,231)
         mat(k,348) = rxt(k,340)*y(k,135)
         mat(k,1171) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1224) = .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375)*y(k,126) &
                      + .530_r8*rxt(k,372)*y(k,198) + .260_r8*rxt(k,373)*y(k,199)
         mat(k,1693) = mat(k,1693) + .700_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
                      *y(k,100) + rxt(k,428)*y(k,115)
         mat(k,1039) = .225_r8*rxt(k,492)*y(k,124) + .125_r8*rxt(k,490)*y(k,199)
         mat(k,1103) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,812) = -(rxt(k,342)*y(k,223))
         mat(k,1682) = -rxt(k,342)*y(k,75)
         mat(k,204) = .650_r8*rxt(k,317)*y(k,223)
         mat(k,961) = .200_r8*rxt(k,341)*y(k,223)
         mat(k,948) = rxt(k,429)*y(k,223)
         mat(k,1839) = rxt(k,440)*y(k,193) + .050_r8*rxt(k,414)*y(k,208) &
                      + .400_r8*rxt(k,454)*y(k,210) + .170_r8*rxt(k,457)*y(k,212) &
                      + .700_r8*rxt(k,460)*y(k,225) + .600_r8*rxt(k,467)*y(k,230) &
                      + .250_r8*rxt(k,433)*y(k,231) + .340_r8*rxt(k,473)*y(k,232) &
                      + .170_r8*rxt(k,476)*y(k,233)
         mat(k,1449) = .050_r8*rxt(k,415)*y(k,208) + .250_r8*rxt(k,434)*y(k,231)
         mat(k,381) = rxt(k,440)*y(k,124)
         mat(k,1287) = .250_r8*rxt(k,430)*y(k,231)
         mat(k,1372) = .100_r8*rxt(k,431)*y(k,231)
         mat(k,2052) = .160_r8*rxt(k,453)*y(k,210) + .070_r8*rxt(k,456)*y(k,212)
         mat(k,1169) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,610) = .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453)*y(k,205)
         mat(k,816) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,205)
         mat(k,1682) = mat(k,1682) + .650_r8*rxt(k,317)*y(k,24) + .200_r8*rxt(k,341) &
                      *y(k,74) + rxt(k,429)*y(k,116)
         mat(k,339) = .700_r8*rxt(k,460)*y(k,124)
         mat(k,622) = .600_r8*rxt(k,467)*y(k,124)
         mat(k,1102) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,638) = .340_r8*rxt(k,473)*y(k,124)
         mat(k,395) = .170_r8*rxt(k,476)*y(k,124)
         mat(k,1766) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,205) + rxt(k,175) &
                      *y(k,135) + rxt(k,178)*y(k,136))
         mat(k,2086) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76)
         mat(k,1430) = -rxt(k,175)*y(k,76)
         mat(k,1576) = -rxt(k,178)*y(k,76)
         mat(k,2110) = rxt(k,297)*y(k,223)
         mat(k,1949) = rxt(k,216)*y(k,77)
         mat(k,1005) = rxt(k,315)*y(k,223)
         mat(k,662) = rxt(k,272)*y(k,77)
         mat(k,1095) = rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73) + rxt(k,170)*y(k,134) &
                      + rxt(k,153)*y(k,219) + rxt(k,179)*y(k,223)
         mat(k,766) = rxt(k,270)*y(k,219)
         mat(k,1342) = rxt(k,247)*y(k,219)
         mat(k,738) = rxt(k,202)*y(k,223)
         mat(k,1914) = rxt(k,170)*y(k,77) + rxt(k,182)*y(k,223)
         mat(k,259) = rxt(k,505)*y(k,223)
         mat(k,592) = rxt(k,510)*y(k,223)
         mat(k,1323) = rxt(k,515)*y(k,223)
         mat(k,1792) = rxt(k,153)*y(k,77) + rxt(k,270)*y(k,81) + rxt(k,247)*y(k,85)
         mat(k,1724) = rxt(k,297)*y(k,42) + rxt(k,315)*y(k,62) + rxt(k,179)*y(k,77) &
                      + rxt(k,202)*y(k,112) + rxt(k,182)*y(k,134) + rxt(k,505) &
                      *y(k,139) + rxt(k,510)*y(k,150) + rxt(k,515)*y(k,152)
         mat(k,1089) = -(rxt(k,153)*y(k,219) + rxt(k,170)*y(k,134) + rxt(k,179) &
                      *y(k,223) + rxt(k,216)*y(k,56) + rxt(k,272)*y(k,73))
         mat(k,1782) = -rxt(k,153)*y(k,77)
         mat(k,1903) = -rxt(k,170)*y(k,77)
         mat(k,1703) = -rxt(k,179)*y(k,77)
         mat(k,1938) = -rxt(k,216)*y(k,77)
         mat(k,657) = -rxt(k,272)*y(k,77)
         mat(k,1756) = rxt(k,172)*y(k,205)
         mat(k,2066) = rxt(k,172)*y(k,76)
         mat(k,505) = -(rxt(k,171)*y(k,134) + rxt(k,180)*y(k,223) + rxt(k,217)*y(k,56))
         mat(k,1889) = -rxt(k,171)*y(k,79)
         mat(k,1656) = -rxt(k,180)*y(k,79)
         mat(k,1928) = -rxt(k,217)*y(k,79)
         mat(k,2032) = 2.000_r8*rxt(k,186)*y(k,205)
         mat(k,1656) = mat(k,1656) + 2.000_r8*rxt(k,185)*y(k,223)
         mat(k,174) = rxt(k,518)*y(k,234)
         mat(k,2144) = rxt(k,518)*y(k,154)
         mat(k,763) = -(rxt(k,263)*y(k,134) + rxt(k,264)*y(k,223) + (rxt(k,269) &
                      + rxt(k,270)) * y(k,219))
         mat(k,1899) = -rxt(k,263)*y(k,81)
         mat(k,1678) = -rxt(k,264)*y(k,81)
         mat(k,1781) = -(rxt(k,269) + rxt(k,270)) * y(k,81)
         mat(k,1349) = rxt(k,250)*y(k,42) + rxt(k,251)*y(k,205)
         mat(k,2098) = rxt(k,250)*y(k,17)
         mat(k,2049) = rxt(k,251)*y(k,17)
         mat(k,1336) = -(rxt(k,234)*y(k,134) + rxt(k,235)*y(k,223) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,219) + (rxt(k,549) + rxt(k,555) + rxt(k,560) &
                      ) * y(k,92) + (rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,60) &
                      + (rxt(k,556) + rxt(k,561)) * y(k,91))
         mat(k,1905) = -rxt(k,234)*y(k,85)
         mat(k,1715) = -rxt(k,235)*y(k,85)
         mat(k,1783) = -(rxt(k,246) + rxt(k,247)) * y(k,85)
         mat(k,757) = -(rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,85)
         mat(k,909) = -(rxt(k,554) + rxt(k,559) + rxt(k,564)) * y(k,85)
         mat(k,703) = -(rxt(k,556) + rxt(k,561)) * y(k,85)
         mat(k,210) = rxt(k,326)*y(k,56)
         mat(k,2101) = rxt(k,214)*y(k,56)
         mat(k,1940) = rxt(k,326)*y(k,28) + rxt(k,214)*y(k,42) + rxt(k,216)*y(k,77) &
                      + rxt(k,217)*y(k,79) + rxt(k,236)*y(k,92) + rxt(k,218)*y(k,205)
         mat(k,1505) = rxt(k,233)*y(k,223)
         mat(k,1090) = rxt(k,216)*y(k,56)
         mat(k,506) = rxt(k,217)*y(k,56)
         mat(k,757) = mat(k,757) + rxt(k,236)*y(k,56)
         mat(k,2077) = rxt(k,218)*y(k,56)
         mat(k,1715) = mat(k,1715) + rxt(k,233)*y(k,59)
         mat(k,103) = -(rxt(k,306)*y(k,223) + rxt(k,314)*y(k,219))
         mat(k,1601) = -rxt(k,306)*y(k,86)
         mat(k,1778) = -rxt(k,314)*y(k,86)
         mat(k,714) = -(rxt(k,307)*y(k,223))
         mat(k,1673) = -rxt(k,307)*y(k,87)
         mat(k,832) = .050_r8*rxt(k,480)*y(k,136)
         mat(k,203) = .350_r8*rxt(k,317)*y(k,223)
         mat(k,454) = .370_r8*rxt(k,319)*y(k,136)
         mat(k,927) = .120_r8*rxt(k,348)*y(k,136)
         mat(k,776) = .110_r8*rxt(k,425)*y(k,136)
         mat(k,1140) = .330_r8*rxt(k,378)*y(k,136)
         mat(k,858) = .050_r8*rxt(k,483)*y(k,136)
         mat(k,1243) = .120_r8*rxt(k,392)*y(k,136)
         mat(k,1833) = rxt(k,310)*y(k,206)
         mat(k,1539) = .050_r8*rxt(k,480)*y(k,6) + .370_r8*rxt(k,319)*y(k,25) &
                      + .120_r8*rxt(k,348)*y(k,29) + .110_r8*rxt(k,425)*y(k,98) &
                      + .330_r8*rxt(k,378)*y(k,105) + .050_r8*rxt(k,483)*y(k,110) &
                      + .120_r8*rxt(k,392)*y(k,111)
         mat(k,2045) = rxt(k,308)*y(k,206)
         mat(k,332) = rxt(k,310)*y(k,124) + rxt(k,308)*y(k,205)
         mat(k,1673) = mat(k,1673) + .350_r8*rxt(k,317)*y(k,24)
         mat(k,655) = rxt(k,272)*y(k,77) + rxt(k,274)*y(k,89) + rxt(k,273)*y(k,234)
         mat(k,1087) = rxt(k,272)*y(k,73)
         mat(k,1734) = rxt(k,274)*y(k,73)
         mat(k,2145) = rxt(k,273)*y(k,73)
         mat(k,1745) = -(rxt(k,211)*y(k,223) + rxt(k,274)*y(k,73))
         mat(k,1723) = -rxt(k,211)*y(k,89)
         mat(k,661) = -rxt(k,274)*y(k,89)
         mat(k,2109) = rxt(k,295)*y(k,126)
         mat(k,983) = rxt(k,328)*y(k,126)
         mat(k,1158) = rxt(k,354)*y(k,126)
         mat(k,913) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,85)
         mat(k,197) = rxt(k,502)*y(k,126)
         mat(k,1341) = (rxt(k,554)+rxt(k,559)+rxt(k,564))*y(k,60)
         mat(k,1990) = rxt(k,210)*y(k,223)
         mat(k,1486) = rxt(k,295)*y(k,42) + rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) &
                      + rxt(k,502)*y(k,67)
         mat(k,1723) = mat(k,1723) + rxt(k,210)*y(k,125)
         mat(k,363) = -(rxt(k,187)*y(k,223))
         mat(k,1640) = -rxt(k,187)*y(k,90)
         mat(k,1963) = rxt(k,208)*y(k,205)
         mat(k,2022) = rxt(k,208)*y(k,125)
         mat(k,702) = -(rxt(k,265)*y(k,134) + (rxt(k,556) + rxt(k,561)) * y(k,85))
         mat(k,1896) = -rxt(k,265)*y(k,91)
         mat(k,1333) = -(rxt(k,556) + rxt(k,561)) * y(k,91)
         mat(k,2123) = rxt(k,257)*y(k,205)
         mat(k,2044) = rxt(k,257)*y(k,19)
         mat(k,756) = -(rxt(k,236)*y(k,56) + rxt(k,237)*y(k,134) + rxt(k,238)*y(k,223) &
                      + (rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,85))
         mat(k,1931) = -rxt(k,236)*y(k,92)
         mat(k,1898) = -rxt(k,237)*y(k,92)
         mat(k,1677) = -rxt(k,238)*y(k,92)
         mat(k,1334) = -(rxt(k,549) + rxt(k,555) + rxt(k,560)) * y(k,92)
         mat(k,1501) = rxt(k,225)*y(k,205)
         mat(k,907) = rxt(k,230)*y(k,223)
         mat(k,2048) = rxt(k,225)*y(k,59)
         mat(k,1677) = mat(k,1677) + rxt(k,230)*y(k,60)
         mat(k,990) = -(rxt(k,371)*y(k,223))
         mat(k,1696) = -rxt(k,371)*y(k,93)
         mat(k,478) = .300_r8*rxt(k,416)*y(k,223)
         mat(k,427) = .500_r8*rxt(k,417)*y(k,223)
         mat(k,1848) = rxt(k,370)*y(k,202) + rxt(k,377)*y(k,209)
         mat(k,471) = rxt(k,370)*y(k,124)
         mat(k,1225) = rxt(k,377)*y(k,124)
         mat(k,1696) = mat(k,1696) + .300_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
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
         mat(k,157) = -(rxt(k,402)*y(k,223))
         mat(k,1609) = -rxt(k,402)*y(k,94)
         mat(k,1007) = -(rxt(k,356)*y(k,223))
         mat(k,1698) = -rxt(k,356)*y(k,95)
         mat(k,479) = .700_r8*rxt(k,416)*y(k,223)
         mat(k,428) = .500_r8*rxt(k,417)*y(k,223)
         mat(k,461) = .500_r8*rxt(k,391)*y(k,223)
         mat(k,1850) = .050_r8*rxt(k,414)*y(k,208) + .220_r8*rxt(k,376)*y(k,209) &
                      + .250_r8*rxt(k,433)*y(k,231)
         mat(k,1463) = .050_r8*rxt(k,415)*y(k,208) + .220_r8*rxt(k,375)*y(k,209) &
                      + .250_r8*rxt(k,434)*y(k,231)
         mat(k,435) = .500_r8*rxt(k,360)*y(k,223)
         mat(k,1291) = .220_r8*rxt(k,372)*y(k,209) + .250_r8*rxt(k,430)*y(k,231)
         mat(k,1380) = .230_r8*rxt(k,373)*y(k,209) + .200_r8*rxt(k,361)*y(k,227) &
                      + .100_r8*rxt(k,431)*y(k,231)
         mat(k,1174) = .050_r8*rxt(k,414)*y(k,124) + .050_r8*rxt(k,415)*y(k,126)
         mat(k,1227) = .220_r8*rxt(k,376)*y(k,124) + .220_r8*rxt(k,375)*y(k,126) &
                      + .220_r8*rxt(k,372)*y(k,198) + .230_r8*rxt(k,373)*y(k,199)
         mat(k,1698) = mat(k,1698) + .700_r8*rxt(k,416)*y(k,99) + .500_r8*rxt(k,417) &
                      *y(k,100) + .500_r8*rxt(k,391)*y(k,109) + .500_r8*rxt(k,360) &
                      *y(k,148)
         mat(k,1076) = .200_r8*rxt(k,361)*y(k,199)
         mat(k,1105) = .250_r8*rxt(k,433)*y(k,124) + .250_r8*rxt(k,434)*y(k,126) &
                      + .250_r8*rxt(k,430)*y(k,198) + .100_r8*rxt(k,431)*y(k,199)
         mat(k,245) = -(rxt(k,403)*y(k,223))
         mat(k,1623) = -rxt(k,403)*y(k,96)
         mat(k,1806) = .870_r8*rxt(k,414)*y(k,208)
         mat(k,1444) = .950_r8*rxt(k,415)*y(k,208)
         mat(k,1283) = rxt(k,410)*y(k,208)
         mat(k,1363) = .750_r8*rxt(k,411)*y(k,208)
         mat(k,1163) = .870_r8*rxt(k,414)*y(k,124) + .950_r8*rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,198) + .750_r8*rxt(k,411)*y(k,199)
         mat(k,70) = -(rxt(k,404)*y(k,223))
         mat(k,1597) = -rxt(k,404)*y(k,97)
         mat(k,579) = .600_r8*rxt(k,427)*y(k,223)
         mat(k,1597) = mat(k,1597) + .600_r8*rxt(k,427)*y(k,103)
         mat(k,777) = -(rxt(k,418)*y(k,126) + rxt(k,425)*y(k,136) + rxt(k,426) &
                      *y(k,223))
         mat(k,1447) = -rxt(k,418)*y(k,98)
         mat(k,1540) = -rxt(k,425)*y(k,98)
         mat(k,1679) = -rxt(k,426)*y(k,98)
         mat(k,476) = -(rxt(k,416)*y(k,223))
         mat(k,1653) = -rxt(k,416)*y(k,99)
         mat(k,1819) = .080_r8*rxt(k,408)*y(k,207)
         mat(k,1195) = .080_r8*rxt(k,408)*y(k,124)
         mat(k,424) = -(rxt(k,417)*y(k,223))
         mat(k,1647) = -rxt(k,417)*y(k,100)
         mat(k,1817) = .080_r8*rxt(k,414)*y(k,208)
         mat(k,1164) = .080_r8*rxt(k,414)*y(k,124)
         mat(k,1126) = -(rxt(k,419)*y(k,198) + rxt(k,420)*y(k,199) + rxt(k,421) &
                      *y(k,205) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126))
         mat(k,1293) = -rxt(k,419)*y(k,101)
         mat(k,1386) = -rxt(k,420)*y(k,101)
         mat(k,2068) = -rxt(k,421)*y(k,101)
         mat(k,1856) = -rxt(k,422)*y(k,101)
         mat(k,1469) = -rxt(k,423)*y(k,101)
         mat(k,780) = rxt(k,418)*y(k,126)
         mat(k,1469) = mat(k,1469) + rxt(k,418)*y(k,98)
         mat(k,299) = -(rxt(k,424)*y(k,223))
         mat(k,1631) = -rxt(k,424)*y(k,102)
         mat(k,1118) = rxt(k,421)*y(k,205)
         mat(k,2012) = rxt(k,421)*y(k,101)
         mat(k,580) = -(rxt(k,427)*y(k,223))
         mat(k,1663) = -rxt(k,427)*y(k,103)
         mat(k,2037) = rxt(k,407)*y(k,207) + rxt(k,412)*y(k,208)
         mat(k,1196) = rxt(k,407)*y(k,205)
         mat(k,1166) = rxt(k,412)*y(k,205)
         mat(k,41) = -(rxt(k,541)*y(k,223))
         mat(k,1591) = -rxt(k,541)*y(k,104)
         mat(k,1142) = -(rxt(k,378)*y(k,136) + rxt(k,379)*y(k,223))
         mat(k,1559) = -rxt(k,378)*y(k,105)
         mat(k,1706) = -rxt(k,379)*y(k,105)
         mat(k,781) = .300_r8*rxt(k,425)*y(k,136)
         mat(k,1857) = .360_r8*rxt(k,408)*y(k,207)
         mat(k,1470) = .400_r8*rxt(k,409)*y(k,207)
         mat(k,1559) = mat(k,1559) + .300_r8*rxt(k,425)*y(k,98)
         mat(k,1294) = .390_r8*rxt(k,405)*y(k,207)
         mat(k,1387) = .310_r8*rxt(k,406)*y(k,207)
         mat(k,1204) = .360_r8*rxt(k,408)*y(k,124) + .400_r8*rxt(k,409)*y(k,126) &
                      + .390_r8*rxt(k,405)*y(k,198) + .310_r8*rxt(k,406)*y(k,199)
         mat(k,214) = -(rxt(k,380)*y(k,223))
         mat(k,1619) = -rxt(k,380)*y(k,106)
         mat(k,2004) = rxt(k,374)*y(k,209)
         mat(k,1222) = rxt(k,374)*y(k,205)
         mat(k,401) = -(rxt(k,389)*y(k,223))
         mat(k,1644) = -rxt(k,389)*y(k,107)
         mat(k,1815) = .800_r8*rxt(k,398)*y(k,192)
         mat(k,877) = .800_r8*rxt(k,398)*y(k,124)
         mat(k,219) = -(rxt(k,390)*y(k,223))
         mat(k,1620) = -rxt(k,390)*y(k,108)
         mat(k,2005) = .800_r8*rxt(k,387)*y(k,213)
         mat(k,571) = .800_r8*rxt(k,387)*y(k,205)
         mat(k,460) = -(rxt(k,391)*y(k,223))
         mat(k,1651) = -rxt(k,391)*y(k,109)
         mat(k,1966) = rxt(k,394)*y(k,211)
         mat(k,1267) = rxt(k,394)*y(k,125)
         mat(k,860) = -(rxt(k,482)*y(k,126) + rxt(k,483)*y(k,136) + rxt(k,484) &
                      *y(k,223))
         mat(k,1451) = -rxt(k,482)*y(k,110)
         mat(k,1543) = -rxt(k,483)*y(k,110)
         mat(k,1685) = -rxt(k,484)*y(k,110)
         mat(k,1250) = -(rxt(k,392)*y(k,136) + rxt(k,393)*y(k,223))
         mat(k,1564) = -rxt(k,392)*y(k,111)
         mat(k,1711) = -rxt(k,393)*y(k,111)
         mat(k,784) = .200_r8*rxt(k,425)*y(k,136)
         mat(k,1862) = .560_r8*rxt(k,408)*y(k,207)
         mat(k,1475) = .600_r8*rxt(k,409)*y(k,207)
         mat(k,1564) = mat(k,1564) + .200_r8*rxt(k,425)*y(k,98)
         mat(k,1299) = .610_r8*rxt(k,405)*y(k,207)
         mat(k,1392) = .440_r8*rxt(k,406)*y(k,207)
         mat(k,1208) = .560_r8*rxt(k,408)*y(k,124) + .600_r8*rxt(k,409)*y(k,126) &
                      + .610_r8*rxt(k,405)*y(k,198) + .440_r8*rxt(k,406)*y(k,199)
         mat(k,734) = -(rxt(k,190)*y(k,124) + (rxt(k,191) + rxt(k,192) + rxt(k,193) &
                      ) * y(k,125) + rxt(k,194)*y(k,135) + rxt(k,202)*y(k,223) &
                      + rxt(k,574)*y(k,222))
         mat(k,1835) = -rxt(k,190)*y(k,112)
         mat(k,1971) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112)
         mat(k,1420) = -rxt(k,194)*y(k,112)
         mat(k,1675) = -rxt(k,202)*y(k,112)
         mat(k,677) = -rxt(k,574)*y(k,112)
         mat(k,1897) = rxt(k,188)*y(k,214) + rxt(k,571)*y(k,217)
         mat(k,1420) = mat(k,1420) + rxt(k,572)*y(k,217)
         mat(k,695) = 1.100_r8*rxt(k,567)*y(k,215) + .200_r8*rxt(k,565)*y(k,216)
         mat(k,420) = rxt(k,188)*y(k,134)
         mat(k,554) = 1.100_r8*rxt(k,567)*y(k,201)
         mat(k,685) = .200_r8*rxt(k,565)*y(k,201)
         mat(k,390) = rxt(k,571)*y(k,134) + rxt(k,572)*y(k,135)
         mat(k,1960) = rxt(k,209)*y(k,126)
         mat(k,1442) = rxt(k,209)*y(k,125)
         mat(k,269) = -(rxt(k,428)*y(k,223))
         mat(k,1627) = -rxt(k,428)*y(k,115)
         mat(k,1117) = .200_r8*rxt(k,420)*y(k,199)
         mat(k,1364) = .200_r8*rxt(k,420)*y(k,101)
         mat(k,950) = -(rxt(k,429)*y(k,223))
         mat(k,1692) = -rxt(k,429)*y(k,116)
         mat(k,1122) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) + rxt(k,419)*y(k,198) &
                      + .800_r8*rxt(k,420)*y(k,199)
         mat(k,1844) = rxt(k,422)*y(k,101)
         mat(k,1457) = rxt(k,423)*y(k,101)
         mat(k,1288) = rxt(k,419)*y(k,101)
         mat(k,1375) = .800_r8*rxt(k,420)*y(k,101)
         mat(k,51) = -(rxt(k,519)*y(k,223))
         mat(k,1593) = -rxt(k,519)*y(k,120)
         mat(k,1875) = -(rxt(k,190)*y(k,112) + rxt(k,199)*y(k,126) + rxt(k,203) &
                      *y(k,205) + rxt(k,204)*y(k,136) + rxt(k,205)*y(k,134) + rxt(k,226) &
                      *y(k,59) + rxt(k,258)*y(k,19) + rxt(k,301)*y(k,199) + rxt(k,310) &
                      *y(k,206) + rxt(k,323)*y(k,195) + rxt(k,334)*y(k,198) + rxt(k,338) &
                      *y(k,204) + rxt(k,351)*y(k,196) + rxt(k,359)*y(k,226) + rxt(k,363) &
                      *y(k,227) + (rxt(k,369) + rxt(k,370)) * y(k,202) + (rxt(k,376) &
                      + rxt(k,377)) * y(k,209) + rxt(k,385)*y(k,211) + rxt(k,388) &
                      *y(k,213) + (rxt(k,398) + rxt(k,399)) * y(k,192) + rxt(k,408) &
                      *y(k,207) + rxt(k,414)*y(k,208) + rxt(k,422)*y(k,101) + rxt(k,433) &
                      *y(k,231) + rxt(k,437)*y(k,191) + rxt(k,440)*y(k,193) + rxt(k,445) &
                      *y(k,194) + rxt(k,447)*y(k,197) + rxt(k,451)*y(k,200) + rxt(k,454) &
                      *y(k,210) + rxt(k,457)*y(k,212) + rxt(k,460)*y(k,225) + rxt(k,467) &
                      *y(k,230) + rxt(k,473)*y(k,232) + rxt(k,476)*y(k,233) + rxt(k,487) &
                      *y(k,218) + rxt(k,492)*y(k,228) + rxt(k,497)*y(k,229) + rxt(k,576) &
                      *y(k,222))
         mat(k,740) = -rxt(k,190)*y(k,124)
         mat(k,1489) = -rxt(k,199)*y(k,124)
         mat(k,2088) = -rxt(k,203)*y(k,124)
         mat(k,1578) = -rxt(k,204)*y(k,124)
         mat(k,1916) = -rxt(k,205)*y(k,124)
         mat(k,1516) = -rxt(k,226)*y(k,124)
         mat(k,2136) = -rxt(k,258)*y(k,124)
         mat(k,1403) = -rxt(k,301)*y(k,124)
         mat(k,333) = -rxt(k,310)*y(k,124)
         mat(k,807) = -rxt(k,323)*y(k,124)
         mat(k,1309) = -rxt(k,334)*y(k,124)
         mat(k,670) = -rxt(k,338)*y(k,124)
         mat(k,750) = -rxt(k,351)*y(k,124)
         mat(k,724) = -rxt(k,359)*y(k,124)
         mat(k,1082) = -rxt(k,363)*y(k,124)
         mat(k,472) = -(rxt(k,369) + rxt(k,370)) * y(k,124)
         mat(k,1238) = -(rxt(k,376) + rxt(k,377)) * y(k,124)
         mat(k,1278) = -rxt(k,385)*y(k,124)
         mat(k,576) = -rxt(k,388)*y(k,124)
         mat(k,889) = -(rxt(k,398) + rxt(k,399)) * y(k,124)
         mat(k,1217) = -rxt(k,408)*y(k,124)
         mat(k,1190) = -rxt(k,414)*y(k,124)
         mat(k,1136) = -rxt(k,422)*y(k,124)
         mat(k,1113) = -rxt(k,433)*y(k,124)
         mat(k,416) = -rxt(k,437)*y(k,124)
         mat(k,383) = -rxt(k,440)*y(k,124)
         mat(k,328) = -rxt(k,445)*y(k,124)
         mat(k,534) = -rxt(k,447)*y(k,124)
         mat(k,652) = -rxt(k,451)*y(k,124)
         mat(k,612) = -rxt(k,454)*y(k,124)
         mat(k,822) = -rxt(k,457)*y(k,124)
         mat(k,341) = -rxt(k,460)*y(k,124)
         mat(k,627) = -rxt(k,467)*y(k,124)
         mat(k,644) = -rxt(k,473)*y(k,124)
         mat(k,398) = -rxt(k,476)*y(k,124)
         mat(k,1069) = -rxt(k,487)*y(k,124)
         mat(k,1049) = -rxt(k,492)*y(k,124)
         mat(k,1030) = -rxt(k,497)*y(k,124)
         mat(k,680) = -rxt(k,576)*y(k,124)
         mat(k,740) = mat(k,740) + 2.000_r8*rxt(k,192)*y(k,125) + rxt(k,194)*y(k,135) &
                      + rxt(k,202)*y(k,223)
         mat(k,1993) = 2.000_r8*rxt(k,192)*y(k,112) + rxt(k,195)*y(k,134) + rxt(k,511) &
                      *y(k,152)
         mat(k,1916) = mat(k,1916) + rxt(k,195)*y(k,125)
         mat(k,1432) = rxt(k,194)*y(k,112) + rxt(k,189)*y(k,214)
         mat(k,1324) = rxt(k,511)*y(k,125)
         mat(k,423) = rxt(k,189)*y(k,135)
         mat(k,1726) = rxt(k,202)*y(k,112)
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
         mat(k,1996) = -((rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,112) + (rxt(k,195) &
                      + rxt(k,197)) * y(k,134) + rxt(k,196)*y(k,136) + rxt(k,208) &
                      *y(k,205) + rxt(k,209)*y(k,126) + rxt(k,210)*y(k,223) + rxt(k,228) &
                      *y(k,59) + rxt(k,259)*y(k,19) + rxt(k,345)*y(k,198) + rxt(k,394) &
                      *y(k,211) + rxt(k,452)*y(k,200) + rxt(k,455)*y(k,210) + rxt(k,458) &
                      *y(k,212) + rxt(k,462)*y(k,143) + rxt(k,465)*y(k,191) + rxt(k,511) &
                      *y(k,152))
         mat(k,742) = -(rxt(k,191) + rxt(k,192) + rxt(k,193)) * y(k,125)
         mat(k,1919) = -(rxt(k,195) + rxt(k,197)) * y(k,125)
         mat(k,1581) = -rxt(k,196)*y(k,125)
         mat(k,2091) = -rxt(k,208)*y(k,125)
         mat(k,1492) = -rxt(k,209)*y(k,125)
         mat(k,1729) = -rxt(k,210)*y(k,125)
         mat(k,1519) = -rxt(k,228)*y(k,125)
         mat(k,2139) = -rxt(k,259)*y(k,125)
         mat(k,1310) = -rxt(k,345)*y(k,125)
         mat(k,1279) = -rxt(k,394)*y(k,125)
         mat(k,653) = -rxt(k,452)*y(k,125)
         mat(k,613) = -rxt(k,455)*y(k,125)
         mat(k,823) = -rxt(k,458)*y(k,125)
         mat(k,362) = -rxt(k,462)*y(k,125)
         mat(k,417) = -rxt(k,465)*y(k,125)
         mat(k,1327) = -rxt(k,511)*y(k,125)
         mat(k,545) = rxt(k,396)*y(k,223)
         mat(k,267) = rxt(k,367)*y(k,126)
         mat(k,2139) = mat(k,2139) + rxt(k,258)*y(k,124)
         mat(k,1519) = mat(k,1519) + rxt(k,226)*y(k,124)
         mat(k,367) = rxt(k,187)*y(k,223)
         mat(k,482) = .700_r8*rxt(k,416)*y(k,223)
         mat(k,1137) = rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126)
         mat(k,1878) = rxt(k,258)*y(k,19) + rxt(k,226)*y(k,59) + rxt(k,422)*y(k,101) &
                      + 2.000_r8*rxt(k,199)*y(k,126) + rxt(k,205)*y(k,134) &
                      + rxt(k,204)*y(k,136) + rxt(k,437)*y(k,191) + rxt(k,398) &
                      *y(k,192) + rxt(k,440)*y(k,193) + rxt(k,445)*y(k,194) &
                      + rxt(k,323)*y(k,195) + rxt(k,351)*y(k,196) + rxt(k,447) &
                      *y(k,197) + rxt(k,334)*y(k,198) + rxt(k,301)*y(k,199) &
                      + rxt(k,451)*y(k,200) + rxt(k,369)*y(k,202) + rxt(k,338) &
                      *y(k,204) + rxt(k,203)*y(k,205) + rxt(k,310)*y(k,206) &
                      + .920_r8*rxt(k,408)*y(k,207) + .920_r8*rxt(k,414)*y(k,208) &
                      + rxt(k,376)*y(k,209) + rxt(k,454)*y(k,210) + rxt(k,385) &
                      *y(k,211) + rxt(k,457)*y(k,212) + rxt(k,388)*y(k,213) &
                      + 1.600_r8*rxt(k,487)*y(k,218) + rxt(k,460)*y(k,225) &
                      + rxt(k,359)*y(k,226) + rxt(k,363)*y(k,227) + .900_r8*rxt(k,492) &
                      *y(k,228) + .800_r8*rxt(k,497)*y(k,229) + rxt(k,467)*y(k,230) &
                      + rxt(k,433)*y(k,231) + rxt(k,473)*y(k,232) + rxt(k,476) &
                      *y(k,233)
         mat(k,1492) = mat(k,1492) + rxt(k,367)*y(k,16) + rxt(k,423)*y(k,101) &
                      + 2.000_r8*rxt(k,199)*y(k,124) + rxt(k,200)*y(k,134) &
                      + rxt(k,198)*y(k,205) + rxt(k,409)*y(k,207) + rxt(k,415) &
                      *y(k,208) + rxt(k,375)*y(k,209) + rxt(k,386)*y(k,211) &
                      + 2.000_r8*rxt(k,488)*y(k,218) + rxt(k,201)*y(k,223) &
                      + rxt(k,434)*y(k,231)
         mat(k,797) = rxt(k,357)*y(k,223)
         mat(k,1919) = mat(k,1919) + rxt(k,205)*y(k,124) + rxt(k,200)*y(k,126)
         mat(k,1581) = mat(k,1581) + rxt(k,204)*y(k,124)
         mat(k,527) = rxt(k,494)*y(k,223)
         mat(k,417) = mat(k,417) + rxt(k,437)*y(k,124)
         mat(k,890) = rxt(k,398)*y(k,124)
         mat(k,384) = rxt(k,440)*y(k,124)
         mat(k,329) = rxt(k,445)*y(k,124)
         mat(k,808) = rxt(k,323)*y(k,124)
         mat(k,751) = rxt(k,351)*y(k,124)
         mat(k,535) = rxt(k,447)*y(k,124)
         mat(k,1310) = mat(k,1310) + rxt(k,334)*y(k,124)
         mat(k,1405) = rxt(k,301)*y(k,124) + .500_r8*rxt(k,485)*y(k,218)
         mat(k,653) = mat(k,653) + rxt(k,451)*y(k,124)
         mat(k,473) = rxt(k,369)*y(k,124)
         mat(k,671) = rxt(k,338)*y(k,124)
         mat(k,2091) = mat(k,2091) + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126)
         mat(k,334) = rxt(k,310)*y(k,124)
         mat(k,1218) = .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126)
         mat(k,1191) = .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126)
         mat(k,1239) = rxt(k,376)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,613) = mat(k,613) + rxt(k,454)*y(k,124)
         mat(k,1279) = mat(k,1279) + rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126)
         mat(k,823) = mat(k,823) + rxt(k,457)*y(k,124)
         mat(k,577) = rxt(k,388)*y(k,124)
         mat(k,1070) = 1.600_r8*rxt(k,487)*y(k,124) + 2.000_r8*rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1729) = mat(k,1729) + rxt(k,396)*y(k,1) + rxt(k,187)*y(k,90) &
                      + .700_r8*rxt(k,416)*y(k,99) + rxt(k,201)*y(k,126) + rxt(k,357) &
                      *y(k,127) + rxt(k,494)*y(k,178)
         mat(k,342) = rxt(k,460)*y(k,124)
         mat(k,725) = rxt(k,359)*y(k,124)
         mat(k,1083) = rxt(k,363)*y(k,124)
         mat(k,1050) = .900_r8*rxt(k,492)*y(k,124)
         mat(k,1031) = .800_r8*rxt(k,497)*y(k,124)
         mat(k,628) = rxt(k,467)*y(k,124)
         mat(k,1114) = rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126)
         mat(k,645) = rxt(k,473)*y(k,124)
         mat(k,399) = rxt(k,476)*y(k,124)
         mat(k,1482) = -(rxt(k,198)*y(k,205) + rxt(k,199)*y(k,124) + rxt(k,200) &
                      *y(k,134) + rxt(k,201)*y(k,223) + rxt(k,209)*y(k,125) + rxt(k,295) &
                      *y(k,42) + rxt(k,328)*y(k,45) + rxt(k,347)*y(k,29) + rxt(k,354) &
                      *y(k,49) + rxt(k,367)*y(k,16) + rxt(k,375)*y(k,209) + rxt(k,386) &
                      *y(k,211) + rxt(k,409)*y(k,207) + rxt(k,415)*y(k,208) + rxt(k,418) &
                      *y(k,98) + rxt(k,423)*y(k,101) + rxt(k,434)*y(k,231) + rxt(k,479) &
                      *y(k,6) + rxt(k,482)*y(k,110) + rxt(k,488)*y(k,218) + rxt(k,499) &
                      *y(k,180) + rxt(k,502)*y(k,67))
         mat(k,2081) = -rxt(k,198)*y(k,126)
         mat(k,1868) = -rxt(k,199)*y(k,126)
         mat(k,1909) = -rxt(k,200)*y(k,126)
         mat(k,1719) = -rxt(k,201)*y(k,126)
         mat(k,1986) = -rxt(k,209)*y(k,126)
         mat(k,2105) = -rxt(k,295)*y(k,126)
         mat(k,981) = -rxt(k,328)*y(k,126)
         mat(k,938) = -rxt(k,347)*y(k,126)
         mat(k,1156) = -rxt(k,354)*y(k,126)
         mat(k,265) = -rxt(k,367)*y(k,126)
         mat(k,1234) = -rxt(k,375)*y(k,126)
         mat(k,1273) = -rxt(k,386)*y(k,126)
         mat(k,1212) = -rxt(k,409)*y(k,126)
         mat(k,1185) = -rxt(k,415)*y(k,126)
         mat(k,787) = -rxt(k,418)*y(k,126)
         mat(k,1132) = -rxt(k,423)*y(k,126)
         mat(k,1110) = -rxt(k,434)*y(k,126)
         mat(k,845) = -rxt(k,479)*y(k,126)
         mat(k,871) = -rxt(k,482)*y(k,126)
         mat(k,1065) = -rxt(k,488)*y(k,126)
         mat(k,923) = -rxt(k,499)*y(k,126)
         mat(k,195) = -rxt(k,502)*y(k,126)
         mat(k,443) = rxt(k,260)*y(k,134)
         mat(k,1944) = rxt(k,227)*y(k,60)
         mat(k,910) = rxt(k,227)*y(k,56) + rxt(k,229)*y(k,134) + rxt(k,230)*y(k,223)
         mat(k,659) = rxt(k,274)*y(k,89)
         mat(k,1741) = rxt(k,274)*y(k,73) + rxt(k,211)*y(k,223)
         mat(k,463) = .500_r8*rxt(k,391)*y(k,223)
         mat(k,1986) = mat(k,1986) + rxt(k,197)*y(k,134) + rxt(k,196)*y(k,136)
         mat(k,1909) = mat(k,1909) + rxt(k,260)*y(k,20) + rxt(k,229)*y(k,60) &
                      + rxt(k,197)*y(k,125)
         mat(k,1571) = rxt(k,196)*y(k,125)
         mat(k,355) = rxt(k,343)*y(k,223)
         mat(k,1719) = mat(k,1719) + rxt(k,230)*y(k,60) + rxt(k,211)*y(k,89) &
                      + .500_r8*rxt(k,391)*y(k,109) + rxt(k,343)*y(k,141)
         mat(k,793) = -(rxt(k,357)*y(k,223))
         mat(k,1680) = -rxt(k,357)*y(k,127)
         mat(k,929) = rxt(k,347)*y(k,126)
         mat(k,425) = .500_r8*rxt(k,417)*y(k,223)
         mat(k,301) = rxt(k,424)*y(k,223)
         mat(k,270) = rxt(k,428)*y(k,223)
         mat(k,947) = rxt(k,429)*y(k,223)
         mat(k,1448) = rxt(k,347)*y(k,29)
         mat(k,1680) = mat(k,1680) + .500_r8*rxt(k,417)*y(k,100) + rxt(k,424)*y(k,102) &
                      + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116)
         mat(k,293) = -(rxt(k,489)*y(k,223))
         mat(k,1630) = -rxt(k,489)*y(k,128)
         mat(k,2011) = rxt(k,486)*y(k,218)
         mat(k,1054) = rxt(k,486)*y(k,205)
         mat(k,1917) = -(rxt(k,167)*y(k,136) + 4._r8*rxt(k,168)*y(k,134) + rxt(k,169) &
                      *y(k,135) + rxt(k,170)*y(k,77) + rxt(k,171)*y(k,79) + rxt(k,176) &
                      *y(k,205) + rxt(k,182)*y(k,223) + (rxt(k,195) + rxt(k,197) &
                      ) * y(k,125) + rxt(k,200)*y(k,126) + rxt(k,205)*y(k,124) &
                      + rxt(k,229)*y(k,60) + rxt(k,231)*y(k,59) + rxt(k,234)*y(k,85) &
                      + rxt(k,237)*y(k,92) + rxt(k,260)*y(k,20) + rxt(k,261)*y(k,19) &
                      + rxt(k,263)*y(k,81) + rxt(k,265)*y(k,91) + rxt(k,296)*y(k,42) &
                      + rxt(k,504)*y(k,139) + (rxt(k,569) + rxt(k,570)) * y(k,215) &
                      + rxt(k,571)*y(k,217))
         mat(k,1579) = -rxt(k,167)*y(k,134)
         mat(k,1433) = -rxt(k,169)*y(k,134)
         mat(k,1097) = -rxt(k,170)*y(k,134)
         mat(k,508) = -rxt(k,171)*y(k,134)
         mat(k,2089) = -rxt(k,176)*y(k,134)
         mat(k,1727) = -rxt(k,182)*y(k,134)
         mat(k,1994) = -(rxt(k,195) + rxt(k,197)) * y(k,134)
         mat(k,1490) = -rxt(k,200)*y(k,134)
         mat(k,1876) = -rxt(k,205)*y(k,134)
         mat(k,914) = -rxt(k,229)*y(k,134)
         mat(k,1517) = -rxt(k,231)*y(k,134)
         mat(k,1344) = -rxt(k,234)*y(k,134)
         mat(k,760) = -rxt(k,237)*y(k,134)
         mat(k,445) = -rxt(k,260)*y(k,134)
         mat(k,2137) = -rxt(k,261)*y(k,134)
         mat(k,768) = -rxt(k,263)*y(k,134)
         mat(k,706) = -rxt(k,265)*y(k,134)
         mat(k,2113) = -rxt(k,296)*y(k,134)
         mat(k,260) = -rxt(k,504)*y(k,134)
         mat(k,559) = -(rxt(k,569) + rxt(k,570)) * y(k,134)
         mat(k,392) = -rxt(k,571)*y(k,134)
         mat(k,1769) = rxt(k,174)*y(k,205)
         mat(k,741) = rxt(k,190)*y(k,124) + rxt(k,191)*y(k,125) + rxt(k,194)*y(k,135) &
                      + rxt(k,574)*y(k,222)
         mat(k,1876) = mat(k,1876) + rxt(k,190)*y(k,112)
         mat(k,1994) = mat(k,1994) + rxt(k,191)*y(k,112)
         mat(k,1433) = mat(k,1433) + rxt(k,194)*y(k,112) + rxt(k,506)*y(k,150) &
                      + rxt(k,512)*y(k,152) + rxt(k,573)*y(k,217) + (rxt(k,156) &
                       +rxt(k,157))*y(k,219) + rxt(k,579)*y(k,224)
         mat(k,593) = rxt(k,506)*y(k,135)
         mat(k,1325) = rxt(k,512)*y(k,135)
         mat(k,700) = rxt(k,565)*y(k,216) + 1.150_r8*rxt(k,566)*y(k,222)
         mat(k,2089) = mat(k,2089) + rxt(k,174)*y(k,76)
         mat(k,689) = rxt(k,565)*y(k,201)
         mat(k,392) = mat(k,392) + rxt(k,573)*y(k,135)
         mat(k,1795) = (rxt(k,156)+rxt(k,157))*y(k,135)
         mat(k,681) = rxt(k,574)*y(k,112) + 1.150_r8*rxt(k,566)*y(k,201)
         mat(k,1727) = mat(k,1727) + 2.000_r8*rxt(k,184)*y(k,223)
         mat(k,518) = rxt(k,579)*y(k,135)
         mat(k,1426) = -(rxt(k,156)*y(k,219) + rxt(k,161)*y(k,220) + rxt(k,169) &
                      *y(k,134) + rxt(k,175)*y(k,76) + rxt(k,189)*y(k,214) + rxt(k,194) &
                      *y(k,112) + rxt(k,340)*y(k,203) + rxt(k,506)*y(k,150) + rxt(k,512) &
                      *y(k,152) + rxt(k,568)*y(k,215) + (rxt(k,572) + rxt(k,573) &
                      ) * y(k,217) + rxt(k,579)*y(k,224))
         mat(k,1786) = -rxt(k,156)*y(k,135)
         mat(k,77) = -rxt(k,161)*y(k,135)
         mat(k,1908) = -rxt(k,169)*y(k,135)
         mat(k,1760) = -rxt(k,175)*y(k,135)
         mat(k,421) = -rxt(k,189)*y(k,135)
         mat(k,736) = -rxt(k,194)*y(k,135)
         mat(k,349) = -rxt(k,340)*y(k,135)
         mat(k,589) = -rxt(k,506)*y(k,135)
         mat(k,1319) = -rxt(k,512)*y(k,135)
         mat(k,556) = -rxt(k,568)*y(k,135)
         mat(k,391) = -(rxt(k,572) + rxt(k,573)) * y(k,135)
         mat(k,517) = -rxt(k,579)*y(k,135)
         mat(k,1352) = rxt(k,252)*y(k,136) + rxt(k,251)*y(k,205)
         mat(k,2128) = 2.000_r8*rxt(k,253)*y(k,19) + (rxt(k,255)+rxt(k,256))*y(k,59) &
                      + rxt(k,261)*y(k,134) + rxt(k,257)*y(k,205)
         mat(k,1943) = rxt(k,220)*y(k,136) + rxt(k,218)*y(k,205)
         mat(k,1508) = (rxt(k,255)+rxt(k,256))*y(k,19) + (2.000_r8*rxt(k,222) &
                       +2.000_r8*rxt(k,223))*y(k,59) + rxt(k,231)*y(k,134) &
                      + rxt(k,225)*y(k,205) + rxt(k,233)*y(k,223)
         mat(k,1760) = mat(k,1760) + rxt(k,178)*y(k,136) + rxt(k,172)*y(k,205)
         mat(k,364) = rxt(k,187)*y(k,223)
         mat(k,736) = mat(k,736) + rxt(k,193)*y(k,125)
         mat(k,1867) = rxt(k,204)*y(k,136) + rxt(k,576)*y(k,222)
         mat(k,1985) = rxt(k,193)*y(k,112) + rxt(k,195)*y(k,134) + rxt(k,196)*y(k,136)
         mat(k,1481) = rxt(k,200)*y(k,134) + rxt(k,198)*y(k,205)
         mat(k,1908) = mat(k,1908) + rxt(k,261)*y(k,19) + rxt(k,231)*y(k,59) &
                      + rxt(k,195)*y(k,125) + rxt(k,200)*y(k,126) &
                      + 2.000_r8*rxt(k,168)*y(k,134) + 2.000_r8*rxt(k,167)*y(k,136) &
                      + rxt(k,176)*y(k,205) + rxt(k,160)*y(k,220) + rxt(k,182) &
                      *y(k,223)
         mat(k,1426) = mat(k,1426) + 2.000_r8*rxt(k,161)*y(k,220)
         mat(k,1570) = rxt(k,252)*y(k,17) + rxt(k,220)*y(k,56) + rxt(k,178)*y(k,76) &
                      + rxt(k,204)*y(k,124) + rxt(k,196)*y(k,125) &
                      + 2.000_r8*rxt(k,167)*y(k,134) + rxt(k,507)*y(k,150) &
                      + rxt(k,513)*y(k,152) + 2.000_r8*rxt(k,177)*y(k,205) &
                      + 2.000_r8*rxt(k,158)*y(k,219) + rxt(k,183)*y(k,223)
         mat(k,589) = mat(k,589) + rxt(k,507)*y(k,136)
         mat(k,1319) = mat(k,1319) + rxt(k,513)*y(k,136)
         mat(k,805) = rxt(k,322)*y(k,205)
         mat(k,748) = rxt(k,350)*y(k,205)
         mat(k,1396) = rxt(k,300)*y(k,205)
         mat(k,2080) = rxt(k,251)*y(k,17) + rxt(k,257)*y(k,19) + rxt(k,218)*y(k,56) &
                      + rxt(k,225)*y(k,59) + rxt(k,172)*y(k,76) + rxt(k,198)*y(k,126) &
                      + rxt(k,176)*y(k,134) + 2.000_r8*rxt(k,177)*y(k,136) &
                      + rxt(k,322)*y(k,195) + rxt(k,350)*y(k,196) + rxt(k,300) &
                      *y(k,199) + 2.000_r8*rxt(k,186)*y(k,205) + rxt(k,181)*y(k,223) &
                      + rxt(k,358)*y(k,226)
         mat(k,1786) = mat(k,1786) + 2.000_r8*rxt(k,158)*y(k,136)
         mat(k,77) = mat(k,77) + rxt(k,160)*y(k,134) + 2.000_r8*rxt(k,161)*y(k,135)
         mat(k,678) = rxt(k,576)*y(k,124)
         mat(k,1718) = rxt(k,233)*y(k,59) + rxt(k,187)*y(k,90) + rxt(k,182)*y(k,134) &
                      + rxt(k,183)*y(k,136) + rxt(k,181)*y(k,205)
         mat(k,722) = rxt(k,358)*y(k,205)
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
         mat(k,1573) = -(rxt(k,158)*y(k,219) + rxt(k,167)*y(k,134) + rxt(k,177) &
                      *y(k,205) + rxt(k,178)*y(k,76) + rxt(k,183)*y(k,223) + rxt(k,196) &
                      *y(k,125) + rxt(k,204)*y(k,124) + rxt(k,220)*y(k,56) + rxt(k,252) &
                      *y(k,17) + rxt(k,319)*y(k,25) + rxt(k,348)*y(k,29) + rxt(k,378) &
                      *y(k,105) + rxt(k,392)*y(k,111) + rxt(k,425)*y(k,98) + rxt(k,463) &
                      *y(k,143) + rxt(k,480)*y(k,6) + rxt(k,483)*y(k,110) + rxt(k,507) &
                      *y(k,150) + rxt(k,513)*y(k,152))
         mat(k,1789) = -rxt(k,158)*y(k,136)
         mat(k,1911) = -rxt(k,167)*y(k,136)
         mat(k,2083) = -rxt(k,177)*y(k,136)
         mat(k,1763) = -rxt(k,178)*y(k,136)
         mat(k,1721) = -rxt(k,183)*y(k,136)
         mat(k,1988) = -rxt(k,196)*y(k,136)
         mat(k,1870) = -rxt(k,204)*y(k,136)
         mat(k,1946) = -rxt(k,220)*y(k,136)
         mat(k,1353) = -rxt(k,252)*y(k,136)
         mat(k,456) = -rxt(k,319)*y(k,136)
         mat(k,939) = -rxt(k,348)*y(k,136)
         mat(k,1147) = -rxt(k,378)*y(k,136)
         mat(k,1256) = -rxt(k,392)*y(k,136)
         mat(k,788) = -rxt(k,425)*y(k,136)
         mat(k,361) = -rxt(k,463)*y(k,136)
         mat(k,846) = -rxt(k,480)*y(k,136)
         mat(k,872) = -rxt(k,483)*y(k,136)
         mat(k,590) = -rxt(k,507)*y(k,136)
         mat(k,1321) = -rxt(k,513)*y(k,136)
         mat(k,1911) = mat(k,1911) + rxt(k,169)*y(k,135)
         mat(k,1428) = rxt(k,169)*y(k,134)
         mat(k,1305) = .150_r8*rxt(k,333)*y(k,205)
         mat(k,2083) = mat(k,2083) + .150_r8*rxt(k,333)*y(k,198) + .150_r8*rxt(k,383) &
                      *y(k,211)
         mat(k,1274) = .150_r8*rxt(k,383)*y(k,205)
         mat(k,229) = -(rxt(k,514)*y(k,152))
         mat(k,1314) = -rxt(k,514)*y(k,138)
         mat(k,2121) = rxt(k,254)*y(k,59)
         mat(k,1500) = rxt(k,254)*y(k,19) + 2.000_r8*rxt(k,224)*y(k,59)
         mat(k,253) = -(rxt(k,504)*y(k,134) + rxt(k,505)*y(k,223))
         mat(k,1885) = -rxt(k,504)*y(k,139)
         mat(k,1625) = -rxt(k,505)*y(k,139)
         mat(k,986) = rxt(k,371)*y(k,223)
         mat(k,1802) = .100_r8*rxt(k,492)*y(k,228)
         mat(k,1611) = rxt(k,371)*y(k,93)
         mat(k,1035) = .100_r8*rxt(k,492)*y(k,124)
         mat(k,352) = -(rxt(k,343)*y(k,223))
         mat(k,1638) = -rxt(k,343)*y(k,141)
         mat(k,1961) = rxt(k,345)*y(k,198)
         mat(k,1284) = rxt(k,345)*y(k,125)
         mat(k,1959) = rxt(k,465)*y(k,191)
         mat(k,412) = rxt(k,465)*y(k,125)
         mat(k,359) = -(rxt(k,462)*y(k,125) + rxt(k,463)*y(k,136))
         mat(k,1962) = -rxt(k,462)*y(k,143)
         mat(k,1533) = -rxt(k,463)*y(k,143)
         mat(k,121) = .070_r8*rxt(k,449)*y(k,223)
         mat(k,1812) = rxt(k,447)*y(k,197)
         mat(k,98) = .060_r8*rxt(k,461)*y(k,223)
         mat(k,150) = .070_r8*rxt(k,477)*y(k,223)
         mat(k,530) = rxt(k,447)*y(k,124)
         mat(k,1639) = .070_r8*rxt(k,449)*y(k,66) + .060_r8*rxt(k,461)*y(k,144) &
                      + .070_r8*rxt(k,477)*y(k,187)
         mat(k,96) = -(rxt(k,461)*y(k,223))
         mat(k,1600) = -rxt(k,461)*y(k,144)
         mat(k,88) = .530_r8*rxt(k,438)*y(k,223)
         mat(k,1600) = mat(k,1600) + .530_r8*rxt(k,438)*y(k,7)
         mat(k,234) = -(rxt(k,464)*y(k,223))
         mat(k,1621) = -rxt(k,464)*y(k,145)
         mat(k,2006) = rxt(k,459)*y(k,225)
         mat(k,337) = rxt(k,459)*y(k,205)
         mat(k,432) = -(rxt(k,360)*y(k,223))
         mat(k,1648) = -rxt(k,360)*y(k,148)
         mat(k,2028) = rxt(k,358)*y(k,226)
         mat(k,718) = rxt(k,358)*y(k,205)
         mat(k,305) = -(rxt(k,364)*y(k,223))
         mat(k,1632) = -rxt(k,364)*y(k,149)
         mat(k,2013) = .850_r8*rxt(k,362)*y(k,227)
         mat(k,1074) = .850_r8*rxt(k,362)*y(k,205)
         mat(k,587) = -(rxt(k,506)*y(k,135) + rxt(k,507)*y(k,136) + rxt(k,510) &
                      *y(k,223))
         mat(k,1416) = -rxt(k,506)*y(k,150)
         mat(k,1537) = -rxt(k,507)*y(k,150)
         mat(k,1664) = -rxt(k,510)*y(k,150)
         mat(k,1317) = -(rxt(k,508)*y(k,19) + rxt(k,509)*y(k,59) + rxt(k,511)*y(k,125) &
                      + rxt(k,512)*y(k,135) + rxt(k,513)*y(k,136) + rxt(k,514) &
                      *y(k,138) + rxt(k,515)*y(k,223))
         mat(k,2125) = -rxt(k,508)*y(k,152)
         mat(k,1504) = -rxt(k,509)*y(k,152)
         mat(k,1981) = -rxt(k,511)*y(k,152)
         mat(k,1424) = -rxt(k,512)*y(k,152)
         mat(k,1567) = -rxt(k,513)*y(k,152)
         mat(k,231) = -rxt(k,514)*y(k,152)
         mat(k,1714) = -rxt(k,515)*y(k,152)
         mat(k,1904) = rxt(k,504)*y(k,139)
         mat(k,1424) = mat(k,1424) + rxt(k,506)*y(k,150)
         mat(k,1567) = mat(k,1567) + rxt(k,507)*y(k,150)
         mat(k,257) = rxt(k,504)*y(k,134)
         mat(k,588) = rxt(k,506)*y(k,135) + rxt(k,507)*y(k,136) + rxt(k,510)*y(k,223)
         mat(k,1714) = mat(k,1714) + rxt(k,510)*y(k,150)
         mat(k,894) = -(rxt(k,517)*y(k,223))
         mat(k,1687) = -rxt(k,517)*y(k,153)
         mat(k,2124) = rxt(k,508)*y(k,152)
         mat(k,1502) = rxt(k,509)*y(k,152)
         mat(k,194) = rxt(k,502)*y(k,126) + (rxt(k,503)+.500_r8*rxt(k,516))*y(k,223)
         mat(k,1974) = rxt(k,511)*y(k,152)
         mat(k,1453) = rxt(k,502)*y(k,67)
         mat(k,1421) = rxt(k,512)*y(k,152)
         mat(k,1545) = rxt(k,513)*y(k,152)
         mat(k,230) = rxt(k,514)*y(k,152)
         mat(k,255) = rxt(k,505)*y(k,223)
         mat(k,1316) = rxt(k,508)*y(k,19) + rxt(k,509)*y(k,59) + rxt(k,511)*y(k,125) &
                      + rxt(k,512)*y(k,135) + rxt(k,513)*y(k,136) + rxt(k,514) &
                      *y(k,138) + rxt(k,515)*y(k,223)
         mat(k,1687) = mat(k,1687) + (rxt(k,503)+.500_r8*rxt(k,516))*y(k,67) &
                      + rxt(k,505)*y(k,139) + rxt(k,515)*y(k,152)
         mat(k,175) = -(rxt(k,518)*y(k,234))
         mat(k,2146) = -rxt(k,518)*y(k,154)
         mat(k,893) = rxt(k,517)*y(k,223)
         mat(k,1613) = rxt(k,517)*y(k,153)
         mat(k,825) = .2202005_r8*rxt(k,535)*y(k,136) + .2202005_r8*rxt(k,536) &
                      *y(k,223)
         mat(k,81) = .0023005_r8*rxt(k,537)*y(k,223)
         mat(k,771) = .0031005_r8*rxt(k,540)*y(k,223)
         mat(k,36) = .2381005_r8*rxt(k,541)*y(k,223)
         mat(k,851) = .0508005_r8*rxt(k,543)*y(k,136) + .0508005_r8*rxt(k,544) &
                      *y(k,223)
         mat(k,1524) = .2202005_r8*rxt(k,535)*y(k,6) + .0508005_r8*rxt(k,543)*y(k,110)
         mat(k,42) = .5931005_r8*rxt(k,545)*y(k,223)
         mat(k,107) = .1364005_r8*rxt(k,546)*y(k,223)
         mat(k,135) = .1677005_r8*rxt(k,547)*y(k,223)
         mat(k,1586) = .2202005_r8*rxt(k,536)*y(k,6) + .0023005_r8*rxt(k,537)*y(k,7) &
                      + .0031005_r8*rxt(k,540)*y(k,98) + .2381005_r8*rxt(k,541) &
                      *y(k,104) + .0508005_r8*rxt(k,544)*y(k,110) &
                      + .5931005_r8*rxt(k,545)*y(k,175) + .1364005_r8*rxt(k,546) &
                      *y(k,183) + .1677005_r8*rxt(k,547)*y(k,185)
         mat(k,826) = .2067005_r8*rxt(k,535)*y(k,136) + .2067005_r8*rxt(k,536) &
                      *y(k,223)
         mat(k,82) = .0008005_r8*rxt(k,537)*y(k,223)
         mat(k,772) = .0035005_r8*rxt(k,540)*y(k,223)
         mat(k,37) = .1308005_r8*rxt(k,541)*y(k,223)
         mat(k,852) = .1149005_r8*rxt(k,543)*y(k,136) + .1149005_r8*rxt(k,544) &
                      *y(k,223)
         mat(k,1525) = .2067005_r8*rxt(k,535)*y(k,6) + .1149005_r8*rxt(k,543)*y(k,110)
         mat(k,43) = .1534005_r8*rxt(k,545)*y(k,223)
         mat(k,108) = .0101005_r8*rxt(k,546)*y(k,223)
         mat(k,136) = .0174005_r8*rxt(k,547)*y(k,223)
         mat(k,1587) = .2067005_r8*rxt(k,536)*y(k,6) + .0008005_r8*rxt(k,537)*y(k,7) &
                      + .0035005_r8*rxt(k,540)*y(k,98) + .1308005_r8*rxt(k,541) &
                      *y(k,104) + .1149005_r8*rxt(k,544)*y(k,110) &
                      + .1534005_r8*rxt(k,545)*y(k,175) + .0101005_r8*rxt(k,546) &
                      *y(k,183) + .0174005_r8*rxt(k,547)*y(k,185)
         mat(k,827) = .0653005_r8*rxt(k,535)*y(k,136) + .0653005_r8*rxt(k,536) &
                      *y(k,223)
         mat(k,83) = .0843005_r8*rxt(k,537)*y(k,223)
         mat(k,773) = .0003005_r8*rxt(k,540)*y(k,223)
         mat(k,38) = .0348005_r8*rxt(k,541)*y(k,223)
         mat(k,853) = .0348005_r8*rxt(k,543)*y(k,136) + .0348005_r8*rxt(k,544) &
                      *y(k,223)
         mat(k,1526) = .0653005_r8*rxt(k,535)*y(k,6) + .0348005_r8*rxt(k,543)*y(k,110)
         mat(k,44) = .0459005_r8*rxt(k,545)*y(k,223)
         mat(k,109) = .0763005_r8*rxt(k,546)*y(k,223)
         mat(k,137) = .086_r8*rxt(k,547)*y(k,223)
         mat(k,1588) = .0653005_r8*rxt(k,536)*y(k,6) + .0843005_r8*rxt(k,537)*y(k,7) &
                      + .0003005_r8*rxt(k,540)*y(k,98) + .0348005_r8*rxt(k,541) &
                      *y(k,104) + .0348005_r8*rxt(k,544)*y(k,110) &
                      + .0459005_r8*rxt(k,545)*y(k,175) + .0763005_r8*rxt(k,546) &
                      *y(k,183) + .086_r8*rxt(k,547)*y(k,185)
         mat(k,828) = .1749305_r8*rxt(k,534)*y(k,126) + .1284005_r8*rxt(k,535) &
                      *y(k,136) + .1284005_r8*rxt(k,536)*y(k,223)
         mat(k,84) = .0443005_r8*rxt(k,537)*y(k,223)
         mat(k,774) = .0590245_r8*rxt(k,538)*y(k,126) + .0033005_r8*rxt(k,539) &
                      *y(k,136) + .0271005_r8*rxt(k,540)*y(k,223)
         mat(k,39) = .0076005_r8*rxt(k,541)*y(k,223)
         mat(k,854) = .1749305_r8*rxt(k,542)*y(k,126) + .0554005_r8*rxt(k,543) &
                      *y(k,136) + .0554005_r8*rxt(k,544)*y(k,223)
         mat(k,1440) = .1749305_r8*rxt(k,534)*y(k,6) + .0590245_r8*rxt(k,538)*y(k,98) &
                      + .1749305_r8*rxt(k,542)*y(k,110)
         mat(k,1527) = .1284005_r8*rxt(k,535)*y(k,6) + .0033005_r8*rxt(k,539)*y(k,98) &
                      + .0554005_r8*rxt(k,543)*y(k,110)
         mat(k,45) = .0085005_r8*rxt(k,545)*y(k,223)
         mat(k,110) = .2157005_r8*rxt(k,546)*y(k,223)
         mat(k,138) = .0512005_r8*rxt(k,547)*y(k,223)
         mat(k,1589) = .1284005_r8*rxt(k,536)*y(k,6) + .0443005_r8*rxt(k,537)*y(k,7) &
                      + .0271005_r8*rxt(k,540)*y(k,98) + .0076005_r8*rxt(k,541) &
                      *y(k,104) + .0554005_r8*rxt(k,544)*y(k,110) &
                      + .0085005_r8*rxt(k,545)*y(k,175) + .2157005_r8*rxt(k,546) &
                      *y(k,183) + .0512005_r8*rxt(k,547)*y(k,185)
         mat(k,829) = .5901905_r8*rxt(k,534)*y(k,126) + .114_r8*rxt(k,535)*y(k,136) &
                      + .114_r8*rxt(k,536)*y(k,223)
         mat(k,85) = .1621005_r8*rxt(k,537)*y(k,223)
         mat(k,775) = .0250245_r8*rxt(k,538)*y(k,126) + .0474005_r8*rxt(k,540) &
                      *y(k,223)
         mat(k,40) = .0113005_r8*rxt(k,541)*y(k,223)
         mat(k,855) = .5901905_r8*rxt(k,542)*y(k,126) + .1278005_r8*rxt(k,543) &
                      *y(k,136) + .1278005_r8*rxt(k,544)*y(k,223)
         mat(k,1441) = .5901905_r8*rxt(k,534)*y(k,6) + .0250245_r8*rxt(k,538)*y(k,98) &
                      + .5901905_r8*rxt(k,542)*y(k,110)
         mat(k,1528) = .114_r8*rxt(k,535)*y(k,6) + .1278005_r8*rxt(k,543)*y(k,110)
         mat(k,46) = .0128005_r8*rxt(k,545)*y(k,223)
         mat(k,111) = .0738005_r8*rxt(k,546)*y(k,223)
         mat(k,139) = .1598005_r8*rxt(k,547)*y(k,223)
         mat(k,1590) = .114_r8*rxt(k,536)*y(k,6) + .1621005_r8*rxt(k,537)*y(k,7) &
                      + .0474005_r8*rxt(k,540)*y(k,98) + .0113005_r8*rxt(k,541) &
                      *y(k,104) + .1278005_r8*rxt(k,544)*y(k,110) &
                      + .0128005_r8*rxt(k,545)*y(k,175) + .0738005_r8*rxt(k,546) &
                      *y(k,183) + .1598005_r8*rxt(k,547)*y(k,185)
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
         mat(k,47) = -(rxt(k,545)*y(k,223))
         mat(k,1592) = -rxt(k,545)*y(k,175)
         mat(k,114) = .100_r8*rxt(k,469)*y(k,223)
         mat(k,140) = .230_r8*rxt(k,471)*y(k,223)
         mat(k,1604) = .100_r8*rxt(k,469)*y(k,183) + .230_r8*rxt(k,471)*y(k,185)
         mat(k,485) = -(rxt(k,493)*y(k,223))
         mat(k,1654) = -rxt(k,493)*y(k,177)
         mat(k,2030) = rxt(k,491)*y(k,228)
         mat(k,1036) = rxt(k,491)*y(k,205)
         mat(k,523) = -(rxt(k,494)*y(k,223))
         mat(k,1658) = -rxt(k,494)*y(k,178)
         mat(k,1821) = .200_r8*rxt(k,487)*y(k,218) + .200_r8*rxt(k,497)*y(k,229)
         mat(k,1367) = .500_r8*rxt(k,485)*y(k,218)
         mat(k,1055) = .200_r8*rxt(k,487)*y(k,124) + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1014) = .200_r8*rxt(k,497)*y(k,124)
         mat(k,370) = -(rxt(k,498)*y(k,223))
         mat(k,1641) = -rxt(k,498)*y(k,179)
         mat(k,2023) = rxt(k,496)*y(k,229)
         mat(k,1013) = rxt(k,496)*y(k,205)
         mat(k,918) = -(rxt(k,499)*y(k,126) + rxt(k,500)*y(k,223))
         mat(k,1455) = -rxt(k,499)*y(k,180)
         mat(k,1690) = -rxt(k,500)*y(k,180)
         mat(k,837) = .330_r8*rxt(k,480)*y(k,136)
         mat(k,863) = .330_r8*rxt(k,483)*y(k,136)
         mat(k,1843) = .800_r8*rxt(k,487)*y(k,218) + .800_r8*rxt(k,497)*y(k,229)
         mat(k,1455) = mat(k,1455) + rxt(k,488)*y(k,218)
         mat(k,1547) = .330_r8*rxt(k,480)*y(k,6) + .330_r8*rxt(k,483)*y(k,110)
         mat(k,524) = rxt(k,494)*y(k,223)
         mat(k,1374) = .500_r8*rxt(k,485)*y(k,218) + rxt(k,495)*y(k,229)
         mat(k,1057) = .800_r8*rxt(k,487)*y(k,124) + rxt(k,488)*y(k,126) &
                      + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1690) = mat(k,1690) + rxt(k,494)*y(k,178)
         mat(k,1017) = .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,199)
         mat(k,968) = -(rxt(k,501)*y(k,223))
         mat(k,1694) = -rxt(k,501)*y(k,181)
         mat(k,838) = .300_r8*rxt(k,480)*y(k,136)
         mat(k,864) = .300_r8*rxt(k,483)*y(k,136)
         mat(k,1846) = .900_r8*rxt(k,492)*y(k,228)
         mat(k,1550) = .300_r8*rxt(k,480)*y(k,6) + .300_r8*rxt(k,483)*y(k,110)
         mat(k,1377) = rxt(k,490)*y(k,228)
         mat(k,1040) = .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,199)
         mat(k,496) = -(rxt(k,468)*y(k,223))
         mat(k,1655) = -rxt(k,468)*y(k,182)
         mat(k,2031) = rxt(k,466)*y(k,230)
         mat(k,617) = rxt(k,466)*y(k,205)
         mat(k,112) = -(rxt(k,469)*y(k,223))
         mat(k,1602) = -rxt(k,469)*y(k,183)
         mat(k,128) = -(rxt(k,435)*y(k,223))
         mat(k,1605) = -rxt(k,435)*y(k,184)
         mat(k,2002) = rxt(k,432)*y(k,231)
         mat(k,1100) = rxt(k,432)*y(k,205)
         mat(k,141) = -(rxt(k,471)*y(k,223))
         mat(k,1607) = -rxt(k,471)*y(k,185)
         mat(k,598) = -(rxt(k,474)*y(k,223))
         mat(k,1665) = -rxt(k,474)*y(k,186)
         mat(k,2038) = rxt(k,472)*y(k,232)
         mat(k,634) = rxt(k,472)*y(k,205)
         mat(k,149) = -(rxt(k,477)*y(k,223))
         mat(k,1608) = -rxt(k,477)*y(k,187)
         mat(k,142) = .150_r8*rxt(k,471)*y(k,223)
         mat(k,1608) = mat(k,1608) + .150_r8*rxt(k,471)*y(k,185)
         mat(k,317) = -(rxt(k,478)*y(k,223))
         mat(k,1634) = -rxt(k,478)*y(k,188)
         mat(k,2015) = rxt(k,475)*y(k,233)
         mat(k,393) = rxt(k,475)*y(k,205)
         mat(k,413) = -(rxt(k,436)*y(k,205) + rxt(k,437)*y(k,124) + rxt(k,465) &
                      *y(k,125))
         mat(k,2027) = -rxt(k,436)*y(k,191)
         mat(k,1816) = -rxt(k,437)*y(k,191)
         mat(k,1964) = -rxt(k,465)*y(k,191)
         mat(k,172) = rxt(k,442)*y(k,223)
         mat(k,1646) = rxt(k,442)*y(k,22)
         mat(k,882) = -(rxt(k,397)*y(k,205) + (rxt(k,398) + rxt(k,399)) * y(k,124))
         mat(k,2054) = -rxt(k,397)*y(k,192)
         mat(k,1841) = -(rxt(k,398) + rxt(k,399)) * y(k,192)
         mat(k,564) = rxt(k,400)*y(k,223)
         mat(k,163) = rxt(k,401)*y(k,223)
         mat(k,1686) = rxt(k,400)*y(k,2) + rxt(k,401)*y(k,15)
         mat(k,379) = -(rxt(k,439)*y(k,205) + rxt(k,440)*y(k,124))
         mat(k,2024) = -rxt(k,439)*y(k,193)
         mat(k,1813) = -rxt(k,440)*y(k,193)
         mat(k,89) = .350_r8*rxt(k,438)*y(k,223)
         mat(k,277) = rxt(k,441)*y(k,223)
         mat(k,1642) = .350_r8*rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8)
         mat(k,325) = -(rxt(k,443)*y(k,205) + rxt(k,445)*y(k,124))
         mat(k,2016) = -rxt(k,443)*y(k,194)
         mat(k,1807) = -rxt(k,445)*y(k,194)
         mat(k,241) = rxt(k,444)*y(k,223)
         mat(k,115) = .070_r8*rxt(k,469)*y(k,223)
         mat(k,143) = .060_r8*rxt(k,471)*y(k,223)
         mat(k,1635) = rxt(k,444)*y(k,23) + .070_r8*rxt(k,469)*y(k,183) &
                      + .060_r8*rxt(k,471)*y(k,185)
         mat(k,802) = -(4._r8*rxt(k,320)*y(k,195) + rxt(k,321)*y(k,199) + rxt(k,322) &
                      *y(k,205) + rxt(k,323)*y(k,124))
         mat(k,1371) = -rxt(k,321)*y(k,195)
         mat(k,2051) = -rxt(k,322)*y(k,195)
         mat(k,1838) = -rxt(k,323)*y(k,195)
         mat(k,249) = .500_r8*rxt(k,325)*y(k,223)
         mat(k,209) = rxt(k,326)*y(k,56) + rxt(k,327)*y(k,223)
         mat(k,1932) = rxt(k,326)*y(k,28)
         mat(k,1681) = .500_r8*rxt(k,325)*y(k,27) + rxt(k,327)*y(k,28)
         mat(k,744) = -(rxt(k,349)*y(k,199) + rxt(k,350)*y(k,205) + rxt(k,351) &
                      *y(k,124))
         mat(k,1369) = -rxt(k,349)*y(k,196)
         mat(k,2047) = -rxt(k,350)*y(k,196)
         mat(k,1836) = -rxt(k,351)*y(k,196)
         mat(k,312) = rxt(k,352)*y(k,223)
         mat(k,58) = rxt(k,353)*y(k,223)
         mat(k,1676) = rxt(k,352)*y(k,30) + rxt(k,353)*y(k,31)
         mat(k,531) = -(rxt(k,446)*y(k,205) + rxt(k,447)*y(k,124))
         mat(k,2034) = -rxt(k,446)*y(k,197)
         mat(k,1822) = -rxt(k,447)*y(k,197)
         mat(k,185) = rxt(k,448)*y(k,223)
         mat(k,1822) = mat(k,1822) + rxt(k,437)*y(k,191)
         mat(k,1536) = rxt(k,463)*y(k,143)
         mat(k,360) = rxt(k,463)*y(k,136)
         mat(k,414) = rxt(k,437)*y(k,124) + .400_r8*rxt(k,436)*y(k,205)
         mat(k,2034) = mat(k,2034) + .400_r8*rxt(k,436)*y(k,191)
         mat(k,1659) = rxt(k,448)*y(k,32)
         mat(k,1301) = -(4._r8*rxt(k,331)*y(k,198) + rxt(k,332)*y(k,199) + rxt(k,333) &
                      *y(k,205) + rxt(k,334)*y(k,124) + rxt(k,345)*y(k,125) + rxt(k,372) &
                      *y(k,209) + rxt(k,405)*y(k,207) + rxt(k,410)*y(k,208) + rxt(k,419) &
                      *y(k,101) + rxt(k,430)*y(k,231))
         mat(k,1394) = -rxt(k,332)*y(k,198)
         mat(k,2076) = -rxt(k,333)*y(k,198)
         mat(k,1864) = -rxt(k,334)*y(k,198)
         mat(k,1980) = -rxt(k,345)*y(k,198)
         mat(k,1232) = -rxt(k,372)*y(k,198)
         mat(k,1210) = -rxt(k,405)*y(k,198)
         mat(k,1183) = -rxt(k,410)*y(k,198)
         mat(k,1130) = -rxt(k,419)*y(k,198)
         mat(k,1108) = -rxt(k,430)*y(k,198)
         mat(k,844) = .060_r8*rxt(k,480)*y(k,136)
         mat(k,979) = rxt(k,328)*y(k,126) + rxt(k,329)*y(k,223)
         mat(k,1155) = rxt(k,354)*y(k,126) + rxt(k,355)*y(k,223)
         mat(k,407) = .500_r8*rxt(k,336)*y(k,223)
         mat(k,785) = .080_r8*rxt(k,425)*y(k,136)
         mat(k,1146) = .100_r8*rxt(k,378)*y(k,136)
         mat(k,870) = .060_r8*rxt(k,483)*y(k,136)
         mat(k,1252) = .280_r8*rxt(k,392)*y(k,136)
         mat(k,1864) = mat(k,1864) + .530_r8*rxt(k,376)*y(k,209) + rxt(k,385)*y(k,211) &
                      + rxt(k,388)*y(k,213) + rxt(k,363)*y(k,227)
         mat(k,1477) = rxt(k,328)*y(k,45) + rxt(k,354)*y(k,49) + .530_r8*rxt(k,375) &
                      *y(k,209) + rxt(k,386)*y(k,211)
         mat(k,1566) = .060_r8*rxt(k,480)*y(k,6) + .080_r8*rxt(k,425)*y(k,98) &
                      + .100_r8*rxt(k,378)*y(k,105) + .060_r8*rxt(k,483)*y(k,110) &
                      + .280_r8*rxt(k,392)*y(k,111)
         mat(k,971) = .650_r8*rxt(k,501)*y(k,223)
         mat(k,1301) = mat(k,1301) + .530_r8*rxt(k,372)*y(k,209)
         mat(k,1394) = mat(k,1394) + .260_r8*rxt(k,373)*y(k,209) + rxt(k,382)*y(k,211) &
                      + .300_r8*rxt(k,361)*y(k,227)
         mat(k,2076) = mat(k,2076) + .450_r8*rxt(k,383)*y(k,211) + .200_r8*rxt(k,387) &
                      *y(k,213) + .150_r8*rxt(k,362)*y(k,227)
         mat(k,1232) = mat(k,1232) + .530_r8*rxt(k,376)*y(k,124) + .530_r8*rxt(k,375) &
                      *y(k,126) + .530_r8*rxt(k,372)*y(k,198) + .260_r8*rxt(k,373) &
                      *y(k,199)
         mat(k,1271) = rxt(k,385)*y(k,124) + rxt(k,386)*y(k,126) + rxt(k,382)*y(k,199) &
                      + .450_r8*rxt(k,383)*y(k,205) + 4.000_r8*rxt(k,384)*y(k,211)
         mat(k,574) = rxt(k,388)*y(k,124) + .200_r8*rxt(k,387)*y(k,205)
         mat(k,1713) = rxt(k,329)*y(k,45) + rxt(k,355)*y(k,49) + .500_r8*rxt(k,336) &
                      *y(k,51) + .650_r8*rxt(k,501)*y(k,181)
         mat(k,1079) = rxt(k,363)*y(k,124) + .300_r8*rxt(k,361)*y(k,199) &
                      + .150_r8*rxt(k,362)*y(k,205)
         mat(k,1395) = -(rxt(k,221)*y(k,59) + (4._r8*rxt(k,298) + 4._r8*rxt(k,299) &
                      ) * y(k,199) + rxt(k,300)*y(k,205) + rxt(k,301)*y(k,124) &
                      + rxt(k,321)*y(k,195) + rxt(k,332)*y(k,198) + rxt(k,349) &
                      *y(k,196) + rxt(k,361)*y(k,227) + rxt(k,373)*y(k,209) + rxt(k,382) &
                      *y(k,211) + rxt(k,406)*y(k,207) + rxt(k,411)*y(k,208) + rxt(k,420) &
                      *y(k,101) + rxt(k,431)*y(k,231) + rxt(k,485)*y(k,218) + rxt(k,490) &
                      *y(k,228) + rxt(k,495)*y(k,229))
         mat(k,1507) = -rxt(k,221)*y(k,199)
         mat(k,2079) = -rxt(k,300)*y(k,199)
         mat(k,1866) = -rxt(k,301)*y(k,199)
         mat(k,804) = -rxt(k,321)*y(k,199)
         mat(k,1302) = -rxt(k,332)*y(k,199)
         mat(k,747) = -rxt(k,349)*y(k,199)
         mat(k,1080) = -rxt(k,361)*y(k,199)
         mat(k,1233) = -rxt(k,373)*y(k,199)
         mat(k,1272) = -rxt(k,382)*y(k,199)
         mat(k,1211) = -rxt(k,406)*y(k,199)
         mat(k,1184) = -rxt(k,411)*y(k,199)
         mat(k,1131) = -rxt(k,420)*y(k,199)
         mat(k,1109) = -rxt(k,431)*y(k,199)
         mat(k,1064) = -rxt(k,485)*y(k,199)
         mat(k,1045) = -rxt(k,490)*y(k,199)
         mat(k,1025) = -rxt(k,495)*y(k,199)
         mat(k,936) = .280_r8*rxt(k,348)*y(k,136)
         mat(k,449) = rxt(k,335)*y(k,223)
         mat(k,282) = .700_r8*rxt(k,303)*y(k,223)
         mat(k,786) = .050_r8*rxt(k,425)*y(k,136)
         mat(k,1131) = mat(k,1131) + rxt(k,419)*y(k,198)
         mat(k,1866) = mat(k,1866) + rxt(k,334)*y(k,198) + .830_r8*rxt(k,451)*y(k,200) &
                      + .170_r8*rxt(k,457)*y(k,212)
         mat(k,1569) = .280_r8*rxt(k,348)*y(k,29) + .050_r8*rxt(k,425)*y(k,98)
         mat(k,1302) = mat(k,1302) + rxt(k,419)*y(k,101) + rxt(k,334)*y(k,124) &
                      + 4.000_r8*rxt(k,331)*y(k,198) + .900_r8*rxt(k,332)*y(k,199) &
                      + .450_r8*rxt(k,333)*y(k,205) + rxt(k,405)*y(k,207) + rxt(k,410) &
                      *y(k,208) + rxt(k,372)*y(k,209) + rxt(k,381)*y(k,211) &
                      + rxt(k,430)*y(k,231)
         mat(k,1395) = mat(k,1395) + .900_r8*rxt(k,332)*y(k,198)
         mat(k,650) = .830_r8*rxt(k,451)*y(k,124) + .330_r8*rxt(k,450)*y(k,205)
         mat(k,2079) = mat(k,2079) + .450_r8*rxt(k,333)*y(k,198) + .330_r8*rxt(k,450) &
                      *y(k,200) + .070_r8*rxt(k,456)*y(k,212)
         mat(k,1211) = mat(k,1211) + rxt(k,405)*y(k,198)
         mat(k,1184) = mat(k,1184) + rxt(k,410)*y(k,198)
         mat(k,1233) = mat(k,1233) + rxt(k,372)*y(k,198)
         mat(k,1272) = mat(k,1272) + rxt(k,381)*y(k,198)
         mat(k,820) = .170_r8*rxt(k,457)*y(k,124) + .070_r8*rxt(k,456)*y(k,205)
         mat(k,1717) = rxt(k,335)*y(k,50) + .700_r8*rxt(k,303)*y(k,53)
         mat(k,1109) = mat(k,1109) + rxt(k,430)*y(k,198)
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
         mat(k,647) = -(rxt(k,450)*y(k,205) + rxt(k,451)*y(k,124) + rxt(k,452) &
                      *y(k,125))
         mat(k,2042) = -rxt(k,450)*y(k,200)
         mat(k,1828) = -rxt(k,451)*y(k,200)
         mat(k,1969) = -rxt(k,452)*y(k,200)
         mat(k,694) = -(rxt(k,565)*y(k,216) + rxt(k,566)*y(k,222) + rxt(k,567) &
                      *y(k,215))
         mat(k,684) = -rxt(k,565)*y(k,201)
         mat(k,676) = -rxt(k,566)*y(k,201)
         mat(k,553) = -rxt(k,567)*y(k,201)
         mat(k,468) = -((rxt(k,369) + rxt(k,370)) * y(k,124))
         mat(k,1818) = -(rxt(k,369) + rxt(k,370)) * y(k,202)
         mat(k,262) = rxt(k,368)*y(k,223)
         mat(k,1652) = rxt(k,368)*y(k,16)
         mat(k,347) = -(rxt(k,340)*y(k,135))
         mat(k,1411) = -rxt(k,340)*y(k,203)
         mat(k,1811) = .750_r8*rxt(k,338)*y(k,204)
         mat(k,665) = .750_r8*rxt(k,338)*y(k,124)
         mat(k,666) = -(rxt(k,337)*y(k,205) + rxt(k,338)*y(k,124))
         mat(k,2043) = -rxt(k,337)*y(k,204)
         mat(k,1829) = -rxt(k,338)*y(k,204)
         mat(k,453) = rxt(k,344)*y(k,223)
         mat(k,1671) = rxt(k,344)*y(k,25)
         mat(k,2092) = -((rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,76) + rxt(k,176) &
                      *y(k,134) + rxt(k,177)*y(k,136) + rxt(k,181)*y(k,223) &
                      + 4._r8*rxt(k,186)*y(k,205) + rxt(k,198)*y(k,126) + rxt(k,203) &
                      *y(k,124) + rxt(k,208)*y(k,125) + (rxt(k,218) + rxt(k,219) &
                      ) * y(k,56) + rxt(k,225)*y(k,59) + rxt(k,251)*y(k,17) + rxt(k,257) &
                      *y(k,19) + rxt(k,294)*y(k,42) + rxt(k,300)*y(k,199) + rxt(k,308) &
                      *y(k,206) + rxt(k,322)*y(k,195) + rxt(k,333)*y(k,198) + rxt(k,337) &
                      *y(k,204) + rxt(k,350)*y(k,196) + rxt(k,358)*y(k,226) + rxt(k,362) &
                      *y(k,227) + rxt(k,374)*y(k,209) + rxt(k,383)*y(k,211) + rxt(k,387) &
                      *y(k,213) + rxt(k,397)*y(k,192) + rxt(k,407)*y(k,207) + rxt(k,412) &
                      *y(k,208) + rxt(k,421)*y(k,101) + rxt(k,432)*y(k,231) + rxt(k,436) &
                      *y(k,191) + rxt(k,439)*y(k,193) + rxt(k,443)*y(k,194) + rxt(k,446) &
                      *y(k,197) + rxt(k,450)*y(k,200) + rxt(k,453)*y(k,210) + rxt(k,456) &
                      *y(k,212) + rxt(k,459)*y(k,225) + rxt(k,466)*y(k,230) + rxt(k,472) &
                      *y(k,232) + rxt(k,475)*y(k,233) + rxt(k,486)*y(k,218) + rxt(k,491) &
                      *y(k,228) + rxt(k,496)*y(k,229))
         mat(k,1772) = -(rxt(k,172) + rxt(k,173) + rxt(k,174)) * y(k,205)
         mat(k,1920) = -rxt(k,176)*y(k,205)
         mat(k,1582) = -rxt(k,177)*y(k,205)
         mat(k,1730) = -rxt(k,181)*y(k,205)
         mat(k,1493) = -rxt(k,198)*y(k,205)
         mat(k,1879) = -rxt(k,203)*y(k,205)
         mat(k,1997) = -rxt(k,208)*y(k,205)
         mat(k,1955) = -(rxt(k,218) + rxt(k,219)) * y(k,205)
         mat(k,1520) = -rxt(k,225)*y(k,205)
         mat(k,1358) = -rxt(k,251)*y(k,205)
         mat(k,2140) = -rxt(k,257)*y(k,205)
         mat(k,2116) = -rxt(k,294)*y(k,205)
         mat(k,1406) = -rxt(k,300)*y(k,205)
         mat(k,335) = -rxt(k,308)*y(k,205)
         mat(k,809) = -rxt(k,322)*y(k,205)
         mat(k,1311) = -rxt(k,333)*y(k,205)
         mat(k,672) = -rxt(k,337)*y(k,205)
         mat(k,752) = -rxt(k,350)*y(k,205)
         mat(k,726) = -rxt(k,358)*y(k,205)
         mat(k,1084) = -rxt(k,362)*y(k,205)
         mat(k,1240) = -rxt(k,374)*y(k,205)
         mat(k,1280) = -rxt(k,383)*y(k,205)
         mat(k,578) = -rxt(k,387)*y(k,205)
         mat(k,891) = -rxt(k,397)*y(k,205)
         mat(k,1219) = -rxt(k,407)*y(k,205)
         mat(k,1192) = -rxt(k,412)*y(k,205)
         mat(k,1138) = -rxt(k,421)*y(k,205)
         mat(k,1115) = -rxt(k,432)*y(k,205)
         mat(k,418) = -rxt(k,436)*y(k,205)
         mat(k,385) = -rxt(k,439)*y(k,205)
         mat(k,330) = -rxt(k,443)*y(k,205)
         mat(k,536) = -rxt(k,446)*y(k,205)
         mat(k,654) = -rxt(k,450)*y(k,205)
         mat(k,614) = -rxt(k,453)*y(k,205)
         mat(k,824) = -rxt(k,456)*y(k,205)
         mat(k,343) = -rxt(k,459)*y(k,205)
         mat(k,629) = -rxt(k,466)*y(k,205)
         mat(k,646) = -rxt(k,472)*y(k,205)
         mat(k,400) = -rxt(k,475)*y(k,205)
         mat(k,1071) = -rxt(k,486)*y(k,205)
         mat(k,1051) = -rxt(k,491)*y(k,205)
         mat(k,1032) = -rxt(k,496)*y(k,205)
         mat(k,848) = .570_r8*rxt(k,480)*y(k,136)
         mat(k,91) = .650_r8*rxt(k,438)*y(k,223)
         mat(k,1358) = mat(k,1358) + rxt(k,250)*y(k,42)
         mat(k,2140) = mat(k,2140) + rxt(k,262)*y(k,223)
         mat(k,207) = .350_r8*rxt(k,317)*y(k,223)
         mat(k,458) = .130_r8*rxt(k,319)*y(k,136)
         mat(k,182) = rxt(k,324)*y(k,223)
         mat(k,943) = .280_r8*rxt(k,348)*y(k,136)
         mat(k,2116) = mat(k,2116) + rxt(k,250)*y(k,17) + rxt(k,214)*y(k,56) &
                      + rxt(k,295)*y(k,126) + rxt(k,296)*y(k,134)
         mat(k,56) = rxt(k,330)*y(k,223)
         mat(k,712) = rxt(k,302)*y(k,223)
         mat(k,1955) = mat(k,1955) + rxt(k,214)*y(k,42) + rxt(k,217)*y(k,79)
         mat(k,1520) = mat(k,1520) + rxt(k,221)*y(k,199) + rxt(k,232)*y(k,223)
         mat(k,1006) = rxt(k,305)*y(k,223)
         mat(k,123) = .730_r8*rxt(k,449)*y(k,223)
         mat(k,198) = .500_r8*rxt(k,516)*y(k,223)
         mat(k,965) = rxt(k,341)*y(k,223)
         mat(k,815) = rxt(k,342)*y(k,223)
         mat(k,1772) = mat(k,1772) + rxt(k,175)*y(k,135)
         mat(k,510) = rxt(k,217)*y(k,56) + rxt(k,171)*y(k,134) + rxt(k,180)*y(k,223)
         mat(k,106) = rxt(k,306)*y(k,223)
         mat(k,716) = rxt(k,307)*y(k,223)
         mat(k,1000) = rxt(k,371)*y(k,223)
         mat(k,1011) = rxt(k,356)*y(k,223)
         mat(k,790) = .370_r8*rxt(k,425)*y(k,136)
         mat(k,483) = .300_r8*rxt(k,416)*y(k,223)
         mat(k,431) = rxt(k,417)*y(k,223)
         mat(k,1138) = mat(k,1138) + rxt(k,422)*y(k,124) + rxt(k,423)*y(k,126) &
                      + rxt(k,419)*y(k,198) + 1.200_r8*rxt(k,420)*y(k,199)
         mat(k,304) = rxt(k,424)*y(k,223)
         mat(k,1150) = .140_r8*rxt(k,378)*y(k,136)
         mat(k,218) = .200_r8*rxt(k,380)*y(k,223)
         mat(k,466) = .500_r8*rxt(k,391)*y(k,223)
         mat(k,874) = .570_r8*rxt(k,483)*y(k,136)
         mat(k,1262) = .280_r8*rxt(k,392)*y(k,136)
         mat(k,274) = rxt(k,428)*y(k,223)
         mat(k,959) = rxt(k,429)*y(k,223)
         mat(k,1879) = mat(k,1879) + rxt(k,422)*y(k,101) + rxt(k,398)*y(k,192) &
                      + rxt(k,440)*y(k,193) + rxt(k,445)*y(k,194) + rxt(k,323) &
                      *y(k,195) + rxt(k,351)*y(k,196) + rxt(k,301)*y(k,199) &
                      + .170_r8*rxt(k,451)*y(k,200) + rxt(k,369)*y(k,202) &
                      + .250_r8*rxt(k,338)*y(k,204) + rxt(k,310)*y(k,206) &
                      + .920_r8*rxt(k,408)*y(k,207) + .920_r8*rxt(k,414)*y(k,208) &
                      + .470_r8*rxt(k,376)*y(k,209) + .400_r8*rxt(k,454)*y(k,210) &
                      + .830_r8*rxt(k,457)*y(k,212) + rxt(k,460)*y(k,225) + rxt(k,359) &
                      *y(k,226) + .900_r8*rxt(k,492)*y(k,228) + .800_r8*rxt(k,497) &
                      *y(k,229) + rxt(k,467)*y(k,230) + rxt(k,433)*y(k,231) &
                      + rxt(k,473)*y(k,232) + rxt(k,476)*y(k,233)
         mat(k,1493) = mat(k,1493) + rxt(k,295)*y(k,42) + rxt(k,423)*y(k,101) &
                      + rxt(k,409)*y(k,207) + rxt(k,415)*y(k,208) + .470_r8*rxt(k,375) &
                      *y(k,209) + rxt(k,201)*y(k,223) + rxt(k,434)*y(k,231)
         mat(k,1920) = mat(k,1920) + rxt(k,296)*y(k,42) + rxt(k,171)*y(k,79)
         mat(k,1436) = rxt(k,175)*y(k,76) + rxt(k,340)*y(k,203)
         mat(k,1582) = mat(k,1582) + .570_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,319) &
                      *y(k,25) + .280_r8*rxt(k,348)*y(k,29) + .370_r8*rxt(k,425) &
                      *y(k,98) + .140_r8*rxt(k,378)*y(k,105) + .570_r8*rxt(k,483) &
                      *y(k,110) + .280_r8*rxt(k,392)*y(k,111) + rxt(k,183)*y(k,223)
         mat(k,100) = .800_r8*rxt(k,461)*y(k,223)
         mat(k,898) = rxt(k,517)*y(k,223)
         mat(k,974) = .200_r8*rxt(k,501)*y(k,223)
         mat(k,118) = .280_r8*rxt(k,469)*y(k,223)
         mat(k,148) = .380_r8*rxt(k,471)*y(k,223)
         mat(k,153) = .630_r8*rxt(k,477)*y(k,223)
         mat(k,891) = mat(k,891) + rxt(k,398)*y(k,124)
         mat(k,385) = mat(k,385) + rxt(k,440)*y(k,124)
         mat(k,330) = mat(k,330) + rxt(k,445)*y(k,124)
         mat(k,809) = mat(k,809) + rxt(k,323)*y(k,124) + 2.400_r8*rxt(k,320)*y(k,195) &
                      + rxt(k,321)*y(k,199)
         mat(k,752) = mat(k,752) + rxt(k,351)*y(k,124) + rxt(k,349)*y(k,199)
         mat(k,1311) = mat(k,1311) + rxt(k,419)*y(k,101) + .900_r8*rxt(k,332)*y(k,199) &
                      + rxt(k,405)*y(k,207) + rxt(k,410)*y(k,208) + .470_r8*rxt(k,372) &
                      *y(k,209) + rxt(k,430)*y(k,231)
         mat(k,1406) = mat(k,1406) + rxt(k,221)*y(k,59) + 1.200_r8*rxt(k,420)*y(k,101) &
                      + rxt(k,301)*y(k,124) + rxt(k,321)*y(k,195) + rxt(k,349) &
                      *y(k,196) + .900_r8*rxt(k,332)*y(k,198) + 4.000_r8*rxt(k,298) &
                      *y(k,199) + rxt(k,406)*y(k,207) + rxt(k,411)*y(k,208) &
                      + .730_r8*rxt(k,373)*y(k,209) + rxt(k,382)*y(k,211) &
                      + .500_r8*rxt(k,485)*y(k,218) + .300_r8*rxt(k,361)*y(k,227) &
                      + rxt(k,490)*y(k,228) + rxt(k,495)*y(k,229) + .800_r8*rxt(k,431) &
                      *y(k,231)
         mat(k,654) = mat(k,654) + .170_r8*rxt(k,451)*y(k,124) + .070_r8*rxt(k,450) &
                      *y(k,205)
         mat(k,474) = rxt(k,369)*y(k,124)
         mat(k,350) = rxt(k,340)*y(k,135)
         mat(k,672) = mat(k,672) + .250_r8*rxt(k,338)*y(k,124)
         mat(k,2092) = mat(k,2092) + .070_r8*rxt(k,450)*y(k,200) + .160_r8*rxt(k,453) &
                      *y(k,210) + .330_r8*rxt(k,456)*y(k,212)
         mat(k,335) = mat(k,335) + rxt(k,310)*y(k,124)
         mat(k,1219) = mat(k,1219) + .920_r8*rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126) &
                      + rxt(k,405)*y(k,198) + rxt(k,406)*y(k,199)
         mat(k,1192) = mat(k,1192) + .920_r8*rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126) &
                      + rxt(k,410)*y(k,198) + rxt(k,411)*y(k,199)
         mat(k,1240) = mat(k,1240) + .470_r8*rxt(k,376)*y(k,124) + .470_r8*rxt(k,375) &
                      *y(k,126) + .470_r8*rxt(k,372)*y(k,198) + .730_r8*rxt(k,373) &
                      *y(k,199)
         mat(k,614) = mat(k,614) + .400_r8*rxt(k,454)*y(k,124) + .160_r8*rxt(k,453) &
                      *y(k,205)
         mat(k,1280) = mat(k,1280) + rxt(k,382)*y(k,199)
         mat(k,824) = mat(k,824) + .830_r8*rxt(k,457)*y(k,124) + .330_r8*rxt(k,456) &
                      *y(k,205)
         mat(k,1071) = mat(k,1071) + .500_r8*rxt(k,485)*y(k,199)
         mat(k,1730) = mat(k,1730) + .650_r8*rxt(k,438)*y(k,7) + rxt(k,262)*y(k,19) &
                      + .350_r8*rxt(k,317)*y(k,24) + rxt(k,324)*y(k,26) + rxt(k,330) &
                      *y(k,47) + rxt(k,302)*y(k,52) + rxt(k,232)*y(k,59) + rxt(k,305) &
                      *y(k,62) + .730_r8*rxt(k,449)*y(k,66) + .500_r8*rxt(k,516) &
                      *y(k,67) + rxt(k,341)*y(k,74) + rxt(k,342)*y(k,75) + rxt(k,180) &
                      *y(k,79) + rxt(k,306)*y(k,86) + rxt(k,307)*y(k,87) + rxt(k,371) &
                      *y(k,93) + rxt(k,356)*y(k,95) + .300_r8*rxt(k,416)*y(k,99) &
                      + rxt(k,417)*y(k,100) + rxt(k,424)*y(k,102) + .200_r8*rxt(k,380) &
                      *y(k,106) + .500_r8*rxt(k,391)*y(k,109) + rxt(k,428)*y(k,115) &
                      + rxt(k,429)*y(k,116) + rxt(k,201)*y(k,126) + rxt(k,183) &
                      *y(k,136) + .800_r8*rxt(k,461)*y(k,144) + rxt(k,517)*y(k,153) &
                      + .200_r8*rxt(k,501)*y(k,181) + .280_r8*rxt(k,469)*y(k,183) &
                      + .380_r8*rxt(k,471)*y(k,185) + .630_r8*rxt(k,477)*y(k,187)
         mat(k,343) = mat(k,343) + rxt(k,460)*y(k,124)
         mat(k,726) = mat(k,726) + rxt(k,359)*y(k,124)
         mat(k,1084) = mat(k,1084) + .300_r8*rxt(k,361)*y(k,199)
         mat(k,1051) = mat(k,1051) + .900_r8*rxt(k,492)*y(k,124) + rxt(k,490)*y(k,199)
         mat(k,1032) = mat(k,1032) + .800_r8*rxt(k,497)*y(k,124) + rxt(k,495)*y(k,199)
         mat(k,629) = mat(k,629) + rxt(k,467)*y(k,124)
         mat(k,1115) = mat(k,1115) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126) &
                      + rxt(k,430)*y(k,198) + .800_r8*rxt(k,431)*y(k,199)
         mat(k,646) = mat(k,646) + rxt(k,473)*y(k,124)
         mat(k,400) = mat(k,400) + rxt(k,476)*y(k,124)
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
         mat(k,331) = -(rxt(k,308)*y(k,205) + rxt(k,310)*y(k,124))
         mat(k,2017) = -rxt(k,308)*y(k,206)
         mat(k,1808) = -rxt(k,310)*y(k,206)
         mat(k,2096) = rxt(k,294)*y(k,205)
         mat(k,2017) = mat(k,2017) + rxt(k,294)*y(k,42)
         mat(k,1206) = -(rxt(k,405)*y(k,198) + rxt(k,406)*y(k,199) + rxt(k,407) &
                      *y(k,205) + rxt(k,408)*y(k,124) + rxt(k,409)*y(k,126))
         mat(k,1297) = -rxt(k,405)*y(k,207)
         mat(k,1390) = -rxt(k,406)*y(k,207)
         mat(k,2072) = -rxt(k,407)*y(k,207)
         mat(k,1860) = -rxt(k,408)*y(k,207)
         mat(k,1473) = -rxt(k,409)*y(k,207)
         mat(k,783) = .600_r8*rxt(k,426)*y(k,223)
         mat(k,1709) = .600_r8*rxt(k,426)*y(k,98)
         mat(k,1179) = -(rxt(k,410)*y(k,198) + rxt(k,411)*y(k,199) + rxt(k,412) &
                      *y(k,205) + rxt(k,414)*y(k,124) + rxt(k,415)*y(k,126))
         mat(k,1296) = -rxt(k,410)*y(k,208)
         mat(k,1389) = -rxt(k,411)*y(k,208)
         mat(k,2071) = -rxt(k,412)*y(k,208)
         mat(k,1859) = -rxt(k,414)*y(k,208)
         mat(k,1472) = -rxt(k,415)*y(k,208)
         mat(k,782) = .400_r8*rxt(k,426)*y(k,223)
         mat(k,1708) = .400_r8*rxt(k,426)*y(k,98)
         mat(k,1230) = -(rxt(k,372)*y(k,198) + rxt(k,373)*y(k,199) + rxt(k,374) &
                      *y(k,205) + rxt(k,375)*y(k,126) + (rxt(k,376) + rxt(k,377) &
                      ) * y(k,124))
         mat(k,1298) = -rxt(k,372)*y(k,209)
         mat(k,1391) = -rxt(k,373)*y(k,209)
         mat(k,2073) = -rxt(k,374)*y(k,209)
         mat(k,1474) = -rxt(k,375)*y(k,209)
         mat(k,1861) = -(rxt(k,376) + rxt(k,377)) * y(k,209)
         mat(k,1144) = .500_r8*rxt(k,379)*y(k,223)
         mat(k,215) = .200_r8*rxt(k,380)*y(k,223)
         mat(k,1249) = rxt(k,393)*y(k,223)
         mat(k,1710) = .500_r8*rxt(k,379)*y(k,105) + .200_r8*rxt(k,380)*y(k,106) &
                      + rxt(k,393)*y(k,111)
         mat(k,609) = -(rxt(k,453)*y(k,205) + rxt(k,454)*y(k,124) + rxt(k,455) &
                      *y(k,125))
         mat(k,2039) = -rxt(k,453)*y(k,210)
         mat(k,1825) = -rxt(k,454)*y(k,210)
         mat(k,1968) = -rxt(k,455)*y(k,210)
         mat(k,1270) = -(rxt(k,381)*y(k,198) + rxt(k,382)*y(k,199) + rxt(k,383) &
                      *y(k,205) + 4._r8*rxt(k,384)*y(k,211) + rxt(k,385)*y(k,124) &
                      + rxt(k,386)*y(k,126) + rxt(k,394)*y(k,125))
         mat(k,1300) = -rxt(k,381)*y(k,211)
         mat(k,1393) = -rxt(k,382)*y(k,211)
         mat(k,2075) = -rxt(k,383)*y(k,211)
         mat(k,1863) = -rxt(k,385)*y(k,211)
         mat(k,1476) = -rxt(k,386)*y(k,211)
         mat(k,1979) = -rxt(k,394)*y(k,211)
         mat(k,1145) = .500_r8*rxt(k,379)*y(k,223)
         mat(k,216) = .500_r8*rxt(k,380)*y(k,223)
         mat(k,1712) = .500_r8*rxt(k,379)*y(k,105) + .500_r8*rxt(k,380)*y(k,106)
         mat(k,817) = -(rxt(k,456)*y(k,205) + rxt(k,457)*y(k,124) + rxt(k,458) &
                      *y(k,125))
         mat(k,2053) = -rxt(k,456)*y(k,212)
         mat(k,1840) = -rxt(k,457)*y(k,212)
         mat(k,1973) = -rxt(k,458)*y(k,212)
         mat(k,572) = -(rxt(k,387)*y(k,205) + rxt(k,388)*y(k,124))
         mat(k,2036) = -rxt(k,387)*y(k,213)
         mat(k,1824) = -rxt(k,388)*y(k,213)
         mat(k,402) = rxt(k,389)*y(k,223)
         mat(k,220) = rxt(k,390)*y(k,223)
         mat(k,1662) = rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108)
         mat(k,419) = -(rxt(k,188)*y(k,134) + rxt(k,189)*y(k,135))
         mat(k,1887) = -rxt(k,188)*y(k,214)
         mat(k,1413) = -rxt(k,189)*y(k,214)
         mat(k,1887) = mat(k,1887) + rxt(k,569)*y(k,215)
         mat(k,690) = .900_r8*rxt(k,567)*y(k,215) + .800_r8*rxt(k,565)*y(k,216)
         mat(k,548) = rxt(k,569)*y(k,134) + .900_r8*rxt(k,567)*y(k,201)
         mat(k,682) = .800_r8*rxt(k,565)*y(k,201)
         mat(k,550) = -(rxt(k,567)*y(k,201) + rxt(k,568)*y(k,135) + (rxt(k,569) &
                      + rxt(k,570)) * y(k,134))
         mat(k,691) = -rxt(k,567)*y(k,215)
         mat(k,1415) = -rxt(k,568)*y(k,215)
         mat(k,1891) = -(rxt(k,569) + rxt(k,570)) * y(k,215)
         mat(k,683) = -(rxt(k,565)*y(k,201))
         mat(k,693) = -rxt(k,565)*y(k,216)
         mat(k,732) = rxt(k,574)*y(k,222)
         mat(k,1831) = rxt(k,576)*y(k,222)
         mat(k,1894) = rxt(k,569)*y(k,215)
         mat(k,1418) = rxt(k,573)*y(k,217)
         mat(k,552) = rxt(k,569)*y(k,134)
         mat(k,389) = rxt(k,573)*y(k,135)
         mat(k,675) = rxt(k,574)*y(k,112) + rxt(k,576)*y(k,124)
         mat(k,386) = -(rxt(k,571)*y(k,134) + (rxt(k,572) + rxt(k,573)) * y(k,135))
         mat(k,1886) = -rxt(k,571)*y(k,217)
         mat(k,1412) = -(rxt(k,572) + rxt(k,573)) * y(k,217)
         mat(k,1061) = -(rxt(k,485)*y(k,199) + rxt(k,486)*y(k,205) + rxt(k,487) &
                      *y(k,124) + rxt(k,488)*y(k,126))
         mat(k,1383) = -rxt(k,485)*y(k,218)
         mat(k,2064) = -rxt(k,486)*y(k,218)
         mat(k,1853) = -rxt(k,487)*y(k,218)
         mat(k,1466) = -rxt(k,488)*y(k,218)
         mat(k,841) = rxt(k,479)*y(k,126)
         mat(k,867) = rxt(k,482)*y(k,126)
         mat(k,1466) = mat(k,1466) + rxt(k,479)*y(k,6) + rxt(k,482)*y(k,110) &
                      + .500_r8*rxt(k,499)*y(k,180)
         mat(k,295) = rxt(k,489)*y(k,223)
         mat(k,922) = .500_r8*rxt(k,499)*y(k,126)
         mat(k,1701) = rxt(k,489)*y(k,128)
         mat(k,1793) = -(rxt(k,153)*y(k,77) + rxt(k,154)*y(k,234) + (rxt(k,156) &
                      + rxt(k,157)) * y(k,135) + rxt(k,158)*y(k,136) + (rxt(k,246) &
                      + rxt(k,247)) * y(k,85) + (rxt(k,269) + rxt(k,270)) * y(k,81) &
                      + rxt(k,275)*y(k,64) + rxt(k,276)*y(k,65) + rxt(k,314)*y(k,86))
         mat(k,1096) = -rxt(k,153)*y(k,219)
         mat(k,2161) = -rxt(k,154)*y(k,219)
         mat(k,1431) = -(rxt(k,156) + rxt(k,157)) * y(k,219)
         mat(k,1577) = -rxt(k,158)*y(k,219)
         mat(k,1343) = -(rxt(k,246) + rxt(k,247)) * y(k,219)
         mat(k,767) = -(rxt(k,269) + rxt(k,270)) * y(k,219)
         mat(k,63) = -rxt(k,275)*y(k,219)
         mat(k,133) = -rxt(k,276)*y(k,219)
         mat(k,105) = -rxt(k,314)*y(k,219)
         mat(k,1431) = mat(k,1431) + rxt(k,189)*y(k,214)
         mat(k,698) = .850_r8*rxt(k,566)*y(k,222)
         mat(k,422) = rxt(k,189)*y(k,135)
         mat(k,679) = .850_r8*rxt(k,566)*y(k,201)
         mat(k,76) = -(rxt(k,160)*y(k,134) + rxt(k,161)*y(k,135))
         mat(k,1883) = -rxt(k,160)*y(k,220)
         mat(k,1409) = -rxt(k,161)*y(k,220)
         mat(k,1883) = mat(k,1883) + rxt(k,164)*y(k,221)
         mat(k,1409) = mat(k,1409) + rxt(k,165)*y(k,221)
         mat(k,1529) = rxt(k,166)*y(k,221)
         mat(k,78) = rxt(k,164)*y(k,134) + rxt(k,165)*y(k,135) + rxt(k,166)*y(k,136)
         mat(k,79) = -(rxt(k,164)*y(k,134) + rxt(k,165)*y(k,135) + rxt(k,166)*y(k,136))
         mat(k,1884) = -rxt(k,164)*y(k,221)
         mat(k,1410) = -rxt(k,165)*y(k,221)
         mat(k,1530) = -rxt(k,166)*y(k,221)
         mat(k,1410) = mat(k,1410) + rxt(k,156)*y(k,219)
         mat(k,1777) = rxt(k,156)*y(k,135)
         mat(k,674) = -(rxt(k,566)*y(k,201) + rxt(k,574)*y(k,112) + rxt(k,576) &
                      *y(k,124))
         mat(k,692) = -rxt(k,566)*y(k,222)
         mat(k,731) = -rxt(k,574)*y(k,222)
         mat(k,1830) = -rxt(k,576)*y(k,222)
         mat(k,1417) = rxt(k,568)*y(k,215) + rxt(k,572)*y(k,217) + rxt(k,579)*y(k,224)
         mat(k,551) = rxt(k,568)*y(k,135)
         mat(k,388) = rxt(k,572)*y(k,135)
         mat(k,513) = rxt(k,579)*y(k,135)
         mat(k,1722) = -(rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) + rxt(k,181)*y(k,205) &
                      + rxt(k,182)*y(k,134) + rxt(k,183)*y(k,136) + (4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185)) * y(k,223) + rxt(k,187)*y(k,90) + rxt(k,201) &
                      *y(k,126) + rxt(k,202)*y(k,112) + rxt(k,210)*y(k,125) + rxt(k,211) &
                      *y(k,89) + rxt(k,230)*y(k,60) + (rxt(k,232) + rxt(k,233) &
                      ) * y(k,59) + rxt(k,235)*y(k,85) + rxt(k,238)*y(k,92) + rxt(k,262) &
                      *y(k,19) + rxt(k,264)*y(k,81) + rxt(k,297)*y(k,42) + rxt(k,302) &
                      *y(k,52) + rxt(k,303)*y(k,53) + (rxt(k,305) + rxt(k,315) &
                      ) * y(k,62) + rxt(k,306)*y(k,86) + rxt(k,307)*y(k,87) + rxt(k,317) &
                      *y(k,24) + rxt(k,324)*y(k,26) + rxt(k,325)*y(k,27) + rxt(k,327) &
                      *y(k,28) + rxt(k,329)*y(k,45) + rxt(k,330)*y(k,47) + rxt(k,335) &
                      *y(k,50) + rxt(k,336)*y(k,51) + rxt(k,341)*y(k,74) + rxt(k,342) &
                      *y(k,75) + rxt(k,343)*y(k,141) + rxt(k,344)*y(k,25) + rxt(k,352) &
                      *y(k,30) + rxt(k,353)*y(k,31) + rxt(k,355)*y(k,49) + rxt(k,356) &
                      *y(k,95) + rxt(k,357)*y(k,127) + rxt(k,360)*y(k,148) + rxt(k,364) &
                      *y(k,149) + rxt(k,365)*y(k,29) + rxt(k,366)*y(k,48) + rxt(k,368) &
                      *y(k,16) + rxt(k,371)*y(k,93) + rxt(k,379)*y(k,105) + rxt(k,380) &
                      *y(k,106) + rxt(k,389)*y(k,107) + rxt(k,390)*y(k,108) + rxt(k,391) &
                      *y(k,109) + rxt(k,393)*y(k,111) + rxt(k,396)*y(k,1) + rxt(k,400) &
                      *y(k,2) + rxt(k,401)*y(k,15) + rxt(k,402)*y(k,94) + rxt(k,403) &
                      *y(k,96) + rxt(k,404)*y(k,97) + rxt(k,416)*y(k,99) + rxt(k,417) &
                      *y(k,100) + rxt(k,424)*y(k,102) + rxt(k,426)*y(k,98) + rxt(k,427) &
                      *y(k,103) + rxt(k,428)*y(k,115) + rxt(k,429)*y(k,116) + rxt(k,435) &
                      *y(k,184) + rxt(k,438)*y(k,7) + rxt(k,441)*y(k,8) + rxt(k,442) &
                      *y(k,22) + rxt(k,444)*y(k,23) + rxt(k,448)*y(k,32) + rxt(k,449) &
                      *y(k,66) + rxt(k,461)*y(k,144) + rxt(k,464)*y(k,145) + rxt(k,468) &
                      *y(k,182) + rxt(k,469)*y(k,183) + rxt(k,471)*y(k,185) + rxt(k,474) &
                      *y(k,186) + rxt(k,477)*y(k,187) + rxt(k,478)*y(k,188) + rxt(k,481) &
                      *y(k,6) + rxt(k,484)*y(k,110) + rxt(k,489)*y(k,128) + rxt(k,493) &
                      *y(k,177) + rxt(k,494)*y(k,178) + rxt(k,498)*y(k,179) + rxt(k,500) &
                      *y(k,180) + rxt(k,501)*y(k,181) + (rxt(k,503) + rxt(k,516) &
                      ) * y(k,67) + rxt(k,505)*y(k,139) + rxt(k,510)*y(k,150) &
                      + rxt(k,515)*y(k,152) + rxt(k,517)*y(k,153) + rxt(k,519) &
                      *y(k,120))
         mat(k,1093) = -rxt(k,179)*y(k,223)
         mat(k,507) = -rxt(k,180)*y(k,223)
         mat(k,2084) = -rxt(k,181)*y(k,223)
         mat(k,1912) = -rxt(k,182)*y(k,223)
         mat(k,1574) = -rxt(k,183)*y(k,223)
         mat(k,366) = -rxt(k,187)*y(k,223)
         mat(k,1485) = -rxt(k,201)*y(k,223)
         mat(k,737) = -rxt(k,202)*y(k,223)
         mat(k,1989) = -rxt(k,210)*y(k,223)
         mat(k,1744) = -rxt(k,211)*y(k,223)
         mat(k,912) = -rxt(k,230)*y(k,223)
         mat(k,1512) = -(rxt(k,232) + rxt(k,233)) * y(k,223)
         mat(k,1340) = -rxt(k,235)*y(k,223)
         mat(k,759) = -rxt(k,238)*y(k,223)
         mat(k,2132) = -rxt(k,262)*y(k,223)
         mat(k,765) = -rxt(k,264)*y(k,223)
         mat(k,2108) = -rxt(k,297)*y(k,223)
         mat(k,711) = -rxt(k,302)*y(k,223)
         mat(k,283) = -rxt(k,303)*y(k,223)
         mat(k,1004) = -(rxt(k,305) + rxt(k,315)) * y(k,223)
         mat(k,104) = -rxt(k,306)*y(k,223)
         mat(k,715) = -rxt(k,307)*y(k,223)
         mat(k,206) = -rxt(k,317)*y(k,223)
         mat(k,181) = -rxt(k,324)*y(k,223)
         mat(k,251) = -rxt(k,325)*y(k,223)
         mat(k,211) = -rxt(k,327)*y(k,223)
         mat(k,982) = -rxt(k,329)*y(k,223)
         mat(k,55) = -rxt(k,330)*y(k,223)
         mat(k,450) = -rxt(k,335)*y(k,223)
         mat(k,409) = -rxt(k,336)*y(k,223)
         mat(k,964) = -rxt(k,341)*y(k,223)
         mat(k,814) = -rxt(k,342)*y(k,223)
         mat(k,356) = -rxt(k,343)*y(k,223)
         mat(k,457) = -rxt(k,344)*y(k,223)
         mat(k,314) = -rxt(k,352)*y(k,223)
         mat(k,59) = -rxt(k,353)*y(k,223)
         mat(k,1157) = -rxt(k,355)*y(k,223)
         mat(k,1010) = -rxt(k,356)*y(k,223)
         mat(k,796) = -rxt(k,357)*y(k,223)
         mat(k,436) = -rxt(k,360)*y(k,223)
         mat(k,308) = -rxt(k,364)*y(k,223)
         mat(k,940) = -rxt(k,365)*y(k,223)
         mat(k,904) = -rxt(k,366)*y(k,223)
         mat(k,266) = -rxt(k,368)*y(k,223)
         mat(k,997) = -rxt(k,371)*y(k,223)
         mat(k,1148) = -rxt(k,379)*y(k,223)
         mat(k,217) = -rxt(k,380)*y(k,223)
         mat(k,405) = -rxt(k,389)*y(k,223)
         mat(k,223) = -rxt(k,390)*y(k,223)
         mat(k,464) = -rxt(k,391)*y(k,223)
         mat(k,1257) = -rxt(k,393)*y(k,223)
         mat(k,544) = -rxt(k,396)*y(k,223)
         mat(k,568) = -rxt(k,400)*y(k,223)
         mat(k,164) = -rxt(k,401)*y(k,223)
         mat(k,160) = -rxt(k,402)*y(k,223)
         mat(k,247) = -rxt(k,403)*y(k,223)
         mat(k,72) = -rxt(k,404)*y(k,223)
         mat(k,480) = -rxt(k,416)*y(k,223)
         mat(k,429) = -rxt(k,417)*y(k,223)
         mat(k,302) = -rxt(k,424)*y(k,223)
         mat(k,789) = -rxt(k,426)*y(k,223)
         mat(k,584) = -rxt(k,427)*y(k,223)
         mat(k,272) = -rxt(k,428)*y(k,223)
         mat(k,955) = -rxt(k,429)*y(k,223)
         mat(k,130) = -rxt(k,435)*y(k,223)
         mat(k,90) = -rxt(k,438)*y(k,223)
         mat(k,279) = -rxt(k,441)*y(k,223)
         mat(k,173) = -rxt(k,442)*y(k,223)
         mat(k,243) = -rxt(k,444)*y(k,223)
         mat(k,186) = -rxt(k,448)*y(k,223)
         mat(k,122) = -rxt(k,449)*y(k,223)
         mat(k,99) = -rxt(k,461)*y(k,223)
         mat(k,237) = -rxt(k,464)*y(k,223)
         mat(k,503) = -rxt(k,468)*y(k,223)
         mat(k,117) = -rxt(k,469)*y(k,223)
         mat(k,147) = -rxt(k,471)*y(k,223)
         mat(k,607) = -rxt(k,474)*y(k,223)
         mat(k,152) = -rxt(k,477)*y(k,223)
         mat(k,321) = -rxt(k,478)*y(k,223)
         mat(k,847) = -rxt(k,481)*y(k,223)
         mat(k,873) = -rxt(k,484)*y(k,223)
         mat(k,296) = -rxt(k,489)*y(k,223)
         mat(k,491) = -rxt(k,493)*y(k,223)
         mat(k,525) = -rxt(k,494)*y(k,223)
         mat(k,374) = -rxt(k,498)*y(k,223)
         mat(k,924) = -rxt(k,500)*y(k,223)
         mat(k,973) = -rxt(k,501)*y(k,223)
         mat(k,196) = -(rxt(k,503) + rxt(k,516)) * y(k,223)
         mat(k,258) = -rxt(k,505)*y(k,223)
         mat(k,591) = -rxt(k,510)*y(k,223)
         mat(k,1322) = -rxt(k,515)*y(k,223)
         mat(k,896) = -rxt(k,517)*y(k,223)
         mat(k,52) = -rxt(k,519)*y(k,223)
         mat(k,847) = mat(k,847) + .630_r8*rxt(k,480)*y(k,136)
         mat(k,206) = mat(k,206) + .650_r8*rxt(k,317)*y(k,223)
         mat(k,457) = mat(k,457) + .130_r8*rxt(k,319)*y(k,136)
         mat(k,251) = mat(k,251) + .500_r8*rxt(k,325)*y(k,223)
         mat(k,940) = mat(k,940) + .360_r8*rxt(k,348)*y(k,136)
         mat(k,2108) = mat(k,2108) + rxt(k,296)*y(k,134)
         mat(k,283) = mat(k,283) + .300_r8*rxt(k,303)*y(k,223)
         mat(k,1947) = rxt(k,219)*y(k,205)
         mat(k,660) = rxt(k,273)*y(k,234)
         mat(k,1764) = rxt(k,178)*y(k,136) + 2.000_r8*rxt(k,173)*y(k,205)
         mat(k,1093) = mat(k,1093) + rxt(k,170)*y(k,134) + rxt(k,153)*y(k,219)
         mat(k,507) = mat(k,507) + rxt(k,171)*y(k,134)
         mat(k,765) = mat(k,765) + rxt(k,263)*y(k,134) + rxt(k,269)*y(k,219)
         mat(k,1340) = mat(k,1340) + rxt(k,234)*y(k,134) + rxt(k,246)*y(k,219)
         mat(k,104) = mat(k,104) + rxt(k,314)*y(k,219)
         mat(k,705) = rxt(k,265)*y(k,134)
         mat(k,759) = mat(k,759) + rxt(k,237)*y(k,134)
         mat(k,789) = mat(k,789) + .320_r8*rxt(k,425)*y(k,136)
         mat(k,584) = mat(k,584) + .600_r8*rxt(k,427)*y(k,223)
         mat(k,1148) = mat(k,1148) + .240_r8*rxt(k,378)*y(k,136)
         mat(k,217) = mat(k,217) + .100_r8*rxt(k,380)*y(k,223)
         mat(k,873) = mat(k,873) + .630_r8*rxt(k,483)*y(k,136)
         mat(k,1257) = mat(k,1257) + .360_r8*rxt(k,392)*y(k,136)
         mat(k,1871) = rxt(k,203)*y(k,205)
         mat(k,1485) = mat(k,1485) + rxt(k,198)*y(k,205)
         mat(k,1912) = mat(k,1912) + rxt(k,296)*y(k,42) + rxt(k,170)*y(k,77) &
                      + rxt(k,171)*y(k,79) + rxt(k,263)*y(k,81) + rxt(k,234)*y(k,85) &
                      + rxt(k,265)*y(k,91) + rxt(k,237)*y(k,92) + rxt(k,176)*y(k,205)
         mat(k,1574) = mat(k,1574) + .630_r8*rxt(k,480)*y(k,6) + .130_r8*rxt(k,319) &
                      *y(k,25) + .360_r8*rxt(k,348)*y(k,29) + rxt(k,178)*y(k,76) &
                      + .320_r8*rxt(k,425)*y(k,98) + .240_r8*rxt(k,378)*y(k,105) &
                      + .630_r8*rxt(k,483)*y(k,110) + .360_r8*rxt(k,392)*y(k,111) &
                      + rxt(k,177)*y(k,205)
         mat(k,436) = mat(k,436) + .500_r8*rxt(k,360)*y(k,223)
         mat(k,130) = mat(k,130) + .500_r8*rxt(k,435)*y(k,223)
         mat(k,415) = .400_r8*rxt(k,436)*y(k,205)
         mat(k,1306) = .450_r8*rxt(k,333)*y(k,205)
         mat(k,651) = .400_r8*rxt(k,450)*y(k,205)
         mat(k,2084) = mat(k,2084) + rxt(k,219)*y(k,56) + 2.000_r8*rxt(k,173)*y(k,76) &
                      + rxt(k,203)*y(k,124) + rxt(k,198)*y(k,126) + rxt(k,176) &
                      *y(k,134) + rxt(k,177)*y(k,136) + .400_r8*rxt(k,436)*y(k,191) &
                      + .450_r8*rxt(k,333)*y(k,198) + .400_r8*rxt(k,450)*y(k,200) &
                      + .450_r8*rxt(k,383)*y(k,211) + .400_r8*rxt(k,456)*y(k,212) &
                      + .200_r8*rxt(k,387)*y(k,213) + .150_r8*rxt(k,362)*y(k,227)
         mat(k,1275) = .450_r8*rxt(k,383)*y(k,205)
         mat(k,821) = .400_r8*rxt(k,456)*y(k,205)
         mat(k,575) = .200_r8*rxt(k,387)*y(k,205)
         mat(k,1790) = rxt(k,153)*y(k,77) + rxt(k,269)*y(k,81) + rxt(k,246)*y(k,85) &
                      + rxt(k,314)*y(k,86) + 2.000_r8*rxt(k,154)*y(k,234)
         mat(k,1722) = mat(k,1722) + .650_r8*rxt(k,317)*y(k,24) + .500_r8*rxt(k,325) &
                      *y(k,27) + .300_r8*rxt(k,303)*y(k,53) + .600_r8*rxt(k,427) &
                      *y(k,103) + .100_r8*rxt(k,380)*y(k,106) + .500_r8*rxt(k,360) &
                      *y(k,148) + .500_r8*rxt(k,435)*y(k,184)
         mat(k,1081) = .150_r8*rxt(k,362)*y(k,205)
         mat(k,2158) = rxt(k,273)*y(k,73) + 2.000_r8*rxt(k,154)*y(k,219)
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
         mat(k,512) = -(rxt(k,579)*y(k,135))
         mat(k,1414) = -rxt(k,579)*y(k,224)
         mat(k,1890) = rxt(k,570)*y(k,215) + rxt(k,571)*y(k,217)
         mat(k,549) = rxt(k,570)*y(k,134)
         mat(k,387) = rxt(k,571)*y(k,134)
         mat(k,338) = -(rxt(k,459)*y(k,205) + rxt(k,460)*y(k,124))
         mat(k,2018) = -rxt(k,459)*y(k,225)
         mat(k,1809) = -rxt(k,460)*y(k,225)
         mat(k,120) = .200_r8*rxt(k,449)*y(k,223)
         mat(k,97) = .140_r8*rxt(k,461)*y(k,223)
         mat(k,235) = rxt(k,464)*y(k,223)
         mat(k,1636) = .200_r8*rxt(k,449)*y(k,66) + .140_r8*rxt(k,461)*y(k,144) &
                      + rxt(k,464)*y(k,145)
         mat(k,719) = -(rxt(k,358)*y(k,205) + rxt(k,359)*y(k,124))
         mat(k,2046) = -rxt(k,358)*y(k,226)
         mat(k,1834) = -rxt(k,359)*y(k,226)
         mat(k,928) = rxt(k,365)*y(k,223)
         mat(k,433) = .500_r8*rxt(k,360)*y(k,223)
         mat(k,1674) = rxt(k,365)*y(k,29) + .500_r8*rxt(k,360)*y(k,148)
         mat(k,1077) = -(rxt(k,361)*y(k,199) + rxt(k,362)*y(k,205) + rxt(k,363) &
                      *y(k,124))
         mat(k,1384) = -rxt(k,361)*y(k,227)
         mat(k,2065) = -rxt(k,362)*y(k,227)
         mat(k,1854) = -rxt(k,363)*y(k,227)
         mat(k,842) = .060_r8*rxt(k,480)*y(k,136)
         mat(k,901) = rxt(k,366)*y(k,223)
         mat(k,868) = .060_r8*rxt(k,483)*y(k,136)
         mat(k,1557) = .060_r8*rxt(k,480)*y(k,6) + .060_r8*rxt(k,483)*y(k,110)
         mat(k,306) = rxt(k,364)*y(k,223)
         mat(k,970) = .150_r8*rxt(k,501)*y(k,223)
         mat(k,1702) = rxt(k,366)*y(k,48) + rxt(k,364)*y(k,149) + .150_r8*rxt(k,501) &
                      *y(k,181)
         mat(k,1042) = -(rxt(k,490)*y(k,199) + rxt(k,491)*y(k,205) + rxt(k,492) &
                      *y(k,124))
         mat(k,1382) = -rxt(k,490)*y(k,228)
         mat(k,2063) = -rxt(k,491)*y(k,228)
         mat(k,1852) = -rxt(k,492)*y(k,228)
         mat(k,1465) = .500_r8*rxt(k,499)*y(k,180)
         mat(k,490) = rxt(k,493)*y(k,223)
         mat(k,921) = .500_r8*rxt(k,499)*y(k,126) + rxt(k,500)*y(k,223)
         mat(k,1700) = rxt(k,493)*y(k,177) + rxt(k,500)*y(k,180)
         mat(k,1020) = -(rxt(k,495)*y(k,199) + rxt(k,496)*y(k,205) + rxt(k,497) &
                      *y(k,124))
         mat(k,1381) = -rxt(k,495)*y(k,229)
         mat(k,2062) = -rxt(k,496)*y(k,229)
         mat(k,1851) = -rxt(k,497)*y(k,229)
         mat(k,840) = rxt(k,481)*y(k,223)
         mat(k,866) = rxt(k,484)*y(k,223)
         mat(k,373) = rxt(k,498)*y(k,223)
         mat(k,1699) = rxt(k,481)*y(k,6) + rxt(k,484)*y(k,110) + rxt(k,498)*y(k,179)
         mat(k,620) = -(rxt(k,466)*y(k,205) + rxt(k,467)*y(k,124))
         mat(k,2040) = -rxt(k,466)*y(k,230)
         mat(k,1826) = -rxt(k,467)*y(k,230)
         mat(k,499) = rxt(k,468)*y(k,223)
         mat(k,116) = .650_r8*rxt(k,469)*y(k,223)
         mat(k,1667) = rxt(k,468)*y(k,182) + .650_r8*rxt(k,469)*y(k,183)
         mat(k,1106) = -(rxt(k,430)*y(k,198) + rxt(k,431)*y(k,199) + rxt(k,432) &
                      *y(k,205) + rxt(k,433)*y(k,124) + rxt(k,434)*y(k,126))
         mat(k,1292) = -rxt(k,430)*y(k,231)
         mat(k,1385) = -rxt(k,431)*y(k,231)
         mat(k,2067) = -rxt(k,432)*y(k,231)
         mat(k,1855) = -rxt(k,433)*y(k,231)
         mat(k,1468) = -rxt(k,434)*y(k,231)
         mat(k,159) = rxt(k,402)*y(k,223)
         mat(k,246) = rxt(k,403)*y(k,223)
         mat(k,71) = rxt(k,404)*y(k,223)
         mat(k,581) = .400_r8*rxt(k,427)*y(k,223)
         mat(k,129) = .500_r8*rxt(k,435)*y(k,223)
         mat(k,1704) = rxt(k,402)*y(k,94) + rxt(k,403)*y(k,96) + rxt(k,404)*y(k,97) &
                      + .400_r8*rxt(k,427)*y(k,103) + .500_r8*rxt(k,435)*y(k,184)
         mat(k,636) = -(rxt(k,472)*y(k,205) + rxt(k,473)*y(k,124))
         mat(k,2041) = -rxt(k,472)*y(k,232)
         mat(k,1827) = -rxt(k,473)*y(k,232)
         mat(k,144) = .560_r8*rxt(k,471)*y(k,223)
         mat(k,600) = rxt(k,474)*y(k,223)
         mat(k,1668) = .560_r8*rxt(k,471)*y(k,185) + rxt(k,474)*y(k,186)
         mat(k,394) = -(rxt(k,475)*y(k,205) + rxt(k,476)*y(k,124))
         mat(k,2025) = -rxt(k,475)*y(k,233)
         mat(k,1814) = -rxt(k,476)*y(k,233)
         mat(k,151) = .300_r8*rxt(k,477)*y(k,223)
         mat(k,318) = rxt(k,478)*y(k,223)
         mat(k,1643) = .300_r8*rxt(k,477)*y(k,187) + rxt(k,478)*y(k,188)
         mat(k,2169) = -(rxt(k,154)*y(k,219) + rxt(k,273)*y(k,73) + rxt(k,518) &
                      *y(k,154))
         mat(k,1801) = -rxt(k,154)*y(k,234)
         mat(k,663) = -rxt(k,273)*y(k,234)
         mat(k,178) = -rxt(k,518)*y(k,234)
         mat(k,213) = rxt(k,327)*y(k,223)
         mat(k,316) = rxt(k,352)*y(k,223)
         mat(k,60) = rxt(k,353)*y(k,223)
         mat(k,2119) = rxt(k,297)*y(k,223)
         mat(k,985) = rxt(k,329)*y(k,223)
         mat(k,905) = rxt(k,366)*y(k,223)
         mat(k,1161) = rxt(k,355)*y(k,223)
         mat(k,451) = rxt(k,335)*y(k,223)
         mat(k,411) = rxt(k,336)*y(k,223)
         mat(k,286) = rxt(k,303)*y(k,223)
         mat(k,1775) = rxt(k,174)*y(k,205)
         mat(k,1099) = rxt(k,179)*y(k,223)
         mat(k,511) = rxt(k,180)*y(k,223)
         mat(k,770) = rxt(k,264)*y(k,223)
         mat(k,1348) = (rxt(k,556)+rxt(k,561))*y(k,91) + (rxt(k,549)+rxt(k,555) &
                       +rxt(k,560))*y(k,92) + rxt(k,235)*y(k,223)
         mat(k,717) = rxt(k,307)*y(k,223)
         mat(k,1755) = rxt(k,211)*y(k,223)
         mat(k,369) = rxt(k,187)*y(k,223)
         mat(k,709) = (rxt(k,556)+rxt(k,561))*y(k,85)
         mat(k,762) = (rxt(k,549)+rxt(k,555)+rxt(k,560))*y(k,85) + rxt(k,238)*y(k,223)
         mat(k,1152) = .500_r8*rxt(k,379)*y(k,223)
         mat(k,53) = rxt(k,519)*y(k,223)
         mat(k,439) = rxt(k,360)*y(k,223)
         mat(k,310) = rxt(k,364)*y(k,223)
         mat(k,2095) = rxt(k,174)*y(k,76) + rxt(k,181)*y(k,223)
         mat(k,1733) = rxt(k,327)*y(k,28) + rxt(k,352)*y(k,30) + rxt(k,353)*y(k,31) &
                      + rxt(k,297)*y(k,42) + rxt(k,329)*y(k,45) + rxt(k,366)*y(k,48) &
                      + rxt(k,355)*y(k,49) + rxt(k,335)*y(k,50) + rxt(k,336)*y(k,51) &
                      + rxt(k,303)*y(k,53) + rxt(k,179)*y(k,77) + rxt(k,180)*y(k,79) &
                      + rxt(k,264)*y(k,81) + rxt(k,235)*y(k,85) + rxt(k,307)*y(k,87) &
                      + rxt(k,211)*y(k,89) + rxt(k,187)*y(k,90) + rxt(k,238)*y(k,92) &
                      + .500_r8*rxt(k,379)*y(k,105) + rxt(k,519)*y(k,120) + rxt(k,360) &
                      *y(k,148) + rxt(k,364)*y(k,149) + rxt(k,181)*y(k,205) &
                      + 2.000_r8*rxt(k,184)*y(k,223)
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
         mat(k, 34) = lmat(k, 34)
         mat(k, 35) = lmat(k, 35)
         mat(k, 41) = mat(k, 41) + lmat(k, 41)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 48) = lmat(k, 48)
         mat(k, 49) = lmat(k, 49)
         mat(k, 50) = lmat(k, 50)
         mat(k, 51) = mat(k, 51) + lmat(k, 51)
         mat(k, 54) = mat(k, 54) + lmat(k, 54)
         mat(k, 57) = mat(k, 57) + lmat(k, 57)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 62) = mat(k, 62) + lmat(k, 62)
         mat(k, 64) = lmat(k, 64)
         mat(k, 65) = lmat(k, 65)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = lmat(k, 69)
         mat(k, 70) = mat(k, 70) + lmat(k, 70)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = lmat(k, 74)
         mat(k, 75) = lmat(k, 75)
         mat(k, 76) = mat(k, 76) + lmat(k, 76)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = lmat(k, 80)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = mat(k, 96) + lmat(k, 96)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = lmat(k, 127)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 134) = mat(k, 134) + lmat(k, 134)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 154) = lmat(k, 154)
         mat(k, 155) = lmat(k, 155)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = lmat(k, 158)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 169) = lmat(k, 169)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 179) = mat(k, 179) + lmat(k, 179)
         mat(k, 183) = mat(k, 183) + lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = lmat(k, 189)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 199) = lmat(k, 199)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = mat(k, 202) + lmat(k, 202)
         mat(k, 208) = mat(k, 208) + lmat(k, 208)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 219) = mat(k, 219) + lmat(k, 219)
         mat(k, 221) = lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = lmat(k, 227)
         mat(k, 228) = lmat(k, 228)
         mat(k, 229) = mat(k, 229) + lmat(k, 229)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = mat(k, 234) + lmat(k, 234)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 238) = lmat(k, 238)
         mat(k, 239) = lmat(k, 239)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 244) = lmat(k, 244)
         mat(k, 245) = mat(k, 245) + lmat(k, 245)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = lmat(k, 254)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 273) = lmat(k, 273)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = mat(k, 279) + lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 287) = lmat(k, 287)
         mat(k, 288) = lmat(k, 288)
         mat(k, 289) = lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 297) = lmat(k, 297)
         mat(k, 298) = lmat(k, 298)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 300) = lmat(k, 300)
         mat(k, 303) = lmat(k, 303)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 307) = lmat(k, 307)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 313) = lmat(k, 313)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 315) = lmat(k, 315)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 338) = mat(k, 338) + lmat(k, 338)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 357) = lmat(k, 357)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 365) = lmat(k, 365)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 368) = lmat(k, 368)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 375) = lmat(k, 375)
         mat(k, 376) = lmat(k, 376)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 394) = mat(k, 394) + lmat(k, 394)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 403) = lmat(k, 403)
         mat(k, 404) = lmat(k, 404)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 430) = lmat(k, 430)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 434) = lmat(k, 434)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 437) = lmat(k, 437)
         mat(k, 438) = lmat(k, 438)
         mat(k, 440) = mat(k, 440) + lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 442) = lmat(k, 442)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 444) = lmat(k, 444)
         mat(k, 446) = lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 452) = mat(k, 452) + lmat(k, 452)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 462) = lmat(k, 462)
         mat(k, 465) = lmat(k, 465)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 481) = lmat(k, 481)
         mat(k, 485) = mat(k, 485) + lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 491) = mat(k, 491) + lmat(k, 491)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 495) = lmat(k, 495)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 497) = lmat(k, 497)
         mat(k, 501) = lmat(k, 501)
         mat(k, 502) = lmat(k, 502)
         mat(k, 503) = mat(k, 503) + lmat(k, 503)
         mat(k, 504) = lmat(k, 504)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 514) = lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 519) = lmat(k, 519)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = lmat(k, 522)
         mat(k, 523) = mat(k, 523) + lmat(k, 523)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = mat(k, 527) + lmat(k, 527)
         mat(k, 528) = lmat(k, 528)
         mat(k, 531) = mat(k, 531) + lmat(k, 531)
         mat(k, 537) = lmat(k, 537)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 545) = mat(k, 545) + lmat(k, 545)
         mat(k, 546) = lmat(k, 546)
         mat(k, 547) = mat(k, 547) + lmat(k, 547)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = lmat(k, 566)
         mat(k, 568) = mat(k, 568) + lmat(k, 568)
         mat(k, 569) = lmat(k, 569)
         mat(k, 570) = lmat(k, 570)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 580) = mat(k, 580) + lmat(k, 580)
         mat(k, 582) = lmat(k, 582)
         mat(k, 583) = lmat(k, 583)
         mat(k, 585) = lmat(k, 585)
         mat(k, 586) = lmat(k, 586)
         mat(k, 587) = mat(k, 587) + lmat(k, 587)
         mat(k, 594) = lmat(k, 594)
         mat(k, 595) = lmat(k, 595)
         mat(k, 596) = lmat(k, 596)
         mat(k, 597) = lmat(k, 597)
         mat(k, 598) = mat(k, 598) + lmat(k, 598)
         mat(k, 602) = lmat(k, 602)
         mat(k, 605) = lmat(k, 605)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 608) = lmat(k, 608)
         mat(k, 609) = mat(k, 609) + lmat(k, 609)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 636) = mat(k, 636) + lmat(k, 636)
         mat(k, 647) = mat(k, 647) + lmat(k, 647)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 658) = lmat(k, 658)
         mat(k, 666) = mat(k, 666) + lmat(k, 666)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 680) = mat(k, 680) + lmat(k, 680)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 702) = mat(k, 702) + lmat(k, 702)
         mat(k, 704) = lmat(k, 704)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 710) = mat(k, 710) + lmat(k, 710)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 719) = mat(k, 719) + lmat(k, 719)
         mat(k, 729) = lmat(k, 729)
         mat(k, 733) = lmat(k, 733)
         mat(k, 734) = mat(k, 734) + lmat(k, 734)
         mat(k, 744) = mat(k, 744) + lmat(k, 744)
         mat(k, 756) = mat(k, 756) + lmat(k, 756)
         mat(k, 759) = mat(k, 759) + lmat(k, 759)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 763) = mat(k, 763) + lmat(k, 763)
         mat(k, 764) = mat(k, 764) + lmat(k, 764)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 777) = mat(k, 777) + lmat(k, 777)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 795) = lmat(k, 795)
         mat(k, 797) = mat(k, 797) + lmat(k, 797)
         mat(k, 798) = lmat(k, 798)
         mat(k, 802) = mat(k, 802) + lmat(k, 802)
         mat(k, 811) = lmat(k, 811)
         mat(k, 812) = mat(k, 812) + lmat(k, 812)
         mat(k, 813) = mat(k, 813) + lmat(k, 813)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 817) = mat(k, 817) + lmat(k, 817)
         mat(k, 834) = mat(k, 834) + lmat(k, 834)
         mat(k, 860) = mat(k, 860) + lmat(k, 860)
         mat(k, 882) = mat(k, 882) + lmat(k, 882)
         mat(k, 894) = mat(k, 894) + lmat(k, 894)
         mat(k, 895) = lmat(k, 895)
         mat(k, 897) = lmat(k, 897)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 902) = lmat(k, 902)
         mat(k, 903) = lmat(k, 903)
         mat(k, 907) = mat(k, 907) + lmat(k, 907)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 910) = mat(k, 910) + lmat(k, 910)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 915) = mat(k, 915) + lmat(k, 915)
         mat(k, 916) = lmat(k, 916)
         mat(k, 918) = mat(k, 918) + lmat(k, 918)
         mat(k, 919) = lmat(k, 919)
         mat(k, 920) = lmat(k, 920)
         mat(k, 925) = lmat(k, 925)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 946) = lmat(k, 946)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 956) = lmat(k, 956)
         mat(k, 958) = lmat(k, 958)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 962) = mat(k, 962) + lmat(k, 962)
         mat(k, 963) = lmat(k, 963)
         mat(k, 965) = mat(k, 965) + lmat(k, 965)
         mat(k, 966) = mat(k, 966) + lmat(k, 966)
         mat(k, 967) = mat(k, 967) + lmat(k, 967)
         mat(k, 968) = mat(k, 968) + lmat(k, 968)
         mat(k, 969) = mat(k, 969) + lmat(k, 969)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 971) = mat(k, 971) + lmat(k, 971)
         mat(k, 974) = mat(k, 974) + lmat(k, 974)
         mat(k, 975) = mat(k, 975) + lmat(k, 975)
         mat(k, 977) = mat(k, 977) + lmat(k, 977)
         mat(k, 978) = lmat(k, 978)
         mat(k, 980) = lmat(k, 980)
         mat(k, 984) = lmat(k, 984)
         mat(k, 987) = lmat(k, 987)
         mat(k, 988) = lmat(k, 988)
         mat(k, 989) = lmat(k, 989)
         mat(k, 990) = mat(k, 990) + lmat(k, 990)
         mat(k, 991) = lmat(k, 991)
         mat(k, 992) = lmat(k, 992)
         mat(k, 994) = lmat(k, 994)
         mat(k, 998) = lmat(k, 998)
         mat(k, 999) = lmat(k, 999)
         mat(k,1000) = mat(k,1000) + lmat(k,1000)
         mat(k,1001) = lmat(k,1001)
         mat(k,1003) = mat(k,1003) + lmat(k,1003)
         mat(k,1007) = mat(k,1007) + lmat(k,1007)
         mat(k,1009) = lmat(k,1009)
         mat(k,1011) = mat(k,1011) + lmat(k,1011)
         mat(k,1012) = lmat(k,1012)
         mat(k,1020) = mat(k,1020) + lmat(k,1020)
         mat(k,1042) = mat(k,1042) + lmat(k,1042)
         mat(k,1061) = mat(k,1061) + lmat(k,1061)
         mat(k,1077) = mat(k,1077) + lmat(k,1077)
         mat(k,1089) = mat(k,1089) + lmat(k,1089)
         mat(k,1106) = mat(k,1106) + lmat(k,1106)
         mat(k,1126) = mat(k,1126) + lmat(k,1126)
         mat(k,1141) = mat(k,1141) + lmat(k,1141)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1150) = mat(k,1150) + lmat(k,1150)
         mat(k,1151) = mat(k,1151) + lmat(k,1151)
         mat(k,1153) = mat(k,1153) + lmat(k,1153)
         mat(k,1154) = mat(k,1154) + lmat(k,1154)
         mat(k,1155) = mat(k,1155) + lmat(k,1155)
         mat(k,1160) = lmat(k,1160)
         mat(k,1162) = lmat(k,1162)
         mat(k,1179) = mat(k,1179) + lmat(k,1179)
         mat(k,1192) = mat(k,1192) + lmat(k,1192)
         mat(k,1206) = mat(k,1206) + lmat(k,1206)
         mat(k,1230) = mat(k,1230) + lmat(k,1230)
         mat(k,1244) = lmat(k,1244)
         mat(k,1246) = mat(k,1246) + lmat(k,1246)
         mat(k,1250) = mat(k,1250) + lmat(k,1250)
         mat(k,1252) = mat(k,1252) + lmat(k,1252)
         mat(k,1253) = lmat(k,1253)
         mat(k,1270) = mat(k,1270) + lmat(k,1270)
         mat(k,1301) = mat(k,1301) + lmat(k,1301)
         mat(k,1315) = lmat(k,1315)
         mat(k,1317) = mat(k,1317) + lmat(k,1317)
         mat(k,1325) = mat(k,1325) + lmat(k,1325)
         mat(k,1336) = mat(k,1336) + lmat(k,1336)
         mat(k,1342) = mat(k,1342) + lmat(k,1342)
         mat(k,1345) = mat(k,1345) + lmat(k,1345)
         mat(k,1351) = mat(k,1351) + lmat(k,1351)
         mat(k,1395) = mat(k,1395) + lmat(k,1395)
         mat(k,1414) = mat(k,1414) + lmat(k,1414)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1419) = lmat(k,1419)
         mat(k,1426) = mat(k,1426) + lmat(k,1426)
         mat(k,1431) = mat(k,1431) + lmat(k,1431)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1481) = mat(k,1481) + lmat(k,1481)
         mat(k,1482) = mat(k,1482) + lmat(k,1482)
         mat(k,1486) = mat(k,1486) + lmat(k,1486)
         mat(k,1489) = mat(k,1489) + lmat(k,1489)
         mat(k,1490) = mat(k,1490) + lmat(k,1490)
         mat(k,1492) = mat(k,1492) + lmat(k,1492)
         mat(k,1510) = mat(k,1510) + lmat(k,1510)
         mat(k,1517) = mat(k,1517) + lmat(k,1517)
         mat(k,1518) = mat(k,1518) + lmat(k,1518)
         mat(k,1529) = mat(k,1529) + lmat(k,1529)
         mat(k,1570) = mat(k,1570) + lmat(k,1570)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1577) = mat(k,1577) + lmat(k,1577)
         mat(k,1579) = mat(k,1579) + lmat(k,1579)
         mat(k,1596) = lmat(k,1596)
         mat(k,1606) = lmat(k,1606)
         mat(k,1716) = mat(k,1716) + lmat(k,1716)
         mat(k,1717) = mat(k,1717) + lmat(k,1717)
         mat(k,1722) = mat(k,1722) + lmat(k,1722)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1730) = mat(k,1730) + lmat(k,1730)
         mat(k,1733) = mat(k,1733) + lmat(k,1733)
         mat(k,1744) = mat(k,1744) + lmat(k,1744)
         mat(k,1745) = mat(k,1745) + lmat(k,1745)
         mat(k,1751) = lmat(k,1751)
         mat(k,1766) = mat(k,1766) + lmat(k,1766)
         mat(k,1776) = mat(k,1776) + lmat(k,1776)
         mat(k,1779) = mat(k,1779) + lmat(k,1779)
         mat(k,1780) = mat(k,1780) + lmat(k,1780)
         mat(k,1782) = mat(k,1782) + lmat(k,1782)
         mat(k,1784) = mat(k,1784) + lmat(k,1784)
         mat(k,1785) = lmat(k,1785)
         mat(k,1786) = mat(k,1786) + lmat(k,1786)
         mat(k,1790) = mat(k,1790) + lmat(k,1790)
         mat(k,1792) = mat(k,1792) + lmat(k,1792)
         mat(k,1793) = mat(k,1793) + lmat(k,1793)
         mat(k,1794) = lmat(k,1794)
         mat(k,1795) = mat(k,1795) + lmat(k,1795)
         mat(k,1796) = mat(k,1796) + lmat(k,1796)
         mat(k,1798) = lmat(k,1798)
         mat(k,1799) = lmat(k,1799)
         mat(k,1831) = mat(k,1831) + lmat(k,1831)
         mat(k,1832) = lmat(k,1832)
         mat(k,1835) = mat(k,1835) + lmat(k,1835)
         mat(k,1875) = mat(k,1875) + lmat(k,1875)
         mat(k,1876) = mat(k,1876) + lmat(k,1876)
         mat(k,1890) = mat(k,1890) + lmat(k,1890)
         mat(k,1895) = lmat(k,1895)
         mat(k,1917) = mat(k,1917) + lmat(k,1917)
         mat(k,1937) = mat(k,1937) + lmat(k,1937)
         mat(k,1940) = mat(k,1940) + lmat(k,1940)
         mat(k,1941) = lmat(k,1941)
         mat(k,1942) = lmat(k,1942)
         mat(k,1953) = mat(k,1953) + lmat(k,1953)
         mat(k,1955) = mat(k,1955) + lmat(k,1955)
         mat(k,1989) = mat(k,1989) + lmat(k,1989)
         mat(k,1990) = mat(k,1990) + lmat(k,1990)
         mat(k,1993) = mat(k,1993) + lmat(k,1993)
         mat(k,1994) = mat(k,1994) + lmat(k,1994)
         mat(k,1996) = mat(k,1996) + lmat(k,1996)
         mat(k,2032) = mat(k,2032) + lmat(k,2032)
         mat(k,2092) = mat(k,2092) + lmat(k,2092)
         mat(k,2099) = mat(k,2099) + lmat(k,2099)
         mat(k,2100) = lmat(k,2100)
         mat(k,2110) = mat(k,2110) + lmat(k,2110)
         mat(k,2117) = mat(k,2117) + lmat(k,2117)
         mat(k,2127) = mat(k,2127) + lmat(k,2127)
         mat(k,2137) = mat(k,2137) + lmat(k,2137)
         mat(k,2142) = mat(k,2142) + lmat(k,2142)
         mat(k,2149) = lmat(k,2149)
         mat(k,2158) = mat(k,2158) + lmat(k,2158)
         mat(k,2160) = lmat(k,2160)
         mat(k,2161) = mat(k,2161) + lmat(k,2161)
         mat(k,2163) = lmat(k,2163)
         mat(k,2169) = mat(k,2169) + lmat(k,2169)
         mat(k, 145) = 0._r8
         mat(k, 146) = 0._r8
         mat(k, 242) = 0._r8
         mat(k, 326) = 0._r8
         mat(k, 327) = 0._r8
         mat(k, 340) = 0._r8
         mat(k, 380) = 0._r8
         mat(k, 382) = 0._r8
         mat(k, 397) = 0._r8
         mat(k, 498) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 532) = 0._r8
         mat(k, 533) = 0._r8
         mat(k, 539) = 0._r8
         mat(k, 540) = 0._r8
         mat(k, 543) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 557) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 567) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 601) = 0._r8
         mat(k, 603) = 0._r8
         mat(k, 604) = 0._r8
         mat(k, 606) = 0._r8
         mat(k, 619) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 623) = 0._r8
         mat(k, 624) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 686) = 0._r8
         mat(k, 687) = 0._r8
         mat(k, 688) = 0._r8
         mat(k, 696) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 707) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 723) = 0._r8
         mat(k, 728) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 735) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 806) = 0._r8
         mat(k, 833) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 843) = 0._r8
         mat(k, 850) = 0._r8
         mat(k, 859) = 0._r8
         mat(k, 861) = 0._r8
         mat(k, 869) = 0._r8
         mat(k, 876) = 0._r8
         mat(k, 880) = 0._r8
         mat(k, 881) = 0._r8
         mat(k, 885) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 937) = 0._r8
         mat(k, 941) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 949) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 953) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 976) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 995) = 0._r8
         mat(k, 996) = 0._r8
         mat(k,1002) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1019) = 0._r8
         mat(k,1021) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1023) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1034) = 0._r8
         mat(k,1043) = 0._r8
         mat(k,1044) = 0._r8
         mat(k,1046) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1059) = 0._r8
         mat(k,1060) = 0._r8
         mat(k,1062) = 0._r8
         mat(k,1063) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1092) = 0._r8
         mat(k,1094) = 0._r8
         mat(k,1111) = 0._r8
         mat(k,1112) = 0._r8
         mat(k,1119) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1124) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1129) = 0._r8
         mat(k,1133) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1159) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1172) = 0._r8
         mat(k,1173) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1180) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1202) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1214) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1216) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1237) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1269) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1282) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1328) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1356) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1370) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1401) = 0._r8
         mat(k,1402) = 0._r8
         mat(k,1408) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1429) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1439) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1459) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1483) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1495) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1511) = 0._r8
         mat(k,1513) = 0._r8
         mat(k,1514) = 0._r8
         mat(k,1515) = 0._r8
         mat(k,1523) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1541) = 0._r8
         mat(k,1544) = 0._r8
         mat(k,1549) = 0._r8
         mat(k,1553) = 0._r8
         mat(k,1554) = 0._r8
         mat(k,1555) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1558) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1563) = 0._r8
         mat(k,1565) = 0._r8
         mat(k,1575) = 0._r8
         mat(k,1585) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1657) = 0._r8
         mat(k,1666) = 0._r8
         mat(k,1669) = 0._r8
         mat(k,1670) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1705) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1739) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1747) = 0._r8
         mat(k,1748) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1752) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1759) = 0._r8
         mat(k,1761) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1773) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1791) = 0._r8
         mat(k,1797) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1872) = 0._r8
         mat(k,1873) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1882) = 0._r8
         mat(k,1892) = 0._r8
         mat(k,1893) = 0._r8
         mat(k,1900) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1935) = 0._r8
         mat(k,1936) = 0._r8
         mat(k,1939) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1967) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1977) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1998) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2019) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2050) = 0._r8
         mat(k,2055) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2057) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2069) = 0._r8
         mat(k,2074) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2103) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2107) = 0._r8
         mat(k,2111) = 0._r8
         mat(k,2112) = 0._r8
         mat(k,2115) = 0._r8
         mat(k,2118) = 0._r8
         mat(k,2126) = 0._r8
         mat(k,2129) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2153) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2159) = 0._r8
         mat(k,2162) = 0._r8
         mat(k,2164) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2166) = 0._r8
         mat(k,2167) = 0._r8
         mat(k,2168) = 0._r8
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
         mat(k, 41) = mat(k, 41) - dti(k)
         mat(k, 47) = mat(k, 47) - dti(k)
         mat(k, 48) = mat(k, 48) - dti(k)
         mat(k, 51) = mat(k, 51) - dti(k)
         mat(k, 54) = mat(k, 54) - dti(k)
         mat(k, 57) = mat(k, 57) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 64) = mat(k, 64) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 70) = mat(k, 70) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 76) = mat(k, 76) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 96) = mat(k, 96) - dti(k)
         mat(k, 101) = mat(k, 101) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 112) = mat(k, 112) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 124) = mat(k, 124) - dti(k)
         mat(k, 128) = mat(k, 128) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 154) = mat(k, 154) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 165) = mat(k, 165) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 171) = mat(k, 171) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 187) = mat(k, 187) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 199) = mat(k, 199) - dti(k)
         mat(k, 202) = mat(k, 202) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 219) = mat(k, 219) - dti(k)
         mat(k, 224) = mat(k, 224) - dti(k)
         mat(k, 229) = mat(k, 229) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 245) = mat(k, 245) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 261) = mat(k, 261) - dti(k)
         mat(k, 269) = mat(k, 269) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 281) = mat(k, 281) - dti(k)
         mat(k, 287) = mat(k, 287) - dti(k)
         mat(k, 293) = mat(k, 293) - dti(k)
         mat(k, 299) = mat(k, 299) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 325) = mat(k, 325) - dti(k)
         mat(k, 331) = mat(k, 331) - dti(k)
         mat(k, 338) = mat(k, 338) - dti(k)
         mat(k, 344) = mat(k, 344) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 352) = mat(k, 352) - dti(k)
         mat(k, 359) = mat(k, 359) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 370) = mat(k, 370) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 394) = mat(k, 394) - dti(k)
         mat(k, 401) = mat(k, 401) - dti(k)
         mat(k, 406) = mat(k, 406) - dti(k)
         mat(k, 413) = mat(k, 413) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 424) = mat(k, 424) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 440) = mat(k, 440) - dti(k)
         mat(k, 448) = mat(k, 448) - dti(k)
         mat(k, 452) = mat(k, 452) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 468) = mat(k, 468) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 485) = mat(k, 485) - dti(k)
         mat(k, 496) = mat(k, 496) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 512) = mat(k, 512) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 523) = mat(k, 523) - dti(k)
         mat(k, 531) = mat(k, 531) - dti(k)
         mat(k, 538) = mat(k, 538) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 561) = mat(k, 561) - dti(k)
         mat(k, 572) = mat(k, 572) - dti(k)
         mat(k, 580) = mat(k, 580) - dti(k)
         mat(k, 587) = mat(k, 587) - dti(k)
         mat(k, 598) = mat(k, 598) - dti(k)
         mat(k, 609) = mat(k, 609) - dti(k)
         mat(k, 620) = mat(k, 620) - dti(k)
         mat(k, 636) = mat(k, 636) - dti(k)
         mat(k, 647) = mat(k, 647) - dti(k)
         mat(k, 656) = mat(k, 656) - dti(k)
         mat(k, 666) = mat(k, 666) - dti(k)
         mat(k, 674) = mat(k, 674) - dti(k)
         mat(k, 683) = mat(k, 683) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 702) = mat(k, 702) - dti(k)
         mat(k, 710) = mat(k, 710) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 719) = mat(k, 719) - dti(k)
         mat(k, 734) = mat(k, 734) - dti(k)
         mat(k, 744) = mat(k, 744) - dti(k)
         mat(k, 756) = mat(k, 756) - dti(k)
         mat(k, 763) = mat(k, 763) - dti(k)
         mat(k, 777) = mat(k, 777) - dti(k)
         mat(k, 793) = mat(k, 793) - dti(k)
         mat(k, 802) = mat(k, 802) - dti(k)
         mat(k, 812) = mat(k, 812) - dti(k)
         mat(k, 817) = mat(k, 817) - dti(k)
         mat(k, 834) = mat(k, 834) - dti(k)
         mat(k, 860) = mat(k, 860) - dti(k)
         mat(k, 882) = mat(k, 882) - dti(k)
         mat(k, 894) = mat(k, 894) - dti(k)
         mat(k, 900) = mat(k, 900) - dti(k)
         mat(k, 908) = mat(k, 908) - dti(k)
         mat(k, 918) = mat(k, 918) - dti(k)
         mat(k, 930) = mat(k, 930) - dti(k)
         mat(k, 950) = mat(k, 950) - dti(k)
         mat(k, 962) = mat(k, 962) - dti(k)
         mat(k, 968) = mat(k, 968) - dti(k)
         mat(k, 977) = mat(k, 977) - dti(k)
         mat(k, 990) = mat(k, 990) - dti(k)
         mat(k,1003) = mat(k,1003) - dti(k)
         mat(k,1007) = mat(k,1007) - dti(k)
         mat(k,1020) = mat(k,1020) - dti(k)
         mat(k,1042) = mat(k,1042) - dti(k)
         mat(k,1061) = mat(k,1061) - dti(k)
         mat(k,1077) = mat(k,1077) - dti(k)
         mat(k,1089) = mat(k,1089) - dti(k)
         mat(k,1106) = mat(k,1106) - dti(k)
         mat(k,1126) = mat(k,1126) - dti(k)
         mat(k,1142) = mat(k,1142) - dti(k)
         mat(k,1154) = mat(k,1154) - dti(k)
         mat(k,1179) = mat(k,1179) - dti(k)
         mat(k,1206) = mat(k,1206) - dti(k)
         mat(k,1230) = mat(k,1230) - dti(k)
         mat(k,1250) = mat(k,1250) - dti(k)
         mat(k,1270) = mat(k,1270) - dti(k)
         mat(k,1301) = mat(k,1301) - dti(k)
         mat(k,1317) = mat(k,1317) - dti(k)
         mat(k,1336) = mat(k,1336) - dti(k)
         mat(k,1351) = mat(k,1351) - dti(k)
         mat(k,1395) = mat(k,1395) - dti(k)
         mat(k,1426) = mat(k,1426) - dti(k)
         mat(k,1482) = mat(k,1482) - dti(k)
         mat(k,1510) = mat(k,1510) - dti(k)
         mat(k,1573) = mat(k,1573) - dti(k)
         mat(k,1722) = mat(k,1722) - dti(k)
         mat(k,1745) = mat(k,1745) - dti(k)
         mat(k,1766) = mat(k,1766) - dti(k)
         mat(k,1793) = mat(k,1793) - dti(k)
         mat(k,1875) = mat(k,1875) - dti(k)
         mat(k,1917) = mat(k,1917) - dti(k)
         mat(k,1953) = mat(k,1953) - dti(k)
         mat(k,1996) = mat(k,1996) - dti(k)
         mat(k,2092) = mat(k,2092) - dti(k)
         mat(k,2117) = mat(k,2117) - dti(k)
         mat(k,2142) = mat(k,2142) - dti(k)
         mat(k,2169) = mat(k,2169) - dti(k)
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
