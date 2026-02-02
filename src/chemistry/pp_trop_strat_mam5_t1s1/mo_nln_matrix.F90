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
         mat(k,667) = -(rxt(k,357)*y(k,218))
         mat(k,1632) = -rxt(k,357)*y(k,1)
         mat(k,1874) = rxt(k,360)*y(k,190)
         mat(k,885) = rxt(k,360)*y(k,123)
         mat(k,633) = -(rxt(k,361)*y(k,218))
         mat(k,1629) = -rxt(k,361)*y(k,2)
         mat(k,884) = rxt(k,358)*y(k,204)
         mat(k,2057) = rxt(k,358)*y(k,190)
         mat(k,1000) = -(rxt(k,440)*y(k,125) + rxt(k,441)*y(k,134) + rxt(k,442) &
                      *y(k,218))
         mat(k,1746) = -rxt(k,440)*y(k,5)
         mat(k,2222) = -rxt(k,441)*y(k,5)
         mat(k,1662) = -rxt(k,442)*y(k,5)
         mat(k,163) = -(rxt(k,399)*y(k,218))
         mat(k,1560) = -rxt(k,399)*y(k,6)
         mat(k,416) = -(rxt(k,402)*y(k,218))
         mat(k,1600) = -rxt(k,402)*y(k,7)
         mat(k,482) = rxt(k,400)*y(k,204)
         mat(k,2040) = rxt(k,400)*y(k,192)
         mat(k,164) = .120_r8*rxt(k,399)*y(k,218)
         mat(k,1561) = .120_r8*rxt(k,399)*y(k,6)
         mat(k,992) = .100_r8*rxt(k,441)*y(k,134)
         mat(k,943) = .100_r8*rxt(k,444)*y(k,134)
         mat(k,2206) = .100_r8*rxt(k,441)*y(k,5) + .100_r8*rxt(k,444)*y(k,109)
         mat(k,1861) = .500_r8*rxt(k,401)*y(k,192) + .200_r8*rxt(k,428)*y(k,224) &
                      + .060_r8*rxt(k,434)*y(k,227)
         mat(k,483) = .500_r8*rxt(k,401)*y(k,123)
         mat(k,721) = .200_r8*rxt(k,428)*y(k,123)
         mat(k,745) = .060_r8*rxt(k,434)*y(k,123)
         mat(k,1855) = .200_r8*rxt(k,428)*y(k,224) + .200_r8*rxt(k,434)*y(k,227)
         mat(k,720) = .200_r8*rxt(k,428)*y(k,123)
         mat(k,743) = .200_r8*rxt(k,434)*y(k,123)
         mat(k,1871) = .200_r8*rxt(k,428)*y(k,224) + .150_r8*rxt(k,434)*y(k,227)
         mat(k,722) = .200_r8*rxt(k,428)*y(k,123)
         mat(k,746) = .150_r8*rxt(k,434)*y(k,123)
         mat(k,1856) = .210_r8*rxt(k,434)*y(k,227)
         mat(k,744) = .210_r8*rxt(k,434)*y(k,123)
         mat(k,227) = -(rxt(k,362)*y(k,218))
         mat(k,1571) = -rxt(k,362)*y(k,14)
         mat(k,991) = .050_r8*rxt(k,441)*y(k,134)
         mat(k,942) = .050_r8*rxt(k,444)*y(k,134)
         mat(k,2205) = .050_r8*rxt(k,441)*y(k,5) + .050_r8*rxt(k,444)*y(k,109)
         mat(k,350) = -(rxt(k,328)*y(k,125) + rxt(k,329)*y(k,218))
         mat(k,1736) = -rxt(k,328)*y(k,15)
         mat(k,1591) = -rxt(k,329)*y(k,15)
         mat(k,1416) = -(rxt(k,212)*y(k,41) + rxt(k,213)*y(k,204) + rxt(k,214) &
                      *y(k,134))
         mat(k,1482) = -rxt(k,212)*y(k,16)
         mat(k,2103) = -rxt(k,213)*y(k,16)
         mat(k,2242) = -rxt(k,214)*y(k,16)
         mat(k,2151) = 4.000_r8*rxt(k,215)*y(k,18) + (rxt(k,216)+rxt(k,217))*y(k,58) &
                      + rxt(k,220)*y(k,123) + rxt(k,223)*y(k,133) + rxt(k,470) &
                      *y(k,150) + rxt(k,224)*y(k,218)
         mat(k,144) = rxt(k,202)*y(k,217)
         mat(k,150) = rxt(k,228)*y(k,217)
         mat(k,469) = 2.000_r8*rxt(k,239)*y(k,55) + 2.000_r8*rxt(k,251)*y(k,217) &
                      + 2.000_r8*rxt(k,240)*y(k,218)
         mat(k,596) = rxt(k,241)*y(k,55) + rxt(k,252)*y(k,217) + rxt(k,242)*y(k,218)
         mat(k,387) = 3.000_r8*rxt(k,246)*y(k,55) + 3.000_r8*rxt(k,229)*y(k,217) &
                      + 3.000_r8*rxt(k,247)*y(k,218)
         mat(k,1996) = 2.000_r8*rxt(k,239)*y(k,40) + rxt(k,241)*y(k,42) &
                      + 3.000_r8*rxt(k,246)*y(k,54)
         mat(k,1713) = (rxt(k,216)+rxt(k,217))*y(k,18)
         mat(k,108) = 2.000_r8*rxt(k,230)*y(k,217)
         mat(k,806) = rxt(k,225)*y(k,133) + rxt(k,231)*y(k,217) + rxt(k,226)*y(k,218)
         mat(k,1913) = rxt(k,220)*y(k,18)
         mat(k,2181) = rxt(k,223)*y(k,18) + rxt(k,225)*y(k,80)
         mat(k,1234) = rxt(k,470)*y(k,18)
         mat(k,1522) = rxt(k,202)*y(k,33) + rxt(k,228)*y(k,34) + 2.000_r8*rxt(k,251) &
                      *y(k,40) + rxt(k,252)*y(k,42) + 3.000_r8*rxt(k,229)*y(k,54) &
                      + 2.000_r8*rxt(k,230)*y(k,77) + rxt(k,231)*y(k,80)
         mat(k,1686) = rxt(k,224)*y(k,18) + 2.000_r8*rxt(k,240)*y(k,40) + rxt(k,242) &
                      *y(k,42) + 3.000_r8*rxt(k,247)*y(k,54) + rxt(k,226)*y(k,80)
         mat(k,2145) = rxt(k,218)*y(k,58)
         mat(k,1707) = rxt(k,218)*y(k,18)
         mat(k,2123) = (rxt(k,531)+rxt(k,536))*y(k,90)
         mat(k,778) = (rxt(k,531)+rxt(k,536))*y(k,84)
         mat(k,2165) = -(4._r8*rxt(k,215)*y(k,18) + (rxt(k,216) + rxt(k,217) + rxt(k,218) &
                      ) * y(k,58) + rxt(k,219)*y(k,204) + rxt(k,220)*y(k,123) &
                      + rxt(k,221)*y(k,124) + rxt(k,223)*y(k,133) + rxt(k,224) &
                      *y(k,218) + rxt(k,470)*y(k,150))
         mat(k,1727) = -(rxt(k,216) + rxt(k,217) + rxt(k,218)) * y(k,18)
         mat(k,2118) = -rxt(k,219)*y(k,18)
         mat(k,1928) = -rxt(k,220)*y(k,18)
         mat(k,1972) = -rxt(k,221)*y(k,18)
         mat(k,2196) = -rxt(k,223)*y(k,18)
         mat(k,1701) = -rxt(k,224)*y(k,18)
         mat(k,1242) = -rxt(k,470)*y(k,18)
         mat(k,1422) = rxt(k,214)*y(k,134)
         mat(k,563) = rxt(k,222)*y(k,133)
         mat(k,810) = rxt(k,232)*y(k,217)
         mat(k,784) = rxt(k,227)*y(k,133)
         mat(k,2196) = mat(k,2196) + rxt(k,222)*y(k,19) + rxt(k,227)*y(k,90)
         mat(k,2257) = rxt(k,214)*y(k,16)
         mat(k,1537) = rxt(k,232)*y(k,80)
         mat(k,557) = -(rxt(k,222)*y(k,133))
         mat(k,2171) = -rxt(k,222)*y(k,19)
         mat(k,2147) = rxt(k,221)*y(k,124)
         mat(k,1940) = rxt(k,221)*y(k,18)
         mat(k,236) = -(rxt(k,403)*y(k,218))
         mat(k,1573) = -rxt(k,403)*y(k,21)
         mat(k,1853) = rxt(k,406)*y(k,194)
         mat(k,434) = rxt(k,406)*y(k,123)
         mat(k,337) = -(rxt(k,405)*y(k,218))
         mat(k,1588) = -rxt(k,405)*y(k,22)
         mat(k,435) = rxt(k,404)*y(k,204)
         mat(k,2034) = rxt(k,404)*y(k,194)
         mat(k,292) = -(rxt(k,277)*y(k,55) + rxt(k,278)*y(k,218))
         mat(k,1977) = -rxt(k,277)*y(k,23)
         mat(k,1582) = -rxt(k,278)*y(k,23)
         mat(k,549) = -(rxt(k,279)*y(k,55) + rxt(k,280)*y(k,134) + rxt(k,305)*y(k,218))
         mat(k,1982) = -rxt(k,279)*y(k,24)
         mat(k,2209) = -rxt(k,280)*y(k,24)
         mat(k,1619) = -rxt(k,305)*y(k,24)
         mat(k,262) = -(rxt(k,285)*y(k,218))
         mat(k,1579) = -rxt(k,285)*y(k,25)
         mat(k,821) = .800_r8*rxt(k,281)*y(k,195) + .200_r8*rxt(k,282)*y(k,199)
         mat(k,1788) = .200_r8*rxt(k,282)*y(k,195)
         mat(k,342) = -(rxt(k,286)*y(k,218))
         mat(k,1589) = -rxt(k,286)*y(k,26)
         mat(k,822) = rxt(k,283)*y(k,204)
         mat(k,2035) = rxt(k,283)*y(k,195)
         mat(k,305) = -(rxt(k,287)*y(k,55) + rxt(k,288)*y(k,218))
         mat(k,1978) = -rxt(k,287)*y(k,27)
         mat(k,1584) = -rxt(k,288)*y(k,27)
         mat(k,1025) = -(rxt(k,308)*y(k,125) + rxt(k,309)*y(k,134) + rxt(k,326) &
                      *y(k,218))
         mat(k,1747) = -rxt(k,308)*y(k,28)
         mat(k,2223) = -rxt(k,309)*y(k,28)
         mat(k,1663) = -rxt(k,326)*y(k,28)
         mat(k,837) = .130_r8*rxt(k,386)*y(k,134)
         mat(k,2223) = mat(k,2223) + .130_r8*rxt(k,386)*y(k,97)
         mat(k,410) = -(rxt(k,313)*y(k,218))
         mat(k,1599) = -rxt(k,313)*y(k,29)
         mat(k,787) = rxt(k,311)*y(k,204)
         mat(k,2039) = rxt(k,311)*y(k,196)
         mat(k,110) = -(rxt(k,314)*y(k,218))
         mat(k,1557) = -rxt(k,314)*y(k,30)
         mat(k,266) = -(rxt(k,409)*y(k,218))
         mat(k,1580) = -rxt(k,409)*y(k,31)
         mat(k,624) = rxt(k,407)*y(k,204)
         mat(k,2029) = rxt(k,407)*y(k,197)
         mat(k,100) = -(rxt(k,201)*y(k,217))
         mat(k,1500) = -rxt(k,201)*y(k,32)
         mat(k,142) = -(rxt(k,202)*y(k,217))
         mat(k,1505) = -rxt(k,202)*y(k,33)
         mat(k,147) = -(rxt(k,228)*y(k,217))
         mat(k,1506) = -rxt(k,228)*y(k,34)
         mat(k,114) = -(rxt(k,203)*y(k,217))
         mat(k,1502) = -rxt(k,203)*y(k,35)
         mat(k,152) = -(rxt(k,204)*y(k,217))
         mat(k,1507) = -rxt(k,204)*y(k,36)
         mat(k,118) = -(rxt(k,205)*y(k,217))
         mat(k,1503) = -rxt(k,205)*y(k,37)
         mat(k,157) = -(rxt(k,206)*y(k,217))
         mat(k,1508) = -rxt(k,206)*y(k,38)
         mat(k,122) = -(rxt(k,207)*y(k,217))
         mat(k,1504) = -rxt(k,207)*y(k,39)
         mat(k,468) = -(rxt(k,239)*y(k,55) + rxt(k,240)*y(k,218) + rxt(k,251)*y(k,217))
         mat(k,1981) = -rxt(k,239)*y(k,40)
         mat(k,1608) = -rxt(k,240)*y(k,40)
         mat(k,1517) = -rxt(k,251)*y(k,40)
         mat(k,1486) = -(rxt(k,176)*y(k,55) + rxt(k,212)*y(k,16) + rxt(k,256)*y(k,204) &
                      + rxt(k,257)*y(k,125) + rxt(k,258)*y(k,133) + rxt(k,259) &
                      *y(k,218))
         mat(k,2000) = -rxt(k,176)*y(k,41)
         mat(k,1418) = -rxt(k,212)*y(k,41)
         mat(k,2107) = -rxt(k,256)*y(k,41)
         mat(k,1773) = -rxt(k,257)*y(k,41)
         mat(k,2185) = -rxt(k,258)*y(k,41)
         mat(k,1690) = -rxt(k,259)*y(k,41)
         mat(k,673) = .400_r8*rxt(k,357)*y(k,218)
         mat(k,1010) = .340_r8*rxt(k,441)*y(k,134)
         mat(k,354) = .500_r8*rxt(k,328)*y(k,125)
         mat(k,553) = rxt(k,280)*y(k,134)
         mat(k,1032) = .500_r8*rxt(k,309)*y(k,134)
         mat(k,614) = .500_r8*rxt(k,297)*y(k,218)
         mat(k,799) = rxt(k,264)*y(k,218)
         mat(k,457) = .300_r8*rxt(k,265)*y(k,218)
         mat(k,1434) = (rxt(k,273)+rxt(k,274))*y(k,217)
         mat(k,1716) = rxt(k,183)*y(k,199)
         mat(k,1099) = .800_r8*rxt(k,302)*y(k,218)
         mat(k,845) = .910_r8*rxt(k,386)*y(k,134)
         mat(k,587) = .300_r8*rxt(k,377)*y(k,218)
         mat(k,1200) = .800_r8*rxt(k,381)*y(k,199)
         mat(k,1215) = .120_r8*rxt(k,339)*y(k,134)
         mat(k,577) = .500_r8*rxt(k,352)*y(k,218)
         mat(k,960) = .340_r8*rxt(k,444)*y(k,134)
         mat(k,1338) = .600_r8*rxt(k,353)*y(k,134)
         mat(k,1917) = .100_r8*rxt(k,359)*y(k,190) + rxt(k,263)*y(k,199) &
                      + .500_r8*rxt(k,330)*y(k,201) + .500_r8*rxt(k,299)*y(k,203) &
                      + .920_r8*rxt(k,369)*y(k,206) + .250_r8*rxt(k,337)*y(k,210) &
                      + rxt(k,346)*y(k,212) + rxt(k,320)*y(k,220) + rxt(k,324) &
                      *y(k,221) + .340_r8*rxt(k,453)*y(k,222) + .320_r8*rxt(k,458) &
                      *y(k,223) + .250_r8*rxt(k,394)*y(k,226)
         mat(k,1773) = mat(k,1773) + .500_r8*rxt(k,328)*y(k,15) + rxt(k,370)*y(k,206) &
                      + .250_r8*rxt(k,336)*y(k,210) + rxt(k,347)*y(k,212)
         mat(k,2246) = .340_r8*rxt(k,441)*y(k,5) + rxt(k,280)*y(k,24) &
                      + .500_r8*rxt(k,309)*y(k,28) + .910_r8*rxt(k,386)*y(k,97) &
                      + .120_r8*rxt(k,339)*y(k,104) + .340_r8*rxt(k,444)*y(k,109) &
                      + .600_r8*rxt(k,353)*y(k,110)
         mat(k,528) = rxt(k,304)*y(k,218)
         mat(k,1064) = .680_r8*rxt(k,462)*y(k,218)
         mat(k,892) = .100_r8*rxt(k,359)*y(k,123)
         mat(k,826) = .700_r8*rxt(k,282)*y(k,199)
         mat(k,791) = rxt(k,310)*y(k,199)
         mat(k,1390) = rxt(k,293)*y(k,199) + rxt(k,366)*y(k,206) + .250_r8*rxt(k,333) &
                      *y(k,210) + rxt(k,342)*y(k,212) + .250_r8*rxt(k,391)*y(k,226)
         mat(k,1825) = rxt(k,183)*y(k,58) + .800_r8*rxt(k,381)*y(k,100) + rxt(k,263) &
                      *y(k,123) + .700_r8*rxt(k,282)*y(k,195) + rxt(k,310)*y(k,196) &
                      + rxt(k,293)*y(k,198) + (4.000_r8*rxt(k,260)+2.000_r8*rxt(k,261)) &
                      *y(k,199) + 1.500_r8*rxt(k,367)*y(k,206) + .750_r8*rxt(k,372) &
                      *y(k,207) + .880_r8*rxt(k,334)*y(k,210) + 2.000_r8*rxt(k,343) &
                      *y(k,212) + .750_r8*rxt(k,446)*y(k,216) + .800_r8*rxt(k,322) &
                      *y(k,221) + .930_r8*rxt(k,451)*y(k,222) + .950_r8*rxt(k,456) &
                      *y(k,223) + .800_r8*rxt(k,392)*y(k,226)
         mat(k,569) = .500_r8*rxt(k,330)*y(k,123)
         mat(k,709) = .500_r8*rxt(k,299)*y(k,123)
         mat(k,2107) = mat(k,2107) + .450_r8*rxt(k,344)*y(k,212) + .150_r8*rxt(k,323) &
                      *y(k,221)
         mat(k,1263) = .920_r8*rxt(k,369)*y(k,123) + rxt(k,370)*y(k,125) + rxt(k,366) &
                      *y(k,198) + 1.500_r8*rxt(k,367)*y(k,199)
         mat(k,1295) = .750_r8*rxt(k,372)*y(k,199)
         mat(k,1316) = .250_r8*rxt(k,337)*y(k,123) + .250_r8*rxt(k,336)*y(k,125) &
                      + .250_r8*rxt(k,333)*y(k,198) + .880_r8*rxt(k,334)*y(k,199)
         mat(k,1358) = rxt(k,346)*y(k,123) + rxt(k,347)*y(k,125) + rxt(k,342)*y(k,198) &
                      + 2.000_r8*rxt(k,343)*y(k,199) + .450_r8*rxt(k,344)*y(k,204) &
                      + 4.000_r8*rxt(k,345)*y(k,212)
         mat(k,1051) = .750_r8*rxt(k,446)*y(k,199)
         mat(k,1526) = (rxt(k,273)+rxt(k,274))*y(k,53)
         mat(k,1690) = mat(k,1690) + .400_r8*rxt(k,357)*y(k,1) + .500_r8*rxt(k,297) &
                      *y(k,50) + rxt(k,264)*y(k,51) + .300_r8*rxt(k,265)*y(k,52) &
                      + .800_r8*rxt(k,302)*y(k,73) + .300_r8*rxt(k,377)*y(k,98) &
                      + .500_r8*rxt(k,352)*y(k,108) + rxt(k,304)*y(k,139) &
                      + .680_r8*rxt(k,462)*y(k,179)
         mat(k,772) = rxt(k,320)*y(k,123)
         mat(k,1134) = rxt(k,324)*y(k,123) + .800_r8*rxt(k,322)*y(k,199) &
                      + .150_r8*rxt(k,323)*y(k,204)
         mat(k,1115) = .340_r8*rxt(k,453)*y(k,123) + .930_r8*rxt(k,451)*y(k,199)
         mat(k,917) = .320_r8*rxt(k,458)*y(k,123) + .950_r8*rxt(k,456)*y(k,199)
         mat(k,1177) = .250_r8*rxt(k,394)*y(k,123) + .250_r8*rxt(k,391)*y(k,198) &
                      + .800_r8*rxt(k,392)*y(k,199)
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
         mat(k,595) = -(rxt(k,241)*y(k,55) + rxt(k,242)*y(k,218) + rxt(k,252)*y(k,217))
         mat(k,1983) = -rxt(k,241)*y(k,42)
         mat(k,1624) = -rxt(k,242)*y(k,42)
         mat(k,1518) = -rxt(k,252)*y(k,42)
         mat(k,126) = -(rxt(k,243)*y(k,218))
         mat(k,1558) = -rxt(k,243)*y(k,43)
         mat(k,1086) = -(rxt(k,289)*y(k,125) + rxt(k,290)*y(k,218))
         mat(k,1751) = -rxt(k,289)*y(k,44)
         mat(k,1667) = -rxt(k,290)*y(k,44)
         mat(k,671) = .800_r8*rxt(k,357)*y(k,218)
         mat(k,353) = rxt(k,328)*y(k,125)
         mat(k,263) = rxt(k,285)*y(k,218)
         mat(k,344) = .500_r8*rxt(k,286)*y(k,218)
         mat(k,1026) = .500_r8*rxt(k,309)*y(k,134)
         mat(k,1328) = .100_r8*rxt(k,353)*y(k,134)
         mat(k,1896) = .400_r8*rxt(k,359)*y(k,190) + rxt(k,284)*y(k,195) &
                      + .270_r8*rxt(k,312)*y(k,196) + rxt(k,330)*y(k,201) + rxt(k,349) &
                      *y(k,214) + rxt(k,320)*y(k,220)
         mat(k,1751) = mat(k,1751) + rxt(k,328)*y(k,15)
         mat(k,2226) = .500_r8*rxt(k,309)*y(k,28) + .100_r8*rxt(k,353)*y(k,110)
         mat(k,890) = .400_r8*rxt(k,359)*y(k,123)
         mat(k,825) = rxt(k,284)*y(k,123) + 3.200_r8*rxt(k,281)*y(k,195) &
                      + .800_r8*rxt(k,282)*y(k,199)
         mat(k,790) = .270_r8*rxt(k,312)*y(k,123)
         mat(k,1806) = .800_r8*rxt(k,282)*y(k,195)
         mat(k,567) = rxt(k,330)*y(k,123)
         mat(k,2086) = .200_r8*rxt(k,348)*y(k,214)
         mat(k,679) = rxt(k,349)*y(k,123) + .200_r8*rxt(k,348)*y(k,204)
         mat(k,1667) = mat(k,1667) + .800_r8*rxt(k,357)*y(k,1) + rxt(k,285)*y(k,25) &
                      + .500_r8*rxt(k,286)*y(k,26)
         mat(k,770) = rxt(k,320)*y(k,123)
         mat(k,366) = -(rxt(k,244)*y(k,55) + rxt(k,245)*y(k,218))
         mat(k,1979) = -rxt(k,244)*y(k,45)
         mat(k,1593) = -rxt(k,245)*y(k,45)
         mat(k,103) = -(rxt(k,291)*y(k,218))
         mat(k,1556) = -rxt(k,291)*y(k,46)
         mat(k,923) = -(rxt(k,327)*y(k,218))
         mat(k,1657) = -rxt(k,327)*y(k,47)
         mat(k,670) = .800_r8*rxt(k,357)*y(k,218)
         mat(k,996) = .520_r8*rxt(k,441)*y(k,134)
         mat(k,352) = .500_r8*rxt(k,328)*y(k,125)
         mat(k,947) = .520_r8*rxt(k,444)*y(k,134)
         mat(k,1889) = .250_r8*rxt(k,359)*y(k,190) + .820_r8*rxt(k,312)*y(k,196) &
                      + .500_r8*rxt(k,330)*y(k,201) + .270_r8*rxt(k,453)*y(k,222) &
                      + .040_r8*rxt(k,458)*y(k,223)
         mat(k,1741) = .500_r8*rxt(k,328)*y(k,15)
         mat(k,2217) = .520_r8*rxt(k,441)*y(k,5) + .520_r8*rxt(k,444)*y(k,109)
         mat(k,1059) = .500_r8*rxt(k,462)*y(k,218)
         mat(k,889) = .250_r8*rxt(k,359)*y(k,123)
         mat(k,789) = .820_r8*rxt(k,312)*y(k,123) + .820_r8*rxt(k,310)*y(k,199)
         mat(k,1800) = .820_r8*rxt(k,310)*y(k,196) + .150_r8*rxt(k,451)*y(k,222) &
                      + .025_r8*rxt(k,456)*y(k,223)
         mat(k,566) = .500_r8*rxt(k,330)*y(k,123)
         mat(k,1657) = mat(k,1657) + .800_r8*rxt(k,357)*y(k,1) + .500_r8*rxt(k,462) &
                      *y(k,179)
         mat(k,1107) = .270_r8*rxt(k,453)*y(k,123) + .150_r8*rxt(k,451)*y(k,199)
         mat(k,914) = .040_r8*rxt(k,458)*y(k,123) + .025_r8*rxt(k,456)*y(k,199)
         mat(k,1222) = -(rxt(k,315)*y(k,125) + rxt(k,316)*y(k,218))
         mat(k,1761) = -rxt(k,315)*y(k,48)
         mat(k,1677) = -rxt(k,316)*y(k,48)
         mat(k,1142) = rxt(k,317)*y(k,218)
         mat(k,1211) = .880_r8*rxt(k,339)*y(k,134)
         mat(k,1331) = .500_r8*rxt(k,353)*y(k,134)
         mat(k,1906) = .170_r8*rxt(k,412)*y(k,200) + .050_r8*rxt(k,375)*y(k,207) &
                      + .250_r8*rxt(k,337)*y(k,210) + .170_r8*rxt(k,418)*y(k,213) &
                      + .400_r8*rxt(k,428)*y(k,224) + .250_r8*rxt(k,394)*y(k,226) &
                      + .540_r8*rxt(k,434)*y(k,227) + .510_r8*rxt(k,437)*y(k,229)
         mat(k,1761) = mat(k,1761) + .050_r8*rxt(k,376)*y(k,207) + .250_r8*rxt(k,336) &
                      *y(k,210) + .250_r8*rxt(k,395)*y(k,226)
         mat(k,859) = rxt(k,318)*y(k,218)
         mat(k,2234) = .880_r8*rxt(k,339)*y(k,104) + .500_r8*rxt(k,353)*y(k,110)
         mat(k,1381) = .250_r8*rxt(k,333)*y(k,210) + .250_r8*rxt(k,391)*y(k,226)
         mat(k,1815) = .240_r8*rxt(k,334)*y(k,210) + .500_r8*rxt(k,322)*y(k,221) &
                      + .100_r8*rxt(k,392)*y(k,226)
         mat(k,762) = .170_r8*rxt(k,412)*y(k,123) + .070_r8*rxt(k,411)*y(k,204)
         mat(k,2095) = .070_r8*rxt(k,411)*y(k,200) + .070_r8*rxt(k,417)*y(k,213)
         mat(k,1288) = .050_r8*rxt(k,375)*y(k,123) + .050_r8*rxt(k,376)*y(k,125)
         mat(k,1311) = .250_r8*rxt(k,337)*y(k,123) + .250_r8*rxt(k,336)*y(k,125) &
                      + .250_r8*rxt(k,333)*y(k,198) + .240_r8*rxt(k,334)*y(k,199)
         mat(k,877) = .170_r8*rxt(k,418)*y(k,123) + .070_r8*rxt(k,417)*y(k,204)
         mat(k,1677) = mat(k,1677) + rxt(k,317)*y(k,94) + rxt(k,318)*y(k,126)
         mat(k,1132) = .500_r8*rxt(k,322)*y(k,199)
         mat(k,730) = .400_r8*rxt(k,428)*y(k,123)
         mat(k,1175) = .250_r8*rxt(k,394)*y(k,123) + .250_r8*rxt(k,395)*y(k,125) &
                      + .250_r8*rxt(k,391)*y(k,198) + .100_r8*rxt(k,392)*y(k,199)
         mat(k,754) = .540_r8*rxt(k,434)*y(k,123)
         mat(k,502) = .510_r8*rxt(k,437)*y(k,123)
         mat(k,685) = -(rxt(k,296)*y(k,218))
         mat(k,1634) = -rxt(k,296)*y(k,49)
         mat(k,1020) = .120_r8*rxt(k,309)*y(k,134)
         mat(k,2211) = .120_r8*rxt(k,309)*y(k,28)
         mat(k,1371) = .100_r8*rxt(k,293)*y(k,199) + .150_r8*rxt(k,294)*y(k,204)
         mat(k,1793) = .100_r8*rxt(k,293)*y(k,198)
         mat(k,2061) = .150_r8*rxt(k,294)*y(k,198) + .150_r8*rxt(k,344)*y(k,212)
         mat(k,1350) = .150_r8*rxt(k,344)*y(k,204)
         mat(k,611) = -(rxt(k,297)*y(k,218))
         mat(k,1626) = -rxt(k,297)*y(k,50)
         mat(k,1370) = .360_r8*rxt(k,294)*y(k,204)
         mat(k,2055) = .360_r8*rxt(k,294)*y(k,198) + .400_r8*rxt(k,344)*y(k,212)
         mat(k,1349) = .400_r8*rxt(k,344)*y(k,204)
         mat(k,798) = -(rxt(k,264)*y(k,218))
         mat(k,1644) = -rxt(k,264)*y(k,51)
         mat(k,1188) = .200_r8*rxt(k,381)*y(k,199)
         mat(k,823) = .300_r8*rxt(k,282)*y(k,199)
         mat(k,1795) = .200_r8*rxt(k,381)*y(k,100) + .300_r8*rxt(k,282)*y(k,195) &
                      + 2.000_r8*rxt(k,261)*y(k,199) + .250_r8*rxt(k,367)*y(k,206) &
                      + .250_r8*rxt(k,372)*y(k,207) + .250_r8*rxt(k,334)*y(k,210) &
                      + .250_r8*rxt(k,446)*y(k,216) + .500_r8*rxt(k,322)*y(k,221) &
                      + .250_r8*rxt(k,451)*y(k,222) + .250_r8*rxt(k,456)*y(k,223) &
                      + .300_r8*rxt(k,392)*y(k,226)
         mat(k,1248) = .250_r8*rxt(k,367)*y(k,199)
         mat(k,1278) = .250_r8*rxt(k,372)*y(k,199)
         mat(k,1306) = .250_r8*rxt(k,334)*y(k,199)
         mat(k,1044) = .250_r8*rxt(k,446)*y(k,199)
         mat(k,1129) = .500_r8*rxt(k,322)*y(k,199)
         mat(k,1105) = .250_r8*rxt(k,451)*y(k,199)
         mat(k,912) = .250_r8*rxt(k,456)*y(k,199)
         mat(k,1168) = .300_r8*rxt(k,392)*y(k,199)
         mat(k,455) = -(rxt(k,265)*y(k,218))
         mat(k,1605) = -rxt(k,265)*y(k,52)
         mat(k,1791) = rxt(k,262)*y(k,204)
         mat(k,2046) = rxt(k,262)*y(k,199)
         mat(k,1431) = -(rxt(k,177)*y(k,55) + rxt(k,233)*y(k,72) + rxt(k,266)*y(k,218) &
                      + (rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,217))
         mat(k,1997) = -rxt(k,177)*y(k,53)
         mat(k,867) = -rxt(k,233)*y(k,53)
         mat(k,1687) = -rxt(k,266)*y(k,53)
         mat(k,1523) = -(rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,53)
         mat(k,1031) = .100_r8*rxt(k,309)*y(k,134)
         mat(k,2243) = .100_r8*rxt(k,309)*y(k,28)
         mat(k,386) = -(rxt(k,229)*y(k,217) + rxt(k,246)*y(k,55) + rxt(k,247)*y(k,218))
         mat(k,1516) = -rxt(k,229)*y(k,54)
         mat(k,1980) = -rxt(k,246)*y(k,54)
         mat(k,1595) = -rxt(k,247)*y(k,54)
         mat(k,2008) = -(rxt(k,176)*y(k,41) + rxt(k,177)*y(k,53) + rxt(k,178)*y(k,76) &
                      + rxt(k,179)*y(k,78) + (rxt(k,180) + rxt(k,181)) * y(k,204) &
                      + rxt(k,182)*y(k,134) + rxt(k,189)*y(k,59) + rxt(k,198)*y(k,91) &
                      + rxt(k,239)*y(k,40) + rxt(k,241)*y(k,42) + rxt(k,244)*y(k,45) &
                      + rxt(k,246)*y(k,54) + rxt(k,287)*y(k,27))
         mat(k,1493) = -rxt(k,176)*y(k,55)
         mat(k,1439) = -rxt(k,177)*y(k,55)
         mat(k,1410) = -rxt(k,178)*y(k,55)
         mat(k,606) = -rxt(k,179)*y(k,55)
         mat(k,2115) = -(rxt(k,180) + rxt(k,181)) * y(k,55)
         mat(k,2254) = -rxt(k,182)*y(k,55)
         mat(k,906) = -rxt(k,189)*y(k,55)
         mat(k,817) = -rxt(k,198)*y(k,55)
         mat(k,472) = -rxt(k,239)*y(k,55)
         mat(k,600) = -rxt(k,241)*y(k,55)
         mat(k,370) = -rxt(k,244)*y(k,55)
         mat(k,390) = -rxt(k,246)*y(k,55)
         mat(k,308) = -rxt(k,287)*y(k,55)
         mat(k,2162) = rxt(k,217)*y(k,58)
         mat(k,102) = 4.000_r8*rxt(k,201)*y(k,217)
         mat(k,146) = rxt(k,202)*y(k,217)
         mat(k,117) = 2.000_r8*rxt(k,203)*y(k,217)
         mat(k,156) = 2.000_r8*rxt(k,204)*y(k,217)
         mat(k,121) = 2.000_r8*rxt(k,205)*y(k,217)
         mat(k,161) = rxt(k,206)*y(k,217)
         mat(k,125) = 2.000_r8*rxt(k,207)*y(k,217)
         mat(k,128) = 3.000_r8*rxt(k,243)*y(k,218)
         mat(k,370) = mat(k,370) + rxt(k,245)*y(k,218)
         mat(k,1724) = rxt(k,217)*y(k,18) + (4.000_r8*rxt(k,184)+2.000_r8*rxt(k,186)) &
                      *y(k,58) + rxt(k,188)*y(k,123) + rxt(k,193)*y(k,133) &
                      + rxt(k,471)*y(k,150) + rxt(k,183)*y(k,199) + rxt(k,194) &
                      *y(k,218)
         mat(k,250) = rxt(k,238)*y(k,217)
         mat(k,246) = rxt(k,253)*y(k,217) + rxt(k,248)*y(k,218)
         mat(k,256) = rxt(k,254)*y(k,217) + rxt(k,249)*y(k,218)
         mat(k,303) = rxt(k,255)*y(k,217) + rxt(k,250)*y(k,218)
         mat(k,2138) = rxt(k,196)*y(k,133) + rxt(k,208)*y(k,217) + rxt(k,197)*y(k,218)
         mat(k,1925) = rxt(k,188)*y(k,58)
         mat(k,2193) = rxt(k,193)*y(k,58) + rxt(k,196)*y(k,84)
         mat(k,1240) = rxt(k,471)*y(k,58)
         mat(k,1833) = rxt(k,183)*y(k,58)
         mat(k,1534) = 4.000_r8*rxt(k,201)*y(k,32) + rxt(k,202)*y(k,33) &
                      + 2.000_r8*rxt(k,203)*y(k,35) + 2.000_r8*rxt(k,204)*y(k,36) &
                      + 2.000_r8*rxt(k,205)*y(k,37) + rxt(k,206)*y(k,38) &
                      + 2.000_r8*rxt(k,207)*y(k,39) + rxt(k,238)*y(k,64) + rxt(k,253) &
                      *y(k,81) + rxt(k,254)*y(k,82) + rxt(k,255)*y(k,83) + rxt(k,208) &
                      *y(k,84)
         mat(k,1698) = 3.000_r8*rxt(k,243)*y(k,43) + rxt(k,245)*y(k,45) + rxt(k,194) &
                      *y(k,58) + rxt(k,248)*y(k,81) + rxt(k,249)*y(k,82) + rxt(k,250) &
                      *y(k,83) + rxt(k,197)*y(k,84)
         mat(k,1976) = rxt(k,189)*y(k,59)
         mat(k,1706) = 2.000_r8*rxt(k,185)*y(k,58)
         mat(k,898) = rxt(k,189)*y(k,55) + (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,84)
         mat(k,2122) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,59) + (rxt(k,524) &
                       +rxt(k,530)+rxt(k,535))*y(k,91)
         mat(k,813) = (rxt(k,524)+rxt(k,530)+rxt(k,535))*y(k,84)
         mat(k,1705) = 2.000_r8*rxt(k,210)*y(k,58)
         mat(k,1719) = -(rxt(k,183)*y(k,199) + (4._r8*rxt(k,184) + 4._r8*rxt(k,185) &
                      + 4._r8*rxt(k,186) + 4._r8*rxt(k,210)) * y(k,58) + rxt(k,187) &
                      *y(k,204) + rxt(k,188)*y(k,123) + rxt(k,190)*y(k,124) + rxt(k,193) &
                      *y(k,133) + (rxt(k,194) + rxt(k,195)) * y(k,218) + (rxt(k,216) &
                      + rxt(k,217) + rxt(k,218)) * y(k,18) + rxt(k,471)*y(k,150))
         mat(k,1828) = -rxt(k,183)*y(k,58)
         mat(k,2110) = -rxt(k,187)*y(k,58)
         mat(k,1920) = -rxt(k,188)*y(k,58)
         mat(k,1964) = -rxt(k,190)*y(k,58)
         mat(k,2188) = -rxt(k,193)*y(k,58)
         mat(k,1693) = -(rxt(k,194) + rxt(k,195)) * y(k,58)
         mat(k,2157) = -(rxt(k,216) + rxt(k,217) + rxt(k,218)) * y(k,58)
         mat(k,1237) = -rxt(k,471)*y(k,58)
         mat(k,2003) = rxt(k,198)*y(k,91) + rxt(k,182)*y(k,134) + rxt(k,181)*y(k,204)
         mat(k,903) = rxt(k,191)*y(k,133)
         mat(k,2133) = rxt(k,209)*y(k,217)
         mat(k,816) = rxt(k,198)*y(k,55) + rxt(k,199)*y(k,133) + rxt(k,200)*y(k,218)
         mat(k,2188) = mat(k,2188) + rxt(k,191)*y(k,59) + rxt(k,199)*y(k,91)
         mat(k,2249) = rxt(k,182)*y(k,55)
         mat(k,329) = rxt(k,476)*y(k,150)
         mat(k,1237) = mat(k,1237) + rxt(k,476)*y(k,136)
         mat(k,2110) = mat(k,2110) + rxt(k,181)*y(k,55)
         mat(k,1529) = rxt(k,209)*y(k,84)
         mat(k,1693) = mat(k,1693) + rxt(k,200)*y(k,91)
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
         mat(k,900) = -(rxt(k,189)*y(k,55) + rxt(k,191)*y(k,133) + rxt(k,192)*y(k,218) &
                      + (rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,84))
         mat(k,1988) = -rxt(k,189)*y(k,59)
         mat(k,2177) = -rxt(k,191)*y(k,59)
         mat(k,1655) = -rxt(k,192)*y(k,59)
         mat(k,2126) = -(rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,59)
         mat(k,1711) = rxt(k,190)*y(k,124)
         mat(k,1949) = rxt(k,190)*y(k,58)
         mat(k,1124) = -(rxt(k,276)*y(k,218))
         mat(k,1670) = -rxt(k,276)*y(k,61)
         mat(k,1005) = .230_r8*rxt(k,441)*y(k,134)
         mat(k,1415) = rxt(k,212)*y(k,41)
         mat(k,295) = .350_r8*rxt(k,278)*y(k,218)
         mat(k,552) = .630_r8*rxt(k,280)*y(k,134)
         mat(k,1027) = .560_r8*rxt(k,309)*y(k,134)
         mat(k,1480) = rxt(k,212)*y(k,16) + rxt(k,176)*y(k,55) + rxt(k,257)*y(k,125) &
                      + rxt(k,258)*y(k,133) + rxt(k,259)*y(k,218)
         mat(k,367) = rxt(k,244)*y(k,55)
         mat(k,1221) = rxt(k,315)*y(k,125) + rxt(k,316)*y(k,218)
         mat(k,1993) = rxt(k,176)*y(k,41) + rxt(k,244)*y(k,45)
         mat(k,981) = rxt(k,303)*y(k,218)
         mat(k,838) = .620_r8*rxt(k,386)*y(k,134)
         mat(k,1209) = .650_r8*rxt(k,339)*y(k,134)
         mat(k,955) = .230_r8*rxt(k,444)*y(k,134)
         mat(k,1329) = .560_r8*rxt(k,353)*y(k,134)
         mat(k,1899) = .170_r8*rxt(k,412)*y(k,200) + .220_r8*rxt(k,337)*y(k,210) &
                      + .400_r8*rxt(k,415)*y(k,211) + .350_r8*rxt(k,418)*y(k,213) &
                      + .225_r8*rxt(k,453)*y(k,222) + .250_r8*rxt(k,394)*y(k,226)
         mat(k,1754) = rxt(k,257)*y(k,41) + rxt(k,315)*y(k,48) + .220_r8*rxt(k,336) &
                      *y(k,210) + .500_r8*rxt(k,395)*y(k,226)
         mat(k,2178) = rxt(k,258)*y(k,41) + rxt(k,465)*y(k,137)
         mat(k,2229) = .230_r8*rxt(k,441)*y(k,5) + .630_r8*rxt(k,280)*y(k,24) &
                      + .560_r8*rxt(k,309)*y(k,28) + .620_r8*rxt(k,386)*y(k,97) &
                      + .650_r8*rxt(k,339)*y(k,104) + .230_r8*rxt(k,444)*y(k,109) &
                      + .560_r8*rxt(k,353)*y(k,110)
         mat(k,361) = rxt(k,465)*y(k,133) + rxt(k,466)*y(k,218)
         mat(k,1061) = .700_r8*rxt(k,462)*y(k,218)
         mat(k,1376) = .220_r8*rxt(k,333)*y(k,210) + .250_r8*rxt(k,391)*y(k,226)
         mat(k,1809) = .110_r8*rxt(k,334)*y(k,210) + .125_r8*rxt(k,451)*y(k,222) &
                      + .200_r8*rxt(k,392)*y(k,226)
         mat(k,761) = .170_r8*rxt(k,412)*y(k,123) + .070_r8*rxt(k,411)*y(k,204)
         mat(k,2089) = .070_r8*rxt(k,411)*y(k,200) + .160_r8*rxt(k,414)*y(k,211) &
                      + .140_r8*rxt(k,417)*y(k,213)
         mat(k,1308) = .220_r8*rxt(k,337)*y(k,123) + .220_r8*rxt(k,336)*y(k,125) &
                      + .220_r8*rxt(k,333)*y(k,198) + .110_r8*rxt(k,334)*y(k,199)
         mat(k,716) = .400_r8*rxt(k,415)*y(k,123) + .160_r8*rxt(k,414)*y(k,204)
         mat(k,876) = .350_r8*rxt(k,418)*y(k,123) + .140_r8*rxt(k,417)*y(k,204)
         mat(k,1670) = mat(k,1670) + .350_r8*rxt(k,278)*y(k,23) + rxt(k,259)*y(k,41) &
                      + rxt(k,316)*y(k,48) + rxt(k,303)*y(k,74) + rxt(k,466)*y(k,137) &
                      + .700_r8*rxt(k,462)*y(k,179)
         mat(k,1111) = .225_r8*rxt(k,453)*y(k,123) + .125_r8*rxt(k,451)*y(k,199)
         mat(k,1172) = .250_r8*rxt(k,394)*y(k,123) + .500_r8*rxt(k,395)*y(k,125) &
                      + .250_r8*rxt(k,391)*y(k,198) + .200_r8*rxt(k,392)*y(k,199)
         mat(k,993) = .270_r8*rxt(k,441)*y(k,134)
         mat(k,1022) = .200_r8*rxt(k,309)*y(k,134)
         mat(k,686) = rxt(k,296)*y(k,218)
         mat(k,612) = .500_r8*rxt(k,297)*y(k,218)
         mat(k,1123) = rxt(k,276)*y(k,218)
         mat(k,1095) = .800_r8*rxt(k,302)*y(k,218)
         mat(k,979) = rxt(k,303)*y(k,218)
         mat(k,929) = rxt(k,268)*y(k,218)
         mat(k,574) = .500_r8*rxt(k,352)*y(k,218)
         mat(k,944) = .270_r8*rxt(k,444)*y(k,134)
         mat(k,1325) = .100_r8*rxt(k,353)*y(k,134)
         mat(k,1883) = rxt(k,295)*y(k,198) + .900_r8*rxt(k,453)*y(k,222)
         mat(k,2213) = .270_r8*rxt(k,441)*y(k,5) + .200_r8*rxt(k,309)*y(k,28) &
                      + .270_r8*rxt(k,444)*y(k,109) + .100_r8*rxt(k,353)*y(k,110)
         mat(k,1058) = 1.800_r8*rxt(k,462)*y(k,218)
         mat(k,1372) = rxt(k,295)*y(k,123) + 4.000_r8*rxt(k,292)*y(k,198) &
                      + .900_r8*rxt(k,293)*y(k,199) + .490_r8*rxt(k,294)*y(k,204) &
                      + rxt(k,366)*y(k,206) + 2.000_r8*rxt(k,342)*y(k,212) &
                      + rxt(k,391)*y(k,226)
         mat(k,1796) = .900_r8*rxt(k,293)*y(k,198) + rxt(k,343)*y(k,212) &
                      + .500_r8*rxt(k,451)*y(k,222)
         mat(k,2072) = .490_r8*rxt(k,294)*y(k,198) + .450_r8*rxt(k,344)*y(k,212)
         mat(k,1249) = rxt(k,366)*y(k,198)
         mat(k,1351) = 2.000_r8*rxt(k,342)*y(k,198) + rxt(k,343)*y(k,199) &
                      + .450_r8*rxt(k,344)*y(k,204) + 4.000_r8*rxt(k,345)*y(k,212)
         mat(k,1645) = rxt(k,296)*y(k,49) + .500_r8*rxt(k,297)*y(k,50) + rxt(k,276) &
                      *y(k,61) + .800_r8*rxt(k,302)*y(k,73) + rxt(k,303)*y(k,74) &
                      + rxt(k,268)*y(k,86) + .500_r8*rxt(k,352)*y(k,108) &
                      + 1.800_r8*rxt(k,462)*y(k,179)
         mat(k,1106) = .900_r8*rxt(k,453)*y(k,123) + .500_r8*rxt(k,451)*y(k,199)
         mat(k,1169) = rxt(k,391)*y(k,198)
         mat(k,239) = -(rxt(k,237)*y(k,217))
         mat(k,1511) = -rxt(k,237)*y(k,63)
         mat(k,143) = rxt(k,202)*y(k,217)
         mat(k,148) = rxt(k,228)*y(k,217)
         mat(k,153) = rxt(k,204)*y(k,217)
         mat(k,119) = 2.000_r8*rxt(k,205)*y(k,217)
         mat(k,158) = 2.000_r8*rxt(k,206)*y(k,217)
         mat(k,123) = rxt(k,207)*y(k,217)
         mat(k,107) = 2.000_r8*rxt(k,230)*y(k,217)
         mat(k,251) = rxt(k,254)*y(k,217) + rxt(k,249)*y(k,218)
         mat(k,298) = rxt(k,255)*y(k,217) + rxt(k,250)*y(k,218)
         mat(k,1511) = mat(k,1511) + rxt(k,202)*y(k,33) + rxt(k,228)*y(k,34) &
                      + rxt(k,204)*y(k,36) + 2.000_r8*rxt(k,205)*y(k,37) &
                      + 2.000_r8*rxt(k,206)*y(k,38) + rxt(k,207)*y(k,39) &
                      + 2.000_r8*rxt(k,230)*y(k,77) + rxt(k,254)*y(k,82) + rxt(k,255) &
                      *y(k,83)
         mat(k,1574) = rxt(k,249)*y(k,82) + rxt(k,250)*y(k,83)
         mat(k,247) = -(rxt(k,238)*y(k,217))
         mat(k,1513) = -rxt(k,238)*y(k,64)
         mat(k,115) = rxt(k,203)*y(k,217)
         mat(k,154) = rxt(k,204)*y(k,217)
         mat(k,243) = rxt(k,253)*y(k,217) + rxt(k,248)*y(k,218)
         mat(k,1513) = mat(k,1513) + rxt(k,203)*y(k,35) + rxt(k,204)*y(k,36) &
                      + rxt(k,253)*y(k,81)
         mat(k,1576) = rxt(k,248)*y(k,81)
         mat(k,195) = -(rxt(k,410)*y(k,218))
         mat(k,1565) = -rxt(k,410)*y(k,65)
         mat(k,189) = .180_r8*rxt(k,430)*y(k,218)
         mat(k,1565) = mat(k,1565) + .180_r8*rxt(k,430)*y(k,181)
         mat(k,283) = -(rxt(k,463)*y(k,125) + (rxt(k,464) + rxt(k,478)) * y(k,218))
         mat(k,1734) = -rxt(k,463)*y(k,66)
         mat(k,1581) = -(rxt(k,464) + rxt(k,478)) * y(k,66)
         mat(k,705) = rxt(k,298)*y(k,204)
         mat(k,2027) = rxt(k,298)*y(k,203)
         mat(k,865) = -(rxt(k,233)*y(k,53) + rxt(k,234)*y(k,76) + rxt(k,235)*y(k,230) &
                      + rxt(k,236)*y(k,88))
         mat(k,1428) = -rxt(k,233)*y(k,72)
         mat(k,1401) = -rxt(k,234)*y(k,72)
         mat(k,2265) = -rxt(k,235)*y(k,72)
         mat(k,1460) = -rxt(k,236)*y(k,72)
         mat(k,149) = rxt(k,228)*y(k,217)
         mat(k,159) = rxt(k,206)*y(k,217)
         mat(k,240) = 2.000_r8*rxt(k,237)*y(k,217)
         mat(k,248) = rxt(k,238)*y(k,217)
         mat(k,1520) = rxt(k,228)*y(k,34) + rxt(k,206)*y(k,38) + 2.000_r8*rxt(k,237) &
                      *y(k,63) + rxt(k,238)*y(k,64)
         mat(k,1097) = -(rxt(k,302)*y(k,218))
         mat(k,1668) = -rxt(k,302)*y(k,73)
         mat(k,583) = .700_r8*rxt(k,377)*y(k,218)
         mat(k,535) = .500_r8*rxt(k,378)*y(k,218)
         mat(k,376) = rxt(k,389)*y(k,218)
         mat(k,1897) = .050_r8*rxt(k,375)*y(k,207) + .530_r8*rxt(k,337)*y(k,210) &
                      + .225_r8*rxt(k,453)*y(k,222) + .250_r8*rxt(k,394)*y(k,226)
         mat(k,1752) = .050_r8*rxt(k,376)*y(k,207) + .530_r8*rxt(k,336)*y(k,210) &
                      + .250_r8*rxt(k,395)*y(k,226)
         mat(k,1375) = .530_r8*rxt(k,333)*y(k,210) + .250_r8*rxt(k,391)*y(k,226)
         mat(k,1807) = .260_r8*rxt(k,334)*y(k,210) + .125_r8*rxt(k,451)*y(k,222) &
                      + .100_r8*rxt(k,392)*y(k,226)
         mat(k,1282) = .050_r8*rxt(k,375)*y(k,123) + .050_r8*rxt(k,376)*y(k,125)
         mat(k,1307) = .530_r8*rxt(k,337)*y(k,123) + .530_r8*rxt(k,336)*y(k,125) &
                      + .530_r8*rxt(k,333)*y(k,198) + .260_r8*rxt(k,334)*y(k,199)
         mat(k,1668) = mat(k,1668) + .700_r8*rxt(k,377)*y(k,98) + .500_r8*rxt(k,378) &
                      *y(k,99) + rxt(k,389)*y(k,114)
         mat(k,1109) = .225_r8*rxt(k,453)*y(k,123) + .125_r8*rxt(k,451)*y(k,199)
         mat(k,1171) = .250_r8*rxt(k,394)*y(k,123) + .250_r8*rxt(k,395)*y(k,125) &
                      + .250_r8*rxt(k,391)*y(k,198) + .100_r8*rxt(k,392)*y(k,199)
         mat(k,980) = -(rxt(k,303)*y(k,218))
         mat(k,1661) = -rxt(k,303)*y(k,74)
         mat(k,294) = .650_r8*rxt(k,278)*y(k,218)
         mat(k,1096) = .200_r8*rxt(k,302)*y(k,218)
         mat(k,1073) = rxt(k,390)*y(k,218)
         mat(k,1892) = rxt(k,401)*y(k,192) + .050_r8*rxt(k,375)*y(k,207) &
                      + .400_r8*rxt(k,415)*y(k,211) + .170_r8*rxt(k,418)*y(k,213) &
                      + .700_r8*rxt(k,421)*y(k,219) + .600_r8*rxt(k,428)*y(k,224) &
                      + .250_r8*rxt(k,394)*y(k,226) + .340_r8*rxt(k,434)*y(k,227) &
                      + .170_r8*rxt(k,437)*y(k,229)
         mat(k,1745) = .050_r8*rxt(k,376)*y(k,207) + .250_r8*rxt(k,395)*y(k,226)
         mat(k,486) = rxt(k,401)*y(k,123)
         mat(k,1373) = .250_r8*rxt(k,391)*y(k,226)
         mat(k,1802) = .100_r8*rxt(k,392)*y(k,226)
         mat(k,2083) = .160_r8*rxt(k,414)*y(k,211) + .070_r8*rxt(k,417)*y(k,213)
         mat(k,1281) = .050_r8*rxt(k,375)*y(k,123) + .050_r8*rxt(k,376)*y(k,125)
         mat(k,715) = .400_r8*rxt(k,415)*y(k,123) + .160_r8*rxt(k,414)*y(k,204)
         mat(k,875) = .170_r8*rxt(k,418)*y(k,123) + .070_r8*rxt(k,417)*y(k,204)
         mat(k,1661) = mat(k,1661) + .650_r8*rxt(k,278)*y(k,23) + .200_r8*rxt(k,302) &
                      *y(k,73) + rxt(k,390)*y(k,115)
         mat(k,450) = .700_r8*rxt(k,421)*y(k,123)
         mat(k,728) = .600_r8*rxt(k,428)*y(k,123)
         mat(k,1170) = .250_r8*rxt(k,394)*y(k,123) + .250_r8*rxt(k,395)*y(k,125) &
                      + .250_r8*rxt(k,391)*y(k,198) + .100_r8*rxt(k,392)*y(k,199)
         mat(k,752) = .340_r8*rxt(k,434)*y(k,123)
         mat(k,501) = .170_r8*rxt(k,437)*y(k,123)
         mat(k,1446) = -((rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,204) + rxt(k,142) &
                      *y(k,134))
         mat(k,2105) = -(rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,75)
         mat(k,2244) = -rxt(k,142)*y(k,75)
         mat(k,1484) = rxt(k,259)*y(k,218)
         mat(k,1432) = rxt(k,273)*y(k,217)
         mat(k,1998) = rxt(k,178)*y(k,76)
         mat(k,868) = rxt(k,234)*y(k,76)
         mat(k,1404) = rxt(k,178)*y(k,55) + rxt(k,234)*y(k,72) + rxt(k,134)*y(k,133) &
                      + rxt(k,125)*y(k,217) + rxt(k,143)*y(k,218)
         mat(k,807) = rxt(k,232)*y(k,217)
         mat(k,2128) = rxt(k,209)*y(k,217)
         mat(k,493) = rxt(k,164)*y(k,218)
         mat(k,2183) = rxt(k,134)*y(k,76) + rxt(k,146)*y(k,218)
         mat(k,363) = rxt(k,466)*y(k,218)
         mat(k,514) = rxt(k,472)*y(k,218)
         mat(k,1235) = rxt(k,477)*y(k,218)
         mat(k,1524) = rxt(k,273)*y(k,53) + rxt(k,125)*y(k,76) + rxt(k,232)*y(k,80) &
                      + rxt(k,209)*y(k,84)
         mat(k,1688) = rxt(k,259)*y(k,41) + rxt(k,143)*y(k,76) + rxt(k,164)*y(k,111) &
                      + rxt(k,146)*y(k,133) + rxt(k,466)*y(k,137) + rxt(k,472) &
                      *y(k,148) + rxt(k,477)*y(k,150)
         mat(k,1402) = -(rxt(k,125)*y(k,217) + rxt(k,134)*y(k,133) + rxt(k,143) &
                      *y(k,218) + rxt(k,178)*y(k,55) + rxt(k,234)*y(k,72))
         mat(k,1521) = -rxt(k,125)*y(k,76)
         mat(k,2180) = -rxt(k,134)*y(k,76)
         mat(k,1685) = -rxt(k,143)*y(k,76)
         mat(k,1995) = -rxt(k,178)*y(k,76)
         mat(k,866) = -rxt(k,234)*y(k,76)
         mat(k,1430) = rxt(k,274)*y(k,217)
         mat(k,1444) = rxt(k,136)*y(k,204)
         mat(k,2102) = rxt(k,136)*y(k,75)
         mat(k,1521) = mat(k,1521) + rxt(k,274)*y(k,53)
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
         mat(k,106) = -(rxt(k,230)*y(k,217))
         mat(k,1501) = -rxt(k,230)*y(k,77)
         mat(k,604) = -(rxt(k,135)*y(k,133) + rxt(k,144)*y(k,218) + rxt(k,179)*y(k,55))
         mat(k,2172) = -rxt(k,135)*y(k,78)
         mat(k,1625) = -rxt(k,144)*y(k,78)
         mat(k,1984) = -rxt(k,179)*y(k,78)
         mat(k,2054) = 2.000_r8*rxt(k,150)*y(k,204)
         mat(k,1625) = mat(k,1625) + 2.000_r8*rxt(k,149)*y(k,218)
         mat(k,257) = rxt(k,479)*y(k,230)
         mat(k,2261) = rxt(k,479)*y(k,152)
         mat(k,805) = -(rxt(k,225)*y(k,133) + rxt(k,226)*y(k,218) + (rxt(k,231) &
                      + rxt(k,232)) * y(k,217))
         mat(k,2174) = -rxt(k,225)*y(k,80)
         mat(k,1646) = -rxt(k,226)*y(k,80)
         mat(k,1519) = -(rxt(k,231) + rxt(k,232)) * y(k,80)
         mat(k,1414) = rxt(k,212)*y(k,41) + rxt(k,213)*y(k,204)
         mat(k,1478) = rxt(k,212)*y(k,16)
         mat(k,2073) = rxt(k,213)*y(k,16)
         mat(k,242) = -(rxt(k,248)*y(k,218) + rxt(k,253)*y(k,217))
         mat(k,1575) = -rxt(k,248)*y(k,81)
         mat(k,1512) = -rxt(k,253)*y(k,81)
         mat(k,252) = -(rxt(k,249)*y(k,218) + rxt(k,254)*y(k,217))
         mat(k,1577) = -rxt(k,249)*y(k,82)
         mat(k,1514) = -rxt(k,254)*y(k,82)
         mat(k,299) = -(rxt(k,250)*y(k,218) + rxt(k,255)*y(k,217))
         mat(k,1583) = -rxt(k,250)*y(k,83)
         mat(k,1515) = -rxt(k,255)*y(k,83)
         mat(k,2140) = -(rxt(k,196)*y(k,133) + rxt(k,197)*y(k,218) + (rxt(k,208) &
                      + rxt(k,209)) * y(k,217) + (rxt(k,524) + rxt(k,530) + rxt(k,535) &
                      ) * y(k,91) + (rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,59) &
                      + (rxt(k,531) + rxt(k,536)) * y(k,90))
         mat(k,2195) = -rxt(k,196)*y(k,84)
         mat(k,1700) = -rxt(k,197)*y(k,84)
         mat(k,1536) = -(rxt(k,208) + rxt(k,209)) * y(k,84)
         mat(k,818) = -(rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,84)
         mat(k,907) = -(rxt(k,529) + rxt(k,534) + rxt(k,539)) * y(k,84)
         mat(k,783) = -(rxt(k,531) + rxt(k,536)) * y(k,84)
         mat(k,309) = rxt(k,287)*y(k,55)
         mat(k,473) = rxt(k,239)*y(k,55)
         mat(k,1495) = rxt(k,176)*y(k,55)
         mat(k,602) = rxt(k,241)*y(k,55)
         mat(k,372) = 2.000_r8*rxt(k,244)*y(k,55)
         mat(k,1441) = rxt(k,177)*y(k,55)
         mat(k,391) = rxt(k,246)*y(k,55)
         mat(k,2010) = rxt(k,287)*y(k,27) + rxt(k,239)*y(k,40) + rxt(k,176)*y(k,41) &
                      + rxt(k,241)*y(k,42) + 2.000_r8*rxt(k,244)*y(k,45) + rxt(k,177) &
                      *y(k,53) + rxt(k,246)*y(k,54) + rxt(k,178)*y(k,76) + rxt(k,179) &
                      *y(k,78) + rxt(k,198)*y(k,91) + rxt(k,180)*y(k,204)
         mat(k,1726) = rxt(k,195)*y(k,218)
         mat(k,1411) = rxt(k,178)*y(k,55)
         mat(k,608) = rxt(k,179)*y(k,55)
         mat(k,818) = mat(k,818) + rxt(k,198)*y(k,55)
         mat(k,2117) = rxt(k,180)*y(k,55)
         mat(k,1700) = mat(k,1700) + rxt(k,195)*y(k,58)
         mat(k,180) = -(rxt(k,267)*y(k,218) + rxt(k,275)*y(k,217))
         mat(k,1563) = -rxt(k,267)*y(k,85)
         mat(k,1509) = -rxt(k,275)*y(k,85)
         mat(k,930) = -(rxt(k,268)*y(k,218))
         mat(k,1658) = -rxt(k,268)*y(k,86)
         mat(k,997) = .050_r8*rxt(k,441)*y(k,134)
         mat(k,293) = .350_r8*rxt(k,278)*y(k,218)
         mat(k,551) = .370_r8*rxt(k,280)*y(k,134)
         mat(k,1024) = .120_r8*rxt(k,309)*y(k,134)
         mat(k,836) = .110_r8*rxt(k,386)*y(k,134)
         mat(k,1208) = .330_r8*rxt(k,339)*y(k,134)
         mat(k,948) = .050_r8*rxt(k,444)*y(k,134)
         mat(k,1326) = .120_r8*rxt(k,353)*y(k,134)
         mat(k,1890) = rxt(k,271)*y(k,205)
         mat(k,2218) = .050_r8*rxt(k,441)*y(k,5) + .370_r8*rxt(k,280)*y(k,24) &
                      + .120_r8*rxt(k,309)*y(k,28) + .110_r8*rxt(k,386)*y(k,97) &
                      + .330_r8*rxt(k,339)*y(k,104) + .050_r8*rxt(k,444)*y(k,109) &
                      + .120_r8*rxt(k,353)*y(k,110)
         mat(k,2081) = rxt(k,269)*y(k,205)
         mat(k,443) = rxt(k,271)*y(k,123) + rxt(k,269)*y(k,204)
         mat(k,1658) = mat(k,1658) + .350_r8*rxt(k,278)*y(k,23)
         mat(k,1426) = rxt(k,233)*y(k,72)
         mat(k,864) = rxt(k,233)*y(k,53) + rxt(k,234)*y(k,76) + rxt(k,236)*y(k,88) &
                      + rxt(k,235)*y(k,230)
         mat(k,1400) = rxt(k,234)*y(k,72)
         mat(k,1459) = rxt(k,236)*y(k,72)
         mat(k,2263) = rxt(k,235)*y(k,72)
         mat(k,1464) = -(rxt(k,173)*y(k,218) + rxt(k,236)*y(k,72))
         mat(k,1689) = -rxt(k,173)*y(k,88)
         mat(k,869) = -rxt(k,236)*y(k,88)
         mat(k,1485) = rxt(k,257)*y(k,125)
         mat(k,1089) = rxt(k,289)*y(k,125)
         mat(k,1224) = rxt(k,315)*y(k,125)
         mat(k,901) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,84)
         mat(k,285) = rxt(k,463)*y(k,125)
         mat(k,2129) = (rxt(k,529)+rxt(k,534)+rxt(k,539))*y(k,59)
         mat(k,1960) = rxt(k,172)*y(k,218)
         mat(k,1772) = rxt(k,257)*y(k,41) + rxt(k,289)*y(k,44) + rxt(k,315)*y(k,48) &
                      + rxt(k,463)*y(k,66)
         mat(k,1689) = mat(k,1689) + rxt(k,172)*y(k,124)
         mat(k,422) = -(rxt(k,151)*y(k,218))
         mat(k,1601) = -rxt(k,151)*y(k,89)
         mat(k,1935) = rxt(k,170)*y(k,204)
         mat(k,2041) = rxt(k,170)*y(k,124)
         mat(k,779) = -(rxt(k,227)*y(k,133) + (rxt(k,531) + rxt(k,536)) * y(k,84))
         mat(k,2173) = -rxt(k,227)*y(k,90)
         mat(k,2124) = -(rxt(k,531) + rxt(k,536)) * y(k,90)
         mat(k,2148) = rxt(k,219)*y(k,204)
         mat(k,2070) = rxt(k,219)*y(k,18)
         mat(k,814) = -(rxt(k,198)*y(k,55) + rxt(k,199)*y(k,133) + rxt(k,200)*y(k,218) &
                      + (rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,84))
         mat(k,1986) = -rxt(k,198)*y(k,91)
         mat(k,2175) = -rxt(k,199)*y(k,91)
         mat(k,1647) = -rxt(k,200)*y(k,91)
         mat(k,2125) = -(rxt(k,524) + rxt(k,530) + rxt(k,535)) * y(k,91)
         mat(k,1709) = rxt(k,187)*y(k,204)
         mat(k,899) = rxt(k,192)*y(k,218)
         mat(k,2074) = rxt(k,187)*y(k,58)
         mat(k,1647) = mat(k,1647) + rxt(k,192)*y(k,59)
         mat(k,1154) = -(rxt(k,332)*y(k,218))
         mat(k,1673) = -rxt(k,332)*y(k,92)
         mat(k,585) = .300_r8*rxt(k,377)*y(k,218)
         mat(k,537) = .500_r8*rxt(k,378)*y(k,218)
         mat(k,1902) = rxt(k,331)*y(k,201) + rxt(k,338)*y(k,210)
         mat(k,568) = rxt(k,331)*y(k,123)
         mat(k,1310) = rxt(k,338)*y(k,123)
         mat(k,1673) = mat(k,1673) + .300_r8*rxt(k,377)*y(k,98) + .500_r8*rxt(k,378) &
                      *y(k,99)
         mat(k,222) = -(rxt(k,363)*y(k,218))
         mat(k,1570) = -rxt(k,363)*y(k,93)
         mat(k,1141) = -(rxt(k,317)*y(k,218))
         mat(k,1672) = -rxt(k,317)*y(k,94)
         mat(k,584) = .700_r8*rxt(k,377)*y(k,218)
         mat(k,536) = .500_r8*rxt(k,378)*y(k,218)
         mat(k,575) = .500_r8*rxt(k,352)*y(k,218)
         mat(k,1901) = .050_r8*rxt(k,375)*y(k,207) + .220_r8*rxt(k,337)*y(k,210) &
                      + .250_r8*rxt(k,394)*y(k,226)
         mat(k,1756) = .050_r8*rxt(k,376)*y(k,207) + .220_r8*rxt(k,336)*y(k,210) &
                      + .250_r8*rxt(k,395)*y(k,226)
         mat(k,544) = .500_r8*rxt(k,321)*y(k,218)
         mat(k,1377) = .220_r8*rxt(k,333)*y(k,210) + .250_r8*rxt(k,391)*y(k,226)
         mat(k,1811) = .230_r8*rxt(k,334)*y(k,210) + .200_r8*rxt(k,322)*y(k,221) &
                      + .100_r8*rxt(k,392)*y(k,226)
         mat(k,1284) = .050_r8*rxt(k,375)*y(k,123) + .050_r8*rxt(k,376)*y(k,125)
         mat(k,1309) = .220_r8*rxt(k,337)*y(k,123) + .220_r8*rxt(k,336)*y(k,125) &
                      + .220_r8*rxt(k,333)*y(k,198) + .230_r8*rxt(k,334)*y(k,199)
         mat(k,1672) = mat(k,1672) + .700_r8*rxt(k,377)*y(k,98) + .500_r8*rxt(k,378) &
                      *y(k,99) + .500_r8*rxt(k,352)*y(k,108) + .500_r8*rxt(k,321) &
                      *y(k,146)
         mat(k,1131) = .200_r8*rxt(k,322)*y(k,199)
         mat(k,1173) = .250_r8*rxt(k,394)*y(k,123) + .250_r8*rxt(k,395)*y(k,125) &
                      + .250_r8*rxt(k,391)*y(k,198) + .100_r8*rxt(k,392)*y(k,199)
         mat(k,347) = -(rxt(k,364)*y(k,218))
         mat(k,1590) = -rxt(k,364)*y(k,95)
         mat(k,1857) = .870_r8*rxt(k,375)*y(k,207)
         mat(k,1735) = .950_r8*rxt(k,376)*y(k,207)
         mat(k,1368) = rxt(k,371)*y(k,207)
         mat(k,1789) = .750_r8*rxt(k,372)*y(k,207)
         mat(k,1274) = .870_r8*rxt(k,375)*y(k,123) + .950_r8*rxt(k,376)*y(k,125) &
                      + rxt(k,371)*y(k,198) + .750_r8*rxt(k,372)*y(k,199)
         mat(k,136) = -(rxt(k,365)*y(k,218))
         mat(k,1559) = -rxt(k,365)*y(k,96)
         mat(k,735) = .600_r8*rxt(k,388)*y(k,218)
         mat(k,1559) = mat(k,1559) + .600_r8*rxt(k,388)*y(k,102)
         mat(k,835) = -(rxt(k,379)*y(k,125) + rxt(k,386)*y(k,134) + rxt(k,387) &
                      *y(k,218))
         mat(k,1738) = -rxt(k,379)*y(k,97)
         mat(k,2214) = -rxt(k,386)*y(k,97)
         mat(k,1649) = -rxt(k,387)*y(k,97)
         mat(k,582) = -(rxt(k,377)*y(k,218))
         mat(k,1622) = -rxt(k,377)*y(k,98)
         mat(k,1870) = .080_r8*rxt(k,369)*y(k,206)
         mat(k,1246) = .080_r8*rxt(k,369)*y(k,123)
         mat(k,533) = -(rxt(k,378)*y(k,218))
         mat(k,1617) = -rxt(k,378)*y(k,99)
         mat(k,1868) = .080_r8*rxt(k,375)*y(k,207)
         mat(k,1275) = .080_r8*rxt(k,375)*y(k,123)
         mat(k,1194) = -(rxt(k,380)*y(k,198) + rxt(k,381)*y(k,199) + rxt(k,382) &
                      *y(k,204) + rxt(k,383)*y(k,123) + rxt(k,384)*y(k,125))
         mat(k,1379) = -rxt(k,380)*y(k,100)
         mat(k,1813) = -rxt(k,381)*y(k,100)
         mat(k,2093) = -rxt(k,382)*y(k,100)
         mat(k,1904) = -rxt(k,383)*y(k,100)
         mat(k,1759) = -rxt(k,384)*y(k,100)
         mat(k,839) = rxt(k,379)*y(k,125)
         mat(k,1759) = mat(k,1759) + rxt(k,379)*y(k,97)
         mat(k,404) = -(rxt(k,385)*y(k,218))
         mat(k,1598) = -rxt(k,385)*y(k,101)
         mat(k,1186) = rxt(k,382)*y(k,204)
         mat(k,2038) = rxt(k,382)*y(k,100)
         mat(k,736) = -(rxt(k,388)*y(k,218))
         mat(k,1639) = -rxt(k,388)*y(k,102)
         mat(k,2066) = rxt(k,368)*y(k,206) + rxt(k,373)*y(k,207)
         mat(k,1247) = rxt(k,368)*y(k,204)
         mat(k,1277) = rxt(k,373)*y(k,204)
         mat(k,75) = -(rxt(k,509)*y(k,218))
         mat(k,1551) = -rxt(k,509)*y(k,103)
         mat(k,1210) = -(rxt(k,339)*y(k,134) + rxt(k,340)*y(k,218))
         mat(k,2233) = -rxt(k,339)*y(k,104)
         mat(k,1676) = -rxt(k,340)*y(k,104)
         mat(k,840) = .300_r8*rxt(k,386)*y(k,134)
         mat(k,1905) = .360_r8*rxt(k,369)*y(k,206)
         mat(k,1760) = .400_r8*rxt(k,370)*y(k,206)
         mat(k,2233) = mat(k,2233) + .300_r8*rxt(k,386)*y(k,97)
         mat(k,1380) = .390_r8*rxt(k,366)*y(k,206)
         mat(k,1814) = .310_r8*rxt(k,367)*y(k,206)
         mat(k,1255) = .360_r8*rxt(k,369)*y(k,123) + .400_r8*rxt(k,370)*y(k,125) &
                      + .390_r8*rxt(k,366)*y(k,198) + .310_r8*rxt(k,367)*y(k,199)
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
         mat(k,311) = -(rxt(k,341)*y(k,218))
         mat(k,1585) = -rxt(k,341)*y(k,105)
         mat(k,2031) = rxt(k,335)*y(k,210)
         mat(k,1305) = rxt(k,335)*y(k,204)
         mat(k,507) = -(rxt(k,350)*y(k,218))
         mat(k,1613) = -rxt(k,350)*y(k,106)
         mat(k,1866) = .800_r8*rxt(k,359)*y(k,190)
         mat(k,883) = .800_r8*rxt(k,359)*y(k,123)
         mat(k,316) = -(rxt(k,351)*y(k,218))
         mat(k,1586) = -rxt(k,351)*y(k,107)
         mat(k,2032) = .800_r8*rxt(k,348)*y(k,214)
         mat(k,677) = .800_r8*rxt(k,348)*y(k,204)
         mat(k,573) = -(rxt(k,352)*y(k,218))
         mat(k,1621) = -rxt(k,352)*y(k,108)
         mat(k,1941) = rxt(k,355)*y(k,212)
         mat(k,1348) = rxt(k,355)*y(k,124)
         mat(k,949) = -(rxt(k,443)*y(k,125) + rxt(k,444)*y(k,134) + rxt(k,445) &
                      *y(k,218))
         mat(k,1743) = -rxt(k,443)*y(k,109)
         mat(k,2219) = -rxt(k,444)*y(k,109)
         mat(k,1659) = -rxt(k,445)*y(k,109)
         mat(k,1333) = -(rxt(k,353)*y(k,134) + rxt(k,354)*y(k,218))
         mat(k,2239) = -rxt(k,353)*y(k,110)
         mat(k,1682) = -rxt(k,354)*y(k,110)
         mat(k,843) = .200_r8*rxt(k,386)*y(k,134)
         mat(k,1910) = .560_r8*rxt(k,369)*y(k,206)
         mat(k,1766) = .600_r8*rxt(k,370)*y(k,206)
         mat(k,2239) = mat(k,2239) + .200_r8*rxt(k,386)*y(k,97)
         mat(k,1385) = .610_r8*rxt(k,366)*y(k,206)
         mat(k,1819) = .440_r8*rxt(k,367)*y(k,206)
         mat(k,1259) = .560_r8*rxt(k,369)*y(k,123) + .600_r8*rxt(k,370)*y(k,125) &
                      + .610_r8*rxt(k,366)*y(k,198) + .440_r8*rxt(k,367)*y(k,199)
         mat(k,492) = -(rxt(k,152)*y(k,123) + (rxt(k,153) + rxt(k,154) + rxt(k,155) &
                      ) * y(k,124) + rxt(k,164)*y(k,218))
         mat(k,1864) = -rxt(k,152)*y(k,111)
         mat(k,1937) = -(rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,111)
         mat(k,1611) = -rxt(k,164)*y(k,111)
         mat(k,184) = -((rxt(k,168) + rxt(k,169)) * y(k,217))
         mat(k,1510) = -(rxt(k,168) + rxt(k,169)) * y(k,112)
         mat(k,491) = rxt(k,153)*y(k,124)
         mat(k,1933) = rxt(k,153)*y(k,111)
         mat(k,1934) = rxt(k,171)*y(k,125)
         mat(k,1733) = rxt(k,171)*y(k,124)
         mat(k,374) = -(rxt(k,389)*y(k,218))
         mat(k,1594) = -rxt(k,389)*y(k,114)
         mat(k,1185) = .200_r8*rxt(k,381)*y(k,199)
         mat(k,1790) = .200_r8*rxt(k,381)*y(k,100)
         mat(k,1074) = -(rxt(k,390)*y(k,218))
         mat(k,1666) = -rxt(k,390)*y(k,115)
         mat(k,1190) = rxt(k,383)*y(k,123) + rxt(k,384)*y(k,125) + rxt(k,380)*y(k,198) &
                      + .800_r8*rxt(k,381)*y(k,199)
         mat(k,1895) = rxt(k,383)*y(k,100)
         mat(k,1750) = rxt(k,384)*y(k,100)
         mat(k,1374) = rxt(k,380)*y(k,100)
         mat(k,1805) = .800_r8*rxt(k,381)*y(k,100)
         mat(k,97) = -(rxt(k,480)*y(k,218))
         mat(k,1555) = -rxt(k,480)*y(k,119)
         mat(k,1923) = -(rxt(k,152)*y(k,111) + rxt(k,161)*y(k,125) + rxt(k,165) &
                      *y(k,204) + rxt(k,166)*y(k,134) + rxt(k,167)*y(k,133) + rxt(k,188) &
                      *y(k,58) + rxt(k,220)*y(k,18) + rxt(k,263)*y(k,199) + rxt(k,271) &
                      *y(k,205) + rxt(k,284)*y(k,195) + rxt(k,295)*y(k,198) + rxt(k,299) &
                      *y(k,203) + rxt(k,312)*y(k,196) + rxt(k,320)*y(k,220) + rxt(k,324) &
                      *y(k,221) + (rxt(k,330) + rxt(k,331)) * y(k,201) + (rxt(k,337) &
                      + rxt(k,338)) * y(k,210) + rxt(k,346)*y(k,212) + rxt(k,349) &
                      *y(k,214) + (rxt(k,359) + rxt(k,360)) * y(k,190) + rxt(k,369) &
                      *y(k,206) + rxt(k,375)*y(k,207) + rxt(k,383)*y(k,100) + rxt(k,394) &
                      *y(k,226) + rxt(k,398)*y(k,189) + rxt(k,401)*y(k,192) + rxt(k,406) &
                      *y(k,194) + rxt(k,408)*y(k,197) + rxt(k,412)*y(k,200) + rxt(k,415) &
                      *y(k,211) + rxt(k,418)*y(k,213) + rxt(k,421)*y(k,219) + rxt(k,428) &
                      *y(k,224) + rxt(k,434)*y(k,227) + rxt(k,437)*y(k,229) + rxt(k,448) &
                      *y(k,216) + rxt(k,453)*y(k,222) + rxt(k,458)*y(k,223))
         mat(k,496) = -rxt(k,152)*y(k,123)
         mat(k,1779) = -rxt(k,161)*y(k,123)
         mat(k,2113) = -rxt(k,165)*y(k,123)
         mat(k,2252) = -rxt(k,166)*y(k,123)
         mat(k,2191) = -rxt(k,167)*y(k,123)
         mat(k,1722) = -rxt(k,188)*y(k,123)
         mat(k,2160) = -rxt(k,220)*y(k,123)
         mat(k,1831) = -rxt(k,263)*y(k,123)
         mat(k,445) = -rxt(k,271)*y(k,123)
         mat(k,829) = -rxt(k,284)*y(k,123)
         mat(k,1394) = -rxt(k,295)*y(k,123)
         mat(k,711) = -rxt(k,299)*y(k,123)
         mat(k,794) = -rxt(k,312)*y(k,123)
         mat(k,774) = -rxt(k,320)*y(k,123)
         mat(k,1137) = -rxt(k,324)*y(k,123)
         mat(k,570) = -(rxt(k,330) + rxt(k,331)) * y(k,123)
         mat(k,1320) = -(rxt(k,337) + rxt(k,338)) * y(k,123)
         mat(k,1362) = -rxt(k,346)*y(k,123)
         mat(k,682) = -rxt(k,349)*y(k,123)
         mat(k,895) = -(rxt(k,359) + rxt(k,360)) * y(k,123)
         mat(k,1267) = -rxt(k,369)*y(k,123)
         mat(k,1299) = -rxt(k,375)*y(k,123)
         mat(k,1204) = -rxt(k,383)*y(k,123)
         mat(k,1181) = -rxt(k,394)*y(k,123)
         mat(k,522) = -rxt(k,398)*y(k,123)
         mat(k,488) = -rxt(k,401)*y(k,123)
         mat(k,439) = -rxt(k,406)*y(k,123)
         mat(k,628) = -rxt(k,408)*y(k,123)
         mat(k,765) = -rxt(k,412)*y(k,123)
         mat(k,717) = -rxt(k,415)*y(k,123)
         mat(k,880) = -rxt(k,418)*y(k,123)
         mat(k,452) = -rxt(k,421)*y(k,123)
         mat(k,732) = -rxt(k,428)*y(k,123)
         mat(k,757) = -rxt(k,434)*y(k,123)
         mat(k,504) = -rxt(k,437)*y(k,123)
         mat(k,1055) = -rxt(k,448)*y(k,123)
         mat(k,1118) = -rxt(k,453)*y(k,123)
         mat(k,920) = -rxt(k,458)*y(k,123)
         mat(k,496) = mat(k,496) + 2.000_r8*rxt(k,154)*y(k,124) + rxt(k,164)*y(k,218)
         mat(k,186) = 2.000_r8*rxt(k,168)*y(k,217)
         mat(k,1967) = 2.000_r8*rxt(k,154)*y(k,111) + rxt(k,157)*y(k,133) + rxt(k,473) &
                      *y(k,150)
         mat(k,2191) = mat(k,2191) + rxt(k,157)*y(k,124)
         mat(k,1238) = rxt(k,473)*y(k,124)
         mat(k,1532) = 2.000_r8*rxt(k,168)*y(k,112)
         mat(k,1696) = rxt(k,164)*y(k,111)
         mat(k,1968) = -((rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,111) + (rxt(k,157) &
                      + rxt(k,159)) * y(k,133) + rxt(k,158)*y(k,134) + rxt(k,170) &
                      *y(k,204) + rxt(k,171)*y(k,125) + rxt(k,172)*y(k,218) + rxt(k,190) &
                      *y(k,58) + rxt(k,221)*y(k,18) + rxt(k,306)*y(k,198) + rxt(k,355) &
                      *y(k,212) + rxt(k,413)*y(k,200) + rxt(k,416)*y(k,211) + rxt(k,419) &
                      *y(k,213) + rxt(k,423)*y(k,141) + rxt(k,426)*y(k,189) + rxt(k,473) &
                      *y(k,150))
         mat(k,497) = -(rxt(k,153) + rxt(k,154) + rxt(k,155)) * y(k,124)
         mat(k,2192) = -(rxt(k,157) + rxt(k,159)) * y(k,124)
         mat(k,2253) = -rxt(k,158)*y(k,124)
         mat(k,2114) = -rxt(k,170)*y(k,124)
         mat(k,1780) = -rxt(k,171)*y(k,124)
         mat(k,1697) = -rxt(k,172)*y(k,124)
         mat(k,1723) = -rxt(k,190)*y(k,124)
         mat(k,2161) = -rxt(k,221)*y(k,124)
         mat(k,1395) = -rxt(k,306)*y(k,124)
         mat(k,1363) = -rxt(k,355)*y(k,124)
         mat(k,766) = -rxt(k,413)*y(k,124)
         mat(k,718) = -rxt(k,416)*y(k,124)
         mat(k,881) = -rxt(k,419)*y(k,124)
         mat(k,466) = -rxt(k,423)*y(k,124)
         mat(k,523) = -rxt(k,426)*y(k,124)
         mat(k,1239) = -rxt(k,473)*y(k,124)
         mat(k,675) = rxt(k,357)*y(k,218)
         mat(k,357) = rxt(k,328)*y(k,125)
         mat(k,2161) = mat(k,2161) + rxt(k,220)*y(k,123)
         mat(k,1723) = mat(k,1723) + rxt(k,188)*y(k,123)
         mat(k,425) = rxt(k,151)*y(k,218)
         mat(k,589) = .700_r8*rxt(k,377)*y(k,218)
         mat(k,1205) = rxt(k,383)*y(k,123) + rxt(k,384)*y(k,125)
         mat(k,1924) = rxt(k,220)*y(k,18) + rxt(k,188)*y(k,58) + rxt(k,383)*y(k,100) &
                      + 2.000_r8*rxt(k,161)*y(k,125) + rxt(k,167)*y(k,133) &
                      + rxt(k,166)*y(k,134) + rxt(k,398)*y(k,189) + rxt(k,359) &
                      *y(k,190) + rxt(k,401)*y(k,192) + rxt(k,406)*y(k,194) &
                      + rxt(k,284)*y(k,195) + rxt(k,312)*y(k,196) + rxt(k,408) &
                      *y(k,197) + rxt(k,295)*y(k,198) + rxt(k,263)*y(k,199) &
                      + rxt(k,412)*y(k,200) + rxt(k,330)*y(k,201) + rxt(k,299) &
                      *y(k,203) + rxt(k,165)*y(k,204) + rxt(k,271)*y(k,205) &
                      + .920_r8*rxt(k,369)*y(k,206) + .920_r8*rxt(k,375)*y(k,207) &
                      + rxt(k,337)*y(k,210) + rxt(k,415)*y(k,211) + rxt(k,346) &
                      *y(k,212) + rxt(k,418)*y(k,213) + rxt(k,349)*y(k,214) &
                      + 1.600_r8*rxt(k,448)*y(k,216) + rxt(k,421)*y(k,219) &
                      + rxt(k,320)*y(k,220) + rxt(k,324)*y(k,221) + .900_r8*rxt(k,453) &
                      *y(k,222) + .800_r8*rxt(k,458)*y(k,223) + rxt(k,428)*y(k,224) &
                      + rxt(k,394)*y(k,226) + rxt(k,434)*y(k,227) + rxt(k,437) &
                      *y(k,229)
         mat(k,1780) = mat(k,1780) + rxt(k,328)*y(k,15) + rxt(k,384)*y(k,100) &
                      + 2.000_r8*rxt(k,161)*y(k,123) + rxt(k,162)*y(k,133) &
                      + rxt(k,160)*y(k,204) + rxt(k,370)*y(k,206) + rxt(k,376) &
                      *y(k,207) + rxt(k,336)*y(k,210) + rxt(k,347)*y(k,212) &
                      + 2.000_r8*rxt(k,449)*y(k,216) + rxt(k,163)*y(k,218) &
                      + rxt(k,395)*y(k,226)
         mat(k,863) = rxt(k,318)*y(k,218)
         mat(k,2192) = mat(k,2192) + rxt(k,167)*y(k,123) + rxt(k,162)*y(k,125)
         mat(k,2253) = mat(k,2253) + rxt(k,166)*y(k,123)
         mat(k,622) = rxt(k,455)*y(k,218)
         mat(k,523) = mat(k,523) + rxt(k,398)*y(k,123)
         mat(k,896) = rxt(k,359)*y(k,123)
         mat(k,489) = rxt(k,401)*y(k,123)
         mat(k,440) = rxt(k,406)*y(k,123)
         mat(k,830) = rxt(k,284)*y(k,123)
         mat(k,795) = rxt(k,312)*y(k,123)
         mat(k,629) = rxt(k,408)*y(k,123)
         mat(k,1395) = mat(k,1395) + rxt(k,295)*y(k,123)
         mat(k,1832) = rxt(k,263)*y(k,123) + .500_r8*rxt(k,446)*y(k,216)
         mat(k,766) = mat(k,766) + rxt(k,412)*y(k,123)
         mat(k,571) = rxt(k,330)*y(k,123)
         mat(k,712) = rxt(k,299)*y(k,123)
         mat(k,2114) = mat(k,2114) + rxt(k,165)*y(k,123) + rxt(k,160)*y(k,125)
         mat(k,446) = rxt(k,271)*y(k,123)
         mat(k,1268) = .920_r8*rxt(k,369)*y(k,123) + rxt(k,370)*y(k,125)
         mat(k,1300) = .920_r8*rxt(k,375)*y(k,123) + rxt(k,376)*y(k,125)
         mat(k,1321) = rxt(k,337)*y(k,123) + rxt(k,336)*y(k,125)
         mat(k,718) = mat(k,718) + rxt(k,415)*y(k,123)
         mat(k,1363) = mat(k,1363) + rxt(k,346)*y(k,123) + rxt(k,347)*y(k,125)
         mat(k,881) = mat(k,881) + rxt(k,418)*y(k,123)
         mat(k,683) = rxt(k,349)*y(k,123)
         mat(k,1056) = 1.600_r8*rxt(k,448)*y(k,123) + 2.000_r8*rxt(k,449)*y(k,125) &
                      + .500_r8*rxt(k,446)*y(k,199)
         mat(k,1697) = mat(k,1697) + rxt(k,357)*y(k,1) + rxt(k,151)*y(k,89) &
                      + .700_r8*rxt(k,377)*y(k,98) + rxt(k,163)*y(k,125) + rxt(k,318) &
                      *y(k,126) + rxt(k,455)*y(k,176)
         mat(k,453) = rxt(k,421)*y(k,123)
         mat(k,775) = rxt(k,320)*y(k,123)
         mat(k,1138) = rxt(k,324)*y(k,123)
         mat(k,1119) = .900_r8*rxt(k,453)*y(k,123)
         mat(k,921) = .800_r8*rxt(k,458)*y(k,123)
         mat(k,733) = rxt(k,428)*y(k,123)
         mat(k,1182) = rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125)
         mat(k,758) = rxt(k,434)*y(k,123)
         mat(k,505) = rxt(k,437)*y(k,123)
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
         mat(k,1777) = -(rxt(k,160)*y(k,204) + rxt(k,161)*y(k,123) + rxt(k,162) &
                      *y(k,133) + rxt(k,163)*y(k,218) + rxt(k,171)*y(k,124) + rxt(k,257) &
                      *y(k,41) + rxt(k,289)*y(k,44) + rxt(k,308)*y(k,28) + rxt(k,315) &
                      *y(k,48) + rxt(k,328)*y(k,15) + rxt(k,336)*y(k,210) + rxt(k,347) &
                      *y(k,212) + rxt(k,370)*y(k,206) + rxt(k,376)*y(k,207) + rxt(k,379) &
                      *y(k,97) + rxt(k,384)*y(k,100) + rxt(k,395)*y(k,226) + rxt(k,440) &
                      *y(k,5) + rxt(k,443)*y(k,109) + rxt(k,449)*y(k,216) + rxt(k,460) &
                      *y(k,178) + rxt(k,463)*y(k,66))
         mat(k,2111) = -rxt(k,160)*y(k,125)
         mat(k,1921) = -rxt(k,161)*y(k,125)
         mat(k,2189) = -rxt(k,162)*y(k,125)
         mat(k,1694) = -rxt(k,163)*y(k,125)
         mat(k,1965) = -rxt(k,171)*y(k,125)
         mat(k,1489) = -rxt(k,257)*y(k,125)
         mat(k,1091) = -rxt(k,289)*y(k,125)
         mat(k,1034) = -rxt(k,308)*y(k,125)
         mat(k,1226) = -rxt(k,315)*y(k,125)
         mat(k,356) = -rxt(k,328)*y(k,125)
         mat(k,1318) = -rxt(k,336)*y(k,125)
         mat(k,1360) = -rxt(k,347)*y(k,125)
         mat(k,1265) = -rxt(k,370)*y(k,125)
         mat(k,1297) = -rxt(k,376)*y(k,125)
         mat(k,847) = -rxt(k,379)*y(k,125)
         mat(k,1202) = -rxt(k,384)*y(k,125)
         mat(k,1179) = -rxt(k,395)*y(k,125)
         mat(k,1012) = -rxt(k,440)*y(k,125)
         mat(k,962) = -rxt(k,443)*y(k,125)
         mat(k,1053) = -rxt(k,449)*y(k,125)
         mat(k,976) = -rxt(k,460)*y(k,125)
         mat(k,287) = -rxt(k,463)*y(k,125)
         mat(k,561) = rxt(k,222)*y(k,133)
         mat(k,2004) = rxt(k,189)*y(k,59)
         mat(k,904) = rxt(k,189)*y(k,55) + rxt(k,191)*y(k,133) + rxt(k,192)*y(k,218)
         mat(k,871) = rxt(k,236)*y(k,88)
         mat(k,1468) = rxt(k,236)*y(k,72) + rxt(k,173)*y(k,218)
         mat(k,579) = .500_r8*rxt(k,352)*y(k,218)
         mat(k,1965) = mat(k,1965) + rxt(k,159)*y(k,133) + rxt(k,158)*y(k,134)
         mat(k,2189) = mat(k,2189) + rxt(k,222)*y(k,19) + rxt(k,191)*y(k,59) &
                      + rxt(k,159)*y(k,124)
         mat(k,2250) = rxt(k,158)*y(k,124)
         mat(k,530) = rxt(k,304)*y(k,218)
         mat(k,1694) = mat(k,1694) + rxt(k,192)*y(k,59) + rxt(k,173)*y(k,88) &
                      + .500_r8*rxt(k,352)*y(k,108) + rxt(k,304)*y(k,139)
         mat(k,858) = -(rxt(k,318)*y(k,218))
         mat(k,1651) = -rxt(k,318)*y(k,126)
         mat(k,1023) = rxt(k,308)*y(k,125)
         mat(k,534) = .500_r8*rxt(k,378)*y(k,218)
         mat(k,406) = rxt(k,385)*y(k,218)
         mat(k,375) = rxt(k,389)*y(k,218)
         mat(k,1071) = rxt(k,390)*y(k,218)
         mat(k,1740) = rxt(k,308)*y(k,28)
         mat(k,1651) = mat(k,1651) + .500_r8*rxt(k,378)*y(k,99) + rxt(k,385)*y(k,101) &
                      + rxt(k,389)*y(k,114) + rxt(k,390)*y(k,115)
         mat(k,392) = -(rxt(k,450)*y(k,218))
         mat(k,1596) = -rxt(k,450)*y(k,127)
         mat(k,2036) = rxt(k,447)*y(k,216)
         mat(k,1042) = rxt(k,447)*y(k,204)
         mat(k,2197) = -(rxt(k,131)*y(k,134) + 4._r8*rxt(k,132)*y(k,133) + rxt(k,134) &
                      *y(k,76) + rxt(k,135)*y(k,78) + rxt(k,140)*y(k,204) + rxt(k,146) &
                      *y(k,218) + (rxt(k,157) + rxt(k,159)) * y(k,124) + rxt(k,162) &
                      *y(k,125) + rxt(k,167)*y(k,123) + rxt(k,191)*y(k,59) + rxt(k,193) &
                      *y(k,58) + rxt(k,196)*y(k,84) + rxt(k,199)*y(k,91) + rxt(k,222) &
                      *y(k,19) + rxt(k,223)*y(k,18) + rxt(k,225)*y(k,80) + rxt(k,227) &
                      *y(k,90) + rxt(k,258)*y(k,41) + rxt(k,465)*y(k,137))
         mat(k,2258) = -rxt(k,131)*y(k,133)
         mat(k,1412) = -rxt(k,134)*y(k,133)
         mat(k,609) = -rxt(k,135)*y(k,133)
         mat(k,2119) = -rxt(k,140)*y(k,133)
         mat(k,1702) = -rxt(k,146)*y(k,133)
         mat(k,1973) = -(rxt(k,157) + rxt(k,159)) * y(k,133)
         mat(k,1785) = -rxt(k,162)*y(k,133)
         mat(k,1929) = -rxt(k,167)*y(k,133)
         mat(k,908) = -rxt(k,191)*y(k,133)
         mat(k,1728) = -rxt(k,193)*y(k,133)
         mat(k,2142) = -rxt(k,196)*y(k,133)
         mat(k,819) = -rxt(k,199)*y(k,133)
         mat(k,564) = -rxt(k,222)*y(k,133)
         mat(k,2166) = -rxt(k,223)*y(k,133)
         mat(k,811) = -rxt(k,225)*y(k,133)
         mat(k,785) = -rxt(k,227)*y(k,133)
         mat(k,1497) = -rxt(k,258)*y(k,133)
         mat(k,365) = -rxt(k,465)*y(k,133)
         mat(k,1456) = rxt(k,138)*y(k,204)
         mat(k,498) = rxt(k,152)*y(k,123) + rxt(k,153)*y(k,124)
         mat(k,1929) = mat(k,1929) + rxt(k,152)*y(k,111)
         mat(k,1973) = mat(k,1973) + rxt(k,153)*y(k,111)
         mat(k,2258) = mat(k,2258) + 2.000_r8*rxt(k,130)*y(k,217)
         mat(k,2119) = mat(k,2119) + rxt(k,138)*y(k,75)
         mat(k,1538) = 2.000_r8*rxt(k,130)*y(k,134)
         mat(k,1702) = mat(k,1702) + 2.000_r8*rxt(k,148)*y(k,218)
         mat(k,2259) = -((rxt(k,129) + rxt(k,130)) * y(k,217) + rxt(k,131)*y(k,133) &
                      + rxt(k,141)*y(k,204) + rxt(k,142)*y(k,75) + rxt(k,147)*y(k,218) &
                      + rxt(k,158)*y(k,124) + rxt(k,166)*y(k,123) + rxt(k,182)*y(k,55) &
                      + rxt(k,214)*y(k,16) + rxt(k,280)*y(k,24) + rxt(k,309)*y(k,28) &
                      + rxt(k,339)*y(k,104) + rxt(k,353)*y(k,110) + rxt(k,386)*y(k,97) &
                      + rxt(k,424)*y(k,141) + rxt(k,441)*y(k,5) + rxt(k,444)*y(k,109) &
                      + rxt(k,469)*y(k,148) + rxt(k,475)*y(k,150))
         mat(k,1539) = -(rxt(k,129) + rxt(k,130)) * y(k,134)
         mat(k,2198) = -rxt(k,131)*y(k,134)
         mat(k,2120) = -rxt(k,141)*y(k,134)
         mat(k,1457) = -rxt(k,142)*y(k,134)
         mat(k,1703) = -rxt(k,147)*y(k,134)
         mat(k,1974) = -rxt(k,158)*y(k,134)
         mat(k,1930) = -rxt(k,166)*y(k,134)
         mat(k,2013) = -rxt(k,182)*y(k,134)
         mat(k,1424) = -rxt(k,214)*y(k,134)
         mat(k,556) = -rxt(k,280)*y(k,134)
         mat(k,1040) = -rxt(k,309)*y(k,134)
         mat(k,1219) = -rxt(k,339)*y(k,134)
         mat(k,1346) = -rxt(k,353)*y(k,134)
         mat(k,850) = -rxt(k,386)*y(k,134)
         mat(k,467) = -rxt(k,424)*y(k,134)
         mat(k,1018) = -rxt(k,441)*y(k,134)
         mat(k,968) = -rxt(k,444)*y(k,134)
         mat(k,517) = -rxt(k,469)*y(k,134)
         mat(k,1244) = -rxt(k,475)*y(k,134)
         mat(k,1398) = .150_r8*rxt(k,294)*y(k,204)
         mat(k,2120) = mat(k,2120) + .150_r8*rxt(k,294)*y(k,198) + .150_r8*rxt(k,344) &
                      *y(k,212)
         mat(k,1366) = .150_r8*rxt(k,344)*y(k,204)
         mat(k,326) = -(rxt(k,476)*y(k,150))
         mat(k,1230) = -rxt(k,476)*y(k,136)
         mat(k,2146) = rxt(k,216)*y(k,58)
         mat(k,1708) = rxt(k,216)*y(k,18) + 2.000_r8*rxt(k,186)*y(k,58)
         mat(k,358) = -(rxt(k,465)*y(k,133) + rxt(k,466)*y(k,218))
         mat(k,2169) = -rxt(k,465)*y(k,137)
         mat(k,1592) = -rxt(k,466)*y(k,137)
         mat(k,1147) = rxt(k,332)*y(k,218)
         mat(k,1852) = .100_r8*rxt(k,453)*y(k,222)
         mat(k,1572) = rxt(k,332)*y(k,92)
         mat(k,1103) = .100_r8*rxt(k,453)*y(k,123)
         mat(k,525) = -(rxt(k,304)*y(k,218))
         mat(k,1616) = -rxt(k,304)*y(k,139)
         mat(k,1939) = rxt(k,306)*y(k,198)
         mat(k,1369) = rxt(k,306)*y(k,124)
         mat(k,1932) = rxt(k,426)*y(k,189)
         mat(k,518) = rxt(k,426)*y(k,124)
         mat(k,464) = -(rxt(k,423)*y(k,124) + rxt(k,424)*y(k,134))
         mat(k,1936) = -rxt(k,423)*y(k,141)
         mat(k,2207) = -rxt(k,424)*y(k,141)
         mat(k,197) = .070_r8*rxt(k,410)*y(k,218)
         mat(k,1862) = rxt(k,408)*y(k,197)
         mat(k,175) = .060_r8*rxt(k,422)*y(k,218)
         mat(k,218) = .070_r8*rxt(k,438)*y(k,218)
         mat(k,625) = rxt(k,408)*y(k,123)
         mat(k,1607) = .070_r8*rxt(k,410)*y(k,65) + .060_r8*rxt(k,422)*y(k,142) &
                      + .070_r8*rxt(k,438)*y(k,185)
         mat(k,173) = -(rxt(k,422)*y(k,218))
         mat(k,1562) = -rxt(k,422)*y(k,142)
         mat(k,165) = .530_r8*rxt(k,399)*y(k,218)
         mat(k,1562) = mat(k,1562) + .530_r8*rxt(k,399)*y(k,6)
         mat(k,331) = -(rxt(k,425)*y(k,218))
         mat(k,1587) = -rxt(k,425)*y(k,143)
         mat(k,2033) = rxt(k,420)*y(k,219)
         mat(k,448) = rxt(k,420)*y(k,204)
         mat(k,541) = -(rxt(k,321)*y(k,218))
         mat(k,1618) = -rxt(k,321)*y(k,146)
         mat(k,2053) = rxt(k,319)*y(k,220)
         mat(k,768) = rxt(k,319)*y(k,204)
         mat(k,398) = -(rxt(k,325)*y(k,218))
         mat(k,1597) = -rxt(k,325)*y(k,147)
         mat(k,2037) = .850_r8*rxt(k,323)*y(k,221)
         mat(k,1128) = .850_r8*rxt(k,323)*y(k,204)
         mat(k,512) = -(rxt(k,469)*y(k,134) + rxt(k,472)*y(k,218))
         mat(k,2208) = -rxt(k,469)*y(k,148)
         mat(k,1614) = -rxt(k,472)*y(k,148)
         mat(k,1233) = -(rxt(k,470)*y(k,18) + rxt(k,471)*y(k,58) + rxt(k,473)*y(k,124) &
                      + rxt(k,475)*y(k,134) + rxt(k,476)*y(k,136) + rxt(k,477) &
                      *y(k,218))
         mat(k,2150) = -rxt(k,470)*y(k,150)
         mat(k,1712) = -rxt(k,471)*y(k,150)
         mat(k,1954) = -rxt(k,473)*y(k,150)
         mat(k,2235) = -rxt(k,475)*y(k,150)
         mat(k,328) = -rxt(k,476)*y(k,150)
         mat(k,1678) = -rxt(k,477)*y(k,150)
         mat(k,2179) = rxt(k,465)*y(k,137)
         mat(k,2235) = mat(k,2235) + rxt(k,469)*y(k,148)
         mat(k,362) = rxt(k,465)*y(k,133)
         mat(k,513) = rxt(k,469)*y(k,134) + rxt(k,472)*y(k,218)
         mat(k,1678) = mat(k,1678) + rxt(k,472)*y(k,148)
         mat(k,852) = -(rxt(k,468)*y(k,218))
         mat(k,1650) = -rxt(k,468)*y(k,151)
         mat(k,2149) = rxt(k,470)*y(k,150)
         mat(k,1710) = rxt(k,471)*y(k,150)
         mat(k,284) = rxt(k,463)*y(k,125) + (rxt(k,464)+.500_r8*rxt(k,478))*y(k,218)
         mat(k,1947) = rxt(k,473)*y(k,150)
         mat(k,1739) = rxt(k,463)*y(k,66)
         mat(k,2215) = rxt(k,475)*y(k,150)
         mat(k,327) = rxt(k,476)*y(k,150)
         mat(k,360) = rxt(k,466)*y(k,218)
         mat(k,1232) = rxt(k,470)*y(k,18) + rxt(k,471)*y(k,58) + rxt(k,473)*y(k,124) &
                      + rxt(k,475)*y(k,134) + rxt(k,476)*y(k,136) + rxt(k,477) &
                      *y(k,218)
         mat(k,1650) = mat(k,1650) + (rxt(k,464)+.500_r8*rxt(k,478))*y(k,66) &
                      + rxt(k,466)*y(k,137) + rxt(k,477)*y(k,150)
         mat(k,258) = -(rxt(k,479)*y(k,230))
         mat(k,2262) = -rxt(k,479)*y(k,152)
         mat(k,851) = rxt(k,468)*y(k,218)
         mat(k,1578) = rxt(k,468)*y(k,151)
         mat(k,985) = .2202005_r8*rxt(k,497)*y(k,134)
         mat(k,936) = .0508005_r8*rxt(k,513)*y(k,134)
         mat(k,1840) = .1279005_r8*rxt(k,496)*y(k,191) + .0097005_r8*rxt(k,501) &
                      *y(k,193) + .0003005_r8*rxt(k,504)*y(k,208) &
                      + .1056005_r8*rxt(k,508)*y(k,209) + .0245005_r8*rxt(k,512) &
                      *y(k,215) + .0154005_r8*rxt(k,518)*y(k,225) &
                      + .0063005_r8*rxt(k,522)*y(k,228)
         mat(k,2200) = .2202005_r8*rxt(k,497)*y(k,5) + .0508005_r8*rxt(k,513)*y(k,109)
         mat(k,44) = .5931005_r8*rxt(k,515)*y(k,218)
         mat(k,50) = .1279005_r8*rxt(k,496)*y(k,123) + .2202005_r8*rxt(k,495)*y(k,204)
         mat(k,56) = .0097005_r8*rxt(k,501)*y(k,123) + .0023005_r8*rxt(k,500)*y(k,204)
         mat(k,2015) = .2202005_r8*rxt(k,495)*y(k,191) + .0023005_r8*rxt(k,500) &
                      *y(k,193) + .0031005_r8*rxt(k,503)*y(k,208) &
                      + .2381005_r8*rxt(k,507)*y(k,209) + .0508005_r8*rxt(k,511) &
                      *y(k,215) + .1364005_r8*rxt(k,517)*y(k,225) &
                      + .1677005_r8*rxt(k,521)*y(k,228)
         mat(k,62) = .0003005_r8*rxt(k,504)*y(k,123) + .0031005_r8*rxt(k,503)*y(k,204)
         mat(k,68) = .1056005_r8*rxt(k,508)*y(k,123) + .2381005_r8*rxt(k,507)*y(k,204)
         mat(k,76) = .0245005_r8*rxt(k,512)*y(k,123) + .0508005_r8*rxt(k,511)*y(k,204)
         mat(k,1541) = .5931005_r8*rxt(k,515)*y(k,173)
         mat(k,82) = .0154005_r8*rxt(k,518)*y(k,123) + .1364005_r8*rxt(k,517)*y(k,204)
         mat(k,88) = .0063005_r8*rxt(k,522)*y(k,123) + .1677005_r8*rxt(k,521)*y(k,204)
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
         mat(k,986) = .2067005_r8*rxt(k,497)*y(k,134)
         mat(k,937) = .1149005_r8*rxt(k,513)*y(k,134)
         mat(k,1841) = .1792005_r8*rxt(k,496)*y(k,191) + .0034005_r8*rxt(k,501) &
                      *y(k,193) + .0003005_r8*rxt(k,504)*y(k,208) &
                      + .1026005_r8*rxt(k,508)*y(k,209) + .0082005_r8*rxt(k,512) &
                      *y(k,215) + .0452005_r8*rxt(k,518)*y(k,225) &
                      + .0237005_r8*rxt(k,522)*y(k,228)
         mat(k,2201) = .2067005_r8*rxt(k,497)*y(k,5) + .1149005_r8*rxt(k,513)*y(k,109)
         mat(k,45) = .1534005_r8*rxt(k,515)*y(k,218)
         mat(k,51) = .1792005_r8*rxt(k,496)*y(k,123) + .2067005_r8*rxt(k,495)*y(k,204)
         mat(k,57) = .0034005_r8*rxt(k,501)*y(k,123) + .0008005_r8*rxt(k,500)*y(k,204)
         mat(k,2016) = .2067005_r8*rxt(k,495)*y(k,191) + .0008005_r8*rxt(k,500) &
                      *y(k,193) + .0035005_r8*rxt(k,503)*y(k,208) &
                      + .1308005_r8*rxt(k,507)*y(k,209) + .1149005_r8*rxt(k,511) &
                      *y(k,215) + .0101005_r8*rxt(k,517)*y(k,225) &
                      + .0174005_r8*rxt(k,521)*y(k,228)
         mat(k,63) = .0003005_r8*rxt(k,504)*y(k,123) + .0035005_r8*rxt(k,503)*y(k,204)
         mat(k,69) = .1026005_r8*rxt(k,508)*y(k,123) + .1308005_r8*rxt(k,507)*y(k,204)
         mat(k,77) = .0082005_r8*rxt(k,512)*y(k,123) + .1149005_r8*rxt(k,511)*y(k,204)
         mat(k,1542) = .1534005_r8*rxt(k,515)*y(k,173)
         mat(k,83) = .0452005_r8*rxt(k,518)*y(k,123) + .0101005_r8*rxt(k,517)*y(k,204)
         mat(k,89) = .0237005_r8*rxt(k,522)*y(k,123) + .0174005_r8*rxt(k,521)*y(k,204)
         mat(k,987) = .0653005_r8*rxt(k,497)*y(k,134)
         mat(k,938) = .0348005_r8*rxt(k,513)*y(k,134)
         mat(k,1842) = .0676005_r8*rxt(k,496)*y(k,191) + .1579005_r8*rxt(k,501) &
                      *y(k,193) + .0073005_r8*rxt(k,504)*y(k,208) &
                      + .0521005_r8*rxt(k,508)*y(k,209) + .0772005_r8*rxt(k,512) &
                      *y(k,215) + .0966005_r8*rxt(k,518)*y(k,225) &
                      + .0025005_r8*rxt(k,522)*y(k,228)
         mat(k,2202) = .0653005_r8*rxt(k,497)*y(k,5) + .0348005_r8*rxt(k,513)*y(k,109)
         mat(k,46) = .0459005_r8*rxt(k,515)*y(k,218)
         mat(k,52) = .0676005_r8*rxt(k,496)*y(k,123) + .0653005_r8*rxt(k,495)*y(k,204)
         mat(k,58) = .1579005_r8*rxt(k,501)*y(k,123) + .0843005_r8*rxt(k,500)*y(k,204)
         mat(k,2017) = .0653005_r8*rxt(k,495)*y(k,191) + .0843005_r8*rxt(k,500) &
                      *y(k,193) + .0003005_r8*rxt(k,503)*y(k,208) &
                      + .0348005_r8*rxt(k,507)*y(k,209) + .0348005_r8*rxt(k,511) &
                      *y(k,215) + .0763005_r8*rxt(k,517)*y(k,225) + .086_r8*rxt(k,521) &
                      *y(k,228)
         mat(k,64) = .0073005_r8*rxt(k,504)*y(k,123) + .0003005_r8*rxt(k,503)*y(k,204)
         mat(k,70) = .0521005_r8*rxt(k,508)*y(k,123) + .0348005_r8*rxt(k,507)*y(k,204)
         mat(k,78) = .0772005_r8*rxt(k,512)*y(k,123) + .0348005_r8*rxt(k,511)*y(k,204)
         mat(k,1543) = .0459005_r8*rxt(k,515)*y(k,173)
         mat(k,84) = .0966005_r8*rxt(k,518)*y(k,123) + .0763005_r8*rxt(k,517)*y(k,204)
         mat(k,90) = .0025005_r8*rxt(k,522)*y(k,123) + .086_r8*rxt(k,521)*y(k,204)
         mat(k,988) = .1749305_r8*rxt(k,494)*y(k,125) + .1284005_r8*rxt(k,497) &
                      *y(k,134)
         mat(k,832) = .0590245_r8*rxt(k,502)*y(k,125) + .0033005_r8*rxt(k,505) &
                      *y(k,134)
         mat(k,939) = .1749305_r8*rxt(k,510)*y(k,125) + .0554005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1843) = .079_r8*rxt(k,496)*y(k,191) + .0059005_r8*rxt(k,501)*y(k,193) &
                      + .0057005_r8*rxt(k,504)*y(k,208) + .0143005_r8*rxt(k,508) &
                      *y(k,209) + .0332005_r8*rxt(k,512)*y(k,215) &
                      + .0073005_r8*rxt(k,518)*y(k,225) + .011_r8*rxt(k,522)*y(k,228)
         mat(k,1731) = .1749305_r8*rxt(k,494)*y(k,5) + .0590245_r8*rxt(k,502)*y(k,97) &
                      + .1749305_r8*rxt(k,510)*y(k,109)
         mat(k,2203) = .1284005_r8*rxt(k,497)*y(k,5) + .0033005_r8*rxt(k,505)*y(k,97) &
                      + .0554005_r8*rxt(k,513)*y(k,109)
         mat(k,47) = .0085005_r8*rxt(k,515)*y(k,218)
         mat(k,53) = .079_r8*rxt(k,496)*y(k,123) + .1284005_r8*rxt(k,495)*y(k,204)
         mat(k,59) = .0059005_r8*rxt(k,501)*y(k,123) + .0443005_r8*rxt(k,500)*y(k,204)
         mat(k,2018) = .1284005_r8*rxt(k,495)*y(k,191) + .0443005_r8*rxt(k,500) &
                      *y(k,193) + .0271005_r8*rxt(k,503)*y(k,208) &
                      + .0076005_r8*rxt(k,507)*y(k,209) + .0554005_r8*rxt(k,511) &
                      *y(k,215) + .2157005_r8*rxt(k,517)*y(k,225) &
                      + .0512005_r8*rxt(k,521)*y(k,228)
         mat(k,65) = .0057005_r8*rxt(k,504)*y(k,123) + .0271005_r8*rxt(k,503)*y(k,204)
         mat(k,71) = .0143005_r8*rxt(k,508)*y(k,123) + .0076005_r8*rxt(k,507)*y(k,204)
         mat(k,79) = .0332005_r8*rxt(k,512)*y(k,123) + .0554005_r8*rxt(k,511)*y(k,204)
         mat(k,1544) = .0085005_r8*rxt(k,515)*y(k,173)
         mat(k,85) = .0073005_r8*rxt(k,518)*y(k,123) + .2157005_r8*rxt(k,517)*y(k,204)
         mat(k,91) = .011_r8*rxt(k,522)*y(k,123) + .0512005_r8*rxt(k,521)*y(k,204)
         mat(k,989) = .5901905_r8*rxt(k,494)*y(k,125) + .114_r8*rxt(k,497)*y(k,134)
         mat(k,833) = .0250245_r8*rxt(k,502)*y(k,125)
         mat(k,940) = .5901905_r8*rxt(k,510)*y(k,125) + .1278005_r8*rxt(k,513) &
                      *y(k,134)
         mat(k,1844) = .1254005_r8*rxt(k,496)*y(k,191) + .0536005_r8*rxt(k,501) &
                      *y(k,193) + .0623005_r8*rxt(k,504)*y(k,208) &
                      + .0166005_r8*rxt(k,508)*y(k,209) + .130_r8*rxt(k,512)*y(k,215) &
                      + .238_r8*rxt(k,518)*y(k,225) + .1185005_r8*rxt(k,522)*y(k,228)
         mat(k,1732) = .5901905_r8*rxt(k,494)*y(k,5) + .0250245_r8*rxt(k,502)*y(k,97) &
                      + .5901905_r8*rxt(k,510)*y(k,109)
         mat(k,2204) = .114_r8*rxt(k,497)*y(k,5) + .1278005_r8*rxt(k,513)*y(k,109)
         mat(k,48) = .0128005_r8*rxt(k,515)*y(k,218)
         mat(k,54) = .1254005_r8*rxt(k,496)*y(k,123) + .114_r8*rxt(k,495)*y(k,204)
         mat(k,60) = .0536005_r8*rxt(k,501)*y(k,123) + .1621005_r8*rxt(k,500)*y(k,204)
         mat(k,2019) = .114_r8*rxt(k,495)*y(k,191) + .1621005_r8*rxt(k,500)*y(k,193) &
                      + .0474005_r8*rxt(k,503)*y(k,208) + .0113005_r8*rxt(k,507) &
                      *y(k,209) + .1278005_r8*rxt(k,511)*y(k,215) &
                      + .0738005_r8*rxt(k,517)*y(k,225) + .1598005_r8*rxt(k,521) &
                      *y(k,228)
         mat(k,66) = .0623005_r8*rxt(k,504)*y(k,123) + .0474005_r8*rxt(k,503)*y(k,204)
         mat(k,72) = .0166005_r8*rxt(k,508)*y(k,123) + .0113005_r8*rxt(k,507)*y(k,204)
         mat(k,80) = .130_r8*rxt(k,512)*y(k,123) + .1278005_r8*rxt(k,511)*y(k,204)
         mat(k,1545) = .0128005_r8*rxt(k,515)*y(k,173)
         mat(k,86) = .238_r8*rxt(k,518)*y(k,123) + .0738005_r8*rxt(k,517)*y(k,204)
         mat(k,92) = .1185005_r8*rxt(k,522)*y(k,123) + .1598005_r8*rxt(k,521)*y(k,204)
         mat(k,49) = -(rxt(k,515)*y(k,218))
         mat(k,1546) = -rxt(k,515)*y(k,173)
         mat(k,190) = .100_r8*rxt(k,430)*y(k,218)
         mat(k,208) = .230_r8*rxt(k,432)*y(k,218)
         mat(k,1566) = .100_r8*rxt(k,430)*y(k,181) + .230_r8*rxt(k,432)*y(k,183)
         mat(k,643) = -(rxt(k,454)*y(k,218))
         mat(k,1630) = -rxt(k,454)*y(k,175)
         mat(k,2058) = rxt(k,452)*y(k,222)
         mat(k,1104) = rxt(k,452)*y(k,204)
         mat(k,618) = -(rxt(k,455)*y(k,218))
         mat(k,1627) = -rxt(k,455)*y(k,176)
         mat(k,1872) = .200_r8*rxt(k,448)*y(k,216) + .200_r8*rxt(k,458)*y(k,223)
         mat(k,1792) = .500_r8*rxt(k,446)*y(k,216)
         mat(k,1043) = .200_r8*rxt(k,448)*y(k,123) + .500_r8*rxt(k,446)*y(k,199)
         mat(k,911) = .200_r8*rxt(k,458)*y(k,123)
         mat(k,475) = -(rxt(k,459)*y(k,218))
         mat(k,1609) = -rxt(k,459)*y(k,177)
         mat(k,2049) = rxt(k,457)*y(k,223)
         mat(k,910) = rxt(k,457)*y(k,204)
         mat(k,970) = -(rxt(k,460)*y(k,125) + rxt(k,461)*y(k,218))
         mat(k,1744) = -rxt(k,460)*y(k,178)
         mat(k,1660) = -rxt(k,461)*y(k,178)
         mat(k,998) = .330_r8*rxt(k,441)*y(k,134)
         mat(k,950) = .330_r8*rxt(k,444)*y(k,134)
         mat(k,1891) = .800_r8*rxt(k,448)*y(k,216) + .800_r8*rxt(k,458)*y(k,223)
         mat(k,1744) = mat(k,1744) + rxt(k,449)*y(k,216)
         mat(k,2220) = .330_r8*rxt(k,441)*y(k,5) + .330_r8*rxt(k,444)*y(k,109)
         mat(k,619) = rxt(k,455)*y(k,218)
         mat(k,1801) = .500_r8*rxt(k,446)*y(k,216) + rxt(k,456)*y(k,223)
         mat(k,1045) = .800_r8*rxt(k,448)*y(k,123) + rxt(k,449)*y(k,125) &
                      + .500_r8*rxt(k,446)*y(k,199)
         mat(k,1660) = mat(k,1660) + rxt(k,455)*y(k,176)
         mat(k,915) = .800_r8*rxt(k,458)*y(k,123) + rxt(k,456)*y(k,199)
         mat(k,1060) = -(rxt(k,462)*y(k,218))
         mat(k,1665) = -rxt(k,462)*y(k,179)
         mat(k,1002) = .300_r8*rxt(k,441)*y(k,134)
         mat(k,953) = .300_r8*rxt(k,444)*y(k,134)
         mat(k,1894) = .900_r8*rxt(k,453)*y(k,222)
         mat(k,2225) = .300_r8*rxt(k,441)*y(k,5) + .300_r8*rxt(k,444)*y(k,109)
         mat(k,1804) = rxt(k,451)*y(k,222)
         mat(k,1108) = .900_r8*rxt(k,453)*y(k,123) + rxt(k,451)*y(k,199)
         mat(k,656) = -(rxt(k,429)*y(k,218))
         mat(k,1631) = -rxt(k,429)*y(k,180)
         mat(k,2059) = rxt(k,427)*y(k,224)
         mat(k,723) = rxt(k,427)*y(k,204)
         mat(k,188) = -(rxt(k,430)*y(k,218))
         mat(k,1564) = -rxt(k,430)*y(k,181)
         mat(k,204) = -(rxt(k,396)*y(k,218))
         mat(k,1567) = -rxt(k,396)*y(k,182)
         mat(k,2028) = rxt(k,393)*y(k,226)
         mat(k,1167) = rxt(k,393)*y(k,204)
         mat(k,209) = -(rxt(k,432)*y(k,218))
         mat(k,1568) = -rxt(k,432)*y(k,183)
         mat(k,694) = -(rxt(k,435)*y(k,218))
         mat(k,1635) = -rxt(k,435)*y(k,184)
         mat(k,2062) = rxt(k,433)*y(k,227)
         mat(k,747) = rxt(k,433)*y(k,204)
         mat(k,217) = -(rxt(k,438)*y(k,218))
         mat(k,1569) = -rxt(k,438)*y(k,185)
         mat(k,210) = .150_r8*rxt(k,432)*y(k,218)
         mat(k,1569) = mat(k,1569) + .150_r8*rxt(k,432)*y(k,183)
         mat(k,428) = -(rxt(k,439)*y(k,218))
         mat(k,1602) = -rxt(k,439)*y(k,186)
         mat(k,2042) = rxt(k,436)*y(k,229)
         mat(k,499) = rxt(k,436)*y(k,204)
         mat(k,519) = -(rxt(k,397)*y(k,204) + rxt(k,398)*y(k,123) + rxt(k,426) &
                      *y(k,124))
         mat(k,2052) = -rxt(k,397)*y(k,189)
         mat(k,1867) = -rxt(k,398)*y(k,189)
         mat(k,1938) = -rxt(k,426)*y(k,189)
         mat(k,237) = rxt(k,403)*y(k,218)
         mat(k,1615) = rxt(k,403)*y(k,21)
         mat(k,888) = -(rxt(k,358)*y(k,204) + (rxt(k,359) + rxt(k,360)) * y(k,123))
         mat(k,2078) = -rxt(k,358)*y(k,190)
         mat(k,1887) = -(rxt(k,359) + rxt(k,360)) * y(k,190)
         mat(k,636) = rxt(k,361)*y(k,218)
         mat(k,228) = rxt(k,362)*y(k,218)
         mat(k,1654) = rxt(k,361)*y(k,2) + rxt(k,362)*y(k,14)
         mat(k,55) = -(rxt(k,495)*y(k,204) + rxt(k,496)*y(k,123))
         mat(k,2020) = -rxt(k,495)*y(k,191)
         mat(k,1845) = -rxt(k,496)*y(k,191)
         mat(k,990) = rxt(k,498)*y(k,218)
         mat(k,1547) = rxt(k,498)*y(k,5)
         mat(k,484) = -(rxt(k,400)*y(k,204) + rxt(k,401)*y(k,123))
         mat(k,2050) = -rxt(k,400)*y(k,192)
         mat(k,1863) = -rxt(k,401)*y(k,192)
         mat(k,166) = .350_r8*rxt(k,399)*y(k,218)
         mat(k,418) = rxt(k,402)*y(k,218)
         mat(k,1610) = .350_r8*rxt(k,399)*y(k,6) + rxt(k,402)*y(k,7)
         mat(k,61) = -(rxt(k,500)*y(k,204) + rxt(k,501)*y(k,123))
         mat(k,2021) = -rxt(k,500)*y(k,193)
         mat(k,1846) = -rxt(k,501)*y(k,193)
         mat(k,162) = rxt(k,499)*y(k,218)
         mat(k,1548) = rxt(k,499)*y(k,6)
         mat(k,436) = -(rxt(k,404)*y(k,204) + rxt(k,406)*y(k,123))
         mat(k,2043) = -rxt(k,404)*y(k,194)
         mat(k,1858) = -rxt(k,406)*y(k,194)
         mat(k,338) = rxt(k,405)*y(k,218)
         mat(k,191) = .070_r8*rxt(k,430)*y(k,218)
         mat(k,211) = .060_r8*rxt(k,432)*y(k,218)
         mat(k,1603) = rxt(k,405)*y(k,22) + .070_r8*rxt(k,430)*y(k,181) &
                      + .060_r8*rxt(k,432)*y(k,183)
         mat(k,824) = -(4._r8*rxt(k,281)*y(k,195) + rxt(k,282)*y(k,199) + rxt(k,283) &
                      *y(k,204) + rxt(k,284)*y(k,123))
         mat(k,1797) = -rxt(k,282)*y(k,195)
         mat(k,2075) = -rxt(k,283)*y(k,195)
         mat(k,1884) = -rxt(k,284)*y(k,195)
         mat(k,343) = .500_r8*rxt(k,286)*y(k,218)
         mat(k,306) = rxt(k,287)*y(k,55) + rxt(k,288)*y(k,218)
         mat(k,1987) = rxt(k,287)*y(k,27)
         mat(k,1648) = .500_r8*rxt(k,286)*y(k,26) + rxt(k,288)*y(k,27)
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
         mat(k,788) = -(rxt(k,310)*y(k,199) + rxt(k,311)*y(k,204) + rxt(k,312) &
                      *y(k,123))
         mat(k,1794) = -rxt(k,310)*y(k,196)
         mat(k,2071) = -rxt(k,311)*y(k,196)
         mat(k,1882) = -rxt(k,312)*y(k,196)
         mat(k,411) = rxt(k,313)*y(k,218)
         mat(k,111) = rxt(k,314)*y(k,218)
         mat(k,1643) = rxt(k,313)*y(k,29) + rxt(k,314)*y(k,30)
         mat(k,626) = -(rxt(k,407)*y(k,204) + rxt(k,408)*y(k,123))
         mat(k,2056) = -rxt(k,407)*y(k,197)
         mat(k,1873) = -rxt(k,408)*y(k,197)
         mat(k,268) = rxt(k,409)*y(k,218)
         mat(k,1873) = mat(k,1873) + rxt(k,398)*y(k,189)
         mat(k,2210) = rxt(k,424)*y(k,141)
         mat(k,465) = rxt(k,424)*y(k,134)
         mat(k,520) = rxt(k,398)*y(k,123) + .400_r8*rxt(k,397)*y(k,204)
         mat(k,2056) = mat(k,2056) + .400_r8*rxt(k,397)*y(k,189)
         mat(k,1628) = rxt(k,409)*y(k,31)
         mat(k,1387) = -(4._r8*rxt(k,292)*y(k,198) + rxt(k,293)*y(k,199) + rxt(k,294) &
                      *y(k,204) + rxt(k,295)*y(k,123) + rxt(k,306)*y(k,124) + rxt(k,333) &
                      *y(k,210) + rxt(k,366)*y(k,206) + rxt(k,371)*y(k,207) + rxt(k,380) &
                      *y(k,100) + rxt(k,391)*y(k,226))
         mat(k,1821) = -rxt(k,293)*y(k,198)
         mat(k,2101) = -rxt(k,294)*y(k,198)
         mat(k,1912) = -rxt(k,295)*y(k,198)
         mat(k,1956) = -rxt(k,306)*y(k,198)
         mat(k,1314) = -rxt(k,333)*y(k,198)
         mat(k,1261) = -rxt(k,366)*y(k,198)
         mat(k,1293) = -rxt(k,371)*y(k,198)
         mat(k,1198) = -rxt(k,380)*y(k,198)
         mat(k,1176) = -rxt(k,391)*y(k,198)
         mat(k,1008) = .060_r8*rxt(k,441)*y(k,134)
         mat(k,1088) = rxt(k,289)*y(k,125) + rxt(k,290)*y(k,218)
         mat(k,1223) = rxt(k,315)*y(k,125) + rxt(k,316)*y(k,218)
         mat(k,613) = .500_r8*rxt(k,297)*y(k,218)
         mat(k,844) = .080_r8*rxt(k,386)*y(k,134)
         mat(k,1214) = .100_r8*rxt(k,339)*y(k,134)
         mat(k,958) = .060_r8*rxt(k,444)*y(k,134)
         mat(k,1335) = .280_r8*rxt(k,353)*y(k,134)
         mat(k,1912) = mat(k,1912) + .530_r8*rxt(k,337)*y(k,210) + rxt(k,346)*y(k,212) &
                      + rxt(k,349)*y(k,214) + rxt(k,324)*y(k,221)
         mat(k,1768) = rxt(k,289)*y(k,44) + rxt(k,315)*y(k,48) + .530_r8*rxt(k,336) &
                      *y(k,210) + rxt(k,347)*y(k,212)
         mat(k,2241) = .060_r8*rxt(k,441)*y(k,5) + .080_r8*rxt(k,386)*y(k,97) &
                      + .100_r8*rxt(k,339)*y(k,104) + .060_r8*rxt(k,444)*y(k,109) &
                      + .280_r8*rxt(k,353)*y(k,110)
         mat(k,1063) = .650_r8*rxt(k,462)*y(k,218)
         mat(k,1387) = mat(k,1387) + .530_r8*rxt(k,333)*y(k,210)
         mat(k,1821) = mat(k,1821) + .260_r8*rxt(k,334)*y(k,210) + rxt(k,343)*y(k,212) &
                      + .300_r8*rxt(k,322)*y(k,221)
         mat(k,2101) = mat(k,2101) + .450_r8*rxt(k,344)*y(k,212) + .200_r8*rxt(k,348) &
                      *y(k,214) + .150_r8*rxt(k,323)*y(k,221)
         mat(k,1314) = mat(k,1314) + .530_r8*rxt(k,337)*y(k,123) + .530_r8*rxt(k,336) &
                      *y(k,125) + .530_r8*rxt(k,333)*y(k,198) + .260_r8*rxt(k,334) &
                      *y(k,199)
         mat(k,1356) = rxt(k,346)*y(k,123) + rxt(k,347)*y(k,125) + rxt(k,343)*y(k,199) &
                      + .450_r8*rxt(k,344)*y(k,204) + 4.000_r8*rxt(k,345)*y(k,212)
         mat(k,680) = rxt(k,349)*y(k,123) + .200_r8*rxt(k,348)*y(k,204)
         mat(k,1684) = rxt(k,290)*y(k,44) + rxt(k,316)*y(k,48) + .500_r8*rxt(k,297) &
                      *y(k,50) + .650_r8*rxt(k,462)*y(k,179)
         mat(k,1133) = rxt(k,324)*y(k,123) + .300_r8*rxt(k,322)*y(k,199) &
                      + .150_r8*rxt(k,323)*y(k,204)
         mat(k,1830) = -(rxt(k,183)*y(k,58) + (4._r8*rxt(k,260) + 4._r8*rxt(k,261) &
                      ) * y(k,199) + rxt(k,262)*y(k,204) + rxt(k,263)*y(k,123) &
                      + rxt(k,282)*y(k,195) + rxt(k,293)*y(k,198) + rxt(k,310) &
                      *y(k,196) + rxt(k,322)*y(k,221) + rxt(k,334)*y(k,210) + rxt(k,343) &
                      *y(k,212) + rxt(k,367)*y(k,206) + rxt(k,372)*y(k,207) + rxt(k,381) &
                      *y(k,100) + rxt(k,392)*y(k,226) + rxt(k,446)*y(k,216) + rxt(k,451) &
                      *y(k,222) + rxt(k,456)*y(k,223))
         mat(k,1721) = -rxt(k,183)*y(k,199)
         mat(k,2112) = -rxt(k,262)*y(k,199)
         mat(k,1922) = -rxt(k,263)*y(k,199)
         mat(k,828) = -rxt(k,282)*y(k,199)
         mat(k,1393) = -rxt(k,293)*y(k,199)
         mat(k,793) = -rxt(k,310)*y(k,199)
         mat(k,1136) = -rxt(k,322)*y(k,199)
         mat(k,1319) = -rxt(k,334)*y(k,199)
         mat(k,1361) = -rxt(k,343)*y(k,199)
         mat(k,1266) = -rxt(k,367)*y(k,199)
         mat(k,1298) = -rxt(k,372)*y(k,199)
         mat(k,1203) = -rxt(k,381)*y(k,199)
         mat(k,1180) = -rxt(k,392)*y(k,199)
         mat(k,1054) = -rxt(k,446)*y(k,199)
         mat(k,1117) = -rxt(k,451)*y(k,199)
         mat(k,919) = -rxt(k,456)*y(k,199)
         mat(k,1035) = .280_r8*rxt(k,309)*y(k,134)
         mat(k,688) = rxt(k,296)*y(k,218)
         mat(k,459) = .700_r8*rxt(k,265)*y(k,218)
         mat(k,1438) = rxt(k,177)*y(k,55) + rxt(k,233)*y(k,72) + rxt(k,272)*y(k,217) &
                      + rxt(k,266)*y(k,218)
         mat(k,2005) = rxt(k,177)*y(k,53)
         mat(k,872) = rxt(k,233)*y(k,53)
         mat(k,848) = .050_r8*rxt(k,386)*y(k,134)
         mat(k,1203) = mat(k,1203) + rxt(k,380)*y(k,198)
         mat(k,1922) = mat(k,1922) + rxt(k,295)*y(k,198) + .830_r8*rxt(k,412)*y(k,200) &
                      + .170_r8*rxt(k,418)*y(k,213)
         mat(k,2251) = .280_r8*rxt(k,309)*y(k,28) + .050_r8*rxt(k,386)*y(k,97)
         mat(k,1393) = mat(k,1393) + rxt(k,380)*y(k,100) + rxt(k,295)*y(k,123) &
                      + 4.000_r8*rxt(k,292)*y(k,198) + .900_r8*rxt(k,293)*y(k,199) &
                      + .490_r8*rxt(k,294)*y(k,204) + rxt(k,366)*y(k,206) + rxt(k,371) &
                      *y(k,207) + rxt(k,333)*y(k,210) + rxt(k,342)*y(k,212) &
                      + rxt(k,391)*y(k,226)
         mat(k,1830) = mat(k,1830) + .900_r8*rxt(k,293)*y(k,198)
         mat(k,764) = .830_r8*rxt(k,412)*y(k,123) + .330_r8*rxt(k,411)*y(k,204)
         mat(k,2112) = mat(k,2112) + .490_r8*rxt(k,294)*y(k,198) + .330_r8*rxt(k,411) &
                      *y(k,200) + .070_r8*rxt(k,417)*y(k,213)
         mat(k,1266) = mat(k,1266) + rxt(k,366)*y(k,198)
         mat(k,1298) = mat(k,1298) + rxt(k,371)*y(k,198)
         mat(k,1319) = mat(k,1319) + rxt(k,333)*y(k,198)
         mat(k,1361) = mat(k,1361) + rxt(k,342)*y(k,198)
         mat(k,879) = .170_r8*rxt(k,418)*y(k,123) + .070_r8*rxt(k,417)*y(k,204)
         mat(k,1531) = rxt(k,272)*y(k,53)
         mat(k,1695) = rxt(k,296)*y(k,49) + .700_r8*rxt(k,265)*y(k,52) + rxt(k,266) &
                      *y(k,53)
         mat(k,1180) = mat(k,1180) + rxt(k,391)*y(k,198)
         mat(k,760) = -(rxt(k,411)*y(k,204) + rxt(k,412)*y(k,123) + rxt(k,413) &
                      *y(k,124))
         mat(k,2068) = -rxt(k,411)*y(k,200)
         mat(k,1880) = -rxt(k,412)*y(k,200)
         mat(k,1944) = -rxt(k,413)*y(k,200)
         mat(k,565) = -((rxt(k,330) + rxt(k,331)) * y(k,123))
         mat(k,1869) = -(rxt(k,330) + rxt(k,331)) * y(k,201)
         mat(k,351) = rxt(k,329)*y(k,218)
         mat(k,1620) = rxt(k,329)*y(k,15)
         mat(k,1854) = .750_r8*rxt(k,299)*y(k,203)
         mat(k,706) = .750_r8*rxt(k,299)*y(k,123)
         mat(k,707) = -(rxt(k,298)*y(k,204) + rxt(k,299)*y(k,123))
         mat(k,2063) = -rxt(k,298)*y(k,203)
         mat(k,1876) = -rxt(k,299)*y(k,203)
         mat(k,550) = rxt(k,305)*y(k,218)
         mat(k,1636) = rxt(k,305)*y(k,24)
         mat(k,2116) = -((rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,75) + rxt(k,140) &
                      *y(k,133) + rxt(k,141)*y(k,134) + rxt(k,145)*y(k,218) &
                      + 4._r8*rxt(k,150)*y(k,204) + rxt(k,160)*y(k,125) + rxt(k,165) &
                      *y(k,123) + rxt(k,170)*y(k,124) + (rxt(k,180) + rxt(k,181) &
                      ) * y(k,55) + rxt(k,187)*y(k,58) + rxt(k,213)*y(k,16) + rxt(k,219) &
                      *y(k,18) + rxt(k,256)*y(k,41) + rxt(k,262)*y(k,199) + rxt(k,269) &
                      *y(k,205) + rxt(k,283)*y(k,195) + rxt(k,294)*y(k,198) + rxt(k,298) &
                      *y(k,203) + rxt(k,311)*y(k,196) + rxt(k,319)*y(k,220) + rxt(k,323) &
                      *y(k,221) + rxt(k,335)*y(k,210) + rxt(k,344)*y(k,212) + rxt(k,348) &
                      *y(k,214) + rxt(k,358)*y(k,190) + rxt(k,368)*y(k,206) + rxt(k,373) &
                      *y(k,207) + rxt(k,382)*y(k,100) + rxt(k,393)*y(k,226) + rxt(k,397) &
                      *y(k,189) + rxt(k,400)*y(k,192) + rxt(k,404)*y(k,194) + rxt(k,407) &
                      *y(k,197) + rxt(k,411)*y(k,200) + rxt(k,414)*y(k,211) + rxt(k,417) &
                      *y(k,213) + rxt(k,420)*y(k,219) + rxt(k,427)*y(k,224) + rxt(k,433) &
                      *y(k,227) + rxt(k,436)*y(k,229) + rxt(k,447)*y(k,216) + rxt(k,452) &
                      *y(k,222) + rxt(k,457)*y(k,223))
         mat(k,1454) = -(rxt(k,136) + rxt(k,137) + rxt(k,138)) * y(k,204)
         mat(k,2194) = -rxt(k,140)*y(k,204)
         mat(k,2255) = -rxt(k,141)*y(k,204)
         mat(k,1699) = -rxt(k,145)*y(k,204)
         mat(k,1782) = -rxt(k,160)*y(k,204)
         mat(k,1926) = -rxt(k,165)*y(k,204)
         mat(k,1970) = -rxt(k,170)*y(k,204)
         mat(k,2009) = -(rxt(k,180) + rxt(k,181)) * y(k,204)
         mat(k,1725) = -rxt(k,187)*y(k,204)
         mat(k,1421) = -rxt(k,213)*y(k,204)
         mat(k,2163) = -rxt(k,219)*y(k,204)
         mat(k,1494) = -rxt(k,256)*y(k,204)
         mat(k,1834) = -rxt(k,262)*y(k,204)
         mat(k,447) = -rxt(k,269)*y(k,204)
         mat(k,831) = -rxt(k,283)*y(k,204)
         mat(k,1396) = -rxt(k,294)*y(k,204)
         mat(k,713) = -rxt(k,298)*y(k,204)
         mat(k,796) = -rxt(k,311)*y(k,204)
         mat(k,776) = -rxt(k,319)*y(k,204)
         mat(k,1139) = -rxt(k,323)*y(k,204)
         mat(k,1322) = -rxt(k,335)*y(k,204)
         mat(k,1364) = -rxt(k,344)*y(k,204)
         mat(k,684) = -rxt(k,348)*y(k,204)
         mat(k,897) = -rxt(k,358)*y(k,204)
         mat(k,1269) = -rxt(k,368)*y(k,204)
         mat(k,1301) = -rxt(k,373)*y(k,204)
         mat(k,1206) = -rxt(k,382)*y(k,204)
         mat(k,1183) = -rxt(k,393)*y(k,204)
         mat(k,524) = -rxt(k,397)*y(k,204)
         mat(k,490) = -rxt(k,400)*y(k,204)
         mat(k,441) = -rxt(k,404)*y(k,204)
         mat(k,630) = -rxt(k,407)*y(k,204)
         mat(k,767) = -rxt(k,411)*y(k,204)
         mat(k,719) = -rxt(k,414)*y(k,204)
         mat(k,882) = -rxt(k,417)*y(k,204)
         mat(k,454) = -rxt(k,420)*y(k,204)
         mat(k,734) = -rxt(k,427)*y(k,204)
         mat(k,759) = -rxt(k,433)*y(k,204)
         mat(k,506) = -rxt(k,436)*y(k,204)
         mat(k,1057) = -rxt(k,447)*y(k,204)
         mat(k,1120) = -rxt(k,452)*y(k,204)
         mat(k,922) = -rxt(k,457)*y(k,204)
         mat(k,1016) = .570_r8*rxt(k,441)*y(k,134)
         mat(k,168) = .650_r8*rxt(k,399)*y(k,218)
         mat(k,1421) = mat(k,1421) + rxt(k,212)*y(k,41)
         mat(k,2163) = mat(k,2163) + rxt(k,224)*y(k,218)
         mat(k,297) = .350_r8*rxt(k,278)*y(k,218)
         mat(k,555) = .130_r8*rxt(k,280)*y(k,134)
         mat(k,265) = rxt(k,285)*y(k,218)
         mat(k,1038) = .280_r8*rxt(k,309)*y(k,134)
         mat(k,1494) = mat(k,1494) + rxt(k,212)*y(k,16) + rxt(k,176)*y(k,55) &
                      + rxt(k,257)*y(k,125) + rxt(k,258)*y(k,133)
         mat(k,601) = rxt(k,241)*y(k,55) + rxt(k,242)*y(k,218)
         mat(k,371) = rxt(k,244)*y(k,55) + rxt(k,245)*y(k,218)
         mat(k,105) = rxt(k,291)*y(k,218)
         mat(k,801) = rxt(k,264)*y(k,218)
         mat(k,1440) = rxt(k,273)*y(k,217)
         mat(k,2009) = mat(k,2009) + rxt(k,176)*y(k,41) + rxt(k,241)*y(k,42) &
                      + rxt(k,244)*y(k,45) + rxt(k,179)*y(k,78)
         mat(k,1725) = mat(k,1725) + rxt(k,183)*y(k,199) + rxt(k,194)*y(k,218)
         mat(k,1126) = rxt(k,276)*y(k,218)
         mat(k,199) = .730_r8*rxt(k,410)*y(k,218)
         mat(k,288) = .500_r8*rxt(k,478)*y(k,218)
         mat(k,1101) = rxt(k,302)*y(k,218)
         mat(k,983) = rxt(k,303)*y(k,218)
         mat(k,607) = rxt(k,179)*y(k,55) + rxt(k,135)*y(k,133) + rxt(k,144)*y(k,218)
         mat(k,183) = rxt(k,267)*y(k,218)
         mat(k,933) = rxt(k,268)*y(k,218)
         mat(k,1164) = rxt(k,332)*y(k,218)
         mat(k,1146) = rxt(k,317)*y(k,218)
         mat(k,849) = .370_r8*rxt(k,386)*y(k,134)
         mat(k,590) = .300_r8*rxt(k,377)*y(k,218)
         mat(k,540) = rxt(k,378)*y(k,218)
         mat(k,1206) = mat(k,1206) + rxt(k,383)*y(k,123) + rxt(k,384)*y(k,125) &
                      + rxt(k,380)*y(k,198) + 1.200_r8*rxt(k,381)*y(k,199)
         mat(k,409) = rxt(k,385)*y(k,218)
         mat(k,1217) = .140_r8*rxt(k,339)*y(k,134)
         mat(k,315) = .200_r8*rxt(k,341)*y(k,218)
         mat(k,581) = .500_r8*rxt(k,352)*y(k,218)
         mat(k,966) = .570_r8*rxt(k,444)*y(k,134)
         mat(k,1344) = .280_r8*rxt(k,353)*y(k,134)
         mat(k,379) = rxt(k,389)*y(k,218)
         mat(k,1084) = rxt(k,390)*y(k,218)
         mat(k,1926) = mat(k,1926) + rxt(k,383)*y(k,100) + rxt(k,359)*y(k,190) &
                      + rxt(k,401)*y(k,192) + rxt(k,406)*y(k,194) + rxt(k,284) &
                      *y(k,195) + rxt(k,312)*y(k,196) + rxt(k,263)*y(k,199) &
                      + .170_r8*rxt(k,412)*y(k,200) + rxt(k,330)*y(k,201) &
                      + .250_r8*rxt(k,299)*y(k,203) + rxt(k,271)*y(k,205) &
                      + .920_r8*rxt(k,369)*y(k,206) + .920_r8*rxt(k,375)*y(k,207) &
                      + .470_r8*rxt(k,337)*y(k,210) + .400_r8*rxt(k,415)*y(k,211) &
                      + .830_r8*rxt(k,418)*y(k,213) + rxt(k,421)*y(k,219) + rxt(k,320) &
                      *y(k,220) + .900_r8*rxt(k,453)*y(k,222) + .800_r8*rxt(k,458) &
                      *y(k,223) + rxt(k,428)*y(k,224) + rxt(k,394)*y(k,226) &
                      + rxt(k,434)*y(k,227) + rxt(k,437)*y(k,229)
         mat(k,1782) = mat(k,1782) + rxt(k,257)*y(k,41) + rxt(k,384)*y(k,100) &
                      + rxt(k,370)*y(k,206) + rxt(k,376)*y(k,207) + .470_r8*rxt(k,336) &
                      *y(k,210) + rxt(k,163)*y(k,218) + rxt(k,395)*y(k,226)
         mat(k,2194) = mat(k,2194) + rxt(k,258)*y(k,41) + rxt(k,135)*y(k,78)
         mat(k,2255) = mat(k,2255) + .570_r8*rxt(k,441)*y(k,5) + .130_r8*rxt(k,280) &
                      *y(k,24) + .280_r8*rxt(k,309)*y(k,28) + .370_r8*rxt(k,386) &
                      *y(k,97) + .140_r8*rxt(k,339)*y(k,104) + .570_r8*rxt(k,444) &
                      *y(k,109) + .280_r8*rxt(k,353)*y(k,110) + rxt(k,147)*y(k,218)
         mat(k,177) = .800_r8*rxt(k,422)*y(k,218)
         mat(k,855) = rxt(k,468)*y(k,218)
         mat(k,1067) = .200_r8*rxt(k,462)*y(k,218)
         mat(k,194) = .280_r8*rxt(k,430)*y(k,218)
         mat(k,216) = .380_r8*rxt(k,432)*y(k,218)
         mat(k,221) = .630_r8*rxt(k,438)*y(k,218)
         mat(k,897) = mat(k,897) + rxt(k,359)*y(k,123)
         mat(k,490) = mat(k,490) + rxt(k,401)*y(k,123)
         mat(k,441) = mat(k,441) + rxt(k,406)*y(k,123)
         mat(k,831) = mat(k,831) + rxt(k,284)*y(k,123) + 2.400_r8*rxt(k,281)*y(k,195) &
                      + rxt(k,282)*y(k,199)
         mat(k,796) = mat(k,796) + rxt(k,312)*y(k,123) + rxt(k,310)*y(k,199)
         mat(k,1396) = mat(k,1396) + rxt(k,380)*y(k,100) + .900_r8*rxt(k,293)*y(k,199) &
                      + rxt(k,366)*y(k,206) + rxt(k,371)*y(k,207) + .470_r8*rxt(k,333) &
                      *y(k,210) + rxt(k,391)*y(k,226)
         mat(k,1834) = mat(k,1834) + rxt(k,183)*y(k,58) + 1.200_r8*rxt(k,381)*y(k,100) &
                      + rxt(k,263)*y(k,123) + rxt(k,282)*y(k,195) + rxt(k,310) &
                      *y(k,196) + .900_r8*rxt(k,293)*y(k,198) + 4.000_r8*rxt(k,260) &
                      *y(k,199) + rxt(k,367)*y(k,206) + rxt(k,372)*y(k,207) &
                      + .730_r8*rxt(k,334)*y(k,210) + rxt(k,343)*y(k,212) &
                      + .500_r8*rxt(k,446)*y(k,216) + .300_r8*rxt(k,322)*y(k,221) &
                      + rxt(k,451)*y(k,222) + rxt(k,456)*y(k,223) + .800_r8*rxt(k,392) &
                      *y(k,226)
         mat(k,767) = mat(k,767) + .170_r8*rxt(k,412)*y(k,123) + .070_r8*rxt(k,411) &
                      *y(k,204)
         mat(k,572) = rxt(k,330)*y(k,123)
         mat(k,713) = mat(k,713) + .250_r8*rxt(k,299)*y(k,123)
         mat(k,2116) = mat(k,2116) + .070_r8*rxt(k,411)*y(k,200) + .160_r8*rxt(k,414) &
                      *y(k,211) + .330_r8*rxt(k,417)*y(k,213)
         mat(k,447) = mat(k,447) + rxt(k,271)*y(k,123)
         mat(k,1269) = mat(k,1269) + .920_r8*rxt(k,369)*y(k,123) + rxt(k,370)*y(k,125) &
                      + rxt(k,366)*y(k,198) + rxt(k,367)*y(k,199)
         mat(k,1301) = mat(k,1301) + .920_r8*rxt(k,375)*y(k,123) + rxt(k,376)*y(k,125) &
                      + rxt(k,371)*y(k,198) + rxt(k,372)*y(k,199)
         mat(k,1322) = mat(k,1322) + .470_r8*rxt(k,337)*y(k,123) + .470_r8*rxt(k,336) &
                      *y(k,125) + .470_r8*rxt(k,333)*y(k,198) + .730_r8*rxt(k,334) &
                      *y(k,199)
         mat(k,719) = mat(k,719) + .400_r8*rxt(k,415)*y(k,123) + .160_r8*rxt(k,414) &
                      *y(k,204)
         mat(k,1364) = mat(k,1364) + rxt(k,343)*y(k,199)
         mat(k,882) = mat(k,882) + .830_r8*rxt(k,418)*y(k,123) + .330_r8*rxt(k,417) &
                      *y(k,204)
         mat(k,1057) = mat(k,1057) + .500_r8*rxt(k,446)*y(k,199)
         mat(k,1535) = rxt(k,273)*y(k,53)
         mat(k,1699) = mat(k,1699) + .650_r8*rxt(k,399)*y(k,6) + rxt(k,224)*y(k,18) &
                      + .350_r8*rxt(k,278)*y(k,23) + rxt(k,285)*y(k,25) + rxt(k,242) &
                      *y(k,42) + rxt(k,245)*y(k,45) + rxt(k,291)*y(k,46) + rxt(k,264) &
                      *y(k,51) + rxt(k,194)*y(k,58) + rxt(k,276)*y(k,61) &
                      + .730_r8*rxt(k,410)*y(k,65) + .500_r8*rxt(k,478)*y(k,66) &
                      + rxt(k,302)*y(k,73) + rxt(k,303)*y(k,74) + rxt(k,144)*y(k,78) &
                      + rxt(k,267)*y(k,85) + rxt(k,268)*y(k,86) + rxt(k,332)*y(k,92) &
                      + rxt(k,317)*y(k,94) + .300_r8*rxt(k,377)*y(k,98) + rxt(k,378) &
                      *y(k,99) + rxt(k,385)*y(k,101) + .200_r8*rxt(k,341)*y(k,105) &
                      + .500_r8*rxt(k,352)*y(k,108) + rxt(k,389)*y(k,114) + rxt(k,390) &
                      *y(k,115) + rxt(k,163)*y(k,125) + rxt(k,147)*y(k,134) &
                      + .800_r8*rxt(k,422)*y(k,142) + rxt(k,468)*y(k,151) &
                      + .200_r8*rxt(k,462)*y(k,179) + .280_r8*rxt(k,430)*y(k,181) &
                      + .380_r8*rxt(k,432)*y(k,183) + .630_r8*rxt(k,438)*y(k,185)
         mat(k,454) = mat(k,454) + rxt(k,421)*y(k,123)
         mat(k,776) = mat(k,776) + rxt(k,320)*y(k,123)
         mat(k,1139) = mat(k,1139) + .300_r8*rxt(k,322)*y(k,199)
         mat(k,1120) = mat(k,1120) + .900_r8*rxt(k,453)*y(k,123) + rxt(k,451)*y(k,199)
         mat(k,922) = mat(k,922) + .800_r8*rxt(k,458)*y(k,123) + rxt(k,456)*y(k,199)
         mat(k,734) = mat(k,734) + rxt(k,428)*y(k,123)
         mat(k,1183) = mat(k,1183) + rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125) &
                      + rxt(k,391)*y(k,198) + .800_r8*rxt(k,392)*y(k,199)
         mat(k,759) = mat(k,759) + rxt(k,434)*y(k,123)
         mat(k,506) = mat(k,506) + rxt(k,437)*y(k,123)
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
         mat(k,442) = -(rxt(k,269)*y(k,204) + rxt(k,271)*y(k,123))
         mat(k,2044) = -rxt(k,269)*y(k,205)
         mat(k,1859) = -rxt(k,271)*y(k,205)
         mat(k,1477) = rxt(k,256)*y(k,204)
         mat(k,2044) = mat(k,2044) + rxt(k,256)*y(k,41)
         mat(k,1257) = -(rxt(k,366)*y(k,198) + rxt(k,367)*y(k,199) + rxt(k,368) &
                      *y(k,204) + rxt(k,369)*y(k,123) + rxt(k,370)*y(k,125))
         mat(k,1382) = -rxt(k,366)*y(k,206)
         mat(k,1816) = -rxt(k,367)*y(k,206)
         mat(k,2096) = -rxt(k,368)*y(k,206)
         mat(k,1907) = -rxt(k,369)*y(k,206)
         mat(k,1763) = -rxt(k,370)*y(k,206)
         mat(k,841) = .600_r8*rxt(k,387)*y(k,218)
         mat(k,1679) = .600_r8*rxt(k,387)*y(k,97)
         mat(k,1289) = -(rxt(k,371)*y(k,198) + rxt(k,372)*y(k,199) + rxt(k,373) &
                      *y(k,204) + rxt(k,375)*y(k,123) + rxt(k,376)*y(k,125))
         mat(k,1383) = -rxt(k,371)*y(k,207)
         mat(k,1817) = -rxt(k,372)*y(k,207)
         mat(k,2097) = -rxt(k,373)*y(k,207)
         mat(k,1908) = -rxt(k,375)*y(k,207)
         mat(k,1764) = -rxt(k,376)*y(k,207)
         mat(k,842) = .400_r8*rxt(k,387)*y(k,218)
         mat(k,1680) = .400_r8*rxt(k,387)*y(k,97)
         mat(k,67) = -(rxt(k,503)*y(k,204) + rxt(k,504)*y(k,123))
         mat(k,2022) = -rxt(k,503)*y(k,208)
         mat(k,1847) = -rxt(k,504)*y(k,208)
         mat(k,834) = rxt(k,506)*y(k,218)
         mat(k,1549) = rxt(k,506)*y(k,97)
         mat(k,73) = -(rxt(k,507)*y(k,204) + rxt(k,508)*y(k,123))
         mat(k,2023) = -rxt(k,507)*y(k,209)
         mat(k,1848) = -rxt(k,508)*y(k,209)
         mat(k,74) = rxt(k,509)*y(k,218)
         mat(k,1550) = rxt(k,509)*y(k,103)
         mat(k,1312) = -(rxt(k,333)*y(k,198) + rxt(k,334)*y(k,199) + rxt(k,335) &
                      *y(k,204) + rxt(k,336)*y(k,125) + (rxt(k,337) + rxt(k,338) &
                      ) * y(k,123))
         mat(k,1384) = -rxt(k,333)*y(k,210)
         mat(k,1818) = -rxt(k,334)*y(k,210)
         mat(k,2098) = -rxt(k,335)*y(k,210)
         mat(k,1765) = -rxt(k,336)*y(k,210)
         mat(k,1909) = -(rxt(k,337) + rxt(k,338)) * y(k,210)
         mat(k,1212) = .500_r8*rxt(k,340)*y(k,218)
         mat(k,312) = .200_r8*rxt(k,341)*y(k,218)
         mat(k,1332) = rxt(k,354)*y(k,218)
         mat(k,1681) = .500_r8*rxt(k,340)*y(k,104) + .200_r8*rxt(k,341)*y(k,105) &
                      + rxt(k,354)*y(k,110)
         mat(k,714) = -(rxt(k,414)*y(k,204) + rxt(k,415)*y(k,123) + rxt(k,416) &
                      *y(k,124))
         mat(k,2064) = -rxt(k,414)*y(k,211)
         mat(k,1877) = -rxt(k,415)*y(k,211)
         mat(k,1943) = -rxt(k,416)*y(k,211)
         mat(k,1355) = -(rxt(k,342)*y(k,198) + rxt(k,343)*y(k,199) + rxt(k,344) &
                      *y(k,204) + 4._r8*rxt(k,345)*y(k,212) + rxt(k,346)*y(k,123) &
                      + rxt(k,347)*y(k,125) + rxt(k,355)*y(k,124))
         mat(k,1386) = -rxt(k,342)*y(k,212)
         mat(k,1820) = -rxt(k,343)*y(k,212)
         mat(k,2100) = -rxt(k,344)*y(k,212)
         mat(k,1911) = -rxt(k,346)*y(k,212)
         mat(k,1767) = -rxt(k,347)*y(k,212)
         mat(k,1955) = -rxt(k,355)*y(k,212)
         mat(k,1213) = .500_r8*rxt(k,340)*y(k,218)
         mat(k,313) = .500_r8*rxt(k,341)*y(k,218)
         mat(k,1683) = .500_r8*rxt(k,340)*y(k,104) + .500_r8*rxt(k,341)*y(k,105)
         mat(k,874) = -(rxt(k,417)*y(k,204) + rxt(k,418)*y(k,123) + rxt(k,419) &
                      *y(k,124))
         mat(k,2077) = -rxt(k,417)*y(k,213)
         mat(k,1886) = -rxt(k,418)*y(k,213)
         mat(k,1948) = -rxt(k,419)*y(k,213)
         mat(k,678) = -(rxt(k,348)*y(k,204) + rxt(k,349)*y(k,123))
         mat(k,2060) = -rxt(k,348)*y(k,214)
         mat(k,1875) = -rxt(k,349)*y(k,214)
         mat(k,508) = rxt(k,350)*y(k,218)
         mat(k,317) = rxt(k,351)*y(k,218)
         mat(k,1633) = rxt(k,350)*y(k,106) + rxt(k,351)*y(k,107)
         mat(k,81) = -(rxt(k,511)*y(k,204) + rxt(k,512)*y(k,123))
         mat(k,2024) = -rxt(k,511)*y(k,215)
         mat(k,1849) = -rxt(k,512)*y(k,215)
         mat(k,941) = rxt(k,514)*y(k,218)
         mat(k,1552) = rxt(k,514)*y(k,109)
         mat(k,1046) = -(rxt(k,446)*y(k,199) + rxt(k,447)*y(k,204) + rxt(k,448) &
                      *y(k,123) + rxt(k,449)*y(k,125))
         mat(k,1803) = -rxt(k,446)*y(k,216)
         mat(k,2084) = -rxt(k,447)*y(k,216)
         mat(k,1893) = -rxt(k,448)*y(k,216)
         mat(k,1748) = -rxt(k,449)*y(k,216)
         mat(k,1001) = rxt(k,440)*y(k,125)
         mat(k,952) = rxt(k,443)*y(k,125)
         mat(k,1748) = mat(k,1748) + rxt(k,440)*y(k,5) + rxt(k,443)*y(k,109) &
                      + .500_r8*rxt(k,460)*y(k,178)
         mat(k,394) = rxt(k,450)*y(k,218)
         mat(k,971) = .500_r8*rxt(k,460)*y(k,125)
         mat(k,1664) = rxt(k,450)*y(k,127)
         mat(k,1527) = -(rxt(k,125)*y(k,76) + rxt(k,126)*y(k,230) + (rxt(k,129) &
                      + rxt(k,130)) * y(k,134) + (rxt(k,168) + rxt(k,169)) * y(k,112) &
                      + rxt(k,201)*y(k,32) + rxt(k,202)*y(k,33) + rxt(k,203)*y(k,35) &
                      + rxt(k,204)*y(k,36) + rxt(k,205)*y(k,37) + rxt(k,206)*y(k,38) &
                      + rxt(k,207)*y(k,39) + (rxt(k,208) + rxt(k,209)) * y(k,84) &
                      + rxt(k,228)*y(k,34) + rxt(k,229)*y(k,54) + rxt(k,230)*y(k,77) &
                      + (rxt(k,231) + rxt(k,232)) * y(k,80) + rxt(k,237)*y(k,63) &
                      + rxt(k,238)*y(k,64) + rxt(k,251)*y(k,40) + rxt(k,252)*y(k,42) &
                      + rxt(k,253)*y(k,81) + rxt(k,254)*y(k,82) + rxt(k,255)*y(k,83) &
                      + (rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,53) + rxt(k,275) &
                      *y(k,85))
         mat(k,1406) = -rxt(k,125)*y(k,217)
         mat(k,2273) = -rxt(k,126)*y(k,217)
         mat(k,2247) = -(rxt(k,129) + rxt(k,130)) * y(k,217)
         mat(k,185) = -(rxt(k,168) + rxt(k,169)) * y(k,217)
         mat(k,101) = -rxt(k,201)*y(k,217)
         mat(k,145) = -rxt(k,202)*y(k,217)
         mat(k,116) = -rxt(k,203)*y(k,217)
         mat(k,155) = -rxt(k,204)*y(k,217)
         mat(k,120) = -rxt(k,205)*y(k,217)
         mat(k,160) = -rxt(k,206)*y(k,217)
         mat(k,124) = -rxt(k,207)*y(k,217)
         mat(k,2131) = -(rxt(k,208) + rxt(k,209)) * y(k,217)
         mat(k,151) = -rxt(k,228)*y(k,217)
         mat(k,388) = -rxt(k,229)*y(k,217)
         mat(k,109) = -rxt(k,230)*y(k,217)
         mat(k,808) = -(rxt(k,231) + rxt(k,232)) * y(k,217)
         mat(k,241) = -rxt(k,237)*y(k,217)
         mat(k,249) = -rxt(k,238)*y(k,217)
         mat(k,470) = -rxt(k,251)*y(k,217)
         mat(k,597) = -rxt(k,252)*y(k,217)
         mat(k,244) = -rxt(k,253)*y(k,217)
         mat(k,254) = -rxt(k,254)*y(k,217)
         mat(k,301) = -rxt(k,255)*y(k,217)
         mat(k,1435) = -(rxt(k,272) + rxt(k,273) + rxt(k,274)) * y(k,217)
         mat(k,181) = -rxt(k,275)*y(k,217)
         mat(k,1692) = -(rxt(k,143)*y(k,76) + rxt(k,144)*y(k,78) + rxt(k,145)*y(k,204) &
                      + rxt(k,146)*y(k,133) + rxt(k,147)*y(k,134) + (4._r8*rxt(k,148) &
                      + 4._r8*rxt(k,149)) * y(k,218) + rxt(k,151)*y(k,89) + rxt(k,163) &
                      *y(k,125) + rxt(k,164)*y(k,111) + rxt(k,172)*y(k,124) + rxt(k,173) &
                      *y(k,88) + rxt(k,192)*y(k,59) + (rxt(k,194) + rxt(k,195) &
                      ) * y(k,58) + rxt(k,197)*y(k,84) + rxt(k,200)*y(k,91) + rxt(k,224) &
                      *y(k,18) + rxt(k,226)*y(k,80) + rxt(k,240)*y(k,40) + rxt(k,242) &
                      *y(k,42) + rxt(k,243)*y(k,43) + rxt(k,245)*y(k,45) + rxt(k,247) &
                      *y(k,54) + rxt(k,248)*y(k,81) + rxt(k,249)*y(k,82) + rxt(k,250) &
                      *y(k,83) + rxt(k,259)*y(k,41) + rxt(k,264)*y(k,51) + rxt(k,265) &
                      *y(k,52) + rxt(k,266)*y(k,53) + rxt(k,267)*y(k,85) + rxt(k,268) &
                      *y(k,86) + rxt(k,276)*y(k,61) + rxt(k,278)*y(k,23) + rxt(k,285) &
                      *y(k,25) + rxt(k,286)*y(k,26) + rxt(k,288)*y(k,27) + rxt(k,290) &
                      *y(k,44) + rxt(k,291)*y(k,46) + rxt(k,296)*y(k,49) + rxt(k,297) &
                      *y(k,50) + rxt(k,302)*y(k,73) + rxt(k,303)*y(k,74) + rxt(k,304) &
                      *y(k,139) + rxt(k,305)*y(k,24) + rxt(k,313)*y(k,29) + rxt(k,314) &
                      *y(k,30) + rxt(k,316)*y(k,48) + rxt(k,317)*y(k,94) + rxt(k,318) &
                      *y(k,126) + rxt(k,321)*y(k,146) + rxt(k,325)*y(k,147) + rxt(k,326) &
                      *y(k,28) + rxt(k,327)*y(k,47) + rxt(k,329)*y(k,15) + rxt(k,332) &
                      *y(k,92) + rxt(k,340)*y(k,104) + rxt(k,341)*y(k,105) + rxt(k,350) &
                      *y(k,106) + rxt(k,351)*y(k,107) + rxt(k,352)*y(k,108) + rxt(k,354) &
                      *y(k,110) + rxt(k,357)*y(k,1) + rxt(k,361)*y(k,2) + rxt(k,362) &
                      *y(k,14) + rxt(k,363)*y(k,93) + rxt(k,364)*y(k,95) + rxt(k,365) &
                      *y(k,96) + rxt(k,377)*y(k,98) + rxt(k,378)*y(k,99) + rxt(k,385) &
                      *y(k,101) + rxt(k,387)*y(k,97) + rxt(k,388)*y(k,102) + rxt(k,389) &
                      *y(k,114) + rxt(k,390)*y(k,115) + rxt(k,396)*y(k,182) + rxt(k,399) &
                      *y(k,6) + rxt(k,402)*y(k,7) + rxt(k,403)*y(k,21) + rxt(k,405) &
                      *y(k,22) + rxt(k,409)*y(k,31) + rxt(k,410)*y(k,65) + rxt(k,422) &
                      *y(k,142) + rxt(k,425)*y(k,143) + rxt(k,429)*y(k,180) + rxt(k,430) &
                      *y(k,181) + rxt(k,432)*y(k,183) + rxt(k,435)*y(k,184) + rxt(k,438) &
                      *y(k,185) + rxt(k,439)*y(k,186) + rxt(k,442)*y(k,5) + rxt(k,445) &
                      *y(k,109) + rxt(k,450)*y(k,127) + rxt(k,454)*y(k,175) + rxt(k,455) &
                      *y(k,176) + rxt(k,459)*y(k,177) + rxt(k,461)*y(k,178) + rxt(k,462) &
                      *y(k,179) + (rxt(k,464) + rxt(k,478)) * y(k,66) + rxt(k,466) &
                      *y(k,137) + rxt(k,468)*y(k,151) + rxt(k,472)*y(k,148) + rxt(k,477) &
                      *y(k,150) + rxt(k,480)*y(k,119))
         mat(k,1407) = -rxt(k,143)*y(k,218)
         mat(k,605) = -rxt(k,144)*y(k,218)
         mat(k,2109) = -rxt(k,145)*y(k,218)
         mat(k,2187) = -rxt(k,146)*y(k,218)
         mat(k,2248) = -rxt(k,147)*y(k,218)
         mat(k,423) = -rxt(k,151)*y(k,218)
         mat(k,1775) = -rxt(k,163)*y(k,218)
         mat(k,495) = -rxt(k,164)*y(k,218)
         mat(k,1963) = -rxt(k,172)*y(k,218)
         mat(k,1467) = -rxt(k,173)*y(k,218)
         mat(k,902) = -rxt(k,192)*y(k,218)
         mat(k,1718) = -(rxt(k,194) + rxt(k,195)) * y(k,218)
         mat(k,2132) = -rxt(k,197)*y(k,218)
         mat(k,815) = -rxt(k,200)*y(k,218)
         mat(k,2156) = -rxt(k,224)*y(k,218)
         mat(k,809) = -rxt(k,226)*y(k,218)
         mat(k,471) = -rxt(k,240)*y(k,218)
         mat(k,598) = -rxt(k,242)*y(k,218)
         mat(k,127) = -rxt(k,243)*y(k,218)
         mat(k,368) = -rxt(k,245)*y(k,218)
         mat(k,389) = -rxt(k,247)*y(k,218)
         mat(k,245) = -rxt(k,248)*y(k,218)
         mat(k,255) = -rxt(k,249)*y(k,218)
         mat(k,302) = -rxt(k,250)*y(k,218)
         mat(k,1488) = -rxt(k,259)*y(k,218)
         mat(k,800) = -rxt(k,264)*y(k,218)
         mat(k,458) = -rxt(k,265)*y(k,218)
         mat(k,1436) = -rxt(k,266)*y(k,218)
         mat(k,182) = -rxt(k,267)*y(k,218)
         mat(k,932) = -rxt(k,268)*y(k,218)
         mat(k,1125) = -rxt(k,276)*y(k,218)
         mat(k,296) = -rxt(k,278)*y(k,218)
         mat(k,264) = -rxt(k,285)*y(k,218)
         mat(k,345) = -rxt(k,286)*y(k,218)
         mat(k,307) = -rxt(k,288)*y(k,218)
         mat(k,1090) = -rxt(k,290)*y(k,218)
         mat(k,104) = -rxt(k,291)*y(k,218)
         mat(k,687) = -rxt(k,296)*y(k,218)
         mat(k,615) = -rxt(k,297)*y(k,218)
         mat(k,1100) = -rxt(k,302)*y(k,218)
         mat(k,982) = -rxt(k,303)*y(k,218)
         mat(k,529) = -rxt(k,304)*y(k,218)
         mat(k,554) = -rxt(k,305)*y(k,218)
         mat(k,413) = -rxt(k,313)*y(k,218)
         mat(k,112) = -rxt(k,314)*y(k,218)
         mat(k,1225) = -rxt(k,316)*y(k,218)
         mat(k,1145) = -rxt(k,317)*y(k,218)
         mat(k,862) = -rxt(k,318)*y(k,218)
         mat(k,546) = -rxt(k,321)*y(k,218)
         mat(k,402) = -rxt(k,325)*y(k,218)
         mat(k,1033) = -rxt(k,326)*y(k,218)
         mat(k,926) = -rxt(k,327)*y(k,218)
         mat(k,355) = -rxt(k,329)*y(k,218)
         mat(k,1159) = -rxt(k,332)*y(k,218)
         mat(k,1216) = -rxt(k,340)*y(k,218)
         mat(k,314) = -rxt(k,341)*y(k,218)
         mat(k,511) = -rxt(k,350)*y(k,218)
         mat(k,320) = -rxt(k,351)*y(k,218)
         mat(k,578) = -rxt(k,352)*y(k,218)
         mat(k,1339) = -rxt(k,354)*y(k,218)
         mat(k,674) = -rxt(k,357)*y(k,218)
         mat(k,641) = -rxt(k,361)*y(k,218)
         mat(k,229) = -rxt(k,362)*y(k,218)
         mat(k,225) = -rxt(k,363)*y(k,218)
         mat(k,349) = -rxt(k,364)*y(k,218)
         mat(k,138) = -rxt(k,365)*y(k,218)
         mat(k,588) = -rxt(k,377)*y(k,218)
         mat(k,539) = -rxt(k,378)*y(k,218)
         mat(k,407) = -rxt(k,385)*y(k,218)
         mat(k,846) = -rxt(k,387)*y(k,218)
         mat(k,741) = -rxt(k,388)*y(k,218)
         mat(k,378) = -rxt(k,389)*y(k,218)
         mat(k,1080) = -rxt(k,390)*y(k,218)
         mat(k,206) = -rxt(k,396)*y(k,218)
         mat(k,167) = -rxt(k,399)*y(k,218)
         mat(k,420) = -rxt(k,402)*y(k,218)
         mat(k,238) = -rxt(k,403)*y(k,218)
         mat(k,340) = -rxt(k,405)*y(k,218)
         mat(k,269) = -rxt(k,409)*y(k,218)
         mat(k,198) = -rxt(k,410)*y(k,218)
         mat(k,176) = -rxt(k,422)*y(k,218)
         mat(k,334) = -rxt(k,425)*y(k,218)
         mat(k,664) = -rxt(k,429)*y(k,218)
         mat(k,193) = -rxt(k,430)*y(k,218)
         mat(k,215) = -rxt(k,432)*y(k,218)
         mat(k,703) = -rxt(k,435)*y(k,218)
         mat(k,220) = -rxt(k,438)*y(k,218)
         mat(k,432) = -rxt(k,439)*y(k,218)
         mat(k,1011) = -rxt(k,442)*y(k,218)
         mat(k,961) = -rxt(k,445)*y(k,218)
         mat(k,396) = -rxt(k,450)*y(k,218)
         mat(k,651) = -rxt(k,454)*y(k,218)
         mat(k,621) = -rxt(k,455)*y(k,218)
         mat(k,480) = -rxt(k,459)*y(k,218)
         mat(k,975) = -rxt(k,461)*y(k,218)
         mat(k,1065) = -rxt(k,462)*y(k,218)
         mat(k,286) = -(rxt(k,464) + rxt(k,478)) * y(k,218)
         mat(k,364) = -rxt(k,466)*y(k,218)
         mat(k,854) = -rxt(k,468)*y(k,218)
         mat(k,515) = -rxt(k,472)*y(k,218)
         mat(k,1236) = -rxt(k,477)*y(k,218)
         mat(k,98) = -rxt(k,480)*y(k,218)
         mat(k,1011) = mat(k,1011) + .630_r8*rxt(k,441)*y(k,134)
         mat(k,296) = mat(k,296) + .650_r8*rxt(k,278)*y(k,218)
         mat(k,554) = mat(k,554) + .130_r8*rxt(k,280)*y(k,134)
         mat(k,345) = mat(k,345) + .500_r8*rxt(k,286)*y(k,218)
         mat(k,1033) = mat(k,1033) + .360_r8*rxt(k,309)*y(k,134)
         mat(k,1488) = mat(k,1488) + rxt(k,258)*y(k,133)
         mat(k,458) = mat(k,458) + .300_r8*rxt(k,265)*y(k,218)
         mat(k,1436) = mat(k,1436) + rxt(k,272)*y(k,217)
         mat(k,2002) = rxt(k,181)*y(k,204)
         mat(k,870) = rxt(k,235)*y(k,230)
         mat(k,1450) = rxt(k,142)*y(k,134) + 2.000_r8*rxt(k,137)*y(k,204)
         mat(k,1407) = mat(k,1407) + rxt(k,134)*y(k,133) + rxt(k,125)*y(k,217)
         mat(k,605) = mat(k,605) + rxt(k,135)*y(k,133)
         mat(k,809) = mat(k,809) + rxt(k,225)*y(k,133) + rxt(k,231)*y(k,217)
         mat(k,2132) = mat(k,2132) + rxt(k,196)*y(k,133) + rxt(k,208)*y(k,217)
         mat(k,182) = mat(k,182) + rxt(k,275)*y(k,217)
         mat(k,781) = rxt(k,227)*y(k,133)
         mat(k,815) = mat(k,815) + rxt(k,199)*y(k,133)
         mat(k,846) = mat(k,846) + .320_r8*rxt(k,386)*y(k,134)
         mat(k,741) = mat(k,741) + .600_r8*rxt(k,388)*y(k,218)
         mat(k,1216) = mat(k,1216) + .240_r8*rxt(k,339)*y(k,134)
         mat(k,314) = mat(k,314) + .100_r8*rxt(k,341)*y(k,218)
         mat(k,961) = mat(k,961) + .630_r8*rxt(k,444)*y(k,134)
         mat(k,1339) = mat(k,1339) + .360_r8*rxt(k,353)*y(k,134)
         mat(k,1919) = rxt(k,165)*y(k,204)
         mat(k,1775) = mat(k,1775) + rxt(k,160)*y(k,204)
         mat(k,2187) = mat(k,2187) + rxt(k,258)*y(k,41) + rxt(k,134)*y(k,76) &
                      + rxt(k,135)*y(k,78) + rxt(k,225)*y(k,80) + rxt(k,196)*y(k,84) &
                      + rxt(k,227)*y(k,90) + rxt(k,199)*y(k,91) + rxt(k,140)*y(k,204)
         mat(k,2248) = mat(k,2248) + .630_r8*rxt(k,441)*y(k,5) + .130_r8*rxt(k,280) &
                      *y(k,24) + .360_r8*rxt(k,309)*y(k,28) + rxt(k,142)*y(k,75) &
                      + .320_r8*rxt(k,386)*y(k,97) + .240_r8*rxt(k,339)*y(k,104) &
                      + .630_r8*rxt(k,444)*y(k,109) + .360_r8*rxt(k,353)*y(k,110) &
                      + rxt(k,141)*y(k,204)
         mat(k,546) = mat(k,546) + .500_r8*rxt(k,321)*y(k,218)
         mat(k,206) = mat(k,206) + .500_r8*rxt(k,396)*y(k,218)
         mat(k,521) = .400_r8*rxt(k,397)*y(k,204)
         mat(k,1391) = .490_r8*rxt(k,294)*y(k,204)
         mat(k,763) = .400_r8*rxt(k,411)*y(k,204)
         mat(k,2109) = mat(k,2109) + rxt(k,181)*y(k,55) + 2.000_r8*rxt(k,137)*y(k,75) &
                      + rxt(k,165)*y(k,123) + rxt(k,160)*y(k,125) + rxt(k,140) &
                      *y(k,133) + rxt(k,141)*y(k,134) + .400_r8*rxt(k,397)*y(k,189) &
                      + .490_r8*rxt(k,294)*y(k,198) + .400_r8*rxt(k,411)*y(k,200) &
                      + .450_r8*rxt(k,344)*y(k,212) + .400_r8*rxt(k,417)*y(k,213) &
                      + .200_r8*rxt(k,348)*y(k,214) + .150_r8*rxt(k,323)*y(k,221)
         mat(k,1359) = .450_r8*rxt(k,344)*y(k,204)
         mat(k,878) = .400_r8*rxt(k,417)*y(k,204)
         mat(k,681) = .200_r8*rxt(k,348)*y(k,204)
         mat(k,1528) = rxt(k,272)*y(k,53) + rxt(k,125)*y(k,76) + rxt(k,231)*y(k,80) &
                      + rxt(k,208)*y(k,84) + rxt(k,275)*y(k,85) + 2.000_r8*rxt(k,126) &
                      *y(k,230)
         mat(k,1692) = mat(k,1692) + .650_r8*rxt(k,278)*y(k,23) + .500_r8*rxt(k,286) &
                      *y(k,26) + .300_r8*rxt(k,265)*y(k,52) + .600_r8*rxt(k,388) &
                      *y(k,102) + .100_r8*rxt(k,341)*y(k,105) + .500_r8*rxt(k,321) &
                      *y(k,146) + .500_r8*rxt(k,396)*y(k,182)
         mat(k,1135) = .150_r8*rxt(k,323)*y(k,204)
         mat(k,2274) = rxt(k,235)*y(k,72) + 2.000_r8*rxt(k,126)*y(k,217)
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
         mat(k,449) = -(rxt(k,420)*y(k,204) + rxt(k,421)*y(k,123))
         mat(k,2045) = -rxt(k,420)*y(k,219)
         mat(k,1860) = -rxt(k,421)*y(k,219)
         mat(k,196) = .200_r8*rxt(k,410)*y(k,218)
         mat(k,174) = .140_r8*rxt(k,422)*y(k,218)
         mat(k,332) = rxt(k,425)*y(k,218)
         mat(k,1604) = .200_r8*rxt(k,410)*y(k,65) + .140_r8*rxt(k,422)*y(k,142) &
                      + rxt(k,425)*y(k,143)
         mat(k,769) = -(rxt(k,319)*y(k,204) + rxt(k,320)*y(k,123))
         mat(k,2069) = -rxt(k,319)*y(k,220)
         mat(k,1881) = -rxt(k,320)*y(k,220)
         mat(k,1021) = rxt(k,326)*y(k,218)
         mat(k,542) = .500_r8*rxt(k,321)*y(k,218)
         mat(k,1642) = rxt(k,326)*y(k,28) + .500_r8*rxt(k,321)*y(k,146)
         mat(k,1130) = -(rxt(k,322)*y(k,199) + rxt(k,323)*y(k,204) + rxt(k,324) &
                      *y(k,123))
         mat(k,1810) = -rxt(k,322)*y(k,221)
         mat(k,2090) = -rxt(k,323)*y(k,221)
         mat(k,1900) = -rxt(k,324)*y(k,221)
         mat(k,1006) = .060_r8*rxt(k,441)*y(k,134)
         mat(k,924) = rxt(k,327)*y(k,218)
         mat(k,956) = .060_r8*rxt(k,444)*y(k,134)
         mat(k,2230) = .060_r8*rxt(k,441)*y(k,5) + .060_r8*rxt(k,444)*y(k,109)
         mat(k,399) = rxt(k,325)*y(k,218)
         mat(k,1062) = .150_r8*rxt(k,462)*y(k,218)
         mat(k,1671) = rxt(k,327)*y(k,47) + rxt(k,325)*y(k,147) + .150_r8*rxt(k,462) &
                      *y(k,179)
         mat(k,1110) = -(rxt(k,451)*y(k,199) + rxt(k,452)*y(k,204) + rxt(k,453) &
                      *y(k,123))
         mat(k,1808) = -rxt(k,451)*y(k,222)
         mat(k,2088) = -rxt(k,452)*y(k,222)
         mat(k,1898) = -rxt(k,453)*y(k,222)
         mat(k,1753) = .500_r8*rxt(k,460)*y(k,178)
         mat(k,648) = rxt(k,454)*y(k,218)
         mat(k,973) = .500_r8*rxt(k,460)*y(k,125) + rxt(k,461)*y(k,218)
         mat(k,1669) = rxt(k,454)*y(k,175) + rxt(k,461)*y(k,178)
         mat(k,913) = -(rxt(k,456)*y(k,199) + rxt(k,457)*y(k,204) + rxt(k,458) &
                      *y(k,123))
         mat(k,1799) = -rxt(k,456)*y(k,223)
         mat(k,2079) = -rxt(k,457)*y(k,223)
         mat(k,1888) = -rxt(k,458)*y(k,223)
         mat(k,995) = rxt(k,442)*y(k,218)
         mat(k,946) = rxt(k,445)*y(k,218)
         mat(k,476) = rxt(k,459)*y(k,218)
         mat(k,1656) = rxt(k,442)*y(k,5) + rxt(k,445)*y(k,109) + rxt(k,459)*y(k,177)
         mat(k,725) = -(rxt(k,427)*y(k,204) + rxt(k,428)*y(k,123))
         mat(k,2065) = -rxt(k,427)*y(k,224)
         mat(k,1878) = -rxt(k,428)*y(k,224)
         mat(k,658) = rxt(k,429)*y(k,218)
         mat(k,192) = .650_r8*rxt(k,430)*y(k,218)
         mat(k,1638) = rxt(k,429)*y(k,180) + .650_r8*rxt(k,430)*y(k,181)
         mat(k,87) = -(rxt(k,517)*y(k,204) + rxt(k,518)*y(k,123))
         mat(k,2025) = -rxt(k,517)*y(k,225)
         mat(k,1850) = -rxt(k,518)*y(k,225)
         mat(k,187) = rxt(k,516)*y(k,218)
         mat(k,1553) = rxt(k,516)*y(k,181)
         mat(k,1174) = -(rxt(k,391)*y(k,198) + rxt(k,392)*y(k,199) + rxt(k,393) &
                      *y(k,204) + rxt(k,394)*y(k,123) + rxt(k,395)*y(k,125))
         mat(k,1378) = -rxt(k,391)*y(k,226)
         mat(k,1812) = -rxt(k,392)*y(k,226)
         mat(k,2092) = -rxt(k,393)*y(k,226)
         mat(k,1903) = -rxt(k,394)*y(k,226)
         mat(k,1758) = -rxt(k,395)*y(k,226)
         mat(k,224) = rxt(k,363)*y(k,218)
         mat(k,348) = rxt(k,364)*y(k,218)
         mat(k,137) = rxt(k,365)*y(k,218)
         mat(k,737) = .400_r8*rxt(k,388)*y(k,218)
         mat(k,205) = .500_r8*rxt(k,396)*y(k,218)
         mat(k,1674) = rxt(k,363)*y(k,93) + rxt(k,364)*y(k,95) + rxt(k,365)*y(k,96) &
                      + .400_r8*rxt(k,388)*y(k,102) + .500_r8*rxt(k,396)*y(k,182)
         mat(k,749) = -(rxt(k,433)*y(k,204) + rxt(k,434)*y(k,123))
         mat(k,2067) = -rxt(k,433)*y(k,227)
         mat(k,1879) = -rxt(k,434)*y(k,227)
         mat(k,212) = .560_r8*rxt(k,432)*y(k,218)
         mat(k,696) = rxt(k,435)*y(k,218)
         mat(k,1640) = .560_r8*rxt(k,432)*y(k,183) + rxt(k,435)*y(k,184)
         mat(k,93) = -(rxt(k,521)*y(k,204) + rxt(k,522)*y(k,123))
         mat(k,2026) = -rxt(k,521)*y(k,228)
         mat(k,1851) = -rxt(k,522)*y(k,228)
         mat(k,207) = rxt(k,520)*y(k,218)
         mat(k,1554) = rxt(k,520)*y(k,183)
         mat(k,500) = -(rxt(k,436)*y(k,204) + rxt(k,437)*y(k,123))
         mat(k,2051) = -rxt(k,436)*y(k,229)
         mat(k,1865) = -rxt(k,437)*y(k,229)
         mat(k,219) = .300_r8*rxt(k,438)*y(k,218)
         mat(k,429) = rxt(k,439)*y(k,218)
         mat(k,1612) = .300_r8*rxt(k,438)*y(k,185) + rxt(k,439)*y(k,186)
         mat(k,2286) = -(rxt(k,126)*y(k,217) + rxt(k,235)*y(k,72) + rxt(k,479) &
                      *y(k,152))
         mat(k,1540) = -rxt(k,126)*y(k,230)
         mat(k,873) = -rxt(k,235)*y(k,230)
         mat(k,261) = -rxt(k,479)*y(k,230)
         mat(k,310) = rxt(k,288)*y(k,218)
         mat(k,415) = rxt(k,313)*y(k,218)
         mat(k,113) = rxt(k,314)*y(k,218)
         mat(k,474) = rxt(k,240)*y(k,218)
         mat(k,1499) = rxt(k,259)*y(k,218)
         mat(k,603) = rxt(k,242)*y(k,218)
         mat(k,129) = rxt(k,243)*y(k,218)
         mat(k,1094) = rxt(k,290)*y(k,218)
         mat(k,373) = rxt(k,245)*y(k,218)
         mat(k,928) = rxt(k,327)*y(k,218)
         mat(k,1229) = rxt(k,316)*y(k,218)
         mat(k,689) = rxt(k,296)*y(k,218)
         mat(k,617) = rxt(k,297)*y(k,218)
         mat(k,460) = rxt(k,265)*y(k,218)
         mat(k,1443) = rxt(k,266)*y(k,218)
         mat(k,1458) = rxt(k,138)*y(k,204)
         mat(k,1413) = rxt(k,143)*y(k,218)
         mat(k,610) = rxt(k,144)*y(k,218)
         mat(k,812) = rxt(k,226)*y(k,218)
         mat(k,304) = rxt(k,250)*y(k,218)
         mat(k,2144) = (rxt(k,531)+rxt(k,536))*y(k,90) + (rxt(k,524)+rxt(k,530) &
                       +rxt(k,535))*y(k,91) + rxt(k,197)*y(k,218)
         mat(k,935) = rxt(k,268)*y(k,218)
         mat(k,1476) = rxt(k,173)*y(k,218)
         mat(k,427) = rxt(k,151)*y(k,218)
         mat(k,786) = (rxt(k,531)+rxt(k,536))*y(k,84)
         mat(k,820) = (rxt(k,524)+rxt(k,530)+rxt(k,535))*y(k,84) + rxt(k,200)*y(k,218)
         mat(k,1220) = .500_r8*rxt(k,340)*y(k,218)
         mat(k,99) = rxt(k,480)*y(k,218)
         mat(k,548) = rxt(k,321)*y(k,218)
         mat(k,403) = rxt(k,325)*y(k,218)
         mat(k,2121) = rxt(k,138)*y(k,75) + rxt(k,145)*y(k,218)
         mat(k,1704) = rxt(k,288)*y(k,27) + rxt(k,313)*y(k,29) + rxt(k,314)*y(k,30) &
                      + rxt(k,240)*y(k,40) + rxt(k,259)*y(k,41) + rxt(k,242)*y(k,42) &
                      + rxt(k,243)*y(k,43) + rxt(k,290)*y(k,44) + rxt(k,245)*y(k,45) &
                      + rxt(k,327)*y(k,47) + rxt(k,316)*y(k,48) + rxt(k,296)*y(k,49) &
                      + rxt(k,297)*y(k,50) + rxt(k,265)*y(k,52) + rxt(k,266)*y(k,53) &
                      + rxt(k,143)*y(k,76) + rxt(k,144)*y(k,78) + rxt(k,226)*y(k,80) &
                      + rxt(k,250)*y(k,83) + rxt(k,197)*y(k,84) + rxt(k,268)*y(k,86) &
                      + rxt(k,173)*y(k,88) + rxt(k,151)*y(k,89) + rxt(k,200)*y(k,91) &
                      + .500_r8*rxt(k,340)*y(k,104) + rxt(k,480)*y(k,119) + rxt(k,321) &
                      *y(k,146) + rxt(k,325)*y(k,147) + rxt(k,145)*y(k,204) &
                      + 2.000_r8*rxt(k,148)*y(k,218)
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
         mat(k, 36) = lmat(k, 36)
         mat(k, 37) = lmat(k, 37)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = lmat(k, 42)
         mat(k, 43) = lmat(k, 43)
         mat(k, 49) = mat(k, 49) + lmat(k, 49)
         mat(k, 55) = mat(k, 55) + lmat(k, 55)
         mat(k, 61) = mat(k, 61) + lmat(k, 61)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 81) = mat(k, 81) + lmat(k, 81)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 93) = mat(k, 93) + lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 100) = mat(k, 100) + lmat(k, 100)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 103) = mat(k, 103) + lmat(k, 103)
         mat(k, 106) = mat(k, 106) + lmat(k, 106)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 110) = mat(k, 110) + lmat(k, 110)
         mat(k, 114) = mat(k, 114) + lmat(k, 114)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 121) = mat(k, 121) + lmat(k, 121)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 125) = mat(k, 125) + lmat(k, 125)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 130) = lmat(k, 130)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = mat(k, 136) + lmat(k, 136)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = lmat(k, 141)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 169) = lmat(k, 169)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 172) = lmat(k, 172)
         mat(k, 173) = mat(k, 173) + lmat(k, 173)
         mat(k, 178) = lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 188) = mat(k, 188) + lmat(k, 188)
         mat(k, 195) = mat(k, 195) + lmat(k, 195)
         mat(k, 200) = lmat(k, 200)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 206) = mat(k, 206) + lmat(k, 206)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 217) = mat(k, 217) + lmat(k, 217)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 223) = lmat(k, 223)
         mat(k, 225) = mat(k, 225) + lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 230) = lmat(k, 230)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = lmat(k, 260)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 266) = mat(k, 266) + lmat(k, 266)
         mat(k, 267) = lmat(k, 267)
         mat(k, 269) = mat(k, 269) + lmat(k, 269)
         mat(k, 270) = lmat(k, 270)
         mat(k, 271) = lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
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
         mat(k, 283) = mat(k, 283) + lmat(k, 283)
         mat(k, 289) = lmat(k, 289)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = mat(k, 292) + lmat(k, 292)
         mat(k, 298) = mat(k, 298) + lmat(k, 298)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 305) = mat(k, 305) + lmat(k, 305)
         mat(k, 311) = mat(k, 311) + lmat(k, 311)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 318) = lmat(k, 318)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = mat(k, 320) + lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 330) = lmat(k, 330)
         mat(k, 331) = mat(k, 331) + lmat(k, 331)
         mat(k, 333) = lmat(k, 333)
         mat(k, 334) = mat(k, 334) + lmat(k, 334)
         mat(k, 335) = lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = mat(k, 337) + lmat(k, 337)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 344) = mat(k, 344) + lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 350) = mat(k, 350) + lmat(k, 350)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 359) = lmat(k, 359)
         mat(k, 361) = mat(k, 361) + lmat(k, 361)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 377) = lmat(k, 377)
         mat(k, 380) = lmat(k, 380)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = mat(k, 386) + lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 393) = lmat(k, 393)
         mat(k, 395) = lmat(k, 395)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 398) = mat(k, 398) + lmat(k, 398)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = lmat(k, 401)
         mat(k, 402) = mat(k, 402) + lmat(k, 402)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 405) = lmat(k, 405)
         mat(k, 408) = lmat(k, 408)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = mat(k, 413) + lmat(k, 413)
         mat(k, 414) = lmat(k, 414)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 419) = lmat(k, 419)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 421) = lmat(k, 421)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 426) = lmat(k, 426)
         mat(k, 428) = mat(k, 428) + lmat(k, 428)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = lmat(k, 431)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 433) = lmat(k, 433)
         mat(k, 436) = mat(k, 436) + lmat(k, 436)
         mat(k, 442) = mat(k, 442) + lmat(k, 442)
         mat(k, 444) = lmat(k, 444)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 456) = lmat(k, 456)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 461) = lmat(k, 461)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = mat(k, 464) + lmat(k, 464)
         mat(k, 468) = mat(k, 468) + lmat(k, 468)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 477) = lmat(k, 477)
         mat(k, 478) = lmat(k, 478)
         mat(k, 479) = lmat(k, 479)
         mat(k, 480) = mat(k, 480) + lmat(k, 480)
         mat(k, 481) = lmat(k, 481)
         mat(k, 484) = mat(k, 484) + lmat(k, 484)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 498) = mat(k, 498) + lmat(k, 498)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 507) = mat(k, 507) + lmat(k, 507)
         mat(k, 509) = lmat(k, 509)
         mat(k, 510) = lmat(k, 510)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 516) = lmat(k, 516)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 538) = lmat(k, 538)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 543) = lmat(k, 543)
         mat(k, 545) = lmat(k, 545)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 547) = lmat(k, 547)
         mat(k, 549) = mat(k, 549) + lmat(k, 549)
         mat(k, 557) = mat(k, 557) + lmat(k, 557)
         mat(k, 558) = lmat(k, 558)
         mat(k, 559) = lmat(k, 559)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = mat(k, 561) + lmat(k, 561)
         mat(k, 562) = lmat(k, 562)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 565) = mat(k, 565) + lmat(k, 565)
         mat(k, 573) = mat(k, 573) + lmat(k, 573)
         mat(k, 576) = lmat(k, 576)
         mat(k, 580) = lmat(k, 580)
         mat(k, 582) = mat(k, 582) + lmat(k, 582)
         mat(k, 586) = lmat(k, 586)
         mat(k, 591) = lmat(k, 591)
         mat(k, 592) = lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = lmat(k, 594)
         mat(k, 595) = mat(k, 595) + lmat(k, 595)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 599) = lmat(k, 599)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 605) = mat(k, 605) + lmat(k, 605)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 612) = mat(k, 612) + lmat(k, 612)
         mat(k, 615) = mat(k, 615) + lmat(k, 615)
         mat(k, 616) = lmat(k, 616)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 619) = mat(k, 619) + lmat(k, 619)
         mat(k, 620) = lmat(k, 620)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 623) = lmat(k, 623)
         mat(k, 626) = mat(k, 626) + lmat(k, 626)
         mat(k, 632) = lmat(k, 632)
         mat(k, 633) = mat(k, 633) + lmat(k, 633)
         mat(k, 637) = lmat(k, 637)
         mat(k, 638) = lmat(k, 638)
         mat(k, 640) = lmat(k, 640)
         mat(k, 641) = mat(k, 641) + lmat(k, 641)
         mat(k, 642) = lmat(k, 642)
         mat(k, 643) = mat(k, 643) + lmat(k, 643)
         mat(k, 644) = lmat(k, 644)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = lmat(k, 650)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = lmat(k, 655)
         mat(k, 656) = mat(k, 656) + lmat(k, 656)
         mat(k, 661) = lmat(k, 661)
         mat(k, 663) = lmat(k, 663)
         mat(k, 664) = mat(k, 664) + lmat(k, 664)
         mat(k, 665) = lmat(k, 665)
         mat(k, 666) = lmat(k, 666)
         mat(k, 667) = mat(k, 667) + lmat(k, 667)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 671) = mat(k, 671) + lmat(k, 671)
         mat(k, 673) = mat(k, 673) + lmat(k, 673)
         mat(k, 675) = mat(k, 675) + lmat(k, 675)
         mat(k, 676) = lmat(k, 676)
         mat(k, 678) = mat(k, 678) + lmat(k, 678)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 690) = lmat(k, 690)
         mat(k, 691) = lmat(k, 691)
         mat(k, 692) = lmat(k, 692)
         mat(k, 693) = lmat(k, 693)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 699) = lmat(k, 699)
         mat(k, 701) = lmat(k, 701)
         mat(k, 703) = mat(k, 703) + lmat(k, 703)
         mat(k, 704) = lmat(k, 704)
         mat(k, 707) = mat(k, 707) + lmat(k, 707)
         mat(k, 714) = mat(k, 714) + lmat(k, 714)
         mat(k, 725) = mat(k, 725) + lmat(k, 725)
         mat(k, 736) = mat(k, 736) + lmat(k, 736)
         mat(k, 738) = lmat(k, 738)
         mat(k, 739) = lmat(k, 739)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = mat(k, 741) + lmat(k, 741)
         mat(k, 742) = lmat(k, 742)
         mat(k, 749) = mat(k, 749) + lmat(k, 749)
         mat(k, 760) = mat(k, 760) + lmat(k, 760)
         mat(k, 769) = mat(k, 769) + lmat(k, 769)
         mat(k, 779) = mat(k, 779) + lmat(k, 779)
         mat(k, 780) = lmat(k, 780)
         mat(k, 781) = mat(k, 781) + lmat(k, 781)
         mat(k, 788) = mat(k, 788) + lmat(k, 788)
         mat(k, 798) = mat(k, 798) + lmat(k, 798)
         mat(k, 802) = lmat(k, 802)
         mat(k, 803) = lmat(k, 803)
         mat(k, 804) = lmat(k, 804)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 807) = mat(k, 807) + lmat(k, 807)
         mat(k, 814) = mat(k, 814) + lmat(k, 814)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 817) = mat(k, 817) + lmat(k, 817)
         mat(k, 824) = mat(k, 824) + lmat(k, 824)
         mat(k, 835) = mat(k, 835) + lmat(k, 835)
         mat(k, 852) = mat(k, 852) + lmat(k, 852)
         mat(k, 853) = lmat(k, 853)
         mat(k, 856) = lmat(k, 856)
         mat(k, 858) = mat(k, 858) + lmat(k, 858)
         mat(k, 860) = lmat(k, 860)
         mat(k, 861) = lmat(k, 861)
         mat(k, 863) = mat(k, 863) + lmat(k, 863)
         mat(k, 865) = mat(k, 865) + lmat(k, 865)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 888) = mat(k, 888) + lmat(k, 888)
         mat(k, 899) = mat(k, 899) + lmat(k, 899)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 901) = mat(k, 901) + lmat(k, 901)
         mat(k, 903) = mat(k, 903) + lmat(k, 903)
         mat(k, 904) = mat(k, 904) + lmat(k, 904)
         mat(k, 905) = lmat(k, 905)
         mat(k, 906) = mat(k, 906) + lmat(k, 906)
         mat(k, 913) = mat(k, 913) + lmat(k, 913)
         mat(k, 923) = mat(k, 923) + lmat(k, 923)
         mat(k, 925) = lmat(k, 925)
         mat(k, 927) = lmat(k, 927)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 949) = mat(k, 949) + lmat(k, 949)
         mat(k, 970) = mat(k, 970) + lmat(k, 970)
         mat(k, 972) = lmat(k, 972)
         mat(k, 974) = lmat(k, 974)
         mat(k, 977) = lmat(k, 977)
         mat(k, 978) = lmat(k, 978)
         mat(k, 980) = mat(k, 980) + lmat(k, 980)
         mat(k, 981) = mat(k, 981) + lmat(k, 981)
         mat(k, 983) = mat(k, 983) + lmat(k, 983)
         mat(k,1000) = mat(k,1000) + lmat(k,1000)
         mat(k,1025) = mat(k,1025) + lmat(k,1025)
         mat(k,1046) = mat(k,1046) + lmat(k,1046)
         mat(k,1058) = mat(k,1058) + lmat(k,1058)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1061) = mat(k,1061) + lmat(k,1061)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1067) = mat(k,1067) + lmat(k,1067)
         mat(k,1070) = lmat(k,1070)
         mat(k,1074) = mat(k,1074) + lmat(k,1074)
         mat(k,1078) = lmat(k,1078)
         mat(k,1083) = lmat(k,1083)
         mat(k,1084) = mat(k,1084) + lmat(k,1084)
         mat(k,1086) = mat(k,1086) + lmat(k,1086)
         mat(k,1087) = lmat(k,1087)
         mat(k,1092) = lmat(k,1092)
         mat(k,1093) = lmat(k,1093)
         mat(k,1097) = mat(k,1097) + lmat(k,1097)
         mat(k,1098) = lmat(k,1098)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1130) = mat(k,1130) + lmat(k,1130)
         mat(k,1141) = mat(k,1141) + lmat(k,1141)
         mat(k,1143) = lmat(k,1143)
         mat(k,1144) = lmat(k,1144)
         mat(k,1146) = mat(k,1146) + lmat(k,1146)
         mat(k,1148) = lmat(k,1148)
         mat(k,1149) = lmat(k,1149)
         mat(k,1150) = lmat(k,1150)
         mat(k,1151) = lmat(k,1151)
         mat(k,1153) = lmat(k,1153)
         mat(k,1154) = mat(k,1154) + lmat(k,1154)
         mat(k,1156) = lmat(k,1156)
         mat(k,1157) = lmat(k,1157)
         mat(k,1158) = lmat(k,1158)
         mat(k,1163) = lmat(k,1163)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1174) = mat(k,1174) + lmat(k,1174)
         mat(k,1194) = mat(k,1194) + lmat(k,1194)
         mat(k,1209) = mat(k,1209) + lmat(k,1209)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1215) = mat(k,1215) + lmat(k,1215)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1221) = mat(k,1221) + lmat(k,1221)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1223) = mat(k,1223) + lmat(k,1223)
         mat(k,1227) = lmat(k,1227)
         mat(k,1231) = lmat(k,1231)
         mat(k,1232) = mat(k,1232) + lmat(k,1232)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1243) = lmat(k,1243)
         mat(k,1257) = mat(k,1257) + lmat(k,1257)
         mat(k,1273) = lmat(k,1273)
         mat(k,1289) = mat(k,1289) + lmat(k,1289)
         mat(k,1301) = mat(k,1301) + lmat(k,1301)
         mat(k,1312) = mat(k,1312) + lmat(k,1312)
         mat(k,1327) = lmat(k,1327)
         mat(k,1329) = mat(k,1329) + lmat(k,1329)
         mat(k,1333) = mat(k,1333) + lmat(k,1333)
         mat(k,1335) = mat(k,1335) + lmat(k,1335)
         mat(k,1341) = lmat(k,1341)
         mat(k,1355) = mat(k,1355) + lmat(k,1355)
         mat(k,1387) = mat(k,1387) + lmat(k,1387)
         mat(k,1402) = mat(k,1402) + lmat(k,1402)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1427) = lmat(k,1427)
         mat(k,1429) = lmat(k,1429)
         mat(k,1430) = mat(k,1430) + lmat(k,1430)
         mat(k,1431) = mat(k,1431) + lmat(k,1431)
         mat(k,1432) = mat(k,1432) + lmat(k,1432)
         mat(k,1434) = mat(k,1434) + lmat(k,1434)
         mat(k,1436) = mat(k,1436) + lmat(k,1436)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1442) = lmat(k,1442)
         mat(k,1443) = mat(k,1443) + lmat(k,1443)
         mat(k,1446) = mat(k,1446) + lmat(k,1446)
         mat(k,1454) = mat(k,1454) + lmat(k,1454)
         mat(k,1464) = mat(k,1464) + lmat(k,1464)
         mat(k,1467) = mat(k,1467) + lmat(k,1467)
         mat(k,1470) = lmat(k,1470)
         mat(k,1480) = mat(k,1480) + lmat(k,1480)
         mat(k,1481) = lmat(k,1481)
         mat(k,1484) = mat(k,1484) + lmat(k,1484)
         mat(k,1486) = mat(k,1486) + lmat(k,1486)
         mat(k,1527) = mat(k,1527) + lmat(k,1527)
         mat(k,1538) = mat(k,1538) + lmat(k,1538)
         mat(k,1692) = mat(k,1692) + lmat(k,1692)
         mat(k,1719) = mat(k,1719) + lmat(k,1719)
         mat(k,1724) = mat(k,1724) + lmat(k,1724)
         mat(k,1728) = mat(k,1728) + lmat(k,1728)
         mat(k,1772) = mat(k,1772) + lmat(k,1772)
         mat(k,1777) = mat(k,1777) + lmat(k,1777)
         mat(k,1779) = mat(k,1779) + lmat(k,1779)
         mat(k,1780) = mat(k,1780) + lmat(k,1780)
         mat(k,1785) = mat(k,1785) + lmat(k,1785)
         mat(k,1830) = mat(k,1830) + lmat(k,1830)
         mat(k,1864) = mat(k,1864) + lmat(k,1864)
         mat(k,1923) = mat(k,1923) + lmat(k,1923)
         mat(k,1929) = mat(k,1929) + lmat(k,1929)
         mat(k,1960) = mat(k,1960) + lmat(k,1960)
         mat(k,1963) = mat(k,1963) + lmat(k,1963)
         mat(k,1967) = mat(k,1967) + lmat(k,1967)
         mat(k,1968) = mat(k,1968) + lmat(k,1968)
         mat(k,1973) = mat(k,1973) + lmat(k,1973)
         mat(k,2008) = mat(k,2008) + lmat(k,2008)
         mat(k,2116) = mat(k,2116) + lmat(k,2116)
         mat(k,2121) = mat(k,2121) + lmat(k,2121)
         mat(k,2128) = mat(k,2128) + lmat(k,2128)
         mat(k,2138) = mat(k,2138) + lmat(k,2138)
         mat(k,2140) = mat(k,2140) + lmat(k,2140)
         mat(k,2151) = mat(k,2151) + lmat(k,2151)
         mat(k,2165) = mat(k,2165) + lmat(k,2165)
         mat(k,2166) = mat(k,2166) + lmat(k,2166)
         mat(k,2197) = mat(k,2197) + lmat(k,2197)
         mat(k,2198) = mat(k,2198) + lmat(k,2198)
         mat(k,2247) = mat(k,2247) + lmat(k,2247)
         mat(k,2258) = mat(k,2258) + lmat(k,2258)
         mat(k,2259) = mat(k,2259) + lmat(k,2259)
         mat(k,2267) = lmat(k,2267)
         mat(k,2270) = lmat(k,2270)
         mat(k,2273) = mat(k,2273) + lmat(k,2273)
         mat(k,2274) = mat(k,2274) + lmat(k,2274)
         mat(k,2284) = lmat(k,2284)
         mat(k,2286) = mat(k,2286) + lmat(k,2286)
         mat(k, 213) = 0._r8
         mat(k, 214) = 0._r8
         mat(k, 253) = 0._r8
         mat(k, 300) = 0._r8
         mat(k, 339) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 438) = 0._r8
         mat(k, 451) = 0._r8
         mat(k, 485) = 0._r8
         mat(k, 487) = 0._r8
         mat(k, 494) = 0._r8
         mat(k, 503) = 0._r8
         mat(k, 627) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 657) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 668) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 695) = 0._r8
         mat(k, 697) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 702) = 0._r8
         mat(k, 708) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 724) = 0._r8
         mat(k, 726) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 729) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 748) = 0._r8
         mat(k, 750) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 771) = 0._r8
         mat(k, 773) = 0._r8
         mat(k, 777) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 792) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 827) = 0._r8
         mat(k, 857) = 0._r8
         mat(k, 886) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 893) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 909) = 0._r8
         mat(k, 916) = 0._r8
         mat(k, 918) = 0._r8
         mat(k, 931) = 0._r8
         mat(k, 934) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 963) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 969) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1007) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1017) = 0._r8
         mat(k,1019) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1047) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1052) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1072) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1079) = 0._r8
         mat(k,1081) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1102) = 0._r8
         mat(k,1112) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1116) = 0._r8
         mat(k,1121) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1152) = 0._r8
         mat(k,1155) = 0._r8
         mat(k,1160) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1178) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1187) = 0._r8
         mat(k,1189) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1228) = 0._r8
         mat(k,1241) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1250) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1260) = 0._r8
         mat(k,1262) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1270) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1283) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1294) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1302) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1315) = 0._r8
         mat(k,1317) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1334) = 0._r8
         mat(k,1336) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1340) = 0._r8
         mat(k,1342) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1367) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1397) = 0._r8
         mat(k,1399) = 0._r8
         mat(k,1403) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1408) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1420) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1437) = 0._r8
         mat(k,1445) = 0._r8
         mat(k,1447) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1451) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1465) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1469) = 0._r8
         mat(k,1471) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1483) = 0._r8
         mat(k,1487) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1530) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1606) = 0._r8
         mat(k,1623) = 0._r8
         mat(k,1637) = 0._r8
         mat(k,1641) = 0._r8
         mat(k,1652) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1675) = 0._r8
         mat(k,1691) = 0._r8
         mat(k,1714) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1717) = 0._r8
         mat(k,1720) = 0._r8
         mat(k,1729) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1742) = 0._r8
         mat(k,1749) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1769) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1774) = 0._r8
         mat(k,1776) = 0._r8
         mat(k,1778) = 0._r8
         mat(k,1781) = 0._r8
         mat(k,1783) = 0._r8
         mat(k,1784) = 0._r8
         mat(k,1786) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1826) = 0._r8
         mat(k,1827) = 0._r8
         mat(k,1829) = 0._r8
         mat(k,1835) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1885) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1918) = 0._r8
         mat(k,1927) = 0._r8
         mat(k,1931) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1961) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1966) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1971) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2001) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2007) = 0._r8
         mat(k,2011) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2030) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2080) = 0._r8
         mat(k,2082) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2091) = 0._r8
         mat(k,2094) = 0._r8
         mat(k,2099) = 0._r8
         mat(k,2104) = 0._r8
         mat(k,2106) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2127) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2134) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2139) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2153) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2158) = 0._r8
         mat(k,2159) = 0._r8
         mat(k,2164) = 0._r8
         mat(k,2167) = 0._r8
         mat(k,2168) = 0._r8
         mat(k,2170) = 0._r8
         mat(k,2176) = 0._r8
         mat(k,2182) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2190) = 0._r8
         mat(k,2199) = 0._r8
         mat(k,2212) = 0._r8
         mat(k,2216) = 0._r8
         mat(k,2221) = 0._r8
         mat(k,2224) = 0._r8
         mat(k,2227) = 0._r8
         mat(k,2228) = 0._r8
         mat(k,2231) = 0._r8
         mat(k,2232) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2237) = 0._r8
         mat(k,2238) = 0._r8
         mat(k,2240) = 0._r8
         mat(k,2245) = 0._r8
         mat(k,2256) = 0._r8
         mat(k,2260) = 0._r8
         mat(k,2264) = 0._r8
         mat(k,2266) = 0._r8
         mat(k,2268) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2271) = 0._r8
         mat(k,2272) = 0._r8
         mat(k,2275) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2281) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2283) = 0._r8
         mat(k,2285) = 0._r8
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
         mat(k, 49) = mat(k, 49) - dti(k)
         mat(k, 55) = mat(k, 55) - dti(k)
         mat(k, 61) = mat(k, 61) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 75) = mat(k, 75) - dti(k)
         mat(k, 81) = mat(k, 81) - dti(k)
         mat(k, 87) = mat(k, 87) - dti(k)
         mat(k, 93) = mat(k, 93) - dti(k)
         mat(k, 94) = mat(k, 94) - dti(k)
         mat(k, 97) = mat(k, 97) - dti(k)
         mat(k, 100) = mat(k, 100) - dti(k)
         mat(k, 103) = mat(k, 103) - dti(k)
         mat(k, 106) = mat(k, 106) - dti(k)
         mat(k, 110) = mat(k, 110) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 118) = mat(k, 118) - dti(k)
         mat(k, 122) = mat(k, 122) - dti(k)
         mat(k, 126) = mat(k, 126) - dti(k)
         mat(k, 130) = mat(k, 130) - dti(k)
         mat(k, 133) = mat(k, 133) - dti(k)
         mat(k, 136) = mat(k, 136) - dti(k)
         mat(k, 139) = mat(k, 139) - dti(k)
         mat(k, 142) = mat(k, 142) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 163) = mat(k, 163) - dti(k)
         mat(k, 169) = mat(k, 169) - dti(k)
         mat(k, 173) = mat(k, 173) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 184) = mat(k, 184) - dti(k)
         mat(k, 188) = mat(k, 188) - dti(k)
         mat(k, 195) = mat(k, 195) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 209) = mat(k, 209) - dti(k)
         mat(k, 217) = mat(k, 217) - dti(k)
         mat(k, 222) = mat(k, 222) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 230) = mat(k, 230) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 236) = mat(k, 236) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 242) = mat(k, 242) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 252) = mat(k, 252) - dti(k)
         mat(k, 258) = mat(k, 258) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 266) = mat(k, 266) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 274) = mat(k, 274) - dti(k)
         mat(k, 280) = mat(k, 280) - dti(k)
         mat(k, 283) = mat(k, 283) - dti(k)
         mat(k, 289) = mat(k, 289) - dti(k)
         mat(k, 292) = mat(k, 292) - dti(k)
         mat(k, 299) = mat(k, 299) - dti(k)
         mat(k, 305) = mat(k, 305) - dti(k)
         mat(k, 311) = mat(k, 311) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 326) = mat(k, 326) - dti(k)
         mat(k, 331) = mat(k, 331) - dti(k)
         mat(k, 337) = mat(k, 337) - dti(k)
         mat(k, 342) = mat(k, 342) - dti(k)
         mat(k, 347) = mat(k, 347) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 358) = mat(k, 358) - dti(k)
         mat(k, 366) = mat(k, 366) - dti(k)
         mat(k, 374) = mat(k, 374) - dti(k)
         mat(k, 380) = mat(k, 380) - dti(k)
         mat(k, 386) = mat(k, 386) - dti(k)
         mat(k, 392) = mat(k, 392) - dti(k)
         mat(k, 398) = mat(k, 398) - dti(k)
         mat(k, 404) = mat(k, 404) - dti(k)
         mat(k, 410) = mat(k, 410) - dti(k)
         mat(k, 416) = mat(k, 416) - dti(k)
         mat(k, 422) = mat(k, 422) - dti(k)
         mat(k, 428) = mat(k, 428) - dti(k)
         mat(k, 436) = mat(k, 436) - dti(k)
         mat(k, 442) = mat(k, 442) - dti(k)
         mat(k, 449) = mat(k, 449) - dti(k)
         mat(k, 455) = mat(k, 455) - dti(k)
         mat(k, 461) = mat(k, 461) - dti(k)
         mat(k, 464) = mat(k, 464) - dti(k)
         mat(k, 468) = mat(k, 468) - dti(k)
         mat(k, 475) = mat(k, 475) - dti(k)
         mat(k, 484) = mat(k, 484) - dti(k)
         mat(k, 492) = mat(k, 492) - dti(k)
         mat(k, 500) = mat(k, 500) - dti(k)
         mat(k, 507) = mat(k, 507) - dti(k)
         mat(k, 512) = mat(k, 512) - dti(k)
         mat(k, 519) = mat(k, 519) - dti(k)
         mat(k, 525) = mat(k, 525) - dti(k)
         mat(k, 533) = mat(k, 533) - dti(k)
         mat(k, 541) = mat(k, 541) - dti(k)
         mat(k, 549) = mat(k, 549) - dti(k)
         mat(k, 557) = mat(k, 557) - dti(k)
         mat(k, 565) = mat(k, 565) - dti(k)
         mat(k, 573) = mat(k, 573) - dti(k)
         mat(k, 582) = mat(k, 582) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 595) = mat(k, 595) - dti(k)
         mat(k, 604) = mat(k, 604) - dti(k)
         mat(k, 611) = mat(k, 611) - dti(k)
         mat(k, 618) = mat(k, 618) - dti(k)
         mat(k, 626) = mat(k, 626) - dti(k)
         mat(k, 633) = mat(k, 633) - dti(k)
         mat(k, 643) = mat(k, 643) - dti(k)
         mat(k, 656) = mat(k, 656) - dti(k)
         mat(k, 667) = mat(k, 667) - dti(k)
         mat(k, 678) = mat(k, 678) - dti(k)
         mat(k, 685) = mat(k, 685) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 707) = mat(k, 707) - dti(k)
         mat(k, 714) = mat(k, 714) - dti(k)
         mat(k, 725) = mat(k, 725) - dti(k)
         mat(k, 736) = mat(k, 736) - dti(k)
         mat(k, 749) = mat(k, 749) - dti(k)
         mat(k, 760) = mat(k, 760) - dti(k)
         mat(k, 769) = mat(k, 769) - dti(k)
         mat(k, 779) = mat(k, 779) - dti(k)
         mat(k, 788) = mat(k, 788) - dti(k)
         mat(k, 798) = mat(k, 798) - dti(k)
         mat(k, 802) = mat(k, 802) - dti(k)
         mat(k, 805) = mat(k, 805) - dti(k)
         mat(k, 814) = mat(k, 814) - dti(k)
         mat(k, 824) = mat(k, 824) - dti(k)
         mat(k, 835) = mat(k, 835) - dti(k)
         mat(k, 852) = mat(k, 852) - dti(k)
         mat(k, 858) = mat(k, 858) - dti(k)
         mat(k, 865) = mat(k, 865) - dti(k)
         mat(k, 874) = mat(k, 874) - dti(k)
         mat(k, 888) = mat(k, 888) - dti(k)
         mat(k, 900) = mat(k, 900) - dti(k)
         mat(k, 913) = mat(k, 913) - dti(k)
         mat(k, 923) = mat(k, 923) - dti(k)
         mat(k, 930) = mat(k, 930) - dti(k)
         mat(k, 949) = mat(k, 949) - dti(k)
         mat(k, 970) = mat(k, 970) - dti(k)
         mat(k, 980) = mat(k, 980) - dti(k)
         mat(k,1000) = mat(k,1000) - dti(k)
         mat(k,1025) = mat(k,1025) - dti(k)
         mat(k,1046) = mat(k,1046) - dti(k)
         mat(k,1060) = mat(k,1060) - dti(k)
         mat(k,1074) = mat(k,1074) - dti(k)
         mat(k,1086) = mat(k,1086) - dti(k)
         mat(k,1097) = mat(k,1097) - dti(k)
         mat(k,1110) = mat(k,1110) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1130) = mat(k,1130) - dti(k)
         mat(k,1141) = mat(k,1141) - dti(k)
         mat(k,1154) = mat(k,1154) - dti(k)
         mat(k,1174) = mat(k,1174) - dti(k)
         mat(k,1194) = mat(k,1194) - dti(k)
         mat(k,1210) = mat(k,1210) - dti(k)
         mat(k,1222) = mat(k,1222) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1257) = mat(k,1257) - dti(k)
         mat(k,1289) = mat(k,1289) - dti(k)
         mat(k,1312) = mat(k,1312) - dti(k)
         mat(k,1333) = mat(k,1333) - dti(k)
         mat(k,1355) = mat(k,1355) - dti(k)
         mat(k,1387) = mat(k,1387) - dti(k)
         mat(k,1402) = mat(k,1402) - dti(k)
         mat(k,1416) = mat(k,1416) - dti(k)
         mat(k,1431) = mat(k,1431) - dti(k)
         mat(k,1446) = mat(k,1446) - dti(k)
         mat(k,1464) = mat(k,1464) - dti(k)
         mat(k,1486) = mat(k,1486) - dti(k)
         mat(k,1527) = mat(k,1527) - dti(k)
         mat(k,1692) = mat(k,1692) - dti(k)
         mat(k,1719) = mat(k,1719) - dti(k)
         mat(k,1777) = mat(k,1777) - dti(k)
         mat(k,1830) = mat(k,1830) - dti(k)
         mat(k,1923) = mat(k,1923) - dti(k)
         mat(k,1968) = mat(k,1968) - dti(k)
         mat(k,2008) = mat(k,2008) - dti(k)
         mat(k,2116) = mat(k,2116) - dti(k)
         mat(k,2140) = mat(k,2140) - dti(k)
         mat(k,2165) = mat(k,2165) - dti(k)
         mat(k,2197) = mat(k,2197) - dti(k)
         mat(k,2259) = mat(k,2259) - dti(k)
         mat(k,2286) = mat(k,2286) - dti(k)
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
