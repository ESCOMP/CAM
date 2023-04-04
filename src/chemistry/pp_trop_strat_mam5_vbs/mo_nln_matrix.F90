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
         mat(k,634) = -(rxt(k,356)*y(k,219))
         mat(k,1630) = -rxt(k,356)*y(k,1)
         mat(k,1875) = rxt(k,359)*y(k,191)
         mat(k,885) = rxt(k,359)*y(k,124)
         mat(k,668) = -(rxt(k,360)*y(k,219))
         mat(k,1633) = -rxt(k,360)*y(k,2)
         mat(k,886) = rxt(k,357)*y(k,205)
         mat(k,2060) = rxt(k,357)*y(k,191)
         mat(k,1001) = -(rxt(k,439)*y(k,126) + rxt(k,440)*y(k,135) + rxt(k,441) &
                      *y(k,219))
         mat(k,1747) = -rxt(k,439)*y(k,6)
         mat(k,2223) = -rxt(k,440)*y(k,6)
         mat(k,1663) = -rxt(k,441)*y(k,6)
         mat(k,164) = -(rxt(k,398)*y(k,219))
         mat(k,1561) = -rxt(k,398)*y(k,7)
         mat(k,417) = -(rxt(k,401)*y(k,219))
         mat(k,1601) = -rxt(k,401)*y(k,8)
         mat(k,483) = rxt(k,399)*y(k,205)
         mat(k,2041) = rxt(k,399)*y(k,193)
         mat(k,165) = .120_r8*rxt(k,398)*y(k,219)
         mat(k,1562) = .120_r8*rxt(k,398)*y(k,7)
         mat(k,993) = .100_r8*rxt(k,440)*y(k,135)
         mat(k,944) = .100_r8*rxt(k,443)*y(k,135)
         mat(k,2207) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,110)
         mat(k,1862) = .500_r8*rxt(k,400)*y(k,193) + .200_r8*rxt(k,427)*y(k,225) &
                      + .060_r8*rxt(k,433)*y(k,228)
         mat(k,484) = .500_r8*rxt(k,400)*y(k,124)
         mat(k,722) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,746) = .060_r8*rxt(k,433)*y(k,124)
         mat(k,1856) = .200_r8*rxt(k,427)*y(k,225) + .200_r8*rxt(k,433)*y(k,228)
         mat(k,721) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,744) = .200_r8*rxt(k,433)*y(k,124)
         mat(k,1872) = .200_r8*rxt(k,427)*y(k,225) + .150_r8*rxt(k,433)*y(k,228)
         mat(k,723) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,747) = .150_r8*rxt(k,433)*y(k,124)
         mat(k,1857) = .210_r8*rxt(k,433)*y(k,228)
         mat(k,745) = .210_r8*rxt(k,433)*y(k,124)
         mat(k,228) = -(rxt(k,361)*y(k,219))
         mat(k,1572) = -rxt(k,361)*y(k,15)
         mat(k,992) = .050_r8*rxt(k,440)*y(k,135)
         mat(k,943) = .050_r8*rxt(k,443)*y(k,135)
         mat(k,2206) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,110)
         mat(k,351) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,219))
         mat(k,1737) = -rxt(k,327)*y(k,16)
         mat(k,1592) = -rxt(k,328)*y(k,16)
         mat(k,1417) = -(rxt(k,211)*y(k,42) + rxt(k,212)*y(k,205) + rxt(k,213) &
                      *y(k,135))
         mat(k,1483) = -rxt(k,211)*y(k,17)
         mat(k,2104) = -rxt(k,212)*y(k,17)
         mat(k,2243) = -rxt(k,213)*y(k,17)
         mat(k,2152) = 4.000_r8*rxt(k,214)*y(k,19) + (rxt(k,215)+rxt(k,216))*y(k,59) &
                      + rxt(k,219)*y(k,124) + rxt(k,222)*y(k,134) + rxt(k,469) &
                      *y(k,151) + rxt(k,223)*y(k,219)
         mat(k,145) = rxt(k,201)*y(k,218)
         mat(k,151) = rxt(k,227)*y(k,218)
         mat(k,470) = 2.000_r8*rxt(k,238)*y(k,56) + 2.000_r8*rxt(k,250)*y(k,218) &
                      + 2.000_r8*rxt(k,239)*y(k,219)
         mat(k,597) = rxt(k,240)*y(k,56) + rxt(k,251)*y(k,218) + rxt(k,241)*y(k,219)
         mat(k,388) = 3.000_r8*rxt(k,245)*y(k,56) + 3.000_r8*rxt(k,228)*y(k,218) &
                      + 3.000_r8*rxt(k,246)*y(k,219)
         mat(k,1997) = 2.000_r8*rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43) &
                      + 3.000_r8*rxt(k,245)*y(k,55)
         mat(k,1714) = (rxt(k,215)+rxt(k,216))*y(k,19)
         mat(k,109) = 2.000_r8*rxt(k,229)*y(k,218)
         mat(k,807) = rxt(k,224)*y(k,134) + rxt(k,230)*y(k,218) + rxt(k,225)*y(k,219)
         mat(k,1914) = rxt(k,219)*y(k,19)
         mat(k,2182) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81)
         mat(k,1235) = rxt(k,469)*y(k,19)
         mat(k,1523) = rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35) + 2.000_r8*rxt(k,250) &
                      *y(k,41) + rxt(k,251)*y(k,43) + 3.000_r8*rxt(k,228)*y(k,55) &
                      + 2.000_r8*rxt(k,229)*y(k,78) + rxt(k,230)*y(k,81)
         mat(k,1687) = rxt(k,223)*y(k,19) + 2.000_r8*rxt(k,239)*y(k,41) + rxt(k,241) &
                      *y(k,43) + 3.000_r8*rxt(k,246)*y(k,55) + rxt(k,225)*y(k,81)
         mat(k,2146) = rxt(k,217)*y(k,59)
         mat(k,1708) = rxt(k,217)*y(k,19)
         mat(k,2124) = (rxt(k,530)+rxt(k,535))*y(k,91)
         mat(k,779) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,2166) = -(4._r8*rxt(k,214)*y(k,19) + (rxt(k,215) + rxt(k,216) + rxt(k,217) &
                      ) * y(k,59) + rxt(k,218)*y(k,205) + rxt(k,219)*y(k,124) &
                      + rxt(k,220)*y(k,125) + rxt(k,222)*y(k,134) + rxt(k,223) &
                      *y(k,219) + rxt(k,469)*y(k,151))
         mat(k,1728) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,19)
         mat(k,2119) = -rxt(k,218)*y(k,19)
         mat(k,1929) = -rxt(k,219)*y(k,19)
         mat(k,1973) = -rxt(k,220)*y(k,19)
         mat(k,2197) = -rxt(k,222)*y(k,19)
         mat(k,1702) = -rxt(k,223)*y(k,19)
         mat(k,1243) = -rxt(k,469)*y(k,19)
         mat(k,1423) = rxt(k,213)*y(k,135)
         mat(k,564) = rxt(k,221)*y(k,134)
         mat(k,811) = rxt(k,231)*y(k,218)
         mat(k,785) = rxt(k,226)*y(k,134)
         mat(k,2197) = mat(k,2197) + rxt(k,221)*y(k,20) + rxt(k,226)*y(k,91)
         mat(k,2258) = rxt(k,213)*y(k,17)
         mat(k,1538) = rxt(k,231)*y(k,81)
         mat(k,558) = -(rxt(k,221)*y(k,134))
         mat(k,2172) = -rxt(k,221)*y(k,20)
         mat(k,2148) = rxt(k,220)*y(k,125)
         mat(k,1941) = rxt(k,220)*y(k,19)
         mat(k,237) = -(rxt(k,402)*y(k,219))
         mat(k,1574) = -rxt(k,402)*y(k,22)
         mat(k,1854) = rxt(k,405)*y(k,195)
         mat(k,435) = rxt(k,405)*y(k,124)
         mat(k,338) = -(rxt(k,404)*y(k,219))
         mat(k,1589) = -rxt(k,404)*y(k,23)
         mat(k,436) = rxt(k,403)*y(k,205)
         mat(k,2035) = rxt(k,403)*y(k,195)
         mat(k,293) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,219))
         mat(k,1978) = -rxt(k,276)*y(k,24)
         mat(k,1583) = -rxt(k,277)*y(k,24)
         mat(k,550) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,135) + rxt(k,304)*y(k,219))
         mat(k,1983) = -rxt(k,278)*y(k,25)
         mat(k,2210) = -rxt(k,279)*y(k,25)
         mat(k,1620) = -rxt(k,304)*y(k,25)
         mat(k,263) = -(rxt(k,284)*y(k,219))
         mat(k,1580) = -rxt(k,284)*y(k,26)
         mat(k,822) = .800_r8*rxt(k,280)*y(k,196) + .200_r8*rxt(k,281)*y(k,200)
         mat(k,1789) = .200_r8*rxt(k,281)*y(k,196)
         mat(k,343) = -(rxt(k,285)*y(k,219))
         mat(k,1590) = -rxt(k,285)*y(k,27)
         mat(k,823) = rxt(k,282)*y(k,205)
         mat(k,2036) = rxt(k,282)*y(k,196)
         mat(k,306) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,219))
         mat(k,1979) = -rxt(k,286)*y(k,28)
         mat(k,1585) = -rxt(k,287)*y(k,28)
         mat(k,1026) = -(rxt(k,307)*y(k,126) + rxt(k,308)*y(k,135) + rxt(k,325) &
                      *y(k,219))
         mat(k,1748) = -rxt(k,307)*y(k,29)
         mat(k,2224) = -rxt(k,308)*y(k,29)
         mat(k,1664) = -rxt(k,325)*y(k,29)
         mat(k,838) = .130_r8*rxt(k,385)*y(k,135)
         mat(k,2224) = mat(k,2224) + .130_r8*rxt(k,385)*y(k,98)
         mat(k,411) = -(rxt(k,312)*y(k,219))
         mat(k,1600) = -rxt(k,312)*y(k,30)
         mat(k,788) = rxt(k,310)*y(k,205)
         mat(k,2040) = rxt(k,310)*y(k,197)
         mat(k,111) = -(rxt(k,313)*y(k,219))
         mat(k,1558) = -rxt(k,313)*y(k,31)
         mat(k,267) = -(rxt(k,408)*y(k,219))
         mat(k,1581) = -rxt(k,408)*y(k,32)
         mat(k,625) = rxt(k,406)*y(k,205)
         mat(k,2030) = rxt(k,406)*y(k,198)
         mat(k,101) = -(rxt(k,200)*y(k,218))
         mat(k,1501) = -rxt(k,200)*y(k,33)
         mat(k,143) = -(rxt(k,201)*y(k,218))
         mat(k,1506) = -rxt(k,201)*y(k,34)
         mat(k,148) = -(rxt(k,227)*y(k,218))
         mat(k,1507) = -rxt(k,227)*y(k,35)
         mat(k,115) = -(rxt(k,202)*y(k,218))
         mat(k,1503) = -rxt(k,202)*y(k,36)
         mat(k,153) = -(rxt(k,203)*y(k,218))
         mat(k,1508) = -rxt(k,203)*y(k,37)
         mat(k,119) = -(rxt(k,204)*y(k,218))
         mat(k,1504) = -rxt(k,204)*y(k,38)
         mat(k,158) = -(rxt(k,205)*y(k,218))
         mat(k,1509) = -rxt(k,205)*y(k,39)
         mat(k,123) = -(rxt(k,206)*y(k,218))
         mat(k,1505) = -rxt(k,206)*y(k,40)
         mat(k,469) = -(rxt(k,238)*y(k,56) + rxt(k,239)*y(k,219) + rxt(k,250)*y(k,218))
         mat(k,1982) = -rxt(k,238)*y(k,41)
         mat(k,1609) = -rxt(k,239)*y(k,41)
         mat(k,1518) = -rxt(k,250)*y(k,41)
         mat(k,1487) = -(rxt(k,175)*y(k,56) + rxt(k,211)*y(k,17) + rxt(k,255)*y(k,205) &
                      + rxt(k,256)*y(k,126) + rxt(k,257)*y(k,134) + rxt(k,258) &
                      *y(k,219))
         mat(k,2001) = -rxt(k,175)*y(k,42)
         mat(k,1419) = -rxt(k,211)*y(k,42)
         mat(k,2108) = -rxt(k,255)*y(k,42)
         mat(k,1774) = -rxt(k,256)*y(k,42)
         mat(k,2186) = -rxt(k,257)*y(k,42)
         mat(k,1691) = -rxt(k,258)*y(k,42)
         mat(k,640) = .400_r8*rxt(k,356)*y(k,219)
         mat(k,1011) = .340_r8*rxt(k,440)*y(k,135)
         mat(k,355) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,554) = rxt(k,279)*y(k,135)
         mat(k,1033) = .500_r8*rxt(k,308)*y(k,135)
         mat(k,615) = .500_r8*rxt(k,296)*y(k,219)
         mat(k,800) = rxt(k,263)*y(k,219)
         mat(k,458) = .300_r8*rxt(k,264)*y(k,219)
         mat(k,1435) = (rxt(k,272)+rxt(k,273))*y(k,218)
         mat(k,1717) = rxt(k,182)*y(k,200)
         mat(k,1100) = .800_r8*rxt(k,301)*y(k,219)
         mat(k,846) = .910_r8*rxt(k,385)*y(k,135)
         mat(k,588) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,1201) = .800_r8*rxt(k,380)*y(k,200)
         mat(k,1216) = .120_r8*rxt(k,338)*y(k,135)
         mat(k,578) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,961) = .340_r8*rxt(k,443)*y(k,135)
         mat(k,1339) = .600_r8*rxt(k,352)*y(k,135)
         mat(k,1918) = .100_r8*rxt(k,358)*y(k,191) + rxt(k,262)*y(k,200) &
                      + .500_r8*rxt(k,329)*y(k,202) + .500_r8*rxt(k,298)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,207) + .250_r8*rxt(k,336)*y(k,211) &
                      + rxt(k,345)*y(k,213) + rxt(k,319)*y(k,221) + rxt(k,323) &
                      *y(k,222) + .340_r8*rxt(k,452)*y(k,223) + .320_r8*rxt(k,457) &
                      *y(k,224) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1774) = mat(k,1774) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,207) &
                      + .250_r8*rxt(k,335)*y(k,211) + rxt(k,346)*y(k,213)
         mat(k,2247) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25) &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,98) &
                      + .120_r8*rxt(k,338)*y(k,105) + .340_r8*rxt(k,443)*y(k,110) &
                      + .600_r8*rxt(k,352)*y(k,111)
         mat(k,529) = rxt(k,303)*y(k,219)
         mat(k,1065) = .680_r8*rxt(k,461)*y(k,219)
         mat(k,893) = .100_r8*rxt(k,358)*y(k,124)
         mat(k,827) = .700_r8*rxt(k,281)*y(k,200)
         mat(k,792) = rxt(k,309)*y(k,200)
         mat(k,1391) = rxt(k,292)*y(k,200) + rxt(k,365)*y(k,207) + .250_r8*rxt(k,332) &
                      *y(k,211) + rxt(k,341)*y(k,213) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1826) = rxt(k,182)*y(k,59) + .800_r8*rxt(k,380)*y(k,101) + rxt(k,262) &
                      *y(k,124) + .700_r8*rxt(k,281)*y(k,196) + rxt(k,309)*y(k,197) &
                      + rxt(k,292)*y(k,199) + (4.000_r8*rxt(k,259)+2.000_r8*rxt(k,260)) &
                      *y(k,200) + 1.500_r8*rxt(k,366)*y(k,207) + .750_r8*rxt(k,371) &
                      *y(k,208) + .880_r8*rxt(k,333)*y(k,211) + 2.000_r8*rxt(k,342) &
                      *y(k,213) + .750_r8*rxt(k,445)*y(k,217) + .800_r8*rxt(k,321) &
                      *y(k,222) + .930_r8*rxt(k,450)*y(k,223) + .950_r8*rxt(k,455) &
                      *y(k,224) + .800_r8*rxt(k,391)*y(k,227)
         mat(k,570) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,710) = .500_r8*rxt(k,298)*y(k,124)
         mat(k,2108) = mat(k,2108) + .450_r8*rxt(k,343)*y(k,213) + .150_r8*rxt(k,322) &
                      *y(k,222)
         mat(k,1264) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) + rxt(k,365) &
                      *y(k,199) + 1.500_r8*rxt(k,366)*y(k,200)
         mat(k,1296) = .750_r8*rxt(k,371)*y(k,200)
         mat(k,1317) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,199) + .880_r8*rxt(k,333)*y(k,200)
         mat(k,1359) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,341)*y(k,199) &
                      + 2.000_r8*rxt(k,342)*y(k,200) + .450_r8*rxt(k,343)*y(k,205) &
                      + 4.000_r8*rxt(k,344)*y(k,213)
         mat(k,1052) = .750_r8*rxt(k,445)*y(k,200)
         mat(k,1527) = (rxt(k,272)+rxt(k,273))*y(k,54)
         mat(k,1691) = mat(k,1691) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,263)*y(k,52) + .300_r8*rxt(k,264)*y(k,53) &
                      + .800_r8*rxt(k,301)*y(k,74) + .300_r8*rxt(k,376)*y(k,99) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,140) &
                      + .680_r8*rxt(k,461)*y(k,180)
         mat(k,773) = rxt(k,319)*y(k,124)
         mat(k,1135) = rxt(k,323)*y(k,124) + .800_r8*rxt(k,321)*y(k,200) &
                      + .150_r8*rxt(k,322)*y(k,205)
         mat(k,1116) = .340_r8*rxt(k,452)*y(k,124) + .930_r8*rxt(k,450)*y(k,200)
         mat(k,918) = .320_r8*rxt(k,457)*y(k,124) + .950_r8*rxt(k,455)*y(k,200)
         mat(k,1178) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,390)*y(k,199) &
                      + .800_r8*rxt(k,391)*y(k,200)
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
         mat(k,596) = -(rxt(k,240)*y(k,56) + rxt(k,241)*y(k,219) + rxt(k,251)*y(k,218))
         mat(k,1984) = -rxt(k,240)*y(k,43)
         mat(k,1625) = -rxt(k,241)*y(k,43)
         mat(k,1519) = -rxt(k,251)*y(k,43)
         mat(k,127) = -(rxt(k,242)*y(k,219))
         mat(k,1559) = -rxt(k,242)*y(k,44)
         mat(k,1087) = -(rxt(k,288)*y(k,126) + rxt(k,289)*y(k,219))
         mat(k,1752) = -rxt(k,288)*y(k,45)
         mat(k,1668) = -rxt(k,289)*y(k,45)
         mat(k,638) = .800_r8*rxt(k,356)*y(k,219)
         mat(k,354) = rxt(k,327)*y(k,126)
         mat(k,264) = rxt(k,284)*y(k,219)
         mat(k,345) = .500_r8*rxt(k,285)*y(k,219)
         mat(k,1027) = .500_r8*rxt(k,308)*y(k,135)
         mat(k,1329) = .100_r8*rxt(k,352)*y(k,135)
         mat(k,1897) = .400_r8*rxt(k,358)*y(k,191) + rxt(k,283)*y(k,196) &
                      + .270_r8*rxt(k,311)*y(k,197) + rxt(k,329)*y(k,202) + rxt(k,348) &
                      *y(k,215) + rxt(k,319)*y(k,221)
         mat(k,1752) = mat(k,1752) + rxt(k,327)*y(k,16)
         mat(k,2227) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,891) = .400_r8*rxt(k,358)*y(k,124)
         mat(k,826) = rxt(k,283)*y(k,124) + 3.200_r8*rxt(k,280)*y(k,196) &
                      + .800_r8*rxt(k,281)*y(k,200)
         mat(k,791) = .270_r8*rxt(k,311)*y(k,124)
         mat(k,1807) = .800_r8*rxt(k,281)*y(k,196)
         mat(k,568) = rxt(k,329)*y(k,124)
         mat(k,2087) = .200_r8*rxt(k,347)*y(k,215)
         mat(k,680) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,205)
         mat(k,1668) = mat(k,1668) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26) &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,771) = rxt(k,319)*y(k,124)
         mat(k,367) = -(rxt(k,243)*y(k,56) + rxt(k,244)*y(k,219))
         mat(k,1980) = -rxt(k,243)*y(k,46)
         mat(k,1594) = -rxt(k,244)*y(k,46)
         mat(k,104) = -(rxt(k,290)*y(k,219))
         mat(k,1557) = -rxt(k,290)*y(k,47)
         mat(k,924) = -(rxt(k,326)*y(k,219))
         mat(k,1658) = -rxt(k,326)*y(k,48)
         mat(k,637) = .800_r8*rxt(k,356)*y(k,219)
         mat(k,997) = .520_r8*rxt(k,440)*y(k,135)
         mat(k,353) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,948) = .520_r8*rxt(k,443)*y(k,135)
         mat(k,1890) = .250_r8*rxt(k,358)*y(k,191) + .820_r8*rxt(k,311)*y(k,197) &
                      + .500_r8*rxt(k,329)*y(k,202) + .270_r8*rxt(k,452)*y(k,223) &
                      + .040_r8*rxt(k,457)*y(k,224)
         mat(k,1742) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,2218) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,110)
         mat(k,1060) = .500_r8*rxt(k,461)*y(k,219)
         mat(k,890) = .250_r8*rxt(k,358)*y(k,124)
         mat(k,790) = .820_r8*rxt(k,311)*y(k,124) + .820_r8*rxt(k,309)*y(k,200)
         mat(k,1801) = .820_r8*rxt(k,309)*y(k,197) + .150_r8*rxt(k,450)*y(k,223) &
                      + .025_r8*rxt(k,455)*y(k,224)
         mat(k,567) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,1658) = mat(k,1658) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,180)
         mat(k,1108) = .270_r8*rxt(k,452)*y(k,124) + .150_r8*rxt(k,450)*y(k,200)
         mat(k,915) = .040_r8*rxt(k,457)*y(k,124) + .025_r8*rxt(k,455)*y(k,200)
         mat(k,1223) = -(rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219))
         mat(k,1762) = -rxt(k,314)*y(k,49)
         mat(k,1678) = -rxt(k,315)*y(k,49)
         mat(k,1143) = rxt(k,316)*y(k,219)
         mat(k,1212) = .880_r8*rxt(k,338)*y(k,135)
         mat(k,1332) = .500_r8*rxt(k,352)*y(k,135)
         mat(k,1907) = .170_r8*rxt(k,411)*y(k,201) + .050_r8*rxt(k,374)*y(k,208) &
                      + .250_r8*rxt(k,336)*y(k,211) + .170_r8*rxt(k,417)*y(k,214) &
                      + .400_r8*rxt(k,427)*y(k,225) + .250_r8*rxt(k,393)*y(k,227) &
                      + .540_r8*rxt(k,433)*y(k,228) + .510_r8*rxt(k,436)*y(k,230)
         mat(k,1762) = mat(k,1762) + .050_r8*rxt(k,375)*y(k,208) + .250_r8*rxt(k,335) &
                      *y(k,211) + .250_r8*rxt(k,394)*y(k,227)
         mat(k,860) = rxt(k,317)*y(k,219)
         mat(k,2235) = .880_r8*rxt(k,338)*y(k,105) + .500_r8*rxt(k,352)*y(k,111)
         mat(k,1382) = .250_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1816) = .240_r8*rxt(k,333)*y(k,211) + .500_r8*rxt(k,321)*y(k,222) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,763) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,205)
         mat(k,2096) = .070_r8*rxt(k,410)*y(k,201) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1289) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1312) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,199) + .240_r8*rxt(k,333)*y(k,200)
         mat(k,878) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1678) = mat(k,1678) + rxt(k,316)*y(k,95) + rxt(k,317)*y(k,127)
         mat(k,1133) = .500_r8*rxt(k,321)*y(k,200)
         mat(k,731) = .400_r8*rxt(k,427)*y(k,124)
         mat(k,1176) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,755) = .540_r8*rxt(k,433)*y(k,124)
         mat(k,503) = .510_r8*rxt(k,436)*y(k,124)
         mat(k,686) = -(rxt(k,295)*y(k,219))
         mat(k,1635) = -rxt(k,295)*y(k,50)
         mat(k,1021) = .120_r8*rxt(k,308)*y(k,135)
         mat(k,2212) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1372) = .100_r8*rxt(k,292)*y(k,200) + .150_r8*rxt(k,293)*y(k,205)
         mat(k,1794) = .100_r8*rxt(k,292)*y(k,199)
         mat(k,2062) = .150_r8*rxt(k,293)*y(k,199) + .150_r8*rxt(k,343)*y(k,213)
         mat(k,1351) = .150_r8*rxt(k,343)*y(k,205)
         mat(k,612) = -(rxt(k,296)*y(k,219))
         mat(k,1627) = -rxt(k,296)*y(k,51)
         mat(k,1371) = .400_r8*rxt(k,293)*y(k,205)
         mat(k,2056) = .400_r8*rxt(k,293)*y(k,199) + .400_r8*rxt(k,343)*y(k,213)
         mat(k,1350) = .400_r8*rxt(k,343)*y(k,205)
         mat(k,799) = -(rxt(k,263)*y(k,219))
         mat(k,1645) = -rxt(k,263)*y(k,52)
         mat(k,1189) = .200_r8*rxt(k,380)*y(k,200)
         mat(k,824) = .300_r8*rxt(k,281)*y(k,200)
         mat(k,1796) = .200_r8*rxt(k,380)*y(k,101) + .300_r8*rxt(k,281)*y(k,196) &
                      + 2.000_r8*rxt(k,260)*y(k,200) + .250_r8*rxt(k,366)*y(k,207) &
                      + .250_r8*rxt(k,371)*y(k,208) + .250_r8*rxt(k,333)*y(k,211) &
                      + .250_r8*rxt(k,445)*y(k,217) + .500_r8*rxt(k,321)*y(k,222) &
                      + .250_r8*rxt(k,450)*y(k,223) + .250_r8*rxt(k,455)*y(k,224) &
                      + .300_r8*rxt(k,391)*y(k,227)
         mat(k,1249) = .250_r8*rxt(k,366)*y(k,200)
         mat(k,1279) = .250_r8*rxt(k,371)*y(k,200)
         mat(k,1307) = .250_r8*rxt(k,333)*y(k,200)
         mat(k,1045) = .250_r8*rxt(k,445)*y(k,200)
         mat(k,1130) = .500_r8*rxt(k,321)*y(k,200)
         mat(k,1106) = .250_r8*rxt(k,450)*y(k,200)
         mat(k,913) = .250_r8*rxt(k,455)*y(k,200)
         mat(k,1169) = .300_r8*rxt(k,391)*y(k,200)
         mat(k,456) = -(rxt(k,264)*y(k,219))
         mat(k,1606) = -rxt(k,264)*y(k,53)
         mat(k,1792) = rxt(k,261)*y(k,205)
         mat(k,2047) = rxt(k,261)*y(k,200)
         mat(k,1432) = -(rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,265)*y(k,219) &
                      + (rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,218))
         mat(k,1998) = -rxt(k,176)*y(k,54)
         mat(k,868) = -rxt(k,232)*y(k,54)
         mat(k,1688) = -rxt(k,265)*y(k,54)
         mat(k,1524) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,54)
         mat(k,1032) = .100_r8*rxt(k,308)*y(k,135)
         mat(k,2244) = .100_r8*rxt(k,308)*y(k,29)
         mat(k,387) = -(rxt(k,228)*y(k,218) + rxt(k,245)*y(k,56) + rxt(k,246)*y(k,219))
         mat(k,1517) = -rxt(k,228)*y(k,55)
         mat(k,1981) = -rxt(k,245)*y(k,55)
         mat(k,1596) = -rxt(k,246)*y(k,55)
         mat(k,2009) = -(rxt(k,175)*y(k,42) + rxt(k,176)*y(k,54) + rxt(k,177)*y(k,77) &
                      + rxt(k,178)*y(k,79) + (rxt(k,179) + rxt(k,180)) * y(k,205) &
                      + rxt(k,181)*y(k,135) + rxt(k,188)*y(k,60) + rxt(k,197)*y(k,92) &
                      + rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43) + rxt(k,243)*y(k,46) &
                      + rxt(k,245)*y(k,55) + rxt(k,286)*y(k,28))
         mat(k,1494) = -rxt(k,175)*y(k,56)
         mat(k,1440) = -rxt(k,176)*y(k,56)
         mat(k,1411) = -rxt(k,177)*y(k,56)
         mat(k,607) = -rxt(k,178)*y(k,56)
         mat(k,2116) = -(rxt(k,179) + rxt(k,180)) * y(k,56)
         mat(k,2255) = -rxt(k,181)*y(k,56)
         mat(k,907) = -rxt(k,188)*y(k,56)
         mat(k,818) = -rxt(k,197)*y(k,56)
         mat(k,473) = -rxt(k,238)*y(k,56)
         mat(k,601) = -rxt(k,240)*y(k,56)
         mat(k,371) = -rxt(k,243)*y(k,56)
         mat(k,391) = -rxt(k,245)*y(k,56)
         mat(k,309) = -rxt(k,286)*y(k,56)
         mat(k,2163) = rxt(k,216)*y(k,59)
         mat(k,103) = 4.000_r8*rxt(k,200)*y(k,218)
         mat(k,147) = rxt(k,201)*y(k,218)
         mat(k,118) = 2.000_r8*rxt(k,202)*y(k,218)
         mat(k,157) = 2.000_r8*rxt(k,203)*y(k,218)
         mat(k,122) = 2.000_r8*rxt(k,204)*y(k,218)
         mat(k,162) = rxt(k,205)*y(k,218)
         mat(k,126) = 2.000_r8*rxt(k,206)*y(k,218)
         mat(k,129) = 3.000_r8*rxt(k,242)*y(k,219)
         mat(k,371) = mat(k,371) + rxt(k,244)*y(k,219)
         mat(k,1725) = rxt(k,216)*y(k,19) + (4.000_r8*rxt(k,183)+2.000_r8*rxt(k,185)) &
                      *y(k,59) + rxt(k,187)*y(k,124) + rxt(k,192)*y(k,134) &
                      + rxt(k,470)*y(k,151) + rxt(k,182)*y(k,200) + rxt(k,193) &
                      *y(k,219)
         mat(k,251) = rxt(k,237)*y(k,218)
         mat(k,247) = rxt(k,252)*y(k,218) + rxt(k,247)*y(k,219)
         mat(k,257) = rxt(k,253)*y(k,218) + rxt(k,248)*y(k,219)
         mat(k,304) = rxt(k,254)*y(k,218) + rxt(k,249)*y(k,219)
         mat(k,2139) = rxt(k,195)*y(k,134) + rxt(k,207)*y(k,218) + rxt(k,196)*y(k,219)
         mat(k,1926) = rxt(k,187)*y(k,59)
         mat(k,2194) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,85)
         mat(k,1241) = rxt(k,470)*y(k,59)
         mat(k,1834) = rxt(k,182)*y(k,59)
         mat(k,1535) = 4.000_r8*rxt(k,200)*y(k,33) + rxt(k,201)*y(k,34) &
                      + 2.000_r8*rxt(k,202)*y(k,36) + 2.000_r8*rxt(k,203)*y(k,37) &
                      + 2.000_r8*rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39) &
                      + 2.000_r8*rxt(k,206)*y(k,40) + rxt(k,237)*y(k,65) + rxt(k,252) &
                      *y(k,82) + rxt(k,253)*y(k,83) + rxt(k,254)*y(k,84) + rxt(k,207) &
                      *y(k,85)
         mat(k,1699) = 3.000_r8*rxt(k,242)*y(k,44) + rxt(k,244)*y(k,46) + rxt(k,193) &
                      *y(k,59) + rxt(k,247)*y(k,82) + rxt(k,248)*y(k,83) + rxt(k,249) &
                      *y(k,84) + rxt(k,196)*y(k,85)
         mat(k,1977) = rxt(k,188)*y(k,60)
         mat(k,1707) = 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,899) = rxt(k,188)*y(k,56) + (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,2123) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60) + (rxt(k,523) &
                       +rxt(k,529)+rxt(k,534))*y(k,92)
         mat(k,814) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85)
         mat(k,1706) = 2.000_r8*rxt(k,209)*y(k,59)
         mat(k,1720) = -(rxt(k,182)*y(k,200) + (4._r8*rxt(k,183) + 4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185) + 4._r8*rxt(k,209)) * y(k,59) + rxt(k,186) &
                      *y(k,205) + rxt(k,187)*y(k,124) + rxt(k,189)*y(k,125) + rxt(k,192) &
                      *y(k,134) + (rxt(k,193) + rxt(k,194)) * y(k,219) + (rxt(k,215) &
                      + rxt(k,216) + rxt(k,217)) * y(k,19) + rxt(k,470)*y(k,151))
         mat(k,1829) = -rxt(k,182)*y(k,59)
         mat(k,2111) = -rxt(k,186)*y(k,59)
         mat(k,1921) = -rxt(k,187)*y(k,59)
         mat(k,1965) = -rxt(k,189)*y(k,59)
         mat(k,2189) = -rxt(k,192)*y(k,59)
         mat(k,1694) = -(rxt(k,193) + rxt(k,194)) * y(k,59)
         mat(k,2158) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,59)
         mat(k,1238) = -rxt(k,470)*y(k,59)
         mat(k,2004) = rxt(k,197)*y(k,92) + rxt(k,181)*y(k,135) + rxt(k,180)*y(k,205)
         mat(k,904) = rxt(k,190)*y(k,134)
         mat(k,2134) = rxt(k,208)*y(k,218)
         mat(k,817) = rxt(k,197)*y(k,56) + rxt(k,198)*y(k,134) + rxt(k,199)*y(k,219)
         mat(k,2189) = mat(k,2189) + rxt(k,190)*y(k,60) + rxt(k,198)*y(k,92)
         mat(k,2250) = rxt(k,181)*y(k,56)
         mat(k,330) = rxt(k,475)*y(k,151)
         mat(k,1238) = mat(k,1238) + rxt(k,475)*y(k,137)
         mat(k,2111) = mat(k,2111) + rxt(k,180)*y(k,56)
         mat(k,1530) = rxt(k,208)*y(k,85)
         mat(k,1694) = mat(k,1694) + rxt(k,199)*y(k,92)
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
         mat(k,901) = -(rxt(k,188)*y(k,56) + rxt(k,190)*y(k,134) + rxt(k,191)*y(k,219) &
                      + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85))
         mat(k,1989) = -rxt(k,188)*y(k,60)
         mat(k,2178) = -rxt(k,190)*y(k,60)
         mat(k,1656) = -rxt(k,191)*y(k,60)
         mat(k,2127) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60)
         mat(k,1712) = rxt(k,189)*y(k,125)
         mat(k,1950) = rxt(k,189)*y(k,59)
         mat(k,1125) = -(rxt(k,275)*y(k,219))
         mat(k,1671) = -rxt(k,275)*y(k,62)
         mat(k,1006) = .230_r8*rxt(k,440)*y(k,135)
         mat(k,1416) = rxt(k,211)*y(k,42)
         mat(k,296) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,553) = .630_r8*rxt(k,279)*y(k,135)
         mat(k,1028) = .560_r8*rxt(k,308)*y(k,135)
         mat(k,1481) = rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56) + rxt(k,256)*y(k,126) &
                      + rxt(k,257)*y(k,134) + rxt(k,258)*y(k,219)
         mat(k,368) = rxt(k,243)*y(k,56)
         mat(k,1222) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219)
         mat(k,1994) = rxt(k,175)*y(k,42) + rxt(k,243)*y(k,46)
         mat(k,982) = rxt(k,302)*y(k,219)
         mat(k,839) = .620_r8*rxt(k,385)*y(k,135)
         mat(k,1210) = .650_r8*rxt(k,338)*y(k,135)
         mat(k,956) = .230_r8*rxt(k,443)*y(k,135)
         mat(k,1330) = .560_r8*rxt(k,352)*y(k,135)
         mat(k,1900) = .170_r8*rxt(k,411)*y(k,201) + .220_r8*rxt(k,336)*y(k,211) &
                      + .400_r8*rxt(k,414)*y(k,212) + .350_r8*rxt(k,417)*y(k,214) &
                      + .225_r8*rxt(k,452)*y(k,223) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1755) = rxt(k,256)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,211) + .500_r8*rxt(k,394)*y(k,227)
         mat(k,2179) = rxt(k,257)*y(k,42) + rxt(k,464)*y(k,138)
         mat(k,2230) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25) &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,98) &
                      + .650_r8*rxt(k,338)*y(k,105) + .230_r8*rxt(k,443)*y(k,110) &
                      + .560_r8*rxt(k,352)*y(k,111)
         mat(k,362) = rxt(k,464)*y(k,134) + rxt(k,465)*y(k,219)
         mat(k,1062) = .700_r8*rxt(k,461)*y(k,219)
         mat(k,1377) = .220_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1810) = .110_r8*rxt(k,333)*y(k,211) + .125_r8*rxt(k,450)*y(k,223) &
                      + .200_r8*rxt(k,391)*y(k,227)
         mat(k,762) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,205)
         mat(k,2090) = .070_r8*rxt(k,410)*y(k,201) + .160_r8*rxt(k,413)*y(k,212) &
                      + .140_r8*rxt(k,416)*y(k,214)
         mat(k,1309) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,199) + .110_r8*rxt(k,333)*y(k,200)
         mat(k,717) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,205)
         mat(k,877) = .350_r8*rxt(k,417)*y(k,124) + .140_r8*rxt(k,416)*y(k,205)
         mat(k,1671) = mat(k,1671) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,258)*y(k,42) &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,75) + rxt(k,465)*y(k,138) &
                      + .700_r8*rxt(k,461)*y(k,180)
         mat(k,1112) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,200)
         mat(k,1173) = .250_r8*rxt(k,393)*y(k,124) + .500_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .200_r8*rxt(k,391)*y(k,200)
         mat(k,994) = .270_r8*rxt(k,440)*y(k,135)
         mat(k,1023) = .200_r8*rxt(k,308)*y(k,135)
         mat(k,687) = rxt(k,295)*y(k,219)
         mat(k,613) = .500_r8*rxt(k,296)*y(k,219)
         mat(k,1124) = rxt(k,275)*y(k,219)
         mat(k,1096) = .800_r8*rxt(k,301)*y(k,219)
         mat(k,980) = rxt(k,302)*y(k,219)
         mat(k,930) = rxt(k,267)*y(k,219)
         mat(k,575) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,945) = .270_r8*rxt(k,443)*y(k,135)
         mat(k,1326) = .100_r8*rxt(k,352)*y(k,135)
         mat(k,1884) = rxt(k,294)*y(k,199) + .900_r8*rxt(k,452)*y(k,223)
         mat(k,2214) = .270_r8*rxt(k,440)*y(k,6) + .200_r8*rxt(k,308)*y(k,29) &
                      + .270_r8*rxt(k,443)*y(k,110) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,1059) = 1.800_r8*rxt(k,461)*y(k,219)
         mat(k,1373) = rxt(k,294)*y(k,124) + 4.000_r8*rxt(k,291)*y(k,199) &
                      + .900_r8*rxt(k,292)*y(k,200) + rxt(k,365)*y(k,207) &
                      + 2.000_r8*rxt(k,341)*y(k,213) + rxt(k,390)*y(k,227)
         mat(k,1797) = .900_r8*rxt(k,292)*y(k,199) + rxt(k,342)*y(k,213) &
                      + .500_r8*rxt(k,450)*y(k,223)
         mat(k,2073) = .450_r8*rxt(k,343)*y(k,213)
         mat(k,1250) = rxt(k,365)*y(k,199)
         mat(k,1352) = 2.000_r8*rxt(k,341)*y(k,199) + rxt(k,342)*y(k,200) &
                      + .450_r8*rxt(k,343)*y(k,205) + 4.000_r8*rxt(k,344)*y(k,213)
         mat(k,1646) = rxt(k,295)*y(k,50) + .500_r8*rxt(k,296)*y(k,51) + rxt(k,275) &
                      *y(k,62) + .800_r8*rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) &
                      + rxt(k,267)*y(k,87) + .500_r8*rxt(k,351)*y(k,109) &
                      + 1.800_r8*rxt(k,461)*y(k,180)
         mat(k,1107) = .900_r8*rxt(k,452)*y(k,124) + .500_r8*rxt(k,450)*y(k,200)
         mat(k,1170) = rxt(k,390)*y(k,199)
         mat(k,240) = -(rxt(k,236)*y(k,218))
         mat(k,1512) = -rxt(k,236)*y(k,64)
         mat(k,144) = rxt(k,201)*y(k,218)
         mat(k,149) = rxt(k,227)*y(k,218)
         mat(k,154) = rxt(k,203)*y(k,218)
         mat(k,120) = 2.000_r8*rxt(k,204)*y(k,218)
         mat(k,159) = 2.000_r8*rxt(k,205)*y(k,218)
         mat(k,124) = rxt(k,206)*y(k,218)
         mat(k,108) = 2.000_r8*rxt(k,229)*y(k,218)
         mat(k,252) = rxt(k,253)*y(k,218) + rxt(k,248)*y(k,219)
         mat(k,299) = rxt(k,254)*y(k,218) + rxt(k,249)*y(k,219)
         mat(k,1512) = mat(k,1512) + rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35) &
                      + rxt(k,203)*y(k,37) + 2.000_r8*rxt(k,204)*y(k,38) &
                      + 2.000_r8*rxt(k,205)*y(k,39) + rxt(k,206)*y(k,40) &
                      + 2.000_r8*rxt(k,229)*y(k,78) + rxt(k,253)*y(k,83) + rxt(k,254) &
                      *y(k,84)
         mat(k,1575) = rxt(k,248)*y(k,83) + rxt(k,249)*y(k,84)
         mat(k,248) = -(rxt(k,237)*y(k,218))
         mat(k,1514) = -rxt(k,237)*y(k,65)
         mat(k,116) = rxt(k,202)*y(k,218)
         mat(k,155) = rxt(k,203)*y(k,218)
         mat(k,244) = rxt(k,252)*y(k,218) + rxt(k,247)*y(k,219)
         mat(k,1514) = mat(k,1514) + rxt(k,202)*y(k,36) + rxt(k,203)*y(k,37) &
                      + rxt(k,252)*y(k,82)
         mat(k,1577) = rxt(k,247)*y(k,82)
         mat(k,196) = -(rxt(k,409)*y(k,219))
         mat(k,1566) = -rxt(k,409)*y(k,66)
         mat(k,190) = .180_r8*rxt(k,429)*y(k,219)
         mat(k,1566) = mat(k,1566) + .180_r8*rxt(k,429)*y(k,182)
         mat(k,284) = -(rxt(k,462)*y(k,126) + (rxt(k,463) + rxt(k,477)) * y(k,219))
         mat(k,1735) = -rxt(k,462)*y(k,67)
         mat(k,1582) = -(rxt(k,463) + rxt(k,477)) * y(k,67)
         mat(k,706) = rxt(k,297)*y(k,205)
         mat(k,2028) = rxt(k,297)*y(k,204)
         mat(k,866) = -(rxt(k,232)*y(k,54) + rxt(k,233)*y(k,77) + rxt(k,234)*y(k,231) &
                      + rxt(k,235)*y(k,89))
         mat(k,1429) = -rxt(k,232)*y(k,73)
         mat(k,1402) = -rxt(k,233)*y(k,73)
         mat(k,2266) = -rxt(k,234)*y(k,73)
         mat(k,1461) = -rxt(k,235)*y(k,73)
         mat(k,150) = rxt(k,227)*y(k,218)
         mat(k,160) = rxt(k,205)*y(k,218)
         mat(k,241) = 2.000_r8*rxt(k,236)*y(k,218)
         mat(k,249) = rxt(k,237)*y(k,218)
         mat(k,1521) = rxt(k,227)*y(k,35) + rxt(k,205)*y(k,39) + 2.000_r8*rxt(k,236) &
                      *y(k,64) + rxt(k,237)*y(k,65)
         mat(k,1098) = -(rxt(k,301)*y(k,219))
         mat(k,1669) = -rxt(k,301)*y(k,74)
         mat(k,584) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,536) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,377) = rxt(k,388)*y(k,219)
         mat(k,1898) = .050_r8*rxt(k,374)*y(k,208) + .530_r8*rxt(k,336)*y(k,211) &
                      + .225_r8*rxt(k,452)*y(k,223) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1753) = .050_r8*rxt(k,375)*y(k,208) + .530_r8*rxt(k,335)*y(k,211) &
                      + .250_r8*rxt(k,394)*y(k,227)
         mat(k,1376) = .530_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1808) = .260_r8*rxt(k,333)*y(k,211) + .125_r8*rxt(k,450)*y(k,223) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,1283) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1308) = .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335)*y(k,126) &
                      + .530_r8*rxt(k,332)*y(k,199) + .260_r8*rxt(k,333)*y(k,200)
         mat(k,1669) = mat(k,1669) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + rxt(k,388)*y(k,115)
         mat(k,1110) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,200)
         mat(k,1172) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,981) = -(rxt(k,302)*y(k,219))
         mat(k,1662) = -rxt(k,302)*y(k,75)
         mat(k,295) = .650_r8*rxt(k,277)*y(k,219)
         mat(k,1097) = .200_r8*rxt(k,301)*y(k,219)
         mat(k,1074) = rxt(k,389)*y(k,219)
         mat(k,1893) = rxt(k,400)*y(k,193) + .050_r8*rxt(k,374)*y(k,208) &
                      + .400_r8*rxt(k,414)*y(k,212) + .170_r8*rxt(k,417)*y(k,214) &
                      + .700_r8*rxt(k,420)*y(k,220) + .600_r8*rxt(k,427)*y(k,225) &
                      + .250_r8*rxt(k,393)*y(k,227) + .340_r8*rxt(k,433)*y(k,228) &
                      + .170_r8*rxt(k,436)*y(k,230)
         mat(k,1746) = .050_r8*rxt(k,375)*y(k,208) + .250_r8*rxt(k,394)*y(k,227)
         mat(k,487) = rxt(k,400)*y(k,124)
         mat(k,1374) = .250_r8*rxt(k,390)*y(k,227)
         mat(k,1803) = .100_r8*rxt(k,391)*y(k,227)
         mat(k,2084) = .160_r8*rxt(k,413)*y(k,212) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1282) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,716) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,205)
         mat(k,876) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1662) = mat(k,1662) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,74) + rxt(k,389)*y(k,116)
         mat(k,451) = .700_r8*rxt(k,420)*y(k,124)
         mat(k,729) = .600_r8*rxt(k,427)*y(k,124)
         mat(k,1171) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,753) = .340_r8*rxt(k,433)*y(k,124)
         mat(k,502) = .170_r8*rxt(k,436)*y(k,124)
         mat(k,1447) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,205) + rxt(k,141) &
                      *y(k,135))
         mat(k,2106) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,76)
         mat(k,2245) = -rxt(k,141)*y(k,76)
         mat(k,1485) = rxt(k,258)*y(k,219)
         mat(k,1433) = rxt(k,272)*y(k,218)
         mat(k,1999) = rxt(k,177)*y(k,77)
         mat(k,869) = rxt(k,233)*y(k,77)
         mat(k,1405) = rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73) + rxt(k,133)*y(k,134) &
                      + rxt(k,125)*y(k,218) + rxt(k,142)*y(k,219)
         mat(k,808) = rxt(k,231)*y(k,218)
         mat(k,2129) = rxt(k,208)*y(k,218)
         mat(k,494) = rxt(k,163)*y(k,219)
         mat(k,2184) = rxt(k,133)*y(k,77) + rxt(k,145)*y(k,219)
         mat(k,364) = rxt(k,465)*y(k,219)
         mat(k,515) = rxt(k,471)*y(k,219)
         mat(k,1236) = rxt(k,476)*y(k,219)
         mat(k,1525) = rxt(k,272)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,231)*y(k,81) &
                      + rxt(k,208)*y(k,85)
         mat(k,1689) = rxt(k,258)*y(k,42) + rxt(k,142)*y(k,77) + rxt(k,163)*y(k,112) &
                      + rxt(k,145)*y(k,134) + rxt(k,465)*y(k,138) + rxt(k,471) &
                      *y(k,149) + rxt(k,476)*y(k,151)
         mat(k,1403) = -(rxt(k,125)*y(k,218) + rxt(k,133)*y(k,134) + rxt(k,142) &
                      *y(k,219) + rxt(k,177)*y(k,56) + rxt(k,233)*y(k,73))
         mat(k,1522) = -rxt(k,125)*y(k,77)
         mat(k,2181) = -rxt(k,133)*y(k,77)
         mat(k,1686) = -rxt(k,142)*y(k,77)
         mat(k,1996) = -rxt(k,177)*y(k,77)
         mat(k,867) = -rxt(k,233)*y(k,77)
         mat(k,1431) = rxt(k,273)*y(k,218)
         mat(k,1445) = rxt(k,135)*y(k,205)
         mat(k,2103) = rxt(k,135)*y(k,76)
         mat(k,1522) = mat(k,1522) + rxt(k,273)*y(k,54)
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
         mat(k,107) = -(rxt(k,229)*y(k,218))
         mat(k,1502) = -rxt(k,229)*y(k,78)
         mat(k,605) = -(rxt(k,134)*y(k,134) + rxt(k,143)*y(k,219) + rxt(k,178)*y(k,56))
         mat(k,2173) = -rxt(k,134)*y(k,79)
         mat(k,1626) = -rxt(k,143)*y(k,79)
         mat(k,1985) = -rxt(k,178)*y(k,79)
         mat(k,2055) = 2.000_r8*rxt(k,149)*y(k,205)
         mat(k,1626) = mat(k,1626) + 2.000_r8*rxt(k,148)*y(k,219)
         mat(k,258) = rxt(k,478)*y(k,231)
         mat(k,2262) = rxt(k,478)*y(k,153)
         mat(k,806) = -(rxt(k,224)*y(k,134) + rxt(k,225)*y(k,219) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,218))
         mat(k,2175) = -rxt(k,224)*y(k,81)
         mat(k,1647) = -rxt(k,225)*y(k,81)
         mat(k,1520) = -(rxt(k,230) + rxt(k,231)) * y(k,81)
         mat(k,1415) = rxt(k,211)*y(k,42) + rxt(k,212)*y(k,205)
         mat(k,1479) = rxt(k,211)*y(k,17)
         mat(k,2074) = rxt(k,212)*y(k,17)
         mat(k,243) = -(rxt(k,247)*y(k,219) + rxt(k,252)*y(k,218))
         mat(k,1576) = -rxt(k,247)*y(k,82)
         mat(k,1513) = -rxt(k,252)*y(k,82)
         mat(k,253) = -(rxt(k,248)*y(k,219) + rxt(k,253)*y(k,218))
         mat(k,1578) = -rxt(k,248)*y(k,83)
         mat(k,1515) = -rxt(k,253)*y(k,83)
         mat(k,300) = -(rxt(k,249)*y(k,219) + rxt(k,254)*y(k,218))
         mat(k,1584) = -rxt(k,249)*y(k,84)
         mat(k,1516) = -rxt(k,254)*y(k,84)
         mat(k,2141) = -(rxt(k,195)*y(k,134) + rxt(k,196)*y(k,219) + (rxt(k,207) &
                      + rxt(k,208)) * y(k,218) + (rxt(k,523) + rxt(k,529) + rxt(k,534) &
                      ) * y(k,92) + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60) &
                      + (rxt(k,530) + rxt(k,535)) * y(k,91))
         mat(k,2196) = -rxt(k,195)*y(k,85)
         mat(k,1701) = -rxt(k,196)*y(k,85)
         mat(k,1537) = -(rxt(k,207) + rxt(k,208)) * y(k,85)
         mat(k,819) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85)
         mat(k,908) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85)
         mat(k,784) = -(rxt(k,530) + rxt(k,535)) * y(k,85)
         mat(k,310) = rxt(k,286)*y(k,56)
         mat(k,474) = rxt(k,238)*y(k,56)
         mat(k,1496) = rxt(k,175)*y(k,56)
         mat(k,603) = rxt(k,240)*y(k,56)
         mat(k,373) = 2.000_r8*rxt(k,243)*y(k,56)
         mat(k,1442) = rxt(k,176)*y(k,56)
         mat(k,392) = rxt(k,245)*y(k,56)
         mat(k,2011) = rxt(k,286)*y(k,28) + rxt(k,238)*y(k,41) + rxt(k,175)*y(k,42) &
                      + rxt(k,240)*y(k,43) + 2.000_r8*rxt(k,243)*y(k,46) + rxt(k,176) &
                      *y(k,54) + rxt(k,245)*y(k,55) + rxt(k,177)*y(k,77) + rxt(k,178) &
                      *y(k,79) + rxt(k,197)*y(k,92) + rxt(k,179)*y(k,205)
         mat(k,1727) = rxt(k,194)*y(k,219)
         mat(k,1412) = rxt(k,177)*y(k,56)
         mat(k,609) = rxt(k,178)*y(k,56)
         mat(k,819) = mat(k,819) + rxt(k,197)*y(k,56)
         mat(k,2118) = rxt(k,179)*y(k,56)
         mat(k,1701) = mat(k,1701) + rxt(k,194)*y(k,59)
         mat(k,181) = -(rxt(k,266)*y(k,219) + rxt(k,274)*y(k,218))
         mat(k,1564) = -rxt(k,266)*y(k,86)
         mat(k,1510) = -rxt(k,274)*y(k,86)
         mat(k,931) = -(rxt(k,267)*y(k,219))
         mat(k,1659) = -rxt(k,267)*y(k,87)
         mat(k,998) = .050_r8*rxt(k,440)*y(k,135)
         mat(k,294) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,552) = .370_r8*rxt(k,279)*y(k,135)
         mat(k,1025) = .120_r8*rxt(k,308)*y(k,135)
         mat(k,837) = .110_r8*rxt(k,385)*y(k,135)
         mat(k,1209) = .330_r8*rxt(k,338)*y(k,135)
         mat(k,949) = .050_r8*rxt(k,443)*y(k,135)
         mat(k,1327) = .120_r8*rxt(k,352)*y(k,135)
         mat(k,1891) = rxt(k,270)*y(k,206)
         mat(k,2219) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25) &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,98) &
                      + .330_r8*rxt(k,338)*y(k,105) + .050_r8*rxt(k,443)*y(k,110) &
                      + .120_r8*rxt(k,352)*y(k,111)
         mat(k,2082) = rxt(k,268)*y(k,206)
         mat(k,444) = rxt(k,270)*y(k,124) + rxt(k,268)*y(k,205)
         mat(k,1659) = mat(k,1659) + .350_r8*rxt(k,277)*y(k,24)
         mat(k,1427) = rxt(k,232)*y(k,73)
         mat(k,865) = rxt(k,232)*y(k,54) + rxt(k,233)*y(k,77) + rxt(k,235)*y(k,89) &
                      + rxt(k,234)*y(k,231)
         mat(k,1401) = rxt(k,233)*y(k,73)
         mat(k,1460) = rxt(k,235)*y(k,73)
         mat(k,2264) = rxt(k,234)*y(k,73)
         mat(k,1465) = -(rxt(k,172)*y(k,219) + rxt(k,235)*y(k,73))
         mat(k,1690) = -rxt(k,172)*y(k,89)
         mat(k,870) = -rxt(k,235)*y(k,89)
         mat(k,1486) = rxt(k,256)*y(k,126)
         mat(k,1090) = rxt(k,288)*y(k,126)
         mat(k,1225) = rxt(k,314)*y(k,126)
         mat(k,902) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,286) = rxt(k,462)*y(k,126)
         mat(k,2130) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60)
         mat(k,1961) = rxt(k,171)*y(k,219)
         mat(k,1773) = rxt(k,256)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) &
                      + rxt(k,462)*y(k,67)
         mat(k,1690) = mat(k,1690) + rxt(k,171)*y(k,125)
         mat(k,423) = -(rxt(k,150)*y(k,219))
         mat(k,1602) = -rxt(k,150)*y(k,90)
         mat(k,1936) = rxt(k,169)*y(k,205)
         mat(k,2042) = rxt(k,169)*y(k,125)
         mat(k,780) = -(rxt(k,226)*y(k,134) + (rxt(k,530) + rxt(k,535)) * y(k,85))
         mat(k,2174) = -rxt(k,226)*y(k,91)
         mat(k,2125) = -(rxt(k,530) + rxt(k,535)) * y(k,91)
         mat(k,2149) = rxt(k,218)*y(k,205)
         mat(k,2071) = rxt(k,218)*y(k,19)
         mat(k,815) = -(rxt(k,197)*y(k,56) + rxt(k,198)*y(k,134) + rxt(k,199)*y(k,219) &
                      + (rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85))
         mat(k,1987) = -rxt(k,197)*y(k,92)
         mat(k,2176) = -rxt(k,198)*y(k,92)
         mat(k,1648) = -rxt(k,199)*y(k,92)
         mat(k,2126) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,92)
         mat(k,1710) = rxt(k,186)*y(k,205)
         mat(k,900) = rxt(k,191)*y(k,219)
         mat(k,2075) = rxt(k,186)*y(k,59)
         mat(k,1648) = mat(k,1648) + rxt(k,191)*y(k,60)
         mat(k,1155) = -(rxt(k,331)*y(k,219))
         mat(k,1674) = -rxt(k,331)*y(k,93)
         mat(k,586) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,538) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,1903) = rxt(k,330)*y(k,202) + rxt(k,337)*y(k,211)
         mat(k,569) = rxt(k,330)*y(k,124)
         mat(k,1311) = rxt(k,337)*y(k,124)
         mat(k,1674) = mat(k,1674) + .300_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100)
         mat(k,223) = -(rxt(k,362)*y(k,219))
         mat(k,1571) = -rxt(k,362)*y(k,94)
         mat(k,1142) = -(rxt(k,316)*y(k,219))
         mat(k,1673) = -rxt(k,316)*y(k,95)
         mat(k,585) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,537) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,576) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,1902) = .050_r8*rxt(k,374)*y(k,208) + .220_r8*rxt(k,336)*y(k,211) &
                      + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1757) = .050_r8*rxt(k,375)*y(k,208) + .220_r8*rxt(k,335)*y(k,211) &
                      + .250_r8*rxt(k,394)*y(k,227)
         mat(k,545) = .500_r8*rxt(k,320)*y(k,219)
         mat(k,1378) = .220_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1812) = .230_r8*rxt(k,333)*y(k,211) + .200_r8*rxt(k,321)*y(k,222) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,1285) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1310) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,199) + .230_r8*rxt(k,333)*y(k,200)
         mat(k,1673) = mat(k,1673) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + .500_r8*rxt(k,351)*y(k,109) + .500_r8*rxt(k,320) &
                      *y(k,147)
         mat(k,1132) = .200_r8*rxt(k,321)*y(k,200)
         mat(k,1174) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,348) = -(rxt(k,363)*y(k,219))
         mat(k,1591) = -rxt(k,363)*y(k,96)
         mat(k,1858) = .870_r8*rxt(k,374)*y(k,208)
         mat(k,1736) = .950_r8*rxt(k,375)*y(k,208)
         mat(k,1369) = rxt(k,370)*y(k,208)
         mat(k,1790) = .750_r8*rxt(k,371)*y(k,208)
         mat(k,1275) = .870_r8*rxt(k,374)*y(k,124) + .950_r8*rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,199) + .750_r8*rxt(k,371)*y(k,200)
         mat(k,137) = -(rxt(k,364)*y(k,219))
         mat(k,1560) = -rxt(k,364)*y(k,97)
         mat(k,736) = .600_r8*rxt(k,387)*y(k,219)
         mat(k,1560) = mat(k,1560) + .600_r8*rxt(k,387)*y(k,103)
         mat(k,836) = -(rxt(k,378)*y(k,126) + rxt(k,385)*y(k,135) + rxt(k,386) &
                      *y(k,219))
         mat(k,1739) = -rxt(k,378)*y(k,98)
         mat(k,2215) = -rxt(k,385)*y(k,98)
         mat(k,1650) = -rxt(k,386)*y(k,98)
         mat(k,583) = -(rxt(k,376)*y(k,219))
         mat(k,1623) = -rxt(k,376)*y(k,99)
         mat(k,1871) = .080_r8*rxt(k,368)*y(k,207)
         mat(k,1247) = .080_r8*rxt(k,368)*y(k,124)
         mat(k,534) = -(rxt(k,377)*y(k,219))
         mat(k,1618) = -rxt(k,377)*y(k,100)
         mat(k,1869) = .080_r8*rxt(k,374)*y(k,208)
         mat(k,1276) = .080_r8*rxt(k,374)*y(k,124)
         mat(k,1195) = -(rxt(k,379)*y(k,199) + rxt(k,380)*y(k,200) + rxt(k,381) &
                      *y(k,205) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126))
         mat(k,1380) = -rxt(k,379)*y(k,101)
         mat(k,1814) = -rxt(k,380)*y(k,101)
         mat(k,2094) = -rxt(k,381)*y(k,101)
         mat(k,1905) = -rxt(k,382)*y(k,101)
         mat(k,1760) = -rxt(k,383)*y(k,101)
         mat(k,840) = rxt(k,378)*y(k,126)
         mat(k,1760) = mat(k,1760) + rxt(k,378)*y(k,98)
         mat(k,405) = -(rxt(k,384)*y(k,219))
         mat(k,1599) = -rxt(k,384)*y(k,102)
         mat(k,1187) = rxt(k,381)*y(k,205)
         mat(k,2039) = rxt(k,381)*y(k,101)
         mat(k,737) = -(rxt(k,387)*y(k,219))
         mat(k,1640) = -rxt(k,387)*y(k,103)
         mat(k,2067) = rxt(k,367)*y(k,207) + rxt(k,372)*y(k,208)
         mat(k,1248) = rxt(k,367)*y(k,205)
         mat(k,1278) = rxt(k,372)*y(k,205)
         mat(k,76) = -(rxt(k,509)*y(k,219))
         mat(k,1552) = -rxt(k,509)*y(k,104)
         mat(k,1211) = -(rxt(k,338)*y(k,135) + rxt(k,339)*y(k,219))
         mat(k,2234) = -rxt(k,338)*y(k,105)
         mat(k,1677) = -rxt(k,339)*y(k,105)
         mat(k,841) = .300_r8*rxt(k,385)*y(k,135)
         mat(k,1906) = .360_r8*rxt(k,368)*y(k,207)
         mat(k,1761) = .400_r8*rxt(k,369)*y(k,207)
         mat(k,2234) = mat(k,2234) + .300_r8*rxt(k,385)*y(k,98)
         mat(k,1381) = .390_r8*rxt(k,365)*y(k,207)
         mat(k,1815) = .310_r8*rxt(k,366)*y(k,207)
         mat(k,1256) = .360_r8*rxt(k,368)*y(k,124) + .400_r8*rxt(k,369)*y(k,126) &
                      + .390_r8*rxt(k,365)*y(k,199) + .310_r8*rxt(k,366)*y(k,200)
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
         mat(k,312) = -(rxt(k,340)*y(k,219))
         mat(k,1586) = -rxt(k,340)*y(k,106)
         mat(k,2032) = rxt(k,334)*y(k,211)
         mat(k,1306) = rxt(k,334)*y(k,205)
         mat(k,508) = -(rxt(k,349)*y(k,219))
         mat(k,1614) = -rxt(k,349)*y(k,107)
         mat(k,1867) = .800_r8*rxt(k,358)*y(k,191)
         mat(k,884) = .800_r8*rxt(k,358)*y(k,124)
         mat(k,317) = -(rxt(k,350)*y(k,219))
         mat(k,1587) = -rxt(k,350)*y(k,108)
         mat(k,2033) = .800_r8*rxt(k,347)*y(k,215)
         mat(k,678) = .800_r8*rxt(k,347)*y(k,205)
         mat(k,574) = -(rxt(k,351)*y(k,219))
         mat(k,1622) = -rxt(k,351)*y(k,109)
         mat(k,1942) = rxt(k,354)*y(k,213)
         mat(k,1349) = rxt(k,354)*y(k,125)
         mat(k,950) = -(rxt(k,442)*y(k,126) + rxt(k,443)*y(k,135) + rxt(k,444) &
                      *y(k,219))
         mat(k,1744) = -rxt(k,442)*y(k,110)
         mat(k,2220) = -rxt(k,443)*y(k,110)
         mat(k,1660) = -rxt(k,444)*y(k,110)
         mat(k,1334) = -(rxt(k,352)*y(k,135) + rxt(k,353)*y(k,219))
         mat(k,2240) = -rxt(k,352)*y(k,111)
         mat(k,1683) = -rxt(k,353)*y(k,111)
         mat(k,844) = .200_r8*rxt(k,385)*y(k,135)
         mat(k,1911) = .560_r8*rxt(k,368)*y(k,207)
         mat(k,1767) = .600_r8*rxt(k,369)*y(k,207)
         mat(k,2240) = mat(k,2240) + .200_r8*rxt(k,385)*y(k,98)
         mat(k,1386) = .610_r8*rxt(k,365)*y(k,207)
         mat(k,1820) = .440_r8*rxt(k,366)*y(k,207)
         mat(k,1260) = .560_r8*rxt(k,368)*y(k,124) + .600_r8*rxt(k,369)*y(k,126) &
                      + .610_r8*rxt(k,365)*y(k,199) + .440_r8*rxt(k,366)*y(k,200)
         mat(k,493) = -(rxt(k,151)*y(k,124) + (rxt(k,152) + rxt(k,153) + rxt(k,154) &
                      ) * y(k,125) + rxt(k,163)*y(k,219))
         mat(k,1865) = -rxt(k,151)*y(k,112)
         mat(k,1938) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,112)
         mat(k,1612) = -rxt(k,163)*y(k,112)
         mat(k,185) = -((rxt(k,167) + rxt(k,168)) * y(k,218))
         mat(k,1511) = -(rxt(k,167) + rxt(k,168)) * y(k,113)
         mat(k,492) = rxt(k,152)*y(k,125)
         mat(k,1934) = rxt(k,152)*y(k,112)
         mat(k,1935) = rxt(k,170)*y(k,126)
         mat(k,1734) = rxt(k,170)*y(k,125)
         mat(k,375) = -(rxt(k,388)*y(k,219))
         mat(k,1595) = -rxt(k,388)*y(k,115)
         mat(k,1186) = .200_r8*rxt(k,380)*y(k,200)
         mat(k,1791) = .200_r8*rxt(k,380)*y(k,101)
         mat(k,1075) = -(rxt(k,389)*y(k,219))
         mat(k,1667) = -rxt(k,389)*y(k,116)
         mat(k,1191) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) + rxt(k,379)*y(k,199) &
                      + .800_r8*rxt(k,380)*y(k,200)
         mat(k,1896) = rxt(k,382)*y(k,101)
         mat(k,1751) = rxt(k,383)*y(k,101)
         mat(k,1375) = rxt(k,379)*y(k,101)
         mat(k,1806) = .800_r8*rxt(k,380)*y(k,101)
         mat(k,98) = -(rxt(k,479)*y(k,219))
         mat(k,1556) = -rxt(k,479)*y(k,120)
         mat(k,1924) = -(rxt(k,151)*y(k,112) + rxt(k,160)*y(k,126) + rxt(k,164) &
                      *y(k,205) + rxt(k,165)*y(k,135) + rxt(k,166)*y(k,134) + rxt(k,187) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,262)*y(k,200) + rxt(k,270) &
                      *y(k,206) + rxt(k,283)*y(k,196) + rxt(k,294)*y(k,199) + rxt(k,298) &
                      *y(k,204) + rxt(k,311)*y(k,197) + rxt(k,319)*y(k,221) + rxt(k,323) &
                      *y(k,222) + (rxt(k,329) + rxt(k,330)) * y(k,202) + (rxt(k,336) &
                      + rxt(k,337)) * y(k,211) + rxt(k,345)*y(k,213) + rxt(k,348) &
                      *y(k,215) + (rxt(k,358) + rxt(k,359)) * y(k,191) + rxt(k,368) &
                      *y(k,207) + rxt(k,374)*y(k,208) + rxt(k,382)*y(k,101) + rxt(k,393) &
                      *y(k,227) + rxt(k,397)*y(k,190) + rxt(k,400)*y(k,193) + rxt(k,405) &
                      *y(k,195) + rxt(k,407)*y(k,198) + rxt(k,411)*y(k,201) + rxt(k,414) &
                      *y(k,212) + rxt(k,417)*y(k,214) + rxt(k,420)*y(k,220) + rxt(k,427) &
                      *y(k,225) + rxt(k,433)*y(k,228) + rxt(k,436)*y(k,230) + rxt(k,447) &
                      *y(k,217) + rxt(k,452)*y(k,223) + rxt(k,457)*y(k,224))
         mat(k,497) = -rxt(k,151)*y(k,124)
         mat(k,1780) = -rxt(k,160)*y(k,124)
         mat(k,2114) = -rxt(k,164)*y(k,124)
         mat(k,2253) = -rxt(k,165)*y(k,124)
         mat(k,2192) = -rxt(k,166)*y(k,124)
         mat(k,1723) = -rxt(k,187)*y(k,124)
         mat(k,2161) = -rxt(k,219)*y(k,124)
         mat(k,1832) = -rxt(k,262)*y(k,124)
         mat(k,446) = -rxt(k,270)*y(k,124)
         mat(k,830) = -rxt(k,283)*y(k,124)
         mat(k,1395) = -rxt(k,294)*y(k,124)
         mat(k,712) = -rxt(k,298)*y(k,124)
         mat(k,795) = -rxt(k,311)*y(k,124)
         mat(k,775) = -rxt(k,319)*y(k,124)
         mat(k,1138) = -rxt(k,323)*y(k,124)
         mat(k,571) = -(rxt(k,329) + rxt(k,330)) * y(k,124)
         mat(k,1321) = -(rxt(k,336) + rxt(k,337)) * y(k,124)
         mat(k,1363) = -rxt(k,345)*y(k,124)
         mat(k,683) = -rxt(k,348)*y(k,124)
         mat(k,896) = -(rxt(k,358) + rxt(k,359)) * y(k,124)
         mat(k,1268) = -rxt(k,368)*y(k,124)
         mat(k,1300) = -rxt(k,374)*y(k,124)
         mat(k,1205) = -rxt(k,382)*y(k,124)
         mat(k,1182) = -rxt(k,393)*y(k,124)
         mat(k,523) = -rxt(k,397)*y(k,124)
         mat(k,489) = -rxt(k,400)*y(k,124)
         mat(k,440) = -rxt(k,405)*y(k,124)
         mat(k,629) = -rxt(k,407)*y(k,124)
         mat(k,766) = -rxt(k,411)*y(k,124)
         mat(k,718) = -rxt(k,414)*y(k,124)
         mat(k,881) = -rxt(k,417)*y(k,124)
         mat(k,453) = -rxt(k,420)*y(k,124)
         mat(k,733) = -rxt(k,427)*y(k,124)
         mat(k,758) = -rxt(k,433)*y(k,124)
         mat(k,505) = -rxt(k,436)*y(k,124)
         mat(k,1056) = -rxt(k,447)*y(k,124)
         mat(k,1119) = -rxt(k,452)*y(k,124)
         mat(k,921) = -rxt(k,457)*y(k,124)
         mat(k,497) = mat(k,497) + 2.000_r8*rxt(k,153)*y(k,125) + rxt(k,163)*y(k,219)
         mat(k,187) = 2.000_r8*rxt(k,167)*y(k,218)
         mat(k,1968) = 2.000_r8*rxt(k,153)*y(k,112) + rxt(k,156)*y(k,134) + rxt(k,472) &
                      *y(k,151)
         mat(k,2192) = mat(k,2192) + rxt(k,156)*y(k,125)
         mat(k,1239) = rxt(k,472)*y(k,125)
         mat(k,1533) = 2.000_r8*rxt(k,167)*y(k,113)
         mat(k,1697) = rxt(k,163)*y(k,112)
         mat(k,1969) = -((rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,112) + (rxt(k,156) &
                      + rxt(k,158)) * y(k,134) + rxt(k,157)*y(k,135) + rxt(k,169) &
                      *y(k,205) + rxt(k,170)*y(k,126) + rxt(k,171)*y(k,219) + rxt(k,189) &
                      *y(k,59) + rxt(k,220)*y(k,19) + rxt(k,305)*y(k,199) + rxt(k,354) &
                      *y(k,213) + rxt(k,412)*y(k,201) + rxt(k,415)*y(k,212) + rxt(k,418) &
                      *y(k,214) + rxt(k,422)*y(k,142) + rxt(k,425)*y(k,190) + rxt(k,472) &
                      *y(k,151))
         mat(k,498) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,125)
         mat(k,2193) = -(rxt(k,156) + rxt(k,158)) * y(k,125)
         mat(k,2254) = -rxt(k,157)*y(k,125)
         mat(k,2115) = -rxt(k,169)*y(k,125)
         mat(k,1781) = -rxt(k,170)*y(k,125)
         mat(k,1698) = -rxt(k,171)*y(k,125)
         mat(k,1724) = -rxt(k,189)*y(k,125)
         mat(k,2162) = -rxt(k,220)*y(k,125)
         mat(k,1396) = -rxt(k,305)*y(k,125)
         mat(k,1364) = -rxt(k,354)*y(k,125)
         mat(k,767) = -rxt(k,412)*y(k,125)
         mat(k,719) = -rxt(k,415)*y(k,125)
         mat(k,882) = -rxt(k,418)*y(k,125)
         mat(k,467) = -rxt(k,422)*y(k,125)
         mat(k,524) = -rxt(k,425)*y(k,125)
         mat(k,1240) = -rxt(k,472)*y(k,125)
         mat(k,642) = rxt(k,356)*y(k,219)
         mat(k,358) = rxt(k,327)*y(k,126)
         mat(k,2162) = mat(k,2162) + rxt(k,219)*y(k,124)
         mat(k,1724) = mat(k,1724) + rxt(k,187)*y(k,124)
         mat(k,426) = rxt(k,150)*y(k,219)
         mat(k,590) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,1206) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126)
         mat(k,1925) = rxt(k,219)*y(k,19) + rxt(k,187)*y(k,59) + rxt(k,382)*y(k,101) &
                      + 2.000_r8*rxt(k,160)*y(k,126) + rxt(k,166)*y(k,134) &
                      + rxt(k,165)*y(k,135) + rxt(k,397)*y(k,190) + rxt(k,358) &
                      *y(k,191) + rxt(k,400)*y(k,193) + rxt(k,405)*y(k,195) &
                      + rxt(k,283)*y(k,196) + rxt(k,311)*y(k,197) + rxt(k,407) &
                      *y(k,198) + rxt(k,294)*y(k,199) + rxt(k,262)*y(k,200) &
                      + rxt(k,411)*y(k,201) + rxt(k,329)*y(k,202) + rxt(k,298) &
                      *y(k,204) + rxt(k,164)*y(k,205) + rxt(k,270)*y(k,206) &
                      + .920_r8*rxt(k,368)*y(k,207) + .920_r8*rxt(k,374)*y(k,208) &
                      + rxt(k,336)*y(k,211) + rxt(k,414)*y(k,212) + rxt(k,345) &
                      *y(k,213) + rxt(k,417)*y(k,214) + rxt(k,348)*y(k,215) &
                      + 1.600_r8*rxt(k,447)*y(k,217) + rxt(k,420)*y(k,220) &
                      + rxt(k,319)*y(k,221) + rxt(k,323)*y(k,222) + .900_r8*rxt(k,452) &
                      *y(k,223) + .800_r8*rxt(k,457)*y(k,224) + rxt(k,427)*y(k,225) &
                      + rxt(k,393)*y(k,227) + rxt(k,433)*y(k,228) + rxt(k,436) &
                      *y(k,230)
         mat(k,1781) = mat(k,1781) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,101) &
                      + 2.000_r8*rxt(k,160)*y(k,124) + rxt(k,161)*y(k,134) &
                      + rxt(k,159)*y(k,205) + rxt(k,369)*y(k,207) + rxt(k,375) &
                      *y(k,208) + rxt(k,335)*y(k,211) + rxt(k,346)*y(k,213) &
                      + 2.000_r8*rxt(k,448)*y(k,217) + rxt(k,162)*y(k,219) &
                      + rxt(k,394)*y(k,227)
         mat(k,864) = rxt(k,317)*y(k,219)
         mat(k,2193) = mat(k,2193) + rxt(k,166)*y(k,124) + rxt(k,161)*y(k,126)
         mat(k,2254) = mat(k,2254) + rxt(k,165)*y(k,124)
         mat(k,623) = rxt(k,454)*y(k,219)
         mat(k,524) = mat(k,524) + rxt(k,397)*y(k,124)
         mat(k,897) = rxt(k,358)*y(k,124)
         mat(k,490) = rxt(k,400)*y(k,124)
         mat(k,441) = rxt(k,405)*y(k,124)
         mat(k,831) = rxt(k,283)*y(k,124)
         mat(k,796) = rxt(k,311)*y(k,124)
         mat(k,630) = rxt(k,407)*y(k,124)
         mat(k,1396) = mat(k,1396) + rxt(k,294)*y(k,124)
         mat(k,1833) = rxt(k,262)*y(k,124) + .500_r8*rxt(k,445)*y(k,217)
         mat(k,767) = mat(k,767) + rxt(k,411)*y(k,124)
         mat(k,572) = rxt(k,329)*y(k,124)
         mat(k,713) = rxt(k,298)*y(k,124)
         mat(k,2115) = mat(k,2115) + rxt(k,164)*y(k,124) + rxt(k,159)*y(k,126)
         mat(k,447) = rxt(k,270)*y(k,124)
         mat(k,1269) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126)
         mat(k,1301) = .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,1322) = rxt(k,336)*y(k,124) + rxt(k,335)*y(k,126)
         mat(k,719) = mat(k,719) + rxt(k,414)*y(k,124)
         mat(k,1364) = mat(k,1364) + rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126)
         mat(k,882) = mat(k,882) + rxt(k,417)*y(k,124)
         mat(k,684) = rxt(k,348)*y(k,124)
         mat(k,1057) = 1.600_r8*rxt(k,447)*y(k,124) + 2.000_r8*rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1698) = mat(k,1698) + rxt(k,356)*y(k,1) + rxt(k,150)*y(k,90) &
                      + .700_r8*rxt(k,376)*y(k,99) + rxt(k,162)*y(k,126) + rxt(k,317) &
                      *y(k,127) + rxt(k,454)*y(k,177)
         mat(k,454) = rxt(k,420)*y(k,124)
         mat(k,776) = rxt(k,319)*y(k,124)
         mat(k,1139) = rxt(k,323)*y(k,124)
         mat(k,1120) = .900_r8*rxt(k,452)*y(k,124)
         mat(k,922) = .800_r8*rxt(k,457)*y(k,124)
         mat(k,734) = rxt(k,427)*y(k,124)
         mat(k,1183) = rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126)
         mat(k,759) = rxt(k,433)*y(k,124)
         mat(k,506) = rxt(k,436)*y(k,124)
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
         mat(k,1778) = -(rxt(k,159)*y(k,205) + rxt(k,160)*y(k,124) + rxt(k,161) &
                      *y(k,134) + rxt(k,162)*y(k,219) + rxt(k,170)*y(k,125) + rxt(k,256) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,211) + rxt(k,346) &
                      *y(k,213) + rxt(k,369)*y(k,207) + rxt(k,375)*y(k,208) + rxt(k,378) &
                      *y(k,98) + rxt(k,383)*y(k,101) + rxt(k,394)*y(k,227) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,110) + rxt(k,448)*y(k,217) + rxt(k,459) &
                      *y(k,179) + rxt(k,462)*y(k,67))
         mat(k,2112) = -rxt(k,159)*y(k,126)
         mat(k,1922) = -rxt(k,160)*y(k,126)
         mat(k,2190) = -rxt(k,161)*y(k,126)
         mat(k,1695) = -rxt(k,162)*y(k,126)
         mat(k,1966) = -rxt(k,170)*y(k,126)
         mat(k,1490) = -rxt(k,256)*y(k,126)
         mat(k,1092) = -rxt(k,288)*y(k,126)
         mat(k,1035) = -rxt(k,307)*y(k,126)
         mat(k,1227) = -rxt(k,314)*y(k,126)
         mat(k,357) = -rxt(k,327)*y(k,126)
         mat(k,1319) = -rxt(k,335)*y(k,126)
         mat(k,1361) = -rxt(k,346)*y(k,126)
         mat(k,1266) = -rxt(k,369)*y(k,126)
         mat(k,1298) = -rxt(k,375)*y(k,126)
         mat(k,848) = -rxt(k,378)*y(k,126)
         mat(k,1203) = -rxt(k,383)*y(k,126)
         mat(k,1180) = -rxt(k,394)*y(k,126)
         mat(k,1013) = -rxt(k,439)*y(k,126)
         mat(k,963) = -rxt(k,442)*y(k,126)
         mat(k,1054) = -rxt(k,448)*y(k,126)
         mat(k,977) = -rxt(k,459)*y(k,126)
         mat(k,288) = -rxt(k,462)*y(k,126)
         mat(k,562) = rxt(k,221)*y(k,134)
         mat(k,2005) = rxt(k,188)*y(k,60)
         mat(k,905) = rxt(k,188)*y(k,56) + rxt(k,190)*y(k,134) + rxt(k,191)*y(k,219)
         mat(k,872) = rxt(k,235)*y(k,89)
         mat(k,1469) = rxt(k,235)*y(k,73) + rxt(k,172)*y(k,219)
         mat(k,580) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,1966) = mat(k,1966) + rxt(k,158)*y(k,134) + rxt(k,157)*y(k,135)
         mat(k,2190) = mat(k,2190) + rxt(k,221)*y(k,20) + rxt(k,190)*y(k,60) &
                      + rxt(k,158)*y(k,125)
         mat(k,2251) = rxt(k,157)*y(k,125)
         mat(k,531) = rxt(k,303)*y(k,219)
         mat(k,1695) = mat(k,1695) + rxt(k,191)*y(k,60) + rxt(k,172)*y(k,89) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,140)
         mat(k,859) = -(rxt(k,317)*y(k,219))
         mat(k,1652) = -rxt(k,317)*y(k,127)
         mat(k,1024) = rxt(k,307)*y(k,126)
         mat(k,535) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,407) = rxt(k,384)*y(k,219)
         mat(k,376) = rxt(k,388)*y(k,219)
         mat(k,1072) = rxt(k,389)*y(k,219)
         mat(k,1741) = rxt(k,307)*y(k,29)
         mat(k,1652) = mat(k,1652) + .500_r8*rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) &
                      + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116)
         mat(k,393) = -(rxt(k,449)*y(k,219))
         mat(k,1597) = -rxt(k,449)*y(k,128)
         mat(k,2037) = rxt(k,446)*y(k,217)
         mat(k,1043) = rxt(k,446)*y(k,205)
         mat(k,2198) = -(rxt(k,130)*y(k,135) + 4._r8*rxt(k,131)*y(k,134) + rxt(k,133) &
                      *y(k,77) + rxt(k,134)*y(k,79) + rxt(k,139)*y(k,205) + rxt(k,145) &
                      *y(k,219) + (rxt(k,156) + rxt(k,158)) * y(k,125) + rxt(k,161) &
                      *y(k,126) + rxt(k,166)*y(k,124) + rxt(k,190)*y(k,60) + rxt(k,192) &
                      *y(k,59) + rxt(k,195)*y(k,85) + rxt(k,198)*y(k,92) + rxt(k,221) &
                      *y(k,20) + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81) + rxt(k,226) &
                      *y(k,91) + rxt(k,257)*y(k,42) + rxt(k,464)*y(k,138))
         mat(k,2259) = -rxt(k,130)*y(k,134)
         mat(k,1413) = -rxt(k,133)*y(k,134)
         mat(k,610) = -rxt(k,134)*y(k,134)
         mat(k,2120) = -rxt(k,139)*y(k,134)
         mat(k,1703) = -rxt(k,145)*y(k,134)
         mat(k,1974) = -(rxt(k,156) + rxt(k,158)) * y(k,134)
         mat(k,1786) = -rxt(k,161)*y(k,134)
         mat(k,1930) = -rxt(k,166)*y(k,134)
         mat(k,909) = -rxt(k,190)*y(k,134)
         mat(k,1729) = -rxt(k,192)*y(k,134)
         mat(k,2143) = -rxt(k,195)*y(k,134)
         mat(k,820) = -rxt(k,198)*y(k,134)
         mat(k,565) = -rxt(k,221)*y(k,134)
         mat(k,2167) = -rxt(k,222)*y(k,134)
         mat(k,812) = -rxt(k,224)*y(k,134)
         mat(k,786) = -rxt(k,226)*y(k,134)
         mat(k,1498) = -rxt(k,257)*y(k,134)
         mat(k,366) = -rxt(k,464)*y(k,134)
         mat(k,1457) = rxt(k,137)*y(k,205)
         mat(k,499) = rxt(k,151)*y(k,124) + rxt(k,152)*y(k,125)
         mat(k,1930) = mat(k,1930) + rxt(k,151)*y(k,112)
         mat(k,1974) = mat(k,1974) + rxt(k,152)*y(k,112)
         mat(k,2120) = mat(k,2120) + rxt(k,137)*y(k,76)
         mat(k,1703) = mat(k,1703) + 2.000_r8*rxt(k,147)*y(k,219)
         mat(k,2260) = -(rxt(k,129)*y(k,218) + rxt(k,130)*y(k,134) + rxt(k,140) &
                      *y(k,205) + rxt(k,141)*y(k,76) + rxt(k,146)*y(k,219) + rxt(k,157) &
                      *y(k,125) + rxt(k,165)*y(k,124) + rxt(k,181)*y(k,56) + rxt(k,213) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,105) + rxt(k,352)*y(k,111) + rxt(k,385)*y(k,98) + rxt(k,423) &
                      *y(k,142) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110) + rxt(k,468) &
                      *y(k,149) + rxt(k,474)*y(k,151))
         mat(k,1540) = -rxt(k,129)*y(k,135)
         mat(k,2199) = -rxt(k,130)*y(k,135)
         mat(k,2121) = -rxt(k,140)*y(k,135)
         mat(k,1458) = -rxt(k,141)*y(k,135)
         mat(k,1704) = -rxt(k,146)*y(k,135)
         mat(k,1975) = -rxt(k,157)*y(k,135)
         mat(k,1931) = -rxt(k,165)*y(k,135)
         mat(k,2014) = -rxt(k,181)*y(k,135)
         mat(k,1425) = -rxt(k,213)*y(k,135)
         mat(k,557) = -rxt(k,279)*y(k,135)
         mat(k,1041) = -rxt(k,308)*y(k,135)
         mat(k,1220) = -rxt(k,338)*y(k,135)
         mat(k,1347) = -rxt(k,352)*y(k,135)
         mat(k,851) = -rxt(k,385)*y(k,135)
         mat(k,468) = -rxt(k,423)*y(k,135)
         mat(k,1019) = -rxt(k,440)*y(k,135)
         mat(k,969) = -rxt(k,443)*y(k,135)
         mat(k,518) = -rxt(k,468)*y(k,135)
         mat(k,1245) = -rxt(k,474)*y(k,135)
         mat(k,1399) = .150_r8*rxt(k,293)*y(k,205)
         mat(k,2121) = mat(k,2121) + .150_r8*rxt(k,293)*y(k,199) + .150_r8*rxt(k,343) &
                      *y(k,213)
         mat(k,1367) = .150_r8*rxt(k,343)*y(k,205)
         mat(k,327) = -(rxt(k,475)*y(k,151))
         mat(k,1231) = -rxt(k,475)*y(k,137)
         mat(k,2147) = rxt(k,215)*y(k,59)
         mat(k,1709) = rxt(k,215)*y(k,19) + 2.000_r8*rxt(k,185)*y(k,59)
         mat(k,359) = -(rxt(k,464)*y(k,134) + rxt(k,465)*y(k,219))
         mat(k,2170) = -rxt(k,464)*y(k,138)
         mat(k,1593) = -rxt(k,465)*y(k,138)
         mat(k,1148) = rxt(k,331)*y(k,219)
         mat(k,1853) = .100_r8*rxt(k,452)*y(k,223)
         mat(k,1573) = rxt(k,331)*y(k,93)
         mat(k,1104) = .100_r8*rxt(k,452)*y(k,124)
         mat(k,526) = -(rxt(k,303)*y(k,219))
         mat(k,1617) = -rxt(k,303)*y(k,140)
         mat(k,1940) = rxt(k,305)*y(k,199)
         mat(k,1370) = rxt(k,305)*y(k,125)
         mat(k,1933) = rxt(k,425)*y(k,190)
         mat(k,519) = rxt(k,425)*y(k,125)
         mat(k,465) = -(rxt(k,422)*y(k,125) + rxt(k,423)*y(k,135))
         mat(k,1937) = -rxt(k,422)*y(k,142)
         mat(k,2208) = -rxt(k,423)*y(k,142)
         mat(k,198) = .070_r8*rxt(k,409)*y(k,219)
         mat(k,1863) = rxt(k,407)*y(k,198)
         mat(k,176) = .060_r8*rxt(k,421)*y(k,219)
         mat(k,219) = .070_r8*rxt(k,437)*y(k,219)
         mat(k,626) = rxt(k,407)*y(k,124)
         mat(k,1608) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,143) &
                      + .070_r8*rxt(k,437)*y(k,186)
         mat(k,174) = -(rxt(k,421)*y(k,219))
         mat(k,1563) = -rxt(k,421)*y(k,143)
         mat(k,166) = .530_r8*rxt(k,398)*y(k,219)
         mat(k,1563) = mat(k,1563) + .530_r8*rxt(k,398)*y(k,7)
         mat(k,332) = -(rxt(k,424)*y(k,219))
         mat(k,1588) = -rxt(k,424)*y(k,144)
         mat(k,2034) = rxt(k,419)*y(k,220)
         mat(k,449) = rxt(k,419)*y(k,205)
         mat(k,542) = -(rxt(k,320)*y(k,219))
         mat(k,1619) = -rxt(k,320)*y(k,147)
         mat(k,2054) = rxt(k,318)*y(k,221)
         mat(k,769) = rxt(k,318)*y(k,205)
         mat(k,399) = -(rxt(k,324)*y(k,219))
         mat(k,1598) = -rxt(k,324)*y(k,148)
         mat(k,2038) = .850_r8*rxt(k,322)*y(k,222)
         mat(k,1129) = .850_r8*rxt(k,322)*y(k,205)
         mat(k,513) = -(rxt(k,468)*y(k,135) + rxt(k,471)*y(k,219))
         mat(k,2209) = -rxt(k,468)*y(k,149)
         mat(k,1615) = -rxt(k,471)*y(k,149)
         mat(k,1234) = -(rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,125) &
                      + rxt(k,474)*y(k,135) + rxt(k,475)*y(k,137) + rxt(k,476) &
                      *y(k,219))
         mat(k,2151) = -rxt(k,469)*y(k,151)
         mat(k,1713) = -rxt(k,470)*y(k,151)
         mat(k,1955) = -rxt(k,472)*y(k,151)
         mat(k,2236) = -rxt(k,474)*y(k,151)
         mat(k,329) = -rxt(k,475)*y(k,151)
         mat(k,1679) = -rxt(k,476)*y(k,151)
         mat(k,2180) = rxt(k,464)*y(k,138)
         mat(k,2236) = mat(k,2236) + rxt(k,468)*y(k,149)
         mat(k,363) = rxt(k,464)*y(k,134)
         mat(k,514) = rxt(k,468)*y(k,135) + rxt(k,471)*y(k,219)
         mat(k,1679) = mat(k,1679) + rxt(k,471)*y(k,149)
         mat(k,853) = -(rxt(k,467)*y(k,219))
         mat(k,1651) = -rxt(k,467)*y(k,152)
         mat(k,2150) = rxt(k,469)*y(k,151)
         mat(k,1711) = rxt(k,470)*y(k,151)
         mat(k,285) = rxt(k,462)*y(k,126) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,219)
         mat(k,1948) = rxt(k,472)*y(k,151)
         mat(k,1740) = rxt(k,462)*y(k,67)
         mat(k,2216) = rxt(k,474)*y(k,151)
         mat(k,328) = rxt(k,475)*y(k,151)
         mat(k,361) = rxt(k,465)*y(k,219)
         mat(k,1233) = rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,125) &
                      + rxt(k,474)*y(k,135) + rxt(k,475)*y(k,137) + rxt(k,476) &
                      *y(k,219)
         mat(k,1651) = mat(k,1651) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,67) &
                      + rxt(k,465)*y(k,138) + rxt(k,476)*y(k,151)
         mat(k,259) = -(rxt(k,478)*y(k,231))
         mat(k,2263) = -rxt(k,478)*y(k,153)
         mat(k,852) = rxt(k,467)*y(k,219)
         mat(k,1579) = rxt(k,467)*y(k,152)
         mat(k,986) = .2202005_r8*rxt(k,497)*y(k,135)
         mat(k,937) = .0508005_r8*rxt(k,513)*y(k,135)
         mat(k,1841) = .1279005_r8*rxt(k,496)*y(k,192) + .0097005_r8*rxt(k,501) &
                      *y(k,194) + .0003005_r8*rxt(k,504)*y(k,209) &
                      + .1056005_r8*rxt(k,508)*y(k,210) + .0245005_r8*rxt(k,512) &
                      *y(k,216) + .0154005_r8*rxt(k,518)*y(k,226) &
                      + .0063005_r8*rxt(k,521)*y(k,229)
         mat(k,2201) = .2202005_r8*rxt(k,497)*y(k,6) + .0508005_r8*rxt(k,513)*y(k,110)
         mat(k,45) = .5931005_r8*rxt(k,515)*y(k,219)
         mat(k,51) = .1279005_r8*rxt(k,496)*y(k,124) + .2202005_r8*rxt(k,495)*y(k,205)
         mat(k,57) = .0097005_r8*rxt(k,501)*y(k,124) + .0023005_r8*rxt(k,500)*y(k,205)
         mat(k,2016) = .2202005_r8*rxt(k,495)*y(k,192) + .0023005_r8*rxt(k,500) &
                      *y(k,194) + .0031005_r8*rxt(k,503)*y(k,209) &
                      + .2381005_r8*rxt(k,507)*y(k,210) + .0508005_r8*rxt(k,511) &
                      *y(k,216) + .1364005_r8*rxt(k,517)*y(k,226) &
                      + .1677005_r8*rxt(k,520)*y(k,229)
         mat(k,63) = .0003005_r8*rxt(k,504)*y(k,124) + .0031005_r8*rxt(k,503)*y(k,205)
         mat(k,69) = .1056005_r8*rxt(k,508)*y(k,124) + .2381005_r8*rxt(k,507)*y(k,205)
         mat(k,77) = .0245005_r8*rxt(k,512)*y(k,124) + .0508005_r8*rxt(k,511)*y(k,205)
         mat(k,1542) = .5931005_r8*rxt(k,515)*y(k,174)
         mat(k,83) = .0154005_r8*rxt(k,518)*y(k,124) + .1364005_r8*rxt(k,517)*y(k,205)
         mat(k,89) = .0063005_r8*rxt(k,521)*y(k,124) + .1677005_r8*rxt(k,520)*y(k,205)
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
         mat(k,987) = .2067005_r8*rxt(k,497)*y(k,135)
         mat(k,938) = .1149005_r8*rxt(k,513)*y(k,135)
         mat(k,1842) = .1792005_r8*rxt(k,496)*y(k,192) + .0034005_r8*rxt(k,501) &
                      *y(k,194) + .0003005_r8*rxt(k,504)*y(k,209) &
                      + .1026005_r8*rxt(k,508)*y(k,210) + .0082005_r8*rxt(k,512) &
                      *y(k,216) + .0452005_r8*rxt(k,518)*y(k,226) &
                      + .0237005_r8*rxt(k,521)*y(k,229)
         mat(k,2202) = .2067005_r8*rxt(k,497)*y(k,6) + .1149005_r8*rxt(k,513)*y(k,110)
         mat(k,46) = .1534005_r8*rxt(k,515)*y(k,219)
         mat(k,52) = .1792005_r8*rxt(k,496)*y(k,124) + .2067005_r8*rxt(k,495)*y(k,205)
         mat(k,58) = .0034005_r8*rxt(k,501)*y(k,124) + .0008005_r8*rxt(k,500)*y(k,205)
         mat(k,2017) = .2067005_r8*rxt(k,495)*y(k,192) + .0008005_r8*rxt(k,500) &
                      *y(k,194) + .0035005_r8*rxt(k,503)*y(k,209) &
                      + .1308005_r8*rxt(k,507)*y(k,210) + .1149005_r8*rxt(k,511) &
                      *y(k,216) + .0101005_r8*rxt(k,517)*y(k,226) &
                      + .0174005_r8*rxt(k,520)*y(k,229)
         mat(k,64) = .0003005_r8*rxt(k,504)*y(k,124) + .0035005_r8*rxt(k,503)*y(k,205)
         mat(k,70) = .1026005_r8*rxt(k,508)*y(k,124) + .1308005_r8*rxt(k,507)*y(k,205)
         mat(k,78) = .0082005_r8*rxt(k,512)*y(k,124) + .1149005_r8*rxt(k,511)*y(k,205)
         mat(k,1543) = .1534005_r8*rxt(k,515)*y(k,174)
         mat(k,84) = .0452005_r8*rxt(k,518)*y(k,124) + .0101005_r8*rxt(k,517)*y(k,205)
         mat(k,90) = .0237005_r8*rxt(k,521)*y(k,124) + .0174005_r8*rxt(k,520)*y(k,205)
         mat(k,988) = .0653005_r8*rxt(k,497)*y(k,135)
         mat(k,939) = .0348005_r8*rxt(k,513)*y(k,135)
         mat(k,1843) = .0676005_r8*rxt(k,496)*y(k,192) + .1579005_r8*rxt(k,501) &
                      *y(k,194) + .0073005_r8*rxt(k,504)*y(k,209) &
                      + .0521005_r8*rxt(k,508)*y(k,210) + .0772005_r8*rxt(k,512) &
                      *y(k,216) + .0966005_r8*rxt(k,518)*y(k,226) &
                      + .0025005_r8*rxt(k,521)*y(k,229)
         mat(k,2203) = .0653005_r8*rxt(k,497)*y(k,6) + .0348005_r8*rxt(k,513)*y(k,110)
         mat(k,47) = .0459005_r8*rxt(k,515)*y(k,219)
         mat(k,53) = .0676005_r8*rxt(k,496)*y(k,124) + .0653005_r8*rxt(k,495)*y(k,205)
         mat(k,59) = .1579005_r8*rxt(k,501)*y(k,124) + .0843005_r8*rxt(k,500)*y(k,205)
         mat(k,2018) = .0653005_r8*rxt(k,495)*y(k,192) + .0843005_r8*rxt(k,500) &
                      *y(k,194) + .0003005_r8*rxt(k,503)*y(k,209) &
                      + .0348005_r8*rxt(k,507)*y(k,210) + .0348005_r8*rxt(k,511) &
                      *y(k,216) + .0763005_r8*rxt(k,517)*y(k,226) + .086_r8*rxt(k,520) &
                      *y(k,229)
         mat(k,65) = .0073005_r8*rxt(k,504)*y(k,124) + .0003005_r8*rxt(k,503)*y(k,205)
         mat(k,71) = .0521005_r8*rxt(k,508)*y(k,124) + .0348005_r8*rxt(k,507)*y(k,205)
         mat(k,79) = .0772005_r8*rxt(k,512)*y(k,124) + .0348005_r8*rxt(k,511)*y(k,205)
         mat(k,1544) = .0459005_r8*rxt(k,515)*y(k,174)
         mat(k,85) = .0966005_r8*rxt(k,518)*y(k,124) + .0763005_r8*rxt(k,517)*y(k,205)
         mat(k,91) = .0025005_r8*rxt(k,521)*y(k,124) + .086_r8*rxt(k,520)*y(k,205)
         mat(k,989) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,497) &
                      *y(k,135)
         mat(k,833) = .0590245_r8*rxt(k,502)*y(k,126) + .0033005_r8*rxt(k,505) &
                      *y(k,135)
         mat(k,940) = .1749305_r8*rxt(k,510)*y(k,126) + .0554005_r8*rxt(k,513) &
                      *y(k,135)
         mat(k,1844) = .079_r8*rxt(k,496)*y(k,192) + .0059005_r8*rxt(k,501)*y(k,194) &
                      + .0057005_r8*rxt(k,504)*y(k,209) + .0143005_r8*rxt(k,508) &
                      *y(k,210) + .0332005_r8*rxt(k,512)*y(k,216) &
                      + .0073005_r8*rxt(k,518)*y(k,226) + .011_r8*rxt(k,521)*y(k,229)
         mat(k,1732) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,502)*y(k,98) &
                      + .1749305_r8*rxt(k,510)*y(k,110)
         mat(k,2204) = .1284005_r8*rxt(k,497)*y(k,6) + .0033005_r8*rxt(k,505)*y(k,98) &
                      + .0554005_r8*rxt(k,513)*y(k,110)
         mat(k,48) = .0085005_r8*rxt(k,515)*y(k,219)
         mat(k,54) = .079_r8*rxt(k,496)*y(k,124) + .1284005_r8*rxt(k,495)*y(k,205)
         mat(k,60) = .0059005_r8*rxt(k,501)*y(k,124) + .0443005_r8*rxt(k,500)*y(k,205)
         mat(k,2019) = .1284005_r8*rxt(k,495)*y(k,192) + .0443005_r8*rxt(k,500) &
                      *y(k,194) + .0271005_r8*rxt(k,503)*y(k,209) &
                      + .0076005_r8*rxt(k,507)*y(k,210) + .0554005_r8*rxt(k,511) &
                      *y(k,216) + .2157005_r8*rxt(k,517)*y(k,226) &
                      + .0512005_r8*rxt(k,520)*y(k,229)
         mat(k,66) = .0057005_r8*rxt(k,504)*y(k,124) + .0271005_r8*rxt(k,503)*y(k,205)
         mat(k,72) = .0143005_r8*rxt(k,508)*y(k,124) + .0076005_r8*rxt(k,507)*y(k,205)
         mat(k,80) = .0332005_r8*rxt(k,512)*y(k,124) + .0554005_r8*rxt(k,511)*y(k,205)
         mat(k,1545) = .0085005_r8*rxt(k,515)*y(k,174)
         mat(k,86) = .0073005_r8*rxt(k,518)*y(k,124) + .2157005_r8*rxt(k,517)*y(k,205)
         mat(k,92) = .011_r8*rxt(k,521)*y(k,124) + .0512005_r8*rxt(k,520)*y(k,205)
         mat(k,990) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,497)*y(k,135)
         mat(k,834) = .0250245_r8*rxt(k,502)*y(k,126)
         mat(k,941) = .5901905_r8*rxt(k,510)*y(k,126) + .1278005_r8*rxt(k,513) &
                      *y(k,135)
         mat(k,1845) = .1254005_r8*rxt(k,496)*y(k,192) + .0536005_r8*rxt(k,501) &
                      *y(k,194) + .0623005_r8*rxt(k,504)*y(k,209) &
                      + .0166005_r8*rxt(k,508)*y(k,210) + .130_r8*rxt(k,512)*y(k,216) &
                      + .238_r8*rxt(k,518)*y(k,226) + .1185005_r8*rxt(k,521)*y(k,229)
         mat(k,1733) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,502)*y(k,98) &
                      + .5901905_r8*rxt(k,510)*y(k,110)
         mat(k,2205) = .114_r8*rxt(k,497)*y(k,6) + .1278005_r8*rxt(k,513)*y(k,110)
         mat(k,49) = .0128005_r8*rxt(k,515)*y(k,219)
         mat(k,55) = .1254005_r8*rxt(k,496)*y(k,124) + .114_r8*rxt(k,495)*y(k,205)
         mat(k,61) = .0536005_r8*rxt(k,501)*y(k,124) + .1621005_r8*rxt(k,500)*y(k,205)
         mat(k,2020) = .114_r8*rxt(k,495)*y(k,192) + .1621005_r8*rxt(k,500)*y(k,194) &
                      + .0474005_r8*rxt(k,503)*y(k,209) + .0113005_r8*rxt(k,507) &
                      *y(k,210) + .1278005_r8*rxt(k,511)*y(k,216) &
                      + .0738005_r8*rxt(k,517)*y(k,226) + .1598005_r8*rxt(k,520) &
                      *y(k,229)
         mat(k,67) = .0623005_r8*rxt(k,504)*y(k,124) + .0474005_r8*rxt(k,503)*y(k,205)
         mat(k,73) = .0166005_r8*rxt(k,508)*y(k,124) + .0113005_r8*rxt(k,507)*y(k,205)
         mat(k,81) = .130_r8*rxt(k,512)*y(k,124) + .1278005_r8*rxt(k,511)*y(k,205)
         mat(k,1546) = .0128005_r8*rxt(k,515)*y(k,174)
         mat(k,87) = .238_r8*rxt(k,518)*y(k,124) + .0738005_r8*rxt(k,517)*y(k,205)
         mat(k,93) = .1185005_r8*rxt(k,521)*y(k,124) + .1598005_r8*rxt(k,520)*y(k,205)
         mat(k,50) = -(rxt(k,515)*y(k,219))
         mat(k,1547) = -rxt(k,515)*y(k,174)
         mat(k,191) = .100_r8*rxt(k,429)*y(k,219)
         mat(k,209) = .230_r8*rxt(k,431)*y(k,219)
         mat(k,1567) = .100_r8*rxt(k,429)*y(k,182) + .230_r8*rxt(k,431)*y(k,184)
         mat(k,644) = -(rxt(k,453)*y(k,219))
         mat(k,1631) = -rxt(k,453)*y(k,176)
         mat(k,2058) = rxt(k,451)*y(k,223)
         mat(k,1105) = rxt(k,451)*y(k,205)
         mat(k,619) = -(rxt(k,454)*y(k,219))
         mat(k,1628) = -rxt(k,454)*y(k,177)
         mat(k,1873) = .200_r8*rxt(k,447)*y(k,217) + .200_r8*rxt(k,457)*y(k,224)
         mat(k,1793) = .500_r8*rxt(k,445)*y(k,217)
         mat(k,1044) = .200_r8*rxt(k,447)*y(k,124) + .500_r8*rxt(k,445)*y(k,200)
         mat(k,912) = .200_r8*rxt(k,457)*y(k,124)
         mat(k,476) = -(rxt(k,458)*y(k,219))
         mat(k,1610) = -rxt(k,458)*y(k,178)
         mat(k,2050) = rxt(k,456)*y(k,224)
         mat(k,911) = rxt(k,456)*y(k,205)
         mat(k,971) = -(rxt(k,459)*y(k,126) + rxt(k,460)*y(k,219))
         mat(k,1745) = -rxt(k,459)*y(k,179)
         mat(k,1661) = -rxt(k,460)*y(k,179)
         mat(k,999) = .330_r8*rxt(k,440)*y(k,135)
         mat(k,951) = .330_r8*rxt(k,443)*y(k,135)
         mat(k,1892) = .800_r8*rxt(k,447)*y(k,217) + .800_r8*rxt(k,457)*y(k,224)
         mat(k,1745) = mat(k,1745) + rxt(k,448)*y(k,217)
         mat(k,2221) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,110)
         mat(k,620) = rxt(k,454)*y(k,219)
         mat(k,1802) = .500_r8*rxt(k,445)*y(k,217) + rxt(k,455)*y(k,224)
         mat(k,1046) = .800_r8*rxt(k,447)*y(k,124) + rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1661) = mat(k,1661) + rxt(k,454)*y(k,177)
         mat(k,916) = .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,200)
         mat(k,1061) = -(rxt(k,461)*y(k,219))
         mat(k,1666) = -rxt(k,461)*y(k,180)
         mat(k,1003) = .300_r8*rxt(k,440)*y(k,135)
         mat(k,954) = .300_r8*rxt(k,443)*y(k,135)
         mat(k,1895) = .900_r8*rxt(k,452)*y(k,223)
         mat(k,2226) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,110)
         mat(k,1805) = rxt(k,450)*y(k,223)
         mat(k,1109) = .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,200)
         mat(k,657) = -(rxt(k,428)*y(k,219))
         mat(k,1632) = -rxt(k,428)*y(k,181)
         mat(k,2059) = rxt(k,426)*y(k,225)
         mat(k,724) = rxt(k,426)*y(k,205)
         mat(k,189) = -(rxt(k,429)*y(k,219))
         mat(k,1565) = -rxt(k,429)*y(k,182)
         mat(k,205) = -(rxt(k,395)*y(k,219))
         mat(k,1568) = -rxt(k,395)*y(k,183)
         mat(k,2029) = rxt(k,392)*y(k,227)
         mat(k,1168) = rxt(k,392)*y(k,205)
         mat(k,210) = -(rxt(k,431)*y(k,219))
         mat(k,1569) = -rxt(k,431)*y(k,184)
         mat(k,695) = -(rxt(k,434)*y(k,219))
         mat(k,1636) = -rxt(k,434)*y(k,185)
         mat(k,2063) = rxt(k,432)*y(k,228)
         mat(k,748) = rxt(k,432)*y(k,205)
         mat(k,218) = -(rxt(k,437)*y(k,219))
         mat(k,1570) = -rxt(k,437)*y(k,186)
         mat(k,211) = .150_r8*rxt(k,431)*y(k,219)
         mat(k,1570) = mat(k,1570) + .150_r8*rxt(k,431)*y(k,184)
         mat(k,429) = -(rxt(k,438)*y(k,219))
         mat(k,1603) = -rxt(k,438)*y(k,187)
         mat(k,2043) = rxt(k,435)*y(k,230)
         mat(k,500) = rxt(k,435)*y(k,205)
         mat(k,520) = -(rxt(k,396)*y(k,205) + rxt(k,397)*y(k,124) + rxt(k,425) &
                      *y(k,125))
         mat(k,2053) = -rxt(k,396)*y(k,190)
         mat(k,1868) = -rxt(k,397)*y(k,190)
         mat(k,1939) = -rxt(k,425)*y(k,190)
         mat(k,238) = rxt(k,402)*y(k,219)
         mat(k,1616) = rxt(k,402)*y(k,22)
         mat(k,889) = -(rxt(k,357)*y(k,205) + (rxt(k,358) + rxt(k,359)) * y(k,124))
         mat(k,2079) = -rxt(k,357)*y(k,191)
         mat(k,1888) = -(rxt(k,358) + rxt(k,359)) * y(k,191)
         mat(k,671) = rxt(k,360)*y(k,219)
         mat(k,229) = rxt(k,361)*y(k,219)
         mat(k,1655) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)
         mat(k,56) = -(rxt(k,495)*y(k,205) + rxt(k,496)*y(k,124))
         mat(k,2021) = -rxt(k,495)*y(k,192)
         mat(k,1846) = -rxt(k,496)*y(k,192)
         mat(k,991) = rxt(k,498)*y(k,219)
         mat(k,1548) = rxt(k,498)*y(k,6)
         mat(k,485) = -(rxt(k,399)*y(k,205) + rxt(k,400)*y(k,124))
         mat(k,2051) = -rxt(k,399)*y(k,193)
         mat(k,1864) = -rxt(k,400)*y(k,193)
         mat(k,167) = .350_r8*rxt(k,398)*y(k,219)
         mat(k,419) = rxt(k,401)*y(k,219)
         mat(k,1611) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)
         mat(k,62) = -(rxt(k,500)*y(k,205) + rxt(k,501)*y(k,124))
         mat(k,2022) = -rxt(k,500)*y(k,194)
         mat(k,1847) = -rxt(k,501)*y(k,194)
         mat(k,163) = rxt(k,499)*y(k,219)
         mat(k,1549) = rxt(k,499)*y(k,7)
         mat(k,437) = -(rxt(k,403)*y(k,205) + rxt(k,405)*y(k,124))
         mat(k,2044) = -rxt(k,403)*y(k,195)
         mat(k,1859) = -rxt(k,405)*y(k,195)
         mat(k,339) = rxt(k,404)*y(k,219)
         mat(k,192) = .070_r8*rxt(k,429)*y(k,219)
         mat(k,212) = .060_r8*rxt(k,431)*y(k,219)
         mat(k,1604) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,182) &
                      + .060_r8*rxt(k,431)*y(k,184)
         mat(k,825) = -(4._r8*rxt(k,280)*y(k,196) + rxt(k,281)*y(k,200) + rxt(k,282) &
                      *y(k,205) + rxt(k,283)*y(k,124))
         mat(k,1798) = -rxt(k,281)*y(k,196)
         mat(k,2076) = -rxt(k,282)*y(k,196)
         mat(k,1885) = -rxt(k,283)*y(k,196)
         mat(k,344) = .500_r8*rxt(k,285)*y(k,219)
         mat(k,307) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,219)
         mat(k,1988) = rxt(k,286)*y(k,28)
         mat(k,1649) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)
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
         mat(k,789) = -(rxt(k,309)*y(k,200) + rxt(k,310)*y(k,205) + rxt(k,311) &
                      *y(k,124))
         mat(k,1795) = -rxt(k,309)*y(k,197)
         mat(k,2072) = -rxt(k,310)*y(k,197)
         mat(k,1883) = -rxt(k,311)*y(k,197)
         mat(k,412) = rxt(k,312)*y(k,219)
         mat(k,112) = rxt(k,313)*y(k,219)
         mat(k,1644) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)
         mat(k,627) = -(rxt(k,406)*y(k,205) + rxt(k,407)*y(k,124))
         mat(k,2057) = -rxt(k,406)*y(k,198)
         mat(k,1874) = -rxt(k,407)*y(k,198)
         mat(k,269) = rxt(k,408)*y(k,219)
         mat(k,1874) = mat(k,1874) + rxt(k,397)*y(k,190)
         mat(k,2211) = rxt(k,423)*y(k,142)
         mat(k,466) = rxt(k,423)*y(k,135)
         mat(k,521) = rxt(k,397)*y(k,124) + .400_r8*rxt(k,396)*y(k,205)
         mat(k,2057) = mat(k,2057) + .400_r8*rxt(k,396)*y(k,190)
         mat(k,1629) = rxt(k,408)*y(k,32)
         mat(k,1388) = -(4._r8*rxt(k,291)*y(k,199) + rxt(k,292)*y(k,200) + rxt(k,293) &
                      *y(k,205) + rxt(k,294)*y(k,124) + rxt(k,305)*y(k,125) + rxt(k,332) &
                      *y(k,211) + rxt(k,365)*y(k,207) + rxt(k,370)*y(k,208) + rxt(k,379) &
                      *y(k,101) + rxt(k,390)*y(k,227))
         mat(k,1822) = -rxt(k,292)*y(k,199)
         mat(k,2102) = -rxt(k,293)*y(k,199)
         mat(k,1913) = -rxt(k,294)*y(k,199)
         mat(k,1957) = -rxt(k,305)*y(k,199)
         mat(k,1315) = -rxt(k,332)*y(k,199)
         mat(k,1262) = -rxt(k,365)*y(k,199)
         mat(k,1294) = -rxt(k,370)*y(k,199)
         mat(k,1199) = -rxt(k,379)*y(k,199)
         mat(k,1177) = -rxt(k,390)*y(k,199)
         mat(k,1009) = .060_r8*rxt(k,440)*y(k,135)
         mat(k,1089) = rxt(k,288)*y(k,126) + rxt(k,289)*y(k,219)
         mat(k,1224) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219)
         mat(k,614) = .500_r8*rxt(k,296)*y(k,219)
         mat(k,845) = .080_r8*rxt(k,385)*y(k,135)
         mat(k,1215) = .100_r8*rxt(k,338)*y(k,135)
         mat(k,959) = .060_r8*rxt(k,443)*y(k,135)
         mat(k,1336) = .280_r8*rxt(k,352)*y(k,135)
         mat(k,1913) = mat(k,1913) + .530_r8*rxt(k,336)*y(k,211) + rxt(k,345)*y(k,213) &
                      + rxt(k,348)*y(k,215) + rxt(k,323)*y(k,222)
         mat(k,1769) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,211) + rxt(k,346)*y(k,213)
         mat(k,2242) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,98) &
                      + .100_r8*rxt(k,338)*y(k,105) + .060_r8*rxt(k,443)*y(k,110) &
                      + .280_r8*rxt(k,352)*y(k,111)
         mat(k,1064) = .650_r8*rxt(k,461)*y(k,219)
         mat(k,1388) = mat(k,1388) + .530_r8*rxt(k,332)*y(k,211)
         mat(k,1822) = mat(k,1822) + .260_r8*rxt(k,333)*y(k,211) + rxt(k,342)*y(k,213) &
                      + .300_r8*rxt(k,321)*y(k,222)
         mat(k,2102) = mat(k,2102) + .450_r8*rxt(k,343)*y(k,213) + .200_r8*rxt(k,347) &
                      *y(k,215) + .150_r8*rxt(k,322)*y(k,222)
         mat(k,1315) = mat(k,1315) + .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335) &
                      *y(k,126) + .530_r8*rxt(k,332)*y(k,199) + .260_r8*rxt(k,333) &
                      *y(k,200)
         mat(k,1357) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,342)*y(k,200) &
                      + .450_r8*rxt(k,343)*y(k,205) + 4.000_r8*rxt(k,344)*y(k,213)
         mat(k,681) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,205)
         mat(k,1685) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,180)
         mat(k,1134) = rxt(k,323)*y(k,124) + .300_r8*rxt(k,321)*y(k,200) &
                      + .150_r8*rxt(k,322)*y(k,205)
         mat(k,1831) = -(rxt(k,182)*y(k,59) + (4._r8*rxt(k,259) + 4._r8*rxt(k,260) &
                      ) * y(k,200) + rxt(k,261)*y(k,205) + rxt(k,262)*y(k,124) &
                      + rxt(k,281)*y(k,196) + rxt(k,292)*y(k,199) + rxt(k,309) &
                      *y(k,197) + rxt(k,321)*y(k,222) + rxt(k,333)*y(k,211) + rxt(k,342) &
                      *y(k,213) + rxt(k,366)*y(k,207) + rxt(k,371)*y(k,208) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,227) + rxt(k,445)*y(k,217) + rxt(k,450) &
                      *y(k,223) + rxt(k,455)*y(k,224))
         mat(k,1722) = -rxt(k,182)*y(k,200)
         mat(k,2113) = -rxt(k,261)*y(k,200)
         mat(k,1923) = -rxt(k,262)*y(k,200)
         mat(k,829) = -rxt(k,281)*y(k,200)
         mat(k,1394) = -rxt(k,292)*y(k,200)
         mat(k,794) = -rxt(k,309)*y(k,200)
         mat(k,1137) = -rxt(k,321)*y(k,200)
         mat(k,1320) = -rxt(k,333)*y(k,200)
         mat(k,1362) = -rxt(k,342)*y(k,200)
         mat(k,1267) = -rxt(k,366)*y(k,200)
         mat(k,1299) = -rxt(k,371)*y(k,200)
         mat(k,1204) = -rxt(k,380)*y(k,200)
         mat(k,1181) = -rxt(k,391)*y(k,200)
         mat(k,1055) = -rxt(k,445)*y(k,200)
         mat(k,1118) = -rxt(k,450)*y(k,200)
         mat(k,920) = -rxt(k,455)*y(k,200)
         mat(k,1036) = .280_r8*rxt(k,308)*y(k,135)
         mat(k,689) = rxt(k,295)*y(k,219)
         mat(k,460) = .700_r8*rxt(k,264)*y(k,219)
         mat(k,1439) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,271)*y(k,218) &
                      + rxt(k,265)*y(k,219)
         mat(k,2006) = rxt(k,176)*y(k,54)
         mat(k,873) = rxt(k,232)*y(k,54)
         mat(k,849) = .050_r8*rxt(k,385)*y(k,135)
         mat(k,1204) = mat(k,1204) + rxt(k,379)*y(k,199)
         mat(k,1923) = mat(k,1923) + rxt(k,294)*y(k,199) + .830_r8*rxt(k,411)*y(k,201) &
                      + .170_r8*rxt(k,417)*y(k,214)
         mat(k,2252) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,98)
         mat(k,1394) = mat(k,1394) + rxt(k,379)*y(k,101) + rxt(k,294)*y(k,124) &
                      + 4.000_r8*rxt(k,291)*y(k,199) + .900_r8*rxt(k,292)*y(k,200) &
                      + .450_r8*rxt(k,293)*y(k,205) + rxt(k,365)*y(k,207) + rxt(k,370) &
                      *y(k,208) + rxt(k,332)*y(k,211) + rxt(k,341)*y(k,213) &
                      + rxt(k,390)*y(k,227)
         mat(k,1831) = mat(k,1831) + .900_r8*rxt(k,292)*y(k,199)
         mat(k,765) = .830_r8*rxt(k,411)*y(k,124) + .330_r8*rxt(k,410)*y(k,205)
         mat(k,2113) = mat(k,2113) + .450_r8*rxt(k,293)*y(k,199) + .330_r8*rxt(k,410) &
                      *y(k,201) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1267) = mat(k,1267) + rxt(k,365)*y(k,199)
         mat(k,1299) = mat(k,1299) + rxt(k,370)*y(k,199)
         mat(k,1320) = mat(k,1320) + rxt(k,332)*y(k,199)
         mat(k,1362) = mat(k,1362) + rxt(k,341)*y(k,199)
         mat(k,880) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1532) = rxt(k,271)*y(k,54)
         mat(k,1696) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,264)*y(k,53) + rxt(k,265) &
                      *y(k,54)
         mat(k,1181) = mat(k,1181) + rxt(k,390)*y(k,199)
         mat(k,761) = -(rxt(k,410)*y(k,205) + rxt(k,411)*y(k,124) + rxt(k,412) &
                      *y(k,125))
         mat(k,2069) = -rxt(k,410)*y(k,201)
         mat(k,1881) = -rxt(k,411)*y(k,201)
         mat(k,1945) = -rxt(k,412)*y(k,201)
         mat(k,566) = -((rxt(k,329) + rxt(k,330)) * y(k,124))
         mat(k,1870) = -(rxt(k,329) + rxt(k,330)) * y(k,202)
         mat(k,352) = rxt(k,328)*y(k,219)
         mat(k,1621) = rxt(k,328)*y(k,16)
         mat(k,1855) = .750_r8*rxt(k,298)*y(k,204)
         mat(k,707) = .750_r8*rxt(k,298)*y(k,124)
         mat(k,708) = -(rxt(k,297)*y(k,205) + rxt(k,298)*y(k,124))
         mat(k,2064) = -rxt(k,297)*y(k,204)
         mat(k,1877) = -rxt(k,298)*y(k,204)
         mat(k,551) = rxt(k,304)*y(k,219)
         mat(k,1637) = rxt(k,304)*y(k,25)
         mat(k,2117) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,76) + rxt(k,139) &
                      *y(k,134) + rxt(k,140)*y(k,135) + rxt(k,144)*y(k,219) &
                      + 4._r8*rxt(k,149)*y(k,205) + rxt(k,159)*y(k,126) + rxt(k,164) &
                      *y(k,124) + rxt(k,169)*y(k,125) + (rxt(k,179) + rxt(k,180) &
                      ) * y(k,56) + rxt(k,186)*y(k,59) + rxt(k,212)*y(k,17) + rxt(k,218) &
                      *y(k,19) + rxt(k,255)*y(k,42) + rxt(k,261)*y(k,200) + rxt(k,268) &
                      *y(k,206) + rxt(k,282)*y(k,196) + rxt(k,293)*y(k,199) + rxt(k,297) &
                      *y(k,204) + rxt(k,310)*y(k,197) + rxt(k,318)*y(k,221) + rxt(k,322) &
                      *y(k,222) + rxt(k,334)*y(k,211) + rxt(k,343)*y(k,213) + rxt(k,347) &
                      *y(k,215) + rxt(k,357)*y(k,191) + rxt(k,367)*y(k,207) + rxt(k,372) &
                      *y(k,208) + rxt(k,381)*y(k,101) + rxt(k,392)*y(k,227) + rxt(k,396) &
                      *y(k,190) + rxt(k,399)*y(k,193) + rxt(k,403)*y(k,195) + rxt(k,406) &
                      *y(k,198) + rxt(k,410)*y(k,201) + rxt(k,413)*y(k,212) + rxt(k,416) &
                      *y(k,214) + rxt(k,419)*y(k,220) + rxt(k,426)*y(k,225) + rxt(k,432) &
                      *y(k,228) + rxt(k,435)*y(k,230) + rxt(k,446)*y(k,217) + rxt(k,451) &
                      *y(k,223) + rxt(k,456)*y(k,224))
         mat(k,1455) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,205)
         mat(k,2195) = -rxt(k,139)*y(k,205)
         mat(k,2256) = -rxt(k,140)*y(k,205)
         mat(k,1700) = -rxt(k,144)*y(k,205)
         mat(k,1783) = -rxt(k,159)*y(k,205)
         mat(k,1927) = -rxt(k,164)*y(k,205)
         mat(k,1971) = -rxt(k,169)*y(k,205)
         mat(k,2010) = -(rxt(k,179) + rxt(k,180)) * y(k,205)
         mat(k,1726) = -rxt(k,186)*y(k,205)
         mat(k,1422) = -rxt(k,212)*y(k,205)
         mat(k,2164) = -rxt(k,218)*y(k,205)
         mat(k,1495) = -rxt(k,255)*y(k,205)
         mat(k,1835) = -rxt(k,261)*y(k,205)
         mat(k,448) = -rxt(k,268)*y(k,205)
         mat(k,832) = -rxt(k,282)*y(k,205)
         mat(k,1397) = -rxt(k,293)*y(k,205)
         mat(k,714) = -rxt(k,297)*y(k,205)
         mat(k,797) = -rxt(k,310)*y(k,205)
         mat(k,777) = -rxt(k,318)*y(k,205)
         mat(k,1140) = -rxt(k,322)*y(k,205)
         mat(k,1323) = -rxt(k,334)*y(k,205)
         mat(k,1365) = -rxt(k,343)*y(k,205)
         mat(k,685) = -rxt(k,347)*y(k,205)
         mat(k,898) = -rxt(k,357)*y(k,205)
         mat(k,1270) = -rxt(k,367)*y(k,205)
         mat(k,1302) = -rxt(k,372)*y(k,205)
         mat(k,1207) = -rxt(k,381)*y(k,205)
         mat(k,1184) = -rxt(k,392)*y(k,205)
         mat(k,525) = -rxt(k,396)*y(k,205)
         mat(k,491) = -rxt(k,399)*y(k,205)
         mat(k,442) = -rxt(k,403)*y(k,205)
         mat(k,631) = -rxt(k,406)*y(k,205)
         mat(k,768) = -rxt(k,410)*y(k,205)
         mat(k,720) = -rxt(k,413)*y(k,205)
         mat(k,883) = -rxt(k,416)*y(k,205)
         mat(k,455) = -rxt(k,419)*y(k,205)
         mat(k,735) = -rxt(k,426)*y(k,205)
         mat(k,760) = -rxt(k,432)*y(k,205)
         mat(k,507) = -rxt(k,435)*y(k,205)
         mat(k,1058) = -rxt(k,446)*y(k,205)
         mat(k,1121) = -rxt(k,451)*y(k,205)
         mat(k,923) = -rxt(k,456)*y(k,205)
         mat(k,1017) = .570_r8*rxt(k,440)*y(k,135)
         mat(k,169) = .650_r8*rxt(k,398)*y(k,219)
         mat(k,1422) = mat(k,1422) + rxt(k,211)*y(k,42)
         mat(k,2164) = mat(k,2164) + rxt(k,223)*y(k,219)
         mat(k,298) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,556) = .130_r8*rxt(k,279)*y(k,135)
         mat(k,266) = rxt(k,284)*y(k,219)
         mat(k,1039) = .280_r8*rxt(k,308)*y(k,135)
         mat(k,1495) = mat(k,1495) + rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56) &
                      + rxt(k,256)*y(k,126) + rxt(k,257)*y(k,134)
         mat(k,602) = rxt(k,240)*y(k,56) + rxt(k,241)*y(k,219)
         mat(k,372) = rxt(k,243)*y(k,56) + rxt(k,244)*y(k,219)
         mat(k,106) = rxt(k,290)*y(k,219)
         mat(k,802) = rxt(k,263)*y(k,219)
         mat(k,1441) = rxt(k,272)*y(k,218)
         mat(k,2010) = mat(k,2010) + rxt(k,175)*y(k,42) + rxt(k,240)*y(k,43) &
                      + rxt(k,243)*y(k,46) + rxt(k,178)*y(k,79)
         mat(k,1726) = mat(k,1726) + rxt(k,182)*y(k,200) + rxt(k,193)*y(k,219)
         mat(k,1127) = rxt(k,275)*y(k,219)
         mat(k,200) = .730_r8*rxt(k,409)*y(k,219)
         mat(k,289) = .500_r8*rxt(k,477)*y(k,219)
         mat(k,1102) = rxt(k,301)*y(k,219)
         mat(k,984) = rxt(k,302)*y(k,219)
         mat(k,608) = rxt(k,178)*y(k,56) + rxt(k,134)*y(k,134) + rxt(k,143)*y(k,219)
         mat(k,184) = rxt(k,266)*y(k,219)
         mat(k,934) = rxt(k,267)*y(k,219)
         mat(k,1165) = rxt(k,331)*y(k,219)
         mat(k,1147) = rxt(k,316)*y(k,219)
         mat(k,850) = .370_r8*rxt(k,385)*y(k,135)
         mat(k,591) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,541) = rxt(k,377)*y(k,219)
         mat(k,1207) = mat(k,1207) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) &
                      + rxt(k,379)*y(k,199) + 1.200_r8*rxt(k,380)*y(k,200)
         mat(k,410) = rxt(k,384)*y(k,219)
         mat(k,1218) = .140_r8*rxt(k,338)*y(k,135)
         mat(k,316) = .200_r8*rxt(k,340)*y(k,219)
         mat(k,582) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,967) = .570_r8*rxt(k,443)*y(k,135)
         mat(k,1345) = .280_r8*rxt(k,352)*y(k,135)
         mat(k,380) = rxt(k,388)*y(k,219)
         mat(k,1085) = rxt(k,389)*y(k,219)
         mat(k,1927) = mat(k,1927) + rxt(k,382)*y(k,101) + rxt(k,358)*y(k,191) &
                      + rxt(k,400)*y(k,193) + rxt(k,405)*y(k,195) + rxt(k,283) &
                      *y(k,196) + rxt(k,311)*y(k,197) + rxt(k,262)*y(k,200) &
                      + .170_r8*rxt(k,411)*y(k,201) + rxt(k,329)*y(k,202) &
                      + .250_r8*rxt(k,298)*y(k,204) + rxt(k,270)*y(k,206) &
                      + .920_r8*rxt(k,368)*y(k,207) + .920_r8*rxt(k,374)*y(k,208) &
                      + .470_r8*rxt(k,336)*y(k,211) + .400_r8*rxt(k,414)*y(k,212) &
                      + .830_r8*rxt(k,417)*y(k,214) + rxt(k,420)*y(k,220) + rxt(k,319) &
                      *y(k,221) + .900_r8*rxt(k,452)*y(k,223) + .800_r8*rxt(k,457) &
                      *y(k,224) + rxt(k,427)*y(k,225) + rxt(k,393)*y(k,227) &
                      + rxt(k,433)*y(k,228) + rxt(k,436)*y(k,230)
         mat(k,1783) = mat(k,1783) + rxt(k,256)*y(k,42) + rxt(k,383)*y(k,101) &
                      + rxt(k,369)*y(k,207) + rxt(k,375)*y(k,208) + .470_r8*rxt(k,335) &
                      *y(k,211) + rxt(k,162)*y(k,219) + rxt(k,394)*y(k,227)
         mat(k,2195) = mat(k,2195) + rxt(k,257)*y(k,42) + rxt(k,134)*y(k,79)
         mat(k,2256) = mat(k,2256) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,98) + .140_r8*rxt(k,338)*y(k,105) + .570_r8*rxt(k,443) &
                      *y(k,110) + .280_r8*rxt(k,352)*y(k,111) + rxt(k,146)*y(k,219)
         mat(k,178) = .800_r8*rxt(k,421)*y(k,219)
         mat(k,856) = rxt(k,467)*y(k,219)
         mat(k,1068) = .200_r8*rxt(k,461)*y(k,219)
         mat(k,195) = .280_r8*rxt(k,429)*y(k,219)
         mat(k,217) = .380_r8*rxt(k,431)*y(k,219)
         mat(k,222) = .630_r8*rxt(k,437)*y(k,219)
         mat(k,898) = mat(k,898) + rxt(k,358)*y(k,124)
         mat(k,491) = mat(k,491) + rxt(k,400)*y(k,124)
         mat(k,442) = mat(k,442) + rxt(k,405)*y(k,124)
         mat(k,832) = mat(k,832) + rxt(k,283)*y(k,124) + 2.400_r8*rxt(k,280)*y(k,196) &
                      + rxt(k,281)*y(k,200)
         mat(k,797) = mat(k,797) + rxt(k,311)*y(k,124) + rxt(k,309)*y(k,200)
         mat(k,1397) = mat(k,1397) + rxt(k,379)*y(k,101) + .900_r8*rxt(k,292)*y(k,200) &
                      + rxt(k,365)*y(k,207) + rxt(k,370)*y(k,208) + .470_r8*rxt(k,332) &
                      *y(k,211) + rxt(k,390)*y(k,227)
         mat(k,1835) = mat(k,1835) + rxt(k,182)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,101) &
                      + rxt(k,262)*y(k,124) + rxt(k,281)*y(k,196) + rxt(k,309) &
                      *y(k,197) + .900_r8*rxt(k,292)*y(k,199) + 4.000_r8*rxt(k,259) &
                      *y(k,200) + rxt(k,366)*y(k,207) + rxt(k,371)*y(k,208) &
                      + .730_r8*rxt(k,333)*y(k,211) + rxt(k,342)*y(k,213) &
                      + .500_r8*rxt(k,445)*y(k,217) + .300_r8*rxt(k,321)*y(k,222) &
                      + rxt(k,450)*y(k,223) + rxt(k,455)*y(k,224) + .800_r8*rxt(k,391) &
                      *y(k,227)
         mat(k,768) = mat(k,768) + .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410) &
                      *y(k,205)
         mat(k,573) = rxt(k,329)*y(k,124)
         mat(k,714) = mat(k,714) + .250_r8*rxt(k,298)*y(k,124)
         mat(k,2117) = mat(k,2117) + .070_r8*rxt(k,410)*y(k,201) + .160_r8*rxt(k,413) &
                      *y(k,212) + .330_r8*rxt(k,416)*y(k,214)
         mat(k,448) = mat(k,448) + rxt(k,270)*y(k,124)
         mat(k,1270) = mat(k,1270) + .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) &
                      + rxt(k,365)*y(k,199) + rxt(k,366)*y(k,200)
         mat(k,1302) = mat(k,1302) + .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,199) + rxt(k,371)*y(k,200)
         mat(k,1323) = mat(k,1323) + .470_r8*rxt(k,336)*y(k,124) + .470_r8*rxt(k,335) &
                      *y(k,126) + .470_r8*rxt(k,332)*y(k,199) + .730_r8*rxt(k,333) &
                      *y(k,200)
         mat(k,720) = mat(k,720) + .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413) &
                      *y(k,205)
         mat(k,1365) = mat(k,1365) + rxt(k,342)*y(k,200)
         mat(k,883) = mat(k,883) + .830_r8*rxt(k,417)*y(k,124) + .330_r8*rxt(k,416) &
                      *y(k,205)
         mat(k,1058) = mat(k,1058) + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1536) = rxt(k,272)*y(k,54)
         mat(k,1700) = mat(k,1700) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,223)*y(k,19) &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,241) &
                      *y(k,43) + rxt(k,244)*y(k,46) + rxt(k,290)*y(k,47) + rxt(k,263) &
                      *y(k,52) + rxt(k,193)*y(k,59) + rxt(k,275)*y(k,62) &
                      + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,477)*y(k,67) &
                      + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,143)*y(k,79) &
                      + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,331)*y(k,93) &
                      + rxt(k,316)*y(k,95) + .300_r8*rxt(k,376)*y(k,99) + rxt(k,377) &
                      *y(k,100) + rxt(k,384)*y(k,102) + .200_r8*rxt(k,340)*y(k,106) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,388)*y(k,115) + rxt(k,389) &
                      *y(k,116) + rxt(k,162)*y(k,126) + rxt(k,146)*y(k,135) &
                      + .800_r8*rxt(k,421)*y(k,143) + rxt(k,467)*y(k,152) &
                      + .200_r8*rxt(k,461)*y(k,180) + .280_r8*rxt(k,429)*y(k,182) &
                      + .380_r8*rxt(k,431)*y(k,184) + .630_r8*rxt(k,437)*y(k,186)
         mat(k,455) = mat(k,455) + rxt(k,420)*y(k,124)
         mat(k,777) = mat(k,777) + rxt(k,319)*y(k,124)
         mat(k,1140) = mat(k,1140) + .300_r8*rxt(k,321)*y(k,200)
         mat(k,1121) = mat(k,1121) + .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,200)
         mat(k,923) = mat(k,923) + .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,200)
         mat(k,735) = mat(k,735) + rxt(k,427)*y(k,124)
         mat(k,1184) = mat(k,1184) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126) &
                      + rxt(k,390)*y(k,199) + .800_r8*rxt(k,391)*y(k,200)
         mat(k,760) = mat(k,760) + rxt(k,433)*y(k,124)
         mat(k,507) = mat(k,507) + rxt(k,436)*y(k,124)
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
         mat(k,443) = -(rxt(k,268)*y(k,205) + rxt(k,270)*y(k,124))
         mat(k,2045) = -rxt(k,268)*y(k,206)
         mat(k,1860) = -rxt(k,270)*y(k,206)
         mat(k,1478) = rxt(k,255)*y(k,205)
         mat(k,2045) = mat(k,2045) + rxt(k,255)*y(k,42)
         mat(k,1258) = -(rxt(k,365)*y(k,199) + rxt(k,366)*y(k,200) + rxt(k,367) &
                      *y(k,205) + rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126))
         mat(k,1383) = -rxt(k,365)*y(k,207)
         mat(k,1817) = -rxt(k,366)*y(k,207)
         mat(k,2097) = -rxt(k,367)*y(k,207)
         mat(k,1908) = -rxt(k,368)*y(k,207)
         mat(k,1764) = -rxt(k,369)*y(k,207)
         mat(k,842) = .600_r8*rxt(k,386)*y(k,219)
         mat(k,1680) = .600_r8*rxt(k,386)*y(k,98)
         mat(k,1290) = -(rxt(k,370)*y(k,199) + rxt(k,371)*y(k,200) + rxt(k,372) &
                      *y(k,205) + rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126))
         mat(k,1384) = -rxt(k,370)*y(k,208)
         mat(k,1818) = -rxt(k,371)*y(k,208)
         mat(k,2098) = -rxt(k,372)*y(k,208)
         mat(k,1909) = -rxt(k,374)*y(k,208)
         mat(k,1765) = -rxt(k,375)*y(k,208)
         mat(k,843) = .400_r8*rxt(k,386)*y(k,219)
         mat(k,1681) = .400_r8*rxt(k,386)*y(k,98)
         mat(k,68) = -(rxt(k,503)*y(k,205) + rxt(k,504)*y(k,124))
         mat(k,2023) = -rxt(k,503)*y(k,209)
         mat(k,1848) = -rxt(k,504)*y(k,209)
         mat(k,835) = rxt(k,506)*y(k,219)
         mat(k,1550) = rxt(k,506)*y(k,98)
         mat(k,74) = -(rxt(k,507)*y(k,205) + rxt(k,508)*y(k,124))
         mat(k,2024) = -rxt(k,507)*y(k,210)
         mat(k,1849) = -rxt(k,508)*y(k,210)
         mat(k,75) = rxt(k,509)*y(k,219)
         mat(k,1551) = rxt(k,509)*y(k,104)
         mat(k,1313) = -(rxt(k,332)*y(k,199) + rxt(k,333)*y(k,200) + rxt(k,334) &
                      *y(k,205) + rxt(k,335)*y(k,126) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,124))
         mat(k,1385) = -rxt(k,332)*y(k,211)
         mat(k,1819) = -rxt(k,333)*y(k,211)
         mat(k,2099) = -rxt(k,334)*y(k,211)
         mat(k,1766) = -rxt(k,335)*y(k,211)
         mat(k,1910) = -(rxt(k,336) + rxt(k,337)) * y(k,211)
         mat(k,1213) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,313) = .200_r8*rxt(k,340)*y(k,219)
         mat(k,1333) = rxt(k,353)*y(k,219)
         mat(k,1682) = .500_r8*rxt(k,339)*y(k,105) + .200_r8*rxt(k,340)*y(k,106) &
                      + rxt(k,353)*y(k,111)
         mat(k,715) = -(rxt(k,413)*y(k,205) + rxt(k,414)*y(k,124) + rxt(k,415) &
                      *y(k,125))
         mat(k,2065) = -rxt(k,413)*y(k,212)
         mat(k,1878) = -rxt(k,414)*y(k,212)
         mat(k,1944) = -rxt(k,415)*y(k,212)
         mat(k,1356) = -(rxt(k,341)*y(k,199) + rxt(k,342)*y(k,200) + rxt(k,343) &
                      *y(k,205) + 4._r8*rxt(k,344)*y(k,213) + rxt(k,345)*y(k,124) &
                      + rxt(k,346)*y(k,126) + rxt(k,354)*y(k,125))
         mat(k,1387) = -rxt(k,341)*y(k,213)
         mat(k,1821) = -rxt(k,342)*y(k,213)
         mat(k,2101) = -rxt(k,343)*y(k,213)
         mat(k,1912) = -rxt(k,345)*y(k,213)
         mat(k,1768) = -rxt(k,346)*y(k,213)
         mat(k,1956) = -rxt(k,354)*y(k,213)
         mat(k,1214) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,314) = .500_r8*rxt(k,340)*y(k,219)
         mat(k,1684) = .500_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,340)*y(k,106)
         mat(k,875) = -(rxt(k,416)*y(k,205) + rxt(k,417)*y(k,124) + rxt(k,418) &
                      *y(k,125))
         mat(k,2078) = -rxt(k,416)*y(k,214)
         mat(k,1887) = -rxt(k,417)*y(k,214)
         mat(k,1949) = -rxt(k,418)*y(k,214)
         mat(k,679) = -(rxt(k,347)*y(k,205) + rxt(k,348)*y(k,124))
         mat(k,2061) = -rxt(k,347)*y(k,215)
         mat(k,1876) = -rxt(k,348)*y(k,215)
         mat(k,509) = rxt(k,349)*y(k,219)
         mat(k,318) = rxt(k,350)*y(k,219)
         mat(k,1634) = rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108)
         mat(k,82) = -(rxt(k,511)*y(k,205) + rxt(k,512)*y(k,124))
         mat(k,2025) = -rxt(k,511)*y(k,216)
         mat(k,1850) = -rxt(k,512)*y(k,216)
         mat(k,942) = rxt(k,514)*y(k,219)
         mat(k,1553) = rxt(k,514)*y(k,110)
         mat(k,1047) = -(rxt(k,445)*y(k,200) + rxt(k,446)*y(k,205) + rxt(k,447) &
                      *y(k,124) + rxt(k,448)*y(k,126))
         mat(k,1804) = -rxt(k,445)*y(k,217)
         mat(k,2085) = -rxt(k,446)*y(k,217)
         mat(k,1894) = -rxt(k,447)*y(k,217)
         mat(k,1749) = -rxt(k,448)*y(k,217)
         mat(k,1002) = rxt(k,439)*y(k,126)
         mat(k,953) = rxt(k,442)*y(k,126)
         mat(k,1749) = mat(k,1749) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,110) &
                      + .500_r8*rxt(k,459)*y(k,179)
         mat(k,395) = rxt(k,449)*y(k,219)
         mat(k,972) = .500_r8*rxt(k,459)*y(k,126)
         mat(k,1665) = rxt(k,449)*y(k,128)
         mat(k,1528) = -(rxt(k,125)*y(k,77) + rxt(k,126)*y(k,231) + rxt(k,129) &
                      *y(k,135) + (rxt(k,167) + rxt(k,168)) * y(k,113) + rxt(k,200) &
                      *y(k,33) + rxt(k,201)*y(k,34) + rxt(k,202)*y(k,36) + rxt(k,203) &
                      *y(k,37) + rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39) + rxt(k,206) &
                      *y(k,40) + (rxt(k,207) + rxt(k,208)) * y(k,85) + rxt(k,227) &
                      *y(k,35) + rxt(k,228)*y(k,55) + rxt(k,229)*y(k,78) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,81) + rxt(k,236)*y(k,64) + rxt(k,237) &
                      *y(k,65) + rxt(k,250)*y(k,41) + rxt(k,251)*y(k,43) + rxt(k,252) &
                      *y(k,82) + rxt(k,253)*y(k,83) + rxt(k,254)*y(k,84) + (rxt(k,271) &
                      + rxt(k,272) + rxt(k,273)) * y(k,54) + rxt(k,274)*y(k,86))
         mat(k,1407) = -rxt(k,125)*y(k,218)
         mat(k,2274) = -rxt(k,126)*y(k,218)
         mat(k,2248) = -rxt(k,129)*y(k,218)
         mat(k,186) = -(rxt(k,167) + rxt(k,168)) * y(k,218)
         mat(k,102) = -rxt(k,200)*y(k,218)
         mat(k,146) = -rxt(k,201)*y(k,218)
         mat(k,117) = -rxt(k,202)*y(k,218)
         mat(k,156) = -rxt(k,203)*y(k,218)
         mat(k,121) = -rxt(k,204)*y(k,218)
         mat(k,161) = -rxt(k,205)*y(k,218)
         mat(k,125) = -rxt(k,206)*y(k,218)
         mat(k,2132) = -(rxt(k,207) + rxt(k,208)) * y(k,218)
         mat(k,152) = -rxt(k,227)*y(k,218)
         mat(k,389) = -rxt(k,228)*y(k,218)
         mat(k,110) = -rxt(k,229)*y(k,218)
         mat(k,809) = -(rxt(k,230) + rxt(k,231)) * y(k,218)
         mat(k,242) = -rxt(k,236)*y(k,218)
         mat(k,250) = -rxt(k,237)*y(k,218)
         mat(k,471) = -rxt(k,250)*y(k,218)
         mat(k,598) = -rxt(k,251)*y(k,218)
         mat(k,245) = -rxt(k,252)*y(k,218)
         mat(k,255) = -rxt(k,253)*y(k,218)
         mat(k,302) = -rxt(k,254)*y(k,218)
         mat(k,1436) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,218)
         mat(k,182) = -rxt(k,274)*y(k,218)
         mat(k,1693) = -(rxt(k,142)*y(k,77) + rxt(k,143)*y(k,79) + rxt(k,144)*y(k,205) &
                      + rxt(k,145)*y(k,134) + rxt(k,146)*y(k,135) + (4._r8*rxt(k,147) &
                      + 4._r8*rxt(k,148)) * y(k,219) + rxt(k,150)*y(k,90) + rxt(k,162) &
                      *y(k,126) + rxt(k,163)*y(k,112) + rxt(k,171)*y(k,125) + rxt(k,172) &
                      *y(k,89) + rxt(k,191)*y(k,60) + (rxt(k,193) + rxt(k,194) &
                      ) * y(k,59) + rxt(k,196)*y(k,85) + rxt(k,199)*y(k,92) + rxt(k,223) &
                      *y(k,19) + rxt(k,225)*y(k,81) + rxt(k,239)*y(k,41) + rxt(k,241) &
                      *y(k,43) + rxt(k,242)*y(k,44) + rxt(k,244)*y(k,46) + rxt(k,246) &
                      *y(k,55) + rxt(k,247)*y(k,82) + rxt(k,248)*y(k,83) + rxt(k,249) &
                      *y(k,84) + rxt(k,258)*y(k,42) + rxt(k,263)*y(k,52) + rxt(k,264) &
                      *y(k,53) + rxt(k,265)*y(k,54) + rxt(k,266)*y(k,86) + rxt(k,267) &
                      *y(k,87) + rxt(k,275)*y(k,62) + rxt(k,277)*y(k,24) + rxt(k,284) &
                      *y(k,26) + rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28) + rxt(k,289) &
                      *y(k,45) + rxt(k,290)*y(k,47) + rxt(k,295)*y(k,50) + rxt(k,296) &
                      *y(k,51) + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,303) &
                      *y(k,140) + rxt(k,304)*y(k,25) + rxt(k,312)*y(k,30) + rxt(k,313) &
                      *y(k,31) + rxt(k,315)*y(k,49) + rxt(k,316)*y(k,95) + rxt(k,317) &
                      *y(k,127) + rxt(k,320)*y(k,147) + rxt(k,324)*y(k,148) + rxt(k,325) &
                      *y(k,29) + rxt(k,326)*y(k,48) + rxt(k,328)*y(k,16) + rxt(k,331) &
                      *y(k,93) + rxt(k,339)*y(k,105) + rxt(k,340)*y(k,106) + rxt(k,349) &
                      *y(k,107) + rxt(k,350)*y(k,108) + rxt(k,351)*y(k,109) + rxt(k,353) &
                      *y(k,111) + rxt(k,356)*y(k,1) + rxt(k,360)*y(k,2) + rxt(k,361) &
                      *y(k,15) + rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364) &
                      *y(k,97) + rxt(k,376)*y(k,99) + rxt(k,377)*y(k,100) + rxt(k,384) &
                      *y(k,102) + rxt(k,386)*y(k,98) + rxt(k,387)*y(k,103) + rxt(k,388) &
                      *y(k,115) + rxt(k,389)*y(k,116) + rxt(k,395)*y(k,183) + rxt(k,398) &
                      *y(k,7) + rxt(k,401)*y(k,8) + rxt(k,402)*y(k,22) + rxt(k,404) &
                      *y(k,23) + rxt(k,408)*y(k,32) + rxt(k,409)*y(k,66) + rxt(k,421) &
                      *y(k,143) + rxt(k,424)*y(k,144) + rxt(k,428)*y(k,181) + rxt(k,429) &
                      *y(k,182) + rxt(k,431)*y(k,184) + rxt(k,434)*y(k,185) + rxt(k,437) &
                      *y(k,186) + rxt(k,438)*y(k,187) + rxt(k,441)*y(k,6) + rxt(k,444) &
                      *y(k,110) + rxt(k,449)*y(k,128) + rxt(k,453)*y(k,176) + rxt(k,454) &
                      *y(k,177) + rxt(k,458)*y(k,178) + rxt(k,460)*y(k,179) + rxt(k,461) &
                      *y(k,180) + (rxt(k,463) + rxt(k,477)) * y(k,67) + rxt(k,465) &
                      *y(k,138) + rxt(k,467)*y(k,152) + rxt(k,471)*y(k,149) + rxt(k,476) &
                      *y(k,151) + rxt(k,479)*y(k,120))
         mat(k,1408) = -rxt(k,142)*y(k,219)
         mat(k,606) = -rxt(k,143)*y(k,219)
         mat(k,2110) = -rxt(k,144)*y(k,219)
         mat(k,2188) = -rxt(k,145)*y(k,219)
         mat(k,2249) = -rxt(k,146)*y(k,219)
         mat(k,424) = -rxt(k,150)*y(k,219)
         mat(k,1776) = -rxt(k,162)*y(k,219)
         mat(k,496) = -rxt(k,163)*y(k,219)
         mat(k,1964) = -rxt(k,171)*y(k,219)
         mat(k,1468) = -rxt(k,172)*y(k,219)
         mat(k,903) = -rxt(k,191)*y(k,219)
         mat(k,1719) = -(rxt(k,193) + rxt(k,194)) * y(k,219)
         mat(k,2133) = -rxt(k,196)*y(k,219)
         mat(k,816) = -rxt(k,199)*y(k,219)
         mat(k,2157) = -rxt(k,223)*y(k,219)
         mat(k,810) = -rxt(k,225)*y(k,219)
         mat(k,472) = -rxt(k,239)*y(k,219)
         mat(k,599) = -rxt(k,241)*y(k,219)
         mat(k,128) = -rxt(k,242)*y(k,219)
         mat(k,369) = -rxt(k,244)*y(k,219)
         mat(k,390) = -rxt(k,246)*y(k,219)
         mat(k,246) = -rxt(k,247)*y(k,219)
         mat(k,256) = -rxt(k,248)*y(k,219)
         mat(k,303) = -rxt(k,249)*y(k,219)
         mat(k,1489) = -rxt(k,258)*y(k,219)
         mat(k,801) = -rxt(k,263)*y(k,219)
         mat(k,459) = -rxt(k,264)*y(k,219)
         mat(k,1437) = -rxt(k,265)*y(k,219)
         mat(k,183) = -rxt(k,266)*y(k,219)
         mat(k,933) = -rxt(k,267)*y(k,219)
         mat(k,1126) = -rxt(k,275)*y(k,219)
         mat(k,297) = -rxt(k,277)*y(k,219)
         mat(k,265) = -rxt(k,284)*y(k,219)
         mat(k,346) = -rxt(k,285)*y(k,219)
         mat(k,308) = -rxt(k,287)*y(k,219)
         mat(k,1091) = -rxt(k,289)*y(k,219)
         mat(k,105) = -rxt(k,290)*y(k,219)
         mat(k,688) = -rxt(k,295)*y(k,219)
         mat(k,616) = -rxt(k,296)*y(k,219)
         mat(k,1101) = -rxt(k,301)*y(k,219)
         mat(k,983) = -rxt(k,302)*y(k,219)
         mat(k,530) = -rxt(k,303)*y(k,219)
         mat(k,555) = -rxt(k,304)*y(k,219)
         mat(k,414) = -rxt(k,312)*y(k,219)
         mat(k,113) = -rxt(k,313)*y(k,219)
         mat(k,1226) = -rxt(k,315)*y(k,219)
         mat(k,1146) = -rxt(k,316)*y(k,219)
         mat(k,863) = -rxt(k,317)*y(k,219)
         mat(k,547) = -rxt(k,320)*y(k,219)
         mat(k,403) = -rxt(k,324)*y(k,219)
         mat(k,1034) = -rxt(k,325)*y(k,219)
         mat(k,927) = -rxt(k,326)*y(k,219)
         mat(k,356) = -rxt(k,328)*y(k,219)
         mat(k,1160) = -rxt(k,331)*y(k,219)
         mat(k,1217) = -rxt(k,339)*y(k,219)
         mat(k,315) = -rxt(k,340)*y(k,219)
         mat(k,512) = -rxt(k,349)*y(k,219)
         mat(k,321) = -rxt(k,350)*y(k,219)
         mat(k,579) = -rxt(k,351)*y(k,219)
         mat(k,1340) = -rxt(k,353)*y(k,219)
         mat(k,641) = -rxt(k,356)*y(k,219)
         mat(k,676) = -rxt(k,360)*y(k,219)
         mat(k,230) = -rxt(k,361)*y(k,219)
         mat(k,226) = -rxt(k,362)*y(k,219)
         mat(k,350) = -rxt(k,363)*y(k,219)
         mat(k,139) = -rxt(k,364)*y(k,219)
         mat(k,589) = -rxt(k,376)*y(k,219)
         mat(k,540) = -rxt(k,377)*y(k,219)
         mat(k,408) = -rxt(k,384)*y(k,219)
         mat(k,847) = -rxt(k,386)*y(k,219)
         mat(k,742) = -rxt(k,387)*y(k,219)
         mat(k,379) = -rxt(k,388)*y(k,219)
         mat(k,1081) = -rxt(k,389)*y(k,219)
         mat(k,207) = -rxt(k,395)*y(k,219)
         mat(k,168) = -rxt(k,398)*y(k,219)
         mat(k,421) = -rxt(k,401)*y(k,219)
         mat(k,239) = -rxt(k,402)*y(k,219)
         mat(k,341) = -rxt(k,404)*y(k,219)
         mat(k,270) = -rxt(k,408)*y(k,219)
         mat(k,199) = -rxt(k,409)*y(k,219)
         mat(k,177) = -rxt(k,421)*y(k,219)
         mat(k,335) = -rxt(k,424)*y(k,219)
         mat(k,665) = -rxt(k,428)*y(k,219)
         mat(k,194) = -rxt(k,429)*y(k,219)
         mat(k,216) = -rxt(k,431)*y(k,219)
         mat(k,704) = -rxt(k,434)*y(k,219)
         mat(k,221) = -rxt(k,437)*y(k,219)
         mat(k,433) = -rxt(k,438)*y(k,219)
         mat(k,1012) = -rxt(k,441)*y(k,219)
         mat(k,962) = -rxt(k,444)*y(k,219)
         mat(k,397) = -rxt(k,449)*y(k,219)
         mat(k,652) = -rxt(k,453)*y(k,219)
         mat(k,622) = -rxt(k,454)*y(k,219)
         mat(k,481) = -rxt(k,458)*y(k,219)
         mat(k,976) = -rxt(k,460)*y(k,219)
         mat(k,1066) = -rxt(k,461)*y(k,219)
         mat(k,287) = -(rxt(k,463) + rxt(k,477)) * y(k,219)
         mat(k,365) = -rxt(k,465)*y(k,219)
         mat(k,855) = -rxt(k,467)*y(k,219)
         mat(k,516) = -rxt(k,471)*y(k,219)
         mat(k,1237) = -rxt(k,476)*y(k,219)
         mat(k,99) = -rxt(k,479)*y(k,219)
         mat(k,1012) = mat(k,1012) + .630_r8*rxt(k,440)*y(k,135)
         mat(k,297) = mat(k,297) + .650_r8*rxt(k,277)*y(k,219)
         mat(k,555) = mat(k,555) + .130_r8*rxt(k,279)*y(k,135)
         mat(k,346) = mat(k,346) + .500_r8*rxt(k,285)*y(k,219)
         mat(k,1034) = mat(k,1034) + .360_r8*rxt(k,308)*y(k,135)
         mat(k,1489) = mat(k,1489) + rxt(k,257)*y(k,134)
         mat(k,459) = mat(k,459) + .300_r8*rxt(k,264)*y(k,219)
         mat(k,1437) = mat(k,1437) + rxt(k,271)*y(k,218)
         mat(k,2003) = rxt(k,180)*y(k,205)
         mat(k,871) = rxt(k,234)*y(k,231)
         mat(k,1451) = rxt(k,141)*y(k,135) + 2.000_r8*rxt(k,136)*y(k,205)
         mat(k,1408) = mat(k,1408) + rxt(k,133)*y(k,134) + rxt(k,125)*y(k,218)
         mat(k,606) = mat(k,606) + rxt(k,134)*y(k,134)
         mat(k,810) = mat(k,810) + rxt(k,224)*y(k,134) + rxt(k,230)*y(k,218)
         mat(k,2133) = mat(k,2133) + rxt(k,195)*y(k,134) + rxt(k,207)*y(k,218)
         mat(k,183) = mat(k,183) + rxt(k,274)*y(k,218)
         mat(k,782) = rxt(k,226)*y(k,134)
         mat(k,816) = mat(k,816) + rxt(k,198)*y(k,134)
         mat(k,847) = mat(k,847) + .320_r8*rxt(k,385)*y(k,135)
         mat(k,742) = mat(k,742) + .600_r8*rxt(k,387)*y(k,219)
         mat(k,1217) = mat(k,1217) + .240_r8*rxt(k,338)*y(k,135)
         mat(k,315) = mat(k,315) + .100_r8*rxt(k,340)*y(k,219)
         mat(k,962) = mat(k,962) + .630_r8*rxt(k,443)*y(k,135)
         mat(k,1340) = mat(k,1340) + .360_r8*rxt(k,352)*y(k,135)
         mat(k,1920) = rxt(k,164)*y(k,205)
         mat(k,1776) = mat(k,1776) + rxt(k,159)*y(k,205)
         mat(k,2188) = mat(k,2188) + rxt(k,257)*y(k,42) + rxt(k,133)*y(k,77) &
                      + rxt(k,134)*y(k,79) + rxt(k,224)*y(k,81) + rxt(k,195)*y(k,85) &
                      + rxt(k,226)*y(k,91) + rxt(k,198)*y(k,92) + rxt(k,139)*y(k,205)
         mat(k,2249) = mat(k,2249) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,141)*y(k,76) &
                      + .320_r8*rxt(k,385)*y(k,98) + .240_r8*rxt(k,338)*y(k,105) &
                      + .630_r8*rxt(k,443)*y(k,110) + .360_r8*rxt(k,352)*y(k,111) &
                      + rxt(k,140)*y(k,205)
         mat(k,547) = mat(k,547) + .500_r8*rxt(k,320)*y(k,219)
         mat(k,207) = mat(k,207) + .500_r8*rxt(k,395)*y(k,219)
         mat(k,522) = .400_r8*rxt(k,396)*y(k,205)
         mat(k,1392) = .450_r8*rxt(k,293)*y(k,205)
         mat(k,764) = .400_r8*rxt(k,410)*y(k,205)
         mat(k,2110) = mat(k,2110) + rxt(k,180)*y(k,56) + 2.000_r8*rxt(k,136)*y(k,76) &
                      + rxt(k,164)*y(k,124) + rxt(k,159)*y(k,126) + rxt(k,139) &
                      *y(k,134) + rxt(k,140)*y(k,135) + .400_r8*rxt(k,396)*y(k,190) &
                      + .450_r8*rxt(k,293)*y(k,199) + .400_r8*rxt(k,410)*y(k,201) &
                      + .450_r8*rxt(k,343)*y(k,213) + .400_r8*rxt(k,416)*y(k,214) &
                      + .200_r8*rxt(k,347)*y(k,215) + .150_r8*rxt(k,322)*y(k,222)
         mat(k,1360) = .450_r8*rxt(k,343)*y(k,205)
         mat(k,879) = .400_r8*rxt(k,416)*y(k,205)
         mat(k,682) = .200_r8*rxt(k,347)*y(k,205)
         mat(k,1529) = rxt(k,271)*y(k,54) + rxt(k,125)*y(k,77) + rxt(k,230)*y(k,81) &
                      + rxt(k,207)*y(k,85) + rxt(k,274)*y(k,86) + 2.000_r8*rxt(k,126) &
                      *y(k,231)
         mat(k,1693) = mat(k,1693) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,264)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,103) + .100_r8*rxt(k,340)*y(k,106) + .500_r8*rxt(k,320) &
                      *y(k,147) + .500_r8*rxt(k,395)*y(k,183)
         mat(k,1136) = .150_r8*rxt(k,322)*y(k,205)
         mat(k,2275) = rxt(k,234)*y(k,73) + 2.000_r8*rxt(k,126)*y(k,218)
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
         mat(k,450) = -(rxt(k,419)*y(k,205) + rxt(k,420)*y(k,124))
         mat(k,2046) = -rxt(k,419)*y(k,220)
         mat(k,1861) = -rxt(k,420)*y(k,220)
         mat(k,197) = .200_r8*rxt(k,409)*y(k,219)
         mat(k,175) = .140_r8*rxt(k,421)*y(k,219)
         mat(k,333) = rxt(k,424)*y(k,219)
         mat(k,1605) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,143) &
                      + rxt(k,424)*y(k,144)
         mat(k,770) = -(rxt(k,318)*y(k,205) + rxt(k,319)*y(k,124))
         mat(k,2070) = -rxt(k,318)*y(k,221)
         mat(k,1882) = -rxt(k,319)*y(k,221)
         mat(k,1022) = rxt(k,325)*y(k,219)
         mat(k,543) = .500_r8*rxt(k,320)*y(k,219)
         mat(k,1643) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,147)
         mat(k,1131) = -(rxt(k,321)*y(k,200) + rxt(k,322)*y(k,205) + rxt(k,323) &
                      *y(k,124))
         mat(k,1811) = -rxt(k,321)*y(k,222)
         mat(k,2091) = -rxt(k,322)*y(k,222)
         mat(k,1901) = -rxt(k,323)*y(k,222)
         mat(k,1007) = .060_r8*rxt(k,440)*y(k,135)
         mat(k,925) = rxt(k,326)*y(k,219)
         mat(k,957) = .060_r8*rxt(k,443)*y(k,135)
         mat(k,2231) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,110)
         mat(k,400) = rxt(k,324)*y(k,219)
         mat(k,1063) = .150_r8*rxt(k,461)*y(k,219)
         mat(k,1672) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,148) + .150_r8*rxt(k,461) &
                      *y(k,180)
         mat(k,1111) = -(rxt(k,450)*y(k,200) + rxt(k,451)*y(k,205) + rxt(k,452) &
                      *y(k,124))
         mat(k,1809) = -rxt(k,450)*y(k,223)
         mat(k,2089) = -rxt(k,451)*y(k,223)
         mat(k,1899) = -rxt(k,452)*y(k,223)
         mat(k,1754) = .500_r8*rxt(k,459)*y(k,179)
         mat(k,649) = rxt(k,453)*y(k,219)
         mat(k,974) = .500_r8*rxt(k,459)*y(k,126) + rxt(k,460)*y(k,219)
         mat(k,1670) = rxt(k,453)*y(k,176) + rxt(k,460)*y(k,179)
         mat(k,914) = -(rxt(k,455)*y(k,200) + rxt(k,456)*y(k,205) + rxt(k,457) &
                      *y(k,124))
         mat(k,1800) = -rxt(k,455)*y(k,224)
         mat(k,2080) = -rxt(k,456)*y(k,224)
         mat(k,1889) = -rxt(k,457)*y(k,224)
         mat(k,996) = rxt(k,441)*y(k,219)
         mat(k,947) = rxt(k,444)*y(k,219)
         mat(k,477) = rxt(k,458)*y(k,219)
         mat(k,1657) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) + rxt(k,458)*y(k,178)
         mat(k,726) = -(rxt(k,426)*y(k,205) + rxt(k,427)*y(k,124))
         mat(k,2066) = -rxt(k,426)*y(k,225)
         mat(k,1879) = -rxt(k,427)*y(k,225)
         mat(k,659) = rxt(k,428)*y(k,219)
         mat(k,193) = .650_r8*rxt(k,429)*y(k,219)
         mat(k,1639) = rxt(k,428)*y(k,181) + .650_r8*rxt(k,429)*y(k,182)
         mat(k,88) = -(rxt(k,517)*y(k,205) + rxt(k,518)*y(k,124))
         mat(k,2026) = -rxt(k,517)*y(k,226)
         mat(k,1851) = -rxt(k,518)*y(k,226)
         mat(k,188) = rxt(k,516)*y(k,219)
         mat(k,1554) = rxt(k,516)*y(k,182)
         mat(k,1175) = -(rxt(k,390)*y(k,199) + rxt(k,391)*y(k,200) + rxt(k,392) &
                      *y(k,205) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126))
         mat(k,1379) = -rxt(k,390)*y(k,227)
         mat(k,1813) = -rxt(k,391)*y(k,227)
         mat(k,2093) = -rxt(k,392)*y(k,227)
         mat(k,1904) = -rxt(k,393)*y(k,227)
         mat(k,1759) = -rxt(k,394)*y(k,227)
         mat(k,225) = rxt(k,362)*y(k,219)
         mat(k,349) = rxt(k,363)*y(k,219)
         mat(k,138) = rxt(k,364)*y(k,219)
         mat(k,738) = .400_r8*rxt(k,387)*y(k,219)
         mat(k,206) = .500_r8*rxt(k,395)*y(k,219)
         mat(k,1675) = rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364)*y(k,97) &
                      + .400_r8*rxt(k,387)*y(k,103) + .500_r8*rxt(k,395)*y(k,183)
         mat(k,750) = -(rxt(k,432)*y(k,205) + rxt(k,433)*y(k,124))
         mat(k,2068) = -rxt(k,432)*y(k,228)
         mat(k,1880) = -rxt(k,433)*y(k,228)
         mat(k,213) = .560_r8*rxt(k,431)*y(k,219)
         mat(k,697) = rxt(k,434)*y(k,219)
         mat(k,1641) = .560_r8*rxt(k,431)*y(k,184) + rxt(k,434)*y(k,185)
         mat(k,94) = -(rxt(k,520)*y(k,205) + rxt(k,521)*y(k,124))
         mat(k,2027) = -rxt(k,520)*y(k,229)
         mat(k,1852) = -rxt(k,521)*y(k,229)
         mat(k,208) = rxt(k,519)*y(k,219)
         mat(k,1555) = rxt(k,519)*y(k,184)
         mat(k,501) = -(rxt(k,435)*y(k,205) + rxt(k,436)*y(k,124))
         mat(k,2052) = -rxt(k,435)*y(k,230)
         mat(k,1866) = -rxt(k,436)*y(k,230)
         mat(k,220) = .300_r8*rxt(k,437)*y(k,219)
         mat(k,430) = rxt(k,438)*y(k,219)
         mat(k,1613) = .300_r8*rxt(k,437)*y(k,186) + rxt(k,438)*y(k,187)
         mat(k,2287) = -(rxt(k,126)*y(k,218) + rxt(k,234)*y(k,73) + rxt(k,478) &
                      *y(k,153))
         mat(k,1541) = -rxt(k,126)*y(k,231)
         mat(k,874) = -rxt(k,234)*y(k,231)
         mat(k,262) = -rxt(k,478)*y(k,231)
         mat(k,311) = rxt(k,287)*y(k,219)
         mat(k,416) = rxt(k,312)*y(k,219)
         mat(k,114) = rxt(k,313)*y(k,219)
         mat(k,475) = rxt(k,239)*y(k,219)
         mat(k,1500) = rxt(k,258)*y(k,219)
         mat(k,604) = rxt(k,241)*y(k,219)
         mat(k,130) = rxt(k,242)*y(k,219)
         mat(k,1095) = rxt(k,289)*y(k,219)
         mat(k,374) = rxt(k,244)*y(k,219)
         mat(k,929) = rxt(k,326)*y(k,219)
         mat(k,1230) = rxt(k,315)*y(k,219)
         mat(k,690) = rxt(k,295)*y(k,219)
         mat(k,618) = rxt(k,296)*y(k,219)
         mat(k,461) = rxt(k,264)*y(k,219)
         mat(k,1444) = rxt(k,265)*y(k,219)
         mat(k,1459) = rxt(k,137)*y(k,205)
         mat(k,1414) = rxt(k,142)*y(k,219)
         mat(k,611) = rxt(k,143)*y(k,219)
         mat(k,813) = rxt(k,225)*y(k,219)
         mat(k,305) = rxt(k,249)*y(k,219)
         mat(k,2145) = (rxt(k,530)+rxt(k,535))*y(k,91) + (rxt(k,523)+rxt(k,529) &
                       +rxt(k,534))*y(k,92) + rxt(k,196)*y(k,219)
         mat(k,936) = rxt(k,267)*y(k,219)
         mat(k,1477) = rxt(k,172)*y(k,219)
         mat(k,428) = rxt(k,150)*y(k,219)
         mat(k,787) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,821) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85) + rxt(k,199)*y(k,219)
         mat(k,1221) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,100) = rxt(k,479)*y(k,219)
         mat(k,549) = rxt(k,320)*y(k,219)
         mat(k,404) = rxt(k,324)*y(k,219)
         mat(k,2122) = rxt(k,137)*y(k,76) + rxt(k,144)*y(k,219)
         mat(k,1705) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31) &
                      + rxt(k,239)*y(k,41) + rxt(k,258)*y(k,42) + rxt(k,241)*y(k,43) &
                      + rxt(k,242)*y(k,44) + rxt(k,289)*y(k,45) + rxt(k,244)*y(k,46) &
                      + rxt(k,326)*y(k,48) + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50) &
                      + rxt(k,296)*y(k,51) + rxt(k,264)*y(k,53) + rxt(k,265)*y(k,54) &
                      + rxt(k,142)*y(k,77) + rxt(k,143)*y(k,79) + rxt(k,225)*y(k,81) &
                      + rxt(k,249)*y(k,84) + rxt(k,196)*y(k,85) + rxt(k,267)*y(k,87) &
                      + rxt(k,172)*y(k,89) + rxt(k,150)*y(k,90) + rxt(k,199)*y(k,92) &
                      + .500_r8*rxt(k,339)*y(k,105) + rxt(k,479)*y(k,120) + rxt(k,320) &
                      *y(k,147) + rxt(k,324)*y(k,148) + rxt(k,144)*y(k,205) &
                      + 2.000_r8*rxt(k,147)*y(k,219)
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
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 118) = mat(k, 118) + lmat(k, 118)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 129) = mat(k, 129) + lmat(k, 129)
         mat(k, 131) = lmat(k, 131)
         mat(k, 132) = lmat(k, 132)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = lmat(k, 136)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = lmat(k, 141)
         mat(k, 142) = lmat(k, 142)
         mat(k, 143) = mat(k, 143) + lmat(k, 143)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 145) = mat(k, 145) + lmat(k, 145)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 160) = mat(k, 160) + lmat(k, 160)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 170) = lmat(k, 170)
         mat(k, 171) = lmat(k, 171)
         mat(k, 172) = lmat(k, 172)
         mat(k, 173) = lmat(k, 173)
         mat(k, 174) = mat(k, 174) + lmat(k, 174)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = lmat(k, 180)
         mat(k, 181) = mat(k, 181) + lmat(k, 181)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 186) = mat(k, 186) + lmat(k, 186)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 201) = lmat(k, 201)
         mat(k, 202) = lmat(k, 202)
         mat(k, 203) = lmat(k, 203)
         mat(k, 204) = lmat(k, 204)
         mat(k, 205) = mat(k, 205) + lmat(k, 205)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 210) = mat(k, 210) + lmat(k, 210)
         mat(k, 218) = mat(k, 218) + lmat(k, 218)
         mat(k, 223) = mat(k, 223) + lmat(k, 223)
         mat(k, 224) = lmat(k, 224)
         mat(k, 226) = mat(k, 226) + lmat(k, 226)
         mat(k, 227) = lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 231) = lmat(k, 231)
         mat(k, 232) = lmat(k, 232)
         mat(k, 233) = lmat(k, 233)
         mat(k, 234) = lmat(k, 234)
         mat(k, 235) = lmat(k, 235)
         mat(k, 236) = lmat(k, 236)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 240) = mat(k, 240) + lmat(k, 240)
         mat(k, 241) = mat(k, 241) + lmat(k, 241)
         mat(k, 243) = mat(k, 243) + lmat(k, 243)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = mat(k, 249) + lmat(k, 249)
         mat(k, 251) = mat(k, 251) + lmat(k, 251)
         mat(k, 252) = mat(k, 252) + lmat(k, 252)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 260) = lmat(k, 260)
         mat(k, 261) = lmat(k, 261)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
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
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = mat(k, 293) + lmat(k, 293)
         mat(k, 299) = mat(k, 299) + lmat(k, 299)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 317) = mat(k, 317) + lmat(k, 317)
         mat(k, 319) = lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 322) = lmat(k, 322)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = lmat(k, 325)
         mat(k, 326) = lmat(k, 326)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 334) = lmat(k, 334)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 336) = lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 338) = mat(k, 338) + lmat(k, 338)
         mat(k, 341) = mat(k, 341) + lmat(k, 341)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 346) = mat(k, 346) + lmat(k, 346)
         mat(k, 347) = lmat(k, 347)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 351) = mat(k, 351) + lmat(k, 351)
         mat(k, 359) = mat(k, 359) + lmat(k, 359)
         mat(k, 360) = lmat(k, 360)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 367) = mat(k, 367) + lmat(k, 367)
         mat(k, 370) = lmat(k, 370)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 378) = lmat(k, 378)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 387) = mat(k, 387) + lmat(k, 387)
         mat(k, 388) = mat(k, 388) + lmat(k, 388)
         mat(k, 393) = mat(k, 393) + lmat(k, 393)
         mat(k, 394) = lmat(k, 394)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 401) = lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 405) = mat(k, 405) + lmat(k, 405)
         mat(k, 406) = lmat(k, 406)
         mat(k, 409) = lmat(k, 409)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 411) = mat(k, 411) + lmat(k, 411)
         mat(k, 413) = lmat(k, 413)
         mat(k, 414) = mat(k, 414) + lmat(k, 414)
         mat(k, 415) = lmat(k, 415)
         mat(k, 417) = mat(k, 417) + lmat(k, 417)
         mat(k, 418) = lmat(k, 418)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 422) = lmat(k, 422)
         mat(k, 423) = mat(k, 423) + lmat(k, 423)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 425) = lmat(k, 425)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 427) = lmat(k, 427)
         mat(k, 429) = mat(k, 429) + lmat(k, 429)
         mat(k, 431) = lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 433) = mat(k, 433) + lmat(k, 433)
         mat(k, 434) = lmat(k, 434)
         mat(k, 437) = mat(k, 437) + lmat(k, 437)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 445) = lmat(k, 445)
         mat(k, 448) = mat(k, 448) + lmat(k, 448)
         mat(k, 450) = mat(k, 450) + lmat(k, 450)
         mat(k, 456) = mat(k, 456) + lmat(k, 456)
         mat(k, 457) = lmat(k, 457)
         mat(k, 458) = mat(k, 458) + lmat(k, 458)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = lmat(k, 464)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 470) = mat(k, 470) + lmat(k, 470)
         mat(k, 476) = mat(k, 476) + lmat(k, 476)
         mat(k, 478) = lmat(k, 478)
         mat(k, 479) = lmat(k, 479)
         mat(k, 480) = lmat(k, 480)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 485) = mat(k, 485) + lmat(k, 485)
         mat(k, 493) = mat(k, 493) + lmat(k, 493)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 499) = mat(k, 499) + lmat(k, 499)
         mat(k, 501) = mat(k, 501) + lmat(k, 501)
         mat(k, 508) = mat(k, 508) + lmat(k, 508)
         mat(k, 510) = lmat(k, 510)
         mat(k, 511) = lmat(k, 511)
         mat(k, 513) = mat(k, 513) + lmat(k, 513)
         mat(k, 514) = mat(k, 514) + lmat(k, 514)
         mat(k, 517) = lmat(k, 517)
         mat(k, 520) = mat(k, 520) + lmat(k, 520)
         mat(k, 526) = mat(k, 526) + lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 528) = lmat(k, 528)
         mat(k, 531) = mat(k, 531) + lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 539) = lmat(k, 539)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 544) = lmat(k, 544)
         mat(k, 546) = lmat(k, 546)
         mat(k, 547) = mat(k, 547) + lmat(k, 547)
         mat(k, 548) = lmat(k, 548)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 559) = lmat(k, 559)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = lmat(k, 561)
         mat(k, 562) = mat(k, 562) + lmat(k, 562)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = mat(k, 564) + lmat(k, 564)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 577) = lmat(k, 577)
         mat(k, 581) = lmat(k, 581)
         mat(k, 583) = mat(k, 583) + lmat(k, 583)
         mat(k, 587) = lmat(k, 587)
         mat(k, 592) = lmat(k, 592)
         mat(k, 593) = lmat(k, 593)
         mat(k, 594) = lmat(k, 594)
         mat(k, 595) = lmat(k, 595)
         mat(k, 596) = mat(k, 596) + lmat(k, 596)
         mat(k, 597) = mat(k, 597) + lmat(k, 597)
         mat(k, 600) = lmat(k, 600)
         mat(k, 605) = mat(k, 605) + lmat(k, 605)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 612) = mat(k, 612) + lmat(k, 612)
         mat(k, 613) = mat(k, 613) + lmat(k, 613)
         mat(k, 616) = mat(k, 616) + lmat(k, 616)
         mat(k, 617) = lmat(k, 617)
         mat(k, 619) = mat(k, 619) + lmat(k, 619)
         mat(k, 620) = mat(k, 620) + lmat(k, 620)
         mat(k, 621) = lmat(k, 621)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 624) = lmat(k, 624)
         mat(k, 627) = mat(k, 627) + lmat(k, 627)
         mat(k, 633) = lmat(k, 633)
         mat(k, 634) = mat(k, 634) + lmat(k, 634)
         mat(k, 637) = mat(k, 637) + lmat(k, 637)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 640) = mat(k, 640) + lmat(k, 640)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 643) = lmat(k, 643)
         mat(k, 644) = mat(k, 644) + lmat(k, 644)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = lmat(k, 646)
         mat(k, 647) = lmat(k, 647)
         mat(k, 648) = lmat(k, 648)
         mat(k, 650) = lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = mat(k, 652) + lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 654) = lmat(k, 654)
         mat(k, 655) = lmat(k, 655)
         mat(k, 656) = lmat(k, 656)
         mat(k, 657) = mat(k, 657) + lmat(k, 657)
         mat(k, 662) = lmat(k, 662)
         mat(k, 664) = lmat(k, 664)
         mat(k, 665) = mat(k, 665) + lmat(k, 665)
         mat(k, 666) = lmat(k, 666)
         mat(k, 667) = lmat(k, 667)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 672) = lmat(k, 672)
         mat(k, 673) = lmat(k, 673)
         mat(k, 675) = lmat(k, 675)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 677) = lmat(k, 677)
         mat(k, 679) = mat(k, 679) + lmat(k, 679)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 691) = lmat(k, 691)
         mat(k, 692) = lmat(k, 692)
         mat(k, 693) = lmat(k, 693)
         mat(k, 694) = lmat(k, 694)
         mat(k, 695) = mat(k, 695) + lmat(k, 695)
         mat(k, 700) = lmat(k, 700)
         mat(k, 702) = lmat(k, 702)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 705) = lmat(k, 705)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 726) = mat(k, 726) + lmat(k, 726)
         mat(k, 737) = mat(k, 737) + lmat(k, 737)
         mat(k, 739) = lmat(k, 739)
         mat(k, 740) = lmat(k, 740)
         mat(k, 741) = lmat(k, 741)
         mat(k, 742) = mat(k, 742) + lmat(k, 742)
         mat(k, 743) = lmat(k, 743)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 761) = mat(k, 761) + lmat(k, 761)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 781) = lmat(k, 781)
         mat(k, 782) = mat(k, 782) + lmat(k, 782)
         mat(k, 789) = mat(k, 789) + lmat(k, 789)
         mat(k, 799) = mat(k, 799) + lmat(k, 799)
         mat(k, 803) = lmat(k, 803)
         mat(k, 804) = lmat(k, 804)
         mat(k, 805) = lmat(k, 805)
         mat(k, 806) = mat(k, 806) + lmat(k, 806)
         mat(k, 807) = mat(k, 807) + lmat(k, 807)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 816) = mat(k, 816) + lmat(k, 816)
         mat(k, 818) = mat(k, 818) + lmat(k, 818)
         mat(k, 825) = mat(k, 825) + lmat(k, 825)
         mat(k, 836) = mat(k, 836) + lmat(k, 836)
         mat(k, 853) = mat(k, 853) + lmat(k, 853)
         mat(k, 854) = lmat(k, 854)
         mat(k, 857) = lmat(k, 857)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 861) = lmat(k, 861)
         mat(k, 862) = lmat(k, 862)
         mat(k, 864) = mat(k, 864) + lmat(k, 864)
         mat(k, 866) = mat(k, 866) + lmat(k, 866)
         mat(k, 875) = mat(k, 875) + lmat(k, 875)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 900) = mat(k, 900) + lmat(k, 900)
         mat(k, 901) = mat(k, 901) + lmat(k, 901)
         mat(k, 902) = mat(k, 902) + lmat(k, 902)
         mat(k, 904) = mat(k, 904) + lmat(k, 904)
         mat(k, 905) = mat(k, 905) + lmat(k, 905)
         mat(k, 906) = lmat(k, 906)
         mat(k, 907) = mat(k, 907) + lmat(k, 907)
         mat(k, 914) = mat(k, 914) + lmat(k, 914)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 926) = lmat(k, 926)
         mat(k, 928) = lmat(k, 928)
         mat(k, 931) = mat(k, 931) + lmat(k, 931)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 971) = mat(k, 971) + lmat(k, 971)
         mat(k, 973) = lmat(k, 973)
         mat(k, 975) = lmat(k, 975)
         mat(k, 978) = lmat(k, 978)
         mat(k, 979) = lmat(k, 979)
         mat(k, 981) = mat(k, 981) + lmat(k, 981)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 984) = mat(k, 984) + lmat(k, 984)
         mat(k,1001) = mat(k,1001) + lmat(k,1001)
         mat(k,1026) = mat(k,1026) + lmat(k,1026)
         mat(k,1047) = mat(k,1047) + lmat(k,1047)
         mat(k,1059) = mat(k,1059) + lmat(k,1059)
         mat(k,1060) = mat(k,1060) + lmat(k,1060)
         mat(k,1061) = mat(k,1061) + lmat(k,1061)
         mat(k,1062) = mat(k,1062) + lmat(k,1062)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1064) = mat(k,1064) + lmat(k,1064)
         mat(k,1065) = mat(k,1065) + lmat(k,1065)
         mat(k,1068) = mat(k,1068) + lmat(k,1068)
         mat(k,1071) = lmat(k,1071)
         mat(k,1075) = mat(k,1075) + lmat(k,1075)
         mat(k,1079) = lmat(k,1079)
         mat(k,1084) = lmat(k,1084)
         mat(k,1085) = mat(k,1085) + lmat(k,1085)
         mat(k,1087) = mat(k,1087) + lmat(k,1087)
         mat(k,1088) = lmat(k,1088)
         mat(k,1093) = lmat(k,1093)
         mat(k,1094) = lmat(k,1094)
         mat(k,1098) = mat(k,1098) + lmat(k,1098)
         mat(k,1099) = lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1102) = mat(k,1102) + lmat(k,1102)
         mat(k,1111) = mat(k,1111) + lmat(k,1111)
         mat(k,1125) = mat(k,1125) + lmat(k,1125)
         mat(k,1131) = mat(k,1131) + lmat(k,1131)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1144) = lmat(k,1144)
         mat(k,1145) = lmat(k,1145)
         mat(k,1147) = mat(k,1147) + lmat(k,1147)
         mat(k,1149) = lmat(k,1149)
         mat(k,1150) = lmat(k,1150)
         mat(k,1151) = lmat(k,1151)
         mat(k,1152) = lmat(k,1152)
         mat(k,1154) = lmat(k,1154)
         mat(k,1155) = mat(k,1155) + lmat(k,1155)
         mat(k,1157) = lmat(k,1157)
         mat(k,1158) = lmat(k,1158)
         mat(k,1159) = lmat(k,1159)
         mat(k,1164) = lmat(k,1164)
         mat(k,1165) = mat(k,1165) + lmat(k,1165)
         mat(k,1175) = mat(k,1175) + lmat(k,1175)
         mat(k,1195) = mat(k,1195) + lmat(k,1195)
         mat(k,1210) = mat(k,1210) + lmat(k,1210)
         mat(k,1211) = mat(k,1211) + lmat(k,1211)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1215) = mat(k,1215) + lmat(k,1215)
         mat(k,1216) = mat(k,1216) + lmat(k,1216)
         mat(k,1218) = mat(k,1218) + lmat(k,1218)
         mat(k,1222) = mat(k,1222) + lmat(k,1222)
         mat(k,1223) = mat(k,1223) + lmat(k,1223)
         mat(k,1224) = mat(k,1224) + lmat(k,1224)
         mat(k,1228) = lmat(k,1228)
         mat(k,1232) = lmat(k,1232)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1244) = lmat(k,1244)
         mat(k,1258) = mat(k,1258) + lmat(k,1258)
         mat(k,1274) = lmat(k,1274)
         mat(k,1290) = mat(k,1290) + lmat(k,1290)
         mat(k,1302) = mat(k,1302) + lmat(k,1302)
         mat(k,1313) = mat(k,1313) + lmat(k,1313)
         mat(k,1328) = lmat(k,1328)
         mat(k,1330) = mat(k,1330) + lmat(k,1330)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1336) = mat(k,1336) + lmat(k,1336)
         mat(k,1342) = lmat(k,1342)
         mat(k,1356) = mat(k,1356) + lmat(k,1356)
         mat(k,1388) = mat(k,1388) + lmat(k,1388)
         mat(k,1403) = mat(k,1403) + lmat(k,1403)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1428) = lmat(k,1428)
         mat(k,1430) = lmat(k,1430)
         mat(k,1431) = mat(k,1431) + lmat(k,1431)
         mat(k,1432) = mat(k,1432) + lmat(k,1432)
         mat(k,1433) = mat(k,1433) + lmat(k,1433)
         mat(k,1435) = mat(k,1435) + lmat(k,1435)
         mat(k,1437) = mat(k,1437) + lmat(k,1437)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1443) = lmat(k,1443)
         mat(k,1444) = mat(k,1444) + lmat(k,1444)
         mat(k,1447) = mat(k,1447) + lmat(k,1447)
         mat(k,1455) = mat(k,1455) + lmat(k,1455)
         mat(k,1465) = mat(k,1465) + lmat(k,1465)
         mat(k,1468) = mat(k,1468) + lmat(k,1468)
         mat(k,1471) = lmat(k,1471)
         mat(k,1481) = mat(k,1481) + lmat(k,1481)
         mat(k,1482) = lmat(k,1482)
         mat(k,1485) = mat(k,1485) + lmat(k,1485)
         mat(k,1487) = mat(k,1487) + lmat(k,1487)
         mat(k,1528) = mat(k,1528) + lmat(k,1528)
         mat(k,1539) = lmat(k,1539)
         mat(k,1693) = mat(k,1693) + lmat(k,1693)
         mat(k,1720) = mat(k,1720) + lmat(k,1720)
         mat(k,1725) = mat(k,1725) + lmat(k,1725)
         mat(k,1729) = mat(k,1729) + lmat(k,1729)
         mat(k,1773) = mat(k,1773) + lmat(k,1773)
         mat(k,1778) = mat(k,1778) + lmat(k,1778)
         mat(k,1780) = mat(k,1780) + lmat(k,1780)
         mat(k,1781) = mat(k,1781) + lmat(k,1781)
         mat(k,1786) = mat(k,1786) + lmat(k,1786)
         mat(k,1831) = mat(k,1831) + lmat(k,1831)
         mat(k,1865) = mat(k,1865) + lmat(k,1865)
         mat(k,1924) = mat(k,1924) + lmat(k,1924)
         mat(k,1930) = mat(k,1930) + lmat(k,1930)
         mat(k,1961) = mat(k,1961) + lmat(k,1961)
         mat(k,1964) = mat(k,1964) + lmat(k,1964)
         mat(k,1968) = mat(k,1968) + lmat(k,1968)
         mat(k,1969) = mat(k,1969) + lmat(k,1969)
         mat(k,1974) = mat(k,1974) + lmat(k,1974)
         mat(k,2009) = mat(k,2009) + lmat(k,2009)
         mat(k,2117) = mat(k,2117) + lmat(k,2117)
         mat(k,2122) = mat(k,2122) + lmat(k,2122)
         mat(k,2129) = mat(k,2129) + lmat(k,2129)
         mat(k,2139) = mat(k,2139) + lmat(k,2139)
         mat(k,2141) = mat(k,2141) + lmat(k,2141)
         mat(k,2152) = mat(k,2152) + lmat(k,2152)
         mat(k,2166) = mat(k,2166) + lmat(k,2166)
         mat(k,2167) = mat(k,2167) + lmat(k,2167)
         mat(k,2198) = mat(k,2198) + lmat(k,2198)
         mat(k,2199) = mat(k,2199) + lmat(k,2199)
         mat(k,2248) = mat(k,2248) + lmat(k,2248)
         mat(k,2259) = mat(k,2259) + lmat(k,2259)
         mat(k,2260) = mat(k,2260) + lmat(k,2260)
         mat(k,2268) = lmat(k,2268)
         mat(k,2271) = lmat(k,2271)
         mat(k,2274) = mat(k,2274) + lmat(k,2274)
         mat(k,2275) = mat(k,2275) + lmat(k,2275)
         mat(k,2285) = lmat(k,2285)
         mat(k,2287) = mat(k,2287) + lmat(k,2287)
         mat(k, 214) = 0._r8
         mat(k, 215) = 0._r8
         mat(k, 254) = 0._r8
         mat(k, 301) = 0._r8
         mat(k, 340) = 0._r8
         mat(k, 438) = 0._r8
         mat(k, 439) = 0._r8
         mat(k, 452) = 0._r8
         mat(k, 486) = 0._r8
         mat(k, 488) = 0._r8
         mat(k, 495) = 0._r8
         mat(k, 504) = 0._r8
         mat(k, 628) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 636) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 658) = 0._r8
         mat(k, 660) = 0._r8
         mat(k, 661) = 0._r8
         mat(k, 663) = 0._r8
         mat(k, 669) = 0._r8
         mat(k, 670) = 0._r8
         mat(k, 674) = 0._r8
         mat(k, 696) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 701) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 709) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 725) = 0._r8
         mat(k, 727) = 0._r8
         mat(k, 728) = 0._r8
         mat(k, 730) = 0._r8
         mat(k, 732) = 0._r8
         mat(k, 749) = 0._r8
         mat(k, 751) = 0._r8
         mat(k, 752) = 0._r8
         mat(k, 754) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 757) = 0._r8
         mat(k, 772) = 0._r8
         mat(k, 774) = 0._r8
         mat(k, 778) = 0._r8
         mat(k, 783) = 0._r8
         mat(k, 793) = 0._r8
         mat(k, 798) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 858) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 888) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 894) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 910) = 0._r8
         mat(k, 917) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 935) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 952) = 0._r8
         mat(k, 955) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 985) = 0._r8
         mat(k, 995) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1004) = 0._r8
         mat(k,1005) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1010) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1016) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1020) = 0._r8
         mat(k,1029) = 0._r8
         mat(k,1030) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1038) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1048) = 0._r8
         mat(k,1049) = 0._r8
         mat(k,1050) = 0._r8
         mat(k,1051) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1113) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1115) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1122) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1141) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1162) = 0._r8
         mat(k,1163) = 0._r8
         mat(k,1166) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1194) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1202) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1219) = 0._r8
         mat(k,1229) = 0._r8
         mat(k,1242) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1251) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1271) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1273) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1286) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1304) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1316) = 0._r8
         mat(k,1318) = 0._r8
         mat(k,1324) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1335) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1338) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1348) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1355) = 0._r8
         mat(k,1358) = 0._r8
         mat(k,1366) = 0._r8
         mat(k,1368) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1400) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1409) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1418) = 0._r8
         mat(k,1420) = 0._r8
         mat(k,1421) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1434) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1446) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1449) = 0._r8
         mat(k,1450) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1453) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1462) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1464) = 0._r8
         mat(k,1466) = 0._r8
         mat(k,1467) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1474) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1484) = 0._r8
         mat(k,1488) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1493) = 0._r8
         mat(k,1497) = 0._r8
         mat(k,1499) = 0._r8
         mat(k,1526) = 0._r8
         mat(k,1531) = 0._r8
         mat(k,1534) = 0._r8
         mat(k,1607) = 0._r8
         mat(k,1624) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1642) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1654) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1716) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1721) = 0._r8
         mat(k,1730) = 0._r8
         mat(k,1731) = 0._r8
         mat(k,1738) = 0._r8
         mat(k,1743) = 0._r8
         mat(k,1750) = 0._r8
         mat(k,1756) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1763) = 0._r8
         mat(k,1770) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1775) = 0._r8
         mat(k,1777) = 0._r8
         mat(k,1779) = 0._r8
         mat(k,1782) = 0._r8
         mat(k,1784) = 0._r8
         mat(k,1785) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1788) = 0._r8
         mat(k,1799) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1824) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1827) = 0._r8
         mat(k,1828) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1836) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1838) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1840) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1915) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1917) = 0._r8
         mat(k,1919) = 0._r8
         mat(k,1928) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1943) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1952) = 0._r8
         mat(k,1953) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1963) = 0._r8
         mat(k,1967) = 0._r8
         mat(k,1970) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2007) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2013) = 0._r8
         mat(k,2015) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2049) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2081) = 0._r8
         mat(k,2083) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2092) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2100) = 0._r8
         mat(k,2105) = 0._r8
         mat(k,2107) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2128) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2138) = 0._r8
         mat(k,2140) = 0._r8
         mat(k,2142) = 0._r8
         mat(k,2144) = 0._r8
         mat(k,2153) = 0._r8
         mat(k,2154) = 0._r8
         mat(k,2155) = 0._r8
         mat(k,2156) = 0._r8
         mat(k,2159) = 0._r8
         mat(k,2160) = 0._r8
         mat(k,2165) = 0._r8
         mat(k,2168) = 0._r8
         mat(k,2169) = 0._r8
         mat(k,2171) = 0._r8
         mat(k,2177) = 0._r8
         mat(k,2183) = 0._r8
         mat(k,2185) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2191) = 0._r8
         mat(k,2200) = 0._r8
         mat(k,2213) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2225) = 0._r8
         mat(k,2228) = 0._r8
         mat(k,2229) = 0._r8
         mat(k,2232) = 0._r8
         mat(k,2233) = 0._r8
         mat(k,2237) = 0._r8
         mat(k,2238) = 0._r8
         mat(k,2239) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2246) = 0._r8
         mat(k,2257) = 0._r8
         mat(k,2261) = 0._r8
         mat(k,2265) = 0._r8
         mat(k,2267) = 0._r8
         mat(k,2269) = 0._r8
         mat(k,2270) = 0._r8
         mat(k,2272) = 0._r8
         mat(k,2273) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2277) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2280) = 0._r8
         mat(k,2281) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2283) = 0._r8
         mat(k,2284) = 0._r8
         mat(k,2286) = 0._r8
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
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 134) = mat(k, 134) - dti(k)
         mat(k, 137) = mat(k, 137) - dti(k)
         mat(k, 140) = mat(k, 140) - dti(k)
         mat(k, 143) = mat(k, 143) - dti(k)
         mat(k, 148) = mat(k, 148) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 158) = mat(k, 158) - dti(k)
         mat(k, 164) = mat(k, 164) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 179) = mat(k, 179) - dti(k)
         mat(k, 181) = mat(k, 181) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 201) = mat(k, 201) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 210) = mat(k, 210) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 223) = mat(k, 223) - dti(k)
         mat(k, 228) = mat(k, 228) - dti(k)
         mat(k, 231) = mat(k, 231) - dti(k)
         mat(k, 234) = mat(k, 234) - dti(k)
         mat(k, 237) = mat(k, 237) - dti(k)
         mat(k, 240) = mat(k, 240) - dti(k)
         mat(k, 243) = mat(k, 243) - dti(k)
         mat(k, 248) = mat(k, 248) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 263) = mat(k, 263) - dti(k)
         mat(k, 267) = mat(k, 267) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 281) = mat(k, 281) - dti(k)
         mat(k, 284) = mat(k, 284) - dti(k)
         mat(k, 290) = mat(k, 290) - dti(k)
         mat(k, 293) = mat(k, 293) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 306) = mat(k, 306) - dti(k)
         mat(k, 312) = mat(k, 312) - dti(k)
         mat(k, 317) = mat(k, 317) - dti(k)
         mat(k, 322) = mat(k, 322) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 332) = mat(k, 332) - dti(k)
         mat(k, 338) = mat(k, 338) - dti(k)
         mat(k, 343) = mat(k, 343) - dti(k)
         mat(k, 348) = mat(k, 348) - dti(k)
         mat(k, 351) = mat(k, 351) - dti(k)
         mat(k, 359) = mat(k, 359) - dti(k)
         mat(k, 367) = mat(k, 367) - dti(k)
         mat(k, 375) = mat(k, 375) - dti(k)
         mat(k, 381) = mat(k, 381) - dti(k)
         mat(k, 387) = mat(k, 387) - dti(k)
         mat(k, 393) = mat(k, 393) - dti(k)
         mat(k, 399) = mat(k, 399) - dti(k)
         mat(k, 405) = mat(k, 405) - dti(k)
         mat(k, 411) = mat(k, 411) - dti(k)
         mat(k, 417) = mat(k, 417) - dti(k)
         mat(k, 423) = mat(k, 423) - dti(k)
         mat(k, 429) = mat(k, 429) - dti(k)
         mat(k, 437) = mat(k, 437) - dti(k)
         mat(k, 443) = mat(k, 443) - dti(k)
         mat(k, 450) = mat(k, 450) - dti(k)
         mat(k, 456) = mat(k, 456) - dti(k)
         mat(k, 462) = mat(k, 462) - dti(k)
         mat(k, 465) = mat(k, 465) - dti(k)
         mat(k, 469) = mat(k, 469) - dti(k)
         mat(k, 476) = mat(k, 476) - dti(k)
         mat(k, 485) = mat(k, 485) - dti(k)
         mat(k, 493) = mat(k, 493) - dti(k)
         mat(k, 501) = mat(k, 501) - dti(k)
         mat(k, 508) = mat(k, 508) - dti(k)
         mat(k, 513) = mat(k, 513) - dti(k)
         mat(k, 520) = mat(k, 520) - dti(k)
         mat(k, 526) = mat(k, 526) - dti(k)
         mat(k, 534) = mat(k, 534) - dti(k)
         mat(k, 542) = mat(k, 542) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 558) = mat(k, 558) - dti(k)
         mat(k, 566) = mat(k, 566) - dti(k)
         mat(k, 574) = mat(k, 574) - dti(k)
         mat(k, 583) = mat(k, 583) - dti(k)
         mat(k, 592) = mat(k, 592) - dti(k)
         mat(k, 596) = mat(k, 596) - dti(k)
         mat(k, 605) = mat(k, 605) - dti(k)
         mat(k, 612) = mat(k, 612) - dti(k)
         mat(k, 619) = mat(k, 619) - dti(k)
         mat(k, 627) = mat(k, 627) - dti(k)
         mat(k, 634) = mat(k, 634) - dti(k)
         mat(k, 644) = mat(k, 644) - dti(k)
         mat(k, 657) = mat(k, 657) - dti(k)
         mat(k, 668) = mat(k, 668) - dti(k)
         mat(k, 679) = mat(k, 679) - dti(k)
         mat(k, 686) = mat(k, 686) - dti(k)
         mat(k, 695) = mat(k, 695) - dti(k)
         mat(k, 708) = mat(k, 708) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 726) = mat(k, 726) - dti(k)
         mat(k, 737) = mat(k, 737) - dti(k)
         mat(k, 750) = mat(k, 750) - dti(k)
         mat(k, 761) = mat(k, 761) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 780) = mat(k, 780) - dti(k)
         mat(k, 789) = mat(k, 789) - dti(k)
         mat(k, 799) = mat(k, 799) - dti(k)
         mat(k, 803) = mat(k, 803) - dti(k)
         mat(k, 806) = mat(k, 806) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 825) = mat(k, 825) - dti(k)
         mat(k, 836) = mat(k, 836) - dti(k)
         mat(k, 853) = mat(k, 853) - dti(k)
         mat(k, 859) = mat(k, 859) - dti(k)
         mat(k, 866) = mat(k, 866) - dti(k)
         mat(k, 875) = mat(k, 875) - dti(k)
         mat(k, 889) = mat(k, 889) - dti(k)
         mat(k, 901) = mat(k, 901) - dti(k)
         mat(k, 914) = mat(k, 914) - dti(k)
         mat(k, 924) = mat(k, 924) - dti(k)
         mat(k, 931) = mat(k, 931) - dti(k)
         mat(k, 950) = mat(k, 950) - dti(k)
         mat(k, 971) = mat(k, 971) - dti(k)
         mat(k, 981) = mat(k, 981) - dti(k)
         mat(k,1001) = mat(k,1001) - dti(k)
         mat(k,1026) = mat(k,1026) - dti(k)
         mat(k,1047) = mat(k,1047) - dti(k)
         mat(k,1061) = mat(k,1061) - dti(k)
         mat(k,1075) = mat(k,1075) - dti(k)
         mat(k,1087) = mat(k,1087) - dti(k)
         mat(k,1098) = mat(k,1098) - dti(k)
         mat(k,1111) = mat(k,1111) - dti(k)
         mat(k,1125) = mat(k,1125) - dti(k)
         mat(k,1131) = mat(k,1131) - dti(k)
         mat(k,1142) = mat(k,1142) - dti(k)
         mat(k,1155) = mat(k,1155) - dti(k)
         mat(k,1175) = mat(k,1175) - dti(k)
         mat(k,1195) = mat(k,1195) - dti(k)
         mat(k,1211) = mat(k,1211) - dti(k)
         mat(k,1223) = mat(k,1223) - dti(k)
         mat(k,1234) = mat(k,1234) - dti(k)
         mat(k,1258) = mat(k,1258) - dti(k)
         mat(k,1290) = mat(k,1290) - dti(k)
         mat(k,1313) = mat(k,1313) - dti(k)
         mat(k,1334) = mat(k,1334) - dti(k)
         mat(k,1356) = mat(k,1356) - dti(k)
         mat(k,1388) = mat(k,1388) - dti(k)
         mat(k,1403) = mat(k,1403) - dti(k)
         mat(k,1417) = mat(k,1417) - dti(k)
         mat(k,1432) = mat(k,1432) - dti(k)
         mat(k,1447) = mat(k,1447) - dti(k)
         mat(k,1465) = mat(k,1465) - dti(k)
         mat(k,1487) = mat(k,1487) - dti(k)
         mat(k,1528) = mat(k,1528) - dti(k)
         mat(k,1693) = mat(k,1693) - dti(k)
         mat(k,1720) = mat(k,1720) - dti(k)
         mat(k,1778) = mat(k,1778) - dti(k)
         mat(k,1831) = mat(k,1831) - dti(k)
         mat(k,1924) = mat(k,1924) - dti(k)
         mat(k,1969) = mat(k,1969) - dti(k)
         mat(k,2009) = mat(k,2009) - dti(k)
         mat(k,2117) = mat(k,2117) - dti(k)
         mat(k,2141) = mat(k,2141) - dti(k)
         mat(k,2166) = mat(k,2166) - dti(k)
         mat(k,2198) = mat(k,2198) - dti(k)
         mat(k,2260) = mat(k,2260) - dti(k)
         mat(k,2287) = mat(k,2287) - dti(k)
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
