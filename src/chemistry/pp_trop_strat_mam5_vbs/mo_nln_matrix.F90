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
         mat(k,550) = -(rxt(k,356)*y(k,219))
         mat(k,1513) = -rxt(k,356)*y(k,1)
         mat(k,1350) = rxt(k,359)*y(k,191)
         mat(k,833) = rxt(k,359)*y(k,124)
         mat(k,561) = -(rxt(k,360)*y(k,219))
         mat(k,1514) = -rxt(k,360)*y(k,2)
         mat(k,834) = rxt(k,357)*y(k,205)
         mat(k,1690) = rxt(k,357)*y(k,191)
         mat(k,815) = -(rxt(k,439)*y(k,126) + rxt(k,440)*y(k,135) + rxt(k,441) &
                      *y(k,219))
         mat(k,1915) = -rxt(k,439)*y(k,6)
         mat(k,1602) = -rxt(k,440)*y(k,6)
         mat(k,1537) = -rxt(k,441)*y(k,6)
         mat(k,117) = -(rxt(k,398)*y(k,219))
         mat(k,1449) = -rxt(k,398)*y(k,7)
         mat(k,300) = -(rxt(k,401)*y(k,219))
         mat(k,1479) = -rxt(k,401)*y(k,8)
         mat(k,402) = rxt(k,399)*y(k,205)
         mat(k,1665) = rxt(k,399)*y(k,193)
         mat(k,118) = .120_r8*rxt(k,398)*y(k,219)
         mat(k,1450) = .120_r8*rxt(k,398)*y(k,7)
         mat(k,812) = .100_r8*rxt(k,440)*y(k,135)
         mat(k,785) = .100_r8*rxt(k,443)*y(k,135)
         mat(k,1590) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,110)
         mat(k,1338) = .500_r8*rxt(k,400)*y(k,193) + .200_r8*rxt(k,427)*y(k,225) &
                      + .060_r8*rxt(k,433)*y(k,228)
         mat(k,403) = .500_r8*rxt(k,400)*y(k,124)
         mat(k,618) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,634) = .060_r8*rxt(k,433)*y(k,124)
         mat(k,1331) = .200_r8*rxt(k,427)*y(k,225) + .200_r8*rxt(k,433)*y(k,228)
         mat(k,617) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,632) = .200_r8*rxt(k,433)*y(k,124)
         mat(k,1347) = .200_r8*rxt(k,427)*y(k,225) + .150_r8*rxt(k,433)*y(k,228)
         mat(k,620) = .200_r8*rxt(k,427)*y(k,124)
         mat(k,635) = .150_r8*rxt(k,433)*y(k,124)
         mat(k,1332) = .210_r8*rxt(k,433)*y(k,228)
         mat(k,633) = .210_r8*rxt(k,433)*y(k,124)
         mat(k,180) = -(rxt(k,361)*y(k,219))
         mat(k,1461) = -rxt(k,361)*y(k,15)
         mat(k,811) = .050_r8*rxt(k,440)*y(k,135)
         mat(k,784) = .050_r8*rxt(k,443)*y(k,135)
         mat(k,1589) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,110)
         mat(k,286) = -(rxt(k,327)*y(k,126) + rxt(k,328)*y(k,219))
         mat(k,1908) = -rxt(k,327)*y(k,16)
         mat(k,1477) = -rxt(k,328)*y(k,16)
         mat(k,1278) = -(rxt(k,210)*y(k,42) + rxt(k,211)*y(k,205) + rxt(k,212) &
                      *y(k,135))
         mat(k,1859) = -rxt(k,210)*y(k,17)
         mat(k,1732) = -rxt(k,211)*y(k,17)
         mat(k,1626) = -rxt(k,212)*y(k,17)
         mat(k,1806) = 4.000_r8*rxt(k,213)*y(k,19) + (rxt(k,214)+rxt(k,215))*y(k,59) &
                      + rxt(k,218)*y(k,124) + rxt(k,221)*y(k,134) + rxt(k,468) &
                      *y(k,151) + rxt(k,222)*y(k,219)
         mat(k,1885) = (rxt(k,214)+rxt(k,215))*y(k,19)
         mat(k,705) = rxt(k,223)*y(k,134) + rxt(k,229)*y(k,218) + rxt(k,224)*y(k,219)
         mat(k,1388) = rxt(k,218)*y(k,19)
         mat(k,1836) = rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81)
         mat(k,1112) = rxt(k,468)*y(k,19)
         mat(k,1412) = rxt(k,229)*y(k,81)
         mat(k,1566) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,81)
         mat(k,1800) = rxt(k,216)*y(k,59)
         mat(k,1879) = rxt(k,216)*y(k,19)
         mat(k,2036) = (rxt(k,530)+rxt(k,535))*y(k,91)
         mat(k,667) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,1815) = -(4._r8*rxt(k,213)*y(k,19) + (rxt(k,214) + rxt(k,215) + rxt(k,216) &
                      ) * y(k,59) + rxt(k,217)*y(k,205) + rxt(k,218)*y(k,124) &
                      + rxt(k,219)*y(k,125) + rxt(k,221)*y(k,134) + rxt(k,222) &
                      *y(k,219) + rxt(k,468)*y(k,151))
         mat(k,1894) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,19)
         mat(k,1741) = -rxt(k,217)*y(k,19)
         mat(k,1397) = -rxt(k,218)*y(k,19)
         mat(k,1992) = -rxt(k,219)*y(k,19)
         mat(k,1845) = -rxt(k,221)*y(k,19)
         mat(k,1575) = -rxt(k,222)*y(k,19)
         mat(k,1118) = -rxt(k,468)*y(k,19)
         mat(k,1284) = rxt(k,212)*y(k,135)
         mat(k,455) = rxt(k,220)*y(k,134)
         mat(k,709) = rxt(k,230)*y(k,218)
         mat(k,671) = rxt(k,225)*y(k,134)
         mat(k,1845) = mat(k,1845) + rxt(k,220)*y(k,20) + rxt(k,225)*y(k,91)
         mat(k,1635) = rxt(k,212)*y(k,17)
         mat(k,1421) = rxt(k,230)*y(k,81)
         mat(k,451) = -(rxt(k,220)*y(k,134))
         mat(k,1826) = -rxt(k,220)*y(k,20)
         mat(k,1802) = rxt(k,219)*y(k,125)
         mat(k,1967) = rxt(k,219)*y(k,19)
         mat(k,189) = -(rxt(k,402)*y(k,219))
         mat(k,1463) = -rxt(k,402)*y(k,22)
         mat(k,1329) = rxt(k,405)*y(k,195)
         mat(k,354) = rxt(k,405)*y(k,124)
         mat(k,265) = -(rxt(k,404)*y(k,219))
         mat(k,1473) = -rxt(k,404)*y(k,23)
         mat(k,355) = rxt(k,403)*y(k,205)
         mat(k,1663) = rxt(k,403)*y(k,195)
         mat(k,221) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,219))
         mat(k,2002) = -rxt(k,276)*y(k,24)
         mat(k,1467) = -rxt(k,277)*y(k,24)
         mat(k,467) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,135) + rxt(k,304)*y(k,219))
         mat(k,2004) = -rxt(k,278)*y(k,25)
         mat(k,1593) = -rxt(k,279)*y(k,25)
         mat(k,1502) = -rxt(k,304)*y(k,25)
         mat(k,200) = -(rxt(k,284)*y(k,219))
         mat(k,1465) = -rxt(k,284)*y(k,26)
         mat(k,712) = .800_r8*rxt(k,280)*y(k,196) + .200_r8*rxt(k,281)*y(k,200)
         mat(k,1750) = .200_r8*rxt(k,281)*y(k,196)
         mat(k,273) = -(rxt(k,285)*y(k,219))
         mat(k,1475) = -rxt(k,285)*y(k,27)
         mat(k,713) = rxt(k,282)*y(k,205)
         mat(k,1664) = rxt(k,282)*y(k,196)
         mat(k,227) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,219))
         mat(k,2003) = -rxt(k,286)*y(k,28)
         mat(k,1468) = -rxt(k,287)*y(k,28)
         mat(k,877) = -(rxt(k,307)*y(k,126) + rxt(k,308)*y(k,135) + rxt(k,325) &
                      *y(k,219))
         mat(k,1919) = -rxt(k,307)*y(k,29)
         mat(k,1606) = -rxt(k,308)*y(k,29)
         mat(k,1542) = -rxt(k,325)*y(k,29)
         mat(k,749) = .130_r8*rxt(k,385)*y(k,135)
         mat(k,1606) = mat(k,1606) + .130_r8*rxt(k,385)*y(k,98)
         mat(k,342) = -(rxt(k,312)*y(k,219))
         mat(k,1485) = -rxt(k,312)*y(k,30)
         mat(k,693) = rxt(k,310)*y(k,205)
         mat(k,1670) = rxt(k,310)*y(k,197)
         mat(k,95) = -(rxt(k,313)*y(k,219))
         mat(k,1446) = -rxt(k,313)*y(k,31)
         mat(k,204) = -(rxt(k,408)*y(k,219))
         mat(k,1466) = -rxt(k,408)*y(k,32)
         mat(k,541) = rxt(k,406)*y(k,205)
         mat(k,1658) = rxt(k,406)*y(k,198)
         mat(k,1870) = -(rxt(k,174)*y(k,56) + rxt(k,210)*y(k,17) + rxt(k,254)*y(k,205) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,134) + rxt(k,257) &
                      *y(k,219))
         mat(k,2028) = -rxt(k,174)*y(k,42)
         mat(k,1286) = -rxt(k,210)*y(k,42)
         mat(k,1743) = -rxt(k,254)*y(k,42)
         mat(k,1953) = -rxt(k,255)*y(k,42)
         mat(k,1847) = -rxt(k,256)*y(k,42)
         mat(k,1577) = -rxt(k,257)*y(k,42)
         mat(k,558) = .400_r8*rxt(k,356)*y(k,219)
         mat(k,829) = .340_r8*rxt(k,440)*y(k,135)
         mat(k,291) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,474) = rxt(k,279)*y(k,135)
         mat(k,888) = .500_r8*rxt(k,308)*y(k,135)
         mat(k,436) = .500_r8*rxt(k,296)*y(k,219)
         mat(k,692) = rxt(k,262)*y(k,219)
         mat(k,372) = .300_r8*rxt(k,263)*y(k,219)
         mat(k,1896) = rxt(k,181)*y(k,200)
         mat(k,912) = .800_r8*rxt(k,301)*y(k,219)
         mat(k,761) = .910_r8*rxt(k,385)*y(k,135)
         mat(k,502) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,1083) = .800_r8*rxt(k,380)*y(k,200)
         mat(k,1097) = .120_r8*rxt(k,338)*y(k,135)
         mat(k,464) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,802) = .340_r8*rxt(k,443)*y(k,135)
         mat(k,1174) = .600_r8*rxt(k,352)*y(k,135)
         mat(k,1399) = .100_r8*rxt(k,358)*y(k,191) + rxt(k,261)*y(k,200) &
                      + .500_r8*rxt(k,329)*y(k,202) + .500_r8*rxt(k,298)*y(k,204) &
                      + .920_r8*rxt(k,368)*y(k,207) + .250_r8*rxt(k,336)*y(k,211) &
                      + rxt(k,345)*y(k,213) + rxt(k,319)*y(k,221) + rxt(k,323) &
                      *y(k,222) + .340_r8*rxt(k,452)*y(k,223) + .320_r8*rxt(k,457) &
                      *y(k,224) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1953) = mat(k,1953) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,207) &
                      + .250_r8*rxt(k,335)*y(k,211) + rxt(k,346)*y(k,213)
         mat(k,1637) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25) &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,98) &
                      + .120_r8*rxt(k,338)*y(k,105) + .340_r8*rxt(k,443)*y(k,110) &
                      + .600_r8*rxt(k,352)*y(k,111)
         mat(k,388) = rxt(k,303)*y(k,219)
         mat(k,937) = .680_r8*rxt(k,461)*y(k,219)
         mat(k,845) = .100_r8*rxt(k,358)*y(k,124)
         mat(k,721) = .700_r8*rxt(k,281)*y(k,200)
         mat(k,701) = rxt(k,309)*y(k,200)
         mat(k,1272) = rxt(k,292)*y(k,200) + rxt(k,365)*y(k,207) + .250_r8*rxt(k,332) &
                      *y(k,211) + rxt(k,341)*y(k,213) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1793) = rxt(k,181)*y(k,59) + .800_r8*rxt(k,380)*y(k,101) + rxt(k,261) &
                      *y(k,124) + .700_r8*rxt(k,281)*y(k,196) + rxt(k,309)*y(k,197) &
                      + rxt(k,292)*y(k,199) + (4.000_r8*rxt(k,258)+2.000_r8*rxt(k,259)) &
                      *y(k,200) + 1.500_r8*rxt(k,366)*y(k,207) + .750_r8*rxt(k,371) &
                      *y(k,208) + .880_r8*rxt(k,333)*y(k,211) + 2.000_r8*rxt(k,342) &
                      *y(k,213) + .750_r8*rxt(k,445)*y(k,217) + .800_r8*rxt(k,321) &
                      *y(k,222) + .930_r8*rxt(k,450)*y(k,223) + .950_r8*rxt(k,455) &
                      *y(k,224) + .800_r8*rxt(k,391)*y(k,227)
         mat(k,485) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,609) = .500_r8*rxt(k,298)*y(k,124)
         mat(k,1743) = mat(k,1743) + .450_r8*rxt(k,343)*y(k,213) + .150_r8*rxt(k,322) &
                      *y(k,222)
         mat(k,1224) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) + rxt(k,365) &
                      *y(k,199) + 1.500_r8*rxt(k,366)*y(k,200)
         mat(k,1154) = .750_r8*rxt(k,371)*y(k,200)
         mat(k,1196) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,199) + .880_r8*rxt(k,333)*y(k,200)
         mat(k,1242) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,341)*y(k,199) &
                      + 2.000_r8*rxt(k,342)*y(k,200) + .450_r8*rxt(k,343)*y(k,205) &
                      + 4.000_r8*rxt(k,344)*y(k,213)
         mat(k,1006) = .750_r8*rxt(k,445)*y(k,200)
         mat(k,1577) = mat(k,1577) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,262)*y(k,52) + .300_r8*rxt(k,263)*y(k,53) &
                      + .800_r8*rxt(k,301)*y(k,74) + .300_r8*rxt(k,376)*y(k,99) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,140) &
                      + .680_r8*rxt(k,461)*y(k,180)
         mat(k,664) = rxt(k,319)*y(k,124)
         mat(k,1020) = rxt(k,323)*y(k,124) + .800_r8*rxt(k,321)*y(k,200) &
                      + .150_r8*rxt(k,322)*y(k,205)
         mat(k,987) = .340_r8*rxt(k,452)*y(k,124) + .930_r8*rxt(k,450)*y(k,200)
         mat(k,967) = .320_r8*rxt(k,457)*y(k,124) + .950_r8*rxt(k,455)*y(k,200)
         mat(k,1060) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,390)*y(k,199) &
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
         mat(k,1024) = -(rxt(k,288)*y(k,126) + rxt(k,289)*y(k,219))
         mat(k,1930) = -rxt(k,288)*y(k,45)
         mat(k,1553) = -rxt(k,289)*y(k,45)
         mat(k,554) = .800_r8*rxt(k,356)*y(k,219)
         mat(k,289) = rxt(k,327)*y(k,126)
         mat(k,201) = rxt(k,284)*y(k,219)
         mat(k,275) = .500_r8*rxt(k,285)*y(k,219)
         mat(k,880) = .500_r8*rxt(k,308)*y(k,135)
         mat(k,1162) = .100_r8*rxt(k,352)*y(k,135)
         mat(k,1377) = .400_r8*rxt(k,358)*y(k,191) + rxt(k,283)*y(k,196) &
                      + .270_r8*rxt(k,311)*y(k,197) + rxt(k,329)*y(k,202) + rxt(k,348) &
                      *y(k,215) + rxt(k,319)*y(k,221)
         mat(k,1930) = mat(k,1930) + rxt(k,327)*y(k,16)
         mat(k,1615) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,111)
         mat(k,839) = .400_r8*rxt(k,358)*y(k,124)
         mat(k,716) = rxt(k,283)*y(k,124) + 3.200_r8*rxt(k,280)*y(k,196) &
                      + .800_r8*rxt(k,281)*y(k,200)
         mat(k,696) = .270_r8*rxt(k,311)*y(k,124)
         mat(k,1772) = .800_r8*rxt(k,281)*y(k,196)
         mat(k,482) = rxt(k,329)*y(k,124)
         mat(k,1720) = .200_r8*rxt(k,347)*y(k,215)
         mat(k,573) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,205)
         mat(k,1553) = mat(k,1553) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26) &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,660) = rxt(k,319)*y(k,124)
         mat(k,89) = -(rxt(k,290)*y(k,219))
         mat(k,1445) = -rxt(k,290)*y(k,47)
         mat(k,847) = -(rxt(k,326)*y(k,219))
         mat(k,1539) = -rxt(k,326)*y(k,48)
         mat(k,553) = .800_r8*rxt(k,356)*y(k,219)
         mat(k,817) = .520_r8*rxt(k,440)*y(k,135)
         mat(k,288) = .500_r8*rxt(k,327)*y(k,126)
         mat(k,790) = .520_r8*rxt(k,443)*y(k,135)
         mat(k,1365) = .250_r8*rxt(k,358)*y(k,191) + .820_r8*rxt(k,311)*y(k,197) &
                      + .500_r8*rxt(k,329)*y(k,202) + .270_r8*rxt(k,452)*y(k,223) &
                      + .040_r8*rxt(k,457)*y(k,224)
         mat(k,1917) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,1604) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,110)
         mat(k,929) = .500_r8*rxt(k,461)*y(k,219)
         mat(k,838) = .250_r8*rxt(k,358)*y(k,124)
         mat(k,695) = .820_r8*rxt(k,311)*y(k,124) + .820_r8*rxt(k,309)*y(k,200)
         mat(k,1761) = .820_r8*rxt(k,309)*y(k,197) + .150_r8*rxt(k,450)*y(k,223) &
                      + .025_r8*rxt(k,455)*y(k,224)
         mat(k,480) = .500_r8*rxt(k,329)*y(k,124)
         mat(k,1539) = mat(k,1539) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,180)
         mat(k,974) = .270_r8*rxt(k,452)*y(k,124) + .150_r8*rxt(k,450)*y(k,200)
         mat(k,952) = .040_r8*rxt(k,457)*y(k,124) + .025_r8*rxt(k,455)*y(k,200)
         mat(k,1100) = -(rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219))
         mat(k,1934) = -rxt(k,314)*y(k,49)
         mat(k,1558) = -rxt(k,315)*y(k,49)
         mat(k,944) = rxt(k,316)*y(k,219)
         mat(k,1089) = .880_r8*rxt(k,338)*y(k,135)
         mat(k,1163) = .500_r8*rxt(k,352)*y(k,135)
         mat(k,1381) = .170_r8*rxt(k,411)*y(k,201) + .050_r8*rxt(k,374)*y(k,208) &
                      + .250_r8*rxt(k,336)*y(k,211) + .170_r8*rxt(k,417)*y(k,214) &
                      + .400_r8*rxt(k,427)*y(k,225) + .250_r8*rxt(k,393)*y(k,227) &
                      + .540_r8*rxt(k,433)*y(k,228) + .510_r8*rxt(k,436)*y(k,230)
         mat(k,1934) = mat(k,1934) + .050_r8*rxt(k,375)*y(k,208) + .250_r8*rxt(k,335) &
                      *y(k,211) + .250_r8*rxt(k,394)*y(k,227)
         mat(k,739) = rxt(k,317)*y(k,219)
         mat(k,1618) = .880_r8*rxt(k,338)*y(k,105) + .500_r8*rxt(k,352)*y(k,111)
         mat(k,1258) = .250_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1776) = .240_r8*rxt(k,333)*y(k,211) + .500_r8*rxt(k,321)*y(k,222) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,651) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,205)
         mat(k,1725) = .070_r8*rxt(k,410)*y(k,201) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1141) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1186) = .250_r8*rxt(k,336)*y(k,124) + .250_r8*rxt(k,335)*y(k,126) &
                      + .250_r8*rxt(k,332)*y(k,199) + .240_r8*rxt(k,333)*y(k,200)
         mat(k,772) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1558) = mat(k,1558) + rxt(k,316)*y(k,95) + rxt(k,317)*y(k,127)
         mat(k,1014) = .500_r8*rxt(k,321)*y(k,200)
         mat(k,627) = .400_r8*rxt(k,427)*y(k,124)
         mat(k,1053) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,643) = .540_r8*rxt(k,433)*y(k,124)
         mat(k,414) = .510_r8*rxt(k,436)*y(k,124)
         mat(k,475) = -(rxt(k,295)*y(k,219))
         mat(k,1503) = -rxt(k,295)*y(k,50)
         mat(k,873) = .120_r8*rxt(k,308)*y(k,135)
         mat(k,1594) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1249) = .100_r8*rxt(k,292)*y(k,200) + .150_r8*rxt(k,293)*y(k,205)
         mat(k,1754) = .100_r8*rxt(k,292)*y(k,199)
         mat(k,1684) = .150_r8*rxt(k,293)*y(k,199) + .150_r8*rxt(k,343)*y(k,213)
         mat(k,1230) = .150_r8*rxt(k,343)*y(k,205)
         mat(k,432) = -(rxt(k,296)*y(k,219))
         mat(k,1498) = -rxt(k,296)*y(k,51)
         mat(k,1248) = .400_r8*rxt(k,293)*y(k,205)
         mat(k,1682) = .400_r8*rxt(k,293)*y(k,199) + .400_r8*rxt(k,343)*y(k,213)
         mat(k,1228) = .400_r8*rxt(k,343)*y(k,205)
         mat(k,689) = -(rxt(k,262)*y(k,219))
         mat(k,1526) = -rxt(k,262)*y(k,52)
         mat(k,1066) = .200_r8*rxt(k,380)*y(k,200)
         mat(k,714) = .300_r8*rxt(k,281)*y(k,200)
         mat(k,1756) = .200_r8*rxt(k,380)*y(k,101) + .300_r8*rxt(k,281)*y(k,196) &
                      + 2.000_r8*rxt(k,259)*y(k,200) + .250_r8*rxt(k,366)*y(k,207) &
                      + .250_r8*rxt(k,371)*y(k,208) + .250_r8*rxt(k,333)*y(k,211) &
                      + .250_r8*rxt(k,445)*y(k,217) + .500_r8*rxt(k,321)*y(k,222) &
                      + .250_r8*rxt(k,450)*y(k,223) + .250_r8*rxt(k,455)*y(k,224) &
                      + .300_r8*rxt(k,391)*y(k,227)
         mat(k,1202) = .250_r8*rxt(k,366)*y(k,200)
         mat(k,1129) = .250_r8*rxt(k,371)*y(k,200)
         mat(k,1179) = .250_r8*rxt(k,333)*y(k,200)
         mat(k,992) = .250_r8*rxt(k,445)*y(k,200)
         mat(k,1011) = .500_r8*rxt(k,321)*y(k,200)
         mat(k,973) = .250_r8*rxt(k,450)*y(k,200)
         mat(k,951) = .250_r8*rxt(k,455)*y(k,200)
         mat(k,1047) = .300_r8*rxt(k,391)*y(k,200)
         mat(k,368) = -(rxt(k,263)*y(k,219))
         mat(k,1488) = -rxt(k,263)*y(k,53)
         mat(k,1753) = rxt(k,260)*y(k,205)
         mat(k,1674) = rxt(k,260)*y(k,200)
         mat(k,2032) = -(rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) + rxt(k,177)*y(k,79) &
                      + (rxt(k,178) + rxt(k,179)) * y(k,205) + rxt(k,180)*y(k,135) &
                      + rxt(k,187)*y(k,60) + rxt(k,196)*y(k,92) + rxt(k,286)*y(k,28))
         mat(k,1874) = -rxt(k,174)*y(k,56)
         mat(k,1043) = -rxt(k,176)*y(k,56)
         mat(k,528) = -rxt(k,177)*y(k,56)
         mat(k,1747) = -(rxt(k,178) + rxt(k,179)) * y(k,56)
         mat(k,1641) = -rxt(k,180)*y(k,56)
         mat(k,862) = -rxt(k,187)*y(k,56)
         mat(k,728) = -rxt(k,196)*y(k,56)
         mat(k,230) = -rxt(k,286)*y(k,56)
         mat(k,1821) = rxt(k,215)*y(k,59)
         mat(k,1900) = rxt(k,215)*y(k,19) + (4.000_r8*rxt(k,182)+2.000_r8*rxt(k,184)) &
                      *y(k,59) + rxt(k,186)*y(k,124) + rxt(k,191)*y(k,134) &
                      + rxt(k,469)*y(k,151) + rxt(k,181)*y(k,200) + rxt(k,192) &
                      *y(k,219)
         mat(k,159) = rxt(k,236)*y(k,218)
         mat(k,2054) = rxt(k,194)*y(k,134) + rxt(k,206)*y(k,218) + rxt(k,195)*y(k,219)
         mat(k,1403) = rxt(k,186)*y(k,59)
         mat(k,1851) = rxt(k,191)*y(k,59) + rxt(k,194)*y(k,85)
         mat(k,1122) = rxt(k,469)*y(k,59)
         mat(k,1797) = rxt(k,181)*y(k,59)
         mat(k,1427) = rxt(k,236)*y(k,65) + rxt(k,206)*y(k,85)
         mat(k,1581) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,85)
         mat(k,2001) = rxt(k,187)*y(k,60)
         mat(k,1878) = 2.000_r8*rxt(k,183)*y(k,59)
         mat(k,853) = rxt(k,187)*y(k,56) + (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,2035) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60) + (rxt(k,523) &
                       +rxt(k,529)+rxt(k,534))*y(k,92)
         mat(k,723) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85)
         mat(k,1877) = 2.000_r8*rxt(k,208)*y(k,59)
         mat(k,1897) = -(rxt(k,181)*y(k,200) + (4._r8*rxt(k,182) + 4._r8*rxt(k,183) &
                      + 4._r8*rxt(k,184) + 4._r8*rxt(k,208)) * y(k,59) + rxt(k,185) &
                      *y(k,205) + rxt(k,186)*y(k,124) + rxt(k,188)*y(k,125) + rxt(k,191) &
                      *y(k,134) + (rxt(k,192) + rxt(k,193)) * y(k,219) + (rxt(k,214) &
                      + rxt(k,215) + rxt(k,216)) * y(k,19) + rxt(k,469)*y(k,151))
         mat(k,1794) = -rxt(k,181)*y(k,59)
         mat(k,1744) = -rxt(k,185)*y(k,59)
         mat(k,1400) = -rxt(k,186)*y(k,59)
         mat(k,1995) = -rxt(k,188)*y(k,59)
         mat(k,1848) = -rxt(k,191)*y(k,59)
         mat(k,1578) = -(rxt(k,192) + rxt(k,193)) * y(k,59)
         mat(k,1818) = -(rxt(k,214) + rxt(k,215) + rxt(k,216)) * y(k,59)
         mat(k,1120) = -rxt(k,469)*y(k,59)
         mat(k,2029) = rxt(k,196)*y(k,92) + rxt(k,180)*y(k,135) + rxt(k,179)*y(k,205)
         mat(k,859) = rxt(k,189)*y(k,134)
         mat(k,2051) = rxt(k,207)*y(k,218)
         mat(k,727) = rxt(k,196)*y(k,56) + rxt(k,197)*y(k,134) + rxt(k,198)*y(k,219)
         mat(k,1848) = mat(k,1848) + rxt(k,189)*y(k,60) + rxt(k,197)*y(k,92)
         mat(k,1638) = rxt(k,180)*y(k,56)
         mat(k,258) = rxt(k,474)*y(k,151)
         mat(k,1120) = mat(k,1120) + rxt(k,474)*y(k,137)
         mat(k,1744) = mat(k,1744) + rxt(k,179)*y(k,56)
         mat(k,1424) = rxt(k,207)*y(k,85)
         mat(k,1578) = mat(k,1578) + rxt(k,198)*y(k,92)
         mat(k,855) = -(rxt(k,187)*y(k,56) + rxt(k,189)*y(k,134) + rxt(k,190)*y(k,219) &
                      + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85))
         mat(k,2011) = -rxt(k,187)*y(k,60)
         mat(k,1832) = -rxt(k,189)*y(k,60)
         mat(k,1540) = -rxt(k,190)*y(k,60)
         mat(k,2039) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60)
         mat(k,1883) = rxt(k,188)*y(k,125)
         mat(k,1976) = rxt(k,188)*y(k,59)
         mat(k,939) = -((rxt(k,265) + rxt(k,275)) * y(k,219))
         mat(k,1547) = -(rxt(k,265) + rxt(k,275)) * y(k,62)
         mat(k,820) = .230_r8*rxt(k,440)*y(k,135)
         mat(k,1277) = rxt(k,210)*y(k,42)
         mat(k,224) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,470) = .630_r8*rxt(k,279)*y(k,135)
         mat(k,878) = .560_r8*rxt(k,308)*y(k,135)
         mat(k,1857) = rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) + rxt(k,255)*y(k,126) &
                      + rxt(k,256)*y(k,134) + rxt(k,257)*y(k,219)
         mat(k,1099) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219)
         mat(k,2013) = rxt(k,174)*y(k,42)
         mat(k,766) = rxt(k,302)*y(k,219)
         mat(k,750) = .620_r8*rxt(k,385)*y(k,135)
         mat(k,1087) = .650_r8*rxt(k,338)*y(k,135)
         mat(k,793) = .230_r8*rxt(k,443)*y(k,135)
         mat(k,1160) = .560_r8*rxt(k,352)*y(k,135)
         mat(k,1371) = .170_r8*rxt(k,411)*y(k,201) + .220_r8*rxt(k,336)*y(k,211) &
                      + .400_r8*rxt(k,414)*y(k,212) + .350_r8*rxt(k,417)*y(k,214) &
                      + .225_r8*rxt(k,452)*y(k,223) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1924) = rxt(k,255)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,211) + .500_r8*rxt(k,394)*y(k,227)
         mat(k,1833) = rxt(k,256)*y(k,42) + rxt(k,464)*y(k,138)
         mat(k,1609) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25) &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,98) &
                      + .650_r8*rxt(k,338)*y(k,105) + .230_r8*rxt(k,443)*y(k,110) &
                      + .560_r8*rxt(k,352)*y(k,111)
         mat(k,281) = rxt(k,464)*y(k,134) + rxt(k,465)*y(k,219)
         mat(k,931) = .700_r8*rxt(k,461)*y(k,219)
         mat(k,1253) = .220_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1766) = .110_r8*rxt(k,333)*y(k,211) + .125_r8*rxt(k,450)*y(k,223) &
                      + .200_r8*rxt(k,391)*y(k,227)
         mat(k,650) = .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410)*y(k,205)
         mat(k,1714) = .070_r8*rxt(k,410)*y(k,201) + .160_r8*rxt(k,413)*y(k,212) &
                      + .140_r8*rxt(k,416)*y(k,214)
         mat(k,1182) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,199) + .110_r8*rxt(k,333)*y(k,200)
         mat(k,613) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,205)
         mat(k,771) = .350_r8*rxt(k,417)*y(k,124) + .140_r8*rxt(k,416)*y(k,205)
         mat(k,1547) = mat(k,1547) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,257)*y(k,42) &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,75) + rxt(k,465)*y(k,138) &
                      + .700_r8*rxt(k,461)*y(k,180)
         mat(k,977) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,200)
         mat(k,1050) = .250_r8*rxt(k,393)*y(k,124) + .500_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .200_r8*rxt(k,391)*y(k,200)
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
         mat(k,111) = -(rxt(k,235)*y(k,218))
         mat(k,1406) = -rxt(k,235)*y(k,64)
         mat(k,156) = -(rxt(k,236)*y(k,218))
         mat(k,1408) = -rxt(k,236)*y(k,65)
         mat(k,144) = -(rxt(k,409)*y(k,219))
         mat(k,1454) = -rxt(k,409)*y(k,66)
         mat(k,138) = .180_r8*rxt(k,429)*y(k,219)
         mat(k,1454) = mat(k,1454) + .180_r8*rxt(k,429)*y(k,182)
         mat(k,233) = -(rxt(k,462)*y(k,126) + (rxt(k,463) + rxt(k,476)) * y(k,219))
         mat(k,1906) = -rxt(k,462)*y(k,67)
         mat(k,1469) = -(rxt(k,463) + rxt(k,476)) * y(k,67)
         mat(k,602) = rxt(k,297)*y(k,205)
         mat(k,1656) = rxt(k,297)*y(k,204)
         mat(k,677) = -(rxt(k,232)*y(k,77) + rxt(k,233)*y(k,231) + rxt(k,234)*y(k,89))
         mat(k,1034) = -rxt(k,232)*y(k,73)
         mat(k,2060) = -rxt(k,233)*y(k,73)
         mat(k,1289) = -rxt(k,234)*y(k,73)
         mat(k,112) = 2.000_r8*rxt(k,235)*y(k,218)
         mat(k,157) = rxt(k,236)*y(k,218)
         mat(k,1409) = 2.000_r8*rxt(k,235)*y(k,64) + rxt(k,236)*y(k,65)
         mat(k,908) = -(rxt(k,301)*y(k,219))
         mat(k,1544) = -rxt(k,301)*y(k,74)
         mat(k,496) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,489) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,296) = rxt(k,388)*y(k,219)
         mat(k,1368) = .050_r8*rxt(k,374)*y(k,208) + .530_r8*rxt(k,336)*y(k,211) &
                      + .225_r8*rxt(k,452)*y(k,223) + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1921) = .050_r8*rxt(k,375)*y(k,208) + .530_r8*rxt(k,335)*y(k,211) &
                      + .250_r8*rxt(k,394)*y(k,227)
         mat(k,1252) = .530_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1764) = .260_r8*rxt(k,333)*y(k,211) + .125_r8*rxt(k,450)*y(k,223) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,1133) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1180) = .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335)*y(k,126) &
                      + .530_r8*rxt(k,332)*y(k,199) + .260_r8*rxt(k,333)*y(k,200)
         mat(k,1544) = mat(k,1544) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + rxt(k,388)*y(k,115)
         mat(k,975) = .225_r8*rxt(k,452)*y(k,124) + .125_r8*rxt(k,450)*y(k,200)
         mat(k,1049) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,765) = -(rxt(k,302)*y(k,219))
         mat(k,1534) = -rxt(k,302)*y(k,75)
         mat(k,223) = .650_r8*rxt(k,277)*y(k,219)
         mat(k,907) = .200_r8*rxt(k,301)*y(k,219)
         mat(k,894) = rxt(k,389)*y(k,219)
         mat(k,1362) = rxt(k,400)*y(k,193) + .050_r8*rxt(k,374)*y(k,208) &
                      + .400_r8*rxt(k,414)*y(k,212) + .170_r8*rxt(k,417)*y(k,214) &
                      + .700_r8*rxt(k,420)*y(k,220) + .600_r8*rxt(k,427)*y(k,225) &
                      + .250_r8*rxt(k,393)*y(k,227) + .340_r8*rxt(k,433)*y(k,228) &
                      + .170_r8*rxt(k,436)*y(k,230)
         mat(k,1913) = .050_r8*rxt(k,375)*y(k,208) + .250_r8*rxt(k,394)*y(k,227)
         mat(k,406) = rxt(k,400)*y(k,124)
         mat(k,1250) = .250_r8*rxt(k,390)*y(k,227)
         mat(k,1760) = .100_r8*rxt(k,391)*y(k,227)
         mat(k,1707) = .160_r8*rxt(k,413)*y(k,212) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1131) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,612) = .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413)*y(k,205)
         mat(k,769) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1534) = mat(k,1534) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,74) + rxt(k,389)*y(k,116)
         mat(k,376) = .700_r8*rxt(k,420)*y(k,124)
         mat(k,624) = .600_r8*rxt(k,427)*y(k,124)
         mat(k,1048) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
         mat(k,640) = .340_r8*rxt(k,433)*y(k,124)
         mat(k,413) = .170_r8*rxt(k,436)*y(k,124)
         mat(k,1304) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,205) + rxt(k,140) &
                      *y(k,135))
         mat(k,1734) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76)
         mat(k,1628) = -rxt(k,140)*y(k,76)
         mat(k,1861) = rxt(k,257)*y(k,219)
         mat(k,2019) = rxt(k,176)*y(k,77)
         mat(k,940) = rxt(k,275)*y(k,219)
         mat(k,680) = rxt(k,232)*y(k,77)
         mat(k,1037) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73) + rxt(k,132)*y(k,134) &
                      + rxt(k,124)*y(k,218) + rxt(k,141)*y(k,219)
         mat(k,706) = rxt(k,230)*y(k,218)
         mat(k,2042) = rxt(k,207)*y(k,218)
         mat(k,325) = rxt(k,162)*y(k,219)
         mat(k,1838) = rxt(k,132)*y(k,77) + rxt(k,144)*y(k,219)
         mat(k,283) = rxt(k,465)*y(k,219)
         mat(k,421) = rxt(k,470)*y(k,219)
         mat(k,1113) = rxt(k,475)*y(k,219)
         mat(k,1414) = rxt(k,124)*y(k,77) + rxt(k,230)*y(k,81) + rxt(k,207)*y(k,85)
         mat(k,1568) = rxt(k,257)*y(k,42) + rxt(k,275)*y(k,62) + rxt(k,141)*y(k,77) &
                      + rxt(k,162)*y(k,112) + rxt(k,144)*y(k,134) + rxt(k,465) &
                      *y(k,138) + rxt(k,470)*y(k,149) + rxt(k,475)*y(k,151)
         mat(k,1035) = -(rxt(k,124)*y(k,218) + rxt(k,132)*y(k,134) + rxt(k,141) &
                      *y(k,219) + rxt(k,176)*y(k,56) + rxt(k,232)*y(k,73))
         mat(k,1411) = -rxt(k,124)*y(k,77)
         mat(k,1834) = -rxt(k,132)*y(k,77)
         mat(k,1554) = -rxt(k,141)*y(k,77)
         mat(k,2015) = -rxt(k,176)*y(k,77)
         mat(k,678) = -rxt(k,232)*y(k,77)
         mat(k,1302) = rxt(k,134)*y(k,205)
         mat(k,1721) = rxt(k,134)*y(k,76)
         mat(k,524) = -(rxt(k,133)*y(k,134) + rxt(k,142)*y(k,219) + rxt(k,177)*y(k,56))
         mat(k,1827) = -rxt(k,133)*y(k,79)
         mat(k,1509) = -rxt(k,142)*y(k,79)
         mat(k,2005) = -rxt(k,177)*y(k,79)
         mat(k,1687) = 2.000_r8*rxt(k,148)*y(k,205)
         mat(k,1509) = mat(k,1509) + 2.000_r8*rxt(k,147)*y(k,219)
         mat(k,195) = rxt(k,478)*y(k,231)
         mat(k,2057) = rxt(k,478)*y(k,153)
         mat(k,704) = -(rxt(k,223)*y(k,134) + rxt(k,224)*y(k,219) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,218))
         mat(k,1829) = -rxt(k,223)*y(k,81)
         mat(k,1528) = -rxt(k,224)*y(k,81)
         mat(k,1410) = -(rxt(k,229) + rxt(k,230)) * y(k,81)
         mat(k,1276) = rxt(k,210)*y(k,42) + rxt(k,211)*y(k,205)
         mat(k,1856) = rxt(k,210)*y(k,17)
         mat(k,1703) = rxt(k,211)*y(k,17)
         mat(k,2055) = -(rxt(k,194)*y(k,134) + rxt(k,195)*y(k,219) + (rxt(k,206) &
                      + rxt(k,207)) * y(k,218) + (rxt(k,523) + rxt(k,529) + rxt(k,534) &
                      ) * y(k,92) + (rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,60) &
                      + (rxt(k,530) + rxt(k,535)) * y(k,91))
         mat(k,1852) = -rxt(k,194)*y(k,85)
         mat(k,1582) = -rxt(k,195)*y(k,85)
         mat(k,1428) = -(rxt(k,206) + rxt(k,207)) * y(k,85)
         mat(k,729) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85)
         mat(k,863) = -(rxt(k,528) + rxt(k,533) + rxt(k,538)) * y(k,85)
         mat(k,674) = -(rxt(k,530) + rxt(k,535)) * y(k,85)
         mat(k,231) = rxt(k,286)*y(k,56)
         mat(k,1875) = rxt(k,174)*y(k,56)
         mat(k,2033) = rxt(k,286)*y(k,28) + rxt(k,174)*y(k,42) + rxt(k,176)*y(k,77) &
                      + rxt(k,177)*y(k,79) + rxt(k,196)*y(k,92) + rxt(k,178)*y(k,205)
         mat(k,1901) = rxt(k,193)*y(k,219)
         mat(k,1044) = rxt(k,176)*y(k,56)
         mat(k,529) = rxt(k,177)*y(k,56)
         mat(k,729) = mat(k,729) + rxt(k,196)*y(k,56)
         mat(k,1748) = rxt(k,178)*y(k,56)
         mat(k,1582) = mat(k,1582) + rxt(k,193)*y(k,59)
         mat(k,132) = -(rxt(k,266)*y(k,219) + rxt(k,274)*y(k,218))
         mat(k,1452) = -rxt(k,266)*y(k,86)
         mat(k,1407) = -rxt(k,274)*y(k,86)
         mat(k,685) = -(rxt(k,267)*y(k,219))
         mat(k,1525) = -rxt(k,267)*y(k,87)
         mat(k,813) = .050_r8*rxt(k,440)*y(k,135)
         mat(k,222) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,469) = .370_r8*rxt(k,279)*y(k,135)
         mat(k,875) = .120_r8*rxt(k,308)*y(k,135)
         mat(k,747) = .110_r8*rxt(k,385)*y(k,135)
         mat(k,1086) = .330_r8*rxt(k,338)*y(k,135)
         mat(k,786) = .050_r8*rxt(k,443)*y(k,135)
         mat(k,1158) = .120_r8*rxt(k,352)*y(k,135)
         mat(k,1358) = rxt(k,270)*y(k,206)
         mat(k,1597) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25) &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,98) &
                      + .330_r8*rxt(k,338)*y(k,105) + .050_r8*rxt(k,443)*y(k,110) &
                      + .120_r8*rxt(k,352)*y(k,111)
         mat(k,1701) = rxt(k,268)*y(k,206)
         mat(k,363) = rxt(k,270)*y(k,124) + rxt(k,268)*y(k,205)
         mat(k,1525) = mat(k,1525) + .350_r8*rxt(k,277)*y(k,24)
         mat(k,676) = rxt(k,232)*y(k,77) + rxt(k,234)*y(k,89) + rxt(k,233)*y(k,231)
         mat(k,1033) = rxt(k,232)*y(k,73)
         mat(k,1288) = rxt(k,234)*y(k,73)
         mat(k,2058) = rxt(k,233)*y(k,73)
         mat(k,1291) = -(rxt(k,171)*y(k,219) + rxt(k,234)*y(k,73))
         mat(k,1567) = -rxt(k,171)*y(k,89)
         mat(k,679) = -rxt(k,234)*y(k,89)
         mat(k,1860) = rxt(k,255)*y(k,126)
         mat(k,1026) = rxt(k,288)*y(k,126)
         mat(k,1102) = rxt(k,314)*y(k,126)
         mat(k,856) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,85)
         mat(k,235) = rxt(k,462)*y(k,126)
         mat(k,2041) = (rxt(k,528)+rxt(k,533)+rxt(k,538))*y(k,60)
         mat(k,1984) = rxt(k,170)*y(k,219)
         mat(k,1943) = rxt(k,255)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) &
                      + rxt(k,462)*y(k,67)
         mat(k,1567) = mat(k,1567) + rxt(k,170)*y(k,125)
         mat(k,318) = -(rxt(k,149)*y(k,219))
         mat(k,1481) = -rxt(k,149)*y(k,90)
         mat(k,1962) = rxt(k,168)*y(k,205)
         mat(k,1667) = rxt(k,168)*y(k,125)
         mat(k,668) = -(rxt(k,225)*y(k,134) + (rxt(k,530) + rxt(k,535)) * y(k,85))
         mat(k,1828) = -rxt(k,225)*y(k,91)
         mat(k,2037) = -(rxt(k,530) + rxt(k,535)) * y(k,91)
         mat(k,1803) = rxt(k,217)*y(k,205)
         mat(k,1700) = rxt(k,217)*y(k,19)
         mat(k,724) = -(rxt(k,196)*y(k,56) + rxt(k,197)*y(k,134) + rxt(k,198)*y(k,219) &
                      + (rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,85))
         mat(k,2009) = -rxt(k,196)*y(k,92)
         mat(k,1830) = -rxt(k,197)*y(k,92)
         mat(k,1530) = -rxt(k,198)*y(k,92)
         mat(k,2038) = -(rxt(k,523) + rxt(k,529) + rxt(k,534)) * y(k,92)
         mat(k,1881) = rxt(k,185)*y(k,205)
         mat(k,854) = rxt(k,190)*y(k,219)
         mat(k,1705) = rxt(k,185)*y(k,59)
         mat(k,1530) = mat(k,1530) + rxt(k,190)*y(k,60)
         mat(k,916) = -(rxt(k,331)*y(k,219))
         mat(k,1545) = -rxt(k,331)*y(k,93)
         mat(k,497) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,490) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,1369) = rxt(k,330)*y(k,202) + rxt(k,337)*y(k,211)
         mat(k,481) = rxt(k,330)*y(k,124)
         mat(k,1181) = rxt(k,337)*y(k,124)
         mat(k,1545) = mat(k,1545) + .300_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100)
         mat(k,175) = -(rxt(k,362)*y(k,219))
         mat(k,1460) = -rxt(k,362)*y(k,94)
         mat(k,943) = -(rxt(k,316)*y(k,219))
         mat(k,1548) = -rxt(k,316)*y(k,95)
         mat(k,498) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,491) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,460) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,1372) = .050_r8*rxt(k,374)*y(k,208) + .220_r8*rxt(k,336)*y(k,211) &
                      + .250_r8*rxt(k,393)*y(k,227)
         mat(k,1925) = .050_r8*rxt(k,375)*y(k,208) + .220_r8*rxt(k,335)*y(k,211) &
                      + .250_r8*rxt(k,394)*y(k,227)
         mat(k,445) = .500_r8*rxt(k,320)*y(k,219)
         mat(k,1254) = .220_r8*rxt(k,332)*y(k,211) + .250_r8*rxt(k,390)*y(k,227)
         mat(k,1767) = .230_r8*rxt(k,333)*y(k,211) + .200_r8*rxt(k,321)*y(k,222) &
                      + .100_r8*rxt(k,391)*y(k,227)
         mat(k,1136) = .050_r8*rxt(k,374)*y(k,124) + .050_r8*rxt(k,375)*y(k,126)
         mat(k,1183) = .220_r8*rxt(k,336)*y(k,124) + .220_r8*rxt(k,335)*y(k,126) &
                      + .220_r8*rxt(k,332)*y(k,199) + .230_r8*rxt(k,333)*y(k,200)
         mat(k,1548) = mat(k,1548) + .700_r8*rxt(k,376)*y(k,99) + .500_r8*rxt(k,377) &
                      *y(k,100) + .500_r8*rxt(k,351)*y(k,109) + .500_r8*rxt(k,320) &
                      *y(k,147)
         mat(k,1012) = .200_r8*rxt(k,321)*y(k,200)
         mat(k,1051) = .250_r8*rxt(k,393)*y(k,124) + .250_r8*rxt(k,394)*y(k,126) &
                      + .250_r8*rxt(k,390)*y(k,199) + .100_r8*rxt(k,391)*y(k,200)
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
         mat(k,270) = -(rxt(k,363)*y(k,219))
         mat(k,1474) = -rxt(k,363)*y(k,96)
         mat(k,1333) = .870_r8*rxt(k,374)*y(k,208)
         mat(k,1907) = .950_r8*rxt(k,375)*y(k,208)
         mat(k,1246) = rxt(k,370)*y(k,208)
         mat(k,1751) = .750_r8*rxt(k,371)*y(k,208)
         mat(k,1125) = .870_r8*rxt(k,374)*y(k,124) + .950_r8*rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,199) + .750_r8*rxt(k,371)*y(k,200)
         mat(k,105) = -(rxt(k,364)*y(k,219))
         mat(k,1447) = -rxt(k,364)*y(k,97)
         mat(k,579) = .600_r8*rxt(k,387)*y(k,219)
         mat(k,1447) = mat(k,1447) + .600_r8*rxt(k,387)*y(k,103)
         mat(k,748) = -(rxt(k,378)*y(k,126) + rxt(k,385)*y(k,135) + rxt(k,386) &
                      *y(k,219))
         mat(k,1912) = -rxt(k,378)*y(k,98)
         mat(k,1599) = -rxt(k,385)*y(k,98)
         mat(k,1533) = -rxt(k,386)*y(k,98)
         mat(k,495) = -(rxt(k,376)*y(k,219))
         mat(k,1506) = -rxt(k,376)*y(k,99)
         mat(k,1346) = .080_r8*rxt(k,368)*y(k,207)
         mat(k,1200) = .080_r8*rxt(k,368)*y(k,124)
         mat(k,487) = -(rxt(k,377)*y(k,219))
         mat(k,1505) = -rxt(k,377)*y(k,100)
         mat(k,1345) = .080_r8*rxt(k,374)*y(k,208)
         mat(k,1126) = .080_r8*rxt(k,374)*y(k,124)
         mat(k,1072) = -(rxt(k,379)*y(k,199) + rxt(k,380)*y(k,200) + rxt(k,381) &
                      *y(k,205) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126))
         mat(k,1256) = -rxt(k,379)*y(k,101)
         mat(k,1774) = -rxt(k,380)*y(k,101)
         mat(k,1723) = -rxt(k,381)*y(k,101)
         mat(k,1379) = -rxt(k,382)*y(k,101)
         mat(k,1932) = -rxt(k,383)*y(k,101)
         mat(k,751) = rxt(k,378)*y(k,126)
         mat(k,1932) = mat(k,1932) + rxt(k,378)*y(k,98)
         mat(k,336) = -(rxt(k,384)*y(k,219))
         mat(k,1484) = -rxt(k,384)*y(k,102)
         mat(k,1064) = rxt(k,381)*y(k,205)
         mat(k,1669) = rxt(k,381)*y(k,101)
         mat(k,580) = -(rxt(k,387)*y(k,219))
         mat(k,1516) = -rxt(k,387)*y(k,103)
         mat(k,1692) = rxt(k,367)*y(k,207) + rxt(k,372)*y(k,208)
         mat(k,1201) = rxt(k,367)*y(k,205)
         mat(k,1128) = rxt(k,372)*y(k,205)
         mat(k,67) = -(rxt(k,509)*y(k,219))
         mat(k,1440) = -rxt(k,509)*y(k,104)
         mat(k,1088) = -(rxt(k,338)*y(k,135) + rxt(k,339)*y(k,219))
         mat(k,1617) = -rxt(k,338)*y(k,105)
         mat(k,1557) = -rxt(k,339)*y(k,105)
         mat(k,752) = .300_r8*rxt(k,385)*y(k,135)
         mat(k,1380) = .360_r8*rxt(k,368)*y(k,207)
         mat(k,1933) = .400_r8*rxt(k,369)*y(k,207)
         mat(k,1617) = mat(k,1617) + .300_r8*rxt(k,385)*y(k,98)
         mat(k,1257) = .390_r8*rxt(k,365)*y(k,207)
         mat(k,1775) = .310_r8*rxt(k,366)*y(k,207)
         mat(k,1210) = .360_r8*rxt(k,368)*y(k,124) + .400_r8*rxt(k,369)*y(k,126) &
                      + .390_r8*rxt(k,365)*y(k,199) + .310_r8*rxt(k,366)*y(k,200)
         mat(k,239) = -(rxt(k,340)*y(k,219))
         mat(k,1470) = -rxt(k,340)*y(k,106)
         mat(k,1660) = rxt(k,334)*y(k,211)
         mat(k,1178) = rxt(k,334)*y(k,205)
         mat(k,438) = -(rxt(k,349)*y(k,219))
         mat(k,1499) = -rxt(k,349)*y(k,107)
         mat(k,1343) = .800_r8*rxt(k,358)*y(k,191)
         mat(k,832) = .800_r8*rxt(k,358)*y(k,124)
         mat(k,244) = -(rxt(k,350)*y(k,219))
         mat(k,1471) = -rxt(k,350)*y(k,108)
         mat(k,1661) = .800_r8*rxt(k,347)*y(k,215)
         mat(k,571) = .800_r8*rxt(k,347)*y(k,205)
         mat(k,459) = -(rxt(k,351)*y(k,219))
         mat(k,1501) = -rxt(k,351)*y(k,109)
         mat(k,1968) = rxt(k,354)*y(k,213)
         mat(k,1229) = rxt(k,354)*y(k,125)
         mat(k,788) = -(rxt(k,442)*y(k,126) + rxt(k,443)*y(k,135) + rxt(k,444) &
                      *y(k,219))
         mat(k,1914) = -rxt(k,442)*y(k,110)
         mat(k,1601) = -rxt(k,443)*y(k,110)
         mat(k,1536) = -rxt(k,444)*y(k,110)
         mat(k,1164) = -(rxt(k,352)*y(k,135) + rxt(k,353)*y(k,219))
         mat(k,1621) = -rxt(k,352)*y(k,111)
         mat(k,1561) = -rxt(k,353)*y(k,111)
         mat(k,754) = .200_r8*rxt(k,385)*y(k,135)
         mat(k,1383) = .560_r8*rxt(k,368)*y(k,207)
         mat(k,1937) = .600_r8*rxt(k,369)*y(k,207)
         mat(k,1621) = mat(k,1621) + .200_r8*rxt(k,385)*y(k,98)
         mat(k,1260) = .610_r8*rxt(k,365)*y(k,207)
         mat(k,1778) = .440_r8*rxt(k,366)*y(k,207)
         mat(k,1212) = .560_r8*rxt(k,368)*y(k,124) + .600_r8*rxt(k,369)*y(k,126) &
                      + .610_r8*rxt(k,365)*y(k,199) + .440_r8*rxt(k,366)*y(k,200)
         mat(k,324) = -(rxt(k,150)*y(k,124) + (rxt(k,151) + rxt(k,152) + rxt(k,153) &
                      ) * y(k,125) + rxt(k,162)*y(k,219))
         mat(k,1334) = -rxt(k,150)*y(k,112)
         mat(k,1963) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112)
         mat(k,1482) = -rxt(k,162)*y(k,112)
         mat(k,1961) = rxt(k,169)*y(k,126)
         mat(k,1905) = rxt(k,169)*y(k,125)
         mat(k,294) = -(rxt(k,388)*y(k,219))
         mat(k,1478) = -rxt(k,388)*y(k,115)
         mat(k,1063) = .200_r8*rxt(k,380)*y(k,200)
         mat(k,1752) = .200_r8*rxt(k,380)*y(k,101)
         mat(k,896) = -(rxt(k,389)*y(k,219))
         mat(k,1543) = -rxt(k,389)*y(k,116)
         mat(k,1068) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) + rxt(k,379)*y(k,199) &
                      + .800_r8*rxt(k,380)*y(k,200)
         mat(k,1367) = rxt(k,382)*y(k,101)
         mat(k,1920) = rxt(k,383)*y(k,101)
         mat(k,1251) = rxt(k,379)*y(k,101)
         mat(k,1763) = .800_r8*rxt(k,380)*y(k,101)
         mat(k,86) = -(rxt(k,479)*y(k,219))
         mat(k,1444) = -rxt(k,479)*y(k,120)
         mat(k,1391) = -(rxt(k,150)*y(k,112) + rxt(k,159)*y(k,126) + rxt(k,163) &
                      *y(k,205) + rxt(k,164)*y(k,135) + rxt(k,165)*y(k,134) + rxt(k,186) &
                      *y(k,59) + rxt(k,218)*y(k,19) + rxt(k,261)*y(k,200) + rxt(k,270) &
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
         mat(k,326) = -rxt(k,150)*y(k,124)
         mat(k,1945) = -rxt(k,159)*y(k,124)
         mat(k,1735) = -rxt(k,163)*y(k,124)
         mat(k,1629) = -rxt(k,164)*y(k,124)
         mat(k,1839) = -rxt(k,165)*y(k,124)
         mat(k,1888) = -rxt(k,186)*y(k,124)
         mat(k,1809) = -rxt(k,218)*y(k,124)
         mat(k,1785) = -rxt(k,261)*y(k,124)
         mat(k,364) = -rxt(k,270)*y(k,124)
         mat(k,717) = -rxt(k,283)*y(k,124)
         mat(k,1267) = -rxt(k,294)*y(k,124)
         mat(k,606) = -rxt(k,298)*y(k,124)
         mat(k,697) = -rxt(k,311)*y(k,124)
         mat(k,661) = -rxt(k,319)*y(k,124)
         mat(k,1016) = -rxt(k,323)*y(k,124)
         mat(k,483) = -(rxt(k,329) + rxt(k,330)) * y(k,124)
         mat(k,1192) = -(rxt(k,336) + rxt(k,337)) * y(k,124)
         mat(k,1237) = -rxt(k,345)*y(k,124)
         mat(k,575) = -rxt(k,348)*y(k,124)
         mat(k,841) = -(rxt(k,358) + rxt(k,359)) * y(k,124)
         mat(k,1219) = -rxt(k,368)*y(k,124)
         mat(k,1149) = -rxt(k,374)*y(k,124)
         mat(k,1079) = -rxt(k,382)*y(k,124)
         mat(k,1056) = -rxt(k,393)*y(k,124)
         mat(k,428) = -rxt(k,397)*y(k,124)
         mat(k,407) = -rxt(k,400)*y(k,124)
         mat(k,358) = -rxt(k,405)*y(k,124)
         mat(k,544) = -rxt(k,407)*y(k,124)
         mat(k,652) = -rxt(k,411)*y(k,124)
         mat(k,614) = -rxt(k,414)*y(k,124)
         mat(k,773) = -rxt(k,417)*y(k,124)
         mat(k,377) = -rxt(k,420)*y(k,124)
         mat(k,628) = -rxt(k,427)*y(k,124)
         mat(k,645) = -rxt(k,433)*y(k,124)
         mat(k,415) = -rxt(k,436)*y(k,124)
         mat(k,1002) = -rxt(k,447)*y(k,124)
         mat(k,983) = -rxt(k,452)*y(k,124)
         mat(k,963) = -rxt(k,457)*y(k,124)
         mat(k,326) = mat(k,326) + 2.000_r8*rxt(k,152)*y(k,125) + rxt(k,162)*y(k,219)
         mat(k,1986) = 2.000_r8*rxt(k,152)*y(k,112) + rxt(k,155)*y(k,134) + rxt(k,471) &
                      *y(k,151)
         mat(k,1839) = mat(k,1839) + rxt(k,155)*y(k,125)
         mat(k,1114) = rxt(k,471)*y(k,125)
         mat(k,1569) = rxt(k,162)*y(k,112)
         mat(k,1997) = -((rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,112) + (rxt(k,155) &
                      + rxt(k,157)) * y(k,134) + rxt(k,156)*y(k,135) + rxt(k,168) &
                      *y(k,205) + rxt(k,169)*y(k,126) + rxt(k,170)*y(k,219) + rxt(k,188) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,305)*y(k,199) + rxt(k,354) &
                      *y(k,213) + rxt(k,412)*y(k,201) + rxt(k,415)*y(k,212) + rxt(k,418) &
                      *y(k,214) + rxt(k,422)*y(k,142) + rxt(k,425)*y(k,190) + rxt(k,471) &
                      *y(k,151))
         mat(k,329) = -(rxt(k,151) + rxt(k,152) + rxt(k,153)) * y(k,125)
         mat(k,1850) = -(rxt(k,155) + rxt(k,157)) * y(k,125)
         mat(k,1640) = -rxt(k,156)*y(k,125)
         mat(k,1746) = -rxt(k,168)*y(k,125)
         mat(k,1956) = -rxt(k,169)*y(k,125)
         mat(k,1580) = -rxt(k,170)*y(k,125)
         mat(k,1899) = -rxt(k,188)*y(k,125)
         mat(k,1820) = -rxt(k,219)*y(k,125)
         mat(k,1274) = -rxt(k,305)*y(k,125)
         mat(k,1244) = -rxt(k,354)*y(k,125)
         mat(k,656) = -rxt(k,412)*y(k,125)
         mat(k,616) = -rxt(k,415)*y(k,125)
         mat(k,777) = -rxt(k,418)*y(k,125)
         mat(k,394) = -rxt(k,422)*y(k,125)
         mat(k,431) = -rxt(k,425)*y(k,125)
         mat(k,1121) = -rxt(k,471)*y(k,125)
         mat(k,559) = rxt(k,356)*y(k,219)
         mat(k,293) = rxt(k,327)*y(k,126)
         mat(k,1820) = mat(k,1820) + rxt(k,218)*y(k,124)
         mat(k,1899) = mat(k,1899) + rxt(k,186)*y(k,124)
         mat(k,322) = rxt(k,149)*y(k,219)
         mat(k,503) = .700_r8*rxt(k,376)*y(k,219)
         mat(k,1085) = rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126)
         mat(k,1402) = rxt(k,218)*y(k,19) + rxt(k,186)*y(k,59) + rxt(k,382)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,126) + rxt(k,165)*y(k,134) &
                      + rxt(k,164)*y(k,135) + rxt(k,397)*y(k,190) + rxt(k,358) &
                      *y(k,191) + rxt(k,400)*y(k,193) + rxt(k,405)*y(k,195) &
                      + rxt(k,283)*y(k,196) + rxt(k,311)*y(k,197) + rxt(k,407) &
                      *y(k,198) + rxt(k,294)*y(k,199) + rxt(k,261)*y(k,200) &
                      + rxt(k,411)*y(k,201) + rxt(k,329)*y(k,202) + rxt(k,298) &
                      *y(k,204) + rxt(k,163)*y(k,205) + rxt(k,270)*y(k,206) &
                      + .920_r8*rxt(k,368)*y(k,207) + .920_r8*rxt(k,374)*y(k,208) &
                      + rxt(k,336)*y(k,211) + rxt(k,414)*y(k,212) + rxt(k,345) &
                      *y(k,213) + rxt(k,417)*y(k,214) + rxt(k,348)*y(k,215) &
                      + 1.600_r8*rxt(k,447)*y(k,217) + rxt(k,420)*y(k,220) &
                      + rxt(k,319)*y(k,221) + rxt(k,323)*y(k,222) + .900_r8*rxt(k,452) &
                      *y(k,223) + .800_r8*rxt(k,457)*y(k,224) + rxt(k,427)*y(k,225) &
                      + rxt(k,393)*y(k,227) + rxt(k,433)*y(k,228) + rxt(k,436) &
                      *y(k,230)
         mat(k,1956) = mat(k,1956) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,101) &
                      + 2.000_r8*rxt(k,159)*y(k,124) + rxt(k,160)*y(k,134) &
                      + rxt(k,158)*y(k,205) + rxt(k,369)*y(k,207) + rxt(k,375) &
                      *y(k,208) + rxt(k,335)*y(k,211) + rxt(k,346)*y(k,213) &
                      + 2.000_r8*rxt(k,448)*y(k,217) + rxt(k,161)*y(k,219) &
                      + rxt(k,394)*y(k,227)
         mat(k,743) = rxt(k,317)*y(k,219)
         mat(k,1850) = mat(k,1850) + rxt(k,165)*y(k,124) + rxt(k,160)*y(k,126)
         mat(k,1640) = mat(k,1640) + rxt(k,164)*y(k,124)
         mat(k,540) = rxt(k,454)*y(k,219)
         mat(k,431) = mat(k,431) + rxt(k,397)*y(k,124)
         mat(k,846) = rxt(k,358)*y(k,124)
         mat(k,410) = rxt(k,400)*y(k,124)
         mat(k,361) = rxt(k,405)*y(k,124)
         mat(k,722) = rxt(k,283)*y(k,124)
         mat(k,702) = rxt(k,311)*y(k,124)
         mat(k,548) = rxt(k,407)*y(k,124)
         mat(k,1274) = mat(k,1274) + rxt(k,294)*y(k,124)
         mat(k,1796) = rxt(k,261)*y(k,124) + .500_r8*rxt(k,445)*y(k,217)
         mat(k,656) = mat(k,656) + rxt(k,411)*y(k,124)
         mat(k,486) = rxt(k,329)*y(k,124)
         mat(k,610) = rxt(k,298)*y(k,124)
         mat(k,1746) = mat(k,1746) + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126)
         mat(k,367) = rxt(k,270)*y(k,124)
         mat(k,1226) = .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126)
         mat(k,1156) = .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126)
         mat(k,1198) = rxt(k,336)*y(k,124) + rxt(k,335)*y(k,126)
         mat(k,616) = mat(k,616) + rxt(k,414)*y(k,124)
         mat(k,1244) = mat(k,1244) + rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126)
         mat(k,777) = mat(k,777) + rxt(k,417)*y(k,124)
         mat(k,578) = rxt(k,348)*y(k,124)
         mat(k,1008) = 1.600_r8*rxt(k,447)*y(k,124) + 2.000_r8*rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1580) = mat(k,1580) + rxt(k,356)*y(k,1) + rxt(k,149)*y(k,90) &
                      + .700_r8*rxt(k,376)*y(k,99) + rxt(k,161)*y(k,126) + rxt(k,317) &
                      *y(k,127) + rxt(k,454)*y(k,177)
         mat(k,380) = rxt(k,420)*y(k,124)
         mat(k,665) = rxt(k,319)*y(k,124)
         mat(k,1021) = rxt(k,323)*y(k,124)
         mat(k,988) = .900_r8*rxt(k,452)*y(k,124)
         mat(k,969) = .800_r8*rxt(k,457)*y(k,124)
         mat(k,631) = rxt(k,427)*y(k,124)
         mat(k,1062) = rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126)
         mat(k,648) = rxt(k,433)*y(k,124)
         mat(k,418) = rxt(k,436)*y(k,124)
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
         mat(k,1955) = -(rxt(k,158)*y(k,205) + rxt(k,159)*y(k,124) + rxt(k,160) &
                      *y(k,134) + rxt(k,161)*y(k,219) + rxt(k,169)*y(k,125) + rxt(k,255) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,211) + rxt(k,346) &
                      *y(k,213) + rxt(k,369)*y(k,207) + rxt(k,375)*y(k,208) + rxt(k,378) &
                      *y(k,98) + rxt(k,383)*y(k,101) + rxt(k,394)*y(k,227) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,110) + rxt(k,448)*y(k,217) + rxt(k,459) &
                      *y(k,179) + rxt(k,462)*y(k,67))
         mat(k,1745) = -rxt(k,158)*y(k,126)
         mat(k,1401) = -rxt(k,159)*y(k,126)
         mat(k,1849) = -rxt(k,160)*y(k,126)
         mat(k,1579) = -rxt(k,161)*y(k,126)
         mat(k,1996) = -rxt(k,169)*y(k,126)
         mat(k,1872) = -rxt(k,255)*y(k,126)
         mat(k,1031) = -rxt(k,288)*y(k,126)
         mat(k,889) = -rxt(k,307)*y(k,126)
         mat(k,1106) = -rxt(k,314)*y(k,126)
         mat(k,292) = -rxt(k,327)*y(k,126)
         mat(k,1197) = -rxt(k,335)*y(k,126)
         mat(k,1243) = -rxt(k,346)*y(k,126)
         mat(k,1225) = -rxt(k,369)*y(k,126)
         mat(k,1155) = -rxt(k,375)*y(k,126)
         mat(k,762) = -rxt(k,378)*y(k,126)
         mat(k,1084) = -rxt(k,383)*y(k,126)
         mat(k,1061) = -rxt(k,394)*y(k,126)
         mat(k,830) = -rxt(k,439)*y(k,126)
         mat(k,803) = -rxt(k,442)*y(k,126)
         mat(k,1007) = -rxt(k,448)*y(k,126)
         mat(k,872) = -rxt(k,459)*y(k,126)
         mat(k,238) = -rxt(k,462)*y(k,126)
         mat(k,457) = rxt(k,220)*y(k,134)
         mat(k,2030) = rxt(k,187)*y(k,60)
         mat(k,860) = rxt(k,187)*y(k,56) + rxt(k,189)*y(k,134) + rxt(k,190)*y(k,219)
         mat(k,683) = rxt(k,234)*y(k,89)
         mat(k,1297) = rxt(k,234)*y(k,73) + rxt(k,171)*y(k,219)
         mat(k,465) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,1996) = mat(k,1996) + rxt(k,157)*y(k,134) + rxt(k,156)*y(k,135)
         mat(k,1849) = mat(k,1849) + rxt(k,220)*y(k,20) + rxt(k,189)*y(k,60) &
                      + rxt(k,157)*y(k,125)
         mat(k,1639) = rxt(k,156)*y(k,125)
         mat(k,389) = rxt(k,303)*y(k,219)
         mat(k,1579) = mat(k,1579) + rxt(k,190)*y(k,60) + rxt(k,171)*y(k,89) &
                      + .500_r8*rxt(k,351)*y(k,109) + rxt(k,303)*y(k,140)
         mat(k,738) = -(rxt(k,317)*y(k,219))
         mat(k,1532) = -rxt(k,317)*y(k,127)
         mat(k,876) = rxt(k,307)*y(k,126)
         mat(k,488) = .500_r8*rxt(k,377)*y(k,219)
         mat(k,338) = rxt(k,384)*y(k,219)
         mat(k,295) = rxt(k,388)*y(k,219)
         mat(k,893) = rxt(k,389)*y(k,219)
         mat(k,1911) = rxt(k,307)*y(k,29)
         mat(k,1532) = mat(k,1532) + .500_r8*rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) &
                      + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116)
         mat(k,312) = -(rxt(k,449)*y(k,219))
         mat(k,1480) = -rxt(k,449)*y(k,128)
         mat(k,1666) = rxt(k,446)*y(k,217)
         mat(k,990) = rxt(k,446)*y(k,205)
         mat(k,1846) = -(rxt(k,129)*y(k,135) + 4._r8*rxt(k,130)*y(k,134) + rxt(k,132) &
                      *y(k,77) + rxt(k,133)*y(k,79) + rxt(k,138)*y(k,205) + rxt(k,144) &
                      *y(k,219) + (rxt(k,155) + rxt(k,157)) * y(k,125) + rxt(k,160) &
                      *y(k,126) + rxt(k,165)*y(k,124) + rxt(k,189)*y(k,60) + rxt(k,191) &
                      *y(k,59) + rxt(k,194)*y(k,85) + rxt(k,197)*y(k,92) + rxt(k,220) &
                      *y(k,20) + rxt(k,221)*y(k,19) + rxt(k,223)*y(k,81) + rxt(k,225) &
                      *y(k,91) + rxt(k,256)*y(k,42) + rxt(k,464)*y(k,138))
         mat(k,1636) = -rxt(k,129)*y(k,134)
         mat(k,1041) = -rxt(k,132)*y(k,134)
         mat(k,527) = -rxt(k,133)*y(k,134)
         mat(k,1742) = -rxt(k,138)*y(k,134)
         mat(k,1576) = -rxt(k,144)*y(k,134)
         mat(k,1993) = -(rxt(k,155) + rxt(k,157)) * y(k,134)
         mat(k,1952) = -rxt(k,160)*y(k,134)
         mat(k,1398) = -rxt(k,165)*y(k,134)
         mat(k,858) = -rxt(k,189)*y(k,134)
         mat(k,1895) = -rxt(k,191)*y(k,134)
         mat(k,2049) = -rxt(k,194)*y(k,134)
         mat(k,726) = -rxt(k,197)*y(k,134)
         mat(k,456) = -rxt(k,220)*y(k,134)
         mat(k,1816) = -rxt(k,221)*y(k,134)
         mat(k,710) = -rxt(k,223)*y(k,134)
         mat(k,672) = -rxt(k,225)*y(k,134)
         mat(k,1869) = -rxt(k,256)*y(k,134)
         mat(k,285) = -rxt(k,464)*y(k,134)
         mat(k,1310) = rxt(k,136)*y(k,205)
         mat(k,328) = rxt(k,150)*y(k,124) + rxt(k,151)*y(k,125)
         mat(k,1398) = mat(k,1398) + rxt(k,150)*y(k,112)
         mat(k,1993) = mat(k,1993) + rxt(k,151)*y(k,112)
         mat(k,1742) = mat(k,1742) + rxt(k,136)*y(k,76)
         mat(k,1576) = mat(k,1576) + 2.000_r8*rxt(k,146)*y(k,219)
         mat(k,1632) = -(rxt(k,128)*y(k,218) + rxt(k,129)*y(k,134) + rxt(k,139) &
                      *y(k,205) + rxt(k,140)*y(k,76) + rxt(k,145)*y(k,219) + rxt(k,156) &
                      *y(k,125) + rxt(k,164)*y(k,124) + rxt(k,180)*y(k,56) + rxt(k,212) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,105) + rxt(k,352)*y(k,111) + rxt(k,385)*y(k,98) + rxt(k,423) &
                      *y(k,142) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,110) + rxt(k,467) &
                      *y(k,149) + rxt(k,473)*y(k,151))
         mat(k,1418) = -rxt(k,128)*y(k,135)
         mat(k,1842) = -rxt(k,129)*y(k,135)
         mat(k,1738) = -rxt(k,139)*y(k,135)
         mat(k,1307) = -rxt(k,140)*y(k,135)
         mat(k,1572) = -rxt(k,145)*y(k,135)
         mat(k,1989) = -rxt(k,156)*y(k,135)
         mat(k,1394) = -rxt(k,164)*y(k,135)
         mat(k,2023) = -rxt(k,180)*y(k,135)
         mat(k,1282) = -rxt(k,212)*y(k,135)
         mat(k,472) = -rxt(k,279)*y(k,135)
         mat(k,885) = -rxt(k,308)*y(k,135)
         mat(k,1095) = -rxt(k,338)*y(k,135)
         mat(k,1171) = -rxt(k,352)*y(k,135)
         mat(k,758) = -rxt(k,385)*y(k,135)
         mat(k,393) = -rxt(k,423)*y(k,135)
         mat(k,827) = -rxt(k,440)*y(k,135)
         mat(k,800) = -rxt(k,443)*y(k,135)
         mat(k,423) = -rxt(k,467)*y(k,135)
         mat(k,1116) = -rxt(k,473)*y(k,135)
         mat(k,1269) = .150_r8*rxt(k,293)*y(k,205)
         mat(k,1738) = mat(k,1738) + .150_r8*rxt(k,293)*y(k,199) + .150_r8*rxt(k,343) &
                      *y(k,213)
         mat(k,1239) = .150_r8*rxt(k,343)*y(k,205)
         mat(k,254) = -(rxt(k,474)*y(k,151))
         mat(k,1108) = -rxt(k,474)*y(k,137)
         mat(k,1801) = rxt(k,214)*y(k,59)
         mat(k,1880) = rxt(k,214)*y(k,19) + 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,278) = -(rxt(k,464)*y(k,134) + rxt(k,465)*y(k,219))
         mat(k,1824) = -rxt(k,464)*y(k,138)
         mat(k,1476) = -rxt(k,465)*y(k,138)
         mat(k,913) = rxt(k,331)*y(k,219)
         mat(k,1328) = .100_r8*rxt(k,452)*y(k,223)
         mat(k,1462) = rxt(k,331)*y(k,93)
         mat(k,971) = .100_r8*rxt(k,452)*y(k,124)
         mat(k,384) = -(rxt(k,303)*y(k,219))
         mat(k,1491) = -rxt(k,303)*y(k,140)
         mat(k,1964) = rxt(k,305)*y(k,199)
         mat(k,1247) = rxt(k,305)*y(k,125)
         mat(k,1960) = rxt(k,425)*y(k,190)
         mat(k,425) = rxt(k,425)*y(k,125)
         mat(k,391) = -(rxt(k,422)*y(k,125) + rxt(k,423)*y(k,135))
         mat(k,1965) = -rxt(k,422)*y(k,142)
         mat(k,1591) = -rxt(k,423)*y(k,142)
         mat(k,146) = .070_r8*rxt(k,409)*y(k,219)
         mat(k,1339) = rxt(k,407)*y(k,198)
         mat(k,129) = .060_r8*rxt(k,421)*y(k,219)
         mat(k,171) = .070_r8*rxt(k,437)*y(k,219)
         mat(k,542) = rxt(k,407)*y(k,124)
         mat(k,1492) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,143) &
                      + .070_r8*rxt(k,437)*y(k,186)
         mat(k,127) = -(rxt(k,421)*y(k,219))
         mat(k,1451) = -rxt(k,421)*y(k,143)
         mat(k,119) = .530_r8*rxt(k,398)*y(k,219)
         mat(k,1451) = mat(k,1451) + .530_r8*rxt(k,398)*y(k,7)
         mat(k,259) = -(rxt(k,424)*y(k,219))
         mat(k,1472) = -rxt(k,424)*y(k,144)
         mat(k,1662) = rxt(k,419)*y(k,220)
         mat(k,374) = rxt(k,419)*y(k,205)
         mat(k,443) = -(rxt(k,320)*y(k,219))
         mat(k,1500) = -rxt(k,320)*y(k,147)
         mat(k,1683) = rxt(k,318)*y(k,221)
         mat(k,657) = rxt(k,318)*y(k,205)
         mat(k,330) = -(rxt(k,324)*y(k,219))
         mat(k,1483) = -rxt(k,324)*y(k,148)
         mat(k,1668) = .850_r8*rxt(k,322)*y(k,222)
         mat(k,1010) = .850_r8*rxt(k,322)*y(k,205)
         mat(k,419) = -(rxt(k,467)*y(k,135) + rxt(k,470)*y(k,219))
         mat(k,1592) = -rxt(k,467)*y(k,149)
         mat(k,1496) = -rxt(k,470)*y(k,149)
         mat(k,1111) = -(rxt(k,468)*y(k,19) + rxt(k,469)*y(k,59) + rxt(k,471)*y(k,125) &
                      + rxt(k,473)*y(k,135) + rxt(k,474)*y(k,137) + rxt(k,475) &
                      *y(k,219))
         mat(k,1805) = -rxt(k,468)*y(k,151)
         mat(k,1884) = -rxt(k,469)*y(k,151)
         mat(k,1980) = -rxt(k,471)*y(k,151)
         mat(k,1619) = -rxt(k,473)*y(k,151)
         mat(k,256) = -rxt(k,474)*y(k,151)
         mat(k,1559) = -rxt(k,475)*y(k,151)
         mat(k,1835) = rxt(k,464)*y(k,138)
         mat(k,1619) = mat(k,1619) + rxt(k,467)*y(k,149)
         mat(k,282) = rxt(k,464)*y(k,134)
         mat(k,420) = rxt(k,467)*y(k,135) + rxt(k,470)*y(k,219)
         mat(k,1559) = mat(k,1559) + rxt(k,470)*y(k,149)
         mat(k,732) = -(rxt(k,477)*y(k,219))
         mat(k,1531) = -rxt(k,477)*y(k,152)
         mat(k,1804) = rxt(k,468)*y(k,151)
         mat(k,1882) = rxt(k,469)*y(k,151)
         mat(k,234) = rxt(k,462)*y(k,126) + (rxt(k,463)+.500_r8*rxt(k,476))*y(k,219)
         mat(k,1973) = rxt(k,471)*y(k,151)
         mat(k,1910) = rxt(k,462)*y(k,67)
         mat(k,1598) = rxt(k,473)*y(k,151)
         mat(k,255) = rxt(k,474)*y(k,151)
         mat(k,280) = rxt(k,465)*y(k,219)
         mat(k,1110) = rxt(k,468)*y(k,19) + rxt(k,469)*y(k,59) + rxt(k,471)*y(k,125) &
                      + rxt(k,473)*y(k,135) + rxt(k,474)*y(k,137) + rxt(k,475) &
                      *y(k,219)
         mat(k,1531) = mat(k,1531) + (rxt(k,463)+.500_r8*rxt(k,476))*y(k,67) &
                      + rxt(k,465)*y(k,138) + rxt(k,475)*y(k,151)
         mat(k,196) = -(rxt(k,478)*y(k,231))
         mat(k,2059) = -rxt(k,478)*y(k,153)
         mat(k,731) = rxt(k,477)*y(k,219)
         mat(k,1464) = rxt(k,477)*y(k,152)
         mat(k,805) = .2202005_r8*rxt(k,497)*y(k,135)
         mat(k,778) = .0508005_r8*rxt(k,513)*y(k,135)
         mat(k,1316) = .1279005_r8*rxt(k,496)*y(k,192) + .0097005_r8*rxt(k,501) &
                      *y(k,194) + .0003005_r8*rxt(k,504)*y(k,209) &
                      + .1056005_r8*rxt(k,508)*y(k,210) + .0245005_r8*rxt(k,512) &
                      *y(k,216) + .0154005_r8*rxt(k,518)*y(k,226) &
                      + .0063005_r8*rxt(k,521)*y(k,229)
         mat(k,1584) = .2202005_r8*rxt(k,497)*y(k,6) + .0508005_r8*rxt(k,513)*y(k,110)
         mat(k,36) = .5931005_r8*rxt(k,515)*y(k,219)
         mat(k,42) = .1279005_r8*rxt(k,496)*y(k,124) + .2202005_r8*rxt(k,495)*y(k,205)
         mat(k,48) = .0097005_r8*rxt(k,501)*y(k,124) + .0023005_r8*rxt(k,500)*y(k,205)
         mat(k,1644) = .2202005_r8*rxt(k,495)*y(k,192) + .0023005_r8*rxt(k,500) &
                      *y(k,194) + .0031005_r8*rxt(k,503)*y(k,209) &
                      + .2381005_r8*rxt(k,507)*y(k,210) + .0508005_r8*rxt(k,511) &
                      *y(k,216) + .1364005_r8*rxt(k,517)*y(k,226) &
                      + .1677005_r8*rxt(k,520)*y(k,229)
         mat(k,54) = .0003005_r8*rxt(k,504)*y(k,124) + .0031005_r8*rxt(k,503)*y(k,205)
         mat(k,60) = .1056005_r8*rxt(k,508)*y(k,124) + .2381005_r8*rxt(k,507)*y(k,205)
         mat(k,68) = .0245005_r8*rxt(k,512)*y(k,124) + .0508005_r8*rxt(k,511)*y(k,205)
         mat(k,1430) = .5931005_r8*rxt(k,515)*y(k,174)
         mat(k,74) = .0154005_r8*rxt(k,518)*y(k,124) + .1364005_r8*rxt(k,517)*y(k,205)
         mat(k,80) = .0063005_r8*rxt(k,521)*y(k,124) + .1677005_r8*rxt(k,520)*y(k,205)
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
         mat(k,806) = .2067005_r8*rxt(k,497)*y(k,135)
         mat(k,779) = .1149005_r8*rxt(k,513)*y(k,135)
         mat(k,1317) = .1792005_r8*rxt(k,496)*y(k,192) + .0034005_r8*rxt(k,501) &
                      *y(k,194) + .0003005_r8*rxt(k,504)*y(k,209) &
                      + .1026005_r8*rxt(k,508)*y(k,210) + .0082005_r8*rxt(k,512) &
                      *y(k,216) + .0452005_r8*rxt(k,518)*y(k,226) &
                      + .0237005_r8*rxt(k,521)*y(k,229)
         mat(k,1585) = .2067005_r8*rxt(k,497)*y(k,6) + .1149005_r8*rxt(k,513)*y(k,110)
         mat(k,37) = .1534005_r8*rxt(k,515)*y(k,219)
         mat(k,43) = .1792005_r8*rxt(k,496)*y(k,124) + .2067005_r8*rxt(k,495)*y(k,205)
         mat(k,49) = .0034005_r8*rxt(k,501)*y(k,124) + .0008005_r8*rxt(k,500)*y(k,205)
         mat(k,1645) = .2067005_r8*rxt(k,495)*y(k,192) + .0008005_r8*rxt(k,500) &
                      *y(k,194) + .0035005_r8*rxt(k,503)*y(k,209) &
                      + .1308005_r8*rxt(k,507)*y(k,210) + .1149005_r8*rxt(k,511) &
                      *y(k,216) + .0101005_r8*rxt(k,517)*y(k,226) &
                      + .0174005_r8*rxt(k,520)*y(k,229)
         mat(k,55) = .0003005_r8*rxt(k,504)*y(k,124) + .0035005_r8*rxt(k,503)*y(k,205)
         mat(k,61) = .1026005_r8*rxt(k,508)*y(k,124) + .1308005_r8*rxt(k,507)*y(k,205)
         mat(k,69) = .0082005_r8*rxt(k,512)*y(k,124) + .1149005_r8*rxt(k,511)*y(k,205)
         mat(k,1431) = .1534005_r8*rxt(k,515)*y(k,174)
         mat(k,75) = .0452005_r8*rxt(k,518)*y(k,124) + .0101005_r8*rxt(k,517)*y(k,205)
         mat(k,81) = .0237005_r8*rxt(k,521)*y(k,124) + .0174005_r8*rxt(k,520)*y(k,205)
         mat(k,807) = .0653005_r8*rxt(k,497)*y(k,135)
         mat(k,780) = .0348005_r8*rxt(k,513)*y(k,135)
         mat(k,1318) = .0676005_r8*rxt(k,496)*y(k,192) + .1579005_r8*rxt(k,501) &
                      *y(k,194) + .0073005_r8*rxt(k,504)*y(k,209) &
                      + .0521005_r8*rxt(k,508)*y(k,210) + .0772005_r8*rxt(k,512) &
                      *y(k,216) + .0966005_r8*rxt(k,518)*y(k,226) &
                      + .0025005_r8*rxt(k,521)*y(k,229)
         mat(k,1586) = .0653005_r8*rxt(k,497)*y(k,6) + .0348005_r8*rxt(k,513)*y(k,110)
         mat(k,38) = .0459005_r8*rxt(k,515)*y(k,219)
         mat(k,44) = .0676005_r8*rxt(k,496)*y(k,124) + .0653005_r8*rxt(k,495)*y(k,205)
         mat(k,50) = .1579005_r8*rxt(k,501)*y(k,124) + .0843005_r8*rxt(k,500)*y(k,205)
         mat(k,1646) = .0653005_r8*rxt(k,495)*y(k,192) + .0843005_r8*rxt(k,500) &
                      *y(k,194) + .0003005_r8*rxt(k,503)*y(k,209) &
                      + .0348005_r8*rxt(k,507)*y(k,210) + .0348005_r8*rxt(k,511) &
                      *y(k,216) + .0763005_r8*rxt(k,517)*y(k,226) + .086_r8*rxt(k,520) &
                      *y(k,229)
         mat(k,56) = .0073005_r8*rxt(k,504)*y(k,124) + .0003005_r8*rxt(k,503)*y(k,205)
         mat(k,62) = .0521005_r8*rxt(k,508)*y(k,124) + .0348005_r8*rxt(k,507)*y(k,205)
         mat(k,70) = .0772005_r8*rxt(k,512)*y(k,124) + .0348005_r8*rxt(k,511)*y(k,205)
         mat(k,1432) = .0459005_r8*rxt(k,515)*y(k,174)
         mat(k,76) = .0966005_r8*rxt(k,518)*y(k,124) + .0763005_r8*rxt(k,517)*y(k,205)
         mat(k,82) = .0025005_r8*rxt(k,521)*y(k,124) + .086_r8*rxt(k,520)*y(k,205)
         mat(k,808) = .1749305_r8*rxt(k,494)*y(k,126) + .1284005_r8*rxt(k,497) &
                      *y(k,135)
         mat(k,744) = .0590245_r8*rxt(k,502)*y(k,126) + .0033005_r8*rxt(k,505) &
                      *y(k,135)
         mat(k,781) = .1749305_r8*rxt(k,510)*y(k,126) + .0554005_r8*rxt(k,513) &
                      *y(k,135)
         mat(k,1319) = .079_r8*rxt(k,496)*y(k,192) + .0059005_r8*rxt(k,501)*y(k,194) &
                      + .0057005_r8*rxt(k,504)*y(k,209) + .0143005_r8*rxt(k,508) &
                      *y(k,210) + .0332005_r8*rxt(k,512)*y(k,216) &
                      + .0073005_r8*rxt(k,518)*y(k,226) + .011_r8*rxt(k,521)*y(k,229)
         mat(k,1903) = .1749305_r8*rxt(k,494)*y(k,6) + .0590245_r8*rxt(k,502)*y(k,98) &
                      + .1749305_r8*rxt(k,510)*y(k,110)
         mat(k,1587) = .1284005_r8*rxt(k,497)*y(k,6) + .0033005_r8*rxt(k,505)*y(k,98) &
                      + .0554005_r8*rxt(k,513)*y(k,110)
         mat(k,39) = .0085005_r8*rxt(k,515)*y(k,219)
         mat(k,45) = .079_r8*rxt(k,496)*y(k,124) + .1284005_r8*rxt(k,495)*y(k,205)
         mat(k,51) = .0059005_r8*rxt(k,501)*y(k,124) + .0443005_r8*rxt(k,500)*y(k,205)
         mat(k,1647) = .1284005_r8*rxt(k,495)*y(k,192) + .0443005_r8*rxt(k,500) &
                      *y(k,194) + .0271005_r8*rxt(k,503)*y(k,209) &
                      + .0076005_r8*rxt(k,507)*y(k,210) + .0554005_r8*rxt(k,511) &
                      *y(k,216) + .2157005_r8*rxt(k,517)*y(k,226) &
                      + .0512005_r8*rxt(k,520)*y(k,229)
         mat(k,57) = .0057005_r8*rxt(k,504)*y(k,124) + .0271005_r8*rxt(k,503)*y(k,205)
         mat(k,63) = .0143005_r8*rxt(k,508)*y(k,124) + .0076005_r8*rxt(k,507)*y(k,205)
         mat(k,71) = .0332005_r8*rxt(k,512)*y(k,124) + .0554005_r8*rxt(k,511)*y(k,205)
         mat(k,1433) = .0085005_r8*rxt(k,515)*y(k,174)
         mat(k,77) = .0073005_r8*rxt(k,518)*y(k,124) + .2157005_r8*rxt(k,517)*y(k,205)
         mat(k,83) = .011_r8*rxt(k,521)*y(k,124) + .0512005_r8*rxt(k,520)*y(k,205)
         mat(k,809) = .5901905_r8*rxt(k,494)*y(k,126) + .114_r8*rxt(k,497)*y(k,135)
         mat(k,745) = .0250245_r8*rxt(k,502)*y(k,126)
         mat(k,782) = .5901905_r8*rxt(k,510)*y(k,126) + .1278005_r8*rxt(k,513) &
                      *y(k,135)
         mat(k,1320) = .1254005_r8*rxt(k,496)*y(k,192) + .0536005_r8*rxt(k,501) &
                      *y(k,194) + .0623005_r8*rxt(k,504)*y(k,209) &
                      + .0166005_r8*rxt(k,508)*y(k,210) + .130_r8*rxt(k,512)*y(k,216) &
                      + .238_r8*rxt(k,518)*y(k,226) + .1185005_r8*rxt(k,521)*y(k,229)
         mat(k,1904) = .5901905_r8*rxt(k,494)*y(k,6) + .0250245_r8*rxt(k,502)*y(k,98) &
                      + .5901905_r8*rxt(k,510)*y(k,110)
         mat(k,1588) = .114_r8*rxt(k,497)*y(k,6) + .1278005_r8*rxt(k,513)*y(k,110)
         mat(k,40) = .0128005_r8*rxt(k,515)*y(k,219)
         mat(k,46) = .1254005_r8*rxt(k,496)*y(k,124) + .114_r8*rxt(k,495)*y(k,205)
         mat(k,52) = .0536005_r8*rxt(k,501)*y(k,124) + .1621005_r8*rxt(k,500)*y(k,205)
         mat(k,1648) = .114_r8*rxt(k,495)*y(k,192) + .1621005_r8*rxt(k,500)*y(k,194) &
                      + .0474005_r8*rxt(k,503)*y(k,209) + .0113005_r8*rxt(k,507) &
                      *y(k,210) + .1278005_r8*rxt(k,511)*y(k,216) &
                      + .0738005_r8*rxt(k,517)*y(k,226) + .1598005_r8*rxt(k,520) &
                      *y(k,229)
         mat(k,58) = .0623005_r8*rxt(k,504)*y(k,124) + .0474005_r8*rxt(k,503)*y(k,205)
         mat(k,64) = .0166005_r8*rxt(k,508)*y(k,124) + .0113005_r8*rxt(k,507)*y(k,205)
         mat(k,72) = .130_r8*rxt(k,512)*y(k,124) + .1278005_r8*rxt(k,511)*y(k,205)
         mat(k,1434) = .0128005_r8*rxt(k,515)*y(k,174)
         mat(k,78) = .238_r8*rxt(k,518)*y(k,124) + .0738005_r8*rxt(k,517)*y(k,205)
         mat(k,84) = .1185005_r8*rxt(k,521)*y(k,124) + .1598005_r8*rxt(k,520)*y(k,205)
         mat(k,41) = -(rxt(k,515)*y(k,219))
         mat(k,1435) = -rxt(k,515)*y(k,174)
         mat(k,139) = .100_r8*rxt(k,429)*y(k,219)
         mat(k,161) = .230_r8*rxt(k,431)*y(k,219)
         mat(k,1455) = .100_r8*rxt(k,429)*y(k,182) + .230_r8*rxt(k,431)*y(k,184)
         mat(k,504) = -(rxt(k,453)*y(k,219))
         mat(k,1507) = -rxt(k,453)*y(k,176)
         mat(k,1685) = rxt(k,451)*y(k,223)
         mat(k,972) = rxt(k,451)*y(k,205)
         mat(k,535) = -(rxt(k,454)*y(k,219))
         mat(k,1511) = -rxt(k,454)*y(k,177)
         mat(k,1348) = .200_r8*rxt(k,447)*y(k,217) + .200_r8*rxt(k,457)*y(k,224)
         mat(k,1755) = .500_r8*rxt(k,445)*y(k,217)
         mat(k,991) = .200_r8*rxt(k,447)*y(k,124) + .500_r8*rxt(k,445)*y(k,200)
         mat(k,950) = .200_r8*rxt(k,457)*y(k,124)
         mat(k,395) = -(rxt(k,458)*y(k,219))
         mat(k,1493) = -rxt(k,458)*y(k,178)
         mat(k,1678) = rxt(k,456)*y(k,224)
         mat(k,949) = rxt(k,456)*y(k,205)
         mat(k,865) = -(rxt(k,459)*y(k,126) + rxt(k,460)*y(k,219))
         mat(k,1918) = -rxt(k,459)*y(k,179)
         mat(k,1541) = -rxt(k,460)*y(k,179)
         mat(k,818) = .330_r8*rxt(k,440)*y(k,135)
         mat(k,791) = .330_r8*rxt(k,443)*y(k,135)
         mat(k,1366) = .800_r8*rxt(k,447)*y(k,217) + .800_r8*rxt(k,457)*y(k,224)
         mat(k,1918) = mat(k,1918) + rxt(k,448)*y(k,217)
         mat(k,1605) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,110)
         mat(k,536) = rxt(k,454)*y(k,219)
         mat(k,1762) = .500_r8*rxt(k,445)*y(k,217) + rxt(k,455)*y(k,224)
         mat(k,993) = .800_r8*rxt(k,447)*y(k,124) + rxt(k,448)*y(k,126) &
                      + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1541) = mat(k,1541) + rxt(k,454)*y(k,177)
         mat(k,953) = .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,200)
         mat(k,930) = -(rxt(k,461)*y(k,219))
         mat(k,1546) = -rxt(k,461)*y(k,180)
         mat(k,819) = .300_r8*rxt(k,440)*y(k,135)
         mat(k,792) = .300_r8*rxt(k,443)*y(k,135)
         mat(k,1370) = .900_r8*rxt(k,452)*y(k,223)
         mat(k,1608) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,110)
         mat(k,1765) = rxt(k,450)*y(k,223)
         mat(k,976) = .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,200)
         mat(k,515) = -(rxt(k,428)*y(k,219))
         mat(k,1508) = -rxt(k,428)*y(k,181)
         mat(k,1686) = rxt(k,426)*y(k,225)
         mat(k,619) = rxt(k,426)*y(k,205)
         mat(k,137) = -(rxt(k,429)*y(k,219))
         mat(k,1453) = -rxt(k,429)*y(k,182)
         mat(k,153) = -(rxt(k,395)*y(k,219))
         mat(k,1456) = -rxt(k,395)*y(k,183)
         mat(k,1657) = rxt(k,392)*y(k,227)
         mat(k,1046) = rxt(k,392)*y(k,205)
         mat(k,162) = -(rxt(k,431)*y(k,219))
         mat(k,1458) = -rxt(k,431)*y(k,184)
         mat(k,591) = -(rxt(k,434)*y(k,219))
         mat(k,1517) = -rxt(k,434)*y(k,185)
         mat(k,1693) = rxt(k,432)*y(k,228)
         mat(k,636) = rxt(k,432)*y(k,205)
         mat(k,170) = -(rxt(k,437)*y(k,219))
         mat(k,1459) = -rxt(k,437)*y(k,186)
         mat(k,163) = .150_r8*rxt(k,431)*y(k,219)
         mat(k,1459) = mat(k,1459) + .150_r8*rxt(k,431)*y(k,184)
         mat(k,348) = -(rxt(k,438)*y(k,219))
         mat(k,1486) = -rxt(k,438)*y(k,187)
         mat(k,1671) = rxt(k,435)*y(k,230)
         mat(k,411) = rxt(k,435)*y(k,205)
         mat(k,426) = -(rxt(k,396)*y(k,205) + rxt(k,397)*y(k,124) + rxt(k,425) &
                      *y(k,125))
         mat(k,1681) = -rxt(k,396)*y(k,190)
         mat(k,1342) = -rxt(k,397)*y(k,190)
         mat(k,1966) = -rxt(k,425)*y(k,190)
         mat(k,190) = rxt(k,402)*y(k,219)
         mat(k,1497) = rxt(k,402)*y(k,22)
         mat(k,837) = -(rxt(k,357)*y(k,205) + (rxt(k,358) + rxt(k,359)) * y(k,124))
         mat(k,1709) = -rxt(k,357)*y(k,191)
         mat(k,1364) = -(rxt(k,358) + rxt(k,359)) * y(k,191)
         mat(k,564) = rxt(k,360)*y(k,219)
         mat(k,181) = rxt(k,361)*y(k,219)
         mat(k,1538) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)
         mat(k,47) = -(rxt(k,495)*y(k,205) + rxt(k,496)*y(k,124))
         mat(k,1649) = -rxt(k,495)*y(k,192)
         mat(k,1321) = -rxt(k,496)*y(k,192)
         mat(k,810) = rxt(k,498)*y(k,219)
         mat(k,1436) = rxt(k,498)*y(k,6)
         mat(k,404) = -(rxt(k,399)*y(k,205) + rxt(k,400)*y(k,124))
         mat(k,1679) = -rxt(k,399)*y(k,193)
         mat(k,1340) = -rxt(k,400)*y(k,193)
         mat(k,120) = .350_r8*rxt(k,398)*y(k,219)
         mat(k,302) = rxt(k,401)*y(k,219)
         mat(k,1494) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)
         mat(k,53) = -(rxt(k,500)*y(k,205) + rxt(k,501)*y(k,124))
         mat(k,1650) = -rxt(k,500)*y(k,194)
         mat(k,1322) = -rxt(k,501)*y(k,194)
         mat(k,116) = rxt(k,499)*y(k,219)
         mat(k,1437) = rxt(k,499)*y(k,7)
         mat(k,356) = -(rxt(k,403)*y(k,205) + rxt(k,405)*y(k,124))
         mat(k,1672) = -rxt(k,403)*y(k,195)
         mat(k,1335) = -rxt(k,405)*y(k,195)
         mat(k,266) = rxt(k,404)*y(k,219)
         mat(k,140) = .070_r8*rxt(k,429)*y(k,219)
         mat(k,164) = .060_r8*rxt(k,431)*y(k,219)
         mat(k,1487) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,182) &
                      + .060_r8*rxt(k,431)*y(k,184)
         mat(k,715) = -(4._r8*rxt(k,280)*y(k,196) + rxt(k,281)*y(k,200) + rxt(k,282) &
                      *y(k,205) + rxt(k,283)*y(k,124))
         mat(k,1758) = -rxt(k,281)*y(k,196)
         mat(k,1704) = -rxt(k,282)*y(k,196)
         mat(k,1360) = -rxt(k,283)*y(k,196)
         mat(k,274) = .500_r8*rxt(k,285)*y(k,219)
         mat(k,228) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,219)
         mat(k,2008) = rxt(k,286)*y(k,28)
         mat(k,1529) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)
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
         mat(k,694) = -(rxt(k,309)*y(k,200) + rxt(k,310)*y(k,205) + rxt(k,311) &
                      *y(k,124))
         mat(k,1757) = -rxt(k,309)*y(k,197)
         mat(k,1702) = -rxt(k,310)*y(k,197)
         mat(k,1359) = -rxt(k,311)*y(k,197)
         mat(k,343) = rxt(k,312)*y(k,219)
         mat(k,96) = rxt(k,313)*y(k,219)
         mat(k,1527) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)
         mat(k,543) = -(rxt(k,406)*y(k,205) + rxt(k,407)*y(k,124))
         mat(k,1689) = -rxt(k,406)*y(k,198)
         mat(k,1349) = -rxt(k,407)*y(k,198)
         mat(k,206) = rxt(k,408)*y(k,219)
         mat(k,1349) = mat(k,1349) + rxt(k,397)*y(k,190)
         mat(k,1595) = rxt(k,423)*y(k,142)
         mat(k,392) = rxt(k,423)*y(k,135)
         mat(k,427) = rxt(k,397)*y(k,124) + .400_r8*rxt(k,396)*y(k,205)
         mat(k,1689) = mat(k,1689) + .400_r8*rxt(k,396)*y(k,190)
         mat(k,1512) = rxt(k,408)*y(k,32)
         mat(k,1264) = -(4._r8*rxt(k,291)*y(k,199) + rxt(k,292)*y(k,200) + rxt(k,293) &
                      *y(k,205) + rxt(k,294)*y(k,124) + rxt(k,305)*y(k,125) + rxt(k,332) &
                      *y(k,211) + rxt(k,365)*y(k,207) + rxt(k,370)*y(k,208) + rxt(k,379) &
                      *y(k,101) + rxt(k,390)*y(k,227))
         mat(k,1782) = -rxt(k,292)*y(k,199)
         mat(k,1731) = -rxt(k,293)*y(k,199)
         mat(k,1387) = -rxt(k,294)*y(k,199)
         mat(k,1982) = -rxt(k,305)*y(k,199)
         mat(k,1189) = -rxt(k,332)*y(k,199)
         mat(k,1216) = -rxt(k,365)*y(k,199)
         mat(k,1146) = -rxt(k,370)*y(k,199)
         mat(k,1076) = -rxt(k,379)*y(k,199)
         mat(k,1054) = -rxt(k,390)*y(k,199)
         mat(k,825) = .060_r8*rxt(k,440)*y(k,135)
         mat(k,1025) = rxt(k,288)*y(k,126) + rxt(k,289)*y(k,219)
         mat(k,1101) = rxt(k,314)*y(k,126) + rxt(k,315)*y(k,219)
         mat(k,433) = .500_r8*rxt(k,296)*y(k,219)
         mat(k,756) = .080_r8*rxt(k,385)*y(k,135)
         mat(k,1092) = .100_r8*rxt(k,338)*y(k,135)
         mat(k,798) = .060_r8*rxt(k,443)*y(k,135)
         mat(k,1166) = .280_r8*rxt(k,352)*y(k,135)
         mat(k,1387) = mat(k,1387) + .530_r8*rxt(k,336)*y(k,211) + rxt(k,345)*y(k,213) &
                      + rxt(k,348)*y(k,215) + rxt(k,323)*y(k,222)
         mat(k,1941) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,211) + rxt(k,346)*y(k,213)
         mat(k,1625) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,98) &
                      + .100_r8*rxt(k,338)*y(k,105) + .060_r8*rxt(k,443)*y(k,110) &
                      + .280_r8*rxt(k,352)*y(k,111)
         mat(k,933) = .650_r8*rxt(k,461)*y(k,219)
         mat(k,1264) = mat(k,1264) + .530_r8*rxt(k,332)*y(k,211)
         mat(k,1782) = mat(k,1782) + .260_r8*rxt(k,333)*y(k,211) + rxt(k,342)*y(k,213) &
                      + .300_r8*rxt(k,321)*y(k,222)
         mat(k,1731) = mat(k,1731) + .450_r8*rxt(k,343)*y(k,213) + .200_r8*rxt(k,347) &
                      *y(k,215) + .150_r8*rxt(k,322)*y(k,222)
         mat(k,1189) = mat(k,1189) + .530_r8*rxt(k,336)*y(k,124) + .530_r8*rxt(k,335) &
                      *y(k,126) + .530_r8*rxt(k,332)*y(k,199) + .260_r8*rxt(k,333) &
                      *y(k,200)
         mat(k,1234) = rxt(k,345)*y(k,124) + rxt(k,346)*y(k,126) + rxt(k,342)*y(k,200) &
                      + .450_r8*rxt(k,343)*y(k,205) + 4.000_r8*rxt(k,344)*y(k,213)
         mat(k,574) = rxt(k,348)*y(k,124) + .200_r8*rxt(k,347)*y(k,205)
         mat(k,1565) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,180)
         mat(k,1015) = rxt(k,323)*y(k,124) + .300_r8*rxt(k,321)*y(k,200) &
                      + .150_r8*rxt(k,322)*y(k,205)
         mat(k,1790) = -(rxt(k,181)*y(k,59) + (4._r8*rxt(k,258) + 4._r8*rxt(k,259) &
                      ) * y(k,200) + rxt(k,260)*y(k,205) + rxt(k,261)*y(k,124) &
                      + rxt(k,281)*y(k,196) + rxt(k,292)*y(k,199) + rxt(k,309) &
                      *y(k,197) + rxt(k,321)*y(k,222) + rxt(k,333)*y(k,211) + rxt(k,342) &
                      *y(k,213) + rxt(k,366)*y(k,207) + rxt(k,371)*y(k,208) + rxt(k,380) &
                      *y(k,101) + rxt(k,391)*y(k,227) + rxt(k,445)*y(k,217) + rxt(k,450) &
                      *y(k,223) + rxt(k,455)*y(k,224))
         mat(k,1893) = -rxt(k,181)*y(k,200)
         mat(k,1740) = -rxt(k,260)*y(k,200)
         mat(k,1396) = -rxt(k,261)*y(k,200)
         mat(k,720) = -rxt(k,281)*y(k,200)
         mat(k,1271) = -rxt(k,292)*y(k,200)
         mat(k,700) = -rxt(k,309)*y(k,200)
         mat(k,1019) = -rxt(k,321)*y(k,200)
         mat(k,1195) = -rxt(k,333)*y(k,200)
         mat(k,1241) = -rxt(k,342)*y(k,200)
         mat(k,1223) = -rxt(k,366)*y(k,200)
         mat(k,1153) = -rxt(k,371)*y(k,200)
         mat(k,1082) = -rxt(k,380)*y(k,200)
         mat(k,1059) = -rxt(k,391)*y(k,200)
         mat(k,1005) = -rxt(k,445)*y(k,200)
         mat(k,986) = -rxt(k,450)*y(k,200)
         mat(k,966) = -rxt(k,455)*y(k,200)
         mat(k,887) = .280_r8*rxt(k,308)*y(k,135)
         mat(k,477) = rxt(k,295)*y(k,219)
         mat(k,371) = .700_r8*rxt(k,263)*y(k,219)
         mat(k,760) = .050_r8*rxt(k,385)*y(k,135)
         mat(k,1082) = mat(k,1082) + rxt(k,379)*y(k,199)
         mat(k,1396) = mat(k,1396) + rxt(k,294)*y(k,199) + .830_r8*rxt(k,411)*y(k,201) &
                      + .170_r8*rxt(k,417)*y(k,214)
         mat(k,1634) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,98)
         mat(k,1271) = mat(k,1271) + rxt(k,379)*y(k,101) + rxt(k,294)*y(k,124) &
                      + 4.000_r8*rxt(k,291)*y(k,199) + .900_r8*rxt(k,292)*y(k,200) &
                      + .450_r8*rxt(k,293)*y(k,205) + rxt(k,365)*y(k,207) + rxt(k,370) &
                      *y(k,208) + rxt(k,332)*y(k,211) + rxt(k,341)*y(k,213) &
                      + rxt(k,390)*y(k,227)
         mat(k,1790) = mat(k,1790) + .900_r8*rxt(k,292)*y(k,199)
         mat(k,655) = .830_r8*rxt(k,411)*y(k,124) + .330_r8*rxt(k,410)*y(k,205)
         mat(k,1740) = mat(k,1740) + .450_r8*rxt(k,293)*y(k,199) + .330_r8*rxt(k,410) &
                      *y(k,201) + .070_r8*rxt(k,416)*y(k,214)
         mat(k,1223) = mat(k,1223) + rxt(k,365)*y(k,199)
         mat(k,1153) = mat(k,1153) + rxt(k,370)*y(k,199)
         mat(k,1195) = mat(k,1195) + rxt(k,332)*y(k,199)
         mat(k,1241) = mat(k,1241) + rxt(k,341)*y(k,199)
         mat(k,776) = .170_r8*rxt(k,417)*y(k,124) + .070_r8*rxt(k,416)*y(k,205)
         mat(k,1574) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,263)*y(k,53)
         mat(k,1059) = mat(k,1059) + rxt(k,390)*y(k,199)
         mat(k,649) = -(rxt(k,410)*y(k,205) + rxt(k,411)*y(k,124) + rxt(k,412) &
                      *y(k,125))
         mat(k,1698) = -rxt(k,410)*y(k,201)
         mat(k,1356) = -rxt(k,411)*y(k,201)
         mat(k,1971) = -rxt(k,412)*y(k,201)
         mat(k,479) = -((rxt(k,329) + rxt(k,330)) * y(k,124))
         mat(k,1344) = -(rxt(k,329) + rxt(k,330)) * y(k,202)
         mat(k,287) = rxt(k,328)*y(k,219)
         mat(k,1504) = rxt(k,328)*y(k,16)
         mat(k,1330) = .750_r8*rxt(k,298)*y(k,204)
         mat(k,603) = .750_r8*rxt(k,298)*y(k,124)
         mat(k,604) = -(rxt(k,297)*y(k,205) + rxt(k,298)*y(k,124))
         mat(k,1694) = -rxt(k,297)*y(k,204)
         mat(k,1352) = -rxt(k,298)*y(k,204)
         mat(k,468) = rxt(k,304)*y(k,219)
         mat(k,1518) = rxt(k,304)*y(k,25)
         mat(k,1739) = -((rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,76) + rxt(k,138) &
                      *y(k,134) + rxt(k,139)*y(k,135) + rxt(k,143)*y(k,219) &
                      + 4._r8*rxt(k,148)*y(k,205) + rxt(k,158)*y(k,126) + rxt(k,163) &
                      *y(k,124) + rxt(k,168)*y(k,125) + (rxt(k,178) + rxt(k,179) &
                      ) * y(k,56) + rxt(k,185)*y(k,59) + rxt(k,211)*y(k,17) + rxt(k,217) &
                      *y(k,19) + rxt(k,254)*y(k,42) + rxt(k,260)*y(k,200) + rxt(k,268) &
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
         mat(k,1308) = -(rxt(k,134) + rxt(k,135) + rxt(k,136)) * y(k,205)
         mat(k,1843) = -rxt(k,138)*y(k,205)
         mat(k,1633) = -rxt(k,139)*y(k,205)
         mat(k,1573) = -rxt(k,143)*y(k,205)
         mat(k,1949) = -rxt(k,158)*y(k,205)
         mat(k,1395) = -rxt(k,163)*y(k,205)
         mat(k,1990) = -rxt(k,168)*y(k,205)
         mat(k,2024) = -(rxt(k,178) + rxt(k,179)) * y(k,205)
         mat(k,1892) = -rxt(k,185)*y(k,205)
         mat(k,1283) = -rxt(k,211)*y(k,205)
         mat(k,1813) = -rxt(k,217)*y(k,205)
         mat(k,1866) = -rxt(k,254)*y(k,205)
         mat(k,1789) = -rxt(k,260)*y(k,205)
         mat(k,365) = -rxt(k,268)*y(k,205)
         mat(k,719) = -rxt(k,282)*y(k,205)
         mat(k,1270) = -rxt(k,293)*y(k,205)
         mat(k,608) = -rxt(k,297)*y(k,205)
         mat(k,699) = -rxt(k,310)*y(k,205)
         mat(k,663) = -rxt(k,318)*y(k,205)
         mat(k,1018) = -rxt(k,322)*y(k,205)
         mat(k,1194) = -rxt(k,334)*y(k,205)
         mat(k,1240) = -rxt(k,343)*y(k,205)
         mat(k,577) = -rxt(k,347)*y(k,205)
         mat(k,843) = -rxt(k,357)*y(k,205)
         mat(k,1222) = -rxt(k,367)*y(k,205)
         mat(k,1152) = -rxt(k,372)*y(k,205)
         mat(k,1081) = -rxt(k,381)*y(k,205)
         mat(k,1058) = -rxt(k,392)*y(k,205)
         mat(k,430) = -rxt(k,396)*y(k,205)
         mat(k,409) = -rxt(k,399)*y(k,205)
         mat(k,360) = -rxt(k,403)*y(k,205)
         mat(k,547) = -rxt(k,406)*y(k,205)
         mat(k,654) = -rxt(k,410)*y(k,205)
         mat(k,615) = -rxt(k,413)*y(k,205)
         mat(k,775) = -rxt(k,416)*y(k,205)
         mat(k,379) = -rxt(k,419)*y(k,205)
         mat(k,630) = -rxt(k,426)*y(k,205)
         mat(k,647) = -rxt(k,432)*y(k,205)
         mat(k,417) = -rxt(k,435)*y(k,205)
         mat(k,1004) = -rxt(k,446)*y(k,205)
         mat(k,985) = -rxt(k,451)*y(k,205)
         mat(k,965) = -rxt(k,456)*y(k,205)
         mat(k,828) = .570_r8*rxt(k,440)*y(k,135)
         mat(k,122) = .650_r8*rxt(k,398)*y(k,219)
         mat(k,1283) = mat(k,1283) + rxt(k,210)*y(k,42)
         mat(k,1813) = mat(k,1813) + rxt(k,222)*y(k,219)
         mat(k,226) = .350_r8*rxt(k,277)*y(k,219)
         mat(k,473) = .130_r8*rxt(k,279)*y(k,135)
         mat(k,203) = rxt(k,284)*y(k,219)
         mat(k,886) = .280_r8*rxt(k,308)*y(k,135)
         mat(k,1866) = mat(k,1866) + rxt(k,210)*y(k,17) + rxt(k,174)*y(k,56) &
                      + rxt(k,255)*y(k,126) + rxt(k,256)*y(k,134)
         mat(k,91) = rxt(k,290)*y(k,219)
         mat(k,691) = rxt(k,262)*y(k,219)
         mat(k,2024) = mat(k,2024) + rxt(k,174)*y(k,42) + rxt(k,177)*y(k,79)
         mat(k,1892) = mat(k,1892) + rxt(k,181)*y(k,200) + rxt(k,192)*y(k,219)
         mat(k,942) = rxt(k,265)*y(k,219)
         mat(k,148) = .730_r8*rxt(k,409)*y(k,219)
         mat(k,237) = .500_r8*rxt(k,476)*y(k,219)
         mat(k,911) = rxt(k,301)*y(k,219)
         mat(k,768) = rxt(k,302)*y(k,219)
         mat(k,526) = rxt(k,177)*y(k,56) + rxt(k,133)*y(k,134) + rxt(k,142)*y(k,219)
         mat(k,135) = rxt(k,266)*y(k,219)
         mat(k,687) = rxt(k,267)*y(k,219)
         mat(k,924) = rxt(k,331)*y(k,219)
         mat(k,947) = rxt(k,316)*y(k,219)
         mat(k,759) = .370_r8*rxt(k,385)*y(k,135)
         mat(k,501) = .300_r8*rxt(k,376)*y(k,219)
         mat(k,494) = rxt(k,377)*y(k,219)
         mat(k,1081) = mat(k,1081) + rxt(k,382)*y(k,124) + rxt(k,383)*y(k,126) &
                      + rxt(k,379)*y(k,199) + 1.200_r8*rxt(k,380)*y(k,200)
         mat(k,340) = rxt(k,384)*y(k,219)
         mat(k,1096) = .140_r8*rxt(k,338)*y(k,135)
         mat(k,243) = .200_r8*rxt(k,340)*y(k,219)
         mat(k,463) = .500_r8*rxt(k,351)*y(k,219)
         mat(k,801) = .570_r8*rxt(k,443)*y(k,135)
         mat(k,1172) = .280_r8*rxt(k,352)*y(k,135)
         mat(k,299) = rxt(k,388)*y(k,219)
         mat(k,903) = rxt(k,389)*y(k,219)
         mat(k,1395) = mat(k,1395) + rxt(k,382)*y(k,101) + rxt(k,358)*y(k,191) &
                      + rxt(k,400)*y(k,193) + rxt(k,405)*y(k,195) + rxt(k,283) &
                      *y(k,196) + rxt(k,311)*y(k,197) + rxt(k,261)*y(k,200) &
                      + .170_r8*rxt(k,411)*y(k,201) + rxt(k,329)*y(k,202) &
                      + .250_r8*rxt(k,298)*y(k,204) + rxt(k,270)*y(k,206) &
                      + .920_r8*rxt(k,368)*y(k,207) + .920_r8*rxt(k,374)*y(k,208) &
                      + .470_r8*rxt(k,336)*y(k,211) + .400_r8*rxt(k,414)*y(k,212) &
                      + .830_r8*rxt(k,417)*y(k,214) + rxt(k,420)*y(k,220) + rxt(k,319) &
                      *y(k,221) + .900_r8*rxt(k,452)*y(k,223) + .800_r8*rxt(k,457) &
                      *y(k,224) + rxt(k,427)*y(k,225) + rxt(k,393)*y(k,227) &
                      + rxt(k,433)*y(k,228) + rxt(k,436)*y(k,230)
         mat(k,1949) = mat(k,1949) + rxt(k,255)*y(k,42) + rxt(k,383)*y(k,101) &
                      + rxt(k,369)*y(k,207) + rxt(k,375)*y(k,208) + .470_r8*rxt(k,335) &
                      *y(k,211) + rxt(k,161)*y(k,219) + rxt(k,394)*y(k,227)
         mat(k,1843) = mat(k,1843) + rxt(k,256)*y(k,42) + rxt(k,133)*y(k,79)
         mat(k,1633) = mat(k,1633) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,98) + .140_r8*rxt(k,338)*y(k,105) + .570_r8*rxt(k,443) &
                      *y(k,110) + .280_r8*rxt(k,352)*y(k,111) + rxt(k,145)*y(k,219)
         mat(k,131) = .800_r8*rxt(k,421)*y(k,219)
         mat(k,735) = rxt(k,477)*y(k,219)
         mat(k,935) = .200_r8*rxt(k,461)*y(k,219)
         mat(k,143) = .280_r8*rxt(k,429)*y(k,219)
         mat(k,169) = .380_r8*rxt(k,431)*y(k,219)
         mat(k,174) = .630_r8*rxt(k,437)*y(k,219)
         mat(k,843) = mat(k,843) + rxt(k,358)*y(k,124)
         mat(k,409) = mat(k,409) + rxt(k,400)*y(k,124)
         mat(k,360) = mat(k,360) + rxt(k,405)*y(k,124)
         mat(k,719) = mat(k,719) + rxt(k,283)*y(k,124) + 2.400_r8*rxt(k,280)*y(k,196) &
                      + rxt(k,281)*y(k,200)
         mat(k,699) = mat(k,699) + rxt(k,311)*y(k,124) + rxt(k,309)*y(k,200)
         mat(k,1270) = mat(k,1270) + rxt(k,379)*y(k,101) + .900_r8*rxt(k,292)*y(k,200) &
                      + rxt(k,365)*y(k,207) + rxt(k,370)*y(k,208) + .470_r8*rxt(k,332) &
                      *y(k,211) + rxt(k,390)*y(k,227)
         mat(k,1789) = mat(k,1789) + rxt(k,181)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,101) &
                      + rxt(k,261)*y(k,124) + rxt(k,281)*y(k,196) + rxt(k,309) &
                      *y(k,197) + .900_r8*rxt(k,292)*y(k,199) + 4.000_r8*rxt(k,258) &
                      *y(k,200) + rxt(k,366)*y(k,207) + rxt(k,371)*y(k,208) &
                      + .730_r8*rxt(k,333)*y(k,211) + rxt(k,342)*y(k,213) &
                      + .500_r8*rxt(k,445)*y(k,217) + .300_r8*rxt(k,321)*y(k,222) &
                      + rxt(k,450)*y(k,223) + rxt(k,455)*y(k,224) + .800_r8*rxt(k,391) &
                      *y(k,227)
         mat(k,654) = mat(k,654) + .170_r8*rxt(k,411)*y(k,124) + .070_r8*rxt(k,410) &
                      *y(k,205)
         mat(k,484) = rxt(k,329)*y(k,124)
         mat(k,608) = mat(k,608) + .250_r8*rxt(k,298)*y(k,124)
         mat(k,1739) = mat(k,1739) + .070_r8*rxt(k,410)*y(k,201) + .160_r8*rxt(k,413) &
                      *y(k,212) + .330_r8*rxt(k,416)*y(k,214)
         mat(k,365) = mat(k,365) + rxt(k,270)*y(k,124)
         mat(k,1222) = mat(k,1222) + .920_r8*rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126) &
                      + rxt(k,365)*y(k,199) + rxt(k,366)*y(k,200)
         mat(k,1152) = mat(k,1152) + .920_r8*rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126) &
                      + rxt(k,370)*y(k,199) + rxt(k,371)*y(k,200)
         mat(k,1194) = mat(k,1194) + .470_r8*rxt(k,336)*y(k,124) + .470_r8*rxt(k,335) &
                      *y(k,126) + .470_r8*rxt(k,332)*y(k,199) + .730_r8*rxt(k,333) &
                      *y(k,200)
         mat(k,615) = mat(k,615) + .400_r8*rxt(k,414)*y(k,124) + .160_r8*rxt(k,413) &
                      *y(k,205)
         mat(k,1240) = mat(k,1240) + rxt(k,342)*y(k,200)
         mat(k,775) = mat(k,775) + .830_r8*rxt(k,417)*y(k,124) + .330_r8*rxt(k,416) &
                      *y(k,205)
         mat(k,1004) = mat(k,1004) + .500_r8*rxt(k,445)*y(k,200)
         mat(k,1573) = mat(k,1573) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,222)*y(k,19) &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,290) &
                      *y(k,47) + rxt(k,262)*y(k,52) + rxt(k,192)*y(k,59) + rxt(k,265) &
                      *y(k,62) + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,476) &
                      *y(k,67) + rxt(k,301)*y(k,74) + rxt(k,302)*y(k,75) + rxt(k,142) &
                      *y(k,79) + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,331) &
                      *y(k,93) + rxt(k,316)*y(k,95) + .300_r8*rxt(k,376)*y(k,99) &
                      + rxt(k,377)*y(k,100) + rxt(k,384)*y(k,102) + .200_r8*rxt(k,340) &
                      *y(k,106) + .500_r8*rxt(k,351)*y(k,109) + rxt(k,388)*y(k,115) &
                      + rxt(k,389)*y(k,116) + rxt(k,161)*y(k,126) + rxt(k,145) &
                      *y(k,135) + .800_r8*rxt(k,421)*y(k,143) + rxt(k,477)*y(k,152) &
                      + .200_r8*rxt(k,461)*y(k,180) + .280_r8*rxt(k,429)*y(k,182) &
                      + .380_r8*rxt(k,431)*y(k,184) + .630_r8*rxt(k,437)*y(k,186)
         mat(k,379) = mat(k,379) + rxt(k,420)*y(k,124)
         mat(k,663) = mat(k,663) + rxt(k,319)*y(k,124)
         mat(k,1018) = mat(k,1018) + .300_r8*rxt(k,321)*y(k,200)
         mat(k,985) = mat(k,985) + .900_r8*rxt(k,452)*y(k,124) + rxt(k,450)*y(k,200)
         mat(k,965) = mat(k,965) + .800_r8*rxt(k,457)*y(k,124) + rxt(k,455)*y(k,200)
         mat(k,630) = mat(k,630) + rxt(k,427)*y(k,124)
         mat(k,1058) = mat(k,1058) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126) &
                      + rxt(k,390)*y(k,199) + .800_r8*rxt(k,391)*y(k,200)
         mat(k,647) = mat(k,647) + rxt(k,433)*y(k,124)
         mat(k,417) = mat(k,417) + rxt(k,436)*y(k,124)
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
         mat(k,362) = -(rxt(k,268)*y(k,205) + rxt(k,270)*y(k,124))
         mat(k,1673) = -rxt(k,268)*y(k,206)
         mat(k,1336) = -rxt(k,270)*y(k,206)
         mat(k,1854) = rxt(k,254)*y(k,205)
         mat(k,1673) = mat(k,1673) + rxt(k,254)*y(k,42)
         mat(k,1214) = -(rxt(k,365)*y(k,199) + rxt(k,366)*y(k,200) + rxt(k,367) &
                      *y(k,205) + rxt(k,368)*y(k,124) + rxt(k,369)*y(k,126))
         mat(k,1262) = -rxt(k,365)*y(k,207)
         mat(k,1780) = -rxt(k,366)*y(k,207)
         mat(k,1729) = -rxt(k,367)*y(k,207)
         mat(k,1385) = -rxt(k,368)*y(k,207)
         mat(k,1939) = -rxt(k,369)*y(k,207)
         mat(k,755) = .600_r8*rxt(k,386)*y(k,219)
         mat(k,1563) = .600_r8*rxt(k,386)*y(k,98)
         mat(k,1142) = -(rxt(k,370)*y(k,199) + rxt(k,371)*y(k,200) + rxt(k,372) &
                      *y(k,205) + rxt(k,374)*y(k,124) + rxt(k,375)*y(k,126))
         mat(k,1259) = -rxt(k,370)*y(k,208)
         mat(k,1777) = -rxt(k,371)*y(k,208)
         mat(k,1726) = -rxt(k,372)*y(k,208)
         mat(k,1382) = -rxt(k,374)*y(k,208)
         mat(k,1936) = -rxt(k,375)*y(k,208)
         mat(k,753) = .400_r8*rxt(k,386)*y(k,219)
         mat(k,1560) = .400_r8*rxt(k,386)*y(k,98)
         mat(k,59) = -(rxt(k,503)*y(k,205) + rxt(k,504)*y(k,124))
         mat(k,1651) = -rxt(k,503)*y(k,209)
         mat(k,1323) = -rxt(k,504)*y(k,209)
         mat(k,746) = rxt(k,506)*y(k,219)
         mat(k,1438) = rxt(k,506)*y(k,98)
         mat(k,65) = -(rxt(k,507)*y(k,205) + rxt(k,508)*y(k,124))
         mat(k,1652) = -rxt(k,507)*y(k,210)
         mat(k,1324) = -rxt(k,508)*y(k,210)
         mat(k,66) = rxt(k,509)*y(k,219)
         mat(k,1439) = rxt(k,509)*y(k,104)
         mat(k,1187) = -(rxt(k,332)*y(k,199) + rxt(k,333)*y(k,200) + rxt(k,334) &
                      *y(k,205) + rxt(k,335)*y(k,126) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,124))
         mat(k,1261) = -rxt(k,332)*y(k,211)
         mat(k,1779) = -rxt(k,333)*y(k,211)
         mat(k,1728) = -rxt(k,334)*y(k,211)
         mat(k,1938) = -rxt(k,335)*y(k,211)
         mat(k,1384) = -(rxt(k,336) + rxt(k,337)) * y(k,211)
         mat(k,1090) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,240) = .200_r8*rxt(k,340)*y(k,219)
         mat(k,1165) = rxt(k,353)*y(k,219)
         mat(k,1562) = .500_r8*rxt(k,339)*y(k,105) + .200_r8*rxt(k,340)*y(k,106) &
                      + rxt(k,353)*y(k,111)
         mat(k,611) = -(rxt(k,413)*y(k,205) + rxt(k,414)*y(k,124) + rxt(k,415) &
                      *y(k,125))
         mat(k,1695) = -rxt(k,413)*y(k,212)
         mat(k,1353) = -rxt(k,414)*y(k,212)
         mat(k,1970) = -rxt(k,415)*y(k,212)
         mat(k,1233) = -(rxt(k,341)*y(k,199) + rxt(k,342)*y(k,200) + rxt(k,343) &
                      *y(k,205) + 4._r8*rxt(k,344)*y(k,213) + rxt(k,345)*y(k,124) &
                      + rxt(k,346)*y(k,126) + rxt(k,354)*y(k,125))
         mat(k,1263) = -rxt(k,341)*y(k,213)
         mat(k,1781) = -rxt(k,342)*y(k,213)
         mat(k,1730) = -rxt(k,343)*y(k,213)
         mat(k,1386) = -rxt(k,345)*y(k,213)
         mat(k,1940) = -rxt(k,346)*y(k,213)
         mat(k,1981) = -rxt(k,354)*y(k,213)
         mat(k,1091) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,241) = .500_r8*rxt(k,340)*y(k,219)
         mat(k,1564) = .500_r8*rxt(k,339)*y(k,105) + .500_r8*rxt(k,340)*y(k,106)
         mat(k,770) = -(rxt(k,416)*y(k,205) + rxt(k,417)*y(k,124) + rxt(k,418) &
                      *y(k,125))
         mat(k,1708) = -rxt(k,416)*y(k,214)
         mat(k,1363) = -rxt(k,417)*y(k,214)
         mat(k,1975) = -rxt(k,418)*y(k,214)
         mat(k,572) = -(rxt(k,347)*y(k,205) + rxt(k,348)*y(k,124))
         mat(k,1691) = -rxt(k,347)*y(k,215)
         mat(k,1351) = -rxt(k,348)*y(k,215)
         mat(k,439) = rxt(k,349)*y(k,219)
         mat(k,245) = rxt(k,350)*y(k,219)
         mat(k,1515) = rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108)
         mat(k,73) = -(rxt(k,511)*y(k,205) + rxt(k,512)*y(k,124))
         mat(k,1653) = -rxt(k,511)*y(k,216)
         mat(k,1325) = -rxt(k,512)*y(k,216)
         mat(k,783) = rxt(k,514)*y(k,219)
         mat(k,1441) = rxt(k,514)*y(k,110)
         mat(k,997) = -(rxt(k,445)*y(k,200) + rxt(k,446)*y(k,205) + rxt(k,447) &
                      *y(k,124) + rxt(k,448)*y(k,126))
         mat(k,1770) = -rxt(k,445)*y(k,217)
         mat(k,1718) = -rxt(k,446)*y(k,217)
         mat(k,1375) = -rxt(k,447)*y(k,217)
         mat(k,1928) = -rxt(k,448)*y(k,217)
         mat(k,822) = rxt(k,439)*y(k,126)
         mat(k,795) = rxt(k,442)*y(k,126)
         mat(k,1928) = mat(k,1928) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,110) &
                      + .500_r8*rxt(k,459)*y(k,179)
         mat(k,314) = rxt(k,449)*y(k,219)
         mat(k,869) = .500_r8*rxt(k,459)*y(k,126)
         mat(k,1551) = rxt(k,449)*y(k,128)
         mat(k,1416) = -(rxt(k,124)*y(k,77) + rxt(k,125)*y(k,231) + rxt(k,128) &
                      *y(k,135) + (rxt(k,206) + rxt(k,207)) * y(k,85) + (rxt(k,229) &
                      + rxt(k,230)) * y(k,81) + rxt(k,235)*y(k,64) + rxt(k,236) &
                      *y(k,65) + rxt(k,274)*y(k,86))
         mat(k,1038) = -rxt(k,124)*y(k,218)
         mat(k,2068) = -rxt(k,125)*y(k,218)
         mat(k,1630) = -rxt(k,128)*y(k,218)
         mat(k,2043) = -(rxt(k,206) + rxt(k,207)) * y(k,218)
         mat(k,707) = -(rxt(k,229) + rxt(k,230)) * y(k,218)
         mat(k,113) = -rxt(k,235)*y(k,218)
         mat(k,158) = -rxt(k,236)*y(k,218)
         mat(k,133) = -rxt(k,274)*y(k,218)
         mat(k,1571) = -(rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) + rxt(k,143)*y(k,205) &
                      + rxt(k,144)*y(k,134) + rxt(k,145)*y(k,135) + (4._r8*rxt(k,146) &
                      + 4._r8*rxt(k,147)) * y(k,219) + rxt(k,149)*y(k,90) + rxt(k,161) &
                      *y(k,126) + rxt(k,162)*y(k,112) + rxt(k,170)*y(k,125) + rxt(k,171) &
                      *y(k,89) + rxt(k,190)*y(k,60) + (rxt(k,192) + rxt(k,193) &
                      ) * y(k,59) + rxt(k,195)*y(k,85) + rxt(k,198)*y(k,92) + rxt(k,222) &
                      *y(k,19) + rxt(k,224)*y(k,81) + rxt(k,257)*y(k,42) + rxt(k,262) &
                      *y(k,52) + rxt(k,263)*y(k,53) + (rxt(k,265) + rxt(k,275) &
                      ) * y(k,62) + rxt(k,266)*y(k,86) + rxt(k,267)*y(k,87) + rxt(k,277) &
                      *y(k,24) + rxt(k,284)*y(k,26) + rxt(k,285)*y(k,27) + rxt(k,287) &
                      *y(k,28) + rxt(k,289)*y(k,45) + rxt(k,290)*y(k,47) + rxt(k,295) &
                      *y(k,50) + rxt(k,296)*y(k,51) + rxt(k,301)*y(k,74) + rxt(k,302) &
                      *y(k,75) + rxt(k,303)*y(k,140) + rxt(k,304)*y(k,25) + rxt(k,312) &
                      *y(k,30) + rxt(k,313)*y(k,31) + rxt(k,315)*y(k,49) + rxt(k,316) &
                      *y(k,95) + rxt(k,317)*y(k,127) + rxt(k,320)*y(k,147) + rxt(k,324) &
                      *y(k,148) + rxt(k,325)*y(k,29) + rxt(k,326)*y(k,48) + rxt(k,328) &
                      *y(k,16) + rxt(k,331)*y(k,93) + rxt(k,339)*y(k,105) + rxt(k,340) &
                      *y(k,106) + rxt(k,349)*y(k,107) + rxt(k,350)*y(k,108) + rxt(k,351) &
                      *y(k,109) + rxt(k,353)*y(k,111) + rxt(k,356)*y(k,1) + rxt(k,360) &
                      *y(k,2) + rxt(k,361)*y(k,15) + rxt(k,362)*y(k,94) + rxt(k,363) &
                      *y(k,96) + rxt(k,364)*y(k,97) + rxt(k,376)*y(k,99) + rxt(k,377) &
                      *y(k,100) + rxt(k,384)*y(k,102) + rxt(k,386)*y(k,98) + rxt(k,387) &
                      *y(k,103) + rxt(k,388)*y(k,115) + rxt(k,389)*y(k,116) + rxt(k,395) &
                      *y(k,183) + rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8) + rxt(k,402) &
                      *y(k,22) + rxt(k,404)*y(k,23) + rxt(k,408)*y(k,32) + rxt(k,409) &
                      *y(k,66) + rxt(k,421)*y(k,143) + rxt(k,424)*y(k,144) + rxt(k,428) &
                      *y(k,181) + rxt(k,429)*y(k,182) + rxt(k,431)*y(k,184) + rxt(k,434) &
                      *y(k,185) + rxt(k,437)*y(k,186) + rxt(k,438)*y(k,187) + rxt(k,441) &
                      *y(k,6) + rxt(k,444)*y(k,110) + rxt(k,449)*y(k,128) + rxt(k,453) &
                      *y(k,176) + rxt(k,454)*y(k,177) + rxt(k,458)*y(k,178) + rxt(k,460) &
                      *y(k,179) + rxt(k,461)*y(k,180) + (rxt(k,463) + rxt(k,476) &
                      ) * y(k,67) + rxt(k,465)*y(k,138) + rxt(k,470)*y(k,149) &
                      + rxt(k,475)*y(k,151) + rxt(k,477)*y(k,152) + rxt(k,479) &
                      *y(k,120))
         mat(k,1039) = -rxt(k,141)*y(k,219)
         mat(k,525) = -rxt(k,142)*y(k,219)
         mat(k,1737) = -rxt(k,143)*y(k,219)
         mat(k,1841) = -rxt(k,144)*y(k,219)
         mat(k,1631) = -rxt(k,145)*y(k,219)
         mat(k,319) = -rxt(k,149)*y(k,219)
         mat(k,1947) = -rxt(k,161)*y(k,219)
         mat(k,327) = -rxt(k,162)*y(k,219)
         mat(k,1988) = -rxt(k,170)*y(k,219)
         mat(k,1294) = -rxt(k,171)*y(k,219)
         mat(k,857) = -rxt(k,190)*y(k,219)
         mat(k,1890) = -(rxt(k,192) + rxt(k,193)) * y(k,219)
         mat(k,2044) = -rxt(k,195)*y(k,219)
         mat(k,725) = -rxt(k,198)*y(k,219)
         mat(k,1811) = -rxt(k,222)*y(k,219)
         mat(k,708) = -rxt(k,224)*y(k,219)
         mat(k,1864) = -rxt(k,257)*y(k,219)
         mat(k,690) = -rxt(k,262)*y(k,219)
         mat(k,370) = -rxt(k,263)*y(k,219)
         mat(k,941) = -(rxt(k,265) + rxt(k,275)) * y(k,219)
         mat(k,134) = -rxt(k,266)*y(k,219)
         mat(k,686) = -rxt(k,267)*y(k,219)
         mat(k,225) = -rxt(k,277)*y(k,219)
         mat(k,202) = -rxt(k,284)*y(k,219)
         mat(k,276) = -rxt(k,285)*y(k,219)
         mat(k,229) = -rxt(k,287)*y(k,219)
         mat(k,1028) = -rxt(k,289)*y(k,219)
         mat(k,90) = -rxt(k,290)*y(k,219)
         mat(k,476) = -rxt(k,295)*y(k,219)
         mat(k,434) = -rxt(k,296)*y(k,219)
         mat(k,910) = -rxt(k,301)*y(k,219)
         mat(k,767) = -rxt(k,302)*y(k,219)
         mat(k,386) = -rxt(k,303)*y(k,219)
         mat(k,471) = -rxt(k,304)*y(k,219)
         mat(k,345) = -rxt(k,312)*y(k,219)
         mat(k,97) = -rxt(k,313)*y(k,219)
         mat(k,1104) = -rxt(k,315)*y(k,219)
         mat(k,946) = -rxt(k,316)*y(k,219)
         mat(k,741) = -rxt(k,317)*y(k,219)
         mat(k,447) = -rxt(k,320)*y(k,219)
         mat(k,333) = -rxt(k,324)*y(k,219)
         mat(k,884) = -rxt(k,325)*y(k,219)
         mat(k,850) = -rxt(k,326)*y(k,219)
         mat(k,290) = -rxt(k,328)*y(k,219)
         mat(k,923) = -rxt(k,331)*y(k,219)
         mat(k,1094) = -rxt(k,339)*y(k,219)
         mat(k,242) = -rxt(k,340)*y(k,219)
         mat(k,442) = -rxt(k,349)*y(k,219)
         mat(k,248) = -rxt(k,350)*y(k,219)
         mat(k,462) = -rxt(k,351)*y(k,219)
         mat(k,1170) = -rxt(k,353)*y(k,219)
         mat(k,556) = -rxt(k,356)*y(k,219)
         mat(k,568) = -rxt(k,360)*y(k,219)
         mat(k,182) = -rxt(k,361)*y(k,219)
         mat(k,178) = -rxt(k,362)*y(k,219)
         mat(k,272) = -rxt(k,363)*y(k,219)
         mat(k,107) = -rxt(k,364)*y(k,219)
         mat(k,500) = -rxt(k,376)*y(k,219)
         mat(k,493) = -rxt(k,377)*y(k,219)
         mat(k,339) = -rxt(k,384)*y(k,219)
         mat(k,757) = -rxt(k,386)*y(k,219)
         mat(k,584) = -rxt(k,387)*y(k,219)
         mat(k,298) = -rxt(k,388)*y(k,219)
         mat(k,902) = -rxt(k,389)*y(k,219)
         mat(k,155) = -rxt(k,395)*y(k,219)
         mat(k,121) = -rxt(k,398)*y(k,219)
         mat(k,304) = -rxt(k,401)*y(k,219)
         mat(k,191) = -rxt(k,402)*y(k,219)
         mat(k,268) = -rxt(k,404)*y(k,219)
         mat(k,207) = -rxt(k,408)*y(k,219)
         mat(k,147) = -rxt(k,409)*y(k,219)
         mat(k,130) = -rxt(k,421)*y(k,219)
         mat(k,262) = -rxt(k,424)*y(k,219)
         mat(k,522) = -rxt(k,428)*y(k,219)
         mat(k,142) = -rxt(k,429)*y(k,219)
         mat(k,168) = -rxt(k,431)*y(k,219)
         mat(k,600) = -rxt(k,434)*y(k,219)
         mat(k,173) = -rxt(k,437)*y(k,219)
         mat(k,352) = -rxt(k,438)*y(k,219)
         mat(k,826) = -rxt(k,441)*y(k,219)
         mat(k,799) = -rxt(k,444)*y(k,219)
         mat(k,316) = -rxt(k,449)*y(k,219)
         mat(k,510) = -rxt(k,453)*y(k,219)
         mat(k,538) = -rxt(k,454)*y(k,219)
         mat(k,399) = -rxt(k,458)*y(k,219)
         mat(k,870) = -rxt(k,460)*y(k,219)
         mat(k,934) = -rxt(k,461)*y(k,219)
         mat(k,236) = -(rxt(k,463) + rxt(k,476)) * y(k,219)
         mat(k,284) = -rxt(k,465)*y(k,219)
         mat(k,422) = -rxt(k,470)*y(k,219)
         mat(k,1115) = -rxt(k,475)*y(k,219)
         mat(k,734) = -rxt(k,477)*y(k,219)
         mat(k,87) = -rxt(k,479)*y(k,219)
         mat(k,826) = mat(k,826) + .630_r8*rxt(k,440)*y(k,135)
         mat(k,225) = mat(k,225) + .650_r8*rxt(k,277)*y(k,219)
         mat(k,471) = mat(k,471) + .130_r8*rxt(k,279)*y(k,135)
         mat(k,276) = mat(k,276) + .500_r8*rxt(k,285)*y(k,219)
         mat(k,884) = mat(k,884) + .360_r8*rxt(k,308)*y(k,135)
         mat(k,1864) = mat(k,1864) + rxt(k,256)*y(k,134)
         mat(k,370) = mat(k,370) + .300_r8*rxt(k,263)*y(k,219)
         mat(k,2022) = rxt(k,179)*y(k,205)
         mat(k,681) = rxt(k,233)*y(k,231)
         mat(k,1306) = rxt(k,140)*y(k,135) + 2.000_r8*rxt(k,135)*y(k,205)
         mat(k,1039) = mat(k,1039) + rxt(k,132)*y(k,134) + rxt(k,124)*y(k,218)
         mat(k,525) = mat(k,525) + rxt(k,133)*y(k,134)
         mat(k,708) = mat(k,708) + rxt(k,223)*y(k,134) + rxt(k,229)*y(k,218)
         mat(k,2044) = mat(k,2044) + rxt(k,194)*y(k,134) + rxt(k,206)*y(k,218)
         mat(k,134) = mat(k,134) + rxt(k,274)*y(k,218)
         mat(k,670) = rxt(k,225)*y(k,134)
         mat(k,725) = mat(k,725) + rxt(k,197)*y(k,134)
         mat(k,757) = mat(k,757) + .320_r8*rxt(k,385)*y(k,135)
         mat(k,584) = mat(k,584) + .600_r8*rxt(k,387)*y(k,219)
         mat(k,1094) = mat(k,1094) + .240_r8*rxt(k,338)*y(k,135)
         mat(k,242) = mat(k,242) + .100_r8*rxt(k,340)*y(k,219)
         mat(k,799) = mat(k,799) + .630_r8*rxt(k,443)*y(k,135)
         mat(k,1170) = mat(k,1170) + .360_r8*rxt(k,352)*y(k,135)
         mat(k,1393) = rxt(k,163)*y(k,205)
         mat(k,1947) = mat(k,1947) + rxt(k,158)*y(k,205)
         mat(k,1841) = mat(k,1841) + rxt(k,256)*y(k,42) + rxt(k,132)*y(k,77) &
                      + rxt(k,133)*y(k,79) + rxt(k,223)*y(k,81) + rxt(k,194)*y(k,85) &
                      + rxt(k,225)*y(k,91) + rxt(k,197)*y(k,92) + rxt(k,138)*y(k,205)
         mat(k,1631) = mat(k,1631) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,140)*y(k,76) &
                      + .320_r8*rxt(k,385)*y(k,98) + .240_r8*rxt(k,338)*y(k,105) &
                      + .630_r8*rxt(k,443)*y(k,110) + .360_r8*rxt(k,352)*y(k,111) &
                      + rxt(k,139)*y(k,205)
         mat(k,447) = mat(k,447) + .500_r8*rxt(k,320)*y(k,219)
         mat(k,155) = mat(k,155) + .500_r8*rxt(k,395)*y(k,219)
         mat(k,429) = .400_r8*rxt(k,396)*y(k,205)
         mat(k,1268) = .450_r8*rxt(k,293)*y(k,205)
         mat(k,653) = .400_r8*rxt(k,410)*y(k,205)
         mat(k,1737) = mat(k,1737) + rxt(k,179)*y(k,56) + 2.000_r8*rxt(k,135)*y(k,76) &
                      + rxt(k,163)*y(k,124) + rxt(k,158)*y(k,126) + rxt(k,138) &
                      *y(k,134) + rxt(k,139)*y(k,135) + .400_r8*rxt(k,396)*y(k,190) &
                      + .450_r8*rxt(k,293)*y(k,199) + .400_r8*rxt(k,410)*y(k,201) &
                      + .450_r8*rxt(k,343)*y(k,213) + .400_r8*rxt(k,416)*y(k,214) &
                      + .200_r8*rxt(k,347)*y(k,215) + .150_r8*rxt(k,322)*y(k,222)
         mat(k,1238) = .450_r8*rxt(k,343)*y(k,205)
         mat(k,774) = .400_r8*rxt(k,416)*y(k,205)
         mat(k,576) = .200_r8*rxt(k,347)*y(k,205)
         mat(k,1417) = rxt(k,124)*y(k,77) + rxt(k,229)*y(k,81) + rxt(k,206)*y(k,85) &
                      + rxt(k,274)*y(k,86) + 2.000_r8*rxt(k,125)*y(k,231)
         mat(k,1571) = mat(k,1571) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,263)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,103) + .100_r8*rxt(k,340)*y(k,106) + .500_r8*rxt(k,320) &
                      *y(k,147) + .500_r8*rxt(k,395)*y(k,183)
         mat(k,1017) = .150_r8*rxt(k,322)*y(k,205)
         mat(k,2069) = rxt(k,233)*y(k,73) + 2.000_r8*rxt(k,125)*y(k,218)
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
         mat(k,375) = -(rxt(k,419)*y(k,205) + rxt(k,420)*y(k,124))
         mat(k,1675) = -rxt(k,419)*y(k,220)
         mat(k,1337) = -rxt(k,420)*y(k,220)
         mat(k,145) = .200_r8*rxt(k,409)*y(k,219)
         mat(k,128) = .140_r8*rxt(k,421)*y(k,219)
         mat(k,260) = rxt(k,424)*y(k,219)
         mat(k,1489) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,143) &
                      + rxt(k,424)*y(k,144)
         mat(k,658) = -(rxt(k,318)*y(k,205) + rxt(k,319)*y(k,124))
         mat(k,1699) = -rxt(k,318)*y(k,221)
         mat(k,1357) = -rxt(k,319)*y(k,221)
         mat(k,874) = rxt(k,325)*y(k,219)
         mat(k,444) = .500_r8*rxt(k,320)*y(k,219)
         mat(k,1523) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,147)
         mat(k,1013) = -(rxt(k,321)*y(k,200) + rxt(k,322)*y(k,205) + rxt(k,323) &
                      *y(k,124))
         mat(k,1771) = -rxt(k,321)*y(k,222)
         mat(k,1719) = -rxt(k,322)*y(k,222)
         mat(k,1376) = -rxt(k,323)*y(k,222)
         mat(k,823) = .060_r8*rxt(k,440)*y(k,135)
         mat(k,848) = rxt(k,326)*y(k,219)
         mat(k,796) = .060_r8*rxt(k,443)*y(k,135)
         mat(k,1614) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,110)
         mat(k,331) = rxt(k,324)*y(k,219)
         mat(k,932) = .150_r8*rxt(k,461)*y(k,219)
         mat(k,1552) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,148) + .150_r8*rxt(k,461) &
                      *y(k,180)
         mat(k,978) = -(rxt(k,450)*y(k,200) + rxt(k,451)*y(k,205) + rxt(k,452) &
                      *y(k,124))
         mat(k,1769) = -rxt(k,450)*y(k,223)
         mat(k,1717) = -rxt(k,451)*y(k,223)
         mat(k,1374) = -rxt(k,452)*y(k,223)
         mat(k,1927) = .500_r8*rxt(k,459)*y(k,179)
         mat(k,509) = rxt(k,453)*y(k,219)
         mat(k,868) = .500_r8*rxt(k,459)*y(k,126) + rxt(k,460)*y(k,219)
         mat(k,1550) = rxt(k,453)*y(k,176) + rxt(k,460)*y(k,179)
         mat(k,956) = -(rxt(k,455)*y(k,200) + rxt(k,456)*y(k,205) + rxt(k,457) &
                      *y(k,124))
         mat(k,1768) = -rxt(k,455)*y(k,224)
         mat(k,1716) = -rxt(k,456)*y(k,224)
         mat(k,1373) = -rxt(k,457)*y(k,224)
         mat(k,821) = rxt(k,441)*y(k,219)
         mat(k,794) = rxt(k,444)*y(k,219)
         mat(k,398) = rxt(k,458)*y(k,219)
         mat(k,1549) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,110) + rxt(k,458)*y(k,178)
         mat(k,622) = -(rxt(k,426)*y(k,205) + rxt(k,427)*y(k,124))
         mat(k,1696) = -rxt(k,426)*y(k,225)
         mat(k,1354) = -rxt(k,427)*y(k,225)
         mat(k,518) = rxt(k,428)*y(k,219)
         mat(k,141) = .650_r8*rxt(k,429)*y(k,219)
         mat(k,1520) = rxt(k,428)*y(k,181) + .650_r8*rxt(k,429)*y(k,182)
         mat(k,79) = -(rxt(k,517)*y(k,205) + rxt(k,518)*y(k,124))
         mat(k,1654) = -rxt(k,517)*y(k,226)
         mat(k,1326) = -rxt(k,518)*y(k,226)
         mat(k,136) = rxt(k,516)*y(k,219)
         mat(k,1442) = rxt(k,516)*y(k,182)
         mat(k,1052) = -(rxt(k,390)*y(k,199) + rxt(k,391)*y(k,200) + rxt(k,392) &
                      *y(k,205) + rxt(k,393)*y(k,124) + rxt(k,394)*y(k,126))
         mat(k,1255) = -rxt(k,390)*y(k,227)
         mat(k,1773) = -rxt(k,391)*y(k,227)
         mat(k,1722) = -rxt(k,392)*y(k,227)
         mat(k,1378) = -rxt(k,393)*y(k,227)
         mat(k,1931) = -rxt(k,394)*y(k,227)
         mat(k,177) = rxt(k,362)*y(k,219)
         mat(k,271) = rxt(k,363)*y(k,219)
         mat(k,106) = rxt(k,364)*y(k,219)
         mat(k,581) = .400_r8*rxt(k,387)*y(k,219)
         mat(k,154) = .500_r8*rxt(k,395)*y(k,219)
         mat(k,1555) = rxt(k,362)*y(k,94) + rxt(k,363)*y(k,96) + rxt(k,364)*y(k,97) &
                      + .400_r8*rxt(k,387)*y(k,103) + .500_r8*rxt(k,395)*y(k,183)
         mat(k,638) = -(rxt(k,432)*y(k,205) + rxt(k,433)*y(k,124))
         mat(k,1697) = -rxt(k,432)*y(k,228)
         mat(k,1355) = -rxt(k,433)*y(k,228)
         mat(k,165) = .560_r8*rxt(k,431)*y(k,219)
         mat(k,593) = rxt(k,434)*y(k,219)
         mat(k,1521) = .560_r8*rxt(k,431)*y(k,184) + rxt(k,434)*y(k,185)
         mat(k,85) = -(rxt(k,520)*y(k,205) + rxt(k,521)*y(k,124))
         mat(k,1655) = -rxt(k,520)*y(k,229)
         mat(k,1327) = -rxt(k,521)*y(k,229)
         mat(k,160) = rxt(k,519)*y(k,219)
         mat(k,1443) = rxt(k,519)*y(k,184)
         mat(k,412) = -(rxt(k,435)*y(k,205) + rxt(k,436)*y(k,124))
         mat(k,1680) = -rxt(k,435)*y(k,230)
         mat(k,1341) = -rxt(k,436)*y(k,230)
         mat(k,172) = .300_r8*rxt(k,437)*y(k,219)
         mat(k,349) = rxt(k,438)*y(k,219)
         mat(k,1495) = .300_r8*rxt(k,437)*y(k,186) + rxt(k,438)*y(k,187)
         mat(k,2081) = -(rxt(k,125)*y(k,218) + rxt(k,233)*y(k,73) + rxt(k,478) &
                      *y(k,153))
         mat(k,1429) = -rxt(k,125)*y(k,231)
         mat(k,684) = -rxt(k,233)*y(k,231)
         mat(k,199) = -rxt(k,478)*y(k,231)
         mat(k,232) = rxt(k,287)*y(k,219)
         mat(k,347) = rxt(k,312)*y(k,219)
         mat(k,98) = rxt(k,313)*y(k,219)
         mat(k,1876) = rxt(k,257)*y(k,219)
         mat(k,1032) = rxt(k,289)*y(k,219)
         mat(k,852) = rxt(k,326)*y(k,219)
         mat(k,1107) = rxt(k,315)*y(k,219)
         mat(k,478) = rxt(k,295)*y(k,219)
         mat(k,437) = rxt(k,296)*y(k,219)
         mat(k,373) = rxt(k,263)*y(k,219)
         mat(k,1315) = rxt(k,136)*y(k,205)
         mat(k,1045) = rxt(k,141)*y(k,219)
         mat(k,530) = rxt(k,142)*y(k,219)
         mat(k,711) = rxt(k,224)*y(k,219)
         mat(k,2056) = (rxt(k,530)+rxt(k,535))*y(k,91) + (rxt(k,523)+rxt(k,529) &
                       +rxt(k,534))*y(k,92) + rxt(k,195)*y(k,219)
         mat(k,688) = rxt(k,267)*y(k,219)
         mat(k,1301) = rxt(k,171)*y(k,219)
         mat(k,323) = rxt(k,149)*y(k,219)
         mat(k,675) = (rxt(k,530)+rxt(k,535))*y(k,85)
         mat(k,730) = (rxt(k,523)+rxt(k,529)+rxt(k,534))*y(k,85) + rxt(k,198)*y(k,219)
         mat(k,1098) = .500_r8*rxt(k,339)*y(k,219)
         mat(k,88) = rxt(k,479)*y(k,219)
         mat(k,450) = rxt(k,320)*y(k,219)
         mat(k,335) = rxt(k,324)*y(k,219)
         mat(k,1749) = rxt(k,136)*y(k,76) + rxt(k,143)*y(k,219)
         mat(k,1583) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31) &
                      + rxt(k,257)*y(k,42) + rxt(k,289)*y(k,45) + rxt(k,326)*y(k,48) &
                      + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50) + rxt(k,296)*y(k,51) &
                      + rxt(k,263)*y(k,53) + rxt(k,141)*y(k,77) + rxt(k,142)*y(k,79) &
                      + rxt(k,224)*y(k,81) + rxt(k,195)*y(k,85) + rxt(k,267)*y(k,87) &
                      + rxt(k,171)*y(k,89) + rxt(k,149)*y(k,90) + rxt(k,198)*y(k,92) &
                      + .500_r8*rxt(k,339)*y(k,105) + rxt(k,479)*y(k,120) + rxt(k,320) &
                      *y(k,147) + rxt(k,324)*y(k,148) + rxt(k,143)*y(k,205) &
                      + 2.000_r8*rxt(k,146)*y(k,219)
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
         mat(k, 53) = mat(k, 53) + lmat(k, 53)
         mat(k, 59) = mat(k, 59) + lmat(k, 59)
         mat(k, 65) = mat(k, 65) + lmat(k, 65)
         mat(k, 67) = mat(k, 67) + lmat(k, 67)
         mat(k, 73) = mat(k, 73) + lmat(k, 73)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 85) = mat(k, 85) + lmat(k, 85)
         mat(k, 86) = mat(k, 86) + lmat(k, 86)
         mat(k, 89) = mat(k, 89) + lmat(k, 89)
         mat(k, 92) = lmat(k, 92)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = mat(k, 95) + lmat(k, 95)
         mat(k, 99) = lmat(k, 99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = lmat(k, 103)
         mat(k, 104) = lmat(k, 104)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 108) = lmat(k, 108)
         mat(k, 109) = lmat(k, 109)
         mat(k, 110) = lmat(k, 110)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 114) = lmat(k, 114)
         mat(k, 115) = lmat(k, 115)
         mat(k, 117) = mat(k, 117) + lmat(k, 117)
         mat(k, 123) = lmat(k, 123)
         mat(k, 124) = lmat(k, 124)
         mat(k, 125) = lmat(k, 125)
         mat(k, 126) = lmat(k, 126)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 132) = mat(k, 132) + lmat(k, 132)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 144) = mat(k, 144) + lmat(k, 144)
         mat(k, 149) = lmat(k, 149)
         mat(k, 150) = lmat(k, 150)
         mat(k, 151) = lmat(k, 151)
         mat(k, 152) = lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 156) = mat(k, 156) + lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 179) = lmat(k, 179)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = lmat(k, 185)
         mat(k, 186) = lmat(k, 186)
         mat(k, 187) = lmat(k, 187)
         mat(k, 188) = lmat(k, 188)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = lmat(k, 193)
         mat(k, 194) = lmat(k, 194)
         mat(k, 196) = mat(k, 196) + lmat(k, 196)
         mat(k, 197) = lmat(k, 197)
         mat(k, 198) = lmat(k, 198)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 204) = mat(k, 204) + lmat(k, 204)
         mat(k, 205) = lmat(k, 205)
         mat(k, 207) = mat(k, 207) + lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = lmat(k, 210)
         mat(k, 211) = lmat(k, 211)
         mat(k, 212) = lmat(k, 212)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = lmat(k, 214)
         mat(k, 215) = lmat(k, 215)
         mat(k, 216) = lmat(k, 216)
         mat(k, 217) = lmat(k, 217)
         mat(k, 218) = lmat(k, 218)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = lmat(k, 220)
         mat(k, 221) = mat(k, 221) + lmat(k, 221)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 246) = lmat(k, 246)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = mat(k, 248) + lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = lmat(k, 250)
         mat(k, 251) = lmat(k, 251)
         mat(k, 252) = lmat(k, 252)
         mat(k, 253) = lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 257) = lmat(k, 257)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 259) = mat(k, 259) + lmat(k, 259)
         mat(k, 261) = lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 263) = lmat(k, 263)
         mat(k, 264) = lmat(k, 264)
         mat(k, 265) = mat(k, 265) + lmat(k, 265)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 270) = mat(k, 270) + lmat(k, 270)
         mat(k, 273) = mat(k, 273) + lmat(k, 273)
         mat(k, 275) = mat(k, 275) + lmat(k, 275)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = mat(k, 278) + lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 286) = mat(k, 286) + lmat(k, 286)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 297) = lmat(k, 297)
         mat(k, 300) = mat(k, 300) + lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 303) = lmat(k, 303)
         mat(k, 304) = mat(k, 304) + lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 306) = lmat(k, 306)
         mat(k, 307) = lmat(k, 307)
         mat(k, 308) = lmat(k, 308)
         mat(k, 309) = lmat(k, 309)
         mat(k, 310) = lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 313) = lmat(k, 313)
         mat(k, 315) = lmat(k, 315)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 318) = mat(k, 318) + lmat(k, 318)
         mat(k, 319) = mat(k, 319) + lmat(k, 319)
         mat(k, 320) = lmat(k, 320)
         mat(k, 321) = lmat(k, 321)
         mat(k, 322) = mat(k, 322) + lmat(k, 322)
         mat(k, 324) = mat(k, 324) + lmat(k, 324)
         mat(k, 326) = mat(k, 326) + lmat(k, 326)
         mat(k, 328) = mat(k, 328) + lmat(k, 328)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = lmat(k, 334)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 337) = lmat(k, 337)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 350) = lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = mat(k, 352) + lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 356) = mat(k, 356) + lmat(k, 356)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 365) = mat(k, 365) + lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 369) = lmat(k, 369)
         mat(k, 370) = mat(k, 370) + lmat(k, 370)
         mat(k, 372) = mat(k, 372) + lmat(k, 372)
         mat(k, 375) = mat(k, 375) + lmat(k, 375)
         mat(k, 381) = lmat(k, 381)
         mat(k, 382) = lmat(k, 382)
         mat(k, 383) = lmat(k, 383)
         mat(k, 384) = mat(k, 384) + lmat(k, 384)
         mat(k, 385) = lmat(k, 385)
         mat(k, 387) = lmat(k, 387)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = lmat(k, 390)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = lmat(k, 397)
         mat(k, 399) = mat(k, 399) + lmat(k, 399)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = lmat(k, 401)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 424) = lmat(k, 424)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 432) = mat(k, 432) + lmat(k, 432)
         mat(k, 434) = mat(k, 434) + lmat(k, 434)
         mat(k, 435) = lmat(k, 435)
         mat(k, 438) = mat(k, 438) + lmat(k, 438)
         mat(k, 440) = lmat(k, 440)
         mat(k, 441) = lmat(k, 441)
         mat(k, 443) = mat(k, 443) + lmat(k, 443)
         mat(k, 446) = lmat(k, 446)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = lmat(k, 448)
         mat(k, 449) = lmat(k, 449)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 452) = lmat(k, 452)
         mat(k, 453) = lmat(k, 453)
         mat(k, 454) = lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 457) = mat(k, 457) + lmat(k, 457)
         mat(k, 458) = lmat(k, 458)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 461) = lmat(k, 461)
         mat(k, 466) = lmat(k, 466)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 492) = lmat(k, 492)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 499) = lmat(k, 499)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 508) = lmat(k, 508)
         mat(k, 510) = mat(k, 510) + lmat(k, 510)
         mat(k, 511) = lmat(k, 511)
         mat(k, 512) = lmat(k, 512)
         mat(k, 513) = lmat(k, 513)
         mat(k, 514) = lmat(k, 514)
         mat(k, 515) = mat(k, 515) + lmat(k, 515)
         mat(k, 516) = lmat(k, 516)
         mat(k, 520) = lmat(k, 520)
         mat(k, 521) = lmat(k, 521)
         mat(k, 522) = mat(k, 522) + lmat(k, 522)
         mat(k, 523) = lmat(k, 523)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 525) = mat(k, 525) + lmat(k, 525)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 533) = lmat(k, 533)
         mat(k, 534) = lmat(k, 534)
         mat(k, 535) = mat(k, 535) + lmat(k, 535)
         mat(k, 536) = mat(k, 536) + lmat(k, 536)
         mat(k, 537) = lmat(k, 537)
         mat(k, 539) = lmat(k, 539)
         mat(k, 540) = mat(k, 540) + lmat(k, 540)
         mat(k, 543) = mat(k, 543) + lmat(k, 543)
         mat(k, 549) = lmat(k, 549)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 553) = mat(k, 553) + lmat(k, 553)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 557) = lmat(k, 557)
         mat(k, 558) = mat(k, 558) + lmat(k, 558)
         mat(k, 559) = mat(k, 559) + lmat(k, 559)
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
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = lmat(k, 589)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 595) = lmat(k, 595)
         mat(k, 598) = lmat(k, 598)
         mat(k, 600) = mat(k, 600) + lmat(k, 600)
         mat(k, 601) = lmat(k, 601)
         mat(k, 604) = mat(k, 604) + lmat(k, 604)
         mat(k, 611) = mat(k, 611) + lmat(k, 611)
         mat(k, 622) = mat(k, 622) + lmat(k, 622)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 649) = mat(k, 649) + lmat(k, 649)
         mat(k, 658) = mat(k, 658) + lmat(k, 658)
         mat(k, 668) = mat(k, 668) + lmat(k, 668)
         mat(k, 669) = lmat(k, 669)
         mat(k, 670) = mat(k, 670) + lmat(k, 670)
         mat(k, 676) = mat(k, 676) + lmat(k, 676)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 682) = lmat(k, 682)
         mat(k, 685) = mat(k, 685) + lmat(k, 685)
         mat(k, 689) = mat(k, 689) + lmat(k, 689)
         mat(k, 694) = mat(k, 694) + lmat(k, 694)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 705) = mat(k, 705) + lmat(k, 705)
         mat(k, 706) = mat(k, 706) + lmat(k, 706)
         mat(k, 715) = mat(k, 715) + lmat(k, 715)
         mat(k, 724) = mat(k, 724) + lmat(k, 724)
         mat(k, 725) = mat(k, 725) + lmat(k, 725)
         mat(k, 728) = mat(k, 728) + lmat(k, 728)
         mat(k, 732) = mat(k, 732) + lmat(k, 732)
         mat(k, 733) = lmat(k, 733)
         mat(k, 736) = lmat(k, 736)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 740) = lmat(k, 740)
         mat(k, 742) = lmat(k, 742)
         mat(k, 743) = mat(k, 743) + lmat(k, 743)
         mat(k, 748) = mat(k, 748) + lmat(k, 748)
         mat(k, 764) = lmat(k, 764)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 768) = mat(k, 768) + lmat(k, 768)
         mat(k, 770) = mat(k, 770) + lmat(k, 770)
         mat(k, 788) = mat(k, 788) + lmat(k, 788)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 837) = mat(k, 837) + lmat(k, 837)
         mat(k, 847) = mat(k, 847) + lmat(k, 847)
         mat(k, 849) = lmat(k, 849)
         mat(k, 851) = lmat(k, 851)
         mat(k, 854) = mat(k, 854) + lmat(k, 854)
         mat(k, 855) = mat(k, 855) + lmat(k, 855)
         mat(k, 856) = mat(k, 856) + lmat(k, 856)
         mat(k, 859) = mat(k, 859) + lmat(k, 859)
         mat(k, 860) = mat(k, 860) + lmat(k, 860)
         mat(k, 861) = lmat(k, 861)
         mat(k, 862) = mat(k, 862) + lmat(k, 862)
         mat(k, 865) = mat(k, 865) + lmat(k, 865)
         mat(k, 866) = lmat(k, 866)
         mat(k, 867) = lmat(k, 867)
         mat(k, 871) = lmat(k, 871)
         mat(k, 877) = mat(k, 877) + lmat(k, 877)
         mat(k, 892) = lmat(k, 892)
         mat(k, 896) = mat(k, 896) + lmat(k, 896)
         mat(k, 900) = lmat(k, 900)
         mat(k, 903) = mat(k, 903) + lmat(k, 903)
         mat(k, 906) = lmat(k, 906)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 909) = lmat(k, 909)
         mat(k, 911) = mat(k, 911) + lmat(k, 911)
         mat(k, 912) = mat(k, 912) + lmat(k, 912)
         mat(k, 914) = lmat(k, 914)
         mat(k, 915) = lmat(k, 915)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 917) = lmat(k, 917)
         mat(k, 918) = lmat(k, 918)
         mat(k, 920) = lmat(k, 920)
         mat(k, 921) = lmat(k, 921)
         mat(k, 922) = lmat(k, 922)
         mat(k, 924) = mat(k, 924) + lmat(k, 924)
         mat(k, 926) = lmat(k, 926)
         mat(k, 927) = lmat(k, 927)
         mat(k, 929) = mat(k, 929) + lmat(k, 929)
         mat(k, 930) = mat(k, 930) + lmat(k, 930)
         mat(k, 931) = mat(k, 931) + lmat(k, 931)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 935) = mat(k, 935) + lmat(k, 935)
         mat(k, 937) = mat(k, 937) + lmat(k, 937)
         mat(k, 939) = mat(k, 939) + lmat(k, 939)
         mat(k, 943) = mat(k, 943) + lmat(k, 943)
         mat(k, 945) = lmat(k, 945)
         mat(k, 947) = mat(k, 947) + lmat(k, 947)
         mat(k, 948) = lmat(k, 948)
         mat(k, 956) = mat(k, 956) + lmat(k, 956)
         mat(k, 978) = mat(k, 978) + lmat(k, 978)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k,1013) = mat(k,1013) + lmat(k,1013)
         mat(k,1023) = lmat(k,1023)
         mat(k,1024) = mat(k,1024) + lmat(k,1024)
         mat(k,1029) = lmat(k,1029)
         mat(k,1030) = lmat(k,1030)
         mat(k,1035) = mat(k,1035) + lmat(k,1035)
         mat(k,1052) = mat(k,1052) + lmat(k,1052)
         mat(k,1072) = mat(k,1072) + lmat(k,1072)
         mat(k,1087) = mat(k,1087) + lmat(k,1087)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1091) = mat(k,1091) + lmat(k,1091)
         mat(k,1092) = mat(k,1092) + lmat(k,1092)
         mat(k,1096) = mat(k,1096) + lmat(k,1096)
         mat(k,1097) = mat(k,1097) + lmat(k,1097)
         mat(k,1099) = mat(k,1099) + lmat(k,1099)
         mat(k,1100) = mat(k,1100) + lmat(k,1100)
         mat(k,1101) = mat(k,1101) + lmat(k,1101)
         mat(k,1105) = lmat(k,1105)
         mat(k,1109) = lmat(k,1109)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1111) = mat(k,1111) + lmat(k,1111)
         mat(k,1119) = lmat(k,1119)
         mat(k,1124) = lmat(k,1124)
         mat(k,1142) = mat(k,1142) + lmat(k,1142)
         mat(k,1152) = mat(k,1152) + lmat(k,1152)
         mat(k,1159) = lmat(k,1159)
         mat(k,1160) = mat(k,1160) + lmat(k,1160)
         mat(k,1164) = mat(k,1164) + lmat(k,1164)
         mat(k,1166) = mat(k,1166) + lmat(k,1166)
         mat(k,1173) = lmat(k,1173)
         mat(k,1187) = mat(k,1187) + lmat(k,1187)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1233) = mat(k,1233) + lmat(k,1233)
         mat(k,1264) = mat(k,1264) + lmat(k,1264)
         mat(k,1278) = mat(k,1278) + lmat(k,1278)
         mat(k,1291) = mat(k,1291) + lmat(k,1291)
         mat(k,1294) = mat(k,1294) + lmat(k,1294)
         mat(k,1298) = lmat(k,1298)
         mat(k,1304) = mat(k,1304) + lmat(k,1304)
         mat(k,1308) = mat(k,1308) + lmat(k,1308)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1391) = mat(k,1391) + lmat(k,1391)
         mat(k,1398) = mat(k,1398) + lmat(k,1398)
         mat(k,1406) = mat(k,1406) + lmat(k,1406)
         mat(k,1408) = mat(k,1408) + lmat(k,1408)
         mat(k,1409) = mat(k,1409) + lmat(k,1409)
         mat(k,1411) = mat(k,1411) + lmat(k,1411)
         mat(k,1412) = mat(k,1412) + lmat(k,1412)
         mat(k,1414) = mat(k,1414) + lmat(k,1414)
         mat(k,1415) = lmat(k,1415)
         mat(k,1416) = mat(k,1416) + lmat(k,1416)
         mat(k,1417) = mat(k,1417) + lmat(k,1417)
         mat(k,1419) = lmat(k,1419)
         mat(k,1420) = lmat(k,1420)
         mat(k,1422) = lmat(k,1422)
         mat(k,1423) = lmat(k,1423)
         mat(k,1427) = mat(k,1427) + lmat(k,1427)
         mat(k,1448) = lmat(k,1448)
         mat(k,1457) = lmat(k,1457)
         mat(k,1566) = mat(k,1566) + lmat(k,1566)
         mat(k,1571) = mat(k,1571) + lmat(k,1571)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1574) = mat(k,1574) + lmat(k,1574)
         mat(k,1581) = mat(k,1581) + lmat(k,1581)
         mat(k,1583) = mat(k,1583) + lmat(k,1583)
         mat(k,1630) = mat(k,1630) + lmat(k,1630)
         mat(k,1632) = mat(k,1632) + lmat(k,1632)
         mat(k,1636) = mat(k,1636) + lmat(k,1636)
         mat(k,1687) = mat(k,1687) + lmat(k,1687)
         mat(k,1739) = mat(k,1739) + lmat(k,1739)
         mat(k,1790) = mat(k,1790) + lmat(k,1790)
         mat(k,1806) = mat(k,1806) + lmat(k,1806)
         mat(k,1815) = mat(k,1815) + lmat(k,1815)
         mat(k,1816) = mat(k,1816) + lmat(k,1816)
         mat(k,1842) = mat(k,1842) + lmat(k,1842)
         mat(k,1846) = mat(k,1846) + lmat(k,1846)
         mat(k,1857) = mat(k,1857) + lmat(k,1857)
         mat(k,1858) = lmat(k,1858)
         mat(k,1861) = mat(k,1861) + lmat(k,1861)
         mat(k,1870) = mat(k,1870) + lmat(k,1870)
         mat(k,1895) = mat(k,1895) + lmat(k,1895)
         mat(k,1897) = mat(k,1897) + lmat(k,1897)
         mat(k,1900) = mat(k,1900) + lmat(k,1900)
         mat(k,1943) = mat(k,1943) + lmat(k,1943)
         mat(k,1945) = mat(k,1945) + lmat(k,1945)
         mat(k,1952) = mat(k,1952) + lmat(k,1952)
         mat(k,1955) = mat(k,1955) + lmat(k,1955)
         mat(k,1956) = mat(k,1956) + lmat(k,1956)
         mat(k,1984) = mat(k,1984) + lmat(k,1984)
         mat(k,1986) = mat(k,1986) + lmat(k,1986)
         mat(k,1988) = mat(k,1988) + lmat(k,1988)
         mat(k,1993) = mat(k,1993) + lmat(k,1993)
         mat(k,1997) = mat(k,1997) + lmat(k,1997)
         mat(k,2013) = mat(k,2013) + lmat(k,2013)
         mat(k,2017) = lmat(k,2017)
         mat(k,2024) = mat(k,2024) + lmat(k,2024)
         mat(k,2025) = lmat(k,2025)
         mat(k,2032) = mat(k,2032) + lmat(k,2032)
         mat(k,2033) = mat(k,2033) + lmat(k,2033)
         mat(k,2042) = mat(k,2042) + lmat(k,2042)
         mat(k,2054) = mat(k,2054) + lmat(k,2054)
         mat(k,2055) = mat(k,2055) + lmat(k,2055)
         mat(k,2062) = lmat(k,2062)
         mat(k,2066) = lmat(k,2066)
         mat(k,2068) = mat(k,2068) + lmat(k,2068)
         mat(k,2069) = mat(k,2069) + lmat(k,2069)
         mat(k,2074) = lmat(k,2074)
         mat(k,2081) = mat(k,2081) + lmat(k,2081)
         mat(k, 166) = 0._r8
         mat(k, 167) = 0._r8
         mat(k, 267) = 0._r8
         mat(k, 357) = 0._r8
         mat(k, 359) = 0._r8
         mat(k, 378) = 0._r8
         mat(k, 405) = 0._r8
         mat(k, 408) = 0._r8
         mat(k, 416) = 0._r8
         mat(k, 517) = 0._r8
         mat(k, 519) = 0._r8
         mat(k, 545) = 0._r8
         mat(k, 546) = 0._r8
         mat(k, 551) = 0._r8
         mat(k, 552) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 562) = 0._r8
         mat(k, 563) = 0._r8
         mat(k, 567) = 0._r8
         mat(k, 592) = 0._r8
         mat(k, 594) = 0._r8
         mat(k, 596) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 599) = 0._r8
         mat(k, 605) = 0._r8
         mat(k, 607) = 0._r8
         mat(k, 621) = 0._r8
         mat(k, 623) = 0._r8
         mat(k, 625) = 0._r8
         mat(k, 626) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 637) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 641) = 0._r8
         mat(k, 642) = 0._r8
         mat(k, 644) = 0._r8
         mat(k, 646) = 0._r8
         mat(k, 659) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 666) = 0._r8
         mat(k, 673) = 0._r8
         mat(k, 698) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 718) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 763) = 0._r8
         mat(k, 787) = 0._r8
         mat(k, 789) = 0._r8
         mat(k, 797) = 0._r8
         mat(k, 804) = 0._r8
         mat(k, 814) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 824) = 0._r8
         mat(k, 831) = 0._r8
         mat(k, 835) = 0._r8
         mat(k, 836) = 0._r8
         mat(k, 840) = 0._r8
         mat(k, 842) = 0._r8
         mat(k, 844) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 879) = 0._r8
         mat(k, 881) = 0._r8
         mat(k, 882) = 0._r8
         mat(k, 883) = 0._r8
         mat(k, 890) = 0._r8
         mat(k, 891) = 0._r8
         mat(k, 895) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 904) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 928) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 954) = 0._r8
         mat(k, 955) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 960) = 0._r8
         mat(k, 961) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 964) = 0._r8
         mat(k, 968) = 0._r8
         mat(k, 970) = 0._r8
         mat(k, 979) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 981) = 0._r8
         mat(k, 982) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 995) = 0._r8
         mat(k, 996) = 0._r8
         mat(k, 998) = 0._r8
         mat(k, 999) = 0._r8
         mat(k,1000) = 0._r8
         mat(k,1001) = 0._r8
         mat(k,1003) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1022) = 0._r8
         mat(k,1027) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1040) = 0._r8
         mat(k,1042) = 0._r8
         mat(k,1055) = 0._r8
         mat(k,1057) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1070) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1073) = 0._r8
         mat(k,1074) = 0._r8
         mat(k,1075) = 0._r8
         mat(k,1077) = 0._r8
         mat(k,1078) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1093) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1117) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1130) = 0._r8
         mat(k,1132) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1138) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1145) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1161) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1169) = 0._r8
         mat(k,1175) = 0._r8
         mat(k,1176) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1205) = 0._r8
         mat(k,1206) = 0._r8
         mat(k,1207) = 0._r8
         mat(k,1208) = 0._r8
         mat(k,1209) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1213) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1217) = 0._r8
         mat(k,1218) = 0._r8
         mat(k,1220) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1227) = 0._r8
         mat(k,1231) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1236) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1266) = 0._r8
         mat(k,1273) = 0._r8
         mat(k,1275) = 0._r8
         mat(k,1279) = 0._r8
         mat(k,1280) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1287) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1295) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1303) = 0._r8
         mat(k,1305) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1311) = 0._r8
         mat(k,1312) = 0._r8
         mat(k,1313) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1390) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1404) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1413) = 0._r8
         mat(k,1425) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1490) = 0._r8
         mat(k,1510) = 0._r8
         mat(k,1519) = 0._r8
         mat(k,1522) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1556) = 0._r8
         mat(k,1570) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1600) = 0._r8
         mat(k,1603) = 0._r8
         mat(k,1607) = 0._r8
         mat(k,1610) = 0._r8
         mat(k,1611) = 0._r8
         mat(k,1612) = 0._r8
         mat(k,1613) = 0._r8
         mat(k,1616) = 0._r8
         mat(k,1620) = 0._r8
         mat(k,1622) = 0._r8
         mat(k,1623) = 0._r8
         mat(k,1624) = 0._r8
         mat(k,1627) = 0._r8
         mat(k,1642) = 0._r8
         mat(k,1643) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1677) = 0._r8
         mat(k,1688) = 0._r8
         mat(k,1706) = 0._r8
         mat(k,1710) = 0._r8
         mat(k,1711) = 0._r8
         mat(k,1712) = 0._r8
         mat(k,1713) = 0._r8
         mat(k,1715) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1727) = 0._r8
         mat(k,1733) = 0._r8
         mat(k,1736) = 0._r8
         mat(k,1759) = 0._r8
         mat(k,1783) = 0._r8
         mat(k,1784) = 0._r8
         mat(k,1786) = 0._r8
         mat(k,1787) = 0._r8
         mat(k,1788) = 0._r8
         mat(k,1791) = 0._r8
         mat(k,1792) = 0._r8
         mat(k,1795) = 0._r8
         mat(k,1798) = 0._r8
         mat(k,1799) = 0._r8
         mat(k,1807) = 0._r8
         mat(k,1808) = 0._r8
         mat(k,1810) = 0._r8
         mat(k,1812) = 0._r8
         mat(k,1814) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1819) = 0._r8
         mat(k,1822) = 0._r8
         mat(k,1823) = 0._r8
         mat(k,1825) = 0._r8
         mat(k,1831) = 0._r8
         mat(k,1837) = 0._r8
         mat(k,1840) = 0._r8
         mat(k,1844) = 0._r8
         mat(k,1853) = 0._r8
         mat(k,1855) = 0._r8
         mat(k,1862) = 0._r8
         mat(k,1863) = 0._r8
         mat(k,1865) = 0._r8
         mat(k,1867) = 0._r8
         mat(k,1868) = 0._r8
         mat(k,1871) = 0._r8
         mat(k,1873) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1887) = 0._r8
         mat(k,1889) = 0._r8
         mat(k,1891) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1909) = 0._r8
         mat(k,1916) = 0._r8
         mat(k,1922) = 0._r8
         mat(k,1923) = 0._r8
         mat(k,1926) = 0._r8
         mat(k,1929) = 0._r8
         mat(k,1935) = 0._r8
         mat(k,1942) = 0._r8
         mat(k,1944) = 0._r8
         mat(k,1946) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1951) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1957) = 0._r8
         mat(k,1958) = 0._r8
         mat(k,1959) = 0._r8
         mat(k,1969) = 0._r8
         mat(k,1972) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1977) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1979) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1987) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1998) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2000) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2007) = 0._r8
         mat(k,2010) = 0._r8
         mat(k,2012) = 0._r8
         mat(k,2014) = 0._r8
         mat(k,2016) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2020) = 0._r8
         mat(k,2021) = 0._r8
         mat(k,2026) = 0._r8
         mat(k,2027) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2040) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2047) = 0._r8
         mat(k,2048) = 0._r8
         mat(k,2050) = 0._r8
         mat(k,2052) = 0._r8
         mat(k,2053) = 0._r8
         mat(k,2061) = 0._r8
         mat(k,2063) = 0._r8
         mat(k,2064) = 0._r8
         mat(k,2065) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2070) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2072) = 0._r8
         mat(k,2073) = 0._r8
         mat(k,2075) = 0._r8
         mat(k,2076) = 0._r8
         mat(k,2077) = 0._r8
         mat(k,2078) = 0._r8
         mat(k,2079) = 0._r8
         mat(k,2080) = 0._r8
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
         mat(k, 53) = mat(k, 53) - dti(k)
         mat(k, 59) = mat(k, 59) - dti(k)
         mat(k, 65) = mat(k, 65) - dti(k)
         mat(k, 67) = mat(k, 67) - dti(k)
         mat(k, 73) = mat(k, 73) - dti(k)
         mat(k, 79) = mat(k, 79) - dti(k)
         mat(k, 85) = mat(k, 85) - dti(k)
         mat(k, 86) = mat(k, 86) - dti(k)
         mat(k, 89) = mat(k, 89) - dti(k)
         mat(k, 92) = mat(k, 92) - dti(k)
         mat(k, 95) = mat(k, 95) - dti(k)
         mat(k, 99) = mat(k, 99) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 114) = mat(k, 114) - dti(k)
         mat(k, 117) = mat(k, 117) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 127) = mat(k, 127) - dti(k)
         mat(k, 132) = mat(k, 132) - dti(k)
         mat(k, 137) = mat(k, 137) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 149) = mat(k, 149) - dti(k)
         mat(k, 153) = mat(k, 153) - dti(k)
         mat(k, 156) = mat(k, 156) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 170) = mat(k, 170) - dti(k)
         mat(k, 175) = mat(k, 175) - dti(k)
         mat(k, 180) = mat(k, 180) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 186) = mat(k, 186) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 192) = mat(k, 192) - dti(k)
         mat(k, 196) = mat(k, 196) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 204) = mat(k, 204) - dti(k)
         mat(k, 208) = mat(k, 208) - dti(k)
         mat(k, 212) = mat(k, 212) - dti(k)
         mat(k, 218) = mat(k, 218) - dti(k)
         mat(k, 221) = mat(k, 221) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 233) = mat(k, 233) - dti(k)
         mat(k, 239) = mat(k, 239) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 249) = mat(k, 249) - dti(k)
         mat(k, 254) = mat(k, 254) - dti(k)
         mat(k, 259) = mat(k, 259) - dti(k)
         mat(k, 265) = mat(k, 265) - dti(k)
         mat(k, 270) = mat(k, 270) - dti(k)
         mat(k, 273) = mat(k, 273) - dti(k)
         mat(k, 278) = mat(k, 278) - dti(k)
         mat(k, 286) = mat(k, 286) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 306) = mat(k, 306) - dti(k)
         mat(k, 312) = mat(k, 312) - dti(k)
         mat(k, 318) = mat(k, 318) - dti(k)
         mat(k, 324) = mat(k, 324) - dti(k)
         mat(k, 330) = mat(k, 330) - dti(k)
         mat(k, 336) = mat(k, 336) - dti(k)
         mat(k, 342) = mat(k, 342) - dti(k)
         mat(k, 348) = mat(k, 348) - dti(k)
         mat(k, 356) = mat(k, 356) - dti(k)
         mat(k, 362) = mat(k, 362) - dti(k)
         mat(k, 368) = mat(k, 368) - dti(k)
         mat(k, 375) = mat(k, 375) - dti(k)
         mat(k, 381) = mat(k, 381) - dti(k)
         mat(k, 384) = mat(k, 384) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 395) = mat(k, 395) - dti(k)
         mat(k, 404) = mat(k, 404) - dti(k)
         mat(k, 412) = mat(k, 412) - dti(k)
         mat(k, 419) = mat(k, 419) - dti(k)
         mat(k, 426) = mat(k, 426) - dti(k)
         mat(k, 432) = mat(k, 432) - dti(k)
         mat(k, 438) = mat(k, 438) - dti(k)
         mat(k, 443) = mat(k, 443) - dti(k)
         mat(k, 451) = mat(k, 451) - dti(k)
         mat(k, 459) = mat(k, 459) - dti(k)
         mat(k, 467) = mat(k, 467) - dti(k)
         mat(k, 475) = mat(k, 475) - dti(k)
         mat(k, 479) = mat(k, 479) - dti(k)
         mat(k, 487) = mat(k, 487) - dti(k)
         mat(k, 495) = mat(k, 495) - dti(k)
         mat(k, 504) = mat(k, 504) - dti(k)
         mat(k, 515) = mat(k, 515) - dti(k)
         mat(k, 524) = mat(k, 524) - dti(k)
         mat(k, 531) = mat(k, 531) - dti(k)
         mat(k, 535) = mat(k, 535) - dti(k)
         mat(k, 543) = mat(k, 543) - dti(k)
         mat(k, 550) = mat(k, 550) - dti(k)
         mat(k, 561) = mat(k, 561) - dti(k)
         mat(k, 572) = mat(k, 572) - dti(k)
         mat(k, 580) = mat(k, 580) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 604) = mat(k, 604) - dti(k)
         mat(k, 611) = mat(k, 611) - dti(k)
         mat(k, 622) = mat(k, 622) - dti(k)
         mat(k, 638) = mat(k, 638) - dti(k)
         mat(k, 649) = mat(k, 649) - dti(k)
         mat(k, 658) = mat(k, 658) - dti(k)
         mat(k, 668) = mat(k, 668) - dti(k)
         mat(k, 677) = mat(k, 677) - dti(k)
         mat(k, 685) = mat(k, 685) - dti(k)
         mat(k, 689) = mat(k, 689) - dti(k)
         mat(k, 694) = mat(k, 694) - dti(k)
         mat(k, 704) = mat(k, 704) - dti(k)
         mat(k, 715) = mat(k, 715) - dti(k)
         mat(k, 724) = mat(k, 724) - dti(k)
         mat(k, 732) = mat(k, 732) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 748) = mat(k, 748) - dti(k)
         mat(k, 765) = mat(k, 765) - dti(k)
         mat(k, 770) = mat(k, 770) - dti(k)
         mat(k, 788) = mat(k, 788) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 837) = mat(k, 837) - dti(k)
         mat(k, 847) = mat(k, 847) - dti(k)
         mat(k, 855) = mat(k, 855) - dti(k)
         mat(k, 865) = mat(k, 865) - dti(k)
         mat(k, 877) = mat(k, 877) - dti(k)
         mat(k, 896) = mat(k, 896) - dti(k)
         mat(k, 908) = mat(k, 908) - dti(k)
         mat(k, 916) = mat(k, 916) - dti(k)
         mat(k, 930) = mat(k, 930) - dti(k)
         mat(k, 939) = mat(k, 939) - dti(k)
         mat(k, 943) = mat(k, 943) - dti(k)
         mat(k, 956) = mat(k, 956) - dti(k)
         mat(k, 978) = mat(k, 978) - dti(k)
         mat(k, 997) = mat(k, 997) - dti(k)
         mat(k,1013) = mat(k,1013) - dti(k)
         mat(k,1024) = mat(k,1024) - dti(k)
         mat(k,1035) = mat(k,1035) - dti(k)
         mat(k,1052) = mat(k,1052) - dti(k)
         mat(k,1072) = mat(k,1072) - dti(k)
         mat(k,1088) = mat(k,1088) - dti(k)
         mat(k,1100) = mat(k,1100) - dti(k)
         mat(k,1111) = mat(k,1111) - dti(k)
         mat(k,1142) = mat(k,1142) - dti(k)
         mat(k,1164) = mat(k,1164) - dti(k)
         mat(k,1187) = mat(k,1187) - dti(k)
         mat(k,1214) = mat(k,1214) - dti(k)
         mat(k,1233) = mat(k,1233) - dti(k)
         mat(k,1264) = mat(k,1264) - dti(k)
         mat(k,1278) = mat(k,1278) - dti(k)
         mat(k,1291) = mat(k,1291) - dti(k)
         mat(k,1304) = mat(k,1304) - dti(k)
         mat(k,1391) = mat(k,1391) - dti(k)
         mat(k,1416) = mat(k,1416) - dti(k)
         mat(k,1571) = mat(k,1571) - dti(k)
         mat(k,1632) = mat(k,1632) - dti(k)
         mat(k,1739) = mat(k,1739) - dti(k)
         mat(k,1790) = mat(k,1790) - dti(k)
         mat(k,1815) = mat(k,1815) - dti(k)
         mat(k,1846) = mat(k,1846) - dti(k)
         mat(k,1870) = mat(k,1870) - dti(k)
         mat(k,1897) = mat(k,1897) - dti(k)
         mat(k,1955) = mat(k,1955) - dti(k)
         mat(k,1997) = mat(k,1997) - dti(k)
         mat(k,2032) = mat(k,2032) - dti(k)
         mat(k,2055) = mat(k,2055) - dti(k)
         mat(k,2081) = mat(k,2081) - dti(k)
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
